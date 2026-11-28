---
title: 幾何光学
description: 光路長の変分から光線を求め、曲面、レンズ、鏡、屈折率分布媒質を近軸行列と光線方程式で解析する。
notation:
  - symbol: 'L'
    meaning: 光路長
  - symbol: 'S'
    meaning: 長さの次元をもつアイコナール
  - symbol: 'y,p'
    meaning: 近軸光線の高さと換算角
  - symbol: 'R'
    meaning: 曲率中心が進行方向側にあれば正の曲率半径
---

## この章の目標

- 光路長の停留と光線方程式を結び付ける。
- 曲面での屈折と反射から結像式を導く。
- 行列の積で近軸光学系を追跡する。
- 屈折率勾配による屈曲と転回を計算する。

## Fermat の原理

### 到達時間の変分
静止した透明媒質で正の屈折率 $n(\vb r)$ が波長に比べてゆっくり変わるとする。曲線 $\vb r(u)$ に沿った :term[光路長]{en="optical path length" yomi="こうろちょう"} は $L=\int n|\dot{\vb r}|\dd u$ である。$n$ が分散を持たなければ $L/c$ は所要時間にも等しい。:term[Fermat の原理]{en="Fermat's principle" yomi="ふぇるまーのげんり"} は、二点を結ぶ実際の光線の $L$ が固定端の変分に対して停留することを述べる。
$$
\delta L=\left[n\hat{\vb t}\cdot\delta\vb r\right]_a^b+
\int\left(|\dot{\vb r}|\grad n-\frac{\dd(n\hat{\vb t})}{\dd u}\right)\cdot\delta\vb r\,\dd u ,
\qquad \hat{\vb t}=\dot{\vb r}/|\dot{\vb r}|.
\label{eq-optics-optical-variation}
$$
第一項は、$\delta|\dot{\vb r}|=\hat{\vb t}\cdot\delta\dot{\vb r}$ を部分積分した境界項である。固定端ではこれが消える。境界の両側が一様なら各区間は直線で、境界上の交点を動かす変分だけが残る。接線に沿った係数は $n_1\hat{\vb t}_1-n_2\hat{\vb t}_2$ の接線成分なので、これをゼロと置くと :eqref[eq-optics-snell] を得る。反射の場合は両区間の屈折率が等しく、法線速度が逆向きの経路を採る。

:::example[平面境界の停留点]
境界からの距離が $a,b>0$ で接線座標が $0,d$ の二点を結ぶ。交点 $x$ に対し $L=n_1\sqrt{a^2+x^2}+n_2\sqrt{b^2+(d-x)^2}$。$L'=0$ は Snell の式であり、$L''$ は二項とも正なので、ここでは停留点が唯一の最小である。
:::
:::example[特別な等光路]
楕円鏡では二焦点から反射点までの距離の和が一定である。鏡面に沿った変分に対し $L'=0$ となるため、一方の焦点から出た光線は他方へ向かう。この場合は反射点を動かしても光路長が変わらない。
:::

### 光路が停留しても最短とは限らない
反射点の接線座標を $y$ とし、凹面鏡の頂点位置を $z=0$、入射側を $z<0$ とする。実物点と実像候補を同じ $z=-a$ に置いて、鏡の位置を $z(y)=-y^2/(2\rho)+O(y^4)$ とする。片道の長さは $\sqrt{(a+z(y))^2+y^2}=a+y^2(1/a-1/\rho)/2+O(y^4)$。往復では二次係数が $1/a-1/\rho$ である。頂点で光路は停留するが、$a>\rho$ のとき近傍では極大である。一様直線区間の長さの凸性を、曲面上の反射点にそのまま適用できない。

:::example[局所極大となる反射経路]
反射鏡の $\rho=0.50\,\mathrm m$ に対して物点と像候補を同じ $a=1.0\,\mathrm m$ に置く。往復の光路長の二次係数は $1/a-1/\rho=-1\,\mathrm{m^{-1}}$。高さを $1\,\mathrm{cm}$ ずらした反射候補の光路は、二次の近似で頂点を通る光路より $0.10\,\mathrm{mm}$ 短い。頂点の経路は対称性で停留するが、光路長はここで局所極大である。
:::

### 平面を横切る光路の一意性と感度

物点 $(0,-a)$ と像候補 $(X,b)$ を平面 $z=0$ の両側に置き、$a,b>0$ とする。光が面を通る点を $(x,0)$ に固定すると、各一様媒質内の最短線は直線である。従って面を一度だけ横切る経路の中では $L(x)=n_1\sqrt{a^2+x^2}+n_2\sqrt{b^2+(X-x)^2}$ を調べればよい。微分すると $L'=n_1x/\sqrt{a^2+x^2}-n_2(X-x)/\sqrt{b^2+(X-x)^2}$。さらに
$$
L''=\frac{n_1 a^2}{(a^2+x^2)^{3/2}}+\frac{n_2 b^2}{[b^2+(X-x)^2]^{3/2}}>0.
$$
左右の無限遠で $L'$ の符号が逆なので、停留点はただ一つで、この範囲で最小になる。平面を曲面に変えると、曲面位置の二階微分から符号不定の項が生じる。Fermat の原理を最短だけに限定できない理由がここにある。

$X/a,X/b$ が小さい場合は $L=L_0+n_1 x^2/(2a)+n_2(X-x)^2/(2b)+O(X^4)$ である。この二次式を平方完成すると $x_*=(n_2/b)X/(n_1/a+n_2/b)$ が最小点で、変位 $\delta x$ に対する光路の増分は $(n_1/a+n_2/b)\delta x^2/2$。停留点の近くで光路が変化しにくいことを、この係数で定量化できる。

:::example[等しい媒質の通過点]
同じ屈折率の二領域に対して、直線経路を変分で再現する。

入力は入射屈折率 $n_{1}=1$、出射屈折率 $n_{2}=1$、物点までの法線距離 $a=\frac{1}{5}\,\mathrm m$、終点までの法線距離 $b=\frac{3}{10}\,\mathrm m$、両点の横ずれ $X=\frac{1}{100}\,\mathrm m$、通過点のずれ $\delta x=\frac{1}{1000}\,\mathrm m$ である。

曲率には面の両側の光路が加わる。

$$
\begin{aligned}
L_2''&=\frac{n_{2}}{b} + \frac{n_{1}}{a},\\
L_2''/(\mathrm{m^{-1}})&= \frac{\left(1\right)}{\left(\frac{3}{10}\right)} + \frac{\left(1\right)}{\left(\frac{1}{5}\right)}=\frac{25}{3}.
\end{aligned}
$$

二つの区間の換算角を等しくして停留点を求める。

$$
\begin{aligned}
x_*&=\frac{X n_{2}}{b L_2''},\\
x_*/(\mathrm m)&= \frac{\left(\frac{1}{100}\right) \cdot \left(1\right)}{\left(\frac{3}{10}\right) \cdot \left(\frac{25}{3}\right)}=\frac{1}{250}.
\end{aligned}
$$

入射区間の傾きを取る。

$$
\begin{aligned}
\theta_i&=\frac{x_*}{a},\\
\theta_i/(\mathrm{rad})&= \frac{\left(\frac{1}{250}\right)}{\left(\frac{1}{5}\right)}=\frac{1}{50}.
\end{aligned}
$$

出射区間で終点へ向かう傾きを取る。

$$
\begin{aligned}
\theta_t&=\frac{X - x_*}{b},\\
\theta_t/(\mathrm{rad})&= \frac{\left(\frac{1}{100}\right) - \left(\frac{1}{250}\right)}{\left(\frac{3}{10}\right)}=\frac{1}{50}.
\end{aligned}
$$

停留経路からの誤差はずれの二乗で増える。

$$
\begin{aligned}
\Delta L_2&=\frac{\delta x^{2} L_2''}{2},\\
\Delta L_2/(\mathrm m)&= \frac{\left(\frac{1}{1000}\right)^{2} \cdot \left(\frac{25}{3}\right)}{2}\simeq 4.1667\times 10^{-6}.
\end{aligned}
$$

両側の傾きは等しく、通過点は直線上にある。停留点を外れると光路が増える。
:::

:::example[透過側の屈折率を増す]
ガラス側の距離を固定して、通過点と角度の変化を見る。

入力は入射屈折率 $n_{1}=1$、出射屈折率 $n_{2}=\frac{3}{2}$、物点までの法線距離 $a=\frac{1}{5}\,\mathrm m$、終点までの法線距離 $b=\frac{3}{10}\,\mathrm m$、両点の横ずれ $X=\frac{1}{100}\,\mathrm m$、通過点のずれ $\delta x=\frac{1}{1000}\,\mathrm m$ である。

曲率には面の両側の光路が加わる。

$$
\begin{aligned}
L_2''&=\frac{n_{2}}{b} + \frac{n_{1}}{a},\\
L_2''/(\mathrm{m^{-1}})&= \frac{\left(\frac{3}{2}\right)}{\left(\frac{3}{10}\right)} + \frac{\left(1\right)}{\left(\frac{1}{5}\right)}=10.
\end{aligned}
$$

二つの区間の換算角を等しくして停留点を求める。

$$
\begin{aligned}
x_*&=\frac{X n_{2}}{b L_2''},\\
x_*/(\mathrm m)&= \frac{\left(\frac{1}{100}\right) \cdot \left(\frac{3}{2}\right)}{\left(\frac{3}{10}\right) \cdot \left(10\right)}=\frac{1}{200}.
\end{aligned}
$$

入射区間の傾きを取る。

$$
\begin{aligned}
\theta_i&=\frac{x_*}{a},\\
\theta_i/(\mathrm{rad})&= \frac{\left(\frac{1}{200}\right)}{\left(\frac{1}{5}\right)}=\frac{1}{40}.
\end{aligned}
$$

出射区間で終点へ向かう傾きを取る。

$$
\begin{aligned}
\theta_t&=\frac{X - x_*}{b},\\
\theta_t/(\mathrm{rad})&= \frac{\left(\frac{1}{100}\right) - \left(\frac{1}{200}\right)}{\left(\frac{3}{10}\right)}=\frac{1}{60}.
\end{aligned}
$$

停留経路からの誤差はずれの二乗で増える。

$$
\begin{aligned}
\Delta L_2&=\frac{\delta x^{2} L_2''}{2},\\
\Delta L_2/(\mathrm m)&= \frac{\left(\frac{1}{1000}\right)^{2} \cdot \left(10\right)}{2}\simeq 5\times 10^{-6}.
\end{aligned}
$$

透過角が小さいため、同じ屈折率だった前例に比べて通過点が終点側へずれる。この設定では通過点の横座標は両端の中点である。
:::

:::example[境界近くに像を置く]
物点を遠くし、終点は面の近くに置く。

入力は入射屈折率 $n_{1}=1$、出射屈折率 $n_{2}=\frac{3}{2}$、物点までの法線距離 $a=\frac{1}{2}\,\mathrm m$、終点までの法線距離 $b=\frac{1}{20}\,\mathrm m$、両点の横ずれ $X=\frac{1}{100}\,\mathrm m$、通過点のずれ $\delta x=\frac{1}{1000}\,\mathrm m$ である。

曲率には面の両側の光路が加わる。

$$
\begin{aligned}
L_2''&=\frac{n_{2}}{b} + \frac{n_{1}}{a},\\
L_2''/(\mathrm{m^{-1}})&= \frac{\left(\frac{3}{2}\right)}{\left(\frac{1}{20}\right)} + \frac{\left(1\right)}{\left(\frac{1}{2}\right)}=32.
\end{aligned}
$$

二つの区間の換算角を等しくして停留点を求める。

$$
\begin{aligned}
x_*&=\frac{X n_{2}}{b L_2''},\\
x_*/(\mathrm m)&= \frac{\left(\frac{1}{100}\right) \cdot \left(\frac{3}{2}\right)}{\left(\frac{1}{20}\right) \cdot \left(32\right)}=\frac{3}{320}.
\end{aligned}
$$

入射区間の傾きを取る。

$$
\begin{aligned}
\theta_i&=\frac{x_*}{a},\\
\theta_i/(\mathrm{rad})&= \frac{\left(\frac{3}{320}\right)}{\left(\frac{1}{2}\right)}=\frac{3}{160}.
\end{aligned}
$$

出射区間で終点へ向かう傾きを取る。

$$
\begin{aligned}
\theta_t&=\frac{X - x_*}{b},\\
\theta_t/(\mathrm{rad})&= \frac{\left(\frac{1}{100}\right) - \left(\frac{3}{320}\right)}{\left(\frac{1}{20}\right)}=\frac{1}{80}.
\end{aligned}
$$

停留経路からの誤差はずれの二乗で増える。

$$
\begin{aligned}
\Delta L_2&=\frac{\delta x^{2} L_2''}{2},\\
\Delta L_2/(\mathrm m)&= \frac{\left(\frac{1}{1000}\right)^{2} \cdot \left(32\right)}{2}\simeq 1.6\times 10^{-5}.
\end{aligned}
$$

終点に近い区間の曲率の寄与が大きい。通過点の同じずれに対して、光路の増分は前例より大きくなる。
:::

## アイコナール方程式と光線方程式

### 波から光線へ
スカラー Helmholtz 方程式 $\nabla^2U+k_0^2n^2U=0$ に $U=Ae^{ik_0 S}$ を入れる。$A,S$ を実とし、十分滑らかで $A\ne0$ とする。積を微分すると
$$
\nabla^2A+ik_0(2\grad A\cdot\grad S+A\nabla^2S)+k_0^2(n^2-|\grad S|^2)A=0 .
$$
$k_0$ が大きい極限で実部の最高次は :term[アイコナール方程式]{en="eikonal equation" yomi="あいこなーるほうていしき"} を与える。虚部は振幅の輸送式になる。
$$
|\grad S|^2=n^2,\qquad \divg(A^2\grad S)=0 .
\label{eq-optics-eikonal}
$$
$S=\mathrm{const.}$ が波面で、光線の接線は $\hat{\vb t}=\grad S/n$ である。弧長を $s$ とすると $\dd(\grad S)/\dd s=(\operatorname{Hess}S)\grad S/n=\grad(|\grad S|^2)/(2n)=\grad n$。従って
$$
\frac{\dd}{\dd s}(n\hat{\vb t})=\grad n .
\label{eq-optics-ray}
$$
この式は :eqref[eq-optics-optical-variation] の積分項からも得られる。輸送式を側面が光線からなる細い管で積分すると、両端面を通る $nA^2$ の流束が等しい。焦点で管の面積がゼロに近づくと、この近似は破綻する。

:::example[球面波]
一定の $n$ で $S=nr$ と取ると $\grad S=n\hat{\vb r}$ である。輸送式は $\dd(r^2nA^2)/\dd r=0$ を与え、外向き振幅の $1/r$ 依存性を再現する。
:::
:::example[曲がる向き]
:eqref[eq-optics-ray] の接線成分を除くと $n\dd\hat{\vb t}/\dd s=\grad n-(\hat{\vb t}\cdot\grad n)\hat{\vb t}$ になる。横方向に屈折率が上がるなら、光線はその側へ曲がる。
:::

### 規約によらない光線の位相増分
アイコナールに沿って光線を追うと $\dd S/\dd s=\grad S\cdot\hat{\vb t}=n$ である。従って $S$ の増分は光路長に等しい。光線を初期点から積分するときは、各光線の初期位相を波面上で設定する必要がある。交差する複数光線の領域では、同じ位置で異なる $\grad S$ をもつ枝が現れるため、すべてを一つの滑らかな $S$ で記述できない。

### 二次元で Hamilton 形式に直す
前へ進む曲線を $y=y(z)$ とし、光路の被積分関数は $\mathcal L=n(y,z)\sqrt{1+(y')^2}$。横の共役変数は $p=\partial\mathcal L/\partial y'=ny'/\sqrt{1+(y')^2}$ であり、近軸の $n\theta$ の厳密版である。$|p|<n$ で $y'=p/\sqrt{n^2-p^2}$ と逆に解ける。Legendre 変換 $\mathcal H=p y'-\mathcal L=-\sqrt{n^2-p^2}$ に対して $\partial_p\mathcal H=y'$、$-\partial_y\mathcal H=n\partial_y n/\sqrt{n^2-p^2}$。後者は Euler–Lagrange 方程式 $p'=\partial_y\mathcal L$ と一致する。媒質が $z$ に依存しないと $\mathcal H'=\partial_y\mathcal H\,y'+\partial_p\mathcal H\,p'=0$。従って軸方向の光学的な運動量 $\sqrt{n^2-p^2}$ が保存する。非近軸でも、平行な方向に媒質が不変であることに対応する保存量がある。

### 光線管の運動と面積保存
近傍の二光線の差 $(\delta y,\delta p)$ は、Hamilton の方程式を一次化した行列
$$
\frac{\dd}{\dd z}\binom{\delta y}{\delta p}
=\begin{pmatrix}\mathcal H_{py}&\mathcal H_{pp}\\-\mathcal H_{yy}&-\mathcal H_{yp}\end{pmatrix}
\binom{\delta y}{\delta p}
$$
で動く。$n$ が二階微分可能なら混合微分が一致し、トレースはゼロである。二つの差ベクトルを並べた行列式の微分は、このトレースと元の行列式の積となる。従って横座標と換算角の位相空間で面積が保存される。開口による光線の除去がないとき、光束を狭い面積に集めるなら一般に角度の広がりが増える。この関係は ABCD 行列の行列式 1 にも現れる。

:::example[波面から曲がりの初期値を読む]
$n=1.5$ の点で、光線の接線に垂直な屈折率成分の勾配が $0.3\,\mathrm{m^{-1}}$ とする。式 :eqref[eq-optics-ray] から初期曲率の大きさは $0.3/1.5=0.2\,\mathrm{m^{-1}}$ で、局所の曲率半径は $5\,\mathrm m$。接線と同方向の勾配は曲率のこの値に入らない。
:::

### 球面から出る光線管の輸送

高周波展開で $U=Ae^{ik_0S}$ と置くと、$k_0$ の一次の項は $2\grad A\cdot\grad S+A\nabla^2S=0$。$A$ を実数とすれば $A$ を掛けて $\divg(A^2\grad S)=0$ と書ける。光線管の側面では $\grad S$ が面に接するので、二断面間で管内の流量が一定になる。球面対称で $S'=n(r)$ の外向き枝を選ぶ場合、球を横切る流量は $4\pi r^2n(r)A^2$。従って $A(r)/A(r_0)=(r_0/r)\sqrt{n(r_0)/n(r)}$ である。これはスカラー Helmholtz 方程式の先頭振幅である。媒質の不均一性が大きいと $\nabla^2 A$ を無視できず、反射も加わる。

$n(r)=n_0+g(r-r_0)$ とすれば位相関数の増分は積分して $n_0\Delta r+g(\Delta r)^2/2$。中心が同じ球の間では方向は変わらない。実際、$\grad n$ が接線と平行なので $\dd\hat{\vb t}/\dd s$ の式の横成分はゼロである。一方、媒質の変化は位相と輸送に現れる。

:::example[均質な球対称媒質]
点源から離れた短い範囲で球を横切る波を評価する。

入力は最初の半径 $r_{0}=\frac{1}{50}\,\mathrm m$、半径の増分 $\Delta r=\frac{1}{100}\,\mathrm m$、最初の屈折率 $n_{0}=\frac{3}{2}$、屈折率勾配 $g=0\,\mathrm{m^{-1}}$、真空波長 $\lambda_0\simeq 6\times 10^{-7}\,\mathrm m$ である。

球の半径そのものを求める。

$$
\begin{aligned}
r&=\Delta r + r_{0},\\
r/(\mathrm m)&= \left(\frac{1}{100}\right) + \left(\frac{1}{50}\right)=\frac{3}{100}.
\end{aligned}
$$

線形分布を終点で評価する。

$$
\begin{aligned}
n&=\Delta r g + n_{0},\\
n&= \left(\frac{1}{100}\right) \cdot \left(0\right) + \left(\frac{3}{2}\right)=\frac{3}{2}.
\end{aligned}
$$

放射状の光線に沿って局所屈折率を積分する。

$$
\begin{aligned}
\Delta S&=\frac{\Delta r^{2} g}{2} + \Delta r n_{0},\\
\Delta S/(\mathrm m)&= \frac{\left(\frac{1}{100}\right)^{2} \cdot \left(0\right)}{2} + \left(\frac{1}{100}\right) \cdot \left(\frac{3}{2}\right)=\frac{3}{200}.
\end{aligned}
$$

真空波長一つ分の光路が位相一周期に相当する。

$$
\begin{aligned}
\Delta S/\lambda_0&=\frac{\Delta S}{\lambda_0},\\
\Delta S/\lambda_0&\simeq \frac{\left(\frac{3}{200}\right)}{\left(6\times 10^{-7}\right)}=25000.
\end{aligned}
$$

面積増加と局所波数の変化を両方入れて輸送式を評価する。

$$
\begin{aligned}
A(r)/A(r_0)&=\frac{r_{0} \sqrt{\frac{n_{0}}{n}}}{r},\\
A(r)/A(r_0)&= \frac{\left(\frac{1}{50}\right) \cdot \sqrt{\frac{\left(\frac{3}{2}\right)}{\left(\frac{3}{2}\right)}}}{\left(\frac{3}{100}\right)}=\frac{2}{3}.
\end{aligned}
$$

振幅比は逆距離則だけで決まる。
:::

:::example[外で屈折率が高い球]
半径とともに屈折率が増す領域を考える。

入力は最初の半径 $r_{0}=\frac{1}{50}\,\mathrm m$、半径の増分 $\Delta r=\frac{1}{100}\,\mathrm m$、最初の屈折率 $n_{0}=\frac{3}{2}$、屈折率勾配 $g=10\,\mathrm{m^{-1}}$、真空波長 $\lambda_0\simeq 6\times 10^{-7}\,\mathrm m$ である。

球の半径そのものを求める。

$$
\begin{aligned}
r&=\Delta r + r_{0},\\
r/(\mathrm m)&= \left(\frac{1}{100}\right) + \left(\frac{1}{50}\right)=\frac{3}{100}.
\end{aligned}
$$

線形分布を終点で評価する。

$$
\begin{aligned}
n&=\Delta r g + n_{0},\\
n&= \left(\frac{1}{100}\right) \cdot \left(10\right) + \left(\frac{3}{2}\right)=\frac{8}{5}.
\end{aligned}
$$

放射状の光線に沿って局所屈折率を積分する。

$$
\begin{aligned}
\Delta S&=\frac{\Delta r^{2} g}{2} + \Delta r n_{0},\\
\Delta S/(\mathrm m)&= \frac{\left(\frac{1}{100}\right)^{2} \cdot \left(10\right)}{2} + \left(\frac{1}{100}\right) \cdot \left(\frac{3}{2}\right)=\frac{31}{2000}.
\end{aligned}
$$

真空波長一つ分の光路が位相一周期に相当する。

$$
\begin{aligned}
\Delta S/\lambda_0&=\frac{\Delta S}{\lambda_0},\\
\Delta S/\lambda_0&\simeq \frac{\left(\frac{31}{2000}\right)}{\left(6\times 10^{-7}\right)}=\frac{77500}{3}.
\end{aligned}
$$

面積増加と局所波数の変化を両方入れて輸送式を評価する。

$$
\begin{aligned}
A(r)/A(r_0)&=\frac{r_{0} \sqrt{\frac{n_{0}}{n}}}{r},\\
A(r)/A(r_0)&= \frac{\left(\frac{1}{50}\right) \cdot \sqrt{\frac{\left(\frac{3}{2}\right)}{\left(\frac{8}{5}\right)}}}{\left(\frac{3}{100}\right)}\simeq 0.6455.
\end{aligned}
$$

光線は放射状のままでも、位相差は増え、振幅は均質媒質の場合より小さくなる。
:::

:::example[外で屈折率が低い球]
逆向きの勾配を与え、同じ距離だけ進ませる。

入力は最初の半径 $r_{0}=\frac{1}{50}\,\mathrm m$、半径の増分 $\Delta r=\frac{1}{100}\,\mathrm m$、最初の屈折率 $n_{0}=\frac{3}{2}$、屈折率勾配 $g=-10\,\mathrm{m^{-1}}$、真空波長 $\lambda_0\simeq 6\times 10^{-7}\,\mathrm m$ である。

球の半径そのものを求める。

$$
\begin{aligned}
r&=\Delta r + r_{0},\\
r/(\mathrm m)&= \left(\frac{1}{100}\right) + \left(\frac{1}{50}\right)=\frac{3}{100}.
\end{aligned}
$$

線形分布を終点で評価する。

$$
\begin{aligned}
n&=\Delta r g + n_{0},\\
n&= \left(\frac{1}{100}\right) \cdot \left(-10\right) + \left(\frac{3}{2}\right)=\frac{7}{5}.
\end{aligned}
$$

放射状の光線に沿って局所屈折率を積分する。

$$
\begin{aligned}
\Delta S&=\frac{\Delta r^{2} g}{2} + \Delta r n_{0},\\
\Delta S/(\mathrm m)&= \frac{\left(\frac{1}{100}\right)^{2} \cdot \left(-10\right)}{2} + \left(\frac{1}{100}\right) \cdot \left(\frac{3}{2}\right)=\frac{29}{2000}.
\end{aligned}
$$

真空波長一つ分の光路が位相一周期に相当する。

$$
\begin{aligned}
\Delta S/\lambda_0&=\frac{\Delta S}{\lambda_0},\\
\Delta S/\lambda_0&\simeq \frac{\left(\frac{29}{2000}\right)}{\left(6\times 10^{-7}\right)}=\frac{72500}{3}.
\end{aligned}
$$

面積増加と局所波数の変化を両方入れて輸送式を評価する。

$$
\begin{aligned}
A(r)/A(r_0)&=\frac{r_{0} \sqrt{\frac{n_{0}}{n}}}{r},\\
A(r)/A(r_0)&= \frac{\left(\frac{1}{50}\right) \cdot \sqrt{\frac{\left(\frac{3}{2}\right)}{\left(\frac{7}{5}\right)}}}{\left(\frac{3}{100}\right)}\simeq 0.69007.
\end{aligned}
$$

局所波数が下がるため、逆距離則だけで求めた値より振幅が大きい。球の全流量は保存する。
:::

## 球面での屈折と近軸近似

### 符号と近軸展開
光を左から右へ進め、光軸との交点を $z=0$、曲率中心を $z=R$ とする。高さ $y$ での曲面位置は $z=y^2/(2R)+O(y^4/R^3)$ で、右向き法線の傾きは $-y/R+O(y^3/R^3)$ である。これは球面式 $(z-R)^2+y^2=R^2$ を微分して得られる。光線角 $\theta$ を $+z$ 軸から $+y$ へ測り、$p=n\theta$ とする。すべての角が十分小さければ、法線との角の正弦を角で置き換えて
$$
p_2=p_1-\Phi y,\qquad \Phi=(n_2-n_1)/R .
\label{eq-optics-surface-power}
$$
$\Phi$ は :term[面屈折力]{en="surface power" yomi="めんくっせつりょく"} で単位は $\mathrm{m^{-1}}$ である。実物点を表面の左 $a>0$、実像を右 $b>0$ と書くと $\theta_1\simeq y/a$、$\theta_2\simeq-y/b$。従って
$$
n_1/a+n_2/b=\Phi .
$$
虚物点、虚像では対応する距離を負にする。物点の高さが $y_o$ のときは入射勾配が $(y-y_o)/a$、出射勾配が $(y_i-y)/b$ なので、横倍率は $y_i/y_o=-n_1b/(n_2a)$ である。

:::example[空気からガラスの球面]
$n_1=1,n_2=1.5,R=0.10\,\mathrm m$ なら $\Phi=5\,\mathrm{m^{-1}}$。平行入射で $b=0.30\,\mathrm m$、$a=0.40\,\mathrm m$ で $b=0.60\,\mathrm m$ となる。
:::
:::example[平面を通して見た深さ]
$\Phi=0$ なら $b=-(n_2/n_1)a$ である。水中 $a=30\,\mathrm{cm}$ の点を空気側から近軸で見ると、$n_1=4/3,n_2=1$ から見かけの深さは $22.5\,\mathrm{cm}$ である。
:::

:::example[球面の近軸誤差]
$R=0.10\,\mathrm m$、$y=0.010\,\mathrm m$ では球面のたわみの二次項は $0.50\,\mathrm{mm}$、次の項 $y^4/(8R^3)$ は $1.25\,\mathrm{\mu m}$ である。光線の角度の見積りには小さくても、光路を波長未満で比較する計算では残す必要がある。
:::

### 面で接線波数を三次元に保存する
近軸近似の前なら、右向き単位法線を $\hat{\vb n}$ とし、入射光線の光学的な運動量を $\vb p_1=n_1\hat{\vb t}_1$ と置く。接線成分は $\vb p_{\parallel}=\vb p_1-(\vb p_1\cdot\hat{\vb n})\hat{\vb n}$。透過の向きを法線の正側に取ると $\vb p_2=\vb p_\parallel+\sqrt{n_2^2-|\vb p_\parallel|^2}\hat{\vb n}$ である。全反射で実平方根がなければ、反射は $\vb p_r=\vb p_1-2(\vb p_1\cdot\hat{\vb n})\hat{\vb n}$ を取る。曲面を通る光線を追跡するときはまず曲面との交点を求め、そこで法線を計算してからこの式を適用する。近軸式の $\Phi y$ は、この曲面則を光軸近くで展開した項である。

### 平面と球面を光線追跡で扱う

球面の中心を $\vb c$、半径の大きさを $\rho=|R|$ とし、追跡中の直線を $\vb r_0+s\hat{\vb t}$ とする。$\vb w=\vb r_0-\vb c$ と置くと、球面との交点は
$$
s=-\vb w\cdot\hat{\vb t}\pm
\sqrt{(\vb w\cdot\hat{\vb t})^2-|\vb w|^2+\rho^2}
$$
である。これは $|\vb w+s\hat{\vb t}|^2=\rho^2$ を解いた結果である。判別式が負なら光線は球面に達しない。実根があるときは、意図する面の枝と開口に入る正の根のうち最小のものを使う。球全体の最初の交点を無条件に採ると、レンズの裏面に相当する枝を選ぶことがある。平面なら $\hat{\vb n}\cdot(\vb r_0+s\hat{\vb t}-\vb r_{\rm plane})=0$ は一次式である。交点で曲面の法線を媒質 1 から 2 へ向け、先ほどのベクトル屈折式を適用する。

:::example[球面までの距離を二つの根から選ぶ]
中心 $z=100\,\mathrm{mm}$、半径 $100\,\mathrm{mm}$ の球に、$z=-10\,\mathrm{mm}$、$y=60\,\mathrm{mm}$ から $+z$ へ進む光線を入れる。平方根は $\sqrt{100^2-60^2}=80\,\mathrm{mm}$ で、交点までの距離は $110\pm80\,\mathrm{mm}$。入射側の球面が光学面なら $30\,\mathrm{mm}$ の根を選び、交点は $z=20\,\mathrm{mm}$ になる。
:::

### プリズムの最小偏角
外側の屈折率を $n_0$、プリズムを $n_1>n_0$ とし、頂角を $\alpha$ とする。光線が二面を横切るとき、内部角 $r_1,r_2$ は $r_1+r_2=\alpha$ を満たす。外の角を $i_1,i_2$ とすれば、各面の向きの変化を足して偏角は $i_1+i_2-\alpha$ である。ここで二面の入射・出射角とも $0\le r_j<\arcsin(n_0/n_1)$ の経路を考える。

$m=n_1/n_0>1$ とすれば Snell 則は $i_j=f(r_j)=\arcsin(m\sin r_j)$ となる。この関数の二階微分は $m(m^2-1)\sin r/(1-m^2\sin^2 r)^{3/2}$ であり、$r>0$ で正である。従って $r_1+r_2$ を固定した和 $f(r_1)+f(r_2)$ は $r_1=r_2=\alpha/2$ で一意に最小になる。対称な経路が伝播域に存在する場合、最小偏角 $\delta_{\min}$ を測れば
$$
\frac{n_1}{n_0}=\frac{\sin[(\alpha+\delta_{\min})/2]}{\sin(\alpha/2)}
$$
で屈折率が得られる。

:::example[波長ごとの偏角を比べる]
空気中で頂角 $60^\circ$ とし、二波長の屈折率が $1.50,1.52$ なら、最小偏角はそれぞれ $2\arcsin(0.75)-60^\circ\simeq37.18^\circ$、$2\arcsin(0.76)-60^\circ\simeq38.93^\circ$ となる。分光では屈折率差を出射角の差として読む。
:::

### 球面での厳密な軸交点

$R>0$ の凸面に左から平行光を入れ、$0<y<R$ の光線を考える。交点は $z=R-\sqrt{R^2-y^2}$。右向き法線と光軸のなす符号付き角は $-\alpha$、$\alpha=\arcsin(y/R)$ である。$n_2>n_1$ なら屈折後の法線との角は $\beta=\arcsin[(n_1/n_2)\sin\alpha]<\alpha$。従って光線の角は $\beta-\alpha$ で、軸交点は $z+y/\tan(\alpha-\beta)$ になる。$y/R\to0$ では $\alpha-\beta\simeq(1-n_1/n_2)y/R$ なので、この式の限界は $n_2/\Phi$。高さを変えながら厳密な交点を求めれば、焦点が一点になるという近軸近似の誤差を分離できる。

面上の高さが一定であっても、面までに進む軸方向距離は二次で変わる。また、角度の代わりに正接を使う補正は三次で効く。厳密な交点と近軸値を比較するときは、面のたわみと角の両方を厳密に取る。

:::example[軸近くの高さで照射]
空気からガラスの球面へ、細い平行光を入れる。

入力は曲率半径 $R=\frac{1}{10}\,\mathrm m$、入射高さ $y=\frac{1}{1000}\,\mathrm m$、入射側と出射側の屈折率比 $n_1/n_2=\frac{2}{3}$ である。

球面の左側の枝である。

$$
\begin{aligned}
z&=R - \sqrt{R^{2} - y^{2}},\\
z/(\mathrm m)&= \left(\frac{1}{10}\right) - \sqrt{\left(\frac{1}{10}\right)^{2} - \left(\frac{1}{1000}\right)^{2}}\simeq 5.0001\times 10^{-6}.
\end{aligned}
$$

高さと半径の比で法線方向を定める。

$$
\begin{aligned}
\alpha&=\operatorname{asin}{\left(\frac{y}{R} \right)},\\
\alpha/(\mathrm{rad})&= \operatorname{asin}{\left(\frac{\left(\frac{1}{1000}\right)}{\left(\frac{1}{10}\right)} \right)}\simeq 0.01.
\end{aligned}
$$

Snell の式で法線に対する角を求める。

$$
\begin{aligned}
\beta&=\operatorname{asin}{\left(\frac{n_1/n_2 y}{R} \right)},\\
\beta/(\mathrm{rad})&= \operatorname{asin}{\left(\frac{\left(\frac{2}{3}\right) \cdot \left(\frac{1}{1000}\right)}{\left(\frac{1}{10}\right)} \right)}\simeq 0.0066667.
\end{aligned}
$$

曲面位置に屈折後の軸方向距離を足す。

$$
\begin{aligned}
z_{\mathrm{axis}}&=\frac{y}{\tan{\left(\alpha - \beta \right)}} + z,\\
z_{\mathrm{axis}}/(\mathrm m)&\simeq \frac{\left(\frac{1}{1000}\right)}{\tan{\left(\left(0.01000016667\right) - \left(0.00666671605\right) \right)}} + \left(5.000125006\times 10^{-6}\right)\simeq 0.29999.
\end{aligned}
$$

高さをゼロに近づけたときの軸交点を使う。

$$
\begin{aligned}
f'&=\frac{R}{1 - n_1/n_2},\\
f'/(\mathrm m)&= \frac{\left(\frac{1}{10}\right)}{1 - \left(\frac{2}{3}\right)}=\frac{3}{10}.
\end{aligned}
$$

軸交点の差を符号付きで比較する。

$$
\begin{aligned}
\Delta z&=z_{\mathrm{axis}} - f',\\
\Delta z/(\mathrm m)&\simeq \left(0.2999933332\right) - \left(\frac{3}{10}\right)\simeq -6.6668\times 10^{-6}.
\end{aligned}
$$

近軸焦点からのずれは小さい。小さな開口の焦点を基準にするのに適した条件である。
:::

:::example[同じ面の周辺を使う]
面そのものを変えず、光線の高さを十倍にする。

入力は曲率半径 $R=\frac{1}{10}\,\mathrm m$、入射高さ $y=\frac{1}{100}\,\mathrm m$、入射側と出射側の屈折率比 $n_1/n_2=\frac{2}{3}$ である。

球面の左側の枝である。

$$
\begin{aligned}
z&=R - \sqrt{R^{2} - y^{2}},\\
z/(\mathrm m)&= \left(\frac{1}{10}\right) - \sqrt{\left(\frac{1}{10}\right)^{2} - \left(\frac{1}{100}\right)^{2}}\simeq 5.0126\times 10^{-4}.
\end{aligned}
$$

高さと半径の比で法線方向を定める。

$$
\begin{aligned}
\alpha&=\operatorname{asin}{\left(\frac{y}{R} \right)},\\
\alpha/(\mathrm{rad})&= \operatorname{asin}{\left(\frac{\left(\frac{1}{100}\right)}{\left(\frac{1}{10}\right)} \right)}\simeq 0.10017.
\end{aligned}
$$

Snell の式で法線に対する角を求める。

$$
\begin{aligned}
\beta&=\operatorname{asin}{\left(\frac{n_1/n_2 y}{R} \right)},\\
\beta/(\mathrm{rad})&= \operatorname{asin}{\left(\frac{\left(\frac{2}{3}\right) \cdot \left(\frac{1}{100}\right)}{\left(\frac{1}{10}\right)} \right)}\simeq 0.066716.
\end{aligned}
$$

曲面位置に屈折後の軸方向距離を足す。

$$
\begin{aligned}
z_{\mathrm{axis}}&=\frac{y}{\tan{\left(\alpha - \beta \right)}} + z,\\
z_{\mathrm{axis}}/(\mathrm m)&\simeq \frac{\left(\frac{1}{100}\right)}{\tan{\left(\left(0.1001674212\right) - \left(0.06671614841\right) \right)}} + \left(5.012562893\times 10^{-4}\right)\simeq 0.29933.
\end{aligned}
$$

高さをゼロに近づけたときの軸交点を使う。

$$
\begin{aligned}
f'&=\frac{R}{1 - n_1/n_2},\\
f'/(\mathrm m)&= \frac{\left(\frac{1}{10}\right)}{1 - \left(\frac{2}{3}\right)}=\frac{3}{10}.
\end{aligned}
$$

軸交点の差を符号付きで比較する。

$$
\begin{aligned}
\Delta z&=z_{\mathrm{axis}} - f',\\
\Delta z/(\mathrm m)&\simeq \left(0.2993320313\right) - \left(\frac{3}{10}\right)\simeq -6.6797\times 10^{-4}.
\end{aligned}
$$

軸交点の前進は先頭では高さの二乗で増える。焦点までの距離を一点の値で扱うとずれる。
:::

:::example[光学面の端に近い点]
さらに高い光線を通し、法線角が大きい場合を比較する。

入力は曲率半径 $R=\frac{1}{10}\,\mathrm m$、入射高さ $y=\frac{1}{20}\,\mathrm m$、入射側と出射側の屈折率比 $n_1/n_2=\frac{2}{3}$ である。

球面の左側の枝である。

$$
\begin{aligned}
z&=R - \sqrt{R^{2} - y^{2}},\\
z/(\mathrm m)&= \left(\frac{1}{10}\right) - \sqrt{\left(\frac{1}{10}\right)^{2} - \left(\frac{1}{20}\right)^{2}}\simeq 0.013397.
\end{aligned}
$$

高さと半径の比で法線方向を定める。

$$
\begin{aligned}
\alpha&=\operatorname{asin}{\left(\frac{y}{R} \right)},\\
\alpha/(\mathrm{rad})&= \operatorname{asin}{\left(\frac{\left(\frac{1}{20}\right)}{\left(\frac{1}{10}\right)} \right)}\simeq 0.5236.
\end{aligned}
$$

Snell の式で法線に対する角を求める。

$$
\begin{aligned}
\beta&=\operatorname{asin}{\left(\frac{n_1/n_2 y}{R} \right)},\\
\beta/(\mathrm{rad})&= \operatorname{asin}{\left(\frac{\left(\frac{2}{3}\right) \cdot \left(\frac{1}{20}\right)}{\left(\frac{1}{10}\right)} \right)}\simeq 0.33984.
\end{aligned}
$$

曲面位置に屈折後の軸方向距離を足す。

$$
\begin{aligned}
z_{\mathrm{axis}}&=\frac{y}{\tan{\left(\alpha - \beta \right)}} + z,\\
z_{\mathrm{axis}}/(\mathrm m)&\simeq \frac{\left(\frac{1}{20}\right)}{\tan{\left(\left(0.5235987756\right) - \left(0.3398369095\right) \right)}} + \left(0.01339745962\right)\simeq 0.28242.
\end{aligned}
$$

高さをゼロに近づけたときの軸交点を使う。

$$
\begin{aligned}
f'&=\frac{R}{1 - n_1/n_2},\\
f'/(\mathrm m)&= \frac{\left(\frac{1}{10}\right)}{1 - \left(\frac{2}{3}\right)}=\frac{3}{10}.
\end{aligned}
$$

軸交点の差を符号付きで比較する。

$$
\begin{aligned}
\Delta z&=z_{\mathrm{axis}} - f',\\
\Delta z/(\mathrm m)&\simeq \left(0.2824191173\right) - \left(\frac{3}{10}\right)\simeq -0.017581.
\end{aligned}
$$

法線に対する角が大きいため、三次までの角展開による補正だけでも不足し得る。
:::

## 薄レンズと結像公式

薄いレンズでは二面での光線高さを同じとみなし、外側媒質を $n_0$、レンズを $n_\ell$ とする。二面の角変化を足すと $\Phi=(n_\ell-n_0)(1/R_1-1/R_2)$ である。距離の符号を前節と同じに取ると
$$
\frac{1}{f}=\frac{\Phi}{n_0},\qquad
\frac{1}{a}+\frac{1}{b}=\frac{1}{f},\qquad M=-b/a .
\label{eq-optics-thin-lens}
$$
$f$ は焦点距離である。平行な光線の交点は $b=f$ であり、$f<0$ なら出射光線を後ろへ延長して交点を求める。物点が $a=f$ の場合は出射が平行になる。さらに $a=f+x,\ b=f+x'$ を代入して通分すると $xx'=f^2$ が得られる。焦点から測るこの積の式が Newton の結像公式である。

:::example[実像と虚像]
$f=10\,\mathrm{cm}$ で $a=30\,\mathrm{cm}$ なら $b=15\,\mathrm{cm},M=-1/2$ である。同じレンズで $a=5\,\mathrm{cm}$ なら $b=-10\,\mathrm{cm},M=2$ である。後者は物体側に正立の虚像を作る。
:::
:::example[液体中のレンズ]
$n_\ell=1.5,n_0=1.5$ とすると両面とも $\Phi=0$ である。同じ形でも周囲の屈折率が変わると焦点距離は変わる。
:::

### 前後の変位と倍率
共役距離が有限である区間で $1/a+1/b=1/f$ を微分すると $\dd b/\dd a=-(b/a)^2$ である。実物点の実験室座標は $-a$ なので、物点を右へ動かしたとき像の実験室座標の変位は同符号で、倍率は横倍率の二乗である。媒質の屈折率が異なる単一曲面では、$n_1/a+n_2/b=\Phi$ から縦倍率が $n_1 b^2/(n_2 a^2)$ になる。横倍率の式を使うと $(n_2/n_1)M^2$ である。

:::example[集光位置の感度]
$f=10\,\mathrm{cm}$、$a=15\,\mathrm{cm}$ の共役で $b=30\,\mathrm{cm}$、$M=-2$ である。物点を右へ $0.10\,\mathrm{mm}$ だけ移すと、一次の像変位は右へ $0.40\,\mathrm{mm}$ である。
:::

:::example[同じ物像間隔を満たす二配置]
距離 $l$ だけ離れた物体面とスクリーンの間に焦点距離 $f>0$ の薄レンズを置く。$a+b=l$、結像式から $ab=fl$。従って位置は $a=(l\pm\sqrt{l(l-4f)})/2$、二配置のレンズ位置差を $\Delta$ とすれば $f=(l^2-\Delta^2)/(4l)$ となる。$l=1.0\,\mathrm m$、$\Delta=0.60\,\mathrm m$ なら $f=0.16\,\mathrm m$。$l<4f$ では実像となる配置がなく、$l=4f$ で一つだけになる。
:::

### 横倍率から開口角の変化を読む

二本の近軸光線に対し $H=y_1p_2-y_2p_1$ を作る。屈折で各 $p_j$ に $-\Phi y_j$ を足すと、二つの補正は打ち消し合う。同じく距離 $d$ の伝播で各 $y_j$ に $dp_j/n$ を足しても $H$ は変わらない。この $H$ を :term[Lagrange の光学不変量]{en="Lagrange optical invariant" yomi="らぐらんじゅのこうがくふへんりょう"} と呼ぶ。物面で一本を軸上点から出る周辺光線 $(0,p_o)$、もう一本を高さ $h_o$ の点から出る主光線 $(h_o,p_c)$ とする。結像面で前者の高さが再びゼロ、後者の高さが $Mh_o$ となるので、同じ不変量は $-p_o h_o=-p_i Mh_o$。従って $p_i=p_o/M$。薄レンズで物側と像側の媒質が同じなら、像が拡大されるほど開口角は小さくなる。光を多く取り込みながら小さく集光したい場合は、その逆で像側に広い開口が要る。

結像が成立する範囲で観測面を $\Delta b$ だけずらすと、軸上点の周辺光線のずれは $|\Delta b\tan\theta_i|$。近軸ではこれが幾何学的なぼけの半径となる。有限の回折や収差があると、正しい焦点にも別の広がりがあるため、この半径は追加の幾何学的ずれを表す。

:::example[等倍の近くで像面を動かす]
小さな物側開口角の像をセンサーに作る。レンズを固定し、センサーだけをずらす。

入力は物体距離 $a=\frac{1}{5}\,\mathrm m$、焦点距離 $f=\frac{1}{10}\,\mathrm m$、物側の周辺光線角 $\theta_o=\frac{1}{50}\,\mathrm{rad}$、像面のずれ $\Delta b=\frac{1}{1000}\,\mathrm m$ である。

薄レンズ式を像距離について解く。

$$
\begin{aligned}
b&=\frac{a f}{a - f},\\
b/(\mathrm m)&= \frac{\left(\frac{1}{5}\right) \cdot \left(\frac{1}{10}\right)}{\left(\frac{1}{5}\right) - \left(\frac{1}{10}\right)}=\frac{1}{5}.
\end{aligned}
$$

負号は実像が倒立することを示す。

$$
\begin{aligned}
M&=- \frac{b}{a},\\
M&= - \frac{\left(\frac{1}{5}\right)}{\left(\frac{1}{5}\right)}=-1.
\end{aligned}
$$

外側が同じ媒質なので、換算角の倍率を角度にも使える。

$$
\begin{aligned}
\theta_i&=\frac{\theta_o}{M},\\
\theta_i/(\mathrm{rad})&= \frac{\left(\frac{1}{50}\right)}{\left(-1\right)}=- \frac{1}{50}.
\end{aligned}
$$

周辺光線の焦点からの横ずれを取る。

$$
\begin{aligned}
r_{\mathrm{blur}}&=\left|{\Delta b \theta_i}\right|,\\
r_{\mathrm{blur}}/(\mathrm m)&= \left|{\left(\frac{1}{1000}\right) \cdot \left(- \frac{1}{50}\right)}\right|\simeq 2\times 10^{-5}.
\end{aligned}
$$

同じ共役条件で物体を軸方向に動かした場合の倍率も評価する。

$$
\begin{aligned}
M_z&=M^{2},\\
M_z&= \left(-1\right)^{2}=1.
\end{aligned}
$$

この共役では像側の開口角の大きさが物側と等しい。
:::

:::example[縮小像でぼけを測る]
物体距離を増やし、物側で同じ角度の光線を通す。

入力は物体距離 $a=\frac{2}{5}\,\mathrm m$、焦点距離 $f=\frac{1}{10}\,\mathrm m$、物側の周辺光線角 $\theta_o=\frac{1}{50}\,\mathrm{rad}$、像面のずれ $\Delta b=\frac{1}{1000}\,\mathrm m$ である。

薄レンズ式を像距離について解く。

$$
\begin{aligned}
b&=\frac{a f}{a - f},\\
b/(\mathrm m)&= \frac{\left(\frac{2}{5}\right) \cdot \left(\frac{1}{10}\right)}{\left(\frac{2}{5}\right) - \left(\frac{1}{10}\right)}=\frac{2}{15}.
\end{aligned}
$$

負号は実像が倒立することを示す。

$$
\begin{aligned}
M&=- \frac{b}{a},\\
M&= - \frac{\left(\frac{2}{15}\right)}{\left(\frac{2}{5}\right)}=- \frac{1}{3}.
\end{aligned}
$$

外側が同じ媒質なので、換算角の倍率を角度にも使える。

$$
\begin{aligned}
\theta_i&=\frac{\theta_o}{M},\\
\theta_i/(\mathrm{rad})&= \frac{\left(\frac{1}{50}\right)}{\left(- \frac{1}{3}\right)}=- \frac{3}{50}.
\end{aligned}
$$

周辺光線の焦点からの横ずれを取る。

$$
\begin{aligned}
r_{\mathrm{blur}}&=\left|{\Delta b \theta_i}\right|,\\
r_{\mathrm{blur}}/(\mathrm m)&= \left|{\left(\frac{1}{1000}\right) \cdot \left(- \frac{3}{50}\right)}\right|\simeq 6\times 10^{-5}.
\end{aligned}
$$

同じ共役条件で物体を軸方向に動かした場合の倍率も評価する。

$$
\begin{aligned}
M_z&=M^{2},\\
M_z&= \left(- \frac{1}{3}\right)^{2}=\frac{1}{9}.
\end{aligned}
$$

像は縮小するが像側開口角は前例より大きい。像面の同じずれが大きな幾何学ぼけになる。
:::

:::example[拡大像で縦倍率を確認する]
物体を前焦点に近づけ、同じ物側角を使う。

入力は物体距離 $a=\frac{1}{8}\,\mathrm m$、焦点距離 $f=\frac{1}{10}\,\mathrm m$、物側の周辺光線角 $\theta_o=\frac{1}{50}\,\mathrm{rad}$、像面のずれ $\Delta b=\frac{1}{1000}\,\mathrm m$ である。

薄レンズ式を像距離について解く。

$$
\begin{aligned}
b&=\frac{a f}{a - f},\\
b/(\mathrm m)&= \frac{\left(\frac{1}{8}\right) \cdot \left(\frac{1}{10}\right)}{\left(\frac{1}{8}\right) - \left(\frac{1}{10}\right)}=\frac{1}{2}.
\end{aligned}
$$

負号は実像が倒立することを示す。

$$
\begin{aligned}
M&=- \frac{b}{a},\\
M&= - \frac{\left(\frac{1}{2}\right)}{\left(\frac{1}{8}\right)}=-4.
\end{aligned}
$$

外側が同じ媒質なので、換算角の倍率を角度にも使える。

$$
\begin{aligned}
\theta_i&=\frac{\theta_o}{M},\\
\theta_i/(\mathrm{rad})&= \frac{\left(\frac{1}{50}\right)}{\left(-4\right)}=- \frac{1}{200}.
\end{aligned}
$$

周辺光線の焦点からの横ずれを取る。

$$
\begin{aligned}
r_{\mathrm{blur}}&=\left|{\Delta b \theta_i}\right|,\\
r_{\mathrm{blur}}/(\mathrm m)&= \left|{\left(\frac{1}{1000}\right) \cdot \left(- \frac{1}{200}\right)}\right|\simeq 5\times 10^{-6}.
\end{aligned}
$$

同じ共役条件で物体を軸方向に動かした場合の倍率も評価する。

$$
\begin{aligned}
M_z&=M^{2},\\
M_z&= \left(-4\right)^{2}=16.
\end{aligned}
$$

像側角が小さくなる一方、物体を前後に動かすと像面は大きく動く。センサーの変位と物体の変位を区別する必要がある。
:::

## 球面鏡と反射系

右側から左側へ戻る反射光を扱うときは進行軸を反転して展開する。鏡の凹面に入射側から見た正半径 $\rho>0$ を割り当てる。入射軸の正方向を鏡へ向かう向きとする。鏡面法線のうち同じ側を向くものを選ぶと、法線角は近軸で $\alpha=y/\rho$ である。反射則から、実験室座標での反射光線の角は $\pi+2\alpha-\theta_{\mathrm{in}}$ となる。出射軸では $z$ の符号だけを反転するので、横座標 $y$ を保った出射角は $\theta_{\mathrm{out}}=\theta_{\mathrm{in}}-2y/\rho$ となる。従って $1/a+1/b=2/\rho$、$f=\rho/2$ である。ここでは実像の $b>0$ は反射光が進む側で測る。凸面では $\rho<0$ となる。

:::example[凹面鏡の曲率中心]
$a=\rho$ を代入すると $b=\rho$ で、倍率は $-1$ である。軸上の物点が曲率中心にあれば、法線入射なので各光線は厳密に同じ軌道を戻る。
:::
:::example[放物面の焦点]
軸上の焦点を $z=f$、鏡面を $z=y^2/(4f)$ とする。焦点までの距離は $\sqrt{(z-f)^2+y^2}=z+f$ である。上方の平面 $z=Z$ から軸に平行に入る光線の全長は $(Z-z)+(z+f)$ で一定である。したがって鏡の全開口で焦点が共通になる。
:::

### 球面鏡の開口を広げる
軸に平行な入射光を考え、入射側へ測った鏡面位置を $s(y)=\rho-\sqrt{\rho^2-y^2}$ とする。反射光が軸へ向かう角の大きさは $2\arcsin(y/\rho)$。軸との交点の位置は $s(y)+y\cot[2\arcsin(y/\rho)]$ である。三角関数を二倍角で展開すると、この位置は $\rho-\rho^2/(2\sqrt{\rho^2-y^2})=\rho/2-y^2/(4\rho)+O(y^4/\rho^3)$ となる。$|y|<\rho/\sqrt2$ の開口では反射光が入射側を向く。この範囲で $|y|$ を増やすと、軸交点の座標は頂点側へ移る。

:::example[縦球面収差]
$\rho=1\,\mathrm m$ で $y=0.10\,\mathrm m$ の光線は、二次近似で近軸焦点より $2.5\,\mathrm{mm}$ 頂点側で軸に交わる。軸上焦点を一つ選んでも幅のある平行束のすべては通らない。
:::

### 球面鏡と放物面鏡のたわみを比較する

頂点を原点、入射側の軸を正とする凹面鏡を考える。球面の半径を $\rho$ とすると、同じ近軸焦点 $\rho/2$ を持つ放物面のたわみは $y^2/(2\rho)$。球面では $s_{\rm sph}=\rho-\sqrt{\rho^2-y^2}$ で、二つの差は
$$
s_{\rm sph}-s_{\rm par}=\frac{y^4}{2\rho(\rho+\sqrt{\rho^2-y^2})^2}.
$$
これは平方根を有理化して二次項を引けば得られる。鏡面を軸方向に $\Delta s$ 動かしたときの同一基準点までの経路長は、往路と復路の寄与を足すと近軸で $-2\Delta s$ だけ変わる。従って球面と放物面の光路差は先頭で $-y^4/(4\rho^3)$。開口の角が小さくても、波長より大きい光路差を生じ得る。鏡の形状誤差を測ったとき、たわみと反射波の位相差を同じ値にしない。

球面の法線角は $\alpha=\arcsin(y/\rho)$ なので、反射軌道の軸への角は $2\alpha$。$y<\rho/\sqrt2$ では鏡の交点 $s_{\rm sph}$ より入射側の軸で交わり、頂点からの距離は $s_{\rm sph}+y/\tan(2\alpha)$ である。

:::example[曲率の緩い鏡の内側]
軸近くを使い、縦変位と波面差を比べる。

入力は凹面鏡の半径 $\rho=\frac{1}{2}\,\mathrm m$、入射高さ $y=\frac{1}{200}\,\mathrm m$、真空波長 $\lambda_0\simeq 5\times 10^{-7}\,\mathrm m$ である。

球面の入射点で法線の傾きを求める。

$$
\begin{aligned}
\alpha&=\operatorname{asin}{\left(\frac{y}{\rho} \right)},\\
\alpha/(\mathrm{rad})&= \operatorname{asin}{\left(\frac{\left(\frac{1}{200}\right)}{\left(\frac{1}{2}\right)} \right)}\simeq 0.01.
\end{aligned}
$$

軸交点までの距離を鏡の頂点から測るために、交点の座標も入れる。

$$
\begin{aligned}
s&=\rho - \sqrt{\rho^{2} - y^{2}},\\
s/(\mathrm m)&= \left(\frac{1}{2}\right) - \sqrt{\left(\frac{1}{2}\right)^{2} - \left(\frac{1}{200}\right)^{2}}\simeq 2.5001\times 10^{-5}.
\end{aligned}
$$

反射後は法線角の二倍だけ内側に向かう。

$$
\begin{aligned}
b&=s + \frac{y}{\tan{\left(2 \alpha \right)}},\\
b/(\mathrm m)&\simeq \left(2.500062503\times 10^{-5}\right) + \frac{\left(\frac{1}{200}\right)}{\tan{\left(2 \cdot \left(0.01000016667\right) \right)}}\simeq 0.24999.
\end{aligned}
$$

全光線で同じ近軸焦点からどれだけ離れるか調べる。

$$
\begin{aligned}
\Delta b&=- \frac{\rho}{2} + b,\\
\Delta b/(\mathrm m)&\simeq - \frac{\left(\frac{1}{2}\right)}{2} + \left(0.2499874991\right)\simeq -1.2501\times 10^{-5}.
\end{aligned}
$$

同じ近軸焦点を持つ放物面との差を波長で割る。

$$
\begin{aligned}
W_4/\lambda_0&=- \frac{y^{4}}{4 \rho^{3} \lambda_0},\\
W_4/\lambda_0&\simeq - \frac{\left(\frac{1}{200}\right)^{4}}{4 \cdot \left(\frac{1}{2}\right)^{3} \cdot \left(5\times 10^{-7}\right)}=- \frac{1}{400}.
\end{aligned}
$$

ここでは波面差は一波長より十分小さい。
:::

:::example[開口を四倍にする]
鏡はそのままで、照射される半径を四倍にする。

入力は凹面鏡の半径 $\rho=\frac{1}{2}\,\mathrm m$、入射高さ $y=\frac{1}{50}\,\mathrm m$、真空波長 $\lambda_0\simeq 5\times 10^{-7}\,\mathrm m$ である。

球面の入射点で法線の傾きを求める。

$$
\begin{aligned}
\alpha&=\operatorname{asin}{\left(\frac{y}{\rho} \right)},\\
\alpha/(\mathrm{rad})&= \operatorname{asin}{\left(\frac{\left(\frac{1}{50}\right)}{\left(\frac{1}{2}\right)} \right)}\simeq 0.040011.
\end{aligned}
$$

軸交点までの距離を鏡の頂点から測るために、交点の座標も入れる。

$$
\begin{aligned}
s&=\rho - \sqrt{\rho^{2} - y^{2}},\\
s/(\mathrm m)&= \left(\frac{1}{2}\right) - \sqrt{\left(\frac{1}{2}\right)^{2} - \left(\frac{1}{50}\right)^{2}}\simeq 4.0016\times 10^{-4}.
\end{aligned}
$$

反射後は法線角の二倍だけ内側に向かう。

$$
\begin{aligned}
b&=s + \frac{y}{\tan{\left(2 \alpha \right)}},\\
b/(\mathrm m)&\simeq \left(4.001601281\times 10^{-4}\right) + \frac{\left(\frac{1}{50}\right)}{\tan{\left(2 \cdot \left(0.04001067435\right) \right)}}\simeq 0.2498.
\end{aligned}
$$

全光線で同じ近軸焦点からどれだけ離れるか調べる。

$$
\begin{aligned}
\Delta b&=- \frac{\rho}{2} + b,\\
\Delta b/(\mathrm m)&\simeq - \frac{\left(\frac{1}{2}\right)}{2} + \left(0.2497997597\right)\simeq -2.0024\times 10^{-4}.
\end{aligned}
$$

同じ近軸焦点を持つ放物面との差を波長で割る。

$$
\begin{aligned}
W_4/\lambda_0&=- \frac{y^{4}}{4 \rho^{3} \lambda_0},\\
W_4/\lambda_0&\simeq - \frac{\left(\frac{1}{50}\right)^{4}}{4 \cdot \left(\frac{1}{2}\right)^{3} \cdot \left(5\times 10^{-7}\right)}=- \frac{16}{25}.
\end{aligned}
$$

近軸角はなお小さいが、位相誤差は半径の四乗で増える。
:::

:::example[小型鏡で同じ開口半径を使う]
より強く曲がる鏡に取り替え、波長と照射半径を固定する。

入力は凹面鏡の半径 $\rho=\frac{1}{5}\,\mathrm m$、入射高さ $y=\frac{1}{50}\,\mathrm m$、真空波長 $\lambda_0\simeq 5\times 10^{-7}\,\mathrm m$ である。

球面の入射点で法線の傾きを求める。

$$
\begin{aligned}
\alpha&=\operatorname{asin}{\left(\frac{y}{\rho} \right)},\\
\alpha/(\mathrm{rad})&= \operatorname{asin}{\left(\frac{\left(\frac{1}{50}\right)}{\left(\frac{1}{5}\right)} \right)}\simeq 0.10017.
\end{aligned}
$$

軸交点までの距離を鏡の頂点から測るために、交点の座標も入れる。

$$
\begin{aligned}
s&=\rho - \sqrt{\rho^{2} - y^{2}},\\
s/(\mathrm m)&= \left(\frac{1}{5}\right) - \sqrt{\left(\frac{1}{5}\right)^{2} - \left(\frac{1}{50}\right)^{2}}\simeq 0.0010025.
\end{aligned}
$$

反射後は法線角の二倍だけ内側に向かう。

$$
\begin{aligned}
b&=s + \frac{y}{\tan{\left(2 \alpha \right)}},\\
b/(\mathrm m)&\simeq \left(0.001002512579\right) + \frac{\left(\frac{1}{50}\right)}{\tan{\left(2 \cdot \left(0.1001674212\right) \right)}}\simeq 0.099496.
\end{aligned}
$$

全光線で同じ近軸焦点からどれだけ離れるか調べる。

$$
\begin{aligned}
\Delta b&=- \frac{\rho}{2} + b,\\
\Delta b/(\mathrm m)&\simeq - \frac{\left(\frac{1}{5}\right)}{2} + \left(0.09949621847\right)\simeq -5.0378\times 10^{-4}.
\end{aligned}
$$

同じ近軸焦点を持つ放物面との差を波長で割る。

$$
\begin{aligned}
W_4/\lambda_0&=- \frac{y^{4}}{4 \rho^{3} \lambda_0},\\
W_4/\lambda_0&\simeq - \frac{\left(\frac{1}{50}\right)^{4}}{4 \cdot \left(\frac{1}{5}\right)^{3} \cdot \left(5\times 10^{-7}\right)}=-10.
\end{aligned}
$$

同じ半径の照射でも小さい曲率半径が波面差を増す。角だけで回折像の質を判定できない。
:::

## 行列光学（ABCD 行列）

### 行列の順序と結像
一様媒質の距離 $d$ の伝播で $y\mapsto y+dp/n$ となり、屈折では $p\mapsto p-\Phi y$ となる。よって :term[光線伝達行列]{en="ray transfer matrix" yomi="こうせんでんたつぎょうれつ"} は
$$
\binom{y_2}{p_2}=M\binom{y_1}{p_1},\qquad
T(d/n)=\begin{pmatrix}1&d/n\\0&1\end{pmatrix},\quad
L(\Phi)=\begin{pmatrix}1&0\\-\Phi&1\end{pmatrix}.
\label{eq-optics-ray-matrix}
$$
先に通る面の行列を右に置く。各行列の行列式が 1 なので積も 1 である。$M=\begin{pmatrix}A&B\\C&D\end{pmatrix}$ で $B=0$ のときだけ、初期値 $p_1$ によらず $y_2=Ay_1$ となる。従ってこの条件が結像を表し、倍率は $A$ である。二本の光線を列に並べた行列の行列式 $y_1p_2-y_2p_1$ も保存される。

:::example[望遠系の間隔]
空気中で凸レンズ二枚を距離 $d$ だけ離す。$L(1/f_2)T(d)L(1/f_1)$ の左下は $-1/f_1-1/f_2+d/(f_1 f_2)$ である。$d=f_1+f_2$ ならゼロとなり、平行入力は平行出力に変わる。
:::
:::example[結像の位置を行列で選ぶ]
レンズの手前 $a$、後ろ $b$ の面を取ると $T(b)L(1/f)T(a)$ の右上は $a+b-ab/f$ である。これをゼロと置けば :eqref[eq-optics-thin-lens] になる。
:::

### 並べたレンズが作る屈折力
同じ媒質 $n$ の薄レンズ二つを間隔 $d$ で置く。:eqref[eq-optics-ray-matrix] を掛けると合成屈折力は $\Phi_1+\Phi_2-(d/n)\Phi_1\Phi_2$。並び順を逆にした場合も合成屈折力は同じだが、$A=1-(d/n)\Phi_1$ と $D=1-(d/n)\Phi_2$ が交換するので、頂点から焦点までの距離は一般に変わる。

:::example[焦点の違い]
空気中で $\Phi_1=4\,\mathrm{m^{-1}},\Phi_2=6\,\mathrm{m^{-1}},d=0.10\,\mathrm m$ なら $\Phi=7.6\,\mathrm{m^{-1}}$。出力で高さが $Ay$、角が $-\Phi y$ となる平行入射を追うと、焦点までの距離は順方向で $0.6/7.6\,\mathrm m$、レンズ順を逆にすると $0.4/7.6\,\mathrm m$ である。
:::

### 望遠系を行列で調整する

外側が空気の二薄レンズを間隔 $d$ で並べる。直前から直後までの写像は $L(1/f_2)T(d)L(1/f_1)$ である。積を取ると $A=1-d/f_1$、$B=d$、$C=(d-f_1-f_2)/(f_1 f_2)$、$D=1-d/f_2$。$d=f_1+f_2$ で $C=0$ となり、軸に平行な入力から再び平行光が出る。このとき横幅の倍率は $A=-f_2/f_1$、入射角の倍率は $D=-f_1/f_2$ で、積が $1$ になる。入力に平行な光線の束は、間隔が $\delta d$ ずれると出力角に $C y=\delta d\,y/(f_1 f_2)$ が加わる。調整では、出力角が入射高さに対して傾かないように $d$ を選ぶ。

結像位置を求める場合には手順が異なる。平行入射に対し出力直後の高さが $A y$、角が $C y$ なので、さらに距離 $-A/C$ 進むと軸へ達する。$C=0$ の近くではこの距離の絶対値が大きくなり、有限の焦点距離を測る方法は間隔誤差に敏感になる。

:::example[正の二レンズを無焦点に合わせる]
両レンズに正の焦点距離を取り、平行光の径を変える。

入力は第一レンズの焦点距離 $f_{1}=\frac{1}{5}\,\mathrm m$、第二レンズの焦点距離 $f_{2}=\frac{1}{20}\,\mathrm m$、レンズ間隔 $d=\frac{1}{4}\,\mathrm m$、入射光線の高さ $y=\frac{1}{250}\,\mathrm m$ である。

第一レンズで角を変え、第二レンズまで進める。

$$
\begin{aligned}
A&=- \frac{d}{f_{1}} + 1,\\
A&= - \frac{\left(\frac{1}{4}\right)}{\left(\frac{1}{5}\right)} + 1=- \frac{1}{4}.
\end{aligned}
$$

第二レンズで再び角を変える。

$$
\begin{aligned}
C&=\frac{d - f_{1} - f_{2}}{f_{1} f_{2}},\\
C/(\mathrm{m^{-1}})&= \frac{\left(\frac{1}{4}\right) - \left(\frac{1}{5}\right) - \left(\frac{1}{20}\right)}{\left(\frac{1}{5}\right) \cdot \left(\frac{1}{20}\right)}=0.
\end{aligned}
$$

高さゼロで斜めに入れた光線の出力角を求める。

$$
\begin{aligned}
D&=- \frac{d}{f_{2}} + 1,\\
D&= - \frac{\left(\frac{1}{4}\right)}{\left(\frac{1}{20}\right)} + 1=-4.
\end{aligned}
$$

入射角をゼロとする。

$$
\begin{aligned}
y'&=A y,\\
y'/(\mathrm m)&= \left(- \frac{1}{4}\right) \cdot \left(\frac{1}{250}\right)=- \frac{1}{1000}.
\end{aligned}
$$

余分な発散または収束を角で表す。

$$
\begin{aligned}
\theta'&=C y,\\
\theta'/(\mathrm{rad})&= \left(0\right) \cdot \left(\frac{1}{250}\right)=0.
\end{aligned}
$$

出力は幅が狭くなる。入力全体を傾けた場合は角度が逆向きに拡大する。
:::

:::example[間隔をわずかに増やす]
組立誤差が出力の平行度にどの程度効くか調べる。

入力は第一レンズの焦点距離 $f_{1}=\frac{1}{5}\,\mathrm m$、第二レンズの焦点距離 $f_{2}=\frac{1}{20}\,\mathrm m$、レンズ間隔 $d=\frac{251}{1000}\,\mathrm m$、入射光線の高さ $y=\frac{1}{250}\,\mathrm m$ である。

第一レンズで角を変え、第二レンズまで進める。

$$
\begin{aligned}
A&=- \frac{d}{f_{1}} + 1,\\
A&= - \frac{\left(\frac{251}{1000}\right)}{\left(\frac{1}{5}\right)} + 1=- \frac{51}{200}.
\end{aligned}
$$

第二レンズで再び角を変える。

$$
\begin{aligned}
C&=\frac{d - f_{1} - f_{2}}{f_{1} f_{2}},\\
C/(\mathrm{m^{-1}})&= \frac{\left(\frac{251}{1000}\right) - \left(\frac{1}{5}\right) - \left(\frac{1}{20}\right)}{\left(\frac{1}{5}\right) \cdot \left(\frac{1}{20}\right)}=\frac{1}{10}.
\end{aligned}
$$

高さゼロで斜めに入れた光線の出力角を求める。

$$
\begin{aligned}
D&=- \frac{d}{f_{2}} + 1,\\
D&= - \frac{\left(\frac{251}{1000}\right)}{\left(\frac{1}{20}\right)} + 1=- \frac{201}{50}.
\end{aligned}
$$

入射角をゼロとする。

$$
\begin{aligned}
y'&=A y,\\
y'/(\mathrm m)&= \left(- \frac{51}{200}\right) \cdot \left(\frac{1}{250}\right)\simeq -0.00102.
\end{aligned}
$$

余分な発散または収束を角で表す。

$$
\begin{aligned}
\theta'&=C y,\\
\theta'/(\mathrm{rad})&= \left(\frac{1}{10}\right) \cdot \left(\frac{1}{250}\right)=\frac{1}{2500}.
\end{aligned}
$$

中心光線が平行でも、周辺光線には角が出る。無焦点の評価には高さの異なる光線が要る。
:::

:::example[負レンズを第二面に使う]
二枚の間に実焦点を置かない長さで平行光を縮める。

入力は第一レンズの焦点距離 $f_{1}=\frac{1}{5}\,\mathrm m$、第二レンズの焦点距離 $f_{2}=- \frac{1}{20}\,\mathrm m$、レンズ間隔 $d=\frac{3}{20}\,\mathrm m$、入射光線の高さ $y=\frac{1}{250}\,\mathrm m$ である。

第一レンズで角を変え、第二レンズまで進める。

$$
\begin{aligned}
A&=- \frac{d}{f_{1}} + 1,\\
A&= - \frac{\left(\frac{3}{20}\right)}{\left(\frac{1}{5}\right)} + 1=\frac{1}{4}.
\end{aligned}
$$

第二レンズで再び角を変える。

$$
\begin{aligned}
C&=\frac{d - f_{1} - f_{2}}{f_{1} f_{2}},\\
C/(\mathrm{m^{-1}})&= \frac{\left(\frac{3}{20}\right) - \left(\frac{1}{5}\right) - \left(- \frac{1}{20}\right)}{\left(\frac{1}{5}\right) \cdot \left(- \frac{1}{20}\right)}=0.
\end{aligned}
$$

高さゼロで斜めに入れた光線の出力角を求める。

$$
\begin{aligned}
D&=- \frac{d}{f_{2}} + 1,\\
D&= - \frac{\left(\frac{3}{20}\right)}{\left(- \frac{1}{20}\right)} + 1=4.
\end{aligned}
$$

入射角をゼロとする。

$$
\begin{aligned}
y'&=A y,\\
y'/(\mathrm m)&= \left(\frac{1}{4}\right) \cdot \left(\frac{1}{250}\right)=\frac{1}{1000}.
\end{aligned}
$$

余分な発散または収束を角で表す。

$$
\begin{aligned}
\theta'&=C y,\\
\theta'/(\mathrm{rad})&= \left(0\right) \cdot \left(\frac{1}{250}\right)=0.
\end{aligned}
$$

同じ倍率の大きさでも $A,D$ は正であり、出力の横座標と角の向きが保存する。
:::

## 屈折率分布媒質と蜃気楼

### 水平に一様な媒質と転回点
$n=n(y)$ のとき :eqref[eq-optics-ray] から $n\dd x/\dd s=P$ は一定である。$P\ne0$ の光線では、単位接線の条件より
$$
(\dd y/\dd x)^2=n(y)^2/P^2-1 .
$$
実光線の領域は $n(y)\ge|P|$ で、等号では水平になる。$P=0$ の光線では $\dd x/\dd s=0$ なので、軌跡を $y(x)$ として表さず鉛直線として追う。$n'(y)\ne0$ なら $\dd^2y/\dd s^2=n'(y)/n(y)$ がその点で成り立つので、光線は向きを変えて許容域へ戻る。熱した地表近くで $n$ が低く、上ほど高い場合、斜め下へ進んできた光は上へ戻る。観察者が直線を仮定して光を逆延長すると、地表近くに空の像が見える。

### 放物型分布
$n(y)=n_0(1-g^2 y^2/2)$、$|gy|\ll1$ とする。近軸で $\dd s\simeq\dd z$ とし、$p\simeq n_0\dd y/\dd z$ を :eqref[eq-optics-ray] に入れると $y''=-g^2y$ である。初期値 $y_0,p_0$ に対し
$$
y(z)=y_0\cos gz+\frac{p_0}{n_0g}\sin gz,\qquad
p(z)=-n_0g y_0\sin gz+p_0\cos gz .
$$
$g$ の単位は $\mathrm{m^{-1}}$ である。

:::example[四分の一周期の板]
$gz=\pi/2$ とすると $p_0=0$ の光線は任意の $y_0$ から出力面の軸上へ集まる。面が平坦でも、内部の勾配がレンズ作用を持つ。
:::
:::example[許容域を選ぶ]
水平成分 $P=1.00025$ の光線は $n=1.00020$ の層へ入れない。上から降りてくる途中で $n=P$ となる高さがあれば、そこで下向きの勾配が反転する。
:::

:::example[正弦的な光線の振幅]
$n_0=1.5,g=100\,\mathrm{m^{-1}}$ で $y_0=0$、$p_0=0.015$ なら、最大の高さは $p_0/(n_0 g)=0.10\,\mathrm{mm}$。$g|y|\le0.01$ なので、屈折率の放物展開に用いた小さい高さの条件を保つ。
:::

### 一次勾配の層で軌跡を積分する
$n(y)=n_b+ay$、$a>0$、水平成分 $P>0$ を考える。上向き枝で $\dd x/\dd y=P/\sqrt{(n_b+ay)^2-P^2}$。転回点の水平座標を $x_t$ とし、$u=(n_b+ay)/P$ で積分すると $x-x_t=(P/a)\operatorname{arcosh}u$。下向き枝も合わせ、光線は
$$
y(x)=\frac{P}{a}\cosh\left(\frac{a(x-x_t)}P\right)-\frac{n_b}{a}.
$$
転回点で傾きはゼロ、曲率の二階微分は $a/P$ である。一次勾配が実際に保たれる高さの範囲でのみこの式を使える。勾配が層を越えて変わるなら、その位置で接線成分をそろえて別の軌跡へつなぐ。

### 勾配媒質の位相空間と光線の周期

$n(y)=n_0[1-g^2 y^2/2+O(y^4)]$ とし、軸方向 $z$ を独立変数に取る。光路積分の被積分関数は $n(y)\sqrt{1+(y')^2}=n_0+n_0((y')^2-g^2y^2)/2+O((y,y')^4)$ である。Euler–Lagrange 方程式は $y''+g^2y=0$。初期条件 $y(0)=y_0$、$y'(0)=\theta_0$ なら $y(z)=y_0\cos(gz)+(\theta_0/g)\sin(gz)$。微分して $y'$ を求めれば、$y^2+(y'/g)^2=y_0^2+(\theta_0/g)^2$ が三角関数の平方和によって保存する。軸を横切る角と最高到達点を一つの量で比較できる。

厚さを $gz=\pi/2$ にすれば高さと角の役割が入れ替わり、異なる高さから光軸と平行に入る光線が、出力面の軸上に集まる。$gz=\pi$ なら $(y,y')$ の符号が両方反転する。式を用いる範囲は $g|y|\ll1$ かつ $|y'|\ll1$ の範囲であり、実際の媒質の外縁が振幅より内側なら光線は途中で外へ出る。

:::example[平行光を軸へ寄せる]
高さを持つ光線を平行に入れ、四分の一周期の厚さで取り出す。

入力は勾配定数 $g=100\,\mathrm{m^{-1}}$、伝播距離 $z\simeq 0.015708\,\mathrm m$、初期高さ $y_{0}=\frac{1}{5000}\,\mathrm m$、初期角 $\theta_0=0\,\mathrm{rad}$ である。

勾配媒質での振動の進み具合を求める。

$$
\begin{aligned}
gz&=g z,\\
gz/(\mathrm{rad})&\simeq \left(100\right) \cdot \left(0.01570796327\right)\simeq 1.5708.
\end{aligned}
$$

初期高さと換算した角を直交成分として合成する。

$$
\begin{aligned}
y_{\max}&=\sqrt{y_{0}^{2} + \frac{\theta_0^{2}}{g^{2}}},\\
y_{\max}/(\mathrm m)&= \sqrt{\left(\frac{1}{5000}\right)^{2} + \frac{\left(0\right)^{2}}{\left(100\right)^{2}}}=\frac{1}{5000}.
\end{aligned}
$$

二つの初期条件が余弦項と正弦項にそれぞれ効く。

$$
\begin{aligned}
y&=y_{0} \cos{\left(gz \right)} + \frac{\theta_0 \sin{\left(gz \right)}}{g},\\
y/(\mathrm m)&\simeq \left(\frac{1}{5000}\right) \cdot \cos{\left(\left(1.570796327\right) \right)} + \frac{\left(0\right) \cdot \sin{\left(\left(1.570796327\right) \right)}}{\left(100\right)}=0.
\end{aligned}
$$

高さを $z$ で微分する。

$$
\begin{aligned}
\theta&=- g y_{0} \sin{\left(gz \right)} + \theta_0 \cos{\left(gz \right)},\\
\theta/(\mathrm{rad})&\simeq - \left(100\right) \cdot \left(\frac{1}{5000}\right) \cdot \sin{\left(\left(1.570796327\right) \right)} + \left(0\right) \cdot \cos{\left(\left(1.570796327\right) \right)}=- \frac{1}{50}.
\end{aligned}
$$

光線が軸に達するとき、初期条件で決まる最大角を取る。

$$
\begin{aligned}
\theta_{\max}&=y_{\max} g,\\
\theta_{\max}/(\mathrm{rad})&= \left(\frac{1}{5000}\right) \cdot \left(100\right)=\frac{1}{50}.
\end{aligned}
$$

高さはゼロ、傾きは負になり、軸を横切る。
:::

:::example[角度を位置へ写す]
軸上の点から斜めに出た光線を前例と同じ板に入れる。

入力は勾配定数 $g=100\,\mathrm{m^{-1}}$、伝播距離 $z=\pi/(2g)\simeq 0.015708\,\mathrm m$、初期高さ $y_{0}=0\,\mathrm m$、初期角 $\theta_0=\frac{1}{50}\,\mathrm{rad}$ である。

勾配媒質での振動の進み具合を求める。

$$
\begin{aligned}
gz&=g z,\\
gz/(\mathrm{rad})&\simeq \left(100\right) \cdot \left(0.01570796327\right)\simeq 1.5708.
\end{aligned}
$$

初期高さと換算した角を直交成分として合成する。

$$
\begin{aligned}
y_{\max}&=\sqrt{y_{0}^{2} + \frac{\theta_0^{2}}{g^{2}}},\\
y_{\max}/(\mathrm m)&= \sqrt{\left(0\right)^{2} + \frac{\left(\frac{1}{50}\right)^{2}}{\left(100\right)^{2}}}=\frac{1}{5000}.
\end{aligned}
$$

二つの初期条件が余弦項と正弦項にそれぞれ効く。

$$
\begin{aligned}
y&=y_{0} \cos{\left(gz \right)} + \frac{\theta_0 \sin{\left(gz \right)}}{g},\\
y/(\mathrm m)&\simeq \left(0\right) \cdot \cos{\left(\pi/2 \right)} + \frac{\left(\frac{1}{50}\right) \cdot \sin{\left(\pi/2 \right)}}{\left(100\right)}=\frac{1}{5000}.
\end{aligned}
$$

高さを $z$ で微分する。

$$
\begin{aligned}
\theta&=- g y_{0} \sin{\left(gz \right)} + \theta_0 \cos{\left(gz \right)},\\
\theta/(\mathrm{rad})&\simeq - \left(100\right) \cdot \left(0\right) \cdot \sin{\left(\pi/2 \right)} + \left(\frac{1}{50}\right) \cdot \cos{\left(\pi/2 \right)}=0.
\end{aligned}
$$

最大の傾きは軸上で生じる。この例では入射面が軸上の点を通り、出力面では高さが最大になる。

$$
\begin{aligned}
\theta_{\max}&=y_{\max} g,\\
\theta_{\max}/(\mathrm{rad})&= \left(\frac{1}{5000}\right) \cdot \left(100\right)=\frac{1}{50}.
\end{aligned}
$$

出力高さが入射角に比例し、出力角はゼロになる。
:::

:::example[半周期の板を使う]
高さと角を両方与え、二倍の厚さへ進める。

入力は勾配定数 $g=100\,\mathrm{m^{-1}}$、伝播距離 $z\simeq 0.031416\,\mathrm m$、初期高さ $y_{0}=\frac{1}{5000}\,\mathrm m$、初期角 $\theta_0=\frac{1}{100}\,\mathrm{rad}$ である。

勾配媒質での振動の進み具合を求める。

$$
\begin{aligned}
gz&=g z,\\
gz/(\mathrm{rad})&\simeq \left(100\right) \cdot \left(0.03141592654\right)\simeq 3.1416.
\end{aligned}
$$

初期高さと換算した角を直交成分として合成する。

$$
\begin{aligned}
y_{\max}&=\sqrt{y_{0}^{2} + \frac{\theta_0^{2}}{g^{2}}},\\
y_{\max}/(\mathrm m)&= \sqrt{\left(\frac{1}{5000}\right)^{2} + \frac{\left(\frac{1}{100}\right)^{2}}{\left(100\right)^{2}}}\simeq 2.2361\times 10^{-4}.
\end{aligned}
$$

二つの初期条件が余弦項と正弦項にそれぞれ効く。

$$
\begin{aligned}
y&=y_{0} \cos{\left(gz \right)} + \frac{\theta_0 \sin{\left(gz \right)}}{g},\\
y/(\mathrm m)&\simeq \left(\frac{1}{5000}\right) \cdot \cos{\left(\left(3.141592654\right) \right)} + \frac{\left(\frac{1}{100}\right) \cdot \sin{\left(\left(3.141592654\right) \right)}}{\left(100\right)}=- \frac{1}{5000}.
\end{aligned}
$$

高さを $z$ で微分する。

$$
\begin{aligned}
\theta&=- g y_{0} \sin{\left(gz \right)} + \theta_0 \cos{\left(gz \right)},\\
\theta/(\mathrm{rad})&\simeq - \left(100\right) \cdot \left(\frac{1}{5000}\right) \cdot \sin{\left(\left(3.141592654\right) \right)} + \left(\frac{1}{100}\right) \cdot \cos{\left(\left(3.141592654\right) \right)}=- \frac{1}{100}.
\end{aligned}
$$

光線が軸に達するとき、初期条件で決まる最大角を取る。

$$
\begin{aligned}
\theta_{\max}&=y_{\max} g,\\
\theta_{\max}/(\mathrm{rad})&\simeq \left(2.236067977\times 10^{-4}\right) \cdot \left(100\right)\simeq 0.022361.
\end{aligned}
$$

出力で両方の符号が反転する。途中の最大高さが媒質の半径より小さければ、この像の関係を利用できる。
:::

## 演習問題

::::exercise{level=basic}
一定の屈折率で光線が直線であることを変分式から示せ。

:::solution
勾配 $\grad n=0$ なので :eqref[eq-optics-optical-variation] から $\dd(n\hat{\vb t})/\dd s=0$。$n>0$ は一定であり $\hat{\vb t}$ も一定となる。
:::
::::

::::exercise{level=basic}
$n_1=1,n_2=1.5,R=0.25\,\mathrm m$ の面の $\Phi$ を求めよ。

:::solution
屈折力は $(1.5-1)/0.25=2\,\mathrm{m^{-1}}$ である。
:::
::::

::::exercise{level=basic}
$f=20\,\mathrm{cm}$ の薄レンズに $a=40\,\mathrm{cm}$ で物を置く。像位置と倍率を求めよ。

:::solution
$1/b=1/20-1/40=1/40\,\mathrm{cm^{-1}}$ なので $b=40\,\mathrm{cm}$、$M=-1$ である。
:::
::::

::::exercise{level=basic}
$\rho=60\,\mathrm{cm}$ の凹面鏡に平行光を当てる。像位置を求めよ。

:::solution
入射光線の傾きがゼロなので、反射後は $-2y/\rho$ の傾きとなる。像位置は $y/(2y/\rho)=\rho/2=30\,\mathrm{cm}$ である。
:::
::::

::::exercise{level=advanced}
外部媒質を空気とみなし、$T(b)L(\Phi)T(a)$ が結像するとき、この行列の $D$ が横倍率の逆数になることを示せ。

:::solution
結像では $B=0$。三行列とも行列式が 1 なので $AD-BC=AD=1$ となる。$A$ が倍率であるから $D=1/A$ である。
:::
::::

::::exercise{level=advanced}
同じ媒質中で $f_1=10\,\mathrm{cm}$、$f_2=20\,\mathrm{cm}$ のレンズの間隔を $30\,\mathrm{cm}$ とする。角倍率を求めよ。

:::solution
出力の角は $C=0$ により初期高さによらない。積の右下は $D=1-d/f_2=1-30/20=-1/2$。同じ媒質では $p$ と角が同じ比例係数を持つので角倍率は $-1/2$ である。
:::
::::

::::exercise{level=advanced}
放物型分布の周期の半分で得られる光線行列を求めよ。

:::solution
$gz=\pi$ では $\sin gz=0,\cos gz=-1$。従って $y=-y_0,p=-p_0$ で、行列は $-I_2$ である。
:::
::::

::::exercise{level=advanced}
$n(y)=n_b+ay$、$a>0$ の層内で、正の水平成分 $P>n_b$ を持つ光線が低屈折率側へ向かう。層内に $n=P$ の高さがあるとき、転回点を求めよ。

:::solution
水平となる点では $n=P$ なので $y_t=(P-n_b)/a>0$ である。許容域は $y\ge y_t$ で、勾配が正なので $y_t$ で上向きに転回する。
:::
::::

::::exercise{level=advanced}
同じ媒質に焦点距離 $f_1,f_2>0$ の二薄レンズを置く。間隔を $f_1+f_2+\delta$、ずれを $|\delta|<f_2$ とする。高さ $y$ で軸に平行に入れた光線の出力角を近軸で求め、$\delta$ の符号によって収束と発散を判定せよ。

:::solution
第一レンズ直後の角は $-y/f_1$、第二レンズ面の高さは $y-(f_1+f_2+\delta)y/f_1=-(f_2+\delta)y/f_1$。第二レンズでその高さを $f_2$ で割った量を角から引くので $\theta'=\delta y/(f_1f_2)$。$|\delta|<f_2$ なら出力高さは入射高さと逆符号である。そのため $\delta>0$ では出力の角が軸へ向かい、$\delta<0$ では軸から離れる。$\delta=0$ なら高さによらず平行に出る。
:::
::::

## 確認テスト

::::quiz
結像条件はどれか。

- [ ] $A=0$
- [x] $B=0$
- [ ] $D=0$

:::explanation
出力高さ $y_2=Ay_1+Bp_1$ を初期角に依存させないため $B=0$ とする。このとき倍率は $A$ である。
:::
::::

::::quiz{answer="30" tol="1e-6"}
空気中で $f=10\,\mathrm{cm}$、$a=15\,\mathrm{cm}$ のとき、$b$ を cm 単位で答えよ。

:::explanation
$1/b=1/10-1/15=1/30\,\mathrm{cm^{-1}}$ なので $b=30\,\mathrm{cm}$ である。
:::
::::

::::quiz
水平成層で保存される量を選べ。

- [x] $n\dd x/\dd s$
- [ ] $n\dd y/\dd s$
- [ ] $n$

:::explanation
水平屈折率勾配がゼロなので、光線方程式の水平成分を積分して $n\dd x/\dd s$ が一定になる。
:::
::::

::::quiz{answer="10" tol="1e-6"}
空気中で距離 $0.50\,\mathrm m$ を傾き $0.02$ の光線が進む。高さの増加を mm 単位で答えよ。

:::explanation
近軸では高さの増加が $d\theta$ である。$0.50(0.02)=0.010\,\mathrm m=10\,\mathrm{mm}$ となる。
:::
::::

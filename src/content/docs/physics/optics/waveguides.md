---
title: 導波光学と光ファイバー
description: 誘電体導波路の境界条件からモードを導き、ファイバーの分散・損失、結合器および周期構造の禁制帯を計算する。
notation:
  - symbol: '\beta'
    meaning: 軸方向の伝播定数
  - symbol: 'n_1,n_2'
    meaning: コアとクラッドの屈折率
  - symbol: 'u,w,V'
    meaning: 横波数と減衰定数およびそれらから決まる規格化周波数
  - symbol: '\kappa'
    meaning: モード結合係数
  - symbol: '\Lambda'
    meaning: 周期構造の一周期
---

## この章の目標

- 平板・円筒導波路でコアとクラッドの場を接続する。
- ファイバーの開口数、規格化周波数、単一モード条件を計算する。
- スペクトル・モード・吸収に由来する伝播の制約を区別する。
- 結合モードと Bloch 波で電力移行と禁制帯を説明する。

## 平板誘電体導波路のモード

### 屈折率の高い層へ場を閉じ込める
コア $|x|<a$ の屈折率を $n_1$、両クラッドを $n_2<n_1$ とする。非磁性、無損失、$y$ 一様で、$+z$ に進むモードを選び、TE モードの電場を $\hat{\vb y}E(x)e^{i\beta z-i\omega t}$ と置く。式 :eqref[eq-optics-em-wave] より各層で $E''+(k_0^2n^2-\beta^2)E=0$。クラッドの場を減衰させ、コアで振動させる条件は $k_0n_2<\beta<k_0 n_1$ である。$u=a\sqrt{k_0^2n_1^2-\beta^2}$、$w=a\sqrt{\beta^2-k_0^2n_2^2}$ と置くと $u^2+w^2=V^2$、$V=k_0 a\sqrt{n_1^2-n_2^2}$。ここで $2a$ が全幅である。

上下対称なので偶は中心で $\cos(ux/a)$、奇は $\sin(ux/a)$ になる。上クラッドで場は $E(a)e^{-w(x-a)/a}$。接線の $E_y,H_z$ が連続で、$H_z=E'/(i\omega\mu_0)$ なので、対数微分の一致によって
$$
u\tan u=w\quad\text{（偶 TE）},\qquad -u\cot u=w\quad\text{（奇 TE）}.
\label{eq-optics-slab-modes}
$$
TM では $H_y$ を未知関数とし、$H_y$ および $\eps^{-1}H_y'$ が連続である。従って二式の右辺が $(n_1^2/n_2^2)w$ になる。偶と奇の右辺は $u$ とともに減少する。左辺は正になる区間で増加し、最初の偶区間では $u=0$ の左辺が 0 で右辺が $V>0$。端で右辺がゼロになるか左辺が発散するので、基本偶モードが必ず一つある。最初の奇モードは $u=\pi/2,w=0$ を境に現れる。カットオフの $w=0$ 自体はクラッドに閉じ込められない。

:::example[開閉の判定]
$V=1$ なら $\pi/2$ に達しないので TE の導波解は最初の偶だけである。$V=2$ では奇も一つ存在する。
:::
:::example[減衰とモード場]
二つのモードで $w=0.1$ と $w=1$ を比べると、クラッドの振幅 $1/e$ 長は $10a$ と $a$ である。カットオフ近くではコアの幅だけでモードサイズを見積もれない。
:::

### 厚さとモード数
対称な平板の TE では、偶モードの左辺は $m\pi<u<m\pi+\pi/2$（$m\ge0$）、奇モードは $(m+1/2)\pi<u<(m+1)\pi$ で正である。各区間で左辺はゼロから増え、右辺 $\sqrt{V^2-u^2}$ は下がる。区間開始が $V$ 未満なら交点はちょうど一つあり、そうでなければ存在しない。TM でも右辺を正の定数倍するだけなので、本数とカットオフは TE と同じである。ただし伝播定数は一般に異なる。従って $V>0$ での束縛 TE モードの個数は $\lceil2V/\pi\rceil$。$V$ がカットオフと等しい場合には、開始点で $w=0$ の解を数えない。この数は $a$ が半幅である約束に依存する。

### 固有モードの直交性
TE の二解を $E_j,E_l$、それぞれの伝播定数を $\beta_j,\beta_l$ とする。片方の微分方程式に他方の共役を掛けて引き、実軸で積分する。二階微分項は境界項 $E_l^*E_j'-E_l'^*E_j$ であり、クラッドの無限遠で消え、各界面では連続性によって相殺する。従って $(\beta_j^2-\beta_l^2)\int E_l^*E_j\dd x=0$。固有値が異なるモードはこの積分で直交する。同じ導波路の横断面で前進電場を与えれば、各導波モードへの展開係数がこの内積で得られる。入口に屈折率の不連続がある場合は、反射と透過を含む境界値問題で、導波路内の前進場を先に求める必要がある。同じ手順で TM では方程式を $(H_y'/\eps)'$ を含む形に書き、重みに $\eps^{-1}$ が入る。導波モードに射影して残った場には、クラッドへ広がる放射成分も含まれる。

:::example[一つの根を寸法に戻す]
偶 TE で $u=\pi/4$ のとき $\tan u=1$ で $w=\pi/4$。従って $V=\pi/(2\sqrt2)$。半幅 $a=1\,\mathrm{\mu m}$ であればクラッドの侵入長は $a/w=4/\pi\,\mathrm{\mu m}\simeq1.27\,\mathrm{\mu m}$。コア半幅より大きい値であり、開口の幾何学的な幅だけで場の広さを表せない。
:::

### モードのパワーをコアとクラッドに配分する

偶 TE のコア電場を $E_0\cos(ux/a)$ とする。磁場の横成分が Faraday 則から $H_x=-\beta E_y/(\omega\mu_0)$ なので、軸方向の平均流は $\beta|E_y|^2/(2\omega\mu_0)$。従ってパワーの比を電場の二乗積分から求められる。コアは $|E_0|^2a[1+\sin(2u)/(2u)]$、クラッドの両側は $|E_0|^2a\cos^2u/w$。後者は境界から指数減衰する場を二つの半直線で積分した結果である。比較する場合は中心振幅の規格化を両領域で同じにする。$w$ が小さくなるほど、モードを規格化したときにコアで使える場の割合が減少する。

:::example[ごく弱い閉じ込め]
小さな正の $u$ から基本枝を作る。

入力は基本偶 TE の横波数の換算値 $u\simeq 0.2618$、コア半厚 $a\simeq 1\times 10^{-6}\,\mathrm m$ である。

境界の対数微分を合わせて外側の減衰を決める。

$$
\begin{aligned}
w&=u \tan{\left(u \right)},\\
w&\simeq \left(0.2617993878\right) \cdot \tan{\left(\left(0.2617993878\right) \right)}\simeq 0.070149.
\end{aligned}
$$

境界値を満たす $V$ を幾何学的に構成する。

$$
\begin{aligned}
V&=\sqrt{u^{2} + w^{2}},\\
V&\simeq \sqrt{\left(0.2617993878\right)^{2} + \left(0.07014893454\right)^{2}}\simeq 0.27103.
\end{aligned}
$$

指数の係数を実際の深さへ換算する。

$$
\begin{aligned}
\ell&=\frac{a}{w},\\
\ell/(\mathrm m)&\simeq \frac{\left(1\times 10^{-6}\right)}{\left(0.07014893454\right)}\simeq 1.4255\times 10^{-5}.
\end{aligned}
$$

中心振幅二乗と半幅を除いた積分を求める。

$$
\begin{aligned}
C_{\rm core}&=1 + \frac{\sin{\left(2 u \right)}}{2 u},\\
C_{\rm core}&\simeq 1 + \frac{\sin{\left(2 \cdot \left(0.2617993878\right) \right)}}{2 \cdot \left(0.2617993878\right)}\simeq 1.9549.
\end{aligned}
$$

境界の振幅と減衰率を用いる。

$$
\begin{aligned}
C_{\rm clad}&=\frac{\cos^{2}{\left(u \right)}}{w},\\
C_{\rm clad}&\simeq \frac{\cos^{2}{\left(\left(0.2617993878\right) \right)}}{\left(0.07014893454\right)}\simeq 13.3.
\end{aligned}
$$

全断面での Poynting 流に対するコアの寄与を求める。

$$
\begin{aligned}
\Gamma_{\rm core}&=\frac{C_{\rm core}}{C_{\rm clad} + C_{\rm core}},\\
\Gamma_{\rm core}&\simeq \frac{\left(1.954929659\right)}{\left(13.30045435\right) + \left(1.954929659\right)}\simeq 0.12815.
\end{aligned}
$$

クラッドへ長い尾を持ち、パワーの多くは外側にある。
:::

:::example[中程度の閉じ込め]
コアの中心から境界までに位相が四分の一周期未満だけ変わる条件を取る。

入力は基本偶 TE の横波数の換算値 $u\simeq 0.7854$、コア半厚 $a\simeq 1\times 10^{-6}\,\mathrm m$ である。

境界の対数微分を合わせて外側の減衰を決める。

$$
\begin{aligned}
w&=u \tan{\left(u \right)},\\
w&\simeq \left(0.7853981634\right) \cdot \tan{\left(\left(0.7853981634\right) \right)}\simeq 0.7854.
\end{aligned}
$$

境界値を満たす $V$ を幾何学的に構成する。

$$
\begin{aligned}
V&=\sqrt{u^{2} + w^{2}},\\
V&\simeq \sqrt{\left(0.7853981634\right)^{2} + \left(0.7853981634\right)^{2}}\simeq 1.1107.
\end{aligned}
$$

指数の係数を実際の深さへ換算する。

$$
\begin{aligned}
\ell&=\frac{a}{w},\\
\ell/(\mathrm m)&\simeq \frac{\left(1\times 10^{-6}\right)}{\left(0.7853981634\right)}\simeq 1.2732\times 10^{-6}.
\end{aligned}
$$

中心振幅二乗と半幅を除いた積分を求める。

$$
\begin{aligned}
C_{\rm core}&=1 + \frac{\sin{\left(2 u \right)}}{2 u},\\
C_{\rm core}&\simeq 1 + \frac{\sin{\left(2 \cdot \left(0.7853981634\right) \right)}}{2 \cdot \left(0.7853981634\right)}\simeq 1.6366.
\end{aligned}
$$

境界の振幅と減衰率を用いる。

$$
\begin{aligned}
C_{\rm clad}&=\frac{\cos^{2}{\left(u \right)}}{w},\\
C_{\rm clad}&\simeq \frac{\cos^{2}{\left(\left(0.7853981634\right) \right)}}{\left(0.7853981634\right)}\simeq 0.63662.
\end{aligned}
$$

全断面での Poynting 流に対するコアの寄与を求める。

$$
\begin{aligned}
\Gamma_{\rm core}&=\frac{C_{\rm core}}{C_{\rm clad} + C_{\rm core}},\\
\Gamma_{\rm core}&\simeq \frac{\left(1.636619772\right)}{\left(0.6366197724\right) + \left(1.636619772\right)}\simeq 0.71995.
\end{aligned}
$$

侵入長はなお大きいが、場の境界振幅が下がりコアの割合は増える。
:::

:::example[強い閉じ込めの基本枝]
基本枝で $u$ を $\pi/2$ に近づける。

入力は基本偶 TE の横波数の換算値 $u\simeq 1.309$、コア半厚 $a\simeq 1\times 10^{-6}\,\mathrm m$ である。

境界の対数微分を合わせて外側の減衰を決める。

$$
\begin{aligned}
w&=u \tan{\left(u \right)},\\
w&\simeq \left(1.308996939\right) \cdot \tan{\left(\left(1.308996939\right) \right)}\simeq 4.8852.
\end{aligned}
$$

境界値を満たす $V$ を幾何学的に構成する。

$$
\begin{aligned}
V&=\sqrt{u^{2} + w^{2}},\\
V&\simeq \sqrt{\left(1.308996939\right)^{2} + \left(4.885243083\right)^{2}}\simeq 5.0576.
\end{aligned}
$$

指数の係数を実際の深さへ換算する。

$$
\begin{aligned}
\ell&=\frac{a}{w},\\
\ell/(\mathrm m)&\simeq \frac{\left(1\times 10^{-6}\right)}{\left(4.885243083\right)}\simeq 2.047\times 10^{-7}.
\end{aligned}
$$

中心振幅二乗と半幅を除いた積分を求める。

$$
\begin{aligned}
C_{\rm core}&=1 + \frac{\sin{\left(2 u \right)}}{2 u},\\
C_{\rm core}&\simeq 1 + \frac{\sin{\left(2 \cdot \left(1.308996939\right) \right)}}{2 \cdot \left(1.308996939\right)}\simeq 1.191.
\end{aligned}
$$

境界の振幅と減衰率を用いる。

$$
\begin{aligned}
C_{\rm clad}&=\frac{\cos^{2}{\left(u \right)}}{w},\\
C_{\rm clad}&\simeq \frac{\cos^{2}{\left(\left(1.308996939\right) \right)}}{\left(4.885243083\right)}\simeq 0.013712.
\end{aligned}
$$

全断面での Poynting 流に対するコアの寄与を求める。

$$
\begin{aligned}
\Gamma_{\rm core}&=\frac{C_{\rm core}}{C_{\rm clad} + C_{\rm core}},\\
\Gamma_{\rm core}&\simeq \frac{\left(1.190985932\right)}{\left(0.01371217296\right) + \left(1.190985932\right)}\simeq 0.98862.
\end{aligned}
$$

境界での振幅も侵入長も小さくなる。別の高次モードを許すかは $V$ から判定する。
:::

### 上下のクラッドが違う平板

コアを $0<x<d$、屈折率を $n_1$、左と右のクラッドをそれぞれ $n_a,n_b<n_1$ とする。TE 場に対し $h=\sqrt{k_0^2n_1^2-\beta^2}$、$q_j=\sqrt{\beta^2-k_0^2n_j^2}$ とし、三つとも正実数の範囲を取る。左界面での条件 $E'(0)/E(0)=q_a$ をコアの余弦 $E=C\cos(hx-\phi_a)$ に課すと $\tan\phi_a=q_a/h$。右界面の条件 $E'(d)/E(d)=-q_b$ から $hd-\phi_a=m\pi+\phi_b$、$\tan\phi_b=q_b/h$ となる。従って二つの反射面の条件は
$$
hd=m\pi+\arctan(q_a/h)+\arctan(q_b/h),\qquad m=0,1,\ldots .
\label{eq-optics-asymmetric-slab}
$$
両クラッドが等しいと $hd-2\phi_a=m\pi$。コア中央での余弦の位相が $m\pi/2$ なので、偶数 $m$ は偶モード、奇数 $m$ は奇モードへ戻る。非対称な場合は中心の偶奇で分類できない。

$n_b>n_a$ とする。$\beta\downarrow k_0 n_b$ で $q_b\to0$ だが
$q_a\to q_{a,c}=k_0\sqrt{n_b^2-n_a^2}>0$。また $h$ は $h_c=k_0\sqrt{n_1^2-n_b^2}$ に近づく。$\beta$ を増すと左辺は減少し、二つの逆正接は増加する。$\beta\uparrow k_0n_1$ では左辺はゼロ、右辺の逆正接の和は $\pi$ に近づく。従って $m=0$ でも束縛モードがある必要十分条件は $d>\arctan(q_{a,c}/h_c)/h_c$ である。等号では右クラッドの場が指数減衰しない。対称導波路では基本モードに有限の厚さの下限がないが、クラッドの非対称性でこの下限が現れる。

:::example[クラッド非対称性の下限を測る]
材料と波長を選んだ結果、カットオフで $q_{a,c}=h_c=2\times10^6\,\mathrm{m^{-1}}$ とする。式 :eqref[eq-optics-asymmetric-slab] の基本モードの下限は $(\pi/4)/(2\times10^6\,\mathrm{m^{-1}})\simeq0.393\,\mathrm{\mu m}$。厚さ $0.30\,\mathrm{\mu m}$ では束縛 TE モードがなく、$0.50\,\mathrm{\mu m}$ では基本モードが一つある。次の枝には右辺に $\pi$ が加わるため、後者でも二本目は存在しない。
:::

## 光ファイバーの導波モードと開口数

### 入射角と全反射
半径 $a$ のステップ屈折率ファイバーの端面へ、外部屈折率 $n_0$ から角 $\theta_0$ で入射するメリディオナル光線を考える。コアで軸と角 $\theta_1$ をなすと、側面の入射角は $\pi/2-\theta_1$。全反射には $\cos\theta_1>n_2/n_1$ が必要で、端面の Snell 則 $n_0\sin\theta_0=n_1\sin\theta_1$ から
$$
\mathrm{NA}=\sqrt{n_1^2-n_2^2},\qquad
V=\frac{2\pi a}{\lambda_0}\mathrm{NA}.
\label{eq-optics-fiber-na}
$$
この :term[開口数]{en="numerical aperture" yomi="かいこうすう"} は、幾何光学で $n_0\sin\theta_0<\mathrm{NA}$ を満たす光線が側面で全反射することを表す。$\mathrm{NA}\le n_0$ なら限界角は $\sin\theta_{0,\max}=\mathrm{NA}/n_0$。$\mathrm{NA}>n_0$ なら、この条件だけでは入射半空間の角を制限しない。厳密な入射結合では端面散乱とモードの重なりを評価するので、光線の許容角だけでは結合効率は決まらない。

### ベクトル場の固有方程式
円筒対称の Maxwell 方程式を $e^{i\ell\varphi+i\beta z-i\omega t}$ で分離すると、縦成分 $E_z,H_z$ は各一様層で $F''+r^{-1}F'+[h^2-\ell^2/r^2]F=0$、$h^2=n^2k_0^2-\beta^2$ を満たす。コアは正則な $J_\ell(ur/a)$、クラッドは指数減衰する $K_\ell(wr/a)$ を両縦成分に取る。ここで $u,w$ は平板と同じ式を使う。回転方程式の横成分を連立して解けば
$$
\begin{aligned}
\vb E_\perp&=\frac{i}{h^2}(\beta\grad_\perp E_z-\omega\mu_0\hat{\vb z}\times\grad_\perp H_z),\\
\vb H_\perp&=\frac{i}{h^2}(\beta\grad_\perp H_z+\omega\eps\hat{\vb z}\times\grad_\perp E_z).
\end{aligned}
$$
例えば $\hat{\vb z}\times(i\beta\vb E_\perp-\grad_\perp E_z)=i\omega\mu_0\vb H_\perp$ と、同形の磁場の式を用いて一方を消去する。$r=a$ で $E_z,H_z,E_\varphi,H_\varphi$ の連続を課すとコア二振幅とクラッド二振幅の同次連立方程式を得る。その行列式がゼロとなる $\beta$ が導波モードである。$\ell=0$ なら角度微分がゼロになり、$E_z=0$ の TE と $H_z=0$ の TM が分かれる。$\ell\ne0$ では両縦成分を持つ混成波が一般に現れる。

:::example[角度とモード幅]
$n_1=1.50,n_2=1.48$ では $\mathrm{NA}\simeq0.244$、空気中の許容半角は約 $14.1^\circ$ である。導波路が単一モードの場合、実際の遠方角幅は基本モードの横断面を変換して求める。
:::
:::example[規格化]
同じ材料でコア半径を二倍にし、真空波長も二倍にすると $V$ は変わらない。屈折率の分散を無視すれば、規格化半径 $r/a$ で見た形は変わらない。
:::

### 縦成分の境界値から固有方程式を作る
$J_\ell(u)\ne0$ のとき、コアとクラッドで共通の境界値を $E_b,H_b$ と置く。各縦成分の対数微分を $D_1=(u/a)J_\ell'(u)/J_\ell(u)$、$D_2=(w/a)K_\ell'(w)/K_\ell(w)$ とする。$h_1^2=u^2/a^2,\ h_2^2=-w^2/a^2$ とし、$B=\beta\ell(h_1^{-2}-h_2^{-2})/a$、$C=D_1/h_1^2-D_2/h_2^2$、$D=\eps_1D_1/h_1^2-\eps_2D_2/h_2^2$ と置く。横場の式から $E_\varphi,H_\varphi$ の差は
$$
\begin{pmatrix}-B&-i\omega\mu_0 C\\i\omega D&-B\end{pmatrix}
\binom{E_b}{H_b}=0 .
$$
従って $B^2=\omega^2\mu_0 CD$。$J_\ell(u)=0$ なら両縦成分の境界値がゼロである。$K_\ell(w)>0$ なのでクラッドの両縦成分がゼロとなり、横場の連続性はコア側の両縦成分の法線微分もゼロとする。$J_\ell$ の正根では微分が非零である。値と微分が同時にゼロなら Bessel の微分方程式の初期値の一意性に反するためである。従って $w>0$ の束縛モードの探索では、この対数微分の極を固有値から除く。$\ell=0$ では $B=0$ であり、$C=0$ が TE、$D=0$ が TM になる。混成解では行列の一行から $E_b/H_b$ まで決まる。

### 多数の横モードを数える近似
$V\gg1$ の弱導波円筒で、横の束縛状態を半径 $a$ の面積の中の局所平面波とみなす。横波数の許容円の面積は $\pi(V/a)^2$。横の実空間面積 $\pi a^2$ と掛けて一状態の位相空間面積 $(2\pi)^2$ で割ると、偏光を一つにした数は $V^2/4$ 程度。二偏光まで数えれば $V^2/2$ 程度になる。境界の位相やクラッドへの侵入を数えていないため、この面積則は先頭の大きい $V$ 近似である。カットオフ近くの本数には個別の Bessel の条件を使う。

:::example[多数モードの計数規約]
$V=20$ の弱導波ステップファイバーでは位相空間から偏光一つ当たり約 100、二偏光を含めて約 200 モードになる。同じ光源を入れても、入口モードとの重なりが小さい高次モードは弱くしか励起されない。固有モードの個数と実際に占有したモードの個数は分けて数える。
:::

### 開口数と単一モード条件を同じ図式で使う

開口数を大きくすると光線としては広い角を受け入れられるが、一定の半径と波長では :eqref[eq-optics-fiber-na] の $V$ も増える。弱導波の円形ステップ形では、高次モードが初めて現れる $V_c\simeq2.4048$ を一定とすれば、カットオフ真空波長は $\lambda_c=2\pi a\,\mathrm{NA}/V_c$。$\lambda_0>\lambda_c$ の波長側で単一空間モードになる。これは材料分散を無視する局所設計式であり、波長範囲が広いなら $\mathrm{NA}(\lambda_0)$ を入れて条件を解く。また単一モードでも外部の全受入角を同じ効率で結合できるわけではない。場の重なりは次節で求める。

:::example[空気から長波長を入れる]
弱導波の小コアを通信波長に近い尺度で調べる。

入力はコア屈折率 $n_{1}=\frac{29}{20}$、クラッド屈折率 $n_{2}=\frac{289}{200}$、入射側の屈折率 $n_{0}=1$、コア半径 $a\simeq 4\times 10^{-6}\,\mathrm m$、真空波長 $\lambda_0\simeq 1.55\times 10^{-6}\,\mathrm m$ である。

端面での角をコアの全反射へつなぐ。

$$
\begin{aligned}
\mathrm{NA}&=\sqrt{n_{1}^{2} - n_{2}^{2}},\\
\mathrm{NA}&= \sqrt{\left(\frac{29}{20}\right)^{2} - \left(\frac{289}{200}\right)^{2}}\simeq 0.12031.
\end{aligned}
$$

外部屈折率で割って角を得る。

$$
\begin{aligned}
\theta_{\max}&=\operatorname{asin}{\left(\min\left(1, \frac{\mathrm{NA}}{n_{0}}\right) \right)},\\
\theta_{\max}/(\mathrm{rad})&\simeq \operatorname{asin}{\left(\min\left(1, \frac{\left(0.1203120942\right)}{\left(1\right)}\right) \right)}\simeq 0.1206.
\end{aligned}
$$

半径を真空波数と開口数で規格化する。

$$
\begin{aligned}
V&=\frac{2 \pi \mathrm{NA} a}{\lambda_0},\\
V&\simeq \frac{2 \cdot \pi \cdot \left(0.1203120942\right) \cdot \left(4\times 10^{-6}\right)}{\left(1.55\times 10^{-6}\right)}\simeq 1.9508.
\end{aligned}
$$

弱導波の最初の高次カットオフで波長を解く。

$$
\begin{aligned}
\lambda_c&=\frac{2\pi \mathrm{NA}a}{V_c},\\
\lambda_c/(\mathrm m)&\simeq \frac{1250 \cdot \pi \cdot \left(0.1203120942\right) \cdot \left(4\times 10^{-6}\right)}{1503}\simeq 1.2574\times 10^{-6}.
\end{aligned}
$$

一より長い側が単一空間モードである。

$$
\begin{aligned}
\lambda_0/\lambda_c&=\frac{\lambda_0}{\lambda_c},\\
\lambda_0/\lambda_c&\simeq \frac{\left(1.55\times 10^{-6}\right)}{\left(1.257390523\times 10^{-6}\right)}\simeq 1.2327.
\end{aligned}
$$

光線の角制限と単一モードの判定を別々に得る。
:::

:::example[同じファイバーを短波長で使う]
屈折率の波長変化を無視し、波長をほぼ半分にする。

入力はコア屈折率 $n_{1}=\frac{29}{20}$、クラッド屈折率 $n_{2}=\frac{289}{200}$、入射側の屈折率 $n_{0}=1$、コア半径 $a\simeq 4\times 10^{-6}\,\mathrm m$、真空波長 $\lambda_0\simeq 8\times 10^{-7}\,\mathrm m$ である。

端面での角をコアの全反射へつなぐ。

$$
\begin{aligned}
\mathrm{NA}&=\sqrt{n_{1}^{2} - n_{2}^{2}},\\
\mathrm{NA}&= \sqrt{\left(\frac{29}{20}\right)^{2} - \left(\frac{289}{200}\right)^{2}}\simeq 0.12031.
\end{aligned}
$$

外部屈折率で割って角を得る。

$$
\begin{aligned}
\theta_{\max}&=\operatorname{asin}{\left(\min\left(1, \frac{\mathrm{NA}}{n_{0}}\right) \right)},\\
\theta_{\max}/(\mathrm{rad})&\simeq \operatorname{asin}{\left(\min\left(1, \frac{\left(0.1203120942\right)}{\left(1\right)}\right) \right)}\simeq 0.1206.
\end{aligned}
$$

半径を真空波数と開口数で規格化する。

$$
\begin{aligned}
V&=\frac{2 \pi \mathrm{NA} a}{\lambda_0},\\
V&\simeq \frac{2 \cdot \pi \cdot \left(0.1203120942\right) \cdot \left(4\times 10^{-6}\right)}{\left(8\times 10^{-7}\right)}\simeq 3.7797.
\end{aligned}
$$

弱導波の最初の高次カットオフで波長を解く。

$$
\begin{aligned}
\lambda_c&=\frac{2\pi \mathrm{NA}a}{V_c},\\
\lambda_c/(\mathrm m)&\simeq \frac{1250 \cdot \pi \cdot \left(0.1203120942\right) \cdot \left(4\times 10^{-6}\right)}{1503}\simeq 1.2574\times 10^{-6}.
\end{aligned}
$$

一より長い側が単一空間モードである。

$$
\begin{aligned}
\lambda_0/\lambda_c&=\frac{\lambda_0}{\lambda_c},\\
\lambda_0/\lambda_c&\simeq \frac{\left(8\times 10^{-7}\right)}{\left(1.257390523\times 10^{-6}\right)}\simeq 0.63624.
\end{aligned}
$$

開口数は同じでも高次の導波条件に入る。
:::

:::example[入射側へ高屈折率媒質を置く]
端面の外側を透明媒質へ変える。

入力はコア屈折率 $n_{1}=\frac{29}{20}$、クラッド屈折率 $n_{2}=\frac{289}{200}$、入射側の屈折率 $n_{0}=\frac{133}{100}$、コア半径 $a\simeq 4\times 10^{-6}\,\mathrm m$、真空波長 $\lambda_0\simeq 1.55\times 10^{-6}\,\mathrm m$ である。

端面での角をコアの全反射へつなぐ。

$$
\begin{aligned}
\mathrm{NA}&=\sqrt{n_{1}^{2} - n_{2}^{2}},\\
\mathrm{NA}&= \sqrt{\left(\frac{29}{20}\right)^{2} - \left(\frac{289}{200}\right)^{2}}\simeq 0.12031.
\end{aligned}
$$

外部屈折率で割って角を得る。

$$
\begin{aligned}
\theta_{\max}&=\operatorname{asin}{\left(\min\left(1, \frac{\mathrm{NA}}{n_{0}}\right) \right)},\\
\theta_{\max}/(\mathrm{rad})&\simeq \operatorname{asin}{\left(\min\left(1, \frac{\left(0.1203120942\right)}{\left(\frac{133}{100}\right)}\right) \right)}\simeq 0.090584.
\end{aligned}
$$

半径を真空波数と開口数で規格化する。

$$
\begin{aligned}
V&=\frac{2 \pi \mathrm{NA} a}{\lambda_0},\\
V&\simeq \frac{2 \cdot \pi \cdot \left(0.1203120942\right) \cdot \left(4\times 10^{-6}\right)}{\left(1.55\times 10^{-6}\right)}\simeq 1.9508.
\end{aligned}
$$

弱導波の最初の高次カットオフで波長を解く。

$$
\begin{aligned}
\lambda_c&=\frac{2\pi \mathrm{NA}a}{V_c},\\
\lambda_c/(\mathrm m)&\simeq \frac{1250 \cdot \pi \cdot \left(0.1203120942\right) \cdot \left(4\times 10^{-6}\right)}{1503}\simeq 1.2574\times 10^{-6}.
\end{aligned}
$$

一より長い側が単一空間モードである。

$$
\begin{aligned}
\lambda_0/\lambda_c&=\frac{\lambda_0}{\lambda_c},\\
\lambda_0/\lambda_c&\simeq \frac{\left(1.55\times 10^{-6}\right)}{\left(1.257390523\times 10^{-6}\right)}\simeq 1.2327.
\end{aligned}
$$

角の許容値は変化するが、内部だけで決まる $V$ は変わらない。
:::

## 単一モードファイバーと LP モード

### 弱導波近似
屈折率差が小さいとき横偏光を一定に近似し、一つのスカラー場 $F(r)e^{i\ell\varphi}$ を考える。Maxwell の連続条件を一次の近似で $F,F'$ の連続に置き換える。$J_\ell$ と $K_\ell$ の対数微分を合わせれば $u J_\ell'(u)/J_\ell(u)=w K_\ell'(w)/K_\ell(w)$。円筒 Bessel 関数の微分関係 $xJ_\ell'=xJ_{\ell-1}-\ell J_\ell$、$xK_\ell'=-xK_{\ell-1}-\ell K_\ell$ で
$$
\frac{uJ_{\ell-1}(u)}{J_\ell(u)}=-\frac{wK_{\ell-1}(w)}{K_\ell(w)}
\label{eq-optics-fiber-lp}
$$
となる。近似的な直線偏光モードの組を :term[LP モード]{en="linearly polarized mode" yomi="えるぴーもーど"} と呼ぶ。

カットオフの $w\downarrow0$ で、$\ell=1$ なら $K_0(w)\sim-\log w$、$K_1(w)\sim1/w$ なので右辺はゼロに近づく。$\ell\ge2$ でも $K_\ell\sim2^{\ell-1}(\ell-1)!w^{-\ell}$ から同じである。従って $\ell=1$ の最初のカットオフは $J_0(V)=0$ の第一正根 $2.4048$。$\ell=0$ では $K_{-1}=K_1$、右辺は $-wK_1/K_0\to0^-$、左辺は $-uJ_1/J_0$ となる。基本 $\ell=0$ 解は $u$ が 0 から $J_0$ の第一根までの枝にあり、中心で正の場を持つ。左辺と右辺の一致を $uJ_1/J_0=wK_1/K_0$ と書くと、左辺はこの区間で 0 から増加し、右辺は $w$ とともに増加する。単調性は Bessel の微分方程式から、左辺の $u$ 微分が $u+y^2/u>0$、右辺の $w$ 微分が $y^2/w-w>0$ となることで分かる。後者で $K_1>K_0$ は正の引数での積分表示 $K_\ell(w)=\int_0^\infty e^{-w\cosh t}\cosh(\ell t)\dd t$ による。従って任意の $V>0$ で基本枝は一つある。追加の $\ell=0$ 枝の最初のカットオフは $J_1$ の第一根で、$2.4048$ より大きい。より高い角次数では、スカラー横固有値問題の遠心項 $\ell^2/r^2$ が増え、変分原理で最初の束縛に必要な $V$ は下がらない。ゆえに $V<2.4048$ が単一空間モード条件である。基本モードには二つの偏光自由度が残る。

:::example[単一モード半径]
$\lambda_0=1.55\,\mathrm{\mu m},\mathrm{NA}=0.12$ では $a<2.4048\lambda_0/(2\pi\mathrm{NA})\simeq4.94\,\mathrm{\mu m}$ である。
:::
:::example[波長の切替え]
同じファイバーで波長を約半分にすれば $V$ は約二倍になる。長波長で単一モードだったファイバーに短波長を入れると、別の横モードが導波され得る。
:::

### LP の角次数とベクトルモードの区別

LP の角次数は支配的な横電場の角次数であり、先の柱座標における縦成分の $\ell$ と区別する。例えば $x$ 偏光を柱座標で書くと $\hat{\vb x}=\hat{\vb r}\cos\varphi-\hat{\boldsymbol\varphi}\sin\varphi$ なので、横電場に角依存のない $\mathrm{LP}_{01}$ でも、柱座標の成分と縦成分には角次数 1 が現れる。これが厳密な $\mathrm{HE}_{11}$ の弱導波極限である。$\mathrm{LP}_{11}$ では角方向に $\cos\varphi,\sin\varphi$ の二成分があり、それぞれ二偏光を選べる。柱座標へ換えると角因子の積には角次数 0 と 2 が現れ、厳密なモードでは $\mathrm{TE}_{01}$、$\mathrm{TM}_{01}$、二つの向きを持つ $\mathrm{HE}_{21}$ に分かれる。これらを一組とするのは弱導波極限で伝播定数の差を省いたためである。

特に角次数 0 のベクトル TE 枝と TM 枝では、先の $C,D$ に $J_0'=-J_1$、$K_0'=-K_1$ を入れると、クラッド側には $K_1(w)/(wK_0(w))$ が現れる。$w\downarrow0$ でこれは発散する。コア側がこれと釣り合うには $J_0(u)$ がゼロに近づく必要があり、最初のしきい値はどちらも $V=2.4048\ldots$ である。高次の混成モードでは、有限の屈折率差があるときはベクトルの固有方程式で個別に判定する。

### 規格化された伝播定数
材料ごとに異なる $\beta$ を比較するため、$b=(\beta^2-k_0^2 n_2^2)/[k_0^2(n_1^2-n_2^2)]=w^2/V^2$ を使う。導波モードでは $0<b<1$ であり、指数減衰の枝のカットオフで $b\to0$。厳密に $w=0$ のしきい値ではクラッドの方程式が変わる。弱導波のスカラー円筒では $\ell\ge1$ の減衰側が $r^{-\ell}$ となり、$\ell\ge2$ では $\int_a^\infty |F|^2r\dd r$ が有限になり得る。この境界例は $w>0$ の指数減衰モードの本数には含めない。$b$ が 1 に近いと有効屈折率がコア値に近づく。一方で幅を推定するには $w=V\sqrt b$ を使う必要があり、同じ $b$ だけから物理的な侵入長 $a/w$ は決められない。

:::example[有効屈折率を求める]
$n_1=1.50,n_2=1.48,b=1/2$ なら $n_{\rm eff}^2=(n_1^2+n_2^2)/2=2.2202$。有効屈折率は約 $1.4900$ である。
:::

### 基本モードへの結合を重なりから求める

基本モードの横断面をガウスで近似し、入射も同じ中心、平面波面、同じ偏光のガウスとする。振幅を $U_j=\sqrt{2/\pi}w_j^{-1}\exp(-r^2/w_j^2)$ とすると $\int|U_j|^2\dd^2r=1$。相互内積は $\int_0^\infty U_1U_2\,2\pi r\dd r=2w_1w_2/(w_1^2+w_2^2)$ である。さらに入射中心が横に $d$ だけずれる場合、二つの指数を平方完成すると振幅の重なりに $e^{-d^2/(w_1^2+w_2^2)}$ が掛かる。従ってパワー評価では指数の二倍を取る。実際の LP モードがガウスから外れる場合には、Bessel 型の固有場との重なりを直接積分する。

:::example[幅と中心が一致する]
偏光を合わせた上で二つのガウス半径をそろえる。

入力は入射ガウスの半径 $w_1\simeq 5\times 10^{-6}\,\mathrm m$、基本モードの等価半径 $w_2\simeq 5\times 10^{-6}\,\mathrm m$、横方向の芯ずれ $d=0\,\mathrm m$、端面に届くパワー $P=\frac{1}{1000}\,\mathrm W$ である。

平方完成に現れる幅の尺度をまとめる。

$$
\begin{aligned}
W&=w_1^{2} + w_2^{2},\\
W/(\mathrm{m^2})&\simeq \left(\left(5\times 10^{-6}\right)\right)^{2} + \left(\left(5\times 10^{-6}\right)\right)^{2}\simeq 5\times 10^{-11}.
\end{aligned}
$$

規格化した二つの関数の積を積分する。

$$
\begin{aligned}
A_0&=\frac{2 w_1 w_2}{W},\\
A_0&\simeq \frac{2 \cdot \left(5\times 10^{-6}\right) \cdot \left(5\times 10^{-6}\right)}{\left(5\times 10^{-11}\right)}=1.
\end{aligned}
$$

横に平行移動したときの積分の指数補正を掛ける。

$$
\begin{aligned}
B&=e^{- \frac{d^{2}}{W}},\\
B&\simeq e^{- \frac{\left(0\right)^{2}}{\left(5\times 10^{-11}\right)}}=1.
\end{aligned}
$$

内積で求める振幅の絶対値を二乗する。

$$
\begin{aligned}
\eta&=A_0^{2} B^{2},\\
\eta&= \left(1\right)^{2} \cdot \left(1\right)^{2}=1.
\end{aligned}
$$

端面の反射を無視した重なりによる結合を求める。

$$
\begin{aligned}
P_{\rm mode}&=P \eta,\\
P_{\rm mode}/(\mathrm W)&= \left(\frac{1}{1000}\right) \cdot \left(1\right)=\frac{1}{1000}.
\end{aligned}
$$

重なりは一になり、反射以外の結合損失がない模型になる。
:::

:::example[幅が広すぎる入射]
中心は同じで入射幅を二倍にする。

入力は入射ガウスの半径 $w_1\simeq 1\times 10^{-5}\,\mathrm m$、基本モードの等価半径 $w_2\simeq 5\times 10^{-6}\,\mathrm m$、横方向の芯ずれ $d=0\,\mathrm m$、端面に届くパワー $P=\frac{1}{1000}\,\mathrm W$ である。

平方完成に現れる幅の尺度をまとめる。

$$
\begin{aligned}
W&=w_1^{2} + w_2^{2},\\
W/(\mathrm{m^2})&\simeq \left(\left(1\times 10^{-5}\right)\right)^{2} + \left(\left(5\times 10^{-6}\right)\right)^{2}\simeq 1.25\times 10^{-10}.
\end{aligned}
$$

規格化した二つの関数の積を積分する。

$$
\begin{aligned}
A_0&=\frac{2 w_1 w_2}{W},\\
A_0&\simeq \frac{2 \cdot \left(1\times 10^{-5}\right) \cdot \left(5\times 10^{-6}\right)}{\left(1.25\times 10^{-10}\right)}=\frac{4}{5}.
\end{aligned}
$$

横に平行移動したときの積分の指数補正を掛ける。

$$
\begin{aligned}
B&=e^{- \frac{d^{2}}{W}},\\
B&\simeq e^{- \frac{\left(0\right)^{2}}{\left(1.25\times 10^{-10}\right)}}=1.
\end{aligned}
$$

内積で求める振幅の絶対値を二乗する。

$$
\begin{aligned}
\eta&=A_0^{2} B^{2},\\
\eta&= \left(\frac{4}{5}\right)^{2} \cdot \left(1\right)^{2}=\frac{16}{25}.
\end{aligned}
$$

端面の反射を無視した重なりによる結合を求める。

$$
\begin{aligned}
P_{\rm mode}&=P \eta,\\
P_{\rm mode}/(\mathrm W)&= \left(\frac{1}{1000}\right) \cdot \left(\frac{16}{25}\right)=\frac{2}{3125}.
\end{aligned}
$$

角の制限を満たすように光を絞っても、幅が一致しなければパワーの一部は導波モードへ入らない。
:::

:::example[幅が合い中心がずれる]
前の損失と異なる原因として芯ずれを入れる。

入力は入射ガウスの半径 $w_1\simeq 5\times 10^{-6}\,\mathrm m$、基本モードの等価半径 $w_2\simeq 5\times 10^{-6}\,\mathrm m$、横方向の芯ずれ $d\simeq 3\times 10^{-6}\,\mathrm m$、端面に届くパワー $P=\frac{1}{1000}\,\mathrm W$ である。

平方完成に現れる幅の尺度をまとめる。

$$
\begin{aligned}
W&=w_1^{2} + w_2^{2},\\
W/(\mathrm{m^2})&\simeq \left(\left(5\times 10^{-6}\right)\right)^{2} + \left(\left(5\times 10^{-6}\right)\right)^{2}\simeq 5\times 10^{-11}.
\end{aligned}
$$

規格化した二つの関数の積を積分する。

$$
\begin{aligned}
A_0&=\frac{2 w_1 w_2}{W},\\
A_0&\simeq \frac{2 \cdot \left(5\times 10^{-6}\right) \cdot \left(5\times 10^{-6}\right)}{\left(5\times 10^{-11}\right)}=1.
\end{aligned}
$$

横に平行移動したときの積分の指数補正を掛ける。

$$
\begin{aligned}
B&=e^{- \frac{d^{2}}{W}},\\
B&\simeq e^{- \frac{\left(\left(3\times 10^{-6}\right)\right)^{2}}{\left(5\times 10^{-11}\right)}}\simeq 0.83527.
\end{aligned}
$$

内積で求める振幅の絶対値を二乗する。

$$
\begin{aligned}
\eta&=A_0^{2} B^{2},\\
\eta&\simeq \left(1\right)^{2} \cdot \left(0.8352702114\right)^{2}\simeq 0.69768.
\end{aligned}
$$

端面の反射を無視した重なりによる結合を求める。

$$
\begin{aligned}
P_{\rm mode}&=P \eta,\\
P_{\rm mode}/(\mathrm W)&\simeq \left(\frac{1}{1000}\right) \cdot \left(0.6976763261\right)\simeq 6.9768\times 10^{-4}.
\end{aligned}
$$

幅の因子が一でも、指数の因子による結合損失が残る。
:::

## ファイバーの分散と損失

### 群遅延を二つに分ける
異なるモード $j$ が同じ周波数で異なる $\beta_j(\omega)$ を持つと、距離 $L$ の群遅延は $L\dd\beta_j/\dd\omega$ で異なる。これが :term[モード分散]{en="modal dispersion" yomi="もーどぶんさん"} である。非分散コア中を角 $\theta_1$ で進む光線の軸方向速度は $(c/n_1)\cos\theta_1$。幾何光学の最大遅延差は $n_1L(c^{-1})(n_1/n_2-1)$ になる。弱導波で $n_1/n_2-1\simeq(n_1-n_2)/n_1$ なら差は $(n_1-n_2)L/c$。放物型の屈折率分布では、経路の長い大角度光線が低屈折率の外側を通るので遅延差を減らせる。

単一モードでも $\beta''(\omega)$ によってパルスが広がる（:eqref[eq-optics-gvd]）。有効屈折率 $n_{\rm eff}=\beta/k_0$ は材料屈折率のほか、波長による閉じ込めの変化にも依存する。単位長の群遅延 $1/v_g=\beta'$ の波長微分を $D$ とすると、$\omega=2\pi c/\lambda_0$ の連鎖律から $D=-2\pi c\,\beta''/\lambda_0^2$ である。$D$ はしばしば $\mathrm{ps/(nm\,km)}$ で与えられ、狭帯域の遅延差は $D L\Delta\lambda_0$ になる。

### 損失のパワー表示
$P(L)=P(0)e^{-\alpha L}$ とすればパワー損失を dB で表示した値は $10\log_{10}(P(0)/P(L))=(10/\ln10)\alpha L$。吸収、微小な密度揺らぎによる散乱、曲げによる放射は別の機構である。放射モードへ移った分は検出する導波モードの減衰として数える。

:::example[光ファイバーの遅延幅]
$D=17\,\mathrm{ps/(nm\,km)}$ のファイバーを $10\,\mathrm{km}$、スペクトル幅を $0.1\,\mathrm{nm}$ とすると、端から端までの群遅延の差は約 $17\,\mathrm{ps}$。
:::
:::example[損失を掛ける]
$0.2\,\mathrm{dB/km}$ で $50\,\mathrm{km}$ 進むと $10\,\mathrm{dB}$ のパワー損失となり、出力は入力の $1/10$。振幅はその平方根である。
:::

### 放物分布で遅延差が小さくなる理由
コアを $n(y)=n_1(1-g^2y^2/2)$ と近似する。光線に沿う光路長を $z$ で積分すると $n\sqrt{1+(y')^2}\simeq n_1[1+((y')^2-g^2y^2)/2]$ である。幾何光学の光線 $y=A\cos(gz+\phi)$ に対し、一周期での $(y')^2$ と $g^2y^2$ の積分は等しい。従ってビームの大きい振動振幅による経路長増加と、外側での速度増加がこの次数では相殺する。有限周期数、材料分散、高次の屈折率項を残すと遅延差は再び現れる。

### 導波路分散だけを取り出す
材料の $n_1,n_2$ が一定でも $V= a\omega\sqrt{n_1^2-n_2^2}/c$ が周波数で動く。$n_{\rm eff}(V)=\sqrt{n_2^2+(n_1^2-n_2^2)b(V)}$ と書けば、$\beta=\omega n_{\rm eff}(V)/c$ の微分で $1/v_g=[n_{\rm eff}+V\,\dd n_{\rm eff}/\dd V]/c$。さらに $\beta''=[2V n_{\rm eff}'+V^2 n_{\rm eff}'']/(c\omega)$ となる。材料に分散がない場合にも横の閉じ込めが周波数で変わるため、$\beta''$ は一般にゼロにならない。設計では材料分散と同符号にも逆符号にもなり得るため、両方を含めた $\beta(\omega)$ の二階微分を取る。

:::example[材料が一定でも群遅延は変わる]
ある動作点で $n_{\rm eff}=1.45$、$V\,\dd n_{\rm eff}/\dd V=0.01$ が得られたとする。群屈折率は $1.46$ であり、位相屈折率だけを使った結果より $100\,\mathrm m$ 当たり約 $3.3\,\mathrm{ns}$ 遅れる。材料分散を無視する場合もこの差は残る。
:::

### 距離と帯域を同時に含めた伝送の見積り

パワー減衰が $\alpha_{\rm dB}$ デシベル毎メートルなら、距離 $L$ の透過割合は $10^{-\alpha_{\rm dB}L/10}$。波長表示の分散係数 $D$ は単位長さ当たりの群遅延の波長微分なので、狭い波長分布の標準偏差 $\sigma_\lambda$ に対応する群遅延の標準偏差は $|D|L\sigma_\lambda$。これは遅延関数を一次に展開して分散をとった結果である。入射時刻と波長に相関がない統計モデルなら、出力の時間分散は入射時間分散とこの値の二乗和になる。既にチャープがあるパルスではこの和に共分散項 $2DL\operatorname{Cov}(t,\lambda)$ が入る。波形の強度だけでは共分散は定まらない。

:::example[短距離で入力幅が支配的]
標準偏差の小さい光源を短いファイバーへ入れる。

入力は波長分散係数 $D=17\,\mathrm{ps/(nm\,km)}$、距離 $L=1\,\mathrm{km}$、波長の標準偏差 $\sigma_\lambda=\frac{1}{10}\,\mathrm{nm}$、入射の時間標準偏差 $\sigma_t=20\,\mathrm{ps}$、パワー損失率 $\alpha_{\rm dB}=\frac{1}{5}\,\mathrm{dB/km}$ である。

実用単位では距離の km と波長幅の nm が分母と相殺する。

$$
\begin{aligned}
\sigma_{\rm g}&=L \sigma_\lambda \left|{D}\right|,\\
\sigma_{\rm g}/(\mathrm{ps})&= \left(1\right) \cdot \left(\frac{1}{10}\right) \cdot \left|{\left(17\right)}\right|=\frac{17}{10}.
\end{aligned}
$$

時刻と波長に初期相関がない場合の分散を加える。

$$
\begin{aligned}
\sigma_{\rm out}&=\sqrt{\sigma_{\rm g}^{2} + \sigma_t^{2}},\\
\sigma_{\rm out}/(\mathrm{ps})&= \sqrt{\left(\frac{17}{10}\right)^{2} + \left(20\right)^{2}}\simeq 20.072.
\end{aligned}
$$

実際の変調時間がどの程度伸びるかを求める。

$$
\begin{aligned}
g&=\frac{\sigma_{\rm out}}{\sigma_t},\\
g&\simeq \frac{\left(20.07211997\right)}{\left(20\right)}\simeq 1.0036.
\end{aligned}
$$

一様なファイバーの長さを掛ける。

$$
\begin{aligned}
A&=L \alpha_{\rm dB},\\
A/(\mathrm{dB})&= \left(1\right) \cdot \left(\frac{1}{5}\right)=\frac{1}{5}.
\end{aligned}
$$

デシベルをパワー比へ戻す。

$$
\begin{aligned}
T&=10^{- \frac{A}{10}},\\
T&= 10^{- \frac{\left(\frac{1}{5}\right)}{10}}\simeq 0.95499.
\end{aligned}
$$

群遅延の広がりは入力幅より小さく、損失も小さい。
:::

:::example[長距離に伸ばす]
光源の特性を変えず減衰と分散の両方を積算する。

入力は波長分散係数 $D=17\,\mathrm{ps/(nm\,km)}$、距離 $L=50\,\mathrm{km}$、波長の標準偏差 $\sigma_\lambda=\frac{1}{10}\,\mathrm{nm}$、入射の時間標準偏差 $\sigma_t=20\,\mathrm{ps}$、パワー損失率 $\alpha_{\rm dB}=\frac{1}{5}\,\mathrm{dB/km}$ である。

実用単位では距離の km と波長幅の nm が分母と相殺する。

$$
\begin{aligned}
\sigma_{\rm g}&=L \sigma_\lambda \left|{D}\right|,\\
\sigma_{\rm g}/(\mathrm{ps})&= \left(50\right) \cdot \left(\frac{1}{10}\right) \cdot \left|{\left(17\right)}\right|=85.
\end{aligned}
$$

時刻と波長に初期相関がない場合の分散を加える。

$$
\begin{aligned}
\sigma_{\rm out}&=\sqrt{\sigma_{\rm g}^{2} + \sigma_t^{2}},\\
\sigma_{\rm out}/(\mathrm{ps})&= \sqrt{\left(85\right)^{2} + \left(20\right)^{2}}\simeq 87.321.
\end{aligned}
$$

実際の変調時間がどの程度伸びるかを求める。

$$
\begin{aligned}
g&=\frac{\sigma_{\rm out}}{\sigma_t},\\
g&\simeq \frac{\left(87.32124598\right)}{\left(20\right)}\simeq 4.3661.
\end{aligned}
$$

一様なファイバーの長さを掛ける。

$$
\begin{aligned}
A&=L \alpha_{\rm dB},\\
A/(\mathrm{dB})&= \left(50\right) \cdot \left(\frac{1}{5}\right)=10.
\end{aligned}
$$

デシベルをパワー比へ戻す。

$$
\begin{aligned}
T&=10^{- \frac{A}{10}},\\
T&= 10^{- \frac{\left(10\right)}{10}}=\frac{1}{10}.
\end{aligned}
$$

パワーが減るだけでなく、時間波形も大きく広がる。
:::

:::example[波長幅を狭める]
入射の時間幅との両立が可能な統計モデルで、波長幅だけを前例より狭める。

入力は波長分散係数 $D=17\,\mathrm{ps/(nm\,km)}$、距離 $L=50\,\mathrm{km}$、波長の標準偏差 $\sigma_\lambda=\frac{7}{100}\,\mathrm{nm}$、入射の時間標準偏差 $\sigma_t=20\,\mathrm{ps}$、パワー損失率 $\alpha_{\rm dB}=\frac{1}{5}\,\mathrm{dB/km}$ である。

実用単位では距離の km と波長幅の nm が分母と相殺する。

$$
\begin{aligned}
\sigma_{\rm g}&=L \sigma_\lambda \left|{D}\right|,\\
\sigma_{\rm g}/(\mathrm{ps})&= \left(50\right) \cdot \left(\frac{7}{100}\right) \cdot \left|{\left(17\right)}\right|=\frac{119}{2}.
\end{aligned}
$$

時刻と波長に初期相関がない場合の分散を加える。

$$
\begin{aligned}
\sigma_{\rm out}&=\sqrt{\sigma_{\rm g}^{2} + \sigma_t^{2}},\\
\sigma_{\rm out}/(\mathrm{ps})&= \sqrt{\left(\frac{119}{2}\right)^{2} + \left(20\right)^{2}}\simeq 62.771.
\end{aligned}
$$

実際の変調時間がどの程度伸びるかを求める。

$$
\begin{aligned}
g&=\frac{\sigma_{\rm out}}{\sigma_t},\\
g&\simeq \frac{\left(62.77141069\right)}{\left(20\right)}\simeq 3.1386.
\end{aligned}
$$

一様なファイバーの長さを掛ける。

$$
\begin{aligned}
A&=L \alpha_{\rm dB},\\
A/(\mathrm{dB})&= \left(50\right) \cdot \left(\frac{1}{5}\right)=10.
\end{aligned}
$$

デシベルをパワー比へ戻す。

$$
\begin{aligned}
T&=10^{- \frac{A}{10}},\\
T&= 10^{- \frac{\left(10\right)}{10}}=\frac{1}{10}.
\end{aligned}
$$

波長分散による幅が下がる一方、同じ長さでのパワー比は変わらない。
:::

## モード結合理論と方向性結合器

### ゆっくり変わるモード振幅
同じ偏光の二局在モードを持つ弱結合系を、スカラー近軸方程式で表す。基準の伝播定数 $\bar\beta$ を除くと $\partial_z U=i \widehat B U$ と書ける。$\widehat B=(\nabla_\perp^2+k_0^2n^2-\bar\beta^2)/(2\bar\beta)$ は損失を無視すればエルミートである。二つの局在モードを直交規格化した基底へこの式を射影すると、各行列成分は $\int u_i^*\widehat B u_j\dd^2r$ となる。元の局在場の重なりを残して直接射影するなら、左辺には重なり行列 $S_{ij}=\int u_i^*u_j\dd^2r$ が掛かり、保存量は $\vb a^\dagger S\vb a$ になる。直交基底では平均の位相を再び除き、基底の位相を調整して一定な非対角成分を実の $\kappa$ に選ぶと
$$
\frac{\dd}{\dd z}\binom{a_1}{a_2}=i
\begin{pmatrix}\Delta/2&\kappa\\\kappa&-\Delta/2\end{pmatrix}\binom{a_1}{a_2}.
\label{eq-optics-coupled-modes}
$$
$\Delta$ は位相不整合に対応する対角差である。右の実行列を $K$ とすれば $K^2=\Omega^2 I$、$\Omega=\sqrt{\kappa^2+(\Delta/2)^2}$ である。指数関数を偶数項と奇数項に分けると $\cos(\Omega z)I+i\sin(\Omega z)K/\Omega$。入力 $(1,0)$ で第 2 のパワーは $(\kappa/\Omega)^2\sin^2(\Omega z)$ である。これが :term[方向性結合器]{en="directional coupler" yomi="ほうこうせいけつごうき"} の基本式である。エルミート性から総パワーが保存される。

:::example[完全移行]
$\Delta=0$、$\kappa=100\,\mathrm{m^{-1}}>0$ なら最初の完全移行長は $\pi/(2\kappa)\simeq15.7\,\mathrm{mm}$、等分配はその半分で得られる。
:::
:::example[不整合の上限]
$\Delta=2\kappa$ では $\Omega^2=2\kappa^2$ なので、長さをどれだけ調整しても移行は最大 $1/2$ である。
:::

### 対称・反対称モードによる見方
$\Delta=0$ では :eqref[eq-optics-coupled-modes] の固有ベクトルは $(1,\pm1)/\sqrt2$、定数の差は $2\kappa$。片側入力はその二モードを等しく重ねたものである。伝播すると位相差が $2\kappa z$ となり、$\pi$ の地点で入力側が打ち消され、他方が強め合う。個々の固有モードの強度分布は伝播中に変わらず、電力移行は両方を重ねたときに見える。

### 非一様結合器で位相を積分する
対角差 $\Delta=0$ のまま $\kappa$ が $z$ で滑らかに変わるとする。全ての位置で行列が同じ $\begin{pmatrix}0&1\\1&0\end{pmatrix}$ の定数倍だから互いに可換である。従って :eqref[eq-optics-coupled-modes] の解の $\kappa z$ を $\int_0^z\kappa(z')\dd z'$ に置けばよい。一般に $\Delta(z)/\kappa(z)$ が変わる場合は微分行列同士が可換でなく、この積分だけでは解を作れない。長手方向を薄い区間に分け、各区間での一定係数の指数を通過順に掛けることで近似できる。

:::example[二つの結合区間の間に空間を置く]
位相整合した二導波路で $\kappa=100\,\mathrm{m^{-1}}$ の区間を二つ、それぞれ長さ $\pi/(400)\,\mathrm m$ とする。間の区間で $\Delta=0,\kappa=0$ なら相対位相は変わらない。合計結合角は $\pi/4+\pi/4=\pi/2$ で、第二導波路に全パワーが移る。
:::

### 離調がある結合器の最短移行長

共通位相を除いた振幅の方程式を $\vb a'=iH\vb a$、$H=\begin{pmatrix}\delta&\kappa\\\kappa&-\delta\end{pmatrix}$ とする。$H^2=(\kappa^2+\delta^2)I$ なので指数の級数の偶数項と奇数項がそれぞれ余弦と正弦を与える。最初に第一導波路だけを励起すれば、第二の強度は $\kappa^2\sin^2(qz)/q^2$、$q=\sqrt{\kappa^2+\delta^2}$。最初の最大は $z=\pi/(2q)$、最大値は $\kappa^2/q^2$ である。長さの調節で補えるのは正弦の因子であり、離調で低下した最大値を取り戻すには伝播定数の差を調整する。$\kappa,\delta$ は実で吸収を無視し、各振幅をパワーで規格化した。

:::example[整合した移行]
離調をゼロにし、完全移行となる長さを選ぶ。

入力は結合係数 $\kappa=100\,\mathrm{m^{-1}}$、伝播定数差の半分 $\delta=0\,\mathrm{m^{-1}}$、結合区間の長さ $L\simeq 0.015708\,\mathrm m$ である。

$H^2$ のスカラー成分の平方根を求める。

$$
\begin{aligned}
q&=\sqrt{\delta^{2} + \kappa^{2}},\\
q/(\mathrm{m^{-1}})&= \sqrt{\left(0\right)^{2} + \left(100\right)^{2}}=100.
\end{aligned}
$$

振幅の非対角成分だけが第二路へパワーを渡す。

$$
\begin{aligned}
\eta_{\max}&=\frac{\kappa^{2}}{q^{2}},\\
\eta_{\max}&= \frac{\left(100\right)^{2}}{\left(100\right)^{2}}=1.
\end{aligned}
$$

正弦二乗が初めて一になる位置を求める。

$$
\begin{aligned}
L_{\rm pk}&=\frac{\pi}{2 q},\\
L_{\rm pk}/(\mathrm m)&= \frac{\pi}{2 \cdot \left(100\right)}\simeq 0.015708.
\end{aligned}
$$

設計した長さで結合振動を評価する。

$$
\begin{aligned}
\eta&=\eta_{\max} \sin^{2}{\left(L q \right)},\\
\eta&\simeq \left(1\right) \cdot \sin^{2}{\left(\left(0.01570796327\right) \cdot \left(100\right) \right)}=1.
\end{aligned}
$$

Hermite 行列からの発展はユニタリなので二路の和は一になる。

$$
\begin{aligned}
R&=1 - \eta,\\
R&= 1 - \left(1\right)=0.
\end{aligned}
$$

上限と実際の移行がともに一になる。
:::

:::example[同じ素子の離調]
長さはそのままで導波路の差を増やす。

入力は結合係数 $\kappa=100\,\mathrm{m^{-1}}$、伝播定数差の半分 $\delta=100\,\mathrm{m^{-1}}$、結合区間の長さ $L\simeq 0.015708\,\mathrm m$ である。

$H^2$ のスカラー成分の平方根を求める。

$$
\begin{aligned}
q&=\sqrt{\delta^{2} + \kappa^{2}},\\
q/(\mathrm{m^{-1}})&= \sqrt{\left(100\right)^{2} + \left(100\right)^{2}}\simeq 141.42.
\end{aligned}
$$

振幅の非対角成分だけが第二路へパワーを渡す。

$$
\begin{aligned}
\eta_{\max}&=\frac{\kappa^{2}}{q^{2}},\\
\eta_{\max}&\simeq \frac{\left(100\right)^{2}}{\left(141.4213562\right)^{2}}=\frac{1}{2}.
\end{aligned}
$$

正弦二乗が初めて一になる位置を求める。

$$
\begin{aligned}
L_{\rm pk}&=\frac{\pi}{2 q},\\
L_{\rm pk}/(\mathrm m)&\simeq \frac{\pi}{2 \cdot \left(141.4213562\right)}\simeq 0.011107.
\end{aligned}
$$

設計した長さで結合振動を評価する。

$$
\begin{aligned}
\eta&=\eta_{\max} \sin^{2}{\left(L q \right)},\\
\eta&\simeq \left(\frac{1}{2}\right) \cdot \sin^{2}{\left(\left(0.01570796327\right) \cdot \left(141.4213562\right) \right)}\simeq 0.31656.
\end{aligned}
$$

Hermite 行列からの発展はユニタリなので二路の和は一になる。

$$
\begin{aligned}
R&=1 - \eta,\\
R&\simeq 1 - \left(0.3165638355\right)\simeq 0.68344.
\end{aligned}
$$

移行の上限が半分となり、その上限を取る長さも前例からずれる。
:::

:::example[離調時の長さを調整]
離調は残したまま最初の最大位置で止める。

入力は結合係数 $\kappa=100\,\mathrm{m^{-1}}$、伝播定数差の半分 $\delta=100\,\mathrm{m^{-1}}$、結合区間の長さ $L\simeq 0.011107\,\mathrm m$ である。

$H^2$ のスカラー成分の平方根を求める。

$$
\begin{aligned}
q&=\sqrt{\delta^{2} + \kappa^{2}},\\
q/(\mathrm{m^{-1}})&= \sqrt{\left(100\right)^{2} + \left(100\right)^{2}}\simeq 141.42.
\end{aligned}
$$

振幅の非対角成分だけが第二路へパワーを渡す。

$$
\begin{aligned}
\eta_{\max}&=\frac{\kappa^{2}}{q^{2}},\\
\eta_{\max}&\simeq \frac{\left(100\right)^{2}}{\left(141.4213562\right)^{2}}=\frac{1}{2}.
\end{aligned}
$$

正弦二乗が初めて一になる位置を求める。

$$
\begin{aligned}
L_{\rm pk}&=\frac{\pi}{2 q},\\
L_{\rm pk}/(\mathrm m)&\simeq \frac{\pi}{2 \cdot \left(141.4213562\right)}\simeq 0.011107.
\end{aligned}
$$

設計した長さで結合振動を評価する。

$$
\begin{aligned}
\eta&=\eta_{\max} \sin^{2}{\left(L q \right)},\\
\eta&\simeq \left(\frac{1}{2}\right) \cdot \sin^{2}{\left(\left(0.01110720735\right) \cdot \left(141.4213562\right) \right)}=\frac{1}{2}.
\end{aligned}
$$

Hermite 行列からの発展はユニタリなので二路の和は一になる。

$$
\begin{aligned}
R&=1 - \eta,\\
R&= 1 - \left(\frac{1}{2}\right)=\frac{1}{2}.
\end{aligned}
$$

上限までは達するが完全移行できない。
:::

## フォトニック結晶入門

### 周期的な層と Bloch 因子
:term[フォトニック結晶]{en="photonic crystal" yomi="ふぉとにっくけっしょう"} は、屈折率などを波長と同程度の周期で変える構造である。一次元で周期 $\Lambda=d_H+d_L$ の二層を繰り返す。左端から右端の状態ベクトルへの伝播行列 $M$ は、式 :eqref[eq-optics-film-matrix] の二層積の逆行列で $\det M=1$。周期点を一つずつ進む固有状態を Bloch 波と呼び、対応する乗数を $e^{iK_B\Lambda}$ と書く。二固有値の積は 1 なので
$$
\cos(K_B\Lambda)=\tfrac12\tr M
=\cos\phi_H\cos\phi_L
-\tfrac12(Y_H/Y_L+Y_L/Y_H)\sin\phi_H\sin\phi_L .
$$
ここで $\phi_j=k_0 n_j d_j$。$|\tr M/2|<1$ では実の Bloch 波数を持つ二根があり、$|\tr M/2|>1$ では二乗数の一方が減衰する。等号では固有値が重なり、固有波が存在しても一般の状態には周期数に比例する項が出ることがある。後者が :term[フォトニックバンドギャップ]{en="photonic band gap" yomi="ふぉとにっくばんどぎゃっぷ"} の一次元垂直入射の例である。周期性を局所的に破ると、両側で指数減衰する解を欠陥の境界で接続できることがある。その周波数が有限多層膜の狭い透過線を作る。

:::example[四分の一波長]
$\phi_H=\phi_L=\pi/2$、$Y_H/Y_L=2$ なら $\tr M/2=-5/4$。固有値は $-2,-1/2$ で、半無限の結晶では振幅が一周期ごとに $-1/2$ 倍の枝を取る。
:::
:::example[一様媒質]
$Y_H=Y_L$ なら右辺は $\cos(\phi_H+\phi_L)$。どこも同じ屈折率であれば、単に周期点を置いても禁制帯は出ない。
:::

### 等光学厚の一次元禁制帯
両層の光学厚が同じなら、分散を無視したとき $\phi_H=\phi_L=\phi$。$Q=(Y_H/Y_L+Y_L/Y_H)/2$ と置くと $\tr M/2=1-(1+Q)\sin^2\phi$。禁制帯の境目は $\sin^2\phi=2/(1+Q)$ である。$Q>1$ のとき $\phi=\pi/2$ の近くで禁制帯となる。$Q=1$ ならこの条件は一点で等号になるだけで幅を持たない。したがって屈折率比を大きくすると中心の減衰だけでなく禁制帯の幅も増す。

### 禁制帯での有限列の透過
停止帯で Bloch 波数が複素数となり、減衰側の一周期乗数の絶対値を $e^{-\kappa_B\Lambda}$ とする。有限な $N$ 周期の全行列には $e^{\pm\kappa_B N\Lambda}$ の二成分がある。両端で外向き条件を合わせると、端で特別な界面共鳴がない限り、透過振幅の分母に増大成分が入り、パワー透過は $e^{-2\kappa_B N\Lambda}$ を尺度として小さくなる。無限周期の禁制帯と有限枚数の反射率は同じ概念ではなく、後者では外側のアドミタンスに応じて前置因子が変わる。設計波長の二層積を使うとこの指数を直接確認できる。

:::example[停止帯での減衰を見積もる]
四分の一波長の高低層でアドミタンス比が 2 なら、中心周波数の減衰固有値は $-1/2$。5 周期の減衰枝の振幅比は $-1/32$、二乗絶対値は $1/1024$ である。有限結晶の透過率ではこれに外部境界の係数が掛かる。
:::

### 有限の四分の一波長膜で境界を接続する

高、低の順で一周期を取り、両側の媒質を同じアドミタンス $Y_0$ とする。設計周波数では各層の位相厚が $\pi/2$ なので、左から右への行列を掛けると一周期は $-\operatorname{diag}(r,r^{-1})$、$r=Y_H/Y_L$ である。入射振幅を 1、反射を $r_N$、透過を $t_N$ と置くと、左端の場は $(1+r_N,Y_0(1-r_N))$、右端は $(t_N,Y_0t_N)$。全行列を左端のベクトルに掛ければ
$$
t_N=(-1)^N r^N(1+r_N)=(-1)^N r^{-N}(1-r_N).
$$
これを解くと $T_N=|t_N|^2=4/(r^N+r^{-N})^2$。先の $r=2$、5 周期の例なら $T_5\simeq0.00390$ となる。一周期の減衰固有値だけから得た $1/1024$ との差は、両端の進行波と Bloch 波を接続した結果である。

### 周期行列のトレースから複素 Bloch 波数を出す

法線入射で無損失な二層の位相厚をともに $\delta$ とし、アドミタンス比を $r>0$ とする。各伝播行列の対角は $\cos\delta$、非対角は $iY^{-1}\sin\delta,iY\sin\delta$ だから、積の跡の半分は $F=\cos^2\delta-(r+r^{-1})\sin^2\delta/2$。行列式は一なので Bloch 固有値は $\xi_\pm=F\pm\sqrt{F^2-1}$。$|F|\le1$ なら絶対値一で伝播できる。$F<-1$ では減衰する根が負実数であり、一周期ごとに符号が反転する。以下では $F<-1$ なら減衰固有波、$|F|<1$ なら二つの Bloch 固有波の片方を選ぶ。選んだ乗数 $\xi$ から、$K_B\Lambda=\arg\xi-i\ln|\xi|$ で複素波数を得る。例えば $F<-1$ の減衰枝では実部を $\pi/\Lambda$、虚部を正に取る。$Q_N$ を選んだ固有波の振幅二乗の比とする。有限膜の入射透過には入射端と出射端のモード接続も加わる。

:::example[四分の一波長の繰返し]
二層の位相厚を $\pi/2$ とする。

入力は二層のアドミタンス比 $r=2$、一層の位相厚 $\delta\simeq 1.5708\,\mathrm{rad}$、周期の数 $N=5$ である。

各層での進行と反射を二層の行列の積にまとめる。

$$
\begin{aligned}
F&=- \frac{\left(r + \frac{1}{r}\right) \sin^{2}{\left(\delta \right)}}{2} + \cos^{2}{\left(\delta \right)},\\
F&\simeq - \frac{\left(\left(2\right) + \frac{1}{\left(2\right)}\right) \cdot \sin^{2}{\left(\left(1.570796327\right) \right)}}{2} + \cos^{2}{\left(\left(1.570796327\right) \right)}=- \frac{5}{4}.
\end{aligned}
$$

固有値の絶対値が一から外れる場合を判定する。

$$
\begin{aligned}
S&=F^{2} - 1,\\
S&= \left(- \frac{5}{4}\right)^{2} - 1=\frac{9}{16}.
\end{aligned}
$$

以下の例は通過帯または $F<-1$ の阻止帯を取る。

$$
\begin{aligned}
\xi&=F + \sqrt{S},\\
\xi&= \left(- \frac{5}{4}\right) + \sqrt{\left(\frac{9}{16}\right)}=- \frac{1}{2}.
\end{aligned}
$$

固有値の絶対値を指数の減衰率へ換算する。

$$
\begin{aligned}
\kappa_B\Lambda&=- \log{\left(\left|{\xi}\right| \right)},\\
\kappa_B\Lambda&= - \log{\left(\left|{\left(- \frac{1}{2}\right)}\right| \right)}\simeq 0.69315.
\end{aligned}
$$

周期ごとの振幅比の積を二乗する。

$$
\begin{aligned}
Q_N&=\left|{\xi}\right|^{2 N},\\
Q_N&= \left|{\left(- \frac{1}{2}\right)}\right|^{2 \cdot \left(5\right)}=\frac{1}{1024}.
\end{aligned}
$$

一周期で符号が反転し、振幅はアドミタンス比の逆数になる。
:::

:::example[中心から離れた阻止帯]
両層の位相厚を同じ比で減らす。

入力は二層のアドミタンス比 $r=2$、一層の位相厚 $\delta\simeq 1.309\,\mathrm{rad}$、周期の数 $N=5$ である。

各層での進行と反射を二層の行列の積にまとめる。

$$
\begin{aligned}
F&=- \frac{\left(r + \frac{1}{r}\right) \sin^{2}{\left(\delta \right)}}{2} + \cos^{2}{\left(\delta \right)},\\
F&\simeq - \frac{\left(\left(2\right) + \frac{1}{\left(2\right)}\right) \cdot \sin^{2}{\left(\left(1.308996939\right) \right)}}{2} + \cos^{2}{\left(\left(1.308996939\right) \right)}\simeq -1.0993.
\end{aligned}
$$

固有値の絶対値が一から外れる場合を判定する。

$$
\begin{aligned}
S&=F^{2} - 1,\\
S&\simeq \left(-1.099278579\right)^{2} - 1\simeq 0.20841.
\end{aligned}
$$

以下の例は通過帯または $F<-1$ の阻止帯を取る。

$$
\begin{aligned}
\xi&=F + \sqrt{S},\\
\xi&\simeq \left(-1.099278579\right) + \sqrt{\left(0.2084133948\right)}\simeq -0.64276.
\end{aligned}
$$

固有値の絶対値を指数の減衰率へ換算する。

$$
\begin{aligned}
\kappa_B\Lambda&=- \log{\left(\left|{\xi}\right| \right)},\\
\kappa_B\Lambda&\simeq - \log{\left(\left|{\left(-0.6427554201\right)}\right| \right)}\simeq 0.44199.
\end{aligned}
$$

周期ごとの振幅比の積を二乗する。

$$
\begin{aligned}
Q_N&=\left|{\xi}\right|^{2 N},\\
Q_N&\simeq \left|{\left(-0.6427554201\right)}\right|^{2 \cdot \left(5\right)}\simeq 0.012035.
\end{aligned}
$$

同じ材料でも減衰率は中心より小さい。
:::

:::example[通過帯へ移る]
位相厚をさらに小さく取る。

入力は二層のアドミタンス比 $r=2$、一層の位相厚 $\delta\simeq 0.7854\,\mathrm{rad}$、周期の数 $N=5$ である。

各層での進行と反射を二層の行列の積にまとめる。

$$
\begin{aligned}
F&=- \frac{\left(r + \frac{1}{r}\right) \sin^{2}{\left(\delta \right)}}{2} + \cos^{2}{\left(\delta \right)},\\
F&\simeq - \frac{\left(\left(2\right) + \frac{1}{\left(2\right)}\right) \cdot \sin^{2}{\left(\left(0.7853981634\right) \right)}}{2} + \cos^{2}{\left(\left(0.7853981634\right) \right)}=- \frac{1}{8}.
\end{aligned}
$$

固有値の絶対値が一から外れる場合を判定する。

$$
\begin{aligned}
S&=F^{2} - 1,\\
S&= \left(- \frac{1}{8}\right)^{2} - 1=- \frac{63}{64}.
\end{aligned}
$$

以下の例は通過帯または $F<-1$ の阻止帯を取る。

$$
\begin{aligned}
\xi&=F + \sqrt{S},\\
\xi&= \left(- \frac{1}{8}\right) + \sqrt{\left(- \frac{63}{64}\right)}\simeq -0.125+0.99216i.
\end{aligned}
$$

固有値の絶対値を指数の減衰率へ換算する。

$$
\begin{aligned}
\kappa_B\Lambda&=- \log{\left(\left|{\xi}\right| \right)},\\
\kappa_B\Lambda&\simeq - \log{\left(\left|{\left(-0.125+0.9921567416i\right)}\right| \right)}=0.
\end{aligned}
$$

周期ごとの振幅比の積を二乗する。

$$
\begin{aligned}
Q_N&=\left|{\xi}\right|^{2 N},\\
Q_N&\simeq \left|{\left(-0.125+0.9921567416i\right)}\right|^{2 \cdot \left(5\right)}=1.
\end{aligned}
$$

固有値は単位円上にある。有限構造の界面反射とは別の条件である。
:::

## 演習問題

::::exercise{level=basic}
対称平板導波路で $V=1$ のとき奇 TE の導波が可能か。

:::solution
最初の奇の正の左辺が出るのは $u>\pi/2$。$u<V=1$ なので条件を満たせない。
:::
::::

::::exercise{level=basic}
$n_1=1.5,n_2=1.2$ のステップ屈折率ファイバーの NA を求めよ。

:::solution
式 :eqref[eq-optics-fiber-na] で $\sqrt{2.25-1.44}=\sqrt{0.81}=0.90$ となる。
:::
::::

::::exercise{level=basic}
$0.5\,\mathrm{dB/km}$ の導波路を $20\,\mathrm{km}$ 進んだパワー比を求めよ。

:::solution
損失が $10\,\mathrm{dB}$ なので、パワー比は $10^{-10/10}=0.1$ である。
:::
::::

::::exercise{level=basic}
$\Delta=0,\kappa>0$ の方向性結合器で半分が初めて移る長さを求めよ。

:::solution
移行率が $\sin^2(\kappa z)$。最初の半分は $\kappa z=\pi/4$ で、$z=\pi/(4\kappa)$ である。
:::
::::

::::exercise{level=advanced}
$V=2$ の円筒ファイバーで LP の空間モード数を弱導波近似で答えよ。

:::solution
最初の高次カットオフが $2.4048$ で、$2<2.4048$。従って空間モードは基本の一つであり、別に偏光自由度がある。
:::
::::

::::exercise{level=advanced}
式 :eqref[eq-optics-coupled-modes] から $\abs{a_1}^2+\abs{a_2}^2$ の保存を示せ。

:::solution
$a'=iKa$ とその共役転置 $a'^\dagger=-ia^\dagger K$ を積の微分に用いると、$(a^\dagger a)'=-ia^\dagger Ka+ia^\dagger Ka=0$ になる。
:::
::::

::::exercise{level=advanced}
$D=20\,\mathrm{ps/(nm\,km)},L=5\,\mathrm{km}$ で波長差 $0.2\,\mathrm{nm}$ の二パルスの遅延差を求めよ。

:::solution
線形近似では $D L\Delta\lambda_0=20(5)(0.2)\,\mathrm{ps}=20\,\mathrm{ps}$。
:::
::::

::::exercise{level=advanced}
$\tr M/2=3/2$ の無損失周期系で一周期当たりの減衰乗数を求めよ。

:::solution
行列式 1 から固有値は $(3\pm\sqrt5)/2$。正で絶対値が 1 未満の枝は $(3-\sqrt5)/2$ である。
:::
::::

::::exercise{level=advanced}
非対称平板の TE モードで、カットオフの横波数が $h_c=4\times10^6\,\mathrm{m^{-1}}$、強く閉じ込める側の減衰率が $q_{a,c}=\sqrt3\,h_c$ とする。厚さ $0.30\,\mathrm{\mu m}$ で基本枝と次の枝のどちらが導波するか求めよ。

:::solution
式 :eqref[eq-optics-asymmetric-slab] のカットオフで、弱く閉じ込める側の減衰率はゼロ、他方の反射位相に対応する逆正接は $\arctan\sqrt3=\pi/3$。基本枝の条件は $h_cd>\pi/3$、次の枝は $h_cd>4\pi/3$ である。与えられた厚さでは $h_cd=1.20$。$\pi/3\simeq1.047<1.20<4\pi/3\simeq4.189$ より、基本枝だけが導波する。
:::
::::

## 確認テスト

::::quiz
誘電体導波コアの束縛モードの $\beta$ の範囲を選べ。

- [x] $k_0 n_2<\beta<k_0 n_1$
- [ ] $\beta<k_0 n_2$
- [ ] $\beta>k_0 n_1$

:::explanation
クラッドの横波数二乗は負でコアの横波数二乗は正である。その両条件を合わせるとこの開区間を得る。
:::
::::

::::quiz{answer="0.01" tol="1e-6"}
長さを変えて損失が $20\,\mathrm{dB}$ となった。出力のパワー比を答えよ。

:::explanation
$20=-10\log_{10}(P/P_0)$ なので $P/P_0=10^{-2}$。
:::
::::

::::quiz
単一空間モードのファイバーでも起こるものを選べ。

- [x] 群速度分散
- [ ] 空間モード間の群遅延差
- [ ] 別モードへの伝播が必須になること

:::explanation
同じ一つのモードでも $\beta'(\omega)$ が周波数で変わり、パルス中の異なる周波数が異なる遅延を受ける。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
$\Delta=2\kappa$ の方向性結合器で他方に移せるパワー比の最大値を答えよ。

:::explanation
式の最大は $\kappa^2/[\kappa^2+(\Delta/2)^2]$。与条件で $1/2$ である。
:::
::::

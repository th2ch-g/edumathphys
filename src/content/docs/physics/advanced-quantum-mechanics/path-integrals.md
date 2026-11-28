---
title: 経路積分
description: 時間発展を時間分割し、実時間の干渉と虚時間の重みを通じて伝播とトンネル効果を計算する。
notation:
  - symbol: 'K(x_b,t_b;x_a,t_a)'
    meaning: 実時間の伝播核
  - symbol: 'S[x]'
    meaning: 実時間の作用
  - symbol: 'S_E[x]'
    meaning: Euclid 作用
  - symbol: '\beta'
    meaning: 逆温度 1/(k_B T)
---

## この章の目標

- 伝播核の初期条件と合成則を説明する。
- 時間分割で自由粒子と調和振動子の核を求める。
- 半古典近似と虚時間を使い、古典作用と量子の揺らぎを分ける。
- 電磁場による相対位相とトンネルの指数因子を計算する。

## 伝播関数と時間発展核

点 $x_a$ を中心とする波束を準備したとき、後で $x_b$ を測る振幅は時間発展が作る核で与えられる。点の状態はデルタ規格化された一般化固有状態なので、絶対値二乗をそのまま確率とは見なさず、正規化した波束に積分して使う。

:::definition[時間発展核]{#def-aqm-propagator}
時間に依存しない自己共役 Hamiltonian $H$ に対し $U(T)=e^{-iHT/\hbar}$ と置く。:term[伝播核]{en="propagator kernel" yomi="でんぱかく"} は $K(x_b,T;x_a,0)=\mel{x_b}{U(T)}{x_a}$ である。
:::

位置の単位分解から、波動関数は $\psi(x_b,T)=\int K(x_b,T;x_a,0)\psi(x_a,0)\dd x_a$ になる。従って一次元の核の単位は長さの逆数である。$U(0)=I$ および $i\hbar\,\partial_T U=HU$ は $K(0)=\delta(x_b-x_a)$ と $i\hbar\,\partial_TK=H_{x_b}K$ を与える。また $U(T_1+T_2)=U(T_2)U(T_1)$ から

$$
K(x_b,T_1+T_2;x_a,0)
=\int K(x_b,T_2;x,0)K(x,T_1;x_a,0)\dd x .
\label{eq-aqm-kernel-composition}
$$

随伴を取ると $K(x_b,-T;x_a,0)=\overline{K(x_a,T;x_b,0)}$ である。遅延 Green 関数は $G_{\mathrm R}(T)=-i\mathbf1_{T>0}K(T)/\hbar$ であり、$(i\hbar\partial_T-H)G_{\mathrm R}=\delta(T)I$ となる。時間 Fourier 変換 $\int G_{\mathrm R}(T)e^{i(E+i\varepsilon)T/\hbar}\dd T$ は $(E+i\varepsilon-H)^{-1}$ であり、前章の外向き resolvent とつながる。

:::example[離散スペクトル]
固有関数 $\phi_n$ が完全なら $K=\sum_n\phi_n(x_b)\overline{\phi_n(x_a)}e^{-iE_nT/\hbar}$ である。この式を $\psi=\phi_m$ と積分すれば、直交性で第 $m$ 項だけが残る。
:::

:::example[経路の加算]
二次元の $H=\hbar\Omega\sigma_x/2$ で $K_{00}(T)=\cos(\Omega T/2)$、$K_{10}(T)=-i\sin(\Omega T/2)$ である。二つの等しい時間幅を $\Omega T/2=\pi/4$ に取ると、途中の $0$ を通る振幅 $1/2$ と $1$ を通る振幅 $-1/2$ が加算されて零になる。途中の位置や状態を実際に測定した場合は、この二振幅を足す実験とは異なる。
:::

### 半直線の境界を像で課す

半直線 $x>0$ の自由粒子に Dirichlet 境界を課す。正の半直線の初期関数を負側へ奇に延長し、全直線の自由核 $K_0(x_b-x_a;T)=\int e^{[ip(x_b-x_a)-ip^2T/(2m)]/\hbar}\dd p/(2\pi\hbar)$ で伝播させる。ここで $m$ は粒子質量であり、この表示は運動量基底で自由時間発展を対角化して得られる。自由 Hamiltonian は反転と可換なので奇性が保たれ、$x=0$ で波動関数は零になる。負側の積分の変数を反転して
$$
K_D(x_b,T;x_a,0)=K_0(x_b-x_a;T)-K_0(x_b+x_a;T)
$$
を得る。初期時刻では第二項のデルタの台が開半直線の外にあり、第一項だけが初期条件を与える。Neumann 境界なら偶延長を使って符号が和に替わる。強連続な半直線のユニタリ性は、奇延長に $1/\sqrt2$ を掛けた写像が全直線の奇部分空間へのユニタリであることから従う。

:::example[壁での反射位相を取り出す]
$x_a=x_b=a>0$ とすると像の経路差は $2a$、作用差は $2ma^2/T$ である。これが $\pi\hbar$ のとき、自由核に対する Dirichlet 核の比は $1-e^{i\pi}=2$、Neumann 核の比は $1+e^{i\pi}=0$ となる。これは点の核の比であり、実際の計数は幅を指定した初期波束と検出器に積分した結果である。虚時間では同じ比の指数が実の負になり、壁へ近づくほど Dirichlet 核は小さくなる。
:::

### 時刻に依存する核と相互作用表示

Hamiltonian $H(t)$ が互いに異なる時刻で可換とは限らない場合、$U(t_b,t_a)$ は $i\hbar\partial_{t_b}U=H(t_b)U$ と $U(t_a,t_a)=I$ で定める。途中の $t$ で $U(t_b,t)U(t,t_a)$ を作ると、$t_b$ の同じ微分方程式を満たし $t_b=t$ で $U(t,t_a)$ になるので、$U(t_b,t_a)$ と等しい。位置の分解を挿入すれば時間非一様な場合にも核の合成則が成り立つ。逆向きは $U(t_a,t_b)=U(t_b,t_a)^\dagger$ である。有限次元で $H(t)$ が連続なら Picard 反復で解は一意であり、有界でノルム連続な $H(t)$ の場合は作用素ノルムでの反復も収束する。非有界の場合は、これらの式を満たす伝播演算子が存在する条件を別に課す。

$H(t)=H_0+V(t)$ とし $U=U_0U_I$ を入れる。$i\hbar\dot U_0=H_0U_0$ を引けば $i\hbar\dot U_I=V_I(t)U_I$、$V_I(t)=U_0(t,t_a)^\dagger V(t)U_0(t,t_a)$ である。初期値の積分式を反復し、$n$ 次の項は $t_a<t_1<\cdots<t_n<t_b$ の領域で $(-i/\hbar)^n V_I(t_n)\cdots V_I(t_1)$ を積分したものになる。有界な $V$ なら各項のノルムは $[\sup\|V\|(t_b-t_a)/\hbar]^n/n!$ 以下なので絶対収束する。$V(t,x)$ が位置で対角なら一次の核は
$$
-\frac{i}{\hbar}\int_{t_a}^{t_b}\dd t\int\dd x\,
K_0(x_b,t_b;x,t)V(t,x)K_0(x,t;x_a,t_a)
$$
となる。短時間の点 $\vb x$ で相互作用を挿入する意味を、固有状態の摂動から直接に読むことができる。

:::example[二つのパルスの時間順序]
二準位で $\pi/2$ の $x$ 回転 $U_x=(I-i\sigma_x)/\sqrt2$ と同じ角の $z$ 回転 $U_z=(I-i\sigma_z)/\sqrt2$ を用いる。$\ket{+z}$ から時刻順に $x,z$ と掛けると、$x$ 回転の下成分の相対位相 $-i$ が次の $z$ 回転で $+1$ になり、終状態は $\ket{+x}$ である。$z,x$ の順なら最初は全体位相だけで、最後は $(1,-i)/\sqrt2$ となる。$\ket{+x}$ に射影した確率は順に $1,1/2$ である。二パルスの積の順序を省略して一つの指数にまとめられるかは、測定量にも関係する。
:::

### 核の合成を異なる時間幅で直接行う

自由粒子の二つの時間幅 $T_1,T_2>0$ と、途中の座標 $x$ を考える。二段の位相の二次式は
$$
\frac{(x-x_a)^2}{T_1}+\frac{(x_b-x)^2}{T_2}
=\left(\frac1{T_1}+\frac1{T_2}\right)(x-\bar x)^2
+\frac{(x_b-x_a)^2}{T_1+T_2},
\quad
\bar x=\frac{T_2x_a+T_1x_b}{T_1+T_2}.
$$
$\bar x$ は全経路のうち時刻 $T_1$ での直線の値である。中間の Gaussian を終えると、二つの前因子の積に $\sqrt{2\pi i\hbar/[m(1/T_1+1/T_2)]}$ が掛かり、一つの時間 $T_1+T_2$ の前因子へ戻る。繰り返すと自由核の一様でない分割でも同じであり、時間分割で収束が必要なのは運動とポテンシャルの非可換な積の部分である。

:::example[中間点と分割の幅を混同しない]
$T_2=2T_1$、$x_a=0,x_b=3d$ なら停留中間点は $d$ である。二つの速度は $d/T_1$ と $2d/(2T_1)$ で同じになる。この点の作用は $3md^2/(2T_1)$ である。単純に座標の平均 $3d/2$ を使うと、上の完成平方の正の項が余って $3md^2/(16T_1)$ だけ作用が増える。停留条件は区間ごとの時間を含めた運動量の一致である。
:::

## 時間分割による経路積分の導出

$H=p^2/(2m)+V(x)$、$T=N\varepsilon$ とする。$V$ が実有界なら Trotter 積は強収束し、運動項の微分が定義される稠密集合で一ステップの差は $o(\varepsilon)$ である。具体的に $\psi$ が $H_0$ の領域にあれば両積を $\psi$ に作用させた一次の微分は $-i(H_0+V)\psi/\hbar$ に一致する。有限時間の $e^{-iHt/\hbar}\psi$ は領域のグラフノルムでコンパクトな軌道をなし、差の $o(\varepsilon)$ はその上で一様である。これは各差分商がグラフノルムで一様に有界であることから従う。積の差を $N$ 項に望遠鏡展開し各因子のユニタリ性を使うと、全誤差は $N o(\varepsilon)\to0$ となる。非有界の振動子は以下で核の直接検証を行う。

各区間で位置の単位分解を挿入し、運動項には $\braket{x}{p}=e^{ipx/\hbar}/\sqrt{2\pi\hbar}$ を使う。各 $p_j$ には指数 $[p_j(x_{j+1}-x_j)-\varepsilon p_j^2/(2m)]i/\hbar$ が掛かる。平方完成し、収束因子を付けた Gaussian 積分を極限にすると

$$
\begin{aligned}
K&=\lim_{N\to\infty}\left(\frac{m}{2\pi i\hbar\varepsilon}\right)^{N/2}
 \int_{\R^{N-1}}\prod_{j=1}^{N-1}\dd x_j\,e^{i S_N/\hbar},\\
S_N&=\sum_{j=0}^{N-1}\varepsilon\left[
\frac{m}{2}\left(\frac{x_{j+1}-x_j}{\varepsilon}\right)^2-V(x_j)\right],
\quad x_0=x_a,\ x_N=x_b .
\end{aligned}
\label{eq-aqm-time-sliced}
$$

平方根は $\sqrt{i}=e^{i\pi/4}$ を取る。:term[経路積分]{en="path integral" yomi="けいろせきぶん"} の記号 $\int\mathcal D x\, e^{iS[x]/\hbar}$ はこのような核の極限を表す。これは古典軌道の確率密度ではない。実時間の積分は振動積分であり、極限は波束に作用させて理解する。

:::example[中点を一つ入れる]
自由粒子で $N=2$ なら $(x_1-x_a)^2+(x_b-x_1)^2=2(x_1-(x_a+x_b)/2)^2+(x_b-x_a)^2/2$ である。積分する変数を平行移動し Gaussian 積分を行うと時間 $2\varepsilon$ の一ステップの核に戻る。
:::

:::example[一定のポテンシャル]
$V=V_0$ なら $S_N$ は自由粒子の作用から $TV_0$ だけずれるため $K=e^{-iTV_0/\hbar}K_0$ となる。$V_0T=\pi\hbar$ なら核全体の符号が反転するが、この系だけの位置確率分布は変わらない。
:::

### 分割の順序と対称積

有界な $A,B$ で時間の係数を吸収して $h=\varepsilon$ と書く。順序積の $h^2$ の係数は $A^2/2+AB+B^2/2$、一回の $e^{h(A+B)}$ では中央の項が $(AB+BA)/2$ である。差は $[A,B]/2$ となる。対称積 $e^{hB/2}e^{hA}e^{hB/2}$ では両順序が各 $1/2$ ずつ現れ、二次の差も消える。作用素ノルムで Taylor 余項を評価し、時間を固定して望遠鏡和を作れば、反 Hermite な $A,B$ では順に全誤差が $\order{\varepsilon}$、$\order{\varepsilon^2}$ になる。非有界な演算子では同じ次数を保証するには、交換子や時間発展した状態の領域に追加条件が要る。

:::example[二準位で順序の差を見る]
$A=-i\sigma_x$、$B=-i\sigma_z$ とすると $[A,B]=2i\sigma_y$ である。$\ket{0}$ に対し順序積と厳密発展の差は $\varepsilon^2 i\sigma_y\ket0+\order{\varepsilon^3}=-\varepsilon^2\ket1+\order{\varepsilon^3}$ となる。順を逆にすると先頭の符号が反転する。対称積では二次まで同じなので、下向きの振幅差は三次から始まる。例えば一歩の $\varepsilon=0.1$ では順序積に約 $10^{-2}$ の振幅差が出る尺度である。
:::

### 分割の内部点に関する運動方程式

離散作用 :eqref[eq-aqm-time-sliced] を中間の $x_j$ で微分すると、隣接する二つの運動項と $V(x_j)$ だけが寄与する。従って
$$
m\frac{x_{j+1}-2x_j+x_{j-1}}{\varepsilon^2}=-V'(x_j).
$$
各区間の運動量を $p_j=m(x_{j+1}-x_j)/\varepsilon$ とすれば $p_j=p_{j-1}-\varepsilon V'(x_j)$ である。連続極限で二階の運動方程式になり、正則な軌道での規格化の離散補正は二次変分で決まる。全経路の積分の前に計算する停留点の方程式と、量子核の極限を与える強収束とは異なる段階である。非二次の $V$ では停留点を一つ求めても有限分割の積分全体は求まらない。

### 運動量積分を正則化してから極限を取る

一歩の運動量積分で $e^{-\eta p^2}$、$\eta>0$ を掛けると、$\int e^{-A p^2+Bp}\dd p=\sqrt{\pi/A}\,e^{B^2/(4A)}$（$\Re A>0$）を使える。ここで $A=\eta+i\varepsilon/(2m\hbar)$、$B=i(x_b-x_a)/\hbar$ である。複素 Gaussian の式は $\Re A>0$ で解析的で、正の実数 $A$ での完成平方と同じ枝に連続につながる。$B$ については $\partial_B$ と部分積分から $\partial_B I=B I/(2A)$ を得て確認できる。$\eta\downarrow0$ で核に必要な $\sqrt i$ の枝が出る。途中の積分を正則化した後に波束に作用させれば、自由の Fourier 乗算 $e^{-\eta p^2-i\varepsilon p^2/(2m\hbar)}$ の極限を有界収束で取れる。これにより切断なしの振動積分を普通の絶対収束積分と同じ扱いにすることを避けられる。

:::example[一歩で位相に効くポテンシャルの尺度]
一歩で一定差 $\Delta V$ を持つ二つの区間を比較すると相対位相は $-\varepsilon\Delta V/\hbar$ である。全体の時間を同じ区間 $N$ 個に細分しても、差が各区間にあれば相対位相は $-T\Delta V/\hbar$ のままである。例えば $\Delta V=1\,\mathrm{meV}$、$T=0.658\,\mathrm{ps}$ では $\hbar\simeq0.658\,\mathrm{meV\,ps}$ より約 $-1$ rad である。一歩の位相が小さくなることと、区間を通じた干渉が小さくなることは区別する。
:::

### 位相空間積分で端点に現れる境界項

運動量と座標の短時間の重なりを使うと、運動量ごとの指数は $p_j(x_{j+1}-x_j)-\varepsilon [p_j^2/(2m)+V(x_j)]$ である。離散的な部分積分は
$$
\sum_{j=0}^{N-1}p_j(x_{j+1}-x_j)
=p_{N-1}x_N-p_0x_0-\sum_{j=1}^{N-1}x_j(p_j-p_{j-1})
$$
である。左では端の位置が固定、右ではもし端の運動量を指定したいなら境界の $p x$ の位相を Fourier 変換で処理する。端点の作用を替えずに運動量表示の核と同じと思うと、この変換の符号を誤る。実際、自由粒子の位置核の両端を Fourier 変換すると、先に $x_a,x_b$ との位相を積分した結果、同じ運動量だけが残り $\delta(p_b-p_a)e^{-ip_a^2T/(2m\hbar)}$ になる。これは運動量固定の経路で作用の位相が $-\int H\dd t$ になることと一致する。

## 自由粒子と調和振動子の経路積分

### 自由核

運動量積分を直接行うと

$$
K_0=\sqrt{\frac{m}{2\pi i\hbar T}}\exp\left(\frac{im(x_b-x_a)^2}{2\hbar T}\right).
\label{eq-aqm-free-kernel}
$$

作用は一定速度の軌道の $m\dot x^2 T/2$ に等しい。収束因子付きの Fourier 積分で $T\downarrow0$ を取るとデルタ関数に収束する。三次元では直交座標の運動量積分が三つに分かれ、前因子が $3/2$ 乗となる。

### 振動子の有限分割

$V=m\omega^2 x^2/2$ とし、$0<T<\pi/\omega$ とする。古典方程式 $\ddot x=-\omega^2x$ に両端条件を入れると

$$
x_{\mathrm{cl}}(t)=\frac{x_a\sin(\omega(T-t))+x_b\sin(\omega t)}{\sin(\omega T)}.
$$

運動方程式上で $m(\dot x^2-\omega^2 x^2)/2=(m/2)\dv{(x\dot x)}t$ なので

$$
S_{\mathrm{cl}}=\frac{m\omega}{2\sin(\omega T)}
 [(x_a^2+x_b^2)\cos(\omega T)-2x_ax_b].
$$

離散積分の内部 Hessian は $m/\varepsilon$ を掛けた三重対角行列で、主対角 $2-\omega^2\varepsilon^2$、副対角 $-1$ である。$d_n$ をこの無次元行列の $n$ 次主小行列式とすると最終行展開から $d_n=(2-\omega^2\varepsilon^2)d_{n-1}-d_{n-2}$、$d_0=1,d_1=2-\omega^2\varepsilon^2$ である。$2\cos\vartheta=2-\omega^2\varepsilon^2$ と置けば $d_n=\sin[(n+1)\vartheta]/\sin\vartheta$ が漸化式を満たす。$N-1$ 個の積分の前因子は $(m/(2\pi i\hbar\varepsilon d_{N-1}))^{1/2}$ で、$\varepsilon d_{N-1}\to\sin(\omega T)/\omega$ である。離散停留点も $x_{j+1}+x_{j-1}=2\cos\vartheta\,x_j$ を満たすので正弦補間に収束する。従って

$$
K_{\mathrm{ho}}=\sqrt{\frac{m\omega}{2\pi i\hbar\sin\omega T}}\ e^{i S_{\mathrm{cl}}/\hbar}.
\label{eq-aqm-oscillator-kernel}
$$

この式は $i\hbar\partial_TK=(-\hbar^2\partial_{x_b}^2/(2m)+m\omega^2x_b^2/2)K$ も満たす。指数の実部の式は Hamilton–Jacobi 方程式、前因子の微分は $\partial_T\log(\sin\omega T)^{-1/2}=-(\omega/2)\cot\omega T$ に対応する。$T\to0$ では自由核に近づくので初期条件が揃う。

:::example[四分の一周期]
$\omega T=\pi/2$ では $S_{\mathrm{cl}}=-m\omega x_ax_b$ となる。従って核は尺度 $m\omega/\hbar$ の Fourier 変換を全位相 $e^{-i\pi/4}$ とともに行う。
:::

:::example[半周期]
スペクトル $E_n=(n+1/2)\hbar\omega$ とパリティ $(-1)^n$ より $U(\pi/\omega)=e^{-i\pi/2}\Pi$ である。そこで核は $e^{-i\pi/2}\delta(x_b+x_a)$ となる。半周期を越えるときは短い時間区間の核を合成して平方根の枝を選ぶ。
:::

### 自由な Gaussian 波束の広がり

初期を $\psi(x,0)=(2\pi s^2)^{-1/4}e^{-(x-a)^2/(4s^2)+ip_0(x-a)/\hbar}$ とする。自由運動の Heisenberg 位置は $x(t)=x+pt/m$ である。$[p^2,x]=-2i\hbar p$ で一回の交換子以降が消えるからである。初期では $\langle p\rangle=p_0$、$\operatorname{Var}(p)=\hbar^2/(4s^2)$ であり、$\langle\{x-a,p-p_0\}\rangle=0$ となる。最後の二つは $(p-p_0)\psi=i\hbar(x-a)\psi/(2s^2)$ を使う。従って平均と分散は
$$
\langle x(t)\rangle=a+\frac{p_0t}{m},\qquad
s_t^2=s^2+\frac{\hbar^2t^2}{4m^2s^2}.
$$
核の積分で直接確認するには $\xi=x_a-a$、$y=x_b-a-p_0t/m$ とし、指数の $\xi$ についての二次項 $-(1-i\,2ms^2/(\hbar t))\xi^2/(4s^2)$ と一次項を平方完成する。結果は初期の規格化係数に $(1+i\hbar t/(2ms^2))^{-1/2}$ を掛け、幅の因子 $e^{-y^2/[4s^2(1+i\hbar t/(2ms^2))]}$ と平面波の位相を掛けたものになる。絶対値二乗で上の $s_t^2$ を得る。

:::example[電子の波束が広がる時刻]
$s=1\,\mathrm{nm}$、電子質量 $m=9.11\times10^{-31}\,\mathrm{kg}$ とする。$t=2ms^2/\hbar\simeq1.73\times10^{-14}\,\mathrm s$ では分散が二倍で、標準偏差は $\sqrt2\,\mathrm{nm}$ となる。$p_0=0$ なら平均位置は動かない。この広がりは初期の運動量分散から生じ、同じ質量で $s$ を十倍にすると分散が二倍になる時刻は百倍になる。
:::

### 振動子の負固有値を数えて枝を選ぶ

時間分割で内部の $j=1,\ldots,N-1$ に Dirichlet 条件を課した行列 $A$ の固有ベクトルに $u^{(s)}_j=\sin(\pi sj/N)$ を入れる。両隣との和は $2\cos(\pi s/N)u^{(s)}_j$ であり、固有値は
$$
\lambda_s=4\sin^2\left(\frac{\pi s}{2N}\right)-\omega^2\varepsilon^2,\qquad s=1,\ldots,N-1
$$
である。これらは異なるので全内部次元を張る。$N\varepsilon=T$ を固定して $N$ を増やすと、負固有値の条件は先頭で $\pi s<\omega T$ になる。$\omega T$ が $\pi$ の整数倍でなければ、負の個数は $\nu=\lfloor\omega T/\pi\rfloor$ に揃う。従って絶対値の $\sin\omega T$ を使った前因子に $e^{-i\pi\nu/2}$ を掛ける。時間が半周期を越えると $\sin$ の符号だけを見て平方根を選ぶ方法では、全周期を越えたときの履歴を保持できない。

:::example[四分の三周期で Fourier 変換の向きが替わる]
$\omega T=3\pi/2$ では $\nu=1$、$S_{\rm cl}=+m\omega x_ax_b$ である。四分の一周期の核と比べて Fourier 指数の符号が逆、前因子の位相はさらに $e^{-i\pi/2}$ 倍である。実際、$U(3\pi/(2\omega))=e^{-i\pi/2}\Pi U(\pi/(2\omega))$ によっても同じである。全周期の $U(2\pi/\omega)=-I$ では分布核は $-\delta(x_b-x_a)$ となり、絶対値発散の式へ直接代入する必要はない。全周期の位置分布は初期と同じでも、経路の片側だけにこの変換を掛ける比較では全位相が残る。
:::

### 固有関数と Gaussian の核を照合する

振動子の自然な尺度を $\ell=\sqrt{\hbar/(m\omega)}$ とする。$0<\omega T<\pi$ の核を基底 $\phi_0(x)=(\pi\ell^2)^{-1/4}e^{-x^2/(2\ell^2)}$ に作用させる。$x_a^2$ の係数は $-A/(2\ell^2)$、$A=1-i\cot\omega T$ で、一次の係数は $-ix_b/(\ell^2\sin\omega T)$ である。平方完成で出る $x_b^2$ の係数は
$$
\frac{1}{2\ell^2}\left[i\cot\omega T-\frac{1}{A\sin^2\omega T}\right]=-\frac1{2\ell^2}.
$$
$A\sin\omega T=\sin\omega T-i\cos\omega T=-i e^{i\omega T}$ を使った。前因子は同じ枝で $(iA\sin\omega T)^{-1/2}=e^{-i\omega T/2}$ に簡約される。従って波形はそのままで零点エネルギーの位相だけが残る。この照合は、作用の符号、規格化、平方根の三つを一つの初期波束で確かめる方法である。

:::example[第一励起の位相も同時に確かめる]
第一励起は $\phi_1(x)=\sqrt2 (x/\ell)\phi_0(x)$ である。上の Gaussian 積分に $x_a$ を挿入すると、線形係数に関する一回の微分となり、平均 $-ix_b/(A\sin\omega T)=x_b e^{-i\omega T}$ が掛かる。したがって位相は $e^{-3i\omega T/2}$ となる。$\omega T=\pi/2$ で $\phi_0+\phi_1$ の規格化した重ね合わせを準備すると、相対位相は $-i$ になる。$\phi_0,\phi_1$ を実に選ぶと、その時点の位置密度の交差項は純虚なので消える。
:::

### 二次系で異なる周波数を持つ方向を分ける

二次元で $V=m(\omega_x^2 x^2+\omega_y^2y^2)/2$ とする。Hamiltonian が各座標の和でその二作用素が可換なため、核は二つの一次元核の積である。停留作用の混合微分行列も対角であり、対角要素は $m\omega_x/\sin\omega_xT$ と $m\omega_y/\sin\omega_yT$ である。焦点の時刻は二方向で別々に現れ、ある方向だけがデルタ核になるとき、他方は位置間の Fourier 型の振幅を保つ。二次作用の空間変数を実直交変換して固有周波数に分ける場合も、測度が変わらないので同じ積の構造を使える。

:::example[片方向だけが反転する時刻]
$\omega_y=2\omega_x$ とし、$T=\pi/(2\omega_x)$ を選ぶ。$x$ は四分の一周期で $S_x=-m\omega_x x_ax_b$、$y$ は半周期で $K_y=-i\delta(y_b+y_a)$ になる。積の全位相は $x$ の平方根の $e^{-i\pi/4}$ と $y$ の $-i$ の積で $e^{-3i\pi/4}$ である。$y$ を初期と同じ位置で測る確率は初期波束が原点近くにどれだけあるかに依存し、$x$ と同じ連続な点の核として評価しない。
:::

## 停留位相近似と半古典極限

振動積分 $\int a(x)e^{iS(x)/\hbar}\dd x$ で $S'$ が零でない領域では、$\partial_x e^{iS/\hbar}=(iS'/\hbar)e^{iS/\hbar}$ を使った部分積分で次数を下げられる。滑らかな実位相 $S$ と滑らかなコンパクト台の振幅を取り、停留点が有限個でいずれも非退化とする。孤立した各停留点 $x_0$ では、$x-x_0$ を $\sqrt{\hbar}$ で拡大すると二次項が有限に残る。符号つき平方根の局所変数 $y=\sgn(x-x_0)\sqrt{2|S(x)-S(x_0)|}$ に替えれば、位相は正確に二次になる。その Jacobian は $1/\sqrt{|S''(x_0)|}$ である。局所振幅を定数と残差に分けると、定数が Fresnel 積分を与え、残差の奇項が消え、偶項の最初は $y^2$ なので一つ高い次数となる。各点の寄与を足すと先頭項は
$$
\sum_{x_0}
a(x_0)\sqrt{\frac{2\pi\hbar}{|S''(x_0)|}}\,
e^{iS(x_0)/\hbar+i\pi\sgn(S''(x_0))/4}
$$
で、余りは $\order{\hbar^{3/2}}$ である。これが :term[停留位相近似]{en="stationary phase approximation" yomi="ていりゅういそうきんじ"} である。

経路の停留条件は端点固定の Euler–Lagrange 方程式 :ref[thm-am-euler-lagrange] である。分割作用での Gaussian 行列式を離散 Jacobi 方程式の解と比べると、一次元では端点の混合微分が前因子になる。具体的に離散停留点で $x_0$ を固定し、$x_1$ を変えて $Y_j=\partial x_j/\partial x_1$ とする。漸化式は $Y_{j+1}=(2-\varepsilon^2 V''(x_j)/m)Y_j-Y_{j-1}$、$Y_0=0,Y_1=1$ であり $Y_N$ は無次元 Hessian の行列式に一致する。初期運動量は $p_a=m(x_1-x_0)/\varepsilon+\varepsilon V'(x_0)$ なので $\partial p_a/\partial x_b=m/(\varepsilon Y_N)$ である。一方 $\partial_{x_a}S_{\mathrm{cl}}=-p_a$ であるから、連続極限では

$$
K_{\mathrm{sc}}\simeq \sum_{\rm paths}
\frac{1}{\sqrt{2\pi i\hbar}}
\left|-\pdv{}{x_b}\pdv{S_{\mathrm{cl}}}{x_a}\right|^{1/2}
 e^{iS_{\mathrm{cl}}/\hbar-i\pi\nu/2}.
$$

$\nu$ は時間分割の Hessian の負固有値数で、各負方向は正方向に比べ Fresnel 因子を $e^{-i\pi/2}$ だけ替える。混合微分が発散する焦線では二次の近似を使えず、一様近似が必要になる。この連続表示は非退化な古典軌道が孤立する領域で用いる。

:::example[自由粒子と振動子の前因子]
自由粒子では $-\partial_{x_b}\partial_{x_a}S=m/T$、振動子では $m\omega/\sin\omega T$ である。$0<\omega T<\pi$ の両者の Hessian は正で、先の厳密な核に一致する。
:::

:::example[二経路の作用差]
同じ絶対値の前因子を持つ二経路で $\nu$ の差が零なら、強度は一経路の $2[1+\cos((S_2-S_1)/\hbar)]$ 倍である。$S_2-S_1=\pi\hbar$ のとき二経路の振幅が打ち消し合う。
:::

### 線形ポテンシャルで停留軌道を検算する

$V=-Fx$ として両端を固定する。$x_{\rm cl}=x_a+(x_b-x_a)t/T+Ft(t-T)/(2m)$ は $m\ddot x=F$ を満たす。作用に入れると直線部分との運動の交差積分は $\int_0^T(2t-T)\dd t=0$ で消える。二次の速度が $F^2T^3/(24m)$、$F x_{\rm cl}$ の曲がりから $-F^2T^3/(12m)$ が出る。従って
$$
S_{\rm cl}=\frac{m(x_b-x_a)^2}{2T}+\frac{FT(x_a+x_b)}2-\frac{F^2T^3}{24m}.
$$
二次変分は自由のものと同じなので前因子は :eqref[eq-aqm-free-kernel] と同じである。線形ポテンシャルの寄与によって端点の和と伝播時間に依存する位相が増えるが、$x_a,x_b$ の混合微分には影響しない。

:::example[端の運動量の差を読む]
$x_b=x_a$ の例では、$-\partial_{x_a}S_{\rm cl}=-FT/2$、$\partial_{x_b}S_{\rm cl}=FT/2$ である。初速度は負で最後に正となり、運動量の増分は $FT$ である。$F$ の符号を替えると両端の運動量も符号を替える一方、$F^2$ の位相項は変わらない。
:::

### Gaussian モーメントで停留位相の次の項を求める

位相を厳密に $S=S_0+\lambda y^2/2$ に直し、実振幅を $b(y)=b(0)+b'(0)y+b''(0)y^2/2+\cdots$ とする。停留点の近くに偶の切断を入れれば奇項は消え、外では位相の微分で部分積分して任意の次数で小さくできる。規格化した Fresnel の二乗モーメントは $i\hbar/\lambda$ である。実際、$I(\lambda)=\sqrt{2\pi\hbar/|\lambda|}e^{i\pi\,\sgn\lambda/4}$ を $\lambda$ で微分すると、被積分関数に $iy^2/(2\hbar)$ が掛かる。従って次の係数は $i\hbar b''(0)/(2\lambda)$ である。$\lambda$ が零に近いと補正に逆冪が現れるので、単純な非退化停留点の近似は一様でなくなる。

:::example[同じ二次位相で奇振幅と偶振幅を入れる]
正則化 $e^{-\eta y^2}$ を付けて振幅 $1+a y+b y^2$ を積分する。$a$ の項は任意の $\eta>0$ で奇関数なので零である。$\eta\to0$ で全体は未変形の Gaussian 因子に $1+ib\hbar/\lambda$ を掛けた値になる。例えば $b\hbar/\lambda=1/5$ なら係数の絶対値二乗は $26/25$、位相は $\arctan(1/5)$ である。補正の一次が純虚なため、この例では強度の変化はその二乗から始まる。
:::

### 停留点の合流を三次の核で表す

二つの停留点が合流する最も単純な模型では位相を $x^3/3-a x$ とする。$a>0$ で停留点は $\pm\sqrt a$、作用差は $4a^{3/2}/3$、二階微分は $\pm2\sqrt a$ である。個別の停留位相を足すと前因子は $a^{-1/4}$ の尺度になり、$a\to0$ まで一様には使えない。一方 $x=\hbar^{1/3}y$ に替えて積分全体を保つと
$$
\int_{\R}e^{i(x^3/3-a x)/\hbar}\dd x
=2\pi\hbar^{1/3}\operatorname{Ai}(-a/\hbar^{2/3}).
$$
右辺の Airy 関数は $\operatorname{Ai}(z)=(2\pi)^{-1}\int e^{i(y^3/3+zy)}\dd y$ を振動積分として定義したものである。だから $a$ を変えるときの一様な尺度は $a/\hbar^{2/3}$ になる。この導入は三次模型の厳密な変数変換であり、一般の焦線が同じ模型で近似できることの証明には局所の変数変換を別に調べる。

:::example[二つの位相を分離できる尺度]
$a=9\hbar^{2/3}$ なら二点の作用差は $36\hbar$ であり、停留点周りの幅 $\sqrt{\hbar/(2\sqrt a)}$ は点間距離 $2\sqrt a$ より小さい。一方 $a=\hbar^{2/3}/4$ なら作用差は $\hbar/6$ で、それぞれの位相を独立に移動する幅の重なりが無視できなくなる。この場合に二つを個別に Gaussian へ置き換えると、合流近くで誤った発散を作る。
:::

### 四次の揺らぎが Gaussian に与える補正

有限次元の位相の一変数で $S(x)=\lambda x^2/2+\zeta x^4/4$、$\lambda>0,\zeta\ge0$ とする。原点が唯一の停留点であり、局所に 1 の滑らかな切断を入れる。$x=\sqrt{\hbar/\lambda}\,y$ で二次位相を正規化すると、四次の指数は $i\zeta\hbar y^4/(4\lambda^2)$ である。二次の正則化積分を $\lambda$ で二回微分すると Fresnel の四乗モーメントは $-3\hbar^2/\lambda^2$ になる。従って指数の一次展開による全体の係数は $1-3i\zeta\hbar/(4\lambda^2)+\order{\hbar^2}$ である。局所で位相を正確に二次へ写す方法でも、変数変換 $y=x\sqrt{1+\zeta x^2/(2\lambda)}$ の Jacobian は $1-3\zeta y^2/(4\lambda)+\cdots$ となり、二乗モーメントを掛けて同じ係数を得る。切断の外の位相は導関数を持つため部分積分で高次にできる。

:::example[四次の作用を持つ一分割の中間点]
同じ幅 $\varepsilon$ 二つ、端点がともに零で、中間点に $V(x)=-b x^4/4$ の局所近似を入れる。作用は $mx^2/\varepsilon+\varepsilon b x^4/4$ であるから $\lambda=2m/\varepsilon,\zeta=\varepsilon b$ である。補正の無次元係数は $-3i\,b\hbar\varepsilon^3/(16m^2)$ となる。例えば $b\hbar\varepsilon^3/m^2=1/10$ なら絶対値の尺度は $3/160$ である。この例は一つの有限分割の局所積分であり、非有界な負の四次ポテンシャルの連続時間の核を一般に構成しているわけではない。
:::

## 虚時間経路積分と分配関数

$T=-i\tau$ とすると、下に有界な $H$ のスペクトル表示から $U(-i\tau)=e^{-\tau H/\hbar}$ となる。時間分割で各 Gaussian は正の実係数になる。実の $V$ に対して重みは $e^{-S_E/\hbar}$ であり

$$
S_E[x]=\int_0^\tau\left(\frac{m}{2}\dot x^2+V(x)\right)\dd t_E .
\label{eq-aqm-euclidean-action}
$$

これを :term[Euclid 作用]{en="Euclidean action" yomi="ゆーくりっどさよう"} と呼ぶ。熱平衡の逆温度 $\beta=1/(k_BT_{\mathrm{th}})$ に対し $\tau=\beta\hbar$ と取る。トレースのために $x_b=x_a$ としその共通値を積分すれば、周期経路の積分で $Z=\tr e^{-\beta H}$ を得る。無限空間の自由粒子ではこのトレースが発散するため、体積を有限にするか単位体積当たりの値を使う。

:::example[振動子の分配関数]
$u=\beta\hbar\omega$ と置く。:eqref[eq-aqm-oscillator-kernel] の解析接続で対角核は $\sqrt{m\omega/(2\pi\hbar\sinh u)}e^{-(m\omega/\hbar)x^2\tanh(u/2)}$ となる。積分すると $Z=1/[2\sinh(u/2)]$ である。$u\gg1$ では $e^{-u/2}$ が支配し、零点エネルギー $\hbar\omega/2$ が残る。
:::

:::example[自由粒子の熱波長]
自由核の対角値は $\sqrt{m/(2\pi\hbar^2\beta)}$ である。辺長 $L$ の周期箱で $L$ が熱波長 $\lambda_T=\sqrt{2\pi\hbar^2\beta/m}$ より十分長ければ、非零の巻き数の寄与は $\exp(-\pi L^2/\lambda_T^2)$ 以下の尺度で小さく、$Z\simeq L/\lambda_T$ となる。
:::

### 熱核を自由な橋の期待値にする

実有界 $V$ について、虚時間分割の各 Gaussian の積から自由核 $K_{0,E}$ を括り出す。残りは両端を固定した多変量 Gaussian の確率密度である。実際、区間ごとの増分を独立な平均零、分散 $\hbar\varepsilon/m$ の正規変数とし、全増分を $x_b-x_a$ に条件付けた密度である。中間位置の平均は $x_a+(x_b-x_a)t/\tau$、共分散は
$$
\frac{\hbar}{m}\big(\min(t,s)-ts/\tau\big).
$$
独立増分の累積の共分散は $(\hbar/m)\min(t,s)$ であり、終点を条件付けると両時刻と終点の共分散の積を終点の分散で割った値を引く。連続な $V$ では経路上の Riemann 和が収束し、有界性で期待値の極限を通せる。従って $K_E=K_{0,E}\mathbb E_{a\to b}[\exp(-\int_0^\tau V(x(t))\dd t/\hbar)]$ となる。これが、この仮定の下での :term[Feynman–Kac 公式]{en="Feynman–Kac formula" yomi="ふぁいんまんかっくこうしき"} である。この表示は標本経路の微分可能性を仮定しない。時間分割の運動項は増分の Gaussian 測度として既に組み込まれており、橋の期待値で積分するのは $V$ の項である。

:::example[橋の中央の分散でポテンシャルを一次まで測る]
$x_a=x_b=0$ の橋で、二次ポテンシャル $V=m\omega^2x^2/2$ の零結合まわりの一次補正を考える。$t=\tau/2$ で位置の分散は $\hbar\tau/(4m)$ である。期待値を全時刻で積分すると $\int_0^\tau\langle V\rangle\dd t/\hbar=\omega^2\tau^2/12$ を得る。非負の二次ポテンシャルでは零結合からの右微分を Gaussian の有限な二次モーメントで支配できる。厳密核と自由核の比 $\sqrt{\omega\tau/\sinh(\omega\tau)}=1-\omega^2\tau^2/12+\cdots$ と一致する。
:::

### 分配関数と虚時間での相関

$\beta$ を固定し $0<\tau<\beta\hbar$ とする。演算子 $X$ の熱相関は $\tr[e^{-(\beta-\tau/\hbar)H}X e^{-\tau H/\hbar}X]/Z$ である。二つの区間にそれぞれ虚時間核を使い、$X=x$ なら間の位置と周期端で位置を因子として挿入する。振動子で $X=\ell(a+a^\dagger)/\sqrt2$ とすると、$e^{\tau H/\hbar}a e^{-\tau H/\hbar}=e^{-\omega\tau}a$ で、反対の演算子 $a^\dagger$ は $e^{\omega\tau}$ を得る。Gibbs 状態では占有を二つ変える期待値が零なので
$$
C_E(\tau)=\frac{\ell^2}{2}\big[(\bar n+1)e^{-\omega\tau}+\bar n e^{\omega\tau}\big],
\quad \bar n=(e^{\beta\hbar\omega}-1)^{-1}.
$$
$\tau\mapsto\beta\hbar-\tau$ で二項が入れ替わる。低温で $\tau$ を固定すれば $\ell^2e^{-\omega\tau}/2$ となり、虚時間減衰の率から励起ギャップ $\hbar\omega$ を読める。$\tau$ を熱周期の半分に近づける場合は第二項も残す。

:::example[熱周期の中央と端を比べる]
$\beta\hbar\omega=\log4$ では $\bar n=1/3$ である。端の $C_E(0)$ は $5\ell^2/6$、中央 $\tau=\beta\hbar/2$ では二項が等しくそれぞれ $2/3$ なので $C_E=2\ell^2/3$ である。比は $4/5$ になる。中央で後者の「増える指数」を捨ててしまうと半分の値になる。これは周期端に近いもう片側からも同じ状態が伝わるためである。
:::

### 熱核から高温の古典分配関数を得る

有界で連続な $V$ と有限な周期箱で、$\beta\downarrow0$ を考える。同じ点へ戻る自由な橋の幅は典型的に $\sqrt{\hbar^2\beta/m}$ である。各点 $x$ に対し短い橋上の $V$ は連続性で $V(x)$ へ近づく。より具体的には橋の時刻の分布を $s=t/\tau$ で固定し、$x+\sqrt{\hbar\tau/m}\,B_s$ と書く。連続な標準橋 $B_s$ は有限区間で有界なので、$\int_0^1V(x+\sqrt{\hbar\tau/m}B_s)\dd s\to V(x)$ である。有界な $V$ を用い、指数を一次まで展開して期待値へ極限を通すと
$$
K_E(x,x;\beta\hbar)
=\sqrt{\frac{m}{2\pi\beta\hbar^2}}\,e^{-\beta V(x)}[1+o(\beta)]
$$
となる。周期箱の巻き経路の補正は短時間で指数的に小さい。したがって $Z$ の先頭は $\int\dd x\,\dd p\,e^{-\beta[p^2/(2m)+V(x)]}/(2\pi\hbar)$ である。非有界な $V$ を無限空間で扱うときは、対角核の点ごとの式だけからトレースへ極限を通せるとは限らない。調和振動子は厳密な $Z$ を持つため直接に照合できる。

:::example[振動子で古典の先頭と最初の補正を比較する]
$u=\beta\hbar\omega$ と置く。古典の積分では $Z_{\rm cl}=1/u$、正確な値は $Z=1/(2\sinh(u/2))$ である。従って $Z/Z_{\rm cl}=1-u^2/24+7u^4/5760+\cdots$ となる。例えば $u=1/2$ では先頭の補正は $-1/96$ であり、古典の値 2 に対し量子は約 1.979 になる。この比の導出には零点エネルギーをあらかじめ差し引いていない。エネルギーの原点を変えれば $Z$ は対応する指数を得る。
:::

## 電磁場中の粒子と Aharonov–Bohm 効果

SI 単位で電荷 $Q$ の Hamiltonian は $(\vb p-Q\vb A)^2/(2m)+Q\phi$ である。短時間の位相空間積分で $\vb p$ を $\vb p-Q\vb A$ に平行移動すると、運動項に加えて $Q\vb A\cdot\Delta\vb x-Q\phi\varepsilon$ が残る。$\vb A$ は中点で評価する。中点規則は Weyl 順序に対応し、$\vb p\cdot\vb A+\vb A\cdot\vb p$ の両順序を同じ重みで与える。連続表示で経路の位相は

$$
\frac{1}{\hbar}\int\left(\frac m2\dot{\vb x}^2+
 Q\,\dot{\vb x}\cdot\vb A-Q\phi\right)\dd t
$$

である。:ref[def-em-gauge] の $\vb A\mapsto\vb A+\grad\Lambda$、$\phi\mapsto\phi-\partial_t\Lambda$ では作用の変化は $Q[\Lambda_b-\Lambda_a]$ であり、核には端点位相だけが掛かる。

粒子が通る領域で磁場が零でも、閉曲線を張る面の内部に磁束があれば $Q\oint\vb A\cdot\dd\vb x/\hbar=Q\Phi/\hbar$ が残る。これが :term[Aharonov–Bohm 効果]{en="Aharonov–Bohm effect" yomi="あはらのふぼーむこうか"} の磁気型である。通過領域で大域的な実一価関数 $\Lambda$ により $\vb A+\grad\Lambda=0$ とできるなら閉曲線積分は零なので、非零の磁束があるとそのような $\Lambda$ は存在しない。ただし量子状態に作用するゲージは位相 $e^{iQ\Lambda/\hbar}$ であり、磁束が $h/Q$ の整数倍なら、この位相だけを一価に保つゲージでリング上の $\vb A$ を消せる。干渉が区別するのは $Q\Phi/\hbar$ を $2\pi$ で割った余りである。

:::example[干渉の磁束周期]
電荷 $Q$ の二経路が一周の閉曲線をなすと、磁束 $h/|Q|$ の増加で相対位相は $2\pi$ 変わる。電子一個では $h/e\simeq4.14\times10^{-15}\,\mathrm{Wb}$ である。磁束がこの半分だけ増えると二経路の干渉項の符号が変わる。
:::

:::example[磁束を囲むリング]
半径 $R$ のリングで接線成分 $A_\theta=\Phi/(2\pi R)$ と選ぶと $H=(-i\hbar\partial_\theta-Q\Phi/(2\pi))^2/(2mR^2)$ である。周期境界の $e^{in\theta}$ に作用させて $E_n=\hbar^2(n-Q\Phi/h)^2/(2mR^2)$ となる。$Q\Phi/h=1/2$ で $n=0,1$ のエネルギーが等しくなる。
:::

### 巻き数とリングの熱核

角度差を $\Delta\theta$ とし、円周を覆う直線の上で $\Delta\theta+2\pi w$ を終点に取る。整数 $w$ の各経路族には磁気作用 $Q\Phi(\Delta\theta+2\pi w)/(2\pi)$ が掛かる。慣性モーメント $I=mR^2$ として、虚時間核を角座標で正規化するなら
$$
K_E=\sqrt{\frac{I}{2\pi\hbar\tau}}\sum_{w\in\Z}
\exp\left[-\frac{I(\Delta\theta+2\pi w)^2}{2\hbar\tau}
+i\frac{Q\Phi}{h}(\Delta\theta+2\pi w)\right].
$$
$\Delta\theta$ に $2\pi$ を加えても和の添字を動かせば同じ値になる。Fourier 係数を取ると全 $w$ の積分区間が実直線を覆い、線形位相の波数は $n-Q\Phi/h$ になる。Gaussian 変換で係数は $(2\pi)^{-1}\exp[-\tau\hbar(n-Q\Phi/h)^2/(2I)]$ であり、固有値による熱核に一致する。

:::example[閉じた経路に残る磁束]
対角では $w=\pm1$ が対を作るため $w=0$ の項に対する補正は $2e^{-2\pi^2 I/(\hbar\tau)}\cos(2\pi Q\Phi/h)$ である。$\tau$ を固定して $Q\Phi/h=0,1/4,1/2$ と変えるとこの補正の係数は順に $2,0,-2$ となる。四分の一磁束では $w=\pm2$ 以降の項まで零になる必要はない。
:::

### 位相空間での中点評価と順序

一成分の位置演算子を $x$、$p=-i\hbar\partial_x$ とする。Weyl 順序の $A(x)p$ に対応する核は $\int (\dd p/2\pi\hbar) e^{ip(x_b-x_a)/\hbar} A((x_b+x_a)/2)p$ である。$p$ を微分表示に替えれば $\delta'(x_b-x_a)$ が出る。試験波動関数と $x_a$ の積を微分すると $A(x_b)\partial_{x_b}\psi +A'(x_b)\psi/2$ であり、$-i\hbar$ を掛ければ $(Ap+pA)\psi/2$ である。従って $\vb p\cdot\vb A$ の中点規則は、量子 Hamiltonian の対称順序を回復する。端点の一方だけで評価すると $A'$ の係数が $0$ または $1$ となり、そのままでは Hermite な最小結合にならない。

:::example[純ゲージを核の端に移す]
一次元で $A(x)=ax$、$\phi=0$、$\Lambda(x)=ax^2/2$ とすれば、粒子の磁場に関する局所的な力は零で、$\int Q A\dot x\dd t=Qa(x_b^2-x_a^2)/2$ である。離散でも中点なら各区間に $Qa(x_{j+1}+x_j)(x_{j+1}-x_j)/2$ が掛かり、完全に差の形で和が端だけに残る。例えば $x_a=0,x_b=d$、$Qa d^2/(2\hbar)=\pi$ なら自由核から符号が替わる。初期状態も同じゲージで変換して確率を求めるので、単一の開経路のこの符号自体を物理的な磁束と扱わない。
:::

### 電位差による干渉と端点位相

二経路で粒子がそれぞれ空間的に一定な電位 $\phi_1(t),\phi_2(t)$ の領域にいる時間を $[t_a,t_b]$ とする。切り替え中も波束が各領域内に収まり、経路の近くの電場と磁場が零で、他の区間で両経路の作用が等しいと仮定する。そのとき電位の項が作る経路 1 の経路 2 に対する位相は $-Q\int_{t_a}^{t_b}(\phi_1-\phi_2)\dd t/\hbar$ である。各経路に沿って $\vb A\cdot\dd\vb x-\phi\dd t$ を積分した差として書けば、全経路の始終点が共通なのでゲージ変換の端点項が相殺される。各領域で局所的に電位を零にするゲージを選ぶ場合、つなぎ目に同じ位相差が移る。

:::example[同じ時間だけ電位差を保つ]
二経路に電位差 $\phi_1-\phi_2=\Delta\phi$ を時間 $t$ だけ与え、残りの作用を揃える。等振幅の出力を合成する一方のポートの確率は $[1+\cos(Q t\Delta\phi/\hbar)]/2$ である。$|Q|t|\Delta\phi|=\pi\hbar$ ならそのポートは暗くなる。電位差が作用した時間を半分にすると確率は $1/2$ になる。
:::

### 環上の磁束と平衡電流

$E_n=\hbar^2(n-\alpha)^2/(2I)$、$\alpha=Q\Phi/h$ とする。各枝の持続電流は $I_n=-\partial E_n/\partial\Phi$ であり、これは電荷 $Q$ が一周する頻度を掛けたものに等しい。角速度 $\dot\theta=\hbar(n-\alpha)/I$ により $I_n=Q\hbar(n-\alpha)/(2\pi I)$ になる。熱平衡での平均は
$$
I_{\rm th}=-\partial_\Phi F
=\sum_n \frac{e^{-\beta E_n}}{\sum_m e^{-\beta E_m}}I_n,\quad
F=-\beta^{-1}\log Z
$$
である。$\alpha\mapsto\alpha+1$ では $n$ を一つずらすと全分配関数が同じになるため $F$ と平均電流は磁束周期を持つ。零温度で一つの基底枝を占有した場合は最小の $|n-\alpha|$ を選ぶので、半整数の所で枝が交替する。そこに正確に乗り Gibbs 極限で二枝を等重みに取る場合には平均は零だが、個別の枝は逆の電流を持つ。

:::example[磁束の半周期の両側で基底を選ぶ]
$Q>0$ とする。$\alpha=1/4$ の基底は $n=0$、$\alpha=3/4$ は $n=1$ である。従って $I_n/[Q\hbar/(2\pi I)]$ は順に $-1/4,+1/4$ になる。エネルギーはいずれも $\hbar^2/(32I)$ だが、傾きは逆である。磁場が粒子の経路の上で零でも、境界条件に入る位相が運動学的な角運動量をずらす。
:::

## トンネル効果とインスタントン

障壁の下では実時間の古典軌道が両側をつなげない。虚時間で :eqref[eq-aqm-euclidean-action] を停留させると $m\ddot x=+V'(x)$ となり、反転したポテンシャルで軌道が得られる。両端が同じ高さの極小で $V=0,\dot x=0$ に近づく軌道では $m\dot x^2/2-V=0$ が保存される。従って単調経路の作用は $\int\sqrt{2mV(x)}\dd x$ になる。この局在した虚時間の遷移軌道を :term[インスタントン]{en="instanton" yomi="いんすたんとん"} と呼ぶ。

対称二重井戸で $V=m\omega^2(x^2-a^2)^2/(8a^2)$ とする。極小の曲率は $m\omega^2$ で、$-a<x<a$ の正方向の軌道は $\dot x=\omega(a^2-x^2)/(2a)$ を満たす。変数分離から $x=a\tanh(\omega(t_E-t_0)/2)$、作用は $S_0=(m\omega/(2a))\int_{-a}^{a}(a^2-x^2)\dd x=2m\omega a^2/3$ である。

### 希薄な遷移の加算

$S_0/\hbar\gg1$ のとき、谷の近くで過ごす時間に比べて遷移の持続時間が短い経路を加算できる。一回の遷移の単位虚時間当たりの重みを $d=Ae^{-S_0/\hbar}$ と置く。$A$ は零モードを時間並進 $t_0$ に置き換えた残りの揺らぎ積分で決まる周波数である。独立な $n$ 回の遷移は順序付き中心時刻の積分 $\tau^n/n!$ を持つ。同じ井戸に戻るには偶数回、反対側に着くには奇数回が必要で、和は $\cosh(d\tau)$ と $\sinh(d\tau)$ になる。共通因子 $e^{-E_0\tau/\hbar}$ を含めて対称・反対称結合に移ると $E_0\mp\hbar d$ であり、分裂は $2\hbar A e^{-S_0/\hbar}$ に比例する。四次ポテンシャルでの前因子は、この節の後半で揺らぎを積分して求める。

:::example[矩形障壁での指数]
幅 $L$、高さ $V_0>E$ の領域では定常方程式が $u''=\kappa^2u$、$\kappa=\sqrt{2m(V_0-E)}/\hbar$ となる。厚い障壁の透過振幅の指数は $e^{-\kappa L}$、透過率は $e^{-2\kappa L}$ に比例する。$V$ からエネルギーを引いた量を障壁上で用いる点と、振幅を二乗する点が重要である。
:::

:::example[質量を変える]
ポテンシャルの形を固定して質量を $4m$ にすると $\int\sqrt{2mV}\dd x$ は二倍になる。したがって分裂の指数因子は元の指数因子の二乗となる。井戸の周波数も変わるため、前因子は別に評価する。
:::

### 二重井戸の揺らぎの前因子

上の四次ポテンシャルで作用の二次変分を質量で割った演算子は

$$
D=-\frac{\mathrm d^2}{\mathrm dt_E^2}
 +\omega^2-\frac{3\omega^2}{2}\operatorname{sech}^2(\omega(t_E-t_0)/2)
$$

である。$V''(x)=m\omega^2(3x^2-a^2)/(2a^2)$ を代入した。遷移中心を $t_0=0$ とし、$y_0(t)=\operatorname{sech}^2(\omega t/2)$ とすると $Dy_0=0$ である。これは運動方程式を時間微分したものに等しい。$w=-y_0'/y_0$ と置けば $D=(-\partial_t+w)(\partial_t+w)$ である。Dirichlet 境界で部分積分すると二次形式は $\|(\partial_t+w)y\|^2$ なので負の固有値がない。独立解は $y_1(t)=y_0(t)I(t)$、$I(t)=\int_0^t y_0(u)^{-2}\dd u$ であり、Wronskian は 1 である。

区間 $[-L,L]$ で Dirichlet 境界を置く。$Dy=0$ を $y(-L)=0,y'(-L)=1$ で解くと、右端では $y(L)=2y_0(L)^2 I(L)$ である。基準演算子 $D_{\rm vac}=-\mathrm d^2/\mathrm dt^2+\omega^2$ では右端の解は $\sinh(2\omega L)/\omega$ になる。両行列式の比はこの右端値の比である。実際、格子間隔 $\varepsilon$ の無次元三重対角行列では $n$ 次行列式が $d_n=(2+\varepsilon^2 U_n)d_{n-1}-d_{n-2}$ に従う。$d_0=1$、$d_{-1}=0$ なので同じ差分方程式の境界解 $y_{n+1}=\varepsilon d_n$ である。比を取って格子極限に移れば対応が得られる。

$y_0(L)\sim4e^{-\omega L}$、$I(L)\sim e^{2\omega L}/(32\omega)$ なので $\det D/\det D_{\rm vac}\sim 2e^{-2\omega L}$ である。この中に時間並進に対応する最小固有値が含まれる。これを計算するには、偶解 $y(t;\lambda)$ を $y(0;\lambda)=1,y'(0;\lambda)=0$ とし、$(D-\lambda)y=0$ とする。ゼロでの微分は定数変化法により

$$
\partial_\lambda y(L;0)
=-y_1(L)\int_0^L y_0(t)^2\dd t
+y_0(L)\int_0^L y_0(t)y_1(t)\dd t
$$

である。第 1 項が先頭で、第 2 項は $L e^{-\omega L}$ の次数である。$\int_0^\infty y_0^2\dd t=4/(3\omega)$ より境界値の零は $\lambda_0\sim 24\omega^2 e^{-2\omega L}$ となる。高次項が小さいことは、定数変化の積分核を反復して各次数に $L$ の冪と $e^{\omega L}$ の上界を得ると分かる。具体的に $0\le u\le t$ で核 $y_1(t)y_0(u)-y_0(t)y_1(u)$ は定数倍の $e^{\omega(t-u)}$ 以下であるため、$n$ 回反復は定数の $n$ 乗と $L^n/n!$ を掛けた $e^{\omega L}$ 以下である。$\lambda=\order{e^{-2\omega L}}$ では二次以上は右端の $y_0(L)$ より小さい。従って零モードを除いた比は $\det{}'D/\det D_{\rm vac}=1/(12\omega^2)$ である。

正規化した並進方向に沿う変位係数 $b$ について、$m\int \dot x_{\rm cl}^2\dd t=S_0$ なので $\sqrt{m}\,|\dd b|=\sqrt{S_0}\,|\dd t_0|$ である。Gaussian の規格化 $\sqrt{m/(2\pi\hbar)}\,\dd b$ が $\sqrt{S_0/(2\pi\hbar)}\,\dd t_0$ に変わる。残りの揺らぎからは行列式比の逆平方根 $\sqrt{12}\omega$ が出る。従ってこの模型の希薄近似の一ループ前因子は

$$
A=\omega\sqrt{\frac{6S_0}{\pi\hbar}}.
$$

この評価は井戸間距離が調和振動の幅より十分広い極限で用いる。例えば $S_0/\hbar=10$ なら $\Delta E/(\hbar\omega)\simeq2\sqrt{60/\pi}\,e^{-10}\simeq4.0\times10^{-4}$ である。準安定井戸の崩壊率を求める場合には、鞍点まわりの負モードの積分路を変形してエネルギーの虚部を定める必要がある。その解析接続の扱いは :cite[coleman-aspects] を参照せよ。

### 偏りのある二井戸の有効模型

左右の局在基底で、トンネル結合の大きさを $J>0$、井戸のエネルギー差の半分を $\varepsilon$ とする。高い準位とのギャップが両者より十分大きければ、$H=E_0I+\varepsilon\sigma_z-J\sigma_x$ を考える。Pauli 行列の反交換で $(H-E_0I)^2=(\varepsilon^2+J^2)I$ なので、固有値差は $2\sqrt{\varepsilon^2+J^2}$ である。時間の指数を偶数項と奇数項に分けると、左から右へ移る確率は
$$
p_{L\to R}(t)=\frac{J^2}{\varepsilon^2+J^2}\sin^2\left(\frac{\sqrt{\varepsilon^2+J^2}\,t}{\hbar}\right).
$$
対称なときにインスタントンで得た $J=\hbar d$ を使える。偏りが障壁の形まで変える場合には結合を再評価する。

:::example[分裂が大きくても移りにくい場合]
$\varepsilon=3J/4$ なら平方根は $5J/4$ で、最大移行確率は $16/25$ である。最初の最大時刻は $2\pi\hbar/(5J)$ になる。対称二井戸では最大確率が 1 で、最初の最大は $\pi\hbar/(2J)$ である。エネルギー分裂が速い振動を与えても、局在基底から両固有状態への射影が等しいとは限らない。
:::

### 障壁の厳密な接続と指数の取り出し

両側で $V=0$、$0<x<L$ で $V_0>E$ の一次元障壁を考える。内部の初期値 $(\psi,\psi')$ を $L$ へ運ぶ行列は
$$
\begin{pmatrix}\cosh\kappa L&\sinh\kappa L/\kappa\\
\kappa\sinh\kappa L&\cosh\kappa L\end{pmatrix}.
$$
各列は $u''=\kappa^2u$ の二つの初期値 $(1,0),(0,1)$ の解である。外部波数を $k=\sqrt{2mE}/\hbar$ とする。左の値と微分は $1+r_A,ik(1-r_A)$、右のそれは $t_Ae^{ikL},ikt_Ae^{ikL}$ である。右側の列に逆行列を掛けて左へ戻すと、左の入射係数 1 は
$$
t_Ae^{ikL}\left[\cosh\kappa L+\frac{i}{2}\left(\frac{\kappa}{k}-\frac k\kappa\right)\sinh\kappa L\right]
$$
となる。従って同じ外部速度での透過率は
$$
\mathcal T=\left[1+\frac{V_0^2\sinh^2(\kappa L)}{4E(V_0-E)}\right]^{-1}.
$$
恒等式 $\cosh^2-\sinh^2=1$ と $k^2+\kappa^2=2mV_0/\hbar^2$ を使った。$\kappa L\gg1$ では $\mathcal T\simeq 16E(V_0-E)e^{-2\kappa L}/V_0^2$ であり、指数だけを使うと境界で決まる係数が失われる。

:::example[障壁の中央エネルギーで幅を変える]
$E=V_0/2$ なら $k=\kappa$ で $\mathcal T=\operatorname{sech}^2(\kappa L)$ である。$\kappa L=\log2$ では $\cosh(\log2)=5/4$ なので率は $16/25$ である。$\kappa L=\log4$ では $64/289$ になる。厚い極限の $4 e^{-2\kappa L}$ は前者で 1、後者で $1/4$ となり、小さい幅での近似誤差が大きい。接続係数まで保持した式を使えば、指数だけを比較してよい範囲が明確になる。
:::

### 零モードの積分を時間で規格化する

インスタントンの移動量を $t_0$ として $x_{\rm cl}(t-t_0)$ を微分すると $-\dot x_{\rm cl}$ が出る。$\int m\dot x_{\rm cl}^2\dd t=S_0$ なので、質量を吸収した揺らぎの規格化ベクトルは $\sqrt m\,\dot x_{\rm cl}/\sqrt{S_0}$ である。時間のずれ $\dd t_0$ による変位をこの方向に射影すると長さの代わりに $\sqrt{S_0}\dd t_0$ が測度に出る。有限区間を先に取る場合はモードの最小固有値を除いた比を用い、中心 $t_0$ の積分と二重に数えない。独立に動かせる中心が $n$ 個あるとき、順序付き領域の体積は区間 $\tau$ の $n$ 次元の立方体を $n!$ 個に分けた一つなので $\tau^n/n!$ である。端に近い遷移や互いに近い遷移は、この体積の全てで独立という近似を破る。密度 $d$ が $\omega$ より小さいとき、典型的な間隔 $1/d$ は中心の幅 $1/\omega$ より大きい。

:::example[虚時間核の非対角から分裂を推定する]
同じ井戸への核と反対側への核の比が、ある長さ $\tau$ で $1/3$ であったとする。希薄近似ではこれが $\tanh(d\tau)$ なので、$d\tau=\operatorname{artanh}(1/3)=\log2/2$ である。従って分裂 $\Delta E=2\hbar d=\hbar\log2/\tau$ となる。時間を $2\tau$ にすれば $\tanh(2d\tau)=2(1/3)/(1+(1/3)^2)=3/5$ と予測する。両時間の結果を比べれば、共通の $d$ の近似がどこまで合うかを判断できる。
:::

### 障壁の指数を WKB と虚時間で照合する

一次元の定常式 $-\hbar^2\psi''/(2m)+V\psi=E\psi$ を、滑らかで $V>E$ の区間で考える。$\psi=A(x)e^{-W(x)/\hbar}$ を代入すると、次数 $\hbar^0$ で $(W')^2=2m(V-E)$ である。減衰枝に $W'>0$ を選ぶ。次の次数は $2W'A'+W''A=0$ なので $A\propto(W')^{-1/2}$ になる。転回点で $W'=0$ のためこの式は破れ、Airy 型の局所接続または別の境界計算が必要である。指数に関しては虚時間の古典保存量 $m\dot x^2/2-V=-E$ を使い、$S_E-E\tau=\int (m\dot x^2/2+V-E)\dd t=\int\sqrt{2m(V-E)}\dd x$ を得る。エネルギーを指定した区間の指数が、二つの方法で同じになる。井戸間をつなぐ基底のインスタントンでは極小のポテンシャルを零に取り、この式の $E=0$ にあたる指数を先頭に使う。零点運動は前因子や高次の評価に入る。

:::example[滑らかな放物障壁の被積分関数]
最大点の近くで $V(x)=V_b-m\Omega^2 x^2/2$、$E<V_b$ とする。転回点は $\pm a$、$a^2=2(V_b-E)/(m\Omega^2)$ である。障壁指数は $\int_{-a}^a m\Omega\sqrt{a^2-x^2}\dd x=m\Omega\pi a^2/2=\pi(V_b-E)/\Omega$ である。深い下の透過率の指数は $e^{-2\pi(V_b-E)/(\hbar\Omega)}$ となる。例えば $(V_b-E)/(\hbar\Omega)=2$ なら指数は $e^{-4\pi}$ である。頂上へ近づくと作用が $\hbar$ の尺度になり、同じ指数の先頭だけでは高精度の率を得られない。
:::

### 滑らかな転回点で減衰波を接続する

$x_t$ で $V(x_t)=E$、$V'(x_t)=F>0$ とする。左は許容領域、右は禁制領域である。ポテンシャルを $V(x)\simeq E+F(x-x_t)$ と線形化し、$\ell_t=(\hbar^2/(2mF))^{1/3}$、$z=(x-x_t)/\ell_t$ を使うと定常方程式は $\psi_{zz}=z\psi$ となる。右で減衰する解 $\operatorname{Ai}(z)$ は、$z\to+\infty$ で $(2\sqrt\pi)^{-1}z^{-1/4}e^{-2z^{3/2}/3}$、$z\to-\infty$ で $\pi^{-1/2}|z|^{-1/4}\sin(2|z|^{3/2}/3+\pi/4)$ に漸近する。振動積分の位相 $y^3/3+zy$ を使うと、$z<0$ では二実停留点の Fresnel 位相が正弦に入る。$z>0$ では積分路を上に移し、鞍点 $i\sqrt z$ の周りで $y=i\sqrt z+s$ とすれば指数は $-2z^{3/2}/3-\sqrt z\,s^2+is^3/3$ となる。Gaussian 積分で右側の係数を得る。積分路の移動は有限の端を持つ積分から行い、垂直の端を無限遠に移すと $e^{-\operatorname{Im}(y)\operatorname{Re}(y)^2}$ の減衰で寄与が消える。

$p(x)=\sqrt{2m(E-V(x))}$、$q(x)=\sqrt{2m(V(x)-E)}$ を各領域で正に取る。線形域で $\int_{x_t}^x q\dd x/\hbar=2z^{3/2}/3$、$q^{-1/2}$ は $z^{-1/4}$ に比例する。左では同じ比例定数で $p^{-1/2}$ と $|z|^{-1/4}$ が対応する。従って右側の減衰枝 $Cq^{-1/2} e^{-\int_{x_t}^x q\dd x/\hbar}$ は、左側の $2Cp^{-1/2}\sin(\int_x^{x_t}p\dd x/\hbar+\pi/4)$ につながる。$\ell_t$ より離れ、ポテンシャルの線形化がまだ成立する重なり領域で係数を比較した。左が禁制の転回点では座標を反転して使う。

たとえば二つの単純転回点 $x_-,x_+$ の間に粒子が束縛され、両側で減衰する場合、許容領域では両側から得る正弦が比例しなければならない。左からの位相を $\theta(x)=\hbar^{-1}\int_{x_-}^{x}p\dd x$、全積分を $I=\hbar^{-1}\int_{x_-}^{x_+}p\dd x$ とすると、二つの正弦は $\sin(\theta+\pi/4)$ と $\sin(I-\theta+\pi/4)$ である。二係数を比較すると $I=(n+1/2)\pi$ が必要である。調和振動子で $p=\sqrt{2mE-m^2\omega^2x^2}$ を積分すると $I=\pi E/(\hbar\omega)$ となり、既知の準位を再現する。一般の井戸ではこの量子化式は半古典の精度で用いる。

## 演習問題

::::exercise{level=basic}
$H\mapsto H+\Delta I$ で核がどう変化するか求めよ。
:::solution
$\Delta I$ は $H$ と可換なので時間発展の指数は積に分かれる。位置表現で $K\mapsto e^{-i\Delta T/\hbar}K$ である。$\Delta T=\pi\hbar/2$ なら係数は $-i$ である。
:::
::::

::::exercise{level=basic}
自由核で $x_b-x_a=d$、伝播時間 $T$ の古典作用と連続端点の混合微分を求めよ。
:::solution
直線解は $x=x_a+dt/T$ で、運動項を積分して $S=md^2/(2T)$ である。$\partial_{x_a}S=-md/T$ をさらに $x_b$ で微分すると $-m/T$ になる。
:::
::::

::::exercise{level=basic}
同じ時間幅を二回使い、自由粒子の中間点の停留値を求めよ。
:::solution
中間点 $x$ の部分は $(x-x_a)^2+(x_b-x)^2$ に比例する。微分は $4x-2(x_a+x_b)$ なので、停留値は $(x_a+x_b)/2$ である。二階微分は正である。
:::
::::

::::exercise{level=basic}
リングで $Q\Phi/h=1/4$ とする。$E_1/E_0$ を求めよ。
:::solution
共通係数を除けば $E_1$ は $(1-1/4)^2=9/16$、$E_0$ は $(0-1/4)^2=1/16$ なので比は $9$ である。
:::
::::

::::exercise{level=advanced}
振動子の $Z$ から平均熱エネルギーを計算せよ。
:::solution
$\langle H\rangle=-\partial_\beta\log Z$ である。$Z=[2\sinh(\beta\hbar\omega/2)]^{-1}$ を微分すると $\langle H\rangle=(\hbar\omega/2)\coth(\beta\hbar\omega/2)$ となる。低温極限は $\hbar\omega/2$、高温では $\coth x\simeq1/x$ より $1/\beta$ である。
:::
::::

::::exercise{level=advanced}
二重井戸の $x=a\tanh(\omega t_E/2)$ が零の Euclid エネルギーを持つことを直接示せ。
:::solution
$\dot x=(a\omega/2)\operatorname{sech}^2(\omega t_E/2)$ である。一方 $x^2-a^2=-a^2\operatorname{sech}^2(\omega t_E/2)$ なので $m\dot x^2/2=V=m a^2\omega^2\operatorname{sech}^4(\omega t_E/2)/8$ となる。差が零なので作用の被積分関数 は $2V$ である。
:::
::::

::::exercise{level=advanced}
実対称な二変数の位相 Hessian の固有値が $2,-8$ である。定数振幅の Fresnel 前因子を求めよ。
:::solution
直交座標変換で測度を保って分離する。各固有値 $\lambda$ の因子は $\sqrt{2\pi\hbar/|\lambda|}e^{i\pi\sgn\lambda/4}$ である。二方向の位相が相殺し、大きさは $2\pi\hbar/\sqrt{16}=\pi\hbar/2$ となる。
:::
::::

::::exercise{level=advanced}
虚時間 $\tau$ が十分長いとき、対称二重井戸で同じ井戸と反対側への核の比 $K_{\mathrm{opp}}/K_{\mathrm{same}}$ を希薄近似で求め、長時間極限を説明せよ。
:::solution
偶数回和と奇数回和の共通因子が消えるので比は $\tanh(d\tau)$ である。$d>0$ なら $\tau\to\infty$ で 1 に収束する。これは両井戸に同じ符号で重みを持つ対称状態が最低エネルギーとなり、虚時間で支配的になることに対応する。
:::
::::

::::exercise{level=advanced}
四つの等間隔の時間区間で振動子を伝播させる。幅を $\varepsilon$ とし、$2-\omega^2\varepsilon^2=2\cos\vartheta$ と置く。内部 Hessian の行列式と、自由な値に対する Gaussian 前因子の比を求めよ。$0<\vartheta<\pi/4$ とする。

:::solution
境界は固定しているので揺らぎの変数は内部の三点だけである。$m/\varepsilon$ を抜くと行列は
$$
A=\begin{pmatrix}
2\cos\vartheta&-1&0\\
-1&2\cos\vartheta&-1\\
0&-1&2\cos\vartheta
\end{pmatrix}.
$$
$z=2\cos\vartheta$ とすれば、先頭行で展開して $\det A=z(z^2-1)-z=z^3-2z$ である。四倍角の公式で $z^3-2z=\sin(4\vartheta)/\sin\vartheta$ になる。自由な場合は $\vartheta\to0$ で行列式は 4 である。三つの固有値は $2\cos\vartheta-2\cos(s\pi/4)$、$s=1,2,3$ である。与えた範囲ならいずれも正なので、自由との比に追加の Fresnel 位相はない。従って前因子の比は $\sqrt{4\sin\vartheta/\sin4\vartheta}$ である。例えば $\vartheta=\pi/6$ では行列式が $\sqrt3$、比は $2/3^{1/4}$ になる。

自由に近いときに比が 1 を超えるのは、振動子のポテンシャル項が内部の二次変分を小さくするためである。$\vartheta=\pi/4$ に近づくと最小固有値が零となり、非退化な Gaussian の公式では発散する。この有限分割での焦点の条件は $4\vartheta=\pi$ であり、連続極限では $\vartheta\simeq\omega\varepsilon$ を使って $\omega T=\pi$ へ移る。
:::
::::

## 確認テスト

::::quiz
虚時間で経路に掛かる因子を選べ。

- [ ] $e^{i S_E/\hbar}$
- [x] $e^{-S_E/\hbar}$
- [ ] $e^{S_E/\hbar}$

:::explanation
時間を $-i\tau$ に回すと発展演算子は $e^{-\tau H/\hbar}$ となる。運動項の Gaussian とポテンシャルの各区間の重みを掛けると指数は $-S_E/\hbar$ になる。
:::
::::

::::quiz{answer="2" tol="1e-9"}
二経路の干渉で両経路の前因子は同じ、作用差は $\pi\hbar/2$ とする。一経路の強度に対する合計強度の比を求めよ。
:::explanation
相対振幅は $e^{i\pi/2}=i$ である。共通の振幅を除くと和は $1+i$、絶対値の二乗は 2 となる。
:::
::::

::::quiz
時間 $T$ が半周期より小さい振動子で、古典作用だけから正しい振幅を求めるために何を追加するか。

- [x] Hessian の Gaussian 積分による前因子
- [ ] 古典経路の運動量に対する確率
- [ ] 電磁場に依存しない磁束の位相

:::explanation
二次の作用では停留点の周囲の積分が厳密に Gaussian となり、$[m\omega/(2\pi i\hbar\sin\omega T)]^{1/2}$ が残る。これが核の次元と合成則に必要である。
:::
::::

::::quiz{answer="0.5" tol="1e-9"}
$Q>0$ の干渉系で、磁束零から相対位相を $\pi$ 変える最小の正の磁束を $h/Q$ の単位で答えよ。
:::explanation
相対位相は $Q\Phi/\hbar$ である。これを $\pi$ にすると $\Phi=\pi\hbar/Q=h/(2Q)$ になるので、指定単位では $1/2$ である。
:::
::::

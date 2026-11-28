---
title: 宇宙論的 Boltzmann 方程式
description: 光子と粒子の分布を角度展開し、散乱で結ばれる宇宙論的な一次方程式と初期条件を導く。
notation:
  - symbol: 'q'
    meaning: 共動運動量の大きさ
  - symbol: '\mu'
    meaning: 波数と進行方向の余弦
  - symbol: '\Theta_\ell'
    meaning: 光子温度の多重極
  - symbol: '\dot\kappa_T'
    meaning: 正の共形散乱率
  - symbol: 'R_b'
    meaning: バリオンと光子の慣性比
---

## この章の目標

- 測地線に沿った分布の変化を温度多重極へ分解する。
- Thomson 散乱で光子とバリオンの運動量が交換されることを示す。
- 自由流のニュートリノと断熱初期値を計算する。

## 光子の Boltzmann 方程式

自然単位、符号 $(-,+,+,+)$、平坦背景の Newton ゲージを採用する。$\eta$ 微分をプライムと書く。局所直交系で測った運動量 $p$ と共動運動量 $q=ap$ を用い、進行方向を $\vb n$、波数との余弦を $\mu=\hat{\vb k}\cdot\vb n$ とする。一次の測地線式は $q'/q=\Phi'-n^i\partial_i\Psi$ である。時間依存の空間尺度の摂動が $\Phi'$、時間計量の勾配による重力赤方偏移が $-n^i\partial_i\Psi$ を与える。方向の微分は背景が等方なため一次 Liouville 項に現れない。$f=f_0(q)-qf_0'(q)\Theta$（ここで $f_0'$ は $q$ 微分）と置けば、温度の摂動 $\Theta$ の式は

$$
\Theta'+ik\mu\Theta=\Phi'-ik\mu\Psi+\mathcal C_\Theta . \label{eq-cosmo-photon-boltzmann}
$$

:eqref[eq-neq-boltzmann] の輸送を曲がった時空で使った形である。角度の :term[多重極]{en="multipole" yomi="たじゅうきょく"} を $\Theta=\sum_\ell(-i)^\ell(2\ell+1)\Theta_\ell P_\ell(\mu)$ で定義し、Legendre 恒等式 $(2\ell+1)\mu P_\ell=(\ell+1)P_{\ell+1}+\ell P_{\ell-1}$ を掛けて射影する。自由流で $\Theta_\ell'=k[\ell\Theta_{\ell-1}-(\ell+1)\Theta_{\ell+1}]/(2\ell+1)$、重力の補正は $\ell=0$ で $\Phi'$、$\ell=1$ で $k\Psi/3$ である。$T^4$ の一次温度変化から $\delta_\gamma=4\Theta_0$ であり、運動量積分から $\theta_\gamma=3k\Theta_1$、跡なし応力から $\sigma_\gamma=2\Theta_2$ となる。

:::example[等方からの自由流]
初期に $\Theta=A$ で重力なしなら $\Theta=Ae^{-ik\mu\Delta\eta}$、従って $\Theta_\ell=A j_\ell(k\Delta\eta)$ である。$j_\ell$ は球 Bessel 関数である。
:::


:::example[短時間の双極子]
$j_1(x)\simeq x/3$ なので $\Theta_1\simeq A k\Delta\eta/3$ であり、四重極は二次から現れる。
:::

### 重力源を固有エネルギーの変化へ分ける

Newton ゲージの時刻面に垂直な単位観測者を $u$ とし、粒子の四元運動量を $P=Eu+p n$ と分解する。ここで $E=-P\cdot u$、$n$ は局所の空間単位ベクトルである。測地線に沿うパラメータ $\lambda$ に対し $\dd E/\dd\lambda=-P^\mu P^\nu\nabla_\mu u_\nu$ となる。$\nabla u=K-u\otimes A$ と分けると $K$ は時刻面の膨張、$A=\grad_{\rm spatial}\log[a(1+\Psi)]$ は法線観測者の加速度である。従って $\dd E/\dd\lambda=-p^2K_{nn}-Ep A_n$ である。計量から一次で $K_{nn}=(\mathcal H-\Phi'-\mathcal H\Psi)/a$、$A_n=n^i\partial_i\Psi/a$、$\dd\eta/\dd\lambda=E/[a(1+\Psi)]$ を得る。$\dd p/\dd\eta=(E/p)\dd E/\dd\eta$ を使うと
$$
q'=q\Phi'-\varepsilon_q\,n^i\partial_i\Psi,\qquad q=ap,\quad \varepsilon_q=aE.
$$
光子では $\varepsilon_q=q$ となる。同じ符号で質量粒子を扱える。方向の変化は一次の方向摂動に作用すると二次になるが、背景が異方的ならそもそもこの簡約はできない。

:::example[静的なポテンシャル]
$\Phi'=0$、粒子が $\Psi$ の高くなる向きへ移動する場合は $q'<0$ となる。真空中の光子がポテンシャルの井戸から出るときにエネルギーを失う符号である。
:::

### 分布の不変性と運動量積分

一粒子分布 $f$ は局所直交系で測った位相空間占有数であり、観測者の座標変換で同じ粒子状態を評価すれば同じ値を持つ。測地線に沿う Liouville 作用素は、質量殻で独立な変数 $(x^i,q,n^i)$ を取ると
$$
\frac{\dd f}{\dd\eta}
=\partial_\eta f+\frac{\dd x^i}{\dd\eta}\partial_i f
  +q'\partial_q f+\frac{\dd n^i}{\dd\eta}\partial_{n^i}f .
$$
背景では $\dd x^i/\dd\eta=n^i$、$q'=0$、$\dd n^i/\dd\eta=0$ である。摂動した空間速度は既に一次の $\partial_i f$ と掛かるため二次になる。同様に角度の変化は等方背景に作用すると零である。従って一次で測地線の補正が残るのは $q'\partial_q f_0$ だけである。$q'\partial_q f_0$ を左辺から移して $-q\partial_q f_0$ で割ると、:eqref[eq-cosmo-photon-boltzmann] の重力源の符号が決まる。これは散乱による周波数再分配を無視し、初期揺らぎが温度変化の形を持つ場合の分離である。

光子のエネルギー密度に対応する角度輝度は $q^3\dd q$ を重みとする。$\delta f=-q(\partial_q f_0)\Theta$ の積分は、部分積分で $\int q^3\delta f\,\dd q=4\Theta\int q^3 f_0\,\dd q$ である。次に角度モーメントを取り、$\langle\mu\Theta\rangle=-i\Theta_1$ と $\langle(\mu^2-1/3)\Theta\rangle=-2\Theta_2/3$ を直交性から得る。物質節の規約 $\theta=ik^iv_i$ と $\sigma_f=-Q^i{}_j\Pi^j{}_i/(\rho+p)$ では
$$
\delta_\gamma=4\Theta_0,\qquad
(\bar\rho_\gamma+\bar p_\gamma)\theta_\gamma=4\bar\rho_\gamma k\Theta_1,\qquad
(\bar\rho_\gamma+\bar p_\gamma)\sigma_\gamma=\frac83\bar\rho_\gamma\Theta_2 .
$$
ここで $Q^i{}_j=\hat k^i\hat k_j-\delta^i_j/3$ である。光子では $\bar\rho+\bar p=4\bar\rho/3$ なので、先ほどの多重極と流体量の対応になる。係数 $3$ と $2$ は輝度の因子 $4$ およびエンタルピーの規格化から出る。

### 視線に沿う無衝突解

固定方向で $X=\Theta+\Psi$ とおく。:eqref[eq-cosmo-photon-boltzmann] に沿った積分因子は $e^{ik\mu\eta}$ である。衝突なしでは
$$
X(\eta,\mu)=e^{-ik\mu(\eta-\eta_i)}X(\eta_i,\mu)
+\int_{\eta_i}^\eta e^{-ik\mu(\eta-u)}[\Phi'(u)+\Psi'(u)]\dd u.
$$
ポテンシャルの時間微分が零なら、初期の有効温度を運ぶだけである。時間依存の重力は連続して新しい球 Bessel モードを作るため、初期に単極子だけあっても、後の各多重極は単一の $j_\ell$ で表せない。

:::example[初期の純粋な双極子が作る単極子]
重力を零とし、初期分布が $\Theta=-3iD\mu$ とする。角度平均では $\partial_x e^{-ix\mu}=-i\mu e^{-ix\mu}$ より $\Theta_0=3D j_0'(x)=-3D j_1(x)$、$x=k(\eta-\eta_i)$ である。短時間では $\Theta_0=-Dk(\eta-\eta_i)+\order{x^3}$ であり、階層の $\Theta_0'=-k\Theta_1$ と一致する。
:::

### 光線の進行方向とテンソルの赤方偏移

テンソルだけの計量 $g_{ij}=a^2(\delta_{ij}+h_{ij})$ を使う。法線観測者は $u^\mu=a^{-1}(1,0,0,0)$ で加速度は零である。混合成分の外部曲率は $K^i{}_j=a^{-1}(\mathcal H\delta^i_j+h^{i\prime}{}_j/2)$ なので、法線方向のエネルギーの変化を使えば
$$
\frac{q'}q=-\frac12 n^in^j h_{ij}' .
$$
ここでも方向の一次変化は等方な背景分布へ作用しない。従ってテンソルによる無衝突の温度輸送は $\Theta'+ik\mu\Theta=-n^in^j h_{ij}'/2$ である。$\vb k$ を第三軸に取り、プラス偏光を $e_{xx}=-e_{yy}=1/\sqrt2$ と選べば、源の角度依存は $(1-\mu^2)\cos(2\varphi)/\sqrt2$ である。スカラーの場合のような $\mu$ だけの Legendre 階層ではこの角度依存を表せない。クロス偏光では $\cos(2\varphi)$ の代わりに $\sin(2\varphi)$ が現れ、二つの方位角の直交性により線形では互いに独立である。

:::example[一定のテンソル計量では生じない赤方偏移]
$h_{ij}'=0$ で初期の分布が等方なら、自由流で一次のテンソル温度は生成されない。逆に短い区間で $k\Delta\eta\ll1$ なら輸送の位相変化を無視でき、テンソルの変化 $\Delta h_{ij}$ に対し $\Delta\Theta\simeq-n^in^j\Delta h_{ij}/2$ である。$x$ 軸と $y$ 軸に進む光で変化の符号が逆になる。
:::

### 初期の角度分布を直接射影する

無衝突で静的な両重力ポテンシャルの下、初期の有効温度を $X_i(\mu)=A-3iD\mu-5Q P_2(\mu)$ とする。自由流により $X=X_i e^{-ix\mu}$、$x=k(\eta-\eta_i)$ である。微分演算 $\partial_x$ が $-i\mu$ を掛けることと、$P_2=(3\mu^2-1)/2$ を使うと、全ての多重極は
$$
X_\ell(x)= A j_\ell(x)+3D j_\ell'(x)
+\frac52 Q[3j_\ell''(x)+j_\ell(x)]
$$
である。特に $x\to0$ では $j_0(0)=1,j_1'(0)=1/3,j_2''(0)=2/15$ なので $X_0=A,X_1=D,X_2=Q$ に戻る。単極子に入る四重極の寄与は、$j_0(x)=1-x^2/6+x^4/120+\cdots$ から $Qx^2/3+\cdots$ である。これは階層の初期で $X_1'=-(2k/3)Q$、$X_0'=-kX_1$ を使っても得られる。実際には時間依存の重力源を各時刻の $x$ で積分するため、一つの $x$ の微分だけではまとめられない。

この式で用いる $X_\ell$ は $X=\Theta+\Psi$ の多重極である。ポテンシャルの $\Psi$ 自体に角度依存がないので $\Theta_\ell=X_\ell$ は $\ell\ge1$ で成り立つが、単極子では $\Theta_0=X_0-\Psi$ である。静的な重力の下で $X_i=0$ と置けば $X$ は後も零で、初期の $\Theta_0=-\Psi$ を維持する。この重力平衡と、初めの $\Theta_0=0$ から自由流で重力の変化を測る場合を区別する。

## Thomson 散乱による衝突項

正の散乱率を $\dot\kappa_T=a n_e\sigma_T$ と書く。物質の静止系で非偏光の散乱核は :eqref[eq-cosmo-compton] に続く Thomson 極限の $3(1+\cos^2\vartheta)/(16\pi)$ である。これを $(4\pi)^{-1}[1+P_2(\cos\vartheta)/2]$ と書けば、回転対称の畳み込みで単極子と四重極だけが入射に寄与する。バリオン速度を $\vb v_b$ として、偏光を捨てた一時的な角度式は
$$
\mathcal C_\Theta=\dot\kappa_T[-\Theta+\Theta_0+\vb n\cdot\vb v_b-\tfrac12 P_2(\mu)\Theta_2].
$$
速度項は電子の静止系と宇宙の座標系の一次 Doppler 変換である。厳密には偏光も生成される。光の二偏光の強度行列を横な三次元行列 $J_{ij}(\vb n)$ とし、$P_{ij}(\vb n)=\delta_{ij}-n_i n_j$ とする。散乱の利得は $\dot\kappa_T(3/8\pi)P(\vb n)[\int J(\vb n')\dd\Omega']P(\vb n)$ である。双極子加速度の場を出射偏光へ射影した振幅を二乗するとこのテンソル積が現れ、等方な非偏光を代入すると損失が相殺する。跡の差の方向依存が線偏光になる。

:::example[等方の温度]
$\Theta=\Theta_0$、$\vb v_b=0$ なら衝突項は零である。
:::


:::example[四重極の減衰]
偏光を捨て、輸送を止めると $\ell=2$ の利得は $\dot\kappa_T\Theta_2/10$ で、正味の減衰率は $0.9\dot\kappa_T$ となる。
:::

### 強度行列から温度と偏光の利得を得る

軸対称なスカラー摂動を波数軸のまわりに展開する。方向の極角と方位角の基底を $(\vb e_\theta,\vb e_\varphi)$ とし、温度と同じ規格化を採る偏光振幅を $P$ とする。横強度行列の対角は $(\Theta+P)/2,(\Theta-P)/2$ である。入射の行列を全方向で積分したものを $A$ とし、$D=A_{zz}-A_{xx}$ と置く。対称性で $A_{xx}=A_{yy}$ である。$A_{zz}=\frac12\int(\Theta+P)(1-\mu^2)\dd\Omega$ と $\tr A=\int\Theta\dd\Omega$ より、Legendre 多項式を $\mathsf P_\ell$、偏光の同じ展開規約での係数を $P_\ell^{\rm pol}$ と書くと
$$
D=\tfrac12\int[-\mathsf P_2(\mu)\Theta+(1-\mathsf P_2(\mu))P]\dd\Omega
=2\pi\Pi,\qquad
\Pi=\Theta_2+P_0^{\rm pol}+P_2^{\rm pol}.
$$
出射方向の横基底へ射影した $A$ の対角は $A_{xx}+D(1-\mu^2),A_{xx}$ である。これに $3/(8\pi)$ を掛けた和と差により、温度と偏光の散乱利得は
$$
\Theta_{\rm gain}=\Theta_0-\tfrac12\mathsf P_2(\mu)\Pi,\qquad
P_{\rm gain}=\tfrac12[1-\mathsf P_2(\mu)]\Pi .
\label{eq-cosmo-thomson-polarization}
$$
各方向から元の $\Theta,P$ を差し引き、温度の式に電子の Doppler 項を加えて散乱項となる。角度平均で温度の利得と損失は打ち消し合う。偏光の単極子は角度に沿って回転する基底での係数であり、全天で同じ向きの偏光を意味しない。

### 衝突項が保つ角度モーメント

電子が静止し、偏光なしの角度分布を考える。核 $[1+P_2(\hat{\vb n}\cdot\hat{\vb n}')/2]/(4\pi)$ は Legendre の加法定理から、球面調和関数 $Y_{\ell m}$ には固有値 $1$（$\ell=0$）、$1/10$（$\ell=2$）、零（その他）で作用する。これは $P_L(\hat{\vb n}\cdot\hat{\vb n}')=4\pi\sum_M Y_{LM}(\hat{\vb n})Y_{LM}^*(\hat{\vb n}')/(2L+1)$ を代入し、直交積分すれば得られる。利得から損失を差し引くと数の角度平均は保存され、高次の異方性は減衰する。動く電子の双極子は静止系へ戻す Doppler 補正を含む。

:::example[二つの異方性の残り方]
$\dot\kappa_T$ 一定、輸送なし、偏光を捨てた模型で、期間の光学的厚さを $10$ とする。四重極の振幅比は $e^{-9}$、八重極は $e^{-10}$ である。どちらも小さいが、四重極は散乱入によって大きい方が残る。
:::

### 温度と偏光の連立階層

$\theta_b=ik^iv_{bi}$ より $\vb n\cdot\vb v_b=-i\mu\theta_b/k$ である。角度式の $-i\mu$ を $\ell=1$ の基底 $-3iP_1$ と比較すると、双極子の利得は $\dot\kappa_T\theta_b/(3k)$ である。また $-P_2/2$ と $\ell=2$ の基底 $-5P_2$ を比較すると、四重極利得は $\dot\kappa_T\Pi/10$ である。温度の最初の式は
$$
\begin{aligned}
\Theta_0'&=-k\Theta_1+\Phi',\\
\Theta_1'&=\frac{k}{3}(\Theta_0-2\Theta_2+\Psi)
  -\dot\kappa_T(\Theta_1-\theta_b/(3k)),\\
\Theta_2'&=\frac{k}{5}(2\Theta_1-3\Theta_3)
  -\dot\kappa_T(\Theta_2-\Pi/10).
\end{aligned}
$$
$\ell\ge3$ では自由流に $-\dot\kappa_T\Theta_\ell$ を加える。偏光では全 $\ell$ が同じ自由流と損失を持ち、利得は $\ell=0$ の $\dot\kappa_T\Pi/2$ と $\ell=2$ の $\dot\kappa_T\Pi/10$ だけである。背景偏光が零なのでスカラーの重力赤方偏移はこの偏光方程式に一次の源を作らない。偏光基底は各方向で波数軸に合わせており、無衝突の平坦背景の光路で一定である。この条件で角度ごとの輸送には温度と同じ $ik\mu$ が作用する。

衝突だけを残すと $(\Theta_2,P_0^{\rm pol},P_2^{\rm pol})$ のベクトルを $X$ として $X'=-\dot\kappa_T[X-u(1,1,1)X]$、$u=(1/10,1/2,1/10)^{\mathsf T}$ である。成分を足すと $\Pi'=-(3/10)\dot\kappa_T\Pi$ となる。和が零の二つの方向は減衰率 $\dot\kappa_T$ を持つ。従って偏光を許した角度消去では、温度四重極だけの減衰率を使うと誤差になる。

:::example[散乱だけで偏光を生成する初期値]
$\Theta_2(0)=A$、$P_0^{\rm pol}(0)=P_2^{\rm pol}(0)=0$ とし、$\dot\kappa_T$ 一定で $u_T=\dot\kappa_T\eta$ とおく。$\Pi=Ae^{-3u_T/10}$ を $P_0^{\rm pol}$ の式に入れて積分すると $P_0^{\rm pol}=(5A/7)(e^{-3u_T/10}-e^{-u_T})$ になる。初期の生成率は $\dot\kappa_T A/2$ であり、輸送がないので十分時間がたつと偏光も消える。
:::

### Thomson 近似が捨てるエネルギー交換

電子の静止系で $q$ が変わらない散乱の利得と損失を考える。偏光なしのスカラー分布なら、入出射方向の関係を同じ核 $K(\vb n,\vb n')=K(\vb n',\vb n)$ で表せる。Bose の増強を含む差は $f(\vb n')[1+f(\vb n)]-f(\vb n)[1+f(\vb n')]=f(\vb n')-f(\vb n)$ である。従ってこの極限で背景占有数が大きくても線形の散乱核を使える。有限の反跳や電子温度で周波数が再配分されると、分布の異なるエネルギーが結び付き、同じ角度核だけでは光子数とエネルギーの変化を同時に計算できない。

動く電子に対して、電子系の温度は一次で $\Theta_{\rm e}(\vb n)=\Theta(\vb n)-\vb n\cdot\vb v_b$ である。$O(v_b)$ の角度のずれは等方な背景へは効果がなく、既に一次の異方性との積は二次になる。電子系で $-\Theta_{\rm e}$ を使った損失から $+\vb n\cdot\vb v_b$ が現れ、利得の単極子と四重極にはその双極子が入らない。これが宇宙座標の温度に現れる速度項の導出である。自由電子の密度の一次揺らぎも、背景の等方分布に作用する利得と損失が相殺するため、一次の温度散乱率には背景の $n_e$ を使う。再結合の揺らぎが散乱に与える効果を求めるには二次の項が必要である。

:::example[電子の静止枠で等方な光子]
宇宙座標の分布が $\Theta=\Theta_0+\vb n\cdot\vb v_b$ で、偏光が零なら電子系では等方である。温度式に代入すると $-\Theta+\Theta_0+\vb n\cdot\vb v_b=0$ であり、光子と電子の運動量は交換されない。この例と宇宙座標で $\Theta=\Theta_0$ の例を比べれば、散乱が駆動するのは座標上の光子速度そのものではなく、電子との相対速度であることが分かる。
:::

### 偏光行列が保つ全強度

行列核の規格化を方向積分で確かめる。固定した入射の横行列を $J$ とし、出射の量 $P(\vb n)JP(\vb n)$ を考える。跡は $\operatorname{tr}(JP)=\operatorname{tr}J-n^iJ_{ij}n^j$ である。等方な方向平均 $\int n^i n^j\dd\Omega=(4\pi/3)\delta^{ij}$ を用いると
$$
\frac{3}{8\pi}\int\operatorname{tr}(PJP)\dd\Omega
=\frac{3}{8\pi}\frac{8\pi}{3}\operatorname{tr}J
=\operatorname{tr}J .
$$
この関係は非偏光の核の積分と異なり、入射偏光の状態を限定しない。光子数を保つ散乱なので、損失した光は別の方向に同じ数だけ現れる。特定方向の温度は再分配されるので、任意の多重極が保たれるわけではない。

例えば非偏光の一本の光線が波数軸方向から来る場合、二つの入射電場はその軸に垂直である。出射角 $\vartheta$ の横基底のうち、散乱面に垂直な電場は射影で大きさを保ち、面内の電場は $\cos\vartheta$ が掛かる。従って出射強度の和は $1+\cos^2\vartheta$、差の絶対値は $1-\cos^2\vartheta$ に比例する。直角散乱に限れば一本の線偏光が残る。実際の原始放射では全方向の入射を加えるため、方向に依らない等方入射では偏光が零になる。

### 光学的厚さと散乱後の分布

方向ごとの一次輸送を $X'+(ik\mu+\dot\kappa_T)X=J$ と書く。$X=\Theta+\Psi$ で、$J$ に重力の時間変化と散乱の利得を集める。積分因子により
$$
X(\eta,\mu)=e^{-ik\mu(\eta-\eta_i)}e^{-\tau(\eta_i,\eta)}X_i(\mu)
+\int_{\eta_i}^\eta e^{-ik\mu(\eta-u)}e^{-\tau(u,\eta)}J(u,\mu)\dd u,
\qquad
\tau(u,\eta)=\int_u^\eta\dot\kappa_T(v)\dd v .
$$
$\tau$ は無次元なので、衝突の時間単位と $a$ を混同すると指数に次元の誤りが出る。電子の静止系で測った時間に対しては $n_e\sigma_T$ が率であり、背景の共形時刻に対しては $a n_e\sigma_T$ が率である。背景に温度や方向での散乱率の違いがなければ、無散乱確率はこの単一の $e^{-\tau}$ である。衝突の利得が未知のモーメントに依存するので、上式だけから光路の全分布は決まらず、階層を解いて $J$ を得る。

:::example[深さの符号を確かめる]
定数率で $\tau(\eta_i,\eta)=3$ なら、初期の非等方分布が無散乱で残る割合は約 $0.0498$ である。ただし光子の数全体がこの倍率に減るわけではない。先の行列核が保つ全強度は、無散乱で残った光と他方向へ散乱された光を両方数えた量である。
:::

## バリオンと暗黒物質の方程式

冷たい暗黒物質は :eqref[eq-cosmo-fluid-perturbations] の $w=\sigma_f=0$ の式に従う。非相対論的バリオンでも $\delta_b'=-\theta_b+3\Phi'$ であるが、運動方程式には音速 $c_b$ と散乱の反作用が入る。光子の双極子に対する衝突項は速度で $\dot\kappa_T(\theta_b-\theta_\gamma)$ となる。運動量密度の重みは光子が $4\bar\rho_\gamma/3$、バリオンが $\bar\rho_b$ なので、$R_b=3\bar\rho_b/(4\bar\rho_\gamma)$ として

$$
\theta_b'=-\mathcal H\theta_b+c_b^2 k^2\delta_b+k^2\Psi+\frac{\dot\kappa_T}{R_b}(\theta_\gamma-\theta_b) \label{eq-cosmo-baryon-euler}
$$

を得る。ここでは圧力摂動を $c_b^2\bar\rho_b\delta_b$ で閉じた。組成一定の理想気体で温度も変数とする場合は $\delta p_b/\bar\rho_b=[k_BT_b/(\mu m_p)](\delta_b+\delta T_b/T_b)$ を使う。$T_b$ はバリオンの温度、$\mu m_p$ は平均粒子質量である。バリオンを電子とイオンの単一流体にするには Coulomb 相互作用が十分速いことを仮定する。$\dot\kappa_T\gg k,\mathcal H$ では速度差は小さいが、衝突項に大きい率が掛かるため片方の速度式から直接消せない。慣性で重み付けして二式を足す。

:::example[速度差の減衰]
散乱項だけを残し $R_b=1/2$ を固定する局所模型 なら $\theta_\gamma-\theta_b$ の減衰率は $3\dot\kappa_T$ である。
:::


:::example[運動量の平均]
散乱項だけを残し $R_b=1$ を固定する模型、初期速度発散が一方 $2$、他方 $0$ なら共通の最終値は $1$ である。
:::

### 速度スリップの先頭次数

$S=\theta_\gamma-\theta_b$ とする。二つの速度式の差から
$$
S'+\dot\kappa_T(1+R_b^{-1})S
=k^2(\delta_\gamma/4-\sigma_\gamma-c_b^2\delta_b)+\mathcal H\theta_b
$$
となる。右辺を $F$、左辺の減衰率を $\gamma_s$ とすると、急な過渡変化が消えた後には $S\simeq F/\gamma_s$ である。この近似の誤差を調べるには $S=F/\gamma_s+\delta S$ を代入し、$\delta S'+\gamma_s\delta S=-(F/\gamma_s)'$ を解けばよい。各係数が $1/\gamma_s$ の時間よりゆっくり変わると誤差が小さい。

:::example[共通の重力が差から消える]
$R_b=1$、$\mathcal H=0$、圧力差の駆動を一定の $F$ とする。$S(0)=0$ からの解は $S=F(1-e^{-2\dot\kappa_T\eta})/(2\dot\kappa_T)$ である。$2\dot\kappa_T\eta=\log 2$ では定常の半分まで達する。同じ重力加速度は両方に掛かるため差には現れない。
:::

### バリオン温度と速度の二つの緩和

運動量の Thomson 抵抗と電子温度の Compton 緩和は別々に計算する。非相対論的で全イオンと電子が同じ温度 $T_b$ を持ち、熱容量を粒子一個当たり $3/2$ とする。自由電子一個当たり、単位固有時の正味の光子へのエネルギー流は $4\sigma_T\rho_\gamma(T_b-T_\gamma)/m_e$ である。これは黒体の反跳と誘導散乱を含むエネルギー交換を線形化した式であり、同温度で零になる。自由電子密度を $n_e$、熱を受け取る全粒子数密度を $n_{\rm tot}$ とすれば
$$
T_b'+2\mathcal H T_b=\Gamma_C(T_\gamma-T_b),\qquad
\Gamma_C=\frac{8a\sigma_T\rho_\gamma}{3m_e}\frac{n_e}{n_{\rm tot}} .
$$
膨張項は理想単原子気体が断熱なら $T_b\propto a^{-2}$ となることから出る。加熱の左辺は $(3/2)n_{\rm tot}$ を掛ければ固有エネルギーの式へ戻る。$n_{\rm tot}$ は中性粒子も含み、組成変化による潜熱や天体からの加熱を無視した。再結合後、$\Gamma_C$ がなお大きければ光子とバリオンの運動量が離れ始めていても、温度はしばらく $T_\gamma$ に近い。

:::example[緩和が強い極限]
$T_\gamma'=-\mathcal H T_\gamma$ とし、$d=T_\gamma-T_b$ の式を作ると $d'+(2\mathcal H+\Gamma_C)d=\mathcal H T_\gamma$ となる。過渡項が消え $\Gamma_C\gg\mathcal H$ なら $d/T_\gamma\simeq\mathcal H/\Gamma_C>0$ である。温度がほぼ同じに保たれるためには、膨張で温度が先に下がる気体へ光子から熱が移る必要がある。
:::

### バリオンと冷たい暗黒物質の相対モード

散乱後でバリオン圧力を無視できる波数を取り、両成分が同じ $\Psi$ を受けるとする。$S_{bc}=\delta_b-\delta_c$、$V_{bc}=\theta_b-\theta_c$ と置く。連続式の差と Euler 式の差は $S_{bc}'=-V_{bc}$、$V_{bc}'+\mathcal HV_{bc}=0$ である。従って $V_{bc}$ は $a^{-1}$ に減衰し、$S_{bc}=S_{bc,i}-a_iV_{bc,i}\int\dd\eta/a$ となる。物質優勢 $a\propto\eta^2$ なら差の密度は定数と $\eta^{-1}$ の和である。両密度の共通成長モードが $a$ と増えるため、相対差の絶対値が一定に残っても、共通の成長に対する比は減る。小尺度でバリオン圧力が再び重要な場合は、この相対方程式に圧力源を戻す。

### 光子とバリオンの系を時間積分する

定数 $R_b$、定数 $\dot\kappa_T=\kappa$ とし、速度以外の右辺を仮に与えた $A_\gamma,A_b$ とする。$V=(\theta_\gamma+R_b\theta_b)/(1+R_b)$ と $S=\theta_\gamma-\theta_b$ へ移ると
$$
V'=\frac{A_\gamma+R_b A_b}{1+R_b},\qquad
S'+\kappa(1+1/R_b)S=A_\gamma-A_b .
$$
従って衝突行列の固有値は $0$ と $-\kappa(1+1/R_b)$ である。$V$ を時間積分すれば衝突なしの合計が進み、変化の速い $S$ は指数積分で解ける。時間積分の一段で $A_\gamma-A_b$ と係数が一定なら
$$
S(\eta+\Delta\eta)=e^{-\gamma_s\Delta\eta}S(\eta)
+\frac{A_\gamma-A_b}{\gamma_s}(1-e^{-\gamma_s\Delta\eta})
$$
が厳密である。陽 Euler 法では同じ均一方程式の増幅率が $1-\gamma_s\Delta\eta$ なので、$0<\gamma_s\Delta\eta<2$ を外れると減衰を再現できない。強い散乱の初期に全多重極を陽的に積分すると、観測したい音の周期よりはるかに短いステップが必要になる。タイトカップリングを用いる理由の一つは、この短い時間尺度の過渡解が既に消えた範囲で遅い組合せを直接進められることである。

率が時間変化する場合は $\int\gamma_s\dd\eta'$ が指数に入る。漸近展開の二次のスリップは $S\simeq F/\gamma_s-\gamma_s^{-1}(F/\gamma_s)'$ である。これをもとの式へ代入すると残差は $-[\gamma_s^{-1}(F/\gamma_s)']'$ となる。率の大きさだけを見て再結合近傍まで二次展開を続けると、$\gamma_s'/\gamma_s^2$ が増えて誤差が大きくなる。切り替えには $k/\kappa,\mathcal H/\kappa,|\kappa'|/\kappa^2$ に加え、切り替え前後の低次の量の連続性を確認する。

:::example[陽 Euler 法の不安定性]
$R_b=1$、$\kappa\Delta\eta=2$、駆動なしとすれば $\gamma_s\Delta\eta=4$ である。厳密な増幅率は $e^{-4}\simeq0.0183$ だが、陽 Euler では $1-4=-3$ となる。初期の物理系が減衰していても、時間離散化が異方性を増幅する。切り替え条件を満たしていても、運動量の和を各ステップで保存しているかを別に検査する。
:::

### 強い散乱から一次の減衰包絡を得る

重力、膨張、バリオンの熱圧を外し $R_b$ と散乱率 $\kappa=\dot\kappa_T$ を一定とする。角度階層の偏光の定常部分は $P_0^{\rm pol}=\Pi/2$、$P_2^{\rm pol}=\Pi/10$ であり、温度の四重極は $(3\kappa/4)\Theta_2=2k\Theta_1/5$ を満たす。従って $\sigma_\gamma=16\theta_\gamma/(45\kappa)$ である。また差 $S=\theta_\gamma-\theta_b$ の方程式の駆動は先頭で $k^2\Theta_0$ なので $S=R_b k^2\Theta_0/[(1+R_b)\kappa]$ となる。この二量を $\kappa^{-1}$ の一次の補正として使う。

全慣性で平均した速度 $V$ の方程式は $V'=k^2(\Theta_0-\sigma_\gamma)/(1+R_b)$ であり、単極子の式は $\Theta_0'=-(V+R_bS/(1+R_b))/3$ である。$\Theta_0'$ を再び微分し、$S'$ に $\Theta_0'$、$\sigma_\gamma$ に $\theta_\gamma=-3\Theta_0'$ を代入すると
$$
\Theta_0''+\Gamma_{\rm ac}\Theta_0'+\frac{k^2}{3(1+R_b)}\Theta_0\simeq0,\qquad
\Gamma_{\rm ac}=\frac{k^2}{3(1+R_b)\kappa}
\left(\frac{R_b^2}{1+R_b}+\frac{16}{15}\right)
$$
を得る。$\Gamma_{\rm ac}$ が音の角周波数より十分小さければ、特性根の実部は $-\Gamma_{\rm ac}/2$ である。従って包絡は $e^{-\Gamma_{\rm ac}\eta/2}$ となる。散乱を増やすと光子・バリオンの速度差と光の四重極を抑えるため、強い結合の範囲では音の散逸が逆に減る。この微分式が次の数値例の比較対象であり、宇宙膨張や重力源がある場合は係数を時間で更新する。

### 有限階層で音響の減衰を検査する

短い時刻区間に背景の $\dot\kappa_T$ と $R_b$ を固定し、重力、膨張、バリオン熱圧を外す。初期は $\Theta_0=1$、速度、偏光、他の多重極が零とする。線形方程式なので振幅を $1$ にしたのは伝達比を測るための規格化である。角度依存の温度と偏光を最大 $L$ で切り、$L+1$ を零にする。高い散乱率により高次の係数が小さい範囲では、$L$ を増やして低い係数が変わらないことを確認できる。散乱が止まる場合は、この零境界のまま長く進めると反射を生じ得る。

次の計算では時間を $k^{-1}$ の単位にし、$k=1$ を用いる。状態は $(\Theta_0,\ldots,\Theta_L,P_0^{\rm pol},\ldots,P_L^{\rm pol},\theta_b)$ の順である。散乱源を含む各時間微分を評価し、古典的な四段の Runge–Kutta 法で一周期を進める。散乱の速い時間尺度でも陽的な積分を安定にするため、$\gamma_s\Delta\eta$ を十分小さく取る。

~~~python
from math import ceil, exp, pi, sqrt

def derivative(state, cutoff, opacity, inertia):
    size = cutoff + 1
    temperature = state[:size]
    polarization = state[size:2 * size]
    baryon = state[-1]
    source = temperature[2] + polarization[0] + polarization[2]

    def stream(moments):
        padded = [0.0, *moments, 0.0]
        return [
            (ell * padded[ell] - (ell + 1) * padded[ell + 2]) / (2 * ell + 1)
            for ell in range(size)
        ]

    temp_rate = stream(temperature)
    pol_rate = stream(polarization)
    for ell in range(1, size):
        temp_rate[ell] -= opacity * temperature[ell]
    for ell in range(size):
        pol_rate[ell] -= opacity * polarization[ell]
    temp_rate[1] += opacity * baryon / 3
    temp_rate[2] += opacity * source / 10
    pol_rate[0] += opacity * source / 2
    pol_rate[2] += opacity * source / 10
    baryon_rate = opacity * (3 * temperature[1] - baryon) / inertia
    return [*temp_rate, *pol_rate, baryon_rate]

def advance(state, step, cutoff, opacity, inertia):
    def evaluate(values):
        return derivative(values, cutoff, opacity, inertia)

    def shifted(values, rates, interval):
        return [value + interval * rate for value, rate in zip(values, rates)]

    first = evaluate(state)
    second = evaluate(shifted(state, first, step / 2))
    third = evaluate(shifted(state, second, step / 2))
    fourth = evaluate(shifted(state, third, step))
    return [
        value + step * (a + 2 * b + 2 * c + d) / 6
        for value, a, b, c, d in zip(state, first, second, third, fourth)
    ]

def period_response(opacity, cutoff=6, step_limit=0.15):
    inertia = 0.6
    sound_speed = 1 / sqrt(3 * (1 + inertia))
    period = 2 * pi / sound_speed
    slip_rate = opacity * (1 + 1 / inertia)
    steps = ceil(period * slip_rate / step_limit)
    state = [1.0] + [0.0] * (2 * cutoff + 2)
    for _ in range(steps):
        state = advance(state, period / steps, cutoff, opacity, inertia)
    diffusion = (
        inertia ** 2 / (1 + inertia) + 16 / 15
    ) / (6 * opacity * (1 + inertia))
    return state[0], exp(-diffusion * period)

for opacity in (16, 32, 64):
    exact, approximation = period_response(opacity)
    finer, _ = period_response(opacity, cutoff=9, step_limit=0.075)
    print(opacity, f"{exact:.6f}", f"{approximation:.6f}", f"{abs(finer-exact):.2e}")
~~~

最後の近似値は、定数係数でスリップと偏光を $\dot\kappa_T^{-1}$ の一次まで求めたときの減衰包絡である。この比較では一周期の温度振幅の差に、過渡解と周波数の高次補正も含まれる。例えば $\kappa=16,32,64$ の一周期後の温度は順に $0.891047,0.943820,0.971484$、近似は $0.890689,0.943763,0.971475$ である。散乱率を上げると近似との差が減り、同じ散乱率で時間刻みと切断次数を変えた差はこの例では $3\times10^{-14}$ 以下である。前者は流体近似の誤差、後者は連立方程式の数値誤差を測る。

## ニュートリノの方程式

質量 $m$ のニュートリノの共形エネルギーを $\varepsilon_q=\sqrt{q^2+a^2m^2}$ とする。速度 $q/\varepsilon_q$ が輸送を、$q'=q\Phi'-\varepsilon_q n^i\partial_i\Psi$ が重力項を決める。$f=f_0(1+\mathcal F)$ と書けば、粒子が衝突の影響なく進む :term[自由流]{en="free streaming" yomi="じゆうりゅう"} の極限では、脱結合後の方程式は

$$
\mathcal F'+ik\mu(q/\varepsilon_q)\mathcal F+\frac{\dd\log f_0}{\dd\log q}\left[\Phi'-ik\mu(\varepsilon_q/q)\Psi\right]=0 . \label{eq-cosmo-nu-boltzmann}
$$

角度展開は光子と同じで、質量があれば各 $q$ の速度が違う。密度は $g a^{-4}\int \varepsilon_q f\dd^3q/(2\pi)^3$、圧力は $\varepsilon_q$ を $q^2/(3\varepsilon_q)$ に替えた積分である。異方性応力は四重極に比例し :eqref[eq-cosmo-einstein-scalar] の $\Phi-\Psi$ を生じる。

:::example[高運動量]
$q\gg am$ なら $q/\varepsilon_q\to1$ となり、質量零の光子の無衝突式と同じ輸送になる。
:::


:::example[有限の質量]
$q=am$ なら速度は $1/\sqrt2$ である。粒子の共動運動量 $q$ が保存されても $a$ が増すと速度が下がる。
:::

### 運動量ごとの多重極から重力源を作る

角度展開を $\mathcal F=\sum_\ell(-i)^\ell(2\ell+1)\mathcal F_\ell P_\ell$ とし、$A=g/(2\pi^2a^4)$ とする。等方測度と応力テンソルの一粒子重みから
$$
\begin{aligned}
\delta\rho&=A\int q^2\varepsilon_q f_0\mathcal F_0\dd q,
&\delta p&=\frac A3\int \frac{q^4}{\varepsilon_q}f_0\mathcal F_0\dd q,\\
(\rho+p)\theta&=Ak\int q^3 f_0\mathcal F_1\dd q,
&(\rho+p)\sigma&=\frac{2A}{3}\int\frac{q^4}{\varepsilon_q}f_0\mathcal F_2\dd q .
\end{aligned}
$$
角度平均では $\langle\mu\mathcal F\rangle=-i\mathcal F_1$、$\langle(\mu^2-1/3)\mathcal F\rangle=-2\mathcal F_2/3$ を使った。質量なしの場合は $q/\varepsilon_q=1$ なので輸送を運動量で積分できる。重力の源にある分布微分は $\int q^4\partial_q f_0\,\dd q=-4\int q^3f_0\,\dd q$ となり、端点項は熱分布で零である。質量がある場合は各運動量の伝播速度が違うので、同じ変換で階層を閉じられない。

:::example[一種類の運動量だけの応力比]
$q=am$ の狭い運動量帯が密度摂動を作り、$\mathcal F_0$ が帯で一定とする。$\delta p/\delta\rho=q^2/(3\varepsilon_q^2)=1/6$ である。これは狭い帯の応答比であり、流体全体の音速を指定する値ではない。
:::

### 自由流の階層を有限に切る

重力源を落とし、初期に等方な一つの運動量ビンを考える。:eqref[eq-cosmo-nu-boltzmann] の解は $\mathcal F=\mathcal F_0(\eta_i)e^{-ik\mu d_q}$、$d_q=\int_{\eta_i}^\eta q\dd\eta'/\varepsilon_q(\eta')$ である。角度展開で $\mathcal F_\ell=\mathcal F_0(\eta_i)j_\ell(kd_q)$ となる。球 Bessel 関数の漸化式 $j_{\ell+1}=(2\ell+1)j_\ell/x-j_{\ell-1}$ を使えば、$x=kd_q>0$ に対しこの解を有限の階層に接続する境界として
$$
\mathcal F_{L+1}\simeq\frac{2L+1}{kd_q}\mathcal F_L-\mathcal F_{L-1}
$$
を使える。純粋な自由流の等方初期条件なら等式であり、重力で継続的に新たな多重極が生じる場合は一般には近似である。異なる生成時刻の成分は異なる $d_q$ を持つからである。

まったく同じ切断でも、速く動くビンほど同じ物理時刻で $kd_q$ が大きくなり、高多重極まで波が及ぶ。数値計算では、切断次数を上げて低次の $\ell=0,1,2$ と応力積分が変わらないかを確かめる。$x$ が極小の間は漸化式の差し引きで桁を失うので、$j_\ell(x)=x^\ell/(2\ell+1)!!+\order{x^{\ell+2}}$ を使うか、十分大きい $L$ で上端の振幅が小さいことを調べる。

:::example[短い伝播での上端の振幅]
$kd_q=0.2$ なら $j_2\simeq0.2^2/15=2.67\times10^{-3}$、$j_3\simeq0.2^3/105=7.62\times10^{-5}$ である。最初の瞬間に八重極が小さくても、伝播時間が延びればこの小波数展開だけでは切断次数を決められない。
:::

### 巨視的な背景保存と運動量ビン

背景の共動分布 $f_0(q)$ は脱結合後に固定される。$\varepsilon_q'=\mathcal H a^2m^2/\varepsilon_q$ と $A'=-4\mathcal HA$ を使うと
$$
\bar\rho'=-4\mathcal H\bar\rho
+\mathcal H A\int q^2 f_0\,\frac{a^2m^2}{\varepsilon_q}\dd q
=-3\mathcal H(\bar\rho+\bar p).
$$
最後の変形は $a^2m^2=\varepsilon_q^2-q^2$ から従う。従って背景保存は各モードを積分してから別に課す条件ではなく、衝突なしで正しく積分した結果である。$am\ll q$ に寄与が集中すれば $\bar p\simeq\bar\rho/3$、逆に $am\gg q$ なら $\bar\rho\simeq ma^{-3}g\int q^2 f_0\dd q/(2\pi^2)$ となる。非相対論的な圧力はさらに $a^{-2}$ を受け、$a^{-5}$ に減る。

$L_q=\dd\log f_0/\dd\log q$ とすれば、運動量ごとの階層は
$$
\begin{aligned}
\mathcal F_0'&=-k(q/\varepsilon_q)\mathcal F_1-L_q\Phi',\\
\mathcal F_1'&=\frac{kq}{3\varepsilon_q}(\mathcal F_0-2\mathcal F_2)
            -\frac{k\varepsilon_q}{3q}L_q\Psi,\\
\mathcal F_\ell'&=\frac{kq}{(2\ell+1)\varepsilon_q}
 [\ell\mathcal F_{\ell-1}-(\ell+1)\mathcal F_{\ell+1}],\qquad \ell\ge2 .
\end{aligned}
$$
$\Phi'$ は角度に依らず、$\Psi$ の項は $i\mu$ の射影で符号が変わる。エネルギー積分で第一式を用いると、$\mathcal F_1$ は $-(\rho+p)\theta$ を与える。$\Phi'$ の係数は $-A\int q^3\varepsilon_q\,\partial_qf_0\dd q$ である。部分積分で $A\int f_0(3q^2\varepsilon_q+q^4/\varepsilon_q)\dd q=3(\rho+p)$ を得る。従って $\delta\rho'+3\mathcal H(\delta\rho+\delta p)=(\rho+p)(3\Phi'-\theta)$ を再現する。運動量への重みを誤ると、質量が効く時期にこの恒等式が破れる。

:::example[非相対論的な熱速度]
$am$ が全ての代表的な $q$ より大きければ、$\langle v^2\rangle=\langle q^2\rangle/(a^2m^2)$ が先頭項になる。尺度因子が三倍になると二乗平均速度は九分の一、圧力密度比 $\bar p/\bar\rho\simeq\langle v^2\rangle/3$ も九分の一である。ただし共動距離 $d_q$ にはそれ以前の相対論的な移動も積算される。
:::

### 分布の運動量保存による照合

前の多重極を使い $J=(\bar\rho+\bar p)\theta=Ak\int q^3 f_0\mathcal F_1\dd q$ を微分する。$A'=-4\mathcal HA$ と $\mathcal F_1'$ の三項を代入すると、圧力の項は $k^2\delta p$、四重極の項は $-k^2(\bar\rho+\bar p)\sigma$ である。重力の項は
$$
-\frac{A k^2\Psi}3\int q^3\varepsilon_q\,\partial_q f_0\dd q
=k^2(\bar\rho+\bar p)\Psi
$$
となる。端点が消えるときの部分積分はエネルギー式と同じである。従って $J'+4\mathcal HJ=k^2[\delta p+(\bar\rho+\bar p)(\Psi-\sigma)]$ を再現する。$\theta$ の式を書く場合は、質量の遷移により $\bar p/\bar\rho$ が変化するので $J$ の分母も微分する。

### 自由流ニュートリノによるテンソルの減衰

相対論的な無衝突種が初期に等方で、テンソル波に後から応答する場合を取る。その角度輝度 $F=4\Theta$ は $F'+ik\mu F=-2n^in^jh_{ij}'$ に従う。プラス偏光の応力を $\Pi_+=\rho_\nu\langle (e^+_{ij}n^in^j)F\rangle_{\vb n}$ とする。方位角の積分で $(e^+_{ij}n^in^j)^2$ の平均は $(1-\mu^2)^2/4$ であり、自由流の積分解を使うと
$$
\Pi_+(\eta)=-4\rho_\nu(\eta)\int_{\eta_i}^{\eta}
\frac{j_2(k(\eta-u))}{[k(\eta-u)]^2}h_+'(u)\dd u .
$$
係数は $\int_{-1}^1(1-\mu^2)^2e^{-i\mu v}\dd\mu=16j_2(v)/v^2$ を用いた。これは左辺を $[1+2\partial_v^2+\partial_v^4](2\sin v/v)$ として微分すれば得られる。核は $v=0$ で $1/15$ に近づくため、上端に発散はない。

テンソルの Einstein 式に入れ、$\rho_\nu=f_\nu\rho_{\rm tot}$、平坦背景の $8\pi G a^2\rho_{\rm tot}=3\mathcal H^2$ を使うと
$$
h_+''+2\mathcal Hh_+'+k^2h_+
=-24 f_\nu\mathcal H^2\int_{\eta_i}^{\eta}
\frac{j_2(k(\eta-u))}{[k(\eta-u)]^2} h_+'(u)\dd u .
$$
クロス偏光でも方位角平均は同じである。応力には過去の時間微分が残る。したがって定数の局所粘性を加えるだけでは全時刻を再現できない。無衝突であることは応力が零であることを意味しない。地平線に進入して $h'$ が生じると、方向ごとに異なるエネルギーへのずれが四重極を作る。

### Fermi 分布と相対論的な背景密度

完全に脱結合して分布の形が $f_0(q)=1/(e^{q/T_{\nu0}}+1)$ と固定された例を取る。$T_{\nu0}$ は $a=1$ に外挿した温度パラメータであり、質量が効き始めた後の運動エネルギーの熱温度と同じではない。質量を無視すれば $\rho_\nu=g T_{\nu0}^4/(2\pi^2a^4)\int_0^\infty x^3/(e^x+1)\dd x$ である。$1/(e^x+1)=\sum_{n=1}^\infty(-1)^{n-1}e^{-nx}$ を積分すると $6(1-2^{-3})\zeta(4)$ になる。対応する Bose 積分は $6\zeta(4)$ なので、自由度と温度が同じなら相対論的な密度は $7/8$ 倍である。電子陽電子の消滅を瞬時の断熱加熱と近似した温度比 $T_\nu/T_\gamma=(4/11)^{1/3}$ を使うと、一味のニュートリノと反ニュートリノを光子の二偏光と比べた密度比は $(7/8)(4/11)^{4/3}$ となる。エネルギー密度をこの比の何倍かで表すのが $N_{\rm eff}$ の規約である。質量の効果が大きい時代にはこの密度を $a^{-4}$ で外挿せず、運動量積分を更新する。

### 放射期に進む自由流の距離

放射優勢で $a=C\eta$、$C>0$ とおき、質量を持つ粒子の運動量ビン $q>0$ を選ぶ。$\eta_{\rm nr}=q/(mC)$ とすれば $q/\varepsilon_q=[1+(\eta/\eta_{\rm nr})^2]^{-1/2}$ である。従って自由流の共動距離は
$$
d_q(\eta)=\eta_{\rm nr}\big[
\operatorname{arsinh}(\eta/\eta_{\rm nr})-
\operatorname{arsinh}(\eta_i/\eta_{\rm nr})\big].
$$
初期には $d_q\simeq\eta-\eta_i$、後期には $\log(2\eta/\eta_{\rm nr})$ に従って増える。固定の遅い時刻でも、運動量が大きいビンはより長く相対論的に動く。このため同じ種の自由流を代表速度だけで計算すると、正確な密度や四重極には誤差が残る。

:::example[相対論的から非相対論的への遷移]
$\eta_i=0$、$\eta=\eta_{\rm nr}$ で $d_q/\eta_{\rm nr}=\log(1+\sqrt2)\simeq0.881$ である。光ならこの比は $1$ である。瞬時に速度を切り替える模型を使わなくても、放射期に限れば各ビンの実際の速度から距離を求められる。
:::

### 非相対論的な応力の積分と桁の比較

運動量単位 $T_{\nu0}$ を使い、$z=am/T_{\nu0}$、$u=q/T_{\nu0}$ とする。固定された Fermi 分布では背景の数密度が $a^{-3}$ に比例し、圧力密度比は
$$
w(z)=\frac{\displaystyle\int_0^\infty
       \frac{u^4}{3\sqrt{u^2+z^2}}\frac{\dd u}{e^u+1}}
     {\displaystyle\int_0^\infty
       u^2\sqrt{u^2+z^2}\frac{\dd u}{e^u+1}}
$$
である。$z\to0$ では分子の積分が分母の三分の一である。$z\gg1$ では分母の $\sqrt{u^2+z^2}$ が先頭で $z$、分子の逆数が $z^{-1}$ なので
$$
w(z)\simeq\frac{1}{3z^2}
\frac{15\zeta(5)}{\zeta(3)} .
$$
分子と分母を別々に展開しているため、寄与する $u$ が $z$ より十分小さくなるまでこの近似は使わない。相対論から非相対論への遷移は、全粒子が同じ運動量で急に切り替わる現象ではない。例えば $\langle q\rangle/T_{\nu0}=\int u^3 f/\int u^2 f=7\pi^4/[180\zeta(3)]\simeq3.151$ である一方、二乗平均の平方根は $\sqrt{15\zeta(5)/\zeta(3)}\simeq3.597$ である。密度と応力で異なる平均速度が現れるため、精密計算で単一の代表運動量を使うと差が出る。

多重極の数値積分は、背景の積分と同じ運動量点を使うと保存式で照合しやすい。$u\gg1$ の尾は指数で抑えられるが、$\delta p$ と $\sigma$ には $\delta\rho$ より高い $u$ の重みがある。積分の上端を変えて背景密度が収束しても、応力の収束は別に調べる。また方程式の $\varepsilon_q/q$ は $q=0$ で形式上発散するが、熱分布の $L_q=q(\partial_q f_0)/f_0$ が $q$ に比例する。先に積に直せば重力源は有限の $\varepsilon_q(\partial_qf_0)/f_0$ となる。小さい $q$ で二係数をばらばらに評価して桁を失うことを避けられる。

:::example[相対論的ニュートリノが担う放射の比]
理想化した三味で $T_\nu/T_\gamma=(4/11)^{1/3}$ とし、他の放射を光子だけとする。$\rho_\nu/\rho_\gamma=3(7/8)(4/11)^{4/3}\simeq0.681$ であり、放射全体に対する分率は $R_\nu=(\rho_\nu/\rho_\gamma)/(1+\rho_\nu/\rho_\gamma)\simeq0.405$ である。この二つを区別しないと、跡なし Einstein の初期条件に使う分率を過大評価する。質量や脱結合の非瞬時性を入れれば、この値をそれぞれの運動量分布で更新できる。
:::

## 断熱的初期条件

:term[断熱的初期条件]{en="adiabatic initial conditions" yomi="だんねつてきしょきじょうけん"} では全成分が同じ局所時間ずれから密度変化を受ける。従って $\delta_i/(1+w_i)$ は共通で、$\delta_b=\delta_c=3\delta_\gamma/4=3\delta_\nu/4$ である。放射期の正則な解で $\mathcal H=1/\eta$、光子とバリオンが強く結合し $R_b\ll1$、$k\eta\ll1$ とし、先頭で $\Phi'=\Psi'=0$ とする。:eqref[eq-cosmo-einstein-scalar] の時間制約により $\delta_\gamma=\delta_\nu=-2\Psi$ となる。圧力を含む放射の速度式は $\theta'=k^2(\Psi+\delta_\gamma/4)+\order{k^4\eta^2}=k^2\Psi/2+\cdots$ だから $\theta=k^2\eta\Psi/2$ である。冷たい物質は $\theta'+\theta/\eta=k^2\Psi$ に従い、正則解が同じ先頭を持つ。ニュートリノの自由流の四重極では $\sigma_\nu'=4\theta_\nu/15+\order{k^4\eta^3}$ なので $\sigma_\nu=(k\eta)^2\Psi/15$ になる。放射のニュートリノ比を $R_\nu$ とすると Einstein の跡なし制約より $\Phi-\Psi=(2R_\nu/5)\Psi$ である。

:::example[密度と速度の次数]
放射密度は $k\eta$ の零次だが $\theta/k$ は一次、ニュートリノの四重極は二次である。
:::


:::example[バリオンの初期過密]
$\Psi=-10^{-5}$ なら $\delta_b=1.5\times10^{-5}$ である。
:::

### 初期振幅を原始曲率に合わせる

放射優勢で断熱な正則モードについて $\delta_r=-2\Psi$、$\Phi=(1+2R_\nu/5)\Psi$ を使う。$\zeta=-\Phi+\delta_r/4$ だから
$$
\zeta=-\left(\frac32+\frac{2R_\nu}5\right)\Psi .
$$
従って与えた原始曲率から初期の $\Psi,\Phi,\delta_i$ を順に定められる。$R_\nu$ が同じでも断熱条件を破って独立な成分の $\zeta_i$ を指定すれば、同じ正規化式は使えない。

:::example[自由流成分を無視した場合との比較]
$R_\nu=2/5$ では係数の絶対値は $83/50$ である。$R_\nu=0$ では $3/2$ なので、同じ $\zeta$ での $|\Psi|$ は $75/83$ 倍になる。固定する初期量を $\Psi$ にするか $\zeta$ にするかで比較の結果が変わる。
:::

### 曲率と成分の相対初期値

独立に保存する成分について、断熱性は $\zeta_i=-\Phi+\delta_i/[3(1+w_i)]$ が共通になる条件である。一方、平坦な放射優勢期に $\zeta_r=0$ と $S_{cr}=3(\zeta_c-\zeta_r)=S$ を指定しても線形初期データを作れる。これを :term[等曲率摂動]{en="isocurvature perturbation" yomi="とうきょくりつせつどう"} と呼ぶ。ここでは光子だけの完全放射と圧力なしの冷たい暗黒物質を取り、$y=\bar\rho_c/\bar\rho_r\ll1$、$k\eta\ll1$ とする。$\Phi=\Psi$、$\delta_r=4\Phi$、$\delta_c=S+3\Phi$ である。全密度の相対摂動は先頭で $4\Phi+yS$ となる。従って時間制約は
$$
\eta\Phi'+3\Phi=-yS/2 .
$$
$y\propto a\propto\eta$ より正則解は $\Phi=-yS/8$ である。同次の $\eta^{-3}$ は初期に発散するので除く。冷たい物質に有限の密度差があっても、初期の重力源はそのエネルギー割合 $y$ で抑えられる。後に物質が増えると全曲率が変わることは、:eqref[eq-cosmo-zeta-conservation] と各成分のエンタルピー平均でも確認できる。自由流ニュートリノを含む場合は、その応力の制約と初期階層も同時に解く。

### 初期条件を数値に落とすときの次数

$\Phi,\Psi$ の先頭は定数、密度は $\order{1}$、速度は $\theta/k=\order{k\eta}$ である。無衝突ニュートリノではこの速度から $\sigma_\nu=\order{(k\eta)^2}$ が出る。$k^2$ で割った Einstein の跡なし制約にこの小さい量を入れるため、$\Phi-\Psi$ は先頭から有限になり得る。一方の多重極を零に初期化しながら両ポテンシャルの差だけ残すと制約に反する。積分開始を早めたときに同じ物理時刻の解が一致するか、時間・運動量・多重極の切断をそれぞれ変えたときに低次の量が収束するかを検証する。

:::example[初期時刻の選択]
$k_{\max}\eta_i=0.02$ なら $(k_{\max}\eta_i)^2=4\times10^{-4}$ である。零次の密度を使うならこの二次の補正が誤差の一つになる。ただし自由流の四重極まで省くと、跡なし制約ではその誤差が $(k\eta_i)^2$ に抑えられるとは限らない。
:::

### 初期の高次項を純放射模型で検査する

光子だけの完全放射、$R_b=R_\nu=0$、平坦背景で $\Phi=\Psi$ とする。定数振幅を $\Phi_i$、$x=k\eta/\sqrt3$ とすれば $\Phi=3\Phi_i(\sin x-x\cos x)/x^3$ である。展開は
$$
\Phi=\Phi_i(1-x^2/10+x^4/280+\order{x^6}).
$$
時間 Einstein 制約から $\delta_\gamma=-2(x^2\Phi+x\partial_x\Phi+\Phi)$ なので $\delta_\gamma=-2\Phi_i[1+7x^2/10-23x^4/280+\order{x^6}]$ である。速度は運動量制約の $4\pi Ga^2(\rho+p)=2/\eta^2$ を使うと
$$
\theta_\gamma=\frac{k^2\eta}{2}(\Phi+x\partial_x\Phi)
=\frac{k^2\eta\Phi_i}2[1-3x^2/10+\order{x^4}]
$$
となる。初期化で $\Phi'$ を零に置くのは零次のポテンシャルに対する近似であり、二次まで一致させるなら $\Phi'$ も一貫して保持する。この模型を解析解と照合すれば、二つの制約の符号と初期時刻を変えた場合の収束を検査できる。ニュートリノなどのある系にこの二次係数をそのまま流用しない。

### 補償された成分間の初期条件

バリオンと冷たい暗黒物質に $\bar\rho_b\delta_b+\bar\rho_c\delta_c=0$ を課し、光子とニュートリノの揺らぎを零にする。全ての速度、異方性応力、計量も零とする。バリオン圧力を無視し、背景の Thomson 率を用いる線形方程式では、重力源は零で速度も生成されず、二つの成分の密度差は残る。従って原始曲率の指定だけで全ての線形初期条件が尽くされるわけではない。有限のバリオン圧力や、局所のバリオン密度が再結合と音速を変える非線形効果を考えると、補償の観測上の見え方が変わる。

:::example[質量分率で重み付けする補償]
$\bar\rho_b/\bar\rho_c=1/5$、$\delta_b=5\times10^{-4}$ とすると、補償する冷たい暗黒物質の $\delta_c=-10^{-4}$ である。分数摂動を単純に $\delta_c=-\delta_b$ とすると、少ない側と多い側の背景密度を同じ重みで足してしまう。
:::

### 自由流する放射を含む断熱解の二次

純粋な放射背景で物質とバリオン慣性を無視する。ただし一定比 $R_\nu$ の無衝突ニュートリノは残す。$x=k\eta$、$\Phi=\Phi_0+\Phi_2 x^2+\cdots$、$\Psi=\Psi_0+\Psi_2x^2+\cdots$ とおく。光子とニュートリノの同じ先頭速度 $\theta/k=\Psi_0 x/2$ を連続式に入れると、両者の密度は $-2\Psi_0+(4\Phi_2-\Psi_0/3)x^2+\cdots$ である。時間制約の $x^2$ の係数は
$$
12\Phi_2+3\Psi_2=-\Phi_0+\Psi_0/2 .
$$
ニュートリノについて $\theta_\nu/k=\Psi_0 x/2+t_3x^3+\cdots$、$\sigma_\nu=s_2x^2+s_4x^4+\cdots$、無次元輝度の三重極を $F_{\nu3}$ とする。速度式で $s_2=\Psi_0/15$ を使えば $3t_3=\Phi_2+\Psi_2-3\Psi_0/20$ である。階層より $\partial_x F_{\nu3}=6s_2x^2/7+\cdots$、従って $F_{\nu3}=2s_2x^3/7+\cdots$ となる。四重極式 $\partial_x\sigma_\nu=(4/15)(\theta_\nu/k)-(3/10)F_{\nu3}$ は $s_4=t_3/15-3s_2/140$ を与える。整理すると
$$
s_4=\frac{\Phi_2+\Psi_2}{45}-\frac{\Psi_0}{210},\qquad
\Phi_2-\Psi_2=6R_\nu s_4 .
$$
後式は跡なし制約の二次である。従って時間制約と合わせて $2\times2$ の連立式で $\Phi_2,\Psi_2$ を決められる。純粋な光子流体なら $\Phi_0=\Psi_0$、$\Phi_2=\Psi_2$ から $\Phi_2=-\Psi_0/30$ であり、前の解析解と一致する。自由流を入れると二次のポテンシャルを得るために四重極の四次までが必要であり、上端の多重極を切る時にこの次数の違いに注意する。

### 自由流の割合を変えた初期値

二次の断熱系に $r=R_\nu$ と $\Phi_0=(1+2r/5)\Psi_0$ を入れる。未知の $(\Phi_2,\Psi_2)$ に対し
$$
4\Phi_2+\Psi_2=-\left(\frac16+\frac{2r}{15}\right)\Psi_0,\qquad
\left(1-\frac{2r}{15}\right)\Phi_2-
\left(1+\frac{2r}{15}\right)\Psi_2=-\frac r{35}\Psi_0 .
$$
右辺の第二式の定数は四重極の $s_4$ に含まれる三重極の戻りを含む。$r$ を固定し、第一式で $\Psi_2$ を消すと
$$
\frac{\Phi_2}{\Psi_0}
=-\frac{(1+2r/15)(1/6+2r/15)+r/35}{5+2r/5}.
$$
$r=0$ では $-1/30$ に戻り、$r=2/5$ では約 $-0.04712$ となる。$\Psi_2/\Psi_0$ は $-1/6-2r/15-4\Phi_2/\Psi_0$ で計算し、$r=2/5$ では約 $-0.03150$ である。固定した原始曲率を比較するなら $\Psi_0/\zeta$ も $r$ とともに変わる。先頭のポテンシャル差を合わせただけでは、有限の初期時刻における $k\eta$ 二次の差まで合わない。

### 共形時刻と観測波数を合わせる

放射優勢で $a=H_0\sqrt{\Omega_{r0}}\eta$ と規格化できる。実際、$a'/a=aH$ と $H=H_0\sqrt{\Omega_{r0}}/a^2$ から $a'=H_0\sqrt{\Omega_{r0}}$ である。自然単位で $H_0$ を共動長さの逆数として扱うため、観測単位を使うなら $H_0/c$ とする。例えば $H_0/c=2.25\times10^{-4}\,\mathrm{Mpc^{-1}}$、$\Omega_{r0}=9\times10^{-5}$ とすれば、$a=10^{-8}$ は $\eta\simeq0.00469\,\mathrm{Mpc}$ である。$k=1\,\mathrm{Mpc^{-1}}$ は $k\eta\simeq0.00469$ なので、先頭の長波長展開が使える。早すぎる時刻では、振動のためというより衝突の時間尺度が短すぎて陽的な積分が難しくなる。波数帯全体で初期化をそろえる際は、最大波数での $k\eta$ と最小平均自由時間の両方を確認する。

## 演習問題

::::exercise{#exr-cosmo-cosmological-boltzmann-1 level=basic}
$\Theta_0=10^{-5}$ の光子密度比を求めよ。
:::solution
放射密度の温度四乗を線形化して $\delta_\gamma=4\times10^{-5}$ となる。
:::
::::

::::exercise{#exr-cosmo-cosmological-boltzmann-2 level=advanced}
初期が等方な自由流の単極子を直接求めよ。
:::solution
$x=k\Delta\eta$ とし、$A e^{-ix\mu}$ を方向平均すると $A\sin x/x$ である。
:::
::::

::::exercise{#exr-cosmo-cosmological-boltzmann-3 level=basic}
$R_b=2$ の場合に速度差のみの減衰係数を散乱率単位で求めよ。
:::solution
光子とバリオンの式の差を取れば $1+1/R_b=3/2$ である。
:::
::::

::::exercise{#exr-cosmo-cosmological-boltzmann-4 level=advanced}
Thomson 核の全立体角における規格化を証明せよ。
:::solution
$2\pi$ を掛けた積分は $(3/8)\int_{-1}^1(1+\mu^2)\dd\mu=(3/8)(8/3)=1$ である。
:::
::::

::::exercise{#exr-cosmo-cosmological-boltzmann-5 level=basic}
$q=3am$ の粒子速度を求めよ。
:::solution
速度 $q/\sqrt{q^2+a^2m^2}$ に代入すると $3/\sqrt{10}$ になる。
:::
::::

::::exercise{#exr-cosmo-cosmological-boltzmann-6 level=advanced}
速度差の衝突項で光子・バリオンの運動量の和が保存されることを示せ。
:::solution
光子の慣性を一単位、バリオンを $R_b$ 単位とすれば $(\theta_b-\theta_\gamma)+R_b(\theta_\gamma-\theta_b)/R_b=0$ となる。
:::
::::

::::exercise{#exr-cosmo-cosmological-boltzmann-7 level=basic}
断熱条件で $\delta_\gamma=8\times10^{-5}$ の $\delta_c$ を求めよ。
:::solution
物質の分数摂動は放射の $3/4$ で $6\times10^{-5}$ である。
:::
::::

::::exercise{#exr-cosmo-cosmological-boltzmann-8 level=advanced}
$R_\nu=0$ の放射優勢で $\zeta/\Psi$ を求めよ。
:::solution
$\Phi=\Psi$ かつ $\delta_\gamma=-2\Psi$ なので $\zeta=-\Phi+\delta_\gamma/4=-3\Psi/2$ である。
:::
::::

::::exercise{#exr-cosmo-cosmological-boltzmann-9 level=advanced}
無衝突ニュートリノの運動量積分で $q^3\varepsilon_q f_0(q)$ の端点が零とする。Euler 式の重力源が、粒子が非相対論的になった後も $(\bar\rho+\bar p)k^2\Psi$ となることを示せ。
:::solution
双極子の重力項は $-k(\varepsilon_q/q)(\dd\log f_0/\dd\log q)\Psi/3$ である。運動量密度を $Ak\int q^3 f_0 \mathcal F_1\dd q$ と測れば、この重力項全体は $-(A k^2\Psi/3)\int q^3\varepsilon_q\partial_q f_0\dd q$ となる。部分積分し、$\partial_q\varepsilon_q=q/\varepsilon_q$ を使うと $A\int (3 q^2\varepsilon_q+q^4/\varepsilon_q)f_0\dd q=3(\bar\rho+\bar p)$ を得る。これで重力係数は述べた値になる。途中で $\varepsilon_q=q$ に置き換えていないため、どの時刻でも成り立つ。
:::
::::

::::exercise{#exr-cosmo-cosmological-boltzmann-10 level=advanced}
偏光を含む散乱だけの系で初期 $\Theta_2=0$、$P_0^{\rm pol}=B$、$P_2^{\rm pol}=0$ とする。$\Theta_2(u_T)$ と長時間の主な減衰率を求めよ。
:::solution
$\Pi=B e^{-3u_T/10}$ であり、$\dd\Theta_2/\dd u_T+\Theta_2=\Pi/10$ を積分する。初期値が零なので $\Theta_2=(B/7)(e^{-3u_T/10}-e^{-u_T})$ である。短時間の傾きは $B/10$、長時間では $e^{-3u_T/10}$ が遅いので率は $(3/10)\dot\kappa_T$ である。温度の四重極が初めからなくても、偏光から散乱利得を経由して作られる。
:::
::::

::::exercise{#exr-cosmo-cosmological-boltzmann-11 level=advanced}
自由流の初期有効温度として $\ell=2$ の振幅 $Q$ だけを与える。波数 $k$ で短い時間 $\Delta\eta$ が過ぎたときの単極子の先頭項を、階層から求めよ。
:::solution
初期では $X_1=0$、$X_2=Q$、$X_0=0$ なので、$\Theta_1$ と同じ双極子の式は $X_1'=k(X_0-2X_2)/3=-2kQ/3$ である。次に $X_0'=-kX_1$ なので、初期で $X_0'=0$、$X_0''=2k^2Q/3$ を得る。従って $X_0=Q(k\Delta\eta)^2/3+\cdots$ である。輸送が一階しか結ばないため、四重極は単極子へ届くまで二つの隣接モードを経る。
:::
::::

## 確認テスト

::::quiz
光子の自由流が結ぶ多重極を選べ。

- [x] 隣接する多重極
- [ ] 同じ多重極だけ
- [ ] 常に単極子だけ

:::explanation
$\mu P_\ell$ が $P_{\ell-1}$ と $P_{\ell+1}$ の線形結合なので隣り合う次数が結合する。
:::
::::

::::quiz{answer="2" tol="1e-6"}
$R_b=1$ で散乱だけによる速度差の減衰率を $\dot\kappa_T$ 単位で答えよ。
:::explanation
速度差の係数は $1+R_b^{-1}=2$ である。
:::
::::

::::quiz{answer="1.16" tol="1e-6"}
$R_\nu=0.4$ の先頭初期値で $\Phi/\Psi$ は何か。
:::explanation
異方性応力から $1+2R_\nu/5=1.16$ となる。
:::
::::

::::quiz
断熱初期条件の関係を選べ。

- [x] $\delta_b=3\delta_\gamma/4$
- [ ] $\delta_\gamma=0$ が必要。
- [ ] $\delta_b=4\delta_\gamma/3$

:::explanation
成分間で $\delta_i/(1+w_i)$ が揃い、圧力なしと放射の比は $1/(4/3)=3/4$ である。
:::
::::

---
title: 宇宙マイクロ波背景放射
description: 背景放射の黒体性と異方性を、重力、音響、拡散、散乱偏光から説明する。
notation:
  - symbol: 'x'
    meaning: 光子エネルギーと温度の比
  - symbol: '\Theta'
    meaning: 温度の一次摂動
  - symbol: 'r_s'
    meaning: 共動音響地平線
  - symbol: 'k_D'
    meaning: 拡散減衰波数
  - symbol: 'C_\ell'
    meaning: 角度パワースペクトル
---

## この章の目標

- 黒体形状と小さな歪みを区別する。
- 重力赤方偏移と音響振動の異方性を導く。
- 温度と偏光の相関から宇宙のパラメータを測る方法を説明する。

## 黒体スペクトルとスペクトル歪み

$c=\hbar=k_B=1$、計量符号 $(-,+,+,+)$ の平坦背景とする。光子の局所エネルギーを $p$、温度を $T$ とすれば、:ref[def-cosmo-planck-function] に対応する占有数は $n(x)=1/(e^x-1)$、$x=p/T$ である。自由流では $p$ と $T$ がともに $a^{-1}$ でずれるため、この形は保たれる。現在の単極子温度は約 $2.725\,\mathrm K$ である。[絶対温度の測定](https://doi.org/10.1088/0004-637X/707/2/916)

Compton 散乱だけが速く光子数を変える過程が遅い区間では、エネルギーを移すと分布は $n=1/[\exp(x+\mu_{\rm d})-1]$ へ近づき得る。この $\mu_{\rm d}$ は無次元歪みで、熱力学的化学ポテンシャルとは符号が逆である。二重 Compton・自由–自由放射が十分速ければ光子数も再調整される。

薄い高温電子の弱い散乱の $y$ 歪みは $\Delta n=y x^{-2}\partial_x[x^4\partial_xn]$ と表せる。軟光子が小さい変化を多数受ける運動量拡散の Kompaneets 作用素 :eqref[eq-cosmo-kompaneets] を黒体に適用し、$T_e\gg T_\gamma$ で残る項である。数の変化は $\int x^2\Delta n\dd x$ が全微分なので零、エネルギーの変化は $\int x^3\Delta n\dd x=-y\int x^4 n'\dd x=4y\int x^3 n\dd x$ である。精密な生成時期と振幅は光子数変化の衝突積分を要し :cite[weinberg-cosmology] を参照する。

:::example[温度変化と歪み]
全体の温度が二倍になっても Planck 形状は変わらず、全エネルギー密度は十六倍になる。
:::


:::example[弱い散乱]
$y=10^{-6}$ なら光子エネルギー密度は一次で $4\times10^{-6}$ 増す。
:::

### 光子数を固定した Bose–Einstein 歪み

光子エネルギーの再分配は速いが、光子数を変える反応が遅いとする。占有数を $n=[\exp(h\nu/(k_BT)+\mu)-1]^{-1}$ と書く。ここで $\mu$ は通常の化学ポテンシャルを $-k_BT$ で割った値である。基準黒体と比べ $T/T_0=1+t$、$|\mu|,|t|\ll1$ とする。積分 $\int_0^\infty x^j/(\exp(x+\mu)-1)\dd x=j!\sum_{n\ge1}e^{-n\mu}/n^{j+1}$ を $\mu\ge0$ で一次展開すると
$$
\frac{\Delta n_\gamma}{n_\gamma}=3t-\frac{\zeta(2)}{\zeta(3)}\mu,\qquad
\frac{\Delta\rho_\gamma}{\rho_\gamma}=4t-\frac{\zeta(3)}{\zeta(4)}\mu .
$$
数を保存すれば $t=\zeta(2)\mu/[3\zeta(3)]$ であり、消去して $\mu\simeq1.401\,\Delta\rho_\gamma/\rho_\gamma$ となる。有限周波数で光子生成も有効な場合、この一つの $\mu$ だけでは歪み全体を記述できない。

:::example[小さな熱の注入]
光子数の変化なしに $\Delta\rho_\gamma/\rho_\gamma=10^{-5}$ なら $\mu\simeq1.40\times10^{-5}$、$t\simeq6.39\times10^{-6}$ である。温度を上げるだけで同じエネルギーを加えれば光子数も増えるので、数固定の条件は満たせない。
:::

### 熱的 Sunyaev–Zeldovich 歪みの零点

高温で薄い電子層に対し、電子温度より低い温度の黒体を背景とする。入射黒体の温度で $x=h\nu/(k_BT_\gamma)$ を定め、:eqref[eq-cosmo-kompaneets] の高温電子による拡散項だけを短時間積分すると
$$
\Delta n=\frac{y}{x^2}\dv{}{x}\!\left(x^4\dv{n_0}{x}\right),\qquad
y=\int n_e\sigma_T\frac{k_BT_e}{m_ec^2}\dd l .
$$
ここで $\dd l$ は電子系の固有距離、$T_e\gg T_\gamma$ とした。黒体の $n_0=(e^x-1)^{-1}$ を微分して
$$
\Delta n=y\,\frac{xe^x}{(e^x-1)^2}[x\coth(x/2)-4]
$$
を得る。低周波で $x\coth(x/2)\to2$ なので $\Delta n/n_0\to-2y$ であり、高周波で括弧が正になる。$g(x)=x\coth(x/2)$ を微分すると、符号は $\sinh x-x>0$ に一致するため、$x>0$ で $g$ は狭義増加である。$g(0)=2$ と $g(x)\to\infty$ から零点は一つである。

:::example[減光と増光の境界]
零点 $x\simeq3.830$ に $T_\gamma=2.725\,\mathrm K$、$k_B/h\simeq20.84\,\mathrm{GHz\,K^{-1}}$ を使うと $\nu\simeq217.5\,\mathrm{GHz}$ である。電子が相対論的になれば拡散係数と反跳の高次項が必要になり、境界の周波数も変化する。
:::

### 異なる温度の黒体を混ぜる

基準温度 $T_0$ に対する対数温度差を $t=\log(T/T_0)$、$D=-x\partial_x$ とする。各方向の占有数は $n_0(xe^{-t})=n_0+tDn_0+t^2D^2n_0/2+\order{t^3}$ である。視野内で平均を取り、平均対数温度の黒体との差を取ると $\operatorname{Var}(t)D^2n_0/2$ が残る。$D^2n_0=x^2n_0''+xn_0'$、$Y=x^{-2}\partial_x(x^4\partial_xn_0)=x^2n_0''+4xn_0'$ なので $D^2n_0=Y+3Dn_0$ である。従ってこの混合には対数温度をさらに $3\operatorname{Var}(t)/2$ だけ変える項と、$y_{\rm mix}=\operatorname{Var}(t)/2$ の歪みがある。平均対数温度を零に取ると、混合した光子数の比は $\langle e^{3t}\rangle=1+9\operatorname{Var}(t)/2+\order{t^3}$ である。これは温度を $3\operatorname{Var}(t)/2$ だけ移した黒体と一致し、歪み自体が光子数を変えないことを確かめられる。

:::example[二方向を同じ重みで足す]
対数温度差を $\pm10^{-3}$ とすれば $y_{\rm mix}=5\times10^{-7}$、温度の二次の対数補正は $1.5\times10^{-6}$ である。各方向で完全な黒体でも、角度平均では歪みが残る。
:::

### 熱化に必要な時間尺度

電子との散乱頻度は $\Gamma_T=n_e\sigma_Tc$ である。非相対論的な熱電子との一回の散乱で光子の典型的な相対エネルギー変化を決めるのは $k_BT_e/(m_ec^2)$ の尺度である。従って角度を再配分する条件 $\Gamma_T/H\gg1$ と、エネルギーを再配分する条件は異なる。後者の目安は $\int [k_BT_e/(m_ec^2)]\Gamma_T\dd t$ である。本章の $y$ はエネルギーの相対増分が小さいとき $4y$ となる規約であり、放射過程の章の平均相対増幅を含む Compton パラメータとは係数 $4$ が異なる。

エネルギーの再配分が進んでも、単一 Compton 散乱は光子数を変えない。二重 Compton と自由–自由放射は光子数の緩和を担い、とくに低周波で化学ポテンシャルを抑える。そのため注入が十分早ければ主に黒体温度の変化、光子数をほぼ保ったまま再分配が進めば主に $\mu$ 型、再分配が少なければ主に $y$ 型に移る。時間帯の境界は連続で、遷移領域ではこの二つの歪みの単純な和に残差が生じ得る。

### 観測する比強度と熱力学的温度

周波数の単位を復元し、$h_{\rm P}=2\pi\hbar$ と書く。両偏光を含む比強度は $I_\nu=(2h_{\rm P}\nu^3/c^2)n$ である。これは等方な場で単位体積・単位周波数のエネルギー密度が $4\pi I_\nu/c$ となり、モード数 $8\pi\nu^2\dd\nu/c^3$ にエネルギー $h_{\rm P}\nu$ と占有数を掛けることから得る。基準黒体の周りでは $x=h_{\rm P}\nu/(k_BT_0)$ とおき、
$$
\delta I_\nu=\frac{2h_{\rm P}\nu^3}{c^2}
\frac{xe^x}{(e^x-1)^2}\frac{\delta T}{T_0}
$$
である。この変換で同じ $\delta T/T_0$ が全周波数で得られるのは黒体の温度摂動である。低周波の Rayleigh–Jeans 温度を任意の周波数で用いると、温度のスペクトルに見かけの形が付く。

:::example[二つの温度単位]
$\delta T_{\rm RJ}=c^2\delta I_\nu/(2k_B\nu^2)$ なら $\delta T_{\rm RJ}/\delta T=x^2e^x/(e^x-1)^2$ である。$x\ll1$ で比は $1$、$x=3$ では約 $0.496$ へ低下する。周波数の異なる測定を比較する場合は、同じ温度単位に直す必要がある。
:::

### 銀河団の運動が作る温度変化

薄い電子雲の全体速度を $\vb v_e$ とし、熱運動はここでは無視する。入射黒体が宇宙の法線枠で等方なら、電子枠では $-\vb n\cdot\vb v_e$ の双極子がある。Thomson 核の利得にこの双極子は入らず、損失はこの双極子を減らす。散乱した光を宇宙枠へ戻すと、温度は $\Delta T/T_0\simeq\tau_T\vb n\cdot\vb v_e$ となる。空の方向 $\hat{\vb s}=-\vb n$ を使えば $-\tau_T\hat{\vb s}\cdot\vb v_e$ である。固有速度で書くと右辺は $c$ で割る。これが運動学的な :term[Sunyaev–Zeldovich 効果]{en="Sunyaev–Zeldovich effect" yomi="すにゃえふぜるどびっちこうか"} である。一次では黒体の温度変化と同じスペクトルを持つため、周波数だけで一次 CMB と分けることはできない。熱的効果は上の $y$ 型の形を持つ。

:::example[遠ざかる銀河団]
$p=\hat{\vb s}\cdot\vb v_e/c=10^{-3}$、$\tau_T=0.01$ なら $\Delta T/T_0=-10^{-5}$、$T_0=2.725\,\mathrm K$ で $\Delta T\simeq-27.3\,\mu\mathrm K$ である。熱的効果はこの結果へ別に加わる。
:::

### 化学ポテンシャルと温度を周波数で見分ける

黒体の小変化を温度と化学ポテンシャルに分け、$t=\delta T/T_0$ とする。$x$ を基準温度で定義し、分布を $[e^{x/(1+t)+\mu}-1]^{-1}$ と展開すると
$$
\Delta n(x)=\frac{e^x}{(e^x-1)^2}(x t-\mu).
$$
数を固定した $\mu$ 型では $t=c_\mu\mu$、$c_\mu=\zeta(2)/[3\zeta(3)]\simeq0.456$ である。従って正のエネルギー注入では $x<1/c_\mu$ で占有数が減り、それより高周波では増える。この近似の零点は $x\simeq2.192$ で、熱的 $y$ 型の零点とは異なる。各周波数で温度に換算すると $\Delta T_{\rm eff}/T_0=t-\mu/x$ なので、同じ $\mu$ 型でも単一の温度に直すことはできない。

小さい $\mu$ の展開は固定された $x>0$ で取る。$x\ll1$ でも $|\mu|\ll x$ が成り立たなければ、黒体に対する線形の相対変化は小さくない。正の $\mu$ では厳密な占有数が $x\to0$ で有限である一方、基準の黒体は $x^{-1}$ に発散する。負の一定 $\mu$ を全ての周波数まで外挿すると零または負の分母が生じる。低周波まで正確に追うときは光子生成も含む周波数依存の衝突式を用いる。数密度・エネルギー密度の線形応答は、必要な積分が収束することを確かめてから取った極限である。

:::example[二つの歪みの零点を比較する]
$T_0=2.725\,\mathrm K$ では $x=1$ は約 $56.8\,\mathrm{GHz}$ である。数を固定した小さい $\mu$ 型の零点は約 $124.5\,\mathrm{GHz}$、熱的 $y$ 型は約 $217.5\,\mathrm{GHz}$ になる。前者では黒体をどの温度とみなすかも歪みの規格化に入っている。温度補正を任意に足せば零点が移るため、比較の際には分解の規約をそろえる。
:::

## Sachs–Wolfe 効果

:term[Sachs–Wolfe 効果]{en="Sachs–Wolfe effect" yomi="ざっくすうぉるふこうか"} は光子の重力赤方偏移が作る温度変動である。:eqref[eq-cosmo-photon-boltzmann] に $\Psi$ の光路全微分を足すと $\dd(\Theta+\Psi)/\dd\eta=\Phi'+\Psi'$（衝突なし）となる。最終散乱を薄い面とし、散乱源の温度と偏光の四重極を無視する。観測者位置の単極子・双極子を除けば

$$
\Theta_o=(\Theta_0+\Psi+\vb n\cdot\vb v_b)_*+\int_{\eta_*}^{\eta_o}(\Phi'+\Psi')\dd\eta . \label{eq-cosmo-sachs-wolfe}
$$

$\vb n$ は光子の進行方向である。右端は積分 Sachs–Wolfe 効果であり、固定のポテンシャル井戸では入射時と出射時の変化が相殺する。物質期の超 Hubble 断熱解では $\delta_m=-2\Psi$、$\Theta_0=\delta_m/3$、$\Phi=\Psi$ なので通常の項は $\Psi/3=-\zeta/5$ となる。

:::example[重力井戸]
$\Psi=-3\times10^{-5}$ の長波長では通常の項は $-10^{-5}$ である。
:::


:::example[時間変化の寄与]
光路に沿って空間勾配を無視でき、各ポテンシャルが $10^{-6}$ 増えた区間では積分項は $2\times10^{-6}$ である。
:::

### 静的な井戸と時間発展する井戸

:eqref[eq-cosmo-sachs-wolfe] の積分に現れるのは、視線に沿った全微分でなく、座標を固定した時間微分である。例えば $\Phi=\Psi=f(\eta)u(\vb x)$ なら積分は $2\int f'(\eta)u(\vb x(\eta))\dd\eta$ である。$u$ の負の領域を通る光にとって、$f'<0$ なら井戸は浅くなり、積分は正となる。$f$ が一定の場合、光が入るときと出るときの静的なエネルギー変化は相殺する。

:::example[広い静的な井戸との比較]
光が通る区間で $u=-1$ とし、$f$ が $10^{-5}$ から $5\times10^{-6}$ に減るとする。積分は $-2(f_{\rm out}-f_{\rm in})=10^{-5}$ である。井戸の外で $f$ が変化しても光路上で $u=0$ ならこの積分には寄与しない。
:::

### 光の進行方向と観測者速度

観測者が見る方向を $\hat{\vb s}$ とすると、到来光の進行方向は $\vb n=-\hat{\vb s}$ である。局所の法線に対し $\vb v_o$ で動く観測者の固有エネルギーは一次で $p_o(1-\vb n\cdot\vb v_o)$ となる。放射体の静止系から法線へは逆向きの変換なので、放射面では $+\vb n\cdot\vb v_b$ が付く。合わせると $\hat{\vb s}\cdot(\vb v_o-\vb v_b)$ である。同じ速度で動く二者の間では一次の相対 Doppler シフトは零になる。視線の向きの取り方を替えるときは、光路と両速度の符号を一緒に変更する。

無摂動の放射温度は $T_\gamma\propto a^{-1}$ なので、放射時刻だけを $\delta\eta$ ずらせば温度は相対的に $-\mathcal H\delta\eta$ だけ変わる。一方、観測までの赤方偏移の比 $a_*/a_o$ は $+\mathcal H\delta\eta$ だけ変わる。背景で両者の積は放射時刻に依存しないので、二つは相殺する。実際の散乱では時間幅にわたって一次の源を積分する。

### 放射・物質の転換による積分項

異方性応力を無視した長波長の断熱解では、放射期と物質期で同じ原始曲率に対する $\Phi/\zeta$ はそれぞれ $-2/3$、$-3/5$ である。従って空間変化が無視できる光路区間がこの転換を丸ごと含む理想化では、$\Phi'+\Psi'=2\Phi'$ の積分は $2(-3/5+2/3)\zeta=2\zeta/15$ である。実際の光は最後の散乱以前には連続的に再散乱するため、転換の全部をこの値で観測することはできない。観測まで無散乱で届く割合 $e^{-\tau}$ と Fourier 位相を掛け、視線に沿って積分する必要がある。

物質が優勢な時期には $\Phi$ の成長モードが一定となる。遅い加速期に物質の集まり方が膨張より遅くなると、Poisson 制約の共動勾配を固定した場合に井戸が浅くなる。従って遅い積分効果は低赤方偏移の密度場と相関する。単一の井戸を横切る例では温度上昇の符号を求められるが、全天の積分には正負の密度、投影、時間変化が全て入る。

### テンソルの重力赤方偏移

テンソルだけの背景摂動 $h_{ij}$ に沿う光子の温度変化は、共動運動量の変化 $q'/q=-h_{ij}'n^in^j/2$ から得られる。衝突なしで初めの温度を等方とする極限では
$$
\Theta_o(\vb k,\vb n)
=-\frac12\int_{\eta_*}^{\eta_0}n^in^jh_{ij}'(\eta,\vb k)
e^{-ik\mu(\eta_0-\eta)}\dd\eta .
$$
波数を第三軸とし、規格化 $e^+_{xx}=1/\sqrt2$、$e^+_{yy}=-1/\sqrt2$ を使えば $n^in^je^+_{ij}=(1-\mu^2)\cos(2\varphi)/\sqrt2$ である。空の Fourier モードを球面調和へ射影すると、方位角積分は $m=\pm2$ のみに残る。例えば $m=2$ では $\int_0^{2\pi}\cos(2\varphi)e^{-2i\varphi}\dd\varphi=\pi$ である。$Y_{\ell2}=N_{\ell2}P_\ell^2(\mu)e^{2i\varphi}$、$N_{\ell2}=\sqrt{(2\ell+1)(\ell-2)!/[4\pi(\ell+2)!]}$ を用い、後の偏光小節で導く角度積分を適用すると
$$
a_{\ell2}^{(+)}(\vb k)
=\frac{\pi N_{\ell2}}{\sqrt2}(-i)^\ell
\frac{(\ell+2)!}{(\ell-2)!}
\int_{\eta_*}^{\eta_0}h_+'(\eta,\vb k)
\frac{j_\ell(k(\eta_0-\eta))}{[k(\eta_0-\eta)]^2}\dd\eta
$$
である。正符号は経路の赤方偏移と角度積分の双方から負符号が一つずつ出ることで決まる。一定の超 Hubble モードに対し $h'=0$ なら、この経路源は零である。実際には光子の散乱による温度と偏光の源も同時に追い、テンソルのパワーを二偏光の規格化と合わせる。

## 音響振動とタイトカップリング近似

:term[タイトカップリング]{en="tight coupling" yomi="たいとかっぷりんぐ"} は散乱率 $\dot\kappa_T$ が波の時間尺度より大きい近似である。:eqref[eq-cosmo-baryon-euler] と光子の双極子式を慣性で足し、$\theta_b=\theta_\gamma=\theta$、バリオン圧力・四重極無視とする。$R_b'= \mathcal H R_b$ から $(1+R_b)\theta'+R_b'\theta=k^2[\Theta_0+(1+R_b)\Psi]$ となる。$\Theta_0'=-\theta/3+\Phi'$ を微分して

$$
\Theta_0''+\frac{R_b'}{1+R_b}\Theta_0'+k^2c_s^2\Theta_0=\Phi''+\frac{R_b'}{1+R_b}\Phi'-\frac{k^2}3\Psi,\qquad c_s^2=\frac1{3(1+R_b)} . \label{eq-cosmo-acoustic}
$$

一定の $R_b,\Psi,\Phi$ では $\Theta_0+(1+R_b)\Psi$ が振動する。可視温度 $\Theta_0+\Psi$ の中点は $-R_b\Psi$ なので圧縮と希薄化のピークが同じ高さにならない。位相は $kr_s$、$r_s(\eta)=\int_0^\eta c_s\dd\eta'$ である。

:::example[バリオンによる音速低下]
$R_b=1$ なら音速は光子だけの値の $1/\sqrt2$ 倍になる。
:::


:::example[位相]
初期速度が無視できる場合、$kr_s=\pi$ で第一の極値、$2\pi$ では振動の中点を基準に逆符号の極値となる。
:::

### 放射期の重力駆動

完全放射流体で異方性応力が零、$R_b=0$ とする。$x=k\eta/\sqrt3$、$\mathcal H=1/\eta$ で :eqref[eq-cosmo-einstein-scalar] の最初の式と Friedmann 式を使うと $\delta_\gamma=-2[x^2\Phi+x\partial_x\Phi+\Phi]$ である。従って
$$
\Theta_0+\Phi=\tfrac12(\Phi-x\partial_x\Phi-x^2\Phi).
$$
放射期の正則解 $\Phi=3\Phi_i(\sin x-x\cos x)/x^3$ を入れる。$x\ll1$ では $\Theta_0+\Phi\to\Phi_i/2$ である。$x\gg1$ では $-x^2\Phi$ が $3\Phi_i\cos x+\order{x^{-1}}$ なので、先頭の振動振幅は $3|\Phi_i|/2$ である。ポテンシャルが音響振動の最中に減るため、初期の有効温度の振幅がそのまま保存されない。

:::example[同じ初期振幅での音響利得]
$\Phi_i=-2\times10^{-5}$ とする。初期の有効温度は $-10^{-5}$ で、後の深い亜地平線では $\cos x=1$ の近くで約 $-3\times10^{-5}$ である。バリオン慣性、ニュートリノ応力、物質成分を入れるとこの係数は変わる。
:::

### 音響地平線の解析積分

$R_b=R_0 a$、$R_0=3\Omega_{b0}/(4\Omega_{\gamma0})$ とする。平坦な放射・物質背景で $a_{\rm eq}=\Omega_{r0}/\Omega_{m0}$、$H=H_0a^{-2}\sqrt{\Omega_{r0}+\Omega_{m0}a}$ を取る。$\dd\eta=\dd a/(a^2H)$ を音響積分へ代入し、$R_{\rm eq}=R_0a_{\rm eq}$ とおけば
$$
r_s(a)=\frac{1}{H_0\sqrt{3\Omega_{m0}R_0}}
\int_0^{R_b(a)}\frac{\dd R}{\sqrt{(R+R_{\rm eq})(R+1)}}
=\frac{2}{H_0\sqrt{3\Omega_{m0}R_0}}
\log\frac{\sqrt{R_b(a)+R_{\rm eq}}+\sqrt{R_b(a)+1}}{\sqrt{R_{\rm eq}}+1}.
$$
積分は対数の引数の分子を微分して得られる。SI 単位で $H_0$ を代入するときは分子を $c$ 倍する。$\Omega_{r0}$ は膨張に寄与する相対論的種を含み、$R_0$ の分母は光子だけである。自由流ニュートリノは背景の膨張に寄与するが、光子とバリオンの音波の圧力と慣性には参加しない。

重力を外した同次方程式で $u=\sqrt{1+R_b}\Theta_0$ とおく。$\gamma=R_b'/(1+R_b)$ を用いると $u''+[k^2c_s^2-\gamma'/2-\gamma^2/4]u=0$ となる。$kc_s$ が $\gamma,\sqrt{|\gamma'|}$ より大きく、周波数が一周期にほぼ一定なら、$u$ の WKB 振幅は $(kc_s)^{-1/2}$、位相は $k\int c_s\dd\eta$ である。従って $\Theta_0$ の振幅は $(1+R_b)^{-1/4}$ に比例する。これは散逸がなくても、膨張による慣性の変化で生じる。

### 圧縮、希薄化、速度の位相

定数 $R_b$ と静的なポテンシャル、初期速度が零の近似で、$X=\Theta_0+\Psi$ は
$$
X=-R_b\Psi+[X_i+R_b\Psi]\cos(kc_s\Delta\eta).
$$
物質期の長波長の断熱値 $X_i=\Psi/3$ を仮に接続する。井戸に対応する $\Psi<0$ では、最初の圧縮 $kc_s\Delta\eta=\pi$ で $X=-(2R_b+1/3)\Psi>0$、次の希薄化で $X=\Psi/3<0$ となる。従って同じ初期振幅の単純な模型で圧縮の絶対振幅が増える。$R_b=1/2$ なら絶対値の比は $4$ である。実際のピーク比は投影と放射の重力駆動、拡散で変わり、この比をそのまま観測の $C_\ell$ の比へ移せない。

連続式で $\theta=-3\Theta_0'$ なので、速度は $\sin(kc_s\Delta\eta)$ に比例する。単極子の音響極値で速度は零となる一方、温度の中点を通る時刻には速度が大きい。速度 Doppler 源は角度射影でも $j_\ell'$ を持つため、単極子成分とは異なる位相と投影を持つ。この二つを一つの振幅として足す前に、符号と視線角を保持する必要がある。偏光の四重極は深いタイトカップリングで速度から作られるため、温度の圧縮とは概ねずれた位相に応答する。

### 重力の時間変化を積分する応答

$\gamma=R_b'/(1+R_b)$ とし、ポテンシャルを含む音響式の源を $F=\Phi''+\gamma\Phi'-k^2\Psi/3$ とする。$u=\sqrt{1+R_b}\Theta_0$ なら $u''+\omega^2u=\sqrt{1+R_b}\,F$、$\omega^2=k^2c_s^2-\gamma'/2-\gamma^2/4$ である。同次解を $u_1,u_2$ とし、$W=u_1u_2'-u_1'u_2$ と置く。方程式に一階微分がないため $W'=0$ である。点源時刻 $u$ からの Green 関数は
$$
G(\eta,u)=\frac{u_1(u)u_2(\eta)-u_2(u)u_1(\eta)}{W}
\quad (\eta>u)
$$
であり、$\eta=u$ で値が零、右微分が $1$ である。従って初期値を運ぶ同次解へ、この関数で重力源を時間積分して足す。重力が音響周期の間に変化すると、積分が余弦モードと正弦モードの両方を作る。重力駆動の効果を単なる一定の振幅倍率で表せるのは、駆動の時期と初期条件が限られた場合である。

:::example[圧力波の位相を途中で押す]
定数 $\omega=kc_s$、$R_b$ の系で、温度振動のある時刻 $\eta_d$ に有効源 $F=f_0\delta(\eta-\eta_d)$ を加える。この項が作る温度は $(f_0/\omega)\sin[\omega(\eta-\eta_d)]$ であり、温度は連続、時間微分が $f_0$ 増える。押す時刻を半周期ずらすと、後の同じ時刻での寄与の符号が反対になる。従って等密度や再結合のタイミングを変えると、変化の強さだけでなく位相も変わり得る。
:::

## 拡散減衰（Silk 減衰）

散乱間隔の有限性により四重極と速度差が残る。四重極は粘性、速度差は熱伝導の効果を持つ。粒子の拡散距離の二乗が自由行程と経過時間の積に増えるため、短波長が先に消える。係数を得るには光子の偏光込み四重極とバリオン Euler 式の速度差を $\dot\kappa_T^{-1}$ の一次まで解き、音響モードの複素周波数の減衰率を積分する。標準の一次結果は

$$
k_D^{-2}=\int_0^\eta\frac{\dd\eta'}{6\dot\kappa_T(1+R_b)}\left[\frac{R_b^2}{1+R_b}+\frac{16}{15}\right],\qquad A(k)\propto e^{-k^2/k_D^2} . \label{eq-cosmo-silk}
$$

$16/15$ の由来を下の小節で計算する。振幅とパワーの指数には二倍の違いがある。最終散乱の有限幅による位相平均も別に細かい構造を弱くする。

:::example[拡散波数での抑制]
$k=k_D$ では振幅が $e^{-1}$、パワーは $e^{-2}$ 倍となる。
:::


:::example[散乱率の効果]
同じ時間区間で $\dot\kappa_T$ を一様に四倍にすると拡散長 $k_D^{-1}$ は半分になる。
:::

### 熱伝導と粘性の係数

散乱率 $\dot\kappa_T$ の逆数を小さい量とし、各振動周期で $R_b$、背景、散乱率を固定する。:eqref[eq-cosmo-thomson-polarization] を角度射影すると、速い散乱で偏光の係数は $P_0^{\rm pol}=\Pi/2$、$P_2^{\rm pol}=\Pi/10$ に近い。従って $\Pi=5\Theta_2/2$、温度四重極の消散は $\dot\kappa_T(\Theta_2-\Pi/10)=(3/4)\dot\kappa_T\Theta_2$ である。輸送による四重極の源 $2k\Theta_1/5$ と釣り合わせ、$\sigma_\gamma=2\Theta_2$、$\theta_\gamma=3k\Theta_1$ を使うと
$$
\sigma_\gamma\simeq\frac{16\theta_\gamma}{45\dot\kappa_T}.
$$

$s=\theta_\gamma-\theta_b$ と置き、重力と膨張の変化を一周期で無視する。差の式では、光子の圧力 $k^2\Theta_0$ に $\dot\kappa_T(1+1/R_b)s$ が釣り合い、$s\simeq R_b k^2\Theta_0/[(1+R_b)\dot\kappa_T]$ となる。慣性で重み付けした速度 $V=(\theta_\gamma+R_b\theta_b)/(1+R_b)$ を取ると $V'=k^2(\Theta_0-\sigma_\gamma)/(1+R_b)$、$\Theta_0'=-[V+R_b s/(1+R_b)]/3$ である。$V$ を消去し、粘性項に $\theta_\gamma=-3\Theta_0'$ を使えば
$$
\Theta_0''+k^2c_s^2\Theta_0+\gamma_{\rm d}\Theta_0'\simeq0,\qquad
\gamma_{\rm d}=\frac{k^2}{3(1+R_b)\dot\kappa_T}
\left[\frac{R_b^2}{1+R_b}+\frac{16}{15}\right].
$$
一周期当たりの減衰が小さいなら振幅は $\exp[-\int\gamma_{\rm d}\dd\eta/2]$ 倍になる。従って :eqref[eq-cosmo-silk] の二項がそれぞれ速度差と四重極から出る。

### 拡散振幅とパワーへの違う効果

狭い時間幅の最後の散乱と小さい角度を仮定する。散乱面までの共動距離を $\chi_*=\eta_0-\eta_*$ とすると、観測の多重極と波数は $\ell\simeq k\chi_*$ の付近で強く結び付く。:eqref[eq-cosmo-silk] の振幅は $\exp[-k^2/k_D^2]$ なので、パワーの減衰は $\exp[-2k^2/k_D^2]$ である。投影を狭くみなすと $\ell_D=k_D\chi_*$ で指数は $-2(\ell/\ell_D)^2$ になる。実際のスペクトルは視線積分を含み、この変換は包絡線の目安である。

:::example[半分になる角度多重極]
パワーを無減衰の半分にするのは $\ell/\ell_D=\sqrt{\log2/2}\simeq0.589$ である。$\ell=\ell_D$ では振幅が $e^{-1}$、パワーは $e^{-2}\simeq0.135$ である。スペクトルに振幅の式を直接当てはめると減衰を過小評価する。
:::

### 有限幅による位相平均

散乱時刻の分布を中心 $\eta_*$、標準偏差 $\sigma_\eta$ の正規化 Gaussian で近似する。局所の音響波が $\cos(kr_s(\eta))$ で、幅の中で $r_s(\eta)\simeq r_s(\eta_*)+c_s(\eta-\eta_*)$ とする。射影の変化を無視すれば畳み込みは
$$
\int\frac{\dd u}{\sqrt{2\pi}\sigma_\eta}e^{-u^2/(2\sigma_\eta^2)}
\cos(kr_s(\eta_*)+kc_su)
=\cos(kr_s(\eta_*))e^{-(kc_s\sigma_\eta)^2/2}
$$
となる。有限幅は異なる放射時刻の位相を平均し、Silk 減衰は各時刻までの拡散が音源自体を変える。視認性が広ければ球 Bessel 関数による射影も時刻とともに変わるため、実際の温度を求めるには視線積分全体を使う。

### 拡散の有効量と光学的厚さ

再結合より十分早い区間で水素が電離し、電子数密度が $n_e\propto a^{-3}$ と仮定する。$\dot\kappa_T\propto a^{-2}$ であり、放射優勢では $a\propto\eta$ である。$R_b\ll1$ で :eqref[eq-cosmo-silk] を使えば $k_D^{-2}=(8/45)\int\dd\eta/\dot\kappa_T\propto\eta^3$ となる。従ってこの理想化では拡散長が $\eta^{3/2}$、音響地平線が $\eta$ に比例する。再結合で自由電子が急減すると、同じ一定電離度のべき則では拡散を追えない。

光学的厚さ $\tau(\eta)=\int_\eta^{\eta_0}\dot\kappa_T\dd\eta'$ とすると $g(\eta)=(e^{-\tau})'=\dot\kappa_T e^{-\tau}$ であり、$\int_{\eta_i}^{\eta_0}g\,\dd\eta=1-e^{-\tau(\eta_i)}$ である。これは自由行程の累積で最後の散乱の時刻を選ぶ分布である。一方、$k_D^{-2}$ は最後の散乱までの多数の散乱の結果を各時刻で積算する。$g$ が大きい時期と $k_D^{-2}$ の変化が大きい時期は近い場合があるが、二量は確率と長さの二乗という異なる次元を持つ。

### 偏光を省いた場合との定量比較

温度の四重極の衝突項で偏光を省くと、利得は温度四重極だけである。従ってタイトカップリングで $(9/10)\dot\kappa_T\Theta_2=2k\Theta_1/5$ となり、$\sigma_\gamma=8\theta_\gamma/(27\dot\kappa_T)$ である。$\theta_\gamma=-3\Theta_0'$ を用いて音響式へ代入すると、:eqref[eq-cosmo-silk] の角括弧の $16/15$ が $8/9$ に替わる。偏光込みの粘性項に対する比は $(8/9)/(16/15)=5/6$ である。この比較は熱伝導項の変化を無視した一次の範囲に限る。

:::example[同じ短波長のパワーへの差]
$R_b=0$ の理想化で、偏光込みの $k=k_D$ を取る。偏光を外すと振幅の指数は $-5/6$ となり、真の一次近似の $-1$ に比べ振幅は $e^{1/6}\simeq1.181$ 倍大きい。パワーでは $e^{1/3}\simeq1.396$ 倍になる。散乱確率のわずかな係数差でも、指数減衰であるため深い短波長側では差が大きくなる。
:::

## 角度パワースペクトル

全天の温度を $\Theta_o(\hat{\vb n})=\sum a_{\ell m}Y_{\ell m}(\hat{\vb n})$ と展開する。統計的等方性により二点関数は $\langle a_{\ell m}a^*_{\ell'm'}\rangle=C_\ell\delta_{\ell\ell'}\delta_{mm'}$ である。各 Fourier モードの視線積分を Rayleigh 展開で球面へ射影し、応答を $\Delta_\ell(k)$ とする。角度の直交性と :eqref[eq-cosmo-primordial-scalar] のパワーの規約で

$$
C_\ell=4\pi\int_0^\infty\mathcal P_\zeta(k)\left|\Delta_\ell(k)\right|^2\frac{\dd k}{k}. \label{eq-cosmo-cmb-cl}
$$

ガウス等方の全天で推定量 $\widehat C_\ell=\sum_m|a_{\ell m}|^2/(2\ell+1)$ は実自由度 $2\ell+1$ 個の二乗の平均である。実 Gaussian の二乗の分散が平均二乗の二倍なので $\operatorname{Var}\widehat C_\ell=2C_\ell^2/(2\ell+1)$ となる。この不確かさを :term[宇宙分散]{en="cosmic variance" yomi="うちゅうぶんさん"} と呼ぶ。観測雑音がなくても空の実現が一つのために残る。

:::example[低多重極]
$\ell=2$ で相対標準偏差は $\sqrt{2/5}\simeq0.632$ である。
:::


:::example[角度と音響]
平坦空間で $\ell\sim k\chi_*$ と見積もり、$\chi_*/r_s=96$ なら音響の間隔尺度は約 $302$ である。
:::

### 視線式の角度射影

$\mathcal T=\Theta+\Psi$ と置く。:eqref[eq-cosmo-photon-boltzmann] と :eqref[eq-cosmo-thomson-polarization] により、固定 $\mu$ で
$$
\mathcal T'+(ik\mu+\dot\kappa_T)\mathcal T=
\Phi'+\Psi'+\dot\kappa_T
[\Theta_0+\Psi-i\mu\,\theta_b/k-\mathsf P_2(\mu)\Pi/2]
$$
となる。初期の光学的厚さが大きいと、積分因子により境界項は消える。各時刻の右辺に $e^{-\tau(\eta)} e^{-ik\mu(\eta_0-\eta)}$ を掛けて積分したものが観測者の $\mathcal T$ である。$x=k(\eta_0-\eta)$ とすると $j_\ell(x)=(i^\ell/2)\int_{-1}^1e^{-ix\mu}\mathsf P_\ell(\mu)\dd\mu$ である。Doppler の $-i\mu$ は $x$ 微分、$\mathsf P_2(\mu)=(3\mu^2-1)/2$ は $-(3\partial_x^2+1)/2$ に置き換えられる。従って一次の温度多重極は
$$
\Theta_{\ell,o}=\int\dd\eta\,\left[
(g(\Theta_0+\Psi)+e^{-\tau}(\Phi'+\Psi'))j_\ell(x)
+g\frac{\theta_b}{k}j_\ell'(x)
+\frac{g\Pi}{4}(3j_\ell''(x)+j_\ell(x))
\right] .
$$
$g=\dot\kappa_T e^{-\tau}$ は :eqref[eq-cosmo-visibility] である。観測地点だけに依存する重力単極子と観測者の運動双極子は除く。各初期曲率の単位振幅でこれを積分すれば :eqref[eq-cosmo-cmb-cl] に入れる $\Delta_\ell$ となる。

:::example[大角度の極限]
物質期の狭い散乱面で Doppler、四重極、時間変化を無視すると $\Delta_\ell=-j_\ell(k\chi_*)/5$ となる。この場合、同じ $k$ は一般に複数の $\ell$ に寄与する。
:::

### 純粋な Sachs–Wolfe 台地

曲率パワー $\mathcal P_\zeta=A_s$ が一定で、瞬間的再結合の温度源が通常の Sachs–Wolfe 項 $-\zeta/5$ だけとする。:eqref[eq-cosmo-cmb-cl] は $C_\ell=(4\pi A_s/25)I_\ell$、$I_\ell=\int_0^\infty j_\ell(x)^2\dd x/x$ となる。$\ell\ge1$ なら $j_\ell$ は両端で零であり、$j_\ell'=j_{\ell-1}-(\ell+1)j_\ell/x$ を $j_\ell$ に掛けて積分すると
$$
(\ell+1)I_\ell=\int_0^\infty j_\ell j_{\ell-1}\dd x .
$$
$\ell\ge2$ では $j_{\ell-1}'=(\ell-1)j_{\ell-1}/x-j_\ell$ を同様に使い、右辺は $(\ell-1)I_{\ell-1}$ と分かる。$\ell=1$ では $j_0'=-j_1$、$j_0(0)=1$ なので $\int j_1j_0\dd x=1/2$、よって $I_1=1/4$ である。帰納的に $I_\ell=1/[2\ell(\ell+1)]$ となる。従って $\ell(\ell+1)C_\ell/(2\pi)=A_s/25$ である。

:::example[四重極の規格化]
$A_s=2.1\times10^{-9}$ なら台地は $8.4\times10^{-11}$ である。$C_2$ はこれに $2\pi/6$ を掛けた約 $8.80\times10^{-11}$ であり、台地の値そのものとは異なる。
:::

### 角相関と観測ビーム

$\xi(\gamma)=\langle\Theta_o(\hat{\vb s})\Theta_o(\hat{\vb s}')\rangle$ は $\cos\gamma=\hat{\vb s}\cdot\hat{\vb s}'$ だけに依存する。加法定理と $C_\ell$ の定義、Legendre 多項式の直交性から
$$
\xi(\gamma)=\sum_{\ell=0}^\infty\frac{2\ell+1}{4\pi}C_\ell P_\ell(\cos\gamma),
\qquad
C_\ell=2\pi\int_{-1}^1\xi(\arccos u)P_\ell(u)\dd u
$$
である。同一点の分散は $\xi(0)=\sum_\ell(2\ell+1)C_\ell/(4\pi)$ であり、$\ell\gg1$ では $\mathcal D_\ell=\ell(\ell+1)C_\ell/(2\pi)$ が対数多重極幅当たりの寄与に近い。観測では単極子と運動双極子を除く。

軸対称ビーム $B(\hat{\vb s}\cdot\hat{\vb s}')$ を単位積分に規格化する。真の地図をこれで畳み込めば、加法定理から各係数は $B_\ell a_{\ell m}$、$B_\ell=2\pi\int_{-1}^1B(u)P_\ell(u)\dd u$ となる。等方な雑音 $N_\ell$ が真の地図と独立なら、測る平均は $B_\ell^2C_\ell+N_\ell$ である。独立な二地図間の相互相関では雑音の平均が消えるが、推定量の分散には双方の雑音が残る。

### Fourier 空間の共分散から角度パワーを作る

温度の単位原始曲率への応答を、位置空間での Fourier 積分と波数軸に対する Legendre 展開で組み合わせる。進行方向を $\hat{\vb n}$ を用いて表示するなら、加法定理から
$$
a_{\ell m}=4\pi(-i)^\ell\int\frac{\dd^3k}{(2\pi)^3}
\zeta_{\vb k}\Delta_\ell(k)Y_{\ell m}^*(\hat{\vb k}).
$$
逆向きの空の座標では全係数に $(-1)^\ell$ が掛かるが、同じ向きで定義する自己スペクトルは変わらない。二つの Fourier 積分を合わせ、$\langle\zeta_{\vb k}\zeta_{\vb k'}^*\rangle=(2\pi)^3\delta^3(\vb k-\vb k')P_\zeta(k)$ を使う。角度積分で $\delta_{\ell\ell'}\delta_{mm'}$、半径積分の係数は $(4\pi)^2/(2\pi)^3=2/\pi$ である。従って $C_\ell=(2/\pi)\int k^2P_\zeta|\Delta_\ell|^2\dd k$ となり、$P_\zeta=2\pi^2\mathcal P_\zeta/k^3$ から :eqref[eq-cosmo-cmb-cl] へ戻る。多重極を作るには中心の波数だけを代入するのではなく、全波数からの寄与を重ねる。

### Doppler 源を球 Bessel の微分から移す

前の視線式で $x=k(\eta_0-\eta)$ とおく。$j_\ell'(x)=-k^{-1}\dd j_\ell(x)/\dd\eta$ なので $g(\theta_b/k)j_\ell'$ の積分は、端点の項を除き $(g\theta_b)'j_\ell/k^2$ になる。同様に $j_\ell''(x)=k^{-2}\dd^2j_\ell(x)/\dd\eta^2$ であり、四重極項は $(g\Pi/4)j_\ell+3(g\Pi)''j_\ell/(4k^2)$ へ移せる。この入れ替えでは $g$ とその微分が端に残る場合や時刻を有限に切る場合は境界項を保持する。端点が消える理想化では全ての物理を時間源 $S_T$ にまとめて $j_\ell$ だけと積分できる。再結合の急な時間変化を数値微分する場合には、元の $j_\ell',j_\ell''$ を使う表示と照合できる。

### Sachs–Wolfe 台地を積分で確かめる

球 Bessel の積分 $I_\ell$ を直接数値積分する場合、$x=0$ の $\sin x$ と $x\cos x$ の差で桁を失う。小さい $x$ では $j_\ell=x^\ell/(2\ell+1)!!$ に級数を掛ける。球 Bessel 方程式に $x^\ell\sum_n a_nx^{2n}$ を入れると $a_n/a_{n-1}=-1/[2n(2\ell+2n+1)]$ である。大きい側では $\sin x,\cos x$ と漸化式を使う。振動の尾は $j_\ell^2/x=\order{x^{-3}}$ であるため、上端 $X$ の積分の打ち切り誤差は $\order{X^{-2}}$ である。刻みを小さくする検査と $X$ を大きくする検査は別に行う。

~~~python
from math import cos, prod, sin

def spherical_bessel(ell, x):
    if x <= 1:
        term = x ** ell / prod(range(1, 2 * ell + 2, 2))
        value = term
        for index in range(1, 15):
            term *= -x * x / (2 * index * (2 * ell + 2 * index + 1))
            value += term
        return value
    previous = sin(x) / x
    if ell == 0:
        return previous
    current = sin(x) / x ** 2 - cos(x) / x
    for order in range(1, ell):
        previous, current = current, (2 * order + 1) * current / x - previous
    return current

def plateau_integral(ell, upper, intervals):
    step = upper / intervals

    def integrand(x):
        return spherical_bessel(ell, x) ** 2 / x if x else 0.0

    odd = sum(integrand(index * step) for index in range(1, intervals, 2))
    even = sum(integrand(index * step) for index in range(2, intervals, 2))
    return step * (4 * odd + 2 * even + integrand(upper)) / 3

for ell in (2, 3, 4):
    analytical = 1 / (2 * ell * (ell + 1))
    first = plateau_integral(ell, upper=80, intervals=16000)
    extended = plateau_integral(ell, upper=160, intervals=32000)
    print(ell, f"{analytical:.8f}", f"{first:.8f}", f"{extended:.8f}")
~~~

この積分では振動の尾が正であるため、上端を広げた値が解析解に下から近づく。例えば $\ell=2$ では上端 $80,160$ に対してそれぞれ $0.08329418,0.08332359$ となり、解析値 $1/12$ との差が約四分の一になる。高い $\ell$ へそのまま上向きの漸化式を適用すると、$x$ が $\ell$ より小さい部分では丸め誤差が増え得る。ここでは低い三つの多重極の解析式を検算し、一般の高速な球 Bessel 評価の代わりには使わない。

### 広い視線源と波数の近似対応

ある三次元の統計的に一様な場 $d$ を $\int W(\chi)d(\chi\hat{\vb s})\dd\chi$ として射影する。$d$ 自体に時刻依存があれば、それが波数によらない範囲で $W$ に含める。球面展開で $C_\ell=(2/\pi)\int k^2P_d(k)[\int W(\chi)j_\ell(k\chi)\dd\chi]^2\dd k$ となる。もし $P_d$ が一定なら、球 Bessel 関数の完全性
$$
\frac2\pi\int_0^\infty k^2 j_\ell(k\chi)j_\ell(k\chi')\dd k
=\frac{\delta_{\rm D}(\chi-\chi')}{\chi^2}
$$
により $C_\ell=P_d\int W^2/\chi^2\dd\chi$ である。この恒等式は平面波表示の三次元デルタ関数に球面調和の加法定理を適用し、$\delta^3(\vb x-\vb x')=\chi^{-2}\delta(\chi-\chi')\delta^2(\hat{\vb s},\hat{\vb s}')$ と係数を比較すれば得る。

多重極が大きく、$W$ が広く滑らかな場合は、$P_d(k)$ が主に $k\chi\simeq\ell+1/2$ の付近で寄与する。そこで $P_d$ をその波数で評価して先に近似的に外へ出すと
$$
C_\ell\simeq\int\frac{W^2(\chi)}{\chi^2}
P_d((\ell+1/2)/\chi)\dd\chi
$$
となる。これを :term[Limber 近似]{en="Limber approximation" yomi="りんばーきんじ"} と呼ぶ。小さい多重極、急な $P_d$、鋭い時間境界を持つ $W$ では波数を一本に絞る近似が悪くなる。瞬間的な散乱面の Sachs–Wolfe 台地で $W^2$ にデルタ関数の二乗を入れることはできず、先ほどの $j_\ell^2$ の積分を使う。広いレンズや銀河の射影と、狭い一次温度の源とで使える近似が違う。

## 偏光と宇宙論パラメータの決定

:term[線偏光]{en="linear polarization" yomi="せんへんこう"} は Thomson 散乱に入る放射の四重極が作る。等方入射では二つの横偏光へ同量が散るため消える。局所基底の回転で Stokes の組 $Q\pm iU$ はスピン $\pm2$ の場となり、全天で$E$ と $B$ に分ける。空間反転で $E_{\ell m}$ は $(-1)^\ell$、$B_{\ell m}$ は $(-1)^{\ell+1}$ を受ける。線形スカラーの一つの Fourier モードは波数と視線を含む面で反射対称なので、向きを反転する $B$ を作れない。テンソルは原始 $B$ を作り得て、重力レンズも異なる方向の $E$ を混ぜて $B$ を生む。再電離で遅く散乱した成分は大角度 $E$ に現れる。

音響角から $r_s/\chi_*$、圧縮と希薄化の比からバリオン慣性、等密度付近のポテンシャル減衰から物質量、振幅と再電離情報から原始振幅を推定する。温度・偏光・レンズを合わせて相関の縮退を減らす。周波数依存する銀河前景と機器応答の同時推定は :cite[dodelson-cosmology] の枠組みを用い、実データではマスクと雑音も扱う。

:::example[再散乱]
光学的厚さ $\tau_{\rm re}=0.05$ の一次温度成分は無散乱率 $e^{-0.05}$ 倍で、パワーは $e^{-0.1}\simeq0.905$ 倍となる。
:::


:::example[角度縮退]
$r_s$ と $\chi_*$ を同じ割合で変えればピーク間隔の尺度は同じであり、間隔だけでは両者を別々に決められない。
:::

### 角尺度が固定するパラメータの組合せ

空間が平坦なら再結合時刻までの共動距離は $\chi_*=\int_0^{z_*}c\,\dd z/H(z)$ である。音響角は $\theta_s=r_s(\eta_*)/\chi_*$、理想的な音響間隔は $\ell_A=\pi/\theta_s$ である。同じ $\theta_s$ を保つ微小変化は $\delta\log r_s=\delta\log\chi_*$ を満たす。例えば初期の膨張を速めると、固定した音速と再結合赤方偏移の近似では $\int c_s\dd a/(a^2H)$ が減る。低赤方偏移の膨張も変えて $\chi_*$ を同じ割合だけ短くすれば、ピーク間隔を保てる。この比較を実際のデータへ使うには、初期成分が再結合赤方偏移、重力駆動、拡散を同時に変えることも入れる。

ピーク間隔を最初のピーク位置と同一視すると、初期位相と投影を落とす。薄い散乱面、一定の初期位相、単極子だけの近似なら $kr_s=n\pi$ と $\ell\simeq k\chi_*$ を重ねて $\ell_n\simeq n\ell_A$ と見積もれる。しかし放射の重力駆動は有効温度を変え、自由流の応力や有限の散乱幅、Doppler との干渉も加わる。従って $\ell_A$ は振動の代表的な間隔であり、最初のピークの $\ell$ をそのまま $\pi/\theta_s$ に代入しても同じ角尺度にはならない。

:::example[角尺度と拡散尺度が異なる感度を持つ]
初期のある区間の $\dd\eta$ を一律に $0.9$ 倍し、$c_s$ と共形散乱率を固定するという近似を考える。音響距離は $0.9$ 倍、拡散距離は $\sqrt{0.9}$ 倍である。$\chi_*$ も $0.9$ 倍すれば音響角は元のままだが、拡散角は $1/\sqrt{0.9}\simeq1.054$ 倍になる。実際に膨張率を変えると電子の電離史も変わるため、この算例では作用する二種類の積分の差だけを取り出している。
:::

### 複数の角度スペクトルの共分散

全空、雑音なしで温度 $T$ と偏光 $E$ が同時に Gaussian かつ等方とする。実球面調和基底の各係数は、固定した $\ell$ について共分散行列の成分が $C_\ell^{XY}$ の二変量正規分布に従う。同じ $\ell$ の二つの推定量の共分散は、独立な $m$ の項だけを数え、Wick 展開の三組のうち積の平均で相殺する組を除くと
$$
\operatorname{Cov}(\hat C_\ell^{XY},\hat C_\ell^{ZW})
=\frac{C_\ell^{XZ}C_\ell^{YW}+C_\ell^{XW}C_\ell^{YZ}}{2\ell+1}
$$
である。$TE$ の分散は $(C_\ell^{TT}C_\ell^{EE}+(C_\ell^{TE})^2)/(2\ell+1)$ であり、相関が零でも積の推定量にはゆらぎが残る。また $\hat C_\ell^{TT},\hat C_\ell^{EE}$ の共分散は $2(C_\ell^{TE})^2/(2\ell+1)$ である。これらを独立とみなすと、同じ空を測った情報を重複して数える。観測領域を切ると異なる $\ell$ も混ざる。

### スカラー偏光の視線積分

偏光基底を波数を極とする $(\vb e_\theta,\vb e_\varphi)$ に取る。軸対称な線形スカラーで $U=0$、$Q=P$ と規格化する。:eqref[eq-cosmo-thomson-polarization] より光路の散乱源は $3(1-\mu^2)\Pi/4$ である。背景偏光が零なら重力は一次では偏光を作らないので、輸送式を積分すると
$$
Q_o(k,\mu)=\frac34(1-\mu^2)\int_0^{\eta_0}g(\eta)\Pi(k,\eta)e^{-ik\mu(\eta_0-\eta)}\dd\eta
$$
となる。スピン重み付き場を $Q\pm iU=-\sum_{\ell m}(E_{\ell m}\pm iB_{\ell m})\,{}_{\pm2}Y_{\ell m}$ と規約する。$m=0$ に対し $\vphantom{Y}{}_2Y_{\ell0}={}_{-2}Y_{\ell0}=[(2\ell+1)(\ell-2)!/(4\pi(\ell+2)!)]^{1/2}P_\ell^2(\mu)$ である。従って二つの展開係数が等しく、$B$ は零である。角度射影に必要な積分は
$$
\frac12\int_{-1}^1e^{-ix\mu}(1-\mu^2)P_\ell^2(\mu)\dd\mu
=-(-i)^\ell\frac{(\ell+2)!}{(\ell-2)!}\frac{j_\ell(x)}{x^2}.
$$
実際、$P_\ell^2=(1-\mu^2)P_\ell''$ とすれば $(1-\mu^2)^2$ とその導関数が端点で零なので、二回部分積分できる。残った $\mu$ のべきを $i\partial_x$ に置き換え、球 Bessel 方程式 $j_\ell''+2j_\ell'/x+[1-\ell(\ell+1)/x^2]j_\ell=0$ を使う。温度と同じ位相規約で
$$
\Delta_\ell^E(k)=\sqrt{\frac{(\ell+2)!}{(\ell-2)!}}\int_0^{\eta_0}
\frac34g(\eta)\Pi(k,\eta)\frac{j_\ell(k(\eta_0-\eta))}{[k(\eta_0-\eta)]^2}\dd\eta
$$
を得る。$\Pi$ は単位原始曲率への応答である。:eqref[eq-cosmo-cmb-cl] の積分で $EE$ は $|\Delta_\ell^E|^2$、$TE$ は $\Delta_\ell^T\Delta_\ell^{E*}$ を使う。$Q$ または $E$ の符号を替える場合、比較時には $TE$ の符号もそろえる。

### レンズによる $B$ 偏光の生成

平坦背景で最終散乱面までの共動距離を $\chi_*$ とする。無摂動の光路上で微小偏向を足す Born 近似で、レンズポテンシャルは
$$
\varphi_{\rm L}(\hat{\vb s})=-\int_0^{\chi_*}\frac{\chi_*-\chi}{\chi_*\chi}
(\Phi+\Psi)(\eta_0-\chi,\chi\hat{\vb s})\dd\chi
$$
である。観測方向のずれは $\nabla_{\hat s}\varphi_{\rm L}$ となる。係数 $(\chi_*-\chi)/\chi_*$ は距離 $\chi$ の偏向が源位置へ及ぼす割合、$1/\chi$ は横座標微分を角度微分へ変える係数である。

小角近似で $P=Q+iU$ とすると、観測偏光は一次の偏向で $P_{\rm obs}=P_{\rm src}+\nabla\varphi_{\rm L}\cdot\nabla P_{\rm src}$ となる。補正の Fourier 積分では $-[\vb L\cdot(\boldsymbol\ell-\vb L)]\varphi_{\rm L}(\vb L)P_{\rm src}(\boldsymbol\ell-\vb L)$ が入る。$E+iB$ へ直す角度位相を掛けると、積分は $e^{2i(\varphi_{\boldsymbol\ell-\vb L}-\varphi_{\boldsymbol\ell})}$ を含む。その正弦は一般には零にならないため、源が $E$ だけでもレンズで $B$ が出る。偏向が空の平面全体で一様なら地図は並進するだけで、この変換は起きない。

### 再電離偏光の時間と角度

観測時刻の上限を $\eta_0$、再電離に対応する区間の光学的厚さを $\tau_{\rm re}$ とする。再電離を通り抜ける一次光子の振幅は $e^{-\tau_{\rm re}}$、新たに散乱する確率は $1-e^{-\tau_{\rm re}}$ である。再結合後にほぼ自由に進んだ温度の単極子 $A$ は、再電離までの間に $\Theta_2\simeq A j_2(k(\eta_{\rm re}-\eta_*))$ の四重極を作る。$k(\eta_{\rm re}-\eta_*)\ll1$ では $j_2$ が二次なので偏光の源も小さい。一方、早く振動する四重極は角度・時間で平均されやすい。これが再電離が追加する $E$ が主に大角度に見える理由である。

低い多重極の新しい $E$ の振幅は薄い再電離の近似で $\tau_{\rm re}$ に比例し、その自己パワーは $A_s\tau_{\rm re}^2$ の尺度を持つ。一次 CMB の小角度パワーは主に $A_s e^{-2\tau_{\rm re}}$ を測る。後者を一定とすれば $\dd\log A_s=2\dd\tau_{\rm re}$ である。従って大角度偏光が加わると、二つのパラメータの縮退を減らせる。この説明では各スペクトルの形や前景も重要なので、二つの単一振幅だけで実際の尤度を置き換えているわけではない。

### パリティと全天の尤度

同時に等方な $T,E,B$ の二点を考える。空間反転の対称性があれば、$T_{\ell m}$ と $E_{\ell m}$ は $(-1)^\ell$、$B_{\ell m}$ は $(-1)^{\ell+1}$ を受ける。同じ $\ell$ の $TB,EB$ は全体で符号を変えるため期待値が零になる。一方 $BB$ は偶であり、レンズの二次の寄与を含めてもこの対称性だけでは消えない。測定で非零の $TB$ や $EB$ を得た場合には、宇宙の非対称性だけでなく偏光角の較正も調べる必要がある。実際、全方向の基底が小角 $\beta$ ずれれば、$E+iB$ に $e^{2i\beta}$ が掛かり、元の $B$ が零でも $B\simeq2\beta E$ となる。

温度だけの全空 Gaussian 場なら、$n=2\ell+1$ 個の実係数の結合確率密度は $C_\ell^{-n/2}\exp[-n\hat C_\ell/(2C_\ell)]$ に比例する。従って $-2\log\mathcal L=n(\log C_\ell+\hat C_\ell/C_\ell)+$ 定数である。最大は $C_\ell=\hat C_\ell$ にあり、近傍の二次微分から分散 $2C_\ell^2/n$ が再現される。低多重極では $\hat C_\ell$ は正で歪んだ $\chi^2_n$ 分布を持つため、対称な誤差棒だけで全尤度を表せない。偏光を同時に使うと各 $m$ の共分散行列が必要になる。

### 再電離の深さと電子の数

水素だけが時刻 $z_{\rm re}$ で一斉に電離し、それ以降の電子が保存するとする。現在に換算した数密度を $n_{e0}$ とすると $n_e(z)=n_{e0}(1+z)^3$ である。固有時間に沿う散乱確率から
$$
\tau_{\rm re}=c\sigma_T\int_0^{z_{\rm re}}\frac{n_e(z)}{(1+z)H(z)}\dd z
$$
となる。分母の $(1+z)$ は $\dd t=-\dd z/[(1+z)H]$ から来る。全区間で物質優勢と近似し、$H(z)=H_0\sqrt{\Omega_{m0}}(1+z)^{3/2}$ を入れると
$$
\tau_{\rm re}=\frac{2c\sigma_T n_{e0}}{3H_0\sqrt{\Omega_{m0}}}
\big[(1+z_{\rm re})^{3/2}-1\big]
$$
である。遅い宇宙定数を無視したため低赤方偏移部分は概算である。二重に電離した He や部分電離の期間を入れると同じ $\tau$ でも $z_{\rm re}$ の割り当てが変わる。大角度偏光の形には電子の全積分だけでなく散乱時刻の分布も現れる。

:::example[一斉に電離する時刻の感度]
同じ $n_{e0}$ と背景で $1+z_{\rm re}$ を $8$ から $10$ に変えると、深さの比は $(10^{3/2}-1)/(8^{3/2}-1)\simeq1.416$ である。薄い再電離で他の量を固定した場合、大角度偏光のパワーはその二乗の約 $2.01$ 倍になる。正確な形の計算では四重極が発達した時間も同時に変わる。
:::

### 複数周波数を用いた成分の推定

各周波数を熱力学的な温度単位へ直した地図を $d_i$ とする。CMB は共通振幅 $s$、ある前景は既知の周波数係数 $f_i$ と空間振幅 $a_f$ を持つと理想化し、$d_i=s+f_i a_f+n_i$ とする。二周波数で $f_1\ne f_2$ なら $\hat s=(f_2 d_1-f_1d_2)/(f_2-f_1)$ として前景を消せる。独立な雑音分散 $N_i$ に対し、推定値の分散は $(f_2^2N_1+f_1^2N_2)/(f_2-f_1)^2$ である。従って係数の差が小さいと、平均で前景が相殺しても雑音が増幅する。

一般に $d=A x+n$、$N=\langle nn^{\mathsf T}\rangle$ で正定値とし、$A$ の列を独立とする。$(d-Ax)^{\mathsf T}N^{-1}(d-Ax)$ を微分すると $(A^{\mathsf T}N^{-1}A)\hat x=A^{\mathsf T}N^{-1}d$ である。$\hat x-x$ を雑音で表し、$N$ を挟んで二乗平均を取れば共分散は $(A^{\mathsf T}N^{-1}A)^{-1}$ となる。これは既知の周波数依存と線形混合のモデルである。銀河の熱塵の周波数依存が場所で変わる場合は、この固定した行列だけでは除けない。テンソルや再電離の偏光を推定するときには、周波数間で共通のビームや較正に直したうえで対象の角度尺度を解析する。

### 観測マスクと角度モードの混合

全天の温度に実数の重み $W_{\rm sky}(\hat{\vb s})$ を掛けた地図を考える。マスクした地図を球面調和で展開すれば
$$
\tilde a_{\ell m}=\sum_{\ell' m'}K_{\ell m,\ell'm'}a_{\ell'm'},\qquad
K_{\ell m,\ell'm'}=\int Y_{\ell m}^*W_{\rm sky}Y_{\ell'm'}\dd\Omega .
$$
等方な入力について二点関数を代入すると $\langle\tilde C_\ell\rangle=\sum_{\ell'}M_{\ell\ell'}C_{\ell'}$、$M_{\ell\ell'}=(2\ell+1)^{-1}\sum_{mm'}|K_{\ell m,\ell'm'}|^2$ となる。マスクが一定なら $K$ は対角であるが、場所に依存すれば異なる $\ell$ が混ざる。単に観測面積の比だけで信号を割ると、狭いピークや急な減衰で混合による偏りが残る。有限幅の帯にまとめ、この混合行列や観測されないモードをモデルに含めれば比較できる。

:::example[双極子の重みが隣の多重極を混ぜる]
$W_{\rm sky}=1+\varepsilon\cos\theta$ とする。$\cos\theta$ は $\ell=1,m=0$ の球面調和に比例するので $m$ を保存し、$m=0$ では Legendre 多項式の漸化式により $\ell$ を $\ell\pm1$ に移す。例えば入力が $Y_{20}$ だけなら、$\cos\theta P_2=3P_3/5+2P_1/5$ により $\ell=1,3$ の偽の振幅が生じる。パワーの混入は $\varepsilon^2$ の次数であり、入力の独立なモードの間にも共分散が生じる。
:::

### レンズを温度の共分散から推定する

小角近似で観測温度を $T(\vb s+\nabla\varphi_{\rm L})$ とする。固定したレンズのもとで、一次の温度補正は
$$
\delta T_{\boldsymbol\ell}
=-\int\frac{\dd^2L}{(2\pi)^2}\,
 [\vb L\cdot(\boldsymbol\ell-\vb L)]
 \varphi_{{\rm L},\vb L}T_{\boldsymbol\ell-\vb L}
$$
である。ここで温度の無レンズの二点は $(2\pi)^2\delta^2(\boldsymbol\ell+\boldsymbol\ell')C_\ell$ とする。$\boldsymbol\ell+\boldsymbol\ell'=\vb L\ne0$ の固定した組に対して、$\langle\delta T_{\boldsymbol\ell}T_{\boldsymbol\ell'}\rangle+\langle T_{\boldsymbol\ell}\delta T_{\boldsymbol\ell'}\rangle$ は
$$
[\vb L\cdot\boldsymbol\ell\, C_\ell+
 \vb L\cdot\boldsymbol\ell'\,C_{\ell'}]\varphi_{{\rm L},\vb L}
$$
である。一方の補正の中の温度波数が $-\boldsymbol\ell'$ となるので、積分の負号が打ち消される。この式は無レンズの一次温度と固定したレンズポテンシャルを独立に近似した応答であり、近い時代の重力温度との相関を別に扱う。

同じ $\vb L$ になる温度モードの組を応答と雑音で重み付けして足すと、温度の二次量でレンズを推定できる。レンズの平均は零でもそのパワーは非零であり、温度四点の連結部分に現れる。無レンズ Gaussian 場からも推定量の分散があるので、四点を使うときは非連結部分を差し引く。温度ピークの平滑化によるレンズの制約と、四点から得る制約には同じ天空が関わるため、合成時に相関も考慮する。

## 演習問題

::::exercise{#exr-cosmo-cmb-1 level=basic}
黒体の温度を三倍にするとエネルギー密度は何倍か。
:::solution
温度の四乗に比例するので $3^4=81$ 倍である。
:::
::::

::::exercise{#exr-cosmo-cmb-2 level=advanced}
$y$ 歪みで光子数が保存される理由を示せ。
:::solution
$x^2\Delta n=y\partial_x(x^4n')$ を積分する。$x\to0$ で $n'\sim-x^{-2}$、無限遠で指数的減衰なので両端が零となる。
:::
::::

::::exercise{#exr-cosmo-cmb-3 level=basic}
物質期の通常 Sachs–Wolfe で $\zeta=5\times10^{-5}$ の温度比を求めよ。
:::solution
$\Theta=-\zeta/5=-10^{-5}$ となる。
:::
::::

::::exercise{#exr-cosmo-cmb-4 level=advanced}
音響式で$k>0$ に対し重力ポテンシャルとバリオン慣性が時間一定なら、音響振動が中心とする有効温度 $\Theta_0+\Psi$ の値を求めよ。
:::solution
定常の解は $\Theta_0=-(1+R_b)\Psi$ なので $\Theta_0+\Psi=-R_b\Psi$ である。
:::
::::

::::exercise{#exr-cosmo-cmb-5 level=basic}
$k=2k_D$ の拡散によるパワー比を求めよ。
:::solution
振幅が $e^{-4}$ であるため、二乗すると $e^{-8}$ になる。
:::
::::

::::exercise{#exr-cosmo-cmb-6 level=advanced}
全天で $\ell=49$ の相対分散を求めよ。
:::solution
独立実自由度が $2\ell+1=99$ なので $2/99$ である。
:::
::::

::::exercise{#exr-cosmo-cmb-7 level=basic}
$R_b=2$ の音速の二乗を求めよ。
:::solution
式 :eqref[eq-cosmo-acoustic] で $1/[3(1+2)]=1/9$ になる。
:::
::::

::::exercise{#exr-cosmo-cmb-8 level=advanced}
再散乱による抑制だけを仮定し、パワーが $1/2$ となる深さを求めよ。
:::solution
パワーは振幅の無散乱率の二乗 $e^{-2\tau_{\rm re}}$ なので $\tau_{\rm re}=(\log2)/2$ である。
:::
::::

::::exercise{#exr-cosmo-cmb-9 level=advanced}
一定 $R_b=1/2$、静的 $\Psi<0$、初期の有効温度 $X_i=\Psi/3$、初期速度零を仮定する。拡散を外し、最初の圧縮と次の希薄化の温度の絶対振幅比を求めよ。両時点の速度も求めよ。
:::solution
$X=\Theta_0+\Psi=-R_b\Psi+(X_i+R_b\Psi)\cos(kc_s\Delta\eta)$ である。位相 $\pi$ で $X=-(2R_b+1/3)\Psi=-(4/3)\Psi$、位相 $2\pi$ で $X=\Psi/3$ に戻る。比は $4$ である。また $\theta=-3\Theta_0'$ は $\sin(kc_s\Delta\eta)$ に比例するので両時点で零である。定数背景の単極子の比較であり、全天の異なる角度多重極でのパワー比ではない。
:::
::::

::::exercise{#exr-cosmo-cmb-10 level=advanced}
全空で雑音を無視し、ある $\ell$ について $C_\ell^{TT}=4C$、$C_\ell^{EE}=C$、$C_\ell^{TE}=C$、$C>0$ とする。$\hat C_\ell^{TE}$ の分散と $\hat C_\ell^{TT},\hat C_\ell^{EE}$ の共分散を求めよ。三つの観測量に同じ誤差を割り当てられるか。
:::solution
$TE$ の分散は $(4C^2+C^2)/(2\ell+1)=5C^2/(2\ell+1)$ である。$TT$ と $EE$ の共分散は $2C^2/(2\ell+1)$ で、非零である。各自己スペクトルの分散はそれぞれ $32C^2/(2\ell+1)$ と $2C^2/(2\ell+1)$ であるため、三者に同じ分散を割り当てると正しくない。共分散があるので、同じ空の三つの推定量を独立として尤度を掛け合わせることもできない。
:::
::::

::::exercise{#exr-cosmo-cmb-11 level=advanced}
小さい角度 $\beta$ で偏光の較正が回転し、真の $B$ が零である。測定する $EB$ と $BB$ の先頭次数を求め、宇宙のパリティとの関係を述べよ。
:::solution
$E+iB=e^{2i\beta}E_{\rm true}$ により、$E=E_{\rm true}\cos(2\beta)$、$B=E_{\rm true}\sin(2\beta)$ となる。従って $C_\ell^{EB}\simeq2\beta C_{\ell,\rm true}^{EE}$、$C_\ell^{BB}\simeq4\beta^2C_{\ell,\rm true}^{EE}$ である。較正の回転でパリティが奇の相関が作られるため、測定した $EB$ のみでは宇宙にパリティの破れがあると判定できない。
:::
::::

::::exercise{#exr-cosmo-cmb-12 level=basic}
小角度の拡散の近似で、$\ell/\ell_D=1/2$ ならパワーの抑制率はいくらか。振幅と区別せよ。
:::solution
振幅は $e^{-(\ell/\ell_D)^2}=e^{-1/4}$ である。パワーはその二乗の $e^{-1/2}\simeq0.607$ 倍である。$\ell_D$ と観測多重極の換算は薄い散乱面での包絡線の近似なので、同じ比が実測の個々のピークに厳密に当てはまるという意味ではない。
:::
::::

## 確認テスト

::::quiz
静的な重力ポテンシャルに対する積分 Sachs–Wolfe 項はどうなるか。

- [x] 零
- [ ] 必ず正
- [ ] 必ず負

:::explanation
時間微分 $\Phi'+\Psi'$ が零なので光路積分も零である。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
バリオン慣性比 $R_b=1/3$ の音速を求めよ。
:::explanation
$c_s^2=1/[3(1+1/3)]=1/4$ なので $c_s=1/2$ である。
:::
::::

::::quiz{answer="0.367879" tol="1e-5"}
$k=k_D$ の拡散による振幅比を求めよ。
:::explanation
指数は $-1$ なので $e^{-1}\simeq0.367879$ である。
:::
::::

::::quiz
線形スカラーが Thomson 散乱で直接生成する偏光はどれか。

- [x] $E$
- [ ] $B$ のみ
- [ ] 円偏光のみ

:::explanation
反射対称性により奇パリティの $B$ は消え、$E$ が生成される。
:::
::::

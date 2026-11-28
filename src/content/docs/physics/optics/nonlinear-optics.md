---
title: 非線形光学
description: 非線形分極から波長変換と Kerr 効果を導き、分散との釣合いおよび光と物質波の相互作用を計算する。
notation:
  - symbol: '\chi^{(j)}'
    meaning: j 次の非線形感受率
  - symbol: 'A_j'
    meaning: 周波数 j の電場包絡線
  - symbol: '\Delta k'
    meaning: 波数不整合
  - symbol: 'n_{\rm K}'
    meaning: 強度規格化した Kerr 係数
  - symbol: '\gamma'
    meaning: 導波モードの非線形係数
  - symbol: '\beta_2'
    meaning: 群速度分散
---

## この章の目標

- 非線形分極から生成周波数と選択則を求める。
- 位相整合、三波結合と利得を計算する。
- Kerr 効果と分散を組み合わせてソリトンを導く。
- 電場、音波、分子や格子の振動と光の相互作用を比較する。

## 非線形分極と非線形感受率

### 電場に比例しない分極
弱い電場では束縛電荷の復元力を変位に比例させられる。しかし場を大きくすると高次の変位も増え、磁気応答を無視した局所・瞬時の模型では
$$
P_i=\eps_0\left(\sum_j\chi^{(1)}_{ij}E_j+\sum_{jk}\chi^{(2)}_{ijk}E_jE_k+
\sum_{jkl}\chi^{(3)}_{ijkl}E_jE_kE_l+\cdots\right)
\label{eq-optics-nonlinear-polarization}
$$
となる。$\chi^{(m)}$ が $m$ 次の :term[非線形感受率]{en="nonlinear susceptibility" yomi="ひせんけいかんじゅりつ"} であり、SI で単位は $(\mathrm{m/V})^{m-1}$ である。反転対称な媒質では極性ベクトルの関係 $\vb P(-\vb E)=-\vb P(\vb E)$ を Taylor 展開すると、バルクの電気双極子応答で偶数次数が消える。表面では反転後に同じ環境に戻らないので、この議論は成り立たない。

分極の遅れを残す場合は、例えば二次項を $\eps_0\int_0^\infty\int_0^\infty h(u,v)E(t-u)E(t-v)\dd u\dd v$ と書く。入力指数を掛けると出力周波数は入力周波数の符号付きの和になり、係数は $h$ の二変数 Fourier 変換で決まる。このため同じ二次でも、光と低周波電場を混ぜる係数と、二つの光を混ぜる係数は一般に異なる。以下では偏光の射影を含めた有効スカラー係数を使う。

:::example[二乗に含まれる二成分]
$E=E_0\cos\omega t$ なら $E^2=E_0^2(1+\cos2\omega t)/2$。直流が光整流、振動成分が第 2 高調波の源である。
:::
:::example[三次で変わる基本波]
$E^3=E_0^3(3\cos\omega t+\cos3\omega t)/4$ である。三次で第三高調波が現れる一方、基本波と同じ周波数の分極も変わる。
:::

### 非調和振動子の摂動計算
束縛電荷の変位について $\ddot x+\Gamma\dot x+\omega_0^2x+\alpha x^2+\zeta x^3=(q/m)E$ とする。$x=x_1+x_2+\cdots$、$x_j$ は場の $j$ 次の量と置くと、一次は $D(\omega)X_1=(q/m)A$、$D(\omega)=\omega_0^2-\omega^2-i\Gamma\omega$。二次は $D(2\omega)X_2=-\alpha X_1^2/2$ となる。ここで $A$ と $X_j$ は実部表示の複素振幅である。$P_2=Nq X_2=\eps_0\chi^{(2)}A^2/2$ と比べれば
$$
\chi^{(2)}(2\omega;\omega,\omega)=
-\frac{Nq^3\alpha}{\eps_0m^2 D(2\omega) D(\omega)^2}.
$$
従って二倍周波数側が共鳴しても、基本波側が共鳴しても分極は強まり得る。同時に減衰による位相・吸収が増えるため、係数の大きさだけから出射強度を推測できない。$\alpha=0$ の反転対称模型では、この機構の二次項はなく、三次の $\zeta$ は残る。

:::example[二周波から差の出現を確かめる]
実電場が $a\cos\omega_a t+b\cos\omega_b t$ のとき、二乗の交差項は $ab[\cos((\omega_a+\omega_b)t)+\cos((\omega_a-\omega_b)t)]$。例えば $\omega_a/(2\pi)=200\,\mathrm{THz}$、$\omega_b/(2\pi)=190\,\mathrm{THz}$ なら和が $390\,\mathrm{THz}$、差が $10\,\mathrm{THz}$ である。両方の振動が分極に出ても、結晶での吸収や位相整合が同じとは限らない。
:::

### 実電場から和・差周波数の複素分極を取り出す

スカラーで瞬時の二次応答を考え、$E=\Re(E_1e^{-i\omega_1t}+E_2e^{-i\omega_2t})$、$\omega_1>\omega_2>0$ とする。以下では生成周波数が互いに異なると仮定する。同じ周波数に複数の項が重なる場合は、その複素振幅を足す。実部を共役との和に展開すると、$E^2$ の $e^{-2i\omega_1t}$ の係数は $E_1^2/4$、$e^{-i(\omega_1+\omega_2)t}$ の係数は $E_1E_2/2$、$e^{-i(\omega_1-\omega_2)t}$ の係数は $E_1 E_2^*/2$。生成した分極も実部を取る約束で複素振幅に戻すなら、これらを二倍して $\eps_0\chi^{(2)}$ を掛ける。従って第 2 高調波には係数 $1/2$、異なる周波数の和と差には係数 $1$ が付く。差周波数では下の周波数の位相が引かれる。実際に出力光が蓄積するかは、それぞれの分極の空間位相と生成波の伝播位相を比較して決まる。

:::example[上の光だけを入れる]
下の周波数を消して直流と高調波の項を比較する。

入力は上の周波数の実振幅 $E_1=1000000\,\mathrm{V/m}$、下の周波数の実振幅 $E_2=0\,\mathrm{V/m}$、下の周波数の位相 $\phi_2=0\,\mathrm{rad}$、二次感受率 $\chi^{(2)}\simeq 2\times 10^{-12}\,\mathrm{m/V}$ である。

ここでは上の波を $E_1\cos\omega_1 t$、下の波を $E_2\cos(\omega_2t-\phi_2)$ とし、$E_1,E_2$ は実数とする。真空誘電率は $\epsilon_0\simeq8.8542\times10^{-12}\,\mathrm{F/m}$ である。

$$
\begin{aligned}
\epsilon_0\chi^{(2)}&=\epsilon_0\cdot\chi^{(2)},\\
\epsilon_0\chi^{(2)}/(\mathrm{F/V})&\simeq (8.854187819\times10^{-12})(2\times10^{-12})\simeq 1.7708\times 10^{-23}.
\end{aligned}
$$

各周波数とその共役の積を足す。

$$
\begin{aligned}
P_{\rm dc}&=\frac{\epsilon_0\chi^{(2)} \left(E_1^{2} + E_2^{2}\right)}{2},\\
P_{\rm dc}/(\mathrm{C/m^2})&\simeq \frac{\left(1.770837564\times 10^{-23}\right) \cdot \left(\left(1000000\right)^{2} + \left(0\right)^{2}\right)}{2}\simeq 8.8542\times 10^{-12}.
\end{aligned}
$$

同じ周波数の二因子から作る項である。

$$
\begin{aligned}
P_{2\omega_1}&=\frac{\epsilon_0\chi^{(2)} E_1^{2}}{2},\\
P_{2\omega_1}/(\mathrm{C/m^2})&\simeq \frac{\left(1.770837564\times 10^{-23}\right) \cdot \left(1000000\right)^{2}}{2}\simeq 8.8542\times 10^{-12}.
\end{aligned}
$$

位相ゼロの上の波と下の波を掛ける。

$$
\begin{aligned}
P_\Sigma&=\epsilon_0\chi^{(2)} E_1 E_2 e^{i \phi_2},\\
P_\Sigma/(\mathrm{C/m^2})&\simeq \left(1.770837564\times 10^{-23}\right) \cdot \left(1000000\right) \cdot \left(0\right) \cdot e^{i \cdot \left(0\right)}=0.
\end{aligned}
$$

差の項では下の波の共役を掛ける。

$$
\begin{aligned}
P_\Delta&=\epsilon_0\chi^{(2)} E_1 E_2 e^{- i \phi_2},\\
P_\Delta/(\mathrm{C/m^2})&\simeq \left(1.770837564\times 10^{-23}\right) \cdot \left(1000000\right) \cdot \left(0\right) \cdot e^{- i \cdot \left(0\right)}=0.
\end{aligned}
$$

二つの入力を要求する和と差の分極は生じない。
:::

:::example[二つの光を同相で入れる]
両周波数の実振幅を等しくする。

入力は上の周波数の実振幅 $E_1=1000000\,\mathrm{V/m}$、下の周波数の実振幅 $E_2=1000000\,\mathrm{V/m}$、下の周波数の位相 $\phi_2=0\,\mathrm{rad}$、二次感受率 $\chi^{(2)}\simeq 2\times 10^{-12}\,\mathrm{m/V}$ である。

ここでは上の波を $E_1\cos\omega_1 t$、下の波を $E_2\cos(\omega_2t-\phi_2)$ とし、$E_1,E_2$ は実数とする。真空誘電率は $\epsilon_0\simeq8.8542\times10^{-12}\,\mathrm{F/m}$ である。

$$
\begin{aligned}
\epsilon_0\chi^{(2)}&=\epsilon_0\cdot\chi^{(2)},\\
\epsilon_0\chi^{(2)}/(\mathrm{F/V})&\simeq (8.854187819\times10^{-12})(2\times10^{-12})\simeq 1.7708\times 10^{-23}.
\end{aligned}
$$

各周波数とその共役の積を足す。

$$
\begin{aligned}
P_{\rm dc}&=\frac{\epsilon_0\chi^{(2)} \left(E_1^{2} + E_2^{2}\right)}{2},\\
P_{\rm dc}/(\mathrm{C/m^2})&\simeq \frac{\left(1.770837564\times 10^{-23}\right) \cdot \left(\left(1000000\right)^{2} + \left(1000000\right)^{2}\right)}{2}\simeq 1.7708\times 10^{-11}.
\end{aligned}
$$

同じ周波数の二因子から作る項である。

$$
\begin{aligned}
P_{2\omega_1}&=\frac{\epsilon_0\chi^{(2)} E_1^{2}}{2},\\
P_{2\omega_1}/(\mathrm{C/m^2})&\simeq \frac{\left(1.770837564\times 10^{-23}\right) \cdot \left(1000000\right)^{2}}{2}\simeq 8.8542\times 10^{-12}.
\end{aligned}
$$

位相ゼロの上の波と下の波を掛ける。

$$
\begin{aligned}
P_\Sigma&=\epsilon_0\chi^{(2)} E_1 E_2 e^{i \phi_2},\\
P_\Sigma/(\mathrm{C/m^2})&\simeq \left(1.770837564\times 10^{-23}\right) \cdot \left(1000000\right) \cdot \left(1000000\right) \cdot e^{i \cdot \left(0\right)}\simeq 1.7708\times 10^{-11}.
\end{aligned}
$$

差の項では下の波の共役を掛ける。

$$
\begin{aligned}
P_\Delta&=\epsilon_0\chi^{(2)} E_1 E_2 e^{- i \phi_2},\\
P_\Delta/(\mathrm{C/m^2})&\simeq \left(1.770837564\times 10^{-23}\right) \cdot \left(1000000\right) \cdot \left(1000000\right) \cdot e^{- i \cdot \left(0\right)}\simeq 1.7708\times 10^{-11}.
\end{aligned}
$$

和と差の源分極は同じ値であり、この例では上の波の高調波の二倍になる。
:::

:::example[一方を四分の一周期進める]
実振幅は前例のまま位相だけを動かす。

入力は上の周波数の実振幅 $E_1=1000000\,\mathrm{V/m}$、下の周波数の実振幅 $E_2=1000000\,\mathrm{V/m}$、下の周波数の位相 $\phi_2\simeq 1.5708\,\mathrm{rad}$、二次感受率 $\chi^{(2)}\simeq 2\times 10^{-12}\,\mathrm{m/V}$ である。

ここでは上の波を $E_1\cos\omega_1 t$、下の波を $E_2\cos(\omega_2t-\phi_2)$ とし、$E_1,E_2$ は実数とする。真空誘電率は $\epsilon_0\simeq8.8542\times10^{-12}\,\mathrm{F/m}$ である。

$$
\begin{aligned}
\epsilon_0\chi^{(2)}&=\epsilon_0\cdot\chi^{(2)},\\
\epsilon_0\chi^{(2)}/(\mathrm{F/V})&\simeq (8.854187819\times10^{-12})(2\times10^{-12})\simeq 1.7708\times 10^{-23}.
\end{aligned}
$$

各周波数とその共役の積を足す。

$$
\begin{aligned}
P_{\rm dc}&=\frac{\epsilon_0\chi^{(2)} \left(E_1^{2} + E_2^{2}\right)}{2},\\
P_{\rm dc}/(\mathrm{C/m^2})&\simeq \frac{\left(1.770837564\times 10^{-23}\right) \cdot \left(\left(1000000\right)^{2} + \left(1000000\right)^{2}\right)}{2}\simeq 1.7708\times 10^{-11}.
\end{aligned}
$$

同じ周波数の二因子から作る項である。

$$
\begin{aligned}
P_{2\omega_1}&=\frac{\epsilon_0\chi^{(2)} E_1^{2}}{2},\\
P_{2\omega_1}/(\mathrm{C/m^2})&\simeq \frac{\left(1.770837564\times 10^{-23}\right) \cdot \left(1000000\right)^{2}}{2}\simeq 8.8542\times 10^{-12}.
\end{aligned}
$$

位相ゼロの上の波と下の波を掛ける。

$$
\begin{aligned}
P_\Sigma&=\epsilon_0\chi^{(2)} E_1 E_2 e^{i \phi_2},\\
P_\Sigma/(\mathrm{C/m^2})&\simeq \left(1.770837564\times 10^{-23}\right) \cdot \left(1000000\right) \cdot \left(1000000\right) \cdot e^{i \cdot \left(1.570796327\right)}\simeq 1.7708\times 10^{-11}i.
\end{aligned}
$$

差の項では下の波の共役を掛ける。

$$
\begin{aligned}
P_\Delta&=\epsilon_0\chi^{(2)} E_1 E_2 e^{- i \phi_2},\\
P_\Delta/(\mathrm{C/m^2})&\simeq \left(1.770837564\times 10^{-23}\right) \cdot \left(1000000\right) \cdot \left(1000000\right) \cdot e^{- i \cdot \left(1.570796327\right)}\simeq -1.7708\times 10^{-11}i.
\end{aligned}
$$

和と差の源の位相は反対向きに動き、直流分極は変わらない。
:::

## 第 2 高調波発生と位相整合

### 生成する場の位相をそろえる
実場を $E_j=\Re[A_j(z)e^{i(k_jz-\omega_jt)}]$ とする。周波数 $\omega_j$ の非線形分極の複素振幅を $\mathcal P_j(z)$ と定義し、そこに空間の速い位相も含める。横 Maxwell 方程式は $(\partial_z^2+k_j^2)(A_j e^{ik_jz})=-\mu_0\omega_j^2\mathcal P_j$。$|A_j''|\ll2k_j|A_j'|$ とし、$k_j=n_j\omega_j/c$ を入れると
$$
A_j'=i\frac{\omega_j}{2n_jc\eps_0}\mathcal P_j e^{-ik_jz}.
\label{eq-optics-nonlinear-envelope}
$$
周波数 $\omega$ の基本波の変化を無視し、$\omega_2=2\omega$ の成分に注目する。二乗の正周波数係数は $A_1^2e^{2ik_1z}/4$ であり、実部表示の複素振幅に戻すと $\mathcal P_2=\eps_0\chi^{(2)} A_1^2 e^{2ik_1z}/2$。$A_2(0)=0$、$\Delta k=2k_1-k_2$ として積分すると
$$
A_2(L)=i\frac{\omega\chi^{(2)}}{2n_2 c}A_1^2 L e^{i\Delta k L/2}\operatorname{sinc}(\Delta k L/2),
\label{eq-optics-shg}
$$
ここで $\operatorname{sinc}x=\sin x/x$、原点の値は 1 である。:term[第 2 高調波発生]{en="second-harmonic generation" yomi="だいにこうちょうははっせい"} の照度は $n_2\eps_0c|A_2|^2/2$ で、ビームの横断面で積分すればパワーになる。$\Delta k=0$ が :term[位相整合]{en="phase matching" yomi="いそうせいごう"} であり、結晶の複屈折を使って基本波と第二高調波の有効屈折率を合わせられる。

位相ずれが $\pi$ になる長さ $L_c=\pi/|\Delta k|$ をコヒーレンス長と呼ぶ。$L_c$ ごとに $\chi^{(2)}$ の符号を変えると、隣の区間からの生成振幅を再び同相にできる。周期 $2L_c$ の矩形の符号関数の第一 Fourier 係数は、二つの半区間を積分して絶対値 $2/\pi$ になる。これが:term[擬似位相整合]{en="quasi-phase matching" yomi="ぎじいそうせいごう"} である。大きく変換する場合は、この小信号の式のままでは基本波のエネルギー減少を扱えない。

:::example[長さの選択]
$\Delta k=200\,\mathrm{m^{-1}}$ なら $L_c\simeq15.7\,\mathrm{mm}$。非補償の結晶で最初の強度零点は $2L_c$ である。
:::
:::example[波長を半分にする]
真空波長 $1064\,\mathrm{nm}$ から生成する第 2 高調波は $532\,\mathrm{nm}$。ただし結晶内では二周波数の屈折率の比も波長比に入る。
:::

### 複屈折による整合角と小信号効率

一軸結晶で基本波は常光、第二高調波は異常光とする。光軸と波数のなす角を $\theta$、第二高調波での主屈折率を $n_o(2\omega),n_e(2\omega)$ とすれば、屈折率楕円体から
$$
\frac{1}{n_{\rm eff}(2\omega,\theta)^2}
=\frac{\cos^2\theta}{n_o(2\omega)^2}+\frac{\sin^2\theta}{n_e(2\omega)^2}
$$
である。従って整合角は右辺を $1/n_o(\omega)^2$ と置いて求める。$n_o(\omega)$ が $n_o(2\omega)$ と $n_e(2\omega)$ の間にある場合に実角が存在する。さらに分極を出力偏光方向へ射影した有効係数が非零でなければ光は生成しない。異常光では波数とエネルギー流がずれることがあるので、細いビームでは位相が整合しても横方向の重なりが長さとともに減る。

$d_{\rm eff}=\chi^{(2)}/2$ とし、平面波で基本波の減少を無視する。式 :eqref[eq-optics-shg] を二乗し、$|A_1|^2=2I_1/(n_1\eps_0 c)$ を代入すれば
$$
I_2=\frac{2\omega^2 |d_{\rm eff}|^2}{n_1^2 n_2\eps_0 c^3}
       I_1^2 L^2\operatorname{sinc}^2(\Delta kL/2).
$$
この式を用いるには、計算した $I_2$ が $I_1$ に比べて十分小さいことも確認する。

### 変換が強いときに残る保存量
損失なし、位相整合、実の $\chi^{(2)}$ とし、$A_1(0)=A>0$、$A_2(0)=0$。$A_1=u$ を実、$A_2=iv$、$v$ を実として、式 :eqref[eq-optics-nonlinear-envelope] は
$$
u'=-a uv,\qquad v'=b u^2,\qquad
a=\frac{\omega\chi^{(2)}}{2n_1c},\quad b=\frac{\omega\chi^{(2)}}{2n_2c}
$$
となる。基本波の源は $A_2A_1^*$ で、異なる周波数の二順序があるため二次分極は $\eps_0\chi^{(2)}A_2A_1^*$ である。従って $bu^2+av^2=bA^2$ が保存される。$a,b>0$ の場合、$v'=bA^2-av^2$ の分離解は $v=A\sqrt{b/a}\tanh(A\sqrt{ab}\,z)$。$u=A\operatorname{sech}(A\sqrt{ab}\,z)$ となり、強度換算では $n_2v^2+n_1u^2=n_1A^2$ である。短い距離で $\tanh x\simeq x$ とすれば小変換の式を回復する。

### パルスの群速度差
基本波と第二高調波の群速度が違うと、位相整合してもパルスが重なる時間が有限になる。逆群速度を $\beta_{11},\beta_{12}$ とし、分散二次以上を無視する。生成波の方程式は $\partial_z A_2+\beta_{12}\partial_t A_2=i b A_1^2e^{i\Delta kz}$、基本波は小変換で $A_1(z,t)=A_{10}(t-\beta_{11}z)$ とする。$T=t-\beta_{12}z$ を一定に積分すれば
$$
A_2(L,T)=ib\int_0^L A_{10}^2(T+(\beta_{12}-\beta_{11})z)e^{i\Delta kz}\dd z .
$$
群遅延差 $|\beta_{12}-\beta_{11}|L$ が入力幅より大きいと、同じ $T$ で全長の振幅は足されない。短パルスの結晶長を決める際は、$\Delta k$ の許容範囲だけでなく、この群遅延差も評価する。

:::example[線形の位相ずれによる強度低下]
$\Delta k L=\pi$ なら、長さと基本波が同じ完全整合例に対する照度比は $\operatorname{sinc}^2(\pi/2)=4/\pi^2\simeq0.405$。相互作用長そのものを半分にして不整合を半減させる場合は、$L^2$ の変化も同時に掛かる。
:::

### 周期反転による擬似位相整合

有効非線形係数を周期 $\Lambda$ で $+d,-d$ に半周期ずつ反転する。$m\ne0$ の複素 Fourier 係数は
$$
d_m=\frac{1}{\Lambda}\left[
\int_0^{\Lambda/2}d e^{-2\pi imz/\Lambda}\dd z-
\int_{\Lambda/2}^\Lambda d e^{-2\pi imz/\Lambda}\dd z
\right]=\frac{d(1-(-1)^m)}{i\pi m}.
$$
平均はゼロである。:eqref[eq-optics-nonlinear-envelope] の分極と伝播波の位相差に $2\pi m/\Lambda$ を加えたものがゼロなら、その項は距離に比例して蓄積する。一次の大きさは $2d/\pi$ で、一様結晶が完全整合する場合に対する小信号パワーの比は $4/\pi^2$。整合からの残差を $\delta k$ とすると有限長での比には $\operatorname{sinc}^2(\delta k L/2)$ がさらに掛かる。この一項の近似は、一周期当たりの変換が弱く、多くの周期を通る場合に使う。

:::example[一次で不整合を打ち消す]
指定の波数差に周期を合わせる。

入力は元の波数の不整合 $\Delta k=200000\,\mathrm{m^{-1}}$、反転周期 $\Lambda\simeq 3.1416\times 10^{-5}\,\mathrm m$、結晶長 $L=\frac{1}{250}\,\mathrm m$、選ぶ正の反転次数 $m=1$ である。

周期反転が持つ基本波数を求める。

$$
\begin{aligned}
K&=\frac{2 \pi}{\Lambda},\\
K/(\mathrm{m^{-1}})&\simeq \frac{2 \cdot \pi}{\left(3.141592654\times 10^{-5}\right)}=200000.
\end{aligned}
$$

整合へ最も近い Fourier 項を選ぶ。

$$
\begin{aligned}
\delta k&=- K m + \Delta k,\\
\delta k/(\mathrm{m^{-1}})&= - \left(200000\right) \cdot \left(1\right) + \left(200000\right)=0.
\end{aligned}
$$

方形波で選んだ次数の係数の絶対値を取る。

$$
\begin{aligned}
d_{\rm eff}/d&=\frac{2 \left|{\sin{\left(\frac{\pi m}{2} \right)}}\right|}{\pi m},\\
d_{\rm eff}/d&= \frac{2 \cdot \left|{\sin{\left(\frac{\pi \cdot \left(1\right)}{2} \right)}}\right|}{\pi \cdot \left(1\right)}\simeq 0.63662.
\end{aligned}
$$

有限長で残差位相が回る積分を二乗する。

$$
\begin{aligned}
q&=\operatorname{sinc}^{2}{\left(\frac{L \delta k}{2} \right)},\\
q&= \operatorname{sinc}^{2}{\left(\frac{\left(\frac{1}{250}\right) \cdot \left(0\right)}{2} \right)}=1.
\end{aligned}
$$

係数のパワー比と残差位相の比を掛ける。

$$
\begin{aligned}
R&=(d_{\rm eff}/d)^2 q,\\
R&\simeq \left(0.6366197724\right)^{2} \cdot \left(1\right)\simeq 0.40528.
\end{aligned}
$$

位相の残差はゼロになるが、矩形の反転パターンにより係数は低下する。
:::

:::example[反転周期がずれた結晶]
同じ波数差で周期を約一パーセント伸ばす。

入力は元の波数の不整合 $\Delta k=200000\,\mathrm{m^{-1}}$、反転周期 $\Lambda\simeq 3.173\times 10^{-5}\,\mathrm m$、結晶長 $L=\frac{1}{250}\,\mathrm m$、選ぶ正の反転次数 $m=1$ である。

周期反転が持つ基本波数を求める。

$$
\begin{aligned}
K&=\frac{2 \pi}{\Lambda},\\
K/(\mathrm{m^{-1}})&\simeq \frac{2 \cdot \pi}{\left(3.17300858\times 10^{-5}\right)}\simeq 1.9802\times 10^{5}.
\end{aligned}
$$

整合へ最も近い Fourier 項を選ぶ。

$$
\begin{aligned}
\delta k&=- K m + \Delta k,\\
\delta k/(\mathrm{m^{-1}})&\simeq - \left(1.98019802\times 10^{5}\right) \cdot \left(1\right) + \left(200000\right)\simeq 1980.2.
\end{aligned}
$$

方形波で選んだ次数の係数の絶対値を取る。

$$
\begin{aligned}
d_{\rm eff}/d&=\frac{2 \left|{\sin{\left(\frac{\pi m}{2} \right)}}\right|}{\pi m},\\
d_{\rm eff}/d&= \frac{2 \cdot \left|{\sin{\left(\frac{\pi \cdot \left(1\right)}{2} \right)}}\right|}{\pi \cdot \left(1\right)}\simeq 0.63662.
\end{aligned}
$$

有限長で残差位相が回る積分を二乗する。

$$
\begin{aligned}
q&=\operatorname{sinc}^{2}{\left(\frac{L \delta k}{2} \right)},\\
q&\simeq \operatorname{sinc}^{2}{\left(\frac{\left(\frac{1}{250}\right) \cdot \left(1980.19802\right)}{2} \right)}\simeq 0.034006.
\end{aligned}
$$

係数のパワー比と残差位相の比を掛ける。

$$
\begin{aligned}
R&=(d_{\rm eff}/d)^2 q,\\
R&\simeq \left(0.6366197724\right)^{2} \cdot \left(0.03400633264\right)\simeq 0.013782.
\end{aligned}
$$

小さな周期差が結晶全体で蓄積し、変換効率を下げる。
:::

:::example[三次の周期を使う]
長い周期に変更し、三次の成分で同じ波数差を補う。

入力は元の波数の不整合 $\Delta k=200000\,\mathrm{m^{-1}}$、反転周期 $\Lambda\simeq 9.4248\times 10^{-5}\,\mathrm m$、結晶長 $L=\frac{1}{250}\,\mathrm m$、選ぶ正の反転次数 $m=3$ である。

周期反転が持つ基本波数を求める。

$$
\begin{aligned}
K&=\frac{2 \pi}{\Lambda},\\
K/(\mathrm{m^{-1}})&\simeq \frac{2 \cdot \pi}{\left(9.424777961\times 10^{-5}\right)}\simeq 66667.
\end{aligned}
$$

整合へ最も近い Fourier 項を選ぶ。

$$
\begin{aligned}
\delta k&=- K m + \Delta k,\\
\delta k/(\mathrm{m^{-1}})&\simeq -3\left(\frac{200000}{3}\right)+200000=0.
\end{aligned}
$$

方形波で選んだ次数の係数の絶対値を取る。

$$
\begin{aligned}
d_{\rm eff}/d&=\frac{2 \left|{\sin{\left(\frac{\pi m}{2} \right)}}\right|}{\pi m},\\
d_{\rm eff}/d&= \frac{2 \cdot \left|{\sin{\left(\frac{\pi \cdot \left(3\right)}{2} \right)}}\right|}{\pi \cdot \left(3\right)}\simeq 0.21221.
\end{aligned}
$$

有限長で残差位相が回る積分を二乗する。

$$
\begin{aligned}
q&=\operatorname{sinc}^{2}{\left(\frac{L \delta k}{2} \right)},\\
q&= \operatorname{sinc}^{2}{\left(\frac{\left(\frac{1}{250}\right) \cdot \left(0\right)}{2} \right)}=1.
\end{aligned}
$$

係数のパワー比と残差位相の比を掛ける。

$$
\begin{aligned}
R&=(d_{\rm eff}/d)^2 q,\\
R&\simeq \left(0.2122065908\right)^{2} \cdot \left(1\right)\simeq 0.045032.
\end{aligned}
$$

残差を消しても結合係数の二乗は一次の九分の一になる。
:::

## 和周波・差周波発生とパラメトリック増幅

### 三つの波の連立式
$\omega_3=\omega_1+\omega_2$、$\omega_1\ne\omega_2$ とする。実場の二乗では、$A_1A_2$ と $A_3A_2^*$ の積がそれぞれ :term[和周波発生]{en="sum-frequency generation" yomi="わしゅうははっせい"} と :term[差周波発生]{en="difference-frequency generation" yomi="さしゅうははっせい"} の源になる。積の二通りの順序と実部表示への換算を行い、$\Delta k=k_1+k_2-k_3$、$c_j=\omega_j\chi^{(2)}/(2n_jc)$ として :eqref[eq-optics-nonlinear-envelope] に入れる。
$$
\begin{aligned}
A_3'&=ic_3A_1A_2e^{i\Delta kz},\\
A_1'&=ic_1A_3A_2^*e^{-i\Delta kz},&
A_2'&=ic_2A_3A_1^*e^{-i\Delta kz}.
\end{aligned}
$$
ここでは同じ実 $\chi^{(2)}$ を使う無損失模型である。$Z=A_3^*A_1A_2e^{i\Delta kz}$ と置けば $\partial_z|A_3|^2=-2c_3\Im Z$、低い二波は $+2c_j\Im Z$。光子束 $\Phi_j=n_j\eps_0c|A_j|^2/(2\hbar\omega_j)$ に直すと $\Phi_1'=\Phi_2'=-\Phi_3'$ となる。これは Manley–Rowe の関係で、高い周波数の光子一つが低い二周波数の光子一つずつと対応する。各エネルギーで重み付けして足すと全強度は保存される。

### ポンプ一定の線形化
$A_3$ を一定のポンプとし、位相整合させる。$A_1'$ を微分して $A_2'^*$ を代入すると $A_1''=G^2 A_1$、$G=\sqrt{c_1c_2}|A_3|$ である。$A_2(0)=0$ なら $A_1'(0)=0$ なので $A_1(z)=A_1(0)\cosh(Gz)$。$A_2$ は $\sinh(Gz)$ に比例する。これが :term[パラメトリック増幅]{en="parametric amplification" yomi="ぱらめとりっくぞうふく"} である。不整合があるとき $\widetilde A_1=A_1e^{i\Delta kz/2}$、$\widetilde A_2^*=A_2^*e^{-i\Delta kz/2}$ と置くと、微分行列の対角は $\pm i\Delta k/2$、非対角の積は $G^2$。従って固有値の二乗は $G^2-(\Delta k/2)^2$ となる。

:::example[もう一つの周波数]
ポンプが $300\,\mathrm{THz}$、信号が $190\,\mathrm{THz}$ ならもう一波は $110\,\mathrm{THz}$ である。
:::
:::example[強度利得]
$GL=\ln2$ の位相整合した増幅器で、種のない側の入力をゼロとする。信号強度比は $\cosh^2(\ln2)=25/16$ である。
:::

### 位相ずれを持つ増幅の解
ポンプ一定とする。前の変換後の二成分を $v$、その微分行列を $M$ とすると $M^2=g^2 I$、$g^2=G^2-(\Delta k/2)^2$ である。指数の偶奇を分けて $v(L)=[\cosh(gL)I+\sinh(gL)M/g]v(0)$。第二入力がゼロのとき信号強度比は
$$
1+\frac{G^2}{g^2}\sinh^2(gL)
$$
である。$g$ が純虚数なら $\sinh(iu)=i\sin u$ を使い、強度比は振動する。$g=0$ では極限 $1+G^2L^2$ を取る。$G$ が一定でも、位相ずれが強ければ長距離の指数的増幅は失われる。

### 位相整合の波長許容幅
ポンプ周波数を固定して信号角周波数 $\omega_1$ を動かすと、もう一波は $\omega_2=\omega_3-\omega_1$ と動く。二波の逆群速度をそれぞれ $k_1',k_2'$ とすると、整合点からの離調による一次の不整合は $\Delta k\simeq(k_1'-k_2')\delta\omega_1$ である。小変換で $\operatorname{sinc}^2(\Delta k L/2)$ の片側第一零点を使えば、離調の大きさは $2\pi/(L|k_1'-k_2'|)$。群速度が一致すると一次項が消えるので、二階微分やポンプ帯域で許容幅が決まる。

:::example[整合から外れて成長が止まる点]
$G=5\,\mathrm{m^{-1}}$、$\Delta k=10\,\mathrm{m^{-1}}$ では $g=0$。もう一波の入力がなく $L=0.10\,\mathrm m$ なら信号利得は $1+G^2L^2=1.25$。同じ条件で $\Delta k$ をこれ以上増すと、長くしても指数増幅は得られない。
:::

### 三波混合で光子流を数える

電力を周波数で割る保存関係を導くため、三波の振幅を光子流で規格化し、高周波を $p$、低い二波を $s,i$ とする。位相整合した無損失の式を $a_p'=ig a_s a_i$、$a_s'=ig a_p a_i^*$、$a_i'=ig a_p a_s^*$ と書く。$g$ を実とし $C=a_p^*a_s a_i$ とおくと、$\partial_z|a_p|^2=-2g\Im C$、残り二つは $+2g\Im C$。従って失われたポンプ光子一個当たり、信号とアイドラーが一個ずつ増える。これは :term[Manley–Rowe 関係]{en="Manley–Rowe relations" yomi="まんりーろーかんけい"} である。パワーに戻すと $\omega_p=\omega_s+\omega_i$ により全パワーも保存する。位相不整合を入れる場合は $C$ に相応の位相因子を付けるが、光子束の増分の関係と、全パワーの保存は同じである。

:::example[縮退した二モードの生成]
ここでは信号とアイドラーを直交する偏光で区別する。周波数が等しくても互いに異なるモードなので、差周波数発生に対応する二つの光子束を別に数えられる。
信号をポンプの半周波数に選ぶ。二波は偏光など別の自由度で区別する。

入力はポンプ真空波長 $\lambda_p\simeq 5.32\times 10^{-7}\,\mathrm m$、信号真空波長 $\lambda_s\simeq 1.064\times 10^{-6}\,\mathrm m$、失ったポンプパワー $\Delta P_p=\frac{1}{10}\,\mathrm W$ である。

エネルギーの和を各真空波長の逆数で表す。

$$
\begin{aligned}
\lambda_i&=\frac{1}{- \frac{1}{\lambda_s} + \frac{1}{\lambda_p}},\\
\lambda_i/(\mathrm m)&\simeq \frac{1}{- \frac{1}{\left(1.064\times 10^{-6}\right)} + \frac{1}{\left(5.32\times 10^{-7}\right)}}\simeq 1.064\times 10^{-6}.
\end{aligned}
$$

失った光パワーを一ポンプ光子のエネルギーで割る。

$$
\begin{aligned}
N&=\frac{\Delta P_p \lambda_p}{hc},\\
N/(\mathrm{s^{-1}})&\simeq \frac{10000000000000000000000000000000000 \cdot \left(\frac{1}{10}\right) \cdot \left(5.32\times 10^{-7}\right)}{1986445857}\simeq 2.6782\times 10^{17}.
\end{aligned}
$$

同じ個数を信号一光子のエネルギーで換算する。

$$
\begin{aligned}
\Delta P_s&=\frac{\Delta P_p \lambda_p}{\lambda_s},\\
\Delta P_s/(\mathrm W)&\simeq \frac{\left(\frac{1}{10}\right) \cdot \left(5.32\times 10^{-7}\right)}{\left(1.064\times 10^{-6}\right)}=\frac{1}{20}.
\end{aligned}
$$

残りの周波数側でも同じ光子数を使う。

$$
\begin{aligned}
\Delta P_i&=\frac{\Delta P_p \lambda_p}{\lambda_i},\\
\Delta P_i/(\mathrm W)&\simeq \frac{\left(\frac{1}{10}\right) \cdot \left(5.32\times 10^{-7}\right)}{\left(1.064\times 10^{-6}\right)}=\frac{1}{20}.
\end{aligned}
$$

等しい光子数であっても、両周波数のパワーは一般に異なる。

$$
\begin{aligned}
R&=\frac{\Delta P_s}{\Delta P_i + \Delta P_s},\\
R&= \frac{\left(\frac{1}{20}\right)}{\left(\frac{1}{20}\right) + \left(\frac{1}{20}\right)}=\frac{1}{2}.
\end{aligned}
$$

信号とアイドラーの波長が一致する。二つのラベルで数える理想模型ではパワーが等分される。
:::

:::example[高周波側の信号を選ぶ]
信号の波長を短くし、別のアイドラーを生成する。

入力はポンプ真空波長 $\lambda_p\simeq 5.32\times 10^{-7}\,\mathrm m$、信号真空波長 $\lambda_s\simeq 8\times 10^{-7}\,\mathrm m$、失ったポンプパワー $\Delta P_p=\frac{1}{10}\,\mathrm W$ である。

エネルギーの和を各真空波長の逆数で表す。

$$
\begin{aligned}
\lambda_i&=\frac{1}{- \frac{1}{\lambda_s} + \frac{1}{\lambda_p}},\\
\lambda_i/(\mathrm m)&\simeq \frac{1}{- \frac{1}{\left(8\times 10^{-7}\right)} + \frac{1}{\left(5.32\times 10^{-7}\right)}}\simeq 1.5881\times 10^{-6}.
\end{aligned}
$$

失った光パワーを一ポンプ光子のエネルギーで割る。

$$
\begin{aligned}
N&=\frac{\Delta P_p \lambda_p}{hc},\\
N/(\mathrm{s^{-1}})&\simeq \frac{10000000000000000000000000000000000 \cdot \left(\frac{1}{10}\right) \cdot \left(5.32\times 10^{-7}\right)}{1986445857}\simeq 2.6782\times 10^{17}.
\end{aligned}
$$

同じ個数を信号一光子のエネルギーで換算する。

$$
\begin{aligned}
\Delta P_s&=\frac{\Delta P_p \lambda_p}{\lambda_s},\\
\Delta P_s/(\mathrm W)&\simeq \frac{\left(\frac{1}{10}\right) \cdot \left(5.32\times 10^{-7}\right)}{\left(8\times 10^{-7}\right)}=\frac{133}{2000}.
\end{aligned}
$$

残りの周波数側でも同じ光子数を使う。

$$
\begin{aligned}
\Delta P_i&=\frac{\Delta P_p \lambda_p}{\lambda_i},\\
\Delta P_i/(\mathrm W)&\simeq \frac{\left(\frac{1}{10}\right) \cdot \left(5.32\times 10^{-7}\right)}{\left(1.588059701\times 10^{-6}\right)}=\frac{67}{2000}.
\end{aligned}
$$

等しい光子数であっても、両周波数のパワーは一般に異なる。

$$
\begin{aligned}
R&=\frac{\Delta P_s}{\Delta P_i + \Delta P_s},\\
R&= \frac{\left(\frac{133}{2000}\right)}{\left(\frac{67}{2000}\right) + \left(\frac{133}{2000}\right)}=\frac{133}{200}.
\end{aligned}
$$

アイドラーより信号に多くのパワーが移るが、増えた光子数は同じである。
:::

:::example[低周波側の信号を選ぶ]
信号をさらに長い波長へ移し、パワーの配分を逆にする。

入力はポンプ真空波長 $\lambda_p\simeq 5.32\times 10^{-7}\,\mathrm m$、信号真空波長 $\lambda_s\simeq 1.6\times 10^{-6}\,\mathrm m$、失ったポンプパワー $\Delta P_p=\frac{1}{10}\,\mathrm W$ である。

エネルギーの和を各真空波長の逆数で表す。

$$
\begin{aligned}
\lambda_i&=\frac{1}{- \frac{1}{\lambda_s} + \frac{1}{\lambda_p}},\\
\lambda_i/(\mathrm m)&\simeq \frac{1}{- \frac{1}{\left(1.6\times 10^{-6}\right)} + \frac{1}{\left(5.32\times 10^{-7}\right)}}\simeq 7.97\times 10^{-7}.
\end{aligned}
$$

失った光パワーを一ポンプ光子のエネルギーで割る。

$$
\begin{aligned}
N&=\frac{\Delta P_p \lambda_p}{hc},\\
N/(\mathrm{s^{-1}})&\simeq \frac{10000000000000000000000000000000000 \cdot \left(\frac{1}{10}\right) \cdot \left(5.32\times 10^{-7}\right)}{1986445857}\simeq 2.6782\times 10^{17}.
\end{aligned}
$$

同じ個数を信号一光子のエネルギーで換算する。

$$
\begin{aligned}
\Delta P_s&=\frac{\Delta P_p \lambda_p}{\lambda_s},\\
\Delta P_s/(\mathrm W)&\simeq \frac{\left(\frac{1}{10}\right) \cdot \left(5.32\times 10^{-7}\right)}{\left(1.6\times 10^{-6}\right)}=\frac{133}{4000}.
\end{aligned}
$$

残りの周波数側でも同じ光子数を使う。

$$
\begin{aligned}
\Delta P_i&=\frac{\Delta P_p \lambda_p}{\lambda_i},\\
\Delta P_i/(\mathrm W)&\simeq \frac{\left(\frac{1}{10}\right) \cdot \left(5.32\times 10^{-7}\right)}{\left(7.970037453\times 10^{-7}\right)}=\frac{267}{4000}.
\end{aligned}
$$

等しい光子数であっても、両周波数のパワーは一般に異なる。

$$
\begin{aligned}
R&=\frac{\Delta P_s}{\Delta P_i + \Delta P_s},\\
R&= \frac{\left(\frac{133}{4000}\right)}{\left(\frac{267}{4000}\right) + \left(\frac{133}{4000}\right)}=\frac{133}{400}.
\end{aligned}
$$

アイドラーがより高い周波数になる。
:::

## 3 次の非線形効果（光 Kerr 効果・自己位相変調）

### 強度依存の屈折率
$E=\Re(Ae^{ikz-i\omega t})$ の三乗を展開し、基本波の複素振幅だけを取ると $3|A|^2A/4$ になる。線形屈折率を $n_0$ とすれば、$n^2$ の補正は $3\chi^{(3)}|A|^2/4$。実部に注目して平方根を展開し $I=n_0\eps_0c|A|^2/2$ を代入すると、
$$
n=n_0+n_{\rm K} I,\qquad n_{\rm K}=\frac{3\Re\chi^{(3)}}{4n_0^2\eps_0 c}.
\label{eq-optics-kerr}
$$
これが :term[光 Kerr 効果]{en="optical Kerr effect" yomi="ひかりかーこうか"} である。$n_{\rm K}$ の単位は $\mathrm{m^2/W}$。パルスの強度が遅延時間 $T$ で変われば、位相 $\phi_{\rm NL}(T)=k_0 n_{\rm K} I(T)L$ も変わる。これを :term[自己位相変調]{en="self-phase modulation" yomi="じこいそうへんちょう"} と呼ぶ。瞬時周波数の偏移は規約 $e^{-i\omega t}$ で $-\partial_T\phi_{\rm NL}$ である。そのため強度包絡線が直ちに変わらなくても、スペクトルが広がり得る。強度分布が横方向で変われば位相にも横勾配がつくため、ビームの集光にもつながる。

同じスカラー導波モード $u(x,y)$ に沿って伝播し材料の $n_{\rm K}$ が横で一定なら、分極をモードへ射影すると重み $\int|u|^4/(\int |u|^2)^2=1/A_{\rm eff}$ を得る。パワー $P$ に対する位相は $\gamma PL$、$\gamma=k_0 n_{\rm K}/A_{\rm eff}$ である。

:::example[量の尺度]
$\gamma=0.002\,\mathrm{W^{-1}m^{-1}}$、$P=10\,\mathrm W$、$L=50\,\mathrm m$ なら $\phi_{\rm NL}=1\,\mathrm{rad}$。
:::
:::example[前縁と後縁]
$n_{\rm K}>0$ では強度の増す前縁で $-\partial_T\phi_{\rm NL}<0$、減る後縁で正である。パルスの中央の強度極大点ではその微分はゼロである。
:::

### 二色の強度が互いの位相を変える
二つの異なる搬送周波数をもつ場を三乗する。第三高調波や他の混合成分が二つの搬送周波数と重ならないとする。このとき第一波の周波数に残る組は $(\omega_1,\omega_1,-\omega_1)$ の三順序と $(\omega_1,\omega_2,-\omega_2)$ の六順序である。周波数の組に依存しない瞬時スカラー $\chi^{(3)}$、同偏光、近い線形屈折率を仮定すると、第一波の位相増分は $k_{01} n_{\rm K}(I_1+2I_2)\dd z$ である。第二波の強度による項が交差位相変調である。周波数の差が群遅延を生むと、同じ時刻の包絡線が全長で重ならないので、この項は各位置の強度で積分する。偏光を変えればテンソルの成分も変わる。

### 四光波混合の波数条件
異なる三つの周波数の電場の三乗には $A_1A_2A_3^*$ があり、$\omega_4=\omega_1+\omega_2-\omega_3$ の正周波数を駆動する。源の空間位相は $k_1+k_2-k_3$ なので、$k_4$ との差を取った量が位相不整合である。基本波二個が同じ周波数の場合は $2\omega_p=\omega_s+\omega_i$ と書け、$A_p^2 A_s^*$ がもう一波を作る。屈折率が強度でも動くため、強ポンプの場合は線形の波数だけでなく、自己位相変調と交差位相変調が出す位相差も積分する。これは同じ三次分極の中の異なる周波数組を連立させた結果である。

:::example[パルスの最も急な部分での周波数偏移]
ガウスパルスの軸上の最大位相が $\phi_0=1\,\mathrm{rad}$、$1/e$ 強度半幅が $T_0=100\,\mathrm{fs}$ のとき、偏移の最大は $\sqrt{2/e}\phi_0/T_0$ である。$|T|e^{-T^2/T_0^2}$ を微分して $|T|=T_0/\sqrt2$ を得る。Hz への換算は $2\pi$ で割るので、最大は約 $1.37\,\mathrm{THz}$。
:::

### Kerr 媒質で近接した側帯が成長する条件

連続波ポンプのパワーを $P$ とし、損失のない単一モードの Kerr 媒質で、群速度分散 $\beta_2$ までを残す。場を $A=e^{i\gamma Pz}(\sqrt{P}+a)$ と書き、分散による位相 $D=\beta_2\Omega^2/2$ を持つ側帯 $\pm\Omega$ の小振幅を $a_\pm$ とする。分極の三次項の線形化では $|\sqrt P+a|^2(\sqrt P+a)=P\sqrt P+2Pa+Pa^*+O(|a|^2)$。共通のポンプ位相を除いた後の式は
$$
\frac{\dd}{\dd z}\binom{a_+}{a_-^*}
=i\begin{pmatrix} D+\gamma P & \gamma P\\ -\gamma P &-D-\gamma P \end{pmatrix}
 \binom{a_+}{a_-^*}.
$$
行列の二乗から、成長指数 $g$ の二乗は $(\gamma P)^2-(D+\gamma P)^2$。例えば $\gamma>0,\ \beta_2<0$ では $-2\gamma P<D<0$ で指数成長し、$D=-\gamma P$ で最大の $g=\gamma P$ に達する。これが :term[変調不安定性]{en="modulation instability" yomi="へんちょうふあんていせい"} である。この結果は連続波近傍の微小摂動に対し、強い成長でポンプが減る前に適用する。

### 損失がある Kerr 媒質で位相を積分する

単一モードを伝わるパルスで分散を無視し、パワー損失率を $\alpha$、非線形係数を $\gamma$ とする。遅延座標で $P(z,t)=P(0,t)e^{-\alpha z}$、位相の増加率は $\gamma P(z,t)$ なので $\phi(t)=\gamma P(0,t)L_{\rm eff}$、$L_{\rm eff}=(1-e^{-\alpha L})/\alpha$。$\alpha\to0$ では $L_{\rm eff}\to L$ である。入射パワーが $P_0e^{-t^2/T^2}$ なら瞬間角周波数のずれ $-\partial_t\phi=2t\phi(0)e^{-t^2/T^2}/T^2$。その絶対値は $|t|=T/\sqrt2$ で最大になり、裾でゼロへ戻る。これは瞬間周波数の範囲であり、スペクトルの半値全幅そのものにはならない。

:::example[短い弱損失区間]
パルスに非線形位相を付け、短距離で損失が効くかを調べる。

入力はパワー損失率 $\alpha=\frac{1}{10}\,\mathrm{m^{-1}}$、伝播長 $L=1\,\mathrm m$、非線形係数 $\gamma=\frac{1}{100}\,\mathrm{W^{-1}m^{-1}}$、入射ピークパワー $P_0=5\,\mathrm W$、パワーの $1/e$ 半幅 $T\simeq 1\times 10^{-11}\,\mathrm s$ である。

伝播中のパワー減衰を出す。

$$
\begin{aligned}
r&=e^{- L \alpha},\\
r&= e^{- \left(1\right) \cdot \left(\frac{1}{10}\right)}\simeq 0.90484.
\end{aligned}
$$

軸に沿うパワーを積分し、入力値で規格化する。

$$
\begin{aligned}
L_{\rm eff}&=\frac{1 - r}{\alpha},\\
L_{\rm eff}/(\mathrm m)&\simeq \frac{1 - \left(0.904837418\right)}{\left(\frac{1}{10}\right)}\simeq 0.95163.
\end{aligned}
$$

中心で最も大きなパワーに有効長を掛ける。

$$
\begin{aligned}
\phi_0&=P_0 L_{\rm eff} \gamma,\\
\phi_0/(\mathrm{rad})&\simeq \left(5\right) \cdot \left(0.9516258196\right) \cdot \left(\frac{1}{100}\right)\simeq 0.047581.
\end{aligned}
$$

Gaussian の傾きが最大となる点を取る。

$$
\begin{aligned}
t_{\max}&=\frac{\sqrt{2} T}{2},\\
t_{\max}/(\mathrm s)&\simeq \frac{\sqrt{2} \cdot \left(1\times 10^{-11}\right)}{2}\simeq 7.0711\times 10^{-12}.
\end{aligned}
$$

時間の位相微分を通常周波数へ換える。

$$
\begin{aligned}
\Delta\nu_{\max}&=\frac{\sqrt{2} \left|{\phi_0}\right|}{2 \pi T e^{\frac{1}{2}}},\\
\Delta\nu_{\max}/(\mathrm{Hz})&\simeq \frac{\sqrt{2} \cdot \left|{\left(0.04758129098\right)}\right|}{2 \cdot \pi \cdot \left(1\times 10^{-11}\right) \cdot e^{\frac{1}{2}}}\simeq 6.4957\times 10^{8}.
\end{aligned}
$$

長さと有効長は近いが、位相はなお有限である。
:::

:::example[同じ区間を長くする]
材料と入力はそのまま伝播長を大きくする。

入力はパワー損失率 $\alpha=\frac{1}{10}\,\mathrm{m^{-1}}$、伝播長 $L=20\,\mathrm m$、非線形係数 $\gamma=\frac{1}{100}\,\mathrm{W^{-1}m^{-1}}$、入射ピークパワー $P_0=5\,\mathrm W$、パワーの $1/e$ 半幅 $T\simeq 1\times 10^{-11}\,\mathrm s$ である。

伝播中のパワー減衰を出す。

$$
\begin{aligned}
r&=e^{- L \alpha},\\
r&= e^{- \left(20\right) \cdot \left(\frac{1}{10}\right)}\simeq 0.13534.
\end{aligned}
$$

軸に沿うパワーを積分し、入力値で規格化する。

$$
\begin{aligned}
L_{\rm eff}&=\frac{1 - r}{\alpha},\\
L_{\rm eff}/(\mathrm m)&\simeq \frac{1 - \left(0.1353352832\right)}{\left(\frac{1}{10}\right)}\simeq 8.6466.
\end{aligned}
$$

中心で最も大きなパワーに有効長を掛ける。

$$
\begin{aligned}
\phi_0&=P_0 L_{\rm eff} \gamma,\\
\phi_0/(\mathrm{rad})&\simeq \left(5\right) \cdot \left(8.646647168\right) \cdot \left(\frac{1}{100}\right)\simeq 0.43233.
\end{aligned}
$$

Gaussian の傾きが最大となる点を取る。

$$
\begin{aligned}
t_{\max}&=\frac{\sqrt{2} T}{2},\\
t_{\max}/(\mathrm s)&\simeq \frac{\sqrt{2} \cdot \left(1\times 10^{-11}\right)}{2}\simeq 7.0711\times 10^{-12}.
\end{aligned}
$$

時間の位相微分を通常周波数へ換える。

$$
\begin{aligned}
\Delta\nu_{\max}&=\frac{\sqrt{2} \left|{\phi_0}\right|}{2 \pi T e^{\frac{1}{2}}},\\
\Delta\nu_{\max}/(\mathrm{Hz})&\simeq \frac{\sqrt{2} \cdot \left|{\left(0.4323323584\right)}\right|}{2 \cdot \pi \cdot \left(1\times 10^{-11}\right) \cdot e^{\frac{1}{2}}}\simeq 5.9021\times 10^{9}.
\end{aligned}
$$

有効長は幾何長より短い。入力パワーを全長に掛けると位相を過大評価する。
:::

:::example[パルスを長くする]
前例の中心パワーを保ち、時間の幅を伸ばす。

入力はパワー損失率 $\alpha=\frac{1}{10}\,\mathrm{m^{-1}}$、伝播長 $L=20\,\mathrm m$、非線形係数 $\gamma=\frac{1}{100}\,\mathrm{W^{-1}m^{-1}}$、入射ピークパワー $P_0=5\,\mathrm W$、パワーの $1/e$ 半幅 $T\simeq 2\times 10^{-11}\,\mathrm s$ である。

伝播中のパワー減衰を出す。

$$
\begin{aligned}
r&=e^{- L \alpha},\\
r&= e^{- \left(20\right) \cdot \left(\frac{1}{10}\right)}\simeq 0.13534.
\end{aligned}
$$

軸に沿うパワーを積分し、入力値で規格化する。

$$
\begin{aligned}
L_{\rm eff}&=\frac{1 - r}{\alpha},\\
L_{\rm eff}/(\mathrm m)&\simeq \frac{1 - \left(0.1353352832\right)}{\left(\frac{1}{10}\right)}\simeq 8.6466.
\end{aligned}
$$

中心で最も大きなパワーに有効長を掛ける。

$$
\begin{aligned}
\phi_0&=P_0 L_{\rm eff} \gamma,\\
\phi_0/(\mathrm{rad})&\simeq \left(5\right) \cdot \left(8.646647168\right) \cdot \left(\frac{1}{100}\right)\simeq 0.43233.
\end{aligned}
$$

Gaussian の傾きが最大となる点を取る。

$$
\begin{aligned}
t_{\max}&=\frac{\sqrt{2} T}{2},\\
t_{\max}/(\mathrm s)&\simeq \frac{\sqrt{2} \cdot \left(2\times 10^{-11}\right)}{2}\simeq 1.4142\times 10^{-11}.
\end{aligned}
$$

時間の位相微分を通常周波数へ換える。

$$
\begin{aligned}
\Delta\nu_{\max}&=\frac{\sqrt{2} \left|{\phi_0}\right|}{2 \pi T e^{\frac{1}{2}}},\\
\Delta\nu_{\max}/(\mathrm{Hz})&\simeq \frac{\sqrt{2} \cdot \left|{\left(0.4323323584\right)}\right|}{2 \cdot \pi \cdot \left(2\times 10^{-11}\right) \cdot e^{\frac{1}{2}}}\simeq 2.951\times 10^{9}.
\end{aligned}
$$

中心の位相は変わらず、瞬間周波数の最大偏移は半分になる。
:::

## 光ソリトン

### 形を保つパルス
$|A(z,T)|^2$ がパワーになるように包絡線を規格化する。式 :eqref[eq-optics-gvd] に Kerr 位相の項 $i\gamma|A|^2A$ を加えると
$$
i\partial_z A-\frac{\beta_2}{2}\partial_T^2 A+\gamma|A|^2A=0.
\label{eq-optics-nlse}
$$
これが :term[非線形 Schrödinger 方程式]{en="nonlinear Schrödinger equation" yomi="ひせんけいしゅれーでぃんがーほうていしき"} である。損失、高次分散、応答の遅れは省いている。$\gamma>0,\beta_2<0$ とし $b=|\beta_2|$、形 $f$ を使って $A=\sqrt{P_0}\,f(T/T_0)e^{iKz}$ と置くと $-K f+(b/2T_0^2)f''+\gamma P_0f^3=0$。$f=\operatorname{sech}$ は $f''=f-2f^3$ を満たすため
$$
P_0=\frac{b}{\gamma T_0^2},\qquad K=\frac{b}{2T_0^2}.
\label{eq-optics-soliton}
$$
この解が基本の明るい :term[光ソリトン]{en="optical soliton" yomi="ひかりそりとん"} である。分散長 $L_D=T_0^2/b$ と非線形長 $L_{\rm NL}=1/(\gamma P_0)$ は等しい。強度積分は $\tanh'=\operatorname{sech}^2$ より $2P_0T_0$。逆に $\beta_2>0,\gamma>0$ では背景上に強度の穴を作る暗い解があり、$A=\sqrt{P_b}\tanh(T/T_b)e^{i\gamma P_b z}$、$T_b^2=\beta_2/(\gamma P_b)$ を入れれば確かめられる。

:::example[単位も合わせる]
$\beta_2=-20\,\mathrm{ps^2/km}$、$\gamma=2\,\mathrm{W^{-1}km^{-1}}$、$T_0=1\,\mathrm{ps}$ なら $P_0=10\,\mathrm W$、エネルギーは $20\,\mathrm{pJ}$。
:::
:::example[幅とエネルギー]
同じ媒質で幅を二倍にする基本解には、四分の一のピークパワーが必要である。エネルギーは幅とピークの積に比例するため半分になる。
:::

### 非線形伝播で保存する包絡線エネルギー
式 :eqref[eq-optics-nlse] に $A^*$ を掛け、共役の式と組み合わせると
$$
\partial_z |A|^2=\beta_2\,\partial_T\Im(A^*\partial_T A).
$$
三次項は絶対値が実なので打ち消し合う。$A$ と $\partial_T A$ が両端で減衰するパルスなら時間積分がゼロとなり、$\int |A|^2\dd T$ は保存される。暗ソリトンのように定数の背景を持つ解では、この積分はそのままでは発散する。強度欠損を積分するなら背景を引いて定義する。

### 幅と位相の一次の釣合いを確かめる
正の Kerr 係数でガウスパルスの後縁は周波数が高い側へずれる。一方、異常群速度分散 $\beta_2<0$ では、高い側の成分ほど単位長当たりの群遅延 $\beta_1+\beta_2\Omega$ が小さい。そのため後縁に生じた成分が中央へ近づく向きに動く。前縁では逆符号の周波数偏移と遅延が組になる。これだけでは任意の形のパルスが保たれる証明にはならず、形が変わらない厳密な条件は式 :eqref[eq-optics-soliton] の幅とピークパワーである。両効果の符号を同じ向きの時間座標で比較することが重要である。

:::example[分散長に合わせて出力を調整する]
$|\beta_2|=20\,\mathrm{ps^2/km},\gamma=2\,\mathrm{W^{-1}km^{-1}}$、$T_0=2\,\mathrm{ps}$ なら分散長は $0.20\,\mathrm{km}$。同じ長さを非線形長にするには $P_0=1/(2\cdot0.20)=2.5\,\mathrm W$ が要る。これを超えるパワーでは基本解の条件から外れる。
:::

### ソリトン次数を二つの長さから決める

異常分散 $\beta_2<0$ と自己集束型 $\gamma>0$、無損失を仮定する。初期包絡線を $\sqrt{P_0}\operatorname{sech}(T/T_0)$ とする。時間を $T_0$、距離を $L_D=T_0^2/|\beta_2|$、振幅を $\sqrt{P_0}$ で換算すると :eqref[eq-optics-nlse] の分散項の係数が $1/2$、三次項の係数が $N^2=\gamma P_0 L_D$ になる。$N=1$ なら $\operatorname{sech} x$ の二階微分が $\operatorname{sech}x-2\operatorname{sech}^3x$ なので、非線形項が形を変える部分を打ち消し、全体位相だけが進む。$N$ が一未満なら分散長の間に非線形位相が一に届かない。$N>1$ では一般に形が変わり、任意の初期パルスを $N$ だけで厳密な高次ソリトンと判定できない。

:::example[分散と非線形が釣り合う]
指定の二係数に対して基本解を作る。

入力は分散係数の絶対値 $\lvert\beta_2\rvert\simeq 2\times 10^{-26}\,\mathrm{s^2/m}$、非線形係数 $\gamma=\frac{1}{500}\,\mathrm{W^{-1}m^{-1}}$、双曲線正割の時間尺度 $T_0\simeq 1\times 10^{-12}\,\mathrm s$、入射ピークパワー $P_0=10\,\mathrm W$ である。

二次の分散がパルスの時間幅に効く距離を求める。

$$
\begin{aligned}
L_D&=\frac{T_0^{2}}{\lvert\beta_2\rvert},\\
L_D/(\mathrm m)&\simeq \frac{\left(\left(1\times 10^{-12}\right)\right)^{2}}{\left(2\times 10^{-26}\right)}=50.
\end{aligned}
$$

ピークで非線形位相が一ラジアンになる距離を求める。

$$
\begin{aligned}
L_{\rm NL}&=\frac{1}{P_0 \gamma},\\
L_{\rm NL}/(\mathrm m)&= \frac{1}{\left(10\right) \cdot \left(\frac{1}{500}\right)}=50.
\end{aligned}
$$

入射形を固定したときのソリトン次数のパラメータである。

$$
\begin{aligned}
N_s&=\sqrt{\frac{L_D}{L_{\rm NL}}},\\
N_s&= \sqrt{\frac{\left(50\right)}{\left(50\right)}}=1.
\end{aligned}
$$

二つの長さが一致するよう入力パワーを解く。

$$
\begin{aligned}
P_1&=\frac{1}{L_D \gamma},\\
P_1/(\mathrm W)&= \frac{1}{\left(50\right) \cdot \left(\frac{1}{500}\right)}=10.
\end{aligned}
$$

双曲線正割の二乗を時刻で積分する。

$$
\begin{aligned}
E_{\rm pulse}&=2 P_0 T_0,\\
E_{\rm pulse}/(\mathrm J)&\simeq 2 \cdot \left(10\right) \cdot \left(1\times 10^{-12}\right)\simeq 2\times 10^{-11}.
\end{aligned}
$$

二つの長さは等しい。この時間形なら一定のパワー分布を保つ。
:::

:::example[ピークが弱い入射]
同じ時間幅と伝播条件でパワーを四分の一にする。

入力は分散係数の絶対値 $\lvert\beta_2\rvert\simeq 2\times 10^{-26}\,\mathrm{s^2/m}$、非線形係数 $\gamma=\frac{1}{500}\,\mathrm{W^{-1}m^{-1}}$、双曲線正割の時間尺度 $T_0\simeq 1\times 10^{-12}\,\mathrm s$、入射ピークパワー $P_0=\frac{5}{2}\,\mathrm W$ である。

二次の分散がパルスの時間幅に効く距離を求める。

$$
\begin{aligned}
L_D&=\frac{T_0^{2}}{\lvert\beta_2\rvert},\\
L_D/(\mathrm m)&\simeq \frac{\left(\left(1\times 10^{-12}\right)\right)^{2}}{\left(2\times 10^{-26}\right)}=50.
\end{aligned}
$$

ピークで非線形位相が一ラジアンになる距離を求める。

$$
\begin{aligned}
L_{\rm NL}&=\frac{1}{P_0 \gamma},\\
L_{\rm NL}/(\mathrm m)&= \frac{1}{\left(\frac{5}{2}\right) \cdot \left(\frac{1}{500}\right)}=200.
\end{aligned}
$$

入射形を固定したときのソリトン次数のパラメータである。

$$
\begin{aligned}
N_s&=\sqrt{\frac{L_D}{L_{\rm NL}}},\\
N_s&= \sqrt{\frac{\left(50\right)}{\left(200\right)}}=\frac{1}{2}.
\end{aligned}
$$

二つの長さが一致するよう入力パワーを解く。

$$
\begin{aligned}
P_1&=\frac{1}{L_D \gamma},\\
P_1/(\mathrm W)&= \frac{1}{\left(50\right) \cdot \left(\frac{1}{500}\right)}=10.
\end{aligned}
$$

双曲線正割の二乗を時刻で積分する。

$$
\begin{aligned}
E_{\rm pulse}&=2 P_0 T_0,\\
E_{\rm pulse}/(\mathrm J)&\simeq 2 \cdot \left(\frac{5}{2}\right) \cdot \left(1\times 10^{-12}\right)\simeq 5\times 10^{-12}.
\end{aligned}
$$

分散長に対する非線形長が伸び、形を保つ基本解から外れる。
:::

:::example[時間幅を増す]
ピークを十ワットに戻し、初期の時間尺度を二倍にする。

入力は分散係数の絶対値 $\lvert\beta_2\rvert\simeq 2\times 10^{-26}\,\mathrm{s^2/m}$、非線形係数 $\gamma=\frac{1}{500}\,\mathrm{W^{-1}m^{-1}}$、双曲線正割の時間尺度 $T_0\simeq 2\times 10^{-12}\,\mathrm s$、入射ピークパワー $P_0=10\,\mathrm W$ である。

二次の分散がパルスの時間幅に効く距離を求める。

$$
\begin{aligned}
L_D&=\frac{T_0^{2}}{\lvert\beta_2\rvert},\\
L_D/(\mathrm m)&\simeq \frac{\left(\left(2\times 10^{-12}\right)\right)^{2}}{\left(2\times 10^{-26}\right)}=200.
\end{aligned}
$$

ピークで非線形位相が一ラジアンになる距離を求める。

$$
\begin{aligned}
L_{\rm NL}&=\frac{1}{P_0 \gamma},\\
L_{\rm NL}/(\mathrm m)&= \frac{1}{\left(10\right) \cdot \left(\frac{1}{500}\right)}=50.
\end{aligned}
$$

入射形を固定したときのソリトン次数のパラメータである。

$$
\begin{aligned}
N_s&=\sqrt{\frac{L_D}{L_{\rm NL}}},\\
N_s&= \sqrt{\frac{\left(200\right)}{\left(50\right)}}=2.
\end{aligned}
$$

二つの長さが一致するよう入力パワーを解く。

$$
\begin{aligned}
P_1&=\frac{1}{L_D \gamma},\\
P_1/(\mathrm W)&= \frac{1}{\left(200\right) \cdot \left(\frac{1}{500}\right)}=\frac{5}{2}.
\end{aligned}
$$

双曲線正割の二乗を時刻で積分する。

$$
\begin{aligned}
E_{\rm pulse}&=2 P_0 T_0,\\
E_{\rm pulse}/(\mathrm J)&\simeq 2 \cdot \left(10\right) \cdot \left(2\times 10^{-12}\right)\simeq 4\times 10^{-11}.
\end{aligned}
$$

分散長が四倍になり、同じピークでも非線形の相対効果は大きい。
:::

## 電気光学効果と音響光学効果

### 電圧で位相を変える
低周波電場 $E_k^{\rm dc}$ を加え、光周波数の逆比誘電率を $\Delta(\eps_r^{-1})_{ij}=\sum_k r_{ijk}E_k^{\rm dc}+\cdots$ と展開する。一次項を :term[Pockels 効果]{en="Pockels effect" yomi="ぽっけるすこうか"} と呼び、$r_{ijk}$ の単位は $\mathrm{m/V}$。反転対称のバルクでは静電場に奇な項が消える。主軸方向の $1/n^2$ を微分して $\Delta n=-n^3 r E^{\rm dc}/2$。長さ $L$、間隔 $d$ の電極で一様に $E^{\rm dc}=V/d$ が生じるとすると、片軸の位相が $\pi$ 動く電圧は $V_\pi=\lambda_0 d/(n^3|r|L)$ である。二腕差や偏光差を使う場合は、その差に対応する係数を用いる。

### 音波が作る動く格子
光弾性による変調が $e^{i(\vb K\cdot\vb r-\Omega t)}$ とその共役を持つとする。入射波に掛けると光は $\vb k\pm\vb K$、$\omega\pm\Omega$ へ散乱される。これが :term[音響光学効果]{en="acousto-optic effect" yomi="おんきょうこうがくこうか"} の周波数変換である。厚い格子で波数整合した一組を選ぶと、光の波数の大きさの変化を無視できる範囲では $|\vb k+\vb K|^2=|\vb k|^2$、すなわち $\vb k\cdot\vb K=-K^2/2$。音波面から測った角の大きさは $\sin\theta_B=\lambda/(2\Lambda)$ になる。$\lambda$ は媒質内の光波長、$\Lambda$ は音波波長である。

:::example[半波電圧]
$n=2,r=10\,\mathrm{pm/V},L=2\,\mathrm{cm},d=0.1\,\mathrm{mm},\lambda_0=1\,\mathrm{\mu m}$ なら、片軸の位相に対し $V_\pi=62.5\,\mathrm V$ である。
:::
:::example[音波に対する傾き]
$\Lambda=10\,\mathrm{\mu m},\lambda=0.5\,\mathrm{\mu m}$ なら $\sin\theta_B=0.025$、角は約 $1.43^\circ$。
:::

### 位相変調を強度変調へ換える
Mach–Zehnder の一腕に Pockels 素子を置き、位相を $\phi_0+\pi V/V_\pi$ とする。分割器と合波器が対称なら出力の一方は総パワーの $(1+\cos(\phi_0+\pi V/V_\pi))/2$。$\phi_0=\pi/2$ で微小電圧に対する一次変化は $-\pi V/(2V_\pi)$。一腕に位相を加えるだけではその腕のパワーは変わらないが、干渉させれば位相を明暗に変換できる。

### 音響格子の二波結合
厚い一様な格子で波数整合する入射と一次回折の二成分だけをとる。周波数差が光周波数に比べて十分小さい範囲で、光子束とパワーの規格化の差を無視する。格子の実屈折率変調を $\Delta n\cos(Kx-\Omega t)$ とすると、二つの指数成分がそれぞれ逆向きの結合を与える。光のスカラー波動方程式で屈折率の二乗の摂動は $2n\Delta n\cos(Kx-\Omega t)$。二波の各空間指数に投影すると余弦の対応する Fourier 係数 $1/2$ が掛かり、縦方向の波数を $k_z$ とすれば結合の大きさは $k_0^2n|\Delta n|/(2k_z)$。従って :eqref[eq-optics-coupled-modes] と同じ二波解により整合時の回折率は $\sin^2(\kappa L)$ になる。この二波近似が成り立つ範囲では、小さい角度ずれは結合式の対角差に入れられる。薄い格子などで三つ以上の次数が同程度に結合する場合は、それらも伝播式に含める。

:::example[電圧に対する微小変調]
$V_\pi=5\,\mathrm V$ の一腕変調器を直交点に置き、$V=0.10\,\mathrm V$ だけ動かす。総入力に対する出力パワーの変化は一次で $-\pi(0.10)/(10)\simeq-0.0314$。出力はもともと約 0.5 なので、背景出力に対する変化は約 $-6.28$ パーセントである。
:::

### 薄い音響格子では多数の回折次数が出る

光が相互作用中に横方向へ回折する効果を無視し、音響格子を薄い位相板として扱う。垂直入射の光は $e^{i\zeta\cos(Kx-\Omega t)}$ を掛けられ、位相の振幅は $\zeta=k_0\Delta nL$ である。Jacobi–Anger 展開 $e^{i\zeta\cos\theta}=\sum_{m\in\Z}i^mJ_m(\zeta)e^{im\theta}$ は、$e^{i\zeta e^{i\theta}/2}e^{i\zeta e^{-i\theta}/2}$ を Taylor 展開し、$e^{im\theta}$ の係数を集めて得る。従って $m$ 次に周波数 $\omega+m\Omega$ の成分が現れ、そのパワー比は近軸で $J_m(\zeta)^2$ となる。位相因子の絶対値は 1 なので、全次数のパワー比の和は Parseval の等式で 1 である。これを :term[Raman–Nath 回折]{en="Raman–Nath diffraction" yomi="らまんなすかいせき"} という。寄与する次数 $m$ に対し、格子を通る間に生じる横回折の位相 $m^2K^2L/(2k)$ が十分小さいときに薄い位相板として扱える。

### 二腕の位相変調を差動出力で測る

半波電圧 $V_\pi$ は差動の位相を $\pi$ 変える電圧とする。等強度の二腕を合波して片方から取り出す規格化パワーは $\tfrac12[1+\cos(\phi_b+\pi V/V_\pi)]$。もう一方には符号が逆の干渉項が現れ、両出力の和は一である。$\phi_b=\pi/2$ ならゼロ電圧の周囲でパワー変化が一次になる。導関数の絶対値はここで最大である。一方 $\phi_b=0$ なら一次の感度がなく、二次の応答が先頭に来る。半波電圧の値を比較するときは一腕の位相と差動位相のどちらで定義したかを明示する。

:::example[微小変調を直交位相で読む]
二つの出力がほぼ等しいバイアスで小電圧を加える。

入力はバイアス位相 $\phi_b\simeq 1.5708\,\mathrm{rad}$、変調電圧 $V=\frac{1}{20}\,\mathrm V$、差動半波電圧 $V_\pi=5\,\mathrm V$ である。

半波電圧を差の位相へ換算する。

$$
\begin{aligned}
\delta\phi&=\frac{\pi V}{V_\pi},\\
\delta\phi/(\mathrm{rad})&= \frac{\pi \cdot \left(\frac{1}{20}\right)}{\left(5\right)}\simeq 0.031416.
\end{aligned}
$$

バイアスと印加電圧の寄与を足す。

$$
\begin{aligned}
\phi&=\phi_b + \delta\phi,\\
\phi/(\mathrm{rad})&\simeq \left(1.570796327\right) + \left(0.03141592654\right)\simeq 1.6022.
\end{aligned}
$$

差の位相から合波のパワーを求める。

$$
\begin{aligned}
T_+&=\frac{\cos{\left(\phi \right)}}{2} + \frac{1}{2},\\
T_+&\simeq \frac{\cos{\left(\left(1.602212253\right) \right)}}{2} + \frac{1}{2}\simeq 0.48429.
\end{aligned}
$$

等分割で無損失としたときの相補出力を使う。

$$
\begin{aligned}
T_-&=1 - T_+,\\
T_-&\simeq 1 - \left(0.4842946205\right)\simeq 0.51571.
\end{aligned}
$$

電圧に対するパワーの一次感度を求める。

$$
\begin{aligned}
\dd T_+/\dd V&=- \frac{\pi \sin{\left(\phi \right)}}{2 V_\pi},\\
\dd T_+/\dd V/(\mathrm{V^{-1}})&\simeq - \frac{\pi \cdot \sin{\left(\left(1.602212253\right) \right)}}{2 \cdot \left(5\right)}\simeq -0.314.
\end{aligned}
$$

パワー差の変化が電圧の一次で現れる。
:::

:::example[最大透過の付近で駆動する]
駆動を同じにして位相バイアスをゼロにする。

入力はバイアス位相 $\phi_b=0\,\mathrm{rad}$、変調電圧 $V=\frac{1}{20}\,\mathrm V$、差動半波電圧 $V_\pi=5\,\mathrm V$ である。

半波電圧を差の位相へ換算する。

$$
\begin{aligned}
\delta\phi&=\frac{\pi V}{V_\pi},\\
\delta\phi/(\mathrm{rad})&= \frac{\pi \cdot \left(\frac{1}{20}\right)}{\left(5\right)}\simeq 0.031416.
\end{aligned}
$$

バイアスと印加電圧の寄与を足す。

$$
\begin{aligned}
\phi&=\phi_b + \delta\phi,\\
\phi/(\mathrm{rad})&\simeq \left(0\right) + \left(0.03141592654\right)\simeq 0.031416.
\end{aligned}
$$

差の位相から合波のパワーを求める。

$$
\begin{aligned}
T_+&=\frac{\cos{\left(\phi \right)}}{2} + \frac{1}{2},\\
T_+&\simeq \frac{\cos{\left(\left(0.03141592654\right) \right)}}{2} + \frac{1}{2}\simeq 0.99975.
\end{aligned}
$$

等分割で無損失としたときの相補出力を使う。

$$
\begin{aligned}
T_-&=1 - T_+,\\
T_-&\simeq 1 - \left(0.9997532802\right)\simeq 2.4672\times 10^{-4}.
\end{aligned}
$$

電圧に対するパワーの一次感度を求める。

$$
\begin{aligned}
\dd T_+/\dd V&=- \frac{\pi \sin{\left(\phi \right)}}{2 V_\pi},\\
\dd T_+/\dd V/(\mathrm{V^{-1}})&\simeq - \frac{\pi \cdot \sin{\left(\left(0.03141592654\right) \right)}}{2 \cdot \left(5\right)}\simeq -0.009868.
\end{aligned}
$$

第一出力の変化は小さく、ゼロ電圧の傾きは消える。
:::

:::example[半波電圧を印加する]
同相から逆相へ切り替える。

入力はバイアス位相 $\phi_b=0\,\mathrm{rad}$、変調電圧 $V=5\,\mathrm V$、差動半波電圧 $V_\pi=5\,\mathrm V$ である。

半波電圧を差の位相へ換算する。

$$
\begin{aligned}
\delta\phi&=\frac{\pi V}{V_\pi},\\
\delta\phi/(\mathrm{rad})&= \frac{\pi \cdot \left(5\right)}{\left(5\right)}\simeq 3.1416.
\end{aligned}
$$

バイアスと印加電圧の寄与を足す。

$$
\begin{aligned}
\phi&=\phi_b + \delta\phi,\\
\phi/(\mathrm{rad})&\simeq \left(0\right) + \left(3.141592654\right)\simeq 3.1416.
\end{aligned}
$$

差の位相から合波のパワーを求める。

$$
\begin{aligned}
T_+&=\frac{\cos{\left(\phi \right)}}{2} + \frac{1}{2},\\
T_+&\simeq \frac{\cos{\left(\left(3.141592654\right) \right)}}{2} + \frac{1}{2}=0.
\end{aligned}
$$

等分割で無損失としたときの相補出力を使う。

$$
\begin{aligned}
T_-&=1 - T_+,\\
T_-&= 1 - \left(0\right)=1.
\end{aligned}
$$

電圧に対するパワーの一次感度を求める。

$$
\begin{aligned}
\dd T_+/\dd V&=- \frac{\pi \sin{\left(\phi \right)}}{2 V_\pi},\\
\dd T_+/\dd V/(\mathrm{V^{-1}})&\simeq - \frac{\pi \cdot \sin{\left(\left(3.141592654\right) \right)}}{2 \cdot \left(5\right)}=0.
\end{aligned}
$$

同じ合波器でも、パワーの出る端子が入れ替わる。
:::

## 誘導 Raman 散乱と誘導 Brillouin 散乱

### 光のうなりで物質を駆動する
ポンプと Stokes 波の差周波数を $\Omega=\omega_p-\omega_s$、差波数を $\vb q=\vb k_p-\vb k_s$ とする。光から媒質へ渡る力が $A_p A_s^*$ に比例すると、媒質の変位の複素振幅は減衰振動子の式から
$$
Q_\Omega=\frac{C A_p A_s^*}{\Omega_m^2-\Omega^2-i\Gamma_m\Omega}
$$
となる。$\Omega_m$ は振動の固有角周波数、$\Gamma_m$ は減衰率。変位が光を変調すると、Stokes 成分の非線形分極は $Q_\Omega^* A_p$ に比例する。実結合かつ正の結合係数の積とすると、応答の正虚部に共役を取り、式 :eqref[eq-optics-nonlinear-envelope] の $i$ を掛けるので正の実利得が出る。ポンプ一定なら $\partial_z I_s=g_s I_p I_s$ であり、$g_s$ の単位は $\mathrm{m/W}$。共鳴付近では分母の二乗絶対値によって Lorentz 型の利得になる。

:term[誘導 Raman 散乱]{en="stimulated Raman scattering" yomi="ゆうどうらまんさんらん"} では分子や格子の内部振動が介在し、周波数差の尺度はその固有振動で決まる。:term[誘導 Brillouin 散乱]{en="stimulated Brillouin scattering" yomi="ゆうどうぶりるあんさんらん"} では音波が介在して $\Omega_m\simeq v_a|\vb q|$ となる。後方散乱では $|\vb q|\simeq2n k_0$ なので差周波数は $2nv_a/\lambda_0$。前方の同方向の二波では差波数は小さく、一般に同じ周波数ではない。反 Stokes 側は $\omega_p+\Omega$ で、媒質の振動からエネルギーを受け取る。

:::example[後方音響散乱]
$n=1.45,v_a=6.0\times10^3\,\mathrm{m/s},\lambda_0=1.55\,\mathrm{\mu m}$ なら後方 Brillouin の周波数差は約 $11.2\,\mathrm{GHz}$。
:::
:::example[一定ポンプの利得]
$g_s I_p L=2$ なら弱い種の強度は $e^2\simeq7.39$ 倍になる。種が厳密にゼロならこの決定論的模型の解はゼロであり、自発散乱を記述するには雑音源が必要である。
:::

### 共鳴の幅と入出力の単位
$\Omega$ が $\Omega_m$ の近くで $|\Omega-\Omega_m|\ll\Omega_m$ なら、分母は $-2\Omega_m[(\Omega-\Omega_m)+i\Gamma_m/2]$ と近似できる。応答の虚部は $(\Gamma_m/2)/[(\Omega-\Omega_m)^2+(\Gamma_m/2)^2]$ に比例し、角周波数の全半値幅が $\Gamma_m$ である。実験の通常周波数での幅は $\Gamma_m/(2\pi)$。Brillouin の後方配置では Stokes 波はポンプと逆向きに進むため、小信号利得の積分座標は Stokes の進行方向に取り直す。ポンプ減衰を含めると、同じ物理座標で二波の伝播式の増減符号は逆向きに現れる。

:::example[線幅を二種類の周波数単位で換算する]
媒質の振動の減衰率が $\Gamma_m=2\pi(30\times10^6)\,\mathrm{s^{-1}}$ なら、Stokes 利得の角周波数全半値幅から Hz に換算した値は $30\,\mathrm{MHz}$ である。ポンプと Stokes の周波数差は固有周波数付近の GHz などの値であり、この狭い利得幅とは別である。
:::

### ポンプと散乱光のエネルギーの行き先

振動の振幅が定常に近い場合、同じ向きの光二波について、光子束の交換量は $\partial_z(I_s/\hbar\omega_s)=-\partial_z(I_p/\hbar\omega_p)$ である。実際、光のうなりに含まれる積 $A_p A_s^*$ はポンプ光子を一つ減らすと同時に Stokes 光子を一つ増やす組に対応する。光の強度で表すと $\partial_z I_p=-(\omega_p/\omega_s)\partial_z I_s$ であり、差に当たる $\hbar\Omega$ は一回の遷移ごとに媒質へ渡る。

一方、反 Stokes 波だけを種とし、そのうなりで同じ型の受動的な振動を駆動すると、振動を駆動する積は $A_{\rm as} A_p^*$、光を変調する積は $Q_\Omega A_p$ となる。分極に応答の共役が現れないので $i$ を掛けた項の実部は負であり、この二波模型では反 Stokes 波は減衰する。媒質の振動が別の光や外部の力で駆動されていれば、その振動から反 Stokes 光へのエネルギー移動も起こり得る。

### 減衰ポンプと離調が散乱利得を下げる

共鳴での小信号係数を $g_0$ とし、角周波数のずれ $\delta$ で $g=g_0/[1+4\delta^2/\Gamma_m^2]$ とする。ポンプは Stokes 波の枯渇を無視して $I_p(z)=I_{p0}e^{-\alpha_p z}$。Stokes とポンプが同じ向きに進む場合、Stokes の独立の損失係数を $\alpha_s$ とすると $\partial_z\ln I_s=gI_p-\alpha_s$。両辺を積分すると正味の対数利得は $g I_{p0}(1-e^{-\alpha_pL})/\alpha_p-\alpha_s L$。離調の絶対値が全半値幅の半分なら小信号の局所利得は二分の一になり、同じ損失に対して正味増幅が維持できない場合もある。

:::example[共鳴で種を増幅する]
小さな Stokes の種をポンプと同方向へ通す振動応答模型を評価する。

入力は共鳴時の小信号係数 $g_0\simeq 1\times 10^{-11}\,\mathrm{m/W}$、入口のポンプ照度 $I_{p0}\simeq 1\times 10^{10}\,\mathrm{W/m^2}$、線幅で測った角離調 $\delta/\Gamma_m=0$、ポンプ損失係数 $\alpha_p=\frac{1}{10}\,\mathrm{m^{-1}}$、Stokes 損失係数 $\alpha_s=\frac{1}{25}\,\mathrm{m^{-1}}$、媒質長 $L=5\,\mathrm m$ である。

減衰振動子の線形状を角離調で評価する。

$$
\begin{aligned}
r&=\frac{1}{1+4(\delta/\Gamma_m)^2},\\
r&= \frac{1}{4 \cdot \left(0\right)^{2} + 1}=1.
\end{aligned}
$$

減衰するポンプの照度を伝播区間で積分する。

$$
\begin{aligned}
L_{\rm eff}&=\frac{1 - e^{- L \alpha_p}}{\alpha_p},\\
L_{\rm eff}/(\mathrm m)&= \frac{1 - e^{- \left(5\right) \cdot \left(\frac{1}{10}\right)}}{\left(\frac{1}{10}\right)}\simeq 3.9347.
\end{aligned}
$$

線幅の補正と有効長を小信号利得へ掛ける。

$$
\begin{aligned}
A&=L_{\rm eff} g_0 I_{p0} r,\\
A&\simeq \left(3.934693403\right) \cdot \left(1\times 10^{-11}\right) \cdot \left(1\times 10^{10}\right) \cdot \left(1\right)\simeq 0.39347.
\end{aligned}
$$

Stokes が独立に受ける損失を引く。

$$
\begin{aligned}
G_{\rm net}&=A - L \alpha_s,\\
G_{\rm net}&\simeq \left(0.3934693403\right) - \left(5\right) \cdot \left(\frac{1}{25}\right)\simeq 0.19347.
\end{aligned}
$$

対数からパワーの比へ戻す。

$$
\begin{aligned}
R&=e^{G_{\rm net}},\\
R&\simeq e^{\left(0.1934693403\right)}\simeq 1.2135.
\end{aligned}
$$

共鳴の誘導利得が損失を上回る。
:::

:::example[線の半値へ移す]
ポンプ照度を固定し、二光波の周波数差だけを変える。

入力は共鳴時の小信号係数 $g_0\simeq 1\times 10^{-11}\,\mathrm{m/W}$、入口のポンプ照度 $I_{p0}\simeq 1\times 10^{10}\,\mathrm{W/m^2}$、線幅で測った角離調 $\delta/\Gamma_m=\frac{1}{2}$、ポンプ損失係数 $\alpha_p=\frac{1}{10}\,\mathrm{m^{-1}}$、Stokes 損失係数 $\alpha_s=\frac{1}{25}\,\mathrm{m^{-1}}$、媒質長 $L=5\,\mathrm m$ である。

減衰振動子の線形状を角離調で評価する。

$$
\begin{aligned}
r&=\frac{1}{1+4(\delta/\Gamma_m)^2},\\
r&= \frac{1}{4 \cdot \left(\frac{1}{2}\right)^{2} + 1}=\frac{1}{2}.
\end{aligned}
$$

減衰するポンプの照度を伝播区間で積分する。

$$
\begin{aligned}
L_{\rm eff}&=\frac{1 - e^{- L \alpha_p}}{\alpha_p},\\
L_{\rm eff}/(\mathrm m)&= \frac{1 - e^{- \left(5\right) \cdot \left(\frac{1}{10}\right)}}{\left(\frac{1}{10}\right)}\simeq 3.9347.
\end{aligned}
$$

線幅の補正と有効長を小信号利得へ掛ける。

$$
\begin{aligned}
A&=L_{\rm eff} g_0 I_{p0} r,\\
A&\simeq \left(3.934693403\right) \cdot \left(1\times 10^{-11}\right) \cdot \left(1\times 10^{10}\right) \cdot \left(\frac{1}{2}\right)\simeq 0.19673.
\end{aligned}
$$

Stokes が独立に受ける損失を引く。

$$
\begin{aligned}
G_{\rm net}&=A - L \alpha_s,\\
G_{\rm net}&\simeq \left(0.1967346701\right) - \left(5\right) \cdot \left(\frac{1}{25}\right)\simeq -0.0032653.
\end{aligned}
$$

対数からパワーの比へ戻す。

$$
\begin{aligned}
R&=e^{G_{\rm net}},\\
R&\simeq e^{\left(-0.003265329856\right)}\simeq 0.99674.
\end{aligned}
$$

媒質の利得は半分になり、損失との収支がわずかに吸収側へ移る。
:::

:::example[共鳴から遠ざける]
同じ長さと光量で線幅一つ分離す。

入力は共鳴時の小信号係数 $g_0\simeq 1\times 10^{-11}\,\mathrm{m/W}$、入口のポンプ照度 $I_{p0}\simeq 1\times 10^{10}\,\mathrm{W/m^2}$、線幅で測った角離調 $\delta/\Gamma_m=1$、ポンプ損失係数 $\alpha_p=\frac{1}{10}\,\mathrm{m^{-1}}$、Stokes 損失係数 $\alpha_s=\frac{1}{25}\,\mathrm{m^{-1}}$、媒質長 $L=5\,\mathrm m$ である。

減衰振動子の線形状を角離調で評価する。

$$
\begin{aligned}
r&=\frac{1}{1+4(\delta/\Gamma_m)^2},\\
r&= \frac{1}{4 \cdot \left(1\right)^{2} + 1}=\frac{1}{5}.
\end{aligned}
$$

減衰するポンプの照度を伝播区間で積分する。

$$
\begin{aligned}
L_{\rm eff}&=\frac{1 - e^{- L \alpha_p}}{\alpha_p},\\
L_{\rm eff}/(\mathrm m)&= \frac{1 - e^{- \left(5\right) \cdot \left(\frac{1}{10}\right)}}{\left(\frac{1}{10}\right)}\simeq 3.9347.
\end{aligned}
$$

線幅の補正と有効長を小信号利得へ掛ける。

$$
\begin{aligned}
A&=L_{\rm eff} g_0 I_{p0} r,\\
A&\simeq \left(3.934693403\right) \cdot \left(1\times 10^{-11}\right) \cdot \left(1\times 10^{10}\right) \cdot \left(\frac{1}{5}\right)\simeq 0.078694.
\end{aligned}
$$

Stokes が独立に受ける損失を引く。

$$
\begin{aligned}
G_{\rm net}&=A - L \alpha_s,\\
G_{\rm net}&\simeq \left(0.07869386806\right) - \left(5\right) \cdot \left(\frac{1}{25}\right)\simeq -0.12131.
\end{aligned}
$$

対数からパワーの比へ戻す。

$$
\begin{aligned}
R&=e^{G_{\rm net}},\\
R&\simeq e^{\left(-0.1213061319\right)}\simeq 0.88576.
\end{aligned}
$$

誘導効果だけなら正でも、透過全体では損失が勝る。
:::

## 演習問題

::::exercise{level=basic}
反転対称の媒質でバルクの電気双極子による $\chi^{(2)}$ はどうなるか。

:::solution
分極は電場を反転すると符号が変わる。二次の単項式は符号が変わらないため、Taylor 係数を比べるとゼロとなる。
:::
::::

::::exercise{level=basic}
$E=E_0\cos\omega t$ の二乗の周波数成分を求めよ。

:::solution
$\cos^2\omega t=(1+\cos2\omega t)/2$ より、直流と二倍周波数を含む。両者の実表示の係数は $E_0^2/2$。
:::
::::

::::exercise{level=basic}
位相整合した小変換の第 2 高調波で結晶を半分の長さにする。生成強度を比較せよ。

:::solution
:eqref[eq-optics-shg] で振幅が長さに比例する。従って振幅は半分、同一媒質での強度は四分の一となる。
:::
::::

::::exercise{level=basic}
ポンプ $300\,\mathrm{THz}$、信号 $180\,\mathrm{THz}$ の差周波を求めよ。

:::solution
周波数の差 $300-180=120\,\mathrm{THz}$ が正の出力周波数となる。
:::
::::

::::exercise{level=advanced}
$\beta_2=-18\,\mathrm{ps^2/km},\gamma=2\,\mathrm{W^{-1}km^{-1}},T_0=3\,\mathrm{ps}$ の基本明ソリトンのピークを求めよ。

:::solution
:eqref[eq-optics-soliton] を使うと $P_0=18/(2\cdot9)=1\,\mathrm W$ となる。
:::
::::

::::exercise{level=advanced}
三波混合で高周波光子束が減るとき、低周波の二つのパワー増加の比を求めよ。

:::solution
Manley–Rowe の関係により二波の光子束の増分は等しい。各波の一光子エネルギーを掛けると、パワー増分の比は $\omega_1:\omega_2$ となる。
:::
::::

::::exercise{level=advanced}
$P(T)=P_0e^{-T^2/T_0^2}$ として、長さ $L$ の Kerr 位相が作る瞬時周波数偏移を求めよ。

:::solution
位相は $\gamma LP(T)$ である。負の時間微分は $2\gamma L T P(T)/T_0^2$ になる。
:::
::::

::::exercise{level=advanced}
同じ屈折率と音速で後方 Brillouin の真空波長を二倍にする。周波数差を比較せよ。

:::solution
差周波数は $2nv_a/\lambda_0$ で、屈折率と音速を一定に保つため、周波数差は半分になる。
:::
::::

::::exercise{level=advanced}
二値反転した二次非線形結晶で、一次の擬似位相整合を用いる。変換が小さいとし、設計値からの残留波数差を $\delta k=\pi/L$ とする。同じ長さの完全位相整合した一様結晶に対する第 2 高調波の強度比を求めよ。

:::solution
反転による一次の有効係数は一様結晶の $2/\pi$ 倍である。さらに残留波数差のため振幅に $\operatorname{sinc}(\delta kL/2)=\operatorname{sinc}(\pi/2)=2/\pi$ が掛かる。従って全体の振幅に掛かる因子は $4/\pi^2$、強度比はこれを二乗した $16/\pi^4$ となる。
:::
::::

## 確認テスト

::::quiz
第 2 高調波を位相整合したとき何が起きるか。

- [x] 各位置で作る振幅が同相で足される。
- [ ] 線形吸収が消える。
- [ ] 基本波の強度への依存が消える。

:::explanation
生成項の因子 $e^{i\Delta kz}$ が 1 になり、場所による位相のずれが作る打ち消しがなくなる。
:::
::::

::::quiz{answer="1" tol="1e-6"}
$\gamma=0.01\,\mathrm{W^{-1}m^{-1}}$、$P=2\,\mathrm W$、$L=50\,\mathrm m$ の Kerr 位相を rad 単位で答えよ。

:::explanation
局所の位相変化は距離当たり $\gamma P$。値が一様なので $\gamma PL=0.01(2)(50)=1$ となる。
:::
::::

::::quiz
$\gamma>0$ で基本の明ソリトンを作る群速度分散の符号を選べ。

- [x] $\beta_2<0$
- [ ] $\beta_2>0$
- [ ] $\beta_2=0$

:::explanation
$\operatorname{sech}$ を :eqref[eq-optics-nlse] に入れると、分散の三次項と Kerr の三次項が相殺されなければならない。$\gamma>0$ では $\beta_2<0$ が必要となる。
:::
::::

::::quiz{answer="532" tol="1e-6"}
$1064\,\mathrm{nm}$ の真空光の第 2 高調波の波長を nm 単位で答えよ。

:::explanation
周波数は二倍、両者の真空位相速度は $c$ なので波長は半分である。
:::
::::

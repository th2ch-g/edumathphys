---
title: 反射と屈折
description: 平面境界での電磁場を解き、Fresnel 係数、全反射、ビームの変位および金属の反射を導く。
notation:
  - symbol: '\theta_i,\theta_t'
    meaning: 入射角および透過角
  - symbol: 'r_s,r_p'
    meaning: 電場の振幅反射係数
  - symbol: 'k_x'
    meaning: 境界に平行な波数
  - symbol: '\kappa'
    meaning: 境界法線方向の減衰定数
---

## この章の目標

- 接線方向の位相整合と境界条件を区別する。
- 偏光ごとの振幅とパワーを計算する。
- 全反射の近接場と反射ビームの変位を関係づける。
- 吸収を含む媒質で波数の枝を選ぶ。

## 境界条件と反射・屈折の法則

### 境界上で位相がそろう条件
境界を $z=0$、入射側を $z<0$ とする。自由な表面電荷と表面電流のない静止境界では、接線 $\vb E,\vb H$ と法線 $\vb D,\vb B$ が連続である。法線の条件は Maxwell の発散式を薄い円柱で積分し、高さをゼロへ近づけると得られる。接線の条件は境界をまたぐ長方形に回転の二式を積分して得られる。面積がゼロへ近づくと、有限な $\partial_t\vb B,\partial_t\vb D$ の面積積分が消え、境界に沿う二辺の周回積分だけが残る。

入射面を $xz$ 面とし、以下では媒質を線形・等方とし、界面の両側でそれぞれ一様であるとする。反射・屈折角の議論では正の実誘電率・透磁率を仮定する。各場の境界での位相は $k_x x-\omega t$ である。異なる $k_x,\omega$ の指数関数は独立なので、任意の $x,t$ で三つの波の境界条件を満たすには、反射と透過で入射と同じ周波数と接線波数を使う。例えば境界に平行な方向に Fourier 変換すると、各接線波数の係数ごとに独立の境界条件を得る。角を法線から測ると
$$
\theta_r=\theta_i,\qquad n_1\sin\theta_i=n_2\sin\theta_t .
\label{eq-optics-snell}
$$
第 2 式を :term[屈折の法則]{en="law of refraction" yomi="くっせつのほうそく"} と呼ぶ。正の実屈折率では $k_{1z}=k_1\cos\theta_i$、反射波の法線波数は $-k_{1z}$ である。透過側には外向きの波だけを採る。

:::example[空気からガラス]
$n_1=1,n_2=1.5,\theta_i=30^\circ$ なら $\sin\theta_t=1/3$ である。位相の平行成分は変わらないが、法線成分は $k_0\sqrt{(3/2)^2-(1/2)^2}=\sqrt2\,k_0$ に増す。
:::

:::example[逆進]
透過した方向を逆にたどる波を媒質 2 から入れると、共通の接線波数の絶対値は同じである。式 :eqref[eq-optics-snell] で 1 と 2 を交換すれば透過波は元の入射方向を逆進する。
:::

### 表面電流を残すとき
界面に自由表面電荷密度 $\sigma_f$、自由表面電流密度 $\vb K_f$ があると、薄い円柱に積分した電荷は $\sigma_f$ と底面積の積で残り、周回の積分でも $\vb K_f$ の項が残る。媒質 1 から 2 への法線を $\hat{\vb n}$ とすると
$$
\hat{\vb n}\cdot(\vb D_2-\vb D_1)=\sigma_f,\qquad
\hat{\vb n}\times(\vb H_2-\vb H_1)=\vb K_f .
$$
磁気単極子や磁気表面電流を置いていないので、$\vb B$ の法線と $\vb E$ の接線は依然連続である。等方な薄い導電シートなら $\vb K_f=G_s\vb E_\parallel$ と置ける。$G_s$ の単位は siemens である。

:::example[シートだけで反射を作る]
両側が真空の静止シートへ垂直入射し、$Y=1/\eta_0$ とする。境界式は $t=1+r$、$Y(1-r)-Yt=G_s t$ なので $r=-G_s/(2Y+G_s)$。実の $G_s=2Y$ なら $r=-1/2,t=1/2$ である。パワー反射と透過は各 $1/4$、残る $1/2$ はシートで仕事をする。
:::

### 境界条件を積分形から導く

界面を横切る高さ $2h$ の直方体に $\divg\vb D=\rho_{\rm f}$ を積分する。面積を $A$ とすれば側面の流束は場が有界なら $O(h)$ であり、蓋と底から $A\hat{\vb z}\cdot(\vb D_2-\vb D_1)$ が残る。自由表面電荷密度を $\sigma_f$ とすると $\hat{\vb z}\cdot(\vb D_2-\vb D_1)=\sigma_f$。磁場の Gauss 則にも同じ操作をすれば $\hat{\vb z}\cdot(\vb B_2-\vb B_1)=0$ になる。次に接線 $\hat{\vb t}$ に沿う辺と法線辺で作った細長い長方形を用いる。Faraday 則の面を通る磁束変化は $O(h)$ で消えるため、電場の接線成分は連続する。Ampère–Maxwell 則では同様に変位電流が消える一方、表面電流 $\vb K_f$ の寄与は有限であり、周回の向きをそろえて $\hat{\vb z}\times(\vb H_2-\vb H_1)=\vb K_f$ が残る。

界面の接線方向と時間に並進対称性があるので、境界方程式は Fourier 成分 $(k_x,k_y,\omega)$ ごとに独立である。二種類の入射があっても周波数や接線波数の異なる成分を混ぜて条件を解く必要はない。例えば面に沿って周期 $d$ の変調を加えると並進対称性は離散になり、同じ位相を与える $k_x+2\pi m/d$ が透過波に含まれ得る。この比較により、平面という形状に加えて面内の一様性も Snell の法則の仮定であることが分かる。

:::example[空気側から水面を見る]
水面を局所的に平らとみなし、斜めに単色波を入れる。

入力は入射側の屈折率 $n_{1}=1$、透過側の屈折率 $n_{2}=\frac{4}{3}$、真空波長 $\lambda_0\simeq 6\times 10^{-7}\,\mathrm m$、入射角 $\theta_i\simeq 0.5236\,\mathrm{rad}$ である。

空間位相を波長から求める。

$$
\begin{aligned}
k_0&=\frac{2 \pi}{\lambda_0},\\
k_0/(\mathrm{m^{-1}})&\simeq \frac{2 \cdot \pi}{\left(6\times 10^{-7}\right)}\simeq 1.0472\times 10^{7}.
\end{aligned}
$$

界面上での位相が揃うように、入射波の接線波数を透過側にも使う。

$$
\begin{aligned}
K&=k_0 n_{1} \sin{\left(\theta_i \right)},\\
K/(\mathrm{m^{-1}})&\simeq \left(1.047197551\times 10^{7}\right) \cdot \left(1\right) \cdot \sin{\left(\left(0.5235987756\right) \right)}\simeq 5.236\times 10^{6}.
\end{aligned}
$$

媒質内の分散関係から接線波数の寄与を引く。

$$
\begin{aligned}
k_{2z}^2&=- K^{2} + k_0^{2} n_{2}^{2},\\
k_{2z}^2/(\mathrm{m^{-2}})&\simeq - \left(\left(5.235987756\times 10^{6}\right)\right)^{2} + \left(\left(1.047197551\times 10^{7}\right)\right)^{2} \cdot \left(\frac{4}{3}\right)^{2}\simeq 1.6754\times 10^{14}.
\end{aligned}
$$

正の実数または正の虚部を持つ枝を選び、半空間の外向き条件を満たす。

$$
\begin{aligned}
k_{2z}&=\sqrt{k_{2z}^2},\\
k_{2z}/(\mathrm{m^{-1}})&\simeq \sqrt{\left(1.675395809\times 10^{14}\right)}\simeq 1.2944\times 10^{7}.
\end{aligned}
$$

界面に沿って測った位相周期を求める。両側とも同じ値である。

$$
\begin{aligned}
\Lambda_x&=\frac{2 \pi}{K},\\
\Lambda_x/(\mathrm m)&\simeq \frac{2 \cdot \pi}{\left(5.235987756\times 10^{6}\right)}\simeq 1.2\times 10^{-6}.
\end{aligned}
$$

透過側の法線波数は正であり、位相が水中へ進む。水中での総波数が増しても、界面上の周期は保存する。
:::

:::example[平行に進む境界]
入射角を選び、透過波の波数が面に沿う条件を作る。

入力は入射側の屈折率 $n_{1}=2$、透過側の屈折率 $n_{2}=1$、真空波長 $\lambda_0\simeq 6\times 10^{-7}\,\mathrm m$、入射角 $\theta_i=\pi/6\,\mathrm{rad}$ である。

空間位相を波長から求める。

$$
\begin{aligned}
k_0&=\frac{2 \pi}{\lambda_0},\\
k_0/(\mathrm{m^{-1}})&\simeq \frac{2 \cdot \pi}{\left(6\times 10^{-7}\right)}\simeq 1.0472\times 10^{7}.
\end{aligned}
$$

界面上での位相が揃うように、入射波の接線波数を透過側にも使う。

$$
\begin{aligned}
K&=k_0 n_{1} \sin{\left(\theta_i \right)},\\
K&=k_0\cdot2\sin(\pi/6)=k_0.
\end{aligned}
$$

媒質内の分散関係から接線波数の寄与を引く。

$$
\begin{aligned}
k_{2z}^2&=- K^{2} + k_0^{2} n_{2}^{2},\\
k_{2z}^2&=-k_0^2+k_0^2=0.
\end{aligned}
$$

臨界角では、伝播側から外向きの根を連続に取ると $k_{2z}=0$ になる。

$$
\begin{aligned}
k_{2z}&=\sqrt{k_{2z}^2},\\
k_{2z}/(\mathrm{m^{-1}})&= \sqrt{\left(0\right)}=0.
\end{aligned}
$$

界面に沿って測った位相周期を求める。両側とも同じ値である。

$$
\begin{aligned}
\Lambda_x&=\frac{2 \pi}{K},\\
\Lambda_x/(\mathrm m)&\simeq \frac{2 \cdot \pi}{\left(1.047197551\times 10^{7}\right)}\simeq 6\times 10^{-7}.
\end{aligned}
$$

法線波数はゼロである。透過側の位相は法線方向には進まない。
:::

:::example[低屈折率側で減衰する条件]
同じ二媒質でさらに角度を増やす。

入力は入射側の屈折率 $n_{1}=2$、透過側の屈折率 $n_{2}=1$、真空波長 $\lambda_0\simeq 6\times 10^{-7}\,\mathrm m$、入射角 $\theta_i\simeq 0.7854\,\mathrm{rad}$ である。

空間位相を波長から求める。

$$
\begin{aligned}
k_0&=\frac{2 \pi}{\lambda_0},\\
k_0/(\mathrm{m^{-1}})&\simeq \frac{2 \cdot \pi}{\left(6\times 10^{-7}\right)}\simeq 1.0472\times 10^{7}.
\end{aligned}
$$

界面上での位相が揃うように、入射波の接線波数を透過側にも使う。

$$
\begin{aligned}
K&=k_0 n_{1} \sin{\left(\theta_i \right)},\\
K/(\mathrm{m^{-1}})&\simeq \left(1.047197551\times 10^{7}\right) \cdot \left(2\right) \cdot \sin{\left(\left(0.7853981634\right) \right)}\simeq 1.481\times 10^{7}.
\end{aligned}
$$

媒質内の分散関係から接線波数の寄与を引く。

$$
\begin{aligned}
k_{2z}^2&=- K^{2} + k_0^{2} n_{2}^{2},\\
k_{2z}^2/(\mathrm{m^{-2}})&\simeq - \left(\left(1.480960979\times 10^{7}\right)\right)^{2} + \left(\left(1.047197551\times 10^{7}\right)\right)^{2} \cdot \left(1\right)^{2}\simeq -1.0966\times 10^{14}.
\end{aligned}
$$

正の実数または正の虚部を持つ枝を選び、半空間の外向き条件を満たす。

$$
\begin{aligned}
k_{2z}&=\sqrt{k_{2z}^2},\\
k_{2z}/(\mathrm{m^{-1}})&\simeq \sqrt{\left(-1.096622711\times 10^{14}\right)}\simeq 1.0472\times 10^{7}i.
\end{aligned}
$$

界面に沿って測った位相周期を求める。両側とも同じ値である。

$$
\begin{aligned}
\Lambda_x&=\frac{2 \pi}{K},\\
\Lambda_x/(\mathrm m)&\simeq \frac{2 \cdot \pi}{\left(1.480960979\times 10^{7}\right)}\simeq 4.2426\times 10^{-7}.
\end{aligned}
$$

法線波数が純虚数になった。接線方向の位相は保存され、法線方向には振幅が減る。
:::

## Fresnel の式

### 接線電場で振幅をそろえる
電場が入射面に垂直なものを :term[s 偏光]{en="s polarization" yomi="えすへんこう"}、電場が入射面に含まれるものを :term[p 偏光]{en="p polarization" yomi="ぴーへんこう"} と呼ぶ。斜入射の振幅には基底の符号が入るため、ここでは両偏光とも接線電場を振幅とする。$s$ では $E_y,H_x$ を取り、前向き波で $H_x=-Y_s E_y$ と書く。$p$ では $E_x,H_y$ を取り、前向き波で $H_y=Y_p E_x$ と書く。式 :eqref[eq-optics-plane] より
$$
Y_s=\frac{k_z}{\omega\mu},\qquad
Y_p=\frac{\omega\eps}{k_z}.
$$
$p$ の二式目は Ampère–Maxwell の式の $x$ 成分からも得られる。反射波では両方とも $Y$ の符号を逆にする。入射振幅を 1 とすると、接線場の連続条件は $1+r=t$ と $Y_1(1-r)=Y_2t$ である。これを解いたものが :term[Fresnel 係数]{en="Fresnel coefficient" yomi="ふれねるけいすう"} である。
$$
r_j=\frac{Y_{1j}-Y_{2j}}{Y_{1j}+Y_{2j}},\qquad
t_j=\frac{2Y_{1j}}{Y_{1j}+Y_{2j}},\qquad j=s,p.
\label{eq-optics-fresnel}
$$
全電場の振幅で $p$ を扱う場合、透過係数には $\cos\theta_i/\cos\theta_t$ が掛かる。反射側の $p$ 基底を接線が負になるように取れば、反射係数には負号も付く。

### パワー係数と保存則
接線振幅 $E_t$ の波が法線方向に運ぶパワー密度は $\Re Y\,|E_t|^2/2$ である。入射側が無損失で伝播波なら $R=|r|^2$、$T=\Re Y_2\,|t|^2/Y_1$ となる。実 $Y_1>0$ に対し
$$
|Y_1+Y_2|^2-|Y_1-Y_2|^2=4Y_1\Re Y_2
$$
であるから $R+T=1$ を得る。媒質 2 が吸収性なら $T$ は境界を横切る割合であり、遠くまで残る割合は伝播損失も含めて求める。

:::example[垂直入射]
非磁性媒質では垂直入射で $Y=n/\eta_0$ である。$n_1=1,n_2=1.5$ なら $r=-1/5,t=4/5$ であり、$R=1/25$、$T=24/25$ である。
:::

:::example[磁気応答を含む整合]
$\eps_r=\mu_r=2$ の媒質に真空から垂直入射すると $n=2$ だが、$\eta=\eta_0$ なので $r=0$ である。屈折率が変わっても垂直入射の反射はインピーダンスで決まる。
:::

### 横磁場を使って符号を検算する
$z=0$ で $p$ 波を考え、入射の $E_x$ と $H_y$ の符号を両方正に取る。反射で $k_z$ が反転すると、$H_y/E_x$ が反転する。透過は外向きなので無損失なら実の $Y_{2p}>0$ である。反射側で進行方向に付随した局所基底を使うときは、その基底で $E_x$ が変わる符号まで反射係数に含める。円偏光や楕円偏光の反射を計算する際は、Fresnel の表の係数だけでなく入出射の基底が同じ約束で与えられているかを確認する。

:::example[Snell 則を直接代入する]
$n_1=1,n_2=\sqrt2$ に $\theta_i=45^\circ$ で入射すると $\sin\theta_t=1/2$。$\eta_0$ を掛けた $s$ の $Y_1,Y_2$ は $1/\sqrt2,\sqrt6/2$、$p$ は $\sqrt2,2\sqrt{2/3}$ である。従って $r_s=(1-\sqrt3)/(1+\sqrt3)$、$r_p=(\sqrt3-2)/(\sqrt3+2)$ となる。
:::

### パワーを同じ尺度で並べる
両側が無損失、伝播波とし、左右からの入射を同時に扱う。各進行波の接線電場に $\sqrt{Y_j}$ を掛けて振幅を規格化すると、その絶対値二乗は法線パワーの共通の尺度となる。左からの入射では反射 $r=(Y_1-Y_2)/(Y_1+Y_2)$、透過 $t_P=2\sqrt{Y_1Y_2}/(Y_1+Y_2)$。逆向きでは反射は $-r$、透過は同じ $t_P$ である。従って入出力の散乱行列は
$$
S=\begin{pmatrix}r&t_P\\ t_P&-r\end{pmatrix}.
$$
$r^2+t_P^2=1$ と非対角の打ち消しから $S^\mathsf T S=I$。両側から同時に光を入れると各出力に干渉項があっても、二出力の総パワーは保存する。この規格化で透過が対称になることは、片側から見た電場の $t$ が入射方向により変わることと矛盾しない。

### 接線成分と偏光単位ベクトルの換算

入射面を $xz$、$s$ 方向を $\hat{\vb y}$ にそろえ、各波の $p$ 単位ベクトルを $\hat{\vb p}=\hat{\vb y}\times\hat{\vb k}$ と決める。これは入射波には $\cos\theta_i\,\hat{\vb x}-\sin\theta_i\,\hat{\vb z}$、反射波には $-\cos\theta_i\,\hat{\vb x}-\sin\theta_i\,\hat{\vb z}$ を与える。従って、この単位ベクトルでの反射振幅は接線電場を比べる :eqref[eq-optics-fresnel] の $r_p$ の負である。透過側の $x$ 成分は $\cos\theta_t$ なので、単位ベクトル規約の透過振幅は接線の $t_p$ に $\cos\theta_i/\cos\theta_t$ を掛けたものになる。両規約でパワーは同じである。実際、無損失で透過の流束は単位ベクトル振幅の二乗に $\cos\theta_t/(2\eta_2)$ を掛けて得られ、接線振幅では同じ式が $|E_x|^2/(2\eta_2\cos\theta_t)$ になる。

一般の非磁性・無損失界面で入射面への分解が重要な理由は、正射影が異なる二成分を別のアドミタンスに送るからである。入射電場を $\cos\psi\,\hat{\vb s}+\sin\psi\,\hat{\vb p}$ と規格化すれば、反射率は $|r_s|^2\cos^2\psi+|r_p|^2\sin^2\psi$。反射パワーの式に交差項がないのは、同方向へ戻る直交偏光の Poynting ベクトルが横振幅のノルムの二乗に比例するためである。検光子で一つの偏光へ射影すれば、両成分の振幅を足してから二乗する。

:::example[法線入射で任意の横方向を選ぶ]
$s,p$ の基底を面上の任意の直交軸に取る。

入力は屈折率比 $m=\frac{3}{2}$、入射角 $\theta_i=0\,\mathrm{rad}$、入射偏光と $s$ 軸の角 $\psi\simeq 1.0472\,\mathrm{rad}$ である。

Snell の式で屈折角を求め、外向きの枝に対応する正の余弦を使う。

$$
\begin{aligned}
\cos\theta_t&=\sqrt{1 - \frac{\sin^{2}{\left(\theta_i \right)}}{m^{2}}},\\
\cos\theta_t&= \sqrt{1 - \frac{\sin^{2}{\left(\left(0\right) \right)}}{\left(\frac{3}{2}\right)^{2}}}=1.
\end{aligned}
$$

$s$ では法線磁束を求める係数が余弦に比例する。

$$
\begin{aligned}
r_s&=\frac{- \cos\theta_t m + \cos{\left(\theta_i \right)}}{\cos\theta_t m + \cos{\left(\theta_i \right)}},\\
r_s&= \frac{- \left(1\right) \cdot \left(\frac{3}{2}\right) + \cos{\left(\left(0\right) \right)}}{\left(1\right) \cdot \left(\frac{3}{2}\right) + \cos{\left(\left(0\right) \right)}}=- \frac{1}{5}.
\end{aligned}
$$

$p$ には余弦の逆数を含むアドミタンスを使い、分母と分子から共通の余弦を払う。

$$
\begin{aligned}
r_p&=\frac{\cos\theta_t - m \cos{\left(\theta_i \right)}}{\cos\theta_t + m \cos{\left(\theta_i \right)}},\\
r_p&= \frac{\left(1\right) - \left(\frac{3}{2}\right) \cdot \cos{\left(\left(0\right) \right)}}{\left(1\right) + \left(\frac{3}{2}\right) \cdot \cos{\left(\left(0\right) \right)}}=- \frac{1}{5}.
\end{aligned}
$$

直交した成分のパワーを入射偏光の比で重み付けして足す。

$$
\begin{aligned}
R&=r_p^{2} \sin^{2}{\left(\psi \right)} + r_s^{2} \cos^{2}{\left(\psi \right)},\\
R&\simeq \left(- \frac{1}{5}\right)^{2} \cdot \sin^{2}{\left(\left(1.047197551\right) \right)} + \left(- \frac{1}{5}\right)^{2} \cdot \cos^{2}{\left(\left(1.047197551\right) \right)}=\frac{1}{25}.
\end{aligned}
$$

反射光を $s$ に透過する検光子で選んだときの割合を求める。

$$
\begin{aligned}
f_s&=\frac{r_s^{2} \cos^{2}{\left(\psi \right)}}{R},\\
f_s&\simeq \frac{\left(- \frac{1}{5}\right)^{2} \cdot \cos^{2}{\left(\left(1.047197551\right) \right)}}{\left(\frac{1}{25}\right)}=\frac{1}{4}.
\end{aligned}
$$

接線振幅の反射係数は等しいので、界面だけでは元の直線偏光は変わらない。
:::

:::example[斜めの入射を混合偏光で測る]
ほぼ等しい二つの偏光成分を持つ光を斜めに入れ、反射光を偏光子で分ける。

入力は屈折率比 $m=\frac{3}{2}$、入射角 $\theta_i\simeq 0.7854\,\mathrm{rad}$、入射偏光と $s$ 軸の角 $\psi\simeq 0.7854\,\mathrm{rad}$ である。

Snell の式で屈折角を求め、外向きの枝に対応する正の余弦を使う。

$$
\begin{aligned}
\cos\theta_t&=\sqrt{1 - \frac{\sin^{2}{\left(\theta_i \right)}}{m^{2}}},\\
\cos\theta_t&\simeq \sqrt{1 - \frac{\sin^{2}{\left(\left(0.7853981634\right) \right)}}{\left(\frac{3}{2}\right)^{2}}}\simeq 0.88192.
\end{aligned}
$$

$s$ では法線磁束を求める係数が余弦に比例する。

$$
\begin{aligned}
r_s&=\frac{- \cos\theta_t m + \cos{\left(\theta_i \right)}}{\cos\theta_t m + \cos{\left(\theta_i \right)}},\\
r_s&\simeq \frac{- \left(0.8819171037\right) \cdot \left(\frac{3}{2}\right) + \cos{\left(\left(0.7853981634\right) \right)}}{\left(0.8819171037\right) \cdot \left(\frac{3}{2}\right) + \cos{\left(\left(0.7853981634\right) \right)}}\simeq -0.30334.
\end{aligned}
$$

$p$ には余弦の逆数を含むアドミタンスを使い、分母と分子から共通の余弦を払う。

$$
\begin{aligned}
r_p&=\frac{\cos\theta_t - m \cos{\left(\theta_i \right)}}{\cos\theta_t + m \cos{\left(\theta_i \right)}},\\
r_p&\simeq \frac{\left(0.8819171037\right) - \left(\frac{3}{2}\right) \cdot \cos{\left(\left(0.7853981634\right) \right)}}{\left(0.8819171037\right) + \left(\frac{3}{2}\right) \cdot \cos{\left(\left(0.7853981634\right) \right)}}\simeq -0.092013.
\end{aligned}
$$

直交した成分のパワーを入射偏光の比で重み付けして足す。

$$
\begin{aligned}
R&=r_p^{2} \sin^{2}{\left(\psi \right)} + r_s^{2} \cos^{2}{\left(\psi \right)},\\
R&\simeq \left(-0.09201336305\right)^{2} \cdot \sin^{2}{\left(\left(0.7853981634\right) \right)} + \left(-0.3033370453\right)^{2} \cdot \cos^{2}{\left(\left(0.7853981634\right) \right)}\simeq 0.05024.
\end{aligned}
$$

反射光を $s$ に透過する検光子で選んだときの割合を求める。

$$
\begin{aligned}
f_s&=\frac{r_s^{2} \cos^{2}{\left(\psi \right)}}{R},\\
f_s&\simeq \frac{\left(-0.3033370453\right)^{2} \cdot \cos^{2}{\left(\left(0.7853981634\right) \right)}}{\left(0.05023991101\right)}\simeq 0.91574.
\end{aligned}
$$

反射光は入射光より強く $s$ に偏る。全反射パワーだけでは偏光の比を決めにくい。
:::

:::example[Brewster 角の少し外側]
ガラス側に光を入れ、$p$ を主成分にして入射表面の影響を減らす。

入力は屈折率比 $m=\frac{3}{2}$、入射角 $\theta_i\simeq 1.0472\,\mathrm{rad}$、入射偏光と $s$ 軸の角 $\psi\simeq 1.309\,\mathrm{rad}$ である。

Snell の式で屈折角を求め、外向きの枝に対応する正の余弦を使う。

$$
\begin{aligned}
\cos\theta_t&=\sqrt{1 - \frac{\sin^{2}{\left(\theta_i \right)}}{m^{2}}},\\
\cos\theta_t&\simeq \sqrt{1 - \frac{\sin^{2}{\left(\left(1.047197551\right) \right)}}{\left(\frac{3}{2}\right)^{2}}}\simeq 0.8165.
\end{aligned}
$$

$s$ では法線磁束を求める係数が余弦に比例する。

$$
\begin{aligned}
r_s&=\frac{- \cos\theta_t m + \cos{\left(\theta_i \right)}}{\cos\theta_t m + \cos{\left(\theta_i \right)}},\\
r_s&\simeq \frac{- \left(0.8164965809\right) \cdot \left(\frac{3}{2}\right) + \cos{\left(\left(1.047197551\right) \right)}}{\left(0.8164965809\right) \cdot \left(\frac{3}{2}\right) + \cos{\left(\left(1.047197551\right) \right)}}\simeq -0.4202.
\end{aligned}
$$

$p$ には余弦の逆数を含むアドミタンスを使い、分母と分子から共通の余弦を払う。

$$
\begin{aligned}
r_p&=\frac{\cos\theta_t - m \cos{\left(\theta_i \right)}}{\cos\theta_t + m \cos{\left(\theta_i \right)}},\\
r_p&\simeq \frac{\left(0.8164965809\right) - \left(\frac{3}{2}\right) \cdot \cos{\left(\left(1.047197551\right) \right)}}{\left(0.8164965809\right) + \left(\frac{3}{2}\right) \cdot \cos{\left(\left(1.047197551\right) \right)}}\simeq 0.042449.
\end{aligned}
$$

直交した成分のパワーを入射偏光の比で重み付けして足す。

$$
\begin{aligned}
R&=r_p^{2} \sin^{2}{\left(\psi \right)} + r_s^{2} \cos^{2}{\left(\psi \right)},\\
R&\simeq \left(0.04244923464\right)^{2} \cdot \sin^{2}{\left(\left(1.308996939\right) \right)} + \left(-0.4202041029\right)^{2} \cdot \cos^{2}{\left(\left(1.308996939\right) \right)}\simeq 0.013509.
\end{aligned}
$$

反射光を $s$ に透過する検光子で選んだときの割合を求める。

$$
\begin{aligned}
f_s&=\frac{r_s^{2} \cos^{2}{\left(\psi \right)}}{R},\\
f_s&\simeq \frac{\left(-0.4202041029\right)^{2} \cdot \cos^{2}{\left(\left(1.308996939\right) \right)}}{\left(0.01350927751\right)}\simeq 0.87555.
\end{aligned}
$$

入射には $s$ が少なくても、反射側では $s$ の比が大きい。反射率を抑えるには角度と偏光の両方を制御する。
:::

## Brewster 角と偏光

無損失、非磁性の異なる二媒質を考える。$p$ 成分の反射が消えるには $Y_{1p}=Y_{2p}$、すなわち $n_1\cos\theta_t=n_2\cos\theta_i$ が必要である。:eqref[eq-optics-snell] から $\sin^2\theta_t=(n_1/n_2)^2\sin^2\theta_i$ である。整合条件の二乗 $n_1^2\cos^2\theta_t=n_2^2\cos^2\theta_i$ に $n_1^2\sin^2\theta_t$ を足すと
$$
n_1^2=n_2^2\cos^2\theta_i+\frac{n_1^4}{n_2^2}\sin^2\theta_i .
$$
$n_1\ne n_2$ を用いて整理すると $\tan^2\theta_i=n_2^2/n_1^2$ である。よって
$$
\tan\theta_B=n_2/n_1,\qquad\theta_t+\theta_B=\pi/2 .
\label{eq-optics-brewster}
$$
第二式は $\sin\theta_t=(n_1/n_2)\sin\theta_B=\cos\theta_B$ による。この角を :term[Brewster 角]{en="Brewster angle" yomi="ぶりゅーすたーかく"} という。$s$ で反射をゼロにする条件は両側の $k_z$ が等しいことである。それと $k_x$ の等しさは $n_1=n_2$ を強制するので、ここでは解がない。

:::example[ガラス面での選別]
$n_2/n_1=3/2$ なら $\theta_B\simeq56.31^\circ$ である。ここで $Y_{1s}/Y_{2s}=n_1^2/n_2^2=4/9$ なので $r_s=-5/13$、$R_s=25/169$ である。等しい平均強度の $s,p$ が入れば反射率は $25/338$ となり、反射光は $s$ 偏光である。
:::

:::example[同じ角の意味]
光を逆方向から入れると、Brewster 角は $\arctan(n_1/n_2)=\pi/2-\theta_B$ になる。これは元の透過角であるため、偏光を失わずに往復する窓を設計できる。
:::

### Brewster 零点の傾き
$\theta=\theta_i$ とし、非磁性の場合の $p$ アドミタンスを $Y_1=n_1/(\eta_0\cos\theta)$、$Y_2=n_2/(\eta_0\cos\theta_t)$ とする。$Y_1=Y_2$ の点で微分すると $r_p'=(\tan\theta-\tan\theta_t\,\theta_t')/2$ である。Snell の角微分から $\theta_t'=(n_1/n_2)\cos\theta/\cos\theta_t$。Brewster 角では $m=n_2/n_1$ として $\tan\theta=m,\tan\theta_t=1/m,\theta_t'=1/m^2$。従って $r_p'=(m-m^{-3})/2$。局所の反射率は $R_p\simeq(r_p')^2(\theta-\theta_B)^2$ である。入射角に有限の広がりがあると、中心が零点でも平均反射は残る。

:::example[角度誤差]
$m=1.5$ なら $r_p'\simeq0.602\,\mathrm{rad^{-1}}$。Brewster 角から $1^\circ$ ずれると、この近似で $R_p\simeq1.10\times10^{-4}$ になる。
:::

### 磁性媒質では二つの零条件を分ける
両側で法線波数が正の実数であるとする。$s$ の無反射条件は $k_{1z}/\mu_1=k_{2z}/\mu_2$、$p$ は $k_{1z}/\eps_1=k_{2z}/\eps_2$ に等しい。そこで $m=n_2/n_1$ とし、$s$ では $b=\mu_2/\mu_1$、$p$ では $b=\eps_2/\eps_1$ と置く。両条件とも $k_{2z}=b k_{1z}$ の形であり、$s_i=\sin^2\theta_i$ とすれば $m^2-s_i=b^2(1-s_i)$ になる。従って $b\ne1$ で候補は $s_i=(b^2-m^2)/(b^2-1)$ である。これが $0\le s_i<\min(1,m^2)$ に入れば、正の法線波数の仮定により元の無反射条件も満たす。$b=1$ のときは $m=1$ ならすべての伝播角で整合し、$m\ne1$ なら整合しない。

両偏光について同じ非零斜入射で整合するなら $\mu_2/\mu_1=\eps_2/\eps_1=k_{2z}/k_{1z}=a$ である。すると $k_2^2=a^2 k_1^2$ かつ $k_{2z}^2=a^2 k_{1z}^2$。共通の $k_x$ を引き算すると $k_x^2=a^2k_x^2$、非零斜入射で正の $a$ なら $a=1$ になる。垂直入射ではこの制約がないため、インピーダンスが等しければ両偏光で無反射にできる。

:::example[磁気応答だけを変える]
$\eps_1=\eps_2$、$\mu_2/\mu_1=2$ なら $n_2/n_1=\sqrt2$。$s$ で $k_{2z}=2k_{1z}$ とすると $2-\sin^2\theta_i=4(1-\sin^2\theta_i)$、従って $\sin^2\theta_i=2/3$ で無反射となる。この角では $p$ のアドミタンスは一致しない。
:::

### 角度分布と偏光の漏れで消え残る反射

Brewster 角近傍で $p$ の反射係数が $a\delta+O(\delta^2)$、角度のずれが $\delta=\theta_i-\theta_B$ であるとする。角度分布の平均が $\delta_0$、分散が $\sigma_\theta^2$ なら、$\langle\delta^2\rangle=\delta_0^2+\sigma_\theta^2$。従ってビーム全体の反射率を角度成分のパワーで平均すると、先頭項は $a^2(\delta_0^2+\sigma_\theta^2)$ である。平均角だけを調整しても、分散の項は残る。Fourier 面で角度ごとのパワーを足す場合、この平均は Parseval の等式から反射ビーム全体のパワーと等価である。

入射偏光に $s$ のパワー比 $f$ が混じるなら、$s$ 側では零点でも $R_{sB}=[(1-m^2)/(1+m^2)]^2$ が残る。狭いビームでは $R_s$ を中心で一定とした一次の見積りにより
$$
\langle R\rangle \simeq (1-f)a^2(\delta_0^2+\sigma_\theta^2)+f R_{sB}.
$$
この式は角度偏差が小さく、偏光の混合が角度によらず、二偏光を区別せずパワーを測る場合に使う。$f$ を下げる対策とビームを広くして角度幅を抑える対策は、右辺のどちらが大きいかで選べる。$p$ の零点のまわりで $R_s$ の変化も必要なら、同じ角分布で Fresnel 係数を直接平均する。

:::example[ビームの発散で決まる残留反射]
偏光の純度を上げて平均角を正しく合わせたとき、角度幅が決める反射を見積もる。

入力は屈折率比 $m=\frac{3}{2}$、平均角のずれ $\delta_0=0\,\mathrm{rad}$、角度の標準偏差 $\sigma_\theta=\frac{3}{1000}\,\mathrm{rad}$、入射での $s$ の割合 $f\simeq 1\times 10^{-7}$ である。

前に求めた $p$ の接線振幅の角微分を入れる。

$$
\begin{aligned}
a&=\frac{m}{2} - \frac{1}{2 m^{3}},\\
a/(\mathrm{rad^{-1}})&= \frac{\left(\frac{3}{2}\right)}{2} - \frac{1}{2 \cdot \left(\frac{3}{2}\right)^{3}}=\frac{65}{108}.
\end{aligned}
$$

Brewster 条件を $s$ の式に代入し、振幅を二乗する。

$$
\begin{aligned}
R_{sB}&=\frac{\left(1 - m^{2}\right)^{2}}{\left(m^{2} + 1\right)^{2}},\\
R_{sB}&= \frac{\left(1 - \left(\frac{3}{2}\right)^{2}\right)^{2}}{\left(\left(\frac{3}{2}\right)^{2} + 1\right)^{2}}=\frac{25}{169}.
\end{aligned}
$$

ビームの平均角のずれと角の分散を足し、入射の $p$ 分に掛ける。

$$
\begin{aligned}
R_{\mathrm{angle}}&=a^{2} \left(1 - f\right) \left(\delta_0^{2} + \sigma_\theta^{2}\right),\\
R_{\mathrm{angle}}&\simeq \left(\frac{65}{108}\right)^{2} \cdot \left(1 - \left(1\times 10^{-7}\right)\right) \cdot \left(\left(0\right)^{2} + \left(\frac{3}{1000}\right)^{2}\right)\simeq 3.26\times 10^{-6}.
\end{aligned}
$$

入射時に含まれる $s$ の分を反射率で重み付けする。

$$
\begin{aligned}
R_{\mathrm{leak}}&=R_{sB} f,\\
R_{\mathrm{leak}}&\simeq \left(\frac{25}{169}\right) \cdot \left(1\times 10^{-7}\right)\simeq 1.4793\times 10^{-8}.
\end{aligned}
$$

反射の総パワーでは直交成分の寄与を足す。

$$
\begin{aligned}
\langle R\rangle&=R_{\mathrm{leak}} + R_{\mathrm{angle}},\\
\langle R\rangle&\simeq \left(1.479289941\times 10^{-8}\right) + \left(3.260030538\times 10^{-6}\right)\simeq 3.2748\times 10^{-6}.
\end{aligned}
$$

角分布の寄与が大きい。偏光をさらに整えても、この反射率はほとんど変わらない。
:::

:::example[偏光の漏れが残る配置]
細い角度分布を用いるが、$s$ 成分の除去に限界がある。

入力は屈折率比 $m=\frac{3}{2}$、平均角のずれ $\delta_0=0\,\mathrm{rad}$、角度の標準偏差 $\sigma_\theta=\frac{1}{1000}\,\mathrm{rad}$、入射での $s$ の割合 $f=\frac{1}{1000}$ である。

前に求めた $p$ の接線振幅の角微分を入れる。

$$
\begin{aligned}
a&=\frac{m}{2} - \frac{1}{2 m^{3}},\\
a/(\mathrm{rad^{-1}})&= \frac{\left(\frac{3}{2}\right)}{2} - \frac{1}{2 \cdot \left(\frac{3}{2}\right)^{3}}=\frac{65}{108}.
\end{aligned}
$$

Brewster 条件を $s$ の式に代入し、振幅を二乗する。

$$
\begin{aligned}
R_{sB}&=\frac{\left(1 - m^{2}\right)^{2}}{\left(m^{2} + 1\right)^{2}},\\
R_{sB}&= \frac{\left(1 - \left(\frac{3}{2}\right)^{2}\right)^{2}}{\left(\left(\frac{3}{2}\right)^{2} + 1\right)^{2}}=\frac{25}{169}.
\end{aligned}
$$

ビームの平均角のずれと角の分散を足し、入射の $p$ 分に掛ける。

$$
\begin{aligned}
R_{\mathrm{angle}}&=a^{2} \left(1 - f\right) \left(\delta_0^{2} + \sigma_\theta^{2}\right),\\
R_{\mathrm{angle}}&= \left(\frac{65}{108}\right)^{2} \cdot \left(1 - \left(\frac{1}{1000}\right)\right) \cdot \left(\left(0\right)^{2} + \left(\frac{1}{1000}\right)^{2}\right)\simeq 3.6186\times 10^{-7}.
\end{aligned}
$$

入射時に含まれる $s$ の分を反射率で重み付けする。

$$
\begin{aligned}
R_{\mathrm{leak}}&=R_{sB} f,\\
R_{\mathrm{leak}}&= \left(\frac{25}{169}\right) \cdot \left(\frac{1}{1000}\right)=\frac{1}{6760}.
\end{aligned}
$$

反射の総パワーでは直交成分の寄与を足す。

$$
\begin{aligned}
\langle R\rangle&=R_{\mathrm{leak}} + R_{\mathrm{angle}},\\
\langle R\rangle&\simeq \left(\frac{1}{6760}\right) + \left(3.618634259\times 10^{-7}\right)\simeq 1.4829\times 10^{-4}.
\end{aligned}
$$

反射パワーは検光子の漏れによって決まる。平均角を掃引しても最小値には床が残る。
:::

:::example[窓の傾きを調整する]
屈折率比が大きい系で角度設定の誤差がある。

入力は屈折率比 $m=2$、平均角のずれ $\delta_0=\frac{1}{100}\,\mathrm{rad}$、角度の標準偏差 $\sigma_\theta=\frac{1}{1000}\,\mathrm{rad}$、入射での $s$ の割合 $f\simeq 1\times 10^{-4}$ である。

前に求めた $p$ の接線振幅の角微分を入れる。

$$
\begin{aligned}
a&=\frac{m}{2} - \frac{1}{2 m^{3}},\\
a/(\mathrm{rad^{-1}})&= \frac{\left(2\right)}{2} - \frac{1}{2 \cdot \left(2\right)^{3}}=\frac{15}{16}.
\end{aligned}
$$

Brewster 条件を $s$ の式に代入し、振幅を二乗する。

$$
\begin{aligned}
R_{sB}&=\frac{\left(1 - m^{2}\right)^{2}}{\left(m^{2} + 1\right)^{2}},\\
R_{sB}&= \frac{\left(1 - \left(2\right)^{2}\right)^{2}}{\left(\left(2\right)^{2} + 1\right)^{2}}=\frac{9}{25}.
\end{aligned}
$$

ビームの平均角のずれと角の分散を足し、入射の $p$ 分に掛ける。

$$
\begin{aligned}
R_{\mathrm{angle}}&=a^{2} \left(1 - f\right) \left(\delta_0^{2} + \sigma_\theta^{2}\right),\\
R_{\mathrm{angle}}&\simeq \left(\frac{15}{16}\right)^{2} \cdot \left(1 - \left(1\times 10^{-4}\right)\right) \cdot \left(\left(\frac{1}{100}\right)^{2} + \left(\frac{1}{1000}\right)^{2}\right)\simeq 8.8761\times 10^{-5}.
\end{aligned}
$$

入射時に含まれる $s$ の分を反射率で重み付けする。

$$
\begin{aligned}
R_{\mathrm{leak}}&=R_{sB} f,\\
R_{\mathrm{leak}}&\simeq \left(\frac{9}{25}\right) \cdot \left(1\times 10^{-4}\right)\simeq 3.6\times 10^{-5}.
\end{aligned}
$$

反射の総パワーでは直交成分の寄与を足す。

$$
\begin{aligned}
\langle R\rangle&=R_{\mathrm{leak}} + R_{\mathrm{angle}},\\
\langle R\rangle&\simeq \left(3.6\times 10^{-5}\right) + \left(8.87606543\times 10^{-5}\right)\simeq 1.2476\times 10^{-4}.
\end{aligned}
$$

零点の傾きが大きく、角設定の改善で最も大きな低減が得られる範囲である。
:::

## 全反射とエバネッセント波

$n_1>n_2$ のとき、$\sin\theta_i>n_2/n_1$ では媒質 2 の法線波数が虚数になる。臨界角は $\theta_c=\arcsin(n_2/n_1)$ である。$z>0$ で増大する解を除き、$k_{2z}=i\kappa$、
$$
\kappa=\sqrt{k_x^2-k_2^2},\qquad e^{ik_x x+ik_{2z}z}=e^{ik_x x-\kappa z}
\label{eq-optics-evanescent}
$$
を採る。これを :term[エバネッセント波]{en="evanescent wave" yomi="えばねっせんとは"} と呼ぶ。$Y_{2s},Y_{2p}$ は純虚数なので、:eqref[eq-optics-fresnel] で $|r_s|=|r_p|=1$ となる。これは :term[全反射]{en="total internal reflection" yomi="ぜんはんしゃ"} である。例えば $s$ 波では $H_z=k_x E_y/(\omega\mu_2)$ なので、境界に平行な平均流 $\langle S_x\rangle=k_x|E_y|^2/(2\omega\mu_2)$ が残る。

非磁性の $s$ 波の反射位相は $\phi_s=-2\arctan(\kappa/k_{1z})$ である。$p$ の接線基底では $\phi_p=2\arctan(\eps_2 k_{1z}/(\eps_1\kappa))$ である。臨界角を越えて二つの位相が別々に動くので、振幅比を変えずに偏光を楕円化できる。

:::example[しみ出しの距離]
$n_1=1.5,n_2=1,\lambda_0=600\,\mathrm{nm},\theta_i=60^\circ$ なら $\kappa=k_0\sqrt{11}/4$ である。振幅の $1/e$ 深さは約 $115\,\mathrm{nm}$ であり、強度の $1/e$ 深さはその半分である。
:::

:::example[第 3 の媒質]
低屈折率層が厚さ $d$ で終わり、その先に伝播できる媒質を置くと、低屈折率層の解は $e^{-\kappa z}$ と $e^{+\kappa z}$ の両方を含む。単一境界の外向き条件はもはや使えない。$\kappa d\gg1$ なら透過パワーには往路の振幅因子の二乗 $e^{-2\kappa d}$ が現れる。
:::

### 低屈折率の薄層を通り抜ける場合
両外側が同じ媒質で $s$ 波を考え、外側の法線波数を $q>0$、間の厚さ $d$ の低屈折率層を $i\kappa$ とする。非磁性なら法線微分と場が連続である。層の右端場を $t$、右端勾配を $iqt$ とし、層の方程式 $E''=\kappa^2 E$ で左端へ戻すと
$$
\binom{E(0)}{E'(0)}
=\begin{pmatrix}\cosh\kappa d&-\sinh(\kappa d)/\kappa\\-\kappa\sinh(\kappa d)&\cosh\kappa d\end{pmatrix}\binom{t}{iqt}.
$$
左端で $E(0)=1+r$、$E'(0)=iq(1-r)$ なので、$E(0)+E'(0)/(iq)=2$。よって $t=1/[\cosh\kappa d+(i/2)(\kappa/q-q/\kappa)\sinh\kappa d]$。パワー透過率は
$$
T=\frac1{1+[(q^2+\kappa^2)/(2q\kappa)]^2\sinh^2\kappa d}
$$
である。半無限にした極限で初めて $T=0$ になる。この形なら厚さゼロで $T=1$ となることも確認できる。

:::example[同じ尺度の波数]
$q=\kappa,\ \kappa d=\ln2$ のとき $T=\cosh^{-2}(\ln2)=16/25$。層内だけを見ればエバネッセント域でも、有限の厚さを介した伝播パワーがある。
:::

### 位相差を遅相器として使う
$k_{1z}=q$、$\eps_1/\eps_2=a$ とすると、全反射の二つの位相は $\phi_s=-2\arctan(\kappa/q)$、$\phi_p=\pi-2\arctan(a\kappa/q)$ である。反射面に固定した接線基底での $\pi$ は、反射側で $p$ の基底の向きを反転すれば取り除ける。後者の約束では実質的な遅相の大きさは $2[\arctan(a\kappa/q)-\arctan(\kappa/q)]$ である。この量は臨界角とすれすれ入射（$\theta_i\to\pi/2$）の極限でゼロになり、中間で最大になる。$x=\kappa/q$ で微分すると零点は $x=1/\sqrt a$。同じ入射面で全反射を二回起こす :term[Fresnel 菱面体]{en="Fresnel rhomb" yomi="ふれねるりょうめんたい"} では、この相対位相が二倍になる。二反射で四分の一周期を得るには一反射の位相差の大きさを $\pi/4$ にする。加法定理から $\tan(\pi/8)=(a-1)x/(1+a x^2)$ であり、$x$ を求めれば $\kappa/q=x$ から入射角が決まる。

:::example[二回の全反射で四分の一周期の位相差を作る]
$n_1=1.5$、$n_2=1$ なら $a=9/4$。$\tau=\tan(\pi/8)=\sqrt2-1$ と置くと $a\tau x^2-(a-1)x+\tau=0$ である。二つの正の根を $\sin^2\theta_i=(x^2+1/a)/(1+x^2)$ に代入すると、入射角は約 $50.23^\circ$ または $53.26^\circ$ となる。等振幅の $s,p$ を入れると、二反射後の位相差は四分の一周期であり、[偏光と結晶光学](/physics/optics/polarization/) で扱う円偏光の条件を満たす。
:::

### 全反射面に沿うパワーとその深さ

$s$ について透過側を $\vb{\mathcal E}=\hat{\vb y}E_b e^{iKx-\kappa z}$ とすると $\vb{\mathcal H}=(-i\kappa\,\hat{\vb x}+K\,\hat{\vb z})E_b e^{iKx-\kappa z}/(\omega\mu_2)$ である。複素 Poynting 流の $x$ 成分は $K|E_b|^2 e^{-2\kappa z}/(2\omega\mu_2)$、$z$ 成分はその式の $K$ を $-i\kappa$ に置き換えたものになる。従って法線の流れは実部が消え、接線の実部は深さとともに減る。場を平面波でモデル化する限りこの平均流は $x$ によらず、法線にも平均流がないので、無損失媒質で流れの発散はゼロである。接線流を深さ方向に積分すると
$$
\int_0^\infty \langle S_x\rangle\dd z=\frac{K|E_b|^2}{4\omega\mu_2\kappa}
$$
となる。単位は $\mathrm{W/m}$ であり、ビームの幅にわたる全パワーを出すには残る座標でも積分する。

厚さ $a$ までに含まれる割合は $1-e^{-2\kappa a}$。厚さ $a$ での値そのものは境界での値の $e^{-2\kappa a}$ であり、積分量と局所量は別の測定に対応する。薄い試料を界面に接近させるセンサーでは、試料が場をほとんど乱さなければその体積で場の二乗を積分して感度を見積もる。吸収のある試料を厚く載せる場合にはこの場自体が変わるため、元の全反射の式を固定せず新しい境界条件を解く。

:::example[臨界角のそばで厚さを測る]
ガラスから空気へ入り、外に置いた薄い試料が場を乱さないとする。

入力は入射側屈折率 $n_{1}=\frac{3}{2}$、外側屈折率 $n_{2}=1$、真空波長 $\lambda_0\simeq 6\times 10^{-7}\,\mathrm m$、入射角 $\theta_i\simeq 0.7854\,\mathrm{rad}$、試料の厚さ $a\simeq 5\times 10^{-8}\,\mathrm m$ である。

界面の接線成分を入射波から取る。

$$
\begin{aligned}
K&=\frac{2 \pi n_{1} \sin{\left(\theta_i \right)}}{\lambda_0},\\
K/(\mathrm{m^{-1}})&\simeq \frac{2 \cdot \pi \cdot \left(\frac{3}{2}\right) \cdot \sin{\left(\left(0.7853981634\right) \right)}}{\left(6\times 10^{-7}\right)}\simeq 1.1107\times 10^{7}.
\end{aligned}
$$

外側で伝播できないため、法線波数の虚部の正値を求める。

$$
\begin{aligned}
\kappa&=\sqrt{K^{2} - \frac{4 \pi^{2} n_{2}^{2}}{\lambda_0^{2}}},\\
\kappa/(\mathrm{m^{-1}})&\simeq \sqrt{\left(\left(1.110720735\times 10^{7}\right)\right)^{2} - \frac{4 \cdot \pi^{2} \cdot \left(1\right)^{2}}{\left(\left(6\times 10^{-7}\right)\right)^{2}}}\simeq 3.7024\times 10^{6}.
\end{aligned}
$$

強度は振幅の二乗なので、指数率は振幅の二倍になる。

$$
\begin{aligned}
\delta_I&=\frac{1}{2 \kappa},\\
\delta_I/(\mathrm m)&\simeq \frac{1}{2 \cdot \left(3.702402448\times 10^{6}\right)}\simeq 1.3505\times 10^{-7}.
\end{aligned}
$$

境界から厚さ $a$ の位置に進んだ強度を評価する。

$$
\begin{aligned}
I(a)/I(0)&=e^{- \frac{a}{\delta_I}},\\
I(a)/I(0)&\simeq e^{- \frac{\left(5\times 10^{-8}\right)}{\left(1.350474474\times 10^{-7}\right)}}\simeq 0.69057.
\end{aligned}
$$

指数減衰を積分した総流から、試料より奥の流れを引く。

$$
\begin{aligned}
f_a&=1 - I(a)/I(0),\\
f_a&\simeq 1 - \left(0.6905684052\right)\simeq 0.30943.
\end{aligned}
$$

試料は侵入長より薄く、接線流のかなりの部分が試料の奥へ分布する。
:::

:::example[角度を増して表面を選ぶ]
同じ厚さの試料を使い、入射角を増やす。

入力は入射側屈折率 $n_{1}=\frac{3}{2}$、外側屈折率 $n_{2}=1$、真空波長 $\lambda_0\simeq 6\times 10^{-7}\,\mathrm m$、入射角 $\theta_i\simeq 1.0472\,\mathrm{rad}$、試料の厚さ $a\simeq 5\times 10^{-8}\,\mathrm m$ である。

界面の接線成分を入射波から取る。

$$
\begin{aligned}
K&=\frac{2 \pi n_{1} \sin{\left(\theta_i \right)}}{\lambda_0},\\
K/(\mathrm{m^{-1}})&\simeq \frac{2 \cdot \pi \cdot \left(\frac{3}{2}\right) \cdot \sin{\left(\left(1.047197551\right) \right)}}{\left(6\times 10^{-7}\right)}\simeq 1.3603\times 10^{7}.
\end{aligned}
$$

外側で伝播できないため、法線波数の虚部の正値を求める。

$$
\begin{aligned}
\kappa&=\sqrt{K^{2} - \frac{4 \pi^{2} n_{2}^{2}}{\lambda_0^{2}}},\\
\kappa/(\mathrm{m^{-1}})&\simeq \sqrt{\left(\left(1.360349523\times 10^{7}\right)\right)^{2} - \frac{4 \cdot \pi^{2} \cdot \left(1\right)^{2}}{\left(\left(6\times 10^{-7}\right)\right)^{2}}}\simeq 8.6829\times 10^{6}.
\end{aligned}
$$

強度は振幅の二乗なので、指数率は振幅の二倍になる。

$$
\begin{aligned}
\delta_I&=\frac{1}{2 \kappa},\\
\delta_I/(\mathrm m)&\simeq \frac{1}{2 \cdot \left(8.682903397\times 10^{6}\right)}\simeq 5.7584\times 10^{-8}.
\end{aligned}
$$

境界から厚さ $a$ の位置に進んだ強度を評価する。

$$
\begin{aligned}
I(a)/I(0)&=e^{- \frac{a}{\delta_I}},\\
I(a)/I(0)&\simeq e^{- \frac{\left(5\times 10^{-8}\right)}{\left(5.758442507\times 10^{-8}\right)}}\simeq 0.41967.
\end{aligned}
$$

指数減衰を積分した総流から、試料より奥の流れを引く。

$$
\begin{aligned}
f_a&=1 - I(a)/I(0),\\
f_a&\simeq 1 - \left(0.4196684267\right)\simeq 0.58033.
\end{aligned}
$$

角度を増すと侵入長が減り、試料の占める割合が増える。境界に生じる場の大きさも比較するには $t_s$ の変化を含める。
:::

:::example[長い波長で同じ角度を使う]
材料の分散を無視し、波長だけを二倍にする。

入力は入射側屈折率 $n_{1}=\frac{3}{2}$、外側屈折率 $n_{2}=1$、真空波長 $\lambda_0\simeq 1.2\times 10^{-6}\,\mathrm m$、入射角 $\theta_i\simeq 1.0472\,\mathrm{rad}$、試料の厚さ $a\simeq 5\times 10^{-8}\,\mathrm m$ である。

界面の接線成分を入射波から取る。

$$
\begin{aligned}
K&=\frac{2 \pi n_{1} \sin{\left(\theta_i \right)}}{\lambda_0},\\
K/(\mathrm{m^{-1}})&\simeq \frac{2 \cdot \pi \cdot \left(\frac{3}{2}\right) \cdot \sin{\left(\left(1.047197551\right) \right)}}{\left(1.2\times 10^{-6}\right)}\simeq 6.8017\times 10^{6}.
\end{aligned}
$$

外側で伝播できないため、法線波数の虚部の正値を求める。

$$
\begin{aligned}
\kappa&=\sqrt{K^{2} - \frac{4 \pi^{2} n_{2}^{2}}{\lambda_0^{2}}},\\
\kappa/(\mathrm{m^{-1}})&\simeq \sqrt{\left(\left(6.801747616\times 10^{6}\right)\right)^{2} - \frac{4 \cdot \pi^{2} \cdot \left(1\right)^{2}}{\left(\left(1.2\times 10^{-6}\right)\right)^{2}}}\simeq 4.3415\times 10^{6}.
\end{aligned}
$$

強度は振幅の二乗なので、指数率は振幅の二倍になる。

$$
\begin{aligned}
\delta_I&=\frac{1}{2 \kappa},\\
\delta_I/(\mathrm m)&\simeq \frac{1}{2 \cdot \left(4.341451698\times 10^{6}\right)}\simeq 1.1517\times 10^{-7}.
\end{aligned}
$$

境界から厚さ $a$ の位置に進んだ強度を評価する。

$$
\begin{aligned}
I(a)/I(0)&=e^{- \frac{a}{\delta_I}},\\
I(a)/I(0)&\simeq e^{- \frac{\left(5\times 10^{-8}\right)}{\left(1.151688501\times 10^{-7}\right)}}\simeq 0.64782.
\end{aligned}
$$

指数減衰を積分した総流から、試料より奥の流れを引く。

$$
\begin{aligned}
f_a&=1 - I(a)/I(0),\\
f_a&\simeq 1 - \left(0.647818205\right)\simeq 0.35218.
\end{aligned}
$$

侵入長は二倍になる。表面に近い部分だけを見たい場合は、波長選択も効く。
:::

## Goos–Hänchen シフト

幅の広いビームは中心接線波数 $K$ の近くに分布 $A(q)$ を持つ平面波の和で表せる。境界上で反射場は $\int A(q)r(K+q)e^{i(K+q)x}\dd q$ に比例する。全反射域内で $\phi(K+q)=\phi(K)+q\phi'(K)+O(q^2)$ と展開すると、一次の項は入射包絡線の引数を $x+\phi'(K)$ へ置き換える。従って :term[Goos–Hänchen シフト]{en="Goos–Hänchen shift" yomi="ぐーすへんへんしふと"} の境界上の変位は
$$
D_j=-\dv{\phi_j}{k_x}\bigg|_{K}.
\label{eq-optics-gh}
$$
微分は $\omega$ 一定で行う。$\kappa'=k_x/\kappa$、$k_{1z}'=-k_x/k_{1z}$ だから、
$$
\dv{}{k_x}\arctan\frac{\kappa}{k_{1z}}
=\frac{(k_x/\kappa)k_{1z}+\kappa k_x/k_{1z}}{k_{1z}^2+\kappa^2}
=\frac{k_x}{\kappa k_{1z}}.
$$
ゆえに $D_s=2K/(\kappa k_{1z})$ である。同様に $a=\eps_1/\eps_2$ とおき $\phi_p=\pi-2\arctan(a\kappa/k_{1z})$ を微分すると、$D_p/D_s=a(k_{1z}^2+\kappa^2)/(k_{1z}^2+a^2\kappa^2)$ となる。

:::example[変位とビーム幅の比較]
前節と同じ角で $K/k_{1z}=\sqrt3$ なので $D_s=2\sqrt3/\kappa\simeq399\,\mathrm{nm}$ である。幅が数十 $\mathrm{\mu m}$ のビームなら変位は幅より十分小さい。
:::

:::example[臨界角近傍の限界]
中心と臨界波数の差 $K-k_2$ がスペクトル幅 $\Delta K$ と同程度なら、一部の成分に $|r|<1$ が現れる。この場合は $r(K+q)$ を積分に残す必要があり、$\kappa\to0$ で上式の変位が発散するという外挿は使えない。
:::

:::example[透過側の法線減衰との比較]
全反射の $s$ 波で進行方向を固定し、$\kappa$ のみが小さいとき、変位は $2\tan\theta_i/\kappa$ である。例えば $\tan\theta_i=1$ なら反射点の変位は振幅の侵入長の二倍になる。一方、反射面に垂直な断面での強度の侵入長は $1/(2\kappa)$ なので、その深さの四倍である。
:::

### ビーム重心の積分式と角度偏向

境界上の反射場を $F(x)=\int B(k)e^{ikx}\dd k$ と書く。$B$ が滑らかに減衰すると、部分積分で $xF(x)=i\int B'(k)e^{ikx}\dd k$。Parseval の等式を用いれば $|F|^2$ で重み付けした位置の平均は $\Re[i\int B^*B'\dd k]/\int|B|^2\dd k$ となる。$B=A r$ とし、$A$ を実数の入射スペクトル、$r=|r|e^{i\phi}$ と置く。$B^*B'$ の虚部は $|B|^2\phi'$ だから、
$$
\bar x=-\frac{\int |A|^2|r|^2\phi'\dd k}{\int |A|^2|r|^2\dd k}
$$
である。反射係数に零点があるときは位相を無理に連続へ延長せず、元の $B^*B'$ の形を使う。全反射域で狭いスペクトルなら平均を中心値で置き換えられ、:eqref[eq-optics-gh] が戻る。ここで定義した重心は境界上の横電場の二乗による。法線流で重み付けする場合は波数ごとの波動アドミタンスも効く。

一方、$|r(k)|$ の傾きは平均の波数を変える。$|A|^2$ が中心 $K$ に関して対称で分散 $\sigma_k^2$ を持つとする。$R(k)=|r|^2$ を一次展開して平均 $k-K$ を求めると、分子の定数項は奇関数で消え、一次項が $R'(K)\sigma_k^2$ を与える。分母の先頭は $R(K)$ なので $\Delta k\simeq\sigma_k^2 R'(K)/R(K)$。全反射域の平面波だけなら $R'=0$ であり、この機構による角度偏向は出ない。$R(K)=0$ 付近では一次式を割って用いず、分子と分母を積分する。

:::example[二偏光を別の検出器に入れる]
二媒質は非磁性で、誘電率の比を $a=9/4$ とする。法線波数 $q=7.5\times10^6\,\mathrm{m^{-1}}$、減衰定数 $\kappa=8.3\times10^6\,\mathrm{m^{-1}}$ を選ぶ。分散式 $K^2+q^2=k_1^2=a k_2^2=a(K^2-\kappa^2)$ から
$$
K=\sqrt{\frac{q^2+a\kappa^2}{a-1}}\simeq 1.30001\times10^7\,\mathrm{m^{-1}}.
$$
ここで $\kappa/q=83/75$ を位相の微分で得た式に入れると
$$
D_s=\frac{2K}{\kappa q}\simeq 417.673\,\mathrm{nm},\qquad
\frac{D_p}{D_s}\simeq 0.695213,\qquad
D_p\simeq 290.372\,\mathrm{nm} .
$$
境界上のビーム幅を $w=100\,\mathrm{\mu m}$ とすると $D_s/w\simeq 0.00417673$ である。反射率は両偏光で一だが、$p$ の重心の変位は $s$ より小さい。
:::

:::example[臨界角に近い反射]
二媒質は非磁性で、誘電率の比を $a=9/4$ とする。法線波数 $q=11\times10^6\,\mathrm{m^{-1}}$、減衰定数 $\kappa=2.5\times10^6\,\mathrm{m^{-1}}$ を選ぶ。分散式 $K^2+q^2=k_1^2=a k_2^2=a(K^2-\kappa^2)$ から
$$
K=\sqrt{\frac{q^2+a\kappa^2}{a-1}}\simeq 1.03947\times10^7\,\mathrm{m^{-1}}.
$$
ここで $\kappa/q=5/22$ を位相の微分で得た式に入れると
$$
D_s=\frac{2K}{\kappa q}\simeq 755.979\,\mathrm{nm},\qquad
\frac{D_p}{D_s}\simeq 1.87573,\qquad
D_p\simeq 1418.01\,\mathrm{nm} .
$$
境界上のビーム幅を $w=200\,\mathrm{\mu m}$ とすると $D_s/w\simeq 0.00377989$ である。$p$ の変位は $s$ より大きい。臨界波数との差 $K-k_2$ は約 $3.05\times10^5\,\mathrm{m^{-1}}$ であり、ビームのスペクトル幅の尺度 $1/w=5.0\times10^3\,\mathrm{m^{-1}}$ より大きい。
:::

:::example[変位比が一になる位置]
二媒質は非磁性で、誘電率の比を $a=9/4$ とする。法線波数 $q=9\times10^6\,\mathrm{m^{-1}}$、減衰定数 $\kappa=6\times10^6\,\mathrm{m^{-1}}$ を選ぶ。分散式 $K^2+q^2=k_1^2=a k_2^2=a(K^2-\kappa^2)$ から
$$
K=\sqrt{\frac{q^2+a\kappa^2}{a-1}}\simeq 1.13842\times10^7\,\mathrm{m^{-1}}.
$$
ここで $\kappa/q=2/3$ を位相の微分で得た式に入れると
$$
D_s=\frac{2K}{\kappa q}\simeq 421.637\,\mathrm{nm},\qquad
\frac{D_p}{D_s}=1,\qquad
D_p\simeq 421.637\,\mathrm{nm} .
$$
境界上のビーム幅を $w=100\,\mathrm{\mu m}$ とすると $D_s/w\simeq 0.00421637$ である。$a(\kappa/q)^2=1$ なので、変位比は厳密に 1 である。
:::

## 金属面での反射

非磁性で吸収を持つ半無限媒質を $\tilde n=n+i b$、$b>0$ と表す。垂直入射では $e^{ik_0\tilde n z}$ が減衰し、斜入射では $k_{2z}^2=k_0^2\tilde n^2-k_x^2$ で $\Im k_{2z}>0$ の枝を選ぶ。垂直入射、真空側 $Y_1=1/\eta_0$ なら
$$
R=\left|\frac{1-\tilde n}{1+\tilde n}\right|^2
=\frac{(1-n)^2+b^2}{(1+n)^2+b^2}.
\label{eq-optics-metal}
$$
$1-R=4n/[(1+n)^2+b^2]$ は媒質内部へ入り熱になる割合である。$b$ が大きいと場の侵入が短くなり、反射係数は $-1$ に近づく。

### 導電率で表す極限
実導電率 $\sigma$ に対し $\rot\vb H=(\sigma-i\omega\eps)\vb E$ なので、実誘電率を使っていた式に $\eps+i\sigma/\omega$ を代入できる。$\sigma\gg\omega\eps$ なら $k_2^2\simeq i\mu\sigma\omega$ であり、
$$
k_2\simeq(1+i)/\delta,\qquad \delta=\sqrt{2/(\mu\sigma\omega)} .
$$
$\delta$ を :term[表皮深さ]{en="skin depth" yomi="ひょうひふかさ"} と呼ぶ。光周波数で一定の直流導電率を仮定できるとは限らず、周波数応答を用いる。

:::example[複素屈折率の計算]
$\tilde n=1+3i$ なら $R=9/13\simeq0.692$、$1-R=4/13$ である。$b$ の値だけでは吸収率は決まらず、界面へどれだけ電力が入るかを求める必要がある。
:::

:::example[周波数の尺度]
導電率が一定の領域内で周波数を 4 倍にすると、$\delta$ は半分になる。反射係数を比較する場合は、さらに $Y_2=k_2/(\omega\mu)$ の周波数依存性を用いる。
:::

### Poynting ベクトルと熱の整合
導電性媒質の垂直入射で、場を $E(z)=E_b e^{ik_z z}$ とする。伝播中の平均 Joule 発熱は $\sigma|E(z)|^2/2$。$\Im k_z=b>0$ なら半無限体の総発熱は単位面積当たり $\sigma|E_b|^2/(4b)$ である。一方、$k_z=a+ib$ を $k_z^2=\omega^2\mu\eps+i\omega\mu\sigma$ へ代入すると $2ab=\omega\mu\sigma$。従って $Y=a/(\omega\mu)+ib/(\omega\mu)$ の実部に $|E_b|^2/2$ を掛けた値と発熱の積分が一致する。表皮厚に達する前の場を無視してはいけない理由が、積分から分かる。

:::example[導電率を大きくする]
良導体近似では $a=b\propto\sqrt{\sigma}$ であり $Y_2\propto\sqrt{\sigma}$。真空側の $Y_1$ が固定なら $t=2Y_1/(Y_1+Y_2)$ の絶対値は $\sigma^{-1/2}$ に比例して減る。そのため総発熱は $\sigma |t|^2/b\propto\sigma^{-1/2}$ となり、完全導体の極限で吸収はゼロになる。
:::

:::example[吸収のある媒質に Brewster 零点は残るか]
真空を媒質 1、非磁性金属の複素比誘電率を $\eps_r$ とし、$p$ のアドミタンスを等しいと仮定する。$\sin^2\theta=s$ とおくと、等式を二乗して $\eps_r-s=\eps_r^2(1-s)$、すなわち $\eps_r\ne1$ で $s=\eps_r/(1+\eps_r)$ となる。$\Im\eps_r>0$ ならこの右辺の虚部は $\Im\eps_r/|1+\eps_r|^2>0$。左辺は実で矛盾する。従って反射率に極小はあり得ても、ここで扱う等方吸収モデルでは $p$ 反射の厳密な零点はない。
:::

### 複素反射の振幅と位相から光学定数を戻す

垂直入射の半無限媒質を考える。真空側で接線成分を基準にすれば $r=(1-\tilde n)/(1+\tilde n)$ なので、$r\ne-1$ で
$$
\tilde n=\frac{1-r}{1+r}
=\frac{1-|r|^2-2i\Im r}{|1+r|^2}.
$$
第二式は分母の共役を掛けて得る。従って吸収で $\Im\tilde n>0$ なら、時間規約 $e^{-i\omega t}$ のもとで $\Im r<0$ である。位相も求めるには基準面が必要であり、界面から距離 $d$ だけ手前へ基準を移すと $r$ に $e^{2ik_0d}$ が掛かる。基準面がずれれば $|r|$ は同じでも逆算した $n,b$ は一般に変わる。

強度だけを測った場合に残る曖昧さは、:eqref[eq-optics-metal] を整理すると見える。$0<R<1$ の等反射率曲線は
$$
\left(n-\frac{1+R}{1-R}\right)^2+b^2
=\frac{4R}{(1-R)^2}
$$
という円である。受動性が選ぶ上半円にも複数の候補があり、厚さや位相の測定が必要になる。複素振幅からの復元でも、$\dd\tilde n/\dd r=-2/(1+r)^2$ なので完全反射の極限に近いほど微小な位相誤差が増幅される。この式は金属膜の薄さによる背面反射が無視できる場合に限って用いる。

:::example[小さな位相遅れから吸収を求める]
基準面が校正できた複素反射測定を使う。

入力は反射係数の実部 $\Re r=- \frac{2}{5}$、反射係数の虚部 $\Im r=- \frac{1}{10}$、真空波長 $\lambda_0\simeq 6\times 10^{-7}\,\mathrm m$ である。

反射波と入射波は同じ真空中にあるため振幅比の二乗を取る。

$$
\begin{aligned}
R&=\Re r^{2} + \Im r^{2},\\
R&= \left(- \frac{2}{5}\right)^{2} + \left(- \frac{1}{10}\right)^{2}=\frac{17}{100}.
\end{aligned}
$$

複素分母を有理化して得られる二乗ノルムを求める。

$$
\begin{aligned}
D_r&=\Im r^{2} + \left(\Re r + 1\right)^{2},\\
D_r&= \left(- \frac{1}{10}\right)^{2} + \left(\left(- \frac{2}{5}\right) + 1\right)^{2}=\frac{37}{100}.
\end{aligned}
$$

有理化した式の実部を入れる。

$$
\begin{aligned}
n&=\frac{1 - R}{D_r},\\
n&= \frac{1 - \left(\frac{17}{100}\right)}{\left(\frac{37}{100}\right)}=\frac{83}{37}.
\end{aligned}
$$

同じ式の虚部を入れる。負の反射位相が正の減衰へつながる。

$$
\begin{aligned}
\kappa_n&=- \frac{2 \Im r}{D_r},\\
\kappa_n&= - \frac{2 \cdot \left(- \frac{1}{10}\right)}{\left(\frac{37}{100}\right)}=\frac{20}{37}.
\end{aligned}
$$

屈折率の虚部から、媒質内でパワーが $1/e$ になる深さを求める。

$$
\begin{aligned}
\ell_I&=\frac{\lambda_0}{4 \pi \kappa_n},\\
\ell_I/(\mathrm m)&\simeq \frac{\left(6\times 10^{-7}\right)}{4 \cdot \pi \cdot \left(\frac{20}{37}\right)}\simeq 8.8331\times 10^{-8}.
\end{aligned}
$$

複素振幅を微分し、その絶対値で一次誤差の増幅を測る。

$$
\begin{aligned}
\lvert\dd\tilde n/\dd r\rvert&=\frac{2}{D_r},\\
\lvert\dd\tilde n/\dd r\rvert&= \frac{2}{\left(\frac{37}{100}\right)}=\frac{200}{37}.
\end{aligned}
$$

位相から正の消衰係数が得られ、反射光だけでは見えない侵入の深さを決められる。
:::

:::example[同じ程度の反射率で位相が違う表面]
反射振幅の実部と虚部の配分が逆に近い試料を測る。

入力は反射係数の実部 $\Re r=- \frac{1}{10}$、反射係数の虚部 $\Im r=- \frac{2}{5}$、真空波長 $\lambda_0\simeq 6\times 10^{-7}\,\mathrm m$ である。

反射波と入射波は同じ真空中にあるため振幅比の二乗を取る。

$$
\begin{aligned}
R&=\Re r^{2} + \Im r^{2},\\
R&= \left(- \frac{1}{10}\right)^{2} + \left(- \frac{2}{5}\right)^{2}=\frac{17}{100}.
\end{aligned}
$$

複素分母を有理化して得られる二乗ノルムを求める。

$$
\begin{aligned}
D_r&=\Im r^{2} + \left(\Re r + 1\right)^{2},\\
D_r&= \left(- \frac{2}{5}\right)^{2} + \left(\left(- \frac{1}{10}\right) + 1\right)^{2}=\frac{97}{100}.
\end{aligned}
$$

有理化した式の実部を入れる。

$$
\begin{aligned}
n&=\frac{1 - R}{D_r},\\
n&= \frac{1 - \left(\frac{17}{100}\right)}{\left(\frac{97}{100}\right)}=\frac{83}{97}.
\end{aligned}
$$

同じ式の虚部を入れる。負の反射位相が正の減衰へつながる。

$$
\begin{aligned}
\kappa_n&=- \frac{2 \Im r}{D_r},\\
\kappa_n&= - \frac{2 \cdot \left(- \frac{2}{5}\right)}{\left(\frac{97}{100}\right)}=\frac{80}{97}.
\end{aligned}
$$

屈折率の虚部から、媒質内でパワーが $1/e$ になる深さを求める。

$$
\begin{aligned}
\ell_I&=\frac{\lambda_0}{4 \pi \kappa_n},\\
\ell_I/(\mathrm m)&\simeq \frac{\left(6\times 10^{-7}\right)}{4 \cdot \pi \cdot \left(\frac{80}{97}\right)}\simeq 5.7893\times 10^{-8}.
\end{aligned}
$$

複素振幅を微分し、その絶対値で一次誤差の増幅を測る。

$$
\begin{aligned}
\lvert\dd\tilde n/\dd r\rvert&=\frac{2}{D_r},\\
\lvert\dd\tilde n/\dd r\rvert&= \frac{2}{\left(\frac{97}{100}\right)}=\frac{200}{97}.
\end{aligned}
$$

前例と $R$ は同じだが、屈折率の虚部と侵入長は異なる。反射率だけでは区別できない二例である。
:::

:::example[ほぼ完全反射する面]
反射係数が $-1$ に近い表面の測定を反転する。

入力は反射係数の実部 $\Re r=- \frac{19}{20}$、反射係数の虚部 $\Im r=- \frac{3}{100}$、真空波長 $\lambda_0\simeq 6\times 10^{-7}\,\mathrm m$ である。

反射波と入射波は同じ真空中にあるため振幅比の二乗を取る。

$$
\begin{aligned}
R&=\Re r^{2} + \Im r^{2},\\
R&= \left(- \frac{19}{20}\right)^{2} + \left(- \frac{3}{100}\right)^{2}=\frac{4517}{5000}.
\end{aligned}
$$

複素分母を有理化して得られる二乗ノルムを求める。

$$
\begin{aligned}
D_r&=\Im r^{2} + \left(\Re r + 1\right)^{2},\\
D_r&= \left(- \frac{3}{100}\right)^{2} + \left(\left(- \frac{19}{20}\right) + 1\right)^{2}=\frac{17}{5000}.
\end{aligned}
$$

有理化した式の実部を入れる。

$$
\begin{aligned}
n&=\frac{1 - R}{D_r},\\
n&= \frac{1 - \left(\frac{4517}{5000}\right)}{\left(\frac{17}{5000}\right)}=\frac{483}{17}.
\end{aligned}
$$

同じ式の虚部を入れる。負の反射位相が正の減衰へつながる。

$$
\begin{aligned}
\kappa_n&=- \frac{2 \Im r}{D_r},\\
\kappa_n&= - \frac{2 \cdot \left(- \frac{3}{100}\right)}{\left(\frac{17}{5000}\right)}=\frac{300}{17}.
\end{aligned}
$$

屈折率の虚部から、媒質内でパワーが $1/e$ になる深さを求める。

$$
\begin{aligned}
\ell_I&=\frac{\lambda_0}{4 \pi \kappa_n},\\
\ell_I/(\mathrm m)&\simeq \frac{\left(6\times 10^{-7}\right)}{4 \cdot \pi \cdot \left(\frac{300}{17}\right)}\simeq 2.7056\times 10^{-9}.
\end{aligned}
$$

複素振幅を微分し、その絶対値で一次誤差の増幅を測る。

$$
\begin{aligned}
\lvert\dd\tilde n/\dd r\rvert&=\frac{2}{D_r},\\
\lvert\dd\tilde n/\dd r\rvert&= \frac{2}{\left(\frac{17}{5000}\right)}=\frac{10000}{17}.
\end{aligned}
$$

微小誤差の倍率は大きい。この場合は基準面を精密に決めてから光学定数へ換算する。
:::

## 演習問題

::::exercise{level=basic}
$n_1=1,n_2=2,\theta_i=30^\circ$ で $\sin\theta_t$ を求めよ。

:::solution
:eqref[eq-optics-snell] から $\sin\theta_t=(1/2)(1/2)=1/4$ である。
:::
::::

::::exercise{level=basic}
真空から非磁性で $n=2$ の媒質へ垂直入射するとき $R,T$ を求めよ。

:::solution
$Y_2/Y_1=2$ なので $r=-1/3,t=2/3$。したがって $R=1/9,T=2(2/3)^2=8/9$ である。
:::
::::

::::exercise{level=basic}
$n_1=2,n_2=1$ の臨界角を求めよ。

:::solution
透過角 $90^\circ$ での $\sin\theta_c=n_2/n_1=1/2$ より、$\theta_c=30^\circ$ である。
:::
::::

::::exercise{level=basic}
$\tilde n=2+i$ の媒質への真空からの垂直入射で反射率を求めよ。

:::solution
:eqref[eq-optics-metal] の分子は $1+1=2$、分母は $9+1=10$ であり、$R=1/5$ となる。
:::
::::

::::exercise{level=advanced}
全反射の $s$ 波で界面の直上の振幅が $t E_0$ である。$z=\kappa^{-1}$ で $\langle S_x(z)\rangle/\langle S_x(0)\rangle$ を求めよ。

:::solution
$H_z$ は $E_y$ に実定数 $k_x/(\omega\mu_2)$ を掛けた量であるから平均接線流は $|E_y|^2$ に比例する。:eqref[eq-optics-evanescent] で評価して $e^{-2}$ を得る。
:::
::::

::::exercise{level=advanced}
両側が伝播域となる接線波数を保ち、無損失の界面に逆向きから入射する。接線電場の規約で、反射振幅が元の負値になることを示せ。

:::solution
同じ接線波数で両側が伝播域なら、入射側と出射側を逆にすると :eqref[eq-optics-fresnel] の二つの正の $Y$ が入れ替わる。分母は同じで分子の符号だけが変わる。
:::
::::

::::exercise{level=advanced}
$D_s$ を $\theta_i$ に関する反射位相の微分で書け。

:::solution
$k_x=k_1\sin\theta_i$ だから $\dd k_x/\dd\theta_i=k_1\cos\theta_i$。連鎖律を :eqref[eq-optics-gh] に用いると $D_s=-(k_1\cos\theta_i)^{-1}\dd\phi_s/\dd\theta_i$ を得る。
:::
::::

::::exercise{level=advanced}
$n_1<n_2$ の非磁性媒質で Brewster 角入射の $r_s$ を屈折率だけで表せ。

:::solution
$\cos\theta_t=\sin\theta_B$ より $Y_{2s}/Y_{1s}=n_2\tan\theta_B/n_1=n_2^2/n_1^2$。従って $r_s=(n_1^2-n_2^2)/(n_1^2+n_2^2)$ である。
:::
::::

::::exercise{level=advanced}
真空から半無限の非磁性媒質に垂直入射すると、界面基準の電場反射係数が $r=-(1+i)/2$ だった。複素屈折率と吸収されるパワーの割合を求めよ。

:::solution
$1+r=(1-i)/2$、$1-r=(3+i)/2$ なので、共役を掛けて $\tilde n=(3+i)/(1-i)=1+2i$ を得る。虚部が正なので透過場は内部で減衰する。反射パワーの比は $|r|^2=1/2$。媒質は半無限で背面への出力がないため、内部へ入る $1-|r|^2=1/2$ が吸収される。
:::
::::

## 確認テスト

::::quiz
静止境界で入射波と透過波に共通の量をすべて選べ。

- [x] 周波数
- [x] 接線波数
- [ ] 法線波数

:::explanation
境界で任意の時刻と接線位置について位相をそろえるので、周波数と接線波数が共通になる。法線波数は各媒質で分散式から決める。
:::
::::

::::quiz{answer="0.25" tol="1e-6"}
垂直入射でアドミタンス比が $Y_2/Y_1=3$ のとき、パワー反射率を答えよ。

:::explanation
振幅反射係数は $(1-3)/(1+3)=-1/2$ なので、同じ媒質へ戻るパワーの割合はその絶対値の二乗で $1/4$ である。
:::
::::

::::quiz
全反射で生じるエバネッセント波について正しいものを選べ。

- [x] 法線の時間平均流はゼロである。
- [ ] 界面上の電場はゼロである。
- [ ] 接線の位相は一定である。

:::explanation
法線波数が純虚数なので $\Re Y_2=0$ である。一方、$t$ は有限であり、位相には $k_x x$ が残る。
:::
::::

::::quiz{answer="2" tol="1e-6"}
$n_1=n_2/2$ で Brewster 角の正接を答えよ。

:::explanation
:eqref[eq-optics-brewster] で $p$ のアドミタンスを整合させる角は $\tan\theta_B=n_2/n_1=2$ を満たす。
:::
::::

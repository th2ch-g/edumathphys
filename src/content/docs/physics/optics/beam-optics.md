---
title: ガウスビームと光共振器
description: 近軸波動方程式からガウスモードを導き、ABCD 則、共振器の安定性、共鳴周波数と損失を計算する。
notation:
  - symbol: 'w,w_0'
    meaning: 強度の 1/e^2 半径とその最小値
  - symbol: 'z_R'
    meaning: Rayleigh 長
  - symbol: 'q'
    meaning: 負の虚部を持つ複素ビームパラメータ
  - symbol: 'g_j'
    meaning: 鏡 j の無次元安定性パラメータ
  - symbol: '\tau_{\mathrm{ph}}'
    meaning: 共振器内エネルギーの 1/e 寿命
---

## この章の目標

- 近軸近似の尺度を確認し、基本モードの幅と位相を求める。
- 直交する高次モードの次数と Gouy 位相を結ぶ。
- ABCD 則でビームと共振器の自己整合を計算する。
- 縦・横モード間隔と光子寿命を見積もる。

## 近軸 Helmholtz 方程式

### 光軸に沿った速い位相を分離する
一様無損失媒質の Helmholtz 方程式に $U=A(x,y,z)e^{ikz}$ を代入すると、$\nabla_\perp^2 A+2ik\partial_zA+\partial_z^2A=0$。縦方向に包絡線が変わる尺度が $k^{-1}$ より十分長いとき最後の項を無視できる。これが :term[近軸 Helmholtz 方程式]{en="paraxial Helmholtz equation" yomi="きんじくへるむほるつほうていしき"} である。
$$
2ik\partial_z A+\nabla_\perp^2 A=0 .
\label{eq-optics-paraxial}
$$
横 Fourier 変換で各成分は $e^{-ik_\perp^2 z/(2k)}$ を受ける。これは式 :eqref[eq-optics-angular-spectrum] を $k_\perp/k\ll1$ で展開しても得られる。距離が長いと位相誤差が蓄積するため、$k_\perp^4 z/(8k^3)\ll1$ も確認する。場が横端で十分減衰すれば、式に $A^*$ を掛けて共役の式との差を積分することで、$\int|A|^2\dd^2\rho$ が一定と分かる。Laplace 項の差は $\divg(A^*\grad A-A\grad A^*)$ の積分で、境界項はゼロである。

:::example[傾いた成分]
$A=e^{ik_xx-ik_x^2z/(2k)}$ は式 :eqref[eq-optics-paraxial] を厳密に満たす。$k_x/k=0.1$ での厳密な縦波数比は $\sqrt{0.99}\simeq0.994987$、近軸式は $0.995$ である。
:::
:::example[平行な成分]
横座標に依存しない $A$ は近軸方程式で一定であり、平面波を再現する。ただし横面の積分は無限大なので、有限パワーの保存則をそのまま適用できない。
:::

### 近軸の核に必要な位相精度
最大傾斜の成分を $k_\perp/k=\epsilon$ とすると、縦波数の厳密値は $k(1-\epsilon^2/2-\epsilon^4/8+\cdots)$。距離 $L$ の最後の項の位相は $kL\epsilon^4/8$ である。例えば $\epsilon=0.1$、$\lambda=500\,\mathrm{nm}$ で $L=1\,\mathrm{mm}$ なら約 $0.16\,\mathrm{rad}$。近軸の縦波数の相対誤差だけを見ると約 $10^{-5}$ でも、位相の精度が必要な距離では積を取って確かめる。

:::example[高い角成分だけ近似から外れる]
$k_x/k=1/2$ の成分で厳密な $k_z/k$ は $\sqrt3/2$、近軸の二次式は $7/8$ である。相対差は一パーセントほどでも $z=20\lambda$ で予測位相が約 $1.13\,\mathrm{rad}$ ずれる。高開口の集光では高角成分を落とせない。
:::

### 傾いた近軸解とビーム中心の移動

$A(x,y,z)$ が $2ikA_z+\nabla_\perp^2A=0$ の解なら、実定数 $K$ に対し $B=e^{iKx-iK^2z/(2k)}A(x-Kz/k,y,z)$ も解である。実際、$2ik\partial_z B$ の移動座標から生じる $-2iK A_x$ が $\partial_x^2 B$ の $2iK A_x$ と打ち消し合い、二次の位相からの $K^2 A$ と $-K^2 A$ も消える。従って小さい横位相勾配を入れると包絡線は傾き $K/k$ で平行移動し、形の発展は元の解が与える。厳密な Helmholtz の中心成分では傾きは $K/\sqrt{k^2-K^2}$ である。多くの横成分を含むビームでは傾きだけでなく、成分ごとの位相の誤差も確認する。

:::example[小さい傾きの短い伝播]
中心波数に小さな横成分を加える。

入力は波長 $\lambda\simeq 1\times 10^{-6}\,\mathrm m$、加えた横波数 $K=100000\,\mathrm{m^{-1}}$、伝播距離 $z=\frac{1}{50}\,\mathrm m$ である。

媒質中の位相速度から波数を固定する。

$$
\begin{aligned}
k&=\frac{2 \pi}{\lambda},\\
k/(\mathrm{m^{-1}})&\simeq \frac{2 \cdot \pi}{\left(1\times 10^{-6}\right)}\simeq 6.2832\times 10^{6}.
\end{aligned}
$$

包絡線の移動速度を縦距離に対する比で求める。

$$
\begin{aligned}
a&=\frac{K}{k},\\
a&\simeq \frac{\left(100000\right)}{\left(6.283185307\times 10^{6}\right)}\simeq 0.015915.
\end{aligned}
$$

傾きを伝播距離へ掛ける。

$$
\begin{aligned}
x_{\rm p}&=a z,\\
x_{\rm p}/(\mathrm m)&\simeq \left(0.01591549431\right) \cdot \left(\frac{1}{50}\right)\simeq 3.1831\times 10^{-4}.
\end{aligned}
$$

厳密な角の正接を使う。

$$
\begin{aligned}
x_{\rm e}&=\frac{a z}{\sqrt{1 - a^{2}}},\\
x_{\rm e}/(\mathrm m)&\simeq \frac{\left(0.01591549431\right) \cdot \left(\frac{1}{50}\right)}{\sqrt{1 - \left(0.01591549431\right)^{2}}}\simeq 3.1835\times 10^{-4}.
\end{aligned}
$$

近軸の位置と厳密な中心波数からの位置を比べる。

$$
\begin{aligned}
\Delta x&=x_{\rm e} - x_{\rm p},\\
\Delta x/(\mathrm m)&\simeq \left(3.183502083\times 10^{-4}\right) - \left(3.183098862\times 10^{-4}\right)\simeq 4.0322\times 10^{-8}.
\end{aligned}
$$

短い伝播では位置の差は小さい。
:::

:::example[長い伝播で位置を読む]
同じ横成分を保ち、距離だけを増やす。

入力は波長 $\lambda\simeq 1\times 10^{-6}\,\mathrm m$、加えた横波数 $K=100000\,\mathrm{m^{-1}}$、伝播距離 $z=20\,\mathrm m$ である。

媒質中の位相速度から波数を固定する。

$$
\begin{aligned}
k&=\frac{2 \pi}{\lambda},\\
k/(\mathrm{m^{-1}})&\simeq \frac{2 \cdot \pi}{\left(1\times 10^{-6}\right)}\simeq 6.2832\times 10^{6}.
\end{aligned}
$$

包絡線の移動速度を縦距離に対する比で求める。

$$
\begin{aligned}
a&=\frac{K}{k},\\
a&\simeq \frac{\left(100000\right)}{\left(6.283185307\times 10^{6}\right)}\simeq 0.015915.
\end{aligned}
$$

傾きを伝播距離へ掛ける。

$$
\begin{aligned}
x_{\rm p}&=a z,\\
x_{\rm p}/(\mathrm m)&\simeq \left(0.01591549431\right) \cdot \left(20\right)\simeq 0.31831.
\end{aligned}
$$

厳密な角の正接を使う。

$$
\begin{aligned}
x_{\rm e}&=\frac{a z}{\sqrt{1 - a^{2}}},\\
x_{\rm e}/(\mathrm m)&\simeq \frac{\left(0.01591549431\right) \cdot \left(20\right)}{\sqrt{1 - \left(0.01591549431\right)^{2}}}\simeq 0.31835.
\end{aligned}
$$

近軸の位置と厳密な中心波数からの位置を比べる。

$$
\begin{aligned}
\Delta x&=x_{\rm e} - x_{\rm p},\\
\Delta x/(\mathrm m)&\simeq \left(0.3183502083\right) - \left(0.3183098862\right)\simeq 4.0322\times 10^{-5}.
\end{aligned}
$$

小さな角度差も距離に比例して位置へ現れる。
:::

:::example[大きな横波数を与える]
同じ短距離で傾きを増し、ずれを比較する。

入力は波長 $\lambda\simeq 1\times 10^{-6}\,\mathrm m$、加えた横波数 $K=1000000\,\mathrm{m^{-1}}$、伝播距離 $z=\frac{1}{50}\,\mathrm m$ である。

媒質中の位相速度から波数を固定する。

$$
\begin{aligned}
k&=\frac{2 \pi}{\lambda},\\
k/(\mathrm{m^{-1}})&\simeq \frac{2 \cdot \pi}{\left(1\times 10^{-6}\right)}\simeq 6.2832\times 10^{6}.
\end{aligned}
$$

包絡線の移動速度を縦距離に対する比で求める。

$$
\begin{aligned}
a&=\frac{K}{k},\\
a&\simeq \frac{\left(1000000\right)}{\left(6.283185307\times 10^{6}\right)}\simeq 0.15915.
\end{aligned}
$$

傾きを伝播距離へ掛ける。

$$
\begin{aligned}
x_{\rm p}&=a z,\\
x_{\rm p}/(\mathrm m)&\simeq \left(0.1591549431\right) \cdot \left(\frac{1}{50}\right)\simeq 0.0031831.
\end{aligned}
$$

厳密な角の正接を使う。

$$
\begin{aligned}
x_{\rm e}&=\frac{a z}{\sqrt{1 - a^{2}}},\\
x_{\rm e}/(\mathrm m)&\simeq \frac{\left(0.1591549431\right) \cdot \left(\frac{1}{50}\right)}{\sqrt{1 - \left(0.1591549431\right)^{2}}}\simeq 0.0032242.
\end{aligned}
$$

近軸の位置と厳密な中心波数からの位置を比べる。

$$
\begin{aligned}
\Delta x&=x_{\rm e} - x_{\rm p},\\
\Delta x/(\mathrm m)&\simeq \left(0.003224195695\right) - \left(0.003183098862\right)\simeq 4.1097\times 10^{-5}.
\end{aligned}
$$

位置ずれは傾きの三次から始まるため、距離が短くても増える。
:::

## ガウスビームと Rayleigh 長・Gouy 位相

### 二次位相を持つ解
半径を $r=\sqrt{x^2+y^2}$ とし、$A=C(z)\exp[ik r^2/(2q(z))]$ と試す。横 Laplace 演算の結果は $(2ik/q-k^2r^2/q^2)A$。縦微分も代入して $r^2$ と定数の係数を比較すると $q'=1,\ C'/C=-1/q$ となる。腰を $z=0$ とし $A(r,0)=A_0e^{-r^2/w_0^2}$ とすると $q_0=-iz_R$、$z_R=kw_0^2/2=\pi w_0^2/\lambda$ である。従って
$$
A=A_0\frac{q_0}{q}e^{ik r^2/(2q)},\qquad q=z-iz_R,\qquad
\frac1q=\frac1{R(z)}+\frac{2i}{kw(z)^2}.
\label{eq-optics-gaussian}
$$
これは :term[ガウスビーム]{en="Gaussian beam" yomi="がうすびーむ"} の基本モードである。実部と虚部を分ければ $w=w_0\sqrt{1+(z/z_R)^2}$、$R=(z^2+z_R^2)/z$。$w$ は強度の $1/e^2$ 半径であり、遠方半角は $\theta=w_0/z_R=\lambda/(\pi w_0)$ である。:term[Rayleigh 長]{en="Rayleigh range" yomi="れいりーちょう"} $z_R$ で半径は $\sqrt2$ 倍となる。前因子 $q_0/q=(1+iz/z_R)^{-1}$ の位相 $-\arctan(z/z_R)$ が :term[Gouy 位相]{en="Gouy phase" yomi="ぐいーいそう"} であり、平面波に対する位相の遅れを表す。$w^2$ の増加と軸上強度の減少が相殺し、パワーは一定である。

:::example[ビームの発散]
$\lambda=1.0\,\mathrm{\mu m}$、$w_0=100\,\mathrm{\mu m}$ なら $z_R\simeq31.4\,\mathrm{mm}$、半角は約 $3.18\,\mathrm{mrad}$。
:::
:::example[腰の前後]
$z=\pm z_R$ の幅は同じになるが、曲率半径は $R=\pm2z_R$ である。収束中の波と発散中の波は同じ強度分布でも曲率の符号で区別できる。
:::

### 腰を直接見ずに求める
ある面で幅 $w$ と波面の曲率半径 $R$ を測れば :eqref[eq-optics-gaussian] の $1/q$ が決まる。腰からの距離は $\Re q$、Rayleigh 長は $-\Im q$ であり、$w_0=\sqrt{\lambda z_R/\pi}$ である。$R$ の符号も測れば、同じ幅の収束側と発散側を区別できる。

:::example[同じ二つの項を持つ逆数]
測定面で $1/R=2/(kw^2)=b>0$ なら $1/q=b(1+i)$、従って $q=(1-i)/(2b)$。腰は測定面の手前 $1/(2b)$ にあり、そこから測定面までの距離が Rayleigh 長に等しい。測定幅は腰の幅の $\sqrt2$ 倍である。
:::

### 幅と発散が下げ合う理由
腰の一次元断面の振幅は $e^{-x^2/w_0^2}$。横波数の振幅は Gauss 積分で $e^{-k_x^2w_0^2/4}$ に比例する。従って横位置の強度の標準偏差は $w_0/2$、波数強度の標準偏差は $1/w_0$ である。その積は $1/2$。一般の正規化振幅 $u$ について平均位置と平均波数をゼロに直すと、$\int x\,\partial_x|u|^2\dd x=-1$。積分の実部の絶対値を Cauchy–Schwarz で抑えると $1\le2\|x u\|_2\|\partial_xu\|_2$ となる。Fourier 側で $\|\partial_xu\|_2$ は波数の標準偏差なので、基本モードはこの幅と波数の下限を達成する。有限開口で形を切れば、単純に腰を狭めるだけでは遠方の発散も同時に小さくできない。

### 強度半径と受光パワーを同時に確認する

円対称のガウスビームで全パワーが $P$ なら、軸上照度は $I_0=2P/(\pi w^2)$ である。強度 $I_0e^{-2r^2/w^2}$ を $2\pi r\dd r$ で積分すれば全パワー $P$ が戻り、半径 $a$ の検出器に入るパワーは $P[1-e^{-2a^2/w^2}]$。従って伝播して光束のパワーが保存されても、固定径の受光器で測るパワーは減り得る。検出器が受光領域を覆い、入射位置による感度差と開口の透過率を校正したとき、二つの受光半径のパワー比からビーム幅を検証できる。

:::example[腰をそのまま受光する]
強度半径と等しい検出半径を使う。

入力は波長 $\lambda\simeq 1\times 10^{-6}\,\mathrm m$、腰の半径 $w_0\simeq 1\times 10^{-4}\,\mathrm m$、腰からの距離 $z=0\,\mathrm m$、全パワー $P=\frac{1}{100}\,\mathrm W$、受光器の半径 $a\simeq 1\times 10^{-4}\,\mathrm m$ である。

腰の波面曲率が変わる軸方向尺度を求める。

$$
\begin{aligned}
z_R&=\frac{\pi w_0^{2}}{\lambda},\\
z_R/(\mathrm m)&\simeq \frac{\pi \cdot \left(\left(1\times 10^{-4}\right)\right)^{2}}{\left(1\times 10^{-6}\right)}\simeq 0.031416.
\end{aligned}
$$

距離に応じたビームの広がりを求める。

$$
\begin{aligned}
w(z)&=w_0 \sqrt{1 + \frac{z^{2}}{z_R^{2}}},\\
w(z)/(\mathrm m)&\simeq \left(1\times 10^{-4}\right) \cdot \sqrt{1 + \frac{\left(0\right)^{2}}{\left(0.03141592654\right)^{2}}}\simeq 1\times 10^{-4}.
\end{aligned}
$$

軸上で平面波に対する位相を計算する。

$$
\begin{aligned}
\phi_G&=- \operatorname{atan}{\left(\frac{z}{z_R} \right)},\\
\phi_G/(\mathrm{rad})&\simeq - \operatorname{atan}{\left(\frac{\left(0\right)}{\left(0.03141592654\right)} \right)}=0.
\end{aligned}
$$

円の面積に対するガウスの規格化を取る。

$$
\begin{aligned}
I_0&=\frac{2 P}{\pi w(z)^{2}},\\
I_0/(\mathrm{W/m^2})&\simeq \frac{2 \cdot \left(\frac{1}{100}\right)}{\pi \cdot \left(\left(1\times 10^{-4}\right)\right)^{2}}\simeq 6.3662\times 10^{5}.
\end{aligned}
$$

指数強度を受光器の円内だけで積分する。

$$
\begin{aligned}
P_a&=P \left(1 - e^{- \frac{2 a^{2}}{w(z)^{2}}}\right),\\
P_a/(\mathrm W)&\simeq \left(\frac{1}{100}\right) \cdot \left(1 - e^{- \frac{2 \cdot \left(\left(1\times 10^{-4}\right)\right)^{2}}{\left(\left(1\times 10^{-4}\right)\right)^{2}}}\right)\simeq 0.0086466.
\end{aligned}
$$

全パワーを取るには、検出半径を腰より大きくする必要がある。
:::

:::example[Rayleigh 長の先で受光する]
同じ受光器を腰から Rayleigh 長だけ離す。

入力は波長 $\lambda\simeq 1\times 10^{-6}\,\mathrm m$、腰の半径 $w_0\simeq 1\times 10^{-4}\,\mathrm m$、腰からの距離 $z\simeq 0.031416\,\mathrm m$、全パワー $P=\frac{1}{100}\,\mathrm W$、受光器の半径 $a\simeq 1\times 10^{-4}\,\mathrm m$ である。

腰の波面曲率が変わる軸方向尺度を求める。

$$
\begin{aligned}
z_R&=\frac{\pi w_0^{2}}{\lambda},\\
z_R/(\mathrm m)&\simeq \frac{\pi \cdot \left(\left(1\times 10^{-4}\right)\right)^{2}}{\left(1\times 10^{-6}\right)}\simeq 0.031416.
\end{aligned}
$$

距離に応じたビームの広がりを求める。

$$
\begin{aligned}
w(z)&=w_0 \sqrt{1 + \frac{z^{2}}{z_R^{2}}},\\
w(z)/(\mathrm m)&\simeq \left(1\times 10^{-4}\right) \cdot \sqrt{1 + \frac{\left(0.03141592654\right)^{2}}{\left(0.03141592654\right)^{2}}}\simeq 1.4142\times 10^{-4}.
\end{aligned}
$$

軸上で平面波に対する位相を計算する。

$$
\begin{aligned}
\phi_G&=- \operatorname{atan}{\left(\frac{z}{z_R} \right)},\\
\phi_G/(\mathrm{rad})&\simeq - \operatorname{atan}{\left(\frac{\left(0.03141592654\right)}{\left(0.03141592654\right)} \right)}\simeq -0.7854.
\end{aligned}
$$

円の面積に対するガウスの規格化を取る。

$$
\begin{aligned}
I_0&=\frac{2 P}{\pi w(z)^{2}},\\
I_0/(\mathrm{W/m^2})&\simeq \frac{2 \cdot \left(\frac{1}{100}\right)}{\pi \cdot \left(\left(1.414213562\times 10^{-4}\right)\right)^{2}}\simeq 3.1831\times 10^{5}.
\end{aligned}
$$

指数強度を受光器の円内だけで積分する。

$$
\begin{aligned}
P_a&=P \left(1 - e^{- \frac{2 a^{2}}{w(z)^{2}}}\right),\\
P_a/(\mathrm W)&\simeq \left(\frac{1}{100}\right) \cdot \left(1 - e^{- \frac{2 \cdot \left(\left(1\times 10^{-4}\right)\right)^{2}}{\left(\left(1.414213562\times 10^{-4}\right)\right)^{2}}}\right)\simeq 0.0063212.
\end{aligned}
$$

軸上照度が半分になり、受光パワーも減る。
:::

:::example[受光径を増す]
前例の位置で受光半径だけを $\sqrt2$ 倍にする。

入力は波長 $\lambda\simeq 1\times 10^{-6}\,\mathrm m$、腰の半径 $w_0\simeq 1\times 10^{-4}\,\mathrm m$、腰からの距離 $z\simeq 0.031416\,\mathrm m$、全パワー $P=\frac{1}{100}\,\mathrm W$、受光器の半径 $a\simeq 1.4142\times 10^{-4}\,\mathrm m$ である。

腰の波面曲率が変わる軸方向尺度を求める。

$$
\begin{aligned}
z_R&=\frac{\pi w_0^{2}}{\lambda},\\
z_R/(\mathrm m)&\simeq \frac{\pi \cdot \left(\left(1\times 10^{-4}\right)\right)^{2}}{\left(1\times 10^{-6}\right)}\simeq 0.031416.
\end{aligned}
$$

距離に応じたビームの広がりを求める。

$$
\begin{aligned}
w(z)&=w_0 \sqrt{1 + \frac{z^{2}}{z_R^{2}}},\\
w(z)/(\mathrm m)&\simeq \left(1\times 10^{-4}\right) \cdot \sqrt{1 + \frac{\left(0.03141592654\right)^{2}}{\left(0.03141592654\right)^{2}}}\simeq 1.4142\times 10^{-4}.
\end{aligned}
$$

軸上で平面波に対する位相を計算する。

$$
\begin{aligned}
\phi_G&=- \operatorname{atan}{\left(\frac{z}{z_R} \right)},\\
\phi_G/(\mathrm{rad})&\simeq - \operatorname{atan}{\left(\frac{\left(0.03141592654\right)}{\left(0.03141592654\right)} \right)}\simeq -0.7854.
\end{aligned}
$$

円の面積に対するガウスの規格化を取る。

$$
\begin{aligned}
I_0&=\frac{2 P}{\pi w(z)^{2}},\\
I_0/(\mathrm{W/m^2})&\simeq \frac{2 \cdot \left(\frac{1}{100}\right)}{\pi \cdot \left(\left(1.414213562\times 10^{-4}\right)\right)^{2}}\simeq 3.1831\times 10^{5}.
\end{aligned}
$$

指数強度を受光器の円内だけで積分する。

$$
\begin{aligned}
P_a&=P \left(1 - e^{- \frac{2 a^{2}}{w(z)^{2}}}\right),\\
P_a/(\mathrm W)&\simeq \left(\frac{1}{100}\right) \cdot \left(1 - e^{- \frac{2 \cdot \left(\left(1.414213562\times 10^{-4}\right)\right)^{2}}{\left(\left(1.414213562\times 10^{-4}\right)\right)^{2}}}\right)\simeq 0.0086466.
\end{aligned}
$$

軸上照度は変わらないが、受光パワーの割合は腰での値に戻る。
:::

## 高次モード（Hermite–Gauss・Laguerre–Gauss）

### 伸縮する無次元変数で分離する
$\xi=\sqrt2x/w,\ \eta=\sqrt2 y/w$、$\psi=\arctan(z/z_R)$ とする。基本解に多項式 $F(\xi,\eta)$ と位相 $e^{-iN\psi}$ を掛ける。近軸式に代入すると、$\xi_z=-(w'/w)\xi$、$w'/w=1/R$ のため変数の伸縮に由来する微分は曲率位相の勾配項と打ち消し合う。残る項を $2/w^2$ で割ると
$$
(\partial_\xi^2+\partial_\eta^2-2\xi\partial_\xi-2\eta\partial_\eta)F+2N F=0
$$
になる。ここで $k\psi'=2/w^2$ を使った。多項式 $H_j$ を Rodrigues の式 $H_j(\xi)=(-1)^j e^{\xi^2}(\dd/\dd\xi)^j e^{-\xi^2}$ で定める。$e^{2t\xi-t^2}=\sum_j H_j(\xi)t^j/j!$ は Taylor 展開から得られ、この母関数へ $\partial_\xi^2-2\xi\partial_\xi+2t\partial_t$ を作用させるとゼロである。係数を比べて $F=H_m(\xi)H_n(\eta),\ N=m+n$ を得る。これが :term[Hermite–Gauss モード]{en="Hermite–Gaussian mode" yomi="えるみーとがうすもーど"} である。二つの母関数の積を $e^{-\xi^2}$ で積分すると $\sqrt\pi e^{2tt'}$ であり、係数比較で異なる次数は直交する。

### 円筒座標の高次解
$s=\xi^2+\eta^2$、方位角を $\varphi$ とし、$F=s^{|\ell|/2}e^{i\ell\varphi}L(s)$ とおく。横演算子の半径部分は $4s\partial_s^2+4\partial_s-4s\partial_s$、角度部分は $\partial_\varphi^2/s$ である。代入して
$$
sL''+(|\ell|+1-s)L'+\frac{N-|\ell|}{2}L=0
$$
を得る。$L=\sum_j a_js^j$ とおくと $(j+1)(j+|\ell|+1)a_{j+1}=[j-(N-|\ell|)/2]a_j$。$N=2p+|\ell|$、$p$ が非負整数なら級数が次数 $p$ で止まる。これが一般化 Laguerre 多項式で、対応するモードを :term[Laguerre–Gauss モード]{en="Laguerre–Gaussian mode" yomi="らげーるがうすもーど"} という。$\ell$ が違えば角度積分で直交する。同じ $\ell$ では微分方程式を $(s^{|\ell|+1}e^{-s}L')'=-p s^{|\ell|}e^{-s}L$ と書き、異なる $p$ の二式を交差して引き、部分積分すれば直交する。

:::example[同じ横次数]
$H_1(\xi)H_0(\eta)=2\xi$ と $p=0,\ell=1$ の $\sqrt{s}e^{i\varphi}=\xi+i\eta$ は共に $N=1$。基本解のものも含めて Gouy 位相は $-2\psi$ になる。
:::
:::example[軸上の節]
$\ell=2,p=0$ では振幅に $r^2$ が掛かるため、軸上の強度はゼロである。位相は軸の周りを一周すると $4\pi$ 変化する。
:::

### 高次モードの二次モーメント
Hermite–Gauss の腰で $W_m(\xi)=H_m(\xi)^2e^{-\xi^2}$ を重みとする。母関数から $\xi H_m=H_{m+1}/2+mH_{m-1}$ を得る。右辺の二項は直交し、それぞれのノルム二乗は $\sqrt\pi\,2^j j!$ である。従って $\langle\xi^2\rangle_m=(2m+1)/2$。横実座標の標準偏差は $(w_0/2)\sqrt{2m+1}$ になる。伝播しても同じ規格化横形状なので、遠方での角の標準偏差も基本モードの $\sqrt{2m+1}$ 倍、幅と角の積は $2m+1$ 倍である。これが各方向で二次モーメントに基づくビーム品質係数 $M_x^2$ の高次モードでの値になる。

:::example[最低横次数の縮退を基底変換する]
腰で規格化した $u_{10},u_{01}$ は同じ次数なので、伝播位相と Gouy 位相が等しい。線形結合 $(u_{10}\pm i u_{01})/\sqrt2$ は軸対称な絶対値と $e^{\pm i\varphi}$ の波面を持つ。従って同じ二次元部分空間を、二つの Cartesian の節線で表すか、逆向きの位相渦の組で表すか選べる。非点収差で二方向の Gouy 位相が異なると、この結合は同じ形を保たない。
:::

### 半径方向の確率分布から渦の輪を測る

径方向次数 $p=0$ の Laguerre–Gauss モードで、$m=|\ell|$ とおく。局所強度は $r^{2m} e^{-2r^2/w^2}$ に比例するから、$m>0$ の最も明るい半径は対数を微分して $r_{\rm pk}=w\sqrt{m/2}$。一方、円環に入るパワーを測るときは面積因子 $r\dd r$ を掛ける。$s=2r^2/w^2$ における規格化密度は $s^me^{-s}/m!$ である。これは $\int_0^\infty s^me^{-s}\dd s=m!$ を部分積分で繰り返して確認できる。もう一つ $s$ を掛けて積分することにより $\langle r^2\rangle=w^2(m+1)/2$ を得る。局所の最大と円環での二乗平均は同じ半径ではない。

:::example[渦を持たない基本モード]
軸上に強度が集まる基準を求める。

入力は方位次数の絶対値 $m=0$、ガウス因子の半径 $w\simeq 1\times 10^{-4}\,\mathrm m$、腰からの距離 $z=\frac{3}{100}\,\mathrm m$、Rayleigh 長 $z_R=\frac{3}{100}\,\mathrm m$ である。

径方向の対数微分の零点を求める。

$$
\begin{aligned}
r_{\rm pk}&=\frac{\sqrt{2} \sqrt{m} w}{2},\\
r_{\rm pk}/(\mathrm m)&\simeq \frac{\sqrt{2} \cdot \sqrt{\left(0\right)} \cdot \left(1\times 10^{-4}\right)}{2}=0.
\end{aligned}
$$

円環の面積で重み付けした分布を使う。

$$
\begin{aligned}
r_{\rm rms}&=w \sqrt{\frac{m}{2} + \frac{1}{2}},\\
r_{\rm rms}/(\mathrm m)&\simeq \left(1\times 10^{-4}\right) \cdot \sqrt{\frac{\left(0\right)}{2} + \frac{1}{2}}\simeq 7.0711\times 10^{-5}.
\end{aligned}
$$

横次数が増すと基本モードの Gouy 位相に整数倍が掛かる。

$$
\begin{aligned}
\phi_G&=\left(- m - 1\right) \operatorname{atan}{\left(\frac{z}{z_R} \right)},\\
\phi_G/(\mathrm{rad})&= \left(- \left(0\right) - 1\right) \cdot \operatorname{atan}{\left(\frac{\left(\frac{3}{100}\right)}{\left(\frac{3}{100}\right)} \right)}\simeq -0.7854.
\end{aligned}
$$

方位角の指数因子から位相が何回転するかを読む。

$$
\begin{aligned}
\lvert\Delta\phi\rvert&=2 \pi m,\\
\lvert\Delta\phi\rvert/(\mathrm{rad})&= 2 \cdot \pi \cdot \left(0\right)=0.
\end{aligned}
$$

最大強度は軸上でも二乗平均半径は有限である。
:::

:::example[一周に一巻きの渦]
同じ腰からの距離とガウス尺度で一方位次数を入れる。

入力は方位次数の絶対値 $m=1$、ガウス因子の半径 $w\simeq 1\times 10^{-4}\,\mathrm m$、腰からの距離 $z=\frac{3}{100}\,\mathrm m$、Rayleigh 長 $z_R=\frac{3}{100}\,\mathrm m$ である。

径方向の対数微分の零点を求める。

$$
\begin{aligned}
r_{\rm pk}&=\frac{\sqrt{2} \sqrt{m} w}{2},\\
r_{\rm pk}/(\mathrm m)&\simeq \frac{\sqrt{2} \cdot \sqrt{\left(1\right)} \cdot \left(1\times 10^{-4}\right)}{2}\simeq 7.0711\times 10^{-5}.
\end{aligned}
$$

円環の面積で重み付けした分布を使う。

$$
\begin{aligned}
r_{\rm rms}&=w \sqrt{\frac{m}{2} + \frac{1}{2}},\\
r_{\rm rms}/(\mathrm m)&\simeq \left(1\times 10^{-4}\right) \cdot \sqrt{\frac{\left(1\right)}{2} + \frac{1}{2}}\simeq 1\times 10^{-4}.
\end{aligned}
$$

横次数が増すと基本モードの Gouy 位相に整数倍が掛かる。

$$
\begin{aligned}
\phi_G&=\left(- m - 1\right) \operatorname{atan}{\left(\frac{z}{z_R} \right)},\\
\phi_G/(\mathrm{rad})&= \left(- \left(1\right) - 1\right) \cdot \operatorname{atan}{\left(\frac{\left(\frac{3}{100}\right)}{\left(\frac{3}{100}\right)} \right)}\simeq -1.5708.
\end{aligned}
$$

方位角の指数因子から位相が何回転するかを読む。

$$
\begin{aligned}
\lvert\Delta\phi\rvert&=2 \pi m,\\
\lvert\Delta\phi\rvert/(\mathrm{rad})&= 2 \cdot \pi \cdot \left(1\right)\simeq 6.2832.
\end{aligned}
$$

局所強度の最大は輪上へ移る。
:::

:::example[三巻きの渦]
方位次数をさらに増やす。

入力は方位次数の絶対値 $m=3$、ガウス因子の半径 $w\simeq 1\times 10^{-4}\,\mathrm m$、腰からの距離 $z=\frac{3}{100}\,\mathrm m$、Rayleigh 長 $z_R=\frac{3}{100}\,\mathrm m$ である。

径方向の対数微分の零点を求める。

$$
\begin{aligned}
r_{\rm pk}&=\frac{\sqrt{2} \sqrt{m} w}{2},\\
r_{\rm pk}/(\mathrm m)&\simeq \frac{\sqrt{2} \cdot \sqrt{\left(3\right)} \cdot \left(1\times 10^{-4}\right)}{2}\simeq 1.2247\times 10^{-4}.
\end{aligned}
$$

円環の面積で重み付けした分布を使う。

$$
\begin{aligned}
r_{\rm rms}&=w \sqrt{\frac{m}{2} + \frac{1}{2}},\\
r_{\rm rms}/(\mathrm m)&\simeq \left(1\times 10^{-4}\right) \cdot \sqrt{\frac{\left(3\right)}{2} + \frac{1}{2}}\simeq 1.4142\times 10^{-4}.
\end{aligned}
$$

横次数が増すと基本モードの Gouy 位相に整数倍が掛かる。

$$
\begin{aligned}
\phi_G&=\left(- m - 1\right) \operatorname{atan}{\left(\frac{z}{z_R} \right)},\\
\phi_G/(\mathrm{rad})&= \left(- \left(3\right) - 1\right) \cdot \operatorname{atan}{\left(\frac{\left(\frac{3}{100}\right)}{\left(\frac{3}{100}\right)} \right)}\simeq -3.1416.
\end{aligned}
$$

方位角の指数因子から位相が何回転するかを読む。

$$
\begin{aligned}
\lvert\Delta\phi\rvert&=2 \pi m,\\
\lvert\Delta\phi\rvert/(\mathrm{rad})&= 2 \cdot \pi \cdot \left(3\right)\simeq 18.85.
\end{aligned}
$$

同じ $w$ でも輪が広がり、Gouy 位相の進行も変わる。
:::

## ABCD 則と q パラメータ

### 曲率と幅を一つの変数で伝える
ここからは一様な媒質内で光線座標を $(y,\theta)$ とする。式 :eqref[eq-optics-ray-matrix] の $(y,n\theta)$ 座標で表した行列とは、成分の規格化が異なる。自由伝播 $d$ は $q\mapsto q+d$。薄レンズは位相 $-k r^2/(2f)$ を加えるので $1/q\mapsto1/q-1/f$、すなわち $q\mapsto q/(1-q/f)$ である。この二つはそれぞれ光線行列の成分を使った分数線形写像である。写像を合成して分母を払えば行列の積と同じ成分が得られるので、
$$
q_{\mathrm{out}}=\frac{Aq_{\mathrm{in}}+B}{Cq_{\mathrm{in}}+D}
\label{eq-optics-abcq}
$$
が成立する。これを :term[ABCD 則]{en="ABCD law" yomi="えーびーしーでぃーそく"} という。$\det M=1$ の実行列では虚部が $\Im q_{\mathrm{in}}/|Cq_{\mathrm{in}}+D|^2$ となり、減衰するガウス解の枝が保たれる。波面のパラメータ $q$ は全体位相を含まないので、共振の計算では Gouy 位相を別に追う。

:::example[平面波面をレンズに置く]
入力 $q=-iz_R$ を焦点距離 $f>0$ のレンズに入れる。出力の逆数は $-1/f+i/z_R$。従ってレンズ直後の波面は収束側に曲がるが、直後の幅は変わらない。
:::
:::example[焦点位置]
この例で $q_{\rm out}=(-f z_R^2-i f^2 z_R)/(f^2+z_R^2)$。次の自由伝播で実部がゼロになる位置は $d=fz_R^2/(f^2+z_R^2)$ で、$z_R\gg f$ のとき幾何光学的な $f$ に近づく。
:::

### 二つの腰のモード整合
同軸で同じ面に腰を置いた二つの規格化基本モードを $u_j=\sqrt{2/\pi}\,w_j^{-1}e^{-r^2/w_j^2}$ とする。実の重なりは二次元 Gauss 積分で $2w_1w_2/(w_1^2+w_2^2)$。射影のパワーはその二乗である。幅が同じでも曲率が異なる場合は、被積分関数に $e^{ik r^2(1/R_1-1/R_2)/2}$ が掛かり、一般に重なりは 1 より小さい。ABCD 則で一方の複素 $q$ を他方に合わせると、曲率と幅の両条件を同時に満たせる。

:::example[腰の径が二倍違う場合]
$w_2=2w_1$ なら重なり振幅は $4/5$、パワーは $16/25$ となる。他モードへの分配と共振器鏡での損失は別に数える。
:::

:::example[角度がずれたビームの重なり]
同じ腰 $w$ の基本モードで片方に $e^{ik\theta x}$ の位相傾きを付ける。腰面の強度重みでこの指数の平均を取ると、振幅重なりは $e^{-(kw\theta)^2/8}$、パワー結合は $e^{-(kw\theta)^2/4}$。$\theta=\lambda/(\pi w)$、すなわち基本ビームの遠方半角だけ傾けるとパワー結合は $e^{-1}$ となる。
:::

### 屈折率が異なる界面の通過

媒質が変わるときは、光線座標 $(y,n\theta)$ に合わせて $Q=q/n$ を使う。真空波数 $k_0$ で式 :eqref[eq-optics-gaussian] を書き直すと $1/Q=n/R+2i/(k_0w^2)$ である。自由伝播距離 $d$ では $Q\mapsto Q+d/n$。平面界面では近軸の Snell 則で $n\theta$ が保存され、界面上のビーム半径も連続なので $Q$ は変わらない。

頂点の右に曲率中心を持つ球面を $R_s>0$ とする。界面の位置は軸近くで $z=y^2/(2R_s)$ であり、その面での入射側と出射側の位相を合わせると $n_2/R_{\rm out}-n_1/R_{\rm in}=-(n_2-n_1)/R_s$ となる。横幅の連続性から $1/Q_{\rm out}=1/Q_{\rm in}-(n_2-n_1)/R_s$ である。これは屈折面の光線行列の下左成分を $-(n_2-n_1)/R_s$ とした分数線形写像であり、伝播や屈折を同じ座標で合成できる。

:::example[平面界面の腰]
腰の面に平面界面があるとき、幅と $Q$ は連続である。空気から屈折率 $1.5$ の媒質に入るなら $q$ と Rayleigh 長は $1.5$ 倍になる。同じ物理的な腰幅でも、媒質内波長が短くなるので発散半角は $1/1.5$ 倍になる。
:::

### レンズが腰の位置と大きさを変える

入力の $q=x-iy$、$y>0$ を焦点距離 $f$ の薄レンズへ入れる。分数線形式に分母の共役を掛けると
$$
\Re q'=\frac{f[x(f-x)-y^2]}{(f-x)^2+y^2},\qquad
-\Im q'=\frac{f^2y}{(f-x)^2+y^2}.
$$
有限な実焦点距離 $f\ne0$ では第二式が正なので、出力は正の腰幅を持つガウスビームになる。レンズ直後から $d=-\Re q'$ の位置に新しい腰があり、$d<0$ なら出射場を逆伝播した像としての腰である。$x$ を正に選べば元の腰はレンズの手前にある。幾何光学での焦点からの差を求める場合、入力の腰とレンズの距離だけでなく $y$ も使う。

:::example[レンズ面が入力の腰]
入射波面は平面で、正レンズで曲率を与える。

入力は入力の腰からレンズの距離 $x=0\,\mathrm m$、入力 Rayleigh 長 $z_R=\frac{1}{25}\,\mathrm m$、レンズの焦点距離 $f=\frac{1}{50}\,\mathrm m$、波長 $\lambda\simeq 1\times 10^{-6}\,\mathrm m$ である。

複素分母にその共役を掛ける。

$$
\begin{aligned}
D&=z_R^{2} + \left(f - x\right)^{2},\\
D/(\mathrm{m^2})&= \left(\frac{1}{25}\right)^{2} + \left(\left(\frac{1}{50}\right) - \left(0\right)\right)^{2}=\frac{1}{500}.
\end{aligned}
$$

レンズ直後での $q'$ の実部を自由伝播で打ち消す。

$$
\begin{aligned}
d&=- \frac{f \left(x \left(f - x\right) - z_R^{2}\right)}{D},\\
d/(\mathrm m)&= - \frac{\left(\frac{1}{50}\right) \cdot \left(\left(0\right) \cdot \left(\left(\frac{1}{50}\right) - \left(0\right)\right) - \left(\frac{1}{25}\right)^{2}\right)}{\left(\frac{1}{500}\right)}=\frac{2}{125}.
\end{aligned}
$$

虚部の符号を変えて幅の尺度を求める。

$$
\begin{aligned}
z'_R&=\frac{f^{2} z_R}{D},\\
z'_R/(\mathrm m)&= \frac{\left(\frac{1}{50}\right)^{2} \cdot \left(\frac{1}{25}\right)}{\left(\frac{1}{500}\right)}=\frac{1}{125}.
\end{aligned}
$$

Rayleigh 長からガウスの半径へ戻す。

$$
\begin{aligned}
w'_0&=\frac{\sqrt{z'_R \lambda}}{\sqrt{\pi}},\\
w'_0/(\mathrm m)&\simeq \frac{\sqrt{\left(\frac{1}{125}\right) \cdot \left(1\times 10^{-6}\right)}}{\sqrt{\pi}}\simeq 5.0463\times 10^{-5}.
\end{aligned}
$$

同じ波長での二つの腰の大きさを比べる。

$$
\begin{aligned}
w'_0/w_0&=\sqrt{\frac{z'_R}{z_R}},\\
w'_0/w_0&= \sqrt{\frac{\left(\frac{1}{125}\right)}{\left(\frac{1}{25}\right)}}\simeq 0.44721.
\end{aligned}
$$

出力の腰はレンズの後ろにあり、位置は焦点距離より近い。
:::

:::example[入力の腰が前焦点面]
腰を前焦点へ移し、出力の腰を測る。

入力は入力の腰からレンズの距離 $x=\frac{1}{50}\,\mathrm m$、入力 Rayleigh 長 $z_R=\frac{1}{25}\,\mathrm m$、レンズの焦点距離 $f=\frac{1}{50}\,\mathrm m$、波長 $\lambda\simeq 1\times 10^{-6}\,\mathrm m$ である。

複素分母にその共役を掛ける。

$$
\begin{aligned}
D&=z_R^{2} + \left(f - x\right)^{2},\\
D/(\mathrm{m^2})&= \left(\frac{1}{25}\right)^{2} + \left(\left(\frac{1}{50}\right) - \left(\frac{1}{50}\right)\right)^{2}=\frac{1}{625}.
\end{aligned}
$$

レンズ直後での $q'$ の実部を自由伝播で打ち消す。

$$
\begin{aligned}
d&=- \frac{f \left(x \left(f - x\right) - z_R^{2}\right)}{D},\\
d/(\mathrm m)&= - \frac{\left(\frac{1}{50}\right) \cdot \left(\left(\frac{1}{50}\right) \cdot \left(\left(\frac{1}{50}\right) - \left(\frac{1}{50}\right)\right) - \left(\frac{1}{25}\right)^{2}\right)}{\left(\frac{1}{625}\right)}=\frac{1}{50}.
\end{aligned}
$$

虚部の符号を変えて幅の尺度を求める。

$$
\begin{aligned}
z'_R&=\frac{f^{2} z_R}{D},\\
z'_R/(\mathrm m)&= \frac{\left(\frac{1}{50}\right)^{2} \cdot \left(\frac{1}{25}\right)}{\left(\frac{1}{625}\right)}=\frac{1}{100}.
\end{aligned}
$$

Rayleigh 長からガウスの半径へ戻す。

$$
\begin{aligned}
w'_0&=\frac{\sqrt{z'_R \lambda}}{\sqrt{\pi}},\\
w'_0/(\mathrm m)&\simeq \frac{\sqrt{\left(\frac{1}{100}\right) \cdot \left(1\times 10^{-6}\right)}}{\sqrt{\pi}}\simeq 5.6419\times 10^{-5}.
\end{aligned}
$$

同じ波長での二つの腰の大きさを比べる。

$$
\begin{aligned}
w'_0/w_0&=\sqrt{\frac{z'_R}{z_R}},\\
w'_0/w_0&= \sqrt{\frac{\left(\frac{1}{100}\right)}{\left(\frac{1}{25}\right)}}=\frac{1}{2}.
\end{aligned}
$$

出力の腰は後焦点面になる。入力の波面曲率が有限であることを $q$ に含めた結果である。
:::

:::example[短い Rayleigh 長で発散する入射]
腰とレンズが近く、入力の発散が強い場合を調べる。

入力は入力の腰からレンズの距離 $x=\frac{1}{100}\,\mathrm m$、入力 Rayleigh 長 $z_R=\frac{1}{500}\,\mathrm m$、レンズの焦点距離 $f=\frac{3}{100}\,\mathrm m$、波長 $\lambda\simeq 1\times 10^{-6}\,\mathrm m$ である。

複素分母にその共役を掛ける。

$$
\begin{aligned}
D&=z_R^{2} + \left(f - x\right)^{2},\\
D/(\mathrm{m^2})&= \left(\frac{1}{500}\right)^{2} + \left(\left(\frac{3}{100}\right) - \left(\frac{1}{100}\right)\right)^{2}\simeq 4.04\times 10^{-4}.
\end{aligned}
$$

レンズ直後での $q'$ の実部を自由伝播で打ち消す。

$$
\begin{aligned}
d&=- \frac{f \left(x \left(f - x\right) - z_R^{2}\right)}{D},\\
d/(\mathrm m)&\simeq - \frac{\left(\frac{3}{100}\right) \cdot \left(\left(\frac{1}{100}\right) \cdot \left(\left(\frac{3}{100}\right) - \left(\frac{1}{100}\right)\right) - \left(\frac{1}{500}\right)^{2}\right)}{\left(4.04\times 10^{-4}\right)}\simeq -0.014554.
\end{aligned}
$$

虚部の符号を変えて幅の尺度を求める。

$$
\begin{aligned}
z'_R&=\frac{f^{2} z_R}{D},\\
z'_R/(\mathrm m)&\simeq \frac{\left(\frac{3}{100}\right)^{2} \cdot \left(\frac{1}{500}\right)}{\left(4.04\times 10^{-4}\right)}=\frac{9}{2020}.
\end{aligned}
$$

Rayleigh 長からガウスの半径へ戻す。

$$
\begin{aligned}
w'_0&=\frac{\sqrt{z'_R \lambda}}{\sqrt{\pi}},\\
w'_0/(\mathrm m)&\simeq \frac{\sqrt{\left(\frac{9}{2020}\right) \cdot \left(1\times 10^{-6}\right)}}{\sqrt{\pi}}\simeq 3.7659\times 10^{-5}.
\end{aligned}
$$

同じ波長での二つの腰の大きさを比べる。

$$
\begin{aligned}
w'_0/w_0&=\sqrt{\frac{z'_R}{z_R}},\\
w'_0/w_0&= \sqrt{\frac{\left(\frac{9}{2020}\right)}{\left(\frac{1}{500}\right)}}\simeq 1.4926.
\end{aligned}
$$

レンズ後でも発散が残り、出力の腰は逆伝播する方向に現れる。
:::

## 光共振器の安定条件

### 往復行列の固有値
二球面鏡の距離を $L$、凹面側が正の半径を $R_1,R_2$ とする。反射は展開座標で $M_j=\begin{pmatrix}1&0\\-2/R_j&1\end{pmatrix}$。往復 $M=M_1T(L)M_2T(L)$ の行列式は 1 である。トレース半分 $s=(A+D)/2$ に対して固有値は $s\pm\sqrt{s^2-1}$。$|s|<1$ なら異なる単位円上の根で、複素対角化によって全光線の繰り返しが有界となる。$|s|>1$ なら一方の固有方向に指数的に増大する。$s=\pm1$ では $M=\pm I$ なら有界だが、非自明な Jordan 部分があればそれが往復数に比例して成長する。

$g_j=1-L/R_j$ と置くと $s=2g_1g_2-1$。開区間の :term[共振器安定条件]{en="resonator stability condition" yomi="きょうしんきあんていじょうけん"} は
$$
0<g_1g_2<1 .
\label{eq-optics-cavity-stability}
$$
Gaussian の自己整合条件は $Cq^2+(D-A)q-B=0$。判別式は $(A+D)^2-4$ なので、安定区間では二共役根のうち負の虚部の根が一意に選ばれる。行列 $M=\pm I$ の境界例では固定点が一意に決まらない。

:::example[平面・凹面]
$R_1=\infty,\ R_2=2L$ なら $g_1g_2=1/2$。平面鏡の曲率と一致するようにそこへ腰を置くと、他端で $R(L)=(L^2+z_R^2)/L=2L$。従って $z_R=L$ となる。
:::
:::example[平面同士]
二枚が平面鏡なら $M=T(2L)$ である。$s=1$ でも $M\ne I$ なので、角度のある光線は毎往復で高さがずれ、有限幅のビームは回折で広がる。
:::

### 安定な二鏡のモード寸法
腰から鏡までの符号付き距離を $z_1,z_2$ とし、$z_1+z_2=L$。曲率の境界条件は $z_R^2=z_j(R_j-z_j)$ である。二式を引いて $z_1=L(R_2-L)/(R_1+R_2-2L)$、$z_2=L(R_1-L)/(R_1+R_2-2L)$ を得る。分母がゼロの対称共焦点のような場合は元の曲率式を使う。有限で非退化な鏡では $w_j^2=(\lambda/(\pi z_R))(z_R^2+z_j^2)$ を求め、鏡の物理的半径と比較する。安定条件は近軸の無限開口での条件であり、実際の開口が小さければ端での回折損失が増える。

:::example[二鏡の積で不安定を判定する]
$R_1=R_2=L/3$ なら $g_1=g_2=-2$ で積は 4、往復トレースの半分は 7 である。大きい方の固有値は $7+4\sqrt3$ なので、わずかな初期ずれでも該当固有方向は一往復で約 13.9 倍になる。
:::

### 平面・凹面共振器の腰から安定性を読む

平面鏡を腰に置く。もう一方の凹面鏡が距離 $L$ にあり曲率半径が $R$ なら、ビーム波面と鏡面の曲率が合う条件は $R=(L^2+z_R^2)/L$。従って $z_R^2=L(R-L)$。有限の腰を持つ解は $0<L<R$ に存在する。$L\to R$ では腰が狭くなるが凹面鏡上のビームは広がる。これは $w^2(L)=\lambda(z_R+L^2/z_R)/\pi$ より分かる。安定境界近くで小開口の鏡を使うと、幾何学的な安定性があってもクリッピング損失が増える。

:::example[中間の鏡間隔]
曲率半径の半分の距離でモードの寸法を求める。

入力は鏡間隔 $L=\frac{1}{20}\,\mathrm m$、凹面鏡の曲率半径 $R=\frac{1}{10}\,\mathrm m$、媒質中の波長 $\lambda\simeq 1\times 10^{-6}\,\mathrm m$、凹面鏡の有効半径 $a=\frac{1}{2000}\,\mathrm m$ である。

平面側の $g$ は一である。

$$
\begin{aligned}
g_1g_2&=- \frac{L}{R} + 1,\\
g_1g_2&= - \frac{\left(\frac{1}{20}\right)}{\left(\frac{1}{10}\right)} + 1=\frac{1}{2}.
\end{aligned}
$$

凹面で波面の曲率が鏡に合うようにする。

$$
\begin{aligned}
z_R&=\sqrt{L \left(- L + R\right)},\\
z_R/(\mathrm m)&= \sqrt{\left(\frac{1}{20}\right) \cdot \left(- \left(\frac{1}{20}\right) + \left(\frac{1}{10}\right)\right)}=\frac{1}{20}.
\end{aligned}
$$

平面波面での強度半径を求める。

$$
\begin{aligned}
w_0&=\frac{\sqrt{\lambda z_R}}{\sqrt{\pi}},\\
w_0/(\mathrm m)&\simeq \frac{\sqrt{\left(1\times 10^{-6}\right) \cdot \left(\frac{1}{20}\right)}}{\sqrt{\pi}}\simeq 1.2616\times 10^{-4}.
\end{aligned}
$$

鏡までの距離に応じた回折の広がりを含める。

$$
\begin{aligned}
w(L)&=w_0 \sqrt{\frac{L^{2}}{z_R^{2}} + 1},\\
w(L)/(\mathrm m)&\simeq \left(1.261566261\times 10^{-4}\right) \cdot \sqrt{\frac{\left(\frac{1}{20}\right)^{2}}{\left(\frac{1}{20}\right)^{2}} + 1}\simeq 1.7841\times 10^{-4}.
\end{aligned}
$$

場を未切断のガウスと近似し、外側の強度を積分する。

$$
\begin{aligned}
L_{\rm clip}&=e^{- \frac{2 a^{2}}{w(L)^{2}}},\\
L_{\rm clip}&\simeq e^{- \frac{2 \cdot \left(\frac{1}{2000}\right)^{2}}{\left(\left(1.784124116\times 10^{-4}\right)\right)^{2}}}\simeq 1.507\times 10^{-7}.
\end{aligned}
$$

凹面側の半径は腰より大きく、有限開口の評価にはその位置の幅が必要になる。
:::

:::example[半球形の境界に近づく]
平面鏡と凹面鏡の距離を増やす。

入力は鏡間隔 $L=\frac{99}{1000}\,\mathrm m$、凹面鏡の曲率半径 $R=\frac{1}{10}\,\mathrm m$、媒質中の波長 $\lambda\simeq 1\times 10^{-6}\,\mathrm m$、凹面鏡の有効半径 $a=\frac{1}{2000}\,\mathrm m$ である。

平面側の $g$ は一である。

$$
\begin{aligned}
g_1g_2&=- \frac{L}{R} + 1,\\
g_1g_2&= - \frac{\left(\frac{99}{1000}\right)}{\left(\frac{1}{10}\right)} + 1=\frac{1}{100}.
\end{aligned}
$$

凹面で波面の曲率が鏡に合うようにする。

$$
\begin{aligned}
z_R&=\sqrt{L \left(- L + R\right)},\\
z_R/(\mathrm m)&= \sqrt{\left(\frac{99}{1000}\right) \cdot \left(- \left(\frac{99}{1000}\right) + \left(\frac{1}{10}\right)\right)}\simeq 0.0099499.
\end{aligned}
$$

平面波面での強度半径を求める。

$$
\begin{aligned}
w_0&=\frac{\sqrt{\lambda z_R}}{\sqrt{\pi}},\\
w_0/(\mathrm m)&\simeq \frac{\sqrt{\left(1\times 10^{-6}\right) \cdot \left(0.009949874371\right)}}{\sqrt{\pi}}\simeq 5.6277\times 10^{-5}.
\end{aligned}
$$

鏡までの距離に応じた回折の広がりを含める。

$$
\begin{aligned}
w(L)&=w_0 \sqrt{\frac{L^{2}}{z_R^{2}} + 1},\\
w(L)/(\mathrm m)&\simeq \left(5.627737892\times 10^{-5}\right) \cdot \sqrt{\frac{\left(\frac{99}{1000}\right)^{2}}{\left(0.009949874371\right)^{2}} + 1}\simeq 5.6277\times 10^{-4}.
\end{aligned}
$$

場を未切断のガウスと近似し、外側の強度を積分する。

$$
\begin{aligned}
L_{\rm clip}&=e^{- \frac{2 a^{2}}{w(L)^{2}}},\\
L_{\rm clip}&\simeq e^{- \frac{2 \cdot \left(\frac{1}{2000}\right)^{2}}{\left(\left(5.627737892\times 10^{-4}\right)\right)^{2}}}\simeq 0.20624.
\end{aligned}
$$

腰が縮んでも、凹面鏡での切り落としが増える。
:::

:::example[鏡の半径を広げる]
同じ安定性の積を保ったまま凹面鏡の有効径だけを増す。

入力は鏡間隔 $L=\frac{99}{1000}\,\mathrm m$、凹面鏡の曲率半径 $R=\frac{1}{10}\,\mathrm m$、媒質中の波長 $\lambda\simeq 1\times 10^{-6}\,\mathrm m$、凹面鏡の有効半径 $a=\frac{1}{1000}\,\mathrm m$ である。

平面側の $g$ は一である。

$$
\begin{aligned}
g_1g_2&=- \frac{L}{R} + 1,\\
g_1g_2&= - \frac{\left(\frac{99}{1000}\right)}{\left(\frac{1}{10}\right)} + 1=\frac{1}{100}.
\end{aligned}
$$

凹面で波面の曲率が鏡に合うようにする。

$$
\begin{aligned}
z_R&=\sqrt{L \left(- L + R\right)},\\
z_R/(\mathrm m)&= \sqrt{\left(\frac{99}{1000}\right) \cdot \left(- \left(\frac{99}{1000}\right) + \left(\frac{1}{10}\right)\right)}\simeq 0.0099499.
\end{aligned}
$$

平面波面での強度半径を求める。

$$
\begin{aligned}
w_0&=\frac{\sqrt{\lambda z_R}}{\sqrt{\pi}},\\
w_0/(\mathrm m)&\simeq \frac{\sqrt{\left(1\times 10^{-6}\right) \cdot \left(0.009949874371\right)}}{\sqrt{\pi}}\simeq 5.6277\times 10^{-5}.
\end{aligned}
$$

鏡までの距離に応じた回折の広がりを含める。

$$
\begin{aligned}
w(L)&=w_0 \sqrt{\frac{L^{2}}{z_R^{2}} + 1},\\
w(L)/(\mathrm m)&\simeq \left(5.627737892\times 10^{-5}\right) \cdot \sqrt{\frac{\left(\frac{99}{1000}\right)^{2}}{\left(0.009949874371\right)^{2}} + 1}\simeq 5.6277\times 10^{-4}.
\end{aligned}
$$

場を未切断のガウスと近似し、外側の強度を積分する。

$$
\begin{aligned}
L_{\rm clip}&=e^{- \frac{2 a^{2}}{w(L)^{2}}},\\
L_{\rm clip}&\simeq e^{- \frac{2 \cdot \left(\frac{1}{1000}\right)^{2}}{\left(\left(5.627737892\times 10^{-4}\right)\right)^{2}}}\simeq 0.0018093.
\end{aligned}
$$

光線行列の安定条件は同じだが、切り落としは指数的に減る。
:::

## 共振器の縦モードと横モード

### 縦の整数と横次数
均一で非分散の共振器で、腰が二鏡間にある場合をまず考える。腰から両端までの距離を $z_1,z_2\ge0$ とすると $z_1+z_2=L$。片道 Gouy 位相の大きさを $\Psi=\arctan(z_1/z_R)+\arctan(z_2/z_R)$ とすると、高次モードの往復位相は $2kL-2(N+1)\Psi$ に二鏡の反射位相の和を加えたものになる。鏡の反射位相の和が $2\pi$ の整数倍であり、全往復位相も $2\pi$ の整数倍となるとき、共鳴周波数は
$$
\nu_{j,N}=\frac{c}{2nL}\left[j+\frac{N+1}{\pi}\Psi\right].
\label{eq-optics-resonator-modes}
$$
$j$ が縦モードを、$N$ が横次数を指定する。例えば対称共振器で $z_1=z_2=L/2$、$R_1=R_2=R$ なら $z_R^2=(L/2)(R-L/2)$。$\cos[2\arctan(L/(2z_R))]=(z_R^2-(L/2)^2)/(z_R^2+(L/2)^2)=1-L/R=g$、従って $\Psi=\arccos g$。$g<0$ の枝では $\Psi>\pi/2$ となる。次数差を求める際に常に $\arccos\sqrt{g_1g_2}$ を使うと、この枝を取り違える。

対称性を仮定しない場合でも枝を決められる。$g_1=(z_R^2-z_1z_2)/(z_R^2+z_1^2)$、$g_2=(z_R^2-z_1z_2)/(z_R^2+z_2^2)$ は鏡の曲率条件を代入したものである。一方、二つの逆正接の余弦の加法公式から
$$
\cos\Psi=\frac{z_R^2-z_1z_2}{\sqrt{(z_R^2+z_1^2)(z_R^2+z_2^2)}} .
$$
従って $0<g_1g_2<1$ では $\Psi=\arccos(\operatorname{sgn}(g_1)\sqrt{g_1g_2})$。符号付きの $z_j$ を用いても、$z_1+z_2=L>0$ なので $0<\Psi<\pi$ となる。

:::example[平面・凹面の横間隔]
前節の $z_R=L$ では $\Psi=\arctan1=\pi/4$。同じ縦番号で横次数を 1 上げると周波数は縦モード間隔の $1/4$ だけ上がる。
:::
:::example[共焦点の縮退]
対称で $R=L$ なら $\Psi=\pi/2$。$N$ を 2 増やして $j$ を 1 減らすと同じ周波数になる。同次数の方位の異なるモードも理想回転対称系では縮退する。
:::

式 :eqref[eq-optics-resonator-modes] の $\Psi$ は、腰が共振器外にある場合にも、距離 $z_j$ を曲率条件で使った符号付きの値とすれば適用できる。あるいは光軸上の固定座標で二鏡間の $\arctan[(z-z_{\rm waist})/z_R]$ の増分をそのまま取る。鏡の位相和 $\phi_m$ が $2\pi$ の整数倍でないときは、同式の角括弧の中へ $-\phi_m/(2\pi)$ を残す。整数 $j$ の値だけを動かして任意の定数位相を除くことはできない。

### 二方向で異なる横モード間隔
光学系が二つの直交軸に分離でき、それぞれの $q$ と片道 Gouy 増分を $\Psi_x,\Psi_y$ とする。高次モードを導いたときの変数分離を一軸ずつ行うと、一次元の基本解の前因子は $(q_0/q)^{1/2}$ になる。従って各軸の基本モードが受ける位相は、その軸の Gouy 増分の半分である。Hermite 次数 $(m,n)$ には片道で $-(m+1/2)\Psi_x-(n+1/2)\Psi_y$ が加わる。従って整数条件の横の項は $[(m+1/2)\Psi_x+(n+1/2)\Psi_y]/\pi$ になる。非点共振器では $m+n$ が同じでもこの和が違い、縮退が解ける。二方向とも安定でなければ横断面全体の正規化モードにならない。

:::example[非点の縦周波数差を測る]
$\Psi_x=\pi/3,\Psi_y=\pi/4$ の系で縦間隔が $600\,\mathrm{MHz}$ なら、同じ縦番号の $(m,n)=(1,0)$ と $(0,1)$ の差は $600(1/3-1/4)\,\mathrm{MHz}=50\,\mathrm{MHz}$ である。基本モードの位置を動かす二方向の共通項はこの差から消える。
:::

### Gouy 位相で横モードの重なりを判定する

以下の例では真空の平面・凹面共振器を用いる。この配置で安定な場合、片道の Gouy 角は $\arctan(L/z_R)=\arccos\sqrt{1-L/R}$ である。鏡の反射位相和を $0$ として一往復の位相条件を $2kL-2(N+1)\psi=2\pi j$ と書く。$N$ は横次数、$j$ は縦次数である。よって縦次数を同じにして横次数を $\Delta N$ 増す周波数差は $\Delta\nu_{\rm FSR}\Delta N\psi/\pi$。この値が整数個の FSR に一致するなら、縦次数の異なるモードが縮退する。実際の共振器で線幅より差が小さいだけなら、周波数の一致とは区別して分離できないという判断をする。

:::example[隣の横モード]
鏡間隔を半径の半分とし、基本モードと一次を比較する。

入力は鏡間隔 $L=\frac{1}{20}\,\mathrm m$、凹面の曲率半径 $R=\frac{1}{10}\,\mathrm m$、比較する二モードの横次数差 $\Delta N=1$ である。

平面側を一として積を求める。

$$
\begin{aligned}
g&=- \frac{L}{R} + 1,\\
g&= - \frac{\left(\frac{1}{20}\right)}{\left(\frac{1}{10}\right)} + 1=\frac{1}{2}.
\end{aligned}
$$

半球形に近いほど角は $\pi/2$ に近づく。

$$
\begin{aligned}
\psi&=\operatorname{acos}{\left(\sqrt{g} \right)},\\
\psi/(\mathrm{rad})&= \operatorname{acos}{\left(\sqrt{\left(\frac{1}{2}\right)} \right)}\simeq 0.7854.
\end{aligned}
$$

空気近似で一往復時間の逆数を求める。

$$
\begin{aligned}
\Delta\nu_{\rm FSR}&=\frac{c}{2L},\\
\Delta\nu_{\rm FSR}/(\mathrm{Hz})&= \frac{149896229}{\left(\frac{1}{20}\right)}=2997924580.
\end{aligned}
$$

横次数だけを変えたときの位相差を周期で数える。

$$
\begin{aligned}
u&=\frac{\Delta N \psi}{\pi},\\
u&\simeq \frac{\left(1\right) \cdot \left(0.7853981634\right)}{\pi}=\frac{1}{4}.
\end{aligned}
$$

実際の周波数軸へ戻す。

$$
\begin{aligned}
\Delta\nu&=\Delta\nu_{\rm FSR} u,\\
\Delta\nu/(\mathrm{Hz})&= \left(2997924580\right) \cdot \left(\frac{1}{4}\right)=749481145.
\end{aligned}
$$

周波数差は一 FSR の四分の一である。
:::

:::example[異なる縦次数との重なり]
同じ共振器で横次数を四つ離す。

入力は鏡間隔 $L=\frac{1}{20}\,\mathrm m$、凹面の曲率半径 $R=\frac{1}{10}\,\mathrm m$、比較する二モードの横次数差 $\Delta N=4$ である。

平面側を一として積を求める。

$$
\begin{aligned}
g&=- \frac{L}{R} + 1,\\
g&= - \frac{\left(\frac{1}{20}\right)}{\left(\frac{1}{10}\right)} + 1=\frac{1}{2}.
\end{aligned}
$$

半球形に近いほど角は $\pi/2$ に近づく。

$$
\begin{aligned}
\psi&=\operatorname{acos}{\left(\sqrt{g} \right)},\\
\psi/(\mathrm{rad})&= \operatorname{acos}{\left(\sqrt{\left(\frac{1}{2}\right)} \right)}\simeq 0.7854.
\end{aligned}
$$

空気近似で一往復時間の逆数を求める。

$$
\begin{aligned}
\Delta\nu_{\rm FSR}&=\frac{c}{2L},\\
\Delta\nu_{\rm FSR}/(\mathrm{Hz})&= \frac{149896229}{\left(\frac{1}{20}\right)}=2997924580.
\end{aligned}
$$

横次数だけを変えたときの位相差を周期で数える。

$$
\begin{aligned}
u&=\frac{\Delta N \psi}{\pi},\\
u&\simeq \frac{\left(4\right) \cdot \left(0.7853981634\right)}{\pi}=1.
\end{aligned}
$$

実際の周波数軸へ戻す。

$$
\begin{aligned}
\Delta\nu&=\Delta\nu_{\rm FSR} u,\\
\Delta\nu/(\mathrm{Hz})&= \left(2997924580\right) \cdot \left(1\right)=2997924580.
\end{aligned}
$$

同じ縦次数では一 FSR 離れるので、隣の縦次数を選べば一致する。
:::

:::example[重なりを外す]
横次数差を保ち、鏡間隔を変えて位相の比を動かす。

入力は鏡間隔 $L=\frac{1}{25}\,\mathrm m$、凹面の曲率半径 $R=\frac{1}{10}\,\mathrm m$、比較する二モードの横次数差 $\Delta N=4$ である。

平面側を一として積を求める。

$$
\begin{aligned}
g&=- \frac{L}{R} + 1,\\
g&= - \frac{\left(\frac{1}{25}\right)}{\left(\frac{1}{10}\right)} + 1=\frac{3}{5}.
\end{aligned}
$$

半球形に近いほど角は $\pi/2$ に近づく。

$$
\begin{aligned}
\psi&=\operatorname{acos}{\left(\sqrt{g} \right)},\\
\psi/(\mathrm{rad})&= \operatorname{acos}{\left(\sqrt{\left(\frac{3}{5}\right)} \right)}\simeq 0.68472.
\end{aligned}
$$

空気近似で一往復時間の逆数を求める。

$$
\begin{aligned}
\Delta\nu_{\rm FSR}&=\frac{c}{2L},\\
\Delta\nu_{\rm FSR}/(\mathrm{Hz})&= \frac{149896229}{\left(\frac{1}{25}\right)}=3747405725.
\end{aligned}
$$

横次数だけを変えたときの位相差を周期で数える。

$$
\begin{aligned}
u&=\frac{\Delta N \psi}{\pi},\\
u&\simeq \frac{\left(4\right) \cdot \left(0.684719203\right)}{\pi}\simeq 0.87181.
\end{aligned}
$$

実際の周波数軸へ戻す。

$$
\begin{aligned}
\Delta\nu&=\Delta\nu_{\rm FSR} u,\\
\Delta\nu/(\mathrm{Hz})&\simeq \left(3747405725\right) \cdot \left(0.8718115663\right)\simeq 3.267\times 10^{9}.
\end{aligned}
$$

FSR に対する比が整数から離れ、前例の重なりが解ける。
:::

## フィネスと光子寿命

### 一往復のエネルギー減衰
無入力で一往復の強度生存率を $0<\mathcal R<1$、往復時間を $t_{\rm rt}$ とする。$m$ 周後のエネルギーは $\mathcal R^m$ 倍なので指数補間して
$$
\tau_{\mathrm{ph}}=-\frac{t_{\rm rt}}{\log\mathcal R}.
\label{eq-optics-photon-lifetime}
$$
$\tau_{\mathrm{ph}}$ を :term[光子寿命]{en="photon lifetime" yomi="こうしじゅみょう"} という。振幅の生存因子は $\sqrt{\mathcal R}$。幾何級数の線形応答で式 :eqref[eq-optics-airy-transmission] の分母の $R$ を $\sqrt{\mathcal R}$ に置き換えれば幅が分かる。高フィネスでは $1-\sqrt{\mathcal R}\simeq t_{\rm rt}/(2\tau_{\rm ph})$ で $\Delta\nu\simeq1/(2\pi\tau_{\rm ph})$。自由スペクトル間隔は $1/t_{\rm rt}$ なので $\mathcal F\simeq2\pi\tau_{\rm ph}/t_{\rm rt}$、品質係数は $Q=\nu_0/\Delta\nu\simeq 2\pi\nu_0\tau_{\rm ph}$ である。エネルギー寿命を振幅寿命と取り違えると幅が二倍違う。

:::example[リングダウン]
$t_{\rm rt}=2\,\mathrm{ns}$、$\mathcal R=0.99$ なら $\tau_{\rm ph}\simeq199\,\mathrm{ns}$、高フィネスで線幅は約 $0.80\,\mathrm{MHz}$。
:::
:::example[損失の積]
鏡二枚の反射率がともに 0.99 で内部損失がないなら、往復生存率は $0.99^2$ である。同じ往復時間で $\mathcal R=0.99$ とした計算に比べ、寿命は半分になる。
:::

### 有限開口の損失を見積もる
半径 $a$ の鏡に強度半径 $w$ の軸対称基本モードを入れる。鏡の外にあるパワー比は $\int_a^\infty r e^{-2r^2/w^2}\dd r/\int_0^\infty r e^{-2r^2/w^2}\dd r=e^{-2a^2/w^2}$。例えば $a=2w$ なら約 $3.35\times10^{-4}$ である。切られた場は厳密には元のモードと同じ形ではないため、この式は小さなクリッピングの損失見積りとして使う。

### 入力を止める実験
定常入力を切った後の出力を、残存エネルギーに比例する一枚の出力鏡で検出する。検出器の応答時間が $\tau_{\rm ph}$ より短く、背景を引ければ、パワーを時刻に対して対数表示した傾きは $-1/\tau_{\rm ph}$。一方、電場振幅を参照光と混ぜて振幅として測るなら傾きは $-1/(2\tau_{\rm ph})$ になる。周波数応答の線幅と比較する際はどちらの信号を記録したかをそろえる。

:::example[低いフィネスで近似を比較する]
エネルギーの一往復生存率が $1/4$ なら振幅の生存率は $1/2$ である。式 :eqref[eq-optics-airy-transmission] より全半値幅は位相で $4\arcsin(1/\sqrt8)$、フィネスは約 $4.35$。寿命からの高フィネス近似 $2\pi/\ln4$ は約 $4.53$ となる。
:::

### リングダウンの時間定数と周波数線幅

一往復時間を $t_r$、パワー保持率を $\mathcal R$ とすれば、$m$ 往復後のエネルギーは $U_0\mathcal R^m=U_0\exp(-mt_r/\tau)$、$\tau=-t_r/\ln\mathcal R$。振幅はその平方根で減衰する。共振の規格化透過を $[1+4\rho\sin^2(\delta/2)/(1-\rho)^2]^{-1}$、$\rho=\sqrt{\mathcal R}$、$\delta=2\pi(\nu-\nu_0)t_r$ とすると、左右の半値点から全半値幅は $2\arcsin[(1-\rho)/(2\sqrt{\rho})]/(\pi t_r)$。引数が一以下、すなわち $\rho\ge(\sqrt2-1)^2$ の場合に限り半値点が存在する。これより損失が大きい場合、隣の共鳴との中間でも透過はピークの半分を上回る。高保持率では $-\ln\mathcal R\simeq2(1-\rho)$ なので、この全幅は $1/(2\pi\tau)$ に近づく。

:::example[保持率が高い共振器]
一往復あたりの損失が小さい場合を調べる。

入力は往復時間 $t_r\simeq 1\times 10^{-9}\,\mathrm s$、往復パワー保持率 $\mathcal R=\frac{99}{100}$ である。

パワーと振幅の減衰率を区別する。

$$
\begin{aligned}
\rho&=\sqrt{\mathcal R},\\
\rho&= \sqrt{\left(\frac{99}{100}\right)}\simeq 0.99499.
\end{aligned}
$$

エネルギーの指数減衰を一往復で合わせる。

$$
\begin{aligned}
\tau&=- \frac{t_r}{\log{\left(\mathcal R \right)}},\\
\tau/(\mathrm s)&\simeq - \frac{\left(1\times 10^{-9}\right)}{\log{\left(\left(\frac{99}{100}\right) \right)}}\simeq 9.9499\times 10^{-8}.
\end{aligned}
$$

規格化した共振透過が二分の一になる条件を出す。

$$
\begin{aligned}
a&=\frac{1 - \rho}{2 \sqrt{\rho}},\\
a&\simeq \frac{1 - \left(0.9949874371\right)}{2 \cdot \sqrt{\left(0.9949874371\right)}}\simeq 0.0025126.
\end{aligned}
$$

離調の左右二点の差に換算する。

$$
\begin{aligned}
\Delta\nu&=\frac{2 \operatorname{asin}{\left(a \right)}}{\pi t_r},\\
\Delta\nu/(\mathrm{Hz})&\simeq \frac{2 \cdot \operatorname{asin}{\left(\left(0.002512586607\right) \right)}}{\pi \cdot \left(1\times 10^{-9}\right)}\simeq 1.5996\times 10^{6}.
\end{aligned}
$$

十分狭い共振の指数寿命から見積もる。

$$
\begin{aligned}
\Delta\nu_{\rm app}&=\frac{1}{2 \pi \tau},\\
\Delta\nu_{\rm app}/(\mathrm{Hz})&\simeq \frac{1}{2 \cdot \pi \cdot \left(9.949916247\times 10^{-8}\right)}\simeq 1.5996\times 10^{6}.
\end{aligned}
$$

二つの幅はよく一致する。
:::

:::example[往復損失が無視できない共振器]
エネルギー保持率を下げ、同じ往復時間で比較する。

入力は往復時間 $t_r\simeq 1\times 10^{-9}\,\mathrm s$、往復パワー保持率 $\mathcal R=\frac{1}{2}$ である。

パワーと振幅の減衰率を区別する。

$$
\begin{aligned}
\rho&=\sqrt{\mathcal R},\\
\rho&= \sqrt{\left(\frac{1}{2}\right)}\simeq 0.70711.
\end{aligned}
$$

エネルギーの指数減衰を一往復で合わせる。

$$
\begin{aligned}
\tau&=- \frac{t_r}{\log{\left(\mathcal R \right)}},\\
\tau/(\mathrm s)&\simeq - \frac{\left(1\times 10^{-9}\right)}{\log{\left(\left(\frac{1}{2}\right) \right)}}\simeq 1.4427\times 10^{-9}.
\end{aligned}
$$

規格化した共振透過が二分の一になる条件を出す。

$$
\begin{aligned}
a&=\frac{1 - \rho}{2 \sqrt{\rho}},\\
a&\simeq \frac{1 - \left(0.7071067812\right)}{2 \cdot \sqrt{\left(0.7071067812\right)}}\simeq 0.17416.
\end{aligned}
$$

離調の左右二点の差に換算する。

$$
\begin{aligned}
\Delta\nu&=\frac{2 \operatorname{asin}{\left(a \right)}}{\pi t_r},\\
\Delta\nu/(\mathrm{Hz})&\simeq \frac{2 \cdot \operatorname{asin}{\left(\left(0.1741553499\right) \right)}}{\pi \cdot \left(1\times 10^{-9}\right)}\simeq 1.1144\times 10^{8}.
\end{aligned}
$$

十分狭い共振の指数寿命から見積もる。

$$
\begin{aligned}
\Delta\nu_{\rm app}&=\frac{1}{2 \pi \tau},\\
\Delta\nu_{\rm app}/(\mathrm{Hz})&\simeq \frac{1}{2 \cdot \pi \cdot \left(1.442695041\times 10^{-9}\right)}\simeq 1.1032\times 10^{8}.
\end{aligned}
$$

寿命が短くなり、線幅が増える。正弦による線形からの補正も現れる。
:::

:::example[広い共振の半値を測る]
保持率をさらに下げる。

入力は往復時間 $t_r\simeq 1\times 10^{-9}\,\mathrm s$、往復パワー保持率 $\mathcal R=\frac{1}{10}$ である。

パワーと振幅の減衰率を区別する。

$$
\begin{aligned}
\rho&=\sqrt{\mathcal R},\\
\rho&= \sqrt{\left(\frac{1}{10}\right)}\simeq 0.31623.
\end{aligned}
$$

エネルギーの指数減衰を一往復で合わせる。

$$
\begin{aligned}
\tau&=- \frac{t_r}{\log{\left(\mathcal R \right)}},\\
\tau/(\mathrm s)&\simeq - \frac{\left(1\times 10^{-9}\right)}{\log{\left(\left(\frac{1}{10}\right) \right)}}\simeq 4.3429\times 10^{-10}.
\end{aligned}
$$

規格化した共振透過が二分の一になる条件を出す。

$$
\begin{aligned}
a&=\frac{1 - \rho}{2 \sqrt{\rho}},\\
a&\simeq \frac{1 - \left(0.316227766\right)}{2 \cdot \sqrt{\left(0.316227766\right)}}\simeq 0.60797.
\end{aligned}
$$

離調の左右二点の差に換算する。

$$
\begin{aligned}
\Delta\nu&=\frac{2 \operatorname{asin}{\left(a \right)}}{\pi t_r},\\
\Delta\nu/(\mathrm{Hz})&\simeq \frac{2 \cdot \operatorname{asin}{\left(\left(0.6079690424\right) \right)}}{\pi \cdot \left(1\times 10^{-9}\right)}\simeq 4.1603\times 10^{8}.
\end{aligned}
$$

十分狭い共振の指数寿命から見積もる。

$$
\begin{aligned}
\Delta\nu_{\rm app}&=\frac{1}{2 \pi \tau},\\
\Delta\nu_{\rm app}/(\mathrm{Hz})&\simeq \frac{1}{2 \cdot \pi \cdot \left(4.342944819\times 10^{-10}\right)}\simeq 3.6647\times 10^{8}.
\end{aligned}
$$

半値点は存在するが、単純な寿命からの近似との差が大きい。
:::

## 演習問題

::::exercise{level=basic}
$\psi=\arctan(z/z_R)$ とする。$z=z_R$ における $w,R,\psi$ と、基本モードが平面波に対して受ける Gouy 位相を求めよ。

:::solution
幅の式に $z/z_R=1$ を代入して $w=\sqrt2w_0$。曲率半径は $R=2z_R$、$\psi=\arctan1=\pi/4$ なので、基本モードが受ける Gouy 位相は $-\pi/4$ である。
:::
::::

::::exercise{level=basic}
$w_0$ が二倍になると $z_R$ と発散半角はどうなるか。

:::solution
$z_R=\pi w_0^2/\lambda$ より四倍、$\theta=\lambda/(\pi w_0)$ より半分になる。
:::
::::

::::exercise{level=basic}
$L=0.30\,\mathrm m$ の真空共振器で縦モードの間隔を求めよ。

:::solution
$c/(2L)\simeq3.0\times10^8/(0.60)=0.50\,\mathrm{GHz}$ である。
:::
::::

::::exercise{level=basic}
$(p,\ell)=(2,-3)$ の Laguerre–Gauss の横次数を求めよ。

:::solution
式で $N=2p+\abs{\ell}=4+3=7$ となり、Gouy 位相には $N+1=8$ が掛かる。
:::
::::

::::exercise{level=advanced}
$q=a-ib,\ b>0$ を実係数で行列式 1 の ABCD 写像へ通す。虚部の符号が保たれることを示せ。

:::solution
出力に分母の共役を掛ける。分子の虚部は $(AD-BC)\Im q$ である。正の分母 $\abs{Cq+D}^2$ で割るから符号は保たれる。
:::
::::

::::exercise{level=advanced}
対称共振器の $R=2L$ のとき安定性と横モード間隔の比を求めよ。

:::solution
$g=1-L/R=1/2$ で積 $1/4$ は安定域内。$\Psi=\arccos(1/2)=\pi/3$ なので横間隔は縦間隔の $1/3$ である。
:::
::::

::::exercise{level=advanced}
$t_{\rm rt}=4\,\mathrm{ns}$、往復生存率 $e^{-0.02}$ の寿命を求めよ。

:::solution
:eqref[eq-optics-photon-lifetime] より $\tau_{\rm ph}=4/0.02\,\mathrm{ns}=200\,\mathrm{ns}$ である。
:::
::::

::::exercise{level=advanced}
$H_2$ を Rodrigues の式から求め、その根の $\xi$ を求めよ。

:::solution
二階微分は $(4\xi^2-2)e^{-\xi^2}$ なので $H_2=4\xi^2-2$。節は $\xi=\pm1/\sqrt2$ である。
:::
::::

::::exercise{level=advanced}
一往復にパワーの割合 $e^{-0.04}$ が残る共振器で、往復時間が $2\,\mathrm{ns}$ とする。パワーと振幅の指数寿命、および高フィネス近似で線幅を求めよ。

:::solution
式 :eqref[eq-optics-photon-lifetime] から $\tau_P=2\,\mathrm{ns}/0.04=50\,\mathrm{ns}$。振幅はパワーの平方根に比例するので $e^{-t/(2\tau_P)}$ で減衰し、寿命は $100\,\mathrm{ns}$。線幅はパワー寿命から $1/(2\pi\tau_P)=(10/\pi)\,\mathrm{MHz}\simeq3.18\,\mathrm{MHz}$ となる。
:::
::::

## 確認テスト

::::quiz
ガウスビームの腰で成り立つものをすべて選べ。

- [x] 幅は最小である。
- [x] 波面の曲率はゼロである。
- [ ] Gouy 位相のずれの大きさは $\pi/2$ である。

:::explanation
$q=-iz_R$ では式 :eqref[eq-optics-gaussian] の $1/q$ は純虚数で、幅の式は最小値 $w_0$ を取る。腰を基準とした $\psi$ はゼロである。
:::
::::

::::quiz{answer="9" tol="1e-6"}
$w_0$ が三倍なら同じ波長での Rayleigh 長は何倍か。

:::explanation
$z_R$ は $w_0^2$ に比例するので $3^2=9$ 倍になる。
:::
::::

::::quiz
安定区間で同じ往復行列を何回も作用させると光線はどうなるか。

- [x] 有界にとどまる。
- [ ] どの光線も指数的に減衰する。
- [ ] どの光線も無限遠へ進む。

:::explanation
二つの固有値が異なり、絶対値が 1 である。対角化した座標の成分は何乗しても有界で、固定の基底へ戻しても有界である。
:::
::::

::::quiz{answer="200" tol="1e-6"}
往復時間 $5\,\mathrm{ns}$ の共振器の自由スペクトル間隔を MHz 単位で答えよ。

:::explanation
間隔は往復時間の逆数なので $1/(5\times10^{-9}\,\mathrm s)=2\times10^8\,\mathrm{Hz}=200\,\mathrm{MHz}$。
:::
::::

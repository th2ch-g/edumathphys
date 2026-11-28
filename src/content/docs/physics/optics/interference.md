---
title: 干渉と干渉計
description: 光路差と境界での位相から干渉縞を計算し、二光束・多光束干渉計と誘電体膜を設計する。
notation:
  - symbol: '\Gamma_{12}'
    meaning: 二光束の相互相関
  - symbol: '\mathcal V'
    meaning: 縞の可視度
  - symbol: '\delta'
    meaning: 二光束または共振器一往復の位相差
  - symbol: 'R,T'
    meaning: パワー反射率と透過率
  - symbol: 'Y'
    meaning: 接線場に対する光学アドミタンス
---

## この章の目標

- 強度と相関から干渉縞の明暗を求める。
- 波面分割と振幅分割の光路を計算する。
- 干渉計で変位、屈折率および周波数を測る原理を説明する。
- 多光束の和と膜の行列から透過スペクトルを求める。

## 2 光束干渉と可干渉性の条件

### 検出器が平均する量
同じ媒体で、ほぼ同方向かつ同偏光の二波を一つの面上で重ねる。電場の複素包絡線を $u_1,u_2$ とし、共通の強度換算係数を振幅に吸収する。片方へ制御位相 $\phi$ を加え、光周期より長い平均を取ると
$$
I(\phi)=I_1+I_2+2\Re[e^{i\phi}\Gamma_{12}],
\quad I_j=\langle\abs{u_j}^2\rangle,\quad \Gamma_{12}=\langle u_1^*u_2\rangle .
\label{eq-optics-two-beam}
$$
これは $\langle|u_1+e^{i\phi}u_2|^2\rangle$ を展開した結果である。二偏光を使うと積は横ベクトルの内積に変わり、完全に直交する二成分では交差項が消える。違う周波数なら交差項は差周波数で振動するので、平均時間がその周期より十分長いと消える。短い平均ならうなりを検出できる。

縞の :term[可視度]{en="visibility" yomi="かしど"} $\mathcal V=(I_{\max}-I_{\min})/(I_{\max}+I_{\min})$ は $2|\Gamma_{12}|/(I_1+I_2)$ である。Cauchy–Schwarz の不等式から $|\Gamma_{12}|\le\sqrt{I_1I_2}$、さらに $(\sqrt{I_1}-\sqrt{I_2})^2\ge0$ なので $\mathcal V\le1$。視認できる強度があるとき、$\mathcal V=1$ となる必要十分条件は $u_2=e^{i\phi_0}u_1$ が平均の範囲でほとんど至る所成り立つことである。実際、Cauchy–Schwarz の等号では二波が一定の複素係数で比例し、相加相乗平均の等号から両波の平均強度が等しくなり、その係数の絶対値は 1 になる。遅延がある場合の相関は [コヒーレンス](/physics/optics/coherence/) で扱う。

:::example[不釣合い]
固定位相で $I_1=9,I_2=1$ とすると $I_{\max}=16,I_{\min}=4$、従って $\mathcal V=3/5$。位相が安定でも暗線はゼロにならない。
:::
:::example[位相の平均]
等強度で位相誤差が $+\pi/2$ と $-\pi/2$ へ等確率で跳ぶ場合、平均位相因子は $(i-i)/2=0$。長い露光では縞が消える。
:::

### 有限の露光でうなりを平均する
二波の周波数差を $\Delta\nu$、露光の中央を $t_0$、長さを $T$ とする。交差項の時間因子の平均は $\frac1T\int_{t_0-T/2}^{t_0+T/2}e^{-2\pi i\Delta\nu t}\dd t=e^{-2\pi i\Delta\nu t_0}\operatorname{sinc}(\pi\Delta\nu T)$。$\operatorname{sinc}x=\sin x/x$ とし原点は 1 である。異なる周波数であっても、短い露光で位相を既知の時刻にそろえれば交差項を測れる。$\Delta\nu T=1$ では丸一周期を平均するのでゼロになる。

:::example[位相を走査する検出器の露光]
同じ偏光の二波の差が $100\,\mathrm{kHz}$ のとき、露光を $10\,\mathrm{\mu s}$ にすると一周期を平均して交差項が消える。$5\,\mathrm{\mu s}$ なら $\Delta\nu T=1/2$ で交差項の振幅の大きさは $2/\pi$ 倍残る。時間窓の中心を同じ位相にそろえることが条件である。
:::

### 直交する偏光が干渉を隠す場合

複素振幅を $a_j\vb u_j$ とし $\vb u_j^\dagger\vb u_j=1$ とする。両方に等しい感度を持つ検出器の二乗平均には、交差項として $2\Re[\langle a_1^*a_2\rangle\vb u_1^\dagger\vb u_2]$ が現れる。従って、位相の関係を固定しても偏光が直交するとこの項がゼロになる。検光子を置いた後には $\vb u_j$ を一つの成分へ射影し、二光束の強度と交差項をその成分で取り直す。例えば $x,y$ の直交入射でも $45^\circ$ の検光子の後では元の強度を両方半分にした干渉になる。位相乱れと偏光の非一致はどちらも縞を薄くするが、検光子を回転させたときの効果が異なる。

:::example[同偏光だが光量が不均衡]
時間的には位相が全く乱れない二光束を重ねる。

入力は第一光束の強度 $I_1=9\,\mathrm{W/m^2}$、第二光束の強度 $I_2=1\,\mathrm{W/m^2}$、二つの直線偏光の差 $\theta=0\,\mathrm{rad}$、時間的な複素相関の絶対値 $\lvert\gamma\rvert=1$ である。

二光束を別々に点灯したときの強度を足す。

$$
\begin{aligned}
S&=I_1 + I_2,\\
S/(\mathrm{W/m^2})&= \left(9\right) + \left(1\right)=10.
\end{aligned}
$$

偏光の内積と時間相関を交差項に掛ける。

$$
\begin{aligned}
I_{\rm mod}&=2 \lvert\gamma\rvert \sqrt{I_1 I_2} \left|{\cos{\left(\theta \right)}}\right|,\\
I_{\rm mod}/(\mathrm{W/m^2})&= 2 \cdot \left(1\right) \cdot \sqrt{\left(9\right) \cdot \left(1\right)} \cdot \left|{\cos{\left(\left(0\right) \right)}}\right|=6.
\end{aligned}
$$

交差項が正に最大となる位相を選ぶ。

$$
\begin{aligned}
I_{\max}&=I_{\rm mod} + S,\\
I_{\max}/(\mathrm{W/m^2})&= \left(6\right) + \left(10\right)=16.
\end{aligned}
$$

位相を半周期変えて交差項を引く。

$$
\begin{aligned}
I_{\min}&=- I_{\rm mod} + S,\\
I_{\min}/(\mathrm{W/m^2})&= - \left(6\right) + \left(10\right)=4.
\end{aligned}
$$

明暗の差を和で割ると変動幅と背景の比になる。

$$
\begin{aligned}
\mathcal V&=\frac{I_{\rm mod}}{S},\\
\mathcal V&= \frac{\left(6\right)}{\left(10\right)}=\frac{3}{5}.
\end{aligned}
$$

縞は残るが完全には暗くならない。
:::

:::example[直交状態を重ねる]
同じ強度で位相関係を一定にし、偏光だけを直交させる。

入力は第一光束の強度 $I_1=4\,\mathrm{W/m^2}$、第二光束の強度 $I_2=4\,\mathrm{W/m^2}$、二つの直線偏光の差 $\theta=\pi/2\,\mathrm{rad}$、時間的な複素相関の絶対値 $\lvert\gamma\rvert=1$ である。

二光束を別々に点灯したときの強度を足す。

$$
\begin{aligned}
S&=I_1 + I_2,\\
S/(\mathrm{W/m^2})&= \left(4\right) + \left(4\right)=8.
\end{aligned}
$$

偏光の内積と時間相関を交差項に掛ける。

$$
\begin{aligned}
I_{\rm mod}&=2 \lvert\gamma\rvert \sqrt{I_1 I_2} \left|{\cos{\left(\theta \right)}}\right|,\\
I_{\rm mod}/(\mathrm{W/m^2})&\simeq 2 \cdot \left(1\right) \cdot \sqrt{\left(4\right) \cdot \left(4\right)} \cdot \left|{\cos(\pi/2)}\right|=0.
\end{aligned}
$$

交差項が正に最大となる位相を選ぶ。

$$
\begin{aligned}
I_{\max}&=I_{\rm mod} + S,\\
I_{\max}/(\mathrm{W/m^2})&= \left(0\right) + \left(8\right)=8.
\end{aligned}
$$

位相を半周期変えて交差項を引く。

$$
\begin{aligned}
I_{\min}&=- I_{\rm mod} + S,\\
I_{\min}/(\mathrm{W/m^2})&= - \left(0\right) + \left(8\right)=8.
\end{aligned}
$$

明暗の差を和で割ると変動幅と背景の比になる。

$$
\begin{aligned}
\mathcal V&=\frac{I_{\rm mod}}{S},\\
\mathcal V&= \frac{\left(0\right)}{\left(8\right)}=0.
\end{aligned}
$$

総強度だけを測ると縞が消える。検光子による射影を加えれば縞を取り出せる。
:::

:::example[角度差と時間変動が共存する]
二束とも入る検出器で角度差と位相の揺らぎを含める。

入力は第一光束の強度 $I_1=4\,\mathrm{W/m^2}$、第二光束の強度 $I_2=4\,\mathrm{W/m^2}$、二つの直線偏光の差 $\theta\simeq 1.0472\,\mathrm{rad}$、時間的な複素相関の絶対値 $\lvert\gamma\rvert=\frac{2}{5}$ である。

二光束を別々に点灯したときの強度を足す。

$$
\begin{aligned}
S&=I_1 + I_2,\\
S/(\mathrm{W/m^2})&= \left(4\right) + \left(4\right)=8.
\end{aligned}
$$

偏光の内積と時間相関を交差項に掛ける。

$$
\begin{aligned}
I_{\rm mod}&=2 \lvert\gamma\rvert \sqrt{I_1 I_2} \left|{\cos{\left(\theta \right)}}\right|,\\
I_{\rm mod}/(\mathrm{W/m^2})&\simeq 2 \cdot \left(\frac{2}{5}\right) \cdot \sqrt{\left(4\right) \cdot \left(4\right)} \cdot \left|{\cos{\left(\left(1.047197551\right) \right)}}\right|=\frac{8}{5}.
\end{aligned}
$$

交差項が正に最大となる位相を選ぶ。

$$
\begin{aligned}
I_{\max}&=I_{\rm mod} + S,\\
I_{\max}/(\mathrm{W/m^2})&= \left(\frac{8}{5}\right) + \left(8\right)=\frac{48}{5}.
\end{aligned}
$$

位相を半周期変えて交差項を引く。

$$
\begin{aligned}
I_{\min}&=- I_{\rm mod} + S,\\
I_{\min}/(\mathrm{W/m^2})&= - \left(\frac{8}{5}\right) + \left(8\right)=\frac{32}{5}.
\end{aligned}
$$

明暗の差を和で割ると変動幅と背景の比になる。

$$
\begin{aligned}
\mathcal V&=\frac{I_{\rm mod}}{S},\\
\mathcal V&= \frac{\left(\frac{8}{5}\right)}{\left(8\right)}=\frac{1}{5}.
\end{aligned}
$$

強度の釣合いを保ったまま、偏光と時間相関の積で可視度が低下する。
:::

## 波面分割干渉（Young の干渉・Lloyd 鏡）

### 二つの開口
距離 $d$ の二開口を同相に照明し、距離 $L$ のスクリーンに投影する。開口から観測点への距離を $r_\pm=\sqrt{L^2+(x\mp d/2)^2}$ とすると $r_--r_+=2xd/(r_-+r_+)$ である。$d,|x|\ll L$ なら $r_--r_+\simeq xd/L$。よって空気中の :term[Young の干渉]{en="Young interference" yomi="やんぐのかんしょう"} の明線は $x_m=m\lambda_0 L/d$、縞間隔は $\lambda_0L/d$ である。遠方角で記述するなら、差は $d\sin\theta$ であり、大角度では $x/L$ で代用しない。各スリットに幅があると単一スリットの回折包絡線が掛かる。

### 鏡が作る第二の光源
:term[Lloyd 鏡]{en="Lloyd's mirror" yomi="ろいどきょう"} は直接波と平面鏡で反射した波を重ねる。光源の鏡面からの高さを $h$ とすると、鏡像は反対側の高さ $-h$ にある。観測点での光路差は Young の式で $d=2h$ とした値になる。ただし反射波には複素反射係数を掛ける。例えば $s$ 偏光を屈折率の高い誘電体面へかすめ入射すると、$r_s\simeq-1$ となる。位相 $\pi$ が加わるため、光路差ゼロでは二波が弱め合う。完全に暗くなるには振幅の大きさも一致する必要がある。実際の位相と振幅は偏光別に :eqref[eq-optics-fresnel] で求める。

:::example[縞幅]
$\lambda_0=500\,\mathrm{nm},d=0.50\,\mathrm{mm},L=2.0\,\mathrm m$ なら縞間隔は $2.0\,\mathrm{mm}$ である。
:::
:::example[半周期の移動]
同じ光路差で一方だけ $\pi$ の反射位相が付くと、明線は $xd/L=(m+1/2)\lambda_0$ に移る。隣接する明線の間隔は変わらない。
:::

:::example[一方だけ薄板を入れる]
平行近軸で屈折率 $n$、厚さ $t$ の板が空気の同じ厚さを置き換えると光路差は $(n-1)t$ 増す。Young の場合、縞がずれる量の絶対値は $|(n-1)t|L/d$。$n=1.50,t=6\,\mathrm{\mu m},\lambda_0=600\,\mathrm{nm}$ なら位相差は $2\pi$ の 5 倍なので 5 縞分のずれである。板による反射損失は可視度の式で別に数える。
:::

### 拡がりを持つ照明で Young 縞が薄まる
二開口への入射方向が角 $\alpha$ だけ異なると、二開口での照明の位相差に $kd\sin\alpha$ が加わる。観測角を $\theta$ とすると小角での合計位相は $kd(\alpha-\theta)$。互いに独立な入射方向が角幅 $\Delta\alpha$ の範囲に一様に分布すると、各方向からの強度を足すため、縞の交差項に $\operatorname{sinc}(kd\Delta\alpha/2)$ が掛かる。同じスクリーン上に少しずつずれた縞が足された結果である。$\Delta\alpha=\lambda_0/d$ でコントラストが初めてゼロとなるので、開口間隔を増すほど照明の角幅の制約が強くなる。

### 有限幅の二スリットで縞と包絡線を分ける

幅 $a$ の同じ二スリットを中心間隔 $d>a$ で置く。均一な平面照明の遠方振幅は、一本の開口の積分 $\int_{-a/2}^{a/2}e^{-iqx}\dd x=a\operatorname{sinc}(qa/2)$ に、二中心の位相因子 $e^{iqd/2}+e^{-iqd/2}=2\cos(qd/2)$ を掛けたもの。従って強度は $I/I(0)=\operatorname{sinc}^2(qa/2)\cos^2(qd/2)$ で、$q=k\sin\theta$。干渉因子の極大方向でも包絡線の零点と重なる縞は見えない。包絡線が傾く位置では、積の真の極大は余弦因子だけの極大からずれる。

:::example[中心の隣の縞]
細いスリットを離して置き、一次の方向を見る。

入力は一つの幅 $a\simeq 2\times 10^{-5}\,\mathrm m$、中心間隔 $d\simeq 1\times 10^{-4}\,\mathrm m$、波長 $\lambda\simeq 6\times 10^{-7}\,\mathrm m$、観測距離 $z=2\,\mathrm m$、干渉因子の極大次数 $m=1$ である。

中心間の光路差が波長の整数倍になる方向を得る。

$$
\begin{aligned}
\sin\theta_m&=\frac{\lambda m}{d},\\
\sin\theta_m&\simeq \frac{\left(6\times 10^{-7}\right) \cdot \left(1\right)}{\left(1\times 10^{-4}\right)}=\frac{3}{500}.
\end{aligned}
$$

遠方方向を、平らなスクリーン上の位置へ換算する。

$$
\begin{aligned}
x_m&=\frac{\sin\theta_m z}{\sqrt{1 - \sin\theta_m^{2}}},\\
x_m/(\mathrm m)&= \frac{\left(\frac{3}{500}\right) \cdot \left(2\right)}{\sqrt{1 - \left(\frac{3}{500}\right)^{2}}}\simeq 0.012.
\end{aligned}
$$

その次数の方向を有限幅の積分へ代入する。

$$
\begin{aligned}
v&=\frac{\pi a m}{d},\\
v&\simeq \frac{\pi \cdot \left(2\times 10^{-5}\right) \cdot \left(1\right)}{\left(1\times 10^{-4}\right)}\simeq 0.62832.
\end{aligned}
$$

余弦因子は最大値一なので、残るのは包絡線である。

$$
\begin{aligned}
I_m/I_0&=\operatorname{sinc}^{2}{\left(v \right)},\\
I_m/I_0&\simeq \operatorname{sinc}^{2}{\left(\left(0.6283185307\right) \right)}\simeq 0.87514.
\end{aligned}
$$

第一暗線の位置と干渉次数の関係を求める。

$$
\begin{aligned}
\sin\theta_{\rm env}&=\frac{\lambda}{a},\\
\sin\theta_{\rm env}&\simeq \frac{\left(6\times 10^{-7}\right)}{\left(2\times 10^{-5}\right)}=\frac{3}{100}.
\end{aligned}
$$

包絡線が一に近く、点状スリットで計算した強度が近い。
:::

:::example[太い開口で外側を見る]
スリット幅を大きくし、外側の干渉次数を選ぶ。

入力は一つの幅 $a\simeq 4\times 10^{-5}\,\mathrm m$、中心間隔 $d\simeq 1\times 10^{-4}\,\mathrm m$、波長 $\lambda\simeq 6\times 10^{-7}\,\mathrm m$、観測距離 $z=2\,\mathrm m$、干渉因子の極大次数 $m=2$ である。

中心間の光路差が波長の整数倍になる方向を得る。

$$
\begin{aligned}
\sin\theta_m&=\frac{\lambda m}{d},\\
\sin\theta_m&\simeq \frac{\left(6\times 10^{-7}\right) \cdot \left(2\right)}{\left(1\times 10^{-4}\right)}=\frac{3}{250}.
\end{aligned}
$$

遠方方向を、平らなスクリーン上の位置へ換算する。

$$
\begin{aligned}
x_m&=\frac{\sin\theta_m z}{\sqrt{1 - \sin\theta_m^{2}}},\\
x_m/(\mathrm m)&= \frac{\left(\frac{3}{250}\right) \cdot \left(2\right)}{\sqrt{1 - \left(\frac{3}{250}\right)^{2}}}\simeq 0.024002.
\end{aligned}
$$

その次数の方向を有限幅の積分へ代入する。

$$
\begin{aligned}
v&=\frac{\pi a m}{d},\\
v&\simeq \frac{\pi \cdot \left(4\times 10^{-5}\right) \cdot \left(2\right)}{\left(1\times 10^{-4}\right)}\simeq 2.5133.
\end{aligned}
$$

余弦因子は最大値一なので、残るのは包絡線である。

$$
\begin{aligned}
I_m/I_0&=\operatorname{sinc}^{2}{\left(v \right)},\\
I_m/I_0&\simeq \operatorname{sinc}^{2}{\left(\left(2.513274123\right) \right)}\simeq 0.054696.
\end{aligned}
$$

第一暗線の位置と干渉次数の関係を求める。

$$
\begin{aligned}
\sin\theta_{\rm env}&=\frac{\lambda}{a},\\
\sin\theta_{\rm env}&\simeq \frac{\left(6\times 10^{-7}\right)}{\left(4\times 10^{-5}\right)}=\frac{3}{200}.
\end{aligned}
$$

二次方向の強度は下がる。位置が明線の条件を満たすだけでは明るさは決まらない。
:::

:::example[欠ける次数]
一つの包絡線の第一零点が干渉因子の極大と揃う配置を作る。

入力は一つの幅 $a\simeq 2\times 10^{-5}\,\mathrm m$、中心間隔 $d\simeq 1\times 10^{-4}\,\mathrm m$、波長 $\lambda\simeq 6\times 10^{-7}\,\mathrm m$、観測距離 $z=2\,\mathrm m$、干渉因子の極大次数 $m=5$ である。

中心間の光路差が波長の整数倍になる方向を得る。

$$
\begin{aligned}
\sin\theta_m&=\frac{\lambda m}{d},\\
\sin\theta_m&\simeq \frac{\left(6\times 10^{-7}\right) \cdot \left(5\right)}{\left(1\times 10^{-4}\right)}=\frac{3}{100}.
\end{aligned}
$$

遠方方向を、平らなスクリーン上の位置へ換算する。

$$
\begin{aligned}
x_m&=\frac{\sin\theta_m z}{\sqrt{1 - \sin\theta_m^{2}}},\\
x_m/(\mathrm m)&= \frac{\left(\frac{3}{100}\right) \cdot \left(2\right)}{\sqrt{1 - \left(\frac{3}{100}\right)^{2}}}\simeq 0.060027.
\end{aligned}
$$

その次数の方向を有限幅の積分へ代入する。

$$
\begin{aligned}
v&=\frac{\pi a m}{d},\\
v&\simeq \frac{\pi \cdot \left(2\times 10^{-5}\right) \cdot \left(5\right)}{\left(1\times 10^{-4}\right)}\simeq 3.1416.
\end{aligned}
$$

余弦因子は最大値一なので、残るのは包絡線である。

$$
\begin{aligned}
I_m/I_0&=\operatorname{sinc}^{2}{\left(v \right)},\\
I_m/I_0&\simeq \operatorname{sinc}^{2}{\left(\left(3.141592654\right) \right)}=0.
\end{aligned}
$$

第一暗線の位置と干渉次数の関係を求める。

$$
\begin{aligned}
\sin\theta_{\rm env}&=\frac{\lambda}{a},\\
\sin\theta_{\rm env}&\simeq \frac{\left(6\times 10^{-7}\right)}{\left(2\times 10^{-5}\right)}=\frac{3}{100}.
\end{aligned}
$$

五次の位置では二スリットとも個別に振幅がゼロになる。
:::

## 振幅分割干渉（薄膜・Newton リング）

### 干渉する平行成分の位相
平行膜の屈折率を $n$、膜厚を $d$、膜内角を $\theta$ とする。入射側で重なる成分を平面波の接線波数一定で比較する。膜を一往復した成分の位相増分は
$$
\delta_{\mathrm{film}}=2 k_0 n d\cos\theta .
\label{eq-optics-film-phase}
$$
斜めの幾何学的距離は $2d/\cos\theta$ だが、膜の上側で出る位置が $2d\tan\theta$ 異なる。面上の位相差 $k_0 n\sin\theta\,(2d\tan\theta)$ を引けば上式になる。反射係数の位相差も加える必要がある。上下が無損失非磁性で垂直入射なら、低屈折率から高屈折率で反射した電場は $\pi$ を受ける。膜内で一往復余分に進む成分の振幅が十分小さければ、表面の直接反射と裏面で一度反射した二光束で近似できる。

### 平行板の等傾角縞と面の形を測る縞
膜厚 $d$ が一定なら、式 :eqref[eq-optics-film-phase] の位相は面上の位置によらず、膜内角に依存する。さまざまな方向から照明し、出てくる平行成分をレンズの後焦点面に集めれば、同じ角度の光が同じ半径に集まる。これが :term[等傾角縞]{en="fringes of equal inclination" yomi="とうけいかくじま"} である。垂直からの外角 $\theta_0$ が小さく $n_0\sin\theta_0=n\sin\theta$ なら、軸上からの位相差の変化は $-k_0 d n_0^2\theta_0^2/n$ である。従って焦点距離 $f$ のレンズで隣接する同種の環の半径二乗の差は $f^2\lambda_0 n/(d n_0^2)$ になる。鏡の反射位相の角度依存は無視した。

これに対して :term[Fizeau 干渉計]{en="Fizeau interferometer" yomi="ふぃぞーかんしょうけい"} は基準面を試料面に近づけ、ほぼ垂直に照らして両面からの反射を重ねる。隙間の厚さ $d(x,y)$ が変わると位相が $2k_0n\,d(x,y)$ だけ変わるので、隙間を結像した面には :term[等厚縞]{en="fringes of equal thickness" yomi="とうこうじま"} が現れる。空気の隙間を真空波長 $600\,\mathrm{nm}$ で測ると、隣接する明線の間で高さが $300\,\mathrm{nm}$ 変わる。小さいくさびを作るとほぼ等間隔の縞が生じ、基準面に対する局所的な高さの変化は、その縞を曲げる。

### 接触面の縞
曲率半径 $\rho$ の球面と平板の間の薄い膜で、軸からの距離を $r$ とする。膜厚は $d(r)=\rho-\sqrt{\rho^2-r^2}\simeq r^2/(2\rho)$。中心で反射の位相が反対で、ほぼ垂直なら暗線の条件は $2nd=m\lambda_0$。したがって :term[Newton リング]{en="Newton rings" yomi="にゅーとんりんぐ"} の反射の暗線は $r_m^2\simeq m\rho\lambda_0/n$。中心の二振幅が同じなら完全に暗い。半径の差でなく二乗の差が一定になる。

:::example[曲率の測定]
空気膜、$\lambda_0=600\,\mathrm{nm}$、隣接暗線の $r^2$ の差が $0.60\,\mathrm{mm^2}$ なら $\rho=\Delta r^2/\lambda_0=1.0\,\mathrm m$。
:::
:::example[液体を入れる]
空隙を $n=4/3$ の液体で満たす。反射位相の組が変わらなければ、同じ暗輪の半径は $\sqrt{3/4}$ 倍になる。
:::

### くさびと等厚干渉
法線入射で膜厚が $d(x)=d_0+\beta x$ と変化する薄いくさびを考える。傾き $\beta$ が小さく膜の各位置を平面層として扱えるなら、反射の位相差は :eqref[eq-optics-film-phase] より $2k_0n(d_0+\beta x)$ に定数の反射位相を加えたものになる。明線が次に現れるには厚さが $\lambda_0/(2n)$ だけ増す必要があるので、位置間隔は $\lambda_0/(2n|\beta|)$。Newton リングと同じく、ほぼ同じ厚さの場所を縞として地図化している。

:::example[傾きの測定]
空気くさびを $600\,\mathrm{nm}$ で照らし、縞間隔 $0.30\,\mathrm{mm}$ を測ったとする。傾きは $\lambda_0/(2\times0.30\,\mathrm{mm})=1.0\times10^{-3}$ である。
:::

### 一枚の薄膜の厳密な反射和
媒質を $0,1,2$ と並べ、$r_{ij},t_{ij}$ を $i$ 側からの接線電場の係数とする。膜の一往復位相を $\delta$ としたとき、反射は直接項 $r_{01}$ と膜を出た項の和である。後者は $t_{01}t_{10}r_{12}e^{i\delta}/(1-r_{10}r_{12}e^{i\delta})$。式 :eqref[eq-optics-fresnel] から $r_{10}=-r_{01}$、$t_{01}t_{10}=1-r_{01}^2$ なので
$$
r_{\mathrm{film}}=\frac{r_{01}+r_{12}e^{i\delta}}{1+r_{01}r_{12}e^{i\delta}}
$$
を得る。一往復あたりの絶対値が 1 未満である場合の幾何級数である。二光束近似は級数の先頭二項 $r_{01}+t_{01}t_{10}r_{12}e^{i\delta}$ を取る。落とした振幅の大きさは、$q=|r_{10}r_{12}e^{i\delta}|<1$ とおくと $|t_{01}t_{10}r_{12}e^{i\delta}|q/(1-q)$ 以下である。厳密式では $\delta=0$ のとき $\frac{r_{01}+r_{12}}{1+r_{01}r_{12}}=r_{02}$ となる。境界を二つ連続に通るが厚さがゼロなら、厚さゼロの膜の媒質 1 が結果から消え、残るのは媒質 0 と基板 2 の直接の境界である。

:::example[四分の一波長で完全に消す]
垂直入射、$n_0=1,n_1=1.5,n_2=2.25$ では $r_{01}=r_{12}=-1/5$。往復位相 $\delta=\pi$ なら分子がゼロで、反射は完全になくなる。二振幅の大きさが一致するのは膜の屈折率を幾何平均に選んだためである。
:::

### 接触しない Newton リングで膜厚を戻す

膜厚 $t(r)=t_0+R-\sqrt{R^2-r^2}$、空気に近い外部から垂直入射、反射で一度だけ位相反転がある配置を考える。膜の屈折率を $n$ とすると暗線条件は $2nt(r)=m\lambda_0$。従って曲面のたわみ $s_m=m\lambda_0/(2n)-t_0$ が $0\le s_m<R$ のとき、球面の式から $r_m^2=2Rs_m-s_m^2$ である。近軸で $s_m^2$ を捨てれば $r_{m+1}^2-r_m^2=R\lambda_0/n$。暗環の間隔は $t_0$ にほぼ依存しないが、絶対次数が分からなければ $t_0$ を $\lambda_0/(2n)$ の整数倍だけずらしても次数を同じ整数だけずらせば同じ環が得られる。この周期的な曖昧さは球面のたわみを厳密に残しても生じる。絶対厚さには接触状態または別波長を使う。

:::example[接触状態で次数を読む]
中心から暗環を数えられる配置を使う。

入力は球面の曲率半径 $R=\frac{1}{5}\,\mathrm m$、膜の屈折率 $n=1$、真空波長 $\lambda_0\simeq 6\times 10^{-7}\,\mathrm m$、中心膜厚 $t_0=0\,\mathrm m$、暗環の絶対次数 $m=10$ である。

往復光路と反射位相の暗条件を使う。

$$
\begin{aligned}
t_m&=\frac{\lambda_0 m}{2 n},\\
t_m/(\mathrm m)&\simeq \frac{\left(6\times 10^{-7}\right) \cdot \left(10\right)}{2 \cdot \left(1\right)}\simeq 3\times 10^{-6}.
\end{aligned}
$$

球面のたわみに相当する量を求める。

$$
\begin{aligned}
s_m&=t_m - t_0,\\
s_m/(\mathrm m)&\simeq \left(3\times 10^{-6}\right) - \left(0\right)\simeq 3\times 10^{-6}.
\end{aligned}
$$

球面の断面から暗環の半径を得る。

$$
\begin{aligned}
r_m&=\sqrt{2 R s_m - s_m^{2}},\\
r_m/(\mathrm m)&\simeq \sqrt{2 \cdot \left(\frac{1}{5}\right) \cdot \left(3\times 10^{-6}\right) - \left(\left(3\times 10^{-6}\right)\right)^{2}}\simeq 0.0010954.
\end{aligned}
$$

二乗半径の四次に相当する補正を省く。

$$
\begin{aligned}
r_{m,\rm par}&=\sqrt{2} \sqrt{R s_m},\\
r_{m,\rm par}/(\mathrm m)&\simeq \sqrt{2} \cdot \sqrt{\left(\frac{1}{5}\right) \cdot \left(3\times 10^{-6}\right)}\simeq 0.0010954.
\end{aligned}
$$

円の直径の測定から曲率を出す際の補正量を評価する。

$$
\begin{aligned}
\Delta r_m^2&=- s_m^{2},\\
\Delta r_m^2/(\mathrm{m^2})&\simeq - \left(\left(3\times 10^{-6}\right)\right)^{2}\simeq -9\times 10^{-12}.
\end{aligned}
$$

近軸の二乗半径は次数に比例する。
:::

:::example[中心が離れた試料]
同じ曲率で中心の隙間が何波長分かある場合を比較する。

入力は球面の曲率半径 $R=\frac{1}{5}\,\mathrm m$、膜の屈折率 $n=1$、真空波長 $\lambda_0\simeq 6\times 10^{-7}\,\mathrm m$、中心膜厚 $t_0\simeq 9\times 10^{-7}\,\mathrm m$、暗環の絶対次数 $m=10$ である。

往復光路と反射位相の暗条件を使う。

$$
\begin{aligned}
t_m&=\frac{\lambda_0 m}{2 n},\\
t_m/(\mathrm m)&\simeq \frac{\left(6\times 10^{-7}\right) \cdot \left(10\right)}{2 \cdot \left(1\right)}\simeq 3\times 10^{-6}.
\end{aligned}
$$

球面のたわみに相当する量を求める。

$$
\begin{aligned}
s_m&=t_m - t_0,\\
s_m/(\mathrm m)&\simeq \left(3\times 10^{-6}\right) - \left(9\times 10^{-7}\right)\simeq 2.1\times 10^{-6}.
\end{aligned}
$$

球面の断面から暗環の半径を得る。

$$
\begin{aligned}
r_m&=\sqrt{2 R s_m - s_m^{2}},\\
r_m/(\mathrm m)&\simeq \sqrt{2 \cdot \left(\frac{1}{5}\right) \cdot \left(2.1\times 10^{-6}\right) - \left(\left(2.1\times 10^{-6}\right)\right)^{2}}\simeq 9.1651\times 10^{-4}.
\end{aligned}
$$

二乗半径の四次に相当する補正を省く。

$$
\begin{aligned}
r_{m,\rm par}&=\sqrt{2} \sqrt{R s_m},\\
r_{m,\rm par}/(\mathrm m)&\simeq \sqrt{2} \cdot \sqrt{\left(\frac{1}{5}\right) \cdot \left(2.1\times 10^{-6}\right)}\simeq 9.1652\times 10^{-4}.
\end{aligned}
$$

円の直径の測定から曲率を出す際の補正量を評価する。

$$
\begin{aligned}
\Delta r_m^2&=- s_m^{2},\\
\Delta r_m^2/(\mathrm{m^2})&\simeq - \left(\left(2.1\times 10^{-6}\right)\right)^{2}\simeq -4.41\times 10^{-12}.
\end{aligned}
$$

同じ絶対次数の半径は小さくなる。中心暗点があっても接触しているとは限らない。
:::

:::example[液体を介した反射測定]
位相反転が一回のまま、膜に屈折率の高い液体を入れる。

入力は球面の曲率半径 $R=\frac{1}{5}\,\mathrm m$、膜の屈折率 $n=\frac{7}{5}$、真空波長 $\lambda_0\simeq 6\times 10^{-7}\,\mathrm m$、中心膜厚 $t_0=0\,\mathrm m$、暗環の絶対次数 $m=10$ である。

往復光路と反射位相の暗条件を使う。

$$
\begin{aligned}
t_m&=\frac{\lambda_0 m}{2 n},\\
t_m/(\mathrm m)&\simeq \frac{\left(6\times 10^{-7}\right) \cdot \left(10\right)}{2 \cdot \left(\frac{7}{5}\right)}\simeq 2.1429\times 10^{-6}.
\end{aligned}
$$

球面のたわみに相当する量を求める。

$$
\begin{aligned}
s_m&=t_m - t_0,\\
s_m/(\mathrm m)&\simeq \left(2.142857143\times 10^{-6}\right) - \left(0\right)\simeq 2.1429\times 10^{-6}.
\end{aligned}
$$

球面の断面から暗環の半径を得る。

$$
\begin{aligned}
r_m&=\sqrt{2 R s_m - s_m^{2}},\\
r_m/(\mathrm m)&\simeq \sqrt{2 \cdot \left(\frac{1}{5}\right) \cdot \left(2.142857143\times 10^{-6}\right) - \left(\left(2.142857143\times 10^{-6}\right)\right)^{2}}\simeq 9.2582\times 10^{-4}.
\end{aligned}
$$

二乗半径の四次に相当する補正を省く。

$$
\begin{aligned}
r_{m,\rm par}&=\sqrt{2} \sqrt{R s_m},\\
r_{m,\rm par}/(\mathrm m)&\simeq \sqrt{2} \cdot \sqrt{\left(\frac{1}{5}\right) \cdot \left(2.142857143\times 10^{-6}\right)}\simeq 9.2582\times 10^{-4}.
\end{aligned}
$$

円の直径の測定から曲率を出す際の補正量を評価する。

$$
\begin{aligned}
\Delta r_m^2&=- s_m^{2},\\
\Delta r_m^2/(\mathrm{m^2})&\simeq - \left(\left(2.142857143\times 10^{-6}\right)\right)^{2}\simeq -4.5918\times 10^{-12}.
\end{aligned}
$$

次数ごとの膜厚の差が短くなり、暗環は内側へ寄る。
:::

## Michelson 干渉計と Mach–Zehnder 干渉計

### 往復測定と一方向測定
:term[Michelson 干渉計]{en="Michelson interferometer" yomi="まいけるそんかんしょうけい"} は一つの分割器から二つの鏡へ光を往復させ、戻る波を重ねる。空気中で片方の鏡を距離 $\Delta L$ 移すと位相差は真空と同じ屈折率とみなせば $2k_0\Delta L$ だけ変わる。従って掃引した縞周期数は $2|\Delta L|/\lambda_0$ である。実際に通過した明線の整数個数には、掃引開始時の位相も関わる。長さ $d$ のセルを一方へ挿入して真空を屈折率 $n$ で置換するなら、外形を固定した往復光路差は $2(n-1)d$ になる。

:term[Mach–Zehnder 干渉計]{en="Mach–Zehnder interferometer" yomi="まっはつぇんだーかんしょうけい"} は分割と合波に別の分割器を使う。一方向のセルは $(n-1)d$ の光路差を与える。対称な無損失分割器を $B=2^{-1/2}\begin{pmatrix}1&i\\i&1\end{pmatrix}$ と取る。$B\diag(e^{i\phi_1},e^{i\phi_2})B$ を入力 $(1,0)^\mathsf T$ へ作用させると出力は $(e^{i\phi_1}-e^{i\phi_2},i[e^{i\phi_1}+e^{i\phi_2}])^\mathsf T/2$。パワーは $\sin^2[(\phi_1-\phi_2)/2]$ と $\cos^2[(\phi_1-\phi_2)/2]$ で、和は 1 である。位相差を $\pi/2$ に設定すると、片方の出力の位相に対する一次感度の大きさが最大になる。

:::example[変位の見積り]
波長 $632.8\,\mathrm{nm}$ の Michelson 干渉計で 100 本の縞が流れれば、鏡の変位は $31.64\,\mathrm{\mu m}$ である。
:::
:::example[セルの感度]
$d=1\,\mathrm{cm}$ の Mach–Zehnder のセルで $n$ が $10^{-5}$ 変わると光路長変化は $100\,\mathrm{nm}$ である。$\lambda_0=600\,\mathrm{nm}$ なら位相変化は $\pi/3$ である。
:::

### 光路差を位相として読み出す
二出力の差を全パワーで規格化すると、和側から差側を引く規約で $\cos(\phi_1-\phi_2)$ である。$\phi_1-\phi_2=\pi/2+\delta$ では出力差は $-\delta+O(\delta^3)$ になり、小さい変位の符号が読める。逆に明線や暗線の中心では微分がゼロなので、同じ出力の微小変化から変位の向きを判別しにくい。大きい掃引では余弦の周期があるので、補助の $\sin$ 信号や縞の通過数を用いて位相を連続に追う。

### 二つの周波数の縞の移り方
Michelson の鏡の移動による二周波数の位相差の差は $4\pi(\nu_1-\nu_2)\Delta L/c$。ゆっくりした可視度のうなりから周波数差を読める。各光の縞周期を別々に測れるなら、位相差の差は合成波長 $\Lambda=c/|\nu_1-\nu_2|$ の縞として解釈できる。鏡の周期変位は $\Lambda/2$ で、単色での $\lambda_0/2$ より長く取れる。小さい変位の分解能と、位相の折り返しのない範囲は同じ尺度ではない。

:::example[位相と変位を較正する]
$\lambda_0=600\,\mathrm{nm}$ の Michelson を直交点に合わせ、規格化出力差が $-0.01$ 動いたとする。微小位相は約 $+0.01\,\mathrm{rad}$、鏡変位は $\lambda_0(0.01)/(4\pi)\simeq0.48\,\mathrm{nm}$ である。直接の出力差で読む場合にも、初期の動作点が直交点にあるかを確かめる。
:::

### 二出力干渉計で入射光の揺らぎを除く

対称な無損失分割器を位相基準を吸収した実行列 $(\begin{smallmatrix}1&1\\1&-1\end{smallmatrix})/\sqrt2$ で表す。二腕の振幅の絶対値を $a,b$、相対位相を $\phi$ とすれば二出力は $(a\pm b e^{i\phi})/\sqrt2$ である。強度の和は $a^2+b^2$、差は $2ab\cos\phi$。差を和で割れば入射パワーが全体を同じ割合で揺らす成分は消える。検出器のゲインが異なる場合は先に校正する。$\phi=\pi/2$ の近くでは差の位相微分の絶対値が最大であり、小変位は差に比例して現れる。明または暗での傾きはゼロなので、その周囲の小さな揺らぎは二次でしか見えない。

:::example[明出力の極値]
両腕が釣り合い同相になる場合を求める。

入力は合流直前の総パワー $P=\frac{1}{500}\,\mathrm W$、一方の腕に残るパワー比 $f=\frac{1}{2}$、位相差 $\phi=0\,\mathrm{rad}$ である。

二腕振幅の積を総パワーで規格化する。

$$
\begin{aligned}
\mathcal V&=2 \sqrt{f \left(1 - f\right)},\\
\mathcal V&= 2 \cdot \sqrt{\left(\frac{1}{2}\right) \cdot \left(1 - \left(\frac{1}{2}\right)\right)}=1.
\end{aligned}
$$

出力の差に残る干渉因子を求める。

$$
\begin{aligned}
r&=\mathcal V \cos{\left(\phi \right)},\\
r&= \left(1\right) \cdot \cos{\left(\left(0\right) \right)}=1.
\end{aligned}
$$

和と差から第一のパワーを戻す。

$$
\begin{aligned}
P_+&=\frac{P \left(r + 1\right)}{2},\\
P_+/(\mathrm W)&= \frac{\left(\frac{1}{500}\right) \cdot \left(\left(1\right) + 1\right)}{2}=\frac{1}{500}.
\end{aligned}
$$

対応する相補的な出力も求める。

$$
\begin{aligned}
P_-&=\frac{P \left(1 - r\right)}{2},\\
P_-/(\mathrm W)&= \frac{\left(\frac{1}{500}\right) \cdot \left(1 - \left(1\right)\right)}{2}=0.
\end{aligned}
$$

位相が少し変わったときの読み出しの傾きを求める。

$$
\begin{aligned}
\dd r/\dd\phi&=- \mathcal V \sin{\left(\phi \right)},\\
\dd r/\dd\phi/(\mathrm{rad^{-1}})&= - \left(1\right) \cdot \sin{\left(\left(0\right) \right)}=0.
\end{aligned}
$$

片側に集中するが、小さな位相偏差の符号は読めない。
:::

:::example[直交位相で変位を検出]
振幅を変えず位相差を四分の一周期にする。

入力は合流直前の総パワー $P=\frac{1}{500}\,\mathrm W$、一方の腕に残るパワー比 $f=\frac{1}{2}$、位相差 $\phi=\pi/2\,\mathrm{rad}$ である。

二腕振幅の積を総パワーで規格化する。

$$
\begin{aligned}
\mathcal V&=2 \sqrt{f \left(1 - f\right)},\\
\mathcal V&= 2 \cdot \sqrt{\left(\frac{1}{2}\right) \cdot \left(1 - \left(\frac{1}{2}\right)\right)}=1.
\end{aligned}
$$

出力の差に残る干渉因子を求める。

$$
\begin{aligned}
r&=\mathcal V \cos{\left(\phi \right)},\\
r&\simeq \left(1\right) \cdot \cos(\pi/2)=0.
\end{aligned}
$$

和と差から第一のパワーを戻す。

$$
\begin{aligned}
P_+&=\frac{P \left(r + 1\right)}{2},\\
P_+/(\mathrm W)&= \frac{\left(\frac{1}{500}\right) \cdot \left(\left(0\right) + 1\right)}{2}=\frac{1}{1000}.
\end{aligned}
$$

対応する相補的な出力も求める。

$$
\begin{aligned}
P_-&=\frac{P \left(1 - r\right)}{2},\\
P_-/(\mathrm W)&= \frac{\left(\frac{1}{500}\right) \cdot \left(1 - \left(0\right)\right)}{2}=\frac{1}{1000}.
\end{aligned}
$$

位相が少し変わったときの読み出しの傾きを求める。

$$
\begin{aligned}
\dd r/\dd\phi&=- \mathcal V \sin{\left(\phi \right)},\\
\dd r/\dd\phi/(\mathrm{rad^{-1}})&\simeq - \left(1\right) \cdot \sin{\left(\left(1.570796327\right) \right)}=-1.
\end{aligned}
$$

二出力は等しいが、差は位相に対して一次で応答する。
:::

:::example[一方の腕に吸収がある]
合波前に一腕のパワーを減らす。

入力は合流直前の総パワー $P=\frac{1}{500}\,\mathrm W$、一方の腕に残るパワー比 $f=\frac{1}{10}$、位相差 $\phi=\pi/2\,\mathrm{rad}$ である。

二腕振幅の積を総パワーで規格化する。

$$
\begin{aligned}
\mathcal V&=2 \sqrt{f \left(1 - f\right)},\\
\mathcal V&= 2 \cdot \sqrt{\left(\frac{1}{10}\right) \cdot \left(1 - \left(\frac{1}{10}\right)\right)}=\frac{3}{5}.
\end{aligned}
$$

出力の差に残る干渉因子を求める。

$$
\begin{aligned}
r&=\mathcal V \cos{\left(\phi \right)},\\
r&\simeq \left(\frac{3}{5}\right) \cdot \cos(\pi/2)=0.
\end{aligned}
$$

和と差から第一のパワーを戻す。

$$
\begin{aligned}
P_+&=\frac{P \left(r + 1\right)}{2},\\
P_+/(\mathrm W)&= \frac{\left(\frac{1}{500}\right) \cdot \left(\left(0\right) + 1\right)}{2}=\frac{1}{1000}.
\end{aligned}
$$

対応する相補的な出力も求める。

$$
\begin{aligned}
P_-&=\frac{P \left(1 - r\right)}{2},\\
P_-/(\mathrm W)&= \frac{\left(\frac{1}{500}\right) \cdot \left(1 - \left(0\right)\right)}{2}=\frac{1}{1000}.
\end{aligned}
$$

位相が少し変わったときの読み出しの傾きを求める。

$$
\begin{aligned}
\dd r/\dd\phi&=- \mathcal V \sin{\left(\phi \right)},\\
\dd r/\dd\phi/(\mathrm{rad^{-1}})&\simeq - \left(\frac{3}{5}\right) \cdot \sin{\left(\left(1.570796327\right) \right)}=- \frac{3}{5}.
\end{aligned}
$$

直交位相での二出力はなお同じであるが、読み出しの感度が低い。
:::

### 時間による縞の位相読み出し

二つの腕を同じ角周波数のまま合流させると、単一画素の出力では位相差 $\phi$ と $-\phi$ を区別できない。そこで、一腕に既知の角周波数差 $\Omega>0$ を与える。第二腕の周波数を第一腕より下げ、合流面で場の位相差が $\Omega t+\phi$ になる約束とする。光周期で平均した二出力の差は、共通の規格で $D(t)=D_0\cos(\Omega t+\phi)$。ここで $D_0=2\sqrt{P_1P_2}$ は二腕が完全に相関するときの振幅である。ビートの周期 $T=2\pi/\Omega$ の整数倍 $T_m$ で積分し、
$$
\begin{aligned}
X&=\frac{2}{T_m}\int_0^{T_m}D(t)\cos(\Omega t)\dd t,\\
Y&=-\frac{2}{T_m}\int_0^{T_m}D(t)\sin(\Omega t)\dd t
\end{aligned}
$$
とする。積和公式から $X=D_0\cos\phi$、$Y=D_0\sin\phi$。従って位相は $X+iY$ の偏角である。このように周波数の違う光を合流させ、差の周波数で位相を取る方法を :term[光ヘテロダイン検出]{en="optical heterodyne detection" yomi="ひかりへてろだいんけんしゅつ"} と呼ぶ。$X,Y$ の積分に使う電気基準も光周波数をずらす信号と位相がそろっている必要がある。

有限の露光を伴う検出器で $D(t)$ を直接サンプリングする場合は、時刻 $t$ を中心に幅 $T_e$ で積分する。積分結果は
$$
\bar D(t)=D_0\operatorname{sinc}\left(\frac{\Omega T_e}{2}\right)
\cos(\Omega t+\phi) .
\label{eq-optics-heterodyne-exposure}
$$
実際、積分区間を $t+u$ に移すと、$u$ の正弦は奇関数なので積分で消え、余弦の積分を $T_e$ で割ると sinc が残る。露光時間を一周期にするとビートが消える。露光の開始時刻を時刻ラベルに使う場合は、その中心との差 $T_e/2$ に起因する位相を校正する。

:::example[順逆の変位を区別する]
真空波長 $600\,\mathrm{nm}$ の Michelson で、腕の長さをビート位相が増える向きに $75\,\mathrm{nm}$ 変える。往復位相の増分は $4\pi(75/600)=\pi/2$。移動前が $(X,Y)=(D_0,0)$ なら、移動後は $(0,D_0)$ である。逆の移動では $(0,-D_0)$ になる。強度の瞬時値を同じ位相だけで読むと、余弦の対称性によってこの二つを区別できない。
:::

:::example[露光幅でビートを失う条件]
ビートが $2\,\mathrm{kHz}$ なら周期は $0.50\,\mathrm{ms}$。露光時間が $0.25\,\mathrm{ms}$ のとき、式 :eqref[eq-optics-heterodyne-exposure] により振幅は $2/\pi$ 倍となる。露光を $0.50\,\mathrm{ms}$ に延ばすと積分内で一周期を掃引し、信号は消える。この比較は、各サンプルから振幅を戻す補正が露光幅で異なることを示す。
:::

:::example[速い位相参照でゆっくり動く鏡を追う]
初期位相をゼロとし、時間 $0.025\,\mathrm s$ の間に鏡が一様に $1.5\,\mathrm{\mu m}$ だけ移動する。真空波長 $600\,\mathrm{nm}$ なら $2(1.5\,\mathrm{\mu m})/(600\,\mathrm{nm})=5$ 周期の変化に相当する。位相を加える方向の移動なら、観測ビートには $5/(0.025\,\mathrm s)=200\,\mathrm{Hz}$ の増分が出る。静止時のビートを $2\,\mathrm{kHz}$ に選んでいれば $2.2\,\mathrm{kHz}$ になり、逆向きの移動なら $1.8\,\mathrm{kHz}$ になる。位相を連続に読む場合も、この周波数変化を標本化で失わないようにする。
:::

### 直交検出のオフセットが位相に効く仕組み

理想的な直交検出では、鏡をゆっくり動かして位相を一周させると $(X,Y)$ は原点を中心とする円を描く。検出のオフセットを $(b_x,b_y)$ とすると、規格化した軌跡は $(\cos\phi,\sin\phi)+(b_x,b_y)/D_0$ である。偏角の微分は $X,Y$ の微小変化に対し $\dd\phi=(X\dd Y-Y\dd X)/(X^2+Y^2)$。これは $X+iY$ の対数を微分して虚部を取れば得られる。従って一次の位相誤差は $(b_y\cos\phi-b_x\sin\phi)/D_0$ である。一定のオフセットでも、位相に対して周期的な誤差になる。

:::example[変位換算した周期誤差]
$X$ だけにビート振幅の $0.02$ 倍のオフセットがあるとする。位相が $\pi/2$ のとき誤差は $-0.02\,\mathrm{rad}$ の一次近似になる。$600\,\mathrm{nm}$ の Michelson では $(-0.02)\lambda_0/(4\pi)\simeq-0.955\,\mathrm{nm}$ に相当する。位相を一周させて軌跡の中心を求め、中心を引いてから偏角を取れば、このモデルのオフセットを補正できる。
:::

## 多光束干渉と Fabry–Pérot 干渉計

### 無限回の反射を足す
二枚の同一な無損失鏡の間に一様媒質を置く。両側の外部モードを同じ規格化で測り、鏡一枚のパワー反射率を $0<R<1$ とする。内部から二回反射する振幅の絶対値は $R$ である。反射位相も加えた往復位相を $\delta$ とすると、初回の透過振幅の大きさは $1-R$、$j$ 往復後はそれに $(Re^{i\delta})^j$ を掛けたものになる。幾何級数を足して
$$
\mathcal T(\delta)=\frac{(1-R)^2}{|1-Re^{i\delta}|^2}
=\frac{1}{1+\mathcal A\sin^2(\delta/2)},\qquad
\mathcal A=\frac{4R}{(1-R)^2}.
\label{eq-optics-airy-transmission}
$$
:term[Fabry–Pérot 干渉計]{en="Fabry–Pérot interferometer" yomi="ふぁぶりぺろーかんしょうけい"} では共鳴 $\delta=2\pi m$ で透過が 1 になる。これは各鏡の透過率より大きい。反射側では、直接の反射と共振器から漏れた振幅が共鳴で相殺する。

### 周波数の間隔と線幅
垂直入射で周波数による鏡位相の変化を無視すれば、$\delta=2n(\nu)(2\pi\nu)L/c+\mathrm{const.}$。$\dd\delta/\dd\nu=4\pi L n_g/c$ で、隣接共鳴間隔は局所的に $\mathrm{FSR}=c/(2Ln_g)$ である。鏡位相の周波数変化を含める場合は、微分に二枚の反射位相の微分を加える。往復位相がその範囲で単調かつほぼ線形なら、共鳴間隔は $2\pi/|\dd\delta/\dd\nu|$ である。半値条件は $\sin^2(\delta/2)=1/\mathcal A$。$\mathcal A\ge1$ のとき位相の全半値幅は $4\arcsin(1/\sqrt{\mathcal A})$。:term[フィネス]{en="finesse" yomi="ふぃねす"} は一周期 $2\pi$ をこの幅で割った値で、高反射率で $\mathcal F\simeq\pi\sqrt R/(1-R)$。$\mathcal A<1$ のときは谷が半値より高く、この半値幅を定義できない。

:::example[共鳴の間隔]
真空で $L=0.15\,\mathrm m$ のとき、$c\simeq3.0\times10^8\,\mathrm{m/s}$ から $\mathrm{FSR}=1.0\,\mathrm{GHz}$ である。
:::
:::example[二枚の鏡の比較]
$R=0.81$ では高反射率近似で $\mathcal F\simeq\pi(0.9)/0.19\simeq14.9$。$R=0.99$ なら約 313 となり、間隔は同じで線幅だけが狭くなる。
:::

### 鏡が異なるときと共振器の内部損失
内部の一方向パワー生存率を $a$、鏡の内部反射率を $R_1,R_2$ とする。往復振幅の絶対値は $a\sqrt{R_1R_2}$ である。入射鏡と出射鏡が無損失なら初回透過のパワーは $(1-R_1)(1-R_2)a$ になる。幾何級数から
$$
T=\frac{(1-R_1)(1-R_2)a}
{1+a^2R_1R_2-2a\sqrt{R_1R_2}\cos\delta}.
$$
$a=1$ の共鳴で分母から分子を引くと $(\sqrt{R_1}-\sqrt{R_2})^2$ になる。従って鏡が異なれば無損失でもピークは 1 を下回る。出射鏡と内部損失を固定し、入射鏡を変えるとき、$b=a\sqrt{R_2}$ と置いて $(1-R_1)/(1-b\sqrt{R_1})^2$ を最大にするには $\sqrt{R_1}=b$ にする。分母から $(1-R_1)(1-b^2)$ を引くと $(\sqrt{R_1}-b)^2$ なので、これが最大の証明にもなる。

:::example[異なる反射率]
$R_1=0.81,R_2=0.64,a=1$ では共鳴透過率は $(0.19)(0.36)/(1-0.72)^2\simeq0.872$ である。
:::

### 非対称な Fabry–Pérot の透過と反射

鏡のパワー反射率を $R_1,R_2$、透明な外部と同じ媒質の内部を一方向に通るパワーの割合を $q$ とする。鏡での位相を共振位相に含める。往復の振幅の大きさは $\rho=q\sqrt{R_1R_2}$、最初の透過のパワーは $(1-R_1)(1-R_2)q$ なので、共振時の透過はこれを $(1-\rho)^2$ で割った値。位相基準を選び、入射鏡の直接反射を $\sqrt{R_1}$、鏡裏からの反射を $-\sqrt{R_1}$ とすれば、反射振幅は $(\sqrt{R_1}-q\sqrt{R_2})/(1-\rho)$ になる。これは直接反射から、透過往復の先頭 $(1-R_1)q\sqrt{R_2}$ を等比和して引いた結果である。共振反射が消える条件は $R_1=q^2R_2$。内部に損失があれば、反射をゼロにしても透過が一にはならない。

:::example[対称で内部損失がない場合]
両鏡の反射率を等しくして共振させる。

入力は入射鏡の反射率 $R_1=\frac{9}{10}$、出力鏡の反射率 $R_2=\frac{9}{10}$、内部の片道透過率 $q=1$ である。

反射の平方根に片道透過率を掛けると往復の振幅率になる。

$$
\begin{aligned}
\rho&=q \sqrt{R_1 R_2},\\
\rho&= \left(1\right) \cdot \sqrt{\left(\frac{9}{10}\right) \cdot \left(\frac{9}{10}\right)}=\frac{9}{10}.
\end{aligned}
$$

直接反射と鏡間から戻る振幅を共振で足す。

$$
\begin{aligned}
r_{\rm res}&=\frac{\sqrt{R_1} - \sqrt{R_2} q}{1 - \rho},\\
r_{\rm res}&= \frac{\sqrt{\left(\frac{9}{10}\right)} - \sqrt{\left(\frac{9}{10}\right)} \cdot \left(1\right)}{1 - \left(\frac{9}{10}\right)}=0.
\end{aligned}
$$

出力が同じ入射端に戻るので絶対値の二乗を取る。

$$
\begin{aligned}
R_{\rm res}&=r_{\rm res}^{2},\\
R_{\rm res}&= \left(0\right)^{2}=0.
\end{aligned}
$$

最初に届いた振幅と繰り返しの和を含める。

$$
\begin{aligned}
T_{\rm res}&=\frac{q \left(1 - R_1\right) \left(1 - R_2\right)}{\left(1 - \rho\right)^{2}},\\
T_{\rm res}&= \frac{\left(1\right) \cdot \left(1 - \left(\frac{9}{10}\right)\right) \cdot \left(1 - \left(\frac{9}{10}\right)\right)}{\left(1 - \left(\frac{9}{10}\right)\right)^{2}}=1.
\end{aligned}
$$

両ポートで測れるパワーを収支から引く。

$$
\begin{aligned}
A&=- R_{\rm res} - T_{\rm res} + 1,\\
A&= - \left(0\right) - \left(1\right) + 1=0.
\end{aligned}
$$

直接反射が打ち消され、全パワーが透過する。
:::

:::example[非対称な二鏡]
片方の反射率を変えて共振透過の高さを見る。

入力は入射鏡の反射率 $R_1=\frac{9}{10}$、出力鏡の反射率 $R_2=\frac{99}{100}$、内部の片道透過率 $q=1$ である。

反射の平方根に片道透過率を掛けると往復の振幅率になる。

$$
\begin{aligned}
\rho&=q \sqrt{R_1 R_2},\\
\rho&= \left(1\right) \cdot \sqrt{\left(\frac{9}{10}\right) \cdot \left(\frac{99}{100}\right)}\simeq 0.94393.
\end{aligned}
$$

直接反射と鏡間から戻る振幅を共振で足す。

$$
\begin{aligned}
r_{\rm res}&=\frac{\sqrt{R_1} - \sqrt{R_2} q}{1 - \rho},\\
r_{\rm res}&\simeq \frac{\sqrt{\left(\frac{9}{10}\right)} - \sqrt{\left(\frac{99}{100}\right)} \cdot \left(1\right)}{1 - \left(0.9439279634\right)}\simeq -0.8258.
\end{aligned}
$$

出力が同じ入射端に戻るので絶対値の二乗を取る。

$$
\begin{aligned}
R_{\rm res}&=r_{\rm res}^{2},\\
R_{\rm res}&\simeq \left(-0.8257973461\right)^{2}\simeq 0.68194.
\end{aligned}
$$

最初に届いた振幅と繰り返しの和を含める。

$$
\begin{aligned}
T_{\rm res}&=\frac{q \left(1 - R_1\right) \left(1 - R_2\right)}{\left(1 - \rho\right)^{2}},\\
T_{\rm res}&\simeq \frac{\left(1\right) \cdot \left(1 - \left(\frac{9}{10}\right)\right) \cdot \left(1 - \left(\frac{99}{100}\right)\right)}{\left(1 - \left(0.9439279634\right)\right)^{2}}\simeq 0.31806.
\end{aligned}
$$

両ポートで測れるパワーを収支から引く。

$$
\begin{aligned}
A&=- R_{\rm res} - T_{\rm res} + 1,\\
A&\simeq - \left(0.6819412569\right) - \left(0.3180587431\right) + 1\simeq 8.2005\times 10^{-143}.
\end{aligned}
$$

内部での損失はゼロでも、反射が残る。
:::

:::example[内部吸収に合わせて鏡を選ぶ]
出力鏡を全反射にし、一周の吸収と入射結合を合わせる。

入力は入射鏡の反射率 $R_1=\frac{81}{100}$、出力鏡の反射率 $R_2=1$、内部の片道透過率 $q=\frac{9}{10}$ である。

反射の平方根に片道透過率を掛けると往復の振幅率になる。

$$
\begin{aligned}
\rho&=q \sqrt{R_1 R_2},\\
\rho&= \left(\frac{9}{10}\right) \cdot \sqrt{\left(\frac{81}{100}\right) \cdot \left(1\right)}=\frac{81}{100}.
\end{aligned}
$$

直接反射と鏡間から戻る振幅を共振で足す。

$$
\begin{aligned}
r_{\rm res}&=\frac{\sqrt{R_1} - \sqrt{R_2} q}{1 - \rho},\\
r_{\rm res}&= \frac{\sqrt{\left(\frac{81}{100}\right)} - \sqrt{\left(1\right)} \cdot \left(\frac{9}{10}\right)}{1 - \left(\frac{81}{100}\right)}=0.
\end{aligned}
$$

出力が同じ入射端に戻るので絶対値の二乗を取る。

$$
\begin{aligned}
R_{\rm res}&=r_{\rm res}^{2},\\
R_{\rm res}&= \left(0\right)^{2}=0.
\end{aligned}
$$

最初に届いた振幅と繰り返しの和を含める。

$$
\begin{aligned}
T_{\rm res}&=\frac{q \left(1 - R_1\right) \left(1 - R_2\right)}{\left(1 - \rho\right)^{2}},\\
T_{\rm res}&= \frac{\left(\frac{9}{10}\right) \cdot \left(1 - \left(\frac{81}{100}\right)\right) \cdot \left(1 - \left(1\right)\right)}{\left(1 - \left(\frac{81}{100}\right)\right)^{2}}=0.
\end{aligned}
$$

両ポートで測れるパワーを収支から引く。

$$
\begin{aligned}
A&=- R_{\rm res} - T_{\rm res} + 1,\\
A&= - \left(0\right) - \left(0\right) + 1=1.
\end{aligned}
$$

反射と透過がともにゼロになり、共振入力がすべて内部で失われる。
:::

## 誘電体多層膜と干渉フィルター

### 層の入力アドミタンス
垂直入射の一様無損失層で、$E=A+B$、$H=Y(A-B)$ を前面の接線電場と磁場とする。後面では $A,B$ に $e^{i\varphi},e^{-i\varphi}$ が掛かる。振幅を消去すると、後面から前面への移送は
$$
\binom{E_{\mathrm f}}{H_{\mathrm f}}=
\begin{pmatrix}\cos\varphi&-i\sin\varphi/Y\\-iY\sin\varphi&\cos\varphi\end{pmatrix}
\binom{E_{\mathrm b}}{H_{\mathrm b}},\qquad \varphi=k_0nd.
\label{eq-optics-film-matrix}
$$
行列は入射側から並べた順で積を作る。基板で $H_{\mathrm b}/E_{\mathrm b}=Y_{\mathrm s}$ なら、積の成分を $M_{jk}$ として $Y_{\mathrm in}=(M_{21}+M_{22}Y_{\mathrm s})/(M_{11}+M_{12}Y_{\mathrm s})$。入射側の $Y_0$ と :eqref[eq-optics-fresnel] を使えば反射を計算できる。

四分の一波長層では $\varphi=\pi/2$ なので $Y_{\mathrm in}=Y^2/Y_{\mathrm s}$。反射防止には $Y=\sqrt{Y_0Y_{\mathrm s}}$ を選ぶ。逆に高、低の二層を順に置くと一組で $Y_{\mathrm s}$ に $(Y_H/Y_L)^2$ が掛かる。$N$ 組なら比の $2N$ 乗になり、$Y_H>Y_L$ の順で $r$ は $-1$ に近づく。この高反射帯の中に半波長の欠陥層を設けると、二つの反射鏡間で特定の往復位相だけが整合し、狭帯域の干渉フィルターになる。実際の通過中心は欠陥層の位相と鏡の位相を併せて評価する。

:::example[反射防止膜]
空気と $n=2.25$ の基板の間では非磁性なら膜の $n=1.50$ を選ぶ。設計波長 $600\,\mathrm{nm}$ で厚さは $600/(4\times1.50)=100\,\mathrm{nm}$。
:::
:::example[設計波長での二層]
$Y_H/Y_L=2$、$Y_{\mathrm s}=Y_0$ の一組では $Y_{\mathrm in}=4Y_0$。よって $R=(3/5)^2=9/25$ になる。
:::

### 設計周波数から外れたとき
一枚の四分の一波長膜では $Y^2=Y_0Y_{\rm s}$ を満たしても、波長が変わると $\varphi=\pi/2$ の条件から外れる。:eqref[eq-optics-film-matrix] の入力式で反射振幅を書くと分子は $(Y_0-Y_{\rm s})\cos\varphi-i(Y_0Y_{\rm s}/Y-Y)\sin\varphi$ である。整合値で第二項は消えるが第一項は残る。このため設計周波数の反射零点の近くでは振幅が位相のずれに一次、パワーが二次に増える。複数層を使う設計では積の $Y_{\rm in}$ の波長依存性を調整して、広帯域で反射を抑えることができる。

### 斜め入射の通過中心
平行な損失の少ない単一欠陥層で鏡の位相を一定とみなし、通過次数を固定する。式 :eqref[eq-optics-film-phase] から設計時の中心波長の比は $\lambda(\theta_0)/\lambda(0)\simeq\cos\theta$。ここで $n_0\sin\theta_0=n\sin\theta$ であり、屈折率の波長依存を無視した。従って外部角が大きくなると中心は短波長へ移る。多層膜の各層で $s,p$ のアドミタンスが異なると、鏡の反射位相と帯域も偏光に依存するため、中心波長を求めるには層ごとの行列を評価する。

:::example[傾けて中心を調整する]
欠陥層の $n=1.5$、空気側の $\theta_0=30^\circ$ では $\sin\theta=1/3$。鏡位相を一定とすれば正面で $600\,\mathrm{nm}$ だった中心は $600\sqrt{8/9}\,\mathrm{nm}\simeq566\,\mathrm{nm}$ になる。
:::

### 周期層で反射帯ができる条件
隣接する媒質の接線電場と接線磁場をそろえれば、斜入射でも :eqref[eq-optics-film-matrix] を使える。具体的には各層の $\varphi$ を $k_z d$、$Y$ を偏光に応じた $Y_s$ または $Y_p$ とする。$s$ では磁場の成分を $-H_x$ と取り、$p$ では $H_y$ と取れば、どちらも前進波で $H=YE$ となる。

屈折率の高い層 $H$ と低い層 $L$ を周期的に並べ、一周期の行列を $M=M_HM_L$ とする。各層の行列式が 1 なので、一周期を逆向きに移る固有値 $\lambda$ は $\lambda^2-(\tr M)\lambda+1=0$ を満たす。無損失で $\varphi_H,\varphi_L$ が実なら
$$
\frac{\tr M}{2}=\cos\varphi_H\cos\varphi_L
-\frac12\left(\frac{Y_H}{Y_L}+\frac{Y_L}{Y_H}\right)\sin\varphi_H\sin\varphi_L .
$$
絶対値が 1 未満なら固有値は単位円上にある。1 より大きければ一方の絶対値が 1 未満、他方が 1 より大きく、半無限層に入る場は減衰する枝を取る。従ってその周波数区間は周期媒質の反射帯となる。有限層では終端があるため、反射率は行列の有限積から求める。

:::example[反射帯の幅]
無分散で垂直入射、設計周波数 $\nu_0$ で両層が四分の一波長とする。共通位相は $\varphi=\pi\nu/(2\nu_0)$。$m=Y_H/Y_L>1$ と置けば、帯域端では跡の半分が $-1$ なので $\sin^2\varphi=4m/(m+1)^2$ である。$\pi/2$ の周囲の二解を取ると相対周波数幅は $(4/\pi)\arcsin[(m-1)/(m+1)]$。$m=2$ なら約 $0.433$ になる。
:::

### 四分の一波長スタックの入射アドミタンス

垂直入射で第 $j$ 層のアドミタンスを $Y_j$ とする。層の厚さが媒質内波長の四分の一なら、電場と磁場の左から右への写像は共通位相を除いて $(E,H)\mapsto(H/Y_j,Y_j E)$。従って右端の終端アドミタンス $Y_s$ を左で見ると $Y_j^2/Y_s$ になる。高・低の順の一周期に二回適用すると $(Y_H/Y_L)^2Y_s$、$N$ 周期なら $(Y_H/Y_L)^{2N}Y_s$。境界で逆向きの場の磁場に負号が付くため、:eqref[eq-optics-fresnel] と同じく $r=(Y_0-Y_{\rm in})/(Y_0+Y_{\rm in})$。層数を増すと設計波長で反射率が一へ近づくが、波長を外した場合の応答には各層の厚さに対応する伝播位相を戻す必要がある。

:::example[一対の層だけを使う]
空気と基板の間に一対の高層・低層を置く。

入力は高い層の屈折率 $n_H=2$、低い層の屈折率 $n_L=\frac{7}{5}$、基板の屈折率 $n_s=\frac{3}{2}$、外側の屈折率 $n_0=1$、高・低の対の数 $N=1$ である。

非磁性媒質でアドミタンスは屈折率に比例する。

$$
\begin{aligned}
u&=\frac{n_H^{2}}{n_L^{2}},\\
u&= \frac{\left(2\right)^{2}}{\left(\frac{7}{5}\right)^{2}}=\frac{100}{49}.
\end{aligned}
$$

終端から入射側へ層の反転を繰り返す。

$$
\begin{aligned}
Y_{\rm in}/Y_0&=\frac{n_s u^{N}}{n_0},\\
Y_{\rm in}/Y_0&= \frac{\left(\frac{3}{2}\right) \cdot \left(\frac{100}{49}\right)^{\left(1\right)}}{\left(1\right)}=\frac{150}{49}.
\end{aligned}
$$

等価な界面一つの Fresnel 問題に置き換える。

$$
\begin{aligned}
r&=\frac{1 - Y_{\rm in}/Y_0}{Y_{\rm in}/Y_0 + 1},\\
r&= \frac{1 - \left(\frac{150}{49}\right)}{\left(\frac{150}{49}\right) + 1}=- \frac{101}{199}.
\end{aligned}
$$

光を同じ側へ戻すので、振幅絶対値の二乗を取る。

$$
\begin{aligned}
R&=r^{2},\\
R&= \left(- \frac{101}{199}\right)^{2}\simeq 0.25759.
\end{aligned}
$$

吸収がないので透過率は残りのパワーである。

$$
\begin{aligned}
T&=1 - R,\\
T&\simeq 1 - \left(0.2575945052\right)\simeq 0.74241.
\end{aligned}
$$

基板だけの界面と異なるアドミタンスが見えるが、まだ透過が残る。
:::

:::example[層数を増やす]
同じ材料と基板で設計波長の反射を増やす。

入力は高い層の屈折率 $n_H=2$、低い層の屈折率 $n_L=\frac{7}{5}$、基板の屈折率 $n_s=\frac{3}{2}$、外側の屈折率 $n_0=1$、高・低の対の数 $N=8$ である。

非磁性媒質でアドミタンスは屈折率に比例する。

$$
\begin{aligned}
u&=\frac{n_H^{2}}{n_L^{2}},\\
u&= \frac{\left(2\right)^{2}}{\left(\frac{7}{5}\right)^{2}}=\frac{100}{49}.
\end{aligned}
$$

終端から入射側へ層の反転を繰り返す。

$$
\begin{aligned}
Y_{\rm in}/Y_0&=\frac{n_s u^{N}}{n_0},\\
Y_{\rm in}/Y_0&= \frac{\left(\frac{3}{2}\right) \cdot \left(\frac{100}{49}\right)^{\left(8\right)}}{\left(1\right)}\simeq 451.36.
\end{aligned}
$$

等価な界面一つの Fresnel 問題に置き換える。

$$
\begin{aligned}
r&=\frac{1 - Y_{\rm in}/Y_0}{Y_{\rm in}/Y_0 + 1},\\
r&\simeq \frac{1 - \left(451.3595323\right)}{\left(451.3595323\right) + 1}\simeq -0.99558.
\end{aligned}
$$

光を同じ側へ戻すので、振幅絶対値の二乗を取る。

$$
\begin{aligned}
R&=r^{2},\\
R&\simeq \left(-0.995578738\right)^{2}\simeq 0.99118.
\end{aligned}
$$

吸収がないので透過率は残りのパワーである。

$$
\begin{aligned}
T&=1 - R,\\
T&\simeq 1 - \left(0.9911770236\right)\simeq 0.008823.
\end{aligned}
$$

入射側に換算したアドミタンスが急増し、反射が強くなる。
:::

:::example[コントラストを下げる]
層数を保ち、二材料の比だけを下げる。

入力は高い層の屈折率 $n_H=\frac{8}{5}$、低い層の屈折率 $n_L=\frac{7}{5}$、基板の屈折率 $n_s=\frac{3}{2}$、外側の屈折率 $n_0=1$、高・低の対の数 $N=8$ である。

非磁性媒質でアドミタンスは屈折率に比例する。

$$
\begin{aligned}
u&=\frac{n_H^{2}}{n_L^{2}},\\
u&= \frac{\left(\frac{8}{5}\right)^{2}}{\left(\frac{7}{5}\right)^{2}}=\frac{64}{49}.
\end{aligned}
$$

終端から入射側へ層の反転を繰り返す。

$$
\begin{aligned}
Y_{\rm in}/Y_0&=\frac{n_s u^{N}}{n_0},\\
Y_{\rm in}/Y_0&= \frac{\left(\frac{3}{2}\right) \cdot \left(\frac{64}{49}\right)^{\left(8\right)}}{\left(1\right)}\simeq 12.705.
\end{aligned}
$$

等価な界面一つの Fresnel 問題に置き換える。

$$
\begin{aligned}
r&=\frac{1 - Y_{\rm in}/Y_0}{Y_{\rm in}/Y_0 + 1},\\
r&\simeq \frac{1 - \left(12.70464139\right)}{\left(12.70464139\right) + 1}\simeq -0.85406.
\end{aligned}
$$

光を同じ側へ戻すので、振幅絶対値の二乗を取る。

$$
\begin{aligned}
R&=r^{2},\\
R&\simeq \left(-0.8540640398\right)^{2}\simeq 0.72943.
\end{aligned}
$$

吸収がないので透過率は残りのパワーである。

$$
\begin{aligned}
T&=1 - R,\\
T&\simeq 1 - \left(0.729425384\right)\simeq 0.27057.
\end{aligned}
$$

同じ層数でも反射率は下がる。層数と屈折率差の両方が効く。
:::

## 演習問題

::::exercise{level=basic}
固定位相の二波で $I_1=I_2=3$。最大と最小強度を求めよ。

:::solution
交差項の振幅は $2\sqrt{3\cdot3}=6$。平均の 6 に足すか引いて、最大 12、最小 0 である。
:::
::::

::::exercise{level=basic}
$\lambda_0=500\,\mathrm{nm},d=1\,\mathrm{mm},L=1\,\mathrm m$ の Young 縞間隔を求めよ。

:::solution
$\lambda_0 L/d=5\times10^{-4}\,\mathrm m=0.50\,\mathrm{mm}$ である。
:::
::::

::::exercise{level=basic}
波長 $600\,\mathrm{nm}$ の Michelson 干渉計で鏡を $6\,\mathrm{\mu m}$ 動かす。縞の周期を何回掃引するか。

:::solution
光路差は $2(6\,\mathrm{\mu m})=12\,\mathrm{\mu m}$。これを $0.6\,\mathrm{\mu m}$ で割ると 20 周期である。
:::
::::

::::exercise{level=basic}
同じ光路長の Mach–Zehnder で $\phi_1=\phi_2$ のとき、出力の規格化パワーを求めよ。

:::solution
差の出力は $\sin^2 0=0$、和の出力は $\cos^2 0=1$ である。
:::
::::

::::exercise{level=advanced}
$R=1/2$ の対称無損失 Fabry–Pérot の共鳴と隣接する谷での透過率を求めよ。

:::solution
:eqref[eq-optics-airy-transmission] で $\mathcal A=8$。共鳴で 1、$\delta=(2m+1)\pi$ の谷では $1/(1+8)=1/9$ である。
:::
::::

::::exercise{level=advanced}
$Y_H/Y_L=2,Y_{\mathrm s}=Y_0$ の四分の一波長層二組の設計波長での反射率を求めよ。

:::solution
入力アドミタンスは $2^4Y_0=16Y_0$。従って接線電場反射振幅は $-15/17$、パワーは $225/289$ となる。
:::
::::

::::exercise{level=advanced}
空気膜で第 4 暗輪の半径が $2\,\mathrm{mm}$ のとき第 1 暗輪の半径を求めよ。

:::solution
$r_m^2$ は $m$ に比例するので、$r_1=r_4/\sqrt4=1\,\mathrm{mm}$ である。
:::
::::

::::exercise{level=advanced}
二光束を合わせる前に強度がともに正のとき、どちらかの腕へ減衰器を入れられるとする。相関度の絶対値を保って可視度を最大にする条件を示せ。

:::solution
$\mathcal V=\abs{\Gamma_{12}}(\sqrt{I_1I_2})^{-1}\,2\sqrt{I_1I_2}/(I_1+I_2)$。相関度の因子が一定なら残りは $(\sqrt{I_1}-\sqrt{I_2})^2\ge0$ より 1 以下である。従って相関度が正のときは、強い方の腕を減衰させて $I_1=I_2>0$ にすると可視度は最大になる。相関度がゼロなら可視度は減衰の前後でゼロである。
:::
::::

::::exercise{level=advanced}
周波数 $8\,\mathrm{kHz}$ のビートを、中心時刻が既知の幅 $62.5\,\mathrm{\mu s}$ の矩形露光で記録する。ビート振幅が真値の何倍になるか求めよ。次に露光を二倍にした場合を求めよ。

:::solution
式 :eqref[eq-optics-heterodyne-exposure] の引数は $\Omega T_e/2=\pi (8000\,\mathrm{Hz})(62.5\,\mathrm{\mu s})=\pi/2$。従って振幅は $2/\pi$ 倍となり、露光の中心を時刻に使えば位相補正は要らない。露光を二倍にすると引数は $\pi$ で、ビート一周期の積分はゼロになる。
:::
::::

## 確認テスト

::::quiz
直交偏光の二波を偏光子なしで重ねるとどうなるか。

- [x] 交差項がゼロとなる。
- [ ] 位相にかかわらず暗くなる。
- [ ] 必ず完全な縞が出る。

:::explanation
複素振幅の交差項は偏光ベクトルの内積に比例する。直交で内積がゼロとなり、検出強度は二強度の和である。
:::
::::

::::quiz{answer="250" tol="1e-6"}
波長 $500\,\mathrm{nm}$ の Michelson で縞が一周期移る鏡変位を nm 単位で答えよ。

:::explanation
往復の光路差変化が波長一つで位相が $2\pi$ 進むため、変位は波長の半分である。
:::
::::

::::quiz
対称無損失な Fabry–Pérot の共鳴で透過率はどうなるか。

- [x] 1 となる。
- [ ] 0 となる。
- [ ] 常に鏡一枚の透過率と同じとなる。

:::explanation
式 :eqref[eq-optics-airy-transmission] で往復位相 $\delta=2\pi m$ なら $\sin^2(\delta/2)=0$ なので透過率は 1 になる。
:::
::::

::::quiz{answer="2" tol="1e-6"}
基板の屈折率が 4 で空気から入射する。四分の一波長の反射防止膜の屈折率を答えよ。

:::explanation
非磁性媒質ではアドミタンスが屈折率に比例する。$n^2/4=1$ を満たす正値は 2 である。
:::
::::

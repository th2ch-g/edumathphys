---
title: 荷電粒子の運動
description: Lorentz 力から回転運動とドリフトを導き、磁気ミラーを調べる。
notation:
  - symbol: '\Omega'
    meaning: 符号付きサイクロトロン角周波数
  - symbol: '\mu'
    meaning: 磁気モーメントの大きさ
  - symbol: '\vb R'
    meaning: 案内中心
  - symbol: '\vb b'
    meaning: 磁場の単位ベクトル
---

## この章の目標

- 磁場中の軌道と旋回方向を求める。
- 平均運動を生む力と適用条件を区別する。
- 断熱不変量から捕捉条件を求める。

## 一様磁場中のサイクロトロン運動

衝突が稀なとき、与えた場における試験粒子の軌道から出発できる。ここでは質量 $m$、符号付き電荷 $q$ の非相対論的粒子を扱う。電場 $\vb E$、磁場 $\vb B$ について
$$
m\dot{\vb v}=q(\vb E+\vb v\times\vb B),\qquad \dot{\vb r}=\vb v
\label{eq-plasma-lorentz-orbit}
$$
であり、速度との内積から $\dv{(mv^2/2)}{t}=q\vb E\cdot\vb v$ となる。

### 速度と位置の解

$\vb E=0$、$\vb B=B\hat{\vb z}$、$B>0$ とする。$\Omega=qB/m$ を符号付き :term[サイクロトロン角周波数]{en="cyclotron angular frequency" yomi="さいくろとろんかくしゅうはすう"} と呼ぶ。$\dot v_x=\Omega v_y,\ \dot v_y=-\Omega v_x,\ \dot v_z=0$ なので $v_x+iv_y=(v_{x0}+iv_{y0})e^{-i\Omega t}$ である。$\vb b=\vb B/B$ とし、
$$
\vb R=\vb r+\frac{\vb v\times\vb b}{\Omega}
\label{eq-plasma-guiding-position}
$$
と置く。定数場で微分すると $(\vb v\times\vb b)\times\vb b=-\vb v_\perp$ より $\dot{\vb R}=v_\parallel\vb b$ である。$\vb R$ を :term[案内中心]{en="guiding center" yomi="あんないちゅうしん"} という。半径 $\rho=v_\perp/\abs\Omega$ と一周期の進み $2\pi v_\parallel/\abs\Omega$ が軌道を決める。

:::example[旋回方向]
初速が $\hat{\vb x}$ 方向なら、正電荷は直後に $-\hat{\vb y}$ へ曲がり、電子は $+\hat{\vb y}$ へ曲がる。いずれも $\dot v_z=0$ である。
:::

:::example[陽子の半径]
$B=1\,\mathrm T$、陽子の垂直速度 $9.58\times10^5\,\mathrm{m\,s^{-1}}$ なら $\Omega\simeq9.58\times10^7\,\mathrm{s^{-1}}$ より $\rho\simeq1.00\,\mathrm{cm}$ である。
:::

### 位相空間とヘリカル軌道

初期点を $(x_0,y_0,z_0)$、初速を $(u,w,h)$ とすると、$z=z_0+ht$、横成分は
$$
\begin{aligned}
x-x_0&=\frac{u\sin\Omega t+w(1-\cos\Omega t)}{\Omega},\\
y-y_0&=\frac{w\sin\Omega t-u(1-\cos\Omega t)}{\Omega}.
\end{aligned}
$$
速度解 $v_x=u\cos\Omega t+w\sin\Omega t$、$v_y=w\cos\Omega t-u\sin\Omega t$ を時間積分した。円の中心は $(x_0+w/\Omega,y_0-u/\Omega)$ である。軌道上の実際の点が回転する向きと、複素表示の偏波ベクトルの固有値の符号を、一つの時間規約のもとで比較する必要がある。

:::example[半周後の位置]
$u=4\,\mathrm{m\,s^{-1}}$、$w=0$、$\Omega=2\,\mathrm{s^{-1}}$ とする。$t=\pi/\abs\Omega$ で $x=x_0$、$y=y_0-4\,\mathrm m$ となる。鏡像の向きに回る同じ半径の粒子では $y=y_0+4\,\mathrm m$ になる。
:::

### 軌道を三点で識別する

磁場 $B\hat{\vb z}$ に対して符号付きの $\Omega=qB/m$ を使う。$\dot v_x=\Omega v_y,\ \dot v_y=-\Omega v_x$ なので、速度平面上の回転角は $-\Omega t$ となる。初期値を $v_x(0)=v_\perp,v_y(0)=0$ とし、原点から出発すれば、積分して $x=(v_\perp/\Omega)\sin\Omega t$、$y=(v_\perp/\Omega)(\cos\Omega t-1)$ を得る。従って旋回中心の $y$ 座標は $-v_\perp/\Omega$ である。$y(t)$ は一般にその付近に留まるが $z=v_\parallel t$ は無限に伸びる。この区別は軌道を断面で観測したときに重要である。負電荷では同時刻の $x$ は同じで $y$ は反転する。以下の例では角周波数に rad を付けず、速度と位置は SI で表す。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $\Omega$ | 符号付き旋回角周波数 | $\mathrm{s^{-1}}$ |
| $v_\perp$ | 初期横速さ | $\mathrm{m\,s^{-1}}$ |
| $v_\parallel$ | 平行速度 | $\mathrm{m\,s^{-1}}$ |
| $t$ | 経過時間 | $\mathrm{s}$ |

角周波数と時間の積を $\theta$ とする。実際の速度の回転角は $-\theta$ である。

$$
\theta=\Omega t
$$

旋回半径を $\rho$ とする。半径を求めるときは電荷の符号を除く。

$$
\rho=\frac{v_{\perp}}{\left|{\Omega}\right|}
$$

初期速度方向の変位を $x$ とする。$v_x$ を積分する。

$$
x=\frac{v_{\perp} \sin{\left(\theta \right)}}{\Omega}
$$

初期速度に直交する横変位を $y$ とする。中心を原点と取り違えないよう初期値を引く。

$$
y=\frac{v_{\perp} \left(\cos{\left(\theta \right)} - 1\right)}{\Omega}
$$

磁場方向の変位を $z$ とする。Lorentz 力に平行成分はない。

$$
z=t v_{\parallel}
$$

一周あたりの進みを $\Delta z$ とする。旋回周期は正の時間である。

$$
\Delta z=\frac{2 \pi v_{\parallel}}{\left|{\Omega}\right|}
$$

:::example[陽電荷の四分の一周期]
原点から出発した後の座標を計算する。

入力値を次のように置く。

$$
\begin{aligned}
\Omega&=2\,\mathrm{s^{-1}} \\v_\perp&=4\,\mathrm{m\,s^{-1}} \\v_\parallel&=1\,\mathrm{m\,s^{-1}} \\t&=\frac{\pi}{4}\,\mathrm{s}
\end{aligned}
$$

これらを順に代入する。

角周波数と時間の積は

$$
\theta=\left(2\right) \left(\frac{\pi}{4}\right)\simeq 1.5708
$$

である。

旋回半径は

$$
\rho=\frac{\left(4\right)}{\left|{\left(2\right)}\right|}=2\,\mathrm{m}
$$

である。

初期速度方向の変位は

$$
x=\frac{\left(4\right) \sin{\left(\left(\frac{\pi}{2}\right) \right)}}{\left(2\right)}=2\,\mathrm{m}
$$

である。

初期速度に直交する横変位は

$$
y=\frac{\left(4\right) \left(\cos{\left(\left(\frac{\pi}{2}\right) \right)} - 1\right)}{\left(2\right)}=-2\,\mathrm{m}
$$

である。

磁場方向の変位は

$$
z=\left(1\right) \left(\frac{\pi}{4}\right)\simeq 0.7854\,\mathrm{m}
$$

である。

一周あたりの進みは

$$
\Delta z=\frac{2 \pi \left(1\right)}{\left|{\left(2\right)}\right|}\simeq 3.1416\,\mathrm{m}
$$

である。

初めに速度が $x$ を向けば Lorentz 力は負の $y$ を向く。
:::

:::example[負電荷の四分の一周期]
大きさの等しい反対符号の角周波数で比べる。

入力値を次のように置く。

$$
\begin{aligned}
\Omega&=-2\,\mathrm{s^{-1}} \\v_\perp&=4\,\mathrm{m\,s^{-1}} \\v_\parallel&=1\,\mathrm{m\,s^{-1}} \\t&=\frac{\pi}{4}\,\mathrm{s}
\end{aligned}
$$

これらを順に代入する。

角周波数と時間の積は

$$
\theta=\left(-2\right) \left(\frac{\pi}{4}\right)\simeq -1.5708
$$

である。

旋回半径は

$$
\rho=\frac{\left(4\right)}{\left|{\left(-2\right)}\right|}=2\,\mathrm{m}
$$

である。

初期速度方向の変位は

$$
x=\frac{\left(4\right) \sin{\left(\left(- \frac{\pi}{2}\right) \right)}}{\left(-2\right)}=2\,\mathrm{m}
$$

である。

初期速度に直交する横変位は

$$
y=\frac{\left(4\right) \left(\cos{\left(\left(- \frac{\pi}{2}\right) \right)} - 1\right)}{\left(-2\right)}=2\,\mathrm{m}
$$

である。

磁場方向の変位は

$$
z=\left(1\right) \left(\frac{\pi}{4}\right)\simeq 0.7854\,\mathrm{m}
$$

である。

一周あたりの進みは

$$
\Delta z=\frac{2 \pi \left(1\right)}{\left|{\left(-2\right)}\right|}\simeq 3.1416\,\mathrm{m}
$$

である。

軸方向の進みと半径は同じで、横方向の旋回だけが反転する。
:::

:::example[平面内の半周]
平行速度がゼロのとき軌道の最遠点を求める。

入力値を次のように置く。

$$
\begin{aligned}
\Omega&=3\,\mathrm{s^{-1}} \\v_\perp&=6\,\mathrm{m\,s^{-1}} \\v_\parallel&=0\,\mathrm{m\,s^{-1}} \\t&=\frac{\pi}{3}\,\mathrm{s}
\end{aligned}
$$

これらを順に代入する。

角周波数と時間の積は

$$
\theta=\left(3\right) \left(\frac{\pi}{3}\right)\simeq 3.1416
$$

である。

旋回半径は

$$
\rho=\frac{\left(6\right)}{\left|{\left(3\right)}\right|}=2\,\mathrm{m}
$$

である。

初期速度方向の変位は

$$
x=\frac{\left(6\right) \sin{\left(\left(\pi\right) \right)}}{\left(3\right)}=0\,\mathrm{m}
$$

である。

初期速度に直交する横変位は

$$
y=\frac{\left(6\right) \left(\cos{\left(\left(\pi\right) \right)} - 1\right)}{\left(3\right)}=-4\,\mathrm{m}
$$

である。

磁場方向の変位は

$$
z=\left(0\right) \left(\frac{\pi}{3}\right)=0\,\mathrm{m}
$$

である。

一周あたりの進みは

$$
\Delta z=\frac{2 \pi \left(0\right)}{\left|{\left(3\right)}\right|}=0\,\mathrm{m}
$$

である。

原点からの距離は半径の二倍に達する。
:::

:::example[後退するらせん]
一周期後の位置とピッチを求める。

入力値を次のように置く。

$$
\begin{aligned}
\Omega&=2\,\mathrm{s^{-1}} \\v_\perp&=1\,\mathrm{m\,s^{-1}} \\v_\parallel&=-2\,\mathrm{m\,s^{-1}} \\t&=\pi\,\mathrm{s}
\end{aligned}
$$

これらを順に代入する。

角周波数と時間の積は

$$
\theta=\left(2\right) \left(\pi\right)\simeq 6.2832
$$

である。

旋回半径は

$$
\rho=\frac{\left(1\right)}{\left|{\left(2\right)}\right|}=\frac{1}{2}\,\mathrm{m}
$$

である。

初期速度方向の変位は

$$
x=\frac{\left(1\right) \sin{\left(\left(2 \pi\right) \right)}}{\left(2\right)}=0\,\mathrm{m}
$$

である。

初期速度に直交する横変位は

$$
y=\frac{\left(1\right) \left(\cos{\left(\left(2 \pi\right) \right)} - 1\right)}{\left(2\right)}=0\,\mathrm{m}
$$

である。

磁場方向の変位は

$$
z=\left(-2\right) \left(\pi\right)\simeq -6.2832\,\mathrm{m}
$$

である。

一周あたりの進みは

$$
\Delta z=\frac{2 \pi \left(-2\right)}{\left|{\left(2\right)}\right|}\simeq -6.2832\,\mathrm{m}
$$

である。

断面上は出発点へ戻るが、三次元の位置は一致しない。
:::


## E×B ドリフトと重力ドリフト

磁場に垂直な一定の外力 $\vb F_\perp$ があるとき、ある一定速度を引けば残りは円運動になる。力を打ち消す速さを $\vb v_F$ とすると、$q\vb v_F\times\vb B=-\vb F_\perp$ から
$$
\vb v_F=\frac{\vb F\times\vb B}{qB^2}.
\label{eq-plasma-force-drift}
$$
$\vb F=q\vb E$ と置いた $\vb v_E=\vb E\times\vb B/B^2$ を :term[電場ドリフト]{en="electric drift" yomi="でんばどりふと"} と呼ぶ。一定一様場では $\vb v-\vb v_E$ に対して前節の横運動が厳密に成り立つ。平行力は $m\dot v_\parallel=F_\parallel$ として別に積分する。

重力なら $\vb F=m\vb g$ である。正負の粒子は逆向きに移動するが、それぞれの電流は $\vb j_{gs}=n_s m_s\vb g\times\vb B/B^2$ で同じ向きに加わる。電場ドリフトの電流は全種で $\rho_q\vb v_E$ となるので、準中性では小さい。

:::example[交差場での移動]
$E_x=20\,\mathrm{V\,m^{-1}}$、$B_z=0.01\,\mathrm T$ なら $\vb v_E=-2000\hat{\vb y}\,\mathrm{m\,s^{-1}}$ である。電荷の符号が逆でも速度は同じである。
:::

:::example[力の仕事]
$\vb F\cdot\vb v_F=0$ なので垂直成分の仕事の一周期平均はゼロである。平行成分があれば仕事率 $F_\parallel v_\parallel(t)$ が残る。電場に平行成分がなければ、周期後の運動エネルギーは元に戻る。
:::

### ドリフト系の保存量

$\vb E\cdot\vb B=0$ かつ両方が一定なら、$\vb w=\vb v-\vb v_E$ は $\dot{\vb w}=(q/m)\vb w\times\vb B$ であり、$m w^2/2$ は一定である。実験室での $mv^2/2$ は $m w^2/2+m v_E^2/2+m\vb w\cdot\vb v_E$ を含むため、旋回位相で振動する。電子とイオンは同じ $\vb v_E$ を持つが、$\vb w$ の旋回向きと半径は異なる。

:::example[重力と電場のつりあい]
静止した陽子の下向き重力 $\vb g=-g\hat{\vb z}$ に対して、$E_z=m_p g/e$ なら陽子は加速されない。一方、電子への電気力は下向きで、自身の重力と加算される。完全に同じ静止条件を要請するなら、各種の圧力勾配も力のつりあいに入れる。
:::

### 交差場における力のつりあいと電流

静的な磁場 $B\hat{\vb z}$ と電場 $E\hat{\vb x}$ のもとで、重力型の加速度 $g\hat{\vb x}$ を加える。平均の横方程式は $q(E+B u_y)+mg=0$ である。$u_y=-E/B-g/\Omega$ に移ると速度の残りは純粋な旋回となるので、この平均は初期速度に依存しない。正負の単一電荷で密度を等しく取ると、電場による速度は同じで電流への寄与が相殺する。一方、重力で生じた電流は $j_y=-n(m_i+m_e)g/B$ である。実際、電流による体積力 $j_yB$ は全質量に加わる外力を打ち消す。以下では $g$ の正負を含め、$u_y$ を負にも取る。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $E$ | 電場の $x$ 成分 | $\mathrm{V\,m^{-1}}$ |
| $B$ | 磁場強度 | $\mathrm{T}$ |
| $g$ | 外部加速度 | $\mathrm{m\,s^{-2}}$ |
| $\Omega_i$ | イオン角周波数 | $\mathrm{s^{-1}}$ |
| $r$ | 質量比 $m_i/m_e$ | 1 |

電場ドリフトを $u_E$ とする。$\vb E\times\vb B$ は負の $y$ 方向である。

$$
u_E=- \frac{E}{B}
$$

イオン平均速度を $u_i$ とする。イオンは正電荷とする。

$$
u_i=u_{E} - \frac{g}{\Omega_{i}}
$$

電子平均速度を $u_e$ とする。電子角周波数は $-r\Omega_i$ である。

$$
u_e=u_{E} + \frac{g}{\Omega_{i} r}
$$

電流を電荷密度で割った量を $j_y/(en)$ とする。電場に共通の速度を差し引く。

$$
j_y/(en)=- u_{e} + u_{i}
$$

電流が生じる力の規格化値を $j_yB/(en)$ とする。磁場との外積は $x$ 向きである。

$$
j_yB/(en)=B [j_y/(en)]
$$

質量平均速度を $u_{\rm m}$ とする。電子の重みは $1/(1+r)$ である。

$$
u_{\rm m}=\frac{r u_{i} + u_{e}}{r + 1}
$$

:::example[電流を運ばない移流]
両種の共通速度と電流を求める。

入力値を次のように置く。

$$
\begin{aligned}
E&=\frac{1}{10}\,\mathrm{V\,m^{-1}} \\B&=\frac{1}{100}\,\mathrm{T} \\g&=0\,\mathrm{m\,s^{-2}} \\\Omega_i&=1000000\,\mathrm{s^{-1}} \\r&=1836
\end{aligned}
$$

これらを順に代入する。

電場ドリフトは

$$
u_E=- \frac{\left(\frac{1}{10}\right)}{\left(\frac{1}{100}\right)}=-10\,\mathrm{m\,s^{-1}}
$$

である。

イオン平均速度は

$$
u_i=\left(-10\right) - \frac{\left(0\right)}{\left(1000000\right)}=-10\,\mathrm{m\,s^{-1}}
$$

である。

電子平均速度は

$$
u_e=\left(-10\right) + \frac{\left(0\right)}{\left(1000000\right) \left(1836\right)}=-10\,\mathrm{m\,s^{-1}}
$$

である。

電流を電荷密度で割った量は

$$
j_y/(en)=0\,\mathrm{m\,s^{-1}}
$$

である。

電流が生じる力の規格化値は

$$
j_yB/(en)=\left(0\right) \left(\frac{1}{100}\right)=0\,\mathrm{N\,C^{-1}}
$$

である。

質量平均速度は

$$
u_{\rm m}=\frac{\left(-10\right) \left(1836\right) + \left(-10\right)}{\left(1836\right) + 1}=-10\,\mathrm{m\,s^{-1}}
$$

である。

準中性のため、この速度で粒子は運ばれても横電流は流れない。
:::

:::example[重力型のドリフトのみ]
重い種の速度と外力をつり合わせる電流を見積もる。

入力値を次のように置く。

$$
\begin{aligned}
E&=0\,\mathrm{V\,m^{-1}} \\B&=\frac{1}{100000}\,\mathrm{T} \\g&=10\,\mathrm{m\,s^{-2}} \\\Omega_i&=1000\,\mathrm{s^{-1}} \\r&=1836
\end{aligned}
$$

これらを順に代入する。

電場ドリフトは

$$
u_E=- \frac{\left(0\right)}{\left(\frac{1}{100000}\right)}=0\,\mathrm{m\,s^{-1}}
$$

である。

イオン平均速度は

$$
u_i=\left(0\right) - \frac{\left(10\right)}{\left(1000\right)}=- \frac{1}{100}\,\mathrm{m\,s^{-1}}
$$

である。

電子平均速度は

$$
u_e=\left(0\right) + \frac{\left(10\right)}{\left(1000\right) \left(1836\right)}=\frac{1}{183600}\,\mathrm{m\,s^{-1}}
$$

である。

電流を電荷密度で割った量は

$$
j_y/(en)=\left(- \frac{1}{100}\right) - \left(\frac{1}{183600}\right)=- \frac{1837}{183600}\,\mathrm{m\,s^{-1}}
$$

である。

電流が生じる力の規格化値は

$$
j_yB/(en)=\left(- \frac{1837}{183600}\right) \left(\frac{1}{100000}\right)\simeq -1.0005\times 10^{-7}\,\mathrm{N\,C^{-1}}
$$

である。

質量平均速度は

$$
u_{\rm m}=\frac{\left(- \frac{1}{100}\right) \left(1836\right) + \left(\frac{1}{183600}\right)}{\left(1836\right) + 1}=- \frac{367}{36720}\,\mathrm{m\,s^{-1}}
$$

である。

イオンのドリフトが主な電流を担う。
:::

:::example[イオンが平均的に静止する条件]
電気力と外力がつり合うイオンを考える。

入力値を次のように置く。

$$
\begin{aligned}
E&=- \frac{1}{10000}\,\mathrm{V\,m^{-1}} \\B&=\frac{1}{100}\,\mathrm{T} \\g&=10000\,\mathrm{m\,s^{-2}} \\\Omega_i&=1000000\,\mathrm{s^{-1}} \\r&=1836
\end{aligned}
$$

これらを順に代入する。

電場ドリフトは

$$
u_E=- \frac{\left(- \frac{1}{10000}\right)}{\left(\frac{1}{100}\right)}=\frac{1}{100}\,\mathrm{m\,s^{-1}}
$$

である。

イオン平均速度は

$$
u_i=\left(\frac{1}{100}\right) - \frac{\left(10000\right)}{\left(1000000\right)}=0\,\mathrm{m\,s^{-1}}
$$

である。

電子平均速度は

$$
u_e=\left(\frac{1}{100}\right) + \frac{\left(10000\right)}{\left(1000000\right) \left(1836\right)}=\frac{1837}{183600}\,\mathrm{m\,s^{-1}}
$$

である。

電流を電荷密度で割った量は

$$
j_y/(en)=\left(0\right) - \left(\frac{1837}{183600}\right)=- \frac{1837}{183600}\,\mathrm{m\,s^{-1}}
$$

である。

電流が生じる力の規格化値は

$$
j_yB/(en)=\left(- \frac{1837}{183600}\right) \left(\frac{1}{100}\right)\simeq -0.00010005\,\mathrm{N\,C^{-1}}
$$

である。

質量平均速度は

$$
u_{\rm m}=\frac{\left(0\right) \left(1836\right) + \left(\frac{1837}{183600}\right)}{\left(1836\right) + 1}=\frac{1}{183600}\,\mathrm{m\,s^{-1}}
$$

である。

電子には両力が同方向に作用するので、同じ静止条件は満たさない。
:::

:::example[陽電子と電子の外力ドリフト]
等質量の対プラズマで質量輸送と電流を比較する。

入力値を次のように置く。

$$
\begin{aligned}
E&=0\,\mathrm{V\,m^{-1}} \\B&=\frac{1}{100}\,\mathrm{T} \\g&=10000\,\mathrm{m\,s^{-2}} \\\Omega_i&=1000000\,\mathrm{s^{-1}} \\r&=1
\end{aligned}
$$

これらを順に代入する。

電場ドリフトは

$$
u_E=- \frac{\left(0\right)}{\left(\frac{1}{100}\right)}=0\,\mathrm{m\,s^{-1}}
$$

である。

イオン平均速度は

$$
u_i=\left(0\right) - \frac{\left(10000\right)}{\left(1000000\right)}=- \frac{1}{100}\,\mathrm{m\,s^{-1}}
$$

である。

電子平均速度は

$$
u_e=\left(0\right) + \frac{\left(10000\right)}{\left(1000000\right) \left(1\right)}=\frac{1}{100}\,\mathrm{m\,s^{-1}}
$$

である。

電流を電荷密度で割った量は

$$
j_y/(en)=\left(- \frac{1}{100}\right) - \left(\frac{1}{100}\right)=- \frac{1}{50}\,\mathrm{m\,s^{-1}}
$$

である。

電流が生じる力の規格化値は

$$
j_yB/(en)=\left(- \frac{1}{50}\right) \left(\frac{1}{100}\right)=- \frac{1}{5000}\,\mathrm{N\,C^{-1}}
$$

である。

質量平均速度は

$$
u_{\rm m}=\frac{\left(- \frac{1}{100}\right) \left(1\right) + \left(\frac{1}{100}\right)}{\left(1\right) + 1}=0\,\mathrm{m\,s^{-1}}
$$

である。

質量平均流はゼロだが、反対方向の速度が電流では加わる。
:::


## 勾配 B ドリフトと曲率ドリフト

回転半径の範囲で磁場が変わると、軌道の片側で曲がり方が強くなる。磁場の向きが曲がる場合には、磁力線に沿う速度も垂直加速度を必要とする。

### 磁場の強度勾配

$\vb B=B(x)\hat{\vb z}$ とし、一旋回での相対変化 $\rho\abs{B'}/B$ を小さいとする。案内中心の周りの変位を $\vb\varrho=\vb b\times\vb v_\perp/\Omega$ とし、磁場を $B+B'\varrho_x$ と展開する。$\varrho_x=-v_y/\Omega$ なので平均横力の $x$ 成分は $qB'\langle \varrho_x v_y\rangle=-qB'v_\perp^2/(2\Omega)$ である。$y$ 成分は $\langle v_xv_y\rangle=0$ で消える。一般の局所勾配に対して横方向へ回転し、:eqref[eq-plasma-force-drift] を使えば
$$
\mu=\frac{mv_\perp^2}{2B},\qquad
\vb v_{\nabla B}=\frac{\mu}{qB}\vb b\times\nabla B.
\label{eq-plasma-gradb}
$$
$\mu$ は :term[磁気モーメント]{en="magnetic moment" yomi="じきもーめんと"} の大きさである。単位は $\mathrm{J\,T^{-1}}$ であり、旋回電流が作る双極子のベクトルは $-\mu\vb b$ である。電荷の符号を変えると旋回の向きも反転するため、電流の向きは共通である。

### 磁力線の曲率

$\vb\kappa=(\vb b\cdot\nabla)\vb b$ とする。$\vb b\cdot\vb\kappa=0$ は $\vb b^2=1$ の微分から従う。曲率加速度 $v_\parallel^2\vb\kappa$ を作るのに $q\vb v_\kappa\times\vb B=mv_\parallel^2\vb\kappa$ が必要であるため
$$
\vb v_\kappa=\frac{mv_\parallel^2}{qB}\vb b\times\vb\kappa .
\label{eq-plasma-curvature}
$$
曲率中心とは逆向きの慣性力を式 :eqref[eq-plasma-force-drift] に入れても同じ式を得る。

:::example[円形磁力線]
円筒座標で $\vb b=\hat{\vb\theta}$ なら $\vb\kappa=-\hat{\vb r}/r$ であり、正電荷の曲率ドリフトは $+\hat{\vb z}$ である。$B\propto1/r$ なら勾配ドリフトも同じ方向を向く。
:::

:::example[等方熱分布]
$\langle mv_\parallel^2\rangle=T$、$\langle mv_\perp^2/2\rangle=T$ である。$B\propto1/r$ の円形磁力線では二つの平均ドリフトの大きさはいずれも $T/(\abs qBr)$ となる。
:::

### 熱粒子に対するドリフト電流

等温で等方な各種について $\langle\mu\rangle=T_s/B$、$\langle v_\parallel^2\rangle=T_s/m_s$ である。種の平均案内中心移流だけから電流を組むと
$$
\vb j_{\rm drift}=\sum_s \frac{n_sT_s}{B}\vb b\times
\left(\nabla\ln B+\vb\kappa\right)
$$
である。等しい符号の圧力を持つ電子とイオンの寄与は加わる。ただし粒子は案内中心周りにも回転しており、場所ごとの旋回電流を合わせる必要がある。例えば一様磁場の密度勾配ではここでの $\vb j_{\rm drift}$ はゼロだが、局所の流体電流がゼロとは限らない。

:::example[ミラーで増える項]
同じ全エネルギーの捕捉粒子が反射点に近づくと、$v_\parallel^2$ とともに曲率ドリフトの係数はゼロに近づく。一方 $\mu$ は残るので、磁場の垂直勾配があれば強度勾配のドリフトは残る。
:::

### 曲率と磁場強度が連動する場

円筒座標で $\vb B=B(R)\hat{\vb\varphi}$、$B(R)=B_0(R/R_0)^{-p}$ とする。$R>0$ では $\vb\kappa=-\hat{\vb R}/R$、$\grad B=-pB\hat{\vb R}/R$ である。$\hat{\vb\varphi}\times(-\hat{\vb R})=\hat{\vb z}$ なので、曲率ドリフトは $v_\parallel^2/(\Omega R)$、強度勾配ドリフトは $p v_\perp^2/(2\Omega R)$ となる。これは各速度の二乗が独立なため、$p<0$ では打ち消しが可能である。軌道を平均するには $\rho/R\ll1$ に加え $\rho|p|/R\ll1$ が必要で、磁場の向きと強さを別に調べる。最低次の横ドリフトが打ち消されると高次項の相対的な重要性が増す。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $R$ | 円形磁力線の半径 | $\mathrm{m}$ |
| $\Omega$ | 局所の旋回角周波数 | $\mathrm{s^{-1}}$ |
| $p$ | 磁場の半径指数 | 1 |
| $v_\parallel$ | 平行速度 | $\mathrm{m\,s^{-1}}$ |
| $v_\perp$ | 横速さ | $\mathrm{m\,s^{-1}}$ |

軌道の曲率尺度比を $\rho/R$ とする。磁力線方向の変化を測る比である。

$$
\rho/R=\frac{v_\perp}{R\abs{\Omega}}
$$

強度変化の尺度比を $\rho/L_B$ とする。$L_B=B/|\partial_R B|$ を使う。

$$
\rho/L_B=(\rho/R)\abs{p}
$$

曲率ドリフトを $u_c$ とする。平行速度の符号によらない。

$$
u_c=\frac{v_{\parallel}^{2}}{R \Omega}
$$

強度勾配ドリフトを $u_g$ とする。磁気モーメントを横速さで表す。

$$
u_g=\frac{p v_{\perp}^{2}}{2 R \Omega}
$$

合成横速度を $u_z$ とする。両者は同じ軸の成分を持つ。

$$
u_z=u_{c} + u_{g}
$$

磁力線一周当たりの軸方向変位を $\Delta z$ とする。一周時間には平行速度の絶対値を使う。

$$
\Delta z=\frac{2 \pi R u_{z}}{\left|{v_{\parallel}}\right|}
$$

:::example[真空トロイダル場]
$B\propto1/R$ で二成分の大きさを比べる。

入力値を次のように置く。

$$
\begin{aligned}
R&=2\,\mathrm{m} \\\Omega&=10000\,\mathrm{s^{-1}} \\p&=1 \\v_\parallel&=100\,\mathrm{m\,s^{-1}} \\v_\perp&=100\,\mathrm{m\,s^{-1}}
\end{aligned}
$$

これらを順に代入する。

軌道の曲率尺度比は

$$
\rho/R=\frac{\left(100\right)}{\left(2\right) \left|{\left(10000\right)}\right|}=\frac{1}{200}
$$

である。

強度変化の尺度比は

$$
\rho/L_B=\left(\frac{1}{200}\right) \left|{\left(1\right)}\right|=\frac{1}{200}
$$

である。

曲率ドリフトは

$$
u_c=\frac{\left(100\right)^{2}}{\left(10000\right) \left(2\right)}=\frac{1}{2}\,\mathrm{m\,s^{-1}}
$$

である。

強度勾配ドリフトは

$$
u_g=\frac{\left(100\right)^{2} \left(1\right)}{2 \left(10000\right) \left(2\right)}=\frac{1}{4}\,\mathrm{m\,s^{-1}}
$$

である。

合成横速度は

$$
u_z=\left(\frac{1}{2}\right) + \left(\frac{1}{4}\right)=\frac{3}{4}\,\mathrm{m\,s^{-1}}
$$

である。

磁力線一周当たりの軸方向変位は

$$
\Delta z=\frac{2 \pi \left(2\right) \left(\frac{3}{4}\right)}{\left|{\left(100\right)}\right|}\simeq 0.094248\,\mathrm{m}
$$

である。

同じ速さなら曲率の寄与が強度勾配の二倍である。
:::

:::example[勾配と曲率が相殺する場]
強度が外側ほど大きい場合に合成を求める。

入力値を次のように置く。

$$
\begin{aligned}
R&=2\,\mathrm{m} \\\Omega&=10000\,\mathrm{s^{-1}} \\p&=-2 \\v_\parallel&=100\,\mathrm{m\,s^{-1}} \\v_\perp&=100\,\mathrm{m\,s^{-1}}
\end{aligned}
$$

これらを順に代入する。

軌道の曲率尺度比は

$$
\rho/R=\frac{\left(100\right)}{\left(2\right) \left|{\left(10000\right)}\right|}=\frac{1}{200}
$$

である。

強度変化の尺度比は

$$
\rho/L_B=\left(\frac{1}{200}\right) \left|{\left(-2\right)}\right|=\frac{1}{100}
$$

である。

曲率ドリフトは

$$
u_c=\frac{\left(100\right)^{2}}{\left(10000\right) \left(2\right)}=\frac{1}{2}\,\mathrm{m\,s^{-1}}
$$

である。

強度勾配ドリフトは

$$
u_g=\frac{\left(-2\right) \left(100\right)^{2}}{2 \left(10000\right) \left(2\right)}=- \frac{1}{2}\,\mathrm{m\,s^{-1}}
$$

である。

合成横速度は

$$
u_z=\left(- \frac{1}{2}\right) + \left(\frac{1}{2}\right)=0\,\mathrm{m\,s^{-1}}
$$

である。

磁力線一周当たりの軸方向変位は

$$
\Delta z=\frac{2 \pi \left(0\right) \left(2\right)}{\left|{\left(100\right)}\right|}=0\,\mathrm{m}
$$

である。

両方のドリフトがあっても一次の合成はゼロである。
:::

:::example[平行運動の向きとドリフト]
負電荷が逆向きに磁力線を通る場合を計算する。

入力値を次のように置く。

$$
\begin{aligned}
R&=1\,\mathrm{m} \\\Omega&=-100000\,\mathrm{s^{-1}} \\p&=2 \\v_\parallel&=-200\,\mathrm{m\,s^{-1}} \\v_\perp&=50\,\mathrm{m\,s^{-1}}
\end{aligned}
$$

これらを順に代入する。

軌道の曲率尺度比は

$$
\rho/R=\frac{\left(50\right)}{\left(1\right) \left|{\left(-100000\right)}\right|}=\frac{1}{2000}
$$

である。

強度変化の尺度比は

$$
\rho/L_B=\left(\frac{1}{2000}\right) \left|{\left(2\right)}\right|=\frac{1}{1000}
$$

である。

曲率ドリフトは

$$
u_c=\frac{\left(-200\right)^{2}}{\left(-100000\right) \left(1\right)}=- \frac{2}{5}\,\mathrm{m\,s^{-1}}
$$

である。

強度勾配ドリフトは

$$
u_g=\frac{\left(2\right) \left(50\right)^{2}}{2 \left(-100000\right) \left(1\right)}=- \frac{1}{40}\,\mathrm{m\,s^{-1}}
$$

である。

合成横速度は

$$
u_z=\left(- \frac{1}{40}\right) + \left(- \frac{2}{5}\right)=- \frac{17}{40}\,\mathrm{m\,s^{-1}}
$$

である。

磁力線一周当たりの軸方向変位は

$$
\Delta z=\frac{2 \pi \left(- \frac{17}{40}\right) \left(1\right)}{\left|{\left(-200\right)}\right|}\simeq -0.013352\,\mathrm{m}
$$

である。

軸向きは電荷で反転し、平行速度の符号だけで反転しない。
:::

:::example[大きな指数の磁場]
旋回半径が曲率半径より小さい場合に強度変化も調べる。

入力値を次のように置く。

$$
\begin{aligned}
R&=1\,\mathrm{m} \\\Omega&=10000\,\mathrm{s^{-1}} \\p&=20 \\v_\parallel&=200\,\mathrm{m\,s^{-1}} \\v_\perp&=100\,\mathrm{m\,s^{-1}}
\end{aligned}
$$

これらを順に代入する。

軌道の曲率尺度比は

$$
\rho/R=\frac{\left(100\right)}{\left(1\right) \left|{\left(10000\right)}\right|}=\frac{1}{100}
$$

である。

強度変化の尺度比は

$$
\rho/L_B=\left(\frac{1}{100}\right) \left|{\left(20\right)}\right|=\frac{1}{5}
$$

である。

曲率ドリフトは

$$
u_c=\frac{\left(200\right)^{2}}{\left(10000\right) \left(1\right)}=4\,\mathrm{m\,s^{-1}}
$$

である。

強度勾配ドリフトは

$$
u_g=\frac{\left(100\right)^{2} \left(20\right)}{2 \left(10000\right) \left(1\right)}=10\,\mathrm{m\,s^{-1}}
$$

である。

合成横速度は

$$
u_z=\left(10\right) + \left(4\right)=14\,\mathrm{m\,s^{-1}}
$$

である。

磁力線一周当たりの軸方向変位は

$$
\Delta z=\frac{2 \pi \left(14\right) \left(1\right)}{\left|{\left(200\right)}\right|}\simeq 0.43982\,\mathrm{m}
$$

である。

強度の尺度比は $0.2$ で、曲率の比より大きい。勾配展開の精度は大きい方の比に制約される。
:::


## 時間変化する場と分極ドリフト

$\vb E_\perp(t)$ の向きや強さが変わると、粒子はその時刻の $\vb v_E$ に少し遅れて追随する。空間一様かつ時間一定の $\vb B$ で、変化の角周波数 $\omega$ が $\abs\Omega$ より小さいとする。

旋回を平均した横速度を $\vb v_E+\delta\vb v$ と書くと、$m\dot{\vb v}_E=q\delta\vb v\times\vb B$ が最初の補正である。$\dot{\delta\vb v}$ はさらに $\omega/\abs\Omega$ だけ小さい。したがって :term[分極ドリフト]{en="polarization drift" yomi="ぶんきょくどりふと"} は
$$
\vb v_{\mathrm{pol}}=\frac{\vb b\times\dot{\vb v}_E}{\Omega}
=\frac{m}{qB^2}\partial_t\vb E_\perp .
\label{eq-plasma-polarization}
$$
空間変化を含む場合には案内中心に沿う微分や速度せん断が加わる。時間変化する磁場を指定する場合は Faraday 則が誘導電場を要求する。

:::example[正弦電場]
$E_x=E_0\cos\omega t$ なら $v_{\mathrm{pol},x}=-m\omega E_0\sin\omega t/(qB^2)$ である。$v_{E,y}$ の振幅との比は $\omega/\abs\Omega$ であり、低周波展開の小ささを直接測れる。
:::

:::example[電流と質量]
種の和を取ると $\vb j_{\mathrm{pol}}=(\rho_m/B^2)\partial_t\vb E_\perp$ である。ここで $\rho_m=\sum_s m_sn_s$ は質量密度である。電子と正イオンのドリフト方向が逆でも電流は加算される。
:::

### 誘導場があるときのエネルギー

空間的に一様な軸磁場 $B(t)$ が有限の円筒領域でゆっくり変わる。軸対称なら半径 $r$ の円に Faraday 則を用いて $2\pi r E_\theta=-\pi r^2\dot B$、従って $E_\theta=-r\dot B/2$ である。軸中心を回る粒子では $q r v_\theta=-mv_\perp^2/B$ が最低次なので、仕事率 $qv_\theta E_\theta=\mu\dot B$ となる。従って $W_\perp$ が $B$ と同じ割合で変わり、$\mu$ は最低次で一定である。同じ誘導電場による半径方向の電場ドリフトは $v_{E,r}=E_\theta/B=-r\dot B/(2B)$ だから、円の面積に磁場を掛けた $r^2B$ が一定になる。

:::example[ゆっくりした磁場の増加]
$B$ を四倍にすると、円の半径は半分、垂直エネルギーは四倍、垂直速度は二倍である。軸方向電場を加えていないので平行速度は変わらない。
:::

### 有限周波数で分極ドリフトを検査する

一様磁場、$E_x=E_0\cos\omega t$、$E_y=0$ とし、自由旋回を除く。$v_x=A\sin\omega t,\ v_y=C\cos\omega t$ と置くと $\omega A=(q/m)E_0+\Omega C$、$\omega C=\Omega A$ となる。したがって $C=-(E_0/B)/(1-\omega^2/\Omega^2)$、$A=(\omega/\Omega)C$ である。低周波で展開すると $C\simeq-E_0/B$ と $A\simeq-(\omega/\Omega)(E_0/B)$ を得る。後者は $(1/\Omega B)\partial_tE_x$ の係数である。平均電力は一周期で $\cos\omega t\sin\omega t$ が積分ゼロなので、非共鳴のこの解では消える。共鳴の分母をまたいで低周波展開を続けることはできない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $V$ | 速度単位 $E_0/B$ | $\mathrm{m\,s^{-1}}$ |
| $r$ | 符号付き周波数比 $\omega/\Omega$ | 1 |
| $\alpha$ | 位相 $\omega t$ | 1 |

有限周波数による振幅倍率を $G$ とする。横の二成分方程式を解いた分母である。

$$
G=\frac{1}{1 - r^{2}}
$$

電場と直交する速度を $v_y$ とする。低周波では静的な電場ドリフトの瞬時値に近づく。

$$
v_y=- G V \cos{\left(\alpha \right)}
$$

電場に平行な速度を $v_x$ とする。電荷符号で $r$ が変わる。

$$
v_x=- G V r \sin{\left(\alpha \right)}
$$

低周波振幅に対する相対補正を $\delta$ とする。位相に依存しない振幅を比較する。

$$
\delta=\left|{G - 1}\right|
$$

規格化横運動エネルギーを $K/(mV^2)$ とする。自由旋回を含めない特解の値である。

$$
K/(mV^2)=\frac{v_{x}^{2} + v_{y}^{2}}{2 V^{2}}
$$

瞬時電力の規格化値を $P/(qE_0 V)$ とする。速度と電場の積なので一周期の平均とは異なる。

$$
P/(qE_0 V)=\frac{v_{x} \cos{\left(\alpha \right)}}{V}
$$

:::example[低周波の四分の一周期]
電場がゼロの時刻の横応答を求める。

入力値を次のように置く。

$$
\begin{aligned}
V&=100\,\mathrm{m\,s^{-1}} \\r&=\frac{1}{20} \\\alpha&=\frac{\pi}{2}
\end{aligned}
$$

これらを順に代入する。

有限周波数による振幅倍率は

$$
G=\frac{1}{1 - \left(\frac{1}{20}\right)^{2}}=\frac{400}{399}
$$

である。

電場と直交する速度は

$$
v_y=- \left(100\right) \left(\frac{400}{399}\right) \cos{\left(\left(\frac{\pi}{2}\right) \right)}=0\,\mathrm{m\,s^{-1}}
$$

である。

電場に平行な速度は

$$
v_x=- \left(100\right) \left(\frac{1}{20}\right) \left(\frac{400}{399}\right) \sin{\left(\left(\frac{\pi}{2}\right) \right)}=- \frac{2000}{399}\,\mathrm{m\,s^{-1}}
$$

である。

低周波振幅に対する相対補正は

$$
\delta=\left|{\left(\frac{400}{399}\right) - 1}\right|=\frac{1}{399}
$$

である。

規格化横運動エネルギーは

$$
K/(mV^2)=\frac{\left(- \frac{2000}{399}\right)^{2} + \left(0\right)^{2}}{2 \left(100\right)^{2}}=\frac{200}{159201}
$$

である。

瞬時電力の規格化値は

$$
P/(qE_0 V)=\frac{\left(- \frac{2000}{399}\right) \cos{\left(\left(\frac{\pi}{2}\right) \right)}}{\left(100\right)}=0
$$

である。

この時刻には電場ドリフトはゼロだが、その変化による速度が残る。
:::

:::example[共鳴に近い周波数]
低周波式の誤差を評価する。

入力値を次のように置く。

$$
\begin{aligned}
V&=100\,\mathrm{m\,s^{-1}} \\r&=\frac{9}{10} \\\alpha&=\frac{\pi}{4}
\end{aligned}
$$

これらを順に代入する。

有限周波数による振幅倍率は

$$
G=\frac{1}{1 - \left(\frac{9}{10}\right)^{2}}=\frac{100}{19}
$$

である。

電場と直交する速度は

$$
v_y=- \left(100\right) \left(\frac{100}{19}\right) \cos{\left(\left(\frac{\pi}{4}\right) \right)}\simeq -372.16\,\mathrm{m\,s^{-1}}
$$

である。

電場に平行な速度は

$$
v_x=- \left(100\right) \left(\frac{100}{19}\right) \left(\frac{9}{10}\right) \sin{\left(\left(\frac{\pi}{4}\right) \right)}\simeq -334.95\,\mathrm{m\,s^{-1}}
$$

である。

低周波振幅に対する相対補正は

$$
\delta=\left|{\left(\frac{100}{19}\right) - 1}\right|=\frac{81}{19}
$$

である。

規格化横運動エネルギーは

$$
K/(mV^2)=\frac{\left(- \frac{4500 \sqrt{2}}{19}\right)^{2} + \left(- \frac{5000 \sqrt{2}}{19}\right)^{2}}{2 \left(100\right)^{2}}=\frac{4525}{361}
$$

である。

瞬時電力の規格化値は

$$
P/(qE_0 V)=\frac{\left(- \frac{4500 \sqrt{2}}{19}\right) \cos{\left(\left(\frac{\pi}{4}\right) \right)}}{\left(100\right)}=- \frac{45}{19}
$$

である。

振幅の差が大きく、粒子が非相対論的であっても低周波の順序が破れている。
:::

:::example[共鳴より高い周波数]
横速度が静的なドリフトと同方向か調べる。

入力値を次のように置く。

$$
\begin{aligned}
V&=100\,\mathrm{m\,s^{-1}} \\r&=2 \\\alpha&=\frac{\pi}{4}
\end{aligned}
$$

これらを順に代入する。

有限周波数による振幅倍率は

$$
G=\frac{1}{1 - \left(2\right)^{2}}=- \frac{1}{3}
$$

である。

電場と直交する速度は

$$
v_y=- \left(- \frac{1}{3}\right) \left(100\right) \cos{\left(\left(\frac{\pi}{4}\right) \right)}\simeq 23.57\,\mathrm{m\,s^{-1}}
$$

である。

電場に平行な速度は

$$
v_x=- \left(- \frac{1}{3}\right) \left(100\right) \left(2\right) \sin{\left(\left(\frac{\pi}{4}\right) \right)}\simeq 47.14\,\mathrm{m\,s^{-1}}
$$

である。

低周波振幅に対する相対補正は

$$
\delta=\left|{\left(- \frac{1}{3}\right) - 1}\right|=\frac{4}{3}
$$

である。

規格化横運動エネルギーは

$$
K/(mV^2)=\frac{\left(\frac{100 \sqrt{2}}{3}\right)^{2} + \left(\frac{50 \sqrt{2}}{3}\right)^{2}}{2 \left(100\right)^{2}}=\frac{5}{36}
$$

である。

瞬時電力の規格化値は

$$
P/(qE_0 V)=\frac{\left(\frac{100 \sqrt{2}}{3}\right) \cos{\left(\left(\frac{\pi}{4}\right) \right)}}{\left(100\right)}=\frac{1}{3}
$$

である。

分母が負に変わり、電場と同相の横速度の符号が逆転する。
:::

:::example[負電荷の分極速度]
正電荷の場合に対して電荷反転の影響を求める。

入力値を次のように置く。

$$
\begin{aligned}
V&=100\,\mathrm{m\,s^{-1}} \\r&=- \frac{1}{20} \\\alpha&=\frac{\pi}{4}
\end{aligned}
$$

これらを順に代入する。

有限周波数による振幅倍率は

$$
G=\frac{1}{1 - \left(- \frac{1}{20}\right)^{2}}=\frac{400}{399}
$$

である。

電場と直交する速度は

$$
v_y=- \left(100\right) \left(\frac{400}{399}\right) \cos{\left(\left(\frac{\pi}{4}\right) \right)}\simeq -70.888\,\mathrm{m\,s^{-1}}
$$

である。

電場に平行な速度は

$$
v_x=- \left(- \frac{1}{20}\right) \left(100\right) \left(\frac{400}{399}\right) \sin{\left(\left(\frac{\pi}{4}\right) \right)}\simeq 3.5444\,\mathrm{m\,s^{-1}}
$$

である。

低周波振幅に対する相対補正は

$$
\delta=\left|{\left(\frac{400}{399}\right) - 1}\right|=\frac{1}{399}
$$

である。

規格化横運動エネルギーは

$$
K/(mV^2)=\frac{\left(- \frac{20000 \sqrt{2}}{399}\right)^{2} + \left(\frac{1000 \sqrt{2}}{399}\right)^{2}}{2 \left(100\right)^{2}}=\frac{40100}{159201}
$$

である。

瞬時電力の規格化値は

$$
P/(qE_0 V)=\frac{\left(\frac{1000 \sqrt{2}}{399}\right) \cos{\left(\left(\frac{\pi}{4}\right) \right)}}{\left(100\right)}=\frac{10}{399}
$$

である。

電場に垂直な $v_y$ は同じで、電場に平行な $v_x$ は反転する。電流を計算するときは電荷も反転させる。
:::


## 断熱不変量と磁気ミラー

磁場が一旋回の間に少ししか変わらないとき、速い旋回を平均して遅い運動を求められる。円運動の作用は力学的運動量とベクトルポテンシャルの両方を含む。旋回の向きに積分すると $\oint m\vb v\cdot\dd\vb r=2\pi m v_\perp\rho$、$\oint q\vb A\cdot\dd\vb r=-\abs qB\pi\rho^2$ である。よって $I=(2\pi)^{-1}\oint (m\vb v+q\vb A)\cdot\dd\vb r=m\mu/\abs q$ を得る。

### 鏡力と反射条件

軸近くで $B_z=B(z)$ と置く。$\nabla\cdot\vb B=0$ と軸対称性から $B_r=-rB'/2$ が最低次である。円運動に沿う磁気力の平行成分は $q(\vb v_\perp\times\vb B_\perp)_z=-mv_\perp^2 B'/(2B)=-\mu B'$ である。場が静的で電場がなければ全運動エネルギー $W$ が保存される。よって
$$
\frac{m v_\parallel^2}{2}=W-\mu B,\qquad
B_m=\frac{B_0}{\sin^2\alpha_0}.
\label{eq-plasma-mirror}
$$
ここで $\sin^2\alpha_0=v_{\perp0}^2/v_0^2$ は基準点のピッチ角、$B_m$ は反射磁場である。$\mu$ の保存は $m v_\parallel\dot v_\parallel=-\mu v_\parallel B'$ を全エネルギー微分と比較して得る。この精度では、磁場増大とともに平行エネルギーが垂直エネルギーへ移る。$B_{\max}$ より大きな $B_m$ を必要とする角度領域を :term[損失円錐]{en="loss cone" yomi="そんしつえんすい"} という。

### 遅い運動に伴う作用

捕捉粒子が二つの反射点を往復するとき、平行作用は $J=2\int_{s_-}^{s_+}\sqrt{2m(W-\mu B(s))}\,\dd s$ である。場の変化時間が往復時間より長い場合の断熱不変量となる。さらに遅い閉じたドリフト軌道がある場合には、その一周に対応する正準作用が保たれる。軌道の長さ尺度を $L$ とすると、力学的運動量の線積分の尺度は $mvL$、ベクトルポテンシャルの寄与は $|q|BL^2$ である。前者と後者の比は旋回半径と $L$ の比の程度である。磁束がこの尺度で有限のドリフト軌道では、Stokes の定理により作用の先頭項は $q$ に磁束を掛けた量になる。従って軌道に囲まれる磁束が最低次の第三不変量である。分離曲線の横断や共鳴で周期の分離ができなくなると、これらの保存は破れる。

:::example[損失角]
$B_{\max}=4B_0$ なら境界角は $\arcsin(1/2)=30^\circ$ である。$\alpha_0=45^\circ$ の粒子は $B=2B_0$ で反射する。
:::

:::example[等方粒子の損失割合]
両側の損失円錐が半角 $\alpha_c$ のとき、一方の立体角は $2\pi(1-\cos\alpha_c)$ である。全立体角に対する割合は $1-\cos\alpha_c$ であり、上の例では約 $0.134$ である。
:::

### 放物形ミラーの作用

反射点が軸近くにあり $B(s)\simeq B_0+b s^2/2$ とする。$W_\parallel=W-\mu B_0>0$ と置くと $s_m=\sqrt{2W_\parallel/(\mu b)}$、往復角周波数は $\omega_b=\sqrt{\mu b/m}$ である。平行作用は $(s,p_\parallel)$ 平面上の楕円の面積になるから
$$
J=\pi s_m\sqrt{2mW_\parallel}=\frac{2\pi W_\parallel}{\omega_b}.
$$
往復周期は $2\pi/\omega_b$ であり、$B_0$ による一定のポテンシャル分は $J$ に入らない。断熱変化が遅い間、$J$ が一定なら $\mu$ が同じ粒子の $W_\parallel$ は $\sqrt b$ に比例する。

:::example[鏡の曲率を変える]
$B_0$ を固定し $b$ を四倍にすると、$\omega_b$ と $W_\parallel$ はともに二倍、反射距離 $s_m$ は $1/\sqrt2$ 倍になる。小振幅の往復という仮定が維持される範囲に限る。
:::

### 作用が断熱的に保たれる理由

一自由度の Hamilton 系 $H(p,s;\lambda(t))$ を考える。固定した $\lambda$ でエネルギー $W$ の周期軌道を持ち、周期 $P(W,\lambda)$ は有限、$\dot s=\partial_p H$ が反射点を除いてゼロにならないとする。作用を $J(W,\lambda)=\oint p\,\dd s$ とする。上下の枝を分けて積分すると、反射点の移動による端の項は両枝の $p$ が一致して相殺する。$H(p,s;\lambda)=W$ を $W$ で微分して $\partial_W p=(\partial_p H)^{-1}$、$\lambda$ で微分して $\partial_\lambda p=-(\partial_\lambda H)/(\partial_p H)$ である。このため
$$
\partial_W J=P,\qquad
\partial_\lambda J=-\oint\partial_\lambda H\,\dd t .
$$
時間変化する軌道では Hamilton 方程式から $\dot W=\partial_\lambda H\,\dot\lambda$ である。以上を $J(W(t),\lambda(t))$ に代入し、一周の間 $\dot\lambda$ と固定した軌道での平均を変えない最低次で評価すると、$P\langle\dot W\rangle+\partial_\lambda J\,\dot\lambda=0$ になる。一次の変化が周期平均で相殺する。この議論には周期と遅い変化の比が小さいこと、閉軌道の種類が変わらないことが必要である。分離軌道で周期が発散すればこの順序が成立しない。

:::example[振動数だけを変える調和振動子]
$H=p^2/(2m)+m\omega(t)^2s^2/2$ の周期軌道で $\langle s^2\rangle=W/(m\omega^2)$ である。$\langle\dot W\rangle=m\omega\dot\omega\langle s^2\rangle=W\dot\omega/\omega$ なので $W/\omega$ は平均で一定となる。楕円の作用 $J=2\pi W/\omega$ と同じ関係である。
:::

### 電位差を含む磁気ミラーの透過条件

静的な場に沿って $K_\parallel=W-\mu B-q\phi$ が負になれないことを使う。入口で総運動エネルギーが $K_0$、ピッチ角の二乗正弦が $s$ なら $\mu=K_0s/B_0$ である。有効ポテンシャル $V_{\rm eff}=\mu B+q\phi$ の区間最大値が $W$ より小さければ通過し、$W$ を超える障壁があれば手前で反射する。両者が等しい軌道では、到達時間 $\int\dd s/\sqrt{2(W-V_{\rm eff})/m}$ も調べる。例えば滑らかな非退化極大と $W$ が等しい場合、速度は極大点までの距離に比例してゼロに近づくため、到達時間は発散する。

入口から端までの $B$ と $q\phi$ がともに単調に増える場合、最大値は端にある。逆向きに変わる場合でも、二つの量が同じ単調な座標の一次関数なら、有効ポテンシャルもその座標の一次関数なので両端だけを比較できる。以下の例では端のエネルギーで判定できるこの条件を採用する。鏡比を $R_m=B_{\rm end}/B_0$、規格化した電位差を $F=q(\phi_{\rm end}-\phi_0)/K_0$ とすると、透過円錐の上限は $s_c=\min(1,\max(0,(1-F)/R_m))$ である。等方な入射分布の前進半球で角度の確率測度は $\sin\alpha\,\dd\alpha$、面を横切る流束の測度は $\sin\alpha\cos\alpha\,\dd\alpha$ となるため、透過粒子数と透過流束の割合は異なる。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $R_m$ | 磁気鏡比 | 1 |
| $F$ | 入口エネルギー単位の静電障壁 $q\Delta\phi/K_0$ | 1 |
| $s$ | 入射粒子の $\sin^2\alpha_0$ | 1 |

許される入射の上限を $s_c$ とする。$0\le s_c\le1$ の場合を扱う。

$$
s_c=\min\!\left(1,\max\!\left(0,\frac{1-F}{R_m}\right)\right)
$$

透過円錐の境界の余弦を $\cos\alpha_c$ とする。ピッチ角は前進半球に限る。

$$
\cos\alpha_c=\sqrt{1 - s_{c}}
$$

粒子数に対する角度分率を $P_N$ とする。等方な分布の立体角を積分する。

$$
P_N=1 - \cos\alpha_{c}
$$

一定速さの入射流束分率を $P_\Gamma$ とする。法線速度の余弦を重みに追加して積分する。

$$
P_\Gamma=s_{c}
$$

障壁端での形式的な平行エネルギーを $K_{\parallel,{\rm end}}/K_0$ とする。負の値はその位置へ進めないことを表す。

$$
K_{\parallel,{\rm end}}/K_0=- F - R_{m} s + 1
$$

障壁端まで到達できたときの横エネルギーを $K_{\perp,{\rm end}}/K_0$ とする。磁気モーメントが一定なので磁場に比例する。

$$
K_{\perp,{\rm end}}/K_0=R_{m} s
$$

:::example[電場のない鏡]
一粒子の透過と等方分布の損失率を比較する。

入力値を次のように置く。

$$
\begin{aligned}
R_m&=4 \\F&=0 \\s&=\frac{1}{8}
\end{aligned}
$$

これらを順に代入する。

許される入射の上限は

$$
s_c=\frac{1 - \left(0\right)}{\left(4\right)}=\frac{1}{4}
$$

である。

透過円錐の境界の余弦は

$$
\cos\alpha_c=\sqrt{1 - \left(\frac{1}{4}\right)}\simeq 0.86603
$$

である。

粒子数に対する角度分率は

$$
P_N=1 - \left(\frac{\sqrt{3}}{2}\right)\simeq 0.13397
$$

である。

一定速さの入射流束分率は

$$
P_\Gamma=\left(\frac{1}{4}\right)=\frac{1}{4}
$$

である。

障壁端での形式的な平行エネルギーは

$$
K_{\parallel,{\rm end}}/K_0=- \left(0\right) - \left(4\right) \left(\frac{1}{8}\right) + 1=\frac{1}{2}
$$

である。

障壁端まで到達できたときの横エネルギーは

$$
K_{\perp,{\rm end}}/K_0=\left(4\right) \left(\frac{1}{8}\right)=\frac{1}{2}
$$

である。

その粒子は透過し、等方分布は粒子数より流束で大きな透過分率を示す。
:::

:::example[正の静電障壁]
磁場を変えず正イオンを静電的に減速させる。

入力値を次のように置く。

$$
\begin{aligned}
R_m&=4 \\F&=\frac{1}{2} \\s&=\frac{1}{8}
\end{aligned}
$$

これらを順に代入する。

許される入射の上限は

$$
s_c=\frac{1 - \left(\frac{1}{2}\right)}{\left(4\right)}=\frac{1}{8}
$$

である。

透過円錐の境界の余弦は

$$
\cos\alpha_c=\sqrt{1 - \left(\frac{1}{8}\right)}\simeq 0.93541
$$

である。

粒子数に対する角度分率は

$$
P_N=1 - \left(\frac{\sqrt{14}}{4}\right)\simeq 0.064586
$$

である。

一定速さの入射流束分率は

$$
P_\Gamma=\left(\frac{1}{8}\right)=\frac{1}{8}
$$

である。

障壁端での形式的な平行エネルギーは

$$
K_{\parallel,{\rm end}}/K_0=- \left(4\right) \left(\frac{1}{8}\right) - \left(\frac{1}{2}\right) + 1=0
$$

である。

障壁端まで到達できたときの横エネルギーは

$$
K_{\perp,{\rm end}}/K_0=\left(4\right) \left(\frac{1}{8}\right)=\frac{1}{2}
$$

である。

端での形式的な平行エネルギーはゼロで、入射はしきい値にある。端への到達時間を判断するには、その近傍で有効ポテンシャルがどのように変化するかを指定する。
:::

:::example[静電加速を含む鏡]
磁場と電位が同じ比で変化する模型で透過条件を求める。

入力値を次のように置く。

$$
\begin{aligned}
R_m&=4 \\F&=- \frac{1}{2} \\s&=\frac{1}{3}
\end{aligned}
$$

これらを順に代入する。

許される入射の上限は

$$
s_c=\frac{1 - \left(- \frac{1}{2}\right)}{\left(4\right)}=\frac{3}{8}
$$

である。

透過円錐の境界の余弦は

$$
\cos\alpha_c=\sqrt{1 - \left(\frac{3}{8}\right)}\simeq 0.79057
$$

である。

粒子数に対する角度分率は

$$
P_N=1 - \left(\frac{\sqrt{10}}{4}\right)\simeq 0.20943
$$

である。

一定速さの入射流束分率は

$$
P_\Gamma=\left(\frac{3}{8}\right)=\frac{3}{8}
$$

である。

障壁端での形式的な平行エネルギーは

$$
K_{\parallel,{\rm end}}/K_0=- \left(- \frac{1}{2}\right) - \left(4\right) \left(\frac{1}{3}\right) + 1=\frac{1}{6}
$$

である。

障壁端まで到達できたときの横エネルギーは

$$
K_{\perp,{\rm end}}/K_0=\left(4\right) \left(\frac{1}{3}\right)=\frac{4}{3}
$$

である。

静電力によって透過円錐が拡大する。一般の電位分布では中間の障壁を確認する。
:::

:::example[反射する斜め入射]
入口で平行速度を持つ粒子が端に到達できるか判定する。

入力値を次のように置く。

$$
\begin{aligned}
R_m&=2 \\F&=\frac{1}{4} \\s&=\frac{1}{2}
\end{aligned}
$$

これらを順に代入する。

許される入射の上限は

$$
s_c=\frac{1 - \left(\frac{1}{4}\right)}{\left(2\right)}=\frac{3}{8}
$$

である。

透過円錐の境界の余弦は

$$
\cos\alpha_c=\sqrt{1 - \left(\frac{3}{8}\right)}\simeq 0.79057
$$

である。

粒子数に対する角度分率は

$$
P_N=1 - \left(\frac{\sqrt{10}}{4}\right)\simeq 0.20943
$$

である。

一定速さの入射流束分率は

$$
P_\Gamma=\left(\frac{3}{8}\right)=\frac{3}{8}
$$

である。

障壁端での形式的な平行エネルギーは

$$
K_{\parallel,{\rm end}}/K_0=- \left(2\right) \left(\frac{1}{2}\right) - \left(\frac{1}{4}\right) + 1=- \frac{1}{4}
$$

である。

障壁端まで到達できたときの横エネルギーは

$$
K_{\perp,{\rm end}}/K_0=\left(2\right) \left(\frac{1}{2}\right)=1
$$

である。

形式的な平行エネルギーが負なので途中で反射する。
:::


## 地球磁気圏の荷電粒子と van Allen 帯

地球磁場の第一近似として磁気双極子をとる。磁気緯度を $\lambda$、地球半径を $R_E$、赤道交点を $L R_E$ と書く。双極子の成分は $B_r:B_\lambda=2\sin\lambda:-\cos\lambda$ である。磁力線方程式 $\dd r/(r\dd\lambda)=-2\tan\lambda$ を積分し、磁場強度をとると
$$
r=L R_E\cos^2\lambda,\qquad
\frac{B(\lambda)}{B(0)}=\frac{\sqrt{1+3\sin^2\lambda}}{\cos^6\lambda}.
$$
粒子は一つの線に沿って往復しながら、勾配と曲率によって経度方向にもドリフトする。これらの運動で捕捉された高エネルギー粒子群が :term[放射線帯]{en="radiation belt" yomi="ほうしゃせんたい"} を構成する。実際の帯では、高エネルギー電子の速度に相対論的補正が必要である。静磁場中では Lorentz 力が仕事をしないので、相対論的因子 $\gamma=1/\sqrt{1-v^2/c^2}$ は一定である。一様磁場で $\vb p=\gamma m\vb v$ の運動式を使うと、角周波数は $qB/(\gamma m)$、旋回半径は $p_\perp/(|q|B)$ になる。旋回の正準作用の式にこの半径を代入すると $p_\perp^2/(2|q|B)$ を得る。静磁場でこれを断熱的に保てる場合、全運動量の大きさも一定なので、反射場の比は非相対論的な場合と同じ $1/\sin^2\alpha_0$ になる。粒子の供給、波との共鳴、損失円錐への散乱を加えなければ、分布の時間発展は求まらない。

:::example[大気との交点]
$L=2$ で $r=R_E$ を大気の境界と近似すると、$\cos^2\lambda=1/2$ より $\lambda=45^\circ$ である。場の比は $8\sqrt{5/2}\simeq12.65$ であり、赤道での損失円錐は約 $16.3^\circ$ となる。
:::

:::example[赤道ピッチ角と反射]
赤道で $\alpha_0=90^\circ$ なら $B_m=B(0)$ である。第一近似では平行速度がゼロで、赤道の鏡力もゼロなので赤道に留まる。少し傾けると二つの反射点が赤道の両側に生じる。
:::

### 双極子の赤道付近の往復

磁力線の $r(\lambda)$ を赤道付近で展開すると $r/r_0=1-\lambda^2+\order{\lambda^4}$、場は $B/B_0=1+(9/2)\lambda^2+\order{\lambda^4}$ である。線に沿う距離は $s=r_0\lambda+\order{\lambda^3 r_0}$ だから $B(s)\simeq B_0+9B_0s^2/(2r_0^2)$ である。赤道に強く捕捉された粒子は $\omega_b\simeq 3\sqrt{\mu B_0/m}/r_0$ で往復する。これは小さな緯度で反射する近似である。

:::example[三つの周期]
局所旋回時間は $\abs\Omega^{-1}$、往復時間は $\order{r_0/v}$、曲率ドリフトが一周する時間は $\order{r_0^2/(v\rho)}$ である。$r_0/\rho=10^3$ なら、係数を除いて隣同士の尺度は千倍離れる。粒子のエネルギーが高くなると $\rho$ が増え、この分離が弱まる。
:::

### 双極子の磁力線と損失円錐

双極子の磁気緯度を $\lambda$ とし、赤道で $r=L R_E$ の磁力線を考える。成分比から $\dd r/(r\dd\lambda)=-2\tan\lambda$ なので $r=L R_E\cos^2\lambda$ である。磁場強度は $r^{-3}\sqrt{1+3\sin^2\lambda}$ に比例する。吸収端を $r=A R_E$ として $A<L$ を選ぶと $\cos^2\lambda_a=A/L$、鏡比 $R_a=(L/A)^3\sqrt{4-3A/L}$ である。端までの距離は赤道面への射影でなく $\dd s=r\sqrt{1+4\tan^2\lambda}\,\dd\lambda$ で積分する。ここでは端での反射判定のみを行う。捕捉の判定から寿命を求める場合には、波や衝突でピッチ角が変わる過程も必要である。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $L$ | 赤道の半径を $R_E$ で割った値 | 1 |
| $A$ | 吸収端の半径を $R_E$ で割った値 | 1 |
| $s$ | 赤道の $\sin^2\alpha_0$ | 1 |

吸収端緯度の二乗正弦を $\sin^2\lambda_a$ とする。磁力線の形を使う。

$$
\sin^2\lambda_a=- \frac{A}{L} + 1
$$

吸収端での鏡比を $R_a$ とする。半径の効果と角因子が掛け合わさる。

$$
R_a=\frac{L^{3} \sqrt{3 \sin^2\lambda_a + 1}}{A^{3}}
$$

赤道損失円錐の境界を $s_c$ とする。$\mu$ と運動エネルギーの保存を使う。

$$
s_c=\frac{1}{R_{a}}
$$

損失円錐の半角を $\alpha_c$ とする。逆正弦を度に変換する。

$$
\alpha_c=\frac{180 \arcsin{\left(\sqrt{s_{c}} \right)}}{\pi}
$$

等方分布で両方向の円錐に入る割合を $F_N$ とする。北と南の二つの円錐を球全体で規格化する。

$$
F_N=1 - \sqrt{1 - s_{c}}
$$

端の形式的平行エネルギーを $K_{\parallel,a}/K_0$ とする。負なら端より赤道に近い位置で反射する。

$$
K_{\parallel,a}/K_0=- R_{a} s + 1
$$

:::example[赤道半径二の磁力線]
一定ピッチ角で大気への到達を調べる。

入力値を次のように置く。

$$
\begin{aligned}
L&=2 \\A&=1 \\s&=\frac{1}{4}
\end{aligned}
$$

これらを順に代入する。

吸収端緯度の二乗正弦は

$$
\sin^2\lambda_a=- \frac{\left(1\right)}{\left(2\right)} + 1=\frac{1}{2}
$$

である。

吸収端での鏡比は

$$
R_a=\frac{\left(2\right)^{3} \sqrt{3 \left(\frac{1}{2}\right) + 1}}{\left(1\right)^{3}}\simeq 12.649
$$

である。

赤道損失円錐の境界は

$$
s_c=\frac{1}{\left(4 \sqrt{10}\right)}\simeq 0.079057
$$

である。

損失円錐の半角は

$$
\alpha_c=\frac{180 \arcsin{\left(\sqrt{\left(\frac{\sqrt{10}}{40}\right)} \right)}}{\pi}\simeq 16.33\,\mathrm{deg}
$$

である。

等方分布で両方向の円錐に入る割合は

$$
F_N=1 - \sqrt{1 - \left(\frac{\sqrt{10}}{40}\right)}\simeq 0.040342
$$

である。

端の形式的平行エネルギーは

$$
K_{\parallel,a}/K_0=- \left(4 \sqrt{10}\right) \left(\frac{1}{4}\right) + 1\simeq -2.1623
$$

である。

吸収端より手前で反射するため、同じ軌道を永続するとした模型では大気に落ちない。
:::

:::example[同じ磁力線での小ピッチ角]
平行に近い粒子の行先を判定する。

入力値を次のように置く。

$$
\begin{aligned}
L&=2 \\A&=1 \\s&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

吸収端緯度の二乗正弦は

$$
\sin^2\lambda_a=- \frac{\left(1\right)}{\left(2\right)} + 1=\frac{1}{2}
$$

である。

吸収端での鏡比は

$$
R_a=\frac{\left(2\right)^{3} \sqrt{3 \left(\frac{1}{2}\right) + 1}}{\left(1\right)^{3}}\simeq 12.649
$$

である。

赤道損失円錐の境界は

$$
s_c=\frac{1}{\left(4 \sqrt{10}\right)}\simeq 0.079057
$$

である。

損失円錐の半角は

$$
\alpha_c=\frac{180 \arcsin{\left(\sqrt{\left(\frac{\sqrt{10}}{40}\right)} \right)}}{\pi}\simeq 16.33\,\mathrm{deg}
$$

である。

等方分布で両方向の円錐に入る割合は

$$
F_N=1 - \sqrt{1 - \left(\frac{\sqrt{10}}{40}\right)}\simeq 0.040342
$$

である。

端の形式的平行エネルギーは

$$
K_{\parallel,a}/K_0=- \left(4 \sqrt{10}\right) \left(\frac{1}{100}\right) + 1\simeq 0.87351
$$

である。

小さな磁気モーメントでは強い端の磁場でも反射に必要な横運動エネルギーに達しない。
:::

:::example[高い吸収高度]
吸収面を地表より上に設定したときの円錐を求める。

入力値を次のように置く。

$$
\begin{aligned}
L&=3 \\A&=\frac{6}{5} \\s&=\frac{1}{20}
\end{aligned}
$$

これらを順に代入する。

吸収端緯度の二乗正弦は

$$
\sin^2\lambda_a=1 - \frac{\left(\frac{6}{5}\right)}{\left(3\right)}=\frac{3}{5}
$$

である。

吸収端での鏡比は

$$
R_a=\frac{\left(3\right)^{3} \sqrt{3 \left(\frac{3}{5}\right) + 1}}{\left(\frac{6}{5}\right)^{3}}\simeq 26.146
$$

である。

赤道損失円錐の境界は

$$
s_c=\frac{1}{\left(\frac{25 \sqrt{70}}{8}\right)}\simeq 0.038247
$$

である。

損失円錐の半角は

$$
\alpha_c=\frac{180 \arcsin{\left(\sqrt{\left(\frac{4 \sqrt{70}}{875}\right)} \right)}}{\pi}\simeq 11.278\,\mathrm{deg}
$$

である。

等方分布で両方向の円錐に入る割合は

$$
F_N=1 - \sqrt{1 - \left(\frac{4 \sqrt{70}}{875}\right)}\simeq 0.01931
$$

である。

端の形式的平行エネルギーは

$$
K_{\parallel,a}/K_0=- \left(\frac{1}{20}\right) \left(\frac{25 \sqrt{70}}{8}\right) + 1\simeq -0.30728
$$

である。

同じ $L$ で $A=1$ とした場合より吸収面での磁場強度が弱く、損失円錐は広い。
:::

:::example[外側の磁力線]
赤道が外側にある場合、二の磁力線と同じ角度で比べる。

入力値を次のように置く。

$$
\begin{aligned}
L&=6 \\A&=1 \\s&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

吸収端緯度の二乗正弦は

$$
\sin^2\lambda_a=- \frac{\left(1\right)}{\left(6\right)} + 1=\frac{5}{6}
$$

である。

吸収端での鏡比は

$$
R_a=\frac{\left(6\right)^{3} \sqrt{3 \left(\frac{5}{6}\right) + 1}}{\left(1\right)^{3}}\simeq 404.1
$$

である。

赤道損失円錐の境界は

$$
s_c=\frac{1}{\left(108 \sqrt{14}\right)}\simeq 0.0024746
$$

である。

損失円錐の半角は

$$
\alpha_c=\frac{180 \arcsin{\left(\sqrt{\left(\frac{\sqrt{14}}{1512}\right)} \right)}}{\pi}\simeq 2.8514\,\mathrm{deg}
$$

である。

等方分布で両方向の円錐に入る割合は

$$
F_N=1 - \sqrt{1 - \left(\frac{\sqrt{14}}{1512}\right)}\simeq 0.0012381
$$

である。

端の形式的平行エネルギーは

$$
K_{\parallel,a}/K_0=- \left(108 \sqrt{14}\right) \left(\frac{1}{100}\right) + 1\simeq -3.041
$$

である。

赤道に対する端の磁場の比が大きいため、二の磁力線で失われる粒子にも反射条件が成り立つ。
:::


## 案内中心近似の系統的導出

軌道を分解する順序を明示すると、何を捨てたかが判断できる。$\rho/L\sim\omega_{\mathrm{slow}}/\abs\Omega=\epsilon\ll1$ とし、磁気零点から離れ、$v_\parallel$ と $v_\perp$ が同程度で、ドリフト速度はその代表速度 $v$ の $\order{\epsilon}$ とする。$L$ には磁場の強度と方向の変化尺度のうち短い方を用いる。場と平行速度を速い位相 $\theta$ に対して固定し、$\vb r=\vb R+\vb\varrho$、$\langle\vb\varrho\rangle_\theta=0$ と置く。

### 一次の平均方程式

最低次の速度 $\vb w$ は $\dot{\vb w}=\Omega\vb w\times\vb b$ を満たす。$\vb\varrho=\vb b\times\vb w/\Omega$ と置けば速い時間微分が $\vb w$ となる。直交基底で $\langle w_iw_j\rangle=w^2(\delta_{ij}-b_i b_j)/2$ である。したがって
$$
q\langle\vb w\times[(\vb\varrho\cdot\nabla)\vb B]\rangle
=-\mu\nabla B .
$$
この恒等式は $\vb b=\hat{\vb z}$ とした局所座標で確認できる。$\langle w_x\varrho_y\rangle=w^2/(2\Omega)$、$\langle w_y\varrho_x\rangle=-w^2/(2\Omega)$ なので、平均力は $q w^2(-\partial_x B_z,-\partial_y B_z,\partial_x B_x+\partial_y B_y)/(2\Omega)$ である。$\nabla\cdot\vb B=0$ を使えば上式となる。

平均加速度の垂直成分は最低次の平均速度 $v_\parallel\vb b$ から $v_\parallel^2\vb\kappa$ である。静磁場、弱い電場の順序では以上を統合すると
$$
\begin{aligned}
\dot{\vb R}&=v_\parallel\vb b+\vb v_E+\vb v_{\nabla B}+\vb v_\kappa+\order{\epsilon^2v},\\
m\dot v_\parallel&=qE_\parallel-\mu\partial_\parallel B+\order{\epsilon^2 m\abs\Omega v}.
\end{aligned}
\label{eq-plasma-gc-first}
$$
垂直電場ドリフトが最低次に大きい順序では、平均加速度に $\vb v_E$ の微分も入る。この式の順序と混ぜない。高次の位相空間 Hamilton 形式では座標変換と位相空間体積を同時に求める。ここでは一次の力の平均に絞るため、その変分導出は :cite[bellan-plasma] を参照する。

:::example[平均の消える項]
$\vb w=w(\cos\theta,\sin\theta,0)$ では $\langle w_xw_y\rangle=0$ である。これに対し $\varrho_y=w_x/\Omega$ なので $\langle w_x\varrho_y\rangle=w^2/(2\Omega)$ は残る。
:::

:::example[勾配が急な場所]
$\rho=2\,\mathrm{mm}$ で $L=1\,\mathrm m$ なら $\rho/L=0.002$、$L=1\,\mathrm{mm}$ なら比は $2$ である。後者では Taylor 展開の二次項を小さいとできず、元の軌道方程式を使う。
:::

### 最低次のエネルギー保存の確認

静場で $E=0$ とする。式 :eqref[eq-plasma-gc-first] に従う横ドリフトの仕事はない。$\mathcal H=mv_\parallel^2/2+\mu B$ を微分すると、平行加速度による $-\mu v_\parallel\partial_\parallel B$ と、平行移流による $+\mu v_\parallel\partial_\parallel B$ が消える。勾配ドリフトは $\nabla B$ と直交して $B$ の変化に寄与しない。曲率ドリフトによる残差は、もともとのエネルギー変化尺度 $\mu B\abs\Omega$ に対し二次である。二次でも保存させるには、二次の $\dot v_\parallel$ と座標の補正を同じ近似から計算する。

:::example[近似の次数の比較]
$\rho/L=0.01$ の静場では、一次ドリフトの代表速度は旋回速度の百分の一である。二次の速度補正は同じ代表速度に対し $10^{-4}$ 程度となる。一次の係数が対称性で消える場所では、この比だけで実際の横速度の相対誤差を見積もれない。
:::

### 直線場での有限半径平均の誤差

一様磁場で、一つの横 Fourier 成分 $F(x)=F_0\cos(kx)$ を案内中心で評価する。零次旋回 $x=X+\rho\cos\vartheta$ を用いると正弦の寄与は $\vartheta\mapsto\vartheta+\pi$ で消え、$\langle F\rangle=F_0\cos(kX)J_0(k\rho)$ となる。$J_0(a)=\langle\cos(a\cos\vartheta)\rangle$ を Taylor 展開し、$\langle\cos^2\vartheta\rangle=1/2$、$\langle\cos^4\vartheta\rangle=3/8$ を使うと $1-a^2/4+a^4/64$ を得る。余弦の Taylor 誤差の点ごとの上界を積分すると、四次までの誤差は $a^6/2304$ 以下である。これは軌道に沿った平均の誤差であり、強い力が軌道そのものを変形させる場合の誤差は別にある。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $k$ | 横波数 | $\mathrm{m^{-1}}$ |
| $\rho$ | 未摂動の旋回半径 | $\mathrm{m}$ |
| $X$ | 案内中心 | $\mathrm{m}$ |

横位相の変動幅を $a$ とする。案内中心から距離 $\rho$ の点で最大になる。

$$
a=\rho k
$$

正確な平均因子を $J_0(a)$ とする。円周積分で決まる。

$$
J_0(a)=\frac{1}{2\pi}\int_0^{2\pi}\cos(a\cos\vartheta)\,\dd\vartheta
$$

二次近似の因子を $J_2^{\rm app}$ とする。$\cos^2\vartheta$ の平均を使う。

$$
J_2^{\rm app}=1 - \frac{a^{2}}{4}
$$

四次近似の因子を $J_4^{\rm app}$ とする。四乗角平均を加える。

$$
J_4^{\rm app}=J_2^{\rm app} + \frac{a^{4}}{64}
$$

四次近似の絶対誤差上界を $\delta_4$ とする。六乗角平均 $5/16$ を用いる。

$$
\delta_4=\frac{a^{6}}{2304}
$$

実際に観測する平均を $\langle F\rangle/F_0$ とする。案内中心の位相が節なら平均はゼロである。

$$
\langle F\rangle/F_0=J_0(a) \cos{\left(X k \right)}
$$

:::example[急な変化のない場]
二次で十分な場合の平均因子を計算する。

入力値を次のように置く。

$$
\begin{aligned}
k&=1\,\mathrm{m^{-1}} \\\rho&=\frac{1}{10}\,\mathrm{m} \\X&=0\,\mathrm{m}
\end{aligned}
$$

これらを順に代入する。

横位相の変動幅は

$$
a=\left(1\right) \left(\frac{1}{10}\right)=\frac{1}{10}
$$

である。

正確な平均因子は

$$
J_0(a)=J_{0}\left(\left(\frac{1}{10}\right)\right)\simeq 0.9975
$$

である。

二次近似の因子は

$$
J_2^{\rm app}=1 - \frac{\left(\frac{1}{10}\right)^{2}}{4}=\frac{399}{400}
$$

である。

四次近似の因子は

$$
J_4^{\rm app}=\frac{\left(\frac{1}{10}\right)^{4}}{64} + \left(\frac{399}{400}\right)=\frac{638401}{640000}
$$

である。

四次近似の絶対誤差上界は

$$
\delta_4=\frac{\left(\frac{1}{10}\right)^{6}}{2304}=\frac{1}{2304000000}
$$

である。

実際に観測する平均は

$$
\langle F\rangle/F_0=\left(J_{0}\left(\frac{1}{10}\right)\right) \cos{\left(\left(0\right) \left(1\right) \right)}\simeq 0.9975
$$

である。

点での評価には約 $0.25$ パーセントの差があるが、二次項でそのほとんどを補正できる。
:::

:::example[旋回半径と波数の逆数が等しい場]
平均を行わず案内中心で評価した値と比較する。

入力値を次のように置く。

$$
\begin{aligned}
k&=2\,\mathrm{m^{-1}} \\\rho&=\frac{1}{2}\,\mathrm{m} \\X&=0\,\mathrm{m}
\end{aligned}
$$

これらを順に代入する。

横位相の変動幅は

$$
a=\left(2\right) \left(\frac{1}{2}\right)=1
$$

である。

正確な平均因子は

$$
J_0(a)=J_{0}\left(\left(1\right)\right)\simeq 0.7652
$$

である。

二次近似の因子は

$$
J_2^{\rm app}=1 - \frac{\left(1\right)^{2}}{4}=\frac{3}{4}
$$

である。

四次近似の因子は

$$
J_4^{\rm app}=\frac{\left(1\right)^{4}}{64} + \left(\frac{3}{4}\right)=\frac{49}{64}
$$

である。

四次近似の絶対誤差上界は

$$
\delta_4=\frac{\left(1\right)^{6}}{2304}=\frac{1}{2304}
$$

である。

実際に観測する平均は

$$
\langle F\rangle/F_0=\left(J_{0}\left(1\right)\right) \cos{\left(\left(0\right) \left(2\right) \right)}\simeq 0.7652
$$

である。

平均値が約四分の一弱だけ減る。低周波でも空間平均の補正が必要である。
:::

:::example[案内中心が節にある場]
円周平均が消えるかどうか、点での評価と比較する。

入力値を次のように置く。

$$
\begin{aligned}
k&=2\,\mathrm{m^{-1}} \\\rho&=\frac{1}{4}\,\mathrm{m} \\X&=\frac{\pi}{4}\,\mathrm{m}
\end{aligned}
$$

これらを順に代入する。

横位相の変動幅は

$$
a=\left(2\right) \left(\frac{1}{4}\right)=\frac{1}{2}
$$

である。

正確な平均因子は

$$
J_0(a)=J_{0}\left(\left(\frac{1}{2}\right)\right)\simeq 0.93847
$$

である。

二次近似の因子は

$$
J_2^{\rm app}=1 - \frac{\left(\frac{1}{2}\right)^{2}}{4}=\frac{15}{16}
$$

である。

四次近似の因子は

$$
J_4^{\rm app}=\left(\frac{15}{16}\right) + \frac{\left(\frac{1}{2}\right)^{4}}{64}=\frac{961}{1024}
$$

である。

四次近似の絶対誤差上界は

$$
\delta_4=\frac{\left(\frac{1}{2}\right)^{6}}{2304}=\frac{1}{147456}
$$

である。

実際に観測する平均は

$$
\langle F\rangle/F_0=\left(J_{0}\left(\frac{1}{2}\right)\right) \cos{\left(\left(2\right) \left(\frac{\pi}{4}\right) \right)}=0
$$

である。

円の正負の寄与が打ち消し合うため、$J_0$ の値にかかわらず平均はゼロである。
:::

:::example[平均の符号が変わる場]
横位相幅が大きい場合の平均を計算する。

入力値を次のように置く。

$$
\begin{aligned}
k&=3\,\mathrm{m^{-1}} \\\rho&=1\,\mathrm{m} \\X&=0\,\mathrm{m}
\end{aligned}
$$

これらを順に代入する。

横位相の変動幅は

$$
a=\left(1\right) \left(3\right)=3
$$

である。

正確な平均因子は

$$
J_0(a)=J_{0}\left(\left(3\right)\right)\simeq -0.26005
$$

である。

二次近似の因子は

$$
J_2^{\rm app}=1 - \frac{\left(3\right)^{2}}{4}=- \frac{5}{4}
$$

である。

四次近似の因子は

$$
J_4^{\rm app}=\left(- \frac{5}{4}\right) + \frac{\left(3\right)^{4}}{64}=\frac{1}{64}
$$

である。

四次近似の絶対誤差上界は

$$
\delta_4=\frac{\left(3\right)^{6}}{2304}=\frac{81}{256}
$$

である。

実際に観測する平均は

$$
\langle F\rangle/F_0=\left(J_{0}\left(3\right)\right) \cos{\left(\left(0\right) \left(3\right) \right)}\simeq -0.26005
$$

である。

正確な平均因子は負である。この場合は四次展開の誤差上界が大きく、近似の符号を根拠に判定できない。
:::


## 演習問題

::::exercise{level=basic}
$\Omega=5\,\mathrm{s^{-1}}$、$v_\perp=10\,\mathrm{m\,s^{-1}}$ の半径を求めよ。

:::solution
速度を角周波数で割ると $\rho=10/5=2\,\mathrm m$ である。
:::
::::

::::exercise{level=basic}
$\vb E=3\hat{\vb y}\,\mathrm{V\,m^{-1}}$、$\vb B=2\hat{\vb z}\,\mathrm T$ の電場ドリフトを求めよ。

:::solution
$\hat{\vb y}\times\hat{\vb z}=\hat{\vb x}$ より $\vb v_E=(3/2)\hat{\vb x}\,\mathrm{m\,s^{-1}}$ である。
:::
::::

::::exercise{level=basic}
$\Omega<0$、$v_y(0)=0$、$v_x(0)>0$ なら初期の $\dot v_y$ の符号は何か。

:::solution
$\dot v_y=-\Omega v_x$ だから正である。
:::
::::

::::exercise{level=basic}
$\sin\alpha_0=1/3$ の粒子の反射磁場を求めよ。

:::solution
平行運動エネルギーは反射でゼロであり、式 :eqref[eq-plasma-mirror] から $B_m=9B_0$ である。
:::
::::

::::exercise{level=advanced}
$E_x=at$、$\vb B=B\hat{\vb z}$ に対する横速度の特解を求めよ。

:::solution
$v_y=-at/B$、$v_x=a/(B\Omega)$ と置く。$\dot v_x=0=q at/m+\Omega v_y$、$\dot v_y=-a/B=-\Omega v_x$ なので厳密な特解である。初期条件に合わせるには同次円運動を足す。
:::
::::

::::exercise{level=advanced}
$B(s)=B_0+bs^2/2$、$b>0$ の軸近くの往復周波数を求めよ。

:::solution
鏡力は $m\ddot s=-\mu b s$ である。従って角周波数は $\sqrt{\mu b/m}$ である。軌道振幅が磁場の二次近似の範囲に入る場合に成り立つ。
:::
::::

::::exercise{level=advanced}
$B_z=B_0+ax>0$ の領域で $q>0$、$a>0$、$v_\perp>0$ の粒子の勾配ドリフトの向きを求めよ。

:::solution
式 :eqref[eq-plasma-gradb] で $\vb b\times\nabla B=a\hat{\vb y}$ となる。$\mu/(qB)>0$ だから $+\hat{\vb y}$ である。
:::
::::

::::exercise{level=advanced}
静電場のミラーでは反射条件がどう変わるか。

:::solution
全エネルギーは $W=mv_\parallel^2/2+\mu B+q\phi$ である。反射点では $\mu B(s_m)+q\phi(s_m)=W$ となる。その点に到達する手前で $W-\mu B-q\phi\ge0$ となる枝を取る。
:::
::::

## 確認テスト

::::quiz
磁場だけがあるとき保存される量を選べ。

- [x] 全運動エネルギー
- [ ] 任意の場における平行速度
- [ ] 速度の全成分

:::explanation
Lorentz 力は速度に垂直なので仕事がゼロである。磁気ミラーでは平行速度と垂直速度は別々には保存されない。
:::
::::

::::quiz{answer="5" tol="1e-6"}
$\sin^2\alpha_0=1/5$ のとき $B_m/B_0$ を答えよ。

:::explanation
$W=\mu B_m$、$\mu B_0=W/5$ なので $B_m/B_0=5$ である。
:::
::::

::::quiz
正負の種に共通のドリフトを選べ。

- [x] 一定の交差電磁場での電場ドリフト
- [ ] 重力ドリフト
- [ ] 勾配ドリフト

:::explanation
$\vb F=q\vb E$ とおけば $\vb F\times\vb B/(qB^2)$ の $q$ が消える。他の二つでは残る。
:::
::::

::::quiz{answer="3" tol="1e-6"}
$v_\perp=30\,\mathrm{m\,s^{-1}}$、$\abs\Omega=10\,\mathrm{s^{-1}}$ の旋回半径を m 単位で答えよ。

:::explanation
$\rho=v_\perp/\abs\Omega=30/10=3\,\mathrm m$ である。
:::
::::

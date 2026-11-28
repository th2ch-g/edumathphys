---
title: 宇宙・天体プラズマ
description: 太陽風、衝撃波、粒子加速、降着円盤を、保存則と局所安定性から調べる。
notation:
  - symbol: 'r_c'
    meaning: Parker 風の臨界半径
  - symbol: 'r'
    meaning: 衝撃波の圧縮比（後半）
  - symbol: '\kappa'
    meaning: 円盤のエピサイクル角周波数
  - symbol: '\Omega'
    meaning: 円盤の回転角周波数
---

## この章の目標

- Parker 風の遷音速解を導く。
- 磁気圏とフレアのエネルギー収支を説明する。
- Fermi 加速のスペクトルと磁気回転不安定性の条件を導く。

## 太陽風と Parker 解

太陽の外層では重力が物質を引く一方、外へ下がる圧力が流れを加速する。:term[太陽風]{en="solar wind" yomi="たいようふう"} の最小模型として、質量 $M$ の中心天体、一様温度で $p=c_s^2\rho$、定常球対称、磁気力と回転を無視した流れを考える。質量流量 $\dot M=4\pi r^2\rho u$ は一定である。運動式 $u u'=-c_s^2(\ln\rho)'-GM/r^2$ に $(\ln\rho)'=-2/r-u'/u$ を入れると
$$
\left(u-\frac{c_s^2}{u}\right)u'=\frac{2c_s^2}{r}-\frac{GM}{r^2}.
\label{eq-plasma-parker}
$$
ここで $G$ は重力定数である。音速 $u=c_s$ を有限の微分で通るには、右辺が同時にゼロとなる $r_c=GM/(2c_s^2)$ でなければならない。両辺をさらに微分すると臨界点で $2(u_c')^2=2c_s^2/r_c^2$ となり、加速解では $u_c'=c_s/r_c$ である。

### 遷音速の枝

$m=u/c_s,\ x=r/r_c$ として式を積分し、$m=x=1$ を課すと
$$
m^2-\ln m^2=4\ln x+4/x-3 .
\label{eq-plasma-parker-integral}
$$
左辺は $m>0$ で $m=1$ が唯一の最小で、値は $1$ である。右辺も $x=1$ が唯一の最小で同じ値である。各 $x\ne1$ で亜音速と超音速の値が一つずつある。内側で亜音速、外側で超音速を選び、上の正の傾きで接続したものが :term[Parker 解]{en="Parker wind solution" yomi="ぱーかーかい"} である。この境界選択は有限の外圧を課した亜音速流と異なる。等温性はエネルギー式の仮定であり、熱伝導、非熱電子、磁場や波を含む太陽風では加熱の機構も別に追う。

:::example[臨界半径の温度依存]
単一イオンと電子が等温なら $c_s^2=(T_e+T_i)/m_i$ である。両温度が二倍になると $c_s^2$ は二倍となり、$r_c$ は半分になる。
:::

:::example[遠方の風]
$x\gg1$ では $4/x$ が小さく、超音速根では $\ln m^2$ が $m^2$ より遅く増える。従って $m\simeq2\sqrt{\ln x}$ が先頭近似である。等温を保つための加熱を仮定した結果なので、無限遠の実際の速度に外挿するにはエネルギー式を確認する。
:::

### 等温でない風のエネルギー

圧力を $p=K\rho^{\Gamma}$、$\Gamma>1$ とし、$K$ が流れに沿って一定のポリトロープを考える。局所音速は $c^2=\dd p/\dd\rho=\Gamma K\rho^{\Gamma-1}$ であり、運動方程式に質量保存を入れると、等温の場合と同様に $(u-c^2/u)u'=2c^2/r-GM/r^2$ となる。ただし $c$ は密度に依存する。運動式を直接積分すれば、比エンタルピー $\int\dd p/\rho=c^2/(\Gamma-1)$ を用いて
$$
\frac{u^2}{2}+\frac{c^2}{\Gamma-1}-\frac{GM}{r}=\mathcal B
$$
となる。無限遠で密度がゼロへ下がり、かつ速度が正の有限値へ近づく枝では $u_\infty=\sqrt{2\mathcal B}$ である。等温風の遠方速度が対数的に増えるのは温度条件の違いによる。滑らかな臨界点で $GM/r_c=2c_c^2$ を使うと $\mathcal B=c_c^2[1/(\Gamma-1)-3/2]$ となる。従ってこの仮定のもとで $\Gamma=5/3$ の解は、同じ臨界点を通って無限遠で正の有限速度を保つことができない。エネルギーの必要条件から得た結論であり、$\mathcal B>0$ の場合の枝選択には内側の境界条件も要る。

### 凍結場の渦巻き

太陽の足元が角速度 $\Omega_\odot$ で回転し、赤道上で風が一定の速さ $u$ で放射状に流れるとする。半径 $r$ までの移動時間は $(r-r_0)/u$ なので、同じ時刻に見た磁力線の足元との方位角差は $-\Omega_\odot(r-r_0)/u$ となる。軌跡の接線をとれば $B_\varphi/B_r=r\dd\varphi/\dd r=-\Omega_\odot r/u$ である。球面で放射状磁束を保存すれば $B_r\propto r^{-2}$、遠方で $B_\varphi\propto r^{-1}$ となる。この運動学的模型は風の加速域と緯度依存を省いている。

:::example[成分が等しくなる半径]
この描像では $r=u/\Omega_\odot$ で $\abs{B_\varphi}=\abs{B_r}$ になる。風の速度を二倍とすれば、同じ回転速度でこの半径も二倍になり、同じ半径でのねじれは半分になる。
:::

### 亜音速風と外側の境界

式 :eqref[eq-plasma-parker-integral] の変数を $y=m^2$、右辺を $G(x)$ とすると $y e^{-y}=e^{-G(x)}$ である。$W(z)e^{W(z)}=z$ の逆関数で書けば、$y=-W_0(-e^{-G})$ が亜音速、$y=-W_{-1}(-e^{-G})$ が超音速である。二枝は $G=1$ で $y=1$ に合流する。積分定数を臨界値からずらし $G=4\ln x+4/x+C$ とし、$C>-3$ なら右辺の最小値は $1$ より大きい。このとき全領域で亜音速枝を選べる。遠方で $y\ll1$ なら $y\simeq e^{-G}\propto x^{-4}$、$u\propto r^{-2}$ である。質量保存により密度が有限値へ近づくので、そのような亜音速の流出には有限の外圧が伴う。

:::example[臨界解と遠方境界]
遷音速解の遠方は $u$ が $\sqrt{\ln r}$ 程度で増し、$\rho=\dot M/(4\pi r^2 u)$ はゼロへ近づく。亜音速枝は $r^2u$ が有限の定数に近づき、密度はゼロにならない。同じ内部の等温方程式でも、外側の圧力条件が解の枝を選ぶ。
:::

### Parker 解を二つの実数枝から復元する

$G=4\ln x+4/x-3$ とし $x=r/r_c$ を指定する。$y=m^2$ について $ye^{-y}=e^{-G}$ は Lambert の $W$ 関数で逆にできる。外向きの遷音速解では $x<1$ で $-W_0(-e^{-G})$、$x>1$ で $-W_{-1}(-e^{-G})$ を使う。局所勾配は :eqref[eq-plasma-parker] から $\dd m/\dd x=2(x-1)/[x^2(m-1/m)]$ である。臨界点では分子と分母が同時に消えるので、正の枝ではその点だけ $\dd m/\dd x=1$ を使う。臨界点での密度を $\rho_c$ とすると $\rho/\rho_c=1/(x^2m)$ である。上流温度が同じでも外側境界が違えば枝が変わる。以下では数値の分枝を明示する。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $x$ | 規格化距離 $r/r_c$。臨界点そのものを除く。 | 1 |
| $b$ | Lambert の実枝番号 | 1 |

Parker の積分定数を含む関数を $G(x)$ とする。最小値は臨界点で一である。

$$
G(x)=4 \ln{\left(x \right)} - 3 + \frac{4}{x}
$$

二乗 Mach 数を $m^2$ とする。枝を選んで実数の根を取る。

$$
m^2=- W_{b}\left(- e^{- G(x)}\right)
$$

Mach 数を $m$ とする。外向きの正根である。

$$
m=\sqrt{m^{2}}
$$

速度の局所勾配を $\dd m/\dd x$ とする。臨界点の有理式の極限は直接代入しない。

$$
\dd m/\dd x=\frac{2 x - 2}{x^{2} \left(m - \frac{1}{m}\right)}
$$

密度比を $\rho/\rho_c$ とする。球面の面積と速度を掛けて質量を保つ。

$$
\rho/\rho_c=\frac{1}{m x^{2}}
$$

動圧の規格化値を $\rho u^2/(\rho_c c_s^2)$ とする。密度の変化と速度の変化を共に入れる。

$$
\rho u^2/(\rho_c c_s^2)=\rho/\rho_{c} m^{2}
$$

:::example[臨界点の内側]
亜音速側の速度と密度を求める。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{2} \\b&=0
\end{aligned}
$$

これらを順に代入する。

Parker の積分定数を含む関数は

$$
G(x)=4 \ln{\left(\left(\frac{1}{2}\right) \right)} - 3 + \frac{4}{\left(\frac{1}{2}\right)}\simeq 2.2274
$$

である。

二乗 Mach 数は

$$
m^2=- W_{\left(0\right)}\left(- e^{- \left(5 - 4 \ln{\left(2 \right)}\right)}\right)\simeq 0.12177
$$

である。

Mach 数は

$$
m=\sqrt{\left(- W\left(- \frac{16}{e^{5}}\right)\right)}\simeq 0.34895
$$

である。

速度の局所勾配は

$$
\dd m/\dd x=\frac{2 \left(\frac{1}{2}\right) - 2}{\left(\frac{1}{2}\right)^{2} \left(\left(\sqrt{- W\left(- \frac{16}{e^{5}}\right)}\right) - \frac{1}{\left(\sqrt{- W\left(- \frac{16}{e^{5}}\right)}\right)}\right)}\simeq 1.5893
$$

である。

密度比は

$$
\rho/\rho_c=\frac{1}{\left(\frac{1}{2}\right)^{2} \left(\sqrt{- W\left(- \frac{16}{e^{5}}\right)}\right)}\simeq 11.463
$$

である。

動圧の規格化値は

$$
\rho u^2/(\rho_c c_s^2)=\left(- W\left(- \frac{16}{e^{5}}\right)\right) \left(\frac{4}{\sqrt{- W\left(- \frac{16}{e^{5}}\right)}}\right)\simeq 1.3958
$$

である。

内側の密度は大きく、風速は音速を下回る。
:::

:::example[臨界点の外側]
同じ解の超音速側を求める。

入力値を次のように置く。

$$
\begin{aligned}
x&=2 \\b&=-1
\end{aligned}
$$

これらを順に代入する。

Parker の積分定数を含む関数は

$$
G(x)=4 \ln{\left(\left(2\right) \right)} - 3 + \frac{4}{\left(2\right)}\simeq 1.7726
$$

である。

二乗 Mach 数は

$$
m^2=- W_{\left(-1\right)}\left(- e^{- \left(-1 + 4 \ln{\left(2 \right)}\right)}\right)\simeq 2.8034
$$

である。

Mach 数は

$$
m=\sqrt{\left(- W_{-1}\left(- \frac{e}{16}\right)\right)}\simeq 1.6743
$$

である。

速度の局所勾配は

$$
\dd m/\dd x=\frac{2 \left(2\right) - 2}{\left(2\right)^{2} \left(\left(\sqrt{- W_{-1}\left(- \frac{e}{16}\right)}\right) - \frac{1}{\left(\sqrt{- W_{-1}\left(- \frac{e}{16}\right)}\right)}\right)}\simeq 0.46421
$$

である。

密度比は

$$
\rho/\rho_c=\frac{1}{\left(2\right)^{2} \left(\sqrt{- W_{-1}\left(- \frac{e}{16}\right)}\right)}\simeq 0.14931
$$

である。

動圧の規格化値は

$$
\rho u^2/(\rho_c c_s^2)=\left(- W_{-1}\left(- \frac{e}{16}\right)\right) \left(\frac{1}{4 \sqrt{- W_{-1}\left(- \frac{e}{16}\right)}}\right)\simeq 0.41859
$$

である。

速度は一より大きいが密度は球面の広がりで低い。
:::

:::example[外側で亜音速枝を使う比較]
外側で亜音速の根を選ぶと勾配がどうなるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=2 \\b&=0
\end{aligned}
$$

これらを順に代入する。

Parker の積分定数を含む関数は

$$
G(x)=4 \ln{\left(\left(2\right) \right)} - 3 + \frac{4}{\left(2\right)}\simeq 1.7726
$$

である。

二乗 Mach 数は

$$
m^2=- W_{\left(0\right)}\left(- e^{- \left(-1 + 4 \ln{\left(2 \right)}\right)}\right)\simeq 0.20948
$$

である。

Mach 数は

$$
m=\sqrt{\left(- W\left(- \frac{e}{16}\right)\right)}\simeq 0.4577
$$

である。

速度の局所勾配は

$$
\dd m/\dd x=\frac{2 \left(2\right) - 2}{\left(2\right)^{2} \left(\left(\sqrt{- W\left(- \frac{e}{16}\right)}\right) - \frac{1}{\left(\sqrt{- W\left(- \frac{e}{16}\right)}\right)}\right)}\simeq -0.28949
$$

である。

密度比は

$$
\rho/\rho_c=\frac{1}{\left(2\right)^{2} \left(\sqrt{- W\left(- \frac{e}{16}\right)}\right)}\simeq 0.54622
$$

である。

動圧の規格化値は

$$
\rho u^2/(\rho_c c_s^2)=\left(- W\left(- \frac{e}{16}\right)\right) \left(\frac{1}{4 \sqrt{- W\left(- \frac{e}{16}\right)}}\right)\simeq 0.11442
$$

である。

負の勾配になり、内側から加速してつながる Parker 風とは異なる。
:::

:::example[十分遠い遷音速風]
加速と動圧が同じ方向に変わるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=20 \\b&=-1
\end{aligned}
$$

これらを順に代入する。

Parker の積分定数を含む関数は

$$
G(x)=4 \ln{\left(\left(20\right) \right)} - 3 + \frac{4}{\left(20\right)}\simeq 9.1829
$$

である。

二乗 Mach 数は

$$
m^2=- W_{\left(-1\right)}\left(- e^{- \left(- \frac{14}{5} + 4 \ln{\left(20 \right)}\right)}\right)\simeq 11.637
$$

である。

Mach 数は

$$
m=\sqrt{\left(- W_{-1}\left(- \frac{e^{\frac{14}{5}}}{160000}\right)\right)}\simeq 3.4113
$$

である。

速度の局所勾配は

$$
\dd m/\dd x=\frac{2 \left(20\right) - 2}{\left(20\right)^{2} \left(\left(\sqrt{- W_{-1}\left(- \frac{e^{\frac{14}{5}}}{160000}\right)}\right) - \frac{1}{\left(\sqrt{- W_{-1}\left(- \frac{e^{\frac{14}{5}}}{160000}\right)}\right)}\right)}\simeq 0.030466
$$

である。

密度比は

$$
\rho/\rho_c=\frac{1}{\left(20\right)^{2} \left(\sqrt{- W_{-1}\left(- \frac{e^{\frac{14}{5}}}{160000}\right)}\right)}\simeq 0.00073285
$$

である。

動圧の規格化値は

$$
\rho u^2/(\rho_c c_s^2)=\left(- W_{-1}\left(- \frac{e^{\frac{14}{5}}}{160000}\right)\right) \left(\frac{1}{400 \sqrt{- W_{-1}\left(- \frac{e^{\frac{14}{5}}}{160000}\right)}}\right)\simeq 0.0085283
$$

である。

風速は増しても動圧は半径の増加とともに小さくなる。
:::


### Parker 渦巻きの磁場とエネルギー密度

一定速さの赤道風で $r_s=u/\Omega_\odot$ を定義する。放射成分 $B_r=B_s(r_s/r)^2$、方位成分 $B_\varphi=-(r/r_s)B_r$ である。質量保存で密度は $\rho_s(r_s/r)^2$ なので、$v_A^2=B^2/(\mu_0\rho)$ は $v_{As,r}^2[1+(r_s/r)^2]$ となる。遠方では方位成分が主となり、この単純模型の $v_A$ は一定値に近づく。方位成分を除いた仮想の放射磁場だけの計算では $v_A\propto1/r$ となり、大きい半径で差が出る。回転軸に近い場では赤道における方位速度をそのまま使えず、ここでの式を角度によらず適用しない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $x$ | ねじれ半径に対する距離 $r/r_s$ | 1 |
| $M_s$ | 風速を $v_{As,r}$ で割った値 | 1 |

放射磁場の比を $B_r/B_s$ とする。球面の磁束保存に従う。

$$
B_r/B_s=\frac{1}{x^{2}}
$$

方位磁場の比を $B_\varphi/B_s$ とする。足元の回転がねじれを作る。

$$
B_\varphi/B_s=- \frac{1}{x}
$$

密度比を $\rho/\rho_s$ とする。同じ一定速度で質量流が球面を通る。

$$
\rho/\rho_s=(B_r/B_s)
$$

規格化磁気圧を $B^2/B_s^2$ とする。正負の成分を二乗して足す。

$$
B^2/B_s^2=(B_\varphi/B_s)^{2} + (B_r/B_s)^{2}
$$

Alfvén 速度の比を $v_A/v_{As,r}$ とする。磁場強度と密度の両方が変わる。

$$
v_A/v_{As,r}=\sqrt{\frac{(B^2/B_s^2)}{\rho/\rho_{s}}}
$$

流れの Alfvén Mach 数を $M_A$ とする。速波の Mach 数とは区別する。

$$
M_A=\frac{M_{s}}{v_A/v_{As,r}}
$$

:::example[放射方向が主な位置]
内側で一定風速の近似を形式的に評価する。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{10} \\M_s&=4
\end{aligned}
$$

これらを順に代入する。

放射磁場の比は

$$
B_r/B_s=\frac{1}{\left(\frac{1}{10}\right)^{2}}=100
$$

である。

方位磁場の比は

$$
B_\varphi/B_s=- \frac{1}{\left(\frac{1}{10}\right)}=-10
$$

である。

密度比は

$$
\rho/\rho_s=\left(100\right)=100
$$

である。

規格化磁気圧は

$$
B^2/B_s^2=\left(-10\right)^{2} + \left(100\right)^{2}=10100
$$

である。

Alfvén 速度の比は

$$
v_A/v_{As,r}=\sqrt{\frac{\left(10100\right)}{\left(100\right)}}\simeq 10.05
$$

である。

流れの Alfvén Mach 数は

$$
M_A=\frac{\left(4\right)}{\left(\sqrt{101}\right)}\simeq 0.39801
$$

である。

放射成分が支配する。この入力で風が亜 Alfvén 的なら磁気力を無視した一定速さは別途確認する。
:::

:::example[二成分が等しい位置]
回転の影響が同じ大きさになる位置を求める。

入力値を次のように置く。

$$
\begin{aligned}
x&=1 \\M_s&=4
\end{aligned}
$$

これらを順に代入する。

放射磁場の比は

$$
B_r/B_s=\frac{1}{\left(1\right)^{2}}=1
$$

である。

方位磁場の比は

$$
B_\varphi/B_s=- \frac{1}{\left(1\right)}=-1
$$

である。

密度比は

$$
\rho/\rho_s=\left(1\right)=1
$$

である。

規格化磁気圧は

$$
B^2/B_s^2=\left(-1\right)^{2} + \left(1\right)^{2}=2
$$

である。

Alfvén 速度の比は

$$
v_A/v_{As,r}=\sqrt{\frac{\left(2\right)}{\left(1\right)}}\simeq 1.4142
$$

である。

流れの Alfvén Mach 数は

$$
M_A=\frac{\left(4\right)}{\left(\sqrt{2}\right)}\simeq 2.8284
$$

である。

二乗磁場には等しい二つの寄与が加わる。
:::

:::example[方位方向が主な位置]
遠方の Alfvén 速度を調べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=10 \\M_s&=4
\end{aligned}
$$

これらを順に代入する。

放射磁場の比は

$$
B_r/B_s=\frac{1}{\left(10\right)^{2}}=\frac{1}{100}
$$

である。

方位磁場の比は

$$
B_\varphi/B_s=- \frac{1}{\left(10\right)}=- \frac{1}{10}
$$

である。

密度比は

$$
\rho/\rho_s=\left(\frac{1}{100}\right)=\frac{1}{100}
$$

である。

規格化磁気圧は

$$
B^2/B_s^2=\left(- \frac{1}{10}\right)^{2} + \left(\frac{1}{100}\right)^{2}=\frac{101}{10000}
$$

である。

Alfvén 速度の比は

$$
v_A/v_{As,r}=\sqrt{\frac{\left(\frac{101}{10000}\right)}{\left(\frac{1}{100}\right)}}\simeq 1.005
$$

である。

流れの Alfvén Mach 数は

$$
M_A=\frac{\left(4\right)}{\left(\frac{\sqrt{101}}{10}\right)}\simeq 3.9801
$$

である。

磁場の大きさと密度の平方根が同じ半径依存になる。
:::

:::example[同じ規格化位置の高速風]
風速で規格化した効果を計算する。

入力値を次のように置く。

$$
\begin{aligned}
x&=1 \\M_s&=8
\end{aligned}
$$

これらを順に代入する。

放射磁場の比は

$$
B_r/B_s=\frac{1}{\left(1\right)^{2}}=1
$$

である。

方位磁場の比は

$$
B_\varphi/B_s=- \frac{1}{\left(1\right)}=-1
$$

である。

密度比は

$$
\rho/\rho_s=\left(1\right)=1
$$

である。

規格化磁気圧は

$$
B^2/B_s^2=\left(-1\right)^{2} + \left(1\right)^{2}=2
$$

である。

Alfvén 速度の比は

$$
v_A/v_{As,r}=\sqrt{\frac{\left(2\right)}{\left(1\right)}}\simeq 1.4142
$$

である。

流れの Alfvén Mach 数は

$$
M_A=\frac{\left(8\right)}{\left(\sqrt{2}\right)}\simeq 5.6569
$$

である。

同じ $r/r_s$ を比較しているので、物理的な距離は速度に応じて変わる。
:::


## 地球磁気圏とバウショック

太陽風の法線速度が速い磁気音波の速度を上回ると、地球磁気圏の前方では小さな摂動が上流へ伝わらない。そこで流れの前方に :term[バウショック]{en="bow shock" yomi="ばうしょっく"} が形成される。衝撃波を通った流れは磁気圏の周囲の磁気シースへ入り、磁気圏界面で内外の法線応力がつりあう。この応力には熱圧と磁気応力も入る。したがって磁気シースはバウショックと磁気圏界面の間にある。

### 磁気圏界面の尺度

昼側の法線応力を $\rho_{\rm sw}U_{\rm sw}^2$、双極子の内側磁場を $B_m=B_*(R_E/r)^3$ と近似する。内部熱圧、太陽風の磁気圧、界面電流による補正を省くと
$$
\frac{r_{\rm mp}}{R_E}
\simeq\left(\frac{B_*^2}{2\mu_0\rho_{\rm sw}U_{\rm sw}^2}\right)^{1/6}.
$$
指数の $1/6$ は双極子の磁気圧が $r^{-6}$ に比例するためである。位置の精密計算では衝撃波通過後の全圧と、磁場の向きによる磁気再結合を考慮する。昼側で開いた磁力線は太陽風に引かれ、夜側に尾部を作る。この流れは全球の磁束と電流の収支に結びつく。昼側の再結合で開く磁束の単位時間当たりの量を $\mathcal V_{\rm d}$、尾部の再結合で閉じる量を $\mathcal V_{\rm n}$ とする。いずれも再結合電場の線積分で与えられる電圧の向きと大きさで表す。北半球側の開いた磁束を $\Phi_{\rm open}$ とすると、Faraday 則により $\dot\Phi_{\rm open}=\mathcal V_{\rm d}-\mathcal V_{\rm n}$ である。例えば $\mathcal V_{\rm d}=60\,\mathrm{kV}$、$\mathcal V_{\rm n}=20\,\mathrm{kV}$ が二分続けば、開いた磁束は $4.8\times10^6\,\mathrm{Wb}$ 増える。この積分収支は、電離圏との応答や境界位置の変化を計算する際の条件になる。

:::example[外圧の変化]
太陽風の動圧が $64$ 倍になれば、単純な圧力模型で $r_{\rm mp}$ は半分になる。バウショックまでの距離は圧力だけでは決まらず、衝撃波の Mach 数にも依存する。
:::

:::example[上流 Mach 数]
垂直な場で $c_s=30,\ v_A=40,\ U_n=100$ を同じ速度単位で与えると速波速度は $50$、法線 Mach 数は $2$ である。したがって上流は速波に対して超速波である。これは粒子反射が必須となる運動論的な臨界 Mach 数の判定とは別である。
:::

:::example[圧縮と磁気シース]
磁場が法線に垂直で、密度比 $r=3$ の MHD 衝撃波を通ると $B_{t2}=3B_{t1}$ である。磁場の圧力は九倍になるため、磁気圏界面で上流の磁気圧をそのまま用いると圧力の分担を誤る。このような差を簡単な動圧の尺度式では省いている。
:::

### 磁気圏界面の圧力収支の感度

外部が動圧 $P_d$ と、それ以外の熱圧・磁気応力をまとめた正の法線圧 $P_t$ でプラズマを押すとする。昼側の局所モデルで内部の熱圧を $P_b$ とし、$P_d+P_t-P_b>0$ なら磁気圧が不足分を支える。双極子の強度を $B_*$ とすると $x=r_{\rm mp}/R_E=[B_*^2/(2\mu_0(P_d+P_t-P_b))]^{1/6}$ である。外側の磁気応力は法線成分があると磁気圧だけにならないので、この式では指定済みの外圧を使う。$P_d$ のみを微小に変える感度は $\partial\ln x/\partial\ln P_d=-P_d/[6(P_d+P_t-P_b)]$ である。動圧の変化を風速に換算する場合は $P_d\propto U^2$ の因子二が加わる。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $D$ | $P_d/P_*$。$P_*=B_*^2/(2\mu_0)$ である。 | 1 |
| $T$ | $P_t/P_*$ | 1 |
| $B$ | $P_b/P_*$ | 1 |

内部磁気圧が支える規格化圧力を $C$ とする。外側全圧から内部熱圧を引く。

$$
C=- B + D + T
$$

界面距離の比を $r_{\rm mp}/R_E$ とする。磁気圧は半径の六乗の逆数である。

$$
r_{\rm mp}/R_E=\frac{1}{\sqrt[6]{C}}
$$

外側の圧力で動圧の担う割合を $F_d$ とする。界面の感度そのものとは異なる。

$$
F_d=\frac{D}{D + T}
$$

動圧だけを変えたときの感度を $\partial\ln x/\partial\ln P_d$ とする。他の二圧力を一定とする。

$$
\partial\ln x/\partial\ln P_d=- \frac{D}{6 C}
$$

風速を変えるときの感度を $\partial\ln x/\partial\ln U$ とする。密度が変わらない場合である。

$$
\partial\ln x/\partial\ln U=2 \partial\ln x/\partial\ln P_{d}
$$

界面磁場の値を $B_m/B_*$ とする。双極子式を圧力収支から逆に求める。

$$
B_m/B_*=\sqrt{C}
$$

:::example[動圧だけの基準]
双極子による界面の距離を計算する。

入力値を次のように置く。

$$
\begin{aligned}
D&=\frac{1}{1000000} \\T&=0 \\B&=0
\end{aligned}
$$

これらを順に代入する。

内部磁気圧が支える規格化圧力は

$$
C=\left(\frac{1}{1000000}\right)=\frac{1}{1000000}
$$

である。

界面距離の比は

$$
r_{\rm mp}/R_E=\frac{1}{\sqrt[6]{\left(\frac{1}{1000000}\right)}}=10
$$

である。

外側の圧力で動圧の担う割合は

$$
F_d=\frac{\left(\frac{1}{1000000}\right)}{\left(0\right) + \left(\frac{1}{1000000}\right)}=1
$$

である。

動圧だけを変えたときの感度は

$$
\partial\ln x/\partial\ln P_d=- \frac{1}{6}
$$

である。

風速を変えるときの感度は

$$
\partial\ln x/\partial\ln U=2 \left(- \frac{1}{6}\right)=- \frac{1}{3}
$$

である。

界面磁場の値は

$$
B_m/B_*=\sqrt{\left(\frac{1}{1000000}\right)}=\frac{1}{1000}
$$

である。

六乗根により圧力に対する距離の依存は緩やかである。
:::

:::example[外側の熱・磁気応力も同程度]
追加の外圧で位置と感度を求める。

入力値を次のように置く。

$$
\begin{aligned}
D&=\frac{1}{1000000} \\T&=\frac{1}{1000000} \\B&=0
\end{aligned}
$$

これらを順に代入する。

内部磁気圧が支える規格化圧力は

$$
C=- \left(0\right) + 2 \left(\frac{1}{1000000}\right)=\frac{1}{500000}
$$

である。

界面距離の比は

$$
r_{\rm mp}/R_E=\frac{1}{\sqrt[6]{\left(\frac{1}{500000}\right)}}\simeq 8.909
$$

である。

外側の圧力で動圧の担う割合は

$$
F_d=\frac{1}{2}
$$

である。

動圧だけを変えたときの感度は

$$
\partial\ln x/\partial\ln P_d=- \frac{\left(\frac{1}{1000000}\right)}{6 \left(\frac{1}{500000}\right)}=- \frac{1}{12}
$$

である。

風速を変えるときの感度は

$$
\partial\ln x/\partial\ln U=2 \left(- \frac{1}{12}\right)=- \frac{1}{6}
$$

である。

界面磁場の値は

$$
B_m/B_*=\sqrt{\left(\frac{1}{500000}\right)}\simeq 0.0014142
$$

である。

全圧は増えるが、風速だけの変動で変わる割合は下がる。
:::

:::example[内部熱圧の補償]
基準と同じ位置で感度も同じか調べる。

入力値を次のように置く。

$$
\begin{aligned}
D&=\frac{1}{1000000} \\T&=\frac{1}{1000000} \\B&=\frac{1}{1000000}
\end{aligned}
$$

これらを順に代入する。

内部磁気圧が支える規格化圧力は

$$
C=\left(\frac{1}{1000000}\right)=\frac{1}{1000000}
$$

である。

界面距離の比は

$$
r_{\rm mp}/R_E=\frac{1}{\sqrt[6]{\left(\frac{1}{1000000}\right)}}=10
$$

である。

外側の圧力で動圧の担う割合は

$$
F_d=\frac{1}{2}
$$

である。

動圧だけを変えたときの感度は

$$
\partial\ln x/\partial\ln P_d=- \frac{1}{6}
$$

である。

風速を変えるときの感度は

$$
\partial\ln x/\partial\ln U=2 \left(- \frac{1}{6}\right)=- \frac{1}{3}
$$

である。

界面磁場の値は

$$
B_m/B_*=\sqrt{\left(\frac{1}{1000000}\right)}=\frac{1}{1000}
$$

である。

この入力では動圧変更に対する感度も基準に戻る。ただし外圧と内圧の各寄与は異なる。
:::

:::example[内圧の補償が強い]
磁気圧に必要な不足分が小さいときの形式値を求める。

入力値を次のように置く。

$$
\begin{aligned}
D&=\frac{1}{1000000} \\T&=0 \\B&=\frac{9}{10000000}
\end{aligned}
$$

これらを順に代入する。

内部磁気圧が支える規格化圧力は

$$
C=\left(0\right) + \left(\frac{1}{1000000}\right) - \left(\frac{9}{10000000}\right)=\frac{1}{10000000}
$$

である。

界面距離の比は

$$
r_{\rm mp}/R_E=\frac{1}{\sqrt[6]{\left(\frac{1}{10000000}\right)}}\simeq 14.678
$$

である。

外側の圧力で動圧の担う割合は

$$
F_d=\frac{\left(\frac{1}{1000000}\right)}{\left(0\right) + \left(\frac{1}{1000000}\right)}=1
$$

である。

動圧だけを変えたときの感度は

$$
\partial\ln x/\partial\ln P_d=- \frac{\left(\frac{1}{1000000}\right)}{6 \left(\frac{1}{10000000}\right)}=- \frac{5}{3}
$$

である。

風速を変えるときの感度は

$$
\partial\ln x/\partial\ln U=2 \left(- \frac{5}{3}\right)=- \frac{10}{3}
$$

である。

界面磁場の値は

$$
B_m/B_*=\sqrt{\left(\frac{1}{10000000}\right)}\simeq 0.00031623
$$

である。

位置の感度が大きくなるため、内部圧を一定と置いたことの影響を点検する。
:::


### バウショックの上流条件を斜め磁場で判定する

上流での法線方向の速い磁気音波速度は、$C=c_s^2/v_A^2$、$S=\cos^2\theta_{Bn}$ として $v_f^2/v_A^2=[1+C+\sqrt{(1+C)^2-4CS}]/2$ である。これは等方の理想 MHD で扱う圧縮波の速い根である。法線流の速波 Mach 数 $M_f=U_n/v_f$ が一を越えると上流は法線方向に超速波となる。$M_f>1$ は速い圧縮衝撃波を議論する出発点で、観測した傾斜層の質量・エネルギー収支や運動論的臨界性をここからは決められない。磁場も流れも斜めなら、流速全体でなく衝撃面の法線成分を使う。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $C$ | 上流音速と Alfvén 速度の二乗比 | 1 |
| $S$ | 磁場と衝撃面法線の角の二乗余弦 | 1 |
| $U_n/v_A$ | 上流法線速度の比 | 1 |

速波と遅波を分ける判別式を $D$ とする。磁場が斜めなら角度によって変わる。

$$
D=- 4 C S + \left(C + 1\right)^{2}
$$

速波の二乗速度比を $v_f^2/v_A^2$ とする。大きい方の根である。

$$
v_f^2/v_A^2=\frac{C}{2} + \frac{\sqrt{D}}{2} + \frac{1}{2}
$$

遅波の二乗速度比を $v_s^2/v_A^2$ とする。二根の積から求める。

$$
v_s^2/v_A^2=\frac{C S}{(v_f^2/v_A^2)}
$$

法線の速波 Mach 数を $M_f$ とする。一との比較で判定する。

$$
M_f=\frac{(U_n/v_A)}{\sqrt{(v_f^2/v_A^2)}}
$$

法線の Alfvén Mach 数を $M_A$ とする。同じ速度でも速波の値とは異なる。

$$
M_A=(U_n/v_A)
$$

速波の超過の指標を $M_f^2-1$ とする。符号で超速波か判断する。

$$
M_f^2-1=M_{f}^{2} - 1
$$

:::example[垂直磁場]
音速と磁化の両方を含む Mach 数を求める。

入力値を次のように置く。

$$
\begin{aligned}
C&=\frac{3}{4} \\S&=0 \\U_n/v_A&=2
\end{aligned}
$$

これらを順に代入する。

速波と遅波を分ける判別式は

$$
D=- 4 \left(0\right) \left(\frac{3}{4}\right) + \left(\left(\frac{3}{4}\right) + 1\right)^{2}=\frac{49}{16}
$$

である。

速波の二乗速度比は

$$
v_f^2/v_A^2=\frac{\left(\frac{3}{4}\right)}{2} + \frac{\sqrt{\left(\frac{49}{16}\right)}}{2} + \frac{1}{2}=\frac{7}{4}
$$

である。

遅波の二乗速度比は

$$
v_s^2/v_A^2=\frac{\left(0\right) \left(\frac{3}{4}\right)}{\left(\frac{7}{4}\right)}=0
$$

である。

法線の速波 Mach 数は

$$
M_f=\frac{\left(2\right)}{\sqrt{\left(\frac{7}{4}\right)}}\simeq 1.5119
$$

である。

法線の Alfvén Mach 数は

$$
M_A=\left(2\right)=2
$$

である。

速波の超過の指標は

$$
M_f^2-1=\left(\frac{4 \sqrt{7}}{7}\right)^{2} - 1=\frac{9}{7}
$$

である。

上流は超速波であり、Alfvén Mach 数の方が大きい。
:::

:::example[平行磁場の境界]
上流速度を Alfvén 速度に合わせる。

入力値を次のように置く。

$$
\begin{aligned}
C&=\frac{3}{4} \\S&=1 \\U_n/v_A&=1
\end{aligned}
$$

これらを順に代入する。

速波と遅波を分ける判別式は

$$
D=- 4 \left(1\right) \left(\frac{3}{4}\right) + \left(\left(\frac{3}{4}\right) + 1\right)^{2}=\frac{1}{16}
$$

である。

速波の二乗速度比は

$$
v_f^2/v_A^2=\frac{\sqrt{\left(\frac{1}{16}\right)}}{2} + \frac{\left(\frac{3}{4}\right)}{2} + \frac{1}{2}=1
$$

である。

遅波の二乗速度比は

$$
v_s^2/v_A^2=\left(\frac{3}{4}\right)=\frac{3}{4}
$$

である。

法線の速波 Mach 数は

$$
M_f=\sqrt{\left(1\right)}=1
$$

である。

法線の Alfvén Mach 数は

$$
M_A=\left(1\right)=1
$$

である。

速波の超過の指標は

$$
M_f^2-1=\left(1\right)^{2} - 1=0
$$

である。

この入力では速波速度と同じで境界にある。
:::

:::example[同じ速度で斜め磁場]
磁場角を変えて境界を越える向きを調べる。

入力値を次のように置く。

$$
\begin{aligned}
C&=\frac{3}{4} \\S&=\frac{1}{2} \\U_n/v_A&=1
\end{aligned}
$$

これらを順に代入する。

速波と遅波を分ける判別式は

$$
D=- 4 \left(\frac{1}{2}\right) \left(\frac{3}{4}\right) + \left(\left(\frac{3}{4}\right) + 1\right)^{2}=\frac{25}{16}
$$

である。

速波の二乗速度比は

$$
v_f^2/v_A^2=\frac{\sqrt{\left(\frac{25}{16}\right)}}{2} + \frac{\left(\frac{3}{4}\right)}{2} + \frac{1}{2}=\frac{3}{2}
$$

である。

遅波の二乗速度比は

$$
v_s^2/v_A^2=\frac{\left(\frac{1}{2}\right) \left(\frac{3}{4}\right)}{\left(\frac{3}{2}\right)}=\frac{1}{4}
$$

である。

法線の速波 Mach 数は

$$
M_f=\frac{\left(1\right)}{\sqrt{\left(\frac{3}{2}\right)}}\simeq 0.8165
$$

である。

法線の Alfvén Mach 数は

$$
M_A=\left(1\right)=1
$$

である。

速波の超過の指標は

$$
M_f^2-1=\left(\frac{\sqrt{6}}{3}\right)^{2} - 1=- \frac{1}{3}
$$

である。

速波速度が増し、流速は速波を下回る。
:::

:::example[音速が大きい上流]
Alfvén Mach 数だけで衝撃波の分類ができないことを確かめる。

入力値を次のように置く。

$$
\begin{aligned}
C&=4 \\S&=\frac{1}{2} \\U_n/v_A&=2
\end{aligned}
$$

これらを順に代入する。

速波と遅波を分ける判別式は

$$
D=- 4 \left(4\right) \left(\frac{1}{2}\right) + \left(\left(4\right) + 1\right)^{2}=17
$$

である。

速波の二乗速度比は

$$
v_f^2/v_A^2=\frac{\sqrt{\left(17\right)}}{2} + \frac{\left(4\right)}{2} + \frac{1}{2}\simeq 4.5616
$$

である。

遅波の二乗速度比は

$$
v_s^2/v_A^2=\frac{\left(4\right) \left(\frac{1}{2}\right)}{\left(\frac{\sqrt{17}}{2} + \frac{5}{2}\right)}\simeq 0.43845
$$

である。

法線の速波 Mach 数は

$$
M_f=\frac{\left(2\right)}{\sqrt{\left(\frac{\sqrt{17}}{2} + \frac{5}{2}\right)}}\simeq 0.93643
$$

である。

法線の Alfvén Mach 数は

$$
M_A=\left(2\right)=2
$$

である。

速波の超過の指標は

$$
M_f^2-1=\left(\frac{2}{\sqrt{\frac{\sqrt{17}}{2} + \frac{5}{2}}}\right)^{2} - 1\simeq -0.12311
$$

である。

Alfvén 的には超過しても速波に対して不足する。
:::


## 太陽フレアと磁気再結合

太陽大気で急に熱と高速粒子のエネルギーが増える現象を :term[太陽フレア]{en="solar flare" yomi="たいようふれあ"} と呼ぶ。大気の一部では低 $\beta$ で磁気力が支配的である。足元の磁束分布を固定した系で、電流を持つ場はポテンシャル場より高いエネルギーを持ちうる。

### 利用可能な磁気エネルギー

滑らかな境界をもつ連結な有界領域で、元の磁場は $\nabla\cdot\vb B=0$ とする。このとき総法線磁束がゼロなので Neumann 問題の整合条件が満たされる。$\vb B_p=-\nabla\Phi$、$\nabla^2\Phi=0$ を解き、境界法線磁場を元の磁場と同じにする。差 $\vb b=\vb B-\vb B_p$ は発散も境界法線成分もゼロである。従って $\int\vb B_p\cdot\vb b\dd V=-\int\nabla\cdot(\Phi\vb b)\dd V=0$ であり
$$
\int\frac{B^2-B_p^2}{2\mu_0}\dd V=\int\frac{b^2}{2\mu_0}\dd V\ge0.
$$
右辺は境界法線磁束を保つ場の間で比較した自由磁気エネルギーである。理想発展では磁力線の連結が制約されるため、この全量が常に利用できるわけではない。再結合で連結が変わると、:eqref[eq-em-poynting] の $\vb j\cdot\vb E$ を通じて場のエネルギーを粒子へ渡せる。流出での衝撃波や波も、熱と非熱成分への配分を変える。

:::example[場の尺度]
余剰成分 $b=0.01\,\mathrm T$ が体積 $10^{21}\,\mathrm{m^3}$ にほぼ一様なら、自由エネルギーは $b^2V/(2\mu_0)\simeq4\times10^{22}\,\mathrm J$ である。これは磁場模型の蓄積量であり、放射エネルギーは配分を別に計算する。
:::

:::example[再結合の時間尺度]
長さ $L$ と上流 Alfvén 速度 $v_A$ の層で、全長を横切る量を代表として流入時間を $L/v_{\rm in}$ と見積もる。$v_{\rm in}/v_A=10^{-2}$ なら $100 L/v_A$ であり、$10^{-4}$ なら $10^4 L/v_A$ である。両率で経過時間は百倍異なる。観測で急に散逸する量を議論するときは、有効層長と流入域を別に特定する。
:::

### フレアの流入 Poynting 流束と電力

案内磁場のない対称な再結合層で、上流の反転する磁場を $B$、密度を $\rho$ とする。流入は磁場と直交するので、理想式 $\vb E=-\vb u\times\vb B$ による流入 Poynting 流束の大きさは $v_{\rm in}B^2/\mu_0$ である。片側の面積を $A$ とすれば両面から供給される電力は $2A v_{\rm in}B^2/\mu_0$ である。係数二は二つの面の和であり、磁気エネルギー密度 $B^2/(2\mu_0)$ を単に移流させた流束とも区別する。流束には磁気応力の仕事が含まれる。定常でもこの入力が全て放射になるとは限らず、流出の運動エネルギーや熱伝導も収支に含める。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $B$ | 上流の反転磁場 | $\mathrm T$ |
| $A$ | 片側の流入面積 | $\mathrm{m^2}$ |
| $v_A$ | 上流の Alfvén 速度 | $\mathrm{m\,s^{-1}}$ |
| $M$ | 流入比 $v_{\rm in}/v_A$ | 1 |

流入速度を $v_{\rm in}$ とする。流出の尺度を $v_A$ に取る。

$$
v_{\rm in}=M v_{A}
$$

再結合電場の大きさを $E_{\rm rec}$ とする。理想的な上流の値である。

$$
E_{\rm rec}=B v_{\rm in}
$$

規格化した片面の電磁流束を $\mu_0 S_{\rm in}$ とする。$\vb E$ と $\vb B$ の外積を取る。

$$
\mu_0 S_{\rm in}=B E_{\rm rec}
$$

両面から供給する電力を $P_{\rm mag}$ とする。$\mu_0\simeq4\pi\times10^{-7}\,\mathrm{H\,m^{-1}}$ とする。

$$
P_{\rm mag}=2A S_{\rm in}
$$

六十秒間の供給エネルギーを $W_{60}$ とする。定常かつ一定面積の形式値である。

$$
W_{60}=60 P_{\rm mag}
$$

上流の磁気エネルギー密度を $U_B$ とする。電力と直接同じ量ではない。

$$
U_B=\frac{B^2}{2\mu_0}
$$

:::example[小さい流入率]
上流の電場と電力を求める。

入力値を次のように置く。

$$
\begin{aligned}
B&=\frac{1}{1000}\,\mathrm T \\A&=1000000000000\,\mathrm{m^2} \\v_A&=1000000\,\mathrm{m\,s^{-1}} \\M&=\frac{1}{1000}
\end{aligned}
$$

これらを順に代入する。

流入速度は

$$
v_{\rm in}=\left(1000000\right) \left(\frac{1}{1000}\right)=1000\,\mathrm{m\,s^{-1}}
$$

である。

再結合電場の大きさは

$$
E_{\rm rec}=\left(1000\right) \left(\frac{1}{1000}\right)=1\,\mathrm{V\,m^{-1}}
$$

である。

規格化した片面の電磁流束は

$$
\mu_0 S_{\rm in}=\left(1\right) \left(\frac{1}{1000}\right)=\frac{1}{1000}\,\mathrm{T^2\,m\,s^{-1}}
$$

である。

両面から供給する電力は

$$
P_{\rm mag}=\frac{5000000 \left(1000000000000\right) \left(\frac{1}{1000}\right)}{\pi}\simeq 1.5915\times 10^{15}\,\mathrm W
$$

である。

六十秒間の供給エネルギーは

$$
W_{60}=60 \left(\frac{5000000000000000}{\pi}\right)\simeq 9.5493\times 10^{16}\,\mathrm J
$$

である。

上流の磁気エネルギー密度は

$$
U_B=\frac{1250000 \left(\frac{1}{1000}\right)^{2}}{\pi}\simeq 0.39789\,\mathrm{J\,m^{-3}}
$$

である。

流出より流入が遅いので、電力に入る速さも小さい。
:::

:::example[一桁速い流入]
同じ磁気エネルギー密度でも供給が増えるか確かめる。

入力値を次のように置く。

$$
\begin{aligned}
B&=\frac{1}{1000}\,\mathrm T \\A&=1000000000000\,\mathrm{m^2} \\v_A&=1000000\,\mathrm{m\,s^{-1}} \\M&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

流入速度は

$$
v_{\rm in}=\left(1000000\right) \left(\frac{1}{100}\right)=10000\,\mathrm{m\,s^{-1}}
$$

である。

再結合電場の大きさは

$$
E_{\rm rec}=\left(10000\right) \left(\frac{1}{1000}\right)=10\,\mathrm{V\,m^{-1}}
$$

である。

規格化した片面の電磁流束は

$$
\mu_0 S_{\rm in}=\left(10\right) \left(\frac{1}{1000}\right)=\frac{1}{100}\,\mathrm{T^2\,m\,s^{-1}}
$$

である。

両面から供給する電力は

$$
P_{\rm mag}=\frac{5000000 \left(1000000000000\right) \left(\frac{1}{100}\right)}{\pi}\simeq 1.5915\times 10^{16}\,\mathrm W
$$

である。

六十秒間の供給エネルギーは

$$
W_{60}=60 \left(\frac{50000000000000000}{\pi}\right)\simeq 9.5493\times 10^{17}\,\mathrm J
$$

である。

上流の磁気エネルギー密度は

$$
U_B=\frac{1250000 \left(\frac{1}{1000}\right)^{2}}{\pi}\simeq 0.39789\,\mathrm{J\,m^{-3}}
$$

である。

時間当たりの供給は十倍で、エネルギー密度は変わらない。
:::

:::example[二倍の磁場]
$v_A$ が同じになるよう密度も調整して比較する。

入力値を次のように置く。

$$
\begin{aligned}
B&=\frac{1}{500}\,\mathrm T \\A&=1000000000000\,\mathrm{m^2} \\v_A&=1000000\,\mathrm{m\,s^{-1}} \\M&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

流入速度は

$$
v_{\rm in}=\left(1000000\right) \left(\frac{1}{100}\right)=10000\,\mathrm{m\,s^{-1}}
$$

である。

再結合電場の大きさは

$$
E_{\rm rec}=\left(10000\right) \left(\frac{1}{500}\right)=20\,\mathrm{V\,m^{-1}}
$$

である。

規格化した片面の電磁流束は

$$
\mu_0 S_{\rm in}=\left(20\right) \left(\frac{1}{500}\right)=\frac{1}{25}\,\mathrm{T^2\,m\,s^{-1}}
$$

である。

両面から供給する電力は

$$
P_{\rm mag}=\frac{5000000 \left(1000000000000\right) \left(\frac{1}{25}\right)}{\pi}\simeq 6.3662\times 10^{16}\,\mathrm W
$$

である。

六十秒間の供給エネルギーは

$$
W_{60}=60 \left(\frac{200000000000000000}{\pi}\right)\simeq 3.8197\times 10^{18}\,\mathrm J
$$

である。

上流の磁気エネルギー密度は

$$
U_B=\frac{1250000 \left(\frac{1}{500}\right)^{2}}{\pi}\simeq 1.5915\,\mathrm{J\,m^{-3}}
$$

である。

電場は二倍で、Poynting 流束は四倍になる。
:::

:::example[流入面積が広い層]
局所電場と全体の入力を比較する。

入力値を次のように置く。

$$
\begin{aligned}
B&=\frac{1}{1000}\,\mathrm T \\A&=10000000000000\,\mathrm{m^2} \\v_A&=1000000\,\mathrm{m\,s^{-1}} \\M&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

流入速度は

$$
v_{\rm in}=\left(1000000\right) \left(\frac{1}{100}\right)=10000\,\mathrm{m\,s^{-1}}
$$

である。

再結合電場の大きさは

$$
E_{\rm rec}=\left(10000\right) \left(\frac{1}{1000}\right)=10\,\mathrm{V\,m^{-1}}
$$

である。

規格化した片面の電磁流束は

$$
\mu_0 S_{\rm in}=\left(10\right) \left(\frac{1}{1000}\right)=\frac{1}{100}\,\mathrm{T^2\,m\,s^{-1}}
$$

である。

両面から供給する電力は

$$
P_{\rm mag}=\frac{5000000 \left(10000000000000\right) \left(\frac{1}{100}\right)}{\pi}\simeq 1.5915\times 10^{17}\,\mathrm W
$$

である。

六十秒間の供給エネルギーは

$$
W_{60}=60 \left(\frac{500000000000000000}{\pi}\right)\simeq 9.5493\times 10^{18}\,\mathrm J
$$

である。

上流の磁気エネルギー密度は

$$
U_B=\frac{1250000 \left(\frac{1}{1000}\right)^{2}}{\pi}\simeq 0.39789\,\mathrm{J\,m^{-3}}
$$

である。

電場は同じだが、系全体で供給される電力が大きくなる。
:::


## 衝撃波統計加速（Fermi 加速）と宇宙線

衝撃波の近くで粒子が上流と下流を何度も往復し、両側の散乱体の相対運動からエネルギーを得る過程を :term[一次 Fermi 加速]{en="first-order Fermi acceleration" yomi="いちじふぇるみかそく"} と呼ぶ。平面定常な非相対論的衝撃波の静止系で上流・下流の流速を $u_1>u_2>0$ とする。粒子速さ $v\gg u_1$、両側の散乱がほぼ等方で弾性的、粒子が流れを変えないと仮定する。

### 一往復の利得と逃げ

座標系の速度差 $\Delta u=u_1-u_2$ による運動量の相対変化は、一次で $\Delta p/p=(\Delta u/v)\mu$ である。これは $\delta E=\Delta u\,p\mu$ と $\dd E/\dd p=v$ を使った。衝撃面を横切る粒子の余弦 $\mu\in(0,1)$ の分布は流束で重み付けされ $2\mu$ であるから $\langle\mu\rangle=2/3$。往路と復路で利得は加算され、一往復当たり $g=4(u_1-u_2)/(3v)$ となる。下流へ運び去られる粒子束は $nu_2$、面を横切る等方な粒子束は $nv/4$ なので、一往復での逃げ確率は $P_{\rm esc}\simeq4u_2/v$ である。

一往復による $\ln p$ の増加は $g$、生存数 $\mathcal N$ の対数の減少は $P_{\rm esc}$ である。$g,P_{\rm esc}\ll1$ の範囲では、$v$ が $p$ に依存しても比 $P_{\rm esc}/g$ は一定である。従って $\dd\ln\mathcal N/\dd\ln p\simeq-P_{\rm esc}/g$ を積分すると積分スペクトルは $\mathcal N(>p)\propto p^{-3/(r-1)}$、$r=u_1/u_2$ である。等方分布の $\mathcal N$ は $4\pi\int_p^\infty p'^2 f(p')\dd p'$ なので、
$$
f(p)\propto p^{-3r/(r-1)}.
\label{eq-plasma-fermi-index}
$$
相対論的な粒子では $E\simeq pc$ だから $r=4$ の微分エネルギースペクトルは $E^{-2}$ となる。衝撃波の有限寿命、粒子の拡散長、放射損失は最大エネルギーを制限する。この局所スペクトルと観測される宇宙線の分布の間には、銀河内の輸送も入る。

:::example[弱い圧縮]
$r=2$ なら $f\propto p^{-6}$ であり、相対論的な粒子の微分エネルギースペクトルは $E^{-4}$ となる。下流へ逃げる効果が利得に比べて大きいほど、分布は急になる。
:::

### 拡散・移流方程式によるスペクトルの導出

同じ試験粒子近似で、等方な分布 $f(x,p)$ は衝撃面を含む一次元の定常拡散・移流方程式
$$
u\partial_x f=\partial_x(\kappa\,\partial_x f)+\frac{\partial_x u}{3}p\partial_p f
$$
を満たす。圧縮項の係数は、局所の等方な粒子集団で各線寸法が同じ比率で縮むとき $p$ が体積の $-1/3$ 乗に比例することから得られる。衝撃面を $x=0$、上流を $x<0$ にとる。注入運動量より上で上流遠方の分布をゼロとし、下流遠方の分布が有界であるとする。両側で $u,\kappa$ を位置によらない量とすれば、面で連続な解は上流で $f_0(p)e^{u_1x/\kappa_1}$、下流で $f_0(p)$ である。面を横切って方程式を積分すると、$u\partial_x f$ の有限な積分は消え、
$$
0=-u_1 f_0+\frac{u_2-u_1}{3}p\dv{f_0}{p}
$$
となる。従って :eqref[eq-plasma-fermi-index] を再び得る。拡散係数がスペクトルの指数から消えるのは、この境界条件と定常性による。

### 拡散で決まる加速時間

各側で衝撃面の法線方向の空間拡散係数を $\kappa_j$、流速を $u_j$ とする。定常な上流の一次元移流・拡散式 $u_1 f_x=\kappa_1 f_{xx}$ の減衰解から、上流への広がりは $\kappa_1/u_1$ である。下流の戻り確率も、距離 $x$ からの逆向き拡散と移流の競合により $e^{-u_2x/\kappa_2}$ で下がる。これは $u_2P'+\kappa_2 P''=0$ に $P(0)=1,P(\infty)=0$ を課した解である。面を渡る束が $nv/4$ であることから、滞在時間の尺度は両側で $4\kappa_j/(vu_j)$ となる。一往復の利得 $g$ で割ると、運動量の対数が一だけ増えるための局所的な時間尺度は
$$
t_{\rm acc}(p)=\frac{3}{u_1-u_2}\left(\frac{\kappa_1(p)}{u_1}+\frac{\kappa_2(p)}{u_2}\right).
$$
拡散が $p$ に依存する場合、$p_0$ から $p$ までの時間は $t_{\rm acc}\dd\ln p$ を積分する。無秩序な散乱体が平均的に静止している二次 Fermi 過程では一次の利得が向きの平均で消え、残る利得は速度差の二次程度である。衝撃波では二方向からの横断が頭合わせとなり、一次の平均が残る。

:::example[等しい拡散係数]
$\kappa_1=\kappa_2=\kappa_d$、$u_1=4u_2$ なら $t_{\rm acc}=20\kappa_d/u_1^2$ である。拡散係数が大きいと衝撃波から遠くまで移動して一往復の時間が長くなる。
:::

### 有限往復利得を保持した Fermi スペクトル

一往復で $p$ を $1+g$ 倍し、確率 $P$ で失う独立な試験粒子模型を考える。整数 $N$ 回の計算では生存率は $(1-P)^N$ である。連続包絡をとれば $\mathcal N(>p)\propto p^{-a}$、$a=-\ln(1-P)/\ln(1+g)$ となる。$g,P$ が小さいと $a\simeq P/g$ である。微分運動量スペクトルは $p^{-(a+1)}$、等方な三次元速度空間での位相密度は $p^{-(a+3)}$ である。以下の有限値での差は離散確率模型の結果であり、衝撃波の一次の $g,P$ の導出を有限 $u/v$ に延長した値ではない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $g$ | 一往復での相対運動量増加 | 1 |
| $P$ | 一往復で逃げる確率 | 1 |
| $n$ | 往復回数 | 1 |

運動量の倍率を $p/p_0$ とする。利得が毎回同じとする。

$$
p/p_0=\left(g + 1\right)^{n}
$$

生じた生存率を $\mathcal N_n/\mathcal N_0$ とする。逃げる確率を乗算する。

$$
\mathcal N_n/\mathcal N_0=\left(1 - P\right)^{n}
$$

積分スペクトルの指数を $a$ とする。回数を消去して求める。

$$
a=- \frac{\ln{\left(1 - P \right)}}{\ln{\left(g + 1 \right)}}
$$

一次近似の指数を $a_{\rm lin}$ とする。対数を小引数で展開する。

$$
a_{\rm lin}=\frac{P}{g}
$$

微分スペクトルの指数を $s_p$ とする。$-\dd\mathcal N/\dd p$ を計算する。

$$
s_p=a + 1
$$

位相密度の指数を $s_f$ とする。$4\pi p^2$ の球殻因子で割る。

$$
s_f=a + 3
$$

:::example[弱い利得と弱い散逸]
一次近似と有限往復の指数を比べる。

入力値を次のように置く。

$$
\begin{aligned}
g&=\frac{1}{100} \\P&=\frac{1}{100} \\n&=100
\end{aligned}
$$

これらを順に代入する。

運動量の倍率は

$$
p/p_0=\left(\left(\frac{1}{100}\right) + 1\right)^{\left(100\right)}\simeq 2.7048
$$

である。

生じた生存率は

$$
\mathcal N_n/\mathcal N_0=\left(1 - \left(\frac{1}{100}\right)\right)^{\left(100\right)}\simeq 0.36603
$$

である。

積分スペクトルの指数は

$$
a=- \frac{\ln{\left(1 - \left(\frac{1}{100}\right) \right)}}{\ln{\left(\left(\frac{1}{100}\right) + 1 \right)}}\simeq 1.0101
$$

である。

一次近似の指数は

$$
a_{\rm lin}=1
$$

である。

微分スペクトルの指数は

$$
s_p=\left(- \frac{\ln{\left(\frac{99}{100} \right)}}{\ln{\left(\frac{101}{100} \right)}}\right) + 1\simeq 2.0101
$$

である。

位相密度の指数は

$$
s_f=\left(- \frac{\ln{\left(\frac{99}{100} \right)}}{\ln{\left(\frac{101}{100} \right)}}\right) + 3\simeq 4.0101
$$

である。

小量が等しい場合でも対数の符号の違いで指数に補正がある。
:::

:::example[一往復の利得が逃げより大きい]
有限回での生存と運動量を求める。

入力値を次のように置く。

$$
\begin{aligned}
g&=\frac{1}{10} \\P&=\frac{1}{100} \\n&=10
\end{aligned}
$$

これらを順に代入する。

運動量の倍率は

$$
p/p_0=\left(\left(\frac{1}{10}\right) + 1\right)^{\left(10\right)}\simeq 2.5937
$$

である。

生じた生存率は

$$
\mathcal N_n/\mathcal N_0=\left(1 - \left(\frac{1}{100}\right)\right)^{\left(10\right)}\simeq 0.90438
$$

である。

積分スペクトルの指数は

$$
a=- \frac{\ln{\left(1 - \left(\frac{1}{100}\right) \right)}}{\ln{\left(\left(\frac{1}{10}\right) + 1 \right)}}\simeq 0.10545
$$

である。

一次近似の指数は

$$
a_{\rm lin}=\frac{\left(\frac{1}{100}\right)}{\left(\frac{1}{10}\right)}=\frac{1}{10}
$$

である。

微分スペクトルの指数は

$$
s_p=\left(- \frac{\ln{\left(\frac{99}{100} \right)}}{\ln{\left(\frac{11}{10} \right)}}\right) + 1\simeq 1.1054
$$

である。

位相密度の指数は

$$
s_f=\left(- \frac{\ln{\left(\frac{99}{100} \right)}}{\ln{\left(\frac{11}{10} \right)}}\right) + 3\simeq 3.1054
$$

である。

生存が多いので積分スペクトルは緩やかになる。
:::

:::example[逃げる確率が大きい]
一次の指数近似からどれだけずれるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
g&=\frac{1}{10} \\P&=\frac{1}{2} \\n&=10
\end{aligned}
$$

これらを順に代入する。

運動量の倍率は

$$
p/p_0=\left(\left(\frac{1}{10}\right) + 1\right)^{\left(10\right)}\simeq 2.5937
$$

である。

生じた生存率は

$$
\mathcal N_n/\mathcal N_0=\left(1 - \left(\frac{1}{2}\right)\right)^{\left(10\right)}=\frac{1}{1024}
$$

である。

積分スペクトルの指数は

$$
a=- \frac{\ln{\left(1 - \left(\frac{1}{2}\right) \right)}}{\ln{\left(\left(\frac{1}{10}\right) + 1 \right)}}\simeq 7.2725
$$

である。

一次近似の指数は

$$
a_{\rm lin}=\frac{\left(\frac{1}{2}\right)}{\left(\frac{1}{10}\right)}=5
$$

である。

微分スペクトルの指数は

$$
s_p=\left(\frac{\ln{\left(2 \right)}}{\ln{\left(\frac{11}{10} \right)}}\right) + 1\simeq 8.2725
$$

である。

位相密度の指数は

$$
s_f=\left(\frac{\ln{\left(2 \right)}}{\ln{\left(\frac{11}{10} \right)}}\right) + 3\simeq 10.273
$$

である。

逃げが大きいので対数の近似は弱い。
:::

:::example[強い衝撃波の一次模型に近い値]
微分指数が二に近づくことを確認する。

入力値を次のように置く。

$$
\begin{aligned}
g&=\frac{1}{1000} \\P&=\frac{1}{1000} \\n&=1000
\end{aligned}
$$

これらを順に代入する。

運動量の倍率は

$$
p/p_0=\left(\left(\frac{1}{1000}\right) + 1\right)^{\left(1000\right)}\simeq 2.7169
$$

である。

生じた生存率は

$$
\mathcal N_n/\mathcal N_0=\left(1 - \left(\frac{1}{1000}\right)\right)^{\left(1000\right)}\simeq 0.3677
$$

である。

積分スペクトルの指数は

$$
a=- \frac{\ln{\left(1 - \left(\frac{1}{1000}\right) \right)}}{\ln{\left(\left(\frac{1}{1000}\right) + 1 \right)}}\simeq 1.001
$$

である。

一次近似の指数は

$$
a_{\rm lin}=1
$$

である。

微分スペクトルの指数は

$$
s_p=\left(- \frac{\ln{\left(\frac{999}{1000} \right)}}{\ln{\left(\frac{1001}{1000} \right)}}\right) + 1\simeq 2.001
$$

である。

位相密度の指数は

$$
s_f=\left(- \frac{\ln{\left(\frac{999}{1000} \right)}}{\ln{\left(\frac{1001}{1000} \right)}}\right) + 3\simeq 4.001
$$

である。

位相密度の指数は四に近く、エネルギーの微分指数と混同しない。
:::


## 降着円盤と磁気回転不安定性

円盤の内側が速く回る場合、弱い磁場が結びつけた二つの流体要素が互いに角運動量を交換する。内側が角運動量を失ってさらに内へ、外側が得てさらに外へ動き、ずれが増幅されうる。これが :term[磁気回転不安定性]{en="magnetorotational instability" yomi="じきかいてんふあんていせい"} である。

### 局所の線形方程式

半径 $R_0$ とともに $\Omega$ で回転する直交座標 $(x,y,z)$ をとる。$x$ は外向き、$y$ は回転方向、$z$ は円盤に垂直である。背景速度を $-q_\Omega\Omega x\hat{\vb y}$、$q_\Omega=-\dd\ln\Omega/\dd\ln R$ とする。密度一定、背景磁場は鉛直、軸対称摂動は $e^{ikz+\gamma t}$ とする。磁場摂動 $\vb b$ は速度単位に換算する。$a=kv_A$ とすれば、線形の運動式と誘導式は
$$
\begin{aligned}
\gamma u_x-2\Omega u_y&=ia b_x,&
\gamma u_y+(2-q_\Omega)\Omega u_x&=ia b_y,\\
\gamma b_x&=ia u_x,&
\gamma b_y&=-q_\Omega\Omega b_x+ia u_y .
\end{aligned}
$$
$-2\Omega u_y$ と $2\Omega u_x$ は Coriolis 力から来る。$-q_\Omega\Omega u_x$ と $-q_\Omega\Omega b_x$ は背景せん断を横切る効果である。誘導式を置換すると $(\gamma^2+a^2)u_x=2\Omega\gamma u_y$、$(\gamma^2+a^2)u_y+(2-q_\Omega)\Omega\gamma u_x-q_\Omega\Omega a^2 u_x/\gamma=0$ となる。これらを消去して
$$
\gamma^4+\gamma^2(\kappa^2+2a^2)+a^2(a^2-2q_\Omega\Omega^2)=0,\qquad
\kappa^2=2(2-q_\Omega)\Omega^2.
\label{eq-plasma-mri}
$$
を得る。$\kappa^2>0$ で磁場のない円盤は軸対称の動径振動に安定である。$q_\Omega>0$ のとき $0<a^2<2q_\Omega\Omega^2$ なら $\gamma^2$ の二根の積が負となり、正の実成長根が一つある。張力がない $a=0$ では結合がなく、張力が強すぎると摂動は安定になる。

:::example[Kepler 円盤]
$\Omega\propto R^{-3/2}$ なら $q_\Omega=3/2$、$\kappa=\Omega$ である。$b=a^2/\Omega^2$ とすれば成長枝の二乗は $\gamma^2/\Omega^2=(-1-2b+\sqrt{1+16b})/2$ である。微分は $-1+4/\sqrt{1+16b}$ なので $b=15/16$ で最大となり、$\gamma_{\max}=3\Omega/4$ を得る。
:::

### 角運動量の流束

局所の速度単位の磁場摂動を用い、平面平均した流れの動径方向の角運動量輸送を考える。せん断応力は $\mathcal T_{xy}=\rho\langle u_x u_y-b_x b_y\rangle$ であり、第 1 項は速度の運動量流束、第 2 項は磁気応力である。周期領域で軸対称摂動の二次エネルギーを積分する。速度式の反対称な Coriolis 項は一対で消え、磁気張力と誘導による仕事も部分積分で相殺する。残るせん断項は $q_\Omega\Omega\int \mathcal T_{xy}\dd V$ である。成長するモードではこの積が正であるため、$q_\Omega>0$ のとき平均応力は角運動量を外へ運ぶ向きである。飽和した乱流の応力の大きさは線形成長率だけからは決まらない。

:::example[許される波長]
高さ $H$ の周期箱なら $k=2\pi j/H$、$j=1,2,\ldots$ である。Kepler 円盤で $2\pi v_A/H\ge\sqrt3\,\Omega$ なら、最低次の非零波数からすでに限界以上なので、この箱に軸対称の鉛直成長モードは入らない。
:::

### 抵抗がある場合の成長境界

運動式はそのままとし、誘導式に磁気拡散係数 $\eta_m>0$ を加えた局所模型を考える。$d=\eta_m k^2$ とすると誘導式左辺の $\gamma$ だけが $\gamma+d$ へ変わる。$b_x=ia u_x/(\gamma+d)$ を先に消去し、残る $b_y$ と $u_y$ を消去すれば
$$
[\gamma(\gamma+d)+a^2]^2+\kappa^2(\gamma+d)^2-2q_\Omega\Omega^2a^2=0
$$
である。$0<q_\Omega<2$ なら $\kappa^2>0$ であり、左辺は実数 $\gamma\ge0$ で単調に増加する。従って $\gamma=0$ での値 $c_0=a^4+\kappa^2d^2-2q_\Omega\Omega^2a^2$ が負のとき、正の実成長根を一つ持つ。

$c_0>0$ では複素成長根もない。実際、四次式の三次・二次・一次係数は $c_3=2d$、$c_2=d^2+2a^2+\kappa^2$、$c_1=2d(a^2+\kappa^2)$ である。Routh 表の最初の列は $1,c_3,(c_3c_2-c_1)/c_3,(c_3c_2c_1-c_1^2-c_3^2c_0)/(c_3c_2-c_1),c_0$ である。$a,d>0$ なら二つの差は $2d(d^2+a^2)$ と $4d^2a^2(d^2+\kappa^2+2q_\Omega\Omega^2)$ で、ともに正である。したがって Routh の安定判定で右半平面の根はない。例えば Kepler 回転で $kv_A=\Omega$ なら $c_0=\Omega^2(d^2-2\Omega^2)$ となり、磁気拡散率を増やして $d=\sqrt2\Omega$ に達するとこの波数の成長は止まる。粘性や Hall 項を足した場合には、この境界を使う前に誘導式と運動式を組み直す必要がある。

### 回転勾配と MRI の最速波数

$q=-\dd\ln\Omega/\dd\ln R$ を $0<q<2$ とする。:eqref[eq-plasma-mri] で $b=(kv_A/\Omega)^2$、$K=2(2-q)$ とおくと、成長枝は $G=\gamma^2/\Omega^2=[-K-2b+\sqrt{K^2+16b}]/2$ である。判別式の $16b$ は $K+2q=4$ から出る。$\partial_bG=-1+4/\sqrt{K^2+16b}$ なので $b_*=q(4-q)/4$、$G_*=q^2/4$ となる。$\partial_b^2G<0$ が極大を保証する。$b=0$ で $G=0$、$b=2q$ で再びゼロである。最速波長が円盤の高さより長くなる場合には局所の最適波数を使えない。抵抗や粘性のある円盤では最短波長側の条件も別に現れる。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $q$ | 回転の対数勾配の符号を反転した量 | 1 |
| $b$ | 鉛直の Alfvén 波数比の二乗 | 1 |

エピサイクル角周波数の二乗を $\kappa^2/\Omega^2$ とする。非磁化の回転の安定性を判定する。

$$
\kappa^2/\Omega^2=4 - 2 q
$$

磁気張力のしきい波数の二乗を $b_c$ とする。$0<b<b_c$ が成長帯である。

$$
b_c=2 q
$$

成長候補枝の二乗値を $\gamma^2/\Omega^2$ とする。負なら振動の枝である。

$$
\gamma^2/\Omega^2=- \frac{\kappa^2/\Omega^2}{2} - b + \frac{\sqrt{\left(\kappa^2/\Omega^2\right)^{2} + 16 b}}{2}
$$

最速波数の二乗値を $b_*$ とする。二乗成長率を波数で微分する。

$$
b_*=\frac{q \left(4 - q\right)}{4}
$$

最速成長率の比を $\gamma_*/\Omega$ とする。連続波数で許されると仮定する。

$$
\gamma_*/\Omega=\frac{q}{2}
$$

最大値を基準にした形式二乗値を $(\gamma/\gamma_*)^2$ とする。安定側では負になる。

$$
(\gamma/\gamma_*)^2=\frac{\gamma^2/\Omega^2}{(\gamma_*/\Omega)^2}
$$

:::example[Kepler 回転の最速波数]
一般回転勾配の式が Kepler の結果に戻るか調べる。

入力値を次のように置く。

$$
\begin{aligned}
q&=\frac{3}{2} \\b&=\frac{15}{16}
\end{aligned}
$$

これらを順に代入する。

エピサイクル角周波数の二乗は

$$
\kappa^2/\Omega^2=4 - 2 \left(\frac{3}{2}\right)=1
$$

である。

磁気張力のしきい波数の二乗は

$$
b_c=2 \left(\frac{3}{2}\right)=3
$$

である。

成長候補枝の二乗値は

$$
\gamma^2/\Omega^2=- \frac{\left(1\right)}{2} - \left(\frac{15}{16}\right) + \frac{\sqrt{\left(1\right)^{2} + 16 \left(\frac{15}{16}\right)}}{2}=\frac{9}{16}
$$

である。

最速波数の二乗値は

$$
b_*=\frac{\left(\frac{3}{2}\right) \left(4 - \left(\frac{3}{2}\right)\right)}{4}=\frac{15}{16}
$$

である。

最速成長率の比は

$$
\gamma_*/\Omega=\frac{\left(\frac{3}{2}\right)}{2}=\frac{3}{4}
$$

である。

最大値を基準にした形式二乗値は

$$
(\gamma/\gamma_*)^2=\frac{\left(\frac{9}{16}\right)}{\left(\frac{3}{4}\right)^{2}}=1
$$

である。

成長率は回転率の四分の三である。
:::

:::example[磁気張力が強い波]
同じ回転で短い波長を調べる。

入力値を次のように置く。

$$
\begin{aligned}
q&=\frac{3}{2} \\b&=4
\end{aligned}
$$

これらを順に代入する。

エピサイクル角周波数の二乗は

$$
\kappa^2/\Omega^2=4 - 2 \left(\frac{3}{2}\right)=1
$$

である。

磁気張力のしきい波数の二乗は

$$
b_c=2 \left(\frac{3}{2}\right)=3
$$

である。

成長候補枝の二乗値は

$$
\gamma^2/\Omega^2=- \frac{\left(1\right)}{2} - \left(4\right) + \frac{\sqrt{\left(1\right)^{2} + 16 \left(4\right)}}{2}\simeq -0.46887
$$

である。

最速波数の二乗値は

$$
b_*=\frac{\left(\frac{3}{2}\right) \left(4 - \left(\frac{3}{2}\right)\right)}{4}=\frac{15}{16}
$$

である。

最速成長率の比は

$$
\gamma_*/\Omega=\frac{\left(\frac{3}{2}\right)}{2}=\frac{3}{4}
$$

である。

最大値を基準にした形式二乗値は

$$
(\gamma/\gamma_*)^2=\frac{\left(- \frac{9}{2} + \frac{\sqrt{65}}{2}\right)}{\left(\frac{3}{4}\right)^{2}}\simeq -0.83355
$$

である。

張力のしきいを超え成長しない。
:::

:::example[平坦な回転曲線]
周速度一定のときの最速値を求める。

入力値を次のように置く。

$$
\begin{aligned}
q&=1 \\b&=\frac{3}{4}
\end{aligned}
$$

これらを順に代入する。

エピサイクル角周波数の二乗は

$$
\kappa^2/\Omega^2=4 - 2 \left(1\right)=2
$$

である。

磁気張力のしきい波数の二乗は

$$
b_c=2 \left(1\right)=2
$$

である。

成長候補枝の二乗値は

$$
\gamma^2/\Omega^2=- \frac{\left(2\right)}{2} - \left(\frac{3}{4}\right) + \frac{\sqrt{\left(2\right)^{2} + 16 \left(\frac{3}{4}\right)}}{2}=\frac{1}{4}
$$

である。

最速波数の二乗値は

$$
b_*=\frac{\left(1\right) \left(4 - \left(1\right)\right)}{4}=\frac{3}{4}
$$

である。

最速成長率の比は

$$
\gamma_*/\Omega=\frac{\left(1\right)}{2}=\frac{1}{2}
$$

である。

最大値を基準にした形式二乗値は

$$
(\gamma/\gamma_*)^2=\frac{\left(\frac{1}{4}\right)}{\left(\frac{1}{2}\right)^{2}}=1
$$

である。

成長率の最大は Kepler の場合より低い。
:::

:::example[ほぼ剛体回転]
回転が速くても勾配が弱い場合を求める。

入力値を次のように置く。

$$
\begin{aligned}
q&=\frac{1}{10} \\b&=1
\end{aligned}
$$

これらを順に代入する。

エピサイクル角周波数の二乗は

$$
\kappa^2/\Omega^2=4 - 2 \left(\frac{1}{10}\right)=\frac{19}{5}
$$

である。

磁気張力のしきい波数の二乗は

$$
b_c=2 \left(\frac{1}{10}\right)=\frac{1}{5}
$$

である。

成長候補枝の二乗値は

$$
\gamma^2/\Omega^2=- \left(1\right) - \frac{\left(\frac{19}{5}\right)}{2} + \frac{\sqrt{16 \left(1\right) + \left(\frac{19}{5}\right)^{2}}}{2}\simeq -0.14138
$$

である。

最速波数の二乗値は

$$
b_*=\frac{\left(\frac{1}{10}\right) \left(4 - \left(\frac{1}{10}\right)\right)}{4}=\frac{39}{400}
$$

である。

最速成長率の比は

$$
\gamma_*/\Omega=\frac{\left(\frac{1}{10}\right)}{2}=\frac{1}{20}
$$

である。

最大値を基準にした形式二乗値は

$$
(\gamma/\gamma_*)^2=\frac{\left(- \frac{29}{10} + \frac{\sqrt{761}}{10}\right)}{\left(\frac{1}{20}\right)^{2}}\simeq -56.551
$$

である。

成長帯の短波長端が移動し、設定した波は安定側となる。
:::


## 演習問題

::::exercise{level=basic}
$c_s$ が二倍になれば Parker 風の $r_c$ は何倍か。

:::solution
$r_c=GM/(2c_s^2)$ なので四分の一である。
:::
::::

::::exercise{level=basic}
動圧が変わらず双極子の基準磁場 $B_*$ が八倍のとき $r_{\rm mp}$ は何倍か。

:::solution
磁気圧は $64$ 倍であり、その六乗根を取るため二倍である。
:::
::::

::::exercise{level=basic}
衝撃波で $r=4$ のとき $f(p)$ の指数を求めよ。

:::solution
式 :eqref[eq-plasma-fermi-index] の指数は $-3r/(r-1)=-12/3=-4$ である。
:::
::::

::::exercise{level=basic}
Kepler 円盤で $kv_A=2\Omega$ の MRI 成長根は存在するか。

:::solution
$a^2=4\Omega^2$ はしきい値 $3\Omega^2$ を上回る。二次式の $\gamma^2$ の根は実で和が負、積が正なのでどちらも負で、成長根はない。
:::
::::

::::exercise{level=advanced}
Parker の積分式の臨界点近くで $m-1\simeq x-1$ を確認せよ。

:::solution
左辺は $1+2(m-1)^2+\order{(m-1)^3}$、右辺は $1+2(x-1)^2+\order{(x-1)^3}$ である。亜音速から超音速へ進む枝で同じ符号を選べば $m-1\simeq x-1$ である。
:::
::::

::::exercise{level=advanced}
$B^2/(2\mu_0)$ の一様な模型で場が半分になり体積が一定なら、放出されるエネルギーの初期値に対する比を求めよ。

:::solution
磁気エネルギーは二乗に比例するので、残りは四分の一である。差は四分の三になる。
:::
::::

::::exercise{level=advanced}
上流速さ $u_1=4u_2$ の Fermi 模型で $P_{\rm esc}/g$ を求めよ。

:::solution
比は $(4u_2/v)/[4(u_1-u_2)/(3v)]=3u_2/(u_1-u_2)=1$ である。従って積分スペクトルは $p^{-1}$ になる。
:::
::::

::::exercise{level=advanced}
剛体回転の局所模型で MRI の成長根がないことを示せ。

:::solution
$q_\Omega=0$、$\kappa^2=4\Omega^2$ である。$\gamma^2$ の判別式は $(4\Omega^2+2a^2)^2-4a^4>0$、根の和は負で積は $a^4>0$ なので $a\ne0$ で両根が負である。$a=0$ は零根を許す。
:::
::::

## 確認テスト

::::quiz
等温 Parker 風が正則に音速を通るために必要な条件はどれか。

- [x] 右辺の圧力と重力の項が臨界点で相殺する。
- [ ] 中心天体の質量がゼロである。
- [ ] 流出密度が臨界点で厳密にゼロになる。

:::explanation
式 :eqref[eq-plasma-parker] で $u=c_s$ では左辺の係数がゼロになる。有限な微分を持つには $2c_s^2/r=GM/r^2$ が必要である。
:::
::::

::::quiz{answer="2" tol="1e-6"}
$r=4$ の衝撃波で相対論的粒子の $\dd N/\dd E\propto E^{-s}$ の $s$ を答えよ。

:::explanation
$f\propto p^{-4}$ に相空間の $p^2$ を掛け、$\dd p/\dd E$ が一定なので指数は $-2$ となる。
:::
::::

::::quiz
磁気回転不安定性で弱い鉛直磁場はどんな役割を持つか。

- [x] 異なる半径の流体要素を結びつけて角運動量を交換させる。
- [ ] 回転角速度の勾配を常にゼロにする。
- [ ] 密度を減らして重力を消す。

:::explanation
磁気張力が半径方向にずれた要素を結び、角速度の差によって交換を続ける。式 :eqref[eq-plasma-mri] では $a\ne0$ と $q_\Omega>0$ が成長範囲を与える。
:::
::::

::::quiz{answer="6" tol="1e-6"}
$\Omega=8\,\mathrm{s^{-1}}$ の Kepler 円盤の局所最大成長率を SI で答えよ。

:::explanation
最大値は $3\Omega/4$ なので $6\,\mathrm{s^{-1}}$ である。
:::
::::

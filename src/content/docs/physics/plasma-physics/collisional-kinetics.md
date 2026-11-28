---
title: 衝突の運動論
description: 小角衝突から速度空間の拡散を導き、電気伝導と輸送の閉包を求める。
notation:
  - symbol: 'A_i'
    meaning: 速度変化の一次係数
  - symbol: 'B_{ij}'
    meaning: 速度変化の二次係数
  - symbol: 'H,G'
    meaning: Rosenbluth ポテンシャル
  - symbol: '\nu_D'
    meaning: 偏向周波数
---

## この章の目標

- 速度の跳びを展開して Fokker–Planck 方程式を得る。
- Landau 演算子の保存と散逸を示す。
- Lorentz 模型から電気伝導度を計算する。

## Fokker–Planck 方程式の導出

短時間に多くの微小な偏向があるとき、速度の飛びを一回ずつ追う代わりに、その条件付き平均と分散を求める。[Fokker–Planck 方程式](/physics/nonequilibrium/fokker-planck/) の導出を速度空間に適用する。$W(\vb v,\vb a)$ を $\vb v$ から $\vb v+\vb a$ への単位時間・単位速度体積当たりの遷移率とする。マスター方程式は
$$
\partial_t f(\vb v)=\int[W(\vb v-\vb a,\vb a)f(\vb v-\vb a)-W(\vb v,\vb a)f(\vb v)]\dd^3a .
$$
利得項を $\vb a$ で二次まで Taylor 展開すると零次は損失と消え、一次と二次から
$$
\partial_t f=-\partial_i(A_i f)+\tfrac12\partial_i\partial_j(B_{ij} f),\qquad
A_i=\int a_i W\dd^3a,\quad B_{ij}=\int a_i a_j W\dd^3a
\label{eq-plasma-fp}
$$
となる。これを :term[Fokker–Planck 方程式]{en="Fokker–Planck equation" yomi="ふぉっかーぷらんくほうていしき"} という。$A$ は $\mathrm{m\,s^{-2}}$、$B$ は $\mathrm{m^2\,s^{-3}}$ である。$h_iB_{ij}h_j=\int(\vb h\cdot\vb a)^2W\dd^3a\ge0$ だから $B$ は半正定値である。有限の跳びの過程では高次の項も残り、小角が対数的に支配する範囲でこの切断を使う。

:::example[一次元での加熱]
$A=0,B=2D_v$ が一定なら $\partial_t f=D_v f''$ である。境界項を落として二度部分積分すると $\dv{\langle v^2\rangle}{t}=2D_v$ であり、速度分布が広がる。
:::

:::example[定常 Gaussian]
一次元で $A=-\nu v,\ B=2\nu T/m$ とする。無流束定常状態では $A f-\frac12\partial_v(B f)=0$ なので $f'/f=-mv/T$ である。積分して $f\propto e^{-mv^2/(2T)}$ を得る。この模型の分散は $T/m$ になり、抗力と拡散の両方が必要である。
:::

### 線形摩擦と速度拡散の初期値解

一定係数の一次元方程式 $\partial_t f=\nu\partial_v(vf)+D\partial_v^2 f$ を考える。確率過程 $\dd v=-\nu v\dd t+\sqrt{2D}\dd W$ に積分因子を掛けて $v(t)=e^{-\nu t}v(0)+\sqrt{2D}\int_0^te^{-\nu(t-s)}\dd W_s$ を得る。最後の項は時間ごとの独立な Gaussian 増分の線形結合だから、分散は $2D\int_0^te^{-2\nu(t-s)}\dd s=D(1-e^{-2\nu t})/\nu$ である。初期値と独立なので分散を加算できる。$\nu>0$ での定常分散 $D/\nu$ は拡散と摩擦のつりあいで決まる。$\nu=0$ では同じ式の極限は $2Dt$ である。速度空間の確率流には :eqref[eq-neq-fp-current] を用いる。その式の拡散行列は本章の二次係数の半分 $D_{ij}=B_{ij}/2$ に当たる。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $x$ | 時間 $\nu t$ | 1 |
| $a$ | 初期平均を定常熱速さ $\sqrt{D/\nu}$ で割った値 | 1 |
| $s$ | 初期分散を $D/\nu$ で割った値 | 1 |

規格化した平均速度を $M$ とする。摩擦で指数的に減る。

$$
M=a e^{- x}
$$

初期分散が残る量を $V_{\rm ini}$ とする。平均と異なり減衰指数は二倍である。

$$
V_{\rm ini}=s e^{- 2 x}
$$

拡散が生んだ分散を $V_{\rm kick}$ とする。増分の分散を積分する。

$$
V_{\rm kick}=1 - e^{- 2 x}
$$

全分散を $V$ とする。初期分布と Wiener 増分を独立とする。

$$
V=V_{\rm ini} + V_{\rm kick}
$$

規格化した運動エネルギーを $K$ とする。単位は $mD/\nu$ である。

$$
K=\frac{M^{2}}{2} + \frac{V}{2}
$$

平衡からの二乗速度超過を $2K-1$ とする。負なら初期状態から加熱される側である。

$$
2K-1=M^{2} + V - 1
$$

:::example[単一速度の冷粒子群]
初期速度も分散もゼロの拡散加熱を求める。

入力値を次のように置く。

$$
\begin{aligned}
x&=1 \\a&=0 \\s&=0
\end{aligned}
$$

これらを順に代入する。

規格化した平均速度は

$$
M=\left(0\right) e^{- \left(1\right)}=0
$$

である。

初期分散が残る量は

$$
V_{\rm ini}=\left(0\right) e^{- 2 \left(1\right)}=0
$$

である。

拡散が生んだ分散は

$$
V_{\rm kick}=1 - e^{- 2 \left(1\right)}\simeq 0.86466
$$

である。

全分散は

$$
V=\left(0\right) + \left(1 - e^{-2}\right)\simeq 0.86466
$$

である。

規格化した運動エネルギーは

$$
K=\frac{\left(0\right)^{2}}{2} + \frac{\left(1 - e^{-2}\right)}{2}\simeq 0.43233
$$

である。

平衡からの二乗速度超過は

$$
2K-1=\left(0\right)^{2} + \left(1 - e^{-2}\right) - 1\simeq -0.13534
$$

である。

平均はゼロのまま、熱運動だけが増える。
:::

:::example[流れの緩和]
平衡分散を持つ初期分布の平均流を調べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=1 \\a&=2 \\s&=1
\end{aligned}
$$

これらを順に代入する。

規格化した平均速度は

$$
M=\left(2\right) e^{- \left(1\right)}\simeq 0.73576
$$

である。

初期分散が残る量は

$$
V_{\rm ini}=\left(1\right) e^{- 2 \left(1\right)}\simeq 0.13534
$$

である。

拡散が生んだ分散は

$$
V_{\rm kick}=1 - e^{- 2 \left(1\right)}\simeq 0.86466
$$

である。

全分散は

$$
V=\left(1 - e^{-2}\right) + \left(e^{-2}\right)=1
$$

である。

規格化した運動エネルギーは

$$
K=\frac{\left(1\right)}{2} + \frac{\left(\frac{2}{e}\right)^{2}}{2}\simeq 0.77067
$$

である。

平衡からの二乗速度超過は

$$
2K-1=\left(1\right) + \left(\frac{2}{e}\right)^{2} - 1\simeq 0.54134
$$

である。

分散は変わらず、流れのエネルギーが散逸する。
:::

:::example[高温の背景摂動]
平衡より大きい分散からの冷却を求める。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{2} \\a&=0 \\s&=4
\end{aligned}
$$

これらを順に代入する。

規格化した平均速度は

$$
M=\left(0\right) e^{- \left(\frac{1}{2}\right)}=0
$$

である。

初期分散が残る量は

$$
V_{\rm ini}=\left(4\right) e^{- 2 \left(\frac{1}{2}\right)}\simeq 1.4715
$$

である。

拡散が生んだ分散は

$$
V_{\rm kick}=1 - e^{- 2 \left(\frac{1}{2}\right)}\simeq 0.63212
$$

である。

全分散は

$$
V=\left(1 - e^{-1}\right) + \left(\frac{4}{e}\right)\simeq 2.1036
$$

である。

規格化した運動エネルギーは

$$
K=\frac{\left(0\right)^{2}}{2} + \frac{\left(1 + \frac{3}{e}\right)}{2}\simeq 1.0518
$$

である。

平衡からの二乗速度超過は

$$
2K-1=\left(0\right)^{2} + \left(1 + \frac{3}{e}\right) - 1\simeq 1.1036
$$

である。

初期分散の過剰は平均流より速い指数で減る。
:::

:::example[近い平均と冷たい分散]
同じ平均二乗速度を持つ非平衡初期条件を追う。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{10} \\a&=1 \\s&=0
\end{aligned}
$$

これらを順に代入する。

規格化した平均速度は

$$
M=\left(1\right) e^{- \left(\frac{1}{10}\right)}\simeq 0.90484
$$

である。

初期分散が残る量は

$$
V_{\rm ini}=\left(0\right) e^{- 2 \left(\frac{1}{10}\right)}=0
$$

である。

拡散が生んだ分散は

$$
V_{\rm kick}=1 - e^{- 2 \left(\frac{1}{10}\right)}\simeq 0.18127
$$

である。

全分散は

$$
V=\left(0\right) + \left(1 - e^{- \frac{1}{5}}\right)\simeq 0.18127
$$

である。

規格化した運動エネルギーは

$$
K=\frac{\left(1 - e^{- \frac{1}{5}}\right)}{2} + \frac{\left(\left(e^{- \frac{1}{10}}\right)\right)^{2}}{2}=\frac{1}{2}
$$

である。

平衡からの二乗速度超過は

$$
2K-1=\left(1 - e^{- \frac{1}{5}}\right) + \left(\left(e^{- \frac{1}{10}}\right)\right)^{2} - 1=0
$$

である。

二乗速度全体は一定でも、流れからランダムな速度へ分配が変わる。
:::


## Landau 衝突演算子

相対速度 $\vb u=\vb v-\vb v'$、$\kappa_{ab}=q_aq_b/(4\pi\varepsilon_0)$、$\Gamma_{ab}=2\pi\kappa_{ab}^2\ln\Lambda/m_a^2$ とする。小角散乱では試験粒子の垂直速度変化の大きさは $(m_r/m_a)u\chi\simeq2\abs{\kappa_{ab}}/(m_a u b)$ である。方位角を平均し、垂直二成分に二乗の半分ずつを割り当てる。散乱率 $u\,2\pi b\dd b\, f_b(\vb v')\dd^3v'$ で積分すると
$$
B_{ij}=2\Gamma_{ab}\int U_{ij}(\vb u) f_b(\vb v')\dd^3v',\qquad
U_{ij}(\vb u)=\frac{\delta_{ij}-u_i u_j/u^2}{u}
$$
である。平行の跳びは $-(m_r/m_a)u(1-\cos\chi)$ であり、同じ断面積積分で $A_i=-2\Gamma_{ab}(1+m_a/m_b)\int u_i u^{-3}f_b\dd^3v'$ になる。ここでは対数の速度依存は無視した。$\partial_{u_j}U_{ij}=-2u_i/u^3$ を用いて式 :eqref[eq-plasma-fp] を流束形に戻し、背景速度で部分積分すれば
$$
C_{ab}=\Gamma_{ab}\partial_{v_i}
\int U_{ij}\left[f_b'\partial_{v_j}f_a-\frac{m_a}{m_b}f_a\partial_{v'_j}f_b'\right]\dd^3v' .
\label{eq-plasma-landau-collision}
$$
$f_b'=f_b(\vb v')$ である。これが :term[Landau 衝突演算子]{en="Landau collision operator" yomi="らんだうしょうとつえんざんし"} である。

### 保存量とエントロピー

速度の境界項を消せば $\int C_{ab}\dd^3v=0$ である。括弧を $m_a f_a f_b'\vb h$、$\vb h=\nabla_v\ln f_a/m_a-\nabla_{v'}\ln f_b'/m_b$ と書く。$a,b$ と $\vb v,\vb v'$ を交換すると $\vb h$ の符号が逆で、$\ln\Lambda$ を粒子対に対して対称に取れば $m_a^2\Gamma_{ab}=m_b^2\Gamma_{ba}$ である。両種を足すと運動量の積分は相殺する。エネルギー重みでは残りが $\vb u\cdot\mathsf U\vb h$ となるが $\mathsf U\vb u=0$ なので消える。$H=\sum_a\int f_a\ln(f_a/f_{a*})\dd^3v$（$f_{a*}>0$ は同じ単位の定数） の時間微分も同様に交換し、
$$
\dot H=-\frac12\sum_{a,b}m_a^2\Gamma_{ab}\int f_a f_b'\,\vb h\cdot\mathsf U\vb h\,\dd^3v\dd^3v'\le0 .
$$
ここで $\mathsf U$ は $\vb u$ に垂直な射影を $u$ で割ったものである。共通の温度・平均速度を持つ Maxwell 分布では $\vb h=-\vb u/T$ となるので各衝突項がゼロである。

:::example[種ごとのエネルギー]
温度の違う二種では $\vb h=-\vb v/T_a+\vb v'/T_b$ が一般に $\vb u$ に平行でない。種別のエネルギーは交換され、両種の和だけが保存される。
:::

### エントロピー生成がゼロになる条件

滑らかで正の分布を仮定する。単一種で $H$ の減少がゼロなら、任意の $\vb v,\vb v'$ で $\nabla\ln f(\vb v)-\nabla\ln f(\vb v')$ が $\vb v-\vb v'$ と平行になる。$\vb v'$ を $\vb v$ に近づけると $\nabla\nabla\ln f$ は任意のベクトルをこれと平行なベクトルへ写す。座標軸と二つの軸の和に適用すると、この行列は $a(\vb v)I$ である。異なる添字 $i,j$ で $\partial_j a=\partial_j\partial_i^2\ln f=\partial_i(\partial_i\partial_j\ln f)=0$ なので $a$ は一定である。二度積分すると $\ln f=a v^2/2+\vb b\cdot\vb v+c$ となり、可積分性から $a<0$ である。これは Maxwell 分布である。多種の場合は各種の結果を種間のゼロ条件へ戻すと、対称な二次項の係数を質量で割ったものと、温度に対する平均速度の比が共通でなければならない。従って温度と平均速度がそろう。

:::example[同じ温度で異なる平均流]
$T_a=T_b=T$ にしても平均速度が別々なら $\vb h=-(\vb v-\vb v')/T+(\vb u_a-\vb u_b)/T$ である。後半の定ベクトルは全ての相対速度に平行にはできない。従って温度が等しくても種間摩擦によるエントロピー生成が残る。
:::

### Newton 核の射影と Landau 散逸

Landau 核は $\mathsf U(\vb w)=(I-\hat{\vb w}\hat{\vb w})/|\vb w|$ である。任意の二体の対数勾配差 $\vb a$ に対し $\vb a\cdot\mathsf U\vb a=(|\vb a|^2-(\hat{\vb w}\cdot\vb a)^2)/|\vb w|\ge0$ で、等号は $\vb a\parallel\vb w$ のときに限る。これは速度交換が相対速度の接平面にあるためである。衝突の弱形式では $f_s f_t$ を掛けて対称化することで、この非負量を積分したエントロピー生成を得る。以下の計算は一つの対の値であり、衝突率は種の密度と速度分布を積分して決まる。$w=0$ をこの点ごとの式へ代入しない。三次元で滑らかな分布なら $1/|\vb w|$ の局所速度積分は有限である。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $w$ | 規格化した相対速さ | 1 |
| $a$ | $\vb a$ の相対速度に平行な成分 | 1 |
| $b$ | 同じ勾配差の横成分の大きさ | 1 |
| $F$ | 二粒子の分布の積の規格化値 | 1 |

勾配差の規格化二乗を $\abs{\vb a}^2$ とする。直交成分を加える。

$$
\abs{\vb a}^2=a^{2} + b^{2}
$$

散逸に寄与する二乗の割合を $R_\perp$ とする。$\vb a\ne0$ とする。

$$
R_\perp=\frac{b^{2}}{\abs{\vb a}^2}
$$

接平面の固有値を $U_\perp$ とする。平行固有値はゼロである。

$$
U_\perp=\frac{1}{w}
$$

二体の非負二次形式を $D_{\rm pair}$ とする。横の勾配差だけが残る。

$$
D_{\rm pair}=U_{\perp} b^{2}
$$

分布で重みを掛けた生成項を $\dot S_{\rm pair}$ とする。全体に共通の正係数は除く。

$$
\dot S_{\rm pair}=D_{\rm pair} F
$$

核を掛けた横の成分を $C_\perp$ とする。散逸の二乗とは区別する。

$$
C_\perp=U_{\perp} b
$$

:::example[勾配が相対速度と平行]
勾配が有限でも散逸が消えるか確かめる。

入力値を次のように置く。

$$
\begin{aligned}
w&=1 \\a&=2 \\b&=0 \\F&=1
\end{aligned}
$$

これらを順に代入する。

勾配差の規格化二乗は

$$
\abs{\vb a}^2=\left(0\right)^{2} + \left(2\right)^{2}=4
$$

である。

散逸に寄与する二乗の割合は

$$
R_\perp=\frac{\left(0\right)^{2}}{\left(4\right)}=0
$$

である。

接平面の固有値は

$$
U_\perp=\frac{1}{\left(1\right)}=1
$$

である。

二体の非負二次形式は

$$
D_{\rm pair}=\left(0\right)^{2} \left(1\right)=0
$$

である。

分布で重みを掛けた生成項は

$$
\dot S_{\rm pair}=\left(0\right) \left(1\right)=0
$$

である。

核を掛けた横の成分は

$$
C_\perp=\left(0\right) \left(1\right)=0
$$

である。

核の零空間にあるため生成項はゼロである。
:::

:::example[直交する勾配]
射影で失われる成分があるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
w&=1 \\a&=0 \\b&=2 \\F&=1
\end{aligned}
$$

これらを順に代入する。

勾配差の規格化二乗は

$$
\abs{\vb a}^2=\left(0\right)^{2} + \left(2\right)^{2}=4
$$

である。

散逸に寄与する二乗の割合は

$$
R_\perp=\frac{\left(2\right)^{2}}{\left(4\right)}=1
$$

である。

接平面の固有値は

$$
U_\perp=\frac{1}{\left(1\right)}=1
$$

である。

二体の非負二次形式は

$$
D_{\rm pair}=\left(1\right) \left(2\right)^{2}=4
$$

である。

分布で重みを掛けた生成項は

$$
\dot S_{\rm pair}=\left(1\right) \left(4\right)=4
$$

である。

核を掛けた横の成分は

$$
C_\perp=\left(1\right) \left(2\right)=2
$$

である。

勾配は全て接平面にある。
:::

:::example[相対速さが大きい対]
同じ横勾配で核の重みを比較する。

入力値を次のように置く。

$$
\begin{aligned}
w&=4 \\a&=1 \\b&=2 \\F&=1
\end{aligned}
$$

これらを順に代入する。

勾配差の規格化二乗は

$$
\abs{\vb a}^2=\left(1\right)^{2} + \left(2\right)^{2}=5
$$

である。

散逸に寄与する二乗の割合は

$$
R_\perp=\frac{\left(2\right)^{2}}{\left(5\right)}=\frac{4}{5}
$$

である。

接平面の固有値は

$$
U_\perp=\frac{1}{\left(4\right)}=\frac{1}{4}
$$

である。

二体の非負二次形式は

$$
D_{\rm pair}=\left(2\right)^{2} \left(\frac{1}{4}\right)=1
$$

である。

分布で重みを掛けた生成項は

$$
\dot S_{\rm pair}=\left(1\right)^{2}=1
$$

である。

核を掛けた横の成分は

$$
C_\perp=\left(2\right) \left(\frac{1}{4}\right)=\frac{1}{2}
$$

である。

相対速さの逆数だけ二次形式が小さくなる。
:::

:::example[速度分布の希薄な領域]
直交成分が残っても分布の積が小さい場合を調べる。

入力値を次のように置く。

$$
\begin{aligned}
w&=1 \\a&=2 \\b&=1 \\F&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

勾配差の規格化二乗は

$$
\abs{\vb a}^2=\left(1\right)^{2} + \left(2\right)^{2}=5
$$

である。

散逸に寄与する二乗の割合は

$$
R_\perp=\frac{\left(1\right)^{2}}{\left(5\right)}=\frac{1}{5}
$$

である。

接平面の固有値は

$$
U_\perp=\frac{1}{\left(1\right)}=1
$$

である。

二体の非負二次形式は

$$
D_{\rm pair}=\left(1\right)^{3}=1
$$

である。

分布で重みを掛けた生成項は

$$
\dot S_{\rm pair}=\left(1\right) \left(\frac{1}{100}\right)=\frac{1}{100}
$$

である。

核を掛けた横の成分は

$$
C_\perp=\left(1\right)^{2}=1
$$

である。

速い尾部での散逸には対の個数も必要である。
:::


## Rosenbluth ポテンシャル

多数の相対速度で同じ核を積分する代わりに、次の二つのスカラーを計算する。
$$
H(\vb v)=\int\frac{f_b(\vb v')}{\abs{\vb v-\vb v'}}\dd^3v',\qquad
G(\vb v)=\int\abs{\vb v-\vb v'} f_b(\vb v')\dd^3v'.
$$
これを :term[Rosenbluth ポテンシャル]{en="Rosenbluth potentials" yomi="ろーぜんぶるーすぽてんしゃる"} という。$U_{ij}=\partial_i\partial_j u$ なので $B_{ij}=2\Gamma_{ab}\partial_i\partial_j G$ である。また $\partial_i (1/u)=-u_i/u^3$ より $A_i=2\Gamma_{ab}(1+m_a/m_b)\partial_i H$ である。半径 $u>0$ で $\nabla^2u=2/u$、原点に追加のデルタがないのは半径 $\epsilon$ の球の勾配流束が $4\pi\epsilon^2$ で消えるためである。一方 $1/u$ の勾配流束は $-4\pi$ である。従って
$$
\nabla_v^2 G=2H,\qquad \nabla_v^2 H=-4\pi f_b.
$$

:::example[止まった背景]
$f_b=n_b\delta^{(3)}(\vb v')$ なら $H=n_b/v$、$G=n_b v$ である。拡散は速度球の接平面に限られ、有限質量の背景では反跳により平均エネルギーも変化する。
:::

:::example[一様な速度球]
背景分布を $v<V$ で $3n_b/(4\pi V^3)$、外でゼロとする。球対称の Poisson 方程式から内部で $H'=-n_b v/V^3$、外で $H'=-n_b/v^2$ である。$H\to0$ の外部条件と $v=V$ の連続性から内部は $H=(n_b/(2V))(3-v^2/V^2)$ となる。摩擦の一次係数は低速で速度に比例し、高速では逆二乗となる。
:::

### 球対称背景の縦・横拡散

球対称な背景では $G=G(v)$ であり、Hessian の固有値は動径方向で $G''$、二つの接線方向で $G'/v$ である。従って速度の二次係数は、接線方向の一成分につき $B_\perp=2\Gamma_{ab}G'/v$、動径方向で $B_\parallel=2\Gamma_{ab}G''$ となる。ここで平行・垂直は磁場ではなく試験粒子の速度を基準とする。$G$ が原点で正則なら $\nabla_v^2G=2H$ は
$$
G'(v)=\frac{2}{v^2}\int_0^v s^2H(s)\dd s,\qquad G''(v)=2H(v)-2G'(v)/v
$$
を与える。これは $\partial_v(v^2G')=2v^2H$ を積分した結果である。大きな速さでの偏向とエネルギー変化は、同じ拡散係数から得られるが、二つの成分の大きさは異なる。

### 等方 Maxwell 背景の Rosenbluth 摩擦

背景粒子の密度 $n$ の分布を $F=n(2\pi v_t^2)^{-3/2}e^{-u^2/(2v_t^2)}$ とする。$H(v)=\int F(\vb u)/|\vb v-\vb u|\dd^3u$ は $\nabla_v^2 H=-4\pi F(v)$ を満たす。球の中心で有限とすれば $H'=-4\pi v^{-2}\int_0^v u^2F(u)\dd u$ である。$x=v/(\sqrt2v_t)$ と置き、$u^2 e^{-u^2}$ の原始関数を部分積分すると内側粒子割合 $C(x)=\operatorname{erf}x-2xe^{-x^2}/\sqrt\pi$ が得られる。$H(\infty)=0$ から $H=n\,\operatorname{erf}x/v$ である。$C(x)\simeq4x^3/(3\sqrt\pi)$ なので中心の摩擦は線形、遠方では逆二乗に近い。実際の試験粒子の減速には場粒子と試験粒子の質量係数も掛ける。

同じ $x$ と $C(x)$ を用いると、$G'(v)=n[\operatorname{erf}x-C(x)/(2x^2)]$ である。これを確かめるには $2\int_0^v s\,\operatorname{erf}[s/(\sqrt2v_t)]\dd s$ を部分積分すればよい。$G''=2H-2G'/v$ より
$$
\frac{G'}v=\frac{n}{v}\left[\operatorname{erf}x-\frac{C(x)}{2x^2}\right],
\qquad
G''=\frac{n}{v}\frac{C(x)}{x^2}.
$$
$x\to0$ では両者とも $2n/(3v_t)\sqrt{2/\pi}$ に近づき、拡散は等方になる。$x\gg1$ なら $C\simeq1$ なので $B_\parallel/B_\perp\simeq x^{-2}$ である。例えば $v=10v_t$ なら先頭近似でこの比は $0.02$ であり、二つの接線方向の拡散が一つの動径方向より大きい。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $x$ | 速度 $v/(\sqrt2v_t)>0$ | 1 |

誤差関数を $G$ とする。Maxwell 積分の原始関数である。

$$
G=\operatorname{erf}{\left(x \right)}
$$

速度球に入る背景粒子の割合を $C(x)$ とする。半径方向の三次元の測度を含める。

$$
C(x)=G - \frac{2 x e^{- x^{2}}}{\sqrt{\pi}}
$$

規格化した Rosenbluth ポテンシャルを $H\sqrt2v_t/n$ とする。無限遠で消える解を選ぶ。

$$
H\sqrt2v_t/n=\frac{G}{x}
$$

摩擦に入る勾配の大きさを ${\cal F}$ とする。単位は $n/(2v_t^2)$ である。

$$
{\cal F}=\frac{C(x)}{x^{2}}
$$

低速近似の値を ${\cal F}_{\rm low}$ とする。原点で摩擦が消えることを表す。

$$
{\cal F}_{\rm low}=\frac{4 x}{3 \sqrt{\pi}}
$$

高速近似の値を ${\cal F}_{\rm high}$ とする。球内の割合を一にする。

$$
{\cal F}_{\rm high}=\frac{1}{x^{2}}
$$

:::example[非常に遅い粒子]
原点近くで逆二乗の近似が使えるか確認する。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{10}
\end{aligned}
$$

これらを順に代入する。

誤差関数は

$$
G=\operatorname{erf}{\left(\left(\frac{1}{10}\right) \right)}\simeq 0.11246
$$

である。

速度球に入る背景粒子の割合は

$$
C(x)=- \frac{2 \left(\frac{1}{10}\right) e^{- \left(\frac{1}{10}\right)^{2}}}{\sqrt{\pi}} + \left(\operatorname{erf}{\left(\frac{1}{10} \right)}\right)\simeq 0.00074776
$$

である。

規格化した Rosenbluth ポテンシャルは

$$
H\sqrt2v_t/n=\frac{\left(\operatorname{erf}{\left(\frac{1}{10} \right)}\right)}{\left(\frac{1}{10}\right)}\simeq 1.1246
$$

である。

摩擦に入る勾配の大きさは

$$
{\cal F}=\frac{\left(- \frac{1}{5 \sqrt{\pi} e^{\frac{1}{100}}} + \operatorname{erf}{\left(\frac{1}{10} \right)}\right)}{\left(\frac{1}{10}\right)^{2}}\simeq 0.074776
$$

である。

低速近似の値は

$$
{\cal F}_{\rm low}=\frac{4 \left(\frac{1}{10}\right)}{3 \sqrt{\pi}}\simeq 0.075225
$$

である。

高速近似の値は

$$
{\cal F}_{\rm high}=\frac{1}{\left(\frac{1}{10}\right)^{2}}=100
$$

である。

低速の線形近似の方が近い。
:::

:::example[熱速さの近く]
内側の粒子の割合と摩擦を求める。

入力値を次のように置く。

$$
\begin{aligned}
x&=1
\end{aligned}
$$

これらを順に代入する。

誤差関数は

$$
G=\operatorname{erf}{\left(\left(1\right) \right)}\simeq 0.8427
$$

である。

速度球に入る背景粒子の割合は

$$
C(x)=- \frac{2 \left(1\right) e^{- \left(1\right)^{2}}}{\sqrt{\pi}} + \left(\operatorname{erf}{\left(1 \right)}\right)\simeq 0.42759
$$

である。

規格化した Rosenbluth ポテンシャルは

$$
H\sqrt2v_t/n=\frac{\left(\operatorname{erf}{\left(1 \right)}\right)}{\left(1\right)}\simeq 0.8427
$$

である。

摩擦に入る勾配の大きさは

$$
{\cal F}=\frac{\left(- \frac{2}{e \sqrt{\pi}} + \operatorname{erf}{\left(1 \right)}\right)}{\left(1\right)^{2}}\simeq 0.42759
$$

である。

低速近似の値は

$$
{\cal F}_{\rm low}=\frac{4 \left(1\right)}{3 \sqrt{\pi}}\simeq 0.75225
$$

である。

高速近似の値は

$$
{\cal F}_{\rm high}=\frac{1}{\left(1\right)^{2}}=1
$$

である。

線形・逆二乗のどちらの先頭項も十分な精度とは限らない。
:::

:::example[高速の粒子]
高速近似が近づく理由を調べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=3
\end{aligned}
$$

これらを順に代入する。

誤差関数は

$$
G=\operatorname{erf}{\left(\left(3\right) \right)}\simeq 0.99998
$$

である。

速度球に入る背景粒子の割合は

$$
C(x)=- \frac{2 \left(3\right) e^{- \left(3\right)^{2}}}{\sqrt{\pi}} + \left(\operatorname{erf}{\left(3 \right)}\right)\simeq 0.99956
$$

である。

規格化した Rosenbluth ポテンシャルは

$$
H\sqrt2v_t/n=\frac{\left(\operatorname{erf}{\left(3 \right)}\right)}{\left(3\right)}\simeq 0.33333
$$

である。

摩擦に入る勾配の大きさは

$$
{\cal F}=\frac{\left(- \frac{6}{\sqrt{\pi} e^{9}} + \operatorname{erf}{\left(3 \right)}\right)}{\left(3\right)^{2}}\simeq 0.11106
$$

である。

低速近似の値は

$$
{\cal F}_{\rm low}=\frac{4 \left(3\right)}{3 \sqrt{\pi}}\simeq 2.2568
$$

である。

高速近似の値は

$$
{\cal F}_{\rm high}=\frac{1}{\left(3\right)^{2}}=\frac{1}{9}
$$

である。

球内が背景密度のほぼ全てを含む。
:::

:::example[さらに低速の粒子]
速度を十分ゼロへ近づけたときの符号を確かめる。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

誤差関数は

$$
G=\operatorname{erf}{\left(\left(\frac{1}{100}\right) \right)}\simeq 0.011283
$$

である。

速度球に入る背景粒子の割合は

$$
C(x)=- \frac{2 \left(\frac{1}{100}\right) e^{- \left(\frac{1}{100}\right)^{2}}}{\sqrt{\pi}} + \left(\operatorname{erf}{\left(\frac{1}{100} \right)}\right)\simeq 7.5221\times 10^{-7}
$$

である。

規格化した Rosenbluth ポテンシャルは

$$
H\sqrt2v_t/n=\frac{\left(\operatorname{erf}{\left(\frac{1}{100} \right)}\right)}{\left(\frac{1}{100}\right)}\simeq 1.1283
$$

である。

摩擦に入る勾配の大きさは

$$
{\cal F}=\frac{\left(- \frac{1}{50 \sqrt{\pi} e^{\frac{1}{10000}}} + \operatorname{erf}{\left(\frac{1}{100} \right)}\right)}{\left(\frac{1}{100}\right)^{2}}\simeq 0.0075221
$$

である。

低速近似の値は

$$
{\cal F}_{\rm low}=\frac{4 \left(\frac{1}{100}\right)}{3 \sqrt{\pi}}\simeq 0.0075225
$$

である。

高速近似の値は

$$
{\cal F}_{\rm high}=\frac{1}{\left(\frac{1}{100}\right)^{2}}=10000
$$

である。

摩擦の大きさは有限で、中心でゼロになる。
:::


## Lorentz 気体模型と電気伝導度

重いイオンを静止散乱体とし、電子同士の衝突を無視する。この :term[Lorentz 気体模型]{en="Lorentz gas model" yomi="ろーれんつきたいもけい"} では速さが保存され、方向だけが拡散する。$m_b\to\infty$ の :eqref[eq-plasma-landau-collision] で $U$ が接平面への射影であることを使えば、球面 Laplacian $\Delta_{\mathbb S^2}$ により $C=(\nu_D(v)/2)\Delta_{\mathbb S^2}f$、$\nu_D=2\Gamma_{ei}n_i/v^3$ となる。軸対称で $\xi=v_z/v$ なら $\Delta_{\mathbb S^2}=\partial_\xi[(1-\xi^2)\partial_\xi]$ であり、$\Delta_{\mathbb S^2}\xi=-2\xi$ である。

### Maxwell 分布の線形伝導度

均一な弱い電場 $E\hat{\vb z}$ に対して $(qE/m)\partial_{v_z}f_0=C[f_1]$ とする。Maxwell 平衡の微分は $-mv_z f_0/T$ なので $f_1=qEv_z f_0/(T\nu_D)$ である。電流積分で $\langle v_z^2\rangle_{\rm angle}=v^2/3$ を使うと
$$
\sigma=\frac{n_e e^2}{3T_e}\left\langle\frac{v^2}{\nu_D(v)}\right\rangle.
\label{eq-plasma-lorentz-conductivity}
$$
$v_t=\sqrt{T_e/m_e}$、$\nu_*=\nu_D(v_t)$ とすれば、Maxwell の速さの分布で $\langle v^5\rangle=48\sqrt{2/\pi}v_t^5$ なので $\sigma=16\sqrt{2/\pi}\, n_e e^2/(m_e\nu_*)$ である。ここでの $\nu_*$ は熱速度 $v_t$ における偏向周波数である。一方、:cite[goldston-rutherford] で用いられる熱平均時間の慣例の一つは $\tau_e=3(4\pi\varepsilon_0)^2\sqrt{m_e}T_e^{3/2}/(4\sqrt{2\pi}n_i Z^2e^4\ln\Lambda)$ である。係数を比較するには周波数の定義をそろえる。実際の水素の伝導度には電子間衝突を含める必要がある。

:::example[電流を担う速度]
式 :eqref[eq-plasma-lorentz-conductivity] の Maxwell 速さ積分の重みは $v^7e^{-v^2/(2v_t^2)}$ である。対数微分をゼロにすると最大は $v=\sqrt7 v_t$ であり、代表速度一つで Drude 係数を決めると差が出る。
:::

:::example[角度モードの緩和]
$\xi^2-1/3$ に角度演算子 $\partial_\xi((1-\xi^2)\partial_\xi)$ を作用させると $2-6\xi^2=-6(\xi^2-1/3)$ である。従って二次の異方成分は $3\nu_D$ で減衰する。一次の流れの成分は $\nu_D$ であり、角度緩和でもモーメントごとに時間が違う。
:::

### Lorentz 気体の角モード緩和

等方固定イオンによる角拡散を $C[f]=(\nu_D(v)/2)\partial_\xi[(1-\xi^2)\partial_\xi f]$ とする。$\xi=v_\parallel/v$ で軸対称な分布を考える。Legendre 多項式 $P_\ell$ の微分方程式は、Rodrigues 公式を微分して $\partial_\xi[(1-\xi^2)P_\ell']=-\ell(\ell+1)P_\ell$ を与える。従って $\ell=0$ は保存し、電流を運ぶ $\ell=1$ は $e^{-\nu_Dt}$、圧力異方性の $\ell=2$ は $e^{-3\nu_Dt}$ で減る。Coulomb 角拡散で $\nu_D\propto v^{-3}$ とすれば全速度域の電流が一つの率で減ることはない。以下では一つの速度殻で角分布を追う。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $x$ | 速度と基準速度の比 | 1 |
| $s$ | 基準散乱率を掛けた時間 | 1 |
| $a_1$ | 初期の $\ell=1$ 係数 | 1 |
| $a_2$ | 初期の $\ell=2$ 係数 | 1 |

角拡散率の比を $\nu_D/\nu_{D0}$ とする。固定イオンの Coulomb 散乱を仮定する。

$$
\nu_D/\nu_{D0}=\frac{1}{x^{3}}
$$

選んだ速度での時間を $u$ とする。該当する速度殻の緩和率で測る。

$$
u=\nu_D/\nu_{D0} s
$$

残る一次の角モードを $a_1(t)$ とする。電流の速度ごとの寄与を運ぶ。

$$
a_1(t)=a_{1} e^{- u}
$$

残る二次の角モードを $a_2(t)$ とする。圧力異方性に寄与する。

$$
a_2(t)=a_{2} e^{- 3 u}
$$

軸正方向の分布を $f(1)/f_0$ とする。$P_1(1)=P_2(1)=1$ とする。

$$
f(1)/f_0=a_1(t) + a_2(t) + 1
$$

軸に垂直な分布を $f(0)/f_0$ とする。$P_1(0)=0$、$P_2(0)=-1/2$ である。

$$
f(0)/f_0=1 - \frac{a_2(t)}{2}
$$

:::example[基準速度での緩和]
電流と圧力異方性の残り方を比較する。

入力値を次のように置く。

$$
\begin{aligned}
x&=1 \\s&=1 \\a_1&=\frac{1}{5} \\a_2&=\frac{1}{5}
\end{aligned}
$$

これらを順に代入する。

角拡散率の比は

$$
\nu_D/\nu_{D0}=\frac{1}{\left(1\right)^{3}}=1
$$

である。

選んだ速度での時間は

$$
u=\left(1\right)^{2}=1
$$

である。

残る一次の角モードは

$$
a_1(t)=\left(\frac{1}{5}\right) e^{- \left(1\right)}\simeq 0.073576
$$

である。

残る二次の角モードは

$$
a_2(t)=\left(\frac{1}{5}\right) e^{- 3 \left(1\right)}\simeq 0.0099574
$$

である。

軸正方向の分布は

$$
f(1)/f_0=\left(\frac{1}{5 e^{3}}\right) + \left(\frac{1}{5 e}\right) + 1\simeq 1.0835
$$

である。

軸に垂直な分布は

$$
f(0)/f_0=1 - \frac{\left(\frac{1}{5 e^{3}}\right)}{2}\simeq 0.99502
$$

である。

二次のモードの方が大きく減衰する。
:::

:::example[高速殻の非平衡]
同じ時間でも異方性が残るか調べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=3 \\s&=1 \\a_1&=\frac{1}{5} \\a_2&=\frac{1}{5}
\end{aligned}
$$

これらを順に代入する。

角拡散率の比は

$$
\nu_D/\nu_{D0}=\frac{1}{\left(3\right)^{3}}=\frac{1}{27}
$$

である。

選んだ速度での時間は

$$
u=\left(1\right) \left(\frac{1}{27}\right)=\frac{1}{27}
$$

である。

残る一次の角モードは

$$
a_1(t)=\left(\frac{1}{5}\right) e^{- \left(\frac{1}{27}\right)}\simeq 0.19273
$$

である。

残る二次の角モードは

$$
a_2(t)=\left(\frac{1}{5}\right) e^{- 3 \left(\frac{1}{27}\right)}\simeq 0.17897
$$

である。

軸正方向の分布は

$$
f(1)/f_0=\left(\frac{1}{5 e^{\frac{1}{27}}}\right) + \left(\frac{1}{5 e^{\frac{1}{9}}}\right) + 1\simeq 1.3717
$$

である。

軸に垂直な分布は

$$
f(0)/f_0=1 - \frac{\left(\frac{1}{5 e^{\frac{1}{9}}}\right)}{2}\simeq 0.91052
$$

である。

高速では角拡散が弱く、初期分布に近い。
:::

:::example[低速殻での等方化]
低速の分布の残りを計算する。

入力値を次のように置く。

$$
\begin{aligned}
x&=\frac{1}{2} \\s&=1 \\a_1&=\frac{1}{5} \\a_2&=\frac{1}{5}
\end{aligned}
$$

これらを順に代入する。

角拡散率の比は

$$
\nu_D/\nu_{D0}=\frac{1}{\left(\frac{1}{2}\right)^{3}}=8
$$

である。

選んだ速度での時間は

$$
u=\left(1\right) \left(8\right)=8
$$

である。

残る一次の角モードは

$$
a_1(t)=\left(\frac{1}{5}\right) e^{- \left(8\right)}\simeq 6.7093\times 10^{-5}
$$

である。

残る二次の角モードは

$$
a_2(t)=\left(\frac{1}{5}\right) e^{- 3 \left(8\right)}\simeq 7.5503\times 10^{-12}
$$

である。

軸正方向の分布は

$$
f(1)/f_0=\left(\frac{1}{5 e^{24}}\right) + \left(\frac{1}{5 e^{8}}\right) + 1\simeq 1.0001
$$

である。

軸に垂直な分布は

$$
f(0)/f_0=1 - \frac{\left(\frac{1}{5 e^{24}}\right)}{2}\simeq 1
$$

である。

電流の一次成分もほとんど消える。
:::

:::example[圧力のみに非等方性]
電流がゼロでも緩和が存在するか調べる。

入力値を次のように置く。

$$
\begin{aligned}
x&=1 \\s&=\frac{1}{3} \\a_1&=0 \\a_2&=\frac{1}{5}
\end{aligned}
$$

これらを順に代入する。

角拡散率の比は

$$
\nu_D/\nu_{D0}=\frac{1}{\left(1\right)^{3}}=1
$$

である。

選んだ速度での時間は

$$
u=\left(1\right) \left(\frac{1}{3}\right)=\frac{1}{3}
$$

である。

残る一次の角モードは

$$
a_1(t)=\left(0\right) e^{- \left(\frac{1}{3}\right)}=0
$$

である。

残る二次の角モードは

$$
a_2(t)=\left(\frac{1}{5}\right) e^{- 3 \left(\frac{1}{3}\right)}\simeq 0.073576
$$

である。

軸正方向の分布は

$$
f(1)/f_0=\left(0\right) + \left(\frac{1}{5 e}\right) + 1\simeq 1.0736
$$

である。

軸に垂直な分布は

$$
f(0)/f_0=1 - \frac{\left(\frac{1}{5 e}\right)}{2}\simeq 0.96321
$$

である。

二次のモードだけが指数で減り、一次は励起されない。
:::


### Lorentz 模型の熱電輸送を同じ積分で求める

静止イオンの中の電子を考え、分布の等方部分を局所 Maxwell 分布 $F_0(n,T)$ とする。平行方向に弱い勾配を与え、角の一次だけを残す。$v_z\partial_z F_0+(qE/m)\partial_{v_z}F_0=-\nu_D f_1$ より
$$
f_1=\frac{v_z F_0}{\nu_D}
\left[X+Y\left(\frac{mv^2}{2T}-\frac32\right)\right],
\qquad
X=\frac{qE}{T}-\partial_z\ln n,\quad Y=-\partial_z\ln T.
$$
$L_0=\int v_z^2 F_0/\nu_D\,\dd^3v$ とおき、$\nu_D\propto v^{-3}$ を使うと速さ積分の重みは $v^7\exp[-mv^2/(2T)]$ となる。$\xi=mv^2/(2T)$ では重みが $\xi^3 e^{-\xi}$ で、規格化した平均 $\langle\xi\rangle_w=4,\ \langle\xi^2\rangle_w=20$ である。これは $\int_0^\infty \xi^j e^{-\xi}\dd\xi=j!$ を部分積分で求めればよい。従って数流束と、運ばれるエネルギーの流束は
$$
\Gamma=L_0(X+5Y/2),\qquad
J_E=T L_0(4X+14Y).
\label{eq-plasma-lorentz-thermoelectric}
$$
$J_E$ は中心化した熱流と区別する。平均速度 $u=\Gamma/n$ に変換すると、一次で $5nTu/2=5T\Gamma/2$ を引くことになる。これは運動エネルギーの移流と圧力仕事である。数流束をゼロにする電場を選ぶと $X=-5Y/2$、熱流は $J_E=4T L_0Y=-4L_0\,\partial_z T$ になる。:eqref[eq-plasma-lorentz-thermoelectric] の係数は固定散乱体のもので、電子間衝突込みの Braginskii 数係数にそのまま使わない。

:::remark[熱流の温度の規約]
上式の $T$ は粒子一個当たりのエネルギーである。ケルビン単位の温度を $T_{\rm K}$ とする場合、$-\partial_z T=-k_B\partial_z T_{\rm K}$ になるので、熱伝導係数を比較する前にこの $k_B$ を含める。
:::

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $X$ | 電気力と密度勾配を合わせた駆動の規格化値 | 1 |
| $Y$ | 温度勾配の駆動の規格化値 | 1 |

規格化数流束を $\Gamma/L_0$ とする。角の一次の速度積分を取る。

$$
\Gamma/L_0=X + \frac{5 Y}{2}
$$

規格化した運動エネルギー流束を $J_E/(TL_0)$ とする。重みにさらに $\xi$ を掛ける。

$$
J_E/(TL_0)=4 X + 14 Y
$$

流れによって運ばれるエンタルピーの規格化値を $5\Gamma/(2L_0)$ とする。平均流への変換で引く項である。

$$
5\Gamma/(2L_0)=\frac{5 \Gamma/L_{0}}{2}
$$

規格化した中心熱流を $q_z/(TL_0)$ とする。一次の中心速度三次モーメントである。

$$
q_z/(TL_0)=\frac{J_E}{TL_0}-\frac{5\Gamma}{2L_0}
$$

数流束を打ち消す駆動を $X_{\Gamma=0}$ とする。電気力で必要な補償を与える。

$$
X_{\Gamma=0}=- \frac{5 Y}{2}
$$

無粒子流での熱流を $J_{E,\Gamma=0}/(TL_0)$ とする。流れを止めても熱は流れる。

$$
J_{E,\Gamma=0}/(TL_0)=4 X_{\Gamma=0} + 14 Y
$$

:::example[等温での電場応答]
数流と熱流の違いを求める。

入力値を次のように置く。

$$
\begin{aligned}
X&=1 \\Y&=0
\end{aligned}
$$

これらを順に代入する。

規格化数流束は

$$
\Gamma/L_0=\frac{5 \left(0\right)}{2} + \left(1\right)=1
$$

である。

規格化した運動エネルギー流束は

$$
J_E/(TL_0)=14 \left(0\right) + 4 \left(1\right)=4
$$

である。

流れによって運ばれるエンタルピーの規格化値は

$$
5\Gamma/(2L_0)=\frac{5 \left(1\right)}{2}=\frac{5}{2}
$$

である。

規格化した中心熱流は

$$
q_z/(TL_0)=\left(4\right) - \left(\frac{5}{2}\right)=\frac{3}{2}
$$

である。

数流束を打ち消す駆動は

$$
X_{\Gamma=0}=- \frac{5 \left(0\right)}{2}=0
$$

である。

無粒子流での熱流は

$$
J_{E,\Gamma=0}/(TL_0)=18 \left(0\right)=0
$$

である。

一次の非平衡補正は高速の粒子を強く重み付けするため、中心熱流も有限である。
:::

:::example[温度勾配のみ]
等密度で電場を補償しない場合を調べる。

入力値を次のように置く。

$$
\begin{aligned}
X&=0 \\Y&=1
\end{aligned}
$$

これらを順に代入する。

規格化数流束は

$$
\Gamma/L_0=\left(0\right) + \frac{5 \left(1\right)}{2}=\frac{5}{2}
$$

である。

規格化した運動エネルギー流束は

$$
J_E/(TL_0)=4 \left(0\right) + 14 \left(1\right)=14
$$

である。

流れによって運ばれるエンタルピーの規格化値は

$$
5\Gamma/(2L_0)=\frac{5 \left(\frac{5}{2}\right)}{2}=\frac{25}{4}
$$

である。

規格化した中心熱流は

$$
q_z/(TL_0)=\left(14\right) - \left(\frac{25}{4}\right)=\frac{31}{4}
$$

である。

数流束を打ち消す駆動は

$$
X_{\Gamma=0}=- \frac{5 \left(1\right)}{2}=- \frac{5}{2}
$$

である。

無粒子流での熱流は

$$
J_{E,\Gamma=0}/(TL_0)=4 \left(- \frac{5}{2}\right) + 14 \left(1\right)=4
$$

である。

粒子とエネルギーの流束がともに生じる。
:::

:::example[数流束を消す電場]
電場で補償しても熱が運ばれることを確かめる。

入力値を次のように置く。

$$
\begin{aligned}
X&=- \frac{5}{2} \\Y&=1
\end{aligned}
$$

これらを順に代入する。

規格化数流束は

$$
\Gamma/L_0=\left(- \frac{5}{2}\right) + \frac{5 \left(1\right)}{2}=0
$$

である。

規格化した運動エネルギー流束は

$$
J_E/(TL_0)=4 \left(- \frac{5}{2}\right) + 14 \left(1\right)=4
$$

である。

流れによって運ばれるエンタルピーの規格化値は

$$
5\Gamma/(2L_0)=\frac{5 \left(0\right)}{2}=0
$$

である。

規格化した中心熱流は

$$
q_z/(TL_0)=- \left(0\right) + \left(4\right)=4
$$

である。

数流束を打ち消す駆動は

$$
X_{\Gamma=0}=- \frac{5 \left(1\right)}{2}=- \frac{5}{2}
$$

である。

無粒子流での熱流は

$$
J_{E,\Gamma=0}/(TL_0)=4 \left(- \frac{5}{2}\right) + 14 \left(1\right)=4
$$

である。

中心熱流と全エネルギー流束が一致する。
:::

:::example[勾配を反転する]
無粒子流の条件を保ったまま熱流を反転する。

入力値を次のように置く。

$$
\begin{aligned}
X&=\frac{5}{2} \\Y&=-1
\end{aligned}
$$

これらを順に代入する。

規格化数流束は

$$
\Gamma/L_0=\frac{5 \left(-1\right)}{2} + \left(\frac{5}{2}\right)=0
$$

である。

規格化した運動エネルギー流束は

$$
J_E/(TL_0)=14 \left(-1\right) + 4 \left(\frac{5}{2}\right)=-4
$$

である。

流れによって運ばれるエンタルピーの規格化値は

$$
5\Gamma/(2L_0)=\frac{5 \left(0\right)}{2}=0
$$

である。

規格化した中心熱流は

$$
q_z/(TL_0)=\left(-4\right) - \left(0\right)=-4
$$

である。

数流束を打ち消す駆動は

$$
X_{\Gamma=0}=- \frac{5 \left(-1\right)}{2}=\frac{5}{2}
$$

である。

無粒子流での熱流は

$$
J_{E,\Gamma=0}/(TL_0)=14 \left(-1\right) + 4 \left(\frac{5}{2}\right)=-4
$$

である。

大きさは同じで向きだけが反転し、熱流は温度が下がる側を向く。
:::


## 逃走電子

小角衝突による抗力の強さは高速で小さくなる。等方な背景の $H(v)$ は $v$ が熱速度より大きいと $n_b/v$ に近い。従って速度の一次係数は $\vb A\simeq-K\vb v/v^3$、$K=2\sum_b\Gamma_{eb}(1+m_e/m_b)n_b$ である。このベクトル平均には向きの散乱が含まれる。速さだけの抗力を求めるには二次係数の補正も必要である。式 :eqref[eq-plasma-fp] を $v=\abs{\vb v}$ と合わせて部分積分すると、平均速さの変化は $\vb A\cdot\hat{\vb v}+\tr[(I-\hat{\vb v}\hat{\vb v})B]/(2v)$ である。高速極限で二次項は $2\sum_b\Gamma_{eb} n_b/v^2$ となり、残る減速は $-K_{\rm sl}/v^2$、$K_{\rm sl}=2\sum_b\Gamma_{eb}(m_e/m_b)n_b$ である。

電気力に沿う高速粒子を考え、角度散乱を一旦無視すると $eE/m_e>K_{\rm sl}/v^2$ のとき加速が勝つ。加速される側へ逸脱する粒子を :term[逃走電子]{en="runaway electron" yomi="とうそうでんし"} という。$E_D\sim m_e K_{\rm sl}/(e v_t^2)$ は熱速度で比較した Dreicer 電場の目安である。分布尾部からの供給率は速度空間の拡散と角度散乱で決まるので、この不等式だけでは出ない。相対論的速度では摩擦則も変わる。相対論的な臨界電場と二次生成率には相対論的衝突演算子が必要であるため、本章では非相対論的な加速の機構までを扱う。 :cite[goldston-rutherford]

:::example[電場を変えたとき]
同じ非相対論的背景で $E$ を四倍にすると、上の釣り合い速度 $\sqrt{m_e K_{\rm sl}/(eE)}$ は半分である。熱速度に近づいたときは高速近似を見直す。
:::

:::example[運動量とエネルギーの違い]
イオンを固定する Lorentz 模型では $m_b\to\infty$、同じ $\Gamma_{eb}n_b$ で $K_{\rm sl}\to0$ である。イオンによる角度散乱の率は有限でも、速い電子の速さそのものを下げる反跳は消える。したがって電気伝導で測る電流緩和率を、そのまま高速粒子のエネルギー損失率に使えない。
:::

### 高速近似での逃走しきい速度

高速・非相対論的な試験粒子の一次元模型では、電子が力を受ける向きに $m_e\dot v=e|E|-C/v^2$ と置ける。$e|E|=C/v_c^2$ の釣り合いは不安定であり、$v=v_c+\delta v$ を代入すると $\dot{\delta v}\simeq[2e|E|/(m_ev_c)]\delta v$ となる。初期速度がしきいより上なら抵抗が減ってさらに加速する。背景の Maxwell 裾にある割合を $\tfrac12\operatorname{erfc}(v_c/\sqrt2v_t)$ で見積もることはできるが、これは方向を一つにとった速度の初期割合である。速度拡散、ピッチ角散乱、相対論的な抵抗、二次生成は入っていない。低速の $v^{-2}$ の発散は模型の外である。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $E/E_*$ | $E_*=C/(ev_t^2)$ を用いた電場比 | 1 |
| $v/v_t$ | 試験粒子の初速度の比 | 1 |

形式的なしきい速度を $v_c/v_t$ とする。電気力と高速摩擦が等しい点である。

$$
v_c/v_t=\frac{1}{\sqrt{E/E_{*}}}
$$

初速での摩擦の大きさを $F_{\rm drag}/(eE_*)$ とする。電子熱速さより十分速い範囲で使う。

$$
F_{\rm drag}/(eE_*)=\frac{1}{(v/v_t)^2}
$$

規格化加速度を $A$ とする。正なら初めに加速する。

$$
A=E/E_{*} - F_{\rm drag}/(eE_*)
$$

しきい速度に対する初速を $v/v_c$ とする。一を越えれば摩擦が電気力より小さい。

$$
v/v_c=\frac{v/v_{t}}{(v_c/v_t)}
$$

しきい近くの線形化率を $G$ とする。単位は $eE_* /(m_ev_t)$ である。

$$
G=\frac{2 E/E_{*}}{(v_c/v_t)}
$$

Maxwell の片側での初期割合を $P_0$ とする。逃走粒子の総生成率ではない。

$$
P_0=\frac{\operatorname{erfc}{\left(\frac{\sqrt{2} (v_c/v_t)}{2} \right)}}{2}
$$

:::example[弱い電場でしきい以下]
初期加速度の符号を計算する。

入力値を次のように置く。

$$
\begin{aligned}
E/E_*&=\frac{1}{100} \\v/v_t&=5
\end{aligned}
$$

これらを順に代入する。

形式的なしきい速度は

$$
v_c/v_t=\frac{1}{\sqrt{\left(\frac{1}{100}\right)}}=10
$$

である。

初速での摩擦の大きさは

$$
F_{\rm drag}/(eE_*)=\frac{1}{\left(5\right)^{2}}=\frac{1}{25}
$$

である。

規格化加速度は

$$
A=\left(\frac{1}{100}\right) - \left(\frac{1}{25}\right)=- \frac{3}{100}
$$

である。

しきい速度に対する初速は

$$
v/v_c=\frac{\left(5\right)}{\left(10\right)}=\frac{1}{2}
$$

である。

しきい近くの線形化率は

$$
G=\frac{2 \left(\frac{1}{100}\right)}{\left(10\right)}=\frac{1}{500}
$$

である。

Maxwell の片側での初期割合は

$$
P_0=\frac{\operatorname{erfc}{\left(\frac{\sqrt{2} \left(10\right)}{2} \right)}}{2}\simeq 7.6199\times 10^{-24}
$$

である。

熱速さの数倍あっても高速摩擦が大きい。
:::

:::example[弱い電場でしきい以上]
同じ電場で初速だけを増やす。

入力値を次のように置く。

$$
\begin{aligned}
E/E_*&=\frac{1}{100} \\v/v_t&=20
\end{aligned}
$$

これらを順に代入する。

形式的なしきい速度は

$$
v_c/v_t=\frac{1}{\sqrt{\left(\frac{1}{100}\right)}}=10
$$

である。

初速での摩擦の大きさは

$$
F_{\rm drag}/(eE_*)=\frac{1}{\left(20\right)^{2}}=\frac{1}{400}
$$

である。

規格化加速度は

$$
A=\left(\frac{1}{100}\right) - \left(\frac{1}{400}\right)=\frac{3}{400}
$$

である。

しきい速度に対する初速は

$$
v/v_c=\frac{\left(20\right)}{\left(10\right)}=2
$$

である。

しきい近くの線形化率は

$$
G=\frac{2 \left(\frac{1}{100}\right)}{\left(10\right)}=\frac{1}{500}
$$

である。

Maxwell の片側での初期割合は

$$
P_0=\frac{\operatorname{erfc}{\left(\frac{\sqrt{2} \left(10\right)}{2} \right)}}{2}\simeq 7.6199\times 10^{-24}
$$

である。

電気力が勝ち、初めに正の加速度となる。
:::

:::example[強い電場としきい]
しきい速度とその不安定性を調べる。

入力値を次のように置く。

$$
\begin{aligned}
E/E_*&=\frac{1}{25} \\v/v_t&=5
\end{aligned}
$$

これらを順に代入する。

形式的なしきい速度は

$$
v_c/v_t=\frac{1}{\sqrt{\left(\frac{1}{25}\right)}}=5
$$

である。

初速での摩擦の大きさは

$$
F_{\rm drag}/(eE_*)=\frac{1}{\left(5\right)^{2}}=\frac{1}{25}
$$

である。

規格化加速度は

$$
A=0
$$

である。

しきい速度に対する初速は

$$
v/v_c=1
$$

である。

しきい近くの線形化率は

$$
G=\frac{2 \left(\frac{1}{25}\right)}{\left(5\right)}=\frac{2}{125}
$$

である。

Maxwell の片側での初期割合は

$$
P_0=\frac{\operatorname{erfc}{\left(\frac{\sqrt{2} \left(5\right)}{2} \right)}}{2}\simeq 2.8665\times 10^{-7}
$$

である。

加速度はゼロだが正の速度ずれは成長する。
:::

:::example[高速近似が崩れる電場]
形式的な Maxwell 割合を解釈する。

入力値を次のように置く。

$$
\begin{aligned}
E/E_*&=1 \\v/v_t&=2
\end{aligned}
$$

これらを順に代入する。

形式的なしきい速度は

$$
v_c/v_t=\frac{1}{\sqrt{\left(1\right)}}=1
$$

である。

初速での摩擦の大きさは

$$
F_{\rm drag}/(eE_*)=\frac{1}{\left(2\right)^{2}}=\frac{1}{4}
$$

である。

規格化加速度は

$$
A=\left(1\right) - \left(\frac{1}{4}\right)=\frac{3}{4}
$$

である。

しきい速度に対する初速は

$$
v/v_c=\frac{\left(2\right)}{\left(1\right)}=2
$$

である。

しきい近くの線形化率は

$$
G=2
$$

である。

Maxwell の片側での初期割合は

$$
P_0=\frac{\operatorname{erfc}{\left(\frac{\sqrt{2} \left(1\right)}{2} \right)}}{2}\simeq 0.15866
$$

である。

しきいが熱速さに等しく、高速近似から生成率は求められない。
:::


## 輸送係数と Braginskii 方程式の概観

衝突時間 $\tau_s$ が流体変化時間より短く、平均自由行程が巨視的長さより短ければ、分布を局所 Maxwell 分布と小さな補正に分けられる。Vlasov の移流と力を Maxwell 分布に作用させ、一次で線形化した衝突演算子と等置する。衝突が保存する粒子数、全運動量、全エネルギーに平行な成分は可解条件として流体の保存式を与える。補正の一次・二次・三次モーメントから摩擦、粘性、熱伝導を得る。これは :term[Braginskii 輸送]{en="Braginskii transport" yomi="ぶらぎんすきーゆそう"} の考え方である。

磁化した種の熱流束を $\vb h=-\kappa_\parallel\vb b\partial_\parallel T-\kappa_\perp\nabla_\perp T-\kappa_\wedge\vb b\times\nabla T$ と分ける。$\kappa_\wedge$ の項と $\nabla T$ の内積はゼロである。輸送の尺度は $\kappa_\parallel\sim n v_t^2\tau$、$\kappa_\perp\sim n\rho_t^2/\tau$ であり、$\kappa_\perp/\kappa_\parallel\sim(\Omega\tau)^{-2}$ となる。ここで $T$ はエネルギーなので $\kappa$ の単位は $\mathrm{m^{-1}s^{-1}}$ である。示した横方向の尺度は $|\Omega|\tau\gg1$ に限る。数値係数と熱力や粘性の全成分は、線形衝突演算子を速度空間のテンソル基底で解く必要があるため :cite[goldston-rutherford] を参照する。

:::example[磁化と熱輸送]
$\Omega\tau=100$ なら、強磁場のオーダー推定では横の熱伝導は平行の $10^{-4}$ 程度である。有限旋回半径が密度や温度の勾配尺度に達すると局所の閉包は使えない。
:::

:::example[局所閉包の適用範囲]
平均自由行程が $1\,\mathrm{cm}$ で温度が $1\,\mathrm m$ にわたって変化するなら、局所展開の空間比は $0.01$ である。辺の温度が $1\,\mathrm{mm}$ で変化すると同じ衝突条件で比は $10$ になる。後者では温度の局所勾配だけで熱流束を与える近似は制御できない。
:::

### 線形化と輸送係数の符号

共通温度・共通平均流の Maxwell 平衡 $f_{aM}$ の近くで $f_a=f_{aM}(1+h_a)$ と書く。式 :eqref[eq-plasma-landau-collision] の括弧を $m_a f_a f_b'\vb h$ で表したとき、平衡の対数微分の差は $\vb u$ に平行で $\mathsf U$ が消す。このため線形の括弧には $\nabla h_a/m_a-\nabla' h_b'/m_b$ だけが残る。線形化した衝突項を $-f_{aM}\mathcal L h_a$ と定義し、内積 $\langle g,h\rangle_M=\sum_a\int f_{aM}g_ah_a\dd^3v$ を用いる。部分積分と種の交換により
$$
\langle g,\mathcal L h\rangle_M
=\frac12\sum_{a,b}m_a^2\Gamma_{ab}\int f_{aM} f_{bM}'
\left(\frac{\nabla g_a}{m_a}-\frac{\nabla' g_b'}{m_b}\right)
\cdot\mathsf U
\left(\frac{\nabla h_a}{m_a}-\frac{\nabla' h_b'}{m_b}\right)
\dd^3v\dd^3v'
$$
である。従って $\mathcal L$ は対称で非負である。粒子数、全運動量、全エネルギーに対応する関数は右辺を消すので、$\mathcal L h=s$ を解くには、それらと $s$ が直交する必要がある。分布の密度・平均速度・温度を Maxwell 部に含めることで、補正の核成分を固定できる。

輸送の一次方程式を磁場なしで $\mathcal L h=s$ とすれば、同じ駆動に対応する輸送率は $\langle s,h\rangle_M=\langle h,\mathcal Lh\rangle_M\ge0$ である。各種の旋回を残す場合、$\mathcal L$ に内積に関して反対称な旋回作用素が加わる。その実二次形式はゼロなので、散逸の源は衝突に残る。これが符号を数値係数から独立に確認する方法である。

:::example[試験関数から得る伝導度の下界]
$\mathcal L h_*=s$ が核の直交補空間で解けるとする。任意の試験関数 $g$ に対し
$$
2\langle s,g\rangle_M-\langle g,\mathcal L g\rangle_M
=\langle s,h_*\rangle_M-\langle g-h_*,\mathcal L(g-h_*)\rangle_M
$$
である。従って左辺を一つの速度関数で最大化して得た値は真の輸送率を越えない。例えば $g=a g_0$、$\langle g_0,\mathcal Lg_0\rangle_M=C>0$、$\langle s,g_0\rangle_M=D$ とすれば最大は $a=D/C$、値は $D^2/C$ である。数値係数を使う場合も、衝突行列と駆動の規格化を先に確認する。
:::

### 平行輸送で熱力と熱流を結ぶ

一価イオンのプラズマで $m_e/m_i\ll1$ とし、$\tau_e$ は Lorentz 伝導度と比較した熱衝突時間の規約を用いる。平行相対流 $u_\parallel=u_{e\parallel}-u_{i\parallel}$ に対し、電子に働く平行摩擦力密度 $R_{e\parallel}$ と、電子流に乗った系での熱流 $h_{e\parallel}$ は
$$
\begin{aligned}
R_{e\parallel}&\simeq
-0.51\,\frac{m_e n_e}{\tau_e}u_\parallel
-0.71\,n_e\partial_\parallel T_e,\\
h_{e\parallel}&\simeq
-3.16\,\frac{n_e T_e\tau_e}{m_e}\partial_\parallel T_e
+0.71\,n_e T_e u_\parallel .
\end{aligned}
$$
数係数は電子同士と電子・イオンの線形衝突問題の解である。等温では電子運動式から平行伝導度 $n_e e^2\tau_e/(0.51 m_e)$ を得る。ここでは粘性テンソル全体の係数を求める数値解は扱わず、平行輸送で二つの駆動がどう結合するかを示す。 :cite[goldston-rutherford]

同じ規格化で Lorentz 模型の $\nu_*\tau_e=3\sqrt{\pi/2}$ を使うと、伝導度は $[32/(3\pi)]\,n_e e^2\tau_e/m_e$ である。したがって電子間衝突を含めた値はこの約 $0.58$ 倍になる。電子間衝突は全電子運動量を保つが、速度ごとの分布補正を変えるため、イオン摩擦で失われる合計電流にも効く。熱力と熱流の $0.71$ の項を落とすと、温度勾配と相対流を同時に与えた場合の相互応答を再現できない。

### 磁場に沿う熱流と横の熱流を分離する

輸送係数を得る完全な線形衝突問題の代わりに、固定の率 $\nu$ を持つ流束緩和式 $\nu\vb q-\Omega\vb q\times\hat{\vb b}=-C\grad T$ を調べる。$\hat{\vb b}=\hat{\vb z}$ で行列を二成分ずつ反転すると、平行は $\kappa_0=C/\nu$、垂直は $\kappa_0/(1+\beta^2)$、直交成分の大きさは $\kappa_0|\beta|/(1+\beta^2)$ である。$\beta=\Omega/\nu$ は符号を含む。温度勾配と直交する成分は $\vb q\cdot\grad T$ には寄与しない。この模型の係数を Braginskii の数係数と同一視しない。実際の係数は衝突率の速度依存と粒子流・熱流の結合を含む。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $\beta$ | 符号付き磁化比 | 1 |
| $X$ | 磁場に垂直な温度勾配の規格化値 | 1 |
| $Z$ | 磁場に平行な温度勾配の規格化値 | 1 |

横の対称係数を $P$ とする。行列の逆の対称部分である。

$$
P=\frac{1}{\beta^{2} + 1}
$$

横の反対称係数を $H$ とする。磁場反転で符号を変える。

$$
H=P \beta
$$

温度勾配に沿う横熱流を $Q_x$ とする。$\kappa_0$ は共通単位に含める。

$$
Q_x=- P X
$$

勾配に直交する横熱流を $Q_y$ とする。流束式の $y$ 成分で $q_y=-\beta q_x$ である。

$$
Q_y=H X
$$

磁場に沿う熱流を $Q_z$ とする。磁化による抑制がない。

$$
Q_z=- Z
$$

温度勾配と熱流の負の内積を $S_T$ とする。直交成分が熱的な散逸を直接担わないことを示す。

$$
S_T=- Q_{x} X - Q_{z} Z
$$

:::example[磁場なし]
二つの勾配の熱流を比較する。

入力値を次のように置く。

$$
\begin{aligned}
\beta&=0 \\X&=1 \\Z&=1
\end{aligned}
$$

これらを順に代入する。

横の対称係数は

$$
P=\frac{1}{\left(0\right)^{2} + 1}=1
$$

である。

横の反対称係数は

$$
H=\left(0\right) \left(1\right)=0
$$

である。

温度勾配に沿う横熱流は

$$
Q_x=- \left(1\right)^{2}=-1
$$

である。

勾配に直交する横熱流は

$$
Q_y=\left(0\right) \left(1\right)=0
$$

である。

磁場に沿う熱流は

$$
Q_z=- \left(1\right)=-1
$$

である。

温度勾配と熱流の負の内積は

$$
S_T=- 2 \left(-1\right) \left(1\right)=2
$$

である。

同じ係数であり、直交熱流が消える。
:::

:::example[強い磁化]
横熱流と磁場方向の熱流を計算する。

入力値を次のように置く。

$$
\begin{aligned}
\beta&=10 \\X&=1 \\Z&=1
\end{aligned}
$$

これらを順に代入する。

横の対称係数は

$$
P=\frac{1}{\left(10\right)^{2} + 1}=\frac{1}{101}
$$

である。

横の反対称係数は

$$
H=\left(10\right) \left(\frac{1}{101}\right)=\frac{10}{101}
$$

である。

温度勾配に沿う横熱流は

$$
Q_x=- \left(1\right) \left(\frac{1}{101}\right)=- \frac{1}{101}
$$

である。

勾配に直交する横熱流は

$$
Q_y=\left(1\right) \left(\frac{10}{101}\right)=\frac{10}{101}
$$

である。

磁場に沿う熱流は

$$
Q_z=- \left(1\right)=-1
$$

である。

温度勾配と熱流の負の内積は

$$
S_T=- \left(- \frac{1}{101}\right) \left(1\right) - \left(-1\right) \left(1\right)=\frac{102}{101}
$$

である。

対称な横伝導は強く抑制される。
:::

:::example[磁場の反転]
熱流のどの成分が変わるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
\beta&=-10 \\X&=1 \\Z&=1
\end{aligned}
$$

これらを順に代入する。

横の対称係数は

$$
P=\frac{1}{\left(-10\right)^{2} + 1}=\frac{1}{101}
$$

である。

横の反対称係数は

$$
H=\left(-10\right) \left(\frac{1}{101}\right)=- \frac{10}{101}
$$

である。

温度勾配に沿う横熱流は

$$
Q_x=- \left(1\right) \left(\frac{1}{101}\right)=- \frac{1}{101}
$$

である。

勾配に直交する横熱流は

$$
Q_y=\left(- \frac{10}{101}\right) \left(1\right)=- \frac{10}{101}
$$

である。

磁場に沿う熱流は

$$
Q_z=- \left(1\right)=-1
$$

である。

温度勾配と熱流の負の内積は

$$
S_T=- \left(- \frac{1}{101}\right) \left(1\right) - \left(-1\right) \left(1\right)=\frac{102}{101}
$$

である。

反対称成分だけが反転し、散逸の指標は同じである。
:::

:::example[垂直勾配だけ]
勾配と直交する熱流が大きい条件を求める。

入力値を次のように置く。

$$
\begin{aligned}
\beta&=2 \\X&=1 \\Z&=0
\end{aligned}
$$

これらを順に代入する。

横の対称係数は

$$
P=\frac{1}{\left(2\right)^{2} + 1}=\frac{1}{5}
$$

である。

横の反対称係数は

$$
H=\left(2\right) \left(\frac{1}{5}\right)=\frac{2}{5}
$$

である。

温度勾配に沿う横熱流は

$$
Q_x=- \left(1\right) \left(\frac{1}{5}\right)=- \frac{1}{5}
$$

である。

勾配に直交する横熱流は

$$
Q_y=\left(1\right) \left(\frac{2}{5}\right)=\frac{2}{5}
$$

である。

磁場に沿う熱流は

$$
Q_z=- \left(0\right)=0
$$

である。

温度勾配と熱流の負の内積は

$$
S_T=- \left(- \frac{1}{5}\right) \left(1\right) - \left(0\right)^{2}=\frac{1}{5}
$$

である。

直交成分が大きくても、勾配方向の成分が散逸を決める。
:::


## 演習問題

::::exercise{level=basic}
ドリフトがゼロで、速度の二次係数 $B=8\,\mathrm{m^2s^{-3}}$ が一定の一次元模型の分散の増加率は何か。

:::solution
式 :eqref[eq-plasma-fp] では $D_v=B/2$ なので $\dv{\langle v^2\rangle}{t}=2D_v=8\,\mathrm{m^2s^{-3}}$ である。平均速度は一定である。
:::
::::

::::exercise{level=basic}
固定散乱体の $\nu_D$ は速さが二倍なら何倍か。

:::solution
$v^{-3}$ のため $1/8$ 倍である。
:::
::::

::::exercise{level=basic}
$\Delta_{\mathbb S^2}\xi$ を求めよ。

:::solution
$\partial_\xi\xi=1$ なので $\partial_\xi(1-\xi^2)=-2\xi$ である。
:::
::::

::::exercise{level=basic}
$\Omega\tau=10^3$ の熱伝導比のオーダーは何か。

:::solution
$(\Omega\tau)^{-2}=10^{-6}$ である。
:::
::::

::::exercise{level=advanced}
Rosenbluth の $H$ が球対称であるとき $H'(v)$ を密度積分で書け。

:::solution
動径 Laplacian から $(v^2H')'=-4\pi v^2 f_b(v)$ である。原点で $v^2H'\to0$ なら $H'(v)=-4\pi v^{-2}\int_0^v s^2 f_b(s)\dd s$ となる。
:::
::::

::::exercise{level=advanced}
固定散乱体極限で速度ベクトルの平均が減っても、各粒子の速さが保存される理由を示せ。

:::solution
$\vb A=-2\Gamma n_b\hat{\vb v}/v^2$、$B=2\Gamma n_b(I-\hat{\vb v}\hat{\vb v})/v$ である。速さの二次補正は $2\Gamma n_b/v^2$ なので動径の一次係数を相殺する。さらに速さの二次係数は $\hat{\vb v}\cdot\mathsf B\hat{\vb v}=0$ であり、速さには拡散もない。
:::
::::

::::exercise{level=advanced}
$f_1$ が $E$ の一次であるとき、電流の符号が電荷の符号に依存しないことを示せ。

:::solution
$f_1=qEv_zf_0/(T\nu_D)$ を $\int qv_z f_1\dd^3v$ に入れると $q^2E$ と非負の積分の積になる。従って $\sigma>0$ である。
:::
::::

::::exercise{level=advanced}
$\Gamma_{ab}$ を保ったまま等方温度を九倍にしたとき Lorentz 伝導度は何倍か。

:::solution
$\nu_*\propto v_t^{-3}\propto T_e^{-3/2}$ なので $\sigma$ は $9^{3/2}=27$ 倍になる。
:::
::::

## 確認テスト

::::quiz
Landau 核がゼロにするベクトルは何か。

- [x] 相対速度 $\vb u$
- [ ] $\vb u$ に垂直な任意の非零ベクトル
- [ ] 全てのベクトル

:::explanation
$\mathsf U$ は $\vb u$ の垂直平面への射影を $u$ で割るため、$\mathsf U\vb u=0$ である。
:::
::::

::::quiz{answer="27" tol="1e-6"}
速さが三倍になったときの $\nu_D$ の逆数の倍率を答えよ。

:::explanation
$\nu_D^{-1}\propto v^3$ なので倍率は $3^3=27$ である。
:::
::::

::::quiz
固定散乱体の Lorentz 演算子はどれを変えるか。

- [x] 電子の速度の向き
- [ ] 電子の速さ
- [ ] 電子数

:::explanation
演算子は一定の速さの球面上だけで作用する。角度分布を変えるが、角度積分はゼロである。
:::
::::

::::quiz{answer="8" tol="1e-6"}
$A=0$、$B=4$ の一次元速度拡散で、2 秒間の分散増加を SI で答えよ。

:::explanation
分散の時間微分が $B=4\,\mathrm{m^2s^{-3}}$ なので増加は $8\,\mathrm{m^2s^{-2}}$ である。
:::
::::

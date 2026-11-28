---
title: 誘電テンソルと波動伝播
description: 磁化プラズマの誘電応答から分散面と伝播条件を求める。
notation:
  - symbol: '\mathsf K'
    meaning: 比誘電テンソル
  - symbol: 'S,D,P,R,L'
    meaning: Stix パラメータ
  - symbol: '\theta'
    meaning: 波数と背景磁場の角度
---

## この章の目標

- 円偏波の応答を直交座標で組み立てる。
- 斜め伝播の根と偏波を求める。
- 非一様媒質の反射と吸収を評価する。

## 冷たいプラズマの誘電テンソルと Stix パラメータ

電流を分極電流 $\vb j=-i\omega\vb P_{\rm pol}$ とみなし、$\vb E+\vb P_{\rm pol}/\varepsilon_0=\mathsf K\vb E$ と書く。したがって $\mathsf K=I+i\mathsf\sigma/(\varepsilon_0\omega)$ である。衝突を無視した :eqref[eq-plasma-rl] の円基底の応答は $R,L$、平行応答は $P=1-\sum_s\omega_{ps}^2/\omega^2$ である。円基底から直交基底へ戻すと
$$
\mathsf K=\begin{pmatrix}S&-iD&0\\iD&S&0\\0&0&P\end{pmatrix},\quad
S=(R+L)/2,\quad D=(R-L)/2 .
\label{eq-plasma-stix}
$$
これらを :term[Stix パラメータ]{en="Stix parameters" yomi="すてぃっくすぱらめーた"} と呼ぶ。極を除く実周波数では行列は Hermite である。したがって $\vb E^\dagger\mathsf K\vb E$ が実となり、時間平均の $\Re(\vb E^\dagger\vb j)/2$ はゼロになる。

:::example[偏波の異なる応答]
$S=3$、$D=1$ では等しい電場強度の R 波と L 波はそれぞれ $4$ と $2$ の誘電値を感じる。座標軸と円偏波の基底を固定したまま背景磁場を反転すると $D$ の符号が変わり、両応答は交換される。
:::

### 摩擦で極を正則化する

一様な種ごとの運動式に $-m_s\nu_s\vb u_s$ を足すと、R 波では $-i(\omega+\Omega_s+i\nu_s)m_s u_{sR}=q_s E_R$ である。このとき $R=1-\sum_s\omega_{ps}^2/[\omega(\omega+\Omega_s+i\nu_s)]$ となり、分極電流への換算に使う $\omega$ 自体は実周波数のままである。実数 $\omega>0$ では各和の虚部が正なので、円偏波の電場から粒子へエネルギーが渡る。

:::example[衝突で広がった共鳴]
イオンを固定し $\omega=\abs{\Omega_e}$ とすれば $\Im R=\omega_{pe}^2/(\omega\nu_e)$ である。同じ電場振幅で $\nu_e$ を半分にすると共鳴点での吸収は二倍になる。吸収を決める $\omega\Im R$ の半値全幅は半分になる。境界からの波では電場振幅自体も $\nu_e$ に依存する。
:::

### 摩擦を含む円偏波の吸収形状

一種の冷たい流体に $\nu\vb u$ の摩擦を加える。$e^{-i\omega t}$ 規約で円偏波を選ぶと慣性と旋回の行列は対角化し、$\sigma_\pm=n q^2/[m(\nu-i\Delta_\pm)]$、$\Delta_\pm=\omega\pm\Omega$ となる。$\Re\sigma_\pm$ の最大値は $\Delta_\pm=0$ にあり、半値全幅は角周波数で $2\nu$ である。$\int_{\R}\nu\,\dd\Delta/(\nu^2+\Delta^2)=\pi$ なので、衝突を弱めてもこの線形応答の吸収の面積は一定である。冷たい系における摩擦幅であり、温度による Doppler 幅や空間方向の共鳴層の幅とは区別する。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $\Delta$ | 選んだ円偏波の離調 | $\mathrm{s^{-1}}$ |
| $\nu$ | 摩擦率 | $\mathrm{s^{-1}}$ |
| $C$ | $n_s q_s^2/m_s$ の係数 | $\mathrm{S\,m^{-1}\,s^{-1}}$ |

線幅単位の離調を $x$ とする。共鳴の左右で符号が変わる。

$$
x=\frac{\Delta}{\nu}
$$

伝導度の実部を $\sigma_R$ とする。$\sigma_R|E|^2/2$ が平均の吸収である。

$$
\sigma_R=\frac{C}{\nu \left(x^{2} + 1\right)}
$$

伝導度の虚部を $\sigma_I$ とする。誘導的な電流成分である。

$$
\sigma_I=\sigma_{R} x
$$

同じ電場での最大加熱に対する比を $P/P_{\rm peak}$ とする。$\sigma_R$ を共鳴の値で割る。

$$
P/P_{\rm peak}=\frac{1}{x^{2} + 1}
$$

吸収線の半値全幅を $W$ とする。正負の半値点の差を取る。

$$
W=2 \nu
$$

全離調で積分した吸収線の面積を $A$ とする。$\nu$ は高さと幅で相殺する。

$$
A=\pi C
$$

:::example[共鳴点]
実部と虚部を分けて計算する。

入力値を次のように置く。

$$
\begin{aligned}
\Delta&=0\,\mathrm{s^{-1}} \\\nu&=2\,\mathrm{s^{-1}} \\C&=10\,\mathrm{S\,m^{-1}\,s^{-1}}
\end{aligned}
$$

これらを順に代入する。

線幅単位の離調は

$$
x=\frac{\left(0\right)}{\left(2\right)}=0
$$

である。

伝導度の実部は

$$
\sigma_R=\frac{\left(10\right)}{\left(2\right) \left(\left(0\right)^{2} + 1\right)}=5\,\mathrm{S\,m^{-1}}
$$

である。

伝導度の虚部は

$$
\sigma_I=\left(0\right) \left(5\right)=0\,\mathrm{S\,m^{-1}}
$$

である。

同じ電場での最大加熱に対する比は

$$
P/P_{\rm peak}=\frac{1}{\left(0\right)^{2} + 1}=1
$$

である。

吸収線の半値全幅は

$$
W=2 \left(2\right)=4\,\mathrm{s^{-1}}
$$

である。

全離調で積分した吸収線の面積は

$$
A=\pi \left(10\right)\simeq 31.416\,\mathrm{S\,m^{-1}\,s^{-1}}
$$

である。

共鳴点では応答は同相である。
:::

:::example[正側の半値点]
同じ摩擦でどの離調が半値になるか確かめる。

入力値を次のように置く。

$$
\begin{aligned}
\Delta&=2\,\mathrm{s^{-1}} \\\nu&=2\,\mathrm{s^{-1}} \\C&=10\,\mathrm{S\,m^{-1}\,s^{-1}}
\end{aligned}
$$

これらを順に代入する。

線幅単位の離調は

$$
x=1
$$

である。

伝導度の実部は

$$
\sigma_R=\frac{\left(10\right)}{\left(2\right) \left(\left(1\right)^{2} + 1\right)}=\frac{5}{2}\,\mathrm{S\,m^{-1}}
$$

である。

伝導度の虚部は

$$
\sigma_I=\left(1\right) \left(\frac{5}{2}\right)=\frac{5}{2}\,\mathrm{S\,m^{-1}}
$$

である。

同じ電場での最大加熱に対する比は

$$
P/P_{\rm peak}=\frac{1}{\left(1\right)^{2} + 1}=\frac{1}{2}
$$

である。

吸収線の半値全幅は

$$
W=2 \left(2\right)=4\,\mathrm{s^{-1}}
$$

である。

全離調で積分した吸収線の面積は

$$
A=\pi \left(10\right)\simeq 31.416\,\mathrm{S\,m^{-1}\,s^{-1}}
$$

である。

この点の虚部と実部の絶対値は一致する。
:::

:::example[負側の半値点]
符号を反転し、発熱が変わるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
\Delta&=-2\,\mathrm{s^{-1}} \\\nu&=2\,\mathrm{s^{-1}} \\C&=10\,\mathrm{S\,m^{-1}\,s^{-1}}
\end{aligned}
$$

これらを順に代入する。

線幅単位の離調は

$$
x=\frac{\left(-2\right)}{\left(2\right)}=-1
$$

である。

伝導度の実部は

$$
\sigma_R=\frac{\left(10\right)}{\left(2\right) \left(\left(-1\right)^{2} + 1\right)}=\frac{5}{2}\,\mathrm{S\,m^{-1}}
$$

である。

伝導度の虚部は

$$
\sigma_I=\left(-1\right) \left(\frac{5}{2}\right)=- \frac{5}{2}\,\mathrm{S\,m^{-1}}
$$

である。

同じ電場での最大加熱に対する比は

$$
P/P_{\rm peak}=\frac{1}{\left(-1\right)^{2} + 1}=\frac{1}{2}
$$

である。

吸収線の半値全幅は

$$
W=2 \left(2\right)=4\,\mathrm{s^{-1}}
$$

である。

全離調で積分した吸収線の面積は

$$
A=\pi \left(10\right)\simeq 31.416\,\mathrm{S\,m^{-1}\,s^{-1}}
$$

である。

発熱は変わらず、虚部だけが反転する。
:::

:::example[幅の狭い吸収線]
衝突率を減らして高さと積分値を比較する。

入力値を次のように置く。

$$
\begin{aligned}
\Delta&=0\,\mathrm{s^{-1}} \\\nu&=\frac{1}{10}\,\mathrm{s^{-1}} \\C&=10\,\mathrm{S\,m^{-1}\,s^{-1}}
\end{aligned}
$$

これらを順に代入する。

線幅単位の離調は

$$
x=\frac{\left(0\right)}{\left(\frac{1}{10}\right)}=0
$$

である。

伝導度の実部は

$$
\sigma_R=\frac{\left(10\right)}{\left(\frac{1}{10}\right) \left(\left(0\right)^{2} + 1\right)}=100\,\mathrm{S\,m^{-1}}
$$

である。

伝導度の虚部は

$$
\sigma_I=\left(0\right) \left(100\right)=0\,\mathrm{S\,m^{-1}}
$$

である。

同じ電場での最大加熱に対する比は

$$
P/P_{\rm peak}=\frac{1}{\left(0\right)^{2} + 1}=1
$$

である。

吸収線の半値全幅は

$$
W=2 \left(\frac{1}{10}\right)=\frac{1}{5}\,\mathrm{s^{-1}}
$$

である。

全離調で積分した吸収線の面積は

$$
A=\pi \left(10\right)\simeq 31.416\,\mathrm{S\,m^{-1}\,s^{-1}}
$$

である。

ピークが高くなる一方で面積は一定である。
:::


## 任意の伝播方向の分散関係

$\vb k$ を $xz$ 平面にとり、$\hat{\vb k}=(\sin\theta,0,\cos\theta)$ とする。式 :eqref[eq-plasma-wave-equation] から $\det(\mathsf K-N^2(I-\hat{\vb k}\hat{\vb k}^{\mathsf T}))=0$ である。$X=N^2$ とすると、行列の対角は $(S-X\cos^2\theta,S-X,P-X\sin^2\theta)$、非対角は $M_{xy}=-iD$、$M_{xz}=X\sin\theta\cos\theta$ とその共役である。第 2 行で展開すると
$$
\begin{aligned}
0&=AX^2-BX+C,\\
A&=S\sin^2\theta+P\cos^2\theta,\\
B&=RL\sin^2\theta+PS(1+\cos^2\theta),\\
C&=PRL .
\end{aligned}
\label{eq-plasma-oblique}
$$
三次項は横射影の階数が二であるため消える。判別式は $(RL-PS)^2\sin^4\theta+4P^2D^2\cos^2\theta$ と整理され、損失のない実係数模型では $N^2$ は実数である。$A=0$、$B\ne0$ では、二次式は線形式になり、有限根は $C/B$ である。$B=0$ も同時に成り立つ場合は、残る定数項と元の行列を調べる。

:::example[二つの特殊方向]
$\theta=0$ では式が $P(X-R)(X-L)=0$ となる。$\theta=\pi/2$ では $(X-P)(SX-RL)=0$ である。前者で $P=0$ は独立の縦モードを与える。
:::

:::example[斜め伝播の数値例]
$R=4,L=2$、$P=1$ とする。このとき $S=3,D=1$ で、$\theta=\pi/4$ なら $A=2,\ B=17/2,\ C=8$ となる。$N^2=(17\pm\sqrt{33})/8$ の二つの正根が得られる。垂直伝播に角度を変えると二根は $P=1$ と $8/3$ になる。
:::

### 縦成分の消去と分枝の数

波数に平行な軸と垂直な二軸に基底を取り替える。$\mathsf K$ をこの基底で
$$
\mathsf K=\begin{pmatrix}A&\vb h^\dagger\\ \vb h&\mathsf K_T\end{pmatrix}
$$
とする。ここでは冷たい無衝突で実周波数、共鳴の極の外を考える。縦行から $AE_\parallel=-\vb h^\dagger\vb E_T$ である。$A\ne0$ ならこれを横行へ戻し、
$$
\left(\mathsf K_T-\frac{\vb h\vb h^\dagger}{A}\right)\vb E_T=N^2\vb E_T
$$
を得る。従って二つの $N^2$ は Hermite な二次元行列の固有値として実数になる。この見方では、$A=0$ へ近づくと $\vb h$ の方向の固有値が大きくなりうることも分かる。$\vb h=0$ ならこの除算は不要であり、縦モードと横モードは結合しない。純粋な平行伝播で $P=0$ の場合がその例である。

:::example[符号の数を変える条件]
実対称な模型で $A=1$、$\mathsf K_T=\mathrm{diag}(2,3)$、$\vb h=(2,0)$ と置く。横の有効行列は $\mathrm{diag}(-2,3)$ であり、$\mathsf K_T$ の二つの固有値が正でも、縦成分と結合した後は一枝が減衰する。$N^2=-2$ の偏波には $E_\parallel=-2E_{T1}$ が伴う。
:::

### 特殊方向を結ぶ行列式の根

無衝突で $\vb k$ を $xz$ 面に取り、$s=\sin^2\theta$ とおく。$N^2$ の多項式の係数は $A=Ss+P(1-s)$、$B=RLs+PS(2-s)$、$C=PRL$ である。これは行列の $xy$ 結合の積が $D^2$ になり、残りの $S^2-D^2$ が $RL$ となるためである。判別式は $[(RL-PS)s]^2+4P^2 D^2(1-s)$ と平方和に整理できる。従って実数の誘電成分を用いると二つの形式根は実数であるが、正であるとは限らない。$A=0$ で $B\ne0$ なら二次式が一次式になるので、有限根を直接 $C/B$ から求める。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $R$ | 右円偏波の誘電率 | 1 |
| $L$ | 左円偏波の誘電率 | 1 |
| $P$ | 平行誘電率 | 1 |
| $s$ | $\sin^2\theta$ | 1 |

対角の横誘電率を $S$ とする。円偏波の平均を取る。

$$
S=\frac{L}{2} + \frac{R}{2}
$$

二次係数を $A$ とする。波数方向への誘電率の射影である。

$$
A=P \left(1 - s\right) + S s
$$

一次係数の符号を反転した量を $B$ とする。展開で残る横の積を使う。

$$
B=L R s + P S \left(2 - s\right)
$$

定数項を $C$ とする。波数ゼロでの誘電行列式に等しい。

$$
C=L P R
$$

正符号の平方根から得る枝を $N_1^2$ とする。$A$ が負なら数値の大小順とは逆になる。

$$
N_1^2=\frac{B + \sqrt{- 4 A C + B^{2}}}{2 A}
$$

もう一つの枝を $N_2^2$ とする。積を使い二つの近い値の減算を避ける。

$$
N_2^2=\frac{C}{A N^{2}_{1}}
$$

:::example[垂直入射の検査]
分枝が $P$ と $RL/S$ に分かれることを確認する。

入力値を次のように置く。

$$
\begin{aligned}
R&=3 \\L&=1 \\P&=2 \\s&=1
\end{aligned}
$$

これらを順に代入する。

対角の横誘電率は

$$
S=\frac{\left(1\right)}{2} + \frac{\left(3\right)}{2}=2
$$

である。

二次係数は

$$
A=\left(1\right) \left(2\right) + \left(2\right) \left(1 - \left(1\right)\right)=2
$$

である。

一次係数の符号を反転した量は

$$
B=\left(1\right)^{2} \left(3\right) + \left(2\right)^{2} \left(2 - \left(1\right)\right)=7
$$

である。

定数項は

$$
C=\left(1\right) \left(2\right) \left(3\right)=6
$$

である。

正符号の平方根から得る枝は

$$
N_1^2=\frac{\left(7\right) + \sqrt{- 4 \left(2\right) \left(6\right) + \left(7\right)^{2}}}{2 \left(2\right)}=2
$$

である。

もう一つの枝は

$$
N_2^2=\frac{\left(6\right)}{\left(2\right)^{2}}=\frac{3}{2}
$$

である。

典型的な O 波と X 波である。
:::

:::example[平行入射の検査]
二根が円偏波に一致するか確かめる。

入力値を次のように置く。

$$
\begin{aligned}
R&=3 \\L&=1 \\P&=2 \\s&=0
\end{aligned}
$$

これらを順に代入する。

対角の横誘電率は

$$
S=\frac{\left(1\right)}{2} + \frac{\left(3\right)}{2}=2
$$

である。

二次係数は

$$
A=\left(0\right) \left(2\right) + \left(2\right) \left(1 - \left(0\right)\right)=2
$$

である。

一次係数の符号を反転した量は

$$
B=\left(0\right) \left(1\right) \left(3\right) + \left(2\right)^{2} \left(2 - \left(0\right)\right)=8
$$

である。

定数項は

$$
C=\left(1\right) \left(2\right) \left(3\right)=6
$$

である。

正符号の平方根から得る枝は

$$
N_1^2=\frac{\left(8\right) + \sqrt{- 4 \left(2\right) \left(6\right) + \left(8\right)^{2}}}{2 \left(2\right)}=3
$$

である。

もう一つの枝は

$$
N_2^2=\frac{\left(6\right)}{\left(2\right) \left(3\right)}=1
$$

である。

縦応答はこの二根とは別であり、ここでは $P$ はゼロでない。
:::

:::example[斜めで一枝のみ伝播]
二つの実数根の符号を求める。

入力値を次のように置く。

$$
\begin{aligned}
R&=3 \\L&=-1 \\P&=2 \\s&=\frac{1}{2}
\end{aligned}
$$

これらを順に代入する。

対角の横誘電率は

$$
S=\frac{\left(-1\right)}{2} + \frac{\left(3\right)}{2}=1
$$

である。

二次係数は

$$
A=\left(1\right) \left(\frac{1}{2}\right) + \left(2\right) \left(1 - \left(\frac{1}{2}\right)\right)=\frac{3}{2}
$$

である。

一次係数の符号を反転した量は

$$
B=\left(-1\right) \left(3\right) \left(\frac{1}{2}\right) + \left(1\right) \left(2\right) \left(2 - \left(\frac{1}{2}\right)\right)=\frac{3}{2}
$$

である。

定数項は

$$
C=\left(-1\right) \left(2\right) \left(3\right)=-6
$$

である。

正符号の平方根から得る枝は

$$
N_1^2=\frac{\left(\frac{3}{2}\right) + \sqrt{- 4 \left(-6\right) \left(\frac{3}{2}\right) + \left(\frac{3}{2}\right)^{2}}}{2 \left(\frac{3}{2}\right)}\simeq 2.5616
$$

である。

もう一つの枝は

$$
N_2^2=\frac{\left(-6\right)}{\left(\frac{1}{2} + \frac{\sqrt{17}}{2}\right) \left(\frac{3}{2}\right)}\simeq -1.5616
$$

である。

正根の偏波だけがこの一様媒質で実波数を持つ。
:::

:::example[共鳴角に近づく]
二次係数がゼロへ近づく側で根を比較する。

入力値を次のように置く。

$$
\begin{aligned}
R&=3 \\L&=1 \\P&=-1 \\s&=\frac{3}{10}
\end{aligned}
$$

これらを順に代入する。

対角の横誘電率は

$$
S=\frac{\left(1\right)}{2} + \frac{\left(3\right)}{2}=2
$$

である。

二次係数は

$$
A=\left(-1\right) \left(1 - \left(\frac{3}{10}\right)\right) + \left(2\right) \left(\frac{3}{10}\right)=- \frac{1}{10}
$$

である。

一次係数の符号を反転した量は

$$
B=\left(-1\right) \left(2\right) \left(2 - \left(\frac{3}{10}\right)\right) + \left(1\right) \left(3\right) \left(\frac{3}{10}\right)=- \frac{5}{2}
$$

である。

定数項は

$$
C=\left(-1\right) \left(1\right) \left(3\right)=-3
$$

である。

正符号の平方根から得る枝は

$$
N_1^2=\frac{\left(- \frac{5}{2}\right) + \sqrt{- 4 \left(- \frac{1}{10}\right) \left(-3\right) + \left(- \frac{5}{2}\right)^{2}}}{2 \left(- \frac{1}{10}\right)}\simeq 1.2639
$$

である。

もう一つの枝は

$$
N_2^2=\frac{\left(-3\right)}{\left(- \frac{1}{10}\right) \left(\frac{25}{2} - \frac{\sqrt{505}}{2}\right)}\simeq 23.736
$$

である。

一方の根だけが大きくなり、もう一方は有限の値へ近づく。
:::


## CMA 図

:term[CMA 図]{en="Clemmow–Mullaly–Allis diagram" yomi="しーえむえーず"} では、冷たいプラズマの波面法線の分散面を密度と磁場の平面上で区分する。ここでいう分散面は、周波数を固定して許される実ベクトル $c\vb k/\omega$ を描いた面である。固定した正の $\omega$ で $X_e=\omega_{pe}^2/\omega^2$、$Y=\abs{\Omega_e}/\omega$ ととる。イオン固定なら
$$
P=1-X_e,\quad R=1-X_e/(1-Y),\quad L=1-X_e/(1+Y),\quad S=1-X_e/(1-Y^2).
$$
図を作るには、$P,R,L,S$ の零曲線と $Y=1$ の極を描き、各領域で $N^2$ の符号と発散角を求める。$PS<0$ なら $A=0$ を満たす角 $\tan^2\theta=-P/S$ がある。根の正負は積 $C/A$ と和 $B/A$ から判定する。イオンを含む図では、同じ手順でイオン共鳴の極と混成共鳴線を追加する。

:::example[図上の領域]
$(X_e,Y)=(1/2,1/4)$ なら $P=1/2$、$R=1/3$、$L=3/5$、$S=7/15$ である。$A,B,C$ が正なので二根とも正である。$P$ と $S$ は同符号で発散角はない。
:::

### 電子のみを動かす CMA 図の領域

$0<Y<1$ に限定すると、密度を増やしたとき零点の順は $X_e=1-Y$（$R$）、$1-Y^2$（$S$）、$1$（$P$）、$1+Y$（$L$）である。後の三つの値がこの順になるのは $Y^2<Y$ と $Y>0$ による。垂直方向を選べば O 波の符号は $P$、X 波の符号は $RL/S$ だけで判定できる。ここで零点と極そのものを除く。

| $X_e$ の範囲 | O 波 | X 波 |
| --- | --- | --- |
| $0<X_e<1-Y$ | 伝播 | 伝播 |
| $1-Y<X_e<1-Y^2$ | 伝播 | 減衰 |
| $1-Y^2<X_e<1$ | 伝播 | 伝播 |
| $1<X_e<1+Y$ | 減衰 | 伝播 |
| $1+Y<X_e$ | 減衰 | 減衰 |

表は垂直方向に限る。斜め方向には $A=0$ の角度が現れうるので、表から波面全体を推定せず式 :eqref[eq-plasma-oblique] を使う。

:::example[垂直波が一つの領域]
$Y=1/2$、$X_e=3/5$ は $1-Y<X_e<1-Y^2$ である。$P=2/5$、$R=-1/5$、$L=3/5$、$S=1/5$ なので O 波だけが垂直に伝播する。
:::

### 角度を連続に変える判定

垂直の二枝と平行の二枝が分かっていても、その間の接続は $A$ と $C$ に依存する。$C\ne0$ のまま角度を変えると、有限の根がゼロを通ることはない。正根の数が変わるには、$A=0$ に達して根が無限大を通る必要がある。$A=0$ の非縮退点で $B\ne0$ なら、一つの根は $C/B$ に近づき、他方は和 $B/A$ から発散する。$A$ は $\sin^2\theta$ に関する一次式なので、$0<\theta<\pi/2$ には発散角が高々一つである。これは CMA 図の各小領域で斜め方向の伝播枝を数える手順を与える。

:::example[発散角をまたぐ枝]
$P>0,S<0$ で $C\ne0,B\ne0$ とする。$\theta=0$ では $A=P>0$、$\theta=\pi/2$ では $A=S<0$ である。$\tan^2\theta=-P/S$ で一枝だけが発散し、その $N^2$ の符号が変わる。交換の前後で伝播枝を結ぶとき、根の大きさ順の名前が保存されるとは限らない。
:::

### 全角度での伝播枝の数

$X_e>0$ とし、端点、極、カットオフ、共鳴円錐上を除く。$Y>1$ では $R,S>0$ であり、$P=0$ の直線が $L=0$ の直線より低密度にある。$PS<0$ の区間では共鳴角 $\theta_c$ の前後で $A$ の符号が反転する。有限根の符号は変わらず、一つの根が発散を経て符号を変える。これを平行と垂直の根に適用すると、正の $N^2$ を持つ枝の数は次のようになる。

| 磁場の範囲 | 密度の範囲 | 正の根の数 |
| --- | --- | --- |
| $0<Y<1$ | $0<X_e<1-Y$ | 2 |
| $0<Y<1$ | $1-Y<X_e<1-Y^2$ | 1 |
| $0<Y<1$ | $1-Y^2<X_e<1$ | $\theta<\theta_c$ で 1、$\theta>\theta_c$ で 2 |
| $0<Y<1$ | $1<X_e<1+Y$ | 1 |
| $0<Y<1$ | $1+Y<X_e$ | 0 |
| $Y>1$ | $0<X_e<1$ | 2 |
| $Y>1$ | $1<X_e<1+Y$ | $\theta<\theta_c$ で 2、$\theta>\theta_c$ で 1 |
| $Y>1$ | $1+Y<X_e$ | $\theta<\theta_c$ で 1、$\theta>\theta_c$ で 0 |

例えば $X_e=3/2,Y=2$ では $P=-1/2,S=3/2$ なので $\theta_c=30^\circ$ である。平行の二根は $R=5/2,L=1/2$ ともに正だが、垂直の O 波の $N^2$ は負で、X 波の $N^2=5/6$ は正である。零点が重なる $Y=0$ では $R=L=S=P$ の等方問題に戻る。

### CMA 座標を密度と磁場から読む

電子のみが動くとし $X=\omega_{pe}^2/\omega^2$、$Y=\abs{\Omega_e}/\omega>0$ とおく。$R=1-X/(1-Y)$、$L=1-X/(1+Y)$、$P=1-X$ である。$Y<1$ なら正密度での R のカットオフは $X=1-Y$、L は $X=1+Y$ となる。$S=0$ は $X=1-Y^2$ なので、$0<Y<1$ では R のカットオフと O 波のカットオフの中間にある。$Y>1$ では R は正であり $S=0$ の正密度解も消える。こうした領域図から特定方向の伝播の有無は読み取れるが、角度を指定しない波数の大きさは分散式で求める。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $X$ | 密度を非磁化の臨界密度で割った値 | 1 |
| $Y$ | 電子旋回角周波数の絶対値と波の周波数の比 | 1 |

R の誘電成分を $R$ とする。電子旋回の極を持つ。

$$
R=- \frac{X}{1 - Y} + 1
$$

L の誘電成分を $L$ とする。正の周波数で電子の極を持たない。

$$
L=- \frac{X}{Y + 1} + 1
$$

O 波の二乗屈折率を $P$ とする。磁場のない横波と同じである。

$$
P=1 - X
$$

X 波の共鳴を支配する成分を $S$ とする。偶関数の磁場依存性を持つ。

$$
S=\frac{L}{2} + \frac{R}{2}
$$

垂直の X 波の二乗屈折率を $N_X^2$ とする。分子と分母の符号を確認する。

$$
N_X^2=\frac{L R}{S}
$$

上部混成の密度位置からのずれを $X-(1-Y^2)$ とする。$Y>1$ ならこの共鳴位置は正密度にない。

$$
X-(1-Y^2)=X + Y^{2} - 1
$$

:::example[低密度側の透過]
全方位の判定を始める前に特殊方向を確認する。

入力値を次のように置く。

$$
\begin{aligned}
X&=\frac{1}{4} \\Y&=\frac{1}{2}
\end{aligned}
$$

これらを順に代入する。

R の誘電成分は

$$
R=- \frac{\left(\frac{1}{4}\right)}{1 - \left(\frac{1}{2}\right)} + 1=\frac{1}{2}
$$

である。

L の誘電成分は

$$
L=- \frac{\left(\frac{1}{4}\right)}{\left(\frac{1}{2}\right) + 1} + 1=\frac{5}{6}
$$

である。

O 波の二乗屈折率は

$$
P=1 - \left(\frac{1}{4}\right)=\frac{3}{4}
$$

である。

X 波の共鳴を支配する成分は

$$
S=\frac{\left(\frac{1}{2}\right)}{2} + \frac{\left(\frac{5}{6}\right)}{2}=\frac{2}{3}
$$

である。

垂直の X 波の二乗屈折率は

$$
N_X^2=\frac{\left(\frac{1}{2}\right) \left(\frac{5}{6}\right)}{\left(\frac{2}{3}\right)}=\frac{5}{8}
$$

である。

上部混成の密度位置からのずれは

$$
X-(1-Y^2)=\left(\frac{1}{2}\right)^{2} + \left(\frac{1}{4}\right) - 1=- \frac{1}{2}
$$

である。

R と L、O と X の二乗屈折率はいずれも正である。
:::

:::example[R のカットオフ直後]
R と X が同時に実波数を失うかを調べる。

入力値を次のように置く。

$$
\begin{aligned}
X&=\frac{3}{5} \\Y&=\frac{1}{2}
\end{aligned}
$$

これらを順に代入する。

R の誘電成分は

$$
R=- \frac{\left(\frac{3}{5}\right)}{1 - \left(\frac{1}{2}\right)} + 1=- \frac{1}{5}
$$

である。

L の誘電成分は

$$
L=- \frac{\left(\frac{3}{5}\right)}{\left(\frac{1}{2}\right) + 1} + 1=\frac{3}{5}
$$

である。

O 波の二乗屈折率は

$$
P=1 - \left(\frac{3}{5}\right)=\frac{2}{5}
$$

である。

X 波の共鳴を支配する成分は

$$
S=\frac{\left(- \frac{1}{5}\right)}{2} + \frac{\left(\frac{3}{5}\right)}{2}=\frac{1}{5}
$$

である。

垂直の X 波の二乗屈折率は

$$
N_X^2=\frac{\left(- \frac{1}{5}\right) \left(\frac{3}{5}\right)}{\left(\frac{1}{5}\right)}=- \frac{3}{5}
$$

である。

上部混成の密度位置からのずれは

$$
X-(1-Y^2)=\left(\frac{1}{2}\right)^{2} + \left(\frac{3}{5}\right) - 1=- \frac{3}{20}
$$

である。

この範囲では R が負かつ $S$ が正なので、X も負になる。
:::

:::example[上部混成より高密度]
X が再び伝播するか調べる。

入力値を次のように置く。

$$
\begin{aligned}
X&=\frac{9}{10} \\Y&=\frac{1}{2}
\end{aligned}
$$

これらを順に代入する。

R の誘電成分は

$$
R=- \frac{\left(\frac{9}{10}\right)}{1 - \left(\frac{1}{2}\right)} + 1=- \frac{4}{5}
$$

である。

L の誘電成分は

$$
L=- \frac{\left(\frac{9}{10}\right)}{\left(\frac{1}{2}\right) + 1} + 1=\frac{2}{5}
$$

である。

O 波の二乗屈折率は

$$
P=1 - \left(\frac{9}{10}\right)=\frac{1}{10}
$$

である。

X 波の共鳴を支配する成分は

$$
S=\frac{\left(- \frac{4}{5}\right)}{2} + \frac{\left(\frac{2}{5}\right)}{2}=- \frac{1}{5}
$$

である。

垂直の X 波の二乗屈折率は

$$
N_X^2=\frac{\left(- \frac{4}{5}\right) \left(\frac{2}{5}\right)}{\left(- \frac{1}{5}\right)}=\frac{8}{5}
$$

である。

上部混成の密度位置からのずれは

$$
X-(1-Y^2)=\left(\frac{1}{2}\right)^{2} + \left(\frac{9}{10}\right) - 1=\frac{3}{20}
$$

である。

R と $S$ がともに負なので X は正となる。
:::

:::example[電子旋回より低い波の周波数]
$Y>1$ の領域で O と X を比べる。

入力値を次のように置く。

$$
\begin{aligned}
X&=\frac{3}{2} \\Y&=2
\end{aligned}
$$

これらを順に代入する。

R の誘電成分は

$$
R=- \frac{\left(\frac{3}{2}\right)}{1 - \left(2\right)} + 1=\frac{5}{2}
$$

である。

L の誘電成分は

$$
L=- \frac{\left(\frac{3}{2}\right)}{\left(2\right) + 1} + 1=\frac{1}{2}
$$

である。

O 波の二乗屈折率は

$$
P=1 - \left(\frac{3}{2}\right)=- \frac{1}{2}
$$

である。

X 波の共鳴を支配する成分は

$$
S=\frac{\left(\frac{1}{2}\right)}{2} + \frac{\left(\frac{5}{2}\right)}{2}=\frac{3}{2}
$$

である。

垂直の X 波の二乗屈折率は

$$
N_X^2=\frac{\left(\frac{1}{2}\right) \left(\frac{5}{2}\right)}{\left(\frac{3}{2}\right)}=\frac{5}{6}
$$

である。

上部混成の密度位置からのずれは

$$
X-(1-Y^2)=\left(2\right)^{2} + \left(\frac{3}{2}\right) - 1=\frac{9}{2}
$$

である。

O のカットオフは磁場に依存しないが、垂直 X の条件は変化する。
:::


## 共鳴とカットオフ

カットオフでは $N^2=0$ となる。式 :eqref[eq-plasma-oblique] からその候補は $P,R,L$ の零点である。共鳴には粒子の旋回と同調して誘電テンソル自体が発散する場合と、テンソルが有限でも $A=0$ で波数が発散する場合がある。後者では大きな $N$ に対し横電場が抑えられ、$\hat{\vb k}\cdot\mathsf K\hat{\vb k}=A=0$ が縦の応答を決める。

:::example[高域混成の見え方]
垂直伝播では $A=S$ である。$S=0$ で $RL\ne0$ なら X 波の $N^2$ が発散する。O 波は $N^2=P$ なので同じ発散を持たない。
:::

数式上の極で実際に発散する振幅を予測することはできない。有限温度、衝突、非一様性により極近傍の模型を更新する必要がある。複数分枝が接近する領域での接続係数の一般形は、各モードの連立境界問題に依存するため :cite[stix-waves] を参照する。

:::example[電子 X 波の二つのカットオフ]
イオンを固定し $w_c=\abs{\Omega_e}$ とする。$R=0$ は $\omega(\omega-w_c)=\omega_{pe}^2$、$L=0$ は $\omega(\omega+w_c)=\omega_{pe}^2$ である。正の根はそれぞれ $(\sqrt{w_c^2+4\omega_{pe}^2}\pm w_c)/2$ であり、その差は $w_c$ である。混成周波数はこれらの間にある。実際、混成周波数の二乗から低い根の二乗を引くと $w_c^2+w_c\omega_L>0$、高い根の二乗から混成周波数の二乗を引くと $w_c(\omega_R-w_c)>0$ である。
:::

### 波束とエネルギーの向き

一様な冷たい無衝突媒質で
$$
\mathcal U=\frac14\left(\varepsilon_0\abs E^2+\frac{\abs{B_1}^2}{\mu_0}
+\sum_s m_sn_{s0}\abs{u_s}^2\right)
$$
は時間平均のエネルギー密度である。円基底で成分ごとの $u_{sR}\propto 1/(\omega+\Omega_s)$、$u_{sL}\propto1/(\omega-\Omega_s)$ と平行成分を使う。実際 $\partial_\omega(\omega R)=1+\sum_s\omega_{ps}^2/(\omega+\Omega_s)^2$ であり、L と P でも同じ操作をすると
$$
\mathcal U=\frac14\left(\varepsilon_0\vb E^\dagger\partial_\omega(\omega\mathsf K)\vb E+\frac{\abs{B_1}^2}{\mu_0}\right)
$$
となる。分散行列を $\mathsf M=\omega^2\mathsf K/c^2+\vb k\vb k^\mathsf T-k^2I$ とする。非縮退枝で $\mathsf M\vb E=0$ を $k_j$ で微分し、左から $\vb E^\dagger$ を掛けると $\partial_{k_j}\omega=-\vb E^\dagger(\partial_{k_j}\mathsf M)\vb E/[\vb E^\dagger(\partial_\omega\mathsf M)\vb E]$ である。分子の符号を反転したベクトルは $2\Re[\vb E^*\times(\vb k\times\vb E)]$、分母は $4\mu_0\omega\mathcal U$ となる。Faraday 則を使えば群速度は時間平均 Poynting 流束を $\mathcal U$ で割った量に等しい。冷たい模型で粒子の圧力や平均移流がないため、粒子のエネルギー流束は別に現れない。

:::example[境界へ運べる電力]
$N^2\to0$ に近づく枝では、有限の電場振幅に対して $B_1=\vb k\times\vb E/\omega$ が小さくなる。したがって同じ振幅でも境界に垂直な電力は小さくなる。非一様問題で入射電力を一定にしたまま転回点に近づくと、ここで固定した振幅の比較は適用できない。
:::

### 共鳴円錐の角度と波数発散

斜めの波の $N^2$ が大きいとき、縦電場の係数は $A=S\sin^2\theta+P\cos^2\theta$ である。$SP<0$ で $A=0$ となる角が存在し、$\tan^2\theta_c=-P/S$ である。$B,C$ を現在の角で評価し、$|AC/B^2|\ll1$、$B\ne0$ とすると、大きい根は $N^2=B/A-C/B+\order{A}$、有限根は $C/B+\order{A}$ となる。二次式に大きな根を入れると $AN^2$ と $B$ が釣り合い、二根の和と積が支配項を決める。この計算は吸収のない冷たい極限であり、発散の近くで有限旋回半径や衝突を落とせるかどうかは別に調べる。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $S$ | 横の対角誘電率 | 1 |
| $P$ | 平行誘電率 | 1 |
| $V$ | 円偏波の誘電率の積 $RL=S^2-D^2$ | 1 |
| $\delta$ | $\sin^2\theta-\sin^2\theta_c$ | 1 |

共鳴角の二乗正弦を $\sin^2\theta_c$ とする。$SP<0$ なら零と一の間である。

$$
\sin^2\theta_c=- \frac{P}{- P + S}
$$

共鳴角付近の二次係数を $A$ とする。$A$ は角の二乗正弦に対して線形である。

$$
A=\delta \left(- P + S\right)
$$

一次係数の符号を反転した値を $B$ とする。現在の角で展開係数を評価する。

$$
B=P S \left(- \delta - \sin^2\theta_c + 2\right) + V \left(\delta + \sin^2\theta_c\right)
$$

発散枝の二項近似を $N_{\rm sing}^2$ とする。二根の和から有限枝の先頭項を引く。

$$
N_{\rm sing}^2=- \frac{P V}{B} + \frac{B}{A}
$$

根の分離を調べる小量を $\abs{AC/B^2}$ とする。有限根と発散根のスケール比になる。

$$
\abs{AC/B^2}=\left|{\frac{A P V}{B^{2}}}\right|
$$

共鳴円錐の半角を $\theta_c$ とする。正の主値を度に変換する。

$$
\theta_c=\frac{180 \arcsin{\left(\sqrt{\sin^2\theta_c} \right)}}{\pi}
$$

:::example[小角側の進行波]
二次係数と一次係数が同符号の側を調べる。

入力値を次のように置く。

$$
\begin{aligned}
S&=1 \\P&=-1 \\V&=-1 \\\delta&=- \frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

共鳴角の二乗正弦は

$$
\sin^2\theta_c=- \frac{\left(-1\right)}{- \left(-1\right) + \left(1\right)}=\frac{1}{2}
$$

である。

共鳴角付近の二次係数は

$$
A=\left(- \frac{1}{100}\right) \left(- \left(-1\right) + \left(1\right)\right)=- \frac{1}{50}
$$

である。

一次係数の符号を反転した値は

$$
B=\left(-1\right) \left(1\right) \left(- \left(- \frac{1}{100}\right) - \left(\frac{1}{2}\right) + 2\right) + \left(-1\right) \left(\left(- \frac{1}{100}\right) + \left(\frac{1}{2}\right)\right)=-2
$$

である。

発散枝の二項近似は

$$
N_{\rm sing}^2=- \frac{\left(-1\right)^{2}}{\left(-2\right)} + \frac{\left(-2\right)}{\left(- \frac{1}{50}\right)}=\frac{201}{2}
$$

である。

根の分離を調べる小量は

$$
\abs{AC/B^2}=\left|{\frac{\left(- \frac{1}{50}\right) \left(-1\right)^{2}}{\left(-2\right)^{2}}}\right|=\frac{1}{200}
$$

である。

共鳴円錐の半角は

$$
\theta_c=\frac{180 \arcsin{\left(\sqrt{\left(\frac{1}{2}\right)} \right)}}{\pi}=45\,\mathrm{deg}
$$

である。

発散枝の $N^2$ は正である。
:::

:::example[大角側の減衰波]
同じ成分の角だけを共鳴の反対側へ移す。

入力値を次のように置く。

$$
\begin{aligned}
S&=1 \\P&=-1 \\V&=-1 \\\delta&=\frac{1}{100}
\end{aligned}
$$

これらを順に代入する。

共鳴角の二乗正弦は

$$
\sin^2\theta_c=- \frac{\left(-1\right)}{- \left(-1\right) + \left(1\right)}=\frac{1}{2}
$$

である。

共鳴角付近の二次係数は

$$
A=\left(\frac{1}{100}\right) \left(- \left(-1\right) + \left(1\right)\right)=\frac{1}{50}
$$

である。

一次係数の符号を反転した値は

$$
B=\left(-1\right) \left(1\right) \left(- \left(\frac{1}{100}\right) - \left(\frac{1}{2}\right) + 2\right) + \left(-1\right) \left(\left(\frac{1}{100}\right) + \left(\frac{1}{2}\right)\right)=-2
$$

である。

発散枝の二項近似は

$$
N_{\rm sing}^2=- \frac{\left(-1\right)^{2}}{\left(-2\right)} + \frac{\left(-2\right)}{\left(\frac{1}{50}\right)}=- \frac{199}{2}
$$

である。

根の分離を調べる小量は

$$
\abs{AC/B^2}=\left|{\frac{\left(-1\right)^{2} \left(\frac{1}{50}\right)}{\left(-2\right)^{2}}}\right|=\frac{1}{200}
$$

である。

共鳴円錐の半角は

$$
\theta_c=\frac{180 \arcsin{\left(\sqrt{\left(\frac{1}{2}\right)} \right)}}{\pi}=45\,\mathrm{deg}
$$

である。

大きな根の符号が反転する。
:::

:::example[偏った共鳴角]
角度が四十五度に固定されないことを確認する。

入力値を次のように置く。

$$
\begin{aligned}
S&=3 \\P&=-1 \\V&=-7 \\\delta&=\frac{1}{50}
\end{aligned}
$$

これらを順に代入する。

共鳴角の二乗正弦は

$$
\sin^2\theta_c=- \frac{\left(-1\right)}{- \left(-1\right) + \left(3\right)}=\frac{1}{4}
$$

である。

共鳴角付近の二次係数は

$$
A=\left(\frac{1}{50}\right) \left(- \left(-1\right) + \left(3\right)\right)=\frac{2}{25}
$$

である。

一次係数の符号を反転した値は

$$
B=\left(-1\right) \left(3\right) \left(- \left(\frac{1}{4}\right) - \left(\frac{1}{50}\right) + 2\right) + \left(-7\right) \left(\left(\frac{1}{4}\right) + \left(\frac{1}{50}\right)\right)=- \frac{177}{25}
$$

である。

発散枝の二項近似は

$$
N_{\rm sing}^2=\frac{\left(- \frac{177}{25}\right)}{\left(\frac{2}{25}\right)} - \frac{\left(-1\right) \left(-7\right)}{\left(- \frac{177}{25}\right)}=- \frac{30979}{354}
$$

である。

根の分離を調べる小量は

$$
\abs{AC/B^2}=\left|{\frac{\left(-1\right) \left(-7\right) \left(\frac{2}{25}\right)}{\left(- \frac{177}{25}\right)^{2}}}\right|=\frac{350}{31329}
$$

である。

共鳴円錐の半角は

$$
\theta_c=\frac{180 \arcsin{\left(\sqrt{\left(\frac{1}{4}\right)} \right)}}{\pi}=30\,\mathrm{deg}
$$

である。

$S/P$ の比から角は三十度になる。
:::

:::example[根の尺度が十分離れない場合]
近似した発散根を信用できるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
S&=1 \\P&=-1 \\V&=1 \\\delta&=\frac{1}{4}
\end{aligned}
$$

これらを順に代入する。

共鳴角の二乗正弦は

$$
\sin^2\theta_c=- \frac{\left(-1\right)}{- \left(-1\right) + \left(1\right)}=\frac{1}{2}
$$

である。

共鳴角付近の二次係数は

$$
A=\left(\frac{1}{4}\right) \left(- \left(-1\right) + \left(1\right)\right)=\frac{1}{2}
$$

である。

一次係数の符号を反転した値は

$$
B=\left(-1\right) \left(1\right) \left(- \left(\frac{1}{2}\right) - \left(\frac{1}{4}\right) + 2\right) + \left(1\right) \left(\left(\frac{1}{2}\right) + \left(\frac{1}{4}\right)\right)=- \frac{1}{2}
$$

である。

発散枝の二項近似は

$$
N_{\rm sing}^2=\frac{\left(- \frac{1}{2}\right)}{\left(\frac{1}{2}\right)} - \frac{\left(-1\right) \left(1\right)}{\left(- \frac{1}{2}\right)}=-3
$$

である。

根の分離を調べる小量は

$$
\abs{AC/B^2}=\left|{\frac{\left(-1\right) \left(1\right) \left(\frac{1}{2}\right)}{\left(- \frac{1}{2}\right)^{2}}}\right|=2
$$

である。

共鳴円錐の半角は

$$
\theta_c=\frac{180 \arcsin{\left(\sqrt{\left(\frac{1}{2}\right)} \right)}}{\pi}=45\,\mathrm{deg}
$$

である。

比 $|AC/B^2|$ が小さくないので、角が近くても元の二次式へ戻る。
:::


## WKB 近似と波の伝播・吸収

定常な一次元の媒質に単一偏波を入れ、局所形 $E''+k(x)^2E=0$ が成り立つとする。$E=a e^{i\int k\dd x}$ を代入すると $a''+2ika'+ik'a=0$ である。例えば振幅の残差から $|k'|/|k|^2\ll1$ と $|k''|/|k|^3\ll1$ を確認する。このように波数が局所波長にわたって滑らかに変化する場合の :term[WKB 近似]{en="WKB approximation" yomi="だぶりゅーけーびーきんじ"} では $a''$ を落とし、$a\propto k^{-1/2}$ とする。実 $k$ に対して $k\abs a^2$ は一定であり、もとの式から $\Im(E^*E')$ も厳密に一定である。

### 非一様媒質での波線

位相 $\Theta(\vb x,t)$ がある領域で滑らかに変わるとし、$\vb k=\nabla\Theta$、$\omega=-\partial_t\Theta$ と置く。局所分散式 $\mathcal D(\omega,\vb k,\vb x,t)=0$ が、$\partial_\omega\mathcal D\ne0$ の一枝で $\omega=\Omega(\vb k,\vb x,t)$ と解けるとする。位相の空間微分により
$$
\partial_t k_i+(\partial_{k_j}\Omega)\partial_{x_i}k_j=-\partial_{x_i}\Omega
$$
となる。$\partial_{x_i}k_j=\partial_{x_j}k_i$ を使えば、特性曲線で
$$
\dot{\vb x}=\nabla_{\vb k}\Omega,\qquad
\dot{\vb k}=-\nabla_{\vb x}\Omega,\qquad \dot\omega=\partial_t\Omega
$$
を得る。最後の式は $\Omega$ の全微分で空間項と波数項が打ち消すことから従う。定常媒質では時間周波数を保ったまま波数が変わる。

例として非磁化の一様な断面内で $\Omega^2=\omega_p(x)^2+c^2(k_x^2+k_y^2)$ とする。$y$ によらないので $k_y$ は保存する。密度上昇に伴って $k_x^2=(\omega^2-\omega_p(x)^2)/c^2-k_y^2$ がゼロになると、法線方向の群速度もゼロになる。ここから先の反射を定量化するには、波線で使った短波長近似を転回点の解に接続する。

### 転回点と減衰

$k^2\simeq -\alpha(x-x_c)$、$\alpha>0$ なら $z=\alpha^{1/3}(x-x_c)$ で $E_{zz}-zE=0$ となる。減衰側で Airy 関数 $\mathrm{Ai}(z)$ を選ぶと、伝播側の漸近形は $\pi^{-1/2}(-z)^{-1/4}\sin(2(-z)^{3/2}/3+\pi/4)$ である。位相を $\Phi(x)=\int_x^{x_c} k(x')\dd x'$ と取る。時間依存 $e^{-i\omega t}$ では $e^{-i\Phi}$ が入射、$e^{i\Phi}$ が反射である。上の正弦をこの二項に分けると反射振幅と入射振幅の比は $-i$ となる。したがって強度は全反射し、位相は $-\pi/2$ ずれる。転回点を含む解には Airy 関数を使う。弱い損失で $k=k_r+ik_i$ のとき、進行波の強度は幾何的振幅を除いて $e^{-2\int k_i\dd x}$ だけ変化する。

:::example[振幅と強度の区別]
$k_i=0.1\,\mathrm{m^{-1}}$ が $10\,\mathrm m$ 続けば、振幅には $e^{-1}$、強度には $e^{-2}$ が掛かる。
:::

### 減衰層の通過

二つの転回点 $x_1<x_2$ の間だけ $k^2<0$ とする。十分厚い減衰層では減衰側の WKB 解が $\exp[-\int_{x_1}^{x}\abs k\dd x]$ を含む。両転回点の接続を行ったときの透過強度の先頭の指数因子は $\exp[-2\int_{x_1}^{x_2}\abs k\dd x]$ である。接続で生じる前置係数と干渉はこの指数評価に入っていない。媒質に散逸がない場合もこの透過率は小さくなり、その分は主に反射に回る。局所的な虚波数だけで吸収を結論できない。

:::example[障壁の厚さ]
層内の主な領域で $\abs k\simeq5\,\mathrm{m^{-1}}$、距離が $1\,\mathrm m$ なら、透過強度の指数因子は $e^{-10}$ である。同じ値の $k_i$ を持つ散逸した進行波と指数は似ているが、エネルギーの行き先は境界の反射で判別できる。
:::

### 転回点の Airy 解と障壁積分

$E''+k_0^2 N^2(x)E=0$ で単純転回点の近くを $k_0^2N^2=-a(x-x_t)$、$a>0$ とする。$\xi=a^{1/3}(x-x_t)$ に写すと $E_{\xi\xi}-\xi E=0$ である。減衰端での条件から Airy の $\operatorname{Ai}$ を選ぶ。侵入長の尺度は $\ell=a^{-1/3}$ で、WKB の相対勾配 $|\kappa'|/\kappa^2=1/(2\xi^{3/2})$ を代入すると $|\xi|\lesssim1$ では近似が破れる。減衰側の作用は $\int_{x_t}^x \sqrt{a(x'-x_t)}\dd x'=(2/3)\xi^{3/2}$ である。以下の指数は作用による抑制を表し、Airy の漸近形の $\xi^{-1/4}$ や境界接続係数を含まない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $a$ | 転回点における $-k_0^2(N^2)'$ | $\mathrm{m^{-3}}$ |
| $d$ | 減衰側の転回点からの距離 | $\mathrm{m}$ |

Airy 層の尺度を $\ell$ とする。二階微分と線形ポテンシャルが競合する長さである。

$$
\ell=\frac{1}{\sqrt[3]{a}}
$$

規格化距離を $\xi$ とする。Airy 方程式の独立変数である。

$$
\xi=\frac{d}{\ell}
$$

局所減衰波数を $\kappa$ とする。実数の正の平方根を使う。

$$
\kappa=\sqrt{a d}
$$

禁制側の作用積分を $S_a$ とする。波数を転回点から積分する。

$$
S_a=\frac{2 \xi^{\frac{3}{2}}}{3}
$$

強度の指数抑制因子を $M_{\rm int}$ とする。振幅指数の二倍である。

$$
M_{\rm int}=e^{- 2 S_{a}}
$$

WKB の局所勾配の指標を $\delta_{\rm WKB}$ とする。波数の逆長さとの比を取る。

$$
\delta_{\rm WKB}=\frac{1}{2 \xi^{\frac{3}{2}}}
$$

:::example[Airy 層の境界付近]
一層の距離で WKB の誤差尺度を出す。

入力値を次のように置く。

$$
\begin{aligned}
a&=8\,\mathrm{m^{-3}} \\d&=\frac{1}{2}\,\mathrm{m}
\end{aligned}
$$

これらを順に代入する。

Airy 層の尺度は

$$
\ell=\frac{1}{\sqrt[3]{\left(8\right)}}=\frac{1}{2}\,\mathrm{m}
$$

である。

規格化距離は

$$
\xi=1
$$

である。

局所減衰波数は

$$
\kappa=\sqrt{\left(8\right) \left(\frac{1}{2}\right)}=2\,\mathrm{m^{-1}}
$$

である。

禁制側の作用積分は

$$
S_a=\frac{2 \left(1\right)^{\frac{3}{2}}}{3}=\frac{2}{3}
$$

である。

強度の指数抑制因子は

$$
M_{\rm int}=e^{- 2 \left(\frac{2}{3}\right)}\simeq 0.2636
$$

である。

WKB の局所勾配の指標は

$$
\delta_{\rm WKB}=\frac{1}{2 \left(1\right)^{\frac{3}{2}}}=\frac{1}{2}
$$

である。

相対勾配が十分小さくないため接続に Airy 解を使う。
:::

:::example[深い禁制側]
層を数倍越えた位置の抑制を求める。

入力値を次のように置く。

$$
\begin{aligned}
a&=8\,\mathrm{m^{-3}} \\d&=2\,\mathrm{m}
\end{aligned}
$$

これらを順に代入する。

Airy 層の尺度は

$$
\ell=\frac{1}{\sqrt[3]{\left(8\right)}}=\frac{1}{2}\,\mathrm{m}
$$

である。

規格化距離は

$$
\xi=\frac{\left(2\right)}{\left(\frac{1}{2}\right)}=4
$$

である。

局所減衰波数は

$$
\kappa=\sqrt{\left(2\right) \left(8\right)}=4\,\mathrm{m^{-1}}
$$

である。

禁制側の作用積分は

$$
S_a=\frac{2 \left(4\right)^{\frac{3}{2}}}{3}=\frac{16}{3}
$$

である。

強度の指数抑制因子は

$$
M_{\rm int}=e^{- 2 \left(\frac{16}{3}\right)}\simeq 2.3309\times 10^{-5}
$$

である。

WKB の局所勾配の指標は

$$
\delta_{\rm WKB}=\frac{1}{2 \left(4\right)^{\frac{3}{2}}}=\frac{1}{16}
$$

である。

作用が大きくなって指数因子が急速に小さくなる。
:::

:::example[勾配の急な層]
同じ距離で勾配だけを急にした場合を比べる。

入力値を次のように置く。

$$
\begin{aligned}
a&=125\,\mathrm{m^{-3}} \\d&=\frac{1}{2}\,\mathrm{m}
\end{aligned}
$$

これらを順に代入する。

Airy 層の尺度は

$$
\ell=\frac{1}{\sqrt[3]{\left(125\right)}}=\frac{1}{5}\,\mathrm{m}
$$

である。

規格化距離は

$$
\xi=\frac{\left(\frac{1}{2}\right)}{\left(\frac{1}{5}\right)}=\frac{5}{2}
$$

である。

局所減衰波数は

$$
\kappa=\sqrt{\left(125\right) \left(\frac{1}{2}\right)}\simeq 7.9057\,\mathrm{m^{-1}}
$$

である。

禁制側の作用積分は

$$
S_a=\frac{2 \left(\frac{5}{2}\right)^{\frac{3}{2}}}{3}\simeq 2.6352
$$

である。

強度の指数抑制因子は

$$
M_{\rm int}=e^{- 2 \left(\frac{5 \sqrt{10}}{6}\right)}\simeq 0.0051412
$$

である。

WKB の局所勾配の指標は

$$
\delta_{\rm WKB}=\frac{1}{2 \left(\frac{5}{2}\right)^{\frac{3}{2}}}\simeq 0.12649
$$

である。

Airy 層が狭まり、同じ距離でも大きい規格化距離に相当する。
:::

:::example[転回点の直近]
強度の指数が一に近い値を示しても WKB を使えるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
a&=8\,\mathrm{m^{-3}} \\d&=\frac{1}{100}\,\mathrm{m}
\end{aligned}
$$

これらを順に代入する。

Airy 層の尺度は

$$
\ell=\frac{1}{\sqrt[3]{\left(8\right)}}=\frac{1}{2}\,\mathrm{m}
$$

である。

規格化距離は

$$
\xi=\frac{\left(\frac{1}{100}\right)}{\left(\frac{1}{2}\right)}=\frac{1}{50}
$$

である。

局所減衰波数は

$$
\kappa=\sqrt{\left(8\right) \left(\frac{1}{100}\right)}\simeq 0.28284\,\mathrm{m^{-1}}
$$

である。

禁制側の作用積分は

$$
S_a=\frac{2 \left(\frac{1}{50}\right)^{\frac{3}{2}}}{3}\simeq 0.0018856
$$

である。

強度の指数抑制因子は

$$
M_{\rm int}=e^{- 2 \left(\frac{\sqrt{2}}{750}\right)}\simeq 0.99624
$$

である。

WKB の局所勾配の指標は

$$
\delta_{\rm WKB}=\frac{1}{2 \left(\frac{1}{50}\right)^{\frac{3}{2}}}\simeq 176.78
$$

である。

相対勾配は大きい。指数だけの平坦な外挿は接続解にならない。
:::


### 二つの転回点を持つ滑らかな障壁

一次元の例として $E''+Q(x)E=0$、$Q(x)=a(x^2-L^2)$、$a,L>0$ とする。両側は進行波領域、中央の $-L<x<L$ は減衰領域である。単独の転回点での減衰を二つ足し合わせるには、両方の Airy 層が重ならないことが必要である。各点の $|Q'|=2aL$ だから層厚は $\ell=(2aL)^{-1/3}$ である。$\ell/L\ll1$ の場合、$x=L\sin\theta$ と置換して
$$
S=\int_{-L}^{L}\sqrt{-Q(x)}\,\dd x
  =\sqrt a L^2\int_{-\pi/2}^{\pi/2}\cos^2\theta\,\dd\theta
  =\frac{\pi}{2}\sqrt a\,L^2
$$
を得る。強度透過の指数部分は $e^{-2S}$ である。これは電力の輸送が保存するスカラー方程式での減衰であり、共鳴吸収のような熱散逸を意味しない。実係数の方程式では $E^*E''-EE^{*\,\prime\prime}=0$ より $W=\Im(E^*E')$ が一定である。左側の入射・反射と右側の透過に $W$ を評価すれば、正しく流束規格化した透過率と反射率の和は一である。指数因子の前に付く量や位相を決めるには接続解が要る。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $a$ | 二次分布の係数 | $\mathrm{m^{-4}}$ |
| $L$ | 二転回点の間の距離の半分 | $\mathrm m$ |

一転回点の Airy 層の長さを $\ell$ とする。微分の傾きから求める。

$$
\ell=\frac{2^{\frac{2}{3}}}{2 \sqrt[3]{L a}}
$$

接続領域の分離の指標を $\ell/L$ とする。一より十分小さい場合に二点を別々に扱える。

$$
\ell/L=\frac{\ell}{L}
$$

障壁中央の減衰波数を $\kappa(0)$ とする。$-Q(0)$ の平方根である。

$$
\kappa(0)=L \sqrt{a}
$$

全作用積分を $S$ とする。中央値を曲線の全幅に掛けただけでは係数が変わる。

$$
S=\frac{\pi L \kappa(0)}{2}
$$

透過の対数の先頭項を $\ln\mathcal T_{\rm exp}$ とする。深い障壁での大きさを比較しやすい。

$$
\ln\mathcal T_{\rm exp}=- 2 S
$$

強度の指数因子を $\mathcal T_{\rm exp}$ とする。厳密な透過率との比較には接続係数も要る。

$$
\mathcal T_{\rm exp}=e^{\ln\mathcal T_{\rm exp}}
$$

:::example[層が重なる境界]
作用に加えて Airy の長さ比を調べる。

入力値を次のように置く。

$$
\begin{aligned}
a&=1\,\mathrm{m^{-4}} \\L&=1\,\mathrm m
\end{aligned}
$$

これらを順に代入する。

一転回点の Airy 層の長さは

$$
\ell=\frac{2^{\frac{2}{3}}}{2 \sqrt[3]{\left(1\right)^{2}}}\simeq 0.7937\,\mathrm m
$$

である。

接続領域の分離の指標は

$$
\ell/L=\frac{\left(\frac{2^{\frac{2}{3}}}{2}\right)}{\left(1\right)}\simeq 0.7937
$$

である。

障壁中央の減衰波数は

$$
\kappa(0)=\left(1\right)^{\frac{3}{2}}=1\,\mathrm{m^{-1}}
$$

である。

全作用積分は

$$
S=\frac{\pi \left(1\right)^{2}}{2}\simeq 1.5708
$$

である。

透過の対数の先頭項は

$$
\ln\mathcal T_{\rm exp}=- 2 \left(\frac{\pi}{2}\right)\simeq -3.1416
$$

である。

強度の指数因子は

$$
\mathcal T_{\rm exp}=e^{\left(- \pi\right)}\simeq 0.043214
$$

である。

作用だけから精密に透過を求めるには層の分離が不足している。
:::

:::example[同じ高さで厚い障壁]
中央で同じ減衰波数を持つ障壁と比べる。

入力値を次のように置く。

$$
\begin{aligned}
a&=\frac{1}{4}\,\mathrm{m^{-4}} \\L&=2\,\mathrm m
\end{aligned}
$$

これらを順に代入する。

一転回点の Airy 層の長さは

$$
\ell=\frac{2^{\frac{2}{3}}}{2 \sqrt[3]{\left(2\right) \left(\frac{1}{4}\right)}}=1\,\mathrm m
$$

である。

接続領域の分離の指標は

$$
\ell/L=\frac{\left(1\right)}{\left(2\right)}=\frac{1}{2}
$$

である。

障壁中央の減衰波数は

$$
\kappa(0)=\left(2\right) \sqrt{\left(\frac{1}{4}\right)}=1\,\mathrm{m^{-1}}
$$

である。

全作用積分は

$$
S=\frac{\pi \left(1\right) \left(2\right)}{2}\simeq 3.1416
$$

である。

透過の対数の先頭項は

$$
\ln\mathcal T_{\rm exp}=- 2 \left(\pi\right)\simeq -6.2832
$$

である。

強度の指数因子は

$$
\mathcal T_{\rm exp}=e^{\left(- 2 \pi\right)}\simeq 0.0018674
$$

である。

最大の局所減衰は同じでも作用は長さに比例して増える。
:::

:::example[同じ幅で高い障壁]
係数だけを変えて作用と転回層を求める。

入力値を次のように置く。

$$
\begin{aligned}
a&=100\,\mathrm{m^{-4}} \\L&=1\,\mathrm m
\end{aligned}
$$

これらを順に代入する。

一転回点の Airy 層の長さは

$$
\ell=\frac{2^{\frac{2}{3}}}{2 \sqrt[3]{\left(100\right) \left(1\right)}}\simeq 0.171\,\mathrm m
$$

である。

接続領域の分離の指標は

$$
\ell/L=\frac{\left(\frac{\sqrt[3]{5}}{10}\right)}{\left(1\right)}\simeq 0.171
$$

である。

障壁中央の減衰波数は

$$
\kappa(0)=\sqrt{\left(100\right)} \left(1\right)=10\,\mathrm{m^{-1}}
$$

である。

全作用積分は

$$
S=\frac{\pi \left(10\right) \left(1\right)}{2}\simeq 15.708
$$

である。

透過の対数の先頭項は

$$
\ln\mathcal T_{\rm exp}=- 2 \left(5 \pi\right)\simeq -31.416
$$

である。

強度の指数因子は

$$
\mathcal T_{\rm exp}=e^{\left(- 10 \pi\right)}\simeq 2.2711\times 10^{-14}
$$

である。

指数抑制が大きくなり、転回層は狭くなる。
:::

:::example[十分分離した転回点]
層厚と半長を比較して指数近似の条件を調べる。

入力値を次のように置く。

$$
\begin{aligned}
a&=4\,\mathrm{m^{-4}} \\L&=3\,\mathrm m
\end{aligned}
$$

これらを順に代入する。

一転回点の Airy 層の長さは

$$
\ell=\frac{2^{\frac{2}{3}}}{2 \sqrt[3]{\left(3\right) \left(4\right)}}\simeq 0.34668\,\mathrm m
$$

である。

接続領域の分離の指標は

$$
\ell/L=\frac{\left(\frac{3^{\frac{2}{3}}}{6}\right)}{\left(3\right)}\simeq 0.11556
$$

である。

障壁中央の減衰波数は

$$
\kappa(0)=\left(3\right) \sqrt{\left(4\right)}=6\,\mathrm{m^{-1}}
$$

である。

全作用積分は

$$
S=\frac{\pi \left(3\right) \left(6\right)}{2}\simeq 28.274
$$

である。

透過の対数の先頭項は

$$
\ln\mathcal T_{\rm exp}=- 2 \left(9 \pi\right)\simeq -56.549
$$

である。

強度の指数因子は

$$
\mathcal T_{\rm exp}=e^{\left(- 18 \pi\right)}\simeq 2.762\times 10^{-25}
$$

である。

この設定では障壁内の WKB 領域が広く、指数因子は極めて小さい。
:::


## 波動加熱と電流駆動

加熱には電場と電流が同位相になる成分が要る。時間平均は $P_{\rm abs}=\Re(\vb E^\dagger\vb j)/2=\varepsilon_0\omega\,\vb E^\dagger(\mathsf K-\mathsf K^\dagger)\vb E/(4i)$ である。平行 Drude 模型で $m(\nu-i\omega)u=qE$ とすれば $\sigma=nq^2/[m(\nu-i\omega)]$ なので $P_{\rm abs}=nq^2\nu\abs E^2/[2m(\nu^2+\omega^2)]$ となる。

### 波が渡す運動量

進行する平行静電波 $\phi=\phi(kz-\omega t)$ では $\partial_t(q\phi)=-(\omega/k)\partial_z(q\phi)$ である。一粒子 Hamiltonian $H=p_z^2/(2m)+q\phi$ は $\dot H=\partial_t(q\phi)$、$\dot p_z=-\partial_z(q\phi)$ を満たすので $\dot p_z=(k/\omega)\dot H$ である。選択的に平行運動量を渡すと、電子とイオンの流れがずれて電流が生じる。定常電流を計算するには、共鳴電子の空間分布と衝突による運動量交換も必要である。電子サイクロトロン加熱では共回転する偏波を選ぶ。案内中心が $z=v_\parallel t$ を進むと、基本波の位相は粒子から角周波数 $\omega-k_\parallel v_\parallel$ に見える。これが $\abs{\Omega_e}$ に一致すると、円運動に対する電場の向きが保たれて正味の仕事をする。斜め入射では旋回による位置の変化も位相に入り、高調波の共鳴が生じる。有限温度と高調波ごとの応答は [磁化プラズマの運動論](/physics/plasma-physics/magnetized-kinetics/) で扱う。

電流には運動量の注入に加えて、速度による摩擦率の違いも効く。例えば二群の電子に等しく逆向きの運動量を単位体積・単位時間当たり $\pm F$ だけ与え、各群の運動量密度 $P_\pm$ が摩擦率 $\nu_\pm$ で減るとする。定常では $P_\pm=\pm F/\nu_\pm$ であり、電子電流は $-(e/m_e)F(1/\nu_+-1/\nu_-)$ となる。総注入力がゼロでも摩擦率が異なれば電流は残る。平行速度の片側を選んで電子を加熱する電流駆動では、この速度選択性を評価する必要がある。

:::example[逆向きの波]
同じ $\omega>0$ で $k$ の符号だけが異なる二つの静電波が同じ電力を吸収させると、上の模型で力密度は相殺する。加熱率は加算される。
:::

:::example[吸収に向く周波数]
固定した $\omega$ に対し、Drude 模型の $\nu/(\nu^2+\omega^2)$ は $\nu=\omega$ で最大である。低衝突側では電子が場とほぼ直交位相で動き、高衝突側では電流振幅が小さくなる。場の空間分布を固定できる場合に限った最大である。
:::

### 波の吸収による電力と平行運動量収支

弱減衰を仮定し、固定した実周波数で電場が $E\propto e^{ik_r x-k_i x}$ と減るとする。波の強度に比例した電力流 $P_{\rm w}$ は $\partial_x P_{\rm w}=-2k_iP_{\rm w}$ に従う。反射と幾何学的な拡がりを無視すれば $P_{\rm abs}=P_{\rm in}(1-e^{-2k_iL})$ である。背景が磁場方向に一様で、波の一つの Fourier 成分が吸収されると、エネルギーと平行運動量の交換比は $k_\parallel/\omega$ である。粒子が位相 $k_\parallel z-\omega t$ の関数である場を見る場合、運動する座標の Hamiltonian $H-(\omega/k_\parallel)p_\parallel$ が保存することから従う。力から定常電流を得るにはさらに衝突と加速の速度分布を解く必要がある。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $P_{\rm in}$ | 入射波の電力 | $\mathrm W$ |
| $k_i$ | 波数の減衰成分 | $\mathrm{m^{-1}}$ |
| $L$ | 一様な吸収層の長さ | $\mathrm m$ |
| $Z$ | 平行運動量とエネルギーの比 $k_\parallel/\omega$ | $\mathrm{s\,m^{-1}}$ |

層の光学厚さを $\tau$ とする。強度で定義するため係数二が付く。

$$
\tau=2 L k_{i}
$$

電力の透過割合を $\mathcal T$ とする。反射を無視している。

$$
\mathcal T=e^{- \tau}
$$

透過電力を $P_{\rm out}$ とする。入射電力に割合を掛ける。

$$
P_{\rm out}=P_{\rm in} \mathcal T
$$

層の吸収電力を $P_{\rm abs}$ とする。失われた波の電力が粒子に渡る。

$$
P_{\rm abs}=P_{\rm in} - P_{\rm out}
$$

吸収による平行力を $F_\parallel$ とする。運動量移行率は電力に $k_\parallel/\omega$ を掛ける。

$$
F_\parallel=P_{\rm abs} Z
$$

層入口での単位長当たり吸収を $S_a$ とする。電力流の減少率を入口で評価する。

$$
S_a=2 P_{\rm in} k_{i}
$$

:::example[長い吸収長を持つ波]
吸収電力と力を計算する。

入力値を次のように置く。

$$
\begin{aligned}
P_{\rm in}&=1000\,\mathrm W \\k_i&=\frac{1}{100}\,\mathrm{m^{-1}} \\L&=1\,\mathrm m \\Z&=\frac{1}{1000000}\,\mathrm{s\,m^{-1}}
\end{aligned}
$$

これらを順に代入する。

層の光学厚さは

$$
\tau=2 \left(1\right) \left(\frac{1}{100}\right)=\frac{1}{50}
$$

である。

電力の透過割合は

$$
\mathcal T=e^{- \left(\frac{1}{50}\right)}\simeq 0.9802
$$

である。

透過電力は

$$
P_{\rm out}=\left(1000\right) \left(e^{- \frac{1}{50}}\right)\simeq 980.2\,\mathrm W
$$

である。

層の吸収電力は

$$
P_{\rm abs}=\left(1000\right) - \left(\frac{1000}{e^{\frac{1}{50}}}\right)\simeq 19.801\,\mathrm W
$$

である。

吸収による平行力は

$$
F_\parallel=\left(1000 - \frac{1000}{e^{\frac{1}{50}}}\right) \left(\frac{1}{1000000}\right)\simeq 1.9801\times 10^{-5}\,\mathrm N
$$

である。

層入口での単位長当たり吸収は

$$
S_a=2 \left(1000\right) \left(\frac{1}{100}\right)=20\,\mathrm{W\,m^{-1}}
$$

である。

吸収割合が小さいので入射電力をそのまま加熱電力にできない。
:::

:::example[光学厚さが一]
通り抜ける電力の割合を求める。

入力値を次のように置く。

$$
\begin{aligned}
P_{\rm in}&=1000\,\mathrm W \\k_i&=\frac{1}{2}\,\mathrm{m^{-1}} \\L&=1\,\mathrm m \\Z&=\frac{1}{1000000}\,\mathrm{s\,m^{-1}}
\end{aligned}
$$

これらを順に代入する。

層の光学厚さは

$$
\tau=2 \left(1\right) \left(\frac{1}{2}\right)=1
$$

である。

電力の透過割合は

$$
\mathcal T=e^{- \left(1\right)}\simeq 0.36788
$$

である。

透過電力は

$$
P_{\rm out}=\left(1000\right) \left(e^{-1}\right)\simeq 367.88\,\mathrm W
$$

である。

層の吸収電力は

$$
P_{\rm abs}=\left(1000\right) - \left(\frac{1000}{e}\right)\simeq 632.12\,\mathrm W
$$

である。

吸収による平行力は

$$
F_\parallel=\left(1000 - \frac{1000}{e}\right) \left(\frac{1}{1000000}\right)\simeq 0.00063212\,\mathrm N
$$

である。

層入口での単位長当たり吸収は

$$
S_a=2 \left(1000\right) \left(\frac{1}{2}\right)=1000\,\mathrm{W\,m^{-1}}
$$

である。

強度は $1/e$ が残り、電場の振幅比とは異なる。
:::

:::example[反対符号の平行波数]
吸収率を固定して力の向きを変える。

入力値を次のように置く。

$$
\begin{aligned}
P_{\rm in}&=1000\,\mathrm W \\k_i&=\frac{1}{2}\,\mathrm{m^{-1}} \\L&=1\,\mathrm m \\Z&=- \frac{1}{1000000}\,\mathrm{s\,m^{-1}}
\end{aligned}
$$

これらを順に代入する。

層の光学厚さは

$$
\tau=2 \left(1\right) \left(\frac{1}{2}\right)=1
$$

である。

電力の透過割合は

$$
\mathcal T=e^{- \left(1\right)}\simeq 0.36788
$$

である。

透過電力は

$$
P_{\rm out}=\left(1000\right) \left(e^{-1}\right)\simeq 367.88\,\mathrm W
$$

である。

層の吸収電力は

$$
P_{\rm abs}=\left(1000\right) - \left(\frac{1000}{e}\right)\simeq 632.12\,\mathrm W
$$

である。

吸収による平行力は

$$
F_\parallel=\left(- \frac{1}{1000000}\right) \left(1000 - \frac{1000}{e}\right)\simeq -0.00063212\,\mathrm N
$$

である。

層入口での単位長当たり吸収は

$$
S_a=2 \left(1000\right) \left(\frac{1}{2}\right)=1000\,\mathrm{W\,m^{-1}}
$$

である。

加熱電力は同じで力だけが反転する。
:::

:::example[垂直な入射]
強い吸収で平行力も大きいか判定する。

入力値を次のように置く。

$$
\begin{aligned}
P_{\rm in}&=1000\,\mathrm W \\k_i&=5\,\mathrm{m^{-1}} \\L&=1\,\mathrm m \\Z&=0\,\mathrm{s\,m^{-1}}
\end{aligned}
$$

これらを順に代入する。

層の光学厚さは

$$
\tau=2 \left(1\right) \left(5\right)=10
$$

である。

電力の透過割合は

$$
\mathcal T=e^{- \left(10\right)}\simeq 4.54\times 10^{-5}
$$

である。

透過電力は

$$
P_{\rm out}=\left(1000\right) \left(e^{-10}\right)\simeq 0.0454\,\mathrm W
$$

である。

層の吸収電力は

$$
P_{\rm abs}=\left(1000\right) - \left(\frac{1000}{e^{10}}\right)\simeq 999.95\,\mathrm W
$$

である。

吸収による平行力は

$$
F_\parallel=\left(0\right) \left(1000 - \frac{1000}{e^{10}}\right)=0\,\mathrm N
$$

である。

層入口での単位長当たり吸収は

$$
S_a=2 \left(1000\right) \left(5\right)=10000\,\mathrm{W\,m^{-1}}
$$

である。

この配置では大きな加熱があっても平行運動量を供給しない。
:::


## 演習問題

::::exercise{level=basic}
$R=5,L=1$ の $S,D$ を求めよ。

:::solution
半和と半差から $S=3,D=2$ である。
:::
::::

::::exercise{level=basic}
$R=L=P=2$ で全ての方向の横波の屈折率を求めよ。

:::solution
$\mathsf K=2I$ なので、横成分の係数は $2-N^2$ である。正の屈折率は $\sqrt2$ である。
:::
::::

::::exercise{level=basic}
$S=1,P=-1$ の共鳴角を求めよ。

:::solution
$A=\sin^2\theta-\cos^2\theta=0$ から、第 1 象限で $\theta=\pi/4$ である。
:::
::::

::::exercise{level=basic}
$k_i L=2$ の強度減衰率を求めよ。

:::solution
振幅が $e^{-2}$ なのでその二乗 $e^{-4}$ である。
:::
::::

::::exercise{level=advanced}
$P=0$、斜め伝播で $S\ne0$ の根を求めよ。

:::solution
$\sin\theta\ne0$ のもとで分散式は $\sin^2\theta X(SX-RL)=0$ であり、$X=0,RL/S$ である。
:::
::::

::::exercise{level=advanced}
同じ実 $\omega$ で $k(x)$ が二倍の場所に達した振幅を求めよ。

:::solution
反射と損失がなく単一の WKB 枝なら $k\abs a^2$ が一定なので $1/\sqrt2$ 倍である。
:::
::::

::::exercise{level=advanced}
Drude 模型で $\nu=\omega$ と $\nu=2\omega$ の加熱率を比較せよ。

:::solution
固定した $E,n,q,m$ で率は $\nu/(\nu^2+\omega^2)$ に比例する。後者は前者の $(2/5)/(1/2)=4/5$ 倍である。
:::
::::

::::exercise{level=advanced}
斜め伝播の波動行列で $E_y$ と $E_x$ の関係を求めよ。

:::solution
第 2 行は $iDE_x+(S-X)E_y=0$ である。$S\ne X$ なら $E_y=-iDE_x/(S-X)$ となる。$S=X$ ではこの行を割らずに使う。
:::
::::

## 確認テスト

::::quiz
冷たい無衝突模型の極から離れた実周波数で平均加熱はどうなるか。

- [x] ゼロ
- [ ] 正定値
- [ ] 周波数とともに単調増加

:::explanation
$\mathsf K=\mathsf K^\dagger$ なので、加熱を与える反 Hermite 部分はゼロである。
:::
::::

::::quiz{answer="4" tol="1e-6"}
$R=6,L=2$ の $S$ を答えよ。

:::explanation
半和 $(6+2)/2=4$ である。
:::
::::

::::quiz
WKB が一般に破れる場所はどこか。

- [x] 単純転回点
- [ ] 実波数が一定の領域
- [ ] 局所波長より十分緩やかに変化する領域

:::explanation
単純転回点で $k^2=0$ となり、$\abs{k'}/\abs k^2$ が発散する。接続には Airy 方程式を使う。
:::
::::

::::quiz{answer="1.5" tol="1e-6"}
$S=2,R=3,L=1$ の垂直 X 波の $N^2$ を答えよ。

:::explanation
式 :eqref[eq-plasma-xmode] より $3/2$ である。
:::
::::

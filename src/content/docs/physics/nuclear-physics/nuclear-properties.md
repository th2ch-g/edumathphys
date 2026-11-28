---
title: 原子核の基本的性質
description: 核図表、電子散乱、結合エネルギー、質量公式、核スピンとアイソスピンを扱う。
notation:
  - symbol: 'Z,N,A'
    meaning: 陽子数、中性子数、質量数
  - symbol: 'F(q)'
    meaning: 電荷形状因子
  - symbol: 'B(A,Z)'
    meaning: 結合エネルギー
  - symbol: 'I^\pi'
    meaning: 核スピンとパリティ
  - symbol: 'T,T_z'
    meaning: アイソスピンと第三成分
---

## この章の目標

- 核種を核図表に位置づける。
- 電子散乱から電荷分布と半径を推定する。
- 結合エネルギーと質量公式から核の安定性を考える。
- 核スピン、電磁モーメント、アイソスピンを区別する。

## 原子核の構成と核図表

陽子数は元素を決め、中性子数は同じ元素の中で安定性や質量を変える。陽子と中性子を合わせて :term[核子]{en="nucleon" yomi="かくし"} という。

:::definition[核種]{#def-nuc-nuclide}
:term[核種]{en="nuclide" yomi="かくしゅ"} を核図表で指定するには、陽子数 $Z$ と中性子数 $N$ を用いる。元素記号を $X$ とすると ${}^{A}_{Z}X$ と書き、$A=Z+N$ を質量数という。同じ $(N,Z)$ に長寿命の励起状態がある場合は、${}^{A\mathrm m}X$ などと書いて区別する。
:::

核図表では通常、横軸に $N$、縦軸に $Z$ を取る。同じ行は同位体、同じ列は同中性子体、傾き $-1$ の直線上は同重体である。一つの点には複数の励起状態があり、寿命が長い励起状態を核異性体という。励起エネルギー $E^*$ は核の静止エネルギーを $E^*$ だけ増す。

:::example[核図表の移動]
${}^{14}\mathrm C$ は $(N,Z)=(8,6)$ にある。中性子が陽子に変わる $\beta^-$ 崩壊後の核は $(7,7)$ の ${}^{14}\mathrm N$ である。$\alpha$ 放出では $N,Z$ がそれぞれ 2 減る。
:::

:::example[半径と密度]
半径を $R=r_0A^{1/3}$ とすると一様球の数密度は $3/(4\pi r_0^3)$ になる。$r_0=1.2\,\mathrm{fm}$ では約 $0.14\,\mathrm{fm}^{-3}$ である。質量中心から測った半径と二乗平均平方根半径は区別する。等質量の二核子が距離 $d$ だけ離れた配置では各核子の中心からの距離は $d/2$ である。
:::

本章では $1\,\mathrm{fm}=10^{-15}\,\mathrm m$、$\hbar c\simeq197.327\,\mathrm{MeV\,fm}$、$1\,\mathrm u\,c^2\simeq931.494\,\mathrm{MeV}$ を用いる。面積の単位 barn は $1\,\mathrm b=100\,\mathrm{fm}^2$ である。

### 核半径に関する三つの尺度

核の広がりは試験方法によって違う重みで測られる。点核子の座標 $\vb r_i$ を質量中心から測るとして、物質半径の二乗を $\frac1A\sum_i\langle r_i^2\rangle$、点陽子半径を $\frac1Z\sum_{p}\langle r_p^2\rangle$ とする。後者は電荷半径に近いが、核子自体の電荷が点に集中しているとした値である。電子散乱で使う核全体の電荷には、陽子の有限サイズや中性子の符号が変わる電荷分布も入る。どの半径を比較しているかを先に揃える。

:::example[中性子の外側への広がり]
点陽子の二乗半径が $9\,\mathrm{fm}^2$、点中性子が $12\,\mathrm{fm}^2$ で、$N/Z=3/2$ とする。物質半径の二乗は $(2/5)9+(3/5)12=10.8\,\mathrm{fm}^2$ で、物質半径と点陽子半径の差は約 $0.29\,\mathrm{fm}$ である。中性子半径と陽子半径の差そのものは約 $0.46\,\mathrm{fm}$ になる。
:::

### 核図表上の移動と平均密度

核反応で電子や光子も出入りする場合、核の陽子数と系全体の電荷を分けて数える。例えば中性子捕獲に続く $\beta^-$ 崩壊では、捕獲で $(N,Z)$ が $(N+1,Z)$ へ移り、崩壊で $(N,Z+1)$ へ移る。電子の電荷が終核の電荷増分を打ち消す。電荷についての同じ収支は、核図表の縦横を入れ替えて描いた場合にも成り立つ。

一様な核子数密度を $\rho_{\rm b}$ とする球では、粒子数は球の体積を使って $A=4\pi R^3\rho_{\rm b}/3$ である。したがって一定密度という仮定が $R=r_0A^{1/3}$ を与え、逆に半径の系統性から $\rho_{\rm b}=3/(4\pi r_0^3)$ が得られる。表面に有限な厚みがあると表面近くの核子の割合は概ね $R^{-1}$ に比例する。このため、ごく軽い核の半径を一定密度の球で再現する精度は重い核より低い。

:::example[核子数密度]
均一球の体積半径を $R=r_0A^{1/3}$ とする。核種を変えても同じ $r_0$ を使う近似での数密度を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
r_{0}=\frac{6}{5}\,\mathrm{fm}
$$

この条件で求める量を ${\rho_{\rm b}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\rho_{\rm b}}&=\frac{3}{4 \cdot \pi \cdot r_{0}^{3}}\\&=\frac{3}{4 \cdot \pi \cdot \left(\frac{6}{5}\right)^{3}}\\&\simeq 0.138155\,\mathrm{fm}^{-3}.\end{aligned}
$$

質量数は半径の三乗とともに消える。点陽子の密度が必要なら、陽子と中性子の形が同じ場合に限り $Z/A$ を掛ける。

半径の係数を少し大きく見積もった場合の誤差を調べる。そこで次の値に変更し、ほかの入力を固定する。
$$
r_{0}=\frac{5}{4}\,\mathrm{fm}
$$

このときの値を ${\rho_{\rm b}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\rho_{\rm b}}_{\mathrm{new}}&=\frac{3}{4 \cdot \pi \cdot \left(\frac{5}{4}\right)^{3}}\\&\simeq 0.122231\,\mathrm{fm}^{-3},\\{\rho_{\rm b}}_{\mathrm{new}}-{\rho_{\rm b}}&\simeq -0.0159243\,\mathrm{fm}^{-3}.\end{aligned}
$$
:::

:::example[投射面積の質量数依存]
同じ密度をもつ球の断面を考える。面積 $\pi R^2$ は寸法を表す幾何学的な尺度である。

量の定義に従い、共通の単位で次の値を用いる。
$$
r_{0}=\frac{6}{5}\,\mathrm{fm},\qquad A=27
$$

この条件で求める量を ${\mathcal A}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\mathcal A}&=\pi \cdot A^{\frac{2}{3}} \cdot r_{0}^{2}\\&=\pi \cdot \left(27\right)^{\frac{2}{3}} \cdot \left(\frac{6}{5}\right)^{2}\\&\simeq 40.715\,\mathrm{fm}^2.\end{aligned}
$$

核子数を増すと断面は増えるが、体積と異なり核子数に比例しない。反応断面積にはさらに相互作用と入射エネルギーが関係する。

別の質量数の球と比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
A=64
$$

このときの値を ${\mathcal A}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\mathcal A}_{\mathrm{new}}&=\pi \cdot \left(64\right)^{\frac{2}{3}} \cdot \left(\frac{6}{5}\right)^{2}\\&\simeq 72.3823\,\mathrm{fm}^2,\\{\mathcal A}_{\mathrm{new}}-{\mathcal A}&\simeq 31.6673\,\mathrm{fm}^2.\end{aligned}
$$
:::

:::example[表面層に入る核子]
一様球で外縁から厚さ $d<R$ 以内にある核子の割合を求める。中心から半径 $R-d$ までを内側と数える。

量の定義に従い、共通の単位で次の値を用いる。
$$
R=4\,\mathrm{fm},\qquad d=1\,\mathrm{fm}
$$

この条件で求める量を ${f_{\rm surf}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{f_{\rm surf}}&=1 - \left(\frac{R - d}{R}\right)^{3}\\&=1 - \frac{\left(- \left(1\right) + \left(4\right)\right)^{3}}{\left(4\right)^{3}}\\&\simeq 0.578125.\end{aligned}
$$

表面領域の割合が大きい場合には、$3d/R$ だけを用いる一次近似の誤差が無視できない。

同じ厚みを保ったまま球を大きくする。そこで次の値に変更し、ほかの入力を固定する。
$$
R=8\,\mathrm{fm}
$$

このときの値を ${f_{\rm surf}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{f_{\rm surf}}_{\mathrm{new}}&=1 - \frac{\left(- \left(1\right) + \left(8\right)\right)^{3}}{\left(8\right)^{3}}\\&\simeq 0.330078,\\{f_{\rm surf}}_{\mathrm{new}}-{f_{\rm surf}}&\simeq -0.248047.\end{aligned}
$$
:::

:::example[中性子過剰度]
中性子の数は $A-Z$ であり、陽子と中性子の数の差を全核子数で割ると無次元の非対称度になる。

量の定義に従い、共通の単位で次の値を用いる。
$$
Z=20,\qquad A=40
$$

この条件で求める量を ${\delta}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\delta}&=1 - \frac{2 \cdot Z}{A}\\&=- \frac{2 \cdot \left(20\right)}{\left(40\right)} + 1\\&=0.\end{aligned}
$$

基準となる $N=Z$ では対称性エネルギーの二次項は消える。質量数の比を使う場合、この量の符号は中性子過剰を正とする。

同じ元素の中性子の多い同位体を考える。そこで次の値に変更し、ほかの入力を固定する。
$$
A=48
$$

このときの値を ${\delta}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\delta}_{\mathrm{new}}&=- \frac{2 \cdot \left(20\right)}{\left(48\right)} + 1\\&\simeq 0.166667,\\{\delta}_{\mathrm{new}}-{\delta}&\simeq 0.166667.\end{aligned}
$$
:::

## 電子散乱と電荷分布・形状因子

### 散乱振幅の比

核の各部分から出る散乱波には位相差がある。この差を利用し、電荷密度を推定する。静的な核に対し電子の入射波と出射波を平面波とし、電気的相互作用の一次までを残す。核反跳と磁気散乱は省く。電子のポテンシャルエネルギー $V$ は、電荷密度 $\rho$ の Coulomb 核との畳み込みに $-e/(4\pi\varepsilon_0)$ を掛けたものである。行列要素の積分で $\vb s=\vb x-\vb r$ と置くと、位相 $e^{i\vb q\cdot\vb x}$ のうち $\vb s$ に依存する部分は点電荷の積分と共通になる。ここで $\vb q$ は入射波と出射波の波数ベクトルの差である。よって :term[電荷形状因子]{en="charge form factor" yomi="でんかけいじょういんし"} は
$$
F(\vb q)=\frac{1}{Ze}\int \rho(\vb r)e^{i\vb q\cdot\vb r}\dd^3r,\qquad
\left(\frac{\dd\sigma}{\dd\Omega}\right)_{\rm el}
=|F(\vb q)|^2\left(\frac{\dd\sigma}{\dd\Omega}\right)_{\rm point}.
\label{eq-nuc-electron-form-factor}
$$
点電荷の相対論的な断面積には電子スピンを含む Mott 断面積を採る。反跳が無視できると $q=2k\sin(\theta/2)$ である。重い核では電子波の Coulomb 歪みが増す。また、電荷密度には核内の陽子の位置分布に加えて核子自体の有限な電荷分布も含まれる。

### 球対称密度

$\vb q$ を極軸にして $\mu=\cos\theta_r$ と置くと $\int_{-1}^1 e^{iqr\mu}\dd\mu=2\sin(qr)/(qr)$ である。$|\rho|$ の四次モーメントが有限なら、Taylor 展開の剰余を $|qr|^4/120$ で抑えて積分できる。したがって
$$
\begin{aligned}
F(q)&=\frac{4\pi}{Ze}\int_0^\infty r^2\rho(r)\frac{\sin(qr)}{qr}\dd r\\
&=1-\frac{q^2}{6}\langle r^2\rangle_{\rm ch}+O(q^4),\qquad
\langle r^2\rangle_{\rm ch}=\frac{4\pi}{Ze}\int_0^\infty r^4\rho(r)\dd r .
\end{aligned}
\label{eq-nuc-charge-radius}
$$
断面積比 $|F|^2$ の $q^2=0$ での勾配は $-\langle r^2\rangle_{\rm ch}/3$ であり、振幅の勾配の 2 倍である。

:::example[一様球]
$\rho=3Ze/(4\pi R^3)$ を $r\le R$ で用いる。$\int_0^R r\sin(qr)\dd r=-R\cos(qR)/q+\sin(qR)/q^2$ なので、$x=qR$ に対し $F=3(\sin x-x\cos x)/x^3$ となる。最初の零点は約 $x=4.493$ である。$R=3\,\mathrm{fm}$ なら $q\simeq1.50\,\mathrm{fm}^{-1}$ となる。積分で求めた二乗半径は $3R^2/5$ である。
:::

:::example[Gaussian 分布]
$\rho=Ze(\pi a^2)^{-3/2}\exp(-r^2/a^2)$ では Fourier 積分が直交座標ごとに分離する。各座標の変換 $\exp(-q_i^2a^2/4)$ を掛けて $F=\exp(-q^2a^2/4)$ が得られる。$a=\sqrt{2/5}R$ なら一様球と同じ二乗半径だが、この形状因子は正のままである。
:::

### 電子スピンの影響と電荷分布の復元

静止点核による散乱で、電子の全エネルギーを $E$、運動量の大きさを $p$ とする。Dirac 平面波の上成分を二成分スピノル $\chi$、下成分を $h\vb\sigma\cdot\hat{\vb p}\,\chi$ とし、$h=pc/(E+m_ec^2)$ と置く。両波の規格化をそろえると電気ポテンシャルのスピン空間の行列は $(1+h^2\vb\sigma\cdot\hat{\vb p}_f\,\vb\sigma\cdot\hat{\vb p}_i)/(1+h^2)$ である。Pauli 行列の積を内積と外積に分け、初期スピン平均と終スピン和を取ると
$$
\frac{(1+h^2\cos\theta)^2+h^4\sin^2\theta}{(1+h^2)^2}
=1-\frac{4h^2}{(1+h^2)^2}\sin^2(\theta/2)
$$
となる。$2h/(1+h^2)=pc/E=v/c$ なので、これが Mott 断面積のスピン補正因子である。点電荷の Fourier 振幅が $q^{-2}$ であり、状態密度と流束により運動学因子が $p^2/v^2$ に比例する。従って補正因子を除く断面積は $[Ze^2/(8\pi\varepsilon_0 pv)]^2/\sin^4(\theta/2)$ である。

球対称で振幅の符号が分かり、逆変換が収束するなら
$$
\rho(r)=\frac{Ze}{2\pi^2}\int_0^\infty q^2 F(q)\frac{\sin(qr)}{qr}\dd q
$$
となる。三次元の逆変換で極軸を $\vb r$ に取った結果である。測定は $|F|^2$ を返すため、振幅の符号は $F(0)=1$、零点の性質、密度の模型などから推定する。

:::example[有限の波数窓]
点電荷の $F=1$ を $0\le q\le Q$ だけで逆変換すると、積分は $Ze[\sin(Qr)-Qr\cos(Qr)]/(2\pi^2r^3)$ となる。点状の電荷でも、有限範囲のデータから得た密度は幅 $1/Q$ の桁に広がり、外側で振動する。
:::

### 核子サイズと中心の形の分離

点陽子密度を $\rho_p(\vb r)$、単一陽子の中心に対する電荷密度を $f_p(\vb s)$ とする。$\int\rho_p=Z,\int f_p=e$ とし、核内で陽子の電荷分布が自由陽子と同じであると仮定する。陽子由来の電荷密度は $\rho_{\rm ch}=\rho_p*f_p$ である。Fourier 積分で $\vb x=\vb r+\vb s$ と置くと $F_{\rm ch}=F_p G_E^p$ となる。$F_p$ は中心分布の $Z$ で規格化した形状因子、$G_E^p$ は陽子の $e$ で規格化した形状因子である。中心のまわりで一次モーメントをゼロに取ると、
$$
\langle r^2\rangle_{\rm ch}=\langle r^2\rangle_{\text{point }p}+\langle r^2\rangle_p
$$
となる。これは $q^2$ の展開を比較しても、座標の二乗の交差項の平均がゼロになることからも得られる。

中性子の電荷の積分はゼロなので、$e$ を単位として規格化した一中性子の形状因子 $G_E^n$ は $q=0$ でゼロである。その補正が $-q^2\langle r^2\rangle_n/6$ から始まるなら、核全体の $\langle r^2\rangle_{\rm ch}$ には $(N/Z)\langle r^2\rangle_n$ が加わる。この中性子の二次モーメントは全電荷 1 の確率分布の分散とは違い、負にもなり得る。スピン電流や相対論補正を省く畳み込み模型の結果である。

:::example[二次モーメントの補正]
点陽子半径 $3\,\mathrm{fm}$、自由陽子に割り当てた半径 $0.84\,\mathrm{fm}$ の模型で中性子の電荷を省くと、二乗半径は $9+0.7056=9.7056\,\mathrm{fm}^2$、半径は約 $3.12\,\mathrm{fm}$ である。二つの長さを直接足すと過大になる。
:::

### 角度と分解能

相対論的な電子で運動エネルギーが静止エネルギーより十分に大きければ $pc\simeq E$ である。大角度の散乱ほど $\hbar q=2p\sin(\theta/2)$ が大きいが、小角度に比べて点核の断面積は小さくなる。また $|F|^2$ の回折最小ではイベント数が減る。電荷半径だけを求めるなら十分低い $q$ の点が重要であり、表面の厚さまで求めるなら零点より高い $q$ の精度も必要になる。

:::example[形状因子の曲率]
一様球の $F$ は $1-x^2/10+x^4/280+O(x^6)$、初期勾配が同じ Gaussian は $e^{-x^2/10}=1-x^2/10+x^4/200+O(x^6)$ である。$x=qR$ とした。二次項が同じでも、四次項の比は $7/5$ である。低角のデータのうち誤差がこの差を上回る範囲だけでは、二つの密度を判別しにくい。
:::

### 表面の厚さを含む電荷分布

安定な中重核の球対称成分を当てはめる模型として、:term[二母数 Fermi 分布]{en="two-parameter Fermi distribution" yomi="にぼすうふぇるみぶんぷ"} $\rho(r)=\rho_0/[1+\exp((r-R_{1/2})/a)]$ を用いる。$R_{1/2}$ は密度が内側の外挿値 $\rho_0$ の半分となる半径、$a>0$ は拡散度である。$\rho/\rho_0=f$ を解くと $r=R_{1/2}+a\ln[(1-f)/f]$ であり、$f=0.9$ と $0.1$ の間の厚さは $2a\ln9$ である。

$a\ll R_{1/2}$ とし、$I_j=\int_0^\infty r^j\rho(r)\dd r/\rho_0$ と置く。鋭い境界との差を $r=R_{1/2}\pm ax$ で内外に分けると、$a\int_0^\infty[(R_{1/2}+ax)^j-(R_{1/2}-ax)^j]/(1+e^x)\dd x$ が得られる。下側の積分を $x>R_{1/2}/a$ まで延長した誤差は指数的に小さい。$\int_0^\infty x/(1+e^x)\dd x=\pi^2/12$ は分母の逆数を指数関数の交代級数に展開して積分すれば得られる。従って
$$
I_j=\frac{R_{1/2}^{j+1}}{j+1}+
\frac{\pi^2a^2}{6}jR_{1/2}^{j-1}+O(a^4R_{1/2}^{j-3}).
$$
$4\pi\rho_0 I_2=Ze$ で規格化すると、二乗半径は $I_4/I_2=3R_{1/2}^2/5+7\pi^2a^2/5+O(a^4/R_{1/2}^2)$ である。同じ二乗半径は $R_{1/2}$ と $a$ の複数の組で得られるため、表面の厚さまで決めるには形状因子の高い $q$ の情報を使う。電子散乱による半径と密度の測定については :cite[krane-nuclear] を参照する。

:::example[拡散度の半径への寄与]
$R_{1/2}=5\,\mathrm{fm}$、$a=0.5\,\mathrm{fm}$ とすると、密度が $90\%$ から $10\%$ まで低下する幅は約 $2.20\,\mathrm{fm}$ である。二乗半径は拡散度の二次までで $15+7\pi^2/20\simeq18.45\,\mathrm{fm}^2$ となり、二乗平均平方根半径は約 $4.30\,\mathrm{fm}$ である。
:::

### 低運動量のモーメントと回折による判定

有限な六次モーメントをもつ球対称分布では、$\sin x/x$ の Taylor 展開とその余項を動径積分に入れられる。余項の絶対値は $|x|^6/7!$ 以下であり、$\rho\ge0$ の場合は積分した誤差を $\langle r^6\rangle |q|^6/7!$ 以下に抑えられる。従って
$$
F(q)=1-\frac{q^2}{6}\langle r^2\rangle+\frac{q^4}{120}\langle r^4\rangle+O(q^6).
$$
同じ二乗半径をもつ二分布は一次の勾配が等しい。:term[回折極小]{en="diffraction minimum" yomi="かいせつきょくしょう"} や四次項は、その勾配だけでは決まらない表面の情報を与える。一様球の場合は分子 $\sin(qR)-qR\cos(qR)$ の最初の正根を観測する。$\tan(qR)=qR$ を使うときは、$\cos(qR)=0$ の点では元の分子がゼロでないことも確認する。

:::example[Gaussian 分布の減衰]
各方向の座標分散を $b^2/2$ とする等方 Gaussian 電荷分布を仮定する。一次元の Fourier 積分を三方向で掛け合わせると振幅は $\exp[-(qb)^2/4]$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
q=\frac{1}{2}\,\mathrm{fm}^{-1},\qquad b=2\,\mathrm{fm}
$$

この条件で求める量を ${F}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{F}&=e^{\frac{\left(-1\right) \cdot \left(b \cdot q\right)^{2}}{4}}\\&=e^{- \frac{\left(2\right)^{2} \cdot \left(\frac{1}{2}\right)^{2}}{4}}\\&\simeq 0.778801.\end{aligned}
$$

同じ条件の散乱断面積の比は振幅の二乗である。この分布では有限な正の $q$ に回折零点をもたない。

より短い長さを分解する運動量で比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
q=1\,\mathrm{fm}^{-1}
$$

このときの値を ${F}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{F}_{\mathrm{new}}&=e^{- \frac{\left(1\right)^{2} \cdot \left(2\right)^{2}}{4}}\\&\simeq 0.367879,\\{F}_{\mathrm{new}}-{F}&\simeq -0.410921.\end{aligned}
$$
:::

:::example[回折零点からの半径]
一様球で $\tan x=x$ の最初の正根を $x_1=4.4934$ とする。計測した零点の波数 $q_1$ から境界半径を推定する。

量の定義に従い、共通の単位で次の値を用いる。
$$
x_{1}=\frac{22467}{5000},\qquad q_{1}=\frac{9}{10}\,\mathrm{fm}^{-1}
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=\frac{x_{1}}{q_{1}}\\&=\frac{\left(\frac{22467}{5000}\right)}{\left(\frac{9}{10}\right)}\\&\simeq 4.99267\,\mathrm{fm}.\end{aligned}
$$

この半径を二乗平均平方根半径に換算するには一様球の係数 $\sqrt{3/5}$ を掛ける。極小が深くないときは背景と分解能も確認する。

極小位置の読み取りを変更する。そこで次の値に変更し、ほかの入力を固定する。
$$
q_{1}=1\,\mathrm{fm}^{-1}
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\frac{\left(\frac{22467}{5000}\right)}{\left(1\right)}\\&\simeq 4.4934\,\mathrm{fm},\\{R}_{\mathrm{new}}-{R}&\simeq -0.499267\,\mathrm{fm}.\end{aligned}
$$
:::

:::example[二成分の二乗半径]
中心部と外側に広い部分の規格化分布を重ねる。規格化総電荷に対する外側の比を $w$ とし、二成分の二乗半径を $b_1^2,b_2^2$ とする。

量の定義に従い、共通の単位で次の値を用いる。
$$
w=\frac{1}{20},\qquad b_{1}=3\,\mathrm{fm},\qquad b_{2}=7\,\mathrm{fm}
$$

この条件で求める量を ${\langle r^2\rangle}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\langle r^2\rangle}&=b_{1}^{2} \cdot \left(1 - w\right) + b_{2}^{2} \cdot w\\&=\left(3\right)^{2} \cdot \left(1 - \left(\frac{1}{20}\right)\right) + \left(7\right)^{2} \cdot \left(\frac{1}{20}\right)\\&=11\,\mathrm{fm}^2.\end{aligned}
$$

少量の外側電荷でも二乗半径には強く効く。振幅自身が密度に線形であるため、この平均には二成分間の交差項が出ない。

外側の割合を増やして影響を見る。そこで次の値に変更し、ほかの入力を固定する。
$$
w=\frac{1}{10}
$$

このときの値を ${\langle r^2\rangle}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\langle r^2\rangle}_{\mathrm{new}}&=\left(3\right)^{2} \cdot \left(1 - \left(\frac{1}{10}\right)\right) + \left(7\right)^{2} \cdot \left(\frac{1}{10}\right)\\&=13\,\mathrm{fm}^2,\\{\langle r^2\rangle}_{\mathrm{new}}-{\langle r^2\rangle}&=2\,\mathrm{fm}^2.\end{aligned}
$$
:::

:::example[反跳なし近似での分解能]
電子の入射と出射の運動量の大きさをともに $p$ とし、散乱角を $\theta$ とする。二つの運動量ベクトルの差の大きさを $\hbar c$ で換算する。

量の定義に従い、共通の単位で次の値を用いる。
$$
p_{c}=200\,\mathrm{MeV},\qquad \theta=\frac{\pi}{6}\,\mathrm{rad},\qquad h_{c}=\frac{197327}{1000}\,\mathrm{MeV\,fm}
$$

この条件で求める量を ${q}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{q}&=\frac{2 \cdot p_{c} \cdot \sin{\left(\frac{\theta}{2} \right)}}{h_{c}}\\&=\frac{2 \cdot \left(200\right) \cdot \sin{\left(\frac{\left(\frac{\pi}{6}\right)}{2} \right)}}{\left(\frac{197327}{1000}\right)}\\&\simeq 0.52465\,\mathrm{fm}^{-1}.\end{aligned}
$$

この見積もりで散乱角を波数に換算できる。核の反跳が大きい場合は出射運動量を弾性散乱の保存則から求める必要がある。

ビーム運動量を固定し大角度側と比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
\theta=\frac{\pi}{2}\,\mathrm{rad}
$$

このときの値を ${q}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{q}_{\mathrm{new}}&=\frac{2 \cdot \left(200\right) \cdot \sin{\left(\frac{\left(\frac{\pi}{2}\right)}{2} \right)}}{\left(\frac{197327}{1000}\right)}\\&\simeq 1.43337\,\mathrm{fm}^{-1},\\{q}_{\mathrm{new}}-{q}&\simeq 0.908721\,\mathrm{fm}^{-1}.\end{aligned}
$$
:::

### 反跳から運動量移行を求める

電子質量を無視し、標的核は静止しているとする。初期の核四元運動量を $P$、電子の初期と終状態を $p_i,p_f$ とする。核は励起しないので $(P+p_i-p_f)^2=P^2$ であり、$p_i^2=p_f^2=0$ である。計量 $(+,-,-,-)$ を選べば
$$
Mc^2(E_i-E_f)=E_iE_f(1-\cos\theta),\qquad
E_f=\frac{E_i}{1+(E_i/Mc^2)(1-\cos\theta)}
$$
を得る。核の空間運動量の大きさは $\sqrt{E_i^2+E_f^2-2E_iE_f\cos\theta}/c$ である。反跳を含む相対論的な四元移行の二乗の符号を反転した量は $2E_iE_f(1-\cos\theta)/c^2$ であり、空間移行の二乗とはエネルギー移行の項だけ違う。静的な形状因子の導出では空間波数を $q$ としていたので、測定から換算した量の定義をそろえる必要がある。

:::example[重い核と軽い核の反跳の比較]
$E_i=200\,\mathrm{MeV}$、$\theta=\pi/2$ のとき、核の静止エネルギーが $4\,\mathrm{GeV}$ なら $E_f=200/(1+0.05)\simeq190.5\,\mathrm{MeV}$ である。$40\,\mathrm{GeV}$ なら $E_f=200/(1+0.005)\simeq199.0\,\mathrm{MeV}$ である。同じ角度でも軽い核の方が弾性ピークの位置は低い。この換算は電子の運動学だけであり、断面積の Coulomb 歪みを補正したものではない。
:::

## 質量と結合エネルギー

:::definition[結合エネルギー]{#def-nuc-binding}
:term[結合エネルギー]{en="binding energy" yomi="けつごうえねるぎー"} は、自由核子と基底状態の核の静止エネルギー差
$$
B(A,Z)=[Zm_p+(A-Z)m_n-M_{\rm nuc}(A,Z)]c^2
\label{eq-nuc-binding}
$$
である。$m_p,m_n$ は陽子と中性子の質量である。
:::

中性原子の質量を $M_{\rm at}$、正の全電子結合エネルギーを $B_e$ とすると $M_{\rm at}c^2=M_{\rm nuc}c^2+Zm_ec^2-B_e$ である。水素原子質量 $m_{\rm H}$ と結合エネルギー $B_{\rm H}$ を使えば、$B=[Zm_{\rm H}+Nm_n-M_{\rm at}]c^2+ZB_{\rm H}-B_e$ となる。核の MeV スケールの概算では最後の差を省けることが多い。電子の個数が反応前後で等しいときは電子の静止質量が相殺する。

核を中性子と残留核へ分ける最小エネルギーは生成物と親核の質量差である。:eqref[eq-nuc-binding] を差し引くと
$$
S_n=B(A,Z)-B(A-1,Z),\qquad S_p=B(A,Z)-B(A-1,Z-1).
\label{eq-nuc-separation}
$$
2 中性子分離エネルギーは隣接する $S_n$ の和である。$B/A$ は全体の平均であり、最後の核子の放出しきい値は分離エネルギーで判定する。

:::example[重陽子]
$m_{\rm H}=1.007825\,\mathrm u$、$m_n=1.008665\,\mathrm u$、$M_{\rm at}=2.014102\,\mathrm u$ から、電子結合差を省いた質量欠損は $0.002388\,\mathrm u$ となる。静止エネルギーに変換すると約 $2.224\,\mathrm{MeV}$ である。
:::

:::example[しきい値]
同じ $Z$ で $B(A-2)=70$、$B(A-1)=73$、$B(A)=79\,\mathrm{MeV}$ なら、$S_n=6$、$S_{2n}=9\,\mathrm{MeV}$ である。$7\,\mathrm{MeV}$ の励起状態からは 1 中性子チャネルが開き、2 中性子チャネルは閉じる。
:::

### 質量超過と誤差

表に $D(A,Z)=[M_{\rm at}(A,Z)-A\,\mathrm u]c^2$ が載ることがある。これは :term[質量超過]{en="mass excess" yomi="しつりょうちょうか"} である。$A$ を保つ反応では $A\,\mathrm u$ が相殺し、発熱量は質量超過の差で得られる。$D_n=(m_n-\mathrm u)c^2$ と置けば、原子質量による中性子分離エネルギーは $D(A-1,Z)+D_n-D(A,Z)$ である。電子束縛補正は要求する精度に応じて入れる。

:::example[質量差の標準偏差]
質量超過の独立な誤差が $3$ と $4\,\mathrm{keV}$ なら、差の分散は $3^2+4^2$ で、標準偏差は $5\,\mathrm{keV}$ である。両測定に同じオフセットが加わる場合、その寄与は差で相殺する。一般の共分散がある場合は分散の和から共分散の 2 倍を引く。
:::

### 質量差に現れる電子の結合と相関誤差

中性原子の全電子の結合エネルギーを正の $B_e$ と書く。電子を無限遠に分けた状態の静止エネルギーは $M_{\rm at}c^2+B_e$ なので、核質量は $M_{\rm at}-Zm_e+B_e/c^2$ である。核反応の両辺で $Z$ の総和が等しいと電子の静止質量は消えるが、$B_e$ は原子ごとに違う。したがって電子の束縛を無視できるかどうかは、質量表の桁数よりも最終的に必要な反応エネルギーの精度で判断する。

差 $D=X-Y$ の分散は、中心化した変数の差を二乗して平均することで $\sigma_D^2=\sigma_X^2+\sigma_Y^2-2\operatorname{Cov}(X,Y)$ になる。同一の参照質量から較正した質量を比較すると共有する誤差が相殺しうる。

:::example[原子と核の質量差]
質量エネルギーの単位を MeV とする。核に含まれない電子の静止エネルギーを引き、電子の結合の分を戻す補正を計算する。

量の定義に従い、共通の単位で次の値を用いる。
$$
Z=20,\qquad m_{ec2}=\frac{511}{1000}\,\mathrm{MeV},\qquad B_{e}=\frac{1}{50}\,\mathrm{MeV}
$$

この条件で求める量を ${\Delta_{\rm at\to nuc}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta_{\rm at\to nuc}}&=B_{e} + - Z \cdot m_{ec2}\\&=- \left(20\right) \cdot \left(\frac{511}{1000}\right) + \left(\frac{1}{50}\right)\\&\simeq -10.2\,\mathrm{MeV}.\end{aligned}
$$

電子一個の静止質量は電子の結合エネルギーとは尺度が異なる。補正の符号は、原子から裸の核へ変換する向きに対応する。

電子結合エネルギーの近似値を変える。そこで次の値に変更し、ほかの入力を固定する。
$$
B_{e}=\frac{1}{40}\,\mathrm{MeV}
$$

このときの値を ${\Delta_{\rm at\to nuc}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta_{\rm at\to nuc}}_{\mathrm{new}}&=- \left(20\right) \cdot \left(\frac{511}{1000}\right) + \left(\frac{1}{40}\right)\\&\simeq -10.195\,\mathrm{MeV},\\{\Delta_{\rm at\to nuc}}_{\mathrm{new}}-{\Delta_{\rm at\to nuc}}&\simeq 0.005\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[二中性子分離と閾値]
親核と二つ中性子の少ない核の結合エネルギーをそれぞれ $B_A,B_{A2}$ と書く。連続状態の二中性子の運動エネルギーは閾値でゼロとする。

量の定義に従い、共通の単位で次の値を用いる。
$$
B_{A}=120,\qquad B_{A2}=108
$$

この条件で求める量を ${S_{2n}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{S_{2n}}&=B_{A} - B_{A2}\\&=- \left(108\right) + \left(120\right)\\&=12\,\mathrm{MeV}.\end{aligned}
$$

これは二中性子を最終状態に置く閾値である。一中性子放出が先に進むかは、中間の核の質量も調べて判断する。

二中性子を失った核の結合が強い模型と比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
B_{A2}=111
$$

このときの値を ${S_{2n}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{S_{2n}}_{\mathrm{new}}&=- \left(111\right) + \left(120\right)\\&=9\,\mathrm{MeV},\\{S_{2n}}_{\mathrm{new}}-{S_{2n}}&=-3\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[同じ基準の質量差]
二つの質量超過に共通標準偏差 $s$ の加法誤差があり、各々に独立な標準偏差 $u,v$ の揺らぎがあるとする。差の標準偏差を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
s=4\,\mathrm{keV},\qquad u=3\,\mathrm{keV},\qquad v=4\,\mathrm{keV}
$$

この条件で求める量を ${\sigma_D}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\sigma_D}&=\sqrt{- 2 \cdot s^{2} + s^{2} + s^{2} + u^{2} + v^{2}}\\&=\sqrt{\left(3\right)^{2} + \left(4\right)^{2}}\\&=5\,\mathrm{keV}.\end{aligned}
$$

加法誤差が完全に共通なら、差の精度を決めるのは二つの独立な成分である。

共通誤差だけの大きさを変更する。そこで次の値に変更し、ほかの入力を固定する。
$$
s=7\,\mathrm{keV}
$$

このときの値を ${\sigma_D}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\sigma_D}_{\mathrm{new}}&=\sqrt{\left(3\right)^{2} + \left(4\right)^{2}}\\&=5\,\mathrm{keV},\\{\sigma_D}_{\mathrm{new}}-{\sigma_D}&=0\,\mathrm{keV}.\end{aligned}
$$
:::

:::example[電荷を保存する捕獲の Q 値]
核質量の代わりに中性原子の質量を使う。中性子捕獲で質量数も電荷も全体で保存するので、MeV 単位の質量超過 $\Delta_t,\Delta_n,\Delta_f$ の差を取ればよい。電子の結合差は無視する。

量の定義に従い、共通の単位で次の値を用いる。
$$
\Delta_{t}=-20\,\mathrm{MeV},\qquad \Delta_{n}=\frac{8071}{1000}\,\mathrm{MeV},\qquad \Delta_{f}=-18\,\mathrm{MeV}
$$

この条件で求める量を ${Q}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{Q}&=- \Delta_{f} + \Delta_{n} + \Delta_{t}\\&=- \left(-18\right) + \left(-20\right) + \left(\frac{8071}{1000}\right)\\&\simeq 6.071\,\mathrm{MeV}.\end{aligned}
$$

質量超過を使うと双方に共通する原子質量単位の整数倍が消える。負の質量超過も符号ごと代入する。

終核の質量評価値を変えた場合を調べる。そこで次の値に変更し、ほかの入力を固定する。
$$
\Delta_{f}=-17\,\mathrm{MeV}
$$

このときの値を ${Q}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{Q}_{\mathrm{new}}&=- \left(-17\right) + \left(-20\right) + \left(\frac{8071}{1000}\right)\\&\simeq 5.071\,\mathrm{MeV},\\{Q}_{\mathrm{new}}-{Q}&=-1\,\mathrm{MeV}.\end{aligned}
$$
:::

## Weizsäcker の質量公式と β 安定線

体積内の核子が近傍と主に相互作用するなら、平均的な結合はほぼ $A$ に比例する。表面では近傍が欠けるので面積 $A^{2/3}$ に比例する減少が生じる。この見方で
$$
B=a_v A-a_s A^{2/3}-a_c Z(Z-1)A^{-1/3}
 -a_a(N-Z)^2/A+\delta
\label{eq-nuc-mass-formula}
$$
と表す。一様な電荷球を内側から組み立てる仕事は、半径 $r$ に既にある電荷 $Ze(r/R)^3$ と追加する電荷 $3Ze r^2\dd r/R^3$ の積を $4\pi\varepsilon_0r$ で割って積分したものである。結果は $3(Ze)^2/(20\pi\varepsilon_0R)$ となる。実際には陽子の自己相互作用を省くため $Z^2$ を $Z(Z-1)$ に置き換える。

対称項の符号は Fermi 運動に由来する部分から確認できる。固定体積で各核子種の運動エネルギーは粒子数の $5/3$ 乗に比例する。$\eta=(N-Z)/A$ に対して $[(1+\eta)^{5/3}+(1-\eta)^{5/3}]/2=1+5\eta^2/9+O(\eta^4)$ となり、不均衡で運動エネルギーが増す。相互作用も $a_a$ に寄与する。代表的な対項は偶偶核で $+a_p/\sqrt A$、奇奇核でその負、奇数 $A$ でゼロである。

固定 $A$ では原子の質量エネルギーは $Am_nc^2-\Delta Z-B$、$\Delta=(m_n-m_{\rm H})c^2$ と近似できる。$\delta$ を省くと $Z$ に関する二次係数が $C=a_c A^{-1/3}+4a_a/A>0$、一次係数が $-(4a_a+a_cA^{-1/3}+\Delta)$ である。平方完成による連続最小は
$$
Z_*=\frac{4a_a+a_cA^{-1/3}+\Delta}{2C}.
\label{eq-nuc-beta-stability}
$$
実際には整数 $Z$ の質量差を取る。特に偶数 $A$ では偶偶と奇奇で対項が違い、二本の質量放物線になる。

:::example[連続最小]
$A=125$、$a_c=0.71\,\mathrm{MeV}$、$a_a=23\,\mathrm{MeV}$、$\Delta=0.782\,\mathrm{MeV}$ とすると $C=0.878\,\mathrm{MeV}$、$Z_*\simeq52.92$ となる。この模型では最小整数は 53 であり、$Z=52$ と $53$ の差は約 $0.74\,\mathrm{MeV}$ である。個々の崩壊の判定には実測質量を用いる。
:::

:::example[中性子の超過]
$\Delta$ と $Z(Z-1)$ の線形補正を省くと $Z_*/A=[2+(a_c/2a_a)A^{2/3}]^{-1}$ である。例えば $A^{2/3}=25$ なら上の係数で $Z_*/A\simeq0.419$ である。重い核ほど陽子間の反発を減らす効果が重要になる。
:::

### 偶数質量数での二重ベータの可能性

偶偶核の質量エネルギーを $F(Z)-p$、奇奇核を $F(Z)+p$ と書く。$F(Z)=M_*c^2+C(Z-Z_*)^2$、$p=a_p/\sqrt A>0$ とする。偶偶の $Z$ から奇奇の $Z+1$ への差は $Q_1=C[-2(Z-Z_*)-1]-2p$、偶偶の $Z+2$ への差は $Q_2=C[-4(Z-Z_*)-4]$ となる。$-1/2-p/C\le Z-Z_*<-1$ なら $Q_1\le0$ かつ $Q_2>0$ である。この区間が存在する条件は $p/C>1/2$ であり、さらに区間内に偶数 $Z$ が必要である。

:::example[崩壊窓]
$C=1\,\mathrm{MeV},p=1\,\mathrm{MeV}$、$Z-Z_*=-5/4$ の模型では $Q_1=-0.5\,\mathrm{MeV}$、$Q_2=1\,\mathrm{MeV}$ である。一段の崩壊は閉じるが、二段の電荷差をもつ終状態へはエネルギーを放出できる。寿命には [弱い遷移](/physics/nuclear-physics/beta-decay/) の二次の振幅が関わる。
:::

### β 安定線の推定を反応に使う

対称性項は核子の陽子・中性子の偏りを抑える。Coulomb 項は電荷が増えるとより強く結合を弱めるので、両者を合わせた質量の極小は一般に $N=Z$ からずれる。頂点の位置を実験と比較するときは、まず奇数 $A$ の組を用いると陽子・中性子の対エネルギーの補正が小さい。

一定 $A$ で中性原子質量を使い、終状態の原子の電離・励起とニュートリノ質量を無視すると、$\beta^-$ のエネルギー窓は $[M_{\rm at}(A,Z)-M_{\rm at}(A,Z+1)]c^2-E^*$ である。ここで $E^*$ は終核の励起エネルギーである。$\beta^+$ と電子捕獲では同じ $Z-1$ の終核を作るが、陽電子放出では中性原子の質量差からさらに $2m_ec^2$ を引く。陽電子自身に加え、親原子の余分な一電子の静止エネルギーも終側に数えるためである。電子捕獲ではこの二電子分が差し引かれず、代わりに捕獲後の空孔などの励起がニュートリノに渡るエネルギーを減らす。それぞれについて窓が正であることは運動学的な必要条件であり、崩壊率には遷移行列要素が入る。詳しくは [β 崩壊と弱い相互作用](/physics/nuclear-physics/beta-decay/) を参照する。

:::example[連続近似の安定線]
核子の質量差と電子質量を無視し、Coulomb 項では $Z(Z-1)$ を $Z^2$ に置く。固定 $A$ で $a_c Z^2/A^{1/3}+a_a(A-2Z)^2/A$ を微分して得た頂点を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
A=27,\qquad a_{c}=\frac{71}{100}\,\mathrm{MeV},\qquad a_{a}=23\,\mathrm{MeV}
$$

この条件で求める量を ${Z_*}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{Z_*}&=\frac{A}{\frac{A^{\frac{2}{3}} \cdot a_{c}}{2 \cdot a_{a}} + 2}\\&=\frac{\left(27\right)}{2 + \frac{\left(27\right)^{\frac{2}{3}} \cdot \left(\frac{71}{100}\right)}{2 \cdot \left(23\right)}}\\&\simeq 12.6232.\end{aligned}
$$

頂点は実数である。安定核を選ぶ段階では整数化し、偶数 $A$ なら対項を入れて隣接核の質量を比較する。

軽い核と重い核の頂点を比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
A=125
$$

このときの値を ${Z_*}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{Z_*}_{\mathrm{new}}&=\frac{\left(125\right)}{\frac{\left(125\right)^{\frac{2}{3}} \cdot \left(\frac{71}{100}\right)}{2 \cdot \left(23\right)} + 2}\\&\simeq 52.3918,\\{Z_*}_{\mathrm{new}}-{Z_*}&\simeq 39.7686.\end{aligned}
$$
:::

:::example[原子質量と核子の質量差]
対項を外し、$Z(Z-1)$ を用いる。$d=(m_{\rm H}-m_n)c^2$ と定めると、質量エネルギーの $Z$ の一次係数には $d$ が加わる。

量の定義に従い、共通の単位で次の値を用いる。
$$
A=125,\qquad a_{c}=\frac{71}{100}\,\mathrm{MeV},\qquad a_{a}=23\,\mathrm{MeV},\qquad d=0\,\mathrm{MeV}
$$

この条件で求める量を ${Z_*}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{Z_*}&=\frac{4 \cdot a_{a} - d + \frac{a_{c}}{A^{1 \cdot \frac{1}{3}}}}{2 \cdot \left(\frac{a_{c}}{A^{1 \cdot \frac{1}{3}}} + \frac{4 \cdot a_{a}}{A}\right)}\\&=\frac{- \left(0\right) + 4 \cdot \left(23\right) + \frac{\left(\frac{71}{100}\right)}{\sqrt[3]{\left(125\right)}}}{\frac{8 \cdot \left(23\right)}{\left(125\right)} + \frac{2 \cdot \left(\frac{71}{100}\right)}{\sqrt[3]{\left(125\right)}}}\\&\simeq 52.4727.\end{aligned}
$$

中性水素の方が軽い $d<0$ の場合、質量差のみの効果は頂点を陽子数の大きい側に動かす。Coulomb 項がもたらすずれとは符号が異なる。

原子質量に合わせて $d$ を含める。そこで次の値に変更し、ほかの入力を固定する。
$$
d=- \frac{391}{500}\,\mathrm{MeV}
$$

このときの値を ${Z_*}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{Z_*}_{\mathrm{new}}&=\frac{- \left(- \frac{391}{500}\right) + 4 \cdot \left(23\right) + \frac{\left(\frac{71}{100}\right)}{\sqrt[3]{\left(125\right)}}}{\frac{8 \cdot \left(23\right)}{\left(125\right)} + \frac{2 \cdot \left(\frac{71}{100}\right)}{\sqrt[3]{\left(125\right)}}}\\&\simeq 52.918,\\{Z_*}_{\mathrm{new}}-{Z_*}&\simeq 0.44533.\end{aligned}
$$
:::

:::example[偶奇差の大きさ]
結合エネルギーの偶奇補正を偶偶で $+a_p A^{-1/2}$、奇奇で $-a_p A^{-1/2}$ とする。同じ $A$ で放物線の中心から等距離の二核を比べる仮想例である。

量の定義に従い、共通の単位で次の値を用いる。
$$
a_{p}=12\,\mathrm{MeV},\qquad A=36
$$

この条件で求める量を ${\Delta_{\rm eo}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta_{\rm eo}}&=\frac{2 \cdot a_{p}}{\sqrt{A}}\\&=\frac{2 \cdot \left(12\right)}{\sqrt{\left(36\right)}}\\&=4\,\mathrm{MeV}.\end{aligned}
$$

結合エネルギーの増加は質量を小さくする。この比較では奇奇核の質量が表示した量だけ上にある。

対補正が質量数とともに弱くなる程度を比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
A=144
$$

このときの値を ${\Delta_{\rm eo}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta_{\rm eo}}_{\mathrm{new}}&=\frac{2 \cdot \left(12\right)}{\sqrt{\left(144\right)}}\\&=2\,\mathrm{MeV},\\{\Delta_{\rm eo}}_{\mathrm{new}}-{\Delta_{\rm eo}}&=-2\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[励起準位への β の窓]
連続曲線の頂点から距離 $x=Z-Z_*$ の親核を置く。娘核は $x+1$ なので、エネルギー差からその励起エネルギー $E_f$ を引く。

量の定義に従い、共通の単位で次の値を用いる。
$$
C=\frac{6}{5}\,\mathrm{MeV},\qquad x=- \frac{3}{2},\qquad E_{f}=1\,\mathrm{MeV}
$$

この条件で求める量を ${Q_f}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{Q_f}&=C \cdot \left(x^{2} - \left(x + 1\right)^{2}\right) - E_{f}\\&=- \left(1\right) + \left(\frac{6}{5}\right) \cdot \left(\left(- \frac{3}{2}\right)^{2} - \left(\left(- \frac{3}{2}\right) + 1\right)^{2}\right)\\&\simeq 1.4\,\mathrm{MeV}.\end{aligned}
$$

負の値は当該終状態へ遷移するエネルギーが足りないことを示す。許された準位でも、選択則から分岐が弱い場合がある。

より高い励起に達する場合を比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
E_{f}=3\,\mathrm{MeV}
$$

このときの値を ${Q_f}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{Q_f}_{\mathrm{new}}&=- \left(3\right) + \left(\frac{6}{5}\right) \cdot \left(\left(- \frac{3}{2}\right)^{2} - \left(\left(- \frac{3}{2}\right) + 1\right)^{2}\right)\\&\simeq -0.6\,\mathrm{MeV},\\{Q_f}_{\mathrm{new}}-{Q_f}&=-2\,\mathrm{MeV}.\end{aligned}
$$
:::

## 核スピン・パリティと電磁モーメント

回転対称な核の状態を角運動量 $\vb I^2,I_z$ の固有値 $\hbar^2 I(I+1),\hbar M$ で分類する（:ref[thm-qm-ladder]）。:term[パリティ]{en="parity" yomi="ぱりてぃ"} $\pi$ は空間反転の固有値である。球面調和関数の反転は $(-1)^\ell$ なので、独立粒子配置の全パリティは $\prod_i(-1)^{\ell_i}$ になる。

核磁子を $\mu_N=e\hbar/(2m_p)$ とする。磁気双極子モーメント $\mu$ は $\ket{I,I}$ における $\hat\mu_z$ の期待値である。同じ $I>0$ 多重項内では階数 1 の演算子が角運動量に比例するので、磁場 $B_0\vb e_z$ の一次エネルギー移動は $-\mu(M/I)B_0$ となる。分光学的四重極モーメントは
$$
Q_s=\bra{I,I}\int(3z^2-r^2)\hat\rho(\vb r)\dd^3r\ket{I,I}
$$
であり、単位は電荷と面積の積である。対角行列要素が存在するには $I\otimes2$ に $I$ が含まれる必要があり、三角条件 $2\le 2I$ が得られる。そのため $I=0,1/2$ では $Q_s$ が消える。核の固有座標における形状は遷移確率からも調べる。

:::example[軌道と符号]
$0^+$ の芯に $\ell=2$ の核子を加えると、角運動量の合成で $3/2^+,5/2^+$ が候補となる。軌道間の順位は相互作用で決まる。
:::

:::example[固定配向の分布]
電荷 $Ze/2$ ずつを $z=\pm a$ に置くと各点で $3z^2-r^2=2a^2$ であり、四重極積分は $2Ze a^2$ となる。$x=\pm a$ に置く場合は $-Ze a^2$ である。量子状態の $Q_s$ ではこの積分を配向について平均する。
:::

### 静電場の多重極からみるモーメント

観測点 $\vb R$ が分布の外側にあるとき、$s=r/R<1$ と $\mu=\hat{\vb R}\cdot\hat{\vb r}$ に対し
$$
|\vb R-\vb r|^{-1}=R^{-1}(1-2s\mu+s^2)^{-1/2}
=R^{-1}\left[1+s\mu+\frac{s^2}{2}(3\mu^2-1)+\cdots\right].
$$
Taylor 展開の各項を電荷密度に掛けて積分する。核を原点に置き反転対称な定常密度を選べば、奇数の次数は角度反転で積分が消える。$\vb R$ が $z$ 軸上なら二次項が $Q_s/(2R^3)$ になる。従って四重極項は単極項より遠方で速く減る。分光学的 $Q_s$ は $\ket{I,I}$ における値なので、核が無配向の試料ではそのまま外部電場に現れない。無偏極の密度行列は多重項の恒等演算子に比例し、非ゼロ階数の成分の平均は回転対称性から消える。

:::example[配向を平均した四重極]
固定配向で正の四重極をもつ二点電荷を考える。軸を等方的に平均すれば、成分 $3z^2-r^2$ の球面平均は $r^2(3/3-1)=0$ である。二乗半径はどの配向でも $r^2$ なので、同じ平均で消えない。
:::

### 磁気双極子と四重極に対する測定の条件

静磁場だけに置いた核では相互作用 $-\hat{\vb\mu}\cdot\vb B$ を多重項に射影する。同じ $I$ の既約表現内の階数一の演算子は一種類なので $\hat{\vb\mu}=g_I\mu_N\vb I/\hbar$ となる。最大射影の期待値から $g_I=\mu/(I\mu_N)$ である。したがって隣接準位の差の大きさは $|g_I|\mu_N B$、その差に対応する振動数はこれを $h=2\pi\hbar$ で割ったものとなる。ゼロ磁場で離れた別の準位との混合は、この一次近似の適用範囲外である。

一様な場の勾配だけを残すと電気エネルギーは $\frac12\sum_{ij}\partial_i\partial_j\phi\int r_i r_j\rho\dd^3r$ である。核の位置に外部電荷がなければ $\nabla^2\phi=0$ なので積分の等方成分は消える。ゆえに四重極結合はモーメントと場の勾配の両方を要する。

:::example[隣接射影間の振動数]
一次 Zeeman 分裂を使う。単位は $b_N=\mu_N/h$ を MHz/T とし、核の $g_I>0$ と磁場 $B$ から振動数を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
g_{I}=2,\qquad b_{N}=\frac{38113}{5000}\,\mathrm{MHz/T},\qquad B=1\,\mathrm{T}
$$

この条件で求める量を ${\nu}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\nu}&=B \cdot b_{N} \cdot g_{I}\\&=\left(1\right) \cdot \left(2\right) \cdot \left(\frac{38113}{5000}\right)\\&\simeq 15.2452\,\mathrm{MHz}.\end{aligned}
$$

これは同一多重項での隣接準位の間隔である。磁気モーメントの符号を決めるには射影とエネルギーの順序も測る。

同じ核を強い磁場に入れた場合と比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
B=3\,\mathrm{T}
$$

このときの値を ${\nu}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\nu}_{\mathrm{new}}&=\left(2\right) \cdot \left(3\right) \cdot \left(\frac{38113}{5000}\right)\\&\simeq 45.7356\,\mathrm{MHz},\\{\nu}_{\mathrm{new}}-{\nu}&\simeq 30.4904\,\mathrm{MHz}.\end{aligned}
$$
:::

:::example[軸対称楕円体の四重極]
半軸 $a,a,c$ の一様な電荷分布をとる。一様単位球では各座標の二乗平均が $1/5$ なので、座標の伸縮後は $x^2,y^2,z^2$ の平均が $a^2/5,a^2/5,c^2/5$ になる。

量の定義に従い、共通の単位で次の値を用いる。
$$
Z=20,\qquad a=3\,\mathrm{fm},\qquad c=4\,\mathrm{fm}
$$

この条件で求める量を ${Q_0/e}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{Q_0/e}&=\frac{2 \cdot Z \cdot \left(- a^{2} + c^{2}\right)}{5}\\&=\frac{2 \cdot \left(20\right) \cdot \left(- \left(3\right)^{2} + \left(4\right)^{2}\right)}{5}\\&=56\,\mathrm{fm}^2.\end{aligned}
$$

符号は固定した対称軸に沿う伸び縮みを区別する。実験室での分光学的値には回転状態による配向平均が加わる。

同じ横半径で扁平な分布に変える。そこで次の値に変更し、ほかの入力を固定する。
$$
c=\frac{5}{2}\,\mathrm{fm}
$$

このときの値を ${Q_0/e}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{Q_0/e}_{\mathrm{new}}&=\frac{2 \cdot \left(20\right) \cdot \left(- \left(3\right)^{2} + \left(\frac{5}{2}\right)^{2}\right)}{5}\\&=-22\,\mathrm{fm}^2,\\{Q_0/e}_{\mathrm{new}}-{Q_0/e}&=-78\,\mathrm{fm}^2.\end{aligned}
$$
:::

:::example[整列と無偏極の違い]
スピン 1 のアンサンブルで軸方向の $M=0$ の占有割合を $w_0$ とし、$M=\pm1$ を等しく準備する。$\langle 3M^2-I(I+1)\rangle$ を計算する。

量の定義に従い、共通の単位で次の値を用いる。
$$
w_{0}=\frac{1}{3}
$$

この条件で求める量を ${\mathcal Q}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\mathcal Q}&=1 - 3 \cdot w_{0}\\&=1 - 3 \cdot \left(\frac{1}{3}\right)\\&=0.\end{aligned}
$$

磁化に比例する $\langle M\rangle$ は両条件でゼロである。二次の整列は磁化を測るだけでは分からない。

中央の射影の占有を除いた場合を比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
w_{0}=0
$$

このときの値を ${\mathcal Q}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\mathcal Q}_{\mathrm{new}}&=1 - 3 \cdot \left(0\right)\\&=1,\\{\mathcal Q}_{\mathrm{new}}-{\mathcal Q}&=1.\end{aligned}
$$
:::

:::example[単極子と四重極の距離則]
軸上で二次のポテンシャル項を単極項で割る。$Q_0/e$ を $q_2$ とし、観測点の距離を $d$ とする。

量の定義に従い、共通の単位で次の値を用いる。
$$
q_{2}=40\,\mathrm{fm}^2,\qquad Z=20,\qquad d=10\,\mathrm{fm}
$$

この条件で求める量を ${u}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{u}&=\frac{q_{2}}{2 \cdot Z \cdot d^{2}}\\&=\frac{\left(40\right)}{2 \cdot \left(10\right)^{2} \cdot \left(20\right)}\\&\simeq 0.01.\end{aligned}
$$

ここでは固定配向で分布の外側を観測している。無偏極の試料では軸方向の二次項をこのまま平均できない。

遠方に移動して相対寄与を調べる。そこで次の値に変更し、ほかの入力を固定する。
$$
d=20\,\mathrm{fm}
$$

このときの値を ${u}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{u}_{\mathrm{new}}&=\frac{\left(40\right)}{2 \cdot \left(20\right)^{3}}\\&\simeq 0.0025,\\{u}_{\mathrm{new}}-{u}&\simeq -0.0075.\end{aligned}
$$
:::

## アイソスピン

:::definition[アイソスピン]{#def-nuc-isospin}
:term[アイソスピン]{en="isospin" yomi="あいそすぴん"} は陽子と中性子を近似対称性に基づく $\SU(2)$ 二重項として扱う量子数である。本科目では核子の $t=1/2$ に対して $t_z(n)=+1/2$、$t_z(p)=-1/2$ とし、核全体では $T_z=(N-Z)/2$ とする。
:::

電磁相互作用や核子の質量差を入れると対称性は破れる。二重項を二つ合成する場合、対称な $\ket{nn},(\ket{np}+\ket{pn})/\sqrt2,\ket{pp}$ が $T=1$ をなし、反対称な $(\ket{np}-\ket{pn})/\sqrt2$ が $T=0$ をなす。後者に上げ下げ演算子を作用させると二項が打ち消し合う。相対軌道 $L$、合成スピン $S$、アイソスピン $T$ の各交換固有値は $(-1)^L,(-1)^{S+1},(-1)^{T+1}$ なので、全反対称化は
$$
(-1)^{L+S+T}=-1
\label{eq-nuc-pauli-isospin}
$$
を要求する。

:::example[二核子の S 波]
$L=0$ では $S+T$ が奇数である。従って $(S,T)=(1,0),(0,1)$ が許される。陽子二つでは $T_z=-1$ なので $T=1$ となり、低エネルギーの S 波は $S=0$ である。
:::

:::example[鏡映核]
${}^{3}\mathrm H$ と ${}^{3}\mathrm{He}$ の基底状態は近似的に $T=1/2$ の多重項をなし、第三成分はそれぞれ $+1/2,-1/2$ である。この規約では $T_-$ が $\beta^-$ の $n\to p$ に対応する。
:::

### アイソスピン射影と類似状態

対称性が成り立つ極限では、一つの $T$ の既約表現の多重項に $2T+1$ 個の射影がある。$T_+$ と $T_-$ のノルムは $\vb T^2=T_z^2+(T_+T_-+T_-T_+)/2$ および交換子 $[T_+,T_-]=2T_z$ で決まる。この二式を連立すると $T_+T_-=\vb T^2-T_z(T_z-1)$ である。最大射影から順に下げると同じ $T$ の類似状態を生成し、最小射影ではこれ以上下げられない。

対称性の破れを一次摂動で考えると、アイソスカラー、ベクトルの第三成分、階数二の第三成分の期待値はそれぞれ定数、$T_z$、$3T_z^2-T(T+1)$ に比例する。最後の形は同じ表現内で階数二の成分を $T_iT_j$ の対称トレースレス部分で表せば得られる。この近似では同じ多重項の質量曲線は二次式となる。状態混合による二次摂動やさらに高階の寄与は、この二次式からのずれを生む。

:::example[Fermi 遷移の強度の比]
核の規約で $T_-$ は $T_z$ を一つ下げる。強い相互作用の対称性を仮定し、同じ $T$ 内での状態のノルム二乗を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
T=2,\qquad m=2
$$

この条件で求める量を ${B_F}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{B_F}&=T \cdot \left(T + 1\right) - m \cdot \left(m - 1\right)\\&=- \left(2\right) \cdot \left(\left(2\right) - 1\right) + \left(2\right) \cdot \left(\left(2\right) + 1\right)\\&=4.\end{aligned}
$$

強度は下げ演算子の振幅ではなく、その二乗である。放出レプトンの位相空間はまだ含めていない。

同じ多重項の別の親準位を比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
m=0
$$

このときの値を ${B_F}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{B_F}_{\mathrm{new}}&=- \left(0\right) \cdot \left(\left(0\right) - 1\right) + \left(2\right) \cdot \left(\left(2\right) + 1\right)\\&=6,\\{B_F}_{\mathrm{new}}-{B_F}&=2.\end{aligned}
$$
:::

:::example[多重項の二次係数]
$T=1$ の類似状態三つの質量超過を $D_-,D_0,D_+$ とし、$D(m)=a+bm+cm^2$ に当てはめる。端の平均から中央を引くと $a,b$ が消える。

量の定義に従い、共通の単位で次の値を用いる。
$$
D_{p}=2\,\mathrm{MeV},\qquad D_{m}=4\,\mathrm{MeV},\qquad D_{0}=\frac{27}{10}\,\mathrm{MeV}
$$

この条件で求める量を ${c}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{c}&=- D_{0} + \frac{D_{m} + D_{p}}{2}\\&=\frac{\left(2\right)}{2} + \frac{\left(4\right)}{2} - \left(\frac{27}{10}\right)\\&\simeq 0.3\,\mathrm{MeV}.\end{aligned}
$$

三点だけでは必ず二次式が存在するため、その精度を実験的に検証するにはさらに射影数の多い多重項が必要である。

中央の読み取り誤差が二次係数にどう伝わるかを見る。そこで次の値に変更し、ほかの入力を固定する。
$$
D_{0}=\frac{14}{5}\,\mathrm{MeV}
$$

このときの値を ${c}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{c}_{\mathrm{new}}&=\frac{\left(2\right)}{2} + \frac{\left(4\right)}{2} - \left(\frac{14}{5}\right)\\&\simeq 0.2\,\mathrm{MeV},\\{c}_{\mathrm{new}}-{c}&\simeq -0.1\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[二準位のアイソスピン不純度]
同じ $T_z$ で異なる $T$ の二準位が $d>0$ だけ離れているとする。非対角摂動 $v$ で下側の状態に入る上側成分は一次の振幅が $-v/d$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
v=\frac{1}{20}\,\mathrm{MeV},\qquad d=1\,\mathrm{MeV}
$$

この条件で求める量を ${P}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{P}&=\left(\frac{v}{d}\right)^{2}\\&=\frac{\left(\frac{1}{20}\right)^{2}}{\left(1\right)^{2}}\\&\simeq 0.0025.\end{aligned}
$$

表示した確率は $v/d$ の二次までである。準位間隔が摂動程度になる場合は二準位を厳密に対角化する。

準位が混ざりやすい小さな間隔と比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
d=\frac{1}{2}\,\mathrm{MeV}
$$

このときの値を ${P}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{P}_{\mathrm{new}}&=\frac{\left(\frac{1}{20}\right)^{2}}{\left(\frac{1}{2}\right)^{2}}\\&\simeq 0.01,\\{P}_{\mathrm{new}}-{P}&\simeq 0.0075.\end{aligned}
$$
:::

:::example[同じ A の中性子過剰度]
規約 $T_z=(N-Z)/2$ から核の射影を求める。偶奇によって許される全 $T$ は整数または半整数で、$T\ge \abs{T_z}$ を満たす。

量の定義に従い、共通の単位で次の値を用いる。
$$
A=48,\qquad Z=20
$$

この条件で求める量を ${T_z}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{T_z}&=\frac{A}{2} - Z\\&=- \left(20\right) + \frac{\left(48\right)}{2}\\&=4.\end{aligned}
$$

$\abs{T_z}$ は全 $T$ の下限である。特定の励起状態の $T$ を決めるときは遷移強度や類似状態の位置も使う。

同じ質量数で陽子数を二つ増やす。そこで次の値に変更し、ほかの入力を固定する。
$$
Z=22
$$

このときの値を ${T_z}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{T_z}_{\mathrm{new}}&=- \left(22\right) + \frac{\left(48\right)}{2}\\&=2,\\{T_z}_{\mathrm{new}}-{T_z}&=-2.\end{aligned}
$$
:::

## 演習問題

::::exercise{level=basic}
${}^{40}_{20}\mathrm{Ca},{}^{42}_{20}\mathrm{Ca},{}^{42}_{22}\mathrm{Ti}$ の $(N,Z)$ と相互の関係を求めよ。
:::solution
$N=A-Z$ より $(20,20),(22,20),(20,22)$ である。最初と二番目は同位体、最初と三番目は同中性子体、二番目と三番目は同重体である。
:::
::::

::::exercise{level=basic}
$R=5\,\mathrm{fm}$ の一様球の二乗半径と $F$ の $q^2$ に対する初期勾配を求めよ。
:::solution
一様球の平均は $3R^2/5=15\,\mathrm{fm^2}$ である。:eqref[eq-nuc-charge-radius] より勾配は $-15/6=-2.5\,\mathrm{fm^2}$ である。
:::
::::

::::exercise{level=basic}
$B(A,Z)=112$、$B(A-1,Z)=106$、$B(A-1,Z-1)=104\,\mathrm{MeV}$ のとき分離エネルギーと $7\,\mathrm{MeV}$ で開く一核子チャネルを求めよ。
:::solution
:eqref[eq-nuc-separation] から $S_n=6$、$S_p=8\,\mathrm{MeV}$ である。$6<7<8$ なので中性子放出だけが開き、放出後の相対運動エネルギーは $1\,\mathrm{MeV}$ である。
:::
::::

::::exercise{level=basic}
$I=3/2$、$\mu=3\mu_N$ の核で $E_{1/2}-E_{3/2}$ を求めよ。
:::solution
$E_M=-(M/I)\mu B_0$ なので、差は $(1/I)\mu B_0=2\mu_N B_0$ である。
:::
::::

::::exercise{level=advanced}
半径 $R$ の薄い球殻に一様に載った電荷の形状因子と初期勾配を求め、一様球と比較せよ。
:::solution
すべての電荷で半径が $R$ なので方位の平均 $F=\tfrac12\int_{-1}^1e^{iqR\mu}\dd\mu=\sin(qR)/(qR)$ を取る。$q^2$ の係数は $-R^2/6$ となる。一様球では $-R^2/10$ なので、勾配の大きさは球殻の方が $5/3$ 倍大きい。
:::
::::

::::exercise{level=advanced}
質量エネルギー $M_*c^2+C(Z-Z_*)^2$ で、$Z\pm1$ のどちらにも質量差が正にならない整数 $Z$ の範囲を導け。
:::solution
$Z+1$ への差は $C[-2(Z-Z_*)-1]$、$Z-1$ への差は $C[2(Z-Z_*)-1]$ である。$C>0$ なので両者が非正である条件は $|Z-Z_*|\le1/2$ である。半整数の頂点では二つの $Z$ が同じ質量をもつ。
:::
::::

::::exercise{level=advanced}
二核子の $P$ 波で許される $(S,T)$ と $J^\pi$ を求めよ。
:::solution
$L=1$ を :eqref[eq-nuc-pauli-isospin] に入れると $S+T$ は偶数である。よって $(0,0)$ なら $J^\pi=1^-$、$(1,1)$ なら $J^\pi=0^-,1^-,2^-$ となる。$J$ は $|L-S|$ から $L+S$ まで、パリティは $(-1)^L$ である。
:::
::::

::::exercise{level=advanced}
$T=2,T_z=2$ の状態に無次元の $T_-^2$ を作用させよ。
:::solution
角運動量の代数により $T_-\ket{T,m}=\sqrt{T(T+1)-m(m-1)}\ket{T,m-1}$ である。$m=2$ からは係数 2、$m=1$ からは $\sqrt6$ となるため、答えは $2\sqrt6\ket{2,0}$ である。
:::
::::

## 確認テスト

::::quiz
$\beta^-$ 崩壊による $(\Delta N,\Delta Z)$ を選べ。

- [ ] $(1,-1)$
- [x] $(-1,1)$
- [ ] $(-2,-2)$

:::explanation
中性子一つが陽子に変わる。最初の組は逆向きの $\beta^+$ または電子捕獲であり、最後は $\alpha$ 放出である。
:::
::::

::::quiz{answer="2.5" tol="1e-6"}
$B=30\,\mathrm{MeV}$、$A=12$ のとき $B/A$ を MeV 単位で求めよ。
:::explanation
核子で等分すると $30/12=2.5\,\mathrm{MeV}$ である。
:::
::::

::::quiz{answer="3" tol="1e-6"}
$q$ の単位を $\mathrm{fm}^{-1}$ として $F=1-1.5q^2+O(q^4)$ である。電荷の二乗平均平方根半径を fm 単位で求めよ。
:::explanation
:eqref[eq-nuc-charge-radius] より $\langle r^2\rangle=6\cdot1.5=9\,\mathrm{fm^2}$ なので $3\,\mathrm{fm}$ である。
:::
::::

::::quiz
$I=1/2$ の静的モーメントについて成り立つものを選べ。

- [x] 磁気双極子は許され、電気四重極子は消える。
- [ ] 四重極子だけが許される。
- [ ] 四重極子が消えれば固有座標での密度も球対称である。

:::explanation
対角行列要素の階数は $2I=1$ 以下である。磁気双極子は階数 1、電気四重極子は階数 2 である。この条件は実験室座標での対角成分を制限する。
:::
::::

::::quiz{answer="-1" tol="1e-6"}
$N=8,Z=10$ の核の $T_z$ を求めよ。
:::explanation
本章の規約で $T_z=(8-10)/2=-1$ である。
:::
::::

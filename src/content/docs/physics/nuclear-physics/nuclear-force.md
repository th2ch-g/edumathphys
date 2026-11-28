---
title: 核力と二核子系
description: 重陽子と低エネルギー散乱から核力のスピン、アイソスピン、距離への依存性を調べる。
notation:
  - symbol: 'S_{12}'
    meaning: 二核子テンソル演算子
  - symbol: 'a,r_e'
    meaning: 散乱長と有効距離
  - symbol: '\kappa'
    meaning: 束縛状態の外部減衰波数
---

## この章の目標

- 重陽子の波動関数を相対座標で表す。
- 位相差のしきい値展開を散乱長と束縛極に結び付ける。
- 核力をスピンとアイソスピンのチャネルに分ける。
- 一中間子交換と短距離・三体相互作用の役割を説明する。

## 重陽子の束縛状態とテンソル力

陽子と中性子を自由にした状態を相対エネルギーの原点とする。質量を $m_p,m_n$、換算質量を $\mu=m_pm_n/(m_p+m_n)$ とすると、相対 Hamiltonian は $-\hbar^2\nabla^2/(2\mu)+V$ である。重陽子の束縛エネルギーは約 $2.2246\,\mathrm{MeV}$、量子数は $J^\pi=1^+$ である。$T=0$ の記述では、:eqref[eq-nuc-pauli-isospin] に偶数 $L$ を入れると $S=1$ となる。角運動量の三角条件から、このスピンと $J=1$ に結合できる偶数 $L$ は $0,2$ である。

:::definition[テンソル演算子]{#def-nuc-tensor-force}
二粒子の Pauli ベクトルを $\vb\sigma_1,\vb\sigma_2$ とすると、:term[テンソル力]{en="tensor force" yomi="てんそるりょく"} の角度依存性は
$$
S_{12}=3(\vb\sigma_1\cdot\hat{\vb r})(\vb\sigma_2\cdot\hat{\vb r})-\vb\sigma_1\cdot\vb\sigma_2
$$
で表せる。空間とスピンの同時回転に対してスカラーであり、別々の回転に対しては一般に不変でない。
:::

$S$ 波と $D$ 波のスピン角関数を $\mathcal Y_S,\mathcal Y_D$ とし、動径部分を $u(r)/r,w(r)/r$ とする。$S_{12}\mathcal Y_S$ は角度平均ゼロで、$J=1$ の許された部分空間のうち $D$ 波だけを含む。三重項で $\hat{\vb r}$ を量子化軸とすると $S_{12}$ の固有値は $2,-4,2$ なので $S_{12}^2=8-2S_{12}$ である。従って $\|S_{12}\mathcal Y_S\|^2=8$、位相を $\mathcal Y_D=S_{12}\mathcal Y_S/\sqrt8$ と取れば $S_{12}\mathcal Y_D=\sqrt8\mathcal Y_S-2\mathcal Y_D$ となる。中心力 $V_C$ とテンソル項 $V_T S_{12}$ だけの場合、
$$
\left[-\frac{\hbar^2}{2\mu}\frac{\dd^2}{\dd r^2}+
\begin{pmatrix}
V_C & \sqrt8 V_T\\
\sqrt8 V_T & V_C-2V_T+6\hbar^2/(2\mu r^2)
\end{pmatrix}\right]\binom{u}{w}
=-B_d\binom{u}{w}.
\label{eq-nuc-deuteron-coupled}
$$
原点で正則、無限遠で減衰する解を取り、$\int_0^\infty(u^2+w^2)\dd r=1$ とする。核力の範囲外では $u\propto e^{-\kappa r}$、$w\propto e^{-\kappa r}[1+3/(\kappa r)+3/(\kappa r)^2]$ であり、$\kappa=\sqrt{2\mu B_d}/\hbar$ である。後者は外部方程式 $w''=(\kappa^2+6/r^2)w$ への代入で確認できる。

:::example[束縛波の長さ]
$\mu c^2=469.5\,\mathrm{MeV}$ とすると $\kappa\simeq0.232\,\mathrm{fm}^{-1}$ で、減衰長は約 $4.3\,\mathrm{fm}$ である。後述の pion の到達距離より長いため、波動関数の外部部分が重要になる。
:::

:::example[純粋な S 波]
陽子を点電荷とすれば、純粋な $L=0$ の状態の空間密度は球対称なので四重極積分は消える。実測の重陽子は非ゼロの四重極モーメントを持つ。二体模型では、$D$ 波との干渉が非ゼロの主要因となる。$D$ 波確率自体は力の模型や短距離の表現に依存する。
:::

### 束縛のしきい値を有限井戸で見る

テンソル力をいったん外し、$r<R$ で $V=-V_0$、外側で $V=0$ の S 波を解く。束縛エネルギー $B$ に対し内側の規則解は $u=A\sin Kr$、$K=\sqrt{2\mu(V_0-B)}/\hbar$ である。外側は $u=A\sin(KR)e^{-\kappa(r-R)}$、$\kappa=\sqrt{2\mu B}/\hbar$ である。$u,u'$ を境界で一致させると $K\cot(KR)=-\kappa$ になる。最初の束縛状態が現れる $B\downarrow0$ の極限では $KR\to\pi/2$ なので、$V_{0,\rm th}=\pi^2\hbar^2/(8\mu R^2)$ が必要である。

最も浅い束縛状態の周辺で $\sqrt{2\mu V_0}R/\hbar=\pi/2+\varepsilon$、$0<\varepsilon\ll1$ とする。$\kappa R$ は小さく、$KR=\sqrt{(\pi/2+\varepsilon)^2-(\kappa R)^2}$ である。境界条件の展開により $\kappa R=(\pi/2)\varepsilon+O(\varepsilon^2)$ となる。一方、演習の散乱長は $\tan(\pi/2+\varepsilon)\simeq-1/\varepsilon$ により $a/R\simeq2/(\pi\varepsilon)$ である。従って $\kappa a\simeq1$ を得る。束縛エネルギーは井戸の深さの余分な分の一次ではなく、しきい値近くで二次から増える。外の尾が長いときに短距離の形への感度が弱まることを示す例である。

:::example[浅い状態の半径]
井戸の外の尾が規格化を支配するなら、$u\simeq\sqrt{2\kappa}e^{-\kappa r}$ として $\langle r^2\rangle\simeq1/(2\kappa^2)$ である。導出は $\int_0^\infty r^2e^{-2\kappa r}\dd r=2/(2\kappa)^3$ による。ここで $r$ は二核子の相対距離である。両質量が等しく点核子なら、一核子の中心から測る二乗距離はその $1/4$ になる。
:::

### 内部座標と漸近規格化

$\vb r=\vb r_1-\vb r_2$、$\vb R=(m_1\vb r_1+m_2\vb r_2)/(m_1+m_2)$ と置く。鎖則を使うと二つの運動エネルギーの交差微分係数は $m_1/[m_1(m_1+m_2)]-m_2/[m_2(m_1+m_2)]=0$ となり、相対運動の係数には $1/\mu=1/m_1+1/m_2$ が現れる。束縛エネルギー $B$ はこの相対運動に対して定義する。

ポテンシャルの範囲外で $x=\kappa r$ と置くと $D$ 波の方程式は $w_{xx}=(1+6/x^2)w$ である。$w=e^{-x}f(x)$ を入れると $f''-2f'=6f/x^2$ となり、$f=1+3/x+3/x^2$ が直接代入で解と分かる。従って $w/u$ の漸近値は尾の係数 $A_D/A_S$ の比となり、全空間の確率 $P_D$ とは異なる。

:::example[束縛と減衰長]
束縛エネルギー $B$ と換算質量の静止エネルギー $m_c$ から外部 $S$ 波の減衰長を得る。

量の定義に従い、共通の単位で次の値を用いる。
$$
h_{c}=\frac{197327}{1000}\,\mathrm{MeV\,fm},\qquad m_{c}=469\,\mathrm{MeV},\qquad B=\frac{89}{40}\,\mathrm{MeV}
$$

この条件で求める量を ${L}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{L}&=\frac{h_{c}}{\sqrt{2 \cdot B \cdot m_{c}}}\\&=\frac{\sqrt{2} \cdot \left(\frac{197327}{1000}\right)}{2 \cdot \sqrt{\left(469\right) \cdot \left(\frac{89}{40}\right)}}\\&\simeq 4.31937\,\mathrm{fm}.\end{aligned}
$$

浅い束縛ほど外部に広がる。内部でこの指数をそのまま使う必要はない。

束縛を弱めた仮想状態と比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
B=1\,\mathrm{MeV}
$$

このときの値を ${L}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{L}_{\mathrm{new}}&=\frac{\sqrt{2} \cdot \left(\frac{197327}{1000}\right)}{2 \cdot \sqrt{\left(1\right) \cdot \left(469\right)}}\\&\simeq 6.44295\,\mathrm{fm},\\{L}_{\mathrm{new}}-{L}&\simeq 2.12359\,\mathrm{fm}.\end{aligned}
$$
:::

:::example[零距離模型の二乗平均平方根半径]
規格化波 $u=\sqrt{2\kappa}e^{-\kappa r}$ を全域で仮定する。$\int r^2u^2\dd r=1/(2\kappa^2)$ であり、等質量なら一核子の中心からの距離は $r/2$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
\kappa=\frac{23}{100}\,\mathrm{fm}^{-1}
$$

この条件で求める量を ${r_{\rm rms}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{r_{\rm rms}}&=1 \cdot \frac{1}{\sqrt{8} \cdot \kappa}\\&=\frac{\sqrt{2}}{4 \cdot \left(\frac{23}{100}\right)}\\&\simeq 1.53719\,\mathrm{fm}.\end{aligned}
$$

点核子の半径である。電荷半径と比較するには核子の電荷分布も含める。

より小さな減衰定数で評価する。そこで次の値に変更し、ほかの入力を固定する。
$$
\kappa=\frac{3}{20}\,\mathrm{fm}^{-1}
$$

このときの値を ${r_{\rm rms}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{r_{\rm rms}}_{\mathrm{new}}&=\frac{\sqrt{2}}{4 \cdot \left(\frac{3}{20}\right)}\\&\simeq 2.35702\,\mathrm{fm},\\{r_{\rm rms}}_{\mathrm{new}}-{r_{\rm rms}}&\simeq 0.819834\,\mathrm{fm}.\end{aligned}
$$
:::

:::example[尾の D 波と S 波の振幅比]
漸近係数比を $\eta=A_D/A_S$ とする。無次元座標 $x$ における局所の比 $w/u$ を調べる。

量の定義に従い、共通の単位で次の値を用いる。
$$
\eta=\frac{1}{40},\qquad x=2
$$

この条件で求める量を ${w/u}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{w/u}&=\eta \cdot \left(1 + \frac{3}{x^{2}} + \frac{3}{x}\right)\\&=\left(\frac{1}{40}\right) \cdot \left(1 + \frac{3}{\left(2\right)} + \frac{3}{\left(2\right)^{2}}\right)\\&\simeq 0.08125.\end{aligned}
$$

この比は $x\to\infty$ で $\eta$ に近づく。小さな $x$ へ漸近形を外挿して内部確率を求めることはできない。

尾のより遠い測定点と比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
x=5
$$

このときの値を ${w/u}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{w/u}_{\mathrm{new}}&=\left(\frac{1}{40}\right) \cdot \left(1 + \frac{3}{\left(5\right)} + \frac{3}{\left(5\right)^{2}}\right)\\&\simeq 0.043,\\{w/u}_{\mathrm{new}}-{w/u}&\simeq -0.03825.\end{aligned}
$$
:::

:::example[混合による固有値の低下]
孤立した $S,D$ 二状態を未混合で $0,d$ とし、非対角要素 $v$ を入れる。固有方程式 $E(E-d)=v^2$ の下側の根を使う。

量の定義に従い、共通の単位で次の値を用いる。
$$
d=4\,\mathrm{MeV},\qquad v=\frac{2}{5}\,\mathrm{MeV}
$$

この条件で求める量を ${E_-}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{E_-}&=\frac{d - \sqrt{d^{2} + 4 \cdot v^{2}}}{2}\\&=\frac{\left(4\right)}{2} - \frac{\sqrt{\left(4\right)^{2} + 4 \cdot \left(\frac{2}{5}\right)^{2}}}{2}\\&\simeq -0.0396078\,\mathrm{MeV}.\end{aligned}
$$

混合の符号を変えても固有値は変わらないが、固有ベクトルの相対位相は変わる。

テンソル混合が大きい場合を比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
v=\frac{4}{5}\,\mathrm{MeV}
$$

このときの値を ${E_-}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{E_-}_{\mathrm{new}}&=\frac{\left(4\right)}{2} - \frac{\sqrt{\left(4\right)^{2} + 4 \cdot \left(\frac{4}{5}\right)^{2}}}{2}\\&\simeq -0.154066\,\mathrm{MeV},\\{E_-}_{\mathrm{new}}-{E_-}&\simeq -0.114458\,\mathrm{MeV}.\end{aligned}
$$
:::

### 四重極モーメントの S 波と D 波の干渉項

陽子を点電荷とし、陽子と中性子を等質量とする。電荷を $e$ で割ると、相対座標での四重極演算子は $r^2P_2(\cos\theta)/2$ になる。これは陽子の重心からの座標が $\vb r/2$ であり、$3z^2-r^2=2r^2 P_2$ だからである。$M_J=1$ の $D$ 波のスピン射影の順を $1,0,-1$ に取ると、そのスピン角関数の成分は
$$
\left(\frac{Y_{20}}{\sqrt{10}},-\sqrt{\frac3{10}}Y_{21},\sqrt{\frac35}Y_{22}\right)^{\mathsf T}
$$
である。この位相は :eqref[eq-nuc-deuteron-coupled] で $S$–$D$ 成分を $+\sqrt{8}V_T$ とするものと一致する。実際、$S$ 波に $S_{12}$ を作用させると、この三成分に $\sqrt8$ を掛けたものが得られる。$Y_{00}$ の成分が $M_S=1$ のみなので、$S$–$D$ の角度積分は $(1/\sqrt{10})\int Y_{00}^*P_2Y_{20}\dd\Omega=1/\sqrt{50}$ になる。一方、$D$–$D$ の三つの対角積分 $\int |Y_{2m}|^2P_2\dd\Omega$ は $m=0,1,2$ について $2/7,1/7,-2/7$ である。これは $Y_{2m}$ の絶対値二乗が極角の余弦 $x$ の多項式 $(3x^2-1)^2,x^2(1-x^2),(1-x^2)^2$ に比例することから、各多項式に $P_2$ を掛けて $[-1,1]$ で積分し、元の多項式の積分で割れば求まる。スピンごとの重みで足すと $-1/10$ になる。従って
$$
\frac{Q_d}{e}=\frac1{20}\int_0^\infty r^2w(r)\,[\sqrt8u(r)-w(r)]\,\dd r.\label{eq-nuc-deuteron-quadrupole}
$$
:eqref[eq-nuc-deuteron-quadrupole] で小さな $D$ 波では、$w^2$ の対角項よりも $uw$ の干渉項が先に効く。したがって四重極モーメントを $P_D$ の制約に使う場合、動径重なりと有効な電荷演算子をモデルに含める。

## 低エネルギー核子-核子散乱と有効距離理論

電気力を無視できる短距離ポテンシャルを考える。位相差の定義と散乱振幅の一般論は [散乱理論](/physics/advanced-quantum-mechanics/scattering/) を参照する。S 波の外部動径解は $\sin(kr+\delta_0)$ であり、$k\to0$ の解を $u_0\propto r-a$ と規格化する。$a$ を :term[散乱長]{en="scattering length" yomi="さんらんちょう"} という。$\tan\delta_0/k\to-a$ から $k\cot\delta_0\to-1/a$ である。

有限範囲の有界なポテンシャルでは、原点からの解の対数微分 $L(k)=u_k'(R)/u_k(R)$ が、$u_k(R)\ne0$ において $k^2$ の解析関数となる。例えば $V=0$ となる範囲より外に $R$ を選び、積分方程式 $u_k(r)=r+\int_0^r(r-s)[2\mu V(s)/\hbar^2-k^2]u_k(s)\dd s$ を反復する。$k^2$ の任意のコンパクト集合で括弧を $M$ 以下に抑えると、第 $n$ 回の補正は $M^n r^{2n+1}/(2n+1)!$ 以下である。微分した級数も収束し、分母が非零の範囲で対数微分は解析的になる。外部解の $L(k)=k\cot(kR+\delta_0)$ に余接の加法公式を使うと
$$
k\cot\delta_0=\frac{k^2+L(k)k\cot(kR)}{k\cot(kR)-L(k)}.
$$
この比が $k=0$ で正則な範囲で、定数項と二次項を分離する。
$$
k\cot\delta_0=-\frac1a+\frac12r_e k^2+O(k^4),\qquad
f_0=\frac{e^{2i\delta_0}-1}{2ik}=\frac1{k\cot\delta_0-ik}.
\label{eq-nuc-effective-range}
$$
$r_e$ が :term[有効距離]{en="effective range" yomi="ゆうこうきょり"} である。区別可能なスカラー粒子では等方的な S 波断面積は $4\pi|f_0|^2$ となる。束縛状態は解析接続した $k=i\kappa$ における極に対応し、二次まででは $\kappa=1/a+r_e\kappa^2/2$ である。

:::example[符号と断面積]
$a=+5\,\mathrm{fm}$ と $a=-5\,\mathrm{fm}$ の場合、ともに $k\to0$ の断面積は $100\pi\,\mathrm{fm^2}=\pi\,\mathrm b$ である。断面積だけでは符号が決まらない。有効距離を無視した極は前者が $+i/5\,\mathrm{fm}^{-1}$、後者が $-i/5\,\mathrm{fm}^{-1}$ にあり、後者の外部波は正規化できない。
:::

:::example[剛体球]
$r=R$ で波がゼロなら $\delta_0=-kR$ と取れる。$k\cot(-kR)=-1/R+Rk^2/3+\cdots$ なので $a=R$、$r_e=2R/3$ である。正の散乱長だけで束縛状態を断定できないことが分かる。
:::

### ゼロエネルギー波から有効距離を計算する

角運動量ゼロの有限範囲ポテンシャルを考える。$a\ne0$ とし、外側解を $u_k(r)=\cos kr+(k\cot\delta_0)\sin(kr)/k$ と規格化する。$k=0$ の解 $u_0$ は外側で $1-r/a$ になる。どちらも原点ではゼロであり、動径方程式の差から
$$
(u_k'u_0-u_0'u_k)'=-k^2u_ku_0
$$
を得る。任意の外側の位置 $R$ まで積分し、$k\to0$ とすれば右辺を $-k^2\int_0^R u_0^2\dd r$ で近似できる。左辺では :eqref[eq-nuc-effective-range] により、外側の $k^2$ の係数が $g(r)=-r^2/2+r^3/(6a)+(r_e/2)r$ となる。従って表面値を計算すると
$$
g'(R)u_0(R)-u_0'(R)g(R)=
\frac{r_e}{2}-R+\frac{R^2}{a}-\frac{R^3}{3a^2}.
$$
これを積分と等置し、
$$
r_e=2\int_0^\infty\left[(1-r/a)^2-u_0(r)^2\right]\dd r
$$
を得る。ポテンシャルの外では被積分関数がゼロである。剛体球なら外部境界の内側を $u_0=0$ として延長でき、結果は $2\int_0^a(1-r/a)^2\dd r=2a/3$ である。

### 三重項の結合チャネル

テンソル力のある $J=1$ の三重項には二つの入射部分波がある。弾性領域では出射振幅を $S$ 波と $D$ 波に分けると、入射振幅からの写像は二行二列のユニタリ行列 $\mathsf S$ になる。時間反転対称性に適合する位相を選ぶと $\mathsf S$ は複素対称行列である。その実部と虚部は実対称で、ユニタリ性から互いに可換である。従って二つを同じ直交変換 $O(\epsilon)$ で対角化でき、
$$
\mathsf S=O(\epsilon)
\begin{pmatrix}e^{2i\delta_S}&0\\0&e^{2i\delta_D}\end{pmatrix}
O(\epsilon)^{\mathsf T},\qquad
O(\epsilon)=\begin{pmatrix}\cos\epsilon&-\sin\epsilon\\\sin\epsilon&\cos\epsilon\end{pmatrix}
$$
と書ける。$\delta_S,\delta_D$ はしきい値でそれぞれ $S,D$ 波へ連続する固有位相差で、$\epsilon$ はこの規約の混合角である。非対角成分は $(e^{2i\delta_S}-e^{2i\delta_D})\sin\epsilon\cos\epsilon$ である。位相差が一致する点では混合角の選び方に任意性が残る。散乱データでテンソル力を検証するには、二つの位相差に加えて混合も調べる。以下で三重項の散乱長を用いるときは $S$ 波に連続するしきい値の成分を指す。

:::example[束縛極の近似を検査する]
三重項の代表値 $a=5.42\,\mathrm{fm},r_e=1.75\,\mathrm{fm}$ を使う。$\kappa=(1-\sqrt{1-2r_e/a})/r_e$ は約 $0.231\,\mathrm{fm}^{-1}$ であり、重陽子の長い減衰長と整合する。二次方程式のもう一つの根は約 $0.912\,\mathrm{fm}^{-1}$ であり、外挿が有効な低運動量域を外れる。打ち切った式の全ての極を物理的な束縛状態に数えると誤る。
:::

### 短距離の力の高い部分波への寄与

範囲 $R$ の中心力について $r=R$ での対数微分 $L_\ell$ を使う。外部解を $u=r[j_\ell(kr)-\tan\delta_\ell\,n_\ell(kr)]$ と書く。$j_\ell,n_\ell$ は球 Bessel 関数と球 Neumann 関数である。$L_\ell u-u'=0$ から、$\tan\delta_\ell$ は $L_\ell rj_\ell-(rj_\ell)'$ と $L_\ell rn_\ell-(rn_\ell)'$ の比として求まる。$kR\ll1$ で $j_\ell(kr)\sim(kr)^\ell/(2\ell+1)!!$、$n_\ell(kr)\sim-(2\ell-1)!!/(kr)^{\ell+1}$ である。内部の対数微分が有限で、分母の主係数がゼロにならないなら比は $k^{2\ell+1}$ から始まる。従って部分波断面積 $(4\pi/k^2)(2\ell+1)\sin^2\delta_\ell$ は $k^{4\ell}$ の桁になる。しきい値に束縛極がある場合などは、分母のゼロによりこの数え方を変える必要がある。

:::example[二つの波の比較]
散乱長がゼロでない通常の S 波と、しきい値に極がない P 波を比べる。$k$ を半分にすると、S 波断面積の主項は同じである。P 波断面積の主項は $1/16$ になる。低エネルギーで角度分布が等方に近づく理由の一つである。
:::

### 弾性ユニタリティと束縛極

弾性の $S$ 波で $S_0=e^{2i\delta_0}$ である。外向波の係数 $f_0=(S_0-1)/(2ik)$ を書き換えると $f_0^{-1}=k\cot\delta_0-ik$ となり、虚部はユニタリティで固定される。有限距離相互作用で $k\cot\delta_0$ を二次まで展開して $k=i\kappa$ を代入すると、極の方程式は $\kappa-1/a-r_e\kappa^2/2=0$ である。有効距離が自然な大きさで $|r_e/a|\ll1$ なら、小さい根は $\kappa=1/a+r_e/(2a^2)+O(r_e^2/a^3)$ となる。大きな根は一般には低エネルギー展開の範囲に入らない。

:::example[有効距離を含む断面積]
区別できる粒子の一つの $S$ 波を考え、逆振幅の実部を $-1/a+r_e k^2/2$ に近似する。断面積は $4\pi$ と振幅の絶対値二乗の積である。

量の定義に従い、共通の単位で次の値を用いる。
$$
a=5\,\mathrm{fm},\qquad r_{e}=\frac{17}{10}\,\mathrm{fm},\qquad k=\frac{1}{20}\,\mathrm{fm}^{-1}
$$

この条件で求める量を ${\sigma}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\sigma}&=\frac{4 \cdot \pi}{k^{2} + \left(\frac{k^{2} \cdot r_{e}}{2} - \frac{1}{a}\right)^{2}}\\&=\frac{4 \cdot \pi}{\left(\frac{1}{20}\right)^{2} + \left(\frac{\left(\frac{17}{10}\right) \cdot \left(\frac{1}{20}\right)^{2}}{2} - \frac{1}{\left(5\right)}\right)^{2}}\\&\simeq 301.681\,\mathrm{fm}^2.\end{aligned}
$$

$ik$ の項は低エネルギーでも有限の $ka$ では無視できない。

異なる波数の入射を比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
k=\frac{3}{20}\,\mathrm{fm}^{-1}
$$

このときの値を ${\sigma}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\sigma}_{\mathrm{new}}&=\frac{4 \cdot \pi}{\left(\frac{3}{20}\right)^{2} + \left(\frac{\left(\frac{17}{10}\right) \cdot \left(\frac{3}{20}\right)^{2}}{2} - \frac{1}{\left(5\right)}\right)^{2}}\\&\simeq 227.587\,\mathrm{fm}^2,\\{\sigma}_{\mathrm{new}}-{\sigma}&\simeq -74.0942\,\mathrm{fm}^2.\end{aligned}
$$
:::

:::example[極からの散乱長]
束縛極 $\kappa$ と有効距離 $r_e$ が与えられたとする。二次近似の極の条件を $a$ について解く。

量の定義に従い、共通の単位で次の値を用いる。
$$
\kappa=\frac{23}{100}\,\mathrm{fm}^{-1},\qquad r_{e}=0\,\mathrm{fm}
$$

この条件で求める量を ${a}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{a}&=1 \cdot \frac{1}{- \frac{\kappa^{2} \cdot r_{e}}{2} + \kappa}\\&=\frac{1}{- \frac{\left(0\right) \cdot \left(\frac{23}{100}\right)^{2}}{2} + \left(\frac{23}{100}\right)}\\&\simeq 4.34783\,\mathrm{fm}.\end{aligned}
$$

束縛状態が同じ浅さでも有限範囲の補正で散乱長は変わる。

有効距離を無視した場合と比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
r_{e}=\frac{17}{10}\,\mathrm{fm}
$$

このときの値を ${a}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{a}_{\mathrm{new}}&=\frac{1}{- \frac{\left(\frac{17}{10}\right) \cdot \left(\frac{23}{100}\right)^{2}}{2} + \left(\frac{23}{100}\right)}\\&\simeq 5.40438\,\mathrm{fm},\\{a}_{\mathrm{new}}-{a}&\simeq 1.05656\,\mathrm{fm}.\end{aligned}
$$

**一次近似の検証。** 有限距離補正の最初の効果は、零距離模型の周りでの微分から見積もれる。

${r_{e}}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${a}$ の一階微分を求めると、
$$
\frac{\mathrm d{a}}{\mathrm dr_{e}}=\frac{2}{\left(\kappa \cdot r_{e} - 2\right)^{2}}.
$$

この導関数を変更前の点で評価して $h=\frac{17}{10}$ を掛ける。ここで $r_{e}$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq 0.85\,\mathrm{fm},\\\Delta-\Delta_{\mathrm{lin}}&\simeq 0.206557\,\mathrm{fm}.\end{aligned}
$$

残差が示すのは、指定した二点間で生じる線形化の誤差である。基準点からさらに離れた値について同じ精度を保証するものではない。
:::

:::example[硬球の S 波位相]
半径 $R$ で動径波をゼロにする条件を選ぶ。低波数の枝では $\delta_0=-kR$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
R=2\,\mathrm{fm},\qquad k=\frac{1}{10}\,\mathrm{fm}^{-1}
$$

この条件で求める量を ${\delta_0}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\delta_0}&=R \cdot \left(- k\right)\\&=- \left(2\right) \cdot \left(\frac{1}{10}\right)\\&\simeq -0.2\,\mathrm{rad}.\end{aligned}
$$

反発的な硬球では負の位相になる。位相は $\pi$ の整数倍の違いを同じ散乱として表す。

エネルギーの高い側を比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
k=\frac{1}{4}\,\mathrm{fm}^{-1}
$$

このときの値を ${\delta_0}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\delta_0}_{\mathrm{new}}&=- \left(2\right) \cdot \left(\frac{1}{4}\right)\\&\simeq -0.5\,\mathrm{rad},\\{\delta_0}_{\mathrm{new}}-{\delta_0}&\simeq -0.3\,\mathrm{rad}.\end{aligned}
$$
:::

:::example[ユニタリティの上限]
弾性 $S$ 波の断面積は $4\pi\sin^2\delta_0/k^2$ である。位相を任意に許したときの最大値を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
k=\frac{1}{10}\,\mathrm{fm}^{-1}
$$

この条件で求める量を ${\sigma_{\rm max}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\sigma_{\rm max}}&=\frac{4 \cdot \pi}{k^{2}}\\&=\frac{4 \cdot \pi}{\left(\frac{1}{10}\right)^{2}}\\&\simeq 1256.64\,\mathrm{fm}^2.\end{aligned}
$$

実測断面積はこの上限以下であり、上限を達成するには共鳴などで位相が適切な値になる必要がある。

波数を二倍にする。そこで次の値に変更し、ほかの入力を固定する。
$$
k=\frac{1}{5}\,\mathrm{fm}^{-1}
$$

このときの値を ${\sigma_{\rm max}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\sigma_{\rm max}}_{\mathrm{new}}&=\frac{4 \cdot \pi}{\left(\frac{1}{5}\right)^{2}}\\&\simeq 314.159\,\mathrm{fm}^2,\\{\sigma_{\rm max}}_{\mathrm{new}}-{\sigma_{\rm max}}&\simeq -942.478\,\mathrm{fm}^2.\end{aligned}
$$
:::

## 核力の荷電独立性とアイソスピン

:ref[def-nuc-isospin] の回転すべてに対して $V$ が不変であることを荷電独立性という。$p$ と $n$ を交換する一つの回転についての不変性は荷電対称性である。二核子のアイソスピン Pauli ベクトルを $\vb\tau_1,\vb\tau_2$ とする。$\vb T=(\vb\tau_1+\vb\tau_2)/2$ の二乗から $\vb\tau_1\cdot\vb\tau_2=2T(T+1)-3$ となる。従って射影子は
$$
P_{T=0}=(1-\vb\tau_1\cdot\vb\tau_2)/4,\qquad
P_{T=1}=(3+\vb\tau_1\cdot\vb\tau_2)/4.
$$
各射影子が固有値 0 または 1 を持つことは $T=0,1$ を代入すれば分かる。$v+u\vb\tau_1\cdot\vb\tau_2$ は荷電独立だが、各 $T$ の力は異なる。

:::example[荷電独立の行列要素]
$v=0,u=2\,\mathrm{MeV}$ なら $T=0$ で $-6\,\mathrm{MeV}$、$T=1$ で $2\,\mathrm{MeV}$ である。後者の値は $T_z=1,0,-1$ で共通である。
:::

:::example[対称性の区別]
$\tau_{1z}\tau_{2z}$ は両方の $z$ 成分を反転しても不変である。$T=1$ では $nn,pp$ に固有値 $+1$、対称な $np$ に $-1$ を与えるので、荷電対称性は満たすが荷電独立性は満たさない。実験で比較するときは陽子同士の Coulomb 散乱を扱う必要がある。
:::

### 多重項内の対称性の破れ

同じ軌道・スピンの状態で、荷電独立な核力は $T=0,1$ の空間にそれぞれ定数として作用する。従って $T_z=0$ の積状態 $\ket{np}$ に戻すと平均 $(V_1+V_0)/2$ が対角項、差 $(V_1-V_0)/2$ が電荷交換項になる。$\ket{nn}$ と $\ket{pp}$ は同じ $V_1$ である。これらの関係は Coulomb 力を含めた測定値にそのまま適用できない。

荷電対称性のみを課すと $T_z\mapsto-T_z$ の対称性は残る。したがって $T=1$ の一次摂動を $A+B T_z+C(3T_z^2-2)$ と書くと $B=0$ だが、$C$ は消えなくてもよい。

:::example[積基底の対角成分]
同じ動径状態の $T=0,1$ における行列要素を $V_0,V_1$ とする。$\ket{np}$ は半分ずつを含む。

量の定義に従い、共通の単位で次の値を用いる。
$$
V_{1}=-4\,\mathrm{MeV},\qquad V_{0}=-8\,\mathrm{MeV}
$$

この条件で求める量を ${V_{np}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{V_{np}}&=\frac{V_{0} + V_{1}}{2}\\&=\frac{\left(-4\right)}{2} + \frac{\left(-8\right)}{2}\\&=-6\,\mathrm{MeV}.\end{aligned}
$$

これはアイソスピン状態だけを積基底に書き直した行列要素である。Pauli 条件は別に確認する。

一重項方向の行列要素を変える。そこで次の値に変更し、ほかの入力を固定する。
$$
V_{0}=-10\,\mathrm{MeV}
$$

このときの値を ${V_{np}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{V_{np}}_{\mathrm{new}}&=\frac{\left(-10\right)}{2} + \frac{\left(-4\right)}{2}\\&=-7\,\mathrm{MeV},\\{V_{np}}_{\mathrm{new}}-{V_{np}}&=-1\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[電荷交換成分]
対称・反対称の固有基底のうち後者の符号が交換操作で反転する。$\ket{np}$ から $\ket{pn}$ への行列要素を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
V_{1}=-4\,\mathrm{MeV},\qquad V_{0}=-8\,\mathrm{MeV}
$$

この条件で求める量を ${V_{np,pn}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{V_{np,pn}}&=\frac{- V_{0} + V_{1}}{2}\\&=\frac{\left(-4\right)}{2} - \frac{\left(-8\right)}{2}\\&=2\,\mathrm{MeV}.\end{aligned}
$$

両作用が等しければ積基底でも対角化されて交換項は消える。

二つのアイソスピンの力を同じにする。そこで次の値に変更し、ほかの入力を固定する。
$$
V_{0}=-4\,\mathrm{MeV}
$$

このときの値を ${V_{np,pn}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{V_{np,pn}}_{\mathrm{new}}&=0\\&=0\,\mathrm{MeV},\\{V_{np,pn}}_{\mathrm{new}}-{V_{np,pn}}&=-2\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[荷電対称性の破れ]
二核子の $T=1$ 多重項で端の行列要素を $V_{nn},V_{pp}$ とする。$B$ はこの差の半分である。

量の定義に従い、共通の単位で次の値を用いる。
$$
V_{nn}=- \frac{41}{10}\,\mathrm{MeV},\qquad V_{pp}=- \frac{39}{10}\,\mathrm{MeV}
$$

この条件で求める量を ${B}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{B}&=\frac{V_{nn} - V_{pp}}{2}\\&=- \frac{\left(- \frac{39}{10}\right)}{2} + \frac{\left(- \frac{41}{10}\right)}{2}\\&\simeq -0.1\,\mathrm{MeV}.\end{aligned}
$$

端だけから分かるのは奇成分であり、中央との関係を求めるには追加の情報が必要である。

端を等しくする。そこで次の値に変更し、ほかの入力を固定する。
$$
V_{pp}=- \frac{41}{10}\,\mathrm{MeV}
$$

このときの値を ${B}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{B}_{\mathrm{new}}&=0\\&=0\,\mathrm{MeV},\\{B}_{\mathrm{new}}-{B}&\simeq 0.1\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[荷電独立性の二次の破れ]
端の平均を $V_e$、中央を $V_c$ とする。端で $3T_z^2-2=1$、中央で $-2$ なので二次係数を求められる。

量の定義に従い、共通の単位で次の値を用いる。
$$
V_{e}=-4\,\mathrm{MeV},\qquad V_{c}=- \frac{43}{10}\,\mathrm{MeV}
$$

この条件で求める量を ${C}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{C}&=\frac{- V_{c} + V_{e}}{3}\\&=- \frac{\left(- \frac{43}{10}\right)}{3} + \frac{\left(-4\right)}{3}\\&\simeq 0.1\,\mathrm{MeV}.\end{aligned}
$$

端が互いに等しくても、この係数が非ゼロなら全回転の対称性は破れている。

三成分が等しい場合と比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
V_{c}=-4\,\mathrm{MeV}
$$

このときの値を ${C}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{C}_{\mathrm{new}}&=0\\&=0\,\mathrm{MeV},\\{C}_{\mathrm{new}}-{C}&\simeq -0.1\,\mathrm{MeV}.\end{aligned}
$$
:::

## 核力のスピン依存性と交換性

スピンでも $\vb\sigma_1\cdot\vb\sigma_2=2S(S+1)-3$ であり、スピン交換 $P_\sigma$ は三重項で $+1$、一重項で $-1$ となる。よって $P_\sigma=(1+\vb\sigma_1\cdot\vb\sigma_2)/2$ である。アイソスピン交換 $P_\tau$ も同形である。空間交換を $P_r$ とすると、反対称化部分空間で $P_rP_\sigma P_\tau=-1$ である。交換演算子を含む力は、軌道やスピンに応じて異なる行列要素を与える。

中心スピン項に加えて、回転やパリティを保つ演算子には $S_{12}$ と $\vb L\cdot\vb S$ がある。後者の固有値は $\hbar^2[J(J+1)-L(L+1)-S(S+1)]/2$ である。テンソル項は同じ $J,S$ の中で $L$ を混ぜ、スピン軌道項は与えた $L,S$ の $J$ を分裂させる。

### 中心交換力の各成分

中心力の強さをスピンとアイソスピンに応じて変えるには、$V_C(r)=A(r)+B(r)P_\sigma+C(r)P_\tau+D(r)P_\sigma P_\tau$ と展開できる。$P_\sigma,P_\tau$ の固有値を $s,t=\pm1$ とすると、チャネルの値は $V_{st}=A+Bs+Ct+Dst$ である。四つの符号の組について足せば $A=\sum_{s,t}V_{st}/4$ であり、$s,t,st$ を掛けてから足せばそれぞれ $B,C,D$ が得られる。これはチャネルごとに異なる中心力と、交換演算子で書く中心力との変換である。

物理的な反対称部分空間では $P_\sigma P_\tau=-P_r$、$P_\tau=-P_rP_\sigma$ なので、同じ力を $1,P_r,P_\sigma,P_rP_\sigma$ の線形結合としても表せる。空間交換 $P_r$ は相対座標を反転し、$L$ 波に $(-1)^L$ を掛ける。従ってこの表示では偶数 $L$ と奇数 $L$ の強さがどの係数で異なるかが見える。符号を比べる際は、交換演算子の定義と反対称化の前後をそろえる。

:::example[三重項 P 波]
$L=S=1$ では $\vb L\cdot\vb S/\hbar^2$ が $J=0,1,2$ に対して $-2,-1,+1$ となる。縮退度 $1,3,5$ で重み付けした和は $-2-3+5=0$ である。
:::

:::example[非偏極散乱]
二つのスピン $1/2$ を偏極せず入射させると、4 状態のうち 1 状態が一重項、3 状態が三重項である。区別可能な $np$ の低エネルギー散乱で両 S 波が寄与するなら $\sigma\simeq(\sigma_s+3\sigma_t)/4$ となる。振幅の重ね合わせでなく、初期スピン状態に対する確率平均である。
:::

:::example[中性子・陽子の低エネルギー断面積]
散乱長の代表値を $a_s=-23.7\,\mathrm{fm}$、$a_t=5.42\,\mathrm{fm}$ とする。しきい値のスピン平均は $\pi(a_s^2+3a_t^2)\simeq2.04\times10^3\,\mathrm{fm}^2$、すなわち約 $20.4\,\mathrm b$ である。一重項の散乱長の絶対値が大きいため、初期状態の重みが $1/4$ でもその寄与が大きい。電磁相互作用を除いた短距離模型の散乱長は、Coulomb 力を含む陽子同士のしきい値に直接代入できない。散乱データと重陽子を組み合わせた核力の決定については :cite[krane-nuclear] を参照する。
:::

### スピン交換の射影と偏極

二つのスピン $1/2$ について $\vb S^2/\hbar^2=3/2+\vb\sigma_1\cdot\vb\sigma_2/2$ である。したがって三重項と一重項への射影子は $(3+\vb\sigma_1\cdot\vb\sigma_2)/4$ と $(1-\vb\sigma_1\cdot\vb\sigma_2)/4$ である。それぞれ対象の固有空間で 1、他方で 0 になることが証明になる。無偏極の積状態の平均は四次元の単位行列を 4 で割ったものなので、スピンの重みは $3/4,1/4$ となる。反対称化した同種粒子の部分波にはこの重みを無条件には使えない。

:::example[スピン平均の中心力]
区別可能な非偏極粒子について、一重項・三重項の断面積を $\sigma_0,\sigma_1$ とする。

量の定義に従い、共通の単位で次の値を用いる。
$$
\sigma_{0}=8\,\mathrm{b},\qquad \sigma_{1}=2\,\mathrm{b}
$$

この条件で求める量を ${\sigma}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\sigma}&=\frac{\sigma_{0}}{4} + \frac{3 \cdot \sigma_{1}}{4}\\&=\frac{3 \cdot \left(2\right)}{4} + \frac{\left(8\right)}{4}\\&\simeq 3.5\,\mathrm{b}.\end{aligned}
$$

平均するのは振幅の平均の二乗ではなく、直交する初期状態ごとの確率である。

三重項が強いエネルギーと比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
\sigma_{1}=4\,\mathrm{b}
$$

このときの値を ${\sigma}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\sigma}_{\mathrm{new}}&=\frac{3 \cdot \left(4\right)}{4} + \frac{\left(8\right)}{4}\\&=5\,\mathrm{b},\\{\sigma}_{\mathrm{new}}-{\sigma}&\simeq 1.5\,\mathrm{b}.\end{aligned}
$$
:::

:::example[反平行の偏極]
独立な二つの偏極ベクトルの軸方向成分を $P_1,P_2$ とする。射影子の期待値で $\langle\vb\sigma_1\cdot\vb\sigma_2\rangle=P_1P_2$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
P_{1}=1,\qquad P_{2}=1
$$

この条件で求める量を ${f_0}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{f_0}&=\frac{1 - P_{1} \cdot P_{2}}{4}\\&=\frac{1}{4} - \frac{\left(1\right)^{2}}{4}\\&=0.\end{aligned}
$$

反平行の積状態は一重項と三重項の $M=0$ の重ね合わせである。

片側だけの偏極を反転する。そこで次の値に変更し、ほかの入力を固定する。
$$
P_{2}=-1
$$

このときの値を ${f_0}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{f_0}_{\mathrm{new}}&=- \frac{\left(-1\right) \cdot \left(1\right)}{4} + \frac{1}{4}\\&\simeq 0.5,\\{f_0}_{\mathrm{new}}-{f_0}&\simeq 0.5.\end{aligned}
$$
:::

:::example[スピン係数の抽出]
中心力を $V_c+V_\sigma\vb\sigma_1\cdot\vb\sigma_2$ と書く。二つの固有値 $V_1,V_0$ を引いて $V_\sigma$ を求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
V_{1}=-5\,\mathrm{MeV},\qquad V_{0}=-1\,\mathrm{MeV}
$$

この条件で求める量を ${V_\sigma}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{V_\sigma}&=\frac{- V_{0} + V_{1}}{4}\\&=- \frac{\left(-1\right)}{4} + \frac{\left(-5\right)}{4}\\&=-1\,\mathrm{MeV}.\end{aligned}
$$

この符号では負の係数が三重項を相対的に低くする。

一重項の引力を強めた場合と比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
V_{0}=-3\,\mathrm{MeV}
$$

このときの値を ${V_\sigma}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{V_\sigma}_{\mathrm{new}}&=- \frac{\left(-3\right)}{4} + \frac{\left(-5\right)}{4}\\&\simeq -0.5\,\mathrm{MeV},\\{V_\sigma}_{\mathrm{new}}-{V_\sigma}&\simeq 0.5\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[P 波のスピン軌道分裂]
$L=S=1$ で係数を $a$、Hamiltonian の項を $a\vb L\cdot\vb S/\hbar^2$ とする。

量の定義に従い、共通の単位で次の値を用いる。
$$
a=\frac{3}{10}\,\mathrm{MeV},\qquad J=0
$$

この条件で求める量を ${E_J}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{E_J}&=\frac{a \cdot \left(J \cdot \left(J + 1\right) - 4\right)}{2}\\&=\frac{\left(\frac{3}{10}\right) \cdot \left(\left(0\right) \cdot \left(\left(0\right) + 1\right) - 4\right)}{2}\\&\simeq -0.6\,\mathrm{MeV}.\end{aligned}
$$

同じ $L,S$ にいても異なる $J$ の位置は分かれる。$2J+1$ で重み付けすると三枝のエネルギーの和はゼロである。

同じ三重項の最大 $J$ と比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
J=2
$$

このときの値を ${E_J}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{E_J}_{\mathrm{new}}&=\frac{\left(\frac{3}{10}\right) \cdot \left(\left(2\right) \cdot \left(\left(2\right) + 1\right) - 4\right)}{2}\\&\simeq 0.3\,\mathrm{MeV},\\{E_J}_{\mathrm{new}}-{E_J}&\simeq 0.9\,\mathrm{MeV}.\end{aligned}
$$
:::

## 湯川理論と一 π 中間子交換ポテンシャル

質量 $m_\pi$ の場の静的伝播を考え、$\lambda_\pi=m_\pi c/\hbar$ と置く。$(-\nabla^2+\lambda_\pi^2)Y=\delta^3(\vb r)$ の解で無限遠に減衰する球対称関数は $Y=e^{-\lambda_\pi r}/(4\pi r)$ である。実際、原点外では $rY$ が定係数の二階方程式を満たす。小球の境界から流出する $-\nabla Y$ の積分が 1 になるよう係数が決まる。

pion は擬スカラーのアイソスピン三重項である。低運動量の各頂点に $\vb\sigma\cdot\vb q$ と $\vb\tau$ を含む微分結合を置くと、運動量空間で長距離ポテンシャルは $-C(\vb\tau_1\cdot\vb\tau_2)(\vb\sigma_1\cdot\vb q)(\vb\sigma_2\cdot\vb q)/(q^2+\lambda_\pi^2)$ となる。$C>0$ は規約に応じた結合係数である。Fourier 変換で $q_iq_j$ を $-\partial_i\partial_j$ に替え、球対称 Hessian を $(Y''-Y'/r)\hat r_i\hat r_j+(Y'/r)\delta_{ij}$ と分けると
$$
V_\pi=\frac C3(\vb\tau_1\cdot\vb\tau_2)
\left[(\vb\sigma_1\cdot\vb\sigma_2)\nabla^2Y+S_{12}(Y''-Y'/r)\right].
$$
原点外では $\nabla^2Y=\lambda_\pi^2Y$、$Y''-Y'/r=(\lambda_\pi^2+3\lambda_\pi/r+3/r^2)Y$ である。原点の接触項と短距離の特異性の扱いには正則化と短距離項が必要である。

:::example[到達距離]
$m_\pi c^2\simeq140\,\mathrm{MeV}$ なら $1/\lambda_\pi\simeq1.41\,\mathrm{fm}$ である。同じ結合形で交換質量を倍にすれば、指数減衰長は半分になる。
:::

:::example[一重項]
スピン一重項では任意の軸で $(\vb\sigma_1\cdot\hat{\vb r})(\vb\sigma_2\cdot\hat{\vb r})=-1$、$\vb\sigma_1\cdot\vb\sigma_2=-3$ なので $S_{12}$ は消える。長距離のスピン・スピン成分は残る。
:::

### 静的な伝播の範囲と接触項

空間 Fourier 変換の規約を $Y(\vb r)=\int e^{i\vb q\cdot\vb r}(q^2+\lambda_\pi^2)^{-1}\dd^3q/(2\pi)^3$ とする。$\nabla^2Y$ の変換は $-q^2/(q^2+\lambda_\pi^2)=-1+\lambda_\pi^2/(q^2+\lambda_\pi^2)$ である。定数 $-1$ は座標空間のデルタ関数 $-\delta^3(\vb r)$ なので、原点外だけの式に落とすと接触項を失う。一方、短距離には他の相互作用もある。低エネルギー有効理論で接触係数を独立に定める場合は、どの正則化と演算子の規約で交換項と接触項を分けたかを揃える。

:::example[滑らかな点源]
規格化した $\rho_a(\vb r)=(\pi a^2)^{-3/2}e^{-r^2/a^2}$ をデルタ関数の代わりに使う。波数空間では $e^{-a^2q^2/4}$ が掛かり、$qa\ll1$ の補正は $-a^2q^2/4$ から始まる。同じ低波数物理を保つには、短距離を変更した分を同程度の次数の接触項に含める。
:::

### 伝播関数と長距離・短距離の分離

$Y(r)=e^{-\lambda r}/(4\pi r)$ は $r>0$ で $\nabla^2Y=\lambda^2Y$ を満たす。原点を囲む球の法線フラックス $4\pi\epsilon^2Y'(\epsilon)$ は $\epsilon\downarrow0$ で $-1$ なので、分布として $(\nabla^2-\lambda^2)Y=-\delta^{(3)}$ である。この式を Fourier 変換すると $\widetilde Y=(q^2+\lambda^2)^{-1}$ となる。スピン頂点から出る二微分を取ると、軌道方向のテンソル部分とラプラシアンのスカラー部分が分かれる。接触項と長距離尾を区別するには原点のフラックスを落とさないことが重要である。

:::example[中間子質量と到達距離]
交換粒子の静止エネルギー $m_c$ から $1/\lambda$ を fm で求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
h_{c}=\frac{197327}{1000}\,\mathrm{MeV\,fm},\qquad m_{c}=140\,\mathrm{MeV}
$$

この条件で求める量を ${L}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{L}&=\frac{h_{c}}{m_{c}}\\&=\frac{\left(\frac{197327}{1000}\right)}{\left(140\right)}\\&\simeq 1.40948\,\mathrm{fm}.\end{aligned}
$$

この距離はポテンシャルが急に切れる半径ではなく、指数因子が $e^{-1}$ になる長さである。

重い交換粒子に変える。そこで次の値に変更し、ほかの入力を固定する。
$$
m_{c}=500\,\mathrm{MeV}
$$

このときの値を ${L}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{L}_{\mathrm{new}}&=\frac{\left(\frac{197327}{1000}\right)}{\left(500\right)}\\&\simeq 0.394654\,\mathrm{fm},\\{L}_{\mathrm{new}}-{L}&\simeq -1.01482\,\mathrm{fm}.\end{aligned}
$$
:::

:::example[同じ距離での尾の比]
結合係数を同じと仮定し、二つの交換の質量差を $\Delta_m>0$ とする。消える共通因子は $1/r$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
\Delta_{m}=360\,\mathrm{MeV},\qquad r=1\,\mathrm{fm},\qquad h_{c}=\frac{197327}{1000}\,\mathrm{MeV\,fm}
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=e^{\frac{- \Delta_{m} \cdot r}{h_{c}}}\\&=e^{- \frac{\left(1\right) \cdot \left(360\right)}{\left(\frac{197327}{1000}\right)}}\\&\simeq 0.161317.\end{aligned}
$$

重い交換が急速に抑えられることを示す。結合定数まで同じとみなす必要は一般にはない。

より外側の距離で比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
r=2\,\mathrm{fm}
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=e^{- \frac{\left(2\right) \cdot \left(360\right)}{\left(\frac{197327}{1000}\right)}}\\&\simeq 0.0260232,\\{R}_{\mathrm{new}}-{R}&\simeq -0.135294.\end{aligned}
$$
:::

:::example[運動量空間の相対振幅]
零運動量で規格化した伝播関数の比を用いる。$q$ と $\lambda$ は逆長さである。

量の定義に従い、共通の単位で次の値を用いる。
$$
\lambda_{0}=\frac{7}{10}\,\mathrm{fm}^{-1},\qquad q=\frac{2}{5}\,\mathrm{fm}^{-1}
$$

この条件で求める量を ${F}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{F}&=\frac{\lambda_{0}^{2}}{\lambda_{0}^{2} + q^{2}}\\&=\frac{\left(\frac{7}{10}\right)^{2}}{\left(\frac{2}{5}\right)^{2} + \left(\frac{7}{10}\right)^{2}}\\&\simeq 0.753846.\end{aligned}
$$

頂点の運動量因子は含めていない。低 $q$ ではこの比を $1-q^2/\lambda_0^2+\cdots$ と展開できる。

逆範囲より大きい波数へ移す。そこで次の値に変更し、ほかの入力を固定する。
$$
q=\frac{7}{5}\,\mathrm{fm}^{-1}
$$

このときの値を ${F}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{F}_{\mathrm{new}}&=\frac{\left(\frac{7}{10}\right)^{2}}{\left(\frac{7}{10}\right)^{2} + \left(\frac{7}{5}\right)^{2}}\\&\simeq 0.2,\\{F}_{\mathrm{new}}-{F}&\simeq -0.553846.\end{aligned}
$$
:::

:::example[スピンテンソルと指数因子]
非接触部分の二階微分で、$Y''-Y'/r$ を $\lambda^2Y$ で割る。$x=\lambda r$ を使う。

量の定義に従い、共通の単位で次の値を用いる。
$$
x=1
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=1 + \frac{3}{x^{2}} + \frac{3}{x}\\&=1 + \frac{3}{\left(1\right)} + \frac{3}{\left(1\right)^{2}}\\&=7.\end{aligned}
$$

遠方で比は 1 に近づく。近接域では二微分がより急な距離依存を作る。

長距離側と比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
x=3
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=1 + \frac{3}{\left(3\right)} + \frac{3}{\left(3\right)^{2}}\\&\simeq 2.33333,\\{R}_{\mathrm{new}}-{R}&\simeq -4.66667.\end{aligned}
$$
:::

## 現実的核力と斥力芯

散乱位相と重陽子の観測量を同時に記述するため、現実的相互作用は長距離の pion 交換と短距離の演算子を組み合わせる。強い短距離斥力をもつ局所ポテンシャルでは、低エネルギー波動関数が小さい $r$ で抑えられる。硬い斥力芯は極限模型であり、位相差からポテンシャルの短距離形が一意に決まるとは限らない。

例えば Argonne $v_{18}$ は、中心、スピン、テンソル、スピン軌道などの演算子に局所的な動径係数を割り当て、荷電独立性を破る成分も持たせる。対して運動量空間で組み立てる有効相互作用では、入射と出射の運動量に別々に依存する非局所項を使うことがある。両者を比較するときは、同じエネルギー範囲での散乱と束縛状態の再現を先に比べる。短距離の相対波動関数を比較するには、ポテンシャルと観測演算子の規約もそろえる。

:::example[芯と境界条件]
半径 $R_c$ の内部に無限斥力を置くと $u(R_c)=0$ が外部問題の境界条件である。さらに外部で $V=0$ なら $u_0=r-R_c$ であり、散乱長はちょうど $R_c$ になる。
:::

:::example[表現と観測]
漸近領域で恒等に戻るユニタリ変換 $U$ で $H'=UHU^\dagger$ とすると、固有値と変換後の漸近散乱波は共通である。内部波動関数は $U$ によって変わる。観測演算子も $O'=UOU^\dagger$ とすれば $\bra{U\psi}O'\ket{U\psi}=\bra\psi O\ket\psi$ となる。
:::

### 位相同値と短距離相関

弾性の位相差は長距離の散乱状態の係数を固定する。短距離だけで働き、漸近領域では恒等となるユニタリ変換 $U$ で $H'=UHU^\dagger$ を作ると、$H\psi=E\psi$ に対して $H'(U\psi)=E(U\psi)$ であり、漸近の係数は保存される。従って同じ散乱データでも内部の波動関数を一意に決められない。内部に敏感な演算子は $O'=UOU^\dagger$ に変換する。このとき $\langle U\psi|O'|U\psi\rangle=\langle\psi|O|\psi\rangle$ である。

:::example[相関孔の規格化]
短距離因子を $f(r)=1-\exp[-r^2/d^2]$ と仮定する。同じ相対座標体積要素の確率比は $f^2$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
d=\frac{4}{5}\,\mathrm{fm},\qquad r=\frac{2}{5}\,\mathrm{fm}
$$

この条件で求める量を ${P_{\rm ratio}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{P_{\rm ratio}}&=\left(1 - e^{- \left(\frac{r}{d}\right)^{2}}\right)^{2}\\&=\left(1 - e^{- \frac{\left(\frac{2}{5}\right)^{2}}{\left(\frac{4}{5}\right)^{2}}}\right)^{2}\\&\simeq 0.0489291.\end{aligned}
$$

分布全体を規格化し直す前の局所的な比である。短距離の抑制は運動量分布にも広がりを生む。

孔の外側と比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
r=\frac{8}{5}\,\mathrm{fm}
$$

このときの値を ${P_{\rm ratio}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{P_{\rm ratio}}_{\mathrm{new}}&=\left(1 - e^{- \frac{\left(\frac{8}{5}\right)^{2}}{\left(\frac{4}{5}\right)^{2}}}\right)^{2}\\&\simeq 0.963704,\\{P_{\rm ratio}}_{\mathrm{new}}-{P_{\rm ratio}}&\simeq 0.914775.\end{aligned}
$$
:::

:::example[接触項の分解能]
Gaussian 正則化の運動量因子を $f(q)=\exp(-q^2/\Lambda^2)$ とする。

量の定義に従い、共通の単位で次の値を用いる。
$$
q=2\,\mathrm{fm}^{-1},\qquad \Lambda=4\,\mathrm{fm}^{-1}
$$

この条件で求める量を ${f}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{f}&=e^{- \left(\frac{q}{\Lambda}\right)^{2}}\\&=e^{- \frac{\left(2\right)^{2}}{\left(4\right)^{2}}}\\&\simeq 0.778801.\end{aligned}
$$

カットオフを変えれば、物理量を保つために係数も同じ規約で再評価する。

低いカットオフで比較する。そこで次の値に変更し、ほかの入力を固定する。
$$
\Lambda=2\,\mathrm{fm}^{-1}
$$

このときの値を ${f}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{f}_{\mathrm{new}}&=e^{-1}\\&\simeq 0.367879,\\{f}_{\mathrm{new}}-{f}&\simeq -0.410921.\end{aligned}
$$

**一次近似の検証。** 運動量を固定してカットオフのみ変えると、正則化を弱める方向の局所的な効果を調べられる。

${\Lambda}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${f}$ の一階微分を求めると、
$$
\frac{\mathrm d{f}}{\mathrm d\Lambda}=\frac{2 \cdot q^{2} \cdot e^{- \frac{q^{2}}{\Lambda^{2}}}}{\Lambda^{3}}.
$$

この導関数を変更前の点で評価して $h=-2$ を掛ける。ここで $\Lambda$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq -0.1947,\\\Delta-\Delta_{\mathrm{lin}}&\simeq -0.216221.\end{aligned}
$$

省いた分が全変化に比べて大きい。従って増分の大きさを予測するには、少なくとも二次の項を調べるか、元の式を使う必要がある。
:::

:::example[二階の摂動によるエネルギー変化]
低い状態が間隔 $d>0$ の高い状態と $v$ で結合する。$v/d$ が小さければ相関エネルギーは $-v^2/d$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
v=1\,\mathrm{MeV},\qquad d=10\,\mathrm{MeV}
$$

この条件で求める量を ${E^{(2)}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{E^{(2)}}&=\frac{\left(-1\right) \cdot v^{2}}{d}\\&=- \frac{\left(1\right)^{2}}{\left(10\right)}\\&\simeq -0.1\,\mathrm{MeV}.\end{aligned}
$$

中間状態のエネルギーが高いほど、同じ行列要素の効果は小さい。

除外した空間がさらに離れた場合と比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
d=20\,\mathrm{MeV}
$$

このときの値を ${E^{(2)}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{E^{(2)}}_{\mathrm{new}}&=- \frac{\left(1\right)^{2}}{\left(20\right)}\\&\simeq -0.05\,\mathrm{MeV},\\{E^{(2)}}_{\mathrm{new}}-{E^{(2)}}&\simeq 0.05\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[硬芯と短距離位相]
硬球境界での零エネルギー外部解を $u\propto r-a$ と書く。節が $R$ なので区別可能な粒子の断面積は $4\pi R^2$ である。

量の定義に従い、共通の単位で次の値を用いる。
$$
R=\frac{2}{5}\,\mathrm{fm}
$$

この条件で求める量を ${\sigma}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\sigma}&=4 \cdot \pi \cdot R^{2}\\&=4 \cdot \pi \cdot \left(\frac{2}{5}\right)^{2}\\&\simeq 2.01062\,\mathrm{fm}^2.\end{aligned}
$$

現実の散乱長には引力による尾も寄与する。模型の半径を直接観測したことにはならない。

硬芯の大きさを変える。そこで次の値に変更し、ほかの入力を固定する。
$$
R=\frac{3}{5}\,\mathrm{fm}
$$

このときの値を ${\sigma}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\sigma}_{\mathrm{new}}&=4 \cdot \pi \cdot \left(\frac{3}{5}\right)^{2}\\&\simeq 4.52389\,\mathrm{fm}^2,\\{\sigma}_{\mathrm{new}}-{\sigma}&\simeq 2.51327\,\mathrm{fm}^2.\end{aligned}
$$
:::

## 三体力の概観

三核子が同時に近づくときのエネルギーを、各対の行列要素の和から独立に選べる項を三体力という。二つの pion を介して三核子が結びつき、中央の核子が中間励起状態を経る模型が代表例である。低エネルギー有効理論では二体の指標だけで全係数が決まらない三体演算子が許される。次数別の全演算子と正則化条件は有効場の理論の導出を要するため、本章では三体項の識別に絞る。

:::example[組合せ]
粒子が 4 個なら対は $\binom42=6$、三つ組は $\binom43=4$ である。全粒子が同じ点にいるという仮想的な配置でも、同じ数え方の係数で二体と三体を取り替えることは一般の $A$ に対してできない。
:::

:::example[二体データで決まらない項]
$W=w\exp[-(r_{12}^2+r_{23}^2+r_{31}^2)/b^2]$ を全三つ組について足す。どれか一粒子を遠ざけるとその三つ組の寄与は消え、二核子系にはそもそも三つ組がない。従って二核子の束縛エネルギーや位相差を固定したまま $w$ を変えられる。
:::

### 多体系に埋め込む際の整合性

三粒子の Hamiltonian を $H=T+\sum_{i<j}V_{ij}$ とする。対ごとの短距離のユニタリ変換 $U=e^X$、$X=\sum_{i<j}X_{ij}$ を全体系に作用させれば、一次の変化は $[X,H]$ である。ここには $[X_{12},V_{23}]$ のように粒子 2 を共有する交換子が含まれる。これは一般に三粒子の座標やスピンに作用する。対をそれぞれ二粒子系で変換して足し直す操作では、この項が入らない。二体散乱で同じ性質を持つ相互作用を使っても、三体系の予測には誘起された多体項を揃える必要がある。

:::example[共有する粒子の交換子]
$X_{12}=i\epsilon(\sigma_{1z}\sigma_{2z})$、$V_{23}=v\,\sigma_{2x}\sigma_{3x}$ とすると、$\sigma_{2z}$ と $\sigma_{2x}$ の交換子は $2i\sigma_{2y}$ である。従って $[X_{12},V_{23}]=-2\epsilon v\,\sigma_{1z}\sigma_{2y}\sigma_{3x}$ となる。空間因子を省いた演算子の例だが、三つの粒子に同時に作用する項が生じることを示す。
:::

### 三体演算子が二体力の変更で生じる理由

ユニタリ変換を $U=e^A$ とし、$A^\dagger=-A$ とする。$e^{tA}He^{-tA}$ の $t=0$ の微分は $[A,H]$ である。二体の $A_{12}$ と二体の $V_{23}$ は粒子 2 の自由度を共有し、一般に可換でない。例えばスピンで $[\sigma_{1x}\sigma_{2x},\sigma_{2y}\sigma_{3y}]=2i\sigma_{1x}\sigma_{2z}\sigma_{3y}$ となる。作用する粒子が三つなので、二体ポテンシャルの変換を多体系に持ち込むだけでは、この寄与が欠ける。

:::example[三つ組の数]
粒子の順序を区別せず $A$ 個から三個を取る。

量の定義に従い、共通の単位で次の値を用いる。
$$
A=3
$$

この条件で求める量を ${n_3}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{n_3}&=\frac{A \cdot \left(A - 1\right) \cdot \left(A - 2\right)}{6}\\&=\frac{\left(3\right) \cdot \left(\left(3\right) - 2\right) \cdot \left(\left(3\right) - 1\right)}{6}\\&=1.\end{aligned}
$$

三つ組の数と三体相関の強度は別の量である。局所力では離れた三粒子の寄与は小さい。

粒子数を増やした場合を数える。そこで次の値に変更し、ほかの入力を固定する。
$$
A=6
$$

このときの値を ${n_3}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{n_3}_{\mathrm{new}}&=\frac{\left(6\right) \cdot \left(\left(6\right) - 2\right) \cdot \left(\left(6\right) - 1\right)}{6}\\&=20,\\{n_3}_{\mathrm{new}}-{n_3}&=19.\end{aligned}
$$
:::

:::example[二体と三体の組合せ比]
すべての組を等しく数えると、三つ組の数を対の数で割った比は $(A-2)/3$ になる。

量の定義に従い、共通の単位で次の値を用いる。
$$
A=5
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=\frac{A - 2}{3}\\&=\frac{\left(5\right)}{3} - \frac{2}{3}\\&=1.\end{aligned}
$$

この比が増えても、短距離だけで働く力のエネルギー比が同様に増えるとはいえない。

粒子数の異なる模型と比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
A=8
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\frac{\left(8\right)}{3} - \frac{2}{3}\\&=2,\\{R}_{\mathrm{new}}-{R}&=1.\end{aligned}
$$
:::

:::example[一様密度の三体寄与]
密度依存のエネルギー密度を $C_3\rho^3$ とする。粒子当たりの値は $\rho$ で割る。

量の定義に従い、共通の単位で次の値を用いる。
$$
C_{3}=100\,\mathrm{MeV\,fm^6},\qquad \rho=\frac{2}{25}\,\mathrm{fm}^{-3}
$$

この条件で求める量を ${e_3}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{e_3}&=C_{3} \cdot \rho^{2}\\&=\left(100\right) \cdot \left(\frac{2}{25}\right)^{2}\\&\simeq 0.64\,\mathrm{MeV}.\end{aligned}
$$

二体の接触寄与が粒子当たり $\rho$ に比例する模型とは密度依存が異なる。

密度を二倍にする。そこで次の値に変更し、ほかの入力を固定する。
$$
\rho=\frac{4}{25}\,\mathrm{fm}^{-3}
$$

このときの値を ${e_3}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{e_3}_{\mathrm{new}}&=\left(100\right) \cdot \left(\frac{4}{25}\right)^{2}\\&\simeq 2.56\,\mathrm{MeV},\\{e_3}_{\mathrm{new}}-{e_3}&\simeq 1.92\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[密度汎関数の三体項の圧力]
粒子当たり $e_3=C_3\rho^2$ の項が加える圧力を $P=\rho^2e_3'$ から求める。

量の定義に従い、共通の単位で次の値を用いる。
$$
C_{3}=100\,\mathrm{MeV\,fm^6},\qquad \rho=\frac{2}{25}\,\mathrm{fm}^{-3}
$$

この条件で求める量を ${P_3}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{P_3}&=2 \cdot C_{3} \cdot \rho^{3}\\&=2 \cdot \left(100\right) \cdot \left(\frac{2}{25}\right)^{3}\\&\simeq 0.1024\,\mathrm{MeV\,fm}^{-3}.\end{aligned}
$$

微分に含まれる密度依存が圧力の係数を決める。一定の一体エネルギーを加えてもこの効果は出ない。

前例と同じ密度の変化で圧力を比べる。そこで次の値に変更し、ほかの入力を固定する。
$$
\rho=\frac{4}{25}\,\mathrm{fm}^{-3}
$$

このときの値を ${P_3}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{P_3}_{\mathrm{new}}&=2 \cdot \left(100\right) \cdot \left(\frac{4}{25}\right)^{3}\\&\simeq 0.8192\,\mathrm{MeV\,fm}^{-3},\\{P_3}_{\mathrm{new}}-{P_3}&\simeq 0.7168\,\mathrm{MeV\,fm}^{-3}.\end{aligned}
$$
:::

## 演習問題

::::exercise{level=basic}
$\mu c^2=500\,\mathrm{MeV}$、$B=1\,\mathrm{MeV}$ の外部減衰長を $\hbar c=200\,\mathrm{MeV\,fm}$ として求めよ。
:::solution
$\kappa=\sqrt{2\mu c^2B}/(\hbar c)=\sqrt{1000}/200\,\mathrm{fm}^{-1}$ なので減衰長は $2\sqrt{10}\,\mathrm{fm}$ である。
:::
::::

::::exercise{level=basic}
$a=10\,\mathrm{fm}$ の区別可能な散乱体のゼロ運動量断面積を b で求めよ。
:::solution
$f_0(0)=-a$ なので $4\pi a^2=400\pi\,\mathrm{fm^2}=4\pi\,\mathrm b$ である。
:::
::::

::::exercise{level=basic}
二核子のスピン交換子の固有値を $S=0,1$ について求めよ。
:::solution
$S=0,1$ における $\vb\sigma_1\cdot\vb\sigma_2$ は $-3,1$ である。$(1+\vb\sigma_1\cdot\vb\sigma_2)/2$ はそれぞれ $-1,+1$ を与える。
:::
::::

::::exercise{level=basic}
核子が 5 個ある場合の対と三つ組の個数を求めよ。
:::solution
添字の順序を区別しないので $\binom52=\binom53=10$ である。
:::
::::

::::exercise{level=advanced}
深さ $V_0>0$、半径 $R$ の S 波井戸の散乱長を求めよ。ゼロエネルギーの内側波数を $K=\sqrt{2\mu V_0}/\hbar$ とする。
:::solution
内側の正則解は $C\sin Kr$、外側は $D(r-a)$ である。境界で値と微分をつなぐと $\tan(KR)/K=R-a$ なので $a=R-\tan(KR)/K$ が得られる。$\cos(KR)=0$ では外部解の微分がゼロとなり、有限な $a$ で規格化できない。
:::
::::

::::exercise{level=advanced}
固有値から、三重項における $S_{12}^2=8-2S_{12}$ を示せ。
:::solution
$\hat{\vb r}$ を量子化軸にすると、スピン積は $m_S=\pm1$ で $+1$、$m_S=0$ で $-1$ である。三重項の $\vb\sigma_1\cdot\vb\sigma_2$ は $+1$ なので固有値は $2,-4,2$ となる。両固有値で $x^2+2x-8=0$ が成り立つので、演算子の二次式もゼロになる。
:::
::::

::::exercise{level=advanced}
$Y=e^{-\lambda r}/(4\pi r)$ の $Y''-Y'/r$ を計算せよ。
:::solution
$Y'=-(\lambda+1/r)Y$、$Y''=[1/r^2+(\lambda+1/r)^2]Y$ である。さらに $-Y'/r=(\lambda/r+1/r^2)Y$ を足すと $(\lambda^2+3\lambda/r+3/r^2)Y$ を得る。
:::
::::

::::exercise{level=advanced}
$\delta_0=-kR$ で :eqref[eq-nuc-effective-range] の係数を求めよ。
:::solution
$\cos x=1-x^2/2+\cdots$、$\sin x=x(1-x^2/6+\cdots)$ なので $k\cot(-kR)=-1/R+k^2R/3+O(k^4)$ である。従って $a=R,r_e=2R/3$ となる。
:::
::::

## 確認テスト

::::quiz
重陽子のテンソル混合で結ばれる二つの軌道を選べ。

- [ ] S 波と P 波
- [x] S 波と D 波
- [ ] P 波と D 波

:::explanation
正パリティなので $L$ は偶数である。$S=1$ と $J=1$ を両立する偶数 $L$ は 0 と 2 である。
:::
::::

::::quiz{answer="2" tol="1e-6"}
剛体球の半径は $3\,\mathrm{fm}$ である。有効距離を fm で求めよ。
:::explanation
$k\cot\delta_0=-1/R+(R/3)k^2+\cdots$ の係数を比較して $r_e=2R/3=2\,\mathrm{fm}$ である。
:::
::::

::::quiz
荷電独立な力について成り立つものを選べ。

- [x] 同じ $T$ の各 $T_z$ で強い相互作用の行列要素が同じである。
- [ ] $T=0$ と $T=1$ の力が必ず等しい。
- [ ] Coulomb 力もアイソスピンで回転できる。

:::explanation
荷電独立性は $\SU(2)$ 回転不変性である。多重項の内部では $T_z$ を選ばないが、異なる $T$ の係数は独立である。Coulomb 力は電荷の違いを区別する。
:::
::::

::::quiz{answer="-3" tol="1e-6"}
$T=0$ の $\vb\tau_1\cdot\vb\tau_2$ の固有値を求めよ。
:::explanation
二核子の $\vb T^2$ を展開すると固有値は $2T(T+1)-3$ であり、$T=0$ で $-3$ である。
:::
::::

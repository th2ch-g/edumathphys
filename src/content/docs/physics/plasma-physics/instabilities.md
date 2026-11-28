---
title: MHD 不安定性と磁気再結合
description: 理想 MHD の不安定性と、抵抗が磁力線の連結を変える過程を解析する。
notation:
  - symbol: '\gamma'
    meaning: 時間成長率
  - symbol: '\Delta^{\prime}'
    meaning: テアリング外部解の対数微分の跳び
  - symbol: 'S_L'
    meaning: Lundquist 数
---

## この章の目標

- 圧力・重力による不安定性と磁気張力を比較する。
- テアリングの外部解と内部層の役割を説明する。
- 磁気再結合の速度を収支から求める。

## 交換型不安定性と Rayleigh–Taylor 不安定性

流体要素を磁場に沿う変化なしに入れ替えると、張力の代償を小さくできる。この型の摂動が平衡からエネルギーを受けて成長する現象を :term[交換型不安定性]{en="interchange instability" yomi="こうかんがたふあんていせい"} と呼ぶ。平面境界 $z=0$ に、上側 $\rho_+$、下側 $\rho_-$ の非圧縮な静止流体を置く。重力は $-g\hat{\vb z}$、磁場 $\vb B_\pm$ は境界に平行である。境界変位を $\xi e^{i\vb k\cdot\vb x+\gamma t}$ とし、上下で摂動が遠方に減衰するものとする。

誘導式は $\vb b_1=i(\vb k\cdot\vb B_\pm)\vb\xi$、運動式は $[\rho_\pm\gamma^2+(\vb k\cdot\vb B_\pm)^2/\mu_0]\vb\xi=-\nabla p_{\rm tot,1}$ である。発散を取れば全圧摂動は調和関数であり、上下で $e^{\mp kz}$ となる。各層で変位の係数が非零なら、境界で $\xi_z=\xi$ なので $p_{\rm tot,1}^\pm=\pm[\rho_\pm\gamma^2+(\vb k\cdot\vb B_\pm)^2/\mu_0]\xi/k$ である。全圧の連続を変位した境界で課し、平衡圧力勾配 $-\rho_\pm g$ の差を加えると
$$
\gamma^2=\frac{gk(\rho_+-\rho_-)-[(\vb k\cdot\vb B_+)^2+(\vb k\cdot\vb B_-)^2]/\mu_0}{\rho_++\rho_-}.
\label{eq-plasma-rt}
$$
係数が消える点では、界面分枝の中立条件をこの式の連続な極限で読む。これは :term[Rayleigh–Taylor 不安定性]{en="Rayleigh–Taylor instability" yomi="れいりーていらーふあんていせい"} の磁気張力つきの式である。

:::example[向きの選択]
上下の磁場がともに $x$ 方向なら $\vb k=k\hat{\vb y}$ が張力項をゼロにする。$\rho_+=3\rho_-$ でその向きを選ぶと $\gamma^2=gk/2$ である。
:::

:::example[短波長側の安定化]
上下で同じ $\vb B=B\hat{\vb x}$ を持ち、重い流体が上にあるとする。$\vb k=k\hat{\vb x}$ なら $\gamma^2=A_g gk-v_T^2k^2$、$A_g=(\rho_+-\rho_-)/(\rho_++\rho_-)$、$v_T^2=2B^2/(\mu_0(\rho_++\rho_-))$ である。最大成長は $k=A_g g/(2v_T^2)$、その二倍の波数が限界である。$\vb k$ を磁場に垂直にした場合、同じ表面模型にはこの短波長側の抑制がない。
:::

### せん断流のある界面

上下に一様な接線流 $\vb U_\pm$ を加え、時間依存を $e^{-i\omega t}$ に変える。各層で線形化する際、時間微分は物質微分 $-i(\omega-\vb k\cdot\vb U_\pm)$ に置き換わる。したがって全圧の接続条件では :eqref[eq-plasma-rt] の $\rho_\pm\gamma^2$ が $-\rho_\pm(\omega-\vb k\cdot\vb U_\pm)^2$ になる。$\rho_\Sigma=\rho_++\rho_-$、$\overline{\vb U}=(\rho_+\vb U_++\rho_-\vb U_-)/\rho_\Sigma$ とおいて平方を完成すると
$$
\begin{aligned}
(\omega-\vb k\cdot\overline{\vb U})^2
&=\frac{(\vb k\cdot\vb B_+)^2+(\vb k\cdot\vb B_-)^2}{\mu_0\rho_\Sigma}
-\frac{gk(\rho_+-\rho_-)}{\rho_\Sigma}\\
&\quad-\frac{\rho_+\rho_-}{\rho_\Sigma^2}[\vb k\cdot(\vb U_+-\vb U_-)]^2 .
\end{aligned}
$$
最後の項が :term[Kelvin–Helmholtz 不安定性]{en="Kelvin–Helmholtz instability" yomi="けるびんへるむほるつふあんていせい"} の駆動である。せん断流も磁場と同じく、波数方向へ射影した成分だけが効く。

:::example[磁化した速度不連続面]
密度と磁場が上下で等しく、$\vb U_\pm=\pm U\hat{\vb x}$、$\vb B_\pm=B\hat{\vb x}$、$\vb k=k\hat{\vb x}$ とする。右辺は $k^2(v_A^2-U^2)$ であり、$U>v_A=B/\sqrt{\mu_0\rho}$ のとき成長する。$U=v_A$ では根が重なるので、指数成長率がゼロという判定と初期値解の有界性を区別する。
:::

### 磁力線の曲率が駆動する交換

低 $\beta$ の真空磁場に近い平衡で、曲率ベクトルを $\vb\kappa=(\vb b\cdot\nabla)\vb b$ とする。このとき $\nabla_\perp B/B=\vb\kappa$ なので、勾配 $B$ ドリフトと曲率ドリフトの和は $(m_s/q_s)(v_\parallel^2+v_\perp^2/2)\vb b\times\vb\kappa/B$ である。これは力によるドリフト $\vb F\times\vb B/(q_sB^2)$ で $\vb F=-m_s(v_\parallel^2+v_\perp^2/2)\vb\kappa$ と置いた形である。等方な分布では $m_s n_s\langle v_\parallel^2\rangle=p_s$ と $m_s n_s\langle v_\perp^2/2\rangle=p_s$ が成り立つ。従って電荷を重みとしたドリフトには、合計で有効な体積力 $-2p\vb\kappa$ に相当する寄与がある。

重い層が有効重力と反対側にあるときに Rayleigh–Taylor の変位が駆動されたのと同じく、$-2p\vb\kappa$ が高圧側から低圧側へ向くと交換を駆動する。悪い曲率では従って $\vb\kappa\cdot\nabla p>0$ である。張力を落とし、圧力の変化長を $L_p=p/|\nabla_\perp p|$、曲率半径を $R_c=1/|\vb\kappa|$ とすると、両者が平行な領域の局所成長率の尺度は $\gamma\sim\sqrt{2p/(\rho L_p R_c)}$ である。これは有効加速度 $2p/(\rho R_c)$ を勾配長 $L_p$ で割って得た尺度であり、有限長の磁力線の固有値を与えない。

例えば軸対称な鏡では、磁束 $\Delta\Phi$ の細い管の面積は $\Delta\Phi/B$ である。半径 $r$ での方位ドリフト角速度はドリフト速度を $r$ で割った値なので、電荷の移動を管全体で平均する際には曲率の符号を保った重み $[p_\parallel+p_\perp]/(r B^2)$ が入る。$B^{-2}$ の一因子はドリフト速度、もう一因子は管断面に由来する。同じ曲率半径で長さと圧力が等しい二領域でも、片側で磁場が二倍ならその領域の重みは四分の一になる。鏡の端と中央で曲率の符号が異なる場合には、符号だけでなくこの重みを考慮する。 :cite[goldston-rutherford]

### 界面波の不安定波数と最速波長

非圧縮な二層の密度を $\rho_+>\rho_-$ とし、重力 $-g\hat{\vb z}$ に逆らって重い方を上に置く。界面変位を $\xi$ とすると各半空間の速度ポテンシャルは $e^{-k|z|}$ で減衰し、運動エネルギーの慣性係数は $(\rho_++\rho_-)/k$ になる。重力による位置エネルギーの二次項は $-(\rho_+-\rho_-)g\xi^2/2$ である。両側の水平磁場による曲げのエネルギーは $(B_{+,\parallel}^2+B_{-,\parallel}^2)k\xi^2/(2\mu_0)$ である。振幅の成長率は $\gamma^2=ak-bk^2$、$a=g(\rho_+-\rho_-)/(\rho_++\rho_-)$、$b=(B_{+,\parallel}^2+B_{-,\parallel}^2)/[\mu_0(\rho_++\rho_-)]$ となる。磁場と波数が直交して $b=0$ ならこの無厚さ模型には最速波数が存在しない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $a$ | 重力の有効加速度 | $\mathrm{m\,s^{-2}}$ |
| $b$ | 曲げを支える速度係数 | $\mathrm{m^2\,s^{-2}}$ |
| $k$ | 水平波数 | $\mathrm{m^{-1}}$ |

重力駆動を $g_k$ とする。変形した層の慣性を含めた値である。

$$
g_k=a k
$$

曲げの復元を $t_k$ とする。波数方向の磁場成分だけが寄与する。

$$
t_k=b k^{2}
$$

成長率の二乗を $\gamma^2$ とする。負の場合は振動する。

$$
\gamma^2=g_{k} - t_{k}
$$

不安定帯の上限を $k_c$ とする。二つの寄与が同じになる。

$$
k_c=\frac{a}{b}
$$

最速成長の波数を $k_*$ とする。$\gamma^2$ を微分する。

$$
k_*=\frac{k_{c}}{2}
$$

最大成長率を $\gamma_*$ とする。最速波数を代入する。

$$
\gamma_*=\frac{a}{2 \sqrt{b}}
$$

:::example[最長側の不安定波]
重力と張力を比べる。

入力値を次のように置く。

$$
\begin{aligned}
a&=4\,\mathrm{m\,s^{-2}} \\b&=1\,\mathrm{m^2\,s^{-2}} \\k&=\frac{1}{10}\,\mathrm{m^{-1}}
\end{aligned}
$$

これらを順に代入する。

重力駆動は

$$
g_k=\left(4\right) \left(\frac{1}{10}\right)=\frac{2}{5}\,\mathrm{s^{-2}}
$$

である。

曲げの復元は

$$
t_k=\left(1\right) \left(\frac{1}{10}\right)^{2}=\frac{1}{100}\,\mathrm{s^{-2}}
$$

である。

成長率の二乗は

$$
\gamma^2=- \left(\frac{1}{100}\right) + \left(\frac{2}{5}\right)=\frac{39}{100}\,\mathrm{s^{-2}}
$$

である。

不安定帯の上限は

$$
k_c=\frac{\left(4\right)}{\left(1\right)}=4\,\mathrm{m^{-1}}
$$

である。

最速成長の波数は

$$
k_*=\frac{\left(4\right)}{2}=2\,\mathrm{m^{-1}}
$$

である。

最大成長率は

$$
\gamma_*=\frac{\left(4\right)}{2 \sqrt{\left(1\right)}}=2\,\mathrm{s^{-1}}
$$

である。

重力が優勢でも波数が小さいため成長は最大ではない。
:::

:::example[最速波数]
成長率が最大値に一致するか調べる。

入力値を次のように置く。

$$
\begin{aligned}
a&=4\,\mathrm{m\,s^{-2}} \\b&=1\,\mathrm{m^2\,s^{-2}} \\k&=2\,\mathrm{m^{-1}}
\end{aligned}
$$

これらを順に代入する。

重力駆動は

$$
g_k=\left(2\right) \left(4\right)=8\,\mathrm{s^{-2}}
$$

である。

曲げの復元は

$$
t_k=\left(1\right) \left(2\right)^{2}=4\,\mathrm{s^{-2}}
$$

である。

成長率の二乗は

$$
\gamma^2=- \left(4\right) + \left(8\right)=4\,\mathrm{s^{-2}}
$$

である。

不安定帯の上限は

$$
k_c=\frac{\left(4\right)}{\left(1\right)}=4\,\mathrm{m^{-1}}
$$

である。

最速成長の波数は

$$
k_*=\frac{\left(4\right)}{2}=2\,\mathrm{m^{-1}}
$$

である。

最大成長率は

$$
\gamma_*=\frac{\left(4\right)}{2 \sqrt{\left(1\right)}}=2\,\mathrm{s^{-1}}
$$

である。

波長をさらに短くすると張力が急増する。
:::

:::example[中立波数]
波数で打ち消しが起こる条件を確かめる。

入力値を次のように置く。

$$
\begin{aligned}
a&=4\,\mathrm{m\,s^{-2}} \\b&=1\,\mathrm{m^2\,s^{-2}} \\k&=4\,\mathrm{m^{-1}}
\end{aligned}
$$

これらを順に代入する。

重力駆動は

$$
g_k=\left(4\right)^{2}=16\,\mathrm{s^{-2}}
$$

である。

曲げの復元は

$$
t_k=\left(1\right) \left(4\right)^{2}=16\,\mathrm{s^{-2}}
$$

である。

成長率の二乗は

$$
\gamma^2=0\,\mathrm{s^{-2}}
$$

である。

不安定帯の上限は

$$
k_c=\frac{\left(4\right)}{\left(1\right)}=4\,\mathrm{m^{-1}}
$$

である。

最速成長の波数は

$$
k_*=\frac{\left(4\right)}{2}=2\,\mathrm{m^{-1}}
$$

である。

最大成長率は

$$
\gamma_*=\frac{\left(4\right)}{2 \sqrt{\left(1\right)}}=2\,\mathrm{s^{-1}}
$$

である。

復元と駆動が等しくなる。
:::

:::example[安定側]
高波数での符号を調べる。

入力値を次のように置く。

$$
\begin{aligned}
a&=4\,\mathrm{m\,s^{-2}} \\b&=1\,\mathrm{m^2\,s^{-2}} \\k&=5\,\mathrm{m^{-1}}
\end{aligned}
$$

これらを順に代入する。

重力駆動は

$$
g_k=\left(4\right) \left(5\right)=20\,\mathrm{s^{-2}}
$$

である。

曲げの復元は

$$
t_k=\left(1\right) \left(5\right)^{2}=25\,\mathrm{s^{-2}}
$$

である。

成長率の二乗は

$$
\gamma^2=\left(20\right) - \left(25\right)=-5\,\mathrm{s^{-2}}
$$

である。

不安定帯の上限は

$$
k_c=\frac{\left(4\right)}{\left(1\right)}=4\,\mathrm{m^{-1}}
$$

である。

最速成長の波数は

$$
k_*=\frac{\left(4\right)}{2}=2\,\mathrm{m^{-1}}
$$

である。

最大成長率は

$$
\gamma_*=\frac{\left(4\right)}{2 \sqrt{\left(1\right)}}=2\,\mathrm{s^{-1}}
$$

である。

$\gamma^2<0$ なので指数成長はない。
:::


## ソーセージ不安定性とキンク不安定性

電流を軸方向に流す円柱で摂動を $e^{i(m\theta+k_z z)}$ と分類する。$m=0$ は断面が周期的に細くなる :term[ソーセージモード]{en="sausage mode" yomi="そーせーじもーど"}、$m=1$ は軸が横へずれる :term[キンクモード]{en="kink mode" yomi="きんくもーど"} である。個々の成長率には内部電流と外部領域の境界条件が要る。

純 Z ピンチの外部真空では、軸電流が一定の軸対称なゆっくりした変形に対し表面磁場は $B_\theta(a)=\mu_0 I/(2\pi a)$ である。従って外部磁気圧の変化は $\delta p_B=-2p_B\delta a/a$ である。細くなる箇所に磁気圧が増えることが駆動の一因となる。内部熱圧の増加や軸磁束の圧縮を同時に入れなければ、これだけで符号は決まらない。例えば一単位長さの物質を独立に断熱圧縮すると $n\propto a^{-2}$、$p\propto a^{-2\gamma_{\rm ad}}$ なので、熱圧も変化する。軸方向の流れを許す実際の摂動では各断面の質量は固定されない。

:::example[偏波と軸の移動]
小さな一様横変位 $\xi_0\hat{\vb x}$ は円柱表面に動径成分 $\xi_0\cos\theta$ を持つ。このため軸の変位は $m=\pm1$ の実結合になる。
:::

:::example[軸磁場の圧縮]
断面の軸磁束 $B_z\pi a^2$ が一定なら、$\delta B_z/B_z=-2\delta a/a$ である。従って軸磁場の圧力変化は $\delta p_{B_z}=-4p_{B_z}\delta a/a$ である。細くなった部分でこの内部圧が増すため、固定した軸電流による外部圧だけを考えた場合より収縮しにくい。
:::

### 表面電流をもつ円柱のソーセージ分枝

半径 $a$、密度 $\rho$ の非圧縮なプラズマを真空が囲むとする。内部磁場は一様な軸磁場 $B_i\hat{\vb z}$、外部磁場は $B_a(a/r)\hat{\vb\theta}$ とし、界面の表面電流が磁場の跳びを担う。平衡圧力は $p_i+(B_i^2/2\mu_0)=B_a^2/(2\mu_0)$ である。$k\ne0$ の軸対称変位 $\xi_r(a)e^{ikz-i\omega t}$ を考える。

外部摂動は発散と回転がゼロであり、軸対称かつ $k\ne0$ なので磁気スカラーポテンシャルで表せる。界面の接線磁場は $\theta$ 方向であり、変位が $\theta$ によらないため、外部摂動の法線成分はゼロである。無限遠で減衰する外部摂動は従ってゼロである。移動した界面での外部全圧変化には、平衡の磁気圧勾配 $-B_a^2/(\mu_0 a)$ に $\xi_r(a)$ を掛けた項が残る。

内部の全圧摂動 $P(r)$ は、非圧縮条件から $P''+P'/r-k^2P=0$ を満たす。軸で正則な解は $P=C I_0(|k|r)$ である。運動式の半径成分は $\partial_r P=(\rho\omega^2-k^2 B_i^2/\mu_0)\xi_r$ なので、界面で $P/\xi_r=(\rho\omega^2-k^2 B_i^2/\mu_0)I_0(|k|a)/[|k|I_1(|k|a)]$ である。これを外部の圧力変化に等置すると
$$
\omega^2=\frac{k^2 B_i^2}{\mu_0\rho}
-\frac{B_a^2}{\mu_0\rho a}|k|\frac{I_1(|k|a)}{I_0(|k|a)} .
$$
圧力勾配が駆動し、内部の軸磁場が張力を与えることが分かる。

:::example[長波長での軸磁場の効果]
$|k|a\ll1$ なら $I_1(|k|a)/I_0(|k|a)\simeq |k|a/2$ である。$B_i=0$ では $\gamma\simeq |k|B_a/\sqrt{2\mu_0\rho}$ となる。$B_i^2>B_a^2/2$ であればこの長波長分枝は振動する。この境界は表面電流の平衡に対する値であり、内部に電流が分布する平衡では内部の圧力方程式も変わる。
:::

### 円柱のモード形を面積と軸変位から区別する

非負の整数 $m$ について円柱断面を $r=a[1+\epsilon\cos(m\theta)]$ とし、ある $z$ 断面で位相を固定する。面積は $\tfrac12\int r^2\dd\theta$ であり、$m\ne0$ では一次変化が消える。$x$ 方向の一次モーメントは $\tfrac13\int r^3\cos\theta\,\dd\theta$ である。角の直交性より $m=1$ だけが一次で重心を移動させる。$m=0$ は一次の面積変化を起こす。これらはソーセージとキンクを識別する幾何学で、成長率には軸方向波数、磁気力、境界条件を含める必要がある。$m=1$ の極座標の形は一次で並進と一致するが、二次までの剛体並進とは異なる。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $m$ | 方位モード番号。非負の整数とする。 | 1 |
| $\epsilon$ | 半径の相対変位 | 1 |
| $\theta$ | 評価する方位角 | 1 |

境界半径の比を $r/a$ とする。実振幅を一つの断面で評価する。

$$
r/a=\epsilon \cos{\left(\theta m \right)} + 1
$$

断面積の比を $A/(\pi a^2)$ とする。余弦と余弦二乗を一周積分する。

$$
A/(\pi a^2)=\frac{\epsilon^{2} \left(\delta_{0 m} + 1\right)}{2} + 2 \epsilon \delta_{0 m} + 1
$$

重心の一次変位を $x_c^{(1)}/a$ とする。一次モーメントの角積分を使う。

$$
x_c^{(1)}/a=\epsilon \delta_{1 m}
$$

面積の一次変化を $\delta A^{(1)}/(\pi a^2)$ とする。線形の圧縮を測る。

$$
\delta A^{(1)}/(\pi a^2)=2 \epsilon \delta_{0 m}
$$

この境界形の二次面積変化を $\delta A^{(2)}/(\pi a^2)$ とする。二次の境界変位は与えていない。

$$
\delta A^{(2)}/(\pi a^2)=A/(\pi a^2) - \delta A^{(1)}/(\pi a^2) - 1
$$

局所変位の二乗を $(r/a-1)^2$ とする。一点で同じ変位が見えても全体形は異なる。

$$
(r/a-1)^2=\left(r/a - 1\right)^{2}
$$

:::example[ソーセージの山]
膨らむ断面の面積を求める。

入力値を次のように置く。

$$
\begin{aligned}
m&=0 \\\epsilon&=\frac{1}{10} \\\theta&=0
\end{aligned}
$$

これらを順に代入する。

境界半径の比は

$$
r/a=\left(\frac{1}{10}\right) \cos{\left(\left(0\right)^{2} \right)} + 1=\frac{11}{10}
$$

である。

断面積の比は

$$
A/(\pi a^2)=\frac{\left(\frac{1}{10}\right)^{2} \left(\delta_{0 \left(0\right)} + 1\right)}{2} + 2 \left(\frac{1}{10}\right) \delta_{0 \left(0\right)} + 1=\frac{121}{100}
$$

である。

重心の一次変位は

$$
x_c^{(1)}/a=\left(\frac{1}{10}\right) \delta_{1 \left(0\right)}=0
$$

である。

面積の一次変化は

$$
\delta A^{(1)}/(\pi a^2)=2 \left(\frac{1}{10}\right) \delta_{0 \left(0\right)}=\frac{1}{5}
$$

である。

この境界形の二次面積変化は

$$
\delta A^{(2)}/(\pi a^2)=\left(\frac{121}{100}\right) - \left(\frac{1}{5}\right) - 1=\frac{1}{100}
$$

である。

局所変位の二乗は

$$
(r/a-1)^2=\left(\left(\frac{11}{10}\right) - 1\right)^{2}=\frac{1}{100}
$$

である。

一次の面積変化があり、重心は変わらない。
:::

:::example[キンクの山]
同じ境界変位で面積と重心を比べる。

入力値を次のように置く。

$$
\begin{aligned}
m&=1 \\\epsilon&=\frac{1}{10} \\\theta&=0
\end{aligned}
$$

これらを順に代入する。

境界半径の比は

$$
r/a=\left(\frac{1}{10}\right) \cos{\left(\left(0\right) \left(1\right) \right)} + 1=\frac{11}{10}
$$

である。

断面積の比は

$$
A/(\pi a^2)=\frac{\left(\frac{1}{10}\right)^{2} \left(\delta_{0 \left(1\right)} + 1\right)}{2} + 2 \left(\frac{1}{10}\right) \delta_{0 \left(1\right)} + 1=\frac{201}{200}
$$

である。

重心の一次変位は

$$
x_c^{(1)}/a=\left(\frac{1}{10}\right) \delta_{1 \left(1\right)}=\frac{1}{10}
$$

である。

面積の一次変化は

$$
\delta A^{(1)}/(\pi a^2)=2 \left(\frac{1}{10}\right) \delta_{0 \left(1\right)}=0
$$

である。

この境界形の二次面積変化は

$$
\delta A^{(2)}/(\pi a^2)=- \left(0\right) + \left(\frac{201}{200}\right) - 1=\frac{1}{200}
$$

である。

局所変位の二乗は

$$
(r/a-1)^2=\left(\left(\frac{11}{10}\right) - 1\right)^{2}=\frac{1}{100}
$$

である。

面積の一次変化はゼロだが重心が移る。
:::

:::example[楕円形に近い変形]
方位二のモードをキンクと区別する。

入力値を次のように置く。

$$
\begin{aligned}
m&=2 \\\epsilon&=\frac{1}{10} \\\theta&=0
\end{aligned}
$$

これらを順に代入する。

境界半径の比は

$$
r/a=\left(\frac{1}{10}\right) \cos{\left(\left(0\right) \left(2\right) \right)} + 1=\frac{11}{10}
$$

である。

断面積の比は

$$
A/(\pi a^2)=\frac{\left(\frac{1}{10}\right)^{2} \left(\delta_{0 \left(2\right)} + 1\right)}{2} + 2 \left(\frac{1}{10}\right) \delta_{0 \left(2\right)} + 1=\frac{201}{200}
$$

である。

重心の一次変位は

$$
x_c^{(1)}/a=\left(\frac{1}{10}\right) \delta_{1 \left(2\right)}=0
$$

である。

面積の一次変化は

$$
\delta A^{(1)}/(\pi a^2)=2 \left(\frac{1}{10}\right) \delta_{0 \left(2\right)}=0
$$

である。

この境界形の二次面積変化は

$$
\delta A^{(2)}/(\pi a^2)=- \left(0\right) + \left(\frac{201}{200}\right) - 1=\frac{1}{200}
$$

である。

局所変位の二乗は

$$
(r/a-1)^2=\left(\left(\frac{11}{10}\right) - 1\right)^{2}=\frac{1}{100}
$$

である。

一次では面積と重心をどちらも変えない。
:::

:::example[キンクの節]
一点の変位がゼロでも重心が動くか確かめる。

入力値を次のように置く。

$$
\begin{aligned}
m&=1 \\\epsilon&=\frac{1}{10} \\\theta&=\frac{\pi}{2}
\end{aligned}
$$

これらを順に代入する。

境界半径の比は

$$
r/a=\left(\frac{1}{10}\right) \cos{\left(\left(1\right) \left(\frac{\pi}{2}\right) \right)} + 1=1
$$

である。

断面積の比は

$$
A/(\pi a^2)=\frac{\left(\frac{1}{10}\right)^{2} \left(\delta_{0 \left(1\right)} + 1\right)}{2} + 2 \left(\frac{1}{10}\right) \delta_{0 \left(1\right)} + 1=\frac{201}{200}
$$

である。

重心の一次変位は

$$
x_c^{(1)}/a=\left(\frac{1}{10}\right) \delta_{1 \left(1\right)}=\frac{1}{10}
$$

である。

面積の一次変化は

$$
\delta A^{(1)}/(\pi a^2)=2 \left(\frac{1}{10}\right) \delta_{0 \left(1\right)}=0
$$

である。

この境界形の二次面積変化は

$$
\delta A^{(2)}/(\pi a^2)=- \left(0\right) + \left(\frac{201}{200}\right) - 1=\frac{1}{200}
$$

である。

局所変位の二乗は

$$
(r/a-1)^2=\left(\left(1\right) - 1\right)^{2}=0
$$

である。

方位測定を一点で行うとモードを見落とす。
:::


## Kruskal–Shafranov 条件

周期長 $L$ の円柱を大アスペクト比トーラスの近似とする。境界 $r=a$ で磁力線が軸方向の一周期に回る角は $L B_\theta/(a B_z)$ である。$q(a)=2\pi a B_z/(L B_\theta)$ は式 :eqref[eq-plasma-safety-factor] の円柱版である。$m=1$ のヘリカル変位に対する線に沿う波数は $k_\parallel=(B_z k_z+B_\theta/a)/B$ であり、$k_z=-2\pi/L$ を選ぶと $q(a)=1$ で消える。これが単純な自由境界電流柱の外部キンクの :term[Kruskal–Shafranov 条件]{en="Kruskal–Shafranov criterion" yomi="くるすかるしゃふらのふじょうけん"} に現れるねじれである。標準の長波長円柱モデルでは $q(a)<1$ 側が不安定になる。しきい値の厳密な変分導出は外部真空磁場と表面項を含める必要があり、固定境界の式 :eqref[eq-plasma-energy-principle] の適用範囲を越えるため :cite[bellan-plasma] を参照する。近接導体壁や線結びされた端を持つ系では条件が変わる。

:::example[周期長の影響]
$L=2\pi a$、$B_z=2B_\theta(a)$ なら $q(a)=2$ で、磁力線の回転は $\pi$ である。$L$ だけ三倍にすると $q(a)=2/3$ になる。
:::

:::example[縦波数を変えたとき]
同じ円柱で $m=1,\ k_z=-4\pi/L$ とすれば、辺の $k_\parallel=0$ は $q(a)=1/2$ で起きる。$q(a)=1$ の候補は最長波長の $k_z=-2\pi/L$ から来る。$k_\parallel=0$ だけで各モードの成長率や境界壁の影響を計算することはできない。
:::

### 自由境界の試験変位

一つの円柱模型で閾値の力学を確かめる。内部を $B_\theta=Ar$、$B_z=$ 一定、$\rho=$ 一定とし、外部真空を $B_\theta=Aa^2/r$ と同じ軸磁場でつなぐ。$A>0$、$\abs k a\ll1$ とする。内部で複素試験変位 $\vb\xi=\xi_0(\hat{\vb x}+i\hat{\vb y})e^{ikz}$ を採ると、動径成分は $\xi_0 e^{i\theta+ikz}$ であり、$\nabla\cdot\vb\xi=0$ である。磁場摂動は $\vb Q=(\vb B\cdot\nabla)\vb\xi-(\vb\xi\cdot\nabla)\vb B=i(kB_z+A)\vb\xi$ となる。実際、$(\vb\xi\cdot\nabla)\vb B=A\hat{\vb z}\times\vb\xi=-iA\vb\xi$ である。

式 :eqref[eq-plasma-energy-principle] の内部積分を複素変位の実二次形式として評価する。電流は $\vb j=2A\hat{\vb z}/\mu_0$、$\vb j\cdot(\vb\xi^*\times\vb Q)=-2A(kB_z+A)|\vb\xi|^2/\mu_0$ であり、$Q^2$ の寄与を足すと内部の係数は $[(kB_z)^2-A^2]/\mu_0$ である。界面の全圧勾配は内外で $-A^2 a/\mu_0$ に一致する。自由境界の二次変分に現れる全圧勾配の跳びの項は、この平衡では消える。

外部真空の磁場摂動を $\nabla\chi$ と置く。界面での動径成分は $i(kB_z+A)\xi_0 e^{i\theta+ikz}$ である。$\abs k a\ll1$ の極限で、半径に関する調和関数の外部解は $\chi(a)(a/r)e^{i\theta+ikz}$ であり、下付きの振幅は $\chi(a)=-ai(kB_z+A)\xi_0$ である。有限 $k$ の外部減衰解は変形 Bessel 関数 $K_1(\abs k r)$ であり、これは境界付近の長波長極限に相当する。真空エネルギーを一周期当たりで評価すると、積分 $\int_a^\infty 2(a/r)^4r\dd r=a^2$ から、内部の変位ノルムで規格化した係数は $(kB_z+A)^2/\mu_0$ である。

両方を足すと符号は $2kB_z(kB_z+A)$ で決まる。$k=-2\pi/L$、$B_z>0$ なら $0<-kB_z<A$ が負の符号を与える。$B_\theta(a)=Aa$ よりこの条件が $q(a)<1$ である。以上は無限遠に壁がある一様電流柱で、長波長の試験変位が負のエネルギーを持つことを示す。一般の電流分布や壁距離では、真空解と表面条件をその平衡で計算し直す。しきい値直近では落とした有限 $ka$ の補正も評価する。

:::example[二つのエネルギーの符号]
$kB_z=-A/2$ なら内部の係数は $-3A^2/(4\mu_0)$、真空の係数は $+A^2/(4\mu_0)$ であり、和は負である。境界の外に磁場を作るエネルギーを加えても、この長波長変位は駆動される。
:::

### 近接する導体壁の効果

同じ長波長の試験変位で、真空の外側に固定完全導体の円筒壁 $r=b>a$ を置く。$|k|b\ll1$ なら $m=1$ の外部ポテンシャルは $\chi=C_1r+C_2/r$ に比例し、壁の法線摂動ゼロ $\chi'(b)=0$ より $C_2=b^2 C_1$ である。真空磁場の二乗積分は Green の恒等式により内側境界での $-\overline\chi(a)\chi'(a)$ に比例する。従って同じ $\chi'(a)$ に対する真空のエネルギーは、壁なしの値の
$$
-\frac{\chi(a)}{a\chi'(a)}=\frac{b^2+a^2}{b^2-a^2}
$$
倍になる。右辺は壁が近づくほど大きい。例えば $b=\sqrt2 a$ なら三倍であり、先の $kB_z=-A/2$ の試験変位では内部の負の項をちょうど打ち消す。これは特定の変位での二次エネルギーの比較なので、他の変位に対する安定性には変分の最小化が要る。

### 円筒の巻き数と有限周期の許可モード

長さ $L$ で周期的な円柱に対し、磁力線のねじれ角は $\Phi=L B_\theta/(aB_z)$ である。$q_a=2\pi/\Phi$ とし、螺旋変位を $e^{i(m\theta-2\pi n z/L)}$ と書く。辺の平行微分係数は $(2\pi B_z/L)(m/q_a-n)$ である。これは円柱の記号上の共鳴であり、外部キンクの完全な安定性判定そのものではない。自由境界の長波長 $m=n=1$ の単純模型では、巻き角が一周を超える条件が Kruskal–Shafranov のしきい値である。導体壁が近い場合、外側真空場の変分が変わるので同じ巻き数から成長率を求めない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $q_a$ | 辺の安全係数。正を仮定する。 | 1 |
| $m$ | 方位モード番号 | 1 |
| $n$ | 軸周期の調波番号 | 1 |

一周期で巻く回数を $\Phi/(2\pi)$ とする。安全係数の逆数である。

$$
\Phi/(2\pi)=\frac{1}{q_{a}}
$$

ねじれ角を $\Phi$ とする。周回数をラジアンへ換算する。

$$
\Phi=2\pi/q_a
$$

平行微分係数の規格化値を $K_\parallel$ とする。軸と方位の変化が競合する。

$$
K_\parallel=\Phi/(2\pi) m - n
$$

平行張力の波数因子を $K_\parallel^2$ とする。共鳴面で小さくなる。

$$
K_\parallel^2=K_{\parallel}^{2}
$$

最長外部キンクの巻き超過を $\Phi/(2\pi)-1$ とする。自由境界かつ長波長のしきい値と比較する。

$$
\Phi/(2\pi)-1=\Phi/(2\pi) - 1
$$

このモードが共鳴する安全係数を $m/n$ とする。$n\ne0$ で $K_\parallel=0$ を解く。

$$
m/n=\frac{m}{n}
$$

:::example[一周より少ない巻き]
最長モードの巻きの尺度を計算する。

入力値を次のように置く。

$$
\begin{aligned}
q_a&=2 \\m&=1 \\n&=1
\end{aligned}
$$

これらを順に代入する。

一周期で巻く回数は

$$
\Phi/(2\pi)=\frac{1}{\left(2\right)}=\frac{1}{2}
$$

である。

ねじれ角は

$$
\Phi=2 \pi \left(\frac{1}{2}\right)\simeq 3.1416
$$

である。

平行微分係数の規格化値は

$$
K_\parallel=\left(1\right) \left(\frac{1}{2}\right) - \left(1\right)=- \frac{1}{2}
$$

である。

平行張力の波数因子は

$$
K_\parallel^2=\left(- \frac{1}{2}\right)^{2}=\frac{1}{4}
$$

である。

最長外部キンクの巻き超過は

$$
\Phi/(2\pi)-1=\left(\frac{1}{2}\right) - 1=- \frac{1}{2}
$$

である。

このモードが共鳴する安全係数は

$$
m/n=1
$$

である。

この自由境界の指標ではしきい値を超えない。
:::

:::example[巻き角が一周]
中立条件と平行微分係数が一致するか調べる。

入力値を次のように置く。

$$
\begin{aligned}
q_a&=1 \\m&=1 \\n&=1
\end{aligned}
$$

これらを順に代入する。

一周期で巻く回数は

$$
\Phi/(2\pi)=\frac{1}{\left(1\right)}=1
$$

である。

ねじれ角は

$$
\Phi=2 \pi \left(1\right)\simeq 6.2832
$$

である。

平行微分係数の規格化値は

$$
K_\parallel=\left(1\right)^{2} - \left(1\right)=0
$$

である。

平行張力の波数因子は

$$
K_\parallel^2=\left(0\right)^{2}=0
$$

である。

最長外部キンクの巻き超過は

$$
\Phi/(2\pi)-1=\left(1\right) - 1=0
$$

である。

このモードが共鳴する安全係数は

$$
m/n=1
$$

である。

この模型でちょうど境界に当たる。
:::

:::example[一周より多い巻き]
最長モードのしきい値超過を求める。

入力値を次のように置く。

$$
\begin{aligned}
q_a&=\frac{2}{3} \\m&=1 \\n&=1
\end{aligned}
$$

これらを順に代入する。

一周期で巻く回数は

$$
\Phi/(2\pi)=\frac{1}{\left(\frac{2}{3}\right)}=\frac{3}{2}
$$

である。

ねじれ角は

$$
\Phi=2 \pi \left(\frac{3}{2}\right)\simeq 9.4248
$$

である。

平行微分係数の規格化値は

$$
K_\parallel=\left(1\right) \left(\frac{3}{2}\right) - \left(1\right)=\frac{1}{2}
$$

である。

平行張力の波数因子は

$$
K_\parallel^2=\left(\frac{1}{2}\right)^{2}=\frac{1}{4}
$$

である。

最長外部キンクの巻き超過は

$$
\Phi/(2\pi)-1=\left(\frac{3}{2}\right) - 1=\frac{1}{2}
$$

である。

このモードが共鳴する安全係数は

$$
m/n=1
$$

である。

巻き数は一より大きく、単純模型の駆動側にある。
:::

:::example[別の有理数]
$m=3$ の共鳴を最長キンクから区別する。

入力値を次のように置く。

$$
\begin{aligned}
q_a&=\frac{3}{2} \\m&=3 \\n&=2
\end{aligned}
$$

これらを順に代入する。

一周期で巻く回数は

$$
\Phi/(2\pi)=\frac{1}{\left(\frac{3}{2}\right)}=\frac{2}{3}
$$

である。

ねじれ角は

$$
\Phi=2 \pi \left(\frac{2}{3}\right)\simeq 4.1888
$$

である。

平行微分係数の規格化値は

$$
K_\parallel=- \left(2\right) + \left(3\right) \left(\frac{2}{3}\right)=0
$$

である。

平行張力の波数因子は

$$
K_\parallel^2=\left(0\right)^{2}=0
$$

である。

最長外部キンクの巻き超過は

$$
\Phi/(2\pi)-1=\left(\frac{2}{3}\right) - 1=- \frac{1}{3}
$$

である。

このモードが共鳴する安全係数は

$$
m/n=\frac{\left(3\right)}{\left(2\right)}=\frac{3}{2}
$$

である。

平行微分は消えるが、$m=1$ の巻きしきい値を超えない。
:::


## バルーニングモード

トーラスでは外側と内側で磁場曲率と圧力勾配の相対方向が変わる。外側で両者が内向きなら、流体変位は局所的な圧力駆動を受ける。磁力線に沿って局在してこの駆動を利用するものを :term[バルーニングモード]{en="ballooning mode" yomi="ばるーにんぐもーど"} と呼ぶ。

### 局在と張力の競合

曲率駆動を $G$、線に沿う張力を $v_A^2$ とした模型の二次形式は $\int (v_A^2\abs{\xi'}^2-G(s)\abs\xi^2)\dd s$ である。$G=G_0>0$ の区間 $-\ell/2<s<\ell/2$ に $\xi=\cos(\pi s/\ell)$ を置き、外でゼロとする。$\int|\xi'|^2/\int|\xi|^2=(\pi/\ell)^2$ なので試験変位の符号は $v_A^2\pi^2/\ell^2-G_0$ で決まる。短い局在ほど張力の負担が増える。トーラスの高モード数極限には磁気シアで変わる垂直波数を組み込む必要があり、上の定数係数モデルに含まれない。幾何学を保持した固有値問題は :cite[goldston-rutherford] を参照する。

:::example[局在長]
$\ell=2\pi v_A/\sqrt{G_0}$ とすれば張力の寄与は $G_0/4$ であり、この試験変位で二次形式は負になる。
:::

:::example[異なる方向の曲率]
同じ線に沿って $G=G_0$ の部分と $G=-G_0$ の部分が等しい長さあるとする。一様な試験変位では駆動の二つの積分が相殺する。外側に局在させれば正の $G$ だけを利用できるが、境界の変位勾配による正の項を払う。
:::

### 磁気シアを含む円形トーラス近似

大アスペクト比の円形磁気面で、外側赤道を $\theta=0$ とする。長い磁力線に沿って局所解を探すときは、極角 $\theta$ を実軸上に延長する。高トロイダルモード数では、近くの二磁気面の巻き角の差により半径方向の波数が線に沿って変わる。トロイダル角を $\zeta$、小半径を $r$ とする。位相を $n[\zeta-q(r)(\theta-\theta_0)]$ とすると、径方向微分は $-nq'(\theta-\theta_0)$、局所の極角方向の波数は $-nq/r$ である。その比は $\hat s(\theta-\theta_0)$、$\hat s=(r/q)\dd q/\dd r$ となる。

円形面の中心の径方向変位も残す標準の $s$–$\alpha$ 近似では、主半径を $R_0$、代表トロイダル磁場を $B_0$ とした規格化圧力勾配 $\alpha=-2\mu_0R_0 q^2(\dd p/\dd r)/B_0^2$ を用い、この比に局所シアの補正 $-\alpha\sin\theta$ が加わる。$\Lambda=\hat s(\theta-\theta_0)-\alpha\sin\theta$ とおくと、垂直波数の二乗は $1+\Lambda^2$ に比例し、法線曲率と測地曲率を波数の向きに射影した駆動は $\cos\theta+\Lambda\sin\theta$ に比例する。ここでは円形・大アスペクト比のモデル係数を採用する。平衡の詳細な幾何からこれらの係数を求めるには [MHD 平衡と磁場閉じ込め](/physics/plasma-physics/equilibrium/) の平衡を用いる。

正の定数倍を省いた二次形式は
$$
\mathcal E[\xi]=\int_\R\Big[(1+\Lambda^2)|\xi'|^2
-\alpha(\cos\theta+\Lambda\sin\theta)|\xi|^2\Big]\dd\theta .
$$
局在する試験関数に対して部分積分すると、中立条件の Euler 方程式は $[(1+\Lambda^2)\xi']'+\alpha(\cos\theta+\Lambda\sin\theta)\xi=0$ である。端で減衰しない解を中立と判定せず、有限の切断域で調べる場合は切断位置と $\theta_0$ も変える。$\alpha=0$ なら二次形式は非負である。有限の $\alpha$ では曲率項と線の曲げがともに変わるので、圧力勾配だけを一項として加えた模型の境界からは判定できない。

### 局在試験関数の曲率仕事と張力

磁力線座標 $s$ に沿って $\delta W=\int[T|\xi'|^2-(D-C s^2)|\xi|^2]\dd s$ と規格化する。$T,C>0$、$D$ は悪い曲率中心での駆動である。$\xi=e^{-s^2/(2\ell^2)}$ なら $\int s^2\xi^2/\int\xi^2=\ell^2/2$、$\int|\xi'|^2/\int\xi^2=1/(2\ell^2)$ である。従って規格化変分は $T/(2\ell^2)+C\ell^2/2-D$ となる。幅を変えて最小にする条件は $\ell^4=T/C$ である。この放物型模型では試験関数から得た下限が厳密である。実際、無限遠で境界項が消える変位について
$$
\int_\R\left|\sqrt{T}\,\xi'+\sqrt{C}\,s\xi\right|^2\dd s
=\delta W+(D-\sqrt{TC})\int_\R|\xi|^2\dd s
$$
となる。交差項は $\sqrt{TC}\int s(\abs\xi^2)'\dd s=-\sqrt{TC}\int\abs\xi^2\dd s$ である。左辺がゼロになる微分方程式の解が先の Gauss 関数である。したがって $D>\sqrt{TC}$ がこの模型の不安定条件である。完全なトロイダルの境界値問題では曲率が周期的に反転する。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $T$ | 張力係数 | 1 |
| $C$ | 局在させる曲率の二次係数 | 1 |
| $D$ | 中心の圧力駆動 | 1 |
| $\ell$ | 試験関数の幅 | 1 |

張力項を $W_T$ とする。狭い関数では微分が大きい。

$$
W_T=\frac{T}{2 \ell^{2}}
$$

外側の曲率による安定化を $W_C$ とする。広い変位では中心から離れた領域も動く。

$$
W_C=\frac{C \ell^{2}}{2}
$$

変分の規格化値を $W$ とする。駆動は負符号で入る。

$$
W=- D + W_{C} + W_{T}
$$

最適な幅を $\ell_*$ とする。張力と外側の寄与を同じにする。

$$
\ell_*=\sqrt[4]{\frac{T}{C}}
$$

最小変分を $W_*$ とする。最適な幅を代入する。

$$
W_*=- D + \sqrt{C T}
$$

幅選択によるエネルギーの過剰値を $W-W_*$ とする。非負になる。

$$
W-W_*=W - W_{*}
$$

:::example[中心だけを鋭く動かす]
駆動が強い場合でも短い変位が安定か調べる。

入力値を次のように置く。

$$
\begin{aligned}
T&=1 \\C&=1 \\D&=2 \\\ell&=\frac{1}{4}
\end{aligned}
$$

これらを順に代入する。

張力項は

$$
W_T=\frac{\left(1\right)}{2 \left(\frac{1}{4}\right)^{2}}=8
$$

である。

外側の曲率による安定化は

$$
W_C=\frac{\left(1\right) \left(\frac{1}{4}\right)^{2}}{2}=\frac{1}{32}
$$

である。

変分の規格化値は

$$
W=- \left(2\right) + \left(8\right) + \left(\frac{1}{32}\right)=\frac{193}{32}
$$

である。

最適な幅は

$$
\ell_*=1
$$

である。

最小変分は

$$
W_*=- \left(2\right) + \sqrt{\left(1\right)^{2}}=-1
$$

である。

幅選択によるエネルギーの過剰値は

$$
W-W_*=- \left(-1\right) + \left(\frac{193}{32}\right)=\frac{225}{32}
$$

である。

試験関数を狭くしすぎると張力で正となり、不安定性を見落とす。
:::

:::example[最適幅]
同じ曲率で負の変分を見つける。

入力値を次のように置く。

$$
\begin{aligned}
T&=1 \\C&=1 \\D&=2 \\\ell&=1
\end{aligned}
$$

これらを順に代入する。

張力項は

$$
W_T=\frac{1}{2 \left(1\right)}=\frac{1}{2}
$$

である。

外側の曲率による安定化は

$$
W_C=\frac{\left(1\right)^{3}}{2}=\frac{1}{2}
$$

である。

変分の規格化値は

$$
W=- \left(2\right) + 2 \left(\frac{1}{2}\right)=-1
$$

である。

最適な幅は

$$
\ell_*=1
$$

である。

最小変分は

$$
W_*=- \left(2\right) + \sqrt{\left(1\right)^{2}}=-1
$$

である。

幅選択によるエネルギーの過剰値は

$$
W-W_*=0
$$

である。

この幅では駆動が安定化二項を上回る。
:::

:::example[駆動が弱い場合]
最低の値でも負になるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
T&=1 \\C&=1 \\D&=\frac{1}{2} \\\ell&=1
\end{aligned}
$$

これらを順に代入する。

張力項は

$$
W_T=\frac{1}{2 \left(1\right)}=\frac{1}{2}
$$

である。

外側の曲率による安定化は

$$
W_C=\frac{\left(1\right)^{3}}{2}=\frac{1}{2}
$$

である。

変分の規格化値は

$$
W=\left(\frac{1}{2}\right)=\frac{1}{2}
$$

である。

最適な幅は

$$
\ell_*=1
$$

である。

最小変分は

$$
W_*=- \left(\frac{1}{2}\right) + \sqrt{\left(1\right)^{2}}=\frac{1}{2}
$$

である。

幅選択によるエネルギーの過剰値は

$$
W-W_*=0
$$

である。

この放物型模型の変分は正である。
:::

:::example[張力を増やす]
許される幅と安定性を比較する。

入力値を次のように置く。

$$
\begin{aligned}
T&=16 \\C&=1 \\D&=2 \\\ell&=2
\end{aligned}
$$

これらを順に代入する。

張力項は

$$
W_T=\frac{\left(16\right)}{2 \left(2\right)^{2}}=2
$$

である。

外側の曲率による安定化は

$$
W_C=\frac{\left(1\right) \left(2\right)^{2}}{2}=2
$$

である。

変分の規格化値は

$$
W=\left(2\right)=2
$$

である。

最適な幅は

$$
\ell_*=\sqrt[4]{\frac{\left(16\right)}{\left(1\right)}}=2
$$

である。

最小変分は

$$
W_*=- \left(2\right) + \sqrt{\left(16\right) \left(1\right)}=2
$$

である。

幅選択によるエネルギーの過剰値は

$$
W-W_*=0
$$

である。

最適な幅は広がるが、最小変分も増える。
:::


## 抵抗性不安定性とテアリングモード

磁気凍結が破れる細い層では、外部の理想変位からは許されない磁気島が成長する。これを :term[テアリングモード]{en="tearing mode" yomi="てありんぐもーど"} という。反転場 $B_y(x)$、波数 $k\hat{\vb y}$、共鳴面 $x=0$ を考える。二次元非圧縮模型で $\vb B_\perp=\hat{\vb z}\times\nabla\psi$、$\vb u=\hat{\vb z}\times\nabla\varphi$ とすると、線形誘導式と渦度式は
$$
\begin{aligned}
\gamma\psi_1&=ik B_y\varphi_1+\eta_m(\psi_1''-k^2\psi_1),\\
\rho\gamma(\varphi_1''-k^2\varphi_1)&=\frac{ik}{\mu_0}[B_y(\psi_1''-k^2\psi_1)-B_y''\psi_1].
\end{aligned}
$$
外部で慣性と抵抗を落とせば $\psi_1''-(k^2+B_y''/B_y)\psi_1=0$ となる。両側の減衰解を同じ $\psi_1(0)$ にそろえ、
$$
\Delta'=\frac{\psi_1'(0^+)-\psi_1'(0^-)}{\psi_1(0)}
\label{eq-plasma-tearing-delta}
$$
を定める。$B_y=B_0\tanh(x/a)$ なら $x>0$ の減衰解は $e^{-kx}[ka+\tanh(x/a)]$ に比例する。$(\tanh X)'=\mathrm{sech}^2X$ を二回適用すれば外部式を満たす。左側は偶対称に取れるので $\Delta'a=2[(ka)^{-1}-ka]$ である。

内層では $B_y\simeq B_0 x/a$ とし、幅 $\delta$ 内の $\psi_1''$ を $\Delta'\psi_1/\delta$ と見積もる。増幅、誘導と拡散のつりあいから $\gamma\sim\eta_m\Delta'/\delta$ であり、渦度式から $\gamma^2\sim (kv_A/a)^2\Delta'\delta^3$ となる。消去すると $\gamma a/v_A\sim S_a^{-3/5}(ka)^{2/5}(\Delta'a)^{4/5}$、$S_a=av_A/\eta_m$ である。この定磁束近似は $\Delta'>0$ かつ $\Delta'\delta\ll1$ で使う。この見積もりで未定の係数は、内層の解を外側の勾配に接続して求める。

:::example[しきい値]
$ka=1/2$ では $\Delta'a=3$、$ka=2$ では $\Delta'a=-3$ である。同じ平衡でも波長で外部の駆動が変わる。
:::

### 島の幅と磁束摂動

共鳴面近傍で平衡磁束を $\psi_0=\psi_0(0)+B_y'(0)x^2/2$ とし、摂動を $\widetilde\psi\cos ky$ とする。$B_y'(0)>0$、$\widetilde\psi>0$ なら鞍点の磁束は $\psi_0(0)+\widetilde\psi$ である。島の中心と同じ $y$ にある分離曲線で $B_y'(0)x^2/2-\widetilde\psi=+\widetilde\psi$ だから半幅は $2\sqrt{\widetilde\psi/B_y'(0)}$ である。線形の磁束成長 $\widetilde\psi\propto e^{\gamma t}$ を幅で測ると、幅には $\gamma/2$ が現れる。幅が線形内層を越えると、内層の線形方程式は成立しない。

:::example[磁束と島幅の倍率]
磁束振幅が十六倍になる間に幅は四倍になる。実験で幅の成長時間を求める場合、どの量を指数で当てはめたかを明確にする。
:::

### Harris 層の外側テアリング解

$\vb B=B_0\tanh(x/a)\hat{\vb y}$ とし、無慣性の外側で磁束摂動が $\psi''-(k^2+B_y''/B_y)\psi=0$ を満たす模型を取る。$X=x/a$、$K=ka>0$ として $B_y''/B_y=-2\operatorname{sech}^2X/a^2$ である。右半空間の減衰解 $\psi/\psi(0)=e^{-KX}(K+\tanh X)/K$ は二回微分すればこの式を満たす。偶対称に左側へ接続すると対数微分の飛びは $\Delta'a=2(1/K-K)$ となる。従って $K<1$ 側が正の駆動である。ここでは外側の条件を求めた。成長率は内側の抵抗性層の接続式が与える。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $K$ | 波数 $ka$ | 1 |
| $X$ | 正側の距離 $x/a$ | 1 |

背景の磁場比 $\tanh X$を $t_X$ とする。軸近くでは距離に比例する。

$$
t_X=\frac{e^{2 X} - 1}{e^{2 X} + 1}
$$

減衰解の規格化磁束を $\psi(X)/\psi(0)$ とする。遠方条件を満たす枝を選ぶ。

$$
\psi(X)/\psi(0)=\frac{\left(K + t_{X}\right) e^{- K X}}{K}
$$

右の対数微分を $a\psi'(0^+)/\psi(0)$ とする。$\tanh$ の原点の導関数は一である。

$$
a\psi'(0^+)/\psi(0)=- K + \frac{1}{K}
$$

対数微分の飛びを $\Delta'a$ とする。偶対称なので左の傾きは逆になる。

$$
\Delta'a=2a\psi'(0^+)/\psi(0)
$$

遠方の減衰長を $\ell_k/a$ とする。遠方での指数で決まる。

$$
\ell_k/a=\frac{1}{K}
$$

飛びの有限な比較尺度を $M$ とする。$|\Delta'a|$ が大きいと一に近い。

$$
M=\frac{\left|{\Delta'a}\right|}{\left|{\Delta'a}\right| + 1}
$$

:::example[長波長の正の飛び]
層幅に対して長い波を計算する。

入力値を次のように置く。

$$
\begin{aligned}
K&=\frac{1}{4} \\X&=1
\end{aligned}
$$

これらを順に代入する。

背景の磁場比 $\tanh X$は

$$
t_X=\frac{e^{2 \left(1\right)} - 1}{e^{2 \left(1\right)} + 1}\simeq 0.76159
$$

である。

減衰解の規格化磁束は

$$
\psi(X)/\psi(0)=\frac{\left(\left(\frac{-1 + e^{2}}{1 + e^{2}}\right) + \left(\frac{1}{4}\right)\right) e^{- \left(1\right) \left(\frac{1}{4}\right)}}{\left(\frac{1}{4}\right)}\simeq 3.1513
$$

である。

右の対数微分は

$$
a\psi'(0^+)/\psi(0)=- \left(\frac{1}{4}\right) + \frac{1}{\left(\frac{1}{4}\right)}=\frac{15}{4}
$$

である。

対数微分の飛びは

$$
\Delta'a=2 \left(\frac{15}{4}\right)=\frac{15}{2}
$$

である。

遠方の減衰長は

$$
\ell_k/a=\frac{1}{\left(\frac{1}{4}\right)}=4
$$

である。

飛びの有限な比較尺度は

$$
M=\frac{\left|{\left(\frac{15}{2}\right)}\right|}{\left|{\left(\frac{15}{2}\right)}\right| + 1}=\frac{15}{17}
$$

である。

遠くでは減衰するが、直近の右解は一度増える。
:::

:::example[中立な層]
飛びのゼロを確かめる。

入力値を次のように置く。

$$
\begin{aligned}
K&=1 \\X&=1
\end{aligned}
$$

これらを順に代入する。

背景の磁場比 $\tanh X$は

$$
t_X=\frac{e^{2 \left(1\right)} - 1}{e^{2 \left(1\right)} + 1}\simeq 0.76159
$$

である。

減衰解の規格化磁束は

$$
\psi(X)/\psi(0)=\frac{\left(\left(1\right) + \left(\frac{-1 + e^{2}}{1 + e^{2}}\right)\right) e^{- \left(1\right)^{2}}}{\left(1\right)}\simeq 0.64805
$$

である。

右の対数微分は

$$
a\psi'(0^+)/\psi(0)=- \left(1\right) + \frac{1}{\left(1\right)}=0
$$

である。

対数微分の飛びは

$$
\Delta'a=2 \left(0\right)=0
$$

である。

遠方の減衰長は

$$
\ell_k/a=\frac{1}{\left(1\right)}=1
$$

である。

飛びの有限な比較尺度は

$$
M=\frac{\left|{\left(0\right)}\right|}{\left|{\left(0\right)}\right| + 1}=0
$$

である。

外側解の傾きは連続になり、この駆動は消える。
:::

:::example[短波長の負の飛び]
正の飛びが維持されるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
K&=2 \\X&=1
\end{aligned}
$$

これらを順に代入する。

背景の磁場比 $\tanh X$は

$$
t_X=\frac{e^{2 \left(1\right)} - 1}{e^{2 \left(1\right)} + 1}\simeq 0.76159
$$

である。

減衰解の規格化磁束は

$$
\psi(X)/\psi(0)=\frac{\left(\left(2\right) + \left(\frac{-1 + e^{2}}{1 + e^{2}}\right)\right) e^{- \left(1\right) \left(2\right)}}{\left(2\right)}\simeq 0.18687
$$

である。

右の対数微分は

$$
a\psi'(0^+)/\psi(0)=- \left(2\right) + \frac{1}{\left(2\right)}=- \frac{3}{2}
$$

である。

対数微分の飛びは

$$
\Delta'a=2 \left(- \frac{3}{2}\right)=-3
$$

である。

遠方の減衰長は

$$
\ell_k/a=\frac{1}{\left(2\right)}=\frac{1}{2}
$$

である。

飛びの有限な比較尺度は

$$
M=\frac{\left|{\left(-3\right)}\right|}{\left|{\left(-3\right)}\right| + 1}=\frac{3}{4}
$$

である。

負に変わるため通常の抵抗性テアリングの条件を満たさない。
:::

:::example[層から離れた観測点]
正の飛びと遠方での減衰を両立できるか調べる。

入力値を次のように置く。

$$
\begin{aligned}
K&=\frac{1}{4} \\X&=4
\end{aligned}
$$

これらを順に代入する。

背景の磁場比 $\tanh X$は

$$
t_X=\frac{e^{2 \left(4\right)} - 1}{e^{2 \left(4\right)} + 1}\simeq 0.99933
$$

である。

減衰解の規格化磁束は

$$
\psi(X)/\psi(0)=\frac{\left(\left(\frac{-1 + e^{8}}{1 + e^{8}}\right) + \left(\frac{1}{4}\right)\right) e^{- \left(4\right) \left(\frac{1}{4}\right)}}{\left(\frac{1}{4}\right)}\simeq 1.8384
$$

である。

右の対数微分は

$$
a\psi'(0^+)/\psi(0)=- \left(\frac{1}{4}\right) + \frac{1}{\left(\frac{1}{4}\right)}=\frac{15}{4}
$$

である。

対数微分の飛びは

$$
\Delta'a=2 \left(\frac{15}{4}\right)=\frac{15}{2}
$$

である。

遠方の減衰長は

$$
\ell_k/a=\frac{1}{\left(\frac{1}{4}\right)}=4
$$

である。

飛びの有限な比較尺度は

$$
M=\frac{\left|{\left(\frac{15}{2}\right)}\right|}{\left|{\left(\frac{15}{2}\right)}\right| + 1}=\frac{15}{17}
$$

である。

有限距離では背景関数の因子も残り、十分遠くで指数則に近づく。
:::


### 斜めのテアリングが選ぶ共鳴面

平衡磁場 $\vb B_t(x)$ が $yz$ 面にあり、摂動の接線波数ベクトルを $\vb k_t=(k_y,k_z)$、その大きさを $k>0$ とする。摂動が変化する方向へ座標を回すと、線形の磁力線微分は $i\vb k_t\cdot\vb B_t(x)=ikB_k(x)$ である。$B_k=\hat{\vb k}_t\cdot\vb B_t$ とすると、非圧縮な二次元射影の誘導式と面内の渦度式には $B_y$ の代わりに $B_k$ が入る。波数と直交する接線磁場の力は面内では磁気圧勾配となり、渦度を取ると消える。従って面内の外側式は $\psi''-[k^2+B_k''/B_k]\psi=0$ であり、共鳴面は $B_k=0$ から求める。

例えば $B_y=B_0\tanh(x/a)$ に一定の案内磁場 $B_z$ を足すと、$k_y\ne0$ かつ $|k_zB_z/(k_yB_0)|<1$ で $x_s=a\,\operatorname{artanh}[-k_zB_z/(k_yB_0)]$ に共鳴面が動く。そこで $B_k''(x_s)$ が一般には非零であるため、外側の左右解は対称とは限らず、勾配に共通の対数項が現れることもある。実際、$x-x_s$ を $X$ とし、$B_k''/B_k\simeq (B_k''(x_s)/B_k'(x_s))/X$ とすると、$\psi''$ の先頭は定数を $X$ で割った形である。積分した左右の微分から共通の $\ln|X|$ を除けば有限部分の跳びが得られる。以下の内層計算は、その対数項が消える反転面の模型で行う。

### 定磁束内層の接続係数

外側の $\Delta'$ から成長率を決めるには、内層が作る磁束勾配の跳びを求める。Furth、Killeen、Rosenbluth の定磁束近似では、この跳びと外部解を使って幅と成長率を同時に定める。磁束摂動を $\sqrt{\mu_0\rho}$ で割った量を $\Psi$、流れ関数の振幅に $i$ を掛けた量を $F$ とする。$B_y$ が原点を正の傾きで横切る向きを選び、$A=kB'_y(0)/\sqrt{\mu_0\rho}>0$ とおく。磁気拡散係数を $\eta_m$、実成長率を $\gamma>0$ とする。内層で $k\delta\ll1$ とし、$B_y''$ の項も外すと
$$
\gamma\Psi=A xF+\eta_m\Psi'',\qquad
\gamma F''=-Ax\Psi'' .
$$
$\Delta'\delta\ll1$ では第一式の左辺で $\Psi$ を $\Psi_0=\Psi(0)$ に置ける。ただし $\Psi''$ は勾配の跳びを担うので残す。速度変数を消さずに $\Psi''$ を消去すると $\gamma\eta_m F''-A^2x^2F=-A x\gamma\Psi_0$ である。そこで
$$
\delta=\left(\frac{\gamma\eta_m}{A^2}\right)^{1/4},
\qquad X=\frac{x}{\delta},
\qquad F=\frac{\gamma\Psi_0}{A\delta}U(X)
$$
とおけば $U''-X^2U=-X$ となる。外側で誘導と移流を合わせる条件は $U(X)\sim1/X$ であり、解は奇関数とする。

:::lemma[内層の奇解]{#lem-plasma-fkr-layer}
$U(X)=X\int_0^{1/2}(1-4s^2)^{-1/4}e^{-sX^2}\,\dd s$ は、$U''-X^2U=-X$ と $U(X)\sim1/X$ を満たす滑らかな奇解である。
:::

:::proof
$g=Xe^{-sX^2}$ とおけば $\partial_X^2g-X^2g=-6s\,g+(1-4s^2)\partial_sg$ である。重み $w=(1-4s^2)^{-1/4}$ を掛けて $s$ で部分積分する。$(w(1-4s^2))'=-6s w$ なので積分の内部は消える。上端で $w(1-4s^2)=0$、下端で $w(1-4s^2)g=X$ だから、境界項は $-X$ である。大きな $|X|$ では重みの $s=0$ 近傍が支配し、$U=X\int_0^\infty e^{-sX^2}\dd s+\order{|X|^{-5}}=1/X+\order{|X|^{-5}}$ となる。特に $1-XU$ は実軸で可積分である。もう一つ同じ境界条件を満たす解があると、その差 $V$ は同次式と減衰条件を満たす。$V V''-X^2V^2$ を積分すると $-\int[(V')^2+X^2V^2]\dd X=0$ となり、差はゼロである。
:::

第一内層式に :ref[lem-plasma-fkr-layer] の解を入れ、$\Psi''=(\gamma\Psi_0/\eta_m)(1-XU)$ を共鳴面にわたって積分する。外側の勾配から求めた :eqref[eq-plasma-tearing-delta] と一致させると、$\Delta'=\gamma\delta C/\eta_m$、$C=\int_{\R}(1-XU)\,\dd X$ である。$X^2\int_0^\infty e^{-sX^2}\dd s=1$ と Gauss 積分を使えば
$$
C=\frac{\sqrt\pi}{2}
\left[
\int_{1/2}^{\infty}s^{-3/2}\,\dd s
+
\int_0^{1/2}\bigl[1-(1-4s^2)^{-1/4}\bigr]s^{-3/2}\,\dd s
\right].
$$
二番目の積分は $s=0$ で差が $s^2$ のオーダーであり、上端の特異性も可積分である。$y=4s^2$ と置くと括弧内は $2^{-1/2}$ 掛ける
$$
4+\int_0^1\bigl[1-(1-y)^{-1/4}\bigr]y^{-5/4}\,\dd y
$$
となる。$g_1=y^{-1/4}(1-y)^{3/4}$ には $g_1'=-(1+2y)y^{-5/4}(1-y)^{-1/4}/4$ が成り立つ。積分の下限を $\epsilon>0$ にそろえて $g_1'$ を積分すると、$\epsilon^{-1/4}$ の発散は差を取った時点で消える。残る値は $2\int_0^1 y^{-1/4}(1-y)^{-1/4}\dd y=2\mathrm B(3/4,3/4)$ である。ここで $\mathrm B$ は Beta 関数である。従って
$$
C=\sqrt{\frac{\pi}{2}}\,\mathrm B(3/4,3/4)
=\sqrt2\,\Gamma(3/4)^2>0 .
$$
$\Gamma$ は Gamma 関数である。$\delta$ を戻し $\gamma$ について解くと
$$
\gamma
=\eta_m^{3/5}A^{2/5}\left(\frac{\Delta'}{C}\right)^{4/5},
\qquad C\simeq 2.124 .
\label{eq-plasma-fkr-match}
$$
この式の次元は $\mathrm{s^{-1}}$ である。$\Delta'>0$ が成長解に必要であることも積分の符号から分かる。:eqref[eq-plasma-fkr-match] は定磁束の最初の次数であり、有限の $\Delta'\delta$、高い粘性、Hall 項を許す接続式とは異なる。外側が与える $\Delta'a$ と内側が与える $\delta/a$ の積を、解を得た後で確認する。

:::example[Harris 層の一波数]
$ka=1/2$ では $\Delta'a=3$ である。$S_a=10^6$ とし、$Aa^2/v_A=ka$ を使うと $\gamma a/v_A=S_a^{-3/5}(ka)^{2/5}(\Delta'a/C)^{4/5}\simeq2.51\times10^{-4}$ となる。幅は $\delta/a=[(\gamma a/v_A)/(S_a(ka)^2)]^{1/4}\simeq5.63\times10^{-3}$ である。従って $\Delta'\delta\simeq0.017$ は一より小さく、接続に使った定磁束の尺度分離がある。
:::

:::example[近い波数での導電性の比較]
$ka$ と平衡を同じにして $\eta_m$ を四分の一にすると、:eqref[eq-plasma-fkr-match] により成長率は $4^{-3/5}\simeq0.435$ 倍となる。$\delta\propto(\gamma\eta_m)^{1/4}$ だから層幅も $4^{-2/5}\simeq0.574$ 倍になる。磁気拡散が弱まると、成長が遅くなるのと同時に接続を担う領域が狭まる。
:::

:::example[長波長へ延長する限界]
Harris 層の $ka\ll1$ で $\Delta'a\simeq2/(ka)$ を代入すると、$\delta/a\propto S_a^{-2/5}(ka)^{-3/5}$ となる。従って $\Delta'\delta\propto S_a^{-2/5}(ka)^{-8/5}$ である。固定の $S_a$ で $ka$ を小さくすると定磁束の仮定が先に破れる。そのため成長率の長波長発散をこの式から結論できない。
:::

### 島が線形内層より大きくなった後

島の半幅 $w$ が線形内層の幅より大きくなっても、$\Delta' w\ll1$ なら外側解は先頭近似で保てる。線形の指数成長から離れ、島内の流体慣性を落とせる抵抗性の段階では、外側から与えられる磁束勾配の跳びは $\Delta'\widetilde\psi$ である。これを幅 $w$ にわたる電流で支えると、磁束の拡散項の尺度は $\eta_m\Delta'\widetilde\psi/w$ となる。他方、半幅と磁束の関係は $\widetilde\psi=B_y'(0)w^2/4$ である。従って磁束の変化率を比較すると $\dot w\sim\eta_m\Delta'$ を得る。この Rutherford 段階では、$\Delta'$ がほぼ一定なら幅は時間にほぼ比例する。係数を得るには島内の磁気面ごとの誘導式を平均し、分離曲線で外側へ接続する。幅が平衡の変化長に近づいた場合は、非線形の外側解も必要である。

## 磁気再結合と Sweet–Parker 模型

反対向きの磁場が接する層へ流れが入り、磁力線がつなぎ替わる現象を :term[磁気再結合]{en="magnetic reconnection" yomi="じきさいけつごう"} という。二次元定常で、長さの半幅 $L$、厚さの半幅 $\delta$ の層を考える。非圧縮なら $v_{\rm in}L=v_{\rm out}\delta$ である。磁気圧を流出の動圧に変えると $v_{\rm out}\sim v_A$ である。

Faraday 則から面に垂直な電場は一定である。層の外では $E\sim v_{\rm in}B$、中心では $E\sim\eta B/(\mu_0\delta)$ である。よって
$$
\frac{\delta}{L}\sim\frac{v_{\rm in}}{v_A}\sim S_L^{-1/2},\qquad S_L=Lv_A/\eta_m.
\label{eq-plasma-sweet-parker}
$$
これが :term[Sweet–Parker 模型]{en="Sweet–Parker model" yomi="すうぃーとぱーかーもけい"} である。層の半幅を両方に使ったので質量式の入口・出口の数が相殺される。

:::example[大きい Lundquist 数]
$S_L=10^8$ なら $v_{\rm in}/v_A\sim10^{-4}$ である。吸い込みが遅いのは流出路が狭く、流入した質量を排出できる速さに上限があるためである。
:::

:::example[層の時間の比較]
Sweet–Parker の幅では $\delta^2/\eta_m\sim L/v_A$ である。横拡散の時間が長さ方向の流出時間に合う。$S_L=10^8$ でも幅の小ささにより、拡散を全長 $L$ に対する拡散時間だけで捨てられない。
:::

### Sweet–Parker 層の時間尺度と局所条件

半長 $L$、半幅 $\delta$ の定常な二次元層で密度一定とする。質量保存は $v_{\rm in} L=v_{\rm out}\delta$、一様な再結合電場は $v_{\rm in}B\simeq\eta_m B/\delta$ である。$v_{\rm out}\simeq v_A$ とおくと $\delta=\sqrt{\eta_mL/v_A}$ であり、両辺の電場は $Bv_A/\sqrt{S}$ となる。$\delta^2/\eta_m=L/v_A$ は層を横切る拡散が流出時間に匹敵することを示す。さらに $\delta/d_i$ を調べ、薄い場合は一流体の抵抗模型の妥当性を評価する。極端に大きな $S$ では層の不安定性で定常単層が維持されない場合もある。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $L$ | 層の半長 | $\mathrm m$ |
| $v_A$ | 上流の Alfvén 速度 | $\mathrm{m\,s^{-1}}$ |
| $\eta_m$ | 磁気拡散係数 | $\mathrm{m^2\,s^{-1}}$ |
| $d_i$ | イオン慣性長 | $\mathrm m$ |

Lundquist 数を $S$ とする。流出時間と全長の拡散時間の比をとる。

$$
S=\frac{L v_{A}}{\eta_{m}}
$$

層の半幅を $\delta$ とする。電場と質量の収支を合わせる。

$$
\delta=\frac{L}{\sqrt{S}}
$$

流入速度を $v_{\rm in}$ とする。流出速度との比は半幅と半長の比である。

$$
v_{\rm in}=\frac{v_{A}}{\sqrt{S}}
$$

流出時間を $t_{\rm out}$ とする。長さを Alfvén 速度で分ける。

$$
t_{\rm out}=\frac{L}{v_{A}}
$$

半幅での拡散時間を $t_{\rm diff}$ とする。流出時間と等しい。

$$
t_{\rm diff}=\frac{\delta^{2}}{\eta_{m}}
$$

一流体性を吟味する長さ比を $\delta/d_i$ とする。これが小さいと Hall 項を無視しにくい。

$$
\delta/d_i=\frac{\delta}{d_{i}}
$$

:::example[十分広い抵抗層]
半幅を計算してイオン尺度と比較する。

入力値を次のように置く。

$$
\begin{aligned}
L&=100\,\mathrm m \\v_A&=100\,\mathrm{m\,s^{-1}} \\\eta_m&=1\,\mathrm{m^2\,s^{-1}} \\d_i&=\frac{1}{100}\,\mathrm m
\end{aligned}
$$

これらを順に代入する。

Lundquist 数は

$$
S=\frac{\left(100\right)^{2}}{\left(1\right)}=10000
$$

である。

層の半幅は

$$
\delta=\frac{\left(100\right)}{\sqrt{\left(10000\right)}}=1\,\mathrm m
$$

である。

流入速度は

$$
v_{\rm in}=\frac{\left(100\right)}{\sqrt{\left(10000\right)}}=1\,\mathrm{m\,s^{-1}}
$$

である。

流出時間は

$$
t_{\rm out}=1\,\mathrm s
$$

である。

半幅での拡散時間は

$$
t_{\rm diff}=\left(1\right)=1\,\mathrm s
$$

である。

一流体性を吟味する長さ比は

$$
\delta/d_i=\frac{\left(1\right)}{\left(\frac{1}{100}\right)}=100
$$

である。

層は細いがイオン慣性長より大きい。
:::

:::example[イオン慣性長と同程度]
大域の $S$ を保って局所条件を変える。

入力値を次のように置く。

$$
\begin{aligned}
L&=100\,\mathrm m \\v_A&=100\,\mathrm{m\,s^{-1}} \\\eta_m&=1\,\mathrm{m^2\,s^{-1}} \\d_i&=1\,\mathrm m
\end{aligned}
$$

これらを順に代入する。

Lundquist 数は

$$
S=\frac{\left(100\right)^{2}}{\left(1\right)}=10000
$$

である。

層の半幅は

$$
\delta=\frac{\left(100\right)}{\sqrt{\left(10000\right)}}=1\,\mathrm m
$$

である。

流入速度は

$$
v_{\rm in}=\frac{\left(100\right)}{\sqrt{\left(10000\right)}}=1\,\mathrm{m\,s^{-1}}
$$

である。

流出時間は

$$
t_{\rm out}=1\,\mathrm s
$$

である。

半幅での拡散時間は

$$
t_{\rm diff}=\left(1\right)=1\,\mathrm s
$$

である。

一流体性を吟味する長さ比は

$$
\delta/d_i=1
$$

である。

定常流入の形式値は同じだが Hall 項を見直す必要がある。
:::

:::example[より低い抵抗]
抵抗が小さくなったときに局所の尺度を求める。

入力値を次のように置く。

$$
\begin{aligned}
L&=100\,\mathrm m \\v_A&=100\,\mathrm{m\,s^{-1}} \\\eta_m&=\frac{1}{100}\,\mathrm{m^2\,s^{-1}} \\d_i&=1\,\mathrm m
\end{aligned}
$$

これらを順に代入する。

Lundquist 数は

$$
S=\frac{\left(100\right)^{2}}{\left(\frac{1}{100}\right)}=1000000
$$

である。

層の半幅は

$$
\delta=\frac{\left(100\right)}{\sqrt{\left(1000000\right)}}=\frac{1}{10}\,\mathrm m
$$

である。

流入速度は

$$
v_{\rm in}=\frac{\left(100\right)}{\sqrt{\left(1000000\right)}}=\frac{1}{10}\,\mathrm{m\,s^{-1}}
$$

である。

流出時間は

$$
t_{\rm out}=1\,\mathrm s
$$

である。

半幅での拡散時間は

$$
t_{\rm diff}=\frac{\left(\frac{1}{10}\right)^{2}}{\left(\frac{1}{100}\right)}=1\,\mathrm s
$$

である。

一流体性を吟味する長さ比は

$$
\delta/d_i=\frac{\left(\frac{1}{10}\right)}{\left(1\right)}=\frac{1}{10}
$$

である。

拡散が弱いため半幅が短くなり、同じ時間での拡散の収支を保つ。
:::

:::example[層が長くなる場合]
半幅と相対幅を別々に比較する。

入力値を次のように置く。

$$
\begin{aligned}
L&=10000\,\mathrm m \\v_A&=100\,\mathrm{m\,s^{-1}} \\\eta_m&=1\,\mathrm{m^2\,s^{-1}} \\d_i&=\frac{1}{10}\,\mathrm m
\end{aligned}
$$

これらを順に代入する。

Lundquist 数は

$$
S=\frac{\left(10000\right) \left(100\right)}{\left(1\right)}=1000000
$$

である。

層の半幅は

$$
\delta=\frac{\left(10000\right)}{\sqrt{\left(1000000\right)}}=10\,\mathrm m
$$

である。

流入速度は

$$
v_{\rm in}=\frac{\left(100\right)}{\sqrt{\left(1000000\right)}}=\frac{1}{10}\,\mathrm{m\,s^{-1}}
$$

である。

流出時間は

$$
t_{\rm out}=\frac{\left(10000\right)}{\left(100\right)}=100\,\mathrm s
$$

である。

半幅での拡散時間は

$$
t_{\rm diff}=\frac{\left(10\right)^{2}}{\left(1\right)}=100\,\mathrm s
$$

である。

一流体性を吟味する長さ比は

$$
\delta/d_i=\frac{\left(10\right)}{\left(\frac{1}{10}\right)}=100
$$

である。

絶対幅は増すが相対幅は減り、流入は遅くなる。
:::


## 高速磁気再結合と Petschek 模型

非理想領域が短く、その外に開いた流出口があれば、長い薄い層による制限を弱められる。:term[Petschek 模型]{en="Petschek model" yomi="ぺちぇっくもけい"} は短い拡散層の外に遅い MHD 衝撃波を置く。衝撃波を斜めに横切る流体は、磁場の接線成分の減少を通じてエネルギーを流出と熱へ変換する。低い流入 Mach 数では、衝撃波の開き角の目安は $v_{\rm in}/v_A$ である。拡散域長を $\ell<L$ として前節の局所収支を取ると $v_{\rm in}/v_A\sim (\eta_m/(\ell v_A))^{1/2}$ となる。ただしこの式だけで $\ell$ は決まらない。

Petschek の対数的な速度見積もりには外部場と衝撃波をつないで磁場の減少を求める必要がある。定数の抵抗率で短い層が維持できるとは限らず、生成される横磁場と境界条件も関わる。外部場と衝撃波の接続を以下で見積もる。大きな系では長い層が複数の磁気島に分かれることがあり、無衝突系では Hall 項、電子圧力テンソルと慣性を区別して評価する。

:::example[局所長の効果]
同じ $\eta_m$ と $v_A$ で $\ell=L/100$ を仮定すると、局所収支の流入速度は長い層の十倍になる。この仮定が保たれるかは外部との接続を別に検証する。
:::

:::example[律速のしきい値]
Sweet–Parker の $\delta$ とイオン慣性長を比較する。$L=10^4\,d_i$、$S_L=10^{10}$ なら $\delta\simeq0.1 d_i$ である。抵抗性一流体の幅がこの尺度に入るなら、Hall 項を最初から無視して流入を見積もる手順を点検する。
:::

### 外部磁場との接続で現れる対数

短い拡散層に接続する外部領域を上半面 $y>0$ とし、係数を捨てた二次元の真空近似で対数を求める。外部の摂動磁場を $\delta\vb B=-\nabla\Phi$、$\nabla^2\Phi=0$ とする。$x$ 方向 Fourier 成分について減衰解は $\Phi(k,y)=\Phi(k,0)e^{-\abs k y}$ である。従って $\delta B_x(k,0)=-i\sgn(k)\delta B_y(k,0)$ となる。この変換を実空間で書くと
$$
\delta B_x(0,0)=\frac{1}{\pi}\operatorname{PV}\int_\R\frac{\delta B_y(x,0)}{-x}\dd x
$$
である。ここで核 $1/(\pi x)$ の Fourier 変換を求めるには、正弦積分に収束因子 $e^{-\alpha |x|}$ を掛ける。$\int_0^\infty e^{-\alpha x}\sin(kx)\dd x/x=\arctan(k/\alpha)$ を $\alpha\to0^+$ とする。

流出口の両側で $\ell<\abs x<L$ にわたって $\delta B_y(x,0)\simeq b_\perp\sgn(x)$ と置けば、中心では $\delta B_x\simeq-(2b_\perp/\pi)\ln(L/\ell)$ である。開き角 $v_{\rm in}/v_A$ の流出口では $b_\perp/B$ はこの角と同程度である。従って中心の上流磁場を大きく減らさずに短い層を維持できる流入速度の尺度には $1/\ln(L/\ell)$ が現れる。$\ell$ と流入速度の関係、衝撃波の跳び、外部境界を同時に解かない限り係数には意味がない。この推定は領域の短縮が許されることを仮定しており、抵抗率が一定の系でその仮定を自動的に満たさない。

:::example[層長を変えたときの違い]
他の条件を保ち $L/\ell$ を $10^2$ から $10^4$ に増やすと、同じ $b_\perp$ による中心磁場の減少量はこの外部模型で二倍となる。局所の Sweet–Parker 収支だけから予想した速度増加には、外部場の変化という制約が加わる。
:::

### 遅い衝撃波の流出と横磁場の供給

局所的に平らな遅い衝撃波を考え、衝撃面に沿う電場がゼロになる de Hoffmann–Teller 系を選ぶ。法線成分を添字 $n$、接線成分を $t$ で表す。理想上下流では $u_t=u_n B_t/B_n$、法線磁束 $B_n$ と質量流束 $m=\rho u_n$ は連続である。接線運動量の保存式は
$$
\left[\left(\frac{m^2}{\rho B_n}-\frac{B_n}{\mu_0}\right)B_t\right]=0
$$
となる。$[A]=A_2-A_1$ とし、反転する接線磁場が下流で消える極限を取る。$B_{t1}\ne0$、$B_n\ne0$ なら $u_{n1}^2=B_n^2/(\mu_0\rho_1)$ である。そのとき上流の接線速度の大きさは $|B_{t1}|/\sqrt{\mu_0\rho_1}$ になる。したがって $|u_{n1}|\ll |u_{t1}|$ の小角極限では、衝撃波は小さな法線磁場で Alfvén 速度に相当する接線速度差を支える。圧力・密度の跳びを求めるには、さらに法線運動量と総エネルギーの保存が必要である。

短い拡散層を保つには、流出が運び去る横磁場を層内で補う必要がある。たとえば層に沿う上流磁場が対称で、長さ $L$ で滑らかに変わる場合を考える。一様抵抗率かつ厚さの変化が小さい層では、長さ $\ell\ll L$ の区間で吸い込む磁束の場所による差は相対的に $\order{(\ell/L)^2}$ である。横磁場の生成率は大きさで $Bv_{\rm in}(\ell/L)^2/\ell$、流出による損失率は $v_A b_\perp/\ell$ である。衝撃波に必要な $b_\perp/B\sim v_{\rm in}/v_A$ を入れると $\ell/L\sim1$ となる。同じ近似のもとで一様抵抗率の短層が維持できない理由である。層の局所化を外から課す場合にはこの供給の見積もりを組み直す。

### 局所拡散層と外側の対数制約

Petschek 型の二次元模型で長い全長 $L$ に対し抵抗層長を $\ell\ll L$ とする。局所の質量・電場収支から $M=v_{\rm in}/v_A\sim(\eta_m/(v_A\ell))^{1/2}$ である。一方、外側の小角の遅い衝撃波が上流場を弱める効果は、局所層から全長までの Green 関数により $\propto M\int_\ell^L \dd r/r=M\ln(L/\ell)$ となる。許す上流磁場の低下に上限を置けば $M\lesssim c/\ln(L/\ell)$ の尺度になる。以下ではこの二つを独立に計算し、整合性を調べる。係数 $c$ は境界条件に依存する。均一抵抗で短い局所層を課すだけでは、その層が自己無撞着に存在するとは言えない。

計算には次の量を指定する。例題では、単位を記した量にその単位を用いる。

| 記号 | 意味 | 単位 |
| --- | --- | --- |
| $S$ | 全長の Lundquist 数 | 1 |
| $x$ | 局所層と全長の比 $\ell/L<1$ | 1 |
| $c$ | 外側の上限を表す係数 | 1 |

局所の Lundquist 数を $S_\ell$ とする。短い層を基準にすると小さくなる。

$$
S_\ell=S x
$$

局所収支からの流入比を $M_{\rm loc}$ とする。流出は Alfvén 速度とした。

$$
M_{\rm loc}=\frac{1}{\sqrt{S_{\ell}}}
$$

外側の上限尺度を $M_{\rm ext}$ とする。対数積分の制約を表す。

$$
M_{\rm ext}=\frac{c}{\ln{\left(\frac{1}{x} \right)}}
$$

外側の制約に対する局所比を $C_{\rm ext}$ とする。一より大きければ指定の局所層と整合しない。

$$
C_{\rm ext}=\frac{M_{\rm loc}}{M_{\rm ext}}
$$

全長を層とした場合の流入比を $M_{\rm SP}$ とする。同じ媒質を基準にする。

$$
M_{\rm SP}=\frac{1}{\sqrt{S}}
$$

層を短くしたときの形式的な増倍率を $M_{\rm loc}/M_{\rm SP}$ とする。制約との比較を済ませて解釈する。

$$
M_{\rm loc}/M_{\rm SP}=\frac{M_{\rm loc}}{M_{\rm SP}}
$$

:::example[外部の制約内にある短層]
局所と外側の両方の比を求める。

入力値を次のように置く。

$$
\begin{aligned}
S&=1000000 \\x&=\frac{1}{100} \\c&=\frac{1}{5}
\end{aligned}
$$

これらを順に代入する。

局所の Lundquist 数は

$$
S_\ell=\left(1000000\right) \left(\frac{1}{100}\right)=10000
$$

である。

局所収支からの流入比は

$$
M_{\rm loc}=\frac{1}{\sqrt{\left(10000\right)}}=\frac{1}{100}
$$

である。

外側の上限尺度は

$$
M_{\rm ext}=\frac{\left(\frac{1}{5}\right)}{\ln{\left(\frac{1}{\left(\frac{1}{100}\right)} \right)}}\simeq 0.043429
$$

である。

外側の制約に対する局所比は

$$
C_{\rm ext}=\frac{\left(\frac{1}{100}\right)}{\left(\frac{1}{5 \ln{\left(100 \right)}}\right)}\simeq 0.23026
$$

である。

全長を層とした場合の流入比は

$$
M_{\rm SP}=\frac{1}{\sqrt{\left(1000000\right)}}=\frac{1}{1000}
$$

である。

層を短くしたときの形式的な増倍率は

$$
M_{\rm loc}/M_{\rm SP}=\frac{\left(\frac{1}{100}\right)}{\left(\frac{1}{1000}\right)}=10
$$

である。

この入力では外部の上限を下回る。ただし接続が実現することの証明にはならない。
:::

:::example[短くしすぎた層]
層を縮め続けることが可能か調べる。

入力値を次のように置く。

$$
\begin{aligned}
S&=1000000 \\x&=\frac{1}{10000} \\c&=\frac{1}{5}
\end{aligned}
$$

これらを順に代入する。

局所の Lundquist 数は

$$
S_\ell=\left(1000000\right) \left(\frac{1}{10000}\right)=100
$$

である。

局所収支からの流入比は

$$
M_{\rm loc}=\frac{1}{\sqrt{\left(100\right)}}=\frac{1}{10}
$$

である。

外側の上限尺度は

$$
M_{\rm ext}=\frac{\left(\frac{1}{5}\right)}{\ln{\left(\frac{1}{\left(\frac{1}{10000}\right)} \right)}}\simeq 0.021715
$$

である。

外側の制約に対する局所比は

$$
C_{\rm ext}=\frac{\left(\frac{1}{10}\right)}{\left(\frac{1}{5 \ln{\left(10000 \right)}}\right)}\simeq 4.6052
$$

である。

全長を層とした場合の流入比は

$$
M_{\rm SP}=\frac{1}{\sqrt{\left(1000000\right)}}=\frac{1}{1000}
$$

である。

層を短くしたときの形式的な増倍率は

$$
M_{\rm loc}/M_{\rm SP}=\frac{\left(\frac{1}{10}\right)}{\left(\frac{1}{1000}\right)}=100
$$

である。

局所の流入が外側の上限を上回り、両条件を同時に満たせない。
:::

:::example[大きな全長の Lundquist 数]
同じ幾何で抵抗を減らした場合を求める。

入力値を次のように置く。

$$
\begin{aligned}
S&=100000000 \\x&=\frac{1}{10000} \\c&=\frac{1}{5}
\end{aligned}
$$

これらを順に代入する。

局所の Lundquist 数は

$$
S_\ell=\left(100000000\right) \left(\frac{1}{10000}\right)=10000
$$

である。

局所収支からの流入比は

$$
M_{\rm loc}=\frac{1}{\sqrt{\left(10000\right)}}=\frac{1}{100}
$$

である。

外側の上限尺度は

$$
M_{\rm ext}=\frac{\left(\frac{1}{5}\right)}{\ln{\left(\frac{1}{\left(\frac{1}{10000}\right)} \right)}}\simeq 0.021715
$$

である。

外側の制約に対する局所比は

$$
C_{\rm ext}=\frac{\left(\frac{1}{100}\right)}{\left(\frac{1}{5 \ln{\left(10000 \right)}}\right)}\simeq 0.46052
$$

である。

全長を層とした場合の流入比は

$$
M_{\rm SP}=\frac{1}{\sqrt{\left(100000000\right)}}=\frac{1}{10000}
$$

である。

層を短くしたときの形式的な増倍率は

$$
M_{\rm loc}/M_{\rm SP}=\frac{\left(\frac{1}{100}\right)}{\left(\frac{1}{10000}\right)}=100
$$

である。

局所比は小さくなるが全長の模型より速いという形式的な比は保たれる。
:::

:::example[長さの分離が弱い層]
対数の数値が小さい場合を検討する。

入力値を次のように置く。

$$
\begin{aligned}
S&=1000000 \\x&=\frac{1}{2} \\c&=\frac{1}{5}
\end{aligned}
$$

これらを順に代入する。

局所の Lundquist 数は

$$
S_\ell=\left(1000000\right) \left(\frac{1}{2}\right)=500000
$$

である。

局所収支からの流入比は

$$
M_{\rm loc}=\frac{1}{\sqrt{\left(500000\right)}}\simeq 0.0014142
$$

である。

外側の上限尺度は

$$
M_{\rm ext}=\frac{\left(\frac{1}{5}\right)}{\ln{\left(\frac{1}{\left(\frac{1}{2}\right)} \right)}}\simeq 0.28854
$$

である。

外側の制約に対する局所比は

$$
C_{\rm ext}=\frac{\left(\frac{\sqrt{2}}{1000}\right)}{\left(\frac{1}{5 \ln{\left(2 \right)}}\right)}\simeq 0.0049013
$$

である。

全長を層とした場合の流入比は

$$
M_{\rm SP}=\frac{1}{\sqrt{\left(1000000\right)}}=\frac{1}{1000}
$$

である。

層を短くしたときの形式的な増倍率は

$$
M_{\rm loc}/M_{\rm SP}=\frac{\left(\frac{\sqrt{2}}{1000}\right)}{\left(\frac{1}{1000}\right)}\simeq 1.4142
$$

である。

外側の対数項だけが支配的という近似が弱いので、上限の値を精密な判定にしない。
:::


## 演習問題

::::exercise{level=basic}
重い流体が下にある無磁場の平面境界で成長するか。

:::solution
$\rho_+-\rho_-<0$ なので式 :eqref[eq-plasma-rt] の $\gamma^2<0$ であり、摂動は振動する。
:::
::::

::::exercise{level=basic}
$ka=1$ の反転場模型の $\Delta'a$ を求めよ。

:::solution
$2(1/1-1)=0$ である。
:::
::::

::::exercise{level=basic}
$S_L=10^6$ の Sweet–Parker 層の $\delta/L$ を求めよ。

:::solution
$S_L^{-1/2}=10^{-3}$ である。
:::
::::

::::exercise{level=basic}
$m=0$ の摂動の方位角依存性を述べよ。

:::solution
$e^{im\theta}=1$ なので全周で変位が同相である。動径成分は軸対称に半径を変える。
:::
::::

::::exercise{level=advanced}
両側に同じ磁場 $B\hat{\vb x}$、$\vb k=k\hat{\vb x}$ のとき重い上側を安定化する最小波数を求めよ。

:::solution
駆動と張力を等置すると $gk(\rho_+-\rho_-)=2k^2B^2/\mu_0$ である。$k>0$ を割って $k_c=\mu_0g(\rho_+-\rho_-)/(2B^2)$ となり、$k\ge k_c$ では成長しない。
:::
::::

::::exercise{level=advanced}
定磁束の二つの釣り合い式から内層幅の指数を求めよ。

:::solution
$\eta_m^2\Delta'^2/\delta^2\sim(kv_A/a)^2\Delta'\delta^3$ なので $\delta^5\sim\eta_m^2\Delta'a^2/(k^2v_A^2)$ である。従って $\delta/a\sim S_a^{-2/5}(ka)^{-2/5}(\Delta'a)^{1/5}$ となる。
:::
::::

::::exercise{level=advanced}
電気抵抗率が四倍のとき Sweet–Parker の流入速度は何倍か。

:::solution
$L,v_A$ 固定では $\eta_m$ は四倍、$S_L$ は四分の一となる。$S_L^{-1/2}$ は二倍である。
:::
::::

::::exercise{level=advanced}
$m,n$ の位相 $m\theta-2\pi n z/L$ が円柱磁力線に沿って一定となる条件を求めよ。

:::solution
$\dd\theta/\dd z=B_\theta/(r B_z)$ である。位相の微分ゼロから $m B_\theta/(r B_z)=2\pi n/L$、$n\ne0$ かつ $B_zB_\theta\ne0$ なら、従って $q=m/n$ となる。$n=0$ の場合は元の式から $m B_\theta=0$ と判定する。
:::
::::

## 確認テスト

::::quiz
磁束凍結が成り立つ滑らかな理想 MHD の解から磁気再結合を説明するには、どの仮定を見直す必要があるか。

- [x] 磁束保存を破る項または理想模型の特異性
- [ ] 大きな磁場だけ
- [ ] 圧力をゼロにすること

:::explanation
滑らかな理想流れでは物質曲面の磁束が保存される。連結を変える過程には非理想層などが必要である。
:::
::::

::::quiz{answer="3" tol="1e-6"}
$ka=1/2$ の $\Delta'a$ を答えよ。

:::explanation
外部解で $2[(ka)^{-1}-ka]=2(2-1/2)=3$ である。
:::
::::

::::quiz
上下の接線磁場が同一直線に沿う平面境界で、磁気張力を避ける波数を選べ。

- [x] 両側の磁場に垂直な波数
- [ ] どちらかの磁場に平行な波数
- [ ] 磁場と無関係な任意の波数

:::explanation
式 :eqref[eq-plasma-rt] の張力は二つの内積の二乗の和なので、双方との内積がゼロのとき消える。上下で磁場が平行でなければ、界面内にそのような非零の波数は存在しない。
:::
::::

::::quiz{answer="0.01" tol="1e-6"}
$S_L=10^4$ の流入 Alfvén Mach 数を Sweet–Parker 近似で答えよ。

:::explanation
式 :eqref[eq-plasma-sweet-parker] から $10^{-2}$ である。
:::
::::

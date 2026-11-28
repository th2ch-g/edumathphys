---
title: 核反応
description: 二体反応の運動学、散乱の流束、複合核、直接反応の解析を扱う。
notation:
  - symbol: 'Q'
    meaning: 入口と出口の静止エネルギー差
  - symbol: 'S_\ell'
    meaning: 部分波の弾性散乱行列要素
  - symbol: '\Gamma_c'
    meaning: チャネル c の部分幅
  - symbol: 'T_c'
    meaning: 統計模型でのチャネル透過係数
---

## この章の目標

- 質量差から反応のしきい値と生成物の運動量を計算する。
- 角運動量と部分波で散乱・反応の断面積を整理する。
- 複合核の共鳴と平均断面積を区別する。
- 核子移行反応と重イオン反応に現れる観測量を説明する。

## 反応の運動学と Q 値

$a+A\to b+B$ を $A(a,b)B$ と書く。大文字は核、小文字は飛跡を測る粒子などを表す記号として使う。核の励起状態を指定するならその静止質量に励起エネルギーを含める。:term[反応の Q 値]{en="reaction Q value" yomi="はんのうのきゅーち"} は $Q=(m_a+m_A-m_b-m_B)c^2$ である。これにより非相対論的な重心系で出口の相対運動エネルギーは $E_{\rm cm}+Q$ であり、出口の運動量は $p_f=\sqrt{2\mu_{bB}(E_{\rm cm}+Q)}$ となる。$Q<0$ の場合、重心エネルギーは少なくとも $-Q$ でなければならない。

標的を実験室で静止させ、入射粒子の運動エネルギーを $K_a$ とする。非相対論では重心運動を引いた相対エネルギーは $E_{\rm cm}=\mu_{aA}v_a^2/2=m_AK_a/(m_a+m_A)$ である。厳密な式は全エネルギー運動量の不変量から得られる。$s=(E_a+m_Ac^2)^2-c^2p_a^2$ をエネルギー二乗の単位で定義すると
$$
s=(m_a+m_A)^2c^4+2m_Ac^2K_a,\qquad
K_{\rm th}=\max\left(0,\frac{[(m_b+m_B)^2-(m_a+m_A)^2]c^2}{2m_A}\right)
\label{eq-nuc-reaction-threshold}
$$
となる。出口二体の重心系での最小エネルギーは $\sqrt s=(m_b+m_B)c^2$ である。この条件を満たす最小の非負の $K_a$ が上式である。吸熱反応で $|Q|$ が静止エネルギーに比べて小さければ $K_{\rm th}\simeq -Q(1+m_a/m_A)$ となる。発熱反応の運動学的なしきい値は 0 だが、障壁で低エネルギーの率が抑えられる。

:::example[標的の反跳]
$m_a/m_A=1/4$、$Q=-2\,\mathrm{MeV}$ の反応では非相対論しきい値は $2.5\,\mathrm{MeV}$ である。実験室で不足する $0.5\,\mathrm{MeV}$ は最終状態全体の移動に使われる。
:::

:::example[欠損質量]
$b$ の全エネルギーと運動量を測れば、$E_B=E_a+m_Ac^2-E_b$、$\vb p_B=\vb p_a-\vb p_b$ と置ける。$M_Bc^2=\sqrt{E_B^2-c^2|\vb p_B|^2}$ から残留核の励起を求める。反応を起こした粒子以外も未検出なら、この量は未検出の全体系の不変質量を測る。
:::

### 二体の相対論的な運動量と逆運動学

$W=\sqrt s$ を全重心エネルギーとする。出口の二粒子の四元運動量を $p_b^\mu,p_B^\mu$、全運動量を $P^\mu$ として $(P-p_b)^2=m_B^2c^2$ を使う。重心系では $P^\mu=(W/c,\vb0)$ だから、$E_b^*=(W^2+m_b^2c^4-m_B^2c^4)/(2W)$ である。従って
$$
p_f^*c=\frac{\sqrt{[s-(m_b+m_B)^2c^4][s-(m_b-m_B)^2c^4]}}{2\sqrt s}
$$
となる。しきい値に近く $W=(m_b+m_B)c^2+\varepsilon$ なら、根号の第一因子が $2(m_b+m_B)c^2\varepsilon$、第二因子が $4m_bm_Bc^4$ から始まるので、$p_f^*\simeq\sqrt{2\mu_{bB}\varepsilon}$ に戻る。

重心系の速さを $V$ とし、実験室の入射方向を $z$ 軸に取る。非相対論の場合、粒子 $b$ の重心系での速さを $u$ とすると $\tan\theta_{\rm lab}=u\sin\theta_*/(V+u\cos\theta_*)$ である。ただし角度の象限は分子と分母の符号で定める。$V>u$ では縦速度が常に正なので全粒子が前方に出る。$\cos\theta_*=-u/V$ で角度が最大になり、$\sin\theta_{\rm lab,max}=u/V$ である。最大角より小さい角度に対応する重心系角が二つある場合、実験室の角断面積にはそれぞれの Jacobian を掛けて足す。

:::example[前方に集中する生成物]
$V=2u$ なら $\theta_{\rm lab,max}=30^\circ$ である。重心系で全方向に出る反応でも、運動する重心から観測すると広い角度まで粒子が届かない場合がある。
:::

非相対論の角度変換の係数を求める。$x=\cos\theta_*$、$w=V/u$ と書くと $\cos\theta_{\rm lab}=(x+w)/\sqrt{1+2wx+w^2}$ である。微分すると

$$
\frac{\dd\cos\theta_{\rm lab}}{\dd x}
=\frac{1+wx}{(1+2wx+w^2)^{3/2}} .
$$

軸周りの方位角は変わらないので、実験室系の微分断面積は、各解での重心系微分断面積に $|\dd x/\dd\cos\theta_{\rm lab}|$ を掛けた和になる。$V>u$ では $x=-1/w$ で変換の分母がゼロになる。実験でこの近傍を観測する場合は、点ごとの値に代えて有限の受容角について積分した計数を用いる。

### 二体運動学の実験室系への変換

非相対論的な二体の終状態で、重心の速度を $\vb V$、粒子 $b$ の重心系速度を $\vb u$ とする。実験室で $\vb v_b=\vb V+\vb u$ だから $T_b=T_0+T_*+2\sqrt{T_0T_*}\cos\theta_*$、ただし $T_0=m_bV^2/2,T_*=m_bu^2/2$ である。角度は $\tan\theta_{\rm lab}=u\sin\theta_* /(V+u\cos\theta_*)$ で、分母の符号まで含めて象限を選ぶ。重心で等方なら $\cos\theta_*$ は一様なので理想的な単一の終状態の $T_b$ は上下端の間で一様になる。実際のピーク幅にはビームの幅、標的での損失、検出器の受容角が関係する。

:::example[重心エネルギー]
静止標的に対して非相対論で、投射体と標的の質量比を $r=m_a/m_A$ とする。

計算には次の値を用いる。
$$
T_{a}=12\,\mathrm{MeV},\qquad r=\frac{1}{4}
$$

この条件で求める量を ${E_*}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{E_*}&=\frac{T_{a}}{r + 1}\\&=\frac{\left(12\right)}{\left(\frac{1}{4}\right) + 1}\\&\simeq 9.6\,\mathrm{MeV}.\end{aligned}
$$

逆運動学では同じ実験室エネルギーでも重心で使える比率が変わる。

質量を交換する。次の値を変え、ほかは固定する。
$$
r=4
$$

このときの値を ${E_*}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{E_*}_{\mathrm{new}}&=\frac{\left(12\right)}{\left(4\right) + 1}\\&\simeq 2.4\,\mathrm{MeV},\\{E_*}_{\mathrm{new}}-{E_*}&\simeq -7.2\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[吸熱反応の実験室しきい値]
負の $Q$ の絶対値を $q$ とする。非相対論近似で $T_{\rm th}=q(1+r)$ である。

計算には次の値を用いる。
$$
q=3\,\mathrm{MeV},\qquad r=\frac{1}{4}
$$

この条件で求める量を ${T_{\rm th}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{T_{\rm th}}&=q \cdot \left(r + 1\right)\\&=\left(3\right) \cdot \left(\left(\frac{1}{4}\right) + 1\right)\\&\simeq 3.75\,\mathrm{MeV}.\end{aligned}
$$

しきいでも終状態の重心に運動が残る。

重い投射体と比較する。次の値を変え、ほかは固定する。
$$
r=4
$$

このときの値を ${T_{\rm th}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{T_{\rm th}}_{\mathrm{new}}&=\left(3\right) \cdot \left(\left(4\right) + 1\right)\\&=15\,\mathrm{MeV},\\{T_{\rm th}}_{\mathrm{new}}-{T_{\rm th}}&\simeq 11.25\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[二体終状態での分配]
重心系で $b,B$ の質量比を $r=m_b/m_B$、相対運動の全エネルギーを $W$ とする。

計算には次の値を用いる。
$$
W=5\,\mathrm{MeV},\qquad r=\frac{1}{4}
$$

この条件で求める量を ${T_b^*}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{T_b^*}&=\frac{W}{r + 1}\\&=\frac{\left(5\right)}{\left(\frac{1}{4}\right) + 1}\\&=4\,\mathrm{MeV}.\end{aligned}
$$

同じ運動量をもつので軽い方が多くの運動エネルギーを受け取る。

同じ質量の二体を比べる。次の値を変え、ほかは固定する。
$$
r=1
$$

このときの値を ${T_b^*}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{T_b^*}_{\mathrm{new}}&=\frac{\left(5\right)}{\left(1\right) + 1}\\&\simeq 2.5\,\mathrm{MeV},\\{T_b^*}_{\mathrm{new}}-{T_b^*}&\simeq -1.5\,\mathrm{MeV}.\end{aligned}
$$
:::

:::example[実験室でのエネルギー幅]
重心が運動していて、前後の角度をすべて観測するとする。上端と下端の差を得る。

計算には次の値を用いる。
$$
T_{0}=\frac{1}{4}\,\mathrm{MeV},\qquad T_{s}=4\,\mathrm{MeV}
$$

この条件で求める量を ${\Delta T}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\Delta T}&=4 \cdot \sqrt{T_{0} \cdot T_{s}}\\&=4 \cdot \sqrt{\left(4\right) \cdot \left(\frac{1}{4}\right)}\\&=4\,\mathrm{MeV}.\end{aligned}
$$

この幅は寿命による固有幅とは区別する。

重心系の相対運動が小さい場合を比較する。次の値を変え、ほかは固定する。
$$
T_{s}=1\,\mathrm{MeV}
$$

このときの値を ${\Delta T}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\Delta T}_{\mathrm{new}}&=4 \cdot \sqrt{\left(1\right) \cdot \left(\frac{1}{4}\right)}\\&=2\,\mathrm{MeV},\\{\Delta T}_{\mathrm{new}}-{\Delta T}&=-2\,\mathrm{MeV}.\end{aligned}
$$
:::

## 断面積と保存則

入射流束を $\Phi$（単位面積・単位時間を横切る入射数）とする。標的一個当たりの事象率は $\Phi\sigma$ であり、:term[断面積]{en="cross section" yomi="だんめんせき"} $\sigma$ は面積の次元を持つ。粒子数 $N_t$ の薄い標的では事象率は $\Phi N_t\sigma$ である。角度を区別するなら $\dd\sigma/\dd\Omega$ を用い、有限角を覆う検出器ではその領域を積分し効率を掛ける。

チャネルは核種、内部状態、相対運動などで指定する。低エネルギーの強い相互作用・電磁相互作用では、エネルギー、運動量、電荷、核子数、全角運動量、パリティを保存する。入口のチャネルスピン $s$ は粒子と標的のスピンを結合したもの、全角運動量は $s$ と相対軌道 $\ell$ の合成である。全パリティは固有パリティの積と $(-1)^\ell$ の積になる。弱いチャネルではパリティをこの選別に使えない。時間反転対称な二体反応では、同じ $\sqrt s$ でスピン和した遷移確率が順逆で等しい。入口のスピン数 $g_a=(2J_a+1)(2J_A+1)$ と出口の $g_b$、相対波数 $k_a,k_b$ に対し、区別できる二体チャネルでは $g_a k_a^2\sigma_{ab}=g_b k_b^2\sigma_{ba}$ となる。ここでは入射状態平均と散乱の流束規格化を反転するときの違いを含めた。

:::example[薄い標的での確率]
面数密度 $10^{19}\,\mathrm{cm}^{-2}$ の薄膜に $\sigma=1\,\mathrm b=10^{-24}\,\mathrm{cm}^2$ の反応を起こす粒子を通す。未衝突の減衰を一次まで展開すると反応確率は $10^{-5}$ である。
:::

:::example[角運動量フィルター]
区別できる二つの $0^+$ 核が相対 $\ell=1$ で入ると $J^\pi=1^-$ である。二つの $0^+$ へ出るなら出口も $\ell=1$ が必要である。一方、片方が $0^-$ に変わった二つのスピンゼロ粒子では同じ $\ell$ によるパリティが正になり、この $J^\pi$ では結合できない。
:::

### 順逆反応の位相空間

同じ全重心エネルギーでの二体反応を考える。非相対論の平面波を体積 $\mathcal V$ で規格化し、時間反転対称な遷移振幅のスピン和を $W_{ab}$ とする。出口の相対状態密度は $\mathcal V p_b^2\dd p_b\dd\Omega/(2\pi\hbar)^3$ であり、エネルギー保存を適用すると $p_b^2/v_b$ に比例する。入射流束が $v_a/\mathcal V$ なので、断面積の運動学因子は $p_b^2/(v_a v_b)$ である。逆反応では入口・出口を入れ替えた $p_a^2/(v_bv_a)$ を得る。両方のスピン状態を足した $W_{ab}=W_{ba}$ から、入射だけを平均すると $g_a p_a^2\sigma_{ab}=g_b p_b^2\sigma_{ba}$ が従う。逆反応の比較に同じ実験室エネルギーを使うと、一般には全重心エネルギーが揃わない。

:::example[入口縮退度を含む比較]
区別可能な二粒子について、$g_a=2,g_b=4$、$p_b=p_a/2$ とする。同じ重心エネルギーで $\sigma_{ba}/\sigma_{ab}=g_ap_a^2/(g_bp_b^2)=2$ である。運動量だけを使うと 4 倍になり、入口スピン平均の違いを失う。
:::

### 流束で規格化した断面積

厚さ $d$ の標的内の標的核密度を $n$ とする。ビームが一つの過程で除去される場合、$\dd I/\dd x=-n\sigma I$ を積分して $I(d)=I(0)e^{-n\sigma d}$ を得る。薄い極限では反応確率が $n\sigma d$ になる。この展開は多重散乱を無視して断面積を求める根拠でもある。微分断面積を角度区間に積分して比較する場合、等方な係数の積分は $2\pi(\cos\theta_l-\cos\theta_h)$ を掛ける。角度が等間隔でも立体角の間隔は等しくない。

:::example[薄い標的の事象数]
一入射粒子当たりの衝突確率を面密度 $D$ と断面積 $\sigma$ の積とし、入射数 $N$ を掛ける。

計算には次の値を用いる。
$$
N=1000000,\qquad D=100000000000000000000\,\mathrm{cm}^{-2},\qquad \sigma=\frac{1}{2000000000000000000000000}\,\mathrm{cm^2}
$$

この条件で求める量を ${Y}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{Y}&=D \cdot N \cdot \sigma\\&=\left(100000000000000000000\right) \cdot \left(1000000\right) \cdot \left(\frac{1}{2000000000000000000000000}\right)\\&=50.\end{aligned}
$$

反応確率が十分小さければ投射体の減少を無視できる。

面密度を二倍にする。次の値を変え、ほかは固定する。
$$
D=200000000000000000000\,\mathrm{cm}^{-2}
$$

このときの値を ${Y}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{Y}_{\mathrm{new}}&=\left(1000000\right) \cdot \left(200000000000000000000\right) \cdot \left(\frac{1}{2000000000000000000000000}\right)\\&=100,\\{Y}_{\mathrm{new}}-{Y}&=50.\end{aligned}
$$
:::

:::example[角度受容の換算]
等方な微分断面積を $s$ とする。軸周りの $\theta_l,\theta_h$ の間で方位をすべて積分する。

計算には次の値を用いる。
$$
s=10\,\mathrm{mb/sr},\qquad t_{l}=0\,\mathrm{rad},\qquad t_{h}=\frac{\pi}{6}\,\mathrm{rad}
$$

この条件で求める量を ${\sigma}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\sigma}&=2 \cdot \pi \cdot s \cdot \left(- \cos{\left(t_{h} \right)} + \cos{\left(t_{l} \right)}\right)\\&=2 \cdot \pi \cdot \left(10\right) \cdot \left(\cos{\left(\left(0\right) \right)} - \cos{\left(\left(\frac{\pi}{6}\right) \right)}\right)\\&\simeq 8.41787\,\mathrm{mb}.\end{aligned}
$$

平面の角度間隔で単純に割ると受容に比べた量が変わる。

円錐を広げる。次の値を変え、ほかは固定する。
$$
t_{h}=\frac{\pi}{3}\,\mathrm{rad}
$$

このときの値を ${\sigma}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\sigma}_{\mathrm{new}}&=2 \cdot \pi \cdot \left(10\right) \cdot \left(\cos{\left(\left(0\right) \right)} - \cos{\left(\left(\frac{\pi}{3}\right) \right)}\right)\\&\simeq 31.4159\,\mathrm{mb},\\{\sigma}_{\mathrm{new}}-{\sigma}&\simeq 22.9981\,\mathrm{mb}.\end{aligned}
$$
:::

:::example[ルジャンドル成分の半球積分]
分布 $\dd\sigma/\dd\Omega=s(1+b\cos\theta)$ で前半球を積分する。

計算には次の値を用いる。
$$
s=10\,\mathrm{mb/sr},\qquad b=0
$$

この条件で求める量を ${\sigma_f}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\sigma_f}&=2 \cdot \pi \cdot s \cdot \left(\frac{b}{2} + 1\right)\\&=2 \cdot \pi \cdot \left(10\right) \cdot \left(\frac{\left(0\right)}{2} + 1\right)\\&\simeq 62.8319\,\mathrm{mb}.\end{aligned}
$$

全立体角では一次項が消え、前後の比には残る。

前方に強い分布と比べる。次の値を変え、ほかは固定する。
$$
b=\frac{3}{5}
$$

このときの値を ${\sigma_f}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\sigma_f}_{\mathrm{new}}&=2 \cdot \pi \cdot \left(10\right) \cdot \left(\frac{\left(\frac{3}{5}\right)}{2} + 1\right)\\&\simeq 81.6814\,\mathrm{mb},\\{\sigma_f}_{\mathrm{new}}-{\sigma_f}&\simeq 18.8496\,\mathrm{mb}.\end{aligned}
$$
:::

:::example[逆反応の運動量因子]
スピン縮退度が両チャネルで同一で、終状態の対称化に差がない場合を取る。時間反転関係で $k_i^2\sigma_{if}=k_f^2\sigma_{fi}$ となる。

計算には次の値を用いる。
$$
k_{i}=\frac{6}{5}\,\mathrm{fm}^{-1},\qquad k_{f}=1\,\mathrm{fm}^{-1}
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=\left(\frac{k_{i}}{k_{f}}\right)^{2}\\&=\frac{\left(\frac{6}{5}\right)^{2}}{\left(1\right)^{2}}\\&\simeq 1.44.\end{aligned}
$$

比較は同じ全重心エネルギーで行う。

終チャネルのしきいに近い場合を比べる。次の値を変え、ほかは固定する。
$$
k_{f}=\frac{1}{2}\,\mathrm{fm}^{-1}
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\frac{\left(\frac{6}{5}\right)^{2}}{\left(\frac{1}{2}\right)^{2}}\\&\simeq 5.76,\\{R}_{\mathrm{new}}-{R}&\simeq 4.32.\end{aligned}
$$
:::

## 光学模型

多くの出口へ逃げた流束を弾性チャネルから見れば、吸収として表せる。:term[光学模型]{en="optical model" yomi="こうがくもけい"} では複素ポテンシャル $U=V-iW$、$W\ge0$ を使う。Schrödinger 方程式に $\psi^*$ を掛け、共役式との差を取ると
$$
\partial_t|\psi|^2+\divg\vb j=-\frac{2W}{\hbar}|\psi|^2,\qquad
\vb j=\frac{\hbar}{2i\mu}(\psi^*\grad\psi-\psi\grad\psi^*)
$$
となる。虚部の符号によって流束が減る。チャネル空間を弾性 $P$ と残り $Q$ に射影すると、出口境界条件で $Q\psi=(E-H_{QQ}+i0)^{-1}H_{QP}P\psi$ である。入口式への代入で $H_{PQ}(E-H_{QQ}+i0)^{-1}H_{QP}$ が追加される。連続状態の分母の虚部は $-\pi\delta(E-H_{QQ})$ なので開いているチャネルが虚部を生む。実際の有効相互作用は一般に非局所でエネルギーに依存し、局所 $U$ はその近似である。

短距離で区別できるスピンなし粒子の球対称散乱の場合、散乱行列の部分波要素を $S_\ell$ とする。散乱振幅は $f(\theta)=\sum_\ell(2\ell+1)(S_\ell-1)P_\ell(\cos\theta)/(2ik)$ である。[部分波展開](/physics/advanced-quantum-mechanics/scattering/) の直交積分で弾性断面積を得る。また、各入射部分波の流束に対して出射の割合が $|S_\ell|^2$ であるため、
$$
\sigma_{\rm el}=\frac\pi{k^2}\sum_\ell(2\ell+1)|1-S_\ell|^2,\qquad
\sigma_{\rm reac}=\frac\pi{k^2}\sum_\ell(2\ell+1)(1-|S_\ell|^2).
\label{eq-nuc-optical-cross-sections}
$$
両式の和を整理すると $\sigma_{\rm tot}=(2\pi/k^2)\sum_\ell(2\ell+1)(1-\Re S_\ell)=4\pi\Im f(0)/k$ になる。これが短距離散乱の光学定理である。長距離 Coulomb の寄与を含めるときは、Coulomb 位相を $\sigma_\ell^C$ として $S_\ell=e^{2i\sigma_\ell^C}S_\ell^N$ と分ける。この分解を振幅に入れると、純 Coulomb 振幅に、$e^{2i\sigma_\ell^C}(S_\ell^N-1)$ を係数とする短距離の補正が加わる。微分弾性断面積ではこの二つの振幅を足してから二乗する。純 Coulomb の前方発散があるため全立体角の弾性積分は発散するが、反応断面積は $|S_\ell|^2=|S_\ell^N|^2$ を用いて求められる。

:::example[完全吸収する部分波]
一つの部分波で $S_\ell=0$ なら、式の弾性と反応の寄与はともに $\pi(2\ell+1)/k^2$ である。透過しない場合にも、入射波を妨げることで弾性の回折が残る。
:::

:::example[位相のみの散乱]
$S_0=-1$、他は $1$ なら、$|S_0|=1$ のため反応断面積はゼロである。弾性は $4\pi/k^2$ となる。同じ絶対値でも $S_0=1$ なら散乱は消える。
:::

### 虚部と流束保存

波動方程式に複素ポテンシャル $U=V-iW$ を入れ、方程式に $\psi^*$、共役方程式に $\psi$ を掛けて差を取ると $\partial_t|\psi|^2+\nabla\cdot\vb j=-2W|\psi|^2/\hbar$ である。$W\ge0$ は選んだ弾性空間からの流出を表す。定常入射で体積積分した消失率を流束で割れば反応断面積が得られる。部分波では入射と出射の流束比が $\abs{S_\ell}^2$ なので各波の消失率は $1-\abs{S_\ell}^2$ に対応する。分配先を個別に知るためにはチャネルに分けた模型を要する。

:::example[一つの部分波の吸収]
$S_\ell=\eta_\ell e^{2i\delta_\ell}$ とし、$k$ と部分波次数から寄与を得る。

計算には次の値を用いる。
$$
k=1\,\mathrm{fm}^{-1},\qquad l=2,\qquad \eta=\frac{9}{10}
$$

この条件で求める量を ${\sigma_\ell}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\sigma_\ell}&=\frac{\pi \cdot \left(1 - \eta^{2}\right) \cdot \left(2 \cdot l + 1\right)}{k^{2}}\\&=\frac{\pi \cdot \left(1 - \left(\frac{9}{10}\right)^{2}\right) \cdot \left(2 \cdot \left(2\right) + 1\right)}{\left(1\right)^{2}}\\&\simeq 2.98451\,\mathrm{fm^2}.\end{aligned}
$$

振幅の非弾性度の二乗を流束に使う。

より吸収の強い部分波と比べる。次の値を変え、ほかは固定する。
$$
\eta=\frac{1}{2}
$$

このときの値を ${\sigma_\ell}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\sigma_\ell}_{\mathrm{new}}&=\frac{\pi \cdot \left(1 - \left(\frac{1}{2}\right)^{2}\right) \cdot \left(2 \cdot \left(2\right) + 1\right)}{\left(1\right)^{2}}\\&\simeq 11.781\,\mathrm{fm^2},\\{\sigma_\ell}_{\mathrm{new}}-{\sigma_\ell}&\simeq 8.79646\,\mathrm{fm^2}.\end{aligned}
$$
:::

:::example[完全吸収模型の和]
$0\le\ell\le L$ で完全吸収、ほかでは消失ゼロとする。$\sum(2\ell+1)=(L+1)^2$ の和を使う。

計算には次の値を用いる。
$$
L=3,\qquad k=1\,\mathrm{fm}^{-1}
$$

この条件で求める量を ${\sigma_R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\sigma_R}&=\frac{\pi \cdot \left(L + 1\right)^{2}}{k^{2}}\\&=\frac{\pi \cdot \left(\left(3\right) + 1\right)^{2}}{\left(1\right)^{2}}\\&\simeq 50.2655\,\mathrm{fm^2}.\end{aligned}
$$

高い部分波が吸収されるなら幾何学的半径の推定も大きくなる。

吸収される波数範囲を広げる。次の値を変え、ほかは固定する。
$$
L=5
$$

このときの値を ${\sigma_R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\sigma_R}_{\mathrm{new}}&=\frac{\pi \cdot \left(\left(5\right) + 1\right)^{2}}{\left(1\right)^{2}}\\&\simeq 113.097\,\mathrm{fm^2},\\{\sigma_R}_{\mathrm{new}}-{\sigma_R}&\simeq 62.8319\,\mathrm{fm^2}.\end{aligned}
$$
:::

:::example[弾性部分波の強さ]
位相 $\delta$ と非弾性度 $\eta$ の波で $\abs{1-S}^2=1+\eta^2-2\eta\cos(2\delta)$ を求める。

計算には次の値を用いる。
$$
\eta=\frac{3}{5},\qquad \delta=0\,\mathrm{rad}
$$

この条件で求める量を ${u}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{u}&=\eta^{2} - 2 \cdot \eta \cdot \cos{\left(2 \cdot \delta \right)} + 1\\&=\left(\frac{3}{5}\right)^{2} - 2 \cdot \left(\frac{3}{5}\right) \cdot \cos{\left(2 \cdot \left(0\right) \right)} + 1\\&\simeq 0.16.\end{aligned}
$$

吸収が同じでも位相が異なると弾性断面積は変わる。

位相を変えた二模型を比べる。次の値を変え、ほかは固定する。
$$
\delta=\frac{\pi}{4}\,\mathrm{rad}
$$

このときの値を ${u}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{u}_{\mathrm{new}}&=\left(\frac{3}{5}\right)^{2} - 2 \cdot \left(\frac{3}{5}\right) \cdot \cos{\left(2 \cdot \left(\frac{\pi}{4}\right) \right)} + 1\\&\simeq 1.36,\\{u}_{\mathrm{new}}-{u}&\simeq 1.2.\end{aligned}
$$
:::

:::example[一定の虚ポテンシャルの時間依存]
運動を省いた有効模型で、選んだチャネルの虚エネルギーを $-iW$ とする。$u=Wt/\hbar$ での生存確率を求める。

計算には次の値を用いる。
$$
u=\frac{1}{5}
$$

この条件で求める量を ${P}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{P}&=e^{- 2 \cdot u}\\&=e^{- 2 \cdot \left(\frac{1}{5}\right)}\\&\simeq 0.67032.\end{aligned}
$$

振幅の減衰率と確率の減衰率は二倍違う。

長い相互作用時間と比べる。次の値を変え、ほかは固定する。
$$
u=\frac{4}{5}
$$

このときの値を ${P}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{P}_{\mathrm{new}}&=e^{- 2 \cdot \left(\frac{4}{5}\right)}\\&\simeq 0.201897,\\{P}_{\mathrm{new}}-{P}&\simeq -0.468424.\end{aligned}
$$
:::

## 複合核反応と Breit–Wigner 共鳴

入射で励起された中間核が入口の情報を多くの自由度に分け、その後に粒子や光子を放出する場合を :term[複合核反応]{en="compound-nucleus reaction" yomi="ふくごうかくはんのう"} と呼ぶ。単一の孤立準位の振幅は $\exp[-iE_R t/\hbar-\Gamma t/(2\hbar)]$ で減衰する。この振幅を $t\ge0$ で Fourier 変換すると分母 $E-E_R+i\Gamma/2$ が現れる。各出口の部分幅 $\Gamma_c=\hbar\lambda_c$ を用いると $\Gamma=\sum_c\Gamma_c$ である。入口の形成と出口への結合の振幅をそれぞれ $\sqrt{\Gamma_a},\sqrt{\Gamma_b}$ で表すと、孤立共鳴の寄与は
$$
\sigma_{a\to b}=\frac{\pi}{k_a^2}
\frac{2J+1}{(2J_a+1)(2J_A+1)}
\frac{\Gamma_a\Gamma_b}{(E-E_R)^2+\Gamma^2/4}.
\label{eq-nuc-breit-wigner}
$$
分子はチャネル結合の積、スピン係数は共鳴の射影数を入口の射影数で割ったものである。幅がエネルギーに依存するしきい値付近では、それを分母・分子に反映する。$b=a$ では弾性の非共鳴振幅との干渉があり、上式を単純に断面積に足せない場合がある。

孤立共鳴の規格化は流束保存から確かめられる。一つの $J^\pi$ と一つの固有準位を取り、各開チャネルを単位流束に規格化する。背景位相を取り去り、各チャネルへの結合を $g_c=\sqrt{\Gamma_c}$ と選ぶ。符号または位相はチャネル状態に吸収できる。$x=E-E_R$、$G_{ba}=g_b g_a$ とすれば、共鳴だけの散乱行列は

$$
S=\mathbf 1-\frac{iG}{x+i\Gamma/2}
$$

である。実際、$G^2=\Gamma G$ なので、$SS^\dagger$ の交差項の和は $-\Gamma G/(x^2+\Gamma^2/4)$、二次項はその逆符号となって $SS^\dagger=\mathbf 1$ を与える。異なるチャネル $b\ne a$ に出る割合は $|S_{ba}|^2$ であり、弾性では断面積に $|S_{aa}-1|^2$ が入る。これらに入口部分波の面積因子 $\pi/k_a^2$ と磁気射影の重みを掛けると :eqref[eq-nuc-breit-wigner] を得る。

中性粒子の短距離チャネルでは、しきい値から離れて一定とした部分幅が、出口直上では変化する。反応領域内の正則な終状態の部分波が $j_\ell(kr)\propto k^\ell$ で始まり、終状態密度をエネルギーで積分すると $k^2/(\dd E/\dd k)\propto k$ が残る。このため行列要素の先頭係数が有限で消えなければ $\Gamma_c\propto k^{2\ell+1}$ となる。例えば $p$ 波の出口で運動量が半分になると、部分幅の先頭項は $1/8$ になる。Coulomb 力がある場合や、しきい値の近くに別の極がある場合は、終状態の短距離振幅もその影響を受ける。

:::example[共鳴の強度]
$\Gamma_a=2\,\mathrm{keV}$、$\Gamma_b=3\,\mathrm{keV}$、$\Gamma=5\,\mathrm{keV}$ の場合、$E_R$ で式の最後の比は $4(2)(3)/25=0.96$ である。両幅が等しく全幅を使い切ると、この比は 1 になる。
:::

:::example[寿命と幅]
$\Gamma=1\,\mathrm{eV}$ の共鳴の平均寿命は $\hbar/\Gamma\simeq6.58\times10^{-16}\,\mathrm s$ である。非常に狭い共鳴なら入射のエネルギー分解能も観測形に効く。
:::

### 孤立共鳴の幅と面積

共鳴の状態を $e^{-iE_rt/\hbar-\Gamma t/(2\hbar)}$ で時間発展させる。$t\ge0$ で Fourier 積分すると $[\Gamma/2-i(E-E_r)]^{-1}$ に比例し、絶対値二乗は Lorentz 形になる。入口と出口の結合の二乗を部分幅 $\Gamma_i,\Gamma_f$ とすれば反応率の分子は $\Gamma_i\Gamma_f$ に比例する。すべての独立な出口への確率が足されるため $\Gamma=\sum_f\Gamma_f$ となる。$\Gamma$ が前因子の変化の尺度より十分小さい場合、エネルギー積分は分母に対する $\int_{-\infty}^\infty\dd E/[(E-E_r)^2+\Gamma^2/4]=2\pi/\Gamma$ を使える。しきい近傍の切断と幅のエネルギー依存が強いと、この面積の近似は成立しない。

:::example[中心の透過強度]
共鳴の入口・出口の幅 $g_i,g_f$ と全幅 $g$ を使い、幾何係数を除いた中心値を求める。

計算には次の値を用いる。
$$
g_{i}=1\,\mathrm{keV},\qquad g_{f}=1\,\mathrm{keV},\qquad g=4\,\mathrm{keV}
$$

この条件で求める量を ${F}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{F}&=\frac{4 \cdot g_{f} \cdot g_{i}}{g^{2}}\\&=\frac{4 \cdot \left(1\right)^{2}}{\left(4\right)^{2}}\\&\simeq 0.25.\end{aligned}
$$

全幅を固定した例では未観測の部分幅が補って変わる。

観測する出口への寄与を増す。次の値を変え、ほかは固定する。
$$
g_{f}=2\,\mathrm{keV}
$$

このときの値を ${F}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{F}_{\mathrm{new}}&=\frac{4 \cdot \left(1\right) \cdot \left(2\right)}{\left(4\right)^{2}}\\&\simeq 0.5,\\{F}_{\mathrm{new}}-{F}&\simeq 0.25.\end{aligned}
$$
:::

:::example[面積の規格化]
ピーク値 $\sigma_0$ と全幅 $g$ の単一線を積分する。

計算には次の値を用いる。
$$
\sigma_{0}=80\,\mathrm{mb},\qquad g=2\,\mathrm{keV}
$$

この条件で求める量を ${A}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{A}&=\frac{\pi \cdot g \cdot \sigma_{0}}{2}\\&=\frac{\pi \cdot \left(2\right) \cdot \left(80\right)}{2}\\&\simeq 251.327\,\mathrm{mb\,keV}.\end{aligned}
$$

ピークが同じでも積分値は線幅に応じて増える。

幅を二倍にした模型を比較する。次の値を変え、ほかは固定する。
$$
g=4\,\mathrm{keV}
$$

このときの値を ${A}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{A}_{\mathrm{new}}&=\frac{\pi \cdot \left(4\right) \cdot \left(80\right)}{2}\\&\simeq 502.655\,\mathrm{mb\,keV},\\{A}_{\mathrm{new}}-{A}&\simeq 251.327\,\mathrm{mb\,keV}.\end{aligned}
$$
:::

:::example[全幅と平均寿命]
分解された共鳴の全幅を MeV で扱う。

計算には次の値を用いる。
$$
\hbar=\frac{164553}{250000000000000000000000000}\,\mathrm{MeV\,s},\qquad g=\frac{1}{100}\,\mathrm{MeV}
$$

この条件で求める量を ${\tau}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\tau}&=\frac{\hbar}{g}\\&=\frac{\left(\frac{164553}{250000000000000000000000000}\right)}{\left(\frac{1}{100}\right)}\\&\simeq 6.58212\times10^{-20}\,\mathrm{s}.\end{aligned}
$$

$\Gamma_f$ だけからは全状態の寿命は決まらない。

幅が一桁大きい線と比べる。次の値を変え、ほかは固定する。
$$
g=\frac{1}{10}\,\mathrm{MeV}
$$

このときの値を ${\tau}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\tau}_{\mathrm{new}}&=\frac{\left(\frac{164553}{250000000000000000000000000}\right)}{\left(\frac{1}{10}\right)}\\&\simeq 6.58212\times10^{-21}\,\mathrm{s},\\{\tau}_{\mathrm{new}}-{\tau}&\simeq -5.92391\times10^{-20}\,\mathrm{s}.\end{aligned}
$$
:::

:::example[中心から外れた高さ]
離調 $d$ と全幅 $g$ で、中心に対する高さの比を求める。

計算には次の値を用いる。
$$
d=1\,\mathrm{keV},\qquad g=2\,\mathrm{keV}
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=1 \cdot \frac{1}{\frac{4 \cdot d^{2}}{g^{2}} + 1}\\&=\frac{1}{\frac{4 \cdot \left(1\right)^{2}}{\left(2\right)^{2}} + 1}\\&\simeq 0.5.\end{aligned}
$$

$d=g/2$ で高さが半分になる。

数線幅だけ離れた位置を見る。次の値を変え、ほかは固定する。
$$
d=3\,\mathrm{keV}
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\frac{1}{1 + \frac{4 \cdot \left(3\right)^{2}}{\left(2\right)^{2}}}\\&\simeq 0.1,\\{R}_{\mathrm{new}}-{R}&\simeq -0.4.\end{aligned}
$$

**一次近似の検証。** 中心線から離れた側で、エネルギー校正のずれが高さにどう伝わるか調べる。

${d}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${R}$ の一階微分を求めると、
$$
\frac{\mathrm d{R}}{\mathrm dd}=- \frac{8 \cdot d \cdot g^{2}}{\left(4 \cdot d^{2} + g^{2}\right)^{2}}.
$$

この導関数を変更前の点で評価して $h=2$ を掛ける。ここで $d$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&=-1,\\\Delta-\Delta_{\mathrm{lin}}&\simeq 0.6.\end{aligned}
$$

この変化では一次近似の誤差が大きいため、増分の計算に元の式を使う。
:::

## Hauser–Feshbach 理論の概観

多くの複合核準位をエネルギー平均すると、各準位の波形を追う代わりに、形成の透過係数と平均の崩壊分岐で記述できる。ただし入口で固定した $J^\pi$ は出口でも守る。透過係数 $T_c$ はチャネルの開き方を表す無次元量である。幅の相関を省く理想化では、形成断面積 $\pi g_J T_a/k_a^2$ に、出口の相対重み $T_b/\sum_c T_c$ を掛けて
$$
\sigma_{ab}^{\rm HF}=\frac{\pi}{k_a^2}\sum_{J,\pi}
g_J\,\frac{T_a^{J\pi}T_b^{J\pi}}{\sum_c T_c^{J\pi}}
\label{eq-nuc-statistical-reaction}
$$
となる。$g_J=(2J+1)/[(2J_a+1)(2J_A+1)]$ であり、複数の相対軌道やチャネルスピンはそれぞれのチャネルに含める。すべての出口について足すと分母と和が相殺し、形成断面積に戻る。幅揺らぎの補正には散乱行列の統計分布が必要であり、単一準位が支配する場合には :eqref[eq-nuc-breit-wigner] を使う。

:::example[競合出口]
透過係数が $T_a=0.2,T_b=0.3,T_c=0.5$ の一つの $J^\pi$ では、形成後に $b$ へ出る重みは $0.3$ である。断面積は $\pi g_J/k_a^2$ の $0.06$ 倍になる。
:::

:::example[開くチャネル]
同じ例で $T_c$ がゼロなら $b$ への重みは $0.3/(0.2+0.3)=0.6$ になる。他の係数が固定なら、競合チャネルが閉じると同じ出口の断面積は 2 倍になる。
:::

### 狭い共鳴の列と統計模型のつながり

同じ $J^\pi$ の平均準位間隔を $D$ とする。各幅が $D$ より小さく、多数の共鳴を含む区間で $k_a$ と幅分布がほぼ一定と仮定する。:eqref[eq-nuc-breit-wigner] を一準位について積分すると、断面積の面積は $\pi g_J/k_a^2$ に $2\pi\Gamma_a\Gamma_b/\Gamma$ を掛けた値である。区間幅で割ると、係数は $(2\pi/D)\langle\Gamma_a\Gamma_b/\Gamma\rangle$ となる。ここで $\langle\cdot\rangle$ は準位にわたる平均である。

幅の相関を無視できる近似では $\langle\Gamma_a\Gamma_b/\Gamma\rangle\simeq\langle\Gamma_a\rangle\langle\Gamma_b\rangle/\sum_c\langle\Gamma_c\rangle$ と扱う。この極限で $T_c\simeq2\pi\langle\Gamma_c\rangle/D$ を代入すると :eqref[eq-nuc-statistical-reaction] の形になる。同じ出口と入口を選ぶと $\Gamma_a^2$ が現れ、平均の二乗と二乗の平均が異なる。その差は幅の揺らぎを省く近似が破れる一例である。

平均散乱行列がチャネルについて対角で、直接過程との混合を無視できるとする。そのとき透過係数の定義は $T_c=1-|\langle S_{cc}\rangle|^2$ である。ここでは振幅をエネルギー平均した後に二乗する。先に二乗した $\langle|S_{cc}|^2\rangle$ には、複合核を経て同じチャネルに戻った揺らぎの寄与も入る。複合核の形成に対応する流束を求めるときに前者を使うのは、この戻りを平均弾性波と区別するためである。

:::example[揺らぎで変わる積の平均]
二準位であるチャネルの幅が $1$ と $3\,\mathrm{keV}$ なら、平均幅は $2\,\mathrm{keV}$、二乗の平均は $5\,\mathrm{keV}^2$ である。平均幅を二乗すると $4\,\mathrm{keV}^2$ であり、幅の分散が $1\,\mathrm{keV}^2$ の差を生む。
:::

### 透過係数による統計的な分岐

固定スピン・パリティの共鳴群で、入口から複合核が形成され、記憶を失って出口へ崩壊する近似を考える。幅揺らぎを無視した平均で出口の比は $T_f/\sum_cT_c$、入口の形成に掛かる因子は $T_i$ である。したがって同じ入口・同じスピンで全出口を足すと、分母が相殺され形成の因子に戻る。少数の準位しか平均できない場合や入口と出口が同じ場合には、部分幅の相関を無視した比が不正確になりうる。異なる $J^\pi$ で透過の和を先に作ると分母の違いを失うため、各組ごとに分岐を求めてから足す。

:::example[二つの出口の競合]
透過重みを $T_f,T_g$ として分岐比を求める。

計算には次の値を用いる。
$$
T_{f}=\frac{1}{5},\qquad T_{g}=\frac{1}{10}
$$

この条件で求める量を ${b_f}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{b_f}&=\frac{T_{f}}{T_{f} + T_{g}}\\&=\frac{\left(\frac{1}{5}\right)}{\left(\frac{1}{10}\right) + \left(\frac{1}{5}\right)}\\&\simeq 0.666667.\end{aligned}
$$

形成率が一定でも新しい出口は個別の分岐を変える。

競合する透過の大きい場合と比較する。次の値を変え、ほかは固定する。
$$
T_{g}=\frac{3}{5}
$$

このときの値を ${b_f}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{b_f}_{\mathrm{new}}&=\frac{\left(\frac{1}{5}\right)}{\left(\frac{1}{5}\right) + \left(\frac{3}{5}\right)}\\&\simeq 0.25,\\{b_f}_{\mathrm{new}}-{b_f}&\simeq -0.416667.\end{aligned}
$$
:::

:::example[入口と出口の積]
同じ $J^\pi$ で入口の形成係数 $T_i$ を掛けた出口確率を求める。分母は開いている全てのチャネルの合計 $T$ である。

計算には次の値を用いる。
$$
T_{i}=\frac{2}{5},\qquad T_{f}=\frac{1}{10},\qquad T=\frac{6}{5}
$$

この条件で求める量を ${P}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{P}&=\frac{T_{f} \cdot T_{i}}{T}\\&=\frac{\left(\frac{1}{10}\right) \cdot \left(\frac{2}{5}\right)}{\left(\frac{6}{5}\right)}\\&\simeq 0.0333333.\end{aligned}
$$

全透過を固定して入口を変えるなら他の成分も補っている。

入口の結合が強い場合を比べる。次の値を変え、ほかは固定する。
$$
T_{i}=\frac{3}{5}
$$

このときの値を ${P}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{P}_{\mathrm{new}}&=\frac{\left(\frac{1}{10}\right) \cdot \left(\frac{3}{5}\right)}{\left(\frac{6}{5}\right)}\\&\simeq 0.05,\\{P}_{\mathrm{new}}-{P}&\simeq 0.0166667.\end{aligned}
$$
:::

:::example[出口の連続近似]
区間幅 $D$ に単位エネルギー当たり $n$ の終準位があり、それぞれの透過の平均を $t$ とする。

計算には次の値を用いる。
$$
n=20\,\mathrm{MeV}^{-1},\qquad D=\frac{1}{2}\,\mathrm{MeV},\qquad t=\frac{1}{50}
$$

この条件で求める量を ${T_{\rm sum}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{T_{\rm sum}}&=D \cdot n \cdot t\\&=\left(20\right) \cdot \left(\frac{1}{2}\right) \cdot \left(\frac{1}{50}\right)\\&\simeq 0.2.\end{aligned}
$$

これは個々の透過が小さい場合でも最終状態が多いと総寄与が増す例である。

終状態の窓を広げる。次の値を変え、ほかは固定する。
$$
D=1\,\mathrm{MeV}
$$

このときの値を ${T_{\rm sum}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{T_{\rm sum}}_{\mathrm{new}}&=\left(1\right) \cdot \left(20\right) \cdot \left(\frac{1}{50}\right)\\&\simeq 0.4,\\{T_{\rm sum}}_{\mathrm{new}}-{T_{\rm sum}}&\simeq 0.2.\end{aligned}
$$
:::

:::example[複数のスピンからの分岐]
二つのスピン形成割合を $w,1-w$、それぞれの出口分岐を $b_1,b_2$ とする。

計算には次の値を用いる。
$$
w=\frac{3}{10},\qquad b_{1}=\frac{1}{5},\qquad b_{2}=\frac{7}{10}
$$

この条件で求める量を ${b}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{b}&=b_{1} \cdot w + b_{2} \cdot \left(1 - w\right)\\&=\left(\frac{1}{5}\right) \cdot \left(\frac{3}{10}\right) + \left(\frac{7}{10}\right) \cdot \left(1 - \left(\frac{3}{10}\right)\right)\\&\simeq 0.55.\end{aligned}
$$

重みは形成時のスピン分布から計算する。

入口でのスピン分布が変わる場合を比べる。次の値を変え、ほかは固定する。
$$
w=\frac{3}{5}
$$

このときの値を ${b}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{b}_{\mathrm{new}}&=\left(\frac{1}{5}\right) \cdot \left(\frac{3}{5}\right) + \left(\frac{7}{10}\right) \cdot \left(1 - \left(\frac{3}{5}\right)\right)\\&\simeq 0.4,\\{b}_{\mathrm{new}}-{b}&\simeq -0.15.\end{aligned}
$$
:::

## 直接反応（ストリッピングとピックアップ）

入射粒子が核の表面を短時間で通り、少ない自由度に運動量や核子を渡す過程を :term[直接反応]{en="direct reaction" yomi="ちょくせつはんのう"} という。$(d,p)$ で入射重陽子から中性子を標的に移すのがストリッピングの例であり、逆の $(p,d)$ で標的から中性子を取り出すのがピックアップの例である。出口の角分布には移した角運動量の情報が残る。

入口と出口の波を光学模型で歪めた $\chi_i^{(+)},\chi_f^{(-)}$ とし、内部状態を $\phi_i,\phi_f$ とする。Hamiltonian をチャネルの基準と差 $V_{\rm tr}$ に分け、転移振幅の一次だけを取ると
$$
T_{fi}^{\rm DWBA}=\bra{\chi_f^{(-)}\phi_f}V_{\rm tr}\ket{\chi_i^{(+)}\phi_i}
$$
である。これを歪曲波 Born 近似という。平面波の極限では運動量移行に対応する位相 $e^{i\vb q\cdot\vb r}$ が残る。これを球面波で展開すると移した軌道の角運動量 $\ell$ に応じて球 Bessel 関数 $j_\ell(qr)$ が現れ、低 $q$ で $(qr)^\ell/(2\ell+1)!!$ となる。実験との強度比較には束縛波の重なり、光学ポテンシャル、反応の他の経路が効く。

:::example[核子数の収支]
${}^{16}\mathrm O(d,p){}^{17}\mathrm O$ では入口の核子数が $16+2=18$、出口は $17+1=18$ である。陽子数も $8+1=8+1$ であり、標的が受け取ったのは中性子である。
:::

:::example[小さい運動量移行]
平面波近似で $s$ 波と $d$ 波の動径積分が規則的なら、振幅はそれぞれ $q^0$ と $q^2$ で始まる。$d$ 波の強度は $q^4$ で立ち上がる。ただし実際の角分布は入口と出口の歪曲を含めて比較する。
:::

### 移行の重なりと選択性

一段階近似の振幅は初期と終状態の歪曲波を $\chi_i,\chi_f$、移行相互作用を $V$ として $\langle\chi_f\phi_f|V|\chi_i\phi_i\rangle$ である。平面波・短距離結合まで単純化すると、運動量移行の Fourier 変換を測る構造になる。角運動量 $L$ の成分では $j_L(qr)$ が現れ、$q\to0$ で振幅は $q^L$ に比例する。ただし動径モーメントが消えないと仮定する。このため小角領域の強度は移行角運動量の判定に役立つ。一方、絶対値から配置の重みを取り出すには、光学模型、反応演算子、重なりの規格化もそろえる。

:::example[運動量移行の最小値]
入出の波数を $k_i,k_f$、なす角を $\theta$ とする。

計算には次の値を用いる。
$$
k_{i}=\frac{7}{5}\,\mathrm{fm}^{-1},\qquad k_{f}=\frac{6}{5}\,\mathrm{fm}^{-1},\qquad \theta=0\,\mathrm{rad}
$$

この条件で求める量を ${q}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{q}&=\sqrt{k_{f}^{2} - 2 \cdot k_{f} \cdot k_{i} \cdot \cos{\left(\theta \right)} + k_{i}^{2}}\\&=\sqrt{\left(\frac{6}{5}\right)^{2} - 2 \cdot \left(\frac{6}{5}\right) \cdot \left(\frac{7}{5}\right) \cdot \cos{\left(\left(0\right) \right)} + \left(\frac{7}{5}\right)^{2}}\\&\simeq 0.2\,\mathrm{fm}^{-1}.\end{aligned}
$$

前方でも波数が違えば $q=0$ にはならない。

散乱角を増やす。次の値を変え、ほかは固定する。
$$
\theta=\frac{\pi}{6}\,\mathrm{rad}
$$

このときの値を ${q}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{q}_{\mathrm{new}}&=\sqrt{\left(\frac{6}{5}\right)^{2} - 2 \cdot \left(\frac{6}{5}\right) \cdot \left(\frac{7}{5}\right) \cdot \cos{\left(\left(\frac{\pi}{6}\right) \right)} + \left(\frac{7}{5}\right)^{2}}\\&\simeq 0.70011\,\mathrm{fm}^{-1},\\{q}_{\mathrm{new}}-{q}&\simeq 0.50011\,\mathrm{fm}^{-1}.\end{aligned}
$$
:::

:::example[小運動量の軌道選択性]
規格化した模型で $x=qR$、同じ動径係数を仮定して振幅の最初の次数の二乗を比較する。

計算には次の値を用いる。
$$
x=\frac{1}{4},\qquad L=1
$$

この条件で求める量を ${R}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{R}&=x^{2 \cdot L}\\&=\left(\frac{1}{4}\right)^{2 \cdot \left(1\right)}\\&\simeq 0.0625.\end{aligned}
$$

一般には $L$ ごとに係数や動径積分も異なる。

高い階数での抑制を比べる。次の値を変え、ほかは固定する。
$$
L=2
$$

このときの値を ${R}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{R}_{\mathrm{new}}&=\left(\frac{1}{4}\right)^{2 \cdot \left(2\right)}\\&\simeq 0.00390625,\\{R}_{\mathrm{new}}-{R}&\simeq -0.0585938.\end{aligned}
$$
:::

:::example[移行強度の規格化]
同じ角度で単一粒子模型の断面積 $s_0$ とデータ $s$ を比べる。一段階で比例関係を仮定する。

計算には次の値を用いる。
$$
s=\frac{4}{5}\,\mathrm{mb/sr},\qquad s_{0}=\frac{6}{5}\,\mathrm{mb/sr}
$$

この条件で求める量を ${C^2S}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{C^2S}&=\frac{s}{s_{0}}\\&=\frac{\left(\frac{4}{5}\right)}{\left(\frac{6}{5}\right)}\\&\simeq 0.666667.\end{aligned}
$$

推定値は反応模型の断面積の規格化に依存する。

歪曲波模型の違いを模擬する。次の値を変え、ほかは固定する。
$$
s_{0}=\frac{8}{5}\,\mathrm{mb/sr}
$$

このときの値を ${C^2S}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{C^2S}_{\mathrm{new}}&=\frac{\left(\frac{4}{5}\right)}{\left(\frac{8}{5}\right)}\\&\simeq 0.5,\\{C^2S}_{\mathrm{new}}-{C^2S}&\simeq -0.166667.\end{aligned}
$$
:::

:::example[有限分解能での二準位]
二つの直交した終状態の角度積分断面積を $s_1,s_2$ とし、窓内にそれぞれ $f_1,f_2$ が入る。

計算には次の値を用いる。
$$
f_{1}=1,\qquad f_{2}=\frac{1}{10},\qquad s_{1}=2\,\mathrm{mb},\qquad s_{2}=3\,\mathrm{mb}
$$

この条件で求める量を ${\sigma_{\rm window}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\sigma_{\rm window}}&=f_{1} \cdot s_{1} + f_{2} \cdot s_{2}\\&=\left(1\right) \cdot \left(2\right) + \left(3\right) \cdot \left(\frac{1}{10}\right)\\&\simeq 2.3\,\mathrm{mb}.\end{aligned}
$$

直交した終状態の未分離ピークでは確率が加わる。

窓を広げる。次の値を変え、ほかは固定する。
$$
f_{2}=\frac{3}{5}
$$

このときの値を ${\sigma_{\rm window}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\sigma_{\rm window}}_{\mathrm{new}}&=\left(1\right) \cdot \left(2\right) + \left(3\right) \cdot \left(\frac{3}{5}\right)\\&\simeq 3.8\,\mathrm{mb},\\{\sigma_{\rm window}}_{\mathrm{new}}-{\sigma_{\rm window}}&\simeq 1.5\,\mathrm{mb}.\end{aligned}
$$
:::

### 遷移演算子の prior 形式と post 形式

Hermite なチャネル Hamiltonian で、入口と出口をそれぞれ $H_i,H_f$ に分ける。内部状態と相対運動を含めたチャネル波を $\Phi_i^{(+)},\Phi_f^{(-)}$ とする。全エネルギーが同じ $E$ の殻上では、$\Phi_i^{(+)}$ に右から $H_i$ を作用させても $\Phi_f^{(-)}$ に左から $H_f$ を作用させても、同じ $E\langle\Phi_f^{(-)}|\Phi_i^{(+)}\rangle$ が得られる。この恒等式は散乱状態を波束で正則化してエネルギーを殻上に取った場合にも成り立つ。従って
$$
\bra{\Phi_f^{(-)}}(H-H_i)\ket{\Phi_i^{(+)}}
-
\bra{\Phi_f^{(-)}}(H-H_f)\ket{\Phi_i^{(+)}}=0
$$
である。前者を :term[prior 形式]{en="prior form" yomi="ぷらいあーけいしき"}、後者を :term[post 形式]{en="post form" yomi="ぽすとけいしき"} と呼ぶ。実際の近似で複素光学場を使う場合は、左の波に随伴方程式を使う必要がある。また負担を減らすために、両方から同じ意味の項を落とせるとは限らない。例えば移行する核子と芯の有限範囲結合を、一方だけで零距離化すればこの差が数値上は残りうる。二形式を比較する場合は境界条件、内部波の規格化、相互作用の分割をそろえる。

## 重イオン反応

核同士が接近すると、遠方での Coulomb 反発と、近距離での核力が競合する。正面に近い衝突、周辺での核子移行、非弾性散乱、融合、破砕などがエネルギーと角運動量に応じて競合する。周辺衝突では二核がそれぞれの形を保つ時間が長く、変形や振動への結合が散乱に現れる。多チャネル結合を無視して剛体球で障壁を表すと、表面振動で複数の障壁に分かれる効果を失う。

古典的な衝突径数 $b$ で完全吸収を $b<R$ に限る模型では、断面積は円環 $2\pi b\dd b$ の積分 $\pi R^2$ となる。量子の :eqref[eq-nuc-optical-cross-sections] で $\ell\le \ell_{\max}$ だけが完全に吸収されると、和は $\pi(\ell_{\max}+1)^2/k^2$ になる。$\ell_{\max}\sim kR$ が大きい極限で古典値に近づく。この模型での $R$ はチャネルの反応が起こる距離を表し、単独の核の電荷半径とは区別する。

:::example[部分波の和]
$\ell=0,1,2,3$ で完全吸収なら $1+3+5+7=16$ で、$\sigma_{\rm reac}=16\pi/k^2$ である。最も高い一つの波だけを足すと、この和を過小評価する。
:::

:::example[幾何学的な尺度]
完全吸収の模型で $R=6\,\mathrm{fm}$ なら $\pi R^2\simeq113\,\mathrm{fm}^2\simeq1.13\,\mathrm b$ となる。低エネルギーでは Coulomb 反発と量子透過でこの見積もりから離れる。
:::

### 直接反応と複合核を見分ける測定

出口の粒子エネルギーと角分布を同時に調べると、同じ生成物へ至る複数の機構を区別しやすい。移行反応で残留核の準位が離れている場合は、欠損質量で各準位を選び、それぞれの角分布を比較する。反応時間が通過時間と同じ桁なら、入口の運動量が出口に強く残る。一方、複合核が多くの自由度へエネルギーを配った後の分岐では、多数の準位が重なって同じ生成物を与え得る。角分布の等方性だけで機構は決まらず、複合核にも全角運動量の配向による異方性が残り得る。

:::example[時間幅の比較]
反応領域を $R=5\,\mathrm{fm}$、粒子の相対速度を $c/10$ とすれば、通過時間 $R/v$ は約 $1.7\times10^{-22}\,\mathrm s$ である。幅 $1\,\mathrm{keV}$ の中間状態の寿命は $\hbar/\Gamma\simeq6.6\times10^{-19}\,\mathrm s$ であり、通過時間より長い。時間の差は二つの機構の模型化を分ける手掛かりになる。
:::

### 重イオンの障壁分布と角運動量

古典的に衝突径数を $b$ とすると重心での軌道角運動量は $pb$ である。引力の核ポテンシャルと Coulomb 反発を半径 $R_b$ の高さ $V_b$ の障壁で近似する。$E=V_b+E b^2/R_b^2$ が古典的に越えられる最大 $b$ を与えるので、$E>V_b$ なら融合断面積の簡単な模型は $\pi R_b^2(1-V_b/E)$ となる。量子模型では下側でもトンネルがあり、また変形や励起への結合で障壁が分布する。1 本の障壁の高さを当てはめる場合は、エネルギー帯を変えたときの当てはめの安定性を確認する。

障壁頂上の近傍を逆放物線で近似し、曲率をエネルギー単位で $W=\hbar\omega_b$ と書く。中に入った流束はすべて吸収されるとする。部分波 $\ell$ の障壁を $V_b+\hbar^2\ell(\ell+1)/(2\mu R_b^2)$ とし、半径と曲率は共通とみなす。Hill–Wheeler の透過率を使うと、部分波和の各項は $(2\ell+1)/[1+\exp((2\pi/W)(V_b-E+\hbar^2\ell(\ell+1)/(2\mu R_b^2)))]$ となる。十分多くの $\ell$ が寄与する場合は和を積分に置き換え、$y=\hbar^2\ell(\ell+1)/(2\mu R_b^2)$ とする。すると $(2\ell+1)\dd\ell=(2\mu R_b^2/\hbar^2)\dd y$ である。前因子 $\pi/k^2$ と合わせ、残りの積分を実行すると

$$
\sigma_{\rm fus}(E)\simeq\frac{W R_b^2}{2E}
\ln\!\left[1+\exp\!\left(\frac{2\pi(E-V_b)}{W}\right)\right].
$$

障壁より十分上では対数の先頭が $2\pi(E-V_b)/W$ なので、古典式に戻る。逆に障壁より深く下ではポテンシャルの非放物線性が効き、寄与する部分波が少ないと積分近似も崩れる。

:::example[古典的な捕獲面積]
障壁より上で、$R_b,V_b,E$ から面積を計算する。

計算には次の値を用いる。
$$
R_{b}=6\,\mathrm{fm},\qquad V_{b}=15\,\mathrm{MeV},\qquad E=16\,\mathrm{MeV}
$$

この条件で求める量を ${\sigma}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\sigma}&=\pi \cdot R_{b}^{2} \cdot \left(1 - \frac{V_{b}}{E}\right)\\&=\pi \cdot \left(6\right)^{2} \cdot \left(- \frac{\left(15\right)}{\left(16\right)} + 1\right)\\&\simeq 7.06858\,\mathrm{fm^2}.\end{aligned}
$$

しきいより下にこの式を延長すると負になるため、古典式はそこでゼロに切る。

しきいから離れて比べる。次の値を変え、ほかは固定する。
$$
E=24\,\mathrm{MeV}
$$

このときの値を ${\sigma}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\sigma}_{\mathrm{new}}&=\pi \cdot \left(6\right)^{2} \cdot \left(- \frac{\left(15\right)}{\left(24\right)} + 1\right)\\&\simeq 42.4115\,\mathrm{fm^2},\\{\sigma}_{\mathrm{new}}-{\sigma}&\simeq 35.3429\,\mathrm{fm^2}.\end{aligned}
$$
:::

:::example[障壁の周りの透過]
逆放物線の Hill–Wheeler 透過 $P=[1+\exp(2\pi(V_b-E)/W)]^{-1}$ を用いる。

計算には次の値を用いる。
$$
V_{b}=15\,\mathrm{MeV},\qquad E=\frac{29}{2}\,\mathrm{MeV},\qquad W=3\,\mathrm{MeV}
$$

この条件で求める量を ${P}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{P}&=1 \cdot \frac{1}{e^{\frac{2 \cdot \pi \cdot \left(- E + V_{b}\right)}{W}} + 1}\\&=\frac{1}{e^{\frac{2 \cdot \pi \cdot \left(\left(15\right) - \left(\frac{29}{2}\right)\right)}{\left(3\right)}} + 1}\\&\simeq 0.259764.\end{aligned}
$$

$W$ は頂上の曲率をエネルギーで表した量である。

中心の上下で比較する。次の値を変え、ほかは固定する。
$$
E=\frac{31}{2}\,\mathrm{MeV}
$$

このときの値を ${P}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{P}_{\mathrm{new}}&=\frac{1}{e^{\frac{2 \cdot \pi \cdot \left(\left(15\right) - \left(\frac{31}{2}\right)\right)}{\left(3\right)}} + 1}\\&\simeq 0.740236,\\{P}_{\mathrm{new}}-{P}&\simeq 0.480473.\end{aligned}
$$
:::

:::example[エネルギーと進行方向]
等方の二体重心エネルギーで、重心の速さと終粒子の速さの比を $r=V/u$ とする。$\theta_*=\pi/2$ の実験室の角度を求める。

計算には次の値を用いる。
$$
r=\frac{1}{2}
$$

この条件で求める量を ${\theta_{\rm lab}}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{\theta_{\rm lab}}&=\arctan{\left(1 \cdot \frac{1}{r} \right)}\\&=\arctan{\left(\frac{1}{\left(\frac{1}{2}\right)} \right)}\\&\simeq 1.10715\,\mathrm{rad}.\end{aligned}
$$

逆運動学の大きな重心速度では横向きの放出も前方へ偏る。

重心がより速い場合を比較する。次の値を変え、ほかは固定する。
$$
r=2
$$

このときの値を ${\theta_{\rm lab}}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{\theta_{\rm lab}}_{\mathrm{new}}&=\arctan{\left(\frac{1}{\left(2\right)} \right)}\\&\simeq 0.463648\,\mathrm{rad},\\{\theta_{\rm lab}}_{\mathrm{new}}-{\theta_{\rm lab}}&\simeq -0.643501\,\mathrm{rad}.\end{aligned}
$$
:::

:::example[二障壁の平均透過]
互いに分離した二つの入口配置で高さ $V_1,V_2$、共通の曲率 $W$ とする。混合割合を $w$ として透過率を平均する。

計算には次の値を用いる。
$$
V_{1}=14\,\mathrm{MeV},\qquad V_{2}=16\,\mathrm{MeV},\qquad E=15\,\mathrm{MeV},\qquad W=3\,\mathrm{MeV},\qquad w=\frac{1}{2}
$$

この条件で求める量を ${P}$ とおく。式に数値を代入すると、
$$
\begin{aligned}{P}&=\frac{w}{e^{\frac{2 \cdot \pi \cdot \left(- E + V_{1}\right)}{W}} + 1} + \frac{1 - w}{e^{\frac{2 \cdot \pi \cdot \left(- E + V_{2}\right)}{W}} + 1}\\&=\frac{\left(\frac{1}{2}\right)}{e^{\frac{2 \cdot \pi \cdot \left(\left(14\right) - \left(15\right)\right)}{\left(3\right)}} + 1} + \frac{1 - \left(\frac{1}{2}\right)}{e^{\frac{2 \cdot \pi \cdot \left(- \left(15\right) + \left(16\right)\right)}{\left(3\right)}} + 1}\\&\simeq 0.5.\end{aligned}
$$

低い障壁の寄与が増すと平均率が上がる。本来の結合チャネルでは入口からの振幅を解く必要がある。

低い障壁の重みを増す。次の値を変え、ほかは固定する。
$$
w=\frac{4}{5}
$$

このときの値を ${P}_{\mathrm{new}}$ と書く。同じ式で計算すると、
$$
\begin{aligned}{P}_{\mathrm{new}}&=\frac{\left(\frac{4}{5}\right)}{e^{\frac{2 \cdot \pi \cdot \left(\left(14\right) - \left(15\right)\right)}{\left(3\right)}} + 1} + \frac{1 - \left(\frac{4}{5}\right)}{e^{\frac{2 \cdot \pi \cdot \left(- \left(15\right) + \left(16\right)\right)}{\left(3\right)}} + 1}\\&\simeq 0.734214,\\{P}_{\mathrm{new}}-{P}&\simeq 0.234214.\end{aligned}
$$

**一次近似の検証。** 各入口の高さと曲率が固定なら、その重みの不確かさは透過の線形誤差になる。

${w}_0$ を変更前の値、$h$ をその入力の増分とする。ほかの入力を固定し、${P}$ の一階微分を求めると、
$$
\frac{\mathrm d{P}}{\mathrm dw}=- \frac{1}{1 + e^{- \frac{2 \cdot \pi \cdot \left(E - V_{2}\right)}{W}}} + \frac{1}{1 + e^{- \frac{2 \cdot \pi \cdot \left(E - V_{1}\right)}{W}}}.
$$

この導関数を変更前の点で評価して $h=\frac{3}{10}$ を掛ける。ここで $w$ の増分は上で用いた単位のまま代入する。変更後の値から変更前の値を引いた厳密な差を $\Delta$ とすると、数値は次のように比較できる。
$$
\begin{aligned}\Delta_{\mathrm{lin}}&\simeq 0.234214,\\\Delta-\Delta_{\mathrm{lin}}&=0.\end{aligned}
$$

残差が示すのは、指定した二点間で生じる線形化の誤差である。基準点からさらに離れた値について同じ精度を保証するものではない。
:::

## 演習問題

::::exercise{level=basic}
$Q=-4\,\mathrm{MeV}$、$m_A=3m_a$ の吸熱反応の非相対論的なしきい値を求めよ。
:::solution
$K_{\rm th}=-Q(1+m_a/m_A)=4(1+1/3)=16/3\,\mathrm{MeV}$ である。標的静止では全体の移動にも運動エネルギーが要る。
:::
::::

::::exercise{level=basic}
$E_{\rm cm}=2\,\mathrm{MeV}$、$Q=1\,\mathrm{MeV}$ の出口相対エネルギーを答えよ。
:::solution
重心系で全運動量はゼロである。入力の相対エネルギーに反応で解放する分を足し、$E'_{\rm cm}=2+1=3\,\mathrm{MeV}$ となる。
:::
::::

::::exercise{level=basic}
S 波のみで $|S_0|=1/2$ のとき、反応断面積を $\pi/k^2$ の倍数で表せ。
:::solution
吸収割合は $1-|S_0|^2=1-1/4=3/4$ であり、縮退度が 1 なので断面積は $(3/4)\pi/k^2$ である。
:::
::::

::::exercise{level=basic}
共鳴の部分幅が $2,3,5\,\mathrm{keV}$ のとき最後のチャネルへの分岐を求めよ。
:::solution
全幅は和の $10\,\mathrm{keV}$ である。部分率と幅は比例するので分岐は $5/10=1/2$ となる。
:::
::::

::::exercise{level=advanced}
$Q<0$ とする。:eqref[eq-nuc-reaction-threshold] の $K_{\rm th}$ を $Q$ で展開し、非相対論式との差を求めよ。
:::solution
$m_b+m_B=m_a+m_A-Q/c^2$ を平方して代入する。$Q$ が負なら $\max$ の第二項が選ばれる。一次の項は $-Q(m_a+m_A)/m_A$、二次の項は $Q^2/(2m_Ac^2)$ である。一次までが非相対論式なので、しきい値は相対論的な補正によって二次の項だけ高くなる。
:::
::::

::::exercise{level=advanced}
$S_\ell=\eta_\ell e^{2i\delta_\ell}$、$0\le\eta_\ell\le1$ とする。固定 $\eta_\ell$ で弾性部分の最大値と最小値を求めよ。
:::solution
$|1-S_\ell|^2=1+\eta_\ell^2-2\eta_\ell\cos(2\delta_\ell)$ である。$\cos(2\delta_\ell)=-1$ のとき最大 $(1+\eta_\ell)^2$、$+1$ のとき最小 $(1-\eta_\ell)^2$ である。いずれも反応部分 $1-\eta_\ell^2$ は共通である。
:::
::::

::::exercise{level=advanced}
共鳴近傍で幅と $k_a$ を固定できるとする。:eqref[eq-nuc-breit-wigner] の比の半値全幅を求めよ。
:::solution
ピークでは分母が $\Gamma^2/4$ である。半値のとき分母はその 2 倍なので $|E-E_R|=\Gamma/2$ となる。左右の間隔は $\Gamma$ である。
:::
::::

::::exercise{level=advanced}
完全吸収する波を $\ell=0,\ldots,L$ とする。反応断面積が $\pi(L+1)^2/k^2$ となることを示せ。
:::solution
各波の吸収割合は 1 であり、縮退度を足すと $\sum_{\ell=0}^L(2\ell+1)=L(L+1)+(L+1)=(L+1)^2$ になる。これを :eqref[eq-nuc-optical-cross-sections] に代入すればよい。
:::
::::

## 確認テスト

::::quiz{answer="6" tol="1e-6"}
$Q=-3\,\mathrm{MeV}$、標的と入射粒子が等質量のとき、非相対論のしきい値を MeV で求めよ。
:::explanation
実験室エネルギーの半分が相対運動に使える。相対運動で $3\,\mathrm{MeV}$ が必要なので実験室では $6\,\mathrm{MeV}$ である。
:::
::::

::::quiz
複素光学ポテンシャルの負の虚部が表すものを選べ。

- [x] 弾性チャネルから他のチャネルへ移る流束。
- [ ] したがって必ず増える弾性流束。
- [ ] 保存されない全電荷。

:::explanation
弾性波だけの連続の式に $-2W|\psi|^2/\hbar$ が現れる。他の出口も含む全体系では、その流束を各チャネルに配分する。
:::
::::

::::quiz{answer="9" tol="1e-6"}
$\ell=0,1,2$ が完全に吸収され、他に吸収がない場合、反応断面積は $\pi/k^2$ の何倍か。
:::explanation
各軌道の縮退度は $1,3,5$ なので和は 9 である。
:::
::::

::::quiz
$(d,p)$ ストリッピングで標的に移る粒子を選べ。

- [ ] 陽子
- [x] 中性子
- [ ] 光子

:::explanation
入射重陽子は陽子と中性子からなり、出口では陽子が観測される。核子数と電荷の差から標的が受け取った粒子は中性子である。
:::
::::

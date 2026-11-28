---
title: 輸送現象
description: 分布関数の運動から電気・熱・熱電応答を導き、散乱機構と磁場依存性を計算する。
notation:
  - symbol: 'f,\delta f'
    meaning: 電子分布と局所平衡からのずれ
  - symbol: 'q'
    meaning: キャリアの符号付き電荷
  - symbol: '\Xi(E),L_j'
    meaning: 輸送分布とそのエネルギーモーメント
  - symbol: 'S,\Pi'
    meaning: Seebeck 係数と Peltier 係数
  - symbol: 'n,p,\mu_e,\mu_h'
    meaning: 電子と正孔の密度および正の移動度
---

## この章の目標

- 半古典的な軌道と散乱から Boltzmann 輸送方程式を組み立てる。
- 同じ輸送分布から電気伝導度、熱伝導度、熱電係数を求める。
- 運動量緩和における散乱角の重みと抵抗の温度依存を説明する。
- 二種類のキャリアによる Hall 係数と磁気抵抗を求める。

エネルギーは SI で扱う。電子は $q=-e$、正孔は $q=+e$ である。磁場がない場合は等方応答を中心に扱うが、速度の積をテンソルで残せば異方的な場合にも同じ変数が使える。

分布 $f$ は一粒子状態の占有確率である。バンドとスピンを区別する添字を省いて書く場合も、電流などの状態和ではそれらを数える。古典的な軌道を使うのは波束の中心の運動であり、占有確率には Fermi–Dirac 分布と Pauli の排他原理を使う。外場は波束の大きさに比べてゆっくり変化し、異なるバンドへの遷移を無視できるとする。

## Boltzmann 輸送方程式と緩和時間近似

### 位相空間の流れ

波束の中心が式 :eqref[eq-ss-semiclassical] に従い、分布に関わる位相の干渉を平均化できるとする。散乱がなければ軌道に沿って占有 $f(\vb r,\vb k,t)$ は一定である。連鎖律で時間微分し、衝突項 $I[f]$ を加えると
$$
 \partial_t f+\vb v\cdot\grad_{\vb r}f+\dot{\vb k}\cdot\grad_{\vb k}f=I[f]
 \label{eq-ss-boltzmann}
$$
である。Berry 曲率のない一様外場では位相空間体積が保存する。具体的に $\vb v$ は $\vb k$ のみの関数であり、$\grad_{\vb k}\cdot(\vb v\times\vb B)=\vb B\cdot(\grad_{\vb k}\times\vb v)=0$ となる。ここで $\vb v$ がバンドエネルギーの勾配であることを使った。

弾性散乱率を $W_{\vb k\vb k'}$ とすると、流入と流出には $f_{\vb k'}(1-f_{\vb k})$ と $f_{\vb k}(1-f_{\vb k'})$ が現れる。対称な率 $W_{\vb k\vb k'}=W_{\vb k'\vb k}$ の場合は積が相殺し、$I[f]=\sum_{\vb k'}W_{\vb k\vb k'}(f_{\vb k'}-f_{\vb k})$ となる。これを全波数について足すと零で、粒子数は保存する。局所平衡 $f_0$ からのずれだけが速さ $1/\tau$ で減ると近似して $I=-\delta f/\tau$ とする。密度の変化も扱うときは $\sum\delta f/\tau=0$ となるよう局所平衡の化学ポテンシャルを選ぶ必要がある。

:::example[散乱なしの時間発展]
一様電場のみで $f(\vb k,0)=f_{\rm init}(\vb k)$ なら、特性曲線 $\vb k(t)=\vb k(0)+q\vb E t/\hbar$ 上で $f$ が一定なので $f(\vb k,t)=f_{\rm init}(\vb k-q\vb E t/\hbar)$ である。
:::

:::example[定常値への接近]
一様系で駆動を一定 $D(\vb k)$ とすると線形化した式は $\partial_t\delta f+D=-\delta f/\tau$ である。初期値零から $\delta f=-\tau D(1-e^{-t/\tau})$ となり、定常値は $-\tau D$ である。
:::

### 粒子数を保存する方程式

通常の位相空間測度が保存されるとき、式 :eqref[eq-ss-boltzmann] の左辺を $\partial_t f+\grad_{\vb r}\cdot(\vb v f)+\grad_{\vb k}\cdot(\dot{\vb k}f)$ と書ける。時間微分の項を全状態で足すと数密度 $n(\vb r,t)=V^{-1}\sum_{\vb k}f$ の変化率になり、第 2 項は粒子流密度 $\vb j_N=V^{-1}\sum_{\vb k}\vb v f$ の発散になる。第 3 項は積分領域の境界の流束である。Bloch バンドでは反対側のゾーン境界を同一視するので二つの項が相殺し、連続放物線分散なら無限遠での分布の減衰を仮定すれば零である。従って粒子数を失わない衝突に対して

$$
\partial_t n+\grad_{\vb r}\cdot\vb j_N=0
\label{eq-ss-boltzmann-continuity}
$$

が成り立つ。電荷 $q$ が一定のキャリアだけを扱う場合、$q$ を掛ければ電荷の連続の式となる。

緩和時間近似で参照分布 $f_0$ を外部の一様平衡に固定すると、全状態で $\sum(f-f_0)/\tau$ が零とは限らない。例えば散乱するだけの空間一様系で粒子を少し加えたとき、$\tau$ が一定なら固定 $f_0$ への緩和項は粒子の余分な部分を失わせる。そこで浴の温度を固定し、各点の $\mu_{\rm loc}$ を $\sum f_0(E;\mu_{\rm loc},T)/\tau=\sum f/\tau$ となるよう選ぶ。有限の状態和、正温度、$0<\tau<\infty$ の下では、左辺は $\mu_{\rm loc}$ に対して連続で狭義単調増加し、$0$ と $\sum1/\tau$ の間を動く。従って右辺がこの区間の内部にあれば解は一つである。外浴がなければエネルギーも保存される場合があり、そのときには温度など他のパラメータも衝突の保存則に合わせる必要がある。

:::example[均一な密度変化は減らない]
一定 $\tau$、一様な無電場系で $f=f_0(E;\mu+\delta\mu,T)$ とする。弾性散乱の各対は同じエネルギー上にあるので $f_{\vb k'}-f_{\vb k}=0$ であり、実際の衝突項は零である。局所化学ポテンシャルに $\mu+\delta\mu$ を用いる緩和時間近似でもずれは零になる。一方、電場によって生じる一次のずれは反転対称な帯では $\vb k$ に関して奇であり、$\tau(\vb k)=\tau(-\vb k)$ なら数密度も $\sum\delta f/\tau$ も変えない。
:::

### 詳細釣り合いと緩和するモード

熱浴とエネルギーを交換する過程では、遷移率自体に熱浴の占有が含まれる。状態を $a,b$、$a$ から $b$ への率を $W_{a\to b}$ とする。温度 $T$ にある浴に対して、平衡では各対の正味の流れが消えると仮定する。これを :term[詳細釣り合い]{en="detailed balance" yomi="しょうさいつりあい"} といい、

$$
P_{ab}=W_{a\to b} f_a^0(1-f_b^0)
=W_{b\to a} f_b^0(1-f_a^0)=P_{ba}
$$

である。$\delta f_a=f_a^0(1-f_a^0)\phi_a$ と置いて、$W_{b\to a}f_b(1-f_a)-W_{a\to b}f_a(1-f_b)$ を一次まで展開する。流出側の相対補正は $(1-f_a^0)\phi_a-f_b^0\phi_b$、流入側は $(1-f_b^0)\phi_b-f_a^0\phi_a$ である。差を取ると $P_{ab}(\phi_b-\phi_a)$ が残る。従って線形化した衝突項とその二次形式は

$$
\begin{aligned}
I_a^{(1)}&=-\sum_bP_{ab}(\phi_a-\phi_b),\\
-\sum_a\phi_a I_a^{(1)}
&=\frac12\sum_{a,b}P_{ab}(\phi_a-\phi_b)^2\ge0 .
\end{aligned}
\label{eq-ss-transport-collision-form}
$$

第 2 式では $a,b$ を交換した和を平均した。従って衝突は状態間の占有のずれをならす。$\phi$ が全状態で一定なら、化学ポテンシャルの微小変化に対応して零モードとなる。弾性散乱だけの場合、エネルギーごとに独立な平均のずれも零モードとなる。これが、全ての分布変化に同じ有限緩和時間を当てはめると保存則を失う理由である。

:::example[二状態の均し方]
同じエネルギーをもつ二状態が互いに同じ率 $\Gamma$ で遷移する。衝突項は $I_1=\Gamma(f_2-f_1)$、$I_2=-I_1$ である。和は一定だが、差は $\partial_t(f_1-f_2)=-2\Gamma(f_1-f_2)$ に従い、緩和時間は $1/(2\Gamma)$ である。一方から出る遷移率だけを反転した $1/\Gamma$ とは因子が異なる。
:::

### 軌道に沿う分布と交流の遅れ

式 :eqref[eq-ss-boltzmann] で衝突項を $-[f-f_0]/\tau$ とし、$\tau$ は一定とする。ある時刻 $t$ の $\vb r,\vb k$ を通る特性軌道を $X(t')=(\vb r(t'),\vb k(t'))$ と書く。軌道に沿った方程式は一階線形なので

$$
f(X(t),t)=e^{-(t-t_i)/\tau}f(X(t_i),t_i)
+\int_{t_i}^t\frac{\dd t'}{\tau}e^{-(t-t')/\tau}f_0(X(t'),t')
\label{eq-ss-transport-characteristic}
$$

となる。初期分布の記憶は指数的に減る一方、現在の分布には過去の軌道の平衡分布が加わる。例えば熱平衡の進む位置が平均自由行程にわたり違うと、この積分を現在位置の勾配だけで置き換える精度が悪くなる。境界を横切る有限試料では、特性が初期時刻より先に境界に当たる場合、入射分布を境界値として用いる。

一様電場が $e^{-i\omega t}$ に比例し、磁場を零にすると、(線形化した) 移流の駆動は $q\vb E\cdot\vb v\,\partial_Ef_0$ である。解は $\delta f=-q\tau\vb E\cdot\vb v\,\partial_E f_0/(1-i\omega\tau)$ となる。$\tau$ のほかの部分が直流と同じなら、各エネルギーの寄与には同じ $1/(1-i\omega\tau)$ が掛かる。$\tau(E)$ が異なる場合はこの因子も状態ごとに異なり、直流の一つの $\tau$ で交流スペクトルまで再現できるとは限らない。

:::example[異なる二つの緩和時間]
二つの独立な伝導成分が等しい直流伝導度 $\sigma_*$ を持ち、緩和時間が $\tau$ と $3\tau$ であるとする。$\omega\tau=1$ での実部は $\sigma_*[1/(1+1)+1/(1+9)]=0.6\sigma_*$ であり、全直流値に対する比は $0.30$ である。虚部は $\sigma_*[1/(1+1)+3/(1+9)]=0.8\sigma_*$ である。低周波の虚部の傾きに合う単一時間は $2\tau$ だが、それを使った同じ周波数での実部の比は $1/(1+4)=0.20$ となる。
:::

## 電気伝導度と熱伝導度

### 分布の線形化

$w=-\partial_E f_0\ge0$ と書く。$f_0(E,\mu,T)$ の変分は $\grad_{\vb r}f_0=w[\grad\mu+(E-\mu)\grad T/T]$ であり、電場の移流項は $q\vb v\cdot\vb E\,\partial_E f_0$ である。弱い駆動で $\delta f$ の移流を二次として落とすと
$$
 \delta f=\tau w\,\vb v\cdot
 \left[q\vb{\mathcal E}-(E-\mu)\frac{\grad T}{T}\right],
 \qquad \vb{\mathcal E}=\vb E-\grad\mu/q .
$$
熱流はエネルギー流から化学ポテンシャルに粒子流を掛けた量を引く。電流と熱流は、それぞれ $q\vb v\delta f$ と $(E-\mu)\vb v\delta f$ の波数和を単位体積で取る。等方体では
$$
 \Xi(E)=\frac1{3V}\sum_{n,\vb k,s}
 \tau_{n\vb k}v_{n\vb k}^2\delta(E-E_{n\vb k}),\quad
 L_j=\int\Xi(E)(E-\mu)^j w(E)\dd E
$$
と置くと
$$
 \vb j=q^2L_0\vb{\mathcal E}-qL_1\grad T/T,\qquad
 \vb j_Q=qL_1\vb{\mathcal E}-L_2\grad T/T .
 \label{eq-ss-transport-moments}
$$
従って $\sigma=q^2L_0$、開回路熱伝導度は $\kappa=[L_2-L_1^2/L_0]/T$ である。正の重み $\Xi w$ について Cauchy–Schwarz を使うと $L_0L_2\ge L_1^2$、従って $\kappa\ge0$ となる。

:::example[Drude 値への帰着]
放物線帯で $\tau$ 一定なら $\Xi(E_F)=g(E_F)(2E_F/m)\tau/3$ である。$g(E_F)=3n/(2E_F)$ を代入すると $n\tau/m$ となり、低温の $\sigma$ は式 :eqref[eq-ss-drude] の直流極限に一致する。
:::

:::example[一つのエネルギーが運ぶ熱]
励起が狭いエネルギー $E_1$ だけで流れる模型では $L_1=(E_1-\mu)L_0$、$L_2=(E_1-\mu)^2L_0$ なので $\kappa=0$ である。全キャリアの熱流が粒子流に比例し、電流を止めると熱流も消える。
:::

### 異方的な伝導と熱流

速度の角度平均をとらず、輸送分布 $\Xi_{ij}(E)$ の状態和で $v^2/3$ を $v_i v_j$ に置き換える。ベクトル $\vb a$ に対して $a_i\Xi_{ij}a_j$ は $\tau(\vb a\cdot\vb v)^2$ のエネルギー分解なので非負である。積分したテンソルを $\mathsf L_j$ とすると、式 :eqref[eq-ss-transport-moments] の $L_j$ を置き換えるだけで電流と熱流が得られる。$\mathsf L_0$ が正定値であると仮定すると、電荷電流をすべての方向で零とした場合に

$$
\boldsymbol\sigma=q^2\mathsf L_0,\qquad
\boldsymbol\kappa
=\frac1T(\mathsf L_2-\mathsf L_1\mathsf L_0^{-1}\mathsf L_1)
\label{eq-ss-transport-tensor}
$$

となる。行列が交換するとは限らないので積の順序を保つ。

この式の正値性は、任意の実ベクトル $\vb a,\vb b$ に対して
$$
\vb a^{\mathsf T}\mathsf L_0\vb a
+2\vb a^{\mathsf T}\mathsf L_1\vb b
+\vb b^{\mathsf T}\mathsf L_2\vb b
=\int w(E)[\vb a+(E-\mu)\vb b]^{\mathsf T}
\boldsymbol\Xi(E)[\vb a+(E-\mu)\vb b]\dd E\ge0
$$
となることから分かる。左辺で $\vb a=-\mathsf L_0^{-1}\mathsf L_1\vb b$ と選ぶと、$\boldsymbol\kappa$ の二次形式に $T$ を掛けた値が残る。ある方向に速度成分を持つ状態がなく $\mathsf L_0$ が特異なら、伝導する方向に制限して電気化学場を消去する。

:::example[有効質量の違い]
楕円放物線分散 $E=\sum_i p_i^2/(2m_i)$、一定 $\tau$ をとる。$v_i=p_i/m_i$ なので $\partial_{p_i}f_0=-w v_i$ である。$\sigma_{ij}=q^2\tau V^{-1}\sum wv_iv_j$ の連続積分で $v_i\partial_{p_j}f_0$ を部分積分すると、$\partial_{p_j}v_i=\delta_{ij}/m_i$ から $\sigma_{ij}=nq^2\tau\delta_{ij}/m_i$ を得る。この結果に低温近似は要らない。例えば $m_y=4m_x$ なら $\sigma_y=\sigma_x/4$ である。同じ電流密度を流すとき、$y$ 方向の電場は $x$ 方向の四倍になる。
:::

### 密度勾配による拡散

一定温度で、数密度を $\mu$ の関数 $n(\mu,T)$ とする。熱力学的な微分 $\chi=(\partial n/\partial\mu)_T>0$ を使えば $\grad\mu=\grad n/\chi$ である。式 :eqref[eq-ss-transport-moments] で熱勾配を零にし $q$ で割ると

$$
\vb j_N=q L_0\vb E-\frac{L_0}{\chi}\grad n,\qquad
D=\frac{L_0}{\chi},\qquad \sigma=q^2D\chi .
\label{eq-ss-generalized-einstein}
$$

従って無電場の小さい密度変動では式 :eqref[eq-ss-boltzmann-continuity] から $\partial_t n=D\nabla^2n$ となる。ただし一定係数とみなせる平衡のまわりで線形化し、電荷変動が新たに作る電場を無視するか、電荷の補償がある場合である。非縮退極限で $n\propto e^{\mu/(k_BT)}$ なら $\chi=n/(k_BT)$、移動度の大きさを $\mu_{\rm mob}=\sigma/(n|q|)$ として $D/\mu_{\rm mob}=k_BT/|q|$ である。一方、縮退三次元の放物線帯では $\chi=g(E_F)$、$L_0\simeq g(E_F)v_F^2\tau/3$ なので $D=v_F^2\tau/3$ になる。二つの極限では熱速度と Fermi 速度のどちらが拡散距離を支配するかが異なる。

:::example[密度の正弦波の減衰]
一次元の空間に $\delta n(x,0)=A\cos(kx)$ の小さな変動を与えると、$\partial_x^2\cos(kx)=-k^2\cos(kx)$ なので $\delta n(x,t)=A e^{-Dk^2t}\cos(kx)$ である。$D=10^{-4}\,\mathrm{m^2\,s^{-1}}$、波長 $1.0\,\mathrm{mm}$ なら時定数 $1/(Dk^2)$ は約 $0.25\,\mathrm{ms}$ になる。波長が平均自由行程と同程度になると、勾配展開で得たこの拡散式の精度は下がる。
:::

### Lorenz 数の低温補正

$t=k_BT$、$x=(E-\mu)/t$ と置く。:term[Lorenz 数]{en="Lorenz number" yomi="ろーれんつすう"} を電気伝導度と電子熱伝導度から $\mathcal L=\kappa/(\sigma T)$ と定める。化学ポテンシャルが帯の端から $t$ より十分遠く、$\Xi$ が滑らかで $\Xi(\mu)>0$ とする。$w(E)\dd E=e^x/(1+e^x)^2\dd x$ は、積分範囲を実軸全体に延長したとき偶関数である。その二次・四次モーメントを $a_2,a_4$ とすると $a_2=\pi^2/3$、$a_4=7\pi^4/15$ である。四次の値は式 :eqref[eq-stm-fermi-gas] の温度四次係数に $4!$ を掛けても得られる。

$\Xi$ を $\mu$ の周りで展開し、奇数次の積分が消えることを用いると

$$
\begin{aligned}
L_0&=\Xi+\tfrac12 a_2 t^2\Xi''+\order{t^4},\\
L_1&=a_2t^2\Xi'+\order{t^4},\\
L_2&=a_2t^2\Xi+\tfrac12 a_4t^4\Xi''+\order{t^6}
\end{aligned}
$$

となる。エネルギー微分は $\mu$ で評価する。この三式を $\mathcal L=(L_2/L_0-L_1^2/L_0^2)/(q^2 T^2)$ に代入して分母の逆数を展開すれば

$$
\mathcal L=\frac{k_B^2}{q^2}\left[
a_2+t^2\left(\frac{a_4-a_2^2}{2}\frac{\Xi''}{\Xi}
-a_2^2\left(\frac{\Xi'}{\Xi}\right)^2\right)
+\order{t^4}\right]
\label{eq-ss-transport-lorenz-correction}
$$

となる。先頭項は式 :eqref[eq-ss-wiedemann-franz] である。第 2 項の最初には輸送分布の曲率が、第 2 項の後者には開回路条件が効く。この補正は電子が共通の弾性的な衝突演算子で緩和する領域の結果である。

:::example[電子と正孔が運ぶ熱]
二キャリアが個別には $\vb j_i=\sigma_i(\vb{\mathcal E}-S_i\grad T)$、$\vb j_{Q,i}=S_iT\vb j_i-\kappa_i\grad T$ に従うとする。電子と正孔の生成・再結合が局所平衡にあり、それぞれの電気化学ポテンシャルを電荷で割った電位が共通である場合を考える。この共通の駆動場と一様な温度勾配のもとで全電流を零とすると、$\vb{\mathcal E}=(\sigma_1 S_1+\sigma_2 S_2)\grad T/(\sigma_1+\sigma_2)$ となる。熱流へ代入すると
$$
\kappa_{\rm tot}=\kappa_1+\kappa_2+
 T\frac{\sigma_1\sigma_2}{\sigma_1+\sigma_2}(S_1-S_2)^2
$$
である。補正項は二つの電流が互いに打ち消し合いながら熱を運ぶ寄与である。$\sigma_1=\sigma_2=10^4\,\mathrm{S\,m^{-1}}$、$S_1=-S_2=200\,\mathrm{\mu V\,K^{-1}}$、$T=300\,\mathrm K$ なら追加分は $0.24\,\mathrm{W\,m^{-1}K^{-1}}$ となる。
:::

### 電気と熱の駆動を同時にかける

温度と化学ポテンシャルがゆっくり変わるとき、局所のエントロピー生成率は $\dot s_{\rm prod}=\vb j\cdot\vb{\mathcal E}/T-\vb j_Q\cdot\grad T/T^2$ である。熱流が運ぶエントロピー流を $\vb j_Q/T$ とし、定常のエネルギー収支 $\grad\cdot\vb j_Q=\vb j\cdot\vb{\mathcal E}$ を使って発散をとれば、この二つの項が得られる。ただし化学反応など別のエントロピー源は省いている。

式 :eqref[eq-ss-transport-moments] から等方体では

$$
T\dot s_{\rm prod}
=\int \Xi(E) w(E)
\left|q\vb{\mathcal E}-(E-\mu)\grad T/T\right|^2\dd E
\ge0
$$

である。電気と熱の応答の交差項は、温度勾配や電場の向きによって正にも負にもなれるが、全生成率は非負である。開回路では電気化学場を変化させて電流を零にするので、同じ温度勾配での二次形式を最小にする。平方完成すると残る係数は $\kappa/T$ になる。従って熱伝導の補正とエントロピーの正値性は同じ積分に由来する。

## 熱電効果（Seebeck 効果・Peltier 効果）

### 二つの境界条件

ここでは開回路の電気化学場で $S=\mathcal E_x/\partial_x T$ と定義する。式 :eqref[eq-ss-transport-moments] に $\vb j=0$ を入れると $S=L_1/(qTL_0)$ である。温度勾配をかけずに電流を流したときの熱流の比が $\Pi=j_Q/j=L_1/(qL_0)$ である。したがって $\Pi=ST$ となる。この等式は両係数が同じ交差モーメントに由来することから出る。異種材料の界面で電流が連続し $\Pi$ が変われば、熱流の差を界面での発熱・吸熱が補う。

低温で $\Xi(E)$ が滑らかなら $\Xi(\mu)+\Xi'(\mu)(E-\mu)$ と展開する。$w$ は $\mu$ 周りで偶関数なので $L_0\simeq\Xi(\mu)$、$L_1\simeq(\pi^2/3)(k_BT)^2\Xi'(\mu)$ である。二次モーメントには式 :eqref[eq-ss-sommerfeld-expansion] の導出を用いた。従って
$$
 S=\frac{\pi^2 k_B^2T}{3q}\left.\dv{\log\Xi}{E}\right|_{E=\mu}
 \label{eq-ss-mott}
$$
が低温の先頭項である。これは散乱のエネルギー依存も含むので符号は電荷だけで一意に決まらない。

:::example[等方放物線帯]
三次元で $\tau\propto E^r$ なら $g\propto E^{1/2}$、$v^2\propto E$ より $\Xi\propto E^{r+3/2}$ である。$r=0$ の電子なら $S=-(\pi^2 k_B^2T)/(2eE_F)$ である。$E_F=5\,\mathrm{eV}$、$T=300\,\mathrm K$ では約 $-2.2\,\mathrm{\mu V\,K^{-1}}$ となる。
:::

:::example[界面の熱]
$T=300\,\mathrm K$、$\Delta S=100\,\mathrm{\mu V\,K^{-1}}$ であれば $\Delta\Pi=0.030\,\mathrm V$ である。$I=0.10\,\mathrm A$ が界面を通ると熱流の差の大きさは $0.0030\,\mathrm W$ となる。向きを逆にするとこの寄与の符号も変わる。
:::

### 測定される起電力と界面の向き

温度差で生じる電圧を測る :term[Seebeck 効果]{en="Seebeck effect" yomi="ぜーべっくこうか"} では、電気化学電位を電荷で割った電位 $\varphi_{\rm ec}=\varphi+\mu/q$ を考える。静電場は $-\grad\varphi$ なので $\vb{\mathcal E}=-\grad\varphi_{\rm ec}$ である。一本道の試料が開回路なら、冷端から熱端までの電気化学電圧は

$$
\varphi_{\rm ec}(T_h)-\varphi_{\rm ec}(T_c)
=-\int_{T_c}^{T_h}S(T)\dd T
$$

となる。計器につなぐリードにも温度勾配があれば、リードの積分が加わる。異種線 A と B で同じ $T_h$ の接点を作り、冷端の温度を共に $T_c$ とすると、A の冷端電位から B の冷端電位を引いた電圧は $\int_{T_c}^{T_h}(S_A-S_B)\dd T$ になる。接点の電気化学電位を共通の基準として、各線を熱端から冷端へ積分した結果である。

:term[Peltier 効果]{en="Peltier effect" yomi="ぺるちぇこうか"} では界面を横切る電流 $I$ と熱流の向きをそろえて符号を取る。A から B に流れる場合、輸送される熱の出ていく分から入る分を引くと $I(\Pi_B-\Pi_A)$ である。接点の温度を保つには、この量を接点へ外部から流し込む必要がある。したがって正なら接点から熱が奪われる向き、負なら接点に熱を与える向きである。Joule 熱は $I^2 R$ と電流反転に対して偶であり、界面の Peltier 熱は $I$ に一次なので奇である。

:::example[熱電対の電圧]
$S_A-S_B=40\,\mathrm{\mu V\,K^{-1}}$ が一定で、$T_h-T_c=50\,\mathrm K$ のとき冷端電圧は $2.0\,\mathrm{mV}$ である。$S_A-S_B$ が温度区間内で $30$ から $50\,\mathrm{\mu V\,K^{-1}}$ まで直線的に増える場合も平均が $40\,\mathrm{\mu V\,K^{-1}}$ なので同じ値になる。最高温度での係数に温度差を掛けると、この場合は過大評価になる。
:::

### Thomson 熱とエネルギー収支

空間方向を一次元に限り、断面を一定とする。定常で電流密度 $j$ は一定である。電気化学場と熱流を $j$ で表すと

$$
\mathcal E_x=j/\sigma+S\partial_xT,\qquad
j_Q=\Pi j-\kappa\partial_xT .
$$

試料の側面から熱が逃げず、電子と格子が同じ局所温度を持つとする。格子が熱を運ぶ場合、$j_Q$ と $\kappa$ にはその寄与も含める。定常で電場がする仕事率は $jE_x$ なので、電子と格子を合わせたエネルギー流は $\partial_xj_U=jE_x$ を満たす。熱流からは化学仕事を引くため、$j_Q=j_U-\mu j/q$ によって $\partial_x j_Q=j\mathcal E_x$ となる。上の式を代入し、均質材で $\Pi=ST$ が温度だけで変わるとすれば

$$
-\partial_x(\kappa\partial_xT)=j^2/\sigma
-T\frac{\dd S}{\dd T}j\,\partial_xT
\label{eq-ss-transport-thomson}
$$

となる。第 2 項が :term[Thomson 効果]{en="Thomson effect" yomi="とむそんこうか"} である。局所の発熱項としての符号は式の通りであり、発熱から吸熱へ切り替わる方向は $j$、$\partial_x T$、$\dd S/\dd T$ の符号で決まる。$S$ が一定のときも接合で $\Pi$ が飛べば Peltier 熱が生じる。

:::example[Thomson 項の桁]
$S(T)=aT$、$a=-0.10\,\mathrm{\mu V\,K^{-2}}$、$T=300\,\mathrm K$ とする。$j=10^4\,\mathrm{A\,m^{-2}}$、$\partial_xT=10^3\,\mathrm{K\,m^{-1}}$ なら式 :eqref[eq-ss-transport-thomson] の Thomson 発熱は $300\,\mathrm{W\,m^{-3}}$ である。例えば $\sigma=10^5\,\mathrm{S\,m^{-1}}$ なら Joule 発熱は $10^3\,\mathrm{W\,m^{-3}}$ である。電流を反転すれば前者だけが $-300\,\mathrm{W\,m^{-3}}$ へ変わる。
:::

### 熱電変換と負荷

長さ $\ell$、断面積 $A$ の均質熱電材料を、両端が $T_h,T_c$ の熱浴につながるように置く。相手側の導線に熱電係数や熱伝導がない理想化を行い、材料の $S,\sigma,\kappa_{\rm tot}$ は一定とする。ここで $\kappa_{\rm tot}$ は電子と格子の熱伝導を合わせた値である。起電力の大きさは $|S|\Delta T$、内部抵抗 $R=\ell/(\sigma A)$、熱コンダクタンス $K=\kappa_{\rm tot}A/\ell$ である。外部抵抗を $R_L$ とすれば $|I|=|S|\Delta T/(R+R_L)$、負荷の受け取る電力は

$$
P_L=I^2R_L=\frac{(S\Delta T)^2R_L}{(R+R_L)^2}
$$

である。$R_L$ で微分すると分子に $R-R_L$ が残るので、最大電力は $R_L=R$ のとき $(S\Delta T)^2/(4R)$ である。一方、熱を漏らす経路があるので、電力最大と熱効率最大は一般に一致しない。

この模型の高温端からの熱流入は、電流の向きを $S I>0$ に取ると $Q_h=S T_h I+K\Delta T-I^2R/2$ である。第 1 項は高温端の Peltier 熱である。第 2 項と第 3 項は、$\kappa_{\rm tot}T''=-j^2/\sigma$ に両端温度を与えて高温端の温度勾配から得られる。例えば座標 $x=0$ を高温端、$x=\ell$ を低温端に置くと、$T'(0)=-\Delta T/\ell+j^2\ell/(2\sigma\kappa_{\rm tot})$ なので $-\kappa_{\rm tot}AT'(0)=K\Delta T-I^2R/2$ となる。従って熱伝導の漏れに対する電気変換の大きさは $Z=S^2/(RK)=S^2\sigma/\kappa_{\rm tot}$ を代表温度で無次元化した $ZT$ で比較される。$Z$ の単位は $\mathrm{K^{-1}}$ である。

:::example[熱電素子の出力と熱漏れ]
$S=200\,\mathrm{\mu V\,K^{-1}}$、$\Delta T=50\,\mathrm K$、$R=0.010\,\mathrm{\Omega}$ とする。負荷を $R$ に合わせると起電力は $0.010\,\mathrm V$、電流は $0.50\,\mathrm A$、出力は $2.5\,\mathrm{mW}$ である。$K=0.010\,\mathrm{W\,K^{-1}}$、$T_h=350\,\mathrm K$ では $Q_h=0.035+0.50-0.00125=0.53375\,\mathrm W$ なので効率は約 $0.47\%$ である。接触抵抗などを省いた値であり、熱漏れ $K\Delta T$ は電力式に現れなくても効率を左右する。
:::

### 最大電力と最大効率

同じ理想模型で $\Delta T=T_h-T_c>0$、$\overline T=(T_h+T_c)/2$、負荷比 $u=R_L/R$ とする。$S\ne0$ の場合、$I=S\Delta T/[R(1+u)]$ を $P_L$ と $Q_h$ に代入し、分母分子に $(1+u)^2/(K\Delta T)$ を掛けると効率は

$$
\eta(u)=\frac{Z\Delta T\,u}
{u^2+(2+ZT_h)u+1+Z\overline T}
$$

となる。定義域 $u\ge0$ で分母は正であり、微分の符号は $1+Z\overline T-u^2$ で決まる。従って最大効率は $u_*=\sqrt{1+Z\overline T}$ のときに得られる。この関係を分母に代入すると $\eta(u_*)=Z\Delta T/(2u_*+2+ZT_h)$ である。$Z=(u_*^2-1)/\overline T$ を入れ、分子と分母を $u_*+1$ で割れば

$$
\eta_{\max}=\frac{T_h-T_c}{T_h}
 \frac{\sqrt{1+Z\overline T}-1}{\sqrt{1+Z\overline T}+T_c/T_h}
\label{eq-ss-transport-efficiency}
$$

を得る。有限の $Z\overline T$ では後ろの因子が $1$ 未満である。$Z\overline T$ が大きくなる極限で、式は Carnot 効率に近づく。ただしこの比較では材料定数の温度依存、接触抵抗、系外への熱漏れを省いた。

:::example[同じ材料で負荷を変える]
$T_h=600\,\mathrm K$、$T_c=300\,\mathrm K$、$Z\overline T=1$ なら、最大効率の負荷は内部抵抗の $\sqrt2$ 倍である。式 :eqref[eq-ss-transport-efficiency] では $\eta_{\max}=(1/2)(\sqrt2-1)/(\sqrt2+1/2)\simeq0.108$ になる。電力最大の $u=1$ を使うと $\eta=2/19\simeq0.105$ となる。最大効率の点では電流が小さくなるため Joule 熱が抑えられる一方、負荷に渡る電力も最大値から下がる。
:::

## 不純物散乱とフォノン散乱

### 散乱角の重み

局在不純物が希薄で配置が独立なら、同じ不純物内での散乱振幅の二乗が濃度に比例して足される。弱い散乱に対する Fermi の黄金律では $W_{\vb k\vb k'}=(2\pi/\hbar)n_i|U_{\vb k'-\vb k}|^2\delta(E_{\vb k'}-E_{\vb k})/V$ となる。ここで $U_{\vb Q}=\int U(\vb r)e^{-i\vb Q\cdot\vb r}\dd^3r$ は単一不純物の Fourier 変換、$n_i$ は密度である。有限体積で平面波の行列要素が $U_{\vb Q}/V$ となり、不純物数が $n_iV$ であることを用いた。散乱が等方帯で回転対称なら $\delta f\propto\vb E\cdot\vb k$ に対する衝突項の平均は $\vb k'$ の横成分を消す。残る縦成分が $\cos\theta$ 倍であるから
$$
 \frac1{\tau_{\rm tr}}=\sum_{\vb k'}W_{\vb k\vb k'}(1-\cos\theta)
$$
を得る。寿命の逆数 $\sum W$ は小さい角度の散乱も一回として数えるため、運動量の緩和率とは異なる。

### 音響フォノンによる散乱

縦波の体積ひずみに電子エネルギーが比例するとし、変形ポテンシャル定数を $D$ とする。変位の規格化が $\sqrt{\hbar/(2\rho V\omega_q)}$、発散が $iq$ を掛けるので $|M_q|^2=D^2\hbar q^2/(2\rho V\omega_q)$ である。吸収は $n_q$、放出は $n_q+1$ で重み付けされ、エネルギー保存のデルタ関数はそれぞれ $E_{\vb k'}-E_{\vb k}=\pm\hbar\omega_q$ を要求する。縮退気体では最終状態の Pauli 禁制も重要である。高温で $n_q\simeq k_BT/(\hbar\omega_q)$、$\omega_q=v_s q$ とすれば $|M_q|^2 n_q$ は $q$ によらず $T$ に比例する。

:::example[前方散乱]
$\theta=0.10$ ラジアンで $1-\cos\theta\simeq0.0050$ であり、後方散乱では重みが $2$ となる。同じ発生率でも伝導への寄与には約四百倍の差がある。
:::

:::example[弾性散乱の残留抵抗]
弱い点状不純物で $U_{\vb Q}$ を一定とすれば、球面で $\cos\theta$ の平均は零である。従って寿命と運動量緩和時間が一致し、固定キャリア密度で残留抵抗率は $n_i$ に比例する。
:::

### 不純物配置を平均する

不純物を $\vb R_1,\dots,\vb R_{N_i}$ に置く。規格化した平面波の間の行列要素は $(U_{\vb Q}/V)\sum_a e^{-i\vb Q\cdot\vb R_a}$ である。その二乗を配置平均したとき、$a=b$ の項は $N_i$ 個あり、$a\ne b$ は配置の独立性により平均が $\left|\langle e^{-i\vb Q\cdot\vb R_a}\rangle\right|^2$ に分かれる。胞の中で位置が一様で、$\vb Q$ が非零の周期境界波数なら後者は零である。この平均化には不純物が集合していないことが効く。$\vb Q=0$ の一様成分はバンド全体のエネルギーをずらすが、その過程では速度が変わらない。

放物線分散でスピンを保存する場合、最終状態の和を $V(2\pi)^{-3}\int k'^2\dd k'\dd\Omega$ に換える。$\dd E_{\vb k'}/\dd k'=\hbar^2 k'/m$ なので、デルタ関数を使った動径積分は $mk/\hbar^2$ である。角度分解した緩和率は従って

$$
\frac1{\tau_{\rm tr}(k)}
=\frac{n_i mk}{4\pi^2\hbar^3}
\int_{\text{球面}}|U_{\vb Q}|^2(1-\cos\theta)\dd\Omega
\label{eq-ss-transport-impurity-rate}
$$

となる。ここで $\dd\Omega$ は立体角であり、遷移は同じスピンに限っているためスピンの縮退数を前に掛けない。点状不純物で $U_{\vb Q}=U_0$ と取れば角度積分は $4\pi|U_0|^2$ となり、$1/\tau_{\rm tr}\propto k$ である。従って $\tau$ を Fermi 面では一定と近似できても、バンド底から上まで一定とは限らない。

:::example[点状不純物と電気伝導の密度依存]
散乱不純物密度と $U_0$ を固定する。三次元で低温なら $\tau_{\rm tr}(k_F)\propto1/k_F$、$k_F\propto n^{1/3}$ なので $\sigma\propto n\tau_{\rm tr}\propto n^{2/3}$ である。密度を八倍にすると電気伝導度は四倍になる。キャリア密度とともに不純物密度も増える状況では、この指数には $n_i(n)$ の変化がさらに入る。
:::

### 角度分布ごとの緩和

等エネルギー面の上で、分布の偏りが $\Phi_2(\hat{\vb k})=(\hat{\vb a}\cdot\hat{\vb k})^2-1/3$ に比例するとする。ここで $\hat{\vb a}$ は固定の単位ベクトルである。分布のこの成分は方向を反転しても同じであり、例えば四極子的な偏りを表す。散乱角を $\theta$ と固定したとき、出射方向 $\hat{\vb k}'$ は入射方向を軸とする円を回る。$\hat{\vb a}\cdot\hat{\vb k}=u$ とすると、出射との内積は $u\cos\theta+\sqrt{1-u^2}\sin\theta\cos\varphi$ と書ける。方位角 $\varphi$ に関して平均すると

$$
\left\langle(\hat{\vb a}\cdot\hat{\vb k}')^2-\tfrac13\right\rangle_\varphi
=\frac{3\cos^2\theta-1}{2}\bigl(u^2-\tfrac13\bigr)
$$

となる。交差項は消え、$\cos^2\varphi$ の平均は $1/2$ である。従ってこの偏りの緩和率では、電流に対する $1-\cos\theta$ の代わりに $1-(3\cos^2\theta-1)/2$ が散乱角の重みになる。

:::example[後方散乱を受ける二つの偏り]
散乱が正確な後方 $\theta=\pi$ に起こる模型で、その率を $W$ とする。電流に関わる奇の偏りは符号を変えるので、緩和率は $2W$ である。一方、$\Phi_2$ は $\hat{\vb k}\mapsto-\hat{\vb k}$ に対して偶なので、上の緩和率は零である。同じ衝突でも、どの分布の偏りを測るかによって緩和時間は変わる。
:::

### 遮蔽された長距離散乱

単一不純物のポテンシャルを $U(r)=C e^{-q_s r}/(4\pi r)$、$q_s>0$ とする。遮蔽の起源は [バンド計算の方法と電子間相互作用](/physics/solid-state/band-calculation/) で見る。この模型では半径方向の Fourier 積分は $U_Q=(C/Q)\int_0^\infty e^{-q_s r}\sin(Qr)\dd r=C/(Q^2+q_s^2)$ である。大きい波数移行ほど振幅が小さくなるため、電流を緩和しにくい小角度の散乱が多い。

$Q=2k\sin(\theta/2)$ を使うと $\dd\Omega=2\pi Q\dd Q/k^2$、$1-\cos\theta=Q^2/(2k^2)$ である。式 :eqref[eq-ss-transport-impurity-rate] に代入すると、$1/\tau_{\rm tr}$ に必要な角度積分は $\pi C^2/k^4$ に

$$
\int_0^{2k}\frac{Q^3\dd Q}{(Q^2+q_s^2)^2}
=\frac12\left[\log(1+b)-\frac{b}{1+b}\right],\qquad b=\frac{4k^2}{q_s^2}
$$

を掛けた値になる。$u=Q^2+q_s^2$ と置けば、被積分関数は $(1/u-q_s^2/u^2)/2$ である。一方、全遷移率 $1/\tau_{\rm q}$ では角度重みを $1$ に置き換えるため、積分は $\int_0^{2k}Q\dd Q/(Q^2+q_s^2)^2=b/[2q_s^2(1+b)]$ である。よって

$$
\frac{\tau_{\rm q}}{\tau_{\rm tr}}
=\frac{2}{b}\left[\frac{1+b}{b}\log(1+b)-1\right]
\label{eq-ss-transport-forward-ratio}
$$

となる。$b\to0$ の強い遮蔽では右辺は $1$ に近づき、$b\to\infty$ の弱い遮蔽では $\log b/b$ 程度に小さくなる。後者では全遷移は多くても方向を反転する機会が少ない。

:::example[二つの寿命を比較する]
$2k/q_s=10$ なら $b=100$ であり、$\tau_{\rm q}/\tau_{\rm tr}\simeq0.0732$ となる。例えば全遷移からの寿命が $10\,\mathrm{fs}$ なら運動量緩和時間は約 $137\,\mathrm{fs}$ である。散乱時間を比較する実験では、線幅が主に全遷移を数える状況と、電流の減衰を測る状況を区別する。
:::

### フォノンによる吸収と放出

平衡フォノンの分布 $n_Q$ について、$n_Q+1=e^{\hbar\omega_Q/(k_BT)}n_Q$ である。この関係と
$$
\frac{f_0(E)[1-f_0(E+\epsilon)]}
 {f_0(E+\epsilon)[1-f_0(E)]}
=e^{\epsilon/(k_BT)}
$$
を合わせると、$\epsilon=\hbar\omega_Q$ の吸収と逆向き放出で $P_{ab}$ が一致する。従って場がない状態では散乱だけで平衡占有が変わらない。温度零で $n_Q=0$ なら吸収はなく、放出可能なエネルギー差があっても最終状態が埋まっていれば放出できない。

縦波に対する変形ポテンシャルには体積ひずみ $\grad\cdot\vb u$ が効く。等方体の横波は $\vb Q\cdot\vb e_Q=0$ のためこの項による一次の結合を持たない。より一般にはテンソル変形ポテンシャルと $\partial_i u_j$ が結合するため、せん断がバンド端や谷のエネルギーを変える材料では横波も関係する。

全ての関係する音響波数で $k_BT\gg\hbar\omega_Q$ とし、さらに $v_s\ll v_F$ なのでエネルギー移行を弾性近似にする。一フォノンの吸収と放出を合わせると、行列要素の重みは $|M_Q|^2(2n_Q+1)\simeq D^2 k_BT/(\rho V v_s^2)$ となる。等方帯では角度の $\cos\theta$ の平均が零となるので、同じスピンの単位体積状態密度を $g_s(E)$ とすると

$$
\tau_{\rm tr}^{-1}(E)
\simeq\frac{2\pi}{\hbar}\frac{D^2 k_BT}{\rho v_s^2}g_s(E).
$$

この式では $Q\le2k_F$ の全散乱角を許すとした。Debye の打切りがそれより小さい場合は許された角度へ積分を制限する。同じ形の $T$ 比例でも係数は変わる。電子の縮退を保って Fermi 面で評価するには $k_BT\ll E_F$ も必要なので、ここでいう高温はフォノンエネルギーに対する高温である。

:::example[光学フォノンのしきい値]
分散のない光学分枝のエネルギーを $\hbar\omega_0=40\,\mathrm{meV}$ とする。ある電子のエネルギーが帯の底から $25\,\mathrm{meV}$ なら、その帯の中への一光学フォノン放出後のエネルギーは $-15\,\mathrm{meV}$ となり、終状態がない。吸収なら最終エネルギーは $65\,\mathrm{meV}$ であるが、ボース重みが掛かる。例えば $k_BT=10\,\mathrm{meV}$ なら $n_0=1/(e^4-1)\simeq0.0187$ である。縮退時には、この判定にさらに終状態の空席も必要である。
:::

## Matthiessen 則と Bloch–Grüneisen 則

### 足せる場合と低温の位相空間

独立な衝突演算子は足し合わされる。各演算子が同じ駆動モードに対して $-\delta f/\tau_i$ となるなら、$1/\tau_{\rm tr}=\sum_i1/\tau_i$ である。同じキャリア密度と質量の Drude 式を使う場合に抵抗率も足せる。これが :term[Matthiessen 則]{en="Matthiessen's rule" yomi="まてぃーせんそく"} である。エネルギー依存が異なる散乱で $\sigma$ が $\tau(E)$ の平均になる場合、率を足すことと最終抵抗率を足すことは同じ操作にならない。

三次元の球状 Fermi 面で弾性に近い音響散乱を考える。波数移行 $Q=2k_F\sin(\theta/2)$ に対し角度積分は $Q\dd Q$、角度の重みは $Q^2/(2k_F^2)$、$|M_Q|^2$ は $Q$ に比例する。さらに線形応答の衝突積分では吸収の粒子因子をエネルギーで積分し、
$$
 \int_{-\infty}^{\infty} f(E)[1-f(E+\epsilon)]\dd E=\frac{\epsilon}{1-e^{-\epsilon/(k_BT)}}
$$
を使う。実際の帯端が熱励起の領域から十分離れているとして、積分を実軸全体に延ばした。$f(E)-f(E+\epsilon)$ は被積分関数に分母 $1-e^{-\epsilon/(k_BT)}$ を掛けた式である。差を $[-A,A]$ で積分すると、変数移動により低エネルギー側の長さ $\epsilon$ の区間の占有から高エネルギー側の同じ長さの区間を引いた値となる。$A\to\infty$ で占有はそれぞれ $1,0$ へ行くため、差の積分は $\epsilon$ となる。吸収のフォノン占有 $n_Q$ と、電流の線形化の $1/(k_BT)$ を掛ければ重みは $x e^x/(e^x-1)^2$、$x=\hbar v_sQ/(k_BT)$ となる。放出も詳細釣り合いで同じ温度依存を与える。従って総合すると
$$
 \rho_{\rm ph}\propto (T/\Theta_{\rm BG})^5
 \int_0^{\Theta_{\rm BG}/T}\frac{x^5 e^x}{(e^x-1)^2}\dd x,\qquad
 \Theta_{\rm BG}=\hbar v_s\min(q_D,2k_F)/k_B .
 \label{eq-ss-bloch-gruneisen}
$$
低温では積分は有限な定数、高温では被積分関数が $x^3$ に近づき上限の四乗になるので、それぞれ $T^5$ と $T$ である。これは三次元、音響変形ポテンシャル、滑らかな Fermi 面、Fermi 速度より十分小さい音速を前提にした :term[Bloch–Grüneisen 則]{en="Bloch–Grüneisen law" yomi="ぶろっほぐりゅーないぜんそく"} である。

:::example[残留抵抗と低温外挿]
$\rho=\rho_0+AT^5$ なら $T^5$ を横軸にした切片が $\rho_0$ である。$T$ を二倍にして増分 $\rho-\rho_0$ を比べれば $32$ 倍になるが、全抵抗率の比は $32$ より小さい。
:::

:::example[小さな Fermi 面]
$q_D>2k_F$ の場合、$k_F$ を半分にすると $\Theta_{\rm BG}$ も半分になる。バンドの密度が変わると比例定数も変わり得るので、同じ温度の抵抗率が一定倍率になるとはこの式だけでは言えない。
:::

### 衝突演算子の加法と測定抵抗

式 :eqref[eq-ss-transport-collision-form] の緩和を $I^{(1)}=-\mathsf A\phi$ と書く。$\mathsf A$ は対称で、二次形式が非負である。一定温度と電場について、移流の駆動の負を $\vb b$ とすると定常解は $\mathsf A\phi=\vb b$ を満たす。有限状態で考え、$\vb b$ は数保存などの零モードに直交するとする。任意の試行分布 $\psi$ について

$$
2\langle\vb b,\psi\rangle-\langle\psi,\mathsf A\psi\rangle
=\langle\phi,\mathsf A\phi\rangle
-\langle\psi-\phi,\mathsf A(\psi-\phi)\rangle
$$

である。両辺を展開し対称性を使えば証明できる。従って左辺は $\psi=\phi$ で最大になる。電場を一方向にかけた場合、$\vb b$ は $qEv$ と平衡の重みの積であり、$\langle\vb b,\phi\rangle$ は一定の正係数を除いて電流の仕事率 $jE$ である。この変分式に特定の速度方向の分布だけを入れると、完全な分布を求める代わりにその輸送率を見積もれる。

独立な衝突機構について $\mathsf A=\mathsf A_1+\mathsf A_2$ である。しかし観測伝導度は定常解を決める $\mathsf A^{-1}$ を含む。二つの機構で駆動が共通の固有分布なら、逆行列を一つの固有値の逆数に置き換えられる。駆動が複数の固有分布を含む場合、各機構が異なる分布を選び、最後の抵抗率は加法的にならない場合がある。

:::example[率を足した結果と抵抗率を足した結果]
二つのエネルギー領域が等しい輸送重みを持ち、同じ単位で緩和時間を $\tau_1,\tau_2$、規格化した伝導度を $(\tau_1+\tau_2)/2$ とする。散乱 A の率が $(1,4)$、散乱 B の率が $(4,1)$ なら、各機構だけの伝導度はともに $(1+1/4)/2=5/8$、抵抗率はともに $8/5$ である。両方を入れると領域ごとの率は $(5,5)$ なので、伝導度は $1/5$ で抵抗率は $5$ である。各機構で得た抵抗率を加えた $16/5$ とは異なる。各機構が弱く散乱させる領域が互いに違うことが、この差を生む。
:::

### Bloch–Grüneisen 積分の両極限

式 :eqref[eq-ss-bloch-gruneisen] は衝突積分を線形化し、電子の試行分布を電場方向の速度に比例させた近似である。式 :eqref[eq-ss-transport-collision-form] の差の二乗は、等しい大きさの二つの速度について方向平均すると $|\vb v-\vb v'|^2/3=2v_F^2(1-\cos\theta)/3$ となる。ここから運動量緩和の角度重みが出る。電子のエネルギーが $\mu$ から数 $k_BT$ 以内にあり、音速が $v_F$ より十分小さいとすれば、散乱角の積分では $k'\simeq k_F$ とできる。一方、占有の積の中には $\hbar\omega_Q$ を残す。占有はエネルギー $k_BT$ の幅で変わるので、$v_s/v_F$ が小さいことからこの項も零とすることはできない。

低温側の定数は級数で求められる。$x>0$ で $e^x/(e^x-1)^2=\sum_{j\ge1}j e^{-jx}$ であり、全項が非負なので積分と和を交換できる。$x^5$ を掛けた各項の積分は $j\int_0^\infty x^5 e^{-jx}\dd x=5!/j^5$ なので

$$
\int_0^\infty\frac{x^5 e^x}{(e^x-1)^2}\dd x=120\zeta(5),
\qquad \zeta(5)=\sum_{j=1}^\infty\frac1{j^5}.
$$

上限を無限大に延長した際の誤差は、打切りを $y=\Theta_{\rm BG}/T$ として $y$ が大きいとき $y^5 e^{-y}$ 程度である。実際、$x\ge y$ では $e^{-x}\le e^x/(e^x-1)^2\le e^{-x}/(1-e^{-y})^2$ であり、$\int_y^\infty x^5 e^{-x}\dd x$ は部分積分で $e^{-y}$ と $y$ の五次多項式の積になる。高温側は $e^x/(e^x-1)^2=1/[4\sinh^2(x/2)]=x^{-2}-1/12+\order{x^2}$ であり、積分上限が $y\ll1$ なら

$$
\int_0^y\frac{x^5 e^x}{(e^x-1)^2}\dd x
=\frac{y^4}{4}-\frac{y^6}{72}+\order{y^8}
$$

となる。先頭の高温抵抗を $\rho_{\rm ph}\simeq AT$ と書くなら、相対補正は $-y^2/18$ である。

:::example[低温域と高温域を分ける]
$v_s=4.0\times10^3\,\mathrm{m\,s^{-1}}$、$\min(q_D,2k_F)=1.0\times10^{10}\,\mathrm{m^{-1}}$ なら $\Theta_{\rm BG}\simeq306\,\mathrm K$ である。$T=10\,\mathrm K$ では上限が約 $31$ なので、積分の打切りによる温度依存は小さい。$T=600\,\mathrm K$ では $y\simeq0.51$、高温近似の相対補正の先頭は約 $-0.014$ である。この見積もりは、温度を上げても音速やキャリア密度を一定とした場合に限る。
:::

### 格子が運動量を受け取る時間

音響フォノンを平衡浴とした式では、電子が渡した運動量は浴側で十分速く失われる。電子とフォノンの運動量の交換が強くても、フォノン側での損失が遅いと、フォノンが電子を引く寄与が生じる。効果を一方向の速度に縮約した模型として、電子とフォノンのドリフトを $u_e,u_p$、不純物や境界への運動量損失係数を $a,b>0$、相互の交換係数を $\gamma>0$ とする。速度は二つの運動量モードを表す変数であり、音速そのものとは異なる。定常で

$$
0=nqE-a u_e-\gamma(u_e-u_p),\qquad
0=-b u_p+\gamma(u_e-u_p)
$$

と書けば、第 2 式から $u_p=\gamma u_e/(b+\gamma)$ である。従って電子に対する有効損失係数は $a+\gamma b/(b+\gamma)$ となる。$b\gg\gamma$ ではフォノン浴に固定する近似の $a+\gamma$ に近づき、$b\to0$ では追加損失が零になる。後者ではフォノンが受け取った運動量を浴へ渡せず、定常状態で電子とともに動く。電子が直接失う運動量に対応する係数 $a$ は残る。

:::example[フォノンが取り去る運動量]
$a=b=\gamma$ なら有効損失は $3a/2$ である。フォノンを動かない浴とみなした $2a$ より小さい。熱勾配がフォノンを先に動かす場合は、同じ相互作用により電子にも駆動が加わる。このフォノンドラッグの寄与は、電子だけの式 :eqref[eq-ss-mott] に含まれない。
:::

## 磁気抵抗と Hall 係数

### 二キャリアの伝導度を加える

密度 $n$ の電子、$p$ の正孔と正の移動度 $\mu_e,\mu_h$ をとる。各種の抵抗行列は式 :eqref[eq-ss-drude-hall] を電荷符号ごとに変えればよい。個別の伝導度行列を逆行列から求めて足すと、$\vb j=\begin{pmatrix}a&b\\-b&a\end{pmatrix}\vb E$ に対し
$$
 a=e\left[\frac{n\mu_e}{1+(\mu_eB)^2}+\frac{p\mu_h}{1+(\mu_hB)^2}\right],
 \quad b=eB\left[\frac{-n\mu_e^2}{1+(\mu_eB)^2}
 +\frac{p\mu_h^2}{1+(\mu_hB)^2}\right].
$$
外から固定する電場に対して両方のキャリアの電流が流れるため、加えるのは伝導度である。開回路横方向では $\rho_{xx}=a/(a^2+b^2)$、$\rho_{yx}=b/(a^2+b^2)$ であり、弱磁場の Hall 係数は
$$
 R_H=\lim_{B\to0}\frac{\rho_{yx}}B
 =\frac{p\mu_h^2-n\mu_e^2}{e(p\mu_h+n\mu_e)^2}
 \label{eq-ss-two-carrier-hall}
$$
である。低磁場展開では
$$
 \frac{\rho_{xx}(B)-\rho_{xx}(0)}{\rho_{xx}(0)}
 =\frac{np\mu_e\mu_h(\mu_e+\mu_h)^2}{(n\mu_e+p\mu_h)^2}B^2
 +\order{B^4}.
$$
これは $a=a_0-a_2B^2$、$b=b_1B$ と展開し $(a_0a_2-b_1^2)/a_0^2$ を整理した形である。

:::example[補償と抵抗]
$n=p$、$\mu_e=\mu_h=\mu$ では $b=0$ が全磁場で成り立ち、$a=2en\mu/(1+\mu^2B^2)$ である。従って $\rho_{xx}(B)/\rho_{xx}(0)=1+\mu^2B^2$ となる。Hall 電圧は零でも磁気抵抗は残る。
:::

:::example[少数キャリアで符号が変わる]
$p=n/4$、$\mu_h=3\mu_e$ なら $p\mu_h^2=(9/4)n\mu_e^2$ なので Hall 係数は正である。一方、零磁場伝導度では $p\mu_h=(3/4)n\mu_e$ だから電子が大きく寄与する。
:::

### 全磁場での二キャリア模型

$N_m=n\mu_e+p\mu_h$、$D_B=N_m^2+(p-n)^2\mu_e^2\mu_h^2B^2$ と置く。移動度は磁場に依存しないとする。$a,b$ を共通分母 $D=(1+\mu_e^2B^2)(1+\mu_h^2B^2)$ にまとめると

$$
\begin{aligned}
\frac{Da}{e}&=N_m+(n\mu_h+p\mu_e)\mu_e\mu_hB^2,\\
\frac{Db}{e}&=B[p\mu_h^2-n\mu_e^2+(p-n)\mu_e^2\mu_h^2B^2].
\end{aligned}
$$

二乗を足して整理すると $a^2+b^2=e^2D_B/D$ である。例えば各キャリアの伝導度を複素数で表せば $a+ib=e[n\mu_e/(1+i\mu_e B)+p\mu_h/(1-i\mu_h B)]$ であり、二項の和の分子は $e[N_m+i(p-n)\mu_e\mu_hB]$ なので絶対値二乗が直ちに出る。従って

$$
\begin{aligned}
\rho_{xx}(B)&=\frac{N_m+(n\mu_h+p\mu_e)\mu_e\mu_hB^2}{eD_B},\\
\rho_{yx}(B)&=\frac{B[p\mu_h^2-n\mu_e^2+(p-n)\mu_e^2\mu_h^2B^2]}{eD_B}.
\end{aligned}
\label{eq-ss-transport-two-carrier-full}
$$

零磁場を引くと $(\rho_{xx}(B)-\rho_{xx}(0))/\rho_{xx}(0)$ の分子は $np\mu_e\mu_h(\mu_e+\mu_h)^2B^2$、分母は $D_B$ である。この形は磁気抵抗が非負になることを示し、$B\to0$ で先の二次近似に戻る。$n\ne p$ で両種が存在する場合、大きい $B$ では縦抵抗率は $(n\mu_h+p\mu_e)/[e(p-n)^2\mu_e\mu_h]$ に飽和し、$\rho_{yx}/B\to1/[e(p-n)]$ となる。$n=p$ では $D_B$ の二次項が消え、磁気抵抗は模型の中では $B^2$ で増え続ける。測定上の高磁場にこの結果を使うには、準位の量子化などで密度や散乱が変わらない範囲であることが必要である。

:::example[Hall 電圧が磁場で符号を変える場合]
$p=n/4$、$\mu_h=3\mu_e$ とする。式 :eqref[eq-ss-transport-two-carrier-full] の角括弧を $n\mu_e^2$ で割ると $5/4-(27/4)(\mu_eB)^2$ である。従って $B>0$ では $\mu_eB=\sqrt{5/27}\simeq0.430$ を境に、$\rho_{yx}$ が正から負へ変わる。弱磁場では移動度の大きい正孔の寄与が優勢だが、強磁場では密度差を電子が支配する。
:::

:::example[ほぼ補償された密度]
$p=0.99n$、移動度は両方とも $\mu$ とする。$D_B$ の二つの項が等しくなる尺度は $(n+p)/(|p-n|\mu)=199/\mu$ である。単に $\mu B$ が $1$ を超えた程度では飽和を見込めない。例えば $\mu B=10$ でも $(p-n)^2\mu^4B^2/N_m^2\simeq2.5\times10^{-3}$ であり、縦抵抗率はほぼ二次で増加する。
:::

### 一つの帯でも緩和率が状態に依存する場合

等方放物線帯で、散乱が弾性的だが $\tau(E)$ が変わる場合を考える。各エネルギー殻の正の零磁場伝導度の重みを $\dd\sigma_0(E)=q^2\Xi(E)w(E)\dd E$ とする。移動度の大きさを $u(E)=|q|\tau(E)/m$ とすると、場の強さが寄与する各殻で $u|B|\ll1$ の範囲で

$$
a=\sigma_0[1-\langle u^2\rangle B^2+\order{B^4}],\qquad
b=\sgn(q)\sigma_0[\langle u\rangle B+\order{B^3}]
$$

である。$\langle\cdot\rangle$ は $\dd\sigma_0/\sigma_0$ で平均したものを表す。殻ごとに一定緩和時間の Drude 式を展開して足した結果である。逆行列をとると

$$
\frac{\rho_{xx}(B)}{\rho_{xx}(0)}
=1+[\langle u^2\rangle-\langle u\rangle^2]B^2+\order{B^4}
\label{eq-ss-transport-mobility-variance}
$$

となる。一帯でも各エネルギーの移動度が異なれば分散が残る。温度零で単一の Fermi エネルギーだけが寄与し、同じ面上で $\tau$ が一定ならこの機構は消える。逆に温度を上げて異なるエネルギーが電流を運ぶと、磁気抵抗にその分布幅が反映される。

同じ分布幅は弱磁場の Hall 係数にも効く。質量 $m$ の等方放物線帯に対して、$\tau$ を除いた輸送重み $(m/n)g(E)v(E)^2 w(E)\dd E/3$ を使う。この重みの積分は、先の放物線帯の部分積分で $\tau=1$ と置けば $1$ と分かる。この測度の平均を $\langle\cdot\rangle_D$ と書くと、零磁場の伝導度は $nq^2\langle\tau\rangle_D/m$ である。$b/B$ に入る殻ごとの重みにはさらに符号付きの $q\tau/m$ が掛かるので、$b/B=nq^3\langle\tau^2\rangle_D/m^2$ となる。従って
$$
 R_H=\frac{r_H}{nq},\qquad
 r_H=\frac{\langle\tau^2\rangle_D}{\langle\tau\rangle_D^2}\ge1
$$
である。:term[Hall 因子]{en="Hall factor" yomi="ほーるいんし"} $r_H$ の不等号は $\tau$ の分散が非負であることによる。単一帯でも、測定された $R_H$ を $1/(nq)$ と置くには $\tau$ のエネルギー依存が無視できるかを確かめる。

:::example[同符号キャリアの正の磁気抵抗]
電気伝導への重みが等しい二つの狭い領域で、$u_1=\mu$、$u_2=3\mu$ とする。$\langle u\rangle=2\mu$、$\langle u^2\rangle=5\mu^2$ なので、縦抵抗率は低磁場で $\rho_{xx}(0)(1+\mu^2B^2)$ まで増える。両者が電子なら Hall 係数の符号は負のままである。符号の反対のキャリアが磁気抵抗に必須ではないことが分かる。
:::

### 過去の軌道を積分する方法

一般の周期的な帯で、磁場を保ち電場だけを一次まで扱う。緩和時間 $\tau$ は一定、$f_0$ は一様であるとする。磁場だけの軌道 $\vb k(t')$ は $\hbar\dot{\vb k}=q\vb v\times\vb B$ を満たし、エネルギーは一定なので $w(E)$ も軌道上で一定である。この軌道上の緩和方程式を積分すると

$$
\delta f(\vb k(0))
=q w(E)\int_{-\infty}^0 e^{t'/\tau}\vb v(\vb k(t'))\cdot\vb E\,\dd t'
$$

となる。現在の速度に掛けて全状態で足せば、一様な電場への伝導度は

$$
\sigma_{ij}(\vb B)
=\frac{q^2}{V}\sum_{\vb k}s_{\vb k}w(E_{\vb k})\,
 v_i(\vb k)\int_0^\infty e^{-t/\tau}
 v_j(\vb k(-t;\vb k,\vb B))\,\dd t
\label{eq-ss-transport-orbit-conductivity}
$$

である。これを :term[Chambers の公式]{en="Chambers formula" yomi="ちぇんばーずのこうしき"} という。$s_{\vb k}$ は表示していない内部縮退数であり、個別に状態和をとれば省ける。$\vb B=0$ では速度が一定なので過去の積分が $\tau v_j$ になり、前節のテンソルを得る。

磁場の符号を反転すると、同じ $\vb k$ を出発する軌道の時間方向が反転する。これは軌道の右辺が $\vb B$ に関して奇であるためである。軌道が波数空間の体積を保つことを使い、各 $t$ で式 :eqref[eq-ss-transport-orbit-conductivity] の積分変数を $\vb k$ から $\vb k(-t;\vb k,\vb B)$ に変えると、$i,j$ が交換されて $\vb B$ が $-\vb B$ となる。従ってこの模型で $\sigma_{ij}(\vb B)=\sigma_{ji}(-\vb B)$ が成り立つ。これは磁場反転と電流・電場の方向交換の関係である。

:::example[円運動の速度相関]
単一の等方放物線帯で $\vb B=B\hat{\vb z}$ とする。符号付き $\omega_c=qB/m$ に対し、過去の速度は $v_x(-t)=v_x(0)\cos(\omega_c t)-v_y(0)\sin(\omega_c t)$、$v_y(-t)=v_x(0)\sin(\omega_c t)+v_y(0)\cos(\omega_c t)$ である。軌道は $\dot v_x=\omega_c v_y$、$\dot v_y=-\omega_c v_x$ の解として得る。角度平均で $\langle v_xv_y\rangle=0$ なので対角の積分は $\int_0^\infty e^{-t/\tau}\cos(\omega_c t)\dd t=\tau/(1+\omega_c^2\tau^2)$、$xy$ 成分には $\int_0^\infty e^{-t/\tau}\sin(\omega_c t)\dd t=\omega_c\tau^2/(1+\omega_c^2\tau^2)$ が入る。指数関数 $e^{(-1/\tau+i\omega_c)t}$ を積分して実部と虚部を分ければよい。これにより前の伝導度行列が軌道の平均として再現される。
:::

### 磁場中の熱電応答

式 :eqref[eq-ss-transport-orbit-conductivity] では、軌道上でエネルギーが一定である。したがって同じ速度相関を使い、$w$ に $(E-\mu)^j$ を掛けた状態和を $\mathsf L_j(\vb B)$ と定義できる。ここでは電荷の前因子 $q^2$ を $\mathsf L_j$ に含めない。粒子の駆動に $-(E-\mu)\grad T/T$ を加え、熱流を求めるときは出力側にも $E-\mu$ を掛ければ、式 :eqref[eq-ss-transport-moments] のテンソル形が得られる。軌道変数を入れ替える同じ証明により $\mathsf L_j(\vb B)=\mathsf L_j^{\mathsf T}(-\vb B)$ である。伝導が有限な方向で $\mathsf L_0$ を逆にできるなら、
$$
 \mathsf S(\vb B)=\frac{1}{qT}\mathsf L_0^{-1}(\vb B)\mathsf L_1(\vb B),
 \qquad
 \boldsymbol\Pi(\vb B)=\frac{1}{q}\mathsf L_1(\vb B)\mathsf L_0^{-1}(\vb B)
$$
となる。積を転置すると順序が逆になるので、磁場中の Kelvin 関係は $\boldsymbol\Pi(\vb B)=T\mathsf S^{\mathsf T}(-\vb B)$ である。

例えば面内が等方的な金属で温度勾配を $x$ 方向、磁場を $z$ 方向に取る。開回路で生じる横の電気化学場 $\mathcal E_y=S_{yx}\partial_xT$ が :term[Nernst 効果]{en="Nernst effect" yomi="ねるんすとこうか"} である。単一放物線帯かつ $\tau$ が全エネルギーで一定なら、角度平均した速度相関の磁場依存はエネルギーに依存しない。実際、円運動の周波数 $qB/m$ と積分の減衰時間 $\tau$ がどの殻でも同じである。従って $\mathsf L_0$ と $\mathsf L_1$ は同じ磁場依存の行列に、それぞれ別のスカラーを掛けた形になる。積 $\mathsf L_0^{-1}\mathsf L_1$ は単位行列に比例し、$S_{yx}$ は零になる。緩和時間がエネルギーで変わる場合はこの相殺が一般に失われる。

## 演習問題

::::exercise{level=basic}
衝突後の分布の第一調和が $e^{-t/\tau}$ で減る。初期値の $1/4$ になる時刻を求めよ。
:::solution
$e^{-t/\tau}=1/4$ の対数をとると $t=\tau\log4$ である。平衡との差についての時刻である。
:::
::::

::::exercise{level=basic}
$\sigma=4\times10^7\,\mathrm{S\,m^{-1}}$、$T=250\,\mathrm K$ で Wiedemann–Franz 則が成立する。$L_0=2.44\times10^{-8}\,\mathrm{W\,\Omega\,K^{-2}}$ として熱伝導度を求めよ。ここで $L_0$ は Lorenz 数である。
:::solution
$\kappa=L_0\sigma T=244\,\mathrm{W\,m^{-1}K^{-1}}$ である。格子の熱伝導分を足す前の値である。
:::
::::

::::exercise{level=basic}
$S=-20\,\mathrm{\mu V\,K^{-1}}$、$T=300\,\mathrm K$ のとき $\Pi$ を求めよ。
:::solution
$\Pi=ST=-6000\,\mathrm{\mu V}=-6.0\,\mathrm{mV}$ である。熱流と電流が反対方向になる符号である。
:::
::::

::::exercise{level=basic}
低温のフォノン抵抗率で温度を三倍にしたとき、残留分を引いた抵抗率の比を求めよ。
:::solution
低温の式 :eqref[eq-ss-bloch-gruneisen] は $T^5$ に比例するから、比は $3^5=243$ である。
:::
::::

::::exercise{level=advanced}
前方と後方の散乱のみが各々同じ率 $W$ で起こる。全寿命と運動量緩和時間の比を求めよ。
:::solution
全散乱率は $W+W=2W$ である。運動量の重みは前方 $0$、後方 $2$ なので緩和率も $2W$ になる。よって二つの時間の比は $1$ である。後方散乱だけなら緩和率は全散乱率の二倍となる。
:::
::::

::::exercise{level=advanced}
$|E-\mu|$ に関して $\Xi$ が対称な場合、状態の積分範囲も対称だとして $S=\Pi=0$ を示せ。
:::solution
$w$ は $\mu$ を中心に偶関数であり、$L_1$ の被積分関数にだけ奇関数因子 $E-\mu$ がある。対称な範囲で積分すると零になる。従って $S=L_1/(qTL_0)$、$\Pi=L_1/(qL_0)$ はともに零である。
:::
::::

::::exercise{level=advanced}
$n=p$、$\mu_h=2\mu_e$ の弱磁場 Hall 係数を求めよ。
:::solution
式 :eqref[eq-ss-two-carrier-hall] の分子は $3n\mu_e^2$、分母は $9en^2\mu_e^2$ なので $R_H=1/(3en)$ である。
:::
::::

::::exercise{level=advanced}
二キャリア模型で $\mu_e=\mu_h=\mu$ とする。$B$ が小さいとき $\rho_{xx}(B)/\rho_{xx}(0)$ の $B^2$ 係数を求めよ。
:::solution
$a=e(n+p)\mu/[1+(\mu B)^2]$、$b=e(p-n)\mu^2B/[1+(\mu B)^2]$ である。これを逆にして比を展開すると $1+\mu^2[1-(p-n)^2/(p+n)^2]B^2+\order{B^4}$ となる。従って係数は $4np\mu^2/(p+n)^2$ である。
:::
::::

::::exercise{level=basic}
拡散方程式で波長 $\lambda$ の密度変動が時定数 $t_0$ で減衰する。別の波長 $\lambda/3$ の時定数を求めよ。ただし拡散係数は一定で境界はないとする。
:::solution
波数が $2\pi/\lambda$ からその三倍へ変わる。正弦波の減衰率は $Dk^2$ なので九倍になり、時定数は $t_0/9$ である。両波長とも拡散近似が使える範囲とした。
:::
::::

::::exercise{level=advanced}
輸送分布 $\Xi(E)=cE$、$c>0$ のバンドを考える。$\mu>0$ で低温とし、Lorenz 数の温度二次の補正が正か負か調べよ。
:::solution
式 :eqref[eq-ss-transport-lorenz-correction] で $\Xi''=0$、$\Xi'/\Xi=1/\mu$ である。従って補正は $-(k_B^2/q^2)(k_BT)^2 a_2^2/\mu^2$ で負である。曲率の項が零なので、開回路で熱流が減少する項だけが残る。$\mu(T)$ の温度変化が二次なら、この補正の係数では $\mu(0)$ を使ってよい。
:::
::::

::::exercise{level=advanced}
遮蔽散乱の式 :eqref[eq-ss-transport-forward-ratio] で $b\to0$ の先頭補正を求めよ。
:::solution
$\log(1+b)=b-b^2/2+b^3/3+\order{b^4}$ を使うと $(1+b)\log(1+b)/b=1+b/2-b^2/6+\order{b^3}$ である。従って比は $1-b/3+\order{b^2}$ となる。$q_s$ が大きく $b$ が小さい場合でも、散乱に角度依存が現れ始めると運動量緩和時間は全遷移の寿命より長くなる。
:::
::::

::::exercise{level=advanced}
二キャリア模型で移動度が両方 $\mu$ であり $p=2n$ とする。$B$ を十分大きくしたときの $\rho_{xx}(B)/\rho_{xx}(0)$ と $\rho_{yx}(B)/B$ の極限を求めよ。
:::solution
式 :eqref[eq-ss-transport-two-carrier-full] で等移動度なら、縦抵抗率の比は $(n+p)^2(1+\mu^2B^2)/[(n+p)^2+(p-n)^2\mu^2B^2]$ である。高磁場極限は $(n+p)^2/(p-n)^2=9$ である。Hall 傾きは同じ式の最高次数について $1/[e(p-n)]=1/(en)$ となる。正孔密度が大きいため符号は正である。
:::
::::

## 確認テスト

::::quiz
開回路の電子熱伝導度に $L_1^2/L_0$ が現れるのはどの条件によるか。

- [x] 温度差で流れようとする電流を電気化学場で打ち消す。
- [ ] 熱流を常に零とする。
- [ ] 散乱率をすべて零とする。

:::explanation
$\vb j=0$ から $\vb{\mathcal E}=L_1\grad T/(qTL_0)$ が必要になる。これを熱流へ代入すると、電流を許した場合の $L_2/T$ から $L_1^2/(TL_0)$ が引かれる。
:::
::::

::::quiz{answer="32" tol="1e-6"}
音響フォノン抵抗の $T^5$ 領域で温度を二倍にすると、抵抗のフォノン成分は何倍か。
:::explanation
式 :eqref[eq-ss-bloch-gruneisen] の積分がほぼ定数となる領域なので $2^5=32$ 倍である。
:::
::::

::::quiz
弱磁場の Hall 係数が零になる条件はどれか。

- [x] $p\mu_h^2=n\mu_e^2$
- [ ] 常に $p=n$ だけで十分
- [ ] 常に $p\mu_h=n\mu_e$ だけで十分

:::explanation
弱磁場では各種の非対角伝導度が電荷符号付きの密度と移動度二乗に比例する。式 :eqref[eq-ss-two-carrier-hall] の分子が零となる条件が上の等式である。
:::
::::

::::quiz{answer="2" tol="1e-6"}
完全補償で移動度が両キャリアとも $\mu$ のとき、$\mu B=1$ における $\rho_{xx}(B)/\rho_{xx}(0)$ を答えよ。
:::explanation
両者の Hall 成分が打ち消し合うので $\rho_{xx}$ は対角伝導度の逆数となる。比は $1+\mu^2B^2=2$ である。
:::
::::

::::quiz
二種類の独立な散乱の運動量緩和率を足したとき、抵抗率も単純に足せる条件はどれか。

- [x] 各散乱が同じ駆動モードを一つの率で緩和し、同じ密度と質量の Drude 式を使える。
- [ ] 各散乱の遷移率がエネルギーごとに足せれば、伝導度をどのように平均しても抵抗率は足せる。
- [ ] 個々に計算した抵抗率が正なら十分である。

:::explanation
駆動が両機構の共通の緩和固有モードなら、そのモードの率は加算され、Drude 抵抗率は合計の率に比例する。異なるエネルギー領域で率の違いがあれば、先に各領域の率を足し、逆数を取って伝導度へ寄与させる。最後に全伝導度の逆数を取ることになるため、一般に各抵抗率の和と一致しない。
:::
::::

::::quiz{answer="0.25" tol="1e-6"}
熱電素子の内部抵抗が $R$、起電力の大きさが $V_0$ のとき、最大出力を $V_0^2/R$ を単位として答えよ。
:::explanation
負荷抵抗 $R_L$ に渡る電力は $V_0^2R_L/(R+R_L)^2$ である。微分すると $R_L=R$ で極大となり、端点 $R_L=0,\infty$ で電力は零なので最大である。その値は $V_0^2/(4R)$ であり、指定単位で $1/4=0.25$ となる。
:::
::::

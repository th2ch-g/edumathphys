---
title: 超伝導
description: 磁束排除と熱力学の測定事実を London・Ginzburg–Landau・BCS の各模型と結び付ける。
notation:
  - symbol: 'q_s,m_s,n_s'
    meaning: 凝縮粒子の電荷、質量、数密度
  - symbol: '\lambda,\xi'
    meaning: 磁場侵入長とコヒーレンス長
  - symbol: 'B_c,B_{c1},B_{c2}'
    meaning: 真空中の印加臨界磁束密度
  - symbol: '\Phi_0'
    meaning: 超伝導磁束量子
  - symbol: '\Delta'
    meaning: BCS の対形成ギャップ
---

## この章の目標

- 零抵抗と Meissner 効果を独立な測定として説明する。
- London 方程式から磁場侵入長、熱力学から凝縮エネルギーを求める。
- Ginzburg–Landau の方程式から秩序変数と渦の基本的な尺度を導く。
- 対形成とトンネル結合の模型を使い、ギャップや Josephson 周波数を求める。

電磁気学は SI 単位とする。バルクの臨界印加場は $H_c$ を $\mathrm{A\,m^{-1}}$ で測り、その真空換算 $B_c=\mu_0H_c$ を $\mathrm T$ で表す。凝縮粒子の符号付き電荷 $q_s$ は電子対なら $-2e$、$e>0$ である。

## 超伝導の実験事実

### 電流の減衰と磁場排除を別々に測る

超伝導ではある温度 $T_c$ 以下でごく小さな電流に対する直流抵抗が測定限界まで下がり、閉じた環に電流が長時間残る。また十分弱い磁場の中で温度を下げても、バルク内部から磁束が排除される。後者を :term[Meissner 効果]{en="Meissner effect" yomi="まいすなーこうか"} という。巨大な伝導度をもつ古典模型だけで $\vb E=0$ と置けば、Faraday の法則 $\partial_t\vb B=-\rot\vb E$ は初期磁場を保存する。この場合、無磁場で冷やしてから場を加える履歴と、場を先に加えて冷やす履歴が異なる。平衡の超伝導では磁束排除が外場下での熱力学的状態として現れる。

比熱には転移で異常があり、トンネル分光では一粒子の励起に有限のギャップが観測される場合がある。ギャップの節をもつ超伝導体もあるため、有限の最小ギャップを普遍的な定義にはできない。低温で残る磁束は渦のピン止めなど非平衡の効果を含み得る。一定の電場が生じる境界電流も、対を壊す限界と渦を動かす限界を区別する。

:::example[持続電流の時定数]
古典的な環の自己インダクタンスを $L$、抵抗を $R$ とすれば $L\dot I+RI=0$、解は $I(t)=I(0)e^{-Rt/L}$ である。$L=1\,\mathrm{\mu H}$ で $10^5\,\mathrm s$ の間に電流低下が $1\%$ 以下なら、この回路模型で $R\le-L\log(0.99)/(10^5\,\mathrm s)\simeq1.0\times10^{-13}\,\mathrm\Omega$ である。
:::

:::example[同位体の比較]
ばね定数が同じでイオン質量が倍になると調和角振動数は $1/\sqrt2$ になる。もし対形成のエネルギー窓だけがこの振動数に比例し、無次元結合が同じなら $T_c$ も $1/\sqrt2$ になる。実際に質量だけを変えたときの $-\partial\log T_c/\partial\log M$ を同位体指数と呼ぶ。
:::

### 試料形状と測定条件

試料の磁場応答を読むには形も指定する。例えば楕円体の主軸方向で反磁場係数を $N_d$ とする。外場 $H_a$ と内部の巨視的な $H_{\rm int}$ は $H_{\rm int}=H_a-N_dM$ の関係にある。磁束を排除した領域では $B=\mu_0(H_{\rm int}+M)=0$ なので $M=-H_a/(1-N_d)$ である。長い針状試料で軸を場にそろえると $N_d$ が零に近く、$M\simeq-H_a$ となる。球では一様に磁化したときの磁気スカラーポテンシャルを求めると $N_d=1/3$ になる。半径 $R$ の球で $\vb M=M\hat{\vb z}$ とし、内部を $\phi_{\rm in}=c r\cos\theta$、外部を $\phi_{\rm out}=cR^3\cos\theta/r^2$ と置けば両方が Laplace 方程式を満たし表面で連続である。磁場は $\vb H_d=-\grad\phi$ で、法線の差 $H_{{\rm out},n}-H_{{\rm in},n}=3c\cos\theta$ は磁束密度の連続性により $M\cos\theta$ と等しい。従って $c=M/3$、$\vb H_d=-\vb M/3$ である。完全な磁束排除を示す球の磁化は $M=-3H_a/2$ となる。ここでは磁場侵入のある表面層が球半径に比べ薄いとした。

二つの温度プロトコルは、試料に束縛された磁束の検出にも使える。零磁場で冷却してから印加する実験は、小さい磁場の変化に対する遮蔽を測る。印加磁場のまま冷却する実験では、磁束が排除されるか残るかを測る。第一種でも試料形状による局所的な磁場の集中で正常領域と超伝導領域が共存することがある。第二種では後述する渦が自由エネルギーの安定な構成要素となる磁場範囲がある。両者は起源と磁束分布の尺度が異なる。例えば渦がピン止めされる第二種で冷却経路による差があっても、測定された全磁化が平衡値を示すとは限らない。

:::example[磁束排除の負の磁化]
完全に排除する体積 $10^{-7}\,\mathrm{m^3}$ の球に $H_a=100\,\mathrm{A\,m^{-1}}$ を加えると、磁気モーメントは $(-3H_a/2)V=-1.5\times10^{-5}\,\mathrm{A\,m^2}$ になる。軸方向の反磁場係数をほぼ零にした同体積の長い試料では、同じ条件で $-1.0\times10^{-5}\,\mathrm{A\,m^2}$ に近い。この差は形状による境界条件から出る。
:::

### 電圧で見る励起のしきい値

トンネル障壁の一方に常伝導金属、他方に完全にギャップの開いた超伝導体を置き、準粒子の弾性トンネルを考える。温度零で正常側の Fermi 準位から超伝導側へ電子を注入するには、接合で結合できる電子追加スペクトルの最低エネルギー $\Delta>0$ を供給する必要がある。従って理想的なしきい値は $\Delta/e$ である。両側が同じ等方ギャップ $\Delta$ の超伝導体なら、片方から電子を除く側も準粒子を作る費用が要る。従ってしきい値は $2\Delta/e$ に移る。有限温度や準位の幅はしきい値を丸める。バイアス零でも起こる対のトンネルは後述の Josephson 電流である。

:::example[二種類のしきい値]
$\Delta=0.50\,\mathrm{meV}$ とすれば、正常・超伝導接合では $0.50\,\mathrm{mV}$、同じ二超伝導体では $1.00\,\mathrm{mV}$ である。原点をまたいで対称なしきい値を観測しても、そのエネルギーを $\Delta$ と解釈するには、接合の両側がどの状態にあるかを確認する必要がある。
:::

## London 方程式と Meissner 効果

### 電流と位相

凝縮粒子が位相 $\theta$ を共有し、緩やかな変化では運動量が $m_s\vb v_s=\hbar\grad\theta-q_s\vb A$ となると仮定する。位相のない通常伝導との違いは、電流が $\vb j_s=n_s q_s\vb v_s$ と共に位相の勾配で拘束されることである。$n_s$ を一定、渦がなく $\rot\grad\theta=0$ とすると
$$
 \rot\vb j_s=-\frac{n_s q_s^2}{m_s}\vb B,\qquad
 \lambda^2=\frac{m_s}{\mu_0 n_s q_s^2}
 \label{eq-ss-london-supercurrent}
$$
を得る。粒子を電場で加速する $m_s\partial_t\vb v_s=q_s\vb E$ は、一定密度の理想化で $\partial_t\vb j_s=n_s q_s^2\vb E/m_s$ を与える。この回転を取って Faraday の式と合わせると、$\partial_t(\rot\vb j_s+n_s q_s^2\vb B/m_s)=0$ である。この関係だけでは括弧内の初期値が任意である。平衡の渦なし位相の仮定を使うと、電流の回転は式 :eqref[eq-ss-london-supercurrent] の値に固定される。

静磁場の Ampère 則 $\rot\vb B=\mu_0\vb j_s$ の回転を取り、$\divg\vb B=0$ を使うと $\nabla^2\vb B=\vb B/\lambda^2$ となる。半無限 $x>0$ で面に平行な $B_z(0)=B_a$ に対し、深部で有界な解は $B_z(x)=B_a e^{-x/\lambda}$ である。$\lambda$ は :term[磁場侵入長]{en="magnetic penetration depth" yomi="じばしんにゅうちょう"} であり、電流は $j_y=-\partial_x B_z/\mu_0$ となる。厚さ $2d$ の板で両端に同じ $B_a$ を課せば $B_z(x)=B_a\cosh(x/\lambda)/\cosh(d/\lambda)$ である。

:::example[磁場が残る割合]
$x=3\lambda$ で磁場は表面値の $e^{-3}\simeq0.050$ である。電流を深さについて積分すると $B_a/\mu_0$ となり、$\lambda$ を変えても全表面電流は一定である。
:::

:::example[薄膜では排除が完全でない]
厚さ $2d$ の板の中心では $B(0)/B_a=1/\cosh(d/\lambda)$ である。$d=\lambda\log2$ なら比は $1/[(2+1/2)/2]=4/5$ であり、半無限試料の指数解を両面にそのまま適用できない。
:::

### 時間依存応答と超流動成分

第 1 London 方程式に時間因子 $e^{-i\omega t}$ を入れると、零でない周波数で $\sigma_s(\omega)=in_s q_s^2/(m_s\omega)$ となる。因果性を課すには、時間積分の核 $\Theta(t)n_s q_s^2/m_s$ を $e^{-\eta t}$ で減衰させ $\eta\downarrow0$ とする。Fourier 変換は $n_s q_s^2/[m_s(\eta-i\omega)]$ であり、実部 $\eta/(\eta^2+\omega^2)$ の全積分が $\pi$ である。虚部は $\omega/(\eta^2+\omega^2)$ であり、原点の左右で符号が逆なので試験関数の偶成分を落とした主値に収束する。従って
$$
 \sigma_s(\omega)=\frac{n_s q_s^2}{m_s}
 \left[\pi\delta(\omega)+i\,\mathop{\rm P}\frac1\omega\right]
 \label{eq-ss-london-optical}
$$
と分布の意味で書ける。虚部が $1/\omega$ であることと零周波数の実部の重みは同じ係数を持つ。非零の周波数で実部の広い吸収が見える場合には、熱励起された準粒子が並行に流れるなどの効果がある。単純な二流体模型では $\vb j=\vb j_s+\vb j_n$ とし、$\vb j_n$ に Drude 型の緩和を与える。波の分散では式 :eqref[eq-ss-drude-permittivity] の導き方と同じく、$\epsilon_r=\epsilon_{\rm bg}+i\sigma/(\epsilon_0\omega)$ を介して $\sigma_s$ が入る。$k^2=\epsilon_r\omega^2/c^2$ に純粋な超流動成分を代入すれば $k^2\simeq-\lambda^{-2}$ となり、低周波の波は長さ $\lambda$ で減衰する。

### 円柱内の London 解

長い円柱の軸に沿って静磁場を入れ、半径を $R$ とする。軸対称の方程式は $(rB_z')'/r=B_z/\lambda^2$ である。$r=0$ で正則な解を $I_0(r/\lambda)$ と書く。級数 $\sum_{j\ge0}(r/(2\lambda))^{2j}/(j!)^2$ を微分すると、$j$ 項に $4j^2/r^2$ が掛かり、前の $j-1$ 項の $\lambda^{-2}$ 倍となるので方程式を満たす。$B_z(R)=B_a$ より $B_z(r)/B_a=I_0(r/\lambda)/I_0(R/\lambda)$ である。半径が小さい場合にはこの比を二次まで展開して $1+(r^2-R^2)/(4\lambda^2)$ となる。断面平均は円環の重み $2r\dd r/R^2$ を使うため、$1-R^2/(8\lambda^2)$ である。平板と円柱では、同じ London 方程式でも形状に合わせて境界を指定する必要がある。

:::example[超流動密度を変えたとき]
凝縮粒子の有効質量と電荷を固定すると、$n_s$ が四分の一になれば $\lambda$ は二倍になる。$x$ を旧侵入長の二倍の位置に固定すると、同じ表面磁場で内部値は $e^{-2}$ 倍から $e^{-1}$ 倍に変わる。式 :eqref[eq-ss-london-optical] の虚部の係数は四分の一になるため、交流応答でも密度変化を確かめられる。
:::

:::example[表面電流の運動エネルギー]
平板の $x>0$ で $B_z=B_a e^{-x/\lambda}$ なら $j_y=B_a e^{-x/\lambda}/(\mu_0\lambda)$ である。流れの運動エネルギー密度 $m_s j_y^2/(2 n_s q_s^2)=\mu_0\lambda^2j_y^2/2$ は、この解では磁場エネルギー密度 $B_z^2/(2\mu_0)$ に等しい。両方を深さ方向で積分すると単位面積当たり $B_a^2\lambda/(2\mu_0)$ になる。これは局所場を保持した配置のエネルギーであり、外部電源の磁気仕事を含む固定印加場の Gibbs エネルギーと比較するときは境界をそろえる。
:::

### 細い導線の運動インダクタンス

断面積 $A$、長さ $\ell$ の導線で、横の寸法が侵入長より小さく電流密度 $j_s\simeq I/A$ が一様だとする。流れの運動エネルギーを全体で積分すると

$$
U_{\rm kin}=\frac{\mu_0\lambda^2\ell}{2A}I^2
=\frac12 L_{\rm kin} I^2,\qquad
L_{\rm kin}=\frac{\mu_0\lambda^2\ell}{A}
$$

となる。これを :term[運動インダクタンス]{en="kinetic inductance" yomi="うんどういんだくたんす"} という。外側を含む磁場のエネルギーで決まる幾何学的なインダクタンスとは、エネルギーの由来が異なる。

同じ線を細い環にし、幾何学的インダクタンスを $L_g$ とする。後の磁束の量子化の節で導く通り、環の向きに合わせて整数を取ると $\Phi_{\rm ext}+(L_g+L_{\rm kin})I=n\Phi_0$ となる。運動の項は $\mu_0\lambda^2(I/A)\ell=L_{\rm kin}I$ に対応する。磁束を与えたときの電流は巻き数ごとに決まり、自己エネルギー $(L_g+L_{\rm kin}) I^2/2$ は外部磁束に対して二次曲線となる。平衡でこの寄与が最小の巻き数を選べるなら、半整数の $\Phi_{\rm ext}/\Phi_0$ で二つの曲線が交差する。実際に巻き数が変わるには振幅が一時的に小さくなる過程を経るため、平衡の枝に直ちに移るとは限らない。

:::example[細線での運動エネルギーの尺度]
$\lambda=200\,\mathrm{nm}$、$\ell=10\,\mathrm{\mu m}$、$A=(100\,\mathrm{nm})^2$ とすると $L_{\rm kin}\simeq50\,\mathrm{pH}$ である。環の $L_g$ が $100\,\mathrm{pH}$ なら一つ隣の巻き数の電流との差は $\Phi_0/(L_g+L_{\rm kin})\simeq14\,\mathrm{\mu A}$ となる。この電流が実現するには、各枝で電流密度が対形成を壊さない範囲にあることも必要である。
:::

### 清浄試料における局所近似の範囲

London 方程式は、電流がその場所のゲージ不変な運動量だけで決まる近似である。清浄で等方的な超伝導体では、Cooper 対の広がりの尺度 $\xi_0$ の範囲にある場も電流に寄与する。この尺度の見積もりは後の BCS 理論で行う。例えば横応答の Fourier 表示で、横成分を $\vb A_\perp(\vb k)$ とすると、非局所模型では $\vb j_\perp(\vb k)=-K_\perp(k)\vb A_\perp(\vb k)/\mu_0$ と書く。$K_\perp(0)=\lambda_L^{-2}$ とし、ここで $\lambda_L$ を局所 London 長と呼ぶ。$k\xi_0\ll1$ の範囲で $K_\perp(k)\simeq K_\perp(0)$ と置いたのが局所模型である。従って $\lambda_L\gg \xi_0$ なら、局所解の主な波数 $k\sim1/\lambda_L$ はこの範囲に入る。

逆に清浄極限で $\xi_0\gg\lambda_L$ のとき、Pippard の模型は $k\xi_0\gg1$ で $K_\perp(k)$ が $1/(\lambda_L^2 k\xi_0)$ の定数倍に下がると与える。これは横の電流応答を Fermi 面の方向にわたって平均する結果であり、実空間の境界条件と正確な係数には非局所の応答核が必要である :cite[tinkham-sc]。Ampère 則で現れる $k^2$ がこの応答と同程度になる波数は $k\sim(\lambda_L^2\xi_0)^{-1/3}$ である。従って特徴的な侵入の尺度は $(\lambda_L^2\xi_0)^{1/3}$ となり、一般の磁場分布を単一の指数で表すことはできない。

:::example[局所近似が崩れる尺度]
$\lambda_L=30\,\mathrm{nm}$、$\xi_0=2.4\,\mathrm{\mu m}$ なら、局所解を仮定した $k\xi_0$ は約 $80$ である。非局所模型の尺度は $(30^2\cdot2400)^{1/3}\,\mathrm{nm}\simeq129\,\mathrm{nm}$ で、その波数でも $k\xi_0$ は $1$ より十分大きい。これは長さの桁の比較であり、半無限面の厳密な解には表面散乱の条件も入る。
:::

## 超伝導の熱力学と臨界磁場

### 磁場下の安定性

反磁場の影響を省ける長い試料で、第一種の平衡を考える。零磁場での単位体積の自由エネルギーを $f_s,f_n$ とし、正常相の磁化は無視する。固定印加場での磁気 Gibbs エネルギーの物質部分は $\partial g/\partial H=-\mu_0M$ を満たす。完全排除では $M_s=-H$ なので、零場から積分すると $g_s(H)-g_n(H)=f_s-f_n+\mu_0H^2/2$ である。従って両相が等しくなる場では
$$
 f_n-f_s=\frac{B_c^2}{2\mu_0}.
 \label{eq-ss-condensation-energy}
$$
これは零場で得られる凝縮エネルギー密度である。表面近くの侵入長の領域の寄与はバルク体積に対して無視した。

$S=-\partial_T f$ とすると $S_n-S_s=-(B_c/\mu_0)\dd B_c/\dd T$ である。$B_c$ が温度と共に減れば正常相のエントロピーの方が大きい。有限磁場で相境界を横切る潜熱密度は $T(S_n-S_s)$ である。零場転移では $B_c(T_c)=0$、傾きが有限ならエントロピーは連続であり、
$$
 C_s-C_n=\frac{T_c}{\mu_0}[B_c'(T_c)]^2
$$
という比熱跳びを得る。$C=-T\partial_T^2 f$ を用いた。

:::example[凝縮エネルギーの桁]
$B_c(0)=0.10\,\mathrm T$ なら $B_c^2/(2\mu_0)\simeq4.0\times10^3\,\mathrm{J\,m^{-3}}$ である。これは一原子ごとの化学結合エネルギーではなく、電子系の二状態の自由エネルギー差である。
:::

:::example[臨界磁場の近似曲線]
$B_c(T)=B_0[1-(T/T_c)^2]$ と近似すれば、$T=T_c/2$ では $B_c=3B_0/4$、凝縮エネルギー密度は零温値の $9/16$ である。比熱跳びは $4B_0^2/(\mu_0 T_c)$ になる。
:::

### 相境界に沿う微分

固定印加場における Gibbs 密度の微分は $\dd g=-s\dd T-\mu_0M\dd H$ である。正常相と超伝導相が同じ $g$ を持つ曲線 $H_c(T)$ では
$$
 s_n-s_s=-\mu_0(M_n-M_s)\dv{H_c}{T}
$$
である。完全排除と $M_n\simeq0$ を入れると $M_n-M_s=H_c$ なので式 :eqref[eq-ss-condensation-energy] からのエントロピー差と一致する。これは磁場で駆動する熱力学的な Clapeyron の関係である。$0<T<T_c$ で $H_c'(T)<0$ のとき、一定温度で磁場を強くして正常相に移ると、その温度を保つために潜熱を吸収する。零磁場の $T=T_c$ では $M_n-M_s=0$ になるため右辺が零になる。実験で場を掃引する際は磁気熱量の可逆な寄与と、渦運動やヒステリシスによる発熱を区別する。

:::example[有限温度の潜熱]
$B_c=B_0(1-T^2/T_c^2)$ として $u=T/T_c$ と置けば、潜熱密度は $T[-B_c B_c'/\mu_0]=2B_0^2u^2(1-u^2)/\mu_0$ である。$u=1/2$ では $3B_0^2/(8\mu_0)$ である。$u^2(1-u^2)$ が $u^2=1/2$ で最大になるので、$u=1/\sqrt2$ で最大値 $B_0^2/(2\mu_0)$ を取る。
:::

:::example[低温の正常電子比熱の反映]
極低温で正常相の電子比熱密度が $\gamma T$、超伝導相の電子比熱が温度に対し指数的に小さいと仮定する。第三法則で $s_n-s_s\simeq\int_0^T\gamma T'\dd T'/T'=\gamma T$ である。式 :eqref[eq-ss-condensation-energy] の微分により $B_c(T)^2\simeq B_c(0)^2-\mu_0\gamma T^2$ となる。経験式の低温二次係数を、比熱から得た $\gamma$ と照合できる。
:::

### 反磁場を含む第一種の中間状態

楕円体の主軸に印加した場を $H_a$、反磁場係数を $0<N_d<1$ とする。表面層と界面のエネルギーを体積項に比べて無視できる極限で、正常相の体積割合を $p$ とする。両相が共存する間は、平面界面に沿う正常領域の場は $H_c$ であり、超伝導領域の磁束は零である。従って平均磁束密度は $pB_c$ となる。正常相の磁化を無視した巨視的な構成式では $H_{\rm int}=H_c$ と置けるので、$\bar B=\mu_0(H_{\rm int}+M)$ から平均磁化は $M=-(1-p)H_c$ となる。楕円体の反磁場式へ代入して

$
 H_a=H_c[1-N_d(1-p)],\qquad
 p=\frac{H_a/H_c-(1-N_d)}{N_d}
$

を得る。従って理想的な中間状態の範囲は $(1-N_d)H_c<H_a<H_c$ である。この計算では正常領域と超伝導領域の細かい形を平均した。個々の領域の形は、正の界面エネルギーと試料表面付近の磁場が決める。

:::example[球で共存が始まる場]
$N_d=1/3$ なら $H_a=2H_c/3$ で正常領域が生じ始める。$H_a=5H_c/6$ では $p=1/2$、$M=-H_c/2$ である。これに対し反磁場がない極限では共存する印加場の幅 $N_dH_c$ が零になる。
:::

## Ginzburg–Landau 理論

### 複素秩序変数の汎関数

転移温度の近傍で緩やかに変化する複素秩序変数 $\psi$ を取り、$|\psi|^2$ を模型内の凝縮粒子密度とする。ゲージ変換 $\vb A\mapsto\vb A+\grad\chi$ と $\psi\mapsto e^{iq_s\chi/\hbar}\psi$ に共変な微分を $P=-i\hbar\grad-q_s\vb A$ と書く。:term[Ginzburg–Landau 理論]{en="Ginzburg–Landau theory" yomi="ぎんずぶるくらんだうりろん"} の自由エネルギーは零印加場を基準に
$$
 F=\int\left[a|\psi|^2+\frac b2|\psi|^4+
 \frac{|P\psi|^2}{2m_s}+\frac{B^2}{2\mu_0}\right]\dd^3r,\quad
 a=a_0(T-T_c),\ b>0 .
$$
一様な場合の極小は $a<0$ で $|\psi_0|^2=-a/b$、エネルギー密度差は $-a^2/(2b)$ である。$\psi^*$ を変分し、境界項が消える条件で部分積分すると
$$
 (a+b|\psi|^2)\psi+\frac{P^2}{2m_s}\psi=0 .
$$
$\vb A$ を変分すると磁場エネルギーは $\rot\vb B/\mu_0$、運動項は $-(q_s/m_s)\Re(\psi^* P\psi)$ を与える。従って
$$
 \vb j_s=\frac{q_s}{m_s}\Re(\psi^*P\psi)
 =\frac{q_s}{m_s}|\psi|^2(\hbar\grad\theta-q_s\vb A).
$$
一定密度なら式 :eqref[eq-ss-london-supercurrent] に戻る。磁場なしで一次元の変化を $\psi=\psi_0 f$ とし尺度を
$$
 \xi=\frac{\hbar}{\sqrt{2m_s|a|}},\qquad
 \lambda=\sqrt{\frac{m_s b}{\mu_0q_s^2|a|}}
 \label{eq-ss-gl-lengths}
$$
とすれば $-\xi^2 f''-f+f^3=0$ である。$\xi$ は秩序変数が空間変化する尺度である。例えば $x=0$ に完全な対破壊境界 $f(0)=0$ を置けば、$f(x)=\tanh[x/(\sqrt2\xi)]$ は直接微分で方程式と $f(+\infty)=1$ を満たす。

:::example[温度依存]
温度差 $T_c-T$ が四分の一になると $|\psi_0|^2$ は四分の一、$\lambda$ と $\xi$ は二倍である。比 $\kappa=\lambda/\xi$ はこの次数の係数の範囲で変わらない。
:::

:::example[凝縮エネルギーの比較]
式 :eqref[eq-ss-condensation-energy] と $a^2/(2b)$ を等置すると $B_c=|a|\sqrt{\mu_0/b}$ になる。式 :eqref[eq-ss-gl-lengths] の積から同じ式を $B_c=\hbar/(\sqrt2|q_s|\lambda\xi)$ と書ける。
:::

### 流れがあるときの振幅と最大電流

細い一様なチャネルで電流の自己磁場を無視し、ゲージ不変な波数を $Q=\partial_x\theta-q_s A_x/\hbar$ とする。振幅を一定として運動エネルギーを加えると、二次係数は $a+\hbar^2Q^2/(2m_s)$ になる。従って非零解は
$$
 |\psi(Q)|^2=|\psi_0|^2(1-\xi^2Q^2),\qquad
 j(Q)=\frac{q_s\hbar|\psi_0|^2}{m_s}Q(1-\xi^2Q^2)
$$
である。最初の式は $|\xi Q|<1$ の範囲で使う。$|j|$ は $|\xi Q|=1/\sqrt3$ で最大であり、電流密度の大きさは
$$
 j_{\rm d}=\frac{2|q_s|\hbar|\psi_0|^2}{3\sqrt3\,m_s\xi}
 =\frac{2\sqrt2}{3\sqrt3}\frac{B_c}{\mu_0\lambda}
 \label{eq-ss-depairing-current}
$$
となる。$|\xi Q|=1$ まで大きくすると振幅も電流も零になる。最大の臨界電流を計算する際に密度を一定とみなすと、この振幅の減少を失う。実際には厚さ、欠陥、渦の侵入によりさらに低い電流で散逸が見えることがある。

### 試料表面での振幅

試料表面に局所的な自由エネルギー $g_s|\psi|^2$ を置く場合、$\psi^*$ による変分の境界項から、$\vb n$ を試料外向き法線として $\vb n\cdot(\grad-iq_s\vb A/\hbar)\psi=-(2m_sg_s/\hbar^2)\psi$ となる。運動エネルギーの部分積分で境界に $\hbar^2\vb n\cdot(\grad-iq_s\vb A/\hbar)\psi/(2m_s)$ が残るためである。絶縁境界で $g_s=0$ なら法線電流は零である。$g_s$ が大きく正なら $\psi$ が表面で小さくなる。局所的な界面でどれだけ凝縮が抑えられるかを表す条件である。

:::example[位相差と最大電流]
$\xi=20\,\mathrm{nm}$ の細い均一なチャネルで、自己磁場なしの最大電流に対応する $|Q|$ は約 $2.9\times10^7\,\mathrm{m^{-1}}$ である。振幅二乗の比 $|\psi(Q)|^2/|\psi_0|^2$ は $2/3$ まで減っている。チャネル長を $200\,\mathrm{nm}$ とすると、位相勾配の長さ積分 $|Q|L$ は約 $5.8$ であり、開いたチャネルでは整数倍に限られない。
:::

:::example[弱い対破壊と表面振幅]
磁場零で $x>0$ に試料を置き、$\vb n=-\hat{\vb x}$、$g_s>0$ とする。無次元振幅の候補 $f(x)=\tanh[(x+x_0)/(\sqrt2\xi)]$ はバルクの微分方程式を満たし、表面で $f'(0)=f(0)/b_s$、$b_s=\hbar^2/(2m_sg_s)$ である。$y=f(0)$ とすれば $(1-y^2)/(\sqrt2\xi)=y/b_s$ であり、$0<y<1$ の根を取る。$b_s=\sqrt2\xi$ なら $1-y^2=y$ により $y=(\sqrt5-1)/2$ である。$b_s\to\infty$ で $y\to1$、$b_s\to0$ で $y\to0$ となる。
:::

## 磁束の量子化と第二種超伝導体

### 位相が一周する条件

:term[磁束量子]{en="flux quantum" yomi="じそくりょうし"} は $\Phi_0=h/|q_s|=h/(2e)$ であり、約 $2.07\times10^{-15}\,\mathrm{Wb}$ である。導出するには、$\psi$ が零でない閉曲線 $C$ に沿い位相の増分を数える。一価性から $\oint_C\grad\theta\cdot\dd\vb l=2\pi n$ である。電流の式を積分すると
$$
 \Phi+\frac{m_s}{n_s q_s^2}\oint_C\vb j_s\cdot\dd\vb l=\frac{h}{q_s}n .
 \label{eq-ss-fluxoid}
$$
$n_s$ は積分路上で一定とした。左辺は磁束だけでなく超電流の循環を含む。十分厚い超伝導領域で $\vb j_s\simeq0$ の経路をとれば、磁束自体が $\Phi_0$ の整数倍となる。渦の中心で位相を一周させるには振幅が零になる芯があり、その半径は $\xi$ の桁、磁場の広がりは $\lambda$ の桁である。

### London 近似における渦の磁場

直線の渦を $z$ 軸上に置き、正の磁束 $\Phi_0$ を通す。芯の外で密度を定数とし、細い芯を線源に置き換える。軸を周回する位相勾配の積分により、分布として $\rot\grad\theta=2\pi n\delta^{(2)}(\vb r_\perp)\hat{\vb z}$ である。実際、軸を含まない積分領域では回転は零であり、軸を一度囲む領域の面積分は Stokes の式で $2\pi n$ になる。磁束の向きを正としたので $hn/q_s=\Phi_0$ とする。Ampère 則を用いて電流を除くと

$$
(1-\lambda^2\nabla_\perp^2)B_z(\vb r_\perp)=\Phi_0\delta^{(2)}(\vb r_\perp)
\label{eq-ss-london-vortex-source}
$$

となる。例えば左辺を大きな円で積分すると、$\lambda^2$ の項はその円周での電流循環に転じ、式 :eqref[eq-ss-fluxoid] と同じ条件になる。

二次元の Fourier 変換では $B_z(\vb k)=\Phi_0/(1+\lambda^2k^2)$ である。分母の逆数を $\int_0^\infty e^{-t(1+\lambda^2k^2)}\dd t$ と書くと、各 $t$ の逆変換は二方向の Gaussian 積分に分かれ、$(4\pi\lambda^2t)^{-1}e^{-r^2/(4\lambda^2t)}$ を生む。従って $u=r/\lambda>0$ について

$$
\begin{aligned}
B_z(r)&=\frac{\Phi_0}{2\pi\lambda^2}K_0(u),\\
K_0(u)&=\frac12\int_0^\infty e^{-t-u^2/(4t)}\frac{\dd t}{t}
       =\int_0^\infty e^{-u\cosh v}\dd v .
\end{aligned}
\label{eq-ss-london-vortex-field}
$$

最後は $t=(u/2)e^v$ と置き、$v$ の偶性を使った。$K_0$ は零次の変形 Bessel 関数である。元の Fourier 積分を分布として定義するには、急減衰する試験関数との積から始めればよい。上の表示で $r>0$ の積分は上下の端で収束する。

$K_1(u)=-K_0'(u)=\int_0^\infty \cosh v\,e^{-u\cosh v}\dd v$ と置く。$y=u\sinh v$ を用いると $uK_1(u)=\int_0^\infty e^{-\sqrt{u^2+y^2}}\dd y$ となり、$e^{-u}\le uK_1(u)\le1$ である。従って $u\downarrow0$ では $K_1(u)=1/u+\order{1}$、これを積分して $K_0(u)=-\log u+\order{1}$ である。逆に $u\to\infty$ では $v=z/\sqrt u$ と置き、$u(\cosh(z/\sqrt u)-1)\to z^2/2$ と $\cosh v\ge1+v^2/2$ を使うと $K_0(u)\sim e^{-u}\sqrt{\pi/(2u)}$ となる。小さい $u$ の対数発散は芯を点にしたためであり、実際には $r\lesssim\xi$ で London 近似が破れる。

磁場の積分は $2\pi\int_0^\infty r B_z(r)\dd r=\Phi_0$ である。積分表示で $\int_0^\infty u e^{-u\cosh v}\dd u=1/\cosh^2 v$ なので、$\int_0^\infty uK_0(u)\dd u=\int_0^\infty\dd v/\cosh^2v=1$ となる。回転の式から渦の周りの電流密度は $j_\phi=-B_z'/\mu_0=\Phi_0 K_1(u)/(2\pi\mu_0\lambda^3)$ であり、$\xi\ll r\ll\lambda$ では $1/r$ に比例して流れる。

:::example[二本の渦が及ぼす力]
London の磁場と運動のエネルギー密度を合わせると $(B^2+\lambda^2|\grad B|^2)/(2\mu_0)$ である。距離 $R\gg\xi$ にある二本の渦を重ね、単独の自己エネルギーを引く。交差項の積分で一回部分積分し、片方に式 :eqref[eq-ss-london-vortex-source] を作用させると、単位長さ当たりの相互作用は $\Phi_0^2 K_0(R/\lambda)/(2\pi\mu_0\lambda^2)$ である。同符号の磁束なら正であり、$K_0'=-K_1<0$ なので引き離す向きに力が働く。逆符号の磁束では交差項の符号が反転する。例えば距離が数 $\lambda$ 以上になると、相互作用の大きさは距離に対し指数的に小さくなる。
:::

芯の半径を $\xi$ で打ち切った場合に、London 領域だけのエネルギーを計算することもできる。遠方で $K_0''+K_0'/u=K_0$ なので $(uK_1)'=-uK_0$ である。従って $(uK_0K_1)'=-u(K_0^2+K_1^2)$ となる。軸から $\xi$ より外側の円環を積分すると

$$
\mathcal U_{\rm out}
=\frac{\Phi_0^2}{4\pi\mu_0\lambda^2}
\frac{\xi}{\lambda} K_0(\xi/\lambda)K_1(\xi/\lambda)
=\frac{\Phi_0^2}{4\pi\mu_0\lambda^2}
 [\log(\lambda/\xi)+\order{1}]
$$

である。下部臨界場に現れる対数は、この外側領域だけで求められる。定数項まで求めるには、振幅が変化する芯のエネルギーと接続しなければならない。

### 芯の近くで変わる振幅

渦の芯まで含めるため、円筒座標の方位角を $\phi$ とし $\psi=\psi_0 f(r)e^{in\phi}$、$\vb A=A_\phi(r)\hat{\vb\phi}$ と置く。GL の運動項では方位角の微分が $n/r$ を与える。従って振幅の式と電流の式は

$$
\begin{aligned}
&-\xi^2(f''+f'/r)+\xi^2(n/r-q_s A_\phi/\hbar)^2 f-f+f^3=0,\\
&j_\phi=\frac{q_s\hbar|\psi_0|^2}{m_s}f^2(n/r-q_s A_\phi/\hbar),\\
&B_z=(r A_\phi)'/r,\qquad -B_z'=\mu_0 j_\phi .
\end{aligned}
$$

磁場が軸上で有界なら、特異なゲージを選ばない限り $A_\phi=\order r$ である。小さい $r$ で微分項と $n^2/r^2$ の項を釣り合わせて $f\propto r^p$ を入れると $p^2=n^2$ になる。振幅が有界で運動エネルギーが有限になる方を選ぶと $f\propto r^{|n|}$ である。例えば $|n|=1$ では $f^2\propto r^2$ なので、電流は軸近くで $r$ に比例し零に近づく。芯の外での $1/r$ の挙動をそのまま内側へ延ばすと逆の結果になる。

:::example[磁束二量子の渦の外側エネルギー]
芯の細い模型で一つの芯に二量子を通すと式 :eqref[eq-ss-london-vortex-source] の源は $2\Phi_0$ となる。磁場と電流は二倍なので、同じ打切りの外側エネルギーは一量子の四倍である。一方、十分離れた二本の一量子渦なら自己エネルギーは概ね二倍で、相互作用は指数的に小さい。大きい $\kappa$ の対数が支配する場合に一量子ずつの渦が有利になる理由である。芯の半径も変わり得るため、同じ打切りで得た四倍は全エネルギーの厳密な比でない。
:::

### 二つの臨界磁場

$\kappa=\lambda/\xi>1/\sqrt2$ の一成分 GL 模型では正常相と超伝導相の平面界面エネルギーが負になり、場の中で渦が入る混合状態が安定になり得る。界面エネルギーの符号は後で固定印加場の汎関数を用いて導く。まず上下の臨界場を求める。

バルクの上部臨界場の近くでは $\psi$ が小さいため GL 方程式を線形化する。表面とスピンの対破壊を省き、正常相の磁束密度を一様とすれば運動項の最低固有値は式 :eqref[eq-ss-landau] と同様に $\hbar|q_s|B/(2m_s)$ である。これが $|a|$ に等しいとき二次係数は零になり、場をさらに下げると正常状態は線形に不安定になる。従って軌道効果によるバルクの上部臨界場は
$$
 B_{c2}=\frac{\Phi_0}{2\pi\xi^2}.
 \label{eq-ss-upper-critical}
$$
他方、$\kappa\gg1$ の孤立渦では $\xi<r<\lambda$ の位相勾配が $1/r$ となる。磁場による相殺をこの対数領域の先頭近似で省くと運動エネルギーは単位長さ当たり $\int_\xi^\lambda (\hbar^2 n_s/2m_s)r^{-2}\,2\pi r\,\dd r$ であり、$\Phi_0^2\log\kappa/(4\pi\mu_0\lambda^2)$ となる。固定印加 $H$ で渦一本を入れる磁気仕事は単位長さ当たり $H\Phi_0$ なので
$$
 B_{c1}\simeq\frac{\Phi_0}{4\pi\lambda^2}\log\kappa
$$
である。対数の外側の寄与と芯のエネルギーは $\log\kappa$ に対して定数の補正になる。混合状態の渦密度は平均内部磁束密度 $\bar B$ に対して $\bar B/\Phi_0$ である。

:::example[上部臨界場の大きさ]
$\xi=10\,\mathrm{nm}$ なら $B_{c2}\simeq3.3\,\mathrm T$ である。$\lambda=200\,\mathrm{nm}$ なら $\kappa=20$ なので、下部臨界場の対数近似は約 $0.012\,\mathrm T$ である。
:::

:::example[渦格子の間隔]
渦を三角格子に並べると一渦の面積は $\sqrt3 a_v^2/2$ である。$\bar B=1\,\mathrm T$ なら $a_v=\sqrt{2\Phi_0/(\sqrt3\bar B)}\simeq49\,\mathrm{nm}$ になる。
:::

### 上部臨界場近くの渦格子

$\kappa\gg1$ で $B_{c2}$ のすぐ下を考える。磁場の空間変動を省き、平均磁束密度 $\bar B$ を固定する。線形運動項が最も低い Landau 準位の中で磁気並進に従う周期配置を取り、$\psi=\sqrt{Q}\,u$、胞内平均を $\langle |u|^2\rangle=1$ に規格化する。一胞に一量子の磁束を通すとその面積は $\Phi_0/\bar B$ である。これは胞境界を一周する際の磁気並進の位相が、秩序変数に対して $2\pi$ の整数倍になる条件から得られる。

形を固定したとき、正常相からの自由エネルギー密度差は

$$
 f_s-f_n=\alpha_B Q+\frac b2\beta_A Q^2,\qquad
 \alpha_B=a+\frac{\hbar|q_s|\bar B}{2m_s},\qquad
 \beta_A=\langle|u|^4\rangle
$$

である。最低準位では二次項が形によらず同じなので、形の選択には四次項を使う。$\alpha_B<0$ なら $Q=-\alpha_B/(b\beta_A)$ で最小になり、その値は $-\alpha_B^2/(2b\beta_A)$ である。従って $\beta_A$ が小さいほど低い。最低 Landau 準位の周期解の密度を逆格子で展開すると、一量子の正方格子では $\beta_A\simeq1.180$、三角格子では $\beta_A\simeq1.160$ となる。周期解の構成と一般の格子形状の最小化は :cite[tinkham-sc] を参照する。この近似では三角格子の凝縮利得の絶対値は正方格子の約 $1.017$ 倍である。異方性や磁場の空間変動を許す場合は、それらのエネルギーも形の比較に加わる。

### 平面界面での平方完成

外場を $H_c=B_c/\mu_0$ に固定し、$x\to-\infty$ で完全排除する超伝導相、$x\to+\infty$ で $B=B_c$ の正常相とする。$X=x/\lambda$、$\psi=\psi_0 f(X)$、$A_y=B_c\lambda u(X)$ と置く。ゲージと位相を選び $\psi_0,f$ は実で、バルクの超伝導側で $u=0$ とする。一枚の平面界面を表すため、$f$ は $1$ から $0$ まで単調に下がる配置を考える。$B/B_c=u'$ であり、ダッシュは $X$ の微分である。場を固定すると磁場の自由エネルギー密度は $B^2/(2\mu_0)-H_cB$ となる。両側で等しいバルク密度 $-B_c^2/(2\mu_0)$ を引き、式 :eqref[eq-ss-gl-lengths] を使って係数を無次元化する。単位面積当たりの界面エネルギーを $B_c^2\lambda/(2\mu_0)$ で割ると

$$
\int_{\R}\left[(u'-1)^2-2f^2+f^4+u^2 f^2+
 \frac{2}{\kappa^2}f'^2\right]\dd X
$$

になる。勾配項と $u^2 f^2$ は、GL の運動項の $x$ 成分と $y$ 成分からそれぞれ出る。この積分の被積分関数は

$$
(2f'+uf)^2+(u'-1+f^2)^2
+\left(\frac{2}{\kappa^2}-4\right)f'^2-2(uf^2)'
\label{eq-ss-gl-interface-squares}
$$

に等しい。展開すると交差項が $4uff'+2u'f^2=2(uf^2)'$ にまとまる。最後の積分は、超伝導側で $u\to0$、正常側で $uf^2\to0$ となる界面では零である。後者は正常側で磁場が有限かつ秩序変数が十分速く減衰する境界条件である。

$\kappa=1/\sqrt2$ では平方の和だけが残り、最小値は負にならない。零を実現するには $2f'=-uf$、$u'=1-f^2$ とすればよい。$f=1,u=0$ から積分して $u^2=2f^2-2-4\log f$ を得る。この式の右辺は $0<f<1$ で正である。従って正の平方根を選び $\dd X=-2\dd f/(uf)$ と分離すれば界面を構成できる。$f\uparrow1$ では $u\sim2(1-f)$ なので $X\to-\infty$、$f\downarrow0$ では $u\sim2\sqrt{-\log f}$ なので $X\to+\infty$ となる。正常側で $u'\to1$、$u f^2\to0$ も従う。従って界面エネルギーはこの値の $\kappa$ で厳密に零となる。

$\kappa>1/\sqrt2$ なら同じ無次元プロファイルを試行配置に使って、式 :eqref[eq-ss-gl-interface-squares] の $f'^2$ の係数だけを負にできる。従って最小の平面界面エネルギーは負になる。$\kappa<1/\sqrt2$ では逆に $f'^2$ に正の係数が残る。この場合に最小値が零へ近づけると仮定すると、$\int f'^2$ と二つの平方の積分が同時に零へ近づく。$f$ が $3/4$ から $1/4$ まで下がる区間で $1/4\le f\le3/4$ を保つ部分を取る。その長さは Cauchy–Schwarz により $1/(4\int f'^2)$ 以上なので発散する。最初の平方が小さいことから、この区間で $\int u^2$ も零へ近づく。一方、長さ $1$ の部分区間に台を持つ非負のなめらかな試験関数を使うと、二番目の平方により $u'$ の試験関数との積分は、正の下限 $7/16$ をもつ $1-f^2$ の積分に近い。部分積分で同じ積分を $u$ と試験関数の微分の積分に変えると零へ近づくので矛盾する。従ってこの側の最小界面エネルギーは正である。

:::example[上部臨界場との比較]
式 :eqref[eq-ss-upper-critical] は $B_{c2}=\hbar/(|q_s|\xi^2)$ である。一方、零場の凝縮エネルギーから得る $B_c$ は $\hbar/(\sqrt2|q_s|\lambda\xi)$ なので、比は $\sqrt2\kappa$ である。界面エネルギーが零の点では二つの場も一致する。例えば $\kappa=2$ の GL 模型では、正常相が線形に不安定になる場 $B_{c2}$ は熱力学的な $B_c$ の約 $2.83$ 倍になる。
:::

## BCS 理論の概観

### 対形成の不安定性

Fermi 面より上のエネルギー幅 $\Lambda$ に一スピン当たり一定の状態密度 $\nu$ があり、反対波数・反対スピン間に有効引力 $-g$ がある模型を考える。二粒子波動関数の係数 $c_k$ は束縛エネルギー $E_b>0$ に対し $(2\xi_k+E_b)c_k=g\sum_{k'}c_{k'}$ を満たす。非零解の条件は $1=g\nu\int_0^\Lambda\dd\xi/(2\xi+E_b)=(g\nu/2)\log(1+2\Lambda/E_b)$ である。任意の $g>0$ に対し根があり、弱結合で $E_b\simeq2\Lambda e^{-2/(g\nu)}$ になる。これは Fermi 面の存在が作る対の不安定性である。格子振動を介した有効引力の定量化には、周波数依存の相互作用と Coulomb 斥力の扱いが必要なので、ここでは有効結合を入力とする。強結合の扱いは :cite[tinkham-sc] を参照する。

### 多数の対と励起エネルギー

異常平均 $\langle c_{-\vb k\downarrow}c_{\vb k\uparrow}\rangle$ を用い、$\Delta=(g/V)\sum_{\vb k}\langle c_{-\vb k\downarrow}c_{\vb k\uparrow}\rangle$ を引力の平均場とする。粒子・正孔演算子を $\Psi_{\vb k}=(c_{\vb k\uparrow},c_{-\vb k\downarrow}^\dagger)^{\mathsf T}$ とすると、その二次形式 $\Psi_{\vb k}^\dagger H_{\vb k}\Psi_{\vb k}$ の非対角成分は対の生成・消滅項である。引力をこの規約で分解すると $H_{\vb k}=\begin{pmatrix}\xi_k&-\Delta\\-\Delta^*&-\xi_k\end{pmatrix}$ を得る。二乗は $(\xi_k^2+|\Delta|^2)I$ なので、準粒子エネルギーは $E_k=\sqrt{\xi_k^2+|\Delta|^2}$ となる。行列の二つの射影 $P_\pm=(I\pm H_{\vb k}/E_k)/2$ に占有重みを掛けると、$\langle(\Psi_{\vb k})_2^\dagger(\Psi_{\vb k})_1\rangle$ は $[f(E_k)P_++f(-E_k)P_-]_{12}=\Delta\tanh[E_k/(2k_BT)]/(2E_k)$ である。これを $\Delta$ の定義に戻し、非零の $\Delta$ で割る。対称な打切り $-\Lambda<\xi_k<\Lambda$ では、被積分関数が偶なので正側だけを積分すると
$$
 1=g\nu\int_0^\Lambda
 \frac{\tanh(\sqrt{\xi^2+|\Delta|^2}/(2k_BT))}
 {\sqrt{\xi^2+|\Delta|^2}}\dd\xi .
 \label{eq-ss-bcs-gap}
$$
状態密度 $\nu$ は単位体積、一スピン当たりとし、$g\nu$ は無次元である。$T=0$ では積分は $\operatorname{arsinh}(\Lambda/|\Delta_0|)$ で、弱結合で $\Delta_0\simeq2\Lambda e^{-1/(g\nu)}$ となる。$T_c$ では $\Delta\to0$ を代入し $\int_0^{\Lambda/(2k_BT_c)}\tanh x\,\dd x/x$ を評価する。次の導出から弱結合では $k_BT_c\simeq(2e^{\gamma_E}/\pi)\Lambda e^{-1/(g\nu)}$、従って $2\Delta_0/(k_BT_c)\simeq3.53$ となる。この比は弱結合で等方的なギャップの結果である。

### 熱で切られる対数の定数

$L=\Lambda/(2k_BT_c)$ とする。積分は $\tanh x$ が $1$ に近い領域から主に対数を得る。定数項まで求めるには部分積分して

$$
\int_0^L\frac{\tanh x}{x}\dd x
=\tanh L\log L-\int_0^L\frac{\log x}{\cosh^2x}\dd x
=\log L-M'(0)+o(1),\qquad
M(s)=\int_0^\infty\frac{x^s}{\cosh^2x}\dd x
$$

とする。原点で $\tanh x\log x\to0$ であり、上端から無限遠までの誤差は指数的に小さい。$M(s)$ は $\Re s>-1$ で収束し、原点近くの $x^s\log x$ も可積分なので導関数を積分の内側で取れる。$x>0$ では幾何級数を微分して $\cosh^{-2}x=4\sum_{j\ge1}(-1)^{j-1}j e^{-2jx}$ である。まず $\Re s>1$ で絶対収束する積分と和を交換すると

$$
M(s)=2^{1-s}\Gamma(1+s)(1-2^{1-s})\zeta(s)
$$

となる。Gamma 積分は変数 $2jx$ で評価し、交代和では偶数を二回引いた。$s=1$ の極は $1-2^{1-s}$ の零で消えるので、:ref[thm-cpx-zeta-functional] の解析接続と一致の定理により右辺は $\Re s>-1$ に延長する。

微分に必要なゼータ関数の値は計算し直せる。$f_s(x)=x^{-s}$ とし、整数区間で台形則の和から積分を引く。各区間の誤差は $\sup|f_s''|$ の定数倍であり、その和の尾部は $\Re s>-1$ で局所一様に収束する。従って

$$
\zeta(s)=\lim_{N\to\infty}
\left[\sum_{n=1}^N n^{-s}-\frac{N^{1-s}}{1-s}-\frac{N^{-s}}2\right]
$$

が成り立つ。$\Re s>1$ ではもとの級数に一致するため、これが原点付近での解析接続である。$s=0$ で値は $-1/2$、微分は $-\log N!+(N+1/2)\log N-N$ の極限である。:ref[thm-cpx-stirling] を使うと $\zeta'(0)=-(1/2)\log(2\pi)$ になる。また :ref[thm-cpx-gamma-extension] の無限積を $s$ の一次まで取ると $\Gamma'(1)=-\gamma_E$ である。従って $M(0)=1$、$M'(0)=-\log2-\gamma_E+\log(\pi/2)$ となる。

以上から

$$
\int_0^L\frac{\tanh x}{x}\dd x
=\log\!\left(\frac{4e^{\gamma_E}}{\pi}L\right)+o(1),\qquad
k_BT_c\simeq\frac{2e^{\gamma_E}}{\pi}\Lambda e^{-1/(g\nu)}
\label{eq-ss-bcs-thermal-log}
$$

である。$\Delta_0\simeq2\Lambda e^{-1/(g\nu)}$ と比を取れば $2\Delta_0/(k_BT_c)\simeq2\pi e^{-\gamma_E}$ となる。指数の部分は同じ対数の発散から出ている。係数の差は、零温のギャップと有限温度の占有による切断の仕方が異なることで生じる。

:::example[ギャップの量]
$T_c=10\,\mathrm K$ の場合、この比を使うと $\Delta_0\simeq1.52\,\mathrm{meV}$ である。二個の準粒子を同時に作る最低エネルギーは $2\Delta_0\simeq3.04\,\mathrm{meV}$ になる。
:::

:::example[トンネルの状態密度]
正常状態の一スピン状態密度を $\nu$ とし、粒子・正孔対称、損失なしを仮定する。$E=\sqrt{\xi^2+\Delta^2}$ の逆変換から $\dd\xi/\dd E=E/\sqrt{E^2-\Delta^2}$ である。二つの $\xi$ の枝の電子重みの和が $1$ になるため、正エネルギーの電子状態密度は $\nu E/\sqrt{E^2-\Delta^2}$ となり、$E<\Delta$ は零である。
:::

### 凝縮エネルギーを平均場から求める

体積 $V$ の対演算子を $b_{\vb k}=c_{-\vb k\downarrow}c_{\vb k\uparrow}$ とし、引力を $-(g/V)\sum_{\vb k,\vb k'}b_{\vb k}^\dagger b_{\vb k'}$ で表す。$\Delta=(g/V)\sum_{\vb k}\langle b_{\vb k}\rangle$ を定めて、各演算子を平均値とずれに分ける。ずれ同士の積を落とすと、引力の項は $-\sum_{\vb k}(\Delta b_{\vb k}^\dagger+\conj{\Delta} b_{\vb k})+V|\Delta|^2/g$ となる。この定数項は平均値を二回数えた分の補正である。

一組の状態 $(\vb k\uparrow,-\vb k\downarrow)$ の偶数占有部分の基底を、空状態と両方入った状態に取る。位相を選んで $\Delta>0$ とすれば、行列は $\begin{pmatrix}0&-\Delta\\-\Delta&2\xi_{\vb k}\end{pmatrix}$ であり、固有値は $\xi_{\vb k}\pm E_{\vb k}$ である。最低状態を $u_{\vb k}\ket{0}+v_{\vb k}b_{\vb k}^\dagger\ket{0}$ とすると、固有方程式から $v_{\vb k}/u_{\vb k}=(E_{\vb k}-\xi_{\vb k})/\Delta$ である。規格化を使えば $u_{\vb k}^2=(1+\xi_{\vb k}/E_{\vb k})/2$、$v_{\vb k}^2=(1-\xi_{\vb k}/E_{\vb k})/2$ となる。奇数占有の二状態はいずれもエネルギー $\xi_{\vb k}$ を持つので、最低状態から $E_{\vb k}$ の励起二つと $2E_{\vb k}$ の励起一つという構造である。

零温での正常状態の同じ組のエネルギーは $\xi_{\vb k}-|\xi_{\vb k}|$ である。従って化学ポテンシャル一定での大ポテンシャル密度差は

$$
\frac{\Omega_s-\Omega_n}{V}
=\frac{\Delta_0^2}{g}+2\nu\int_0^\Lambda
 \bigl[\xi-\sqrt{\xi^2+\Delta_0^2}\bigr]\dd\xi
=-\nu\Lambda\bigl(\sqrt{\Lambda^2+\Delta_0^2}-\Lambda\bigr)
\simeq-\frac{\nu\Delta_0^2}{2}
\label{eq-ss-bcs-condensation}
$$

である。積分したときの $\Delta_0^2\operatorname{arsinh}(\Lambda/\Delta_0)$ は式 :eqref[eq-ss-bcs-gap] の零温の関係によって $\Delta_0^2/g$ と相殺した。最後は $\Delta_0/\Lambda\ll1$ の展開である。状態密度が対称なこの模型では、$\xi$ と $-\xi$ の二組で $v^2(\xi)+v^2(-\xi)=1$ なので粒子数が正常状態と同じである。従って固定密度での自由エネルギー差にも同じ式が使える。式 :eqref[eq-ss-condensation-energy] と比べると $B_c(0)\simeq\Delta_0\sqrt{\mu_0\nu}$ となる。両スピンの密度を使う場合には $\nu$ をその半分とする。

:::example[比熱と凝縮エネルギーの照合]
正常相の電子比熱を単位体積当たり $C_n=\gamma_n T$ とすると、両スピンの密度が $2\nu$ なので自由電子型の係数は $\gamma_n=2\pi^2\nu k_B^2/3$ である。従って零温の凝縮エネルギーを $\gamma_nT_c^2$ で割った値は $3[\Delta_0/(k_BT_c)]^2/(4\pi^2)\simeq0.236$ である。比熱と磁化のそれぞれから独立に求められる量を、この近似により比較できる。
:::

### 転移直下のギャップと比熱跳び

$t=k_BT$、$\omega_j=(2j+1)\pi t$ とする。:ref[thm-cpx-sine-product] の余接の部分分数展開で $z=1/2+i E/(2\pi t)$ と置き、虚部を取る。$a=E/(2\pi t)$ とすると、整数についての和の虚部は $-a\sum_{n\in\Z}[(1/2-n)^2+a^2]^{-1}$、余接側は $-\pi\tanh(\pi a)$ である。これを $j\ge0$ に折り返すと

$$
\frac{\tanh(E/(2t))}{E}=4t\sum_{j=0}^\infty\frac{1}{\omega_j^2+E^2}
$$

である。式 :eqref[eq-ss-bcs-gap] の積分を $I(T,\Delta)$ とし、$I(T,\Delta)-I(T,0)$ にこの展開を使う。差では大きな $\xi$ の寄与が十分減るので、弱結合で上限を無限大に延ばせる。$\Delta^2$ の一次係数には $\int_0^\infty(\omega_j^2+\xi^2)^{-2}\dd\xi=\pi/(4\omega_j^3)$ を用いる。これは $\xi=\omega_j\tan u$ で $\omega_j^{-3}\int_0^{\pi/2}\cos^2u\,\dd u$ となる。$\sum_{j\ge0}(2j+1)^{-3}=(1-2^{-3})\zeta(3)$ なので、$\Delta/t\ll1$ で

$$
I(T,\Delta)-I(T,0)=-C_\Delta\frac{\Delta^2}{t^2}
+\order{(\Delta/t)^4},\qquad C_\Delta=\frac{7\zeta(3)}{8\pi^2}.
$$

級数での次の係数も $\sum_j\omega_j^{-5}$ が収束するため、展開の余りを積分後に抑えられる。有限の $\Lambda$ を無限大に延ばすと、$\xi\gg t$ における積分の差が $\order{\Delta^2/\xi^3}$ なので、$\order{\Delta^2/\Lambda^2}$ の補正が入る。式 :eqref[eq-ss-bcs-thermal-log] から $I(T,0)-I(T_c,0)\simeq\log(T_c/T)$ であり、ギャップ方程式より $\Delta(T)^2\simeq k_B^2 T_c(T_c-T)/C_\Delta$ となる。

有限温度で各組の励起分配関数は $(1+e^{-E_{\vb k}/t})^2$ である。従って大ポテンシャルは零温の式に $-2t\sum_{\vb k}\log(1+e^{-E_{\vb k}/t})$ を足した形である。$p=|\Delta|^2$ で微分し、$\partial_p E_{\vb k}=1/(2E_{\vb k})$ を使うと単位体積当たりの微分が $1/g-\nu I(T,\Delta)$ となる。$p=0$ を基準に上の展開を積分すると

$$
\frac{\Omega(p)-\Omega(0)}{V}
\simeq \nu\left[\log(T/T_c)p+\frac{C_\Delta}{2t^2}p^2\right].
$$

最低点の値は $-\nu k_B^2(T_c-T)^2/(2C_\Delta)$ の先頭項を持つ。温度の二階微分から比熱差を求め、正常相の $C_n(T_c)=2\pi^2\nu k_B^2T_c/3$ と比べると

$$
\frac{C_s(T_c^-)-C_n(T_c^+)}{C_n(T_c^+)}
\simeq \frac{12}{7\zeta(3)}\simeq1.43
\label{eq-ss-bcs-heat-jump}
$$

となる。この結果は平均場と弱結合の近似に依存する。エントロピーが連続である一方で比熱に跳びが出る点は、臨界磁場を用いた導出と一致する。

### 対の広がりとスピンによる臨界磁場

弱結合 BCS の基底状態で、反対運動量の二つの電子の異常平均は $u_{\vb k}v_{\vb k}=\Delta_0/(2E_{\vb k})$ である。この積が Fermi 面からのエネルギー差 $|\xi_k|\sim\Delta_0$ の範囲で変化することは、$\xi_k\simeq\hbar v_F(k-k_F)$ から波数幅 $\Delta_0/(\hbar v_F)$ に対応する。従って対相関の距離の尺度は Fourier 変換により $\hbar v_F/\Delta_0$ である。慣用的な清浄 BCS コヒーレンス長を $\xi_0=\hbar v_F/(\pi\Delta_0)$ と定める。これは弱結合の零温ギャップから定める長さであり、転移の近くで秩序変数の勾配エネルギーを表す式 :eqref[eq-ss-gl-lengths] の $\xi(T)$ とは温度依存も定義も異なる。

磁場は軌道運動のほか、電子スピンにも作用する。軌道効果を省き、温度零で一重項の一様な $s$ 波超伝導を考える。電子の $g$ 因子を $2$ とし、正常相のスピン磁化率を $\chi_P=2\mu_0\mu_B^2\nu$ とする。Fermi 面で一スピンずつのエネルギーが $\mp\mu_B B$ 動くため、正常相は $-\chi_P B^2/(2\mu_0)=-\nu\mu_B^2B^2$ の自由エネルギーを得る。完全にギャップの開いた一重項相では $\mu_B B<\Delta_0$ のとき正の準粒子エネルギーが残り、零温でスピン磁化は生じない。従って式 :eqref[eq-ss-bcs-condensation] の凝縮利得と正常相のスピン利得は

$
 B_P=\frac{\Delta_0}{\sqrt2\,\mu_B}
$

で等しくなる。これが一様な弱結合模型の :term[Pauli 限界]{en="Pauli limit" yomi="ぱうりげんかい"} である。ここでは磁場を $\mu_0H$ とみなし、スピン軌道結合や空間的に変調した対形成を省いた。実際の臨界場を求めるには、軌道効果とスピン効果が同時に働く場合の安定性を調べる。

:::example[二つの磁場尺度]
$T_c=10\,\mathrm K$ では弱結合の $\Delta_0\simeq1.52\,\mathrm{meV}$ から $B_P\simeq18.6\,\mathrm T$ を得る。一方、GL の式で $\xi=10\,\mathrm{nm}$ のときの軌道磁場は約 $3.3\,\mathrm T$ である。GL の値を低温まで外挿する前に、その $\xi$ をどの温度域で求めたかを確認する。
:::

## Josephson 効果と SQUID

### 位相差と電流

絶縁層で弱く結合した二つの凝縮体を考える。最低次の位相依存のエネルギーは $-E_J\cos\delta$ であり、$\delta$ は $\theta_2-\theta_1-(q_s/\hbar)\int_1^2\vb A\cdot\dd\vb l$ である。これはゲージ変換で両端の位相差と線積分の変化が相殺する。線積分に対するエネルギーの負の微分は $q_sE_J\sin\delta/\hbar$ となる。電流と電圧の向き、および位相差の符号をそろえ $\varphi$ を選ぶと
$$
 I=I_c\sin\varphi,\qquad I_c=2eE_J/\hbar,\qquad
 \dot\varphi=2eV/\hbar .
 \label{eq-ss-josephson}
$$
最後の式は、凝縮粒子一個の電気化学エネルギー差 $2eV$ が、時間発展位相の相対速度に $\hbar$ を掛けた値であることによる。向きを反転すれば $\varphi,I,V$ が同時に符号を変える。電圧零でも $|\sin\varphi|\le1$ の範囲で電流が流れ、一定電圧なら振動周波数は $2e|V|/h$ である。マイクロ波に同期して $\dot\varphi$ の平均が $2\pi f$ の整数倍となると、平均電圧は $V_n=nhf/(2e)$ となる。

### 二つの接合の干渉

自己インダクタンスを無視した環に二つの同じ接合を置く。磁束を $\Phi$ とすると閉路の位相一価性から $\varphi_1-\varphi_2=2\pi\Phi/\Phi_0$ を $2\pi$ を法として得る。二つの枝の電流を加えると
$$
 I=2I_c\sin[(\varphi_1+\varphi_2)/2]\cos(\pi\Phi/\Phi_0).
$$
最大絶対電流は $2I_c|\cos(\pi\Phi/\Phi_0)|$ であり、磁束の周期は $\Phi_0$ である。これが直流 SQUID の原理である。接合の臨界電流が異なると打ち消し点でも有限な電流が残る。

:::example[電圧を周波数に換算する]
$V=1\,\mathrm{\mu V}$ で周波数は $V/\Phi_0\simeq4.84\times10^8\,\mathrm{Hz}$ である。磁束量子は $\mathrm{V\,s}$ と同じ単位であり、この換算が成り立つ。
:::

:::example[半磁束量子での打ち消し]
$\Phi=\Phi_0/2$ では二接合の位相差が $\pi$ であり、同じ臨界電流なら和は零である。もし臨界電流が $I_a=2I_b$ なら和は $(I_a-I_b)\sin\varphi_1$ となり、最大値は $I_b$ である。
:::

### 抵抗と容量を含めた接合

理想な接合の電流に、同じ電圧で流れる正常成分 $V/R$ と、電荷の蓄積による変位電流 $C\dot V$ を加える。ここで $R>0$ は抵抗、$C$ は容量であり、素子を集中定数で扱う。この並列等価回路が :term[RCSJ 模型]{en="resistively and capacitively shunted junction model" yomi="あーるしーえすじぇーもけい"} である。$A_J=\hbar/(2e)=\Phi_0/(2\pi)$ と置く。式 :eqref[eq-ss-josephson] より $V=A_J\dot\varphi$ なので、外部から加える電流 $I_{\rm b}$ に対して

$$
CA_J\ddot\varphi+\frac{A_J}{R}\dot\varphi+I_c\sin\varphi=I_{\rm b}
\label{eq-ss-junction-rcsj}
$$

となる。$I_{\rm b}$ を一定にしたとき、式に $A_J$ を掛けると慣性係数 $CA_J^2$ の粒子が $U(\varphi)=-E_J\cos\varphi-A_J I_{\rm b}\varphi$ で運動する式になる。さらに $\dot\varphi$ を掛けると

$$
\frac{\dd}{\dd t}\left[\tfrac12 CA_J^2\dot\varphi^2+U(\varphi)\right]
=-\frac{A_J^2}{R}\dot\varphi^2=-V^2/R
$$

である。慣性項のエネルギーは容量に蓄えられる $CV^2/2$ であり、右辺は抵抗で失われる電力である。傾きのある $U$ には電流源がする仕事も含まれる。

$|I_{\rm b}|<I_c$ では $\sin\varphi_0=I_{\rm b}/I_c$ の解のうち $\cos\varphi_0>0$ が静止の安定点である。$C>0$ として小さい変位を考えれば

$$
\ddot{\delta\varphi}+\frac{1}{RC}\dot{\delta\varphi}+
\omega_{\rm Jp}^2\delta\varphi=0,\qquad
\omega_{\rm Jp}^2=\frac{I_c}{CA_J}\sqrt{1-(I_{\rm b}/I_c)^2}
$$

となる。減衰を無視した $\omega_{\rm Jp}$ は接合の小振幅のプラズマ角振動数であり、一定電圧で回り続ける位相の交流 Josephson 周波数とは異なる。減衰が弱ければこの振動を観測でき、減衰が強ければ静止点へ単調に戻る場合がある。容量に十分なエネルギーが蓄えられると位相が隣の極小を通り越すことがあり、その場合は電流掃引の上りと下りで電圧の現れ方が違う。

:::example[接合の共鳴の桁]
$I_c=2.0\,\mathrm{\mu A}$、$C=100\,\mathrm{fF}$ とする。無バイアスでの $\omega_{\rm Jp}$ は $\sqrt{I_c/(CA_J)}\simeq2.47\times10^{11}\,\mathrm{s^{-1}}$、普通振動数では約 $39\,\mathrm{GHz}$ である。$I_{\rm b}=0.60 I_c$ にすると、無減衰の振動数は $(1-0.60^2)^{1/4}\simeq0.894$ 倍に下がる。$I_{\rm b}$ が $I_c$ に近づくと曲率が小さくなるため、同じ $R,C$ でも減衰を無視しにくくなる。
:::

### 強減衰の電流–電圧特性

容量項が無視でき、一定バイアスが $I_{\rm b}>I_c$ の場合には $\dot\varphi=(R/A_J)(I_{\rm b}-I_c\sin\varphi)>0$ である。位相が一周する時間 $T_\varphi$ は

$$
T_\varphi=\frac{A_J}{R}\int_0^{2\pi}
\frac{\dd\varphi}{I_{\rm b}-I_c\sin\varphi}
=\frac{2\pi A_J}{R\sqrt{I_{\rm b}^2-I_c^2}} .
$$

積分の評価では角度をずらして $\cos\varphi$ を使い、$-\pi<\varphi<\pi$ で $u=\tan(\varphi/2)$ とする。すると積分は $\int_{\R}2\dd u/[(I_{\rm b}-I_c)+(I_{\rm b}+I_c)u^2]$ であり、変数を正の定数で伸縮すれば $\int_{\R}\dd u/(1+u^2)=\pi$ になる。平均電圧は一周の位相差 $2\pi$ に $A_J/T_\varphi$ を掛けるため

$$
\overline V=R\sqrt{I_{\rm b}^2-I_c^2}\qquad (I_{\rm b}>I_c)
\label{eq-ss-junction-rsj-voltage}
$$

である。負電流の場合は $\varphi\mapsto-\varphi$ により電圧の符号が反転する。$|I_{\rm b}|\le I_c$ の定常な静止解では $\overline V=0$ である。これは抵抗が一定で熱揺らぎがない回路模型の結果である。

:::example[しきい値のすぐ上の電圧]
$I_{\rm b}=1.25I_c$ なら $\overline V/(RI_c)=\sqrt{1.25^2-1}=0.75$ である。$I_c=20\,\mathrm{\mu A}$、$R=10\,\mathrm\Omega$ なら平均電圧は $150\,\mathrm{\mu V}$ となる。正常成分だけの電圧 $RI_{\rm b}=250\,\mathrm{\mu V}$ との差は、周回中も正弦の超電流が流れ、位相が各角度に同じ時間を費やさないことから生じる。
:::

### 交流に同期した位相

一定の電圧と単一の交流電圧 $V(t)=V_0+V_1\cos\omega t$ を接合に与えた場合を考える。位相は $\varphi(t)=\varphi_0+(V_0/A_J)t+x\sin\omega t$、$x=V_1/(A_J\omega)$ である。$e^{ix\sin\omega t}$ を Fourier 展開し、その係数を $J_j(x)$ と書く。この $J_j$ は第一種 Bessel 関数である。係数は $e^{(x/2)e^{i\omega t}}e^{-(x/2)e^{-i\omega t}}$ の二つの指数級数を掛けて、同じ $e^{ij\omega t}$ の項を集めれば得られる。位相の指数の平均で零周波数が出るには $V_0/A_J+j\omega=0$ が必要である。

従って $V_0=n A_J\omega= n\Phi_0 f$ なら、超電流に $I_cJ_{-n}(x)\sin\varphi_0$ の直流成分が残る。$f=\omega/(2\pi)$ である。電流バイアスに交流を重ねる測定では、位相が同期できる範囲の電流で電圧がこの値に留まり、階段状の特性を示すことがある。段の位置は周波数の整数比で決まるが、段の幅を計算するには回路の動力学と交流の駆動条件も必要である。

:::example[周波数と段間隔]
$f=10.0\,\mathrm{GHz}$ での電圧間隔は $\Phi_0 f\simeq20.7\,\mathrm{\mu V}$ であり、第 $3$ 段では約 $62.0\,\mathrm{\mu V}$ である。段の間隔は $I_c$ を式に含まないので、電圧と周波数の換算として使える。
:::

### 接合の幅と磁束の干渉

二つの接合の臨界電流が $I_a,I_b$ で異なる場合でも、自己インダクタンスを無視できれば位相差の差は $\delta=2\pi\Phi/\Phi_0$ である。全電流は複素数 $(I_a e^{i\delta/2}+I_b e^{-i\delta/2})e^{i(\varphi_1+\varphi_2)/2}$ の虚部であるから、最大は $\sqrt{I_a^2+I_b^2+2I_aI_b\cos\delta}$ となる。この値は半磁束量子で $|I_a-I_b|$ へ下がる。大きい自己インダクタンスでは全磁束が外部磁束に環の循環電流の寄与を足した値になるので、外部磁束をそのまま $\delta$ に代入できない。

幅 $W$ の単一接合でも場所ごとに位相が異なれば干渉が起こる。電流が作る磁場を無視し、一様な接合の幅方向 $y$ にわたって有効な磁束が $\Phi$ だとする。接合を横切る経路を幅方向にずらしたとき、二経路を結ぶ閉路の磁束差から $\varphi(y)=\varphi_0+2\pi\Phi y/(W\Phi_0)$ となる。有効な磁束には電極表面の侵入層までの磁束も含まれ得る。一様な局所臨界電流密度で電流を積分し、$\varphi_0$ を調整して絶対値を最大化すると

$$
I_{\rm max}(\Phi)=I_{\rm max}(0)
\left|\frac{\sin(\pi\Phi/\Phi_0)}{\pi\Phi/\Phi_0}\right|
$$

となる。実際、$y\in[-W/2,W/2]$ で $\sin\varphi(y)$ の積分のうち奇の部分が消え、$\sin\varphi_0$ に有限幅の余弦の積分が掛かる。一接合の幅での連続的な位相の差なので、二接合の周期的な干渉とは包絡の形も零点も異なる。例えば同じ磁束 $\Phi_0/2$ で、幅の一様な単一接合は $2/\pi$ 倍、等しい二接合の環は零になる。

## 演習問題

::::exercise{level=basic}
半無限の超伝導体で磁場が表面の $1/10$ になる深さを求めよ。
:::solution
$e^{-x/\lambda}=1/10$ の対数を取ると $x=\lambda\log10$ である。
:::
::::

::::exercise{level=basic}
$B_c$ を半分にしたときの凝縮エネルギー密度の比を求めよ。
:::solution
式 :eqref[eq-ss-condensation-energy] は $B_c^2$ に比例するので比は $1/4$ である。
:::
::::

::::exercise{level=basic}
$\lambda=150\,\mathrm{nm}$、$\xi=30\,\mathrm{nm}$ で GL の $\kappa$ を求めよ。
:::solution
比は $150/30=5$ であり、第二種の境界 $1/\sqrt2$ より大きい。
:::
::::

::::exercise{level=basic}
一定直流電圧 $V=4\,\mathrm{\mu V}$ と $V=1\,\mathrm{\mu V}$ の Josephson 周波数比を求めよ。
:::solution
式 :eqref[eq-ss-josephson] の周波数は $|V|/\Phi_0$ なので、比は $4$ である。
:::
::::

::::exercise{level=advanced}
磁束 $\Phi=\Phi_0/3$ の等しい二接合の SQUID で、零磁束における最大電流に対する比を求めよ。
:::solution
最大電流は $2I_c|\cos(\pi\Phi/\Phi_0)|$ である。零場では $2I_c$、$\Phi_0/3$ では $2I_c\cos(\pi/3)=I_c$ なので比は $1/2$ である。
:::
::::

::::exercise{level=advanced}
GL の一次元界面解 $f(x)=\tanh(x/(\sqrt2\xi))$ を方程式に代入して確認せよ。
:::solution
$y=x/(\sqrt2\xi)$ とすれば $\partial_y f=1-f^2$、$\partial_y^2 f=-2f(1-f^2)$ である。従って $-\xi^2 f''=f(1-f^2)$ が $-f+f^3$ と相殺する。$x=0$ で零、遠方で $1$ になる。
:::
::::

::::exercise{level=advanced}
BCS の零温積分を用い、弱結合近似の前に $\Delta_0$ を求めよ。
:::solution
$\xi=|\Delta_0|\sinh u$ で分母と Jacobian が相殺する。式 :eqref[eq-ss-bcs-gap] は $1/(g\nu)=\operatorname{arsinh}(\Lambda/|\Delta_0|)$ となるので、$|\Delta_0|=\Lambda/\sinh[1/(g\nu)]$ である。$g\nu\ll1$ で分母の先頭は $e^{1/(g\nu)}/2$ となる。
:::
::::

::::exercise{level=advanced}
有限厚の板の平均磁束密度を求めよ。厚さは $2d$、面に平行な表面の磁束密度は $B_a$ とする。
:::solution
$B(x)/B_a=\cosh(x/\lambda)/\cosh(d/\lambda)$ を $-d$ から $d$ まで積分し、$2d$ で割る。従って $\bar B/B_a=(\lambda/d)\tanh(d/\lambda)$ である。薄い極限で $1$、厚い極限で $\lambda/d$ に近づく。
:::
::::

::::exercise{level=basic}
抵抗と容量を含む Josephson 接合で、無バイアスの $I_c$ を一定に保ち、容量を四倍にすると小振幅のプラズマ角振動数は何倍になるか。減衰は無視する。
:::solution
式 :eqref[eq-ss-junction-rcsj] を $I_{\rm b}=0,\varphi=0$ のまわりで線形化すると、$\omega_{\rm Jp}^2=I_c/(CA_J)$ となる。従って $\omega_{\rm Jp}\propto C^{-1/2}$ であり、容量を四倍にすると二分の一になる。
:::
::::

::::exercise{level=advanced}
一量子渦で $r>\xi$ の方位角電流が $j_\phi=-B_z'/\mu_0$ のとき、式 :eqref[eq-ss-london-vortex-field] による半径 $R>\xi$ までの磁束と、式 :eqref[eq-ss-fluxoid] の電流循環項の和を求めよ。芯も含めて点線源模型の磁束を積分する。
:::solution
$u=R/\lambda$ とすると磁束は $\Phi_0\int_0^u vK_0(v)\dd v=\Phi_0[1-uK_1(u)]$ である。$(uK_1)'=-uK_0$ と原点の $uK_1\to1$ を用いた。電流側は $\mu_0\lambda^2\cdot2\pi R\cdot\Phi_0 K_1(u)/(2\pi\mu_0\lambda^3)=\Phi_0uK_1(u)$ である。和は $\Phi_0$ で、磁束単独がこの値へ近づくには $R\gg\lambda$ が必要である。
:::
::::

::::exercise{level=advanced}
BCS の零温凝縮エネルギーの厳密な打切り模型の式 :eqref[eq-ss-bcs-condensation] を用い、先頭の $-\nu\Delta_0^2/2$ に対する次の補正を求めよ。
:::solution
$z=\Delta_0^2/\Lambda^2$ と置くと式は $-\nu\Lambda^2(\sqrt{1+z}-1)$ である。平方根を $1+z/2-z^2/8+\order{z^3}$ と展開すれば、先頭に続く補正は $+\nu\Delta_0^4/(8\Lambda^2)$ である。負の凝縮自由エネルギー差の絶対値を小さくする補正である。
:::
::::

::::exercise{level=advanced}
容量を無視した抵抗付き接合で、平均電圧が $RI_c$ となる正のバイアス電流を求めよ。
:::solution
式 :eqref[eq-ss-junction-rsj-voltage] に $\overline V=RI_c$ を入れ、$R$ で割って二乗すると $I_c^2=I_{\rm b}^2-I_c^2$ になる。正の枝を取れば $I_{\rm b}=\sqrt2 I_c$ であり、式の前提 $I_{\rm b}>I_c$ を満たす。
:::
::::

## 確認テスト

::::quiz
零抵抗を示す古典導体の式 $\vb E=0$ から、場を先にかけた場合に導けるものはどれか。

- [x] Faraday の法則により初期磁場が保存する。
- [ ] 必ず熱平衡で内部磁場が零に排除される。
- [ ] 磁束量子が必ず $h/(2e)$ になる。

:::explanation
$\partial_t\vb B=-\rot\vb E$ によれば、$\vb E=0$ では時間変化が零である。超伝導の平衡 Meissner 効果を導くには、凝縮位相に由来する電流と磁場の関係が必要になる。
:::
::::

::::quiz{answer="4" tol="1e-6"}
$\xi$ を半分にしたときの $B_{c2}$ の倍率を答えよ。
:::explanation
式 :eqref[eq-ss-upper-critical] では磁束量子を固定すると $B_{c2}\propto1/\xi^2$ である。従って倍率は $4$ である。
:::
::::

::::quiz
磁束と超電流の循環を足した量が量子化される理由はどれか。

- [x] 秩序変数の位相が閉曲線を回って $2\pi$ の整数倍だけ変わる。
- [ ] 任意の位置で電流が零だから。
- [ ] すべての閉曲線の磁束が同じだから。

:::explanation
位相勾配の積分が $2\pi n$ である。ゲージ共変な運動量式を積分すると、磁束と電流循環がその整数に対応する。電流が無視できる経路でのみ磁束単独の量子化に還元される。
:::
::::

::::quiz{answer="0" tol="1e-6"}
同じ二接合の SQUID に $\Phi_0/2$ を通したとき、最大超電流を一接合の $I_c$ 単位で答えよ。自己インダクタンスは無視する。
:::explanation
二枝の位相差が $\pi$ であるため、一枝の正弦ともう一枝の正弦が打ち消し合う。従って最大絶対電流も $0$ である。
:::
::::

::::quiz
一成分 GL 模型で $\kappa=1/\sqrt2$ のとき、外場を $B_c/\mu_0$ に固定した平面界面について正しいものはどれか。

- [x] 境界項が消え、二つの平方を同時に零にする配置の界面エネルギーは零である。
- [ ] 磁場が常伝導側でも完全に排除される。
- [ ] 二つの平方の前の係数はどちらも負になる。

:::explanation
式 :eqref[eq-ss-gl-interface-squares] では $2/\kappa^2-4=0$ であり、両端の $uf^2$ が零なので積分は二つの平方だけになる。$2f'=-uf$ と $u'=1-f^2$ を満たす界面ではその平方がともに零になる。正常側の $B/B_c=u'$ は $1$ に近づく。
:::
::::

::::quiz{answer="1.43" tol="0.01"}
弱結合の等方 BCS 模型について、転移での電子比熱の跳びを正常相の電子比熱で割った値を小数で答えよ。
:::explanation
式 :eqref[eq-ss-bcs-heat-jump] の比は $12/[7\zeta(3)]$ である。$\zeta(3)\simeq1.20206$ を代入すると $1.426\ldots$ であり、小数第二位までで $1.43$ になる。跳びには転移直下の超伝導相から直上の正常相を引いた値を使う。
:::
::::

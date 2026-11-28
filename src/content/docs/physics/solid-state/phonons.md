---
title: 格子振動とフォノン
description: 結晶の基準振動を導出して量子化し、散乱実験、比熱、熱膨張、熱伝導を説明する。
notation:
  - symbol: '\Phi_{s\alpha,t\beta}(\vb R)'
    meaning: 格子の力定数行列
  - symbol: 'D(\vb q)'
    meaning: 質量規格化した力学行列
  - symbol: '\omega_{\vb q j}'
    meaning: 波数 q と分枝 j に属する角振動数
  - symbol: '\gamma_{\vb q j}'
    meaning: モードごとの Grüneisen パラメータ
---

## この章の目標

- 力定数から力学行列を作り、単原子鎖と二原子鎖の分散を求める。
- 格子振動を量子化し、一フォノン散乱の運動学と強度を導く。
- フォノンの状態密度から Einstein 比熱と Debye 比熱を計算する。
- 非調和性による膨張と熱抵抗を区別し、それぞれの近似を説明する。

実格子と逆格子には :ref[def-ss-bravais]、:ref[def-ss-reciprocal] の規約を使う。波数の単位は $\mathrm{m^{-1}}$、角振動数の単位は $\mathrm{s^{-1}}$ とする。普通振動数 $f$ とは $\omega=2\pi f$ の関係にある。

## 調和近似と格子の運動方程式

### 力定数の意味

平衡位置 $\vb R+\vb\tau_s$ から第 $s$ 原子が $\vb u_{\vb R s}$ だけずれたとする。電子の状態はイオン位置に追随し、イオンに対して滑らかなポテンシャルを与えるとする。平衡では一階微分が零になるので、:term[調和近似]{en="harmonic approximation" yomi="ちょうわきんじ"} では変位の二次まで残す。ここでは $U$ をその近似ポテンシャル、$\Phi$ を平衡点での二階微分の行列とする。従って
$$
 U-U_0=\frac12\sum_{\vb R,\vb R',s,t,\alpha,\beta}
 u_{\vb R s\alpha}\Phi_{s\alpha,t\beta}(\vb R'-\vb R)u_{\vb R't\beta}
$$
である。二階微分の成分を :term[力定数]{en="force constant" yomi="ちからていすう"} といい、単位は $\mathrm{J\,m^{-2}}=\mathrm{N\,m^{-1}}$ である。結晶の周期性により胞の番号をすべて同じだけずらしても $U$ は同じなので、力定数は胞の差だけに依存する。混合微分の交換から $\Phi_{s\alpha,t\beta}(\vb R)=\Phi_{t\beta,s\alpha}(-\vb R)$ である。外部の固定点がなく、エネルギーが剛体並進に対して不変なら、変位に一定のベクトルを加えても $U$ は同じである。任意の変位について二次形式の交差項が消える条件は $\sum_{\vb R',t}\Phi_{s\alpha,t\beta}(\vb R')=0$ である。これを :term[音響和則]{en="acoustic sum rule" yomi="おんきょうわそく"} という。

運動エネルギーは $\sum_{\vb R,s}M_s|\dot{\vb u}_{\vb R s}|^2/2$ である。微分して Newton の式に代入すると $M_s\ddot u_{\vb R s\alpha}=-\sum_{\vb R',t,\beta}\Phi_{s\alpha,t\beta}(\vb R'-\vb R)u_{\vb R't\beta}$ となる。平衡が局所最小なら力定数の二次形式は非負である。負の固有値が出る配置では時間依存が単振動でなく指数増大になるため、安定な平衡からの展開になっていない。

:::example[結合ばねから力定数を読む]
一次元の結合 $K(u_{n+1}-u_n)^2/2$ を全結合について足すと、$\Phi(0)=2K$、$\Phi(\pm a)=-K$ である。和が $2K-K-K=0$ なのは一様な変位でばねが伸びないためである。各サイトに固定壁とのばね $K_0$ を加えれば $\Phi(0)$ に $K_0$ が加わり、この和則は破れる。
:::

:::example[調和領域の見積もり]
単一自由度で $U=Ku^2/2+Au^3/3!$ とする。調和力は $-Ku$、三次項由来の力は $-Au^2/2$ である。従って力の相対補正は $|Au|/(2K)$ であり、例えば $A/K=1\,\mathrm{\mathring A^{-1}}$、$|u|=0.10\,\mathrm{\mathring A}$ では $0.05$ になる。
:::

### 二体ポテンシャルから力定数を作る

力定数の非対角成分の符号を理解するため、平衡時に位置の差が $\vb d$ である一組を考える。長さを $d=|\vb d|$、方向を $\hat{\vb d}=\vb d/d$、変位の差を $\vb v$ とする。$\sqrt{d^2+2\vb d\cdot\vb v+|\vb v|^2}$ の展開から

$$
 |\vb d+\vb v|-d
 =\hat{\vb d}\cdot\vb v+
 \frac{|\vb v|^2-(\hat{\vb d}\cdot\vb v)^2}{2d}
 +\order{|\vb v|^3}
$$

となる。距離にだけ依存するポテンシャル $V$ では、この伸びを $V$ の Taylor 展開へ代入して二次の項を集めればよい。結合の Hessian は

$$
 H_{\alpha\beta}(\vb d)=V''(d)\hat d_\alpha\hat d_\beta+
 \frac{V'(d)}{d}(\delta_{\alpha\beta}-\hat d_\alpha\hat d_\beta)
 \label{eq-ss-pair-hessian}
$$

である。結合のエネルギーは $\frac12(\vb u_2-\vb u_1)^{\mathsf T}H(\vb u_2-\vb u_1)$ だから、力定数の対角ブロックへ $H$ を加え、二原子間のブロックへ $-H$ を加える。個々の対ごとに音響和則が成立する。多体相互作用でも剛体並進に不変なら和則は成立するので、対ごとに組み立てられるという仮定は必要ない。

:::example[伸びと横ずれの力定数]
長さ $d=0.30\,\mathrm{nm}$ の結合が $x$ 軸方向を向き、$V''(d)=30\,\mathrm{N\,m^{-1}}$、$V'(d)=3.0\times10^{-9}\,\mathrm N$ とする。式 :eqref[eq-ss-pair-hessian] は $H=\diag(30,10,10)\,\mathrm{N\,m^{-1}}$ を与える。初期張力が零なら横ずれの固有値も零になる。実際に張力のない結合を横へ変位させると伸びが変位の二次から始まり、エネルギーは四次から始まる。
:::

:::example[剛体並進を含む自由系]
質量がともに $M$ の二原子を一本のばね $K$ で結ぶ。座標 $Q_0=\sqrt{M/2}(u_1+u_2)$、$Q_1=\sqrt{M/2}(u_2-u_1)$ を使うと運動エネルギーは $(\dot Q_0^2+\dot Q_1^2)/2$、弾性エネルギーは $KQ_1^2/M$ である。よって並進は $Q_0(t)=Q_0(0)+\dot Q_0(0)t$、内部振動は $\omega=\sqrt{2K/M}$ である。ゼロ固有値を外す際には、結晶全体の質量と運動量を別に扱っている。
:::

### 計算した力から行列を求める

解析的なエネルギーの式が得られない場合も、力 $F_i=-\partial_i U$ が分かれば、座標 $j$ を平衡から正負に $\delta$ だけずらして

$$
 \Phi_{ij}=-\frac{F_i(+\delta\vb e_j)-F_i(-\delta\vb e_j)}{2\delta}+\order{\delta^2}
$$

と推定できる。$F_i(\pm\delta\vb e_j)$ を三次まで展開すると、偶数次の項が差で消え、三次の項を $2\delta$ で割った残差が二次の誤差となる。実際の計算には、境界をまたぐ対の変位まで含めて与えた並進周期が入る。力定数がその胞より遠くへ及ぶ場合は、異なる格子ベクトルの成分が同じ周期像に重なる。したがって変位幅に加え、計算胞の大きさも変えて分散の収束を確かめる。

例えば $F(\delta)=-K\delta-A\delta^2/2-B\delta^3/6$ なら推定する曲率は $K+B\delta^2/6$ である。$\delta$ を半分にするとこの残差は四分の一になる。一方、各力に大きさ $\varepsilon_F$ の数値誤差があると差分の誤差は最大 $\varepsilon_F/\delta$ なので、変位を小さくし続けると逆に悪化し得る。

## 1 次元単原子鎖と二原子鎖

### 単原子鎖の分散

質量 $M$ の原子を間隔 $a$ で置き、ばね定数 $K$ の最近接ばねを仮定する。$M\ddot u_n=K(u_{n+1}+u_{n-1}-2u_n)$ に $u_n=u e^{i(qna-\omega t)}$ を入れると
$$
 \omega(q)=2\sqrt{K/M}|\sin(qa/2)|,\qquad
 v_g=\dv{\omega}{q}.
 \label{eq-ss-mono-phonon}
$$
周期境界で $q=2\pi j/(Na)$、逆格子の差だけなら同じ変位なので第一ゾーンには $N$ 個の独立な波数がある。$|qa|\ll1$ では $\omega\simeq a\sqrt{K/M}|q|$ である。群速度の意味は、近い波数の重ね合わせで位相を $\omega(q+\delta q)\simeq\omega(q)+v_g\delta q$ と展開すると、包絡が $x-v_gt$ に依存することから分かる。

一つの波を指定すると位相点は $qx-\omega t=\mathrm{const}$ に沿って動く。その速さ $\omega/q$ を :term[位相速度]{en="phase velocity" yomi="いそうそくど"} という。一方、近い波数の干渉で作る包絡が動く速さ $\dv{\omega}{q}$ が :term[群速度]{en="group velocity" yomi="ぐんそくど"} である。分散が原点近くで直線なら両者は一致するが、曲率のある分枝では差が出る。調和鎖のエネルギー輸送でも、この違いは局所的なエネルギー保存則から確認できる。

第 $n$ サイトへ運動エネルギーと左右のばねのエネルギーの半分を割り振り、$h_n=M\dot u_n^2/2+K[(u_{n+1}-u_n)^2+(u_n-u_{n-1})^2]/4$ とする。$h_n$ を微分して運動方程式を代入すると、右向きの流れを

$$
 J_{n+1/2}=-\frac K2(\dot u_{n+1}+\dot u_n)(u_{n+1}-u_n)
$$

と置いたとき $\dot h_n=J_{n-1/2}-J_{n+1/2}$ が得られる。実波 $u_n=A\cos(qna-\omega t)$ では $u_{n+1}-u_n=-2A\sin(qna-\omega t+qa/2)\sin(qa/2)$ であり、速度の和は $2A\omega\sin(qna-\omega t+qa/2)\cos(qa/2)$ となる。時間平均を取ると $\overline J=KA^2\omega\sin(qa)/2$ である。時間平均の運動エネルギーは $M A^2\omega^2/4$ で、ばねからの寄与も運動方程式により同じ値になる。従って $\overline h/a=M A^2\omega^2/(2a)$ であり、$\overline J/(\overline h/a)=a K\sin(qa)/(M\omega)=\dv{\omega}{q}$ となる。熱伝導の式でエネルギーに群速度を掛けるのは、この保存則と整合する。

:::example[境界と長波長]
$q=\pi/a$ では隣接変位が反対位相で、$\omega=2\sqrt{K/M}$ だが群速度は零である。$q=\pi/(2a)$ では $\omega=\sqrt{2K/M}$、正の群速度は $a\sqrt{K/(2M)}$ である。
:::

:::example[固定点を設けると零振動数が消える]
各サイトを固定点へ定数 $K_0$ のばねで結ぶと、運動方程式の右辺に $-K_0u_n$ が入る。分散の二乗には一様に $K_0/M$ が加わる。$K_0=K$ なら最低角振動数は $\sqrt{K/M}$、最高角振動数は $\sqrt{5K/M}$ である。$q=0$ の一様変位でも全ての固定ばねが伸びるので、もとの零モードは消える。
:::

### 二原子鎖

一胞に $M_1,M_2$ の二原子を置き、胞間隔を $a$、隣の原子との距離を $a/2$ とする。振幅に胞内位置の位相を含めれば最近接ばねの二つの寄与は $2K\cos(qa/2)$ となる。変位振幅 $u,v$ に対し
$$
 \begin{pmatrix}2K-M_1\omega^2&-2K\cos(qa/2)\\
 -2K\cos(qa/2)&2K-M_2\omega^2\end{pmatrix}
 \begin{pmatrix}u\\v\end{pmatrix}=0.
$$
行列式を展開して $M_1M_2\omega^4-2K(M_1+M_2)\omega^2+4K^2\sin^2(qa/2)=0$ を得る。従って
$$
 \omega_\pm^2=K\left(\frac1{M_1}+\frac1{M_2}\right)
 \pm K\sqrt{\left(\frac1{M_1}+\frac1{M_2}\right)^2-\frac{4\sin^2(qa/2)}{M_1M_2}}.
 \label{eq-ss-diatomic-phonon}
$$
同じ行列の解析は式 :eqref[eq-am-diatomic-chain] にも現れる。$q=0$ で下側の音響分枝は $u=v$、上側の光学分枝は $M_1u+M_2v=0$ となる。小さい $q$ で定数項と $\omega^2$ 項を釣り合わせると $v_s=a\sqrt{K/[2(M_1+M_2)]}$ を得る。

:::example[質量差による振動数の間隔]
$M_2=4M_1$ なら境界で $\omega_-=\sqrt{K/(2M_1)}$、$\omega_+=\sqrt{2K/M_1}$ である。境界では非対角要素が零なので低い方は重い原子、高い方は軽い原子だけが動く。$M_2=M_1$ にすると、力学方程式が $a/2$ の並進で不変になる。従って分散は間隔 $a/2$ の単原子鎖を二倍の胞へ折り返した形になり、境界で二枝が一致する。原子の種類が異なっても、質量とこの模型のばね定数が同じなら、同じ分散が得られる。
:::

### 一次元状態密度と実際の音速

状態密度 $g(\omega)$ は、区間 $\dd\omega$ に含まれるモード数を $g(\omega)\dd\omega$ と表したものである。長さ $L=Na$ の周期鎖で $\dd q$ に含まれる波数の数は $L\dd q/(2\pi)$ である。正の $q$ に対する逆関数は $q=(2/a)\arcsin(\omega/\omega_{\max})$、$\omega_{\max}=2\sqrt{K/M}$ である。負の $q$ も同じ振動数に寄与するので両者を数えると

$$
 g(\omega)=\frac{2N}{\pi\sqrt{\omega_{\max}^2-\omega^2}},
 \qquad 0<\omega<\omega_{\max}.
 \label{eq-ss-chain-dos}
$$

$\omega=\omega_{\max}\sin\theta$ と置けば積分は $(2N/\pi)\int_0^{\pi/2}\dd\theta=N$ となる。境界で群速度が零になるため密度が発散するが、その近傍でも積分値は有限である。実際、$\omega_{\max}-\omega=\delta$ なら分母は $\sqrt{(2\omega_{\max}-\delta)\delta}$ であり、幅 $\delta$ の領域のモード数は $\sqrt{\delta}$ の大きさで減る。

:::example[振動数の低い側にあるモードの割合]
$\omega<\omega_{\max}/2$ にあるモードの割合は $(2/\pi)\arcsin(1/2)=1/3$ である。従って帯域の半分より下にあるモードは半数より少ない。この違いは端に向かって群速度が下がることから生じる。
:::

隣から二つ離れたサイトまで相互作用があるときは、ばね定数 $K_1,K_2$ を使って、力に $K_2(u_{n+2}+u_{n-2}-2u_n)$ を加える。平面波を代入して

$$
 M\omega^2=4K_1\sin^2(qa/2)+4K_2\sin^2(qa),
 \qquad
 v_s=a\sqrt{(K_1+4K_2)/M}
$$

を得る。第 2 項は境界で消えるが、小さい $q$ では距離の二乗に対応する因子 $4$ が付く。例えば $K_2=K_1/4$ なら境界の振動数は変わらず、音速は $\sqrt2$ 倍になる。反対に境界だけを測ると、この第 2 近接力定数は推定できない。

### 交互に異なる結合の鎖

質量が等しくても、結合が交互に異なると二原子胞が必要になる。第 $n$ 胞内の変位を $u_n,v_n$ とし、$u_n,v_n$ を $K_1$、$v_n,u_{n+1}$ を $K_2$ で結ぶ。胞番号にだけ位相 $e^{iqna}$ を含めると、$u$ へ働く胞間の力は $K_2(v e^{-iqa}-u)$ であり、同じ胞の力は $K_1(v-u)$ である。$v$ への力も同様に足せば

$$
 M\omega^2
 \begin{pmatrix}u\\v\end{pmatrix}
 =
 \begin{pmatrix}
 K_1+K_2 & -K_1-K_2e^{-iqa}\\
 -K_1-K_2e^{iqa} & K_1+K_2
 \end{pmatrix}
 \begin{pmatrix}u\\v\end{pmatrix}.
$$

共通の対角成分を引くと対角が零になるため、残った行列の固有値は非対角成分の絶対値の正負である。したがって

$$
 M\omega_\pm^2=K_1+K_2\pm\sqrt{(K_1+K_2)^2-4K_1K_2\sin^2(qa/2)}.
$$

この式では低い方の小波数展開が $\omega_-^2=a^2q^2K_1K_2/[2M(K_1+K_2)]+\order{q^4}$ になる。静的に一本の胞を伸ばしたときも、二本のばねに同じ力が通り、それぞれの伸びが足し合わさる。従って有効ばね定数は $K_{\rm eff}=K_1K_2/(K_1+K_2)$ である。質量 $2M$ の胞がこのばねでつながる長波長の運動を考えると、同じ音速が得られる。

:::example[二量体と弱い結合]
$K_1=3K$、$K_2=K$ なら境界で二つの固有値は $\omega_-^2=2K/M$ と $\omega_+^2=6K/M$ で、音速は $a\sqrt{3K/(8M)}$ である。$K_2$ を零に近づけると二量体の重心を移動させる費用がなくなるため、音響分枝は全波数で零へ近づく。内部振動の分枝は $\sqrt{2K_1/M}$ に残る。
:::

## 3 次元格子の基準振動と音響モード・光学モード

### 力学行列の固有値

胞内原子数を $r$、胞の総数を $N$ とする。胞内位相は波数に依存する基底の取り方で変えられる。ここでは
$u_{\vb R s\alpha}=M_s^{-1/2}e_{s\alpha}e^{i(\vb q\cdot\vb R-\omega t)}$ を選ぶ。運動方程式に代入して
$$
 D_{s\alpha,t\beta}(\vb q)=\sum_{\vb R}
 \frac{\Phi_{s\alpha,t\beta}(\vb R)}{\sqrt{M_sM_t}}e^{i\vb q\cdot\vb R},
 \qquad D(\vb q)\vb e_{\vb q j}=\omega_{\vb q j}^2\vb e_{\vb q j}
 \label{eq-ss-dynamical-matrix}
$$
となる。添字を交換し $\vb R\mapsto-\vb R$ とすれば $D^\dagger=D$ である。よって各波数で $3r$ 個の直交固有ベクトルを取れる。実変位を作るときは $\vb q$ と $-\vb q$ の振幅を共役にする。これは一般の基準振動の :ref[thm-am-normal-modes] を胞間 Fourier 変換でブロック化したものである。

### 長波長での区別

音響和則から $\vb e_{s}=\sqrt{M_s}\vb d$ は $\vb q=0$ における固有値零の固有ベクトルである。独立な $\vb d$ が三つあるので、他に零モードのない安定な結晶では音響分枝が三本、残りの光学分枝が $3r-3$ 本になる。長波長の音響モードは弾性変位場に対応し、音速は :eqref[eq-ss-elastic-energy] を用いた連続体の計算と一致する。光学モードの変位は質量規格化で並進モードと直交するため、$\sum_s M_s\vb u_s=0$ を満たす。イオンの有効電荷が異なると双極子が振動し、赤外電場に結合できる。全ての光学分枝が赤外活性になるわけではなく、分枝の対称性が行列要素を決める。

:::example[モードを数える]
ダイヤモンド型は原始胞に二原子なので六本のうち三本が光学分枝である。慣用立方胞の八原子を独立な原始胞とした計算では二十四本出るが、余分な本数は小さなゾーンへの折り返しで説明される。
:::

:::example[並進と内部運動の固有ベクトル]
一方向に動く二原子で $M_2=4M_1$ とする。$\vb q=0$ で質量規格化した並進固有ベクトルは $(1,2)/\sqrt5$ であり、直交する内部モードは $(2,-1)/\sqrt5$ である。後者の実変位比は $u_2/u_1=-1/4$ となる。
:::

### 有限格子の直交性と偏極の取り方

各方向に $N_i$ 胞を含む周期境界をとり、$N=N_1N_2N_3$ とする。許される波数は $\vb q=\sum_i (m_i/N_i)\vb b_i$ であり、$m_i$ は $N_i$ を法として一つずつ選ぶ。実格子の表示 $\vb R=\sum_i n_i\vb a_i$ を使えば、二つの平面波の内積は

$$
 \frac1N\sum_{\vb R} e^{i(\vb q-\vb q')\cdot\vb R}
 =\prod_i\frac1{N_i}\sum_{n_i=0}^{N_i-1}e^{2\pi i(m_i-m_i')n_i/N_i}
 =\delta_{\vb q,\vb q'}.
 \label{eq-ss-phonon-orthogonality}
$$

最後の式では幾何級数を用いた。分枝ベクトルを $\sum_{s,\alpha} e^*_{\vb q js\alpha}e_{\vb q j's\alpha}=\delta_{jj'}$ で規格化すると、$e_{\vb q js\alpha}e^{i\vb q\cdot\vb R}/\sqrt N$ は質量規格化した変位空間の直交系になる。列の数は $3rN$ で空間の次元と同じだから完全である。$\vb q$ の偏極の位相を変える自由度は、展開係数へ逆の位相を掛ければ物理的変位に影響しない。縮退した部分空間では、偏極をユニタリ混合しても同じ分枝空間を表す。

胞内座標を平面波に含める場合、振幅は $\tilde e_{s\alpha}=e^{-i\vb q\cdot\vb\tau_s}e_{s\alpha}$ と変換される。$U_{s\alpha,t\beta}=\delta_{st}\delta_{\alpha\beta}e^{-i\vb q\cdot\vb\tau_s}$ とすると $\tilde D=UDU^\dagger$ であり、固有値は同じである。この違いは散乱振幅を計算するときに現れる。例えば $\vb K=\vb G+\vb q$ での位相因子 $e^{i\vb K\cdot\vb\tau_s}e^*_{s\alpha}$ は、胞内位相を含めた偏極では $e^{i\vb G\cdot\vb\tau_s}\tilde e^*_{s\alpha}$ となる。

:::example[周期境界の波数点の数]
$N_1=2,N_2=3,N_3=4$ なら波数は $24$ 点ある。二原子胞で各点に六分枝があるから全モード数は $144$ で、$48$ 個ある原子の変位自由度と一致する。$q_z$ 方向の $m_3=2$ は $-2$ と合同であり、境界の二つの表示を別々に数えない。
:::

### 弾性波を離散模型から取り出す

単原子 Bravais 格子で Hessian の二次モーメントが絶対収束するとする。力定数の和則と反転対称性を使うと $D(\vb q)=M^{-1}\sum_{\vb d}' H(\vb d)[1-\cos(\vb q\cdot\vb d)]$ である。ここでプライムは原点を除く。したがって小波数で

$$
 D_{\alpha\beta}(\vb q)
 =\frac1{2M}\sum_{\vb d}'H_{\alpha\beta}(\vb d)
 (\vb q\cdot\vb d)^2+o(q^2).
$$

この展開は、和を取る各項で $(1-\cos x)/x^2\le1/2$ となるため、二次モーメントを支配関数とする収束から正当化できる。連続場で結合の変位差を $d_j\partial_j u_\alpha$ と近似すると、単位体積当たりのエネルギーは $(4v_c)^{-1}\sum_{\vb d}'H_{\alpha\beta}d_jd_l\,\partial_j u_\alpha\partial_l u_\beta$ となる。分母の $4$ は各結合の二次形式の係数 $1/2$ と、両向きの結合を数えた補正 $1/2$ の積である。この式の変分に $\rho=M/v_c$ を用いれば、平面波の固有値は上の $D(\vb q)$ と一致する。

:::example[軸方向のばねだけでは横波を支えられない]
単純立方格子の六最近接結合について $V'=0$、$V''=K$ とする。縦ばねだけがあるので $D_{xx}=2K[1-\cos(q_xa)]/M$、$D_{yy},D_{zz}$ は巡回置換である。$\vb q=(q,0,0)$ の横偏極の固有値は二つとも零である。この欠落を埋めるには、面対角の結合や結合角のエネルギーなど、せん断で二次から変化する相互作用を加える。
:::

:::example[三角格子の縦波と横波]
面内変位だけを許す二次元の単原子三角格子で、六最近接結合の長さを $a$、縦のばね定数を $K$ とする。代表三方向は $\vb n_1=(1,0)$、$\vb n_2=(1/2,\sqrt3/2)$、$\vb n_3=(-1/2,\sqrt3/2)$ である。張力が零なら結合の Hessian は $K\vb n_i\vb n_i^{\mathsf T}$ である。小波数で各方向とその反対方向をまとめると

$$
 D(\vb q)
 =\frac{Ka^2}{M}\sum_{i=1}^3
 \vb n_i\vb n_i^{\mathsf T}(\vb q\cdot\vb n_i)^2+o(q^2)
 =\frac{3Ka^2}{8M}
 \left(q^2I+2\vb q\vb q^{\mathsf T}\right)+o(q^2).
$$

最後の等式は、例えば $xx$ 成分で $q_x^2$ の係数が $1+1/16+1/16=9/8$、$q_y^2$ の係数が $3/16+3/16=3/8$、交差項が零となることから始めて各成分を足せば得られる。$xy$ 成分は $3q_xq_y/4$、$yy$ 成分は $3q_x^2/8+9q_y^2/8$ である。$\vb q$ に垂直な偏極には第一項だけが働き、平行な偏極には全体が働く。従って $v_T^2=3Ka^2/(8M)$、$v_L^2=9Ka^2/(8M)$ である。単純立方格子の軸方向ばねと異なり、斜めのばねが横ずれでも伸びて復元力を生む。
:::

多原子胞では平均変位がゆっくり変わるときに胞内原子も緩和する。平面波の弾性エネルギーを並進部分 $u$ と内部部分 $w$ に分け、力学行列のブロックを $D_{aa},D_{ao},D_{oa},D_{oo}$ と書く。光学モードにギャップがあれば小波数で $D_{oo}$ は可逆である。固有方程式の内部成分から $w=-(D_{oo}-\omega^2 I)^{-1}D_{oa}u$ を得る。音響分枝の $\omega^2$ が $q^2$ 程度、$D_{oa}$ が $q$ 程度であれば、二次までは逆行列の中で $\omega^2$ を無視できる。その結果、並進成分の行列は $D_{aa}-D_{ao}D_{oo}^{-1}D_{oa}$ である。第二項は非負の二次形式だから、内部緩和を許すと固定した胞内配置に対する剛性以下になる。

### 極性結晶の縦光学・横光学モード

正負の有効電荷 $\pm Z^*e$ を持つ二原子胞の長波長極限では、相対変位 $\vb w$ が分極 $\vb P_{\rm ion}=Z^*e\vb w/v_c$ を作る。短距離力だけでの相対振動数を $\omega_{\rm TO}$、換算質量を $\mu=M_1M_2/(M_1+M_2)$、電子の比誘電率を $\epsilon_\infty$ とする。運動方程式は $\mu\ddot{\vb w}=-\mu\omega_{\rm TO}^2\vb w+Z^*e\vb E$ である。自由電荷なしの静電近似では $\vb q\cdot(\epsilon_0\epsilon_\infty\vb E+\vb P_{\rm ion})=0$、$\vb q\times\vb E=0$ であるから、$\vb E=-\hat{\vb q}(\hat{\vb q}\cdot\vb P_{\rm ion})/(\epsilon_0\epsilon_\infty)$ となる。垂直偏極には電場が生じず、平行偏極には復元力が加わる。従って

$$
 \omega_{\rm LO}^2-\omega_{\rm TO}^2
 =\frac{(Z^*e)^2}{\epsilon_0\epsilon_\infty\mu v_c}.
 \label{eq-ss-lo-to}
$$

この長距離寄与は $\vb q$ の方向を含むので、$\vb q\to0$ の力学行列は一般に方向によって異なる。$\vb q=0$ に短距離の力定数をそのまま入れる計算では、電気的な境界条件を別途そろえる必要がある。例えば式 :eqref[eq-ss-lo-to] の右辺が $\omega_{\rm TO}^2/4$ なら、縦と横の振動数の比は $\sqrt5/2$ である。同じ有効電荷の二原子が逆向きに動く光学モードでも、双極子の寄与が相殺すればこの機構による分裂はない。[誘電性と光学的性質](/physics/solid-state/optical-dielectric/) では同じ復元力を誘電関数の零点として求める。

## フォノンの量子化

### 基準座標と Bose 統計

有限結晶で規格化した実固有ベクトルを質量規格化座標に使う。運動とポテンシャルの二次形式が対角になるため、並進零モードを除けば各モードは $H_\lambda=(P_\lambda^2+\omega_\lambda^2 Q_\lambda^2)/2$ である。$\comm{Q_\lambda}{P_{\lambda'}}=i\hbar\delta_{\lambda\lambda'}$ とし、
$$
 a_\lambda=\sqrt{\omega_\lambda/(2\hbar)}Q_\lambda+
 \frac{iP_\lambda}{\sqrt{2\hbar\omega_\lambda}}
$$
と置く。交換関係から $\comm{a_\lambda}{a_{\lambda'}^\dagger}=\delta_{\lambda\lambda'}$ となり、積を展開すると $H=\sum_\lambda\hbar\omega_\lambda(a_\lambda^\dagger a_\lambda+1/2)$ が得られる。量子一個を :term[フォノン]{en="phonon" yomi="ふぉのん"} という。振動子の構成は式 :eqref[eq-qm-oscillator] に従う。$\vb q$ と $-\vb q$ の縮退する実モードから進行波を作るユニタリ変換により、独立な生成演算子 $a^\dagger_{\vb q j}$ で記述してもよい。

温度 $T$ で占有数 $n$ の重みは $e^{-n x}$、$x=\hbar\omega/(k_BT)$ である。幾何級数 $Z=e^{-x/2}/(1-e^{-x})$ を微分すると $\bar n=1/(e^x-1)$ となる。これは :ref[def-stm-bose-gas] の分布で化学ポテンシャルを零とした形である。フォノン数が固体の粒子数のような固定量でないためである。変位は $a,a^\dagger$ に一次なので遷移行列要素は $\sqrt n,\sqrt{n+1}$ に比例する。

:::example[零点変位]
一質量 $M$ の振動子なら $\langle u^2\rangle=\hbar(2\bar n+1)/(2M\omega)$ である。これは $u=\sqrt{\hbar/(2M\omega)}(a+a^\dagger)$ の二乗で対角項だけを取った結果である。$k_BT\gg\hbar\omega$ では $k_BT/(M\omega^2)$ となり、古典的な等分配と一致する。
:::

:::example[結晶運動量の合同類]
一方向に $q=3\pi/(4a)$ のフォノンを二つ加えると合計は $3\pi/(2a)$ で、還元すると $-\pi/(2a)$ となる。これでエネルギーも保存されるかは個々の分散で別に判定する。
:::

### 変位演算子と交換関係

並進零モードをいったん固定し、正の振動数の部分を考える。式 :eqref[eq-ss-phonon-orthogonality] の規格化を使うと、変位演算子は

$$
 \hat u_{\vb R s\alpha}=\frac1{\sqrt{NM_s}}\sum_{\vb q,j}'\sqrt{\frac{\hbar}{2\omega_{\vb q j}}}
 \left(e_{\vb q js\alpha}e^{i\vb q\cdot\vb R}a_{\vb q j}
       +e_{\vb q js\alpha}^*e^{-i\vb q\cdot\vb R}a_{\vb q j}^\dagger\right).
 \label{eq-ss-displacement-operator}
$$

プライムは固定した並進モードを除くことを意味する。実の基準座標からこの展開へ移るには、$\vb q$ と $-\vb q$ が別点なら、二つの複素共役な基準ベクトルの和と差を $\sqrt2$ と $i\sqrt2$ で割って二つの実ベクトルを作る。二つの実振動子から逆のユニタリ変換で進行波の消滅演算子を作ると上式になる。$\vb q$ と $-\vb q$ が合同なら力学行列は実であり、一つずつ実の固有ベクトルを選べばよい。

共役運動量は $M_s\dot{\hat u}_{\vb R s\alpha}$ である。時間に関して $a(t)=a(0)e^{-i\omega t}$、$a^\dagger(t)=a^\dagger(0)e^{i\omega t}$ だから、交換関係では振動数とその逆数が打ち消し、残る和は完全系の射影である。零モードを除いた結果は内部変位に射影した交換関係である。並進運動を含めるなら結晶の重心位置 $\hat{\vb X}$ を全変位に加え、各運動量に $M_s\hat{\vb P}_{\rm tot}/(N\sum_t M_t)$ を加える。この加算の交換子は、欠けた一様並進成分をちょうど補う。

微小な時間差を比較するときは、変位どうしの相関が便利である。一モードの熱平衡では $\langle a^\dagger a\rangle=\bar n$、$\langle a a^\dagger\rangle=\bar n+1$ で、数を二つ変える項の平均は零になる。従って質量 $M$ の振動子では

$$
 \left\langle \hat u(t)\hat u(0)\right\rangle
 =\frac{\hbar}{2M\omega}\left[(\bar n+1)e^{-i\omega t}
 +\bar n e^{i\omega t}\right].
$$

正・負の周波数に非対称な重みが現れるのは、位置演算子が異なる時刻に可換でないためである。古典高温極限で $\bar n\gg1$ なら両者の重みは近くなり、相関は $[k_BT/(M\omega^2)]\cos\omega t$ に近づく。この非対称性は次節の中性子散乱で直接測れる。

:::example[零点変位の長さ]
$M=4.0\times10^{-26}\,\mathrm{kg}$、$\omega=1.0\times10^{13}\,\mathrm{s^{-1}}$ の振動子を考える。$\hbar\simeq1.055\times10^{-34}\,\mathrm{J\,s}$ とすれば零点の二乗平均平方根は $\sqrt{\hbar/(2M\omega)}\simeq1.15\times10^{-11}\,\mathrm m$ である。温度零でも原子は定まった位置の固有状態にない。$\bar n=1$ ならこの長さは $\sqrt3$ 倍になる。
:::

:::example[質量による振動数と変位の変化]
ばね定数を維持して質量を二倍にすると $\omega$ は $1/\sqrt2$ 倍になる。零点変位は $(M\omega)^{-1/2}$ なので $2^{-1/4}$ 倍になる。一方、高温古典極限の二乗平均は $k_BT/(M\omega^2)=k_BT/K$ であり、質量を変えても同じである。
:::

### 結晶運動量が表す対称性

実格子ベクトル $\vb a$ に対応する離散並進作用を $\hat{\mathcal T}_{\vb a}\hat u_{\vb R s}\hat{\mathcal T}_{\vb a}^{-1}=\hat u_{\vb R+\vb a,s}$ と定める。式 :eqref[eq-ss-displacement-operator] の両辺にこの作用を施すと、右辺の $a_{\vb qj}$ には $e^{i\vb q\cdot\vb a}$、$a_{\vb qj}^\dagger$ には $e^{-i\vb q\cdot\vb a}$ が掛かる。従って並進に不変な真空の上で生成したフォノンは、$\hat{\mathcal T}_{\vb a}$ の固有値 $e^{-i\vb q\cdot\vb a}$ を持つ。$\hbar\vb q$ を :term[結晶運動量]{en="crystal momentum" yomi="けっしょううんどうりょう"} という。$\vb q$ に逆格子ベクトルを加えても、すべての実格子ベクトルに対する固有値は変わらない。

相互作用の演算子 $\hat H_{\rm int}$ がこの離散並進と可換なら、固有位相が異なる状態間で行列要素は零である。実際、$\mel{f}{\hat H_{\rm int}}{i}$ の両側に $\hat{\mathcal T}_{\vb a}$ とその逆を挿入すると、行列要素は $e^{i(\vb q_f-\vb q_i)\cdot\vb a}$ 倍されても同じ値を持たなければならない。非零ならこの位相は全ての $\vb a$ について $1$ であり、$\vb q_f-\vb q_i$ は逆格子になる。多フォノン状態では生成演算子の個数だけ位相を掛けるため、比較すべき量は波数の和である。

質量を持つ結晶全体の実運動量は別の保存量である。式 :eqref[eq-ss-displacement-operator] を時間微分して $M_s$ を掛け、全ての $\vb R,s$ にわたって足す。$\vb q$ が零でなければ胞の位相和が消え、$\vb q=0$ の光学モードも並進固有ベクトルとの直交性により消える。したがって内部モードの共役運動量の総和は零である。自由な結晶の反跳を考えるときは、除いた重心自由度に実運動量が入る。

## 中性子非弾性散乱によるフォノン分散の測定

### 生成と消滅の条件

試料へ渡す波数を $\vb K=\vb k_i-\vb k_f=-\vb Q$、渡すエネルギーを $E=\hbar^2(k_i^2-k_f^2)/(2m_n)$ とする。核の散乱長を $b_s$ とすれば、核密度の散乱演算子は $\sum_{\vb R,s}b_s e^{i\vb K\cdot(\vb R+\vb\tau_s+\vb u_{\vb R s})}$ に比例する。$e^{i\vb K\cdot\vb u}=1+i\vb K\cdot\vb u+\cdots$ と展開すると定数項は弾性、一次項は一フォノン過程に対応する。

生成項では変位演算子に $a^\dagger_{\vb q j}e^{-i\vb q\cdot\vb R}$ が入るので胞間位相の和から $\vb K=\vb G+\vb q$ となる。エネルギー保存は $E=\hbar\omega_{\vb q j}$ である。消滅項では $\vb q$ の符号が逆になり $\vb K=\vb G-\vb q$、$E=-\hbar\omega_{\vb q j}$ となる。変位規格化因子は $\sqrt{\hbar/(2NM_s\omega_{\vb q j})}$ なので、一フォノン生成の強度は共通の運動学因子とエネルギー保存のデルタ関数を除き
$$
 \frac{\bar n+1}{\omega_{\vb q j}}
 \left|\sum_s\frac{b_s}{\sqrt{M_s}} e^{i\vb K\cdot\vb\tau_s}
 \vb K\cdot\conj{\vb e_{\vb q js}}\right|^2
$$
に比例する。熱変位をまとめて扱う場合は各項に Debye–Waller 因子を掛ける。消滅では $\bar n$ が入り、逆向き $\vb K$ と反転したエネルギー移行の強度比の熱因子は $\bar n/(\bar n+1)=e^{-E/(k_BT)}$ である。

:::example[散乱で見えない横波]
単原子で $\vb K=(K,0,0)$、偏極が $(0,1,0)$ なら強度は零である。測りたい $\vb q$ を変えずに別の $\vb G$ を選んで $\vb K$ に $y$ 成分を与えると可視になることがある。
:::

:::example[エネルギー損失と温度]
$\hbar\omega=k_BT\log2$ では $\bar n=1$ なので、対応する消滅と生成の熱因子の比は $1/2$ である。温度零では消滅が起きない一方、生成は中性子がエネルギーを渡せる範囲で可能である。
:::

### 動的構造因子と詳細釣り合い

散乱の強度は試料の固有振動数だけでは決まらない。入射波から出射波へ移るときの遷移振幅と、初期状態の熱分布も必要である。非磁性の単一散乱とし、核の散乱長 $b_s$ を実数とする。核の演算子を $\hat\rho_{\vb K}=\sum_{\vb R,s}b_s e^{i\vb K\cdot(\vb R+\vb\tau_s+\hat{\vb u}_{\vb R s})}$ と置く。試料の初期・終状態を $\ket{m},\ket{n}$、エネルギーを $E_m,E_n$、初期確率を $p_m=e^{-\beta E_m}/Z$ とすると、エネルギー $E$ に関する :term[動的構造因子]{en="dynamic structure factor" yomi="どうてきこうぞういんし"} を

$$
 S(\vb K,E)=\sum_{m,n}p_m|\mel{n}{\hat\rho_{\vb K}}{m}|^2\delta(E-E_n+E_m)
$$

と定義する。この規格化では、断面積 $\dd^2\sigma/(\dd\Omega\,\dd E_f)$ は $(k_f/k_i)S(\vb K,E_i-E_f)$ である。これを確かめるには、散乱ポテンシャルを $(2\pi\hbar^2/m_n)\sum_{\vb R,s}b_s\delta^{(3)}(\vb r-\vb R-\vb\tau_s-\hat{\vb u}_{\vb R s})$ とし、体積 $\mathcal V$ の平面波を代入する。散乱行列要素は $(2\pi\hbar^2/m_n\mathcal V)\mel{n}{\hat\rho_{\vb K}}{m}$ となる。Fermi の黄金律において終状態の状態密度 $\mathcal V m_n k_f/[(2\pi)^3\hbar^2]$ を掛け、入射粒子流 $\hbar k_i/(m_n\mathcal V)$ で割ると、係数は $k_f/k_i$ にまとめられる。$S$ の次元は面積をエネルギーで割ったものである。

同じ量は実時間の相関でも求められる。演算子 $\hat\rho_{-\vb K}(t)=e^{i\hat Ht/\hbar}\hat\rho_{-\vb K}e^{-i\hat Ht/\hbar}$ を展開し、間に完全系を入れると

$$
 S(\vb K,E)=\frac1{2\pi\hbar}\int_{-\infty}^{\infty}
 e^{iEt/\hbar}\langle\hat\rho_{-\vb K}(t)\hat\rho_{\vb K}(0)\rangle\dd t
 \label{eq-ss-phonon-dynamic-structure}
$$

を得る。指数因子の時間積分は $\delta(E-E_n+E_m)$ であり、$\hat\rho_{-\vb K}=\hat\rho_{\vb K}^\dagger$ だから行列要素の積は絶対値の二乗となる。$S(-\vb K,-E)$ の和で $m,n$ を入れ替えると、同じ二乗とデルタ関数に $p_n$ が掛かる。デルタ関数の台では $p_n=p_m e^{-\beta E}$ なので、:term[詳細釣り合い]{en="detailed balance" yomi="しょうさいつりあい"} は

$$
 S(-\vb K,-E)=e^{-\beta E}S(\vb K,E)
 \label{eq-ss-phonon-detailed-balance}
$$

となる。結晶の反転対称性を仮定しなくても、この関係では $\vb K$ と $E$ の両方を反転すればよい。実測値を比にする際には、試料の強度に掛かる $k_f/k_i$ と装置の効率を先に補正する。

### Debye–Waller 因子と多フォノン散乱

式 :eqref[eq-ss-phonon-dynamic-structure] の展開を調和結晶で整理する。原子 $l,m$ の平衡位置を $\vb r_l,\vb r_m$ と書き、$X=\vb K\cdot\hat{\vb u}_l(t)$、$Y=\vb K\cdot\hat{\vb u}_m(0)$ と置く。どちらも生成・消滅演算子の一次式なので交換子は数であり、熱平衡では Gaussian の特性関数を使える。$\langle e^{i(Y-X)}\rangle=\exp[-\langle(Y-X)^2\rangle/2]$ は、$Y-X$ を独立な振動子の線形結合に分けて各振動子の特性関数を掛ければ得られる。さらに $e^{-iX}e^{iY}=\exp[i(Y-X)+[X,Y]/2]$ により、指数に残る交差項は順序付きの $\langle XY\rangle$ となる。従って

$$
 \langle e^{-iX}e^{iY}\rangle
 =e^{-W_l-W_m}\exp(\langle XY\rangle),\qquad
 W_l=\frac12\left\langle(\vb K\cdot\hat{\vb u}_l)^2\right\rangle .
$$

$e^{-W_l}$ は振幅に掛かる :term[Debye–Waller 因子]{en="Debye–Waller factor" yomi="でばいわらーいんし"} である。Gaussian 性は正規順序からも確かめられる。一モードについて $e^{z a^\dagger-z^* a}=e^{-|z|^2/2}e^{z a^\dagger}e^{-z^*a}$ と変形する。右辺を展開し、数状態で平均すると同じ次数の生成と消滅の積だけが残る。幾何分布の下降階乗モーメントは、生成関数 $(1-e^{-x})/(1-te^{-x})$ を $t=1$ で $p$ 回微分して $p!\bar n^p$ となるので、平均は $e^{-(\bar n+1/2)|z|^2}$ である。

$\exp(\langle XY\rangle)$ を展開すると、定数項は時間因子を持たず弾性散乱になる。一次項は変位相関に含まれる $\pm\omega$ を選び、一フォノン散乱になる。二次以上は複数の振動数の和や差を含むので、多フォノン散乱になる。大きい $|\vb K|$ や高温で変位相関の大きさが増すと、同じ次数までの展開は精度が下がる。

散乱長が同位体や核スピンによってサイト間で独立にばらつく場合は $\overline{b_l b_m}=\overline b_l\,\overline b_m+\delta_{lm}\operatorname{Var}(b_l)$ と分かれる。平均値の積は別々の原子の位相を含む干渉性成分を与える。分散の項は一つのサイトの相関しか残さないため、胞の和に由来する $\vb K=\vb G\pm\vb q$ の制限がなくなる。等方的な粉末で方向平均すれば $|\vb K\cdot\vb e|^2$ の平均は $K^2|\vb e|^2/3$ であり、この成分は分枝を横断して質量と散乱長で重み付けした状態密度を反映する。

:::example[二原子胞の消滅則]
結合が交互に異なる一次元鎖を考え、等質量の原子を $\tau_1=0,\tau_2=a/2$ に置く。両原子の $b_s/\sqrt{M_s}$ を同じ値 $b$ とし、$W$ も等しいとする。$q=0$ の光学分枝は並進モードに直交するので、偏極を $(1,-1)/\sqrt2$ と取れる。一フォノン振幅に含まれる和は $bK[1-(-1)^h]/\sqrt2$ となる。ただし $K=2\pi h/a$ である。偶数 $h$ では消え、奇数 $h$ では二つの原子の寄与が同位相になる。同じ光学分枝でも測定する逆格子点で強度が変わる。
:::

### 散乱角と分解能

散乱角を $\theta$ とすると、固定した入射エネルギーで到達できる範囲は $K^2=k_i^2+k_f^2-2k_ik_f\cos\theta$ により決まる。$E=E_i-E_f$ を固定したときは $|k_i-k_f|\le K\le k_i+k_f$ である。この範囲外にある場合は、試料を回すだけでは目的の点に届かず、入射波長も変更する。三軸型では入射・出射の波数を結晶による反射で選び、試料の向きと組み合わせて $(\vb K,E)$ を移動する。飛行時間型では飛行時間と距離から出射速度を求め、検出器の方向と組み合わせて同じ座標に変換する。

:::example[運動学から決まる測定範囲]
中性子について $\hbar^2/(2m_n)\simeq2.072\,\mathrm{meV\,\mathring A^2}$ とする。$E_i=20\,\mathrm{meV}$ の中性子が $E=5\,\mathrm{meV}$ を試料に渡し、$\theta=60^\circ$ で検出されたとする。$k_i=\sqrt{20/2.072}\simeq3.11\,\mathrm{\mathring A^{-1}}$、$k_f=\sqrt{15/2.072}\simeq2.69\,\mathrm{\mathring A^{-1}}$ だから、$K\simeq2.92\,\mathrm{\mathring A^{-1}}$ である。このとき $K=0.20\,\mathrm{\mathring A^{-1}}$ には到達できない。両波数の差だけで約 $0.42\,\mathrm{\mathring A^{-1}}$ になるためである。
:::

有限の分解能を持つ測定は、理想的な $S$ を周囲の波数・エネルギーにわたって重み付き平均した値である。分散が局所的に $\delta E=\hbar\vb v_g\cdot\delta\vb q$ の直線で近似でき、装置の波数誤差とエネルギー誤差が独立な Gaussian なら、観測線のエネルギー分散は装置のエネルギー分散に $\hbar^2\vb v_g^{\mathsf T}\Sigma_q\vb v_g$ を加えたものになる。ここで $\Sigma_q$ は波数誤差の共分散である。従って急な音響分枝では、振動自体に減衰がなくても有限幅の線が現れる。例えば $\hbar v_g=30\,\mathrm{meV\,\mathring A}$、分散に沿う波数分解能 $0.02\,\mathrm{\mathring A^{-1}}$、エネルギー分解能 $0.80\,\mathrm{meV}$ を標準偏差で与えると、合成した幅は $1.00\,\mathrm{meV}$ である。

### 線幅と寿命の換算

弱い減衰により単一モードの振幅相関が $e^{-\Gamma|t|}e^{-i\omega_0t}$ で近似されるとする。$\Gamma$ は角振動数と同じ次元を持つ。式 :eqref[eq-ss-phonon-dynamic-structure] の時間積分を正負に分けて指数を積分すると、エネルギーの関数として

$$
 L(E)=\frac{1}{\pi}\frac{\hbar\Gamma}{(E-\hbar\omega_0)^2+(\hbar\Gamma)^2}
$$

が得られる。$E$ での全積分は $1$ で、分母の二項が等しいところが半値であるから、全半値幅は $2\hbar\Gamma$ である。振幅が $e^{-\Gamma t}$ で減るモードのエネルギーがその二乗に比例するなら、エネルギーの寿命は $\tau=1/(2\Gamma)$ である。この約束で全半値幅は $\hbar/\tau$ となる。振幅と占有数のどちらを寿命と呼ぶかによって係数が異なるので、線幅の値とともに約束を明記する。

:::example[装置の幅を除いて寿命を読む]
単一の Lorentz 型線について装置応答を分離した後の全半値幅が $0.20\,\mathrm{meV}$ なら、エネルギーの寿命は $\hbar/(0.20\,\mathrm{meV})\simeq3.29\,\mathrm{ps}$ である。群速度が $3000\,\mathrm{m\,s^{-1}}$ なら、その間に進む距離は約 $9.9\,\mathrm{nm}$ である。ただし熱流を緩和する散乱率には、運動方向がどれだけ変わるかも関わる。線幅が決める単一モードの寿命から熱伝導率まで求められるかは、散乱の種類の仮定に依存する。
:::

## 格子比熱（Einstein 模型と Debye 模型）

### 一振動数と連続分散

一モードの内部エネルギーは $\hbar\omega(\bar n+1/2)$ なので、固定体積で微分すると $k_B x^2e^x/(e^x-1)^2$ である。以下、胞の総数 $N$ と区別するため、原子の総数を $N_{\rm at}=rN$ と書く。Einstein 模型で $3N_{\rm at}$ モードを同じ $\omega_E$ とすれば、一モードの比熱を $3N_{\rm at}$ 倍する。$x\to0$ では $e^x-1\simeq x$ なので $C_V\to3N_{\rm at}k_B$、$x\to\infty$ では $C_V\simeq3N_{\rm at}k_B x^2e^{-x}$ である。ここで $x=\hbar\omega_E/(k_BT)$ である。

Debye 模型では三つの音響モードを等方的な $\omega=v_j q$ で近似する。体積 $V$ の周期境界の波数点密度は $V/(2\pi)^3$ なので、球殻 $4\pi q^2\dd q$ と $\dd q=\dd\omega/v_j$ により、一分枝の状態密度は $V\omega^2/(2\pi^2 v_j^3)$ である。等価な速さ $3/v_D^3=1/v_L^3+2/v_T^3$ を導入し、$3N_{\rm at}$ モードへ規格化する周波数遮断を $\omega_D=v_D(6\pi^2 N_{\rm at}/V)^{1/3}$ とすると
$$
 C_V=9N_{\rm at}k_B(T/\Theta_D)^3
 \int_0^{\Theta_D/T}\frac{x^4e^x}{(e^x-1)^2}\dd x,\qquad \Theta_D=\hbar\omega_D/k_B.
 \label{eq-ss-debye}
$$
多原子胞では低温係数は音速で決まり、全 $3N_{\rm at}$ を線形分散で近似する遮断は模型上の処方である。高温では被積分関数が $x^2$ に近づき $3N_{\rm at}k_B$ になる。これは一原子の三方向にそれぞれ運動と弾性の $k_BT/2$ がある :term[Dulong–Petit 則]{en="Dulong–Petit law" yomi="でゅろんぷてぃそく"} に一致する。低温では上限を無限大にできる。$e^x/(e^x-1)^2=\sum_{n\ge1}ne^{-nx}$ を項別積分すると $\sum_n 4!/n^4=4\pi^4/15$ である。各項は非負であり、$\int_0^\infty x^4e^{-nx}\dd x=4!/n^5$ は四回の部分積分で求まる。従って $C_V\simeq(12\pi^4/5)N_{\rm at}k_B(T/\Theta_D)^3$ となる。

級数の値 $\sum_{n\ge1}n^{-4}=\pi^4/90$ は Fourier 級数で計算できる。$x^2$ を区間 $[-\pi,\pi]$ で展開すると、定数項は $\pi^2/3$、第 $n$ 余弦係数は二回の部分積分により $4(-1)^n/n^2$、正弦係数は零になる。Parseval の等式を適用すると

$$
 \frac{1}{2\pi}\int_{-\pi}^{\pi}x^4\dd x
 =\frac{\pi^4}{9}+8\sum_{n\ge1}\frac1{n^4}.
$$

左辺は $\pi^4/5$ なので、整理すると必要な値が出る。また $\sum_{n\ge1}n^{-2}=\pi^2/6$ は $x$ の正弦係数 $2(-1)^{n+1}/n$ に同じ等式を適用して得られる。Parseval の等式は :ref[thm-four-circle-parseval] で示されている。

:::example[低温での温度比]
同じ結晶で $T,\ 2T\ll\Theta_D$ なら比熱の比は $8$ になる。有限結晶で最小の非零振動数が $k_BT/\hbar$ を超えるほど低温では、連続化が破れてこの比は成り立たない。
:::

:::example[二次元の音響分枝]
面積 $A$ の線形分散一分枝では球殻の代わりに $2\pi q\,\dd q$ が現れ、状態密度は $A\omega/(2\pi v^2)$ となる。従って低温比熱の積分で $\omega$ を $k_BT x/\hbar$ に変えると $T^2$ が出る。薄膜の曲げ分枝が二次分散なら、同じ数え方で状態密度は定数になり比熱は $T$ に比例する。
:::

### 分散の次数と低温比熱

$d$ 次元の等方的な一分枝で、小波数の分散を $\omega=a_s q^s$ と仮定する。$a_s>0$、$s>0$ とし、系の $d$ 次元体積を $V_d$、単位 $(d-1)$ 球面の面積を $\Omega_{d-1}$ とする。$q\le(\omega/a_s)^{1/s}$ の波数球の体積を数え、$\omega$ で微分すれば

$$
 g(\omega)=
 \frac{\Omega_{d-1}V_d}{(2\pi)^d s a_s^{d/s}}\omega^{d/s-1}
$$

となる。低温で他の極小や分枝からの寄与が小さいなら、$C_V$ の積分で $x=\hbar\omega/(k_BT)$ と変換して

$$
 \frac{C_V}{V_d}
 \simeq
 \frac{\Omega_{d-1}k_B}{(2\pi)^d s}
 \left(\frac{k_BT}{\hbar a_s}\right)^{d/s}
 \Gamma(2+d/s)\zeta(1+d/s)
$$

を得る。ここで $\Gamma(z)=\int_0^\infty x^{z-1}e^{-x}\dd x$、$\zeta(z)=\sum_{n\ge1}n^{-z}$ である。係数は積分 $\int_0^\infty x^{d/s+1}\sum_{n\ge1}ne^{-nx}\dd x$ で、各項の積分に $y=nx$ を使えば出る。低周波で $g$ が発散する場合も、$d/s>0$ なら比熱への寄与は積分可能である。例えば $d=2,s=2$ の曲げ分枝では $C_V/A\simeq\pi k_B^2 T/(12\hbar a_2)$ となる。ただし張力で分散の最低波数側が直線へ変わると、十分低温の次数も変わる。

一方、最低振動数 $\omega_0>0$ の光学分枝を考える。$\vb q_0$ が孤立した最小点で、その近くで $\omega\simeq \omega_0+A|\vb q-\vb q_0|^2$、$A>0$ とする。十分低温なら一モードの比熱は $k_B(\hbar\omega/k_BT)^2 e^{-\hbar\omega/(k_BT)}$ に近い。最低点の近くでは前因子の $\omega$ を $\omega_0$ に置き換え、$d$ 次元の波数を Gaussian で積分すると

$$
 \frac{C_{V,\rm opt}}{V_d}\simeq
 \frac{k_B}{(2\pi)^d}
 \left(\frac{\hbar\omega_0}{k_BT}\right)^2
 \left(\frac{\pi k_BT}{\hbar A}\right)^{d/2}
 e^{-\hbar\omega_0/(k_BT)}
$$

となる。熱励起される領域の波数幅が $\sqrt{k_BT/(\hbar A)}$ であることが $T^{d/2}$ を生む。平らな光学分枝を同じ $\omega_0$ の振動子とした Einstein 値では、この温度依存は現れない。どちらの場合も指数因子のため十分低温では音響分枝が優勢になる。

### 実在の分散から比熱を計算する

音響分枝と光学分枝の両方の振動数が分かる場合は、Debye 模型を介さずに $g(\omega)=\sum_{\vb q j}'\delta(\omega-\omega_{\vb q j})$ を作る。任意の関数 $f$ について $\int_0^\infty f(\omega)g(\omega)\dd\omega=\sum_{\vb q j}'f(\omega_{\vb q j})$ である。従って比熱は $C_V=\int_0^\infty k_B[x/(2\sinh(x/2))]^2 g(\omega)\dd\omega$ となる。高温では角括弧の二乗は $1-x^2/12+\order{x^4}$ なので、並進零モードの熱力学極限で

$$
 C_V=3N_{\rm at}k_B-\frac{\hbar^2}{12 k_BT^2}
 \int_0^\infty \omega^2 g(\omega)\dd\omega+\order{T^{-4}}
$$

を得る。固体の調和分散には上限があるから、$T$ を上げると $x$ は全分枝で一様に小さくなり、剰余を和の外に評価できる。従って高温の $3N_{\rm at}k_B$ への近づき方は、単一の遮断周波数よりも分散全体の二次モーメントを反映する。

異方的な結晶では音速が方向によって変わる。$\omega=v_j(\hat{\vb q})q$ と置き、球座標で $\dd^3q=q^2\dd q\,\dd\Omega$ を使うと

$$
 \frac{g(\omega)}V=\frac{\omega^2}{(2\pi)^3}
 \sum_{j=1}^{3}\int_{S^2}\frac{\dd\Omega}{v_j(\hat{\vb q})^3}
 +o(\omega^2)
$$

である。ここでは全方向で $v_j>0$ とし、線形項の剰余を方向について一様に制御できる場合を考えた。$v_j$ の逆三乗が入るため、遅い方向や遅い分枝は状態密度を大きく増やす。等方体に戻すと、一分枝の係数は $1/(2\pi^2v_j^3)$ になる。例えば $v_L=2v_T$ なら、横波二分枝の低温比熱は縦波一分枝の $16$ 倍である。

:::example[比熱から Debye 温度を読む]
原子一モル当たり $C_V/T^3\simeq1.94\times10^{-3}\,\mathrm{J\,mol^{-1}K^{-4}}$ と測定されたとする。低温を仮定し、気体定数を $R=8.314\,\mathrm{J\,mol^{-1}K^{-1}}$ とすると、$\Theta_D=[(12\pi^4/5)R/(1.94\times10^{-3}\,\mathrm{J\,mol^{-1}K^{-4}})]^{1/3}\simeq100\,\mathrm K$ である。化合物一モルの値を使う場合は、式量一単位中の原子数だけ $R$ を倍にする。金属で電子比熱 $C_{\rm el}=\gamma T$ もある場合は、低温の $C/T$ を $T^2$ に対して描く。切片が $\gamma$、傾きが格子の三乗係数になる。
:::

:::example[零点エネルギーと比熱を分ける]
Debye 模型の零点エネルギーを求めるには、$g(\omega)=9N_{\rm at}\omega^2/\omega_D^3$ に $\hbar\omega/2$ を掛けて $0$ から $\omega_D$ まで積分する。結果は $9N_{\rm at}\hbar\omega_D/8$ である。固定体積で温度を変えてもこの項は比熱に寄与しない。しかし $\omega_D$ が体積によって変わる場合、零点エネルギーは圧力に寄与する。次節ではこの体積微分を扱う。
:::

### 平均二乗変位と低周波の寄与

比熱は低周波の各モードから高々 $k_B$ の寄与を受けるが、変位の重みは別である。質量 $M$ の単原子格子に対し、式 :eqref[eq-ss-displacement-operator] を二乗して占有数で平均すると、全方向を足した一原子の二乗変位は

$$
 \langle|\hat{\vb u}|^2\rangle
 =\frac{\hbar}{2MN}\int_0^\infty \frac{g(\omega)}{\omega}
 \coth\!\left(\frac{\hbar\omega}{2k_BT}\right)\dd\omega .
$$

三次元の等方 Debye 模型では $T=0$ の値が $9\hbar/(4M\omega_D)$、高温での主要項が $9k_BT/(M\omega_D^2)$ となる。零温では積分に $\int_0^{\omega_D}\omega\,\dd\omega$、高温では $\coth y\simeq1/y$ により $\int_0^{\omega_D}\dd\omega$ が現れる。一方、線形分散の $d$ 次元格子では、非零温度で小振動数の被積分関数が $\omega^{d-3}$ に比例する。二次元では最小波数を $2\pi/L$ とした有限系で $\log L$ に比例し、一次元では $L$ に比例する。この発散は周期境界の零モードを取り去った後にも残る。従って、低次元で並進秩序を考える場合には、温度と系サイズをどの順序で扱うかが重要になる。

## 非調和効果と熱膨張・Grüneisen パラメータ

### 二つの導き方

一結合の $U=Ku^2/2+Au^3/6$ を安定な高次項のある局所展開とみなす。古典分布を $A$ の一次まで展開すると $\langle u\rangle=-(A/6k_BT)\langle u^4\rangle_0=-Ak_BT/(2K^2)$ である。最後に Gaussian の四次モーメント $3(k_BT/K)^2$ を使った。典型的な $A<0$ では伸びる側で斜面が緩く、平均位置が伸びる。調和近似だけでは対称分布なのでこの効果は出ない。

準調和近似では各体積で振動は調和的だが $\omega_\lambda=\omega_\lambda(V)$ とする。:term[Grüneisen パラメータ]{en="Grüneisen parameter" yomi="ぐりゅーないぜんぱらめーた"} は $\gamma_\lambda=-\partial\log\omega_\lambda/\partial\log V$ である。振動自由エネルギー
$$
 F_{\rm vib}=\sum_\lambda\left[\hbar\omega_\lambda/2+
 k_BT\log(1-e^{-\hbar\omega_\lambda/(k_BT)})\right]
$$
を微分すると $\partial F_\lambda/\partial\omega_\lambda=\hbar(\bar n_\lambda+1/2)$ である。従って $P_{\rm vib}=-\partial_V F_{\rm vib}=\sum_\lambda\gamma_\lambda U_\lambda/V$ となる。全圧力の微分 $0=(\partial_T P)_V\dd T-(B_T/V)\dd V$ から
$$
 \alpha_V=\frac1V\left(\pdv{V}{T}\right)_P
 =\frac1{B_T V}\sum_\lambda\gamma_\lambda C_{V,\lambda}
 \label{eq-ss-gruneisen}
$$
を得る。これは周波数の温度依存を体積経由で扱った式であり、固定体積でも起こる非調和な周波数ずれや寿命は別扱いとなる。

:::example[膨張係数の符号]
二分枝の比熱が同じで $\gamma_1=2,\gamma_2=-3$ なら和は負である。負の $\gamma$ の分枝の熱励起が優勢な温度域では体積は加熱により縮む。
:::

### 体積微分と実験量

モードの体積微分では、波数の絶対値でなく、逆格子基底に対する座標を固定する。結晶を一様に膨張しても同じ整数ラベルの波は系内で同じ回数だけ振動するためである。例えば単原子の等方的な音響波で $v\propto V^p$ と仮定する。固定ラベルでは $q\propto V^{-1/3}$ だから $\omega=vq\propto V^{p-1/3}$、従って $\gamma=1/3-p$ である。絶対値 $q$ を固定して微分すると、この式の $1/3$ が失われる。

非縮退分枝では、式 :eqref[eq-ss-dynamical-matrix] の固有方程式を体積で微分すると

$$
 \gamma_\lambda=-\frac{V}{2\omega_\lambda^2}\vb e_\lambda^\dagger
 \left(\pdv{D}{V}\right)\vb e_\lambda
$$

となる。実際、微分した式に左から $\vb e_\lambda^\dagger$ を掛けると、$\partial_V\vb e_\lambda$ に掛かる $D-\omega_\lambda^2 I$ が消える。縮退点では同じ式の右辺の行列を縮退部分空間で対角化して、体積変化により分かれる分枝を定める。変形前後の振動数を機械的に低い順に対応付けると、分枝が交差するところで微分を誤る場合がある。

熱容量の測定は一定圧力で行われる場合が多い。自由エネルギー $F(T,V)$ について $\dd F=-S\dd T-P\dd V$ だから $\partial_V S=\partial_T P$ である。従って

$$
 C_P-C_V
 =T\left(\pdv{S}{V}\right)_T\left(\pdv{V}{T}\right)_P
 =TVB_T\alpha_V^2 .
 \label{eq-ss-phonon-heat-capacity-difference}
$$

一つ目の等号は $S(T,V(T,P))$ を温度で微分した結果である。二つ目には $(\partial_T P)_V=B_T\alpha_V$ と $B_T=-V(\partial_VP)_T$ を用いた。膨張が負の場合でも差は非負になる。$\gamma_\lambda$ が低周波で有限、$B_T$ も有限の安定な三次元絶縁体なら、低温では :eqref[eq-ss-gruneisen] により $\alpha_V=\order{T^3}$ となり、差は $\order{T^7}$ である。この場合、低温の比熱から三乗係数を読む際に $C_P$ と $C_V$ の差は主要項に効かない。

:::example[測定された定圧比熱を補正する]
$T=300\,\mathrm K$、モル体積 $10^{-5}\,\mathrm{m^3\,mol^{-1}}$、$B_T=100\,\mathrm{GPa}$、$\alpha_V=3.0\times10^{-5}\,\mathrm{K^{-1}}$ とする。式 :eqref[eq-ss-phonon-heat-capacity-difference] は、モル当たりの差 $0.27\,\mathrm{J\,mol^{-1}K^{-1}}$ を与える。立方晶で測定値が線膨張係数 $\alpha_L$ なら、$V=L^3$ を微分して得る $\alpha_V=3\alpha_L$ を用いる。線膨張係数をそのまま代入すると差を九分の一に見積もる。
:::

### 異方的な熱膨張

結晶では方向ごとの伸びが異なり、温度変化により体積とともに形も変わる。主軸を保つ微小な垂直ひずみだけを扱い、軸 $i$ のひずみを $\epsilon_i$、引張応力を $\sigma_i$ とする。微分の基準を現在の平衡状態に取り、外部応力を一定に保つ。モードごとのひずみ応答を $\gamma_{\lambda i}=-\partial_{\epsilon_i}\log\omega_\lambda$ とすれば、振動自由エネルギーからの応力は $-\sum_\lambda \gamma_{\lambda i}U_\lambda/V$ である。これは $V^{-1}\partial_{\epsilon_i}F_{\rm vib}$ に、既に求めた $\partial_{\omega_\lambda}F_{\rm vib}$ を代入した式である。全自由エネルギーの等温の弾性行列を $C_{ij}=(\partial_{\epsilon_j}\sigma_i)_T$ とすると、固定応力で

$$
 0=\dd\sigma_i=\sum_j C_{ij}\dd\epsilon_j
 -\frac1V\sum_\lambda\gamma_{\lambda i} C_{V,\lambda}\dd T
$$

である。ここでは垂直方向とせん断の成分が対称性で分離する場合を考えた。逆行列 $S=C^{-1}$ を使えば、方向 $i$ の線膨張係数は $\dd\epsilon_i/\dd T=V^{-1}\sum_{j,\lambda}S_{ij}\gamma_{\lambda j}C_{V,\lambda}$ である。従って直角方向どうしの弾性結合を通して、他方向の熱応力も伸びに影響する。

:::example[方向によって符号が異なる膨張]
軸 $1,2$ の弾性部分を $C=\begin{pmatrix}2c&c\\c&2c\end{pmatrix}$、$\sum_\lambda \gamma_{\lambda i}C_{V,\lambda}/V=(f,0)_i$、$c,f>0$ とする。逆行列は $(3c)^{-1}\begin{pmatrix}2&-1\\-1&2\end{pmatrix}$ で、膨張係数は $(2f/(3c),-f/(3c))$ となる。二方向目に直接の熱応力がなくても、一方向目の膨張に伴って二方向目が縮む。
:::

### 固定体積での非調和性

固定体積でも振動の大きさに応じて振動数が変わることがある。一自由度のポテンシャルを $U=Ku^2/2+Bu^4/4!$ とし、$B>0$、振幅を $A$ とする。運動方程式は $M\ddot u+Ku+Bu^3/6=0$ である。$u=A\cos\omega t+\order{BA^3/K}$ を代入し、$\cos^3 t=(3\cos t+\cos3t)/4$ と分解する。基本周波数の成分を消すと $M\omega^2=K+BA^2/8+\order{B^2A^4/K}$ となる。三倍の周波数の成分は、小さい第三高調波を変位に含めれば消える。実際、$u$ に $C\cos3\omega t$ を加えると、その調和力の係数は $(K-9M\omega^2)C$ なので、一次では $C=BA^3/(192K)$ である。従って典型的な振幅が温度とともに増せば、体積を固定しても振動数の変化が観測され得る。

量子力学では、調和振動子の長さを $\ell_0=\sqrt{\hbar/(2M\omega_0)}$、$\omega_0=\sqrt{K/M}$ として、$u=\ell_0(a+a^\dagger)$ と書く。$(a+a^\dagger)^2\ket{n}$ のうち、$\ket{n-2},\ket n,\ket{n+2}$ の係数はそれぞれ $\sqrt{n(n-1)},2n+1,\sqrt{(n+1)(n+2)}$ である。このベクトルのノルムの二乗が $\langle u^4\rangle_n/\ell_0^4=3(2n^2+2n+1)$ になる。従って一次のエネルギー補正は $B\ell_0^4(2n^2+2n+1)/8$、$n$ から $n+1$ の間隔の補正は $B\ell_0^4(n+1)/2$ である。質量が小さく零点振動が大きいときは、温度零からの励起でも補正が残る。

:::example[硬化した振動数]
古典模型で $BA^2/K=0.16$ なら、基本周波数の補正は一次で $\omega/\omega_0-1=0.01$ である。量子模型で $B\ell_0^4/(\hbar\omega_0)=0.02$ なら、基底状態からの遷移の補正も一次で $0.01\omega_0$ である。前者は外から与えた振幅で決まり、後者には零点振動の長さが入る。
:::

## フォノン熱伝導と Umklapp 過程

### 群速度と緩和時間

局所温度に対応する平衡占有 $\bar n_\lambda[T(\vb r)]$ へ時間 $\tau_\lambda$ で戻ると仮定する。定常で小さい温度勾配なら移流と緩和から $\delta n_\lambda=-\tau_\lambda(\vb v_\lambda\cdot\grad T)\partial_T\bar n_\lambda$ となる。熱流密度は $V^{-1}\sum_\lambda\hbar\omega_\lambda\vb v_\lambda\delta n_\lambda$ なので
$$
 \kappa_{ij}=\frac1V\sum_\lambda C_{V,\lambda}v_{\lambda i}v_{\lambda j}\tau_\lambda.
$$
角度平均が等方なら $v_i v_j$ は $v^2\delta_{ij}/3$ に置き換わる。全モードで同一の $v,\ell=v\tau$ を取れば $\kappa=C_V v\ell/(3V)$ である。境界散乱で $\ell$ が温度非依存なら、低温の Debye 領域で $\kappa\propto T^3$ になる。

### 正常過程と逆格子への移行

三次相互作用の Fourier 展開で、二フォノンから一フォノンを作る項には $\sum_{\vb R}e^{i(\vb q_1+\vb q_2-\vb q_3)\cdot\vb R}$ が乗る。従って波数の条件は $\vb q_1+\vb q_2=\vb q_3+\vb G$ であり、同時に $\omega_1+\omega_2=\omega_3$ が必要である。$\vb G=0$ を正常過程、$\vb G\ne0$ を :term[Umklapp 過程]{en="Umklapp process" yomi="うむくらっぷかてい"} という。正常過程だけでは非折り返しの結晶運動量は緩和しない。熱流と結晶運動量が同じ量である必要はないので、正常過程も運動量を保つ分布の範囲で熱流を再配分する。

:::example[波数だけでは過程は決まらない]
一次元で $q_1=q_2=3\pi/(4a)$ なら $q_3=-\pi/(2a)$、$G=2\pi/a$ が可能である。しかし式 :eqref[eq-ss-mono-phonon] では $2\sin(3\pi/8)>\sin(\pi/4)$ なのでこの合体は起こらない。
:::

:::example[境界で制限された熱伝導]
$C_V/V=300\,\mathrm{J\,m^{-3}K^{-1}}$、$v=3000\,\mathrm{m\,s^{-1}}$、$\ell=0.10\,\mathrm{mm}$ なら $\kappa=30\,\mathrm{W\,m^{-1}K^{-1}}$ である。試料幅を半分にして平均自由行程も半分になるなら、同温度の値は $15\,\mathrm{W\,m^{-1}K^{-1}}$ になる。
:::

### 衝突積分の釣り合いと波数の緩和

衝突ごとに振動数と波数が保存される条件を調べると、どの分布が緩和するかが分かる。三つの異なるモード $1,2,3$ を選び、三次のエネルギーに $A_{123} a_3^\dagger a_1 a_2+A_{123}^* a_1^\dagger a_2^\dagger a_3$ が含まれるとする。前向きの行列要素は $A_{123}\sqrt{n_1n_2(n_3+1)}$ である。各モードの占有が相関しないと近似し、平均占有を $f_i$ とすると、黄金律による正味生成率は

$$
 I_{12\to3}=\frac{2\pi}{\hbar}|A_{123}|^2
 \delta(\hbar\omega_3-\hbar\omega_1-\hbar\omega_2)
 \bigl[f_1f_2(1+f_3)-(1+f_1)(1+f_2)f_3\bigr].
 \label{eq-ss-phonon-three-mode-rate}
$$

ここでは $A_{123}$ に波数の選択則を含めた。二つの入力が同じモードの場合は、独立な占有二つの積の代わりに $n_1(n_1-1)$ を使い、行列要素と組の数え方をそろえる。通常の熱平衡では $f_i/(1+f_i)=e^{-\beta\hbar\omega_i}$ である。式 :eqref[eq-ss-phonon-three-mode-rate] の角括弧を $(1+f_1)(1+f_2)(1+f_3)$ で割ると、エネルギー保存の下で零になる。

正常過程では $\vb q_3=\vb q_1+\vb q_2$ なので、平衡分布の指数を $\beta(\hbar\omega_i-\hbar\vb q_i\cdot\vb u)$ に置き換えても同じ釣り合いが成立する。$\vb u$ は速度の次元を持ち、分母が正となるほど小さいとする。置き換えた分布は、全結晶運動量が零でないまま正常過程に対して定常になり得る。対して $\vb q_1+\vb q_2-\vb q_3=\vb G\ne0$ のとき、釣り合いに $\exp(\beta\hbar\vb G\cdot\vb u)$ が残る。従って各衝突を静止した平衡分布へ個別に戻すように扱うと、正常過程だけの結晶運動量保存が失われる。モード間の緩和を含む衝突演算子全体の解法と、境界近くで流れが分布する場合の輸送理論は、空間依存の分布と保存量を同時に扱う必要があるため本章では導かない。結晶の熱輸送におけるこの区別は :cite[ziman-solids] を参照。

低温では $\vb G\ne0$ の条件が過程を制約する。三つの波数がすべて $|\vb q_i|<G_{\min}/3$ なら、三角不等式によりそれらの符号付き和の長さは $G_{\min}$ に届かない。ここで $G_{\min}$ は最短の非零逆格子ベクトルの長さである。零振動数は原点の近くにしかないと仮定すると、Umklapp の合体は少なくとも一つの高いエネルギーを伴う。波数が大きいのが出力側であっても、エネルギー保存により入力エネルギーの和が大きくなる。そのため低温の熱平衡で寄与する合体率には指数抑制が入る。実際のしきい値は振動数保存と偏極依存の行列要素を同時に満たす領域から決まり、逆格子の長さだけで決まるわけではない。

### 高温における三フォノン散乱

三フォノン散乱が熱抵抗を支配するときの温度依存を調べる。式 :eqref[eq-ss-phonon-three-mode-rate] で $f_i=f_i^0+\delta f_i$、$\delta f_i=f_i^0(1+f_i^0)\psi_i$ と置く。平衡では $P=f_1^0f_2^0(1+f_3^0)=(1+f_1^0)(1+f_2^0)f_3^0$ である。生成と消滅の積を一次まで展開して引くと、$\delta f_i$ の係数は入力側が $P/[f_i^0(1+f_i^0)]$、出力側がその負になる。従って角括弧の一次項は $P(\psi_1+\psi_2-\psi_3)$ である。

高温で主要なモードが $f_i^0\simeq k_BT/(\hbar\omega_i)$ を満たし、振動数と三次の行列要素を一定と近似する。すると $P$ は $T^3$、$f_i^0(1+f_i^0)$ は $T^2$ に比例する。従って $\delta f_i$ に作用する線形化衝突演算子の係数は $T$ に比例する。Umklapp が熱流に重なる保存モードを解消し、衝突演算子を解いた熱伝導率が有限になるなら、一定の温度勾配による $\delta f_i$ は $1/T$ に比例する。移流項に入る $\partial_T f_i^0$ が高温で一定だからである。この仮定の下で $\kappa\propto 1/T$ となる。境界や四フォノン散乱が無視できない場合や、温度による振動数の変化が大きい場合には、同じ次数は保証されない。

### 同位体散乱と試料境界

一様な調和結晶では、異なる基準振動は結合しない。質量分布が乱れると調和近似のままでも波が散乱され、進行方向が変わる。平均質量 $M$ の単原子格子に、各サイトで独立な質量偏差 $\delta M_l$ があるとする。変位の共役運動量を $\hat{\vb p}_l$ とすると、運動エネルギーは $(M+\delta M_l)^{-1}\hat{\vb p}_l^2/2$ である。一次まで展開すると散乱項は $-\sum_l\delta M_l\hat{\vb p}_l^2/(2M^2)$ となる。運動量展開には各振動子から $\sqrt{\omega}$ が入るので、一フォノンから別の一フォノンへの行列要素は $\sqrt{\omega\omega'}\,N^{-1}\sum_l(\delta M_l/M)e^{i(\vb q-\vb q')\cdot\vb R_l}$ に、偏極の内積と $\hbar$ を掛けた形になる。

平均をとると、質量偏差の積のうち同じサイトの積だけが残る。従って行列要素の二乗は $\omega\omega'\operatorname{Var}(\delta M/M)/N$ に比例する。静的な乱れなのでエネルギー保存により $\omega'=\omega$ であり、黄金律で終状態を足すと状態密度 $g(\omega)\propto\omega^2$ が掛かる。この積から三次元の低周波では $\tau_{\rm iso}^{-1}\propto\omega^4$ となる。これは弱い無相関の質量乱れと線形な音響分散を仮定した :term[Rayleigh 散乱]{en="Rayleigh scattering" yomi="れいりーさんらん"} の次数である。この散乱だけを全周波数に適用すると、熱伝導の積分は低周波で $\int_0 \omega^2\tau_{\rm iso}\dd\omega\propto\int_0\omega^{-2}\dd\omega$ と発散する。実際の有限試料では長波長の波の伝播も境界などで制限される。

:::example[質量の分散を調整する]
同位体の質量を $M_a,M_b$、存在割合を $p,1-p$ とすると、分散は $p(1-p)(M_a-M_b)^2$ である。質量差が平均質量の $0.10$ 倍で $p=1/2$ なら、無次元の質量分散は $0.0025$ である。同じ周波数・偏極で比較し、他の条件を保てば、$p$ を $0$ に近づけるにつれてこの項による散乱は弱くなる。
:::

境界が完全に拡散反射する場合は、試料の厚さの程度を平均自由行程とする模型が使える。波が境界から境界へほぼ鏡面反射する場合は、境界が運動方向を乱す程度が異なる。温度を下げて内部散乱が弱まると、測定した熱伝導率が試料寸法に依存することがある。三次元で境界の平均自由行程を一定とする Debye 領域では $T^3$ が見込まれる一方、さらに低温で横方向が量子化され、内部散乱が端から端まで起きない場合には、熱伝導率より端間の熱コンダクタンスが直接の量になる。

### 弾道輸送の一分枝

両端を温度 $T_L,T_R$ の熱浴に接続した一次元の伝播路を考える。単調な上向きの分散を持つギャップのない一分枝で、反射と内部散乱を無視する。右向きの状態の占有は $n(\omega,T_L)$、左向きの占有は $n(\omega,T_R)$ である。長さ当たりの右向き状態密度 $\dd q/(2\pi)$ に群速度 $v=\dd\omega/\dd q$ を掛けると、速度は変数変換で消える。従って熱流は

$$
 J_Q=\int_0^{\omega_{\max}}\frac{\dd\omega}{2\pi}\hbar\omega
 [n(\omega,T_L)-n(\omega,T_R)] .
$$

$k_BT_{L,R}\ll\hbar\omega_{\max}$ では上限を無限大としてよい。$\int_0^\infty x/(e^x-1)\dd x=\sum_{m\ge1}m^{-2}=\pi^2/6$ を使うと、$J_Q=\pi k_B^2(T_L^2-T_R^2)/(12\hbar)$ となる。左右の平均温度 $T$ の近くで一次化して、分枝一つ当たりの熱コンダクタンスは

$$
 G_Q=\frac{\pi k_B^2T}{6\hbar}
 \label{eq-ss-phonon-ballistic-conductance}
$$

である。単位は $\mathrm{W\,K^{-1}}$ で、拡散熱伝導率の単位 $\mathrm{W\,m^{-1}K^{-1}}$ とは異なる。例えば $T=1\,\mathrm K$ の理想的な一分枝なら $G_Q\simeq9.46\times10^{-13}\,\mathrm{W\,K^{-1}}$ である。複数分枝が開けばそれぞれの流れを足すが、分枝にギャップがあると低温での寄与は指数的に小さい。

## 演習問題

::::exercise{level=basic}
単原子鎖の $qa=\pi/3$ で $\omega/\sqrt{K/M}$ を求めよ。
:::solution
式 :eqref[eq-ss-mono-phonon] に代入すれば $2\sin(\pi/6)=1$ である。
:::
::::

::::exercise{level=basic}
三原子の原始胞をもつ安定結晶で、光学分枝の本数を求めよ。
:::solution
三方向に動けるので全体は $3\cdot3=9$ 本である。一様並進につながる三本を引くと $6$ 本である。
:::
::::

::::exercise{level=basic}
$\bar n=3$ のモードで消滅と生成の熱因子の比を求めよ。
:::solution
消滅には $\bar n$、生成には $\bar n+1$ が入るので $3/4$ となる。逆向きの波数移行を比べ、他の運動学因子は除く。
:::
::::

::::exercise{level=basic}
Debye 領域において温度を $1/3$ にしたときの比熱の比を求めよ。
:::solution
低温では式 :eqref[eq-ss-debye] の積分上限を無限大にでき、残る $T^3$ の比は $1/27$ である。
:::
::::

::::exercise{level=advanced}
二原子鎖で重い原子を $M_2=9M_1$ とする。振動数が最大の光学モードの $\omega^2$ を $K/M_1$ 単位で求めよ。
:::solution
根号は $\sin^2(qa/2)$ が最小のとき最大である。$q=0$ で $\omega_+^2=2K(1/M_1+1/M_2)=(20/9)K/M_1$ となる。
:::
::::

::::exercise{level=advanced}
力定数が実で式 :eqref[eq-ss-dynamical-matrix] を満たすとき、振動数の集合が $\vb q$ と $-\vb q$ で同じことを証明せよ。
:::solution
$D(-\vb q)=\conj{D(\vb q)}$ である。固有方程式を共役にし、$\omega^2$ が実であることを用いると、共役ベクトルが $-\vb q$ で同じ固有値を持つ。共役は逆変換でもあるため重複度も保たれる。
:::
::::

::::exercise{level=advanced}
全ての振動数が $V^{-2/3}$ に比例するとする。$B_T=50\,\mathrm{GPa}$、$C_V/V=1.5\times10^6\,\mathrm{J\,m^{-3}K^{-1}}$ のとき $\alpha_V$ を求めよ。
:::solution
対数微分により各モードの $\gamma=2/3$ である。式 :eqref[eq-ss-gruneisen] は $(2/3)C_V/(VB_T)$ となり、$2.0\times10^{-5}\,\mathrm{K^{-1}}$ である。
:::
::::

::::exercise{level=advanced}
線形分散の三次元一分枝で、緩和時間が $\tau=A/\omega$ とする。低温で $\kappa$ の温度依存を求めよ。
:::solution
状態密度は $\omega^2$ に比例し、速さは定数だから、積分は $\int_0^\infty \omega\,k_Bx^2e^x/(e^x-1)^2\,\dd\omega$ に比例する。$\omega=(k_BT/\hbar)x$ とすると $\kappa\propto T^2$ となる。無次元積分は零付近で被積分関数が $x$ に比例し、無限遠で指数減衰する。
:::
::::

::::exercise{level=advanced}
単原子鎖を弱く斜めに変形しても最近接のばね定数 $K_1$ は変わらず、第二近接のばね定数 $K_2$ だけ変わるとする。分散を測ると $q=\pi/a$ で変化がなく、音速が第二近接を含まないときの $3/2$ 倍になった。$K_2/K_1$ を求めよ。
:::solution
第二近接の寄与 $4K_2\sin^2(qa)/M$ は $q=\pi/a$ で零になるので、境界に変化がないことは条件と一致する。一方、音速の比は $\sqrt{(K_1+4K_2)/K_1}$ である。両辺を二乗して $1+4K_2/K_1=9/4$、したがって $K_2/K_1=5/16$ となる。
:::
::::

::::exercise{level=advanced}
入射エネルギー $25\,\mathrm{meV}$ の測定で、同じ大きさのエネルギー移行 $\pm9\,\mathrm{meV}$ と互いに反対の波数移行を選ぶ。$9\,\mathrm{meV}/(k_BT)=\log2$ とする。検出効率と分解能の違いを除き、エネルギー利得側の散乱断面積を損失側で割った比を求めよ。
:::solution
式 :eqref[eq-ss-phonon-detailed-balance] から $S(-\vb K,-E)/S(\vb K,E)=e^{-\log2}=1/2$ である。断面積の比には出射波数の比も掛かる。出射エネルギーは利得側で $25+9=34\,\mathrm{meV}$、損失側で $25-9=16\,\mathrm{meV}$ であるから、波数の比は $\sqrt{34/16}$ となる。従って求める値は $\sqrt{34}/8\simeq0.729$ である。動的構造因子の比だけを用いると出射側の状態密度の違いが抜ける。
:::
::::

## 確認テスト

::::quiz
調和結晶の音響和則は何に対応するか。

- [x] 全原子を同じ方向へ同じ距離だけずらしても内部エネルギーは不変である。
- [ ] 一原子だけをずらしたとき力が常に零になる。
- [ ] 各胞の光学モードが必ず赤外活性になる。

:::explanation
一様並進に対して二次形式が不変であるため、一定変位を力定数に掛けた結果が零になる。個々の原子を動かすと相対位置は変わる。
:::
::::

::::quiz{answer="2" tol="1e-6"}
単原子最近接鎖の最大角振動数を $\sqrt{K/M}$ 単位で答えよ。
:::explanation
第一ゾーンで $|\sin(qa/2)|$ の最大は境界での $1$ なので、係数は $2$ である。
:::
::::

::::quiz
Umklapp 過程について成り立つ条件を全て選べ。

- [x] エネルギー保存を満たす。
- [x] 波数の和は非零の逆格子ベクトルだけずれる。
- [ ] 波数を折り返さなくても和が等しい。

:::explanation
時間に依存しない格子相互作用によりエネルギーが保たれる。離散並進が要求するのは逆格子を法とした波数保存であり、非零の逆格子ベクトルが必要な場合を Umklapp と呼ぶ。
:::
::::

::::quiz{answer="8" tol="1e-6"}
十分低温の三次元 Debye 模型で、温度が二倍になると比熱は何倍になるか。
:::explanation
$C_V$ の低温極限は $T^3$ に比例するから $2^3=8$ 倍である。高温での極限は $3Nk_B$ なので、温度を二倍にする前後とも低温であることを用いた。
:::
::::

::::quiz{answer="4" tol="1e-6"}
完全透過の一次元音響分枝が同じ低温で四本開いているとき、熱コンダクタンスは一本のときの何倍か。
:::explanation
各分枝の右向きと左向きの流れの差を独立に足す。すべてギャップがなく同じ低温の極限にあるなら式 :eqref[eq-ss-phonon-ballistic-conductance] が分枝ごとに成り立つ。従って合計は四倍である。分枝の群速度は右向き状態密度と掛けたときに消えるので、同じ速度である必要はない。
:::
::::

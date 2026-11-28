---
title: 結晶中の電子の運動
description: バンド内の波束の運動から、Bloch 振動、サイクロトロン共鳴、Landau 量子化と磁気振動を導く。
notation:
  - symbol: 'E_n(\vb k)'
    meaning: バンド n のエネルギー
  - symbol: '(m^{*-1})_{ij}'
    meaning: 逆有効質量テンソル
  - symbol: 'A(E,k_\parallel)'
    meaning: 磁場に垂直な等エネルギー面の断面積
  - symbol: '\ell_B'
    meaning: 磁気長
---

## この章の目標

- 結晶中の波束の速度と、外場による波数の変化を導く。
- バンドの曲率から電子と正孔の慣性を求める。
- 電場中の Bloch 振動と磁場中の閉軌道を計算する。
- Landau 準位の縮退度と磁気振動の周期から Fermi 面の断面積を求める。

SI 単位を使い、電子の電荷を $q=-e$、$e>0$ とする。正孔との比較では $q=+e$ とする。特に断らない限り、スピンとバンド間遷移を省き、Berry 曲率や軌道磁気モーメントによる補正も無視する。

## 半古典的運動方程式

### 波束とゆっくり変わる外場

単一 Bloch 状態は結晶全体に広がるため、位置をもつ粒子として扱うには近い波数を重ねる。波束を波数 $\vb k_c$ の近傍 $\Delta k$ に集める。$\Delta k$ がゾーンの尺度より小さく、外場の変化尺度が波束幅 $1/\Delta k$ より大きいとする。エネルギーを $\vb k_c$ の周りで一次に展開すると、位相 $e^{i(\delta\vb k\cdot\vb r-E_n(\vb k)t/\hbar)}$ の包絡は $\vb r-\hbar^{-1}\grad_{\vb k}E_n\,t$ に依存する。従って速度は $\vb v_n=\hbar^{-1}\grad_{\vb k}E_n$ である。

弱い外場に対する中心の Lagrangian を
$$
 L=(\hbar\vb k+q\vb A)\cdot\dot{\vb r}-E_n(\vb k)-q\phi
$$
とする。外場なしでは $\hbar\vb k$ が波束の位置に共役であり、外場は中心の位置で評価する。この先頭次数の近似を、次の小節で波束の位相から確かめる。$\vb k$ の変分から群速度が出る。$\vb r$ の変分では $\dv{}{t}(q A_i)$ と $q\,\partial_i A_j\dot r_j$ の差が入り、$\vb E=-\grad\phi-\partial_t\vb A$、$\vb B=\rot\vb A$ を用いると
$$
 \dot{\vb r}=\frac1\hbar\grad_{\vb k}E_n,\qquad
 \hbar\dot{\vb k}=q(\vb E+\dot{\vb r}\times\vb B).
 \label{eq-ss-semiclassical}
$$
弱い場で孤立帯に留まることが前提である。$\dot E_n=\hbar\vb v\cdot\dot{\vb k}=q\vb v\cdot\vb E$ で、磁場は仕事をしない。

:::example[放物線による検算]
$E=\hbar^2k^2/(2m)$ なら $\vb v=\hbar\vb k/m$ となり、式 :eqref[eq-ss-semiclassical] は $m\dot{\vb v}=q(\vb E+\vb v\times\vb B)$ と一致する。周期帯では $\hbar\vb k$ を $m\vb v$ とみなせない。
:::

:::example[時間に依存する電場]
磁場なしで $E_x(t)=E_0\cos\omega t$、初期 $k_x(0)=0$ とする。積分して $k_x(t)=qE_0\sin(\omega t)/(\hbar\omega)$ である。この振幅がゾーン幅に比べて小さいかどうかにより、局所放物線近似を使える範囲が決まる。
:::

### 波束の中心と広がり

連続波数の記法を使い、孤立バンド $n$ の状態を $\Psi=\int C(\vb k)\psi_{n\vb k}\dd^3k$ とする。$\psi_{n\vb k}$ はデルタ関数で規格化し、$\int|C|^2\dd^3k=1$ とする。波束の分布がゾーン境界にかかる場合には、逆格子並進で同じ局所領域へ移して考える。係数の位相 $e^{-i\vb k\cdot\vb r_0}$ は包絡を $\vb r_0$ に移す。しかし周期部分 $u_{n\vb k}$ も波数とともに変化するので、波束の中心を厳密に求めるには、その変化も含める必要がある。

式 :eqref[eq-ss-bloch-form] を $\vb k$ で微分すると、位置演算子を Bloch 波に作用させた式は

$$
 \vb r\psi_{n\vb k}
 =-i\grad_{\vb k}\psi_{n\vb k}
  +i e^{i\vb k\cdot\vb r}\grad_{\vb k}u_{n\vb k}
$$

となる。第 1 項は波束の積分について部分積分できる。第 2 項は Bloch 位相が同じ波数にだけ射影され、その係数は胞内積になる。従って一帯へ射影した位置は係数に対し $i\grad_{\vb k}+\vb{\mathcal A}_n(\vb k)$ と作用する。ここで $\vb{\mathcal A}_n=i\langle u_{n\vb k},\grad_{\vb k}u_{n\vb k}\rangle_c$ であり、周期部分を胞内で規格化した内積を使う。波束の中心は

$$
 \vb r_c=\int \conj{C(\vb k)}
 (i\grad_{\vb k}+\vb{\mathcal A}_n)C(\vb k)\dd^3k.
 \label{eq-ss-packet-center}
$$

$u_{n\vb k}\mapsto e^{i\chi(\vb k)}u_{n\vb k}$ としたとき、同じ波束を表すには $C\mapsto e^{-i\chi}C$ とする。第 1 項が $\grad_{\vb k}\chi$ だけ増え、第 2 項は同量だけ減るので、中心は不変である。位置に現れる $\vb{\mathcal A}_n$ の微分は波数が時間変化する場合に速度を補正する。

波束の作用に現れる $i\hbar\langle\Psi,\partial_t\Psi\rangle$ を、実の包絡 $a(\vb k-\vb k_c)$ と位相 $e^{-i(\vb k-\vb k_c)\cdot\vb r_0}$ から求める。規格化により包絡の実微分の積分は零である。また $\int |C|^2(\vb k-\vb k_c)\dd^3k=0$ なので、位相の微分から残る項は $-\hbar\vb r_0\cdot\dot{\vb k}_c$ である。式 :eqref[eq-ss-packet-center] から狭い波束では $\vb r_c\simeq\vb r_0+\vb{\mathcal A}_n(\vb k_c)$ である。曲率が零の局所領域では $\vb{\mathcal A}_n$ を位相の選択で消去できる。従って全時間微分を加えて作用を書き直せば、位置に共役な項は $\hbar\vb k_c\cdot\dot{\vb r}_c$ となる。外場が緩やかなら、波束内ではベクトルポテンシャルを $\vb A(\vb r_c,t)$ として最小結合を使える。局所の平面波位相を指定する正準波数 $\vb K_c$ と分散に入る波数の関係は $\hbar\vb K_c=\hbar\vb k_c+q\vb A(\vb r_c,t)$ である。そこで作用の運動項を正準波数で書き、中心でのエネルギーを $E_n(\vb k_c)+q\phi(\vb r_c,t)$ とすれば、上の Lagrangian が得られる。これは場の変化と波束の幅についての先頭次数であり、内部の軌道磁気モーメントを加えたエネルギーは扱っていない。

外場がないとき $C(\vb k,t)=C(\vb k,0)e^{-iE_n(\vb k)t/\hbar}$ なので、式 :eqref[eq-ss-packet-center] の時間微分は $\langle\grad_{\vb k}E_n\rangle/\hbar$ である。中心波数で速度を代表させる近似は、包絡の波数幅にわたる速度変化を無視している。二次の分散は、中心速度を定めるだけでなく包絡を広げる。例えば一次元で $E(k)=\hbar^2 k^2/(2m^*)$、初期の $|C|^2$ が平均 $k_0$、分散 $\sigma_k^2$ の Gaussian で、位相に二次項がないとする。$u_k$ の変化が無視できる場合、初期の $\sigma_x=1/(2\sigma_k)$ である。自由発展下の位置は係数に対して $i\partial_k+\hbar kt/m^*$ となり、Gaussian の初期の対称共分散が零なので

$$
 \sigma_x^2(t)=\sigma_x^2(0)+
 \left(\hbar t\sigma_k/m^*\right)^2
$$

である。この式では包絡の広がりが変化しても、速度の一次補正は中心波数左右で打ち消し合う。

:::example[波束と外場の三つの長さ]
格子定数 $a=0.50\,\mathrm{nm}$、位置の標準偏差 $10\,\mathrm{nm}$ の最小不確定 Gaussian をとる。波数の標準偏差は $5.0\times10^7\,\mathrm{m^{-1}}$ で、第一ゾーン半幅 $\pi/a$ の約 $0.008$ 倍である。外場の変化尺度が $1\,\mathrm{\mu m}$ なら、位置分布の一標準偏差にわたる場の相対変化は $0.01$ 程度となる。この例では胞、包絡、場の変化という三つの長さが分離している。一方、場が表面の一胞で急変する場所には同じ近似を適用できない。
:::

### 一様電場による加速則の量子力学的な由来

空間一様な電場 $\vb E(t)$ に対して $\phi=0$、$\vb A(t)=-\int^t \vb E(t')\dd t'$ を選ぶ。このゲージでは $\vb A$ は空間に依存しないので、時間発展中も格子並進の固有値は保存する。保存される波数ラベルを $\vb K$ とし、状態を $e^{i\vb K\cdot\vb r}u(\vb r,t)$ と書く。最小結合した Hamiltonian をこの形に作用させると周期部分は

$$
 i\hbar\partial_tu=H_{\vb k(t)}u,\qquad
 H_{\vb k}=\frac{(-i\hbar\grad+\hbar\vb k)^2}{2m}+V,\qquad
 \vb k(t)=\vb K-\frac{q}{\hbar}\vb A(t)
$$

を満たす。格子並進の固有値が固定でも、瞬時の帯エネルギーに入る引数は変わる。式を微分すれば $\hbar\dot{\vb k}=q\vb E$ である。十分ゆっくり変化して孤立した帯の瞬時の固有状態を追う場合、速度は式 :eqref[eq-ss-bloch-velocity] になる。

断熱性を点検するには、$H_{\vb k}u_n=E_nu_n$ を時間で微分し、$m\ne n$ の状態で射影する。直交性から

$$
 \langle u_m,\partial_t u_n\rangle_c=
 \frac{\langle u_m,(\partial_t H_{\vb k})u_n\rangle_c}{E_n-E_m}
$$

を得る。左辺の絶対値がエネルギー差に対応する振動数 $|E_n-E_m|/\hbar$ に比べて小さいことが、一遷移の局所的な断熱性の目安である。結晶の場が弱くても、非常に狭いギャップではこの判定が破れ得る。有限次元に打ち切った Hermite 行列に対する一様な評価は :ref[thm-aqm-adiabatic] を参照する。

### 外場のエネルギーと位相空間の体積

波束に外場を入れる Lagrangian の項 $q\vb A\cdot\dot{\vb r}-q\phi$ は、電磁ゲージを $\vb A\mapsto\vb A+\grad\Lambda$、$\phi\mapsto\phi-\partial_t\Lambda$ と変えると $q\dd\Lambda/\dd t$ だけ変化する。作用の端点を固定した変分ではこの差が消える。従って式 :eqref[eq-ss-semiclassical] の $\vb k$ はゲージに依存しない量であり、位置の正準運動量 $\hbar\vb k+q\vb A$ とは区別する。

Berry 曲率を無視し、$E_n$ が $\vb k$ だけに依存する場合、位相空間での流れは体積を保存する。実際、$\grad_{\vb r}\cdot\dot{\vb r}=0$ であり、$\vb k$ についての発散は $(q/\hbar)\grad_{\vb k}\cdot(\vb v\times\vb B)$ に限られる。ここで $\vb B$ は固定した $\vb r$ では $\vb k$ に依存せず、$\rot_{\vb k}\vb v=0$ である。成分で書けば対称な Hessian と反対称な Levi-Civita 記号の縮約が零となる。衝突がなければ軌道に沿って動く微小体積内の粒子数も保存する。粒子数と体積がともに一定なので、分布 $f(\vb r,\vb k,t)$ は軌道上で一定である。この結果は [輸送現象](/physics/solid-state/transport/) の運動方程式の左辺に使う。

:::example[二つのゲージで仕事の式を確かめる]
静電ポテンシャルを $\phi(x)=-E_0x$ とする。全エネルギー $E_n(k)-qE_0x$ の時間微分には、式 :eqref[eq-ss-semiclassical] から $qE_0v_x$ と $-qE_0v_x$ が現れて相殺する。時間依存するベクトルポテンシャルで同じ電場を表すと、瞬時の帯の波数 $k$ はそのポテンシャルから直接求まる。この場合も $\dot E_n=qE_0v_x$ となる。
:::

### 投影した位置の交換関係

式 :eqref[eq-ss-packet-center] で使った一帯内の位置を $X_i=i\partial_{k_i}+\mathcal A_i$ とする。波数領域の内部に台を持つ滑らかな係数に作用させると、二階微分と $\mathcal A_i\mathcal A_j$ の項が交換子で相殺し、$[X_i,X_j]=i(\partial_i\mathcal A_j-\partial_j\mathcal A_i)=i\Omega_{ij}$ が残る。実空間の位置成分は可換だが、ほかの帯への成分を除いた演算子では可換性が失われ得る。曲率の定義と多帯への拡張は [トポロジカル物質入門](/physics/solid-state/topological/) に譲る。

一様な静電場をこの帯へ投影した Hamiltonian は $E_n(\vb k)-q\sum_j E_j X_j$ である。$[X_i,E_n]=i\partial_i E_n$ と上の交換子を代入すると、

$$
 \dot X_i=\frac1\hbar\partial_iE_n-\frac{q}{\hbar}\sum_j E_j\Omega_{ij},
 \qquad \hbar\dot k_i=qE_i .
$$

初めの項は群速度であり、曲率の項は電場に垂直である。ここでの結果は、一帯の幾何が電場下の速度に入ることを、波束の位置を定義する同じ演算子から示している。磁気モーメントや帯間遷移を含めるには投影の近似を拡張する必要がある。

:::example[投影と位置の順番]
$\Omega_{xy}=\Omega_0$ を一定とみなせる小領域で、電子に正の $x$ 向き電場を加える。$[X_y,X_x]=-i\Omega_0$ から $\dot Y$ の曲率による補正は $-eE_x\Omega_0/\hbar$ である。例えば $\Omega_0=10\,\mathrm{\mathring A^2}$、$E_x=10^5\,\mathrm{V\,m^{-1}}$ なら大きさは約 $15\,\mathrm{m\,s^{-1}}$ である。$\Omega_0$ の単位は面積なので、電場を掛けて電子の電荷と $\hbar$ を入れると速度になる。
:::

## 有効質量テンソルとホール

### 曲率と正孔

$\dot v_i=\sum_j(\partial v_i/\partial k_j)\dot k_j$ より、:term[逆有効質量テンソル]{en="inverse effective mass tensor" yomi="ぎゃくゆうこうしつりょうてんそる"} は
$$
 (m^{*-1})_{ij}=\frac1{\hbar^2}\partial_{k_i}\partial_{k_j}E_n
 \label{eq-ss-effective-mass}
$$
である。曲率の主軸ではこの行列は対角となる。主曲率が零でない方向では $F_i=m_i^*\dot v_i$ と書ける。主曲率が零の方向では、逆質量テンソルを使うと、その方向への加速度が零であることが分かる。バンド上端で電子の曲率が負なら、ほぼ満たされた帯の応答を空いた状態、すなわち :term[正孔]{en="hole" yomi="せいこう"} で記述する。抜いた電子の波数を $\vb k_e$ とし $\vb k_h=-\vb k_e$、$E_h(\vb k_h)=E_v-E_e(-\vb k_h)$ と置く。すると $\vb v_h=\vb v_e$、$\hbar\dot{\vb k}_h=+e(\vb E+\vb v_h\times\vb B)$ である。充填帯の電流が零なので、電子の欠落による電流は $+e\vb v_h$ となる。また、上端に対応する点で正孔の曲率は電子と逆符号になる。

一定の質量テンソルと緩和時間を用い、キャリア密度を $n_c$、電荷を $q$ とすると、定常運動は $0=q(\vb E+\vb v\times\vb B)-m^*\vb v/\tau$ である。$\vb j=q n_c\vb v$ を代入すれば $\vb E=m^*\vb j/(n_cq^2\tau)-\vb j\times\vb B/(n_cq)$ となる。主軸で縦に流して横電流を止めると $R_H=1/(n_cq)$ である。一般の非放物線帯ではこの一定質量の導出を適用できない。

:::example[斜め向きの加速]
主質量が $m_x=m,\ m_y=2m$ で $F_x=F_y=F$ とすると加速度は $(F/m,F/(2m))$ である。力の方向と加速度の方向は一致しない。
:::

:::example[帯底と帯頂の曲率]
$E=-2t\cos ka$、$t>0$ の底では $m_e^*=\hbar^2/(2ta^2)$ で、頂上の電子の質量はその負である。頂上付近に空きが少数あるなら正孔の質量は $+\hbar^2/(2ta^2)$ である。
:::

### バンド間の混合から曲率を求める

逆質量テンソルを $E_n$ の差分だけで求める代わりに、ある波数 $\vb k_0$ の固有状態から計算することもできる。局所スカラーポテンシャルを仮定すると、胞 Hamiltonian は $\vb k=\vb k_0+\delta\vb k$ において

$$
 H_{\vb k}=H_{\vb k_0}
  +\frac{\hbar}{m}\delta\vb k\cdot
    (-i\hbar\grad+\hbar\vb k_0)
  +\frac{\hbar^2|\delta\vb k|^2}{2m}
$$

と厳密に展開できる。$\vb P=-i\hbar\grad+\hbar\vb k_0$ とし、$P^i_{nm}=\mel{u_n}{P_i}{u_m}$ を胞積分で定める。孤立した非縮退エネルギーのまわりで [二次の非縮退摂動](/physics/solid-state/band-theory/) を適用すると、線形項の二次補正は $\frac{\hbar^2}{m^2}\sum_{m'\ne n}|\sum_i\delta k_i P^i_{m'n}|^2/(E_n-E_{m'})$ である。これと明示的な二次項をそれぞれ二階微分して

$$
 (m^{*-1})_{ij}
 =\frac{\delta_{ij}}{m}
 +\frac{2}{m^2}\Re\sum_{m'\ne n}
 \frac{P^i_{nm'}P^j_{m'n}}{E_n-E_{m'}}
 \label{eq-ss-kp-mass}
$$

を得る。$m$ は裸の電子質量、$m'$ は帯の添字である。上の帯が分母に現れる項は負、下の帯が分母に現れる項は正の符号を持つ。例えば伝導帯の上端側に曲率を下げる寄与があっても、すぐ下の価電子帯との結合は逆質量を増やす方向に働き得る。縮退点ではこの和を分母零のまま使わず、縮退帯を含む小さな行列を先に対角化する。

:::example[二帯を対角化して曲率を検算する]
二帯行列を $H(k)=\begin{pmatrix}\Delta/2&\beta k\\\beta k&-\Delta/2\end{pmatrix}+\hbar^2 k^2 I/(2m)$、$\Delta>0$ とする。正の帯の平方根は $\sqrt{\Delta^2/4+\beta^2 k^2}=\Delta/2+(\beta^2/\Delta)k^2+\order{k^4}$ だから、逆質量は $1/m+2\beta^2/(\Delta\hbar^2)$ である。式 :eqref[eq-ss-kp-mass] でも、下の帯との結合は $P=\beta m/\hbar$、分母は $\Delta$ となり一致する。
:::

### テンソルの軸を測定座標へ移す

二階微分が連続なら逆質量は実対称で、回転行列 $R$ によって対角化できる。主曲率を $\lambda_i=1/m_i^*$ とすると、単位方向 $\vb n$ に力をかけたとき、その方向への加速度は $F\,\vb n^{\mathsf T}m^{*-1}\vb n$ である。すなわち軸とずれた測定の逆質量は $\sum_i\lambda_i(\vb n\cdot\hat{\vb e}_i)^2$ になる。曲率の符号が混在する鞍点ではこの値が零となる方向がある。そのときも一般には加速度ベクトル全体が零にはならない。

:::example[鞍点で力と直交する加速度]
$E=\hbar^2(k_x^2-k_y^2)/(2m)$ の点で、力を $\vb F=F(1,1)/\sqrt2$ とする。加速度は $F(1,-1)/(\sqrt2 m)$ で、加速度の力方向への射影は零になる。方向微分に基づく一個のスカラー質量だけではこの状態を表せない。
:::

### Hall 測定で奇成分を取り出す

放物線状の単一キャリア帯で、質量テンソルが一定、緩和時間がスカラー一定であるとする。このとき抵抗率は $m^*/(n_cq^2\tau)$ の対称部分と $-(\vb j\times\vb B)/(n_cq)$ の反対称部分に分かれる。測定軸が主軸とずれていれば、磁場がなくても横方向に電圧が出る。磁場を反転させた二つの測定の差を半分にすると、磁場偶の対称成分が消え、磁場奇の Hall 成分が残る。式 :eqref[eq-ss-drude-hall] はこの単一帯の等方の場合に当たる。

ポケットが複数ある場合、各バンドで同じ電場を受け、電流を足す。抵抗率を足すのではなく、伝導度テンソルを足してから逆行列を取る理由はこの共通電場にある。電子密度 $n$、正孔密度 $p$、正の移動度 $\mu_e,\mu_h$ をもつ等方的な二帯では、弱磁場で $\sigma_{xx}=e(n\mu_e+p\mu_h)$ である。式 :eqref[eq-ss-drude-hall] の各帯の反転を磁場一次まで展開すれば $\sigma_{yx}=eB(n\mu_e^2-p\mu_h^2)$ となる。$j_y=0$ から $E_y=-\sigma_{yx}E_x/\sigma_{xx}$ であり、

$$
 R_H=\frac{p\mu_h^2-n\mu_e^2}{e(n\mu_e+p\mu_h)^2}
$$

を得る。密度を相互に打ち消す場合でも移動度が異なると Hall 係数は残る。強磁場の式は [輸送現象](/physics/solid-state/transport/) で扱う。

:::example[密度が同じ二種類のキャリア]
$p=n$、$\mu_h=2\mu_e$ とすると分子は $3n\mu_e^2$、分母は $9en^2\mu_e^2$ となる。従って $R_H=1/(3en)>0$ である。この符号は電子が存在しないことを表すのではなく、正孔の移動度二乗に比例する寄与が大きいことを表す。
:::

### 非放物線帯の弱磁場 Hall 応答

単一帯でも、軌道上で曲率の大きさや符号が変わると一定質量の導出は使えない。ここでは二次元、時間反転対称な非磁性帯、面に垂直な $B>0$ を仮定し、Berry 曲率の効果を無視する。$f_0(E)$ を平衡占有、$w=-\partial_E f_0$、$\tau(\vb k)>0$ を緩和時間とする。$\vb E$ に一次の定常分布の変化 $\delta f$ は、緩和時間近似で

$$
 (\tau^{-1}+\mathcal D_B)\delta f=q w\,\vb E\cdot\vb v,\qquad
 \mathcal D_B=\frac{qB}{\hbar}(v_y\partial_{k_x}-v_x\partial_{k_y})
$$

を満たす。$\mathcal D_BE=0$ なので $\mathcal D_Bw=0$ である。$\mathcal D_B$ が $\tau^{-1}$ に比べて小さい範囲で二回代入し、$\vb j=q\int\vb v\delta f\,\dd^2k/(2\pi)^2$ とすれば、一スピンの Hall 項は

$$
 \sigma_{yx}^{(1)}
 =-q^2\int\frac{\dd^2k}{(2\pi)^2}\,w\,(\tau v_y)\mathcal D_B(\tau v_x).
$$

温度零の $w=\delta(E-\mu)$ では、$\nabla_{\vb k}E$ が零でない閉軌道に沿って積分できる。磁場中の波数の速度は $|\dot{\vb k}|=|q|B|\nabla_{\vb k}E|/\hbar^2$ なので、弧長 $\dd l$ に対して $\dd l/|\nabla_{\vb k}E|=|q|B\,\dd t/\hbar^2$ である。これと軌道上の $\mathcal D_B=\dd/\dd t$ を代入する。$\boldsymbol\ell=\tau\vb v$ を平均自由行程ベクトルとすると

$$
 \sigma_{yx}^{(1)}
 =\frac{q^2|q| B}{(2\pi)^2\hbar^2}
 \sum_{\rm orbits}\oint \ell_x\,\dd\ell_y .
 \label{eq-ss-hall-mean-free-path}
$$

部分積分で $-\oint\ell_y\dd\ell_x=\oint\ell_x\dd\ell_y$ とした。周回は磁場中で実際に進む向きで取る。右辺は $\boldsymbol\ell$ の先端が描く符号付き面積の和である。曲線が自分と交差する場合は、周回の向きの逆な部分が差し引かれる。結晶軸を主な伝導軸に選べば、弱場で $\rho_{yx}^{(1)}=-\sigma_{yx}^{(1)}/(\sigma_{xx}^{(0)}\sigma_{yy}^{(0)})$ である。従って Hall 電圧の符号をキャリア密度だけから決められるのは、帯や散乱に追加の仮定がある場合である。

:::example[円形電子ポケットの符号と係数]
$q=-e$、$E=\hbar^2k^2/(2m)$ の電子軌道は、$B>0$ で $+k_x$ 軸から $+k_y$ 方向に回る。一定 $\tau$ なら $\oint\ell_x\dd\ell_y=\pi(\tau\hbar k_F/m)^2$ である。式 :eqref[eq-ss-hall-mean-free-path] は一スピン密度 $n=k_F^2/(4\pi)$ に対し $\sigma_{yx}^{(1)}=ne^3 B\tau^2/m^2$ を与える。$\sigma_{xx}^{(0)}=\sigma_{yy}^{(0)}=ne^2\tau/m$ を使うと、抵抗率では $R_H=-1/(ne)$ である。
:::

## 一様電場中の Bloch 振動

### ゾーンを一周する時間

一次元で $E_{\rm field}$ を一定とすると $k(t)=k_0+qE_{\rm field}t/\hbar$ である。一バンドの速度は逆格子周期 $2\pi/a$ をもつため :term[Bloch 振動]{en="Bloch oscillation" yomi="ぶろっほしんどう"} の周期と角振動数は
$$
 T_B=\frac{2\pi\hbar}{|q|\,|E_{\rm field}|a},\qquad
 \omega_B=\frac{|q|\,|E_{\rm field}|a}{\hbar}.
$$
式 :eqref[eq-ss-semiclassical] の仕事の式を積分すれば $x(t)-x(0)=[E_n(k(t))-E_n(k_0)]/(qE_{\rm field})$ となる。エネルギーは有界なので変位も有界である。余弦帯の最大値と最小値の差は $4t$ だから、最大振幅の半幅は $2t/(|q|\,|E_{\rm field}|)$ になる。

:::example[時間と長さの尺度]
$a=0.50\,\mathrm{nm}$、$|E_{\rm field}|=10^6\,\mathrm{V\,m^{-1}}$ なら $e|E_{\rm field}|a=0.00050\,\mathrm{eV}$ であり、$T_B\simeq8.27\,\mathrm{ps}$ である。通常の弾性散乱時間がこれより短ければ一周期を見る前に波束が緩和する。
:::

:::example[振幅の見積もり]
帯幅 $4t=0.10\,\mathrm{eV}$、$|E_{\rm field}|=10^7\,\mathrm{V\,m^{-1}}$ では最大と最小の位置差は $10\,\mathrm{nm}$ である。場を強めると周期は短くなるが、近接帯への非断熱遷移が無視できるかを別に確認する必要がある。
:::

### 三次元で戻る条件

一様な場で波数の軌道は $\vb k(t)=\vb k_0+q\vb E t/\hbar$ である。波数を逆格子で同一視するため、軌道が閉じる条件はある $T>0$ で $q\vb E T/\hbar=\sum_i n_i\vb b_i$、$n_i\in\Z$ となることである。実格子の基本並進と内積を取れば $\omega_i T=2\pi n_i$、$\omega_i=q\vb E\cdot\vb a_i/\hbar$ となる。従って非零の $\omega_i$ の比がすべて有理数なら共通周期があり、無理数を含む場合には波数軌道に厳密な周期がない。それでも速度の一成分だけを見ると周期的に見えることがある。

分散を Fourier 展開して $E_n(\vb k)=\sum_{\vb R}E_{\vb R}e^{i\vb k\cdot\vb R}$ と書く。滑らかな帯を考え、項別微分と積分ができるとする。群速度には $\sum_{\vb R}(i/\hbar)\vb R E_{\vb R}e^{i\vb k_0\cdot\vb R}e^{iq\vb E\cdot\vb R\,t/\hbar}$ が現れる。$\vb E\cdot\vb R\ne0$ の項は時間積分で振動を与える。$\vb E\cdot\vb R=0$ の項は速度に定数を加え、その $\vb R$ は場に垂直である。仕事の式によって束縛される変位は $\vb E$ 方向の成分であり、直交成分にはドリフトがあり得る。

:::example[二つの異なる振動数]
正方格子の分散を $-2t(\cos k_xa+\cos k_ya)$ とし、$|E_y|=2|E_x|$ とする。$y$ の振動周期は $x$ の半分なので、波数軌道の共通周期は $x$ の周期である。一方、$E_y=0$ なら $k_y$ は変わらず、例えば $k_ya=\pi/2$ の波束は横へ速さ $2ta/\hbar$ で進み続ける。$x$ は有界に振動する。
:::

### Wannier–Stark ラダーの厳密な一帯解

一次元の直交した Wannier 軌道 $\ket{j}$ に対し、隣接ホッピングを $-t_0$、同じ胞のエネルギーを $\epsilon_0$ とする。中心位置 $ja$ の軌道に対して、位置演算子の対角項だけを残す近似をとる。電場 $E_{\rm field}$ による力を $F=q E_{\rm field}$ とすれば胞ごとのポテンシャルは $-Faj$ である。無限鎖上の固有方程式は

$$
 -t_0(c_{j+1}+c_{j-1})+(\epsilon_0-Faj)c_j=\mathcal E c_j.
$$

$\widetilde c(k)=\sum_j c_j e^{-ikja}$ と定めると $ja$ は $i\partial_k$ に、ホッピングは $-2t_0\cos ka$ に写る。そこで

$$
 -iF\partial_k\widetilde c
 =(\mathcal E-\epsilon_0+2t_0\cos ka)\widetilde c
$$

を解くと $\widetilde c(k)\propto\exp\{i[(\mathcal E-\epsilon_0)k+(2t_0/a)\sin ka]/F\}$ になる。$k\mapsto k+2\pi/a$ で同じ値となるためには $(\mathcal E-\epsilon_0)/(Fa)$ が整数でなければならない。中心胞を $j_0$ として符号を選べば、エネルギー列は

$$
 \mathcal E_{j_0}=\epsilon_0-Fa j_0,\qquad j_0\in\Z.
 \label{eq-ss-stark-ladder}
$$

この等間隔列を :term[Wannier–Stark ラダー]{en="Wannier–Stark ladder" yomi="わにえしゅたるくらだー"} と呼ぶ。固有関数の $k$ 表示は $e^{-ikaj_0}e^{iz\sin ka}$、$z=2t_0/(Fa)$ である。実空間の係数はこの Fourier 係数であり、$e^{iz\sin x}=\sum_{s\in\Z}J_s(z)e^{isx}$ の展開で定義される Bessel 関数 $J_s$ を用いると $c_{j_0+\ell}=J_{-\ell}(z)$ となる。直交性は、すべての状態に共通する位相 $e^{iz\sin ka}$ がユニタリ変換であることから従う。平均位置と位置分散も $i\partial_k$ を使って求まる。$i\partial_k$ の期待値は $aj_0-za\cos ka$ のゾーン平均 $aj_0$、位置分散は $(za)^2/2=2t_0^2/F^2$ である。

同じ Hamiltonian による時間発展で、時刻 $T_B=2\pi\hbar/(|F|a)$ のあと、各固有状態の位相は $e^{-i\epsilon_0T_B/\hbar}$ にそろう。従って波束全体がこの共通位相を除いて戻る。半古典論の周期は、量子論では一帯のラダー間隔を $\hbar$ で割った角振動数に対応する。場で別の帯に遷移できる場合や、端のある有限鎖ではこの厳密な式の前提が失われる。

:::example[場を強めるとラダー状態が狭くなる]
$t_0=0.050\,\mathrm{eV}$、一胞でのポテンシャル差 $|Fa|=0.025\,\mathrm{eV}$ なら $|z|=4$ である。固有状態の位置の標準偏差は $|z|a/\sqrt2=2\sqrt2a$ になる。電場を二倍にするとラダー間隔は二倍、標準偏差と再帰時間は半分になる。ここでの幅は一つの定常状態の幅であり、前節で求めた移動波束の最大と最小の位置差とは異なる。
:::

### 近い帯がある場合の時間尺度

ゾーン境界で二つの自由波が混ざる場合、二波行列から共通エネルギーを引いて $h(t)=\begin{pmatrix}\alpha t&U\\U&-\alpha t\end{pmatrix}$ と近似する。ここで $U>0$ は位相を選んで実数にできる結合であり、$t=0$ が交差の中心である。二本の自由分散の差の時間微分が $2\alpha$ である。固有ベクトルの混合角を $\theta$ とすれば、$\cos\theta=\alpha t/\sqrt{\alpha^2t^2+U^2}$、$\sin\theta=U/\sqrt{\alpha^2t^2+U^2}$ である。角度を微分すると $|\dot\theta|=|\alpha|U/(\alpha^2t^2+U^2)$、分離は $2\sqrt{\alpha^2t^2+U^2}$ になる。この二状態の瞬時基底の非対角な時間微分は $|\dot\theta|/2$ なので、境界中心での無次元の断熱指標は $\hbar|\alpha|/(4U^2)$ である。

式 :eqref[eq-ss-nearly-free-gap] の自由波を使えば $\alpha=(\hbar^2G/(2m))\dot k$、$\hbar\dot k=qE_{\rm field}$ である。従って指標は場の強さに比例し、ギャップ $2U$ の二乗に反比例する。また一周期を緩和時間 $\tau$ より短くするには、より強い場が必要になる。実験で Bloch 振動を分解するには、この二つの条件が両立する領域を探す。上の局所的な指標は正確な遷移確率ではない。開始と終了時刻が有限の場合には、境界から離れた区間での位相も影響する。

### 単一交差の Landau–Zener 遷移

二準位近似を交差の前後へ理想的に延長すると、瞬時の固有帯を乗り換える確率を求められる。結合を $U>0$、$\alpha>0$ とし、$h(t)=\alpha t\sigma_z+U\sigma_x$ が $-\infty<t<\infty$ で成立するとする。反対符号の掃引は準位を交換すればよい。固定基底の第一成分 $a$、第二成分 $b$ に対し、$b=(i\hbar\dot a-\alpha t a)/U$ を第二の Schrödinger 方程式へ代入すると

$$
 \ddot a+\left(\frac{\alpha^2t^2+U^2}{\hbar^2}+\frac{i\alpha}{\hbar}\right)a=0 .
$$

$z=e^{i\pi/4}\sqrt{2\alpha/\hbar}\,t$、$\nu=-iU^2/(2\hbar\alpha)$ と置くと $a''+(\nu+1/2-z^2/4)a=0$ になる。ここでプライムは $z$ 微分を表す。第一成分だけが出射する解を得るため、$\arg z=\pi/4$ に沿って $D_\nu(z)\sim z^\nu e^{-z^2/4}$ となる解を取る。この位相の先頭は $e^{-i\alpha t^2/(2\hbar)}$ なので、$b$ を復元したときに $\alpha t$ の項が打ち消され、$b$ は $1/t$ の桁になる。

:::lemma[交差で同じ固定成分に戻る振幅]
$D_\nu(z)$ を負の時刻側で $D_\nu(-z)$ と $D_{-\nu-1}(iz)$ の一次結合で表すと、前者の係数は $e^{-i\pi\nu}$ である。
:::

:::proof
$\Re\nu<0$ では
$$
 D_\nu(z)=\frac{e^{-z^2/4}}{\Gamma(-\nu)}
 \int_0^\infty e^{-zs-s^2/2}s^{-\nu-1}\dd s
$$
を使う。積分を微分して $s$ に関する全微分を積分すれば上の微分方程式を満たす。$z=0$ で $s^2/2$ に変数を替え、Gamma 関数の倍角公式を使うと
$$
 D_\nu(0)=\frac{2^{\nu/2}\sqrt\pi}{\Gamma((1-\nu)/2)},\qquad
 D_\nu'(0)=-\frac{2^{(\nu+1)/2}\sqrt\pi}{\Gamma(-\nu/2)}
$$
である。$s$ の積分を零の近くで展開すると、$\arg z=\pm\pi/4$ では先頭項が $z^\nu e^{-z^2/4}$ になる。その際、$e^{-s^2/2}=1+\order{s^2}$ であり、$\Re z\ge |z|/\sqrt2$ が残りの積分を抑える。これらの式は $\Re\nu=0$ へ、部分積分から得る漸化式 $D_\nu=zD_{\nu-1}-(\nu-1)D_{\nu-2}$ によって延長できる。

負の時刻では $-z$ は偏角 $\pi/4$、$iz$ は偏角 $-\pi/4$ の半直線にある。第一の解 $D_\nu(-z)$ は第一成分が有限、第二の解 $D_{-\nu-1}(iz)$ は逆向きの動的位相を持ち第一成分が $1/|t|$ の桁となる。従って二解は独立である。$Y=D_\nu$、$W=D_{-\nu-1}$ と書くと、$z=0$ での値と微分の二元方程式から、第一の解の係数は
$$
 \frac{Y(0)iW'(0)-Y'(0)W(0)}{Y(0)iW'(0)+Y'(0)W(0)}
$$
となる。Gamma 関数の反射公式で $Y(0)iW'(0)=-i\cos(\pi\nu/2)$、$Y'(0)W(0)=\sin(\pi\nu/2)$ となるので、この比は $e^{-i\pi\nu}$ である。
:::

入射側 $D_\nu(-z)$ と出射側 $D_\nu(z)$ は、同じ偏角 $\pi/4$ で第一成分の漸近規格化が等しい。従って出射が第一成分だけの解を逆にたどったとき、入射で第一成分にいる確率は $|e^{-i\pi\nu}|^2$ である。時間発展のユニタリ性により、順方向で第一成分から第一成分へ行く確率も同じである。第一成分は負の時刻には下の固有帯、正の時刻には上の固有帯に漸近する。この :term[Landau–Zener 遷移]{en="Landau–Zener transition" yomi="らんだうつぇなーせんい"} の確率は

$$
 P_{\rm jump}=\exp\!\left(-\frac{\pi U^2}{\hbar|\alpha|}\right).
 \label{eq-ss-landau-zener}
$$

実際の結晶でこの式を使うには、交差から離れた位置で二帯の混成が小さく、遷移の起こる間は自由分散の差を時間に線形とみなせることが必要である。繰り返して交差を通ると二つの遷移経路の位相が干渉し得るので、相干性が保たれている場合は確率だけを各回独立に掛けられない。

:::example[ギャップと電場の効果]
$\pi U^2/(\hbar|\alpha|)=2$ の掃引では $P_{\rm jump}=e^{-2}\simeq0.135$ である。電場を二倍にして $|\alpha|$ を二倍にすれば $e^{-1}\simeq0.368$ になる。同じ電場のまま $U$ を二倍にすれば元の例の指数は $8$ となり、確率は約 $3.35\times10^{-4}$ に減る。
:::

## 磁場中の半古典運動とサイクロトロン共鳴

### 断面積のエネルギー微分

$\vb E=0$、$\vb B=B\hat{\vb z}$、$B>0$ なら $E_n$ と $k_z$ は保存する。従って運動は等エネルギー面と平面 $k_z=\mathrm{const}$ の交わりに沿う。閉じた曲線の弧長を $\dd l$ とすると $\dd l/\dd t=|q|B|\grad_\perp E|/\hbar^2$ である。曲線をエネルギーで微小に変えると法線変位の絶対値は $\dd E/|\grad_\perp E|$ なので、断面積の変化を積分して
$$
 T_c=\frac{\hbar^2}{|q|B}\left|\pdv{A}{E}\right|,\qquad
 m_c=\frac{\hbar^2}{2\pi}\left|\pdv{A}{E}\right|,\qquad
 \omega_c=\frac{|q|B}{m_c}
 \label{eq-ss-cyclotron-mass}
$$
を得る。交流電場が軌道の回転と同じ周波数および向きを持てば繰り返し仕事をし、吸収が増える。これが :term[サイクロトロン共鳴]{en="cyclotron resonance" yomi="さいくろとろんきょうめい"} である。有限 $\tau$ の円軌道上では速度を円偏極成分に分けると、応答の分母が $1/\tau-i(\omega\mp\omega_c)$ となる。共鳴幅の桁は $1/\tau$ である。

:::example[楕円体の断面]
$E=\sum_i\hbar^2k_i^2/(2m_i)$ では $k_z$ 固定の楕円の面積は $2\pi\sqrt{m_xm_y}(E-\hbar^2 k_z^2/(2m_z))/\hbar^2$ である。従って $m_c=\sqrt{m_xm_y}$ で、$m_z$ はこの向きの共鳴周波数に出ない。
:::

:::example[運動の幾何]
$\hbar \dot{\vb k}=qB\dot{\vb r}\times\hat{\vb z}$ を平面内で積分すれば $\vb r_\perp(t)-\vb r_\perp(0)=\hbar\hat{\vb z}\times[\vb k_\perp(t)-\vb k_\perp(0)]/(qB)$ となる。円形 $k$ 軌道の半径 $k_\perp$ は実空間で $r_c=\hbar k_\perp/(|q|B)$ に対応する。
:::

### 軌道を一周したときの縦方向の移動

軌道を囲む領域がエネルギーとともに広がる場合を $s_A=+1$、狭まる場合を $s_A=-1$ とする。断面の境界ではエネルギーの勾配の外向き成分が $s_A|\grad_\perp E_n|$ である。$k_\parallel$ を微小にずらしてエネルギーを固定すると、境界の外向き変位は

$$
 \delta l_\perp=-\frac{\partial E_n/\partial k_\parallel}
 {s_A|\grad_\perp E_n|}\,\delta k_\parallel
$$

となる。これは $E_n(\vb k)$ を境界で一次に展開した式である。周の弧長で積分すると $(\partial A/\partial k_\parallel)_E$ が得られる。一方、一周期の縦の変位は、$v_\parallel=\partial_{k_\parallel}E_n/\hbar$ を $\dd t=\hbar^2\dd l/(|q|B|\grad_\perp E_n|)$ で積分した値になる。二つを比較すると

$$
 \Delta r_\parallel=-s_A\frac{\hbar}{|q|B}
  \left(\pdv{A}{k_\parallel}\right)_E,\qquad
 \overline v_\parallel
 =-\frac{1}{\hbar}\frac{(\partial A/\partial k_\parallel)_E}
  {(\partial A/\partial E)_{k_\parallel}}.
 \label{eq-ss-orbit-longitudinal}
$$

断面積が極値の軌道では縦の変位も零になる。この結果は閉曲線の近くで形と囲む領域を連続に選べる場合に成り立つ。鞍点を通る軌道のように勾配が零になる場合には周期の式自体の発散を点検する必要がある。

:::example[円筒と球の違い]
二次元的な分散が $k_z$ に依存しなければ $A(E,k_z)$ は一定で、すべての断面の縦速度は零になる。一方、自由電子の球では $A(E,k_z)=2\pi mE/\hbar^2-\pi k_z^2$ である。式 :eqref[eq-ss-orbit-longitudinal] による周期平均は $\hbar k_z/m$ で、実際の一定な $v_z$ と一致する。球の赤道 $k_z=0$ だけは縦移動がない。
:::

### 楕円体の傾いた切断

主軸質量が正の放物線の等エネルギー面を楕円体とする。半径 $K_i=\sqrt{2m_i E}/\hbar$ を使い、単位球を $\vb k=D\vb p$、$D=\diag(K_1,K_2,K_3)$ で移す。向き $\vb n$ に垂直な断面は $\vb n\cdot\vb k=k_\parallel$ であり、球では $\vb p\cdot D\vb n=k_\parallel$ になる。$r_n=|D\vb n|$ とすれば、球の断面は面積 $\pi(1-k_\parallel^2/r_n^2)$ の円である。面積の変換係数は $\det D/r_n$ である。実際、断面の厚みは球で $\dd k_\parallel/r_n$、楕円体で $\dd k_\parallel$ となり、三次元の体積は $\det D$ 倍になるためである。従って中心断面の面積とサイクロトロン質量は

$$
 A_{\rm cent}=\frac{2\pi E}{\hbar^2}
  \sqrt{\frac{\det m^*}{\vb n^{\mathsf T}m^*\vb n}},
 \qquad
 m_c(\vb n)=\sqrt{\frac{\det m^*}{\vb n^{\mathsf T}m^*\vb n}}.
 \label{eq-ss-ellipsoid-cyclotron}
$$

同じ $\vb n$ で $k_\parallel\ne0$ の断面は、中心断面から $E$ に依存しない項を引いたものになる。これは $r_n^2=2E\,\vb n^{\mathsf T}m^*\vb n/\hbar^2$ により分かる。従って $m_c$ はこの放物線模型では切断の高さにも依存しない。全状態数から得る有効質量と、特定の断面で得るサイクロトロン質量は異なる平均である。式 :eqref[eq-ss-band-edge-dos] の三次元の状態密度は $\sqrt{\det m^*}$ に比例するので、等方模型に換算した質量は $(\det m^*)^{1/3}$ である。

:::example[結晶軸から傾けた共鳴]
$m_x=m_y=m$、$m_z=4m$ の楕円体で $\vb n=(1,0,1)/\sqrt2$ とする。このとき $\det m^*=4m^3$、$\vb n^{\mathsf T}m^*\vb n=5m/2$ であり、$m_c=\sqrt{8/5}\,m$ である。$\vb n=\hat{\vb z}$ では $m_c=m$、$\vb n=\hat{\vb x}$ では $m_c=2m$ なので、傾けた結果はその間にある。
:::

### 円偏光に対する伝導度

放物線帯で $\vb B=B\hat{\vb z}$ とし、有限の緩和時間を入れた平均速度に $m\dot{\vb v}=q(\vb E+\vb v\times\vb B)-m\vb v/\tau$ を使う。$v_\pm=v_x\pm i v_y$、$E_\pm=E_x\pm iE_y$ と置くと $(\vb v\times\hat{\vb z})_\pm=\mp i v_\pm$ である。周波数因子を $e^{-i\omega t}$ として

$$
 \sigma_\pm(\omega)=\frac{n_cq^2/m}{\tau^{-1}-i(\omega\mp qB/m)}
 \label{eq-ss-cyclotron-optical}
$$

を得る。正周波数の吸収は $qB>0$ では $+$、$qB<0$ では $-$ の成分に現れる。ここでの偏光の符号は $v_x\pm i v_y$ の定義に従う。偏光を指定せずに共鳴周波数だけを測る場合、分かるのは質量の絶対値である。

分母を実数化すると $\Re\sigma_\pm=(n_cq^2\tau/m)/[1+(\omega\mp qB/m)^2\tau^2]$ である。半値幅は $1/\tau$、全幅は $2/\tau$ になる。非放物線帯では、閉軌道上の速度 $\vb v(t)$ が正弦波とは限らない。周期関数を Fourier 展開すると $\vb v(t)=\sum_p\vb v_p e^{ip\omega_c t}$ であるから、仕事 $q\vb v\cdot\vb E$ を一周で平均すると、場の周波数が $p\omega_c$ と一致する項が残る。各係数の有無は軌道の対称性に依存する。

:::example[周波数から質量を逆算する]
一価キャリアで $B=2.0\,\mathrm T$ において周期の周波数 $56\,\mathrm{GHz}$ が共鳴位置であるとする。角振動数 $2\pi\cdot56\times10^9\,\mathrm{s^{-1}}$ を使うと $m_c=eB/\omega_c\simeq9.1\times10^{-31}\,\mathrm{kg}$ となる。Hz で示された値をそのまま $\omega_c$ に代入すると、質量を $2\pi$ 倍大きく見積もる。
:::

### 開軌道と実空間の横移動

半古典方程式は、常に閉じた曲線を与えるわけではない。例えば二次元正方格子で $E=-2t_x\cos k_xa-2t_y\cos k_ya$、$t_x>t_y>0$、$E=0$ とする。等エネルギー線は $k_xa=\pm\arccos[-(t_y/t_x)\cos k_ya]$ という二本である。それぞれの枝で $\sin k_xa$ は零にならず、$\dot k_y$ の符号は途中で変わらない。$k_y$ が $2\pi/a$ 進むごとに $k_x$ は元に戻る。従ってゾーンを同一視すれば同じ波数に戻るが、逆格子を引かず追跡すると軌道は $k_y$ 方向に延びる。この意味で開軌道である。

変位の式 $\Delta\vb r_\perp=\hbar\hat{\vb z}\times\Delta\vb k_\perp/(qB)$ は、逆格子を引かない $\Delta\vb k$ を使う。上の例で $\Delta k_y=2\pi/a$ なら $\Delta r_x=-2\pi\hbar/(qBa)$ である。開軌道を閉軌道に見立てて面積量子化をしてはならない。開軌道に沿った実空間の移動は磁気抵抗の異方性にも現れ得る。

### 共鳴の高調波を軌道から求める

閉軌道上の時間座標を $s$ とする。$\omega_c=2\pi/T_c$、一定緩和時間 $\tau$ を仮定し、電場 $\vb E_\omega e^{-i\omega t}$ を一次で入れる。平衡分布からの変化を $\delta f(s)e^{-i\omega t}$ と書くと、上の弱場 Hall の式と同じ運動項を用いて

$$
 (\tau^{-1}-i\omega+\partial_s)\delta f
 =q w\,\vb E_\omega\cdot\vb v(s)
$$

を得る。今回は $\partial_s$ を小さいとは仮定しない。$v_i(s)=\sum_{p\in\Z}v_{ip}e^{ip\omega_c s}$ と展開すると各係数が独立に決まり、電流の一周平均に現れる項は

$$
 q^2 w\sum_{p\in\Z}
 \frac{v_{i,-p}v_{jp}}{\tau^{-1}-i(\omega-p\omega_c)}E_{\omega,j}.
$$

軌道上の状態数の測度は $\dd l/|\nabla_\perp E|$ に比例し、これは $\dd s$ に比例する。従って Fourier モードを一周平均で直交させてよい。同じ方向の吸収では $v_{i,-p}v_{ip}=|v_{ip}|^2$ なので各共鳴の重みは非負である。軌道の速度が半周期で反転する場合、積分を二つの半周期に分けると偶数 $p$ の係数は零になる。円で正弦波しかない場合はさらに $p=\pm1$ だけとなり、円偏光の式 :eqref[eq-ss-cyclotron-optical] と整合する。

:::example[第三高調波の吸収]
$v_x(s)=v_1\cos\omega_c s+(v_1/5)\cos 3\omega_c s$ で、同じ軌道と $\tau$ を比較する。共鳴成分の重みはそれぞれ $v_1^2/4$ と $v_1^2/100$ なので、第三高調波の共鳴成分は基本波の $1/25$ である。ピークを区別するには幅 $\tau^{-1}$ が間隔 $2\omega_c$ より十分小さいことも必要になる。
:::

## Landau 準位

### 一様磁場による横方向の量子化

自由電子で $\vb A=(0,Bx,0)$ とする。最小結合 $(\vb p-q\vb A)^2/(2m)$ では $p_y,p_z$ は可換な保存演算子である。$e^{i(k_yy+k_zz)}$ を代入すると
$$
 H=\frac{p_x^2}{2m}+\frac12m\omega_c^2(x-x_0)^2+
 \frac{\hbar^2k_z^2}{2m},\qquad x_0=\frac{\hbar k_y}{qB}.
$$
従って :term[Landau 準位]{en="Landau level" yomi="らんだうじゅんい"} は
$$
 E_{j,k_z}=\hbar\omega_c(j+1/2)+\frac{\hbar^2 k_z^2}{2m},\qquad j=0,1,\ldots
 \label{eq-ss-landau}
$$
である。$y$ に周期 $L_y$ を課すと許される中心の間隔は $2\pi\hbar/(|q|BL_y)$ になる。$x$ 方向の幅 $L_x$ に中心があるものを数えれば、一スピンで縮退数は $L_xL_y |q|B/(2\pi\hbar)$ である。端から磁気長 $\ell_B=\sqrt{\hbar/(|q|B)}$ 程度の中心は境界により変形するが、面積当たりの熱力学極限では寄与は消える。スピンを入れる場合は電子の磁気モーメントによる Zeeman エネルギー $\pm g\mu_BB/2$ を加える。

一般の閉軌道の量子化に移る前に、自由電子の円では $A=2\pi mE_\perp/\hbar^2$ である。式 :eqref[eq-ss-landau] を代入すると $A=2\pi eB(j+1/2)/\hbar$ が出る。円でない軌道に対しても面積が量子化される理由は、横の二つの波数成分を共役な座標へ写すことで分かる。

:::example[磁場の単位と縮退密度]
一スピンでは面積当たり一準位に $eB/h$ 個の状態がある。$B=1.0\,\mathrm T$ なら約 $2.42\times10^{14}\,\mathrm{m^{-2}}$ である。$\mu_B B\simeq5.79\times10^{-5}\,\mathrm{eV}$、自由電子の準位間隔はその二倍程度である。
:::

:::example[異方的な放物線]
$x,y$ の質量が $m_x,m_y$ なら振動子の質量は $m_x$、ばねは $q^2B^2/m_y$ になる。従って準位間隔は $\hbar |q|B/\sqrt{m_xm_y}$ で、式 :eqref[eq-ss-cyclotron-mass] と一致する。
:::

### 運動量の交換関係と案内中心

運動学的運動量 $\Pi_i=p_i-qA_i$ の交換関係は、$[p_i,A_j]=-i\hbar\partial_iA_j$ によって $[\Pi_x,\Pi_y]=i\hbar qB$ となる。$s_B=\operatorname{sgn}(qB)$ と置いて

$$
 a_B=\frac{\ell_B}{\sqrt2\hbar}(\Pi_x+i s_B\Pi_y)
$$

を作ると、$[a_B,a_B^\dagger]=1$ である。積 $a_B^\dagger a_B$ では交差項が $i s_B[\Pi_x,\Pi_y]=- \hbar|qB|$ になる。従って横方向の Hamiltonian は $\hbar\omega_c(a_B^\dagger a_B+1/2)$ となり、特定のゲージを選ばずに式 :eqref[eq-ss-landau] を再現できる。

同じ準位内の自由度を表す :term[案内中心]{en="guiding center" yomi="あんないちゅうしん"} を

$$
 R_x=x+\Pi_y/(qB),\qquad R_y=y-\Pi_x/(qB)
$$

とする。半古典的な円の中心に対応する。例えば $[R_x,\Pi_x]=[x,\Pi_x]+[\Pi_y,\Pi_x]/(qB)=0$ であり、他の組も同様に可換である。従って $R_x,R_y$ は横の Hamiltonian と可換である。しかし中心同士の交換関係は零でなく、

$$
 [R_x,R_y]=-i\hbar/(qB)
$$

である。$[x,-\Pi_x]/(qB)$、$[\Pi_y,y]/(qB)$、$[\Pi_y,-\Pi_x]/(qB)^2$ の三項を足して得られる。$R_x$ を固有値で指定した状態が $R_y$ まで一点に局在するわけではない。$\vb A=(0,Bx,0)$ では $R_x=p_y/(qB)$ となるので、固有値は前節の $x_0$ に一致する。

### 磁束と縮退数

長方形 $L_x\times L_y$ の向かい合う辺をつなげて境界のない系を考える。この面に正の磁場 $B$ を貫く場合、$\vb A=(0,Bx,0)$ は $x$ 方向に周期的でない。$x\mapsto x+L_x$ で $\vb A$ の変化は $\grad(BL_x y)$ なので、最小結合を保つには波動関数に

$$
 \psi(x+L_x,y)=e^{iqBL_x y/\hbar}\psi(x,y),\qquad
 \psi(x,y+L_y)=\psi(x,y)
$$

を課す。二つの移動を逆の順序で行って同じ値になるには $qB L_xL_y/(2\pi\hbar)$ が整数でなければならない。$y$ の運動量の間隔は $2\pi\hbar/L_y$ で、$x$ の境界条件はその添字を $qB L_xL_y/(2\pi\hbar)$ だけ移す。対応する中心は $L_x$ だけずれて同じ状態の周期的な像になる。従って独立な中心の数は添字を絶対値を取ったその整数で割った剰余類の個数である。面積を $S$ と書けば、一準位一スピンの縮退は

$$
 D_B=\frac{S}{2\pi\ell_B^2}=\frac{BS}{h/|q|}
 \label{eq-ss-landau-degeneracy}
$$

となる。$h/|q|$ は一状態に対応する磁束である。面積の周期境界を課さない場合にも、端の寄与を除く熱力学極限で単位面積当たりの結果は一致する。

### 三次元の状態密度と量子極限

一スピンで $E_j=\hbar\omega_c(j+1/2)$ と置く。$z$ 軸方向の状態数は単位長さ当たり $\dd k_z/(2\pi)$ である。二方向の縮退密度を掛けて $E$ へ変数変換すると、$j$ 番目の底からの状態密度の寄与は

$$
 g_j(E)=\frac{\Theta(E-E_j)}{2\pi^2\ell_B^2\hbar}
 \sqrt{\frac{m}{2(E-E_j)}}.
$$

$E>E_j$ では同じエネルギーに $k_z=\pm\sqrt{2m(E-E_j)}/\hbar$ の二根があり、各根の Jacobian が $1/|\partial E/\partial k_z|$ である。従って滑らかな三次元の自由電子の状態密度は、場の中では一次元の底が多数重なった形になる。場が弱く底の間隔より粗く平均すると、和を積分に置き換えられる。$E_j$ の間隔 $\hbar\omega_c$ で割り、$\int_0^E (E-u)^{-1/2}\dd u=2\sqrt E$ を使えば、一スピンの自由電子の $g(E)=m\sqrt{2mE}/(2\pi^2\hbar^3)$ を得る。

強い場で固定密度 $n$ の一スピンの電子が最低準位だけを占めると、残った方向の Fermi 波数は

$$
 n=\frac{1}{2\pi\ell_B^2}\frac{2k_{F,z}}{2\pi},
 \qquad k_{F,z}=2\pi^2 n\ell_B^2
$$

である。従って最低準位の底から測った $\hbar^2 k_{F,z}^2/(2m)$ は $B^{-2}$ に比例する。これが準位間隔より小さいことが、最低準位だけを使う量子極限の自己整合条件である。スピンが分離しない場合は縮退に二を掛ける。場を変えても化学ポテンシャルが固定だと仮定した計算とは、粒子数を固定するこの計算を区別する。

:::example[最低準位での縦方向の状態数]
一スピンの電子密度を $n=10^{22}\,\mathrm{m^{-3}}$、磁場を $10\,\mathrm T$ とする。$\ell_B^2\simeq6.58\times10^{-17}\,\mathrm{m^2}$ なので $k_{F,z}\simeq1.30\times10^7\,\mathrm{m^{-1}}$ である。裸の電子質量なら底からのエネルギーは約 $6.4\times10^{-6}\,\mathrm{eV}$ であり、準位間隔の約 $1.16\times10^{-3}\,\mathrm{eV}$ より小さい。この比較ではスピン分離を仮定した。
:::

### 閉軌道の面積量子化

結晶の帯を孤立した有効 Hamiltonian $E_n(\vb k)$ で表し、磁気長が格子間隔に比べて大きい場合を考える。磁場を入れる最低次の近似では $\hbar\vb k$ を、互いの積の順序に注意しながら運動学的運動量 $\vb\Pi$ へ置き換える。これは位置に緩やかに変わるベクトルポテンシャルを用いた一帯近似であり、一般の結晶で $\hbar\vb k$ が裸の電子の運動量になったことを意味しない。この近似では $[k_x,k_y]=i qB/\hbar$ となる。そこで

$$
 Q=k_x,\qquad P=\frac{\hbar^2}{qB}k_y
$$

と置けば $[Q,P]=i\hbar$ である。磁場に平行な波数を固定し、単純な閉曲線に二つの滑らかな転回点があり、一つの軌道から他の軌道へのトンネルを無視できるとする。式 :eqref[eq-qm-wkb-quantization] の往復の作用は $|\oint P\,\dd Q|=2\pi\hbar(j+1/2)$ である。左辺は $k$ 空間で囲む面積 $A$ の $\hbar^2/(|q|B)$ 倍なので

$$
 A(E_j,k_\parallel)=\frac{2\pi |q|B}{\hbar}(j+\gamma),
 \qquad\gamma=\frac12
 \label{eq-ss-orbit-quantization}
$$

が得られる。半整数のずれは二つの転回点をつないだ位相から来る。Berry 位相や軌道磁気モーメントによるずれが無視できない帯では、位相とエネルギーの補正を加える必要がある。両者の議論は [トポロジカル物質入門](/physics/solid-state/topological/) に譲る。自由電子の円ではこの式が低い $j$ まで厳密解と一致するが、一般の曲線では小さい場で $j$ が大きい領域の近似である。

:::example[分散が放物線でない閉軌道]
二次元で $E=\alpha k^4$、$\alpha>0$ とし、幾何学的な位相補正がない場合を考える。面積は $A=\pi\sqrt{E/\alpha}$ なので、式 :eqref[eq-ss-orbit-quantization] から $E_j=4\alpha(j+1/2)^2/\ell_B^4$ となる。差 $E_{j+1}-E_j=8\alpha(j+1)/\ell_B^4$ は一定でない。大きい $j$ では $A$ の差を $\partial_E A$ で割って求めた $\hbar\omega_c$ と一致する。これは面積とエネルギーとで等間隔の意味が異なる例である。
:::

### 案内中心の並進と一周する位相

磁場中の連続平面では、通常の並進にゲージの補正を組み合わせると Hamiltonian と可換になる。案内中心から $K_x=-qBR_y$、$K_y=qBR_x$ を作り、変位 $\vb d$ に対し $T(\vb d)=\exp[-i\vb d\cdot\vb K/\hbar]$ と置く。$[R_i,K_j]=i\hbar\delta_{ij}$ なので $T(\vb d)^\dagger\vb R T(\vb d)=\vb R+\vb d$ である。一方で $K$ は運動学的運動量と可換なので横方向のエネルギーを変えない。これが :term[磁気並進]{en="magnetic translation" yomi="じきへいしん"} である。

$[K_x,K_y]=-i\hbar qB$ を使うと、$-i\vb d\cdot\vb K/\hbar$ と $-i\vb d'\cdot\vb K/\hbar$ の交換子は定数 $iqB(\vb d\times\vb d')_z/\hbar$ となる。定数は二つの指数と可換なので、この場合の指数の交換にはその交換子の指数だけが残る。よって

$$
 T(\vb d)T(\vb d')=
 e^{iqB(\vb d\times\vb d')_z/\hbar}T(\vb d')T(\vb d).
$$

二方向の順序を入れ替えると、囲んだ磁束に応じた位相差が出る。面積 $L_xL_y$ の周期境界で両周期の並進が共通の固有状態を持てる条件は、この位相が $1$ であることになる。これは波動関数の境界条件から得た磁束量子化と同じ条件である。周期格子に弱い磁場を入れる場合も位相の交換関係が残り、格子並進を扱うときの胞の選び方に関わる。

:::example[半磁束と二倍の周期]
一原始胞の磁束が $h/(2|q|)$ の場合、二つの基本並進は位相 $-1$ を挟むので同時に固有値を指定できない。一方の並進を二回行えば位相は $+1$ に戻るため、一方向を二倍とした磁気単位胞を選べる。二倍胞の面積を貫く磁束は一磁束量子である。
:::

## de Haas–van Alphen 効果と Fermi 面の測定

### 準位が化学ポテンシャルを通る周期

閉じた軌道の断面積を $A(E,k_z)$ とする。化学ポテンシャル $\mu$ で準位条件を解くと $F/B=j+\gamma$、$F=\hbar A(\mu,k_z)/(2\pi e)$ である。隣接準位の通過は
$$
 \Delta(1/B)=\frac{2\pi e}{\hbar A(\mu,k_z)}
 \label{eq-ss-onsager}
$$
を満たす。固定化学ポテンシャルで一準位の大きな縮退が占有に加わると、熱力学ポテンシャルの磁場依存に折れが生じ、$M=-\partial_B\Omega$ に振動が現れる。これを :term[de Haas–van Alphen 効果]{en="de Haas–van Alphen effect" yomi="ではーすふぁんあるふぇんこうか"} という。

三次元では $k_z$ を積分する。振動成分の位相は $2\pi F(k_z)/B$ であり、$F'$ が零でない区間は部分積分により $B$ の次数で小さくなる。内部の非退化停留点では $F=F_0+F''(k_z-k_0)^2/2+\cdots$ と展開し、幅 $\sqrt{B/|F''|}$ の範囲が主要な寄与となる。従って通常は極値断面が観測の周波数を決める。温度による丸まりはエネルギー $E$ ごとの零温応答を $-f'(E)$ で平均する効果である。$F(E)$ を一次化すると位相の傾きは $2\pi/(\hbar\omega_c)$ であり、幅 $k_BT$ の平均によって $k_BT\gg\hbar\omega_c$ では振動が見えにくくなる。散乱による準位の幅が間隔より大きい場合も同様である。

:::example[周期を断面積に直す]
$F=100\,\mathrm T$ なら $\Delta(1/B)=0.010\,\mathrm{T^{-1}}$ である。円断面なら $k_F=\sqrt{2eF/\hbar}\simeq5.5\times10^8\,\mathrm{m^{-1}}$ になる。
:::

:::example[磁場を回す]
球状の Fermi 面なら極値断面積は方向によらない。半径 $k_x,k_y,k_z$ の楕円体では $z$ 軸方向の場に対し $\pi k_x k_y$、$x$ 軸方向の場に対し $\pi k_y k_z$ である。測定周波数の比 $F_x/F_z$ は $k_z/k_x$ になる。
:::

### 二次元自由電子の磁気振動を和から求める

温度零、固定化学ポテンシャル $\mu$ の二次元自由電子を考える。スピンは一種類とし、$b=\hbar\omega_c=\hbar eB/m$、一準位の単位面積当たりの状態数を $g_B=eB/h$ とする。最低準位より $\mu$ が大きいとき、占有する準位の添字は $0$ から $N=\lfloor x\rfloor$、$x=\mu/b-1/2$ までである。等号の準位を足しても、以下の熱力学ポテンシャルには影響しない。面積を $S$ とすると、粒子一個当たり $E-\mu$ を足して

$$
 \frac{\Omega_0}{S}
 =-g_B b\sum_{j=0}^{N}(x-j)
 =-\frac{g_B b}{2}\bigl(x^2+x+f(1-f)\bigr),
 \qquad f=x-N
$$

を得る。最後は有限和 $(N+1)x-N(N+1)/2$ に $N=x-f$ を代入した結果である。$x^2+x=(\mu/b)^2-1/4$ なので、$B$ に対して急に変化する部分は $f$ のみに含まれる。$B_2(f)=f^2-f+1/6$ を周期関数に延長する。平均が零であり、$p\ne0$ の Fourier 係数は

$$
 \int_0^1 B_2(f)e^{-2\pi i p f}\dd f
 =\frac1{2\pi^2 p^2}
$$

である。一回目の部分積分では $B_2(0)=B_2(1)$ により境界項がなく、二回目では $B_2'(1)-B_2'(0)=2$ が残る。従って振動部分は

$$
 \frac{\Omega_{\rm osc}}{S}
 =\frac{g_B b}{2\pi^2}
 \sum_{p=1}^{\infty}\frac{1}{p^2}
 \cos\!\left[2\pi p(F/B-1/2)\right],\qquad
 F=\frac{m\mu}{e\hbar}
 \label{eq-ss-oscillation-omega}
$$

となる。円断面の $A(\mu)=2\pi m\mu/\hbar^2$ を用いれば $F$ は式 :eqref[eq-ss-onsager] の値と同じである。平坦な二次元模型で鋸歯状になるのは、その磁化が占有準位の個数の変化を直接反映するためである。$M/S=-\partial_B(\Omega/S)$ を計算すると、位相を微分した項は $B$ が小さい範囲で $1/p$ の振幅を持つ。これは $B$ とともに緩やかに変わる前因子の微分より大きい。丸まりを入れない零温度極限では、磁化が跳ぶ点の微分を片側極限で考える。

:::example[準位が横切る磁場を並べる]
$F=100\,\mathrm T$ とする。自由電子では $F/B=j+1/2$ だから、$j=9$ と $10$ の横切りは $B=200/19\,\mathrm T$ と $200/21\,\mathrm T$ である。磁場の差は約 $1.00\,\mathrm T$ だが、隣の組では差が変わる。逆数の差は $21/200-19/200=1/100\,\mathrm{T^{-1}}$ であり、こちらが一定である。
:::

### 三次元で極値断面が測られる理由

三次元の放物線帯では、$k_\parallel$ を固定すると二次元の結果を $k_\parallel$ ごとの化学ポテンシャルで使える。より一般の閉軌道では、式 :eqref[eq-ss-orbit-quantization] を化学ポテンシャルの近くで一次に展開すると、準位の局所的な間隔は

$$
 \frac{2\pi eB}{\hbar|\partial_E A|}
 =\frac{\hbar e B}{m_c}
$$

である。以下ではこの間隔と $m_c$ のエネルギー依存が一準位幅で緩やかであり、化学ポテンシャルの近くに多数の準位がある場合を考える。第 $p$ 高調波の $k_\parallel$ 積分には、緩やかな振幅に掛かった $e^{2\pi i pF(k_\parallel)/B}$ が現れる。

まず $F'$ が零にならない閉区間を考える。指数関数を微分すると $(2\pi i p/B)F'e^{2\pi i pF/B}$ なので、振幅を $F'$ で割って部分積分すれば積分は $B/p$ 程度に抑えられる。次に $F'(k_0)=0$、$F''(k_0)\ne0$ の点の周囲では、$F=F_0+F''(k_0)(k_\parallel-k_0)^2/2+\cdots$ と展開する。幅 $\sqrt{B/(p|F''|)}$ の周囲で位相がほぼ同じであり、振幅はそこでの値で置き換えられる。積分をこの尺度で拡大すると、先頭項には

$$
 \int_{-\infty}^{\infty}
  e^{i\pi pF''(k_0)u^2/B}\dd u
 =\sqrt{\frac{B}{p|F''(k_0)|}}\,
 e^{i\,\operatorname{sgn}(F''(k_0))\pi/4}
$$

が掛かる。収束因子 $e^{-\varepsilon u^2}$ を入れた Gaussian 積分を行い、$\varepsilon\downarrow0$ とすれば右辺が得られる。非停留区間の $B/p$ より停留点の $\sqrt{B/p}$ の方が弱磁場で大きい。これが極値の周波数が主に残る理由である。最大断面では位相補正が $-\pi/4$、最小断面では $+\pi/4$ となる。孤立した二次の停留点でない場合は、この位相と振幅則をそのまま使えない。

放物線帯の例では一断面当たりの式 :eqref[eq-ss-oscillation-omega] の振幅は $B^2/p^2$ に比例する。三次元の単位体積では $\dd k_\parallel/(2\pi)$ を掛けて積分するので、極値断面当たりの振幅は $B^{5/2}/p^{5/2}$ となる。磁化では位相微分の因子 $pF/B^2$ が掛かり、先頭は $B^{1/2}/p^{3/2}$ である。一方、完全な円筒では $F''=0$ で全断面が同位相になる。係数に円筒の高さが掛かり、同じ三次元の体積を持つ場合でも $B^{1/2}$ の因子は現れない。

:::example[球での停留位相の符号]
球状ポケットでは $A(\mu,k_z)=\pi(k_F^2-k_z^2)$ なので $F''=-\hbar/e$ である。赤道は最大断面であり、熱力学ポテンシャルの余弦には $-\pi/4$ が加わる。$k_z=\pm k_F$ では面積が零なので、端点の位相から有限の振動周波数は生じない。
:::

### 温度、散乱、スピンの振幅因子

有限温度の一準位の熱力学ポテンシャルは $-k_BT\log(1+e^{(\mu-E)/(k_BT)})$ である。零温度ポテンシャルを化学ポテンシャル $E'$ の関数として $-(E'-E)\Theta(E'-E)$ と書く。これを

$$
 w_T(E'-\mu)=\frac1{4k_BT\cosh^2((E'-\mu)/(2k_BT))}
$$

で平均すると有限温度の式が再現される。$w_T(E'-\mu)=-\partial_{E'}(1+e^{(E'-\mu)/(k_BT)})^{-1}$ を使って $E'\ge E$ の一次式を部分積分し、残る Fermi 分布を積分すればよい。この関係を全準位に適用できる。$k_BT$ の範囲で振幅を一定とみなし、周波数だけ $F(E')\simeq F(\mu)+F'(\mu)(E'-\mu)$ と近似する。すると振動の平均には

$$
 \int_{-\infty}^{\infty}
 \frac{e^{i y u}}{4\cosh^2(u/2)}\dd u,\qquad
 y=\frac{2\pi p k_BT F'(\mu)}{B}
$$

が掛かる。$y>0$ について、この積分を高さ $2\pi i$ の長方形に沿って評価する。縦辺は実部を無限に離すと消え、上辺は下辺の $-e^{-2\pi y}$ 倍になる。内部の極 $u=i\pi$ で $4\cosh^2(u/2)\simeq-(u-i\pi)^2$ なので、留数は $-iy e^{-\pi y}$ である。従って積分は $2\pi y e^{-\pi y}/(1-e^{-2\pi y})=\pi y/\sinh(\pi y)$ である。重みが偶関数なので負の $y$ でも同じ式になる。式 :eqref[eq-ss-cyclotron-mass] から $|F'(\mu)|=m_c/(e\hbar)$ なので、

$$
 R_T=\frac{X}{\sinh X},\qquad
 X=\frac{2\pi^2 p k_BT}{\hbar\omega_c}
 \label{eq-ss-lk-thermal}
$$

を得る。$X\ll1$ では補正が小さく、$X\gg1$ では $2Xe^{-X}$ となる。温度を変えた振幅の比から $m_c$ を推定できるのは、$X$ の温度依存に質量だけが残るためである。$p=2$ の高調波は $p=1$ より温度に敏感になる。

準位の線幅を半値半幅 $\Gamma$ の Lorentz 関数 $\Gamma/[\pi(\Gamma^2+x^2)]$ で模型化する。同じエネルギー平均をその分布で行うと、角周波数 $2\pi p/(\hbar\omega_c)$ に対して $e^{-2\pi p\Gamma/(\hbar\omega_c)}$ が掛かる。実際、正の周波数で上半平面の極 $i\Gamma$ の留数を取ると、指数は $e^{-\Gamma (2\pi p/\hbar\omega_c)}$ になる。振幅の時間減衰を $e^{-t/(2\tau_q)}$ としたとき $\Gamma=\hbar/(2\tau_q)$ なので、:term[Dingle 因子]{en="Dingle factor" yomi="でぃんぐるいんし"} は

$$
 R_D=e^{-\pi p/(\omega_c\tau_q)}
$$

となる。この $\tau_q$ はエネルギー準位を広げる時間であり、前方散乱の重みを除く伝導の緩和時間と一般には異なる。

軌道が二スピンで同形、Zeeman エネルギーを $\pm g\mu_B B/2$、$\mu_B=e\hbar/(2m_0)$ とする。$m_0$ は裸の電子質量、$g$ はこの向きでの有効な因子である。それぞれのスピン系列に対して $\mu\mapsto\mu\mp g\mu_B B/2$ とすると、局所的に位相が $\mp 2\pi pF'(\mu)g\mu_B/2$ だけずれる。余弦を二つ足してスピン縮退の因子二で割ると

$$
 R_s=\cos\!\left(\frac{\pi p g m_c}{2m_0}\right)
$$

となる。値が零になると同じ極値面積に由来する二スピンの振動が打ち消し合う。以上は閉軌道、弱磁場、滑らかな化学ポテンシャルを仮定した :term[Lifshitz–Kosevich の振幅因子]{en="Lifshitz–Kosevich amplitude factors" yomi="りふしっつこせゔぃちのしんぷくいんし"} である。磁化に加え電気抵抗も振動し得る。抵抗の振動は Shubnikov–de Haas 効果と呼び、周波数が共通でも抵抗には散乱と電流の重みが入る。

:::example[有限温度の基本波と第二高調波]
$m_c=0.20m_0$、磁場 $10\,\mathrm T$、温度 $4.0\,\mathrm K$ とする。基本波の $X$ は約 $1.18$ で $R_T$ は約 $0.80$ である。第二高調波では $X\simeq2.35$、$R_T\simeq0.45$ となる。もとの零温度振幅も高調波ほど小さいので、温度を上げると第二高調波がさらに見えにくくなる。
:::

:::example[スピンによる基本波の消失]
$g=2$、$m_c=m_0/2$ のとき基本波の $R_s=\cos(\pi/2)=0$ である。一方で第二高調波では $\cos\pi=-1$ となり、位相が反転して残る。基本波が見えないことだけで軌道が存在しないと結論することはできない。磁場の向きを変えた場合は、$g$ と $m_c$ の角度依存も点検する。
:::

### 周波数の角度依存から分かる形

式 :eqref[eq-ss-ellipsoid-cyclotron] より、放物線楕円体の中心断面では $F(\vb n)\propto1/\sqrt{\vb n^{\mathsf T}m^*\vb n}$ である。主軸から傾けた測定を含めれば、軸方向の三つの値だけで予測した角度依存を確認できる。ただし交換可能な複数のポケットがある場合、各角度で同じ周波数枝を追跡しなければならない。

完全な円筒で軸が $z$、磁場が軸と角度 $\theta$ をなす場合を考える。断面を $xy$ 平面に投影すると円になり、平面上の面積は投影面積の $1/|\cos\theta|$ 倍である。従って断面が閉じている範囲で

$$
 F(\theta)=F(0)/|\cos\theta|
$$

となる。側面に小さい凹凸があれば太い部分と細い部分が別の周波数を与え、基本波が近い範囲では重ね合わせのうなりとして観測される。各振幅が同じとき、差 $\Delta F$ の二つの余弦は $\cos(\pi\Delta F/B)$ の包絡を持つ。包絡の零は $1/B$ で $1/|\Delta F|$ おきに現れる。どの断面が得られたかは場の向きとポケットの形に依存し、一方向の周波数だけで三次元の全形状が定まるわけではない。

:::example[傾斜と近接周波数]
円筒の軸方向で $120\,\mathrm T$ の周波数が、$\theta=60^\circ$ で $240\,\mathrm T$ となる。別に、同一方向で $110\,\mathrm T$ と $130\,\mathrm T$ の基本波が同振幅で重なるなら、包絡の隣り合う零の間隔は $1/(20\,\mathrm T)=0.050\,\mathrm{T^{-1}}$ である。前者は一つの枝の角度依存、後者は周波数差を測っている。
:::

### 二次元で粒子数を固定した場合

式 :eqref[eq-ss-oscillation-omega] は固定化学ポテンシャルである。孤立した清浄二次元系で総密度 $n$ を固定するなら、温度零では最後の準位を部分充填できる。$g_B=eB/h$、$b=\hbar eB/m$ として $\nu=n/g_B=N+f$、$N=\lfloor\nu\rfloor$、$0\le f<1$ と置く。最低から $N$ 本を満たし、添字 $N$ を割合 $f$ だけ満たす。全エネルギー $\mathcal E$ は

$$
 \frac{\mathcal E}{S}
 =g_B b\left[\frac{N^2}{2}+f(N+1/2)\right]
 =\frac{g_B b}{2}[\nu^2+f(1-f)]
$$

である。第一項は $N$ 本の零点項を含む有限和 $\sum_{j=0}^{N-1}(j+1/2)=N^2/2$ から得た。$\nu\propto B^{-1}$ なので $g_B b\nu^2/2$ は磁場一定の項である。整数を横切らない区間では $\dd f/\dd B=-\nu/B$ だから、面積当たりの磁気モーメントは

$$
 -\partial_B(\mathcal E/S)
 =\frac{g_B b}{2B}\left[\nu(1-2f)-2f(1-f)\right].
$$

この場合、端数 $f$ が $0$ に戻るときに鋸歯状波が跳ぶ。固定 $\mu$ では準位の中心が $\mu$ を通る $\mu/b-1/2\in\Z$ で跳ぶため、両条件で逆磁場の周期が同じでも跳ぶ位置はずれる。三次元でほかの帯も粒子を受け持つ場合には、ひとつのポケットだけで固定密度を課す必要がないことがある。

:::example[二つの測定条件の位相差]
密度を $nh/e=100\,\mathrm T$ とする一スピン系では、固定 $n$ の跳びの位置は $1/B=N/(100\,\mathrm T)$ である。無磁場で同じ密度になる $\mu$ を固定した場合、$\mu=m^{-1}2\pi\hbar^2n$ なので $F=m\mu/(e\hbar)=100\,\mathrm T$ である。跳びの位置は $1/B=(j+1/2)/(100\,\mathrm T)$ に移る。どちらも逆磁場周期は $0.010\,\mathrm{T^{-1}}$ となる。
:::

## 演習問題

::::exercise{level=basic}
$E=\hbar^2(k_x^2+4k_y^2)/(2m)$ の主質量を求めよ。
:::solution
二階微分を $\hbar^2$ で割ると $m^{-1}$ と $4m^{-1}$ が得られる。従って主質量は $m,m/4$ である。
:::
::::

::::exercise{level=basic}
格子定数と一様電場を両方二倍にしたとき、Bloch 周期は何倍か。
:::solution
$\omega_B=e|E_{\rm field}|a/\hbar$ は四倍になるので、周期は $1/4$ 倍になる。
:::
::::

::::exercise{level=basic}
電子で $m_x=m$、$m_y=9m$ のとき $z$ 軸方向の場の $m_c$ を求めよ。
:::solution
式 :eqref[eq-ss-cyclotron-mass] の楕円放物線に対する結果を使い、$\sqrt{m_xm_y}=3m$ である。
:::
::::

::::exercise{level=basic}
Landau 準位の $j=3$ と $j=1$ のエネルギー差を $\hbar\omega_c$ 単位で求めよ。$k_z$ は同じとする。
:::solution
式 :eqref[eq-ss-landau] の並進運動と零点項が消えるため、差は $(3-1)\hbar\omega_c$ である。
:::
::::

::::exercise{level=advanced}
電場がないとき一次元実空間での閉じた Bloch 振動は維持されるか。式 :eqref[eq-ss-semiclassical] に基づいて考えよ。
:::solution
一次元で横方向の運動はなく、電場を零にすれば波数は時間に依存しない。したがって一定の群速度で進むか、その波数が極値なら停止する。周期運動は波数がゾーンを移動して初めて起こる。
:::
::::

::::exercise{level=advanced}
等方二次元分散 $E=\hbar v k$、$E>0$ のサイクロトロン質量を求めよ。
:::solution
$A=\pi[E/(\hbar v)]^2$ なので $A'=2\pi E/(\hbar v)^2$ である。式 :eqref[eq-ss-cyclotron-mass] から $m_c=E/v^2$ となる。分散の動径二階微分は零でも、閉軌道の周期は有限となる。
:::
::::

::::exercise{level=advanced}
電子の円軌道で $k_\perp\ell_B=\sqrt{2j+1}$ を導け。
:::solution
$E_\perp=\hbar^2 k_\perp^2/(2m)=\hbar(eB/m)(j+1/2)$ と置く。両辺を $2m\ell_B^2/\hbar^2$ 倍すれば $k_\perp^2\ell_B^2=2j+1$ である。$k_\perp$ は量子状態の平面波波数でなく、同エネルギーの古典軌道半径を表す。
:::
::::

::::exercise{level=advanced}
磁気振動の周波数が $50\,\mathrm T$ と $200\,\mathrm T$ の二つで、両軌道が円形だとする。半径比を求めよ。
:::solution
式 :eqref[eq-ss-onsager] で面積は周波数に比例する。大きい方の面積が四倍なので半径は二倍である。
:::
::::

::::exercise{level=advanced}
立方格子の基本並進を $a\hat{\vb x},a\hat{\vb y},a\hat{\vb z}$ とし、一様電場の成分を $E_0(2,3,0)$、$E_0>0$ とする。波数軌道の最短の周期を求めよ。$z$ 向きの群速度が非零のまま一定になるバンドの例も挙げよ。
:::solution
波数が戻る条件は、$u=|q|E_0 aT/(2\pi\hbar)$ と置くと $2u,3u\in\Z$ である。この差 $u$ も整数で、最小の正の値は $1$ なので $T=2\pi\hbar/(|q|E_0a)$ になる。バンドに $-2t_z\cos k_z a$ の項を加え、他の項を $k_x,k_y$ のみに依存させる。$k_za=\pi/2$ の初期状態なら $z$ の場が零のため $k_z$ は一定であり、$v_z=2t_za/\hbar$ が保たれる。波数が戻っても実空間で戻らない例である。
:::
::::

::::exercise{level=advanced}
主軸質量 $(m,4m,9m)$ の楕円体で、$(1,0,1)/\sqrt2$ に磁場を向けた場合の $m_c$ を求めよ。同じエネルギーで $z$ 軸方向の場と比較した極値断面に由来する周波数比も求めよ。
:::solution
$\det m^*=36m^3$、$\vb n^{\mathsf T}m^*\vb n=5m$ である。式 :eqref[eq-ss-ellipsoid-cyclotron] より $m_c=6m/\sqrt5$ となる。$z$ 軸方向では $m_c=\sqrt{m\cdot4m}=2m$ である。同じエネルギーでは中心の断面積が $m_c$ に比例するので、傾けた方向と $z$ 軸方向の周波数比は $3/\sqrt5$ になる。
:::
::::

::::exercise{level=advanced}
一つの軌道について $R_D=e^{-\pi/(\omega_c\tau_q)}$ を使う。一スピンで基本波を測り、温度と幾何学に由来する振幅因子を取り除いた比が $B_1=5\,\mathrm T$ と $B_2=10\,\mathrm T$ に対し $R_D(B_1)/R_D(B_2)=1/2$ であった。$m_c$ が既知なら $\tau_q$ を求めよ。
:::solution
$\omega_c=eB/m_c$ なので
$$
 \log\frac{R_D(B_1)}{R_D(B_2)}
 =-\frac{\pi m_c}{e\tau_q}\left(\frac1{B_1}-\frac1{B_2}\right)
$$
である。左辺は $-\log2$、逆磁場の差は $0.10\,\mathrm{T^{-1}}$ である。従って $\tau_q=\pi m_c(0.10\,\mathrm{T^{-1}})/(e\log2)$ を得る。$\pi m_c/(e\tau_q)$ は磁場の単位を持ち、指数は無次元になる。
:::
::::

::::exercise{level=advanced}
半古典量子化が成り立つ円形の電子軌道で、同じ $B$ の一連の Landau 準位を調べる。$\partial_E A$ が一定でない場合にも、隣接準位の間に $m_c(E)$ を使って準位差を表せ。
:::solution
式 :eqref[eq-ss-orbit-quantization] を隣り合う準位で引くと $A(E_{j+1})-A(E_j)=2\pi eB/\hbar$ である。$A$ はエネルギーに対して単調増加とする。微分可能なら平均値の定理から区間内の $E_*$ で $A'(E_*)(E_{j+1}-E_j)=2\pi eB/\hbar$ となる。そこで $A'(E_*)=2\pi m_c(E_*)/\hbar^2$ を使えば
$$
 E_{j+1}-E_j=\frac{\hbar eB}{m_c(E_*)}
$$
である。$m_c$ の変化が小さい範囲では $E_*$ を一方の端のエネルギーで置き換えられる。
:::
::::

::::exercise{level=advanced}
閉軌道が四回回転対称で、速度の複素成分が $v_+(s+T_c/4)=iv_+(s)$ を満たすとする。$v_+(s)$ の Fourier 係数が非零となり得る整数 $p$ を求めよ。
:::solution
各モード $e^{ip\omega_c s}$ は $s$ を $T_c/4$ 進めると $e^{ip\pi/2}$ 倍される。Fourier 係数の一意性から、非零なら $e^{ip\pi/2}=i$ を満たす必要がある。従って $p=1+4j$、$j\in\Z$ である。共役成分 $v_-$ では逆の合同類になるため、偏光によって現れ得る高調波が異なる。
:::
::::

## 確認テスト

::::quiz
静磁場だけの下で半古典電子について保存するものを全て選べ。

- [x] バンドエネルギー
- [x] 磁場と平行な波数成分
- [ ] 波数ベクトル全体

:::explanation
磁場による力は速度に垂直なので仕事をしない。また $\dot{\vb k}$ は磁場にも垂直なので、平行成分は保存する。一般に垂直成分は軌道に沿って変化する。
:::
::::

::::quiz{answer="3" tol="1e-6"}
$t>0$ の余弦帯で $t$ だけを三倍にしたとき、固定電場で Bloch 振動の位置の最大値と最小値の差は何倍か。
:::explanation
位置の差は帯幅 $4t$ を電気力の大きさで割った値である。力が同じなので帯幅の比 $3$ が位置の差の比となる。
:::
::::

::::quiz
磁気振動の主な周波数に対応する幾何量はどれか。

- [x] Fermi 面の、磁場に垂直な極値断面積
- [ ] 実空間の原始胞の面積
- [ ] Fermi 面全体の表面積

:::explanation
個々の断面の位相は $A/B$ に比例する。三次元で平行波数を積分すると停留位相に対応する極値断面が残り、その $A$ が周波数を決める。
:::
::::

::::quiz{answer="0.02" tol="1e-6"}
$F=50\,\mathrm T$ の磁気振動の周期を $\mathrm{T^{-1}}$ 単位で答えよ。
:::explanation
準位の通過は $F/B$ が一ずつ変わると起こる。したがって $1/B$ の周期は $1/F=0.02\,\mathrm{T^{-1}}$ である。
:::
::::

::::quiz{answer="0.25" tol="1e-6"}
単一の Landau–Zener 交差で、あるギャップのとき帯を乗り換える確率は $1/2$ である。掃引速度を変えずギャップを $\sqrt2$ 倍にすると新しい確率はいくらか。
:::explanation
最小ギャップは $2U$ で、式 :eqref[eq-ss-landau-zener] の指数は $U^2$ に比例する。ギャップを $\sqrt2$ 倍すると指数が二倍になるため、新しい確率は $(1/2)^2=1/4$ である。
:::
::::

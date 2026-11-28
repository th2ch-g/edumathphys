---
title: "アノマリー"
description: "軸性カレントの量子補正を三角図と測度の変化から求め、ゲージと尺度のアノマリーを区別する。"
notation:
  - symbol: 'j_5^\mu'
    meaning: 軸性カレント
  - symbol: '\widetilde F^{\mu\nu}'
    meaning: 場の強さの双対
  - symbol: '\mathcal A(\omega)'
    meaning: ゲージ変換に対する有効作用の変化
---

## この章の目標

- 古典のカイラル変換と量子測度の破れを区別する。
- 四次元の軸性アノマリーの係数を二通りで求める。
- 二次元の質量生成とゲージアノマリーの相殺条件を説明する。

## カイラル対称性と軸性カレント

### 二つのカイラル成分の位相

$\gamma^5$ を Dirac 章の規約で取り、$\psi'=e^{i\alpha\gamma^5}\psi$ とする。$\gamma^0$ と $\gamma^5$ が反交換するため $\bar\psi'=\bar\psi e^{i\alpha\gamma^5}$ である。従って運動項は定数 $\alpha$ で不変で、質量項は $\delta\mathcal L_m=-2im\alpha\bar\psi\gamma^5\psi$ となる。局所的な $\alpha(x)$ を用いると運動項に $-(\partial_\mu\alpha)\bar\psi\gamma^\mu\gamma^5\psi$ が現れる。よって $j_5^\mu=\bar\psi\gamma^\mu\gamma^5\psi$ の古典式は $\partial_\mu j_5^\mu=2im\bar\psi\gamma^5\psi$ である。質量ゼロでは右・左の位相を独立に回せる。量子論では作用の変分に加え、正則化された測度がこの変換を保つかを調べる。

### フレーバーの軸性変換と質量行列

同じ電荷を持つ複数の Dirac 場を列 $\Psi$ にまとめ、Hermite 質量行列を $M$ とする。定数の Hermite 行列 $T$ に対して $\delta\Psi=i\alpha T\gamma^5\Psi$、$\delta\bar\Psi=i\alpha\bar\Psi T\gamma^5$ と置く。運動項ではガンマ行列と $\gamma^5$ を移すと定数部分が相殺する。質量の変分は $-i\alpha\bar\Psi\{T,M\}\gamma^5\Psi$。従って古典式は $\partial_\mu j_{5,T}^\mu=i\bar\Psi\{T,M\}\gamma^5\Psi$ となる。例えば $M=\diag(m,0)$、$T=\diag(0,1)$ では質量項が消える。$T=\diag(1,-1)$ なら質量項は $2im\bar\Psi_1\gamma^5\Psi_1$ であり、二つの電荷が等しいと量子の電磁アノマリーは逆符号で相殺する。

### 左右の射影と保存則の独立性

射影 $P_{R,L}=(1\pm\gamma^5)/2$ は $(\gamma^5)^2=1$ から $P_RP_L=0$、$P_R+P_L=1$ を満たす。ただし Dirac 共役について $\overline{P_L\psi}=\bar\psi P_R$ である。この関係から $j^\mu=j_R^\mu+j_L^\mu$、$j_5^\mu=j_R^\mu-j_L^\mu$、$j_{R,L}^\mu=\overline{\psi}_{R,L}\gamma^\mu\psi_{R,L}$ と分解できる。質量は $m(\bar\psi_L\psi_R+\bar\psi_R\psi_L)$ のように左右を結ぶ。従って質量を固定した軸回転と、質量も変換する形式的な対称性は異なる。例えば実の $m$ の代わりに複素数 $m_c$ を $\bar\psi_L\psi_R$ の係数に使うと、軸回転でこの積は $e^{2i\alpha}$ を受ける。$m_c\mapsto e^{-2i\alpha}m_c$ まで含めると古典作用は不変だが、物理的な $m_c$ を固定した Noether 電流の保存を意味しない。

ベクトルと軸性の源を区別して $\mathcal L_{\rm source}=-v_\mu j^\mu-a_\mu j_5^\mu$ とする。質量ゼロで右成分が見る源は $v+a$、左成分が見る源は $v-a$。$a=0$ で同じ Dirac 場の測度をベクトルゲージ変換に適合させると、保存すべき Ward 恒等式は $v$ の変換である。軸変換の二つの測度の Jacobian は同符号で足し合わさるのに対し、ベクトル変換での $\psi$ と $\bar\psi$ は逆の位相を受ける。外部源としての軸変換の破れと、$v_\mu$ を動的に積分する理論のゲージ不変性の破れでは、量子論への影響が異なる。

:::example[複数電荷で保存する軸方向を求める]
電荷 $q_1=e,q_2=2e$ の質量ゼロ Dirac 場二つを考える。$T=\diag(t_1,t_2)$ に対し、二つのベクトル頂点の電荷は $q_i^2$ なので軸係数は $e^2(t_1+4t_2)$ に比例する。例えば $T=\diag(4,-1)$ の軸方向なら係数は零。$T=\diag(1,-1)$ はフレーバーのトレースが零でも電荷の二乗で重みを付けると $-3e^2$ となる。さらに第一成分だけに質量 $m$ があれば、$T=\diag(4,-1)$ でも質量の寄与 $8im\bar\Psi_1\gamma^5\Psi_1$ は残る。
:::

### 局所対称性の Ward 式における接触点

電磁場を固定したまま、$O=\psi(y)\bar\psi(z)$ を量子相関へ挿入する。局所軸変換で $\delta O=i\alpha(y)\gamma^5 O+i\alpha(z)O\gamma^5$。作用の変分は部分積分により $\delta S=\int\dd^4x\,\alpha(x)[\partial_\mu j_5^\mu-2imP]$。Minkowski の変数変換で $0=\langle\delta O+iO\delta S+O\delta\log\mathcal J\rangle$ と書けば、分離点で得たアノマリーを $\mathcal X(x)$ として
$$
\begin{aligned}
\langle T(\partial_\mu j_5^\mu-2imP-\mathcal X)(x)\,\psi(y)\bar\psi(z)\rangle
={}&-\delta^4(x-y)\gamma^5\langle T\psi(y)\bar\psi(z)\rangle\\
&-\delta^4(x-z)\langle T\psi(y)\bar\psi(z)\rangle\gamma^5 .
\end{aligned}
$$
ここでは左辺の発散を時間順序相関の $x$ 微分として読む。右辺は外脚の変換から生じ、二点を離せば消える。例えば外部場が零で $m=0$ でも、$x$ と挿入が接触するところまで発散が零になるわけではない。ベクトル回転では $\bar\psi$ の変分の符号が逆なので、対応する二つの接触項も相対的に逆になる。この違いがベクトル電荷の Ward 式と軸電荷の式の最初の比較になる。

### 電荷行列とフレーバー変換

電荷をフレーバーの Hermite 行列 $Q$ として、$T$ が $Q$ と可換な軸変換を考える。古典の運動項で背景結合を含めて不変なのはこの可換性による。$[Q,T]\ne0$ なら、$\alpha$ が定数でも運動項に背景場を含む追加の変分がある。例えば $Q=\diag(q_1,q_2)$ に $T=\tau^1/2$ を取れば $q_1\ne q_2$ でその変分は零でない。逆に $[Q,T]=0$ なら、二つのベクトル頂点と軸頂点のフレーバートレースは $\tr(TQ^2)$。質量行列の $\{T,M\}$ と独立に、この量がゲージ背景で保存する軸方向を指定する。大域フレーバーの軸回転が動的なゲージ場との結合で破れることと、ゲージ変換自身が保てないことは、この二つの Ward 式を別々に取ると区別できる。

## 三角ダイアグラムと Adler–Bell–Jackiw アノマリー

### ベクトル電流を保存する正則化

電荷 $q$ の Dirac 場を外部電磁場へ結合する。軸性一つとベクトル二つの挿入を持つ三角図は表面的に一次発散するので、各図を個別に運動量シフトして三つの Ward 恒等式を同時に課せない。Pauli–Villars 場を同じベクトル電荷と大きい質量 $M$ で導入し、ループを差し引く。ベクトル挿入は隣接逆伝播関数の差になるため、差し引いた収束積分の変数を移して二つとも保存する。軸性挿入は $\gamma^5$ の反交換により逆伝播関数の和に加え、質量挿入 $2imP$ を生む。正則化した発散は従って物理場の質量項から重い場の同型の項を引いたものである。

重い場の擬スカラー頂点の二つの順序を足すと、紫外の残りは場の強さの二乗から作る擬スカラーになる。以下で分子と積分を別々に求める。$\epsilon^{0123}=1$、$\gamma^5=i\gamma^0\gamma^1\gamma^2\gamma^3$ とし、重い場を引く符号を含めると
$$
\partial_\mu j_5^\mu=
2im\bar\psi\gamma^5\psi-\frac{q^2}{16\pi^2}\epsilon^{\mu\nu\rho\sigma}F_{\mu\nu}F_{\rho\sigma}.
\label{eq-qft-abj}
$$
この残りを :term[軸性アノマリー]{en="axial anomaly" yomi="じくせいあのまりー"} という。$\widetilde F^{\mu\nu}=\epsilon^{\mu\nu\rho\sigma}F_{\rho\sigma}/2$ として $F\widetilde F=-4\vb E\cdot\vb B$ なので、電磁項は $q^2\vb E\cdot\vb B/(2\pi^2)$ である。重い調整場が消えても、この項は $M^2$ と分母積分の $M^{-2}$ が相殺するため残る。

### 局所な移し替えと二つのベクトル Ward 恒等式

二つの同種ベクトル挿入に入る運動量を $k_1,k_2$ とする。一次発散の三角振幅のルーティングの差は局所な一次多項式になり得る。軸性添字を $\lambda$ とし、Bose 交換に不変な擬テンソルとして $\Delta T^{\lambda\mu\nu}=c\epsilon^{\lambda\mu\nu\rho}(k_1-k_2)_\rho$ が許される。実際、二つのベクトル添字の交換と運動量差の交換がともに負号を出す。$k_{1\mu}$ を縮約すれば $k_{1\mu}k_{1\rho}$ の項は消えるが $k_{1\mu}k_{2\rho}$ は残る。従ってベクトル Ward 恒等式を満たしている振幅には、任意の $c$ を加えられない。軸性側を保存するためにこの項を選び直すと、一般にベクトル側で保存が破れる。

### 調整場の三点積分を最後まで評価する

この計算では $D_\mu=\partial_\mu+iqA_\mu$ とし、光子の背景への相互作用は $-qA_\mu j^\mu$。二光子核は $A_\mu(k_1),A_\nu(k_2)$ で背景期待値を微分したものとする。外場の Fourier 波を $e^{-ikx}$ と書き、三つの内部運動量を $a=\ell+k_1,b=\ell,c=\ell-k_2$ と取る。擬スカラー挿入 $P=\bar\psi\gamma^5\psi$ を先頭にした一順序の分子は
$$
N^{\mu\nu}=\tr[\gamma^5(\not a+M)\gamma^\mu(\not b+M)\gamma^\nu(\not c+M)] .
$$
四ガンマだけが残り、それぞれの伝播分子から質量を選ぶ三通りがある。$\epsilon(\mu,\nu,u,v)=\epsilon^{\mu\nu\rho\sigma}u_\rho v_\sigma$ とおく。$a,c$ を含む項のガンマ行列の順は $(a,\mu,\nu,c)$、$b,c$ の項は $(\mu,b,\nu,c)$、$a,b$ の項は $(a,\mu,b,\nu)$。従って合計は $-4iM[\epsilon(\mu,\nu,a,c)-\epsilon(\mu,\nu,b,c)-\epsilon(\mu,\nu,a,b)]$。角括弧は $-\epsilon(\mu,\nu,a-b,b-c)$ に等しいので
$$
N^{\mu\nu}=4iM\,\epsilon^{\mu\nu\rho\sigma}k_{1\rho}k_{2\sigma}.
$$
分子から $\ell$ が消えるのは三つの項を足した結果である。二光子を入れ替えた順序は $\mu,\nu$ と $k_1,k_2$ がともに交換され、同じ符号を与える。以下で擬スカラーの挿入への運動量の向きも含め、分子と局所場の Fourier 微分を同じ規約で比較する。

背景中の Dirac Green 関数は $i(\not p-M-q\not A)^{-1}$。$\bar\psi\gamma^5\psi$ の一点関数では、縮約を $\langle\psi\bar\psi\rangle$ の順へ戻す負号が出る。背景二次は $-i q^2$ を上のトレースと分母積分に掛け、$2iM P$ ではさらに $2iM$ を掛ける。三分母積分は
$$
\int\frac{\dd^4\ell}{(2\pi)^4}
\frac1{(a^2-M^2+i0)(b^2-M^2+i0)(c^2-M^2+i0)}
=-\frac{i}{32\pi^2 M^2}+\mathcal O(M^{-4}).
$$
Feynman パラメータで $u+v+w=1$、$u,v,w\ge0$ として平方完成すると、有効質量二乗は $M^2-uvk_1^2-vwk_2^2-uw(k_1+k_2)^2$ である。外部運動量固定の極限では先頭で $M^2$ を使え、パラメータ領域の面積 $1/2$ は結合公式の係数 $2$ と消える。三つの Minkowski 分母の負号と時間積分の回転の $i$ から $-i I_3$ となる。Euclid の $I_3$ は $2\pi^2\int_0^\infty r^3\dd r/(2\pi)^4(r^2+M^2)^3$ であり、$u=r^2/M^2$ で $\int_0^\infty u(1+u)^{-3}\dd u=1/2$ となる。

以上から二つの順序の $2iM P$ の核は $16q^2 M^2 I_3\,\epsilon(\mu,\nu,k_1,k_2)$ の正符号を持つ。

同じ背景で $\epsilon F F$ を二つの $A$ で微分すると $8\,\epsilon(\mu,\nu,k_1,k_2)$。係数八は各 $F$ の反対称部分で二、両 $F$ の交換で二、もう一つの反対称部分で二である。二つの $-i$ と Levi–Civita の並べ替えの負号が打ち消し合う。よって重い場では $\lim_{M\to\infty}\langle2iMP\rangle=(q^2/(16\pi^2))\epsilon F F$。Pauli–Villars の差ではこれを引くので :eqref[eq-qft-abj] の符号になる。異なる符号の調整場を追加してさらに高い紫外次数を改善しても、ベクトルを保った Ward 恒等式が同じ低エネルギーの局所項を固定する。

### 軸性 Ward 恒等式の分子

因子 $i$ を外した逆フェルミオン線を $S_m^{-1}(p)=\not p-m$ とする。$k=p'-p$ なら
$$
\not k\gamma^5=S_m^{-1}(p')\gamma^5+\gamma^5 S_m^{-1}(p)+2m\gamma^5 .
$$
最初の二項は隣の伝播関数を約し、二つの二分母積分の差を作る。順序を交換した図も合わせて Pauli–Villars の差を取ると、軸頂点でのこの差はループ変数をずらして消せる。最後の項が質量付き擬スカラー挿入になる。ベクトル側の $\not k=S_m^{-1}(p')-S_m^{-1}(p)$ にはそもそも質量の残りがない。従って一つの正則化の中で、二方のベクトル側と軸側の違いを同時に確認できる。

### 中性擬スカラーから光子二個への低エネルギー結合

低エネルギーに軸電流と重なる擬スカラー $\pi$ があり、近似対称性で $j_5^\mu=f\partial^\mu\pi+\cdots$ となるとする。光子背景における発散を $a F_{\mu\nu}\widetilde F^{\mu\nu}$ と書けば、角方向の方程式がこれを再現するための局所項は $(a/f)\pi F_{\mu\nu}\widetilde F^{\mu\nu}$。実際 $\mathcal L=(\partial\pi)^2/2+(a/f)\pi F\widetilde F$ なら $f\Box\pi=aF\widetilde F$ となる。四次元で $[\pi]=[f]=1$ なので、係数 $a/f$ は次元 $-1$。有質量の擬 Goldstone として $\pi$ を扱う場合は $\partial_\mu j_5^\mu=-fm_\pi^2\pi+aF\widetilde F+\cdots$ となる。例えば一つのフレーバー方向の係数が $\sum_f T_f Q_f^2=0$ なら、この経路による先頭の二光子結合は零。二光子崩壊が高い次数や他の演算子からも全て禁じられることまでは従わない。

### 局所対項でどの Ward 式が変わるか

外部の Abel ベクトル源 $v$ と軸源 $a$ を一次形式として書く。例えば $B=c\int a\wedge v\wedge\dd v$ は四次元の局所な対項である。境界項を捨てられるとし、$\delta_\beta v=\dd\beta$、$\delta_\alpha a=\dd\alpha$ とすると
$$
\delta_\beta B=c\int\beta\,\dd a\wedge\dd v,\qquad
\delta_\alpha B=-c\int\alpha\,\dd v\wedge\dd v .
$$
最初の符号は $\dd(a\beta\dd v)=\beta\dd a\wedge\dd v-a\wedge\dd\beta\wedge\dd v$ から得る。後者では $\dd(\alpha v\wedge\dd v)=\dd\alpha\wedge v\wedge\dd v+\alpha\dd v\wedge\dd v$ を使う。従って軸発散の $vv$ 項の係数を移すと、一般の軸源に対してベクトル Ward 恒等式の $av$ 項も動く。$a=0$ に制限して二点だけを計算するとこの変化が見えないが、$a$ で微分して三点を取ると現れる。動的にゲージ化するベクトル対称性の Ward 式を保つには、含める局所対項をこの条件に合わせる。

### 二光子結合と崩壊幅

二つの等しい力学的質量を無視できるフレーバーに $T=\diag(1/2,-1/2)$、電荷を $Q_u=2/3,Q_d=-1/3$ とする。色三を掛けたトレースは $3[(1/2)(4/9)-(1/2)(1/9)]=1/2$。:eqref[eq-qft-abj] を $\epsilon F F=2 F\widetilde F$ で書くと、この電流の電磁係数は $a=-e^2/(16\pi^2)$。これを $\langle0|j_5^\mu|\pi(p)\rangle$ の大きさを $fp^\mu$ と取る低エネルギーの表示へ移す。$\mathcal L_{\pi\gamma\gamma}=-(g_{\pi\gamma\gamma}/4)\pi F\widetilde F$ なら $g_{\pi\gamma\gamma}=e^2/(4\pi^2 f)$。擬スカラーの電荷と質量の規約はここで固定し、物理的な $f$ は電流の行列要素で定義する。

二光子の反対向きの運動量を $k_1,k_2$ とすると、振幅の大きさは $g_{\pi\gamma\gamma}$ 倍の $\epsilon^{\mu\nu\rho\sigma}\epsilon_{1\mu}\epsilon_{2\nu}k_{1\rho}k_{2\sigma}$。静止系で二つの直交する横線形偏極の組だけが寄与し、各組で絶対値 $g_{\pi\gamma\gamma}m_\pi^2/2$。従って偏極の和は $g_{\pi\gamma\gamma}^2m_\pi^4/2$。二体位相空間 $1/(8\pi)$、崩壊率の $1/(2m_\pi)$ と同種光子の $1/2!$ を掛けると
$$
\Gamma_{\pi\to\gamma\gamma}
=\frac{g_{\pi\gamma\gamma}^2m_\pi^3}{64\pi}.
$$
例えば $f$ を固定して小さい明示的破れで $m_\pi$ が二倍になれば、この先頭幅は八倍になる。ここでは電流の低エネルギー極でアノマリーを飽和する先頭近似を使っており、有限質量や高次の強い相互作用による修正まで同じ係数とする必要はない。

### 一ループの厳密性の範囲

ベクトル Ward 恒等式を保つくりこみで、正規化した軸電流のアノマリーの係数は対応する場の強さの正規化に合わせれば一ループで定まる。これは :term[Adler–Bardeen の定理]{en="Adler–Bardeen theorem" yomi="あどらーばーでぃーんのていり"} の内容である。全次数の証明には複合電流の有限規格化と局所な Ward 恒等式の同時制御が必要なので :cite[weinberg-qft2] を参照する。正則化した三角図と Fujikawa の計算で一致を確認しても、それだけで一般の時間順序三点関数の横部分や、軸電流の任意の基底での異常次元まで消えるわけではない。例えば $\epsilon F F$ と電流の発散が局所演算子のくりこみで混ざる処方では、係数を移す有限な再規格化と定理の適用条件を明示する必要がある。

## Fujikawa の経路積分による導出

### 固有モード上の Jacobian

Euclid の Hermite 行列 $\gamma_E^\mu$ を取り、$\epsilon_E^{1234}=1$、$\gamma_*=-\gamma_E^1\gamma_E^2\gamma_E^3\gamma_E^4$ とする。$\not D_E$ の直交固有モードで二つの Grassmann 測度を展開する。変換の行列はどちらも無限小で $1+i\alpha\gamma_*$ なので、:term[Fujikawa の方法]{en="Fujikawa method" yomi="ふじかわのほうほう"} では $\log\mathcal J=-2i\operatorname{Tr}[\alpha\gamma_*e^{\not D_E^2/M^2}]$ と正則化する。逆行列式を二回使ったためである。$\not D_E^2=D_E^2+\gamma_E^\mu\gamma_E^\nu\Omega_{\mu\nu}/2$、$\Omega=-igF$ と書く。$\tr_{\rm spin}\gamma_*$ と二ガンマのトレースは零なので最初に残るのは曲率スピン項の二乗である。その項は $(-ig)^2/(8M^4)$ に四ガンマのトレース $-4\epsilon_E^{\mu\nu\rho\sigma}$ を掛ける。平面波で $\int\dd^4k\,e^{-k^2/M^2}/(2\pi)^4=M^4/(16\pi^2)$。従って
$$
\lim_{M\to\infty}\tr\bigl[\gamma_*e^{\not D_E^2/M^2}\bigr](x,x)
=\frac{g^2}{32\pi^2}\epsilon_E^{\mu\nu\rho\sigma}\tr_R(F_{\mu\nu}F_{\rho\sigma}).
\label{eq-qft-fujikawa}
$$
である。軌道の熱核にだけ現れる曲率二乗はスピン単位行列なので寄与しない。$t=-i\tau$、$\gamma_E^j=-i\gamma^j$、$\gamma_E^4=\gamma^0$ では $\gamma_*=\gamma^5$、$\partial_E\cdot j_{5,E}=-i\partial\cdot j_5$。$F^E_{i4}=iF_{0i}$ から $\epsilon_E F_E F_E=i\epsilon F F$。この関係と二測度の係数二を使うと :eqref[eq-qft-abj] が得られる。積分した熱核のトレースでは非零固有値が $\gamma_*$ によって反対符号の固有空間と組になるので、カイラルな零モード数の差だけが残る。

### 平行な電場と磁場で係数を確かめる

一定の $\vb B=B\widehat{\vb z}$ に対し、二成分の質量ゼロ場の Hamiltonian は右手で $\boldsymbol\sigma\cdot\boldsymbol\pi$、左手でその負である。運動学的運動量には $[\pi_x,\pi_y]=iqB$ が成り立つので
$$
(\boldsymbol\sigma\cdot\boldsymbol\pi)^2
 =\pi_x^2+\pi_y^2+p_z^2-qB\sigma_z .
$$
横振動子のエネルギー二乗は $(2n+1)|qB|$、スピン項は $-qB\sigma_z$。$n=0$ かつ $\sigma_z=\sgn(qB)$ の成分だけが横ギャップ零である。従って最下準位の分散は右手で $\sgn(qB)p_z$、左手でその負。面積あたりの縮退 $|qB|/(2\pi)$ は横振動子の中心 $p_y/(qB)$ を箱内に取ったときのモード数から出る。平行な電場 $E$ は $\dot p_z=qE$ を与えるので、軸電荷の密度変化率は二分枝を数えて
$$
\frac{2\,\sgn(qB)\,qE\,|qB|}{(2\pi)^2}
=\frac{q^2 EB}{2\pi^2}.
$$
高い Landau 準位には左右に零を跨ぐ分枝がない。従ってこの一様背景の計算は :eqref[eq-qft-abj] の符号と係数を独立に与える。

### 曲率の積分と正負の零モード

Euclid の $\not D_E$ を反 Hermite とし、$i\not D_Eu_\lambda=\lambda u_\lambda$ とする。$\gamma_*u_\lambda$ の固有値は $-\lambda$。$\lambda\ne0$ なら直交性から $\inner{u_\lambda}{\gamma_*u_\lambda}=0$。従って有限体積で固有値和を取った $\Tr(\gamma_* e^{\not D_E^2/M^2})$ は調整尺度によらず $n_+-n_-$ となる。$\lambda=0$ の空間では $\gamma_*$ と同時対角化し、カイラリティの固有値ごとに数えた。:eqref[eq-qft-fujikawa] を積分して $SU(N)$ の基本表現では $g^2\int F^a\widetilde F^a/(32\pi^2)$ が得られる。例えば自己双対の電荷一の背景では $n_+-n_-=1$ となり、二つの核がすべて対で相殺することはない。ここではコンパクト化または同等の境界条件を仮定し、境界の寄与を含む場合はスペクトルの扱いを分ける。

### 局所な変換と正則化した跡

完全系 $u_n$ を取り $\psi=\sum_n a_nu_n$、$\bar\psi=\sum_n\bar b_n u_n^\dagger$ とする。無限小の係数変換は $a_n'=\sum_m(\delta_{nm}+iA_{nm})a_m$、$A_{nm}=\int u_n^\dagger\alpha(x)\gamma_*u_m$。$\bar b$ は右から同じ行列で変換される。有限モードで $\int\dd a\,a=1$ を保つには $\dd(ka)=k^{-1}\dd a$ なので、双方から $\det(1+iA)^{-1}$。線形に展開して $\delta\log\mathcal J=-2i\Tr A$ を得る。点ごとの $\tr\gamma_*$ が零でも、モード和と局所極限の交換は無条件にできない。$\not D_E$ の固有値を使うカットオフでは正則化核がゲージ変換で共役され、$\alpha$ が群の単位行列なら挿入の色トレースがゲージ不変になる。

Euclid の変数変換で、質量を含む作用の古典変化は
$$
\delta S_E=i\int\dd^4x\,\alpha(-\partial_E\cdot j_{5,E}+2m\bar\psi\gamma_*\psi)
$$
と書ける。ここでは $S_E=\int\bar\psi(\not D_E+m)\psi$ である。挿入が変換で不変か、挿入との接触点を除くと $0=\langle \delta\log\mathcal J-\delta S_E\rangle$。従って $\partial_E\cdot j_{5,E}=2m\bar\psi\gamma_*\psi+2\tr[\gamma_*e^{\not D_E^2/M^2}](x,x)$ の極限を取る。Minkowski へ戻すと軸性発散には全体の $i$ が掛かり、Euclid の $\epsilon_E F_E F_E$ も $i\epsilon F F$ なので曲率二乗に負号が付く。ベクトル回転では $\bar b$ が逆行列で変換され、二つの調整した跡が消える。

:::example[曲率を二平面に置く]
Euclid の $F_{12}=B_1,F_{34}=B_2$ の Abel 背景では $\epsilon_E F_EF_E=8B_1B_2$。:eqref[eq-qft-fujikawa] の局所密度は $q^2 B_1B_2/(4\pi^2)$。辺長を $L_i$ とする四次元トーラスで量子化されたフラックス $q B_1L_1L_2=2\pi n_1$、$q B_2L_3L_4=2\pi n_2$ を選ぶと、体積積分は $n_1n_2$。各平面の横 Landau 状態数 $|n_i|$ と、両磁場の向きで決まるカイラリティの積に対応する。一方の $B_i$ が零なら局所密度は消える。これは $\tr\gamma_*$ に必要な四つの異なるガンマが揃わないことにも対応する。
:::

### 質量ゼロ極限と零モードのある積分

有限のコンパクト Euclid 空間で $\det(\not D_E+m)$ を考える。非零の対は固有値 $i\lambda,-i\lambda$ をもち、積は $\lambda^2+m^2$。零モード一つは因子 $m$ だけを与える。従って固定した背景の行列式は $m^{n_++n_-}$ に比例し、$n_++n_->0$ のセクターの分配関数の重みは質量ゼロ極限で消える。しかし $\langle\psi(x)\bar\psi(y)\rangle_A=(\not D_E+m)^{-1}(x,y)$ は零モードから $u_0(x)u_0^\dagger(y)/m$ を含む。行列式と合わせると、零モードが一つのときはこの挿入が有限になり得る。Grassmann 積分で見ると、零モード対 $a_0,\bar b_0$ を作用の質量項か挿入で一度ずつ飽和する必要があるためである。例えば独立な質量ゼロ Dirac 場が二つあり両方で零モードが一つなら、両フレーバーの零モードを含む四フェルミオン挿入が先頭になる。質量ゼロの真空行列式が零という計算だけで、そのセクターの全ての相関を零とすることはできない。

### 熱核の曲率展開を固定した背景で照合する

一つの色成分の Euclid 接続を可換にし、$F_{12}=B_1,F_{34}=B_2$ を定数とする。ここでは $D=\partial+iqA$ として $\Omega_{12}=iqB_1,\Omega_{34}=iqB_2$。そのスピン部分は $S=iq B_1\gamma^1_E\gamma^2_E+iq B_2\gamma^3_E\gamma^4_E$。$H_1=i\gamma^1_E\gamma^2_E,H_2=i\gamma^3_E\gamma^4_E$ とすると $H_i^2=1$、$[H_1,H_2]=0$、$H_1H_2=\gamma_*$。四次元の各組 $(H_1,H_2)=(\pm1,\pm1)$ の次元は、射影 $(1\pm H_1)(1\pm H_2)/4$ の跡が一なので一である。よって
$$
\tr_{\rm spin}\gamma_*e^{tS}
=4\sinh(tqB_1)\sinh(tqB_2).
$$
一方軌道熱核は二つの Landau 平面の積で、
$$
(4\pi t)^{-2}\prod_{i=1}^2\frac{tqB_i}{\sinh(tqB_i)} .
$$
両者の積は $q^2B_1B_2/(4\pi^2)$ となる。$B_i$ の符号を逆にするとスピンの超跡の符号だけが変わり、軌道の密度に絶対値が入ることと整合する。短時間の一般展開でスピン曲率二乗だけが残ることを、この定数背景では有限の $t$ でも直接確認できる。

### 自己双対と零モードのカイラリティ

$D$ を反 Hermite、各 $\gamma_E$ を Hermite とし、$-\not D_E^2=-D^2-S$ とする。$\gamma_*=-\gamma_E^1\gamma_E^2\gamma_E^3\gamma_E^4$ の規約では、例えば $\gamma_E^1\gamma_E^2\gamma_*=\gamma_E^3\gamma_E^4$。従って自己双対の曲率のスピン作用 $S$ は、$\gamma_*=-1$ の部分で消える。同じ関係は他の二平面の組にも成り立つ。例えば $F_{12}=F_{34}$ なら $\gamma_E^1\gamma_E^2+\gamma_E^3\gamma_E^4=\gamma_E^1\gamma_E^2(1+\gamma_*)$。正規化可能な負カイラリティの零モードがあれば
$$
0=\langle u,-\not D_E^2u\rangle=\langle u,-D^2u\rangle
=\sum_\mu\|D_\mu u\|^2
$$
となるので共変に定数である。$\R^4$ の二乗可積分条件ではそのノルム密度が一定なら零しかない。従って、その境界条件で自己双対の電荷一の背景は正カイラリティの零モードを一個持つ。反自己双対では符号が逆。この結論は索引の差だけを計算した場合より強く、もう一方が零であることに曲率の自己双対性を使っている。

## 2 次元のアノマリーと Schwinger 模型

### 二次元では電場が有質量の波になる

$1+1$ 次元で $\eta=(+,-)$、$\gamma^5=\gamma^0\gamma^1$、$E=F_{01}$ とし、この節では物理電荷 $q$ が Maxwell 方程式の源 $qj^\mu$ になるよう $D_\mu=\partial_\mu+iqA_\mu$ を選ぶ。質量ゼロの円周 $L$ で一様電場をゆっくり変えると、運動量モードの間隔は $2\pi/L$、運動学的運動量の変化率は $qE$ である。右向き分枝の準位が零を跨ぐと軸電荷が一つ増え、左向き分枝の穴も軸電荷を一つ増やす。従って軸電荷の変化率を長さで割ると $\partial_\mu j_5^\mu=qE/\pi$ となる。有限の $L$ では準位交差は離散であるが、大きな $L$ の単位長さあたりの流れがこの式を与える。

質量ゼロの QED を :term[Schwinger 模型]{en="Schwinger model" yomi="しゅうぃんがーもけい"} という。ガンマ積から $j_5^0=j^1$、$j_5^1=j^0$。$\partial_\mu F^{\mu\nu}=qj^\nu$ より $j^0=\partial_x E/q$、$j^1=-\partial_t E/q$ なので、軸性発散は $-\Box E/q$。アノマリーと合わせると $(\Box+q^2/\pi)E=0$ となる。二次元で $[q]=1$ であり、ゲージ不変な電場の励起に質量 $|q|/\sqrt{\pi}$ が現れる。四次元の Maxwell 場の偏極数をこの模型に移すことはできない。

### 静電源の遮蔽

Schwinger 模型に外部保存電流 $J^\mu$ を加え $\partial_\mu F^{\mu\nu}=q j^\nu+J^\nu$ とする。静的な電荷密度 $J^0=\rho$、$J^1=0$ では $\partial_x E=q j^0+\rho$。軸性発散 $\partial_x j^0=qE/\pi$ を代入して $(-\partial_x^2+m_\gamma^2)E=-\rho'$、$m_\gamma^2=q^2/\pi$ となる。$E=-\partial_x A^0$ とし、無限遠で $A^0$ を零に取れば $(-\partial_x^2+m_\gamma^2)A^0=\rho$。例えば $\rho=Q\delta(x)$ の解は $A^0=Qe^{-m_\gamma|x|}/(2m_\gamma)$。原点を除く指数解を連続に貼り、$-\partial_x A^0$ の跳びを $Q$ とすることから係数が決まる。無質量の裸 Maxwell 場での一定電場とは異なり、この模型の点源応答は減衰する。

### 保存則から電場の二次作用を読む

$1+1$ 次元では保存するベクトル電流を局所的に $j^0=\partial_x\varphi/\sqrt\pi$、$j^1=-\partial_t\varphi/\sqrt\pi$ と書ける。ここでは電流の二点関数を再現する正準スカラー変数を用い、零モードには境界条件を別に指定する。このとき $j_5^0=j^1,j_5^1=j^0$ から $\partial_\mu j_5^\mu=-\Box\varphi/\sqrt\pi$。Minkowski で
$$
\mathcal L_{\rm bos}=\frac12(\partial\varphi)^2+\frac12 E^2-\frac{q}{\sqrt\pi}E\varphi
$$
と書けば、$\Box\varphi=-qE/\sqrt\pi$ がアノマリーを与える。最後の相互作用は $-q A_\mu j^\mu$ を部分積分して得る。$A_\mu$ の運動方程式は $E-q\varphi/\sqrt\pi$ が時空定数であることを意味する。定数を固定した電束セクターで $E=q\varphi/\sqrt\pi$ の部分を代入すると、有効ポテンシャルは $q^2\varphi^2/(2\pi)$。従って正準なモードの質量が $|q|/\sqrt\pi$ と分かる。ここで使った局所電流の表示は、フェルミオン演算子を含む全ての相関の同値性の証明までは含まない。

:::example[二つの静電荷間のエネルギー]
外部の二電荷を $\rho(x)=Q[\delta(x)-\delta(x-R)]$ とし、ゲージ結合 $q$ の動的フェルミオンを含める。静的な線形応答の Green 関数は $G(x)=e^{-m_\gamma|x|}/(2m_\gamma)$。電荷を零から一定比率で加える仕事は $\frac12\int\dd x\,\rho A^0$ なので
$$
V(R)=Q^2\bigl(G(0)-G(R)\bigr)=
\frac{Q^2}{2m_\gamma}(1-e^{-m_\gamma R}) .
$$
$R\ll 1/m_\gamma$ では $Q^2R/2+\mathcal O(R^2)$、$R\gg 1/m_\gamma$ では独立な遮蔽電荷の自己エネルギーの和に飽和する。短距離の線形性だけを見て長距離にも閉じ込めを結論することはできない。
:::

### フレーバー数を増やした場合

質量ゼロの $N_f$ 種で物理電荷を $q_i$ とする。各 $j_{5,i}$ の発散は $q_iE/\pi$。時間・空間で Maxwell の二式をそれぞれ微分してまとめると $-\Box E=\sum_i q_i\partial_\mu j_{5,i}^\mu$ であり、質量二乗は $\sum_iq_i^2/\pi$。正準スカラーの局所表示では結合するのは $\sum_iq_i\varphi_i$ だけだから、質量行列は $q_i q_j/\pi$。例えば $q_1=q_2=q$ なら $\varphi_+=(\varphi_1+\varphi_2)/\sqrt2$ が質量 $\sqrt2|q|/\sqrt\pi$ を持ち、直交する差はこの二次表示で質量ゼロである。軸性の差のカレントには両者が逆符号で入るので、同じアノマリーが相殺する。

### 二次元真空偏極を使った質量の照合

Euclid 運動量 $Q^2>0$ における純粋な一フェルミオンループを取る。電荷の規約はこの節と同じで、真空偏極を $(Q^2\delta_{\mu\nu}-Q_\mu Q_\nu)\Pi_E(Q^2)$ と置く。正則化はベクトル Ward 恒等式を保つものを選ぶ。Dirac 行列の二次元の跡は $\tr1=2$、四行列の跡は計量の三つの積にこの二を掛ける。二分母をパラメータ化し、分子の $\ell$ を $r-xQ$ に移すと、$Q_\mu Q_\nu$ の部分は $-4x(1-x)Q_\mu Q_\nu$。従って横テンソルのスカラー係数は
$$
\Pi_E(Q^2)=4q^2\int_0^1\dd x\,x(1-x)
\int\frac{\dd^2r}{(2\pi)^2}
 \frac{1}{[r^2+x(1-x)Q^2]^2}
=\frac{q^2}{\pi Q^2}.
$$
最後の積分は $1/(4\pi x(1-x)Q^2)$ であり、二つの端点でもパラメータの積は一定へ延びる。残りの計量項は Ward 恒等式で $Q^2$ 倍の形に固定される。従って横の逆伝播は $Q^2+q^2/\pi$。軸性発散から求めたゲージ不変な電場の質量と一致する。ベクトル Ward 式を破る切断を用いると余分な局所な計量項が入るため、この証明で $Q_\mu Q_\nu$ の係数だけから全二点を定める手順が使えない。

質量 $m$ を入れたループで同じ $Q_\mu Q_\nu$ の係数を取ると、分母に $m^2$ が加わる。$Q^2\ll m^2$ で $\Pi_E\simeq q^2/(6\pi m^2)$、すなわちこの一ループの横の逆伝播の補正は $Q^2$ に比例する。同時に軸性の Ward 式には $2imP$ が加わるので、電場だけの自由な有質量方程式は閉じない。$m>0$ の強結合の物理スペクトルをこの低運動量の一ループだけで決めるものではない。例えば $m\to0$ を先に取った上の $\Pi_E$ は $Q^2=0$ で非解析になり、二つの極限の順序によって近似が異なる。

### フラックスの挿入と局所保存の対応

円周上の物理電荷 $q$ と一定の $A_1(t)$ を取ると、正準運動量の量子数 $p_n=2\pi(n+\nu)/L$ は一定に保たれ、運動学的運動量は $p_n+qA_1$。ここでは $A_1$ は下添字で、$E=F_{01}=\partial_t A_1$。一周期の大きなゲージ変換で $q A_1 L$ は $2\pi$ だけ動くため、変換後のスペクトル全体は添字を一つずらせば同じになる。しかし占有を時間に沿って追うと、片方の分枝の粒子と他方の穴が残る。このため古典の背景の初期と終期がゲージ同値であることと、軸電荷が二だけ変化することが両立する。ベクトル数は二つの変化で相殺する。$\nu=1/2$ を選べば境界に零モードが来る場合を避けられるが、整数回の全交差数は $\nu$ の一般値でも同じである。

## 非可換ゲージ理論のアノマリーと Wess–Zumino の無矛盾条件

### ゲージ変換そのものの整合性

Dirac 場のベクトル結合では左右の異符号の寄与が打ち消し合う。カイラルに結合すると三角図の対称な色部分は $\tr_R[T^a\{T^b,T^c\}]$ となる。二つのベクトル挿入の順序を交換して足すためである。右手 Weyl 場は左手場と逆符号なので、ゲージ理論の局所摂動アノマリーが消える条件は各 $a,b,c$ について左手の和と右手の和が一致することである。Abel では電荷の三乗の和になる。例えば左電荷 $1,1,-2$ は電荷の和が零でも三乗の和が $-6$ で消えない。大域的なゲージアノマリーの条件は別であり、局所三角図の相殺だけでは保証されない。

$\delta_\omega\Gamma=\mathcal A(\omega)$ と定義し、変換の符号を $[\delta_{\omega_1},\delta_{\omega_2}]=\delta_{[\omega_1,\omega_2]}$ となるよう取る。両辺を $\Gamma$ に作用させると
$$
\delta_{\omega_1}\mathcal A(\omega_2)-\delta_{\omega_2}\mathcal A(\omega_1)
=\mathcal A([\omega_1,\omega_2])
$$
を得る。これが :term[Wess–Zumino の無矛盾条件]{en="Wess–Zumino consistency condition" yomi="うぇすずみののむむじゅんじょうけん"} である。局所な対項 $B[A]$ は $\mathcal A$ を $\delta_\omega B$ だけ変えるので、消去可能な項をこの形で判別する。以下では Chern–Weil の推移式から代表元を求める。局所解の完全な分類は、すべての局所汎関数での BRST コホモロジーを要するため :cite[weinberg-qft2] を参照する。

### 電荷三乗と重力を含む一次和

すべて左 Weyl として電荷を数えると、Abel の局所三角条件は $\sum_i q_i^3=0$。背景重力との混合アノマリーの電荷因子は $\sum_i q_i$ である。重力への結合ではフレーバーに単位行列が作用し、ゲージ生成子だけをフレーバーでトレースするためである。例えば五つの電荷 $(-9,-5,-1,7,8)$ は一次和が $-15+15=0$、三次和が $-729-125-1+343+512=0$。対になった逆電荷を持たないので、この条件の成立は Dirac 質量を全場に許すことと同じではない。実際、素朴な二場のゲージ不変質量には和が零の電荷が必要だが、ここにはその対がない。

### 降下で無矛盾な局所代表を作る

反 Hermite の行列値接続一形式を $A_{\rm a}$、曲率を $F_{\rm a}=\dd A_{\rm a}+A_{\rm a}^2$ とする。積は外積と行列積を合わせたものとし、$\delta_\alpha A_{\rm a}=D\alpha$ を使う。Chern–Weil の推移式 :ref[thm-dtop-chern-weil] に $A_t=tA_{\rm a}$ を入れると $F_t=t\,\dd A_{\rm a}+t^2 A_{\rm a}^2$、$\dd Q_5=\tr_R F_{\rm a}^3$ となる。ただし
$$
Q_5=3\int_0^1\dd t\,\tr_R(A_{\rm a} F_t^2)
=\tr_R\left(A_{\rm a}(\dd A_{\rm a})^2+
\tfrac32 A_{\rm a}^3\,\dd A_{\rm a}+\tfrac35 A_{\rm a}^5\right).
$$
例えば混合項の $3/2$ は $3\int_0^1 2t^3\dd t$ である。ここで全体の跡の中では二形式の巡回に次数符号は付かない。曲率は無限小変換で共役されるため $\delta_\alpha \tr F_{\rm a}^3=0$。その一段下の式は
$$
\delta_\alpha Q_5=\dd Q_4^{(1)}(\alpha,A_{\rm a}),\qquad
Q_4^{(1)}=\tr_R\!\left[\alpha\,\dd\!\left(A_{\rm a}\dd A_{\rm a}+\tfrac12 A_{\rm a}^3\right)\right]
$$
であり、$Q_4^{(1)}$ には完全形式を加える自由がある。変分の計算は、跡の巡回で $[A_{\rm a},\alpha]$ による同時共役分を消し、$\delta(\dd A_{\rm a})=\dd(\delta A_{\rm a})$ の $\dd\alpha$ の項を合わせて行う。具体的に $B=\dd A_{\rm a}$、$u=\dd\alpha$ と置くと
$$
\delta A_{\rm a}=[A_{\rm a},\alpha]+u,\qquad
\delta B=[B,\alpha]-A_{\rm a}u-u A_{\rm a}.
$$
同時共役を除く変分の項は $\tr u[B^2+\tfrac12(BA_{\rm a}^2-A_{\rm a}BA_{\rm a}+A_{\rm a}^2B)]$。これは $\tr u\,\dd(A_{\rm a}B+A_{\rm a}^3/2)$ だから、上の外微分が得られる。$B^2$ は二つの背景場に二微分を持ち、$A_{\rm a}^3$ による部分が非可換の整合性を満たすために加わる。

$\int Q_4^{(1)}$ は有効作用の変分として使う :term[整合アノマリー]{en="consistent anomaly" yomi="せいごうあのまりー"} の代表である。規格化には六形式 $[\tr_R\exp(iF_{\rm a}/(2\pi))]_{6}=(i/(2\pi))^3\tr_R F_{\rm a}^3/3!$ の降下を用い、左右カイラリティで符号を逆にする。四形式の局所積分の全体位相を決めるには $\Gamma$ を Euclid の $-\log Z$ とするか Minkowski の有効作用とするかを指定する。曲率だけで共変に変換する電流の発散は :term[共変アノマリー]{en="covariant anomaly" yomi="きょうへんあのまりー"} と呼ぶ。例えば $\mathcal A_{\rm cov}(\alpha)=k\int\tr(\alpha F_{\rm a}^2)$ なら $\delta_{\alpha_1}\mathcal A_{\rm cov}(\alpha_2)=k\int\tr([\alpha_1,\alpha_2] F_{\rm a}^2)$。二つの変換の差は右辺の二倍になるので、一般の非可換ではこれをそのまま有効作用の変分にできない。二種類の電流の差は局所な多項式の電流で記述できるが、その多項式が常に一つの局所作用の変分とは限らない。

### 一世代の電荷で全ての局所条件を点検する

例として群 $SU(3)\times SU(2)\times U(1)$ で、すべてを左 Weyl として表現を
$(3,2)_{1/6}$、$(\bar3,1)_{-2/3}$、$(\bar3,1)_{1/3}$、$(1,2)_{-1/2}$、$(1,1)_1$ とする。$SU(3)^3$ の係数は基本三重項の一個分で $2-1-1=0$。$SU(2)^3$ は基本二重項の対称三トレースが零のため消える。混合 $SU(3)^2U(1)$ は $T(F)[2(1/6)-2/3+1/3]=0$、$SU(2)^2U(1)$ は $T(F)[3(1/6)-1/2]=0$。単純群の生成子一つだけの混合はその無跡性で零である。Abel の三乗と混合重力は
$$
\begin{aligned}
&6(1/6)^3+3(-2/3)^3+3(1/3)^3+2(-1/2)^3+1^3=0,\\
&6(1/6)+3(-2/3)+3(1/3)+2(-1/2)+1=0 .
\end{aligned}
$$
多重度を落とすと、粒子の成分が見る電荷でなく表現の名前だけを数えることになる。さらに $SU(2)$ 基本二重項は色を含めて $3+1=4$ 個。通常のスピン構造を持つ四次元のこの表現では、奇数個の基本 Weyl 二重項で起きる符号の大域アノマリーにも抵触しない。後者の一般証明は大きな変換の経路での Dirac 固有値の流れを要するため :cite[weinberg-qft2] を参照する。

### 実表現と共役表現

ある表現 $R$ の共役の生成子は $T_{\bar R}^a=-(T_R^a)^{\mathsf T}$ である。従って三生成子のトレースは共役表現で全体の負号を取り、積の向きが逆になる。反交換子により二順序を足しているため反転後も同じ対称トレースである。このため $R\oplus\bar R$ は左手場だけで書いても局所ゲージアノマリーが消える。$R$ と $\bar R$ が同型ならその対称トレースは自身の負に等しい。例えば $SU(2)$ の基本表現では $\{T^b,T^c\}=\delta^{bc}1/2$、$\tr T^a=0$ なので局所三角項は零になる。これで検査したのは Lie 代数の無限小変換に対する整合性である。

### Abel の整合係数と共変係数の比

非可換の降下の規格化を、同じ電荷の Abel 成分で確かめる。一組の Abel 背景に絞ると $A_{\rm a}^2=0$、$Q_4^{(1)}=\tr(\alpha\,\dd A_{\rm a}\,\dd A_{\rm a})$。この積分は有効作用から三点頂点を三回微分した Bose 対称な配置に対応する。二つの源の微分で得る Ward 式には三つの電流が同等に扱われる。例えば三角図の局所項を、二つのベクトル Ward 式を保った軸側へ全て送る場合と、単一カイラル電流の三挿入の各頂点へ等しく送る場合を比較すると、一頂点あたりは三分の一になる。降下でも $\tr F_{\rm a}^3/3!$ を一段降ろした整合項は、発散させる一つの電流の共変なトレースに比べこの因子を持つ。

より直接には、整合する Abel の三形式の電流を $J_{\rm cons}$、その発散を $\dd J_{\rm cons}=\kappa\,\dd A_{\rm a}\,\dd A_{\rm a}$ と表す。共変な局所電流との差は $J_{\rm cov}-J_{\rm cons}=2\kappa A_{\rm a}\dd A_{\rm a}$ となる規約で、外微分は $2\kappa(\dd A_{\rm a})^2$。この差が有効作用の Abel 対項の微分にならないことは、一つの場だけの $A_{\rm a}^2\dd A_{\rm a}$ が反対称性で零になる例からも見える。絶対値の比較を四次元の :eqref[eq-qft-abj] に合わせるなら、右 Weyl の共変な粒子数は Dirac 軸電流の半分、ゲージ電流にはもう一つの $q$ が入る。従って同じ向きで整合するゲージ発散の $|\epsilon F F|$ 当たりの係数は $|q|^3/(96\pi^2)$、共変な係数は $|q|^3/(32\pi^2)$ である。左右で符号が逆になるので、局所消去条件はどちらの係数から読んでも同じになる。

### 反対表現と大域対称性に残る異常

同じ群の表現の直和では、カイラルな対称トレースは和になる。例えば $SU(3)$ の基本と共役を一組にすれば三次のゲージ異常は打ち消せるが、同じ二場に独立な大域電荷を割り当てれば、その大域電流と二つのゲージ挿入の混合トレースは一般に消えない。二つのカイラルな場に同じ大域電荷を与えると、二つのゲージ生成子の跡は共役表現で同符号なので足し合わさる。このように、ゲージの一貫性を保った後に大域軸電流の非保存が残る例を作れる。

群を背景にだけ結合して大域対称性の Ward 式に異常がある場合、その局所対項で消せない部分は有効理論の低エネルギーでも一致する。背景源の変数変換の位相は、重い自由度を正しく積分しても同じ生成汎関数に残るためである。この条件を :term[アノマリー整合]{en="anomaly matching" yomi="あのまりーせいごう"} という。例えば真空で大域対称性が破れ、適切な角場があれば、上の $\pi F\widetilde F$ のように角場が変換する局所項で位相を再現できる。対称性が破れず、自明で一意なギャップだけを持つ低エネルギー理論でこの位相が再現できるかは、許す局所対項や大域的な応答を含めて判定する。これだけで強結合相の励起の内容が一意に決まるとは限らない。

## トレースアノマリー

質量ゼロの理論では古典の尺度変換に保存流 $x_\nu\Theta^{\mu\nu}$ が対応する。その発散は $\Theta^\mu{}_\mu$ である。くりこみで尺度 $\mu$ が入り、尺度変換に伴う無次元結合の変化は $\beta_g$ で与えられる。従って平坦空間のゲージ不変相関では運動方程式と接触項を除き、トレースは $\beta_g\,\partial\mathcal L/\partial g$ に対応する。これを :term[トレースアノマリー]{en="trace anomaly" yomi="とれーすあのまりー"} という。接続を $\mathcal A=gA$、曲率を $\mathcal F=gF$ として固定すれば $\mathcal L=-\mathcal F^2/(4g^2)$、$\partial_g\mathcal L=\mathcal F^2/(2g^3)$ なので
$$
\Theta^\mu{}_\mu=\frac{\beta_g}{2g}F_{\rho\sigma}^a F^{a\rho\sigma}
\label{eq-qft-trace-anomaly}
$$
となる。曲率や演算子の積は同じ規約でくりこむ。質量があれば $\mu\,\dd m/\dd\mu=-\gamma_m m$ の規約で $(1+\gamma_m)m\bar\psi\psi$ を加える。これは質量の古典次元一と発展による追加分を合わせた係数である。例えば純 $SU(3)$ の先頭 $\beta_g=-11g^3/(16\pi^2)$ ならゲージ項の係数は $-11g^2/(32\pi^2)$ である。

### 次元正則化で古典トレースと極を掛ける

曲率を $\mathcal F=gF$ とする $d=4-2\epsilon$ の規約で、裸の純ゲージのトレースは $(d-4)\mathcal F_0^2/(4g_0^2)$。これは $T_{\mu\nu}=-\mathcal F_{\mu\rho}\mathcal F_\nu{}^\rho/g_0^2+\eta_{\mu\nu}\mathcal F^2/(4g_0^2)$ を縮約して得る。背景場の規約で $1/g_0^2=\mu^{-2\epsilon}[1/g^2+b_0/(16\pi^2\bar\epsilon)+\cdots]$ なので $d-4$ を掛けた有限な残りは $-b_0\mathcal F^2/(32\pi^2)$。これは :eqref[eq-qft-ym-beta] の四次元係数を :eqref[eq-qft-trace-anomaly] へ入れた式と一致する。ここで曲率挿入はゲージ不変な背景に対応づけた。独立な局所演算子として使う際には、場の方程式に比例する挿入や接触項を別に考える。

### 尺度 Ward 恒等式と改善項

平坦空間で対称なエネルギー運動量テンソル $\Theta^{\mu\nu}$ を取り、$\partial_\mu\Theta^{\mu\nu}=0$ とする。$\partial_\mu(x_\nu\Theta^{\mu\nu})=\Theta^\mu{}_\mu$ なので、尺度の Ward 恒等式ではこの局所挿入を調べる。ただし $\Theta$ は保存を変えない改善項を加えられる。実スカラーの一例として
$$
\Delta\Theta_{\mu\nu}=\frac16(\eta_{\mu\nu}\Box-\partial_\mu\partial_\nu)\phi^2
$$
を加える。発散は微分の可換性で零、四次元のトレースは $\Box\phi^2/2=(\partial\phi)^2+\phi\Box\phi$。正準トレースは微分項 $-(\partial\phi)^2$ を含むので、改善によって自由質量ゼロの殻上では零になる。$V=\lambda\phi^4/4!$ なら改善後の古典トレースは $\phi\Box\phi+4V$。運動方程式 $\Box\phi=-V'$ と斉次性 $\phi V'=4V$ で相殺する。これで、古典の改善可能な全微分項と量子の $\beta$ 項を区別できる。

固定裸パラメータで RG を動かすと、相関の結合微分は作用の $\int \partial_g\mathcal L$ を挿入する。物理的な尺度変換では、座標を同じ比率で縮めた分を $\mu$ の変更と結合の変化で補う。局所化した Ward 恒等式でこれを $\beta_g\,\partial_g\mathcal L$ と同定するが、場の再規格化に伴う運動方程式の挿入や演算子の衝突による接触項は残り得る。質量項では、次元一の古典スケーリングと $\mu\dd m/\dd\mu=-\gamma_m m$ を合わせた $(1+\gamma_m)m$ が掛かる。例えば自由で質量 $m$ の Dirac 場なら $\beta=0,\gamma_m=0$ でも $\Theta^\mu{}_\mu=m\bar\psi\psi$ が残る。従ってトレースの全てを量子的なアノマリーと数えない。

### 曲がった背景で許される別の局所項

四次元で背景計量を導入し局所 Weyl 変換 $\delta g_{\mu\nu}=2\omega(x)g_{\mu\nu}$ を調べると、平坦空間では零になる曲率二乗も挿入に入る。固定点で、パリティ偶かつ境界を除いた代表は
$$
\langle\Theta^\mu{}_\mu\rangle_g=
\frac{1}{16\pi^2}
\left(c\,C_{\mu\nu\rho\sigma}C^{\mu\nu\rho\sigma}
-a E_4+b\Box R\right),\qquad
E_4=R_{\mu\nu\rho\sigma}^2-4R_{\mu\nu}^2+R^2 .
$$
$C$ は Weyl 曲率、$R$ はスカラー曲率である。曲率二乗スカラーの基底は $C^2,E_4,R^2$ と取れる。局所 Weyl 変換は互いに可換なので、アノマリーの二重変分も対称でなければならない。$\delta_\omega\sqrt g=4\omega\sqrt g$、$\delta_\omega R=-2\omega R-6\Box\omega$ から、$\delta_{\omega_2}\int\sqrt g\,\omega_1 R^2=-12\int\sqrt g\,\omega_1 R\Box\omega_2$。$R$ が変化する背景では、部分積分しても一般には $\omega_1,\omega_2$ に対称にならない。従って独立な $R^2$ はトレースに入れられない。

一方 $\sqrt g C^2$ は Weyl 不変である。また $E_4=C^2-2R_{\mu\nu}^2+2R^2/3$ に $\delta_\omega R_{\mu\nu}=-2\nabla_\mu\nabla_\nu\omega-g_{\mu\nu}\Box\omega$ を入れると $\delta_\omega(\sqrt g E_4)=8\sqrt g (R^{\mu\nu}-g^{\mu\nu}R/2)\nabla_\mu\nabla_\nu\omega$。Bianchi 恒等式により括弧内の発散は零である。よって二つの Weyl パラメータを掛けた変分は一度部分積分すると対称になる。最後の $\Box R$ は $\delta_\omega\int\sqrt g R^2=-12\int\sqrt g\,\omega\Box R$ から生じ、局所対項の変分なので整合する。有限な $R^2$ 項で $b$ が動くが、$a,c$ の値は個々の理論で決まる。例えば最大対称で $R$ が一定なら $C=0,\Box R=0$、$R_{\mu\nu\rho\sigma}^2=R^2/6,R_{\mu\nu}^2=R^2/4$ なので $E_4=R^2/6$。平坦で $\beta=0$ でも、この背景では $a$ を含むトレースが残り得る。

### 重い電荷を積分したときのトレース

質量 $M$ の重い Dirac フェルミオン一種を持つ可換ゲージ理論を考え、低運動量で重い場を除いた有効理論に照合する。背景の接続 $\mathcal A=qA$ を同じ規約で固定すると、フル理論には $\beta_H F^2/(2q)$ と $M\bar\psi\psi$ の両方がある。低エネルギーで他の軽い荷電場がなければ、次元四の Maxwell 理論の $\beta_L$ は零である。そこで $\Theta^\mu{}_\mu$ を低運動量の二光子行列要素で一致させるには
$$
M\bar\psi\psi\ \longrightarrow\
\frac{\beta_L-\beta_H}{2q}\,F^2
=-\frac{q^2}{24\pi^2}F^2
$$
が必要となる。最後は一ループ $\beta_H=q^3/(12\pi^2)$ である。高い次元の項は外部運動量の逆質量展開から生じ、質量の異常次元の補正はこの表示で高次である。この結果は、質量演算子を単独で有効理論から捨てると、重い場のベータ関数への寄与が物理トレースに誤って残ることを示す。

照合は同じ背景の閾値対数を質量で微分しても確かめられる。Euclid の一ループ作用の非定数の低運動量部分は、$\mathcal F=qF$ で $\delta(1/q^2)\mathcal F^2/4$。Dirac の $b_0=-4/3$ から $\delta(1/q^2)=-(1/(12\pi^2))\log(M^2/\mu^2)$。従って $M\partial_M$ は $-\mathcal F^2/(24\pi^2)$ を与える。Euclid で $M\partial_M S_E$ は正の符号で $M\bar\psi\psi$ の挿入に対応し、ここでは同じ演算子の規格化で Minkowski の表示と比較した。対数の有限定数は質量微分で消える。発散する軸電流の重い場による局所項と、トレースの重い場による質量挿入は、どちらも測定する背景と局所演算子の整合が重要である。

### トレースを場の方程式と比較する

自由実スカラーの質量項を $V=m^2\phi^2/2$ に取ると、改善後の古典トレース $\phi\Box\phi+4V$ は、$\Box\phi=-m^2\phi$ により $m^2\phi^2$ になる。これは質量二乗の結合を使えば、古典次元が二なので $-2m^2\,\partial\mathcal L/\partial m^2=+m^2\phi^2$ である。これに対し無次元結合の量子項は $\beta_\lambda\,\partial_\lambda\mathcal L$。$\mathcal L_{\rm int}=-\lambda\phi^4/4!$ の規約でこれは $-\beta_\lambda[\phi^4]/4!$ である。例えば同じ量を $V$ で書いてから $\beta_\lambda\partial_\lambda V$ とすると、Minkowski の作用との符号を取り違える。挿入が運動方程式に比例する部分は、分離した相関では Schwinger–Dyson の変数変換で消えるものの、接触点の恒等式では現れる。これを除いた式が、物理相関の尺度依存と対応する。

次元正則化でもこの負号を確認できる。$d=4-2\epsilon$ で改善した裸の質量ゼロのトレースは、運動方程式の項を除き $2\epsilon\lambda_0\phi_0^4/4!$。くりこんだ局所挿入を結合微分で定めると、その $\phi_0^4$ 部分には $\partial_\lambda\lambda_0$ が掛かる。固定裸量で $0=2\epsilon\lambda_0+\beta^{(d)}_\lambda\partial_\lambda\lambda_0$ が成り立つので、裸のトレースをこの挿入へ移した係数は $-\beta^{(d)}_\lambda$。例えば $\lambda_0=\mu^{2\epsilon}(\lambda+b\lambda^2/\epsilon+\cdots)$ なら $2\epsilon\lambda_0/(\partial_\lambda\lambda_0)=2\epsilon\lambda-2b\lambda^2+\cdots$ となる。$\lambda_0$ だけを展開して裸の $\phi_0^4$ を有限として扱うと、有限項の符号を誤る。曲率二乗にも結合微分で定めた演算子を使う理由は同じである。

## 演習問題

::::exercise{#exr-qft-anomalies-1 level=basic}
質量ゼロで Dirac の古典軸性発散を求めよ。

:::solution
局所軸回転で運動項の変化は $-j_5^\mu\partial_\mu\alpha$、質量項の変化は零である。部分積分と運動方程式で $\partial_\mu j_5^\mu=0$ を得る。
:::
::::

::::exercise{#exr-qft-anomalies-2 level=basic}
純電場背景で :eqref[eq-qft-abj] の電磁項を求めよ。

:::solution
$\vb B=0$ なので $F\widetilde F=-4\vb E\cdot\vb B=0$、従って零である。
:::
::::

::::exercise{#exr-qft-anomalies-3 level=basic}
質量を $M$ 倍した三分母 Euclid 積分の重い極限で $M^2 I_3(M^2)$ を求めよ。

:::solution
:eqref[eq-qft-loop-master] で $d=4,n=3$ を取ると $I_3=1/(32\pi^2 M^2)$、積は $1/(32\pi^2)$ となる。
:::
::::

::::exercise{#exr-qft-anomalies-4 level=basic}
左 Weyl の電荷 $2,-2$ で Abel の局所ゲージアノマリーを求めよ。

:::solution
三乗和は $8-8=0$ で相殺する。左右を分けず左だけで符号つき表現を数えた結果である。
:::
::::

::::exercise{#exr-qft-anomalies-5 level=advanced}
$F_{01}=E$、$F_{23}=-B$ だけが独立な場合、$\epsilon F F$ を計算せよ。

:::solution
$(01)(23)$ の代表項は $-EB$。各組の反転は $\epsilon$ と $F$ の二負号を作り、二組の交換は偶置換である。全八項から $-8EB$ となる。
:::
::::

::::exercise{#exr-qft-anomalies-6 level=advanced}
Euclid 熱核の非零固有モードの寄与が積分で消えることを示せ。

:::solution
$\{\gamma_*,\not D_E\}=0$ なので非零固有値 $i\lambda$ の固有ベクトルは $\gamma_*$ で $-i\lambda$ の空間へ移る。直交性でその対角要素は零である。零固有空間だけで $\gamma_*$ の固有値 $\pm1$ を数えればよい。
:::
::::

::::exercise{#exr-qft-anomalies-7 level=advanced}
Schwinger 模型の電場の静止励起の周波数を求めよ。

:::solution
$\partial_xE=0$ なら $\partial_t^2E+(q^2/\pi)E=0$。正の周波数は $|q|/\sqrt\pi$ である。
:::
::::

::::exercise{#exr-qft-anomalies-8 level=advanced}
一定なベータ関数係数 $b_0$ を使い、$\beta_g=-b_0g^3/(16\pi^2)$ によるトレース項を求めよ。

:::solution
:eqref[eq-qft-trace-anomaly] の係数は $\beta_g/(2g)=-b_0g^2/(32\pi^2)$ なので、これに $F^2$ を掛ける。
:::
::::

::::exercise{#exr-qft-anomalies-9 level=advanced}
円周 $L$ の質量ゼロの Dirac 場に一様な電場を期間 $T$ だけ作用させ、$q\int_0^T E(t)\dd t=2\pi n/L$、$n\in\Z$ とする。初期と最終のスペクトルを比べ、軸電荷とベクトル電荷の変化を求めよ。初期・最終で零モードが境界にちょうど載る曖昧さは、例えば空間反周期境界で避けてよい。

:::solution
運動学的運動量は隣接する準位の間隔 $2\pi/L$ の $n$ 倍だけ移る。$n>0$ なら右向き枝で負から正へ $n$ 準位、左向き枝で正から負へ $n$ 準位が移る。充填した Dirac 海と比べ、前者は右粒子、後者は左の穴を与える。従ってベクトル電荷の差は $n-n=0$、右数から左数を引いた軸電荷の差は $n+n=2n$。$n<0$ なら交差の向きが逆なので式は同じである。場の方程式の形でも $\Delta Q_5=L(q/\pi)\int_0^T E\,\dd t=2n$ となる。ベクトル電荷の保存は二つの枝を両方数えて得られる。
:::
::::

::::exercise{#exr-qft-anomalies-10 level=advanced}
二次元の電荷 $q,2q$ の質量ゼロ Dirac 場に対し、ゲージ背景で保存する軸カレントを一つ求めよ。正準スカラーの局所表示で有質量な方向と直交方向を取り、電場の遮蔽長を求めよ。

:::solution
発散は各々 $qE/\pi,2qE/\pi$ なので $2j_{5,1}-j_{5,2}$ が保存する。正準な場の結合は $E(q\varphi_1+2q\varphi_2)/\sqrt\pi$ に比例するから、有質量な方向は $(\varphi_1+2\varphi_2)/\sqrt5$、直交方向は $(2\varphi_1-\varphi_2)/\sqrt5$。電場に結合する組の質量二乗は $5q^2/\pi$、遮蔽長は $\sqrt\pi/(\sqrt5|q|)$。四次元の電磁軸性アノマリーで係数に電荷の二乗が現れることと区別する。この問題の保存する方向は、二次元での一次の電荷の係数から求めた。一方で力学的な質量二乗には、Maxwell 方程式でもう一度電荷が掛かるため二乗が現れる。
:::
::::

::::exercise{#exr-qft-anomalies-11 level=advanced}
$SU(2)\times U(1)$ の左 Weyl 場に二重項を二つ、電荷 $+1,-1$ で入れる。さらに一重項の電荷を $2,-1,-1$ とする。局所ゲージ条件、混合重力条件、基本二重項の大域条件をそれぞれ調べよ。

:::solution
$SU(2)^3$ の対称トレースは基本で零。$SU(2)^2U(1)$ は $T(F)(1-1)=0$、$SU(2)U(1)^2$ は無跡性で零である。一重項はこれらに寄与しない。$U(1)^3$ では二重項の二成分を足して $2(1^3+(-1)^3)=0$、一重項から $2^3-1-1=6$。従って局所ゲージ異常は残る。一方混合重力の一次和は $2(1-1)+2-1-1=0$。基本二重項は偶数の二個で、その大域的な符号の異常は相殺する。例えば一重項をすべて逆電荷の組で置き換えれば $U(1)^3$ は消えるが、奇数個の二重項を加えた場合の大域的条件は別途確認する必要がある。
:::
::::

::::exercise{#exr-qft-anomalies-12 level=advanced}
Euclid の四次元トーラスに、電荷 $q$ が見る $12,34$ 平面のフラックスを $2\pi(2),2\pi(-3)$ とする。局所熱核の積分、零モードのカイラリティ差、調整質量の位相の効果を求めよ。

:::solution
定数背景の局所密度を積分すると $q^2 B_1B_2L_1L_2L_3L_4/(4\pi^2)=(2)(-3)=-6$。従って $n_+-n_-=-6$。二つの Landau 問題の最下準位は各々二個と三個であり、反対向きの磁場を選んだのでカイラリティの積は負。従ってこの定数背景では六つの負カイラリティの零モードがある。$\det(\not D_E+m e^{i\phi\gamma_*})$ はこれらから $m^6e^{-6i\phi}$ を受ける。非零固有値の対では $\gamma_*$ の作用で二次行列にでき、その行列式は $\lambda^2+m^2$ なので質量位相は相殺する。従って零モードの位相は、軸回転で曲率一項の位相が移ることの有限体積の例になる。
:::
::::

## 確認テスト

::::quiz
Dirac の三角図でベクトル Ward 恒等式を保つ差引きが残すものを選べ。

- [x] 軸性発散の局所項
- [ ] ベクトル電荷の非保存
- [ ] 通常の質量項の消失

:::explanation
ベクトル頂点の縮約はシフト可能な有限積分の差になる。軸性頂点では調整場の質量挿入があり、その分子の質量と積分の逆質量が相殺して局所項が残る。
:::
::::

::::quiz{answer="-6" tol="1e-6"}
左 Weyl 電荷 $1,1,-2$ の三乗和を答えよ。

:::explanation
$1^3+1^3+(-2)^3=1+1-8=-6$ である。
:::
::::

::::quiz
二次元の Schwinger 模型で新しい極が現れるゲージ不変演算子を選べ。

- [x] 電場 $E$
- [ ] ゲージ変換のパラメータ
- [ ] 反ゴーストの場だけ

:::explanation
Maxwell 方程式と軸性アノマリーを合わせると電場が $(\Box+q^2/\pi)E=0$ を満たす。従って電場の励起にその質量が現れる。
:::
::::

::::quiz{answer="4" tol="1e-6"}
$\epsilon F F=2F\widetilde F$ の規約で、$-F\widetilde F= c\vb E\cdot\vb B$ の $c$ を答えよ。

:::explanation
平行な電場と磁場を選ぶと $\epsilon F F=-8 EB$ である。従って $-F\widetilde F=4EB$。任意の向きでも反対称な縮約は同じ内積を与える。
:::
::::

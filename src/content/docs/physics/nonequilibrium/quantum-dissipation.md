---
title: 量子散逸系
description: 調和浴を消去し、量子の摩擦、雑音、マスター方程式およびトンネルの抑制を導く。
notation:
  - symbol: 'J(\omega)'
    meaning: 正周波数での浴のスペクトル密度
  - symbol: '\gamma(t)'
    meaning: 摩擦の記憶核
  - symbol: '\xi(t)'
    meaning: 浴の初期演算子による力
  - symbol: '\rho_S'
    meaning: 系の縮約密度演算子
  - symbol: '\beta'
    meaning: 逆温度
---

## この章の目標

- 調和浴と反項を含む Hamiltonian を解く。
- 演算子の力の相関を対称化と非対称化に分ける。
- Born–Markov と永年近似が行う操作を導く。
- 虚時間の散逸作用から半古典的なトンネルを評価する。

## Caldeira–Leggett 模型

### 浴の座標をずらして結合する

:term[Caldeira–Leggett 模型]{en="Caldeira–Leggett model" yomi="かるでいられげっともけい"} は、質量 $m$ の座標 $x$ と独立な調和振動子の浴を

$$
H=\frac{p^2}{2m}+U(x)+\sum_j\left[\frac{p_j^2}{2m_j}+
\frac{m_j\omega_j^2}{2}\left(q_j-\frac{c_jx}{m_j\omega_j^2}\right)^2\right]
\label{eq-neq-caldeira}
$$

で結ぶ。:eqref[eq-neq-caldeira] の共役な各対には $[x,p]=i\hbar$、$[q_j,p_k]=i\hbar\delta_{jk}$ を課し、異なる対の交換子は零である。平方の中には相互作用 $-x\sum_j c_jq_j$ と反項 $x^2\sum_j c_j^2/(2m_j\omega_j^2)$ がある。固定した $x$ で浴座標を平衡位置に移しても浴の振動数は変わらない。反項を落とすと浴による静的なポテンシャル変化 $-x^2\sum_j c_j^2/(2m_j\omega_j^2)$ が残る。

浴のスペクトル密度を $J(\omega)=(\pi/2)\sum_j c_j^2/(m_j\omega_j)\delta(\omega-\omega_j)$、$\omega>0$ と定める。連続浴では $J$ の低周波の形と高周波の切断が別々の物理を表す。$J\simeq\gamma_0\omega$ を Ohmic と呼び、$\gamma_0$ は $\mathrm{kg\,s^{-1}}$ の次元をもつ。

:::example[単一の浴モード]
固定した粒子の座標を二倍にすると浴平衡位置のずれも二倍になる。周波数はそのままである。反項なしなら系のポテンシャルの曲率が $c_j^2/(m_j\omega_j^2)$ だけ減る。
:::

:::example[Drude 型の切断]
$J(\omega)=\gamma_0\omega\Omega^2/(\omega^2+\Omega^2)$ は $\omega\ll\Omega$ で線形であり、$\Omega$ は浴の速い尺度を指定する。反項の係数に使う $(2/\pi)\int_0^\infty J(\omega)\dd\omega/\omega$ は $\gamma_0\Omega$ と有限である。
:::

### 反項が静的な安定性を保つ

浴の一振動子を展開した結合は $-cxQ+c^2x^2/(2m_b\omega_b^2)$ を含む。$x$ を固定して浴のエネルギーを最小にする座標は $Q_*=cx/(m_b\omega_b^2)$ である。二次項まで入れると浴の最小エネルギーは $x$ に依存しない。一方、二次項を省くと最小値は $-c^2x^2/(2m_b\omega_b^2)$ である。これは調和な系の曲率を変え、十分強い結合では不安定化する。量子の場合も固定した $x$ に対する浴 Hamiltonian は座標の並進で元の振動子とユニタリ同値である。ただし $x$ 自体が演算子なら、系の運動量を含む全 Hamiltonian の同値性まで従うわけではない。

:::example[同じ散逸結合で静的曲率が変わる]
$m_b=\omega_b=1$、$c=2$、系の曲率を $\kappa=3$ とする。二次項を省いた場合は最小化後の曲率が $\kappa-c^2=-1$ になる。二次項を入れればポテンシャルは $3x^2/2+(Q-2x)^2/2$ であり、$x,Q$ の両方向で正定値である。この差は浴の近似や時間発展を用いずに見える。
:::

### 浴の分布をスペクトル密度にまとめる

浴の座標を無次元の生成消滅演算子から
$$
q_j=\sqrt{\frac{\hbar}{2m_j\omega_j}}(b_j+b_j^\dagger),\qquad
p_j=-i\sqrt{\frac{\hbar m_j\omega_j}{2}}(b_j-b_j^\dagger)
$$
と書く。$[b_j,b_k^\dagger]=\delta_{jk}$ を代入すると $[q_j,p_k]=i\hbar\delta_{jk}$ になり、自由浴の Hamiltonian は $\sum_j\hbar\omega_j(b_j^\dagger b_j+1/2)$ である。結合の力を $B=\sum_j c_jq_j$ とすれば、各モードの一量子を作る行列要素の大きさは $c_j\sqrt{\hbar/(2m_j\omega_j)}$ である。その二乗を振動数の分布として足した量は $\hbar J(\omega)/\pi$ に等しい。従ってスペクトル密度 $J$ は単なる浴の状態数ではなく、どのモードが系の座標にどれだけ結合するかを含む。例えば多くの振動子が存在しても $c_j=0$ のものは $J$ に入らない。

振動数の密度を $g(\omega)$、その範囲での結合の平均を $\overline{c^2/(m_b\omega)}$ とすれば、連続模型では $J=(\pi/2)g(\omega)\overline{c^2/(m_b\omega)}$ と書ける。低周波の $J\propto\omega^s$ に対し $s=1$ を Ohmic、$0<s<1$ を劣 Ohmic、$s>1$ を超 Ohmic と呼ぶ。抵抗の静的な近似に関わるのは $J(\omega)/\omega$ の零周波に近い挙動である。反項の強さは逆に $(2/\pi)\int J(\omega)/\omega\dd\omega$ によって全帯域の寄与を受ける。例えば $\omega$ の一次に比例する式を全帯域で使うと反項が発散するため、有限の模型の定義にはカットオフまたは繰り込んだ曲率の指定が必要である。

:::example[同じ数のモードで効果が違う浴]
同じ質量 $m_b$ の二つの振動子の振動数を $\omega_0,2\omega_0$、結合係数を同じ $c$ とする。$J$ の積分強度は一つ目が二つ目の二倍である。反項の係数は $c^2/(m_b\omega_j^2)$ に比例するので、一つ目が二つ目の四倍である。浴の二つの振動数を見つけただけでは、雑音や曲率変化のどちらに強く寄与するかは決まらない。
:::

### 調和系と単一モードを正規座標へ移す

系も $U=\kappa x^2/2$ の調和形とし、浴を一つとする。質量で規格化した位置を $X=\sqrt m x$、$Q=\sqrt{m_b}q$ とすれば、運動エネルギーは両座標で単位質量の形になる。ポテンシャルの二次行列は
$$
M=\begin{pmatrix}
\kappa/m+c^2/(m m_b\omega_b^2)&-c/\sqrt{m m_b}\\
-c/\sqrt{m m_b}&\omega_b^2
\end{pmatrix}.
$$
$\det M=(\kappa/m)\omega_b^2>0$ で跡も正なので二つの振動数平方は正である。実対称な $M$ を直交行列で対角化し、共役運動量にも同じ回転を使えば交換関係を保ったまま二つの独立な振動子になる。系の座標 $x$ は二つの正規モードの線形結合であり、一つの平均の位置応答を測っても二本の周波数が現れることがある。有限浴でエネルギーが一旦別の座標へ移る現象を、連続浴での不可逆な散逸と区別する。

:::example[共鳴した二モードが分裂する]
$m=m_b=\kappa=\omega_b=1$、$c=1$ とする。$M=\begin{pmatrix}2&-1\\-1&1\end{pmatrix}$ であり、振動数平方は $(3\pm\sqrt5)/2$ である。元の二つの振動数はともに $1$ だが、結合を入れると異なる二本が現れる。初期変位を系の座標だけに与えても両モードが励起され、変位には二つの余弦が残る。反項を外したなら $\det M=0$ になるため、同じ結合強さの比較には Hamiltonian の静的部分の規約が必要である。
:::

### 座標の平衡分布と量子の基底状態

同じ二次 Hamiltonian でも、古典の位置の平衡分布と量子の縮約状態は異なることがある。古典の正準積分では $q_j-c_jx/(m_j\omega_j^2)$ を独立な積分変数とし、$x$ の周辺密度は正確に $e^{-\beta U(x)}$ に比例する。量子では位置の平行移動で浴の平方を簡単にできるが、全 Hamiltonian の運動項も変わる。実際、一浴モードで $\alpha=c/(m_b\omega_b^2)$ とし、位置だけを $X=x,Q=q-\alpha x$ と変換する。正準一形式 $p\,\dd x+p_b\dd q=(p+\alpha p_b)\dd X+p_b\dd Q$ から、新しい運動量は $P=p+\alpha p_b$、$P_Q=p_b$ である。従って系の運動項は $(P-\alpha P_Q)^2/(2m)$ となり、位置の結合を消しても運動量の結合が残る。この線形正準変換は交換関係も保つ。

調和系の場合、質量規格化した位置の二次行列 $M$ の正規振動数を $\Omega_i$、直交固有ベクトル行列を $O$ とすれば、基底状態の位置共分散は $\hbar O\diag(\Omega_i^{-1})O^\mathsf T/2=\hbar M^{-1/2}/2$ である。各正規モードの基底状態の二乗平均を足した結果である。同じ条件の古典では共分散は $k_{\mathrm B}T M^{-1}$ となり、行列の平方根の有無に差がある。弱結合の近似で縮約系を裸の Gibbs 状態に置く場合は、結合がこれらの共分散をどの程度変えたかを確認する。

:::example[静的曲率は同じでも零点分散が変わる]
$m=m_b=\kappa=\omega_b=c=1$ の二モードで、$M=\begin{pmatrix}2&-1\\-1&1\end{pmatrix}$ である。$\det M=1$、$\Tr M=3$ なので Cayley–Hamilton の関係により $(M+I)^2=5M$ である。従って正の平方根は $(M+I)/\sqrt5$ であり、$M^{-1/2}$ の左上は $2/\sqrt5$ になる。系の基底状態での二乗位置は $\hbar/\sqrt5$ であり、裸の振動子の $\hbar/2$ と違う。一方 $M^{-1}$ の左上は $1$ であり、古典での静的な位置分散は $k_{\mathrm B}T$ のままである。
:::

### 離散浴の復帰を初期条件から読む

単一または有限個の調和浴と線形に結合した系では、全ての正規モードが離散である。実の時間発展を $x(t)=\sum_j[a_j\cos(\Omega_jt)+b_j\sin(\Omega_jt)]$ と書けるので、全振幅が時間無限大で指数的に消える結論は一般には従わない。実際、異なる振動数での長時間平均の二乗は $\frac12\sum_j(a_j^2+b_j^2)$ であり、係数の一つでも非零なら正である。交差項は三角関数を和と差の周波数に分けて積分すると消える。同じ振動数のモードがあるならまず係数を足す。連続浴に移ると積分の中で位相が相殺し、局所の観測量が減衰しうる。この二つは有限の系と連続極限をとる順序が違う。

:::example[二本の線がうなりを作る]
二本の周波数 $\Omega_1=1$、$\Omega_2=1.2$ の余弦を同じ係数で足すと $2\cos(1.1t)\cos(0.1t)$ である。包絡が零になる時刻があるが、その後に最大の振幅は戻る。観測の短い窓だけを指数で合わせると減衰と読めることがある。振動数の分離 $0.2$ を測れる時間まで記録すると、エネルギーが完全に失われる模型とは区別できる。
:::

### 実時間の二経路から浴の影響を求める

浴を自由 Hamiltonian の Gibbs 状態 $\rho_B$ として系と独立に用意する。系の密度行列を伝播させる経路積分では、順側の経路 $x_+$ と共役側の経路 $x_-$ がある。相互作用表示の自由浴の力を $B(t)=\sum_j c_j q_j(t)$、摩擦核を $\gamma(t)=\sum_j c_j^2\cos(\omega_jt)/(m_j\omega_j^2)$ とする。経路 $x$ を指定した浴の発展を
$$
U_x=\mathcal T\exp\!\left[\frac{i}{\hbar}\int_0^\tau x(t)B(t)\dd t
-\frac{i\gamma(0)}{2\hbar}\int_0^\tau x(t)^2\dd t\right]
$$
と書く。第二項は反項である。浴を跡にとった :term[影響汎関数]{en="influence functional" yomi="えいきょうはんかんすう"} は $\mathcal F[x_+,x_-]=\Tr(U_{x_+}\rho_B U_{x_-}^\dagger)$ となる。

自由浴の $q_j(t)=q_j\cos\omega_jt+p_j\sin\omega_jt/(m_j\omega_j)$ と共役交換関係から、$[B(t),B(s)]=i\hbar\dot\gamma(t-s)$ である。交換子は数だから、時間順序積を短い刻みの指数の積で表して合成すると、指数には一次積分と二項の交換子の半分だけが残る。跡の巡回性で $U_{x_-}^\dagger U_{x_+}$ を作り、再び指数を合成する。$y=x_+-x_-$、$r=(x_++x_-)/2$ とすれば、線形演算子は $(i/\hbar)\int yB$ となり、続く交換子の係数は時刻 $t>s$ で $-y(t)r(s)/\hbar^2$ となる。熱浴の一つの線形演算子 $B_y=\int yB$ は、各モードの直交位相成分の線形和である。熱状態でのその特性関数は Gaussian であり、$\log\langle e^{i B_y/\hbar}\rangle_B=-\langle B_y^2\rangle_B/(2\hbar^2)$ である。これは一振動子の熱状態の位置表示の対角カーネルが Gaussian で、その位置を位相回転させても熱状態が不変であることから従う。積分域を交換すると交換子の反対称な部分は消えるので、$N(t-s)=\langle\{B(t),B(s)\}\rangle_B/2$ を使って
$$
\log\mathcal F=
-\frac{1}{2\hbar^2}\int_0^\tau\!\!\int_0^\tau y(t)N(t-s)y(s)\dd s\dd t
-\frac{i}{\hbar}\int_0^\tau y(t)
\left[\gamma(0)r(t)+\int_0^t\dot\gamma(t-s)r(s)\dd s\right]\dd t
$$
を得る。角括弧内を部分積分すると $\gamma(t)r(0)+\int_0^t\gamma(t-s)\dot r(s)\dd s$ である。実時間の位相からも、速度との畳み込みと同じ符号の初期すべりが得られる。

:::example[離れた二経路の干渉が弱くなる]
両経路の差が一定 $y=a$ で、単一浴モードが $N(t)=N_0\cos\omega_b t$ を与えるとする。二重積分は $a^2N_0|\int_0^\tau e^{i\omega_b t}\dd t|^2=2a^2N_0(1-\cos\omega_b\tau)/\omega_b^2$ である。従って $|\mathcal F|=\exp[-a^2N_0(1-\cos\omega_b\tau)/(\hbar^2\omega_b^2)]$ である。一周期後には絶対値が $1$ に戻り、有限の単一モードの復帰を確認できる。
:::

## 量子 Langevin 方程式

### 浴振動子の運動を厳密に消去する

Heisenberg 方程式で $m_j\ddot q_j+m_j\omega_j^2q_j=c_jx(t)$ である。斉次解を $q_j^{\rm h}(t)=q_j(0)\cos\omega_jt+p_j(0)\sin\omega_jt/(m_j\omega_j)$ とすれば

$$
q_j(t)=q_j^{\rm h}(t)+\frac{c_j}{m_j\omega_j}\int_0^t\sin[\omega_j(t-s)]x(s)\dd s .
$$

積分を $s$ で部分積分すると、上端の $c_j x(t)/(m_j\omega_j^2)$、下端の $-c_j x(0)\cos\omega_jt/(m_j\omega_j^2)$、速度を余弦で畳み込んだ負の項が出る。系の方程式の反項が上端を消すので

$$
m\ddot x+U'(x)+\int_0^t\gamma(t-s)\dot x(s)\dd s=\xi(t)-\gamma(t)x(0),\quad
\gamma(t)=\frac2\pi\int_0^\infty\frac{J(\omega)}{\omega}\cos\omega t\,\dd\omega .
\label{eq-neq-quantum-langevin}
$$

ここで $\xi(t)=\sum_j c_jq_j^{\rm h}(t)$ である。最後の項は初期に浴を変位していない基準から用意したために残る。初期浴の斉次解をずれた座標から定義し直せば右辺を一つの力にまとめられる。その場合は初期状態の系と浴の相関を確認する。演算子の方程式を古典雑音で記述するにはさらに極限が要る。例えば浴の力の交換子は

$$
[\xi(t),\xi(s)]=-i\hbar\sum_j \frac{c_j^2}{m_j\omega_j}\sin[\omega_j(t-s)]
=i\hbar\,\dot\gamma(t-s)
$$

であり、一般に零でない。

:::example[一つの余弦の記憶]
単一モードは $\gamma(t)=c^2\cos(\omega t)/(m_b\omega^2)$ を与える。積分の減衰がないため、熱浴として長時間の局所摩擦に置き換えることはできない。
:::

:::example[Drude 型の記憶]
前節の $J$ を余弦変換すると $\gamma(t)=\gamma_0\Omega e^{-\Omega|t|}$ である。正の半直線での積分は $\gamma_0$ であり、系の周波数が $\Omega$ より小さいときに速度に比例する抵抗が得られる。
:::

### Drude 型カットオフと指数記憶

連続浴のスペクトルを $J(\omega)=\gamma_0\omega\Omega^2/(\Omega^2+\omega^2)$ とする。:eqref[eq-neq-quantum-langevin] の摩擦核 $\gamma(t)=(2/\pi)\int_0^\infty[J(\omega)/\omega]\cos(\omega t)\dd\omega$ は前の例の余弦変換から、$\gamma(t)=\gamma_0\Omega e^{-\Omega t}$（$t\ge0$）である。これは $e^{-\Omega|t|}$ の Fourier 変換が $2\Omega/(\Omega^2+\omega^2)$ になることから得られる。因果変換は $\widehat\gamma(\omega)=\gamma_0\Omega/(\Omega-i\omega)$ なので、低周波で摩擦係数 $\gamma_0$、高周波で実部が小さくなる。初期時刻から積分する式では、固定した非零時刻でのみ $\Omega\to\infty$ の極限を普通の摩擦として扱える。初期すべりの寄与は浴の初期状態の取り方にも依存する。

:::example[白色近似の周波数域]
$\omega=\Omega$ なら $\Re\widehat\gamma=\gamma_0/2$ である。摩擦を $\gamma_0$ と固定した近似は、この周波数で吸収に効く実部を二倍に見積もる。一方 $\omega=\Omega/10$ なら実部は $100\gamma_0/101$ であり、差は約 $1\%$ である。
:::

### 初期すべりの有無を積分表示で確かめる

調和拘束 $U(x)=\kappa x^2/2$、$\kappa>0$ を持つ系について、初期時刻以後の Laplace 変換を考える。系の座標・速度を $x_0,v_0$ とし、浴を変位させていない基準で測った自由力を $\xi$ とする。:eqref[eq-neq-quantum-langevin] から
$$
(mz^2+\kappa+z\widetilde\gamma(z))\,\widetilde x(z)
=m(zx_0+v_0)+\widetilde\xi(z)
$$
である。摩擦の積分を変換した $\widetilde\gamma(z)(z\widetilde x-x_0)$ と右辺の $-\widetilde\gamma(z)x_0$ が相殺した。一方、浴を $x_0$ に応じてシフトした座標で初期化し、力を $\xi_{\rm sh}=\xi-\gamma(t)x_0$ と書けば、右辺は $m(zx_0+v_0)+\widetilde\gamma x_0+\widetilde\xi_{\rm sh}$ になる。演算子の関係としては同じ式だが、独立な浴の分布として選ぶならどの力を中心化したかで平均の初期応答は変わる。

Drude の核 $\gamma(t)=\gamma_0\Omega e^{-\Omega t}$ で $\Omega$ を大きくすると、$-\gamma(t)x_0$ は $t>0$ で点ごとに零へ行く。しかし正の半直線の積分は $-\gamma_0x_0$ であり、一定の力積として残る。初期層の幅 $O(\Omega^{-1})$ の間に系の復元や遅い力が作る力積は小さいと仮定すれば、速度の変化に $-\gamma_0x_0/m$ が残る。従って点ごとの極限で初期項を捨てた後、真の初期値をそのまま使うと過渡を誤ることがある。

:::example[静止した座標に見える短い力積]
$m=\gamma_0=1$、$x_0=2$、$v_0=0$ とする。浴をずらさず平均零に用意した場合、初期すべりの平均の力積は $-2$ である。短い時間窓に着目すれば速度は負に動き始める。座標が固定値で最初からその条件に合った浴を準備すれば、この平均の力積は出ない。量子状態での $x_0$ と浴の相関は別途指定する必要があり、この数値例は固定した初期平均の力の計算である。
:::

### 調和応答の極を浴の時間尺度と照合する

Drude 核の Laplace 変換 $\widetilde\gamma=\gamma_0\Omega/(z+\Omega)$ を使うと、調和応答の分母を払った多項式は
$$
m z^3+m\Omega z^2+(\kappa+\gamma_0\Omega)z+\kappa\Omega
$$
である。三次の根を個別に求めなくても、決定論的な拡張方程式で安定性を判定できる。補助力 $y$ を $\dot y=-\Omega y+\gamma_0\Omega v$、$m\dot v=-\kappa x-y$、$\dot x=v$ とすればよい。二次の量 $\mathcal E=(m|v|^2+\kappa|x|^2+|y|^2/(\gamma_0\Omega))/2$ の微分は $-|y|^2/\gamma_0$ となる。交差する $v,y$ の項は相殺するためである。$m,\kappa,\gamma_0,\Omega>0$ なら、実部が零の固有モードで $y=0$ となり、式から $v=x=0$ も従う。従って全てのモードが厳密に減衰する。この証明は外力と雑音のない線形応答の安定性に関するものであり、量子浴の力を古典の白色乱数で置換するものではない。

:::example[連続浴に加えた速い尺度を極で読む]
$m=\kappa=\gamma_0=1$、$\Omega=2$ とすると多項式は $z^3+2z^2+3z+2=(z+1)(z^2+z+2)$ である。根は $-1$ と $-1/2\pm i\sqrt7/2$ である。局所摩擦 $\gamma_0$ だけでは $z^2+z+1$ の二根であり、虚部は $\sqrt3/2$ である。同じ零周波の摩擦でも、減衰する振動の周期や追加の極が変わる。
:::

### 周波数依存の反応から力の位相差を得る

安定な調和系の定常の平均座標に、力 $f_0 e^{-i\omega t}$ を加える。畳み込みの変換は $\widehat\gamma(\omega)=\int_0^\infty e^{i\omega s}\gamma(s)\dd s$ であり、平均の応答は $\chi_x f_0$、$\chi_x^{-1}=\kappa-m\omega^2-i\omega\widehat\gamma$ である。これを実・虚に分けると
$$
\chi_x^{-1}=\kappa-m\omega^2+\omega\Im\widehat\gamma-i\omega\Re\widehat\gamma .
$$
従って摩擦の位相遅れのある部分の $\Im\widehat\gamma$ は実の復元に入り、散逸を直接決めるのは $\Re\widehat\gamma$ である。時間遅れのある力を単一の定数に置けば、両方を同時には再現しない。周期の平均入力は $\omega|f_0|^2\Im\chi_x/2$ であり、$\Re\widehat\gamma\ge0$ なら非負である。これは応答に関する計算であり、雑音の分布が古典か量子かとは別に成立する。

:::example[ピークの近くで実の補正が残る]
Drude 核で $\gamma_0=\Omega=1$、$m=\kappa=1$、$\omega=1$ とすると $\widehat\gamma=(1+i)/2$、$\chi_x^{-1}=(1-i)/2$ なので $\chi_x=1+i$ になる。局所の $\widehat\gamma=1$ なら $\chi_x=i$ である。両者でこの周波数の吸収の虚部は同じだが、位置の同位相の応答は違う。一点で吸収が一致しても、摩擦の形は同じとは限らない。
:::

### 波数のない一座標で交換関係を照合する

調和系の正規座標を使い、$x=\sum_j a_j X_j$、単位質量の $X_j$ に対して $\dot X_j=P_j$ とする。質量で規格化した全空間の直交変換から $\sum_j a_j^2=1/m$ である。Heisenberg 発展により $[X_j(t),X_j(0)]=-i\hbar\sin(\Omega_jt)/\Omega_j$ であり、
$$
[x(t),x(0)]=-i\hbar\sum_j a_j^2\frac{\sin(\Omega_jt)}{\Omega_j}.
$$
$t=0$ の微分は $-i\hbar/m$ であり、$[x,p]=i\hbar$ と一致する。この関係から応答の二重交換子の総和則に $\sum_j a_j^2$ が現れる。自由な浴の力を古典乱数として扱う近似は、このような交換子を直接再現しない。調和平衡の対称化した二次の相関だけを再現する古典的な確率模型と、演算子の時間発展を記述する模型を区別する。

:::example[正規モードの重みが保存する値]
質量 $m=2$ で二本の正規モードの位置への重み $a_j^2$ が $1/8$ と $3/8$ なら、周波数がどうであっても初期の交換子の微分は $-i\hbar/2$ である。低い周波数に属する線は $\sin(\Omega t)/\Omega$ の振幅が大きくなるが、短時間の微分では周波数が相殺される。測定の帯域で片方の線だけを見た場合、総和則から得る質量を過大に見積もる。
:::

## 量子マスター方程式と Born–Markov 近似

### 密度演算子を二次まで展開する

[開放量子系](/physics/advanced-quantum-mechanics/open-systems/) の記法に沿い、相互作用表示で $V(t)=g A(t)\otimes B(t)$ とする。$A,B$ を Hermite とし、$\rho_B$ は自由浴の定常状態で $\operatorname{Tr}(\rho_B B)=0$ を満たすと仮定する。初期状態を $\rho_S(0)\otimes\rho_B$ とする。$\dot\rho_I=-i[V(t),\rho_I]/\hbar$ を一度積分し、右辺に代入して浴を跡に取れば

$$
\dot\rho_S(t)=-\hbar^{-2}\int_0^t
\operatorname{Tr}_B[V(t),[V(s),\rho_I(s)]]\dd s
$$

である。ここで $\rho_S(t)$ も相互作用表示である。Born 近似では積分内を $\rho_I(s)\simeq\rho_S(s)\otimes\rho_B$ と置く。Markov 近似では浴相関の時間 $\tau_B$ が系の緩和時間より小さいとして、遅い $\rho_S(t-\tau)$ を $\rho_S(t)$ に置き、$\tau$ の積分上限を $\infty$ に延ばす。

$[H_S,A_\omega]=-\hbar\omega A_\omega$ と分解し、$A(t)=\sum_\omega e^{-i\omega t}A_\omega$ とする。Born–Markov の積分から各組に $e^{i(\omega'-\omega)t}$ が現れる。異なる周波数差の周期より長く、緩和より短い時間窓で平均する永年近似により、同じ周波数の組だけを残す。半直線変換 $G_\omega=g^2\hbar^{-2}\int_0^\infty e^{i\omega t}\langle B(t)B(0)\rangle\dd t$ の実部と虚部に二重交換子を分けると、Schrödinger 表示で

$$
\dot\rho_S=-\frac{i}{\hbar}[H_S+H_{\rm LS},\rho_S]
+\sum_\omega 2\Re G_\omega
\left(A_\omega\rho_S A_\omega^\dagger-\tfrac12\{A_\omega^\dagger A_\omega,\rho_S\}\right),\qquad
H_{\rm LS}=\hbar\sum_\omega(\Im G_\omega) A_\omega^\dagger A_\omega .
$$

ここで $\{C,D\}=CD+DC$ である。$2\Re G_\omega$ は浴の両側スペクトルに比例する。エネルギー固有基底で相関を展開すると非負な重みの和になる。縮退する周波数は同じ $A_\omega$ にまとめる。Born–Markov だけの非永年の式には、同じ完全正値性の保証がない。

:::example[二準位の占有率]
下降率を $d$、上昇率を $u$ とすれば $\dot p_e=u-(u+d)p_e$ であり、定常値は $u/(u+d)$ である。浴が平衡なら次節の詳細つり合いから $u/d=e^{-\beta\hbar\omega}$ となる。
:::

:::example[エネルギーを変えない散逸]
$A_0=\operatorname{diag}(1,-1)$ だけがあるとき、散逸部は $r(A_0\rho A_0-\rho)$ である。非対角成分は $e^{-2rt}$ で消え、対角成分は変化しない。
:::

### 二準位の緩和と位相緩和を区別する

基底 $\ket g$ と励起状態 $\ket e$ を結ぶ Lindblad 飛躍を $L_-=\sqrt{\gamma_\downarrow}\ketbra{g}{e}$、$L_+=\sqrt{\gamma_\uparrow}\ketbra{e}{g}$ とする。各散逸項は $L\rho L^\dagger-\{L^\dagger L,\rho\}/2$ である。対角成分では $\dot\rho_{ee}=\gamma_\uparrow(1-\rho_{ee})-\gamma_\downarrow\rho_{ee}$、非対角成分では散逸によって $\dot\rho_{eg}=-(\gamma_\uparrow+\gamma_\downarrow)\rho_{eg}/2$ となる。各ジャンプ項は対角成分にしか直接入り込まず、非対角の減衰は反交換子から出る。さらに $L_z=\sqrt{\gamma_\phi/2}(\ketbra{e}{e}-\ketbra{g}{g})$ を加えると、対角成分は不変、非対角に $-\gamma_\phi\rho_{eg}$ が加わる。

:::example[熱平衡と緩和の二つの時間]
$\gamma_\uparrow=1$、$\gamma_\downarrow=3$、$\gamma_\phi=2$ を毎秒の率とする。励起状態の定常確率は $1/4$、占有確率の緩和率は $4\,\mathrm{s^{-1}}$ である。位相を含む非対角項の減衰率は $2+2=4\,\mathrm{s^{-1}}$ になる。純粋位相緩和を取り去れば占有確率を変えず、非対角の減衰率だけが $2\,\mathrm{s^{-1}}$ になる。
:::

### 二次近似で不可逆な項が生まれる場所

相互作用表示で $V_I(t)=A(t)\otimes B(t)$、$\langle B\rangle_B=0$ とし、結合強さは $B$ に含める。全状態は $\dot\rho_{\rm tot,I}=-(i/\hbar)[V_I,\rho_{\rm tot,I}]$ である。積分式を一度代入し、浴の跡を取れば
$$
\dot\rho_I(t)=-\hbar^{-2}\int_0^t\Tr_B
[V_I(t),[V_I(t-s),\rho_{\rm tot,I}(t-s)]]\,\dd s
$$
となる。最初の一次項が消えるのは、初期状態が積で浴の一階平均が零だからである。Born 近似は積分の中の全状態を $\rho_I(t-s)\otimes\rho_B$ と置く操作である。$[H_B,\rho_B]=0$ なら浴相関は時差 $s$ だけに依存する。相関が消える時間内で $\rho_I$ の変化が小さいとき、相互作用表示で $\rho_I(t-s)$ を $\rho_I(t)$ に置き換え、十分遅い時刻では積分を無限へ延長する。これらが Markov 近似である。Schrödinger 表示では状態に速い自由回転が入るため、そのまま同じ置き換えを行うと遷移の周波数を失う。

$A(t)=\sum_\omega e^{-i\omega t}A_\omega$ と展開し、$A_{-\omega}=A_\omega^\dagger$ とする。浴の相関を $C(s)=\Tr_B[\rho_B B(s)B(0)]$ とすれば、二重交換子の四項の浴の跡は順序に応じて $C(s)$ または $C(-s)$ になる。Markov 近似を施した :term[Redfield 方程式]{en="Redfield equation" yomi="れっどふぃーるどほうていしき"} は
$$
\dot\rho_I(t)=-\hbar^{-2}\int_0^\infty
\bigl(C(s)[A(t),A(t-s)\rho_I(t)]
+C(-s)[\rho_I(t)A(t-s),A(t)]\bigr)\dd s .
$$
第一の交換子の二項は $B(t)$ が $B(t-s)$ より左にある順序でまとめた。第二では逆順であり、$C(-s)=C(s)^*$ なので二つは互いに Hermite 共役である。周波数展開を代入した各項には $\int_0^\infty C(s)e^{i\omega s}\dd s$ またはその共役が入る。異なる周波数の $\omega,\omega'$ の積には相互作用表示で振動因子 $e^{i(\omega'-\omega)t}$ が残る。時間窓 $\Delta t$ で平均すると規格化した絶対値は
$$
\left|\frac{\sin[(\omega'-\omega)\Delta t/2]}{(\omega'-\omega)\Delta t/2}\right|
$$
である。従ってこれらを落とす永年近似には、浴の相関時間より長く、異なる周波数の逆差よりも長く、状態の散逸時間より短い窓が必要である。縮退した周波数は平均で消えないため同じブロックに残す。

:::example[近い二つの遷移は分離できない]
相関時間を $10^{-12}\,\mathrm s$、散逸時間を $10^{-7}\,\mathrm s$ とする。二つの角周波数の差が $10^9\,\mathrm{s^{-1}}$ なら、例えば $\Delta t=10^{-8}\,\mathrm s$ では差と窓の積が $10$ になる。差が $10^6\,\mathrm{s^{-1}}$ なら同じ窓で積は $0.01$ であり、振動はほぼ減らない。浴を短記憶として扱う条件が満たされても、すべての周波数を別の散逸チャネルとして分けられるわけではない。
:::

### 永年近似後の行列の正値性を確かめる

複数の Hermite な結合 $A_\alpha\otimes B_\alpha$ を許す。周波数 $\omega$ での率行列は $\Gamma_{\alpha\beta}(\omega)=\hbar^{-2}\int e^{i\omega t}\langle B_\alpha(t)B_\beta(0)\rangle\dd t$ である。浴の固有状態で展開すれば、任意の複素列 $c$ に対する $c^\dagger\Gamma c$ の重みは
$$
\frac{2\pi}{\hbar^2}\sum_{n,m}p_n\left|\sum_\alpha c_\alpha (B_\alpha)_{mn}\right|^2
\delta(\omega-(E_m-E_n)/\hbar)
$$
である。従って正半定値であり、適切な連続浴で通常の関数として率が得られるとき、各周波数で固有分解して正の係数をもつ飛躍演算子に直せる。有限の系で生成器が $-i[H,\rho]/\hbar+\sum_j(L_j\rho L_j^\dagger-\{L_j^\dagger L_j,\rho\}/2)$ となった場合、短い刻みを $K_0=I-iHh/\hbar-\frac h2\sum_jL_j^\dagger L_j$、$K_j=\sqrt h\,L_j$ の作用の和で近似できる。$\sum K_j^\dagger K_j=I+O(h^2)$ であり、右からこの和の逆平方根を掛けて厳密に規格化すると、一次の生成器は変わらない。刻みを繰り返して有限時間の極限をとれば、完全正値で跡を保つ写像の極限になる。従って率の正値性が必要なのは、単に占有が正かどうかの検査に留まらない。

:::example[同じ浴が二つの遷移を駆動する場合]
同一周波数の二つの結合で率行列を $\Gamma=\begin{pmatrix}1&r\\r&1\end{pmatrix}$ とする。固有値は $1+r,1-r$ であり、$|r|\le1$ である。$r=1$ なら和の遷移だけが係数 $2$ で残り、差の組の係数は零である。二本を独立に率 $1$ として扱うと、浴が同時に駆動できない組にも減衰を付けてしまう。
:::

### 詳細つり合いから Gibbs 状態を確認する

単一の結合の永年化した生成器で、$[H_S,A_\omega]=-\hbar\omega A_\omega$、$A_{-\omega}=A_\omega^\dagger$ とする。正準状態 $\rho_\beta\propto e^{-\beta H_S}$ では、固有基底での行列要素を比較して $\rho_\beta A_\omega=e^{\beta\hbar\omega}A_\omega\rho_\beta$ である。また $A_\omega^\dagger A_\omega$ と $A_\omega A_\omega^\dagger$ は $H_S$ と可換である。このため下降と上昇の散逸項に $\rho_\beta$ を代入すると、率の関係 $\Gamma(-\omega)=e^{-\beta\hbar\omega}\Gamma(\omega)$ により二つの寄与が相殺する。例えば $A_\omega\rho_\beta A_\omega^\dagger=e^{-\beta\hbar\omega}\rho_\beta A_\omega A_\omega^\dagger$ である。$\omega=0$ の Hermite な飛躍も $\rho_\beta$ と可換であり、散逸項は零である。Lamb シフトも $H_S$ と可換なので、$\rho_\beta$ は定常である。保存するブロックがある場合は他の定常状態があり得るため、これだけで一意性は結論できない。

:::example[縮退した遷移に励起が残る]
三準位で励起 $\ket{e_1},\ket{e_2}$ のエネルギーが等しく、一つの下降演算子 $L=\sqrt\gamma\ket g(\bra{e_1}+\bra{e_2})$ だけがあるとする。差の状態 $\ket d=(\ket{e_1}-\ket{e_2})/\sqrt2$ は $L\ket d=0$ なので、その占有は散逸で変わらない。一方、和の状態には $\sqrt{2\gamma}$ の振幅で下降があるため、占有の減衰率は $2\gamma$ である。零温度の Gibbs 状態は定常でも、結合しない励起が準備できる場合はすべての初期状態がそこへ向かうわけではない。
:::

### 調和系の占有確率を正逆の跳躍で求める

調和振動子の消滅演算子を $a$ とし、飛躍を $\sqrt{\gamma(\bar n+1)}\,a$ と $\sqrt{\gamma\bar n}\,a^\dagger$ とする。占有確率 $p_n=\langle n|\rho|n\rangle$ では $a\ket n=\sqrt n\ket{n-1}$ から、下降率が $\gamma(\bar n+1)n$、上昇率が $\gamma\bar n(n+1)$ である。定常の隣接比は $\bar n/(\bar n+1)=e^{-\beta\hbar\omega_0}$ なので、基底からの占有は幾何分布になる。平均占有の変化は一回の増分 $\pm1$ に率を掛けて足し、$\dv{\langle n\rangle}{t}=-\gamma(\langle n\rangle-\bar n)$ である。系の位置の平均は、随伴の散逸を $a$ に作用させると $-\gamma a/2$ となるので、率 $\gamma/2$ で減る。計算では $[a,a^\dagger]=1$ を使うと、下降の寄与が $-\gamma(\bar n+1)a/2$、上昇の寄与が $+\gamma\bar n a/2$ になる。

:::example[振幅とエネルギーを同じ時間で比べる]
$\gamma=2\,\mathrm{s^{-1}}$、$\bar n=1$ とする。初期の平均占有が $5$ なら $\langle n(t)\rangle=1+4e^{-2t}$ であり、平衡からのずれが半分になるのは $(\log2)/2$ 秒である。平均の変位の包絡は $e^{-t}$ で減るため、同じ時刻では $1/\sqrt2$ が残る。エネルギーの減衰とコヒーレントな振幅の減衰から同じ率を読み取るには、この二倍の違いを考慮する。
:::

### ガウス型の位相緩和は指数と限らない

二準位のエネルギー差に、古典的な平均零の定常 Gaussian 雑音 $\xi(t)$ が線形に入る極限を考える。Hamiltonian を $H(t)=[\hbar\omega_0+\xi(t)]\sigma_z/2$ とする。異なる時刻の Hamiltonian が可換なので、固有基底の非対角成分は各実現で $e^{-i\omega_0t-(i/\hbar)\int_0^t\xi(s)\dd s}$ を掛けた値である。積分した雑音は正規であり、特性関数を平均すると
$$
\frac{\langle\rho_{eg}(t)\rangle}{\rho_{eg}(0)}
=e^{-i\omega_0t}
\exp\left[-\hbar^{-2}\int_0^t(t-s)C_\xi(s)\dd s\right]
$$
となる。これは対角に結合した外部の雑音の模型であり、量子浴の交換子は含めていない。$C_\xi(s)$ が原点で連続なら初期は $t^2$ で減る。$C_\xi$ が可積分なら指数の長い時間の傾きは $-\hbar^{-2}\int_0^\infty C_\xi$ である。従って指数の緩和模型は長い記憶のない範囲で使うとよいが、そのまま短時間の曲率を再現しない。

:::example[同じ位相緩和で初期だけが Gaussian になる]
$C_\xi(s)=C_0e^{-|s|/\tau_b}$ とすれば位相減衰の指数は $-C_0\tau_b[t-\tau_b(1-e^{-t/\tau_b})]/\hbar^2$ である。$t=\tau_b$ のとき、この指数の絶対値は長時間の率を時刻零から使った値の $1/e$ 倍である。占有は全時刻で変わらないので、平均エネルギーだけの測定ではこの緩和を検出できない。
:::

### 高温の局所摩擦を密度行列で表す

調和な系 $H_S=p^2/(2m)+\kappa x^2/2$ で、系の周波数が浴のカットオフより低く、高温の雑音を古典白色とみなせる範囲を考える。初期層を過ぎた Kramers 型の近似は、Wigner 関数 $W(x,p)$ に対し
$$
\partial_t W=-\frac{p}{m}\partial_xW+\kappa x\partial_pW
+\frac{\gamma_0}{m}\partial_p(pW)+\gamma_0k_{\mathrm B}T\,\partial_p^2W
$$
である。調和な自由発展の Wigner 方程式はこの古典の移流と一致する。Wigner 変換で $[x,\rho]$ は $i\hbar\partial_p W$、$\{p,\rho\}$ は $2pW$ に対応する。前者は密度行列の二つの位置の差が Fourier 積分下で $p$ 微分になること、後者は二つの位置の差に関する微分が部分積分で $p$ との積になることから従う。従って対応する高温の Caldeira–Leggett マスター方程式は
$$
\dot\rho=-\frac{i}{\hbar}[H_S,\rho]
-\frac{i\gamma_0}{2m\hbar}[x,\{p,\rho\}]
-\frac{\gamma_0 k_{\mathrm B}T}{\hbar^2}[x,[x,\rho]]
$$
である。

この近似を任意の初期状態へ外挿すると、正値性を破り得る。平均零で交差共分散が零の最小不確定 Gaussian 状態をとり、位置と運動量の分散を $V_{xx},V_{pp}$ とする。そのとき $V_{xx}V_{pp}=\hbar^2/4$ である。上の Kramers 方程式によれば初期の $V_{xx}$ の微分は零、$V_{pp}$ の微分は $-2\gamma_0 V_{pp}/m+2\gamma_0k_{\mathrm B}T$ である。従って $V_{pp}>mk_{\mathrm B}T$ なら直後の共分散の行列式は $\hbar^2/4$ を下回る。物理的な密度演算子ならこの不確定性の下限を満たすため、初期層や量子的な拡散を省く近似の適用範囲がここで現れる。

## 量子揺動散逸定理

### 順序による二つのスペクトル

熱平衡の Hermite 演算子 $A$ に対し $S^>(\omega)=\int e^{i\omega t}\langle A(t)A(0)\rangle\dd t$ とする。固有状態 $n,m$ の寄与は $2\pi p_n|A_{nm}|^2\delta(\omega-(E_m-E_n)/\hbar)$ である。$n,m$ を交換して $p_m/p_n=e^{-\beta\hbar\omega}$ を用いると

$$
S^>(-\omega)=e^{-\beta\hbar\omega} S^>(\omega).
$$

:eqref[eq-neq-kubo] で与えた量子応答の虚部は、交換子の Fourier 変換を $2\hbar$ で割った値である。実際、$f(t)=\langle[A(t),A(0)]\rangle$ は純虚数値の奇関数であり、$2\hbar\Im\chi(\omega)=2i\int_0^\infty f(t)\sin\omega t\,\dd t=S^>(\omega)-S^>(-\omega)$ である。従って対称化された $S_{\rm sym}=[S^>(\omega)+S^>(-\omega)]/2$ に対し

$$
S_{\rm sym}(\omega)=\hbar\coth(\beta\hbar\omega/2)\chi''(\omega)
\label{eq-neq-quantum-fdt}
$$

となる。零周波数で保存成分のデルタ線がある場合は、交換子の関係だけでその線を決められない。正の周波数で $\beta\hbar\omega\ll1$ なら :eqref[eq-neq-classical-fdt] に近づく。

浴を自由 Hamiltonian の Gibbs 状態で用意すると、各振動子で $\langle q_j^2\rangle=\hbar(2n_j+1)/(2m_j\omega_j)$、$n_j=(e^{\beta\hbar\omega_j}-1)^{-1}$ である。従って基準を変位させていない自由力 $\xi$ の対称化スペクトルは $\hbar J(|\omega|)\coth(\beta\hbar|\omega|/2)$ である。零温度でも対称化スペクトルには零点の寄与が残る。正の励起エネルギーに対する上向き遷移に必要な $S^>(-\omega)$ は、同じ極限で零になる。

:::example[温かい浴]
$\beta\hbar\omega=\log3$ なら逆向きのスペクトルは正向きの三分の一である。上昇と下降の率の比がこの値なら定常励起率は $1/4$ である。
:::

:::example[低周波の Ohmic 雑音]
$J(\omega)=\gamma_0\omega$ とすれば $S_{\xi,\rm sym}\simeq2k_{\mathrm B}T\gamma_0$ になる。これは :eqref[eq-neq-brownian-fdt] の古典白色極限である。
:::

### 正負周波数を二準位の検出器で読む

検出器に励起エネルギー $\hbar\Omega$ があり、浴の Hermite 演算子 $B$ と弱く結合するとする。浴を状態 $\ket n$ から $\ket m$ へ移す遷移のエネルギー保存は、検出器が上がるとき $E_m-E_n=-\hbar\Omega$、下がるとき $E_m-E_n=\hbar\Omega$ である。$S^>_{BB}(\omega)=\int e^{i\omega t}\langle B(t)B(0)\rangle\dd t$ の固有状態展開では $E_m-E_n=\hbar\omega$ に重みがある。従って励起率は $S^>(-\Omega)$、緩和率は $S^>(+\Omega)$ に比例する。熱平衡での比は $e^{-\beta\hbar\Omega}$ である。対称化スペクトルだけを同じ率として両方向に入れると、有限温度でも比を $1$ にしてしまう。

:::example[温度を率の比から逆算する]
$\gamma_\uparrow/\gamma_\downarrow=1/4$ なら $\beta\hbar\Omega=\log4$ である。定常検出器は励起確率 $1/5$ をもつ。零温度では励起率が消えるが、浴が検出器からエネルギーを受け取る緩和率は、対応するモードがあれば残る。
:::

### 交換子から得た応答と浴の雑音を照合する

調和な系の復元力を $-\kappa x$ とし、初期層が減衰した定常状態を考える。:eqref[eq-neq-quantum-langevin] を $e^{-i\omega t}$ 規約の周波数に移すと、位置の感受率は
$$
\chi_x(\omega)=\frac1{\kappa-m\omega^2-i\omega\widehat\gamma(\omega)}.
$$
因果核の実部が非負なら $\Im\chi_x(\omega)=\omega|\chi_x|^2\Re\widehat\gamma$ である。浴の定義から正の周波数では $\Re\widehat\gamma=J(\omega)/\omega$ である。これは偶延長した核の Fourier 変換が因果核の実部の二倍になるためである。従って対称な位置スペクトルが線形変換で $|\chi_x|^2$ 倍した浴の力のスペクトルになることと :eqref[eq-neq-quantum-fdt] を使えば、$S^{\rm sym}_{xx}(\omega)=\hbar\coth(\beta\hbar\omega/2)\Im\chi_x(\omega)$ を得る。初期条件の寄与が消えるという仮定は、無減衰の離散モードが残る場合には使えない。

:::example[応答を保ったまま温度だけを変える]
$\widehat\gamma=\gamma_0$ の局所摩擦を近似として使い、正の $\omega=\sqrt{\kappa/m}$ をとると、$\Im\chi_x=1/(\gamma_0\omega)$ である。この周波数の対称スペクトルは高温で $2k_{\mathrm B}T/(\gamma_0\omega^2)$ に近づき、零温度で $\hbar/(\gamma_0\omega)$ に近づく。温度で雑音は変わるが、この調和模型の決定論的な応答は同じである。
:::

### 一つの振動子の順序付き相関を計算する

振動子の座標を $q=q_{\rm zp}(b+b^\dagger)$、$q_{\rm zp}^2=\hbar/(2m_b\omega_b)$ とする。自由発展では $b(t)=b e^{-i\omega_b t}$、$b^\dagger(t)=b^\dagger e^{i\omega_bt}$ である。熱平衡の占有数は $\bar n=(e^{\beta\hbar\omega_b}-1)^{-1}$ であり、$\langle bb^\dagger\rangle=\bar n+1$、$\langle b^\dagger b\rangle=\bar n$ である。二つの消滅または二つの生成の期待値は粒子数基底で非対角なので零になる。従って順序付き相関は
$$
\langle q(t)q(0)\rangle=q_{\rm zp}^2[(\bar n+1)e^{-i\omega_b t}+\bar n e^{i\omega_b t}]
$$
である。実部は対称化した相関であり、虚部に $i$ を掛けた量が交換子の半分に等しい。交換子は $-2i q_{\rm zp}^2\sin\omega_b t$ であり温度によらないが、実部は $2\bar n+1$ に比例する。従って測定の順序で変わる部分と、どちらの順序にも共通の熱の変化を分けられる。

:::example[熱により生じる線と零点の線]
$\beta\hbar\omega_b=\log2$ では $\bar n=1$ であり、正の線と負の線の重みは $2q_{\rm zp}^2$ と $q_{\rm zp}^2$ に比例する。対称化するとどちらも $3q_{\rm zp}^2/2$ に比例し、交換子の正の重みは差の $q_{\rm zp}^2$ に比例する。零温度では負の線が消えるが、正の線とその対称化の半分は残る。
:::

### 同じ温度でも帯域で古典近似が変わる

量子と古典の対称雑音の比は $r(x)=x\coth x$、$x=\beta\hbar|\omega|/2$ である。$x=0$ での正則な極限は $1$ であり、$\sinh x=x+x^3/6+\cdots$、$\cosh x=1+x^2/2+\cdots$ から $r=1+x^2/3-x^4/45+\cdots$ となる。$x\gg1$ では $r=x[1+O(e^{-2x})]$ であり、雑音は古典式に比べて強くなる。観測帯域の上限まで古典近似を用いるなら、その上限で $\beta\hbar|\omega|\ll1$ が必要である。温度の高低を述べるには、同時にどの振動数と比較しているかを指定する。

:::example[二つのモードで近似の精度を比べる]
同じ温度で $\beta\hbar|\omega_1|=0.2$ と $\beta\hbar|\omega_2|=4$ を比較する。前者の比は $0.1\coth0.1\simeq1.0033$、後者は $2\coth2\simeq2.075$ である。低い方のモードで雑音がほぼ古典値でも、高い方のモードの対称雑音は古典の二倍以上になる。実験で対称化したスペクトルを温度へ換算するときは、その周波数域での係数を使う。
:::

### 応答と雑音の周波数積分を等時刻へ戻す

減衰する線形な位置の応答で $S^{\rm sym}_{xx}=\hbar\coth(\beta\hbar\omega/2)\chi_x''$ とする。平均零の等時刻分散は $1/(2\pi)$ 倍した全周波数の積分である。高温近似が積分に寄与する帯域で一様に使え、極限と積分を交換できるなら、
$$
\langle x^2\rangle_{\rm cl}=\frac{2k_{\mathrm B}T}{\pi}
\int_0^\infty\frac{\chi_x''(\omega)}{\omega}\dd\omega
=k_{\mathrm B}T\,\chi_x(0).
$$
最後の等号は瞬時の定数項を持たない位置感受率の静的分散関係である。反項を入れた調和系なら $\chi_x(0)=1/\kappa$ なので古典の等分配に戻る。低温では $\coth$ を $2/(\beta\hbar\omega)$ で置き換えられないため、同じ静的感受率でも全周波数の吸収分布が分散に入る。積分に現れる帯域が切断で抑えられなければ、温度を上げた近似の一様性や微分した観測量の分散をさらに確認する。

:::example[位置と速度で高い周波数の重みが違う]
理想的な Ohmic の局所核で $\kappa,m,\gamma_0>0$ とし、有限温度で高周波を考える。$|\chi_x|^2\sim1/(m^2\omega^4)$、力の対称スペクトルは $\hbar\gamma_0|\omega|$ なので位置のスペクトルは $|\omega|^{-3}$ に比例する。位置分散の高周波積分は収束する。一方、速度には $\omega^2$ が掛かり、同じ近似では $|\omega|^{-1}$ の尾となって速度分散が対数発散する。位置の計算に支障が見えない場合でも、速度の観測にはカットオフの指定が必要になる。
:::

### 熱浴が二つあると一つの温度に直せない

二つの独立な浴で散逸核を $\gamma_1,\gamma_2$、温度を $T_1,T_2$ とする。力の対称スペクトルは独立性で足し算になり、応答の分母には $\widehat\gamma_1+\widehat\gamma_2$ が入る。従って $\omega>0$ における位置スペクトルを $\hbar\Im\chi_x$ で割った値は
$$
\frac{\sum_i \Re\widehat\gamma_i(\omega)\coth(\hbar\omega/(2k_{\mathrm B}T_i))}
{\sum_i\Re\widehat\gamma_i(\omega)}
$$
である。応答の虚部を $\omega|\chi_x|^2\sum_i\Re\widehat\gamma_i$ と書いた。両温度が等しければ平衡の一因子に戻る。古典の周波数域でだけ、重み付き平均温度 $\sum_i T_i\Re\widehat\gamma_i/\sum_i\Re\widehat\gamma_i$ という形にできる。異なる浴の周波数依存が異なれば、その平均自体が周波数で変わる。

:::example[異なる帯域を異なる浴が駆動する]
古典の低周波で $T_1=2T_2$ とする。ある周波数で $\Re\widehat\gamma_1=\Re\widehat\gamma_2$ なら応答と雑音の比から読んだ温度は $3T_2/2$ である。他の周波数で両抵抗の比が $1:3$ なら $5T_2/4$ である。二つの測定はどちらも同じ線形モデルに従うが、一つの平衡温度では同時に説明できない。量子の周波数域では $\coth$ の平均を温度へそのまま移すこともできない。
:::

## 散逸のある量子トンネル

### 虚時間における経路のコスト

[経路積分](/physics/advanced-quantum-mechanics/path-integrals/) を使い、虚時間 $0\le\tau\le\beta\hbar$ の周期経路を $x(\tau)=\sum_n x_n e^{i\nu_n\tau}$、$\nu_n=2\pi n/(\beta\hbar)$ と展開する。調和浴の各 Fourier 成分は $m_j(\nu_n^2+\omega_j^2)|q_{j,n}|^2/2$ と $x_n$ による一次項を含む。平方完成で $q_{j,n}$ を積分すると、行列式は $x$ によらず、反項との差を残す。残りの作用は

$$
S_{\rm diss}[x]=\frac{\beta\hbar}{2}\sum_n |x_n|^2 K_n,\qquad
K_n=\sum_j\frac{c_j^2}{m_j\omega_j^2}\frac{\nu_n^2}{\nu_n^2+\omega_j^2}
=\frac2\pi\int_0^\infty\frac{J(\omega)}{\omega}\frac{\nu_n^2}{\nu_n^2+\omega^2}\dd\omega .
$$

$n=0$ の成分は零であり、これは反項が静的ポテンシャルを保つことと対応する。理想 Ohmic で正の $\nu=|\nu_n|$ を用いると $\int_0^\infty \nu^2/(\nu^2+\omega^2)\dd\omega=\pi\nu/2$ なので $K_n=\gamma_0|\nu_n|$ である。

準安定井戸の最低半古典近似では、井戸に留まる経路との差の作用 $S_{\rm b}$ をもつ虚時間の停留経路で率の指数 $e^{-S_{\rm b}/\hbar}$ を評価する。散逸結合を弱く加えると、非散逸の停留経路 $x_{\rm b}$ での $S_{\rm diss}[x_{\rm b}]$ が作用の一次補正となる。経路の変化による非散逸作用の一次変分は停留性で零だからである。非定数ならこの補正は非負で、先頭の指数を抑える。強結合の率の前因子、二重井戸での多重経路の相互作用や局在転移は、この一経路での摂動から決まらない。一般化摩擦と障壁越えの背景は :cite[zwanzig-noneq] を参照する。

:::example[単一の正弦変位]
$x(\tau)=a\cos\nu_1\tau$ なら $x_{\pm1}=a/2$ で、散逸作用は $\beta\hbar a^2 K_1/4$ である。一定のオフセットを加えても変わらない。
:::

:::example[Drude 切断での高周波]
Drude 型の $J$ を代入して分母を部分分数にすると $K_n=\gamma_0\Omega|\nu_n|/(\Omega+|\nu_n|)$ となる。低周波では Ohmic、高周波では $\gamma_0\Omega$ へ飽和する。
:::

### 三次ポテンシャルのトンネル指数

零温度の非散逸極限で $U(x)=m\omega_0^2x^2(1-x/a)/2$、$a>0$ を考える。$x=0$ が準安定極小で、$x=a$ がエネルギー零の折り返し点である。虚時間の保存式 $m\dot x^2/2-U(x)=0$ を使う。$x_{\rm b}(\tau)=a\,\operatorname{sech}^2[\omega_0(\tau-\tau_0)/2]$ と置くと、$\dot x_{\rm b}^2=\omega_0^2x_{\rm b}^2(1-x_{\rm b}/a)$ であり、端点 $\tau\to\pm\infty$ で井戸底へ戻る。井戸底に留まる経路との差の作用は

$$
S_0=2\int_0^a\sqrt{2mU(x)}\,\dd x
=2m\omega_0 a^2\int_0^1 y\sqrt{1-y}\,\dd y
=\frac{8}{15}m\omega_0a^2 .
$$

往復経路であるため、片道の運動量積分に二が掛かる。散逸を弱く加えたとき、零温度の Ohmic 補正は $\frac{\gamma_0}{2}\int|\omega||\widetilde x_{\rm b}(\omega)|^2\dd\omega/(2\pi)$ である。端点での値が同じなので定数の成分は寄与せず、厳密に正の補正が得られる。$x_{\rm b}(\tau)=a\,g(\omega_0\tau)$ と尺度を変換すると補正は $\gamma_0a^2$ に比例し、元の作用に対する比は $\gamma_0/(m\omega_0)$ に比例すると分かる。実際の有限帯域の浴では $K_n$ の周波数依存を保つ。

### 虚時間核の分母に浴振動数が残る理由

虚時間の周期を $\beta\hbar$ とし、Fourier 係数を $x_n$、Matsubara 周波数を $\nu_n=2\pi n/(\beta\hbar)$ とする。一振動子に由来する二次作用は各 $n$ で
$$
\frac{\beta\hbar}2\left[m_b(\nu_n^2+\omega_b^2)|Q_n|^2
-2c\Re(Q_n^*x_n)+\frac{c^2}{m_b\omega_b^2}|x_n|^2\right]
$$
である。$Q_n$ を平方完成した後、Gaussian 積分は $x_n$ に依存しない因子を残す。$x_n$ の係数は $\frac{c^2}{m_b\omega_b^2}\frac{\nu_n^2}{\nu_n^2+\omega_b^2}$ である。$\nu_n=0$ の静的移動は 反項と相殺し、高周波では振動子が移動に追随できず静的二次項に近づく。多数の振動子なら各係数を足す。

:::example[同じ振幅の二つの虚時間モード]
$m_b=c=\omega_b=1$ とする。$\nu=1$ の係数は $1/2$、$\nu=2$ の係数は $4/5$ である。二次作用に占める比は $8/5$ となる。経路の平均位置だけが違う場合はこの項は変わらず、同じ平均でも時刻方向の形が違えば変わる。
:::

### 虚時間作用を変分して障壁の静的解を調べる

散逸作用を含む $S_E=\int_0^{\beta\hbar}[m\dot x^2/2+U(x)]\dd\tau+S_{\rm diss}[x]$ を考える。周期条件で運動項を部分積分すると、停留条件の Fourier 成分は
$$
(m\nu_n^2+K_n)x_n+[U'(x)]_n=0
$$
である。$U$ が非線形なら最後の項は複数の Fourier 成分を混ぜる。障壁頂上 $b$ に留まる静的な解の近くで $U'(b+\delta x)\simeq-m\omega_b^2\delta x$ とすると、二次の各変分係数は $m(\nu_n^2-\omega_b^2)+K_n$ である。$n=0$ の負の方向は静的な障壁の不安定性である。温度を下げて $|\nu_1|$ が小さくなり、非零の最小周波数で係数が零になると、静的な解に虚時間方向の追加の不安定性が出る。この条件が熱活性化から虚時間で変化する経路への転移の目安になる。ただし非線形の作用の比較により転移が別の温度で起きる場合があるため、この二次条件だけで率の転移点を一般には確定できない。

Ohmic 核 $K_n=\gamma_0|\nu_n|$ では境目の角周波数は $\nu_*=(\sqrt{\gamma_0^2+4m^2\omega_b^2}-\gamma_0)/(2m)$ であり、対応する温度は $T_*=\hbar\nu_*/(2\pi k_{\mathrm B})$ である。$\gamma_0$ を増やすと $\nu_*$ が下がり、局所の不安定性が高い温度では現れにくい。これは実時間での障壁の不安定成長率と同じ二次式だが、こちらでは最低の非零 Matsubara 周波数との比較で温度へ換算する。

:::example[摩擦がない場合とある場合の比較]
$\gamma_0/(m\omega_b)=3/2$ とすれば $\nu_*=\omega_b/2$ である。従って二次安定性が変わる温度は非散逸の $\hbar\omega_b/(2\pi k_{\mathrm B})$ の半分になる。計算は頂上の曲率と低い非零周波数の核を使う。深いトンネルの作用の全体を求めるには、頂上だけでなく井戸底から折り返し点までの経路を決める。
:::

### 三次の障壁で作用を高さに換算する

三次ポテンシャル $U(x)=m\omega_0^2x^2(1-x/a)/2$ の頂上は、$U'=m\omega_0^2 x(1-3x/(2a))=0$ から $b=2a/3$ である。障壁高さは $\Delta U=2m\omega_0^2a^2/27$、曲率は $U''(b)=-m\omega_0^2$ である。従って零温度の非散逸作用は $S_0=36\Delta U/(5\omega_0)$ と書き直せる。熱による指数 $\beta\Delta U$ と $S_0/\hbar$ の比較では、同じ高さでも井戸の曲率が大きければトンネルの無次元作用が小さくなる。ポテンシャル形状を固定しない一般の障壁では、同じ高さと局所の曲率から全範囲の作用を決められるわけではない。

:::example[高さだけを固定して曲率を変える]
この三次の形で $\Delta U/\hbar=10\omega_0$ なら $S_0/\hbar=72$ である。障壁高さを同じにした模型で $\omega_0$ を二倍にすれば作用は $36$ になり、指数には大きな違いが出る。その場合は $a$ が半分に変わるため、障壁の幅も変わっている。率の前因子や有限温度の寄与は、この作用の比だけでは決まらない。
:::

### 零温度の Ohmic 作用を時刻差の積分にする

滑らかで速く減衰する実の経路 $x$ について、Parseval の等式により
$$
\int_\R |x(\tau+s)-x(\tau)|^2\dd\tau
=\int_\R\frac{\dd\omega}{2\pi}\,2(1-\cos\omega s)|\widetilde x(\omega)|^2.
$$
これを正の $s$ で $s^{-2}$ の重みを掛けて積分する。非負なので積分順を交換でき、$\int_0^\infty(1-\cos\omega s)/s^2\dd s=\pi|\omega|/2$ を使う。後者は部分積分して $\omega\int_0^\infty\sin(\omega s)/s\dd s$ とし、$e^{-\epsilon s}$ を掛けた積分が $\arctan(\omega/\epsilon)$ になることから従う。従って
$$
S_{\rm diss}[x]=\frac{\gamma_0}{4\pi}
\int_\R\int_\R\frac{(x(\tau)-x(\tau'))^2}{(\tau-\tau')^2}\,\dd\tau\dd\tau'
$$
である。差が現れるため定数の移動は作用に寄与しない。局所の摩擦仕事のように一時刻の速度の二乗だけで決まらず、異なる二時刻の位置がどれだけずれたかを積分する。微分可能で速く減衰する経路では対角近くの分子が二乗で消え、遠くでも積分できる。段差を無限に急にした経路の作用を求める場合には、高周波のカットオフが必要になることがある。

:::example[時間を伸ばしても Ohmic の補正が変わらない形]
$x_T(\tau)=a\,g(\tau/T)$ とし、$g$ を固定して速く減衰する滑らかな形とする。二重積分で $\tau=T u,\tau'=T v$ と置くと、測度の $T^2$ と分母の $T^2$ が消える。従って $S_{\rm diss}$ は $a^2$ に比例し、$T$ には依存しない。運動項は反対に $a^2/T$、二次のポテンシャル項は $a^2T$ に比例する。この尺度の違いがあるため、散逸の項だけから停留経路の時間幅は決められない。
:::

### 三次の往復経路の散逸係数を求める

井戸底から折り返して戻る非散逸経路 $x_{\rm b}=a\,\operatorname{sech}^2(\omega_0\tau/2)$ について、時刻の中心を零とする。$u=e^{\omega_0\tau}$ と置くと $\operatorname{sech}^2(\omega_0\tau/2)=4u/(1+u)^2$ なので
$$
\widetilde x_{\rm b}(\omega)=\frac{4a}{\omega_0}
\int_0^\infty\frac{u^{i\omega/\omega_0}}{(1+u)^2}\dd u
=\frac{4\pi a\omega}{\omega_0^2\sinh(\pi\omega/\omega_0)}
$$
となる。積分は Euler の Beta 積分により $\Gamma(1+iy)\Gamma(1-iy)$、$y=\omega/\omega_0$ であり、:ref[thm-cpx-gamma-reflection] の相補公式と $\Gamma(1+iy)=iy\Gamma(iy)$ から $\pi y/\sinh(\pi y)$ になる。$\omega\to0$ で $4a/\omega_0$ となり、元の経路の面積と一致する。この変換を Ohmic の作用に代入し、$y=\pi|\omega|/\omega_0$ とおくと
$$
S_{\rm diss}[x_{\rm b}]=\frac{8\gamma_0 a^2}{\pi^3}\int_0^\infty\frac{y^3}{\sinh^2 y}\dd y.
$$
$1/\sinh^2 y=4e^{-2y}/(1-e^{-2y})^2=4\sum_{n\ge1}n e^{-2ny}$ の各項は非負であるため積分順を交換できる。$\int_0^\infty y^3e^{-2ny}\dd y=6/(2n)^4$ なので、積分は $(3/2)\sum_{n\ge1}n^{-3}$ となる。従って $S_{\rm diss}=12\zeta(3)\gamma_0a^2/\pi^3$ である。ここで $\zeta(3)=\sum_{n\ge1}n^{-3}$ である。非散逸の作用で割ると $[45\zeta(3)/(2\pi^3)]\gamma_0/(m\omega_0)$ であり、角括弧の係数は約 $0.872$ である。この値は非散逸の経路を固定して評価した結合の一次補正である。

:::example[小さい相対作用でも率に差が出る]
$S_0/\hbar=30$、$\gamma_0/(m\omega_0)=0.02$ とする。作用の相対増加は約 $0.0174$ だが、無次元の作用の増加は約 $0.523$ である。前因子の変化を無視した指数の比は約 $e^{-0.523}=0.59$ になる。元の作用に比べて小さい補正でも、トンネル率そのものの相対変化は小さいとは限らない。
:::

### 有限温度でエネルギーの異なる経路を比べる

非散逸の一次元の半古典近似で、井戸底をエネルギー零とし、障壁を越える途中の各エネルギー $E$ について転回点 $x_-(E),x_+(E)$ を取る。障壁に入る前の状態を熱平衡で選ぶと確率の指数に $\beta E$ が入る。一方、禁止領域の片道の波動関数の減衰は近似的に $\exp[-\int_{x_-}^{x_+}\sqrt{2m(U-E)}\,\dd x/\hbar]$ であり、確率はその二乗である。これを得るには Schrödinger の式に $\psi=e^{-S/\hbar}$ を入れ、$\hbar S''$ を $S'^2$ に比べて小さいとして $S'=\sqrt{2m(U-E)}$ を選ぶ。転回点の接続は前因子に効くが、ここでは先頭の指数だけを求める。

従って各エネルギーの支配指数は
$$
\mathcal B(E)=\beta E+S(E)/\hbar,\qquad
S(E)=2\int_{x_-}^{x_+}\sqrt{2m(U(x)-E)}\,\dd x .
$$
$S'(E)=-2\int_{x_-}^{x_+}m/\sqrt{2m(U-E)}\,\dd x$ であり、端点で被積分関数が零なので端の寄与はない。右辺の負号を取ったものは虚時間で往復する周期である。従って内部の停留エネルギーは虚時間周期 $-S'(E)=\beta\hbar$ を満たす。端点に最小があれば、底のトンネルと頂上からの熱活性化のどちらかが支配する。複数の停留点があればそれぞれの指数を比較する。実際の散逸浴がある場合、全てを一変数 $E$ の WKB 積分だけで表すことはできず、虚時間の非局所作用を使う。

:::example[放物線の頂上で熱とトンネルを比較する]
頂上近傍を $U(x)=\Delta U-m\omega_b^2 x^2/2$ とする。両端の振幅は $\sqrt{2(\Delta U-E)/(m\omega_b^2)}$ であり、半円の面積として平方根を積分すると $S(E)=2\pi(\Delta U-E)/\omega_b$ となる。従って $\mathcal B(E)$ の傾きは $\beta-2\pi/(\hbar\omega_b)$ である。高温側でこの傾きが負なら、頂上付近では $E$ を上げるほど指数が小さくなり、熱活性化が有利である。低温側では傾きが正になり、頂上より低いエネルギーでの経路を考える必要がある。これは局所の放物線の計算であり、井戸の底まで同じ曲率を延長する近似にはなっていない。
:::

### 二重井戸でトンネルの対が相互作用する

二つの安定な井戸の間隔を $d$ とし、障壁が高いため虚時間の経路がほとんどの時刻でどちらかの井戸に留まる場合を考える。幅 $\tau_c$ の短い遷移で井戸を移るとし、往復の時刻差を $T_{\rm b}\gg\tau_c$ とする。井戸の一方を基準にして、往復を高さ $d$、長さ $T_{\rm b}$ のパルスで近似すると、低周波での Fourier 変換の絶対値の二乗は $2d^2(1-\cos\omega T_{\rm b})/\omega^2$ である。有限の遷移幅を $|\omega|\lesssim\tau_c^{-1}$ の切断で表せば、Ohmic の作用は
$$
\frac{S_{\rm diss}}{\hbar}
\simeq \frac{\gamma_0d^2}{\pi\hbar}
\int_0^{1/\tau_c}\frac{1-\cos(\omega T_{\rm b})}{\omega}\dd\omega
=2\alpha_{\rm diss}\log(T_{\rm b}/\tau_c)+O(1),\qquad
\alpha_{\rm diss}=\frac{\gamma_0d^2}{2\pi\hbar}.
$$
$\omega T_{\rm b}<1$ の積分は二次の余りから有界である。それ以上の周波数では $1/\omega$ の積分が対数を作り、余弦の積分は部分積分により有界である。従って遠く離れた往復の重みは $(T_{\rm b}/\tau_c)^{-2\alpha_{\rm diss}}$ となる。

遷移が $N$ 回あり、時刻 $t_i$ の跳びを $d s_i$、$s_i=\pm1$ とする。周期経路では $\sum_i s_i=0$ である。位置の Fourier 変換に $d\sum_i s_i e^{i\omega t_i}/\omega$ が入るため、絶対値の二乗を展開し、$\sum_i s_i=0$ から定数を引くと、遠く離れた対の作用は $-2\alpha_{\rm diss}\sum_{i<j}s_i s_j\log(|t_i-t_j|/\tau_c)$ となる。各遷移に無次元の小さい重み $y$ を割り当てる。時間の尺度を $b>1$ 倍すると時刻積分は $b^N$ 倍となる。一方、$\sum_{i<j}s_i s_j=-N/2$ なので、重みの作用部分は $b^{-N\alpha_{\rm diss}}$ 倍である。従って希薄な遷移の先頭の尺度則は $y\mapsto b^{1-\alpha_{\rm diss}}y$ となる。強い散逸で低い周波数のトンネルが抑えられる理由が分かる。ここで使ったのは遷移の幅より長い時間の対数相互作用と小さい $y$ の展開であり、接近した対による係数の繰り込みまでは含めていない。

:::example[往復の間隔と重み]
$\alpha_{\rm diss}=1/2$ なら、長い往復の間隔を十倍にしたときに散逸による重みは先頭の次数で十分の一になる。一方、遷移一回の重みの尺度則は $b^{1/2}$ であり、時刻積分も含めた繰り込みでの増減と、固定した二時刻の経路の抑制は異なる。
:::

## 演習問題

::::exercise{level=basic}
浴モードの静的な変位を求めよ。
:::solution
$q_j$ による平方を最小にする値は $c_jx/(m_j\omega_j^2)$ である。これを代入すると浴のポテンシャルは零になり、系には $U(x)$ が残る。
:::
::::
::::exercise{level=basic}
$\gamma(t)=8e^{-4t}$、$t>0$ の局所摩擦係数を求めよ。
:::solution
$\int_0^\infty \gamma(t)\dd t=8/4=2$ である。この近似が使えるのは、速度の変化が $1/4$ の時間より遅い場合である。
:::
::::
::::exercise{level=basic}
上昇率と下降率の比が $1/7$ の二準位系の定常励起確率を求めよ。
:::solution
励起の流入 $u(1-p_e)$ と流出 $dp_e$ を等置し、$p_e=u/(u+d)=1/8$ を得る。
:::
::::
::::exercise{level=basic}
位相緩和の率係数 $r=3$ のとき非対角成分の初期値に対する比を時刻 $t$ の関数として求めよ。
:::solution
$A_0\rho A_0$ によって非対角成分の符号が変わり、元を引くと $-2\rho_{12}$ となる。従って比は $e^{-6t}$ である。
:::
::::
::::exercise{level=advanced}
一つの調和浴の力の交換子を求めよ。
:::solution
$q^{\rm h}=q\cos\omega t+p\sin\omega t/(m_b\omega)$ を使う。交差項から $i\hbar c^2[\cos\omega t\sin\omega s-\sin\omega t\cos\omega s]/(m_b\omega)$ を得る。従って $-i\hbar c^2\sin[\omega(t-s)]/(m_b\omega)$ となる。
:::
::::
::::exercise{level=advanced}
$\beta\hbar\omega\ll1$ で量子 FDT と古典 FDT の比の最初の補正を求めよ。
:::solution
$z=\beta\hbar\omega/2$ とすれば、係数の比は $z\coth z=1+z^2/3+O(z^4)$ である。従って補正は $(\beta\hbar\omega)^2/12$ である。
:::
::::
::::exercise{level=advanced}
一つの浴モードの $K_n$ の周波数極限を比較せよ。
:::solution
$K_n=[c^2/(m_b\omega_b^2)]\,\nu_n^2/(\nu_n^2+\omega_b^2)$ である。低周波の先頭は $c^2\nu_n^2/(m_b\omega_b^4)$、高周波の極限は $c^2/(m_b\omega_b^2)$ である。零モードが消える条件も満たす。
:::
::::
::::exercise{level=advanced}
二準位の遷移と位相緩和を同時に入れる。非対角成分の減衰率を求めよ。
:::solution
上下の飛び移りの項 $A_\omega\rho A_\omega^\dagger$ は対角成分だけに寄与する。二つの反交換子は非対角成分にそれぞれ $-u/2$、$-d/2$ を掛ける。さらに位相緩和が $-2r$ なので全減衰率は $(u+d)/2+2r$ である。
:::
::::

## 確認テスト

::::quiz
反項を含む調和浴の模型で、浴を静的に平衡位置へ動かした場合の系のポテンシャルを選べ。

- [x] $U(x)$
- [ ] $U(x)+x^2\sum_j c_j^2/(2m_j\omega_j^2)$
- [ ] $U(x)-x^2\sum_j c_j^2/(2m_j\omega_j^2)$

:::explanation
平衡位置では $q_j=c_jx/(m_j\omega_j^2)$ であり、Hamiltonian に含まれる浴の平方が零になる。
:::
::::
::::quiz{answer="2" tol="1e-6"}
$u=1,d=3,r=0$ の二準位系で、非対角成分の減衰率を答えよ。

:::explanation
上昇と下降の反交換子の二項を足すと係数は $(u+d)/2=2$ である。
:::
::::
::::quiz{answer="0.5" tol="1e-6"}
$\beta\hbar\omega=\log2$ のとき $S^>(-\omega)/S^>(\omega)$ を答えよ。

:::explanation
固有状態の Boltzmann 重みの比が $e^{-\beta\hbar\omega}=1/2$ を与える。
:::
::::
::::quiz
熱浴の零温度極限について正しいものを選べ。

- [x] 正周波数の対称化雑音には零点の寄与が残り得る。
- [ ] 上昇率と下降率は必ず一致する。
- [ ] 浴の力はすべての時刻で可換になる。

:::explanation
対称化相関は二つの演算子順序を平均するため、正周波数の寄与が残る。一方、正の励起エネルギーに対応する負周波数は熱重みで抑えられる。
:::
::::

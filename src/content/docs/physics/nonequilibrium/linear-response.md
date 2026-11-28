---
title: 線形応答理論
description: 外力に対する因果的な応答を平衡相関で表し、感受率、揺動散逸定理と輸送係数を導く。
notation:
  - symbol: 'H_0-f(t)B'
    meaning: 外場 f によって摂動された Hamiltonian
  - symbol: 'R_{AB}(t)'
    meaning: B に共役な場に対する A の応答核
  - symbol: '\chi_{AB}(\omega)'
    meaning: 複素感受率
  - symbol: 'S_{AA}(\omega)'
    meaning: 中心化した観測量の両側スペクトル
---

## この章の目標

- 因果的な応答を畳み込みとして書く。
- 古典と量子の久保公式を導く。
- 感受率の虚部を揺らぎおよび散逸と結び付ける。
- 輸送係数の相関積分と時間反転の対称性を計算する。

## 力学的摂動と応答関数

### 摂動と因果性

平衡にある系に $H(t)=H_0-f(t)B$ という摂動を加える。$fB$ はエネルギーの次元をもち、$A$ 自身には $f$ への明示的な依存がないと仮定する。初期時刻より前に $f=0$ とし、逆温度 $\beta=(k_{\mathrm B}T)^{-1}$ での $H_0$ の正準平衡から始める。一次の変化 $\delta\langle A(t)\rangle$ は $f$ の線形汎関数である。入力を短いパルスに分けると、平衡ではパルスからの時差だけが応答を決め、因果性により未来の場は寄与しない。従って

$$
\delta\langle A(t)\rangle=\int_{-\infty}^t R_{AB}(t-s)f(s)\,\dd s .
\label{eq-neq-response}
$$

:term[線形応答]{en="linear response" yomi="せんけいおうとう"} とはこの一次の項である。半直線 Fourier 変換を $\chi(\omega)=\int_0^\infty R(t)e^{i\omega t}\dd t$ と定めると、$f(t)\propto e^{-i\omega t}$ に対し $\delta A(t)=\chi(\omega)f(t)$ である。瞬時応答がある場合も因果的な分布として扱え、後の分散関係ではその定数項を別に引く。

:::example[段差応答]
$R(t)=a e^{-\lambda t}$、$t>0$ に対し、時刻零からの一定場は $a f(1-e^{-\lambda t})/\lambda$ を与える。段差応答の時間微分を $f$ で割れば力積応答になる。
:::

:::example[移動度と感受率]
:eqref[eq-neq-brownian-velocity] による速度応答は $e^{-t/\tau}/m$ で、積分は $1/\gamma$ である。変位を出力に選べばこの核を時間積分するため、自由粒子の静的な変位感受率は有限にならない。
:::

### 観測量そのものが外場に依存するとき

Hamiltonian を $H(f)=H_0-fB+f^2C/2+O(f^3)$ とし、外場に共役な観測量を $A(f)=-\partial_fH(f)$ とする。外場による期待値の一次変化には、状態が変化して $B$ の期待値を変える項と、観測量が $-fC$ だけ変わる項がある。同じ平衡状態で $C$ を平均すれば後者は $-\langle C\rangle_0 f(t)$ である。従って周波数表示では $\chi_{BB}(\omega)-\langle C\rangle_0$ となる。この定数項を含む感受率は高周波で零に近づくとは限らない。分散関係を適用するときは、定数項を物質の時間遅れのある応答から引く。

:::example[平衡で消える静的応答]
古典調和系 $H(f)=p^2/(2m)+\kappa(x-f)^2/2$ では $B=\kappa x$、$C=\kappa$ である。静的に $x$ の平均は $f$ なので、状態の変化は $\kappa f$ を与える。観測する共役力 $\kappa(x-f)$ の平均は、観測量に含まれる $-\kappa f$ と相殺して零になる。瞬間的に中心をずらした直後は状態がまだ移動しておらず、共役力は $-\kappa f$ だけ変わる。
:::

### 線形性を時刻の記録から検証する

時間一様な基準状態で、摂動を時刻 $s$ に入れたときの応答を $R(t-s)$ とする。入力の大きさを $\varepsilon$ 倍した期待値を $\langle A\rangle_\varepsilon$ と書き、$\delta A(t)=\left.\partial_\varepsilon\langle A(t)\rangle_\varepsilon\right|_0$ を一次応答とする。この微分は入力の振幅を変える操作であり、時刻についての微分とは違う。有限時間で滑らかに依存する系では、入力を短い矩形に分け、応答の和を取ると最後に畳み込みになる。各矩形の幅を $h$ とし、その面積を $f(s)h$ とすれば、和 $\sum_s f(s)R(t-s)h$ は過去の時刻だけを含む。核が可積分なら、有限サポートの連続な入力に対してこの和は積分へ収束する。応答に瞬時項があるときは、その部分を積分核の連続な部分とは別に $c f(t)$ として加える。

段差入力 $f(t)=f_0\mathbf1_{\{t\ge0\}}$ に対して、瞬時項を除いた応答は $f_0\int_0^tR(s)\dd s$ である。測定の時間微分から $R$ を復元できるが、長い時間の誤差や瞬時誤差も微分で強調される。幅 $w$ の矩形入力なら、無限に鋭い力積を作らなくても $f_0\int_{\max(0,t-w)}^t R(s)\dd s$ を検査できる。有限幅が応答より長ければ面積だけで決まらないため、振幅と幅を独立に指定する。

:::example[同じ面積の矩形でもピークが違う]
$R(t)=e^{-at}$（$t>0$）とし、力積を $I$ とする。幅 $w$ の矩形入力を高さ $I/w$ で加えると、停止時点の応答は $I(1-e^{-aw})/(aw)$ である。$aw=\log2$ のときは $I/(2\log2)\simeq0.721I$ となる。幅を零へ近づけた極限の $I$ より小さいのは、最初に加えた部分が停止時までに減衰しているためである。停止直後からは同じ係数に $e^{-a(t-w)}$ を掛ける。
:::

### 複数の場を同じ周波数で加える

Hamiltonian に $-\sum_j f_j B_j$ を加え、各力に共役な観測量 $B_i$ を測る。複素振幅を列 $F$ とすれば $B$ の一次応答は $\Re(\chi F e^{-i\omega t})$ である。入力仕事率の周期平均は
$$
-\sum_i\overline{\dot f_i\,\delta B_i}
=\frac{\omega}{2}\Im(F^\dagger\chi F)
=\frac\omega2 F^\dagger\frac{\chi-\chi^\dagger}{2i} F
$$
である。第二式は $F^\dagger\chi^\dagger F=(F^\dagger\chi F)^*$ を使った。すべての位相と振幅の組で仕事率が非負であるための条件は、この Hermite 行列が半正定値であることである。単一成分の虚部が正でも、二つの場が同時に変化する場合の交差項まで正とは限らない。

:::example[場の相対位相が吸収を変える]
$\omega=1$ で吸収行列を $\begin{pmatrix}2&i\\-i&2\end{pmatrix}$ とする。固有値は $1,3$ で両方とも正である。$F=(1,i)^\mathsf T$ では行列を掛けると $F$ 自身になるため、仕事率は $1$ である。$F=(1,-i)^\mathsf T$ では $3F$ となり、仕事率は $3$ である。各場の絶対値は同じであるが、二つの場の位相関係で吸収が変わる。
:::

### 有限のパルスがする仕事を周波数で判定する

外場 $f(t)$ を実で滑らか、時間にコンパクトな台を持つ関数とし、摂動を $-f(t)B$ とする。仕事の二次項は $-\int \dot f(t)\delta B(t)\dd t$ である。場が始点と終点で零なので、平衡時の一定の $\langle B\rangle_0$ による一次の仕事は消えている。$\widehat f(\omega)=\int e^{i\omega t}f(t)\dd t$ とすると $\widehat{\dot f}=-i\omega\widehat f$ であり、応答の変換は $\chi(\omega)\widehat f(\omega)$ である。Parseval の等式から仕事は
$$
W^{(2)}=\frac1{2\pi}\int_\R \omega\,\chi''(\omega)|\widehat f(\omega)|^2\dd\omega .
$$
実の因果応答なら $\chi''$ は奇関数であり、$\omega\chi''(\omega)$ は偶関数である。従って正周波数で $\chi''\ge0$ なら、どの有限パルスでも二次の正味仕事は非負である。複数の場では $\chi''$ を吸収の Hermite 行列で置き換え、パルスの周波数ごとの複素振幅を二次形式に入れる。保存的な瞬時応答が実定数なら、上の積分に寄与しない。時間表示でも $-\int \dot f\,\chi_\infty f=-(\chi_\infty/2)[f^2]_{-\infty}^\infty=0$ である。

:::example[指数緩和でパルス幅を変える]
$R(t)=Ae^{-at}$（$t>0$、$A,a>0$）とし、$f_0$ を時刻 $0$ から $w$ までだけ加える。矩形入力の跳躍は滑らかな切替の極限で読む。立上がり前の $\delta B$ は零、立下がり直前の応答は $f_0 A(1-e^{-aw})/a$ なので、正味の二次仕事は $f_0^2 A(1-e^{-aw})/a$ である。同じ強さで $w$ を長くしてもこの値は $f_0^2A/a$ に飽和する。一方、面積 $I=f_0w$ を固定したまま幅を大きくすると $f_0$ が小さくなり、仕事は $I^2 A/(aw^2)$ へ減る。面積固定と振幅固定では同じパルスの比較にならない。
:::

### 複数時刻の観測から応答を同定する

離散の刻み $h$ で入力 $f_n$ と出力 $y_n$ を記録し、過去 $m+1$ 点で閉じる近似
$$
y_n=\sum_{j=0}^m r_j f_{n-j}
$$
を考える。係数 $r_j$ には刻みを吸収した。無雑音なら、入力の過去の列を並べた行列 $F_{nj}=f_{n-j}$ を使い、線形方程式 $Fr=y$ を解けばよい。解を一意にするための必要十分条件は $F$ の列が線形独立なことである。実際、二つの解の差は $F$ の核に入り、核のベクトルを任意の解に加えても出力は変わらない。例えば一つの周波数だけなら、すべての列は $\sin$ と $\cos$ の二列が張る空間に入るため、三つ以上の係数は一意に決まらない。周波数応答が一点で決まることと、時間核が決まることは区別される。

:::example[周期入力だけが残す自由度]
$m=2$、入力を $f_n=(-1)^n$ とすると $y_n=(-1)^n(r_0-r_1+r_2)$ である。$r=(1,0,0)$ と $r=(0,0,1)$ はこの入力に対して同じ応答を返す。一方、時刻零だけに単位入力を与えれば、三時刻の応答は $(1,0,0)$ と $(0,0,1)$ であり、遅れの差をそのまま観測できる。実験でインパルスに近い入力を作れない場合は、異なる周波数や切替時刻を組み合わせて列を独立にする。
:::

## 久保公式

### 相互作用表示による量子応答

時間発展の相互作用表示で $B_I(s)=e^{iH_0s/\hbar}B e^{-iH_0s/\hbar}$ とする。状態の von Neumann 方程式は $i\hbar\dot\rho_I=[-f(s)B_I(s),\rho_I]$ である。右辺で $\rho_I$ を :ref[prop-stm-quantum-gibbs-minimum] の平衡密度行列 $\rho_0=e^{-\beta H_0}/Z$ に置き換えると、
$\delta\rho_I(t)=(i/\hbar)\int^t f(s)[B_I(s),\rho_0]\dd s$ である。$A_I(t)$ と跡をとり、巡回性と平衡の定常性を使うと :term[久保公式]{en="Kubo formula" yomi="くぼこうしき"} を得る。

$$
R_{AB}(t)=\frac{i}{\hbar}\Theta(t)\langle[A(t),B(0)]\rangle_0. \label{eq-neq-kubo}
$$

ここで $\Theta(t)$ は $t>0$ に対するステップ関数である。摂動 Hamiltonian の符号を反転すれば応答も反転する。

古典では :eqref[eq-neq-liouville] を一次まで積分する。摂動が $-fB$ なので密度への源は $f\{\rho_0,B\}=f\beta\rho_0\{B,H_0\}$ となる。源を無摂動で伝播し $A$ と積分すると

$$
R_{AB}(t)=\beta\Theta(t)\langle A(t)\dot B(0)\rangle_0
=-\beta\Theta(t)\partial_t\langle\delta A(t)\delta B(0)\rangle_0.
\label{eq-neq-classical-kubo}
$$

最後の等号は時刻の平行移動に対する不変性から得る。相関が長時間で零に向かうという条件の下で $\chi_{AB}(0)=\beta\langle\delta A\delta B\rangle_0$ になる。保存量の成分が残るときには長時間極限を引く必要がある。

:::example[調和振動子への力積]
$A=B=x$ とし、無摂動の $x(t)=x(0)\cos\omega_0t+p(0)\sin\omega_0t/(m\omega_0)$ を入れる。交換子が $-i\hbar\sin\omega_0t/(m\omega_0)$ なので $R(t)=\Theta(t)\sin\omega_0t/(m\omega_0)$ である。運動方程式の力積による初速度の増加と一致する。
:::

:::example[調和井戸での静的変位]
過減衰で $C_{xx}(t)=(k_{\mathrm B}T/\kappa)e^{-\kappa |t|/\gamma}$ とすれば、$R_{xx}(t)=e^{-\kappa t/\gamma}/\gamma$ である。積分が $1/\kappa$ になり、力とばねの釣り合いから求めた変位と合う。
:::

### 固有状態での応答の符号

$H_0\ket n=E_n\ket n$、$p_n=e^{-\beta E_n}/Z$ とする。行列要素を $A_{nm}=\langle n|A|m\rangle$、遷移周波数を $\omega_{mn}=(E_m-E_n)/\hbar$ と書く。:eqref[eq-neq-kubo] の二項を固有状態で展開し、第二項の添字を交換すると
$$
R_{AB}(t)=\frac{i}{\hbar}\Theta(t)\sum_{n,m}(p_n-p_m)A_{nm}B_{mn}e^{-i(E_m-E_n)t/\hbar}
$$
である。収束因子 $e^{-\epsilon t}$ を加えて積分すれば、各項は $(p_n-p_m)A_{nm}B_{mn}/[\hbar(\omega_{mn}-\omega-i\epsilon)]$ になる。$A=B$ が Hermite なら分子の行列要素は絶対値の二乗である。$\epsilon\downarrow0$ で虚部は共鳴周波数のデルタ測度となる。正の周波数では $p_n\ge p_m$ なので、平衡系の吸収は非負である。

:::example[二準位系の静的応答]
エネルギー $0,\hbar\omega_0$ の二準位を、行列要素 $b$ の非対角観測量で結ぶ。正負の周波数の項を足すと
$$
\chi(\omega)=\frac{(p_0-p_1)|b|^2}{\hbar}
\left(\frac1{\omega_0-\omega-i0}+\frac1{\omega_0+\omega+i0}\right).
$$
$p_0-p_1=1/2$、$b=1$ なら $\chi(0)=1/(\hbar\omega_0)$ である。負周波数の項を捨てると静的感受率を半分にしてしまう。正周波数の吸収を近似する場合と、全周波数を使う場合で必要な項が異なる。
:::

### 定常 Markov 過程から応答を求める

有限状態の連鎖の前進生成行列を列ベクトルに作用させ、外場 $f(t)$ のもとで $W(t)=W_0+f(t)W_1+O(f^2)$ とする。$W_0\pi=0$、$\pi_i>0$ とし、遠い過去で $\pi$ にあるように準備する。一次変分 $\delta p$ は $\partial_t\delta p=W_0\delta p+f(t)W_1\pi$ であり、積分因子から
$$
\delta\langle A(t)\rangle=\int_{-\infty}^t f(s) A^\mathsf T e^{(t-s)W_0}W_1\pi\,\dd s .
$$
従って $B_i=(W_1\pi)_i/\pi_i$ とおけば応答核は定常相関 $\langle A(t)B(0)\rangle_\pi$ である。列和が零だから $\langle B\rangle_\pi=0$ である。$B$ は外場が平衡エネルギーを変えるかどうかに関係なく定義できる。

各固定 $f$ で $\pi(f)\propto e^{-\beta(E_i-fC_i)}$ が定常であれば $W_1\pi+W_0\pi'=0$ である。$\pi'=\beta(C-\langle C\rangle_\pi)\pi$ を代入すると核は $-\beta$ 倍した $A$ と $C$ の時間相関の微分になる。非平衡定常状態ではこの $\pi'$ の形が一般には使えず、同じ $W_0$ でも $W_1$ の選び方で応答が変わる。

:::example[二種類の場で静的応答を区別する]
二状態で元の $0\to1,1\to0$ の率をともに $1$ とする。両率を同じ $1+f$ 倍にすれば $W_1\pi=0$ であり、定常占有からの一次応答は零である。一方 $0\to1$ の率だけを $1+f$ にすれば、$p_1$ の一次式は $\partial_t\delta p_1=-2\delta p_1+f/2$ となり、一定場に対して $\delta p_1=f/4$ である。元の二状態の相関が同じでも、制御がどの率を変えたかによって感受率は変わる。
:::

### 平衡感受率を虚時間の積分で求める

有限次元で $H(f)=H_0-fB$ とする。非可換な演算子の指数を微分するため、$X,Y$ を有限行列として $F(s)=e^{(1-s)X}e^{s(X+\epsilon Y)}$ と置く。$s$ で微分すると $\epsilon e^{(1-s)X}Y e^{s(X+\epsilon Y)}$ であり、$0$ から $1$ の差を $\epsilon$ で割って極限を取れば
$$
\left.\partial_\epsilon e^{X+\epsilon Y}\right|_0=\int_0^1 e^{(1-s)X}Ye^{sX}\dd s .
$$
従って $X=-\beta H_0$、$Y=\beta B$ として、$A$ が場に依存しなければ正準平均の微分は
$$
\chi_{AB}^{\rm iso}=\int_0^\beta
\langle \delta B(-i\hbar u)\delta A(0)\rangle_0\,\dd u,
\qquad B(-i\hbar u)=e^{uH_0}Be^{-uH_0}
$$
となる。規格化 $Z(f)$ を微分して引く項が $\delta A,\delta B$ を作る。$B$ が $H_0$ と可換なら虚時間に依存せず、$\beta$ 倍の等時刻共分散に戻る。可換でない場合には積の順序を保った虚時間積分が必要である。

これを固有状態で展開すると、同じエネルギーの組の積分は $\beta$、エネルギー差がある組では $(e^{u(E_n-E_m)})$ の積分になる。Boltzmann 重みを掛けると後者の係数は $(p_n-p_m)/(E_m-E_n)$ である。これは力学的な久保公式の零周波数に現れる項と同じである。一方、エネルギーの等しい組は交換子の平衡重みの差が零で、力学的な応答には入らない。従って孤立系で無限にゆっくり場を変える応答と、各場で熱浴により再平衡化した応答を比較すると、準位占有が変わるかどうかが差を生みうる。温度固定の導関数を実験の静的応答と同定するには、浴による再分配の時間を確認する。

:::example[可換な二準位で応答の種類が分かれる]
$B=\diag(-1,1)$、$H_0=0$ とする。$H(f)=-fB$ での平衡平均は $\tanh(\beta f)$ であり、温度固定の係数は $\beta$ である。一方、孤立した時間発展では $B$ が $H(f(t))$ とすべての時刻で可換なので、$\langle B\rangle$ は初期値のままである。力学的応答の交換子も零になる。外場がエネルギー差を作っても、占有を移す機構を入れなければ温度固定の値へは進まない。
:::

### 古典粒子で久保の符号を計算する

Hamiltonian への摂動を $-f(t)B$ とする。密度の Liouville 式の一次の源は $f(t)\{\rho_0,B\}$ である。$\rho_0\propto e^{-\beta H_0}$ に対して $\{\rho_0,B\}=-\beta\rho_0\{H_0,B\}=\beta\rho_0\dot B$ である。無摂動の流れでこの源を時刻 $s$ から $t$ に伝播させ、積分で観測量を初期点へ移せば $R_{AB}(t-s)=\beta\langle A(t-s)\dot B(0)\rangle_0$ となる。定常性を用いて $\langle A(t)\dot B(0)\rangle=-\partial_t\langle A(t)B(0)\rangle$ である。このため $\langle A(t)B(0)\rangle$ が時間につれて減る場合、場と同じ符号の応答が出る。

:::example[調和振動子で変位の応答を再現する]
$B=A=x$、$H_0=p^2/(2m)+m\Omega^2x^2/2$ とすると、相関は $k_{\mathrm B}T\cos\Omega t/(m\Omega^2)$ である。従って $R_{xx}(t)=\sin\Omega t/(m\Omega)$（$t>0$）である。時刻零に力積 $I$ を入れた場合、正準方程式から $\delta p=I$、その後の変位は $I\sin\Omega t/(m\Omega)$ となり、応答核と一致する。無減衰模型なので、周波数変換には収束因子が必要である。
:::

### 二次までの応答で時刻の順序を確かめる

相互作用表示で $\dot\rho_I=(i/\hbar) f(t)[B_I(t),\rho_I]$ とする。有限次元で場と演算子が有界な場合、積分方程式を二回反復すると
$$
\rho_I(t)=\rho_0+\frac{i}{\hbar}\int_0^t f(s)[B_I(s),\rho_0]\dd s
+\left(\frac{i}{\hbar}\right)^2
\int_{0<u<s<t}f(s)f(u)[B_I(s),[B_I(u),\rho_0]]\dd u\dd s+\cdots .
$$
三次以降の項は、交換子のノルムを $2\|B\|$ 倍で抑え、時間領域の体積が $t^n/n!$ であることを使えば、固定時刻で場の振幅の三次以上となる。観測量を掛けて跡を取ると、$\Tr(A[B,\rho])=\Tr([A,B]\rho)$ により二次の核は $\langle [[A_I(t),B_I(s)],B_I(u)]\rangle_0$ である。入れ子の順序は新しい摂動 $s$ が先、古い摂動 $u$ が後になる。一般に演算子を交換できないため、確率のない単純な平方で二次の応答を表せない。

:::example[線形では消える変化を二次で数える]
二準位で初期の密度行列を $\rho_0=(I+p\sigma_z)/2$ とし、$B=\sigma_x$ の短いパルスを面積 $I_f$ だけ与える。無摂動の発展を無視できる幅の極限では、$U=e^{iI_f\sigma_x/\hbar}$ である。$\sigma_x^2=I$ なので $U=\cos u\,I+i\sin u\,\sigma_x$、$u=I_f/\hbar$ である。その後の $\langle\sigma_z\rangle=p\cos2u$ だから変化は $-2pu^2+O(u^4)$ となる。交換子の期待値は初期の対角状態で零であり、二重交換子 $[[\sigma_z,\sigma_x],\sigma_x]=4\sigma_z$ が二次の係数を与える。一つの短いパルスでは時間順序の領域は正方形の半分なので $I_f^2/2$ が掛かり、展開の符号と大きさが合う。
:::

### 量子の時間応答を二準位の歳差で確かめる

$H_0=(\hbar\Omega/2)\sigma_z$ とし、場は $-\!f(t)\sigma_x$ に結合する。$\sigma_x(t)=\sigma_x\cos\Omega t-\sigma_y\sin\Omega t$ は $i[H_0,\sigma_x]/\hbar=-\Omega\sigma_y$ とその次の微分から得られる。従って $[\sigma_x(t),\sigma_x]=2i\sigma_z\sin\Omega t$ である。正準平均 $\langle\sigma_z\rangle_0=-\tanh(\beta\hbar\Omega/2)$ を入れると
$$
R_{xx}(t)=\frac2\hbar\Theta(t)\tanh(\beta\hbar\Omega/2)\sin\Omega t
$$
となる。零温度で基底状態から上の準位へ向かう遷移は残り、高温では上下の準位の重みが近づいて交換子の期待値は小さくなる。調和振動子の古典的な応答と似た正弦の時間依存でも、振幅は温度に依存する。

:::example[同じ準位差でも占有の差が半分になる温度]
$\beta\hbar\Omega=\log3$ なら $\tanh(\beta\hbar\Omega/2)=1/2$ である。時刻 $t=\pi/(2\Omega)$ の応答核は $1/\hbar$ となり、零温度の値 $2/\hbar$ の半分である。場の形や周波数を変えなくても温度だけで吸収強度が変わる。ここでの減衰のない正弦核の静的極限は、収束因子を付けた周波数変換から読む。
:::

## 複素感受率と Kramers–Kronig の関係

### 上半平面での解析性

実の因果核が可積分なら $\chi(z)=\int_0^\infty R(t)e^{izt}\dd t$ は $\Im z>0$ で解析的である。任意の内点近くで $\Im z$ の正の下限を取ると、$t e^{-t\Im z}$ が有界になり積分を微分できる。実軸では $\chi(-\omega)=\overline{\chi(\omega)}$ である。上半円で $\chi$ が十分に減衰し、境界値が滑らかと仮定する。Cauchy の積分公式を実軸上の点の直上で評価し、そこへの極限をとる。小半円の寄与は $i\pi\chi(\omega)$ であり、

$$
\chi'(\omega)=\frac1\pi\operatorname{PV}\int_{\R}\frac{\chi''(u)}{u-\omega}\,\dd u,\qquad
\chi''(\omega)=-\frac1\pi\operatorname{PV}\int_{\R}\frac{\chi'(u)}{u-\omega}\,\dd u
$$

が従う。ここで $\chi',\chi''$ は微分でなく実部と虚部である。この :term[Kramers–Kronig の関係]{en="Kramers–Kronig relations" yomi="くらまーすくろーにひのかんけい"} は、虚部の測定から実部を求めるのに使える。高周波で一定値に収束する場合は $\chi-\chi(\infty)$ に適用する。

$A=B$、$f(t)=f_0\cos\omega t$ なら仕事率は $-\dot f(t)\langle B(t)\rangle$ である。応答を $\langle B\rangle=f_0(\chi'\cos\omega t+\chi''\sin\omega t)$ と書くと周期平均は $\omega f_0^2\chi''/2$ である。従ってこの約束では正周波数の正の虚部が吸収を示す。

:::example[一つの緩和時間]
$R=a e^{-\lambda t}$ のとき $\chi=a/(\lambda-i\omega)$ で、極は $-i\lambda$ にある。$\omega=\lambda$ では実部と虚部がともに $a/(2\lambda)$ になる。
:::

:::example[静的な総和]
$\omega=0$ の分散関係と虚部の奇性から $\chi'(0)=2\int_0^\infty\chi''(u)\dd u/(\pi u)$ である。前例の右辺は $2a\int_0^\infty(\lambda^2+u^2)^{-1}\dd u/\pi=a/\lambda$ となる。
:::

### 静的分散関係と減衰の時間遅れ

実の因果応答では $\chi(-\omega)=\chi(\omega)^*$ である。零周波数の分散関係で正負を組にすると
$$
\chi(0)-\chi_\infty=\frac2\pi\int_0^\infty\frac{\Im\chi(\omega)}{\omega}\,\dd\omega
$$
を得る。ここで $\chi_\infty$ は実の瞬時項であり、引いた後の関数に上半平面での解析性と半円上の減衰を仮定する。従って吸収の周波数積分から静的応答を推定できる。$1/\omega$ の重みがあるため、観測窓から外れた遅い緩和に注意が要る。

:::example[同じ静的応答で異なる吸収帯]
$\chi-\chi_\infty=A/(1-i\omega\tau)$ とすると、虚部は $A\omega\tau/(1+\omega^2\tau^2)$ である。$0<\omega<\Omega$ の寄与は $(2A/\pi)\arctan(\Omega\tau)$ になる。$\Omega=1/\tau$ なら全体の半分である。$\tau$ を二倍にするとピークは半分の周波数へ移るが、全積分は $A$ のまま変わらない。
:::

### 高周波の項を短時間の核から得る

$\chi(\omega)=\int_0^\infty e^{i\omega t}R(t)\dd t$ とし、$R$ と必要な階数の導関数が無限遠で減衰し絶対可積分とする。積分を部分積分すると
$$
\chi(\omega)=-\frac{R(0^+)}{i\omega}
-\frac1{i\omega}\int_0^\infty e^{i\omega t}R'(t)\dd t
$$
である。さらに反復すれば、$R'(0^+)$ の項は $R'(0^+)/(i\omega)^2$、$R''(0^+)$ の項は $-R''(0^+)/(i\omega)^3$ になる。最後に残る Fourier 変換は、可積分なら Riemann–Lebesgue の補題により零へ行く。瞬時応答を含む場合の定数項は、この滑らかな核の展開とは別に先に足す。従って応答の開始の仕方が、高周波で最初に残るべきを決める。

:::example[一次の位置応答と慣性応答を比べる]
$R(t)=Ae^{-at}$ なら $R(0^+)=A$ であり、高い正の周波数で $\chi\sim iA/\omega$ である。一方、質量 $m$ の調和粒子の位置に力積を加えた直後は変位が零で、速度だけが $1/m$ 増える。従って位置の核は $R(0^+)=0$、$R'(0^+)=1/m$ であり、$\chi\sim-1/(m\omega^2)$ となる。低周波で同じ静的ばね定数をもつ二模型でも、高周波では先頭の次数から違う。
:::

### 分散関係で低周波側の欠測を見積もる

複数の指数緩和を $R(t)=\sum_j A_j \tau_j^{-1}e^{-t/\tau_j}$、$A_j>0$ とする。吸収の重み $\chi''(\omega)/\omega$ の各寄与は $A_j\tau_j/(1+\omega^2\tau_j^2)$ である。測定の下限を $\omega_{\min}$、上限を $\omega_{\max}$ とすれば、静的係数のうち範囲に入る部分は
$$
\sum_j\frac{2A_j}{\pi}
[\arctan(\omega_{\max}\tau_j)-\arctan(\omega_{\min}\tau_j)]
$$
である。各項の全積分は $A_j$ であるが、$1/\tau_j$ が測定域の外にあれば重みの大半が欠測になる。特に遅い緩和では、吸収のピークが測定下限のはるか下でも、静的係数には $A_j$ がそのまま入る。

:::example[長い緩和が応答の面積に隠れる]
同じ強度 $A_1=A_2=1$ で $\tau_1=1\,\mathrm s$、$\tau_2=100\,\mathrm s$ とし、$\omega_{\min}=0.1\,\mathrm{s^{-1}}$ から無限大まで測る。第一の重みは $1-(2/\pi)\arctan0.1\simeq0.937$、第二は $1-(2/\pi)\arctan10\simeq0.063$ である。合計の測定値はちょうど $1$ だが、真の静的係数は $2$ である。データの高周波側で精度を上げても、下限を変えなければこの不足は消えない。
:::

### 因果性を時間領域の変換から検査する

複素感受率の実部と虚部がそろっている場合、逆 Fourier 変換
$$
R(t)=\frac1{2\pi}\int_\R e^{-i\omega t}\chi(\omega)\dd\omega
$$
で時間核を復元できる。十分な減衰や収束因子を仮定する。因果的であれば $t<0$ は零である。例えば $\chi=1/(a-i\omega)$、$a>0$ は下半平面に極 $\omega=-ia$ を持つ。$t>0$ では $e^{-izt}$ が下半平面で減衰するので極を囲み、$R(t)=e^{-at}$ となる。$t<0$ では上半平面を使い極がないので零である。分母を複素共役した $1/(a+i\omega)$ は絶対値が同じだが、上半平面に極があり、対応する核は負の時間にある。振幅の絶対値だけを測って位相を捨てるとこの違いを判定できない。

:::example[純粋な遅延は絶対値を変えない]
因果核 $R(t)$ を $R_d(t)=R(t-d)$、$d>0$ として遅らせる。変数変換から $\chi_d(z)=e^{izd}\chi(z)$ である。実周波数上では絶対値は同じだが、位相は $\omega d$ だけ増える。上半平面では $|e^{izd}|=e^{-d\Im z}$ なので因果性と矛盾しない。例えば $d=1$ 秒、$\omega=\pi/2\,\mathrm{s^{-1}}$ なら位相差は $\pi/2$ である。同じ吸収のピークの幅を持つように見える振幅の記録でも、時間遅れは位相を測らなければ見えない。
:::

### 正則成分に分けて零周波のピークを読む

電流の核に一定の成分 $D_0$ があるとする。正の $\epsilon$ を掛けて変換すると、定数成分は $\chi_\epsilon=D_0/(\epsilon-i\omega)$ であり、
$$
\Re\chi_\epsilon=\frac{D_0\epsilon}{\epsilon^2+\omega^2},
\qquad\Im\chi_\epsilon=\frac{D_0\omega}{\epsilon^2+\omega^2}
$$
となる。前者は全周波数で積分すると $\pi D_0$ であり、幅が狭くなるにつれて $\pi D_0\delta(\omega)$ に収束する。後者を連続で滑らかな試験関数に掛けると、原点の近くで試験関数の偶部分が消え、奇部分は $\omega$ に比例する。そのため極限は $D_0\,\operatorname{PV}(1/\omega)$ である。点ごとの有限な値だけを残すと前者の面積を失う。滑らかな成分とは別にこの対を計上すれば、実部と虚部は同じ因果応答から導かれたままである。

:::example[分解能を二倍にしても面積は一定]
$\epsilon$ を半分にすれば、$\omega=0$ での実部の高さは二倍、半値半幅は半分になる。固定した $\Omega>0$ の窓に入る面積は $2D_0\arctan(\Omega/\epsilon)$ であり、$\epsilon\to0$ で全体の面積 $\pi D_0$ に近づく。逆に高さだけから散逸の積分強度を読むと、分解能を変えたときに異なる値を報告してしまう。
:::

## 揺動散逸定理（第一種と第二種）

### 応答と自発的な揺らぎ

$A=B$ を実観測量とし、古典平衡で自己共分散 $C$ とその正の半直線での導関数が可積分で、$C(t)\to0$ とする。:eqref[eq-neq-classical-kubo] を正の半直線で部分積分すると、境界項は実であり、虚部は $\beta\omega\int_0^\infty C(t)\cos(\omega t)\dd t$ である。:ref[thm-neq-wiener-khinchin] の両側規約から

$$
\chi''(\omega)=\frac{\beta\omega}{2}S_{AA}(\omega)
\label{eq-neq-classical-fdt}
$$

となる。これが :term[揺動散逸定理]{en="fluctuation–dissipation theorem" yomi="ようどうさんいつていり"} の第一種の形である。

第二種では、摩擦の記憶核とランダムな力の相関を結ぶ。例として平衡の線形模型 $m\dot v=-\kappa x-\int_0^\infty K(s)v(t-s)\dd s+\xi(t)+F(t)$、$\kappa>0$ を考える。安定で、平衡の応答と相関に孤立した非減衰モードがないとする。$\widehat K_+(\omega)$ を半直線 Fourier 変換すると、$\omega\ne0$ における速度移動度は $\mu=(\widehat K_+-i\omega m+i\kappa/\omega)^{-1}$ となる。平衡位置には第一種を適用でき、$\mu=-i\omega\chi_{xx}$、$S_{vv}=\omega^2 S_{xx}$ から $S_{vv}=2k_{\mathrm B}T\Re\mu$ を得る。線形方程式のスペクトルは $S_{vv}=|\mu|^2S_{\xi\xi}$ であり、慣性とばねは $1/\mu$ の虚部にしか入らないため

$$
S_{\xi\xi}(\omega)=2k_{\mathrm B}T\Re\widehat K_+(\omega),
\qquad C_{\xi\xi}(t)=k_{\mathrm B}T K(|t|).
$$

最後の式では通常の可積分な $K$ を仮定した。この力相関は $\kappa$ を含まないので、同じ浴で拘束を外した速度のスペクトルに使える。自由粒子は平衡位置分布を持たないため、その位置相関に第一種を直接適用する必要はない。白色極限では半直線上の摩擦の面積 $\gamma$ を保つと、両側の力の相関は $2k_{\mathrm B}T\gamma\delta(t)$ になる。

:::example[Brown 粒子の速度]
移動度の実部 $\gamma/(\gamma^2+\omega^2m^2)$ を用いると $S_{vv}=2k_{\mathrm B}T\gamma/(\gamma^2+\omega^2m^2)$ となる。指数相関を Fourier 変換した値に等しい。
:::

:::example[色付きの力]
$K(s)=(\gamma/\tau)e^{-s/\tau}$、$s>0$ なら $S_{\xi\xi}=2k_{\mathrm B}T\gamma/(1+\omega^2\tau^2)$ である。$\omega\tau\ll1$ でだけ白色雑音の強度に近い。
:::

### 量子相関の重み

量子平衡での古典の相関公式に対応する量を

$$
C^{\mathrm K}_{BA}(t)=\frac1\beta\int_0^\beta
\langle\delta B(-i\hbar\lambda)\delta A(t)\rangle_0\,\dd\lambda
$$

と定める。$B(-i\hbar\lambda)=e^{\lambda H_0}B e^{-\lambda H_0}$ である。$\partial_t\langle B(-i\hbar\lambda)A(t)\rangle=-\langle\dot B(-i\hbar\lambda)A(t)\rangle$ は、跡の巡回性と $\rho_0$ が $H_0$ と可換なことから従う。また $\partial_\lambda B(-i\hbar\lambda)=[H_0,B(-i\hbar\lambda)]=-i\hbar\dot B(-i\hbar\lambda)$ である。ここで $\dot B=(i/\hbar)[H_0,B]$ を用いた。$\lambda$ について積分し、$\langle B(-i\hbar\beta)A(t)\rangle_0=\langle A(t)B(0)\rangle_0$ を使うと

$$
-\beta\partial_t C^{\mathrm K}_{BA}(t)=\frac{i}{\hbar}\langle[A(t),B(0)]\rangle_0.
$$

従ってこの相関に対しては量子の応答も :eqref[eq-neq-classical-kubo] と同じ微分の形をとる。$[B,H_0]=0$ なら虚時間の積分は不要になるが、この場合応答の交換子の期待値も零である。保存量の等時刻分散があっても力学的に応答する成分とは限らない点は、古典の長時間極限と対応する。

### スペクトルを用いた一般化摩擦の検査

定常な速度が単位質量の式 $\dot v=-\int_{-\infty}^tK(t-s)v(s)\dd s+\eta+f$ に従うとする。絶対可積分な因果核の変換 $\widehat K(\omega)$ を用い、速度の移動度を $\mu(\omega)=[\widehat K(\omega)-i\omega]^{-1}$ とする。第二種の揺動散逸関係 $\langle\eta(t)\eta(0)\rangle=k_{\mathrm B}T K(|t|)$ のもとでは $S_\eta=2k_{\mathrm B}T\Re\widehat K$ である。従って
$$
S_v=|\mu|^2S_\eta=2k_{\mathrm B}T|\mu|^2\Re\widehat K
=2k_{\mathrm B}T\Re\mu.
$$
最後の等号は複素数の逆数の実部から得られる。平衡の測定で速度相関と強制応答の両方を測れば、記憶核を個別に復元せずにこの関係を検査できる。

:::example[時間尺度のある摩擦]
$K(t)=ab e^{-bt}$ とする。$\widehat K=ab/(b-i\omega)$ である。$a=b=\omega=1$ では $\widehat K=(1+i)/2$、$\mu=1+i$ であり、$S_v/(k_{\mathrm B}T)=2$ となる。同じ面積 $a=1$ の瞬時摩擦では $\mu=1/(1-i)=(1+i)/2$ なので、同じ周波数で $S_v/(k_{\mathrm B}T)=1$ になる。零周波数ではどちらも $S_v/(k_{\mathrm B}T)=2/a$ である。
:::

### KMS 条件を固有状態の重みから示す

有限次元で $C^>_{AB}(t)=\Tr(\rho_0 A(t)B)$、$C^<_{AB}(t)=\Tr(\rho_0 B A(t))$ とする。$C^>(t)$ の固有状態表示は $\sum_{nm}p_n A_{nm}B_{mn}e^{-i(E_m-E_n)t/\hbar}$ である。$t$ を $t-i\beta\hbar$ に移すと、指数に $e^{-\beta(E_m-E_n)}$ が掛かり、$p_n$ が $p_m$ に変わる。これは $C^<(t)$ と一致する。従って $C^>(t-i\beta\hbar)=C^<(t)$ であり、これが有限系での :term[KMS 条件]{en="KMS condition" yomi="けーえむえすじょうけん"} である。無限系では複素時刻の解析性や境界値の意味を指定する必要があるが、有限系は有限個の指数の和である。

両側変換を $e^{i\omega t}$ とすると、同じ遷移 $(n,m)$ の $S^<(\omega)$ は $S^>(\omega)$ の $e^{-\beta\hbar\omega}$ 倍である。$A=B$ が Hermite のとき $\chi''=(S^>-S^<)/(2\hbar)$、対称化したスペクトルは $(S^>+S^<)/2$ である。後者を前者で表す比 $(1+e^{-x})/(1-e^{-x})=\coth(x/2)$（$x=\beta\hbar\omega$）が量子の係数になる。零周波数に保存量のデルタ線がある場合、交換子でその成分は消える。式をデルタ線の割り算として零周波にそのまま延長せず、正則部分と保存成分を分けて考える。

:::example[量子の温度をスペクトル線の比で読む]
$A=B$ の遷移が $\Omega>0$ の線を作り、正周波数の $S^>$ の重みが負周波数の三倍とする。固有状態を交換した二つの重みの比から $\beta\hbar\Omega=\log3$ である。正周波数で対称化重みを交換子の差の半分で割ると $(1+1/3)/(1-1/3)=2$ である。両線の高さを比較する場合は同じ周波数分解能が必要であり、デルタ線のままなら高さの代わりに積分強度を使う。
:::

### 量子の力雑音と摩擦の実部

$A=B$ の正則な周波数成分について、直前の KMS 条件は
$$
S_{AA}^{\rm sym}(\omega)=\hbar\coth(\beta\hbar\omega/2)\,\chi''_{AA}(\omega)
$$
を与える。$S^{\rm sym}$ は反交換積の半分から作る両側スペクトルである。$\beta\hbar|\omega|\ll1$ では $\hbar\coth(\beta\hbar\omega/2)\sim2/(\beta\omega)$ なので :eqref[eq-neq-classical-fdt] に戻る。

量子の線形な浴で、位置の応答と雑音の両方が同じ線形伝達関数で結ばれ、上の平衡条件を満たす場合を考える。速度移動度は古典と同じ形であり、$S_{vv}^{\rm sym}=\omega^2 S_{xx}^{\rm sym}=\hbar\omega\coth(\beta\hbar\omega/2)\Re\mu$ となる。もう一方で $S_{vv}^{\rm sym}=|\mu|^2S_{\xi\xi}^{\rm sym}$ だから、
$$
S_{\xi\xi}^{\rm sym}(\omega)
=\hbar\omega\coth(\beta\hbar\omega/2)\Re\widehat K_+(\omega).
$$
これは演算子としての力の対称化相関である。交換子を落として、任意の古典確率過程と同じだと解釈することはできない。$T\downarrow0$ で正の周波数の係数は $\hbar\omega$ に近づく。従って零温度でも力の対称化相関は残る。

### 白色雑音の近似が許す周波数

古典浴の一般化摩擦核 $K(t)=K_0e^{-|t|/\tau}$ に対し、偶の力相関を $k_{\mathrm B}T K(t)$ とする。対称な雑音スペクトルは $2k_{\mathrm B}T K_0\tau/(1+\omega^2\tau^2)$ である。$K_0\tau=\gamma$ を固定して $\tau$ を小さくすると、有限の $\omega$ では定数 $2k_{\mathrm B}T\gamma$ に近づく。一方 $|\omega|\sim1/\tau$ の周波数も一緒に追うと係数の差は残る。量子浴なら、摩擦核を短くしても $\beta\hbar|\omega|$ が小さくなければ $\coth(\beta\hbar\omega/2)$ を古典の係数へ置き換えられない。同じ白色近似に見える操作でも、浴の緩和時間と熱エネルギーの条件の両方がある。

:::example[二つの独立な小さい量]
$\omega\tau=0.1$ なら古典の指数核の白色近似の差は $1-1/(1+0.01)\simeq0.0099$ である。一方、同じ周波数で $\beta\hbar\omega=2$ なら量子と古典の係数の比は $\coth(1)\simeq1.313$ である。記憶の補正が小さいという観測だけで、量子による雑音の補正も小さいとは判断できない。
:::

### 局所の相関から有色の浴を構成する

古典の正定値な指数相関 $\langle\eta(t)\eta(s)\rangle=C_0 e^{-\lambda|t-s|}$ を考える。これは $\dd\eta=-\lambda\eta\dd t+\sqrt{2\lambda C_0}\dd W$ を定常から始めれば実現できる。単位質量の変数 $v$ に対して第二種の揺動散逸を課すなら、過去の速度に掛かる核は $K(t)=C_0 e^{-\lambda t}/(k_{\mathrm B}T)$（$t\ge0$）である。従って同じ雑音分散でも、$\lambda$ を変えると核の面積 $C_0/(\lambda k_{\mathrm B}T)$ が変わる。摩擦の面積を保ちたい場合は $C_0\propto\lambda$ にする。短い記憶極限で瞬時の分散が発散しても、時間平均の雑音に有限な極限が残る。

:::example[幅が短くなる分を強度で補う]
$C_0=\gamma k_{\mathrm B}T\lambda$ とする。観測窓 $\Delta$ の雑音の積分の分散は
$$
2\gamma k_{\mathrm B}T\left[\Delta-\frac{1-e^{-\lambda\Delta}}\lambda\right]
$$
である。二重積分して相関の指数を時差で評価した。$\lambda\Delta=\log2$ では白色モデルの $2\gamma k_{\mathrm B}T\Delta$ に対する比が $1-1/(2\log2)\simeq0.279$ である。$\lambda\Delta\gg1$ になるまで、同じ面積を持つ雑音でも短い時間窓では揺らぎが異なる。
:::

### 相関と時間の規約から温度を照合する

両側の古典自己相関を $S(\omega)=\int_\R e^{i\omega t}C(t)\dd t$ とする。定常な中心化実過程では $C(t)$ は偶であり、$S=2\int_0^\infty C(t)\cos\omega t\,\dd t$ である。従って周波数分解した $\chi''$ と比較する際に、一側の余弦積分だけを $S$ と呼ぶ流儀では係数が二倍違う。時間領域でも、雑音の偶相関がデルタに近づくとき、片側の積分に入る摩擦の面積は全域の半分である。どちらの規約を使っても、次元をそろえ、同じ係数を応答と雑音の両方で使えば温度の判定は変わらない。

:::example[二つの測定から熱浴の温度を読む]
減衰率 $a=2\,\mathrm{s^{-1}}$、速度相関の等時刻値 $C(0)=3$ の指数相関を測ったとする。両側スペクトルは $S(\omega)=12/(4+\omega^2)$ である。速度の移動度が $\mu=1/[m(a-i\omega)]$ なら $2\Re\mu=4/[m(4+\omega^2)]$ なので、比 $S/(2\Re\mu)$ は $3m$ である。これは $m\langle v^2\rangle=k_{\mathrm B}T$ と一致する。もし測定した比が周波数で変われば、単一の平衡浴や同じ観測量の比較になっているかを調べる必要がある。
:::

## Green–Kubo 公式と輸送係数

### 保存流の輸送係数

場が位置のモーメント $B$ に共役で、$\dot B=J$ が全流なら、古典の :eqref[eq-neq-classical-kubo] で出力を $A=J$ と選ぶと、流に対する応答核は $R_{JB}(t)=\beta\Theta(t)\langle J(t)J(0)\rangle$ になる。外場は $B$ に共役で、時間微分 $\dot B=J$ によって流の自己相関が現れる。体積 $V$ で割り、定常場に対する電流密度の係数を読めば $\sigma=\beta V^{-1}\int_0^\infty\langle J_e(t)J_e(0)\rangle\dd t$ である。電荷のモーメントは周期箱では巻き戻さない座標で定義する。

同じ関係は積分流の分散からも求まる。$G(t)=\int_0^t J(s)\dd s$ の分散を :eqref[eq-neq-brownian-green-kubo] と同様に展開すると、$\lim_{t\to\infty}\operatorname{Var}G(t)/(2t)=\int_0^\infty C_J$ である。平衡で定義した全せん断応力 $J_{xy}$ とせん断速度に対し $\eta=(k_{\mathrm B}TV)^{-1}\int C_{xy}$、対流分を引いた全熱流 $J_q$ と温度勾配に対し $\kappa=(k_{\mathrm B}T^2V)^{-1}\int C_q$ となる。熱の力が $\nabla(1/T)=-\nabla T/T^2$ であるため、熱伝導では $T$ の因子が一つ増える。これらを :term[Green–Kubo 公式]{en="Green–Kubo formula" yomi="ぐりーんくぼこうしき"} という。相関の可積分性と、流のうち保存量に重なる成分を取り除く条件は個別に確かめる。

:::example[指数相関の面積]
$C_J(t)=C_0e^{-t/\tau}$ なら面積は $C_0\tau$ である。同じ等時刻の分散でも、緩和が二倍遅いと輸送係数は二倍になる。
:::

:::example[保存成分の分離]
$J(t)=a P+J_\perp(t)$ で $P$ は保存され、$J_\perp$ はそれに直交するとする。$C_J$ に $a^2\langle P^2\rangle$ が残り、積分した流の分散は $t^2$ の項をもつ。有限の拡散係数を求めるには対象となる流の定義からこの成分を分離する。
:::

### 積分運動量のモーメントで応力を得る

周期箱を考える前に、対力が働く有限粒子系で $G_{xy}=\sum_a x_a p_{ay}$ を微分する。対力が作用反作用 $F_{ab}=-F_{ba}$ を満たせば

$$
\dot G_{xy}=\sum_a\frac{p_{ax}p_{ay}}{m_a}
+\sum_{a<b}(x_a-x_b)F_{ab,y}
$$

である。第二項は $\sum_a x_a\dot p_{ay}$ の対ごとの組合せである。これが全運動量流のせん断成分であり、単位はエネルギーである。体積で割ると応力の単位になる。$G_{xy}$ に共役な場 $h$ の次元は時間の逆数である。応答の規約で $h$ をせん断速度に合わせれば、:eqref[eq-neq-classical-kubo] において $B=G_{xy}$、$A=J_{xy}$ として $\eta=\beta V^{-1}\int_0^\infty C_{xy}(t)\dd t$ を得る。境界を含めた運動量束と場の符号は、流をどちら向きに正としたかに合わせる。

境界のない並進一様系の輸送を得るには、有限の箱で時間無限大を先に取る操作に注意する。閉鎖箱の有界な位置モーメントでは、長時間に積分流の分散が飽和し得る。周期箱では位置の巻き戻りと力を周期像に取る補正を含め、$G(t)-G(0)$ を測定した全流の積分として構成する。相関表現を使う場合も、輸送を測る時刻範囲で境界や保存モードの影響を確かめる。

### 有限時間で測った輸送係数

平均零の定常な実電流 $J$ に対し、$G_t=\int_0^t J(s)\dd s$ とする。二重積分を時差で整理すると
$$
\frac{\operatorname{Var}G_t}{2t}
=\int_0^t (1-s/t)C_J(s)\dd s.
$$
これは相関を時刻 $t$ まで単純に積分した量と異なり、端に向かって線形に重みが落ちる。$C_J\ge0$ なら右辺は単純な打切り積分より小さい。一方、負の相関も含む電流では差の符号は相関の一階時間モーメントに依存する。

:::example[単一緩和の有限窓補正]
$C_J(s)=C_0e^{-s/\tau}$ なら比は $C_0\tau[1-(\tau/t)(1-e^{-t/\tau})]$ となる。$t=\tau$ では無限時間の輸送積分の $1/e$ である。その時点の単純な打切り積分は $1-1/e$ 倍となる。前者は区間ごとの積分量から得られる値、後者は各遅れの相関を推定してから積分した値である。
:::

### 局所平衡の準備から熱伝導率を導く

時間反転対称な単一成分の流体を考える。平衡の一粒子当たりのエンタルピーを $h_0$、内部エネルギーと粒子数の密度をそれぞれ $\mathcal E(x)$、$\mathcal N(x)$ とし、$\mathcal Q=\mathcal E-h_0\mathcal N$ を用いる。対応する全熱流を $J_q$ とする。一方向について $B_q=\int x\,\mathcal Q(x)\dd x$ と書けば、境界項を除いた連続式は $\dot B_q=J_q$ を与える。有限の試料に対してこの恒等式を使うか、境界の効果を避けるために長波長の摂動の極限として扱う。

小さい逆温度勾配 $g=\partial_x\beta$ を持つ局所平衡から出発する。エントロピーに共役な係数を位置に沿って変えると、正準重みには一次で $-g(B_q-\langle B_q\rangle_0)$ が掛かる。粒子数の係数も同時に選び、エンタルピーを運ぶ平均流を除いた摂動にしている。そこから無摂動の力学で発展させれば $\delta\langle J_q(t)\rangle=-g\langle J_q(t)\delta B_q(0)\rangle_0$ となる。等時刻の相関は時間反転で奇なので零であり、相関の時間微分は定常性と $\dot B_q=J_q$ によって $-\langle J_q(t)J_q(0)\rangle_0$ である。従って熱流密度は
$$
\frac{\delta\langle J_q(t)\rangle}{V}
=\frac{\partial_x(1/T)}{k_{\mathrm B}V}
\int_0^t\langle J_q(s)J_q(0)\rangle_0\,\dd s .
$$
相関が減衰して積分が収束し、巨視的な温度分布がまだ大きく変わらない時間窓を取れば、$j_q=-\kappa\partial_x T$ により $\kappa=(k_{\mathrm B}T^2 V)^{-1}\int_0^\infty\langle J_q(t)J_q(0)\rangle_0\dd t$ となる。温度勾配が消えるまで待つ実験では、孤立した試料の熱流は最終的に零へ戻る。

### 粒子流から熱流へ基準を取り替える

エネルギー流 $J_E$ と粒子流 $J_N$ がともにあるとする。熱として測りたい流れを $J_Q=J_E-h_0J_N$ とし、$h_0$ を平衡の一粒子当たりのエンタルピーとする。積分した流の共分散を使う Green–Kubo の形では、基準を変えた熱への係数は
$$
L_{QQ}=L_{EE}-h_0(L_{EN}+L_{NE})+h_0^2L_{NN}
$$
となる。右辺は $\operatorname{Var}(G_E-h_0G_N)/(2t)$ の極限を展開したものである。積分した流の共分散で定められるのは輸送行列の対称部分であり、磁場中の Hall 応答のような反対称部分には片側の応答積分が要る。このような線形変換では、同じ発熱・散逸の収支を記述するために共役な力も逆転置で変換する。実験で粒子の移流を取り除いた熱輸送を求めるには、基準として引くエンタルピーと、粒子流を零にする境界条件を区別する。前者は観測量の変換であり、後者は力を調整する実験条件である。

:::example[移流の寄与を引くと相関が変わる]
規格化した係数を $L_{EE}=9,L_{EN}=L_{NE}=3,L_{NN}=2$、$h_0=1$ とする。変換後は $L_{QQ}=9-6+2=5$、$L_{QN}=3-2=1$ である。さらに粒子流が零という条件で自由な力を消去すると、熱の係数は $5-1^2/2=9/2$ になる。$5$ と $9/2$ の差は、熱と粒子の相関が変換後にも残ることから出る。
:::

### 保存流の平均が残ると輸送積分が発散する

平衡の電流 $J$ が保存量 $P$ に成分を持つ場合を考える。中心化し、$J=aP+J_\perp$、$a=\langle JP\rangle/\langle P^2\rangle$、$\langle J_\perp P\rangle=0$ とする。$P(t)=P(0)$ なので交差相関はどちらも時刻にかかわらず零であり、
$$
C_J(t)=a^2\langle P^2\rangle+C_{J_\perp}(t)
$$
である。$G_t=\int_0^tJ(s)\dd s$ には $atP$ が入るため、分散には $a^2t^2\langle P^2\rangle$ が残る。輸送係数を有限の拡散として測るには、この保存成分を固定するか、物理的に測る流をその成分と直交させる必要がある。その条件を課した後にも、残りの相関が可積分かどうかを調べる。

:::example[箱の重心移動が粒子流に含まれる場合]
質量 $m$ が同じ粒子では、全粒子流 $\sum_i v_i$ は全運動量 $P$ を $m$ で割った量である。外力のない周期箱で $P$ が保存され、その分散が正なら相関は定数である。$\operatorname{Var}(G_t)/(2t)=t\langle P^2\rangle/(2m^2)$ は時間に比例して増える。重心運動を引くことは、この例で粒子流全体を引くことに相当するが、熱や相対的な成分の流れでは引いた後に非零の輸送が残る。
:::

### 密度の応答で波長と時間の極限を比べる

単独で拡散する保存密度 $n$ を考える。例えば一様な背景を固定して、標識粒子を追う場合である。外場 $f(x,t)$ は $-\int f(x,t)n(x)\dd x$ で結合し、一様な基準状態の等温感受率を $\chi_n=\partial n/\partial\mu_{\rm chem}>0$ とする。ここで $\mu_{\rm chem}$ は粒子一個当たりの化学ポテンシャルである。局所平衡からのずれでは $\delta\mu_{\rm chem}=\delta n/\chi_n-f$ であり、流れ $j=-D\chi_n\nabla\delta\mu_{\rm chem}$ と連続式を組み合わせると
$$
\partial_t\delta n=D\nabla^2(\delta n-\chi_n f)
$$
となる。従って波数 $k$ と角周波数 $\omega$ の応答は
$$
\chi(k,\omega)=\frac{\chi_n Dk^2}{Dk^2-i\omega}.
$$
$k\ne0$ を固定して $\omega\to0$ とすれば、平衡の密度変化 $\delta n=\chi_n f$ が得られる。逆に $\omega\ne0$ を固定して $k\to0$ とすると、空間一様な場では粒子が再配分されないため応答は零に近づく。ここでは閉じた系の調和成分を考えており、一様モードの粒子数を外部と交換させる条件では別の応答を測る。

:::example[波数を半分にすると応答のピークが移る]
虚部は $\chi_n Dk^2\omega/[(Dk^2)^2+\omega^2]$ であり、正の周波数 $\omega=Dk^2$ で高さ $\chi_n/2$ を取る。波数を半分にすればピーク周波数は四分の一になる。静的な変化の大きさが同じでも、粒子が再配分される時間は長くなる。
:::

### 電流の緩和から Drude の伝導を求める

希薄な担体の電荷を $q$、数密度を $n$、質量を $m$ とし、速度が率 $1/\tau$ で緩和するとする。一方向の各平均速度は $\dot v=-v/\tau+qE(t)/m$ であり、電流密度は $j=nqv$ である。従って一定な周波数の伝導率は
$$
\sigma(\omega)=\frac{nq^2\tau/m}{1-i\omega\tau}.
$$
平衡において各担体の速度相関を $(k_{\mathrm B}T/m)e^{-|t|/\tau}$ とし、担体間の交差相関を無視する。体積 $V$ で全電流を $J=q\sum_i v_i$ とすれば $\langle J(t)J(0)\rangle=Vnq^2(k_{\mathrm B}T/m)e^{-|t|/\tau}$ になる。これを $\beta/V$ 倍した片側の Fourier 積分へ入れると同じ $\sigma(\omega)$ を得る。電場による決定論的な計算と平衡のゆらぎをつなぐ際、全電流で始めるなら体積で割る必要がある。

:::example[違う緩和時間で同じ重みを持つ伝導]
$\sigma(0)$ は $\tau$ に比例するが、$\Re\sigma$ の正周波数側の面積は $(\pi/2)nq^2/m$ であり、$\tau$ に依存しない。$\tau$ を四倍にするとピークの高さは四倍、半値半幅は四分の一である。例えば $\omega\tau=1$ では実部が静的係数の半分、虚部の大きさも同じになる。$\tau\to\infty$ で零周波に面積が集中し、時間領域の緩和が消える極限に対応する。
:::

### 流の時間積分に全微分を加えた場合

輸送係数を $\lim_{t\to\infty}\operatorname{Var}G_t/(2t)$（必要な熱力学的係数を除いた形）で測る。観測量の定義を $G_t'=G_t+U_t-U_0$ と変え、$U_t$ は定常で二次モーメントが有界であるとする。差 $H_t=U_t-U_0$ の分散は $4\operatorname{Var}U_0$ 以下である。従って $\operatorname{Var}G_t=O(t)$ なら $\operatorname{Cov}(G_t,H_t)=O(\sqrt t)$ であり、
$$
\frac{\operatorname{Var}G_t'-\operatorname{Var}G_t}{2t}
=\frac{\operatorname{Var}H_t+2\operatorname{Cov}(G_t,H_t)}{2t}\to0.
$$
これは局所のエネルギーの分け方を変えて流に全微分が付いたときの検算に使える。加えた量が無制限に拡散する、または分散が体積・時間の極限に対して一様に抑えられない場合には、この証明はそのまま適用できない。

:::example[有限時間でだけ見える端点のゆらぎ]
独立な $G_t$ と端点差の組で、$\operatorname{Var}G_t=2Dt$、$U_t$ の分散を $s^2$、相関を $s^2e^{-at}$ とする。$\operatorname{Var}G_t'/(2t)=D+s^2(1-e^{-at})/t$ である。$D=s^2=a=1$ なら $t=1$ で約 $1.632$、$t=10$ で約 $1.100$ となる。長時間の値が同じであっても、計算で使う時間が有限なら端点の揺らぎの大きさは観測できる。
:::

## 総和則と Onsager–Casimir の対称性

### 短時間と時間反転

因果核 $R$ とその導関数が十分に減衰するとき、部分積分を二回行うと

$$
\chi(z)=-\frac{R(0^+)}{iz}+\frac{R'(0^+)}{(iz)^2}+\cdots
$$

である。端点値は応答の高周波を支配する。位置の応答では $R(0^+)=0$、$R'(0^+)=1/m$ である。分散関係を上半平面で書いた $\chi(i y)=(2/\pi)\int_0^\infty u\chi''(u)/(u^2+y^2)\dd u$ に $y^2$ を掛け、$y\to\infty$ とすると $\int_0^\infty u\chi''(u)\dd u=\pi/(2m)$ を得る。ここでは $u\chi''$ が可積分な場合を仮定した。この種の積分制約が総和則である。

古典で時間反転 $\Theta$ に対する $A,B$ の符号を $\epsilon_A,\epsilon_B$ とし、外部磁場を $\mathcal B$ とする。平衡分布が $\Theta$ により $-\mathcal B$ の分布へ移り、軌道が逆にたどられるなら

$$
\langle A(t)B(0)\rangle_{\mathcal B}
=\epsilon_A\epsilon_B\langle B(t)A(0)\rangle_{-\mathcal B}.
$$

変数変換で両観測量に符号を付け、定常性で $-t$ を相手側の $t$ に移した。各相関が可積分で、共役な流と力を同じ規約で選んでいる場合、これを流の積分に適用すると :term[Onsager–Casimir の対称性]{en="Onsager–Casimir reciprocity" yomi="おんさがーかしみーるのたいしょうせい"} $L_{ij}(\mathcal B)=\epsilon_i\epsilon_j L_{ji}(-\mathcal B)$ を得る。流を生成する二つの変数の符号を使っても、微分で両方の符号が変わるので積は同じである。

:::example[磁場反転と Hall 成分]
二つの電流はともに時間反転で奇である。従って $\sigma_{xy}(\mathcal B)=\sigma_{yx}(-\mathcal B)$ である。これらが反対称な成分を作る場合、その成分は磁場に関して奇である。
:::

:::example[速さの尺度]
質量 $m=2\,\mathrm{kg}$ の Drude 型の拘束粒子で $\chi_{xx}\sim-1/(m\omega^2)$ が成立すれば、正周波数側の $\omega\chi''$ の積分は $\pi/4$ に対応する。摩擦を変えると吸収の分布は変わるが、この重み付き面積は質量で決まる。
:::

### 二重交換子による吸収の一次モーメント

離散準位をもち、以下の絶対収束が成立する $A=B$ の応答を考える。正負を含むスペクトルを固有状態で展開すると、$\int \omega\chi''(\omega)\dd\omega/\pi$ は $\sum_{n,m}(p_n-p_m)(E_m-E_n)|A_{nm}|^2/\hbar^2$ になる。添字交換で $p_m$ の項を $p_n$ の項へまとめれば二倍となる。これを演算子で書くと
$$
\int_{\R}\frac{\omega\chi''(\omega)}\pi\,\dd\omega
=\frac{\langle[A,[H_0,A]]\rangle_0}{\hbar^2}.
$$
実際、二重交換子の $n$ での対角成分は $2\sum_m(E_m-E_n)|A_{nm}|^2$ である。正の周波数だけの積分はこの半分である。

:::example[位置の総和則]
$H_0=p^2/(2m)+U(x)$、$A=x$ とする。$[H_0,x]=-i\hbar p/m$、さらに $[x,[H_0,x]]=\hbar^2/m$ である。従って両側の一次モーメントは $1/m$ となる。$U$ によって吸収線の位置は変わるが、重みをすべて足した値は変わらない。
:::

### 磁場反転を相関の時刻の組で表す

古典平衡で、場 $\mathcal B$ における軌道を時間反転写像 $\Theta$ で反転すると場 $-\mathcal B$ の軌道になるとする。平衡分布も $\Theta$ で写り、観測量 $A_i(\Theta z)=\epsilon_i A_i(z)$ とする。$C_{ij}(t;\mathcal B)=\langle A_i(t) A_j(0)\rangle_{\mathcal B}$ を初期点での積分として書き、積分変数を $z=\Theta z'$ と変える。反転した軌道上では第一の観測が時刻 $-t$ なので
$$
C_{ij}(t;\mathcal B)=\epsilon_i\epsilon_j C_{ij}(-t;-\mathcal B)
=\epsilon_i\epsilon_j C_{ji}(t;-\mathcal B)
$$
である。最後では定常性で時刻を $t$ だけ移し、古典変数の積の可換性を使った。相関の時間積分が有限なら、同じ等式が輸送係数に移る。量子では積の交換を同じようにはできないため、Kubo の内積や交換子の順序を保って時間反転の作用を追う。

:::example[横の流れは逆場で符号を変える]
同じ時間反転符号の二成分に対し、ある磁場で $L=\begin{pmatrix}2&1\\-1&2\end{pmatrix}$ とする。逆場の $L$ は転置であり、横の成分の符号が反転する。力 $(1,0)^\mathsf T$ を加えると横の流れは $-1$ と $+1$ になるが、生成 $X^\mathsf T LX$ は両方で $2$ である。横応答の符号の変化を散逸の変化と読み違えると、相反関係の内容を見失う。
:::

### 密度の総和則で質量を検査する

質量 $m$ の同種粒子が $H_0=\sum_i p_i^2/(2m)+U(x_1,\ldots,x_N)$ に従うとする。実で滑らかな関数 $a(x)$ に対して $A=\sum_i a(x_i)$ とすると、位置だけの $U$ は $A$ と可換である。$[p_{i\alpha},a(x_i)]=-i\hbar\partial_\alpha a(x_i)$ を二回使うと
$$
[A,[H_0,A]]=\frac{\hbar^2}{m}\sum_i|\nabla a(x_i)|^2
$$
となる。具体的には $[p^2,a]=p[p,a]+[p,a]p$ を使い、外側の $A$ との交換では $[A,p]=i\hbar\nabla a$ だけが残る。よって吸収の両側の一次モーメントは $m^{-1}\sum_i\langle|\nabla a|^2\rangle_0$ である。ポテンシャルは交換子に直接入らないが、密度が不均一なら期待値の重みに入る。

:::example[均一な媒質の一つの密度波]
$a(x)=\cos(kx)$、領域はその周期の整数倍、平衡密度は一様とする。$\langle\sin^2 kx_i\rangle=1/2$ なので、両側の一次モーメントは $Nk^2/(2m)$ である。波数を二倍にすればこのモーメントは四倍になる。$k=0$ では $A=N$ が保存され、有限周波数の吸収は消える。構造の違いで個々の線が移動しても、前提を満たす範囲ではその重み付きの和を比較できる。
:::

### 量子の相反性を反ユニタリ写像で確かめる

磁場 $\mathcal B$ と $-\mathcal B$ の有限系を結ぶ反ユニタリ作用素を $\Theta$ とし、$\Theta H(\mathcal B)\Theta^{-1}=H(-\mathcal B)$ とする。観測量 $A_i$ は Hermite で $\Theta A_i(\mathcal B)\Theta^{-1}=\epsilon_i A_i(-\mathcal B)$ とする。場 $\mathcal B$ での固有ベクトル $\ket n$ を $\Theta\ket n$ へ写せばエネルギーと Boltzmann 重みは同じである。反ユニタリ性から $A_{i,nm}(-\mathcal B)=\epsilon_i A_{i,nm}(\mathcal B)^*$ となる。従って二つの行列要素の積は
$$
A_{i,nm}(\mathcal B)A_{j,mn}(\mathcal B)
=\epsilon_i\epsilon_j A_{j,nm}(-\mathcal B)A_{i,mn}(-\mathcal B)
$$
となる。複素共役の後に各演算子の Hermite 性を使った。固有状態の応答の分母と重みは同じなので、$\chi_{ij}(\omega;\mathcal B)=\epsilon_i\epsilon_j\chi_{ji}(\omega;-\mathcal B)$ となる。瞬時項を含む場合は、その観測量の変化も同じ対称性を満たすか確認する。

:::example[反対の偶奇を持つ交差応答]
磁場を零にし、$A=x$ を偶、$B=p$ を奇とする。調和振動子では $x(t)=x\cos\Omega t+[p/(m\Omega)]\sin\Omega t$、$p(t)=p\cos\Omega t-m\Omega x\sin\Omega t$ である。従って $[x(t),p]=i\hbar\cos\Omega t$、$[p(t),x]=-i\hbar\cos\Omega t$ となる。久保公式の係数を掛けると、$R_{xp}=-\cos\Omega t$、$R_{px}=+\cos\Omega t$（$t>0$）であり、偶奇の積が $-1$ であることに対応する。
:::

## 演習問題

::::exercise{level=basic}
$R(t)=3e^{-2t}$、$t>0$ の静的感受率を求めよ。
:::solution
$\chi(0)=\int_0^\infty 3e^{-2t}\dd t=3/2$ である。段差応答の長時間の増分を場で割っても同じ値になる。
:::
::::
::::exercise{level=basic}
$R(t)=e^{-t}$ の $\omega=1$ における実部と虚部を求めよ。
:::solution
$\chi(1)=(1-i)^{-1}=(1+i)/2$ なのでともに $1/2$ である。
:::
::::
::::exercise{level=basic}
$S_{AA}(\omega)=6$、$\beta=2$、$\omega=3$ のとき古典の $\chi''$ を求めよ。
:::solution
:eqref[eq-neq-classical-fdt] より $\beta\omega S/2=18$ である。
:::
::::
::::exercise{level=basic}
$\chi''=2$、$\omega=3$、$f_0=1$ の平均仕事率を求めよ。
:::solution
$-\dot f=\omega f_0\sin\omega t$ と応答の正弦成分を掛け、$\sin^2$ を平均すると $\omega f_0^2\chi''/2=3$ である。
:::
::::
::::exercise{level=advanced}
$\langle\delta A(t)\delta B(0)\rangle=C_\infty+c e^{-t/\tau}$ とする。静的な力学的応答を求めよ。
:::solution
:eqref[eq-neq-classical-kubo] で $R=\beta c e^{-t/\tau}/\tau$ である。その積分は $\beta c$ で、等時刻の共分散から長時間に残る $C_\infty$ を引いた量に比例する。
:::
::::
::::exercise{level=advanced}
速度が $m\dot v+\gamma v+\kappa x=F$ に従う。$e^{-i\omega t}$ 規約の位置感受率と吸収の符号を求めよ。
:::solution
$v=-i\omega x$ だから $\chi=(\kappa-m\omega^2-i\gamma\omega)^{-1}$ である。虚部は $\gamma\omega/[(\kappa-m\omega^2)^2+\gamma^2\omega^2]$ であり、$\gamma>0,\omega>0$ で正となる。
:::
::::
::::exercise{level=advanced}
二つの流の対称係数を $L_{ij}=\lim_{t\to\infty}\operatorname{Cov}(G_i,G_j)/(2t)$ と定義する。この行列の非負性を証明せよ。
:::solution
任意の実ベクトル $a$ について $a^\top L a=\lim_{t\to\infty}\operatorname{Var}(\sum_i a_iG_i)/(2t)\ge0$ である。極限は有限個の成分で取っているので和と交換できる。
:::
::::
::::exercise{level=advanced}
磁場零で反対の時間反転符号をもつ二つの変数に対し、交差輸送の寄与が散逸から消えることを示せ。
:::solution
$L_{12}=-L_{21}$ なので、力 $X_1,X_2$ を掛けた二次形式の交差項は $(L_{12}+L_{21})X_1X_2=0$ となる。個々の交差流が零である必要はない。
:::
::::

## 確認テスト

::::quiz
$e^{-i\omega t}$ の外場に対し、正周波数での吸収を表す量を選べ。

- [ ] $\chi'$ の符号
- [x] $\chi''$ の符号
- [ ] $|\chi|$ の符号

:::explanation
平均仕事率は $-\overline{\dot f\,\delta B}=\omega f_0^2\chi''/2$ である。$\omega>0$ なら吸収の符号は $\chi''$ で決まる。
:::
::::
::::quiz{answer="2" tol="1e-6"}
$R(t)=4e^{-2t}$ の静的感受率を答えよ。

:::explanation
時間積分で $4/2=2$ である。ここでは瞬時の応答はない。
:::
::::
::::quiz{answer="3" tol="1e-6"}
$C_J(t)=6e^{-2|t|}$ の正の時間側の積分を答えよ。

:::explanation
$\int_0^\infty 6e^{-2t}\dd t=3$ である。半直線を用いるため両側積分の半分になる。
:::
::::
::::quiz
磁場を反転して比較する Onsager–Casimir の式で、二つの電流がともに時間反転で奇の場合を選べ。

- [x] $L_{ij}(\mathcal B)=L_{ji}(-\mathcal B)$
- [ ] $L_{ij}(\mathcal B)=-L_{ji}(-\mathcal B)$
- [ ] 各係数が必ず零である。

:::explanation
二つの符号の積が $(-1)^2=1$ である。時間反転後の相関で変数の順を交換して積分すると最初の式になる。
:::
::::

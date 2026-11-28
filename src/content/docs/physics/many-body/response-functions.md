---
title: "線形応答と相関関数"
description: "外場に対する線形応答を相関関数から求め、伝導度の総和則と Ward 恒等式を導く。"
notation:
  - symbol: '\chi^R_{AB}'
    meaning: 演算子 A の演算子 B に対する遅延応答
  - symbol: 'S(\vb q,\omega)'
    meaning: 単位体積あたりの動的構造因子
  - symbol: 'D_{ij}'
    meaning: 反磁性接触項
  - symbol: '\Gamma^\mu'
    meaning: 粒子数電流の完全頂点
---

## この章の目標

- 外場と観測量の符号をそろえて久保公式を導く。
- 虚時間相関と動的構造因子から遅延応答を計算する。
- 電気伝導度に含まれる接触項と光学総和則を求める。
- 連続の式から頂点の Ward 恒等式を導く。

## 久保公式と遅延相関関数

### 外場を入れたときの時間発展

[線形応答理論](/physics/nonequilibrium/linear-response/) の演算子公式を、電子の密度と電流に適用する。弱い場 $f(t)$ を $H'(t)=f(t)B$ として加え、$B$ は Hermite とする。相互作用表示の発展演算子は $1-i\int_{-\infty}^t\dd t'\,f(t')B(t')+O(f^2)$ である。その随伴を $A(t)$ の左に、発展演算子を右に置いて一次まで展開する。得られる :term[久保公式]{en="Kubo formula" yomi="くぼこうしき"} は
$$
\delta\langle A(t)\rangle=\int\dd t'\,\chi^R_{AB}(t-t')f(t'),\qquad
\chi^R_{AB}(t)=-i\theta(t)\langle[A(t),B(0)]\rangle .
\label{eq-mb-kubo}
$$
ここで初期状態は無摂動 Hamiltonian と可換である。摂動を $-fB$ と書く規約なら右辺の符号が変わる。測る演算子自体が外場に依存するときは、その微分で得られる接触項も足す。

:::example[二準位の静的応答]
$H_0=\operatorname{diag}(0,\Delta)$、$A=B=\sigma_x$、基底状態を取る。$\langle[A(t),A]\rangle=e^{-i\Delta t}-e^{i\Delta t}$ なので $\chi^R(z)=(z-\Delta)^{-1}-(z+\Delta)^{-1}$ である。$\Delta=2\,\mathrm{eV}$ なら $\chi^R(0)=-1\,\mathrm{eV^{-1}}$ である。基底エネルギーの二次補正 $-f^2/\Delta$ を $f$ で微分しても同じ応答を得る。
:::

:::example[保存量への場]
$[Q,H_0]=0$ なら $\chi^R_{QQ}=0$ である。熱浴と再平衡した系では $\partial_f\langle Q\rangle_{H_0+fQ}=-\beta\operatorname{Var}(Q)$ であり、閉じた系の動的応答と測り方が異なる。
:::

### 二準位で Kubo 公式の符号を決める

$H_0=\Delta\ketbra e e$、$A=B=\ketbra e g+\ketbra g e$ とし、基底のエネルギーを零に取る。摂動 $+uB$ に対する $\chi^R(t)=-i\theta(t)\langle[A(t),B]\rangle$ は $-2\theta(t)\sin(\Delta t)$ である。実際、二つの積の基底期待値はそれぞれ $e^{-i\Delta t}$ と $e^{+i\Delta t}$ である。これを Fourier 変換すると
$$
\chi^R(z)=\frac1{z-\Delta}-\frac1{z+\Delta},\qquad \Im z>0 .
$$
静的極限は $-2/\Delta$ である。一方、厳密に $H_0+uB$ を対角化すると基底エネルギーは $(\Delta-\sqrt{\Delta^2+4u^2})/2$ である。一階微分は基底の $\langle B\rangle$ に等しく、その $u$ の一次は $-2u/\Delta$ である。従って場の符号と静的応答が一致する。

:::example[急に外場を加える場合]
$\Delta=2\,\mathrm{meV}$ の系へ時間零に弱い一定の $u_0$ を加える。因果核の積分は $\delta\langle A(t)\rangle=-2u_0(1-\cos\Delta t)/\Delta$ となる。$u_0=0.1\,\mathrm{meV}$ なら長時間平均は $-0.1$、最初の最小は $-0.2$ である。散逸がないため定常値へ単調に近づかず、外場を断熱的に加えた結果とは振動成分が異なる。
:::

### 実時間の源が二時刻の核を生む

$H(t)=H_0+\lambda(t)B$ とし、初期密度行列 $\rho_0$ は $H_0$ と可換であるとする。相互作用表示で
$$
U_I(t,t_0)=I-i\int_{t_0}^t\lambda(t')B_I(t')\dd t'+O(\lambda^2)
$$
を使うと、観測量 $A$ の平均は $\Tr[\rho_0U_I^\dagger A_I(t)U_I]$ である。左の $U_I^\dagger$ からは $+i\lambda B_I A_I$、右の $U_I$ からは $-i\lambda A_I B_I$ が出るので、一次は $-i\int\lambda(t')\langle[A_I(t),B_I(t')]\rangle_0\dd t'$ となる。時刻が $t'>t$ の摂動は発展に入らず、遅延関数の $\theta(t-t')$ はこの積分領域から生じる。初期状態が定常なら両方の時刻の平行移動をトレースの循環性で消せるため、核は時刻差だけで決まる。非定常な初期状態でも一次の二時刻式は成立するが、時刻差だけへ縮約できない。

正の源の代わりに Hamiltonian へ $-\lambda(t)B$ を入れた場合、同じ $-i\theta\langle[A,B]\rangle$ を定義すれば応答は $-\chi^R\lambda$ となる。例えば結合 $+un$ では正の $u$ が占有を減らすので密度感受率が負であり、結合 $-hS^z$ では正の $h$ が上向きの磁化を増やす。どの量を外場と呼ぶかを式の中で確認すると、二つの教科書で異なる符号規約にも対応できる。

:::example[二準位で交差の遅延成分を調べる]
$H_0=(\Delta/2)\sigma_z$、$A=\sigma_y$、$B=\sigma_x$ とし、$\sigma_z=-1$ の基底を取る。Heisenberg 表示で $\sigma_y(t)=\sigma_y\cos\Delta t+\sigma_x\sin\Delta t$ である。従って $[\sigma_y(t),\sigma_x]=-2i\sigma_z\cos\Delta t$、その基底期待値は $2i\cos\Delta t$ である。$+\lambda\sigma_x$ に対する $\chi_{yx}^R(t)$ は $2\theta(t)\cos\Delta t$ となり、$\chi_{xy}^R(t)$ はその負になる。遅延関数を Fourier 変換すると $\chi_{yx}^R(z)=i[(z-\Delta)^{-1}+(z+\Delta)^{-1}]$ である。静的に $\lambda$ を置いた基底は実の固有ベクトルを持つので $\langle\sigma_y\rangle$ は零であり、$\chi_{yx}^R(0)=0$ と整合する。交差成分では正の周波数の虚部だけから吸収の符号を判断せず、実際に加える場の組を指定する。
:::

### 周期駆動で散逸するエネルギーの符号

摂動 $u(t)B$ を $u(t)=\Re[u_0e^{-i\omega t}]$ とする。時間に依存する全 Hamiltonian のエネルギーの変化は $\langle\partial_tH\rangle=\dot u\langle B\rangle$ である。これは密度行列の発展で出る $[H,\rho]$ のトレースが循環性で零となるためである。定常な線形応答では $\delta B(t)=\Re[\chi_{BB}^R(\omega)u_0e^{-i\omega t}]$ なので、一周期の積の平均は
$$
\overline P=\frac12\Re[(-i\omega u_0)\overline{\chi_{BB}^R(\omega)u_0}]
=-\frac{\omega |u_0|^2}{2}\Im\chi_{BB}^R(\omega).
$$
従って正の $\omega$ で熱平衡からエネルギーを吸収するためには $\Im\chi_{BB}^R(\omega)\le0$ である。これは Lehmann 表示の各正周波数の重み $p_r-p_s\ge0$ と一致する。離散系の厳密な周期駆動で共鳴に長時間置く場合、線形応答を無制限に用いることはできない。分布の極限として吸収線を解釈するか、有限時間の分解能を持つ源で積分する。

:::example[同じ振幅で異なる位相の応答]
$\chi^R=-2-i/2$ をある正の周波数 $\omega=4$ で取る。$|u_0|=1/10$ なら一周期あたりの平均仕事率は $4(1/100)(1/2)/2=1/100$ である。実部の $-2$ は瞬間の観測量を変えるが、周期平均の仕事には現れない。純粋に実の $\chi^R$ に置き換えると、場と応答の積には振動が残っても一周期の仕事は零になる。極の狭い分光線を測るときは位相遅れの成分が直接エネルギー吸収を決める。
:::

### 静的な Kubo の値を基底の摂動と照合する

離散の非縮退基底 $\ket0$、励起 $\ket m$ のギャップ $\Delta_m=E_m-E_0>0$ を取る。$H_0+uB$ の基底の一階の補正は $u\sum_{m>0} \ket m\mel mB0/(E_0-E_m)$ である。固有方程式を $u$ で微分し、$\ket m$ を左から掛ければこの係数が出る。$\langle0|\partial_u0\rangle=0$ に位相を選べるので、観測量 $A$ の一次は
$$
\partial_u\langle A\rangle=-\sum_{m>0}
\frac{A_{0m}B_{m0}+B_{0m}A_{m0}}{\Delta_m}.
$$
一方、遅延関数の Lehmann 分母は $z-\Delta_m$ と $z+\Delta_m$ で、後者の分子には交換子の負号がある。従って $z=0$ では同じ式になる。$A=B$ ならこれは $-2\sum_m |B_{m0}|^2/\Delta_m$ で非正である。有限の基底縮退、熱占有、連続体がある場合には再平衡や極限の順序を指定する必要がある。

:::example[二つの励起ギャップの重みを比較する]
基底からエネルギー $2$ と $5\,\mathrm{meV}$ の二つの励起への行列要素をともに $1$ とする。静的な感受率は $-2(1/2+1/5)=-7/5\,\mathrm{meV^{-1}}$ である。遷移のスペクトルの各線の積分重みは等しいが、静的応答では低い方の寄与が大きい。逆にエネルギー重みの総和は $2+5=7\,\mathrm{meV}$ なので、高い線が大きい。静的な勾配だけで全ての分光線を復元できないのは、測るモーメントが違うためである。
:::

### 短い駆動で吸収する確率を計算する

非縮退の基底 $\ket0$ から時刻 $0<t<T$ の間だけ $u_0B$ を加え、終わった後の $H_0$ のエネルギーを測る。一次の遷移振幅は $-iu_0 B_{m0}\int_0^T e^{i\Delta_m t}\dd t$ である。従って励起の確率は二次で
$$
P_m=4u_0^2|B_{m0}|^2\,\frac{\sin^2(\Delta_m T/2)}{\Delta_m^2}.
$$
基底の生存確率はユニタリ性により $1-\sum_m P_m$ である。全エネルギーの増加は $\sum_m\Delta_m P_m\ge0$ になる。静的なポテンシャルを途中で置いても、開始と終了で外場が変化するため有限周波数の成分を持ち、励起が残り得る。$T$ が非常に短ければ各 $\sin(\Delta_mT/2)\simeq\Delta_mT/2$ となるので、低エネルギーの分母が相殺される。この場合、遷移の重みは静的な応答の逆ギャップ重みと異なる。

:::example[同じ場の強さで終端時刻を変える]
一つの励起 $\Delta$、$|B_{10}|=1$ に対し、$T=\pi/\Delta$ なら確率は $4u_0^2/\Delta^2$、$T=2\pi/\Delta$ なら二次の値は零である。開始から終了までの場の時間積分は後者の方が大きいが、遷移振幅の位相が一周して相殺される。例えば $u_0/\Delta=1/20$ なら最初の確率は $1/100$ であり、場が弱い近似を使える。より強い場では時間発展を全次数で取るため、展開の零点は変わり得る。
:::

## 松原形式による応答関数の計算

### 熱重みの差とバブル

$\chi_{AB}(\tau)=-\langle T_\tau A(\tau)B\rangle_c$ と定義する。:eqref[eq-mb-thermal-spectral] の Bose 型と同様に固有状態を入れ、$\Delta=E_s-E_r$ とすると、$\Delta\ne0$ の各項は $(p_r-p_s)A_{rs}B_{sr}/(i\nu_m-\Delta)$ になる。接続 $i\nu_m\to\omega+i0$ で :eqref[eq-mb-kubo] の周波数表示が得られる。$\Delta=0$ の連結成分は $\nu_m=0$ のみに寄与し、遅延交換子のスペクトルからは消える。

自由な並進対称 Fermi 系で $n_{\vb q}=\sum_{\vb k\sigma} c_{\vb k\sigma}^\dagger c_{\vb k+\vb q,\sigma}$ と取る。以下、密度応答は観測点に $n_{\vb q}$、源に $n_{-\vb q}$ を置き、相関を体積 $V$ で割る。Wick 分解で連結項は交換縮約一つである。相関の定義の負符号が縮約の符号を消すので、単位体積あたり
$$
\chi_0(\vb q,z)=\frac1V\sum_{\vb k\sigma}\frac{f(\xi_{\vb k})-f(\xi_{\vb k+\vb q})}{z+\xi_{\vb k}-\xi_{\vb k+\vb q}} .
\label{eq-mb-free-density}
$$
周波数和は :eqref[eq-mb-matsubara-two-poles] と Bose 外部周波数の周期性で求めた。表示式は動的な解析関数である。零松原点を求めるときは、$\xi_{\vb k}=\xi_{\vb k+\vb q}$ の項に $f'(\xi_{\vb k})$ を入れる。この値は二本の松原分母が一致する場合の和であり、非零の周波数から接続した等エネルギーの項は零である。

:::example[極限の順序]
波数が厳密に $\vb q=0$ で周波数 $z\ne0$ なら各分子は零である。一方、静的に先に $z=0$ とすると、零温度で $\lim_{\vb q\to0}\chi_0=-D(\mu)$ となる。$D$ は単位体積あたりの状態密度であり、$f'=-\delta(\xi)$ から得た。一定密度の Fermi 海の静的な圧縮性と粒子数保存が両立するのは、この順序の違いによる。
:::

:::example[有限波数の二つの準位]
$\xi_{\vb k}=0$、$\xi_{\vb k+\vb q}=\log3/\beta$ とした遷移では占有は $1/2,1/4$ である。バブルのこの項は体積因子を除いて $(1/4)/(i\nu_m-\log3/\beta)$ となる。反対の遷移を含める密度相関では、負の周波数にも極の項を足す。
:::

### 静的な接続が失う熱重み

$B$ を Hermite 演算子とし、エネルギー固有値 $E$ に射影する $P_E$ を取る。$B_0=\sum_E P_EBP_E$ は $K$ と可換であり、時間に依存しない。$B-B_0$ の行列要素にはエネルギー差があるから、$B_0$ との交差項は熱トレースで零である。従って $-\langle T_\tau\delta B(\tau)\delta B(0)\rangle$ の定数部分は $-\operatorname{Var}_\beta B_0$ であり、周波数では $-\beta\operatorname{Var}_\beta B_0\,\delta_{m0}$ になる。遅延交換子に同じ項を入れると、等エネルギーでは熱重みの差が零なので消える。静的な熱力学微分を動的接続から復元するとき、この寄与を別途計算する。

:::example[保存量にかかる準位変化]
二準位 $\pm h$ の Hamiltonian $K=h\sigma_z$ と $B=\sigma_z$ を取る。$\langle B\rangle=-\tanh(\beta h)$ だから、静的な微分は $-\beta\operatorname{sech}^2(\beta h)$ である。$B$ の交換子は常に零で動的な $\chi^R$ は零である。$\beta h=\log3$ なら $\tanh(\beta h)=4/5$ なので、零松原点にのみ $-9\beta/25$ が出る。
:::

### 自由な多バンド応答に入る行列要素

軌道基底の密度型の場を $B=\sum_{ij}B_{ij}c_i^\dagger c_j$ とし、一体 Hamiltonian の固有基底の行列を $B_{ab}$ とする。連結な虚時間相関のうち二つの頂点をまたぐ縮約だけを取ると、$0<\tau<\beta$ で $-\langle B(\tau)B(0)\rangle_c=-\sum_{ab} |B_{ab}|^2 e^{-(\xi_b-\xi_a)\tau} f_a(1-f_b)$ である。両スピンやバンドの状態を全て含め、行列要素を別に掛ける。$\xi_a\ne\xi_b$ の対には積分で
$$
\chi_{BB}(i\nu_m)=\sum_{ab}\frac{(f_a-f_b)|B_{ab}|^2}{i\nu_m+\xi_a-\xi_b}
$$
の項が現れる。等エネルギーの対では $\nu_m\ne0$ の係数は零だが、$\nu_m=0$ には $-\beta f_a(1-f_a)|B_{ab}|^2$ を加える。この部分は縮退した部分空間の基底を回転しても $\sum|B_{ab}|^2$ が同じであるため、基底の選び方に依存しない。行列要素を省略してバンドの占有差を足すだけでは、同じ分散を持つ異なる観測量を区別できない。

:::example[同じ二準位に一様場と差の場を加える]
$h=-t\sigma_x$ で一電子の基底は $\sigma_x=+1$ とする。一様場 $B=I$ はどの非対角行列要素も持たず、固定粒子数での応答は零である。差の場 $B=\sigma_z$ には二固有状態の間で行列要素 $1$ があり、ギャップは $2t$ なので静的感受率は $-1/t$ である。$t=2\,\mathrm{meV}$ なら $-1/2\,\mathrm{meV^{-1}}$ となり、一様場と差の場で応答が異なる。同じ二本のバンドでも、位置の形状が行列要素を通じて選択則を決める。
:::

### 等温応答の全体を Gibbs の微分から求める

$K(u)=K_0+uB$ を考える。熱平均の微分は
$$
\partial_u\langle B\rangle_\beta
=-\int_0^\beta\langle\delta B(\tau)\delta B(0)\rangle_\beta\dd\tau .
$$
これは指数の微分を $\tau$ の積分で書き、分母 $Z$ の微分を引けば得られる。固有基底で $B_0=\sum_E P_EBP_E$ と $B_\perp=B-B_0$ を分ける。$B_0$ の平均と分散は等エネルギーから来て、応答に $-\beta\operatorname{Var}_\beta B_0$ を与える。$B_\perp$ は $\Delta\ne0$ の遷移のみを持ち、各積分は $(p_r-p_s)/(E_s-E_r)$ の形になる。従って等温の静的微分は動的な遅延関数を零周波数に近づけた値と、この対角のゆらぎの和である。有限系で熱を交換しない測定では後者が現れない場合がある。浴の緩和をモデルに入れるなら、緩和時間と周波数のどちらを先に極限に取るかで静的値が変わる。

:::example[斜めの外場を同じ二準位に入れる]
$K_0=h\sigma_z$、$B=\alpha\sigma_z+b\sigma_x$ とする。$\sigma_z$ と $\sigma_x$ の交差した熱トレースは零である。従って等温静的応答は $-\beta\alpha^2\operatorname{sech}^2(\beta h)-b^2\tanh(\beta h)/h$ である。孤立系の遅延関数の零周波数極限は第二項だけになる。例えば $h=1\,\mathrm{meV}$、$\alpha=b=1$、$\beta h=\log3$ なら動的極限は $-4/(5\,\mathrm{meV})$、等温の差は $-9\beta/25$ である。外場の向きを傾けると準位の分布を動かす成分と準位を混ぜる成分が同時に入る。
:::

### 外部周波数の零と連続スペクトルの極限

有限次元で $\Delta=E_s-E_r$ の非零遷移を先に保持する。熱重みは $p_s=e^{-\beta\Delta}p_r$ なので静的な分子・分母の比は $(p_r-p_s)/(-\Delta)=-p_r(1-e^{-\beta\Delta})/\Delta$ であり、$\Delta\to0$ で $-\beta p_r$ に近づく。一方、非零の周波数で先に $\Delta\to0$ とすると分子が零になり、遅延の寄与は消える。従って有限次元でも、二つの準位を縮退に近づけたとき何を零にしたかに注意する。熱力学極限では $\Delta$ が任意に小さい遷移を持つため、この入れ替えは積分の一様収束を確認して行う。輸送での拡散極や Goldstone モードは、外部波数と周波数の極限を保ったまま解析する理由になる。

:::example[縮退に近い二軌道を混ぜる場]
$K=(\Delta/2)\sigma_z$、$B=\sigma_x$ とする。対角のゆらぎは $\Delta>0$ ではないため、静的な値は $-2\tanh(\beta\Delta/2)/\Delta$ である。$\Delta\to0$ の値は $-\beta$ である。逆に $\Delta=0$ を最初に置けば $B$ は保存であり、$\chi^R(z)=0$ であるが等温応答は $-\beta\langle B^2\rangle=-\beta$ のままになる。固有基底がどのように縮退するかを考えると、静的重みが遅延で零になる過程を追える。
:::

## 密度応答と動的構造因子

### 遷移確率と詳細釣り合い

$S(\vb q,\omega)=V^{-1}\sum_{rs}p_r|\mel r{n_{\vb q}}s|^2\delta(\omega-E_s+E_r)$ を :term[動的構造因子]{en="dynamic structure factor" yomi="どうてきこうぞういんし"} と呼ぶ。$p_s=e^{-\beta\omega}p_r$ と $n_{\vb q}^\dagger=n_{-\vb q}$ を使い、状態を入れ替えると $S(-\vb q,-\omega)=e^{-\beta\omega}S(\vb q,\omega)$ である。両向きを差し引いた交換子のスペクトルから
$$
-\frac1\pi\Im\chi^R_{nn}(\vb q,\omega)=(1-e^{-\beta\omega})S(\vb q,\omega).
\label{eq-mb-fdt}
$$
弾性線 $\omega=0$ は交換子に出ない。$\int S\,\dd\omega$ は等時刻の密度積を体積で割ったものである。

### 密度のエネルギー重み

反転対称な系で二重交換子を固有状態で展開すると $\langle[n_{-\vb q},[H,n_{\vb q}]]\rangle/V=2\int\omega S(\vb q,\omega)\dd\omega$ となる。連続空間で $H=\sum_j p_j^2/(2m)+U(\{\vb r\})$ とし、$n_{\vb q}=\sum_j e^{-i\vb q\cdot\vb r_j}$ とする。位置ポテンシャルは密度と可換であり、$e^{i\vb q\cdot\vb r}\vb p e^{-i\vb q\cdot\vb r}=\vb p-\vb q$ である。二重交換子に入る二つの運動エネルギー差を足すと $q^2/m$ である。したがって全周波数の一次モーメントは $nq^2/(2m)$ であり、$n=N/V$ である。SI のエネルギー変数では $q^2$ を $\hbar^2q^2$ に替える。

:::example[一つの遷移の熱重み]
$\beta\Delta=\log 3$、密度の遷移行列要素の二乗と $1/V$ の積が $2$ の二準位では、線の重みは正側 $3/2$、負側 $1/2$ である。$\omega=\Delta$ での $-\Im\chi^R/\pi$ の重みは $(1-1/3)(3/2)=1$ である。
:::

:::example[エネルギー重みを取る]
零温度において、ある波数で正の二本の線が $\omega=2,4$、重み $1/2,1/4$ を持つとする。一次は $2(1/2)+4(1/4)=2$、零次は $3/4$ なので平均励起は $8/3$ である。連続系の総和則が $nq^2/(2m)=2$ なら、この二本で一次モーメントを満たす。
:::

### スペクトル重みの移り方

二つのエネルギー準位 $\Delta>0$ の間だけを結ぶ Hermite な演算子を考え、行列要素の大きさを $b$ とする。スペクトル $S$ の規約を $S(\omega)=\sum_{rs}p_r|B_{rs}|^2\delta(\omega-E_s+E_r)$ とすれば、正の線は $p_g b^2$、負の線は $p_e b^2$ である。積分は $b^2$ で温度に依存しないが、一次モーメントは $\Delta b^2(p_g-p_e)$ で減少する。交換子スペクトルは $S(\omega)-S(-\omega)$ なので、正の遷移が増えるかどうかだけでなく逆遷移との重みの差が応答を決める。

:::example[詳細釣り合いを積分重みで測る]
$\beta\Delta=\log 3$、$b=2$ とすると $S$ の二本の積分重みは $3$ と $1$ である。一方で正周波数の $-\Im\chi^R/\pi$ の重みは $2$ であり、負周波数側では $-2$ である。静的構造因子に相当する総重み $4$ と、遅延の正周波数の重み $2$ は別の量である。
:::

### 構造因子を時間相関から作る

密度の記号を $n_{\vb q}$ にそろえる。$\langle n_{\vb q}(t)n_{-\vb q}(0)\rangle/V$ に固有状態を挿入し、$e^{i\omega t}$ を掛けて $\dd t/(2\pi)$ で積分すると、
$$
S(\vb q,\omega)=\frac1V\sum_{rs}p_r
|\mel r{n_{\vb q}}s|^2\delta(\omega-E_s+E_r)
$$
となる。$S$ を $\omega$ で積分すると $V^{-1}\langle n_{\vb q}n_{-\vb q}\rangle$ である。詳細釣り合いは添字 $r,s$ を交換すると $S(-\vb q,-\omega)=e^{-\beta\omega}S(\vb q,\omega)$ になる。逆の波数が現れるのは、随伴で $n_{\vb q}$ が $n_{-\vb q}$ へ移るためである。反転対称性などが $S(-\vb q,\omega)=S(\vb q,\omega)$ を保証する場合は、同じ波数の正負の周波数を比較できる。

遅延の交換子では二つの順序の差が入る。従って実軸境界の $-\Im\chi^R/\pi$ は $S(\vb q,\omega)-S(-\vb q,-\omega)=(1-e^{-\beta\omega})S(\vb q,\omega)$ である。正の周波数で非負であり、負の周波数では負となる。等エネルギー遷移の $\omega=0$ に独立な弾性線がある場合、係数 $1-e^{-\beta\omega}$ がこれを消すため、逆の割り算で弾性線は復元できない。この点は回折で現れる弾性散乱と、外場からのエネルギー吸収を区別するために必要である。

:::example[弾性線を加えても散逸が変わらない]
$S(\omega)=A\delta(\omega-\Delta)+Ae^{-\beta\Delta}\delta(\omega+\Delta)+C\delta(\omega)$、$A,C\ge0$ とする。交換子の正の線は $A(1-e^{-\beta\Delta})$、負の線はその負であり、$C$ の線は消える。$\beta\Delta=\log4$、$A=4$ なら正の線は $3$ であり、$C=0$ と $C=10$ の両方で同じである。静的な等時刻構造因子は $5+C$ となるので、同じ遅延の散逸スペクトルでも構造因子の全強度は違う。
:::

:::example[波数の反転で励起エネルギーが変わる環]
一粒子を三サイトの周期環に入れ、Peierls 位相 $\theta=\pi/6$ の分散 $\epsilon(k)=-2t\cos(k-\theta)$、$t>0$ を取る。格子定数を $1$ とすると、$k=0,2\pi/3,-2\pi/3$ のエネルギーは $-\sqrt3t,0,+\sqrt3t$ である。基底は $k=0$ であり、$n_q=\sum_k c_k^\dagger c_{k+q}$ に対して $|\langle 0|n_q|q\rangle|=1$ である。従って $q=2\pi/3$ では $S(q,\omega)$ は $\omega=\sqrt3t$ の線、$S(-q,\omega)$ は $\omega=2\sqrt3t$ の線を持つ。いずれも重みは一粒子の計算を体積で割った $1/V$ である。波数の符号を固定するには、最初にどの演算子を時間発展させたかを指定する。
:::

### 密度の二重交換子に相互作用が直接入らない条件

$H=\sum_j\vb p_j^2/(2m)+V(\vb r_1,\ldots,\vb r_N)$、$\rho_{\vb q}=\sum_j e^{-i\vb q\cdot\vb r_j}$ とする。$V$ は位置にだけ依存するので $\rho$ と可換である。運動項では $e^{i\vb q\cdot\vb r}\vb p e^{-i\vb q\cdot\vb r}=\vb p-\hbar\vb q$ であり、異なる粒子の因子は可換である。従って
$$
[\rho_{-\vb q},[H,\rho_{\vb q}]]
=\sum_j\frac{(\vb p_j-\hbar\vb q)^2+(\vb p_j+\hbar\vb q)^2-2\vb p_j^2}{2m}
=\frac{N\hbar^2q^2}{m}.
$$
ここでは周波数をエネルギーとして $\omega=E_s-E_r$ の単位に取る。Lehmann 表示を入れると、二重交換子は $V\int \omega[S(\vb q,\omega)+S(-\vb q,\omega)]\dd\omega$ である。反転対称なら二つの積分が等しく、単位体積あたりの :term[密度の f 総和則]{en="density f-sum rule" yomi="みつどのえふそうわそく"} を得る。
$$
\int_{-\infty}^{\infty}\omega S(\vb q,\omega)\dd\omega
=\frac{n\hbar^2q^2}{2m}.
\label{eq-mb-density-fsum}
$$
運動量依存相互作用やバンドへの射影は $\rho$ との可換性を変えるため、同じ二重交換子から計算し直す。位置だけに依存する相互作用は密度の二重交換子に寄与しないが、スペクトルでは各ピークと連続体への重みの分配を変え得る。

:::example[一つのモードに強度が移る場合]
二つの遷移のエネルギーを $E_1=1,E_2=3$、総和則の右辺を $M_1=2$ とし、零温度で正の重みのみを考える。低い方の重みが $1/2$ なら高い方は $(2-1/2)/3=1/2$ である。同じ $M_1$ を低いモードだけが運ぶ場合、そのモードの重みは $2$ になる。強度の総和は $1$ から $2$ へ変わるが、エネルギーを掛けた和は変わらない。密度の総和則で制約されるのは、この重み付きの組合せである。
:::

### 対称化したゆらぎと散逸を結び付ける

Hermite な $B$ の連結相関で、二つの演算子の順序を足して $C_{\rm sym}(t)=\langle\{\delta B(t),\delta B(0)\}\rangle/2$ とする。スペクトルを $\int C(t)e^{i\omega t}\dd t/(2\pi)$ の規格化で取ると、これは $(S(\omega)+S(-\omega))/2$ である。交換子のスペクトルが $S(\omega)-S(-\omega)$ であるため、詳細釣り合いを用いて非零周波数では
$$
\frac{S(\omega)+S(-\omega)}2
=-\frac{\coth(\beta\omega/2)}{2\pi}\Im\chi_{BB}^R(\omega)
$$
となる。$|\beta\omega|\ll1$ では $\coth(\beta\omega/2)\simeq2/(\beta\omega)$ であり、散逸の虚部が $\omega$ に比例して消える系では、対称ゆらぎは有限の定数を持つ。$\beta\to\infty$ では $\coth$ の値は $\sgn\omega$ になり、零温度でも正負両側に対称な量子ゆらぎが残る。$\omega=0$ の弾性のデルタは交換子から消えるため、この非零周波数の式では別に扱う。

:::example[ゆらぎの差と和を同じ二準位で比較する]
遷移の正の線の重みを $3$、負の線の重みを $1$ とする。正の周波数で対称ゆらぎの線は $2$、交換子の線も $2$ になる。この一致は $\beta\Delta=\log3$ のとき $\coth(\beta\Delta/2)=2$ だからである。零温度で正の重みを $3$ に保ち負の線を零にすれば、対称ゆらぎの線は $3/2$、交換子は $3$ になる。順序を足した測定と、吸収・放出の差を測る応答の関係は温度によって変わる。
:::

## 電流応答と電気伝導度の久保公式

### ベクトルポテンシャルと接触項

電荷 $e_*$ の粒子に対し最小結合を $(\vb p-e_*\vb A)^2/(2m)$ とする。一様な $\vb A$ なら $H(\vb A)=H(0)-J^p_i A_i+VD_{ij}A_iA_j/2$、$J^p_i=e_*\sum p_i/m$、$D_{ij}=ne_*^2\delta_{ij}/m$ である。電流密度 $j_i=-V^{-1}\partial H/\partial A_i$ は $J_i^p/V-D_{ij}A_j$ になる。式 :eqref[eq-mb-kubo] に $f=-A_j$ を入れ、$\chi_{ij}^R=-i\theta(t)\langle[J_i^p(t),J_j^p]\rangle/V$ と書くと $\delta j_i=-(\chi^R_{ij}+D_{ij})A_j$ である。$E=i\omega A$ を使うと
$$
\sigma_{ij}(\omega)=\frac{i}{\omega+i0}\,[D_{ij}+\chi^R_{ij}(\omega)].
\label{eq-mb-kubo-conductivity}
$$
格子で結合位相を入れると同様に $D_{ij}=V^{-1}\langle\partial_{A_i}\partial_{A_j}H\rangle$ となり、単一帯なら $(e_*^2/V)\sum_k n_k\partial_{k_i}\partial_{k_j}\epsilon_k$ となる。$\hbar$ を戻すとこの格子式に $\hbar^{-2}$ が付く。

:::example[並進不変な放物線分散]
位置の差だけに依存する相互作用では $J^p=e_*\vb P/m$ が保存されるため、同じ成分の $\chi^R$ は零である。$\sigma=iD/(\omega+i0)$ から $\Re\sigma=\pi D\delta(\omega)$ を得る。運動量緩和がない条件がこのデルタ線の由来である。
:::

:::example[接触項がギャップ内に与える制約]
絶縁体で一様な静的ベクトルポテンシャルに電流が出ない場合、$\chi^R(0)=-D$ が必要である。例えば $D=2$、$\chi^R(z)=8/(z^2-4)$ の二準位の模型なら静的な和は零である。小さい $z$ では和が $-z^2/2$ から始まり、伝導度は $-iz/2+O(z^3)$ となる。
:::

### 有限緩和をもつ Drude 応答

密度 $n$、電荷 $Q$、質量 $m$ の粒子の速度が $\partial_t v+v/\tau=QE/m$ に従うと仮定する。遅延境界で解くと $v(t)=Qm^{-1}\int_{-\infty}^t e^{-(t-t')/\tau}E(t')\dd t'$ である。$j=nQv$ を使い Fourier 変換すると
$$
\sigma(\omega)=\frac{nQ^2/m}{\tau^{-1}-i\omega}.
$$
虚部と実部を有理化すると、実部は Lorentz 型で、$\omega>0$ の虚部は正である。直流値は $nQ^2\tau/m$ になり、$\tau$ が長いほど幅 $\tau^{-1}$ は狭まる。$\tau$ は運動量の緩和時間として置いたもので、一粒子のスペクトルの寿命と区別する必要がある。

:::example[半値周波数と位相]
$\tau=0.2\,\mathrm{ps}$ なら角周波数での半値は $5\,\mathrm{ps^{-1}}$、エネルギーでは $\hbar/\tau\simeq3.29\,\mathrm{meV}$ である。この周波数では $\Re\sigma=\Im\sigma=\sigma_{\rm dc}/2$ となる。$e^{-i\omega t}$ の規約で $\sigma$ の位相が $\pi/4$ なので、実電流の最大は電場の最大より $\pi/(4\omega)$ 遅れる。
:::

### 電流演算子を外場の微分で定義する

均一なベクトルポテンシャル $\vb A$ に対し、粒子の電荷 $Q$ を含む運動エネルギーは $\sum_j(\vb p_j-Q\vb A)^2/(2m)$ である。単位体積で $\vb j=-V^{-1}\partial_{\vb A}H$ と定義すると $\vb j=\vb J_p/V-(nQ^2/m)\vb A$、$\vb J_p=(Q/m)\sum_j\vb p_j$ を得る。第一項は運動する状態の応答、第二項は演算子自体が外場に依存する部分である。一次の Hamiltonian の摂動は $-\vb A\cdot\vb J_p$ なので、$-i\theta\langle[J_{p,i}(t),J_{p,j}(0)]\rangle/V$ と定義した遅延関数を $\chi_{ij}^R$ とすれば
$$
\delta j_i=-\sum_j(\chi_{ij}^R+D\delta_{ij})A_j,
\qquad D=nQ^2/m.
$$
$\vb E=-\partial_t\vb A$ で $\omega\ne0$ には $E_j=i\omega A_j$ である。因果的な積分で $\omega$ の代わりに $z=\omega+i0$ を使うと、$\sigma_{ij}=i(\chi_{ij}^R+D\delta_{ij})/z$ を得る。ここで $\chi_{ij}^R$ は常磁性電流の相関で、接触項はまだ含まれない。この二つをまとめたものが外場に対する物理電流の応答である。

空間依存の $\vb A$ を入れるには $(\vb p-Q\vb A(\vb r))^2$ の左右を展開する。外場に一次の項は $-Q(\vb p\cdot\vb A+\vb A\cdot\vb p)/(2m)$ であり、行列要素では入射運動量と出射運動量の平均を使う。これは座標依存の場に対して Hermite 性を保つためである。縦方向の接触項と相関の相殺を調べるには、この平均の頂点と同じ外場の結合を使わなければならない。

:::example[並進対称な連続自由系の加速]
位置の差だけに依存する相互作用では全運動量 $\sum_j\vb p_j$ は保存され、均一電流の同じ方向の遅延交換子は零である。このとき $\sigma(z)=iD/z$ であり、時間の核は $D\theta(t)$ である。一定の電場 $E_0$ を $t=0$ に加えると $j(t)=DE_0t$ になる。$\omega=0$ の極は有限の直流値ではなく、運動量を失う機構がないことで現れる。例えば合計 $D=3$ の単位で $E_0=2$ を時間 $1/2$ だけ加えると、場を切る直前の電流は $3$ である。場を切ってもその後に散乱がなければこの電流は保存される。
:::

:::example[接触項の規模を見積もる]
自由電子の密度 $n$ を二倍にし質量を固定すると $D$ は二倍である。一方、電荷の符号を反転すると加速する速度の符号も反転し、電流への電荷の因子も反転するため、$D$ は同じになる。従ってこの応答から電荷の符号を直接読み取ることはできない。Hall の横応答には磁場で曲がる向きが入るので、この対角の式とは異なる符号情報を持つ。
:::

### 閉じ込められた粒子で接触項の相殺を計算する

一粒子を一次元の調和ポテンシャル $m\Omega^2x^2/2$ に置く。単位体積の因子を省き、$\hbar=1$ とする。常磁性電流は $J_p=Qp/m$ であり、$p=i\sqrt{m\Omega/2}(a^\dagger-a)$ と書ける。基底から結ぶのは一量子だけで、$|\mel1{J_p}0|^2=Q^2\Omega/(2m)=D\Omega/2$ である。従って
$$
\chi_{JJ}^R(z)=\frac{D\Omega}2
 \left(\frac1{z-\Omega}-\frac1{z+\Omega}\right)
 =\frac{D\Omega^2}{z^2-\Omega^2}
$$
である。この値に接触項 $D$ を足すと $D+\chi=D z^2/(z^2-\Omega^2)$、伝導度は $\sigma(z)=iDz/(z^2-\Omega^2)$ となる。従って $z=0$ の Drude 極は消え、正周波数の $\Omega$ に $\pi D/2$ の吸収線が現れる。静的な均一 $\vb A$ はこの非周期の閉じ込めで純ゲージとして消せるため、零周波数で相関と接触項が相殺することを別に確認できる。周期環では同じ均一 $\vb A$ が境界の磁束を変え得るので、静的な相殺を全ての境界へ仮定しない。

:::example[同じ総重みで共鳴の位置だけが動く]
$D=4$、$\Omega=2$ なら $\chi_{JJ}^R(0)=-4$ であり、$0<\omega<2$ の実吸収は零、$\omega=2$ の線の重みは $2\pi$ である。$\Omega$ を二倍にしても正の線の積分重みは $2\pi$ のままで、位置だけが $4$ へ移る。$\Omega\downarrow0$ で両側の吸収線が原点へ合流すると、全周波数での重みは $\pi D$ の Drude デルタに戻る。測定の正半軸で原点の線を半分数える規約を使えば、総和則はこの極限でも連続につながる。
:::

### 有限の減衰を入れたまま振動子の応答を解く

古典的または線形化した調和振動子で、位置の平均が $\ddot x+\gamma\dot x+\Omega^2x=QE/m$ を満たすとする。$j=nQ\dot x$ を用い、周波数を Fourier 変換すると
$$
\sigma(\omega)=\frac{-iD\omega}{\Omega^2-\omega^2-i\gamma\omega},\qquad
\Re\sigma(\omega)=\frac{D\gamma\omega^2}{(\Omega^2-\omega^2)^2+\gamma^2\omega^2}.
$$
$\gamma>0$ なら実部は非負である。$\Omega=0$ では $\omega$ の共通因子を消して $D/(\gamma-i\omega)$ に戻るが、$\Omega>0$ を固定した静的な応答は $\Re\sigma\propto\omega^2$ で零へ向かう。これは位置が有限の値へ変わり、定常電流が残らないという運動方程式と一致する。周波数を先に零にするか、閉じ込めを先に外すかで直流値が違う。

:::example[共鳴と直流を別々に測る]
$\Omega=3,\gamma=1,D=2$ なら共鳴点 $\omega=\Omega$ の実部は $D/\gamma=2$、$\omega=0$ の値は零である。一方、閉じ込めを外した $\Omega=0$ の場合は直流が $2$ であり、$\omega=3$ の実部は $2/(1+9)=1/5$ になる。散乱係数 $\gamma$ をどちらでも同じに置いても、吸収が現れる周波数範囲は変わる。
:::

## 光学伝導度と f 総和則

### 光学スペクトルの積分

対角成分の遅延応答が上半平面で解析的で十分速く減衰し、$\sigma(z)=iD/z+O(z^{-2})$ とする。上半平面の半円を反時計回りに閉じると円弧積分は $-\pi D$ になる。実軸に沿う積分の実部は従って $\pi D$ である。実時間の応答が実であり $\sigma(-\omega)=\sigma(\omega)^*$ なので
$$
\int_0^\infty \Re\sigma_{ii}(\omega)\dd\omega=\frac{\pi D_{ii}}{2}.
\label{eq-mb-optical-sum}
$$
原点のデルタ線は半分の重みで数える。有限散逸を入れても高周波の係数が同じなら総面積は変わらない。放物線分散で $n=10^{28}\,\mathrm{m^{-3}}$、電子質量なら、SI の $D=ne^2/m$ は約 $2.8\times10^{20}\,\mathrm{S\,m^{-1}\,s^{-1}}$ である。

:::example[Drude 形の面積]
$\sigma(\omega)=D/(\gamma-i\omega)$、$\gamma>0$ なら実部は $D\gamma/(\omega^2+\gamma^2)$ である。$\omega=\gamma x$ と置くと正の半軸の面積は $D\pi/2$ であり、$\gamma$ には依存しない。一方で直流値は $D/\gamma$ である。
:::

:::example[二つの緩和率を足す]
$D_1=1,D_2=3$ の独立な Drude 成分を幅 $\gamma_1,\gamma_2$ で足すと高周波の係数は $i(1+3)/z$ である。正半軸の面積は $2\pi$ である。例えば幅がそれぞれ $1,3$ なら直流値は $1+3/3=2$ であり、面積と直流値の依存性が異なる。
:::

### 格子上で総和則に入る運動エネルギー

最近接の一次元バンド $\epsilon(k)=-2t\cos(ka)$ を考え、空間一様な電磁ベクトルポテンシャルで $k\mapsto k-QA/\hbar$ とする。接触項は $D=Q^2(\hbar^2 L)^{-1}\sum_k n_k\,\partial_k^2\epsilon(k)$ である。実際、Hamiltonian の $A^2/2$ の係数を長さで割ればこの量になる。$\partial_k^2\epsilon=-a^2\epsilon$ なので、この場合 $D$ は負の運動エネルギー密度と比例する。自由なスピンなしのバンドでは $A=0$ の電流と $H$ が可換なので全光学重みが零周波数にあり、正周波数側に半分取る規約で $\int_0^\infty\Re\sigma\,\dd\omega=\pi D/2$ となる。

:::example[空帯と満帯の両方で和が零になる]
熱力学極限で占有範囲を $[-k_F,k_F]$ とすると、$L^{-1}\sum_k n_k\cos(ka)=\sin(k_Fa)/(\pi a)$ である。よって $D=2t aQ^2\sin(k_Fa)/(\pi\hbar^2)$ となる。$k_Fa=\pi/2$ では最大で、$k_Fa=\pi$ の満帯では零である。満帯で個々の電子が消えるのではなく、帯域の正と負の曲率を全て足すためである。
:::

### 因果性と高周波の係数から正半軸の重みを求める

実軸から上半平面に解析的で、境界値に実の時間応答の条件 $\sigma(-\omega)=\overline{\sigma(\omega)}$ があるとする。さらに大きな上半平面の半円上で $\sigma(z)=iD/z+O(z^{-2})$ とする。実軸にある極を上から避け、上半円の積分を評価する。$z=Re^{i\theta}$、$0<\theta<\pi$ で $\sigma(z)\dd z\simeq -D\,\dd\theta$ であり、半円は $-\pi D$ を与える。内部に極がないので実軸の積分は $\pi D$ である。実部を取ると偶関数で、$\int_0^\infty\Re\sigma(\omega)\dd\omega=\pi D/2$ になる。原点のデルタ線は実軸全体の半分を正側へ割り当てる。境界値の分布を含めるには、まず実軸より $\delta>0$ だけ上に水平な線を取り、半円積分を閉じた後に $\delta\downarrow0$ とする。

Kubo の式で $\chi_{JJ}^R(z)=O(z^{-1})$ 以下なら $\sigma$ の先頭は接触項 $iD/z$ であり、これが重みの総量を決める。応答の再配分で Drude のデルタから有限周波数の線や連続体に重みが移っても、連続の運動項が同じなら $D$ は保存される。格子では接触項がバンドの曲率の平均なので、温度や相互作用により $D$ 自身が変わり得る。従って同じ密度で測った二つの格子系の全光学重みが常に等しいとは言えない。

:::example[幅を広げて積分を照合する]
$\sigma(\omega)=D/(\gamma-i\omega)$ の実部は $D\gamma/(\omega^2+\gamma^2)$ である。$x=\omega/\gamma$ と置けば正半軸の積分は $D\int_0^\infty(1+x^2)^{-1}\dd x=\pi D/2$ である。$\gamma$ を四倍にするとピークの高さは四分の一、半幅は四倍になり、積分重みは同じである。積分上限を観測可能な $\omega_{\max}$ に止めると比は $(2/\pi)\arctan(\omega_{\max}/\gamma)$ になるため、幅の変化だけでも有限範囲の見かけの重みは変わる。
:::

### 密度の高周波展開を光学の極限と比べる

反転対称な系の密度感受率は $S(\vb q,E)$ の正と負の線の差であり、零温度の $z$ の展開の先頭は $2m_1(\vb q)/z^2$ である。実際、
$$
\frac1{z-E}-\frac1{z+E}=\frac{2E}{z^2-E^2}.
$$
モーメントが有限なら $z=iy$ で $z^2/(z^2-E^2)\to1$ を $E S$ の積分の内側で取れる。連続系の $m_1=n q^2/(2m)$ に対し $\chi_{nn}\sim n q^2/(m z^2)$ である。連続の式で縦伝導度を $\sigma_L=iQ^2z\chi_{nn}/q^2$ と結び付けると、先頭は $i nQ^2/(m z)$ になる。従って密度の重みと光学の重みは同じ運動項の交換子を測っている。縦の $q\ne0$ を保った式で $\omega$ の極限を先に取り、その後に $q\to0$ とすると静電遮蔽を含む別の応答になるため、均一光学の測定と取り違えない。

:::example[密度の線が一つだけの場合]
$\chi_{nn}(q,z)=nq^2/[m(z^2-\omega_q^2)]$ とする。これは低エネルギーの密度線が全ての重みを担う模型である。$\sigma_L=iDz/(z^2-\omega_q^2)$ なので、正の吸収線の重みは $\pi D/2$ である。音波の $\omega_q=cq$ なら $q\downarrow0$ で線は原点へ、三次元 Coulomb の $\omega_q\to\omega_p>0$ なら有限の周波数へ残る。両方で高周波の先頭は同じである。ここで密度の応答は実際に加えた外場を基準とし、誘起ポテンシャルを含む全場に対する応答とは区別する。
:::

### 格子の二重交換子を自由バンドで確かめる

スピンを省き、$H_0=\sum_k\epsilon_k c_k^\dagger c_k$ とし、$\rho_q=\sum_k c_k^\dagger c_{k+q}$ とする。$[H_0,\rho_q]=\sum_k(\epsilon_k-\epsilon_{k+q})c_k^\dagger c_{k+q}$ である。さらに $\rho_{-q}$ と交換すると $[c_{k+q}^\dagger c_k,c_k^\dagger c_{k+q}]=n_{k+q}-n_k$ なので
$$
\langle[\rho_{-q},[H_0,\rho_q]]\rangle
=\sum_k \langle n_k\rangle(\epsilon_{k+q}+\epsilon_{k-q}-2\epsilon_k)
$$
となる。第二項に $k$ の平行移動を使った。最近接の $\epsilon_k=-2t\cos ka$ では括弧は $4t(1-\cos qa)\cos ka$ であり、$q^2$ の係数は $2ta^2\cos ka=\partial_k^2\epsilon_k$ である。格子上の位置の密度と可換な Hubbard 相互作用なら、この二重交換子に直接は入らない。相互作用は $\langle n_k\rangle$ を変えるため、重みの数値にも影響する。

:::example[帯域が埋まったときの有限波数の打ち消し]
最近接鎖でスピンなしの全波数を占有すると $\sum_k\cos ka=0$ であり、二重交換子は任意の非零 $q$ でも零である。実際に $\rho_q$ は占有済みの軌道へ粒子を移すため、完全に埋まった状態を消す。半充填で $\cos ka>0$ の範囲を占有する場合は同じ波数で二重交換子が正になる。格子の総和則に裸質量を入れず、バンドと占有を合わせて積分する必要がある。
:::

## ゲージ不変性と Ward 恒等式

### 連続の式を時間順序積に入れる

$\psi$ の粒子数電流 $j^\mu=(n,\vb j)$ が $\partial_t n+\divg\vb j=0$ を満たすとする。三点関数 $\langle T j^\mu(x)\psi(1)\psi^\dagger(2)\rangle$ に発散を取ると、演算子の微分は連続の式で消える。時間順序の階段関数の微分から、等時刻交換子 $[n(\vb x),\psi(\vb r)]=-\delta(\vb x-\vb r)\psi(\vb r)$、生成に対する逆符号の項が出る。これらは :eqref[eq-mb-ccr-car] を密度に適用すれば得られる。二つの接触項を Fourier 変換すると二点関数の差になる。三点関数の両端を $G$ で切断し、自由頂点を $\Gamma^0=1$ とする規約では
$$
q_0\Gamma^0(k+q,k)-\vb q\cdot\vb\Gamma(k+q,k)
=G^{-1}(k+q)-G^{-1}(k).
\label{eq-mb-ward}
$$
これを :term[Ward 恒等式]{en="Ward identity" yomi="わーどこうとうしき"} と呼ぶ。ゲージ変換に対する応答核では電流に接触項を含める。その項を含む核に四元運動量を縮約すると、純ゲージの外場が物理量を変えないことに対応して零になる。頂点を相互作用と整合させる方法は :eqref[eq-mb-dyson] の外場微分である。

:::example[自由粒子の頂点]
自由な放物線分散では $\Gamma^0=1$、$\vb\Gamma=(\vb k+\vb q/2)/m$ であり、内積は $[(\vb k+\vb q)^2-k^2]/(2m)$ に等しい。相互作用で $\Sigma$ が周波数依存なら、$\vb q=0$ の $q_0\to0$ で $\Gamma^0=1-\partial_\omega\Sigma$ が要る。$\partial_\omega\Sigma=-1/2$ なら $\Gamma^0=3/2$ であり、極の留数 $Z=2/3$ と掛けて $1$ となる。
:::

### 縦電流と密度の応答を照合する

連続系の粒子数電流を $j^\parallel_{\vb q}$ とし、$\dot n_{\vb q}=-iqj^\parallel_{\vb q}$ とする。$E_s-E_r=\Delta$ の遷移に対して $-i\Delta\mel r{n_{\vb q}}s=-iq\mel r{j^\parallel_{\vb q}}s$ であり、電流の行列要素の絶対値の二乗は密度のそれの $\Delta^2/q^2$ 倍である。Lehmann 表示に入れると $z^2-\Delta^2=(z-\Delta)(z+\Delta)$ なので、$z^2\chi_{nn}-q^2\chi_{j^\parallel j^\parallel}$ では分母が消える。残る $z$ の項は $[n_{\vb q},n_{-\vb q}]=0$ なので消え、$\Delta$ の項は二重交換子の期待値 $nq^2/m$ になる。従って
$$
z^2\chi_{nn}(\vb q,z)=q^2\left[\chi_{j^\parallel j^\parallel}(\vb q,z)+\frac nm\right].
$$
ここで双方の応答は単位体積あたりであり、電流は電荷を掛けていない。式を電荷電流へ直すと $n/m$ の項が :eqref[eq-mb-kubo-conductivity] の接触項になる。$\vb q\ne0$ で遅延関数を静的に評価し、$\chi_{nn}$ が有限なら縦の和は零になる。

:::example[二つの遷移に対する検算]
密度の交換子スペクトルが $\omega=\pm\Delta$ に重み $\pm w$ を持つとする。密度応答は $2w\Delta/(z^2-\Delta^2)$、電流応答はその $\Delta^2/q^2$ 倍である。従って差は $2w\Delta$ である。例えば $w=1/4,\Delta=2$ なら値は $1$ となり、連続系のモデルと整合させるには $nq^2/m=1$ を満たす必要がある。
:::

### 格子の有限差分での頂点関係

一バンドの逆伝播 $G_0^{-1}(k,z)=z-\epsilon(k)$ を使う。有限差の密度頂点を $\Gamma^0=1$ と決めれば、Ward 恒等式の縦成分は $q\Gamma^L=\epsilon(k+q)-\epsilon(k)$ である。バンド $\epsilon=-2t\cos ka$ なら差は $4t\sin(qa/2)\sin(ka+qa/2)$ である。結合上の電流の格子発散には $2\sin(qa/2)/a$ が掛かるので、格子頂点を $2ta\sin(ka+qa/2)$ とすれば同じ差を再現する。$q\to0$ では $\partial_k\epsilon$ に戻る。有限 $q$ で単に $q\,\partial_k\epsilon$ とすると曲率の項を落とす。

:::example[帯域底で速度が零でも有限差は残る]
$ka=0$、$qa=\pi/2$ のとき差は $2t$ である。格子発散係数は $\sqrt2/a$、結合頂点は $\sqrt2ta$ なので積は $2t$ になる。一方 $\partial_k\epsilon(0)=0$ なので、点での速度だけを使うと連続の式を満たさない。有限波数の密度計算では、外場を結合に入れてこの位相を保持する。
:::

### 二時刻の Ward 恒等式の接触項を確認する

密度 $\rho(\vb r,t)$ と粒子流 $\vb j$ が演算子として $\partial_t\rho+\nabla\cdot\vb j=0$ を満たすとする。時間順序積に $\rho(1)$ と Fermi 場 $\psi(2),\psi^\dagger(3)$ を挿入する。$t_1$ の微分を分布として取ると、$\rho$ の微分、さらに $\rho$ が各場と同時刻になるときの左と右の差が出る。密度は Fermi 偶なので、後者の順序の差は $[\rho(\vb r_1,t_1),\psi(\vb r_2,t_1)]=-\delta(\vb r_1-\vb r_2)\psi(\vb r_2,t_1)$ と $[\rho,\psi^\dagger]=+\delta\psi^\dagger$ である。従って
$$
\begin{aligned}
\partial_{t_1}\langle T\,\rho(1)\psi(2)\psi^\dagger(3)\rangle
+\divg_1\langle T\,\vb j(1)\psi(2)\psi^\dagger(3)\rangle
={}&[-\delta(1-2)+\delta(1-3)]
 \langle T\,\psi(2)\psi^\dagger(3)\rangle .
\end{aligned}
$$
ここで $\delta(1-2)$ は時空間のデルタである。密度の期待値と二点の積を引いた連結部分に移しても、定常な密度と連続の式のため右辺は変わらない。Fourier 変換では二つのデルタが外脚の運動量を異なる向きへずらし、三点関数の両側に逆伝播を掛けることで $G^{-1}(K+Q)-G^{-1}(K)$ が出る。頂点を切断するときのこの二つの項は、各外場挿入の電荷の符号から生じる。時間微分を順序の内側だけで取れば消えてしまい、保存則の予言と合わない。

:::example[自由粒子の全頂点をその場で検査する]
$G_0^{-1}(k,z)=z-k^2/(2m)+\mu$ で $\Gamma^0=1$、$\Gamma^L=(k+q/2)/m$ とする。差は $\nu-(2kq+q^2)/(2m)$ であり、$\nu\Gamma^0-q\Gamma^L$ と一致する。例えば $k=1,q=2,m=2,\nu=5$ なら $\Gamma^L=1$、差は $3$ である。電流を点 $k$ の速度 $k/m$ だけで置くと左辺は $4$ となる。有限の転送では入射と出射の中点の速度を使うことが必要である。
:::

### 静的な縦場が純ゲージであることを使う

一次元の無限または開いた区間で、時間に依存しない $A(x)=\partial_x\lambda(x)$ を考え、粒子の境界条件を変えない滑らかな $\lambda$ を取る。波動関数に $e^{iQ\lambda(x)/\hbar}$ を掛けると、$(p-Q\partial_x\lambda)e^{iQ\lambda/\hbar}\psi=e^{iQ\lambda/\hbar}p\psi$ である。従ってこの $\vb A$ をもつ Hamiltonian と $\vb A=0$ はユニタリ同値であり、物理電流の平均は同じである。線形応答では静的な縦の常磁性相関が接触項を打ち消す。この議論は静的な縦場であり、横場の磁気応答や、円周で消せない全磁束の感度を消すものではない。時間依存のゲージ変換を使うと作用の時間微分からスカラーポテンシャルも変わり、$-\partial_t\vb A-\grad\varphi$ の電場が同じであることが必要になる。

:::example[円周上の一定のベクトルポテンシャル]
長さ $L$ の一粒子の輪で $A=\Phi/L$ とする。$\lambda=Ax$ で場を消そうとすると波動関数の位相は一周で $e^{iQ\Phi/\hbar}$ だけ変わる。従って $\Phi$ が $h/Q$ の整数倍でなければ、元の周期境界のままでは純ゲージとして消すことができない。固有エネルギーは $(2\pi\hbar j/L-Q\Phi/L)^2/(2m)$ であり、$\Phi$ を $h/Q$ だけ変えると $j$ の名前を一つずらして同じ集合に戻る。有限輪の静的な流れを連続無限系の縦の相殺と区別する理由になる。
:::

### 保存則を近似の検査に使う

自己エネルギーが周波数に依存し $G^{-1}(z,k)=z-\xi_k-\Sigma(z,k)$ とする。転送の空間成分を零に取れば
$$
\Gamma^0(z+\nu,z;k)=1-\frac{\Sigma(z+\nu,k)-\Sigma(z,k)}{\nu}.
$$
$\nu\ne0$ での式であり、滑らかな場合は $\nu\to0$ で $1-\partial_z\Sigma$ となる。従ってドレスした伝播を使いながら裸の密度頂点 $1$ のみを保つと、$\Sigma$ の変化がある場合には Ward 恒等式を満たさない。系の密度が保存されるなら、一様な動的な密度応答は零でなければならない。これは Lehmann のレベルでは $\rho_{\vb q=0}=N$ に非対角の遷移がないことから従うので、計算結果の直接の検査になる。静的な粒子数の揺らぎは動的な一様応答とは異なり、再平衡なら有限でも矛盾しない。

:::example[密度を運ぶ規格化の確認]
$\Sigma(z)=\alpha-z$ とする狭い極の局所モデルでは $Z=1/(1-\partial_z\Sigma)=1/2$ である。動的な頂点は $\Gamma^0=2$ なので、規格化した準粒子が運ぶ粒子数の因子 $Z\Gamma^0$ は $1$ である。裸の頂点をそのまま使うと $1/2$ になる。この $\Sigma$ は低エネルギーでの線形化であり、全周波数で $G(z)\sim 1/z$ を満たす物理的な自己エネルギーを仮定したわけではない。
:::

### 拡散での連続の式と Einstein 関係

粒子数を保存し、長波長で密度の勾配による拡散が支配的な系を取る。密度の局所平衡感受率を $\kappa_n=\partial_\mu n>0$、拡散定数を $D_{\rm diff}>0$ とする。外部の粒子一個あたりのエネルギー $u$ に対し、局所平衡の密度変化は $-\kappa_n u$ である。従って流れを
$$
\vb j_n=-D_{\rm diff}\grad(\delta n+\kappa_n u)
$$
とすれば平衡で流れは消える。連続の式 $\partial_t\delta n+\divg\vb j_n=0$ を Fourier 変換すると
$$
\chi_{nn}(q,z)=-\frac{\kappa_n D_{\rm diff}q^2}{D_{\rm diff}q^2-iz}.
$$
静的な $\omega=0$ を $q\ne0$ で先に取れば $-\kappa_n$、逆順の動的な一様極限は零である。電荷 $Q$ で $\sigma_L=iQ^2z\chi_{nn}/q^2$ とすれば、$q\to0$ を先に取った低周波数の値は $Q^2\kappa_n D_{\rm diff}$ である。これが Einstein 関係であり、拡散の速さと静的な加えやすさの積が電流を決める。導出に使った流れの勾配展開は散乱の時間と距離より長い範囲に限るため、この式を高周波の光学総和則の上限まで外挿しない。

:::example[拡散の極と静的応答を比較する]
$D_{\rm diff}=0.01\,\mathrm{m^2/s}$、$q=10^6\,\mathrm{m^{-1}}$ なら緩和率は $10^{10}\,\mathrm{s^{-1}}$ である。角周波数がこの値と同じなら $\chi_{nn}=-\kappa_n(1+i)/2$ であり、静的値の実部の半分である。$q$ を半分にすると緩和率は四分の一になる。長い波長ほど密度を均す距離が長いため、同じ小さな外場でも再平衡のための時間が増える。
:::

## 演習問題

::::exercise{level=basic}
$\Delta=4\,\mathrm{eV}$ の二準位の $\sigma_x$ の静的応答を求めよ。

:::solution
二本の遷移項を $z=0$ で引くと $-2/\Delta=-1/2\,\mathrm{eV^{-1}}$ である。
:::
::::

::::exercise{level=basic}
正の遷移の構造因子の重みが $9$、$\beta\omega=\log 3$ なら逆向き遷移の重みを求めよ。

:::solution
詳細釣り合いの係数は $e^{-\beta\omega}=1/3$ なので $3$ である。
:::
::::

::::exercise{level=basic}
$D=6,\gamma=2$ の Drude 伝導度の直流値を求めよ。

:::solution
$\sigma(0)=D/\gamma=3$ である。
:::
::::

::::exercise{level=basic}
光学総和則で $D=4$ の正半軸の面積を求めよ。

:::solution
係数 $\pi/2$ を掛けるので $2\pi$ である。原点の線がある場合はその半分を積分に含める。
:::
::::

::::exercise{level=advanced}
自由分散 $\xi_k=k^2/(2m)-\mu$ で、$k=2,q=1,m=2$ の電流頂点とエネルギー差を求めよ。

:::solution
対称な頂点は $(k+q/2)/m=5/4$ である。分散差は $(3^2-2^2)/(2m)=5/4$ であり、$q$ 倍した頂点に一致する。
:::
::::

::::exercise{level=advanced}
単一粒子で連続空間の密度の二重交換子を計算せよ。

:::solution
$[H,e^{-iqx}]=e^{-iqx}((p-q)^2-p^2)/(2m)$ である。さらに $e^{iqx}$ と交換すると $(p-q)^2+(p+q)^2-2p^2$ を $2m$ で割ったものが残り、$q^2/m$ になる。位置だけのポテンシャルの交換子は零である。
:::
::::

::::exercise{level=advanced}
連続系の $D=3$、$\gamma=1$ の Drude 形で、$0\le\omega\le1$ の面積の全正半軸面積に対する比を求めよ。

:::solution
積分は $3[\arctan\omega]_0^1=3\pi/4$ であり、全正半軸では $3\pi/2$ である。従って比は $1/2$ である。
:::
::::

::::exercise{level=advanced}
自己エネルギーが $\Sigma(z)=a+bz$ なら $\vb q=0$ の密度頂点と、狭い極での $Z\Gamma^0$ を求めよ。$1-b>0$ とする。

:::solution
逆伝播関数の差は $(1-b)q_0$ であるため、Ward 恒等式から $\Gamma^0=1-b$ となる。伝播関数の極の留数は $1/(1-b)$ であるから積は $1$ である。
:::
::::

::::exercise{level=advanced}
有限の離散系で Hermite な観測量 $B$ とする。零温度の非縮退基底に対し、$\chi_{BB}^R(0)\le0$ を示し、等号条件を述べよ。外場は $+uB$ とする。

:::solution
Lehmann 表示に挿入された励起状態 $\ket m$ にはギャップ $E_m-E_0>0$ がある。正の極の分子は $|B_{m0}|^2$、負の極の分子はその負である。従って $\chi^R(0)=-2\sum_{m>0}|B_{m0}|^2/(E_m-E_0)\le0$ である。各分母が正なので、等号には全ての $B_{m0}=0$ が必要十分である。完全性によりこれは $B\ket0=B_{00}\ket0$ と同値である。$B$ が Hamiltonian と可換ならこの条件を満たすが、基底を固有状態とするだけで十分である。
:::
::::

::::exercise{level=advanced}
調和振動子の電流応答で $D>0$、$\Omega>0$ とする。式 $\sigma(z)=iDz/(z^2-\Omega^2)$ を部分分数に分け、実軸境界の実部と虚部を求めよ。

:::solution
$z/(z^2-\Omega^2)=[(z-\Omega)^{-1}+(z+\Omega)^{-1}]/2$ である。$1/(x+i0)=\operatorname{PV}(1/x)-i\pi\delta(x)$ を使うと
$$
\begin{aligned}
\Re\sigma(\omega)&=\frac{\pi D}2[\delta(\omega-\Omega)+\delta(\omega+\Omega)],\\
\Im\sigma(\omega)&=\frac D2\left[\operatorname{PV}\frac1{\omega-\Omega}
+\operatorname{PV}\frac1{\omega+\Omega}\right].
\end{aligned}
$$
実部は偶で正半軸の積分は $\pi D/2$、虚部は奇である。$\Omega\downarrow0$ で二つのデルタは $\pi D\delta(\omega)$ へ合流し、正半軸の規約ではその半分を取る。閉じ込められた場合に零周波数のデルタがないことも、この分解から分かる。
:::
::::

::::exercise{level=advanced}
長波長の拡散模型で $\chi_{nn}(q,z)=-\kappa_n Dq^2/(Dq^2-iz)$ とする。$q\ne0$ での静的な密度応答と、一様な縦伝導度の低周波数値を求め、極限の順序を示せ。

:::solution
$z=0$ を先に取れば $\chi_{nn}=-\kappa_n$ である。非零の $z$ から $q\to0$ を先に取れば、密度応答の極限は零になる。一方、縦伝導度は連続の式から $\sigma_L=iQ^2z\chi_{nn}/q^2=-iQ^2z\kappa_n D/(Dq^2-iz)$ である。従って $q\to0$ を先に取れば $\sigma_L=Q^2\kappa_n D$ となり、その後の $z\to0$ でも同じである。$q\ne0$ で $z=0$ を先に取った伝導度は零であり、静的な密度の勾配によって外場が釣り合った状態を表す。
:::
::::

## 確認テスト

::::quiz
密度の虚時間応答の外部周波数はどれか。

- [x] Bose の松原周波数
- [ ] Fermi の松原周波数
- [ ] 零周波数に限る

:::explanation
密度は偶数個の Fermi 場からなるので、時間境界が周期的になる。従って整数倍の $2\pi/\beta$ が許される。
:::
::::

::::quiz{answer="3" tol="1e-6"}
$D=6$、$\gamma=2$ の Drude 形の直流値を答えよ。

:::explanation
直流では $\omega=0$ なので分母は $\gamma=2$ となり、値は $3$ である。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
$\beta\omega=\log2$ の遷移について、負周波数側と正周波数側の重みの比を求めよ。反転対称とする。

:::explanation
反転対称なので波数の符号によって正側の重みは変わらない。詳細釣り合いで比は $e^{-\log2}=1/2$ となる。
:::
::::

::::quiz
電流応答について正しいものを選べ。

- [x] ベクトルポテンシャルに対する電流の明示的な微分が入る。
- [x] 原点のデルタ線は正半軸の総和則で半分を数える。
- [ ] 自己エネルギーを入れても密度頂点は常に $1$ である。

:::explanation
最小結合の二次項を微分すると電流の接触項が出る。原点の線は正負の周波数の対称性で等分する。周波数依存の自己エネルギーは式 :eqref[eq-mb-ward] により密度頂点を変える。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
Drude 伝導度 $\sigma(\omega)=D/(\gamma-i\omega)$、$\gamma>0$ を考える。正半軸の光学重みのうち $0\le\omega\le\gamma$ に含まれる割合を答えよ。

:::explanation
実部は $D\gamma/(\gamma^2+\omega^2)$ である。$x=\omega/\gamma$ とすると分子の積分は $D\int_0^1\dd x/(1+x^2)=D\pi/4$ であり、全正半軸は $D\int_0^\infty\dd x/(1+x^2)=D\pi/2$ である。比は $1/2$ となる。幅の位置で積分を止めたため、この割合自体は $\gamma$ によらない。
:::
::::

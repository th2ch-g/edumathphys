---
title: "対称性の自発的破れ"
description: "縮退した真空での対称性と質量を求め、大域対称性の破れとゲージ理論の摂動展開を結ぶ。"
notation:
  - symbol: 'v'
    meaning: 選んだ背景の大きさ
  - symbol: 'h,\chi'
    meaning: 動径方向と角方向の場
  - symbol: 'm_A'
    meaning: ゲージ場の質量
---

## この章の目標

- 真空の選択に必要な極限の順序を説明する。
- 線形 σ 模型と Goldstone の定理から質量ゼロの励起を得る。
- Higgs 相の二次作用を対角化し、高エネルギーでの等価定理を使う。

## 自発的対称性の破れと真空の縮退

### 対称な作用で向きを選ぶ

作用が対称でも、エネルギーを最小にする配置が複数あり、対称変換がその一つを別の配置へ写すことがある。無限体積で選んだ真空がそのような変換で不変でないことを :term[自発的対称性の破れ]{en="spontaneous symmetry breaking" yomi="じはつてきたいしょうせいのやぶれ"} という。例えば実スカラーの $V=\lambda(\phi^2-v^2)^2/4$ は $\phi\mapsto-\phi$ に不変だが、$\lambda,v>0$ では $\phi=\pm v$ に極小がある。有限体積の正則化で一意な基底状態を持つなら、対称性によりその状態の奇な期待値は零である。二重井戸の半古典近似では、この基底状態はトンネル効果による対称な組合せである。真空を一つ選ぶには源 $+h\int\phi$ を重みに加え、先に体積を無限にしてから $h\to0^+$ とする。二つの一様な配置の相対重みは $e^{2hv\mathcal V_4}$ なので、この二配置近似で $+v$ が選ばれる。先に $h=0$ なら期待値は対称性により零である。

連続対称性 $G$ では、選んだ真空を固定する部分群 $H$ と、真空を動かす方向を分ける。古典的な極小が一軌道なら、その集合は $G/H$ である。離散対称性の例では二極小をつなぐ連続な零エネルギー方向はない。

### 極限の順序を式で比較する

二相の重みだけを残した模型で $\langle\phi\rangle_{h,\mathcal V}=v\tanh(hv\mathcal V)$。固定 $h>0$ で $\mathcal V\to\infty$ を先に取れば $+v$、固定 $\mathcal V$ で $h\to0$ を先に取れば零になる。微分 $\partial_h\langle\phi\rangle_{h,\mathcal V}|_{h=0}=v^2\mathcal V$ は、極限が $h=0$ で一様でないことを示す。例えば $hv\mathcal V=\log3/2$ の曲線に沿う極限では $\tanh(\log3/2)=1/2$ なので $v/2$ になる。純相を定義するには体積と源の順序を指定する必要がある。

### 安定化部分群と真空の接空間

対称性を持つ理論を展開する前に、極小集合が単一の群軌道かを調べる。同じエネルギーを持つ二つの点が対称変換で結ばれる保証はない。基準点 $\Phi_0$ を選び、$H=\{g\in G\mid g\Phi_0=\Phi_0\}$ とする。$g_1\Phi_0=g_2\Phi_0$ なら $g_2^{-1}g_1\in H$ であり、逆も成り立つ。従って $gH\mapsto g\Phi_0$ は軌道を $G/H$ と同定する。滑らかな有限次元の作用では、その微分は $X\mapsto X\Phi_0$。核は $H$ の Lie 代数であるため、接空間の次元は $\dim G-\dim H$ となる。

:::example[対称性で結ばれない平坦方向]
二つの実ベクトル $\Phi,\Psi$ に $V=\lambda(\Phi^2+\Psi^2-v^2)^2/4$ を取る。運動項も等しく取ると実際の対称性は全成分を回す直交群まで拡大し、一つの球面に見える。一方、$\Phi,\Psi$ の別々の回転を保つ高次微分相互作用を加えれば、一様な $V$ は変化しなくても全成分を混ぜる回転は作用の対称性でなくなる。例えば $(\partial_\mu\Phi\cdot\partial^\mu\Phi)^2$ は $\Phi=0$ を含む一様な極小にエネルギーを加えない。しかし $\Phi^2$ の異なる二極小は、残した対称群の各軌道に分かれる。この場合、対称性から保証される平坦方向と、選んだポテンシャルで偶然平坦な方向を区別する。
:::

真空を選ぶときは、時間の長さを増して状態を射影する操作と、空間の体積を増す操作を区別する。空間体積を固定して Euclid 時間だけを無限にすると、源を零とした理論の基底状態へ射影するだけである。例えば局在した二状態を $\ket+,\ket-$ とし、その基底で
$$
H=E_0 1-\kappa\sigma^1-hvL^3\sigma^3,\qquad \kappa>0
$$
と近似する。$h=0$ では固有値は $E_0\mp\kappa$ で、基底状態の $\sigma^3$ の期待値は零。$h>0$ では $u=hvL^3$ とおくと、基底状態の射影は $(1+(\kappa\sigma^1+u\sigma^3)/\sqrt{\kappa^2+u^2})/2$。よって場の空間平均は $v u/\sqrt{\kappa^2+u^2}$ となる。固定 $L$ の $h\to0$ で零だが、$\kappa(L)/L^3\to0$ と仮定して $L\to\infty$ を先に取れば $+v$。二配置の有限 Euclid 体積の式は、時間発展でのトンネル混合を無視した別の近似である。

### クラスター性で純相を見分ける

無限体積で $\langle\phi\rangle_+=v_{\rm R}$、$\langle\phi\rangle_-=-v_{\rm R}$ となる二相を仮定し、それぞれが離れた挿入を因数分解するとする。対称な混合状態 $\omega=\tfrac12(\omega_++\omega_-)$ では $\omega(\phi)=0$ だが、大きく空間的に離した二点では
$$
\lim_{|\vb x|\to\infty}\omega(\phi(\vb x)\phi(0))
=\tfrac12(v_{\rm R}^2+(-v_{\rm R})^2)=v_{\rm R}^2 .
$$
従ってこの混合はクラスター性を持たない。ここで $v_{\rm R}$ は選んだくりこみ済み演算子の期待値であり、古典のポテンシャルの $v$ に等しいと仮定する必要はない。例えば $O=\phi^2$ だけを測ると両相の期待値は等しく、同じクラスター性の破れは見えない。相を区別する演算子を :term[秩序変数]{en="order parameter" yomi="ちつじょへんすう"} と呼ぶ。

### 有限体積で連続な向きを量子化する

線形 $O(N)$ 模型で動径を固定し、空間の体積 $L^3$ に一様な単位ベクトル $n(t)$ だけを残す。運動項を積分すると $I\dot n^2/2$、$I=v^2L^3$。従って自由回転子の Hamiltonian は球面 $S^{N-1}$ の $-\Delta_{S^{N-1}}/(2I)$ である。例えば $N=2$ なら $n=(\cos\theta,\sin\theta)$、波動関数の周期条件は $2\pi$。固有関数 $e^{i\ell\theta}$ のエネルギーは $\ell^2/(2I)$、基底状態 $\ell=0$ は向きに一様なので $\langle n\rangle=0$。第一励起のギャップ $1/(2v^2L^3)$ は $L\to\infty$ で消える。スピン波などの非一様なモードとは別に、真空の向きを選ぶモードが有限体積で量子化されることが分かる。

源 $-hvL^3 n_N$ を加えると、源の一つの行列要素は $\mathcal O(hvL^3)$、無源の回転子のギャップは $\mathcal O(1/(v^2L^3))$。従って有限 $L$ の小さい $h$ の線形応答を固定 $h$ の無限体積へ一様には延長できない。例えば $N=2$ で $n_N=\cos\theta$ を選ぶと、$\ell=\pm1$ との行列要素は各 $1/2$、一次摂動の期待値は $2hvL^3(2\times \frac14)(2I)=2hvL^3 I$ である。この値が一よりはるかに小さい範囲でのみ線形近似が使える。$L$ を固定せずにこの応答をそのまま積分すると、$|\langle n_N\rangle|\le1$ の条件にも抵触する。

## 線形 σ 模型

$N$ 成分の実場 $\Phi=(\pi_1,\ldots,\pi_{N-1},v+\sigma)$ を取り、$\mathcal L=(\partial\Phi)^2/2-\lambda(\Phi^2-v^2)^2/4$。$O(N)$ のうちこの真空を固定するのは $O(N-1)$。$\Phi^2-v^2=2v\sigma+\sigma^2+\pi^2$ だから
$$
V=\lambda v^2\sigma^2+\lambda v\sigma(\sigma^2+\pi^2)+\tfrac{\lambda}{4}(\sigma^2+\pi^2)^2.
\label{eq-qft-sigma-expansion}
$$
質量は $m_\sigma^2=2\lambda v^2$、$m_\pi^2=0$。$v=2\,\mathrm{GeV}$、$\lambda=1/2$ なら $m_\sigma=2\,\mathrm{GeV}$ である。例えば第 $i$ 方向と最後の成分の回転で $j_i^\mu=(v+\sigma)\partial^\mu\pi_i-\pi_i\partial^\mu\sigma$ となる。先頭項 $v\partial^\mu\pi_i$ は真空から一つの $\pi_i$ を作る。$m_\sigma$ より十分低いエネルギーで動径場を古典的に消すと、$\Phi^2=v^2$ 上の運動項が角度の先頭作用になり、相互作用に微分が付く。

### 小さな明示的破れで角方向に質量を作る

線形模型に $V_h=V-h\Phi_N$ を加え、$\Phi_N=f>0$ の極小を選ぶ。停留条件は $\lambda(f^2-v^2)f=h$。横方向の Hessian は $\lambda(f^2-v^2)$ なので $m_\pi^2=h/f$ になる。元のカレントの発散も $h$ を含み、$j_i^\mu=\Phi_N\partial^\mu\pi_i-\pi_i\partial^\mu\Phi_N$ に運動方程式を入れると $\partial_\mu j_i^\mu=-h\pi_i$ が残る。線形に $j_i^\mu\simeq f\partial^\mu\pi_i$ なので、これは $(\Box+m_\pi^2)\pi_i=0$ と同じである。例えば $h/f=0.04\,\mathrm{GeV}^2$ では角方向の小振動の質量は $0.2\,\mathrm{GeV}$。$h\downarrow0$ を取ると Goldstone 極に連続につながる。

### Goldstone 散乱で非微分項が消える

一種類の $\pi$ の $\pi\pi\to\pi\pi$ を木で計算する。:eqref[eq-qft-sigma-expansion] から四点は $-6\lambda$、各 $\sigma$ 交換は $-4\lambda^2v^2/(s-m_\sigma^2)$ など。低運動量では各交換の定数が $2\lambda$ なので三チャネルを足して四点を消す。次の項は $2\lambda(s+t+u)/m_\sigma^2$、質量ゼロの殻上ではこれも零。従って同じ一種類の角モードに残る先頭は $2\lambda(s^2+t^2+u^2)/m_\sigma^4$ である。例えば $s\ll m_\sigma^2$ と $t=u=-s/2$ なら $3\lambda s^2/m_\sigma^4$。複数の種類では色添字の δ によってチャネルの組合せが異なり、一般には二微分の振幅が残る。

### 非線形表示と動径場の消去

$\Phi=(v+r)n$、$n^2=1$ と置く。$n\cdot\partial_\mu n=0$ なので交差項がなく、
$$
\mathcal L=\frac12(\partial r)^2+
\frac12(v+r)^2(\partial n)^2-\frac{\lambda}{4}(2vr+r^2)^2 .
$$
$r$ の運動方程式を小運動量で解くと $(\Box+m_\sigma^2)r=v(\partial n)^2+\mathcal O(r(\partial n)^2,r^2)$。従って $r$ の最初の応答は $v(\partial n)^2/m_\sigma^2$ である。これを作用で代入すると、次数四までの角モードの木の作用は
$$
\mathcal L_{\rm low}=\frac{v^2}{2}(\partial n)^2+
\frac{v^2}{2m_\sigma^2}\bigl[(\partial n)^2\bigr]^2+\mathcal O(\partial^6).
$$
第一項は球面の誘導計量を用いた :term[非線形 σ 模型]{en="nonlinear sigma model" yomi="ひせんけいしぐまもけい"} である。第二項の係数の半分は $vr(\partial n)^2$ と $-m_\sigma^2r^2/2$ を合わせた結果である。動径運動項は $r=\mathcal O(\partial^2)$ から六微分以上に寄与する。ここで $\mathcal O(\partial^6)$ は振幅のすべてのループを含む評価ではなく、動径場を古典的に消す微分展開の次数を表す。

例えば $n=(\pi/v,\sqrt{1-\pi^2/v^2})$ とすると
$$
\frac{v^2}2(\partial n)^2
=\frac12(\partial\pi)^2+\frac{(\pi\cdot\partial_\mu\pi)^2}{2(v^2-\pi^2)} .
$$
$\pi=0$ で運動項は正準で、四場の相互作用は $1/(2v^2)$ を係数とする二微分である。$N=2$ では角度そのものを使う座標 $n=(\sin(\varphi/v),\cos(\varphi/v))$ があり、二微分ラグランジアンが厳密に $(\partial\varphi)^2/2$ になる。従ってこの場合の最初の局所的な四点振幅は四微分から来る。その係数は $\mathcal L$ の $(\partial\varphi)^4/(2v^2m_\sigma^2)$ に由来し、質量ゼロの四脚を縮約すると $(s^2+t^2+u^2)/(v^2m_\sigma^2)$。元の線形模型で三つの交換を足した値と合う。

### 複数の角モードによる低エネルギー散乱

すべての外場の運動量を入射として $s=(p_1+p_2)^2$、$t=(p_1+p_3)^2$、$u=(p_1+p_4)^2$ とする。外場添字を $a,b,c,d$ として、$(\pi^2)^2$ を四回微分した結果の頂点は $-2i\lambda(\delta_{ab}\delta_{cd}+\delta_{ac}\delta_{bd}+\delta_{ad}\delta_{bc})$。$\sigma\pi_a\pi_b$ の頂点は $-2i\lambda v\delta_{ab}$。従って $s$ チャネルに付く添字部分の係数は
$$
-2\lambda-\frac{4\lambda^2v^2}{s-m_\sigma^2}
=\frac{s}{v^2}\frac{1}{1-s/m_\sigma^2}.
$$
同様に $t,u$ と対応する添字を交換する。先頭は
$$
\mathcal M_{ab\to cd}=
\frac{\delta_{ab}\delta_{cd}s+\delta_{ac}\delta_{bd}t+
\delta_{ad}\delta_{bc}u}{v^2}
+\mathcal O(p^4).
$$
例えば $\pi_1\pi_1\to\pi_2\pi_2$ では $s/v^2$、$\pi_1\pi_2\to\pi_1\pi_2$ では $t/v^2$。全ての添字が同じ場合は $s+t+u=0$ を使う。外脚を一つ軟らかくすると、残りの脚も質量殻に置いた運動学で $s,t,u$ は零へ近づく。この消え方は :term[Adler の零点]{en="Adler zero" yomi="あどらーのれいてん"} の木の例である。$h$ により質量を加えたときは同じ零をそのまま課せない。

### 明示的な破れを含む振幅の低エネルギー展開

源 $h>0$ を保って $\Phi_N=f+\sigma$、$\pi$ を横方向の場とする。$V_h$ の三次・四次の非微分頂点は前と同じ形で $v$ を $f$ に替える。一方 $\Phi^2-v^2$ の定数 $f^2-v^2$ が二次に効き、
$$
m_\pi^2=\lambda(f^2-v^2),\qquad
m_\sigma^2=m_\pi^2+2\lambda f^2 .
$$
例えば $V_h$ の二次を $\lambda[(f^2-v^2)(\sigma^2+\pi^2)/2+f^2\sigma^2]$ と展開すれば直接分かる。$s$ チャネルの $\delta_{ab}\delta_{cd}$ の係数は
$$
-2\lambda+\frac{4\lambda^2 f^2}{m_\sigma^2-s}
=\frac{s-m_\pi^2}{f^2}\frac{m_\sigma^2-m_\pi^2}{m_\sigma^2-s}.
$$
従って $|s|,m_\pi^2\ll m_\sigma^2$ で先頭は $(s-m_\pi^2)/f^2$。交差した $t,u$ の項にも同じ式を使う。物理的な二体閾値 $s=4m_\pi^2$ では $\pi_1\pi_1\to\pi_2\pi_2$ の先頭は $3m_\pi^2/f^2$、一方このチャネルの解析的な零は $s=m_\pi^2$ で閾値の下である。$h\to0$ で質量が消えると、無源のソフト極限に戻る。

:::example[二フレーバーのカイラル模型と同じ四成分]
$M=\sigma 1_2+i\pi_a\tau^a$ とし、Pauli 行列を $\tau^a$ とする。$\tfrac12\tr M^\dagger M=\sigma^2+\pi^2$、$\tfrac14\tr(\partial M^\dagger\partial M)=((\partial\sigma)^2+(\partial\pi)^2)/2$。$SU(2)_L\times SU(2)_R$ の作用 $M\mapsto LMR^\dagger$ はこの二つを保つ。$M=f1_2$ の安定化部分群は $L=R$、角方向は三つである。例えば $L=e^{i\alpha\tau^b/2}$、$R=e^{-i\alpha\tau^b/2}$ では $\delta M=i\alpha f\tau^b$ となり、真空が $\pi_b$ に動く。この表示を使うと、三本の角方向と発散が小さい軸性カレントの組を線形模型の中で具体化できる。これは強い相互作用を持つ微視的な真空をこの古典ポテンシャルだけで計算した結果ではない。
:::

## Goldstone の定理

### 保存流の相関に現れる極

連続な大域対称性の電流 $j^\mu$ が局所保存し、局所場 $O$ に対し $v_O=\langle\delta O\rangle\ne0$ とする。$\delta O=i[Q,O]$ の規約を取る。$F^\mu(p)=\int\dd^4x\,e^{ipx}\langle Tj^\mu(x)O(0)\rangle$ とすれば、発散は等時刻の交換子だけから来る。局所性により等時刻での台は $\vb x=0$。積分した係数は $\langle[Q,O]\rangle=-iv_O$ であり、Fourier 変換した微分の $-ip_\mu$ で割ると $p_\mu F^\mu(p)=v_O$ に局所な微分接触項が加わる。$O$ が Lorentz スカラーで真空が Lorentz 不変なら、共変な局所項は $p^\mu$ と $p^2$ の多項式であり、縮約には少なくとも $p^2$ が付く。従って縦の係数には $v_O/(p^2+i0)$ の極が必要である。質量ゼロの連続スペクトルがあり得る場合も、極に対応するスペクトルのデルタと連続な寄与は区別する。これはこの仮定での :term[Goldstone の定理]{en="Goldstone theorem" yomi="ごーるどすとーんのていり"} である。複数の生成子で $\langle\delta_a O_b\rangle$ のランクが $r$ なら質量ゼロの寄与を外積の和で表したとき少なくとも $r$ 項が必要になる。ゲージ固定に依存する場に同じ結論をそのまま適用するには、電荷と正定値空間の仮定を満たす必要がある。

### 極の重みと破れの行列のランク

質量ゼロの状態 $\ket{\pi_i(p)}$ に $\mel0{j_a^\mu(0)}{\pi_i(p)}=-i F_{ai}p^\mu$ と位相を選び、$\mel{\pi_i(p)}{O_b(0)}0=G_{ib}$ と置く。時間順序二点関数の極は $p^\mu\sum_i F_{ai}G_{ib}/(p^2+i0)$ となる。従って接触項を与える行列は $FG$ であり、$\rank(FG)\le\min(\rank F,\rank G)$。例えば $O(3)$ の真空を第三方向に取り $O_b=\pi_b$ とすれば、二つのカレントの線形項は $v\partial^\mu\pi_a$。木では $F=v1_2,G=1_2$ だからランクは $2$。零運動量の接触項が非零のとき、電流と演算子の重なりをどちらも持つ質量ゼロ成分が必要である。

### Wightman スペクトルによる質量ゼロ寄与の特定

時間順序を取らない二点関数の Fourier 変換を $W^\mu(p)=\int\dd^4x\,e^{ipx}\langle j^\mu(x)O(0)\rangle$ とする。真空の一項は引き、Lorentz 共変性と未来スペクトルを仮定する。非零運動量で唯一のベクトルは $p^\mu$ なので、$W^\mu$ は $p^\mu$ を掛けた質量のスペクトルに書ける。混合スペクトルは自己相関の正のスペクトル二つを用いて偏極恒等式を適用でき、局所有限な複素測度になると仮定する。この仮定ではデルタの微分を重みとする項はない。

局所保存から $p_\mu W^\mu=0$ である。$s=p^2$ の測度を $\rho(s)$ と書けば $s\rho(s)=0$。任意のコンパクトな台を持つ試験関数 $f$ が $s=0$ から離れているとき、$f/s$ も試験関数なので $\int f(s)\rho(s)\dd s=0$。従って測度は $s=0$ に集中し、$p^\mu\theta(p^0)\delta(p^2)$ 型が残る。真空の定数寄与は電流とスカラーの非零の交換子を作らないため、Ward 式の非零定数があればこの光的スペクトルの係数が非零である。スカラー補間場でこの光的軌道の重みを取り出すと、質量ゼロの状態が電流に結合する。これが混合相関での証明である。

:::example[小さな源を入れたときの極の移動]
$j^\mu\simeq v\partial^\mu\pi$、$\partial_\mu j^\mu=-h\pi$ を満たす先頭の線形模型では $h/v=m_\pi^2$。従って混合振幅は $F^\mu=v p^\mu/(p^2-m_\pi^2+i0)$、縮約は $v+v m_\pi^2/(p^2-m_\pi^2+i0)$ である。第一項は等時刻の接触項、第二項は明示的な非保存の項に対応する。$h=0$ のとき後者は消えるが、$F^\mu$ の質量ゼロ極は残る。
:::

### 局所電荷と Ward 恒等式の接触項

空間の滑らかなカットオフ $f_R(\vb x)$ を、原点付近で $1$、半径 $2R$ の外で零とする。電荷 $Q_R(t)=\int\dd^3x\,f_R(\vb x)j^0(t,\vb x)$ の微分は保存則から $\int\dd^3x\,\partial_i f_R\,j^i$。局所な $O(0)$ との交換子で考えると、固定した $t$ について十分大きい $R$ では $\partial_i f_R$ の台は $O$ から空間的に離れている。微視的因果性によりこの交換子は零になる。従って $\lim_{R\to\infty}\langle[Q_R(t),O]\rangle$ は $t$ に依存しない。この極限は局所場との交換子についての主張であり、破れた相の Hilbert 空間上で $Q_R$ 自体が収束することを要しない。

$T$ 積を時間で微分すると、二つの階段関数の微分が $\delta(x^0)\langle[j^0(x),O]\rangle$ を与える。空間積分に定数を使えば $\partial_i\delta^3(\vb x)$ などの接触項は零であり、電荷の作用が残る。相関に対称性を保つ局所対項を加えても、Lorentz スカラー $O$ の場合に許されるベクトル多項式は $p^\mu P(p^2)$ である。その縮約は $p^2P(p^2)$。従って Ward 恒等式の定数 $v_O$ は局所多項式だけで作れない。有限の $\langle\delta O\rangle$ があることを確認してはじめてこの極の議論が使える。

:::example[明示的な質量がある場合の電荷]
一成分の $\pi$ で $j^\mu=f\partial^\mu\pi$ と近似できるとする。$\pi$ の正準交換関係を使うと $[Q_R,\pi(0)]=-if$ は $f_R(0)=1$ のため変わらない。一方 $\partial_\mu j^\mu=-fm_\pi^2\pi$ なので、時間をずらした交換子は一般に一定でない。零空間運動量の自由場の交換子は $\cos(m_\pi t)$ に比例する。実際 $\partial_\mu j^\mu$ に由来する極が Ward 恒等式に加わり、$p_\mu F^\mu=f p^2/(p^2-m_\pi^2+i0)$ となる。この例では等時刻の非零の交換子だけから質量ゼロを結論すると誤る。
:::

### 低次元で仮定が成立するかを調べる

$1+1$ 次元の Lorentz 不変で正定値な場の理論を考え、連続な内部対称性を局所保存カレントで生成し、局所 Lorentz スカラーを秩序変数とする。通常の自発的破れがこの条件で起きないことを :term[Coleman の定理]{en="Coleman's theorem" yomi="こーるまんのていり"} という。実際、前項の Ward 式で秩序変数が非零の接触項を持つなら、混合相関には非零の光的スペクトルが要る。スペクトル射影 $E(B)$ と正定値性を使うと $|\langle j|E(B)|O\rangle|^2\le\langle j|E(B)|j\rangle\langle O|E(B)|O\rangle$ であり、$O$ 自身のスペクトルにも光的な台が必要になる。しかし局所 Lorentz スカラーの正の Wightman スペクトルは運動量空間の局所有限な測度である。$1+1$ 次元の非零の光的軌道上で正かつ Lorentz 不変な測度は、各枝で定数倍の $\dd p^1/|p^1|$。これは $p^1=0$ に近い各区間で積分が発散するので、非零の係数を持てない。離散の真空寄与 $p=0$ は先の電流の混合相関の接触項を再現しない。よって非零の秩序変数を仮定したことと矛盾する。

同じ赤外の問題は非線形模型の自由角の近似にも現れる。二次元 Euclid で $S=f^2\int(\partial\theta)^2/2$ を取り、零モードを分離して UV 正則化すると
$$
\langle(\theta(x)-\theta(0))^2\rangle
=\frac{2}{f^2}\int\frac{\dd^2p}{(2\pi)^2}
\frac{1-\cos(p\cdot x)}{p^2}
=\frac{1}{\pi f^2}\log(|x|/a)+\mathcal O(1)
$$
となる。最後は角度で平均した分子が $p|x|\ll1$ では二次から始まり、$1/|x|\lesssim p\lesssim1/a$ の $\dd p/p$ が対数を作ることから出る。Gaussian の指数挿入の相関は $\exp[-\langle(\theta(x)-\theta(0))^2\rangle/2]\propto |x|^{-1/(2\pi f^2)}$。例えば $f^2=1/(2\pi)$ ならこの近似では $1/|x|$。長距離で非零の秩序変数に対応する定数は残らない。この計算は障害を示す一例であり、非相対論的な理論、有限密度、局所電流の仮定を満たさない系の Goldstone モード数をこの式で決めることはできない。

### 壊れた方向の個数と非線形な実現

相の対称群を $H$ とし、秩序変数 $\Phi$ が $G$ の有限次元表現にある場合、$v_{ab}=\langle \delta_a\Phi_b\rangle$ は $T_a\Phi_0$ の行列である。その写像の核は $\mathfrak h$。従って $r=\dim G-\dim H$ 個の独立方向について、スカラー秩序変数の Ward 恒等式に質量ゼロの重みが必要である。線形模型では実際に $r$ 本の正準な角場を取れた。一般の場の理論でこの下限を論じる際は、真空の Lorentz 不変性、正定値性、連続で局所保存する電流を用いている。例えば有限密度で優先する時間方向があると、一つの極が複数の破れた生成子に対応する場合があるので、ここでの共変なスペクトル表示をそのまま使えない。

角方向を $G/H$ の局所座標 $\pi^a$ に取れば、$\Phi=\exp(\pi^a T_a/v)\Phi_0$ の原点で破れた変換は $\delta\pi^a=v\alpha^a+\mathcal O(\pi)$ となる。これは $\delta\Phi|_{\pi=0}=\alpha^aT_a\Phi_0$ と $\partial_{\pi^a}\Phi|_0=T_a\Phi_0/v$ を比較した。同じ変換でも固定する生成子は原点でゼロになり、$\pi$ の一次から作用する。この違いが、角多重項の作用に局所質量項を独立には加えられない理由である。例えば定数 $\pi$ を一つ選んだときのポテンシャルは、その軌道上で元と同じ値である。明示的破れがあればこの変換の定数項に対し作用が不変でなく、先の $h/f$ の質量を持てる。

## 有効ポテンシャルによる解析

連結生成汎関数の Legendre 変換で一様な有効ポテンシャルを求める。:eqref[eq-qft-effective-ward] から $V_i(T\Phi)^i=0$。非零の停留背景で微分すると $V_{ji}(T\Phi)^i=0$ なので、対称性で動く方向の零固有値が残る。運動項の係数が有限で非零なら逆二点関数の零運動量の零は質量ゼロに対応する。一次近似の曲率は零運動量であり、一般の有質量粒子の極を求めるには運動量依存の自己エネルギーも必要である。

例えば $Z_2$ の二極小を持つ理論の厳密な Legendre ポテンシャルは凸である。同じ低い値を持つ端点の間で、凸性は値が端点を超えないことを意味し、端点が大域最小なので区間は平坦になる。無限体積の一相を選んだ応答は $J\to0^+$ など片側の極限で扱う。外部源で相を選ぶ順序が重要になる理由である。

### 曲率と極位置の違い

単一揺らぎの Euclid 逆二点を $\Gamma^{(2)}(p)=a+Zp^2+cp^4/M^2+\cdots$ とする。ポテンシャルの曲率は $a$ だが、極は $p^2=-m_{\rm pole}^2$ で $a-Zm_{\rm pole}^2+cm_{\rm pole}^4/M^2=0$。$a/(ZM^2)\ll1$ なら先頭は $a/Z$、一次の補正は $ca^2/(Z^3M^2)$。例えば $a=2\,\mathrm{GeV}^2$、$Z=2$、$c=0$ なら極質量は $1\,\mathrm{GeV}$、曲率の平方根は $\sqrt2\,\mathrm{GeV}$。一方で対称性から $a=0$ となる角方向では、$Z$ が有限で非零なら零の位置が動かない。

### 凸性を生成汎関数から得る

Euclid で正の測度を持つ実スカラーの理論を考え、一様な実源 $J$ と $\overline\Phi=\mathcal V_4^{-1}\int\Phi$ を用いる。$w(J)=\mathcal V_4^{-1}\log Z[J]$ の微分は $\dd w/\dd J=\langle\overline\Phi\rangle_J$、二階微分は $\mathcal V_4\operatorname{Var}_J(\overline\Phi)\ge0$。複数成分のときは同様に共分散行列が半正定値である。Legendre–Fenchel 変換 $V_{\rm eff}(\varphi)=\sup_J[J\varphi-w(J)]$ はアフィン関数の上限なので凸となる。有限体積で $w''>0$ なら $w'(J)=\varphi$ を解けて $V_{\rm eff}''(\varphi)=1/w''(J)$。例えば二相の近似で $w''(0)=v^2\mathcal V_4$ なら対称な平均場での曲率は $1/(v^2\mathcal V_4)$ となる。

$\mathcal V_4\to\infty$ で二相の近似から $w(J)-w(0)\to v|J|$ を得る。$|\varphi|\le v$ なら $J\varphi-v|J|\le0$ で $J=0$ が上限を与えるため、$V_{\rm eff}$ は区間で一定である。$|\varphi|>v$ ではこの単純な近似の変換は無限大となる。実際には一相の中での分極があるので、その範囲にも有限なポテンシャルを持ち得る。木の非凸な二重井戸は局所背景のエネルギーを表す。無限体積の厳密な Legendre 変換の平坦部は、異なる相を混ぜて一様な平均場を実現できることに対応する。

### 大域模型の一ループ背景と停留条件

線形 $O(N)$ 模型で背景ノルムを $\varphi$ とする。二次ポテンシャルの Hessian は $\lambda(\varphi^2-v^2)\delta_{ij}+2\lambda\varphi_i\varphi_j$。動径方向で $M_r^2=\lambda(3\varphi^2-v^2)$、直交方向で $M_\perp^2=\lambda(\varphi^2-v^2)$ である。両質量二乗が正の領域で、:eqref[eq-qft-loop-master] から Gaussian の一ループ寄与は
$$
V_1=\frac{1}{64\pi^2}
\left[(M_r^2)^2\left(\log\frac{M_r^2}{\mu^2}-\frac32\right)
 +(N-1)(M_\perp^2)^2\left(\log\frac{M_\perp^2}{\mu^2}-\frac32\right)\right]
$$
となる。これは $\overline{\rm MS}$ で場に依存しない定数を省いた式である。各行列式 $\tfrac12\int\dd^dk/(2\pi)^d\,\log(k^2+M^2)$ を $M^2$ で微分し、極を引いた結果 $M^2(\log(M^2/\mu^2)-1)/(32\pi^2)$ を積分すれば係数と $3/2$ を再現できる。負の固有値の一様背景では Gaussian が安定せず、対数の解析接続を厳密な凸ポテンシャルの値と解釈できない。

$V_{\rm eff}=U(\Phi^2)$ と書けるなら $\partial_iV=2U'\Phi_i$、$V_{ij}=2U'\delta_{ij}+4U''\Phi_i\Phi_j$。従って角曲率は $V'(\varphi)/\varphi$。木の極小の $v$ で一ループだけの角曲率を評価すると、一般に $V_1'(v)/v$ が残る。正しい次数の停留点 $\varphi_*=v+\delta v$ を使い、$m_\sigma^2\delta v+V_1'(v)=0$ とすれば木の角曲率の変化がちょうどこれを打ち消す。この $V_1'$ の Goldstone 項は $M_\perp^2\log M_\perp^2\to0$ で極限を持つ。一方 $V_1''$ には $\log M_\perp^2$ が残り得る。角方向の Ward 恒等式と、動径方向の零運動量の赤外特異性は別の問題である。

### 停留点をずらす数値例と有効展開

$\varphi=v$ で $M_r^2=m_\sigma^2$ とおき、$d M_r^2/d\varphi=6\lambda\varphi$ を使う。角部分の $V_1'$ の極限は零なので
$$
V_1'(v)=\frac{6\lambda v\,m_\sigma^2}{32\pi^2}
 \left(\log\frac{m_\sigma^2}{\mu^2}-1\right),
\qquad
\delta v=-\frac{3\lambda v}{16\pi^2}
 \left(\log\frac{m_\sigma^2}{\mu^2}-1\right).
$$
例えば $\mu=m_\sigma$ に選ぶと $\delta v/v=3\lambda/(16\pi^2)$。$\lambda=1/2$ なら約 $0.0095$。この例は質量・結合を同じ $\overline{\rm MS}$ パラメータで比較し、非凸な領域の揺らぎを積分した値ではない。角曲率の木の移動は $2\lambda v\delta v$ であり、$2\lambda v=m_\sigma^2/v$ を使うと一ループの角曲率 $V_1'(v)/v$ を消す。この検査では停留条件を次数ごとに展開しているため、数値的な曲率を単一の切断ポテンシャルで最小化して高い次数を混ぜる必要はない。

真空に近い微分展開で $\Gamma$ の二次の運動項を $Z(\varphi)(\partial\Phi)^2/2$ と書くと、動径・角方向で同じ係数になる保証はない。例えば大域 $O(N)$ 不変な追加項 $Y(\Phi^2)(\Phi\cdot\partial_\mu\Phi)^2/2$ は、背景 $(0,\dots,v)$ で動径の $Z$ だけに $Y(v^2)v^2$ を加える。従って大域 Ward 恒等式が角の零運動量の曲率を固定しても、波動関数の規格化まで同じにはしない。電流の留数 $f$ を測るには、正準規格化した励起と電流の行列要素を使う。角方向の運動項が $Z_\perp$ で、線形変換が $v$ の移動を持つとき、その先頭の電流は $Z_\perp v\partial\pi$。正準角場 $\pi_{\rm c}=\sqrt{Z_\perp}\pi$ に対する係数は $\sqrt{Z_\perp}v$ となる。

## Higgs 機構

### 角方向をゲージ場に組み込む

複素スカラー $\phi=(v+h)e^{i\theta/v}/\sqrt2$ と $D_\mu=\partial_\mu-igA_\mu$ を使う。$U=e^{ig\omega}$ で $\theta\mapsto\theta+gv\omega$、$A_\mu\mapsto A_\mu+\partial_\mu\omega$。従って $\theta$ の微分と $A$ は
$$
|D\phi|^2=\tfrac12(\partial h)^2+
\tfrac12(v+h)^2\left(\frac{\partial_\mu\theta}{v}-gA_\mu\right)^2
\label{eq-qft-higgs-kinetic}
$$
になる。摂動真空の近傍では $\theta=0$ を選べる。すると $m_A=gv$ の Proca 質量、$g^2v h A_\mu A^\mu$ の三点、$g^2h^2A_\mu A^\mu/2$ の四点が得られる。二つの質量ゼロ横偏極と二つの実スカラーの計四自由度は、三偏極の有質量ベクトルと一つの動径スカラーになる。この仕組みを :term[Higgs 機構]{en="Higgs mechanism" yomi="ひっぐすきこう"} と呼ぶ。局所ゲージ変換は記述の冗長性であり、ゲージに依存する $\langle\phi\rangle$ 単独を物理的な秩序変数とみなすことはできない。

非可換群で $\Phi_0$ が背景なら、複素場の質量行列は $(M^2)_{ab}=2g^2\Re[(T^a\Phi_0)^\dagger(T^b\Phi_0)]$。二次形式が作用する生成子のノルムなので半正定値であり、核は $\Phi_0$ を固定する生成子である。例として $SU(2)$ の実随伴場を第三方向に置けば $T^3$ だけが背景を固定し、他の二方向が質量を得る。

### 基本表現の例と随伴表現の違い

$SU(2)$ の基本複素二成分を $\Phi_0=(0,v/\sqrt2)^{\mathsf T}$ に置き、$T^a=\sigma^a/2$ とする。$\Re[(T^a\Phi_0)^\dagger(T^b\Phi_0)]=\Phi_0^\dagger\{T^a,T^b\}\Phi_0/2=v^2\delta^{ab}/8$。従って $M^2_{ab}=g^2v^2\delta^{ab}/4$ であり、三つすべてが質量 $gv/2$ を持つ。随伴場の第三成分を選んだ前の例と違い、$T^3\Phi_0$ も零でないためである。さらに別の $U(1)$ もゲージ化した場合、二つの対角生成子の線形結合が $\Phi_0$ を消し得る。その判定には群の名称でなく、質量行列を背景上で実際に作用させる。

### 二つの対角ゲージ場の混合

$SU(2)\times U(1)$ で複素二重項の $U(1)$ 電荷を $1/2$、背景を $(0,v/\sqrt2)^\mathsf T$ とする。下成分の $T^3=-1/2$ により中性部分の質量項は $v^2(gA_\mu^3-g'B_\mu)^2/8$。従って非零の固有値は $(g^2+g'^2)v^2/4$、零固有方向は $(g',g)$ である。荷電の二方向には $g^2v^2/4$。例えば $g=4/5$、$g'=3/5$ なら中性の質量は $v/2$、荷電の質量は $2v/5$、比は $5/4$ である。この計算は一つの二重項の正準運動項だけを仮定した木の結果である。

### ゲージ不変な変数で自由度を追う

Abel 模型で $\phi\ne0$ の局所領域に限り、$C_\mu=A_\mu-\partial_\mu\theta/(gv)$ と置ける。$F(A)=F(C)$、:eqref[eq-qft-higgs-kinetic] は $g^2(v+h)^2C_\mu C^\mu/2+(\partial h)^2/2$。二次の $C$ の方程式 $\partial_\nu F^{\nu\mu}+m_A^2C^\mu=0$ を微分すると $m_A^2\partial_\mu C^\mu=0$。$m_A^2>0$ では四成分の一つが制約を受け、残る三成分が質量殻を持つ。例えば静止系 $k^\mu=(m_A,\vb0)$ で $C^0=0$、三つの空間偏極が独立である。角変数の記述は $\phi=0$ の点や巻き付きのある場を覆い尽くさないため、この変数だけから大域的なゲージ固定の存在を仮定できない。

ゲージ不変な局所演算子を作るには、例えば $\phi^\dagger\phi$ と
$$
\mathcal J_\mu=i[\phi^\dagger D_\mu\phi-(D_\mu\phi)^\dagger\phi]
=g(v+h)^2 C_\mu
$$
を取る。展開の先頭で $\phi^\dagger\phi=v^2/2+vh+\cdots$、$\mathcal J_\mu=gv^2 C_\mu+\cdots$ だから、連結二点に動径粒子とベクトルの極が現れる。ゲージ固定で得た質量を、このようにゲージ不変相関の極と照合できる。例えば有限格子で一つのサイトに独立な Abel 変換を行い Haar 測度で平均すると、荷電場の単独の期待値は $\int\dd\alpha\,e^{i\alpha}/(2\pi)=0$。$\phi^\dagger\phi$ と $\mathcal J$ に相当するゲージ不変な挿入は同じ平均で消えない。この局所平均が、ゲージ固定前の荷電場を大域対称性の秩序変数として使えない具体的な理由である。

### 一般の実スカラーと残るゲージ群

実表現を反対称行列 $t^a$ で表し、$D_\mu\Phi=\partial_\mu\Phi-gA_\mu^at^a\Phi$、運動項を $(D\Phi)^2/2$ とする。質量行列の成分は $g^2(t^a\Phi_0)\cdot(t^b\Phi_0)$。複素表示との係数二の違いは運動項の $1/2$ に由来する。実随伴 $SU(2)$ の $t^a_{ij}=-\epsilon_{aij}$ で $\Phi_0=(0,0,v)$ を取ると、二方向のノルムは $v$、第三方向は零。従って固有値は $g^2v^2,g^2v^2,0$ である。$v=0$ の極限では三つとも零へ戻る。基本複素場の例で第三方向も質量を持つことと比較すれば、破れる生成子を Lie 代数の次元だけで決められないことが分かる。

二重項の例では $s_W=g'/\sqrt{g^2+g'^2}$、$c_W=g/\sqrt{g^2+g'^2}$ とおく。正準運動項を保つ直交回転 $Z_\mu=c_W A_\mu^3-s_W B_\mu$、$A_\mu^{\rm em}=s_W A_\mu^3+c_W B_\mu$ により、$A_\mu^{\rm em}$ の結合は $g s_W T^3+g' c_W Y=e(T^3+Y)$、$e=g s_W=g'c_W$ になる。真空の下成分では $T^3+Y=0$。例えば同じ二重項の上成分は $T^3+Y=1$ なので質量ゼロ場に電荷 $e$ で結合する。この例の $m_W/m_Z=c_W$ は正準規格化と二重項の背景から従う木の関係である。

### Higgs 機構に電荷や真空を追加する例

Abel の二場に電荷 $q_1,q_2$、実の背景 $v_1/\sqrt2,v_2/\sqrt2$ を取る。位相を正準角場 $\theta_1,\theta_2$ で表すと、二次の混合は $-g A^\mu\partial_\mu(q_1v_1\theta_1+q_2v_2\theta_2)$、ベクトル質量二乗は $m_A^2=g^2(q_1^2v_1^2+q_2^2v_2^2)$。従って角場のうちゲージ変換で動く組は $\chi_\parallel=(q_1v_1\theta_1+q_2v_2\theta_2)/v_q$、$v_q=(q_1^2v_1^2+q_2^2v_2^2)^{1/2}$。直交する $\chi_\perp=(-q_2v_2\theta_1+q_1v_1\theta_2)/v_q$ は不変である。例えば $q_1=q_2=1$、$v_1=v_2=v$ のとき動くのは和、残る差をゲージ変換で除けない。もし作用が位相差の大域的シフトを保つなら、この残りに Goldstone の機構が働く。

同じ模型で $-\kappa(\phi_1^\dagger\phi_2+\mathrm{h.c.})$ をポテンシャルに加えると、位相差のシフトが明示的に破れる。$\kappa>0$ とし、背景の大きさは別の半径ポテンシャルで選んだものを使う。位相部分は $-\kappa v^2\cos[(\theta_2-\theta_1)/v]$。二次は $\kappa(\theta_2-\theta_1)^2/2=\kappa\chi_\perp^2$、角差の質量二乗は $2\kappa$ となる。和にはポテンシャルのこの二次項がなく、$m_A$ は運動項で決まる。この例では局所ゲージの角と大域的な物理角を一つずつ持ち、二つを同じ「消える角」として数えない。

### 非可換ゲージ固定で質量行列を対角化する

実成分のスカラー揺らぎを $\delta\Phi^i$ とし、混合係数の行列を $K_{ia}=g(t^a\Phi_0)_i$ とする。ベクトル質量は $M_A^2=K^{\mathsf T}K$。直交変換で $K$ の特異値 $m_a>0$ を対角形にすれば、対応するスカラーの基底は $K e_a/m_a$ となり、内積は $\delta_{ab}$ である。したがって混合を $-m_a A_\mu^a\partial^\mu\chi^a$ の形にできる。未破れ方向は $K e_a=0$ なので、そもそも相手となる $\chi^a$ を持たない。例えば実随伴の $SU(2)$ では $K$ のランクは二で、二つだけの $\chi$ が角に対応する。一般の $R_\xi$ の設計では、この質量基底の各非零ペアに $F^a=\partial\cdot A^a+\xi m_a\chi^a$ を用いる。

## R_ξ ゲージと Goldstone ボソン等価定理

線形座標 $\phi=(v+h+i\chi)/\sqrt2$ では二次の混合が $-m_A A^\mu\partial_\mu\chi$。従って
$$
\mathcal L_{\rm gf}=-\frac1{2\xi}(\partial\cdot A+\xi m_A\chi)^2
\label{eq-qft-rxi}
$$
の交差項が部分積分後の混合 $+m_A\chi\partial\cdot A$ を消す。$\chi$ の質量二乗は $\xi m_A^2$ となる。また $\delta\chi=g(v+h)\omega$ なのでゴースト演算子は $\Box+\xi m_Ag(v+h)$。Abel でもこのゲージではゴーストが $h$ と結合する。自由ベクトル演算子では横固有値が $-(k^2-m_A^2)$、縦固有値が $-(k^2-\xi m_A^2)/\xi$ であり、逆に $i$ を掛けて伝播関数を得る。
$$
D_{\mu\nu}(k)=-i\left[
\frac{\eta_{\mu\nu}-k_\mu k_\nu/k^2}{k^2-m_A^2+i0}
+\frac{\xi k_\mu k_\nu/k^2}{k^2-\xi m_A^2+i0}\right].
$$

高エネルギーで縦偏極を角方向へ置き換えるには、偏極ベクトルの近似に加えて振幅の Ward 恒等式が必要である。以下ではまず $\xi=1$ の木で両外脚の質量殻を一致させ、位相を含めて導く。一般の $\xi$ やループでは外脚の規格化と運動量の指定が必要である。

### 縦偏極の高エネルギー展開

$k=(E,0,0,p)$ なら $\epsilon_L=(p,0,0,E)/m_A$。差は
$$
\epsilon_L-\frac{k}{m_A}
=\frac{E-p}{m_A}(-1,0,0,1)
=\frac{m_A}{E+p}(-1,0,0,1).
$$
例えば $E=5m_A$ では係数は $1/(5+\sqrt{24})\simeq0.101$。振幅での補正はこの差を打ち切った $\mathcal M_\mu$ に縮約したものである。具体的に $k\cdot\mathcal M=0$ で角モード振幅が零なら、縦偏極の振幅はこの差から始まる。その場合、角モードとの比による誤差評価はできないが、どの成分が抑制を受けるかはこの式から読める。

### ゴースト質量と非物理極の打消し

線形場の BRST は $s\chi=g(v+h)c$、$sh=-g\chi c$、$sA_\mu=\partial_\mu c$。従って $sF=(\Box+\xi m_A^2+\xi m_A gh)c$ であり、$\mathcal L_{\rm gh}=-\bar c\,sF$ からゴーストの極も $k^2=\xi m_A^2$ にある。例えば $\xi=1$ では $D_{\mu\nu}=-i\eta_{\mu\nu}/(k^2-m_A^2+i0)$。$\xi=0$ は $B$ で $\partial\cdot A=0$ を課す極限として扱い、$\chi,c$ は木で質量ゼロになる。$\xi$ を有限の正値に固定することと、最初から $\xi\to\infty$ を取ることでは、大運動量での伝播関数の振る舞いが異なる。

自由な $C_\mu=A_\mu-\partial_\mu\chi/m_A$ を取ると、$A$ と $\chi$ の混合は消えているので相関は二項の和である。Goldstone 伝播関数 $i/(k^2-\xi m_A^2+i0)$ に $k_\mu k_\nu/m_A^2$ が掛かり、
$$
\begin{aligned}
D_{\mu\nu}^C
&=\frac{-iP^T_{\mu\nu}}{k^2-m_A^2+i0}
+\frac{i(-\xi+k^2/m_A^2) P^L_{\mu\nu}}{k^2-\xi m_A^2+i0}\\
&=\frac{-i}{k^2-m_A^2+i0}
\left(\eta_{\mu\nu}-\frac{k_\mu k_\nu}{m_A^2}\right)
\end{aligned}
$$
となる。$k^2=\xi m_A^2$ の極は分子の零により消え、縦射影の係数 $i/m_A^2$ は物理ベクトルの殻には極を持たない。例えば静止系の殻で $C_0C_0$ の留数は零、空間成分だけが残る。この打消しは自由なゲージ不変挿入の具体例であり、相互作用のある一般の外脚には ST 恒等式を用いる。

### 等価定理で使う一脚の恒等式

消去前のゲージ固定関数 $F=\partial\cdot A+\xi m_A\chi$ は $F=-\xi s\bar c$（補助場の方程式を使う挿入式）である。その他の外脚を物理的な漸近状態へ簡約すれば、離れた $F$ の行列要素は零。接触項は最終的な物理外脚の留数を取る際に扱い、以下では木の一脚の極を取る。外部 $x$ から運動量 $k$ が外へ向かうとき、Fourier 位相 $e^{-ikx}$ では微分が $-ik$ を与える。打ち切った振幅を $T^\nu,T_\chi$ とし、$\xi=1$ で $D_{\mu\nu}=-i\eta_{\mu\nu}/(k^2-m_A^2+i0)$、$D_\chi=i/(k^2-m_A^2+i0)$ を掛けると $F$ の留数は
$$
-k_\nu T^\nu+i m_A T_\chi .
$$
ゆえに殻上で $k_\nu T^\nu/m_A=iT_\chi$。外脚の向きや $\chi$ の定義でこの位相は変わる。$v^\mu(k)=\epsilon_L^\mu-k^\mu/m_A$ とおくと
$$
\mathcal M(A_L)=i\mathcal M(\chi)+v_\mu(k)T^\mu .
$$
最後の項は $v^\mu=\mathcal O(m_A/E)$ のテンソル縮約であり、$\mathcal M(\chi)$ が零や小さい場合に比として評価できない。複数の縦脚では、他の脚を非物理な $\chi$ に置いた一脚の Ward 式を独立に反復するだけでは証明できない。BRST 完全な漸近挿入を組み合わせた恒等式を用いると、先頭の置換項と、少なくとも一つの $v^\mu$ を含む残項に分かれる。複数脚の等価定理は :cite[peskin-schroeder] を参照する。ループまでの式は混合する外脚の留数の解析も要する。

一般の有限 $\xi$ で $k_\mu D^{\mu\nu}$ と $D_\chi$ を用いると、$F$ 挿入から直接読む不物理極は $k^2=\xi m_A^2$ である。物理ベクトルの殻は $k^2=m_A^2$ なので、同じ残差の式に二種類の殻を混ぜることはできない。一般の定理ではベクトルの運動量で評価した打切りスカラー挿入と有限な係数を用いる。例えば位相を吸収して比較する木の $\xi=1$ では両殻が一致し、上の一脚の導出でその係数の絶対値は一になる。

### 二本の縦偏極が出る木の崩壊

Abel Higgs ポテンシャルを $\lambda(|\phi|^2-v^2/2)^2$ とし、$m_h^2=2\lambda v^2$、$m_A=gv$ とする。$hAA$ は $2ig^2v\,\eta_{\mu\nu}$、$h\chi\chi$ は $-2i\lambda v$。$h$ の静止系で、二本のベクトルを $k_{1,2}=(E,0,0,\pm p)$、$E=m_h/2$ とする。縦偏極は $\epsilon_{1,2}=(p,0,0,\pm E)/m_A$ なので $\epsilon_1\cdot\epsilon_2=(p^2+E^2)/m_A^2=(m_h^2-2m_A^2)/(2m_A^2)$。従って $i$ を共通に落とした振幅は
$$
\mathcal M(h\to A_L A_L)=\frac{m_h^2-2m_A^2}{v},\qquad
\mathcal M(h\to\chi\chi)=-\frac{m_h^2}{v}.
$$
木の $\xi=1$ なら二種の位相空間も同じであり、二本の外向きの $i$ は積で $-1$。$m_h\gg m_A$ では先頭が一致し、この過程での相対振幅補正は $2m_A^2/m_h^2$。例えば $m_h=5m_A$ なら補正は $8\%$、同じ位相空間で振幅の二乗の比は $(23/25)^2$。比較している物理的な振幅は縦偏極二本の部分であり、横偏極も足した全崩壊幅の等式ではない。

運動学が十分高くても、特定の偏極や励起の先頭項が対称性で消える場合には、残差項を計算しなければならない。例えば $hAA$ の結合で一方を縦、もう一方を同じ軸に垂直な偏極に取れば、この静止系の計量縮約は零になる。$\chi$ の散乱だけから任意の偏極の振幅を復元できるわけではない。等価定理は縦の外脚に限ることが、偏極ベクトルの近似とこの例の双方から分かる。

## 自発的に破れたゲージ理論のくりこみ可能性

有限の $\xi$ のベクトル・スカラー・ゴースト伝播関数は高運動量で $1/k^2$ として減衰する。元のゲージ作用は質量次元四以下であり、線形場を $v$ だけずらしても新しい項の次元は増えない。従って :eqref[eq-qft-power-counting] の局所対項の数え上げを適用できる。一方、ユニタリゲージの :eqref[eq-qft-proca-propagator] では縦分子が $k^2/m_A^2$ で増えるので、同じ個別グラフの数え上げを適用できない。物理振幅の関係は別のゲージで BRST 恒等式を保って示す。

:::theorem[破れたゲージ理論の摂動的なくりこみ]{#thm-qft-broken-renormalization}
四次元で次元四以下の Yang–Mills・物質作用を真空の周りで展開し、ゲージアノマリーが無いと仮定する。$R_\xi$ ゲージの摂動論では、局所な対項を選んで Slavnov–Taylor 恒等式を各次数で満たせる。
:::

:::proof
数え上げと量子作用原理により、低い次数で恒等式を回復した後に残る破れは局所な挿入となる。以下で示す帰納法では、この挿入が線形化した BRST 微分の像に属することが決定的である。アノマリーの分類と量子作用原理の一般証明には全次数での局所挿入のくりこみが必要なので :cite[kugo-gauge] と :cite[weinberg-qft2] を参照する。
:::

### 真空の移動と対項の連動

Abel 模型で $\phi=(v+h+i\chi)/\sqrt2$ とする。場と結合のくりこみに加え $v\mapsto v+\delta v$ とすると、同じ $|D\phi|^2$ の平方から $A^2$ の項は $g^2v\,\delta v A^2$、$hA^2$ の項は $g^2\delta v\,hA^2$ を受ける。例えば $g$ を固定した微小変化では $m_A^2=g^2v^2$ の相対変化は $2\delta v/v$ に対し、三点係数 $g^2v$ は $\delta v/v$。頂点ごとにこれらを独立に選ぶと、もとのゲージ不変な運動項へ戻せない。真空を定める条件は一次の $h$ の係数を消す条件として課しつつ、同じパラメータ変化を二次以上の頂点にも入れる。

### 線形化した恒等式を対項で回復する

非線形 BRST 変換には源を入れ、補助場は消す前に扱う。古典作用を $S$、古典作用の周りで線形化した Slavnov 演算子を $\mathcal B_S$ とする。:eqref[eq-qft-slavnov-taylor] と $s^2=0$ により、場・源対を含めた局所汎関数上で $\mathcal B_S^2=0$。この性質は Slavnov 汎関数を場と源の反括弧の自己積および補助場の線形部分で書き、Jacobi 恒等式により二重作用を $\mathcal S(S)$ に比例させると得られる。そこで $n-1$ ループまで恒等式を満たし、$n$ で
$$
\mathcal S(\Gamma)=\hbar^n\Delta_n+\mathcal O(\hbar^{n+1})
$$
となったとする。Slavnov 汎関数にさらに線形化演算子を作用させる同じ恒等式から、最低次数で $\mathcal B_S\Delta_n=0$。$\Delta_n$ はゴースト数一の局所汎関数である。ゲージアノマリーが無く、許された対称性・作用の条件の下でこの閉項が $\Delta_n=\mathcal B_S X_n$ と書ければ、対項 $-\hbar^n X_n$ を加えた差は $-\hbar^n\mathcal B_S X_n$。よって $n$ 次でも回復できる。二次以上の積は $n$ 次で寄与せず、帰納が閉じる。

回復後に残る対項の自由度はゴースト数零で $\mathcal B_S$ に閉じる。ゲージ不変な物理パラメータの変更、場と変換源の再定義、ゲージ固定の変更に対応する項がある。例えば二重項の湯川結合を持つ模型で、背景のシフトから出るフェルミオン質量だけを独立に増やすと、元のゲージ不変な $\bar\psi_L\Phi\psi_R$ の変分に必要な三点対項が欠ける。場をシフトする前の不変な演算子で対項を組み立てれば、質量項と角場の結合が同時に変わる。この連動は有限の $\xi$ で伝播関数が $1/k^2$ に落ちるという事実だけからは従わない。

### ゲージ依存な座標と物理量

ゲージ固定・ゴーストを $s\Psi_\xi$ と書き、BRST 変換と局所測度が保たれる範囲で $\xi$ を微小に変える。普通の場の源を入れた積分で $s(\partial_\xi\Psi_\xi)$ の挿入を変数変換により取り除くと、$\partial_\xi\Gamma$ は $\Gamma$ の場微分と、BRST 変換・$\partial_\xi\Psi_\xi$ の混合挿入の積に比例する。一様背景一変数の場合、符号を係数 $C$ に含めて
$$
\partial_\xi V(\varphi,\xi)+C(\varphi,\xi)\,\partial_\varphi V(\varphi,\xi)=0
$$
となる。これを :term[Nielsen 恒等式]{en="Nielsen identity" yomi="にーるせんこうとうしき"} と呼ぶ。導出では源付き積分の変数変換を用いているので、ゲージに依存する背景の値そのものを固定して比較する議論ではない。

非退化な停留点 $\varphi_*(\xi)$ で $C$ が有限なら、点の値の全微分は $0$ だが、停留条件を $\xi$ で微分すると $0=\partial_\xi V'+V''\dd\varphi_*/\dd\xi$。恒等式を $\varphi$ で微分して $\partial_\xi V'=-C V''$ を代入すると $\dd\varphi_*/\dd\xi=C$ となり、位置は動き得る。もう一度 $\varphi$ 微分すると、停留点に沿う曲率の微分は $-2C' V''$ を含む。従って単純な曲率をゲージ不変な質量とすることはできない。例えば $V(\varphi,\xi)=U(\varphi-a\xi)$ なら $C=a$、極小の値は一定で位置だけが $a\xi$ だけ移る。複数の場が混合する物理的な極の議論では二点関数全体の恒等式を使う。

### Nielsen 恒等式を源の変分から求める

符号の例を一つ固定する。普通の偶な場 $\Phi_i$ の源を $\int J_i\Phi_i$、$\int\mathcal D\Phi\,e^{i(S+J\Phi)}=e^{iW}$ とし、$\Psi_\xi$ は奇とする。$\partial_\xi W=\langle s(\partial_\xi\Psi_\xi)\rangle_J$。源のある積分で $(\partial_\xi\Psi_\xi)e^{i(S+J\Phi)}$ に $s$ を作用させて積分すると、第二項には奇の Leibniz 則の負号が付く。従って
$$
\partial_\xi W=i\langle(\partial_\xi\Psi_\xi)J_i\,s\Phi_i\rangle_J
$$
となる。Legendre 変換では古典場固定の $\partial_\xi\Gamma=\partial_\xi W$ であり、$J_i=-\delta\Gamma/\delta\Phi_i$。よって古典場の微分を係数とする恒等式が得られる。全ての場に源を入れるときは奇源の順を固定し、物理的な一様背景に制限したところでポテンシャルの式を得る。局所な混合挿入のくりこみが必要な場合は、BRST を保つ源の対項も含めてこの変数変換に合わせる。

### 高運動量の数え上げを二つのゲージで照合する

$R_\xi$ でスカラー・ベクトル・ゴーストの次元を一、Dirac を $3/2$ とする。頂点 $v$ が場と微分を合わせて次元 $d_v\le4$ なら、連結 1PI グラフの表面的次数は
$$
\omega=4-E_A-E_\phi-E_c-\tfrac32 E_\psi+\sum_v(d_v-4)
$$
以下である。$\Phi\mapsto v+\delta\Phi$ は一つの場を定数に替えるので、誘導される頂点の $d_v$ は増えない。例えば $hA^2$ の場の次元は三、$h^2A^2$ は四。前者のみを含む六点図に高次の局所対項が必要になることはなく、部分発散の方を低点の対項で除く。一方、:eqref[eq-qft-proca-propagator] の縦成分は大運動量で $iP^L/m_A^2$ だから、その一線につきこの数え上げより二次数高くなる。ユニタリゲージの個々のオフシェル相関で必要な引き算と、$R_\xi$ で定義し BRST を保った物理振幅のくりこみの主張は同じではない。

Nielsen 恒等式を一般の有効作用に使うと、停留背景と共に移動した二点の行列 $\Gamma^{(2)}$ のゲージ微分は、他の脚への混合が無い場合 $\partial_\xi\Gamma^{(2)}=-C'^{\mathsf T}\Gamma^{(2)}-\Gamma^{(2)} C'$ の形になる。これは $\partial_\xi\Gamma=-\Gamma_i C^i$ を二回微分し、第三微分の項を背景の移動 $\dd\Phi_*^i/\dd\xi=C^i$ で相殺した式である。極の近傍で $C'$ が正則なら、これは正則な行・列の変換なので物理的な零の位置は変わらない。質量殻の付近に混合した特異性がある場合は、適切な物理部分空間と留数を同時に扱う必要がある。

## 演習問題

::::exercise{#exr-qft-symmetry-breaking-1 level=basic}
線形 $O(4)$ 模型の非零真空で角方向の個数を求めよ。

:::solution
真空の単位方向を固定すると残る回転は $O(3)$。真空球面 $S^3$ の接空間は三次元なので三つである。
:::
::::

::::exercise{#exr-qft-symmetry-breaking-2 level=basic}
$\lambda=1/2$、$v=2\,\mathrm{GeV}$ の動径質量を求めよ。

:::solution
:eqref[eq-qft-sigma-expansion] より $m_\sigma=\sqrt{2\lambda v^2}=2\,\mathrm{GeV}$ である。
:::
::::

::::exercise{#exr-qft-symmetry-breaking-3 level=basic}
$g=1/2$、$v=4\,\mathrm{GeV}$ のゲージ場の質量を求めよ。

:::solution
:eqref[eq-qft-higgs-kinetic] の質量項は $g^2v^2 A^2/2$ なので $m_A=gv=2\,\mathrm{GeV}$ である。
:::
::::

::::exercise{#exr-qft-symmetry-breaking-4 level=basic}
$\xi=4$、$m_A=2\,\mathrm{GeV}$ で $\chi$ の質量を求めよ。

:::solution
質量二乗は $\xi m_A^2=16\,\mathrm{GeV}^2$ なので質量は $4\,\mathrm{GeV}$ である。
:::
::::

::::exercise{#exr-qft-symmetry-breaking-5 level=advanced}
有効ポテンシャルの Ward 恒等式を微分し零固有値を導け。

:::solution
$V_i(T\Phi)^i=0$ の微分は $V_{ji}(T\Phi)^i+V_iT^i{}_j=0$。停留点では第二項が零なので、そこを動かす $T\Phi$ が Hessian の核に入る。
:::
::::

::::exercise{#exr-qft-symmetry-breaking-6 level=advanced}
:eqref[eq-qft-rxi] の交差項と運動項の混合が消えることを示せ。

:::solution
固定項の交差は $-m_A\chi\partial\cdot A$。運動項の $-m_A A^\mu\partial_\mu\chi$ を部分積分すると $+m_A\chi\partial\cdot A$ なので消える。
:::
::::

::::exercise{#exr-qft-symmetry-breaking-7 level=advanced}
質量行列が半正定値であることを示せ。

:::solution
実係数 $u_a$ に対し $u_a(M^2)_{ab}u_b=2g^2\|\sum_a u_aT^a\Phi_0\|^2$ である。従って負の固有値はなく、核は生成子が背景を消す方向である。
:::
::::

::::exercise{#exr-qft-symmetry-breaking-8 level=advanced}
$\epsilon_L^\mu-k^\mu/m_A$ の各非零成分の高エネルギー次数を求めよ。

:::solution
$p$ を $z$ 方向とすれば差は $((p-E)/m_A,0,0,(E-p)/m_A)$。$E-p=m_A^2/(E+p)$ なので絶対値は $m_A/(E+p)$、すなわち $\mathcal O(m_A/E)$ である。
:::
::::

::::exercise{#exr-qft-symmetry-breaking-9 level=advanced}
複素場 $\phi=(v+h+i\chi)/\sqrt2$ の Abel Higgs 模型で、ゲージ不変なポテンシャルを $\lambda(|\phi|^2-v^2/2)^2$ とする。質量二乗 $m_h^2,m_A^2,m_\chi^2$ を $R_\xi$ ゲージで求めよ。$\lambda=1/8,g=1/2,\xi=4,v=2\,\mathrm{GeV}$ では値はいくらか。

:::solution
$|\phi|^2-v^2/2=vh+(h^2+\chi^2)/2$ なのでポテンシャルの二次は $\lambda v^2h^2$、従って $m_h^2=2\lambda v^2$。:eqref[eq-qft-higgs-kinetic] の二次から $m_A^2=g^2v^2$。$\chi$ に元のポテンシャルからの二次はなく、:eqref[eq-qft-rxi] の二次項が $m_\chi^2=\xi g^2v^2$ を与える。指定値ではそれぞれ $1,1,4\,\mathrm{GeV}^2$。$h$ と物理ベクトルの質量はこの木の例では等しいが、$\chi$ の極は $\xi$ を変えると動く。ゲージ固定に伴うスカラー極を観測可能な質量と数えると自由度の数も変わってしまう。
:::
::::

::::exercise{#exr-qft-symmetry-breaking-10 level=advanced}
電荷が $q_1=1,q_2=2$、背景がともに $v/\sqrt2$ である二場の Abel 模型を考える。$\phi_1^2\phi_2^\dagger$ はゲージ不変であることを示し、ポテンシャルに $-\kappa(\phi_1^2\phi_2^\dagger+\mathrm{h.c.})$、$\kappa>0$ を加えたとき角方向の物理質量を求めよ。半径ポテンシャルはこの背景が停留するよう選ぶ。

:::solution
電荷は二つの $\phi_1$ が $2$、$\phi_2^\dagger$ が $-2$ なので和が零。正準角場を $\theta_1,\theta_2$ とすると、相互作用は一様背景で
$$
-\frac{\kappa v^3}{\sqrt2}\cos[(2\theta_1-\theta_2)/v].
$$
ベクトルに混合する角は $(\theta_1+2\theta_2)/\sqrt5$、物理角は $(2\theta_1-\theta_2)/\sqrt5$ と取れる。余弦の二次に代入すると物理角の二次は $5\kappa v\,\chi_\perp^2/(2\sqrt2)$、従って $m_\perp^2=5\kappa v/\sqrt2$。一方 $m_A^2=5g^2v^2$。$\kappa$ は次元一、$g$ は無次元で、両式は次元二である。$\kappa=0$ では大域位相の自由度が復活し、物理角は木で質量ゼロになる。ベクトルの質量は同じ $v$ を保つ限りこの極限でも変わらない。
:::
::::

::::exercise{#exr-qft-symmetry-breaking-11 level=advanced}
動径場のポテンシャルを線形 $O(3)$ 模型とし、小さい源で角質量を $m_\pi$、動径質量を $m_\sigma=4m_\pi$ にしたとする。木で $\pi_1\pi_1\to\pi_2\pi_2$ の閾値振幅を $m_\pi^2/f^2$ の単位で求め、先頭の二微分の結果と比べよ。

:::solution
この添字では $\delta_{ab}\delta_{cd}=1$、交差した二つの積は零なので $s$ チャネルの係数だけが入る。$s=4m_\pi^2$、$m_\sigma^2=16m_\pi^2$ を明示的破れの振幅へ入れると
$$
\mathcal M=\frac{(4-1)m_\pi^2}{f^2}\frac{16-1}{16-4}
=\frac{15}{4}\frac{m_\pi^2}{f^2}.
$$
二微分の近似は $3m_\pi^2/f^2$ なので、厳密な木の振幅がその $5/4$ 倍。この差はループから来たものではなく、有限の $s/m_\sigma^2$ を動径伝播に残したためである。$m_\sigma/m_\pi\to\infty$ なら比は一に近づく。$s$ が $m_\sigma^2$ に近い場合は、さらに打切りを上げるだけでは共鳴の留数や幅を一様に扱えない。
:::
::::

## 確認テスト

::::quiz
連続対称性を持つ作用で自発的な破れが生じたときの説明を選べ。

- [x] 選んだ真空が対称変換で別の真空に写る
- [ ] 作用の対称性が明示的に除かれる
- [ ] 全ての極小が同じ場値である

:::explanation
作用は変換を保つが、秩序変数の期待値を選んだ真空は同じ変換で別の期待値へ移る。
:::
::::

::::quiz{answer="2" tol="1e-6"}
$O(3)$ 線形 σ 模型で非零真空の角方向の個数を答えよ。

:::explanation
三成分から動径方向を一つ除くので二つである。
:::
::::

::::quiz
高エネルギーの縦偏極で先頭に現れる四元ベクトルを選べ。

- [x] $k^\mu/m_A$
- [ ] $(1,\vb0)$
- [ ] $m_A k^\mu/E^2$

:::explanation
$\epsilon_L^\mu=(p,0,0,E)/m_A$ で、$E-p=m_A^2/(E+p)$ と展開すると $k^\mu/m_A$ が先頭になる。
:::
::::

::::quiz{answer="9" tol="1e-6"}
$\xi=9$ のとき $\chi$ の質量二乗を $m_A^2$ の倍数で答えよ。

:::explanation
ゲージ固定項の $\chi^2$ の係数は $-\xi m_A^2/2$ なので九倍である。
:::
::::

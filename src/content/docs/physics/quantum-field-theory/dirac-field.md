---
title: "Dirac 場の量子化"
description: "Dirac 方程式の自由解を規格化し、反交換関係とスピノル伝播関数を導く。"
notation:
  - symbol: '\not p'
    meaning: ガンマ行列と運動量の縮約
  - symbol: '\bar\psi'
    meaning: Dirac 共役
  - symbol: 'u_s,v_s'
    meaning: 粒子と反粒子のスピノル
---

## この章の目標

- Dirac 方程式と双一次形式の共変性を示す。
- 自由解のスピン和と反交換関係を求める。
- Dirac 場と Weyl 場、Majorana 場の自由度を区別する。

## Dirac 方程式とその Lagrangian

### 一階の時間発展

相対論的エネルギーの二乗は $E^2=\vb p^2+m^2$ である。この式を一次の演算子に因数分解するには、互いに反交換する行列が必要になる。:term[Clifford 代数]{en="Clifford algebra" yomi="くりふぉーどだいすう"} を
$$
\{\gamma^\mu,\gamma^\nu\}=2\eta^{\mu\nu},\quad
(\gamma^0)^\dagger=\gamma^0,\quad(\gamma^i)^\dagger=-\gamma^i
\label{eq-qft-clifford}
$$
で定め、$\not p=\gamma^\mu p_\mu$ と書く。:ref[prop-qft-lorentz-irrep] の二つの Weyl 空間の直和で、$\gamma^\mu$ の非対角ブロックを $\sigma^\mu,\bar\sigma^\mu$ に取るとこの関係を満たす。下添字つきの運動量は $p_\mu=(E,-\vb p)$ である。

$\bar\psi=\psi^\dagger\gamma^0$ として $\mathcal L=\bar\psi(i\not \partial-m)\psi$ を用いる。$\bar\psi$ の変分で $(i\not \partial-m)\psi=0$、$\psi$ の変分と部分積分で $i(\partial_\mu\bar\psi)\gamma^\mu+m\bar\psi=0$ が得られる。両式で発散を計算すると $j^\mu=\bar\psi\gamma^\mu\psi$ は保存し、$j^0=\psi^\dagger\psi$ である。Dirac 演算子に $i\not \partial+m$ を掛けると $-(\Box+m^2)$ になる。自由スピノルの各成分は従って Klein–Gordon 方程式にも従う。行列の反交換関係は :ref[def-aqm-clifford] も参照する。

### Hermite な密度と境界条件

Dirac 密度の微分部分を $i\bar\psi\gamma^\mu\partial_\mu\psi$ から $i[\bar\psi\gamma^\mu\partial_\mu\psi-(\partial_\mu\bar\psi)\gamma^\mu\psi]/2$ に変えると、差は $-i\partial_\mu j^\mu/2$ である。:eqref[eq-qft-clifford] から $(\gamma^\mu)^\dagger=\gamma^0\gamma^\mu\gamma^0$ が従うので、後者は局所的に Hermite である。境界項を固定した作用は同じ場の方程式を与える。箱の空間境界で $h=-i\boldsymbol\alpha\cdot\nabla+\beta m$ の Hermite 性を調べると
$$
\int_V\dd^3x\,(\psi^\dagger h\chi-(h\psi)^\dagger\chi)
=-i\int_{\partial V}\dd S\,\psi^\dagger(\boldsymbol\alpha\cdot\widehat{\vb n})\chi.
$$
周期境界なら対面の項が相殺する。有限領域で別の条件を使うときは、この双一次境界項が許される任意の二つの場で消えることを確認する。

### 連続の式と一粒子の局所速度

自由な Dirac 方程式を $i\partial_t\psi=h\psi$、$h=-i\boldsymbol\alpha\cdot\nabla+\beta m$ と書く。$\alpha^i=\gamma^0\gamma^i$ と $\beta=\gamma^0$ は Hermite。$\partial_t(\psi^\dagger\psi)$ に方程式とその随伴を代入すると質量の二項が相殺し、$-\partial_i(\psi^\dagger\alpha^i\psi)$ が残る。従って規格化した古典複素スピノルの波動関数では $\int\dd^3x\,\psi^\dagger\psi$ が時間に依らない。固定した実単位ベクトル $\widehat{\vb n}$ に対し $(\boldsymbol\alpha\cdot\widehat{\vb n})^2=1$ なので、その行列の固有値は $\pm1$。Cauchy–Schwarz または固有基底での展開から
$$
|\widehat{\vb n}\cdot\vb j|\le j^0,\qquad j^0=\psi^\dagger\psi,\quad j^i=\psi^\dagger\alpha^i\psi .
$$
$\widehat{\vb n}$ を $\vb j$ の方向に選べば $|\vb j|\le j^0$。従ってこの古典的な電流は未来向きの因果的なベクトルである。量子場の $:\bar\psi\gamma^\mu\psi:$ は粒子と反粒子の逆符号の電荷を数えるため、正規順序化した演算子の電荷密度を同じ各点で正の確率とみなさない。

:::example[正規化しない静止波と移動波の電流]
$m=3\,\mathrm{GeV}$、$E=5\,\mathrm{GeV}$ の正エネルギー状態 $u=\sqrt8(1,0,1/2,0)^\mathsf T$ を取る。$j^0=10\,\mathrm{GeV}$、$\alpha^3$ は該当する上下の成分を交換するため $j^3=2\cdot\sqrt8\cdot\sqrt2\,\mathrm{GeV}=8\,\mathrm{GeV}$。従って比は $4/5$ で、分散関係から求めた群速度 $p_z/E$ と一致する。$p_z=0$ の静止解では下成分が零だから空間電流も零となる。ノルムを一へ直してもこの比は変わらない。
:::

### 低速展開で Pauli の運動項を得る

外部電磁ポテンシャルに $D_\mu=\partial_\mu+iq A_\mu$ で結合する。$A^\mu=(A^0,\vb A)$、$\boldsymbol\pi=-i\nabla-q\vb A$ とすると Hamiltonian は $h=\beta m+qA^0+\boldsymbol\alpha\cdot\boldsymbol\pi$。$\psi=e^{-imt}(\varphi,\chi)^\mathsf T$ とし、Dirac 基底の二成分へ分ければ
$$
(i\partial_t-qA^0)\varphi=(\boldsymbol\sigma\cdot\boldsymbol\pi)\chi,\qquad
(2m+i\partial_t-qA^0)\chi=(\boldsymbol\sigma\cdot\boldsymbol\pi)\varphi .
$$
残留エネルギーと $qA^0$ が $m$ より小さく、運動量も $m$ より小さい状態を取る。第二式の逆を $1/(2m)$ の先頭で展開し $\chi=(\boldsymbol\sigma\cdot\boldsymbol\pi)\varphi/(2m)+\cdots$。$\pi_i$ の交換子は $[-i\partial_i,-qA_j]+[-qA_i,-i\partial_j]=iq(\partial_iA_j-\partial_jA_i)$。従って
$$
(\boldsymbol\sigma\cdot\boldsymbol\pi)^2
=\boldsymbol\pi^2+
\tfrac i2\epsilon_{ijk}\sigma_k[\pi_i,\pi_j]
=\boldsymbol\pi^2-q\,\boldsymbol\sigma\cdot\vb B .
$$
こうして
$$
i\partial_t\varphi=\left(qA^0+\frac{\boldsymbol\pi^2}{2m}
-\frac{q}{2m}\boldsymbol\sigma\cdot\vb B\right)\varphi+\cdots
$$
を得る。:term[Pauli 方程式]{en="Pauli equation" yomi="ぱうりほうていしき"} では $\vb S=\boldsymbol\sigma/2$ なので、磁気モーメントを $\boldsymbol\mu=g(q/(2m))\vb S$ と書くとこの次数で $g=2$ である。外場が急に変わると $i\partial_t$ と $\boldsymbol\pi$ の交換子が高次に入るため、低速の先頭式だけで時間依存を任意に扱うことはできない。

:::example[磁場の向きとスピンのエネルギー差]
$B>0$ を $z$ 方向に取り、同じ軌道状態を二つのスピンで比べる。$q>0$ ならスピン項は $\sigma_z=+1$ で $-qB/(2m)$、$-1$ で $+qB/(2m)$。上と下の差の大きさは $|qB|/m$。例えば $m=2\,\mathrm{GeV}$、$qB=0.04\,\mathrm{GeV}^2$ なら差は $0.02\,\mathrm{GeV}$ であり、$qB/m^2=0.01$ は先頭近似の条件と整合する。$q<0$ にすればどちらのスピンが低いかが交換される。
:::

### 自由 Hamiltonian を正負のエネルギーに対角化する

固定した $\vb p$ で $h=\beta m+K$、$K=\boldsymbol\alpha\cdot\vb p$ と書く。$K^\dagger=K$、$K^2=|\vb p|^2$、$\{\beta,K\}=0$。$E=\sqrt{m^2+|\vb p|^2}$ とし、$m>0$ で
$$
U(\vb p)=\frac{E+m+\beta K}{\sqrt{2E(E+m)}}
$$
と置く。$\beta K$ は反 Hermite で二乗が $-|\vb p|^2$。従って分子とその随伴の積は $(E+m)^2+|\vb p|^2=2E(E+m)$ となり、$U$ はユニタリである。$A=E+m$ と略記すれば
$$
(A+\beta K)(m\beta+K)
=(Am+|\vb p|^2)\beta+(A-m)K
=EA\beta+EK
=E\beta(A+\beta K).
$$
従って $UhU^{-1}=E\beta$。これが自由な :term[Foldy–Wouthuysen 変換]{en="Foldy–Wouthuysen transformation" yomi="ふぉるでぃうーずほいせんへんかん"} の運動量表示である。場に作用する Fourier の変換としては $U$ が $\vb p$ に依存するので非局所になる。外部場があると $\boldsymbol\pi$ 同士やポテンシャルと $\boldsymbol\pi$ が可換でなく、同じ式に $\vb p\mapsto\boldsymbol\pi$ と代入しただけでは対角化できない。

:::example[正エネルギーのスピノルを上二成分へ移す]
$m=3,p_z=4,E=5\,\mathrm{GeV}$ とする。$u=\sqrt8(1,0,1/2,0)^\mathsf T$ に $U=(8+4\gamma^3)/\sqrt{80}$ を作用させると、上の第一成分は $\sqrt8(8+2)/\sqrt{80}=\sqrt{10}$、下の第三成分は $\sqrt8(4-4)/\sqrt{80}=0$。通常のノルム二乗はユニタリ性で $10\,\mathrm{GeV}$ のままである。この変換は実ラピディティの Lorentz スピノル変換とは違い、固定した $\vb p$ で正定値内積を保つ基底変更である。
:::

## 自由粒子解とスピン和

### 静止系から任意の運動量へ

$u_s(p)e^{-ipx},v_s(p)e^{ipx}$ に対し $(\not p-m)u_s=0$、$(\not p+m)v_s=0$ である。$m>0$ とし、静止系で $\gamma^0 u=u$、$\gamma^0 v=-v$。固有空間は各二次元である。$\gamma^0=\diag(1,1,-1,-1)$、$\gamma^i=\begin{pmatrix}0&\sigma_i\\-\sigma_i&0\end{pmatrix}$ の基底では
$$
u_s=\sqrt{E+m}\binom{\xi_s}{\frac{\boldsymbol\sigma\cdot\vb p}{E+m}\xi_s},\qquad
v_s=\sqrt{E+m}\binom{\frac{\boldsymbol\sigma\cdot\vb p}{E+m}\eta_s}{\eta_s}
$$
である。$\xi_s,\eta_s$ をそれぞれ正規直交基底に取り、$(\boldsymbol\sigma\cdot\vb p)^2=E^2-m^2$ を使うと $u_s^\dagger u_r=v_s^\dagger v_r=2E\delta_{sr}$、$\bar u_su_r=2m\delta_{sr}$、$\bar v_sv_r=-2m\delta_{sr}$。例えば $E=5,m=3,p_z=4$ を GeV 単位で取ると、$\sigma_3$ の正の固有スピノルの $u$ で下と上の非零成分の比は $1/2$ である。

スピン和は $\sum_s\xi_s\xi_s^\dagger=1$ を各ブロックに代入すれば求まる。例えば $u\bar u$ の左上は $(E+m)1$、右下は $-(E-m)1$、非対角は $-\boldsymbol\sigma\cdot\vb p$ とその負号つきである。従って
$$
\sum_su_s(p)\bar u_s(p)=\not p+m,\qquad
\sum_sv_s(p)\bar v_s(p)=\not p-m .
\label{eq-qft-spin-sums}
$$
$m=0$ でも $E>0$ なら上の式で平方根の分母が消えないので同じスピン和が成り立つ。

### スピンを指定した外部脚

$\gamma^5=i\gamma^0\gamma^1\gamma^2\gamma^3$ とする。$p^2=m^2>0$ 上で四元スピン $s^\mu$ を $p\cdot s=0$、$s^2=-1$ として選ぶ。静止系で $s=(0,\widehat{\vb n})$ なら、正エネルギー部分の $\gamma^5\not s$ は $\boldsymbol\sigma\cdot\widehat{\vb n}$ である。$p\cdot s=0$ と Clifford 関係から $[\gamma^5\not s,\not p]=0$、$(\gamma^5\not s)^2=1$。従って正射影をスピンにも掛ければ
$$
u(p,s)\bar u(p,s)=\tfrac12(1+\gamma^5\not s)(\not p+m)
$$
となる。静止系で右辺の上二成分が $m(1+\boldsymbol\sigma\cdot\widehat{\vb n})$、下二成分が零になるため規格化も一致する。任意の $p$ への移行は共変性で従う。偏極を測定する散乱では :eqref[eq-qft-spin-sums] の和の代わりにこの式を使える。

### Hamiltonian の射影と負周波数

行列 $h(\vb p)=\boldsymbol\alpha\cdot\vb p+\beta m$ は $h^2=E^2$ だから、$\Lambda_\pm(\vb p)=(1\pm h/E)/2$ は直交射影である。$\tr h=0$ より両者とも階数二である。ここで $v(p)e^{ipx}$ の空間波数は $-\vb p$ なので、$h(\vb p)$ の負固有空間を張るのは $v(E,-\vb p)$ である。実際、:eqref[eq-qft-spin-sums] の $u\bar u$ に $\gamma^0$ を右から掛けると $E+h(\vb p)$、$v(E,-\vb p)\bar v(E,-\vb p)$ に同じ操作をすると $E-h(\vb p)$ になる。

:::example[正・負の固有ベクトルの直交]
$m=3$、$p_z=4\,\mathrm{GeV}$ では $u_\uparrow=\sqrt8(1,0,1/2,0)^\mathsf T$、$v_\uparrow(E,-\vb p)=\sqrt8(-1/2,0,1,0)^\mathsf T$ である。各ノルム二乗は $10\,\mathrm{GeV}$、内積は $8(-1/2+1/2)=0$。$\vb p$ を反転しない $v(p)$ と比較するとこの条件は一般に成立しない。
:::

### 共変射影と Hamiltonian 射影の使い分け

正質量殻で $\not p^2=m^2$ だから $\Pi_\pm=(1\pm\not p/m)/2$ は冪等で、積が零で和は一である。例えば $\Pi_+u=u$、$\Pi_-v=v$。任意の連結 Lorentz 変換に対し $\not{\Lambda p}=S\not pS^{-1}$ なので、$\Pi_\pm$ は共変な射影になる。一方 $\not p^\dagger=\gamma^0\not p\gamma^0$ であり、$\vb p\ne0$ では一般に $\Pi_\pm^\dagger=\Pi_\pm$ ではない。従って $\Pi_\pm$ は通常のユークリッド型の四成分内積で直交射影だとみなせない。正定値の内積をもつ Hamiltonian の射影 $\Lambda_\pm(\vb p)$ と比べると、両者は同じスピノルを異なる運動量の割付けで分類する。

具体的に、$p^0=E$ の上で $(\not p-m)u=0$ に左から $\gamma^0$ を掛ければ $(E-h(\vb p))u=0$。$(\not p+m)v=0$ では $h(-\vb p)v=-Ev$。従って $\Pi_-$ の像は $h(\vb p)$ の負固有空間とは $\vb p$ の符号だけ違う。混合直交性もこの区別を反映する。同じ $p$ の $u,v$ で $\bar u v=0$ は、$\bar u\not p=m\bar u$ と $\not p v=-mv$ を比較して $2m\bar u v=0$ から従う。一方 $u^\dagger v$ は一般に零ではなく、$u^\dagger(\vb p)v(-\vb p)=0$ が Hermite な $h(\vb p)$ の直交性になる。

:::example[二種類の直交性を成分で確認する]
$m=3,E=5,p_z=4$ で $u=\sqrt8(1,0,1/2,0)^\mathsf T$、$v(p)=\sqrt8(1/2,0,1,0)^\mathsf T$。$u^\dagger v(p)=8$ だが $\bar uv(p)=8(1/2-1/2)=0$。負枝の波数を合わせた $v(E,-\vb p)=\sqrt8(-1/2,0,1,0)^\mathsf T$ なら $u^\dagger v(E,-\vb p)=0$。従ってスピン和で運動量の符号を変える必要がある箇所は、どの積で直交させるかから判断できる。
:::

### Gordon 恒等式で電流の中心運動とスピンを分ける

正の質量殻上の $p,p'$ に対し、両方の Dirac 方程式を使うと $2m\bar u'\gamma^\mu u=\bar u'(\not p'\gamma^\mu+\gamma^\mu\not p)u$。各積を反交換子と交換子の半分に分けると
$$
\not p'\gamma^\mu+\gamma^\mu\not p
=(p'+p)^\mu+i\sigma^{\mu\nu}(p'-p)_\nu .
$$
$\sigma^{\mu\nu}=i[\gamma^\mu,\gamma^\nu]/2$ の規約を用いた。従って
$$
\bar u'\gamma^\mu u
=\frac1{2m}\bar u'\left[(p'+p)^\mu+i\sigma^{\mu\nu}(p'-p)_\nu\right]u .
\label{eq-qft-gordon}
$$
これが :term[Gordon 恒等式]{en="Gordon identity" yomi="ごーどんこうとうしき"} である。$q_\mu=p'_\mu-p_\mu$ を入れると最初の項は $p'^2-p^2=0$、第二項は $q_\mu q_\nu$ の対称性で零。従って電流保存を二項それぞれで確認できる。この分解には $m>0$ と外脚の殻方程式を使っている。内部線を任意の運動量で置いた行列恒等式としては使えない。

:::example[前方散乱と小さい運動量移行]
$p'=p$ ではスピン項が消え、同じスピンに対して電流は $(p^\mu/m)\bar u(p)u(p)=2p^\mu$。例えば $p=(5,0,0,4)\,\mathrm{GeV}$ なら $(10,0,0,8)\,\mathrm{GeV}$。移行が小さいとき、中心運動の項は $p'+p$、スピン項は $p'-p$ で別のベクトルに比例する。$\vb p'=-\vb p$ の重心に近い配置では空間の中心運動の項が消えるが、スピン項まで零になるとは限らない。両項を分けると磁気結合の抽出に使える。
:::

### ヘリシティ基底でスピン和を確かめる

$\vb p\ne0$ の場合、$\boldsymbol\sigma\cdot\widehat{\vb p}$ の固有値 $\lambda=\pm1$ の固有ベクトル $\xi_\lambda$ を使う。$(\boldsymbol\sigma\cdot\widehat{\vb p})^2=1$ と Hermite 性から、射影は $R_\lambda=(1+\lambda\boldsymbol\sigma\cdot\widehat{\vb p})/2$。各固有値の空間は一次元なので $\xi_\lambda\xi_\lambda^\dagger=R_\lambda$。Dirac 表示で $u_\lambda=(\sqrt{E+m}\,\xi_\lambda,\lambda\sqrt{E-m}\,\xi_\lambda)^\mathsf T$ となる。$|\vb p|/(E+m)=\sqrt{(E-m)/(E+m)}$ を使った。この基底では四成分スピン $\boldsymbol\Sigma=\diag(\boldsymbol\sigma,\boldsymbol\sigma)$ が上下に同じ作用をするから、ヘリシティは $\lambda/2$。さらに $\boldsymbol\Sigma\cdot\widehat{\vb p}$ は $h(\vb p)$ と可換であり、$\boldsymbol\sigma\cdot\vb p$ と同じ軸の二行列が可換であることから分かる。

例えば $\widehat{\vb p}=(\sin\theta\cos\phi,\sin\theta\sin\phi,\cos\theta)$ に対し、$\xi_+=(\cos(\theta/2),e^{i\phi}\sin(\theta/2))^\mathsf T$ を取る。行列 $\begin{pmatrix}\cos\theta&e^{-i\phi}\sin\theta\\e^{i\phi}\sin\theta&-\cos\theta\end{pmatrix}$ を掛けると半角の加法定理で同じ $\xi_+$ になる。直交する $\xi_-=(-e^{-i\phi}\sin(\theta/2),\cos(\theta/2))^\mathsf T$ は固有値 $-1$。$u_+\bar u_+$ と $u_-\bar u_-$ を足すと、対角ブロックでは $R_++R_-=1$、非対角では $R_+-R_-=\boldsymbol\sigma\cdot\widehat{\vb p}$ が出る。従って :eqref[eq-qft-spin-sums] のスピン和は軸の選び方に依らない。

:::example[運動量を $x$ 方向に変えたスピン状態]
$m=3\,\mathrm{GeV}$、$\vb p=(4,0,0)\,\mathrm{GeV}$ なら $\theta=\pi/2,\phi=0$。正ヘリシティの二成分は $\xi_+=(1,1)^\mathsf T/\sqrt2$ であり、Dirac スピノルは $u_+=\sqrt8(1,1,1/2,1/2)^\mathsf T/\sqrt2$。$\Sigma_x u_+=u_+$ である。$\Sigma_z$ の期待値は上と下の二成分とも零なので、$z$ 方向のスピン測定では二つの符号を等しい確率で得る。運動量方向の射影と固定軸での射影は同じ観測量ではない。
:::

## Dirac 行列と双一次形式

Lorentz 変換 $(x,\psi)\mapsto(\Lambda x,S\psi)$ の共変性には $S^{-1}\gamma^\mu S=\Lambda^\mu{}_\nu\gamma^\nu$ を用いる。$[\gamma^\mu,\gamma^\nu]$ と $\gamma^\rho$ の交換子を :eqref[eq-qft-clifford] で整理すると $4(\eta^{\nu\rho}\gamma^\mu-\eta^{\mu\rho}\gamma^\nu)$ である。$\Lambda^\mu{}_\nu=\delta^\mu_\nu+\omega^\mu{}_\nu$、$\omega_{\mu\nu}=-\omega_{\nu\mu}$ と書くと、この恒等式により $S=1+X$ の $X$ を求められる。実際、$X=\omega_{\mu\nu}[\gamma^\mu,\gamma^\nu]/8$ を選べば、
$$
[\gamma^\rho,X]=\tfrac12\omega_{\mu\nu}
 (\eta^{\rho\mu}\gamma^\nu-\eta^{\rho\nu}\gamma^\mu)
=\omega^\rho{}_\nu\gamma^\nu .
$$
従って $S^{-1}\gamma^\rho S=\gamma^\rho+[\gamma^\rho,X]$ により無限小の共変性が従う。有限変換ではこの関係を指数化して、連結群の変換を組み合わせればよい。さらに $S^\dagger\gamma^0=\gamma^0S^{-1}$ が成り立つ。実際、無限小生成子 $X$ は $X^\dagger\gamma^0=-\gamma^0 X$ を満たし、指数級数に適用すると有限変換の関係になる。

$\gamma^5=i\gamma^0\gamma^1\gamma^2\gamma^3$ は各 $\gamma^\mu$ と反交換し、二乗が一である。従って偶数個のガンマ行列からなる Lorentz 生成子と可換である。スピノル二個を縮約する :term[双一次形式]{en="spinor bilinear" yomi="そういちじけいしき"} の例は $\bar\psi\psi$、$i\bar\psi\gamma^5\psi$、$\bar\psi\gamma^\mu\psi$、$\bar\psi\gamma^\mu\gamma^5\psi$、$\bar\psi\sigma^{\mu\nu}\psi$（$\sigma^{\mu\nu}=i[\gamma^\mu,\gamma^\nu]/2$）であり、それぞれスカラー、擬スカラー、ベクトル、軸性ベクトル、反対称テンソルである。連結群では「擬」の有無で同じ変換をし、空間反転で区別される。行列としては $1+1+4+4+6=16$ 個である。Clifford 基底の異なる添字集合の積はトレースが零である。積が奇数階なら含まないガンマ行列、偶数階なら含むガンマ行列で共役して符号が逆になり、トレースの巡回性で従う。各基底要素の二乗は非零定数なのでこの 16 個は独立である。

### Clifford 基底への係数の射影

完全な基底を $1,\gamma^5,\gamma^\mu,\gamma^\mu\gamma^5,\sigma^{\mu\nu}$ とする。トレースの対 $(A,B)\mapsto\tr(AB)$ では $\tr(\gamma^\mu\gamma^\nu)=4\eta^{\mu\nu}$、$\tr(\gamma^\mu\gamma^5\gamma^\nu\gamma^5)=-4\eta^{\mu\nu}$。さらに四つのガンマ行列の反交換で $\tr(\sigma^{\mu\nu}\sigma^{\rho\sigma})=4(\eta^{\mu\rho}\eta^{\nu\sigma}-\eta^{\mu\sigma}\eta^{\nu\rho})$。従って $M$ の展開のベクトル係数は $\tr(\gamma_\mu M)/4$、軸性ベクトル係数は $-\tr(\gamma_\mu\gamma^5 M)/4$ となる。テンソルを全添字で足す場合は $M\supset \sigma^{\mu\nu}\tr(\sigma_{\mu\nu}M)/8$ である。これは同じ反対称成分を二度数えるためである。例えば $M=\gamma^0\gamma^1$ は $-i\sigma^{01}$ だけを持つ。

### スピノルの生成子から有限変換を求める

$S^{-1}\gamma^\mu S=\Lambda^\mu{}_\nu\gamma^\nu$ とする。$z$ 方向のブーストでは $\Lambda^0{}_0=\Lambda^3{}_3=\cosh\chi$、$\Lambda^0{}_3=\Lambda^3{}_0=\sinh\chi$ とする。$\alpha^3=\gamma^0\gamma^3$ は二乗が単位だから
$$
S=\exp(\chi\alpha^3/2)=\cosh(\chi/2)+\alpha^3\sinh(\chi/2).
$$
左に逆、右に $S$ を掛け、$\gamma^0\alpha^3=\gamma^3$、$\alpha^3\gamma^0=-\gamma^3$ を使うと $S^{-1}\gamma^0S=\gamma^0\cosh\chi+\gamma^3\sinh\chi$。$\gamma^3$ にも同様に $\gamma^0$ が混じる。横の二つのガンマ行列は $\alpha^3$ と可換なのでそのまま。これで有限の変換則が直接確かめられる。$S^\dagger=S$ だが $\gamma^0 S^{-1}=S\gamma^0$ なので Dirac 共役との双一次は正しく変換する。

:::example[静止スピノルをブーストで規格化する]
$m=3\,\mathrm{GeV}$、$\cosh\chi=5/3$、$\sinh\chi=4/3$ とする。半角の公式で $\cosh(\chi/2)=2/\sqrt3$、$\sinh(\chi/2)=1/\sqrt3$。静止の $\sqrt6(1,0,0,0)^\mathsf T$ に作用させると $\sqrt2(2,0,1,0)^\mathsf T$ になる。通常のノルム二乗は $6$ から $10\,\mathrm{GeV}$ に増すが、Dirac 双一次は $2(4-1)=6\,\mathrm{GeV}$ を保つ。二つの規格化量が異なる変換をするためである。
:::

### 奇数ガンマのトレースと Levi–Civita テンソル

$\gamma^5\gamma^\mu=-\gamma^\mu\gamma^5$ と $(\gamma^5)^2=1$ を用いる。奇数個のガンマ積 $M$ に対して $\tr M=\tr(\gamma^5 M\gamma^5)=-\tr M$ なので零。四ガンマの積に $\gamma^5$ を一つ入れる場合、同じ添字があれば反交換でその対を隣に移せ、残りは二ガンマの積のトレースになる。このとき $\gamma^5$ は少なくとも一つの使っていないガンマで共役して符号が変わる一方、二つの因子は符号を二つ受ける。従ってこのトレースは零である。四つの添字が全て異なる場合は反交換で順を $(0,1,2,3)$ にそろえ、その置換の符号が付く。$\epsilon^{0123}=+1$ とすれば
$$
\tr(\gamma^5\gamma^\mu\gamma^\nu\gamma^\rho\gamma^\sigma)
=-4i\epsilon^{\mu\nu\rho\sigma}.
\label{eq-qft-gamma5-trace}
$$
実際 $\gamma^0\gamma^1\gamma^2\gamma^3=-i\gamma^5$ だから、順をそろえた積のトレースは $-4i$。この公式は四次元の Clifford 代数に依存する。連続次元のループで同じ反交換と巡回性をすべて無条件に要求すると、軸性電流の正則化で齟齬が起こり得る。

:::example[反対称な縮約に同方向の運動量を入れる]
四つのベクトルを $a,b,c,d$ とすれば上式で $\tr(\gamma^5\not a\not b\not c\not d)=-4i\epsilon^{\mu\nu\rho\sigma}a_\mu b_\nu c_\rho d_\sigma$。$d=a+b$ なら右辺は線形性で同じベクトルを二つ持つ項の和なので零。これは別々の質量殻を仮定しなくても成り立つ。一方、$a=e_0,b=e_1,c=e_2,d=e_3$ と反変の基底を取ると下添字の空間三成分に負号が付き、結果は $+4i$。上添字のガンマ積に対する式の $-4i$ と矛盾しない。
:::

### Clifford 基底の完全性から Fierz の形を得る

$A$ を十六個の Clifford 基底の一つ、$A^\vee$ を $\tr(A^\vee B)=4\delta_{AB}$ の双対とする。例えば $\gamma^\mu$ の双対は $\gamma_\mu$、$\gamma^\mu\gamma^5$ の双対は $-\gamma_\mu\gamma^5$、$\sigma^{\mu\nu}$ の双対は独立な $\mu<\nu$ に対して $\sigma_{\mu\nu}$。任意の四次行列 $M$ の展開は $M=\sum_A A\,\tr(A^\vee M)/4$ である。成分で $\tr(A^\vee M)=(A^\vee)_{lk}M_{kl}$ と書き、$M_{kl}$ の係数を比較して
$$
\delta_{ik}\delta_{jl}=\frac14\sum_A A_{ij}(A^\vee)_{lk}
$$
を得る。添字を適切に交換すれば、四つのスピノルの間で縮約の相手を変えられる。これが :term[Fierz 変換]{en="Fierz rearrangement" yomi="ふぃえるつへんかん"} の基礎である。

例えば通常の数としてのスピノル $u,v,w,z$ について、上の式の $k$ と $l$ の割付けを交換し
$$
(\bar u v)(\bar w z)
=\frac14\sum_A(\bar u A z)(\bar w A^\vee v)
$$
となる。Grassmann 奇な場では左辺の順序 $\bar u,v,\bar w,z$ から右辺の $\bar u,z,\bar w,v$ へ入れ替える際に一回の奇置換があるので、右辺全体に負号を付ける。従ってスピン和を計算する通常のスピノルの恒等式と、局所四フェルミオン演算子を並べ替える式で符号を区別する。同一点の量子の積でこれを使う場合は、共通の正則化と演算子の規約の下で行う。

:::example[単一の行列単位で完全性を検算する]
$M_{kl}$ を $M=e_1e_1^\mathsf T$ とし、Dirac 基底で $e_1=(1,0,0,0)^\mathsf T$。双対トレースでは $\tr M=1$、$\tr(\gamma^0M)=1$、$\tr(\gamma^3\gamma^5 M)=1$、$\tr(\sigma^{12}M)=1$ が非零である。独立成分の展開では $M=(1+\gamma^0+\gamma^3\gamma^5+\sigma^{12})/4$。実際後二つは対角で $\gamma^3\gamma^5=\diag(1,-1,-1,1)$、$\sigma^{12}=\diag(1,-1,1,-1)$ なので和の対角は $(4,0,0,0)$。軸性部分では $\gamma_3=-\gamma^3$ と双対の全体の負号が打ち消し合う。
:::

### 双一次形式の Hermite 性と整数階の共変量

$\bar\psi\Gamma\psi$ の随伴は $\bar\psi(\gamma^0\Gamma^\dagger\gamma^0)\psi$ である。随伴を取ると二つの場の並びは反転するが、もともとの $\bar\psi=\psi^\dagger\gamma^0$ を使えば同じ $\psi^\dagger(\cdots)\psi$ の順序に戻るので、ここに独立なフェルミの交換符号を加えない。従って $\Gamma$ が Dirac 随伴で自身と等しいときに、この双一次は Hermite になる。例えば $\gamma^5$ は通常の随伴で Hermite だが $\gamma^0\gamma^5\gamma^0=-\gamma^5$ なので $\bar\psi\gamma^5\psi$ は反 Hermite。$i\gamma^5$ を用いると、$i$ の共役と $\gamma^5$ の負号が相殺して Hermite になる。$\gamma^\mu\gamma^5$ では通常の随伴で順序が反転するが、Dirac 随伴で $\gamma^5$ に付く負号をもう一度反交換で打ち消す。従って軸性ベクトルも Hermite である。

例として自由スピノルのスピン射影 $u\bar u=(1+\gamma^5\not s)(\not p+m)/2$ を使う。$\tr(\gamma^\mu\gamma^5u\bar u)$ を取り、$\gamma^5$ のない項や奇数ガンマの項を除くと
$$
\bar u(p,s)\gamma^\mu\gamma^5u(p,s)
=\frac m2\tr(\gamma^\mu\not s)=2m s^\mu .
$$
従って運動量を測るベクトル電流 $2p^\mu$ と異なり、軸性ベクトルにはスピン四元ベクトルが現れる。例えば静止系で $s=(0,0,0,1)$、$m=3\,\mathrm{GeV}$ なら軸性電流は $(0,0,0,6)\,\mathrm{GeV}$。$z$ 方向へ $p_z=4\,\mathrm{GeV}$ までブーストすれば時間成分が $8\,\mathrm{GeV}$、空間成分が $10\,\mathrm{GeV}$ に変わる。二乗は $8^2-10^2=-36\,\mathrm{GeV}^2$ で、静止系の長さを保つ。

## Dirac 場の量子化と反交換関係

### 正エネルギーと反粒子

$\dd\Pi_p=\dd^3p/[(2\pi)^3 2E_{\vb p}]$ はスカラー場と同じ測度である。場を
$$
\psi(x)=\sum_s\int\dd\Pi_p\bigl[b_s(p)u_s(p)e^{-ipx}+d_s^\dagger(p)v_s(p)e^{ipx}\bigr]
\label{eq-qft-dirac-expansion}
$$
と置く。:term[正準反交換関係]{en="canonical anticommutation relations" yomi="せいじゅんはんこうかんかんけい"} は $\{b_r(\vb p),b_s^\dagger(\vb q)\}=\{d_r(\vb p),d_s^\dagger(\vb q)\}=2E(2\pi)^3\delta_{rs}\delta^3(\vb p-\vb q)$、他の組は零である。:eqref[eq-qft-spin-sums] に右から $\gamma^0$ を掛け、負周波数項の $\vb p$ を反転して足すと $2E$ 倍の単位行列が残る。従って等時刻で $\{\psi_\alpha(\vb x),\psi_\beta^\dagger(\vb y)\}=\delta_{\alpha\beta}\delta^3(\vb x-\vb y)$ となる。

一粒子 Hamiltonian $h=-i\gamma^0\gamma^i\partial_i+m\gamma^0$ は $u e^{-ipx}$ に $+E$、$v e^{ipx}$ に $-E$ を与える。場の $H=\int\psi^\dagger h\psi\,\dd^3x$ を計算すると交差項は固有空間の直交性で消え、箱で $H=\sum E(b^\dagger b-d d^\dagger)$ となる。$dd^\dagger=1-d^\dagger d$ なので真空定数を引いた $H$ は $\sum E(b^\dagger b+d^\dagger d)$ である。電荷 $Q=\int:\psi^\dagger\psi:\,\dd^3x$ は同じ計算で $\sum(b^\dagger b-d^\dagger d)$ になる。$b^\dagger$ が作る粒子の電荷を $+1$ と数えると $d^\dagger$ の反粒子は $-1$ である。

### 正準運動量と一次の作用

運動項 $i\psi^\dagger\dot\psi$ では、Grassmann 奇な $\dot\psi$ を右端に置いたときの係数は $i\psi^\dagger$ であり、$\dot\psi^\dagger$ の係数は零である。従って二つの速度を共役量から独立に解くことはできない。代わりに一次作用の係数が定める対を使う。直交モードに制限して一つの複素成分を $c$ と書けば $L=i c^*\dot c-Ec^*c$。$c^*$ の変分で $i\dot c=Ec$、$c$ の変分で $-i\dot c^*=Ec^*$ が出る。量子化した $H=Ec^\dagger c$ と $\{c,c^\dagger\}=1$ からは $[c,H]=Ec$ となり同じ時間発展を再現する。実際 $cc^\dagger c=(1-c^\dagger c)c=c$ であり、$(c)^2=0$ を使った。

:::example[二つのモードと交換の符号]
正規化した $c_1^\dagger c_2^\dagger\ket0$ に $c_1$ を作用させると $c_2^\dagger\ket0$、$c_2$ を作用させると $-c_1^\dagger\ket0$ になる。$E_1=3\,\mathrm{GeV}$、$E_2=5\,\mathrm{GeV}$ なら初状態のエネルギーは $8\,\mathrm{GeV}$ である。消滅後に残る位相の符号とエネルギーの符号は別の量である。
:::

### モード展開から反交換関係を逆に取り出す

規格化した箱内で、等時刻の場からモードを取り出す。$\{u_s(E,\vb p),v_s(E,-\vb p)\}$ を合わせた $h(\vb p)$ の正規直交基底は、各ベクトルを $\sqrt{2E}$ で割れば得られる。従って
$$
\begin{aligned}
b_s(\vb p)&=\frac{1}{\sqrt{2EV}}\int_V\dd^3x\,
                         u_s(p)^\dagger e^{ipx}\psi(x),\\
d^\dagger_s(\vb p)&=\frac{1}{\sqrt{2EV}}\int_V\dd^3x\,
                         v_s(p)^\dagger e^{-ipx}\psi(x)
\end{aligned}
$$
となる。ここで箱の反交換子は波数とスピンについて Kronecker のデルタで、規格化係数を連続の測度のものから変えた。第一式に展開を戻すと粒子側は $u^\dagger u=2E$ と空間積分 $V$ が前置係数の積を相殺する。反粒子側では運動量の符号が反転するので $u^\dagger(\vb p)v(-\vb p)=0$。第二式も同様である。

逆に等時刻の $\{\psi_\alpha(\vb x),\psi_\beta^\dagger(\vb y)\}$ をデルタと仮定し、このモード積分へ代入する。空間デルタで二つの積分の点が一致し、その後の Fourier 積分で $\vb p=\vb q$。スピノルの内積から $\delta_{sr}$ が出るので $\{b_s(\vb p),b_r^\dagger(\vb q)\}=\delta_{sr}\delta_{\vb p\vb q}$。$\{d_s,d_r^\dagger\}$ も同じ。混合項に同じ空間指数が出る組では $u,v$ が直交し、残りの組は $\{\psi,\psi\}=0$ またはその随伴を使う。従って等時刻の局所な反交換関係とモードの反交換関係は、この完全性の下で同値である。

:::example[二つの異なる波数で同じスピンを選ぶ]
一辺 $L=2\pi\,\mathrm{GeV}^{-1}$ の箱で $p_z=1\,\mathrm{GeV}$ と $q_z=2\,\mathrm{GeV}$ を取る。時間部分やスピノル成分に非零の内積があっても、$\int_0^L \exp[i(q_z-p_z)z]\dd z=0$ なので $b(p)$ と $b^\dagger(q)$ の反交換子は零である。逆に波数が同じでスピンが異なる場合には空間積分は $V$ だが、$\xi_\uparrow^\dagger\xi_\downarrow=0$ により同じ結論になる。デルタは二種類の直交性から出る。
:::

### 反対称 Fock 空間と波束のノルム

一粒子波束の空間を $\mathcal H_1$ とする。$n$ 粒子には置換の符号を付けた射影 $A_n=\sum_{\sigma\in\mathfrak S_n}\sgn(\sigma)U_\sigma/n!$ を使う。随伴で置換を逆にしても符号は同じだから $A_n^\dagger=A_n$。二乗では各積 $\sigma\tau$ に対して係数 $\sgn(\sigma)\sgn(\tau)=\sgn(\sigma\tau)$ が共通で、同じ積が $n!$ 回出るので $A_n^2=A_n$。従って反対称部分への直交射影である。生成写像は $c^\dagger(f)\psi_n=\sqrt{n+1}\,A_{n+1}(f\otimes\psi_n)$、消滅写像はその随伴になる。先頭の波束を巡回させると交換回数に応じて符号が付くので、二つの生成演算子は互いに反交換する。

モードの反交換関係から、波束の生成演算子を $c^\dagger(f)=\sum_s\int\dd\Pi_p f_s(p)c_s^\dagger(p)$ とすれば $\{c(f),c^\dagger(g)\}=\inner f g$。二粒子のノルム二乗は $\|f\|^2\|g\|^2-|\inner f g|^2$。実際、$c(f)$ を二つの生成因子の右へ進めると直接項は $\|f\|^2 c^\dagger(g)\ket0$、交換項は $-\inner f g\,c^\dagger(f)\ket0$。さらに $c(g)$ を作用させると式を得る。一般の $n$ でも最初の因子を消す順序の符号が行列式の余因子展開と一致するので、ノルム二乗は波束の Gram 行列式になる。線形従属な波束しかない場合にはノルムが零で状態が消える。

:::example[重なる波束を二つ入れた場合]
規格化した直交モード $f,h$ から $g=(f+\sqrt3 h)/2$ と取る。$c^\dagger(f)c^\dagger(g)\ket0=(\sqrt3/2)c^\dagger(f)c^\dagger(h)\ket0$ であり、同じ $f$ の二重占有は零。従ってノルム二乗は $3/4$ で、Gram の式 $1-1/4$ と一致する。状態を正規化した後に占有数を測ると、$f$ と $h$ に一個ずつある。ボースの場合に同じ係数の $f,g$ を使うと、$f$ の二重占有は消えず、規格化も変わる。
:::

### 運動量を測る演算子と電荷を測る演算子

正規順序化した空間運動量は $P^i=\int\dd^3x:\psi^\dagger(-i\partial^i_{\rm sp})\psi:$ と書ける。ここで $\partial^i_{\rm sp}$ は通常の空間座標に関する微分 $\partial/\partial x^i$ であり、計量で上げた四次元微分ではない。展開すると $b^\dagger b$ の空間指数は $e^{i\vb p\cdot\vb x}$ なので $+p^i$、$dd^\dagger$ の指数は $e^{-i\vb p\cdot\vb x}$ なので $-p^i$。後者を正規順序にするとさらに負号が付き、結果は両種とも $+p^i$ になる。時間成分も同じように両種へ $E$ を掛ける。他方、電荷密度には空間微分がないため、正規順序の負号が相殺されず $Q=N_b-N_d$ が残る。粒子と反粒子の運動量はどちらも同じ未来殻に置けるが、電荷は逆である。

:::example[反粒子を二粒子状態に一つ加える]
箱で $m=3\,\mathrm{GeV}$、粒子が $\vb p=(0,0,4)\,\mathrm{GeV}$ と $\vb q=\vb0$ の二モード、反粒子が $\vb p'=(0,0,-4)\,\mathrm{GeV}$ の一モードを占めるとする。正規順序の全電荷は $2-1=1$。全運動量の空間部分は零、全エネルギーは $5+3+5=13\,\mathrm{GeV}$。もし反粒子の運動量を負周波数の指数の向きだけで数えれば $8\,\mathrm{GeV}$ の空間成分を誤って得る。物理粒子の演算子 $d^\dagger$ の変換を正規順序の $P$ から読む必要がある。
:::

## スピンと統計の関係

反交換関係により同一モードの $(b^\dagger)^2=0$ であり、多粒子空間は反対称テンソル積の直和である。自由 Dirac 場で普通の正ノルムのボース交換関係を使うと $-Edd^\dagger=-Ed^\dagger d-E$ なので、反粒子を増やすほどエネルギーが下がる。逆に反粒子の交換子の符号を変えてこの問題を避けようとすると、一粒子ノルムが負になる。

さらに自由場の反交換子は
$$
\{\psi_\alpha(x),\bar\psi_\beta(y)\}
=[(i\not \partial_x+m)C(x-y)]_{\alpha\beta}
$$
である。$C$ は :eqref[eq-qft-scalar-commutator] である。微分は分布の台を広げないので空間的分離でこの反交換子が消える。偶数個の場からなる局所可観測量同士では、因子を移す符号が偶数回現れて交換子が消える。:term[スピンと統計の関係]{en="spin–statistics connection" yomi="すぴんととうけいのかんけい"} の一般形では、Poincaré 共変でスペクトルが閉じた未来光円錐にある理論を、正定値 Hilbert 空間と巡回真空の上で考える。場とその共役が空間的分離で一様な交換符号を持つとき、整数スピンには交換、半整数スピンには反交換が対応する。逆の符号を課すと、その場は真空から励起を作れない。一般の証明は Wightman 分布の複素領域での解析接続を要するため :cite[weinberg-qft1] に譲る。

### 偶の観測量の局所性

空間的に分離した $x,y$ で、二つの双一次形式の積 $\bar\psi(x)\Gamma\psi(x)\bar\psi(y)M\psi(y)$ を交換する。$y$ の最初の場は $x$ の二因子を越えるので符号は $(-1)^2$ であり、$y$ の次の場でも同じである。従って全体では $+1$。この計算では相互反交換子が零であることを用いた。例えば電荷密度を空間的に分離した二領域に塗りならせば交換子は零である。一方、同一点の積には正則化が要るので、この議論で接触項を否定することはできない。

### 交換関係を誤ると空間的分離で何が残るか

正の一粒子ノルムを保ったまま、仮に Dirac 展開で反交換関係を普通のボース交換関係へ置き換える。$\psi$ の $b$ と $\bar\psi$ の $b^\dagger$ の交換子は正、$\psi$ の $d^\dagger$ と $\bar\psi$ の $d$ は負である。スピン和を使うと
$$
[\psi(x),\bar\psi(0)]_{\rm B}
=\int\dd\Pi_p\big[(\not p+m)e^{-ipx}-(\not p-m)e^{ipx}\big]
=(i\not\partial+m)\bigl(D^+(x)+D^+(-x)\bigr).
$$
フェルミの反交換子では右の括弧が差になるため、空間的な分離で消える。一方、ここで得た和は等時刻の $r>0$ でも一般に非零。$m>0$ ならスピノル添字のトレースは $4m[D^+(x)+D^+(-x)]$ であり、微分の寄与は $\tr\gamma^\mu=0$ で消える。等時刻相関 $D^+(0,\vb x)>0$ だからトレースが非零となる。従ってエネルギーの不安定性に加え、ボース型の局所交換も成り立たない。

逆向きの練習として、複素スカラーを粒子と反粒子の正ノルムを保って反交換させると、$\{\varphi(x),\varphi^\dagger(y)\}=D^+(x-y)+D^+(y-x)$。やはり空間的分離で非零となる。例えば四次元の質量ゼロで等時刻、距離 $2\,\mathrm{GeV}^{-1}$ ならこの和は $1/(8\pi^2)\,\mathrm{GeV}^2$。整数スピンと半整数スピンで必要な符号の違いは、自由モードの正と負周波数を組み合わせたときにも直接見える。この自由場での反例は、一般の相互作用理論の定理の証明とは別の計算である。

## Dirac 伝播関数

時間順序ではフェルミオンを交換するたびに負号を付ける。したがって $t>0$ の $u\bar u$ の項と、$t<0$ の $-v\bar v$ の項は合わせて
$$
S_F(x)=\mel{0}{T\psi(x)\bar\psi(0)}{0}
=(i\not \partial+m)\Delta_F(x)
=\int\frac{\dd^4p}{(2\pi)^4}\frac{i(\not p+m)e^{-ipx}}{p^2-m^2+i0}
\label{eq-qft-dirac-feynman}
$$
となる。微分を $\theta(t)$ に作用させた二つの境界項は $\Delta_F$ の左右極限が同じため相殺する。$(i\not \partial-m)S_F=i\delta^4(x)$ は :eqref[eq-qft-scalar-feynman] と演算子の因数分解から従う。同じ式を運動量空間で確かめると、分子積は $(p^2-m^2)1$ である。

### 等時刻の跳びと Green 方程式

正の時間側の極限は $\int\dd\Pi_p(\not p+m)e^{i\vb p\cdot\vb x}$。負の時間側の極限は $-\int\dd\Pi_p(\not p-m)e^{-i\vb p\cdot\vb x}$。負の時間側で $\vb p\mapsto-\vb p$ として差を取ると $2E\gamma^0$ だけが残る。従って $S_F(0^+,\vb x)-S_F(0^-,\vb x)=\gamma^0\delta^3(\vb x)$。一次演算子の時間微分を階段関数にも作用させると、$i\gamma^0\gamma^0\delta^4(x)=i\delta^4(x)$ が得られる。例えば $\vb p=0$ に制限した正の時間部分は $e^{-imt}(1+\gamma^0)/2$ であり、負の時間部分は $-e^{imt}(1-\gamma^0)\gamma^0/2$ である。

### 時間処方を替えたときの Dirac 核

遅延 Klein–Gordon 核 $G_R$ が $(\Box+m^2)G_R=\delta^4$ を満たすとする。一次の遅延逆核を $S_R=-(i\not\partial+m)G_R$ と置くと
$(i\not\partial-m)S_R=(\Box+m^2)G_R=\delta^4$。微分分布の台は元の分布の台に含まれるので、$S_R$ も未来の光円錐外で零である。運動量空間では
$$
S_R(p)=\frac{\not p+m}{(p^0+i0)^2-E_{\vb p}^2}.
$$
分母の両極が下半面にあり、$t<0$ の輪郭がどちらも含まない。:eqref[eq-qft-dirac-feynman] は負エネルギー極が上半面なので、同じ分子でもこの性質は持たない。ここでは $S_R$ を $i$ なしの逆核として定義しており、時間順序積 $S_F$ の逆方程式の右辺には $i$ が付く。

:::example[静止モードへの瞬時入力]
$\vb p=0$ の一次方程式 $(i\gamma^0\partial_t-m)\psi=\eta_0\delta(t)$ で、$t<0$ の場を零とする。時刻零を横切って積分すると $i\gamma^0\psi(0^+)=\eta_0$、従って $\psi(t)=-i\theta(t)e^{-im\gamma^0t}\gamma^0\eta_0$。$\eta_0$ が $\gamma^0$ の正固有成分だけなら位相は $e^{-imt}$、負固有成分なら $e^{imt}$。遅延性は二つの周波数を排除する条件ではなく、ともに源の入った後に応答させる境界条件である。
:::

### 残留行列から粒子と反粒子の射影を取り出す

:eqref[eq-qft-dirac-feynman] の運動量表示で、$p^0$ の正極 $E$ の留数は $i(\not p_++m)/(2E)$、$p_+=(E,\vb p)$。負極 $-E$ の留数は $-i(\not p_-+m)/(2E)$、$p_-=(-E,\vb p)$。後者で $q=-p_-=(E,-\vb p)$ と置けば $-i(-\not q+m)/(2E)=i(\not q-m)/(2E)$。正・負極の分子はそれぞれ $u\bar u,v\bar v$ の和になり、$t<0$ の積分路の向きが時間順序のフェルミ符号を作る。例えば $\vb p=0$ なら正極の分子は $m(1+\gamma^0)$、負極側の $v\bar v$ は $m(\gamma^0-1)$。負極の残留行列の符号だけを見て反粒子ノルムが負だと判断してはいけない。状態の正定値ノルムでは Dirac 共役に含まれる $\gamma^0$ をもう一度掛けて $v v^\dagger$ を使う。

:::example[時間の Fourier 積分を Hamiltonian の分解で求める]
右から $\gamma^0$ を掛けた $S_F(t,\vb p)$ を考える。$t>0$ では $u u^\dagger/(2E)=\Lambda_+(\vb p)$ の和に $e^{-iEt}$ が掛かる。$t<0$ では $-v(E,-\vb p)v(E,-\vb p)^\dagger/(2E)=-\Lambda_-(\vb p)$ の和に $e^{iEt}$ が掛かる。従って
$$
S_F(t,\vb p)\gamma^0=\theta(t)e^{-iEt}\Lambda_+
                         -\theta(-t)e^{iEt}\Lambda_- .
$$
差の時刻零の極限は $\Lambda_++\Lambda_-=1$。一つの符号の時間指数だけを使うとこの完全性を失い、デルタの右辺が全四成分に作用する一次方程式の逆にならない。
:::

### 質量ゼロの Dirac 伝播関数の特異性

空間的な分離 $x^2<0$ では質量ゼロのスカラー核は $\Delta_F(x)=1/[4\pi^2(-x^2)]$。この領域で $S_F=i\not\partial\Delta_F$ を取ると
$$
S_F(x)=\frac{i\not x}{2\pi^2(x^2)^2},\qquad x^2<0 .
$$
微分では $\partial_\mu x^2=2x_\mu$ を使い、$\not x=\gamma^\mu x_\mu$ とした。ノルムの尺度を変え $x\mapsto\lambda x$ とすると $\lambda^{-3}$。四次元で二つの Dirac 場の正準次元の和が $3$ であることと一致する。光円錐上ではこの通常関数の式だけで分布を定められないので、:eqref[eq-qft-dirac-feynman] の $i0$ と導関数の定義を使う。この例でも時系列相関は空間的分離で一般に非零であり、局所性を説明するのはその値でなく反交換子である。

:::example[同時刻の一軸で尺度の冪を確認する]
$x=(0,0,0,r)$、$r>0$ なら $\not x=-r\gamma^3$、$S_F=-i\gamma^3/(2\pi^2r^3)$。$r=1\,\mathrm{GeV}^{-1}$ から $2\,\mathrm{GeV}^{-1}$ へ移すと係数は $1/8$。同じ二点でスカラー核なら比は $1/4$ である。Dirac 伝播関数の余分な逆長さは、分子の一次微分から出る。
:::

## Weyl 場と Majorana 場

### カイラル射影と共役条件

$P_{L,R}=(1\mp\gamma^5)/2$ は直交する冪等射影で、$\gamma^\mu P_L=P_R\gamma^\mu$。質量ゼロでは Dirac 方程式を射影すると各カイラル成分だけの一次方程式に分離する。Weyl 基底で $\gamma^5=\diag(-1,-1,1,1)$ なら、左成分は $i\bar\sigma^\mu\partial_\mu\xi=0$ を満たす。正周波数解で $(E+\boldsymbol\sigma\cdot\vb p)\xi=0$ なので粒子のヘリシティは $-1/2$。質量項は $\bar\psi_L\psi_R+\bar\psi_R\psi_L$ となって両成分を結ぶ。$\bar\psi_L=\bar\psi P_R$ であることに注意する。

荷電共役を $\psi^c=C\bar\psi^{\mathsf T}$、$C\gamma^{\mu\mathsf T}C^{-1}=-\gamma^\mu$ で定義する。上で用いた Dirac 基底では $C=i\gamma^2\gamma^0$ と取れるので $\psi^c=i\gamma^2\psi^*$、$(\psi^c)^c=\psi$ である。:term[Majorana 場]{en="Majorana field" yomi="まよらなば"} は $\psi^c=\psi$ を満たす場である。この条件は粒子と反粒子のモードを同一視するので独立な生成演算子の組が半分になる。連続位相回転で $\psi^c$ は逆の位相を受けるから、条件を保つ連続的な $U(1)$ 電荷は零である。Majorana の作用には Dirac 型の密度に $1/2$ を付け、共役条件で二重に数える成分を除く。

### 質量があるときのカイラリティの混合

Weyl 基底で $u=(u_L,u_R)^{\mathsf T}$ とし、$\boldsymbol\sigma\cdot\widehat{\vb p}$ の固有値を $\lambda=\pm1$ とする。殻方程式は $(E+\lambda p)u_L=mu_R$、$(E-\lambda p)u_R=mu_L$ である。$u^\dagger u=2E$ で規格化すると $\|u_L\|^2=E-\lambda p$、$\|u_R\|^2=E+\lambda p$ になる。例えば $\lambda=+1$ なら左成分が占める割合は $(E-p)/(2E)=m^2/[2E(E+p)]$。高エネルギーで抑制されるが、静止に近い場合は左右が同程度である。従ってカイラリティは場の Lorentz 表現の区別であり、ヘリシティは指定した運動量方向に対する状態のスピン射影である。

### 軸性密度とカイラル射影

$\bar\psi_L=\bar\psi P_R$、$P_R\gamma^\mu=\gamma^\mu P_L$ なので $j_L^\mu=\bar\psi\gamma^\mu P_L\psi$、$j_R^\mu=\bar\psi\gamma^\mu P_R\psi$。和が $j^\mu$、差が $j_5^\mu=\bar\psi\gamma^\mu\gamma^5\psi$ である。自由場の方程式を両因子に作用させると
$$
\partial_\mu j_5^\mu
 =im\bar\psi\gamma^5\psi+im\bar\psi\gamma^5\psi .
$$
第二項に負号が付かないのは、$\gamma^\mu$ と $\gamma^5$ を入れ替えるためである。従って質量ゼロでは自由な左右の電流が別々に保存する。量子の相互作用でこの結論を使う場合は [アノマリー](/physics/quantum-field-theory/anomalies/) の正則化を確認する。

### Weyl 質量と Grassmann の反対称性

左 Weyl スピノルを $\xi=(\xi_1,\xi_2)^\mathsf T$ とする。$\varepsilon_{12}=1=-\varepsilon_{21}$ と選ぶと $\xi^\mathsf T\varepsilon\xi=\xi_1\xi_2-\xi_2\xi_1=2\xi_1\xi_2$ である。Grassmann 奇なので、この収縮は同じ場でも零にならない。$A^\mathsf T\varepsilon A=\varepsilon$ より Lorentz 不変である。従って $-(m\xi^\mathsf T\varepsilon\xi+\text{Hermite 共役})/2$ は中性の Weyl 場に質量を与えられる。電荷 $q$ の場であればこの双一次は位相 $e^{-2iq\alpha}$ を受け、非零の $q$ では対称性を保たない。逆符号の電荷をもつ二つの Weyl 場なら、同じ収縮で電荷零の Dirac 質量を作れる。

### 二成分で Majorana 条件を解く

Weyl 基底を $\gamma^\mu=\begin{pmatrix}0&\sigma^\mu\\\bar\sigma^\mu&0\end{pmatrix}$、$\sigma^\mu=(1,\boldsymbol\sigma)$、$\bar\sigma^\mu=(1,-\boldsymbol\sigma)$ と定める。$\gamma^5=\diag(-1_2,1_2)$ で、$\varepsilon=i\sigma^2$ は実で $\varepsilon^2=-1_2$。この基底の荷電共役を
$$
\Psi=\binom{\xi}{\eta}\quad\mapsto\quad
\Psi^c=i\gamma^2\Psi^*
=\binom{\varepsilon\eta^*}{-\varepsilon\xi^*}
$$
と取る。二度作用させると上は $-\varepsilon^2\xi=\xi$、下も同様で元に戻る。また左の射影をしてから共役すると結果は右の空間に入る。Majorana 条件 $\Psi^c=\Psi$ を解くと $\eta=-\varepsilon\xi^*$、上の条件は $\varepsilon\eta^*=-\varepsilon^2\xi=\xi$ で自動的に成り立つ。従って独立な複素二成分は $\xi$ 一つである。Dirac 場は両成分を独立に持つので、この条件を付けない。

自由 Majorana 場では質量殻にある同じ演算子の生成と消滅を正・負周波数に使う。$u$ と $v=i\gamma^2u^*$ を一組とすれば
$$
\Psi_M(x)=\sum_s\int\dd\Pi_p\,[c_s(p)u_s(p)e^{-ipx}+c_s^\dagger(p)v_s(p)e^{ipx}]
$$
となる。共役を取ると二項が入れ替わるので条件を満たす。独立な $c$ だけを数えるため、同じ質量を持つ Dirac 場の $b,d$ の二組と比べて粒子種の数は半分である。例えば静止系にスピン状態が二つあっても、反粒子に別の二つを追加する必要はない。$\langle T\Psi_M(x)\bar\Psi_M(y)\rangle$ の二つの時間枝では、通常と同じ $\{c,c^\dagger\}$ とスピン和が出るので :eqref[eq-qft-dirac-feynman] と同じ核になる。独立な自由度の数が半分という条件を、伝播関数全体を半分にする操作と混同しない。

### 中性場二つと荷電した Dirac 場一つの関係

同じ質量の二つの独立な Majorana 場を $\Psi_1,\Psi_2$ とする。$\Psi_D=(\Psi_1+i\Psi_2)/\sqrt2$ を作ると $\Psi_D^c=(\Psi_1-i\Psi_2)/\sqrt2$。従って $\Psi_D$ と $\Psi_D^c$ の和差から二つの Majorana 場へ戻れる。二つの実な場の間の $SO(2)$ 回転を $(\Psi_1,\Psi_2)\mapsto(\cos\alpha\,\Psi_1-\sin\alpha\,\Psi_2,\sin\alpha\,\Psi_1+\cos\alpha\,\Psi_2)$ と取ると、$\Psi_D$ の作用は $e^{i\alpha}\Psi_D$。同じ質量では自由作用がこの回転を保つため、Dirac 表現で連続な電荷になる。質量が異なる場合、質量行列はこの回転に可換でなくなり、同じ構成をしても電荷は保存しない。

:::example[Majorana のベクトル電流が零になること]
同一点の Grassmann 奇な場の局所積を用いる。荷電共役した Dirac ベクトル電流は $\overline{\Psi^c}\gamma^\mu\Psi^c=-\bar\Psi\gamma^\mu\Psi$ である。実際 $\overline{\Psi^c}=-\Psi^\mathsf T C^{-1}$ を代入して二場を交換し、$(C^{-1}\gamma^\mu C)^\mathsf T=-\gamma^\mu$ を使う。$\Psi=\Psi^c$ なら電流が自身の負に等しいので零。二つの独立な Majorana 場を使った $\Psi_D$ の電流では、対角の $\bar\Psi_i\gamma^\mu\Psi_i$ がそれぞれ零でも、交差項が残る。$\bar\Psi_1\gamma^\mu\Psi_2=-\bar\Psi_2\gamma^\mu\Psi_1$ なので、$\bar\Psi_D\gamma^\mu\Psi_D=i\bar\Psi_1\gamma^\mu\Psi_2$。一場のベクトル電流が零という結論は、二場の遷移電流までは消さない。
:::

### 複数の Weyl 場の質量を固有値で読む

左 Weyl 場 $\xi_i$ の質量項を $-\tfrac12(M_{ij}\xi_i^\mathsf T\varepsilon\xi_j+\text{Hermite 共役})$ とする。$\xi_i^\mathsf T\varepsilon\xi_j=\xi_j^\mathsf T\varepsilon\xi_i$ であり、二成分の反対称性と奇変数の交換の負号が相殺する。従って $M$ の対称部分だけが作用する。ここでは実対称の $M$ に限り、実直交行列で対角化する。この変換は運動項を保つ。固有値が負なら該当する Weyl 場を $i$ 倍する再定義でその質量項の符号が逆になるため、物理質量は固有値の絶対値である。

:::example[Dirac の質量と対角の Majorana 質量を比較する]
$M=\begin{pmatrix}0&m_D\\m_D&M_R\end{pmatrix}$、$m_D,M_R>0$ とする。固有値は $[M_R\pm\sqrt{M_R^2+4m_D^2}]/2$。$M_R=0$ なら $\pm m_D$ で絶対値が等しい。対角の質量項がない元の作用で二つの Weyl 場に正反対の位相を掛けても質量項が変わらないので、一つの Dirac 粒子として扱える。$M_R>0$ なら二つの絶対値が異なり、その位相対称性は破れる。例えば $M_R=3\,\mathrm{GeV}$、$m_D=2\,\mathrm{GeV}$ なら固有値は $4,-1\,\mathrm{GeV}$、従って質量は $4,1\,\mathrm{GeV}$。$M_R\gg m_D$ では軽い方は $m_D^2/M_R+\order{m_D^4/M_R^3}$。これは二次方程式の平方根を小さい比で展開して求める。
:::

### カイラル射影を伝播関数の両端に掛ける

$\psi_L=P_L\psi$ だが $\bar\psi_L=\bar\psi P_R$ であることを使う。同じ左成分の二点 $\langle T\psi_L(x)\bar\psi_L(y)\rangle$ の運動量核は $P_L S_F(p)P_R$。$\not p$ は $\gamma^5$ と反交換して左右を交換するため $P_L\not pP_R=\not pP_R$、一方 $mP_LP_R=0$。従って
$$
P_L S_F(p)P_R=\frac{i\not pP_R}{p^2-m^2+i0}.
$$
これに対し $\bar\psi_R=\bar\psi P_L$ を右に置く相関では、$P_L\not pP_L=0$、$P_LmP_L=mP_L$ なので
$$
P_L S_F(p)P_L=\frac{imP_L}{p^2-m^2+i0}.
$$
前者は質量ゼロでも残り、左 Weyl 場単独の伝播に対応する。後者は左右の成分を結ぶが、質量を零にすると消える。従って質量ゼロでのカイラル電流の独立性は、波動方程式からも内部線の分子からも確認できる。

:::example[二成分の核の逆を直接求める]
左 Weyl 場の自由逆核は $\bar\sigma^\mu p_\mu$。$(\bar\sigma^\mu p_\mu)(\sigma^\nu p_\nu)=p^2 1_2$ は Pauli 関係で従うので、Feynman の二成分核は $i\sigma^\nu p_\nu/(p^2+i0)$ である。例えば $p=(E,0,0,k)$ で $\bar\sigma^\mu p_\mu=E+k\sigma_3$。分子 $E-k\sigma_3$ を掛けて $E^2-k^2$ を得る。正の光的 $E=k>0$ に近づくと分子は $\diag(0,2E)$ であり、極の留数は下の一成分だけに入る。これは $z$ 方向の左粒子のヘリシティが $-1/2$ であることと同じである。
:::

### ヘリシティとカイラリティの混合を状態の観測値で測る

有質量のヘリシティ $\lambda=\pm1$ の正エネルギー $u_\lambda$ について、Weyl 成分のノルムは $\|u_R\|^2=E+\lambda p$、$\|u_L\|^2=E-\lambda p$。従って規格化した $\gamma^5$ の期待値は
$$
\frac{u_\lambda^\dagger\gamma^5u_\lambda}{u_\lambda^\dagger u_\lambda}
=\lambda \frac pE .
$$
これは $\lambda/2$ というヘリシティ固有値とは異なる。例えば $m=3\,\mathrm{GeV}$、$p=4\,\mathrm{GeV}$、$\lambda=+1$ なら右成分の割合は $9/10$、左は $1/10$、期待値は $4/5$。静止極限では方向のヘリシティの定義自体が退化するが、任意の固定方向からその極限へ近づけると左右の割合はともに $1/2$。高エネルギー極限では $\lambda=+1$ の左成分の割合が $m^2/(4E^2)$ の程度になる。

この関係は外部脚にカイラルな頂点を掛けるときに使える。一例として、軸に沿う状態の左成分への射影そのもののノルムは $\sqrt{(E-p)/(2E)}\simeq m/(2E)$。しかし振幅の大きさを決めるには、他の脚や Lorentz 指標との縮約も必要である。射影したノルムの比を、あらゆる散乱過程の断面積の比とみなさない。

:::example[荷電共役が二成分のスピン射影を変えること]
$\xi=(1,0)^\mathsf T$ は $\sigma_3$ の正固有ベクトルである。一方 $-\varepsilon\xi^*=(0,1)^\mathsf T$ には負の固有値が対応する。一般に $\varepsilon\sigma_i^*\varepsilon^{-1}=-\sigma_i$ だから、任意の実方向 $\widehat{\vb n}$ に対して $\xi$ の固有値を $\lambda$ とすると $-\varepsilon\xi^*$ の固有値は $-\lambda$ である。これは二成分の数学的なスピン射影の比較である。粒子と反粒子の物理状態のヘリシティを読む際は、負周波数の波数と、生成演算子に移したスピン指標の変換を合わせて扱う。例えば左 Weyl 場の二成分に含まれる正周波数は左粒子を消し、負周波数は右ヘリシティの反粒子を作る。場の成分がともに左 Weyl 表現に属することと、二つの物理状態が逆のヘリシティを持つことは両立する。
:::

### 複素対称な Weyl 質量行列の正準化

複素の場合も、正準な運動項を持つ左 Weyl 場の質量行列 $M$ は $M^{\mathsf T}=M$ としてよい。一般の相似対角化では質量項に入る二つの場に逆行列を掛けるわけではないから、固有値の絶対値をそのまま使えない。必要なのはユニタリ行列 $V$ が与える合同変換 $V^{\mathsf T} M V$ である。

有限次元の複素空間で反線形作用 $B v=M\bar v$ を考える。$B^2=M\overline M=MM^\dagger$ は Hermite で半正定値である。$B$ はこの二乗と可換で、$B^2$ の固有値 $m_i^2\ge0$ は実なので各固有空間を保つ。$m_i>0$ の空間では $B/m_i$ の二乗は一で、内積は
$$
\langle Bv,Bw\rangle
=\overline{\langle v,MM^\dagger w\rangle}
=m_i^2\langle w,v\rangle
$$
を満たす。第一式は $M^\dagger M=\overline{MM^\dagger}$ を使った。従って $B/m_i$ は反ユニタリな対合である。任意のベクトルを $v=(v+C v)/2+i(v-C v)/(2i)$ と書くと、$C=B/m_i$ の不動ベクトル二つの複素線形結合になる。不動空間の内積は実なので、実 Gram–Schmidt 法で得た正規直交な実基底が元の複素空間の基底にもなる。$m_i=0$ の空間では $\|Bv\|^2=\langle v,MM^\dagger v\rangle=0$ より $B=0$。これで全部の基底をまとめると $M\overline U=U\diag(m_i)$、$U$ はユニタリである。従って $V=\overline U$ を用いれば $V^{\mathsf T} M V=\diag(m_i)$。この :term[Takagi 分解]{en="Takagi factorization" yomi="たかぎぶんかい"} により、物理質量は $MM^\dagger$ の固有値の非負平方根になる。

:::example[複素質量の通常の固有値と物理質量]
対称行列 $M=\begin{pmatrix}1&i\\i&0\end{pmatrix}$ を共通のエネルギー単位で取る。通常の固有値は $(1\pm i\sqrt3)/2$ で絶対値はともに一。一方 $MM^\dagger=\begin{pmatrix}2&-i\\ i&1\end{pmatrix}$ の固有値は $(3\pm\sqrt5)/2$ なので、二つの物理質量は $(\sqrt5\pm1)/2$ である。運動項を保つユニタリな合同変換で求めた量は、通常の相似変換で得た固有値と異なる。
:::

## 演習問題

::::exercise{#exr-qft-dirac-field-1 level=basic}
$(\not p)^2$ を求めよ。

:::solution
$p_\mu p_\nu$ が対称なのでガンマ積を反交換子の半分に置ける。よって $p_\mu p_\nu\eta^{\mu\nu}=p^2$ である。
:::
::::

::::exercise{#exr-qft-dirac-field-2 level=basic}
$m=3,p_z=4,E=5$ を GeV 単位とする。$u$ の上下成分の比を $\sigma_3=+1$ で求めよ。

:::solution
式の比は $p_z/(E+m)=4/8=1/2$ である。
:::
::::

::::exercise{#exr-qft-dirac-field-3 level=basic}
$P_LP_R$ と $P_L+P_R$ を求めよ。

:::solution
$(\gamma^5)^2=1$ から積は $(1-1)/4=0$、和は $1$ である。
:::
::::

::::exercise{#exr-qft-dirac-field-4 level=basic}
スピンが指定された一反粒子モードに粒子が一個ある。$dd^\dagger$ の固有値を求めよ。

:::solution
$dd^\dagger=1-d^\dagger d$ なので占有数一では零である。
:::
::::

::::exercise{#exr-qft-dirac-field-5 level=advanced}
$\bar u(\not p-m)=0$ を示せ。

:::solution
$(\not p-m)u=0$ の随伴を取る。$(\gamma^\mu)^\dagger=\gamma^0\gamma^\mu\gamma^0$ を代入して右に $\gamma^0$ を掛けると $\bar u(\not p-m)=0$ になる。
:::
::::

::::exercise{#exr-qft-dirac-field-6 level=advanced}
$\not p+m$ の正質量殻でのランクを求めよ。

:::solution
静止系で $m(\gamma^0+1)$ は二つの固有値 $2m$ と二つの零を持つ。任意の未来質量殻上の $p$ は静止系からの Lorentz 変換で得られ、対応する行列は相似だからランクは二である。
:::
::::

::::exercise{#exr-qft-dirac-field-7 level=advanced}
Dirac 伝播関数が逆演算子になることを運動量空間で検算せよ。

:::solution
$(\not p-m)i(\not p+m)=i(p^2-m^2)$ なので極の処方を保って分母で割ると分布の極限は $i$ である。逆 Fourier 変換では $i\delta^4(x)$ になる。
:::
::::

::::exercise{#exr-qft-dirac-field-8 level=advanced}
Majorana 条件のもとで $\psi\mapsto e^{i\alpha}\psi$ が許される位相を求めよ。

:::solution
共役は $e^{-i\alpha}\psi^c$ である。変換後も両者が等しい条件は $e^{2i\alpha}=1$ なので、位相は $\pm1$ だけである。
:::
::::

::::exercise{#exr-qft-dirac-field-9 level=advanced}
自由粒子の $h(\vb p)$ を用い、正エネルギー部分空間への射影の中では $\boldsymbol\alpha$ の各成分が $p^i/E$ 倍の恒等作用になることを示せ。$m=5\,\mathrm{GeV}$、$\vb p=(0,0,12)\,\mathrm{GeV}$ では電流の時間成分に対する $z$ 成分の比はいくらか。

:::solution
$\{\alpha^i,h\}=2p^i1$ は $\{\alpha^i,\alpha^j\}=2\delta^{ij}$ と $\{\alpha^i,\beta\}=0$ から得られる。左右に同じ射影 $\Lambda_+=(1+h/E)/2$ を掛けると、$h\Lambda_+=\Lambda_+h=E\Lambda_+$ により左辺は $2E\Lambda_+\alpha^i\Lambda_+$。従って $\Lambda_+\alpha^i\Lambda_+=(p^i/E)\Lambda_+$ である。正エネルギースピノルで $j^i=u^\dagger\alpha^iu$、$j^0=u^\dagger u$ なので比は $p^i/E$。指定の運動量では $E=\sqrt{5^2+12^2}\,\mathrm{GeV}=13\,\mathrm{GeV}$ だから比は $12/13$ である。射影を外した $\alpha^i$ の固有値は $\pm1$ だが、この計算はその固有値ではなく、正エネルギーに制限した電流を評価している。
:::
::::

::::exercise{#exr-qft-dirac-field-10 level=advanced}
加速していない粒子の正エネルギーと負エネルギーの波束の重ね合わせを、固定した $\vb p$ の四成分の計算で考える。$h(\vb p)$ の固有値を $\pm E$ とし、時刻零に両成分 $\psi_+,\psi_-$ があるとき、演算子 $\alpha^i$ の期待値にどの周波数が現れるか。純粋な正エネルギーならその成分が消えることを示せ。

:::solution
時間発展は $\psi(t)=e^{-iEt}\psi_++e^{iEt}\psi_-$。$\alpha^i$ は Hermite なので
$$
\psi(t)^\dagger\alpha^i\psi(t)
=\psi_+^\dagger\alpha^i\psi_+
+\psi_-^\dagger\alpha^i\psi_-
+2\Re\bigl[e^{2iEt}\psi_+^\dagger\alpha^i\psi_-\bigr].
$$
最後の項に角周波数 $2E$ がある。純粋な正エネルギーでは $\psi_-=0$ だから消える。固定 $\vb p$ で正負成分の通常の内積は零だが、間に $\alpha^i$ を挿入すれば零とは限らない。例えば静止 $m=3\,\mathrm{GeV}$ で $\psi_+=e_1/\sqrt2$、$\psi_-=e_3/\sqrt2$ とし、$i=3$ を選ぶ。$\alpha^3$ はこの二成分を交換するので、交差の係数は $1/2$、対角項は零、期待値は $\cos(6\,\mathrm{GeV}\,t)$。一粒子では正エネルギーに射影する一方、場の演算子では負周波数側を反粒子の生成として実装する。
:::
::::

## 確認テスト

::::quiz
Dirac 場を反交換関係で量子化した場合、反粒子の励起のエネルギーはどれか。

- [x] $+E$
- [ ] $-E$
- [ ] $0$

:::explanation
負周波数の寄与 $-Edd^\dagger$ を並べ直すと $+Ed^\dagger d-E$ になる。定数を引いた後は励起に $+E$ が伴う。
:::
::::

::::quiz{answer="4" tol="1e-6"}
$m=2\,\mathrm{GeV}$ の静止スピノルの $\bar u u$ を GeV 単位で答えよ。

:::explanation
静止系で $u$ は $\gamma^0=+1$ の固有ベクトルで、規格化から $\bar uu=2m=4\,\mathrm{GeV}$ である。
:::
::::

::::quiz
左 Weyl 場の正エネルギー粒子のヘリシティはどれか。

- [x] $-1/2$
- [ ] $0$
- [ ] $+1/2$

:::explanation
$(E+\boldsymbol\sigma\cdot\vb p)\xi=0$ と $E=|\vb p|$ から、$\boldsymbol\sigma\cdot\widehat{\vb p}$ の固有値は $-1$、スピン射影はその半分になる。
:::
::::

::::quiz{answer="16" tol="1e-6"}
四成分スピノル上の複素行列空間の次元を答えよ。

:::explanation
四行四列に独立な成分が一つずつあるので次元は $4^2=16$ である。
:::
::::

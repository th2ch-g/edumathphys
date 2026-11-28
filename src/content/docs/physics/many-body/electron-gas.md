---
title: "電子ガス"
description: "一様な電子ガスの交換、誘電応答、遮蔽、集団励起と高密度相関を計算する。"
notation:
  - symbol: 'g_C'
    meaning: Coulomb 結合 e の二乗を 4π ε で割った量
  - symbol: 'k_F'
    meaning: Fermi 波数
  - symbol: 'r_s'
    meaning: Wigner–Seitz 半径を有効 Bohr 半径で割った量
  - symbol: '\chi_0,\varepsilon_{\rm RPA}'
    meaning: 自由密度応答と RPA 誘電関数
  - symbol: 'W(q,z)'
    meaning: 遮蔽した Coulomb 相互作用
---

## この章の目標

- 一様背景と交換孔がエネルギーに与える効果を求める。
- Lindhard 関数の小波数極限と特異点を導く。
- 遮蔽とプラズモンを同じ誘電関数から計算する。
- リング和から相関エネルギーの高密度対数を得る。

## Jellium 模型

### 一様背景による中性化

[バンド計算](/physics/solid-state/band-calculation/) ではイオンの周期ポテンシャルを含める。一様な電子ガスではイオンを平均密度だけで表し、電子間相互作用を調べる。三次元の体積 $V$ の箱にスピン二成分の電子を置き、同じ密度の固定正電荷を入れる。これを :term[Jellium 模型]{en="jellium model" yomi="じぇりうむもけい"} と呼ぶ。媒質の誘電率を $\epsilon$、$g_C=e^2/(4\pi\epsilon)$、$v_q=4\pi g_C/q^2$ とする。逆 Laplace 演算子の Fourier 核に $4\pi g_C$ を掛ければ $g_C/r$ が得られる。零波数の電荷は電子と背景が打ち消し、電子の同一粒子項を引くと
$$
H=\sum_{\vb k\sigma}\frac{k^2}{2m}n_{\vb k\sigma}
+\frac1{2V}\sum_{\vb q\ne0} v_q(n_{\vb q}n_{-\vb q}-\hat N).
\label{eq-mb-jellium}
$$
背景を入れてから熱力学極限を取る。$n=k_F^3/(3\pi^2)$ は Fermi 球の体積を数えた結果である。有効 Bohr 半径は $a_B=\hbar^2/(m g_C)$、長さの比は $r_s=(3/(4\pi n))^{1/3}/a_B$ であり、$k_Fa_B=(9\pi/4)^{1/3}/r_s$ になる。

体積 $1/n$ の球の半径 $r_s a_B$ を :term[Wigner–Seitz 半径]{en="Wigner–Seitz radius" yomi="うぃぐなーさいつはんけい"} と呼ぶ。以下の交換と応答は、特に断らない限り三次元、零温度、非偏極の一様な相について計算する。例えば交換エネルギーを求める際には一様な Slater 行列式を使い、密度やスピンが空間変調する試行状態のエネルギーとは別に扱う。

:::example[密度とエネルギーの桁]
真空誘電率と電子質量で $a_B\simeq0.0529\,\mathrm{nm}$ である。$n=10^{29}\,\mathrm{m^{-3}}$ なら $k_F\simeq14.4\,\mathrm{nm^{-1}}$、$r_s\simeq2.5$、$\hbar^2 k_F^2/(2m)\simeq7.9\,\mathrm{eV}$ である。高密度漸近式の小さい変数は $r_s$ である。
:::

:::example[密度を八倍にする]
$k_F\propto n^{1/3}$ なので波数は二倍になる。$r_s\propto n^{-1/3}$ は半分になり、Fermi エネルギーは四倍になる。同じ質量と誘電率の比較であり、波数が二倍なら背景による中性条件はそのまま新しい電子密度に合わせる。
:::

### 密度を変えると二つのエネルギー尺度が分かれる

三次元で $4\pi(r_s a_B)^3n/3=1$ と定義し、$a_B=\hbar^2/(mg_C)$ とする。非偏極の Fermi 球の半径は、許された波数を状態あたり $(2\pi)^3/V$ で数えて $k_F=(3\pi^2n)^{1/3}$ である。従って $k_Fa_B=(9\pi/4)^{1/3}/r_s$ であり、運動エネルギー $E_F=\hbar^2 k_F^2/(2m)$ は $r_s^{-2}$、平均間隔での Coulomb エネルギー $g_C/(r_sa_B)$ は $r_s^{-1}$ に比例する。両者の比は $2r_s/(9\pi/4)^{2/3}$ である。中性化によって巨視的な一様場を除いても、この比は各電子の運動と局所的な相互作用を比較する尺度として残る。

:::example[密度八倍のときの尺度]
$n$ を八倍にすると $r_s$ は半分になる。従って $E_F$ は四倍、平均間隔の Coulomb エネルギーは二倍になる。エネルギーの絶対値は双方とも増えるが、相互作用を運動エネルギーで割った比は半分になる。
:::

### 電荷の零モードと電子の自己項

一辺 $L$ の周期箱で $\vb q=2\pi\vb l/L$、$\vb l\in\Z^3$ とし、電子の数密度の Fourier 成分を $\hat n_{\vb q}=\sum_{\vb k\sigma}c^\dagger_{\vb k\sigma}c_{\vb k+\vb q,\sigma}$ とする。これは第一量子化の $\sum_j e^{-i\vb q\cdot\vb r_j}$ に対応する。電子密度と背景密度の差の零成分は、背景の粒子数を $N$ に合わせた $N$ 粒子部門では零である。Coulomb 核の周期的な定義は、まず $\ell>0$ で $v_{q,\ell}=4\pi g_C/(q^2+\ell^2)$ とし、正負の電荷を合わせてから $\ell\downarrow0$ とする。零成分が $\ell^{-2}$ に比例するので、電子項だけに先にこの極限を適用すると有限のエネルギーを得られない。

位置表示の積 $\hat n_{\vb q}\hat n_{-\vb q}$ には、同じ電子を二度選んだ位相積 $e^{-i\vb q\cdot\vb r_j}e^{i\vb q\cdot\vb r_j}=1$ が電子ごとに含まれる。従って電子対の Fourier 和は $\hat n_{\vb q}\hat n_{-\vb q}-\hat N$ である。演算子でも $c_{\vb k+\vb q}c_{\vb p}^\dagger=\delta_{\vb k+\vb q,\vb p}-c_{\vb p}^\dagger c_{\vb k+\vb q}$ を使えば、デルタ項が $\hat N$ を作り、残りは二体の正規順序積となる。背景は固定した連続密度なので、同じ粒子を選ぶ自己項を持たない。周期 Coulomb 核の空間平均を零に定めれば、この二体の項の $\vb q\ne0$ だけを足した :eqref[eq-mb-jellium] になる。点電荷の自己項の減算は、例えば上限のある Fourier 和の中で先に行う。

:::example[一個の粒子の Coulomb エネルギー]
$N=1$ の部門では任意の $\vb q\ne0$ に対して $\hat n_{\vb q}\hat n_{-\vb q}=1$ である。従って零モードを除き自己項を引いた本章の Hamiltonian では、Coulomb の二体項は零となる。二電子なら自己項を引いた積は $2\cos[\vb q\cdot(\vb r_1-\vb r_2)]$ である。これに係数 $1/(2V)$ を掛けると同一の対を一度数えた値になる。
:::

### 密度を変える座標変換と基底エネルギー

一粒子あたりの体積が $4\pi(r_s a_B)^3/3$ になるように箱の辺を変え、粒子数 $N$ を固定する。長さを $r_s a_B$ で割った変数 $\vb x_j=\vb r_j/(r_s a_B)$ を使うと、無次元の箱の辺 $L/(r_s a_B)=(4\pi N/3)^{1/3}$ は一定である。$\grad_{\vb r_j}=(r_s a_B)^{-1}\grad_{\vb x_j}$ と $g_C/(r_sa_B)=E_{\rm H}/r_s$ から、Hamiltonian を
$$
H=E_{\rm H}\left(\frac{\mathcal T}{r_s^2}+\frac{\mathcal V}{r_s}\right)
$$
と書ける。$\mathcal T=-\tfrac12\sum_j\nabla_{\vb x_j}^2$、$\mathcal V$ は無次元箱の中性化した Coulomb 項である。背景密度も同じ座標変換を受けるため、電子と背景の項と背景どうしの項にも同じ $1/r_s$ が掛かる。この表示では $H/(E_{\rm H}/r_s^2)=\mathcal T+r_s\mathcal V$ となるので、高密度で弱くなるのは運動項に対する相互作用の比である。

基底エネルギーを粒子一個あたり $e(r_s)$ とし、微分可能な基底枝を取る。基底の固有方程式を $r_s$ で微分して基底の左ベクトルを掛けると、基底ベクトルの微分を含む項が相殺する。従って
$$
r_s e'(r_s)=-2e_{\rm kin}-e_{\rm pot},\qquad
e_{\rm kin}=-r_s e'(r_s)-e,\qquad
e_{\rm pot}=2e+r_s e'(r_s)
$$
となる。ここで $e=e_{\rm kin}+e_{\rm pot}$ である。同じ量は、体積を変えると $r_s\propto V^{1/3}$ となることから圧力 $P=-nr_s e'(r_s)/3$ と結び付く。背景も電子とともに中性のまま圧縮しているので、この圧力は固定背景に電子だけを一様に加えたときの応答とは異なる。

:::example[圧力と平均運動エネルギーの照合]
自由な Fermi 海では $e_{\rm kin}=3E_F/5$、$e_{\rm pot}=0$ である。$E_F\propto r_s^{-2}$ なので $r_se'=-2e$、従って $P=2nE_F/5$ となる。例えば $E_F=10\,\mathrm{eV}$、$n=10^{29}\,\mathrm{m^{-3}}$ なら、この自由電子の圧力は $4\times10^{29}\,\mathrm{eV\,m^{-3}}\simeq6.4\times10^{10}\,\mathrm{Pa}$ である。
:::

## Hartree–Fock 近似と交換エネルギー

### Fermi 球の交換積分

一様な Slater 行列式では背景で直接項が相殺され、:eqref[eq-mb-hf-density] の交換項だけが残る。相対変数を $\vb q=\vb k-\vb p$ に変えると二つの球が交わる体積 $B(q)$ を積分できる。$0\le q\le2k_F$ では、底から高さ $k_F-q/2$ の球冠が二つある。円板面積を積分すると $B(q)=\pi(4k_F+q)(2k_F-q)^2/12$ である。スピンの和は $2$ であり、
$$
\frac{E_x}{V}=-\frac{2g_C}{8\pi^4}\int_0^{2k_F}B(q)\dd q
=-\frac{g_Ck_F^4}{4\pi^3},\quad
\frac{E_x}{N}=-\frac{3g_Ck_F}{4\pi}.
\label{eq-mb-jellium-exchange}
$$
運動エネルギーは $k^2$ を球内で平均するため $3E_F/5$ である。Hartree 単位 $E_{\rm H}=g_C/a_B$ では、運動項が $r_s^{-2}$、交換項が $r_s^{-1}$ に比例する。

:::example[交換の比率]
$k_Fa_B=2$ なら $E_{\rm kin}/N=6E_{\rm H}/5$、$E_x/N=-3E_{\rm H}/(2\pi)$ である。スピンが偏極すれば各スピンで別の $k_F$ を使う。全ての電子を一つのスピンに入れ、全密度を保つと $k_F$ が $2^{1/3}$ 倍になるので、交換の絶対値は同じ倍率になる。
:::

:::example[二球の中間の重なり]
$q=k_F$ の重なり体積は $5\pi k_F^3/12$ であり、一つの球の体積 $4\pi k_F^3/3$ に対し $5/16$ である。$q=0$ で重なりが全体積に戻ることも同じ球冠の式から分かる。
:::

### 偏極を入れた運動と交換の競合

全密度を一定に保ち、$\zeta=(n_\uparrow-n_\downarrow)/n$ とする。各スピンで $n_\sigma=k_\sigma^3/(6\pi^2)$ なので $k_\pm=k_F(1\pm\zeta)^{1/3}$ となる。運動エネルギーの密度は $k_\sigma^5$ に比例する。一方、交換の二つの運動量積分は $k_\sigma^6$、Coulomb 頂点は $k_\sigma^{-2}$ に比例し、異スピン交換はない。従って一粒子あたりでは
$$
\begin{aligned}
e_{\rm kin}(\zeta)&=\frac{3E_F}{5}\frac{(1+\zeta)^{5/3}+(1-\zeta)^{5/3}}2,\\
e_{\rm x}(\zeta)&=e_{\rm x}(0)\frac{(1+\zeta)^{4/3}+(1-\zeta)^{4/3}}2 .
\end{aligned}
$$
$\zeta^2$ の係数はそれぞれ $E_F/3$ と $2e_{\rm x}(0)/9$ で、後者は負である。完全偏極と非偏極のエネルギーが交差する条件には $\zeta^4$ 以上の項も入るので、二次の曲率が零になる点と一致するとは限らない。またこの比較は Slater 行列式での交換近似に限る。

:::example[完全偏極での交換の変化率]
$\zeta=1$ とすると交換の大きさは非偏極の $2^{1/3}$ 倍、運動エネルギーは $2^{2/3}$ 倍になる。例えば $3E_F/5=|e_{\rm x}(0)|$ なら完全偏極の正味の差は $(2^{2/3}-2^{1/3})|e_{\rm x}(0)|>0$ である。
:::

### 交換孔を位置空間で積分する

非偏極の自由 Fermi 海で、一スピンの一体密度行列の相対位置 $\vb r$ 成分を $\gamma_\sigma(r)$ とする。球の角積分を取ると
$$
\begin{aligned}
\gamma_\sigma(r)
&=\int_{k<k_F}\frac{\dd^3k}{(2\pi)^3}e^{i\vb k\cdot\vb r}
=\frac{1}{2\pi^2r}\int_0^{k_F}k\sin(kr)\dd k\\
&=\frac{n}{2}f(k_Fr),\qquad
f(x)=\frac{3(\sin x-x\cos x)}{x^3}.
\end{aligned}
$$
$x\to0$ の Taylor 展開で $f(0)=1$ である。:eqref[eq-mb-hf-density] から、原点にスピン $\sigma$ の電子を見つけた条件の下で、残りの電子の数密度は $n+h_x(r)$、$h_x(r)=-|\gamma_\sigma(r)|^2/(n/2)=-(n/2)f(k_Fr)^2$ になる。異なるスピンはこの行列式では独立なので、欠損は同じスピンの密度にだけ現れる。この負の密度を :term[交換孔]{en="exchange hole" yomi="こうかんこう"} と呼ぶ。

規格化は軌道射影の二乗が同じ射影であることから従う。実際、$\int|\gamma_\sigma(r)|^2\dd^3r=\int_{k<k_F}\dd^3k/(2\pi)^3=n/2$ なので $\int h_x(r)\dd^3r=-1$ である。原点に電子を一つ置いた条件では、その電子自身を除くため、残りの密度を全空間で積分した数が $N-1$ になる。背景による直接項の相殺後の交換エネルギーは $E_x/N=\tfrac12\int (g_C/r)h_x(r)\dd^3r$ である。Fourier の積で $|\gamma_\sigma|^2$ を表せば二つの占有球が交わる体積 $B(q)$ を得るので、:eqref[eq-mb-jellium-exchange] と一致する。

:::example[交換孔の原点値と積分値]
原点で同スピンの条件付き密度は $(n/2)[1-f(0)^2]=0$ で、異スピンの条件付き密度は $n/2$ である。全スピンで平均した対分布を $g_0(r)=1+h_x(r)/n$ とすると、$g_0(0)=1/2$、$n\int[g_0(r)-1]\dd^3r=-1$ になる。独立な二粒子を仮定して $g_0=1$ と置くと、原点での同スピンの排除と、電子一個分の規格化の両方が失われる。
:::

### Fock 自己エネルギーの運動量依存

スピンごとの Fock 自己エネルギーは占有した同スピン球で
$$
\Sigma_x(k)=-\int_{p<k_F}\frac{\dd^3p}{(2\pi)^3}\frac{4\pi g_C}{|\vb k-\vb p|^2}
=-\frac{g_C}{\pi k}\int_0^{k_F}p\log\left|\frac{k+p}{k-p}\right|\dd p
$$
となる。角積分の分母は $k^2+p^2-2kpu$ なので、$u$ の積分は $\log|(k+p)/(k-p)|/(kp)$ である。半径積分では、$p^2/2$ を積分して部分積分し、内部に $p=k$ があるときは左右を同じ距離で切ってから極限を取る。得られる式は
$$
\Sigma_x(k)=-\frac{g_Ck_F}{\pi}
\left[1+\frac{1-x^2}{2x}\log\left|\frac{1+x}{1-x}\right|\right],
\qquad x=\frac{k}{k_F}.
\label{eq-mb-jellium-fock}
$$
ここで $-\frac{g_Ck_F}{\pi}$ は角括弧全体に掛かる。$k=0$ の極限は $-2g_Ck_F/\pi$、$k=k_F$ では $-g_Ck_F/\pi$ である。$x=1+t$ の近くでは $\Sigma_x(k)=-g_Ck_F/\pi+(g_Ck_F/\pi)t\log(2/|t|)+O(t^2\log|t|)$ であり、傾きは正に対数発散する。従って裸の Coulomb 頂点の Hartree–Fock 分散から $\rho(E)=k^2/(\pi^2|\partial_k E|)$ を作ると Fermi 面で状態密度が零へ向かう。この特異性は裸の小転送相互作用を残したことに由来し、後で導く動的な遮蔽を含む分散をこの式の微分だけで評価することはできない。

:::example[化学ポテンシャルを二通りで求める]
交換エネルギー密度は $ne_x=-3ng_Ck_F/(4\pi)$ である。$k_F\propto n^{1/3}$ を用いると $\partial_n(ne_x)=-g_Ck_F/\pi=\Sigma_x(k_F)$ となる。また $\partial_n(n\,3E_F/5)=E_F$ である。中性を保って粒子を加える Hartree–Fock 化学ポテンシャルは $\mu_{\rm HF}=E_F-g_Ck_F/\pi$ となり、最後に占有した Fock 軌道のエネルギーと一致する。$k_Fa_B=2$ なら $\mu_{\rm HF}/E_{\rm H}=2-2/\pi$ である。
:::

## Lindhard 関数

### 自由密度応答

:eqref[eq-mb-free-density] の三次元零温度形を :term[Lindhard 関数]{en="Lindhard function" yomi="りんとはるとかんすう"} と呼ぶ。$D_F=mk_F/\pi^2$ とする。分子の第二項で運動量をずらし、$\vb q$ を極軸に取ると静的な積分は
$$
\chi_0(q,0)=-\frac{m}{\pi^2q}\int_0^{k_F} k\log\left|\frac{q+2k}{q-2k}\right|\dd k
=-D_F\left[\frac12+\frac{1-x^2}{4x}\log\left|\frac{1+x}{1-x}\right|\right],\quad x=\frac q{2k_F}.
\label{eq-mb-lindhard-static}
$$
最初の等号では $\int_{-1}^1\dd u/(a+bu)=b^{-1}\log|(a+b)/(a-b)|$ の主値を用いた。次は $y=2k/q$ と置き、原始関数 $y+(y^2-1)\log|(1+y)/(1-y)|/2$ を使った。$x\to0$ の角括弧は $1$、$x=1$ では $1/2$ であり、後者の近くでは微分に対数特異点がある。

小波数で $s=z/(v_Fq)$ を保つと、占有差は $v_Fq u\,\delta(\xi_k)$ に近づき、分母は $v_Fq(s-u)$ である。従って
$$
\chi_0(q,z)\simeq D_F\int_{-1}^1\frac{u}{2(s-u)}\dd u
=D_F\left[\frac s2\log\frac{s+1}{s-1}-1\right].
\label{eq-mb-lindhard-dynamic}
$$
対数の枝は上半平面から取り、$0<s<1$ の実軸では虚部が $-\pi D_Fs/2$ になる。これは粒子と正孔へエネルギーを渡せる領域である。虚周波数 $z=i v_Fq u$、$u>0$ では $\chi_0=-D_F[1-u\arctan(1/u)]$ である。

:::example[小波数での二つの順序]
$s=0$ では $\chi_0=-D_F$ であり、$|s|\gg1$ では対数の展開により $\chi_0=D_F/(3s^2)+O(s^{-4})$ になる。後者の極限では $D_Fv_F^2/3=n/m$ なので $\chi_0\simeq nq^2/(mz^2)$ である。
:::

### 有限波数での粒子–正孔連続体

放物線分散では $a=q^2/(2m)$、$b=v_Fq$ と置き、分子の第二項で運動量をずらすと占有球上の二項の差を取れる。$k_z=k_Fu$ を固定した断面積は $\pi k_F^2(1-u^2)$ であり、スピン和を含めて
$$
\chi_0(q,z)=\frac{k_F^3}{4\pi^2 b}\left[
F\left(\frac{z-a}b\right)-F\left(\frac{z+a}b\right)\right],\qquad
F(s)=2s+(1-s^2)\log\frac{s+1}{s-1}.
$$
実際に $F(s)=\int_{-1}^1(1-u^2)\dd u/(s-u)$ であり、被積分関数は $s+u+(1-s^2)/(s-u)$ に分かれる。正の周波数の粒子–正孔過程は $\omega=k_zq/m+a$ を満たし、終状態が球の外にある必要がある。上の二項はこの Pauli 条件を実軸の虚部でも実装する。$\Im F(s+i0)=-\pi(1-s^2)$ が $-1<s<1$ の区間にだけ出るので、$\omega>a+b$ では虚部が零となる。$q>2k_F$ では $a-b>0$ となり、そこより下にも励起がない。

:::example[周波数帯の端]
$q=3k_F$ なら $a=9E_F$、$b=6E_F$ である。正の励起は $3E_F$ から $15E_F$ の範囲に入る。一方で $q=k_F$ なら $a=E_F,b=2E_F$ である。正の微小な周波数で二項の両方に虚部が出るため、差の重みは周波数とともに連続に零へ近づく。
:::

### 高運動量極限での反跳

$x=q/(2k_F)>1$ とし、$\log[(1+x)/(x-1)]=2\sum_{l\ge0}x^{-2l-1}/(2l+1)$ を :eqref[eq-mb-lindhard-static] に代入する。最初の定数は $1/2$ と対数の先頭項で相殺し、角括弧は $1/(3x^2)+1/(15x^4)+O(x^{-6})$ となる。従って $\chi_0(q,0)\sim-4mn/q^2$ である。これは占有球を動かして得る二分母に反跳 $q^2/2m$ が支配し、各分母から $-2mn/q^2$ が出るという解釈と一致する。

:::example[高波数の展開を途中で打ち切る]
$q=4k_F$ なら $x=2$、無次元の角括弧の厳密値は $1/2-3\log3/8\simeq0.0880$ である。先頭だけでは $1/12\simeq0.0833$、次まで入れると $7/80=0.0875$ となる。応答の符号はどの段階でも負である。
:::

### 連続体の内部での吸収強度

有限波数の式から、$\omega>0$、$q>0$ の虚部をそのまま求められる。$s_\pm=(\omega\pm a)/b$、$a=q^2/(2m)$、$b=v_Fq$ とする。$\Im F(s+i0)=-\pi(1-s^2)\theta(1-|s|)$ なので
$$
-\Im\chi_0(q,\omega+i0)=\frac{k_F^3}{4\pi b}
\left[(1-s_-^2)\theta(1-|s_-|)
      -(1-s_+^2)\theta(1-|s_+|)\right].
$$
二つの区間が重なると $s_+^2-s_-^2=4\omega a/b^2$ が残る。従って $q<2k_F$ かつ $0<\omega<b-a$ では、$-\Im\chi_0=m^2\omega/(2\pi q)$ である。$b-a<\omega<b+a$ では第一項だけが寄与し、上端で直線的に零へ向かう。$q>2k_F$ では第二項の区間は正周波数に入らず、第一項が $a-b<\omega<a+b$ にのみ残る。二項の差を取る操作は、両方の球の内側にある占有済みの終状態を除くことに対応する。

:::example[吸収強度の二つの式がつながる点]
$q=k_F$ とすると $a=E_F$、$b=2E_F$ なので境界は $\omega=E_F$ である。ここで $s_+=1$、$s_-=0$ となり、上側の式は $k_F^3/(4\pi b)=mk_F/(4\pi)$ である。下側の式も $m^2 E_F/(2\pi k_F)=mk_F/(4\pi)$ となる。連続体の上端 $\omega=3E_F$ では両項が消える。
:::

### 動的応答から自由気体の静的構造因子を求める

零温度では正周波数の $S_0(q,\omega)=-\Im\chi_0/\pi$ が粒子–正孔の数を測る。積分すると $\mathcal S_0(q)=n^{-1}\int_0^\infty S_0(q,\omega)\dd\omega$ は、各スピンで初期の Fermi 球の内部にあり、移動後の球の外に出る占有状態の割合となる。$q>0$ でこの割合は $1-B(q)/B(0)$ である。従って
$$
\mathcal S_0(q)=
\begin{cases}
3q/(4k_F)-q^3/(16k_F^3),&0<q<2k_F,\\
1,&q\ge2k_F .
\end{cases}
$$
球の重なりを数える代わりに、:eqref[eq-mb-jellium] の密度積を Wick 分解しても、$\langle \hat n_{\vb q}\hat n_{-\vb q}\rangle/N=(1/N)\sum_{\vb k\sigma}n_k(1-n_{\vb k+\vb q})$ となる。自己項を含む密度積をそのまま使うと高波数で $1$ に近づく。交換エネルギーを計算するときはこの $1$ を引くので、$E_x/N=\frac12\int v_q[\mathcal S_0(q)-1]\dd^3q/(2\pi)^3$ である。$q=0$ での数演算子の自明な弾性線は、$q\downarrow0$ の連結な構造因子とは区別する。

:::example[交換エネルギーを構造因子で積分する]
$q=k_F$ では $\mathcal S_0=11/16$ であり、欠ける割合は $5/16$ である。また球対称の交換エネルギーを動径積分に直すと $g_C/\pi$ 倍の $\int_0^{2k_F}[\mathcal S_0(q)-1]\dd q$ となる。多項式の三つの項の積分は順に $3k_F/2,-k_F/4,-2k_F$ なので、和は $-3k_F/4$ である。従って前節と同じ $-3g_Ck_F/(4\pi)$ が得られる。
:::

### 二次元の静的応答と角積分

電子を二次元に制限し、分散は引き続き $k^2/(2m)$、スピン縮退は $2$ とする。有限波数の静的応答を二つの占有円板へ移すと
$$
\chi_{0,\mathrm{2D}}(q,0)=-4\int_{k<k_F}\frac{k\dd k\,\dd\theta}{(2\pi)^2}
\operatorname{PV}\frac{1}{a+(qk/m)\cos\theta}.
$$
角積分は $a>qk/m$ なら $2\pi/\sqrt{a^2-(qk/m)^2}$、$0<a<qk/m$ なら主値が零である。これを確かめるには、$w=e^{i\theta}$ の単位円上に $a+i\delta$ を用いた分母を写す。二根の積は $1$ であり、内部の根の留数から角積分は $2\pi/\sqrt{(a+i\delta)^2-(qk/m)^2}$ となる。$\delta\downarrow0$ で後者の範囲は純虚数なので、実部を取ると零になる。残りの動径積分を $k\le\min(k_F,q/2)$ で行うと
$$
\chi_{0,\mathrm{2D}}(q,0)=-\frac{m}{\pi}
\left[1-\theta(q-2k_F)\sqrt{1-(2k_F/q)^2}\right].
$$
ここで平方根を使うのは $q>2k_F$ の枝だけである。有限な低波数区間で応答が一定なのは、放物線分散の二次元状態密度が一定であることと整合する。しきい値の上で平方根が立ち上がるため、三次元の場合とは $2k_F$ の非解析性の次数が違う。

:::example[二次元と三次元で同じ波数比を取る]
$q=2k_F$ では二次元の応答は $-m/\pi$ のままである。一方、三次元では $q\downarrow0$ の値の半分となる。$q=4k_F$ なら二次元の低波数値に対する比は $1-\sqrt3/2\simeq0.134$、三次元では $1/2-3\log3/8\simeq0.088$ である。二つの応答を数で比較するときは、面積あたりと体積あたりの違いを除くため、それぞれの状態密度で割っている。
:::

## 乱雑位相近似（RPA）と誘電関数

### 自己無撞着なポテンシャル

外から与えた電子一個あたりのポテンシャルエネルギーを $U_{\rm ext}$ とする。誘起密度による Hartree 項を足すと $U_{\rm tot}=U_{\rm ext}+v_q\delta n$ である。$\delta n=\chi_0U_{\rm tot}$ で閉じる :term[乱雑位相近似]{en="random phase approximation" yomi="らんざついそうきんじ"} から
$$
\varepsilon_{\rm RPA}=1-v_q\chi_0,\qquad
\chi_{\rm RPA}=\frac{\chi_0}{1-v_q\chi_0},\qquad
W=\frac{v_q}{\varepsilon_{\rm RPA}}.
\label{eq-mb-rpa-dielectric}
$$
これは自由バブルの鎖を全次数で足したものに等しい。分母が零に近いとき、有限個のバブルで打ち切ると集団極が失われる。交換や頂点補正はこの閉じ方から外れる。

:::example[静的な分母]
$v_q D_F=3$ かつ $\chi_0=-D_F$ なら $\varepsilon_{\rm RPA}=4$ となる。$W=v_q/4$、$\chi_{\rm RPA}=-D_F/4$ であり、正の外部ポテンシャルは密度を低下させる。
:::

:::example[正の静的ポテンシャルへの密度変化]
$v_q=2$、$\chi_0=-1/2$、$U_{\rm ext}=3$ を同じ単位系で使うと、誘電関数は $2$、全ポテンシャルは $3/2$、密度変化は $-3/4$ である。誘起項 $v_q\delta n=-3/2$ を足しても同じ全ポテンシャルに戻る。
:::

### 虚軸上の安定性と損失関数

零温度で $\chi_0(z)=\int_0^\infty S_0(E)[1/(z-E)-1/(z+E)]\dd E$ と書ける波数を取る。$S_0\ge0$ なら $\chi_0(i\nu)=-\int_0^\infty 2E S_0(E)/(\nu^2+E^2)\dd E\le0$ である。従って $v_q>0$ では虚軸上の $\varepsilon\ge1$ であり、リング積分の対数は実数になる。実軸で $\varepsilon=\varepsilon_1+i\varepsilon_2$ なら $-\Im(1/\varepsilon)=\varepsilon_2/|\varepsilon|^2$ であり、$\omega>0$ の $\varepsilon_2$ は $-v_q\Im\chi_0\ge0$ である。連続体外の単純零点では、有限の $\Im z$ を保ってから境界値を取ることで離散の損失ピークも含まれる。

:::example[静的応答が有限でも単純零点を持つ]
$\chi_0(z)=2E S/(z^2-E^2)$、$E,S>0$ とすると $\varepsilon^{-1}=1+2v_qES/(z^2-E^2-2v_qES)$ である。$E=2$、$v_qS=3$ なら、正の極は $4$、留数は $3/2$ である。$\varepsilon(0)=4$ の有限な遮蔽と、動的な離散の損失が同時に現れる。
:::

### 分極と応答で外場を区別する

電子の密度変化を $\delta n$、外部の電位を $\varphi_{\rm ext}$ とする。電子の電荷は $-e$ なので $U_{\rm ext}=-e\varphi_{\rm ext}$ である。誘起電荷 $-e\delta n$ に Poisson 方程式を適用すると $-e\varphi_{\rm ind}=v_q\delta n$ になる。従って $U_{\rm tot}=U_{\rm ext}+v_q\delta n$ は符号も含めて電気的な全場である。相互作用する系でも $\delta n=\Pi U_{\rm tot}$ により :term[既約分極]{en="proper polarization" yomi="きやくぶんきょく"} $\Pi$ を定義できる。$\chi=\delta n/U_{\rm ext}$ とすると、$\chi=\Pi/(1-v_q\Pi)$、$\Pi=\chi/(1+v_q\chi)$ である。図の言葉では、$\Pi$ の寄与は Coulomb 線一本で二つの密度応答に分けられないものを集める。RPA はこの $\Pi$ に自由値 $\chi_0$ を使う。

熱平衡で外から加える密度の一様変化を背景も打ち消すなら、バルクの中性気体の熱力学的な圧縮率は $\kappa_T=n^{-2}\partial n/\partial\mu$ である。有限波数の静的な試験場への応答では、電子だけが動き、背景は固定しているので Coulomb の $v_q$ が残る。このため RPA では $\Pi(q,0)\to-D_F$ で有限でも $\chi_{\rm RPA}(q,0)\sim -q^2/(4\pi g_C)$ になる。自由気体で $\partial n/\partial\mu=D_F$ を使えば、圧縮率の情報は RPA の分母を外した $\Pi$ にあると分かる。相互作用を含めた圧縮率と比較するには、分極に同じ近似の自己エネルギーと頂点を使う必要がある。

:::example[全場を逆算する]
$v_qD_F=9$、$\Pi=-D_F$ の静的なモードでは $\chi=-D_F/10$ である。逆にこの結果から既約分極を得ると、分母 $1+v_q\chi=1/10$ で割るため $-D_F$ に戻る。全場での $\delta n$ の比を直接測る場合と外から加えた場だけに対する比を測る場合で、十倍の違いが出る。
:::

### 斥力の RPA に上半平面の極がないことを確かめる

反転対称な安定な自由気体で、$S_0(q,E)$ を $E>0$ の規格化済み構造因子とする。$\chi_0(q,z)=\int_0^\infty 2E S_0(q,E)/(z^2-E^2)\dd E$ である。$z=x+iy$、$y>0$ とすると、$x\ne0$ での虚部は
$$
\Im \chi_0(q,z)=-2xy\int_0^\infty
\frac{2E S_0(q,E)}{|z^2-E^2|^2}\dd E.
$$
スペクトルの正周波数の重みが非零であれば符号は $-x$ の符号であり、$1/v_q$ のような正の実数に一致できない。$x=0$ なら $\chi_0(iy)\le0$ なのでやはり一致しない。重みが全て零なら $\varepsilon_{\rm RPA}=1$ である。従って $v_q>0$ では上半平面に誘電関数の零がなく、遅延応答の RPA 和はそこで解析的である。等比級数を各点で収束させる条件 $|v_q\chi_0|<1$ は、この解析的な逆関数の存在より強い。小波数では静的な $|v_q\chi_0|$ が大きくなり、自己無撞着な式から閉じた形を使う必要がある。

:::example[符号を反転したときの静的な極]
玩具模型 $\chi_0(z)=D_0\Omega_0^2/(z^2-\Omega_0^2)$、$D_0>0$ を取る。$1-v\chi_0=0$ は $z^2=\Omega_0^2(1+vD_0)$ を与える。斥力 $v>0$ なら極は上へ移る実振動である。引力 $v<0$ を仮定して $1+vD_0<0$ とすれば、根が虚軸に分かれ、上半平面の根は線形化した静止状態の不安定性を表す。この根を安定な系の遅延吸収線として扱うことはできない。
:::

## 遮蔽と Friedel 振動

### 指数遮蔽と長距離の振動

小さい $q$ で $\chi_0(q,0)\simeq-D_F$ とすると $W(q,0)=4\pi g_C/(q^2+\kappa^2)$、$\kappa^2=4\pi g_C D_F$ である。逆 Fourier 変換では角積分が $\sin(qr)/(qr)$ を与え、極 $i\kappa$ の留数を使うと $W(r)=g_Ce^{-\kappa r}/r$ となる。これは :term[Thomas–Fermi 遮蔽]{en="Thomas–Fermi screening" yomi="とーますふぇるみしゃへい"} である。

一方、$Q=2k_F$ の付近では :eqref[eq-mb-lindhard-static] の非解析部分が $-D_F(q-Q)\log|q-Q|/(2Q)$ である。滑らかな弱い点外場の $\delta n(q)=\chi_{\rm RPA}(q,0)U_{\rm ext}(q)$ にも $C(q-Q)\log|q-Q|$ が残る。逆変換では $q\sin(qr)/(2\pi^2r)$ が掛かる。$y=q-Q$ の奇関数部分に $\sin(yr)$ が掛かり、分布として $\int_\R\cos(yr)\log|y|\dd y=-\pi/r$（$r>0$）を微分して $\int_\R y\sin(yr)\log|y|\dd y=-\pi/r^2$ を得る。最初の恒等式は $\int_0^\infty e^{-\delta y}y^\alpha e^{iry}\dd y=\Gamma(1+\alpha)(\delta-ir)^{-1-\alpha}$ を $\alpha=0$ で微分して実部を取れば示せる。局所化した特異項から $\cos(2k_Fr)/r^3$ が得られ、これを :term[Friedel 振動]{en="Friedel oscillation" yomi="ふりーでるしんどう"} と呼ぶ。

:::example[距離と周期]
$k_F=10\,\mathrm{nm^{-1}}$ なら振動の周期 $\pi/k_F$ は約 $0.314\,\mathrm{nm}$ である。$\kappa=5\,\mathrm{nm^{-1}}$ の指数遮蔽長 $0.2\,\mathrm{nm}$ と由来が異なる。
:::

:::example[指数遮蔽の比]
$r=2/\kappa$ でのポテンシャルは、同じ距離の裸の $g_C/r$ に対し $e^{-2}\simeq0.135$ 倍である。$r=4/\kappa$ では比が $e^{-4}$ になる。距離どうしで全ポテンシャルを比較する場合には、さらに $1/r$ の違いを含める。
:::

### 外部電荷を積分したときの中性化

Thomas–Fermi 近似で点電荷 $Q$ の全電位を SI で $\varphi(r)=Qe^{-\kappa r}/(4\pi\epsilon_0 r)$ とする。$r>0$ で放射状の Laplacian を取ると $\nabla^2\varphi=\kappa^2\varphi$ なので、誘起電荷密度は $\rho_{\rm ind}=-\epsilon_0\kappa^2\varphi$ である。半径 $R$ 内の誘起電荷は
$$
Q_{\rm ind}(R)=-Q\int_0^R\kappa^2 r e^{-\kappa r}\dd r=-Q[1-(1+\kappa R)e^{-\kappa R}] .
$$
$R\to\infty$ で合計は $-Q$ になる。等価な確認は Fourier の $q\to0$ で $\rho_{\rm ind}(\vb q)=-Q\kappa^2/(q^2+\kappa^2)$ を評価することである。一様成分を最初から除く有限箱では、局所応答からこの極限を取る順序を明記する。

:::example[遮蔽長内に集まる電荷]
$R=1/\kappa$ では誘起電荷の絶対値は $|Q|(1-2/e)\simeq0.264|Q|$ である。$R=2/\kappa$ で約 $0.594|Q|$ になる。半径で積分すると $r^2$ の体積重みがあるため、電位が $1/e$ へ減る距離と全電荷が集まる距離は一致しない。
:::

### 点状の試験電荷がつくる振動の係数

$Q_0=2k_F$ とし、外部の電子一個あたりのポテンシャル $U_{\rm ext}(q)$ は球対称で $Q_0$ の近くで滑らかとする。$\chi_0(q,0)$ の特異項は $c(q-Q_0)\log(|q-Q_0|/Q_0)$、$c=-D_F/(2Q_0)$ である。対数の中で使う波数の単位を変えても、差は $q-Q_0$ に比例する滑らかな項になる。$\chi_{\rm RPA}=\chi_0/(1-v_q\chi_0)$ を $\chi_0$ で微分すると $1/(1-v_q\chi_0)^2$ なので、誘起密度の Fourier 成分 $\delta n=\chi_{\rm RPA}U_{\rm ext}$ の特異項の係数は
$$
C=-\frac{D_F\,U_{\rm ext}(Q_0)}
{2Q_0\,\varepsilon_{\rm RPA}(Q_0,0)^2}.
$$
動径逆変換は $(2\pi^2r)^{-1}\int_0^\infty q\sin(qr)\delta n(q)\dd q$ である。$q=Q_0+y$ とすると、特異項の先頭では $q$ の前因子を $Q_0$ に置ける。$\sin(qr)=\sin(Q_0r)\cos(yr)+\cos(Q_0r)\sin(yr)$ であり、$y\log|y|$ は奇なので第二項だけが寄与する。前に求めた変換 $-\pi/r^2$ を使って
$$
\delta n(r)=
\frac{D_F U_{\rm ext}(Q_0)}{4\pi\,\varepsilon_{\rm RPA}(Q_0,0)^2}
\frac{\cos(Q_0 r)}{r^3}+o(r^{-3})
$$
を得る。この式では外場の原点周辺にある局在部を除き、$Q_0$ の特異点以外の滑らかな部分は十分速く減衰すると仮定した。例えば $U_{\rm ext}(q)$ が点電荷の形で $q=0$ に $1/q^2$ を持つときも、$\chi_{\rm RPA}\sim-q^2/(4\pi g_C)$ がこれを打ち消し、誘起密度の動径積分の被積分関数は原点の近傍で $q$ の滑らかな奇関数となる。そこで境界項まで含めて部分積分すると、原点から $r^{-3}$ の項は生じない。

:::example[ポテンシャルの振動と密度の振動]
電子一個を試験源として $U_{\rm ext}=v_q$ とする。全ポテンシャル $W=v_q+v_q\delta n(q)$ の $Q_0$ における特異項には、密度の場合にさらに $v_{Q_0}$ が掛かる。従って $W(r)$ の振動の係数は $D_F v_{Q_0}^2/[4\pi\varepsilon(Q_0,0)^2]$ であり、$U_{\rm ext}$ の符号を逆に取った正の試験電荷では誘起された電子数密度の振動の位相が $\pi$ だけ変わる。小波数の Yukawa 形は距離が十分長いと指数で減るため、振動の節を除く遠方の包絡はこの代数項が担う。
:::

### 二次元で減衰のべきが変わる理由

二次元で $Q_0=2k_F$ とすると、先に求めた $\chi_{0,\mathrm{2D}}$ の非解析部分は $+(m/\pi)\sqrt{2(q-Q_0)/Q_0}\,\theta(q-Q_0)$ である。角積分は $J_0(qr)=(2\pi)^{-1}\int_0^{2\pi}e^{iqr\cos\theta}\dd\theta$ なので、逆変換の測度は $q J_0(qr)\dd q/(2\pi)$ となる。$qr\gg1$ で位相の停留点は $\theta=0,\pi$ の二つであり、各点で $\cos\theta$ を二次まで展開して振動 Gauss 積分を取ると $J_0(qr)\sim\sqrt{2/(\pi qr)}\cos(qr-\pi/4)$ を得る。二点の位相はそれぞれ $qr-\pi/4$ とその逆である。

次に $y=q-Q_0\ge0$ とし、滑らかな切断を掛けた $y^{1/2}$ を逆変換する。減衰因子 $e^{-\delta y}$ を持つ積分は $\int_0^\infty y^{1/2}e^{(ir-\delta)y}\dd y=\Gamma(3/2)(\delta-ir)^{-3/2}$ であり、$\delta\downarrow0$ で先頭に $e^{3i\pi/4}r^{-3/2}$ を与える。特異点を滑らかに局所化しても、その近傍の先頭係数は同じである。角積分の $r^{-1/2}$ と合わせて振動は $r^{-2}$ となり、位相は $\cos(Q_0r-\pi/4+3\pi/4)=-\sin(Q_0r)$ である。連続体の小波数での非振動項は別にあり、ここで求めたのは $2k_F$ に由来する寄与である。

:::example[距離を二倍にして包絡を比べる]
同じ $k_F$ で振動の包絡を見ると、距離を二倍にして二次元では $1/4$、三次元では $1/8$ になる。どちらも波数 $2k_F$ なので、距離を変えた二点の値そのものは位相によってこの比から外れ得る。包絡の次元依存は、角積分と $2k_F$ の特異点の次数が両方変わることで生じる。
:::

## プラズモン

### 誘電関数の零点

高周波の $\chi_0\simeq nq^2/(mz^2)$ を :eqref[eq-mb-rpa-dielectric] に入れると、零点は $\omega_p^2=4\pi g_C n/m$ になる。この密度の集団振動を :term[プラズモン]{en="plasmon" yomi="ぷらずもん"} と呼ぶ。$q\to0$ で $\omega_p$ が有限なので、粒子–正孔連続体の $\omega\lesssim qv_F$ から離れる。次の係数は :eqref[eq-mb-lindhard-dynamic] の $D_F/(5s^4)$ から出て、$\omega^2=\omega_p^2+3v_F^2q^2/5+O(q^4)$ となる。厳密な放物線 Lindhard 関数での反跳項はここでは $q^4$ 以降に寄与する。実際、占有球上で分母の差を展開すると $2(q^2/2m)z^{-2}$ と $6(q^2/2m)(\vb k\cdot\vb q/m)^2z^{-4}$ が最初の二項であり、球平均 $\langle k_z^2\rangle=k_F^2/5$ が上の係数を与える。

:::example[プラズマエネルギー]
真空誘電率で $n=10^{29}\,\mathrm{m^{-3}}$ とすると SI の $\sqrt{ne^2/(m\epsilon_0)}$ は約 $1.8\times10^{16}\,\mathrm{s^{-1}}$、$\hbar\omega_p$ は約 $11.7\,\mathrm{eV}$ である。
:::

:::example[小さい波数の分散補正]
$v_Fq=\omega_p/2$ の点に小波数展開を代入すると $\omega^2/\omega_p^2=1+3/(5\cdot4)=23/20$ となる。一次の周波数補正としては $\omega/\omega_p\simeq1+3/40$ である。二乗からの平方根を厳密に取っても、もとの $q^4$ の項は復元されない。
:::

### 二次元の長距離 Coulomb 則

電子が平面に閉じ込められ、電場が三次元へ広がるとする。平面 Fourier 変換後、点電荷の電位は平面外で $(\partial_z^2-q^2)\varphi=0$ を満たし、遠方で消えるので $e^{-q|z|}$ に比例する。電荷平面で Poisson 方程式を積分すると $\partial_z\varphi(0^+)-\partial_z\varphi(0^-)=-Q/\epsilon_0$ である。よって $v_q=e^2/(2\epsilon_0 q)$ となる。電子の連続の式と運動方程式、または密度の総和則から $\chi_0\simeq n_{\rm 2D}q^2/(m\omega^2)$ なので、
$$
\omega_{\rm 2D}^2\simeq \frac{e^2 n_{\rm 2D} q}{2m\epsilon_0}.
$$
$\,\omega_{\rm 2D}\propto\sqrt q$ であり、長波長で零へ近づく。それでも粒子–正孔の上端は $q$ に比例するため、十分小さい $q$ では集団モードは上端の外にある。

:::example[波長による周波数比]
長波長の先頭次数で比べる。同じ平面密度で波長を四倍にすると二次元プラズモンの周波数は半分になる。粒子–正孔の上端の先頭項は四分の一になるので、上端に対する周波数比は二倍となる。面密度 $n_{\rm 2D}$ を体積密度へ割り戻して三次元公式に入れると、この波数依存を失う。
:::

### 長波長分散の次の係数と反跳

占有球の運動量積分を $\langle\cdot\rangle_{\rm ball}$ で規格化し、$d=\vb k\cdot\vb q/m$ と書く。二つの分母に運動量移動を吸収すると
$$
\chi_0(q,z)=n\left\langle \frac1{z-d-a}-\frac1{z-d+a}\right\rangle_{\rm ball},
\qquad a=\frac{q^2}{2m}.
$$
$|z|>v_F q+a$ で両分母を幾何級数で展開する。偶対称な球平均を使うと、$z^{-2},z^{-4},z^{-6}$ の分子はそれぞれ $2a$、$6a\langle d^2\rangle+2a^3$、$10a\langle d^4\rangle+20a^3\langle d^2\rangle+2a^5$ となる。球の断面積に比例する密度 $3(1-u^2)/4$ で $u=k_\parallel/k_F$ を平均すると $\langle u^2\rangle=1/5$、$\langle u^4\rangle=3/35$ である。従って $z$ を $\omega_p$ の近くに保った展開は
$$
\chi_0(q,z)=\frac{nq^2}{mz^2}\left[
1+\frac{3v_F^2q^2/5+q^4/(4m^2)}{z^2}
+\frac{3v_F^4q^4}{7z^4}+O(q^6)
\right].
$$
$v_q\chi_0(q,\omega_q)=1$ の両辺に $\omega_q^2$ を掛け、$\omega_q^2=\omega_p^2+Aq^2+Bq^4+O(q^6)$ と仮定する。右辺の $1/\omega_q^2$ を二次まで、$1/\omega_q^4$ を零次まで展開して係数を合わせると $A=3v_F^2/5$、$B=1/(4m^2)+(3/7-9/25)v_F^4/\omega_p^2$ になる。よって
$$
\omega_q^2=\omega_p^2+\frac35 v_F^2 q^2+
\left(\frac{1}{4m^2}+\frac{12 v_F^4}{175\omega_p^2}\right)q^4+O(q^6).
$$
$q=0$ で $1-\omega_p^2/z^2$ の零点の微分は非零なので、十分小さい $q$ でこの展開に対応する実の根が存在する。$1/(4m^2)$ は反跳に由来する。長波長で $z/(v_Fq)$ を固定した Lindhard 近似からはこの項が出ないため、四次まで求めるには二つの分母に $a$ を残す。

:::example[反跳と Fermi 速度の補正を見分ける]
$E_F=mv_F^2/2$ を使うと、$\omega_q^2$ の $q^4$ 係数における Fermi 速度由来の項と反跳項の比は $192 E_F^2/(175\omega_p^2)$ である。例えば $\omega_p=2E_F$ なら比は $48/175$ になる。いずれの項を足す場合も展開自体には $v_Fq+q^2/(2m)\ll\omega_p$ が必要である。
:::

### 離散プラズモンが連続体の端に達する波数

$q>0$ で連続体の上端を $E_+(q)=a+b$ とする。$z>E_+(q)$ の実軸では $\chi_0(q,z)=\int_0^\infty 2E S_0(q,E)/(z^2-E^2)\dd E$ が正で、$z$ を増やすと厳密に減少して零へ向かう。従って $\varepsilon=0$ による上端の外の根がある条件は $v_q\chi_0(q,E_++0)>1$ であり、その場合には根は一つである。有限波数の Lindhard 式で $F(1)=2$ だから、端の値は
$$
v_q\chi_0(q,E_++0)
=v_q\frac{k_F^3}{4\pi^2b}[F(1)-F(1+q/k_F)]
=\frac{g_C m k_F}{\pi q^2}
\int_{-1}^1\frac{1+u}{1+q/k_F-u}\dd u .
$$
最後の式は $F$ を $(1-u^2)/(s-u)$ の積分に戻し、分母の差から $q/k_F$ を取り出したものである。$q$ が増えると $q^{-2}$ と積分の両方が減る。$q\downarrow0$ では発散し、$q\to\infty$ では零に近づくので、端の値が $1$ となる波数 $q_*$ は一つである。$q<q_*$ では離散の RPA プラズモンが上端の外にあり、$q>q_*$ ではこの枝に対応する上端の外の実極は消える。連続体に入った周波数では粒子–正孔対への :term[Landau 減衰]{en="Landau damping" yomi="らんだうげんすい"} によりスペクトルに幅が生じる。

周波数の実部 $\omega_*$ で $\varepsilon_1(\omega_*)=0$、傾き $s=\partial_\omega\varepsilon_1(\omega_*)>0$ とし、幅の範囲で $\varepsilon_2$ がほぼ一定の弱減衰の場合を取る。$\varepsilon(\omega)\simeq s(\omega-\omega_*)+i\varepsilon_2(\omega_*)$ なので、損失関数の半値半幅は $\varepsilon_2/s$ である。$\varepsilon_2=-v_q\Im\chi_0$ は正周波数で非負であり、これが対の生成による減衰の符号を決める。しきい値直上や重なった広いピークでは、この一次近似が必要とする滑らかさを確かめる。

:::example[境界を一つの波数で判定する]
$q=k_F$ では積分は $3\log3-2$ になる。このとき $v_q\chi_0(E_++0)=(3\log3-2)/(\pi k_Fa_B)$ である。例えば $k_Fa_B=2$ では約 $0.206$ なので、端の値の単調性により $q_*<k_F$ となる。$\omega_p$ をそのまま水平な分散とみなして $v_Fq=\omega_p$ を解く方法では、反跳と有限波数での Lindhard の形が使われない。
:::

### 密度スペクトルの重みと損失の重み

連続体外に単純な零 $\omega_q$ を持つ場合、$\chi_{\rm RPA}=\chi_0/\varepsilon$ の正の極の留数は $Z_n=1/[v_q\,\varepsilon'(q,\omega_q)]$ である。$\chi_0(q,\omega_q)=1/v_q$ を用いた。逆誘電関数の留数は $v_q Z_n$ であり、遮蔽相互作用 $W=v_q/\varepsilon$ の留数は $v_q^2Z_n$ になる。測定量によって同じ周波数のピークでも積分重みが違う。

$q\downarrow0$ では $\varepsilon'(q,\omega_q)\to2/\omega_p$ なので
$$
Z_n=\frac{nq^2}{2m\omega_p}+O(q^4),\qquad
\omega_q Z_n=\frac{nq^2}{2m}+O(q^4).
$$
ここで最初の式に $v_q=4\pi g_C/q^2$ と $\omega_p^2=4\pi g_Cn/m$ を使った。一方、$\chi_{\rm RPA}$ の高周波での先頭は $\chi_0$ と同じ $nq^2/(m z^2)$ である。従って 式 :eqref[eq-mb-density-fsum] の一次モーメントの先頭は、長波長の RPA ではプラズモン一つで担われる。$\int_0^\infty S(q,\omega)\dd\omega/n$ の先頭が $q^2/(2m\omega_p)$ になることは、連続体の重みの見積もりからも分かる。実際、$q\ll k_F$ では $0\le\omega\le a+b$ において $\chi_0$ は $D_F$ 程度かそれ以上の絶対値を持つ。従って $|\varepsilon|^{-2}=O(q^4)$ により $S_{\rm RPA}\le O(q^4)S_0$ となり、自由な連続体の総重み $O(q)$ は $O(q^5)$ 以下になる。ここでは小波数の動的 Lindhard 式の実部と虚部が同時に零にならず、上端近くでは実部が大きいことを用いた。

:::example[密度ピークが弱くても損失ピークが残る]
波数を半分にすると長波長の $Z_n$ は約四分の一になる。一方で $v_q$ は四倍になり、逆誘電関数の留数 $v_q Z_n$ は $\omega_p/2$ へ近づく。電荷の密度に結合した分光と、誘電関数の逆数を通じたエネルギー損失では、この $v_q$ の因子が異なる。
:::

## 相関エネルギーと Gell-Mann–Brueckner の結果

### リング和と対数領域

相関エネルギーは同じ密度の運動エネルギーと Hartree–Fock 交換を基底エネルギーから引いた値である。結合定数を $\lambda v_q$ に変え、非縮退の基底枝の固有方程式を $\lambda$ で微分する。基底射影を左から掛けるとエネルギーの微分は相互作用の期待値になり、:eqref[eq-mb-jellium] により $\frac12\int v_q[\int_0^\infty S_\lambda(q,E)\dd E-n]\dd^3q/(2\pi)^3$ である。ここで $S_\lambda$ は単位体積あたりで、$q\ne0$ の一様基底の密度期待値は零とした。自由値で同じ期待値を求めると一次の交換項になるため、これを引くと自己項 $-n$ も消える。

虚軸では $\chi_\lambda(q,i\nu)=-\int_0^\infty 2E S_\lambda(q,E)/(\nu^2+E^2)\dd E$ である。$2E/(\nu^2+E^2)$ の $\dd\nu/(2\pi)$ での積分は $1$ なので、連結密度積は虚周波数の応答の負の積分になる。これに $\chi_\lambda=\chi_0/(1-\lambda v_q\chi_0)$ を入れて自由値を引き、$\lambda$ を積分する。$\int_0^1 [-v_q\chi_0/(1-\lambda v_q\chi_0)+v_q\chi_0]\dd\lambda=\log(1-v_q\chi_0)+v_q\chi_0$ なので
$$
\frac{E_{\rm ring}}{V}=\frac12\int\frac{\dd^3q}{(2\pi)^3}
\int_{-\infty}^{\infty}\frac{\dd\nu}{2\pi}
\{\log[1-v_q\chi_0(q,i\nu)]+v_q\chi_0(q,i\nu)\}.
\label{eq-mb-ring-energy}
$$
括弧を展開すると $-[(v_q\chi_0)^2/2+(v_q\chi_0)^3/3+\cdots]$ であり、一次は除かれている。二次だけでは小さい $q$ の積分が発散する。

$\kappa\ll q\ll k_F$ で $v_q\chi_0=-\kappa^2R(u)/q^2$、$u=|\nu|/(v_Fq)$、$R(u)=1-u\arctan(1/u)$ と置く。二次項、球の測度、周波数の Jacobian を掛けると粒子一個あたりの対数部分は
$$
-\frac{v_F\kappa^4}{8\pi^3 n}\log(k_F/\kappa)
\int_0^\infty R(u)^2\dd u .
$$
$R(u)=\int_0^1 x^2/(u^2+x^2)\dd x$ であり、部分分数で $u$ を先に積分すると $\frac\pi2\int_0^1\int_0^1 xy/(x+y)\dd x\dd y$ となる。領域を $y<x$ に分けて $y=xt$ を代入すると $\pi(1-\log2)/3$ を得る。$\kappa/k_F\propto\sqrt{r_s}$ なので
$$
\frac{E_c}{N E_{\rm H}}=\frac{1-\log2}{\pi^2}\log r_s+B+o(1).
\label{eq-mb-gell-mann-brueckner}
$$
これが非偏極三次元電子ガスに対する Gell-Mann–Brueckner の高密度結果である。定数 $B$ には直接リングの有限部分と二次交換図が入る。リング以外の高次の図がここで示した対数と定数を変えないことの証明は、遮蔽した構成で複数の内部運動量が同時に軟化する領域の次数評価を必要とするため、:cite[fetter-walecka] の高密度展開を参照する。本章では、これら二つの寄与を積分で計算する。

:::example[対数項の大きさ]
$r_s=0.1$ なら対数項は約 $-0.0716\,E_{\rm H}$ である。これに定数と残りが加わり、式 :eqref[eq-mb-jellium-exchange] の交換項とも別に数える。
:::

:::example[リング積分の符号]
一つの虚周波数で $v_q\chi_0=-1$ なら、リングの括弧は $\log2-1<0$ である。一般に $x=-v_q\chi_0\ge0$ に対し $\log(1+x)-x=-\int_0^x t/(1+t)\dd t\le0$ なので、各周波数で同じ符号を持つ。
:::

### 定数項と二次の交換

定数項は Hartree 単位で $B\simeq-0.047$ である。直接リングの有限部分は約 $-0.071$、二次交換の有限部分は $B_{\rm x}^{(2)}=\log2/6-3\zeta(3)/(4\pi^2)\simeq+0.024179$ であり、両者を足して $B$ を得る。$\zeta(3)=\sum_{l\ge1}l^{-3}$ である。二次交換の閉じた式の解析的評価は積分領域の分割と特殊関数の恒等式を要するため :cite[fetter-walecka] を参照する。以下では直接リングと二次交換の両方について数値評価できる積分を導く。Rydberg をエネルギー単位に使えば、対数と定数の両方の係数が二倍になる。

二次交換がどこから出るかは自由 Fermi 海で追える。摂動を一度作用させた励起状態の振幅は、二粒子二正孔の反対称化行列要素を正の励起費用 $\Delta$ で割ったものの負である。これをもう一度相互作用で基底へ戻すと、二次エネルギーは $-|\bar v|^2/\Delta$ を励起について足したものになる。占有の二添字と空の二添字を順序付きで足すなら各励起を四回数えるため係数は $1/4$ になる。占有運動量を $\vb k,\vb p$、空の運動量を $\vb k+\vb q,\vb p-\vb q$ とすれば
$$
\Delta=\frac{q^2+(\vb k-\vb p)\cdot\vb q}{m},\qquad
v_{\rm d}=v_{\vb q},\quad v_{\rm x}=v_{\vb k+\vb q-\vb p}.
$$
スピン和では直接項の平方が $4v_{\rm d}^2+4v_{\rm x}^2$、交差項は両占有スピンが同じ場合の $-4v_{\rm d}v_{\rm x}$ である。従って正の交換エネルギー密度は
$$
\frac{E_{\rm x}^{(2)}}V=
\int\frac{\dd^3k\,\dd^3p\,\dd^3q}{(2\pi)^9}
\frac{n_k n_p(1-n_{\vb k+\vb q})(1-n_{\vb p-\vb q})}{\Delta}\,
v_{\rm d}v_{\rm x}.
$$
ここで $n_k$ は一スピンの占有である。正の符号は、反対称化した振幅の積の負符号と、二次補正の負符号が相殺するためである。二次交換では直接と交換の二つの転送が異なるため、直接リングの赤外対数をそのまま適用できない。下で占有条件の境界を積分端点に移し、有限の定数を取り出す。

:::example[対数と有限部分を足す]
$r_s=0.1$ なら対数と定数を足した相関エネルギーは約 $-0.119\,E_{\rm H}$ である。ここには $r_s\to0$ で消える高次項を含めていない。Hartree–Fock の交換は既に基準値に入っているため、二次交換の正の項と区別する。
:::

### 二次リングの赤外発散を有限箱で追う

小波数領域での二次の直接項は $\int_{q_{\rm min}}^{q_{\rm max}}\dd q/q$ を含む。箱の辺を $L$、周期境界とすると最低の非零波数は $q_{\rm min}=2\pi/L$ であり、二次項だけなら粒子一個あたりの値に $\log L$ が残る。一方、和を取った括弧では $x=-v_q\chi_0$ として $\log(1+x)-x$、$x\gg1$ での主要項は $-x$ である。固定した $u=|\nu|/(qv_F)$ で $x=O(q^{-2})$ なので、測度 $q^2\dd q\,q\dd u$ は低波数で $\dd q\,q$ にする。$u$ 積分を含めても $0\le x-\log(1+x)\le x$、$\int_0^\infty R(u)\dd u=\pi/4$ で可積分である。最後の積分は $R=\int_0^1 x^2/(u^2+x^2)\dd x$ に対し $u$ を先に積分して得る。従って赤外の切断は $2\pi/L$ から、相互作用が分母を変える $\kappa$ の周辺に置き換わる。

:::example[同じ打ち切りの長さ依存性]
二次項だけで $L$ を二倍にすると積分の下端の変化は $\log2$ である。和を取った後で $2\pi/L\ll\kappa$ なら、低波数側の欠損は上の評価から $O(L^{-2})$ 以下になる。この上界は赤外の連続積分を切る効果に対するもので、離散格子和の全誤差を与えるものではない。
:::

### リング積分の有限部分を取り出す

$y=q/k_F$、$u=\nu/(v_Fq)\ge0$、$h=\kappa/k_F$ として $Q(y,u)=-\chi_0(q,i\nu)/D_F$ を導入する。有限波数の Lindhard 式では $F(-\bar s)=-\overline{F(s)}$ なので
$$
Q(y,u)=\frac{\Re F(y/2+iu)}{2y},\qquad
Q(0,u)=R(u)=1-u\arctan(1/u).
$$
対数の枝は $u>0$ で上半平面から評価し、$u=0$ はその極限で読む。$\nu$ の正負が等しく、球の測度が $q^2\dd q/(2\pi^2)$ であることから、:eqref[eq-mb-ring-energy] は
$$
\frac{E_{\rm ring}}{N E_{\rm H}}
=\frac{12}{\pi^3 h^4}\int_0^\infty y^3\dd y\int_0^\infty\dd u\,
\left[\log\left(1+\frac{h^2 Q(y,u)}{y^2}\right)-\frac{h^2 Q(y,u)}{y^2}\right]
$$
となる。前因子では $v_F\kappa^4/(4\pi^3 n E_{\rm H})=12/\pi^3$ を用いた。自由な気体の $R$ だけを全波数に使うと、有限部分が任意の上限に依存するため、上の $Q$ が必要である。

$h\ll y_0\ll1$ に分け、硬い側 $y\ge y_0$ では $h\to0$ を先に取る。$0\le x-\log(1+x)\le x^2/2$ により、$h^4$ で割った被積分関数は $Q(y,u)^2/(2y)$ で支配される。$y\ge4$ では $F(s)=\int_{-1}^1(1-t^2)/(s-t)\dd t$ を使うと $Q(y,u)\le C/(y^2+u^2)$ になる。$y$ が $y_0$ と $4$ の間でも、$Q$ は $u$ が大きいと $O(u^{-2})$ であり、有限周波数の範囲で有界である。従って硬い側は極限で $-\frac12\int_{y_0}^\infty\dd y\,I(y)/y$、$I(y)=\int_0^\infty Q(y,u)^2\dd u$ となる。

軟らかい側 $y<y_0$ では $Q$ を $R$ に替える。実際、$F'(iu)=4R(u)$、$F'''(s)=-8/(s^2-1)^2$ であり、$|\Re s|\le y_0/2$ では三階微分を $C/(1+u^2)^2$ で抑えられる。上半平面の境界で $\Re F(iu)=0$、$\Re F''(iu)=0$ を使うと、$|Q-R|\le C y^2/(1+u^2)^2$ である。関数 $\log(1+x)-x$ の微分の絶対値は $x/(1+x)\le x$ であり、$R$ と $Q$ は $C/(1+u^2)$ 以下である。したがって二つの括弧の差に $y^3/h^4$ を掛けた絶対値は $Cy/(1+u^2)^3$ 以下になる。$u$ と $0<y<y_0$ を積分すると $O(y_0^2)$ なので、$h\to0$ の後で $y_0\to0$ として消せる。

軟らかい側で $p=y/h$、$t=p^2$ に置き換える。固定した $R>0$ に対する $p$ の積分の原始関数は
$$
\frac{t^2}{4}\log(1+R/t)-\frac{Rt}{4}-\frac{R^2}{4}\log(1+t/R).
$$
$t=0$ で零となり、$t=L^2=(y_0/h)^2$ で $-\frac12 R^2\log L+\frac14R^2\log R-\frac18 R^2+O(R^3/L^2)$ となる。$0<R\le1$ なので、この余りは $L\to\infty$ で $u$ 積分した後も消える。硬い側からは $I(0)=I_0=\pi(1-\log 2)/3$ を $0<y<1$ だけ引いて、残りを $y_0=0$ から積分できる。二領域の $\log y_0$ を相殺し、$h^2=4r_s/(\pi b_0)$、$b_0=(9\pi/4)^{1/3}$ を代入すると
$$
\begin{aligned}
B_{\rm dir}=\frac{12}{\pi^3}\Bigg[
&\frac{I_0}{4}\log\frac{4}{\pi b_0}
+\frac14\int_0^\infty R(u)^2\log R(u)\dd u-\frac{I_0}{8}\\
&-\frac12\int_0^\infty\frac{I(y)-I_0\theta(1-y)}{y}\dd y\Bigg].
\end{aligned}
\label{eq-mb-ring-constant}
$$
この表示では $I(y)-I_0=O(y^2)$ なので、最後の積分の左端も収束する。数値計算では $0<y<1$ と $1<y<\infty$ を分け、各無限区間には $x=t/(1-t)$ などの変換を用いる。小さい $y$ で二つの近い値を直接差し引くと丸め誤差が出るため、必要なら $Q=R+O(y^2)$ の差の方を $F'''$ から評価する。この積分を数値評価すると $B_{\rm dir}\simeq-0.07110$ となる。Gell-Mann–Brueckner の定数は、これに次の二次交換の値を足して得る。

### 二次交換を動径と縦成分に分ける

転送 $\vb q$ を極軸とし、以後の波数を $k_F$ で割る。$x=q/k_F$、初期状態の縦成分を $k_\parallel/k_F=u-x/2$、$p_\parallel/k_F=x/2-v$ とする。横成分の二乗をそれぞれ $R_1,R_2$ とする。粒子が球の内から外へ動くためには $u,v>0$ であり、各縦成分 $w=u,v$ の範囲は
$$
w_{\min}=\max(0,x/2-1)<w<x/2+1=w_{\max}.
$$
初期が球内で終状態が球外という二つの条件を横成分の範囲に直すと
$$
R_-(w)=\max[0,1-(w+x/2)^2]<R<R_+(w)=1-(w-x/2)^2
$$
である。二励起の費用は $\Delta=k_F^2x(u+v)/m$ で、交換側の転送の縦成分は $k_F(u+v)$ になる。横成分の二方向がなす角を $\theta$ とすると、交換転送の二乗を $k_F^2$ で割った分母は $(u+v)^2+R_1+R_2-2\sqrt{R_1R_2}\cos\theta$ である。角積分は、二次元の静的応答で使った Gauss の角積分と同じ形である。これを評価し、二つの横成分の動径測度 $\rho\dd\rho=\dd R/2$ を使うと
$$
\begin{aligned}
B_{\rm x}^{(2)}=\frac{3}{8\pi^2}
\int_0^\infty\frac{\dd x}{x}\int_{w_{\min}}^{w_{\max}}\dd u
\int_{w_{\min}}^{w_{\max}}\frac{\dd v}{u+v}
\int_{R_-(u)}^{R_+(u)}\dd R_1\int_{R_-(v)}^{R_+(v)}
\frac{\dd R_2}{\sqrt{[(u+v)^2+R_1+R_2]^2-4R_1R_2}}.
\end{aligned}
\label{eq-mb-exchange-constant}
$$
前因子を得るには、$\vb q$ の全立体角 $4\pi$、二つの横角の積分の係数 $(2\pi)^2$、各 Coulomb 核の $4\pi g_C$ を掛け、$(2\pi)^9$ と二つの動径の $4$ で割る。$m g_C^2=E_{\rm H}$、$n=k_F^3/(3\pi^2)$ を使うと $3/(8\pi^2)$ になる。内側の $R_2$ 積分は $t=u+v>0$、$R_1>0$ に対し、原始関数
$$
\operatorname{arsinh}\left(\frac{R_2+t^2-R_1}{2t\sqrt{R_1}}\right)
$$
を端で評価すればよい。$R_1=0$ では $\log(t^2+R_2)$ を使える。従って数値積分は $x,u,v,R_1$ の四変数にできる。$R_-$ の式が変わる $w=1-x/2$ と、$w_{\min}$ の式が変わる $x=2$ で区間を分ければ、各区間の境界が明示される。

この式が赤外で有限な理由も横成分の範囲から読める。$x\to0$ では $R_+(w)-R_-(w)\le2xw$ であり、交換核の逆平方根は $1/(u+v)^2$ 以下である。従って固定した $x$ の $u,v$ の積分は、定数倍の $x^2\int_0^{1+x/2}\int_0^{1+x/2}uv/(u+v)^3\,\dd u\dd v$ で抑えられ、両方が零に近い部分でも有限である。外の $\dd x/x$ を掛けても積分できる。一方、$x\gg1$ では $u,v$ の積分の幅は $2$、$u+v$ は $x$ 程度であり、交換核は $O(x^{-2})$ になるので $x$ の積分は $x^{-4}$ で収束する。直接二次の転送二乗の発散とこの積分の違いは、交換転送が独立な角度に依存する点にある。

:::example[相関の運動部分とポテンシャル部分]
$e_c/E_{\rm H}=A\log r_s+B$ の先頭二項を使う。座標変換から導いた式を相関部分に適用すると、運動の補正は $e_{{\rm kin},c}/E_{\rm H}=-A\log r_s-B-A$、ポテンシャルの補正は $e_{{\rm pot},c}/E_{\rm H}=2A\log r_s+2B+A$ である。ここでは高密度級数を項別に微分できる次数まで用いた。例えば $r_s=0.1$、$A\simeq0.0311$、$B\simeq-0.047$ なら運動は約 $+0.088$、ポテンシャルは約 $-0.206$ となる。和が約 $-0.119$ になるのは、相関が運動エネルギーを増やしても、電子間のポテンシャルをそれ以上に下げるからである。
:::

## 演習問題

::::exercise{level=basic}
非偏極で $k_F=2/a_B$ の $E_F/E_{\rm H}$ を求めよ。

:::solution
$E_F=k_F^2/(2m)$、$E_{\rm H}=1/(ma_B^2)$ より比は $2^2/2=2$ である。
:::
::::

::::exercise{level=basic}
静的 Lindhard 関数の $q=2k_F$ の値を $D_F$ で表せ。

:::solution
$(1-x^2)\log|(1+x)/(1-x)|$ は $x\to1$ で零になるので、角括弧は $1/2$ である。値は $-D_F/2$ である。
:::
::::

::::exercise{level=basic}
$\kappa=4\,\mathrm{nm^{-1}}$ の遮蔽長を求めよ。

:::solution
指数 $e^{-\kappa r}$ が $e^{-1}$ になる距離 $r=1/\kappa=0.25\,\mathrm{nm}$ である。
:::
::::

::::exercise{level=basic}
$m,g_C$ を保って密度を四倍にするとプラズマ周波数は何倍になるか。

:::solution
$\omega_p^2$ は密度に比例するため $\sqrt4=2$ 倍になる。
:::
::::

::::exercise{level=advanced}
球の交換体積を $q$ で積分し、式 :eqref[eq-mb-jellium-exchange] を確認せよ。

:::solution
$B(q)=\pi(16k_F^3-12k_F^2q+q^3)/12$ である。$0$ から $2k_F$ の積分は $\pi(32-24+4) k_F^4/12=\pi k_F^4$ である。$-g_C/(4\pi^4)$ を掛け、$n=k_F^3/(3\pi^2)$ で割ると $-3g_Ck_F/(4\pi)$ になる。
:::
::::

::::exercise{level=advanced}
静的 RPA 密度応答の $q\to0$ の先頭項を求めよ。

:::solution
$\chi_0\simeq-D_F$ なので $\chi_{\rm RPA}\simeq-D_F q^2/(q^2+\kappa^2)$ である。$\kappa^2=4\pi g_CD_F$ を入れると $-q^2/(4\pi g_C)$ が先頭である。
:::
::::

::::exercise{level=advanced}
小波数でプラズモンが粒子–正孔励起の外側にあることを示せ。

:::solution
占有球の $\vb k$ から $\vb k+\vb q$ に励起するエネルギーは $\vb k\cdot\vb q/m+q^2/(2m)\le v_Fq+q^2/(2m)$ であり、上限は $q\to0$ で零に行く。プラズモンの周波数は同じ極限で $\omega_p>0$ なので、十分小さい $q$ では上限より高い。
:::
::::

::::exercise{level=advanced}
対数係数に現れる $\int_0^1\int_0^1 xy/(x+y)\dd x\dd y$ を計算せよ。

:::solution
対称性で $y<x$ の積分を二倍する。$y=xt$ で内側は $x^2\int_0^1 t/(1+t)\dd t=x^2(1-\log2)$ となる。最後に $x$ を積分して $2(1-\log2)/3$ を得る。
:::
::::

## 確認テスト

::::quiz
Jellium の平行移動対称な Hartree 項はどうなるか。

- [x] 一様な背景により相殺される。
- [ ] 交換項と同じ値になる。
- [ ] 正負の電荷で加算される。

:::explanation
一様電子密度と背景は零波数で電荷が反対であり、電荷の和が零になる。交換は同スピン電子の反対称化に由来し、背景との相殺の対象ではない。
:::
::::

::::quiz{answer="0.25" tol="1e-6"}
$\kappa=4\,\mathrm{nm^{-1}}$ の長さを $\mathrm{nm}$ 単位で答えよ。

:::explanation
減衰係数の逆数なので $1/4=0.25$ である。
:::
::::

::::quiz{answer="2" tol="1e-6"}
他を固定して密度を四倍にしたときの $\omega_p$ の倍率を答えよ。

:::explanation
二乗が密度に比例するため倍率は $\sqrt4=2$ である。
:::
::::

::::quiz
零温度の三次元電子ガスで Friedel 振動を生むのはどれか。

- [ ] 自由応答の $q=0$ での一定項だけ
- [x] 自由応答の $q=2k_F$ での非解析項
- [ ] 背景電荷の運動

:::explanation
$2k_F$ の対数微分を持つ項を位置へ逆変換すると振動と代数的な減衰を得る。小波数の定数だけを残した近似からは指数遮蔽が得られる。
:::
::::

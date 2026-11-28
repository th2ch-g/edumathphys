---
title: "第二量子化と多体 Hamiltonian"
description: "占有数基底から多体 Hamiltonian を組み立て、Hartree–Fock 近似と Wick 展開を導く。"
notation:
  - symbol: '\mathcal F_\pm'
    meaning: 対称または反対称 Fock 空間
  - symbol: 'a_i^\dagger,a_i'
    meaning: モード i の生成・消滅演算子
  - symbol: '\gamma'
    meaning: 一粒子密度行列
  - symbol: '\eta'
    meaning: Bose 系では +1、Fermi 系では −1
---

## この章の目標

- 占有数基底で生成・消滅演算子の係数と符号を計算する。
- 連続空間と格子上で一体・二体の Hamiltonian を立てる。
- Slater 行列式に対する変分から Hartree–Fock 方程式を導く。
- 縮約を定義し、Wick 展開の各項の符号を求める。

## Fock 空間と生成消滅演算子の復習

多体系で一粒子軌道を取り替えるとき、各粒子に名前をつけて追跡する必要はない。直交軌道を選び、それぞれを何個の粒子が占めるかを指定すればよい。[同種粒子](/physics/quantum-mechanics/identical-particles/) の :ref[def-qm-occupation] を出発点に、粒子の出入りも扱える空間を作る。粒子数部門の直和は :ref[def-aqm-fock]{full} に対応する。

:::definition[Fock 空間]{#def-mb-fock-space}
一粒子 Hilbert 空間を $\mathcal H$ とし、内積は ket 側について線形とする。$\eta=+1$ を Bose 系、$\eta=-1$ を Fermi 系に割り当てる。$\mathcal H_\eta^{(N)}$ を、粒子の任意の交換に対して Bose 系では対称、Fermi 系では反対称なテンソルの空間とする。$\mathcal H_\eta^{(0)}=\C\ket{0}$ とし、Hilbert 直和 $\mathcal F_\eta=\bigoplus_{N=0}^\infty\mathcal H_\eta^{(N)}$ を :term[Fock 空間]{en="Fock space" yomi="ふぉっくくうかん"} と呼ぶ。生成・消滅演算子の積は、まず有限個の占有数基底ベクトルが張る共通の定義域で扱う。
:::

### 占有数基底と交換関係

直交規格化した一粒子基底 $\{u_i\}$ に順序を固定する。積 $u_{i_1}\otimes\cdots\otimes u_{i_N}$ を対称化・反対称化して規格化すると、占有数基底 $\ket{\vb n}$ が得られる。Bose 系では $n_i$ が非負整数、Fermi 系では $0,1$ である。実際、異なる占有数をもつ積の軌道列は一つも一致しないため、その対称化の内積は零である。対称化は射影であり、元の積基底が張る空間を射影すれば像全体を張る。したがってこの基底は完全である。

生成演算子は占有数を一つ増やし、消滅演算子はその随伴である。$\vb e_i$ を第 $i$ 成分だけが $1$ の列とすると、非零成分は
$$
\begin{aligned}
a_i^\dagger\ket{\vb n}&=\sqrt{n_i+1}\ket{\vb n+\vb e_i},&
a_i\ket{\vb n}&=\sqrt{n_i}\ket{\vb n-\vb e_i}
&&(\eta=+1),\\
a_i^\dagger\ket{\vb n}&=(-1)^{\sum_{j<i}n_j}\ket{\vb n+\vb e_i}\quad(n_i=0),&
a_i\ket{\vb n}&=(-1)^{\sum_{j<i}n_j}\ket{\vb n-\vb e_i}\quad(n_i=1)
&&(\eta=-1).
\end{aligned}
$$
Bose の係数は、規格化した対称化積に $u_i$ を挿入して $\sqrt{N+1}$ 倍すると得られる。挿入後の占有数の積の階乗 $\prod_j(n_j+\delta_{ij})!$ と挿入前の積の比は $n_i+1$ である。Fermi 系では、左端に挿入した $u_i$ を所定の位置へ動かす交換の回数が $\sum_{j<i}n_j$ である。

:::proposition[正準交換・反交換関係]{#prop-mb-canonical-algebra}
有限粒子ベクトル上で
$$
a_i a_j^\dagger-\eta a_j^\dagger a_i=\delta_{ij},\qquad
a_i a_j-\eta a_j a_i=0,\qquad
a_i^\dagger a_j^\dagger-\eta a_j^\dagger a_i^\dagger=0.
\label{eq-mb-ccr-car}
$$
:::

:::proof
同じモードの Bose 系では $aa^\dagger$ と $a^\dagger a$ の固有値がそれぞれ $n+1,n$ である。同じモードの Fermi 系ではそれぞれ $1-n,n$ であり、同種の演算子を二度作用させると零になる。異なるモードの Bose 系では一方の作用が他方の係数を変えない。Fermi 系で $i<j$ とすると、$i$ の占有数が変わると $j$ の符号だけが反転する。非零の作用に限れば到達先の占有数はどちらの順でも同じであり、零の作用は二つの順で同時に生じる。これで全ての式が従う。
:::

:::example[Bose 増強と Fermi 符号]
Bose 系で $a_2^\dagger a_1\ket{2,3}=\sqrt 2\sqrt 4\ket{1,4}$ であり、この遷移の振幅の絶対値の二乗は $8$ である。Fermi 系では $a_2 a_1\ket{1,1}=\ket{0,0}$、$a_1a_2\ket{1,1}=-\ket{0,0}$ である。
:::

### 粒子数と場演算子

$\hat N=\sum_i a_i^\dagger a_i$ は $\hat N\ket{\vb n}=(\sum_i n_i)\ket{\vb n}$ を満たす。この固有値を生成前後で比較すると $\comm{\hat N}{a_i^\dagger}=a_i^\dagger$、$\comm{\hat N}{a_i}=-a_i$ である。$e^{i\theta\hat N}$ による共役で生成演算子は $e^{i\theta}$ 倍されるため、生成と消滅が同数の単項式は粒子数を保存する。

位置・スピンをまとめて $x$ と書き、$\psi(x)=\sum_i u_i(x)a_i$ と置く。$\eta$ を用いた交換関係の右辺は完全性から $\sum_i u_i(x)u_i(y)^*=\delta(x-y)$ となる。連続空間ではこの関係を試験関数で積分した演算子の等式として読む。例えば体積 $V$ の箱の平面波一個に対して $\expval{\psi^\dagger(x)\psi(x)}=1/V$ であり、箱内で積分すると粒子数 $1$ を得る。

### 有限基底での状態数と基底変換

$M$ モードに $N$ 個の Bose 粒子を置くと、$N$ 個の粒子記号の間または端に $M-1$ 本の区切りを入れた列と占有数が対応する。従って次元は $\binom{M+N-1}{N}$ である。Fermi 系では占有する $N$ モードを選ぶので次元は $\binom MN$、全 Fock 空間の次元は $2^M$ になる。

:::example[同じ粒子数での次元]
$M=4,N=3$ なら Bose 部門は $\binom63=20$ 次元、Fermi 部門は $\binom43=4$ 次元である。Fermi の $N=0,\ldots,4$ の次元を足すと $1+4+6+4+1=16$ となる。Hamiltonian が粒子数を保存するときは、全空間を一度に対角化せず各部門へ制限できる。
:::

軌道を $w_j=\sum_i u_iU_{ij}$ へ変えると、$b_j^\dagger=\sum_i U_{ij}a_i^\dagger$、$b_j=\sum_i U_{ij}^*a_i$ である。:eqref[eq-mb-ccr-car] の右辺は $\sum_i U_{ij}^*U_{il}=\delta_{jl}$ となるので交換関係を保つ。例えば $w_\pm=(u_1\pm u_2)/\sqrt2$ の二つを占める Fermi 状態は $b_+^\dagger b_-^\dagger\ket0=-a_1^\dagger a_2^\dagger\ket0$ であり、変換行列の行列式だけが掛かる。一般の $N$ 個の占有軌道について、積を展開した各係数は対応する $N$ 行を選んだ小行列式である。置換の符号は生成演算子を元の順序に戻す反交換から出る。

### 粒子を区別しない空間と占有数の規格化

一粒子空間を有限次元の $\mathcal H_1$ とし、その直交規格化基底を $(u_i)_{i=1}^M$ とする。添字付きのテンソル積では、例えば $u_1\otimes u_2$ と $u_2\otimes u_1$ は直交する。しかし、同種粒子の交換では観測量が変わらない。従って Bose と Fermi の $N$ 粒子空間はそれぞれ対称化と反対称化で取り出す。

置換 $\pi\in S_N$ に対応するユニタリを $U_\pi$、符号を $\chi_+(\pi)=1$、$\chi_-(\pi)=\sgn(\pi)$ とし、$P_\pm=(N!)^{-1}\sum_\pi \chi_\pm(\pi)U_\pi$ と置く。$U_\pi^\dagger=U_{\pi^{-1}}$ かつ $\chi_\pm(\pi^{-1})=\chi_\pm(\pi)$ なので $P_\pm^\dagger=P_\pm$ である。二乗では各積 $\rho=\pi\sigma$ が $N!$ 回現れ、符号は常に $\chi_\pm(\rho)$ である。従って $P_\pm^2=P_\pm$ である。また $U_\sigma P_\pm=\chi_\pm(\sigma)P_\pm$ であり、逆にこの交換条件を満たすベクトルには $P_\pm$ が恒等的に作用する。

Bose の占有数 $(n_1,\ldots,n_M)$ に対応する一つの積ベクトルを $u_{\vb n}$ とする。$P_+u_{\vb n}$ のノルム二乗は $\inner{u_{\vb n}}{P_+u_{\vb n}}$ であり、内積が $1$ になる置換は同じ軌道の中の $\prod_i n_i!$ 個である。規格化状態は $\sqrt{N!/\prod_i n_i!}\,P_+u_{\vb n}$ である。Fermi の場合、重複軌道があればその交換で項が打ち消し合う。重複のない場合に内積へ残る置換は恒等置換だけなので、規格化には $\sqrt{N!}$ を掛ける。占有数の異なる二つの積は、どの置換を入れても直交する。従ってこれらの規格化状態は正規直交基底である。

:::example[同じ軌道に置くと規格化が変わる]
二つの軌道 $u,v$ に二個の Bose 粒子を置く。$\ket{2,0}$ は $u\otimes u$ であり、$\ket{1,1}$ は $(u\otimes v+v\otimes u)/\sqrt2$ である。例えば $\phi=(u+v)/\sqrt2$ の積状態をこの基底で展開すると
$$
\phi\otimes\phi=\frac12\ket{2,0}+\frac1{\sqrt2}\ket{1,1}+\frac12\ket{0,2}
$$
となる。従って各軌道に一個ずつ見つかる確率は $1/2$ であり、測定前の状態は占有数の固有状態を一つ選んだものではない。
:::

### Fock 演算子の定義域

有限個の Bose 軌道でも全 Fock 空間は無限次元である。例えば一つの軌道の $a^\dagger$ は $\ket n$ のノルムを $\sqrt{n+1}$ 倍するから有界ではない。演算子積を計算するときは、まず有限個の占有数ベクトルの線形結合に制限する。この部分空間は生成・消滅で閉じており、:eqref[eq-mb-ccr-car] の両辺をそこで評価できる。一般のノルム収束列に関係式を移す場合は、粒子数のモーメントなどを用いて定義域を確かめる必要がある。例えば一モード状態 $\sum_{n\ge0}d_n\ket n$ について $\|a^\dagger\psi\|^2=\sum_{n\ge0}(n+1)|d_n|^2$ であり、この和が有限であることが生成演算子の自然な定義域を与える。

### Fock 空間上で基底変換を実装する

一粒子空間のユニタリ $U$ を全粒子空間へ持ち上げる。$N$ 粒子空間では $U^{\otimes N}$ とし、真空では $1$ とする。全ての粒子に同じ変換を施すため $U^{\otimes N}$ は各 $U_\pi$ と可換であり、対称・反対称部分空間を保つ。各部門の直和を $\mathscr U$ とすると、直和のノルムは部門のノルム二乗の和だから、$\mathscr U$ は全 Fock 空間でユニタリである。

規格化した積への生成演算子の作用を使うと、$a^\dagger(u)$ は $N$ 粒子状態に対し $\sqrt{N+1}$ を掛けてから $u$ を添え、対称化または反対称化する操作である。従って
$$
\mathscr U a^\dagger(u)\mathscr U^\dagger=a^\dagger(Uu),\qquad
\mathscr U a(u)\mathscr U^\dagger=a(Uu)
$$
となる。消滅の等式は随伴を取ったものである。$U(\lambda)=e^{-i\lambda h}$ が一体 Hermite 行列 $h$ で生成されるなら、$N$ 粒子で $\partial_\lambda U(\lambda)^{\otimes N}|_0=-i\sum_{j=1}^Nh_j$ である。従ってその Fock 空間の生成子が $\widehat h$ である。これは一体演算子の持ち上げを、変換した生成演算子の作用から確認する方法でもある。Bose の無限直和で生成子は有界とは限らないが、有限粒子の部分空間上では微分を部門ごとに行える。

:::example[二つのモードを混合したときの粒子数測定]
$u_1,u_2$ の二軌道で $Uu_1=(\sqrt3 u_1+u_2)/2$ とし、最初は $u_1$ に Bose 粒子三個を置く。変換後は
$$
\mathscr U\ket{3,0}
=\sum_{j=0}^3\sqrt{\binom 3j}
\left(\frac{\sqrt3}{2}\right)^{3-j}\left(\frac12\right)^j\ket{3-j,j}.
$$
二項展開の係数に $(a_1^\dagger)^{3-j}(a_2^\dagger)^j\ket0=\sqrt{(3-j)!j!}\ket{3-j,j}$ を入れた結果である。従って軌道 $2$ の個数 $j=0,1,2,3$ の確率は $27/64,27/64,9/64,1/64$ であり、平均は $3/4$、分散は $9/16$ である。粒子が一つずつ決まった軌道を選んだ混合状態を仮定しなくても、この占有数測定には二項分布が現れる。
:::

:::example[二電子を入れた空間の位相]
二軌道を両方占めた Fermi 状態に同じ一粒子の $U$ を作用させる。$(Uu_1)\wedge(Uu_2)=(U_{11}U_{22}-U_{21}U_{12})u_1\wedge u_2$ であるため、二軌道を混ぜても占有は変わらず、全体の位相 $\det U$ だけが残る。例えば実の回転角 $\pi/6$ では $\det U=1$ であり、軌道を交換する行列では $-1$ になる。これは同じ二電子状態の基底変換であり、一粒子を別の準位へ散乱した結果と区別する。
:::

## 一体演算子と二体演算子

### 一体演算子の持ち上げ

一粒子演算子 $h$ が各粒子に作用する場合の和 $\sum_{p=1}^N h(p)$ を考える。$h_{ij}=\mel{u_i}{h}{u_j}$ とすると、その Fock 空間での表示は
$$
\hat h=\sum_{ij}h_{ij}a_i^\dagger a_j .
\label{eq-mb-one-body}
$$
これを確かめるには、積の原子となる $\ket{j_1,\dots,j_N}=a_{j_1}^\dagger\cdots a_{j_N}^\dagger\ket{0}$ に作用させればよい。:eqref[eq-mb-ccr-car] により、$a_j$ は第 $p$ 因子を $\eta^{p-1}\delta_{j j_p}$ とともに除く。左から $a_i^\dagger$ を掛けて第 $p$ 位置に戻すと、同じ $\eta^{p-1}$ が掛かる。積が $1$ になるため、$h$ が一粒子を置換する項の和を得る。

### 二体演算子と反対称化

交換対称な二粒子演算子 $v$ を用い、$v_{ij;kl}=\mel{u_i\otimes u_j}{v}{u_k\otimes u_l}$ とする。$\sum_{p<q}v(p,q)$ の表示は
$$
\hat V=\tfrac12\sum_{ijkl}v_{ij;kl}a_i^\dagger a_j^\dagger a_l a_k .
\label{eq-mb-two-body}
$$
二度の消滅を順に追うと符号を決められる。$a_k$ が第 $p$ 因子、$a_l$ が第 $q$ 因子を除くとき、除去の符号は二つの軌道を順に左端へ移す置換の符号に等しい。その逆置換で $i,j$ を元の位置へ挿入するので、符号の積は $1$ である。和は順序付きの粒子対を数えるため係数 $1/2$ が要る。Fermi 系では $\bar v_{ij;kl}=v_{ij;kl}-v_{ij;lk}$ と定めると、$a_l a_k=-a_k a_l$ により係数 $1/4$ と $\bar v$ を用いる表式に移れる。

位置表示で局所的な $v(x,y)=v(y,x)$ なら $\hat V=\tfrac12\int\dd x\dd y\,v(x,y)\psi^\dagger(x)\psi^\dagger(y)\psi(y)\psi(x)$ となる。$x$ にスピンを含めた積分はスピン和も表す。密度の積に書き換えると交換関係のデルタ項が生じるため、自己相互作用を除くには正規順序が必要である。

:::example[Bose の接触相互作用]
一つの軌道へ射影した相互作用を $\tfrac U2 a^\dagger a^\dagger aa$ とする。二度の消滅で係数 $\sqrt{n(n-1)}$ が出るので、エネルギーは $U n(n-1)/2$ である。$U=2\,\mathrm{meV}$、$n=4$ なら $12\,\mathrm{meV}$ であり、$n=1$ では零である。
:::

:::example[二電子の直接項と交換項]
異なる二軌道 $1,2$ を占めると $\expval{\hat V}=v_{12;12}-v_{12;21}$ である。例えば二つの積分が $5$ と $2\,\mathrm{eV}$ なら期待値は $3\,\mathrm{eV}$ である。$v$ がスピンに作用せず二軌道のスピンが直交するとき、交換積分には互いに直交するスピンの内積が入り零になる。
:::

### 一体密度行列の許される範囲

Fermi 粒子の一体密度行列 $\gamma_{ij}=\langle a_j^\dagger a_i\rangle$ を考える。任意の列 $u$ に対して $a(u)=\sum_i u_i^*a_i$ とすれば $u^\dagger\gamma u=\langle a(u)^\dagger a(u)\rangle\ge0$ である。反交換関係から $\|u\|^2-u^\dagger\gamma u=\langle a(u)a(u)^\dagger\rangle\ge0$ でもある。よって固有値は $0$ から $1$ に収まる。一方、Bose の可換関係では後者の式が和になるため、固有値の上限には粒子数を指定する必要がある。$N$ 粒子状態では両統計で $\Tr\gamma=N$ が成立する。

$N$ 粒子の純粋 Fermi 状態では、$\gamma^2=\gamma$ と単一の Slater 行列式であることが同値である。実際、Slater 行列式の軌道基底では $\gamma$ の対角要素は占有数 $0,1$ である。逆に $\gamma$ が射影なら、これを対角化する軌道基底で各 $n_i$ の分散は $p_i-p_i^2=0$ である。従って状態は各 $n_i$ の固有状態であり、指定された占有ビット列に対応する空間は一次元だから単一の行列式になる。

:::example[四軌道での純粋相関状態]
二電子状態 $(\ket{1,1,0,0}+\ket{0,0,1,1})/\sqrt2$ を考える。一体演算子 $a_j^\dagger a_i$ は高々一つの占有軌道を別の軌道へ動かせる。二つの成分には占有軌道が共通しないため、交差する行列要素は消える。従って $\gamma=I_4/2$ であり、$\gamma^2=I_4/4$ となる。位相をもつ純粋状態でも、占有が二分されると一体密度だけでは Slater 行列式へ戻せない。
:::

### Slater–Condon 則の添字の追跡

占有集合 $I$ から $J$ へ移るとき、二体演算子は二個の占有を除いて二個を加えるので、$|I\setminus J|>2$ なら行列要素は零である。ちょうど二つ $i,j$ を $a,b$ に移す場合は、行列式の位相を $\ket J=a_a^\dagger a_b^\dagger a_j a_i\ket I$ と選べる。消滅添字と生成添字を入れ替えた四つの項に対し、演算子の交換符号を行列要素へ移すと $(v_{ab;ij}-v_{ab;ji}-v_{ba;ij}+v_{ba;ji})/2=\bar v_{ab;ij}$ になる。最後の等式では粒子交換の対称性 $v_{ba;ji}=v_{ab;ij}$ を使った。

一つの軌道だけを $i$ から $a$ へ移す場合は、$\ket J=a_a^\dagger a_i\ket I$ と位相を選ぶ。二体項では片方の電子が共通の占有軌道 $j$ に戻る。その各項にも同じ四通りの添字順序があり、残る行列要素は $h_{ai}+\sum_{j\in I\setminus\{i\}}\bar v_{aj;ij}$ である。対角行列要素では二つの占有軌道の組が順序付き和に二度現れるため、$E=\sum_i h_{ii}+\frac12\sum_{ij}\bar v_{ij;ij}$ を得る。この三つの式により、行列全体を作らず、参照配置と結合する配置を特定できる。

:::example[二重励起を一体ポテンシャルで結べない理由]
占有集合 $\{1,2\}$ と $\{3,4\}$ を取る。例えば $a_3^\dagger a_1$ を最初の配置に作用させると、軌道 $2$ が残るので重なりは零である。これに対し $a_3^\dagger a_4^\dagger a_2 a_1$ は順序まで合った $\ket{3,4}$ を生成する。$v_{34;12}=3\,\mathrm{meV}$、$v_{34;21}=1\,\mathrm{meV}$ のとき、二体行列要素は $2\,\mathrm{meV}$ である。
:::

### 連続基底・離散基底間での二体行列要素

体積 $V$ の周期箱でスピンを省略し、平面波を $u_k(\vb r)=V^{-1/2}e^{i\vb k\cdot\vb r}$ とする。相対位置にだけ依存する作用を $v(\vb r-\vb r')$ として、周期関数の Fourier 係数を $v_q=\int_V v(\vb r)e^{-i\vb q\cdot\vb r}\dd^dr$ と置く。一体軌道から作る直積での行列要素は
$$
\begin{aligned}
v_{k_1 k_2;k_3 k_4}
&=V^{-2}\int_V\dd^dr\,\dd^dr'\,
e^{-i(\vb k_1-\vb k_3)\cdot\vb r-i(\vb k_2-\vb k_4)\cdot\vb r'}v(\vb r-\vb r')\\
&=V^{-1}\delta_{\vb k_1+\vb k_2,\vb k_3+\vb k_4}\,v_{\vb k_1-\vb k_3}.
\end{aligned}
$$
ここで $\vb s=\vb r-\vb r'$ で相対座標を積分し、$\vb r'$ の積分で全運動量の Kronecker のデルタを得た。スピンに依存しない相互作用なら直接項には $\delta_{\sigma_1\sigma_3}\delta_{\sigma_2\sigma_4}$ が掛かり、交換項には $\delta_{\sigma_1\sigma_4}\delta_{\sigma_2\sigma_3}$ が掛かる。分母の $V$ は $v_q$ の単位がエネルギーかける体積であるために必要である。

:::example[接触結合の同スピン消失]
$v(\vb s)=g\delta(\vb s)$ の行列要素は運動量保存の下で $g/V$ である。同じスピンの二電子に対して直接と交換はどちらも $g/V$ なので差は零となる。位置表示でも二電子の反対称な波動関数は $\vb r=\vb r'$ で零になり、接触で行列要素が消える。一方、$\uparrow,\downarrow$ をそれぞれ同じスピンの終状態へ動かす直接過程には $g/V$ が残る。この議論は接触頂点の行列要素についてであり、三次元の散乱へ適用するときは紫外の規格化を別に行う。
:::

位置の二体密度に関しても同じ直接・交換が現れる。まずスピンなしの Slater 状態で $\hat n(\vb r)=\psi^\dagger(\vb r)\psi(\vb r)$、$\gamma(\vb r,\vb r')=\langle\psi^\dagger(\vb r')\psi(\vb r)\rangle$ とする。正規順序の二点密度 $P_2(\vb r,\vb r')=\langle\psi^\dagger(\vb r)\psi^\dagger(\vb r')\psi(\vb r')\psi(\vb r)\rangle$ は
$$
P_2(\vb r,\vb r')=n(\vb r)n(\vb r')-|\gamma(\vb r,\vb r')|^2
$$
である。通常の積 $\langle \hat n(\vb r)\hat n(\vb r')\rangle$ では $\psi(\vb r)$ と $\psi^\dagger(\vb r')$ の反交換で得る $\delta(\vb r-\vb r')n(\vb r)$ も加わる。二体作用に入れるときは同じ粒子を数えない $P_2$ を用いる。射影性から $\int|\gamma(\vb r,\vb r')|^2\dd^dr'=\gamma(\vb r,\vb r)$ であるため、$\int P_2(\vb r,\vb r')\dd^dr'=(N-1)n(\vb r)$ となる。固定した場所の電子の周りに他の $N-1$ 個を数える関係であり、交換の負の項を捨てれば $N$ 個を数えてしまう。

スピンを持つ粒子の数密度では $n(\vb r)=\sum_\sigma\gamma_{\sigma\sigma}(\vb r,\vb r)$ を用いる。各スピン対に上の式を適用して和を取ると、直接項は $n(\vb r)n(\vb r')$、交換項は $-\sum_{\sigma\sigma'}|\gamma_{\sigma\sigma'}(\vb r,\vb r')|^2$ となる。例えばスピンごとに軌道を占めた状態では $\gamma_{\uparrow\downarrow}=0$ なので、異スピン間の交換孔はない。

:::example[三軌道での交換による同時占有の低下]
行列式が $u_1$ と $(u_2+u_3)/\sqrt2$ を占めるとする。軌道 $2,3$ のそれぞれの占有は $1/2$ だが、そこにいるのは同じ一電子なので同時占有は零である。軌道 $2,3$ に制限した密度行列は全成分 $1/2$ なので、対の行列式 $1/4-1/4=0$ となる。これに対し $1,2$ の成分は対角 $1,1/2$、非対角零だから同時占有は $1/2$ である。固定粒子数の関係として $\langle n_2n_1\rangle+\langle n_2n_3\rangle=(2-1)\langle n_2\rangle$ が成り立つ。
:::

## 格子模型（強束縛模型と Hubbard 模型）

### 格子 Fourier 変換

各格子点 $j$ に局在軌道を一つ選び、電子のスピンを $\sigma$ とする。局在軌道の間の遷移を一体項に取る方法を :term[強束縛近似]{en="tight-binding approximation" yomi="きょうそくばくきんじ"} と呼ぶ。格子間隔 $a$、$L$ 点の周期鎖の最近接模型は
$$
H_t=-t\sum_{j,\sigma}(c^\dagger_{j+1,\sigma}c_{j,\sigma}+
c^\dagger_{j,\sigma}c_{j+1,\sigma}),\quad
c_{j,\sigma}=L^{-1/2}\sum_k e^{ikja}c_{k,\sigma}.
$$
ここで $k=2\pi m/(La)$ を一つの Brillouin 域から選ぶ。サイト和が $L\delta_{kk'}$ であることは等比級数で示せる。これを代入すると $H_t=\sum_{k\sigma}\epsilon(k)n_{k\sigma}$、$\epsilon(k)=-2t\cos(ka)$ となる。状態の時間依存は $e^{-i\epsilon t_{\rm time}/\hbar}$ であり、速度は $\hbar^{-1}\partial_k\epsilon$ である。$t=1\,\mathrm{eV}$、$a=0.3\,\mathrm{nm}$ なら最大速度 $2ta/\hbar$ はおよそ $9.1\times10^5\,\mathrm{m\,s^{-1}}$ である。

### Hubbard 相互作用

同じサイトの二重占有にだけエネルギー $U$ を課す :term[Hubbard 模型]{en="Hubbard model" yomi="はばーどもけい"} は
$$
H=H_t+U\sum_j n_{j\uparrow}n_{j\downarrow}.
\label{eq-mb-hubbard}
$$
軌道を一つに絞ると同スピンの二重占有は反交換関係で消える。Fourier 変換では、四つの演算子が与える $1/L^2$ とサイト和の $L$ により係数は $U/L$ となる。具体的には
$$
H_U=\frac UL\sum_{k,p,q}
c^\dagger_{k+q,\uparrow}c^\dagger_{p-q,\downarrow}c_{p,\downarrow}c_{k,\uparrow}.
$$
総結晶運動量は逆格子ベクトルを法として保存する。単純な一粒子分散をもつ系でも各 $n_{k\sigma}$ は、$q\ne0$ の散乱により一般には保存されない。

:::example[原子極限と四点鎖]
一サイトで $t=0$ とし、化学ポテンシャルを $\mu$ とすると、空、上スピン、下スピン、二重占有の $H-\mu\hat N$ の固有値は $0,-\mu,-\mu,U-2\mu$ である。$U=4\,\mathrm{eV}$、$\mu=2\,\mathrm{eV}$ なら一粒子状態が最低である。一方、$U=0$ の四点鎖では、許容波数に対する一スピンの固有値は $-2t,0,2t,0$ である。
:::

:::example[ホップと相互作用の単位]
$t=0.25\,\mathrm{eV}$ の周期鎖の帯域幅は $4t=1\,\mathrm{eV}$ である。$U=4\,\mathrm{eV}$ を入れると局所二重占有の費用は帯域幅の四倍となる。原子極限の計算から始める場合、小さい比には $t/U=1/16$ を用いる。
:::

### Peierls 位相と格子上の流れ

電荷を $Q$ とし、サイト $j$ から $i$ へ向く結合に $\theta_{ij}=Q\int_j^i\vb A\cdot\dd\vb r/\hbar$ を付ける。変換 $\vb A\mapsto\vb A+\grad\lambda$ は $\theta_{ij}$ に $Q(\lambda_i-\lambda_j)/\hbar$ を加える。演算子を $c_i\mapsto e^{iQ\lambda_i/\hbar}c_i$ と再位相化すると、$e^{i\theta_{ij}}c_i^\dagger c_j$ の位相は相殺される。従って枝の位相は基底に依存し、閉路上の符号付き和だけが磁束を記録する。

一つの結合の Hamiltonian を $H_{ij}=-t(e^{i\theta_{ij}}c_i^\dagger c_j+e^{-i\theta_{ij}}c_j^\dagger c_i)$ とする。$[n_i,c_i^\dagger c_j]=c_i^\dagger c_j$ と $[n_i,c_j^\dagger c_i]=-c_j^\dagger c_i$ を使うと、結合から $i$ への粒子流は
$$
\dot n_i\big|_{ij}=\frac{it}{\hbar}(e^{i\theta_{ij}}c_i^\dagger c_j-e^{-i\theta_{ij}}c_j^\dagger c_i)
=-\frac1\hbar\partial_{\theta_{ij}}H_{ij}.
$$
$\dot n_j$ は逆符号になるため、全サイトを足すと結合ごとに消える。電流を得るときにはこれに $Q$ を掛ける。局所 Hubbard 項は各 $n_i$ と可換であり、この流れの演算子に直接は寄与しない。

:::example[四サイト環を貫く半磁束量子]
$L=4$ の環で $i=j+1$ 向きの位相を一様に $\theta=\pi/4$ とする。平面波係数 $e^{ikj}$ を作用させると固有値は $-2t\cos(k-\theta)$、$k=0,\pi/2,\pi,3\pi/2$ である。順に $-\sqrt2t,-\sqrt2t,+\sqrt2t,+\sqrt2t$ となり、磁束が零での $-2t,0,+2t,0$ から変化する。隣接サイトの再位相化で三本から位相を除いても、最後の一本の位相に $\pi$ が残るので固有値は同じである。
:::

### 粒子・正孔変換と半充填

二部格子を $A,B$ とし、$\varepsilon_i=\pm1$ を部分格子で反転させる。両スピンについて $c_{i\sigma}\mapsto\varepsilon_i c_{i\sigma}^\dagger$ と置くと、反交換関係は保たれ、$n_{i\sigma}\mapsto1-n_{i\sigma}$ となる。実のホッピングで $i,j$ が異なる部分格子なら $c_i^\dagger c_j\mapsto-\varepsilon_i\varepsilon_j c_j^\dagger c_i=c_j^\dagger c_i$ であり、両向きの和は不変である。サイトの相互作用と化学ポテンシャルは
$$
Un_{i\uparrow}n_{i\downarrow}-\mu(n_{i\uparrow}+n_{i\downarrow})
\mapsto Un_{i\uparrow}n_{i\downarrow}-(U-\mu)(n_{i\uparrow}+n_{i\downarrow})+U-2\mu .
$$
従って $\mu=U/2$ で大正準 Hamiltonian が対称であり、有限系の熱状態では $\langle n_{i\uparrow}+n_{i\downarrow}\rangle=1$ である。基底縮退がある場合は対称な混合状態に同じ議論を適用する。第二近接ホッピングには同じ部分格子が入るので、この対称性から半充填を決めることはできない。

### 二サイト Hubbard 模型の対称性で行列を小さくする

二サイトに電子を二つ置き、ホッピングを $-t$、同サイトの反発を $U\ge0$ とする。スピンの基底を上・下とし、軌道順序を $(1\uparrow,1\downarrow,2\uparrow,2\downarrow)$ に固定する。二つの同サイト対を $\ket{D_i}=c_{i\uparrow}^\dagger c_{i\downarrow}^\dagger\ket0$ とし、別サイト一重項を
$$
\ket S=(c_{1\uparrow}^\dagger c_{2\downarrow}^\dagger-c_{1\downarrow}^\dagger c_{2\uparrow}^\dagger)\ket0/\sqrt2
$$
とする。$c_{2\uparrow}^\dagger c_{1\uparrow}\ket{D_1}=-c_{1\downarrow}^\dagger c_{2\uparrow}^\dagger\ket0$、$c_{2\downarrow}^\dagger c_{1\downarrow}\ket{D_1}=c_{1\uparrow}^\dagger c_{2\downarrow}^\dagger\ket0$ である。従って $H_t\ket{D_1}=-\sqrt2t\ket S$ となり、$\ket{D_2}$ にも同じ値が出る。$\ket{D_\pm}=(\ket{D_1}\pm\ket{D_2})/\sqrt2$ とすると、$D_-$ は作用されず固有値 $U$ を持つ。$S,D_+$ の二次元ブロックは
$$
\begin{pmatrix}0&-2t\\-2t&U\end{pmatrix}.\label{eq-mb-hubbard-dimer-singlet}
$$
三つの三重項は $S^z$ が最大の状態で二重占有が禁じられ、スピン回転で残りも同じ固有値 $0$ を持つ。従って $t\ne0$ では下の固有値 $E_-=(U-\sqrt{U^2+16t^2})/2$ が基底である。$U\gg t$ で $E_-=-4t^2/U+O(t^4/U^3)$ となる。ホップを一度入れた中間の二重占有に $U$ の費用があるため、直接ホップが消える部門にも二次の分裂が残る。

:::example[二重占有の重みを二通りで求める]
$t=1\,\mathrm{meV}$、$U=3\,\mathrm{meV}$ なら固有値は $-1,4\,\mathrm{meV}$ である。基底の固有ベクトルは $(2,1)^\top/\sqrt5$ であり、合計の二重占有の確率は $1/5$ である。別の計算では $\partial_U E_-=(1-U/\sqrt{U^2+16t^2})/2=1/5$ になる。固有値方程式を $U$ で微分して基底の左ベクトルを掛けると $\partial_U H=n_{1\uparrow}n_{1\downarrow}+n_{2\uparrow}n_{2\downarrow}$ の期待値だけが残るので、両手順は一致する。各サイトの二重占有確率は交換対称性から $1/10$ である。
:::

## 電子-フォノン系の Hamiltonian

### 格子振動の量子化

質量 $M$ のイオンが $L$ 個並ぶ鎖で縦変位を $u_j$、共役運動量を $p_j$ とする。調和近似の Hamiltonian は $\sum_j[p_j^2/(2M)+\kappa(u_{j+1}-u_j)^2/2]$ である。格子 Fourier 変換で差が $e^{iqa}-1$ 倍になるため、固有角振動数は $\omega_q=2\sqrt{\kappa/M}|\sin(qa/2)|$ となる。$q=0$ の重心運動を固定し、残りを Bose 粒子の生成消滅演算子で表す。
$$
u_j=\frac1{\sqrt L}\sum_{q\ne0}\frac{e^{iqja}}{\sqrt{2M\omega_q}}(b_q+b_{-q}^\dagger),\qquad
H_{\rm ph}=\sum_{q\ne0}\omega_q(b_q^\dagger b_q+\tfrac12).
$$
$p_j$ は $\dot u_j=p_j/M$ で決める。すなわち各 $q$ の係数を $-i\sqrt{M\omega_q/2}(b_q-b^\dagger_{-q})$ とすれば、交換関係は $[u_j,p_l]=i(\delta_{jl}-1/L)$ となる。$H_{\rm ph}$ への代入では、運動エネルギーとばねの項から $b_qb_{-q}$ とその随伴が反対符号で出て打ち消し合う。SI の変位係数には分子に $\sqrt{\hbar}$ が入り、各量子のエネルギーは $\hbar\omega_q$ となる。

### 線形結合

電子が局所的なひずみを感じ、イオンの変位に比例してサイトエネルギーが変わると仮定する。$H_{\rm ep}=D\sum_j (u_{j+1}-u_j)n_j/a$ なら
$$
H_{\rm ep}=\frac1{\sqrt L}\sum_{kq\sigma}g_q c^\dagger_{k+q,\sigma}c_{k,\sigma}(b_q+b^\dagger_{-q}),\qquad
g_q=\frac{D(e^{iqa}-1)}{a\sqrt{2M\omega_q}} .
\label{eq-mb-ep-lattice}
$$
$D$ はエネルギー、$g_q$ もエネルギーである。$g_{-q}=g_q^*$ により随伴を取っても同じ式になる。$b_q$ の項では電子が $q$ を得て格子が $q$ を失い、$b_{-q}^\dagger$ の項では格子が $-q$ を得る。音響分枝の小さい $q$ では分子が $q$、分母が $|q|^{1/2}$ なので $|g_q|\propto |q|^{1/2}$ である。光学フォノンの局所変位に比例する模型では $\omega_q=\omega_0$ と $g_q=g$ を取ることもできる。

:::example[原子の変位と結合の桁]
$M=40$ 原子質量単位、$\hbar\omega_0=20\,\mathrm{meV}$ の独立振動子なら、零点変位 $\sqrt{\hbar/(2M\omega_0)}$ はおよそ $0.0051\,\mathrm{nm}$ である。サイトエネルギーの傾きが $1\,\mathrm{eV\,nm^{-1}}$ なら、その零点変位に対応する結合は約 $5.1\,\mathrm{meV}$ である。
:::

:::example[フォノン放出の行列要素]
一電子が波数 $k$ にあり、波数 $-q$ のフォノンが三個あるとする。波数 $k+q$ への遷移と同時にフォノンを一個作る振幅は $2g_q/\sqrt L$ である。その過程がエネルギーを保存するには $\epsilon(k)-\epsilon(k+q)=\omega_{-q}$ が必要である。
:::

### 並進が結合しない条件

イオン $\ell$ の平衡位置を $\vb R_\ell$、変位を $\vb u_\ell$ とする。電子とイオンの相対位置だけに依存する $V(\vb r-\vb R_\ell-\vb u_\ell)$ を一次で展開すると、摂動は $-\sum_\ell\vb u_\ell\cdot\grad V(\vb r-\vb R_\ell)$ である。全イオンを同じ $\vb u$ だけ動かした場合、電子 Hamiltonian は座標の並進演算子 $e^{-i\vb u\cdot\vb p/\hbar}$ で共役されたものになる。一次展開すると $\delta h=-(i/\hbar)[\vb u\cdot\vb p,h]$ である。従って固有状態の対角行列要素は零であり、固有エネルギーは一様な平行移動では変わらない。変形ポテンシャルではこの事実を満たすよう、音響モードに対して変位そのものより $\divg\vb u$ を用いる。

一本の縦音響分枝について密度への結合を $D\divg\vb u$ とし、系の体積を $V$、質量密度を $\rho$ とする。フォノン一個での変位の規格化は $\sqrt{\hbar/(2\rho V\omega_{\vb q})}$ であり、勾配の Fourier 因子から
$$
g_{\vb q}=iD(\vb q\cdot\vb e_{\vb q})\sqrt{\frac{\hbar}{2\rho V\omega_{\vb q}}}
$$
を得る。縦偏極 $\vb e_{\vb q}$ の位相を $\vb e_{-\vb q}=\vb e_{\vb q}^*$ と取る。勾配の $\vb q$ が反転するため $g_{-\vb q}=g_{\vb q}^*$ が成り立つ。$\omega_{\vb q}=c_s|\vb q|$ なら $|g_{\vb q}|^2$ は $|\vb q|$ に比例する。$D$ はエネルギー、平方根は長さ、$|\vb q|$ は長さの逆数なので、$g_{\vb q}$ はエネルギーの次元になる。

:::example[音響波の波長を変えたときの行列要素]
線形分散の範囲で波長 $60\,\mathrm{nm}$ を $15\,\mathrm{nm}$ に変え、質量密度と箱は固定する。波数は四倍、変位の零点振幅は半分になる。従って結合の大きさは二倍、黄金律で使う行列要素の二乗は四倍になる。最終状態の密度も変わり得るため、この四倍は結合の部分に対する比較である。
:::

### 固定した電荷と緩和後の電荷

Holstein 型の局所モードで $\omega_0>0$ とし、固定した電子占有 $n$ に対し $B=b+gn/\omega_0$ と置く。$H=\omega_0 B^\dagger B-g^2n^2/\omega_0$ なので、振動の基底を選んだ各電子配置は $-g^2n^2/\omega_0$ だけ緩和する。スピン $1/2$ では $n^2=n+2n_\uparrow n_\downarrow$ であり、一体項の変化と局所二体項 $-2g^2n_\uparrow n_\downarrow/\omega_0$ を区別できる。局所振動子と電荷の同時対角化はホップが零のとき厳密である。ホップがある場合、電荷配置の変更に付随する振動子の重なりが残る。

:::example[局所反発の符号が変わる条件]
$\omega_0=40\,\mathrm{meV}$、$g=20\,\mathrm{meV}$ とする。電荷の二次の差 $E_2+E_0-2E_1$ への寄与は $-2g^2/\omega_0=-20\,\mathrm{meV}$ である。例えば裸の $U=15\,\mathrm{meV}$ なら緩和後の差は $-5\,\mathrm{meV}$ となる。ここで比較したのは振動子が各占有に対して基底まで緩和したときのエネルギーである。
:::

### フォノンの吸収・放出と運動量保存

$H_{\rm ep}=\sum_{\vb k,\vb q}g_{\vb q}c_{\vb k+\vb q}^\dagger c_{\vb k}(b_{\vb q}+b_{-\vb q}^\dagger)$ とし、箱の規格化は $g$ に含める。全運動量 $\vb P=\sum_{\vb k}\hbar\vb k\,n_{\vb k}+\sum_{\vb q}\hbar\vb q\,b_{\vb q}^\dagger b_{\vb q}$ の電子部分は $c_{\vb k+\vb q}^\dagger c_{\vb k}$ に対し $+\hbar\vb q$ の交換子を与える。振動部分は $b_{\vb q}$ に対しても $b_{-\vb q}^\dagger$ に対しても $-\hbar\vb q$ を与える。従って両過程で全運動量が保たれる。格子では $\vb k+\vb q$ を第一 Brillouin ゾーンに折り返せるため、保存は逆格子ベクトルを法として読む。$H_{\rm ep}$ の随伴で $\vb q\mapsto-\vb q$、$\vb k\mapsto\vb k+\vb q$ と添字を動かすと、Hermite 性には $g_{-\vb q}=g_{\vb q}^*$ が必要だと分かる。

:::example[吸収と放出で終状態が異なる]
初期の電子波数を $k=0$ とし、$q=\pi/(3a)$ の頂点を取る。$b_q$ を含む項では $q$ のフォノンを除き、電子は $k=q$ へ変わる。$b_{-q}^\dagger$ を含む項では、同じ電子の終状態へ動く一方でフォノンを $-q$ に作る。初期に両フォノンが空なら吸収振幅は零、放出に対するフォノンの行列要素は $1$ である。初期に $q$ に三個、$-q$ に二個あれば、吸収と放出の行列要素はどちらも $\sqrt3$ になる。散乱率では、これらにエネルギー保存と電子の終状態の空きをさらに掛ける。
:::

### 一振動子に急に電子を加える

局所振動子のエネルギーを $\Omega=\omega_0$ と書き、電子の占有 $n$ に対する変位後の真空を $\ket{0_n}$ とする。$D(\alpha)=e^{\alpha b^\dagger-\bar\alpha b}$ は $D^\dagger(\alpha)b D(\alpha)=b+\alpha$ を満たす。交換子 $[b,\alpha b^\dagger-\bar\alpha b]=\alpha$ が数なので、共役の展開は一次で止まる。従って $\ket{0_n}=D(-gn/\Omega)\ket0$ と書ける。二つの占有の差 $\Delta n$ による基底間の確率は $|\braket{0_{n+\Delta n}}{0_n}|^2=\exp[-g^2(\Delta n)^2/\Omega^2]$ である。電子の作用が振動より速ければ、その直後の振動状態は切り替え前のままなので、基底だけに移るとは限らない。

:::example[緩和エネルギーと垂直遷移]
$g=\Omega/2$ とする。一電子を加えた後の新しい振動基底は旧真空より $\Omega/4$ 低い。旧真空が新基底に射影される確率は $e^{-1/4}$ であり、残りは振動励起へ入る。電子を入れた直後に旧真空のままで Hamiltonian を評価すれば $\langle b+b^\dagger\rangle=0$ であるため、変位結合の平均は零である。最低の遷移エネルギーと多数の最終状態の平均に現れるエネルギーの差は、どの振動状態を測定したかで決まる。
:::

### ホッピングを残した変位変換

各サイトに同じ $\Omega>0$ の振動子を置き、結合を $g\sum_i n_i(b_i+b_i^\dagger)$ とする。:term[Lang–Firsov 変換]{en="Lang–Firsov transformation" yomi="らんぐふぃるそふへんかん"} は $W=\exp[\alpha\sum_i n_i(b_i^\dagger-b_i)]$、$\alpha=g/\Omega$ による共役である。$Wb_iW^\dagger=b_i-\alpha n_i$ が成り立つので、振動と結合の和には $\Omega\sum_i b_i^\dagger b_i-(g^2/\Omega)\sum_i n_i^2$ が残る。一方、$[n_i,c_{j\sigma}]=-\delta_{ij}c_{j\sigma}$ により、$Wc_{i\sigma}W^\dagger=c_{i\sigma}X_i$、$X_i=e^{-\alpha(b_i^\dagger-b_i)}$ である。実際、反復交換のたびに $-\alpha(b_i^\dagger-b_i)$ が一つ掛かり、これを指数にまとめられる。従って $c_{i\sigma}^\dagger c_{j\sigma}$ は $c_{i\sigma}^\dagger c_{j\sigma}X_i^\dagger X_j$ となり、ホップと振動の結合が残る。

:::example[振動を励起しないホップ]
異なる二サイト $i,j$ の振動真空上に射影すると、$\langle X_i^\dagger X_j\rangle=e^{-\alpha^2}$ となる。各サイトの真空重なりが $e^{-\alpha^2/2}$ であり、積を取るためである。例えば $g=\Omega$ では、振動真空どうしのホップ行列要素は元の $e^{-1}$ 倍になる。振動励起を含む状態を省いて有効 Hamiltonian とするには、その部門への結合を摂動として扱える必要がある。より詳しい散乱過程は [電子-フォノン相互作用](/physics/many-body/electron-phonon/) で扱う。
:::

## 平均場近似と Hartree–Fock 近似

### Slater 行列式でのエネルギー

任意の直交する $N$ 軌道を占める Slater 行列式 $\ket{\Phi}$ のみを試行状態とする。これを :term[Hartree–Fock 近似]{en="Hartree–Fock approximation" yomi="はーとりーふぉっくきんじ"} と呼ぶ。$\gamma_{ij}=\expval{a_j^\dagger a_i}_\Phi$ と定義すると、占有軌道を基底に取った際の $\gamma$ は $N$ 個の $1$ と残りの $0$ を持つ直交射影である。二体項の期待値は
$$
\expval{a_i^\dagger a_j^\dagger a_l a_k}_\Phi
=\gamma_{ki}\gamma_{lj}-\gamma_{li}\gamma_{kj}.
\label{eq-mb-hf-density}
$$
占有軌道を基底にすると、非零の項は相異なる二つの占有軌道を除いた後に同じ二軌道を戻す場合だけである。戻す順序を交換したときの符号を比較して右辺を得る。両辺は基底変換の下で同じ四重線形の変換をするため、任意の基底で成り立つ。よって
$$
E[\gamma]=\Tr(h\gamma)+\tfrac12\sum_{ijkl}\bar v_{ij;kl}\gamma_{ki}\gamma_{lj},\qquad
F_{ik}=h_{ik}+\sum_{jl}\bar v_{ij;kl}\gamma_{lj}.
\label{eq-mb-hf-fock}
$$

:::proposition[変分の停留条件]{#prop-mb-hf-stationarity}
有限基底において $\gamma$ が階数 $N$ の直交射影の集合上で $E[\gamma]$ の停留点であるための必要十分条件は $\comm{F}{\gamma}=0$ である。
:::

:::proof
階数 $N$ の直交射影は一つの射影をユニタリ共役したもの全体である。$\gamma(\epsilon)=e^{\epsilon X}\gamma e^{-\epsilon X}$、$X^\dagger=-X$ と置けば接ベクトルは $\comm{X}{\gamma}$ である。逆に射影条件を微分すると接ベクトルは占有・非占有の非対角ブロックだけをもち、任意のそのような Hermite 行列は $\comm{X}{\gamma}$ と表せる。二体項に $\bar v_{ij;kl}=\bar v_{ji;lk}$ を用いて微分すると $\delta E=\Tr(F\delta\gamma)$ だから、停留条件は全ての反 Hermite 行列 $X$ に対して $\Tr(\comm{\gamma}{F}X)=0$ となる。$X=\comm{\gamma}{F}^\dagger$ を選ぶと交換子の Hilbert–Schmidt ノルムの二乗が零になり、必要性が示される。十分性は同じ微分式から従う。
:::

### 自己無撞着性と二重計数

連立方程式は、選んだ占有軌道から $F$ を作り、その固有軌道で新しい $\gamma$ を作って解く。停留条件だけでは最小値かどうかが決まらないため、候補ごとの変分エネルギーを比較する。$F$ の占有固有値を $\epsilon_a$ とすると $\sum_a\epsilon_a=\Tr(h\gamma)+2E_{\rm int}$ なので、全エネルギーは $\sum_a\epsilon_a-E_{\rm int}$ である。Hartree 近似は粒子が平均密度の作る直接ポテンシャルを見る近似であり、Hartree–Fock ではさらに同じスピンの交換を含む。

:::example[二重占有の平均場]
$U n_\uparrow n_\downarrow$ で $m_\sigma=\expval{n_\sigma}$ と置く。演算子を $n_\sigma=m_\sigma+\delta n_\sigma$ と展開して揺らぎの積を落とすと $U(m_\downarrow n_\uparrow+m_\uparrow n_\downarrow-m_\uparrow m_\downarrow)$ である。$U=4\,\mathrm{eV}$、$m_\uparrow=m_\downarrow=1/2$ なら一体シフトは各 $2\,\mathrm{eV}$、定数項は $-1\,\mathrm{eV}$ であり、この密度における期待値は $1\,\mathrm{eV}$ である。
:::

:::example[交換による自己相互作用の消去]
一軌道だけを占めると $\gamma=\ketbra{u}{u}$ である。二体エネルギーの直接項と交換項はともに $\mel{u\otimes u}{v}{u\otimes u}$ であり、差は零である。これは粒子対の数が零であることと一致する。
:::

### 二サイトの自己無撞着方程式

半充填の Hubbard 二サイトで、$\gamma$ をスピン対角とし、各スピンの粒子を一個とする。$n_{1\uparrow}=n_{2\downarrow}=(1+m)/2$、$n_{2\uparrow}=n_{1\downarrow}=(1-m)/2$ の形を試す。上スピンの Fock 行列から定数 $U/2$ を引くと $-t\sigma_x-(Um/2)\sigma_z$ になる。下スピンでは $\sigma_z$ の符号が逆である。この行列の二乗は $E^2I$、$E=\sqrt{t^2+(Um/2)^2}$ なので、下の固有状態の射影は $\tfrac12[I+(t\sigma_x+Um\sigma_z/2)/E]$ である。従って密度差の自己無撞着条件は $m=Um/(2E)$ である。$t>0$ なら $m=0$ があり、$U>2t$ では $m^2=1-(2t/U)^2$ の解もある。Hartree–Fock 全エネルギーは定数を戻した占有固有値の和 $U-2E$ から、$E_{\rm int}=U(1-m^2)/2$ を一回引いたものになる。

:::example[二つの停留点の比較]
$t=1,U=4$ で、対称解 $m=0$ は $E=1$ で全エネルギーが $0$ である。$m=\sqrt3/2$ の解は $E=2$、$E_{\rm int}=1/2$ なので全エネルギー $-1/2$ を得る。これは指定した行列式の集合でエネルギーを下げる。有限二サイトでの厳密な基底のスピンはさらに [強相関電子系と量子磁性](/physics/many-body/strong-correlation/) で調べる。
:::

### 停留条件と安定条件の違い

占有軌道をひとつ $i$、空軌道をひとつ $a$ とし、軌道 $i$ を $i\cos\vartheta+a\sin\vartheta$ へ置き換える。他の占有軌道は固定する。新しい軌道の自己相互作用は交換で消えるため、この変分で変わるエネルギーは $h$ と不変な他の軌道から作った有効一体演算子 $F^{(i)}$ の期待値である。従って
$$
E(\vartheta)-E(0)=
(F^{(i)}_{aa}-F^{(i)}_{ii})\sin^2\vartheta
+2\Re F^{(i)}_{ai}\cos\vartheta\sin\vartheta .
$$
非対角成分 $F^{(i)}_{ai}$ は完全な $F$ の $ai$ と等しい。実際、除いた $i$ との反対称化行列要素には同じ消滅軌道が二つ並ぶため零である。同じ理由で $F^{(i)}_{ii}=F_{ii}$ だが、$F^{(i)}_{aa}=F_{aa}-\bar v_{ai;ai}$ となる。従って、停留点で一方向の曲率は $2(F_{aa}-F_{ii}-\bar v_{ai;ai})$ である。複数の軌道を同時に回せば交差する曲率も現れる。従って停留条件だけで局所最小かどうかを判定できない。

:::example[対称軌道の占有を外す効果]
$F_{aa}-F_{ii}=3\,\mathrm{eV}$、$\bar v_{ai;ai}=4\,\mathrm{eV}$ の停留点を考える。方向 $i\to a$ の厳密な差は $-1\,\mathrm{eV}\sin^2\vartheta$ なので、$\vartheta=\pi/6$ では $-0.25\,\mathrm{eV}$ になる。占有を固定した Fock 固有値には軌道 $a$ と、これから除く軌道 $i$ の相互作用も含まれている。この寄与を外すと変分でエネルギーが下がる。
:::

### Brillouin の条件と相関エネルギー

停留した行列式 $\ket\Phi$ と一重励起 $\ket{\Phi_i^a}=a_a^\dagger a_i\ket\Phi$ をとる。上の一重励起行列要素を使えば $\mel{\Phi_i^a}{H}{\Phi}=F_{ai}=0$ である。これを :term[Brillouin の定理]{en="Brillouin theorem" yomi="ぶりるあんのていり"} という。この関係は行列式の微小変分を許したときの一次の停留条件であり、二重励起との結合を禁じない。例えば二配置の行列
$$
\begin{pmatrix} E_\Phi & V^*\\ V&E_\Phi+\Delta \end{pmatrix},\qquad \Delta>0
$$
では下の固有値は $E_\Phi+(\Delta-\sqrt{\Delta^2+4|V|^2})/2$ である。$|V|\ll\Delta$ での差は $-|V|^2/\Delta+O(|V|^4/\Delta^3)$ となる。離れた軌道へ対で移る配置が、平均場の停留エネルギーを下げる仕組みである。

:::example[二配置で変分上限を改善する]
$\Delta=3\,\mathrm{eV}$、$V=2\,\mathrm{eV}$ とすると平方根は $5\,\mathrm{eV}$ であり、厳密な低下は $1\,\mathrm{eV}$ である。二次近似では $4/3\,\mathrm{eV}$ となる。低い固有状態は $(2,-1)^\top/\sqrt5$ なので、励起配置の確率は $1/5$ である。従ってこの例ではエネルギー分母は大きなパラメータとみなせない。
:::

### 二サイトの非制限平均場を厳密解と比べる

半充填の二サイトで $t>0$ とし、上スピンの占有軌道を $(\cos(\vartheta/2),\sin(\vartheta/2))^\top$、下スピンを $(\sin(\vartheta/2),\cos(\vartheta/2))^\top$ とする。スピンが直交するためこれらの空間部は互いに直交する必要がない。ホップの期待値は各スピンで $-t\sin\vartheta$、相互作用の期待値は $U\sum_i n_{i\uparrow}n_{i\downarrow}=U\sin^2\vartheta/2$ である。従って $s=\sin\vartheta\in[0,1]$ と置くと
$$
E_{\rm var}(s)=-2ts+\frac U2s^2.
$$
$U\le2t$ では境界 $s=1$ が最小で、両スピンは同じ対称軌道を占める。$U>2t$ では $s=2t/U$ が最小で、エネルギーは $-2t^2/U$ になる。$\cos\vartheta=\pm\sqrt{1-s^2}$ の二つはサイトとスピンを入れ替えた解である。いずれも $S^z=0$ だが、厳密な一重項の条件は一般に満たさない。厳密解の $-4t^2/U$ に比べて二次の係数が小さいのは、一つの行列式では交互のスピン配置を独立な振幅で重ねられないためである。

:::example[有限系での対称性を照合する]
$t=1\,\mathrm{meV}$、$U=3\,\mathrm{meV}$ なら制限した $s=1$ のエネルギーは $-1/2\,\mathrm{meV}$、$s=2/3$ の解は $-2/3\,\mathrm{meV}$、厳密な一重項は $-1\,\mathrm{meV}$ である。行列式の各サイトの上・下の占有差は $\pm\sqrt5/3$ だが、厳密一重項ではスピン回転の平均が零である。変分エネルギーの改善と、有限系の厳密な対称性は別々に比較する。
:::

## 演算子形式の Wick の定理

### 正規順序と縮約

真空に対し、:term[正規順序]{en="normal ordering" yomi="せいきじゅんじょ"} $:X_1\cdots X_m:$ は生成演算子を左に置き、Fermi 系では並べ替えの符号を掛ける操作である。ここで $X_i$ は生成・消滅演算子の線形結合であり、操作は線形に拡張する。二因子の :term[縮約]{en="contraction" yomi="しゅくやく"} は $C_{ij}=X_iX_j-:X_iX_j:=\mel{0}{X_iX_j}{0}$（$i<j$）である。例えば $C(a_i,a_j^\dagger)=\delta_{ij}$ であり $C(a_i^\dagger,a_j)=0$ である。Slater 行列式を基準にするときは、空の軌道には元の演算子を使い、占有軌道では生成と消滅の役割を取り替える。すると新しい消滅演算子が $\ket{\Phi}$ を消すため、同じ議論を使える。

:::theorem[演算子形式の Wick の定理]{#thm-mb-operator-wick}
$X_1\cdots X_m$ は、互いに素な添字対を全て選び、それぞれの縮約の積に残りの因子の正規順序積を掛けたものの和である。対のない項も含める。各符号は、各対を元の向き $i<j$ のまま隣り合わせにし、残った因子の元の相対順序を保つ置換の符号のうち、Fermi 系のものを取る。対どうしの順序は符号に影響しない。
:::

:::proof
左端 $X_1$ を生成部と消滅部に分ける。残りの正規順序積に左から掛けると生成部はそのまま正規順序を保つ。消滅部を右へ運ぶ各交換では、符号を付けた交換項と、数としての交換子または反交換子が得られる。後者が $X_1$ とその位置の因子との縮約である。各交換で数を選んだ項からは $X_1$ が消えるため、$X_1$ は高々一度しか縮約されない。縮約しない項は $X_1$ を含む正規順序積に戻る。$m=1$ からこの操作を帰納的に適用すると、残りだけの全ての組合せに対して、$X_1$ を残すか、残ったいずれか一因子と結ぶ場合がちょうど一度ずつ得られる。対の結合に必要な置換の符号と、残りを正規順序化する符号を掛けたものが各項の符号である。
:::

全縮約だけが基準状態の期待値に寄与する。四因子では期待値が $C_{12}C_{34}+\eta C_{13}C_{24}+C_{14}C_{23}$ となる。三つ目の対の選び方では、第 $4$ 因子を第 $1$ 因子の隣まで移すのに二回の交換を使う。最初に二体密度行列を求めた :eqref[eq-mb-hf-density] はこの四因子の式の具体例でもある。

:::example[真空での四因子]
$a_1a_2a_2^\dagger a_1^\dagger$ の真空期待値は $1$ である。第 $3,4$ 因子を入れ替えると、Bose 系では $1$ のまま、Fermi 系では $-1$ になる。縮約を使うと、前者の非零対は $(1,4),(2,3)$ であり、後者は $(1,3),(2,4)$ である。
:::

:::example[局所占有の相関]
二つのモードについて $\gamma_{11}=\gamma_{22}=1/2$、$\gamma_{12}=\gamma_{21}=1/2$ とする。これは一粒子の対称軌道を占めた射影である。$\expval{n_1 n_2}=1/4-1/4=0$ であり、占有数の共分散は $-1/4$ となる。密度の平均だけから積の期待値を求めると、この交換による相関が消える。
:::

### Slater 状態に対する正規順序

粒子・正孔による正規順序を作るには、空軌道 $a$ に対して $b_a=c_a$、占有軌道 $i$ に対して $b_i=c_i^\dagger$ と置く。すると全ての $b$ が $\ket\Phi$ を消す。正準反交換関係も保たれるので、空の真空に関する Wick の証明を $\ket\Phi$ へ適用できる。例えば縮約 $\langle c_i^\dagger c_j\rangle=\delta_{ij}$ は占有領域で現れ、$\langle c_a c_b^\dagger\rangle=\delta_{ab}$ は空領域で現れる。従って同じ頂点でも何を参照真空に選ぶかで零の縮約が異なる。

$H$ をこの順序で並べると、完全に縮約した項は $\langle H\rangle_\Phi$ である。二体演算子の一対の縮約は $\sum_j\bar v_{ij;kj}\langle n_j\rangle_\Phi$ のため一体項に $F$ が現れる。一対も縮約しない部分は四演算子の正規順序積となる。ゆえに
$$
H=E_\Phi+\sum_{ij}F_{ij}:c_i^\dagger c_j:_\Phi
+\frac14\sum_{ijkl}\bar v_{ij;kl}:c_i^\dagger c_j^\dagger c_l c_k:_\Phi .
$$
この式の一体項で $\ket\Phi$ から一重励起へ出る振幅は $F_{ai}$ であり、四演算子項からは一重励起へ出ない。四演算子を残したまま $\ket\Phi$ に作用させると、$b$ があれば消え、全てが $b^\dagger$ なら二粒子と二正孔を作るためである。

:::example[占有で正規順序を変える]
参照状態の局所占有を $n_{\uparrow,0}=1$、$n_{\downarrow,0}=0$ とする。$n_\uparrow=1+:n_\uparrow:_\Phi$、$n_\downarrow=:n_\downarrow:_\Phi$ なので、
$$
Un_\uparrow n_\downarrow=U:n_\downarrow:_\Phi+U:n_\uparrow n_\downarrow:_\Phi .
$$
右辺の四演算子項は $U(n_\uparrow-1)n_\downarrow$ である。実際に $\ket{0,1}$ で一体項と四演算子項は $U$ と $-U$ になり、元の相互作用は零になる。
:::

### 行列式にまとめるときの順序

数を保存する Slater 状態について、$p$ 本の生成と $p$ 本の消滅を
$c_{i_1}^\dagger\cdots c_{i_p}^\dagger c_{j_p}\cdots c_{j_1}$ の順に並べる。$i_1$ をどの $j$ と縮約するかで Wick の和を分けると、残った成分の符号は $(-1)^{1+l}$ になる。消滅添字 $j_l$ を $c_{j_1}$ に対応する位置へ移すと $l-1$ 回の交換を要し、その後、同じ順序の $p-1$ 対が残るからである。従って $p=1$ からの帰納法で
$$
\langle c_{i_1}^\dagger\cdots c_{i_p}^\dagger c_{j_p}\cdots c_{j_1}\rangle
=\det(\gamma_{j_l i_k})_{k,l=1}^p .
$$
消滅演算子を逆に並べた場合は $(-1)^{p(p-1)/2}$ が必要である。密度の高次積を計算するときは、まず同一軌道での $n_i^2=n_i$ を使い、異なる軌道の積に対してこの行列式を用いる。

:::example[二軌道の同時占有]
一電子の軌道 $(u_1+u_2)/\sqrt2$ では $\gamma_{ij}=1/2$ である。従って $\langle n_1n_2\rangle$ は $1/4-1/4=0$ である。各軌道に電子がいる確率の積は $1/4$ だが、交換縮約を引くことで固定した粒子数と整合する。
:::

### 反交換子から三対の符号を追う

直交した三軌道 $1,2,3$ を使う。真空平均 $\langle c_i c_j c_k c_3^\dagger c_2^\dagger c_1^\dagger\rangle$ は、消滅列が $(i,j,k)=(1,2,3)$ のとき一である。実際、右から $c_k$ が $c_3^\dagger$ を消し、残りが $c_2^\dagger c_1^\dagger$ に対して同様に働く。消滅列を一回交換するたびに符号が反転し、重複添字があれば零になる。行列式の表式では $\det(\delta_{\text{消滅添字},\text{生成添字}})$ を生成添字 $1,2,3$ の順に取る。生成の列を逆順に記す理由は、基準となる $(1,2,3)$ に正の符号を与えるためである。

:::example[重ね合わせを消滅させる]
$c'_1=(c_1+c_2)/\sqrt2$、$c'_2=(c_1-c_2)/\sqrt2$、$c'_3=c_3$ を取る。縮約行列は
$$
\begin{pmatrix}
1/\sqrt2&1/\sqrt2&0\\
1/\sqrt2&-1/\sqrt2&0\\
0&0&1
\end{pmatrix}
$$
で、その行列式は $-1$ である。非零の積は $-1/2$ と $1/2$ である。前者は恒等置換、後者は生成添字 $1,2$ の交換なので、置換の符号を入れた和は $-1/2-1/2=-1$ である。実際に $c'_1c'_2=(c_1c_1-c_1c_2+c_2c_1-c_2c_2)/2=-c_1c_2$ と直接展開すると同じ負号が出る。
:::

場の Wick 展開を使うと、正規順序化していない積の等時刻のデルタ項も系統的に残せる。二軌道の異なる場に対し $n_i n_j=c_i^\dagger c_j^\dagger c_j c_i$ であるが、$i=j$ の場合は右辺が零で左辺は $n_i$ になる。従って添字を全て含む式は $n_i n_j=c_i^\dagger c_j^\dagger c_j c_i+\delta_{ij}n_i$ である。$N=\sum_i n_i$ の二乗ではデルタ項を含めて $N^2=N+\sum_{ij}c_i^\dagger c_j^\dagger c_j c_i$ になる。粒子対の演算子は $N(N-1)/2$ であり、同じ軌道を二回数える項は入らない。

:::example[粒子数の分散と対の数]
独立な二軌道を確率 $1/4,3/4$ で占有し、他を空とする混合を取る。全数の平均は $1$、対の平均は $(1/4)(3/4)=3/16$ である。従って $\langle N^2\rangle=1+2(3/16)=11/8$、分散は $3/8$ である。等時刻のデルタ項を外して $\langle N^2\rangle=3/8$ とすると、平均を引いた分散は負になってしまう。
:::

## 演習問題

::::exercise{level=basic}
Bose 基底 $\ket{3,1}$ に $a_2^\dagger a_1$ を作用させよ。

:::solution
右の因子が第一モードを $3$ から $2$ に減らし $\sqrt3$ を出す。次の因子は第二モードを $1$ から $2$ に増やし $\sqrt2$ を出す。したがって $\sqrt6\ket{2,2}$ である。
:::
::::

::::exercise{level=basic}
三つの Fermi 軌道の順序を $1,2,3$ とする。$a_2^\dagger a_3\ket{1,0,1}$ を求めよ。

:::solution
$a_3$ は左に占有軌道を一つ見るので $-\ket{1,0,0}$ を与える。$a_2^\dagger$ も軌道 $1$ のため符号が $-1$ である。積は $+\ket{1,1,0}$ となる。
:::
::::

::::exercise{level=basic}
$U=3\,\mathrm{meV}$ の一モード Bose 系で、$n=2$ から $3$ に増やす相互作用エネルギーを求めよ。

:::solution
各エネルギーは $Un(n-1)/2$ である。$E_3-E_2=U(3-1)=6\,\mathrm{meV}$ となる。この差に一体エネルギーを加えると粒子追加に必要な全エネルギーになる。
:::
::::

::::exercise{level=basic}
四点周期鎖で $t=2\,\mathrm{eV}$ とする。$k=\pi/a$ の一粒子エネルギーと群速度を求めよ。

:::solution
分散に $\cos\pi=-1$ を代入すると $4\,\mathrm{eV}$ である。速度は $2ta\sin(ka)/\hbar$ であり、$\sin\pi=0$ のため零である。
:::
::::

::::exercise{level=advanced}
一体演算子の Fock 表示を $\widehat h$ と書く。$\comm{\widehat h}{\widehat g}=\widehat{\comm{h}{g}}$ を証明せよ。

:::solution
$\comm{a_i^\dagger a_j}{a_k^\dagger a_l}$ を展開し、各積で中央の消滅と生成を :eqref[eq-mb-ccr-car] で交換すると $\delta_{jk}a_i^\dagger a_l-\delta_{il}a_k^\dagger a_j$ が残る。四因子項は生成どうしと消滅どうしを交換する二つの符号が相殺するため消える。$h_{ij}g_{kl}$ で重みを付けると二つのデルタがそれぞれ $hg$ と $gh$ の行列積を作る。
:::
::::

::::exercise{level=advanced}
二サイト間で $h=-t\begin{pmatrix}0&1\\1&0\end{pmatrix}$、$t>0$ とする。一電子を置く場合の Hartree–Fock 最小状態とエネルギーを求めよ。

:::solution
階数一の射影では直接項と交換項が相殺するので、$E[u]=\mel{u}{h}{u}$ のみである。正規化した対称・反対称軌道の固有値はそれぞれ $-t,+t$ である。$u=\alpha u_++\beta u_-$ と展開すれば $E=-t+2t|\beta|^2$ なので、$\beta=0$ のとき最小値 $-t$ を取る。
:::
::::

::::exercise{level=advanced}
Fermi 真空で $\mel{0}{a_3 a_1 a_2 a_1^\dagger a_2^\dagger a_3^\dagger}{0}$ を求めよ。

:::solution
非零の縮約は $(1,6),(2,4),(3,5)$ だけである。第 $6$ 因子を第 $1$ 因子の次へ移すのに四交換を要する。残った列で元の第 $4$ 因子を第 $2$ 因子の次へ移すのに一交換を要する。縮約は全て $1$ であるため値は $(-1)^5=-1$ である。
:::
::::

::::exercise{level=advanced}
独立振動子をもつ一サイトで $H=\omega_0 b^\dagger b+gn(b+b^\dagger)$、$\omega_0>0$、$g\in\R$ とする。固定した $n$ のスペクトルを求めよ。

:::solution
$B=b+gn/\omega_0$ と置くと $[B,B^\dagger]=1$ であり、平方完成により $H=\omega_0 B^\dagger B-g^2n^2/\omega_0$ となる。$B$ の真空は $b$ の固有値 $-gn/\omega_0$ をもつ規格化状態で、元の占有数基底で係数は $e^{-g^2n^2/(2\omega_0^2)}(-gn/\omega_0)^m/\sqrt{m!}$ である。変位はユニタリ $e^{(gn/\omega_0)(b-b^\dagger)}$ で生成でき、$B^\dagger B$ の固有値は $0,1,\ldots$ である。よって $E_m=m\omega_0-g^2n^2/\omega_0$ である。
:::
::::

::::exercise{level=advanced}
二つの軌道 $u_1,u_2$ に二個の Bose 粒子を置き、$\ket{1,1}$ を考える。基底を $\phi_\pm=(u_1\pm u_2)/\sqrt2$ に変えたとき、出力の占有数とその確率を求めよ。

:::solution
新しい生成演算子を $b_\pm^\dagger$ とすると逆変換は $a_1^\dagger=(b_+^\dagger+b_-^\dagger)/\sqrt2$、$a_2^\dagger=(b_+^\dagger-b_-^\dagger)/\sqrt2$ である。積の交差項は Bose の可換性で相殺されるので、
$$
a_1^\dagger a_2^\dagger\ket0
=\tfrac12[(b_+^\dagger)^2-(b_-^\dagger)^2]\ket0
=(\ket{2,0}_b-\ket{0,2}_b)/\sqrt2 .
$$
従ってどちらかの軌道に二個集まる確率は各 $1/2$ である。一個ずつに分かれる確率は零である。
:::
::::

::::exercise{level=advanced}
二サイト Hubbard 模型の一重項ブロック :eqref[eq-mb-hubbard-dimer-singlet] について $t>0$ を固定し、$U/t\to\infty$ で基底における合計の二重占有確率を $t/U$ の四次まで求めよ。

:::solution
$E_-=(U-\sqrt{U^2+16t^2})/2$ を使う。$t>0$ では下の固有値は単純であり、固有値方程式の微分により確率は
$$
P_D=\partial_U E_-=\frac12\left(1-\frac{U}{\sqrt{U^2+16t^2}}\right).
$$
$x=16t^2/U^2$ とすると $(1+x)^{-1/2}=1-x/2+3x^2/8+O(x^3)$ なので、$P_D=x/4-3x^2/16+O(x^3)=4t^2/U^2-48t^4/U^4+O((t/U)^6)$ である。これと独立に、ブロックの固有方程式から成分比 $D_+/S=-E_-/(2t)$ が出る。その二乗を $1+(E_-/(2t))^2$ で割って同じ展開を得る。先頭の確率が振幅の二乗の次数になるため、エネルギー差にも仮想過程の二次が現れる。
:::
::::

::::exercise{level=basic}
同じ一粒子エネルギー $\epsilon$ をもつ三つの軌道を考える。二体相互作用のない Bose または Fermi の二粒子部門の分配関数を求めよ。粒子を区別できる場合とも比較せよ。

:::solution
エネルギーはいずれの配置でも $2\epsilon$ なので、分配関数は部門の次元に $e^{-2\beta\epsilon}$ を掛けたものである。Bose では三軌道のうち一つに二個置く配置が三つ、異なる二軌道を一つずつ占める配置が $\binom32=3$ 個あり、合計六状態である。Fermi では後者だけなので三状態である。区別できる粒子はそれぞれ独立に三軌道を選ぶので九状態である。Bose と Fermi では異なる二軌道を選んだ状態を入れ替えても独立な二状態にはならず、対称化または反対称化された一状態だけを数える。
:::
::::

## 確認テスト

::::quiz
二つの Fermi 軌道について正しいものを選べ。

- [x] $a_1a_2=-a_2a_1$ である。
- [ ] $(a_1^\dagger)^2$ は粒子数 $2$ の状態を作る。
- [x] $a_1^\dagger a_2$ は全粒子数を保存する。

:::explanation
正準反交換関係から一つ目が成り立つ。同じ生成演算子の反交換子が零なので二つ目の作用は零である。三つ目では生成と消滅が一回ずつのため、粒子数演算子との交換子が零になる。
:::
::::

::::quiz{answer="18" tol="1e-6"}
$U=3\,\mathrm{meV}$、$n=4$ の一軌道 Bose 相互作用エネルギーを $\mathrm{meV}$ で答えよ。

:::explanation
粒子対は $4\cdot3/2=6$ 個である。したがって $U n(n-1)/2=18\,\mathrm{meV}$ となる。
:::
::::

::::quiz{answer="8" tol="1e-6"}
$t=2\,\mathrm{eV}$ の最近接鎖で、熱力学極限の帯域幅を $\mathrm{eV}$ で答えよ。

:::explanation
$\epsilon(k)=-2t\cos(ka)$ は最低値 $-2t$、最高値 $2t$ を取る。差は $4t=8\,\mathrm{eV}$ である。
:::
::::

::::quiz
Hartree–Fock エネルギーを計算するとき、占有した Fock 固有値の和にはどの補正が要るか。

- [x] 相互作用エネルギーを一回引く。
- [ ] 相互作用エネルギーを一回足す。
- [ ] 任意の相互作用に対して補正なしでよい。

:::explanation
Fock 演算子の相互作用成分を占有軌道で和を取ると粒子対を二度数える。よって $\sum_a\epsilon_a=\Tr(h\gamma)+2E_{\rm int}$ から $E_{\rm int}$ を一回引けば $E=\Tr(h\gamma)+E_{\rm int}$ になる。
:::
::::

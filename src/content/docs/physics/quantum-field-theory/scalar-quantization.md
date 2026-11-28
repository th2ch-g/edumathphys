---
title: "スカラー場の正準量子化"
description: "Klein–Gordon 場の正準量子化から粒子、反粒子、因果性と伝播関数を導く。"
notation:
  - symbol: 'E_{\vb p}'
    meaning: 正の質量殻エネルギー
  - symbol: '\dd\Pi_p'
    meaning: 正質量殻上の不変測度
  - symbol: '\Delta_F,G_R'
    meaning: Feynman 伝播関数と遅延 Green 関数
---

## この章の目標

- 場と運動量の交換関係から粒子モードを構成する。
- エネルギー、運動量と電荷を占有数で表す。
- 空間的分離での交換子と Green 関数の境界条件を区別する。

## Klein–Gordon 場の正準量子化

:eqref[eq-qft-euler] の実スカラー場で $J=0$ とする。この自由方程式 $(\Box+m^2)\phi=0$ を :term[Klein–Gordon 方程式]{en="Klein–Gordon equation" yomi="くらいんごるどんほうていしき"} と呼ぶ。各 Fourier モードは周波数 $E_{\vb p}=\sqrt{\vb p^2+m^2}$ の調和振動子なので、各モードに振動子の量子化を適用できる。$\hbar=c=1$、以下 $px=E_{\vb p}t-\vb p\cdot\vb x$ とし、点の場は試験関数と積分して定める演算子値分布である。

:::postulate[等時刻正準交換関係]{#post-qft-scalar-ccr}
$\pi=\partial_t\phi$ とし、等時刻で
$$
[\phi(\vb x),\pi(\vb y)]=i\delta^3(\vb x-\vb y),\qquad
[\phi,\phi]=[\pi,\pi]=0
\label{eq-qft-scalar-ccr}
$$
と課す。
:::

$\dd\Pi_p=\dd^3p/[(2\pi)^3 2E_{\vb p}]$ と書く。量子化の際には古典解の正・負周波数の係数を互いに随伴な演算子とする。これは実数値という条件を Hermite 性へ移したものである。
$$
\phi(x)=\int\dd\Pi_p\,[a(\vb p)e^{-ipx}+a^\dagger(\vb p)e^{ipx}].
\label{eq-qft-scalar-expansion}
$$
逆式は $a(\vb p)=i\int\dd^3x\,e^{ipx}(\dot\phi-iE_{\vb p}\phi)$ である。空間積分で同じ $\vb p$ の正周波数項は $2E_{\vb p}$ を生じ、逆向きの負周波数項は係数が零になる。両逆式を :eqref[eq-qft-scalar-ccr] に代入すると
$$
[a(\vb p),a^\dagger(\vb q)]=(2\pi)^3 2E_{\vb p}\delta^3(\vb p-\vb q),\qquad
[a,a]=[a^\dagger,a^\dagger]=0.
\label{eq-qft-scalar-modes}
$$
逆にこの関係を場へ戻せば、$[\phi,\pi]$ の二項が等時刻に Fourier デルタの半分ずつを与える。

### 解空間の内積と真空の選択

複素化した自由解について $(u,v)_{\rm KG}=i\int\dd^3x\,u^*\overleftrightarrow{\partial_t}v$ を置く。時刻微分は $i\int(u^*\ddot v-\ddot u^*v)\dd^3x$。両方の方程式を入れると質量項が消え、空間 Laplacian の差が表面項となる。従ってこの積は時刻と無関係である。箱で $u_{\vb p}=e^{-ipx}/\sqrt{2VE_{\vb p}}$ とすると $(u_{\vb p},u_{\vb q})_{\rm KG}=\delta_{\vb p\vb q}$、$(u_{\vb p}^*,u_{\vb q}^*)_{\rm KG}=-\delta_{\vb p\vb q}$。混合積では Fourier 選択条件 $\vb q=-\vb p$ の下で周波数差が零となる。従って正周波数部分に制限すると正定値の一粒子内積になる。一般の背景で時間並進がない場合は、どの解を正周波数と呼ぶかも追加で指定する必要がある。周期箱の質量ゼロの零モードは調和振動子に含まれず、有限の質量を置くか零モードを別に扱う。

### 零モードと質量ゼロ極限

周期箱で空間平均を $q_0=V^{-1/2}\int\phi\,\dd^3x$、$p_0=V^{-1/2}\int\pi\,\dd^3x$ とする。等時刻の交換関係を両点で積分すると $[q_0,p_0]=i$。運動項は $p_0^2/2$、勾配は零なので、Hamiltonian は $(p_0^2+m^2q_0^2)/2$ である。従って $m>0$ の基底波動関数は $\exp(-mq_0^2/2)$ に比例する。Gauss 積分を微分して $\expval{q_0^2}=1/(2m)$ を得る。$m=0$ では自由粒子 $p_0^2/2$ となり、非コンパクトな場の値を取る限り正規化可能な基底状態は存在しない。

:::example[箱の大きさとスペクトル間隔]
一辺 $L=2\pi\,\mathrm{GeV}^{-1}$ では $\vb p=(n_1,n_2,n_3)\,\mathrm{GeV}$。$m=2\,\mathrm{GeV}$ なら零モードと $\vb n=(1,2,0)$ の角周波数は $2$ と $3\,\mathrm{GeV}$ である。$m=0$ にしても後者の周波数は $\sqrt5\,\mathrm{GeV}$ であり、零モードの規格化の問題はこのモードに波及しない。
:::

### 実場の Fourier モードを正準変数として数える

一辺 $L$ の周期箱で $V=L^3$、$\vb k=2\pi\vb n/L$ とする。零でない $\vb k$ の対 $\{\vb k,-\vb k\}$ から片方だけを取る集合を $\mathcal K$ とし、実基底を使う。
$$
\phi(t,\vb x)=\frac{q_0(t)}{\sqrt V}
+\sqrt{\frac2V}\sum_{\vb k\in\mathcal K}
\bigl[q_{\vb k,c}(t)\cos(\vb k\cdot\vb x)+q_{\vb k,s}(t)\sin(\vb k\cdot\vb x)\bigr].
$$
複素指数の積を箱に積分すると、運動量の和が零のときだけ $V$ が残る。従ってこの実基底は直交してノルム一。勾配項は部分積分して $-\nabla^2$ に移すと各モードの固有値 $|\vb k|^2$ が出る。作用は零モード一つと、各 $\{\vb k,-\vb k\}$ に二つの実振動子の和となる。運動量 $\pi=\dot\phi$ も同じ基底で展開すれば、モードの $p=\dot q$ は $\int\pi\,\delta\phi\,\dd^3x=\sum p\,\delta q$ を満たす。これで周期条件での等時刻交換関係は、各実振動子の正準交換関係を場へ戻したものになる。

:term[定在波]{en="standing wave" yomi="ていざいは"} の二つの座標から、運動量を指定した :term[進行波]{en="traveling wave" yomi="しんこうは"} の消滅演算子を作る。周波数を $\omega>0$ とし、$a_{c,s}=\sqrt{\omega/2}\,q_{c,s}+i p_{c,s}/\sqrt{2\omega}$。すると
$$
a_{\vb k}=\frac{a_c-i a_s}{\sqrt2},\qquad
a_{-\vb k}=\frac{a_c+i a_s}{\sqrt2}.
$$
この二行の係数行列はユニタリなので、両進行波の交換関係は互いに独立である。空間 Fourier の係数 $q_{-\vb k}=q_{\vb k}^*$ は古典場が実である条件であり、$a_{-\vb k}=a_{\vb k}^\dagger$ という条件ではない。実際、正周波数側で $a_{\vb k}$ が $e^{i\vb k\cdot\vb x}$、$a_{-\vb k}$ が $e^{-i\vb k\cdot\vb x}$ を持ち、その随伴がそれぞれ負周波数側を作る。

:::example[定在波の一量子と進行波の確率]
余弦モードの一量子は $a_c^\dagger\ket0=(a_{\vb k}^\dagger+a_{-\vb k}^\dagger)\ket0/\sqrt2$。運動量が異なる二つの状態の直交性から、運動量を測ると $+\vb k$ と $-\vb k$ がそれぞれ確率 $1/2$ で出る。一方、両状態の角周波数は同じなので、エネルギーはこの測定で確率的に分かれない。例えば $m=3\,\mathrm{GeV}$、$|\vb k|=4\,\mathrm{GeV}$ ならエネルギーは $5\,\mathrm{GeV}$、運動量の平均は零である。平均が零という条件は、零運動量モードに一量子がいる条件と異なる。
:::

### Klein–Gordon 内積が初期面に依らないこと

二つの複素解 $u,v$ に対して $J^\mu(u,v)=i(u^*\partial^\mu v-(\partial^\mu u^*)v)$ と定める。発散は $i(u^*\Box v-\Box u^*v)$ なので、両方が同じ実の質量 $m$ を持てば零。二つの空間的 Cauchy 面 $\Sigma_1,\Sigma_2$ で、$J$ の無限遠の流束が消える波束を考える。両面の間の四次元領域に発散定理を使うと
$$
(u,v)_{\Sigma_2}-(u,v)_{\Sigma_1}=0,\qquad
(u,v)_\Sigma=\int_\Sigma\dd\Sigma_\mu\,J^\mu(u,v).
$$
ここで $\dd\Sigma_\mu$ は未来向き法線を使った面要素である。Cauchy 面とは、両向きに延長できる限り延長した任意の因果的曲線を一度だけ切る初期面を指す。$t=\text{定数}$ なら成分は $\dd\Sigma_0=\dd^3x$。ブーストで固定時刻面を別の面へ移しても値は保たれる。これは時間発展での保存を、Lorentz 共変な面積分として書く方法である。

平面波を有限箱で正規化し $u_{\vb p}(\vb x,t)=e^{-iEt+i\vb p\cdot\vb x}/\sqrt{2EV}$ とすると $(u_{\vb p},u_{\vb q})=\delta_{\vb p\vb q}$、$(u_{\vb p}^*,u_{\vb q}^*)=-\delta_{\vb p\vb q}$。異なる周波数符号の内積は空間デルタで $\vb p=-\vb q$ となったとき分子が $E_{\vb p}-E_{\vb q}=0$ なので零。従って解全体では不定値であるが、未来質量殻に台をもつ正周波数解の空間では正になる。連結 Lorentz 変換が未来殻を保つため、慣性系のブーストでこの選択を失わない。場の自由解にこの内積を使うと $a_{\vb p}=(u_{\vb p},\phi)$、$a_{\vb p}^\dagger=-(u_{\vb p}^*,\phi)$ が得られる。負周波数のノルムの符号が、二本目の式に負号を付ける理由である。

:::example[負周波数を含む解の内積]
正規化した一つの $u$ とその共役に対し $f=2u+u^*$、$g=u+2u^*$ とする。$(f,f)=4-1=3$、$(g,g)=1-4=-3$、$(f,g)=2-2=0$。さらに実解 $u+u^*$ は非零だが自身との積は零になる。従ってこの積を古典的な実解のノルムと解釈できない。量子化で一粒子 Hilbert 空間に使うのは選択した正周波数の部分である。
:::

### 質量殻の測度と二つのモード規格化

Lorentz 不変な積分の出発点を $\dd^4p\,\theta(p^0)\delta(p^2-m^2)$ とする。$m>0$ または $\vb p\ne0$ では $p^0$ の根が $\pm E_{\vb p}$ の二つで、それぞれの微分の絶対値は $2E_{\vb p}$。従って
$$
\int \dd p^0\,\theta(p^0)\delta((p^0)^2-|\vb p|^2-m^2)f(p^0)
=\frac{f(E_{\vb p})}{2E_{\vb p}} .
$$
連結 Lorentz 変換は四体積と $p^2$ に加え、未来殻でのエネルギーの符号を保つ。従って残りの三次元測度も、$\vb p$ を殻の上で変換すると不変である。別の流儀では $\phi$ を $\int \dd^3p/(2\pi)^3 [\widetilde a(p)e^{-ipx}+\widetilde a^\dagger(p)e^{ipx}]/\sqrt{2E}$ として、$\widetilde a$ の交換関係を $(2\pi)^3\delta^3$ に取る。:eqref[eq-qft-scalar-expansion] の $a$ とは $a(p)=\sqrt{2E}\,\widetilde a(p)$。この係数を交換関係でも変えれば両規約は同じである。

:::example[測度の不変性を軸方向の Jacobian で確かめる]
$z$ 軸のブーストでは $p_z'=\cosh\chi\,p_z+\sinh\chi\,E$、横二成分は変わらない。従って殻上の変数変換の行列式は $\partial p_z'/\partial p_z=\cosh\chi+\sinh\chi\,p_z/E=E'/E$。$m=3\,\mathrm{GeV}$、$p_z=0$、$\cosh\chi=5/3,\sinh\chi=4/3$ なら $E'=5\,\mathrm{GeV}$、$p_z'=4\,\mathrm{GeV}$ で比は $5/3$。$\dd^3p/(2E)$ では分母にも $5/3$ が掛かって相殺する。$\widetilde a$ の規格化のまま平面波の状態を変換する場合は、逆にこの比の平方根を状態の変換則に含める。
:::

### 真空二点を正定値の分布として使う

正周波数の二点分布を $D^+(x-y)=\langle0|\phi(x)\phi(y)|0\rangle$ とし、滑らかで急減少する複素試験関数 $f$ に対し $\phi(f)=\int\dd^4x\,f(x)\phi(x)$ とする。随伴は $\phi(f^*)$ である。$\widehat f(p)=\int\dd^4x\,e^{ipx}f(x)$ の規約では
$$
\langle0|\phi(f)^\dagger\phi(f)|0\rangle
=\int\dd\Pi_p|\widehat f(p)|^2 .
$$
右辺は非負であり、一粒子ノルムになる。時空の関数として非零でも正殻上の Fourier 係数が零の $f$ はこのノルムで零になる。例えば $f=(\Box+m^2)g$ なら $\widehat f=(m^2-p^2)\widehat g$ なので、$\phi(f)=0$ という演算子の方程式を満たす。従って試験関数の空間全体をそのまま一粒子空間とせず、このような零の方向を除いた後に完備化する。

:::example[質量殻に沿う試験関数と外れた試験関数]
$\widehat f(p)$ が正殻から離れた領域に台をもつなら $\phi(f)\ket0=0$。一方、質量 $m>0$ で $\widehat g(p)=e^{-[(p^0-m)^2+|\vb p|^2]/(2\sigma^2)}$ を取ると、静止点 $(m,\vb0)$ では値が一で、その周りに正の測度の領域がある。従って $\|\phi(g)\ket0\|>0$。$\sigma$ を小さくすると波束は静止近くの運動量へ集中する。このとき空間の幅は Fourier の幅の逆で大きくなる。時空の関数を一点へ絞る極限と、一粒子の運動量を一点へ絞る極限は同じ操作ではない。
:::

## 生成消滅演算子と Fock 空間

### 一振動子から対称テンソル積へ

箱の体積を $V=L^3$ とし、$\vb p=2\pi\vb n/L$ とする。場の各モードを $1/\sqrt{2VE_{\vb p}}$ で規格化すれば $[a_{\vb p},a_{\vb q}^\dagger]=\delta_{\vb p\vb q}$ となる。真空 $\ket0$ を単位ベクトルで $a_{\vb p}\ket0=0$ とする。一振動子では $a(a^\dagger)^n=(a^\dagger)^na+n(a^\dagger)^{n-1}$ なので、帰納法で $(a^\dagger)^n\ket0$ のノルム二乗は $n!$。規格化した状態に $a$ は $\sqrt n$、$a^\dagger$ は $\sqrt{n+1}$ を掛けて占有数を一つ変える。例えば二粒子状態に $a$ を作用させると $\sqrt2\ket1$ である。

:::definition[ボース Fock 空間]{#def-qft-scalar-fock}
一粒子空間を $\mathscr H_1=L^2(\dd\Pi_p)$ とし、対称化した $n$ 粒子テンソル積の Hilbert 直和 $\bigoplus_{n\ge0}\operatorname{Sym}^n\mathscr H_1$ を :term[Fock 空間]{en="Fock space" yomi="ふぉっくくうかん"} と呼ぶ。$n=0$ の成分は $\C\ket0$ である。
:::

$a^\dagger(f)=\int\dd\Pi_p\,f(p)a^\dagger(\vb p)$ とすれば $[a(f),a^\dagger(g)]=\int\dd\Pi_p\,f^*g$ である。従って規格化波束一個のノルムは $\int\dd\Pi_p\,|f|^2$。生成演算子同士は可換であり、粒子ラベルを入れ替えても同じ状態となる。完全な対称化はこの性質を波動関数の立場から表したものである。

### 波束の重なりと二粒子のノルム

$a^\dagger(f)=\int\dd\Pi_p\,f(p)a^\dagger(p)$ とし、$a(f)=(a^\dagger(f))^\dagger$ とする。:eqref[eq-qft-scalar-modes] により $[a(f),a^\dagger(g)]=\int\dd\Pi_p f(p)^*g(p)$。正規化した $f,g$ の二粒子状態 $a^\dagger(f)a^\dagger(g)\ket0$ に $a(f)$ を作用させると $a^\dagger(g)\ket0+\inner{f}{g}a^\dagger(f)\ket0$ となる。残る $a(g)$ を作用させてノルム二乗 $1+|\inner{f}{g}|^2$ を得る。これは内積を第 2 引数について線形とする規約である。

:::example[部分的に重なる二つのモード]
$f,h$ が正規直交し、$g=(f+\sqrt3h)/2$ なら二粒子ノルム二乗は $5/4$。実際、二粒子状態は $(a^\dagger(f))^2\ket0/2+\sqrt3a^\dagger(f)a^\dagger(h)\ket0/2$ である。二つの項の占有数が違うため干渉せず、規格化後に $f$ が二個見つかる確率は $(2/4)/(5/4)=2/5$ となる。
:::

### Fock 空間での生成消滅作用

一粒子空間を $\mathcal H_1$ とする。$n$ 粒子の対称射影は $S_n=\sum_{\sigma\in\mathfrak S_n}U_\sigma/n!$ である。置換作用の随伴は逆置換だから $S_n^\dagger=S_n$、二重和では各積の置換が $n!$ 回現れるため $S_n^2=S_n$。従って $S_n$ は対称部分への直交射影になる。有限粒子数のベクトルに対し $a^\dagger(f)\psi_n=\sqrt{n+1}\,S_{n+1}(f\otimes\psi_n)$ と置く。消滅 $a(f)$ をこの随伴とすれば、対称波動関数の一変数を $f$ で内積した写像に $\sqrt n$ が付く。これは $a^\dagger$ の定義を $n$ 粒子の対称状態と内積し、$S_n$ をそちらへ移せば分かる。

$a(f)a^\dagger(g)\psi_n$ では、消す粒子と作った粒子が同じ場合に $\inner f g\psi_n$ が残る。他の $n$ 項はもともとの粒子の一つを消して $g$ を作る。対称化の $1/(n+1)$ は生成と消滅の $\sqrt{n+1}$ 二つに相殺され、後者は $a^\dagger(g)a(f)\psi_n$ に一致する。従って $[a(f),a^\dagger(g)]=\inner f g$。作る演算子同士は最後の対称射影で因子の順番を交換できるため可換である。有限粒子数の稠密部分はこれらの積を定義する共通領域であり、数演算子などの非有界な作用を Fock 空間の全ベクトルに定義したと仮定しない。

:::example[三粒子の規格化に出る重複]
相異なる直交モード $f,g$ に二粒子と一粒子を置く。規格化した状態は $(a^\dagger(f))^2a^\dagger(g)\ket0/\sqrt2$。消滅演算子を左から移せば、$a(f)$ は $a^\dagger(f)$ の二因子と交換し、$a(g)$ は $g$ の一因子とだけ交換する。従って $f$ を一つ消したノルムは $\sqrt2$、$g$ を一つ消したノルムは $1$。この区別は、粒子総数が三であることだけからは決まらない。
:::

### Bogoliubov 変換と真空

一振動子の $[a,a^\dagger]=1$ に対し $b=u a+v a^\dagger$ とする。$[b,b^\dagger]=|u|^2-|v|^2$ なので $b$ も正準変数にする条件は $|u|^2-|v|^2=1$。例えば実係数なら $u=\cosh r,v=\sinh r$ と書ける。$b$ の消滅真空を $\ket{0_b}=\sum_n c_n\ket{n_a}$ とすると
$$
u\sqrt{n+1}\,c_{n+1}+v\sqrt n\,c_{n-1}=0,\qquad c_{-1}=0.
$$
$n=0$ から $c_1=0$、以後奇数成分はすべて零。偶数成分は
$$
c_{2n}=c_0(-v/u)^n\frac{\sqrt{(2n)!}}{2^n n!}.
$$
ノルムの和は $\sum_n\binom{2n}n(|v/u|^2/4)^n=(1-|v/u|^2)^{-1/2}$。従って $|c_0|^2=\sqrt{1-|v/u|^2}=1/|u|$。逆変換 $a=u^* b-v b^\dagger$ によって $\bra{0_b}a^\dagger a\ket{0_b}=|v|^2$ となる。これを :term[Bogoliubov 変換]{en="Bogoliubov transformation" yomi="ぼごりゅーぼふへんかん"} と呼ぶ。正準交換関係を保っても、どのモードを消滅と呼ぶかを変えると真空の粒子数も変わる。無限個の振動子の変換では、一モードの規格化が可能でもその積が元の Fock 空間に残るとは限らない。

:::example[周波数を瞬時に変えた振動子]
質量を一にした振動子で、切替時に座標 $q$ と運動量 $p$ を同じ正準変数のまま、周波数を $\omega_a$ から $\omega_b$ へ変える。$a=\sqrt{\omega_a/2}q+i p/\sqrt{2\omega_a}$ を解いて $q,p$ を $b$ の定義に代入すると
$$
u=\tfrac12(\sqrt{\omega_b/\omega_a}+\sqrt{\omega_a/\omega_b}),\qquad
v=\tfrac12(\sqrt{\omega_b/\omega_a}-\sqrt{\omega_a/\omega_b})
$$
となる。$\omega_b/\omega_a=4$ なら $u=5/4,v=3/4$。切替前の真空の新粒子数は $v^2=9/16$ であり、零点を引いた新エネルギーは $9\omega_b/16$。時間に依存しない自由場のブーストではこのような混合は起こらない。この例では Hamiltonian を時間に依存させたため、前後の正周波数を別々に定めた。
:::

### 非直交の多粒子状態と permanent

$f_1,\ldots,f_n$ と $g_1,\ldots,g_n$ の二つの波束の組を用い、真空上で生成演算子を積む。最初の消滅演算子を右へ移すと、$n$ 個の生成演算子の一つと交換して内積を出す項が $n$ 個ある。最後まで交換子を取らない項は真空を消す。残りに同じ操作を反復するので
$$
\bra0 a(f_n)\cdots a(f_1)a^\dagger(g_1)\cdots a^\dagger(g_n)\ket0
=\sum_{\sigma\in\mathfrak S_n}\prod_{i=1}^n\inner{f_i}{g_{\sigma(i)}} .
$$
右辺は Gram 型の行列の :term[パーマネント]{en="permanent" yomi="ぱーまねんと"} である。各置換の符号を付けずに足す。フェルミ場では演算子を交換する符号が残って行列式になる。ボースとフェルミの違いは、独立なモードの規格化だけでなく、波束が重なるときに現れる。

:::example[二つが同じ波束で第三だけが異なる場合]
$f$ と $g$ が規格化され、$c=\inner f g$ とする。状態 $(a^\dagger(f))^2a^\dagger(g)\ket0$ のノルム二乗は $2+4|c|^2$。実際、$g=c f+\sqrt{1-|c|^2}\,h$ と直交成分に分ければ、三個の $f$ の項のノルム重みは $6|c|^2$、二個の $f$ と一個の $h$ の重みは $2(1-|c|^2)$。$|c|^2=1/4$ では和は $3$、規格化後に三個とも $f$ と測る確率は $(6/4)/3=1/2$ である。
:::

### 状態の数と生成演算子の順序

$m$ 個の直交モードを使い、総粒子数が $n$ の部分空間を考える。占有数 $n_i\ge0$、$\sum_i n_i=n$ が一つの正規直交基底を決める。座標にモードを順に割り当てる数ではなく、この整数方程式を数える。$n$ 個の印と $m-1$ 個の仕切りを横一列に並べると、各区画の印の数が $n_i$ になる。従って次元は $\binom{n+m-1}{m-1}$。これは対称テンソルの次元と一致する。生成演算子の順番を変えても可換であり、どのモードがいくつかだけが基底指標になる。

:::example[三モードに二粒子を置く]
$(2,0,0),(0,2,0),(0,0,2)$ と、$(1,1,0),(1,0,1),(0,1,1)$ の六通り。前半の基底には $(a_i^\dagger)^2/\sqrt2$、後半には $a_i^\dagger a_j^\dagger$ を使う。例えばモードのエネルギーが $1,2,3\,\mathrm{GeV}$ なら、二粒子で全エネルギー $4\,\mathrm{GeV}$ の基底は $(0,2,0)$ と $(1,0,1)$ の二つである。同じエネルギーでも占有数は異なり、番号演算子を測れば区別できる。
:::

## 粒子解釈とエネルギー運動量演算子

### Hamiltonian の対角化

$H=\int\dd^3x\,[\pi^2+|\nabla\phi|^2+m^2\phi^2]/2$ にモードを代入する。$aa$ 項は空間積分で $\vb q=-\vb p$ になり、その係数は $-E_{\vb p}^2+\vb p^2+m^2=0$ である。$a^\dagger a$ 項では $\vb q=\vb p$ で係数は $2E_{\vb p}^2$ になる。従って箱では $H=\sum_{\vb p}E_{\vb p}(a^\dagger_{\vb p}a_{\vb p}+1/2)$ である。:eqref[eq-qft-canonical-tensor] の空間成分では異なる順の積が $\vb p$ と $-\vb p$ で反対符号の真空項を与える。対称なカットオフを使うと
$$
:P^\mu:=\int\dd\Pi_p\,p^\mu a^\dagger(\vb p)a(\vb p)
\label{eq-qft-scalar-momentum}
$$
となる。波束生成演算子との交換子は $[P^\mu,a^\dagger(f)]=a^\dagger(p^\mu f)$ であり、狭い波束のエネルギー・運動量は中心 $p^\mu$ に近づく。例えば $m=3\,\mathrm{GeV}$、$p_z=4\,\mathrm{GeV}$ のモードに二量子あれば、その寄与は $(10,0,0,8)\,\mathrm{GeV}$ である。

### 局在した源が作る自由粒子

自由場に古典源を加え、$H_I(t)=-\int J(t,\vb x)\phi(t,\vb x)\dd^3x$ とする。$J$ は実の滑らかな関数とし、時空にコンパクトな台をもつと仮定する。この条件なら正エネルギー質量殻に制限した Fourier 変換も一粒子測度について二乗可積分である。$[H_I(t),H_I(s)]$ は :eqref[eq-qft-scalar-commutator] により数なので、積の時間順序を通常の指数へ移すと数の位相だけが出る。これは各時間片の指数に Baker–Campbell–Hausdorff 公式を適用して得られる。残る指数は
$$
\int\dd\Pi_p\,(\alpha_p a^\dagger(p)-\alpha_p^*a(p)),\qquad
\alpha_p=i\int\dd^4x\,e^{ipx}J(x).
$$
消滅演算子を右へ移すと、真空に作用した状態は $e^{-\|\alpha\|^2/2}\exp(a^\dagger(\alpha))\ket0$ である。従って $\expval{N}=\int\dd\Pi_p|\alpha_p|^2$、$\expval{H}=\int\dd\Pi_p E_p|\alpha_p|^2$。源の Fourier 変換のうち、実粒子に変わるのは正の質量殻にある成分である。

:::example[一モードの粒子数の分布]
規格化した一モードで $|\alpha|=2$ とすると、$n$ 粒子の係数は $e^{-2}\alpha^n/\sqrt{n!}$。従って粒子数の分布は $e^{-4}4^n/n!$ で、平均 $4$ である。モードの角周波数が $3\,\mathrm{GeV}$ なら正規順序化したエネルギーの平均は $12\,\mathrm{GeV}$ である。
:::

### 四元運動量が場の並進を実装すること

正規順序化した $P^\mu=\int\dd\Pi_p\,p^\mu a^\dagger(p)a(p)$ に対し :eqref[eq-qft-scalar-modes] を用いる。例えば積の交換子 $[a^\dagger(p)a(p),a(k)]=[a^\dagger(p),a(k)]a(p)$ によって $[P^\mu,a(k)]=-k^\mu a(k)$。生成演算子では符号が逆になる。これらを平面波展開へ入れると
$$
[P^\mu,\phi(x)]=-i\partial^\mu\phi(x).
$$
実際、$-i\partial^\mu$ は正周波数の $e^{-ipx}$ に $-p^\mu$、負周波数の $e^{ipx}$ に $+p^\mu$ を掛ける。従って $e^{ia_\mu P^\mu}\phi(x)e^{-ia_\mu P^\mu}$ を $a$ で微分すれば $\phi(x+a)$ の方向微分と一致する。零点エネルギーの定数を加えても随伴作用では相殺される。粒子状態の並進の位相と場の引数の移動で正負の指数を混同しないよう、ここでは随伴作用の順を明記した。

:::example[二粒子波束の全エネルギーと不変質量]
二つの箱モードに一粒子ずつ置き、両方の質量を $3\,\mathrm{GeV}$、運動量の大きさを $4\,\mathrm{GeV}$ とする。平行な場合は全ベクトルが $(10,0,0,8)\,\mathrm{GeV}$ なので不変質量 $6\,\mathrm{GeV}$。逆向きなら全ベクトルは $(10,\vb0)\,\mathrm{GeV}$ で不変質量 $10\,\mathrm{GeV}$。粒子数と実験室の全エネルギーは同じでも、多粒子状態の $P^2$ は運動量の配置によって違う。散乱の閾値では全不変質量を使う。
:::

### 空間運動量をモードから直接求める

実場の正準テンソルから $P^i=\int\dd^3x\,\pi\,\partial^i\phi$ と書く。箱の展開を使うと $aa$ の項は空間デルタで $\vb q=-\vb p$、両周波数が一致する。このとき係数は $p^i$ に比例して奇なので $\vb p$ と $-\vb p$ を足せば消える。$a^\dagger a^\dagger$ も同様。混合項では $\vb q=\vb p$ が選ばれ、$P^i=\sum_{\vb p}p^i(a^\dagger_{\vb p}a_{\vb p}+a_{\vb p}a^\dagger_{\vb p})/2$。交換関係で $N_{\vb p}+1/2$ となる。対称な正則化では $\sum_{\vb p}p^i/2=0$ なので、エネルギーの零点項とは違い全運動量の定数は零である。正則化前の奇の和を無条件に零と置くより、$\pm\vb p$ を同時に入れることを先に指定する。

:::example[定在波は運動量の分散を持つ]
余弦の一量子では二つの進行波が同じ重みをもつ。例えば波数 $\vb k=(0,0,4)\,\mathrm{GeV}$ では $P^z$ の平均は零、二乗の平均は $16\,\mathrm{GeV}^2$。二つの異なる進行波に一粒子ずつ入れれば全運動量演算子で厳密に零となり、その分散も零である。両者は「左右の運動量を含む」という説明だけでは同じに見えるが、粒子数と運動量の分布を指定すると区別できる。
:::

## 正規順序と零点エネルギー

:term[正規順序]{en="normal ordering" yomi="せいきじゅんじょ"} $:\!X\!:$ は、モードの単項式で生成演算子を左、消滅演算子を右に配置する操作である。従って $:aa^\dagger:=a^\dagger a$ であり、元の積との差は交換子である。正規順序化した単項式の自由真空期待値は、場の因子が一個以上あれば零である。二次の場の積では、元の積から自由真空期待値を引く操作と一致する。四次以上では、真空定数に加えて少ない場を持つ項も差し引く。例えば $:\phi^2:=\phi^2-\mel{0}{\phi^2}{0}$ は正則化した演算子での式である。

質量ゼロで $|\vb p|\le\Lambda$ の三運動量カットオフを入れると単位体積あたりの零点エネルギーは
$$
\rho_0=\frac12\int_{|\vb p|\le\Lambda}\frac{\dd^3p}{(2\pi)^3}|\vb p|
       =\frac{\Lambda^4}{16\pi^2}.
$$
例えば $\Lambda=1\,\mathrm{GeV}$ なら約 $6.3\times10^{-3}\,\mathrm{GeV}^4$。このカットオフは計算をした静止系を特別扱いする。非重力の時間発展では $H\mapsto H+C$ の効果が全状態の共通位相なので差し引ける。異なる境界条件間のエネルギー差や、重力に結合する絶対値を議論するときは、許される減算を別に指定する必要がある。

### Casimir エネルギーの有限な差

一次元の区間 $0<x<L$ の質量ゼロ場に Dirichlet 条件を課すと、固有関数 $\sin(n\pi x/L)$ の周波数は $n\pi/L$ である。長さの単位をもつ $\epsilon>0$ を使って真空和を減衰させると
$$
E_\epsilon=\frac{\pi}{2L}\sum_{n=1}^{\infty}n e^{-n\pi\epsilon/L}
=\frac{\pi}{2L}\frac{e^{-h}}{(1-e^{-h})^2},\qquad h=\pi\epsilon/L .
$$
幾何級数を $e^{-h}$ で微分した。$e^{-h}/(1-e^{-h})^2=1/[4\sinh^2(h/2)]=h^{-2}-1/12+\order{h^2}$。一方で同じ減衰をかけた自由空間の密度は $(2\pi)^{-1}\int_0^\infty k e^{-\epsilon k}\dd k=1/(2\pi\epsilon^2)$。区間長を掛けて引くと有限部分は $-\pi/(24L)$。$L$ に依存しない端のエネルギーを除いた力は $-\partial_LE=-\pi/(24L^2)$ である。ここで計算したのは一次元模型であり、三次元の導体板では横運動量と電磁場の偏極を含める。

### 真空のエネルギー密度と圧力を同じ正則化で求める

:term[真空エネルギー]{en="vacuum energy" yomi="しんくうえねるぎー"} の扱いは、$H$ の加法定数と局所テンソルで区別する。実スカラーの最小結合テンソルは $\Theta^{\mu\nu}=\partial^\mu\phi\partial^\nu\phi-\eta^{\mu\nu}[(\partial\phi)^2-m^2\phi^2]/2$。各縮約に同じ空間カットオフ $\Lambda$ を置けば
$$
\left\langle\partial^\mu\phi\,\partial^\nu\phi\right\rangle_\Lambda
=\int_{|\vb p|<\Lambda}\frac{\dd^3p}{(2\pi)^3 2E}\,p^\mu p^\nu,\qquad
\langle\phi^2\rangle_\Lambda=\int_{|\vb p|<\Lambda}\frac{\dd^3p}{(2\pi)^3 2E}.
$$
質量殻 $p^2=m^2$ なので密度の括弧は縮約後に零であり、$\Theta^{00}$ は上の $p^0p^0$ の積。回転対称なため空間ブロックの非対角は消え、$\langle\Theta^{ij}\rangle=P_\Lambda\delta^{ij}$ となる。
$$
\rho_\Lambda=\frac1{4\pi^2}\int_0^\Lambda\dd p\,p^2\sqrt{p^2+m^2},\qquad
P_\Lambda=\frac1{12\pi^2}\int_0^\Lambda\dd p\,\frac{p^4}{\sqrt{p^2+m^2}} .
$$
例えば $m=0$ では $P_\Lambda=\rho_\Lambda/3$ であり、切り取った正周波数モードの方向平均から出た式である。

これに対し、Lorentz 共変にくりこんだ真空テンソルにおいて不変な二階対称テンソルは $\eta^{\mu\nu}$ の定数倍に限られる。空間回転で $T^{0i}=0$、$T^{ij}=P\delta^{ij}$ とし、$z$ 方向のブーストを掛けると $T'^{0z}=(\rho+P)\sinh\chi\cosh\chi$。全ての $\chi$ で零であるには $P=-\rho$ が必要である。空間カットオフでは異なる慣性系で選ぶモードの集合が違うため、この議論の不変性が仮定できない。従って $\rho_\Lambda$ だけでなく、正則化と局所対項でテンソル全体がどう変わるかを指定してから真空の応力を議論する。

:::example[カットオフ付きの質量ゼロ真空]
$\Lambda=2\,\mathrm{GeV}$ なら $\rho_\Lambda=1/\pi^2\,\mathrm{GeV}^4$、$P_\Lambda=1/(3\pi^2)\,\mathrm{GeV}^4$。$v=3/5$ のブーストを、同じテンソル成分の数値に作用させれば $T'^{0z}=(4\rho_\Lambda/3)(5/4)(3/4)=5\rho_\Lambda/4$。新しい系で球カットオフを最初から課すと回転対称性により $T^{0z}=0$ である。二つの手順の差は、運動量空間の球がブーストで別の領域へ移ることに対応する。
:::

### 周期境界と固定端の Casimir 力を比べる

一周 $L$ の円周で質量ゼロ場を取り、零モードを別に扱う。非零モードには $\pm2\pi n/L$、$n\ge1$ の二つの進行波がある。従って減衰を入れた和は
$$
E_\epsilon^{\rm per}=\frac{2\pi}{L}\sum_{n\ge1}n e^{-2\pi n\epsilon/L}
  =\frac{L}{2\pi\epsilon^2}-\frac{\pi}{6L}+\order{\epsilon^2/L^3}.
$$
二つの進行波と各振動子の $1/2$ が相殺して最初の係数を与える。無限直線の密度の減算は固定端の場合と同じ $L/(2\pi\epsilon^2)$。残りの有限部分は $-\pi/(6L)$、長さに対する力は $-\pi/(6L^2)$ である。例えば $L=2\,\mathrm{GeV}^{-1}$ なら周期境界の有限エネルギーは $-\pi/12\,\mathrm{GeV}$、同じ $L$ の固定端は $-\pi/48\,\mathrm{GeV}$。主発散が同じでも有限な力は周波数間隔と多重度によって変わる。コンパクトな位相場などでは零モードや巻き数状態も変わるので、周期条件だけから同じ全スペクトルとみなさない。

### 真空を変えると正規順序の定数も変わる

一モードで座標 $q$ と運動量 $p$ を固定し、基準周波数を $\omega$ とする。真空における $\langle q^2\rangle_\omega=1/(2\omega)$ を使うと $:q^2:_\omega=q^2-1/(2\omega)$。同じ物理演算子 $q^2$ を別の周波数 $\Omega$ の基準で正規順序化すれば
$$
:q^2:_\omega=:q^2:_\Omega+\frac1{2\Omega}-\frac1{2\omega}.
$$
四次では全ての縮約を数えると $:q^4:_\omega=q^4-6D_\omega q^2+3D_\omega^2$、$D_\omega=1/(2\omega)$。$\Omega$ での式を代入して
$$
:q^4:_\omega=:q^4:_\Omega
+6(D_\Omega-D_\omega):q^2:_\Omega+3(D_\Omega-D_\omega)^2.
$$
違う基準を選ぶと四次から二次の局所項も変わる。二つの規約の係数を同じままにすると、同じ Hamiltonian を書いたことにはならない。

:::example[小さな四次項が二つの規約で作る質量項]
周波数 $\omega=1,\Omega=2$ を同じ単位で取り、$\lambda:q^4:_\omega/24$ を $\Omega$ の規約に移す。$D_\Omega-D_\omega=1/4-1/2=-1/4$ なので二次の係数は $-\lambda/16$、定数は $\lambda/128$。従って新しい正規順序の四次項だけを使うと、元の演算子が持っていた二次と定数を失う。場の理論では $D$ を同一点の二点分布の正則化値で置き、くりこまれた質量や真空エネルギーの条件と合わせる。
:::

## 複素スカラー場と反粒子

複素場では $\pi_\varphi=\dot\varphi^\dagger$、$\pi_{\varphi^\dagger}=\dot\varphi$ である。独立な二実場を :eqref[eq-qft-scalar-expansion] で展開し、二つの消滅演算子の和と差を $a=(a_1+ia_2)/\sqrt2$、$b=(a_1-ia_2)/\sqrt2$ と置くと
$$
\varphi(x)=\int\dd\Pi_p\,[a(\vb p)e^{-ipx}+b^\dagger(\vb p)e^{ipx}]
\label{eq-qft-complex-scalar}
$$
となる。$a,b$ はそれぞれ :eqref[eq-qft-scalar-modes] を満たし、交差交換子は零である。古典電荷 :eqref[eq-qft-uone-current] へ代入したとき $a^\dagger b^\dagger$ と $ba$ は空間積分で同じエネルギーの差を持つため消える。正規順序で $Q=\int\dd\Pi_p(a^\dagger a-b^\dagger b)$ となる。従って $[Q,\varphi]=-\varphi$、$e^{i\alpha Q}\varphi e^{-i\alpha Q}=e^{-i\alpha}\varphi$ である。$b^\dagger$ が作る :term[反粒子]{en="antiparticle" yomi="はんりゅうし"} は負電荷であり、Hamiltonian には $a^\dagger$ と同じ正エネルギーを寄与する。

:::example[電荷を保つ局所的な対生成]
正規順序積 $:\varphi^\dagger\varphi:$ の中には $a^\dagger b^\dagger$ があり、真空から電荷の反対な二粒子を生成する。$[Q,a^\dagger]=a^\dagger$、$[Q,b^\dagger]=-b^\dagger$ だから $[Q,a^\dagger b^\dagger]=0$ である。例えば箱で同じ絶対値 $4\,\mathrm{GeV}$ の逆向き運動量を持つ質量 $3\,\mathrm{GeV}$ の対では全運動量は零、全エネルギーは $10\,\mathrm{GeV}$ となる。
:::

### 電荷の符号から複素場の二種類の相関を読む

箱の粒子を $a^\dagger$、反粒子を $b^\dagger$ が作るとし、真空は $Q\ket0=0$。演算子 $O$ が $[Q,O]=q_O O$ を満たせば、$\bra0[Q,O]\ket0=0$ だから $q_O\ne0$ の期待値は零である。例えば $[Q,\varphi]=-\varphi$ から $[\!Q,\varphi(x)\varphi(y)\!]=-2\varphi(x)\varphi(y)$。従って $\langle T\varphi(x)\varphi(y)\rangle=0$。一方 $\varphi(x)\varphi^\dagger(y)$ は電荷零であり、$t_x>t_y$ で残る縮約は $a(p)a^\dagger(q)$、反対の順序では $b(p)b^\dagger(q)$ の縮約である。時間差が正なら粒子、負なら反粒子の周波数が正の測度から出て、全体は実場と同じ $\Delta_F(x-y)$ になる。複素場の反粒子は負エネルギーを伝える状態ではなく、同じ正エネルギーで反対の電荷を持つ。

:::example[Hermite な挿入と荷電した挿入の差]
$\varphi_1=(\varphi+\varphi^\dagger)/\sqrt2$、$\varphi_2=(\varphi-\varphi^\dagger)/(i\sqrt2)$ とする。電荷を持つ二つの縮約が零で、電気的に中性の二つが等しいから $\langle T\varphi_1\varphi_1\rangle=\Delta_F$、$\langle T\varphi_1\varphi_2\rangle=0$。従って複素場を二実場とみなしても、粒子と反粒子に分けても同じ自由相関が得られる。ただし $\varphi_1$ が真空へ作用した一量子は二電荷の重ね合わせで、$\varphi$ が作用した一量子は反粒子の電荷固有状態である。
:::

### 電荷を明示的に破ると粒子の基底が変わる

複素場の自由ポテンシャルに $\mu^2(\varphi^2+\varphi^{*2})/2$ を足す。実場へ分けるとポテンシャルは $(m^2+\mu^2)\varphi_1^2/2+(m^2-\mu^2)\varphi_2^2/2$ となる。$m^2>|\mu^2|$ として各二乗を正に取ると、それぞれ質量 $m_1,m_2$ の実場を独立に量子化できる。もとの位相回転は二つを混ぜるから、$\mu^2\ne0$ でその電荷は Hamiltonian と可換でなくなる。異なる質量の二つの伝播関数を $\Delta_1,\Delta_2$ とすれば
$$
\langle T\varphi(x)\varphi^\dagger(y)\rangle=\tfrac12(\Delta_1+\Delta_2),\qquad
\langle T\varphi(x)\varphi(y)\rangle=\tfrac12(\Delta_1-\Delta_2).
$$
後者の存在は電荷を持つ二つの挿入が相関を持つことを示し、電荷が保存する場合には $m_1=m_2$ で消える。例えば $m^2=5\,\mathrm{GeV}^2$、$\mu^2=4\,\mathrm{GeV}^2$ なら二質量は $3,1\,\mathrm{GeV}$。荷電相関の二極の留数は $i/2,-i/2$ であり、正のスペクトル表示を課す自己随伴の二点とは異なる。

### 化学ポテンシャルは反粒子に逆符号で作用する

正規順序化した $H$ と $Q$ で、熱平衡の重みを $e^{-\beta(H-\mu Q)}$ とする。箱の一つの運動量で $H=E(N_a+N_b)$、$Q=N_a-N_b$ だから、粒子と反粒子の重みは $e^{-\beta(E-\mu)N_a}$ と $e^{-\beta(E+\mu)N_b}$。幾何級数を足すとモードの分配関数は
$$
Z_{\vb p}=\frac1{(1-e^{-\beta(E-\mu)})(1-e^{-\beta(E+\mu)})}.
$$
これは $E\pm\mu>0$ の両条件で収束する。$\partial_\mu\log Z_{\vb p}/\beta$ から
$$
\langle Q\rangle_{\vb p}
=\frac1{e^{\beta(E-\mu)}-1}
-\frac1{e^{\beta(E+\mu)}-1}.
$$
従って最小の $E$ が $m$ なら、有限箱の全モードを含む通常のトレースには $|\mu|<m$ が必要である。$\mu=m$ で零運動量粒子の重みが一定になり、任意の占有数を足して発散する。有限の総電荷を固定した計算では、深い赤外の占有数を熱積分と別に扱う必要がある。

:::example[一モードにおける電荷の熱平均]
温度とエネルギーを $\beta E=\frac32\log2$、化学ポテンシャルを $\beta\mu=\frac12\log2$ に選ぶ。粒子の指数は $e^{\beta(E-\mu)}=2$、反粒子は $e^{\beta(E+\mu)}=4$。平均占有数は $1,1/3$ なので、電荷の平均は $2/3$ である。$\mu$ の符号を反転するとこの平均は $-2/3$ に変わる一方、全粒子数の平均の $4/3$ は同じである。
:::

## 因果律と微視的因果性

### 交換子で因果性を調べる

局所可観測量を空間的な間隔に置いたとき、測定の順番は座標の選択で逆になる。そこで空間的分離で局所演算子が可換という条件を :term[微視的因果性]{en="microcausality" yomi="びしてきいんがせい"} と呼ぶ。自由実場では :eqref[eq-qft-scalar-modes] によって
$$
C(x)=[\phi(x),\phi(0)]=\int\dd\Pi_p\,(e^{-ipx}-e^{ipx})
\label{eq-qft-scalar-commutator}
$$
となる。測度は Lorentz 不変なので、$x^2<0$ なら $x^0=0$ に移せる。その系で $\vb p\mapsto-\vb p$ が二項を交換するため零である。この議論は分布の台が閉じた光円錐に含まれることを述べる。二点期待値の各項は空間的分離でも非零になり得る。例えば質量ゼロで $r>0$ の等時刻に $e^{-\epsilon|\vb p|}$ を入れると、角積分に続いて $\int_0^\infty\dd p\,e^{-\epsilon p}\sin(pr)=r/(r^2+\epsilon^2)$ となるので、期待値は $\epsilon\downarrow0$ で $1/(4\pi^2r^2)$ になる。

### 真空相関の空間的な減衰

自由な質量付き場で空間的な間隔を等時刻へ移し、距離を $r>0$ とする。等時刻では二点関数の空間 Fourier 係数が $1/(2E_{\vb p})$ である。これは $\int\dd p_4/(2\pi)(p_4^2+E_{\vb p}^2)^{-1}$ に等しい。従って等時刻の相関を $1/(p_E^2+m^2)$ の四次元 Euclid Fourier 積分として計算できる。逆数を $\int_0^\infty\dd s\,e^{-s(p_E^2+m^2)}$ として Gaussian 積分を実行すると
$$
D(r)=\frac{1}{16\pi^2}\int_0^\infty\frac{\dd s}{s^2}
     \exp\left[-m^2s-\frac{r^2}{4s}\right]
=\frac{m}{4\pi^2r}K_1(mr) .
$$
最後は $s=(r/(2m))e^u$ を入れ、積分の偶部分を取ると $K_1(z)=\int_0^\infty e^{-z\cosh u}\cosh u\,\dd u$ に一致する。$m=0$ では $r^2/(4s)$ を変数にして $1/(4\pi^2r^2)$ を得る。$mr\gg1$ では $\cosh u\simeq1+u^2/2$ により $K_1(mr)\simeq\sqrt{\pi/(2mr)}e^{-mr}$。従って質量の逆数が相関の減衰長を与える。一方、交換子が空間的な間隔で零であることは質量によらない。

### 局所多項式と応答

空間的に分離した二領域に台をもつ試験関数 $f,g$ について、場を $\phi(f)=\int f(x)\phi(x)\dd^4x$ として塗りならすと $[\phi(f),\phi(g)]=\int f(x)g(y)C(x-y)\dd^4x\dd^4y=0$ である。積の交換子 $[AB,C]=A[B,C]+[A,C]B$ を反復すると、各領域の塗りならした場の多項式も可換である。例えば 実の $f$ に対して局所ユニタリ $U=e^{-i\lambda\phi(f)}$ を作用させても $\expval{U^\dagger\phi(g)U}=\expval{\phi(g)}$。二領域の真空期待値の相関が非零でも、一方の局所操作では他方の観測分布を変えられない。

### 非局所な正周波数部分と局所性

場を $\phi=\phi^{(+)}+\phi^{(-)}$ に分ける。正周波数は消滅部だけ、負周波数は生成部だけである。従って $[\phi^{(+)}(x),\phi^{(-)}(y)]=D^+(x-y)$ で、$D^+(z)=\int\dd\Pi_p e^{-ipz}$。$[\phi^{(-)}(x),\phi^{(+)}(y)]=-D^+(y-x)$。空間的分離で二つの $D^+$ が同じ値をもつため、全交換子は消える。一方、$\phi^{(+)}$ は単独では局所な可観測量ではない。時刻固定の初期値から得るには
$$
\phi^{(+)}(0,\vb x)
=\tfrac12\left[\phi(0,\vb x)+i(-\nabla^2+m^2)^{-1/2}\pi(0,\vb x)\right]
$$
を使う。これはモードごとに $a$ の係数を解いた式である。平方根の逆はフーリエ空間で非多項式の $1/E_{\vb p}$ を掛けるため、実空間では微分の有限和で表せず、離れた初期値も読む。正周波数だけで信号を表そうとすると、この射影の非局所性を持ち込む。

:::example[同時刻での個別の縮約]
質量ゼロの $t=0$, $r=2\,\mathrm{GeV}^{-1}$ で $D^+(0,r)=1/(16\pi^2)\,\mathrm{GeV}^2$。正周波数部分と負周波数部分の交換子はこの非零値を取るが、引数を逆にしたときも同じなので全場の交換子は零である。二点関数の減衰距離を使って測定の因果性を判定すると、この差が分からない。測定の順序を変える操作は全場の交換子に対応するためである。
:::

### 自由な交換子を初期値問題として求める

$C(t,\vb x)=[\phi(t,\vb x),\phi(0,\vb0)]$ は先頭の場が Klein–Gordon 方程式を満たすので、同じ微分方程式を満たす。等時刻で $C(0,\vb x)=0$、$\partial_t C(0,\vb x)=-i\delta^3(\vb x)$。$\vb x$ を Fourier 変換すると各モードの一般解は $A\cos Et+B\sin Et$、初期値から $A=0,B=-i/E$。したがって
$$
C(t,\vb x)=-i\int\frac{\dd^3p}{(2\pi)^3}\frac{\sin(E t)}{E}e^{i\vb p\cdot\vb x}.
$$
この式は二つの質量殻指数の差から得た :eqref[eq-qft-scalar-commutator] と一致する。局在した初期値に対する領域依存性は古典場と同じであり、デルタを小球に台をもつ平滑関数で近似し、分布極限を取ると閉光円錐の外で零になる。この証明は、ローレンツ変換で等時刻へ移す証明とは違い、初期値の伝播に基づく。境界がある場で Lorentz 対称性が破れても、局所微分方程式の領域依存性が保たれることがある。

:::example[一次元の質量ゼロでは円錐内部に尾が残る]
時間と空間が各一成分で $\Box=\partial_t^2-\partial_x^2$ とする。$C$ の時刻零のデータを上と同じに取れば、D'Alembert の式で $t>0$ の解は $C(t,x)=-(i/2)\int_{x-t}^{x+t}\delta(y)\dd y$。従って $|x|<t$ で $-i/2$、外で零。遅延核は円錐内部に一定値 $1/2$ を持つ。四次元の質量ゼロでは境界上のデルタだけが残るので、光速で信号の先頭が届くという結論は同じでも、その後の応答は時空次元で異なる。
:::

### 質量ゼロの二点分布を光円錐の近くで読む

$t\mapsto t-i\epsilon$、$\epsilon>0$ として正周波数積分に収束を付ける。$m=0,r>0$ で角積分を行えば
$$
D_\epsilon^+(t,r)=\frac1{4\pi^2 r}
\int_0^\infty\dd p\,e^{-(\epsilon+it)p}\sin(pr)
=\frac{1}{4\pi^2[r^2-(t-i\epsilon)^2]} .
$$
最後の等式は $\sin(pr)=(e^{ipr}-e^{-ipr})/(2i)$ とし、二つの指数を $1/(\epsilon+i(t\mp r))$ に積分して引けば得る。空間的分離では $r^2-t^2>0$ なので極限は実関数である。時間的分離でも光円錐から離れれば有限だが、$t=\pm r$ の境界には分布の処方が残る。例えば $t>0$ の光円錐の近くでは分母の虚部は正、$t<0$ では負である。従って二点の差から光円錐上のデルタの係数が決まる。

:::example[質量ゼロの交換子の光円錐項]
$r>0$ で $s=t^2-r^2$ とすると、$D^+(t,r)-D^+(-t,r)$ の主値部分は相殺し、デルタは $-i\,\sgn(t)\delta(s)/(2\pi)$ になる。従って $i\theta(t)C=\theta(t)\delta(s)/(2\pi)$。$t>0$ の根で $\delta(s)=\delta(t-r)/(2r)$ となるため、遅延核は $\delta(t-r)/(4\pi r)$。実の真空相関の主値部分はそのまま応答の尾として残らない。自由場で因果応答を作るのは正と負の周波数の差である。
:::

## Feynman 伝播関数と遅延 Green 関数

### 極の処方から境界条件を読む

:term[Feynman 伝播関数]{en="Feynman propagator" yomi="ふぁいんまんでんぱかんすう"} を真空の時間順序積 $\Delta_F(x)=\mel{0}{T\phi(x)\phi(0)}{0}$ と定める。$t>0$ では正周波数、$t<0$ では負周波数が残る。$\epsilon>0$ に対し二極をそれぞれ下半面・上半面にずらし、その極限を記号 $+i0$ で書けば
$$
\Delta_F(x)=\int\frac{\dd^4p}{(2\pi)^4}\frac{i\,e^{-ipx}}{p^2-m^2+i0},\qquad
(\Box+m^2)\Delta_F=-i\delta^4(x).
\label{eq-qft-scalar-feynman}
$$
$t>0$ は下半面を時計回りに閉じ、留数 $i/(2E_{\vb p})$ に $-i$ を掛ける。$t<0$ は上半面を反時計回りに閉じ、留数 $-i/(2E_{\vb p})$ に $i$ を掛ける。いずれも時間順序積と一致する。

一方、外部源に対する応答には :term[遅延 Green 関数]{en="retarded Green function" yomi="ちえんぐりーんかんすう"} $G_R=i\theta(t)C(x)$ を用いる。空間 Fourier 係数は $\theta(t)\sin(E_{\vb p}t)/E_{\vb p}$ である。二階微分で $\theta$ の境界から $\delta(t)$ が出るので $(\Box+m^2)G_R=\delta^4(x)$ となる。両極を下半面に置く Fourier 表現は $-1/[(p^0+i0)^2-E_{\vb p}^2]$ である。質量ゼロでは角積分の後に $\int_0^\infty\sin(pr)\sin(pt)\dd p=\pi[\delta(t-r)-\delta(t+r)]/2$ を使い、$r>0$ で $G_R=\delta(t-r)/(4\pi r)$ を得る。

### 質量によって光円錐の内側に残る応答

$r=|\vb x|$、$s=t^2-r^2$ とする。四次元の遅延解を実空間で書くと
$$
G_R(t,\vb x)=\theta(t)\left[\frac{\delta(s)}{2\pi}
-\theta(s)\frac{m J_1(m\sqrt{s})}{4\pi\sqrt{s}}\right]
\qquad (r>0).
$$
$J_1$ は原点で $J_1(z)=z/2+\mathcal O(z^3)$ の Bessel 解である。この式の検算は $t>0$ で光円錐を横切る分布を調べればよい。$f=f(s)$ に対し $\Box f=4s f''+8f'$ だから、内側の $f$ は $4s f''+8f'+m^2f=0$ を満たす。$u=\sqrt s$ と置くと $f=h(mu)/u$ は Bessel 方程式 $z^2h''+zh'+(z^2-1)h=0$ に従う。$\theta(s)f(s)$ の $\Box$ には $4f(0)\delta(s)$ が残る。実際、$s\delta'= -\delta$ を使って $4s$ と $8$ の項を足すとこの係数になる。$f(0)=-m^2/(8\pi)$ なので、この項が光円錐上のデルタの質量項 $m^2\delta(s)/(2\pi)$ を消す。時刻零の跳びは質量ゼロと同じ光円錐項から生じ、内側の有限密度の項は積分領域が $t^3$ で縮むため跳びを変えない。従って質量があれば衝撃の後も光円錐の内側に応答が残る。

### 遅延解の一意性

源と初期値が与えられた二つの解の差 $u$ は自由方程式を満たし、源が入る前の $u,\dot u$ が零である。空間積分した $E_u=\int(\dot u^2+|\nabla u|^2+m^2u^2)\dd^3x/2$ の時間微分は、部分積分で $\int\dot u(\Box+m^2)u\,\dd^3x$ となる。無限遠の境界項を消せば $E_u=0$ が保たれ、零の初期条件から $u=0$ が従う。従って同じ微分演算子の逆でも、Feynman の境界条件を指定した解を遅延応答に使うことはできない。

:::example[有限時間だけ作用する源]
質量ゼロで原点の源を $J(t,\vb x)=J_0\delta^3(\vb x)[\theta(t)-\theta(t-T)]$ とする。遅延核を畳み込むと、$r>0$ で $\phi=J_0[\theta(t-r)-\theta(t-r-T)]/(4\pi r)$ である。$r=2\,\mathrm{GeV}^{-1}$、$T=3\,\mathrm{GeV}^{-1}$ では応答がある時間域は $2<t<5$ である。質量が非零なら光円錐内側の項も畳み込むため、$t>r+T$ にも尾が残り得る。
:::

### 先進・遅延・時間順序の核を分布として比較する

$E=\sqrt{|\vb p|^2+m^2}>0$ とし、時間 Fourier 変換の逆を $\int\dd p^0\,e^{-ip^0t}/(2\pi)$ とする。:term[先進 Green 関数]{en="advanced Green function" yomi="せんしんぐりーんかんすう"} を $G_A=-i\theta(-t)C(x)$ と定める。各空間モードは $-\theta(-t)\sin(Et)/E$ であり、時刻零の右側と左側の一階微分の差は $1$。従って $(\Box+m^2)G_A=\delta^4(x)$ である。$G_A$ は過去の光円錐に台を持ち、源の時刻より後に応答がない。遅延核では両極が下半面、先進核では両極が上半面である。
$$
G_R(p)=-\frac1{(p^0+i0)^2-E^2},\qquad
G_A(p)=-\frac1{(p^0-i0)^2-E^2}.
$$
同じ演算子を掛けると右辺は同じになるため、両者の差は自由方程式を満たす。実際、$G_R-G_A=iC$。逆核を定めるには微分方程式だけでなく境界条件が必要である。

主値とデルタを分けるには $1/(u+i0)=\operatorname{PV}(1/u)-i\pi\delta(u)$ を使う。これを示すには $1/(u+i\epsilon)=u/(u^2+\epsilon^2)-i\epsilon/(u^2+\epsilon^2)$ と書く。滑らかな試験関数の原点値を引いた実部は対称極限で主値に移る。虚部は $u=\epsilon v$ と置けば重み $1/(1+v^2)$ の積分が $\pi$ となり、試験関数の原点値を選ぶ。これで
$$
\begin{aligned}
G_R(p)&=-\operatorname{PV}\frac1{p^2-m^2}
+i\pi\sgn(p^0)\delta(p^2-m^2),\\
G_A(p)&=-\operatorname{PV}\frac1{p^2-m^2}
-i\pi\sgn(p^0)\delta(p^2-m^2),\\
\Delta_F(p)&=i\operatorname{PV}\frac1{p^2-m^2}+\pi\delta(p^2-m^2)
\end{aligned}
$$
となる。遅延・先進式の $\sgn(p^0)$ は、正・負エネルギー極で $2p^0 i0$ の向きが異なることを示す。これに対し時間順序は $p^2-m^2+i0$ とするので、どちらの枝でもデルタの係数は同符号である。例えば $\Delta_F=-iG_R+2\pi\theta(-p^0)\delta(p^2-m^2)$。負エネルギー殻の項を足すことで、遅延の境界条件から真空の時間順序の境界条件へ移る。

:::example[同じ瞬時入力に対する二つの応答]
周波数 $E=2\,\mathrm{GeV}$ の一モードで $J=A\delta(t)$ とする。遅延応答は $A\theta(t)\sin(Et)/E$、先進応答は $-A\theta(-t)\sin(Et)/E$。例えば $t=\pi/(4\,\mathrm{GeV})$ では遅延が $A/(2\,\mathrm{GeV})$ で先進は零。$t=-\pi/(4\,\mathrm{GeV})$ では先進が同じ値で遅延が零である。差の $A\sin(Et)/E$ は $t=0$ で一階微分の跳びを持たず、同次方程式を満たす。
:::

:::example[時間順序核の一階微分の跳び]
空間運動量を固定すると $\Delta_F(t,\vb p)=e^{-iE|t|}/(2E)$。左右の一階微分は $-i/2,+i/2$ なので、差は $-i$。従って $(\partial_t^2+E^2)\Delta_F=-i\delta(t)$。時間順序核を直接実の源と畳み込むと、この右辺と境界条件の両方が遅延応答とは違う。例えば一モードで $\Delta_F(0)=1/(2E)$ に対し $G_R(0)$ の連続な値は零である。
:::

### 低周波の源と静的ポテンシャル

外部源が十分ゆっくり入って定常になり、自由な過渡振動は選ばない境界条件を考える。このとき $p^0=0$ の遅延核は $1/(|\vb p|^2+m^2)$。静的平衡の方程式 $(-\nabla^2+m^2)\phi=J$ の逆核に一致する。三次元で点源 $J=q\delta^3(\vb x)$ の場合、$r>0$ の回転対称な解を $\phi(r)=u(r)/r$ と書くと $u''=m^2u$ になる。無限遠で有界な解を選び、$m>0$ なら $u=Ae^{-mr}$。原点の小球に方程式を積分すれば、質量項は $r^2$ の程度で消え、$-\int_{\partial B_r}\partial_r\phi\,\dd S\to4\pi A=q$。従って静的な核は $e^{-mr}/(4\pi r)$ である。有限の広がりを持つ源ならこの分布と畳み込む。

:::example[二倍の距離での質量による抑制]
$m=1\,\mathrm{GeV}$、同じ点源に対し $r_1=1\,\mathrm{GeV}^{-1}$ と $r_2=2\,\mathrm{GeV}^{-1}$ の振幅比は $\phi(r_2)/\phi(r_1)=e^{-1}/2$。質量ゼロでの比 $1/2$ に比べてさらに $e^{-1}$ が掛かる。等時刻の真空相関は四次元の Euclid 逆核であり、ここで求めた三次元の静的応答とは違う。両者で指数の減衰長は $1/m$ だが、距離への冪や次元は異なる。
:::

### 負の虚時間へ続けて正周波数を選ぶ

$t>0$ の時間順序核は正周波数の積分である。$t=-i\tau$、$\tau>0$ へ続ければ $D_E(\tau,\vb p)=e^{-E\tau}/(2E)$。$t<0$ の枝では負周波数が選ばれるため、Euclid の負時間では $e^{-E|\tau|}/(2E)$ になる。この核は $(-\partial_\tau^2+E^2)D_E=\delta(\tau)$ を満たす。実際、正時間の導関数は $-1/2$、負時間は $+1/2$ であり、$-\partial_\tau^2$ が正のデルタを作る。Fourier 変換では $1/(\omega_E^2+E^2)$ となる。境界値としての Minkowski 核の $i0$ は、どの指数がどの虚時間の向きで減衰するかと対応する。

正周波数の二点関数自体を続けるときは $D^+(t,\vb p)=e^{-iEt}/(2E)$ を $\Im t<0$ へ移す。この関数は平坦時空で $E>0$ のためその半面で減衰する。例えば二種類の正のエネルギー $E_1<E_2$ の成分があるなら $D_E(\tau)=c_1e^{-E_1\tau}+c_2e^{-E_2\tau}$、$c_i>0$。$\tau$ が大きいと低い成分だけが残り、$-\partial_\tau\log D_E$ は $E_1$ に近づく。逆に時間順序を取っていない負周波数の積分を同じ向きに回すと指数が成長する。

:::example[二つの質量を混ぜた Euclid 相関から低い質量を読む]
零空間運動量で $D_E(\tau)=e^{-2\tau}+3e^{-4\tau}$ を単位 $\mathrm{GeV}$ に基づく数値で考える。実効エネルギーは $(2+12e^{-2\tau})/(1+3e^{-2\tau})$。$\tau=(\log3)/2$ では二つの重みが等しく、値は $3\,\mathrm{GeV}$。$\tau\to\infty$ では $2\,\mathrm{GeV}$ に近づく。時刻の違う二点関数の比から質量を推定する際、有限の時間では高い状態が混ざることをこのモデルで確認できる。
:::

:::example[四次元の遅延核を時刻固定で全空間に積分する]
質量ゼロで $t>0$ の $\delta(t-r)/(4\pi r)$ を球座標で積分すると $\int_0^\infty r\,\delta(t-r)\dd r=t$。一方、空間運動量を零にしたときの一モードの遅延核は $\lim_{E\to0}\theta(t)\sin(Et)/E=\theta(t)t$。両方は同じであり、空間の積分が零運動量 Fourier 成分を選ぶことと整合する。例えば空間全体に一様な瞬時源を置くと、その後の場は一様に時間に比例する。これは各点が局所の源から受けた応答であり、離れた点の瞬時の相互作用を仮定していない。周期箱の量子論で零モードの真空が定まらないことと、古典の遅延初期値が定まることも別である。
:::

## 演習問題

::::exercise{#exr-qft-scalar-quantization-1 level=basic}
一モードで $a^\dagger\ket3$ のノルムを求めよ。

:::solution
二乗ノルムは $\bra3aa^\dagger\ket3=3+1=4$ なのでノルムは $2$ である。
:::
::::

::::exercise{#exr-qft-scalar-quantization-2 level=basic}
$m=5\,\mathrm{GeV}$、$|\vb p|=12\,\mathrm{GeV}$ の一粒子のエネルギーを求めよ。

:::solution
$E_{\vb p}=\sqrt{5^2+12^2}\,\mathrm{GeV}=13\,\mathrm{GeV}$ である。
:::
::::

::::exercise{#exr-qft-scalar-quantization-3 level=basic}
$a$ 粒子三個と $b$ 粒子一個の正規順序電荷を求めよ。

:::solution
$Q=N_a-N_b$ なので $Q=3-1=2$ である。
:::
::::

::::exercise{#exr-qft-scalar-quantization-4 level=basic}
質量ゼロの零点密度で $\Lambda$ を倍にすると何倍になるか。

:::solution
球対称積分は $\int_0^\Lambda p^3\dd p$ に比例するので $2^4=16$ 倍である。
:::
::::

::::exercise{#exr-qft-scalar-quantization-5 level=advanced}
相異なる直交モード $f,g$ に一粒子ずつ置く。$a(f)a(g)$ の作用を求めよ。

:::solution
両モード間の交換子は零で、各モードの交換子は一である。よって $a(f)a(g)a^\dagger(f)a^\dagger(g)\ket0=a(f)a^\dagger(f)\ket0=\ket0$ となる。
:::
::::

::::exercise{#exr-qft-scalar-quantization-6 level=advanced}
$[Q,b^\dagger]=-b^\dagger$ と $[H,b^\dagger]=E b^\dagger$ を番号演算子で示せ。

:::solution
$[b^\dagger b,b^\dagger]=b^\dagger[b,b^\dagger]=b^\dagger$ である。この番号は $Q$ に負、$H$ に正の $E$ の係数で現れる。
:::
::::

::::exercise{#exr-qft-scalar-quantization-7 level=advanced}
$C$ の時間微分を等時刻で求めよ。

:::solution
時刻微分で各指数項に $-iE_{\vb p}$ または $+iE_{\vb p}$ が掛かる。二項を運動量反転して足すと $\partial_t C(0,\vb x)=-i\int\dd^3p/(2\pi)^3\,e^{i\vb p\cdot\vb x}=-i\delta^3(\vb x)$ となる。
:::
::::

::::exercise{#exr-qft-scalar-quantization-8 level=advanced}
時刻零に一モードを瞬間的に駆動する。源 $J_p(t)=A\delta(t)$ に対する遅延解を求めよ。

:::solution
零より前で場は零とする。方程式の微小時刻区間での積分により $\dot\phi_p(0^+)=A$ であり、場自体は連続なので $\phi_p(0^+)=0$。従って $\phi_p=A\theta(t)\sin(Et)/E$ である。
:::
::::

::::exercise{#exr-qft-scalar-quantization-9 level=advanced}
周波数 $\omega>0$ の一つの自由モードの古典源 $J(t)=A[\delta(t)+\delta(t-T)]$ を考える。真空から出発した最終状態の平均粒子数を求めよ。二つのパルスが打ち消し合う最小の $T>0$ はいくらか。

:::solution
$H_I=-Jq_I$、$q_I(t)=(a e^{-i\omega t}+a^\dagger e^{i\omega t})/\sqrt{2\omega}$ と書く。二つのパルスの指数の交換子は数なので、位相を除いて一つの変位にまとめられる。生成演算子の係数は $\alpha=i A(1+e^{i\omega T})/\sqrt{2\omega}$。従って $\langle N\rangle=|\alpha|^2=A^2(1+\cos\omega T)/\omega$ である。$A\ne0$ なら最初の零は $T=\pi/\omega$。反対に $T=2\pi/\omega$ なら平均は $2A^2/\omega$ で、一パルスの平均 $A^2/(2\omega)$ の四倍である。独立な確率として足すとこの値は得られない。各パルスが同じモードに残す複素振幅を先に足す必要がある。
:::
::::

::::exercise{#exr-qft-scalar-quantization-10 level=advanced}
箱モード一つの周波数を $\omega$ から $4\omega$ へ切り替える。切替前の真空を切替後の粒子数基底で測るとき、零粒子と二粒子の確率、および平均粒子数を求めよ。

:::solution
係数は $u=(2+1/2)/2=5/4$、$v=(2-1/2)/2=3/4$。古い消滅演算子は新しい $b$ で $a=u^*b-vb^\dagger$ と表せる。従って初期真空が $b$ 基底で満たす係数の漸化式は、前に求めた式で $v\mapsto-v$ としたもの。絶対値は同じなので零粒子は $P_0=1/|u|=4/5$。二粒子は
$$
P_2=P_0\,\frac{2!}{2^2(1!)^2}\left|\frac vu\right|^2
=\frac45\cdot\frac12\cdot\frac9{25}=\frac{18}{125}.
$$
従って二つの確率を足しても $118/125$ で、残りは四以上の偶数個にある。平均は $b=ua+va^\dagger$ を直接使う。$a\ket{0_a}=0$ だから $b\ket{0_a}=v a^\dagger\ket{0_a}$、そのノルム二乗が $\langle N_b\rangle=|v|^2=9/16$。零粒子から変わる確率 $1/5$ とは異なり、一度に二個以上が見つかる重みも平均に入る。
:::
::::

## 確認テスト

::::quiz
空間的分離で自由実スカラー場に成り立つものを選べ。

- [x] 場の交換子が消える
- [ ] 各順序の二点期待値が消える
- [ ] 場の反交換子が恒等的に消える

:::explanation
等時刻に移せば交換子の二項は運動量反転で同じ積分となる。個々の積分にその相殺は起こらない。
:::
::::

::::quiz{answer="3" tol="1e-6"}
一振動子で $\bra2aa^\dagger\ket2$ を求めよ。

:::explanation
$aa^\dagger=N+1$ なので固有値は $2+1=3$ である。
:::
::::

::::quiz
遅延 Green 関数の両極の位置を選べ。

- [ ] 両方とも上半面
- [x] 両方とも下半面
- [ ] 正のエネルギーだけ下半面

:::explanation
負の時刻には上半面を閉じるので、遅延条件には上半面に極を置かない。両極を下半面に置くと正の時刻に応答を得る。
:::
::::

::::quiz{answer="-2" tol="1e-6"}
$a$ 粒子一個、$b$ 粒子三個の電荷を答えよ。

:::explanation
$N_a-N_b=1-3=-2$ である。
:::
::::

::::quiz{answer="0.5625" tol="1e-6"}
一つのモードで $b=(5a+3a^\dagger)/4$ とする。$a$ の真空における $b^\dagger b$ の期待値を数値で答えよ。

:::explanation
$a$ が真空を消すので $b\ket{0_a}=(3/4)a^\dagger\ket{0_a}$。一粒子状態 $a^\dagger\ket{0_a}$ のノルム二乗は交換関係により一である。従って $\bra{0_a}b^\dagger b\ket{0_a}=9/16=0.5625$。
:::
::::

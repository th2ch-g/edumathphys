---
title: "古典場の理論"
description: "作用から場の方程式と保存則を導き、相互作用の質量次元を求める。"
notation:
  - symbol: '\mathcal L,\mathcal H'
    meaning: Lagrangian 密度と Hamiltonian 密度
  - symbol: '\Pi^\mu_a,\pi_a'
    meaning: 共役密度とその時間成分
  - symbol: 'T^{\mu\nu}'
    meaning: エネルギー運動量テンソル
---

## この章の目標

- 作用を変分し、境界条件と運動方程式を求める。
- 場の共役運動量を定義し、Hamilton 形式に移る。
- 並進と位相回転に対応する保存流を構成する。
- 自然単位系で質量次元を計算する。

## 場の Lagrangian 形式と作用原理

### 配置を時空上の関数で表す

弦を記述するとき、各時刻に一点の座標を測るだけでは弦の形を定められない。弦の各点に変位を割り当てる必要がある。三次元空間でも、位置ごとに値をもつ関数 $\phi^a(t,\vb x)$ を使う。$a$ は成分の番号である。弦の連続極限については [連続体と場の解析力学](/physics/analytical-mechanics/field-mechanics/) を参照する。

:::definition[局所作用]{#def-qft-local-action}
同じ点での $\phi$ とその微分を組み合わせた :term[Lagrangian 密度]{en="Lagrangian density" yomi="らぐらんじあんみつど"} を $\mathcal L$ とし、$S[\phi]=\int_\Omega\dd^4x\,\mathcal L$ を :term[作用]{en="action" yomi="さよう"} という。以下では $\mathcal L$ は一階微分までを含む。
:::

$\hbar=c=1$、$\eta=\diag(1,-1,-1,-1)$ とする。:ref[def-sr-fourvector] の計量とは全体の符号が逆である。添字を $\eta$ で上下し、$\Box=\partial_\mu\partial^\mu=\partial_t^2-\nabla^2$ と書く。実スカラー場 $\phi$ の例は
$$
\mathcal L=\tfrac12\partial_\mu\phi\partial^\mu\phi-\tfrac12m^2\phi^2+J\phi
$$
である。源 $J$ もスカラーとして変換すれば、Lorentz 変換で $\mathcal L$ と $\dd^4x$ が保たれる。有限領域なら領域自体と境界条件を一緒に写す。例えば $m=2\,\mathrm{GeV}$、一様で静止した瞬間の $\phi=1\,\mathrm{GeV}$、$J=0$ なら $\mathcal L=-2\,\mathrm{GeV}^4$ である。

### 変分に対する境界条件

経路 $\phi+s h$ に沿って $s$ で微分する。時間端点で $h=0$ とし、空間方向は周期境界、または無限遠で境界積分が消える条件を採る。$\mathcal L$ に全微分を足すと作用の差は境界積分である。従って変分でその積分も固定していれば内部の運動方程式が一致する。自由境界では逆に、その境界で残る変分の係数を零にする必要がある。

### 同じ内部方程式から異なる端の条件が出る例

一次元の半直線を扱い、$\mathcal L=(\dot\phi^2-\phi_x^2-m^2\phi^2)/2$ とする。$x=0$ の値を固定すれば許される変分は $h(0)=0$ である。端に作用を足さず $\phi(0)$ を自由にすると、$\int_0^\infty-\phi_x h_x\,\dd x$ の端から $+\phi_x(0)h(0)$ が出る。従って端の条件は $\phi_x(0)=0$ となる。境界作用 $-\int\dd t\,b\phi(0,t)^2/2$ を加えれば $\phi_x(0)-b\phi(0)=0$ である。例えば端へ向かう波と反射波を $e^{-i\omega t}(e^{-ikx}+Re^{ikx})$ と取ると $ik(R-1)=b(1+R)$。従って $R=(ik+b)/(ik-b)$ であり、実数 $b$ に対し $|R|=1$ である。$b=0$ では同位相、$|b|\to\infty$ では逆位相になる。

### 連成振動子から連続体へ

長さ $L$ の円周を格子間隔 $a=L/N$ で分ける。変位を $q_j$、線密度を $\rho$、張力を $\tau$ とし、$q_N=q_0$ とする。各質点の質量は $\rho a$ である。隣接区間の勾配を $(q_{j+1}-q_j)/a$ と近似すると
$$
L_a=a\sum_{j=0}^{N-1}\left[\frac{\rho}{2}\dot q_j^2
-\frac{\tau}{2}\left(\frac{q_{j+1}-q_j}{a}\right)^2-\frac{\rho m^2}{2}q_j^2\right].
$$
$q_j$ を含む二つのばねをそれぞれ微分すると、$\rho\ddot q_j=\tau(q_{j+1}-2q_j+q_{j-1})/a^2-\rho m^2q_j$ である。滑らかな $q(t,x)$ を格子に制限したとき、二階差分の Taylor 展開は $q_{xx}+a^2q_{xxxx}/12+\order{a^4}$ となる。従って固定波長で $a\to0$ とすると $q_{tt}-(\tau/\rho)q_{xx}+m^2q=0$ を得る。これに対し波長を格子間隔に比例させると、この展開には小さいパラメータがない。

:::example[連続近似が使える波数]
$q_j=Qe^{i(kja-\omega t)}$ を差分式に代入すると $\omega^2=m^2+4(\tau/\rho)\sin^2(ka/2)/a^2$ である。$m=0$、$\tau/\rho=1$、$ka=\pi/3$ では格子上の $\omega=1/a$ に対し、連続近似は $\pi/(3a)$ である。格子上の群速度は $\cos(ka/2)=\sqrt3/2$ である。従って周波数の相対誤差と伝播速度の相対誤差は別々に見積もる必要がある。
:::

:::example[停留解に沿う二次変分]
周期箱の一様モードで $\phi(t)=s h(t)$ とし、両端で $h=0$ とする。$J=0$ なら $s=0$ は停留点で、二次変分は $V\int_0^T(\dot h^2-m^2h^2)\dd t$ である。$h(t)=\sin(n\pi t/T)$ のとき $\int\sin^2=\int\cos^2=T/2$ から $VT[(n\pi/T)^2-m^2]/2$ となる。$mT=2\pi$ では $n=1$ と $n=3$ で符号が逆である。従って作用原理の「停留」には鞍点も含まれる。
:::

### 二階微分を含む作用で固定する端のデータ

一階微分の密度では、部分積分は変分の微分を一度だけ除く。二階の微分を含むと、変分そのものに加えてその一階微分も端に残る。例えば時間だけの作用
$$
S[q]=\int_0^T\dd t\left(\frac{\dot q^2}{2}
-\frac{\omega^2q^2}{2}-\frac{a}{2}\ddot q^2\right)
$$
では
$$
\delta S=\int_0^T\dd t\,(-\omega^2 q-\ddot q-aq^{(4)})h
+[(\dot q+a q^{(3)})h-a\ddot q\,\dot h]_0^T .
$$
従って $h=\dot h=0$ で固定するか、$a\ne0$ で $h=0$ だけなら端で $\ddot q=0$ を自然条件として課す。一般に $\mathcal L(\phi,\partial\phi,\partial^2\phi)$ からは、コンパクト台の変分を用いて
$$
\frac{\partial\mathcal L}{\partial\phi}
-\partial_\mu\frac{\partial\mathcal L}{\partial(\partial_\mu\phi)}
+\partial_\mu\partial_\nu\frac{\partial\mathcal L}{\partial(\partial_\mu\partial_\nu\phi)}=0
$$
を得る。二回の部分積分で最後は正符号である。高階微分を基本の運動方程式として解くときと、低エネルギーで先頭方程式への小さい補正として扱うときでは、許す初期データが違う。

:::example[微小係数が増やす振動枝]
上の $a>0$ で $q=e^{-i\Omega t}$ を代入すると $a\Omega^4-\Omega^2+\omega^2=0$。$4a\omega^2<1$ なら $\Omega_\pm^2=[1\pm\sqrt{1-4a\omega^2}]/(2a)$。$a\to0$ で小さい枝は $\omega^2$、大きい枝は $1/a$ の程度となる。例えば $a\omega^2=3/16$ では比 $\Omega_+^2/\Omega_-^2=3$。高周波数枝は元の二階運動方程式にはなく、誘導された高次項を非摂動的に解くときだけ現れる。
:::

### 異なる境界値を結ぶ作用の変分と古典解の端点

作用原理では内部の経路を比べるときに端点を固定する。一方、既に求めた古典解を端点の異なる古典解へ移すと、その作用の変化は境界だけで決まる。ここでは時間区間を $[t_i,t_f]$ と固定し、空間は周期箱とする。両端の場の配置を $\phi_i,\phi_f$ とし、その近くで古典解が一意で境界値に滑らかに依存すると仮定する。古典解に作用を代入した $S_{\rm cl}[\phi_f,\phi_i]$ を変分すると、内部の項は $E_a=0$ で消え、
$$
\delta S_{\rm cl}=\int\dd^3x\,\bigl[\pi_f\,\delta\phi_f-\pi_i\,\delta\phi_i\bigr].
$$
右辺の符号の違いは境界の法線の向きによる。従って、両端の場を引数とする $S_{\rm cl}$ から最終と初期の共役運動量が求まる。この変分公式は境界値問題の解が分岐する付近では、選んだ枝ごとに適用する。

:::example[一モードの端点作用を求める]
自由な Fourier モードを実座標 $q$ に取り、$t_i=0,t_f=T$、$\omega>0$ とする。$\sin\omega T\ne0$ なら与えられた $(q_i,q_f)$ を結ぶ解は
$$
q(t)=q_i\frac{\sin[\omega(T-t)]}{\sin\omega T}
     +q_f\frac{\sin(\omega t)}{\sin\omega T}.
$$
二階微分すれば $-\omega^2 q$、両端の代入で境界値を満たす。端点作用では $\dot q^2-\omega^2q^2=\dd(q\dot q)/\dd t$ だから、途中の積分をせず
$$
S_{\rm cl}(q_f,q_i)=
\frac{\omega}{2\sin\omega T}
\left[(q_f^2+q_i^2)\cos\omega T-2q_iq_f\right]
$$
を得る。例えば $\omega T=\pi/2$、$q_i=2,q_f=3$ なら $S_{\rm cl}=-6\omega$。終点での微分は $\partial_{q_f}S_{\rm cl}=-2\omega$、始点での負の微分は $-\partial_{q_i}S_{\rm cl}=3\omega$ になる。$q(t)$ を直接微分しても $\dot q(T)=-2\omega$、$\dot q(0)=3\omega$ である。数値の座標は規格化したモードの単位である。運動項を $1/2$ にしたモードであれば $q$ の質量次元は $-1/2$ で、$\omega q^2$ の作用は無次元である。
:::

:::example[端点解が一意でない時刻]
同じ振動子に対し $\omega T=\pi$ とする。初期速度を $p_i$ とすれば $q(T)=q_i\cos\pi+(p_i/\omega)\sin\pi=-q_i$ となる。従って $q_f\ne-q_i$ なら古典解が存在せず、$q_f=-q_i$ なら任意の $p_i$ が同じ位置に到達する。この場合、初期値問題の発展は依然一意である。二点の値から初期速度を逆算する段階だけが特異になっている。作用の公式で $\sin\omega T=0$ を分母から消すことは、存在する解の範囲と余分な自由度を確認してから行う。
:::

### 薄い界面の作用と波の接続条件

媒質を一次元の二領域に分け、界面 $x=0$ で場 $\phi$ の値を連続とする。曲がりのある関数 $\phi$ を許すが、各半直線では滑らかと仮定する。運動密度を $\mathcal L_\pm=(\rho_\pm\dot\phi^2-\tau_\pm\phi_x^2)/2$、$\rho_\pm,\tau_\pm>0$ とし、界面の作用に $-\int\kappa\phi(0,t)^2\dd t/2$ を入れる。左側の部分積分の上端は $-\tau_-\partial_x\phi(0^-)\delta\phi$、右側の下端は $+\tau_+\partial_x\phi(0^+)\delta\phi$。界面での変分値も連続なので一つの係数にまとめられ、停留条件は
$$
\tau_+\partial_x\phi(0^+)-\tau_-\partial_x\phi(0^-)=\kappa\phi(0).
$$
$\kappa=0$ では張力を掛けた勾配が連続する。$\rho$ と $\tau$ が場所で異なる場合、単に $\phi_x$ の連続性を課すと界面の運動量の流れを誤って扱う。

:::example[密度が異なる弦での反射と透過]
角周波数 $\omega$ を固定し、左からの波を $e^{-i\omega t}(e^{ik_-x}+Re^{-ik_-x})$、右側を $e^{-i\omega t}\mathcal T e^{ik_+x}$ とする。$k_\pm=\omega\sqrt{\rho_\pm/\tau_\pm}$。$\kappa=0$ のとき場の連続性から $\mathcal T=1+R$、勾配条件から $Z_+\mathcal T=Z_-(1-R)$、$Z_\pm=\sqrt{\rho_\pm\tau_\pm}$。二式を解くと
$$
R=\frac{Z_--Z_+}{Z_-+Z_+},\qquad
\mathcal T=\frac{2Z_-}{Z_-+Z_+}.
$$
エネルギー密度 $u=(\rho\dot\phi^2+\tau\phi_x^2)/2$ の方程式は $\partial_tu=\partial_x(\tau\dot\phi\,\phi_x)$。従って右向き波の時間平均の流束は $\omega^2 Z|\mathcal A|^2/2$ である。例えば $Z_+=2Z_-$ なら $R=-1/3,\mathcal T=2/3$。流束の比は反射が $1/9$、透過が $2(2/3)^2=8/9$ で、和が一となる。境界で連続なのは振幅の和であり、透過率には媒質の $Z$ も入る。
:::

## Euler–Lagrange 方程式と場の Hamilton 形式

:::theorem[場の Euler–Lagrange 方程式]{#thm-qft-field-euler}
滑らかな $\mathcal L$ と二階微分可能な場に対し、固定境界での停留条件は
$$
E_a\defeq\pdv{\mathcal L}{\phi^a}-\partial_\mu\Pi^\mu_a=0,\qquad
\Pi^\mu_a\defeq\pdv{\mathcal L}{(\partial_\mu\phi^a)}
\label{eq-qft-euler}
$$
と同値である。
:::

:::proof
連鎖律と積の微分により $\delta\mathcal L=E_ah^a+\partial_\mu(\Pi^\mu_a h^a)$ である。境界項は消える。ある点で $E_a>0$ なら、その近傍に台をもつ非負の $h^a$ を選び $\delta S>0$ にできる。負の場合も符号を変えればよい。従って停留なら $E_a=0$ であり、逆も積分式から従う。
:::

### Legendre 変換

$\pi_a=\Pi^0_a$ を :term[共役運動量密度]{en="canonical momentum density" yomi="きょうやくうんどうりょうみつど"} とする。速度についての Hessian が可逆なら、逆関数定理で局所的に $\dot\phi$ を $(\phi,\pi,\partial_i\phi)$ で表せる。その範囲で $H=\int\dd^3x\,(\pi_a\dot\phi^a-\mathcal L)$ とする。変分の $\delta\dot\phi^a$ の係数が相殺するので
$$
\delta H=\int\dd^3x\,\left[\dot\phi^a\delta\pi_a+
\left(-\pdv{\mathcal L}{\phi^a}+\partial_i\Pi^i_a\right)\delta\phi^a\right].
$$
汎関数微分を変分の係数として定義すると、:eqref[eq-qft-euler] から $\dot\phi^a=\delta H/\delta\pi_a$、$\dot\pi_a=-\delta H/\delta\phi^a$ を得る。従って :term[Poisson 括弧]{en="Poisson bracket" yomi="ぽあそんかっこ"}
$$
\{F,G\}=\int\dd^3x\,\left(\fdv{F}{\phi^a}\fdv{G}{\pi_a}-\fdv{F}{\pi_a}\fdv{G}{\phi^a}\right)
\label{eq-qft-field-poisson}
$$
によって $\dot F=\partial_t F+\{F,H\}$ となる。点で評価した場の括弧は試験関数で積分すると $\{\phi^a(\vb x),\pi_b(\vb y)\}=\delta^a_b\delta^3(\vb x-\vb y)$ である。Maxwell 場は $\dot A_0$ が密度に現れず Hessian が特異になる。その場合は [ベクトル場と電磁場の量子化](/physics/quantum-field-theory/vector-fields/) で扱う。

:::example[実スカラー場の発展]
$\pi=\dot\phi$、$\mathcal H=\pi^2/2+|\nabla\phi|^2/2+m^2\phi^2/2-J\phi$ である。従って $\dot\pi=\nabla^2\phi-m^2\phi+J$ となり、$(\Box+m^2)\phi=J$ に戻る。自由平面波の $m=3\,\mathrm{GeV}$、$|\vb{k}|=4\,\mathrm{GeV}$ なら $\omega=\sqrt{k^2+m^2}=5\,\mathrm{GeV}$ である。
:::

### 初期波束の時間発展

$\widehat\phi(\vb k,t)$ を空間 Fourier 変換とし、初期値 $\widehat f,\widehat g$ を時刻零の場と速度とする。$\omega_k^2=|\vb k|^2+m^2$ に対し一モードの方程式は $\ddot{\widehat\phi}+\omega_k^2\widehat\phi=\widehat J$。従って
$$
\widehat\phi(\vb k,t)=\widehat f(\vb k)\cos\omega_k t
 +\widehat g(\vb k)\frac{\sin\omega_k t}{\omega_k}
 +\int_0^t\dd s\,\frac{\sin[\omega_k(t-s)]}{\omega_k}\widehat J(\vb k,s).
$$
微分すると最初の境界項は $\sin0=0$、二回目には $\cos0=1$ により $\widehat J(\vb k,t)$ が出る。これで方程式と二つの初期条件を満たす。$m=0$ かつ $\vb k=0$ では $\sin(\omega_k t)/\omega_k$ の極限は $t$ である。例えば $\widehat f=0$、$\widehat g=0$、一定の $\widehat J=J_0$ なら $\omega_k>0$ で $\widehat\phi=J_0(1-\cos\omega_k t)/\omega_k^2$。$\omega_k\to0$ で $J_0t^2/2$ に移る。

### 正準変換と境界項

$\mathcal L'=\mathcal L+\partial_t f(\phi,t)$ とする。$f$ が場の空間微分を含まなければ $\pi'_a=\pi_a+\partial f/\partial\phi^a$ である。時刻を固定した正準一形式の差は $\delta\int f\,\dd^3x$ だから、二度目の外微分は零である。従って $\int\delta\pi'_a\wedge\delta\phi^a=\int\delta\pi_a\wedge\delta\phi^a$ となる。また $H'=H(\phi,\pi'-f_\phi)-\int f_t\,\dd^3x$ である。ここで $f_t$ は $\phi$ 固定での微分である。この項は時間に依存する正準変換で必要になる。

:::example[境界項を加えた振動子]
$\mathcal L=(\dot\phi^2-\omega^2\phi^2)/2$、$f=b(t)\phi^2/2$ なら $\pi'=\dot\phi+b\phi$、$\mathcal H'=(\pi'-b\phi)^2/2+(\omega^2-\dot b)\phi^2/2$ である。Hamilton 方程式は $\dot\phi=\pi'-b\phi$、$\dot\pi'=b(\pi'-b\phi)-(\omega^2-\dot b)\phi$ を与える。最初の式を微分して後の式を代入すると、$b\dot\phi$ と $\dot b\phi$ がそれぞれ相殺し、$\ddot\phi=-\omega^2\phi$ に戻る。
:::

:::example[連成した二成分の正規モード]
$\Phi=(\phi_1,\phi_2)^\mathsf T$ の運動項を標準形とし、質量二乗行列を $M^2=\begin{pmatrix}5&4\\4&5\end{pmatrix}\,\mathrm{GeV}^2$ とする。$\phi_\pm=(\phi_1\pm\phi_2)/\sqrt2$ はそれぞれ固有値 $9,1$ の方向である。同じ直交変換を $\pi$ に掛けると $\pi^\mathsf T\delta\Phi=\pi_+\delta\phi_++\pi_-\delta\phi_-$ である。従って波数 $4\,\mathrm{GeV}$ での二つの正準振動子の角周波数は $5\,\mathrm{GeV}$ と $\sqrt{17}\,\mathrm{GeV}$ になる。
:::

### 非正準の運動項を正準化する

二次の多成分場で $\mathcal L=\frac12\partial_\mu\Phi^\mathsf T K\partial^\mu\Phi-\frac12\Phi^\mathsf T M^2\Phi$、$K=K^\mathsf T>0$、$M^2=(M^2)^\mathsf T$ とする。$\Psi=K^{1/2}\Phi$ を使うと運動項は単位行列、質量は $K^{-1/2}M^2K^{-1/2}$。後者は実対称なので直交行列で対角化できる。元の変数から見ると固有値の式は $M^2v=\lambda K v$ であり、固有ベクトルは $v_i^\mathsf T K v_j=\delta_{ij}$ に選ぶ。実際、異なる固有値をもつ $i,j$ について $v_i^\mathsf T M^2 v_j$ を二通りに計算すれば $(\lambda_i-\lambda_j)v_i^\mathsf T K v_j=0$ が従う。共役運動量は $\pi=K\dot\Phi$ であり、$\pi'=K^{-1/2}\pi$ とすれば $\pi'^\mathsf T\delta\Psi=\pi^\mathsf T\delta\Phi$ である。

:::example[運動項と質量項が可換でない二成分場]
$K=\diag(4,1)$、$M^2=\begin{pmatrix}20&8\\8&5\end{pmatrix}\,\mathrm{GeV}^2$ とする。$K^{-1/2}=\diag(1/2,1)$ を両側から掛けると $\begin{pmatrix}5&4\\4&5\end{pmatrix}\,\mathrm{GeV}^2$ になる。従って質量は $3,1\,\mathrm{GeV}$。元の二成分の比は $v_\pm=(1/2,\pm1)^\mathsf T/\sqrt2$ である。Euclid の内積だけを取ると両者の内積は $-3/8$ だが、運動項 $K$ を挟んだ積では零になる。
:::

### エネルギー法による領域依存性

自由な Klein–Gordon 方程式の解の差を $u$ とする。一時刻で半径 $R$ の球内のデータが一致したなら、その内部に未来の円錐が入っている点では解も一致する。これをエネルギー収支から示す。時刻 $t$ で半径 $R-t$ の球を $B_t$ とし、
$$
\mathcal E(t)=\int_{B_t}\dd^3x\,
    \tfrac12[(\partial_t u)^2+|\nabla u|^2+m^2u^2],\qquad 0\le t<R
$$
とする。積の微分と方程式から $\partial_t$ を密度に作用させると $\divg(\partial_tu\,\nabla u)$ となる。積分領域の境界は速さ $1$ で内へ動くため、密度の流出も差し引く。単位外法線を $\vb n$、接線の勾配を $\nabla_Tu$ とすると
$$
\begin{aligned}
\dot{\mathcal E}
&=\int_{\partial B_t}\dd S\left[
\partial_t u\,\partial_nu-\tfrac12((\partial_t u)^2+|\nabla u|^2+m^2u^2)\right]\\
&=-\frac12\int_{\partial B_t}\dd S\,
\bigl[(\partial_tu-\partial_nu)^2+|\nabla_Tu|^2+m^2u^2\bigr]\le0 .
\end{aligned}
$$
$\mathcal E(0)=0$ なら全時刻で零である。$m>0$ では密度の各項の非負性から $u=0$。$m=0$ でも時間と空間の微分が零であり、初期の値が零なのでその定数も零である。こうして源や初期値の違いが、光より速く内部へ入らないことが分かる。$m^2\ge0$ と滑らかさを仮定した証明であり、変係数の理論に移すときは収支をその密度で作り直す。

:::example[球の外から変えた初期値の影響]
時刻零で中心の周り $R=4\,\mathrm{GeV}^{-1}$ の球に二つの同じ初期値を与え、外側では異なるとする。$t=1.5\,\mathrm{GeV}^{-1}$ の時点で、中心から $2.5\,\mathrm{GeV}^{-1}$ 未満の点では二解の差が零である。逆に、時刻 $3\,\mathrm{GeV}^{-1}$ に中心から $2\,\mathrm{GeV}^{-1}$ の点はこの証明の領域に含まれない。この点で差が実際に非零になることまでは、領域依存性だけからは結論できない。
:::

### 陽に時間に依存する駆動と初期値

有限箱の一モードで $H=\pi^2/2+\omega^2q^2/2-J(t)q$ とする。初期時刻を $0$、$q(0)=\pi(0)=0$ とし、$J(t)=J_0\cos\nu t$ を $t\ge0$ で作用させる。$\nu\ne\omega$ のとき特解を $A\cos\nu t$ とすれば $A=J_0/(\omega^2-\nu^2)$。初期変位を消すために $-A\cos\omega t$ を加え、
$$
q(t)=\frac{J_0}{\omega^2-\nu^2}(\cos\nu t-\cos\omega t)
$$
となる。$\nu\to\omega$ を極限で取ると、分子と分母を $\nu$ で微分して $q(t)=J_0t\sin\omega t/(2\omega)$ を得る。従って共鳴では、自由振動を加えて初期値を合わせた後も時間に比例する包絡線が残る。時刻 $t$ までの駆動源とのエネルギー移動は $\dd(\pi^2/2+\omega^2q^2/2)/\dd t=J\dot q$。源と場の相互作用まで入れた Hamiltonian の微分は $-\dot Jq$ であり、両式の差は $\dd(Jq)/\dd t$ で一致する。

:::example[共鳴に近い二つの角周波数]
$\omega=3\,\mathrm{GeV}$、$\nu=2\,\mathrm{GeV}$ とし、$t=\pi\,\mathrm{GeV}^{-1}$ を取る。$\cos\nu t=1$、$\cos\omega t=-1$ だから変位は $2J_0/(5\,\mathrm{GeV}^2)$。同時刻で両正弦が零となり速度は零である。もし同じ源の $\nu$ を $\omega$ に合わせるとこの時刻で変位は零だが、速度は $-J_0t/2$ になる。従って変位だけを比較すると、共鳴で蓄えられたエネルギーを見落とす。
:::

### 特異な Legendre 変換を補助変数で扱う

速度をもたない実変数 $u$ も含め、単一モードの作用を
$$
L(q,\dot q,u)=\tfrac12\dot q^2-\tfrac12\omega^2q^2+\tfrac12M^2u^2-g u q,
\qquad M^2>0
$$
とする。速度の Hessian は $\diag(1,0)$ なので $p_q=\dot q$ だが $p_u=0$ である。$u$ の運動方程式は $M^2u=gq$ であり、初期値として任意に選べない。これを $L$ に入れると有効な一変数の密度は $\dot q^2/2-(\omega^2+g^2/M^2)q^2/2$。$u$ の部分を平方完成しても同じ式を得る。代入後の変分では $\delta u$ が生じるが、元の $u$ の方程式を満たすため係数が零なので、元の $q$ の方程式と同値である。

Hamilton 形式で $u$ も座標として保つなら $H_{\rm tot}=p_q^2/2+\omega^2q^2/2-M^2u^2/2+guq+\lambda p_u$ と取る。$\lambda$ は $p_u=0$ を課す未定乗数。Hamilton 方程式の $\dot p_u=-\partial H_{\rm tot}/\partial u=M^2u-gq$ より、一次の条件の時間保存に $M^2u-gq=0$ が要る。これを時間微分し $\dot u=\lambda,\dot q=p_q$ を使うと $\lambda=gp_q/M^2$ が決まる。従って独立な初期データは $q,p_q$ の一組である。ここでは $\{p_u,M^2u-gq\}=-M^2$ が可逆なので、二つの条件は同じゲージ軌道を表す冗長性ではない。$u$ のような非伝播変数の消去と、Maxwell のゲージ自由度の消去は、速度 Hessian が特異という共通点をもつが、その後の拘束の括弧で区別する。

:::example[補助変数を含む系の周波数]
規格化した単一モードで $\omega=3\,\mathrm{GeV}$、$M=2\,\mathrm{GeV}$、$g=8\,\mathrm{GeV}^2$ とする。$u=2q$、$p_u=0$ だから、初期値 $q(0)=q_0,p_q(0)=0$ の時間発展は $q(t)=q_0\cos(5\,\mathrm{GeV}\,t)$、$u(t)=2q(t)$ である。$H_{\rm tot}$ の $-M^2u^2/2$ だけを見ると下に無限に広がるが、同じ拘束を課した全体が与えるエネルギーは $p_q^2/2+25\,\mathrm{GeV}^2q^2/2$。従って物理的な初期データ上では安定である。
:::

## Noether の定理と保存カレント

### 定数の変換を局所化する

:ref[thm-am-field-noether] は弦を含む連続系で保存流を与える。ここでは全微分まで作用を保つ変換を扱い、後で量子論に移せる形で書く。固定座標で $\delta\phi^a=\epsilon R^a$ とし、定数 $\epsilon$ に対し $\delta\mathcal L=\epsilon\partial_\mu K^\mu$ と仮定する。

:::theorem[Noether の恒等式]{#thm-qft-noether-current}
$j^\mu=\Pi^\mu_a R^a-K^\mu$ と置けば
$$
\partial_\mu j^\mu=-E_a R^a
\label{eq-qft-noether}
$$
である。運動方程式と無流束の境界条件の下で $Q=\int\dd^3x\,j^0$ は時間に依存しない。
:::

:::proof
第一式は $\delta\mathcal L=\epsilon E_aR^a+\epsilon\partial_\mu(\Pi^\mu_aR^a)$ を仮定と比較すれば従う。第二式は有限領域 $V$ で $\dot Q=-\int_{\partial V}\vb j\cdot\dd\vb S$ と書いて流束を零にすればよい。
:::

$\epsilon(x)$ を試験関数にすると $\mathcal L$ の変分は $\partial_\mu(\epsilon K^\mu)+j^\mu\partial_\mu\epsilon$ になる。微分が $\epsilon$ に作用する箇所を残したためである。また $j^\mu$ に反対称テンソルの発散 $\partial_\nu U^{[\nu\mu]}$ を加えても保存式は変わらない。二階微分が対称であるため、電荷の差も境界項になる。

:::example[シフト対称性]
自由な質量ゼロの場では $R=1$、$K=0$ だから $j^\mu=\partial^\mu\phi$ である。周期箱の一様解 $\phi=A+Bt$ の電荷は $VB$ である。質量項を加えるとシフトに対し $\delta\mathcal L=-\epsilon m^2\phi$ なので、解の上で流の発散が $-m^2\phi$ になる。
:::

### 非可換な内部対称性を持つ実場

$N$ 成分の実場 $\Phi$ に $\mathcal L=(\partial\Phi)^2/2-V(\Phi^{\mathsf T}\Phi)$ を取る。反対称実行列 $A$ で $R=A\Phi$ とすると、$\delta(\Phi^{\mathsf T}\Phi)=2\Phi^{\mathsf T}A\Phi=0$。運動項も各 $\partial_\mu\Phi$ について同じ理由で不変である。従って流は $j_A^\mu=(\partial^\mu\Phi)^{\mathsf T}A\Phi$ であり、$Q_A=\int\pi^{\mathsf T}A\Phi\,\dd^3x$ である。汎関数微分は $\delta Q_A/\delta\pi=A\Phi$、$\delta Q_A/\delta\Phi=A^{\mathsf T}\pi$。:eqref[eq-qft-field-poisson] により $Q_A$ は場と運動量に同じ回転を生成する。また
$$
\{Q_A,Q_B\}=Q_{[A,B]}
$$
となる。括弧の第一項は $\int\pi^{\mathsf T}AB\Phi$、第二項は $\int\pi^{\mathsf T}BA\Phi$ である。例えば三成分で $A\Phi=\vb a\times\Phi$、$B\Phi=\vb b\times\Phi$ とすれば、$[A,B]\Phi=(\vb a\times\vb b)\times\Phi$ となる。こうして正準変数で計算した電荷の括弧が回転の合成則を再現する。

### 局所変換と Noether の第二定理

:ref[thm-qft-noether-current] では変換のパラメータは定数であった。これに対し任意のコンパクトな台をもつ $\alpha(x)$ に対して作用が不変なら、運動方程式を使う前から変分微分の間に恒等式がある。$\delta\phi^a=R^a\alpha+R^{a\mu}\partial_\mu\alpha$ とすると
$$
0=\delta S=\int\dd^4x\,\alpha\left[E_aR^a-\partial_\mu(E_aR^{a\mu})\right].
$$
境界値の固定と部分積分を用いた。$\alpha$ の任意性から角括弧が恒等的に零である。これを :term[Noether の第二定理]{en="Noether's second theorem" yomi="のーたーのだいにていり"} という。この恒等式は運動方程式の独立な成分数を数えるのに使う。

:::example[Maxwell 方程式と保存源]
$E^\nu=\partial_\mu F^{\mu\nu}-j^\nu$ とする。$\delta A_\nu=\partial_\nu\alpha$ なので上式の要請は $\partial_\nu E^\nu=0$ である。反対称性から $\partial_\nu\partial_\mu F^{\mu\nu}=0$ であり、従って外部源には $\partial_\nu j^\nu=0$ が必要になる。時間だけで変わる一様な電荷密度を、空間電流を零としたまま外部源に指定すると、この条件に反する。
:::

:::example[電荷から位相空間上の変換を復元する]
一様な二成分振動子で $Q=q_1p_2-q_2p_1$ とする。$\{q_i,Q\}=\partial Q/\partial p_i$ より $\dot q_1=-q_2,\dot q_2=q_1$ であり、ここで点は変換パラメータについての微分である。さらに $\dot p_1=-p_2,\dot p_2=p_1$ となる。この四式を積分すると $q,p$ が同じ角度だけ回る。$q=(2,1)$、$p=(3,-1)$ の $Q=-5$ は、$\pi/2$ 回転後の $q'=(-1,2)$、$p'=(1,3)$ でも $-5$ である。
:::

### Noether 流の保存と外からの駆動

内部回転 $\delta\phi^a=\alpha R^a$ に対して $\delta\mathcal L=\alpha\partial_\mu K^\mu+\alpha B$ とする。変分の公式を二通りに書いて引くと $\partial_\mu j^\mu=B-E_aR^a$。$B$ は明示的破れであり、場の方程式を課しても残る。例えば $O(2)$ の二実場に外部源 $J_1\phi_1+J_2\phi_2$ を足し、$\delta\phi_1=-\alpha\phi_2$、$\delta\phi_2=\alpha\phi_1$ とすると $B=J_2\phi_1-J_1\phi_2$。場の電荷 $Q=\int(\phi_1\pi_2-\phi_2\pi_1)\dd^3x$ の時間微分は、流束に加えて $\int B$ を持つ。源も同じ回転で変換する規則を与えれば、密度の内積は共変になる。ただし、固定した源のもとで場だけに計算した電荷には上の駆動が入る。

:::example[空間的に一様な駆動の初期トルク]
一様に $\Phi=(2,1)\,\mathrm{GeV}$、$J=(3,5)\,\mathrm{GeV}^3$ と置く。初期時刻での電荷密度の変化率は $5\cdot2-3\cdot1=7\,\mathrm{GeV}^4$。同じ値は $Q/V=\phi_1\dot\phi_2-\phi_2\dot\phi_1$ を時間微分し、時間発展の二式 $\ddot\phi_i=J_i-m^2\phi_i$ を代入して得る。質量による二項は逆符号で相殺する。従って初期の運動量が零でも、源が場と平行でなければ電荷が発生する。
:::

### 正準変換として Noether の電荷を導く

場に関する変換 $\delta\phi^a=R^a(\phi,t)$ を考え、ここでは速度や空間微分を含まないとする。定数パラメータで $\delta\mathcal L=\partial_t K^0(\phi,t)+\partial_i K^i$ と書けると仮定する。$\partial\delta\mathcal L/\partial\dot\phi^a$ を計算する。左辺の変分は $\mathcal L_{\phi^b}R^b+\pi_b\partial_t R^b+\Pi^i_b\partial_iR^b$。従って $\dot\phi^a$ による微分は、$\pi_a$ が場と速度の変化を受ける分 $\delta\pi_a$ と、$\pi_b R^b{}_{,\phi^a}$ の和になる。ただし $K^i$ に速度依存性がない場合、右辺から出るのは $K^0{}_{,\phi^a}$ だけである。従って
$$
\delta\pi_a=-\pi_b R^b{}_{,\phi^a}+K^0{}_{,\phi^a}.
$$
一方、電荷 $Q=\int(\pi_aR^a-K^0)\dd^3x$ の $\pi$ 微分は $R^a$、$\phi$ 微分の負は上式である。従って :eqref[eq-qft-field-poisson] で $\{\phi^a,Q\}=\delta\phi^a,\{\pi_a,Q\}=\delta\pi_a$ となる。$\partial_t R^a$ を全微分として計算したので、$\delta\pi$ は恒等的な座標変換に伴う運動量の変化を含む。明示的な時間依存がある場合も $\partial_t Q+\{Q,H\}=0$ が保存を表す。

:::example[自由な一様モードの時間依存シフト]
$L=\dot q^2/2$ の変換 $q\mapsto q+\alpha t$ では、一次の $\delta L=\alpha\dot q$ が全微分である。$R=t,K^0=q$ なので $Q=tp-q$。実際 $\partial_tQ=p$、$\{Q,H\}=\{-q,p^2/2\}=-p$ が相殺する。さらに $\{q,Q\}=t$、$\{p,Q\}=1$ なので、速度も一定量だけ増える。例えば $q(0)=2,p(0)=3$ なら $Q=-2$ で、$t=1$ で $q=5,p=3$ を代入しても同じ値になる。$K^0$ を除くと $tp$ だけになり、運動量の変換も保存も再現しない。
:::

:::example[二つの大域変換が同じ向きを保つ場合]
同じ自由振動子で質量項がなく、定数シフト $q\mapsto q+\epsilon$ の電荷は $P=p$。時間依存シフトの電荷 $Q=tp-q$ との括弧は $\{P,Q\}=1$。両方の変換を $q,p$ に作用させた順序を入れ替えても結果は同じである。実際、括弧の右辺が定数なら任意の観測量との括弧は零である。このような定数は電荷の代数には現れるが古典的な正準変換の作用には現れない。
:::

## エネルギー運動量テンソルと Belinfante テンソル

### 並進の保存流

固定座標で $\delta\phi^a=\epsilon\partial^\nu\phi^a$ とすると、明示的な座標依存がない場合に $\delta\mathcal L=\epsilon\partial^\nu\mathcal L$ である。従って :eqref[eq-qft-noether] から
$$
T^{\mu\nu}=\Pi^\mu_a\partial^\nu\phi^a-\eta^{\mu\nu}\mathcal L,\qquad
P^\nu=\int\dd^3x\,T^{0\nu}
\label{eq-qft-canonical-tensor}
$$
を得る。特に $T^{00}=\mathcal H$ である。座標依存が明示的にある場合は、上の計算で $\partial^\nu\mathcal L$ からその明示的な微分を引くので $\partial_\mu T^{\mu\nu}=-(\partial^\nu\mathcal L)_{\rm explicit}$ となる。

Lorentz 変換の固定座標での作用を $x^\rho\partial^\sigma-x^\sigma\partial^\rho+\Sigma^{\rho\sigma}$ と書く。このとき成分を回す行列 $\Sigma$ が寄与する流を $S^{\mu\rho\sigma}=\Pi^\mu_a(\Sigma^{\rho\sigma}\phi)^a$ と置く。軌道部分を足した保存流は $x^\rho T^{\mu\sigma}-x^\sigma T^{\mu\rho}+S^{\mu\rho\sigma}$ である。その発散を展開すると $\partial_\mu S^{\mu\rho\sigma}=T^{\sigma\rho}-T^{\rho\sigma}$ が従う。

:::proposition[Belinfante の対称化]{#prop-qft-belinfante}
$$
\begin{aligned}
B^{\lambda\mu\nu}&=\tfrac12(S^{\lambda\mu\nu}+S^{\mu\nu\lambda}+S^{\nu\mu\lambda}),\\
\Theta^{\mu\nu}&=T^{\mu\nu}+\partial_\lambda B^{\lambda\mu\nu}
\end{aligned}
\label{eq-qft-belinfante}
$$
は解の上で対称かつ保存し、表面項が消えれば $P^\nu$ を変えない。
:::

:::proof
$S$ は最後の二添字で反対称だから $B^{\lambda\mu\nu}=-B^{\mu\lambda\nu}$ かつ $B^{\lambda\mu\nu}-B^{\lambda\nu\mu}=S^{\lambda\mu\nu}$ である。前者で発散は変わらず、$B^{00\nu}=0$ のため電荷の差は空間発散に限られる。後者と Lorentz 流の保存式で $\Theta^{\mu\nu}-\Theta^{\nu\mu}=0$ になる。さらに $x^\rho\Theta^{\mu\sigma}-x^\sigma\Theta^{\mu\rho}$ と元の全角運動量流の差は
$$
\partial_\lambda(x^\rho B^{\lambda\mu\sigma}-x^\sigma B^{\lambda\mu\rho})
$$
である。積の微分から出る $-B^{\rho\mu\sigma}+B^{\sigma\mu\rho}$ は $B^{\mu\rho\sigma}-B^{\mu\sigma\rho}=S^{\mu\rho\sigma}$ であり、元のスピン流と相殺する。この表面項も消える境界条件なら、全角運動量も同じである。
:::

:::example[電磁場の対称テンソル]
$F_{\mu\nu}=\partial_\mu A_\nu-\partial_\nu A_\mu$、$\mathcal L=-F^2/4$ のとき $\partial\mathcal L/\partial(\partial_\mu A_\rho)=-F^{\mu\rho}$ である。$T^{\mu\nu}$ に $\partial_\rho(F^{\mu\rho}A^\nu)$ を足すと、自由方程式 $\partial_\rho F^{\rho\mu}=0$ により $\Theta^{\mu\nu}=-F^{\mu\rho}F^\nu{}_\rho+\eta^{\mu\nu}F^2/4$ となる。例えば純電場なら $F^{0\rho}F^0{}_\rho=-|\vb E|^2$ と $F^2=-2|\vb E|^2$ から $\Theta^{00}=|\vb E|^2/2$ である。
:::

### 計量を変分して得るテンソル

平坦時空を一度背景計量 $g_{\mu\nu}$ に置き換え、$S=\int\dd^4x\sqrt{-g}\,\mathcal L$ とする。この計量符号では計量変分 $\delta S=\frac12\int\sqrt{-g}\,\Theta_{\mu\nu}\delta g^{\mu\nu}\dd^4x$ を Hilbert のテンソルの定義とする。逆計量の変分では $\delta\sqrt{-g}=-\frac12\sqrt{-g}\,g_{\mu\nu}\delta g^{\mu\nu}$。最小結合スカラーの場合、$\partial\mathcal L/\partial g^{\mu\nu}=\partial_\mu\phi\partial_\nu\phi/2$ だから $\Theta_{\mu\nu}=\partial_\mu\phi\partial_\nu\phi-g_{\mu\nu}\mathcal L$ である。平坦極限で :eqref[eq-qft-canonical-tensor] と一致する。Maxwell 場の場合は $F_{\mu\nu}$ の下添字を固定し、$F^2$ の二つの逆計量をそれぞれ変分すると同じ $F_{\mu\rho}F_\nu{}^\rho$ が出る。この結果は上で改良した Maxwell のテンソルと一致する。

:::example[波のエネルギー流]
$\phi=A\cos(\omega t-kz)$、$\omega^2=k^2+m^2$ を考える。時間平均で $\langle T^{00}\rangle=A^2(\omega^2+k^2+m^2)/4=A^2\omega^2/2$。$\langle T^{0z}\rangle=A^2\omega k/2$ なのでエネルギーの流速は $k/\omega$ となる。分散関係を微分した群速度 $\dd\omega/\dd k$ と一致する。
:::

### 保存則を保つ改良とトレース

平坦時空では任意の滑らかな局所スカラー $F$ に対し $I^{\mu\nu}=(\eta^{\mu\nu}\Box-\partial^\mu\partial^\nu)F$ が対称である。$\partial_\mu I^{\mu\nu}=\partial^\nu\Box F-\Box\partial^\nu F=0$ は運動方程式に依存しない。さらに $I^{00}=-\nabla^2 F$、$I^{0i}=\partial_t\partial_i F$ である。従って $\nabla F,\partial_t F$ の対応する表面積分が消える周期境界や十分な減衰のもとでは、全エネルギーと全運動量を変えない。一方で局所トレースには $(d-1)\Box F$ が加わる。

:::example[質量ゼロの四次理論の改良テンソル]
四次元で $V=\lambda\phi^4/4!$ とする。最小結合のテンソルのトレースは $-(\partial\phi)^2+4V$。$\Box\phi^2=2(\partial\phi)^2+2\phi\Box\phi=2(\partial\phi)^2-8V$ である。従って $F=\phi^2/6$ による改良を加えると殻上でトレースは零になる。さらに $\partial_\mu(x_\nu\Theta^{\mu\nu})=\Theta^\mu{}_\mu$ だから、この改良を使ってスケール変換の保存流が得られる。$V=m^2\phi^2/2$ なら同じ計算をして $m^2\phi^2$ が残る。
:::

:::example[外部源とのエネルギー移動]
$J(t)\phi$ の項を含むスカラー場で体積平均を取り、空間的に一様とする。場だけのエネルギー密度 $u=(\dot\phi^2+m^2\phi^2)/2$ は $\dot u=\dot\phi(\ddot\phi+m^2\phi)=J\dot\phi$ を満たす。源との結合を含む正準エネルギー密度は $u-J\phi$ なのでその微分は $-\dot J\phi$ である。例えば $J$ を一定に保って場が動く間、前者は変化し得るが後者は保存する。
:::

### エネルギーの重心とブースト電荷

対称な $\Theta^{\mu\nu}$ が保存し、全エネルギー $E>0$ と一次の空間モーメントが有限とする。$X^i=E^{-1}\int x^i\Theta^{00}\dd^3x$ と置く。$\dot E=0$、$\partial_t\Theta^{00}=-\partial_j\Theta^{j0}$ を用い、境界積分が消えると
$$
E\dot X^i=\int \Theta^{i0}\dd^3x=\int\Theta^{0i}\dd^3x=P^i .
$$
最後は対称性を使った。従って $EX^i-tP^i$ が保存する。これは軌道部分の Lorentz 流 $x^\rho\Theta^{\mu\sigma}-x^\sigma\Theta^{\mu\rho}$ のブースト成分と符号を除いて同じである。内部のスピン流を対称テンソルへ吸収すると、重心の運動も同じ $\Theta$ の成分で記述できる。

:::example[離れた二つのパケットの重心速度]
二つの波束が離れていて干渉を無視でき、一方のエネルギー $E_1=5\,\mathrm{GeV}$ と運動量 $P_1^z=4\,\mathrm{GeV}$、他方の $E_2=3\,\mathrm{GeV}$ と $P_2^z=-1\,\mathrm{GeV}$ を足すと、重心速度は $(4-1)/(5+3)=3/8$。二つの速度を単純平均した $(4/5-1/3)/2=7/30$ とは異なる。エネルギーの重心では二波束のエネルギーの違いが重みに入る。
:::

### エネルギー運動量テンソルの保存を Maxwell 方程式だけで示す

外部保存源を加えた電磁場で $\Theta_{\rm em}^{\mu\nu}=-F^{\mu\rho}F^\nu{}_\rho+\eta^{\mu\nu}F^2/4$ とする。積の微分で
$$
\partial_\mu\Theta_{\rm em}^{\mu\nu}
=-(\partial_\mu F^{\mu\rho})F^\nu{}_\rho
-F^{\mu\rho}\partial_\mu F^\nu{}_\rho+\tfrac14\partial^\nu F^2 .
$$
Bianchi 恒等式 $\partial_\mu F_{\lambda\rho}+\partial_\lambda F_{\rho\mu}+\partial_\rho F_{\mu\lambda}=0$ に $F^{\mu\rho}$ を掛ける。第一項と第三項は $\mu,\rho$ の入れ替えにより等しく、中央項は $-F^{\mu\rho}\partial_\lambda F_{\mu\rho}=-\partial_\lambda F^2/2$。従って $F^{\mu\rho}\partial_\mu F^\nu{}_\rho=\partial^\nu F^2/4$。運動方程式 $\partial_\mu F^{\mu\rho}=j^\rho$ を使うと
$$
\partial_\mu\Theta_{\rm em}^{\mu\nu}=-F^\nu{}_\rho j^\rho .
$$
これは源が電磁場から得る力密度と逆符号である。電磁場だけのエネルギーと運動量は源があると保存しないが、動的な荷電物質のテンソルに $+F^\nu{}_\rho j^\rho$ があれば両者の和が保存する。源を外部に固定した計算では後者を模型の外に置く。

:::example[電磁場のエネルギーが電流に移る符号]
$\nu=0$ なら $F^{0i}=-E^i$、$F^0{}_i=E^i$ である。従って $\partial_t\Theta^{00}_{\rm em}+\partial_i\Theta^{i0}_{\rm em}=-\vb E\cdot\vb j$。例えば局所で $\vb E=(2,0,0)\,\mathrm{GeV}^2$、$\vb j=(3,0,0)\,\mathrm{GeV}^3$ なら、電流への仕事率密度は $6\,\mathrm{GeV}^5$。逆方向の電流なら符号が反転し、仕事を受けるのは電磁場の側になる。ここで $\Theta^{00}$ だけの時間変化を求めるには、さらに空間の流束の発散を知る必要がある。
:::

### 古典のスケールカレントと改良の係数

$d>2$ 次元で正準実スカラーの無限小スケール変換を $\delta\phi=\epsilon(\Delta_\phi+x^\mu\partial_\mu)\phi$、$\Delta_\phi=(d-2)/2$ とする。運動項は $2(1+\Delta_\phi)=d$ なので、その変分は $\epsilon\partial_\mu[x^\mu(\partial\phi)^2/2]$。ポテンシャルが $V(\phi)=\lambda\phi^{2d/(d-2)}$ など、考えている領域で $\Delta_\phi\phi V'=d V$ を満たせば、同じ全微分の形になる。Noether 流は $D^\mu=x_\nu T^{\mu\nu}+\Delta_\phi\phi\partial^\mu\phi$。実際、トレースと第二項の発散を合わせると
$$
T^\mu{}_\mu+\Delta_\phi\bigl[(\partial\phi)^2+\phi\Box\phi\bigr]
=d V-\Delta_\phi\phi V'
$$
であり、運動方程式 $\Box\phi=-V'$ を用いた。この式はポテンシャルに質量項が入った場合のスケール対称性の破れも与える。

改良 $\Theta^{\mu\nu}=T^{\mu\nu}+\xi(\eta^{\mu\nu}\Box-\partial^\mu\partial^\nu)\phi^2$ を考える。トレースの追加は $\xi(d-1)\Box\phi^2$。上式の追加項を再現するには $2\xi(d-1)=\Delta_\phi$、すなわち $\xi=(d-2)/[4(d-1)]$ が必要である。従って $\partial_\mu(x_\nu\Theta^{\mu\nu})=dV-\Delta_\phi\phi V'$。二つのスケール流の差は局所全微分である。これを確認するには $I^{\mu\nu}=(\eta^{\mu\nu}\Box-\partial^\mu\partial^\nu)F$ に対し
$$
x_\nu I^{\mu\nu}=
(d-1)\partial^\mu F+
\partial_\lambda(x^\mu\partial^\lambda F-x^\lambda\partial^\mu F)
$$
と展開する。ここで $F=\xi\phi^2$ とすれば第一項が $\Delta_\phi\phi\partial^\mu\phi$ で、第二項は反対称二階テンソルの発散である。

:::example[四次元と三次元の改良係数を比べる]
四次元では $\xi=1/6$。三次元では $\Delta_\phi=1/2$、$\xi=1/8$ である。例えば $d=3$ で $V=a\phi^6+b\phi^4$ の場合、$dV-\Delta_\phi\phi V'=3(a\phi^6+b\phi^4)-(6a\phi^6+4b\phi^4)/2=b\phi^4$。六次項の寄与が消え、次元一の係数 $b$ があることによる破れが残る。これは古典の式であり、量子論では無次元と数えられた結合にも [くりこみ群](/physics/quantum-field-theory/renormalization-group/) による尺度依存が加わり得る。
:::

## 内部対称性と U(1) カレント

### 位相回転と場の電荷

二つの実場を $\varphi=(\phi^1+i\phi^2)/\sqrt2$ にまとめる。$\mathcal L=\partial_\mu\varphi^*\partial^\mu\varphi-V(|\varphi|^2)$ では $\varphi$ と $\varphi^*$ の変分を実部と虚部の変分の独立な線形結合として取れる。$\varphi\mapsto e^{-i\alpha}\varphi$ は時空点を移さず場の成分を回す :term[内部対称性]{en="internal symmetry" yomi="ないぶたいしょうせい"} である。:eqref[eq-qft-noether] に $R=-i\varphi$、$R^*=i\varphi^*$ を代入すると
$$
j^\mu=i(\varphi^*\partial^\mu\varphi-\partial^\mu\varphi^*\,\varphi)
\label{eq-qft-uone-current}
$$
になる。発散を直接取ると交差項が相殺し、$\Box\varphi=-V'(|\varphi|^2)\varphi$ とその共役により残りも消える。自由場で $\varphi=Ae^{-ipx}$ なら $j^\mu=2p^\mu|A|^2$ となり、振動の向きを逆にすると電荷の符号が反転する。

位相を位置に依存させると、微分が位相に作用する。そこで $D_\mu=\partial_\mu+iqA_\mu$、$\varphi'=e^{-iq\alpha}\varphi$、$A_\mu'=A_\mu+\partial_\mu\alpha$ と置けば、$\partial_\mu\alpha$ を含む二項が打ち消し合い $D_\mu'\varphi'=e^{-iq\alpha}D_\mu\varphi$ となる。この局所変換は [量子電磁力学](/physics/quantum-field-theory/qed/) で使う。

### 極座標で見た電荷の輸送

$\varphi=\rho e^{-i\theta}/\sqrt2$ と書ける領域 $\rho>0$ では
$$
\mathcal L=\tfrac12(\partial\rho)^2+\tfrac12\rho^2(\partial\theta)^2-V(\rho^2/2),
\qquad j^\mu=\rho^2\partial^\mu\theta .
$$
位相の変分で $\partial_\mu(\rho^2\partial^\mu\theta)=0$、動径の変分で $\Box\rho-\rho(\partial\theta)^2+\rho V'=0$ となる。従って位相勾配を外から固定しても、上の動径方程式を満たさなければ全体の解にはならない。$\rho=0$ の点では位相座標を使えないので、零点を含む配置では元の複素場で方程式を課す。

:::example[一定の大きさで回転する配置]
$V(s)=m^2s+\lambda s^2/2$、$\rho=\rho_0>0$、$\theta=\omega t-kz$ を代入すると $\omega^2-k^2=m^2+\lambda\rho_0^2/2$ である。$m=1\,\mathrm{GeV}$、$\rho_0=2\,\mathrm{GeV}$、$\lambda=4$、$k=4\,\mathrm{GeV}$ なら $\omega=5\,\mathrm{GeV}$ の枝で $j^0=20\,\mathrm{GeV}^3$、$j^z=16\,\mathrm{GeV}^3$ である。電流が一定なので連続の式を満たす。
:::

:::example[位相対称性を破る混合質量]
$\Delta\mathcal L=-\kappa(\varphi^2+\varphi^{*2})/2$ を加える。$\delta\varphi=-i\alpha\varphi$ で定数 $\alpha$ の変分は $\delta\mathcal L=i\kappa\alpha(\varphi^2-\varphi^{*2})$。従って殻上で $\partial_\mu j^\mu=-2\kappa\Im(\varphi^2)$ である。実成分では二つの質量二乗にそれぞれ $+\kappa,-\kappa$ が足される。例えば $\varphi=(1+i)\,\mathrm{GeV}$、$\kappa=3\,\mathrm{GeV}^2$ なら発散は $-12\,\mathrm{GeV}^4$ であり、位相回転の電荷は保存しない。
:::

### 背景ポテンシャルと物質の保存電流

複素場に固定した外部ポテンシャルを結合し、$\mathcal L=(D_\mu\varphi)^*D^\mu\varphi-V(s)$、$s=|\varphi|^2$、$D_\mu=\partial_\mu+iqA_\mu$ とする。$D_\mu^*=\partial_\mu-iqA_\mu$ と書くと、場の方程式は $D_\mu D^\mu\varphi+V'(s)\varphi=0$ である。実際、$(D_\mu^*\delta\varphi^*)D^\mu\varphi$ の通常微分を部分積分すると、二項を合わせて $-\delta\varphi^*D_\mu D^\mu\varphi$。$\varphi^*$ の式はこれの複素共役である。定数位相に関する正準流は
$$
j_A^\mu=i\left(\varphi^*D^\mu\varphi-(D^\mu\varphi)^*\varphi\right)
=j^\mu-2q A^\mu|\varphi|^2 .
$$
$\partial_\mu$ で発散を取ると、各積は電荷零なので共変微分を積の法則で使える。一次微分の交差項は逆符号で消え、残りは $i[\varphi^*D^2\varphi-(D^2\varphi)^*\varphi]$ である。両方に運動方程式を入れると零。別に $A_\mu$ を変分すると、運動項の変分は $-qj_A^\mu\delta A_\mu$。従って動的な電磁場を加えたとき Maxwell の源は $qj_A^\mu$ である。背景がある場合に自由場の $j^\mu$ だけを使うと、後半の項を失う。

:::example[一定の電位を変えたときの見かけの周波数]
$\vb A=0,A^0=A_c$ を一定とし、$\varphi=a e^{-i\omega t+i\vb k\cdot\vb x}$ の自由な質量項だけを残す。$D_0\varphi=-i(\omega-qA_c)\varphi$ なので分散式は $(\omega-qA_c)^2=|\vb k|^2+m^2$。従って正の運動学的周波数の枝で $j_A^0=2\sqrt{|\vb k|^2+m^2}|a|^2$。局所変換 $\alpha=-A_c t$ によって $A_c$ を零へ移すと $\varphi$ の時間振動数も $\omega-qA_c$ へ移り、電流は同じになる。例えば $m=3\,\mathrm{GeV}$、$|\vb k|=4\,\mathrm{GeV}$、$q A_c=2\,\mathrm{GeV}$ なら $\omega=7\,\mathrm{GeV}$ だが $j_A^0=10\,\mathrm{GeV}\,|a|^2$。ゲージ変換を制限する境界や周期条件がある場合は、一定ポテンシャルを同じ手順で除けるかを別途調べる。
:::

### 固定電荷で一様な場を最小化する

位相変数を $\varphi=\rho e^{-i\theta}/\sqrt2$ とし、空間平均だけを扱う。$\rho>0$、位相の電荷密度を $n=\rho^2\dot\theta$ とする。固定した $n$ では
$$
u=\tfrac12\dot\rho^2+\frac{n^2}{2\rho^2}+V(\rho^2/2)
$$
となる。これは $\dot\theta$ を電荷で置き換えたエネルギーであり、Lagrangian に代入してそのまま変分すると同じ式にならない。固定電荷で座標を消去するには Routh の変換 $\mathcal R=\mathcal L-n\dot\theta$ を使う。実際、$\delta_\rho\mathcal R$ で生じる $\delta_\rho\dot\theta$ の係数は $\partial\mathcal L/\partial\dot\theta-n=0$。従って $\mathcal R=\dot\rho^2/2-n^2/(2\rho^2)-V$ の動径方程式は、元の動径方程式に $\dot\theta=n/\rho^2$ を代入した式と一致する。

:::example[自由な質量付き場で電荷に比例するエネルギー]
$V=m^2\rho^2/2$、$m>0,n\ne0$ とする。静的動径の極小は $-n^2/\rho^3+m^2\rho=0$、すなわち $\rho^2=|n|/m$。$n^2/(2\rho^2)$ と $m^2\rho^2/2$ が等しく、極小エネルギー密度は $m|n|$ である。また $\dot\theta=\sgn(n)m$。例えば $m=2\,\mathrm{GeV},n=8\,\mathrm{GeV}^3$ なら $\rho=2\,\mathrm{GeV}$、$u=16\,\mathrm{GeV}^4$。$\rho$ だけを零へ近づけてエネルギーを下げようとすると、固定電荷の項が発散する。
:::

## 自然単位系と次元解析

### 作用の位相と質量次元

作用は $S/\hbar$ が位相に現れるので自然単位系で無次元である。質量次元を $[X]$ と書き、$[x]=-1$、$[\partial]=1$ とすると $d$ 次元時空で $[\mathcal L]=d$ である。係数が無次元の運動項から $[\phi]=(d-2)/2$、一次微分の Dirac 場は $[\psi]=(d-1)/2$、Maxwell 場は $[A]=(d-2)/2$ になる。

:::example[高次の相互作用]
四次元で $\lambda_n\phi^n/n!$ の係数は $[\lambda_n]=4-n$ である。四次項は無次元、六次項は $\mathrm{GeV}^{-2}$ などで測る。係数が $1/M^2$ で振幅にエネルギーの二乗が掛かる場合、$E/M=0.1$ なら無次元の補正因子は $0.01$ である。
:::

例えば $[T^{00}]=4$ と空間体積の次元 $-3$ から $[H]=1$ になる。$\hbar c\simeq0.1973\,\mathrm{GeV\,fm}$ なので $1\,\mathrm{GeV}^{-1}$ を長さに戻すと約 $0.197\,\mathrm{fm}$ である。エネルギー密度 $1\,\mathrm{GeV}^4$ を体積あたりの単位で表すと約 $130\,\mathrm{GeV/fm^3}$ となる。

### 無次元化して近似の条件を読む

$m>0$ の四次理論を考える。$\widehat x=mx$、$\varphi=\phi/m$、$\widehat J=J/m^3$ とすると、作用の測度 $m^{-4}$ が密度の $m^4$ を打ち消す。その結果、方程式は $(\widehat\Box+1)\varphi+\lambda\varphi^3/6=\widehat J$ になる。一様場で非線形項と質量項の絶対値の比は $\abs{\lambda}\varphi^2/6$ である。従って $\lambda$ の値だけでは線形化の誤差は判断できず、場の振幅も併せて指定する。

:::example[局所演算子の係数]
四次元で $\phi^2(\partial\phi)^2$ は次元 $6$ なので係数の次元は $-2$ である。背景 $\phi_0$ に小さな揺らぎ $\eta$ を加えると、$c\phi^2(\partial\phi)^2/M^2$ の二次部分は $c\phi_0^2(\partial\eta)^2/M^2$ である。$\phi_0/M=1/5$、$c=2$ なら、標準の運動項の係数 $1/2$ からの相対変化は $4/25$ である。勾配が小さくても背景が大きいとこの補正は抑えられない。
:::

:::example[次元を変えたときの四次結合]
$d=3$ では $[\phi]=1/2$ だから $[\lambda_4]=1$ である。運動量尺度 $\mu$ で使う無次元の比は $\lambda_4/\mu$ であり、$\lambda_4=0.2\,\mathrm{GeV}$ なら $\mu=2\,\mathrm{GeV}$ で $0.1$、$\mu=0.1\,\mathrm{GeV}$ で $2$ となる。これは次元だけで決まる変化である。量子補正による尺度依存は [くりこみ群](/physics/quantum-field-theory/renormalization-group/) で区別する。
:::

### 相互作用の次元と比較する尺度を決める

$d$ 次元で局所項 $c\,\mathcal O$ を密度に加え、$\mathcal O$ が正準次元 $\Delta$ を持つとする。唯一の物理尺度が $E>0$ で、場の振幅や微分もその尺度で数えられるなら、無次元変数を $x=\widehat x/E$、$\phi=E^{[\phi]}\widehat\phi$ と置いて
$$
\int \dd^dx\,c\,\mathcal O(x)
=\bigl(cE^{\Delta-d}\bigr)\int\dd^d\widehat x\,\widehat{\mathcal O}(\widehat x)
$$
となる。この係数が小さいとき、その作用への寄与を展開できる候補になる。ただし複数の尺度がある場合は、微分の $E$ と大きな一定背景の値を別々に代入する必要がある。例えば $\phi^2(\partial\phi)^2/M^2$ の四脚過程では微分二つが $E^2/M^2$ の抑制を与える。一様背景の上の二点過程ではこの演算子の二つの $\phi$ が背景に置き替わるため、標準運動項との比は背景二乗で決まる。

:::example[次元五の相互作用の比較]
四次元でスカラー $\phi$ と Dirac 場 $\psi$ を持つとする。$\phi\,\bar\psi\psi$ の次元は $1+3=4$、$\phi^2\bar\psi\psi$ は $2+3=5$。後者に $c/M$ を掛ければ正しい密度になる。共通尺度 $E$ では後者の無次元係数は $cE/M$。例えば $c=2$、$M=10\,\mathrm{GeV}$、$E=0.5\,\mathrm{GeV}$ なら $0.1$。一方、$\phi$ を一定の $v=2\,\mathrm{GeV}$ の周りで展開した場合、質量に入る補正は $cv^2/M=0.8\,\mathrm{GeV}$ である。低い外部運動量での減衰だけから、背景での質量効果まで小さいとは言えない。
:::

:::example[作用の全体係数と量子展開の区別]
$d=4$ で $\lambda>0$、$\mathcal L=(\partial\phi)^2/2-m^2\phi^2/2-\lambda\phi^4/4!$ とする。$\phi=m\chi/\sqrt\lambda$、$\widehat x=mx$ と取れば $S=\lambda^{-1}\int\dd^4\widehat x[(\widehat\partial\chi)^2/2-\chi^2/2-\chi^4/4!]$。古典方程式では全体係数 $1/\lambda$ が消え、$\widehat\Box\chi+\chi+\chi^3/6=0$ となる。一方、量子論の位相 $e^{iS}$ では $1/\lambda$ が残る。例えば $\lambda=1/100$ なら、同じ無次元配置 $\chi$ の作用差は括弧の積分の $100$ 倍になる。小さい $\lambda$ に対して古典の非線形効果が消えるかどうかは、どちらの場の振幅を固定したかで決まる。
:::

### 長さと時間を SI 単位へ戻す

自然単位系で質量とエネルギーの次元を同一視する。エネルギー $E$ に対し時間発展の位相は $Et$ であるため、SI の時間は $\hbar/E$ を尺度とする。空間の位相は自然単位系で $|\vb p|x$。相対論的な高運動量で $E=c|\vb p|$ を戻すと長さは $\hbar c/E$。$m$ を質量ではなく静止エネルギーの単位で与えた場合の Compton 長も $\hbar c/m$ である。例えば $E=2\,\mathrm{GeV}$ に対し $1/E=0.5\,\mathrm{GeV}^{-1}$ の長さは、$\hbar c\simeq0.1973\,\mathrm{GeV\,fm}$ から $0.0987\,\mathrm{fm}$。$t$ を戻すときは同じ式に長さの単位を代入せず、$\hbar\simeq6.58\times10^{-25}\,\mathrm{GeV\,s}$ を使って $3.29\times10^{-25}\,\mathrm{s}$ である。長さを時間に換算した $\ell/c$ とも一致する。

:::example[電荷密度と全電荷]
四次元で内部 $U(1)$ の $j^0$ は $\varphi^*\partial_t\varphi$ より質量次元三。自然単位系で $j^0=2\,\mathrm{GeV}^3$、体積が $3\,\mathrm{GeV}^{-3}$ の一様配置では $Q=6$。長さの換算係数を $L_0=0.1973\,\mathrm{fm}$ とすると、密度は $2/L_0^3\,\mathrm{fm}^{-3}$、体積は $3L_0^3$ であり、積は同じ $6$。これを電気電荷に直すには一粒子あたりの電荷を掛ける。密度だけを換算して体積の数値をそのまま使うと、保存量自体を変えてしまう。
:::

### 同じ質量次元でも項の大きさが変わる例

二つの次元六のスカラー項 $c_1\phi^6/M^2$ と $c_2\phi^2(\partial\phi)^2/M^2$ を比較する。振幅 $\Phi$、変化の長さ $L$ を別の尺度として固定すれば、前者は $c_1\Phi^6/M^2$、後者は $c_2\Phi^4/(M^2L^2)$ の程度である。従って比は $c_1\Phi^2L^2/c_2$ であり、両項が同じ次元だからといって同程度ではない。四次ポテンシャルとの比較では前者に $\Phi^2/M^2$、運動項との比較では後者にも $\Phi^2/M^2$ が掛かる。ところが固定外脚の散乱や一様背景の小振動では、同じ演算子に入る微分の位置が異なるため、比較する二つの基準項を先に決める必要がある。

:::example[緩い勾配と小さい振幅を独立に変える]
$c_1=c_2=1$、$\Phi=0.1M$ とする。$1/L=0.1M$ なら二項の比は一、$1/L=0.01M$ なら比は $100$ になる。一方、いずれの場合も標準の運動項に対する二微分の次元六項の係数比は、$1/2$ を含めて $2\Phi^2/M^2=0.02$ である。従って勾配を弱めて二つの補正間の比を大きくしても、運動項に対する二微分補正の抑制は同じである。
:::

:::note
ここでの質量次元は自然単位系における単位の換算である。量子相関の尺度変換で現れる異常次元は、くりこみ条件を定めてから計算する。
:::

## 演習問題

::::exercise{#exr-qft-classical-fields-1 level=basic}
$\mathcal L=(\partial\phi)^2/2-m^2\phi^2/2-\lambda\phi^4/4!$ の運動方程式を求めよ。

:::solution
$\Pi^\mu=\partial^\mu\phi$、$\partial\mathcal L/\partial\phi=-m^2\phi-\lambda\phi^3/6$ なので $\Box\phi+m^2\phi+\lambda\phi^3/6=0$ である。
:::
::::

::::exercise{#exr-qft-classical-fields-2 level=basic}
$m=4\,\mathrm{GeV}$、$|\vb k|=3\,\mathrm{GeV}$ の自由モードの正の周波数を求めよ。

:::solution
平面波の $\Box$ 固有値は $-\omega^2+k^2$ なので $\omega=\sqrt{3^2+4^2}\,\mathrm{GeV}=5\,\mathrm{GeV}$ である。
:::
::::

::::exercise{#exr-qft-classical-fields-3 level=basic}
$\varphi=Ae^{-i\omega t}$ の電荷密度を求めよ。

:::solution
$\dot\varphi=-i\omega\varphi$ と共役を :eqref[eq-qft-uone-current] に代入すると $j^0=2\omega|A|^2$ になる。
:::
::::

::::exercise{#exr-qft-classical-fields-4 level=basic}
$d=3$ のスカラー理論で $\phi^6$ の係数の次元を求めよ。

:::solution
$[\phi]=1/2$ であり、積の次元は $3$ である。従って係数は無次元である。
:::
::::

::::exercise{#exr-qft-classical-fields-5 level=advanced}
自由な質量ゼロのスカラー場で $T^{\mu\nu}+a(\eta^{\mu\nu}\Box-\partial^\mu\partial^\nu)\phi^2$ のトレースを零にする $a$ を四次元で求めよ。

:::solution
トレースは $-(\partial\phi)^2+3a\Box\phi^2$ である。$\Box\phi=0$ と積の微分から $\Box\phi^2=2(\partial\phi)^2$ なので $a=1/6$ である。改良項の発散は微分の可換性で消える。
:::
::::

::::exercise{#exr-qft-classical-fields-6 level=advanced}
外部源 $J$ のある実場で $\dot H$ を求めよ。

:::solution
連鎖律で $(\delta H/\delta\phi)\dot\phi+(\delta H/\delta\pi)\dot\pi=-\dot\pi\dot\phi+\dot\phi\dot\pi=0$ になる。明示的依存は $-J\phi$ だけなので $\dot H=-\int\dot J\phi\,\dd^3x$ である。
:::
::::

::::exercise{#exr-qft-classical-fields-7 level=advanced}
$P^i$ が $\pi$ を並進することを示せ。

:::solution
$P^i=\int\pi\partial^i\phi\,\dd^3x$ の $\phi$ 変分を部分積分すると $-\partial^i\pi$ が得られる。従って :eqref[eq-qft-field-poisson] から $\{\pi,P^i\}=\partial^i\pi$ である。
:::
::::

::::exercise{#exr-qft-classical-fields-8 level=advanced}
$\Delta\mathcal L=-\mu^2(\varphi^2+\varphi^{*2})/2$ を加えたとき $\partial_\mu j^\mu$ を求めよ。

:::solution
定数変換に対して $\delta(\Delta\mathcal L)=i\alpha\mu^2(\varphi^2-\varphi^{*2})$ である。一般の変分式より解の上では $\delta\mathcal L=\alpha\partial_\mu j^\mu$ なので、発散は $i\mu^2(\varphi^2-\varphi^{*2})$ である。
:::
::::

::::exercise{#exr-qft-classical-fields-9 level=advanced}
時空を $1+1$ 次元、空間を半直線 $x\ge0$ とする。$m>0$ の実スカラー場に境界条件 $\partial_x\phi(0)=b\phi(0)$ を課す。$b<0$ で指数的に局在する正規モードの周波数を求め、線形に不安定になる条件を示せ。$m=3\,\mathrm{GeV}$、$b=-4\,\mathrm{GeV}$ での増幅率も求めよ。

:::solution
時間依存を $e^{-i\omega t}$ とし、空間側を $e^{-\kappa x}$、$\kappa>0$ とする。境界条件は $-\kappa=b$ なので $b<0$ が必要である。一方、Klein–Gordon 方程式には $-\omega^2-\kappa^2+m^2=0$ と入る。従って $\omega^2=m^2-b^2$。$|b|<m$ なら実周波数、$|b|=m$ なら零周波数であり、このとき時間方程式には定数解と時間に比例する解がある。$|b|>m$ なら指数関数 $e^{\pm\sqrt{b^2-m^2}\,t}$ が現れ、成長側が線形不安定性を与える。指定した値では増幅率は $\sqrt7\,\mathrm{GeV}$ である。

この結論は Hamiltonian からも確かめられる。$\phi=q(t)e^{bx}$ と置くと $\int_0^\infty e^{2bx}\dd x=-1/(2b)$。運動項の係数と同じ重みで $\dot q^2/2+(b^2+m^2)q^2/2$ が生じ、さらに境界の Hamiltonian $bq^2/2$ が加わる。これを同じ重みの中に入れれば $(m^2-b^2)q^2/2$ となる。境界エネルギーを落とすと不安定な範囲を読み違える。
:::
::::

## 確認テスト

::::quiz
固定した境界値で作用が停留するための条件はどれか。

- [ ] $\mathcal L=0$
- [x] $E_a=0$
- [ ] $\Pi^\mu_a=0$

:::explanation
境界変分は消え、内部に台を持つ任意の変分の係数が $E_a$ になる。その係数が零であることが必要十分である。
:::
::::

::::quiz{answer="2" tol="1e-6"}
四次元での $m^2$ の質量次元を答えよ。

:::explanation
$[\phi]=1$ なので $[m^2]+2=4$ から $[m^2]=2$ となる。
:::
::::

::::quiz
表面流束が零のとき、局所保存式から何が従うか。

- [ ] 各点の $j^0$ が一定
- [x] 空間積分した $j^0$ が一定
- [ ] 各点の $\vb j$ が零

:::explanation
連続の式を空間で積分すると $\dot Q$ が境界流束の負号になる。流束が零なので $Q$ が一定である。
:::
::::

::::quiz{answer="6" tol="1e-6"}
GeV を基準にした自然単位で $\varphi=e^{-3it+3ix}$ の $j^0$ の数値を答えよ。

:::explanation
$\dot\varphi=-3i\varphi$ であり、$j^0=2\cdot3|\varphi|^2=6$ である。
:::
::::

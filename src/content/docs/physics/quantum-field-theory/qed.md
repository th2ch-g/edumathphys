---
title: "量子電磁力学"
description: "QED のゲージ対称性から Feynman 則を作り、対生成、光子散乱、同種粒子散乱を計算する。"
notation:
  - symbol: '\alpha'
    meaning: 微細構造定数 e²/(4π)
  - symbol: 'D_\mu'
    meaning: 共変微分
  - symbol: 'S_F,D_F'
    meaning: フェルミオンと光子の伝播関数
---

## この章の目標

- ゲージ変換に対する電流と相互作用の関係を求める。
- 対生成と Compton 散乱の断面積を計算する。
- 脚の交換、束縛状態、スカラー電荷に応じた規則を使う。
- 光子の偏極を運動量で置き換えたときの恒等式を示す。

自然単位系 $\hbar=c=1$、計量 $\eta=\diag(+1,-1,-1,-1)$ を用いる。電荷の大きさを $e>0$ とし、電子の電荷は $-e$、$\alpha=e^2/(4\pi)$ とする。

## QED の Lagrangian とゲージ不変性

### 局所位相と電磁ポテンシャル

各 Dirac 場の符号つき電荷を $q_f$ とし、$D_\mu=\partial_\mu+iq_f A_\mu$ を用いる。
$$
\mathcal L_{\rm QED}=-\tfrac14F_{\mu\nu}F^{\mu\nu}+\sum_f\bar\psi_f(i\not D-m_f)\psi_f
\label{eq-qft-qed-lagrangian}
$$
で $A_\mu'=A_\mu+\partial_\mu\omega$、$\psi_f'=e^{-iq_f\omega}\psi_f$ とすれば $F$ は不変で $D_\mu'\psi_f'=e^{-iq_f\omega}D_\mu\psi_f$。$D$ の微分と $A$ の変化の $\partial_\mu\omega$ が相殺する。展開すると相互作用は $-\sum_f q_f\bar\psi_f\gamma^\mu\psi_f A_\mu$ である。$A_\mu$ の式は $\partial_\nu F^{\nu\mu}=j^\mu$ であり、$j^\mu=\sum_f q_f\bar\psi_f\gamma^\mu\psi_f$ である。電流の保存条件 $\partial_\mu j^\mu=0$ は反対称性と Dirac 方程式のいずれからも従う。

### 局所電荷と外部源

$A_\mu$ の源を $J^\mu$ とし、作用に $-\int J^\mu A_\mu\dd^4x$ を加える。コンパクトな台をもつ $\omega$ なら変分は $+\int\omega\partial_\mu J^\mu\dd^4x$ なので、任意の $\omega$ に対する不変性と $\partial_\mu J^\mu=0$ は同値である。例えば $k\ne0$ で $\rho=\rho_0\cos(kz-\omega t)$ に対し $J^z=(\omega/k)\rho$、他の空間成分が零なら、時間微分と縦方向の微分が相殺する。$\omega,k$ を独立に固定して $J^z=0$ とする源は $\omega=0$ でなければ使えない。

### Dirac 方程式と保存則の独立な導出

物質場とその共役を独立に変分すると、$(i\gamma^\mu\partial_\mu-q_f\gamma^\mu A_\mu-m_f)\psi_f=0$ と $i(\partial_\mu\bar\psi_f)\gamma^\mu+q_f\bar\psi_f\gamma^\mu A_\mu+m_f\bar\psi_f=0$ を得る。後式に右から $\psi_f$、前式に左から $\bar\psi_f$ を掛けて足すと、質量項と $A_\mu$ を含む項が消え、$\partial_\mu(\bar\psi_f\gamma^\mu\psi_f)=0$ である。従って各フレーバーの電流は物質方程式だけから保存する。一方、Maxwell 方程式に発散を取ると、対称な二微分と反対称な $F^{\mu\nu}$ の縮約により全電流の保存が必要と分かる。この整合性は、局所位相と $A_\mu$ の変換が連動することの帰結である。

:::example[ゲージ変換した自由な平面波]
$A_\mu=0$ の自由解を $\psi=e^{-ipx}u(p)$ とする。滑らかな $\omega$ から $A_\mu'=\partial_\mu\omega$ と $\psi'=e^{-iq\omega}e^{-ipx}u(p)$ を作ると、$iD'_\mu\psi'=p_\mu\psi'$ である。通常微分の固有値は位置によって $p_\mu+q\partial_\mu\omega$ になるが、運動方程式で使うのは共変微分である。$F_{\mu\nu}'=0$ なのでこの変換によって物理的な電磁場は生じない。
:::

### Gauss 則から静電項を取り出す

Coulomb ゲージ $\boldsymbol\nabla\cdot\vb A=0$ では、電荷密度を $\rho=j^0$ として Gauss 則は $-\nabla^2 A^0=\rho$。境界で $A^0$ が減衰する解を選ぶと $A^0(\vb x)=\int\dd^3y\,\rho(\vb y)/(4\pi|\vb x-\vb y|)$ である。静電エネルギーは $(1/2)\int(\nabla A^0)^2=(1/2)\int\rho A^0$。従って二つの異なる点電荷の交差項は
$$
V=\frac{q_1q_2}{4\pi r}.
$$
二つの組の順番がそれぞれ $1/2$ を受けて一つにまとまる。自己項は点電荷では発散し、二つの電荷を離すことによるエネルギーの差には含めない。例えば $q_1=-e,q_2=Ze$ なら $V=-Z\alpha/r$ で、静電場だけから交換振幅と同じ符号を得る。

光子の自由場に質量を手で加えると、この導出の Poisson 方程式が逆二乗の Coulomb 核から質量付きの核に移る。一方、局所位相の変換で $A_\mu A^\mu$ は変わるので、元の QED ではそのまま許される項ではない。長距離の $1/r$ とゲージ不変な質量ゼロの運動項は、同じ自由光子の小運動量での性質を別の観測量から見たものである。

## QED の Feynman 則

頂点で $\bar\psi,\psi,A^\mu$ を一つずつ縮約すると $-iq_f\gamma^\mu$ が残る。内部線は :eqref[eq-qft-dirac-feynman] と :eqref[eq-qft-photon-propagator] を使い、外部光子は入射なら $\epsilon^\mu$、出射なら $\epsilon^{*\mu}$ を掛ける。閉フェルミオン線ごとに負号を付ける。各頂点で運動量が保存し、独立な内部運動量は $\dd^4\ell/(2\pi)^4$ で積分する。例えば一光子交換の電荷依存性は異なる場の符号つき電荷の積である。静的な非相対論的極限なら $\widetilde V=q_1q_2/\vb q^2$、すなわち同符号電荷は斥力を受ける。

### 一光子交換を二つのゲージで計算する

電子線の端の電流は $J^\mu=\bar u(p')\gamma^\mu u(p)$。$k=p'-p$ の縮約は $\bar u(p')(\not p'-\not p)u(p)=0$ となる。ここで左からの殻方程式と右からの殻方程式が同じ質量 $m$ を与えることを使った。従って :eqref[eq-qft-photon-propagator] の縦部分は各電流の向きに依らず消える。例えば静止した二粒子を極低速で近似すると $J^0\simeq2m$、空間電流は運動量に比例して抑制される。$k^2=-|\vb k|^2$ で残る振幅を $-4m_1m_2$ で割ると $\widetilde V=q_1q_2/|\vb k|^2$。逆変換で $V=q_1q_2/(4\pi r)$ になる。電荷の単位は Heaviside–Lorentz 規約である。

### 運動量の向きと対生成電流の保存

一つの頂点に入る電子と陽電子が $p_1,p_2$ を運ぶと、光子の運動量はその和 $k=p_1+p_2$ である。粒子の殻方程式 $\not p_1u_1=m u_1$ と反粒子の殻方程式 $\bar v_2\not p_2=-m\bar v_2$ を使うと
$$
k_\mu\bar v_2\gamma^\mu u_1
=\bar v_2(\not p_1+\not p_2)u_1=0.
$$
散乱電流では運動量の差と同じ質量符号を使ったが、対消滅では和と逆の質量符号を使う。従って、対生成でも内部光子の縦方向は消える。例えば実の光子一個に終状態を制限した場合、この電流が作る四元運動量は $k^2\ge4m^2>0$ なので終状態の $k^2=0$ と両立しない。仮想光子を通じた二粒子生成ではこの運動学の制約はない。

### 外部線から見た伝播関数の極

外部光子に横偏極を掛ける規則は、横方向の伝播関数の一粒子極の係数を選ぶ操作である。自由光子について、$k=(\omega,0,0,\omega)$ の横基底を $\epsilon_1=(0,1,0,0),\epsilon_2=(0,0,1,0)$ と選ぶ。ゲージパラメータに比例する分子は $k\cdot\epsilon_r=0$ で落ち、$-\eta$ と偏極の縮約は $\delta_{rs}$ になる。従って同じ極に二つの独立な偏極があり、入射無偏極ではその二つを平均する。出射では測定しない偏極を和にする。閉線の符号と同様に、入射平均は振幅を組み立てた後、二乗に対して行う。

### 電流を保つ二つのゲージの明示的な比較

運動量 $k=(k^0,\vb k)$ の交換に二つの保存電流 $J,J'$ を結ぶ。$\vb k\ne0$ とし、$\vb J=\vb J_T+J_L\widehat{\vb k}$ に分ける。保存則により $J_L=k^0 J^0/|\vb k|$。Feynman ゲージの縮約は
$$
\frac{-i}{k^2+i0}(J^0J'^0-\vb J\cdot\vb J')
=\frac{i}{|\vb k|^2}J^0J'^0+
\frac{i}{k^2+i0}\vb J_T\cdot\vb J'_T .
$$
第一項の分子の $(1-(k^0)^2/|\vb k|^2)$ は $-k^2/|\vb k|^2$ であり、分母の $k^2$ と分布として相殺する。これで Coulomb ゲージの瞬時項と横光子項の和を得た。例えば $k^0=0$、空間電流が零なら第一項だけである。一方、$J^0=J'^0=0$ で横電流を流す源には第二項だけが寄与する。この二つの例は同じ規格化での異なる交換を表す。保存しない任意の源についてこの分解の二つの側を等置することはできない。

## e⁺e⁻→μ⁺μ⁻ 散乱

入射を $e^-(p_1),e^+(p_2)$、出射を $\mu^-(p_3),\mu^+(p_4)$ とする。両電荷を $-e$ とすれば $\mathcal M=e^2[\bar v_2\gamma^\mu u_1][\bar u_3\gamma_\mu v_4]/s$。$s\gg m_\mu^2,m_e^2$ で質量を落とすと、入射平均つきの二乗は
$$
\overline{|\mathcal M|^2}=\frac{e^4}{4s^2}
\tr(\not p_2\gamma^\mu\not p_1\gamma^\nu)
\tr(\not p_3\gamma_\mu\not p_4\gamma_\nu)
=2e^4\frac{t^2+u^2}{s^2}.
$$
:eqref[eq-qft-four-gamma-trace] の各トレースを $4(a^\mu b^\nu+a^\nu b^\mu-\eta^{\mu\nu}a\cdot b)$ と置くと積は $32[(p_2\cdot p_3)(p_1\cdot p_4)+(p_2\cdot p_4)(p_1\cdot p_3)]$ である。$t=-s(1-\cos\theta)/2$、$u=-s(1+\cos\theta)/2$ を使い
$$
\frac{\dd\sigma}{\dd\Omega}=\frac{\alpha^2}{4s}(1+\cos^2\theta),\qquad
\sigma=\frac{4\pi\alpha^2}{3s}.
\label{eq-qft-muon-cross-section}
$$
例えば同じ近似で $\sqrt{s}$ を二倍にすると全断面積は四分の一になる。入射質量だけ落とす場合は末脚のトレースで質量項 $-4m_\mu^2\eta_{\mu\nu}$ が加わる。$\beta=\sqrt{1-4m_\mu^2/s}$ とすると断面積は $\alpha^2\beta[1+\cos^2\theta+(1-\beta^2)\sin^2\theta]/(4s)$ であり、閾値では位相空間の $\beta$ が効く。

### 閾値から連続につながる全断面積

$\int\cos^2\theta\dd\Omega=4\pi/3$、$\int\sin^2\theta\dd\Omega=8\pi/3$ を使うと
$$
\sigma(s)=\frac{4\pi\alpha^2}{3s}\,\beta\left(1+\frac{2m_\mu^2}{s}\right)
$$
を得る。終状態の運動量 $p_f=\sqrt{s}\beta/2$ は閾値で消えるが、$1+2m_\mu^2/s$ は $3/2$ に近づく。例えば $s=16m_\mu^2$ なら質量ゼロの全断面積に対する比は $(\sqrt3/2)(9/8)=9\sqrt3/16$。微分断面積ではこの平均比を各角度に一律に掛けることはできない。$\theta=0$ では比は $\beta$、$\theta=\pi/2$ では $\beta(2-\beta^2)$ である。

### 入射と出射の質量を両方保つ

電子の質量を $m_e$、ミュー粒子の質量を $m_\mu$ とし、$\sqrt s>2\max(m_e,m_\mu)$ を仮定する。重心系では $E=\sqrt s/2$、$\beta_i=\sqrt{1-m_e^2/E^2}$、$\beta_f=\sqrt{1-m_\mu^2/E^2}$ である。四つのスピンを平均した振幅二乗を計算するため、各トレースを $4$ で割ったものを $L_i^{\mu\nu}$ と $L_f^{\mu\nu}$ とする。入射側は
$$
L_i^{\mu\nu}=p_2^\mu p_1^\nu+p_2^\nu p_1^\mu
-\eta^{\mu\nu}(p_1\cdot p_2+m_e^2).
$$
重心系で $p_1\cdot p_2+m_e^2=2E^2$。時間成分は全電流の保存から零であり、空間成分は $2E^2(\delta^{jk}-\beta_i^2 \widehat z^j\widehat z^k)$。出射側は $\widehat z$ を $\widehat{\vb p}_3$ に、$\beta_i$ を $\beta_f$ に替えた形になる。二つを縮約して $4e^4/s^2$ を掛けると、振幅二乗は $e^4[3-\beta_i^2-\beta_f^2+\beta_i^2\beta_f^2\cos^2\theta]$。従って
$$
\frac{\dd\sigma}{\dd\Omega}
=\frac{\alpha^2}{4s}\frac{\beta_f}{\beta_i}
 [3-\beta_i^2-\beta_f^2+\beta_i^2\beta_f^2\cos^2\theta].
$$
例えば $\beta_i\to1$ で前節の有限な $m_\mu$ の結果に戻る。$\beta_f\to0$ では角度項が消え、位相空間が終状態の率を抑える。ただし、電荷を持つ終状態の速度が $\alpha$ と同程度なら Coulomb 交換を繰り返す寄与が大きくなり、木レベルの閾値則の近傍では別の近似が必要である。

:::example[角度分布から速度因子を分ける]
$\beta_i=1$、$\beta_f^2=1/2$ の場合、前方と直角方向の微分断面積の比は $2/(2-\beta_f^2)=4/3$。共通の位相空間の比 $\beta_f/\beta_i$ はここで消える。一方、絶対値を求めるにはその比が必要である。角分布の形の検査と規格化の検査を別々に行うと、位相空間の欠落を見つけやすい。
:::

### 仮想光子の偏極を前方と直角で比べる

入射電子の質量を落とし、ビームを $z$ 軸とすると $L_i^{zz}=0$、$L_i^{xx}=L_i^{yy}=s/2$。出射軸を $xz$ 平面に置く。$\theta=0$ では $L_f^{xx}=L_f^{yy}=s/2$ で、ミュー粒子質量の効果は位相空間だけに残る。$\theta=\pi/2$ では $L_f^{xx}=(s/2)(1-\beta^2)$、$L_f^{yy}=s/2$。従って角分布の直角値には $1-\beta^2=4m_\mu^2/s$ に比例する別の成分が加わる。例えば $s=8m_\mu^2$ の直角値は質量を無視した断面積の $(3/2)/\sqrt2$ 倍であり、前方値は $1/\sqrt2$ 倍。質量があると全率が減る場合でも、特定の角度で対応する無質量の値を超えることがある。

立体角ではなく $t$ を測るなら、$t=m_e^2+m_\mu^2-s/2+(s/2)\beta_i\beta_f\cos\theta$。従って方位角を積分した分布は $\dd\sigma/\dd t=(4\pi/(s\beta_i\beta_f))\,\dd\sigma/\dd\Omega$。終状態が閾値へ近づくと、$\dd\sigma/\dd t$ で見えなくなった一つの $\beta_f$ は、許される $t$ 区間の幅 $s\beta_i\beta_f$ に移る。変数を換えた後も全断面積の閾値の冪は同じである。

## ヘリシティ構造と角分布

### カイラルな電流と角度

質量ゼロのベクトル電流は $\bar\psi_L\gamma^\mu\psi_L+\bar\psi_R\gamma^\mu\psi_R$ と分かれる。同じカイラル成分の場が作る粒子・反粒子のヘリシティは逆であるため、同一のヘリシティを持つ入射対は寄与しない。ビーム軸を $z$ に取ると、寄与する対の角運動量射影は $m=\pm1$。終状態の軸に対しても $m'=\pm1$ である。スピン一の回転は横円偏極の内積で計算でき、$y$ 軸周りの回転で $\vb e_\pm^*\cdot R_y(\theta)\vb e_\pm=(1+\cos\theta)/2$、反対の偏極との内積は位相の選択を除き $(1-\cos\theta)/2$。従って平方の和は $(1+\cos\theta)^2+(1-\cos\theta)^2=2(1+\cos^2\theta)$ となり、トレースの角度依存性を再現する。例えば前方では符号の異なる射影の振幅が零である。

### 角度を限った検出器の受容率

両ビームに対して対称な角度領域 $|\cos\theta|<a$ を考える。$x=\cos\theta$ と置くと寄与する比は $\int_{-a}^a(1+x^2)\dd x/\int_{-1}^1(1+x^2)\dd x=(3a+a^3)/4$。例えば $a=1/2$ なら $13/32$ となる。等方分布なら比は $a$ なので、同じ立体角でも検出割合が違う。さらに片側のヘリシティ振幅だけを選べるなら $(1+x)^2$ の奇の部分は前後差を与える。$0<x<1$ の積分は $7/3$、$-1<x<0$ は $1/3$ だから前後非対称度は $3/4$ である。逆のヘリシティを同じ重みで足すと符号が反転して相殺する。

### 回転行列の係数と振幅の規格化

$\vb e_\pm=(\vb e_x\pm i\vb e_y)/\sqrt2$ とする。$R_y(\theta)\vb e_y=\vb e_y$、$R_y(\theta)\vb e_x=\cos\theta\,\vb e_x-\sin\theta\,\vb e_z$ なので、同じ円偏極同士の内積は $(1+\cos\theta)/2$、逆の円偏極では $(\cos\theta-1)/2$ となる。回転の行列要素に現れる半分はこの偏極の規格化に由来する。質量ゼロの四つの非零ヘリシティ振幅の絶対値を $C(1+\cos\theta)$ が二つ、$C(1-\cos\theta)$ が二つと書く。ベクトル相互作用は左右の係数を等しくするため、無偏極平均は $C^2(1+\cos^2\theta)$。前節の二乗が $e^4(1+\cos^2\theta)$ だから $C=e^2$ である。スピン軸の位相の選び方は個々の振幅の位相を変えるが、ここでの絶対値は変えない。

:::example[初期ヘリシティを選んだ測定]
初期の一つの非零ヘリシティ関係を固定し、終状態を和にすると二乗は $e^4[(1+x)^2+(1-x)^2]$。従って無偏極結果の二倍である。初期の四通りのうち二つは零であり、非零な二つの角分布が左右で等しいためである。終状態のヘリシティも一方に固定すると、前後の奇成分が残る。
:::

## Compton 散乱と Klein–Nishina 公式

初期電子 $p=(m,\vb0)$、入射光子 $k$、出射光子 $k'$、出射電子 $p'$ とする。$p'^2=m^2$ から $m(\omega-\omega')=\omega\omega'(1-\cos\theta)$。二つの光子の挿入順序を足すと振幅は
$$
\mathcal M=-e^2\bar u'\left[
\not \epsilon'{}^*\frac{\not p+\not k+m}{2p\cdot k}\not \epsilon
+\not \epsilon\frac{\not p-\not k'+m}{-2p\cdot k'}\not \epsilon'{}^*
\right]u .
\label{eq-qft-compton-amplitude}
$$
偏極和を計算する前に縦成分を調べる。$\epsilon\mapsto k$ とし、$Q(r)=\not r-m$ とおく。第一図では $\not k=Q(p+k)-Q(p)$ と $Q(p)u=0$ により $\bar u'\not\epsilon'{}^*u$ が残る。第二図では $\not k=Q(p')-Q(p-k')$ と $\bar u'Q(p')=0$ により同じ鎖が逆符号で残る。出射側は $k\leftrightarrow-k'$ と二つの頂点の交換で同様に消える。従って物理偏極和と $-\eta$ の差に含まれる光子運動量は、両図を足した振幅との縮約で消える。$a=p\cdot k$、$b=p\cdot k'$ とし、この置換で得たトレースを二つの自己積と交差積に分けると
$$
\begin{aligned}
\overline{|\mathcal M|^2}_{ss}/(2e^4)&=b/a+m^2/a+m^4/a^2,\\
\overline{|\mathcal M|^2}_{uu}/(2e^4)&=a/b-m^2/b+m^4/b^2,\\
\overline{|\mathcal M|^2}_{su+us}/(2e^4)&=m^2(1/a-1/b)-2m^4/(ab)
\end{aligned}
$$
となる。この整理では $\gamma^\mu\not r\gamma_\mu=-2\not r$ と $\gamma^\mu\not r\not v\gamma_\mu=4r\cdot v$ を使い、$p'\cdot k=b$、$p'\cdot k'=a$ に置いた。三項の和を取り、$m^2(1/a-1/b)=\cos\theta-1$ を用いると、平方は $2e^4[a/b+b/a-\sin^2\theta]$ である。実験室系で電子の運動量をデルタで積分すると、残るエネルギーデルタの微分は $(m+\omega(1-\cos\theta))/E'$。これを :eqref[eq-qft-rates] に入れると
$$
\frac{\dd\sigma}{\dd\Omega}=\frac{\alpha^2}{2m^2}
\left(\frac{\omega'}{\omega}\right)^2
\left(\frac{\omega'}{\omega}+\frac{\omega}{\omega'}-\sin^2\theta\right)
\label{eq-qft-klein-nishina}
$$
となる。これが :term[Klein–Nishina 公式]{en="Klein–Nishina formula" yomi="くらいんにしなこうしき"} である。$\omega/m\to0$ では $(1+\cos^2\theta)\alpha^2/(2m^2)$、全断面積 $8\pi\alpha^2/(3m^2)$ の Thomson 極限になる。

### 反跳とエネルギーのJacobian

$r=\omega'/\omega$ とする。空間デルタで $\vb p'=\vb k-\vb k'$ を代入した後、$f(\omega')=m+\omega-\omega'-\sqrt{m^2+\omega^2+\omega'^2-2\omega\omega'\cos\theta}$ である。根では $E'=m+\omega-\omega'$ なので $|f'|=(m+\omega(1-\cos\theta))/E'$。従って残る測度の $\omega'/[E'|f'|]$ は $\omega'/(m+\omega(1-\cos\theta))=\omega r^2/m$ である。これを流束 $4m\omega$ で割ると :eqref[eq-qft-klein-nishina] の $r^2$ が出る。例えば $\omega=m$、$\theta=\pi/2$ では $r=1/2$ であり、微分断面積の Thomson 値に対する比は $r^2(r+r^{-1}-1)=3/8$ である。

### 全 Klein–Nishina 断面積を積分する

$x=\omega/m>0$、$u=1+x(1-\cos\theta)$ とする。これは $\omega/\omega'$ であり、$\theta=0,\pi$ で $1,1+2x$ となる。$\dd\Omega=2\pi\,\dd u/x$ を使うと :eqref[eq-qft-klein-nishina] の角積分は
$$
\sigma=\frac{\pi\alpha^2}{m^2x}\int_1^{1+2x}\dd u\,
\left[\frac1{u^3}+\frac1u-\frac{2(u-1)}{xu^2}+\frac{(u-1)^2}{x^2u^2}\right].
$$
最後の二項は $\sin^2\theta=2(u-1)/x-(u-1)^2/x^2$ を代入したものである。被積分関数を $u^{-3}$、$u^{-2}$、$u^{-1}$、定数の冪に整理すると
$$
u^{-3}+(2/x+1/x^2)u^{-2}
 +(1-2/x-2/x^2)u^{-1}+x^{-2}.
$$
原始関数は $-1/(2u^2)-(2/x+1/x^2)/u+(1-2/x-2/x^2)\log u+u/x^2$。従って $U=1+2x$ とすれば
$$
\sigma=\frac{\pi\alpha^2}{m^2x}
\left[\tfrac12(1-U^{-2})+(2/x+1/x^2)(1-U^{-1})+
(1-2/x-2/x^2)\log U+(U-1)/x^2\right].
$$

:::example[入射光子のエネルギーが電子質量と等しい場合]
$x=1$ では $U=3$ で、角括弧の有理部分は $4/9+2+2=40/9$、対数部分は $-3\log3$。従って Thomson 値に対する全断面積の比は $(3/8)(40/9-3\log3)\simeq0.431$ である。$\theta=\pi/2$ だけで求めた $3/8$ とは異なる。全断面積には前方での反跳の小ささも含まれる。
:::

:::example[低エネルギー展開と高エネルギー展開]
$x\ll1$ なら $r=1-x(1-\cos\theta)+\order{x^2}$、$r+r^{-1}=2+\order{x^2}$。従って微分断面積の一次補正は $-2x(1-\cos\theta)$ 倍であり、奇の $\cos\theta$ 部分は全角積分で消える。全断面積は $\sigma=\sigma_{\rm Th}(1-2x+\order{x^2})$。一方、$x\gg1$ では括弧内の先頭が $\log(2x)+1/2$ なので $\sigma\simeq(\pi\alpha^2/(m^2x))[\log(2x)+1/2]$。高エネルギー値は Thomson 値から大きく離れ、角度を固定して低エネルギー展開する方法は使えない。
:::

### Compton の自己積を直接縮約する

中間運動量を $r=p+k$ とし、$N=\not r+m$ とする。$S=2a$ は内部線の分母である。非物理成分を含む $-\eta$ で二つの偏極を縮約した $s$ 図の自己積は、
$$
\frac{e^4}{4S^2}\tr[(\not p'+m)\gamma^\nu N\gamma^\mu
  (\not p+m)\gamma_\mu N\gamma_\nu]
=\frac{e^4}{4S^2}\tr[(-2\not p'+4m)N(-2\not p+4m)N].
$$
二つの順序を加えた振幅でこの偏極和が許されるため、ここでは自己積をその和の一項として整理している。最後のトレースを $4$ で割ると
$$
4[2(p'\cdot r)(p\cdot r)-(p'\cdot p)r^2]
+4m^2(p'\cdot p)-16m^2(p'\cdot r+p\cdot r)+16m^2(r^2+m^2).
$$
ここで $p'\cdot p=m^2+a-b$、$p'\cdot r=p\cdot r=m^2+a$、$r^2=m^2+2a$。代入すると $8(ab+m^2a+m^4)$ になり、自己積は $2e^4(b/a+m^2/a+m^4/a^2)$。$k\mapsto-k'$、$k'\mapsto-k$ は $a\mapsto-b$、$b\mapsto-a$ に対応し、逆順の自己積も同じトレースで求まる。この計算では $m^4/a^2$ などの項が自己積に出る。低エネルギー極限でこれらを残したまま一図の確率と解釈すると、二図の干渉を落とすことになる。

### Thomson 極限での偏極の選別

電子を初期静止系に置き、物理光子の偏極を $\epsilon^0=\epsilon'{}^0=0$ とする。殻方程式と $p\cdot\epsilon=p\cdot\epsilon'=0$ から、分子の $\not p+m$ は右端の偏極を通すと消える。従って二つの分子にはそれぞれ $\not k$ と $-\not k'$ だけが残る。$\omega\ll m$ で出射スピノルも静止スピノルで近似でき、三つの空間ガンマの積は上下成分をつなぐので寄与しない。時間成分を用いた和は Pauli 行列の反交換子 $\{\boldsymbol\sigma\cdot\boldsymbol\epsilon'{}^*,\boldsymbol\sigma\cdot\boldsymbol\epsilon\}=2\boldsymbol\epsilon'{}^*\cdot\boldsymbol\epsilon$ に比例する。従って全体位相を除く先頭の振幅は $2e^2\boldsymbol\epsilon'{}^*\cdot\boldsymbol\epsilon$ であり、線偏極を固定した微分断面積は $\alpha^2|\boldsymbol\epsilon'{}^*\cdot\boldsymbol\epsilon|^2/m^2$。散乱面に垂直な入射偏極では出射の偏極和が一、散乱面内の入射偏極では $\cos^2\theta$ になる。二つを平均して Thomson 分布を得る。

### 交差した二光子による対生成

Compton の解析式を用いて $\gamma(k_1)\gamma(k_2)\to e^-(p_3)e^+(p_4)$ を求める。Compton の $p,k,k',p'$ をそれぞれ $-p_4,k_1,-k_2,p_3$ に取る。旧 $a=p\cdot k$ は $-A$、旧 $b=p\cdot k'$ は $B$、ここで $A=p_4\cdot k_1>0$、$B=p_4\cdot k_2>0$。電子を陽電子へ交差した外脚では、物理的な和 $\not p_4-m$ は解析式の $\not p+m=-\not p_4+m$ と逆符号になる。光子を交差しても偏極和は同じで、入射の平均はどちらの反応でも $1/4$。従って物理的な振幅二乗は
$$
\overline{|\mathcal M_{\gamma\gamma}|^2}
=2e^4\left[
 \frac AB+\frac BA
 +2m^2\left(\frac1A+\frac1B\right)
 -m^4\left(\frac1A+\frac1B\right)^2
 \right].
$$
電子の外場を交差したことによる負号を省くと、この二乗に全体の負号が付いてしまう。終状態の電子と陽電子は別種であり、位相空間に $1/2!$ は不要である。入射の二光子には異なる運動量を指定したため、入射にもその因子を入れない。

重心で $x=\cos\theta$、$\theta$ は $k_1$ と $p_3$ の角度。$\beta=\sqrt{1-4m^2/s}$ として $A=s(1+\beta x)/4$、$B=s(1-\beta x)/4$。角括弧は
$$
\frac{4(2-\beta^2)}{1-\beta^2x^2}-2
-\frac{4(1-\beta^2)^2}{(1-\beta^2x^2)^2}.
$$
第一項と第三項の角積分には、$L=\log[(1+\beta)/(1-\beta)]$ として
$$
\int_{-1}^1\frac{\dd x}{1-\beta^2x^2}=\frac L\beta,\qquad
\int_{-1}^1\frac{\dd x}{(1-\beta^2x^2)^2}
=\frac1{1-\beta^2}+\frac L{2\beta}
$$
を用いる。二番目は $x/(1-\beta^2x^2)$ の微分と一番目を加えることで得られる。位相空間の $\beta$ と流束を含めた全断面積は
$$
\sigma_{\gamma\gamma}
=\frac{4\pi\alpha^2}{s}
\left[\frac{3-\beta^4}{2}L-\beta(2-\beta^2)\right].
$$
例えば閾値直上では $L=2\beta+O(\beta^3)$ なので $\sigma_{\gamma\gamma}\simeq 4\pi\alpha^2\beta/s$。一方 $s\gg m^2$ では $\sigma_{\gamma\gamma}\simeq4\pi\alpha^2[\log(s/m^2)-1]/s$。高エネルギーの式を得るときに質量を最初から零とすると $x=\pm1$ で角積分が発散するため、質量を残して積分した後に極限を取る必要がある。

:::example[逆反応の断面積]
$\gamma\gamma\leftrightarrow e^+e^-$ の木の反応を同じ $s$ で比べる。行列要素の二乗は時間反転で同じであり、入射の内部状態数も四で等しい。生成では $p_{\rm out}/p_{\rm in}=\beta$、消滅では $1/\beta$。後者に二光子の $1/2!$ が加わるので $\sigma(e^+e^-\to2\gamma)=\sigma_{\gamma\gamma}/(2\beta^2)$。閾値で相対速度 $v_{\rm rel}=2\beta$ を掛けると $\pi\alpha^2/m^2$。一重項だけが消滅に寄与する静止極限では、無偏極平均を四倍して先に用いた $4\pi\alpha^2/m^2$ を回復する。
:::

### 入射の線偏極を指定した Klein–Nishina 分布

入射電子の静止系で $\epsilon^0=0$、$\epsilon^2=-1$、$k\cdot\epsilon=0$ とする。電子のスピンだけを入射で平均し、出射光子を和にする。自由な偏極方向に対する二乗は、入射偏極について二次である。散乱面で反転する成分と変わらない成分の干渉は、無偏極電子の反射対称性で消える。従って面内と面外の二つを求めれば任意の直線偏極が決まる。

面外では $\epsilon\cdot k'=0$ であり、$\epsilon\cdot p'=\epsilon\cdot p=0$ でもある。$P=\not p$、$P'=\not p'$、$R=P+\not k$、$U=P-\not k'$、$E=\not\epsilon$ と書く。出射偏極を $-\eta$ で縮約した $s$ 図の自己積のトレースでは $E(P+m)E=P-m$。中間分子に挟むと $(R+m)(P-m)(R+m)=2a\not k$ になる。ここでは Clifford 関係と $R^2=m^2+2a$ を用いた。さらに $\gamma^\nu\not k\gamma_\nu=-2\not k$、$p'\cdot k=b$ だから、電子の平均と二乗の分母を入れて $s$ 図は $2e^4 b/a$。$u$ 図はトレースを巡回して $E(P'+m)E=P'-m$ を使う。$U=P'-\not k$ なので $(U+m)(P'-m)(U+m)=2b\not k$、$\tr[(P+m)\not k]=4a$ から $2e^4 a/b$ となる。

干渉の一方のトレースは $\tr[(P'+m)\gamma^\nu(R+m)E(P+m)\gamma_\nu(U+m)E]$。ガンマの反交換により中央の縮約は $2P\not k E$。さらに $E(U+m)E=U-m$ なので、残りは $2\tr[(P'+m)P\not k(U-m)]$ となる。四ガンマの式を入れて $8$ をくくると $(p'\cdot p-m^2)(a+b)+b^2-a^2$。$p'\cdot p-m^2=a-b$ なのでこれは零。従って面外の総和は $2e^4(b/a+a/b)$ になる。この分解で共変な偏極和を用いたのは出射側だけである。

二つの入射偏極の無偏極平均で角括弧は $r+r^{-1}-\sin^2\theta$ である。面外の値が $r+r^{-1}$ なので、面内は $r+r^{-1}-2\sin^2\theta$。入射の偏極と散乱面のなす角を $\phi$ とすれば、面内に重み $\cos^2\phi$、面外に $\sin^2\phi$ が付く。従って
$$
\frac{\dd\sigma_{\rm lin}}{\dd\Omega}
=\frac{\alpha^2 r^2}{2m^2}
 \bigl[r+r^{-1}-2\sin^2\theta\cos^2\phi\bigr].
$$
例えば $\theta=\pi/2$ では低エネルギーで面内が零、面外が Thomson の最大値となる。一方 $\omega=m$ なら $r=1/2$ なので、同じ角度でも面内と面外の比は $(1/2+2-2)/(1/2+2)=1/5$。反跳が無視できないと直角方向の面内散乱にも非零の寄与が残る。

## Bhabha 散乱と Møller 散乱

Bhabha 散乱 $e^-e^+\to e^-e^+$ では消滅 $s$ と交換 $t$、Møller 散乱 $e^-e^-\to e^-e^-$ では $t$ と $u$ を足す。同種出射の Møller の鎖にはフェルミオン交換による負号を付ける。高エネルギーで質量を無視すると
$$
\begin{aligned}
\overline{|\mathcal M_{\rm B}|^2}&=2e^4\!\left[\frac{u^2+t^2}{s^2}+\frac{u^2+s^2}{t^2}+\frac{2u^2}{st}\right],\\
\overline{|\mathcal M_{\rm M}|^2}&=2e^4\!\left[\frac{s^2+u^2}{t^2}+\frac{s^2+t^2}{u^2}+\frac{2s^2}{tu}\right].
\end{aligned}
$$
対角項は対生成の二トレースを脚の交差で置き換えたものである。Møller の交差では一つのトレース
$\tr(\not p_3\gamma^\mu\not p_1\gamma^\nu\not p_4\gamma_\mu\not p_2\gamma_\nu)$
を使う。$\gamma^\mu\gamma^\alpha\gamma^\nu\gamma^\beta\gamma_\mu=-2\gamma^\beta\gamma^\nu\gamma^\alpha$ と二行列の縮約により $-32(p_1\cdot p_2)(p_3\cdot p_4)=-8s^2$。入射平均 $1/4$ と相対符号 $-2$ で交差項が $4e^4s^2/(tu)$ になる。Bhabha は $s,u$ の交差で従う。例えば $t=u=-s/2$ の Møller では $36e^4$ である。終状態を全立体角で積分するなら同種の Møller には $1/2$ を用いる。

### 角度の二重計数を除く

Møller 散乱では重心系で一方の電子を角度 $\theta$ に置くと、他方は $\pi-\theta$ にある。$x=\cos\theta$ に対する振幅二乗は $t=-s(1-x)/2$ と $u=-s(1+x)/2$ の交換で不変である。従って対称な角度領域で全事象を数えるには $-a<x<a$ の積分に $1/2$ を掛けるか、$0<x<a$ だけを積分すればよい。例えば $t=u=-s/2$ の Bhabha 振幅二乗では角括弧の三項が $1/2,5,-1$ となり、結果は $9e^4$。終状態の電子と陽電子にはこの半分の係数を掛けない。

### 同種粒子散乱で角度カットを実際に積分する

質量を無視した Møller 散乱で $x=\cos\theta$ とする。$t=-s(1-x)/2$ と $u=-s(1+x)/2$ を代入すると角括弧は
$$
\frac{4+(1+x)^2}{(1-x)^2}
+\frac{4+(1-x)^2}{(1+x)^2}
+\frac8{1-x^2}=\frac{2(3+x^2)^2}{(1-x^2)^2}.
$$
従って二重計数を除くため正の $x$ だけを選ぶと、$0<x<a<1$ の断面積は
$$
\sigma_{\rm M}(|x|<a)=\frac{2\pi\alpha^2}{s}
  \int_0^a\dd x\,\frac{(3+x^2)^2}{(1-x^2)^2}.
$$
積分の分子を $(1-x^2)^2+8(1+x^2)$ と分けると、被積分関数は $1+8\,\dd[x/(1-x^2)]/\dd x$。従って積分は $a+8a/(1-a^2)$ である。例えば $a=1/2$ では $35/6$。$a\uparrow1$ でこれが発散するのは前方およびそれと同じ事象を表す後方の長距離散乱による。電子質量を戻しても質量ゼロの媒介光子の前方極は残る。

### 電子質量を保った二つのチャネル

Møller の $t$ 図のスピン和を先に計算する。$p_i^2=m^2$ とし、$L^{\mu\nu}(p_a,p_b)=\tr[(\not p_a+m)\gamma^\mu(\not p_b+m)\gamma^\nu]/4$。$p_a,p_b$ が $t$ で結ばれれば $L^{\mu\nu}=p_a^\mu p_b^\nu+p_a^\nu p_b^\mu+t\eta^{\mu\nu}/2$。縮約の四つの運動量項は $2[(p_3\cdot p_4)(p_1\cdot p_2)+(p_3\cdot p_2)(p_1\cdot p_4)]$。一つまたは二つの計量の項を足すと $t(p_3\cdot p_1+p_4\cdot p_2)+t^2=2m^2t$。よって $t$ 図の二乗は
$$
\frac{2e^4}{t^2}\left[(s-2m^2)^2+(u-2m^2)^2+4m^2t\right].
$$
$u$ 図の自己積は $t,u$ の交換で従う。

干渉項では
$$
X=\tr[(\not p_3+m)\gamma^\mu(\not p_1+m)\gamma^\nu
       (\not p_4+m)\gamma_\mu(\not p_2+m)\gamma_\nu]
$$
を計算する。四つの $\not p$ の項は既に求めたように $-32(p_1\cdot p_2)(p_3\cdot p_4)$。全て $m$ の項は $\gamma^\mu\gamma^\nu\gamma_\mu=-2\gamma^\nu$ から $-32m^4$。残る六組のうち $\not p_3,\not p_1$ を残す項では、先に $\gamma^\nu\gamma_\mu\gamma_\nu=-2\gamma_\mu$ を使い、次に $\gamma^\mu\not p_1\gamma_\mu=-2\not p_1$。従って $16m^2\,p_3\cdot p_1$ になる。非隣接の $\not p_3,\not p_4$ を残す項は $\gamma^\mu\gamma^\nu\not p_4\gamma_\mu=4p_4^\nu$ から $16m^2\,p_3\cdot p_4$。他の項はトレースの巡回と $\mu,\nu$ の交換でこの二形へ移る。従って二つの $\not p$ の和は $16m^2\sum_{i<j}p_i\cdot p_j=32m^2(s-m^2)$ である。全てを足して $X=-8(s-2m^2)(s-6m^2)$。入射平均の $1/4$ と交換の $-2$ を掛ければ、干渉項は $4e^4(s-2m^2)(s-6m^2)/(tu)$ となる。

Bhabha は一つの入射と一つの出射を反粒子へ交差し、$s,u$ を入れ替える。交差した二つのスピン和にはそれぞれ負号が出るので、二乗では符号は打ち消し合う。例えば Bhabha の干渉項は $4e^4(u-2m^2)(u-6m^2)/(ts)$。物理領域で $s>0,t<0$ のため、同じ式の形でも Møller と干渉の符号が必ず同じとは限らない。

:::example[非相対論的な同種粒子で干渉の符号を見る]
$|\vb p|\ll m$ では $s=4m^2+O(\vb p^2)$、$t,u=O(\vb p^2)$。そのとき $t$ 図と $u$ 図の分子はそれぞれ $16e^4m^4$、干渉項の分子は $-16e^4m^4$ である。従って二乗の先頭は $16e^4m^4(t^{-2}+u^{-2}-(tu)^{-1})$。これは低速の直接・交換のスピン作用素 $t^{-1}1-u^{-1}P$ を二乗して四スピンで平均した結果と一致する。実際 $P^2=1$、$\tr P=2$ なので交差項に $\tr P/4=1/2$ が付く。$\pi/2$ で三重項が相殺しても、平均には一重項の散乱が残る。
:::

## 束縛状態とポジトロニウム

電子・陽電子を低速で扱うと Coulomb ポテンシャル $-\alpha/r$ と換算質量 $\mu=m/2$ を持つ Schrödinger 問題になる。この束縛状態を :term[ポジトロニウム]{en="positronium" yomi="ぽじとろにうむ"} という。基底波動関数 $N e^{-r/a}$ を放射状方程式へ入れると $1/r$ の係数の消失で $a=1/(\mu\alpha)$、定数から $E_1=-\mu\alpha^2/2=-m\alpha^2/4$。規格化は $\int\dd^3x\,e^{-2r/a}=\pi a^3$ なので $|\psi(0)|^2=1/(\pi a^3)$ である。電子の質量を $511\,\mathrm{keV}$、$\alpha\simeq1/137$ とすれば基底の束縛は約 $6.8\,\mathrm{eV}$。

粒子と反粒子の交換で空間波動関数に $(-1)^L$、スピン交換に $(-1)^{S+1}$ が付き、生成演算子の入れ替えでも負号が付くので $C=(-1)^{L+S}$ である。$S$ 波の一重項は $C=+$ で二光子が最低次数、三重項は $C=-$ で三光子が最低次数である。単一光子では静止した親の時間的な四元運動量を保存できない。消滅の距離が $1/m$ 程度で $a$ より短いので、先頭の $S$ 波の幅は短距離の $\sigma v_{\rm rel}$ に $|\psi(0)|^2$ を掛ける形になる。

### 一重項の二光子幅

閾値で相対速度を零へ近づけ、入射 $e^-$ と $e^+$ を $p=(m,\vb0)$ に置く。出射を $k_{1,2}=m(1,0,0,\pm1)$、二偏極を $x,y$ 方向に選ぶ。二つの挿入の電子分母はいずれも $-2m^2$ で、Dirac 基底の静止スピノルは $u=\sqrt{2m}(\xi,0)^{\mathsf T}$、$v=\sqrt{2m}(0,\eta)^{\mathsf T}$。外側の空間ガンマは非対角であるから、中間分子の $m$ の項は上下成分をつなげず零である。$\pm m\gamma^3$ の項を乗じて二順序を足すと、直線偏極が平行なら消え、直交なら位相を除き $2e^2\eta^\dagger\xi$。二つの直交する偏極配置と二つずつの入射スピンを足して四で割ると $\overline{|\mathcal M|^2}=4e^4$ となる。

電荷共役により閾値の $S$ 波で二光子へ行くのは一重項のみなので、一重項の二乗は無偏極平均の四倍である。入射流束を $v_{\rm rel}$ で割ると $4m^2$ に近づく。終状態の同一性の $1/2$ と二体位相空間 $1/(8\pi)$ を掛けると、$\sigma_{\rm singlet}v_{\rm rel}=4\pi\alpha^2/m^2$。従って先頭幅は
$$
\Gamma_{\rm para}
 =|\psi(0)|^2\frac{4\pi\alpha^2}{m^2}
 =\frac{m\alpha^5}{2}.
$$
$m\simeq511\,\mathrm{keV}$、$\alpha\simeq1/137$ では幅が約 $5.3\times10^{-6}\,\mathrm{eV}$、$\hbar\simeq6.58\times10^{-16}\,\mathrm{eV\,s}$ を戻した寿命は約 $1.25\times10^{-10}\,\mathrm{s}$ である。これは束縛時の速度補正と輻射補正を落とした結果である。

### 接触相互作用による超微細分裂の一部

低速 QED の一光子交換には二つの磁気モーメントから来るスピン接触項 $V_{\rm ss}=8\pi\alpha\,\vb S_1\cdot\vb S_2\,\delta^3(\vb r)/(3m^2)$ がある。相対運動量に対する磁気電流は $\vb q\times\vb S_i/m$ に比例する。二電流の縮約から $\vb S_1\cdot\vb S_2-(\vb q\cdot\vb S_1)(\vb q\cdot\vb S_2)/|\vb q|^2$ が出て、$S$ 波では角平均で $2\vb S_1\cdot\vb S_2/3$ となる。この係数に $4\pi\alpha/m^2$ を掛けたものが上の接触項である。$\vb S_1\cdot\vb S_2=[S(S+1)-3/2]/2$ なので三重項と一重項の差は $1$。基底の $|\psi(0)|^2$ を使うと交換による差は $m\alpha^4/3$ である。電子と陽電子では仮想一光子への消滅も短距離に寄与するため、この交換の値だけでは全分裂は得られない。

### 一光子消滅の選択則と先頭超微細分裂

重心の閾値 $p_1=p_2=(m,\vb0)$ で時間成分 $\bar v\gamma^0u$ は零、空間成分は $2m\,\eta^\dagger\sigma^i\xi$ である。ここで $u=\sqrt{2m}(\xi,0)^\mathsf T$、$v=\sqrt{2m}(0,\eta)^\mathsf T$。二つの低速スピンを、係数行列 $M=\xi\eta^\dagger$ の空間で表す。回転は $M\mapsto U M U^\dagger$。従って一重項は $1/\sqrt2$、三重項の直交基底は $\sigma^a/\sqrt2$ である。行列要素 $\eta^\dagger\sigma^i\xi=\tr(\sigma^i M)$ を評価すると、一重項は零、三重項では $\sqrt2\,\delta^{ia}$。これで空間ベクトル電流が三重項だけを消滅させる。

仮想光子の運動量は $P=(2m,\vb0)$ で、分母は $4m^2$。両端の電流の係数 $2m$ を掛けると相殺し、空間の計量で $-\delta_{ij}$ が残る。従って三重項の閾値振幅は $-2e^2$、一重項は零である。外脚を低速に規格化して $\widetilde V=-\mathcal M/(4m^2)$ と対応させると、三重項への接触項は $2\pi\alpha\,\delta^3(\vb r)/m^2$。一重項との差は $2\pi\alpha|\psi(0)|^2/m^2=m\alpha^4/4$ である。先の交換項と合わせて
$$
\Delta E_{\rm hfs}^{(0)}=m\alpha^4\left(\frac13+\frac14\right)=\frac7{12}m\alpha^4 .
$$
これは実光子への崩壊幅とは異なり、仮想一光子の実の局所効果を求めたもの。速度展開の高次と輻射補正は含めていない。

### 離散対称性で崩壊の候補を絞る

束縛状態の空間反転では軌道の球面調和関数が $(-1)^L$ を与え、Dirac 粒子と反粒子の固有パリティの積は $-1$ である。従って $P=(-1)^{L+1}$。電荷共役の結果 $C=(-1)^{L+S}$ と合わせると、基底一重項は $J^{PC}=0^{-+}$、三重項は $1^{--}$ になる。電磁場が電荷共役で反転するため $n$ 個の光子の状態は $C=(-1)^n$。これで一重項の先頭消滅は二光子、三重項の先頭は三光子と分かる。三光子幅の次数は、閾値で $\sigma v_{\rm rel}$ の次元が $m^{-2}$、三つの QED 頂点を二乗して $\alpha^3$、原点の確率密度が $m^3\alpha^3$ なので $m\alpha^6$。厳密な係数には三体位相空間と六つの挿入順序の干渉を含む積分が要る。

:::example[超微細分裂と幅の尺度を比較する]
上で求めた二光子幅と基底の超微細分裂の比は $\Gamma_{\rm para}/\Delta E_{\rm hfs}^{(0)}=6\alpha/7$。従って $\alpha\ll1$ では幅の小ささを分裂の大きさから区別できる。これは同じ一光子相互作用から出る桁の見積もりでも、束縛状態における局所消滅率とエネルギー補正では $\alpha$ の冪が異なる例である。
:::

### 三重項の三光子幅

基底三重項でも硬い消滅の部分は静止した電子と陽電子の散乱振幅で計算できる。ここでは電子と陽電子を $p=(m,\vb0)$ に置き、束縛エネルギーを位相空間に入れる補正を落とす。光子 $i$ のエネルギーを $\omega_i>0$、方向を $\widehat{\vb k}_i$、時間成分が零の実横偏極を $\vb e_i$ とする。$E_i=\boldsymbol\sigma\cdot\vb e_i$、$N_i=\boldsymbol\sigma\cdot\widehat{\vb k}_i$、$D_i=\boldsymbol\sigma\cdot(\widehat{\vb k}_i\times\vb e_i)$ と置く。すると $N_i E_i=iD_i$、$E_i N_i=-iD_i$ である。これは Pauli 行列の積と偏極の横性から従う。

電子端から順に $i,j,k$ の光子を出す場合、中間運動量は $p-k_i$ と $k_k-p$。分母は $-2m\omega_i$ と $-2m\omega_k$ である。Dirac 基底では $\not\epsilon_i$ の非対角ブロックは上が $-E_i$、下が $E_i$。右端から $\not\epsilon_i u$ を一つ目の分子に通すと、$(N_iE_i\xi,E_i\xi)^\mathsf T$ に $\sqrt{2m}\,\omega_i$ が掛かる。次に $\not\epsilon_j$ を通すと上の成分は $-E_jE_i\xi$、下は $E_jN_iE_i\xi$ になる。二つ目の分子の上の行のブロックは $(\omega_k,-\omega_kN_k)$。さらに左端の $\not\epsilon_k$ と $\bar v$ を掛けると、$\eta^\dagger E_k(E_jE_i+N_kE_jN_iE_i)\xi$ に $2m\omega_i\omega_k$ が掛かる。従ってこの順序の振幅は、共通位相を除いて
$$
\mathcal M_{ijk}=\frac{e^3}{2m}\eta^\dagger
   (E_kE_jE_i+D_kE_jD_i)\xi .
$$
分母の小さなエネルギー因子は、この静止近似の分子で相殺する。順序 $(i,j,k)$ と逆順 $(k,j,i)$ を足すと、三つの Pauli 行列の積のスカラー項が反転して消える。実際、任意の実三ベクトルに対し
$$
(\boldsymbol\sigma\cdot a)(\boldsymbol\sigma\cdot b)(\boldsymbol\sigma\cdot c)
+(\boldsymbol\sigma\cdot c)(\boldsymbol\sigma\cdot b)(\boldsymbol\sigma\cdot a)
=2\boldsymbol\sigma\cdot[(a\cdot b)c+(b\cdot c)a-(a\cdot c)b].
$$
これは最初の二つを $(a\cdot b)1+i\boldsymbol\sigma\cdot(a\times b)$ と置いて得る。よって三つの逆順の組を加えた振幅は $(e^3/m)\eta^\dagger(\boldsymbol\sigma\cdot\vb V)\xi$ であり、$\vb d_i=\widehat{\vb k}_i\times\vb e_i$ と書けば
$$
\vb V=\sum_{\text{巡回 }i,j,k}\left[
 (\vb e_j\cdot\vb e_k)\vb e_i
 +(\vb d_i\cdot\vb e_j)\vb d_k
 +(\vb d_k\cdot\vb e_j)\vb d_i
 -(\vb d_i\cdot\vb d_k)\vb e_j\right].
$$
三重項のスピンを平均する。スピン行列の基底は $\sigma^a/\sqrt2$ なので、$\eta^\dagger(\boldsymbol\sigma\cdot\vb V)\xi$ の二乗平均は $\sum_a|\tr(\boldsymbol\sigma\cdot\vb V\,\sigma^a/\sqrt2)|^2/3=2|\vb V|^2/3$。一重項ではトレースが零になり、前の電荷共役の選択則と一致する。

三光子の空間運動量は和が零なので、一つの平面にある。その平面の法線を $\vb n$、$\vb t_i=\vb n\times\widehat{\vb k}_i$ とし、偏極の二択を $\vb n,\vb t_i$ に取る。$\vb e_i=\vb n$ なら $\vb d_i=-\vb t_i$、$\vb e_i=\vb t_i$ なら $\vb d_i=\vb n$。$d_i=1-\widehat{\vb k}_j\cdot\widehat{\vb k}_k$ と書く。三つが法線偏極なら $\vb V=(d_1+d_2+d_3)\vb n$。一つだけ法線でそれが $i$ なら $\vb V=(-d_i+d_j+d_k)\vb n$ になる。全て面内なら $\vb V=-\sum_i d_i\vb t_i$、一つだけ面内でそれが $i$ なら $\vb V=d_i\vb t_i-d_j\vb t_j-d_k\vb t_k$。これらは上の $\vb V$ の内積で $\vb t_j\cdot\vb t_k=\widehat{\vb k}_j\cdot\widehat{\vb k}_k$ を使えば得られる。法線の四通りを二乗して足すと相互の積が消え $4\sum_i d_i^2$。面内の四通りでも内積の前の符号が均等になるので同じ結果である。従って光子偏極を和にして初期三重項を平均した相対論的な振幅二乗は
$$
\overline{|\mathcal M|^2}_{S=1}
=\frac{16e^6}{3m^2}\sum_{i=1}^3 d_i^2.
$$

次に $x_i=\omega_i/m$ とする。$P=(2m,\vb0)$ として $(k_j+k_k)^2=(P-k_i)^2=4m^2(1-x_i)$ だから $d_i=2(1-x_i)/(x_j x_k)$。各 $x_i$ は $0<x_i<1$ で $\sum_i x_i=2$。従って $x_1,x_2$ の積分領域は $0<x_1<1$、$1-x_1<x_2<1$ である。三体位相空間では運動量のデルタが第三の運動量を決める。第一の立体角は $4\pi$、第二の方位角は $2\pi$。残る極角のエネルギーデルタの微分は $\omega_1\omega_2/\omega_3$ なので、測度に残る同じ因子が消える。従ってラベル付きの測度は $\dd\Phi_3=\dd\omega_1\dd\omega_2/(32\pi^3)=m^2\dd x_1\dd x_2/(32\pi^3)$。これに $1/3!$、閾値の規格化 $|\psi(0)|^2/(4m^2)$、振幅二乗を掛けると
$$
\frac{\dd^2\Gamma_{\rm ortho}}{\dd x_1\dd x_2}
=\frac{2m\alpha^6}{9\pi}
 \sum_{\text{巡回 }i,j,k}\frac{(1-x_i)^2}{x_j^2x_k^2}.
$$
運動量を全てラベル付けした領域を積分するため、ここでは $1/3!$ が必要である。実験でどの光子か特定しない単光子分布を作るには、この一光子のエネルギー密度の三倍を数える。

全幅のため、一つの循環項に $z=1-x_i$ を適用する。他方の $x_j$ を $y$ と置くと $z<y<1$、$x_k=1+z-y$。部分分数分解から、
$$
\begin{aligned}
\int_z^1 \frac{z^2\,\dd y}{y^2(1+z-y)^2}
&=\frac{2z(1-z)}{(1+z)^2}-\frac{4z^2\log z}{(1+z)^3},\\
\int_0^1 \frac{2z(1-z)}{(1+z)^2}\dd z&=-4+6\log 2.
\end{aligned}
$$
対数項では $z^2/(1+z)^3=(1+z)^{-1}-2(1+z)^{-2}+(1+z)^{-3}$。第一の対数積分 $\int_0^1\log z\,\dd z/(1+z)$ は $-\pi^2/12$ である。幾何級数を積分すると $-\sum_{n\ge1}(-1)^{n-1}/n^2$ となり、偶数項と奇数項を分けた値を使った。残り二つは、原点で零になる原始関数を選んで部分積分すると $-\log2$、$-(\log2+1/2)/2$。従って一項の積分は $\pi^2/3-3$。三項を足すと
$$
\Gamma_{\rm ortho}^{(0)}=\frac{2(\pi^2-9)}{9\pi}m\alpha^6
$$
である。例えば $m=511\,\mathrm{keV}$、$\alpha=1/137$ で $\hbar$ を戻した寿命は約 $1.4\times10^{-7}\,\mathrm{s}$。二光子の幅よりさらに $\alpha$ の一冪が小さい。ここでのスペクトルは短距離の静止近似であり、一光子のエネルギーが束縛の尺度 $m\alpha^2$ 以下になると、途中の束縛準位を保った計算を要する。硬い領域が全幅の先頭を与えるため、全幅の次数と軟端の点ごとの精度は区別する。

## スカラー QED

複素スカラーを $D_\mu=\partial_\mu+iqA_\mu$ で結合すると
$$
(D_\mu\varphi)^*D^\mu\varphi
=\partial_\mu\varphi^*\partial^\mu\varphi-qA_\mu j^\mu+q^2A_\mu A^\mu\varphi^*\varphi
$$
である。$j^\mu$ は :eqref[eq-qft-uone-current]。粒子線に沿って運動量 $p$ が $p'$ へ変わる三点頂点は $-iq(p+p')^\mu$ であり、二光子の局所頂点は $2iq^2\eta^{\mu\nu}$。後者の $2$ は二つの光子を $A^2$ へ結ぶ選択数である。三点の運動量縮約は $k_\mu(p+p')^\mu=p'^2-p^2$ なので、質量殻上で消え、内部線なら逆伝播関数の差になる。スカラーの Compton 散乱ではこの差が作る局所項を二光子頂点が打ち消す。

### スカラー Compton 散乱の接触項を検算する

電子を複素スカラーへ替え、二つの逐次挿入で入射光子の偏極を $k$ に置き換える。第一図では $k\cdot(2p+k)=(p+k)^2-m^2$ なので内部線が消え、出射光子側の頂点の $p'+p+k$ が残る。第二図では $k\cdot(2p-2k'+k)=m^2-(p-k')^2$ なので、一線が消えて逆符号の $p'+p-k$ が残る。逐次挿入の共通因子は $(-iq)^2i=-iq^2$ だから、二つの差は $-2iq^2k^\nu$。二光子頂点 $2iq^2\eta^{\mu\nu}$ を同じ運動量で縮約すると $+2iq^2k^\nu$ になり和が消える。従って局所二光子項は、正規化のためでなくゲージ変換に対する振幅の相殺に必要である。

### 一光子によるスカラー対の生成

質量ゼロの入射電子・陽電子に対し、終状態のスカラーと反スカラーの質量を $m$ とする。振幅は全体符号を除いて $q_e q_\varphi\,\bar v_2\gamma^\mu u_1(p_3-p_4)_\mu/s$。重心系で $p_3-p_4=(0,2\vb p_f)$ であり、入射スピン平均の対称テンソル $p_2^\mu p_1^\nu+p_2^\nu p_1^\mu-\eta^{\mu\nu}s/2$ の非零空間成分は $xx,yy$ の $s/2$ だけである。従って振幅二乗は $q_e^2 q_\varphi^2\beta^2\sin^2\theta/2$、$|q_e|=|q_\varphi|=e$ なら $\dd\sigma/\dd\Omega=\alpha^2\beta^3\sin^2\theta/(8s)$ となる。$\sin^2\theta$ の角積分は $8\pi/3$ なので全断面積は $\pi\alpha^2\beta^3/(3s)$。閾値での $\beta^3$ は、電流が運動量差に比例することと位相空間の $\beta$ の積である。

### 荷電スカラーの保存電流には場の項がある

スカラーの密度を $|D\varphi|^2-m^2|\varphi|^2$ とすると、$-\partial\mathcal L/\partial A_\mu$ で定義する電流は
$$
j_{\rm em}^\mu=q\,i[\varphi^*D^\mu\varphi-(D^\mu\varphi)^*\varphi]
=q\,j^\mu-2q^2 A^\mu|\varphi|^2 .
$$
$D_\mu D^\mu\varphi+m^2\varphi=0$ と共役を用いると、$\partial_\mu j_{\rm em}^\mu=qi[\varphi^*D^2\varphi-(D^2\varphi)^*\varphi]=0$。積の微分で、左の普通微分に入っていない接続の項は中性な積の中で相殺する。この電流を $A_\mu$ でさらに微分すれば二光子頂点に対応する局所応答が出る。例えば $\varphi$ が局所的に一定であっても $A_\mu\ne0$ なら、共変微分と電流は一般に零でない。

### 二光子から荷電スカラー対を作る

複素スカラーの Compton 振幅を交差して $\gamma(k_1)\gamma(k_2)\to\varphi(p_3)\varphi^*(p_4)$ を考える。二つの逐次挿入に加えて、同じ結合の二光子接触項を含める。重心枠で $k_1,k_2$ は $z$ 方向と逆方向、$p_3$ は $xz$ 面内とする。各光子の偏極基底を $x,y$ 方向に取れば、逐次の各頂点の偏極との縮約は $2p_3\cdot\epsilon_i$ または $-2p_4\cdot\epsilon_i$。$\vb p_4=-\vb p_3$ なので二つの値は同じである。内部線の分母は $-2p_4\cdot k_1$ と $-2p_4\cdot k_2$。従って異なる偏極では両方の逐次挿入と接触項が零、$yy$ では接触項の $-2q^2$ だけが残る。

$x=\cos\theta$、$\beta=\sqrt{1-4m^2/s}$ とする。$xx$ の逐次分子は二図とも $4|\vb p_3|^2(1-x^2)=s\beta^2(1-x^2)$。分母の逆数の和は $-4/[s(1-\beta^2x^2)]$ で、二つの三点頂点から来る負号を掛けると逐次部分は正。従って
$$
\mathcal M_{xx}=2q^2\left[-1+\frac{2\beta^2(1-x^2)}{1-\beta^2x^2}\right],
\quad
\mathcal M_{yy}=-2q^2,\quad
\mathcal M_{xy}=\mathcal M_{yx}=0 .
$$
入射平均で $q^4[1+(1-2T)^2]$、$T=\beta^2(1-x^2)/(1-\beta^2x^2)$。例えば閾値では交換図の分子が零となり、接触項だけで $\overline{|\mathcal M|^2}=2q^4$。$|q|=e$ なら $\dd\sigma/\dd\Omega=\alpha^2\beta/(2s)+O(\beta^3)$。高エネルギーでは端点以外で $T\to1$ であり、$0\le T\le1$ から収束の優関数がある。角積分して $\sigma\to2\pi\alpha^2/s$ の先頭を得る。質量を落とすと端点で対数が出た Dirac 対の生成と比べると、スカラーの場合は交換頂点の横運動量が端点の特異性を抑える。

## Ward 恒等式

### 縦偏極を入れて相殺を確かめる

:term[Ward 恒等式]{en="Ward identity" yomi="わーどこうとうしき"} のオンシェル形は、外部光子偏極をその運動量に置き換えた振幅が零になるというものである。:eqref[eq-qft-compton-amplitude] で $\epsilon\to k$ とする。$S^{-1}(r)=\not r-m$ と置き、$\not k=S^{-1}(p+k)-S^{-1}(p)$ を第一項へ入れると $S^{-1}(p)u=0$ により $+\bar u'\not \epsilon'{}^*u$ が残る。第二項では $\not k=S^{-1}(p')-S^{-1}(p-k')$ で、$\bar u'S^{-1}(p')=0$ により同じ鎖が負号つきで残る。よって和は零である。内部への挿入でも隣接逆伝播関数の差になるため、線に沿って全挿入点を足すと端だけが残る。この端の寄与を含む相関の式は [くりこみ](/physics/quantum-field-theory/renormalization/) で用いる。

### 偏極和の基準ベクトルを取り替える

光子の振幅を $\epsilon_\mu\mathcal M^\mu$ とする。偏極和は $-\eta^{\mu\nu}$ に $k^\mu$ または $k^\nu$ を含む項を足したものである。$k_\mu\mathcal M^\mu=0$ を使えば、それらの項は $\mathcal M_\mu\mathcal M_\nu^*$ との縮約で零になる。例えば $k=(\omega,0,0,\omega)$ なら $\mathcal M^0=\mathcal M^3$ なので $-\mathcal M_\mu\mathcal M^{*\mu}=|\mathcal M^1|^2+|\mathcal M^2|^2$。$-\eta$ を使えるのは全振幅で縦成分が相殺するためであり、Compton の一方の挿入順序だけの自己積にこの議論を適用することはできない。

### 任意の位置への光子挿入を足す

木の一つの開いた電子線に、中性の外部挿入を固定した順序で並べる。さらに光子運動量 $k$ を持つ頂点を、その線上の許される全ての位置へ加える。$i$ を除いた内部線を $G(r)=(\not r-m)^{-1}$ とすると、挿入を運動量で縮約した部分は
$$
G(r+k)\not k G(r)=G(r)-G(r+k)
$$
である。左の逆核との差を分けた二項を隣の挿入位置の同じ表式と比べると、中間の位置で互いに逆符号となる。従って和は両端にのみ残る。両端の場が同じ電荷を運ぶ場合、殻方程式を用いた外部の挿入は直前の内部項を打ち消す。外部脚の代わりに源を固定した相関では、端が物理的な状態に射影されないので差が残る。この差の相関関数での形は [くりこみ](/physics/quantum-field-theory/renormalization/) で導く。

同じ運動量恒等式を閉線の全位置へ適用すると、端の代わりにループ運動量を $k$ だけずらした差が残る。収束積分または並進を保つ正則化なら、この差の積分は零である。逆に発散積分をそのままシフトすると境界項の扱いが未定義になるため、運動量恒等式と積分路の規約を両方確認する必要がある。例えば二点の閉ループの横性は、この縮約の最も簡単な応用である。

### 偏極の基準ベクトルを明示する

光的 $k$ と時間的 $n$ に対して $k\cdot n\ne0$ とする。基準 $n$ の静止系で横偏極を選び、$\epsilon\cdot k=\epsilon\cdot n=0$ を課す。二本の偏極に対する和は
$$
\sum_{r=1}^2\epsilon_r^\mu\epsilon_r^{*\nu}
=-\eta^{\mu\nu}
+\frac{k^\mu n^\nu+n^\mu k^\nu}{k\cdot n}
-\frac{n^2 k^\mu k^\nu}{(k\cdot n)^2}.
$$
右辺は $k$ と $n$ との縮約が零であり、両者に直交する二次元部分で $-\eta^{\mu\nu}$ と一致する。従って横偏極の和と同じ作用をする。例えば $n=(1,0,0,0)$、$k=(\omega,0,0,\omega)$ なら空間の $xx,yy$ 成分だけが一で、他は零。$n$ が異なる観測者の速度を表しても、両端に保存電流を縮約すると $k$ を含む項が零となるため、偏極の和の結果は変わらない。個々の偏極の方向を指定した測定では、基準と偏極ベクトルを変換した上で比較する。

### 外部の縦偏極とオフシェルの線を区別する

木レベルで $G^{-1}(p)=\not p-m$ を用いると $k_\mu\gamma^\mu=G^{-1}(p+k)-G^{-1}(p)$。この式を一般の二つの行列の源の間に挟んでも零ではない。例えば四元運動量が両方とも粒子の殻から外れており、$k=(k^0,\vb0)$ なら差は $k^0\gamma^0$ である。外部の物理的な電子を表す $u,\bar u$ に当たると、それぞれの逆核が消えるため縦偏極の寄与が消える。この区別は内部光子の縦成分がどの一図でも個別に消えるという推論を排除する。例えば Compton 散乱の二図では、逆核が中間線に当たって局所的な鎖を残すので、両方の順序を含める必要がある。

ゲージ変換を $A\mapsto A+\partial\omega$ とする外場の振幅では、摂動の一階は $\int j^\mu\partial_\mu\omega$。境界項がなければ $-\int\omega\,\partial_\mu j^\mu$ になり、保存する外部粒子電流では零である。相関の時間順序積にこの部分積分を行うと、時間順序の階段関数も微分される。従ってオンシェルで消える式を相関の任意のオフシェルの値へ延長する場合には両端の接触項を保つ。相関の Ward–Takahashi 恒等式は [くりこみ](/physics/quantum-field-theory/renormalization/) でこの接触項から導く。

## 演習問題

::::exercise{#exr-qft-qed-1 level=basic}
$\theta=\pi/2$ の対生成角分布を前方値と比較せよ。

:::solution
:eqref[eq-qft-muon-cross-section] の角度因子は一、前方では二なので比は $1/2$ である。
:::
::::

::::exercise{#exr-qft-qed-2 level=basic}
実験室で $\omega=m$、$\theta=\pi/2$ の Compton 散乱の $\omega'/\omega$ を求めよ。

:::solution
$1/\omega'-1/\omega=(1-\cos\theta)/m$ なので $\omega'/\omega=1/2$ である。
:::
::::

::::exercise{#exr-qft-qed-3 level=basic}
スカラー QED の三点頂点を光子運動量で縮約せよ。

:::solution
$k=p'-p$ とすれば積は $-iq(p'^2-p^2)$ となる。二つの脚が同じ質量殻なら零である。
:::
::::

::::exercise{#exr-qft-qed-4 level=basic}
$L=0$、$S=1$ のポジトロニウムの $C$ を求めよ。

:::solution
空間交換は正、スピンの交換も正、生成演算子の交換は負なので $C=-1$ である。
:::
::::

::::exercise{#exr-qft-qed-5 level=advanced}
無偏極質量ゼロの対生成で一方の半球への断面積を求めよ。

:::solution
$\cos\theta=c$ とし、$2\pi\alpha^2/(4s)\int_0^1(1+c^2)\dd c$ は $2\pi\alpha^2/(3s)$ である。
:::
::::

::::exercise{#exr-qft-qed-6 level=advanced}
低エネルギー Compton 散乱の全断面積を求めよ。

:::solution
$\omega'/\omega\to1$ なので角度因子は $1+\cos^2\theta$。その立体角積分は $16\pi/3$、前の係数は $\alpha^2/(2m^2)$ なので $8\pi\alpha^2/(3m^2)$ になる。
:::
::::

::::exercise{#exr-qft-qed-7 level=advanced}
ポジトロニウムの基底波動関数が与える $|\psi(0)|^2$ を $m,\alpha$ で表せ。

:::solution
$\dd^3x$ の球積分で $4\pi\int r^2 e^{-2r/a}\dd r=\pi a^3$ なので $|\psi(0)|^2=1/(\pi a^3)$。$a=2/(m\alpha)$ を代入し $m^3\alpha^3/(8\pi)$ を得る。
:::
::::

::::exercise{#exr-qft-qed-8 level=advanced}
Møller 散乱の振幅二乗に $\theta=\pi/2$ を入れよ。

:::solution
$t=u=-s/2$ で対角項が $5$ ずつ、交差項が $8$ である。$2e^4(5+5+8)=36e^4$ となる。
:::
::::

::::exercise{#exr-qft-qed-9 level=advanced}
入射電子の質量を無視した $\mu^+\mu^-$ 生成について、$|\cos\theta|<a$ の角度領域への受容率を有限な終状態質量で求めよ。$s=16m_\mu^2$、$a=1/2$ で数値を得よ。

:::solution
$\beta=\sqrt{1-4m_\mu^2/s}$ とし、角度に依らない因子を除いた分布を $f(x)=1+x^2+(1-\beta^2)(1-x^2)=2-\beta^2+\beta^2x^2$ と書く。$\dd\Omega=2\pi\dd x$ だから受容率は
$$
\frac{(2-\beta^2)a+\beta^2 a^3/3}{2-\beta^2+\beta^2/3}.
$$
指定の $s$ では $\beta^2=3/4$。分子は $5/8+1/32=21/32$、分母は $3/2$、従って $7/16$ である。閾値へ近づく極限では分布が等方的になり、比は $a$ になる。一方、質量を落とした値は $a=1/2$ で $13/32$。分解能の効果を議論する場合は全断面積の質量補正と受容率の補正を別々に計算する。
:::
::::

## 確認テスト

::::quiz
Compton 散乱の木で電子線への光子の挿入順序はいくつあるか。

- [ ] 一つ
- [x] 二つ
- [ ] 三つ

:::explanation
入射光子を吸収する挿入と出射光子を放射する挿入のどちらを先に置くかで二つある。外部電流を保存するには両方を足す。
:::
::::

::::quiz{answer="0.3333333333" tol="1e-6"}
$\omega=m$、$\theta=\pi$ の $\omega'/\omega$ を答えよ。

:::explanation
$\omega'/\omega=[1+(\omega/m)(1-\cos\theta)]^{-1}=[1+2]^{-1}=1/3$ である。
:::
::::

::::quiz
Møller 散乱で二つの出射線を入れ替えた寄与の相対符号を選べ。

- [ ] 正
- [x] 負
- [ ] どちらでもよい

:::explanation
同種のフェルミオンの外部状態は生成演算子の積であり、二つの出射生成演算子を交換すると状態の符号が反転する。
:::
::::

::::quiz{answer="2" tol="1e-6"}
スカラー QED の二光子頂点の $iq^2\eta^{\mu\nu}$ に掛かる数を答えよ。

:::explanation
密度の $q^2 A_\rho A^\rho\varphi^*\varphi$ へ二つの外部光子を結ぶ対応が二通りあるので係数は二である。
:::
::::

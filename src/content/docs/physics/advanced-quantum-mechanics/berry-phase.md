---
title: 断熱近似と Berry 位相
description: ギャップのある時間発展を推定し、固有状態の幾何学的位相と非可換な平行移動を求める。
notation:
  - symbol: 'P(s)'
    meaning: 注目する固有空間への射影
  - symbol: '\mathcal A'
    meaning: Berry 接続
  - symbol: '\mathcal F'
    meaning: Berry 曲率
  - symbol: '\gamma_{\mathrm B}'
    meaning: Berry 位相
---

## この章の目標

- スペクトルギャップのもとで断熱誤差を評価する。
- 固有ベクトルの局所位相と閉路に沿う干渉位相を区別する。
- 二準位系の Berry 曲率を計算し、磁気単極子の形を導く。
- 核の遅い運動や縮退固有空間で幾何学的な効果を記述する。

## 断熱定理

時刻 $t$ における準位が他から離れているとき、Hamiltonian を十分ゆっくり変えると、状態は初めに占めた準位の瞬間固有空間に近いまま動く。ギャップは遷移先との位相の回転を速くし、異なる時刻の遷移振幅を打ち消し合わせる。以下では有限次元の閉じた系を扱い、$s=t/T\in[0,1]$、$\varepsilon=\hbar/T$ とする。従って方程式は $i\varepsilon U'=H(s)U$ である。

:::theorem[ギャップ付き断熱定理]{#thm-aqm-adiabatic}
$H(s)$ は $C^2$ 級の Hermite 行列族とし、$U(0)=I$ とする。連続な固有値の枝 $E(s)$ の多重度が一定で、他の固有値との距離が一様に $g>0$ 以上とする。この固有値への射影を $P(s)$ とすると、$\|(I-P(s))U(s)P(0)\|=\order{\varepsilon}$ が $s$ について一様に成り立つ。
:::

:::proof
$P(s)$ は $E(s)$ のみを囲む小さな閉曲線上で resolvent $(zI-H(s))^{-1}$ を $1/(2\pi i)$ 倍して積分すれば得られる。曲線を局所的に固定でき、逆行列を二度微分して $C^2$ 級である。$E(s)$ も $\tr(HP)/\rank P$ から $C^2$ 級である。$Q=I-P$、$K=[P',P]$ と置く。$P^2=P$ の微分から $PP'P=0$、$QP'Q=0$ であり、$[K,P]=P'$ となる。$K^\dagger=-K$ なので $U_{\rm A}'=(-iH/\varepsilon+K)U_{\rm A}$、$U_{\rm A}(0)=I$ の解はユニタリである。$U_{\rm A}^\dagger P U_{\rm A}$ の微分に $[H,P]=0$ と $P'=[K,P]$ を用いると零になる。従って $U_{\rm A}(s)P(0)=P(s)U_{\rm A}(s)$ である。

$R=Q(H-E)^{-1}Q$ と書く。ただし逆は $Q$ 上だけで取り、$P$ 上で零とする。$\|R\|\le1/g$ である。$X=RP'P+PP'R$ と置くと $[H,X]=QP'P-PP'Q=K$ になる。この $X$ は $C^1$ 級である。実際 $R=Q(H-E+gP)^{-1}Q$ として普通の逆行列を微分できる。すると

$$
\begin{aligned}
(U^\dagger U_{\rm A})'&=U^\dagger K U_{\rm A},\\
(U^\dagger X U_{\rm A})'
&=(i/\varepsilon)U^\dagger K U_{\rm A}
+U^\dagger(X'+XK)U_{\rm A}.
\end{aligned}
$$

第二の式を積分し第一に代入すると、$\|U-U_{\rm A}\|\le\varepsilon[2\sup\|X\|+\int_0^1(\|X'\|+\|XK\|)\dd s]$ である。$Q(s)U_{\rm A}(s)P(0)=0$ と組み合わせて結論を得る。
:::

この評価の定数には $R$ やその微分が入るため、ギャップが閉じる極限とは独立に取れない。$P'$ の大きさは $H'$ とギャップから見積もれる。$HP=EP$ を微分し左から $Q$、右から $P$ を掛けると $QP'P=-RH'P$ である。従って離れた準位への遷移が小さいための目安に $\hbar\|H'\|/(Tg^2)$ が現れる。

:::example[固有ベクトルが固定される場合]
$H(s)=\Delta(s)\sigma_z/2$ なら $P'=0$ で $K=X=0$ である。初期の $\sigma_z$ 固有状態は任意の速さでそのままの成分にとどまり、時間依存する位相しか得ない。
:::

:::example[一定の速さで回転する磁場]
$H(t)=(\Delta/2)(\cos(\Omega t)\sigma_z+\sin(\Omega t)\sigma_x)$ とする。$R=e^{-i\Omega t\sigma_y/2}$ で $\psi=R\varphi$ と置くと $H_{\mathrm{rot}}=\Delta\sigma_z/2-\hbar\Omega\sigma_y/2$ である。初期状態が $+$ なら反対側の瞬間固有状態の確率は

$$
\frac{\hbar^2\Omega^2}{\Delta^2+\hbar^2\Omega^2}
\sin^2\left(\frac{\sqrt{\Delta^2+\hbar^2\Omega^2}}{2\hbar}t\right)
$$

となる。$H_{\mathrm{rot}}^2$ がスカラーなので指数を正弦と余弦に分けて得る。$\hbar\Omega\ll|\Delta|$ なら確率は二次の小量である。
:::

### 補助項で瞬間固有空間を厳密に運ぶ

実時間で $P=P(t)$ とする。補助 Hamiltonian を $H_{\rm cd}=i\hbar[\dot P,P]$ と置くと、交換子が反 Hermite なのでこれは Hermite である。$H+H_{\rm cd}$ の時間発展 $W$ で $W^\dagger PW$ を微分すると、括弧の中は $\dot P-[[\dot P,P],P]=0$ となる。従ってこの発展は $T$ を大きくする極限に頼らず瞬間空間を厳密に運ぶ。二準位で $P=(I+\vb n\cdot\boldsymbol\sigma)/2$ なら $[\dot P,P]=(i/2)(\dot{\vb n}\times\vb n)\cdot\boldsymbol\sigma$ なので $H_{\rm cd}=(\hbar/2)(\vb n\times\dot{\vb n})\cdot\boldsymbol\sigma$ である。$\dot{\vb n}$ に垂直な成分が、基底の回転を打ち消すために必要となる。

:::example[回転系で落ちる項を加える]
本文の $\vb n=(\sin\Omega t,0,\cos\Omega t)$ では $\vb n\times\dot{\vb n}=(0,\Omega,0)$ である。補助項は $+\hbar\Omega\sigma_y/2$、これを回転座標へ移しても変わらないので、$-\hbar\Omega\sigma_y/2$ と相殺する。そこで回転系では $+\Delta\sigma_z/2$ のみとなり、初期 $+$ の遷移確率は厳密に零である。例えば $\hbar\Omega=\Delta/3$ では、補助項を除いた最大遷移確率は $1/10$ である。同じ経路でも制御によって漏れを取り除ける。
:::

### 漸近評価のギャップと端点の役割

非縮退の二準位で、$s$ 上の固有ベクトルの位相を $\braket{n}{n'}=0$ になるように選ぶ。内積が純虚なので、位相の微分をそれに合わせれば区間内で実現できる。$\psi=\sum_n c_n(s)e^{-i\int_0^s E_n/\varepsilon}\ket{n(s)}$ を代入して成分を取ると
$$
c_m'=-\sum_{n\ne m} c_n(s)\braket{m}{n'}\,
e^{i\int_0^s (E_m-E_n)\dd v/\varepsilon}.
$$
ギャップがある場合に $\braket{m}{n'}= \mel m{H'}n/(E_n-E_m)$ が成り立つ。これは $H\ket n=E_n\ket n$ を微分して $m$ へ射影した式で、同じ準位への対角微分とは異なる。初期 $n$ から他方 $m$ への先頭振幅を調べるため、第一反復で $c_n=1$ と置く。$\Delta=E_m-E_n$、$g=-\braket m{n'}$、$S(s)=\int_0^s\Delta$ とすると
$$
\begin{aligned}
\int_0^1 g(s)e^{iS(s)/\varepsilon}\dd s
&=\frac{\varepsilon}{i}\left[\frac{g(s)}{\Delta(s)}e^{iS(s)/\varepsilon}\right]_0^1\\
&\quad-\frac{\varepsilon}{i}\int_0^1
\left(\frac g\Delta\right)'e^{iS(s)/\varepsilon}\dd s.
\end{aligned}
$$
これで境界項に $\varepsilon\,H'/\Delta^2$ が出る。さらに $H$ が $C^3$ 級で、両端で $g=0$ なら、終点では厳密な遷移振幅も $\order{\varepsilon^2}$ である。実際、厳密な式では上の積分の $g$ を $gc_n$ に替える。:ref[thm-aqm-adiabatic] で $c_m=\order{\varepsilon}$ が一様に成り立つため、二準位の係数方程式より $c_n'=\order{\varepsilon}$ である。一度部分積分すると境界項が消え、$(g/\Delta)c_n'$ を含む積分は外の $\varepsilon$ と合わせてすでに二次である。残った $(g/\Delta)'c_n$ を振幅とする積分も、ギャップで割ってもう一度部分積分すれば $\order{\varepsilon}$ である。振幅とその一階微分が一様有界だからである。例えば両端で $H'=0$ ならこの条件を満たす。途中の時刻では上端の境界項が残り得るため、終点だけの評価と区別する。

:::example[同じギャップで始動の仕方を比較する]
$H=\tfrac\Delta2(\cos\theta\,\sigma_z+\sin\theta\,\sigma_x)$、$\theta(s)=\theta_f F(s)$ とする。実の固有ベクトルを使うと $|\braket{-}{+'}|=|\theta'|/2$ である。$F=s$ では端点の導関数はどちらも 1、$F=3s^2-2s^3$ ではどちらも零である。$\kappa=\Delta/\varepsilon$ として第一反復の振幅には、前者で $\int_0^1e^{i\kappa s}\dd s=(e^{i\kappa}-1)/(i\kappa)$、後者で $\int_0^1 6s(1-s)e^{i\kappa s}\dd s$ が掛かる。後者は部分積分二回で $6(e^{i\kappa}+1)/(i\kappa)^2-12(e^{i\kappa}-1)/(i\kappa)^3$ になる。例えば $\kappa=5\pi$ なら第一項は消え、絶対値は $24/(5\pi)^3$ になる。一定速度の $2/(5\pi)$ と比べ、同じ終点まで動かしても時間の形で終点の励起が変わる。
:::

### 円錐状の磁場を有限速度で一周する

$H(t)=\Delta[\sin\theta(\cos\Omega t\,\sigma_x+\sin\Omega t\,\sigma_y)+\cos\theta\,\sigma_z]/2$ とする。$R(t)=e^{-i\Omega t\sigma_z/2}$ で $\psi=R\varphi$ と置くと、回転系の一定のベクトルは $\vb b=(\Delta\sin\theta,0,\Delta\cos\theta-\hbar\Omega)$ である。初期の上枝の軸 $\vb n_0=(\sin\theta,0,\cos\theta)$ への垂直成分の大きさは $|\hbar\Omega\sin\theta|$ である。この成分だけが反対枝を作るので
$$
p_-(t)=
\frac{(\hbar\Omega\sin\theta)^2}{\Delta^2-2\Delta\hbar\Omega\cos\theta+(\hbar\Omega)^2}
\sin^2\left(\frac{|\vb b|t}{2\hbar}\right).
$$
$\vb b=0$ のときは回転系の発展が恒等演算子となり、確率を零とする。$\theta=0$ では初期の固有ベクトル自体が固定で確率は零である。$\theta\ne0$ でも $\Delta>0$ を固定して $\Omega$ を十分小さくすると二次の漏れとなる。$\Omega>0$ として $T=2\pi/\Omega$ では $R(T)=-I$ である。初期枝に戻る複素振幅は、回転系の指数の対角成分にこの $-1$ を掛けたものになる。

:::example[一周の振幅から断熱位相を抜き出す]
$\hbar\Omega/\Delta\to0$ で $\Delta>0$ とする。初期上枝の回転系の重みは $1-\order{(\hbar\Omega/\Delta)^2}$ なので、戻りの振幅は先頭で $-e^{-i|\vb b|T/(2\hbar)}$ である。$|\vb b|=\Delta-\hbar\Omega\cos\theta+\order{\Omega^2}$ だから動力学的因子 $e^{-i\Delta T/(2\hbar)}$ を割ると、残りは $-e^{i\pi\cos\theta}$ になる。偏角は $-\pi(1-\cos\theta)$ と同じであり、Berry 位相の球面からの結果を動力学的な解で再現する。有限速度では正枝に戻る確率が 1 であると仮定せず、振幅の絶対値も確認する。
:::

### 避けられた交差と Landau–Zener の確率

ギャップが狭い領域を通過するとき、初めに選んだ瞬間の固有空間から別の枝へ移り得る。傾き $v>0$、最小ギャップ $\Delta>0$ の線形模型
$$
H(t)=\frac12\begin{pmatrix}vt&\Delta\\\Delta&-vt\end{pmatrix},
\qquad -\infty<t<\infty
$$
を考える。$v$ の単位はエネルギーを時間で割った値である。時間に依存しない基底を :term[透熱基底]{en="diabatic basis" yomi="とうねつきてい"} と呼ぶ。この模型の遠い両端では瞬間固有状態が二つの透熱基底に近づくが、例えば上の成分のエネルギーは時間負の側と正の側で順位が逆になる。従って同じ透熱成分に残る確率は、反対側の断熱枝に移る確率でもある。

:::proposition[線形通過における Landau–Zener の確率]{#prop-aqm-landau-zener}
線形模型で時間を $-\infty$ から $+\infty$ へ通過するとき、同じ透熱成分に出る確率は
$$
p_{\rm d}=\exp\left(-\frac{\pi\Delta^2}{2\hbar v}\right)
$$
である。
:::

:::proof
上と下の振幅を $a,b$ とする。$b=2(i\hbar a'-vt\,a/2)/\Delta$ を第二の式に入れると
$$
a''+\left(\frac{v^2t^2+\Delta^2}{4\hbar^2}+\frac{iv}{2\hbar}\right)a=0
$$
になる。$z=e^{i\pi/4}\sqrt{v/\hbar}\,t$、$\delta=\Delta^2/(4\hbar v)$、$\nu=-i\delta$ とすれば
$$
a_{zz}+\left(\nu+\frac12-\frac{z^2}4\right)a=0.
$$

この方程式の解を積分から作る。$\Re\nu<0$ について
$$
D_\nu(z)=\frac{e^{-z^2/4}}{\Gamma(-\nu)}
\int_0^\infty e^{-zu-u^2/2}u^{-\nu-1}\dd u
$$
と置く。$F=e^{z^2/4}D_\nu$ とすれば方程式は $F''-zF'+\nu F=0$ になる。積分での被積分関数は $(u^2+zu+\nu)e^{-zu-u^2/2}u^{-\nu-1}=-\partial_u(e^{-zu-u^2/2}u^{-\nu})$ で、両端の境界値が零になる。$\Re z>0$ の固定された方向で $u=y/|z|$ と置けば、指数の $e^{-u^2/2}$ が 1 に近づき、Laplace 積分から $D_\nu(z)\sim e^{-z^2/4}z^\nu$ である。微分の漸近も積分の因子を増やして得られる。

$\Re\nu=0$ で必要な解は $D_\nu=zD_{\nu-1}-(\nu-1)D_{\nu-2}$ により定義できる。両辺が積分で定義できる領域では $u^{-\nu}e^{-zu-u^2/2}$ の部分積分で従うので、これが同じ解の継続となる。右辺の二項を微分しても同じ方程式を満たす。漸近の第一項は $z\cdot e^{-z^2/4}z^{\nu-1}$、第二項は二次数下なので同じ先頭を使える。原点では $y=u^2/2$ とする Gamma 積分と Gamma の倍角公式から
$$
D_\nu(0)=\frac{2^{\nu/2}\sqrt\pi}{\Gamma((1-\nu)/2)},
\qquad
D_\nu'(0)=-\frac{2^{(\nu+1)/2}\sqrt\pi}{\Gamma(-\nu/2)}
$$
になる。微分の式では $z=0$ で外の指数の微分が消え、積分に $-u$ が掛かる。漸化式による延長でも、Gamma の漸化式で同じ原点値を得る。

$t>0$ の遠方で $a=D_\nu(z)$ を選ぶ。$e^{-z^2/4}=e^{-ivt^2/(4\hbar)}$ であり、$2i\hbar a'$ の先頭と $vt\,a$ が相殺するので、下の成分 $b$ は $1/t$ の尺度で小さい。負の $t$ では $-z$ と $iz$ の偏角はそれぞれ $\pi/4,-\pi/4$ なので、$D_\nu(-z)$ と $B(z)=D_{-\nu-1}(iz)$ の漸近を使える。後者も同じ微分方程式を満たし、負の側で時間位相が逆、$a$ 自体は $1/|t|$ の尺度になる。従って前者の係数だけが $t\to-\infty$ の上成分を与える。

$D_\nu(z)=A D_\nu(-z)+C B(z)$ の係数 $A$ を原点の値と微分から決める。$k=D_\nu'(0)/D_\nu(0)$、$l=B'(0)/B(0)$ とすると、値の式で $C B(0)=(1-A)D_\nu(0)$、微分で $k=-Ak+(1-A)l$ である。従って $A=(1-k/l)/(1+k/l)$ となる。原点の式を入れ、Gamma の反射公式 $\Gamma(z)\Gamma(1-z)=\pi/\sin\pi z$ を使うと $k/l=i\tan(\pi\nu/2)$ である。よって $A=e^{-i\pi\nu}=e^{-\pi\delta}$ となる。$+\infty$ の $D_\nu(z)$ と $-\infty$ の $D_\nu(-z)$ の $z^\nu$ の絶対値の規格化は、偏角が同じなので揃う。最後に上成分だけを置いて逆に伝播した時の上成分の重みは $|A|^2=e^{-2\pi\delta}$ である。有限の左右の時刻でユニタリの逆の対角成分は順伝播の対角成分の共役である。遠方の極限でも両者の絶対値が同じなので、これが順伝播の滞在確率になる。
:::

:::example[最小ギャップの二乗が効く尺度]
$\Delta^2/(\hbar v)=2\log2/\pi$ とすれば、:ref[prop-aqm-landau-zener] の確率は $1/2$ になる。傾き $v$ を保ち $\Delta$ を二倍にすると指数が四倍となり、同じ透熱成分の確率は $1/16$ になる。有限時刻で準備・測定する場合は、初期・終点の固有ベクトルに加え、その有限区間の時間発展も必要になる。上で求めた独立な二解から二成分の解を復元して列に並べた行列を $F(t)$ とすれば、有限区間の発展は $F(t_f)F(t_i)^{-1}$ であり、その初期・終状態間の行列要素を使う。
:::

## Berry 位相

非縮退の固有状態 $\ket{n(s)}$ を滑らかに選ぶ。$\psi_{\rm A}(s)=e^{-i\int_0^s E(v)\dd v/\varepsilon}e^{i\gamma(s)}\ket{n(s)}$ と書く。断熱発展の補助項 $K$ は $PKP=0$ なので、瞬間固有状態を bra として方程式に掛ければ

$$
\gamma'(s)=i\braket{n(s)}{n'(s)}
$$

となる。内積の実部は規格化の微分から零であり、右辺は実数である。閉路で端点の射線が同じなら

$$
\gamma_{\mathrm B}=
i\int_0^1\braket{n(s)}{n'(s)}\dd s+\arg\braket{n(0)}{n(1)}
\pmod{2\pi}.
\label{eq-aqm-berry-phase}
$$

これが :term[Berry 位相]{en="Berry phase" yomi="べりーいそう"} である。端点の位相を揃えた表示なら第 2 項は零である。$\ket n\mapsto e^{i\chi}\ket n$ で積分は $-[\chi(1)-\chi(0)]$ だけ変わるが、第 2 項は同じ量を加えて相殺する。このため動力学的位相を差し引いた干渉で、閉路の位相を比較できる。

有限個の点 $s_j$ で重なりを取ると、各 $j$ で $\braket{n(s_j)}{n(s_{j+1})}=1+\braket{n}{n'}\Delta s+\order{(\Delta s)^2}$ である。最後に $\braket{n(1)}{n(0)}$ を掛けた積の偏角に負号を付けると、細分の極限が :eqref[eq-aqm-berry-phase] になる。有限分割でも位相因子は隣接する積の間で厳密に消える。

:::example[状態を回さず代表だけ変える]
$\ket{n(s)}=e^{2\pi i s}\ket{n_0}$ なら積分は $-2\pi$、端点位相は $0$ を法として零である。従って Berry 位相は零で、射線が定数であることと一致する。
:::

:::example[逆向きにたどる]
閉路 $C$ とその逆 $C^{-1}$ を同じ速度の大きさでたどると $\int E\,\dd t$ は同じで、接続の積分は逆符号である。重ね合わせで比較すると相対位相は $2\gamma_{\mathrm B}(C)$ になる。この操作で力学的な位相を共通因子として除ける。
:::

### 非閉路の端点と離散的な閉路

始点と終点の重なりが非零の開経路では、$\arg\braket{n(0)}{n(1)}+i\int\braket n{n'}\dd s$ はゲージ不変である。積分のゲージ変化は端点位相差の負、第一項はその正だからである。閉路の場合と異なり、端点の射線が一致することを前提にしない。この位相は端点の比較を含む干渉で使える。直交端点なら重なりが零であり、位相を定義するには別の比較経路または中間状態を指定する。

:::example[三つの射線の位相を重なりだけで決める]
$\ket{+z}=(1,0)$、$\ket{+x}=(1,1)/\sqrt2$、$\ket{+y}=(1,i)/\sqrt2$ をこの順につなぎ、最後に $+z$ に戻る。三つの重なりの積は $(1+i)/4$ で、負の偏角は $-\pi/4$ である。隣接する二状態を結ぶ短い大円上で、代表を実正の重なりを持つように取ると、その区間の接続は零になる。実際、端点の張る実二次元平面内に $\cos s\,u+\sin s\,v$ と書けば、規格化と $u,v$ の直交性で内積が零である。従って閉路の位相は区間の継ぎ目の位相だけで、先の積から得られる。球面の八分の一の面積 $\pi/2$ に $-1/2$ を掛けた値とも一致する。
:::

### 参照経路と同じ準位の位相を比較する

二つの経路を制御自由度 $\ket0,\ket1$ で分け、内部状態がともに同じ最終 $\ket{n_0}$ へ戻るとする。各経路の動力学的位相を $\delta_0,\delta_1$、幾何学的位相を $\gamma_0,\gamma_1$ と書く。終状態は $\ket{n_0}$ に共通で、制御の相対位相は $\delta_1-\delta_0+\gamma_1-\gamma_0$ である。追加の位相 $\zeta$ を一枝に入れて $\ket+$ へ射影すると
$$
P_+(\zeta)=\frac12[1+\cos(\zeta+\delta_1-\delta_0+\gamma_1-\gamma_0)]
$$
となる。戻りの軸が違う場合には、二つの最終の内部状態の重なりも干渉項に掛かる。したがって幾何学的な位相を比べたいときに、単に内部の始点だけを揃えて経路を分けるだけでは足りない。終点の状態、位相の基準、動力学的位相のどれが比較で共通になるかを指定する。

:::example[逆向きの経路で位相が二倍になる]
上枝の $\theta=\pi/3$ の緯線は立体角 $2\pi(1-\cos\theta)=\pi$ なので $\gamma=-\pi/2$ である。参照を同じ準位で静止させ、同じ時間と同じ固有値を持たせれば $\delta_1-\delta_0=0$ である。$\zeta=0$ の出力は $1/2$、$\zeta=\pi/2$ なら 1 になる。一方、参照を逆向きに一周させれば幾何学的な差は $-\pi$ となり $\zeta=0$ の出力が零になる。逆向き同士の比較だけでは二倍の差が $2\pi$ を法とするため、単一経路の位相は $\pi$ ずれた候補を区別できない。
:::

:::example[同じ閉路を二つのゲージで積分する]
$\ket{u_N}=(\cos(\theta/2),e^{i\phi}\sin(\theta/2))$ を、$\ket{u_S}=e^{-i\phi}\ket{u_N}$ に替える。$\theta=\pi/3$ の円では北の接続成分は $-1/4$、南では $3/4$ である。$0\le\phi\le2\pi$ を積分した値は $-\pi/2,3\pi/2$ で、指数はいずれも $-i$ である。この例で両ゲージの状態は区間の初めと終わりで周期的だが、位相関数 $-\phi$ を実数として見た増分は $-2\pi$ である。接続の線積分の数値そのものが一価と仮定する必要はない。
:::

### 一次元の二成分バンドで一周の符号を得る

格子の無次元波数 $k$ が $[-\pi,\pi]$ を一周するとし、$H(k)=[t_1+t_2\cos k]\sigma_x+t_2\sin k\,\sigma_y$、$t_1,t_2>0$ とする。$t_1\ne t_2$ なら零のギャップはない。$d_x+i d_y=t_1+t_2e^{ik}=|d|e^{i\varphi}$ と書けば、上枝の局所代表は $(1,e^{i\varphi})/\sqrt2$ で、接続は $-\partial_k\varphi/2$ である。一周の位相は $\gamma=-\pi w$、$w$ は複素数の曲線が原点を回る回数である。$t_2<t_1$ なら曲線は原点を含まず、$t_2\to0$ まで切れずに縮められるので $w=0$ である。$t_2>t_1$ なら中心を原点へ動かしても原点を横切らず、$w=1$ になる。$t_2=t_1$ では $k=\pi$ に縮退があり、この位相の変化を妨げるギャップが消える。位置や胞のとり方を固定した二成分の基底で比較する。$k$ に依存して基底自体を変える場合は、基底の接続も観測の位置の規約に入る。

:::example[波数を一周せず一部だけ積分する]
$t_1=0$ の極限なら $\varphi(k)=k$ であり、$0$ から $\pi/2$ の積分だけでは $-\pi/4$ になる。ただしこれはこの代表での接続の値である。開経路の位相を両端の同一基底の比較で定めるなら、重なりは $(1+i)/2$ で偏角 $+\pi/4$ が加わり、和は零である。全周の積分では端点の代表が一致して偏角の補正は零、結果は $-\pi$ になる。閉路と開経路で同じ線積分をそのまま観測値とみなせない。
:::

## Berry 接続と Berry 曲率

パラメータを $\vb R=(R^1,\ldots,R^d)$ とする。局所的に規格化した状態を選び、:term[Berry 接続]{en="Berry connection" yomi="べりーせつぞく"} $\mathcal A_i=i\braket n{\partial_i n}$ を定める。ゲージ変換で $\mathcal A_i\mapsto\mathcal A_i-\partial_i\chi$ である。この場合のゲージは固有ベクトルの位相の自由度を指す。:term[Berry 曲率]{en="Berry curvature" yomi="べりーきょくりつ"} は

$$
\mathcal F_{ij}=\partial_i\mathcal A_j-\partial_j\mathcal A_i
=i(\braket{\partial_i n}{\partial_j n}-\braket{\partial_j n}{\partial_i n})
=-2\Im\braket{\partial_i n}{\partial_j n}.
\label{eq-aqm-berry-curvature}
$$

二階微分は差で消える。$\mathcal F$ は位相に依存せず、曲面全体で滑らかなゲージがあるとき Stokes の定理により $\oint\mathcal A=\int\mathcal F$ となる。複数のゲージが必要なら継ぎ目の位相差も積分し、結果を $2\pi$ を法として扱う。

### エネルギー分母による曲率

固有方程式を微分し、異なる固有状態 $\ket m$ と内積を取ると $\braket m{\partial_i n}=\mel m{\partial_i H}n/(E_n-E_m)$ である。完全性を :eqref[eq-aqm-berry-curvature] に入れると $m=n$ の項は実で、虚部に寄与しない。よって

$$
\mathcal F_{ij}=-2\Im\sum_{m\ne n}
\frac{\mel n{\partial_i H}m\mel m{\partial_j H}n}{(E_n-E_m)^2}.
$$

従って他の準位とのギャップが小さいところでは曲率が大きくなり得る。テンソル $\braket{\partial_i n}{Q\partial_j n}$ の実部は、隣り合う射線の距離の二次項を与える。実際 $1-|\braket{n(\vb R)}{n(\vb R+\Delta\vb R)}|^2=\|Q n(\vb R+\Delta\vb R)\|^2$ を展開すればよい。

:::example[実の固有状態]
実対称 Hamiltonian の非縮退領域で実の規格化固有ベクトルを滑らかに取れるなら、規格化の微分から $\braket n{\partial_i n}=0$ で曲率も零である。しかし閉路全体で同じ符号に戻らないとき、端点の重なり $-1$ から位相 $\pi$ を得ることがある。
:::

:::example[純粋な位相の変化]
$\mathcal A$ に $-\dd\chi$ を加えても、成分表示で $\partial_i\partial_j\chi-\partial_j\partial_i\chi=0$ となる。小長方形の四辺の和を展開すれば、先頭は $\mathcal F_{ij}\Delta R^i\Delta R^j$ であり、一辺の接続からは読めない閉路効果である。
:::

### 曲率と量子計量を同じ行列で求める

$Q_{ij}=\braket{\partial_i n}{(I-P)\partial_j n}$ を :term[量子幾何テンソル]{en="quantum geometric tensor" yomi="りょうしきかてんそる"} とする。ゲージを変えて微分したときに生じる $\partial_i\chi\ket n$ の項は $I-P$ で消え、$Q$ は不変である。Hermite 性から $g_{ij}=\Re Q_{ij}$ は対称、$\mathcal F_{ij}=-2\Im Q_{ij}$ は反対称である。二点の重なりを二次まで展開すると $1-|\braket{n(\vb R)}{n(\vb R+\dd\vb R)}|^2=g_{ij}\dd R^i\dd R^j+\order{|\dd\vb R|^3}$ になる。一次の実部は規格化の微分で零、二次では $\Re\braket n{\partial_i\partial_j n}=-\Re\braket{\partial_i n}{\partial_j n}$ を使う。残った規格化方向の二乗が $P$ の引き算を与える。

:::example[二つの角を測る距離]
$\ket n=(\cos(\theta/2),e^{i\phi}\sin(\theta/2))$ なら $Q_{\theta\theta}=1/4$、$Q_{\phi\phi}=\sin^2\theta/4$、$Q_{\theta\phi}=i\sin\theta/4$ となる。最後は $\braket{\partial_\theta n}{n}=0$ と $\braket{\partial_\theta n}{\partial_\phi n}=i\sin\theta/4$ を使った。赤道で方位角を $\delta$ だけ変えた重なりは $(1+e^{i\delta})/2$ なので、確率の差は $\sin^2(\delta/2)=\delta^2/4+\cdots$ となる。同じ位置で曲率は $\mathcal F_{\theta\phi}=-1/2$ であり、距離の係数と違う成分から得る。
:::

### 小さな閉路の重なりから曲率を読む

点 $R$ の周りで $x,y$ の二座標をそれぞれ $\delta x,\delta y$ だけ進む長方形を正向きに回る。各辺の状態の重なりの偏角の負を順に足した値は、微小分割の極限で $\oint\mathcal A$ である。小さな一辺で一つの重なりを使う場合も、辺の中点の接続の値に辺の長さを掛けた評価に二次まで一致する。実際、$u(t+\delta)=u+\delta u'+\delta^2 u''/2+\cdots$ の重なりの対数の虚部を取ると、$(\braket u{u'})^2$ は実であり、虚部の二次は $\Im\braket u{u''}/2=\partial_t\Im\braket u{u'}/2$ である。従って四辺の総和の先頭は
$$
\left(\partial_x\mathcal A_y-\partial_y\mathcal A_x\right)\delta x\delta y
$$
である。各頂点の任意の位相選択は積の中で厳密に打ち消し合う。格子間隔が大きく、どこかの重なりが零に近づくと偏角が不安定なので、曲率の符号を読む前に間隔を減らす。

:::example[接続が零でも計量が残る実模型]
$H=X\sigma_z+Y\sigma_x$、$X,Y$ の原点を除く。局所の実単位ベクトルでは $\braket n{\partial_i n}=0$ なので局所曲率は零である。ところが $\varphi=\arg(X+iY)$ に対し $g_{\varphi\varphi}=1/4$ である。方位を $\delta$ だけ動かすと重なり二乗は $\cos^2(\delta/2)$ で、たとえば $\delta=\pi/3$ では $3/4$ となる。さらに原点を一周する大きな閉路ではゲージの継ぎ目を通り $-1$ が残る。局所曲率が零であることは、原点を除いた領域で一周の位相が零であることを意味しない。
:::

### 時間反転のあるバンドでの曲率の関係

波数の反転 $k\mapsto-k$ で $H(-k)=\Theta H(k)\Theta^{-1}$、$\Theta$ は波数に依存しない反ユニタリとする。非縮退な一枝が自身へ移る場合、$P(-k)=\Theta P(k)\Theta^{-1}$ である。$-k$ の変数に関する微分は元の微分の符号を一回ずつ替えるので、二つの交換子では相殺する。曲率の $i\tr P[\partial_iP,\partial_jP]$ による表示で、反ユニタリな共役はトレースを複素共役にする。$P$ は Hermite、交換子は反 Hermite なのでこのトレースは純虚で、最後に $\mathcal F_{ij}(-k)=-\mathcal F_{ij}(k)$ となる。積分領域が反転対称で二次元全域を覆うなら総 Chern 数は零になる。この仮定は枝が個別に対応する場合であり、時間反転で二枝が交換したり、点で縮退する場合に個々の枝へそのまま適用できない。

:::example[反転点で二つの制約がかかる二準位]
$H(k_x,k_y)=\cos k_x\,\sigma_x+\sin k_y\,\sigma_y+2\sigma_z$ とする。$K$ の共役で $\sigma_y$ だけの符号が変わるので時間反転条件を満たす。上枝の $(x,y)$ 曲率の分子 $\vb d\cdot(\partial_x\vb d\times\partial_y\vb d)$ は $-2\sin k_x\cos k_y$ である。従って $k_x=0$ で零、$k_x=\pm\pi/2,k_y=0$ で逆の値 $\pm1/8$ になる。局所の曲率が全て零でなくても、反転する二領域の寄与が総和で打ち消し合う。
:::

## スピンの Berry 位相と磁気単極子

$H=\vb R\cdot\boldsymbol\sigma$ とし、$R=|\vb R|>0$、$\theta,\varphi$ を $\vb R$ の球座標とする。上枝の北極で滑らかな状態は $\ket{+}_N=(\cos(\theta/2),e^{i\varphi}\sin(\theta/2))^\mathsf T$ である。Pauli 行列を掛けると固有値は $R$ である。$\theta$ 微分の内積は零、$\varphi$ 微分は $i\sin^2(\theta/2)$ なので

$$
\mathcal A_N=-\frac{1-\cos\theta}{2}\dd\varphi,\qquad
\mathcal F=-\frac12\sin\theta\,\dd\theta\wedge\dd\varphi .
\label{eq-aqm-spin-curvature}
$$

南極では $\ket{+}_S=e^{-i\varphi}\ket{+}_N$ を取り、$\mathcal A_S=\mathcal A_N+\dd\varphi$ となる。接続は二領域で違うが曲率は共通である。球全体の流束は $-2\pi$ なので、一つの滑らかな接続を球全体に置くことは Stokes の定理に反する。半径 $R$ の球面の面積要素が $R^2\sin\theta\dd\theta\dd\varphi$ であるため、曲率の双対ベクトルは $-\hat{\vb R}/(2R^2)$ となる。この :term[磁気単極子]{en="magnetic monopole" yomi="じきたんきょくし"} の形はパラメータ空間の原点の縮退に由来する。荷電粒子の [Aharonov–Bohm 位相](/physics/advanced-quantum-mechanics/path-integrals/) でも、局所接続の線積分を閉路の位相として読む。ここでの接続はスピンの固有状態から決まり、積分路は実空間の粒子軌道ではなく Hamiltonian のパラメータの曲線である。

スピン $j$ の $m$ については $\ket{n}=e^{-i\varphi J_z/\hbar}e^{-i\theta J_y/\hbar}\ket{jm}$ を使う。$\mathcal A_\varphi=\langle n|J_z|n\rangle/\hbar=m\cos\theta$、$\mathcal A_\theta=0$ である。この代表には $2\pi$ 回転の符号があるため、閉路全体では端点位相を加える。曲率は $-m\sin\theta\dd\theta\wedge\dd\varphi$ で、立体角 $\Omega$ の閉路では $\gamma_{\mathrm B}=-m\Omega\pmod{2\pi}$ である。曲面の取り方を変えて $\Omega$ に $4\pi$ を足しても $2m$ が整数なので結果は同じである。

:::example[余緯度が一定の閉路]
上枝の $\theta=\pi/3$ を正向きに一周すると $\Omega=2\pi(1-\cos\theta)=\pi$ で、位相は $-\pi/2$ である。赤道では $-\pi$ になり、両枝を比べる相対位相は $2\pi$ を法として零になる。
:::

:::example[結合の符号と枝]
$H=-\vb R\cdot\boldsymbol\sigma$ の基底は $\vb R\cdot\boldsymbol\sigma$ の上枝である。従ってその曲率は :eqref[eq-aqm-spin-curvature] の符号を保つ。$H$ の基底という語だけから符号を決めず、スピンの向きと Hamiltonian の符号を指定する。
:::

### 射影で表す曲率と Chern 数

固有ベクトルの位相を選ばずに使う式は

$$
\mathcal F_{ij}=i\,\tr(P[\partial_iP,\partial_jP]).
$$

実際、$P=\ketbra n n$ を微分すると $(\partial_iP)\ket n=(I-P)\partial_i\ket n$ となる。右辺の積は $\braket{\partial_i n}{(I-P)\partial_j n}$ から $i,j$ を交換したものを引く。$P$ 成分の寄与は積が実なので差から消え、:eqref[eq-aqm-berry-curvature] に一致する。

球面に縮退点がなく、球面上で非縮退の準位を考える。北・南の二つの半球のそれぞれで滑らかな固有状態を取れる。例えば半球の内部の極から全点へ向かう縮約経路に沿って $[P',P]$ で一つの状態を移せば、極を含む滑らかな規格化が得られる。赤道で $\ket n_S=e^{i\chi}\ket n_N$ と書くと $\mathcal A_S=\mathcal A_N-\dd\chi$ である。北半球の境界と南半球の境界は逆向きなので

$$
\int_{S^2}\mathcal F
=\oint_{\rm eq}(\mathcal A_N-\mathcal A_S)
=\oint_{\rm eq}\dd\chi\in2\pi\Z .
$$

$\chi$ そのものが一周の前後で同じ実数である必要はないが、位相 $e^{i\chi}$ が同じ値に戻るので増分は $2\pi$ の整数倍である。$C=(1/(2\pi))\int\mathcal F$ がこの線束の :term[Chern 数]{en="Chern number" yomi="ちゃーんすう"} である。ギャップが開いたまま連続に Hamiltonian を変えると積分は連続に変わるが整数値を取るので一定となる。スピン $1/2$ の上枝では $C=-1$、下枝では $C=+1$ である。閉路に応じて連続的に変わる Berry 位相と、閉曲面全体で量子化される積分を使い分ける。

### 異方性を持つ二準位の縮退点

$H(\vb R)=\vb d(\vb R)\cdot\boldsymbol\sigma$、$\vb d=M\vb R$ で $M$ を実正則行列とする。上枝の曲率は $\vb d$ 座標での単極子の引き戻しなので
$$
\mathcal F_{ij}=-\frac{\vb d\cdot(\partial_i\vb d\times\partial_j\vb d)}{2|\vb d|^3}
=-\frac{\det(M)\,\epsilon_{ijk}R^k}{2|M\vb R|^3}.
$$
第二式は三重積の行列式を比較した結果である。原点を囲む球の像は楕円体で、向きは $\sgn\det M$ で替わる。$\vb d$ 空間では原点を除きベクトル場 $-\vb d/(2|\vb d|^3)$ の発散が零なので、楕円体の束は内側の小球と同じである。従って上枝の Chern 数は $-\sgn\det M$ となる。

:::example[北極の局所曲率と全束の符号を調べる]
$M=\diag(2,1,-1)$ とする。$R_x=R_y=0,R_z=R>0$ では $\mathcal F_{xy}=1/R^2$ である。$M=I$ での $-1/(2R^2)$ に対して符号も大きさも変わる。一方、全束から得る Chern 数は $+1$ であり、局所曲率を球面積の一定値で外挿して求めるものではない。$\det M$ を負から正へ連続に変える間には $\det M=0$ となるため、そのままでは球面のギャップを保てない。
:::

### スピンの任意の磁気数に対する立体角

$\ket{n(\theta,\phi)}=e^{-i\phi J_z/\hbar}e^{-i\theta J_y/\hbar}\ket{jm}$ とする。$\theta$ 微分で出る $J_y$ の期待値は零である。$\phi$ 微分で出る $J_z$ の期待値は $\hbar m\cos\theta$ である。回した $J_z$ の横成分は梯子を持つので対角が零、縦成分は $\cos\theta$ を掛けた値になる。従ってこのゲージの接続は $m\cos\theta\,\dd\phi$、曲率は $-m\sin\theta\,\dd\theta\wedge\dd\phi$ である。$2\pi$ の回転が $(-1)^{2j}$ を与えるので緯線での端点位相も足す。$j-m$ が整数であることから、一周の位相は $-2\pi m(1-\cos\theta)$ を法 $2\pi$ で得る。球面を積分した Chern 数は $-2m$ で整数になる。

:::example[同じスピン 1 の三枝を比較する]
緯線 $\theta=\pi/3$ では立体角が $\pi$ であり、$m=1,0,-1$ の幾何位相は順に $-\pi,0,+\pi$ である。$m=\pm1$ の両相の指数はこの円では同じ $-1$ になるが、曲率と Chern 数は逆符号である。例えば立体角を $\pi/2$ の円に替えると二つの指数は $-i,+i$ と区別できる。閉路の一つで位相が同じであることから、各枝の曲率まで同じといえない。
:::

### 格子で複数の単極子の寄与を数える

二次元の周期波数で
$$
\vb d(k_x,k_y)=(\sin k_x,\sin k_y,m+\cos k_x+\cos k_y)
$$
とする。$m$ を変えると $\vb d$ が零となるのは $m=-2,0,2$ のいずれかである。$\sin k_x=\sin k_y=0$ の四点だけが候補であり、その第三成分が順に $m+2,m,m,m-2$ となるからである。上枝の北ゲージの接続は $-\frac12(1-d_z/|\vb d|)\dd\varphi$、$\varphi=\arg(d_x+i d_y)$ である。北極に対応する点では規則的で、南極に写る点では $\mathcal A_N\simeq-\dd\varphi$ となる。波数トーラスから南極への各逆像のまわりの小円を除き、Stokes の定理を使う。境界は各小円を時計回りにたどるので、曲率積分の $2\pi$ で割った値は、各小円を反時計回りに見た $\varphi$ の巻き数の和である。

小円の先頭で $(d_x,d_y)$ は $(\cos k_x^0\,\delta k_x,\cos k_y^0\,\delta k_y)$ である。よって巻き数は $\sgn(\cos k_x^0\cos k_y^0)$ で、四点の順を $(0,0),(\pi,0),(0,\pi),(\pi,\pi)$ とすれば $(+1,-1,-1,+1)$ である。$m>2$ では南極がなく $C_+=0$、$0<m<2$ では最後の点だけで $C_+=1$、$-2<m<0$ では中央二つと最後を足し $-1$ である。$m<-2$ では四点全部で和が零になる。曲率に局所的な $-\frac12$ の係数があるにもかかわらず、この閉曲面の積分は整数になる。$m$ がギャップを閉じる点ではこの上枝の線束を全域に定義できない。

:::example[一つと二つの閉じ方を比較する]
$m$ を $3$ から $1$ へ変えると、$(\pi,\pi)$ の $d_z$ が符号を替え、上枝の $C$ は $0$ から $1$ になる。続けて $m=-1$ へ動かすと中央の二点が同時に符号を替え、$C$ は $1$ から $-1$ へ二つ変わる。例えば二点の第三成分に小さく逆方向のずれが入るよう $\cos k_x-\cos k_y$ に比例する項を足せば、二つの交差の $m$ が分離する。各閉じ方が局所に持つ符号を数えることで、ジャンプの個数まで説明できる。
:::

## Born–Oppenheimer 近似と幾何学的位相

分子では核が電子より重く、核座標 $\vb R$ を固定した電子 Hamiltonian $H_e(\vb R)$ を先に解ける。核質量を $M$ とし、全波動関数を $\Psi(\vb R)=\chi(\vb R)\phi(\vb R)$、$H_e\phi=E_e(\vb R)\phi$ と近似する。電子について射影し、核の運動項 $-\hbar^2\nabla_R^2/(2M)$ の積微分を整理すると

$$
H_{\mathrm{eff}}=\frac{(-i\hbar\grad_R-\hbar\mathcal A)^2}{2M}
   +E_e(\vb R)+\Phi_{\mathrm{BH}},\qquad
\Phi_{\mathrm{BH}}=\frac{\hbar^2}{2M}\sum_i
 \braket{\partial_i\phi}{(I-P)\partial_i\phi}.
\label{eq-aqm-bo-effective}
$$

実際 $\braket{\phi}{\partial_i\phi}=-i\mathcal A_i$ と $\braket{\phi}{\partial_i^2\phi}=\partial_i\braket{\phi}{\partial_i\phi}-\|\partial_i\phi\|^2$ を使うと、$\partial_i-i\mathcal A_i$ の二乗との差は $\|(I-P)\partial_i\phi\|^2$ になる。正の補正 $\Phi_{\mathrm{BH}}$ は対角 Born–Huang 項である。電子位相 $\phi\mapsto e^{i\chi_0}\phi$ を替えると核振幅は $e^{-i\chi_0}\chi$ に変わり、共変微分の形が保たれる。この射影近似は電子間のギャップと核運動の尺度が分離する領域で使う。

:::example[円錐交差を囲む核運動]
実二準位で $H_e=X\sigma_z+Y\sigma_x$ とする。$(X,Y)=R(\cos\varphi,\sin\varphi)$ の上枝は $(\cos(\varphi/2),\sin(\varphi/2))^\mathsf T$ で、局所接続は零だが一周すると符号を変える。全波動関数を周期的にするには核の $\chi$ も反周期にする。核を角度自由度だけに制限すると $-i\partial_\varphi$ の固有値は整数 $n$ に $1/2$ を足した値となる。
:::

:::example[同じ円での対角補正]
上の実固有状態で $\braket{\phi}{\partial_\varphi\phi}=0$、$\|\partial_\varphi\phi\|^2=1/4$ なので、慣性モーメント $I$ の角運動だけを残した Born–Huang 項は $\hbar^2/(8I)$ である。反周期境界の運動項の最低値も $\hbar^2/(8I)$ なので、電子エネルギーを除いた模型の最低値は $\hbar^2/(4I)$ になる。
:::

### 円環上の有効模型と厳密模型を比べる

電子二準位と核の角運動を同時に残し $H=-\hbar^2\partial_\varphi^2/(2I)+\Delta(\cos\varphi\,\sigma_z+\sin\varphi\,\sigma_x)/2$ とする。$\Psi=e^{-i\varphi\sigma_y/2}\eta$ と置くと電子項は固定の $\Delta\sigma_z/2$ になる。周期的な $\Psi$ に対応して $\eta$ は反周期である。運動項は $\hbar^2(-i\partial_\varphi-\sigma_y/2)^2/(2I)$ になる。$m\in\Z+1/2$ の角モードで $J=\hbar^2/(2I)$ と書けば
$$
E_{m,\pm}=J(m^2+1/4)\pm\sqrt{\Delta^2/4+J^2m^2}.
$$
大きい $\Delta>0$ で上枝を選ぶと $+\Delta/2+J(m^2+1/4)+J^2m^2/\Delta+\cdots$ である。$J/4$ は Born–Huang 項で、その次が電子の他枝との混合から出る。射影近似で保持するものと、混合が小さい条件 $J|m|/\Delta\ll1$ をこの例から読める。

:::example[有効模型の誤差を最小角運動量で測る]
$\Delta=10J$、$m=1/2$ なら、上枝の平方根は $\sqrt{101}\,J/2$ であり、射影模型との差は $(\sqrt{101}-10)J/2\simeq0.02494J$ である。先頭の混合補正 $J^2m^2/\Delta=0.025J$ と近い。$m$ を大きくすると $J|m|$ に比例する非対角成分が増えるので、同じ電子ギャップでも精度は低下する。
:::

### 核と電子の積を全ての枝で展開する

固定した各 $\vb R$ の電子基底を列に持つユニタリ $N(\vb R)$ を考え、全波動関数を $\Psi=N\chi$ とする。核の微分を作用させると $N^\dagger\partial_i(N\chi)=\partial_i\chi+N^\dagger(\partial_i N)\chi=(\partial_i-i\mathcal A_i)\chi$ である。全基底が揃っていれば同じ置換を二階でも行え、核運動項は $-\hbar^2\sum_i(\partial_i-i\mathcal A_i)^2/(2M)$ となる。単一枝 $a$ の射影を最後に行うと、対角の共変微分に加え、二回の微分の間に別の枝 $b\ne a$ を通った項 $\hbar^2\sum_{i,b\ne a}|(\mathcal A_i)_{ba}|^2/(2M)$ が残る。これは :eqref[eq-aqm-bo-effective] の Born–Huang 項に等しい。最初から電子行列を対角に切り詰めた後で運動の二乗を取るとこの項が抜ける。

完全な基底の $\mathcal A=iN^\dagger\dd N$ には、$\partial_iN^\dagger=-N^\dagger(\partial_iN)N^\dagger$ により $\partial_i\mathcal A_j-\partial_j\mathcal A_i=i[\mathcal A_i,\mathcal A_j]$ が成り立つ。すなわち行列の曲率は零である。単一枝への射影でこの交換子の相手が減り、射影の接続は一般に非零の曲率を持つ。全空間で基底を変える自由度と、ある固有枝だけを選んだ状態の曲率を区別できる。

:::example[円錐の半径方向と角方向で対角補正を比べる]
$H_e=X\sigma_z+Y\sigma_x$ の実の上枝は半径 $R$ に依存しない。従って $R$ 微分からの Born–Huang 項は零である。Cartesian 微分の二乗和は $\|\partial_R\phi\|^2+R^{-2}\|\partial_\varphi\phi\|^2=1/(4R^2)$ なので、対角補正は $\hbar^2/(8MR^2)$ である。$R$ を半分にするとこのエネルギーは四倍、電子のギャップ $2R$ は半分であり、その比は八倍となる。局所の電子ギャップが正でも、交差へ近づくと単一枝の近似が急に悪くなり得る。
:::

### 核の境界条件とゲージを一つのスペクトルで照合する

円錐交差の周りの実の上枝 $\phi_r$ は反周期で $\mathcal A_r=0$ である。電子代表を $\phi_p=e^{i\varphi/2}\phi_r$ に替えれば、$\phi_p$ は周期的になる。核の係数は $\chi_p=e^{-i\varphi/2}\chi_r$ なので周期的で、接続は $\mathcal A_p=-1/2$ である。このゲージで運動項は $J(-i\partial_\varphi+1/2)^2$、$J=\hbar^2/(2I)$ である。周期固有関数 $e^{in\varphi}$ には $J(n+1/2)^2$ が掛かる。実の代表における反周期の $e^{i(n+1/2)\varphi}$ と同じエネルギーである。さらに単一枝の対角項 $J/4$ は $\partial_\varphi\phi$ のゲージ不変な直交部分のノルムから出るので、両表示で同じである。

:::example[エネルギー零の取り方と比を区別する]
運動項だけなら最小の二種 $|n+1/2|=1/2,3/2$ のエネルギーは $J/4,9J/4$ であり、比は 9 である。対角 Born–Huang の定数を加えた全模型では $J/2,5J/2$ になるので比は 5 になる。両者の差はどちらも $2J$ であり、定数の移動に依存しない。比を指定する問題では、電子の基準値や対角補正を入れたエネルギーを使うのか、それとも核の角運動エネルギーだけを使うのかを分ける必要がある。
:::

## Aharonov–Anandan 位相と非可換 Berry 位相

断熱近似で使ったのは、固有状態の射線が閉路を作ることである。一般の Schrödinger 発展でも $\ket{\psi(T)}=e^{i\theta}\ket{\psi(0)}$ なら、全位相から動力学的位相を除いた

$$
\gamma_{\rm AA}=\theta+\frac1\hbar\int_0^T\langle\psi(t)|H(t)|\psi(t)\rangle\dd t
$$

が定まる。これを :term[Aharonov–Anandan 位相]{en="Aharonov–Anandan phase" yomi="あはらのあなんだんいそう"} という。$i\braket{\psi}{\dot\psi}=\langle H\rangle/\hbar$ なので、:eqref[eq-aqm-berry-phase] の式を状態の実際の射線の閉路に適用したものと一致する。従って時刻ごとの位相の取り方に依存しない。

### 縮退した固有空間

$P(s)$ の階数が $r>1$ で、正規直交枠を列に並べた $N(s)$ を取る。断熱解を $\psi=e^{-i\int E/\varepsilon}N(s)c(s)$ とすると、$N^\dagger K N=0$ により $c'=i\mathcal A_s c$、$\mathcal A=iN^\dagger\dd N$ である。各点で $\mathcal A$ は Hermite だが、異なる点では可換とは限らない。従って輸送は微小ユニタリを時間順に掛けた $W=\mathcal P\exp(i\int\mathcal A)$ である。$N\mapsto NG$ に対し $c\mapsto G^\dagger c$ なので

$$
\mathcal A\mapsto G^\dagger\mathcal A G+iG^\dagger\dd G,\qquad
W\mapsto G(1)^\dagger W G(0).
$$

閉路の端点に同じ枠を選べばホロノミーは共役変換し、固有値が不変である。小長方形の四つの指数を順に二次まで掛けると $W_{\square}=I+i(\partial_i\mathcal A_j-\partial_j\mathcal A_i-i[\mathcal A_i,\mathcal A_j])\Delta R^i\Delta R^j+\cdots$ となる。行列の交換子項は閉路上で積の順序が効くことを示す。

:::example[静磁場での AA 位相]
$H=\hbar\omega\sigma_z/2$、初期の Bloch 球の余緯度を $\theta$ とする。周期 $2\pi/\omega$ で $U=-I$ となり $\theta_{\rm total}=\pi$ である。$\langle H\rangle=\hbar\omega\cos\theta/2$ なので $\gamma_{\rm AA}=\pi(1+\cos\theta)\equiv-\pi(1-\cos\theta)$ となり、同じ射線閉路の立体角に一致する。
:::

:::example[積の順序]
二つの輸送が $W_1=i\sigma_x$、$W_2=i\sigma_y$ なら $W_2W_1=i\sigma_z$、$W_1W_2=-i\sigma_z$ である。係数空間で二つの順序を比較すると符号が異なる。閉路の結果を比較するときは同じ始点の枠にそろえる。
:::

### 枠を選ばずに閉路の輸送を表す

射影の曲線が $P(1)=P(0)$ を満たすとする。$K=[P',P]$ で発展する $V'=KV$、$V(0)=I$ は、固有空間を運ぶユニタリである。$v(s)=V(s)v(0)$ が $P(s)$ の像にあるとき $P(s)v'(s)=PKPv(s)=0$ なので、微小変化は固有空間に直交する。逆に $Pv=v$ と $Pv'=0$ を満たす曲線は、$Pv=v$ の微分から $v'=P'v=Kv$ となる。従ってこの条件が平行移動を一意に定める。$v=N(s)c(s)$ と書けば $Pv'=0$ は $c'=i\mathcal A_s c$ である。端点に違う枠を使った場合、閉路の作用を始点の枠で比較する行列は $N(0)^\dagger N(1)W$ になる。例えば階数 1 ならこの端点の重なりを掛けることが、Berry 位相の端点項に当たる。

### 縮退部分空間の離散輸送

隣接する枠に対して $M_j=N_{j+1}^\dagger N_j$ と置く。細分が十分細かければ $M_j$ は可逆であり、極分解のユニタリ部分 $V_j=M_j(M_j^\dagger M_j)^{-1/2}$ を使える。$M_j=I-\Delta s\,N^\dagger N'+\order{\Delta s^2}=I+i\Delta s\,\mathcal A_s+\order{\Delta s^2}$ で、$M_j^\dagger M_j=I+\order{\Delta s^2}$ だからユニタリ化しても一次は変わらない。各 $V_j$ は第 $j$ 枠の係数から第 $j+1$ 枠の係数へ移す。$N_j\mapsto N_jG_j$ に対して $M_j\mapsto G_{j+1}^\dagger M_jG_j$ であり、逆平方根の共役性から $V_j$ も同じ変換をする。時間順に積を取れば中間の $G_j$ が相殺し、閉路で端点の枠を合わせた積の固有位相が残る。

:::example[一定の部分空間で生じる見かけの回転]
二次元部分空間そのものを固定し、枠だけ $N(s)=N(0)e^{-i\pi s\sigma_y/2}$ と回す。このとき $\mathcal A_s=\pi\sigma_y/2$、係数の輸送は $W=e^{i\pi\sigma_y/2}$ である。終点の物理ベクトルは $N(1)Wc(0)=N(0)c(0)$ であり、枠の回転と係数の逆回転が消える。閉路で端の枠を揃えずに $W$ の固有値 $\pm i$ だけを見ても幾何学的位相は得られない。
:::

### 行列の曲率を補空間で確かめる

$P=NN^\dagger$ とする。$\mathcal A_i=iN^\dagger\partial_i N$ を微分して引けば、二階微分は消えて $i[(\partial_iN^\dagger)\partial_jN-(\partial_jN^\dagger)\partial_iN]$ である。ここに $I=P+(I-P)$ を挿入する。$P$ を通る部分では $(\partial_iN^\dagger)N=-N^\dagger\partial_iN$ により $i[\mathcal A_i,\mathcal A_j]$ が出る。これを引いた残りは
$$
\mathcal F_{ij}
=i\big[(\partial_iN^\dagger)(I-P)(\partial_jN)
 -(\partial_jN^\dagger)(I-P)(\partial_iN)\big].
$$
さらに $(\partial_iP)N=(I-P)\partial_iN$ なので、$iN^\dagger[\partial_iP,\partial_jP]N$ と等しい。基底を $NG$ に替えると $\partial_iG$ を含む成分は補空間で消えるから、曲率は $G^\dagger\mathcal F_{ij}G$ と変換する。接続自体と異なり、余分な微分項は伴わない。

:::example[三次元の補空間を通って偏極が回る]
実の単位ベクトル $\vb n$ が選ぶ一方向に $H=\Delta\ketbra{\vb n}{\vb n}$、$\Delta>0$ を取り、直交する二次元の固有空間を運ぶ。球面上の接線の単位ベクトル $(\vb e_\theta,\vb e_\phi)$ を枠にすると、$\partial_\phi\vb e_\theta=\cos\theta\,\vb e_\phi$、$\partial_\phi\vb e_\phi=-\cos\theta\,\vb e_\theta-\sin\theta\,\vb n$ である。従って $\mathcal A_\phi=i\cos\theta\begin{pmatrix}0&-1\\1&0\end{pmatrix}$、$\mathcal A_\theta=0$ である。緯線の輸送では行列が一定なので、角 $-2\pi\cos\theta$ の実回転になる。$\theta=\pi/3$ では $-I_2$ になる。これは縮退の輸送の例だが、この枠で全接続成分が同じ生成子に比例し、非可換な順序効果は出ない。縮退しているだけで二つの閉路の積が必ず非可換になるわけではない。
:::

### 縮退の空間で本当に順序が効く模型

三次元の基底を $\ket1,\ket2,\ket3$ とし、実の小パラメータ $x,y,z$ で
$$
\ket{n(x,y,z)}=\frac{x\ket1+(y+iz)\ket2+\ket3}{\sqrt{1+x^2+y^2+z^2}}
$$
とする。Hamiltonian $\Delta\ketbra nn$ の零の固有空間は $P=I-\ketbra nn$、階数 2 でギャップは $\Delta>0$ の一定である。原点では枠を $\ket1,\ket2$ とする。$\partial_iP=-\ket{v_i}\bra3-\ket3\bra{v_i}$ であり、$v_x=\ket1,v_y=\ket2,v_z=i\ket2$ である。二つの微分の交換子を $P$ で挟むと $\ketbra{v_i}{v_j}-\ketbra{v_j}{v_i}$ となる。従って曲率の原点での値は $\mathcal F_{xy}=-\sigma_y$、$\mathcal F_{xz}=\sigma_x$ である。二つの小長方形の輸送は、それぞれ面積を $A_1,A_2$ として $W_1=I-iA_1\sigma_y+\cdots$、$W_2=I+iA_2\sigma_x+\cdots$ になる。始点の枠を揃えた上で差 $W_1W_2-W_2W_1$ を取ると、先頭は $-2i A_1 A_2\sigma_z$ である。従って縮退していても単に二つのスカラー位相を別々に足す記述ではこの差を表せない。辺の長さをともに $\varepsilon$ とした場合、一閉路の誤差は $\order{\varepsilon^3}$ だが、差では恒等成分との交換が消え、次の誤差は $\order{\varepsilon^5}$ である。

:::example[枠を変えても順序効果を比較する]
始点の枠を $G$ で回せば差全体は $G^\dagger(W_1W_2-W_2W_1)G$ になる。その Hilbert–Schmidt ノルムは共役で不変で、先頭は $2\sqrt2|A_1A_2|$ である。たとえば両面積を同じ $A$ とすれば差の大きさは $A^2$ の次数になる。一つずつの固有位相だけを記録すると、二つの閉路の回転軸の相対方向が分からない。順序を変えた輸送を測るには同じ始点の準備を用いて、その相対方向も保持する。
:::

### 射線が閉じるが瞬間の固有線をたどらない例

静的な有限次元 Hamiltonian の二固有状態 $\ket{E_0},\ket{E_1}$ を $|\psi(0)\rangle=\sqrt{1-p}\ket{E_0}+\sqrt p\ket{E_1}$ とし、$E_1>E_0$ とする。両成分の相対位相が一周する $T=2\pi\hbar/(E_1-E_0)$ で、全体の位相は $-E_0T/\hbar$ を法 $2\pi$ で得る。時間積分した期待エネルギーは $[(1-p)E_0+pE_1]T$ なので、AA 位相は $2\pi p$ となる。この結果は Hamiltonian の基底そのものが時間変化しない場合にも非零になり、断熱の固有線の接続とは違う射線を使っている。二次元の部分空間を Bloch 球に対応させれば、緯線の向きと面積からも同じ指数が出る。

:::example[力学の基準を変えても残る位相]
$p=1/4$ なら AA 位相は $\pi/2$ である。$H\mapsto H+\Lambda I$ にすると全位相が $-\Lambda T/\hbar$ だけ余分にずれるが、期待エネルギーの積分が $+\Lambda T/\hbar$ を加えるので $\pi/2$ は保たれる。全位相だけを幾何学的な量として記録すると、エネルギーの原点に依存する数になってしまう。
:::

## 演習問題

::::exercise{level=basic}
局所状態を $\ket n\mapsto e^{i\chi(\vb R)}\ket n$ と変えたとき接続の変化を求めよ。
:::solution
微分で $\partial_i\ket n$ に加えて $i(\partial_i\chi)\ket n$ が出る。bra の位相は相殺し、先頭の $i$ と掛けて $\mathcal A_i\mapsto \mathcal A_i-\partial_i\chi$ となる。
:::
::::

::::exercise{level=basic}
上枝スピン $1/2$ の $\theta=\pi/2$ の閉路の位相を求めよ。
:::solution
北極ゲージで $\mathcal A_\varphi=-(1-\cos\theta)/2=-1/2$ である。$0\le\varphi\le2\pi$ で積分して $-\pi\pmod{2\pi}$ となる。
:::
::::

::::exercise{level=basic}
スピン $1$ の $m=1$ が立体角 $\pi/2$ の閉路をたどる。幾何学的位相を求めよ。
:::solution
曲率の積分は $-m\Omega$ である。従って位相は $-\pi/2$、振幅の因子は $e^{-i\pi/2}=-i$ である。
:::
::::

::::exercise{level=basic}
$H$ の固有状態で時間発展する閉じた射線の AA 位相を求めよ。
:::solution
初期の固有値が $E$ で一定なら $\theta=-ET/\hbar$ である。期待値の積分は $ET/\hbar$ なので合計は零である。
:::
::::

::::exercise{level=advanced}
断熱定理の $K=[P',P]$ が射影を運び、対角ブロックが零となることを示せ。
:::solution
$P'=P'P+PP'$ から $PP'P=0$、$QP'Q=0$ である。すると $PKP=QKQ=0$ であり、$KP-PK=P'P+PP'=P'$ となる。$U'=KU$ に対し $(U^\dagger P U)'=U^\dagger(P'-[K,P])U=0$ なので射影が輸送される。
:::
::::

::::exercise{level=advanced}
$H=\vb R\cdot\boldsymbol\sigma$ の上枝で $\vb R=(0,0,R)$、$R>0$ における $\mathcal F_{xy}$ をエネルギー分母で計算せよ。
:::solution
上枝は $\ket+$、下枝は $\ket-$、ギャップは $2R$ である。$\mel +{\sigma_x}-=1$、$\mel -{\sigma_y}+=i$ なので和は $i/(4R^2)$ である。虚部に $-2$ を掛けると $-1/(2R^2)$ となる。これは北極方向の曲率ベクトルと一致する。
:::
::::

::::exercise{level=advanced}
反周期な核の角度関数で、最も低い二つの異なる運動エネルギーの比を求めよ。
:::solution
固有波は $e^{i(n+1/2)\varphi}$ であり、絶対値の最小は $1/2$、次は $3/2$ である。運動エネルギーは二乗に比例するので、高い方を低い方で割ると $9$ となる。各値で逆向きの角運動量が縮退する。
:::
::::

::::exercise{level=advanced}
縮退空間の輸送 $c'=i\mathcal A_s c$ がノルムを保つことと、枠を変えた場合の輸送の変換式を示せ。
:::solution
$\mathcal A_s^\dagger=\mathcal A_s$ なので $(c^\dagger c)'=-ic^\dagger\mathcal A_s c+ic^\dagger\mathcal A_s c=0$ である。枠を $N(s)G(s)$ に替えると $c_{\rm new}(s)=G(s)^\dagger c(s)$ である。$c(1)=Wc(0)$ を代入して $W_{\rm new}=G(1)^\dagger WG(0)$ を得る。
:::
::::

::::exercise{level=advanced}
上枝のスピン $1/2$ を移動させる。経路 $C$ を北極から方位角零の子午線で赤道まで進み、赤道を方位角 $\varphi=\pi/3$ まで進んで、その方位の子午線を北極へ戻るものとする。北ゲージでの接続の線積分、曲率による面積積分、三頂点の重なりの積から位相をそれぞれ求めよ。

:::solution
北ゲージで $\ket{n}=(\cos(\theta/2),e^{i\phi}\sin(\theta/2))^\mathsf T$ とする。接続は $\mathcal A=-\tfrac12(1-\cos\theta)\dd\phi$ であり、二つの子午線では $\dd\phi=0$ で寄与しない。赤道では $\cos\theta=0$ なのでそこから $-\varphi/2=-\pi/6$ を得る。北極で代表は方位に依らず $(1,0)$ なので、端点位相は零である。

曲率は $-\sin\theta\,\dd\theta\wedge\dd\phi/2$ である。向きは経路の順を採ると領域 $0<\theta<\pi/2,0<\phi<\varphi$ の正向きの境界になり、積分は $-\frac12\varphi\int_0^{\pi/2}\sin\theta\dd\theta=-\varphi/2$ である。立体角は $\varphi$ である。

各頂点の代表を $u_0=(1,0)$、$u_1=(1,1)/\sqrt2$、$u_2=(1,e^{i\varphi})/\sqrt2$ とすれば、$\braket{u_0}{u_1}\braket{u_1}{u_2}\braket{u_2}{u_0}=(1+e^{i\varphi})/4$ である。この偏角は $\varphi/2$ なので、負号を付けて $-\pi/6$ である。三辺はいずれも短い大円弧であり、一辺の輸送位相を端点の重なりで置き換えられる。辺を別の曲線に変えた場合は一般に頂点だけでは位相が決まらない。
:::
::::

::::exercise{level=advanced}
線形 Landau–Zener 模型で、最小ギャップ $\Delta$、駆動の傾き $v$ とする。同じ透熱成分に残る確率が $1/4$ のとき、傾きを二倍にした場合の確率と、同じ傾きでギャップを二倍にした場合の確率を求めよ。速い駆動で断熱枝の滞在がどう変わるか説明せよ。

:::solution
:ref[prop-aqm-landau-zener] より $\pi\Delta^2/(2\hbar v)=\log4$ である。傾きを二倍にすると指数の絶対値が半分になり、確率は $e^{-(\log4)/2}=1/2$ である。傾きを元のまま、ギャップだけ二倍にすると指数の絶対値は四倍で $e^{-4\log4}=1/256$ である。第一の変更は通過を速くし、第二は準位を離すため、結果は逆向きに変化する。

時刻を負の遠方から正の遠方へ渡ると、同じ透熱成分のエネルギーの順位は入れ替わる。従って初めの瞬間固有値の上または下の順位を最後まで保つ確率は $1-p_{\rm d}$ である。元の値は $3/4$、速い駆動では $1/2$、大きいギャップでは $255/256$ となる。この比較は線形通過で遠方の初期・終条件を使った場合であり、初期の状態が両成分を重ねたものなら相対位相も必要になる。
:::
::::

## 確認テスト

::::quiz
ギャップ付きの非縮退断熱近似で、閉路の位相のうち $\int E(t)\dd t$ によって決まるものはどれか。

- [x] 動力学的位相
- [ ] Berry 曲率
- [ ] 固有ベクトルのゲージ変換そのもの

:::explanation
エネルギーによる時間発展の位相は $-\int E(t)\dd t/\hbar$ である。Berry 位相には固有状態の重なりに沿う接続の積分と端点の位相が入る。
:::
::::

::::quiz{answer="-0.5" tol="1e-9"}
$\vb R\cdot\boldsymbol\sigma$ の上枝で北極方向、半径 $R$ の点の $R^2\mathcal F_{xy}$ を求めよ。
:::explanation
曲率ベクトルは $-\hat{\vb R}/(2R^2)$ なので北極で $\mathcal F_{xy}=-1/(2R^2)$ である。$R^2$ を掛けた値は $-1/2$ となる。
:::
::::

::::quiz
Berry 接続の曲率の局所値を、固有ベクトルの位相選択を替えた場合について答えよ。

- [x] 変わらない
- [ ] 必ず符号が替わる
- [ ] 任意の値にできる

:::explanation
接続の変化は $-\partial_i\chi$ で、曲率はその反対称な微分を取る。滑らかな $\chi$ の混合微分は可換なので、その差は零となる。
:::
::::

::::quiz{answer="9" tol="1e-9"}
半整数の角運動量で角運動エネルギーがその二乗に比例する模型において、絶対値が $3/2$ と $1/2$ の二状態のエネルギー比を求めよ。
:::explanation
共通係数は消える。比は $(3/2)^2/(1/2)^2=9$ である。電子の固定エネルギーなど角運動以外の加算項はこの問題では含めない。
:::
::::

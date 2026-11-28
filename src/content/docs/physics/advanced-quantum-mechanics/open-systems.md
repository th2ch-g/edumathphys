---
title: 密度行列と開放量子系
description: 密度演算子と量子チャネルを定義し、環境との相互作用からマルコフ型の緩和と位相緩和を導く。
notation:
  - symbol: '\rho'
    meaning: 密度演算子
  - symbol: '\tr_B'
    meaning: 部分系 B についての部分トレース
  - symbol: '\mathcal E'
    meaning: 量子チャネル
  - symbol: '\mathcal D[L]'
    meaning: 演算子 L に対応する散逸項
---

## この章の目標

- 密度行列で統計混合と部分系の状態を表す。
- Schmidt 分解と純粋化を具体的な行列から構成する。
- 完全正値性を Kraus 表現と環境模型から判定する。
- Born–Markov 近似の条件を明示し、緩和とデコヒーレンスを計算する。

## 純粋状態と混合状態

状態 $\ket{\psi_i}$ を確率 $p_i$ で準備する。測定が Hermite 行列 $A$ に対応するなら、全平均は $\sum_i p_i\mel{\psi_i}A{\psi_i}$ になる。複数の準備法が同じ測定結果を与えるか比較するには、測定に先立ってこの和を演算子としてまとめる。

:::definition[密度演算子]{#def-aqm-density}
有限次元 Hilbert 空間の :term[密度演算子]{en="density operator" yomi="みつどえんざんし"} は $\rho=\rho^\dagger\ge0$、$\tr\rho=1$ を満たす演算子である。準備の統計混合に対して $\rho=\sum_i p_i\ketbra{\psi_i}{\psi_i}$、期待値は $\langle A\rangle=\tr(\rho A)$ である。
:::

規格化状態の射影は正値でトレース 1 なので、統計混合も条件を満たす。逆に任意の密度演算子は固有基底で展開すれば固有値を確率とする混合になる。測定確率は $A$ を測定結果への射影に取って得られる。すべての射影で同じ確率なら、両密度演算子の差の任意の単位ベクトルにおける期待値が零である。差の固有ベクトルを使えば差は零である。

純粋状態は階数 1 の密度演算子である。固有値 $\lambda_i\ge0$、$\sum_i\lambda_i=1$ なので、$\tr\rho^2=\sum_i\lambda_i^2=1$ となるには一つの固有値が 1 で残りが零でなければならない。従って $\rho^2=\rho$、$\tr\rho^2=1$、階数 1 は同値である。$d$ 次元で Cauchy–Schwarz を使うと $\tr\rho^2\ge1/d$ で、等号は $\rho=I/d$ のときである。$S(\rho)=-k_B\tr(\rho\log\rho)$ は von Neumann エントロピーであり、$0\log0=0$ とする。固有値表示から純粋なら零、一様なら $k_B\log d$ になる。また $-x\log x$ の凹性から $\sum_i(-\lambda_i\log\lambda_i)\le d[-(1/d)\log(1/d)]$ なので、一様状態はエントロピーを最大にする。

閉じた系で $\rho(t)=U(t)\rho(0)U(t)^\dagger$ を微分すると $i\hbar\dot\rho=[H,\rho]$ となる。ユニタリ共役は固有値を保つので、純度もエントロピーも一定である。

純粋状態は密度演算子全体の凸集合の端点でもある。$\ketbra{\psi}{\psi}=q\rho_1+(1-q)\rho_2$、$0<q<1$ なら、$\psi$ に直交する $v$ での期待値は零である。両方の密度が正なので $\|\sqrt{\rho_j}v\|^2=0$、従って $\rho_jv=0$ となる。規格化から $\rho_1=\rho_2=\ketbra{\psi}{\psi}$ である。階数が 2 以上なら正の固有値二つに十分小さな量を正・負に振り分けた二つの異なる密度を作り、その等混合として元の密度を表せる。

:::example[異なる準備が同じ密度行列を作る]
$\ket0,\ket1$ を等確率で混ぜると $I/2$ である。$\ket\pm=(\ket0\pm\ket1)/\sqrt2$ を等確率で混ぜても、非対角成分が相殺して $I/2$ になる。従って単体系を測ってどちらの準備法かを判定できない。
:::

:::example[Bloch ベクトル]
$\rho=(I+\vb r\cdot\boldsymbol\sigma)/2$ の固有値は $(1\pm|\vb r|)/2$ なので $|\vb r|\le1$ である。$r_z=1/2,r_x=r_y=0$ なら確率は $3/4,1/4$、純度は $9/16+1/16=5/8$ になる。
:::

### 二準位で密度行列を再構成する

二準位の Hermite 行列は $\rho=(I+\vb r\cdot\boldsymbol\sigma)/2$ と書ける。トレースで $I$ の係数が固定され、$r_i=\tr(\rho\sigma_i)$ は実である。$(\vb r\cdot\boldsymbol\sigma)^2=|\vb r|^2I$ とトレース零から固有値は $(1\pm|\vb r|)/2$ になる。よって正値性は $|\vb r|\le1$ と同値であり、純度は $(1+|\vb r|^2)/2$ である。各軸の射影測定の $+1$ の頻度 $p_i$ を母集団から得れば、$r_i=2p_i-1$ で再構成できる。有限標本で $|\vb r|>1$ が出れば、その推定値をそのまま状態と解釈せず、統計誤差を含む推定を行う。

:::example[同じ対角成分で固有値が変わる場合]
$(p_x,p_y,p_z)=(3/4,1/2,3/4)$ なら $\vb r=(1/2,0,1/2)$ で、固有値は $(1\pm1/\sqrt2)/2$、純度は $3/4$ である。$z$ の測定だけなら、同じ対角成分を持つ $\diag(3/4,1/4)$ と区別できない。後者の純度は $5/8$ であり、違いの $1/8$ は非対角成分の絶対値二乗を二回足した値である。
:::

### 混合の分解と測定可能な距離

二準位の任意の密度行列 $\rho=(I+r\hat{\vb n}\cdot\boldsymbol\sigma)/2$ は、軸 $\hat{\vb n}$ の二つの純粋射影を重み $(1\pm r)/2$ で混合したものである。別に、互いに逆の横成分を持つ二単位ベクトル $\vb n_\pm=r\hat{\vb n}\pm\sqrt{1-r^2}\hat{\vb m}$、$\hat{\vb m}\perp\hat{\vb n}$ を取ると、対応する二射影を等確率で混ぜても同じ $\rho$ になる。期待値は跡に対して線形なので、どちらの準備かを系だけの測定で判定できない。混合の重みの Shannon エントロピーも分解ごとに異なる。このため量子状態の von Neumann エントロピーは分解ではなく、固有値から $-\sum_i\lambda_i\log\lambda_i$ と定める。

二密度 $\rho,\sigma$ の差 $D=\rho-\sigma$ は跡零の Hermite 行列である。効果 $0\le E\le I$ による一結果の確率差は $\tr(ED)$ である。$D=D_+-D_-$ を正負のスペクトル部分に分けると、それぞれの正行列の跡は等しく $\frac12\tr|D|$ に等しい。不等式 $0\le\tr(E D_\pm)\le\tr D_\pm$ から
$$
\max_{0\le E\le I}|\tr(E(\rho-\sigma))|
=\frac12\tr|\rho-\sigma|
$$
であり、達成するには $D_+$ の台への射影を選べばよい。右辺を :term[トレース距離]{en="trace distance" yomi="とれーすきょり"} という。二候補が等確率なら、$E$ で $\rho$ と判定し $I-E$ で $\sigma$ と判定する成功確率は $(1+\tr(ED))/2$ なので、最適値は $(1+\frac12\tr|D|)/2$ になる。

:::example[二つの Bloch ベクトルを一回の測定で区別する]
候補を $\rho=(I+\sigma_z)/2$ と $\sigma=(I+\sigma_x)/2$ とする。$D=(\sigma_z-\sigma_x)/2$ の固有値は $\pm1/\sqrt2$ である。正の軸 $(\hat{\vb z}-\hat{\vb x})/\sqrt2$ で測ると、成功率は $(1+1/\sqrt2)/2\simeq0.854$ になる。$z$ で測る方法では、$\rho$ を確実に同定する一方で $\sigma$ の半分を誤るため、等事前確率で $3/4$ である。二候補それぞれの測定を完全にする軸と、二候補の差を最大にする軸は同じとは限らない。
:::

### 純粋状態までの距離を固有値で読む

$\rho$ の純度は、最大固有値 $\lambda_{\max}$ に対して $\lambda_{\max}^2\le\tr\rho^2\le\lambda_{\max}$ を満たす。下側は固有値の二乗の和の一項を取ったもの、上側は $\lambda_i^2\le\lambda_{\max}\lambda_i$ を足したものである。また純粋な射影 $P=\ketbra{\psi}{\psi}$ への射影の確率 $\tr(P\rho)$ は最大で $\lambda_{\max}$ である。固有基底で各確率が固有値の凸結合になるからである。二準位では、Bloch ベクトルの長さが $r$ なら最も近い純粋状態を同じ軸に選べて、トレース距離は $(1-r)/2$ である。実際、差の行列の二固有値は $\pm|\vb r-\vb n|/2$ であり、単位 $\vb n$ との距離の最小は $1-r$ である。

:::example[純度が同じ二状態を区別する]
$\rho_z=(I+\sigma_z/2)/2$、$\rho_x=(I+\sigma_x/2)/2$ では固有値が同じ $(3/4,1/4)$ で純度は $5/8$ である。最も近い純粋状態はそれぞれ $\ket{+z},\ket{+x}$ で距離は $1/4$ である。一方、$\rho_z$ と $\rho_x$ の距離は $\sqrt2/4$ であり、等事前確率での最適な判定の成功率は $(1+\sqrt2/4)/2$ になる。純度は一つの状態の固有値の情報であり、二状態を区別するには相対的な基底のずれを使う。
:::

## 縮約密度行列と量子もつれ

全体が $\rho_{AB}$ で、実験で $A$ の観測量だけを測る場合を考える。$\braket{\vb b_j}{\vb b_k}=\delta_{jk}$ となる $B$ の基底を選ぶと

$$
\rho_A=\tr_B\rho_{AB}
=\sum_j(I_A\otimes\bra{b_j})\rho_{AB}(I_A\otimes\ket{b_j}).
\label{eq-aqm-partial-trace}
$$

これを :term[縮約密度行列]{en="reduced density matrix" yomi="しゅくやくみつどぎょうれつ"} という。$A$ の基底も入れてトレースを展開すると $\tr(\rho_A X)=\tr(\rho_{AB}(X\otimes I_B))$ を得る。この関係で $\rho_A$ は一意に決まるので $B$ の基底に依存しない。正値性は各項が正値であることから、規格化は $X=I_A$ から従う。

純粋二体系で $\ket{\Psi}$ が積 $\ket a\otimes\ket b$ と書けない場合を :term[量子もつれ]{en="quantum entanglement" yomi="りょうしもつれ"} という。混合状態では積状態の密度演算子の確率混合で書けるものを分離可能と呼ぶ。全体が分離可能でも部分系の結果は統計的に相関し得る。局所の密度行列だけでは全体の分離可能性を判定できない。

:::example[Bell 状態の縮約]
$\ket{\Phi^+}=(\ket{00}+\ket{11})/\sqrt2$ では交差項の $B$ の内積が零なので $\rho_A=\rho_B=I/2$ である。$\langle\sigma_z\otimes\sigma_z\rangle=1$ だが各単独期待値は零であり、相関は局所密度に残らない。
:::

:::example[古典相関との比較]
$\rho_{\rm cl}=(\ketbra{00}{00}+\ketbra{11}{11})/2$ も同じ縮約と $zz$ 相関を持つ。$\sigma_x\otimes\sigma_x$ は $\ket{00}$ と $\ket{11}$ を交換するので $\ket{\Phi^+}$ では期待値 1、$\rho_{\rm cl}$ では 0 になる。
:::

### 部分転置で可分性を検査する

有限次元で $\rho_{AB}=\sum_i p_i\rho_{A,i}\otimes\rho_{B,i}$ と書ける状態を可分とする。基底を固定した部分転置 $\rho_{AB}^{T_B}$ は、行列単位としての $B$ 成分を転置したものである。可分なら $\rho_{AB}^{T_B}=\sum_i p_i\rho_{A,i}\otimes\rho_{B,i}^{\mathsf T}$ である。転置した正行列は固有値を保ち、正行列のテンソル積と正の和も正である。従って部分転置に負の固有値が出れば、可分ではない。ここで必要なのは十分な検出法としての向きである。任意次元で正の部分転置から可分性までは従わない。

:::example[二準位の混合の比較]
$\rho=p\ketbra{\Phi^+}{\Phi^+}+(1-p)I_4/4$、$\ket{\Phi^+}=(\ket{00}+\ket{11})/\sqrt2$ とする。純粋部分の転置は交換演算子の半分である。例えば $\ket{00}\bra{11}$ は $\ket{01}\bra{10}$ に移る。交換演算子は対称空間の三状態で $+1$、反対称空間の一状態で $-1$ なので、固有値は $(1+p)/4$ が三つと $(1-3p)/4$ が一つである。従って $p>1/3$ ではもつれている。$p=1/2$ で負の値は $-1/8$ である。どの $p$ でも部分系は $I_2/2$ なので、局所の純度ではこの検査を置き換えられない。

$0\le p\le1/3$ では可分な分解を作れる。$i=x,z$ では $\sigma_i$ の同符号の固有状態を両側に並べた二つの積の等混合が $(I_4+\sigma_i\otimes\sigma_i)/4$ である。$i=y$ では逆符号の積の等混合が $(I_4-Y\otimes Y)/4$ になる。この三つの密度の等混合が $p=1/3$ の密度である。$0\le p<1/3$ ではその密度を確率 $3p$、積状態 $I_4/4$ を確率 $1-3p$ で混ぜればよい。
:::

### 部分系への作用を添字の縮約で確かめる

積基底で $(\rho_A)_{ij}=\sum_\mu(\rho_{AB})_{i\mu,j\mu}$ とする。$B$ の基底 $\{\ket\nu\}$ を $\ket{\tilde\mu}=\sum_\nu U_{\nu\mu}\ket\nu$ へ替えると、和に $\sum_\mu\overline{U_{\nu\mu}}U_{\nu'\mu}=\delta_{\nu\nu'}$ が入る。従って $\rho_A$ は変わらない。$B$ にチャネル $X\mapsto\sum_j B_jXB_j^\dagger$ が作用する場合、任意の $A$ の観測量の平均は
$$
\sum_j\tr[(A\otimes I)(I\otimes B_j)\rho(I\otimes B_j^\dagger)]
=\tr[(A\otimes \sum_j B_j^\dagger B_j)\rho]
=\tr[(A\otimes I)\rho]
$$
になる。したがって跡保存の局所変換は他方の密度を変えない。ただし選択された結果一つだけを残すと、$\sum B_j^\dagger B_j$ が恒等作用になる条件を使えず、他方の条件付き密度は変わり得る。

:::example[位相を消した後の三つの測定]
$\ket{\Phi^+}$ の $B$ に確率 $p$ の位相反転を掛ける。出力は $(1-p)\ketbra{\Phi^+}{\Phi^+}+p\ketbra{\Phi^-}{\Phi^-}$ であり、局所はともに $I/2$ のままである。$Z\otimes Z$ はどちらの Bell 状態でも $+1$、$X\otimes X$ はそれぞれ $+1,-1$ なので、平均は順に $1,1-2p$ である。$p=1/2$ では $X$ の相関が消えるが、$Z$ の相関は変わらない。この時点の状態は $\ket{00},\ket{11}$ の等混合である。別々の測定の平均だけでは相関の変化が見えない。
:::

### 可分な混合でも相関は現れる

二変数の共分散を $\langle A\otimes B\rangle-\langle A\otimes I\rangle\langle I\otimes B\rangle$ と定義する。積状態では跡の積の性質から零だが、可分な混合では $\rho_i=\rho_{A,i}\otimes\rho_{B,i}$ ごとの平均が同じ混合添字に依存しているため、共分散は一般に零でない。逆に全ての局所演算子 $A,B$ について共分散が零なら、$\rho_{AB}$ と $\rho_A\otimes\rho_B$ の行列要素は全て一致する。実際、部分空間ごとの行列単位のテンソル積が全行列を張り、その平均を比較すればよい。従って相関を一つ見つけることで積状態は排除できるが、もつれを判定するには混合も考える必要がある。

### Bell 基底で相関を一つずつ読む

$\ket{\Phi^\pm}=(\ket{00}\pm\ket{11})/\sqrt2$、$\ket{\Psi^\pm}=(\ket{01}\pm\ket{10})/\sqrt2$ を取る。$Z\otimes Z$ は $\Phi$ に $+1$、$\Psi$ に $-1$ を掛け、$X\otimes X$ は上付きの符号を固有値に持つ。$Y\otimes Y=-(X\otimes X)(Z\otimes Z)$ である。これを使って Bell 基底で対角の状態を $\rho=(I_4+t_x X\otimes X+t_yY\otimes Y+t_z Z\otimes Z)/4$ と書ける。各固有値は $(1+s_x t_x+s_y t_y+s_z t_z)/4$ で、$s_i$ は対応する基底ベクトルでの $\sigma_i\otimes\sigma_i$ の固有値である。例えば $\Phi^+$ では $(s_x,s_y,s_z)=(1,-1,1)$ となる。すべての Bell 射影の局所縮約が $I/2$ なので、三係数は二体系の相関だけにある。部分転置は $Y^\mathsf T=-Y$ により $t_y$ の符号を替えるため、同じ固有値の計算で正値性を確かめられる。

:::example[相関の三方向で可分性の候補を検査する]
$t_x=t_z=1/2,t_y=-1/2$ とすれば $\Phi^+$ の重みは $5/8$、残り三つは $1/8$ である。部分転置では $\Psi^-$ の符号が $(-1,-1,-1)$ なので、そこに $(1-3/2)/4=-1/8$ が出る。対角の測定だけを見れば $ZZ$ の平均は $1/2$ であり、古典的な混合でも作れる範囲である。複数の軸の相関を合わせると、一つの状態に必要な可分性の条件に反する。
:::

## Schmidt 分解と純粋化

二体系の係数行列を $C$ とすると、$A$ 側の密度行列は $CC^\dagger$ になる。これを対角化すれば両側に直交する基底を選べる。

:::theorem[Schmidt 分解]{#thm-aqm-schmidt}
任意の規格化純粋状態 $\Psi\in\mathcal H_A\otimes\mathcal H_B$ は

$$
\ket\Psi=\sum_{\alpha=1}^r\sqrt{\lambda_\alpha}\ket{u_\alpha}\otimes\ket{v_\alpha},
\qquad \lambda_\alpha>0,\quad\sum_\alpha\lambda_\alpha=1
$$

と表せる。両側のベクトルは正規直交し、$\lambda_\alpha$ は両縮約の非零固有値である。この :term[Schmidt 分解]{en="Schmidt decomposition" yomi="しゅみっとぶんかい"} が一項のとき、かつそのときに限り純粋状態は積である。
:::

:::proof
$\rho_A$ を対角化して $\rho_A u_\alpha=\lambda_\alpha u_\alpha$ とする。$w_\alpha=(\bra{u_\alpha}\otimes I)\ket\Psi$ と置くと $\braket{w_\beta}{w_\alpha}=\mel{u_\alpha}{\rho_A}{u_\beta}$ である。よって零固有値の $w_\alpha$ は零で、それ以外を $v_\alpha=w_\alpha/\sqrt{\lambda_\alpha}$ と規格化すれば正規直交する。$A$ の基底展開から式を得る。片側を部分トレースすれば反対側にも同じ固有値が出る。積なら縮約が階数 1 であり、階数 1 なら表示が一項なので積である。
:::

これを逆向きに用いると、$\rho_A=\sum_i\lambda_i\ketbra{u_i}{u_i}$ に補助空間の直交基底 $\ket i_B$ を割り当てた $\ket\Psi=\sum_i\sqrt{\lambda_i}\ket{u_i}\ket i_B$ が得られる。この :term[純粋化]{en="purification" yomi="じゅんすいか"} の補助空間には少なくとも $\rank\rho_A$ 次元が必要である。二つの純粋化は $A$ の固有基底を同じに取れば、$B$ の正規直交な係数ベクトルの組だけが異なる。同じ次元の補助空間なら両組を基底まで延長して、片方からもう片方へのユニタリ変換を作れる。

:::example[非等重みのもつれ]
$\sqrt{3/4}\ket{00}+\sqrt{1/4}\ket{11}$ の Schmidt 係数は $\sqrt3/2,1/2$ で、縮約純度は $5/8$ である。全体の純度は 1 であるため、全体と部分の混合を区別できる。
:::

:::example[熱状態の純粋化]
$\rho_\beta=e^{-\beta H}/Z$ の純粋化は $\ket\Psi=Z^{-1/2}\sum_i e^{-\beta E_i/2}\ket{i,i}$ である。$B$ の内積で異なる $i$ が消え、確率に $e^{-\beta E_i}$ が残る。係数の指数は確率の指数の半分である。
:::

### 純粋化の自由度と最小の補助次元

$\rho_A$ の正の固有値を $\lambda_i$ とし、固有ベクトル $\ket i$ を固定する。任意の純粋化 $\ket\Psi$ の $\ket i$ 成分を $\sqrt{\lambda_i}\ket{\alpha_i}$ と書けば、部分トレースから $\braket{\alpha_j}{\alpha_i}=\delta_{ij}$ になる。零固有値に対応する成分はノルムが零なので消える。従って補助系は少なくとも $\rank\rho_A$ 次元要り、その次元で標準的な純粋化 $\sum_i\sqrt{\lambda_i}\ket{i,i}$ が作れる。同じ次元の補助系で二種類の純粋化を取ると、対応する直交列を互いに移すユニタリを正規直交基底へ延長できる。したがって両者の差は補助系だけのユニタリである。

:::example[同じ混合を異なる測定で準備する]
$\rho_A=\diag(3/4,1/4)$ の純粋化 $\sqrt3\ket{00}/2+\ket{11}/2$ を作る。補助系を $z$ 基底で測ると $\ket0,\ket1$ が確率 $3/4,1/4$ で残る。$x$ 基底で測ると結果は等確率で、条件付き状態は $(\sqrt3\ket0\pm\ket1)/2$ である。後者の二射影を等しく混ぜると非対角が相殺し、同じ $\rho_A$ に戻る。補助の測定基底は分解を変えるが、結果を捨てた密度を変えない。
:::

### 任意の純粋状態分解を補助系で実現する

$\rho_A=\sum_j p_j\ketbra{\psi_j}{\psi_j}$ に対し、$p_j>0$ の項を並べた $\sum_j\sqrt{p_j}\ket{\psi_j}\otimes\ket j$ は $\rho_A$ の純粋化である。補助空間の次元をこの項数以上にしておけば、標準純粋化を補助側のユニタリだけでこのベクトルへ移せる。そこで補助を $\ket j$ 基底で測ると、結果 $j$ は確率 $p_j$ で出て $A$ には $\ket{\psi_j}$ が残る。項数が $\rank\rho_A$ より大きくても、補助系に空の次元を加えれば同じ議論が使える。

### 係数行列に特異値分解を施す

係数を $C_{ij}$ として $\ket\Psi=\sum_{ij}C_{ij}\ket i_A\ket j_B$ とする。$\rho_A=CC^\dagger$、$\rho_B=\transp C\bar C$ である。$CC^\dagger$ の正の固有値を $s_k^2$、固有ベクトルを $u_k$ とし、$v_k=C^\dagger u_k/s_k$ を作ると $v_k^\dagger v_l=\delta_{kl}$ である。逆に $C v_k=s_k u_k$ であり、零固有空間に対する $C^\dagger$ は零である。従って $C=\sum_k s_k u_k v_k^\dagger$、ket では $B$ の成分に $\bar v_k$ が使われる。これで :ref[thm-aqm-schmidt] の係数を基底が複素の場合にも求められる。$\sum s_k^2=1$ で、$\rho_B$ の正の固有値も同じである。

:::example[複素の係数で共役を確認する]
$\ket\Psi=(\ket{00}+i\ket{01}+2\ket{10})/\sqrt6$ とする。$\rho_A$ は $\begin{pmatrix}2&2\\2&4\end{pmatrix}/6$、$\rho_B$ は $\begin{pmatrix}5&-i\\i&1\end{pmatrix}/6$ である。$B$ の $01$ 成分には $\sum_i C_{i0}\overline{C_{i1}}=-i/6$ が入る。両方のトレースと行列式は $1,1/9$ であり、固有値は $(3\pm\sqrt5)/6$ である。純度は $1-2/9=7/9$ で、全体が純粋でも各部分は混合である。$\rho_B$ で $C^\dagger C$ をそのまま使うと固有値は一致しても、$y$ 軸で読む期待値の符号を誤る。
:::

### 最大重なりから Schmidt 係数を読む

規格化した積状態 $\ket{a,b}$ と $\ket\Psi=\sum_k s_k\ket{u_k,w_k}$ との重なりの絶対値は
$$
\left|\sum_k s_k\braket a{u_k}\braket b{w_k}\right|
\le s_{\max}\sqrt{\sum_k|\braket a{u_k}|^2\sum_k|\braket b{w_k}|^2}
\le s_{\max}
$$
である。最初は Cauchy–Schwarz、次は二つの直交列への射影のノルムを使った。等号は最大係数に対応する一組の積で達する。したがって最大の Schmidt 重み $s_{\max}^2$ が純粋状態を積状態で近似できる最大の確率になる。純粋な二体系でこの値が 1 なら階数 1、それより小さければもつれている。

:::example[同じ純度で積近似を評価する]
$\sqrt{3/4}\ket{00}+\sqrt{1/4}\ket{11}$ では積状態の最大重なり二乗は $3/4$ で、$\ket{00}$ で達する。$\ket+\ket+$ との重なり二乗は $(2+\sqrt3)/8\simeq0.467$ であり、固有成分を均等に重ねても大きくならない。局所純度は $5/8$ で、局所状態の混合の程度と、積状態へ射影する最大確率は別の数である。
:::

### 純粋化の重なりを二つの密度から計算する

同じ $d$ 次元の補助系に非規格化の $\ket\Omega=\sum_i\ket{i,i}$ を置く。密度 $\rho$ の標準純粋化を $(\sqrt\rho\otimes I)\ket\Omega$ とすればノルムは $\tr\rho=1$ である。$\sigma$ の純粋化を補助側のユニタリ $U$ で動かしたとき、重なりは $\tr(\sqrt\rho\sqrt\sigma\,\transp U)$ となる。これは $\bra\Omega A\otimes B\ket\Omega=\sum_{ij}A_{ij}B_{ij}=\tr(A\transp B)$ から分かる。$A=\sqrt\rho\sqrt\sigma$ の特異値分解を $V D W^\dagger$ とすると、任意のユニタリ $B$ について
$$
|\tr(AB)|=|\tr(D W^\dagger B V)|\le\sum_i D_{ii}
$$
である。$W^\dagger B V$ の対角成分の絶対値は 1 以下であり、恒等行列になる $B$ を選べば等号に達する。前節の純粋化の自由度により、最大重なりは $\tr\sqrt{\sqrt\rho\,\sigma\sqrt\rho}$ である。ここでの :term[忠実度]{en="fidelity" yomi="ちゅうじつど"} はその二乗とする。純粋な $\rho=\ketbra\psi\psi$ では $\sqrt\rho\,\sigma\sqrt\rho=\mel\psi\sigma\psi\,\rho$ なので忠実度は射影の確率に戻る。

:::example[古典的な対角分布の比較]
$\rho=\diag(3/4,1/4)$、$\sigma=\diag(1/4,3/4)$ に対し $\sqrt\rho\sqrt\sigma=\diag(\sqrt3/4,\sqrt3/4)$ である。最大重なりは $\sqrt3/2$、忠実度は $3/4$ になる。補助の $\sigma_x$ で純粋化を入れ替えるとこの重なりは零になる。二つの密度が同じなら最大は 1 になるが、選んだ純粋化の組の重なりが常に 1 であるわけではない。
:::

## 量子チャネルと Kraus 表現

閉じた系のユニタリ発展は、環境を観測しない場合には部分系だけのユニタリ発展に帰着しない。例えば環境に状態の情報が移れば、局所の密度行列の固有値も変わる。部分系から部分系への線形写像でその変化を表す。

:::definition[量子チャネル]{#def-aqm-channel}
有限次元の :term[量子チャネル]{en="quantum channel" yomi="りょうしちゃねる"} はトレース保存で完全正値な線形写像である。:term[完全正値性]{en="complete positivity" yomi="かんぜんせいちせい"} とは任意の補助次元 $d$ で $\mathcal E\otimes\operatorname{id}_d$ が正値演算子を正値演算子に写すことをいう。
:::

:::theorem[Kraus 表現]{#thm-aqm-kraus}
$\mathcal E$ が量子チャネルであることと、ある演算子の組 $(K_j)$ により $\mathcal E(X)=\sum_j K_j X K_j^\dagger$、$\sum_jK_j^\dagger K_j=I$ と書けることは同値である。
:::

:::proof
表示があれば補助系で $K_j\otimes I$ を使うので完全正値であり、巡回性でトレースを保つ。逆には入力基底 $\ket i$ と非規格化ベクトル $\ket\Omega=\sum_i\ket{i,i}$ を用い、$C=(\mathcal E\otimes I)(\ketbra{\Omega}{\Omega})$ を作る。完全正値性により $C\ge0$ なので $C=\sum_j\ketbra{v_j}{v_j}$ と、固有値の平方根を係数に吸収して書ける。$\ket{v_j}=\sum_i(K_j\ket i)\otimes\ket i$ と読めば、固定した入力添字 $i,l$ のブロックは $\mathcal E(\ketbra i l)=\sum_j K_j\ketbra i l K_j^\dagger$ になる。行列単位を線形に足して任意の $X$ に拡張する。トレース保存から $\tr[X(\sum K_j^\dagger K_j-I)]=0$ が全行列で成り立つため、括弧内は零である。
:::

$V=\sum_j K_j\otimes\ket j$ は $V^\dagger V=I$ の等長写像である。$\tr_B(V X V^\dagger)$ が Kraus 和となる。入出力が同じ次元なら $V$ の像の正規直交基底を拡張して $V\psi=U(\psi\otimes\ket0)$ となるユニタリを作れる。逆に、初期環境が純粋で相関がないなら $K_j=(I\otimes\bra j)U(I\otimes\ket0)$ である。初期環境が混合の場合も純粋化して同じ構成を使える。初期に系と環境が相関している場合は、系の初期密度だけから全入力に共通のチャネルを構成できるとは限らない。

:::example[完全正値でない転置]
転置 $X\mapsto X^\mathsf T$ は単体系の正値性を保つ。ところが $\ket{\Phi^+}$ の片側だけを転置すると、得られる演算子は交換演算子 $F_{AB}$ の半分である。$\ket{01}-\ket{10}$ に対する固有値は $-1/2$ なので完全正値ではない。
:::

:::example[位相反転のランダム適用]
$K_0=\sqrt{1-p}I$、$K_1=\sqrt p\,\sigma_z$ なら正規化条件を満たす。対角は保たれ、非対角には $1-2p$ が掛かる。$p=1/2$ では位相コヒーレンスが完全に消える。
:::

### 最小 Kraus 数と転置の反例

Choi 行列は $\sum_j\ketbra{\operatorname{vec}K_j}{\operatorname{vec}K_j}$ である。従ってその階数は Kraus 演算子の本数以下である。逆に Choi 行列の正の固有値だけを固有ベクトルの外積として展開すれば、その階数と同じ本数で表せる。これは有限次元での最小本数を与える。等しい本数を揃えた二つの分解を列行列 $A,B$ にまとめると $AA^\dagger=BB^\dagger$ である。その正の固有空間で逆平方根を左に掛ければ両行列の行が正規直交になり、列空間のユニタリによって移り合う。従って Kraus の添字のユニタリ混合だけで同じチャネルの表現が得られる。

:::example[転置はなぜ量子チャネルにならないか]
単独の密度の転置は固有値を保つので正で跡保存である。だが $\ket{\Phi^+}$ に片側だけ転置を掛けると交換演算子の半分になる。反対称ベクトルの期待値は $-1/2$ であるから、二次元の補助だけで完全正値性に反する。Choi の規約で非規格化の $\ket\Omega=\ket{00}+\ket{11}$ を使うと固有値はこれを二倍し、負の値は $-1$ である。
:::

### チャネルをユニタリと環境に分ける

Kraus 表現を $\mathcal E(X)=\sum_j K_j X K_j^\dagger$ とする。補助基底 $\ket j_E$ を用い $V\ket\psi=\sum_j K_j\ket\psi\otimes\ket j_E$ と置く。$V^\dagger V=\sum_j K_j^\dagger K_j=I$ なので $V$ は等長である。環境を跡消去すると二つの環境添字の和から $j=j'$ だけが残り、元のチャネルになる。入出力の系が同じ有限次元なら $\ket\psi\ket0_E\mapsto V\ket\psi$ の等長作用を、直交補の基底を補うことで全体のユニタリに延長できる。別の出力次元を持つ場合にも等長作用の表示までは同様である。

:::example[振幅減衰を一つの環境量子に写す]
$0\le\eta\le1$ とし、$\ket{g,0}\mapsto\ket{g,0}$、$\ket{e,0}\mapsto\sqrt\eta\ket{e,0}+\sqrt{1-\eta}\ket{g,1}$ と取る。二つの像は直交しノルム 1 なので、残りの空間への回転を補える。$\eta=1/4$、入力 $(\ket g+\ket e)/\sqrt2$ では励起占有 $1/8$、基底 $7/8$、両者の非対角係数 $1/4$ になる。環境の 1 を見た試行の確率は $3/8$ であり、その系は $\ket g$ に固定される。環境の 0 の場合は初期の相対位相が振幅比 $1:1/2$ で残る。環境を跡消去と観測のどちらで扱うかで質問が変わる。
:::

### トレース距離は跡保存の完全正値発展で増えない

双対 $\mathcal E^\dagger(A)=\sum_j K_j^\dagger A K_j$ は正で単位を保つ。$0\le E\le I$ なら $\mathcal E^\dagger(E)\ge0$ かつ $I-\mathcal E^\dagger(E)=\mathcal E^\dagger(I-E)\ge0$ である。測定の確率差 $\tr(E\mathcal E(D))$ は、入力で効果 $\mathcal E^\dagger(E)$ を測るのと等しい。出力の効果の全体に対する最大化は、入力の全効果に対する最大化を超えない。従ってトレース距離は縮小する。補助を含めた距離にも同じ証明を適用できるのは $\mathcal E\otimes I$ がチャネルだからである。

:::example[振幅減衰で軸によって距離の減少が違う]
$\ket e,\ket g$ の密度の差は $\ketbra e e-\ketbra g g$ なので、減衰後は $\eta$ 倍になり距離は $\eta$ である。$\ket{\pm}=(\ket g\pm\ket e)/\sqrt2$ では対角が元から同じであり、差は非対角だけなので距離は $\sqrt\eta$ である。$\eta=1/4$ なら順に $1/4,1/2$ となる。Bloch 球の中心が基底側へ移ることは、二状態の差をとると相殺される。
:::

### Bloch 球を縮める Pauli チャネルの条件

$\mathcal E(\rho)=p_0\rho+p_x X\rho X+p_yY\rho Y+p_zZ\rho Z$、$\sum_a p_a=1$ とする。$X,Y,Z$ は Pauli 行列である。共役で同じ軸の Pauli 行列は保たれ、残り二軸は符号が替わる。従って $\lambda_x=p_0+p_x-p_y-p_z$ などが各軸の倍率である。逆に $\mathcal E(I)=I$、$\mathcal E(\sigma_i)=\lambda_i\sigma_i$ の線形写像には、この四式を解いた $p_a$ がある。Choi 行列では $\operatorname{vec}(I),\operatorname{vec}(X),\operatorname{vec}(Y),\operatorname{vec}(Z)$ が直交してノルム二乗 2 なので、固有値は $2p_a$ である。したがってこの写像が完全正値であるための必要十分条件は、各 $p_a$ が非負である。$\max_i|\lambda_i|\le1$ は一系の Bloch 球を保つ条件だが、Choi の条件より弱い。

:::example[等方的な逆向きの写像]
全軸の倍率を同じ $\lambda$ とすると $p_0=(1+3\lambda)/4$、他は $(1-\lambda)/4$ である。従って完全正値は $-1/3\le\lambda\le1$ である。$\lambda=-1$ は一系で Bloch 球の各点を対蹠点へ移すが、$p_0=-1/2$ なので量子チャネルではない。$\lambda=-1/3$ なら三つの Pauli 共役を等確率で掛ければ実現できる。純粋な入力と反対の射影への成功確率は $(1-\lambda)/2=2/3$ である。
:::

## Born–Markov 近似と Lindblad 方程式

系 $S$ と浴 $B$ の Hamiltonian を $H_S+H_B+H_I$ とし、$H_I=\sum_\alpha A_\alpha\otimes B_\alpha$ を Hermite な因子で表す。浴は $H_B$ と可換な $\rho_B$ にあり $\tr(\rho_B B_\alpha)=0$ とする。平均が非零なら一体の項に吸収できる。相互作用表示の式 $\dot{\widetilde\rho}_{SB}=-(i/\hbar)[\widetilde H_I,\widetilde\rho_{SB}]$ を一度積分して元の右辺に代入する。初期積状態で一次項は跡を取ると零で、厳密に二重交換子の記憶積分が残る。

弱結合で浴の変化を無視し、積分中の全密度を $\widetilde\rho_S(t-\tau)\otimes\rho_B$ に替えるのが Born 近似である。浴の相関が $\tau_B$ で減衰し、系の緩和が $\tau_R\gg\tau_B$ なら $\widetilde\rho_S(t-\tau)$ を $\widetilde\rho_S(t)$ に替え、$t\gg\tau_B$ の時間帯で上端を無限に延ばす。これが Markov 近似であり、

$$
\dot{\widetilde\rho}_S(t)=-\frac1{\hbar^2}\int_0^\infty\tr_B
[\widetilde H_I(t),[\widetilde H_I(t-\tau),\widetilde\rho_S(t)\otimes\rho_B]]
\,\dd\tau
$$

となる。この段階だけでは完全正値性が保証されない。

### 周波数ごとの散逸項

$A_\alpha(\omega)=\sum_{E'-E=\hbar\omega}P_E A_\alpha P_{E'}$ と置けば $\widetilde A_\alpha(t)=\sum_\omega e^{-i\omega t}A_\alpha(\omega)$ である。二重交換子を展開すると、$A_\beta(\omega)\rho A_\alpha(\omega')^\dagger$ と、$\rho$ の同じ側に二つの $A$ を掛ける項が出る。異なる周波数の項には $e^{i(\omega'-\omega)t}$ が掛かる。緩和時間に比べ短くかつ異周波数の逆差より長い窓で平均すると、それらは消える。この永年近似では同じ周波数を持つ縮退遷移は残す。残る片側 Fourier 積分を

$$
G_{\alpha\beta}(\omega)=\hbar^{-2}\int_0^\infty
e^{i\omega\tau}\tr[\rho_B B_\alpha(\tau)B_\beta(0)]\,\dd\tau
$$

と置く。異なる Bohr 周波数の逆差が緩和時間に近いと、この平均の窓は取れない。以下では周波数差を分解できる弱結合の領域に限る。同一周波数の四項は

$$
\sum_{\alpha\beta}\big[
G_{\alpha\beta}(A_\beta\rho A_\alpha^\dagger-A_\alpha^\dagger A_\beta\rho)
+G_{\beta\alpha}^*(A_\beta\rho A_\alpha^\dagger-\rho A_\alpha^\dagger A_\beta)
\big]
$$

になる。よって Hermite 部 $\gamma=G+G^\dagger$ が散逸、$\hbar(G-G^\dagger)/(2i)$ が $A_\alpha^\dagger A_\beta$ の係数となって Hamiltonian をずらす。定常性から $\gamma$ は相関の両側 Fourier 変換である。浴のエネルギー基底を挟むと、各デルタ関数 $\delta(\omega-(E_b-E_a)/\hbar)$ に $p_a(B_\alpha)_{ab}(B_\beta)_{ba}$ が掛かる。これは各 $(a,b)$ について一つのベクトルの外積を正の重みで足した形だから、$\gamma$ は正半定値である。連続スペクトルでは和を測度に替える。

$\gamma$ を対角化し非負固有値をジャンプ演算子 $L_j$ に平方根で吸収すると

$$
\dot\rho=-\frac i\hbar[H_{\mathrm{eff}},\rho]+\sum_j\mathcal D[L_j]\rho,\qquad
\mathcal D[L]\rho=L\rho L^\dagger-\frac12\{L^\dagger L,\rho\}
\label{eq-aqm-lindblad}
$$

になる。これが :term[Lindblad 方程式]{en="Lindblad equation" yomi="りんどぶらっどほうていしき"} であり、$L_j$ の単位は $\mathrm{s}^{-1/2}$ である。短時間 $\Delta t$ で $K_j=\sqrt{\Delta t}L_j$、$K_0=I-\Delta t(iH_{\rm eff}/\hbar+\sum L_j^\dagger L_j/2)$ とすると Kraus 和は $I+\Delta t\mathcal L+\order{\Delta t^2}$ の時間発展を与える。$M=\sum_j K_j^\dagger K_j=I+\order{\Delta t^2}$ は短時間で正則なので、各 $K_j$ の右から $M^{-1/2}$ を掛ければ、一次を保ったまま厳密にトレース保存になる。このチャネルを $t/\Delta t$ 回合成し極限を取ると、有限次元でノルム収束して $e^{t\mathcal L}$ は完全正値かつ跡保存となる。

:::example[跳躍なしの確率]
短時間の $j$ 番目の跳躍の確率は $\Delta t\tr(L_j^\dagger L_j\rho)$ である。残りは $K_0$ の確率で、一次の減少はそれらの総和に等しい。$K_0$ に $L^\dagger L$ の項が必要な理由である。
:::

:::example[自発的な位相緩和と両側 Fourier 成分]
浴相関を実かつ偶の $C(\tau)=C_0e^{-|\tau|/\tau_B}$ と近似すると、$\gamma(\omega)=2C_0\tau_B/[\hbar^2(1+\omega^2\tau_B^2)]$ になる。正負の $\tau$ の寄与は共役であり、$\int_0^\infty e^{-(1/\tau_B-i\omega)\tau}\dd\tau$ の実部から得られる。
:::

### 量子力学的半群の生成子

有限次元で時間並進不変の発展が $\mathcal E_{t+s}=\mathcal E_t\mathcal E_s$、$\mathcal E_0=\operatorname{id}$ を満たし、各 $\mathcal E_t$ がチャネルで $t=0$ で微分可能とする。このとき生成子 $\mathcal L=\dot{\mathcal E}_0$ は :eqref[eq-aqm-lindblad] の形で書ける。この逆向きを示せば、特定の浴の近似から得た式が、有限次元の完全正値な半群に共通する形であると分かる。

行列空間で Hilbert–Schmidt 正規直交基底 $F_0=I/\sqrt d,F_a$ を選ぶ。$a>0$ のとき $\tr F_a=0$ である。任意の線形写像の Choi 行列を $\operatorname{vec}(F_a)=\sum_i (F_a\ket i)\otimes\ket i$ の基底で展開できる。$\mathcal L$ の Choi 行列の成分を $c_{ab}$ とすると、:ref[thm-aqm-kraus] のブロックの対応と同じ計算により $\mathcal L(X)=\sum_{ab}c_{ab}F_aXF_b^\dagger$ である。$\mathcal E_t$ は Hermite 性を保つため、生成子の Choi 行列も Hermite である。

恒等写像の Choi 行列は $\ketbra{\Omega}{\Omega}=d\ketbra{\operatorname{vec}(F_0)}{\operatorname{vec}(F_0)}$ である。これに直交する任意の $w$ について、$\mathcal E_t$ の Choi 行列の期待値は $t\langle w|C_{\mathcal L}|w\rangle+o(t)$ である。$t\ge0$ で非負なので、$a,b>0$ のブロック $(c_{ab})$ は正半定値である。そこだけを対角化すると非負の固有値を取り込んだ $L_j$ により $\sum_j L_j X L_j^\dagger$ が得られる。他の項は片側に $F_0$ を持つので $C X+X C^\dagger$ と書ける。$a=b=0$ の実項は両側に半分ずつ分ければよい。トレース保存を微分すると、任意の $X$ に対して $\tr[(C+C^\dagger+\sum L_j^\dagger L_j)X]=0$ である。従って $C=-\sum L_j^\dagger L_j/2-iH_{\rm eff}/\hbar$ と分けられ、$H_{\rm eff}$ は Hermite である。これで生成子の形が示された。定数倍の恒等演算子を $H_{\rm eff}$ に足しても、生成子の交換子は変わらない。

### 熱浴での詳細釣り合い

$\rho_B=e^{-\beta H_B}/Z_B$ なら、浴エネルギーの基底で正周波数の係数に現れる条件は $E_b-E_a=\hbar\omega$ である。添字 $a,b$ を交換して負周波数を作ると、Boltzmann 重みが $p_b=e^{-\beta\hbar\omega}p_a$ に替わり、二つの浴演算子の行列要素の添字も逆になる。従って $\gamma_{\alpha\beta}(-\omega)=e^{-\beta\hbar\omega}\gamma_{\beta\alpha}(\omega)$ である。例えば二準位系が一つの熱浴とエネルギーを交換するなら $\Gamma_\uparrow/\Gamma_\downarrow=e^{-\beta\hbar\omega}$ となる。$\beta\hbar\omega=\log3$ では上向き率は下向き率の $1/3$、定常占有は $1/(1+3)=1/4$ になる。一般の系で全ての準位が浴と遷移でつながっているとは限らず、定常 Gibbs 状態の一意性は別の条件である。

### 単位を保つチャネルでのエントロピー

$\mathcal E$ がトレース保存に加えて $\mathcal E(I)=I$ を満たすとする。入力と出力の密度の固有基底を $u_j,v_i$ とし、固有値をそれぞれ $\lambda_j,\mu_i$ と書くと $\mu_i=\sum_j T_{ij}\lambda_j$、$T_{ij}=\mel{v_i}{\mathcal E(\ketbra{u_j}{u_j})}{v_i}\ge0$ である。トレース保存から各列、単位を保つことから各行の和は 1 である。各行に $-x\log x$ の凹性を適用すると $-\mu_i\log\mu_i\ge-\sum_j T_{ij}\lambda_j\log\lambda_j$ となる。行について足し、各列の和を用いると $S(\mathcal E(\rho))\ge S(\rho)$ となる。例えば等確率の位相反転が $\ket{+x}$ に作用すると $I/2$ を生じ、エントロピーは零から $k_B\log2$ へ増える。振幅減衰は単位を保たず、混合状態を純粋な基底へ近づけることができる。

### 平衡分布と連結性

Gibbs 密度 $\rho_\beta=e^{-\beta H_S}/Z_S$ が定常であることは、縮退を含めて確かめられる。$A_\alpha(\omega)\rho_\beta=e^{-\beta\hbar\omega}\rho_\beta A_\alpha(\omega)$ であり、$A_\alpha(\omega)^\dagger=A_\alpha(-\omega)$ である。また $A_\alpha^\dagger(\omega)A_\beta(\omega)$ は $H_S$ と可換なので、散逸の反交換子はその積と $\rho_\beta$ の通常の積に帰着する。一方、跳躍項に $\rho_\beta$ を左へ移して周波数と二つの添字を反転すると、その係数は $e^{\beta\hbar\omega}\gamma_{\beta\alpha}(-\omega)=\gamma_{\alpha\beta}(\omega)$ となり、反交換子の項を打ち消す。Hamiltonian のずれも同じ周波数内の積で作るため、$[H_{\rm eff},\rho_\beta]=0$ である。

有限準位の占有が閉じた方程式 $\dot p_i=\sum_{j\ne i}(W_{ij}p_j-W_{ji}p_i)$ を満たすとする。全準位に正の $p_i^\ast$ があり $W_{ij}p_j^\ast=W_{ji}p_i^\ast$ とする。$q_i=p_i/p_i^\ast$、$c_{ij}=W_{ij}p_j^\ast=c_{ji}$ とおけば $p_i^\ast\dot q_i=\sum_{j\ne i}c_{ij}(q_j-q_i)$ となる。定常解に $q_i$ を掛けて全和を取ると $-\sum_{i<j}c_{ij}(q_i-q_j)^2=0$ である。従って遷移率が正の辺で結ばれた連結成分の中では $q$ が一定となる。全体が連結なら規格化により $q=1$ で、定常占有は一意である。これは対角占有の主張であり、縮退空間のコヒーレンスが残らないことは別に確かめる。

:::example[浴に結合しない準位を含む場合]
$0\rightleftarrows1$ の上向き率 1、下向き率 3 とし、別の準位 2 を孤立させる。定常状態の占有は $p_2=w$ を固定して $(p_0,p_1)=(3(1-w)/4,(1-w)/4)$ である。初期が全部準位 2 ならそのままであり、初期が全部準位 0 なら $(3/4,1/4,0)$ へ緩和する。詳細釣り合いを満たす一つの Gibbs 分布があっても、孤立した準位があれば全体の初期条件が消えない。
:::

### 純度の時間微分が示す散逸の向き

Lindblad 方程式で $\dv{\tr\rho^2}{t}=2\tr(\rho\dot\rho)$ であり、Hamiltonian の項は巡回性で零である。$L=L^\dagger$ の場合、散逸は $-\tfrac12[L,[L,\rho]]$ と書ける。従って寄与は $-\tr([L,\rho]^\dagger[L,\rho])$ となり、零以下である。これは $[L,\rho]$ が反 Hermite であることと、トレース内で交換子をもう一方へ移すことから従う。$L$ が非 Hermite な一般の過程では純度が増えることもある。例えば減衰では環境へエネルギーと情報を捨てて、どの初期状態でも最終的に基底へ近づき得る。

:::example[純度が一度下がってから増える]
初期が励起 $\ket e$ の零温度減衰なら、$\rho=\diag(\eta,1-\eta)$、$\eta=e^{-\Gamma t}$ である。純度は $\eta^2+(1-\eta)^2=1-2\eta(1-\eta)$ なので、初期 1 から $\eta=1/2$ で $1/2$ まで下がり、長時間で 1 に戻る。この最小は $t=\log2/\Gamma$ である。純度の増減のみで Markov 近似の成否を判定できない。
:::

### 永年近似で残す縮退の結合

浴との結合が同じ周波数を持つ二遷移 $A_1,A_2$ に掛かり、その率行列が $\gamma=\gamma_0\begin{pmatrix}1&\chi\\\bar\chi&1\end{pmatrix}$、$\gamma_0>0$ とする。正値性は $|\chi|\le1$ と同値であり、固有率は $\gamma_0(1\pm|\chi|)$ である。同じ周波数なので $A_1\rho A_2^\dagger$ などの交差項の相対回転はない。率を独立と仮定した $L_1,L_2$ で置き換えると、この固有率は再現されない。

:::example[同じ浴に見える二つの励起]
$\ket{e_1},\ket{e_2}$ が縮退し、一つの $\ket g$ に落ちる $A_j=\ket g\bra{e_j}$ とする。ここでは Hamiltonian が後述の対称・反対称の両成分を混ぜず、逆向きの励起もないとする。$\chi=1$ なら固有率は $2\gamma_0,0$、ジャンプは $\sqrt{2\gamma_0}\ket g\bra{e_+}$ だけである。ただし $\ket{e_\pm}=(\ket{e_1}\pm\ket{e_2})/\sqrt2$ である。従って $\ket{e_-}$ は線形結合の打ち消しで崩壊せず、$\ket{e_+}$ は率 $2\gamma_0$ で崩壊する。$\ket{e_1}$ を準備した場合は両成分が各半分であるため、励起全占有は $(1+e^{-2\gamma_0t})/2$ となる。励起準位が分裂して永年平均で交差項が落ちれば、同じ計算は使えなくなる。
:::

### 時間局所の率と時間区間の完全正値性

対角を固定して非対角を実の $g(t)\ne0$ 倍にする二準位チャネルを考える。初期からの写像がチャネルであるには $|g(t)|\le1$ が必要十分であり、$p=(1-g(t))/2$ の位相反転の確率で実現できる。任意の中間 $s$ から遷移を合成する場合、$\mathcal E_t=\mathcal V_{t,s}\mathcal E_s$ で $\mathcal E_s$ が可逆なため、$\mathcal V_{t,s}$ の倍率は $g(t)/g(s)$ と一意に決まる。従って $|g|$ が増える区間ではその中間写像は完全正値ではない。$g>0$ を微分できるなら $\gamma(t)=-\dot g/g$ が時間局所な減衰率であり、$\gamma(t)\ge0$ の条件は倍率が単調に下がる条件に一致する。初期から各時刻への完全正値性と、任意時刻間でチャネルを分けられる性質は別である。

:::example[二準位の環境との可逆な回転]
環境の $\ket0$ へ、系の $\ket0$ では恒等、$\ket1$ では $e^{-i\Omega t\sigma_y}$ を掛ける。初期と回した環境の重なりは $\cos\Omega t$ なので $g(t)=\cos\Omega t$ である。$t_1=\pi/(3\Omega)$ で $|g|=1/2$、$t_2=5\pi/(6\Omega)$ で $\sqrt3/2$ になる。各時点では場を跡消去した完全正値な状態だが、二つの時点を結ぶ非対角の倍率は $-\sqrt3$ になる。$\ket{+x},\ket{-x}$ の距離が増えるので、その区間全体を一つの系だけのチャネルで説明できない。
:::

## デコヒーレンス

環境が系の異なる成分を区別できると相対位相を局所測定で読み出せなくなる。例えば同じ初期環境が $\ket0,\ket1$ に応じて $\ket{e_0},\ket{e_1}$ に写れば、全体は $a\ket{0,e_0}+b\ket{1,e_1}$ である。:eqref[eq-aqm-partial-trace] により非対角成分は $a\bar b\braket{e_1}{e_0}$ となる。環境状態が直交に近いほど、局所で干渉を測りにくい。この過程を :term[デコヒーレンス]{en="decoherence" yomi="でこひーれんす"} という。他方、全体系を逆ユニタリで戻せる場合には、局所で消えた位相を回復できる。

二準位で $H=\hbar\omega\sigma_z/2$、$L=\sqrt{\gamma_\varphi/2}\,\sigma_z$ とする。$\mathcal D[L]\rho=(\gamma_\varphi/2)(\sigma_z\rho\sigma_z-\rho)$ は対角を変えず、非対角に $-\gamma_\varphi$ を掛ける。従って $\rho_{01}(t)=e^{-(i\omega+\gamma_\varphi)t}\rho_{01}(0)$ である。エネルギー基底で占有を変えず位相だけを失うので純粋位相緩和という。

:::example[環境の重なりを測る]
$a=b=1/\sqrt2$、$\braket{e_1}{e_0}=1/2$ なら $\rho_A=(I+\sigma_x/2)/2$ で、$\sigma_x$ の $+1$ 確率は $3/4$ である。重なりが 0 なら確率は $1/2$、重なりが 1 なら $1$ となる。
:::

:::example[位相緩和の純度]
初期が $\ket+$ で Hamiltonian を取り除いた表示なら $\rho(t)=(I+e^{-\gamma_\varphi t}\sigma_x)/2$ である。純度は $(1+e^{-2\gamma_\varphi t})/2$ となり、長時間で $1/2$ になる。
:::

### 静的なばらつきと反転パルス

同じ試行中は一定だが試行ごとに異なる $\delta$ を持つ $H_\delta=\hbar\delta\sigma_z/2$ を考える。各試行の非対角は $e^{-i\delta t}$ 倍になり、平均は $\delta$ の特性関数である。平均零、分散 $s^2$ の Gaussian なら $e^{-s^2t^2/2}$ である。時刻 $t/2$ に $\sigma_x$ を挟むと $\sigma_x H_\delta\sigma_x=-H_\delta$ なので、$e^{-iH_\delta t/(2\hbar)}\sigma_x e^{-iH_\delta t/(2\hbar)}=\sigma_x$ となり、全試行で同じ出力を得る。環境が時刻の前半と後半で変動すると、残る位相は二つの区間の積分の差になる。

:::example[指数減衰と Gaussian 減衰を同じ時刻で合わせる]
両モデルの包絡を $T$ で $e^{-1}$ に揃える。指数は $e^{-t/T}$、静的 Gaussian は $e^{-t^2/T^2}$ である。$t=2T$ ではそれぞれ $e^{-2},e^{-4}$ になる。静的模型の理想エコーは 1 に戻る。一方、$\sigma_z\rho\sigma_z-\rho$ の時間一定の散逸は $\sigma_x$ で共役しても符号が変わらず、同じ指数包絡を持つ。パルス後の包絡の違いで二つの機構を区別できる。
:::

### デコヒーレンスの重なりと環境での識別

系の二成分が等振幅で、環境の二つの純粋状態の重なりを $\xi=\braket{e_1}{e_0}$ とする。系の局所状態は対角 $1/2$、非対角 $\xi/2$ である。相対位相 $\phi$ の干渉計で片方の出力確率は $(1+\Re(e^{i\phi}\xi))/2$ となり、最大と最小の差は $|\xi|$ である。環境だけで二純粋状態を等事前確率で区別する場合、差の行列の固有値は $\pm\sqrt{1-|\xi|^2}$ である。これは $\ket{e_1}$ を $\ket{e_0}$ と直交補の二成分に分けると、トレース零で行列式が $-(1-|\xi|^2)$ になるためである。従って干渉の大きさの二乗と環境の識別距離の二乗の和は 1 になる。この等式には純粋な二分岐と等事前確率を用いている。

:::example[環境を強く見分けられる場合の干渉]
$|\xi|=3/5$ とすると出力確率は位相を走査して $1/5$ から $4/5$ の間を動く。環境での識別距離は $4/5$ なので、最適な一回の判定の成功確率は $9/10$ である。この場合でも環境の二状態は直交しておらず、干渉は完全には消えない。$|\xi|=0$ にすると環境の識別が完全になり、系の二出力はどの位相でも等確率になる。
:::

### 独立な環境との衝突から指数を得る

短い時刻 $\delta t$ ごとに新しい補助状態を用意し、系の非対角を $g_{\delta t}$ 倍にする同じ相互作用を掛ける。衝突した環境を毎回跡消去するため、固定 $t=N\delta t$ の倍率は $g_{\delta t}^N$ である。$g_{\delta t}=\cos\sqrt{2\gamma\delta t}=1-\gamma\delta t+\order{\delta t^2}$ と取れば、$N\log g_{\delta t}\to-\gamma t$ である。微小結合が $\sqrt{\delta t}$ に比例するのは一衝突の振幅に対して情報の低下がその二乗から出るためである。もし結合角を $\delta t$ に比例させると、この極限では指数の率が零になる。同じ環境を毎回使う場合には重なりの積を取れず、位相緩和の回復が起こり得る。

:::example[二回の衝突を独立に扱うかを数で区別する]
一衝突の角が $\pi/4$ なら重なりは $1/\sqrt2$ である。新しい補助を使う二回の倍率は $1/2$、同じ補助で制御回転を二回重ねた場合は $\cos(\pi/2)=0$ になる。四回では前者が $1/4$、後者は $\cos\pi=-1$ で絶対値が 1 に戻る。実験の $X$ の平均は等振幅の初期状態に対してこの $g$ になるので、初期位相を基準に符号も判別できる。
:::

## 例：減衰する二準位系と調和振動子

### 二準位の振幅減衰

励起 $\ket e$ から基底 $\ket g$ へ $L=\sqrt{\Gamma}\ketbra g e$ が作用し、$H=\hbar\omega\ketbra e e$ とする。$L^\dagger L=\Gamma\ketbra e e$ なので :eqref[eq-aqm-lindblad] は

$$
\dot\rho_{ee}=-\Gamma\rho_{ee},\qquad
\dot\rho_{eg}=-(i\omega+\Gamma/2)\rho_{eg},\qquad
\dot\rho_{gg}=\Gamma\rho_{ee}
$$

を与える。従って $\rho_{ee}(t)=e^{-\Gamma t}\rho_{ee}(0)$、$\rho_{eg}(t)=e^{-(i\omega+\Gamma/2)t}\rho_{eg}(0)$ となる。追加の純粋位相緩和があれば $T_1^{-1}=\Gamma$、$T_2^{-1}=\Gamma/2+\gamma_\varphi$ である。熱浴による上向き率 $\Gamma_\uparrow$ と下向き率 $\Gamma_\downarrow$ があれば $\dot\rho_{ee}=\Gamma_\uparrow(1-\rho_{ee})-\Gamma_\downarrow\rho_{ee}$ なので、緩和率は和、平衡励起は $\Gamma_\uparrow/(\Gamma_\uparrow+\Gamma_\downarrow)$ になる。

:::example[有限時間 Kraus 作用素]
$q=e^{-\Gamma t}$、Hamiltonian の位相を外した表示で $K_0=\ketbra g g+\sqrt q\ketbra e e$、$K_1=\sqrt{1-q}\ketbra g e$ とする。直接掛け合わせると占有は $q$ 倍、非対角は $\sqrt q$ 倍で、微分方程式と一致する。
:::

### 調和振動子の減衰

真空浴で $H=\hbar\omega a^\dagger a$、$L=\sqrt{\kappa}\,a$ とする。観測量の散逸の随伴は $\mathcal D^\dagger[L]X=L^\dagger X L-\{L^\dagger L,X\}/2$ である。これはトレースの巡回性で $\tr(X\mathcal D[L]\rho)=\tr((\mathcal D^\dagger[L]X)\rho)$ とすれば得られる。$a$ に対し $a^\dagger a^2-\{N,a\}/2=[N,a]/2=-a/2$、$N$ に対し $a^\dagger N a-N^2=-N$ である。よって $\langle a\rangle_t=e^{-(i\omega+\kappa/2)t}\langle a\rangle_0$、$\langle N\rangle_t=e^{-\kappa t}\langle N\rangle_0$ となる。無限次元での式は有限の数状態、または必要なモーメントが有限な状態に適用する。

対角要素の式は $\dot p_n=\kappa[(n+1)p_{n+1}-np_n]$ である。初期が数状態 $n_0$ のとき、母関数 $G(z,t)=\sum p_n z^n$ は $\partial_tG=\kappa(1-z)\partial_zG$ になる。$G(z,t)=(1+(z-1)e^{-\kappa t})^{n_0}$ が初期と方程式を満たすので、残存数は二項分布になる。熱浴なら $L_1=\sqrt{\kappa(\bar n+1)}a$、$L_2=\sqrt{\kappa\bar n}a^\dagger$ とする。上向き遷移は $\kappa\bar n(n+1)$、下向きは $\kappa(\bar n+1)n$ なので、隣接準位で流れを釣り合わせると $p_{n+1}/p_n=\bar n/(\bar n+1)$ となる。

:::example[二量子からの減衰]
$n_0=2$、$e^{-\kappa t}=1/2$ なら $p_0,p_1,p_2=(1/4,1/2,1/4)$ であり、平均は 1 である。決まった個数が失われるのではなく、放出の確率が時間とともに増す。
:::

### 真空振動子の非対角成分

透過率 $\eta=e^{-\kappa t}$ の二モード分岐で、環境の真空を加える。数状態を $a^\dagger\mapsto\sqrt\eta\,a^\dagger+\sqrt{1-\eta}\,b^\dagger$ で展開すると $\ket n$ から環境へ $r$ 個出る係数は $\sqrt{\binom nr}\eta^{(n-r)/2}(1-\eta)^{r/2}$ である。環境を跡消去すれば左右で同じ $r$ だけが残るので
$$
\mathcal E_\eta(\ket n\bra m)
=\sum_{r=0}^{\min(n,m)}\sqrt{\binom nr\binom mr}
  (1-\eta)^r\eta^{(n+m)/2-r}\ket{n-r}\bra{m-r}.
$$
$n=m$ では先の二項分布である。$1-\eta=\kappa t+\order{t^2}$ を入れると $r=0$ は $-\kappa(n+m)/2$、$r=1$ は $+\kappa\sqrt{nm}\ket{n-1}\bra{m-1}$ を与え、生成子と一致する。真空の基底を二回続けて加えると残存の振幅が $\sqrt{\eta_1\eta_2}$ になるため、同じ族の合成は透過率の積で表せる。

:::example[非対角と対角の減衰を別々に読む]
初期を $(\ket0+\ket2)/\sqrt2$ とし、$\eta=1/2$ まで減衰させる。$\ket2$ からの確率は $(1/4,1/2,1/4)$ なので、初期の真空も足した出力の対角は $(5/8,1/4,1/8)$ である。非対角の $\ket0\bra2$ は $r=0$ だけで $1/2$ 倍され、係数は $1/4$ になる。出力の平均数は $1/4+2/8=1/2$ で、初期の 1 から半分になったことを確認できる。
:::

### 熱浴で振動子の平均がどこへ緩和するか

有限のモーメントを持つ振動子で $L_1=\sqrt{\kappa(\bar n+1)}a$、$L_2=\sqrt{\kappa\bar n}a^\dagger$ とする。$\mathcal D^\dagger[a]N=-N$、$\mathcal D^\dagger[a^\dagger]N=N+1$ である。後者は $aNa^\dagger$ を数状態に作用させ、$(n+1)^2$ から $\{N+1,N\}/2$ の $n(n+1)$ を引いて得る。従って
$$
\frac{\dd}{\dd t}\langle N\rangle=-\kappa(\langle N\rangle-\bar n),
\qquad
\langle N\rangle_t=\bar n+(\langle N\rangle_0-\bar n)e^{-\kappa t}.
$$
同じ方法で $\mathcal D^\dagger[a^\dagger]a=a/2$ だから、量子を上げる率があっても $\langle a\rangle$ の減衰率は $\kappa/2$ のままである。振幅の平均が零に向かう一方、エネルギーは熱占有に向かう。

:::example[浴が振動子を加熱する初期値]
$\bar n=2$、初期真空とすると $t=\log2/\kappa$ で平均は 1 になる。同じ浴で初期平均 4 ならこの時点で 3 である。どちらも 2 に向かうが、瞬間の上向きと下向きの率は分布によって異なる。例えば真空から初期の上向き率は $2\kappa$、$\ket4$ の下向き率は $12\kappa$、上向き率は $10\kappa$ であり、平均の変化は導いた方程式に一致する。
:::

### 減衰中のコヒーレント入力と重ね合わせ

真空への線形な漏出の等長作用で、生成演算子の変換と指数表示を使うと
$$
V_\eta\ket{\alpha,0}=\ket{\sqrt\eta\,\alpha,\sqrt{1-\eta}\,\alpha}
$$
である。実際、変位指数に $a^\dagger\mapsto\sqrt\eta\,a^\dagger+\sqrt{1-\eta}\,b^\dagger$ を入れると、異なるモードの指数は可換なので二つの変位の積に分かれる。従って入力が一つのコヒーレント状態なら減衰しても出力は純粋である。一方、入力が二つのコヒーレント状態の重ね合わせなら、非対角の $\ket\alpha\bra\beta$ の出力には環境の重なり
$$
\exp\left[-\frac{1-\eta}{2}(|\alpha|^2+|\beta|^2)
 +(1-\eta)\bar\beta\alpha\right]
$$
が掛かる。絶対値は $e^{-(1-\eta)|\alpha-\beta|^2/2}$ である。純粋な一分岐の保ちやすさと、離れた二分岐の干渉が弱まる速さを同じ等長作用から計算できる。

:::example[少しの漏出で大きな振幅差が見える]
$\alpha=2,\beta=-2$ を実に取り、$\eta=3/4$ とする。二分岐の環境振幅は $1,-1$ で、重なりは $e^{-2}$ である。出力の各コヒーレント振幅は $\pm\sqrt3$ であり、各枝の平均数は 4 から 3 に下がるだけである。しかし非対角の係数の絶対値は約 $0.135$ 倍になる。状態の全ノルムを計算するとき、出力の二つのコヒーレント状態が直交するとは仮定しない。実際、それらの重なりは $e^{-6}$ であり、環境の $e^{-2}$ と掛けて初期の $e^{-8}$ を再現する。
:::

## 演習問題

::::exercise{level=basic}
$\rho=(I+0.6\sigma_z)/2$ の固有値と純度を求めよ。
:::solution
$\sigma_z$ の固有値 $\pm1$ を用いると $\lambda_\pm=0.8,0.2$ である。二乗和は $0.64+0.04=0.68$ であり、混合状態である。
:::
::::

::::exercise{level=basic}
$\ket\Psi=\sqrt{1/3}\ket{00}+\sqrt{2/3}\ket{10}$ の Schmidt 階数を求めよ。
:::solution
$B$ の因子が共通の $\ket0$ なので $(\sqrt{1/3}\ket0+\sqrt{2/3}\ket1)\otimes\ket0$ である。両因子は規格化され階数は 1 となる。
:::
::::

::::exercise{level=basic}
位相反転チャネルで $p=1/4$ のとき非対角成分の係数を求めよ。
:::solution
恒等作用の重みは $3/4$、$\sigma_z$ 共役は非対角の符号を替えて重み $1/4$ である。合計は $3/4-1/4=1/2$ である。
:::
::::

::::exercise{level=basic}
振幅減衰で $\Gamma t=\log 4$ のとき励起占有と非対角成分の大きさの割合を求めよ。
:::solution
占有は $e^{-\Gamma t}=1/4$ 倍、非対角の大きさは $e^{-\Gamma t/2}=1/2$ 倍である。エネルギーによる位相は絶対値に寄与しない。
:::
::::

::::exercise{level=advanced}
任意のチャネル $\mathcal E_B$ が局所 $B$ に作用しても $\rho_A$ が変わらないことを示せ。
:::solution
任意の $X$ を用い、$\tr[(X\otimes I)\sum_j(I\otimes K_j)\rho_{AB}(I\otimes K_j^\dagger)]$ を巡回させると $\tr[(X\otimes\sum K_j^\dagger K_j)\rho_{AB}]$ になる。Kraus の規格化で元の期待値に等しく、すべての $X$ で成り立つため $\rho_A$ が同じである。ただし個々の測定結果に条件づければ状態は一般に変わる。
:::
::::

::::exercise{level=advanced}
純粋状態 $\Psi$ の部分系を考え、$\rho_A$ と $\rho_B$ の von Neumann エントロピーが等しいことを示せ。
:::solution
Schmidt 分解から非零固有値 $\lambda_\alpha$ が一致する。余分な零固有値は $0\log0=0$ と定めるので寄与しない。従って両者とも $-k_B\sum_\alpha\lambda_\alpha\log\lambda_\alpha$ である。
:::
::::

::::exercise{level=advanced}
一定の $\Gamma_\uparrow=2\,\mathrm{s^{-1}},\Gamma_\downarrow=6\,\mathrm{s^{-1}}$ で、基底から始めた励起占有を求めよ。
:::solution
占有の式は $\dot p=2-8p$ で、定常値 $1/4$ との差が $e^{-8t}$ で減る。$t$ を秒で表せば $p(t)=(1-e^{-8t})/4$ である。
:::
::::

::::exercise{level=advanced}
真空浴で初期が三量子の振動子について $e^{-\kappa t}=1/3$ のとき二量子が残る確率を求めよ。
:::solution
母関数 $(2/3+z/3)^3$ の $z^2$ の係数を取る。$\binom32(1/3)^2(2/3)=2/9$ である。
:::
::::

::::exercise{level=advanced}
$\mathcal E(\rho)=\tfrac12\rho+\tfrac14 X\rho X+\tfrac14 Z\rho Z$ の Pauli 基底での固有値と Choi 行列の階数を求めよ。$\ket{+x},\ket{-x}$ を同頻度で入力し、一回の出力からどちらを入力したかを推定するときの最良の成功確率を求めよ。

:::solution
まず $\mathcal E(I)=I$ である。$X$ に対し $X$ の共役は $X$、$Z$ の共役は $-X$ なので $\lambda_x=1/2+1/4-1/4=1/2$ となる。$Z$ は逆の順に二つの係数が寄与して $\lambda_z=1/2$ である。$Y$ はどちらの共役でも符号が替わるので $\lambda_y=1/2-1/4-1/4=0$ である。従って Bloch 球では $y$ の方向を失い、$x,z$ は半分になる。

非規格化 Choi 行列では四つの Pauli 行列をベクトルにしたもののノルム二乗は 2 であり、確率係数を掛けた固有値は $1,1/2,0,1/2$ である。階数は 3 であり、最小の Kraus 表現は三本を必要とする。一方、入力の差の演算子は $X$ で、出力の差は $X/2$ である。トレース距離は $1/2$、等事前確率での最良の成功率は $(1+1/2)/2=3/4$ となる。実際 $x$ 軸で測ると、各入力で元と同じ符号になる確率は $3/4$ である。$y$ 軸の二つの入力なら同じ密度へ写るため成功率は $1/2$ になる。Kraus の本数だけからどの情報が残るかは決まらない。
:::
::::

## 確認テスト

::::quiz
全体が純粋な二体系で部分系 $A$ が混合なら、全体についてどれが言えるか。

- [x] もつれている
- [ ] 直積である
- [ ] 全体の密度行列の階数が 2 以上である

:::explanation
部分系が混合であることは縮約の階数が 2 以上であることを意味する。Schmidt 分解に二項以上あるので全体は直積ではなく、全体の純粋性は保たれる。
:::
::::

::::quiz{answer="0.5" tol="1e-9"}
最大混合な二準位状態の純度を求めよ。
:::explanation
二つの固有値は $1/2,1/2$ である。純度は二乗和なので $1/4+1/4=1/2$ となる。
:::
::::

::::quiz
一般の初期密度行列について純粋位相緩和で保存されるものを選べ。

- [x] 位相緩和を定義した基底での対角占有
- [ ] 非対角成分の大きさ
- [ ] 純度

:::explanation
$\mathcal D[\sqrt{\gamma_\varphi/2}\sigma_z]\rho$ の対角は零である。非対角は $\gamma_\varphi>0$ で指数減衰し、最初に非対角があれば純度も低下する。
:::
::::

::::quiz{answer="0.25" tol="1e-9"}
二準位の零温度振幅減衰で非対角成分の大きさが半分になった時刻の励起占有の割合を求めよ。
:::explanation
非対角成分は $e^{-\Gamma t/2}$、占有はその二乗の $e^{-\Gamma t}$ 倍である。よって占有は $(1/2)^2=1/4$ 倍になる。
:::
::::

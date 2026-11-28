---
title: 量子回路
description: 量子回路のゲート分解と普遍性を調べ、もつれを使う通信と測定型計算を構成する。
notation:
  - symbol: 'R_{\boldsymbol n}(\theta)'
    meaning: 軸 n の周りの量子ビット回転
  - symbol: 'H,S,T'
    meaning: Hadamard ゲートと二種類の位相ゲート
  - symbol: '\mathrm{CX}'
    meaning: 制御 NOT
---

## この章の目標

- 1 量子ビットゲートと Bloch 球上の回転の対応を導く。
- 制御ゲートによる計算と普遍ゲート集合による近似を区別する。
- テレポーテーション、超密度符号化、測定型計算を基底状態で追跡する。

## 1 量子ビットゲートと Bloch 球上の回転

回路では左から右へ時刻が進むとし、行列積では右端の演算が先に働く。計算基底の相対位相を変える操作も、別の基底で後から測れば読み出しに影響する。

:::definition[量子ゲート]{#def-qi-unitary-gate}
閉じた有限次元系へのユニタリ変換を :term[量子ゲート]{en="quantum gate" yomi="りょうしげーと"}という。単位軸 $\boldsymbol n$ に対し $R_{\boldsymbol n}(\theta)=e^{-i\theta\boldsymbol n\cdot\boldsymbol\sigma/2}$ とする。
:::

:::proposition[回転の随伴作用]{#prop-qi-bloch-rotation}
$R_{\boldsymbol n}(\theta)\rho R_{\boldsymbol n}(\theta)^\dagger$ の Bloch ベクトルは

$$
\boldsymbol r\cos\theta+(\boldsymbol n\times\boldsymbol r)\sin\theta+\boldsymbol n(\boldsymbol n\cdot\boldsymbol r)(1-\cos\theta)

$$

である。
:::

:::proof
$(\boldsymbol n\cdot\boldsymbol\sigma)^2=I$ なので指数を偶数項と奇数項に分けると、$R=cI-is\boldsymbol n\cdot\boldsymbol\sigma$、$c=\cos(\theta/2),s=\sin(\theta/2)$ となる。共役積の交差項は :eqref[eq-qi-pauli-product] から $2cs(\boldsymbol n\times\boldsymbol r)\cdot\boldsymbol\sigma$ である。$s^2$ の項は、反交換関係により $s^2(2(\boldsymbol n\cdot\boldsymbol r)\boldsymbol n-\boldsymbol r)\cdot\boldsymbol\sigma$ となる。$c^2$ の項と合算する。
:::

任意の 1 量子ビットユニタリの行列式を大域位相で 1 に直すと $\begin{pmatrix}a&-b^*\\b&a^*\end{pmatrix}$ と書ける。第 1 列が $(a,b)^\top$ なら、直交性と行列式 1 が第 2 列をこの形に決める。この行列は $R_z(\alpha)R_y(\beta)R_z(\gamma)$ で表せる。実際、第 1 列は $(e^{-i(\alpha+\gamma)/2}\cos(\beta/2),e^{i(\alpha-\gamma)/2}\sin(\beta/2))^\top$ なので絶対値と位相を別々に合わせられる。

:::example[基底と位相の変換]
$H=(X+Z)/\sqrt2$ は $\ket0,\ket1$ を $\ket+,\ket-$ に移す。$S=\diag(1,i)$、$T=\diag(1,e^{i\pi/4})$ とすると $T^2=S$、$S^2=Z$ である。$S\ket+=(\ket0+i\ket1)/\sqrt2$ は $X$ 軸上の状態を $Y$ 軸へ移したものになる。
:::

### Euler 角から行列を復元する
$R_j(\theta)=e^{-i\theta\sigma_j/2}$ とする。$U\in U(2)$ から $\det U=e^{2i\delta}$ となる $\delta$ を用いて位相を除くと、列の直交性と行列式 1 から $e^{-i\delta}U=[a,-\bar b;b,\bar a]$、$|a|^2+|b|^2=1$ と書ける。一方、行列を掛けると
$$
R_z(\alpha)R_y(\beta)R_z(\gamma)=
\begin{pmatrix}
e^{-i(\alpha+\gamma)/2}\cos(\beta/2)&-e^{-i(\alpha-\gamma)/2}\sin(\beta/2)\\
e^{i(\alpha-\gamma)/2}\sin(\beta/2)&e^{i(\alpha+\gamma)/2}\cos(\beta/2)
\end{pmatrix}.
$$
$\beta\in[0,\pi]$ を $\cos(\beta/2)=|a|$ で決め、非零の成分について $\alpha+\gamma=-2\arg a$、$\alpha-\gamma=2\arg b$ と選べば積が一致する。成分が零ならその位相条件は不要である。例えば $\beta=0$ では $\alpha+\gamma$ だけが行列に残る。このため Euler 角を個別に復元する問題は、極で一意性を失う。

:::example[Hadamard の回転を分解する]
$H=[1,1;1,-1]/\sqrt2$ に対して $-iH$ は行列式 1 で、$a=b=-i/\sqrt2$ である。$\beta=\pi/2$、$\alpha+\gamma=\pi$、$\alpha-\gamma=-\pi$ なので $\alpha=0$、$\gamma=\pi$ を取れる。したがって $H=iR_y(\pi/2)R_z(\pi)$ である。$\ket0$ への作用では $R_z(\pi)$ が $-i$ を掛け、$R_y(\pi/2)$ が $\ket+$ を作るので、大域位相を戻せば $H\ket0=\ket+$ となる。$\ket1$ に対しては $i(-\ket0+\ket1)/\sqrt2$ に $i$ を掛けて $\ket-$ を得る。
:::

:::example[同じ角度で順序だけを入れ替える]
$\ket0$ の初期 Bloch 軸は $\boldsymbol e_z$ である。$R_y(\pi/2)$ を先に掛け、次に $R_z(\pi/2)$ を掛ければ、軸は $\boldsymbol e_z\to\boldsymbol e_x\to\boldsymbol e_y$ と移る。逆順では第 1 ゲートが $\boldsymbol e_z$ を固定し、第 2 ゲートで $\boldsymbol e_x$ に移る。両方とも $Z$ の確率は $1/2$ だが、$Y$ の $+1$ 確率は順に $1,1/2$ である。
:::

### ゲート積を Bloch 回転と振幅の両方で検査する

入力を $\ket0$ とし、各ゲートを掛けた後で ket の座標を出す。出力の Bloch 成分は Pauli 期待値であり、座標からも $r_x=2\Re(\bar a b)$、$r_y=2\Im(\bar a b)$、$r_z=|a|^2-|b|^2$ で計算できる。ゲートの順番が違えば、同じ二つの回転を使っても最終軸が変わる場合がある。

:::example[Y 回転の後に Z 回転]

第 1 ゲートを $U$、第 2 を $V$ とする。積 $VU$ を適用し、各行列を作用の順に追う。

$$
\begin{aligned}
U&=\left(\begin{matrix}\frac{\sqrt{3}}{2} & - \frac{1}{2}\\\frac{1}{2} & \frac{\sqrt{3}}{2}\end{matrix}\right)\\
V&=\left(\begin{matrix}\frac{\sqrt{2} \left(1 - i\right)}{2} & 0\\0 & \frac{\sqrt{2} \left(1 + i\right)}{2}\end{matrix}\right)\\
U\ket0&=\left(\begin{matrix}\frac{\sqrt{3}}{2}\\\frac{1}{2}\end{matrix}\right)\\
VU\ket0&=\left(\begin{matrix}\frac{\sqrt{6} \left(1 - i\right)}{4}\\\frac{\sqrt{2} \left(1 + i\right)}{4}\end{matrix}\right)\\
\boldsymbol r_{\mathrm{out}}&=(0,\frac{\sqrt{3}}{2},\frac{1}{2})\\
(p_X(+),p_Y(+),p_Z(+))&=(\frac{1}{2},\frac{\sqrt{3}}{4} + \frac{1}{2},\frac{3}{4})
\end{aligned}
$$

第一段で正の $X$ 成分が生まれ、第二段でその成分が正の $Y$ へ移る。計算基底での確率差は第二段の $Z$ 回転では変わらない。
:::

:::example[前後を逆にする]

第 1 ゲートを $U$、第 2 を $V$ とする。積 $VU$ を適用し、各行列を作用の順に追う。

$$
\begin{aligned}
U&=\left(\begin{matrix}\frac{\sqrt{2} \left(1 - i\right)}{2} & 0\\0 & \frac{\sqrt{2} \left(1 + i\right)}{2}\end{matrix}\right)\\
V&=\left(\begin{matrix}\frac{\sqrt{3}}{2} & - \frac{1}{2}\\\frac{1}{2} & \frac{\sqrt{3}}{2}\end{matrix}\right)\\
U\ket0&=\left(\begin{matrix}\frac{\sqrt{2} \left(1 - i\right)}{2}\\0\end{matrix}\right)\\
VU\ket0&=\left(\begin{matrix}\frac{\sqrt{6} \left(1 - i\right)}{4}\\\frac{\sqrt{2} \left(1 - i\right)}{4}\end{matrix}\right)\\
\boldsymbol r_{\mathrm{out}}&=(\frac{\sqrt{3}}{2},0,\frac{1}{2})\\
(p_X(+),p_Y(+),p_Z(+))&=(\frac{\sqrt{3}}{4} + \frac{1}{2},\frac{1}{2},\frac{3}{4})
\end{aligned}
$$

第一段の入力は $Z$ の正固有状態なので、先に掛けた位相は全体に掛かる。第二段では相対位相が実の回転となり、最終の $Y$ 成分は零のままである。
:::

:::example[X 回転から Y 回転]

第 1 ゲートを $U$、第 2 を $V$ とする。積 $VU$ を適用し、各行列を作用の順に追う。

$$
\begin{aligned}
U&=\left(\begin{matrix}\frac{\sqrt{2}}{2} & - \frac{\sqrt{2} i}{2}\\- \frac{\sqrt{2} i}{2} & \frac{\sqrt{2}}{2}\end{matrix}\right)\\
V&=\left(\begin{matrix}\frac{\sqrt{3}}{2} & - \frac{1}{2}\\\frac{1}{2} & \frac{\sqrt{3}}{2}\end{matrix}\right)\\
U\ket0&=\left(\begin{matrix}\frac{\sqrt{2}}{2}\\- \frac{\sqrt{2} i}{2}\end{matrix}\right)\\
VU\ket0&=\left(\begin{matrix}\frac{\sqrt{6}}{4} + \frac{\sqrt{2} i}{4}\\\frac{\sqrt{2}}{4} - \frac{\sqrt{6} i}{4}\end{matrix}\right)\\
\boldsymbol r_{\mathrm{out}}&=(0,-1,0)\\
(p_X(+),p_Y(+),p_Z(+))&=(\frac{1}{2},0,\frac{1}{2})
\end{aligned}
$$

第一段が状態を負の $Y$ 軸へ移すので、その方向は第二段の回転軸と平行である。従って第二段は位相を変えうるが、最終状態の Bloch ベクトルを変えない。
:::

### SU(2) の二重被覆を状態の軸から見る

Pauli 展開から $U\in\mathrm{SU}(2)$ は実数 $a,\boldsymbol b$ を用いて $U=aI-i\boldsymbol b\cdot\boldsymbol\sigma$ と書け、$a^2+\|\boldsymbol b\|^2=1$ である。実際、列の直交性と行列式 1 を使って $[u,-\bar v;v,\bar u]$ と書き、$u,v$ を実部と虚部に分ければよい。$\boldsymbol b$ が零でなければ、$a=\cos(\theta/2),\boldsymbol b=\boldsymbol n\sin(\theta/2)$ と置ける。$\boldsymbol b=0$ なら $U=\pm I$ であり、回転軸は任意である。共役作用は行列の積を保つので、これに対応する Bloch 回転の写像は群準同型である。全ての三次元回転には軸と角度があるので全射である。核では $X,Y,Z$ の全てと可換である。$Z$ との可換性で対角、$X$ との可換性で両対角成分が等しくなり、行列式 1 から $\pm I$ を得る。従って同じ Bloch 回転を表す $\mathrm{SU}(2)$ の元は二つあり、符号だけが異なる。

:::example[一周と二周を回路の比較で見分ける]
$R_z(2\pi)=-I$、$R_z(4\pi)=I$ であり、どちらも単独で Bloch ベクトルを固定する。しかし補助ビットを $\ket+$ として、$1$ の枝だけに $R_z(2\pi)$ を掛けると、$\ket+$ は $\ket-$ に移る。対して $4\pi$ の制御では不変である。軸回転として見える周期と、制御回路で振幅の符号を比較した周期が異なる。制御を作る場合に大域位相を追う必要性は、この二重被覆にも表れる。
:::

### 回転を二つ重ねた軸を求める

$U_j=a_jI-i\boldsymbol b_j\cdot\boldsymbol\sigma$ と書くと、Pauli 積により
$$
U_2 U_1=(a_2a_1-\boldsymbol b_2\cdot\boldsymbol b_1)I
-i(a_2\boldsymbol b_1+a_1\boldsymbol b_2+\boldsymbol b_2\times\boldsymbol b_1)\cdot\boldsymbol\sigma .
$$
従って積の回転軸には二つの元の軸の平面に直交する成分が現れる。ゲートの順序を逆にすると外積の符号が変わる。これは可換でないことの具体的な計算である。

:::example[二つの四分の一回転が三分の一周になる]
$U_1=R_x(\pi/2)$、$U_2=R_z(\pi/2)$ では $a_1=a_2=1/\sqrt2$、$\boldsymbol b_1=\boldsymbol e_x/\sqrt2$、$\boldsymbol b_2=\boldsymbol e_z/\sqrt2$ である。積の $a$ は $1/2$、$\boldsymbol b=(1,1,1)/2$ となる。従って軸は $(1,1,1)/\sqrt3$、角は $2\pi/3$ である。ベクトルの回転で照合すると $\boldsymbol e_x\to\boldsymbol e_y$、$\boldsymbol e_y\to\boldsymbol e_z$、$\boldsymbol e_z\to\boldsymbol e_x$ で、向きも一致する。順序を逆にすると $\boldsymbol b=(1,-1,1)/2$ になるため同じ角度でも異なる操作である。
:::

## 制御ゲートと多量子ビットゲート

:term[制御ゲート]{en="controlled gate" yomi="せいぎょげーと"}を $C(U)=\ketbra00\otimes I+\ketbra11\otimes U$ とする。$C(X)$ を $\mathrm{CX}$ と書く。制御レジスタの値が重ね合わせであっても、基底上の作用を線形に広げればよい。

:::example[Bell 状態の準備]
$\mathrm{CX}(H\otimes I)\ket{00}=(\ket{00}+\ket{11})/\sqrt2$ である。一般の $\ket{ab}$ は $(\ket{0,b}+(-1)^a\ket{1,1\oplus b})/\sqrt2$ に移る。したがって逆回路 $(H\otimes I)\mathrm{CX}$ は四つの Bell 状態から $(a,b)$ を読み出す。ここで各ゲートは自己逆なので、順番を逆にすればよい。
:::

:::example[交換と可逆論理]
$\mathrm{CX}_{1\to2}\mathrm{CX}_{2\to1}\mathrm{CX}_{1\to2}$ は $(a,b)\to(a,a\oplus b)\to(b,a\oplus b)\to(b,a)$ と進むので SWAP である。Toffoli ゲートは $\ket{a,b,c}\mapsto\ket{a,b,c\oplus ab}$ と作用する。空白の第 3 ビットを 0 とすれば AND を計算し、入力を残すので逆演算も実行できる。
:::

制御版を作る際には位相に注意する。例えば $I$ と $-I$ は単体では同じ状態変換だが $C(-I)=Z\otimes I$ なので、制御ビットの相対位相を変える。

### ゲートの大域位相が制御で相対位相に変わる

$C(U)=\ketbra00\otimes I+\ketbra11\otimes U$ とする。$U$ を $e^{i\delta}U$ へ置き換えると、制御ビットに $P_\delta=\diag(1,e^{i\delta})$ を掛けた回路になる。標的を作用させる分枝の位相だけが変わるので、制御された二つのゲートは一般に異なる。$C(U)$ が局所ゲートの積になるのは $U$ がスカラーの場合に限る。逆を示すには、$U$ がスカラーでなければ異なる固有値を持つ二本の固有ベクトルを足し、標的 $\ket v$ とする。$U\ket v$ と $\ket v$ は平行でない。制御 $\ket+$ との積を入れた後の係数行列ではこの二ベクトルが独立な行をなすので、出力は積でない。

:::example[制御された負号]
$U=-I$ なら単独で状態を変えない。ところが $C(-I)=Z\otimes I$ であるから、入力 $\ket+\ket0$ を $\ket-\ket0$ に移す。制御を $X$ で測れば $-1$ が確実に出る。$C(I)$ では結果は $+1$ である。
:::

:::example[逆向きの制御 NOT]
$\mathrm{CX}_{1\to2}=(I\otimes H)\mathrm{CZ}(I\otimes H)$ であり、同様に $\mathrm{CX}_{2\to1}=(H\otimes I)\mathrm{CZ}(H\otimes I)$ である。$\mathrm{CZ}$ は $\ket{11}$ にだけ $-1$ を掛けるので二量子ビットを交換しても変わらない。二つの等式を合わせると $(H\otimes H)\mathrm{CX}_{1\to2}(H\otimes H)=\mathrm{CX}_{2\to1}$ となる。例えば $\ket{01}$ は右辺で $\ket{11}$ に移り、$\ket{10}$ は固定される。どちらが制御であるかは回路の見た目だけでなく各基底の像から確認できる。
:::

### 制御ビットで条件付けた行列を計算する

制御ビットが重ね合わせのとき、$U$ を作用させる枝には、標的を変えることと制御との相対位相を変えることの二効果がある。入力 $\ket+\ket0$ に制御 $U$ を掛けた後、制御だけの縮約を作ると、非対角成分は $\braket{U0}{0}/2$ である。標的側を無視するときにこの内積が残る。

:::example[標的が直交方向へ移る制御 NOT]

制御が第 1、標的が第 2 の順で、入力は $\ket+\ket0$ である。制御の $X,Y$ の平均値も追う。

$$
\begin{aligned}
U&=\left(\begin{matrix}0 & 1\\1 & 0\end{matrix}\right)\\
\operatorname{C}(U)&=\left(\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 0 & 1\\0 & 0 & 1 & 0\end{matrix}\right)\\
\ket{\psi_{\mathrm{out}}}&=\left(\begin{matrix}\frac{\sqrt{2}}{2}\\0\\0\\\frac{\sqrt{2}}{2}\end{matrix}\right)\\
\rho_{\mathrm{control}}&=\Tr_{\mathrm{target}}\ketbra{\psi_{\mathrm{out}}}{\psi_{\mathrm{out}}}=\left(\begin{matrix}\frac{1}{2} & 0\\0 & \frac{1}{2}\end{matrix}\right)\\
(\langle X\rangle,\langle Y\rangle)&=(0,0)
\end{aligned}
$$

標的の二枝が直交するので、制御の縮約は完全混合である。にもかかわらず共同状態は純粋であり、個別の枝を捨てたことが局所の混合を生む。
:::

:::example[制御回転に角度を残す]

制御が第 1、標的が第 2 の順で、入力は $\ket+\ket0$ である。制御の $X,Y$ の平均値も追う。

$$
\begin{aligned}
U&=\left(\begin{matrix}\frac{\sqrt{3}}{2} & - \frac{1}{2}\\\frac{1}{2} & \frac{\sqrt{3}}{2}\end{matrix}\right)\\
\operatorname{C}(U)&=\left(\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & \frac{\sqrt{3}}{2} & - \frac{1}{2}\\0 & 0 & \frac{1}{2} & \frac{\sqrt{3}}{2}\end{matrix}\right)\\
\ket{\psi_{\mathrm{out}}}&=\left(\begin{matrix}\frac{\sqrt{2}}{2}\\0\\\frac{\sqrt{6}}{4}\\\frac{\sqrt{2}}{4}\end{matrix}\right)\\
\rho_{\mathrm{control}}&=\Tr_{\mathrm{target}}\ketbra{\psi_{\mathrm{out}}}{\psi_{\mathrm{out}}}=\left(\begin{matrix}\frac{1}{2} & \frac{\sqrt{3}}{4}\\\frac{\sqrt{3}}{4} & \frac{1}{2}\end{matrix}\right)\\
(\langle X\rangle,\langle Y\rangle)&=(\frac{\sqrt{3}}{2},0)
\end{aligned}
$$

二つの標的状態の重なりは $\sqrt3/2$ であり、制御の $X$ 平均はその値になる。角度が大きくなると、標的にどちらの枝かを表す情報が増える。
:::

:::example[標的は不変で枝の位相だけが変化]

制御が第 1、標的が第 2 の順で、入力は $\ket+\ket0$ である。制御の $X,Y$ の平均値も追う。

$$
\begin{aligned}
U&=\left(\begin{matrix}i & 0\\0 & i\end{matrix}\right)\\
\operatorname{C}(U)&=\left(\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & i & 0\\0 & 0 & 0 & i\end{matrix}\right)\\
\ket{\psi_{\mathrm{out}}}&=\left(\begin{matrix}\frac{\sqrt{2}}{2}\\0\\\frac{\sqrt{2} i}{2}\\0\end{matrix}\right)\\
\rho_{\mathrm{control}}&=\Tr_{\mathrm{target}}\ketbra{\psi_{\mathrm{out}}}{\psi_{\mathrm{out}}}=\left(\begin{matrix}\frac{1}{2} & - \frac{i}{2}\\\frac{i}{2} & \frac{1}{2}\end{matrix}\right)\\
(\langle X\rangle,\langle Y\rangle)&=(0,1)
\end{aligned}
$$

標的はどの枝でも同じ。このため共同状態は積であり、制御の $Y$ の期待値が 1 になる。$U$ を単独で使う場合には全体に掛かる位相が、制御に現れる。
:::

### 制御した演算を位相から構成する

$U=e^{i\delta}A X B X C$、$ABC=I$ と分解できるとき、標的へ $C$、制御 NOT、$B$、制御 NOT、$A$ の順に作用させる。制御 0 の枝は $ABC=I$、1 の枝は $AXBXC=e^{-i\delta}U$ になる。最後に制御へ $\diag(1,e^{i\delta})$ を掛けると目的の制御ゲートになる。行列積は回路の時系列と逆順に書くので、先頭の因子を最初に実装しない点に注意する。前節の Euler 分解で $A,B,C$ を求めれば、どの $U$ でもこの構成を適用できる。

:::example[制御位相を二つの CX で実装する]
$U=\diag(1,e^{i\theta})$ は $e^{i\theta/2}R_z(\theta)$ である。$A=R_z(\theta/2)$、$B=R_z(-\theta/2)$、$C=I$ とすれば $AXBX=R_z(\theta)$ になる。従って制御側に位相 $\theta/2$ を掛ける。制御 0 の枝では二回の $R_z$ が消える。制御 1 の枝では $\ket0$ の $-\theta/2$ の位相が制御側の $\theta/2$ と消え、$\ket1$ には合計 $\theta$ が残る。例えば $\theta=\pi/3$ とすれば、計算基底の最後の成分だけに $e^{i\pi/3}$ を掛ける。標的だけでの回転と比較して係数 2 が現れないよう、制御側の位相まで含めて確かめる。
:::

### SWAP を分解したときの中間ラベル

二列 $(a,b)$ を第 1 制御の $\mathrm{CX}$、第 2 制御の $\mathrm{CX}$、再び第 1 制御の $\mathrm{CX}$ に通す。中間は $(a,a\oplus b)$、$(b,a\oplus b)$ で、最後は $(b,a)$ になる。すべての計算基底で一致するため、線形性で任意の重ね合わせに作用する。近接位置だけにゲートを掛けられる装置では、遠い位置へデータを移す手段に使える。ただし中間に測定を入れると、途中のラベルが作業記録に残るため、入れ替えられた量子状態を同じようには得られない。

:::example[第三系との相関を残しながら位置を移す]
データ $A$ が参照 $R$ と $\Phi^+$ を共有し、補助 $B$ は $\ket0$ とする。$A,B$ の SWAP の前は $(\ket{0,0,0}+\ket{1,1,0})/\sqrt2$ を $R,A,B$ の順で持つ。第 1 ゲートの後は $(\ket{000}+\ket{111})/\sqrt2$、第 2 の後は $(\ket{000}+\ket{101})/\sqrt2$ になる。最後のゲートでは両成分の $A=0$ なので変わらない。従って $R,B$ が Bell 対を共有し、$A$ は零である。もし最初のゲート後に $B$ を測って結果を捨てた場合、交差項が消え、$R,B$ の最終状態は $\ket{00},\ket{11}$ の古典混合になる。
:::

## 普遍ゲート集合と Solovay–Kitaev の定理

選んだ二つの基底状態だけに作用する二準位ユニタリを使うと、大きな行列を消去できる。列の二成分 $a,b$ がともに零ならその組を飛ばし、それ以外では行列 $\begin{pmatrix}\bar a&\bar b\\-b&a\end{pmatrix}/\sqrt{|a|^2+|b|^2}$ を作用させて第 2 成分を零にする。$d\times d$ 行列の第 $j$ 列では第 $j$ 行と各 $k>j$ 行を組に選ぶ。これらの行の第 $j-1$ 列までの成分は既にすべて零なので、操作後も零である。第 $d-1$ 列まで処理すると上三角のユニタリを得る。第 1 列には先頭の絶対値 1 の成分しかなく、直交性により他の列の先頭は零になる。残ったブロックに同じ議論を繰り返すと、行列は対角である。各基底に残る位相は、その基底のラベルだけに条件付けて実装する。

任意の二つのビット列は、異なるビットを一つずつ反転する経路で結べる。経路の最後の辺まで一方の基底状態を隣接交換で移し、最後の辺に条件付き 1 量子ビットゲートを適用して交換を戻せば、元の二準位ユニタリになる。他のビットの AND を Toffoli で一時レジスタに計算し、制御に使ってから逆計算すれば多重制御になる。0 による制御は前後の $X$ で 1 制御に変換できる。

Toffoli 自体も、$V^2=X$ となる平方根 $V$ を使い、標的へ順に $b$ 制御 $V$、$a$ 制御で $b$ を反転、$b$ 制御 $V^\dagger$、$a$ 制御で $b$ を復元、$a$ 制御 $V$ とすれば実現できる。$V$ の指数は $b-(a\oplus b)+a=2ab$ である。任意の $C(U)$ は $U=e^{i\delta}R_z(\alpha)R_y(\beta)R_z(\gamma)$ と書き、$A=R_z(\alpha)R_y(\beta/2)$、$B=R_y(-\beta/2)R_z(-(\alpha+\gamma)/2)$、$C=R_z((\gamma-\alpha)/2)$ に分けると $ABC=I$、$AXBXC=e^{-i\delta}U$ となる。従って二つの $\mathrm{CX}$ と 1 量子ビットゲート、制御側の位相で実装できる。

### 有限集合 $\{H,T,\mathrm{CX}\}$ の普遍性

任意の 1 量子ビット演算を正確に使える仮定を、有限ゲート集合での近似へ置き換える。$T$ は位相を除き $R_z(\pi/4)$ であり、$HTH$ は $R_x(\pi/4)$ である。従って $U=R_x(\pi/4)R_z(\pi/4)$ が位相を除き生成される。$\mathrm{SU}(2)$ の固有位相を $\pm\omega$ とすると、Pauli の積の跡を取って $2\cos\omega=2\cos^2(\pi/8)=1+\sqrt2/2$ である。

$\omega/\pi$ が有理なら $e^{i\omega}$ は 1 の冪根であり、その逆との和 $2\cos\omega$ は代数的整数となる。一方、$1+\sqrt2/2$ とその共役の積は $1/2$ である。代数的整数の共役と積も代数的整数であり、有理数の代数的整数は整数に限るので矛盾する。従って $U$ の冪はその回転軸を持つ 1 パラメータ部分群に稠密である。円周上の無理数回転の稠密性は、鳩の巣原理で任意に短い非零の差を整数倍の間に見付け、その差の反復で円周を近似すれば従う。

$U$ の軸は $x,z$ の成分をともに持つ。これは $R_xR_z$ を展開したとき、両軸の一次項の係数がともに $\cos(\pi/8)\sin(\pi/8)$ になるためである。この軸を $R_z(\pi/4)$ で回せば平行でない軸を得る。閉包にはこの二軸の全角度の回転が入る。一方の軸の周りに他方の軸を $\pi/2$ 回した第 3 軸も共役で得られ、三つの軸は一次独立である。三軸の回転を順に掛ける角度写像を原点で微分すると軸が列に並び、逆関数定理で恒等演算の近傍が生成される。連結な $\mathrm{SU}(2)$ の中で近傍を含む部分群は全群である。実際、各剰余類が開になるので、真部分群なら群を互いに素な二つの開集合へ分離してしまう。こうして $H,T$ が任意の 1 量子ビットゲートを大域位相を除いて近似できることが分かる。

:::theorem[Solovay–Kitaev の定理]{#thm-qi-solovay-kitaev}
固定次元 $d\ge2$ で、逆元を含む $\mathrm{SU}(d)$ の有限ゲート集合が稠密部分群を生成するとする。任意の標的 $U\in\mathrm{SU}(d)$ を作用素ノルム誤差 $\eps$ 以下で近似する長さ $O(\log^c(1/\eps))$ のゲート列を求められる。$c=\log_{3/2}5$ と取れる。$U(d)$ のゲートを扱うときは位相を除いて $\mathrm{SU}(d)$ に直してからこの定理を使う。
:::

### 固定次元での再帰の証明

以後、ノルムは作用素ノルムである。近恒等ゲートを、二つの小さいゲートの交換子で近似する。小さいゲート自体の合成が少しずれても、積全体のずれをさらに小さくできることが要点である。

:::lemma[近恒等ゲートの交換子近似]{#lem-qi-balanced-commutator}
固定次元 $d\ge2$ の $W\in\mathrm{SU}(d)$ が $\norm{W-I}\le\eps$ を満たすとする。十分小さい $\eps$ に対し、$V_1,V_2\in\mathrm{SU}(d)$ が存在し、$\norm{V_j-I}=O(\sqrt\eps)$ かつ

$$

\norm{W-V_1V_2V_1^\dagger V_2^\dagger}=O(\eps^{3/2})

$$

となる。
:::

:::proof
固有位相の主値で $W=e^K$ と書くと、$K^\dagger=-K$、$\norm K=O(\eps)$ である。$\det W=1$ より固有位相の和は $2\pi$ の整数倍であり、固定次元で $\eps$ が十分小さければ和は零なので $\Tr K=0$ である。$K$ を対角化し、対角行列を全成分の絶対値が $1/\sqrt d$ の Fourier 行列で共役すると、対角成分がすべて $\Tr K/d=0$ となる。以後この基底で計算する。

$K=0$ なら $V_1=V_2=I$ とすればよい。$K\ne0$ なら $u=\sqrt{\norm K}$ とし、$P$ を対角成分 $iu(j-(d+1)/2)$ の行列とする。$Q$ の対角を零、非対角成分を $Q_{jk}=K_{jk}/(iu(j-k))$ とする。$Q_{kj}=-\overline{Q_{jk}}$ なので両者は跡零の反自己共役行列である。成分ごとに $[P,Q]=K$ となり、$\norm P,\norm Q=O_d(u)$ である。

$V_1=e^P,V_2=e^Q$ と置く。積 $e^Pe^Qe^{-P}e^{-Q}$ の展開で定数項は $I$、一次項は消え、二次項は $PQ-QP$ である。他の項のノルム和は、各指数の絶対収束級数を掛け合わせて $O((\norm P+\norm Q)^3)$ と評価できる。一方 $W=I+K+O(\norm K^2)$ なので主張を得る。最後に共役を戻しても同じノルム評価が成り立つ。
:::

:::lemma[交換子の安定性]{#lem-qi-commutator-stability}
$V_1,V_2$ とその近似 $\tilde V_1,\tilde V_2$ はユニタリとする。各近似誤差は $\eta$ 以下、元の各ゲートの恒等からの距離は $\delta$ 以下とする。このとき交換子の誤差は高々 $4\eta\delta+2\eta^2$ である。
:::

:::proof
$C(A,B)=ABA^\dagger B^\dagger$ と書く。第 1 因子だけを替えるときは、$ABA^\dagger=I+A(B-I)A^\dagger$ なので $\norm{C(A,B)-C(A',B)}\le2\norm{A-A'}\norm{B-I}$ である。第 2 因子だけの差は $C(A,B)^\dagger=C(B,A)$ と同じ見積りから得る。初めに第 1 因子を替え、次に第 2 因子を替えると、二つの上界はそれぞれ $2\eta\delta$、$2\eta(\delta+\eta)$ である。
:::

稠密性とコンパクト性から、固定した小さい $\eps_0$ で全群を覆う有限個のゲート列を事前に用意できる。その最大長を $L_0$ とする。レベル $j$ で標的 $U$ の近似 $U_j$ が得られ、誤差が $\eps_j$ 以下とする。残差 $UU_j^\dagger$ は恒等の $\eps_j$ 近傍にある。:ref[lem-qi-balanced-commutator] で残差に対応する二ゲートを選び、それぞれをレベル $j$ の手続で近似する。逆元を含むゲート集合なので、近似回路の逆を同じ長さで作れる。:ref[lem-qi-commutator-stability] に $\eta=\eps_j,\delta=O(\sqrt{\eps_j})$ を入れると、近似交換子を $U_j$ の左に掛けた誤差は $\eps_{j+1}\le K\eps_j^{3/2}$ である。長さは $L_{j+1}\le5L_j$ である。

$\eps_0<K^{-2}$ に取れば、$K^2\eps_j$ は高々 $(K^2\eps_0)^{(3/2)^j}$ である。従って精度 $\eps$ までの段数は $\log_{3/2}\log(1/\eps)+O(1)$ であり、長さは $O(\log^{\log_{3/2}5}(1/\eps))$ となる。一段の再帰では粗い近似と二つの交換子因子のために再帰を三回呼ぶ。固定次元で厳密な行列演算を一操作として数えると、再帰木の節点数は $O(\log^{\log_{3/2}3}(1/\eps))$ である。回路記述を出力する時間には列の長さも含める。$H,T$ の場合は対応する $\mathrm{SU}(2)$ の代表元と逆元に定理を適用し、物理ゲートの列では $H^{-1}=H$、$T^{-1}=T^7$ を使える。その列に付く全体位相は単独のゲートとして用いる限り状態の統計を変えない。

:::example[回路全体の誤差予算]
100 個の 1 量子ビット演算をそれぞれノルム誤差 $10^{-5}$ 以下で近似する。元の積との差を一因子ずつ替える望遠和で書くと、各差の両側にある因子はユニタリなのでその項のノルムを変えない。従って全体は $100\cdot10^{-5}=10^{-3}$ 以下である。標的が制御ゲートで用いられる場合は、条件付きの位相もゲート分解に含めてから同じ誤差を適用する。
:::

### 厳密合成と近似合成

有限種類のゲートを有限個並べて得る演算の集合は可算である。一方、例えば $P_\theta=\diag(1,e^{i\theta})$ は全体位相を無視しても $\theta\in[0,2\pi)$ に非可算個ある。そのため有限集合の普遍性では、一般には誤差を指定して近似する。作用素ノルムで $\norm{U_j-V_j}\le\eps_j$ なら
$$
U_m\cdots U_1-V_m\cdots V_1
=\sum_{j=1}^m U_m\cdots U_{j+1}(U_j-V_j)V_{j-1}\cdots V_1 .
$$
全外側因子がユニタリであれば全体の誤差は $\sum_j\eps_j$ 以下である。大域位相を各 $V_j$ に合わせてからこの評価を使えるが、ゲートが制御下にある場合は、その位相も制御側に移して一致させる。

:::example[位相の誤差を直接測る]
$P_\theta$ を $P_\phi$ で近似すると $\norm{P_\theta-P_\phi}=|e^{i\theta}-e^{i\phi}|=2|\sin((\theta-\phi)/2)|$ である。$\theta-\phi=\pi/6$ なら $\sqrt{2-\sqrt3}$ である。一方、両方に $\ket+$ を入力したときの二出力のトレース距離は $|\sin((\theta-\phi)/2)|$ で、ちょうどその半分になる。同じ位相誤差を二連続で掛ければ位相差は $\pi/3$ になり、距離は $1/2$ となる。ノルムの和による上界は、誤差の符号が相殺する場合も含む上界である。
:::

:::example[二準位回転を四次元で選択する]
$\ket{00},\ket{11}$ にだけ $R_y(\beta)$ を作用させる。$\mathrm{CX}_{1\to2}$ により両基底は $\ket{00},\ket{10}$ に移る。この間は第 2 ビットが 0 で共通なので、第 2 ビットの 0 を条件に第 1 ビットを回転させる。最後に $\mathrm{CX}_{1\to2}$ をもう一度掛ける。$\beta=\pi/2$ で $\ket{00}$ は $\Phi^+$ になり、$\ket{01},\ket{10}$ は途中で第 2 ビットが 1 となり回転を受けない。基底状態の順序を逆にした場合には回転の非対角成分の符号も確認する必要がある。
:::

### 位相を除かずに有限ゲート列を組み立てる

位相ゲートの積は計算基底上で指数を加える。$T=\diag(1,e^{i\pi/4})$、$S=T^2$ とする。回路の置換でどの基底に位相を選択的に付けるかを変えられる。ここでは各ゲート列の積を直接作り、目的の作用と位相まで一致するか調べる。

:::example[位相ゲートを X 軸へ共役する]

回路は $U_{\mathrm{seq}}=HTH$ であり、右端から作用する。二次元の計算基底で、全体の位相も保持する。

$$
\begin{aligned}
U_{\mathrm{seq}}&=\left(\begin{matrix}\frac{\sqrt{2}}{4} + \frac{1}{2} + \frac{\sqrt{2} i}{4} & - \frac{\sqrt{2}}{4} + \frac{1}{2} - \frac{\sqrt{2} i}{4}\\- \frac{\sqrt{2}}{4} + \frac{1}{2} - \frac{\sqrt{2} i}{4} & \frac{\sqrt{2}}{4} + \frac{1}{2} + \frac{\sqrt{2} i}{4}\end{matrix}\right)\\
U^\dagger_{\mathrm{seq}}U_{\mathrm{seq}}&=\left(\begin{matrix}1 & 0\\0 & 1\end{matrix}\right)\\
\det U_{\mathrm{seq}}&=\frac{\sqrt{2} \left(1 + i\right)}{2}\\
U_{\mathrm{seq}}\ket0&=\left(\begin{matrix}\frac{\sqrt{2}}{4} + \frac{1}{2} + \frac{\sqrt{2} i}{4}\\- \frac{\sqrt{2}}{4} + \frac{1}{2} - \frac{\sqrt{2} i}{4}\end{matrix}\right)\\
U_{\mathrm{seq}}\ket1&=\left(\begin{matrix}- \frac{\sqrt{2}}{4} + \frac{1}{2} - \frac{\sqrt{2} i}{4}\\\frac{\sqrt{2}}{4} + \frac{1}{2} + \frac{\sqrt{2} i}{4}\end{matrix}\right)
\end{aligned}
$$

$H$ は $Z$ を $X$ に移すので、位相を除いた部分は $R_x(\pi/4)$ になる。全体の位相は行列式から読み取れる。
:::

:::example[T の繰り返しで S を作る]

回路は $U_{\mathrm{seq}}=T^2$ である。前の例と同じ二次元の基底で行列を計算する。

$$
\begin{aligned}
U_{\mathrm{seq}}&=\left(\begin{matrix}1 & 0\\0 & i\end{matrix}\right)\\
U^\dagger_{\mathrm{seq}}U_{\mathrm{seq}}&=\left(\begin{matrix}1 & 0\\0 & 1\end{matrix}\right)\\
\det U_{\mathrm{seq}}&=i\\
U_{\mathrm{seq}}\ket0&=\left(\begin{matrix}1\\0\end{matrix}\right)\\
U_{\mathrm{seq}}\ket1&=\left(\begin{matrix}0\\i\end{matrix}\right)
\end{aligned}
$$

$T$ の二回の作用は両基底で $\diag(1,i)$ になる。$Z$ を作るには四回、恒等には八回である。従って有限集合に逆元を入れる場合、$T^\dagger$ は $T^7$ で代用できる。
:::

:::example[CX で偶奇の位相を局所へ移す]

回路は $U_{\mathrm{seq}}=\mathrm{CX}(I\otimes S)\mathrm{CX}$ である。ここでは四次元の計算基底を辞書順に並べる。

$$
\begin{aligned}
U_{\mathrm{seq}}&=\left(\begin{matrix}1 & 0 & 0 & 0\\0 & i & 0 & 0\\0 & 0 & i & 0\\0 & 0 & 0 & 1\end{matrix}\right)\\
U^\dagger_{\mathrm{seq}}U_{\mathrm{seq}}&=\left(\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right)\\
\det U_{\mathrm{seq}}&=-1\\
U_{\mathrm{seq}}\ket0&=\left(\begin{matrix}1\\0\\0\\0\end{matrix}\right)\\
U_{\mathrm{seq}}\ket{11}&=\left(\begin{matrix}0\\0\\0\\1\end{matrix}\right)
\end{aligned}
$$

標的に書かれた偶奇だけが位相を受け、最後の $\mathrm{CX}$ がビット列を元に戻す。標的値だけを条件にせず、入力 01 と 10 の両方で位相が出る。
:::

:::example[偶奇の位相と局所位相で CZ になる]

回路は $U_{\mathrm{seq}}=(S\otimes S)\mathrm{CX}(I\otimes S^\dagger)\mathrm{CX}$ である。偶奇の位相を前の例と逆にし、その後で両ビットに $S$ を掛ける。

$$
\begin{aligned}
U_{\mathrm{seq}}&=\left(\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & -1\end{matrix}\right)\\
U^\dagger_{\mathrm{seq}}U_{\mathrm{seq}}&=\left(\begin{matrix}1 & 0 & 0 & 0\\0 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right)\\
\det U_{\mathrm{seq}}&=-1\\
U_{\mathrm{seq}}\ket0&=\left(\begin{matrix}1\\0\\0\\0\end{matrix}\right)\\
U_{\mathrm{seq}}\ket{11}&=\left(\begin{matrix}0\\0\\0\\-1\end{matrix}\right)
\end{aligned}
$$

計算基底 $ab$ の総位相は $\exp(i\pi[a+b-(a\oplus b)]/2)=(-1)^{ab}$ である。位相を取り出すための作業ビットを別に準備せず、標的を一時的な記録に使って復元できる。
:::

### 普遍性と必要なゲート数の下限

連続な局所ゲートを許す場合でも、任意の $n$ 量子ビットユニタリを実現する費用は指数的に増える。$d=2^n$ とすると、$\mathrm{SU}(d)$ は実次元 $d^2-1$ である。恒等の近傍での接空間は跡零の反 Hermite 行列であり、対角が $d-1$ 個、非対角の各対が二実数を持つことを数えればよい。各 1 量子ビットゲートを位相を除いて三実数で指定し、固定した $\mathrm{CX}$ の配置と順序を許す。$m$ 個の局所ゲートが与える積は、高々 $3m$ 次元の滑らかなパラメータ写像の像である。$3m<d^2-1$ なら、その像は全体を埋められない。局所座標で有界な微分を持つ領域を一辺 $\delta$ の立方体で覆うと、$O(\delta^{-3m})$ 個必要であり、各像は半径 $O(\delta)$ に入る。全体の $d^2-1$ 次元体積の上界は $O(\delta^{d^2-1-3m})$ で、$\delta\downarrow0$ で零になる。パラメータ領域を可算個の有界な部分へ分けても零集合の和である。有限長で配線の選択肢を変えても有限和なので、同じ次元の下界が成り立つ。

:::example[二量子ビットと五量子ビットでパラメータを数える]
二量子ビットでは $d=4$、全体の実次元は 15 である。連続パラメータを完全に使っても 1 量子ビットゲートは少なくとも 5 個必要である。五量子ビットでは $d=32$、次元は 1023 で、少なくとも 341 個必要である。この下界は $\mathrm{CX}$ の必要数を直接与えない。個々の回路に要る数でもなく、任意のユニタリを表す能力についての下界である。多項式の回路族は、全ユニタリから特定用途の構造を選んで使う。
:::

## 量子テレポーテーション

送信者の未知状態を $\ket\psi$ とし、送受信者は $\Phi^+$ の対を一つ共有する。:term[量子テレポーテーション]{en="quantum teleportation" yomi="りょうしてれぽーてーしょん"}では、未知状態と送信者の共有ビットを Bell 測定する。

:::proposition[テレポーテーションの回路恒等式]{#prop-qi-teleportation}
$\mathrm{CX}_{1\to2}$ の後に $H_1$ を適用すると初期状態 $\ket\psi_1\ket{\Phi^+}_{23}$ は

$$

\frac12\sum_{a,b=0}^1\ket{ab}_{12}X^bZ^a\ket\psi_3
\label{eq-qi-teleportation}

$$

となる。
:::

:::proof
入力を $\sum_x\psi_x\ket x$ とすると、共有対のラベル $j$ に対し $\mathrm{CX}$ 後は $\ket{x,j\oplus x,j}$ である。$H$ で第 1 ビットが $a$ になる振幅は $(-1)^{ax}/\sqrt2$ である。$b=j\oplus x$ を固定すると第 3 ビットは $\ket{x\oplus b}$ であり、位相と合わせると $X^bZ^a$ になる。
:::

各結果が $1/4$ なので、二つの古典ビットからは入力を推定できない。受信者は $Z^aX^b$ を適用すれば元に戻せる。入力が他の系と相関している場合も、式の線形性によりその相関を保つ。古典ビットが届く前の受信状態は、一様 Pauli 共役で平均されるため $I/2$ である。

### 雑音がある共有対

共有資源が $\Phi^+$ に $I\otimes P_j$ を掛けた四つの Bell 状態の混合であり、各確率を $p_j$ とする。共有対の受信側に作用する $P_j$ は、送信側の Bell 測定と可換である。分枝 $(a,b)$ の後で補正すると受信側の作用は $Z^aX^b P_j X^bZ^a=\pm P_j$ になる。大域符号は出力密度行列から消えるので、実装されたチャネルは $\sum_jp_jP_j\rho P_j^\dagger$ である。参照系と相関する入力でも同じである。共有対の誤りを、通信された状態への Pauli 雑音として読める。

共有資源 $\tau$ が Bell 基底で非対角でも、標準手順で測定ラベルを捨てた量子出力には対角成分だけが残る。$\ket{\beta_j}=(I\otimes P_j)\ket{\Phi^+}$ とし、$\tau$ の行列単位 $\ketbra{\beta_j}{\beta_k}$ を入力資源とした場合を考える。補正後の枝 $(a,b)$ では $P_j\rho P_k^\dagger$ に $1/4$ と二つの共役符号の積が掛かる。四つの $(a,b)$ に対し、$P_j=I,X,Y,Z$ の符号はそれぞれ $1,(-1)^a,(-1)^{a+b},(-1)^b$ である。異なる二つの符号列の積は平均すると零になり、同じ符号列なら平均は 1 になる。従って $p_j=\mel{\beta_j}{\tau}{\beta_j}$ を用いた Pauli 混合が得られる。

純粋入力を Bloch 球面から一様に選ぶと、入力と出力の射影重なりの平均は $p_0+(1-p_0)/3$ である。非自明な Pauli $P_j$ による枝の射影重なりは $|\mel{\psi}{P_j}{\psi}|^2=r_j^2$ となる。回転対称性から三座標の二乗の平均は等しく、その和は 1 なのでそれぞれ $1/3$ である。例えば $\tau$ が $\Phi^+$ と $\Phi^-$ の等確率混合なら、$\ket0,\ket1$ はそのまま伝わる一方で、球面での平均重なりは $2/3$ になる。

:::example[位相の混ざった Bell 資源]
$\Phi^+$ を確率 $4/5$、$\Phi^-$ を $1/5$ で共有する。入力 $\ket+$ は出力で $\ket+$ と $\ket-$ の $4/5,1/5$ の混合になるので、$X$ の平均は $3/5$ である。入力 $\ket0$ では $Z$ が $\ket0$ を固定し、出力は純粋なままである。情報を運ぶのに必要な共有対の品質は、送る状態の基底に依存して観測される。
:::

:::example[一方の古典ビットを落とす]
共有対は理想とし、送信者の $b$ だけを受け取った場合、受信者は $X^b$ を掛けられる。その後に残るのは $Z^a$ で、$a$ は公平だからチャネルは $(\rho+Z\rho Z)/2$ になる。$\ket+$ は $I/2$ に移るが、$\diag(3/4,1/4)$ は保存される。$a$ だけ受信した場合は $X$ の平均化となって役割が入れ替わる。
:::

### 測定枝を明示して補正を照合する

入力の一量子ビットを $\psi$ とし、受信した二ビットが $a,b$ なら受信系の枝は $X^b Z^a\psi/2$ である。以下では全三ビットの回路を直接掛け、各分枝の振幅から確率と規格化状態を求める。いずれの枝でも補正後は同じ入力に戻る。

:::example[計算基底状態を送る]

Alice の測定順を $(a,b)$ とし、非規格化枝 $\xi_{ab}$ を区別して書く。

$$
\begin{aligned}
\psi&=\left(\begin{matrix}0\\1\end{matrix}\right)\\
\xi_{00}&=\left(\begin{matrix}0\\\frac{1}{2}\end{matrix}\right)\\
\xi_{01}&=\left(\begin{matrix}\frac{1}{2}\\0\end{matrix}\right)\\
\xi_{10}&=\left(\begin{matrix}0\\- \frac{1}{2}\end{matrix}\right)\\
\xi_{11}&=\left(\begin{matrix}- \frac{1}{2}\\0\end{matrix}\right)\\
(p_{00},p_{01},p_{10},p_{11})&=(\frac{1}{4},\frac{1}{4},\frac{1}{4},\frac{1}{4})\\
Z^a X^b\, (2\xi_{ab})&=\left(\begin{matrix}0\\1\end{matrix}\right)
\end{aligned}
$$

受信側の位相符号を変える枝はあるが、単独の基底状態では観測できない。入力が基底だけでも測定結果は一様であり、送信者の二ビットを読むだけでは入力を推定できない。
:::

:::example[実の二つの係数を持つ信号]

Alice の測定順を $(a,b)$ とし、非規格化枝 $\xi_{ab}$ を区別して書く。

$$
\begin{aligned}
\psi&=\left(\begin{matrix}\frac{3}{5}\\\frac{4}{5}\end{matrix}\right)\\
\xi_{00}&=\left(\begin{matrix}\frac{3}{10}\\\frac{2}{5}\end{matrix}\right)\\
\xi_{01}&=\left(\begin{matrix}\frac{2}{5}\\\frac{3}{10}\end{matrix}\right)\\
\xi_{10}&=\left(\begin{matrix}\frac{3}{10}\\- \frac{2}{5}\end{matrix}\right)\\
\xi_{11}&=\left(\begin{matrix}- \frac{2}{5}\\\frac{3}{10}\end{matrix}\right)\\
(p_{00},p_{01},p_{10},p_{11})&=(\frac{1}{4},\frac{1}{4},\frac{1}{4},\frac{1}{4})\\
Z^a X^b\, (2\xi_{ab})&=\left(\begin{matrix}\frac{3}{5}\\\frac{4}{5}\end{matrix}\right)
\end{aligned}
$$

一部の枝では二係数の入れ替え、他では符号の反転が生じる。補正順で誤ると Pauli 同士の順序差に大域符号が出る。ここでは定義に従い、元の ket まで一致させた。
:::

:::example[虚数成分のある信号]

Alice の測定順を $(a,b)$ とし、非規格化枝 $\xi_{ab}$ を区別して書く。

$$
\begin{aligned}
\psi&=\left(\begin{matrix}\frac{\sqrt{2}}{2}\\\frac{\sqrt{2} i}{2}\end{matrix}\right)\\
\xi_{00}&=\left(\begin{matrix}\frac{\sqrt{2}}{4}\\\frac{\sqrt{2} i}{4}\end{matrix}\right)\\
\xi_{01}&=\left(\begin{matrix}\frac{\sqrt{2} i}{4}\\\frac{\sqrt{2}}{4}\end{matrix}\right)\\
\xi_{10}&=\left(\begin{matrix}\frac{\sqrt{2}}{4}\\- \frac{\sqrt{2} i}{4}\end{matrix}\right)\\
\xi_{11}&=\left(\begin{matrix}- \frac{\sqrt{2} i}{4}\\\frac{\sqrt{2}}{4}\end{matrix}\right)\\
(p_{00},p_{01},p_{10},p_{11})&=(\frac{1}{4},\frac{1}{4},\frac{1}{4},\frac{1}{4})\\
Z^a X^b\, (2\xi_{ab})&=\left(\begin{matrix}\frac{\sqrt{2}}{2}\\\frac{\sqrt{2} i}{2}\end{matrix}\right)
\end{aligned}
$$

位相が虚数の場合も $Z$ の符号を掛ける位置が変わるだけである。例えば $a=1,b=0$ の枝は $Y$ の負固有状態になり、補正前の $Y$ の結果は入力とは逆になる。
:::

### 一般次元のテレポーテーションと補正順序

$d$ 次元で添字を法 $d$ とし、$X\ket j=\ket{j+1}$、$Z\ket j=\omega^j\ket j$、$\omega=e^{2\pi i/d}$ とする。共有対は $\sum_j\ket{jj}/\sqrt d$ とする。送信者の入力を $\sum_x c_x\ket x$ とし、送信者の二レジスタで第 2 から第 1 を引く。その後、第一レジスタに $\ket x\mapsto d^{-1/2}\sum_a\omega^{-ax}\ket a$ の逆 Fourier 変換を掛け、二出力を $(a,b)$ として読む。固定した $b=j-x$ に対し、受信者は $\ket{j}=\ket{x+b}$ を得て位相 $\omega^{-ax}$ が付く。従って条件付きの非規格化演算は $d^{-1}X^b Z^{-a}$ であり、確率は全結果で $1/d^2$ である。受信者は $Z^a X^{-b}$ を掛け、元の状態に戻す。二次元では正負のビットが同じ剰余になるが、一般次元では符号を残す。

:::example[三準位の未知位相を受け取る]
$d=3$、入力 $(\ket0+i\ket2)/\sqrt2$、測定結果 $(a,b)=(1,2)$ とする。条件付き状態は $(\ket2+i\omega^{-2}\ket1)/\sqrt2$ である。$X^{-2}$ は添字を $2\to0,1\to2$ に戻し、$Z$ の $\omega^2$ が第二項の位相と相殺する。結果は元の状態に一致する。古典ラベルが届く前に全九分枝を平均すると、$Z^{-a}$ の平均で非対角項が消え、$X^b$ の平均で対角成分が一様になり、 $I/3$ になる。
:::

### 受信結果の記録を誤った場合

二つの測定ビットを古典通信で送り、それぞれの反転指示を $e_a,e_b$ とする。受信者が誤ラベルに基づく補正をすると、大域位相を除いて $Z^{e_a}X^{e_b}$ が残る。二反転が独立で、ともに確率 $q$ なら、結果チャネルの $I,X,Y,Z$ の重みは $((1-q)^2,q(1-q),q^2,q(1-q))$ である。共有 Bell 対が理想でも、補正を指定する古典通信が乱れると量子出力は雑音を受ける。

:::example[通信誤りの軸の違い]
位相に対応するビット $a$ だけが確率 $q=1/8$ で誤るとする。入力 $\ket+$ の出力は $\ket+$ を確率 $7/8$、$\ket-$ を $1/8$ で混ぜ、$X$ の期待値は $3/4$ になる。入力 $\ket0$ は常に同じである。逆に $b$ だけ誤れば $\ket0$ の出力で $Z$ が $3/4$ に減り、$\ket+$ は変わらない。共有対の雑音と古典チャネルの雑音を別に測らない場合、受信状態ではこのような Pauli 型の誤りが合成されて見える。
:::

## 超密度符号化

:term[超密度符号化]{en="superdense coding" yomi="ちょうみつどふごうか"}では送信者が共有 Bell 対の片方に $X^bZ^a$ を適用し、その 1 量子ビットを受信者へ送る。受信者は Bell 基底で二ビットを読む。四つの状態は、$b=0$ で $\ket{00}+(-1)^a\ket{11}$、$b=1$ で $\ket{10}+(-1)^a\ket{01}$ の規格化である。従って異なる $b$ の支持は直交し、同じ $b$ で異なる $a$ も内積が零である。

:::example[二つのビットの読み出し]
$(a,b)=(1,0)$ なら $\Phi^-$ を送り、逆 Bell 回路で $\ket{10}$ が得られる。共有対がまだ配られていないときは、その配布にも独立の通信資源が要る。
:::

### 局所演算の直交性と高次元版

$d$ 次元で $\omega=e^{2\pi i/d}$、$X\ket j=\ket{j+1\bmod d}$、$Z\ket j=\omega^j\ket j$ とする。共有対 $\ket{\Omega_d}=d^{-1/2}\sum_j\ket{jj}$ に送信側で $X^bZ^a$ を掛けると
$$
\ket{\Omega_{ab}}=\frac1{\sqrt d}\sum_j\omega^{aj}\ket{j+b,j}
$$
である。異なる $b$ は積基底の支持が重ならず、同じ $b$ の異なる $a$ は幾何級数の和が零なので直交する。受け取った順に第 1 レジスタから第 2 レジスタを引けば $\ket b\otimes d^{-1/2}\sum_j\omega^{aj}\ket j$ に分離し、第 2 レジスタへ逆 Fourier 変換を掛けて $a$ を読む。従って一回の $d$ 次元系の送信で $d^2$ 通りのラベルを区別できる。受信者がもともと持つ半分の縮約は、送信前のどの局所ユニタリでも $I/d$ のままである。

:::example[三準位での復号]
$d=3,(a,b)=(2,1)$ なら送信状態は $(\ket{10}+\omega^2\ket{21}+\omega\ket{02})/\sqrt3$ である。引き算で $\ket1$ が第 1 レジスタに出て、第 2 レジスタに $(\ket0+\omega^2\ket1+\omega\ket2)/\sqrt3$ が残る。これは $\ket2$ の Fourier 像だから逆変換で $\ket2$ を得る。受信状態の全空間が九次元なので、直交した九通りを超える完全識別は不可能である。
:::

:::example[最大もつれから外れた共有対]
$\sqrt p\ket{00}+\sqrt{1-p}\ket{11}$ に同じ四つの Pauli 符号を等確率で掛ける。各 $b$ に対し $a=0,1$ の内積は $2p-1$ となる。例えば $p=3/4$ では重なり $1/2$ なので、$b$ を確実に識別してから :ref[thm-qi-helstrom] で $a$ を最適二択する成功率は $(1+\sqrt3/2)/2$ である。各 $b$ の入力状態は対応する直交部分空間にあるため、最初にその空間へ射影しても、任意の POVM で得る統計を再現できる。$b$ が既知なら対応する二択に集中できるので、この値は等事前四択でも最適である。
:::

### 復号行列で位相と値のラベルを区別する

Bell の四状態に対し、量子回路の復号 $\mathrm{CX}$ の後で $H$ を触れる位置は相対位相を持つ第 1 ビットである。以下では送信側の誤りや準備の偏りが復号後の確率にどう現れるかを、二量子ビットの振幅で求める。

:::example[位相だけを送る]

符号ラベルは $(a,b)=(1,0)$ とする。共有対から順に二量子ビットの係数を並べる。

$$
\begin{aligned}
\psi_{\mathrm{shared}}&=\left(\begin{matrix}\frac{\sqrt{2}}{2}\\0\\0\\\frac{\sqrt{2}}{2}\end{matrix}\right)\\
\psi_{\mathrm{sent}}&=(X^bZ^a\otimes I)\psi_{\mathrm{shared}}=\left(\begin{matrix}\frac{\sqrt{2}}{2}\\0\\0\\- \frac{\sqrt{2}}{2}\end{matrix}\right)\\
\psi_{\mathrm{CX}}&=\left(\begin{matrix}\frac{\sqrt{2}}{2}\\0\\- \frac{\sqrt{2}}{2}\\0\end{matrix}\right)\\
\psi_{\mathrm{read}}&=\left(\begin{matrix}0\\0\\1\\0\end{matrix}\right)\\
(p_{00},p_{01},p_{10},p_{11})&=(0,0,1,0)
\end{aligned}
$$

入力の二つの計算基底振幅の間の負号が、復号後の第 1 ビット 1 になる。第 2 ビットは同じラベルを二度排他的論理和するので 0 である。
:::

:::example[値と位相を両方送る]

符号ラベルは $(a,b)=(1,1)$ とする。共有対から順に二量子ビットの係数を並べる。

$$
\begin{aligned}
\psi_{\mathrm{shared}}&=\left(\begin{matrix}\frac{\sqrt{2}}{2}\\0\\0\\\frac{\sqrt{2}}{2}\end{matrix}\right)\\
\psi_{\mathrm{sent}}&=(X^bZ^a\otimes I)\psi_{\mathrm{shared}}=\left(\begin{matrix}0\\- \frac{\sqrt{2}}{2}\\\frac{\sqrt{2}}{2}\\0\end{matrix}\right)\\
\psi_{\mathrm{CX}}&=\left(\begin{matrix}0\\- \frac{\sqrt{2}}{2}\\0\\\frac{\sqrt{2}}{2}\end{matrix}\right)\\
\psi_{\mathrm{read}}&=\left(\begin{matrix}0\\0\\0\\-1\end{matrix}\right)\\
(p_{00},p_{01},p_{10},p_{11})&=(0,0,0,1)
\end{aligned}
$$

第 2 レジスタは第 1 と異なるラベルであることから 1 になる。第 1 は二基底の差の組合せから 1 を得る。符号の作用順序を変えると、結果の ket 全体の符号だけが変わる。
:::

:::example[偏った共有対を理想 Bell 復号へ送る]

符号ラベルは $(a,b)=(0,0)$ とする。共有対から順に二量子ビットの係数を並べる。

$$
\begin{aligned}
\psi_{\mathrm{shared}}&=\left(\begin{matrix}\frac{\sqrt{3}}{2}\\0\\0\\\frac{1}{2}\end{matrix}\right)\\
\psi_{\mathrm{sent}}&=(X^bZ^a\otimes I)\psi_{\mathrm{shared}}=\left(\begin{matrix}\frac{\sqrt{3}}{2}\\0\\0\\\frac{1}{2}\end{matrix}\right)\\
\psi_{\mathrm{CX}}&=\left(\begin{matrix}\frac{\sqrt{3}}{2}\\0\\\frac{1}{2}\\0\end{matrix}\right)\\
\psi_{\mathrm{read}}&=\left(\begin{matrix}\frac{\sqrt{2}}{4} + \frac{\sqrt{6}}{4}\\0\\- \frac{\sqrt{2}}{4} + \frac{\sqrt{6}}{4}\\0\end{matrix}\right)\\
(p_{00},p_{01},p_{10},p_{11})&=(\frac{\sqrt{3}}{4} + \frac{1}{2},0,\frac{1}{2} - \frac{\sqrt{3}}{4},0)
\end{aligned}
$$

偶奇は確実に 0 だが、位相ラベルは完全には識別できない。ここでは理想資源用の復号を固定している。共有対の非最大性を許すときには、Bell の完全な四択を根拠に確率 1 としてはいけない。
:::

### Weyl 演算子の直交性とラベルの上限

$d$ 次元の最大もつれ状態に片側から $U$ を掛けると、二つの符号状態の内積は $\Tr(U^\dagger V)/d$ である。このことは $\sum_{ij}\bra{ii}(U^\dagger V\otimes I)\ket{jj}=\sum_i\mel i{U^\dagger V}i$ から分かる。従って Hilbert–Schmidt 内積で直交する $d^2$ 個のユニタリを選べば完全なラベルを作れる。一方、誤りなしの識別が可能な純粋信号は互いに直交する必要があるので、合成の受信空間が $d^2$ 次元である場合はそれ以上のラベルを送れない。$X^bZ^a$ の直交性では、$b$ が異なると差の積に非零のシフトが残って跡が零になる。同じ $b$ で $a$ が異なると、跡は 1 の根の和になって零である。

:::example[受信側だけの観測に符号が現れないことを成分で示す]
$d=2$ の四符号で、受信者が元から持つ片側はいつも $I/2$ である。一方、送信した量子ビットが届いても、二つを計算基底だけで測ると $b$ は二ビットの排他的論理和として分かるが、$a$ は分からない。例えば $b=0$ の $\Phi^\pm$ はどちらも $00,11$ を半々で返す。Bell 測定を使うには先に $\mathrm{CX}$ で偶奇を第 2 ビットに移し、$H$ で位相を第 1 ビットへ移す。この順序なら $(a,b)$ が両方出る。
:::

## 測定型量子計算の概観

:term[測定型量子計算]{en="measurement-based quantum computation" yomi="そくていがたりょうしけいさん"}では、もつれた資源状態を作っておき、局所測定の順序と基底を選んで演算を進める。

一段を計算する。$\ket\psi=u\ket0+v\ket1$ と $\ket+$ に制御 $Z$ を作用させると $u\ket0\ket++v\ket1\ket-$ である。第 1 ビットを $(\ket0+(-1)^m e^{i\alpha}\ket1)/\sqrt2$ に射影した第 2 ビットの非規格化状態は

$$

\frac1{\sqrt2}H Z^m\diag(1,e^{-i\alpha})\ket\psi
=\frac1{\sqrt2}X^m H\diag(1,e^{-i\alpha})\ket\psi

$$

である。従って各結果は $1/2$ であり、$X^m$ を追跡すれば所望の回転を得る。次段に未補正の $X$ が入ると、$R_z(\alpha)X=XR_z(-\alpha)$ により測定角の符号を切り替える。この古典的な依存関係を :term[適応測定]{en="adaptive measurement" yomi="てきおうそくてい"}と呼ぶ。


:::example[出力側を結ぶ二本の経路]
データ $\ket u,\ket v$ に補助 $\ket+,\ket+$ を一つずつ付け、各データと対応補助を制御 $Z$ で結ぶ。さらに二つの補助の間にも制御 $Z$ を掛ける。後者はデータ測定と可換なので、測定後に作用したと考えてよい。両データを $X$ 基底で測ると、補助への条件付き入力は

$$

\mathrm{CZ}(X^mH\otimes X^nH)\ket{u,v}.
$$

$\mathrm{CZ}$ と各 $X$ の交換で、反対側に $Z$ が付く。従って既知の Pauli 補正を追跡すれば、出力は $\mathrm{CZ}(H\otimes H)\ket{u,v}$ である。例えば $u=v=0$ で測定結果 $m=n=0$ なら、出力は $(\ket{00}+\ket{01}+\ket{10}-\ket{11})/2$ である。係数行列の行列式は $-1/2$ なので、測定と二本の経路の結合によりもつれを生成したことが分かる。
:::

### Pauli 補正を測定基底へ送る

理想ゲートを $J(\alpha)=HR_z(-\alpha)$ とし、測定の入力が理想状態の $X^xZ^z$ 補正付きだとする。測定角を $\alpha'=(-1)^x\alpha$ と選ぶと、実際の出力は大域位相を除いて
$$
X^m J(\alpha') X^x Z^z\ket\psi
=X^{m\oplus z}Z^x J(\alpha)\ket\psi
$$
になる。$R_z(-\alpha')X^x=X^xR_z(-\alpha)$ と、$H XH=Z,H ZH=X$ を用いた。従って記録 $(x,z)$ は $(m\oplus z,x)$ に更新するだけでよい。この規則は、測定のランダムな結果を装置上の量子的な訂正ゲートに毎回直す必要がないことを示す。

:::example[角度が切り替わる二段]
最初に $\alpha_1=0$ で結果 $m_1=1$ を得ると、理想の $H$ の後に $X$ が残る。次に理想ゲートを $J(\pi/4)$ としたければ $\alpha_2'=-\pi/4$ を選ぶ。$m_2=0$ なら出力の余分な Pauli は $Z$、$m_2=1$ なら $XZ$ である。入力 $\ket0$ に対する理想列は $HR_z(-\pi/4)\ket+$ で、理想の $Z$ の期待値は $\cos(\pi/4)=1/\sqrt2$ である。$m_2=1$ の枝では $X$ が $Z$ の符号を変えるので、測定結果の記録の符号を反転すれば両分枝で同じ分布に戻る。
:::

:::example[測定だけで任意の一量子ビット演算へ進む]
$J(0)J(\alpha)=R_z(-\alpha)$、$J(\alpha)J(0)=R_x(-\alpha)$ となる。従って $R_z(\pi/2)R_x(\pi/4)R_z(-\pi/2)=R_y(\pi/4)$ は、この分解を使えば六段の鎖で構成できる。$R_z(\pi/2)$ が $X$ を $Y$ に共役で移すためである。隣り合う $J(0)$ の積が $I$ になる箇所は省略できる。
:::

### ゲートの列を資源グラフへ写す

単一線の鎖で $J(\alpha)$ を実装し、二線の出力を結んで $\mathrm{CZ}(H\otimes H)$ を実装できる。後者の直前に両線で $J(0)=H$ を一度ずつ実装すれば、理想演算の積は $\mathrm{CZ}$ になる。実装途中の Pauli 補正は、単一線では測定角と結果ラベルに反映させる。二線を結ぶ段では $X$ が $\mathrm{CZ}$ を通過すると反対側に $Z$ を残すので、その値も補正記録に加える。Euler 分解と $\mathrm{CZ}$ による $\mathrm{CX}$ の構成から、これらの段で任意の量子回路を置き換えられる。

各段で新しい頂点を用意し、出力を次の段の入力と同じ頂点にする。辺に対応する全 $\mathrm{CZ}$ は互いに可換である。まだ測定していない頂点だけを結ぶ後段の辺は、前段の測定演算とも可換である。そのため、各辺の $\mathrm{CZ}$ はその両端の測定よりも前にまとめて作用させてよい。残る操作は一量子ビット測定と古典的な補正記録である。ここでは回路に応じて辺を選んだ。固定した二次元格子を使う場合は配線の構成が別途必要であり、具体的な測定配置は :cite[preskill-notes]{loc="Cluster states（補講ノート）"} を参照する。

### 枝の Pauli を残したまま回転角を調整する

測定角 $\alpha$ の一段に対し、非規格化枝は $X^mH\diag(1,e^{-i\alpha})\psi/\sqrt2$ である。下の例では ket の複素係数まで求めるため、$J_\alpha=H\diag(1,e^{-i\alpha})=e^{-i\alpha/2}J(\alpha)$ を用いる。両規約の差は一段ごとの大域位相であり、補正後の確率には影響しない。既知の $X$ が前段に残っていれば次の角を反転し、$Z$ が残っていれば出力の記録の $X$ ビットを反転する。

:::example[ゼロ角で二回進める]

二段の理想角を $\alpha=0$、$\beta=0$ とし、実際の測定結果は $(m_1,m_2)=(1,0)$ とする。各枝は規格化して書く。

$$
\begin{aligned}
\psi&=\left(\begin{matrix}\frac{\sqrt{3}}{2}\\\frac{1}{2}\end{matrix}\right)\\
\psi_1&=\left(\begin{matrix}- \frac{\sqrt{2}}{4} + \frac{\sqrt{6}}{4}\\\frac{\sqrt{2}}{4} + \frac{\sqrt{6}}{4}\end{matrix}\right)\\
\beta_{\mathrm{adapt}}&=(-1)^{m_1}\beta=0\\
\psi_2&=\left(\begin{matrix}\frac{\sqrt{3}}{2}\\- \frac{1}{2}\end{matrix}\right)\\
X^{m_2}Z^{m_1}&=\left(\begin{matrix}1 & 0\\0 & -1\end{matrix}\right)\\
\psi_{\mathrm{ideal}}&=\left(\begin{matrix}\frac{\sqrt{3}}{2}\\\frac{1}{2}\end{matrix}\right)\\
|\braket{\psi_{\mathrm{ideal}}}{Z^{m_1}X^{m_2}\psi_2}|^2&=1
\end{aligned}
$$

理想では $HH=I$ で入力が移るだけである。第一測定の $X$ は第二の $H$ を通ると $Z$ に移り、第二測定の結果 0 は新しい $X$ を付けない。
:::

:::example[角度を反転しないと向きが変わる]

二段の理想角を $\alpha=0$、$\beta=\frac{\pi}{2}$ とし、実際の測定結果は $(m_1,m_2)=(1,0)$ とする。各枝は規格化して書く。

$$
\begin{aligned}
\psi&=\left(\begin{matrix}\frac{\sqrt{3}}{2}\\\frac{1}{2}\end{matrix}\right)\\
\psi_1&=\left(\begin{matrix}- \frac{\sqrt{2}}{4} + \frac{\sqrt{6}}{4}\\\frac{\sqrt{2}}{4} + \frac{\sqrt{6}}{4}\end{matrix}\right)\\
\beta_{\mathrm{adapt}}&=(-1)^{m_1}\beta=- \frac{\pi}{2}\\
\psi_2&=\left(\begin{matrix}\frac{\left(-1 + \sqrt{3} i\right) \left(1 - i\right)}{4}\\\frac{\left(1 - i\right) \left(\sqrt{3} - i\right)}{4}\end{matrix}\right)\\
X^{m_2}Z^{m_1}&=\left(\begin{matrix}1 & 0\\0 & -1\end{matrix}\right)\\
\psi_{\mathrm{ideal}}&=\left(\begin{matrix}\frac{\left(1 - i\right) \left(\sqrt{3} + i\right)}{4}\\\frac{\left(1 - i\right) \left(1 + \sqrt{3} i\right)}{4}\end{matrix}\right)\\
|\braket{\psi_{\mathrm{ideal}}}{Z^{m_1}X^{m_2}\psi_2}|^2&=1
\end{aligned}
$$

第一の結果が 1 なので、第二角を $-\pi/2$ にする。反転しない場合は相対位相が逆になる。ここでは最後の $Z$ 補正によって同じ理想回転に戻る。
:::

:::example[両測定の結果が負になる]

二段の理想角を $\alpha=\frac{\pi}{3}$、$\beta=\frac{\pi}{2}$ とし、実際の測定結果は $(m_1,m_2)=(1,1)$ とする。各枝は規格化して書く。

$$
\begin{aligned}
\psi&=\left(\begin{matrix}\frac{\sqrt{3}}{2}\\\frac{1}{2}\end{matrix}\right)\\
\psi_1&=\left(\begin{matrix}- \frac{\sqrt{2}}{8} + \frac{\sqrt{6} \left(2 + i\right)}{8}\\\frac{\sqrt{2}}{8} + \frac{\sqrt{6} \left(2 - i\right)}{8}\end{matrix}\right)\\
\beta_{\mathrm{adapt}}&=(-1)^{m_1}\beta=- \frac{\pi}{2}\\
\psi_2&=\left(\begin{matrix}\frac{\left(1 - i\right) \left(\sqrt{3} - i\right)}{8}\\\frac{\left(-1 + 3 \sqrt{3} i\right) \left(1 - i\right)}{8}\end{matrix}\right)\\
X^{m_2}Z^{m_1}&=\left(\begin{matrix}0 & -1\\1 & 0\end{matrix}\right)\\
\psi_{\mathrm{ideal}}&=\left(\begin{matrix}\frac{\left(1 - i\right) \left(3 \sqrt{3} + i\right)}{8}\\\frac{\left(1 - i\right) \left(1 + \sqrt{3} i\right)}{8}\end{matrix}\right)\\
|\braket{\psi_{\mathrm{ideal}}}{Z^{m_1}X^{m_2}\psi_2}|^2&=1
\end{aligned}
$$

最後の副作用に $X$ と $Z$ の両方が入る。Pauli の順序を固定して追えば、最終測定の符号を必要な軸だけで反転できる。
:::

### グラフ状態の検査を測定後も追跡する

無向単純グラフ $G$ の頂点に $\ket+$ を準備し、辺ごとに $\mathrm{CZ}$ を掛けた状態を :term[グラフ状態]{en="graph state" yomi="ぐらふじょうたい"}という。初期には各 $X_v$ が状態を固定する。$\mathrm{CZ}_{uv}$ で共役すると、$X_v$ に $Z_u$ が掛かる。計算基底で $X_v$ が $v$ を反転し、制御位相の前後差が隣接する $u$ の値だからである。従って最終状態は $X_v\prod_{u\in N(v)}Z_u$ を全頂点で $+1$ にする。これらは元の $X_v$ の共役なので可換であり、独立である。頂点 $v$ を $Z$ で測って値 $s$ を得ると、$\mathrm{CZ}$ の辺因子のうち $v$ に触れる分は各隣接頂点への $Z^s$ になる。残りの辺因子は同じなので、$v$ を削除したグラフ状態へ隣接 $Z$ の補正が付く。この計算により、グラフの一部を使わない場合にどの測定で取り除けるか分かる。

:::example[鎖の中央を一つ測る]
五頂点の鎖 $1-2-3-4-5$ を用意する。頂点 3 を $Z$ で測り $s=1$ なら、残る辺は $1-2$ と $4-5$ で、$Z$ が 2 と 4 へ付く。二つの部分グラフは独立なグラフ状態である。従って結果 $s$ が分かれば、両方の局所補正を行って二本の独立な線に分けられる。$s$ を忘れると両側に同じ符号が付くため、条件付きでは積でも、無条件には古典相関が残りうる。
:::

## 演習問題

::::exercise{level=basic}
$HXH$ を求めよ。

:::solution
$HX\ket+=H\ket+=\ket0$ と $HX\ket-= -H\ket-= -\ket1$ であるから、両計算基底上の作用を合わせて $HXH=Z$ である。
:::
::::

::::exercise{level=basic}
$R_z(\pi/2)\ket+$ の Bloch ベクトルを求めよ。

:::solution
:ref[prop-qi-bloch-rotation] で $Z$ 軸の周りに $X$ 軸を正に $\pi/2$ 回すので $(0,1,0)$ となる。
:::
::::

::::exercise{level=basic}
$\mathrm{CX}(\ket+\otimes\ket-)$ を求めよ。

:::solution
第 2 系は $X$ の固有値 $-1$ の状態である。第 1 ビットが 1 の成分だけに $-1$ が掛かり、結果は $\ket-\otimes\ket-$ である。
:::
::::

::::exercise{level=basic}
Toffoli の入力 $110$ の出力を求めよ。

:::solution
制御が両方 1 なので第 3 ビットを反転し、$111$ となる。
:::
::::

::::exercise{level=advanced}
$C(Z)=(I\otimes H)\mathrm{CX}(I\otimes H)$ を示せ。

:::solution
制御 0 では標的に $HH=I$ が働き、制御 1 では $HXH=Z$ が働く。基底上の指定が線形演算を定める。
:::
::::

::::exercise{level=advanced}
式 :eqref[eq-qi-teleportation] で $(a,b)=(1,1)$ のとき、未補正の受信状態が $\ket\psi=\ket+$ に対して何になるか求めよ。

:::solution
$Z\ket+=\ket-$、$X\ket-=-\ket-$ なので $-\ket-$ である。順序を逆にした補正 $ZX$ を施すと $\ket+$ に戻る。
:::
::::

::::exercise{level=advanced}
超密度符号化の四状態の等確率混合を求めよ。

:::solution
四状態は正規直交基底なので、射影の和は $I_4$ である。平均は $I_4/4$ となる。
:::
::::

::::exercise{level=advanced}
測定型計算の一段で $\alpha=0$、入力 $\ket1$、結果 $m=1$ の条件付き出力を求めよ。

:::solution
出力は $XH\ket1=X\ket-=-\ket-$ である。$X$ で補正すると $\ket-$ となる。
:::
::::

## 確認テスト

::::quiz
$\mathrm{CX}$ を $\ket+\ket0$ に作用させた結果はどれか。

- [x] Bell 状態
- [ ] $\ket+\ket+$

:::explanation
制御 0 の枝では $\ket{00}$、制御 1 の枝では $\ket{11}$ が出る。二枝の振幅は同じである。
:::
::::

::::quiz{answer="0.25" tol="1e-6"}
標準テレポーテーションで Alice が $00$ を得る確率はいくらか。

:::explanation
式 :eqref[eq-qi-teleportation] の当該枝は $\ket\psi/2$ なので、ノルムの二乗は $1/4$ である。
:::
::::

::::quiz{answer="2" tol="1e-6"}
標準超密度符号で共有対を配布した後、1 量子ビットの送信で何古典ビットを伝えるか。

:::explanation
四つの直交状態を区別するので $\log_2 4=2$ ビットを伝える。
:::
::::

::::quiz
制御ゲートを回路へ移すときに正しい記述を選べ。

- [x] $U$ に掛けた大域位相が制御側の相対位相に現れる。
- [ ] 任意の大域位相を制御版でも必ず消せる。

:::explanation
制御 1 の枝だけに位相が掛かるので、制御が重ね合わせの場合は二枝の干渉が変わる。
:::
::::

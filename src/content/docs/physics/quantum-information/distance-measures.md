---
title: 量子状態の距離と識別
description: 状態とチャネルの識別限界を導き、有限回の測定で状態を推定する方法を扱う。
notation:
  - symbol: 'D(\rho,\sigma)'
    meaning: トレース距離
  - symbol: 'F(\rho,\sigma)'
    meaning: 二乗しない忠実度
  - symbol: '\|\Phi\|_\diamond'
    meaning: 写像のダイヤモンドノルム
---

## この章の目標

- 誤判定率を最小にする測定を計算する。
- 忠実度を純粋化間の内積と結び付ける。
- 未知の非直交状態の完全識別と複製が不可能であることを示す。
- 補助系を使うチャネル識別と、反復測定による状態推定を区別する。

## トレース距離と Helstrom 測定

状態を一つ受け取り、準備者が用いた二つの候補のどちらかを答える。測定の結果を二種類の答えへ振り分ければ、一般の測定も二効果 $E,I-E$ にまとめられる。

:::definition[トレース距離]{#def-qi-trace-distance}
$\norm A_1=\Tr\sqrt{A^\dagger A}$ をトレースノルムとし、:term[トレース距離]{en="trace distance" yomi="とれーすきょり"}を $D(\rho,\sigma)=\norm{\rho-\sigma}_1/2$ と定める。
:::

:::theorem[Helstrom の公式]{#thm-qi-helstrom}
事前確率 $p,1-p$ の二状態 $\rho,\sigma$ の最適成功率は $(1+\norm{p\rho-(1-p)\sigma}_1)/2$ である。
:::

:::proof
$\Delta=p\rho-(1-p)\sigma$ とおくと、$\rho$ と答える効果を $E$ とした成功率は $1-p+\Tr(E\Delta)$ である。固有基底で $\Tr(E\Delta)=\sum_j\lambda_jE_{jj}$ かつ $0\le E_{jj}\le1$ なので、正固有値の和が上限である。その固有空間への射影で上限を達成する。正固有値の和は $(\norm\Delta_1+\Tr\Delta)/2$ であり、$\Tr\Delta=2p-1$ を代入する。
:::

等確率なら成功率は $(1+D)/2$ である。さらに $\Delta=\rho-\sigma$ とすると $\Tr\Delta=0$ なので、$D=\max_{0\le E\le I}\Tr(E\Delta)$ である。すべての測定結果について確率差が正となる効果を足すと一つの $E$ になる。従って測定後の古典分布の全変動距離の最大値も $D$ である。双対表示と三角不等式 $\Tr(E(\rho-\tau))=\Tr(E(\rho-\sigma))+\Tr(E(\sigma-\tau))$ により $D$ は距離である。トレース保存正写像の随伴は $0\le E\le I$ を保つので、同じ最大化から状態距離はその写像で増加しない。

:::example[対角状態と事前確率]
$\rho=\diag(3/4,1/4)$、$\sigma=\diag(1/4,3/4)$ を等確率で受け取れば成功率は $3/4$ である。$p=4/5$ に変更すると $\Delta=\diag(11/20,1/20)\ge0$ となるため、測定せず常に $\rho$ と答えると最適である。
:::

:::example[球上の弦長]
Bloch 表示の差は $(\boldsymbol r-\boldsymbol s)\cdot\boldsymbol\sigma/2$ で、固有値は $\pm\norm{\boldsymbol r-\boldsymbol s}/2$ である。従って $D=\norm{\boldsymbol r-\boldsymbol s}/2$ であり、$\ket0,\ket+$ では $D=1/\sqrt2$ である。
:::

### 偏りを含む二つの純粋信号

事前確率の差があると、識別装置は等事前確率の場合とは別の方向を測る。$\Delta=p\ketbra\psi\psi-(1-p)\ketbra\phi\phi$ を、$\ket\psi$ とそれに直交するベクトルを基底にして書く。位相を選んで $\ket\phi=c\ket\psi+\sqrt{1-c^2}\ket{\psi^\perp}$ とできる。$0<p<1$、$0\le c<1$ の場合、$\Tr\Delta=2p-1$、$\det\Delta=-p(1-p)(1-c^2)$ である。固有値は正負一つずつなので、トレースノルムは二根の差 $\sqrt{1-4p(1-p)c^2}$ である。$p=0,1$ や $c=1$ でも連続性から同じ式を使える。したがって
$$
P_{\mathrm{succ}}=\frac{1+\sqrt{1-4p(1-p)c^2}}{2}.
$$
同一信号では事前確率の大きい方を常に答える確率 $\max(p,1-p)$ に戻る。

:::example[偏りが最適測定軸を変える]
$\ket\psi=\ket0,\ket\phi=\ket+$ とする。$\Delta=\{(2p-1)I+pZ-(1-p)X\}/2$ なので、正固有値の効果は軸 $(-(1-p),0,p)/\sqrt{p^2+(1-p)^2}$ の $+1$ 射影に対応する。$p=3/4$ の場合、軸は $(-1,0,3)/\sqrt{10}$ である。この測定は $\ket0$ を正しく答える条件付き確率 $(1+3/\sqrt{10})/2$、$\ket+$ を正しく答える確率 $(1+1/\sqrt{10})/2$ を与える。事前確率で重みを付けると $(1+\sqrt{10}/4)/2$ となる。$p=1/4$ へ反転すると、軸の $X$ 成分と $Z$ 成分の絶対値も入れ替わる。
:::

:::example[測定するより頻度で答える場合]
$\rho=\diag(3/5,2/5),\sigma=\diag(1/5,4/5)$、$p=3/4$ とする。$\Delta=\diag(2/5,1/10)$ は正なので、どちらの基底結果が出ても $\rho$ の重み付き尤度の方が大きい。例えば第 2 成分では $p(2/5)=3/10$、$(1-p)(4/5)=1/5$ である。したがって測定装置が状態の差を検出できても、二択の正答率は常に $\rho$ と答えた場合と同じ $3/4$ である。$p=1/2$ なら $\Delta$ の二根は $\pm1/5$ になり、基底ラベルで答えを切り替えると $7/10$ を得る。
:::

### 状態距離、測定統計、連続した変換

トレース距離はノルムに由来するから、対称性、三角不等式、距離零と一致の同値性を満たす。測定で得られる古典分布の全変動距離の最大値でもある。任意の古典事象 $S$ は効果 $\sum_{y\in S}E_y$ に対応し、その確率差は Helstrom の変分式により高々 $D$ である。逆に $\rho-\sigma$ の正固有空間への射影とその補射影の二結果を取れば $D$ を達成する。したがって距離 $\eps$ は「どの事象を記録しても、その確率の差は $\eps$ 以下」という意味を持つ。

写像 $\Phi$ が正かつトレース保存であれば、随伴 $\Phi^\dagger$ は効果を効果に送る。実際、$\Tr[\Phi^\dagger(E)\ketbra vv]=\Tr[E\Phi(\ketbra vv)]\ge0$ であり、$\Phi^\dagger(I)=I$ なので $I-\Phi^\dagger(E)\ge0$ も同様である。このため $D(\Phi(\rho),\Phi(\sigma))\le D(\rho,\sigma)$ が成り立つ。参照系がある場合にも同じ議論を使うためには、$\id\otimes\Phi$ が正である必要がある。

:::example[情報を失う測定と保存する測定]
$\ket+$ と $\ket-$ は直交するので距離 1 である。両者から $Z$ 測定のラベルだけを取り出すと、分布は $(1/2,1/2)$ で一致し距離 0 になる。$X$ 測定なら分布は $(1,0),(0,1)$ で距離 1 を保つ。さらに各記録の符号を確率 $q$ で反転させると、二分布は $(1-q,q)$ と $(q,1-q)$ となり距離は $\abs{1-2q}$ になる。$q=1/4$ なら最適判定の成功率は $3/4$ である。
:::

### 条件付き正答率から Helstrom 値を再構成する

$\Delta$ の異なる固有値が $\lambda_+>0>\lambda_-$ の二つだけなら、$E_0=(\Delta-\lambda_-I)/(\lambda_+-\lambda_-)$ は正部分への射影である。実際、固有ベクトルごとに適用するとそれぞれ 1、0 を与える。零固有空間がある場合は、その直交補に式を適用し、零固有空間上の効果は自由に選べる。以下の二次元の例ではこの式で条件付き正答率を求め、事前確率で平均する。

:::example[等確率でコヒーレンスの符号を読む]

第一仮説は $\rho_0=\ketbra00$ で事前確率 $p=\frac{1}{2}$、第二仮説は $\rho_1$ で事前確率 $1-p$ である。$q_i$ を入力 $i$ での正答率とする。

$$
\begin{aligned}
\rho_1&=\left(\begin{matrix}\frac{9}{10} & - \frac{3}{10}\\- \frac{3}{10} & \frac{1}{10}\end{matrix}\right)\\
\Delta&=\frac{1}{2}\ketbra00-\left(1-\frac{1}{2}\right)\rho_1=\left(\begin{matrix}\frac{1}{20} & \frac{3}{20}\\\frac{3}{20} & - \frac{1}{20}\end{matrix}\right)\\
\lambda_+-\lambda_-&=\sqrt{(\Tr\Delta)^2-4\det\Delta}=\frac{\sqrt{10}}{10}\\
E_0&=(\Delta-\lambda_- I)/(\lambda_+-\lambda_-)=\left(\begin{matrix}\frac{\sqrt{10}}{20} + \frac{1}{2} & \frac{3 \sqrt{10}}{20}\\\frac{3 \sqrt{10}}{20} & \frac{1}{2} - \frac{\sqrt{10}}{20}\end{matrix}\right)\\
q_0&=\Tr(E_0\ketbra00)=\frac{\sqrt{10}}{20} + \frac{1}{2}\\
q_1&=\Tr[(I-E_0)\rho_1]=\frac{\sqrt{10}}{20} + \frac{1}{2}\\
P_{\mathrm{succ}}&=pq_0+(1-p)q_1=\frac{\sqrt{10}}{20} + \frac{1}{2}
\end{aligned}
$$

二状態はともに $Z$ の正方向に近いので、差を読むには $X$ を主成分とする軸が必要である。この軸に射影した二つの Bloch 成分は反対符号になり、二条件の正答率は一致する。
:::

:::example[混合信号と純粋信号の偏りを補う]

第一仮説は $\rho_0=\ketbra00$ で事前確率 $p=\frac{2}{3}$、第二仮説は $\rho_1$ で事前確率 $1-p$ である。$q_i$ を入力 $i$ での正答率とする。

$$
\begin{aligned}
\rho_1&=\left(\begin{matrix}\frac{1}{2} & \frac{1}{4}\\\frac{1}{4} & \frac{1}{2}\end{matrix}\right)\\
\Delta&=\frac{2}{3}\ketbra00-\left(1-\frac{2}{3}\right)\rho_1=\left(\begin{matrix}\frac{1}{2} & - \frac{1}{12}\\- \frac{1}{12} & - \frac{1}{6}\end{matrix}\right)\\
\lambda_+-\lambda_-&=\sqrt{(\Tr\Delta)^2-4\det\Delta}=\frac{\sqrt{17}}{6}\\
E_0&=(\Delta-\lambda_- I)/(\lambda_+-\lambda_-)=\left(\begin{matrix}\frac{2 \sqrt{17}}{17} + \frac{1}{2} & - \frac{\sqrt{17}}{34}\\- \frac{\sqrt{17}}{34} & \frac{1}{2} - \frac{2 \sqrt{17}}{17}\end{matrix}\right)\\
q_0&=\Tr(E_0\ketbra00)=\frac{2 \sqrt{17}}{17} + \frac{1}{2}\\
q_1&=\Tr[(I-E_0)\rho_1]=\frac{\sqrt{17}}{68} + \frac{1}{2}\\
P_{\mathrm{succ}}&=pq_0+(1-p)q_1=\frac{\sqrt{17}}{12} + \frac{1}{2}
\end{aligned}
$$

混合信号が少数でも、差に負固有値があるため、常に多数ラベルを返すよりも改善できる。$q_1$ は混合候補に対する条件付き率であり、全試行に対する寄与にはさらに $1-p$ を掛ける。
:::

:::example[少数派に大きな分離があるとき]

第一仮説は $\rho_0=\ketbra00$ で事前確率 $p=\frac{3}{10}$、第二仮説は $\rho_1$ で事前確率 $1-p$ である。$q_i$ を入力 $i$ での正答率とする。

$$
\begin{aligned}
\rho_1&=\left(\begin{matrix}\frac{1}{10} & 0\\0 & \frac{9}{10}\end{matrix}\right)\\
\Delta&=\frac{3}{10}\ketbra00-\left(1-\frac{3}{10}\right)\rho_1=\left(\begin{matrix}\frac{23}{100} & 0\\0 & - \frac{63}{100}\end{matrix}\right)\\
\lambda_+-\lambda_-&=\sqrt{(\Tr\Delta)^2-4\det\Delta}=\frac{43}{50}\\
E_0&=(\Delta-\lambda_- I)/(\lambda_+-\lambda_-)=\left(\begin{matrix}1 & 0\\0 & 0\end{matrix}\right)\\
q_0&=\Tr(E_0\ketbra00)=1\\
q_1&=\Tr[(I-E_0)\rho_1]=\frac{9}{10}\\
P_{\mathrm{succ}}&=pq_0+(1-p)q_1=\frac{93}{100}
\end{aligned}
$$

混合状態には $0$ の成分が残る。基底の読み取りが理想の場合でも、誤りはその成分に由来する。測定ノイズと信号自体の重なりは、条件付き率から分けて考える。
:::

### 三つ以上の仮説の最小誤り識別

状態 $\rho_j$ が事前確率 $p_j$ で与えられるとき、ラベル $j$ と答える効果を $E_j$ とすれば成功率は $\sum_j\Tr(p_j\rho_jE_j)$ である。この量を $E_j\ge0$、$\sum_j E_j=I$ のもとで最大化する。二状態のときには一つの効果を消去できるが、三つ以上の場合は全効果の正性を同時に課す。上界を求めるには、すべての $j$ に対して $\Gamma\ge p_j\rho_j$ を満たす Hermite 行列 $\Gamma$ を選ぶ。各 $\Tr[(\Gamma-p_j\rho_j)E_j]\ge0$ なので、成功率は $\Tr\Gamma$ 以下である。ある測定とある $\Gamma$ が同じ値を与えれば、両者の比較だけで測定の最適性を証明できる。

:::example[赤道上の三状態の最適測定]
$\ket{\psi_j}=(\ket0+e^{2\pi ij/3}\ket1)/\sqrt2$、$j=0,1,2$ を等確率で準備する。射影を $P_j$ とすると $\sum_jP_j=3I/2$ なので $E_j=2P_j/3$ は測定を成す。どの入力でも正答率は $2/3$ である。一方、$\frac13P_j\le I/3$ だから、$\Gamma=I/3$ は最適成功率が高々 $2/3$ であることを示す。三本は二次元で線形従属なので、後述する非曖昧識別のように全てのラベルを誤りなしで答えることはできない。
:::

## 忠実度と Uhlmann の定理

状態を純粋化して内積を取るには、参照系の基底を揃える必要がある。例えば全く同じ混合状態でも、異なる参照系ラベルを使えば純粋化は直交しうる。参照系の自由度を最適化した値が状態同士の量になる。

:::definition[忠実度]{#def-qi-fidelity}
:term[忠実度]{en="fidelity" yomi="ちゅうじつど"}を $F(\rho,\sigma)=\norm{\sqrt\rho\sqrt\sigma}_1$ とする。本章では $F$ を二乗しない。二乗を忠実度と呼ぶ文献の数値と照合するときは平方根を取る。
:::

:::theorem[Uhlmann の定理]{#thm-qi-uhlmann}
状態空間と同じ次元以上の参照系を用いると、純粋化の重なり $\abs{\braket{\Psi_\rho}{\Psi_\sigma}}$ の最大値は $F(\rho,\sigma)$ である。
:::

:::proof
$\ket\Omega=\sum_j\ket{jj}$ として標準純粋化を $(\sqrt\rho\otimes I)\ket\Omega$ と $(\sqrt\sigma\otimes I)\ket\Omega$ に取る。:ref[thm-qi-purification] により自由度は参照系ユニタリである。同じ次元なら重なりは $\Tr(\sqrt\rho\sqrt\sigma U^\top)$ となる。大きい参照系では $U^\top$ の対応する圧縮が現れ、その作用素ノルムは高々 1 である。$A=\sqrt\rho\sqrt\sigma=W\Sigma V^\dagger$ を特異値分解すると、$\abs{\Tr(AC)}=\abs{\sum_j\Sigma_{jj}(V^\dagger CW)_{jj}}\le\sum_j\Sigma_{jj}$ である。$C=VW^\dagger$ はユニタリで等号を達成し、大きい参照系でも残りに恒等演算を付けられる。
:::

従って $0\le F\le1$ であり、$F=1$ なら等号を与える二つの純粋化は大域位相を除き一致するので、縮約も等しい。$\rho=\ketbra\psi\psi$ のとき $\sqrt\rho\sigma\sqrt\rho=\mel\psi\sigma\psi\,\rho$ なので $F=\sqrt{\mel\psi\sigma\psi}$ である。

:::proposition[距離と忠実度の関係]{#prop-qi-fuchs-vdg}
$1-F\le D\le\sqrt{1-F^2}$ である。
:::

:::proof
二純粋状態の差を張る二次元空間に制限すると、トレースは 0 で二乗のトレースは $2(1-\abs{\braket\psi\phi}^2)$ である。従ってその $D$ は $\sqrt{1-\abs{\braket\psi\phi}^2}$ となる。Uhlmann の最適純粋化を縮約してトレース距離の単調性を使うと上界を得る。

下界には、まず $\rho,\sigma>0$ として $M=\rho^{-1/2}(\sqrt\rho\sigma\sqrt\rho)^{1/2}\rho^{-1/2}$ と取る。$M\rho M=\sigma$ なので、$M$ の固有射影 $P_j$ に対し $q_j=m_j^2 p_j$ となる。よって $\sum_j\sqrt{p_jq_j}=\Tr(M\rho)=F$ である。古典分布では $\min(p_j,q_j)\le\sqrt{p_jq_j}$ であり、$1-\sum_j\min(p_j,q_j)$ が全変動距離である。測定後の距離は $D$ 以下なので下界が出る。特異な状態は両者に一様状態を微量混ぜて極限を取る。固有値の連続性と平方根の一様多項式近似により必要な量は連続である。
:::

:::example[対角分布と純粋状態]
$\diag(3/4,1/4)$ と $\diag(1/4,3/4)$ の忠実度は $\sqrt3/2$、距離は $1/2$ で、上界が等号となる。$\ket0,\ket+$ では忠実度と距離がともに $1/\sqrt2$ である。
:::

### 量子ビットの閉形式と距離の選び方

正の二次行列 $A=\sqrt\rho\sigma\sqrt\rho$ の固有値を $a,b$ とすれば、$(\Tr\sqrt A)^2=a+b+2\sqrt{ab}$ である。$\Tr A=\Tr(\rho\sigma)$、$\det A=(\det\rho)(\det\sigma)$ だから
$$
F(\rho,\sigma)^2=\Tr(\rho\sigma)+2\sqrt{\det\rho\,\det\sigma}
=\frac{1+\boldsymbol r\cdot\boldsymbol s+
\sqrt{(1-\norm{\boldsymbol r}^2)(1-\norm{\boldsymbol s}^2)}}2 .
$$
最後の等式では Bloch 表現と $\det\rho=(1-\norm{\boldsymbol r}^2)/4$ を用いた。二つの Bloch ベクトルの差だけでトレース距離は決まるが、忠実度にはそれぞれの長さも必要である。

:::example[同じ球内距離でも忠実度が違う]
$\rho=\diag(1,0)$ と $\sigma=I/2$ の距離は $1/2$、忠実度は $1/\sqrt2$ である。比較として $\rho'=\diag(3/4,1/4)$ と $\sigma'=\diag(1/4,3/4)$ を取ると距離は同じ $1/2$ だが、忠実度は $\sqrt3/2$ になる。純粋状態を含む最初の組で $D$ の上界は $1/\sqrt2$、後の組では $1/2$ である。距離の上界と実際の距離を比較すれば、境界のどれが鋭いかを個別に判断できる。
:::

:::example[四つの基底状態と下界の飽和]
$\rho=\diag(1/3,2/3,0,0)$、$\sigma=\diag(1/3,0,0,2/3)$ を取る。重なるラベルは第 1 成分だけなので $F=1/3$ であり、差の対角成分を足せば $D=2/3$ である。結果が第 1 成分なら識別には事前確率しか使えず、第 2 または第 4 成分なら入力が確定する。事前確率が等しければ平均成功率は $(1/3)/2+2/3=5/6$ になる。
:::

### 忠実度の単調性と純粋化による計量

:term[量子チャネル]{en="quantum channel" yomi="りょうしちゃねる"} $\Phi$ はトレースを保存し、任意の参照系を付けた $\id_R\otimes\Phi$ が正である線形写像とする。入力の基底で $\ket\Omega=\sum_i\ket{ii}$ と置くと、$J=(\id\otimes\Phi)(\ketbra\Omega\Omega)$ は正である。従って固有値の平方根をベクトルに吸収して $J=\sum_\ell\ketbra{v_\ell}{v_\ell}$ と書ける。$\ket{v_\ell}=\sum_i\ket i\otimes K_\ell\ket i$ で $K_\ell$ を定めると、参照系の $(i,j)$ ブロックの比較から $\Phi(\ketbra ij)=\sum_\ell K_\ell\ketbra ij K_\ell^\dagger$ である。トレース保存より $\sum_\ell K_\ell^\dagger K_\ell=I$ が従う。よって $V\ket x=\sum_\ell K_\ell\ket x\otimes\ket\ell$ は等長であり、$\Tr_R(V\rho V^\dagger)=\Phi(\rho)$ となる。ここで $R$ は $\ket\ell$ を基底に持つ出力補助系である。この表現の自由度は [量子チャネル](/physics/quantum-information/channels/) で扱う。

Uhlmann の定理で重なり $F(\rho,\sigma)$ を達成する純粋化を取り、両方に同じ $V$ を作用させると重なりは変わらない。出力の純粋化として許される組を一つ選んだことになるので、$F(\Phi(\rho),\Phi(\sigma))\ge F(\rho,\sigma)$ である。等長写像だけで状態を移す場合は、$V\sqrt\rho V^\dagger=\sqrt{V\rho V^\dagger}$ が固有分解から従い、積の非零特異値も保存されるので等号になる。

$d_{\mathrm B}(\rho,\sigma)=\sqrt{2-2F(\rho,\sigma)}$ を Bures 距離と呼ぶ。Uhlmann の定理からこれは、純粋化同士のベクトル距離の最小値である。内積が実非負になるよう片方の大域位相を選べば二乗距離が $2-2\abs{\langle\Psi|\Phi\rangle}$ となる。三角不等式を証明するには、中間状態 $\tau$ の純粋化を一つ固定する。同じ参照空間で $\rho,\sigma$ の純粋化をそれぞれ $\tau$ から最小距離になるよう選べる。通常のベクトルの三角不等式を使い、$\rho,\sigma$ 間の純粋化距離はその特定の組の距離以下であることを用いる。距離零なら最適な純粋化の二つが等しく、その縮約も等しい。

:::example[コピー数が忠実度に及ぼす効果]
$\ket0$ と $(3\ket0+4\ket1)/5$ は $F=3/5$ を持つ。独立に三個受け取る場合は内積が $(3/5)^3=27/125$、距離は $\sqrt{1-(27/125)^2}$ である。例えば誤り率 $P_{\mathrm{err}}$ は $\frac12(1-\sqrt{1-c^{2n}})$ だから $n=3$ で約 $0.0118$ となる。ここで使ったのは三個をまとめた最適測定である。各コピーを個別に測る場合は、その古典記録の分布から成功率を別に計算する。
:::

### 固有基底に沿った測定と横切る測定

量子ビットでは平方忠実度を行列の積の跡と行列式で計算できる。測定後の分布を対角行列で表すと、同じ公式は二結果の Bhattacharyya 係数の二乗になる。以下では測定軸を変えると保存される情報も変わることを確認する。

:::example[固有基底で測定する二つの混合]

二入力を $\boldsymbol r=(0,0,\frac{3}{5})$、$\boldsymbol s=(0,0,- \frac{3}{5})$ で指定する。測定軸は Z とする。

$$
\begin{aligned}
\rho&=\left(\begin{matrix}\frac{4}{5} & 0\\0 & \frac{1}{5}\end{matrix}\right)\\
\sigma&=\left(\begin{matrix}\frac{1}{5} & 0\\0 & \frac{4}{5}\end{matrix}\right)\\
D(\rho,\sigma)&=\|\boldsymbol r-\boldsymbol s\|/2=\frac{3}{5}\\
\Tr(\rho\sigma)&=(1+\boldsymbol r\cdot\boldsymbol s)/2=\frac{8}{25}\\
F(\rho,\sigma)^2&=\Tr(\rho\sigma)+2\sqrt{\det\rho\det\sigma}=\frac{16}{25}\\
(p_+,q_+)&=(\frac{4}{5},\frac{1}{5})\\
D_{\mathrm{obs}}&=|p_+-q_+|=\frac{3}{5}\\
F_{\mathrm{obs}}^2&=\bigl(\sqrt{p_+q_+}+\sqrt{(1-p_+)(1-q_+)}\bigr)^2=\frac{16}{25}
\end{aligned}
$$

両方がこの基底で対角なので測定で失う相対位相はない。忠実度も距離も変わらず、測定前の比較を古典計算に置き換えられる。
:::

:::example[互いに直交する軸を持つ二つの混合]

二入力を $\boldsymbol r=(\frac{4}{5},0,0)$、$\boldsymbol s=(0,0,\frac{3}{5})$ で指定する。測定軸は Z とする。

$$
\begin{aligned}
\rho&=\left(\begin{matrix}\frac{1}{2} & \frac{2}{5}\\\frac{2}{5} & \frac{1}{2}\end{matrix}\right)\\
\sigma&=\left(\begin{matrix}\frac{4}{5} & 0\\0 & \frac{1}{5}\end{matrix}\right)\\
D(\rho,\sigma)&=\|\boldsymbol r-\boldsymbol s\|/2=\frac{1}{2}\\
\Tr(\rho\sigma)&=(1+\boldsymbol r\cdot\boldsymbol s)/2=\frac{1}{2}\\
F(\rho,\sigma)^2&=\Tr(\rho\sigma)+2\sqrt{\det\rho\det\sigma}=\frac{37}{50}\\
(p_+,q_+)&=(\frac{1}{2},\frac{4}{5})\\
D_{\mathrm{obs}}&=|p_+-q_+|=\frac{3}{10}\\
F_{\mathrm{obs}}^2&=\bigl(\sqrt{p_+q_+}+\sqrt{(1-p_+)(1-q_+)}\bigr)^2=\frac{9}{10}
\end{aligned}
$$

一方の状態の $X$ 成分はこの測定で全て落ちる。それでも $Z$ の成分は残るので、二つの測定分布は同じにはならない。
:::

:::example[同じ対角成分と異なるコヒーレンス]

二入力を $\boldsymbol r=(\frac{3}{5},0,\frac{3}{5})$、$\boldsymbol s=(0,\frac{3}{5},\frac{3}{5})$ で指定する。測定軸は Z とする。

$$
\begin{aligned}
\rho&=\left(\begin{matrix}\frac{4}{5} & \frac{3}{10}\\\frac{3}{10} & \frac{1}{5}\end{matrix}\right)\\
\sigma&=\left(\begin{matrix}\frac{4}{5} & - \frac{3 i}{10}\\\frac{3 i}{10} & \frac{1}{5}\end{matrix}\right)\\
D(\rho,\sigma)&=\|\boldsymbol r-\boldsymbol s\|/2=\frac{3 \sqrt{2}}{10}\\
\Tr(\rho\sigma)&=(1+\boldsymbol r\cdot\boldsymbol s)/2=\frac{17}{25}\\
F(\rho,\sigma)^2&=\Tr(\rho\sigma)+2\sqrt{\det\rho\det\sigma}=\frac{41}{50}\\
(p_+,q_+)&=(\frac{4}{5},\frac{4}{5})\\
D_{\mathrm{obs}}&=|p_+-q_+|=0\\
F_{\mathrm{obs}}^2&=\bigl(\sqrt{p_+q_+}+\sqrt{(1-p_+)(1-q_+)}\bigr)^2=1
\end{aligned}
$$

$Z$ 基底での対角成分が同じなので、測定分布は一致する。状態は異なり、量子状態の距離は正である。後状態の忠実度 1 は入力の一致を意味しない。
:::

:::example[純粋と混合の間で測定を選ぶ]

二入力を $\boldsymbol r=(0,0,1)$、$\boldsymbol s=(\frac{3}{5},0,0)$ で指定する。測定軸は X とする。

$$
\begin{aligned}
\rho&=\left(\begin{matrix}1 & 0\\0 & 0\end{matrix}\right)\\
\sigma&=\left(\begin{matrix}\frac{1}{2} & \frac{3}{10}\\\frac{3}{10} & \frac{1}{2}\end{matrix}\right)\\
D(\rho,\sigma)&=\|\boldsymbol r-\boldsymbol s\|/2=\frac{\sqrt{34}}{10}\\
\Tr(\rho\sigma)&=(1+\boldsymbol r\cdot\boldsymbol s)/2=\frac{1}{2}\\
F(\rho,\sigma)^2&=\Tr(\rho\sigma)+2\sqrt{\det\rho\det\sigma}=\frac{1}{2}\\
(p_+,q_+)&=(\frac{1}{2},\frac{4}{5})\\
D_{\mathrm{obs}}&=|p_+-q_+|=\frac{3}{10}\\
F_{\mathrm{obs}}^2&=\bigl(\sqrt{p_+q_+}+\sqrt{(1-p_+)(1-q_+)}\bigr)^2=\frac{9}{10}
\end{aligned}
$$

純粋状態の測定分布が一様でも、もう一方の $X$ 測定の結果は一様でない。従って忠実度は 1 までは上がらない。どの軸上の情報が残るかが測定後の距離を決める。
:::

### あらゆる測定で得る古典忠実度の下限

二つの分布の古典忠実度を $\sum_a\sqrt{p_aq_a}$ とする。量子状態に POVM $E_a$ を作用させると、古典忠実度は必ず元の $F$ 以上である。等長実装を使わずに示すには、特異値分解から $F=\abs{\Tr(U\sqrt\rho\sqrt\sigma)}$ を実現するユニタリ $U$ を取る。$\sqrt\rho$ と $\sqrt\sigma$ の間に $\sum_a E_a$ を入れると
$$
F\le\sum_a\abs{\Tr(U\sqrt\rho E_a\sqrt\sigma)}
\le\sum_a \sqrt{\Tr(\rho E_a)\Tr(\sigma E_a)}.
$$
Cauchy–Schwarz で二行列 $\sqrt{E_a}\sqrt\rho U^\dagger$ と $\sqrt{E_a}\sqrt\sigma$ を使った。前者の二乗ノルムはユニタリが跡の中で消えて $\Tr(\rho E_a)$ になる。正定値の場合には、前述の $M\rho M=\sigma$ を満たす正行列 $M$ の固有基底を測ることで等号を得る。

特異な状態の場合は $\rho_\eps=(1-\eps)\rho+\eps I/d$、$\sigma_\eps=(1-\eps)\sigma+\eps I/d$ を考える。各 $\eps>0$ に対し等号を与える正規直交基底をユニタリ行列 $U_\eps$ で表す。有限次元のユニタリ群は閉かつ有界なので、$\eps\to0$ の列から収束する部分列を取れる。射影と確率と平方根は連続なので、極限の基底は元の二状態の忠実度を達成する。従って $F$ は全 POVM の古典忠実度の最小値であり、有限次元なら射影測定で実現できる。

:::example[距離を最大にする軸と忠実度を最小にする軸]
$\rho=\diag(4/5,1/5)$、$\sigma=\ketbra++$ では $F=1/\sqrt2$ である。$X$ を測れば $\sigma$ で確率が $(1,0)$、$\rho$ で $(1/2,1/2)$ であり同じ古典忠実度になる。一方、トレース距離の最適軸は差の Bloch ベクトル $(-1,0,3/5)$ に沿う。量子距離は $\sqrt{34}/10$ だが、忠実度を最小にした $X$ 測定の距離は $1/2$ である。二つの最適化を同じ測定で同時に解けるとは限らない。対照的に二状態が可換であれば共通固有基底を測ることで、確率の差も重なりも元の行列の値を反映する。
:::

### 高確率の結果への条件付け

:::lemma[穏やかな測定]{#lem-qi-gentle-measurement}
状態 $\rho$ と効果 $E$ に対し $p=\Tr(E\rho)>0$ とする。結果 $E$ の測定演算子が $\sqrt E$ であるとき、規格化した後状態 $\rho'=\sqrt E\rho\sqrt E/p$ は $D(\rho,\rho')\le\sqrt{1-p}$ を満たす。
:::

:::proof
$\rho$ の純粋化を $\ket\psi$ とすると、$\ket{\psi'}=(\sqrt E\otimes I)\ket\psi/\sqrt p$ は $\rho'$ の純粋化である。$0\le E\le I$ の固有値は $[0,1]$ にあるので $\sqrt E\ge E$ である。したがって $\braket{\psi}{\psi'}=\Tr(\rho\sqrt E)/\sqrt p\ge\sqrt p$ であり、Uhlmann の定理から $F(\rho,\rho')\ge\sqrt p$ である。:ref[prop-qi-fuchs-vdg] の上界を適用する。
:::

:::example[失敗しにくい射影と位相の擾乱]
$\ket\psi=(\sqrt{99}\ket0+\ket1)/10$ を $Z$ 基底で測り、$0$ を得た場合を考える。$p=99/100$ で後状態は $\ket0$ であるため、距離 $\sqrt{1-p}=1/10$ は上界と一致する。一方、成功演算子を $X\ketbra00$ とすれば同じ確率で $0$ を記録するが、後状態は $\ket1$ で、入力との距離は $\sqrt{99}/10$ である。条件付き状態への上界は確率に加えて平方根による更新を仮定している。
:::

## 非直交状態の識別と複製不可能定理

誤りを許す Helstrom 測定とは別に、「分からない」と答えることで、答えた結果の誤りをゼロにする設計がある。これを :term[非曖昧識別]{en="unambiguous discrimination" yomi="ひあいまいしきべつ"}という。

二純粋状態の内積の絶対値を $c\in(0,1)$ とする。相手の状態と直交する単位ベクトルに沿う射影を $Q_\psi,Q_\phi$ とし、二状態を張る平面で $E_j=Q_j/(1+c)$ と取る。二射影の和の固有値は $1\pm c$ である。これはトレース 2 と行列式 $1-c^2$ から得られる。従って $E_?=I-\sum_j E_j\ge0$ である。正答確率は $(1-c^2)/(1+c)=1-c$ となる。最適性を示すには、誤答に対応する効果が入力を消す事実を使う。例えば $\mel\phi{E_\psi}{\phi}=0$ なら $\sqrt{E_\psi}\ket\phi=0$ である。従って $\braket\psi\phi=\mel\psi{E_?}{\phi}$ であり、失敗確率 $q_\psi,q_\phi$ は $c\le\sqrt{q_\psi q_\phi}$ を満たす。等事前確率では失敗率は相加相乗平均により $c$ 以上である。

:::theorem[複製不可能定理]{#thm-qi-no-cloning}
固定した空白状態を使う等長演算は、相異なる非直交純粋状態を両方とも完全に複製できない。補助系を捨てる実装でも結論は同じである。
:::

:::proof
入力二つの内積の絶対値を $c\in(0,1)$ とする。出力の複製部分が純粋なら、全体はその純粋ベクトルと残りの積である。これは :ref[thm-qi-purification] の階数 1 の場合である。従って出力同士の内積の絶対値は $c^2\abs{\braket{e_\psi}{e_\phi}}\le c^2$ になる。等長演算は内積を保つため $c\le c^2$ となって矛盾する。
:::

:::example[直交する入力だけを複製する]
$\ket{x,y}\mapsto\ket{x,x\oplus y}$ は $\ket{0,0},\ket{1,0}$ をそれぞれ $\ket{00},\ket{11}$ に送る。だが $\ket+, \ket0$ の組を入れると $(\ket{00}+\ket{11})/\sqrt2$ となり、$\ket+\otimes\ket+$ に含まれる $\ket{01},\ket{10}$ の成分は現れない。
:::

### 混合状態の放送と可換性

二つの出力の各縮約がもとの状態 $\rho$ に一致する操作を、$\rho$ の :term[放送]{en="broadcasting" yomi="ほうそう"}という。複製では全出力を $\rho\otimes\rho$ にするが、放送では出力同士に相関が残ってよい。例えば $\rho=\sum_i p_i\ketbra ii$ に対して、基底ラベルを別系へ書けば $\sum_i p_i\ketbra{ii}{ii}$ が得られる。二つの縮約はともに $\rho$ である。

:::theorem[放送可能性と可換性]{#thm-qi-no-broadcasting}
二状態 $\rho,\sigma$ をともに放送する一つのチャネルが存在する必要十分条件は $[\rho,\sigma]=0$ である。
:::

:::proof
可換なら共通固有基底があり、その基底のラベルを二出力に記録すればよい。逆に放送できると仮定し、まず $\rho,\sigma$ は正定値とする。放送の出力を $\widetilde\rho,\widetilde\sigma$ とする。チャネルと部分トレースの忠実度単調性から
$F(\rho,\sigma)\le F(\widetilde\rho,\widetilde\sigma)\le F(\rho,\sigma)$ であり、いずれも等号になる。

$M=\rho^{-1/2}(\sqrt\rho\sigma\sqrt\rho)^{1/2}\rho^{-1/2}$ とし、その異なる固有値 $m_j$ に対応する射影を $P_j$ とする。両出力のいずれか一方で $(P_j)$ を測ると、結果を $p_j,q_j$ として $q_j=m_j^2p_j$ になる。この古典忠実度はもとの量子忠実度に等しい。$U$ を $\Tr(U\sqrt{\widetilde\rho}\sqrt{\widetilde\sigma})=F(\widetilde\rho,\widetilde\sigma)>0$ を満たすユニタリとする。前節の測定に対する忠実度の証明を、出力上の射影 $Q_j=P_j\otimes I$ に適用する。三角不等式と Cauchy–Schwarz の両方で等号だから、$Q_j\sqrt{\widetilde\sigma}$ は $Q_j\sqrt{\widetilde\rho}U^\dagger$ の実非負倍である。二つのノルムの比は $\sqrt{q_j/p_j}=m_j$ である。同じ $U$ を使い、$Q_j=I\otimes P_j$ に対しても同じ関係を得る。

$j$ について足すと $\sqrt{\widetilde\sigma}=(M\otimes I)\sqrt{\widetilde\rho}U^\dagger=(I\otimes M)\sqrt{\widetilde\rho}U^\dagger$ である。従って $\widetilde\rho$ の支持は $M\otimes I-I\otimes M$ の核に入る。核は $P_j\otimes P_j$ の像の直和であり、第 2 系を部分トレースすると異なる $j$ の間のブロックが消える。従って縮約 $\rho$ は各 $P_j$ と可換であり、$\sigma=M\rho M$ も $\rho$ と可換である。

特異な場合は $\rho'=(2\rho+\sigma)/3$、$\sigma'=(\rho+2\sigma)/3$ を取る。線形性で両状態も放送され、両者は $\rho+\sigma$ の支持 $S$ 上で正定値である。出力は両縮約の支持が $S$ なので $S\otimes S$ 上にある。実際、それぞれの直交補の射影を測る確率は零であり、正の出力行列はその直交補を消す。従って前の議論を $S$ 上で適用できる。$[\rho',\sigma']=[\rho,\sigma]/3$ だから、元の二状態も可換である。
:::

:::example[放送できる混合状態とできない混合状態]
$\rho=(I+Z/2)/2$ と $\sigma=(I-Z/2)/2$ はともに計算基底で対角である。放送写像を $\ket j\mapsto\ket{jj}$ とすれば、$\rho$ からは $(3/4)\ketbra{00}{00}+(1/4)\ketbra{11}{11}$ が得られ、計算基底の二結果は必ず一致する。一方、二つの独立な $\rho$ で結果が一致する確率は $9/16+1/16=5/8$ である。$\sigma$ を $\sigma'=(I+X/2)/2$ に取り替えると、$[\rho,\sigma']=iY/8\ne0$ であり、どちらが入力でも放送できる一つの写像は存在しない。
:::

### 事前確率が異なる非曖昧識別

二純粋状態を $\ket{\psi_0},\ket{\psi_1}$、事前を $p,1-p$、内積の絶対値を $c\in(0,1)$ とする。失敗確率 $(q_0,q_1)$ に必要な条件は $0\le q_j\le1$ と $q_0q_1\ge c^2$ である。十分性を示すため、まず入力の位相を選び込み、内積を $c$ にする。行き先として
$$
\ket{\psi_j}\longmapsto \sqrt{1-q_j}\ket{\mathrm{ok},j}
     +\sqrt{q_j}\ket{\mathrm{fail}}\ket{e_j}
$$
を取る。成功ラベルは直交させ、失敗部分の単位ベクトルは $\langle e_0|e_1\rangle=c/\sqrt{q_0q_1}$ になるよう選ぶ。これは右辺が $[0,1]$ に入るので二次元で実現できる。二入力の Gram 行列が保たれるため、線形に延長して等長写像を得る。フラグを読めば、指定した失敗確率の非曖昧識別になる。

$0<p<1$ なら最適点で積は $c^2$ である。よって $q_0\in[c^2,1]$ について $f(q_0)=pq_0+(1-p)c^2/q_0$ を最小にすればよい。$f''(q_0)=2(1-p)c^2/q_0^3>0$ なので、停留点 $c\sqrt{(1-p)/p}$ を区間へ切り詰めたものが最適な $q_0$ であり、$q_1=c^2/q_0$ である。

:::example[両状態を非曖昧に答える]
$c=1/2$、$p=2/3$ なら、最適点は $q_0=1/(2\sqrt2)$、$q_1=1/\sqrt2$ である。重み付き失敗率は $2c\sqrt{p(1-p)}=\sqrt2/3$ になる。例えば $\ket{\psi_0}=\ket0$、$\ket{\psi_1}=(\ket0+\sqrt3\ket1)/2$ では、結果 0 の効果を $\ket{\psi_1^\perp}$ への射影の $(1-q_0)/(1-c^2)$ 倍に取る。入力 0 での射影確率が $1-c^2$ だからである。結果 1 の効果も $\ket1$ への射影を $(1-q_1)/(1-c^2)$ 倍したもので、上の等長構成により二つを引いた残りも正である。
:::

:::example[まれな状態を一度も断定しない設定]
$c=1/2$、$p=9/10$ なら、停留点は $1/6$ で許容区間 $[1/4,1]$ の外に出る。したがって $q_0=1/4,q_1=1$ となり、成功率は $(9/10)(3/4)=27/40$ である。単一の射影 $\ket{\psi_1^\perp}$ を検出したときだけ 0 と答え、その他を失敗とする構成でこの値を実現できる。二択を強制すると最適成功率は $(1+\sqrt{91}/10)/2$ に変わるが、そこには誤答する分枝がある。
:::

### 確率的複製と完全なコピーの区別

成功を告げるフラグ付きの複製を考える。両入力で条件付き成功率が同じ $p$ のとき、出力の成功枝の内積の絶対値は高々 $pc^2$ で、失敗枝では高々 $1-p$ である。内積保存から $c\le pc^2+1-p$、従って $p\le1/(1+c)$ である。位相を選んで内積を正にし、成功枝の補助状態を同一に、失敗枝の全状態も同一に取る。$p=1/(1+c)$ なら出力の Gram 行列は入力と一致するので等長構成が可能である。この上界は事前確率に依存しない。

:::example[三分の二で二コピーを得る]
$c=1/2$ では上界は $2/3$ であり、成功枝の内積への寄与は $(2/3)(1/4)=1/6$、失敗枝は $1/3$ で、和は元の $1/2$ である。先に非曖昧識別を行い、答えられたときだけ二つを再準備する方法では成功率は $1/2$ になる。複製の成功枝に二つの非直交出力を残すことを認めれば、入力のラベルを確定する必要がなくなる。
:::

### 三つ以上の純粋信号を非曖昧に識別する

誤答する確率を零にするには、確定効果の支持が他候補の線形包の直交補に入る必要がある。正効果 $E$ において $\mel{v}{E}{v}=0$ なら $Ev=0$ であるためである。各候補を正の確率で答えるには、全候補の線形独立性が必要である。独立なら双対ベクトル $w_j$ を $\braket{w_j}{\psi_k}=\delta_{jk}$ で選ぶ。各状態を張る空間に効果を圧縮しても入力確率は変わらず、その中で $E_j$ の支持は $w_j$ の方向に限られる。従って全候補で同一の成功率 $t$ を要求すると、$E_j=t\ketbra{w_j}{w_j}$ を使い、$t\sum_j\ketbra{w_j}{w_j}\le I$ となる最大値を取ればよい。

:::example[二つの準位を共有する入力]

各入力を行列 $V$ の列で表す。双対ベクトルの列を $W$ と書き、$W^\dagger V=I$ から効果を得る。

$$
\begin{aligned}
V&=\left(\begin{matrix}1 & \frac{3}{5} & 0\\0 & \frac{4}{5} & 0\\0 & 0 & 1\end{matrix}\right)\\
G&=V^\dagger V=\left(\begin{matrix}1 & \frac{3}{5} & 0\\\frac{3}{5} & 1 & 0\\0 & 0 & 1\end{matrix}\right)\\
W&=(V^{-1})^\dagger=\left(\begin{matrix}1 & 0 & 0\\- \frac{3}{4} & \frac{5}{4} & 0\\0 & 0 & 1\end{matrix}\right)\\
W W^\dagger&=\left(\begin{matrix}1 & - \frac{3}{4} & 0\\- \frac{3}{4} & \frac{17}{8} & 0\\0 & 0 & 1\end{matrix}\right)\\
t&=\lambda_{\max}(WW^\dagger)^{-1}=\frac{2}{5}\\
E_{?}&=I-tWW^\dagger=\left(\begin{matrix}\frac{3}{5} & \frac{3}{10} & 0\\\frac{3}{10} & \frac{3}{20} & 0\\0 & 0 & \frac{3}{5}\end{matrix}\right)\\
P(\mathrm{ok}|j)&=t|\braket{w_j}{\psi_j}|^2=\frac{2}{5}
\end{aligned}
$$

最初の二入力は同じ平面にあり、第三入力はその直交軸にある。第三軸だけを測ると第三候補を断定できるが、残り二候補は区別できない。二つの確定効果を平面の中で相手と直交する方向に取れば、三候補に答えを出せる。
:::

:::example[重なりを大きくして退化に近付ける]

各入力を行列 $V$ の列で表す。双対ベクトルの列を $W$ と書き、$W^\dagger V=I$ から効果を得る。

$$
\begin{aligned}
V&=\left(\begin{matrix}1 & \frac{4}{5} & 0\\0 & \frac{3}{5} & 0\\0 & 0 & 1\end{matrix}\right)\\
G&=V^\dagger V=\left(\begin{matrix}1 & \frac{4}{5} & 0\\\frac{4}{5} & 1 & 0\\0 & 0 & 1\end{matrix}\right)\\
W&=(V^{-1})^\dagger=\left(\begin{matrix}1 & 0 & 0\\- \frac{4}{3} & \frac{5}{3} & 0\\0 & 0 & 1\end{matrix}\right)\\
W W^\dagger&=\left(\begin{matrix}1 & - \frac{4}{3} & 0\\- \frac{4}{3} & \frac{41}{9} & 0\\0 & 0 & 1\end{matrix}\right)\\
t&=\lambda_{\max}(WW^\dagger)^{-1}=\frac{1}{5}\\
E_{?}&=I-tWW^\dagger=\left(\begin{matrix}\frac{4}{5} & \frac{4}{15} & 0\\\frac{4}{15} & \frac{4}{45} & 0\\0 & 0 & \frac{4}{5}\end{matrix}\right)\\
P(\mathrm{ok}|j)&=t|\braket{w_j}{\psi_j}|^2=\frac{1}{5}
\end{aligned}
$$

平面での重なりが大きくなると Gram 行列の最小固有値が下がる。従って各入力に同じ成功率を要求した場合、曖昧な結果が多くなる。
:::

:::example[直交する三候補への極限]

各入力を行列 $V$ の列で表す。双対ベクトルの列を $W$ と書き、$W^\dagger V=I$ から効果を得る。

$$
\begin{aligned}
V&=\left(\begin{matrix}1 & 0 & 0\\0 & 1 & 0\\0 & 0 & 1\end{matrix}\right)\\
G&=V^\dagger V=\left(\begin{matrix}1 & 0 & 0\\0 & 1 & 0\\0 & 0 & 1\end{matrix}\right)\\
W&=(V^{-1})^\dagger=\left(\begin{matrix}1 & 0 & 0\\0 & 1 & 0\\0 & 0 & 1\end{matrix}\right)\\
W W^\dagger&=\left(\begin{matrix}1 & 0 & 0\\0 & 1 & 0\\0 & 0 & 1\end{matrix}\right)\\
t&=\lambda_{\max}(WW^\dagger)^{-1}=1\\
E_{?}&=I-tWW^\dagger=\left(\begin{matrix}0 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{matrix}\right)\\
P(\mathrm{ok}|j)&=t|\braket{w_j}{\psi_j}|^2=1
\end{aligned}
$$

双対ベクトルが元の三入力と一致する。総和は恒等であり、判定不能の効果は零になる。非曖昧識別の式が通常の射影測定を含む例である。
:::

### 反復観測で得られる誤り指数

二純粋状態の重なりの絶対値を $c\in(0,1)$ とすると、独立な $n$ コピーの等事前二択で
$$
P_e(n)=\frac{1-\sqrt{1-c^{2n}}}{2}
       =\frac{c^{2n}}{2(1+\sqrt{1-c^{2n}})}
$$
である。従って $c^{2n}/4\le P_e(n)\le c^{2n}/2$ であり、$-\log_2 P_e(n)/n$ は $-2\log_2 c$ に収束する。定数係数は $n$ で割ると消える。ここでは最適な共同測定を許しているが、同じ指数だけなら各コピーで $\ketbra{\psi_0}{\psi_0}$ の測定を行い、すべてが $\psi_0$ の結果なら 0、1 回でも外れたら 1 と答えて達成できる。入力 0 では正答し、入力 1 では $n$ 回とも通る確率が $c^{2n}$ だからである。共同測定は有限の $n$ で前の定数を改善する。

:::example[同じ指数を持つ二つの識別手順]
$c=1/2,n=4$ の共同測定では、$P_e=(1-\sqrt{255}/16)/2$ で約 $9.78\times10^{-4}$ である。個別測定で一方だけを検出する手順の誤りは $1/512\simeq1.95\times10^{-3}$ である。共同測定の下界 $1/1024$ はこの段階でも真値に近い。有限コピーで要求する確率と、コピー数を増やしたときの指数を区別すると、どこに測定の能力差があるかが分かる。
:::

## ダイヤモンドノルムとチャネルの距離

状態変換を線形写像 $\Phi$ とし、補助系 $R$ はそのまま残す。入力を対象だけに制限する場合と、対象を参照系ともつれさせる場合では区別の精度が異なる。

:::definition[ダイヤモンドノルム]{#def-qi-diamond}
:term[ダイヤモンドノルム]{en="diamond norm" yomi="だいやもんどのるむ"}を

$$
\norm\Phi_\diamond=\sup_{k\ge1}\ \sup_{\norm T_1=1}\norm{(\id_k\otimes\Phi)(T)}_1

$$

とする。ここで $T$ は参照系と入力系の上の行列である。チャネルの距離は $\norm{\Phi-\Psi}_\diamond/2$ とする。
:::

:::proposition[Hermite 保存写像での最適入力]{#prop-qi-diamond-input}
Hermite 保存写像に対し、上の最大化は、参照系次元が入力次元に等しい純粋状態に制限できる。
:::

:::proof
任意の $T$ をブロック Hermite 行列 $\begin{pmatrix}0&T\\T^\dagger&0\end{pmatrix}/2$ に埋め込む。特異値分解から、そのノルムは $\norm T_1$ である。写像を作用させた後も Hermite 保存性により同じ形で、ノルムは元の出力のノルムに等しい。よって参照系を拡大して Hermite 入力だけで最適化できる。次に Hermite 行列を固有射影の符号付き和にすると、係数の絶対値の和は 1 である。三角不等式により、ノルムが元の値以上の純粋状態が一つある。その純粋状態を入力系の基底で展開すれば、参照系の係数ベクトルは入力次元以下の個数しかない。その張る空間へ圧縮しても、等長埋め込みにより出力の非零特異値は変わらない。この有限次元の単位球面はコンパクトなので最大値を持つ。
:::

従って二チャネルが事前確率 $p,1-p$ で与えられる一回の実験の最適成功率は $(1+\norm{p\Phi-(1-p)\Psi}_\diamond)/2$ である。各純粋入力で :ref[thm-qi-helstrom] を適用し、Hermite 保存写像 $p\Phi-(1-p)\Psi$ について最大化すればよい。等確率の場合は $(1+\norm{\Phi-\Psi}_\diamond/2)/2$ になる。

一般の線形写像でも入力のトレースノルムを最適化する参照系の次元は入力次元 $d$ で足りる。任意の $T$ を特異値分解すると、$\norm T_1=1$ のとき $T$ は $\ketbra uv$ の凸結合になる。従って三角不等式により階数 1 の入力だけを調べればよい。それぞれの単位ベクトル $u,v$ は、入力系の基底で展開すると高々 $d$ 本の参照ベクトルを使う。それぞれを等長写像 $V_u,V_v$ で $d$ 次元から元の参照空間に埋め込めば、出力の左右に $V_u\otimes I$ と $V_v^\dagger\otimes I$ を掛けた形になり、非零特異値は変わらない。元の参照空間が $d$ 未満の場合は先に零方向を付け足してよい。この一般の場合、左右のベクトルを異なるものに取る必要があり、密度演算子だけで最大化できるとは限らない。

:::example[補助系で改善する識別]
量子ビットの恒等写像と $\mathcal D(T)=\Tr(T)I/2$ を考える。補助系なしでは任意の純粋状態の出力差の固有値は $\pm1/2$ で、距離は $1/2$ である。混合入力は凸性によりこれを超えない。入力を $\ket{\Phi^+}$ に拡張すると、出力差は $\ketbra{\Phi^+}{\Phi^+}-I_4/4$ である。固有値が $3/4,-1/4,-1/4,-1/4$ なので距離 $3/4$ を達成する。
:::


### Choi 行列から距離を求める

ここでは $\Phi$ を Hermite 保存写像とする。入力と同次元の参照系を取り、$\ket\Omega=\sum_j\ket{jj}$、$J=(\id\otimes\Phi)(\ketbra\Omega\Omega)$ とする。最大化する任意の単位ベクトルは $(A\otimes I)\ket\Omega$ と書け、$\Tr(A^\dagger A)=1$ である。極分解で $A=U\sqrt\rho$、$\rho=A^\dagger A$ と書けば、$U$ は出力のノルムを変えない。従って

$$

\norm\Phi_\diamond
=\max_{\rho\ge0,\ \Tr\rho=1}\norm{(\sqrt\rho\otimes I)J(\sqrt\rho\otimes I)}_1.

$$

出力行列は自己共役であり、符号演算子とのトレースでノルムが得られる。具体的には、固有分解から $\norm B_1=\max_{-I\le K\le I}\Tr(KB)$ である。$W=(\sqrt\rho\otimes I)K(\sqrt\rho\otimes I)$ と置けば、距離を求めるための半正定値計画は

$$

\max\ \Tr(JW),\qquad
\rho\ge0,\quad\Tr\rho=1,\quad
-\rho\otimes I\le W\le\rho\otimes I

$$

となる。逆に最後の二つの行列不等式を満たす $W$ は $\rho\otimes I$ の支持に収まり、その上で逆平方根を掛けると $-I\le K\le I$ となる。支持の主張は $\rho\otimes I$ の核上で $\rho\otimes I\pm W$ の二次形式が零になり、両正演算子が核を消すことから従う。従って変数変換による緩和はない。

### 量子ビットの完全混合化の厳密値

$X,Y,Z$ の共役と恒等を一様に混ぜると、どの Pauli 係数も相殺されて $\mathcal D(T)=\Tr(T)I/2$ となる。したがって $\id-\mathcal D$ は、$\id$ と各非自明 Pauli 共役の差の $1/4$ ずつの和である。どのチャネルも、純粋入力と補助系に対する出力のトレースノルムは 1 である。:ref[prop-qi-diamond-input] により、二チャネルの差のダイヤモンドノルムは高々 2 である。よって $\norm{\id-\mathcal D}_\diamond\le3/2$ であり、先の Bell 入力が等号を達成する。補助系による改善は距離で $1/2$ から $3/4$ への変化である。

### Pauli チャネル間の距離を確率から計算する

$\mathcal P_p(T)=\sum_{j=0}^3 p_j P_j T P_j^\dagger$ とし、$(P_0,P_1,P_2,P_3)=(I,X,Y,Z)$ とする。二つの確率分布 $p,q$ に対して $\norm{\mathcal P_p-\mathcal P_q}_\diamond=\sum_j|p_j-q_j|$ である。上界はチャネルのノルムが 1 であることと三角不等式から出る。下界には規格化した Bell 入力を使う。$(I\otimes P_j)\ket{\Phi^+}$ の内積は $\Tr(P_i^\dagger P_j)/2=\delta_{ij}$ なので、出力の差はこの正規直交基底で対角である。その対角成分は $p_j-q_j$ に等しく、上界を達成する。この式では距離は半分、ノルムはそのままである。

:::example[ビット反転と位相反転を見分ける]
確率 $b=1/5$ の $X$ 反転と確率 $b$ の $Z$ 反転を比較する。分布の差は $(0,b,0,-b)$ で、ダイヤモンドノルムは $2/5$、等事前の最適成功率は $3/5$ である。この組では補助系を使わなくても、$\ket0$ を入力して $Z$ を測れば距離 $b$ を達成する。$X$ 反転では $1$ が確率 $b$、$Z$ 反転では $1$ が確率 0 だからである。Bell 測定を使った最適化は、常に補助系が必須であるという意味ではない。
:::

### ユニタリの相対スペクトルで区別する

チャネルが $U,V$ で共役する写像なら、純粋入力の後の二出力の重なりは $\Tr(\rho U^\dagger V)$ である。ここで $\rho$ は入力系の縮約である。$U^\dagger V$ を対角化すると、このトレースはその複素固有値の凸結合になる。任意の係数を達成するには、入力として対応する固有ベクトルに、その係数の正平方根を掛けて足せばよい。したがって補助系は不要で、固有値の凸包と原点の最短距離 $\nu$ を求めればダイヤモンドノルムは $2\sqrt{1-\nu^2}$ となる。

:::example[位相ゲートを一回で判定する]
$U=I$、$V=\diag(1,i)$ では固有値を結ぶ線分の中点 $(1+i)/2$ が原点から最も近く、$\nu=1/\sqrt2$、ノルムは $\sqrt2$ である。$\ket+$ を入力すれば、二出力 $\ket+$ と $(\ket0+i\ket1)/\sqrt2$ がこの重なりを持つ。一方、$V=\diag(1,-1)$ では固有値の中点が原点なので、同じ $\ket+$ が完全識別を与える。なお $V=e^{i\theta}I$ なら固有値は一点で $\nu=1$、ダイヤモンドノルムは 0 である。
:::

### 合成された実験の誤差を足し合わせる

写像を二回合成したとき、$\norm{\Phi\circ\Lambda}_\diamond\le\norm\Phi_\diamond\norm\Lambda_\diamond$ である。任意の参照系と任意のトレースノルム 1 の入力に対して、内側の写像を適用した結果をそのノルムで割り、外側の写像の定義を適用すればよい。内側の出力が零の場合の不等式は自明である。これとチャネルのノルム 1 を用いると、$m$ 段の二系列の差は各段の差の和で抑えられる。途中の段を一つずつ交換し、両側に残る写像との合成に劣乗法性を使う。

:::example[繰り返す位相反転]
確率 $1/10$ で $Z$ が起こるチャネルを独立に二回使う。二反転は打ち消し合うので、奇数回起こる確率は $2(1/10)(9/10)=9/50$ である。従って恒等写像からのダイヤモンド距離は $9/50$ で、段ごとの和による上界 $1/5$ 以下である。一般の $m$ 回では奇偶の期待値が $(1-2/10)^m=(4/5)^m$ になり、距離は $(1-(4/5)^m)/2$ と計算できる。
:::

### Bell プローブの出力を対角化する

Pauli 分布を $I,X,Y,Z$ の順に書く。差の作用を $RB$ の密度行列上で計算してからノルムを取れば、ダイヤモンドノルムの半分と二択成功率を混同しにくい。Bell 基底の固有値と計算基底での差行列の両方を確認する。

:::example[二種類の対称な Pauli 混合]

二つのチャネルが等しい事前確率で与えられる。Bell 入力を通した差を $\Delta_{RB}$ とする。

$$
\begin{aligned}
p&=(\frac{5}{8},\frac{1}{8},\frac{1}{8},\frac{1}{8})\\
q&=(\frac{1}{4},\frac{1}{4},\frac{1}{4},\frac{1}{4})\\
\Delta_{RB}&=(\id\otimes(\mathcal P_p-\mathcal P_q))(\ketbra{\Phi^+}{\Phi^+})=\left(\begin{matrix}\frac{1}{8} & 0 & 0 & \frac{1}{4}\\0 & - \frac{1}{8} & 0 & 0\\0 & 0 & - \frac{1}{8} & 0\\\frac{1}{4} & 0 & 0 & \frac{1}{8}\end{matrix}\right)\\
\lambda(\Delta_{RB})&=(p_j-q_j)_{j=0}^3=(\frac{3}{8},- \frac{1}{8},- \frac{1}{8},- \frac{1}{8})\\
\|\mathcal P_p-\mathcal P_q\|_\diamond&=\sum_j|\lambda_j(\Delta_{RB})|=\frac{3}{4}\\
P_{\mathrm{succ}}&=\tfrac12+\|\mathcal P_p-\mathcal P_q\|_\diamond/4=\frac{11}{16}
\end{aligned}
$$

第 1 チャネルの Bloch 倍率は $1/2$、第 2 は 0 なので、補助なしの純粋状態での距離は $1/4$ である。Bell プローブでは三つの誤りラベルの違いも個別に読める。
:::

:::example[異なる回転成分を同じ確率で含む混合]

二つのチャネルが等しい事前確率で与えられる。Bell 入力を通した差を $\Delta_{RB}$ とする。

$$
\begin{aligned}
p&=(\frac{3}{5},\frac{1}{5},0,\frac{1}{5})\\
q&=(\frac{3}{5},0,\frac{1}{5},\frac{1}{5})\\
\Delta_{RB}&=(\id\otimes(\mathcal P_p-\mathcal P_q))(\ketbra{\Phi^+}{\Phi^+})=\left(\begin{matrix}0 & 0 & 0 & 0\\0 & 0 & \frac{1}{5} & 0\\0 & \frac{1}{5} & 0 & 0\\0 & 0 & 0 & 0\end{matrix}\right)\\
\lambda(\Delta_{RB})&=(p_j-q_j)_{j=0}^3=(0,\frac{1}{5},- \frac{1}{5},0)\\
\|\mathcal P_p-\mathcal P_q\|_\diamond&=\sum_j|\lambda_j(\Delta_{RB})|=\frac{2}{5}\\
P_{\mathrm{succ}}&=\tfrac12+\|\mathcal P_p-\mathcal P_q\|_\diamond/4=\frac{3}{5}
\end{aligned}
$$

共通のラベルを読んだ分枝では事前に頼るしかない。$\ket+$ を使うと $X$ は不変で $Y$ は直交状態へ移す。この例では補助なしでも同じ距離に達する。
:::

:::example[二つのチャネルに共通する誤りを引き去る]

二つのチャネルが等しい事前確率で与えられる。Bell 入力を通した差を $\Delta_{RB}$ とする。

$$
\begin{aligned}
p&=(\frac{1}{2},\frac{1}{4},\frac{1}{4},0)\\
q&=(\frac{1}{8},\frac{1}{4},\frac{1}{4},\frac{3}{8})\\
\Delta_{RB}&=(\id\otimes(\mathcal P_p-\mathcal P_q))(\ketbra{\Phi^+}{\Phi^+})=\left(\begin{matrix}0 & 0 & 0 & \frac{3}{8}\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\\frac{3}{8} & 0 & 0 & 0\end{matrix}\right)\\
\lambda(\Delta_{RB})&=(p_j-q_j)_{j=0}^3=(\frac{3}{8},0,0,- \frac{3}{8})\\
\|\mathcal P_p-\mathcal P_q\|_\diamond&=\sum_j|\lambda_j(\Delta_{RB})|=\frac{3}{4}\\
P_{\mathrm{succ}}&=\tfrac12+\|\mathcal P_p-\mathcal P_q\|_\diamond/4=\frac{11}{16}
\end{aligned}
$$

最適測定では二つの共通ラベルの割当てを変更しても成功率が同じである。差の固有値零の空間にどのような効果を置いても Helstrom 値が変わらないためである。
:::

## 量子状態トモグラフィー

:term[量子状態トモグラフィー]{en="quantum state tomography" yomi="りょうしじょうたいともぐらふぃー"}では、同じ状態を独立に準備できると仮定する。一つの標本で複数の非可換な観測量を同時に読む代わりに、標本を測定設定ごとに分ける。

1 量子ビットで $X,Y,Z$ の測定を各 $m$ 回行い、$\pm1$ の標本平均を $\hat r_j$ とする。期待値は $r_j$、分散は独立性より $(1-r_j^2)/m$ である。従って $\mathbb E\norm{\hat{\boldsymbol r}-\boldsymbol r}^2\le3/m$ となり、Markov の不等式で長さ $\eps$ 以上の誤差確率は $3/(m\eps^2)$ 以下となる。

:::example[有限標本で球をはみ出す場合]
$X,Y,Z$ の平均が $(1,1,0)$ となると線形逆変換は固有値 $(1\pm\sqrt2)/2$ を返す。$\hat{\boldsymbol r}$ を長さ 1 に縮めれば正の推定値が得られる。一般に推定ベクトル $t\boldsymbol u$ が $t>1$ のとき、任意の真の $\boldsymbol r$ に対し距離の二乗の減少分は $(t-1)(t+1-2\boldsymbol u\cdot\boldsymbol r)\ge0$ である。
:::

:::example[多数の量子ビット]
$P\in\{I,X,Y,Z\}^{\otimes n}$ は積トレースの公式から $\Tr(PQ)=2^n\delta_{P,Q}$ を満たす。$4^n$ 本の直交基底なので $\rho=2^{-n}\sum_P\Tr(P\rho)P$ である。例えば Bell 状態の非零係数は $II,XX,YY,ZZ$ に対し $1,1,-1,1$ である。局所平均だけを測っても $XX$ などの積の係数が欠ける。
:::


### 一つの測定設定による再構成

:::example[正四面体 POVM]
$\boldsymbol v_j$ を $(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)$ のそれぞれを $\sqrt3$ で割ったものとし、$E_j=(I+\boldsymbol v_j\cdot\boldsymbol\sigma)/4$ とする。各効果は固有値 $0,1/2$ を持つ。$\sum_j\boldsymbol v_j=0$ なので和は $I$ である。測定確率は $p_j=(1+\boldsymbol v_j\cdot\boldsymbol r)/4$ となる。符号を足し上げると $\sum_j\boldsymbol v_j\boldsymbol v_j^\top=(4/3)I_3$ なので

$$

\boldsymbol r=3\sum_j p_j\boldsymbol v_j
$$

である。例えば真の状態が $\boldsymbol r=\boldsymbol v_1$ なら確率は $1/2,1/6,1/6,1/6$ であり、式に代入すれば元のベクトルに戻る。この測定は単一の設定でも全三座標を識別できる。各結果の頻度を確率に代入したときの推定値は有限標本で球外に出うるので、必要なら球へ射影する。
:::

### 測定行列の階数と線形逆変換

有限 POVM が情報完全であるとは、測定確率が密度行列を一意に定めることである。同値な条件は、効果の実線形包が Hermite 行列の全空間であることである。包が全空間なら、確率が等しい差はトレース内積で全空間に直交する。逆に張らないなら、全効果に直交する非零の Hermite 行列 $D$ を取る。$\sum_a E_a=I$ より $\Tr D=0$ である。$\eps\in(0,1/(d\norm D_\infty))$ なら $I/d\pm\eps D$ が異なる状態であり、測定確率は等しい。したがって情報完全な POVM の結果数は少なくとも $d^2$ である。

:::example[一つの軸が不足した測定]
$X,Z$ の射影測定を半々で選び、四つのラベルを記録する。効果は四つあっても、張るのは $I,X,Z$ の三次元である。$\rho=(I+Y/2)/2$、$\sigma=(I-Y/2)/2$ なら四ラベルはすべて $1/4$ で、測定統計から区別できない。一方、$Y$ を測れば $+1$ 確率は $3/4,1/4$ と異なる。結果数と測定行列の階数を両方調べる必要がある。
:::

### 最尤推定と誤差の位置依存性

効果 $E_a$ に対する観測度数が $n_a$ なら、独立試行の対数尤度は定数を除いて $\ell(\rho)=\sum_a n_a\log\Tr(E_a\rho)$ である。対数は正の範囲で凹であり、確率は状態のアフィン関数なので $\ell$ は状態集合上で凹である。度数が正の結果に確率零を割り当てる状態では $\ell=-\infty$ とする。正の度数を持つ効果が全状態を区別し、有限の尤度を持つなら二点を結ぶ線分上で少なくとも一項が狭義に凹なので最大解は一意である。

:::example[最尤推定と単純な球への射影]
$X,Y$ を各 20 回測って全て $+1$、$Z$ を 20 回測って $\pm1$ を 10 回ずつ得たとする。$r_z$ の尤度因子 $(1-r_z^2)^{10}$ は $r_z=0$ で最大になり、球の制約も緩む。したがって最尤推定では $r_z=0$ である。対称性と凹性により $r_x=r_y$ と取れ、$(1+r_x)^{20}(1+r_y)^{20}$ はその値が大きいほど増える。よって最尤解は $(1/\sqrt2,1/\sqrt2,0)$ である。

今度は $X$ の測定数だけ 40 にする。他の記録は同じで、$X$ はなお全て $+1$ とする。線形逆変換は依然 $(1,1,0)$ だが、最尤解は同じ点にない。境界を $(\cos t,\sin t,0)$ と書けば、$t=\pi/4$ における対数尤度の導関数は $-20/(1+\sqrt2)<0$ である。角 $t$ を少し減らすと尤度が上がるので、球への最短射影とは異なる。
:::

:::example[正四面体測定での二乗誤差]
正四面体の $n$ 回の独立な結果を $J_1,\dots,J_n$ とする。$\hat{\boldsymbol r}=3n^{-1}\sum_i\boldsymbol v_{J_i}$ は平均が $\boldsymbol r$ に等しく、一回の標本 $T_i=3\boldsymbol v_{J_i}$ は常に長さ 3 である。独立性から $\mathbb E\norm{\hat{\boldsymbol r}-\boldsymbol r}^2=(9-\norm{\boldsymbol r}^2)/n$ である。$n=400$、真の状態が球面上なら値は $1/50$ となる。比較として Pauli 三軸を各 $m$ 回測る場合は $(3-\norm{\boldsymbol r}^2)/m$ であり、総標本数 $3m$ に直すと同じ球面上で $6/(3m)$ となる。この比較は二乗誤差についてであり、設定を切り替えるための実験時間を含める場合の優劣は別に調べる。
:::

### 観測度数から正当な行列を組み立てる

三つの別々の設定で正結果の度数を得たとき、各 $\hat r_i=2n_{i,+}/m_i-1$ をまとめる。球への射影推定では長さが 1 を超えた場合にだけ長さ 1 まで縮める。次の計算で、行列の固有値と純度を併せて確認する。

:::example[零の平均と有限分散]

観測の $X,Y,Z$ の順で、正結果の数と試行数の組を $(10,20)$、$(10,20)$、$(10,20)$ とする。

$$
\begin{aligned}
\hat{\boldsymbol r}&=(0,0,0)\\
\|\hat{\boldsymbol r}\|^2&=\hat r_x^2+\hat r_y^2+\hat r_z^2=0\\
\hat\rho_{\mathrm{lin}}&=(I+\hat{\boldsymbol r}\cdot\boldsymbol\sigma)/2=\left(\begin{matrix}\frac{1}{2} & 0\\0 & \frac{1}{2}\end{matrix}\right)\\
\lambda_{\mathrm{lin}}&=(\frac{1}{2},\frac{1}{2})\\
\hat\rho_{\mathrm{proj}}&=\left(\begin{matrix}\frac{1}{2} & 0\\0 & \frac{1}{2}\end{matrix}\right)\\
\Tr(\hat\rho_{\mathrm{proj}}^2)&=(1+\|\boldsymbol r_{\mathrm{proj}}\|^2)/2=\frac{1}{2}
\end{aligned}
$$

線形逆変換と射影推定は同じになる。三標本平均が零でも、標本の数が有限なので真の状態が厳密に一様であるとは結論できない。
:::

:::example[非物理的な高偏りを球面へ戻す]

観測の $X,Y,Z$ の順で、正結果の数と試行数の組を $(10,10)$、$(9,10)$、$(5,10)$ とする。

$$
\begin{aligned}
\hat{\boldsymbol r}&=(1,\frac{4}{5},0)\\
\|\hat{\boldsymbol r}\|^2&=\hat r_x^2+\hat r_y^2+\hat r_z^2=\frac{41}{25}\\
\hat\rho_{\mathrm{lin}}&=(I+\hat{\boldsymbol r}\cdot\boldsymbol\sigma)/2=\left(\begin{matrix}\frac{1}{2} & \frac{1}{2} - \frac{2 i}{5}\\\frac{1}{2} + \frac{2 i}{5} & \frac{1}{2}\end{matrix}\right)\\
\lambda_{\mathrm{lin}}&=(\frac{1}{2} + \frac{\sqrt{41}}{10},\frac{1}{2} - \frac{\sqrt{41}}{10})\\
\hat\rho_{\mathrm{proj}}&=\left(\begin{matrix}\frac{1}{2} & \frac{\sqrt{41} \left(5 - 4 i\right)}{82}\\\frac{\sqrt{41} \left(5 + 4 i\right)}{82} & \frac{1}{2}\end{matrix}\right)\\
\Tr(\hat\rho_{\mathrm{proj}}^2)&=(1+\|\boldsymbol r_{\mathrm{proj}}\|^2)/2=1
\end{aligned}
$$

長さへの割り算は二つの非零成分の比を変えない。非物理的な推定値の純度は 1 を越える一方で、射影後は球面上の純粋状態に対応する。
:::

:::example[非対角成分の虚部を読む]

観測の $X,Y,Z$ の順で、正結果の数と試行数の組を $(39,60)$、$(24,60)$、$(42,60)$ とする。

$$
\begin{aligned}
\hat{\boldsymbol r}&=(\frac{3}{10},- \frac{1}{5},\frac{2}{5})\\
\|\hat{\boldsymbol r}\|^2&=\hat r_x^2+\hat r_y^2+\hat r_z^2=\frac{29}{100}\\
\hat\rho_{\mathrm{lin}}&=(I+\hat{\boldsymbol r}\cdot\boldsymbol\sigma)/2=\left(\begin{matrix}\frac{7}{10} & \frac{3}{20} + \frac{i}{10}\\\frac{3}{20} - \frac{i}{10} & \frac{3}{10}\end{matrix}\right)\\
\lambda_{\mathrm{lin}}&=(\frac{\sqrt{29}}{20} + \frac{1}{2},\frac{1}{2} - \frac{\sqrt{29}}{20})\\
\hat\rho_{\mathrm{proj}}&=\left(\begin{matrix}\frac{7}{10} & \frac{3}{20} + \frac{i}{10}\\\frac{3}{20} - \frac{i}{10} & \frac{3}{10}\end{matrix}\right)\\
\Tr(\hat\rho_{\mathrm{proj}}^2)&=(1+\|\boldsymbol r_{\mathrm{proj}}\|^2)/2=\frac{129}{200}
\end{aligned}
$$

負の $Y$ 平均は、密度行列の上の非対角成分の虚部を正にする。ここで虚部の符号を読み違えると $Y$ の確率だけが入れ替わる。
:::

:::example[設定ごとに試行数が異なる場合]

観測の $X,Y,Z$ の順で、正結果の数と試行数の組を $(8,15)$、$(19,40)$、$(40,50)$ とする。

$$
\begin{aligned}
\hat{\boldsymbol r}&=(\frac{1}{15},- \frac{1}{20},\frac{3}{5})\\
\|\hat{\boldsymbol r}\|^2&=\hat r_x^2+\hat r_y^2+\hat r_z^2=\frac{1321}{3600}\\
\hat\rho_{\mathrm{lin}}&=(I+\hat{\boldsymbol r}\cdot\boldsymbol\sigma)/2=\left(\begin{matrix}\frac{4}{5} & \frac{1}{30} + \frac{i}{40}\\\frac{1}{30} - \frac{i}{40} & \frac{1}{5}\end{matrix}\right)\\
\lambda_{\mathrm{lin}}&=(\frac{\sqrt{1321}}{120} + \frac{1}{2},\frac{1}{2} - \frac{\sqrt{1321}}{120})\\
\hat\rho_{\mathrm{proj}}&=\left(\begin{matrix}\frac{4}{5} & \frac{1}{30} + \frac{i}{40}\\\frac{1}{30} - \frac{i}{40} & \frac{1}{5}\end{matrix}\right)\\
\Tr(\hat\rho_{\mathrm{proj}}^2)&=(1+\|\boldsymbol r_{\mathrm{proj}}\|^2)/2=\frac{4921}{7200}
\end{aligned}
$$

この場合は三方向の標本数が違う。各座標はその軸での試行数で割る必要がある。尤度を最大にする推定では、異なる試行数は各座標の重みも変える。
:::

### 古典的な標本誤差の同時評価

三軸の二値測定を独立に $m$ 回ずつ行う。値が $\pm1$ の標本 $S$ と平均 $\mu$ に対し $\log\mathbb E e^{t(S-\mu)}$ を $f(t)$ とすると、$f(0)=f'(0)=0$ である。二階微分は、もとの重みに $e^{tS}$ を掛けて規格化した分布での $S$ の分散であり、高々 1 である。従って $f(t)\le t^2/2$ である。標本平均のずれが $a>0$ 以上となる確率に Markov の不等式を適用し、$(t^2/2-ta)m$ を $t=a$ で最小化すると $e^{-ma^2/2}$ である。負方向も同じなので両側はその 2 倍になる。

三軸の長さの誤差が $\eps$ 以上なら、少なくとも一座標で $\eps/\sqrt3$ 以上のずれがある。合併上界より失敗率は $6e^{-m\eps^2/6}$ 以下である。Bloch ベクトルの長さの半分がトレース距離なので、球へ射影した推定について $D>\eps/2$ となる確率もこの値以下である。例えば目標距離 $0.05$、失敗率 $0.01$ なら $\eps=0.1$ を取り、$m\ge600\ln600$、すなわち各軸 3839 回でこの評価を満たす。これは状態に依らない十分条件であり、実際の推定の分散に応じた必要回数と一致するとは限らない。

## 演習問題

::::exercise{level=basic}
$\ket0,\ket1$ を事前確率 $1/3,2/3$ で受け取るときの最適成功率を求めよ。

:::solution
$\Delta=\diag(1/3,-2/3)$ のトレースノルムは 1 である。Helstrom の公式で成功率は 1 であり、計算基底測定が達成する。
:::
::::

::::exercise{level=basic}
$D(\diag(1,0),I/2)$ を求めよ。

:::solution
差の固有値は $\pm1/2$ なので絶対値の和は 1、従って距離は $1/2$ である。
:::
::::

::::exercise{level=basic}
$F(\ketbra00,I/2)$ を求めよ。

:::solution
第 1 状態が純粋なので $F=\sqrt{\mel0{I/2}0}=1/\sqrt2$ である。
:::
::::

::::exercise{level=basic}
$X$ 測定 10 回のうち $+1$ が 7 回、$-1$ が 3 回なら $r_x$ の推定値はいくらか。

:::solution
結果の合計を回数で割るので $(7-3)/10=2/5$ である。
:::
::::

::::exercise{level=advanced}
純粋状態の完全識別が可能であることと、直交性が同値であることを示せ。

:::solution
直交なら対応する射影で識別できる。逆に完全識別なら $E$ が一方の状態に確率 1、他方に確率 0 を与える。正性より $E$ は一方を固定し、他方を消す。Hermite 性から二状態の内積は零である。
:::
::::

::::exercise{level=advanced}
内積が $1/2$ の二純粋状態の等確率非曖昧識別の最適成功率を求めよ。

:::solution
直交補への効果を $1/(1+1/2)=2/3$ 倍すれば、$(2/3)(1-1/4)=1/2$ を達成する。失敗確率の相加相乗平均による下限も $1/2$ なので最適である。
:::
::::

::::exercise{level=advanced}
任意の二状態のテンソル積について $F(\rho_1\otimes\rho_2,\sigma_1\otimes\sigma_2)=F(\rho_1,\sigma_1)F(\rho_2,\sigma_2)$ を示せ。

:::solution
平方根は固有分解からテンソル積を保つ。さらに二行列の特異値分解をテンソル積にすると、積の全特異値は各因子の特異値の全組の積になる。その和も積になる。
:::
::::

::::exercise{level=advanced}
どのチャネルでもダイヤモンドノルムが 1 であることを示せ。

:::solution
チャネルは任意の補助系と合わせても正性とトレースを保つ。:ref[prop-qi-diamond-input] によって純粋入力だけを調べれば、出力は密度行列でトレースノルム 1 である。
:::
::::

::::exercise{level=advanced}
2 種類の量子ビット Pauli チャネルが等確率で与えられる。$I,X,Y,Z$ の順で確率が $(1/2,1/2,0,0)$ と $(1/2,0,1/2,0)$ であるとき、最適識別成功率を求めよ。参照を使わない入力でもその値を達成できるか。

:::solution
Pauli チャネルの差は $\frac12(\mathcal X-\mathcal Y)$ なので、各チャネルのダイヤモンドノルムが 1 であることと三角不等式からノルムは高々 1 である。Bell 入力で $X$ と $Y$ のラベルは直交し、確率差が $\pm1/2$ なのでこのノルムを達成する。従って最適成功率は $(1+1/2)/2=3/4$ である。

補助を使わず $\ket+$ を入れれば、$I$ と $X$ は $\ket+$ を固定する。一方 $Y$ は位相を除いて $\ket-$ へ移す。従って前者の出力は $\ketbra++$、後者は $I/2$ である。測定 $X$ で $-1$ のときだけ後者と答えると、第一仮説では常に正答、第二仮説では半分正答する。両者の平均は再び $3/4$ であり、参照なしでも最適値になる。
:::
::::

## 確認テスト

::::quiz
二状態の事前確率が等しいときの最適な二択成功率はどれか。

- [ ] $D$
- [x] $(1+D)/2$
- [ ] $(1-D)/2$

:::explanation
:ref[thm-qi-helstrom] に $p=1/2$ を代入すると $\norm{\Delta}_1=D$ となる。残る $(1-D)/2$ は誤り率である。
:::
::::

::::quiz{answer="0.75" tol="1e-6"}
Bell 入力で恒等写像と完全混合への置換写像を区別する例において、出力間距離を求めよ。

:::explanation
差の固有値は $3/4$ が一つ、$-1/4$ が三つである。従って距離は $(3/4+3/4)/2=3/4$ である。
:::
::::

::::quiz
非曖昧識別について正しいものを選べ。

- [x] 判定不能の結果を許す。
- [ ] 判定した結果の誤りを最小にするため、必ず二択で答える。

:::explanation
非曖昧識別では答えを出したときの誤りを零にする。二つの純粋信号で内積の絶対値が $c>0$ なら、判定不能の効果に対する Cauchy–Schwarz の不等式から $q_0q_1\ge c^2$ となる。従って判定不能は正の確率で発生する。
:::
::::

::::quiz{answer="3" tol="1e-6"}
量子ビットの任意の状態を表す Bloch ベクトルには実座標がいくつあるか。

:::explanation
トレース 1 の Hermite 行列は $I/2$ に Pauli 行列の実線形結合を加えた形であり、$X,Y,Z$ の 3 係数を持つ。
:::
::::

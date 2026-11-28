---
title: 漸近理論
description: '収束の計算を整理し、推定量の正規近似、局所代替仮説、漸近有効性を導く。'
notation:
  - symbol: '\to_p'
    meaning: 確率収束
  - symbol: '\Rightarrow'
    meaning: 分布収束
  - symbol: 'o_P(1),O_P(1)'
    meaning: 確率収束で無視できる項と確率的に有界な項
  - symbol: '\Delta_n'
    meaning: 局所尤度展開の線形項
---

## この章の目標

- 確率収束と分布収束を区別し、変形後の極限を求める。
- 最尤推定量と M・Z 推定量の一致性と正規近似を組み立てる。
- 超有効性が局所代替で崩れる仕組みを説明する。
- 局所尤度比を使って別の分布の下での極限を求める。

## 確率収束と分布収束の計算と Slutsky の補題

推定量が真値に近づくことと、拡大した誤差の分布が一定の形に近づくことは異なる。標準誤差の推定値を代入するとき、前者と後者を一緒に扱う。

:::definition[収束の記号]{#def-stat-convergence}
全 $\epsilon>0$ で $P(\norm{X_n-X}>\epsilon)\to0$ なら $X_n$ は $X$ に :term[確率収束]{en="convergence in probability" yomi="かくりつしゅうそく"} するといい、$X_n\to_p X$ と書く。有界連続関数 $f$ のすべてについて $E f(X_n)\to E f(X)$ なら :term[分布収束]{en="convergence in distribution" yomi="ぶんぷしゅうそく"} といい、$X_n\Rightarrow X$ と書く。$X_n=o_P(1)$ は $X_n\to_p0$、$X_n=O_P(1)$ は任意の $\epsilon>0$ に対し十分大きい固定 $M$ と十分大きい全 $n$ で $P(\norm{X_n}>M)<\epsilon$ を意味する。
:::

:::theorem[Slutsky の補題]{#thm-stat-slutsky}
$X_n\Rightarrow X$、$Y_n\to_p c$（$c$ は定数）なら $(X_n,Y_n)\Rightarrow(X,c)$ である。従って連続関数を適用した和、積、分母が非零の場合の商も対応する極限を持つ。
:::

:::proof
分布収束から $X_n$ は確率的に有界である。実際、大きい球の外で $1$、小さい球の中で $0$ の連続な関数で裾確率を抑えられる。任意の有界連続 $f$ について $X_n$ が球内で $Y_n$ が $c$ に近いとき、一様連続性により $f(X_n,Y_n)-f(X_n,c)$ は一様に小さい。補集合の期待値は関数の有界性と確率で抑えられる。従って差の期待値はゼロに近づき、$f(X_n,c)$ の期待値は分布収束で扱える。商は $c$ の周りで分母をゼロから離した連続関数に置き換えられ、その置換は確率 $1$ に近づく範囲で元と一致する。
:::

:::example[標準誤差を推定する]
有限な正分散 $\sigma^2$ の独立標本では、:ref[thm-prob-iid-clt] により $\sqrt n(\bar X-\mu)\Rightarrow N(0,\sigma^2)$ である。大数の法則を $X_i,X_i^2$ に適用すると $S^2\to_p\sigma^2$ なので、$S$ で割った誤差は $N(0,1)$ に近づく。$\sigma^2=4$ なら分子の漸近標準偏差は $2$ である。
:::

:::example[確率的な二項の結合]
毎回同じ $Z\sim N(0,1)$ に対し $X_n=Z$、$Y_n=(-1)^n Z$ とする。各周辺分布は同じ標準正規だが、和は偶数で $2Z$、奇数で $0$ である。二つの周辺分布収束だけでは和の法則を得られない。
:::

### 分布関数と確率収束の関係

分布収束の定義では、すべての有界連続関数で平均を比較する。標本統計量が実数値なら、次の判定によって累積確率の計算だけで済む。

:::proposition[分布関数による判定]{#prop-stat-cdf-convergence}
実数値確率変数 $Y_n,Y$ の分布関数を $F_n,F$ とする。$Y_n\Rightarrow Y$ と、$F$ のすべての連続点 $t$ で $F_n(t)\to F(t)$ であることは同値である。特に $Y_n\Rightarrow c$ が定数への収束なら $Y_n\to_p c$ である。
:::

:::proof
$t$ を連続点とする。$\mathbf1_{(-\infty,t]}$ を、区間 $[t-\delta,t]$ で $1$ から $0$ になる折線で下から、$[t,t+\delta]$ で $1$ から $0$ になる折線で上から挟む。分布収束で折線の期待値の極限を取り、$\delta\downarrow0$ とすると上下とも $F(t)$ になる。

逆に、有界連続関数 $f$ を取る。両端を $F$ の連続点に選び、区間外の $Y$ の確率を任意に小さくできる。端点の累積確率の収束から、$Y_n$ の区間外の確率も十分大きい $n$ で小さい。区間内では $f$ は一様連続なので、分割点を $F$ の連続点に取った細かい分割上の階段関数で近似できる。各区間の確率は両端の累積確率の差で収束する。近似幅と区間外の確率をゼロに近づければ期待値の収束を得る。定数の場合には $c-\epsilon,c+\epsilon$ が連続点なので、この区間の外へ出る確率はゼロに近づく。
:::

$Y_n\to_p Y$ なら $Y_n\Rightarrow Y$ である。実際、有界連続 $f$ は $Y$ が大きい球内にあり $\norm{Y_n-Y}\le 1$ のとき、一段大きい球で一様連続である。この球で差を抑え、球外と $\norm{Y_n-Y}>\epsilon$ の寄与を関数の上界で処理する。逆に、定数以外への分布収束は同じ確率空間上の差が消えることまで要求しない。例えば $Y\sim N(0,1)$ に対して $Y_n=-Y$ とおくと $F_n=F$ だが、$P(|Y_n-Y|>2)=P(|Y|>1)>0$ である。

### 確率的な剰余とモーメント

正の数列 $a_n$ に対し $Y_n=O_P(a_n)$ を $Y_n/a_n=O_P(1)$、$Y_n=o_P(a_n)$ を $Y_n/a_n=o_P(1)$ の略記とする。$X_n=O_P(a_n)$ と $Y_n=O_P(b_n)$ なら
$$
X_n+Y_n=O_P(a_n+b_n),\qquad X_nY_n=O_P(a_nb_n)
$$
である。実際、各正規化量を同じ $M$ で抑える事象の補集合を和事象で評価すれば、左辺の正規化量はそれぞれ $M,M^2$ 以下である。一方の因子が $o_P$ なら、もう一方をまず固定の $M$ で抑え、残りを確率収束させると積は $o_P(a_nb_n)$ になる。$Y_n\to_p c\ne0$ の逆数は、$|Y_n-c|\le |c|/2$ で
$$
\left|\frac1{Y_n}-\frac1c\right|\le \frac{2|Y_n-c|}{c^2}
$$
と抑えられる。逆数の定義に使えない $Y_n=0$ の事象では任意の有限値を与えても極限は同じである。

:::example[標本分散を二つの剰余に分ける]
$E X^4<\infty$、$E X=\mu$ とし、分母を $n$ とした標本分散を $V_n=n^{-1}\sum_i(X_i-\bar X)^2$ とする。恒等式
$$
V_n-\sigma^2=\frac1n\sum_i\{(X_i-\mu)^2-\sigma^2\}-(\bar X-\mu)^2
$$
の右端は $O_P(n^{-1})$ である。最初の項に中心極限定理を使うと、$\sqrt n(V_n-\sigma^2)$ の極限分散は $E(X-\mu)^4-\sigma^4$ になる。$X-\mu\sim N(0,4)$ なら $E(X-\mu)^4=48$ なので、極限分散は $32$ である。分母を $n-1$ に替えても差は $V_n/(n-1)=O_P(n^{-1})$ である。
:::

:::example[弱収束で二乗リスクを代用できない例]
$Y_n=n$ を確率 $1/n$、$0$ を確率 $1-1/n$ で取るとする。$P(|Y_n|>\epsilon)=1/n$ はゼロに近づく。ところが $E Y_n=1$、$E Y_n^2=n$ である。大きい誤差がまれに出る推定量を比べるときには、漸近法則に加えて裾を評価する必要がある。
:::

この差を埋める一つの条件が :term[一様可積分性]{en="uniform integrability" yomi="いちようかせきぶんせい"} である。確率変数列 $H_n$ について $\lim_{K\to\infty}\sup_n E[|H_n|\mathbf1_{|H_n|>K}]=0$ と定める。$H_n\Rightarrow H$ かつこの条件が成立すれば $E H_n\to E H$ である。証明には $[-K,K]$ で切った連続関数を使う。極限側の裾も、連続関数 $(|x|-K)_+\wedge M$ の期待値を収束させ、$M\to\infty$ とすれば $\sup_n E[|H_n|\mathbf1_{|H_n|>K}]$ で抑えられる。例えばある $\delta>0$ に対し $\sup_n E|H_n|^{1+\delta}<\infty$ なら、裾は $K^{-\delta}\sup_n E|H_n|^{1+\delta}$ 以下なので条件を満たす。二乗誤差に使う場合は $H_n$ に拡大誤差の二乗を入れる。

## デルタ法

推定対象を変えると誤差の単位も変わる。滑らかな写像の局所的な傾きが、この変換の近似を与える。

:::theorem[デルタ法]{#thm-stat-delta}
$\sqrt n(T_n-t)\Rightarrow Z$ とし、$g:\R^p\to\R^q$ が $t$ で微分可能とする。$\sqrt n(g(T_n)-g(t))\Rightarrow Dg(t)Z$ である。
:::

:::proof
分布収束による有界性から $T_n-t=o_P(1)$ である。微分の定義による剰余は $\norm{T_n-t}\,r(T_n-t)$ で $\norm{r(v)}\to0$ である。$\sqrt n$ を掛けると $O_P(1)o_P(1)=o_P(1)$ となる。後者は大きい固定の上界を確率的に有界な因子に課してから、小さい因子を確率収束させれば分かる。一次項は分布収束し、:ref[thm-stat-slutsky] で結論を得る。
:::

:::example[比率の対数オッズ]
$0<p<1$ の Bernoulli 標本平均に $g(p)=\log[p/(1-p)]$ を適用する。微分は $1/[p(1-p)]$ なので漸近分散は $1/[p(1-p)]$ になる。$p=1/4$ では $16/3$ である。推定率が $0,1$ のときは変換を任意の有限値に定めても、その事象の確率は固定した内点でゼロに収束する。
:::

:::example[一次導関数がゼロの点]
$\theta=0$ の正規平均の二乗では一次の極限が退化する。$n(\bar X)^2=(\sqrt n\bar X)^2$ にスケールを変えると、分散 $4$ の場合は $4\chi_1^2$ となる。一般に $g'(t)=0$、$g$ が二回微分可能なら二次の Taylor 展開で $n(g(T_n)-g(t))\Rightarrow g''(t)Z^2/2$ を得る。
:::

### 多変量の変換と共分散の伝播

推定対象が二つ以上なら、出力の誤差は共通の入力から生じる。従って精度を求めるには入力の共分散も残す。$\sqrt n(T_n-t)\Rightarrow N_p(0,\Sigma)$ のとき、:ref[thm-stat-delta] から
$$
\sqrt n\{g(T_n)-g(t)\}\Rightarrow
N_q(0,Dg(t)\Sigma Dg(t)^{\mathsf T})
\label{eq-stat-delta-covariance}
$$
となる。正規ベクトル $Z$ に行列 $Dg$ を掛けた共分散は $E[(DgZ)(DgZ)^{\mathsf T}]$ である。$Dg$ の階数が $q$ より小さければ極限分布はある部分空間に集中するため、共分散の逆行列を使う操作には階数の確認が要る。

:::example[同じ標本から得る二つの平均の比]
一観測が $(X,Y)=(1,2),(3,4)$ を等確率で取るとする。平均ベクトルは $(2,3)$、共分散は全成分が $1$ の $2\times2$ 行列である。比 $g(x,y)=x/y$ の微分は $(1/3,-2/9)$ なので、拡大した比の誤差の漸近分散は $(1/3-2/9)^2=1/81$ になる。実際 $Y=X+1$ だから比は $\bar X/(1+\bar X)$ であり、$\bar X=2$ での傾きは $1/9$ である。周辺の分散だけで計算すると交差項 $-4/27$ が落ちる。同じ計算を正の平均 $\mu_X,\mu_Y$ の一般のペアに行うと対数比の分散は
$$
\frac{\sigma_X^2}{\mu_X^2}+\frac{\sigma_Y^2}{\mu_Y^2}
-\frac{2\sigma_{XY}}{\mu_X\mu_Y}
$$
となる。
:::

:::example[標準偏差と変動係数を同時に推定する]
$X\sim N(2,4)$ とする。$T_n=(\bar X,V_n)$ の漸近共分散は $\diag(4,32)$ である。実際、分散の一次項と平均の一次項の共分散は中心三次積率 $0$ である。変換 $g(m,v)=(\sqrt v,\sqrt v/m)$ の $(2,4)$ での行列は
$$
Dg(2,4)=\begin{pmatrix}0&1/4\\-1/2&1/8\end{pmatrix}.
$$
したがって $\sqrt n$ を掛けた二つの誤差の共分散行列は
$$
\begin{pmatrix}2&1\\1&3/2\end{pmatrix}.
$$
正規標本では平均と分散が独立でも、変動係数は両者を使う。このため標準偏差と変動係数の極限には相関が残る。
:::

### 二次変換と微分できない変換

二回連続微分可能な実数値関数 $g$ について $Dg(t)=0$ なら、二次のデルタ法は
$$
n\{g(T_n)-g(t)\}\Rightarrow \tfrac12 Z^{\mathsf T}D^2g(t)Z
$$
を与える。積分形 Taylor 展開で、二次形式の係数は線分上の Hessian を重み付きで平均したものである。連続性と $T_n\to_p t$ から、その平均は真値の Hessian の半分に確率収束する。$\sqrt n(T_n-t)=O_P(1)$ なので、係数の差を掛けた剰余は消える。$Z\sim N(0,\Sigma)$ なら、$Z=\Sigma^{1/2}W$、$W\sim N(0,I)$ と書ける。:ref[thm-linalg-real-spectral] による対角化から、二次形式の極限は $\tfrac12\Sigma^{1/2}D^2g(t)\Sigma^{1/2}$ の固有値を重みにした独立な $\chi_1^2$ の和になる。

:::example[二つの方向で符号が異なる二次近似]
独立な二成分を持つ $T_n$ について $\sqrt n T_n\Rightarrow N(0,\diag(2,1))$ とする。$g(x,y)=x^2-3y^2$ なら Hessian は $\diag(2,-6)$ なので $n g(T_n)\Rightarrow 2U-3V$ である。ただし $U,V$ は独立な $\chi_1^2$ である。極限の平均は $-1$、分散は $2(2^2+3^2)=26$ となる。負の重みがあるため、この法則は通常のカイ二乗分布の定数倍では表せない。
:::

非滑らかな場合にも、局所で一次同次な近似があれば結果が得られる。連続な写像 $d$ があり、任意の $K$ について
$$
\sup_{\norm h\le K}\left\|\frac{g(t+rh)-g(t)}r-d(h)\right\|\to0
\quad(r\downarrow0)
$$
とする。$r=n^{-1/2}$、$h=\sqrt n(T_n-t)$ を代入し、$\norm h\le K$ の外を緊密性で処理すると、拡大誤差と $d(h)$ の差は $o_P(1)$ である。従って連続写像と Slutsky の補題で極限は $d(Z)$ になる。

:::example[二案の良い方を報告する]
$g(x,y)=\max(x,y)$、$t=(0,0)$ とする。局所写像は $d(h)=\max(h_1,h_2)$ そのものである。入力の正規極限が独立な標準正規二成分なら、拡大誤差の分布関数は $P(Z_1\le z,Z_2\le z)=\Phi(z)^2$ になる。特に $0$ 以下の確率は $1/4$ であり、対称な正規近似を当てることはできない。真値の一方が他方より大きい場合には、真値の近傍でその成分を選ぶ線形写像になる。
:::

## 最尤推定量の一致性

最尤推定量が真値に近づくには、期待対数尤度の最大点だけでなく、標本上での対数尤度が同時に近いことが要る。

:::theorem[最大点の一致性]{#thm-stat-argmax}
任意の $\epsilon>0$ に対し $M(\theta_0)-\sup_{\norm{\theta-\theta_0}\ge\epsilon}M(\theta)>0$ とする。$\sup_\theta|M_n(\theta)-M(\theta)|\to_p0$ かつ $M_n(\widehat\theta_n)\ge M_n(\theta_0)-o_P(1)$ なら $\widehat\theta_n\to_p\theta_0$ である。
:::

:::proof
$M(\theta_0)-M(\widehat\theta_n)$ は高々 $2\sup|M_n-M|+o_P(1)$ である。一方、推定値が $\epsilon$ 以上離れる事象では左辺は正定数以上である。右辺がその定数を超える確率はゼロに収束する。
:::

各 $\theta$ で $E_{\theta_0}|\log p_\theta(X)|<\infty$ と仮定し、$M(\theta)=E_{\theta_0}\log p_\theta(X)$ とおけば差 $M(\theta_0)-M(\theta)$ は相対エントロピーである。$-\log u\ge1-u$ により非負で、共通台と識別可能性のもとで $\theta=\theta_0$ に限りゼロになる。コンパクトな $\Theta$ で $M$ が連続なら任意の閉補集合上に正の差がある。例えば $|\log p_\theta-\log p_{\vartheta}|\le L(X)\norm{\theta-\vartheta}$、$EL<\infty$ なら、一様大数の法則は有限の $\delta$ 網で示せる。網上の誤差は有限個の大数の法則でゼロに近づき、網からの誤差は $\delta(P_n L+E L)$ 以下なので $\delta$ を小さくできる。

:::example[コンパクト区間上の Bernoulli]
$p\in[1/4,3/4]$ なら各 $\log p_\theta(X)$ の導関数は絶対値 $4$ 以下である。標本率 $s/n=1/8$ のとき非制約最大点は区間外なので最尤値は端の $1/4$ となる。一致性には最大点が標本ごとに内点であることは不要である。
:::

:::example[端点モデルの直接計算]
一様分布の最大値 $M$ では $P(|M-\theta|>\epsilon)=(1-\epsilon/\theta)^n$（$0<\epsilon<\theta$）である。例えば $\epsilon=\theta/2$ なら $2^{-n}$ になる。この場合は対数尤度が台の外で $-\infty$ になるため、上の一様近似の条件を直接使わずに証明できる。
:::

### 一様大数の法則を積分可能な包絡から得る

連続な目的関数を多くのパラメータで同時に平均するとき、各点での大数の法則だけでは最大化に必要な精度を制御できない。パラメータ集合を有限個の近傍で覆い、その中での関数の変動も平均できるようにする。

:::lemma[コンパクトな族の一様大数の法則]{#lem-stat-compact-ulln}
$\Theta\subset\R^d$ はコンパクト、$m_\theta(x)$ は $x$ について可測で、ほとんどすべての $x$ で $\theta$ について連続とする。$\sup_{\theta\in\Theta}|m_\theta(X)|\le H(X)$ かつ $EH<\infty$ とする。独立同分布標本について
$$
\sup_{\theta\in\Theta}|P_nm_\theta-Pm_\theta|\longrightarrow0
\quad\text{ほぼ確実}
$$
である。ここで $Pm_\theta=E m_\theta(X)$ である。
:::

:::proof
$\omega_\delta(x)=\sup_{\norm{\theta-\vartheta}\le\delta}|m_\theta(x)-m_\vartheta(x)|$ とおく。連続性から上限は添字の集合の可算稠密部分でも求められ、可測である。各 $x$ で一様連続性により $\omega_\delta(x)\to0$、また $\omega_\delta(x)\le2H(x)$ なので $P\omega_\delta\to0$ である。

$\delta$ 網 $\theta_1,\ldots,\theta_J$ を取り、各 $\theta$ を距離 $\delta$ 以下の網の点に対応させる。三角不等式を適用すると
$$
\sup_\theta |(P_n-P)m_\theta|
\le \max_{1\le j\le J}|(P_n-P)m_{\theta_j}|+P_n\omega_\delta+P\omega_\delta.
$$
独立同分布の強大数の法則で右辺の上極限は $2P\omega_\delta$ 以下である。例えば $\delta=1/k$ の可算列についてこの事象を交差させても確率 $1$ である。そこで $k\to\infty$ とすると結論を得る。
:::

:::example[指数分布でコンパクトな制約を使う]
率 $\lambda$ を $[1,4]$ に制限し、$m_\lambda(x)=\log\lambda-\lambda x$ とする。$x\ge0$ では $|m_\lambda(x)|\le \log4+4x$ であり、真の率 $\lambda_0>0$ の下で積分可能である。$M(\lambda)=\log\lambda-\lambda/\lambda_0$ の導関数は $1/\lambda-1/\lambda_0$、二階導関数は $-1/\lambda^2$ である。例えば真の率が $2$ なら $M(2)-M(1)=\log2-1/2>0$ で、最大点は $2$ である。一様大数の法則と :ref[thm-stat-argmax] により、制約付き最尤推定量は $2$ に収束する。直接には $1/\bar X$ を $[1,4]$ に射影しても同じ結論を得る。
:::

### 最大点の分離と非コンパクトなパラメータ空間

局所的な一意性と、遠い点で最大値に近づけないことは別の条件である。例えば目的関数 $M(\theta)=-\theta^2/(1+\theta^4)$ は $0$ で唯一最大となるが、$|\theta|\to\infty$ でも $M(\theta)\to0$ である。$M_n=M$、$\widehat\theta_n=n$ とすれば近似最大化の誤差は $O(n^{-2})$ だが、推定値は発散する。従って近似最適化と各点での最大点の一意性から直ちに一致性は従わない。

非コンパクトな場合にはまず $\widehat\theta_n$ を高確率で固定コンパクトに閉じ込める。例えば各 $\eta>0$ で真値を含むコンパクト $K_\eta$ があり、$\limsup_n P(\widehat\theta_n\notin K_\eta)\le\eta$ とする。各 $K_\eta$ の上だけで :ref[thm-stat-argmax] の一様収束と分離を仮定すれば、誤差事象の上極限は $\eta$ 以下になる。実際、$\widehat\theta_n\in K_\eta$ 上ではその定理の証明と同じ不等式が成立し、外の確率は $\eta$ で抑えられる。$\eta$ は任意なので一致性が従う。

:::example[正規の未知分散を上下から抑える]
正規モデルで平均 $\mu$ と分散 $v>0$ を未知とする。$\bar X$ と $V_n=n^{-1}\sum_i(X_i-\bar X)^2$ を用いると、標本平均対数尤度は
$$
-\tfrac12\log(2\pi v)-\frac{V_n+(\bar X-\mu)^2}{2v}
$$
である。各 $v$ で $\mu=\bar X$ が最大なので、残る一変数の微分は $(V_n-v)/(2v^2)$ となる。$V_n>0$ なら $v=V_n$ が大域的最大である。真値が $(\mu_0,v_0)=(2,4)$ なら大数の法則で $(\bar X,V_n)$ は確率 $1$ に近づいて $[1,3]\times[2,6]$ に入る。$n\ge2$ の連続正規標本では $V_n=0$ の確率はゼロであるが、同じ値が並ぶデータを一般の数値入力として扱う場合、$v\downarrow0$ で尤度が発散する。この場合は最大点の存在を主張できない。
:::

### モデルが真の分布を含まないとき

独立標本の真の分布を $P$ とし、真の密度 $p$ が候補 $p_\theta$ の族に属さないこともある。$P\log p$ と $P\log p_\theta$ が有限なら
$$
P\log p_\theta=P\log p-D(P\|P_\theta),
\qquad D(P\|Q)=\int p\log(p/q)\,\dd\nu
$$
である。従って一様収束と分離が成立する限り、最尤推定値は相対エントロピーを最小化する点へ近づく。この点を :term[擬似真値]{en="pseudo-true parameter" yomi="ぎじしんち"} という。候補について期待対数密度が定義できるなら、$P$ 自体に密度がない場合でも最大化の議論は使える。

:::example[一様分布を正規分布で近似する]
真の分布を $\operatorname{Unif}(0,4)$ として、正規密度 $N(\mu,v)$ の期待対数密度を最大化する。$E X=2$、$E(X-2)^2=4/3$ なので、$E(X-\mu)^2=4/3+(\mu-2)^2$ である。まず $\mu=2$ が全 $v$ で最大を与え、次に $v$ の微分 $(4/3-v)/(2v^2)$ から $v=4/3$ が最大を与える。標本での最尤値が $\bar X,V_n$ になる代数は正規性に依存しないので、両者は大数の法則でも同じ点へ近づく。このモデルで得た正規分布は平均と分散を再現するが、観測の台 $[0,4]$ は再現しない。
:::

## 最尤推定量の漸近正規性

最大点が内点で真値に近いなら、スコア方程式を真値の周囲で展開できる。一次項の分布と二次微分の確率極限を別々に求める。

:::theorem[滑らかな最尤推定量の漸近正規性]{#thm-stat-mle-normal}
独立同分布モデルで真値 $\theta_0$ は内点とする。一致する最尤推定量が確率 $1$ に近づく事象で内点にあり、真値の近傍で確率密度は正の共通台を持ち、その台で対数密度は二回連続微分可能とする。スコアは平均ゼロ、共分散 $I>0$ とし、標本平均 Hessian が近傍で期待 Hessian に一様確率収束し、期待 Hessian は真値で連続で値は $-I$ とする。このとき
$$
\sqrt n(\widehat\theta_n-\theta_0)=I^{-1}n^{-1/2}\sum_i s_{\theta_0}(X_i)+o_P(1)\Rightarrow N(0,I^{-1}).
\label{eq-stat-mle-expansion}
$$
:::

:::proof
内点では $\nabla\ell_n(\widehat\theta_n)=0$ である。線分に沿う積分形の Taylor 展開により $0=n^{-1/2}\nabla\ell_n(\theta_0)+H_n\sqrt n(\widehat\theta_n-\theta_0)$ となる。$H_n$ はこの線分上の標本平均 Hessian の平均である。一致性、一様収束、期待値の連続性から $H_n\to_p-I$ であり、正則性より $H_n^{-1}\to_p-I^{-1}$ である。中心極限定理と Slutsky の補題で結論を得る。
:::

:::example[Poisson 推定量]
$\widehat\lambda=\bar X$、$I=1/\lambda$ なので漸近分散は $\lambda$ である。$\lambda=4$ では $n^{1/2}(\widehat\lambda-4)$ の極限標準偏差は $2$ である。
:::

:::example[境界のある正規平均]
$\theta\ge0$ の正規標本で最尤値は $\max(0,\bar X)$ である。真値 $0$ では $\sqrt n\widehat\theta=\max(0,Z)$、ただし $Z$ は標準正規である。従って$0$ に質量 $1/2$ を持つ。内点条件を除くと正規近似は変わる。
:::

### スコアの共分散と負の期待 Hessian

:ref[thm-stat-mle-normal] の行列 $I$ が二通りの式で現れる理由を確認する。ある近傍で密度が正の共通の台を持ち、その台で $p_\theta$ を二回微分できるとする。さらに各偏微分を積分可能な関数で一様に抑えられるとする。すると $\int p_\theta\,\dd\nu=1$ を一回微分し、$s_\theta=\nabla\log p_\theta$ と置くと
$$
0=\int\nabla p_\theta\,\dd\nu=E_\theta s_\theta.
$$
二回微分については
$$
\nabla^2p_\theta
=p_\theta\{\nabla^2\log p_\theta+s_\theta s_\theta^{\mathsf T}\}
$$
なので、積分を行い $-E_\theta\nabla^2\log p_\theta=E_\theta s_\theta s_\theta^{\mathsf T}$ を得る。期待 Hessian の連続性と標本平均の一様収束は、例えば :ref[lem-stat-compact-ulln] の条件を成分ごとに課せば保証される。動く台を持つ模型では、台の端から微分の項が出る場合がある。従って密度が正の点で微分式を書けることだけでは、期待値と微分を交換できない。

:::example[平均と分散を別々に推定する]
$N(\mu,v)$ の一観測のスコアは、$U=X-\mu$ として
$$
s_\mu=U/v,\qquad s_v=(U^2-v)/(2v^2)
$$
である。$EU^2=v$、$EU^3=0$、$EU^4=3v^2$ から共分散は $\diag(1/v,1/(2v^2))$ になる。真値が $(2,4)$ なら逆行列は $\diag(4,32)$ である。推定量の一階展開から、平均の誤差は $n^{-1}\sum U_i$、分散の誤差は $n^{-1}\sum(U_i^2-v)+o_P(n^{-1/2})$ となる。前者は厳密な恒等式であり、後者では引かれる項 $(\bar X-\mu)^2$ が $O_P(n^{-1})$ である。スコアからの展開と標本積率からの直接計算は一致する。
:::

### 観測情報と期待情報を標準誤差に使う

標本平均の負の Hessian $\widehat I_{\mathrm{obs}}=-n^{-1}\nabla^2\ell_n(\widehat\theta_n)$ を :term[観測情報量]{en="observed information" yomi="かんそくじょうほうりょう"} と呼ぶことがある。期待情報 $I(\widehat\theta_n)$ と区別するため、ここではどちらも一観測当たりとする。期待情報が真値で連続なら、後者は一致性で $I(\theta_0)$ に収束する。前者は Hessian の一様確率収束と連続性から同じ極限を持つ。正定値行列の逆と正定値平方根は正定値集合で連続である。逆は $A^{-1}=(\det A)^{-1}\operatorname{adj}A$ で確認でき、平方根の連続性は固有値を含む固定正区間上で平方根関数を多項式で一様近似すれば従う。従って
$$
\widehat I_{\mathrm{obs}}^{1/2}\sqrt n(\widehat\theta_n-\theta_0)
\Rightarrow N_d(0,I_d)
$$
となる。区間を作る場合は $n^{-1}\widehat I_{\mathrm{obs}}^{-1}$ の対角要素が個々の係数の分散の近似になる。

:::example[幾何分布の成功確率]
$P_p(X=x)=p(1-p)^{x-1}$、$x=1,2,\ldots$ とする。$\ell_n(p)=n\log p+(\sum_i X_i-n)\log(1-p)$ なので、$\bar X>1$ なら $\widehat p=1/\bar X$ である。負の一観測 Hessian は $1/p^2+(X-1)/(1-p)^2$ である。$E_p(X-1)=(1-p)/p$ を入れると $I(p)=1/[p^2(1-p)]$ となる。$\bar X=4$ なら $\widehat p=1/4$、観測情報と代入期待情報はいずれも $64/3$ である。$n=64$ での標準誤差は $\sqrt3/64$ となる。全観測が $1$ なら $\widehat p=1$ は閉区間での境界解であり、この事象の確率は $n=64$ なら $p^{64}$ であり、固定した $p<1$ のもとで標本数を増やすとゼロへ向かう。
:::

### 攪乱パラメータを含むときの情報の逆

$\theta=(\psi,\lambda)$ で $\psi$ を推定したいとする。情報行列を
$$
I=\begin{pmatrix}I_{\psi\psi}&I_{\psi\lambda}\\ I_{\lambda\psi}&I_{\lambda\lambda}\end{pmatrix}
$$
と分ける。展開 :eqref[eq-stat-mle-expansion] の連立方程式から $\lambda$ の行を使って $\lambda$ の誤差を消去すると
$$
\sqrt n(\widehat\psi-\psi_0)
=J^{-1}\left\{\Delta_{\psi,n}
-I_{\psi\lambda}I_{\lambda\lambda}^{-1}\Delta_{\lambda,n}\right\}+o_P(1),\qquad
J=I_{\psi\psi}-I_{\psi\lambda}I_{\lambda\lambda}^{-1}I_{\lambda\psi}
$$
となる。括弧の極限共分散を展開すると $J$ になるため、対象だけの極限共分散は $J^{-1}$ である。$\lambda$ 既知の場合の $I_{\psi\psi}^{-1}$ と比べると、独立に動かせる $\lambda$ の方向を除いた分だけ情報が減る。$I_{\psi\lambda}=0$ なら両方は同じである。

:::example[二つの平均の線形結合を推定する]
一観測で独立な $(X,Y)$ が $N(\psi+\lambda,1)$ と $N(\lambda,1)$ に従うとする。スコアは $(U,U+V)$、ただし $U=X-\psi-\lambda,V=Y-\lambda$ と書ける。$I_{\psi\psi}=1$、$I_{\lambda\lambda}=2$、$I_{\psi\lambda}=1$ なので、$J=1/2$ である。従って $\psi$ の最尤推定量の極限分散は $2$ となる。実際 $\widehat\psi=\bar X-\bar Y$ の有限標本分散は $2/n$ である。$\lambda$ を知っていれば $\bar X-\lambda$ の分散は $1/n$ で済む。差の推定を追加した代償が行列の引き算に対応する。
:::

## M 推定量と Z 推定量の漸近論

尤度以外の損失を最小にした場合も、微分で方程式が得られる。その方程式の変動と傾きが、一般には別の行列になる。

:::theorem[滑らかな Z 推定量]{#thm-stat-z-estimator}
$X_1,\ldots,X_n$ は独立同分布で、$P_n g_\theta=n^{-1}\sum_i g_\theta(X_i)$ とする。真値 $\theta_0$ の凸開近傍で $g_\theta$ は $\theta$ についてほぼ確実に連続微分可能とする。$\widehat\theta$ が一致し、$P_ng_{\widehat\theta}=o_P(n^{-1/2})$、$E g_{\theta_0}=0$、$B=E g_{\theta_0}g_{\theta_0}^{\mathsf T}$ が有限とする。$\theta_0$ と $\widehat\theta$ を結ぶ線分上の標本平均 Jacobian を線分に沿って積分した行列が正則な $A$ に確率収束するなら、$\sqrt n(\widehat\theta-\theta_0)\Rightarrow N(0,A^{-1}B A^{-\mathsf T})$ である。$A^{-\mathsf T}$ は $A^{-1}$ の転置を表す。
:::

:::proof
$\sqrt nP_ng_{\widehat\theta}=\sqrt nP_ng_{\theta_0}+A_n\sqrt n(\widehat\theta-\theta_0)$ と積分形で展開する。$A_n$ の逆を掛けて、中心極限定理と :ref[thm-stat-slutsky] を適用する。一致性を持つ滑らかな内点 M 推定量では、$g_\theta$ に損失の勾配を選ぶとこの式が使える。
:::

:::example[重み付き平均]
$g_\theta=W(X-\theta)$ とし、$W$ は $X$ と独立で $P(W=1)=P(W=3)=1/2$、$\operatorname{Var}X=1$ とする。$A=-2,B=5$ なので漸近分散は $5/4$ である。重みが定数でないため、単純平均の漸近分散 $1$ より大きい。
:::

:::example[Huber 位置の分散]
閾値 $c>0$ の Huber スコア $\psi_c(u)=\max(-c,\min(c,u))$ を用い、$g_\theta(X)=\psi_c(X-\theta)$ とおく。独立同分布の観測は真の位置 $\theta_0$ に関して対称で、連続な分布を持つとする。$\psi_c$ は奇関数なので真値での期待スコアはゼロになる。差商は $1$ で有界で、$|X-\theta|=c$ は零確率だから、支配収束により期待スコアの傾きは $-P(|X-\theta|<c)$ である。$P(|X-\theta_0|<c)>0$ なら、$\theta_0$ のわずかに左と右で期待スコアは正と負になる。標本平均スコアは連続かつ非増加であり、大数の法則から高確率でこれら二点の間にすべての根が入る。両点を任意に真値へ近づけられるので、根は一致する。

標本の勾配は折点で変わるが、各スコアは絶対連続なので線分上の積分による展開は成立する。根が $\theta_0$ から $\epsilon$ 以内なら、線分上の負の勾配の平均は $|X_i-\theta_0|<c-\epsilon$ と $|X_i-\theta_0|\le c+\epsilon$ の標本頻度で挟める。大数の法則と $\epsilon\downarrow0$ により傾きの極限が得られる。スコアの分散は $E\min((X-\theta_0)^2,c^2)$ である。例えば一様誤差 $[-2,2]$ で $c=1$ なら分子は $1/6+1/2=2/3$、分母は $(1/2)^2$ なので漸近分散は $8/3$ である。
:::

### 一段階推定

計算しやすい一致推定値 $\widetilde\theta$ がすでに $\sqrt n$ の収束率を持つなら、推定式を一度線形化するだけで、その根と同じ一次展開を得られる。$A_n$ を $\theta_0$ と $\widetilde\theta$ の間の線分に沿って積分した標本平均 Jacobian とし、$A_n\to_p A$ とする。$\widehat A\to_p A$ は反復に使う傾きである。$A$ が正則なら :term[一段階推定量]{en="one-step estimator" yomi="いちだんかいすいていりょう"} を
$$
\theta^{(1)}=\widetilde\theta-\widehat A^{-1}P_ng_{\widetilde\theta}
$$
と定められる。逆行列が存在しない事象では任意の有限値を与える。線分に沿って $P_ng_{\widetilde\theta}$ を展開すると
$$
\sqrt n(\theta^{(1)}-\theta_0)
=-\widehat A^{-1}\sqrt nP_ng_{\theta_0}
+(I_d-\widehat A^{-1}A_n)\sqrt n(\widetilde\theta-\theta_0)
$$
となる。右端の行列は $o_P(1)$ であり、それを $\sqrt n$ の収束率を持つ誤差に掛けた項は消える。従って :ref[thm-stat-z-estimator] の根と一次展開が一致する。スコアを使うなら $A=-I$ であり、傾きに負の期待情報を代入できる。

:::example[指数分布の率を一度の補正で求める]
率 $\lambda$ の指数標本で、標本中央値 $m_n$ から $\widetilde\lambda=\log 2/m_n$ を作る。中央値の正規近似は次項で示す。この推定値は $\sqrt n$ の収束率を持つ。スコアの標本平均は $1/\lambda-\bar X$、傾きは $-1/\lambda^2$ なので、一段階推定値は $2\widetilde\lambda-\widetilde\lambda^2\bar X$ になる。有限標本で正でない場合は元の $\widetilde\lambda$ を使う。一致性からこの置換の確率はゼロに近づく。実際、通常の最尤値 $1/\bar X$ との差は $-\bar X(\widetilde\lambda-1/\bar X)^2=O_P(n^{-1})$ である。両方の推定値を $\sqrt n$ 倍して比較すると、その差は消える。
:::

### 標本分位点の漸近正規性

:::proposition[標本分位点の正規極限]{#prop-stat-sample-quantile-normal}
連続な分布関数 $F$ からの独立同分布標本を取り、$F$ は $q$ で微分可能、$F(q)=p\in(0,1)$、$f(q)=F'(q)>0$ とする。$k_n=\lceil np\rceil$ と置くと、$\widehat q=X_{(k_n)}$ に対し
$$
\sqrt n(\widehat q-q)\Rightarrow N(0,p(1-p)/f(q)^2)
$$
である。
:::

:::proof
固定 $t$ に対し、$\sqrt n(\widehat q-q)\le t$ は閾値 $q+t/\sqrt n$ 以下の観測個数 $S_n$ が $k_n$ 以上であることと同値である。成功確率 $p_n=F(q+t/\sqrt n)=p+f(q)t/\sqrt n+o(n^{-1/2})$ の二項変数なので、分散は $np_n(1-p_n)$ である。各指示変数は有界であり、$p_n$ は内点に近づくので :ref[thm-prob-lindeberg-feller] の条件が成立する。標準化した境界 $(k_n-np_n)/\sqrt{np_n(1-p_n)}$ は $-f(q)t/\sqrt{p(1-p)}$ に収束する。標準正規の分布関数を $\Phi$ とすると、その境界以上の確率は $\Phi(f(q)t/\sqrt{p(1-p)})$ に収束する。これは主張の正規分布の分布関数である。
:::

損失が絶対値で滑らかでない標本中央値も、上の二項確率の計算で扱える。

:::example[正規中央値と平均の精度]
$N(\mu,1)$ の中央値では $p=1/2$、$f(\mu)=1/\sqrt{2\pi}$ なので漸近分散は $\pi/2$ である。平均の漸近分散は $1$ なので、正規分布内で中央値を使うと効率は $2/\pi$ になる。一方、極端な観測一つを変えた場合の中央値の応答は、観測順位の中央付近に限られる。効率の評価と頑健性の評価には異なる比較が要る。
:::

### 過不足のある方程式と平滑化の役割

$P_ng_\theta=0$ を解くとき、未知成分と方程式が同数でも厳密な解が存在するとは限らない。:ref[thm-stat-z-estimator] の残差条件は、数値計算による近似値が統計誤差より小さい残差を持つことを意味する。例えば $g_\theta(x)=x-\theta$ なら、残差は $\bar X-\widehat\theta$ そのものである。固定された許容誤差で近似を打ち切ると、その誤差は $n$ を増やしても小さくならない。漸近分布が必要であれば、反復後の残差を $o_P(n^{-1/2})$ に抑える。逆に残差が $b/\sqrt n$ に漸近するなら、定理の展開から極限中心が $A^{-1}b$ だけずれる。

方程式の方が多い場合には、重み付き二乗で残差を最小にする。$g_\theta(X)\in\R^k$、$\theta\in\R^d$、$k\ge d$ とし、独立同分布標本について $Pg_{\theta_0}=0$、$P\norm{g_{\theta_0}}^2<\infty$ とする。$A=P[\nabla_\theta g_{\theta_0}]$ は列フルランクとする。$W_n\to_p W>0$ とし、内点にある一致した最小点で $(P_n\nabla g_\theta)^{\mathsf T} W_nP_ng_\theta=0$ が成立すると仮定する。評価点と真値を結ぶ線分の平均 Jacobian と、最小点での Jacobian はともに $A$ に確率収束するとする。すると
$$
\sqrt n(\widehat\theta-\theta_0)=
-(A^{\mathsf T}WA)^{-1}A^{\mathsf T}W\sqrt nP_ng_{\theta_0}+o_P(1).
$$
実際、方程式の後半の因子を $P_ng_{\theta_0}+A_n(\widehat\theta-\theta_0)$ に展開し、左側の行列を掛けて解く。$A^{\mathsf T}WA$ は正定値なので逆行列は確率収束する。さらに標本平均の中心極限定理で極限法則を得る。$B=\operatorname{Cov}(g_{\theta_0})>0$ なら、$W=B^{-1}$ の選択がこの形の共分散を最小にする。実際 $D=(A^{\mathsf T}WA)^{-1}A^{\mathsf T}W$、$D_*=(A^{\mathsf T}B^{-1}A)^{-1}A^{\mathsf T}B^{-1}$ とおくと $DA=D_*A=I$ である。この制約により $(D-D_*)BD_*^{\mathsf T}=0$ となるため
$$
DBD^{\mathsf T}-D_*BD_*^{\mathsf T}=(D-D_*)B(D-D_*)^{\mathsf T}\ge0
$$
である。これが :term[一般化モーメント法]{en="generalized method of moments" yomi="いっぱんかもーめんとほう"} の重みの選び方である。

:::example[精度の違う二つの測定を組み合わせる]
同一の平均 $\theta$ を持つ独立な測定 $(X,Y)$ の分散を $1,4$ とする。$g=(X-\theta,Y-\theta)^{\mathsf T}$ なので $A=(-1,-1)^{\mathsf T}$、$B=\diag(1,4)$ である。$W=I_2$ なら $\widehat\theta=(\bar X+\bar Y)/2$ で漸近分散 $5/4$ である。$W=B^{-1}$ では $\widehat\theta=(4\bar X+\bar Y)/5$ となり、分散は $16/25+4/25=4/5$ になる。この例では重みが固定なので、両方の分散は有限標本でも $1/n$ 倍で得られる。
:::

### 局所二次展開による M 推定

損失の勾配が標本ごとに存在しない場合でも、目的関数を局所で二次式に近似し、最大点を比較できる。$M_n$ を最大化し、$U_n=\sqrt n(\widehat\theta_n-\theta_0)=O_P(1)$ とする。固定した任意の $K$ に対し $\norm h\le K$ で一様に
$$
n\{M_n(\theta_0+h/\sqrt n)-M_n(\theta_0)\}
=h^{\mathsf T}G_n-\tfrac12h^{\mathsf T}Ah+o_P(1),
\quad G_n\Rightarrow N(0,B),\quad A>0
$$
とする。真値に開近傍があり、最大化の誤差は $o_P(1/n)$ と仮定する。このとき $U_n=A^{-1}G_n+o_P(1)$ である。証明には、$U_n$ と二次式の最大点 $h_n^*=A^{-1}G_n$ の両方を高確率で同じ球へ入れる。二次式での最大値と $U_n$ での値の差は平方完成で $\tfrac12(U_n-h_n^*)^{\mathsf T}A(U_n-h_n^*)$ である。一様展開と最大化誤差でこの差は $o_P(1)$ となる。$A$ の最小固有値が正なので結論を得る。この方法でも、先に拡大誤差の緊密性を証明する条件は別に要る。

:::example[正規の分散モデルを誤指定した場合のサンドイッチ分散]
真の分布を $\operatorname{Unif}(0,4)$ とし、正規対数尤度で $(\mu,v)$ を推定する。擬似真値は $(2,4/3)$ である。$U=X-2$ とおくと $EU^4=\int_{-2}^2u^4\,\dd u/4=16/5$ である。正規モデルの最尤点はやはり $\bar X,V_n$ なので、実際の漸近共分散は
$$
\diag\!\left(\tfrac43,\tfrac{16}5-(\tfrac43)^2\right)
=\diag\!\left(\tfrac43,\tfrac{64}{45}\right)
$$
となる。一方、モデル内の Fisher 情報の逆に代入すると分散成分は $2(4/3)^2=32/9$ である。ずれの原因は、真の分布で計算したスコア共分散 $B$ と負の期待 Hessian が一致しないことである。独立標本ならスコアの外積を標本平均し、Hessian の逆を両側に掛けた推定で実際の漸近分散を推定できる。必要な一様収束は、:ref[lem-stat-compact-ulln] をスコアの外積と Hessian に適用して確認する。
:::

## 漸近有効性と超有効性

一つの真値に絞れば、そこへの縮小によって誤差を極端に小さくできる。その利得が $1/\sqrt n$ だけ真値を動かしても残るかを調べる。

:::definition[正則性と漸近有効性]{#def-stat-regular-estimator}
推定量列 $T_n$ は、各固定 $h$ に対し $\theta_n=\theta_0+h/\sqrt n$ のもとで $\sqrt n(T_n-\theta_n)$ が $h$ に依存しない法則に収束すれば、$\theta_0$ で :term[正則]{en="regular" yomi="せいそく"} という。正則で極限が $N(0,I^{-1})$ のとき漸近有効という。
:::

:::theorem[畳み込み下界]{#thm-stat-convolution-bound}
真値で二次平均微分可能な独立同分布モデルとし、$\theta_0$ は内点、二次平均微分のスコアの共分散 $I$ は正定値とする。二次平均微分可能とは平方根密度の $L^2$ 展開の線形項が $h^{\mathsf T}s_{\theta_0}\sqrt{p_{\theta_0}}/2$ で、残差ノルムが $o(\norm h)$ となることである。$\theta$ 自体を推定する正則推定量の極限法則は $N(0,I^{-1})$ とある確率分布との畳み込みである。証明は最終節で、二次平均微分から得る尤度比の展開と第三補題を用いて与える。
:::

:::example[Hodges の推定量]
$X_i\sim N(\theta,1)$ で $|\bar X|\le n^{-1/4}$ なら $T_n=0$、他で $T_n=\bar X$ とする。固定した $\theta=0$ では $P(T_n=0)\to1$ なので拡大誤差は $0$ に確率収束する。固定した $\theta\ne0$ では確率 $1$ に近づいて $\bar X$ と一致する。一方 $\theta_n=h/\sqrt n$ では縮小される確率が $P(|Z+h|\le n^{1/4})\to1$ なので中心化誤差は $-h$ に収束する。例えば $h=2$ では極限値が $-2$ であり、正則性を失う。
:::

:::example[正則な正規標本平均]
同じモデルで $\sqrt n(\bar X-\theta_n)$ は各 $n,h$ で $N(0,1)$ である。$I=1$ なので正則かつ漸近有効である。定数 $a/\sqrt n$ を加えると正則性は保つが、中心が $a$ の正規極限に変わる。
:::

### 漸近線形な推定量の下界を直接求める

:eqref[eq-stat-lan] の条件の下で、真値において $\sqrt n(T_n-\theta_0)=n^{-1/2}\sum_i a(X_i)+o_P(1)$、$E_{\theta_0}a=0$、$E_{\theta_0}\norm a^2<\infty$ とする。多変量中心極限定理でこの誤差と $\Delta_n$ は結合正規極限を持つ。局所対立に移ると、第三補題により平均が $Ch$ だけ動く。ここで $C=E[as_{\theta_0}^{\mathsf T}]$ である。$\theta_0+h/\sqrt n$ で中心化すると平均は $(C-I_d)h$ になる。従って正則性には $C=I_d$ が必要である。$E[(a-I^{-1}s_{\theta_0})(a-I^{-1}s_{\theta_0})^{\mathsf T}]=\operatorname{Cov}(a)-I^{-1}$ は半正定値なので、情報量の逆がこのクラスでの漸近分散の下界になる。等号は $a=I^{-1}s_{\theta_0}$ がほぼ確実に成り立つ場合である。この議論では畳み込み定理のうち漸近線形な規則を直接扱っている。

:::example[非線形な誤差成分を加えた推定量]
$N(\theta,1)$ で $T_n=\bar X+c(S^2-1)$ とする。$S^2$ を真の平均の周りで展開すると、一次の項は $a(X)=X-\theta+c((X-\theta)^2-1)$ である。スコア $X-\theta$ との積の期待値は $1$、$a$ の分散は $1+2c^2$ である。局所対立での平均は正しく移動するが、$c\ne0$ なら情報量の逆より大きい分散が残る。
:::

### 関数を推定する場合の有効な一次項

$\psi(\theta)\in\R^k$ が真値で微分可能とし、その微分を $D$ とする。正則性は中心を $\psi(\theta_0+h/\sqrt n)$ とした極限の $h$ からの独立性で定義する。漸近線形な推定量の一次項を $a(X)$ とすれば、局所代替での極限平均の移動は $E[a s^{\mathsf T}]h$ である。一方で中心の移動は $Dh$ であるため、正則性には $E[a s^{\mathsf T}]=D$ が必要である。従って $a_{\mathrm{eff}}=DI^{-1}s$ を引くと
$$
\operatorname{Cov}(a)-DI^{-1}D^{\mathsf T}
=E[(a-a_{\mathrm{eff}})(a-a_{\mathrm{eff}})^{\mathsf T}]
$$
となる。この計算では両方の交差項に $E[a s^{\mathsf T}]=D$ を使っている。:ref[thm-stat-mle-normal] のもとで $\psi(\widehat\theta)$ はデルタ法により一次項が $a_{\mathrm{eff}}$ になり、局所代替での移送が成立する場合にはこの下界を達成する。

:::example[対数率の有効な分散]
$X\sim\operatorname{Pois}(\lambda)$ で $\psi(\lambda)=\log\lambda$ とする。$I=1/\lambda$、$D=1/\lambda$ なので有効な極限分散は $1/\lambda$ である。影響関数 $DI^{-1}s$ は $(X-\lambda)/\lambda$ である。例えば $\lambda=4$ では極限標準偏差は $1/2$ となる。$h=2$ の局所対立で元の $\log4$ から中心化した誤差の平均は $h/\lambda=1/2$ だけ移動する。$P(\bar X=0)=e^{-n\lambda}$ なので、その事象で有限値を任意に置いても極限は変わらない。
:::

### 正則性と局所最大リスク

畳み込み下界は正則な推定量の極限法則を比較する。一方、真値付近の最大リスクを使えば、正則性を失う推定量も評価できる。まず $X_i\sim N(\theta,\sigma^2)$ の場合を考える。$\theta=\theta_0+h/\sqrt n$ と置き、$Z_n=\sqrt n(\bar X-\theta_0)$ とすると $Z_n\sim N(h,\sigma^2)$ である。残差は $Z_n$ と独立で、分布は $h$ に依存しない。二乗損失なら、最大リスクを増やさずに残差を平均化できる。実際、残差を条件付きで積分した決定値の二乗リスクは Jensen の不等式で小さくなる。

$h$ に $N(0,\tau^2)$ の事前分布を置けば、$Z_n$ を知った後の分散は平方完成で $(\tau^{-2}+\sigma^{-2})^{-1}$ となる。事後平均は事後二乗損失を最小にするので、
$$
\sup_{h\in\R}E_h(\widehat h-h)^2
\ge E_{\pi_\tau}E_h(\widehat h-h)^2
\ge(\tau^{-2}+\sigma^{-2})^{-1}
$$
である。$\tau\to\infty$ とすると右辺は $\sigma^2$ になり、決定値 $Z_n$ が全 $h$ でこのリスクを達成する。従ってこの正規模型における局所表示では、どの推定量も全方向の最大二乗リスクを $\sigma^2$ 未満にできない。

同じ結論は大きい有限区間の極限でも得られる。$|h|\le H$ に制限した問題では、任意の規則を $[-H,H]$ に射影すると二乗リスクが減る。従って $|h|>H$ での損失は $(|h|+H)^2\le4h^2$ と抑えてよい。上の正規事前について
$$
\begin{aligned}
(\tau^{-2}+\sigma^{-2})^{-1}
&\le\int E_h(\widehat h-h)^2\,\pi_\tau(\dd h)\\
&\le\sup_{|h|\le H}E_h(\widehat h-h)^2+
4\int_{|h|>H}h^2\,\pi_\tau(\dd h).
\end{aligned}
$$
固定 $\tau$ で $H\to\infty$ としてから $\tau\to\infty$ とすれば、有限な局所区間での最小最大リスクも $\sigma^2$ へ近づく。この順序により、分散を大きくした事前の裾が固定した区間からはみ出す問題を避けられる。

:::example[Hodges 推定量の局所リスクを直接比較する]
分散 $1$、$\theta_0=0$ の例で、局所対立 $h/\sqrt n$ を考える。$\sqrt n$ 倍した中心化誤差は、縮小が起きると定数 $-h$、起きないと元の標本平均の誤差 $Z$ である。従って拡大二乗リスクは
$$
h^2 P(|Z+h|\le n^{1/4})
+E[Z^2\mathbf1_{\{|Z+h|>n^{1/4}\}}]\longrightarrow h^2
$$
である。二項目は $Z^2$ で上から抑えられるので支配収束を使える。$h=0$ ではリスクが $0$ に収束し、$h=2$ では $4$ に収束する。標本平均の拡大リスクは常に $1$ であるため、$h=0$ での利得と $h=2$ での損失を同じ尺度で比較できる。
:::

## 局所漸近正規性と Le Cam の第三補題

局所対立は標本数とともに真値へ近づく。標本全体の尤度比を展開すると、その変化が非退化な正規の問題として残る。

:::theorem[局所漸近正規性]{#thm-stat-lan}
:ref[thm-stat-mle-normal] のスコアと Hessian の仮定が成り立つとする。$\Delta_n=n^{-1/2}\sum_i s_{\theta_0}(X_i)$ とすれば、任意の固定 $K$ について $\norm h\le K$ で一様に
$$
\ell_n(\theta_0+h/\sqrt n)-\ell_n(\theta_0)
=h^{\mathsf T}\Delta_n-\tfrac12 h^{\mathsf T}Ih+o_P(1),
\quad \Delta_n\Rightarrow N(0,I).
\label{eq-stat-lan}
$$
この性質を :term[局所漸近正規性]{en="local asymptotic normality" yomi="きょくしょぜんきんせいきせい"} という。
:::

:::proof
二次の積分形 Taylor 展開を用いる。線分上の標本平均 Hessian は $\norm h\le K$ で一様に $-I+o_P(1)$ なので、二次項の剰余は高々 $K^2o_P(1)$ である。一次項の極限は多変量中心極限定理である。
:::

:::theorem[Le Cam の第三補題の尤度比形]{#thm-stat-le-cam-third}
$Q_n\ll P_n$ とし、$L_n=\dd Q_n/\dd P_n$ とする。$P_n$ のもとで $(T_n,L_n)\Rightarrow (T,L)$、$L\ge0$、$EL=1$ なら、$Q_n$ のもとで $T_n$ は $E[f(T)L]$ を期待値に持つ法則へ分布収束する。
:::

:::proof
有界連続 $f$ に対して $E_{Q_n}f(T_n)=E_{P_n}[f(T_n)L_n]$ である。$L_n$ を $M$ で切れば結合収束で極限へ移せる。切った差の期待値は $\norm f_\infty[1-E_{P_n}\min(L_n,M)]$ 以下で、上極限は $\norm f_\infty[1-E\min(L,M)]$ である。$M\to\infty$ で消えるので結論を得る。
:::

:::example[正規モデルの位置移動]
$N(\theta,4)$ で $\theta_0=0$ とすれば $\Delta_n=\sqrt n\bar X/4$、$I=1/4$ である。局所対立 $h=2$ の比は $\exp(2\Delta_n-1/2)$ となる。$\Delta_n$ の $N(0,1/4)$ 密度に掛けて平方完成すると平均は $1/2$ になり、分散は $1/4$ のままである。
:::

:::example[相関による極限の変化]
$(T,\Delta)$ が中心正規で $\operatorname{Var}\Delta=I$、$\operatorname{Cov}(T,\Delta)=C$ とする。重み $e^{h^{\mathsf T}\Delta-h^{\mathsf T}Ih/2}$ を掛けた $T$ の積率母関数は、同時正規の母関数を $(t,h)$ で評価して $e^{t^{\mathsf T}Ch+t^{\mathsf T}\operatorname{Cov}(T)t/2}$ となる。従って極限平均は $Ch$、分散は不変である。$C=0$ ならこのスケールの局所対立に対する極限分布は変わらない。
:::

### 接近する分布の間で確率ゼロの極限を移す

$Q_n$ が $P_n$ に :term[連接]{en="contiguity" yomi="れんせつ"} するとは、どの事象列 $A_n$ についても $P_n(A_n)\to0$ なら $Q_n(A_n)\to0$ となることである。局所対立への尤度比 $L_n$ が $P_n$ の下で $L\ge0$ に分布収束し $EL=1$ なら、連接性が従う。実際、$Q_n(A_n)=E_{P_n}L_n\mathbf1_{A_n}\le M P_n(A_n)+1-E_{P_n}\min(L_n,M)$ である。固定 $M$ で極限を取り、次に $M\to\infty$ とすれば右辺はゼロになる。さらに $L>0$ がほぼ確実なら、逆向きも従う。$P_n(A_n)\le P_n(L_n\le c)+Q_n(A_n)/c$ であり、連続点 $c>0$ をゼロに近づけると第一項も消える。

:::example[最尤展開を局所対立へ移す]
:eqref[eq-stat-lan] の極限比 $e^{h^{\mathsf T}\Delta-h^{\mathsf T}Ih/2}$ は正で平均 $1$ である。従って :eqref[eq-stat-mle-expansion] の $o_P(1)$ は固定 $h$ の局所対立でも確率ゼロへ収束する。一方、第三補題から $\Delta_n$ は対立の下で $N(Ih,I)$ に収束する。よって $\sqrt n(\widehat\theta_n-\theta_0-h/\sqrt n)$ は中心 $0$、共分散 $I^{-1}$ の正規分布へ収束する。分散 $4$ の正規位置で $h=2$ なら、元の真値で中心化した平均の極限平均は $2$、局所の真値で中心化すると $0$ になる。
:::

### 二次平均微分からの局所尤度展開

滑らかな対数尤度の仮定を外し、平方根密度の微分だけから :eqref[eq-stat-lan] の固定方向版を得る。$p_0=p_{\theta_0}$、$p_{n,h}=p_{\theta_0+h/\sqrt n}$ とし、$p_0>0$ 上で $t_n=\sqrt{p_{n,h}/p_0}-1$ とする。二次平均微分の定義から
$$
t_n(X)=\frac{h^{\mathsf T}s(X)}{2\sqrt n}+r_n(X),\qquad
nE_0r_n^2\to0.
$$
$p_0=0$ 上での $p_{n,h}$ の質量も $o(1/n)$ である。これはその領域では平方根密度の定数項と線形項がともにゼロになるためである。正規化すると $2E_0t_n+E_0t_n^2=o(1/n)$ なので、$\sqrt n E_0t_n\to0$ である。一次展開を代入して $h^{\mathsf T}E_0s=0$ を得る。方向 $h$ を任意に選べるため $E_0s=0$ である。

標本上では $\sum_i r_n(X_i)^2\to_p0$ が Markov の不等式で従う。さらに大数の法則と Cauchy–Schwarz で $\sum_i t_n(X_i)^2\to_p h^{\mathsf T}Ih/4$ となる。$\max_i|t_n(X_i)|\to_p0$ も成り立つ。実際、剰余は二乗和で上から抑えられる。線形項については固定 $\epsilon>0$ に対し $nP(|h^{\mathsf T}s|>\epsilon\sqrt n)\le E[(h^{\mathsf T}s)^2\mathbf1_{\{|h^{\mathsf T}s|>\epsilon\sqrt n\}}]/\epsilon^2\to0$ である。従って対数展開の三次剰余の和は、最大値と二乗和の積で確率ゼロに収束する。

中心化した $\sum_i r_n(X_i)$ の分散は $n\operatorname{Var}r_n=o(1)$ である。さらに正規化関係と $nE_0t_n^2\to h^{\mathsf T}Ih/4$ から $2nE_0 t_n\to-h^{\mathsf T}Ih/4$ である。以上をまとめると
$$
2\sum_i\log(1+t_n(X_i))=
h^{\mathsf T}\Delta_n-\tfrac12h^{\mathsf T}Ih+o_P(1).
$$
対立の積法則のうち基準積法則の台の外にある質量は高々 $n$ 倍の単一観測の台外質量なので $o(1)$ である。第三補題の証明ではこの部分を除いても、有界な試験関数の積分の誤差は $o(1)$ である。従って固定方向の局所代替への移送にもこの展開を使える。一方、コンパクト上の一様版が必要なときは、ここで証明した固定方向版と区別する。

### 畳み込み下界の証明

:ref[thm-stat-convolution-bound] の正則推定量について $U_n=\sqrt n(T_n-\theta_0)$、$W_n=I^{-1}\Delta_n$ とする。基準の下で $U_n$ は正則性の極限法則 $L$ に、$W_n$ は $N(0,I^{-1})$ に弱収束する。周辺の緊密性から組も緊密であり、:ref[thm-prob-prokhorov] で結合極限 $(U,W)$ を持つ部分列をとれる。任意の $h\in\R^d$ で、この同じ部分列に固定方向の局所展開と :ref[thm-stat-le-cam-third] を使う。正則性により局所代替での $U_n$ は $L$ を $h$ だけ移動した法則に近づく。従って任意の $t$ について
$$
E[e^{it^{\mathsf T}U}e^{h^{\mathsf T}IW}]
=\varphi_L(t)\exp\!\left(it^{\mathsf T}h+\tfrac12h^{\mathsf T}Ih\right).
$$
ここで $\varphi_L$ は $L$ の特性関数である。$W$ の周辺が正規なので、左辺は $h$ を任意の固定実方向に沿って複素化しても正則である。コンパクト集合上で $e^{c|v^{\mathsf T}IW|}$ とその多項式倍を積分できるためである。一致の定理で、実数 $h$ での等式を $h=i I^{-1}(v-t)$ に移せる。よって
$$
E[e^{it^{\mathsf T}(U-W)+iv^{\mathsf T}W}]
=\varphi_L(t)e^{t^{\mathsf T}I^{-1}t/2}\,
e^{-v^{\mathsf T}I^{-1}v/2}.
$$
$v=0$ と $t=0$ をそれぞれ入れると、右辺の二因子は $U-W$ と $W$ の特性関数である。従って両者は独立で、$U=(U-W)+W$ の法則が主張の畳み込みになる。有限な共分散を持つ場合、独立な和の共分散は足し算なので、余分な成分の共分散は分散を減らせない。

### 二次平均微分の条件を密度から調べる

平方根密度の点ごとの導関数だけで二次平均微分を結論できない。パラメータを少しずらしたときの導関数が、共通の二乗可積分関数で抑えられるかを確認する。一次元で $q_\theta=\sqrt{p_\theta}$ の点ごとの導関数 $\dot q_\theta$ が真値近傍で連続で、$|\dot q_\theta(x)|\le M(x)$、$\int M^2\,\dd\nu<\infty$ とする。このとき
$$
q_{\theta+h}(x)-q_\theta(x)-h\dot q_\theta(x)
=h\int_0^1\{\dot q_{\theta+uh}(x)-\dot q_\theta(x)\}\,\dd u
$$
の右辺の被積分関数は、$h$ で割ると点ごとにゼロへ近づき $2M$ で抑えられる。支配収束から二乗ノルムは $o(h^2)$ である。実数パラメータの内点で $q_\theta(x)=0$ なら微分可能な非負関数の最小点なので $\dot q_\theta(x)=0$ である。従って $2\dot q_\theta/q_\theta$ を正密度上で定義すれば、これが二次平均微分のスコアになる。

:::example[指数分布の平方根密度]
率 $\lambda>0$ について、$x\ge0$ で $q_\lambda(x)=\sqrt\lambda e^{-\lambda x/2}$ である。
$$
\dot q_\lambda(x)=\tfrac12(1/\lambda-x)\sqrt\lambda e^{-\lambda x/2}.
$$
真値 $\lambda_0=2$ の周囲 $1\le\lambda\le3$ では絶対値が $\sqrt3(1+x)e^{-x/2}/2$ 以下なので、この上界は二乗可積分である。情報は $\operatorname{Var}_{2}(X)=1/4$ であり、方向 $h=2$ の局所対立への対数尤度比の極限は、基準の下で平均 $-1/2$、分散 $1$ の正規分布になる。
:::

:::example[端点の移動が平方根密度に残る場合]
$P_\theta=\operatorname{Unif}(0,\theta)$ とする。$h>0$ のとき $[\theta,\theta+h]$ 上では $p_\theta=0$、$p_{\theta+h}=1/(\theta+h)$ である。二次平均微分の線形項は $p_\theta=0$ の領域ではゼロなので、残差ノルムの二乗は少なくとも $h/(\theta+h)$ である。これは $o(h^2)$ ではない。さらに $h=c/n$、$c>0$ として
$$
P_{\theta+c/n}^{\otimes n}(\max X_i>\theta)
=1-\left(\frac{\theta}{\theta+c/n}\right)^n
\longrightarrow1-e^{-c/\theta}.
$$
基準法則ではこの事象の確率はゼロであるため、$1/n$ の方向でさえ対立法則は基準法則に連接しない。局所化の速さを変えても、共通台を前提とする正規移動の計算には戻らない。
:::

### 局所対立での検出力

一次元で :eqref[eq-stat-lan] が成立し $I>0$ とする。基準法則では $Z_n=I^{-1/2}\Delta_n\Rightarrow N(0,1)$ である。局所方向 $h$ のもとでは重みの指数は $h\sqrt I Z_n-h^2I/2$ に近づく。標準正規密度 $\phi(z)$ に極限重みを掛けると
$$
\phi(z)e^{h\sqrt I z-h^2 I/2}
=\frac{1}{\sqrt{2\pi}}\exp[-(z-h\sqrt I)^2/2]
$$
である。よって $Z_n>z_{1-\alpha}$ を棄却域とする検定の極限検出力は $1-\Phi(z_{1-\alpha}-h\sqrt I)$ となる。ここで $\Phi$ は標準正規の分布関数、$z_\gamma$ はその $\gamma$ 分位点である。帰無側の境界の一点は極限法則のもとで零確率なので、分布収束から棄却確率を得られる。

:::example[Bernoulli の局所対立と固定対立の違い]
$p_0=1/2$ で $I=4$、$Z_n=2\sqrt n(\bar X-1/2)$ である。$p_n=1/2+1/(2\sqrt n)$ なら $h=1/2$ であり、$Z_n$ の極限は $N(1,1)$ となる。棄却閾値を $1$ にすれば、基準での極限棄却確率は $1-\Phi(1)$、この局所対立では $1/2$ である。一方、固定した $p>1/2$ では $Z_n$ が確率的に正の無限大へ発散する。実際 $Z_n=2\sqrt n(p-1/2)+2\sqrt n(\bar X-p)$ で、後半は確率的に有界である。この場合の極限検出力は $1$ となる。
:::

### 局所尤度比の質量を確認する

:ref[thm-stat-le-cam-third] の $E L=1$ は、極限法則の重みが全確率を残すための条件である。例えば基準法則 $P_n$ の下で $A_n$ の確率を $1/n$、対立法則 $Q_n$ の下でその確率を $1/2$ とする。各事象内の条件付き法則は同じとし、$n\ge2$ とする。尤度比は $A_n$ で $n/2$、補集合で $1/[2(1-1/n)]$ である。基準法則で $L_n\Rightarrow1/2$ だが、この極限の平均は $1/2$ である。重みだけで移送すると半分の質量が失われ、$P_n(A_n)\to0$ でも $Q_n(A_n)=1/2$ だから連接性は成立しない。対照的に正規の局所尤度比の極限では、指数の定数項が正規の積率母関数をちょうど打ち消す。

### 複数の方向で同じ統計量を使う場合

$h_1,\ldots,h_k$ を固定すると、それぞれの対数尤度比の極限は同じ $\Delta$ に $h_j^{\mathsf T}$ を掛けて得る。有限個の剰余の最大は、和事象の上界から $o_P(1)$ のままである。従ってこれらの尤度比をまとめたベクトルは、正規移動モデルの尤度比ベクトルへ収束する。局所方向ごとに独立な $Z$ を用意すると、方向間の相関 $h_i^{\mathsf T}Ih_j$ を失う。

:::example[逆向きの二つの局所対立]
$I=2$、方向を $h_+=1,h_-=-1$ とする。基準法則の極限では対数尤度比は $\Delta-1,-\Delta-1$ なので和が $-2$、共分散は $-2$ である。両対立を同じ重みで混ぜた法則の尤度比は $e^{-1}\cosh\Delta$ となる。従って尤度比の大きい順に棄却するなら $|\Delta|$ の大きい両裾を使う。片側で検出力を求めるときと異なり、正負の対立に対して同じ確率で棄却する。
:::

### 局所漸近最小最大下界

二乗損失では、正則性を仮定せずに局所リスクの下界を示せる。真値が内点で、正定値情報行列 $I$ を持つ固定方向の局所展開と、基準法則の台の外の質量が $o(1)$ となる条件を仮定する。これは上で示した二次平均微分可能な独立同分布モデルを含む。任意の推定量列 $T_n$ に対して
$$
\lim_{H\to\infty}\liminf_{n\to\infty}
\sup_{\norm h\le H}E_{\theta_0+h/\sqrt n}
\bigl\|\sqrt n(T_n-\theta_0)-h\bigr\|^2
\ \ge\ \tr(I^{-1})
$$
となる。まず有限個の局所方向で平均したリスクを下から抑え、次いでそれらの方向を増やして正規位置模型のリスクに近づける。

半径 $H$ の球内の点 $h_j$ に確率 $w_j>0$ を置く。球に射影された行動 $a$ だけを比べればよい。球内のどの $h_j$ に対しても射影は二乗距離を増やさないためである。基準法則のもとで各点の尤度比を $L_{nj}$ とすると、平均リスクは
$$
E_0\, b(L_{n1},\ldots,L_{nk}),\qquad
b(l)=\min_{\norm a\le H}\sum_{j=1}^k w_j l_j\norm{a-h_j}^2
$$
以上である。基準の台外の寄与は損失が非負なので、この下界では落としてよい。$b$ は連続で $0\le b(l)\le 4H^2\sum_jw_jl_j$ を満たす。各 $L_{nj}$ の期待値は $1$ に近づき、結合極限の各成分 $L_j=\exp(h_j^{\mathsf T}\Delta-h_j^{\mathsf T}Ih_j/2)$ も期待値 $1$ を持つ。第三補題の証明と同じ切断により $\sum_j L_{nj}$ は一様可積分である。従ってこの下界の極限は $E b(L_1,\ldots,L_k)$ である。これは $Z\sim N(h,I^{-1})$ を観測する模型の有限事前に対する最小平均リスクである。実際、基準の観測は $Z=I^{-1}\Delta$ であり、点ごとの最適行動は $h_j$ の事後重み付き平均なので球内に入る。

この正規位置模型で $h\sim N(0,\tau^2 I_d)$ としたとき、平方完成による事後共分散は $(\tau^{-2}I_d+I)^{-1}$ であり、最小平均リスクはその跡である。半径 $H$ の外の事前を除く影響は、球内の行動に対して $4E[\norm h^2\mathbf1_{\{\norm h>H\}}]$ 以下である。従って $\tau$ を固定して $H$ を大きく取れば、球に条件付けた事前の最小平均リスクも元の値に任意に近い下界を持つ。さらに球を有限の網で分け、各点を距離 $\delta$ 以内の網の一点に置き換える。球内の行動では二乗損失の変化は高々 $4H\delta$ であり、正規密度の平行移動の $L^1$ 連続性により観測法則を入れ替えた誤差も $\delta\downarrow0$ で一様に消える。このため有限事前でも同じ下界に任意に近づける。最後に $\tau\to\infty$ とすれば所望の跡を得る。ここで比較したのは各局所方向での極限分散ではなく、局所領域全体の二乗リスクである。

## 演習問題

::::exercise{#exr-stat-asymptotics-1 level=basic}
$\sqrt n(T_n-2)\Rightarrow N(0,9)$ のとき $\sqrt n(T_n^2-4)$ の極限を求めよ。

:::solution
$x^2$ の $2$ での傾きは $4$ なので、極限は分散 $4^2\cdot9=144$ の中心正規である。
:::
::::

::::exercise{#exr-stat-asymptotics-2 level=basic}
$A_n\to_p0$ と $B_n=O_P(1)$ から $A_nB_n\to_p0$ を示せ。

:::solution
$P(|A_nB_n|>\epsilon)\le P(|B_n|>M)+P(|A_n|>\epsilon/M)$ とする。前者は大きい $M$ で任意に小さくでき、固定 $M$ で後者はゼロに近づく。
:::
::::

::::exercise{#exr-stat-asymptotics-3 level=basic}
独立な $X_i\sim\operatorname{Pois}(9)$ で $\sqrt n(\sqrt{\bar X}-3)$ の極限分散を求めよ。

:::solution
微分は $1/(2\sqrt9)=1/6$ で平均の漸近分散は $9$ なので $9/36=1/4$ である。
:::
::::

::::exercise{#exr-stat-asymptotics-4 level=basic}
$A=-3,B=18$ の Z 推定量の漸近分散を求めよ。

:::solution
一次元の式 $B/A^2$ に代入して $18/9=2$ となる。
:::
::::

::::exercise{#exr-stat-asymptotics-5 level=advanced}
$M(\theta)=-\theta^2$、$\sup|M_n-M|\le n^{-1/2}$ とする。最大点の誤差を評価せよ。

:::solution
$M_n(\widehat\theta)\ge M_n(0)$ と一様誤差を使うと $\widehat\theta^2\le2n^{-1/2}$ になる。従って $|\widehat\theta|\le\sqrt2 n^{-1/4}$ である。
:::
::::

::::exercise{#exr-stat-asymptotics-6 level=advanced}
一様端点模型で $n(1-M/\theta)$ の極限を求めよ。

:::solution
固定 $t\ge0$ で右裾は $P(M/\theta<1-t/n)=(1-t/n)^n$（$n>t$）である。極限は $e^{-t}$ なので率 $1$ の指数分布である。
:::
::::

::::exercise{#exr-stat-asymptotics-7 level=advanced}
Hodges 推定量の閾値を $n^{-a}$ とする。$0<a<1/2$ で局所対立における極限を求めよ。

:::solution
$\sqrt n\bar X=Z+h$、閾値は $\sqrt n$ 倍すると $n^{1/2-a}\to\infty$ である。従って $T_n=0$ となる確率が $1$ に近づき、$\sqrt n(T_n-\theta_n)\to_p-h$ である。
:::
::::

::::exercise{#exr-stat-asymptotics-8 level=advanced}
中心二変量正規で $\operatorname{Var}\Delta=2$、$\operatorname{Cov}(T,\Delta)=3$、$\operatorname{Var}T=5$ とする。$h=2$ の尤度比重みでの $T$ の極限を求めよ。

:::solution
母関数の指数は $(5t^2+2\cdot3\cdot2t+2\cdot2^2)/2-2\cdot2^2/2=5t^2/2+6t$ となる。従って $N(6,5)$ である。
:::
::::

## 確認テスト

::::quiz
Slutsky の補題が直接扱う組を選べ。

- [x] 一方が分布収束、他方が定数へ確率収束する。
- [ ] 二つの周辺分布収束だけから同時分布収束が得られる。
- [ ] ゼロへ確率収束する分母で割っても必ず有限極限になる。

:::explanation
定数への確率収束があるので、同時分布でその成分を定数に置き換えられる。商には分母の極限定数が非零であることが要る。
:::
::::

::::quiz{answer="144" tol="1e-6"}
$\sqrt n(T_n-2)$ の極限分散が $9$ の場合、二乗に対するデルタ法の極限分散を求めよ。

:::explanation
二乗の導関数は $2t$ なので $t=2$ で $4$ になる。分散には導関数の二乗を掛けて $144$ を得る。
:::
::::

::::quiz
Hodges 推定量が正則性を失う理由を選べ。

- [ ] どの真値でも一致しない。
- [x] 局所代替下での中心化誤差の法則が方向 $h$ に依存する。
- [ ] 推定量が観測から計算できない。

:::explanation
局所代替 $h/\sqrt n$ でも閾値内に入りやすい。推定量はゼロに設定されるので真値で中心化した誤差は $-h$ に近づく。
:::
::::

::::quiz{answer="2" tol="1e-6"}
一次元 Z 推定量で $A=-3,B=18$ の場合の漸近分散を答えよ。

:::explanation
推定式の標本変動 $B$ に、傾きの逆の二乗 $1/A^2$ を掛けて $2$ となる。
:::
::::

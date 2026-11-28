---
title: "無限分解可能分布と安定分布"
description: "畳み込みの根、Lévy 測度、正規分布以外の和の極限を扱う。"
notation:
  - symbol: '\mu^{*n}'
    meaning: "測度 μ の n 重畳み込み"
  - symbol: '\nu'
    meaning: "Lévy 測度"
  - symbol: '\psi'
    meaning: "連続な特性指数"
---

## この章の目標

- 無限分解可能性と複合 Poisson 分布を具体例から判定する。
- Lévy–Khintchine 表現を証明し、三角配列の候補極限を求める。
- 安定分布を Lévy 測度の拡大則で特徴付ける。
- 裾のべき減衰から和の正規化を選ぶ。

## 無限分解可能分布の定義と例

同じ分布の独立な項を足して一つの観測値を作る場合、項数を何通りにも選べることがある。項の分布も項数とともに変えてよいとすると、正規分布以外にもこの性質をもつ分布が得られる。

:::definition[無限分解可能性]{#def-prob-infinitely-divisible}
確率測度 $\mu$ が :term[無限分解可能]{en="infinitely divisible" yomi="むげんぶんかいかのう"} であるとは、任意の整数 $n\ge1$ に対し、ある確率測度 $\mu_n$ が $\mu=\mu_n^{*n}$ を満たすことをいう。
:::

正規分布 $N(m,v)$ の根は $N(m/n,v/n)$ である。Poisson 分布 $\operatorname{Pois}(\lambda)$ の根は $\operatorname{Pois}(\lambda/n)$ である。いずれも特性関数を $n$ 乗して確かめられる。例えば $N(6,9)$ の三乗根は $N(2,3)$ である。また無限分解可能な分布の畳み込みでは、二つの $n$ 乗根を畳み込めばよい。

:::example[有界な分布の制約]
$\mu$ の台が長さ $L<\infty$ の区間に入るとする。根の本質的上限・下限が有限であることは、一項が任意に大きい値をとり得るとき残りを固定した有界区間に入れると和も任意に大きくなることから従う。さらにすべての根の変数が上限（下限）のいくらでも近くに来る確率は独立性より正なので、根の台の幅は $L/n$ 以下である。従って $\operatorname{Var}(\mu)=n\operatorname{Var}(\mu_n)\le L^2/n$ である。$n\to\infty$ とすると点質量に限られる。例えば公平な Bernoulli 分布は無限分解可能でない。
:::

### 根を揃えると時間を添字にできる

:::proposition[無限分解可能分布の畳み込み半群]{#prop-prob-convolution-semigroup}
無限分解可能な $\mu$ について $\mu_0=\delta_0$、$\mu_1=\mu$、$\mu_s*\mu_t=\mu_{s+t}$ を満たす弱連続な確率測度族 $(\mu_t)_{t\ge0}$ が一意に存在する。
:::

:::proof
:ref[thm-prob-levy-khintchine] が与える指数を $\psi$ とする。全 $t\ge0$ で、三つ組を $t$ 倍した特性関数 $e^{t\psi}$ は確率測度を定める。積公式により半群性が成立し、$t_n\to t$ なら特性関数が点ごとに $e^{t\psi}$ へ向かうので弱連続である。

一意性には、まず任意の正整数 $n$ に対して $\phi_{\mu_{1/n}}^n=\phi_\mu$ を使う。$\phi_\mu$ は零を持たない。商 $\phi_{\mu_{1/n}}/e^{\psi/n}$ は連続で $n$ 乗が $1$ であり、実直線が連結で原点の値が $1$ なので恒等的に $1$ である。半群性で非負有理数時刻の法則が決まり、弱連続性で残りの時刻も決まる。
:::

:::example[率と時間を分けて管理する]
$\mu=\operatorname{Pois}(6)$ なら $\mu_{1/3}=\operatorname{Pois}(2)$、$\mu_{2/3}=\operatorname{Pois}(4)$ である。前半の区間でゼロ回、後半で二回となる確率は、それらを独立増分として実現すれば $e^{-2}(4^2e^{-4}/2)=8e^{-6}$ となる。このように時間を分けたとき、同じ分解を整数時刻と分数時刻で使える。
:::

:::example[Gamma 分布の根]
形状 $\alpha$、率 $\lambda$ の Gamma 分布の密度を積分すると Laplace 変換は $(1+s/\lambda)^{-\alpha}$（$s\ge0$）である。形状が $\alpha/n$ の独立な $n$ 個の和では変換が元と等しい。あるいは :ref[exm-prob-gamma-convolution] を繰り返せば密度の畳み込みで同じ分布を得る。例えば形状 $3$、率 $2$ の分布は、三つの率 $2$ の指数変数の和であり、二乗根は形状 $3/2$、率 $2$ である。
:::

### 弱極限を取っても分解可能性は残る

:::proposition[無限分解可能分布の弱閉性]{#prop-prob-id-weak-closure}
$\mu_j$ が無限分解可能で $\mu_j\Rightarrow\mu$ ならば、$\mu$ も無限分解可能である。
:::

:::proof
$f_j,f$ を特性関数とする。緊密性により $f_j\to f$ はコンパクト区間上一様である。まず $f$ が零を持たないことを示す。$f_j$ の連続対数を $\psi_j$ とすると、$n$ 乗根の特性関数の絶対値の二乗は $\exp(2\Re\psi_j/n)$ である。:ref[thm-prob-levy-khintchine] の証明にある余弦の不等式をこの根に使い、$n$ 倍して $n\to\infty$ とすれば
$$
-\log |f_j(t)|\le m^2(-\log |f_j(t/m)|)
$$
になる。固定した $t$ に対して $m$ を十分大きく取ると、$f(t/m)$ は零でない。従って右辺は $j$ について有界であり、$|f(t)|>0$ が従う。

$f$ の連続対数 $\psi$ を原点で零となるよう取る。任意の原点を含むコンパクト区間で $f$ は零から離れるため、$f_j/f$ は一様に $1$ へ向かう。この比の $1$ の近くの対数は原点で零であり、連続対数の一意性から $\psi_j-\psi$ に等しい。従って固定 $n$ について $e^{\psi_j/n}\to e^{\psi/n}$ となる。極限は原点で連続なので Lévy の連続性定理を適用すると、$\mu$ の $n$ 乗根が得られる。
:::

:::example[有限跳躍だけの分布から正規成分が出る]
$X_j=(N_j^+-N_j^-)/\sqrt j$ とし、$N_j^\pm$ は互いに独立な $\operatorname{Pois}(j/2)$ とする。指数は $j(\cos(t/\sqrt j)-1)$ であり、$-t^2/2$ へ収束する。各 $j$ で正規成分は零、跳躍の強度は両符号合わせて $j$ である。しかし幅 $1/\sqrt j$ の二乗を掛けた全強度は $1$ を保つため、極限では正規分散 $1$ に移る。跳躍測度を零から離れた集合でだけ比較すると、この分散の寄与を失う。
:::

## 複合 Poisson 分布

:::definition[複合 Poisson 分布]{#def-prob-compound-poisson}
$N\sim\operatorname{Pois}(\lambda)$ と独立同分布の跳躍 $J_k\sim Q$ が独立なとき、$S=\sum_{k=1}^NJ_k$（空和は零）の分布を :term[複合 Poisson 分布]{en="compound Poisson distribution" yomi="ふくごうぽあそんぶんぷ"} という。
:::

$N=n$ の場合の特性関数は $\phi_Q(t)^n$ なので、Poisson の質量で平均して $\phi_S(t)=\exp(\lambda(\phi_Q(t)-1))$ となる。強度を $\lambda/n$ に変えると根が得られる。$J$ に二次積率があれば、$n$ 個の和の平均 $n\mathbb EJ$ と分散 $n\operatorname{Var}(J)$ を使って $\mathbb ES=\lambda\mathbb EJ$、$\operatorname{Var}(S)=\lambda\mathbb EJ^2$ を得る。$N$ がゆらぐため、第 2 式には跳躍の中心積率と平均の二乗がともに入る。

:::example[二方向の跳躍]
$\lambda=3$ で $J$ が $-1,2$ を等確率で取ると、指数は $\tfrac32(e^{-it}+e^{2it}-2)$ である。従って分布は独立な $\operatorname{Pois}(3/2)$ の $N_1,N_2$ を使って $-N_1+2N_2$ とも表せる。平均は $3/2$、分散は $15/2$ である。両向きの跳躍が打ち消し合うことがあるので、和が零の確率には $N=0$ 以外の寄与もある。
:::

### 跳躍の種類ごとの Poisson 回数

:::proposition[複合 Poisson の分割]{#prop-prob-compound-thinning}
跳躍の法則が有限混合 $Q=\sum_j p_jQ_j$ なら、その複合 Poisson 和の分布は強度 $\lambda p_j$、跳躍法則 $Q_j$ の独立な複合 Poisson 和を足した分布である。$p_j=0$ の項は零とする。
:::

:::proof
候補の独立和の指数は $\sum_j\lambda p_j(\phi_{Q_j}(t)-1)$ である。$\sum_jp_j=1$ と特性関数の線形性により、これは元の指数 $\lambda(\phi_Q(t)-1)$ になる。一意性定理により分布は一致する。
:::

:::example[零の跳躍を数えない]
強度 $5$ で跳躍が確率 $2/5$ で $0$、確率 $3/5$ で $2$ とする。非零ジャンプの強度は $3$ なので、和は $2N$、$N\sim\operatorname{Pois}(3)$ と同分布である。和が $0$ の確率は $e^{-3}$ で、元の回数が $0$ の確率 $e^{-5}$ より大きい。平均は $6$、分散は $12$ であり、跳躍に零が入っていても法則は非零部分の強度で決まる。
:::

:::example[非負跳躍の和の小さな点を数える]
強度 $4$ で跳躍が $1,2$ をそれぞれ確率 $3/4,1/4$ で取る。$S=N_1+2N_2$、$N_1\sim\operatorname{Pois}(3)$、$N_2\sim\operatorname{Pois}(1)$ と表せる。$S=2$ は $(N_1,N_2)=(2,0),(0,1)$ の二通りなので、確率は $(9/2+1)e^{-4}=(11/2)e^{-4}$ である。
:::

### 跳躍の型に合わせて母関数と密度を使い分ける

:::proposition[整数値跳躍の質量の漸化式]
$N$ は率 $\lambda$ の Poisson 変数とし、独立な跳躍は正整数 $j$ を確率 $q_j$ で取るとする。和の質量を $p_m=\mathbb P(S=m)$ とすると
$$
p_0=e^{-\lambda},\qquad
m p_m=\lambda\sum_{j=1}^m j q_j p_{m-j}\quad(m\ge1)
$$
である。
:::

:::proof
$S=0$ となるのは $N=0$ だけである。$|z|<1$ なら $Q(z)=\sum_{j\ge1}q_jz^j$ は絶対収束し、$N$ で条件付けると和の母関数は $P(z)=\exp(\lambda(Q(z)-1))$ になる。コンパクトな部分円板上で各べき級数を項別微分でき、$P'(z)=\lambda Q'(z)P(z)$ である。両辺で $z^{m-1}$ の係数を比較すれば式が出る。和は $j\le m$ だけを使うため、大きな跳躍の平均が無限であっても有限段階ごとの漸化式は使える。
:::

:::example[漸化式で二種類の跳躍を扱う]
$\lambda=3$、$q_1=2/3,q_2=1/3$ とする。$p_0=e^{-3}$、$p_1=2e^{-3}$、$2p_2=2p_1+2p_0$ なので $p_2=3e^{-3}$ である。次に $3p_3=2p_2+2p_1$ なので $p_3=(10/3)e^{-3}$ である。独立な率 $2,1$ の二つの個数 $N_1,N_2$ を使うと $S=N_1+2N_2$ になり、$S=3$ の組は $(3,0),(1,1)$ である。その確率の和 $(2^3/3!+2)e^{-3}$ と一致する。
:::

:::example[連続な跳躍を足しても零に原子が残る]
強度 $\lambda=1$、跳躍は率 $2$ の指数分布とする。$S$ には零の質量 $e^{-1}$ があり、条件 $N=n\ge1$ の下では形状 $n$、率 $2$ の Gamma 密度である。従って $s>0$ の連続成分の密度は
$$
e^{-1-2s}\sum_{n\ge1}\frac{2^n s^{n-1}}{n!(n-1)!}
$$
である。各項の積分は $e^{-1}/n!$ なので、和の全質量は $1-e^{-1}$ である。平均は $\lambda\mathbb EJ=1/2$、分散は $\lambda\mathbb EJ^2=1/2$ である。密度を積分しても $1$ にならない分は、跳躍が一つもない事象に対応する。
:::

## Lévy–Khintchine の表現定理

小さな跳躍が無限個ある場合、跳躍の絶対値の積分は発散し得る。一次項を引いて二次の余りだけを積分すると、こうした極限も表現できる。

:::theorem[Lévy–Khintchine の表現]{#thm-prob-levy-khintchine}
$\mu$ が無限分解可能であることと、その特性関数が次の形をもつことは同値である。
$$
\phi_\mu(t)=e^{\psi(t)},\quad
\psi(t)=it\gamma-\frac{v t^2}{2}
 +\int_{\R\setminus\{0\}}\left(e^{itx}-1-it\frac{x}{1+x^2}\right)\nu(\dd x).
$$
ここで $\gamma\in\R$、$v\ge0$、$\nu$ は $\int (1\wedge x^2)\,\nu(\dd x)<\infty$ を満たす正測度である。三つ組 $(\gamma,v,\nu)$ は一意であり、$\nu$ を :term[Lévy 測度]{en="Lévy measure" yomi="れびそくど"} という。
:::

:::proof
まず表示があれば、$\nu$ を $|x|>\eps$ に制限して積分を有限測度の複合 Poisson 指数と定数移動に分ける。正規因子と掛ければ特性関数である。$\eps\downarrow0$ で積分は収束する。実際、小さい $x$ では被積分関数は $O_t(x^2)$、$|x|\ge1$ では有界であり、$t$ が有界な範囲でこれらの評価は一様である。従って極限は原点で連続であり、:ref[thm-prob-levy-continuity] が使える。三つ組を $1/n$ 倍しても同じ構成が可能であり、根を得る。

逆に根の特性関数を $\phi_n$ とする。$\phi$ は原点近くで零を持たない。整数 $m\ge1$ について $1-\cos(mu)\le m^2(1-\cos u)$、および $|\phi_n(t)|^2$ は独立な差の特性関数であることから
$$
1-|\phi_n(t)|^2\le m^2(1-|\phi_n(t/m)|^2).
$$
固定 $t$ で $t/m$ を原点近くに取れば右辺は $m^2(1-|\phi(t/m)|^{2/n})\to0$ であり、$\phi(t)=0$ を排除できる。従って連続な偏角を選んで原点で零となる連続対数 $\psi$ を取れる。連続性により $n\log\phi_n=\psi$ である。従って $n(\phi_n-1)\to\psi$ は有界 $t$ 区間で一様に成立する。

$w(x)=x^2/(1+x^2)$、$b(x)=x/(1+x^2)$ とし、$\rho_n(\dd x)=nw(x)\mu_n(\dd x)$ と置く。$1-\sin x/x\ge c w(x)$ となる $c>0$ がある。これは原点での比の正の極限と、原点以外での厳密な正値性と、無限遠での正の極限から分かる。$-\Re n(\phi_n(t)-1)$ を $[-1,1]$ で平均すると、この評価により $\sup_n\rho_n(\R)<\infty$ となる。同じ平均を $[-u,u]$ で取ると $|x|\ge2/u$ で $1-\sin(ux)/(ux)\ge1/2$ である。従って $\limsup_n \rho_n(|x|\ge2/u)$ は $2\sup_{|t|\le u}|\Re\psi(t)|$ 以下となり、$u\downarrow0$ で零になる。初期有限項を合わせれば緊密である。有限測度の全質量が収束する部分列を取り、極限が正なら正規化して :ref[thm-prob-helly-selection] を使う。極限が零なら有限測度として零に収束する。こうして弱収束部分列 $\rho_n\Rightarrow\rho$ を得る。

$g_t(x)=(e^{itx}-1-itb(x))/w(x)$ は $x=0$ で値 $-t^2/2$ を補うと有界連続である。$n(\phi_n(t)-1)=it\gamma_n+\int g_t\,\dd\rho_n$、$\gamma_n=n\int b\,\dd\mu_n$ と書ける。$t=1$ の虚部を比較すれば同じ部分列で $\gamma_n\to\gamma$ となる。$v=\rho(\{0\})$、$\nu=w^{-1}\rho|_{\R\setminus\{0\}}$ と置けば表示を得る。

一意性には $\psi(t+h)+\psi(t-h)-2\psi(t)$ を使う。これは有限正測度 $vh^2\delta_0+2(1-\cos(hx))\nu(\dd x)$ の Fourier 変換の負である。有限測度は全質量で正規化すると :ref[thm-prob-fourier-inversion] で一意なので、$h\ne0$ から $v$ が分かる。$h=1,\,\sqrt2$ の二つを使えば $x\ne0$ で両方の重みが零にはならず、合算した重みで割ると $\nu$ が分かる。最後に一次項から $\gamma$ が定まる。
:::

跳躍が $2$ に固定され率が $3$ なら $\nu=3\delta_2$ である。ここで切断関数が $b$ なので、補正項 $it\,3b(2)$ を戻すため $\gamma=6/5$ が必要になる。$\gamma$ は切断関数の選び方に依存し、$v$ と $\nu$ は変わらない。

### 切断関数の変更と有限二次積率

:::proposition[切断の変更と平均・分散]{#prop-prob-levy-cutoff-moments}
:ref[thm-prob-levy-khintchine] の $b(x)=x/(1+x^2)$ を $\tilde b(x)=x\mathbf1_{\{|x|\le1\}}$ に変えると、ドリフトは $\tilde\gamma=\gamma+\int(\tilde b-b)\,\dd\nu$ となる。また、法則が有限二次積率をもつことと $\int x^2\,\nu(\dd x)<\infty$ は同値である。この条件の下で
$$
\mathbb EX=\gamma+\int(x-b(x))\,\nu(\dd x),\qquad
\operatorname{Var}(X)=v+\int x^2\,\nu(\dd x)
$$
である。
:::

:::proof
零付近では $\tilde b-b=O(|x|^3)$、零から離れた場所で有界なので、その積分は絶対収束する。新しい表示で一次項を比較すれば $\tilde\gamma$ の式が得られる。

二次積率の条件の下では $x-b(x)$ は零付近で $O(|x|^3)$、$|x|>1$ で $O(|x|)$ であり、全体で可積分である。$\mu_\eps$ を $|x|>\eps$ に切断した表現の法則とする。これは正規、定数移動、および複合 Poisson の独立和なので、平均は $\gamma+\int_{|x|>\eps}(x-b(x))\,\nu(\dd x)$、分散は $v+\int_{|x|>\eps}x^2\,\nu(\dd x)$ である。これらの極限を法則に移すため、切断を $\eps=2^{-j}$ に限る。各環状部分で中心化した複合 Poisson 和を独立に取ると、差の二次積率は追加した部分の $\int x^2\dd\nu$ である。従って累積和は $L^2$ Cauchy であり、正規項と収束する平均を加えれば全体が $L^2$ で収束する。その弱極限は切断した特性関数の極限で元の $\mu$ なので、平均と分散の式が従う。

逆に $X$ が二乗可積分なら、独立な複製との差の特性関数は $|\phi_\mu(t)|^2$ だから $|\phi_\mu(t)|^2=1-\operatorname{Var}(X)t^2+o(t^2)$ である。従って $-\Re\psi(t)/t^2\to\operatorname{Var}(X)/2$ である。一方 $-\Re\psi(t)=vt^2/2+\int(1-\cos(tx))\,\nu(\dd x)$ である。非負な被積分関数に Fatou を適用すると $\int x^2\,\nu(\dd x)\le\operatorname{Var}(X)-v<\infty$ を得る。
:::

:::example[切断の違いを数値で確認する]
$X=2N$、$N\sim\operatorname{Pois}(3)$ とする。滑らかな切断での $\gamma$ は $6/5$ である。$\tilde b(2)=0$ なので $\tilde\gamma=6/5-3\cdot2/5=0$ となる。一方、どちらの表示でも実際の平均は $6$、分散は $12$ である。ドリフトの数値だけで平均を読めるのは、残った跳躍積分の一次部分が零の場合に限る。
:::

:::example[Gamma の Lévy 測度は無限だが大きな跳躍は有限]
$\nu(\dd x)=\alpha e^{-\lambda x}\mathbf1_{\{x>0\}}\,\dd x/x$ とする。零の近くでは $\int\nu=\infty$ だが、$\int_0^1x\,\nu(\dd x)<\infty$ である。従って補正をまとめて、正の半直線上で Laplace 指数を $\int_0^\infty(e^{-sx}-1)\nu(\dd x)$ とできる。$s$ で微分すると $-\alpha/(\lambda+s)$ であり、$s=0$ の値が $0$ だから指数は $-\alpha\log(1+s/\lambda)$ である。切断した正の複合 Poisson 和は、失われる小跳躍の期待値 $\int_0^\eps x\,\nu(\dd x)\to0$ によって $L^1$ で構成できる。得られる変換は Gamma 分布のものと一致する。二つの非負変数の Laplace 変換が一致すれば、$e^{-X}$ の $[0,1]$ 上での全積率が一致するため、積率決定性で元の法則も一致する。
:::

### 小跳躍の級数として特性指数を実現する

:::proposition[環状領域から補償付き小跳躍を足す]
$b(x)=x/(1+x^2)$ とし、Lévy 測度 $\nu$ のみを考える。$A_0=\{|x|>1\}$、$A_j=\{2^{-j}<|x|\le2^{-j+1}\}$（$j\ge1$）とする。各 $A_j$ で $\nu|_{A_j}$ に対応する複合 Poisson 和 $J_j$ を独立に選ぶ。このとき
$$
J_0-\int_{A_0}b\,\dd\nu+
\sum_{j\ge1}\left(J_j-\int_{A_j}b\,\dd\nu\right)
$$
は概確実に収束し、表示定理の跳躍積分を特性指数に持つ。小跳躍の級数は $L^2$ でも収束する。
:::

:::proof
各 $A_j$ の強度は有限である。$A_0$ では $b$ が有界なので、補正の積分は有限であり、跳躍の個数も概確実に有限である。大跳躍の平均や分散は仮定しない。$j\ge1$ では跳躍和は二乗可積分で、平均は $\int_{A_j}x\,\dd\nu$、分散は $\int_{A_j}x^2\,\dd\nu$ である。従って $J_j-\mathbb E J_j$ の分散の和は $\int_{0<|x|\le1}x^2\,\dd\nu<\infty$ である。Kolmogorov の最大不等式で末尾の Cauchy 条件を調べる :ref[thm-prob-three-series] の議論により概収束し、直交性から $L^2$ Cauchy でもある。残りの定数の和は $\int_{0<|x|\le1}(x-b(x))\,\dd\nu$ であり、$|x-b(x)|\le|x|^3$ だから絶対収束する。有限部分和の特性指数を足し、$|e^{itx}-1-itb(x)|\le C_t x^2$（$|x|\le1$）で積分の極限を取る。概収束した確率変数の特性関数は有界収束で同じ極限となる。
:::

:::example[各環で跳躍の絶対値は大きく足し上がる]
$\nu$ が正の $2^{-j}$（$j\ge1$）に質量 $2^j$ を置くとする。跳躍の平均絶対値の和は $\sum_j1=\infty$ だが、分散の和は $\sum_j2^{-j}=1$ である。$N_j\sim\operatorname{Pois}(2^j)$ とすると、補償付き和 $\sum_{j\ge1}2^{-j}(N_j-2^j)$ は概確実かつ $L^2$ で収束し、平均零、分散 $1$ を持つ。これに対して $\sum 2^{-j}N_j$ の各部分和は補償付き部分和に項数を足した値なので、概確実に無限大へ向かう。この例では一次補償が本当に必要である。
:::

### 測度の台の条件から片側に値を取る和を作る

:::proposition[正の有限変動の跳躍と Laplace 指数]
Lévy 測度 $\nu$ が $(0,\infty)$ に集中し、$\int_0^1 x\nu(\dd x)<\infty$ とする。$d\ge0$ に対し、非負確率変数 $Y$ で
$$
-\log\mathbb E e^{-sY}
=ds+\int_{(0,\infty)}(1-e^{-sx})\,\nu(\dd x),\qquad s\ge0
$$
を満たすものが存在する。その法則は無限分解可能である。
:::

:::proof
零から離れた環ごとに独立な正の複合 Poisson 和を作る。$x>1$ の強度は有限だから、その部分は有限個の有限な値の和である。$0<x\le1$ の全ての環の和は非負で、期待値を足すと $\int_0^1x\,\nu(\dd x)$ である。従って和は概確実に有限となる。$d$ を足した和を $Y$ とすれば、有限個の環では積公式によって指数は右辺の制限された積分になる。$1-e^{-sx}\le sx$ と大跳躍の有限強度により全体で積分でき、環を増やすと単調収束で右辺を得る。

さらに全ての環の強度と $d$ を $1/n$ にした同じ構成が可能である。独立な $n$ コピーの和は元の $Y$ と同じ Laplace 変換を持つ。非負変数の Laplace 変換が等しければ、指数を取った $e^{-Y}$ の全積率が等しく、台が有界なので分布が一致する。従って $n$ 乗根が得られる。
:::

:::example[離散的な小跳躍の和]
$\nu=\sum_{j\ge1}j\,\delta_{2^{-j}}$ とする。全強度は無限だが、小跳躍の平均は $\sum_{j\ge1}j2^{-j}=2$ で有限である。独立な $N_j\sim\operatorname{Pois}(j)$ なら $Y=\sum_{j\ge1}2^{-j}N_j$ は概確実に有限である。分散の和は $\sum_{j\ge1}j4^{-j}=4/9$。この式は有限部分和の直交性と $L^2$ Cauchy 性で極限へ渡せる。二つの和は $\sum_{j\ge1}jr^j=r/(1-r)^2$ に $r=1/2,1/4$ を入れて求める。この例は全強度が無限であることと、小跳躍に補償が必要であることが別の条件であることを示す。
:::

:::proposition[非負な無限分解可能分布の特徴付け]{#prop-prob-nonnegative-id}
無限分解可能な法則が $[0,\infty)$ に集中するための必要十分条件は、$v=0$、$\nu((-\infty,0))=0$、$\int_{(0,1]}x\,\nu(\dd x)<\infty$、かつ $d=\gamma-\int b\,\dd\nu\ge0$ である。ここで $b(x)=x/(1+x^2)$ である。その Laplace 変換は上の正の跳躍の公式で与えられる。
:::

:::proof
十分性は、正の跳躍だけの構成で得た $Y$ の特性指数を比較すれば従う。補償なしの跳躍積分に一次項 $itd$ を足すと、もとの切断 $b$ に戻したときのドリフトは $d+\int b\,\dd\nu=\gamma$ である。

必要性では、$n$ 乗根の法則 $\mu_n$ は非負である。もし $\mu_n((-\infty,0))>0$ なら、独立な $n$ 項が全て負となる確率が正であり、和の非負性に反する。表示定理の証明で使った $\rho_n=nw\mu_n$ の任意の弱極限は従って非負の半直線に集中し、$\nu$ は負の側に質量を置かない。

$\ell(s)=\mathbb Ee^{-sX}$ とする。根の Laplace 変換は $\ell(s)^{1/n}$ なので $n(\ell(s)^{1/n}-1)\to\log\ell(s)$ である。表示定理で $\rho_n,\gamma_n$ が収束する部分列に沿い、そこでの特性関数の計算と同様に非負の半直線上で
$$
g_s(x)=\frac{e^{-sx}-1+s b(x)}{w(x)}
$$
を使う。$g_s$ は $x=0$ で $s^2/2$ を補えば有界連続であり、
$$
-\frac{\log\ell(s)}s
=\gamma-\frac{v s}{2}+
\int_{(0,\infty)}\left(\frac{1-e^{-sx}}{s}-b(x)\right)\nu(\dd x)
$$
を得る。左辺は $s>0$ で非負である。積分の被積分関数は $s$ とともに減少して $-b(x)$ に向かい、$s=1$ で絶対可積分である。$s=1$ の関数との差に単調収束を使うと、積分の極限は $-\int b\,\dd\nu$ である。従って $v=0$、$\int b\,\dd\nu\le\gamma$ が必要である。$0<x\le1$ では $b(x)\ge x/2$ なので、小跳躍の一次積率も有限となる。
:::

## 三角配列の極限定理

:::proposition[同分布行の反復畳み込み]{#prop-prob-infinitesimal-array}
行内で独立同分布の $Y_{n,1},\ldots,Y_{n,k_n}$ を考え、$k_n\to\infty$、$Y_{n,1}\to0$ が確率収束するとする。和の分布が弱収束するとき、その極限は無限分解可能である。逆に無限分解可能な任意の分布はこの形で表せる。
:::

:::proof
行の一項の特性関数を $f_n$、和の極限の特性関数を $f$ とする。任意のコンパクトな $t$ 区間で $f_n\to1$ は一様である。実際、確率収束で $\mathbb E\min(2,T|Y_{n,1}|)\to0$ となる。和の特性関数 $f_n^{k_n}$ も緊密性から $f$ へコンパクト上で一様収束する。これは質量の大きい有界区間では $e^{itx}$ が $t$ について一様 Lipschitz となるためである。$\phi_n=f_n$ に前節の余弦の不等式を使う。$t/m$ で $f$ が零でないなら $f_n(t/m)^{k_n}\to f(t/m)$ から $k_n(1-|f_n(t/m)|^2)$ が有界になり、$k_n(1-|f_n(t)|^2)$ も有界である。よって $f(t)\ne0$ である。

コンパクト区間上で $(f_n^{k_n})/f\to1$ が一様だから、原点で零の連続対数について $k_n\log f_n\to\log f$ となる。従って $\lfloor k_n/m\rfloor$ 項の和の特性関数は $e^{(\log f)/m}$ に収束する。この極限は原点で連続であり、連続性定理により $m$ 乗根を構成できる。逆向きは :ref[thm-prob-levy-khintchine] の三つ組を $k_n$ で割ればよい。根は特性関数の連続性定理で $\delta_0$ へ弱収束し、定数への弱収束は確率収束と同値である。
:::

例えば各項が確率 $\lambda/n$ で $Q$ に従う跳躍を生じ、残りで零なら、行和の特性関数は $(1+\lambda(\phi_Q-1)/n)^n$ である。極限は率 $\lambda$ の複合 Poisson になる。行内で分布が異なる一様微小配列の一般極限定理では、切断して中心化した各項の分散と Lévy 測度の収束を別々に判定する。各項の中心化も含む一般形は、この節の :ref[thm-prob-infinitesimal-array-general] で示す。

### 行内で異なる配列の十分条件

:::proposition[重みを付けた跳躍測度による極限の十分条件]{#prop-prob-array-weighted-convergence}
各行の独立な $Y_{n,k}$ について $f_{n,k}(t)=\mathbb E e^{itY_{n,k}}$ とする。各固定 $t$ で $\sum_k|f_{n,k}(t)-1|^2\to0$ とする。さらに
$$
\rho_n=\sum_k\frac{x^2}{1+x^2}\,\mathcal L(Y_{n,k})(\dd x),\qquad
\gamma_n=\sum_k\mathbb E\frac{Y_{n,k}}{1+Y_{n,k}^2}
$$
とし、有限測度 $\rho_n$ が有界連続関数の積分で $\rho$ へ収束し、$\gamma_n\to\gamma$ と仮定する。すると行和は $v=\rho(\{0\})$、$\nu(\dd x)=(1+x^2)x^{-2}\rho(\dd x)$（$x\ne0$）の Lévy–Khintchine 表現をもつ分布へ弱収束する。
:::

:::proof
二乗和の仮定から最大の $|f_{n,k}(t)-1|$ も零へ向かう。従って例えば半径 $1/2$ の $1$ の近傍で対数を取り、和における二次の余りを $C\sum_k|f_{n,k}(t)-1|^2$ で抑えられる。$w(x)=x^2/(1+x^2)$、$b(x)=x/(1+x^2)$ と書くと一次部分は
$$
\sum_k(f_{n,k}(t)-1)=it\gamma_n+
\int\frac{e^{itx}-1-itb(x)}{w(x)}\,\rho_n(\dd x)
$$
である。積分の関数は $x=0$ に $-t^2/2$ を補うと有界連続であり、有限測度の収束を使える。零の質量は正規分散を与え、零以外の部分は $\nu=\rho/w$ を与える。$(1\wedge x^2)/w(x)\le2$ なので積分条件を満たし、表現定理で極限は特性関数となる。連続性定理で弱収束が従う。
:::

:::example[小さな正規項と固定幅の跳躍を混ぜる]
行 $n$ は互いに独立な $n$ 個の $N(0,2/n)$ と、$n$ 個の $3B_{n,j}$（$B_{n,j}\sim\operatorname{Bern}(1/n)$）から成る。正規項の特性関数は $e^{-t^2/n}$、跳躍項は $1+(e^{3it}-1)/n$ である。いずれの二乗和も $O_t(1/n)$ なので、極限指数は $-t^2+e^{3it}-1$ となる。従って $v=2$、$\nu=\delta_3$、$\gamma=b(3)=3/10$ である。正規項はすべて原点近くへ集まり、固定幅 $3$ の項は回数だけが有限の尺度に残る。
:::

### 一様微小配列の一般的な判定

行の決定論的な移動を大きな跳躍として数えると、極限の分散を誤って評価する。例えば各項が確率 $1$ で $n^{-1/2}$ なら、二乗の和は $1$ だが、和から $\sqrt n$ を引けば零である。そのため、各項の小さな位置のずれを先に取り除く。

:::definition[一様微小性と切断中心]
独立な行 $(X_{n,k})_{k\le k_n}$ について、任意の $\eps>0$ で $\max_k\mathbb P(|X_{n,k}|>\eps)\to0$ を :term[一様微小性]{en="uniform asymptotic negligibility" yomi="いちようびしょうせい"} という。$h:\R\to\R$ を有界な Lipschitz 関数で $|x|\le1$ では $h(x)=x$ とする。$c_{n,k}=\mathbb Eh(X_{n,k})$、$Z_{n,k}=X_{n,k}-c_{n,k}$、$w(x)=x^2/(1+x^2)$ と置く。
:::

:::theorem[一様微小配列の収束判定]{#thm-prob-infinitesimal-array-general}
行が一様微小であるとする。定数 $d_n$ を足した和 $d_n+\sum_k X_{n,k}$ が弱収束するための必要十分条件は、次の二つがともに収束することである。
$$
\rho_n(\dd x)=\sum_k w(x)\mathcal L(Z_{n,k})(\dd x),\qquad
\gamma_n=d_n+\sum_k c_{n,k}+\sum_k\mathbb Eh(Z_{n,k}).
$$
有限測度 $\rho_n$ の収束は全ての有界連続関数の積分で取る。極限を $\rho,\gamma$ とすると、和の特性指数は
$$
it\gamma-\tfrac12 t^2\rho(\{0\})+
\int_{x\ne0}(e^{itx}-1-ith(x))w(x)^{-1}\rho(\dd x)
$$
である。
:::

:::proof
まず一様微小性と $h$ の有界性から $\max_k|c_{n,k}|\to0$ であり、$Z_{n,k}$ も一様微小である。十分後ろの行で $|c_{n,k}|\le1/2$ とする。$\mathbb Eh(X_{n,k})=c_{n,k}$ だから
$$
\mathbb Eh(Z_{n,k})=\mathbb E[h(X_{n,k}-c_{n,k})-h(X_{n,k})+c_{n,k}].
$$
被積分関数は $|Z_{n,k}|\le1/2$ では零、他では $C|c_{n,k}|$ 以下の絶対値をもつ。従って $|\mathbb Eh(Z_{n,k})|\le C'|c_{n,k}|\int w\,\dd\mathcal L(Z_{n,k})$ である。

$\sup_n\rho_n(\R)<\infty$ が成立する場合を考える。$g_t(x)=(e^{itx}-1-ith(x))/w(x)$ は零で $-t^2/2$ を補うと有界連続である。したがって $f_{n,k}(t)=\mathbb Ee^{itZ_{n,k}}$ に対し $\sum_k|f_{n,k}(t)-1|\le C_t\rho_n(\R)$ である。一様微小性で $\max_k|f_{n,k}(t)-1|\to0$ だから、差の二乗の和は零へ向かう。この議論はコンパクトな $t$ の範囲で一様にも行える。そこで対数を一次まで展開すると、全体の特性関数は、固定 $t$ で
$$
\exp\left(it\gamma_n+\int g_t\,\dd\rho_n+o(1)\right)
$$
と書ける。$\rho_n,\gamma_n$ がともに収束すれば、切断を $h$ に変えた Lévy–Khintchine 表現が得られ、連続性定理で十分性が従う。

必要性のため和が弱収束するとする。$X_{n,k}'$ を独立な複製、$D_{n,k}=X_{n,k}-X_{n,k}'$ とし、$\widetilde\rho_n=\sum_k w\,\mathcal L(D_{n,k})$ と置く。積の絶対値に定数移動は効かない。和の特性関数 $\Phi_n$ は緊密性により原点付近で一様に極限へ向かい、十分小さい $u>0$ なら $[-u,u]$ で零から離れる。この範囲では
$$
\sum_k\mathbb E[1-\cos(tD_{n,k})]
=\sum_k(1-|\mathbb Ee^{itX_{n,k}}|^2)
\le-2\log|\Phi_n(t)|.
$$
$[-u,u]$ で平均すると左辺は各項の $1-\sin(uD_{n,k})/(uD_{n,k})$ の期待値の和となる。この被積分関数は $c_u w(D_{n,k})$ 以上である。比は原点で正の極限、無限遠で $1$ をもち、その間でも正であるからである。従って $\widetilde\rho_n$ の全質量は有界である。さらに $|D|\ge2/u$ で被積分関数は $1/2$ 以上であり、$\sup_{|t|\le u}(-2\log|\Phi_n(t)|)$ の上極限は $u\downarrow0$ で零に向かう。よって $\widetilde\rho_n$ は緊密である。

この評価を $\rho_n$ に移す。$X=X_{n,k}$、$X'=X_{n,k}'$ と略記する。$h$ は有界 Lipschitz なので $(h(X)-h(X'))^2\le Cw(X-X')$ であり、$\operatorname{Var}(h(X))\le(C/2)\mathbb Ew(X-X')$ である。$|X|\le1$ では $Z=h(X)-\mathbb Eh(X)$、その外では $w(Z)\le1$ である。十分後ろでは $\mathbb P(|X'|\le1/2)\ge1/2$ であるから、$\mathbb P(|X|>1)\le2\mathbb P(|D|>1/2)\le C'\mathbb Ew(D)$ となる。これで $\rho_n$ の質量は有界となる。$R\ge2$ で $|Z|>R$ かつ $|X'|\le1/2$ なら $|D|>R-1$ である。そのため $\mathbb P(|Z|>R)\le2\mathbb P(|D|>R-1)$ となり、右辺を足して $\widetilde\rho_n$ の尾で抑えると $\rho_n$ の緊密性も得られる。

任意の部分列から、有限測度 $\rho_n$ が弱収束する部分列を選ぶ。全質量が収束するよう取り、極限が正なら正規化して Helly の定理を使い、零なら積分はすべて零へ向かう。この部分列で和から $\gamma_n$ を引くと、上の特性関数の展開と連続性定理により弱収束する。元の和も緊密なので、決定論的な差 $\gamma_n$ は有界である。実際、同じ標本空間で、両変数が共通の有界区間に入る確率をともに $3/4$ 以上とすれば、その共通部分は正確率である。従って $\gamma_n$ の収束部分列も取れる。

残るのはこんな取り方で極限が変わらないことを示すことである。$\rho$ の原点の質量は正規分散、零以外は $w\nu$ であるため、:ref[thm-prob-levy-khintchine] の一意性が $\rho$ を決める。$h-b$ は零付近で $O(|x|^3)$、離れたところで有界だから、その $\nu$ による積分も定まり、切断変更の式で $\gamma$ が一意に決まる。全部分列で同じ収束部分列が得られるので、元の列もこれらの極限に収束する。
:::

:::example[どちらの列が本当に分散を作るか]
一行の $n$ 項を全て $X_{n,k}=1/\sqrt n+\xi_k/\sqrt n$（公平符号 $\xi_k$）とする。十分大きい $n$ で切断中心は $c_{n,k}=1/\sqrt n$、中心化した項は $\pm1/\sqrt n$ である。$\rho_n$ は全質量 $n/(n+1)$ を $\pm1/\sqrt n$ に等分し、$\delta_0$ に弱収束する。$d_n=-\sqrt n$ と選べば $\gamma_n=0$ なので極限は標準正規となる。一方 $\xi_k$ を全て零にすれば $\rho_n=0$ であり、同じ $d_n$ の下で和は零となる。元の変数を中心化せず二乗だけ比較すると、後者でも総和 $1$ が現れてしまう。
:::

:::example[二種類の稀な跳躍を同じ行で足す]
偶数 $n$ の行に独立な $n$ 項を置き、各項に固定の $1/n$ を足す。前半の $n/2$ 項は確率 $1/n$ でさらに $-2$ を、後半の $n/2$ 項は確率 $2/n$ でさらに $1$ を足す。行和の特性関数は
$$
e^{it}\left(1+\frac{e^{-2it}-1}{n}\right)^{n/2}
\left(1+\frac{2(e^{it}-1)}{n}\right)^{n/2}.
$$
対数の一次展開で極限は $\exp(it+\tfrac12(e^{-2it}-1)+(e^{it}-1))$ となる。従って極限は $1-2N_1+N_2$ であり、$N_1,N_2$ は独立な率 $1/2,1$ の Poisson 変数である。平均は $1-2(1/2)+1=1$、分散は $4(1/2)+1=3$ である。負の跳躍と正の跳躍は一次積率では相殺しても、分散には両方が加わる。各行で異なる跳躍を区別しておくと、平均と分散が同じ別の分布にも取り違えない。
:::

## 安定分布とその特性関数

:::definition[安定分布]{#def-prob-stable}
非退化な $X$ の独立な複製 $X_1,X_2$ に対し、任意の $a,b>0$ で $aX_1+bX_2\stackrel d=cX+d$ となる $c>0,d\in\R$ が存在するとき、分布が :term[安定]{en="stable" yomi="あんてい"} であるという。いつも $d=0$ にできるものを狭義安定という。
:::

### 拡大係数を決める

:::lemma[独立な複製の個数と拡大係数]{#lem-prob-stable-scaling}
非退化な $X$ に対し、各 $n\ge1$ で独立な複製の和が $c_nX+d_n$ と同分布であるとする。ここで $c_n>0$ とする。このときある $H>0$ が存在し、$c_n=n^H$ である。
:::

:::proof
非退化な分布から相異なる値を取る二つの分位点を選べる。正の倍率で拡大しても、対応する二つの分位点の差はその倍率倍になる。従って $X$ とそれをアフィン変換した分布から、倍率と平行移動が一意に決まる。

$mn$ 個の複製を $n$ 個ずつの $m$ 組に分けると、全体の倍率は $c_mc_n$ である。一意性から $c_{mn}=c_mc_n$ となる。また、独立な複製 $X'$ を使うと $|\phi_X(t)|^2=\mathbb E\cos(t(X-X'))$ である。$0<|X-X'|\le M$ となる確率が正であるよう $M$ を選ぶと、$0<|t|<2\pi/M$ では $|\phi_X(t)|<1$ となる。従って十分小さい $s>0$ で
$$
r(s)=\inf\{t>0:|\phi_X(t)|\le e^{-s}\}
$$
は正かつ有限である。特性関数の絶対値を取れば $|\phi_X(c_nt)|=|\phi_X(t)|^n$ なので $r(s)/c_n=r(s/n)$ である。$r$ は非減少だから $c_n$ も非減少であり、連続性と原点の近くの厳密不等式から $r(s/n)\downarrow0$ なので $c_n\to\infty$ となる。

$a,b\ge2$ に対し $a^{\lfloor k\log b/\log a\rfloor}\le b^k<a^{\lfloor k\log b/\log a\rfloor+1}$ である。非減少性と積の公式で $c_b^k$ を両側から $c_a$ のべきで挟み、対数を $k$ で割って極限を取れば $\log c_b/\log b=\log c_a/\log a$ である。この共通値は $c_n\to\infty$ により正である。
:::

:::proposition[安定分布の Lévy 測度]{#prop-prob-stable-levy}
安定分布には指数 $\alpha\in(0,2]$ があり、独立和の係数は $c=(a^\alpha+b^\alpha)^{1/\alpha}$ となる。$0<\alpha<2$ の特性指数は :ref[thm-prob-levy-khintchine] で $v=0$、
$$
\nu(\dd x)=(c_+\mathbf1_{\{x>0\}}+c_-\mathbf1_{\{x<0\}})|x|^{-1-\alpha}\dd x
$$
としたものであり、$c_+,c_-\ge0$、$c_++c_->0$ である。逆にこの形の測度と任意の $\gamma$ は指数 $\alpha$ の安定分布を与える。指数 $2$ の安定分布は分散が正の正規分布である。
:::

:::proof
安定性を繰り返すと同じ分布の $n$ 項の和が $c_nX+d_n$ と書ける。:ref[lem-prob-stable-scaling] により $c_n=n^{1/\alpha}$ と書ける（ここでは $\alpha>0$）。この畳み込みの式を $X$ について解けば、独立同分布の $n$ 項で $X$ が実現できる。従って分布は無限分解可能である。

一次の項以外を比較すると、正規分散は $nv=n^{2/\alpha}v$ となる。Lévy 測度の一意性から $\nu(rB)=r^{-\alpha}\nu(B)$ が $r=n^{1/\alpha}$ で成り立ち、比を取って $r=(m/n)^{1/\alpha}$ に拡張できる。片側の尾測度は単調で、稠密な $r$ の両側から挟めば $\nu((x,\infty))=C_+x^{-\alpha}$、$\nu((-\infty,-x))=C_-x^{-\alpha}$（$x>0$）になる。微分すれば所望の密度になる。$\alpha\ge2$ では零付近の二次積分が発散するため $C_\pm=0$ であり、非退化性には $v>0$ が必要となる。分散の等式で $\alpha=2$ となる。$0<\alpha<2$ なら分散の等式で $v=0$ になる。

逆にこの測度では $aX_1+bX_2$ と $(a^\alpha+b^\alpha)^{1/\alpha}X$ の正規分散と Lévy 測度が一致する。従って特性指数の差は一次項だけであり、二つの分布は定数移動で一致する。正規の場合も分散から同じ倍率が出る。倍率は分位点の差で一意なので、これで $c$ が決まる。証明には任意の $a,b$ での安定性のうち整数個の等式だけを使ったため、それらを満たす非退化分布も既に安定である。
:::

対称な場合は $c_+=c_-=c$ であり、虚部の跳躍積分が消える。$u=|t|x$ の置換から $\phi(t)=\exp(it\gamma-A|t|^\alpha)$、$A=2c\int_0^\infty(1-\cos u)u^{-1-\alpha}\dd u>0$ となる。特に $\alpha=1$ では $\exp(-A|t|)$ の逆変換が $A/(\pi(A^2+x^2))$ であり、対称 Cauchy 分布を与える。$\alpha=1$ で非対称な場合は切断項の拡大による一次補正が必要であり、対称例から $d=0$ を一般化できない。

### 有限平均と片側の安定分布

:::proposition[有限の二次積率をもつ安定分布]
非退化な安定分布で分散が有限なら指数は $2$ であり、分布は正規である。また有限平均をもつ安定分布は、平均を引くと狭義安定になる。
:::

:::proof
独立な二つの和 $X_1+X_2$ の分散は $2\operatorname{Var}(X)$ であり、一コピーの拡大では $2^{2/\alpha}\operatorname{Var}(X)$ である。正の分散で割ると $\alpha=2$ となり、:ref[prop-prob-stable-levy] から正規分布である。平均 $m$ が有限なら $aX_1+bX_2\stackrel d=cX+d$ の平均を比べると $d=(a+b-c)m$ である。両辺から $(a+b)m$ を引けば、$X-m$ の平行移動が消える。
:::

:::example[正の安定分布の Laplace 指数]
$0<\alpha<1$ として $\nu(\dd x)=c x^{-1-\alpha}\mathbf1_{\{x>0\}}\,\dd x$ とする。小跳躍の絶対値は可積分であり、大跳躍は有限強度なので、補正なしの正の和が作れる。その Laplace 指数は $-c\int_0^\infty(1-e^{-sx})x^{-1-\alpha}\dd x$ である。$y=sx$ と置き、部分積分すると $\int_0^\infty(1-e^{-y})y^{-1-\alpha}\dd y=\Gamma(1-\alpha)/\alpha$ である。境界では零付近が $O(y^{1-\alpha})$、無限遠が $O(y^{-\alpha})$ なので境界項は消える。例えば $\alpha=1/2$、$c=1/(2\sqrt\pi)$ なら指数は $-\sqrt s$ となる。二コピーを足すと指数は $-2\sqrt s=-\sqrt{4s}$ であり、一コピーの $4$ 倍と同分布である。
:::

### 指数が一の非対称な場合の対数項

:::proposition[指数一で拡大したときの移動量]
$\alpha=1$ の安定分布で、切断を $h(x)=x\mathbf1_{\{|x|\le1\}}$ とし、Lévy 密度の係数を $c_+,c_-$ とする。特性指数 $\psi$ は $r>0$ に対して
$$
\psi(rt)=r\psi(t)-irt(c_+-c_-)\log r
$$
を満たす。
:::

:::proof
$y=rx$ と置くと $\nu$ は $r$ 倍で現れる。二つの積分を引くと、一次の切断部分の差だけが残り、
$$
\psi(rt)-r\psi(t)=-irt\int y(\mathbf1_{\{|y|\le r\}}-\mathbf1_{\{|y|\le1\}})\,\nu(\dd y)
$$
となる。これは零から離れた環状領域での積分である。$r\ge1$ なら正の側は $c_+\int_1^r y^{-1}\dd y$、負の側は $-c_-\int_1^r y^{-1}\dd y$ である。$r<1$ では区間の向きが反転し、同じ $\log r$ を得る。ドリフトの項はちょうど相殺する。
:::

:::example[二コピーの和を二倍に戻す]
$c_+-c_-=1$ の指数一の分布では、$X_1+X_2$ の指数 $2\psi(t)$ は $\psi(2t)+2it\log2$ に等しい。従って $X_1+X_2\stackrel d=2X+2\log2$ である。$X$ をどの定数だけ平行移動しても、和と二倍の平行移動はともにその定数の二倍なので、$2\log2$ の差は残る。対称 Cauchy のように $c_+=c_-$ なら差は零となる。
:::

### 安定分布の尺度と非対称性を特性関数で表す

$\alpha\ne1$ では、Lévy 測度の係数を使って特性指数を
$$
\psi(t)=imt-A|t|^\alpha
 \left(1-i\beta\operatorname{sgn}(t)\tan(\pi\alpha/2)\right)
$$
と書ける。ただし $0<\alpha<2$、$\beta=(c_+-c_-)/(c_++c_-)$、$A>0$ であり、$m$ は切断に応じて決まる位置の母数である。実際、$\alpha<1$ では正負それぞれの積分で一次補正を外し、$\alpha>1$ では一次補正を全域の $itx$ に変えてよい。変更を積分した値は有限で、一次項に吸収できる。

正の側の積分を計算するため、まず実数 $s>0$ に対し $e^{-sx}$ を使う。$\alpha<1$ では $\int_0^\infty(e^{-sx}-1)x^{-1-\alpha}\dd x$、$\alpha>1$ では $\int_0^\infty(e^{-sx}-1+sx)x^{-1-\alpha}\dd x$ を考える。$u=sx$ と置き、前者は一回、後者は二回部分積分すれば、どちらも $C_\alpha s^\alpha$ となる。ただし
$$
C_\alpha=\frac{\Gamma(2-\alpha)}{\alpha(\alpha-1)}.
$$
第一の境界項は、零付近では $e^{-u}-1=O(u)$、無限遠では有界であるため消える。第二では $e^{-u}-1+u=O(u^2)$（零付近）、$O(u)$（無限遠）であり、一回目の積分後も $1-e^{-u}=O(u)$（零付近）、$O(1)$（無限遠）なので境界項は消える。

これらの等式は $\Re s>0$ で正則性と一致の定理により成り立つ。右辺のべきは偏角を $(-\pi/2,\pi/2)$ に取る。積分は $s$ が有界な閉右半平面上で、零付近の Taylor 評価と無限遠の評価によって優収束を使える。従って $s$ を $-it$ に近づけると、正の側の Fourier 積分は $C_\alpha|t|^\alpha e^{-i\operatorname{sgn}(t)\pi\alpha/2}$ である。負の側では符号を反転する。二つを足して実部と虚部に分けると上の式を得る。ここで $A=-C_\alpha(c_++c_-)\cos(\pi\alpha/2)>0$ である。

$\alpha=1$ では、前節の拡大公式に $t=1,-1$ を入れて、$t\ne0$ で
$$
\psi(t)=im t-A|t|-i(c_+-c_-)t\log|t|
$$
と書ける。ここで $A=-\Re\psi(1)=(c_++c_-)\pi/2$ である。最後の積分定数は $\int_0^\infty(1-\cos u)u^{-2}\dd u$ を部分積分して $\int_0^\infty\sin(u)/u\,\dd u$ に帰着すれば、区間反転公式の証明の通り $\pi/2$ になる。$t\log|t|$ は $t=0$ で零とする。従って指数一の非対称性では、定数移動で吸収できない対数項が残る。

### 正の指数二分の一の安定分布を反転する

:::example[正規変数の逆二乗から作る安定分布]
$a>0$ とし、標準正規変数 $Z$ から $T_a=a^2/Z^2$ を作る。$x>0$ の密度は変数変換で $a(2\pi x^3)^{-1/2}e^{-a^2/(2x)}$ である。Laplace 変換を調べるため $c=a\sqrt{2s}$ と置くと、
$$
\mathbb E e^{-sT_a}
=\sqrt{\frac2\pi}\int_0^\infty\exp\!\left[-\frac{z^2+c^2/z^2}{2}\right]\dd z.
$$
この積分を $J(c)$ とする。$c>0$ で積分内を微分し、$w=c/z$ と置換すると
$$
J'(c)=-\int_0^\infty\frac{c}{z^2}\exp\!\left[-\frac{z^2+c^2/z^2}{2}\right]\dd z=-J(c)
$$
となる。微分は $c$ の正のコンパクト区間で、零付近の $e^{-\mathrm{const}/z^2}$ と無限遠の $e^{-z^2/2}$ により正当化できる。$c\downarrow0$ では優収束で $J(c)\to\sqrt{\pi/2}$ なので、変換は $e^{-a\sqrt{2s}}$ である。

独立な $T_a,T_b'$ の Laplace 変換を掛けると $T_{a+b}$ の変換になる。例えば独立な $T_1$ 二つの和は $T_2\stackrel d=4T_1$ と同分布である。密度の無限遠での減衰は $x^{-3/2}$ の定数倍なので、平均は無限である。この分布が到達時刻として現れる例は [Brown 運動](/math/probability/brownian-motion/) で扱う。
:::

## 吸引領域と一般化中心極限定理

:::definition[吸引領域]{#def-prob-domain-of-attraction}
独立同分布な $X_k$ について、定数 $a_n>0,b_n$ で $(\sum_{k\le n}X_k-b_n)/a_n$ が非退化な分布 $\mu$ へ弱収束するとき、$X_1$ の分布は $\mu$ の :term[吸引領域]{en="domain of attraction" yomi="きゅういんりょういき"} に属するという。
:::

:::theorem[べきの裾をもつ非正規極限]{#thm-prob-stable-attraction}
$0<\alpha<2$ とし、$x\to\infty$ で $x^\alpha\mathbb P(X>x)\to A_+$、$x^\alpha\mathbb P(X<-x)\to A_-$、$A_++A_->0$ とする。$a_n=n^{1/\alpha}$、$b_n=a_nn\mathbb E[b(X/a_n)]$、$b(x)=x/(1+x^2)$ とすると、中心化した和は $\gamma=v=0$、$c_\pm=\alpha A_\pm$ の安定分布へ弱収束する。
:::

:::proof
$\nu_n=n\,\mathbb P(X/a_n\in\cdot)$ とおく。零から離れた二つの端点の尾の差を取ると、対応する区間の $\nu_n$ の質量は候補 $\nu$ の質量に収束する。区間分割で、零を含まないコンパクト集合に台をもつ連続関数の積分にも拡張できる。残りを抑える。$\nu_n(|x|>R)\le CR^{-\alpha}$ であり、$\int_{|x|\le\eps}x^2\,\dd\nu_n\le 2n\int_0^\eps y\,\mathbb P(|X|>a_ny)\dd y\le C'\eps^{2-\alpha}$ である。後者には全 $u>0$ で $\mathbb P(|X|>u)\le C u^{-\alpha}$ となる上界を使った。従って $g(x)=e^{itx}-1-itb(x)$ の積分は $\nu$ の積分へ収束する。

さらに $|f_n(t)-1|\le C_t\mathbb E(1\wedge|X|/a_n)$ である。尾積分を $1/a_n$ のところで分ければ、この上界は $\alpha<1$ で $O(n^{-1})$、$\alpha=1$ で $O(n^{-1}\log n)$、$\alpha>1$ で $O(n^{-1/\alpha})$ である。$\alpha<2$ なので $n|f_n(t)-1|^2\to0$ となり、$n\log f_n$ と $n(f_n-1)$ の差も零に収束する。$b_n$ で一次項を差し引くと極限の特性指数が得られ、連続性定理で結論となる。
:::

例えば $X$ が対称で $\mathbb P(|X|>x)=x^{-3/2}$（$x\ge1$）なら平均は零であり、正規化は $n^{2/3}$ である。$\sqrt n$ で割ると正規極限は得られない。

### 独立同分布和の極限が安定になる理由

:::lemma[正規化の収束]{#lem-prob-convergence-of-types}
$Y_n\Rightarrow Y$、$r_nY_n+d_n\Rightarrow Z$ で $r_n>0$、$Y,Z$ は非退化とする。このとき $(r_n,d_n)$ は $(r,d)$（$r>0$）へ収束し、$Z\stackrel d=rY+d$ である。
:::

:::proof
$Y$ の分位点関数の連続点 $u<v$ を、$Q_Y(u)<Q_Y(v)$ となるよう選ぶ。$Z$ の分位点関数の不連続点は高々可算なので、$u,v$ はその連続点でもあるよう選べる。:ref[thm-prob-skorokhod-real] の証明にある分位点収束により、$Y_n$ の二点の分位点の差は正の有限値へ収束する。$r_nY_n+d_n$ の分位点の差はその $r_n$ 倍なので、$r_n$ は非負の有限値 $r$ へ収束する。また一方の分位点の式から $d_n$ も有限値 $d$ へ収束する。$Y_n$ は緊密であり、仮に $r=0$ なら第二の極限は $\delta_d$ となるため $r>0$ である。各有界連続関数に対し、共通のコンパクト区間上で $f(r_nx+d_n)\to f(rx+d)$ が一様となり、その外の積分を緊密性で抑えられる。よって $Z\stackrel d=rY+d$ である。
:::

:::theorem[一般化中心極限定理の必要条件]{#thm-prob-generalized-clt}
独立同分布列の部分和のアフィン正規化が非退化な分布に弱収束するなら、極限分布は安定である。
:::

:::proof
$Y_n=(S_n-b_n)/a_n\Rightarrow Y$ とする。$m$ を固定し、$mn$ 項を長さ $n$ の $m$ ブロックに分けると、$(S_{mn}-mb_n)/a_n$ は独立な $Y_n$ の $m$ 個の和と等分布である。特性関数を掛け合わせると、その極限は独立な $Y$ の $m$ 個の和である。この和は非退化である。もし和が定数なら $|\phi_Y|^m=1$ となるが、:ref[lem-prob-stable-scaling] の証明にある小さい非零 $t$ での厳密不等式に反する。

一方 $(S_{mn}-mb_n)/a_n=(a_{mn}/a_n)Y_{mn}+(b_{mn}-mb_n)/a_n$ であり、$Y_{mn}\Rightarrow Y$ である。従って :ref[lem-prob-convergence-of-types] により、独立な $Y$ の $m$ 項の和は $Y$ のアフィン変換と同分布である。全 $m$ に対するこの等式は :ref[prop-prob-stable-levy] の証明の仮定を満たすため、$Y$ は安定である。
:::

### 吸引領域の一般的な判定

正の関数 $L$ が全 $c>0$ で $L(cx)/L(x)\to1$ を満たすとき、無限遠で緩変動という。$0<\alpha<2$ の吸引領域の必要十分条件は、$\mathbb P(|X|>x)=x^{-\alpha}L(x)$ がこの $L$ をもち、$\mathbb P(X>x)/\mathbb P(|X|>x)$ が極限 $p\in[0,1]$ をもつことである。このとき $n\mathbb P(|X|>a_n)\to1$ となる $a_n$ を選び、極限の Lévy 測度の左右の尾の係数は $p,1-p$ になる。指数 $2$ では尾自体でなく、切断二次積率 $\mathbb E[X^2\mathbf1_{\{|X|\le x\}}]$ の緩変動性が正規吸引領域の条件になる。指数が $2$ 未満の必要十分性は :ref[thm-prob-stable-attraction-regular-variation] で示す。指数 $2$ では有限分散の場合は通常の中心極限定理、分散が無限の場合は :ref[prop-prob-normal-attraction-truncated-variance] が十分性を与える。必要性は :ref[prop-prob-normal-attraction-necessity] で示す。

### 正規の吸引領域で分散が無限になる例

:::example[対数補正が必要となる正規極限]
$X$ の符号が独立に公平で、$\mathbb P(|X|>x)=x^{-2}$（$x\ge1$）とする。$|X|$ の密度は $2x^{-3}$（$x>1$）なので $\mathbb E[X^2\mathbf1_{\{|X|\le y\}}]=2\log y$（$y\ge1$）であり、分散は無限である。$a_n=\sqrt{n\log n}$ とし、行で各 $X_k$ を $\{|X_k|\le a_n\}$ で切断する。どれかが変わる確率は $n/a_n^2=1/\log n\to0$ である。

切断和を $a_n$ で割ると、分散は $n\,2\log a_n/a_n^2=1+\log\log n/\log n\to1$ である。固定 $0<\eps<1$ で $\eps a_n>1$ なら、Lindeberg の尾二次積率和は $n\cdot2\log(1/\eps)/a_n^2\to0$ である。従って :ref[thm-prob-lindeberg-feller] により正規極限を得る。切断による変更確率も零へ向かうので、元の和 $S_n/a_n$ も標準正規へ収束する。ここでは正規という極限の形は保たれるが、切断二次積率の増加が正規化に現れる。
:::

### 緩変動を含む非正規吸引領域の判定

:ref[thm-prob-stable-attraction] では尾に $x^\alpha$ を掛けた値が有限の係数へ収束する場合を扱った。尾に対数因子などがあるときは、べきの指数を保ちつつ、規格化の長さを尾の分位点に合わせる。

:::theorem[指数が二未満の吸引領域]{#thm-prob-stable-attraction-regular-variation}
$H(x)=\mathbb P(|X|>x)$ とし、$X$ は非退化とする。$0<\alpha<2$ の安定分布への吸引領域に属するための必要十分条件は、十分大きい $x$ で $H(x)>0$ であり、次の二つの極限が成り立つことである。
$$
H(cx)/H(x)\longrightarrow c^{-\alpha}\quad(c>0),\qquad
\mathbb P(X>x)/H(x)\longrightarrow p\in[0,1]
$$
十分性では $nH(a_n)\to1$ となる $a_n$ を選ぶと、極限の Lévy 測度の尾は $p u^{-\alpha}$ と $(1-p)u^{-\alpha}$（$u>0$）である。
:::

:::proof[尾の条件からの十分性]
まず $a_n$ として $H(x)\le1/n$ になる最初の水準を取る。右連続性により $H(a_n)\le1/n$ であり、任意の $0<r<1$ で $H(ra_n)>1/n$ が十分大きい $n$ で成立する。$a_n\to\infty$ であり、$H(a_n)/H(ra_n)\to r^\alpha$ だから $\liminf nH(a_n)\ge r^\alpha$ である。$r\uparrow1$ で規格化の式を得る。

$X/a_n$ の法則の $n$ 倍を $\nu_n$ とする。固定 $u>0$ で、例えば
$$
\nu_n((u,\infty))
=nH(a_n)\frac{H(ua_n)}{H(a_n)}
 \frac{\mathbb P(X>ua_n)}{H(ua_n)}
\longrightarrow p u^{-\alpha}.
$$
負の側も同じ計算になる。区間の差から、零から離れたコンパクト台上での連続関数の積分が候補の Lévy 測度へ収束する。候補の尾は連続で、有限個の区間への近似で境界の問題を避けられる。

尾を一様に抑えるため $0<\eta<\min(\alpha,2-\alpha)$ を取る。十分大きな $x_0$ と $x\ge x_0$ について
$$
2^{-\alpha-\eta}\le H(2x)/H(x)\le2^{-\alpha+\eta}
$$
である。反復して、$y\ge1$ では $nH(a_ny)\le C y^{-\alpha+\eta}$ となる。詳しくは $2^j\le y<2^{j+1}$ として $H(a_n y)\le H(2^j a_n)$ を使う。$0<y<1$ については、$ya_n$ と $a_n$ の間を倍にしながら比較し、$ya_n\ge x_0$ なら $nH(a_ny)\le C y^{-\alpha-\eta}$ となる。最後の一段が $a_n$ をまたいでも、$nH(a_n)$ の有界性と一回分の定数で吸収できる。さらに同じ下側の二進評価を $x_0$ から $a_n$ まで繰り返すと $H(a_n)\ge c a_n^{-\alpha-\eta}$ である。従って $n\le C'a_n^{\alpha+\eta}$ であり、$ya_n<x_0$ では $H\le1$ と合わせて、定数を変えれば同じ上界を得る。初期の有限行は極限には不要である。

この評価から $\nu_n(|x|>R)\le CR^{-\alpha+\eta}$ となる。また
$$
\int_{|x|\le\eps}x^2\,\nu_n(\dd x)
\le 2n\int_0^\eps yH(a_ny)\,\dd y
\le C'\eps^{2-\alpha-\eta}.
$$
従って $g_t(x)=e^{itx}-1-itx/(1+x^2)$ を積分するとき、原点近くと無限遠の誤差を一様に小さくできる。中間部では先に得た区間による収束を使えるから、$\int g_t\,\dd\nu_n\to\int g_t\,\dd\nu$ である。

最後に、$\alpha/2<\beta<\min(\alpha,1)$ を取る。必要なら $\eta$ をさらに小さくして $\alpha+\eta<2\beta$、$\beta<\alpha-\eta$ とする。二進評価の上側から $H(x)\le Cx^{-\alpha+\eta}$ が得られ、尾積分で $\mathbb E|X|^\beta<\infty$ である。$f_n(t)=\mathbb Ee^{itX/a_n}$ に対し $|f_n(t)-1|\le C_t a_n^{-\beta}\mathbb E|X|^\beta$ なので、$n|f_n(t)-1|^2\to0$ である。従って $n\log f_n(t)-n(f_n(t)-1)\to0$ となる。中心化を $b_n=a_nn\,\mathbb E[(X/a_n)/(1+(X/a_n)^2)]$ と選べば一次項が消える。極限の特性指数は上の $\nu$ を持ち、正規成分は零である。Lévy の連続性定理と :ref[prop-prob-stable-levy] で安定な極限分布が得られる。
:::

:::proof[安定極限からの必要性]
$(S_n-b_n)/a_n$ が非退化な指数 $\alpha$ の安定分布へ収束するとする。まず $a_n\to\infty$ を示す。もし $a_n$ が有界な部分列があれば、その上で独立なコピーとの差 $S_n-S_n'$ は緊密である。実際、標準化した二つの変数は緊密であり、差を有界な $a_n$ 倍しても緊密性が残る。一方その特性関数は $|\phi_X(t)|^{2n}$ である。非退化性により、原点に十分近い非零の $t$ では $|\phi_X(t)|<1$ であるから、この部分列の極限は原点で不連続になる。緊密性による部分列の弱収束に反する。従って $X/a_n\to0$ が確率収束する。

次に $(S_{n+1}-b_n)/a_n$ は $(S_n-b_n)/a_n$ と一項 $X_{n+1}/a_n$ の和なので、同じ弱極限を持つ。これは $a_{n+1}/a_n$ 倍した $(S_{n+1}-b_{n+1})/a_{n+1}$ に定数 $(b_{n+1}-b_n)/a_n$ を足したものでもある。:ref[lem-prob-convergence-of-types] の倍率と移動の一意性により $a_{n+1}/a_n\to1$ となる。

一様微小配列の必要条件を行内同分布の配列 $X/a_n$ に用いる。切断中心 $c_n$ は零へ向かい、その中心を引いた分布の $n$ 倍に $w$ を掛けた有限測度は、極限の $w\nu$ と原点の成分へ弱収束する。固定 $u>0$ では $(u,\infty)$ 上で $1/w$ は有界であり、境界 $u$ に極限測度の原子はない。よって $n\mathbb P(X/a_n-c_n>u)$ は $\nu((u,\infty))=C_+u^{-\alpha}$ に収束する。これは有界連続な上・下近似を境界の小区間で使って確かめられる。$c_n\to0$ なので $u\pm\eps$ で挟み、$\eps\downarrow0$ とすれば $c_n$ を外せる。負の側にも同じ操作を行い、
$$
nH(a_nu)\to(C_++C_-)u^{-\alpha},\qquad
n\mathbb P(X>a_nu)\to C_+u^{-\alpha}
$$
を得る。非退化な指数 $\alpha<2$ の安定分布では $C_++C_->0$ である。

最後に任意の大きい $x$ に対し、$a_n\ge x$ となる最初の $n$ を選ぶ。$a_n\to\infty$ なのでこの $n$ は存在し、$x\to\infty$ で $n\to\infty$ である。$a_{n-1}<x\le a_n$ と $a_{n-1}/a_n\to1$ から $x/a_n\to1$ となる。任意の固定 $u>0$ について $ux$ を $(u-\eps)a_n$ と $(u+\eps)a_n$ で挟み、尾の単調性と上の極限を用いる。$\eps\downarrow0$ で $nH(ux)\to(C_++C_-)u^{-\alpha}$ が得られ、正の尾も対応する極限を持つ。$u=1$ の全尾で割ると二つの必要条件が従う。
:::

:::example[対数因子に合わせて正規化を変更する]
対称な $X$ で、十分大きい $x$ に対し $H(x)=x^{-3/2}\log x$ とする。この関数は十分大きい $x$ で単調減少し、尾の残りを補って確率分布を作れる。固定 $c>0$ で $\log(cx)/\log x\to1$ なので指数は $3/2$、尾の比は $p=1/2$ である。$a_n=(\tfrac23 n\log n)^{2/3}$ を取ると
$$
nH(a_n)
=\frac{(2/3)\log((2/3)n\log n)}{(2/3)\log n}
\longrightarrow1.
$$
従ってこの尺度で対称安定極限を得る。$n^{2/3}$ だけを使うと $nH(n^{2/3})=(2/3)\log n$ が発散し、幅を固定した大きい跳躍が期待値の上で増え続ける。べき指数の情報だけでは規格化を決められない。
:::

### 正規吸引領域の切断二次積率による十分性

:::proposition[緩変動な切断二次積率からの正規極限]{#prop-prob-normal-attraction-truncated-variance}
非退化な $X$ の切断二次積率 $V(x)=\mathbb E[X^2\mathbf1_{\{|X|\le x\}}]$ が無限遠で緩変動し、$V(x)\to\infty$ とする。このとき $nV(a_n)/a_n^2\to1$ を満たす $a_n\to\infty$ を選べ、$X$ の独立同分布和について
$$
\frac{S_n-n\mathbb E X}{a_n}\Rightarrow N(0,1)
$$
である。
:::

:::proof
十分大きい $x$ で $V(x)>0$ である。$0<\eta<1$ を固定すると $V(2x)\le2^\eta V(x)$ が十分大きい $x$ で成立する。従って $y\ge x$ で $V(y)\le C_\eta (y/x)^\eta V(x)$ である。尾を $[2^jx,2^{j+1}x]$ の環に分けると、
$$
\frac{x^2\mathbb P(|X|>x)}{V(x)}
\le\sum_{j\ge0}2^{-2j}
   \frac{V(2^{j+1}x)-V(2^jx)}{V(x)}.
$$
各固定 $j$ の項は緩変動性で零へ向かう。右辺は $C_\eta\sum_j2^{-(2-\eta)j}$ で支配できるので、左辺は零へ向かう。この尾評価と $V(x)=O(x^\eta)$ により、尾積分で $\mathbb E|X|<\infty$ となる。

$V_*(x)=\mathbb E(X^2\wedge x^2)=V(x)+x^2\mathbb P(|X|>x)$ は $V(x)$ に漸近する。$V_*(x)/x^2=\mathbb E((X^2/x^2)\wedge1)$ は連続非増加で零へ向かう。$X$ が恒等的に零ではないので、任意の有限 $x$ で正である。そのため大きい $n$ で $nV_*(a_n)/a_n^2=1$ の解を取れ、$a_n\to\infty$ である。従って $nV(a_n)/a_n^2\to1$ となる。尾評価から $n\mathbb P(|X|>a_n)\to0$ であり、各行を $a_n$ で切断しても和が変わる確率は零へ向かう。

$Y_n=X\mathbf1_{\{|X|\le a_n\}}$ とする。さらに固定 $\eps>0$ に対し
$$
\frac{n}{a_n^2}\bigl(V(a_n)-V((\eps/2)a_n)\bigr)\longrightarrow0
$$
である（$\eps<2$ を選べば十分である）。これは緩変動性による。$\mathbb E Y_n$ は有界なので、中心化した変数の尾 $|Y_n-\mathbb E Y_n|>\eps a_n$ は大きい $n$ では $|Y_n|>(\eps/2)a_n$ に含まれる。二乗は $2Y_n^2+2(\mathbb E Y_n)^2$ 以下であり、第二項を尾の上で積分しても $|Y_n|>(\eps/2)a_n$ では第一項の定数倍で抑えられる。従って規格化した Lindeberg の尾和は零へ向かう。

仮定 $V(x)\to\infty$ により $n(\mathbb E Y_n)^2/a_n^2\to0$ なので、切断和の分散比は $1$ に収束する。

$V(x)\to\infty$ の場合、中心の入れ替えを確かめる。任意の $\delta>0$ で十分大きい $y$ に対し $\mathbb P(|X|>y)\le\delta V(y)/y^2$ である。$y\ge a_n$ で二進評価を使い、$\mathbb E[|X|\mathbf1_{\{|X|>a_n\}}]$ を $a_nH(a_n)+\int_{a_n}^\infty H(y)\dd y$ と書く。すると $\delta C_\eta V(a_n)/a_n$ 以下になる。これに $n/a_n$ を掛けた値の上極限は任意の $\delta>0$ で抑えられるから零である。従って切断平均を真の平均に戻しても極限は変わらない。Lindeberg–Feller と切断の誤差の評価で結論を得る。
:::

:::example[尾に対数がある指数二の例]
対称な $X$ で十分大きい $x$ において $\mathbb P(|X|>x)=x^{-2}/\log x$ とする。尾積分と端点の補正から $V(x)=2\log\log x+O(1)$ である。実際、$2\int^x yH(y)\dd y$ の可変部分は $2\log\log x$ であり、$x^2H(x)=1/\log x$ を引いても主項は同じである。従って $a_n=\sqrt{2n\log\log n}$ を取ると、$\log\log a_n/\log\log n\to1$ なので $nV(a_n)/a_n^2\to1$ である。平均は対称性により零であり、この尺度で標準正規極限を得る。尾が $x^{-2}$ の場合の $\sqrt{n\log n}$ と比べると、対数因子が小さいため必要な正規化も小さい。
:::

### 正規極限から切断二次積率の条件を導く

:::proposition[正規吸引領域の必要条件]{#prop-prob-normal-attraction-necessity}
独立同分布な非退化実変数の和が、決定論的な中心化と正の尺度により非退化な正規分布へ弱収束するなら、$V(x)=\mathbb E[X^2\mathbf1_{\{|X|\le x\}}]$ は無限遠で緩変動する。
:::

:::proof
極限の中心と分散を定数で調整し、標準正規を極限とする。:ref[thm-prob-stable-attraction-regular-variation] の必要性の証明の前半は、極限の安定指数に依存しない。その議論から $a_n\to\infty$ と $a_{n+1}/a_n\to1$ である。独立な複製との差 $D=X-X'$ に対しては和を $a_n$ で割ると $N(0,2)$ へ弱収束する。$D$ は対称なので、:ref[thm-prob-infinitesimal-array-general] の切断 $h$ を奇関数に選べば中心 $c_{n,k}$ は零である。この定理から、有限測度 $n w(x)\mathcal L(D/a_n)(\dd x)$ は $2\delta_0$ へ弱収束する。固定 $u>0$ で、関数 $(1+x^2)\mathbf1_{\{|x|\le u\}}$ は有界で、極限の質量を持つ原点で連続であるため積分も収束する。従って $V_D(x)=\mathbb E[D^2\mathbf1_{\{|D|\le x\}}]$ に対し
$$
\frac{nV_D(ua_n)}{a_n^2}\to2\qquad (u>0)
$$
である。$x\to\infty$ に対し初めて $a_n\ge x$ となる $n$ を取り、$a_{n-1}<x\le a_n$ から $x/a_n\to1$ を使う。$ux$ を $(u\pm\eps)a_n$ で挟み、$V_D$ の単調性と上の収束を使えば $V_D(ux)/V_D(x)\to1$ である。

前の十分性の証明にある環状領域の評価を $V_D$ に使うと、$x^2\mathbb P(|D|>x)/V_D(x)\to0$ が従う。この評価は $V_D$ が無限へ発散すると仮定せずに成立する。また独立性から
$$
\mathbb P(|X|>2x)\mathbb P(|X'|\le x)\le\mathbb P(|D|>x)
$$
なので、$x^2\mathbb P(|X|>x)/V_D(x)\to0$ である。$V_D(x)=O(x^\eta)$（$0<\eta<1$）を同じ二進評価から得れば、尾積分で $\mathbb E|X|<\infty$ となる。

$\widetilde X=X-\mathbb EX$ と置く。$\mathbb P(|\widetilde X|>x)\le\mathbb P(|X|>x/2)$ が十分大きい $x$ で成り立つため、この尾にも同じ評価が使える。$W(x)=\mathbb E[\widetilde X^2\mathbf1_{\{|\widetilde X|\le x\}}]$ と置く。$|\widetilde X|,|\widetilde X'|\le x/2$ なら $|D|\le x$ である。その領域で $D^2$ の積分を二乗展開し、独立性を使うと
$$
2\mathbb P(|\widetilde X|\le x/2)W(x/2)
-2\left(\mathbb E[\widetilde X\mathbf1_{\{|\widetilde X|\le x/2\}}]\right)^2
$$
となる。残りの $|D|\le x$ の領域の積分は $2x^2\mathbb P(|\widetilde X|>x/2)=o(V_D(x))$ 以下である。第一の係数は $2$ に向かい、中心化した一次積率の切断は零へ向かう。$V_D$ は最終的に正の定数以上なので、$W(x/2)\sim V_D(x)/2$ である。従って $W$ も緩変動である。

二次積率が有限なら $V(x)\to\mathbb EX^2>0$ で結論は自明である。無限なら $W(x)\to\infty$ である。定数 $m=\mathbb EX$ に対し $|X^2-\widetilde X^2|\le2|m||\widetilde X|+m^2$ は可積分であり、$|X|\le x$ の範囲は $\{|\widetilde X|\le x-|m|\}$ と $\{|\widetilde X|\le x+|m|\}$ の間にある。よって $V(x)$ は $W(x-|m|)-O(1)$ と $W(x+|m|)+O(1)$ で挟まれる。$W$ は非減少で緩変動なので固定した定数のずれでは $W(x\pm|m|)/W(x)\to1$ となる。従って $V(x)\sim W(x)$ となり、緩変動性が従う。
:::

## 演習問題

::::exercise{level=basic}
$N(6,12)$ の四乗根を求めよ。

:::solution
平均と分散の加法性で $N(3/2,3)$ の四重畳み込みは $N(6,12)$ になる。
:::
::::

::::exercise{level=basic}
$\operatorname{Pois}(8)$ の四乗根を求めよ。

:::solution
特性指数 $8(e^{it}-1)$ を四で割ると $\operatorname{Pois}(2)$ の指数になる。
:::
::::

::::exercise{level=basic}
強度 $2$ で常に跳躍が $3$ の複合 Poisson 変数の平均と分散を求めよ。

:::solution
変数は $3N$ で $N\sim\operatorname{Pois}(2)$ である。従って平均は $6$、分散は $9\cdot2=18$ である。
:::
::::

::::exercise{level=basic}
$\phi(t)=e^{-2|t|}$ に対応する密度を求めよ。

:::solution
反転で $(1/\pi)\int_0^\infty e^{-2t}\cos(tx)\,\dd t$ となる。積分値は $2/(4+x^2)$ なので密度は $2/(\pi(4+x^2))$ である。
:::
::::

::::exercise{level=advanced}
$\nu(\dd x)=\mathbf1_{\{x>0\}}x^{-3/2}\dd x$ が Lévy 測度であることを示せ。

:::solution
零から一までは二乗を掛けた $\int_0^1x^{1/2}\dd x=2/3$、一から先は $\int_1^\infty x^{-3/2}\dd x=2$ であり、ともに有限である。
:::
::::

::::exercise{level=advanced}
$\mu$ が無限分解可能なら、任意の $a\in\R$ に対し $aX$ の分布もそうなることを示せ。

:::solution
$X$ を $n$ 個の独立同分布変数の和として実現し、全項を $a$ 倍する。独立同分布性は可測な同じ関数の適用で保たれ、和は $aX$ の分布になる。
:::
::::

::::exercise{level=advanced}
各項が確率 $2/n$ で $-3$、それ以外で零を取る $n$ 項の独立和の極限を求めよ。

:::solution
特性関数は $(1+(2/n)(e^{-3it}-1))^n$ であり、極限は $\exp(2(e^{-3it}-1))$ である。従って $N\sim\operatorname{Pois}(2)$ に対する $-3N$ の分布である。
:::
::::

::::exercise{level=advanced}
$\mathbb P(X>x)\sim2x^{-3/2}$、$\mathbb P(X<-x)\sim x^{-3/2}$ のとき、前節の極限の $a_n,c_+,c_-$ を求めよ。

:::solution
$\alpha=3/2$ を代入すると $a_n=n^{2/3}$、$c_+=(3/2)\cdot2=3$、$c_-=(3/2)\cdot1=3/2$ である。
:::
::::

## 確認テスト

::::quiz
無限分解可能性の意味はどれか。

- [x] 任意の正整数について畳み込みの根をもつ。
- [ ] 和が常に有界な独立変数で書ける。
- [ ] すべての積率が有限である。

:::explanation
根の存在を項数ごとに要求する。特に Cauchy 分布のように一次絶対積率のない例も含まれる。
:::
::::

::::quiz{answer="16" tol="1e-6"}
強度 $4$ で跳躍が常に $2$ の複合 Poisson 分布の分散を求めよ。

:::explanation
$N\sim\operatorname{Pois}(4)$ による $2N$ なので、分散は $2^2\cdot4=16$ である。
:::
::::

::::quiz{answer="2" tol="1e-6"}
$N(0,10)$ の五乗根の分散を求めよ。

:::explanation
独立和の正規分散は足される。五つの同じ分散の和が $10$ なので根の分散は $10/5=2$ である。
:::
::::

::::quiz
Lévy 測度に対して零付近で課す条件はどれか。

- [x] 二乗を掛けた積分が有限。
- [ ] 測度自体の全質量が有限。
- [ ] 密度が原点で有界。

:::explanation
二次の余りを積分するため $\int_{|x|\le1}x^2\,\nu(\dd x)<\infty$ を課す。密度 $|x|^{-3/2}$ のように全質量が無限でも満たせる。
:::
::::

::::quiz{answer="12" tol="1e-6"}
$N\sim\operatorname{Pois}(3)$、$X=2N$ とする。$X$ の分散を答えよ。

:::explanation
Poisson 変数は分散が率と同じなので $\operatorname{Var}(N)=3$ である。定数倍では分散を倍率の二乗倍するから $\operatorname{Var}(X)=2^2\cdot3=12$ である。跳躍測度を使えば $\int x^2\,3\delta_2(\dd x)=12$ と計算できる。
:::
::::

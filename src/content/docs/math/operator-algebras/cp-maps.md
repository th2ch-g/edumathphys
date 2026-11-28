---
title: 完全正写像
description: 行列レベルの正性から膨張と拡張を導き、条件付き期待値を特徴付ける。
notation:
  - symbol: '\Phi^{(n)}'
    meaning: Φ を各成分に適用した写像
  - symbol: 'C_\Phi'
    meaning: Φ の Choi 行列
  - symbol: 'E:A\to B'
    meaning: 部分環 B への条件付き期待値
---

## この章の目標

- 正値性と完全正値性の違いを具体的に判定する。
- Stinespring の膨張と Choi 行列から写像を構成する。
- Arveson の拡張を正汎関数の延長に帰着する。
- 条件付き期待値の圧縮と両側加群性を使う。

## 作用素系と完全正写像

:::definition[作用素系と行列順序]{#def-oa-complete-positivity}
単位的 C*-環 $A$ の $1$ を含む自己共役線形部分空間 $S$ を :term[作用素系]{en="operator system" yomi="さようそけい"} という。$M_n(S)$ の正錐は $M_n(A)_+$ との共通部分とする。作用素系 $S,T$ の間の線形写像 $\Phi:S\to T$ に対し $\Phi^{(n)}([a_{ij}])=[\Phi(a_{ij})]$ と置く。全 $n\ge1$ で正錐を保つものを :term[完全正写像]{en="completely positive map" yomi="かんぜんせいしゃぞう"}（CP）という。単位を保てば UCP と略記する。
:::

:::example[転置の失敗]
転置 $a\mapsto a^{\mathsf T}$ は正行列の共役に一致するので正値である。だが $M_2(M_2)$ の正な階数一行列 $[E_{ij}]$ は転置増幅で $[E_{ji}]$ になる。これはテンソルの二因子を交換する作用素で、$e_1\otimes e_2-e_2\otimes e_1$ に固有値 $-1$ を与える。従って転置は $2$ 正ですらない。
:::

$T^*\pi(a)T$ という形は CP である。各増幅で *-準同型と $1_n\otimes T$ による共役を行うためである。この構成が全 CP 写像を与える。

作用素系の行列順序から行列ノルムも復元できる。$x\in M_n(S)$ と $r>0$ に対し
$$
\begin{pmatrix}r1_n&x\\x^*&r1_n\end{pmatrix}\ge0
\quad\Longleftrightarrow\quad \norm x\le r
$$
である。忠実表現でブロック行列の二次形式を平方完成すると、残る条件は $r1_n-r^{-1}x^*x\ge0$ となる。従って UCP 写像は全行列段で縮小的であり、行列順序を逆写像も保つ単位的な同型は全行列段で等長である。このノルム判定は $x^*x$ が $M_n(S)$ に属するとは限らない作用素系でも使える。

### k 正性の階層

$\Phi$ の $k$ 次増幅が正であるとき $k$ 正という。$l<k$ の正行列を零ブロックと直和して $k$ 次に入れ、像の左上角を取れば $l$ 正性も従う。一方 $k$ 正性からすべての増幅を制御できるとは限らない。

:::proposition[還元型の写像の正性]
$1\le k<n$ とする。$R_k:M_n\to M_n$、$R_k(a)=k\tr(a)I_n-a$ は $k$ 正だが $(k+1)$ 正ではない。
:::

:::proof
任意の正行列は階数一正元の和なので、$yy^*$ だけを調べればよい。$y\in\C^k\otimes\C^n$ の係数行列を特異値分解すると、二つの因子の正規直交系により $y=\sum_{j=1}^{r}\sqrt{\lambda_j}e_j\otimes f_j$ と書ける。ここで $r\le k$、$\lambda_j>0$ である。$G=\sum_j\lambda_jE_{jj}$、$P$ をその台射影、$\Omega_r=\sum_{j=1}^re_j\otimes f_j$ とすれば、
$$
(\id_k\otimes R_k)(yy^*)
=(G^{1/2}\otimes I)(kP\otimes I-\Omega_r\Omega_r^*)(G^{1/2}\otimes I)
$$
となる。中央の括弧は $P\otimes I$ の像では $kI$ から固有値 $r$ の正の階数一作用素を引く形で、その直交補では零である。従って非負である。これに対し $\Omega_{k+1}$ を試せば、$G=I_{k+1}$ で $\Omega_{k+1}$ の向きの固有値は $k-(k+1)=-1$ になる。
:::

### Schwarz 不等式を二正性だけから導く

:::proposition[単位的二正写像の Schwarz 不等式]{#prop-oa-two-positive-schwarz}
単位的二正写像 $\Phi:A\to B(H)$ は $\Phi(a)^*\Phi(a)\le\Phi(a^*a)$ を満たす。
:::

:::proof
$\begin{pmatrix}a^*a&a^*\\a&1\end{pmatrix}=\begin{pmatrix}a^*\\1\end{pmatrix}\begin{pmatrix}a&1\end{pmatrix}$ は $M_2(A)$ の正元である。その像は下右が $I$ の正行列になる。正写像が対合を保つことは、自己共役元を二つの正元の差に分ければ従う。従って像の二次形式に $(\xi,-\Phi(a)\xi)$ を代入すれば
$$
\inner{(\Phi(a^*a)-\Phi(a)^*\Phi(a))\xi}{\xi}\ge0
$$
となる。全ての $\xi$ で成立するため主張を得る。
:::

:::example[正写像でも Schwarz 不等式が破れる]
$M_2$ の転置 $\Phi$ と $a=E_{12}$ を取る。転置は正行列の共役であり、正性を保つ。一方 $\Phi(a^*a)=E_{22}$、$\Phi(a)^*\Phi(a)=E_{11}$ である。差 $E_{22}-E_{11}$ の第 1 座標の値は $-1$ なので、Schwarz 不等式を満たさない。単位的な正写像でも二正性までは自動的には従わないことを、一つの行列単位で検出できる。
:::

:::example[二正だが完全正でない写像の Schwarz 差]
$M_3$ の $R_2(a)=2\tr(a)I-a$ を $5$ で割って $\Phi$ とする。$\Phi(1)=I$ で、二正性は還元型の写像の命題による。$a=E_{12}$ なら $\Phi(a)=-E_{12}/5$、$\Phi(a^*a)=(2I-E_{22})/5$ である。差は $\diag(2/5,4/25,2/5)$ で正である。一方 $3$ 次の増幅で $\Omega_3\Omega_3^*$ を試すと $\Omega_3$ に沿って固有値 $-1/5$ が現れる。従って一つの Schwarz 不等式は完全正性を判定できない。
:::

:::example[正だが完全正でない二次行列の写像]
$M_d$ 上の写像 $\Phi_p(x)=p x+(1-p)\tr(x)I/d$ を調べる。$p$ が負なら二つの完全正写像の凸結合という説明は使えないので、階数一の入力で各増幅を検査する。

$k$ 次増幅では任意の正行列が階数一正行列の和である。Schmidt 階数 $r\le\min(k,d)$ のベクトルから作る入力を考えると、還元型の写像と同様に、合同の間に残る行列は $(1-p)(P\otimes I)/d+p\Omega_r\Omega_r^*$ である。

$$
\begin{aligned}
(d,k,p)&=\left( 2, \  1, \  - \frac{1}{2}\right)\\
(1-p)/d&=\frac{3}{4}
\end{aligned}
$$

中央行列は台の中で $\Omega_r$ とその直交補に分かれる。$\Omega_r$ の二乗ノルムは $r$ なので、固有値は $(1-p)/d+pr$ と $(1-p)/d$ になる。$p<0$ では前者が最小になる最大の $r$ を試せばよい。

$$
\begin{aligned}
r_{\max}&=1\\
(1-p)/d+p r_{\max}&=\frac{1}{4}
\end{aligned}
$$

完全正性の検査では $r=d$ を取る。これは $d$ 個の標準基底を対応させた最大の Schmidt 階数であり、その入力の像は Choi 行列になる。

$$
\begin{aligned}
(1-p)/d+pd&=- \frac{1}{4}\\
(1-p)/d&=\frac{3}{4}
\end{aligned}
$$

ここでは指定した増幅までは非負だが、Choi 行列に負の固有値があるので完全正ではない。
:::

:::example[二次行列の完全正性の境界]
$M_d$ 上の写像 $\Phi_p(x)=p x+(1-p)\tr(x)I/d$ を調べる。$p$ が負なら二つの完全正写像の凸結合という説明は使えないので、階数一の入力で各増幅を検査する。

$k$ 次増幅では任意の正行列が階数一正行列の和である。Schmidt 階数 $r\le\min(k,d)$ のベクトルから作る入力を考えると、還元型の写像と同様に、合同の間に残る行列は $(1-p)(P\otimes I)/d+p\Omega_r\Omega_r^*$ である。

$$
\begin{aligned}
(d,k,p)&=\left( 2, \  2, \  - \frac{1}{3}\right)\\
(1-p)/d&=\frac{2}{3}
\end{aligned}
$$

中央行列は台の中で $\Omega_r$ とその直交補に分かれる。$\Omega_r$ の二乗ノルムは $r$ なので、固有値は $(1-p)/d+pr$ と $(1-p)/d$ になる。$p<0$ では前者が最小になる最大の $r$ を試せばよい。

$$
\begin{aligned}
r_{\max}&=2\\
(1-p)/d+p r_{\max}&=0
\end{aligned}
$$

完全正性の検査では $r=d$ を取る。これは $d$ 個の標準基底を対応させた最大の Schmidt 階数であり、その入力の像は Choi 行列になる。

$$
\begin{aligned}
(1-p)/d+pd&=0\\
(1-p)/d&=\frac{2}{3}
\end{aligned}
$$

このパラメータでは両方の固有値が非負なので完全正である。
:::

:::example[二正だが三正でない三次行列の写像]
$M_d$ 上の写像 $\Phi_p(x)=p x+(1-p)\tr(x)I/d$ を調べる。$p$ が負なら二つの完全正写像の凸結合という説明は使えないので、階数一の入力で各増幅を検査する。

$k$ 次増幅では任意の正行列が階数一正行列の和である。Schmidt 階数 $r\le\min(k,d)$ のベクトルから作る入力を考えると、還元型の写像と同様に、合同の間に残る行列は $(1-p)(P\otimes I)/d+p\Omega_r\Omega_r^*$ である。

$$
\begin{aligned}
(d,k,p)&=\left( 3, \  2, \  - \frac{1}{6}\right)\\
(1-p)/d&=\frac{7}{18}
\end{aligned}
$$

中央行列は台の中で $\Omega_r$ とその直交補に分かれる。$\Omega_r$ の二乗ノルムは $r$ なので、固有値は $(1-p)/d+pr$ と $(1-p)/d$ になる。$p<0$ では前者が最小になる最大の $r$ を試せばよい。

$$
\begin{aligned}
r_{\max}&=2\\
(1-p)/d+p r_{\max}&=\frac{1}{18}
\end{aligned}
$$

完全正性の検査では $r=d$ を取る。これは $d$ 個の標準基底を対応させた最大の Schmidt 階数であり、その入力の像は Choi 行列になる。

$$
\begin{aligned}
(1-p)/d+pd&=- \frac{1}{9}\\
(1-p)/d&=\frac{7}{18}
\end{aligned}
$$

ここでは指定した増幅までは非負だが、Choi 行列に負の固有値があるので完全正ではない。
:::

:::example[三次行列の完全正性の境界]
$M_d$ 上の写像 $\Phi_p(x)=p x+(1-p)\tr(x)I/d$ を調べる。$p$ が負なら二つの完全正写像の凸結合という説明は使えないので、階数一の入力で各増幅を検査する。

$k$ 次増幅では任意の正行列が階数一正行列の和である。Schmidt 階数 $r\le\min(k,d)$ のベクトルから作る入力を考えると、還元型の写像と同様に、合同の間に残る行列は $(1-p)(P\otimes I)/d+p\Omega_r\Omega_r^*$ である。

$$
\begin{aligned}
(d,k,p)&=\left( 3, \  3, \  - \frac{1}{8}\right)\\
(1-p)/d&=\frac{3}{8}
\end{aligned}
$$

中央行列は台の中で $\Omega_r$ とその直交補に分かれる。$\Omega_r$ の二乗ノルムは $r$ なので、固有値は $(1-p)/d+pr$ と $(1-p)/d$ になる。$p<0$ では前者が最小になる最大の $r$ を試せばよい。

$$
\begin{aligned}
r_{\max}&=3\\
(1-p)/d+p r_{\max}&=0
\end{aligned}
$$

完全正性の検査では $r=d$ を取る。これは $d$ 個の標準基底を対応させた最大の Schmidt 階数であり、その入力の像は Choi 行列になる。

$$
\begin{aligned}
(1-p)/d+pd&=0\\
(1-p)/d&=\frac{3}{8}
\end{aligned}
$$

このパラメータでは両方の固有値が非負なので完全正である。
:::

:::example[恒等写像と完全混合の平均]
$M_d$ 上の写像 $\Phi_p(x)=p x+(1-p)\tr(x)I/d$ を調べる。$p$ が負なら二つの完全正写像の凸結合という説明は使えないので、階数一の入力で各増幅を検査する。

$k$ 次増幅では任意の正行列が階数一正行列の和である。Schmidt 階数 $r\le\min(k,d)$ のベクトルから作る入力を考えると、還元型の写像と同様に、合同の間に残る行列は $(1-p)(P\otimes I)/d+p\Omega_r\Omega_r^*$ である。

$$
\begin{aligned}
(d,k,p)&=\left( 2, \  2, \  \frac{1}{2}\right)\\
(1-p)/d&=\frac{1}{4}
\end{aligned}
$$

中央行列は台の中で $\Omega_r$ とその直交補に分かれる。$\Omega_r$ の二乗ノルムは $r$ なので、固有値は $(1-p)/d+pr$ と $(1-p)/d$ になる。$p<0$ では前者が最小になる最大の $r$ を試せばよい。

$$
\begin{aligned}
r_{\max}&=2\\
(1-p)/d+p r_{\max}&=\frac{5}{4}
\end{aligned}
$$

完全正性の検査では $r=d$ を取る。これは $d$ 個の標準基底を対応させた最大の Schmidt 階数であり、その入力の像は Choi 行列になる。

$$
\begin{aligned}
(1-p)/d+pd&=\frac{5}{4}\\
(1-p)/d&=\frac{1}{4}
\end{aligned}
$$

このパラメータでは両方の固有値が非負なので完全正である。
:::

## Stinespring の膨張定理

:::theorem[Stinespring の膨張]{#thm-oa-stinespring}
単位的 C*-環 $A$ と Hilbert 空間 $H$ に対し、$\Phi:A\to B(H)$ が CP なら、ある単位的 *-表現 $\pi$ と有界作用素 $V:H\to K$ により $\Phi(a)=V^*\pi(a)V$ と書ける。$\overline{\spn\pi(A)VH}=K$ を満たす模型はユニタリ同値を除き一意である。
:::
:::proof
記号 $a\otimes\xi$ の線形包を双線形関係で割った代数的テンソル積を $A\odot H$ と書く。その上で
$$
\inner{\sum_i a_i\otimes\xi_i}{\sum_j b_j\otimes\eta_j}
=\sum_{i,j}\inner{\Phi(b_j^*a_i)\xi_i}{\eta_j}
$$
とする。$[a_j^*a_i]_{j,i}$ は行列積で正なので、完全正値性によりこの形式は半正定値である。零空間で割り完備化する。左乗法の有界性は $[a_j^*c^*ca_i]_{j,i}\le\norm c^2[a_j^*a_i]_{j,i}$ を増幅した不等式から従う。積は左乗法で保たれ、随伴は形式の定義で保たれる。$V\xi=[1\otimes\xi]$ とすると $\norm{V\xi}^2=\inner{\Phi(1)\xi}{\xi}$ であり、$\Phi(a)=V^*\pi(a)V$ が全行列係数で成立する。生成した空間は定義により稠密である。他の最小模型にも $\sum_i\pi(a_i)V\xi_i$ の内積が同じ式で計算されるため、対応は表現と $V$ を保つユニタリに延びる。
:::

$\norm{\Phi^{(n)}}=\norm{\Phi(1)}$ が全ての $n$ で成り立つ。膨張を $I_n\otimes V$ による圧縮へ増幅すると上界は $\norm V^2$ であり、$I_n\otimes1$ に代入して等号となる。従って増幅のノルムの上限である完全有界ノルム $\norm{\Phi}_{\rm cb}$ も $\norm{\Phi(1)}$ に等しい。UCP なら $V$ は等長であり、
$$
\Phi(a^*a)-\Phi(a)^*\Phi(a)=V^*\pi(a)^*(1-VV^*)\pi(a)V\ge0
$$
である。$a^*a$ と $aa^*$ に関する二つの等号が成立すると、$\pi(a)$ と $\pi(a)^*$ が $VH$ を保つ。従って圧縮を挟んでも $\Phi(ab)=\Phi(a)\Phi(b)$、$\Phi(ba)=\Phi(b)\Phi(a)$ となる。このような元の全体は $VH$ を還元する元の逆像なので部分環であり、乗法領域と呼ぶ。

:::example[端の座標への圧縮]
$\Phi:M_2\to\C$、$\Phi(a)=a_{11}$ は $V1=e_1$ による圧縮で UCP である。$a=E_{21}$ に対し $\Phi(a^*a)=1$、$\Phi(a)^*\Phi(a)=0$ であり、等号に必要な不変性は破れる。
:::

:::proposition[正規な完全正写像の最小膨張]{#prop-oa-normal-stinespring}
$M$ が von Neumann 環なら、$\Phi:M\to B(H)$ の正規性とその最小 Stinespring 表現 $\pi$ の正規性は同値である。
:::

:::proof
$0\le x_i\uparrow x$ を有界増大ネットとする。$\Phi$ が正規なら、有限和 $u=\sum_j\pi(a_j)V\xi_j$ で
$$
\inner{\pi(x_i)u}{u}=\sum_{j,k}\inner{\Phi(a_k^*x_i a_j)\xi_j}{\xi_k}
$$
は $x_i$ を $x$ に替えた値へ収束する。実際 $x_i\to x$ は σ 弱収束であり、固定元による左右の積と $\Phi$ は σ 弱連続である。$\pi(x_i)$ は一様有界なので、最小性により全 $u$ へ延長できる。従って $\pi(x_i)\uparrow\pi(x)$ である。逆に $\pi$ がこの上限を保てば、$V$ による圧縮も行列係数で上限を保つ。
:::

:::corollary[単位元のない場合の Stinespring 表示]{#cor-oa-nonunital-stinespring}
非単位的 C*-環 $A$ の完全正写像 $\Phi:A\to B(H)$ にも、非退化表現 $\pi$ と有界作用素 $V$ による $\Phi(a)=V^*\pi(a)V$ がある。
:::

:::proof
まず $\Phi$ は有界である。正縮小元上で有界でなければ、$\norm{\Phi(a_n)}\ge n2^n$ となる正縮小元 $a_n$ を選べる。$a=\sum_n 2^{-n}a_n$ は $a\ge2^{-n}a_n$ なので、正値性により $\norm{\Phi(a)}\ge n$ が全 $n$ で成り立ってしまう。自己共役元の正負分解と一般元の自己共役分解で有界性は全空間に延びる。

:ref[thm-oa-stinespring] と同じ半内積で $A\odot H$ を完備化する。不等式 $\norm{[a\otimes\xi]}\le\norm{\Phi}^{1/2}\norm a\norm\xi$ により、近似単位元 $(e_\lambda)$ の左乗法は生成ベクトル上で恒等に強収束し、表現は非退化である。固定した $\xi$ に対して $\varphi_\xi(a)=\inner{\Phi(a)\xi}{\xi}$ は正値汎関数である。ネット $[e_\lambda\otimes\xi]$ の二項間の距離の二乗は $\varphi_\xi((e_\lambda-e_\mu)^2)$ である。この距離は $\varphi_\xi$ を規格化し単位化した GNS 空間での $[e_\lambda]-[e_\mu]$ の距離の定数倍になる。両方とも $[1]$ に近づくので Cauchy であり、零汎関数の場合は各項が零である。$V\xi=\lim_\lambda[e_\lambda\otimes\xi]$ は線形で $\norm V\le\norm{\Phi}^{1/2}$ である。$ae_\lambda\to a$ により $\pi(a)V\xi=[a\otimes\xi]$ である。これと $V\eta$ の内積を取ると $\inner{\Phi(e_\lambda a)\xi}{\eta}$ の極限となり、所要の式を得る。
:::

### 乗法域

単位的な完全正写像 $\Phi(a)=V^*\pi(a)V$ に対して $V$ は等長である。$Q=VV^*$ とすると
$$
\Phi(a^*a)-\Phi(a)^*\Phi(a)=V^*\pi(a)^*(1-Q)\pi(a)V\ge0
$$
であり、これを Schwarz 不等式と呼ぶ。等号なら $(1-Q)\pi(a)V=0$ であり、$a^*$ にも等号を仮定すれば $\pi(a)$ が $VH$ を還元する。逆向きも式から分かる。従って二つの等号が成り立つ集合は、写像 $\pi$ で $Q$ と可換になる元の逆像なので C*-部分環である。これを :term[乗法域]{en="multiplicative domain" yomi="じょうほういき"} という。その元 $a$ に対して $\Phi(ab)=\Phi(a)\Phi(b)$、$\Phi(ba)=\Phi(b)\Phi(a)$ が成り立つ。実際、積の途中に $Q$ を挿入できるからである。

:::example[ユニタリの像で乗法域を検査する]
単位的完全正写像がユニタリ $u$ をユニタリへ写せば、$\Phi(u)^*\Phi(u)=\Phi(u^*u)=1$ と逆順の等式から $u$ は乗法域に属する。例えば対角期待値は対角ユニタリ $\diag(1,-1)$ を固定するので、これを任意の元に左右から掛ける操作と交換する。一方 $E_{12}+E_{21}$ はユニタリだが期待値で零になるため乗法域に属さない。
:::

### CP の大小を膨張空間で表す

:::proposition[完全正写像の Radon–Nikodym 表示]{#prop-oa-cp-radon-nikodym}
単位的 C*-環 $A$ から $B(H)$ への完全正写像 $\Phi$ を取り、$(\pi,K,V)$ を最小 Stinespring 表示とする。$\Psi$ と $\Phi-\Psi$ がともに完全正であることは、ある一意な $0\le T\le I_K$ が $\pi(A)'$ にあって $\Psi(a)=V^*T\pi(a)V$ となることと同値である。
:::

:::proof
最小性により $\sum_i\pi(a_i)V\xi_i$ の線形空間は $K$ に稠密である。この形の二つのベクトル $u,v$ 上に
$$
B_\Psi(u,v)=\sum_{i,j}\inner{\Psi(b_j^*a_i)\xi_i}{\eta_j},
\qquad v=\sum_j\pi(b_j)V\eta_j
$$
を定める。まず代数的テンソルの表示で評価すると、完全正性から対角値は非負である。また $\Phi-\Psi$ にもこれを使うと $0\le B_\Psi(u,u)\le\norm u^2$ となる。正半定値形式の Cauchy–Schwarz によりノルム零の表示は他の任意の表示との組合せで零を与える。従って商に降りてよく定義され、$|B_\Psi(u,v)|\le\norm u\norm v$ で $K$ 全体に延長する。

Riesz 表現により $B_\Psi(u,v)=\inner{Tu}{v}$ の正縮小元 $T$ が存在する。任意の $x\in A$ で $B_\Psi(\pi(x)u,v)=B_\Psi(u,\pi(x^*)v)$ が成り立つ。実際、両辺の式は同じ $b_j^*xa_i$ を用いる。従って $T\pi(x)=\pi(x)T$ である。$a_1=a$、$b_1=1$ の表示では $V^*T\pi(a)V=\Psi(a)$ となる。逆に可換子の正縮小元 $T$ があれば、$T^{1/2}V$ による圧縮で $\Psi$ が、$(I-T)^{1/2}V$ による圧縮で $\Phi-\Psi$ が得られ、両方とも完全正である。

最後に $T$ の候補は $\Psi(b^*a)$ から全ての稠密ベクトル $\pi(a)V\xi,\pi(b)V\eta$ 上の行列係数が決まる。よって一意である。
:::

:::example[部分測定の重みを可換子で選ぶ]
$\Phi:M_2\to M_2$ を $\Phi(a)=(a+ZaZ)/2$、$Z=\diag(1,-1)$ とする。最小表示は $\pi(a)=a\oplus a$、$V\xi=(\xi,Z\xi)/\sqrt2$ である。例えば $T=\diag(I/3,2I/3)$ は可換子の正縮小元なので、$\Psi(a)=a/6+ZaZ/3$ は $\Phi$ に CP 順序で支配される。実際 $\Phi(a)-\Psi(a)=a/3+ZaZ/6$ であり、両方とも正の重みの和である。$T$ を可換子の非対角方向に選べば、二つの Kraus 作用素の交差項も使える。例えば $T=\tfrac12\begin{pmatrix}I&I\\I&I\end{pmatrix}$ は対角方向の射影で、このとき $\Psi(a)=(1+Z)a(1+Z)/4=E_{11}aE_{11}$ となる。元の全測定 $\Phi$ から一つの座標だけを取る操作がこの $T$ で実現される。
:::

:::example[コヒーレンスを完全に消す]
$d=\diag(1,-1)$、$p=(1+k)/2$、$q=1-p$ とし、$\Phi(a)=pa+qdad$ を考える。有限次元の表示で完全正性と非対角成分の変化を同時に調べる。

代数の表現を $\pi(a)=a\oplus a$ とし、$V\xi=(\sqrt p\,\xi,\sqrt q\,d\xi)$ と取る。$V^*V=(p+q)I=I$ なので、これは等長による圧縮になる。

$$
\begin{aligned}
(k,p,q)&=\left( 0, \  \frac{1}{2}, \  \frac{1}{2}\right)\\
V&=\begin{pmatrix}\frac{\sqrt{2}}{2} & 0\\0 & \frac{\sqrt{2}}{2}\\\frac{\sqrt{2}}{2} & 0\\0 & - \frac{\sqrt{2}}{2}\end{pmatrix}\\
a&=\begin{pmatrix}2 & 1 + i\\1 - i & 3\end{pmatrix}
\end{aligned}
$$

$V^*\pi(a)V$ の二項は、同じブロックからの寄与 $pa$ と $qdad$ である。$dad$ では対角は同じで、非対角だけ符号が変わる。従って非対角には $p-q=k$ が掛かる。

$$
\begin{aligned}
dad&=\begin{pmatrix}2 & -1 - i\\-1 + i & 3\end{pmatrix}\\
\Phi(a)&=\begin{pmatrix}2 & 0\\0 & 3\end{pmatrix}
\end{aligned}
$$

この表示を行列増幅しても $I_m\otimes V$ を使った圧縮になるため正性を保つ。端点 $p=0$ または $q=0$ では対応するブロックが不要になる。それ以外では $d$ と $I$ は一次独立なので、行列単位を適用して生成される空間は両ブロック全体である。

$$
\begin{aligned}
\dim H_{\min}&=4
\end{aligned}
$$

この場合、対角の確率を全て保っても非対角の内積は変化し得る。写像が状態を変えるかどうかを対角だけから判定できない。
:::

:::example[コヒーレンスを半分にする]
$d=\diag(1,-1)$、$p=(1+k)/2$、$q=1-p$ とし、$\Phi(a)=pa+qdad$ を考える。有限次元の表示で完全正性と非対角成分の変化を同時に調べる。

代数の表現を $\pi(a)=a\oplus a$ とし、$V\xi=(\sqrt p\,\xi,\sqrt q\,d\xi)$ と取る。$V^*V=(p+q)I=I$ なので、これは等長による圧縮になる。

$$
\begin{aligned}
(k,p,q)&=\left( \frac{1}{2}, \  \frac{3}{4}, \  \frac{1}{4}\right)\\
V&=\begin{pmatrix}\frac{\sqrt{3}}{2} & 0\\0 & \frac{\sqrt{3}}{2}\\\frac{1}{2} & 0\\0 & - \frac{1}{2}\end{pmatrix}\\
a&=\begin{pmatrix}2 & 1 + i\\1 - i & 3\end{pmatrix}
\end{aligned}
$$

$V^*\pi(a)V$ の二項は、同じブロックからの寄与 $pa$ と $qdad$ である。$dad$ では対角は同じで、非対角だけ符号が変わる。従って非対角には $p-q=k$ が掛かる。

$$
\begin{aligned}
dad&=\begin{pmatrix}2 & -1 - i\\-1 + i & 3\end{pmatrix}\\
\Phi(a)&=\begin{pmatrix}2 & \frac{1}{2} + \frac{i}{2}\\\frac{1}{2} - \frac{i}{2} & 3\end{pmatrix}
\end{aligned}
$$

この表示を行列増幅しても $I_m\otimes V$ を使った圧縮になるため正性を保つ。端点 $p=0$ または $q=0$ では対応するブロックが不要になる。それ以外では $d$ と $I$ は一次独立なので、行列単位を適用して生成される空間は両ブロック全体である。

$$
\begin{aligned}
\dim H_{\min}&=4
\end{aligned}
$$

この場合、対角の確率を全て保っても非対角の内積は変化し得る。写像が状態を変えるかどうかを対角だけから判定できない。
:::

:::example[コヒーレンスの符号を反転する]
$d=\diag(1,-1)$、$p=(1+k)/2$、$q=1-p$ とし、$\Phi(a)=pa+qdad$ を考える。有限次元の表示で完全正性と非対角成分の変化を同時に調べる。

代数の表現を $\pi(a)=a\oplus a$ とし、$V\xi=(\sqrt p\,\xi,\sqrt q\,d\xi)$ と取る。$V^*V=(p+q)I=I$ なので、これは等長による圧縮になる。

$$
\begin{aligned}
(k,p,q)&=\left( - \frac{1}{2}, \  \frac{1}{4}, \  \frac{3}{4}\right)\\
V&=\begin{pmatrix}\frac{1}{2} & 0\\0 & \frac{1}{2}\\\frac{\sqrt{3}}{2} & 0\\0 & - \frac{\sqrt{3}}{2}\end{pmatrix}\\
a&=\begin{pmatrix}2 & 1 + i\\1 - i & 3\end{pmatrix}
\end{aligned}
$$

$V^*\pi(a)V$ の二項は、同じブロックからの寄与 $pa$ と $qdad$ である。$dad$ では対角は同じで、非対角だけ符号が変わる。従って非対角には $p-q=k$ が掛かる。

$$
\begin{aligned}
dad&=\begin{pmatrix}2 & -1 - i\\-1 + i & 3\end{pmatrix}\\
\Phi(a)&=\begin{pmatrix}2 & - \frac{1}{2} - \frac{i}{2}\\- \frac{1}{2} + \frac{i}{2} & 3\end{pmatrix}
\end{aligned}
$$

この表示を行列増幅しても $I_m\otimes V$ を使った圧縮になるため正性を保つ。端点 $p=0$ または $q=0$ では対応するブロックが不要になる。それ以外では $d$ と $I$ は一次独立なので、行列単位を適用して生成される空間は両ブロック全体である。

$$
\begin{aligned}
\dim H_{\min}&=4
\end{aligned}
$$

この場合、対角の確率を全て保っても非対角の内積は変化し得る。写像が状態を変えるかどうかを対角だけから判定できない。
:::

:::example[恒等写像に一致する]
$d=\diag(1,-1)$、$p=(1+k)/2$、$q=1-p$ とし、$\Phi(a)=pa+qdad$ を考える。有限次元の表示で完全正性と非対角成分の変化を同時に調べる。

代数の表現を $\pi(a)=a\oplus a$ とし、$V\xi=(\sqrt p\,\xi,\sqrt q\,d\xi)$ と取る。$V^*V=(p+q)I=I$ なので、これは等長による圧縮になる。

$$
\begin{aligned}
(k,p,q)&=\left( 1, \  1, \  0\right)\\
V&=\begin{pmatrix}1 & 0\\0 & 1\\0 & 0\\0 & 0\end{pmatrix}\\
a&=\begin{pmatrix}2 & 1 + i\\1 - i & 3\end{pmatrix}
\end{aligned}
$$

$V^*\pi(a)V$ の二項は、同じブロックからの寄与 $pa$ と $qdad$ である。$dad$ では対角は同じで、非対角だけ符号が変わる。従って非対角には $p-q=k$ が掛かる。

$$
\begin{aligned}
dad&=\begin{pmatrix}2 & -1 - i\\-1 + i & 3\end{pmatrix}\\
\Phi(a)&=\begin{pmatrix}2 & 1 + i\\1 - i & 3\end{pmatrix}
\end{aligned}
$$

この表示を行列増幅しても $I_m\otimes V$ を使った圧縮になるため正性を保つ。端点 $p=0$ または $q=0$ では対応するブロックが不要になる。それ以外では $d$ と $I$ は一次独立なので、行列単位を適用して生成される空間は両ブロック全体である。

$$
\begin{aligned}
\dim H_{\min}&=2
\end{aligned}
$$

この場合、対角の確率を全て保っても非対角の内積は変化し得る。写像が状態を変えるかどうかを対角だけから判定できない。
:::

:::example[対角ユニタリ共役に一致する]
$d=\diag(1,-1)$、$p=(1+k)/2$、$q=1-p$ とし、$\Phi(a)=pa+qdad$ を考える。有限次元の表示で完全正性と非対角成分の変化を同時に調べる。

代数の表現を $\pi(a)=a\oplus a$ とし、$V\xi=(\sqrt p\,\xi,\sqrt q\,d\xi)$ と取る。$V^*V=(p+q)I=I$ なので、これは等長による圧縮になる。

$$
\begin{aligned}
(k,p,q)&=\left( -1, \  0, \  1\right)\\
V&=\begin{pmatrix}0 & 0\\0 & 0\\1 & 0\\0 & -1\end{pmatrix}\\
a&=\begin{pmatrix}2 & 1 + i\\1 - i & 3\end{pmatrix}
\end{aligned}
$$

$V^*\pi(a)V$ の二項は、同じブロックからの寄与 $pa$ と $qdad$ である。$dad$ では対角は同じで、非対角だけ符号が変わる。従って非対角には $p-q=k$ が掛かる。

$$
\begin{aligned}
dad&=\begin{pmatrix}2 & -1 - i\\-1 + i & 3\end{pmatrix}\\
\Phi(a)&=\begin{pmatrix}2 & -1 - i\\-1 + i & 3\end{pmatrix}
\end{aligned}
$$

この表示を行列増幅しても $I_m\otimes V$ を使った圧縮になるため正性を保つ。端点 $p=0$ または $q=0$ では対応するブロックが不要になる。それ以外では $d$ と $I$ は一次独立なので、行列単位を適用して生成される空間は両ブロック全体である。

$$
\begin{aligned}
\dim H_{\min}&=2
\end{aligned}
$$

この場合、対角の確率を全て保っても非対角の内積は変化し得る。写像が状態を変えるかどうかを対角だけから判定できない。
:::

## Choi の定理と Choi 行列

:::theorem[Choi 行列による判定]{#thm-oa-choi}
$\Phi:M_n\to M_m$ について、完全正値性、$n$ 正性、$C_\Phi=[\Phi(E_{ij})]_{i,j=1}^n\ge0$、および有限個の $m\times n$ 行列 $K_r$ で $\Phi(a)=\sum_r K_raK_r^*$ と書けることは同値である。
:::
:::proof
完全正なら定義により $n$ 正である。$\Omega=\sum_i e_i\otimes e_i$ とすると $\Omega\Omega^*=[E_{ij}]$ は正なので、$n$ 正性から Choi 行列は正である。正なら $C_\Phi=\sum_r v_rv_r^*$ と分解し、$v_r=\sum_i e_i\otimes w_{ri}$ と書く。第 $i$ 列が $w_{ri}$ の行列 $K_r$ では $C_\Phi$ の $(i,j)$ ブロックは $\sum_r K_rE_{ij}K_r^*$ である。行列単位は基底なので等式は全 $a$ に延びる。この表示の各項は共役によって全増幅で正性を保つ。
:::

:::example[対角化と恒等写像]
対角成分だけを残す写像は $K_1=E_{11},K_2=E_{22}$ で表されるので CP である。Choi 行列は $\diag(1,0,0,1)$ となり階数 $2$ である。恒等写像の Choi 行列は $\Omega\Omega^*$ で階数 $1$、非零固有値は $2$ である。一般に最小の Kraus 行列の個数は、階数一分解に必要なベクトル数から Choi 行列の階数に一致する。
:::

### Kraus 作用素の取り替え

最小な Kraus 表示と、項がそれより多い表示との関係を Choi 行列から読める。$C_\Phi=FF^*$ と書き、$F$ の列を最小表示の Kraus 行列を縦に並べたベクトルとする。$F$ の列数は $r=\rank C_\Phi$ であり、列は一次独立である。別の表示の列を並べた行列 $G$ は $GG^*=FF^*$ を満たすので、$\Img G=\Img F$ である。従って $G=FW$ と一意に書ける。$F$ の左逆を $FF^*=FWW^*F^*$ の両側に適用すると $WW^*=I_r$ となる。逆にこの式を満たす $W$ を掛ければ同じ Choi 行列を与える。二つの表示がともに最小なら $W$ は正方ユニタリである。

### Kraus 分解の正規化条件

$\Phi(a)=\sum_lK_laK_l^*$（$K_l:\C^n\to\C^m$）なら、単位的である条件は $\sum_lK_lK_l^*=I_m$ である。一方、行列トレースを保つ条件は $\sum_lK_l^*K_l=I_n$ である。後者は $\tr(\Phi(a))=\tr(a\sum_lK_l^*K_l)$ と、行列単位で係数を検査すれば従う。二つの和の行列サイズは一般に異なる。

:::example[二準位から一準位への緩和]
$0\le\gamma\le1$ とし $K_0=\diag(1,\sqrt{1-\gamma})$、$K_1=\sqrt{\gamma}E_{12}$ とする。$\sum K_l^*K_l=I$ だが $\sum K_lK_l^*=\diag(1+\gamma,1-\gamma)$ である。従って $\gamma>0$ のとき写像はトレースを保つが単位的ではない。$\gamma=3/4$ なら $E_{22}$ は $(3/4)E_{11}+(1/4)E_{22}$ に、$E_{12}$ は $(1/2)E_{12}$ に写る。随伴写像 $a\mapsto\sum K_l^*aK_l$ は単位的であり、観測量の側の写像を与える。
:::

:::example[階数一の Choi 行列を持つ場合]
成分ごとの積で定めた Schur 写像 $\Phi(x)=C\circ x$ を考える。Choi 行列を作り、完全正性と必要な Kraus 作用素の個数を求める。

行列単位 $E_{ij}$ の像は $C_{ij}E_{ij}$ である。従って Choi 行列では対応する $e_i\otimes e_i$ の二方向だけに成分がある。

$$
\begin{aligned}
C&=\begin{pmatrix}\frac{1}{2} & \frac{1}{2}\\\frac{1}{2} & \frac{1}{2}\end{pmatrix}\\
J_\Phi&=\begin{pmatrix}\frac{1}{2} & 0 & 0 & \frac{1}{2}\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\\frac{1}{2} & 0 & 0 & \frac{1}{2}\end{pmatrix}
\end{aligned}
$$

この二方向への圧縮は $C$ 自身なので、正性は対角の非負性と行列式の非負性で判定できる。行列式が零なら非零な圧縮の階数は一、正なら二である。

$$
\begin{aligned}
\det C&=0\\
\rank J_\Phi&=1
\end{aligned}
$$

$C=\sum_l v_lv_l^*$ と正行列の固有分解で表せば $K_l=\diag(v_l)$ が Kraus 作用素になる。実際 $(K_lxK_l^*)_{ij}=(v_l)_i\overline{(v_l)_j}x_{ij}$ なので全項を足すと $C_{ij}x_{ij}$ になる。

$$
\begin{aligned}
x&=\begin{pmatrix}1 & 2 + i\\2 - i & 4\end{pmatrix}\\
\Phi(x)&=\begin{pmatrix}\frac{1}{2} & 1 + \frac{i}{2}\\1 - \frac{i}{2} & 2\end{pmatrix}
\end{aligned}
$$

最小個数は Choi 行列の階数である。なお $\Phi(I)$ は $C$ の対角だけを持つので、この例では写像は単位的ではない。
:::

:::example[非対角ブロックが零の場合]
成分ごとの積で定めた Schur 写像 $\Phi(x)=C\circ x$ を考える。Choi 行列を作り、完全正性と必要な Kraus 作用素の個数を求める。

行列単位 $E_{ij}$ の像は $C_{ij}E_{ij}$ である。従って Choi 行列では対応する $e_i\otimes e_i$ の二方向だけに成分がある。

$$
\begin{aligned}
C&=\begin{pmatrix}\frac{1}{2} & 0\\0 & \frac{1}{2}\end{pmatrix}\\
J_\Phi&=\begin{pmatrix}\frac{1}{2} & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & \frac{1}{2}\end{pmatrix}
\end{aligned}
$$

この二方向への圧縮は $C$ 自身なので、正性は対角の非負性と行列式の非負性で判定できる。行列式が零なら非零な圧縮の階数は一、正なら二である。

$$
\begin{aligned}
\det C&=\frac{1}{4}\\
\rank J_\Phi&=2
\end{aligned}
$$

$C=\sum_l v_lv_l^*$ と正行列の固有分解で表せば $K_l=\diag(v_l)$ が Kraus 作用素になる。実際 $(K_lxK_l^*)_{ij}=(v_l)_i\overline{(v_l)_j}x_{ij}$ なので全項を足すと $C_{ij}x_{ij}$ になる。

$$
\begin{aligned}
x&=\begin{pmatrix}1 & 2 + i\\2 - i & 4\end{pmatrix}\\
\Phi(x)&=\begin{pmatrix}\frac{1}{2} & 0\\0 & 2\end{pmatrix}
\end{aligned}
$$

最小個数は Choi 行列の階数である。なお $\Phi(I)$ は $C$ の対角だけを持つので、この例では写像は単位的ではない。
:::

:::example[正だが境界にある非対角係数]
成分ごとの積で定めた Schur 写像 $\Phi(x)=C\circ x$ を考える。Choi 行列を作り、完全正性と必要な Kraus 作用素の個数を求める。

行列単位 $E_{ij}$ の像は $C_{ij}E_{ij}$ である。従って Choi 行列では対応する $e_i\otimes e_i$ の二方向だけに成分がある。

$$
\begin{aligned}
C&=\begin{pmatrix}\frac{2}{3} & \frac{\sqrt{2}}{3}\\\frac{\sqrt{2}}{3} & \frac{1}{3}\end{pmatrix}\\
J_\Phi&=\begin{pmatrix}\frac{2}{3} & 0 & 0 & \frac{\sqrt{2}}{3}\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\\frac{\sqrt{2}}{3} & 0 & 0 & \frac{1}{3}\end{pmatrix}
\end{aligned}
$$

この二方向への圧縮は $C$ 自身なので、正性は対角の非負性と行列式の非負性で判定できる。行列式が零なら非零な圧縮の階数は一、正なら二である。

$$
\begin{aligned}
\det C&=0\\
\rank J_\Phi&=1
\end{aligned}
$$

$C=\sum_l v_lv_l^*$ と正行列の固有分解で表せば $K_l=\diag(v_l)$ が Kraus 作用素になる。実際 $(K_lxK_l^*)_{ij}=(v_l)_i\overline{(v_l)_j}x_{ij}$ なので全項を足すと $C_{ij}x_{ij}$ になる。

$$
\begin{aligned}
x&=\begin{pmatrix}1 & 2 + i\\2 - i & 4\end{pmatrix}\\
\Phi(x)&=\begin{pmatrix}\frac{2}{3} & \frac{\sqrt{2} \left(2 + i\right)}{3}\\\frac{\sqrt{2} \left(2 - i\right)}{3} & \frac{4}{3}\end{pmatrix}
\end{aligned}
$$

最小個数は Choi 行列の階数である。なお $\Phi(I)$ は $C$ の対角だけを持つので、この例では写像は単位的ではない。
:::

:::example[二つの正固有値を持つ場合]
成分ごとの積で定めた Schur 写像 $\Phi(x)=C\circ x$ を考える。Choi 行列を作り、完全正性と必要な Kraus 作用素の個数を求める。

行列単位 $E_{ij}$ の像は $C_{ij}E_{ij}$ である。従って Choi 行列では対応する $e_i\otimes e_i$ の二方向だけに成分がある。

$$
\begin{aligned}
C&=\begin{pmatrix}\frac{2}{3} & \frac{1}{3}\\\frac{1}{3} & \frac{1}{3}\end{pmatrix}\\
J_\Phi&=\begin{pmatrix}\frac{2}{3} & 0 & 0 & \frac{1}{3}\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\\frac{1}{3} & 0 & 0 & \frac{1}{3}\end{pmatrix}
\end{aligned}
$$

この二方向への圧縮は $C$ 自身なので、正性は対角の非負性と行列式の非負性で判定できる。行列式が零なら非零な圧縮の階数は一、正なら二である。

$$
\begin{aligned}
\det C&=\frac{1}{9}\\
\rank J_\Phi&=2
\end{aligned}
$$

$C=\sum_l v_lv_l^*$ と正行列の固有分解で表せば $K_l=\diag(v_l)$ が Kraus 作用素になる。実際 $(K_lxK_l^*)_{ij}=(v_l)_i\overline{(v_l)_j}x_{ij}$ なので全項を足すと $C_{ij}x_{ij}$ になる。

$$
\begin{aligned}
x&=\begin{pmatrix}1 & 2 + i\\2 - i & 4\end{pmatrix}\\
\Phi(x)&=\begin{pmatrix}\frac{2}{3} & \frac{2}{3} + \frac{i}{3}\\\frac{2}{3} - \frac{i}{3} & \frac{4}{3}\end{pmatrix}
\end{aligned}
$$

最小個数は Choi 行列の階数である。なお $\Phi(I)$ は $C$ の対角だけを持つので、この例では写像は単位的ではない。
:::

:::example[複素位相を持つ場合]
成分ごとの積で定めた Schur 写像 $\Phi(x)=C\circ x$ を考える。Choi 行列を作り、完全正性と必要な Kraus 作用素の個数を求める。

行列単位 $E_{ij}$ の像は $C_{ij}E_{ij}$ である。従って Choi 行列では対応する $e_i\otimes e_i$ の二方向だけに成分がある。

$$
\begin{aligned}
C&=\begin{pmatrix}\frac{1}{3} & \frac{i}{3}\\- \frac{i}{3} & \frac{2}{3}\end{pmatrix}\\
J_\Phi&=\begin{pmatrix}\frac{1}{3} & 0 & 0 & \frac{i}{3}\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\- \frac{i}{3} & 0 & 0 & \frac{2}{3}\end{pmatrix}
\end{aligned}
$$

この二方向への圧縮は $C$ 自身なので、正性は対角の非負性と行列式の非負性で判定できる。行列式が零なら非零な圧縮の階数は一、正なら二である。

$$
\begin{aligned}
\det C&=\frac{1}{9}\\
\rank J_\Phi&=2
\end{aligned}
$$

$C=\sum_l v_lv_l^*$ と正行列の固有分解で表せば $K_l=\diag(v_l)$ が Kraus 作用素になる。実際 $(K_lxK_l^*)_{ij}=(v_l)_i\overline{(v_l)_j}x_{ij}$ なので全項を足すと $C_{ij}x_{ij}$ になる。

$$
\begin{aligned}
x&=\begin{pmatrix}1 & 2 + i\\2 - i & 4\end{pmatrix}\\
\Phi(x)&=\begin{pmatrix}\frac{1}{3} & \frac{i \left(2 + i\right)}{3}\\- \frac{i \left(2 - i\right)}{3} & \frac{8}{3}\end{pmatrix}
\end{aligned}
$$

最小個数は Choi 行列の階数である。なお $\Phi(I)$ は $C$ の対角だけを持つので、この例では写像は単位的ではない。
:::

## Arveson の拡張定理

:::theorem[Arveson の拡張]{#thm-oa-arveson-extension}
単位的 C*-環 $A$ の作用素系 $S$ 上の CP 写像 $\Phi:S\to B(H)$ は $A$ 上の CP 写像へ延長できる。単位元での値は保たれる。
:::
:::proof
まず $H=\C^d$ とする。$s_\Phi([a_{ij}])=\sum_{i,j}\inner{\Phi(a_{ij})e_j}{e_i}$ とおくと $M_d(S)$ 上の正値汎関数である。逆にこの式で与えられる汎関数の正性から $\Phi$ の CP 性が従う。$\xi_j=\sum_l c_{lj}e_l$ として $[a_{ij}]\in M_n(S)_+$ の増幅の二次形式を展開すると、それは $s_\Phi(C^*[a_{ij}]C)$ となる。ただし $C$ の $(j,l)$ 成分は $c_{lj}$ である。スカラー圧縮は正なので二次形式が非負である。任意の $s$ から $\inner{\Phi(a)e_j}{e_i}=s(E_{ij}\otimes a)$ として写像が復元される。

作用素系の正値汎関数 $s$ は実自己共役部分で $\abs{s(x)}\le s(1)\norm x$ を満たす。実 Hahn–Banach で全自己共役部分へ同ノルムで延長する。$0\le x\le1$ なら $\norm{1-x}\le1$ なので延長値は非負である。複素化すると正値汎関数の延長が得られ、それを $M_d(A)$ に適用して CP 写像を構成する。

一般の $H$ では各有限次元部分空間 $F$ への圧縮を上の方法で延長し、像を $F^\perp$ 上で零として $B(H)$ に入れる。Stinespring の定理によりノルムは $\norm{\Phi(1)}$ 以下である。各 $a$ について作用素ノルム有界球の WOT コンパクト性を使い、その積空間でクラスター点を取る。線形性と各増幅の二次形式の非負性はいずれも WOT の極限で保たれる。$S$ の値は $P_F\Phi(s)P_F\to\Phi(s)$ なので延長である。
:::

:::example[延長の自由度]
$S=\spn\{I_2,E_{12},E_{21}\}$ 上で $\varphi(1)=1$、$\varphi(E_{12})=\varphi(E_{21})=3/10$ とする。$M_2$ の状態としての全延長は、密度の対角を $t,1-t$ として
$$
D_t=\begin{pmatrix}t&3/10\\3/10&1-t\end{pmatrix},\qquad 1/10\le t\le9/10
$$
で与えられる。非対角成分は指定された値で決まり、正値性は対角非負と $t(1-t)\ge9/100$ に等しい。従って延長は存在し、範囲の内点で異なる $t$ を選ぶと第 1 座標の確率だけが異なる。いずれも状態なので、スカラー値の Stinespring 表示により完全正である。
:::

### CP 拡張を有限行列の変数として検査する

:::example[作用素系上の制約を Choi 行列に移す]
$S=\spn\{I,E_{12},E_{21}\}\subset M_2$ とし、スカラー値の単位的写像で $\Phi(E_{12})=z$、$\Phi(E_{21})=\bar z$ を要求する。拡張の Choi 行列は
$$
\begin{pmatrix}t&z\\\bar z&1-t\end{pmatrix}
$$
であり、正性と単位性の条件は $0\le t\le1$ および $t(1-t)\ge|z|^2$ になる。従って $|z|\le1/2$ なら少なくとも $t=1/2$ を選べる。$z=1/2$ では $t=1/2$ だけが可能であり、$z=0$ では閉区間全体を動ける。拡張定理は存在を保証するが、制約によって延長全体の形は変わる。
:::

:::example[非対角を測らない指定]
$S=\spn\{I,E_{12},E_{21}\}$ 上の単位的な汎関数を $E_{12}\mapsto z$ で指定する。$M_2$ への正値延長の集合を Choi 行列で求める。

対合を保つため $E_{21}$ の像は $\bar z$ であり、単位元の像を $1$ にするため対角値は $t,1-t$ となる。選択できるパラメータは $t$ だけである。

$$
\begin{aligned}
z&=0\\
J_{\Phi,\,t=1/2}&=\begin{pmatrix}\frac{1}{2} & 0\\0 & \frac{1}{2}\end{pmatrix}
\end{aligned}
$$

正性には $t\in[0,1]$ と $t(1-t)-|z|^2\ge0$ が必要十分になる。平方を完成して両端の値を求める。

$$
\begin{aligned}
|z|^2&=0\\
\sqrt{1-4|z|^2}&=1\\
\text{許される }t&=\left[0, 1\right]
\end{aligned}
$$

中央の $t=1/2$ を選ぶと密度行列は Choi 行列の転置である。これで指定通り $\Phi(E_{12})=z$ となる。任意の $x$ への値は密度との積のトレースである。

$$
\begin{aligned}
x&=\begin{pmatrix}1 & 2 + i\\2 - i & 4\end{pmatrix}\\
\Phi(x)&=\frac{5}{2}
\end{aligned}
$$

区間が一点に縮むと延長は一意である。区間に幅がある場合は対角行列単位に対する値を変えられる。
:::

:::example[小さい実の非対角指定]
$S=\spn\{I,E_{12},E_{21}\}$ 上の単位的な汎関数を $E_{12}\mapsto z$ で指定する。$M_2$ への正値延長の集合を Choi 行列で求める。

対合を保つため $E_{21}$ の像は $\bar z$ であり、単位元の像を $1$ にするため対角値は $t,1-t$ となる。選択できるパラメータは $t$ だけである。

$$
\begin{aligned}
z&=\frac{1}{4}\\
J_{\Phi,\,t=1/2}&=\begin{pmatrix}\frac{1}{2} & \frac{1}{4}\\\frac{1}{4} & \frac{1}{2}\end{pmatrix}
\end{aligned}
$$

正性には $t\in[0,1]$ と $t(1-t)-|z|^2\ge0$ が必要十分になる。平方を完成して両端の値を求める。

$$
\begin{aligned}
|z|^2&=\frac{1}{16}\\
\sqrt{1-4|z|^2}&=\frac{\sqrt{3}}{2}\\
\text{許される }t&=\left[\frac{1}{2} - \frac{\sqrt{3}}{4}, \frac{\sqrt{3}}{4} + \frac{1}{2}\right]
\end{aligned}
$$

中央の $t=1/2$ を選ぶと密度行列は Choi 行列の転置である。これで指定通り $\Phi(E_{12})=z$ となる。任意の $x$ への値は密度との積のトレースである。

$$
\begin{aligned}
x&=\begin{pmatrix}1 & 2 + i\\2 - i & 4\end{pmatrix}\\
\Phi(x)&=\frac{7}{2}
\end{aligned}
$$

区間が一点に縮むと延長は一意である。区間に幅がある場合は対角行列単位に対する値を変えられる。
:::

:::example[虚の非対角指定]
$S=\spn\{I,E_{12},E_{21}\}$ 上の単位的な汎関数を $E_{12}\mapsto z$ で指定する。$M_2$ への正値延長の集合を Choi 行列で求める。

対合を保つため $E_{21}$ の像は $\bar z$ であり、単位元の像を $1$ にするため対角値は $t,1-t$ となる。選択できるパラメータは $t$ だけである。

$$
\begin{aligned}
z&=\frac{i}{3}\\
J_{\Phi,\,t=1/2}&=\begin{pmatrix}\frac{1}{2} & \frac{i}{3}\\- \frac{i}{3} & \frac{1}{2}\end{pmatrix}
\end{aligned}
$$

正性には $t\in[0,1]$ と $t(1-t)-|z|^2\ge0$ が必要十分になる。平方を完成して両端の値を求める。

$$
\begin{aligned}
|z|^2&=\frac{1}{9}\\
\sqrt{1-4|z|^2}&=\frac{\sqrt{5}}{3}\\
\text{許される }t&=\left[\frac{1}{2} - \frac{\sqrt{5}}{6}, \frac{\sqrt{5}}{6} + \frac{1}{2}\right]
\end{aligned}
$$

中央の $t=1/2$ を選ぶと密度行列は Choi 行列の転置である。これで指定通り $\Phi(E_{12})=z$ となる。任意の $x$ への値は密度との積のトレースである。

$$
\begin{aligned}
x&=\begin{pmatrix}1 & 2 + i\\2 - i & 4\end{pmatrix}\\
\Phi(x)&=\frac{5}{2} - \frac{i \left(2 - i\right)}{3} + \frac{i \left(2 + i\right)}{3}
\end{aligned}
$$

区間が一点に縮むと延長は一意である。区間に幅がある場合は対角行列単位に対する値を変えられる。
:::

:::example[一意に決まる境界の指定]
$S=\spn\{I,E_{12},E_{21}\}$ 上の単位的な汎関数を $E_{12}\mapsto z$ で指定する。$M_2$ への正値延長の集合を Choi 行列で求める。

対合を保つため $E_{21}$ の像は $\bar z$ であり、単位元の像を $1$ にするため対角値は $t,1-t$ となる。選択できるパラメータは $t$ だけである。

$$
\begin{aligned}
z&=\frac{1}{2}\\
J_{\Phi,\,t=1/2}&=\begin{pmatrix}\frac{1}{2} & \frac{1}{2}\\\frac{1}{2} & \frac{1}{2}\end{pmatrix}
\end{aligned}
$$

正性には $t\in[0,1]$ と $t(1-t)-|z|^2\ge0$ が必要十分になる。平方を完成して両端の値を求める。

$$
\begin{aligned}
|z|^2&=\frac{1}{4}\\
\sqrt{1-4|z|^2}&=0\\
\text{許される }t&=\left\{\frac{1}{2}\right\}
\end{aligned}
$$

中央の $t=1/2$ を選ぶと密度行列は Choi 行列の転置である。これで指定通り $\Phi(E_{12})=z$ となる。任意の $x$ への値は密度との積のトレースである。

$$
\begin{aligned}
x&=\begin{pmatrix}1 & 2 + i\\2 - i & 4\end{pmatrix}\\
\Phi(x)&=\frac{9}{2}
\end{aligned}
$$

区間が一点に縮むと延長は一意である。区間に幅がある場合は対角行列単位に対する値を変えられる。
:::

:::example[絶対値が境界未満の複素指定]
$S=\spn\{I,E_{12},E_{21}\}$ 上の単位的な汎関数を $E_{12}\mapsto z$ で指定する。$M_2$ への正値延長の集合を Choi 行列で求める。

対合を保つため $E_{21}$ の像は $\bar z$ であり、単位元の像を $1$ にするため対角値は $t,1-t$ となる。選択できるパラメータは $t$ だけである。

$$
\begin{aligned}
z&=\frac{1}{4} + \frac{i}{3}\\
J_{\Phi,\,t=1/2}&=\begin{pmatrix}\frac{1}{2} & \frac{1}{4} + \frac{i}{3}\\\frac{1}{4} - \frac{i}{3} & \frac{1}{2}\end{pmatrix}
\end{aligned}
$$

正性には $t\in[0,1]$ と $t(1-t)-|z|^2\ge0$ が必要十分になる。平方を完成して両端の値を求める。

$$
\begin{aligned}
|z|^2&=\frac{25}{144}\\
\sqrt{1-4|z|^2}&=\frac{\sqrt{11}}{6}\\
\text{許される }t&=\left[\frac{1}{2} - \frac{\sqrt{11}}{12}, \frac{\sqrt{11}}{12} + \frac{1}{2}\right]
\end{aligned}
$$

中央の $t=1/2$ を選ぶと密度行列は Choi 行列の転置である。これで指定通り $\Phi(E_{12})=z$ となる。任意の $x$ への値は密度との積のトレースである。

$$
\begin{aligned}
x&=\begin{pmatrix}1 & 2 + i\\2 - i & 4\end{pmatrix}\\
\Phi(x)&=\frac{5}{2} + \left(\frac{1}{4} - \frac{i}{3}\right) \left(2 - i\right) + \left(\frac{1}{4} + \frac{i}{3}\right) \left(2 + i\right)
\end{aligned}
$$

区間が一点に縮むと延長は一意である。区間に幅がある場合は対角行列単位に対する値を変えられる。
:::

## 条件付き期待値と Tomiyama の定理

:::definition[条件付き期待値]{#def-oa-conditional-expectation}
単位元を共有する $B\subset A$ に対し、$E:A\to B$ が正で $E|_B=\id$、$E(bac)=bE(a)c$ を全 $b,c\in B$ で満たすとき :term[条件付き期待値]{en="conditional expectation" yomi="じょうけんつききたいち"} という。
:::

:::theorem[Tomiyama の定理]{#thm-oa-tomiyama}
単位元を共有する C*-部分環 $B\subset A$ への縮小的線形射影 $E$ は条件付き期待値で、完全正である。
:::
:::proof
まず $B$ の射影 $p$ に対し $x=(1-p)a$、$y=pE(x)$ と置く。$x^*y=0$ だから $\norm{x+ty}^2\le\norm x^2+t^2\norm y^2$ である。他方 $pE(x+ty)=(1+t)y$ だから縮小性から $(2t+1)\norm y^2\le\norm x^2$ となる。$t\to+\infty$ で $pE((1-p)a)=0$ を得る。$1-p$ にも適用すれば $E(pa)=pE(a)$ である。右側は $xy^*=0$ として同様に証明できる。$B$ が von Neumann 環なら自己共役元を有限スペクトル近似できるので、これで両側加群性が全ての $b,c$ に延びる。

一般には :ref[prop-oa-enveloping-vn] の包絡環 $A^{**},B^{**}$ に移る。$B^{**}\subset A^{**}$ は単位を共有する正規な *-埋め込みである。$E^{**}$ は同ノルムで $B^{**}$ に値を取る射影なので、射影での結論を $A,B$ に制限する。

$E(1)=1$ であり、像の状態と合成すると単位元で $1$ を取る縮小汎関数となる。:ref[lem-oa-state-extension] の正値性判定によりこれは正値である。状態は正錐を分離するので $E$ は正である。実際、状態による値が全て実数なら虚部は零であり、全て非負ならスペクトルの下端を取る状態から正性が従う。

最後の完全正性は任意の条件付き期待値に成り立つ。$X=[x_{ij}]\ge0$、列 $b=(b_j)\in B^n$ に対し $b^*E^{(n)}(X)b=E(b^*Xb)\ge0$ である。もし自己共役 $Y=E^{(n)}(X)$ に負のスペクトルがあれば、負の部分に台を持つ連続関数 $f$ を選んで $u=f(Y)$ とし、$u^*Yu$ を非零な負元にできる。正元 $-u^*Yu$ の全対角成分が零なら、その正平方根の各列の二乗和が零となり、その元自体が零になる。従って $u$ の少なくとも一つの列 $b$ に対して $b^*Yb$ は非零な負元となり、矛盾する。これで完全正性が分かり、Stinespring の表示から $\norm E=1$ も出る。
:::

:::example[部分トレース]
$A=M_n\otimes M_m$、$B=M_n\otimes I_m$ とする。$E(a\otimes b)=a\otimes I_m\,\tr(b)/m$ は第 2 因子の各座標への圧縮の平均を使うので CP である。第 1 因子の左右からの積を取り出せる。特に $E(E_{11}\otimes\diag(1,3))=2E_{11}\otimes I_2$ である。
:::

### トレースを保存する期待値

有限 von Neumann 環の忠実正規トレース $\tau$ を保存する $N$ への期待値 $E$ があるとする。$b\in N$ について $\tau(b^*E(a))=\tau(E(b^*a))=\tau(b^*a)$ である。従って $E$ と $F$ が同じ条件を満たすなら、$b=E(a)-F(a)$ をこの式に代入して $\tau(b^*b)=0$ を得る。忠実性により $b=0$ なので一意である。この等式は $L^2$ 内積で $E(a)$ が $N$ への直交射影であることも示す。

:::proposition[有限トレースを保存する期待値の存在]{#prop-oa-tracial-expectation}
$M$ が忠実な有限正規トレース $\tau$ を持ち、$N\subset M$ が単位を共有する von Neumann 部分環なら、$\tau$ を保存する $N$ への期待値 $E$ が存在する。$E$ は忠実かつ正規である。
:::

:::proof
$P$ を $L^2(M,\tau)$ から $N$ の閉包への直交射影とする。まず $0\le x\le1$ に対して $Px$ を近似する $N$ の自己共役元の列 $y_j$ を取る。これは対合が $L^2$ で反線形等長であり、$P$ と可換だからである。連続関数カルキュラスで $z_j=\max(0,\min(y_j,1))$ とする。$y=y_j,z=z_j$ に対し $y-z=(y-1)_+-(-y)_+$ である。前者の台で $z=1$、後者の台で $z=0$ なので
$$
\Re\tau((y-z)(x-z))
=\tau((y-1)_+(x-1))-\tau((-y)_+x)\le0 .
$$
各項の符号は、トレース性により一方の正因子の平方根を他方の両側に掛ければ分かる。従って $\norm{x-z_j}_2^2\le\norm{x-y_j}_2^2$ であり、直交分解により $z_j\to Px$ が $L^2$ で成り立つ。

$N$ の正縮小元の集合は σ 弱コンパクトなので、列 $(z_j)$ の σ 弱集積点 $z$ を取れる。$\tau(b^*z_j)$（$b\in N$）の極限は $Px$ との $L^2$ 内積であり、正規性から $\tau(b^*z)$ と一致する。従って $Px=z$ である。任意の元は正縮小元の線形結合なので、$P$ は $M$ を $N$ へ写し、そこで $E$ を定める。$E$ は正で $N$ を固定する。$L^2$ 上で $N$ の左右の作用はその完備化を還元するため $P$ と可換となり、両加群性を得る。また $1\in N$ との内積を取れば $\tau$ が保たれる。

$E(x^*x)=0$ なら $\tau(x^*x)=0$ から $x=0$ となる。さらに $x_i\uparrow x$ に対して $y=\sup_i E(x_i)\le E(x)$ と置くと、$\tau(y)=\sup_i\tau(x_i)=\tau(x)=\tau(E(x))$ である。忠実性により $y=E(x)$ なので、$E$ は正規である。
:::

:::example[異なる重みで期待値を作る]
$M=M_n\otimes M_m$、$N=M_n\otimes I_m$ とする。第 2 因子の状態 $\omega_D(a)=\tr(Da)$ に対して $\id\otimes\omega_D$ に単位元をテンソルした写像は $N$ を固定し $N$ 両加群的である。$\omega_D$ の Kraus 表現を増幅すれば完全正になり、従って期待値である。第 2 因子上の $I_n\otimes E_{jj}$ の像は $D_{jj}I$ なので、$D=\diag(3/4,1/4)$ とすれば規格化トレースで $1/2$ だった値が $3/4$ または $1/4$ へ変わる。この期待値は規格化トレースを保存しない。
:::

:::example[入れ子の部分環に期待値を取る]
$M=M_2\otimes M_2$、$N=M_2\otimes1$、$D\subset N$ を第 1 因子の対角部分とする。規格化トレースを保つ期待値 $E_N,E_D$ は $E_D E_N=E_D$ を満たす。実際、$E_D E_N$ は $D$ を固定する UCP 射影で、トレースを保存するため一意性が適用できる。例えば $x=\begin{pmatrix}2&1\\1&0\end{pmatrix}\otimes\diag(1,3)$ に対し $E_N(x)=\begin{pmatrix}4&2\\2&0\end{pmatrix}\otimes1$ であり、次に対角を取ると $E_D(x)=\diag(4,0)\otimes1$ になる。$E_N E_D=E_D$ も、$E_D$ の像が既に $N$ に入ることから出る。二段の近似の測定誤差を分けるときに用いる。
:::

:::example[二つの座標をまとめる期待値]
座標をいくつかのブロックに分け、各ブロックへの射影 $p_j$ を取る。写像 $E(x)=\sum_jp_jxp_j$ で測定して、何が保たれるか求める。

この写像は各射影による圧縮の和で、$\sum_jp_j=1$ だから単位的完全正である。出力では同じブロックに属する行と列の組だけを残す。

$$
\begin{aligned}
\text{各座標のブロック番号}&=\left( 0, \  0, \  1\right)\\
a&=\begin{pmatrix}1 & 2 & 3\\2 & 3 & 1\\3 & 1 & 0\end{pmatrix}\\
E(a)&=\begin{pmatrix}1 & 2 & 0\\2 & 3 & 0\\0 & 0 & 0\end{pmatrix}
\end{aligned}
$$

両加群性は、ブロック対角な $b,c$ が全 $p_j$ と可換であるため $E(bxc)=bE(x)c$ と計算できる。さらに $E(x)$ は既にブロック対角なので $E(E(x))=E(x)$ である。

$$
\begin{aligned}
a-E(a)&=\begin{pmatrix}0 & 0 & 3\\0 & 0 & 1\\3 & 1 & 0\end{pmatrix}\\
\tr(a)=\tr(E(a))&=4
\end{aligned}
$$

行列単位の直交性を使うと $E(a)$ と $a-E(a)$ はトレース内積で直交する。実際、片方はブロック内の位置だけ、もう片方はブロック間だけに成分を持つ。

$$
\begin{aligned}
\tr(E(a)^*(a-E(a)))&=0\\
\tr((a-E(a))^*(a-E(a)))&=20
\end{aligned}
$$

この期待値は各ブロック内の非対角成分も保つ。対角だけを取り出す期待値になるのは、各ブロックが一次元の場合である。
:::

:::example[離れた座標が同じブロックになる]
座標をいくつかのブロックに分け、各ブロックへの射影 $p_j$ を取る。写像 $E(x)=\sum_jp_jxp_j$ で測定して、何が保たれるか求める。

この写像は各射影による圧縮の和で、$\sum_jp_j=1$ だから単位的完全正である。出力では同じブロックに属する行と列の組だけを残す。

$$
\begin{aligned}
\text{各座標のブロック番号}&=\left( 0, \  1, \  0\right)\\
a&=\begin{pmatrix}2 & 0 & 1\\0 & 4 & 3\\1 & 3 & 3\end{pmatrix}\\
E(a)&=\begin{pmatrix}2 & 0 & 1\\0 & 4 & 0\\1 & 0 & 3\end{pmatrix}
\end{aligned}
$$

両加群性は、ブロック対角な $b,c$ が全 $p_j$ と可換であるため $E(bxc)=bE(x)c$ と計算できる。さらに $E(x)$ は既にブロック対角なので $E(E(x))=E(x)$ である。

$$
\begin{aligned}
a-E(a)&=\begin{pmatrix}0 & 0 & 0\\0 & 0 & 3\\0 & 3 & 0\end{pmatrix}\\
\tr(a)=\tr(E(a))&=9
\end{aligned}
$$

行列単位の直交性を使うと $E(a)$ と $a-E(a)$ はトレース内積で直交する。実際、片方はブロック内の位置だけ、もう片方はブロック間だけに成分を持つ。

$$
\begin{aligned}
\tr(E(a)^*(a-E(a)))&=0\\
\tr((a-E(a))^*(a-E(a)))&=18
\end{aligned}
$$

この期待値は各ブロック内の非対角成分も保つ。対角だけを取り出す期待値になるのは、各ブロックが一次元の場合である。
:::

:::example[対角部分環への期待値]
座標をいくつかのブロックに分け、各ブロックへの射影 $p_j$ を取る。写像 $E(x)=\sum_jp_jxp_j$ で測定して、何が保たれるか求める。

この写像は各射影による圧縮の和で、$\sum_jp_j=1$ だから単位的完全正である。出力では同じブロックに属する行と列の組だけを残す。

$$
\begin{aligned}
\text{各座標のブロック番号}&=\left( 0, \  1, \  2\right)\\
a&=\begin{pmatrix}1 & 2 & 1\\2 & 3 & 4\\1 & 4 & 9\end{pmatrix}\\
E(a)&=\begin{pmatrix}1 & 0 & 0\\0 & 3 & 0\\0 & 0 & 9\end{pmatrix}
\end{aligned}
$$

両加群性は、ブロック対角な $b,c$ が全 $p_j$ と可換であるため $E(bxc)=bE(x)c$ と計算できる。さらに $E(x)$ は既にブロック対角なので $E(E(x))=E(x)$ である。

$$
\begin{aligned}
a-E(a)&=\begin{pmatrix}0 & 2 & 1\\2 & 0 & 4\\1 & 4 & 0\end{pmatrix}\\
\tr(a)=\tr(E(a))&=13
\end{aligned}
$$

行列単位の直交性を使うと $E(a)$ と $a-E(a)$ はトレース内積で直交する。実際、片方はブロック内の位置だけ、もう片方はブロック間だけに成分を持つ。

$$
\begin{aligned}
\tr(E(a)^*(a-E(a)))&=0\\
\tr((a-E(a))^*(a-E(a)))&=42
\end{aligned}
$$

この期待値は各ブロック内の非対角成分も保つ。対角だけを取り出す期待値になるのは、各ブロックが一次元の場合である。
:::

:::example[一つのブロックが全体を占める]
座標をいくつかのブロックに分け、各ブロックへの射影 $p_j$ を取る。写像 $E(x)=\sum_jp_jxp_j$ で測定して、何が保たれるか求める。

この写像は各射影による圧縮の和で、$\sum_jp_j=1$ だから単位的完全正である。出力では同じブロックに属する行と列の組だけを残す。

$$
\begin{aligned}
\text{各座標のブロック番号}&=\left( 0, \  0, \  0\right)\\
a&=\begin{pmatrix}3 & 2 & 1\\2 & 1 & 2\\1 & 2 & 0\end{pmatrix}\\
E(a)&=\begin{pmatrix}3 & 2 & 1\\2 & 1 & 2\\1 & 2 & 0\end{pmatrix}
\end{aligned}
$$

両加群性は、ブロック対角な $b,c$ が全 $p_j$ と可換であるため $E(bxc)=bE(x)c$ と計算できる。さらに $E(x)$ は既にブロック対角なので $E(E(x))=E(x)$ である。

$$
\begin{aligned}
a-E(a)&=\begin{pmatrix}0 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{pmatrix}\\
\tr(a)=\tr(E(a))&=4
\end{aligned}
$$

行列単位の直交性を使うと $E(a)$ と $a-E(a)$ はトレース内積で直交する。実際、片方はブロック内の位置だけ、もう片方はブロック間だけに成分を持つ。

$$
\begin{aligned}
\tr(E(a)^*(a-E(a)))&=0\\
\tr((a-E(a))^*(a-E(a)))&=0
\end{aligned}
$$

この期待値は各ブロック内の非対角成分も保つ。対角だけを取り出す期待値になるのは、各ブロックが一次元の場合である。
:::

:::example[複素成分を含む場合]
座標をいくつかのブロックに分け、各ブロックへの射影 $p_j$ を取る。写像 $E(x)=\sum_jp_jxp_j$ で測定して、何が保たれるか求める。

この写像は各射影による圧縮の和で、$\sum_jp_j=1$ だから単位的完全正である。出力では同じブロックに属する行と列の組だけを残す。

$$
\begin{aligned}
\text{各座標のブロック番号}&=\left( 0, \  1, \  1\right)\\
a&=\begin{pmatrix}1 & i & 2\\- i & 4 & 1\\2 & 1 & 2\end{pmatrix}\\
E(a)&=\begin{pmatrix}1 & 0 & 0\\0 & 4 & 1\\0 & 1 & 2\end{pmatrix}
\end{aligned}
$$

両加群性は、ブロック対角な $b,c$ が全 $p_j$ と可換であるため $E(bxc)=bE(x)c$ と計算できる。さらに $E(x)$ は既にブロック対角なので $E(E(x))=E(x)$ である。

$$
\begin{aligned}
a-E(a)&=\begin{pmatrix}0 & i & 2\\- i & 0 & 0\\2 & 0 & 0\end{pmatrix}\\
\tr(a)=\tr(E(a))&=7
\end{aligned}
$$

行列単位の直交性を使うと $E(a)$ と $a-E(a)$ はトレース内積で直交する。実際、片方はブロック内の位置だけ、もう片方はブロック間だけに成分を持つ。

$$
\begin{aligned}
\tr(E(a)^*(a-E(a)))&=0\\
\tr((a-E(a))^*(a-E(a)))&=10
\end{aligned}
$$

この期待値は各ブロック内の非対角成分も保つ。対角だけを取り出す期待値になるのは、各ブロックが一次元の場合である。
:::

## 演習問題

::::exercise{level=basic}
$V:\C\to\C^2$ を $Vz=ze_2$ とする。$a\mapsto V^*aV$ を計算せよ。
:::solution
$V^*$ は第 2 座標を取り出すので $\Phi(a)=a_{22}$ である。$V^*V=1$ で、恒等表現の圧縮だから UCP である。
:::
::::
::::exercise{level=basic}
$\Phi(a)=\tfrac12 a$（$a\in M_2$）の Choi 行列の固有値を求めよ。
:::solution
$C_\Phi=\Omega\Omega^*/2$ で $\norm{\Omega}^2=2$ である。$\Omega$ 上で固有値 $1$、直交補で $0$ となる。
:::
::::
::::exercise{level=basic}
行列 $K=\diag(1,2)$ による $\Phi(a)=KaK^*$ の $\Phi(E_{12})$ と $\norm\Phi$ を求めよ。
:::solution
両側の積から $\Phi(E_{12})=2E_{12}$ である。$\Phi(1)=\diag(1,4)$ なので $\norm\Phi=4$ である。
:::
::::
::::exercise{level=basic}
$M_2$ の対角期待値に $a=\begin{pmatrix}2&1\\1&3\end{pmatrix}$ を入れ、$E(a^2)-E(a)^2$ を求めよ。
:::solution
$a^2$ の対角は $5,10$ である。$E(a)^2$ の対角は $4,9$ なので差は $I_2$ である。
:::
::::
::::exercise{level=advanced}
$\Phi$ が UCP のとき、ユニタリ $u$ がユニタリに写れば乗法領域に入ることを示せ。
:::solution
$\Phi(u^*u)=\Phi(uu^*)=1$ で、$\Phi(u)^*\Phi(u)=\Phi(u)\Phi(u)^*=1$ でもある。膨張で $\pi(u)$ と $\pi(u)^*$ はともに $VH$ を保つから積を左右で取り出せる。
:::
::::
::::exercise{level=advanced}
CP 写像の合成が CP であることを示し、UCP についても同様か答えよ。
:::solution
各 $n$ で合成の増幅は増幅の合成である。正元は一つ目で正元になり、二つ目でも正元になる。両者が単位を保てば合成も単位を保つ。
:::
::::
::::exercise{level=advanced}
$M_n$ 上の $a\mapsto\tr(a)I_n/n$ の Choi 行列と最小 Kraus 個数を求めよ。
:::solution
$\Phi(E_{ij})=\delta_{ij}I_n/n$ だから $C_\Phi=I_{n^2}/n$ となり、階数 $n^2$ が最小個数である。実際 $K_{ij}=E_{ij}/\sqrt n$ が表示を与える。
:::
::::
::::exercise{level=advanced}
完全正写像 $\Phi:A\to B(H)$ に対し $\Phi(1)=0$ なら $\Phi=0$ を示せ。
:::solution
Stinespring の表示で $\Phi(1)=V^*V$ である。$\norm{V\xi}^2=0$ が全 $\xi$ で成り立つため $V=0$ で、表示から $\Phi(a)=0$ である。
:::
::::

## 確認テスト

::::quiz{answer="2" tol="1e-6"}
$M_2$ の恒等写像の Choi 行列の最大固有値を答えよ。
:::explanation
階数一行列 $\Omega\Omega^*$ の非零固有値は $\inner{\Omega}{\Omega}=2$ である。
:::
::::
::::quiz
行列の転置が満たす性質はどれか。

- [x] 正値
- [ ] $2$ 正値
- [ ] 完全正値

:::explanation
正行列 $a^*a$ の転置は $a^{\mathsf T}\bar a$ で正である。それでも増幅で正の $[E_{ij}]$ は交換作用素に移り、反対称ベクトルに負の固有値を持つ。
:::
::::
::::quiz{answer="9" tol="1e-6"}
$M_3$ 上の $a\mapsto\tr(a)I_3/3$ の最小 Kraus 個数を答えよ。
:::explanation
Choi 行列は $I_9/3$ であり、階数 $9$ が最小個数である。
:::
::::
::::quiz
UCP 写像の Stinespring 表示で $V$ が満たす条件はどれか。

- [x] $V^*V=1$
- [ ] いつでも $VV^*=1$
- [ ] $V=0$

:::explanation
$\pi$ は単位的なので $\Phi(1)=V^*V=1$ である。行列のベクトル状態では $\C$ から $\C^2$ への等長埋め込みとなるので全射性は一般にはない。
:::
::::

---
title: "総和法と良い核"
description: "正の近似単位で周期関数を復元し、係数の減衰、近似精度、通常部分和の発散を比較する。"
notation:
  - symbol: 'f*g'
    meaning: "正規化測度に関する円周畳み込み"
  - symbol: 'F_N'
    meaning: "Fejér 核"
  - symbol: 'P_r'
    meaning: "円板の Poisson 核"
  - symbol: '\sigma_N f'
    meaning: "第 N Cesàro 平均"
  - symbol: '\omega_f'
    meaning: "連続度"
---

## この章の目標

- 円周上の畳み込みと近似単位の収束を証明する。
- Fejér 核と Poisson 核の係数を求め、近似と境界値の復元に使う。
- 滑らかさと係数の評価を両方向から調べる。
- Lebesgue 定数の増大を利用して発散する連続関数の存在を示す。

## 畳み込みと良い核

周波数ごとの係数を掛け合わせる操作を、元の変数では移動した関数の重み付き平均として表せる。重みを正に選ぶと、誤差の見積もりに打ち消しを仮定する必要がなくなる。以下では $dm(t)=dt/(2\pi)$ とする。

:::definition[円周畳み込み]{#def-four-circle-convolution}
$f,g\in L^1(\mathbb T)$ に対し、:term[畳み込み]{en="convolution" yomi="たたみこみ"} を

$$
(f*g)(x)=\int_{\mathbb T}f(x-t)g(t)\,dm(t)
$$

と定める。被積分関数の絶対値の二重積分は $\|f\|_1\|g\|_1$ なので、積分は a.e. で定義できる。
:::

:::proposition[畳み込みと係数]{#prop-four-circle-product}
$f,g,h\in L^1$ では $f*g=g*f$、$(f*g)*h=f*(g*h)$、$\widehat{f*g}(n)=\widehat f(n)\widehat g(n)$ である。$1\le p\le\infty$ で $f\in L^p$ なら $\|f*g\|_p\le\|f\|_p\|g\|_1$ である。
:::

:::proof
可換性は $u=x-t$ の置換による。結合性は $\iint f(x-s-t)g(s)h(t)\,dm(s)dm(t)$ の順序交換による。絶対値の $x,s,t$ 積分は三つのノルムの積なので交換できる。係数の計算では $x=u+t$ を代入して指数を $e^{-inu}e^{-int}$ に分ける。ノルム評価は積分版 Minkowski の不等式（:ref[thm-meas-minkowski-integral]）と移動でノルムが保存されることに従う。$p=\infty$ では Fubini により a.e. の $x$ について $|f(x-t)|\le\|f\|_\infty$ が a.e. の $t$ で成り立つ。
:::

:::example[一つの周波数を抽出する核]
$e_k(t)=e^{ikt}$ とする。$f*e_k=e_k\widehat f(k)$ であり、直交性によって $(f*e_k)*e_k=f*e_k$ となる。$g=\cos t$ を用いれば、例えば $f=4\sin t+3\cos(2t)$ に対して $f*g=2\sin t$ となる。
:::

:::definition[良い核]{#def-four-good-kernel}
$K_j\in L^1(\mathbb T)$ が $\int K_j\,dm=1$、$\sup_j\|K_j\|_1<\infty$、任意の $0<\delta<\pi$ について $\int_{|t|\ge\delta}|K_j(t)|\,dm(t)\to0$ を満たすとき、:term[近似単位]{en="approximate identity" yomi="きんじたんい"} または良い核という。連続パラメータでも同じ条件を対応する極限で課す。
:::

:::theorem[近似単位による復元]{#thm-four-good-kernel}
$f\in C(\mathbb T)$ なら $f*K_j\to f$ が一様に成り立つ。$f\in L^p$、$1\le p<\infty$ なら $L^p$ ノルムでも収束する。$f\in L^\infty$ が $x$ で連続な代表元を持てば、その $x$ で $f*K_j(x)\to f(x)$ である。
:::

:::proof
$B=\sup_j\|K_j\|_1$ と置く。差は

$$
f*K_j-f=\int K_j(t)(\tau_t f-f)\,dm(t).
$$

$|t|<\delta$ の $L^p$ ノルムは $B\sup_{|t|<\delta}\|\tau_tf-f\|_p$ で抑えられる。連続関数には一様連続性を用いる。有限 $p$ では区間指示関数が移動に連続であり、それらの有限結合の稠密性を使う。稠密性は単関数近似と測度の正則性に従う。残りの積分のノルムは $2\|f\|_p\int_{|t|\ge\delta}|K_j|\,dm$ 以下である。先に $\delta$ を小さくし、次に $j$ を大きくする。一点での収束は同じ分割で連続性を $|f(x-t)-f(x)|$ に適用する。
:::

:::example[区間平均]
$0<\eps<\pi$ で $K_\eps(t)=\pi\eps^{-1}1_{[-\eps,\eps]}(t)$ と置く。質量は $1$ で、$\eps<\delta$ で外側の質量はゼロである。$f=\cos(2t)$ では $f*K_\eps=(\sin(2\eps)/(2\eps))f$ になる。有限 $\eps$ での振幅の低下は、平均が一点評価と異なることを表す。
:::

### 正則化と移動の交換

$t$ に関する小さな移動を $\tau_h f(x)=f(x-h)$ と定めると、畳み込みには $\tau_h(f*g)=(\tau_h f)*g=f*(\tau_hg)$ が成り立つ。最後の等式は積分で $t$ を $t-h$ へ移すと分かる。とくに $f\in L^\infty,g\in L^1$ なら

$$
\sup_x|(f*g)(x-h)-(f*g)(x)|\le\|f\|_\infty\|\tau_hg-g\|_1
$$

なので、各点で積分により定義した $f*g$ は一様連続になる。$f\in L^1$ で $g$ が連続な場合にも、$f$ と $g$ の役割を交換して同じ結論を得る。さらに $g$ が $C^k$ 級なら $(f*g)^{(j)}=f*g^{(j)}$（$0\le j\le k$）である。実際、$f*g(x)=\int f(t)g(x-t)\,dm(t)$ の差商は、$|f(t)|$ と $g'$ の上限の積で支配される。積分の中で微分でき、これを繰り返せばよい。

:::example[矩形パルスを局所平均で正則化する]
$0<\eps<a<\pi/4$ とし、$f=1_{[-a,a]}$ を $K_\eps$ で平均する。$|x|\le a-\eps$ では区間 $[x-\eps,x+\eps]$ の全体がパルスに入るため平均は $1$ である。$a-\eps<|x|<a+\eps$ では重なりの長さが $a+\eps-|x|$ なので平均は $(a+\eps-|x|)/(2\eps)$ になる。残りは $0$ である。$p<\infty$ の誤差は二つの幅 $2\eps$ の境界層から来る。各境界の両側で誤差の大きさは $s/(2\eps)$（$0\le s\le\eps$）を取るから、

$$
\|f*K_\eps-f\|_p^p=\frac4{2\pi}\int_0^\eps\left(\frac{s}{2\eps}\right)^p ds
=\frac{2^{1-p}\eps}{\pi(p+1)} .
$$

一方、境界にいくらでも近い正の長さの弧上で絶対誤差が $1/2$ に近づくため、$L^\infty$ の誤差は $1/2$ である。
:::

### 正値性と濃度が果たす別々の役割

積分が $1$ の非負関数 $K$ は、入力の値が閉区間 $[a,b]$ に入れば出力もそこへ入れる。これは $aK(t)\le f(x-t)K(t)\le bK(t)$ を積分した結果である。より一般に、複素数値の出力は入力の値域の閉凸包に入る。ある閉半平面が入力の値域を含めば、半平面を定める実一次式に同じ積分の議論を適用できる。そのような半平面の共通部分が閉凸包である。濃度の条件は、平均する範囲をある点へ近づける。正値のまま質量を集中できる核なら、平均の値域を保ちながら復元できる。

:::example[符号付き核でも一様収束する]
$K_\eps$ を区間平均の核として $L_\eps=2K_\eps-K_{2\eps}$ を考える。質量は $1$、ノルムは高々 $3$ で、台は $|t|\le2\eps$ に入る。したがって良い核である。具体的に内側 $|t|\le\eps$ で $L_\eps=3\pi/(2\eps)$、外側 $\eps<|t|\le2\eps$ で $L_\eps=-\pi/(2\eps)$ となり、絶対質量は $2$ である。$f(t)=\cos t$ の乗数は $2\sin\eps/\eps-\sin(2\eps)/(2\eps)=1+\eps^2/3+O(\eps^4)$ なので、十分小さい $\eps>0$ で、原点の出力は入力の最大値 $1$ を超える。近似単位であっても、有限の尺度での値域保存には正値性が効く。
:::

### 良い核の条件と各点収束の限界

近似単位であるだけでは、有界でない $f$ の一点での収束や、全 Lebesgue 点での収束は保証できない。核が原点から離れたところにごく小さな質量を置き、その位置で入力が大きくなる場合がある。例えば $a=\pi/2$ とし、$f(x)=|x+a|^{-1/2}$ を $|x+a|<1/4$ で取り、それ以外を $0$ とする。点 $-a$ 自身では値を $0$ と決める。$\eps\downarrow0$ で

$$
L_\eps=(1-\sqrt\eps)K_\eps+\sqrt\eps\,\tau_a K_\eps
$$

は正の良い核である。実際、後ろの項の質量は $\sqrt\eps$ で消え、前の項は原点に集中する。しかし $f*L_\eps(0)=2$ である。前の項は積分でゼロになり、後ろの項は $\sqrt\eps(2\eps)^{-1}\int_{-\eps}^{\eps}|u|^{-1/2}\,du=2$ を与える。$f$ は原点近傍で恒等的にゼロなので、原点は連続点であり Lebesgue 点でもある。各点収束に Fejér 核や Poisson 核を使う際には、このような遠い尖峰を排除する点ごとの核の上界が活躍する。

## Fejér 核と Cesàro 和

負の重みを含む $D_N$ を直接用いる代わりに、次数の異なる部分和を同じ重さでまとめる。次数を増やすごとに高周波へ急に切り替える操作が、線形な重みに置き換わる。

:::definition[Cesàro 平均]{#def-four-cesaro}
$\sigma_N f=N^{-1}\sum_{j=0}^{N-1}S_j f$（$N\ge1$）を :term[Cesàro 平均]{en="Cesàro mean" yomi="ちぇざろへいきん"} とする。核 $F_N=N^{-1}\sum_{j=0}^{N-1}D_j$ を :term[Fejér 核]{en="Fejér kernel" yomi="ふぇいぇーるかく"} という。
:::

:::proposition[Fejér 核の三つの表示]{#prop-four-fejer-kernel}
$$
F_N(t)=\sum_{|n|<N}\left(1-\frac{|n|}{N}\right)e^{int}
=\frac1N\left|\sum_{j=0}^{N-1}e^{ijt}\right|^2
=\frac{\sin^2(Nt/2)}{N\sin^2(t/2)}.
$$

商は $0$ で $N$ と補う。$F_N$ は偶関数の良い核であり、$\sigma_N f=f*F_N$ である。
:::

:::proof
$D_j$ の和では第 $n$ 項が $N-|n|$ 回出る。絶対値二乗の展開でも差が $n$ の添字対は同じ個数である。最後の式は等比級数で得られる。正値性と定数項 $1$ から $\|F_N\|_1=1$ である。$|t|\ge\delta$ で $F_N(t)\le (N\sin^2(\delta/2))^{-1}$ なので質量が集中する。部分和の核表示（:ref[prop-four-series-dirichlet]）を平均すれば最後の主張が従う。
:::

:::example[低周波と高周波の重み]
$N=4$ では $F_4=1+(3/2)\cos t+\cos(2t)+(1/2)\cos(3t)$ である。$f=3+4\cos t+4\cos(3t)+4\cos(4t)$ なら $\sigma_4f=3+3\cos t+\cos(3t)$ である。周波数 $3$ は残るが振幅は $1/4$ になる。
:::

左右の極限 $L_\pm$ がある有界関数では $\sigma_N f(x)\to(L_++L_-)/2$ となる。実際、偶核なので各片側に質量 $1/2$ があり、近傍で片側の値との差は小さい。外側の寄与は質量集中で消える。例えば $\sgn(\sin t)$ の Fejér 平均は偶核と奇性から原点でゼロである。また $-1\le f\le1$ なら同じ範囲に平均値が入る。したがって :ref[thm-four-series-gibbs] 型の値域外への超過は生じない。

:::example[点に集まる高さと幅]
$N=8$ では $F_8(0)=8$、$t=\pi/4$ で $F_8(\pi/4)=0$ である。この零点は $\sin(8t/2)=0$ から最初の正の零点である。一般に $N\ge2$ なら最初の正の零点は $2\pi/N$ であり、中心の高さは $N$ になる。$N=1$ では $F_1=1$ なので零点はない。どの $N$ でも核の質量は $1$ である。
:::

### Cesàro 平均を繰り返すと重みはどう変わるか

より高い階の :term[Cesàro 総和法]{en="Cesàro summation method" yomi="ちぇざろそうわほう"} では、$A_j^{(\alpha)}=\binom{j+\alpha}{j}$ を $\alpha$ が非負整数の場合に使う。整数 $r\ge1$ に対して $N$ 次の平均を

$$
\sigma_N^{(r)} f=\frac1{A_N^{(r)}}\sum_{j=0}^N A_{N-j}^{(r-1)}S_jf
$$

と置く。ここでは次数が $N$ 以下であり、既出の $\sigma_{N+1}$ が $r=1$ に対応する。$A_0^{(r)}=1$ と Pascal の関係 $A_k^{(r)}-A_{k-1}^{(r)}=A_k^{(r-1)}$ を使うと、重みの総和は $A_N^{(r)}$ になる。また $n$ 項の乗数は $j\ge|n|$ の重みを足した $A_{N-|n|}^{(r)}/A_N^{(r)}$ であり、$|n|>N$ ならゼロである。

:::proposition[整数階の Cesàro 核の正値性]{#prop-four-cesaro-higher}
$r\ge1$ の上記の平均は質量 $1$ の非負な核による畳み込みで表せる。固定した $r$ について、$N\to\infty$ でこの核は近似単位になる。
:::

:::proof
$r=1$ は Fejér 核である。$r\ge2$ では $S_j=(j+1)\sigma_{j+1}-j\sigma_j$ を代入し、隣の項と結合すると

$$
A_N^{(r)}\sigma_N^{(r)} f=(N+1)\sigma_{N+1}f+
\sum_{j=0}^{N-1}(j+1)A_{N-j}^{(r-2)}\sigma_{j+1}f .
$$

各係数は非負であり、定数関数を入れれば正規化後に和が $1$ になる。固定された $j$ に対する重みは $O(N^{r-2}/N^r)$ で消える。従って $\delta$ の外側の質量を評価する際に、有限個の小さい $j$ の寄与は消え、大きい $j$ の寄与は Fejér 核の外側質量で一様に小さくできる。
:::

:::example[二階平均のモード重み]
$N=3,r=2$ なら、分母は $A_3^{(2)}=10$ である。周波数 $0,1,2,3$ の乗数は順に $1,6/10,3/10,1/10$ となる。同じ最高次数を持つ Fejér 核 $F_4$ では $1,3/4,1/2,1/4$ である。二階平均が各モードをより強く抑えることは、$r=2$ の乗数が $(N-|n|+1)(N-|n|+2)/((N+1)(N+2))$ の二因子に分かれることからも分かる。
:::

### Fejér 平均で見た跳躍の拡大極限

Gibbs 現象と同じ倍率で跳躍を観察すると、核を替えたときの差がわかる。$v=Nt$ と置き、$F_N(t)\,dm(t)$ を実数上へ移す。密度を区間 $|v|\le\pi N$ で $F_N(v/N)/(2\pi N)$、外側でゼロとすると、各 $v$ で

$$
k(v)=\frac1{2\pi}\left(\frac{\sin(v/2)}{v/2}\right)^2
$$

に収束する。$F_N(t)\le\min(N,\pi^2/(Nt^2))$ から、密度は $C\min(1,|v|^{-2})$ で一様に抑えられる。優収束により $k$ の質量は $1$、$L^1$ で密度は $k$ に近づく。$f$ が有界で原点に片側極限 $L_\pm$ を持つなら、固定した $y$ について $v\ne y$ で $f((y-v)/N)$ は $v<y$ のとき $L_+$、$v>y$ のとき $L_-$ に近づく。再度優収束を適用すると

$$
\lim_{N\to\infty}\sigma_N f(y/N)=L_-+(L_+-L_-)\int_{-\infty}^y k(v)\,dv .
$$

:::example[段差の中点と単調な立ち上がり]
左で $0$、右で $2$ の波なら、拡大極限は $2\int_{-\infty}^y k(v)\,dv$ である。原点の値は偶性から $1$、原点の傾きは $2k(0)=1/\pi$ である。$k$ は非負で零点が離散的なので、極限波形は狭義増加し、どの有限の $y$ でも $0$ と $2$ の間にある。一方、対称部分和の拡大極限は $\sin y/y$ を積分したものであり、導関数の符号が変わる。この差は平均に使った核の符号に対応する。
:::

## Abel 和と Poisson 核

Cesàro 平均では次数に有限の切断があった。$0\le r<1$ で各周波数に $r^{|n|}$ を掛けると、全周波数を含んだまま高い周波数を抑えられる。特に $r$ を円板の半径と読むと各モードは調和関数になる。

:::definition[Abel 平均]{#def-four-abel}
$A_rf(x)=\sum_{n\in\Z}r^{|n|}\widehat f(n)e^{inx}$ を :term[Abel 平均]{en="Abel mean" yomi="あーべるへいきん"} とする。$r<1$ では係数の有界性により一様絶対収束する。
:::

:::proposition[Poisson 核]{#prop-four-circle-poisson}
$$
P_r(t)=\sum_{n\in\Z}r^{|n|}e^{int}=\frac{1-r^2}{1-2r\cos t+r^2}.
$$

これは $r\uparrow1$ の良い核であり、$A_rf=f*P_r$ である。さらに $u(re^{ix})=A_rf(x)$ と置くと $u$ は円板内で調和的である。
:::

:::proof
正負の等比級数を足すと商を得る。分母は $(1-r)^2+2r(1-\cos t)>0$ なので $P_r>0$ である。項別積分で質量は $1$ となる。$|t|\ge\delta$ では $r\ge1/2$ のとき分母が $1-\cos\delta>0$ 以上なので外側の質量はゼロに近づく。級数を積分の中に入れることで畳み込み表示が従う。$r^{|n|}e^{inx}$ は $n\ge0$ なら $z^n$、$n<0$ なら $\bar z^{|n|}$ である。任意の小円板 $|z|\le r_0<1$ で二階導関数の級数も $\sum_n |n|^2r_0^{|n|-2}\|f\|_1$ により一様収束する。各項が調和的なので和も調和的である。
:::

:::example[境界データから円板内の多項式を作る]
$f(t)=2+3\cos(2t)-\sin(3t)$ なら $u(re^{it})=2+3r^2\cos(2t)-r^3\sin(3t)$ である。$z=i/2$ では $u=2-3/4+1/8=11/8$ である。最大値原理により、連続な境界値を取る調和関数は高々一つである。従って :ref[thm-four-good-kernel] と合わせると連続データの円板 Dirichlet 問題を解ける。
:::

:::proposition[二つの総和と通常の極限]{#prop-four-summability-consistency}
固定した $x$ で $S_N f(x)\to L$ なら $\sigma_N f(x),A_rf(x)\to L$ である。さらに $\sigma_N f(x)\to L$ なら $A_rf(x)\to L$ である。
:::

:::proof
通常の収束から Cesàro 収束は有限個の初項と小さい誤差の尾に和を分ければ従う。$A_rf=(1-r)\sum_{j\ge0}r^jS_j f$ は $r^{|n|}=(1-r)\sum_{j\ge|n|}r^j$ から従う。さらに $S_j=(j+1)\sigma_{j+1}-j\sigma_j$ を代入すると

$$
A_rf=(1-r)^2\sum_{N\ge1}Nr^{N-1}\sigma_N f.
$$

重みの和は両表示とも $1$ で、各有限部分の重み和は $r\uparrow1$ でゼロになる。収束する列の誤差を有限部と尾部に分ければ各結論を得る。和の交換は $|S_jf|\le(2j+1)\|f\|_1$ で正当化できる。
:::

$P_r*P_s=P_{rs}$ も成り立つ。係数が $(rs)^{|n|}$ になるため、:ref[thm-four-series-uniqueness] と連続性で従う。例えば $r=s=1/2$ の二回の平滑化は $r=1/4$ の一回と同じである。

:::example[Poisson 核の反対側の値]
$r=1/3$ とする。$P_r(0)=(1+r)/(1-r)=2$、$P_r(\pi)=(1-r)/(1+r)=1/2$ である。畳み込みを周波数で見ると、定数の重みは $1$、第 $2$ モードの重みは $1/9$ である。核の点の値と周波数の乗数を混同すると、この振幅の評価を誤る。
:::

### Abel 平均から通常の和を取り戻す条件

:ref[prop-four-summability-consistency] は通常の収束から総和法への向きを与える。逆向きには係数が十分小さいという追加条件が必要である。各項の大きさが $1/n$ より速く減るとき、指数で重み付けした和と切断した和との差を積分によらず比較できる。

:::theorem[Tauber の十分条件]{#thm-four-abel-tauber}
複素数列 $a_n$（$n\ge0$）が $na_n\to0$ を満たし、$r\uparrow1$ で $\sum_{n\ge0}a_nr^n\to L$ ならば $\sum_{n=0}^Na_n\to L$ である。
:::

:::proof
$N\ge2$ に対して $r_N=1-1/N$ と置く。差を
$$
\sum_{n=0}^Na_n-\sum_{n\ge0}a_nr_N^n
=\sum_{n\le N}(1-r_N^n)a_n-\sum_{n>N}r_N^na_n
$$
と書く。$\eps>0$ を固定し、$n\ge n_0$ で $|a_n|\le\eps/n$ とする。前半の $n<n_0$ の有限個の項は $r_N\uparrow1$ で消える。残りの前半は $1-r_N^n\le n/N$ から $\eps$ 以下である。後半は $n>N$ なので $1/n\le1/N$ を使えば $\eps N^{-1}r_N^{N+1}/(1-r_N)\le\eps$ となる。差の上極限は $2\eps$ 以下であり、任意の $\eps$ に対する評価なのでゼロである。仮定の Abel 極限を合わせて結論を得る。
:::

:::example[Tauber 条件を満たさない交代列]
$\sum_{n\ge0}(-1)^n$ の Abel 和は $\lim_{r\uparrow1}(1+r)^{-1}=1/2$ である。一方通常部分和は $1,0$ と交代する。ここでは $na_n$ はゼロに近づかない。Fourier 級数の一地点に適用する場合は $a_0=\widehat f(0)$、$a_n=\widehat f(n)e^{inx}+\widehat f(-n)e^{-inx}$ と置く。各組の打ち消しも条件に使える。
:::

### 円板で弧が持つ重み

連続な定数データを分割する代わりに、一つの弧を $1$、残りを $0$ とした境界データを考える。$x=0$、$0<a<\pi$ に対して $\int_{-a}^aP_r(t)\,dm(t)$ を計算するには、$s=\tan(t/2)$ と置く。$dt=2\,ds/(1+s^2)$、$1-2r\cos t+r^2=((1-r)^2+(1+r)^2s^2)/(1+s^2)$ なので、

$$
\int_{-a}^aP_r(t)\,dm(t)
=\frac2\pi\arctan\left(\frac{1+r}{1-r}\tan\frac a2\right)
$$

となる。これは円板上の点 $r$ における弧の :term[調和測度]{en="harmonic measure" yomi="ちょうわそくど"} の具体形である。$r=0$ なら値は $a/\pi$ であり、中心からの距離に向きを指定しないため、円周の長さの比だけが現れる。$r\uparrow1$ では $a>0$ に対し値は $1$ へ近づく。

:::example[半円を高く保つ境界]
境界で $|\theta|<\pi/2$ の温度を $3$、残りを $1$ に保つ。円板の定常解は $1+2\int_{-\pi/2}^{\pi/2}P_r(x-t)\,dm(t)$ である。半径 $r=1/2$、偏角 $x=0$ の点 では $1+(4/\pi)\arctan 3$ となり、反対側 $x=\pi$ では $1+(4/\pi)\arctan(1/3)$ である。後者は、$\pi$ 回転で半円と補集合を入れ替え、両点の値の和が $4$ となることからも従う。境界の二つの端点では偶性により両側の重みが半分ずつで、各半径上で温度は $2$ になる。
:::

### Abel 平均の各点収束と半群

$1/2\le r<1$ なら $\eps=1-r$ と置いて $P_r(t)\le C\eps/(\eps^2+t^2)$（$|t|\le\pi$）である。分子は $1-r^2\le2\eps$、分母は $\eps^2+4r\sin^2(t/2)$ である。$x$ が $f\in L^1$ の Lebesgue 点なら $h(t)=|f(x-t)-f(x)|$ の原点付近の積分に $\int_{|t|\le a}h\,dt\le2\eta a$ を使える。$|t|\le\eps$ の寄与は $C\eta$、$2^j\eps<|t|\le2^{j+1}\eps$ の寄与は $C\eta 2^{-j}$ 以下になる。$\delta$ を固定してその外側を $C\eps\delta^{-2}\|h\|_1$ で抑え、最後に $\eta\downarrow0$ とすれば $A_rf(x)\to f(x)$ である。Lebesgue 点は a.e. 存在するので、この議論で a.e. 収束も得られる。

半径の代わりに $r=e^{-s}$ と書くと $s>0$ で $A_{e^{-s}}A_{e^{-t}}=A_{e^{-(s+t)}}$ となる。$B_s=A_{e^{-s}}$ と書く。これは :term[半群]{en="semigroup" yomi="はんぐん"} であり、$s=0$ に恒等写像を置くと $L^p$（$p<\infty$）と $C(\mathbb T)$ でノルム連続な初期値を取る。三角多項式 $Q$ に対し $|D|Q=\sum |n|\widehat Q(n)e^{int}$ と置くと

$$
\partial_s B_sQ=-|D|B_sQ,\qquad
\partial_s^2 B_sQ+\partial_x^2 B_sQ=0 .
$$

初めの式から半群の一次の発展方程式を、次の式から円柱座標 $(s,x)$ における調和性を読める。一般の $f\in L^1$ でも $s$ が正の下限を持てば、$|n|^k e^{-s|n|}$ の総和が有限なので任意の階で項別微分できる。

:::example[熱半群と Poisson 半群の比較]
$f(t)=\cos t+2\cos(3t)$ に対し $s=\log 2$ では $B_s f=\frac12\cos t+\frac14\cos(3t)$ である。時間 $s$ の熱半群なら、同じ $f$ への出力は $\frac12\cos t+\frac1{256}\cos(3t)$ となる。熱では指数が $-sn^2$、Poisson では $-s|n|$ である。両方とも平均を保ち、正の時間で滑らかにするが、第 $3$ 項の減衰を計算するにはどちらの半群か指定する必要がある。
:::

### 通常和を戻す別の十分条件

Cesàro 総和でも、項が小さいと第 $N$ 部分和と平均との差を制御できる。$s_N=\sum_{j=0}^N a_j$ とすれば

$$
s_N-\frac1{N+1}\sum_{k=0}^Ns_k=\frac1{N+1}\sum_{j=1}^Nj a_j .
$$

これは $a_j$ が平均に $N+1-j$ 回現れることから従う。よって $ja_j\to0$ なら Cesàro の極限から通常の極限が従う。$a_j\ge0$ の場合は減衰条件なしで Abel 極限から戻せる。実際、$0<r<1$ に対して各 $a_jr^j$ は $r\uparrow1$ で単調に $a_j$ へ増える。数え上げ測度に単調収束を使えば Abel 極限は $\sum a_j$ そのものである。このように総和法からの逆の推論には、符号の情報または項の減衰が必要になる。

### Littlewood の Tauber 定理

$na_n\to0$ を $a_n=O(1/n)$ に弱めると、単一の指数重みと切断和の差を絶対値だけで小さくすることはできない。そこで異なる $r$ の Abel 平均を有限個組み合わせ、切断の重み自体に近づける。係数の和は $1$ に保つので、組み合わせた平均も同じ Abel 極限を持つ。

:::theorem[Littlewood の Tauber 条件]{#thm-four-littlewood-tauber}
複素数列 $a_n$ が $|a_n|\le C/n$（$n\ge1$）を満たし、$\lim_{r\uparrow1}\sum_{n=0}^\infty a_n r^n=L$ ならば $\lim_{N\to\infty}\sum_{n=0}^N a_n=L$ である。
:::

:::proof
$g(x)=1_{[e^{-1},1]}(x)$ とする。$g(e^{-n/N})=1$ は $n\le N$ と同値である。$0<\delta<1/2$ に対し $g_\delta$ を、$x\le e^{-(1+\delta)}$ で $0$、$x\ge e^{-(1-\delta)}$ で $1$ とし、間を線形につなぐ連続関数とする。$0\le g_\delta\le1$ で、$g$ との差は $(1-\delta)N\le n\le (1+\delta)N$ の添字に限る。よって

$$
\limsup_{N\to\infty}\sum_{n\ge1}\frac{|g(e^{-n/N})-g_\delta(e^{-n/N})|}{n}
\le\frac{2\delta}{1-\delta}.
$$

次に $H_\delta(x)=(g_\delta(x)-x)/(x(1-x))$ を $0<x<1$ で定める。両端に連続に延長できる。なぜなら $x=0$ の近くでは $-1/(1-x)$、$x=1$ の近くでは $1/x$ に一致するからである。$\eta>0$ を固定し、$H_\delta$ を一様誤差 $\eta$ 以下の多項式 $Q$ で近似する。この近似は次の二項分布による構成で得られる。任意の連続な $H$ に対し $\sum_{j=0}^m H(j/m)\binom mj x^j(1-x)^{m-j}$ を取る。重みの和は $1$、$j/m$ の重み付き平均は $x$、分散は $x(1-x)/m$ である。これらは二項定理と $j\binom mj=m\binom{m-1}{j-1}$ および $j(j-1)\binom mj=m(m-1)\binom{m-2}{j-2}$ から従う。$|j/m-x|<h$ の項の誤差は一様連続度で抑え、外側の重みの和は分散を $h^2$ で割って抑える。$h$ を選び、次に $m$ を増やせば一様近似が得られる。

$P(x)=x+x(1-x)Q(x)$ と置く。$P(0)=0,P(1)=1$ で、$|g_\delta(x)-P(x)|\le\eta x(1-x)$ である。従って

$$
\begin{aligned}
\sum_{n\ge1}\frac{|g_\delta(e^{-n/N})-P(e^{-n/N})|}{n}
&\le\eta\sum_{n\ge1}\frac{e^{-n/N}-e^{-2n/N}}{n}\\
&=\eta\log(1+e^{-1/N})\le\eta\log2 .
\end{aligned}
$$

ここでは $-\log(1-z)=\sum_{n\ge1}z^n/n$ を $z=e^{-1/N}$ と $z=e^{-2/N}$ に使った。$P(x)=\sum_{k=1}^d c_kx^k$ かつ $\sum c_k=1$ なので、$\sum_{n\ge0}a_n P(e^{-n/N})=\sum_{k=1}^dc_k\sum_{n\ge0}a_ne^{-kn/N}\to L$ である。すべての内側の級数は $|a_n|\le C/n$ により絶対収束する。したがって切断和と $L$ との差の上極限は $C(2\delta/(1-\delta)+\eta\log2)$ 以下である。$\eta$ と $\delta$ を任意に小さくして結論を得る。
:::

:::example[滑らかさを仮定せず跳躍点の通常和へ移る]
$f$ が円周上で有界変動なら、$a_n=\widehat f(n)e^{inx}+\widehat f(-n)e^{-inx}$ は $O(1/n)$ である。Abel 平均は有界な関数の片側極限がある点で $(f(x+)+f(x-))/2$ へ近づく。偶核の片側質量が $1/2$ で、外側質量が消えることを使えばよい。Littlewood の定理を適用すれば、例えば左から $-2$、右から $3$ に近づく有界変動な入力で $S_N f(x)\to1/2$ となる。これは Dirichlet–Jordan の定理を総和法から得る別の道筋である。
:::

### Poisson 平均の連続度による誤差

$0<\alpha<1$、$\omega_f(t)\le A t^\alpha$ なら、$\eps=1-r$ に対し $\|A_r f-f\|_\infty\le C_\alpha A\eps^\alpha$ である。実際、$\eps\le1/2$ では $P_r(t)\le C\min(\eps^{-1},\eps/t^2)$ を使い、

$$
\int_0^\pi t^\alpha P_r(t)\,dt
\le C\left(\eps^{-1}\int_0^\eps t^\alpha\,dt+\eps\int_\eps^\pi t^{\alpha-2}\,dt\right)
\le C_\alpha\eps^\alpha .
$$

$\alpha=1$ では後半の積分が $\log(\pi/\eps)$ となる。例えば $f(t)=|t|$ では原点で $A_rf(0)=\frac4\pi\sum_{n\ge1,\ n\ \mathrm{odd}}(1-r^n)/n^2$ である。$1\le n\le1/\eps$ なら $1-r^n\ge1-e^{-n\eps}\ge n\eps/2$ なので、下界にも $\eps\log(1/\eps)$ が現れる。上界の対数はこの例で必要である。一方、$f(t)=\cos(kt)$ の誤差は $1-r^k\le k\eps$ であり、単一モードでは対数は現れない。

## Fejér の定理と三角多項式による近似

:::theorem[Fejér の定理]{#thm-four-fejer}
連続な周期関数は $\sigma_N f$ に一様近似される。また $1\le p<\infty$ で $f\in L^p$ は $\sigma_N f$ に $L^p$ 近似される。
:::

:::proof
:ref[prop-four-fejer-kernel] と :ref[thm-four-good-kernel] を組み合わせる。$\sigma_N f$ は次数 $N-1$ 以下の三角多項式なので、同時に各空間で三角多項式が稠密である。
:::

:::proposition[Fejér 近似の誤差]{#prop-four-fejer-rate}
$\omega_f(t)=\sup_{x,|h|\le t}|f(x+h)-f(x)|$ とする。$0<\alpha<1$ で $\omega_f(t)\le A t^\alpha$ なら $\|\sigma_N f-f\|_\infty\le C_\alpha A N^{-\alpha}$ である。$\alpha=1$ なら上界は $CA(1+\log N)/N$ となる。
:::

:::proof
$\int \omega_f(|t|)F_N(t)\,dm$ を $0<t<1/N$ と $1/N\le t\le\pi$ に分ける。最初の区間で $F_N\le N$、後半で $\sin(t/2)\ge t/\pi$ より $F_N\le\pi^2/(Nt^2)$ である。よって上界は定数を除いて

$$
A N\int_0^{1/N}t^\alpha\,dt+\frac{A}{N}\int_{1/N}^{\pi}t^{\alpha-2}\,dt
$$

となる。初項は $AN^{-\alpha}/(1+\alpha)$ であり、後半は $\alpha<1$ なら $AN^{-\alpha}/(1-\alpha)$ 以下、$\alpha=1$ なら $A\log(\pi N)/N$ である。
:::

:::example[三角多項式に対する誤差]
$f=\sin(3t)$、$N>3$ なら $\sigma_N f-f=-(3/N)f$ であり、上限ノルムの誤差は $3/N$ となる。角のない一つのモードでは、核全体の絶対値で評価した対数損失は現れない。
:::

:::proposition[平頂な近似核]{#prop-four-de-la-vallee}
$V_N=2F_{2N}-F_N$ は $\|V_N\|_1\le3$ を満たし、$|n|\le N$ で係数 $1$、$|n|\ge2N$ で係数 $0$ を持つ。次数 $N$ 以下の三角多項式の集合を $\mathcal P_N$ とすると

$$
\|f-f*V_N\|_p\le4\inf_{Q\in\mathcal P_N}\|f-Q\|_p
$$

が $1\le p\le\infty$ で成り立つ。
:::

:::proof
係数表示を差し引けば主張の係数を得る。$Q*V_N=Q$ なので誤差は $(f-Q)-(f-Q)*V_N$ である。:ref[prop-four-circle-product] と三角不等式で評価し、下限を取る。
:::

:::example[既知の周波数を保存する]
$N=3$ に対して $\cos(3t)*V_3=\cos(3t)$ である。一方、$\cos(5t)*V_3=(1/3)\cos(5t)$、$\cos(6t)*V_3=0$ である。境界 $|n|=N$ のモードはちょうど保存される。
:::

### Lebesgue 点での Fejér 収束

連続性がなくても、平均した差が半径とともに小さくなる点では先ほどの核を使える。これは一点の値の連続性を課すより弱い条件である。

:::proposition[Fejér 平均の各点収束]{#prop-four-fejer-lebesgue}
$f\in L^1(\mathbb T)$ とし、$x$ が $f$ の Lebesgue 点、つまり $\frac1{2r}\int_{-r}^{r}|f(x-t)-f(x)|dt\to0$ となる点とする。このとき $\sigma_N f(x)\to f(x)$ である。
:::

:::proof
$h(t)=|f(x-t)-f(x)|$ とする。$\eps>0$ に対し、$r\le2\delta$ で $\int_{-r}^{r}h\le 2r\eps$ となる $\delta>0$ を選ぶ。$|t|\le1/N$ では $F_N\le N$ なので積分は $C\eps$ 以下である。$2^j/N<|t|\le2^{j+1}/N$ かつ $|t|<\delta$ では $F_N\le\pi^2N2^{-2j}$ で、$h$ の積分は外側の球の評価で $C\eps 2^j/N$ 以下になる。外半径が $\delta$ をまたぐ最後の環も $2\delta$ 以内である。従って全内側の積分は $C'\eps$ 以下である。外側では $F_N$ の上限が $O_\delta(N^{-1})$ なので $h\in L^1$ から消える。最後に $\eps\downarrow0$ とする。
:::

:::example[一点の値を変えた連続関数]
$f$ が $\sin t$ と原点以外で一致し、$f(0)=1$ とする。$0$ は指定された値に関する Lebesgue 点ではない。一方、積分は常に $\sin t$ の値に基づくので、原点の Fejér 平均は任意の次数でゼロになる。可積分関数の代表元では、Lebesgue 点の条件も選んだ値で述べる必要がある。
:::

### Jackson 核と Lipschitz 近似

Fejér 核そのものでは Lipschitz 誤差の上界に対数が現れた。核を二乗すると外側での減衰が速くなり、低い周波数を近似する際にその対数を除ける。

:::proposition[Jackson の一様近似]{#prop-four-jackson}
$J_N=F_N^2/\int F_N^2dm$ とする。$0<\alpha\le1$、$\omega_f(t)\le A t^\alpha$ なら $\|f*J_N-f\|_\infty\le C_\alpha A N^{-\alpha}$ である。$f*J_N$ の次数は $2N-2$ 以下である。
:::

:::proof
直交性で $\int F_N^2dm=1+2N^{-2}\sum_{k=1}^{N-1}k^2=(2N^2+1)/(3N)$ である。この値は $2N/3$ 以上である。$F_N\le N$ と $F_N\le\pi^2/(Nt^2)$ から
$$
0\le J_N(t)\le C\min(N,N^{-3}|t|^{-4})
$$
となる。$|t|=1/N$ で分けると $\int|t|^\alpha J_N(t)dm\le C_\alpha N^{-\alpha}$ を得る。ここでは後半の積分が $N^{-3}\int_{1/N}^\pi t^{\alpha-4}dt$ で、$\alpha\le1$ なら下端からそのべきが得られる。畳み込みの差を連続度で抑えればよい。$F_N$ の次数が $N-1$ なので二乗の次数は高々 $2N-2$ である。
:::

:::example[$J_2$ の係数]
$F_2=1+\cos t$ なので $F_2^2=3/2+2\cos t+(1/2)\cos(2t)$ となる。従って $J_2=1+(4/3)\cos t+(1/3)\cos(2t)$ であり、$\cos t*J_2=(2/3)\cos t$ である。$F_2$ では同じモードの乗数は $1/2$ である。
:::

### 一般の連続度での Jackson 評価

$f\in C(\mathbb T)$ の連続度 $\omega_f(h)$ は増加関数であり、$h,s>0$ に対して $\omega_f(s)\le(1+s/h)\omega_f(h)$ となる。実際、長さが高々 $s$ の弧を $\lceil s/h\rceil$ 個以内に分割し、隣接点の値の差を足せばよい。従って $J_N$ の質量が $1$、一次の絶対モーメントが $C/N$ 以下であることから

$$
\|f-f*J_N\|_\infty\le C'\omega_f(1/N)
$$

を得る。この結果では特定の Hölder 指数を仮定していない。$E_M(f)=\inf_{Q\in\mathcal P_M}\|f-Q\|_\infty$ と置き $N=\lfloor M/2\rfloor+1$ を選ぶと、$E_M(f)\le C''\omega_f(1/M)$ である。$M=1$ も含めて、尺度の比が有界なことと連続度の分割の評価を使った。

### 高階の差分を使う近似

一次の連続度には、定数でない滑らかな関数でも一般に $h$ より速い上界を置けない。多項式近似の誤差を高いべきで抑えるには、定数だけでなく低い次数の Taylor 項を消す差分を使う。

:::proposition[高階の Jackson 評価]{#prop-four-higher-jackson}
整数 $k\ge1$ と $N\ge1$ に対し、$f\in C^k(\mathbb T)$ なら、次数が高々 $(k+1)(N-1)$ の三角多項式 $Q_N$ で

$$
\|Q_N-f\|_\infty\le C_kN^{-k}\|f^{(k)}\|_\infty
$$

を満たすものが存在する。
:::

:::proof
$J=F_N^{k+1}/\int F_N^{k+1}\,dm$ と置く。$|t|\le1/N$ では $\sin(Nt/2)$ の絶対値が $N|t|/\pi$ 以上、$\sin(t/2)$ の絶対値が $|t|/2$ 以下なので、$F_N\ge (4/\pi^2)N$ である。従って分母は $c_kN^k$ 以上である。さらに $F_N\le\min(N,\pi^2/(Nt^2))$ を使うと $J\le C_k\min(N,N^{-2k-1}|t|^{-2k-2})$ となる。$|t|=1/N$ で分割して積分すると $\int|t|^kJ\,dm\le C'_kN^{-k}$ である。

$Q_N(x)=\int J(t)\sum_{j=1}^k (-1)^{j-1}\binom{k}{j} f(x+jt)\,dm(t)$ と置く。各 $j$ に対する積分の第 $n$ 係数は $\widehat f(n)\widehat J(-jn)$ である。$J$ の次数は $(k+1)(N-1)$ 以下なので、$|n|$ がそれを超えるとこれらの係数はすべてゼロである。連続関数の有限帯域の復元により $Q_N$ は所定の多項式である。一方、$f-Q_N$ は $\Delta_t^kf(x)=\sum_{j=0}^k(-1)^{k-j}\binom{k}{j} f(x+jt)$ の積分に符号を掛けたものである。基本定理を $k$ 回使えば $\Delta_t^kf(x)=t^k\int_{[0,1]^k}f^{(k)}(x+t(s_1+\cdots+s_k))\,ds$ となる。絶対値を取って $J$ の $k$ 次モーメントで積分すれば評価が従う。
:::

:::example[同じ近似次数で一階と二階の情報を使う]
$f(x)=\cos x+\cos(4x)/16$ では $\|f'\|_\infty\le5/4$、$\|f''\|_\infty\le 2$ である。$k=2$ の構成では $Q_N=\int J(t)(2f(x+t)-f(x+2t))\,dm(t)$ となる。第 $n$ 乗数は $2\widehat J(n)-\widehat J(2n)$ である。$t$ が小さいと $2e^{int}-e^{2int}=1+n^2t^2+O(|nt|^3)$ となり、定数は保存され、一次の誤差は打ち消される。例えば $N=2$ では $J=(1+\cos t)^3/\int(1+\cos t)^3\,dm$ で、その平均は分子で $5/2$ である。$\widehat J(1)=3/4$、$\widehat J(2)=3/10$ を用いると $Q_2=(6/5)\cos x$ になる。差分を使った近似では核の正値性があっても、出力全体が入力の値域に入るとは限らない。
:::

### 近似の速さから連続度を戻す

低い次数で正確に近似できるとき、細かい尺度での $f$ の差を、近似多項式の導関数から抑える。まず次数と導関数を結ぶ評価を得る。

:::lemma[三角多項式の導関数の次数評価]{#lem-four-trig-bernstein}
整数 $N\ge1$ と $Q\in\mathcal P_N$ に対して $\|Q'\|_\infty\le C N\|Q\|_\infty$ が絶対定数 $C$ で成り立つ。
:::

:::proof
$Q=Q*V_N$ なので $Q'=Q*V_N'$ である。$\|F_M'\|_1\le CM$ を示せばよい。$|t|\le1/M$ では係数から $|F_M'|\le\sum_{|n|<M}|n|\le M^2$ であり、その区間上の積分は $O(M)$ である。外側で正弦の商を微分する。分子の導関数は絶対値が $M$ 以下で、その項は $C|t|^{-2}$ に抑えられる。分母の微分からの項は $CM^{-1}|t|^{-3}$ 以下である。いずれも $[1/M,\pi]$ 上の積分は $O(M)$ である。$V_N'=2F_{2N}'-F_N'$ の $L^1$ 評価と畳み込みの不等式を合わせる。
:::

:::proposition[近似誤差による Hölder 連続性]{#prop-four-jackson-converse}
$0<\alpha<1$ とする。$f\in C(\mathbb T)$ で $E_N(f)\le A N^{-\alpha}$ が全 $N\ge1$ について成り立つなら、$f$ は $\alpha$-Hölder 連続である。
:::

:::proof
$2^j$ 次以下の $Q_j$ で $\|f-Q_j\|_\infty\le2 A2^{-\alpha j}$ を選ぶ。$A=0$ なら有限次元の $\mathcal P_1$ が閉じているため $f\in\mathcal P_1$ である。$A>0$ では誤差の上界より、$R_j=Q_j-Q_{j-1}$ は $\|R_j\|_\infty\le C_\alpha A 2^{-\alpha j}$ を満たす。$f=Q_0+\sum_{j\ge1}R_j$ は一様収束である。$2^{-J-1}<|h|\le2^{-J}$ に対し $j\le J$ では :ref[lem-four-trig-bernstein] から差を $C |h|2^j\|R_j\|_\infty$ で抑える。その和は $C_\alpha A|h|2^{(1-\alpha)J}\le C_\alpha A|h|^\alpha$ である。$j>J$ では $2\|R_j\|_\infty$ を足せば $C_\alpha A2^{-\alpha J}\le C'_\alpha A|h|^\alpha$ になる。$Q_0$ は $O(|h|)$ の連続度を持つので結論が従う。
:::

この証明で $\alpha=1$ とすると、低周波の各 $j$ の上界が $C A|h|$ となり、個数 $J$ の因子が残る。Lipschitz 連続性と次数 $N^{-1}$ の近似誤差を同値と考えると、この端点で違いが出る。例えば $f(t)=\sum_{j=0}^\infty 2^{-j}\sin(2^jt)$ では尾の総和によって $E_N(f)\le2/N$ である。$h=2^{-J}$ の場合、$j\le J$ で $\sin(2^jh)\ge (2^jh)/2$ が成り立つので、$f(h)/h\ge(J+1)/2-1$ となる。これは $j>J$ の尾の絶対値が $2^{-J}=h$ 以下であることを使った。従ってこの連続な $f$ は原点で Lipschitz 条件を満たさない。

### Lipschitz 関数に対する Fejér 近似の対数損失

Fejér 核の上界に現れた対数は、実際に達成される。$f(t)=|t|$ を周期延長し、$x=0$ で評価する。係数は $\widehat f(n)=((-1)^n-1)/(\pi n^2)$（$n\ne0$）であり、絶対収束から $f(0)=\sum_n\widehat f(n)=0$ である。したがって

$$
\sigma_N f(0)=\frac4{\pi N}\sum_{\substack{1\le n<N\\ n\ \mathrm{odd}}}\frac1n+
\frac4{\pi}\sum_{\substack{n\ge N\\ n\ \mathrm{odd}}}\frac1{n^2}.
$$

第 2 項は $O(N^{-1})$、第 1 項の調和和は $\frac12\log N+O(1)$ である。後者は全添字の調和和から偶数項を引けばよい。従って $\sigma_N f(0)=(2/\pi)(\log N)/N+O(N^{-1})$ であり、Fejér 平均そのものでは Lipschitz 関数の誤差を一律に $O(N^{-1})$ にはできない。一方、次数 $N$ 以下の近似多項式の中で最適な一つに対する誤差は $E_N(f)=O(N^{-1})$ である。実際、この $f$ では通常部分和の係数尾部を直接絶対値で足しても $O(N^{-1})$ になる。

### 有限測度の Fourier 係数

有限複素 Radon 測度 $\mu$ に対して $\widehat\mu(n)=\int_{\mathbb T}e^{-int}\,d\mu(t)$ とする。$\mu=f\,dm$ なら通常の係数に一致する。これに対し点質量 $\delta_a$ の係数は $e^{-ina}$ であり、一般の有限測度には Riemann–Lebesgue の補題を適用できない。

:::proposition[有限測度の Fejér 近似と一意性]{#prop-four-fejer-measures}
$u_N(x)=\int_{\mathbb T} F_N(x-t)\,d\mu(t)$ と置くと、$u_N$ は $\widehat u_N(n)=(1-|n|/N)_+\widehat\mu(n)$ を持つ三角多項式であり、$\|u_N\|_1\le|\mu|(\mathbb T)$ である。任意の連続関数 $g$ について $\int g u_N\,dm\to\int g\,d\mu$ となる。従って $\widehat\mu$ は $\mu$ を一意に定める。
:::

:::proof
有限展開に積分を施せば $u_N$ の係数が得られる。$F_N\ge0$ と $\int F_Ndm=1$ から、$F_N(x-t)$ の $dm(x)d|\mu|(t)$ 積分は $|\mu|(\mathbb T)$ である。従ってノルム評価が成り立ち、さらに Fubini と $F_N$ の偶性によって $\int g u_N\,dm=\int (g*F_N)\,d\mu$ となる。:ref[thm-four-fejer] の一様収束を右辺に使う。二つの測度の係数が一致するならすべての $u_N$ が共通であるため、連続関数に対する積分も一致する。有限 Radon 測度の表現の一意性（:ref[thm-meas-c-zero-dual]）から測度は等しい。
:::

例えば $\mu=\delta_a$ の近似は $u_N=F_N(\cdot-a)$ である。その $L^1$ ノルムは常に $1$ だが、$L^1$ 極限があれば上の一意性によって $dm$ に関する密度を持つ点質量が得られてしまう。従ってこの近似が保証する収束は連続関数で試験した積分の収束である。

## Fourier 係数の減衰と関数の滑らかさ

積分で係数を求める場合、微分するたびに周波数の因子を分離できる。逆に係数を足して関数を作る場合には、微分後の級数まで一様に制御できるかが問題になる。

:::theorem[微分による係数評価と逆方向の十分条件]{#thm-four-coeff-smooth}
整数 $k\ge1$ とする。$f,\ldots,f^{(k-1)}$ が周期的な絶対連続代表元を持ち、$f^{(k)}\in L^1$ なら $n\ne0$ で

$$
\widehat f(n)=(in)^{-k}\widehat{f^{(k)}}(n)=o(|n|^{-k}).
$$

逆に $f\in L^1$ について $\sum_n(1+|n|)^k|\widehat f(n)|<\infty$ なら $f$ は $C^k$ 級の代表元を持つ。
:::

:::proof
絶対連続関数の部分積分を繰り返す。境界項は各導関数の周期性によりゼロである。最後の係数に :ref[lem-four-series-rl] を適用する。逆方向では $0\le j\le k$ の各導関数の級数に Weierstrass の判定法を使う。有限和の $j-1$ 階導関数を原点から積分表示し、一様収束で積分と極限を交換すれば極限の微分が第 $j$ の和であることが分かる。得られた関数の係数は一様収束から $f$ と一致するので :ref[thm-four-series-uniqueness] で同じ代表元となる。
:::

:::proposition[Hölder 条件と係数]{#prop-four-holder-coeff}
$f\in L^1(\mathbb T)$ とする。$0<\alpha\le1$ で $|f(x+h)-f(x)|\le A|h|^\alpha$ が全点で成り立てば $|\widehat f(n)|\le A\pi^\alpha/(2|n|^\alpha)$（$n\ne0$）である。逆に $0<\alpha<1$ で全ての $n\ne0$ に対して $|\widehat f(n)|\le A|n|^{-1-\alpha}$ なら $f$ は $\alpha$-Hölder 連続な代表元を持つ。
:::

:::proof
先の方向では $h=\pi/n$ とすると移動で係数が $-\widehat f(n)$ になり、差の係数を積分の絶対値で抑えられる。逆方向では係数が絶対可和なので連続代表元を Fourier 級数で表せる。差の第 $n$ 項に $|e^{inh}-1|\le\min(2,|nh|)$ を使い、$M=\lfloor1/|h|\rfloor$ で分ける。$0<|h|\le1/2$ なら

$$
|\tau_h f(x)-f(x)|
\le 2A|h|\sum_{n=1}^M n^{-\alpha}
 +4A\sum_{n>M}n^{-1-\alpha}.
$$

前半は積分比較で $C_\alpha A |h|M^{1-\alpha}$、後半は $4AM^{-\alpha}/\alpha$ 以下である。$M$ は $1/|h|$ の定数倍の範囲に入る。大きな $|h|$ では一様有界性で評価する。
:::

:::example[端点での導関数の跳躍]
$f(t)=|t|$ を周期延長すると、$n\ne0$ で $\widehat f(n)=(({-1})^n-1)/(\pi n^2)$ である。これは一度の部分積分と $\int_0^\pi\sin(nt)\,dt=(1-(-1)^n)/n$ から分かる。$f$ の級数は絶対可和だが、微分級数の絶対値の和は奇数添字で調和級数に比例して発散する。
:::

:::example[係数の公比と正則性]
$\widehat f(n)=3^{-|n|}$ なら任意の $k$ で $\sum_n(1+|n|)^k3^{-|n|}<\infty$ なので $f$ は $C^\infty$ である。実際 :ref[prop-four-circle-poisson] から $f(t)=4/(5-3\cos t)$ であり、分母は実軸上で $2$ 以上である。
:::

### 平方平均を使った絶対収束判定

各係数を別々に評価すると、$\alpha$-Hölder 連続でも $|\widehat f(n)|\le C|n|^{-\alpha}$ までしか分からない。この評価をそのまま足しても、$0<\alpha\le1$ では絶対可和性は得られない。係数の集まりを同時に見積もると、異なる周波数の直交性を利用できる。

:::proposition[Hölder 条件と絶対収束]{#prop-four-bernstein-absolute}
$\alpha>1/2$ の Hölder 連続な $f$ は絶対収束する Fourier 級数を持つ。
:::

:::proof
$2^j\le |n|<2^{j+1}$ とし、$h=\pi/2^{j+1}$ を選ぶ。この範囲で $\pi/2\le |nh|<\pi$ なので $|e^{inh}-1|\ge\sqrt2$ である。$g=f(\cdot+h)-f$ の有限個の周波数への射影を $P$ とする。直交性から $\|g\|_2^2=\|g-P\|_2^2+\sum|\widehat g(n)|^2$ である。従って Hölder 条件 $\|g\|_2\le A|h|^\alpha$ より

$$
\sum_{2^j\le |n|<2^{j+1}}|\widehat f(n)|^2
\le A^2 |h|^{2\alpha}/2 .
$$

添字の個数は $2^{j+1}$ なので Cauchy–Schwarz で同じ範囲の絶対値の和は $C_\alpha A 2^{(1/2-\alpha)j}$ 以下で、$j$ について総和可能である。定数項を足した後、一意性で絶対収束級数の和は $f$ となる。
:::

:::example[点ごとの減衰よりまとまりで見る]
$f(t)=|t|^{3/4}$ を一周期から延長すると $\alpha=3/4$ の Hölder 条件を満たす。原点をまたぐ場合は $||t|^{3/4}-|s|^{3/4}|\le|t-s|^{3/4}$、周期の継ぎ目をまたぐ場合は円周上で原点から測った距離を使えば同じ不等式になる。係数ごとの粗い評価は $|n|^{-3/4}$ で総和できないが、:ref[prop-four-bernstein-absolute] の評価では $2^j$ 付近の絶対値の和は $O(2^{-j/4})$ である。この評価は個々の係数の符号を必要としない。
:::

### 指数減衰と複素変数への延長

実変数での Taylor 展開の収束半径に対応して、周期関数では Fourier 係数の指数減衰が複素平面の帯の幅と対応する。

:::proposition[帯上の正則性と係数の指数減衰]{#prop-four-analytic-strip}
$a>0$ とする。連続な周期関数 $f$ が帯 $|\operatorname{Im}z|<a$ へ $2\pi$ 周期の正則関数 $F$ として延長できるなら、各 $b<a$ で $|\widehat f(n)|\le M_b e^{-b|n|}$ である。ただし $M_b=\max_{0\le x\le2\pi,\,|y|\le b}|F(x+iy)|$ とする。逆に各 $b<a$ についてこの型の評価があるなら、$f$ は同じ帯へ正則に延長できる。
:::

:::proof
$g(z)=F(z)e^{-inz}$ と置く。$I(y)=\int_0^{2\pi}g(x+iy)\,dx$ の導関数は $i\int_0^{2\pi}\partial_xg(x+iy)\,dx=0$ である。微分はコンパクトな長方形で正当化され、最後は周期性を使う。$n>0$ なら $y=-b$、$n<0$ なら $y=b$ を選ぶと、$|e^{-in(x+iy)}|=e^{-b|n|}$ である。これが前半の評価を与える。

逆に、$b'<b<a$ とすると $|\operatorname{Im}z|\le b'$ 上で $|\widehat f(n)e^{inz}|\le M_b e^{-(b-b')|n|}$ である。導関数の和も $|n|$ を掛けた同じ評価で一様収束する。分割した線分上で有限和の微積分の基本定理を使って極限を取れば、この和は正則である。実軸上では絶対収束し、一意性で $f$ に一致する。各小さい帯を合わせて延長を得る。
:::

:::example[Poisson 核の帯の幅]
$P_{1/2}(t)=3/(5-4\cos t)$ は係数 $2^{-|n|}$ を持つため $|\operatorname{Im}z|<\log2$ で正則である。実際 $\cos(iy)=\cosh y$ であり、$\cosh(\log2)=5/4$ で分母がゼロになる。分母の導関数は $4\sin(i\log 2)\ne0$ で、分子は $3$ なので両側に単純極がある。従って両側の極を含む広い帯へは延長できない。係数の公比 $1/2$ と帯の幅 $\log2$ は同じ情報を別の形で表す。
:::

### 滑らかさの逆向きが失敗する端点

$|\widehat f(n)|=O(|n|^{-2})$ は Lipschitz 条件の十分条件ではない。例えば $f(t)=\sum_{n=1}^\infty\sin(nt)/n^2$ は絶対収束によって連続で、$|\widehat f(\pm n)|=1/(2n^2)$ である。$0<h<1/2$、$M=\lfloor h^{-1}\rfloor$ とすると、$n\le M$ で $\sin(nh)\ge nh/2$ なので、低い項は合計で $\frac{h}{2}\sum_{n=1}^M 1/n$ 以上を与える。高い項の絶対値の和は $1/M\le2h$ である。したがって $f(h)/h$ は非有界である。一方、$|e^{inh}-1|\le\min(2,|nh|)$ と同じ分割から $\omega_f(h)\le C h(1+\log(1/h))$ である。係数条件の逆方向で $\alpha<1$ を課すのは、$\alpha=1$ のとき低い項に調和級数が現れるためである。

### 絶対収束する Fourier 級数の代数

連続関数 $f$ で $\|f\|_A=\sum_n|\widehat f(n)|<\infty$ となるものの全体を :term[Wiener 代数]{en="Wiener algebra" yomi="うぃーなーだいすう"} $A(\mathbb T)$ という。任意の $\ell^1(\Z)$ の列は一様収束する級数を定め、項別積分で元の列を係数として持つ。従って係数写像は $\ell^1$ への等長全単射であり、$\|\cdot\|_A$ は完備ノルムである。

$f,g\in A(\mathbb T)$ なら $\sum_{k,l}|\widehat f(k)\widehat g(l)|=\|f\|_A\|g\|_A$ なので、積の級数を周波数ごとにまとめて $\widehat{fg}(n)=\sum_k\widehat f(k)\widehat g(n-k)$ を得る。これを絶対値で足すと $\|fg\|_A\le\|f\|_A\|g\|_A$ となる。さらに $f$ が円周上でゼロを持たなければ、:ref[thm-fa-wiener] により $1/f$ も $A(\mathbb T)$ に属する。例えば $f=1-ae^{it}$、$|a|<1$ なら $\sum_{k\ge0}a^ke^{ikt}$ が逆数を表し、逆数の $A$ ノルムは $1/(1-|a|)$ である。

## Fourier 級数が発散する連続関数

良い核には絶対値の質量の一様有界性があった。Dirichlet 核では中心以外にも正負の振動が続き、その絶対値を積み上げると対数的に増大する。

:::proposition[Lebesgue 定数]{#prop-four-lebesgue-constant}
$N\ge2$ に対し $c\log(N+1)\le\|D_N\|_1\le C\log(N+1)$ が絶対定数 $c,C>0$ で成り立つ。$C(\mathbb T)$ 上の線形汎関数 $f\mapsto S_Nf(0)$ のノルムは $\|D_N\|_1$ である。
:::

:::proof
上界には $|D_N(t)|\le\min(2N+1,\pi/|t|)$ を使い、$1/(2N+1)$ で積分を分ける。下界では $\sin(t/2)\le t/2$ を使い、

$$
\|D_N\|_1\ge\frac2\pi\int_0^{(N+1/2)\pi}\frac{|\sin v|}{v}\,dv
\ge\frac4{\pi^2}\sum_{j=1}^{N}\frac1j.
$$

最後は各区間 $[(j-1)\pi,j\pi]$ で $1/v\ge1/(j\pi)$ と $\int|\sin v|=2$ を用いた。作用素ノルムは積分の三角不等式で上から抑えられる。下からは連続関数 $f_\eps(t)=D_N(-t)/\sqrt{D_N(-t)^2+\eps^2}$ を代入し、優収束で $S_N f_\eps(0)\to\|D_N\|_1$ を得る。
:::

:::theorem[指定点での発散]{#thm-four-continuous-divergence}
任意の $x_0\in\mathbb T$ に対し、$|S_N f(x_0)|$ が非有界となる実数値連続関数 $f$ が存在する。この性質を持つ関数の集合は $C(\mathbb T;\R)$ で稠密な $G_\delta$ である。
:::

:::proof
移動で $x_0=0$ に帰着する。各整数 $m\ge1$ で $U_m=\bigcup_N\{f:|S_N f(0)|>m\}$ は開集合である。これが稠密でないと、ある開球 $B(f_0,r)$ 上で全ての $|S_N f(0)|$ が $m$ 以下となる。$\|h\|_\infty\le r/2$ に対し $f_0+h,f_0$ の差を取り $|S_N h(0)|\le2m$ を得る。従って各汎関数のノルムが $4m/r$ 以下となり :ref[prop-four-lebesgue-constant] に矛盾する。完備空間の Baire の範疇定理により $\bigcap_m U_m$ は稠密であり、これは求める集合である。Baire の範疇定理そのものは :ref[thm-fa-baire]、一様有界性原理による発散の導出は :ref[prop-fa-fourier-divergence] を参照。
:::

:::example[有限の次数での増幅]
$D_1=1+2\cos t$ の負の領域は $2\pi/3<|t|\le\pi$ である。そこでの積分を差し引くと

$$
\|D_1\|_1=1+\frac{2}{\pi}\left(\sqrt3-\frac{\pi}{3}\right)>1.
$$

従って上限ノルム $1$ の連続関数の第 $1$ 部分和でも原点の値が $1$ を超え得る。一方、同じ制約で $\sigma_N f(0)$ や $A_rf(0)$ は絶対値 $1$ 以下である。
:::

### 連続度による通常部分和の一様収束

:::corollary[Hölder 連続な関数の通常部分和]{#cor-four-holder-uniform-series}
ある $\alpha>0$ で $f$ が $\alpha$-Hölder 連続なら $S_N f\to f$ は一様である。また $0<\alpha\le1$ では誤差は $O(N^{-\alpha}\log(N+1))$ である。
:::

:::proof
$M=\lfloor N/2\rfloor+1$ とし $Q=f*J_M$ と置く。次数 $2M-2$ は $N$ 以下なので $S_NQ=Q$ である。差を $S_N(f-Q)-(f-Q)$ と書き、:ref[prop-four-lebesgue-constant] と :ref[prop-four-jackson] を適用する。$\alpha>1$ の仮定からも $\alpha=1$ の条件が従うので一様収束を得る。
:::

:::example[連続度と係数の総和を別々に判断する]
$x=0$ に尖点を持つ $|\sin(x/2)|^\alpha$（$0<\alpha<1$）は $\alpha$-Hölder 連続である。実際、$a,b\ge0$ で $|a^\alpha-b^\alpha|\le|a-b|^\alpha$ である。これは $a=b+h$ として $(b+h)^\alpha-b^\alpha$ が $b$ の非増加関数になることから従う。さらに正弦の Lipschitz 条件を使えばよい。この近似の議論は係数の絶対値の和を計算せずに一様収束を与える。
:::

### Lebesgue 定数の主項

上界と下界が同じ次数であるだけでなく、対数の係数も求められる。$0<t\le\pi$ で $1/\sin(t/2)=2/t+O(t)$ である。原点近くでは $\sin(t/2)=t/2+O(t^3)$ から従い、原点から離れた閉区間では比が有界である。これを Dirichlet 核の積分に使うと

$$
\|D_N\|_1=\frac2\pi \int_0^{(N+1/2)\pi}\frac{|\sin v|}{v}\,dv+O(1)
$$

である。$|\sin v|$ の周期は $\pi$、平均は $2/\pi$ である。$B(v)=\int_0^v(|\sin s|-2/\pi)\,ds$ は周期的かつ有界である。$1$ から $R$ まで部分積分して

$$
\int_1^R\frac{|\sin v|-2/\pi}{v}\,dv
=\frac{B(R)}R-B(1)+\int_1^R\frac{B(v)}{v^2}\,dv=O(1)
$$

となる。$0<v<1$ の元の積分は有限なので $\|D_N\|_1=(4/\pi^2)\log (N+1)+O(1)$ である。この係数は $dm=dt/(2\pi)$ とした正規化を含む。

### 三角多項式のブロックで発散例を作る

前の証明で得た発散例は、どの周波数がいつ誤差を大きくするかを指定しない。次の構成では周波数の間に間隔を取り、各ブロックの前半だけを切り出すときに大きな値を作る。

:::lemma[有限正弦和の一様有界性]{#lem-four-sine-block}
任意の $M\ge1$ で $\|\sum_{k=1}^M\sin(kt)/k\|_\infty\le C$ が絶対定数で成り立つ。
:::

:::proof
奇性と周期性から $0<t\le\pi$ で示せば十分である。$J=\min(M,\lfloor 1/t\rfloor)$ とする。$k\le J$ の総和は $Jt\le1$ である。残りがあれば $m=J+1$ とし、$\sum_{k=m}^{l}\sin(kt)$ の上界 $2/|1-e^{it}|=1/\sin(t/2)$ を等比級数から得る。部分和を $A_l$ と書けば $\sum_{k=m}^M\sin(kt)/k=A_M/M+\sum_{l=m}^{M-1} A_l(1/l-1/(l+1))$ である。従って絶対値は $1/(m\sin(t/2))$ 以下で、$mt>1$ と $\sin(t/2)\ge t/\pi$ により $\pi$ 以下である。
:::

:::proposition[原点で非有界になる具体的な連続関数]{#prop-four-explicit-continuous-divergence}
$M_j=2^{2^{2j}}$（$j\ge1$）とし、

$$
f(t)=\sum_{j\ge1}2^{-j}\sin(2M_jt)\sum_{k=1}^{M_j}\frac{\sin(kt)}{k}
$$

と置く。$f$ は連続で、$S_{2M_j}f(0)\to+\infty$ となる。一方 $S_{3M_j}f(0)=f(0)=0$ である。
:::

:::proof
:ref[lem-four-sine-block] によって各ブロックは一様有界であるから、外の級数は一様絶対収束する。第 $j$ ブロックの周波数は

$$
\sin(2M_jt)\sum_{k=1}^{M_j}\frac{\sin(kt)}k
=\frac12\sum_{k=1}^{M_j}\frac{\cos((2M_j-k)t)-\cos((2M_j+k)t)}k
$$

より $M_j$ 以上 $3M_j$ 以下にある。$3M_j<M_{j+1}$ なので区間は互いに交わらない。従って $S_{2M_j}f$ には過去のブロックが全体で入り、未来のブロックは入らない。過去のブロックは $t=0$ で各対が打ち消し合う。第 $j$ ブロックでは前半だけが残るため $S_{2M_j}f(0)=2^{-j-1}\sum_{k=1}^{M_j}1/k$ である。調和和は $\log(M_j+1)$ 以上であり、$2^{-j-1}\log M_j=2^{j-1}\log 2$ は無限大に向かう。$S_{3M_j}$ に切り替えると第 $j$ ブロックも全体で入るから、値はゼロである。
:::

:::example[一つの小さなブロックでの計算]
$M=3$ のブロックは
$[\cos(5t)-\cos(7t)]/2+[\cos(4t)-\cos(8t)]/4+[\cos(3t)-\cos(9t)]/6$
である。$S_6$ を取って原点に代入すると $1/2+1/4+1/6=11/12$ である。$S_9$ では各対が入り、原点の値はゼロに戻る。この小さなブロックでは有限の誤差しか作らないが、一様有界性を保ったまま調和和を大きくできることが無限の構成で効く。
:::

### 近似誤差と通常部分和の境界

任意の $Q\in\mathcal P_N$ に対し $S_NQ=Q$ なので、

$$
\|f-S_Nf\|_\infty\le(1+\|D_N\|_1) E_N(f)
$$

である。従って $E_N(f)\log N\to0$ なら一様収束する。一般の Jackson 評価を合わせると $\omega_f(1/N)\log N\to0$ でも十分である。これは :term[Dini–Lipschitz 条件]{en="Dini–Lipschitz condition" yomi="でぃにりぷしっつじょうけん"} と呼ばれる形である。積分 $\int_0^{1/2}\omega_f(t)\,dt/t$ が収束する場合にもこの条件は成り立つ。実際、$\omega_f$ は増加関数なので $\frac12\omega_f(h)\log(1/h)\le\int_h^{\sqrt h}\omega_f(t)\,dt/t$ であり、右辺は $h\downarrow0$ でゼロに近づく。判定に使うのは連続度の積分または近似誤差であり、連続という仮定だけでは対数因子を打ち消せない。

## 演習問題

::::exercise{level=basic}
$N=3$ の $F_N$ を余弦で表し、$F_3(0)$ と $F_3(\pi)$ を求めよ。
:::solution
係数は定数項 $1$、$\pm1$ で $2/3$、$\pm2$ で $1/3$ である。従って $F_3=1+(4/3)\cos t+(2/3)\cos(2t)$ である。$t=0$ では $3$、$t=\pi$ では $1/3$ を得る。
:::
::::

::::exercise{level=basic}
$f=1+2\sin(3t)$ の $A_{1/2}f$ を求めよ。
:::solution
定数は残り、$\pm3$ の係数には $(1/2)^3=1/8$ が掛かる。従って $A_{1/2}f=1+(1/4)\sin(3t)$ である。
:::
::::

::::exercise{level=basic}
$\cos(2t)*\sin(2t)$ を求めよ。
:::solution
$\cos(2t)$ の $\pm2$ 係数はいずれも $1/2$ である。$\sin(2t)$ の係数は $\pm2$ にしかないので、積を取ると全ての係数が $\sin(2t)/2$ と一致する。一意性と連続性から畳み込みは $\sin(2t)/2$ である。
:::
::::

::::exercise{level=basic}
$F_4$ と $P_{1/2}$ の畳み込みにおける $\cos(2t)$ の振幅を求めよ。
:::solution
第 $2$ 係数の積は $(1-2/4)(1/2)^2=1/8$ である。両核は実偶関数なので第 $-2$ 係数も同じである。二項をまとめると振幅は $1/4$ となる。
:::
::::

::::exercise{level=advanced}
$|n|\ge1$ で $|\widehat f(n)|\le A/n^2$ なら、連続代表元の連続度が $\omega_f(h)\le CAh(1+|\log h|)$（$0<h\le1/2$）で抑えられることを示せ。
:::solution
係数が絶対可和なので級数で連続代表元を取る。$M=\lfloor1/h\rfloor$ まで $|e^{int}-1|\le nh$、それ以降 $2$ で抑える。誤差は $2Ah\sum_{n\le M}n^{-1}+4A\sum_{n>M}n^{-2}$ 以下である。積分比較から和は $1+\log M$ と $1/M$ 以下であり、$M\ge1/(2h)$ を使う。
:::
::::

::::exercise{level=advanced}
$f\in C(\mathbb T)$ に対し $r\uparrow1$ で $\|f*P_r-f\|_\infty$ を $\omega_f(\delta)$ と $(1-r)/\delta$ で評価せよ。
:::solution
$r\ge1/2$ とする。$0<|t|\le\pi$ で $1-\cos t=2\sin^2(t/2)\ge2t^2/\pi^2$ なので $P_r(t)\le C(1-r)/t^2$ である。内側は質量 $1$ 以下なので $\omega_f(\delta)$ 以下、外側は $2\|f\|_\infty$ に $\int_\delta^\pi C(1-r)t^{-2}\,dt$ の定数倍を掛ければよい。よって和は $\omega_f(\delta)+C'\|f\|_\infty(1-r)/\delta$ 以下である。
:::
::::

::::exercise{level=advanced}
$f\ge0$、$f\in L^1$、$f\ne0$ とする。$0<r<1$ で $A_rf(x)>0$ が全ての $x$ で成り立つことを示せ。
:::solution
$P_r(t)$ の分母は $(1+r)^2$ 以下なので $P_r\ge(1-r)/(1+r)>0$ である。畳み込みを評価すると $A_rf(x)\ge(1-r)\|f\|_1/(1+r)>0$ である。
:::
::::

::::exercise{level=advanced}
可算集合 $E\subset\mathbb T$ の全点で $|S_N f(x)|$ が非有界となる連続関数が存在することを示せ。
:::solution
各 $x\in E$ に対し :ref[thm-four-continuous-divergence] の集合は開稠密集合の可算共通部分である。$E$ が可算なので全点に対する共通部分も開稠密集合の可算共通部分である。Baire の範疇定理によりこの共通部分は稠密であり、特に空でない。
:::
::::

::::exercise{level=basic}
$N=4$、$r=2$ の整数階 Cesàro 平均を $f(t)=6+5\sin t+3\cos(2t)+15\cos(4t)$ に施せ。
:::solution
分母は $A_4^{(2)}=15$、周波数 $0,1,2,4$ の分子は $15,10,6,1$ である。正負の添字は同じ乗数なので、出力は $6+\frac{10}{3}\sin t+\frac65\cos(2t)+\cos(4t)$ となる。最高次数 $4$ の乗数はゼロではなく $1/15$ である。
:::
::::

::::exercise{level=advanced}
$f$ が連続、$s>0$、$|\widehat f(n)|\le(1+|n|)^{-2}$ とする。全部分和が一様に収束することを示し、Poisson 平均の誤差 $\|B_s f-f\|_\infty$ を $s\downarrow0$ で評価せよ。
:::solution
係数が絶対可和だから級数は一様収束し、一意性で和は $f$ である。$B_s$ は第 $n$ 項に $e^{-s|n|}$ を掛ける。従って差は高々 $2\sum_{n\ge1}\min(sn,1)/(1+n)^2$ である。$0<s<1/2$ で $M=\lfloor1/s\rfloor$ と置けば、前半は $2s\sum_{n\le M}1/n\le2s(1+\log M)$、後半は $2/M\le4s$ である。従って $O(s(1+\log(1/s)))$ の誤差となる。平均項は差に現れない。
:::
::::

## 確認テスト

::::quiz{answer="4" tol="1e-6"}
$F_4(0)$ を求めよ。
:::explanation
和 $\sum_{j=0}^3e^{ij0}$ は $4$ であり、二乗して $4$ で割ると $4$ になる。
:::
::::

::::quiz
良い核による畳み込みの収束について正しいのはどれか。

- [x] $f$ が連続なら一様収束する。
- [x] $f\in L^1$ なら $L^1$ 収束する。
- [ ] $f\in L^\infty$ なら常に $L^\infty$ 収束する。

:::explanation
前二者は :ref[thm-four-good-kernel] に従う。例えば跳躍のある矩形波を Fejér 核で畳み込むと連続関数になる。跳躍幅が $2$ のとき連続関数との本質上限距離は少なくとも $1$ である。左右からの極限が連続性で同じ値になる一方、その値は両端から距離 $1$ 未満には入れないためである。
:::
::::

::::quiz{answer="3" tol="1e-6"}
$P_{1/2}(0)$ の値を求めよ。
:::explanation
$P_r(0)=(1+r)/(1-r)$ に $r=1/2$ を代入すると $3$ となる。
:::
::::

::::quiz
ある点で $\sigma_N f(x)$ が $L$ に収束する。必ず従うのはどれか。

- [ ] $S_N f(x)$ が $L$ に収束する。
- [x] $A_rf(x)$ が $L$ に収束する。

:::explanation
$A_r$ は重み $(1-r)^2Nr^{N-1}$ で $\sigma_N$ を平均する。この重みの質量は $1$ で有限部分が消えるから極限は $L$ である。連続関数には $\sigma_N$ が一様収束するものの、指定点で $S_N$ が非有界な例がある。
:::
::::

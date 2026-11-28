---
title: "超関数"
description: "試験関数上の連続線形汎関数を用いて微分、台、畳み込み、基本解を定義する。"
notation:
  - symbol: '\mathcal D(\Omega)'
    meaning: "開集合 Ω 上のコンパクト台付き滑らかな試験関数"
  - symbol: '\mathcal D''(\Omega)'
    meaning: "Ω 上の超関数"
  - symbol: '\mathcal E''(\Omega)'
    meaning: "Ω 上のコンパクト台を持つ超関数"
  - symbol: '\langle T,\phi\rangle'
    meaning: "試験関数に対する超関数の線形作用"
  - symbol: '\delta_a'
    meaning: "a で評価する Dirac のデルタ"
---

## この章の目標

- 試験関数の位相と超関数の連続性を局所評価で記述する。
- 跳躍、点質量、主値を超関数として計算する。
- 台を使って超関数の積と畳み込みの定義域を決める。
- 基本解の作り方を具体例と定数係数作用素で説明する。

## 試験関数の空間とその位相

関数を一点で評価する代わりに、狭い範囲の滑らかな重みと掛けて積分する。そこで重みの台と導関数を制御すれば、関数より特異な対象にも同じ操作を定義できる。

:::definition[試験関数と帰納極限位相]{#def-four-test-space}
$\Omega\subset\R^n$ を開集合とする。コンパクト集合 $K\subset\Omega$ に対し $\mathcal D_K=\{\phi\in C^\infty(\Omega):\operatorname{supp}\phi\subset K\}$ と置く。台は非零点集合の閉包である。$\mathcal D_K$ の位相は

$$
p_m(\phi)=\max_{|\alpha|\le m}\sup_{x\in K}|\partial^\alpha\phi(x)|\quad(m=0,1,\ldots)
$$

で定める。$\mathcal D(\Omega)=\bigcup_{K\Subset\Omega}\mathcal D_K$ を :term[試験関数の空間]{en="test function space" yomi="しけんかんすうのくうかん"} といい、各包含が連続となる最も細かい局所凸位相を入れる。
:::

$\mathcal D_K$ は距離化可能で完備である。距離は $\sum_m2^{-m-1}\min(1,p_m(\phi-\psi))$ を使える。Cauchy 列の各導関数は一様極限を持つ。小直方体で座標線に沿った微積分の基本定理を極限に移すと、各極限は下位の極限の導関数になる。関数の極限は $K$ の外でゼロなので $\mathcal D_K$ に残る。局所凸帰納極限の一般論は :ref[def-fa-inductive-limit] を参照。

:::proposition[列の収束条件]{#prop-four-test-convergence}
$\phi_j\to0$ が $\mathcal D$ で成り立つためには、ある一つのコンパクト集合 $K\Subset\Omega$ に全ての台が含まれ、各 $m$ で $p_m(\phi_j)\to0$ となることが必要十分である。
:::

:::proof
十分性は包含の連続性である。必要性のため、コンパクトな汲み尽くし $K_1\subset\operatorname{int}K_2\subset\cdots$ を取る。台が共通のコンパクト集合に入らない場合、増加する添字 $j_l$ と $x_l\notin K_l$ を $\phi_{j_l}(x_l)\ne0$ となるように選べる。有限個の初項の台を含むよう $K_l$ の添字を毎回大きくして選べばよい。すると $q(\phi)=\sup_l|\phi(x_l)|/|\phi_{j_l}(x_l)|$ は有限のセミノルムであり、各 $\mathcal D_K$ 上では有限個の項しか寄与しない。帰納極限位相の定義から $q$ は連続であるが、$q(\phi_{j_l})\ge1$ に矛盾する。共通台が得られれば、大域的な導関数の上限セミノルムも各 $\mathcal D_K$ 上連続なので、その値はゼロに近づく。
:::

例えば $\phi_j(x)=j^{-1}\psi(x)$ は固定した $\psi\in\mathcal D$ の倍数なのでゼロへ収束する。一方 $\phi_j(x)=e^{-j^2}\psi(x-j)$ は $\R$ で全ての導関数が一様にゼロへ近づくが、$\psi\ne0$ なら台が無限遠へ逃げるため $\mathcal D$ では収束しない。

任意の $K\Subset U\subset\Omega$（$U$ は開集合）について、$K$ の近傍で $1$ となる $\chi\in C_c^\infty(U)$ を選べる。実際、$K$ の小さい有界近傍の指示関数を、それよりさらに狭い台を持つ滑らかな非負の質量 $1$ の関数で畳み込む。台と境界からの距離をそれぞれ確保すると $K$ の近傍で積分の全質量が残る。このような $\chi$ を切断関数として用いる。

:::example[境界へ台が近づく列]
$\Omega=(0,1)$、空でない台を $[-1,1]$ に持つ $\psi$ を取り、$\phi_j(x)=e^{-j^2}\psi(4j(x-1/j))$ とする。$j\ge2$ なら台は $[3/(4j),5/(4j)]$ に入る。任意の固定階数の導関数の上限は $C_m(4j)^me^{-j^2}\to0$ である。しかし台は端点 $0$ へ集まるので、$\Omega$ 内の一つのコンパクト集合にすべてを入れられない。従って :ref[prop-four-test-convergence] の収束は起きない。
:::

:::example[固定台での微分の違い]
原点の近くで $1$ となる $\chi\in\mathcal D(\R)$ に対し $\phi_j(x)=j^{-1}\chi(x)\sin(jx)$ とする。台は固定でき、$\|\phi_j\|_\infty\le\|\chi\|_\infty/j$ である。一方 $\phi'_j(0)=1$ なので $\|\phi_j'\|_\infty\ge1$ である。従って $p_0$ の収束だけでは試験関数の収束を決められない。
:::

### 切断関数を具体的に作る

試験関数の存在を使う際は、どこでゼロで、どこで $1$ になるかを同時に指定することがある。切断には境界に近い場所での全ての導関数をそろえる必要があるので、一次式をつないだ関数をそのまま使うことはできない。

$h(t)=e^{-1/t}$（$t>0$）、$h(t)=0$（$t\le0$）と置く。$t>0$ での $m$ 階導関数は $P_m(1/t)e^{-1/t}$（$P_m$ は多項式）になる。実際 $u=1/t$ と書けば $d/dt=-u^2d/du$ なので、この形は微分で保たれる。任意の $j\ge0$ について $u^je^{-u}\to0$ であり、同じことは各階の導関数をさらに $t$ で割った場合にも成り立つ。このため $t=0$ で導関数が全てゼロとなる滑らかな延長である。従って $x\mapsto h(1-|x|^2)$ は台が閉単位球で、内部で正の試験関数である。積分で割れば、先に使った非負の質量 $1$ の核になる。

平らな部分を持つ切断には
$$
b(t)=\frac{h(t)}{h(t)+h(1-t)}
$$
が使える。分母は全ての実数で正であり、$b=0$（$t\le0$）、$b=1$（$t\ge1$）である。$x\mapsto b((R^2-|x|^2)/(R^2-r^2))$ は $0<r<R$ のとき、閉球 $|x|\le r$ で $1$、$|x|\ge R$ で $0$ である。これを平行移動すれば任意の中心の周囲に切断を置ける。

:::example[重ねた切断から台を制御する]
$\Omega=(-4,4)$、$K=[-2,-1]\cup[1,2]$ とする。中心を $-3/2$、$3/2$ に取ってそれぞれ $r=5/8$、$R=3/4$ の切断を作ると、第一の台は $[-9/4,-3/4]$、第二の台は $[3/4,9/4]$ に入る。二つの支持が交わらないため、その和は $0$ から $1$ の間にあり、$K$ の近傍で $1$ になる。この切断を掛ければ、$K$ の二つの成分の近くで与えられた関数の値を保ったまま、中央の区間 $(-3/4,3/4)$ でゼロにできる。
:::

### 帰納極限と有界集合

試験関数の有界性は、振幅の大小に加えて台が同じコンパクト集合に入るかどうかを含む。これは後で、超関数の連続な族に共通の評価を求めるときに使う。

:::proposition[試験関数族の有界性]{#prop-four-test-bounded-sets}
$B\subset\mathcal D(\Omega)$ が有界集合であることと、ある $K\Subset\Omega$ に対し $B\subset\mathcal D_K$ かつ全ての $m$ について $\sup_{\phi\in B} p_m(\phi)<\infty$ であることは同値である。ここで有界とは、任意の連続セミノルムが集合上で有界であることをいう。
:::

:::proof
各包含が連続なので、共通台とセミノルムの上界があれば任意の連続セミノルムが $B$ 上で有界になる。逆に共通台がないと仮定し、コンパクトな汲み尽くし $(K_j)$ を取る。各 $j$ に対し $\phi_j\in B$ と $x_j\notin K_j$ を $\phi_j(x_j)\ne0$ となるように選ぶ。
$$
q(\phi)=\sup_{j\ge1}j\frac{|\phi(x_j)|}{|\phi_j(x_j)|}
$$
は各 $\mathcal D_K$ 上で有限個の評価の最大値に等しくなるため、そこで連続である。各包含が連続な最も細かい局所凸位相は、各 $\mathcal D_K$ への制限が連続な全セミノルムで生成される。このため $q$ も連続である。しかし $q(\phi_j)\ge j$ なので有界性に反する。共通台の存在が分かった後、各 $\max_{|\alpha|\le m}\|\partial^\alpha\phi\|_\infty$ を全空間上の連続セミノルムとみなせば残りの評価も従う。
:::

:::example[振動と増幅で有界性を分ける]
固定した $\chi$ が原点近くで $1$ なら、$B=\{j^{-4}\chi(x)e^{ijx}:j\ge1\}$ は台が固定しており、$m\le4$ で導関数の上限が有界である。しかし $x=0$ での五階微分は $i^5j$ であるから、$B$ は $\mathcal D$ では有界でない。一方 $\{e^{-j^2}\chi(x)e^{ijx}:j\ge1\}$ は、任意の固定階数の評価が $C_m e^{-j^2}(1+j)^m$ 以下なので有界であり、その列はゼロへ収束する。
:::

:::example[台が縮む速さと導関数の速さ]
$\psi$ の台を $[-1,1]$ とし、$\phi_\eps(x)=\eps^2\psi(x/\eps)$ とする。$0<\eps\le1$ なら台は $[-1,1]$ に入り、$m$ 階導関数は $\eps^{2-m}\psi^{(m)}(x/\eps)$ である。$m=2$ では一般に上限がゼロへ行かず、$m=3$ では $\psi'''$ が非零なら上限が発散する。そのため単なる振幅の縮小と試験関数の収束は区別する必要がある。これを $\phi_\eps=e^{-1/\eps}\psi(x/\eps)$ に変えると各固定の $m$ で $e^{-1/\eps}\eps^{-m}\to0$ であり、共通台もあるので試験関数としてゼロに収束する。
:::

### 局所計算をつなぐ分割

$K$ が有限個の開集合 $U_i\subset\Omega$ に覆われているとする。各 $x\in K$ についてどれか一つの $U_i$ に閉包が入る小球を選べる。そこでさらに小さい球で $1$ になる非負の切断を置き、小さい球の有限部分被覆を取る。同じ $U_i$ に入る切断を足して $b_i$ と置けば、$B=\sum_i b_i$ は $K$ の近傍で正である。この領域に台を持ち $K$ の近傍で $1$ になる別の切断 $\chi$ を掛け、$\theta_i=\chi b_i/B$ をその領域の外でゼロに延長する。すると $\theta_i$ は $U_i$ にコンパクト台を持ち、$\sum_i\theta_i=1$ が $K$ の近くで成り立つ。従って $\operatorname{supp}\phi\subset K$ なら $\phi=\sum_i\theta_i\phi$ と分解できる。$\chi$ の台が $B>0$ の領域内にあるため、$B$ がゼロへ近づく位置で商を滑らかに延ばせるという仮定は不要である。

:::example[重なりを使って区間の試験を割る]
$K=[-1,1]$ を $U_1=(-2,1/2)$ と $U_2=(-1/2,2)$ で覆う。$\theta_1(x)=1-b(x+1/2)$、$\theta_2(x)=b(x+1/2)$ と置くと、$\theta_2$ は $x\le-1/2$ でゼロ、$\theta_1$ は $x\ge1/2$ でゼロである。台を開区間の内部に入れるには遷移をより狭くし、$b(2x+1/2)$ に変える。すると二つの境界は $-1/4,1/4$ に動き、$K$ に台を持つ $\phi$ について $\theta_1\phi$ の台は $[-1,1/4]$、$\theta_2\phi$ は $[-1/4,1]$ に入る。両端が開被覆の境界に触れないことが、各片をその開集合の試験関数として扱うために必要である。
:::

## 超関数の定義と例

:::definition[超関数]{#def-four-distribution}
$\mathcal D(\Omega)$ 上の連続な複素線形汎関数を :term[超関数]{en="distribution" yomi="ちょうかんすう"} と呼び、空間を $\mathcal D'(\Omega)$ と書く。作用 $\langle T,\phi\rangle$ は $\phi$ に関して線形であり、ここでは共役を取らない。
:::

:::proposition[連続性の局所評価]{#prop-four-distribution-estimate}
線形写像 $T$ が超関数であることと、各 $K\Subset\Omega$ についてある $m,C$ が存在して $|\langle T,\phi\rangle|\le Cp_m(\phi)$ が $\mathcal D_K$ 上で成り立つこととは同値である。
:::

:::proof
連続なら各 $\mathcal D_K$ 上でも連続である。$p_m$ が増加するため、$|\langle T,\phi\rangle|<1$ となる零近傍に $p_m(\phi)<\eps$ という形の近傍を含められる。斉次性で評価を得る。逆にその評価があればセミノルム $\phi\mapsto|\langle T,\phi\rangle|$ は各 $\mathcal D_K$ 上連続なので、帰納極限位相でも連続である。
:::

:::example[局所可積分関数と測度]{#exm-four-regular-distribution}
$f\in L^1_{\mathrm{loc}}$ なら $\langle T_f,\phi\rangle=\int f\phi$ で超関数を得る。局所評価は $p_0(\phi)\int_K|f|$ である。$T_f=0$ なら、任意のコンパクト集合の近くで $f$ を切断した後、滑らかな近似単位と畳み込むとその集合上で値がゼロになる。近似単位は $L^1$ で元へ戻るので $f=0$ a.e. である。従ってこれは単射である。例えば $\R$ の $|x|^{-1/2}$ は原点近くで積分でき、この方法で超関数になる。局所有限な Radon 測度 $\mu$ も $\langle\mu,\phi\rangle=\int\phi\,d\mu$、定数 $|\mu|(K)$ を用いて定義できる。
:::

:::example[Dirac のデルタと主値]{#exm-four-delta-pv}
$a\in\Omega$ に対し $\langle\delta_a,\phi\rangle=\phi(a)$ と定める。これが関数 $f$ から来るとすれば、$a$ を避けた試験関数で作用がゼロなので $f=0$ が $\Omega\setminus\{a\}$ 上でほとんど至る所成り立つ。しかし $\phi(a)=1$ なら作用は $1$ であり矛盾する。

$\R$ では $\langle\operatorname{pv}(1/x),\phi\rangle=\lim_{\eps\downarrow0}\int_{|x|>\eps}\phi(x)/x\,dx$ と定める。台を $[-R,R]$ に入れると積分は $\int_\eps^R(\phi(t)-\phi(-t))/t\,dt$ である。被積分関数は $2p_1(\phi)$ 以下で、極限と連続性の評価 $2Rp_1(\phi)$ を得る。例えば偶関数の試験関数では作用はゼロである。
:::

### 超関数の階数と測度

局所評価に現れる $m$ は、試験関数の値だけで作用が決まるか、導関数の情報も必要かを表す。ある $m$ が全てのコンパクト $K$ に対して使えるとき、超関数の :term[階数]{en="order of a distribution" yomi="かいすう"} は高々 $m$ であるという。ただし評価定数 $C_K$ は $K$ に依存できる。また各 $K$ に対してある $m_K$ が使えるのは超関数の定義そのものであり、この条件から全領域に共通の $m$ があるとは限らない。

:::proposition[階数ゼロと正の超関数]{#prop-four-order-zero-measure}
階数が高々ゼロの超関数は、局所有限な複素 Radon 測度により一意に表される。また $\phi\ge0$ の全ての実試験関数に対して $\langle T,\phi\rangle\ge0$ となる超関数は階数が高々ゼロで、表す測度は正である。
:::

:::proof
階数が高々ゼロなら、コンパクト台の連続関数を滑らかに近似し、作用を延長する。具体的には台の周りの開近傍の閉包 $K$ を $\Omega$ 内のコンパクトに取り、連続関数を $\R^n$ にゼロで延長して半径がゼロへ行く滑らかな近似単位で畳み込む。畳み込みの台は十分小さい半径で $K$ に含まれ、元の連続関数は一様連続なので一様収束する。零階の評価により作用の列は Cauchy であり、近似方法や $K$ を変えても、二つの近似の共通台を使えば極限は一致する。従って $C_c(\Omega)$ 上に、各コンパクト台で上限ノルムにより有界な作用を得る。

$\Omega$ を相対コンパクトな開集合 $V$ で汲み尽くす。各 $V$ の $C_c(V)$ では台を $\overline V$ に入れた一つの評価が使えるので、上限ノルムで閉包を取った $C_0(V)$ へ延長できる。:ref[thm-meas-c-zero-dual] により有限複素 Radon 測度が各 $V$ に得られる。重なりでは $C_c$ 上の作用が一致するため、同じ定理の一意性で測度の制限も一致する。これらを貼り合わせると局所有限な測度となる。元の試験関数に対する作用を持つ測度は一様近似によって連続試験関数への作用も同じであり、一意である。

正の超関数の場合は、$K$ 近くで $1$ となる非負切断 $\chi$ を選ぶ。$K$ に台を持つ実 $\phi$ について $\|\phi\|_\infty\chi\pm\phi\ge0$ なので $|\langle T,\phi\rangle|\le\|\phi\|_\infty\langle T,\chi\rangle$ である。複素 $\phi$ は実部と虚部を別々に抑えればよい。非負の連続関数の近似には非負の核を用いられるから、延長も正であり、正の表現定理 :ref[thm-meas-rmk] により測度は正になる。
:::

:::example[主値超関数を測度で表せるか]
原点近くで $1$ の非負な偶切断 $\chi$ を選び、$\phi_\eps(x)=\chi(x)x/\sqrt{x^2+\eps^2}$ とする。上限ノルムは $\|\chi\|_\infty$ 以下だが
$$
\langle\operatorname{pv}(1/x),\phi_\eps\rangle
=\int_\R \frac{\chi(x)}{\sqrt{x^2+\eps^2}}dx
$$
である。$\chi=1$ の範囲 $|x|\le a$ で、$\eps\le x\le a$ の被積分関数は $1/(\sqrt2 x)$ 以上なので作用は発散する。共通台の零階評価が成り立たないため測度では表せない。元の主値の定義から一階の評価は成り立つので、階数は $1$ である。
:::

### 弱収束から共通の局所評価へ

各試験関数に対して値が有界である族は、固定したコンパクト台では一つの微分階数で抑えられる。このことにより、作用の点ごとの極限が再び超関数になることも保証できる。

:::proposition[点ごとに有界な超関数族]{#prop-four-distribution-uniform-bound}
$\mathcal F\subset\mathcal D'(\Omega)$ について各 $\phi$ で $\sup_{T\in\mathcal F}|\langle T,\phi\rangle|<\infty$ とする。各 $K\Subset\Omega$ についてある $m_K,C_K$ が存在し、
$$
\sup_{T\in\mathcal F}|\langle T,\phi\rangle|\le C_K p_{m_K}(\phi)
\quad(\phi\in\mathcal D_K)
$$
である。従って超関数列の作用が全試験関数で収束すれば、その極限は超関数を定める。
:::

:::proof
$\mathcal D_K$ は完備距離空間である。$j\ge1$ を整数として $F_j=\{\phi\in\mathcal D_K:\sup_{T\in\mathcal F}|\langle T,\phi\rangle|\le j\}$ と置く。各 $F_j$ は閉集合の交わりで閉であり、仮定より和集合は $\mathcal D_K$ 全体である。Baire の定理から、ある $F_j$ は内点 $\phi_0$ を持つ。$p_m(\phi)<\delta$ なら $\phi_0+\phi$ が $F_j$ に入るように選び、$\phi_0$ の値を引くと、全ての $T$ で絶対値は $2j$ 以下になる。線形性によって $\phi$ をこの近傍へ縮小すれば、例えば $C_K=4j/\delta$ で主張の評価を得る。$p_m(\phi)=0$ なら $p_0$ もゼロなので $\phi=0$ である。作用の極限にこの評価を移せば、線形性と :ref[prop-four-distribution-estimate] により極限も超関数である。
:::

## 超関数の微分と演算

通常の関数を微分して試験関数と積分すると、部分積分で微分を試験関数に移せる。この右辺だけを用いれば、関数としての導関数を持たない場合にも定義が残る。

:::definition[微分と滑らかな関数倍]{#def-four-distribution-derivative}
$T\in\mathcal D'(\Omega)$、$a\in C^\infty(\Omega)$ とし、

$$
\langle\partial^\alpha T,\phi\rangle=(-1)^{|\alpha|}\langle T,\partial^\alpha\phi\rangle,\qquad
\langle aT,\phi\rangle=\langle T,a\phi\rangle
$$

と定める。連続性は :ref[prop-four-distribution-estimate] と積の微分公式に従う。
:::

:::proposition[積の微分]{#prop-four-distribution-leibniz}
$\partial_j(aT)=(\partial_j a)T+a\partial_jT$ である。また、超関数の混合偏微分は交換する。
:::

:::proof
右辺の作用を展開すると $\langle T,(\partial_j a)\phi\rangle-\langle T,\partial_j(a\phi)\rangle=-\langle T,a\partial_j\phi\rangle$ であり左辺に一致する。混合偏微分は符号が同じで、試験関数の混合偏微分が交換する。
:::

:::example[デルタの導関数は必要な階数を持つ]
$\R$ 上で $a\in\R$、整数 $k\ge1$ とする。$\delta_a^{(k)}$ は定義から階数が高々 $k$ である。一方、$\phi_\eps(x)=\eps^k\psi((x-a)/\eps)$ とし、$\psi^{(k)}(0)=1$ の固定した試験関数を選ぶ。例えば原点近くで $1$ の切断 $\chi$ に対し $\psi(t)=t^k\chi(t)/k!$ でよい。$\eps\le1$ では台を固定でき、$p_{k-1}(\phi_\eps)=O(\eps)$ だが作用は $(-1)^k$ である。よって階数 $k-1$ の評価は存在せず、階数はちょうど $k$ である。
:::

:::example[全体で階数を固定できない超関数]
$T=\sum_{j\ge1}\delta_j^{(j)}$ とする。任意のコンパクト台には有限個の点 $j$ しか入らないので、定義する和は局所有限であり、:ref[prop-four-distribution-estimate] を満たす。一方、$j$ の近くの、他の整数を含まない区間では $T=\delta_j^{(j)}$ である。従って任意の候補 $m$ に対して $j>m$ を選ぶと、その区間に台を持つ試験関数だけで階数 $m$ の評価が破れる。
:::

:::example[跳躍の微分]{#exm-four-step-derivative}
$H=1_{(0,\infty)}$ とすると $\langle H',\phi\rangle=-\int_0^\infty\phi'=\phi(0)$ であり、$H'=\delta_0$ となる。一般に区分的に $C^1$ で有限個の跳躍を持ち、片側で $C^1$ が端点まで延長する関数では、区間ごとに部分積分すると、通常微分に $\sum_j(f(a_j+)-f(a_j-))\delta_{a_j}$ が加わる。
:::

:::example[デルタの導関数と座標倍]
$x\delta_0=0$ なので積の微分から $x\delta_0'=-\delta_0$ である。実際 $\langle x\delta_0',\phi\rangle=-(x\phi)'(0)=-\phi(0)$ となる。一方 $x\operatorname{pv}(1/x)=1$ は、積分の中で $x$ が約分され、除外した幅がゼロになることで従う。
:::

二つの任意の超関数の積は、この定義からは得られない。例えば $\rho$ を台の狭い滑らかな非負関数で $\int\rho=1$ とすると、$\rho_\eps(x)=\eps^{-1}\rho(x/\eps)$ は作用の意味で $\delta_0$ に収束するが、$\langle\rho_\eps^2,\phi\rangle=\eps^{-1}\int\rho^2(t)\phi(\eps t)dt$ は $\phi(0)>0$ のとき発散する。

### 滑らかな近似と弱収束

試験関数に対する値が収束することを :term[超関数の弱収束]{en="distributional convergence" yomi="ちょうかんすうのじゃくしゅうそく"} と呼ぶ。極限は点ごとの極限を指定する必要がない。むしろ、小さな領域に集中する列でも積分後の極限が残る。

:::example[差分でデルタの導関数を作る]
$T_\eps=(\delta_\eps-\delta_0)/\eps$（$\eps>0$）とする。作用の絶対値は $\sup_{[0,\eps]}|\phi'|$ 以下であり、極限は $\phi'(0)$ である。従って $T_\eps\to-\delta_0'$ となる。各項は階数 $0$ の符号付き測度だが、全変動は $2/\eps$ なので、二点を含む固定コンパクト近傍で階数 $0$ の共通評価は得られない。一階の評価なら $\eps$ に依存せず、この極限を扱える。
:::

:::proposition[超関数の正則化]{#prop-four-distribution-mollification}
$\rho\in\mathcal D(\R^n)$、$\int\rho=1$ とし $\rho_\eps=\eps^{-n}\rho(\cdot/\eps)$ と置く。$T\in\mathcal D'(\R^n)$ に対し
$$
u_\eps(x)=\langle T,\rho_\eps(x-\cdot)\rangle
$$
は滑らかで、$\partial^\alpha u_\eps(x)=\langle T,\partial^\alpha\rho_\eps(x-\cdot)\rangle$ である。さらに $u_\eps\to T$ が超関数の弱収束で成り立つ。
:::

:::proof
$x$ をコンパクト集合内で動かすと試験関数の台は共通のコンパクト集合に含まれる。$x$ に関する差商は全ての $y$ 導関数の一様収束を伴うので、:ref[prop-four-distribution-estimate] により差商を作用の外に出せる。これを繰り返せば滑らかさと導関数の式を得る。

試験関数 $\phi$ に対し、$\int \phi(x)u_\eps(x)dx=\langle T,\rho_\eps(-\cdot)*\phi\rangle$ である。交換は $x$ 積分を台のコンパクト性と、$x$ に関する連続性から全試験関数ノルムで Riemann 和に近似して正当化できる。右辺の試験関数は $\phi$ の各導関数の移動平均であり、移動の一様連続性と台の共通性から $\phi$ に収束する。
:::

:::example[点質量を滑らかにする]
$T=\delta_a$ なら $u_\eps(x)=\rho_\eps(x-a)$ であり、$T=\delta_a'$ なら $u_\eps(x)=\rho_\eps'(x-a)$ である。$a$ 以外の固定点でこれらはやがてゼロだが、試験関数に対する極限の作用はそれぞれ $\phi(a)$、$-\phi'(a)$ になる。
:::

:::example[振動して弱く収束する関数]
$f_j(x)=\cos(jx)$ とする。各 $\phi\in\mathcal D(\R)$ について部分積分すると $\int f_j\phi=-j^{-1}\int\sin(jx)\phi'(x)dx$ となり、絶対値は $\|\phi'\|_1/j$ 以下である。従って超関数としてゼロに収束する。さらに微分の定義により $\partial_x^k f_j\to0$ も弱収束として成り立つ。導関数の振幅が $j^k$ と増えることは、この弱収束に反しない。
:::

### 一般の開集合での滑らかな近似

任意の $T\in\mathcal D'(\Omega)$ は、$\Omega$ にコンパクト台を持つ滑らかな関数の列で弱近似できる。実際、コンパクトな汲み尽くし $K_j$ を取り、$K_j$ の近くで $1$ になる $\chi_j\in\mathcal D(\Omega)$ を選ぶ。$\chi_j T$ は、$\R^n$ の試験関数を $\chi_j$ 倍してから $T$ に渡すことで $\R^n$ の超関数になり、その台は $\operatorname{supp}\chi_j$ に含まれる。台が単位球に入る近似核 $\rho$ を選び、$\eps_j\le1/j$ かつ $\eps_j<\operatorname{dist}(\operatorname{supp}\chi_j,\R^n\setminus\Omega)$ として $u_j=(\chi_jT)*\rho_{\eps_j}$ を作る。$\Omega=\R^n$ では第二の条件は不要である。畳み込みの積分表示から $u_j$ の台は $\operatorname{supp}\chi_j$ の閉 $\eps_j$ 近傍に入るので $u_j\in\mathcal D(\Omega)$ である。

$\phi\in\mathcal D(\Omega)$ を固定する。十分大きな $j$ に対し、$\phi$ と $\rho_{\eps_j}(-\cdot)*\phi$ の台は、一つの $K_\ell$ に含まれる。$j\ge\ell$ では $\chi_j$ はこの台で $1$ であり、従って $\int u_j\phi=\langle T,\rho_{\eps_j}(-\cdot)*\phi\rangle$ となる。右辺の試験関数が共通台で全ての導関数について $\phi$ に収束するため、$u_j\to T$ が従う。

### 座標変換と曲面上のデルタ

通常関数 $f$ の変数変換 $f\circ F$ を積分で調べると、試験関数を逆写像で引き戻してから Jacobian を掛ける。従って元の関数が超関数であっても、滑らかな微分同相写像なら同じ式が使える。

:::definition[微分同相写像による引き戻し]{#def-four-distribution-pullback}
$F:U\to V$ を開集合の間の $C^\infty$ 級微分同相写像、$T\in\mathcal D'(V)$ とする。$F^*T\in\mathcal D'(U)$ を
$$
\langle F^*T,\phi\rangle
=\left\langle T,(\phi\circ F^{-1})\,|\det DF^{-1}|\right\rangle
$$
で定める。
:::

右辺の試験関数の台は $F(\operatorname{supp}\phi)$ に入る。台をコンパクト $K$ に固定すると、$F^{-1}$ と Jacobian の導関数は $F(K)$ 上で有界なので、連鎖律と積の微分から各階の上限は $K$ 上での有限階の上限で抑えられる。ゆえに右辺は連続であり、通常関数の引き戻しとの一致は積分の変数変換から従う。二つの微分同相写像を続けて使えば、逆写像の合成と行列式の連鎖律により $(F\circ G)^*=G^*F^*$ となる。

:::example[点質量を伸縮と移動で変える]
一次元で $F(x)=3x-1$、$T=\delta_2$ とする。$F^{-1}(2)=1$、$|(F^{-1})'(2)|=1/3$ であるから $F^*\delta_2=\delta_1/3$ である。$T=\delta_2'$ なら、試験関数 $\phi((y+1)/3)/3$ を $y=2$ で微分して符号を反転するため $\delta_1'/9$ になる。点の位置だけの移動なら係数はそのままだが、伸縮では各階で追加の因子が現れる。
:::

:::example[グラフ上の評価]
$g\in C^\infty(\R^{n-1})$ とする。$T\in\mathcal D'(\R^n)$ を $\langle T,\Phi\rangle=\int_{\R^{n-1}}\Phi(x',0)dx'$ で定める。$F(x',s)=(x',s-g(x'))$ は微分同相写像で逆写像の Jacobian は $1$ である。従って $F^*T$ は $\int\phi(x',g(x'))dx'$ と作用し、グラフに台を持つ。これを $\delta(s-g(x'))$ と書く場合、積分は $dx'$ に対するものである。通常の曲面積分では $\sqrt{1+|\nabla g(x')|^2}$ が掛かるので、両者を区別する。例えば $n=2$、$g(x)=2x$ のグラフでは、曲線の長さに対する測度を $1/\sqrt5$ 倍したものがこのデルタである。
:::

### 跳躍面をまたぐ微分

体積の関数の勾配は、領域の境界をまたぐとき表面上の測度を含み得る。符号を決めるには、外向き法線に向かってどちらの値へ変わるかを確認する。

:::proposition[滑らかな境界での微分]{#prop-four-distribution-interface}
$D\subset\R^n$ を有界な $C^1$ 級領域、$\nu$ を外向き単位法線とする。$f_-,f_+\in C^1(\R^n)$、$f=f_-1_D+f_+1_{\R^n\setminus D}$ と置くと、
$$
\partial_j f=(\partial_j f_-)1_D+(\partial_j f_+)1_{\R^n\setminus D}
 +(f_+-f_-)\nu_j\,dS\big|_{\partial D}.
$$
右端は $\phi$ に対して $\int_{\partial D}(f_+-f_-)\nu_j\phi\,dS$ と作用する。
:::

:::proof
$f=f_++(f_--f_+)1_D$ と書く。$\langle\partial_j 1_D,\phi\rangle=-\int_D\partial_j\phi=-\int_{\partial D}\nu_j\phi\,dS$ が発散定理から成り立つ。$f_--f_+$ はここでは $C^1$ なので、同じ発散定理を $(f_--f_+)\phi$ に使って積の微分式を直接得られる。領域内部の項を整理すると通常の各偏微分になり、境界項の係数は $(f_+-f_-)\nu_j$ である。
:::

:::example[直線を境に変わる関数]
$\R^2$ の $u(x,y)=1_{\{y>2x\}}$ では $y$ 積分を半直線で部分積分して $\partial_yu=\delta(y-2x)$ である。$x<(y/2)$ の範囲を先に積分すると
$$
\langle\partial_xu,\phi\rangle=-\int_\R \phi(y/2,y)dy
=-2\int_\R \phi(x,2x)dx
$$
となり $\partial_xu=-2\delta(y-2x)$ である。勾配は $( -2,1)$ に向き、領域の内向き単位法線は $(-2,1)/\sqrt5$ である。グラフのデルタと表面積の因子を合わせると、$1_D$ の勾配が外向き法線の負になる符号も確認できる。
:::

### 一次元で任意の超関数の原始関数を作る

原始関数を求めるとき、試験関数を微分して得られる関数の積分は常にゼロである。その一つの制約を固定した試験関数で調整すると、任意の超関数から原始関数を構成できる。

:::proposition[区間上の原始超関数]{#prop-four-distribution-primitive}
開区間 $I$ 上の任意の $T\in\mathcal D'(I)$ に対し、$U'=T$ を満たす $U\in\mathcal D'(I)$ が存在する。その二つの解の差は定数である。
:::

:::proof
$\eta\in\mathcal D(I)$、$\int\eta=1$ を固定する。$\phi$ に対して
$$
A\phi(x)=\int_{-\infty}^x\left(\phi(t)-\eta(t)\int_I\phi\right)dt
$$
と置く。積分の被積分関数は全積分がゼロなので、$A\phi$ は $\phi,\eta$ の台の凸包に台を持つ。$\phi$ の台が固定のコンパクトで動くと、この凸包も一つのコンパクトな区間に入る。$A\phi$ の零階の上限は被積分関数の上限と区間長の積で抑えられ、正の階数では被積分関数の一つ下の階数の上限で抑えられる。よって $\langle U,\phi\rangle=-\langle T,A\phi\rangle$ は超関数を定める。$\int_I\phi'=0$ なので $A(\phi')=\phi$ となり、$\langle U',\phi\rangle=\langle T,\phi\rangle$ である。

逆に $V'=0$ なら積分ゼロの任意の試験関数 $\psi$ は $(A\psi)'$ だから作用がゼロである。任意の $\phi$ から $\eta\int\phi$ を引いたものの積分はゼロであり、従って $\langle V,\phi\rangle=\langle V,\eta\rangle\int\phi$ である。これは定数を表す。
:::

:::example[原始関数の定数は連続な形で残る]
$T=\delta_2-\delta_5$ を $\R$ 上で考える。$U=1_{(2,5)}$ とすれば $-\int_2^5\phi'=\phi(2)-\phi(5)$ なので解である。:ref[prop-four-distribution-primitive] より一般解は $1_{(2,5)}+c$ である。$U$ が $(-\infty,2)$ 上でゼロである条件を追加すれば $c=0$ になる。この条件は端点の関数値を指定するものではなく、その開集合に台を持つ全試験関数で作用がゼロであることを意味する。
:::

## 超関数の台とコンパクト台の超関数

:::definition[超関数の台]{#def-four-distribution-support}
開集合 $U$ 上で $T$ がゼロであるとは、$U$ にコンパクト台を持つ全試験関数に作用がゼロであることをいう。$T$ がゼロである最大の開集合の補集合を :term[台]{en="support" yomi="だい"} と呼ぶ。最大の開集合は、試験関数の台を有限個の局所切断で分割できるため存在する。実際、有限部分被覆の中で非負の切断関数の和を台の近傍で正にし、その和で割って試験関数を分解する。
:::

台の補集合で $T$ に試験関数を作用させるとゼロであるから、微分や滑らかな関数倍で台は広がらない。例えば $\operatorname{supp}H=[0,\infty)$ だが $\operatorname{supp}H'=\{0\}$ である。

:::proposition[コンパクト台の超関数の延長]{#prop-four-compact-distribution}
$T$ の台がコンパクトなら $C^\infty(\Omega)$ に連続に一意に延長できる。台の近傍で $1$ となる切断 $\chi$ を用い $\langle T,\phi\rangle=\langle T,\chi\phi\rangle$ とする。逆に $C^\infty(\Omega)$ の通常のコンパクト上一様微分の位相で連続な線形汎関数はコンパクト台を持つ。
:::

:::proof
異なる切断の差は台の近傍でゼロなので値は等しい。$p_m(\chi\phi)\le C\max_{|\alpha|\le m}\sup_{\operatorname{supp}\chi}|\partial^\alpha\phi|$ から連続性が従う。$C^\infty$ の任意の元に、汲み尽くすコンパクト集合の近くで $1$ となる切断を掛ければ、$\mathcal D$ の列が元へ収束する。ゆえに延長は一意である。逆に連続性を有限個のセミノルムで評価し、そのコンパクト集合の和集合を取れば、そこから離れた試験関数では作用がゼロである。
:::

:::theorem[一点に台を持つ超関数]{#thm-four-point-supported}
$\operatorname{supp}T\subset\{a\}$ なら $T$ は有限個の $\partial^\alpha\delta_a$ の線形結合である。
:::

:::proof
原点に移し、固定した小閉球上の局所評価の次数を $m$ とする。$\partial^\alpha\phi(0)=0$ が $|\alpha|\le m$ で成り立つと仮定する。原点近くで $1$ となる $\chi(x/\eps)$ を掛けても作用は変わらない。Taylor の定理で $\partial^\beta\phi(x)=o(|x|^{m-|\beta|})$（$|\beta|\le m$）である。積の微分で $\partial^\gamma\chi(x/\eps)=O(\eps^{-|\gamma|})$ を掛けると、次数 $k\le m$ の各導関数の上限は $o(\eps^{m-k})$ に収まる。従って作用はゼロである。一般の $\phi$ から原点の Taylor 多項式に固定切断を掛けたものを引けば、作用が有限個の $\partial^\alpha\phi(0)$ の線形結合だと分かる。定義によりそれらは $(-1)^{|\alpha|}\partial^\alpha\delta_0$ で表せる。
:::

:::example[二点に台がある場合]
$\operatorname{supp}T\subset\{-1,2\}$ とする。二点の互いに素な近傍にそれぞれ台を持ち、その点の近くで $1$ となる $\chi_1,\chi_2$ を選ぶと $T=\chi_1T+\chi_2T$ である。差の乗数が台の近傍でゼロだからである。各項は一点に台を持つので :ref[thm-four-point-supported] で有限ジェットの和になる。例えば $\phi$ が $-1$ 近傍でゼロなら、残るのは $2$ での有限階の導関数だけである。
:::

:::example[コンパクト台の作用を多項式へ延ばす]
$T=\delta_0''+2\delta_1$ とし、$\phi(x)=x^2+3$ を考える。:ref[prop-four-compact-distribution] による作用は $\phi''(0)+2\phi(1)=2+8=10$ である。切断を用いても両点の近傍で $1$ なら、それぞれで必要な値と導関数は変わらない。
:::

### 制限と貼り合わせ

超関数を開部分集合に制限するとき、その部分集合の試験関数は境界近くで既にゼロである。従って外側をゼロにする延長は滑らかであり、その延長に元の超関数を作用させれば制限を定められる。一方、複数の領域で別々に得た作用をつなぐ際には、重なりでの一致を確認する。

:::proposition[制限された超関数の貼り合わせ]{#prop-four-distribution-gluing}
$\Omega=\bigcup_{i\in I}U_i$ を開被覆とし、$T_i\in\mathcal D'(U_i)$ が全ての重なり $U_i\cap U_j$ で一致するとする。全ての $i$ について $T|_{U_i}=T_i$ であるような $T\in\mathcal D'(\Omega)$ が一意に存在する。
:::

:::proof
試験関数 $\phi$ の台 $K$ を覆う有限個の $U_i$ を選び、そこで作った切断 $\theta_i$ により $\phi=\sum_i\theta_i\phi$ とする。$\langle T,\phi\rangle=\sum_i\langle T_i,\theta_i\phi\rangle$ と置く。別の分割 $\eta_j$ を使ったときは、各 $\theta_i\phi$ に $\sum_j\eta_j=1$ を掛けると、比較する二つの和は $\theta_i\eta_j\phi$ 上の作用に分解される。その試験関数は $U_i\cap U_j$ に台を持つため二つの作用は一致する。異なる $\phi$ に対する線形性も、二つの台の合併を同じ有限分割で覆って確かめられる。

$K$ を固定すれば分割は固定でき、各 $\theta_i\phi$ の台は $\operatorname{supp}\theta_i$ と $K$ の交わりに入る。有限個の $T_i$ の局所評価を取り、その中の最大の階数と各切断の導関数の上界を用いると $T$ の局所評価が成り立つ。$\phi$ の台が $U_i$ 内にある場合は $\theta_i$ 一つで分割できるので制限は $T_i$ である。全制限が一致する二つの作用の差は同じ分割で全ての試験関数を消すため、一意性も従う。
:::

:::example[開被覆の重なりだけで定数を決める]
$U_1=(-3,1)$、$U_2=(-1,3)$ とする。$T_1$ を $U_1$ 上の $1_{(0,\infty)}$、$T_2$ を $U_2$ 上の $\delta_0$ とすると、両方とも原点に関係するが重なりで一致しない。実際、$(1/4,1/2)$ に台を持つ積分 $1$ の試験関数に対し、第一の値は $1$、第二の値は $0$ である。$T_2$ を $1_{(0,\infty)}+2\delta_2$ に替えると、重なりでは $\delta_2$ が消えるので貼り合わせられる。得られる超関数は $(-3,3)$ 上の $1_{(0,\infty)}+2\delta_2$ であり、台は相対閉集合 $[0,3)$ である。
:::

### 滑らかさを失う点を記録する

超関数 $T$ が滑らかである最大の開集合の補集合を :term[特異台]{en="singular support" yomi="とくいだい"} という。ここで開集合上で滑らかとは、ある $C^\infty$ 関数がその制限を表すことをいう。二つの滑らかな関数が重なり上で同じ超関数を表せば a.e. に一致し、両者の連続性から各点で一致する。このため局所代表元を貼り合わせて最大の開集合が得られる。

例えば $H$ の台は $[0,\infty)$ で、特異台は $\{0\}$ である。正の半直線では $1$ を、負の半直線では $0$ を代表元にできる。しかし原点を含む開区間では、一つの連続関数が両側の定数を表すことはできない。一方、滑らかな非零関数の台は空でなくても特異台は空である。微分で特異台は広がらない。元の超関数が滑らかな開集合では、微分を行うとその滑らかな代表元の導関数が新しい代表元となるためである。

### 台と零点を用いる方程式

$aT=0$（$a$ は滑らかな関数）なら、$a$ が零でない開集合で、任意の試験関数を $a$ で割って作用させられる。そのため $\operatorname{supp}T\subset\{a=0\}$ である。この包含は必要条件であるが、$T$ がデルタの導関数を含む場合には $a$ の微分も作用に影響する。

:::proposition[有限個の零点を持つ乗数の核]{#prop-four-distribution-multiplier-kernel}
一次元で $a\in C^\infty(I)$ が有限個の零点 $a_1,\ldots,a_N$ を持ち、それぞれ有限な位数 $m_1,\ldots,m_N$ で消えるとする。$aT=0$ の全解は
$$
T=\sum_{j=1}^N\sum_{r=0}^{m_j-1} c_{j,r}\delta_{a_j}^{(r)}
$$
である。
:::

:::proof
上の包含から $T$ の台は零点集合に入り、互いに素な近傍で切断して一点ずつに分けられる。各点の近くで Taylor の積分剰余により $a(x)=(x-a_j)^{m_j}b_j(x)$、$b_j(a_j)\ne0$ と書ける。$b_j$ はその近傍で滑らかに可逆であるから、その近傍に台を持つ試験関数によって $aT=0$ を調べることは $(x-a_j)^{m_j}T=0$ を調べることと同値である。

:ref[thm-four-point-supported] によって一点で $T$ を $\delta_{a_j}^{(r)}$ の有限和にできる。Leibniz の公式より $(x-a_j)\delta_{a_j}^{(r)}=-r\delta_{a_j}^{(r-1)}$（$r\ge1$）であり、$r=0$ ならゼロである。従って $m_j$ 回掛けると $r<m_j$ の項は消え、$r\ge m_j$ の項は非零定数倍の $\delta_{a_j}^{(r-m_j)}$ になる。これらの次数は異なり、切断付きの単項式で各ジェットを独立に指定できるので係数は全てゼロである。残りの $r<m_j$ の組合せは実際に $aT=0$ を満たす。
:::

:::example[零点の位置と重複度を分ける]
$a(x)=(x-1)^2(x+2)$ なら全解は $c_0\delta_1+c_1\delta_1'+d_0\delta_{-2}$ である。$\delta_1''$ を加えると、$a(1)=a'(1)=0$ であっても $a''(1)=6$ なので $a\delta_1''=6\delta_1$ が残る。$-2$ では一階の導関数が $9$ なので、$a\delta_{-2}'=-9\delta_{-2}$ である。点で $a$ がゼロになるだけではデルタの微分の係数を自由に選べないことが分かる。
:::

:::example[台は離散点の集合そのものとは限らない]
$T=\sum_{j\ge1}2^{-j}\delta_{1/j}$ は $\R$ 上の超関数であり、作用は $\|\phi\|_\infty$ 以下である。台は $\{1/j:j\ge1\}\cup\{0\}$ である。実際、点列とその極限を避ける開集合では作用がゼロであり、どの $1/j$ の小さい近傍にも正の作用がある試験関数を取れる。さらに $0$ のどの近傍にもある $1/j$ が入り、同様にゼロとはならない。$0$ 自身の点質量はゼロだが、台の定義は各近傍での作用を調べるため $0$ も台に入る。
:::

## 超関数の畳み込み

二つの超関数を足し算 $x+y$ に沿ってまとめるには、$\phi(x+y)$ を二変数の試験関数として扱う必要がある。この関数の台は一般には非有界である。そこで、二つの台の上で $x+y$ が有界なら $x,y$ も有界になる条件を課す。

:::proposition[テンソル積と畳み込み]{#prop-four-distribution-convolution}
$S,T\in\mathcal D'(\R^n)$ に対し、$\Phi\in\mathcal D(\R^{2n})$ への作用

$$
\langle S\otimes T,\Phi\rangle=\langle S_x,\langle T_y,\Phi(x,y)\rangle\rangle
$$

は超関数を定め、作用の順序を交換できる。さらに加法写像が $\operatorname{supp}S\times\operatorname{supp}T$ に制限して固有である場合、つまりコンパクト集合の逆像がコンパクトである場合、

$$
\langle S*T,\phi\rangle=\langle S\otimes T,\chi(x,y)\phi(x+y)\rangle
$$

で畳み込みを定められる。$\chi$ は台の積と $\{x+y\in\operatorname{supp}\phi\}$ の交わりの近傍で $1$ の切断である。一方の台がコンパクトなら条件を満たす。
:::

:::proof
$\Phi$ の台の射影はコンパクトである。$x$ による差商の収束は、$y$ に関する任意次数の一様微分でも成り立つため、内側の作用は滑らかで $x$ 射影に台を持つ。各作用を局所評価すると、全体は $\Phi$ の有限次数の導関数の上限で抑えられる。有限和の積 $\Phi(x,y)=\sum_j a_j(x)b_j(y)$ では順序交換が成り立つ。一般の $\Phi$ を大きな直方体に収めて周期延長し、部分積分で急減少する Fourier 係数の級数を取ると、一様に全階の導関数まで収束する。各座標群に固定切断を掛ければ、積の有限和で共通台の近似を得る。これにより交換を一般に移せる。

$S\otimes T$ の台は $\operatorname{supp}S\times\operatorname{supp}T$ に含まれる。積の台から離れた点は、少なくとも一方の因子で作用がゼロとなる開直方体を持つ。試験関数を有限個のそのような直方体に分ければ作用はゼロである。従って上の切断は選択によらない。$\phi$ の台を固定すれば同じ $\chi$ を使え、その微分が局所評価を与える。一方の台が有界なら $x+y$ の有界性がもう一方の変数の有界性を与える。逆像は閉でもあるので固有性が従う。
:::

定義から $S*T=T*S$ である。$\partial^\alpha(S*T)=(\partial^\alpha S)*T=S*(\partial^\alpha T)$ は、試験関数の $x+y$ に関する微分が $x$ 微分にも $y$ 微分にも等しいことから従う。切断の導関数が現れる項は、台の積と対象の逆像の近傍でゼロなので消える。同様に三つの台の積で $x+y+z$ の逆像がコンパクトで、各段階の畳み込みが定義できれば、三変数の作用に帰着して結合できる。

:::example[点による平行移動と半直線上の畳み込み]
$\langle\delta_a*T,\phi\rangle=\langle T,\phi(a+\cdot)\rangle$ であり、これは $T$ を $a$ だけ右に移す操作である。$H*H$ では両変数が非負で和の上限が変数の上限になる。通常の積分で $(H*H)(x)=0$（$x<0$）、$x$（$x>0$）である。対して定数 $1$ 同士では $x+y=0$ を満たす組が非有界で、通常の畳み込み積分も発散する。
:::

:::example[台をずらしてから微分する]
$T=1_{[0,1]}$ として $\delta_2'*T$ を求める。まず $\delta_2*T=1_{[2,3]}$、その微分が $\delta_2-\delta_3$ である。畳み込みと微分の交換が成り立つので、先に $\delta_2$ を微分しても同じ点質量の差になる。
:::

### 台の制御と結合則

和の逆像がコンパクトになる条件を使うと、畳み込みの台を個々の台の和で制御できる。ただし、台の和に入るという情報だけでは、その中の各点で作用が残るとは限らない。

:::proposition[台の和と半直線上の結合則]{#prop-four-convolution-support-associative}
$A=\operatorname{supp}S$、$B=\operatorname{supp}T$ 上で加法写像が固有なら、$A+B$ は閉であり $\operatorname{supp}(S*T)\subset A+B$ である。一次元で三つの超関数の台が全て下に有界なら、両方の順序の三重畳み込みが存在し、結合則が成り立つ。
:::

:::proof
$z_j=a_j+b_j$ が $z$ へ収束するなら、集合 $\{z\}\cup\{z_j:j\ge1\}$ はコンパクトである。固有性により $(a_j,b_j)$ はコンパクト集合内に入るので適当な部分列の極限 $(a,b)$ を持ち、$a\in A,b\in B,a+b=z$ となる。このため $A+B$ は閉である。台を $A+B$ から外れた場所に持つ試験関数を取ると、$\phi(x+y)$ は $A\times B$ の近くでゼロなので作用はゼロである。

三つの台がそれぞれ $[a,\infty),[b,\infty),[c,\infty)$ に入るとする。二変数でも三変数でも合計が上から抑えられれば、各変数は他の下限を引いた数で上から抑えられる。よって和写像の逆像は閉かつ有界で、必要な各段の畳み込みが定義される。固定した試験関数 $\phi$ の台に対して、三変数の逆像全体とその中間和の射影を含むコンパクト集合を取る。各変数の切断と中間和の切断を、この集合の各射影の近傍で $1$ となるように選ぶ。例えば $(S*T)*U$ の作用を展開すると、三重のテンソル積が $x+y+z$ のところで $\phi$ に作用する式に、これらの切断が掛かる。$S*(T*U)$ の方でも同じ式になり、切断の違いは三つの台の積と逆像の近くでゼロである。よって作用は一致する。
:::

:::example[半直線上の繰り返し積分]
$H=1_{(0,\infty)}$ の $k$ 重畳み込みは $x>0$ で $x^{k-1}/(k-1)!$、$x<0$ でゼロである。$k=1$ で定義と一致し、$k$ 重の式を $H$ と畳み込むと $x>0$ で $\int_0^x y^{k-1}/(k-1)!dy=x^k/k!$ になる。$x<0$ では積分範囲が空である。例えば三重畳み込みの連続代表元は $x=2$ で $2$ になる。結合則によって、どの順番で二つずつ計算しても同じ結果が得られる。
:::

:::example[台の和で可能だった点が相殺する]
$S=\delta_0-\delta_1$、$T=\delta_0+\delta_1$ とする。台の和は $\{0,1,2\}$ であるが、畳み込みは $\delta_0+\delta_1-\delta_1-\delta_2=\delta_0-\delta_2$ となる。$x=1$ へ至る二組の寄与が逆符号だからである。非負測度同士の場合は、任意の $a\in A,b\in B$ の近傍に正の質量があるので、$a+b$ の十分小さい任意の近傍にも積から正の質量が入る。その場合は固有性の下で台が和に一致する。
:::

### 局所的に連続関数の微分で表す

有限階の導関数だけで評価できる作用に対し、十分な回数の積分を行えば連続関数へ移せる。超関数の微分はこれを逆にたどる操作なので、連続関数から任意の局所的な超関数を表す方法になる。

:::proposition[連続関数による局所表示]{#prop-four-distribution-local-structure}
$T\in\mathcal D'(\Omega)$、$U\Subset\Omega$ を開集合とする。ある整数 $r\ge2$ と連続関数 $V$ が存在して、$U$ 上で
$$
T=\partial_1^r\cdots\partial_n^r V
$$
となる。
:::

:::proof
$\overline U$ 近くで $1$ の切断を掛け、$\Omega$ の外にゼロで延ばしたコンパクト台超関数を $S$ とする。固定切断で試験関数を台の近くに局所化できるため、$S$ の作用は固定コンパクト集合で $m$ 階までの上限によって評価される。よって $C^m$ 級関数にも作用を延長できる。実際、切断付きの $C^m$ 関数を近似単位で平滑化すると各導関数が $m$ 階まで一様に収束し、共通台に入る。その作用の極限は評価により一意である。

$r=m+2$ と置き、$L_r(t)=(\max(t,0))^{r-1}/(r-1)!$、$K(x)=\prod_{j=1}^n L_r(x_j)$ とする。$L_r$ は $C^m$ 級であり、$r-1$ 回微分すると $H$、さらに微分すると $\delta_0$ になる。従って積 $K$ では $\partial_1^r\cdots\partial_n^r K=\delta_0$ である。これは各変数で試験関数との積分を行い、定義の符号を含めて順に導関数を移せば分かる。

$V(x)=\langle S_y,K(x-y)\rangle$ と置く。$x$ がコンパクト集合内にあるとき、$K$ の $m$ 階までの全導関数は $x-y$ の走るコンパクト集合上で一様連続である。従って $x\mapsto K(x-\cdot)$ は作用を評価する $C^m$ の上限ノルムで連続であり、$V$ も連続である。$V$ は $S*K$ を表す。実際、$x$ を試験関数の台に制限し、$y$ に固定切断を掛ければ、$\int K(x-y)\phi(x)dx$ は $y$ の $m$ 階までの上限ノルムで Riemann 和から得られる。これに $S$ を作用させれば、畳み込みの定義と一致する。微分と畳み込みの交換によって $\partial_1^r\cdots\partial_n^r V=S*\delta_0=S$ を得る。$U$ では $S=T$ なので結論が従う。
:::

:::example[階数と積分回数は異なる数である]
$T=\delta_a'$ は階数 $1$ である。上の構成では $r=3$、$K(t)=t_+^2/2$ とでき、その作用は $V(x)=K'(x-a)=(x-a)_+$ である。実際 $V'=H(x-a)$、$V''=\delta_a$、$V'''=\delta_a'$ となる。$T$ の作用を上限ノルムで評価するために必要な導関数は一階だが、連続関数の三階微分で表している。これらの階数を同じと考える必要はない。
:::

## 基本解と Malgrange–Ehrenpreis の定理（概観）

:::definition[基本解]{#def-four-fundamental-solution}
線形定数係数微分作用素 $P(\partial)$ に対し、$P(\partial)E=\delta_0$ を満たす超関数 $E$ を :term[基本解]{en="fundamental solution" yomi="きほんかい"} という。コンパクト台のデータ $F$ に対し $u=E*F$ とすると、微分の交換で $P(\partial)u=F$ となる。
:::

:::example[一次元と三次元の基本解]
一次元で $E(x)=|x|/2$ とすれば $E'=\sgn x/2$、$E''=\delta_0$ である。三次元の $-\Delta$ に対しては $E(x)=1/(4\pi|x|)$ である。原点外で $\Delta E=0$ であり、半径 $\eps$ の球を除いて Green の公式を使うと

$$
\int_{|x|>\eps} E(-\Delta\phi)\,dx
=\int_{|x|=\eps}(-E\,\partial_\nu\phi+\phi\,\partial_\nu E)\,dS .
$$

内側境界の外向き法線 $\nu$ は $-x/|x|$ である。第 1 項は $O(\eps)$、第 2 項は球面上の $\phi$ の平均で $\phi(0)$ に収束する。$\eps$ 球内での $E\Delta\phi$ の積分は $O(\eps^2)$ なので結論を得る。
:::

:::theorem[Malgrange–Ehrenpreis]{#thm-four-malgrange-ehrenpreis}
非零の定数係数多項式 $P$ で定まる $P(\partial)$ は $\R^n$ 上の超関数の基本解を持つ。
:::

:::proof
定数作用素では $E=P^{-1}\delta_0$ を取れる。次数 $m>0$ とし、最高次斉次部分が消えない実方向 $v$ を選ぶ。$p(\xi)=P(2\pi i\xi)$ とすると、$\zeta\mapsto p(\xi+\zeta v)$ は次数 $m$ で最高次係数 $c\ne0$ が $\xi$ によらない。半径 $1,2,\ldots,m+1$ のうち、全ての根の絶対値と距離 $1/3$ 以上離れたものがある。各根は高々一つの半径を除外するからである。その円上で $|p|\ge |c|3^{-m}$ となる。円上の最小値がこの下界以上になる最初の半径を $r(\xi)$ として選ぶと、各円の最小値は連続なので $r$ は可測である。

$\check\phi(z)=\int\phi(x)e^{2\pi i x\cdot z}dx$ と置き、

$$
\langle E,\phi\rangle=\int_{\R^n}\frac1{2\pi}\int_0^{2\pi}
\frac{\check\phi(\xi+r(\xi)e^{it}v)}
{p(\xi+r(\xi)e^{it}v)}\,dt\,d\xi
$$

で定める。固定台 $K$ では指数の虚部は有界な範囲を動き、$e^{-2\pi x\cdot\Im z}$ とその導関数は $x\in K$ 上で一様有界である。$x$ で部分積分を繰り返し、$1-\Delta_x$ を指数の実部に対して使うと $|\check\phi(z)|\le C_Kp_{2l}(\phi)(1+|\Re z|^2)^{-l}$ となる。$2l>n$ とすれば $\Re z-\xi$ が有界なので積分は絶対収束し、連続性の評価を与える。

$P(-\partial)\phi$ の逆変換は部分積分により $p(z)\check\phi(z)$ である。従って $\langle P(\partial)E,\phi\rangle$ の分母は消え、$\zeta\mapsto\check\phi(\xi+\zeta v)$ の円平均は整関数の平均値公式で $\check\phi(\xi)$ になる（:ref[thm-cpx-cauchy-formula]）。反転公式で $\int\check\phi(\xi)d\xi=\phi(0)$ であるから結論を得る。
:::

### 基本解の選択と初期値

:::proposition[基本解の差]{#prop-four-fundamental-difference}
一つの基本解 $E_0$ が分かっているとき、全ての基本解は $E_0+V$（$P(\partial)V=0$）の形である。
:::

:::proof
二つの基本解の差を取ると右辺のデルタが消える。逆に斉次方程式を満たす $V$ を足すと、線形性で右辺は変わらない。
:::

:::example[前向きと後ろ向きの原始関数]
$\partial$ の基本解として $H$ と $-H(-x)$ の両方が使える。後者の微分でも、反射の符号と前置符号が打ち消され $\delta_0$ になる。差は定数 $1$ である。コンパクト台の可積分関数 $F$ と畳み込むと、それぞれ $u(x)=\int_{-\infty}^x F(s)ds$ と $u(x)=-\int_x^\infty F(s)ds$ になる。一般のコンパクト台超関数 $F$ に対しては畳み込みを試験関数への作用で定める。二つの解の差は定数 $\langle F,1\rangle$ であり、:ref[prop-four-compact-distribution] によりこの作用も定義できる。
:::

### 一次元の前向き基本解を組み立てる

多項式を一次式に分解すると、定数係数常微分作用素に対して台が正の半直線にある基本解を作れる。この台は、右辺を与えた点より前に解へ寄与させない条件である。

:::proposition[定数係数常微分作用素の前向き基本解]{#prop-four-causal-ode}
$P(z)=c\prod_{j=1}^m(z-\lambda_j)$、$c\ne0$ とする。$m\ge1$ なら
$$
E=c^{-1}(e^{\lambda_1 x}H)*\cdots*(e^{\lambda_m x}H)
$$
が $P(\partial)$ の基本解であり、$\operatorname{supp} E\subset[0,\infty)$ を満たす基本解はこれに限る。$m=0$ の場合は $E=c^{-1}\delta_0$ である。
:::

:::proof
各 $E_j=e^{\lambda_j x}H$ は $E_j'=\lambda_j E_j+\delta_0$ を満たす。実際積の微分のデルタ項では滑らかな乗数を原点で評価し、値は $1$ である。全ての台が $[0,\infty)$ に含まれるので :ref[prop-four-convolution-support-associative] により多重畳み込みが定義でき、結合順序によらない。$(\partial-\lambda_j)$ を $j$ 番目の因子へ移すと $\delta_0$ になり、畳み込みの単位元であるから除ける。式の全因子にこれを行うと $P(\partial)E=\delta_0$ である。

一意性では、同じ台を持つ二つの解の差 $V$ に $P(\partial)V=0$ を適用する。$W=\prod_{j=2}^m(\partial-\lambda_j)V$ と置けば $(\partial-\lambda_1)W=0$ である。微分と乗算で台は広がらないから $W$ の台も $[0,\infty)$ に含まれる。一方 $(e^{-\lambda_1x}W)'=0$ なので :ref[prop-four-distribution-primitive] より $e^{-\lambda_1x}W$ は定数である。その定数は負の半直線でゼロなので全体でゼロである。残る因子について同じ議論を繰り返して $V=0$ を得る。
:::

:::example[同じ根の重複と異なる二根]
二因子が同じ $\lambda$ なら、正の $x$ で $\int_0^x e^{\lambda(x-y)}e^{\lambda y}dy=xe^{\lambda x}$ となる。一般に同じ根を $m$ 回取れば $e^{\lambda x}x^{m-1}H(x)/(m-1)!$ である。一方 $\lambda\ne\mu$ なら二因子の畳み込みは $(e^{\lambda x}-e^{\mu x})H(x)/(\lambda-\mu)$ である。例えば $P(\partial)=\partial^2+4$ では $\lambda=2i,\mu=-2i$ なので $E=(\sin(2x)/2)H(x)$ になる。原点で連続かつ右導関数が $1$、左導関数が $0$ なので、$(\partial^2+4)E$ には原点のデルタが一つだけ残る。
:::

### 低次元の Laplace 作用素

三次元の計算は $n\ge3$ に拡張できる。単位球面の面積を $\sigma_{n-1}$ とし、$E(x)=|x|^{2-n}/((n-2)\sigma_{n-1})$ と置く。放射状関数の Laplace 作用素は $\Delta u=r^{1-n}(r^{n-1}u'(r))'$ なので、$E$ は原点を除いて調和である。半径 $\eps$ の内側境界では $\partial_\nu E=\eps^{1-n}/\sigma_{n-1}$ である。そのため三次元の境界積分と同じく第二項は球面平均になり、第一項は球面積と $E$ の積が $O(\eps)$ になるので消える。除いた球内での $|E\Delta\phi|$ の積分も $\int_0^\eps r\,dr$ の定数倍以下である。従ってこの $E$ について $-\Delta E=\delta_0$ となる。

$n=2$ では放射状の調和方程式が $(ru'(r))'=0$ になり、非定数の解が対数になる。そこで
$$
E(x)=-\frac1{2\pi}\log|x|\qquad (x\in\R^2\setminus\{0\})
$$
と置く。原点では局所可積分である。実際 $\int_0^\eps r|\log r|dr=O(\eps^2(1+|\log\eps|))$ である。$|x|>\eps$ では調和であり、内側境界の外向き法線に対して $\partial_\nu E=1/(2\pi\eps)$ となる。従って Green の公式で得られる境界項は、円周上の $\phi$ の平均と、絶対値が $O(\eps|\log\eps|)$ の項の和である。前者は $\phi(0)$ に収束し、後者と除いた円板上の誤差はゼロに収束する。よって $-\Delta E=\delta_0$ である。$n=1$ では $-\partial_x^2$ の基本解は $-|x|/2$ である。一階微分は $-\sgn(x)/2$、その微分は両側の跳躍 $-1$ により $-\delta_0$ だからである。

:::example[二つの点源を重ねる]
$\R^2$ で $a=(1,0)$、$b=(-1,0)$ とする。右辺 $\delta_a-\delta_b$ に対し、$u(x)=(2\pi)^{-1}\log(|x-b|/|x-a|)$ は基本解の平行移動の差なので $-\Delta u=\delta_a-\delta_b$ を満たす。原点では二つの距離が等しいから $u(0)=0$ であり、$x=(2,0)$ では距離の比が $3$ なので $u(2,0)=\log3/(2\pi)$ である。各点源の基本解の対数増大は差で打ち消され、$|x|\to\infty$ では $|x-b|/|x-a|\to1$ なので $u(x)\to0$ となる。
:::

### 基本解を右辺に作用させる

非零の定数係数作用素 $P(\partial)$ の基本解を $E$ とし、コンパクト台超関数 $F$ に対して $u=E*F$ とする。一方の台がコンパクトなので畳み込みが定義でき、$P(\partial)u=(P(\partial)E)*F=\delta_0*F=F$ である。$F$ がコンパクト台の滑らかな関数なら、$E$ の正則化と同じ差商の議論で $u(x)=\langle E_y,F(x-y)\rangle$ は滑らかである。従ってこの場合は等式を古典微分でも読める。右辺がコンパクト台でない場合は、この畳み込みをそのまま定義できるかどうかを台の固有性で改めて判定する。例えば $F$ と $E$ の台がともに時間の下側で有界なら、一次元の前向き解にも同じ畳み込みを使える。

:::example[遅れて入る二つの衝撃]
$\R$ 上で $P(\partial)=\partial^2+4$、$F=\delta_1+2\delta_3$ とする。前向き基本解の畳み込みは
$$
u(x)=\tfrac12 H(x-1)\sin(2(x-1))+H(x-3)\sin(2(x-3)).
$$
$x<1$ でゼロであり、$1<x<3$ では第一項だけが動く。$x=1,3$ で $u$ は連続だが、導関数にそれぞれ $1,2$ の跳躍がある。二階微分のデルタ項が右辺を与え、各開区間で $u''+4u=0$ である。先に斉次解を加えれば同じ右辺を満たしても $x<1$ での条件が変わるので、この条件によって解を選んでいる。
:::

## 演習問題

::::exercise{level=basic}
$1_{[1,3]}'$ を超関数として求めよ。
:::solution
$-\int_1^3\phi'(x)dx=\phi(1)-\phi(3)$ なので $\delta_1-\delta_3$ である。
:::
::::

::::exercise{level=basic}
$x^2\delta_0''$ を求めよ。
:::solution
$\langle x^2\delta_0'',\phi\rangle=(x^2\phi)''(0)=2\phi(0)$ である。従って $2\delta_0$ である。
:::
::::

::::exercise{level=basic}
$(xH)''$ を求めよ。
:::solution
積の微分で $(xH)'=H+x\delta_0=H$ である。さらに一度微分して $\delta_0$ となる。
:::
::::

::::exercise{level=basic}
$\delta_2*\delta_{-3}$ を求めよ。
:::solution
畳み込みでの作用は $\phi(2+(-3))=\phi(-1)$ なので $\delta_{-1}$ である。
:::
::::

::::exercise{level=advanced}
$\R$ 上で $T'=0$ なら $T$ は定数関数の超関数であることを示せ。
:::solution
$\int\phi=0$ なら $\psi(x)=\int_{-\infty}^x\phi(s)ds$ はコンパクト台で $\psi'=\phi$ なので $\langle T,\phi\rangle=0$ である。積分 $1$ の試験関数 $\eta$ を固定し、任意の $\phi$ から $\eta\int\phi$ を引く。差の積分がゼロなので $\langle T,\phi\rangle=\langle T,\eta\rangle\int\phi$ である。
:::
::::

::::exercise{level=advanced}
$\R$ 上の $xT=0$ の全解を求めよ。
:::solution
原点を避ける試験関数 $\phi$ には $\phi/x$ も試験関数なので、$T$ の作用はゼロである。:ref[thm-four-point-supported] により $T=\sum_{k=0}^m c_k\delta_0^{(k)}$ である。Leibniz 則から $x\delta_0^{(k)}=-k\delta_0^{(k-1)}$（$k\ge1$）である。各ジェットは切断付き単項式で独立に指定できるので $c_1,\ldots,c_m$ はゼロである。逆に $c_0\delta_0$ は解になる。
:::
::::

::::exercise{level=advanced}
$\partial^3$ の基本解を一つ求めよ。
:::solution
$E=x^2H/2$ を取る。$x^2\delta_0=x\delta_0=0$ により微分の順に $xH,H,\delta_0$ が得られる。
:::
::::

::::exercise{level=advanced}
$(\partial-a)E=\delta_0$ を満たす、台が $[0,\infty)$ に含まれる解を求め、その条件での一意性を示せ。
:::solution
$E=e^{ax}H$ と置くと $E'-aE=e^{ax}\delta_0=\delta_0$ である。差 $V$ は $(e^{-ax}V)'=0$ を満たすので前の定数の特徴付けから $e^{-ax}V$ は定数である。負の半直線でゼロになる定数はゼロなので $V=0$ である。
:::
::::

::::exercise{level=basic}
$F(x)=4x+3$ による $\delta_7''$ の引き戻しを求めよ。
:::solution
$F^{-1}(y)=(y-3)/4$ である。:ref[def-four-distribution-pullback] により、$\delta_7''$ が作用する試験関数は $\phi((y-3)/4)/4$ である。二階微分でさらに $4^{-2}$ を掛け、$y=7$ に代入すると $\phi''(1)/64$ になる。二階なので微分の符号は正で、結果は $\delta_1''/64$ である。
:::
::::

::::exercise{level=advanced}
$x^2(x-3)T=0$ を満たす超関数 $T$ の全体を求めよ。
:::solution
零点は $0,3$ で、位数はそれぞれ $2,1$ である。:ref[prop-four-distribution-multiplier-kernel] より $T=c_0\delta_0+c_1\delta_0'+d_0\delta_3$ である。これらは実際に、乗数が両点でゼロであり原点では一次導関数もゼロであるため消える。一方 $\delta_0''$ では乗数の二階導関数 $-6$ が原点で残るので、単独で加えられない。
:::
::::

::::exercise{level=advanced}
$(\partial-1)^3E=\delta_0$ を満たし、台が $[0,\infty)$ に含まれる解を求めよ。
:::solution
:ref[prop-four-causal-ode] で同じ根 $1$ に対応する三因子を畳み込む。各積分では指数因子 $e^x$ が出て、残りは $H$ の三重畳み込みなので $E=e^x x_+^2/2$ である。確認には $(\partial-1)(e^xV)=e^x V'$ を三回用いる。$V=x_+^2/2$ の三階微分は $\delta_0$ であり、$e^x\delta_0=\delta_0$ である。一意性は同命題の台を使う議論に従う。
:::
::::

## 確認テスト

::::quiz{answer="-1" tol="1e-6"}
$\chi$ は原点近くで $1$ の切断とする。$\langle\delta_0',x\chi(x)\rangle$ を求めよ。
:::explanation
作用は $-(x\chi)'(0)=-\chi(0)=-1$ である。
:::
::::

::::quiz
超関数の微分で台はどうなるか。

- [x] 元の台に含まれる。
- [ ] 常に元の台と一致する。
- [ ] 常に空になる。

:::explanation
台の外にある試験関数を微分しても台は広がらず、元の超関数の作用はゼロになる。例えば $H$ の微分では台が半直線から一点に縮む。$\delta_0'$ では一点が残る。
:::
::::

::::quiz{answer="2" tol="1e-6"}
$H*H$ の $x=2$ における連続代表元の値を求めよ。
:::explanation
二つの因子がともに $1$ となる $y$ の範囲は $0<y<2$ である。従って積分値は $2$ である。
:::
::::

::::quiz
$\R$ 上の任意の超関数 $T$ に対し、常に定義できる積はどれか。

- [x] $xT$
- [ ] $T\delta_0$
- [ ] $T^2$

:::explanation
$x$ は滑らかなので $xT$ は試験関数に $x$ を掛けて作用させれば定義できる。$T=\delta_0$ の場合、残りの二つはいずれもデルタ同士の積となり、一般の超関数に対するこの積は定義されていない。
:::
::::

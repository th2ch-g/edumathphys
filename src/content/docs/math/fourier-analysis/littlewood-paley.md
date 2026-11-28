---
title: "Littlewood–Paley 理論と Fourier 乗算子"
description: "周波数を尺度ごとに分け、二乗関数と乗算子の Lᵖ 評価を導く。"
notation:
  - symbol: 'T_m'
    meaning: "Fourier 乗算子 m の作用素"
  - symbol: '\Delta_j'
    meaning: "第 j 周波数帯の滑らかな切断"
  - symbol: 'S(f)'
    meaning: "Littlewood–Paley 二乗関数"
  - symbol: 'B^s_{p,q}'
    meaning: "Besov 空間"
  - symbol: 'F^s_{p,q}'
    meaning: "Triebel–Lizorkin 空間"
---

## この章の目標

- 有界な周波数乗数と $L^p$ 有界性の関係を整理する。
- 周波数の分割を作り、二乗関数によるノルム評価を証明する。
- Mikhlin–Hörmander と一次元 Marcinkiewicz の条件を使い分ける。
- Besov 空間と Triebel–Lizorkin 空間で尺度の和を取る順序を比較する。

## Fourier 乗算子

:ref[def-four-riesz-transform] では周波数ごとに有界な数を掛けた。この操作を一般化すると、微分方程式を解いたり、特定の周波数帯を残したりできる。$L^2$ では乗数の大きさがそのまま作用素ノルムになる。

:::definition[Fourier 乗算子]{#def-four-multiplier}
可測な $m\in L^\infty(\R^n)$ に対し $T_mf=\mathcal F^{-1}(m\widehat f)$ を $L^2$ 上で定める。$T_m$ が $L^p$ に有界に延長されるとき、$m$ を :term[Fourier 乗算子]{en="Fourier multiplier" yomi="ふーりえじょうさんし"} と呼ぶ。どの空間に対する乗算子かは指数 $p$ で指定する。
:::

:::proposition[二乗ノルムと合成]{#prop-four-multiplier-l2}
$\|T_m\|_{2\to2}=\|m\|_\infty$、$T_m^*=T_{\bar m}$、$T_mT_b=T_{mb}$ である。
:::

:::proof
Plancherel で $T_m$ を積作用素へ移す。上界は $|m|\le\|m\|_\infty$ による。任意の $\eps>0$ について $|m|>\|m\|_\infty-\eps$ の集合から正の有限測度の部分集合を取り、その指示関数を周波数側の入力とすればノルム比は $\|m\|_\infty-\eps$ 以上である。随伴と合成は積作用素上で計算できる。
:::

:::example[移動と帯域制限]
$m=e^{-2\pi ia\cdot\xi}$ なら $T_m=\tau_a$ で、全ての $p$ の等長作用素である。$m=1_{[\alpha,\beta]}$（一次元）は $L^2$ の直交射影である。$1<p<\infty$ の有界性は一次元 Riesz 変換からも従う。実際、$P_+=(I+iR_1)/2$ は正半軸への射影であり、変調で端点を移した二つの半軸射影の差が区間射影になる。端点の値は零測度なのでどちらでも同じである。
:::

:::example[二つの射影の積と重なり]
一次元で $m=1_{[1,3]}$、$b=1_{[2,4]}$ とすると $T_mT_b=T_{1_{[2,3]}}$ である。これは順番を替えても同じである。これに対して互いに交わらない二つの区間なら積作用素はゼロである。作用素を空間の核で順に積分するよりも、乗数の台を見るだけで区別できる。
:::

### 移動と可換な二乗有界作用素

乗算子から作った作用素は、周波数では平行移動と同じく積作用素になる。反対に、空間の全ての平行移動と交換する二乗有界作用素を考えると、周波数の異なる部分を混ぜられない。この性質から乗数そのものを復元する。

:::proposition[移動不変な $L^2$ 作用素の表示]{#prop-four-multiplier-translation-invariant}
$T:L^2(\R^n)\to L^2(\R^n)$ が有界な線形作用素であるとする。$T$ が全ての平行移動と可換であることと、ある $m\in L^\infty$ に対し $T=T_m$ であることは同値である。$m$ は a.e. 一意であり $\|T\|_{2\to2}=\|m\|_\infty$ である。
:::

:::proof
乗算子の向きは $\widehat{\tau_a f}=e^{-2\pi ia\cdot\xi}\widehat f$ により従う。逆に $A=\mathcal F T\mathcal F^{-1}$ とすると、任意の $a$ について $e^{-2\pi ia\cdot\xi}$ の乗算と可換である。$\psi\in C_c^\infty$ なら $\psi(\xi)=\int \widehat\psi(-a)e^{-2\pi ia\cdot\xi} da$ と書ける。$\widehat\psi\in L^1$ なので、$L^2$ ノルムで積分と $A$ を交換でき、$A$ は $\psi$ の乗算とも可換である。

有限測度の有界可測集合 $E$ に対し、$0\le\psi_j\le1$ の $C_c^\infty$ 関数で $\psi_j\to1_E$ a.e. とできる。実際、拡大した箱内で $1_E$ を平滑化して $L^1$ 収束させ、その部分列を取ればよい。有界性と優収束で乗算作用素の強収束が得られ、$A$ は $1_E$ の乗算と可換である。箱 $Q_k=[-k,k]^n$ 上で $m_k=A1_{Q_k}$ と置くと、$m_k$ は $Q_k$ に台を持つ。$l\ge k$ なら $1_{Q_k}m_l=A1_{Q_k}=m_k$ なので、全空間上の可測関数 $m$ に貼り合わせられる。任意の $E\subset Q_k$ に対し $A1_E=1_E m$ となり、$\int_E|m|^2\le\|A\|^2|E|$ である。$|m|>\|A\|+\eps$ の集合を $Q_k$ 内で試せばその測度はゼロになる。従って $m\in L^\infty$ である。有限台の単関数で $Af=mf$ が成立し、それらは $L^2$ で稠密であるため全体で成立する。ノルムと一意性は :ref[prop-four-multiplier-l2] による。
:::

:::example[有限個の遅れを組み合わせる]
一次元で $Tf(x)=2f(x)-f(x-1)-f(x+1)$ とする。これは移動と可換で、乗数は $m(\xi)=2-2\cos(2\pi\xi)=4\sin^2(\pi\xi)$ である。従って $\|T\|_{2\to2}=4$ である。$\xi=k$ では $m=0$ になるため整数付近の帯域を抑え、$\xi=k+1/2$ 付近では最大の大きさで通す。実空間での二階差分を連続微分と同一視できるのは、$\xi$ が原点に近く $m(\xi)/(2\pi\xi)^2\to1$ となる範囲である。
:::

### 可積分核と有限測度の端点

$L^2$ での乗数は有界でさえあればよい。一方、核が有限測度なら全ての $L^p$ で上界を持ち、$L^1$ で移動と可換な有界作用素はこの方法で尽くされる。ここでは関数としての核に点質量を加えることが自然である。

:::proposition[$L^1$ の移動不変作用素と有限測度]{#prop-four-multiplier-measure-kernel}
$T$ が $L^1(\R^n)$ 上で有界かつ全ての平行移動と可換なら、一意な有限複素 Radon 測度 $\mu$ があり $Tf=\mu*f$ と書ける。このとき $\|T\|_{1\to1}=|\mu|(\R^n)$ であり、$1\le p\le\infty$ で $\|\mu*f\|_p\le|\mu|(\R^n)\|f\|_p$ が成り立つ。
:::

:::proof
測度を与えた場合の上界は $\mu$ の全変動と積分形の Minkowski から従う。逆に非負で積分 $1$ の滑らかなコンパクト台近似単位 $u_j$ を取り、$Tu_j$ を $C_0$ 上の積分汎関数と見る。$\|Tu_j\|_1\le\|T\|$ なので双対単位球の弱 $*$ コンパクト性と測度の表現定理により、ある部分ネットが有限測度 $\mu$ に弱 $*$ 収束し $|\mu|(\R^n)\le\|T\|$ を満たす。任意の $f\in C_c^\infty$ に対し、$T(u_j*f)=(Tu_j)*f$ である。実際、畳み込みを $u_j$ の平行移動の積分と書き、$L^1$ で積分と $T$ を交換する。左辺は $Tf$ に $L^1$ 収束する。一方、任意の $h\in C_c^\infty$ を掛けて積分した右辺は、$y\mapsto\int f(x-y)h(x)dx\in C_c$ を $Tu_j$ に作用させた値である。このため部分ネット上の極限は $\int(\mu*f)h$ となる。従って $Tf=\mu*f$ a.e. であり、稠密性で全ての $L^1$ に移る。

測度の作用がゼロなら $\mu*u_j=0$ である。$h\in C_c$ と組み合わせると $\int h(x)(\mu*u_j)(x)dx\to\int h\,d\mu$ が $h$ の一様連続性で従うため、$\mu=0$ となる。これで一意性を得る。二つのノルム不等式を合わせると等号になる。
:::

:::example[可積分核がない等長作用素]
$m(\xi)=e^{-2\pi ia\cdot\xi}$ では $\mu=\delta_a$ が核となり、$T_mf(x)=f(x-a)$ である。$\delta_a$ を $L^1$ 関数で置き換えることはできない。実際、$L^1$ の Fourier 変換は無限遠でゼロへ向かうが $|m|=1$ である。従って「核が $L^1$ なら端点で有界」という十分条件から「端点で有界なら核が $L^1$」という逆は導けない。
:::

:::example[区間射影の端点での失敗]
$m=1_{[-1,1]}$ は $1<p<\infty$ で有界であると前に見た。もし $L^1$ でも有界なら、:ref[prop-four-multiplier-measure-kernel] と Fubini により $m$ は有限測度の連続な変換と a.e. に一致する。実際、変換が全点で非零の Schwartz 関数を入力に取れば二つの乗数を比較できる。しかし区間の内部と外部で別の定数を取るため、そのような連続代表元は存在しない。$L^\infty$ 側でも、$L^\infty\cap L^2$ 上で一様評価があれば、$L^2$ 上の自己随伴性と $L^1$–$L^\infty$ の双対ノルムで $L^1$ の評価が従い矛盾する。
:::

### 変数の変更と乗算子ノルム

空間で変調や拡大をした前後の作用を比較すると、周波数の位置と尺度を変更しても乗算子の評価は変わらない。このため一つの区間射影の評価を、全ての位置と長さに使える。

:::proposition[変調・拡大・双対]{#prop-four-multiplier-conjugations}
$1<p<\infty$ とし、$m$ が $L^p$ 乗算子であるとする。$a\in\R^n$、$b>0$ について $m(\xi-a)$ と $m(b\xi)$ は $m$ と同じ作用素ノルムを持つ。さらに $T_{\overline m}$ は $L^{p'}$ 上で $T_m$ と同じノルムを持つ。
:::

:::proof
$M_a f(x)=e^{2\pi ia\cdot x}f(x)$ と置くと $\widehat{M_af}(\xi)=\widehat f(\xi-a)$ なので $M_aT_mM_{-a}=T_{m(\cdot-a)}$ である。$D_b f(x)=b^{n/p} f(bx)$ は $L^p$ の等長写像で、$D_b^{-1}T_m D_b=T_{m(b\cdot)}$ を変数変換で確認できる。どちらも逆の共役を使えば両方向のノルム不等式が成り立つ。随伴の式は Schwartz 関数の間で Plancherel から従う。双方で稠密性を用い、双対作用素のノルムが等しいことを使う。
:::

:::example[位置と長さの異なる周波数窓]
$[0,1]$ への射影の乗数を $m$ とする。$[3,7]$ の射影は $m((\xi-3)/4)$ であり、変調と拡大で得られるため $L^p$ ノルムは $[0,1]$ の射影と等しい。$L^2$ ではどちらも直交射影でノルムは $1$ である。一方、$m((\xi-3)/4)+m((\xi+7)/4)$ のように窓を二つ足す場合には、共役だけで一つの窓のノルムへ帰着できない。複数の帯域で一様評価を出すために二乗関数を使う。
:::

### $L^p$ で必要な乗数の大きさ

$L^2$ では、乗数が大きい集合を周波数側でそのまま入力の台として使えた。一般の $L^p$ では、狭い周波数帯に集中させた滑らかな入力で代用できる。ただし大きさを制御するだけで全ての $p$ の有界性が得られるわけではない。

:::proposition[$L^p$ 作用素ノルムによる乗数の下界]{#prop-four-lp-multiplier-lower-bound}
$1<p<\infty$ で $T_m:L^p\to L^p$ が有界なら、$\|m\|_\infty\le\|T_m\|_{p\to p}$ である。
:::

:::proof
$u(x)=e^{-\pi|x|^2}$、$v(x)=u(x)^{p-1}$ とすると $\int u v=\|u\|_p\|v\|_{p'}$ である。$\xi_0$ を $m$ の Lebesgue 点とし、$\eps>0$ に対し
$$
f_\eps(x)=\eps^{n/p}u(\eps x)e^{2\pi i x\cdot\xi_0},\qquad
h_\eps(x)=\eps^{n/p'}v(\eps x)e^{2\pi i x\cdot\xi_0}
$$
とする。$L^p$ と $L^{p'}$ のノルムはそれぞれ $u,v$ のノルムに等しい。Plancherel で内積を計算し $\xi=\xi_0+\eps\eta$ と置くと
$$
\langle T_m f_\eps,h_\eps\rangle
=\int m(\xi_0+\eps\eta)\widehat u(\eta)\overline{\widehat v(\eta)}\,d\eta
$$
となる。$\widehat u\overline{\widehat v}$ は非負の Gauss 関数で、その積分は $\int u v$ である。これを拡大縮小した核による平均は、Gauss 核の尾の減衰と Lebesgue 点の局所平均により $m(\xi_0)\int u v$ に収束する。例えば $\eta$ を固定球内と外に分け、外側は $2\|m\|_\infty$ と Gauss 核の尾で抑え、内側は $m$ の平均振動と核の上限で評価できる。Hölder により絶対値は $\|T_m\|_{p\to p}\|u\|_p\|v\|_{p'}$ 以下なので、共通因子を割ると $|m(\xi_0)|\le\|T_m\|_{p\to p}$ となる。a.e. の $\xi_0$ で成り立つため結論が従う。
:::

:::example[絶対値が $1$ の二次位相でも二乗空間以外では非有界]
$m(\xi)=e^{-i\pi|\xi|^2}$ とする。:ref[prop-four-multiplier-l2] から $T_m$ は $L^2$ の等長作用素である。一方、$f_a(x)=e^{-\pi a|x|^2}$、$a>0$ に作用させると :ref[thm-four-schrodinger] の Gauss 計算により
$$
T_m f_a(x)=(1+ia)^{-n/2}\exp\left(-\frac{\pi a|x|^2}{1+ia}\right)
$$
となる。$|1+ia|=\sqrt{1+a^2}$、$\Re[a/(1+ia)]=a/(1+a^2)$ なので、Gauss 積分から
$$
\frac{\|T_m f_a\|_p}{\|f_a\|_p}
=(1+a^2)^{-n/4+n/(2p)}.
$$
$1<p<2$ なら指数が正なので $a\to\infty$ でノルム比が発散する。$p>2$ については $L^{p'}$ 上の随伴 $T_{\overline m}$ を考える。上の式で $i$ を $-i$ に替えても絶対値の計算は同じなので、$p'<2$ で非有界となる。従って元の作用素も $p>2$ で非有界である。
:::

## Littlewood–Paley 分解

原点付近の低周波を残す切断を拡大縮小すると、隣り合う切断の差が一つの環を切り取る。全部の環を足せば望遠和で元へ戻る。

:::definition[二進環の分割]{#def-four-dyadic-frequency}
実数値で球対称の $\eta\in C_c^\infty$ を $|\xi|\le1$ で $1$、$|\xi|\ge2$ で $0$、半径方向に非増加に取る。$\varphi(\xi)=\eta(\xi)-\eta(2\xi)$、$\varphi_j(\xi)=\varphi(2^{-j}\xi)$、$\Delta_j=T_{\varphi_j}$ とする。この :term[Littlewood–Paley 分解]{en="Littlewood–Paley decomposition" yomi="りとるうっどぺいりーぶんかい"} では $\varphi_j$ の台は $2^{j-1}\le|\xi|\le2^{j+1}$ に含まれる。
:::

:::proposition[復元と二乗評価]{#prop-four-lp-reproduce}
$f\in L^2$ に対し $\sum_{j\in\Z}\Delta_jf=f$ が $L^2$ で成り立つ。さらに

$$
\frac13\|f\|_2^2\le\sum_j\|\Delta_jf\|_2^2\le\|f\|_2^2 .
$$
:::

:::proof
有限和 $\sum_{j=a}^b\varphi_j=\eta(2^{-b}\xi)-\eta(2^{1-a}\xi)$ は $\xi\ne0$ で $1$ に近づき、$0\le\sum\varphi_j\le1$ である。Plancherel と優収束で復元する。各点で高々三つの帯域が重なり、非負の係数和が $1$ なので $1/3\le\sum_j|\varphi_j|^2\le1$ となる。これを $\xi$ 積分に入れる。
:::

:::example[周波数の大きさと帯域番号]
$\widehat f$ が $5<|\xi|<6$ に台を持つなら、$j=2,3$ の帯域だけが交わる。この場合 $\Delta_2 f+\Delta_3f=f$ である。各核は $2^{jn}\check\varphi(2^jx)$ であるから、第 $3$ 帯の空間尺度は第 $2$ 帯の半分になる。
:::

:::example[同じ周波数を二帯域で分担する]
例えば $\eta(3/2)=1/2$ と選んだ場合、$|\xi|=3\cdot2^{j-1}$ では $\varphi_j(\xi)=1/2$、$\varphi_{j+1}(\xi)=1/2$ になる。係数の和は $1$、二乗和は $1/2$ である。境界近くの波は二つの出力へ分かれるので、二乗和を使う場合でも各帯域が厳密に直交すると仮定できない。
:::

### 核の尺度と消えるモーメント

環の切断を実空間に戻すと Schwartz 関数になる。切断が原点の近くでゼロなので、核の積分だけでなく高いモーメントも全て消える。この性質により、滑らかな入力の低次数の Taylor 項を振り落とせる。

:::lemma[帯域核の消失モーメントと微分評価]{#lem-four-lp-band-moments}
$K=\check\varphi$、$K_j(x)=2^{jn}K(2^j x)$ とする。任意の多重指数 $\alpha$ に対して $\int x^\alpha K_j(x)dx=0$ である。$m\ge1$ を整数とし、$f$ の $m$ 階までの導関数が連続で、$m$ 階が $L^p$ に属する場合、$f\in L^p$ なら
$$
\|\Delta_j f\|_p\le C_m 2^{-jm}\sum_{|\alpha|=m}\|\partial^\alpha f\|_p
\qquad(1\le p\le\infty)
$$
である。
:::

:::proof
Schwartz 関数の積分を微分すると $\partial^\alpha\varphi(0)=\int(-2\pi i x)^\alpha K(x)dx$ である。$\varphi$ は原点の近くでゼロなので、微分は全てゼロである。$x=2^{-j}y$ で変数変換すれば $K_j$ のモーメントも消える。

$f(x-y)$ を $t\mapsto f(x-ty)$ の $t=0$ における Taylor 公式で、$m-1$ 階まで展開する。多項式部分を $K_j(y)$ で積分するとモーメントが打ち消す。積分形剰余は、$|\alpha|=m$ の導関数の移動 $\partial^\alpha f(x-ty)$ に $y^\alpha$ と $C_m(1-t)^{m-1}$ を掛けた項の和である。$L^p$ で平行移動のノルムが変わらないので Minkowski により
$$
\|\Delta_j f\|_p\le C_m\int |y|^m|K_j(y)|dy
 \sum_{|\alpha|=m}\|\partial^\alpha f\|_p .
$$
積分は $2^{-jm}\int|z|^m|K(z)|dz$ であり有限である。
:::

:::example[一次の打ち消しから二次の減衰へ]
$f(x)=e^{-\pi x^2}$ では $f''(x)=(4\pi^2x^2-2\pi)e^{-\pi x^2}$ が有界である。従って $m=2$ を用いると $\|\Delta_jf\|_\infty\le C2^{-2j}\|f''\|_\infty$ となる。$j$ を $3$ から $4$ へ増やすと、この上界の尺度因子は $1/4$ になる。実際の Gauss 関数は周波数側でも Gauss 関数なので減衰はさらに速い。Taylor による上界の利点は、Gauss 変換を知らない入力にも同じ方法を使えることである。
:::

### Bernstein の不等式

周波数の台を固定した入力には、それより少し広い切断を掛けても値が変わらない。その切断の核を使えば、帯域の上端から微分による増加を、環の下端から微分で戻せる範囲を評価できる。

:::proposition[Bernstein の帯域評価]{#prop-four-lp-bernstein}
$u\in\mathcal S'$ が $L^p$ に属し、$\operatorname{supp}\widehat u\subset\{|\xi|\le A2^j\}$ とする。$1\le p\le q\le\infty$、多重指数 $\alpha$ に対し
$$
\|\partial^\alpha u\|_q\le C 2^{j(|\alpha|+n(1/p-1/q))}\|u\|_p
$$
である。さらに $0<a<A$ で台が $a2^j\le|\xi|\le A2^j$ に含まれれば、整数 $m\ge1$ に対し $\|u\|_p\le C 2^{-2jm}\|(-\Delta)^m u\|_p$ である。
:::

:::proof
半径 $A$ の球の近傍で $1$ になる $\rho\in C_c^\infty$ を固定する。$\rho(2^{-j}\cdot)\widehat u=\widehat u$ なので $u=2^{jn}\check\rho(2^j\cdot)*u$ である。これは超関数の等式であり、右辺は $L^p$ 関数と Schwartz 関数の畳み込みだから、微分は核に移せる。$1+1/q=1/p+1/r$ とすると、微分した核の $L^r$ ノルムは $2^{j(|\alpha|+n(1-1/r))}\|\partial^\alpha\check\rho\|_r$ である。Young の不等式で最初の評価を得る。

環に台がある場合には、$a\le|\xi|\le A$ の近傍で $1$ となる環上の $\rho$ を取る。$\theta(\xi)=\rho(\xi)/(4\pi^2|\xi|^2)^m$ は原点付近でゼロと定めれば滑らかでコンパクト台を持つ。周波数の積を取ると $u=2^{-2jm}T_{\theta(2^{-j}\cdot)}(-\Delta)^m u$ である。$\check\theta$ を拡大した核の $L^1$ ノルムは $j$ に依存しないので、Young の不等式から逆の評価が従う。
:::

:::example[二乗ノルムから点の大きさを見積もる]
一次元で $\widehat u=1_{[-4,4]}$ とする。反転で $u(x)=\sin(8\pi x)/(\pi x)$、$u(0)=8$ であり、$|u(x)|\le8$ なので $\|u\|_\infty=8$ である。一方、$\|u\|_2=\sqrt8$ である。同じ計算を $[-2^j,2^j]$ に替えると $\|u\|_\infty/\|u\|_2=2^{(j+1)/2}$ である。この例では $p=2,q=\infty$ の指数 $n/2=1/2$ が実際に現れる。台の長さと Cauchy–Schwarz を使って同じ指数を直接確かめられる。
:::

:::example[低周波を含むと逆の微分評価が壊れる]
$m=1$ とし、$u_\eps(x)=u(\eps x)$ で $\widehat u$ を単位球に台を持つ非零の Schwartz 関数に取る。$0<\eps<1$ なら $\widehat u_\eps$ の台は単位球に含まれるが、$\|\Delta u_\eps\|_p/\|u_\eps\|_p=\eps^2\|\Delta u\|_p/\|u\|_p\to0$ である。従って同じ球に台を持つ入力から、$\|u_\eps\|_p\le C\|\Delta u_\eps\|_p$ という一様な逆評価は得られない。環の内側を原点から離した仮定が、逆乗数で $|\xi|^2$ を割れる理由である。
:::

### 原点に集中する超関数

$L^2$ の復元では、周波数の原点が零測度なので捨てられた。一方、超関数では一点の台に情報が集中できる。このため斉次分割では多項式を区別できなくなり、超関数全体で復元するには低周波の切断を残す。

:::proposition[斉次帯域の共通核]{#prop-four-lp-homogeneous-kernel}
$f\in\mathcal S'$ に対し、全ての $j\in\Z$ で $\Delta_j f=0$ であることと、$f$ が多項式であることは同値である。また $\widehat f$ に $\xi=0$ を含まないコンパクト台の試験関数を作用させる場合、有限部分和 $\sum_{j=a}^b\Delta_j f$ は $a\to-\infty,b\to+\infty$ で $f$ の作用を復元する。
:::

:::proof
原点を避ける試験関数 $\phi$ の台には有限個の環しか交わらず、$\phi=\sum_j\varphi_j\phi$ と有限和で書ける。従って全ての帯域がゼロなら $\widehat f$ はこうした試験関数を消し、台は原点に含まれる。有限部分和についても、添字範囲がこの有限個の環を含めば、試験関数に作用した値は $\widehat f$ と一致する。

:ref[thm-four-point-supported] により原点に台を持つ超関数は Dirac のデルタの導関数の有限和である。この等式は $C_c^\infty$ が $\mathcal S$ で稠密であることから $\mathcal S'$ の等式でもある。これを逆変換すると多項式になる。反対に多項式の変換はそのようなデルタの導関数であり、原点付近で恒等的にゼロの $\varphi_j$ を掛ければ全ての作用が消える。
:::

:::example[一次関数と点質量の違い]
一次元で $f(x)=x$ なら $\widehat f=(i/(2\pi))\delta_0'$ であり、各 $\Delta_j f=0$ である。これに対し空間側で $f=\delta_0$ とすると $\widehat f=1$ なので $\Delta_j f=K_j$ となる。空間側の点質量は各帯域に広がり、周波数側の点質量は斉次分割で消える。どちらの変数に台があるかで振る舞いが逆になる。
:::

## 二乗関数と Lᵖ ノルムの同値性

各帯域の出力を絶対値の和で測ると、同じ関数を多数の帯域に分けたときの評価が悪くなる。二乗和を使えば、周波数を無作為に符号反転した結果の平均から $L^p$ を制御できる。

:::lemma[独立符号の二乗和]{#lem-four-khintchine}
独立で等確率の符号 $\eps_j\in\{-1,1\}$ と有限個の複素数 $a_j$ に対し、$0<p<\infty$ で

$$
c_p\Bigl(\sum_j|a_j|^2\Bigr)^{1/2}
\le (\mathbb E|\sum_j\eps_j a_j|^p)^{1/p}
\le C_p\Bigl(\sum_j|a_j|^2\Bigr)^{1/2}.
$$
:::

:::proof
実係数の和 $X$ で分散を $V=\sum a_j^2$ とすると $\mathbb E e^{tX}=\prod\cosh(ta_j)\le e^{t^2V/2}$ である。最後は冪級数で $(2k)!\ge 2^k k!$ と比較した。Markov の不等式を $t=u/V$ と $-t$ に使うと $\mathbb P(|X|>u)\le 2e^{-u^2/(2V)}$ となる。尾の積分で任意の $p$ の上界を得る。複素係数では実部と虚部に適用し、$|X|>u$ なら片方が $u/\sqrt2$ を超えることを用いる。二乗平均は複素係数でも $V=\sum|a_j|^2$ である。四乗平均が $CV^2$ 以下なので

$$
V\le V/2+(\mathbb E|X|^4)^{1/2}\mathbb P(|X|^2>V/2)^{1/2}
$$

からその確率は定数倍以上である。その集合で $|X|^p\ge(V/2)^{p/2}$ なので下界を得る。$V=0$ は直接従う。
:::

:::lemma[滑らかな環の符号和]{#lem-four-signed-annuli}
任意の固定された滑らかな環上の関数 $\psi$ に対し、有限個の係数 $|\eps_j|\le1$ で作る $m=\sum_j\eps_j\psi(2^{-j}\cdot)$ は、$1<p<\infty$ で一様な乗算子ノルムを持つ。
:::

:::proof
各 $\xi\ne0$ に重なる台は有限個で個数が共通に有界なので $L^2$ の上界を得る。核の各項は $2^{jn}\check\psi(2^j x)$ で、その微分も Schwartz 減衰を持つ。$2^j|x|\le1$ の項の和には $\sum 2^{jn}\le C|x|^{-n}$ を使い、$2^j|x|>1$ では $C2^{jn}(2^j|x|)^{-n-2}$ の等比級数を使う。微分では $2^{jn}$ を $2^{j(n+1)}$ に替えれば $C|x|^{-n-1}$ となる。この評価は係数や項数によらない。:ref[prop-four-kernel-smoothness] と :ref[thm-four-cz] で結論を得る。
:::

:::theorem[Littlewood–Paley の二乗関数]{#thm-four-lp-square}
$1<p<\infty$ で :term[二乗関数]{en="square function" yomi="にじょうかんすう"} $S(f)=(\sum_j|\Delta_j f|^2)^{1/2}$ は $\|S(f)\|_p\asymp\|f\|_p$ を満たす。ここで $\asymp$ は両方向に $f$ によらない正定数の評価があることをいう。
:::

:::proof
まず有限個の $j$ に対し :ref[lem-four-khintchine] を点ごとに使い、Tonelli で積分と確率平均を交換する。:ref[lem-four-signed-annuli] が各符号和のノルムを $C\|f\|_p$ で抑えるため、単調収束で $\|S(f)\|_p\le C\|f\|_p$ を得る。

逆向きでは、$\varphi$ の台の近傍で $1$ である実数値の滑らかな環上の切断 $\psi$ を取り、$\widetilde\Delta_j=T_{\psi(2^{-j}\cdot)}$ とする。$f,g\in\mathcal S$ について Plancherel と有限和の極限から

$$
\langle f,g\rangle=\sum_j\langle\Delta_jf,\widetilde\Delta_jg\rangle.
$$

点ごとの Cauchy–Schwarz と Hölder で右辺の絶対値は $\|S(f)\|_p\|(\sum_j|\widetilde\Delta_jg|^2)^{1/2}\|_{p'}$ 以下である。上側の証明を $\psi$ にも使えるので、$g$ を $L^{p'}$ の単位球の Schwartz 関数に動かし双対ノルムを取る。一般の $f$ へは稠密性で移せる。実際、上界が差にも成り立つので帯域列が $L^p(\ell^2)$ で収束し、各成分は畳み込みで定まる $\Delta_jf$ に一致する。
:::

:::example[一つの帯域列と二つの帯域列]
$f$ の変換が一つの $\varphi_j=1$ の部分に含まれ、他の切断と交わらなければ $S(f)=|f|$ となる。二帯域で出力が $u,v$ なら $S(f)=(|u|^2+|v|^2)^{1/2}$ になる。例えば同一点で $u=3,v=4i$ なら値は $5$ である。位相の違いは二乗関数の値に直接は入らない。
:::

:::example[符号平均の四乗の具体値]
$a_1=1,a_2=2$ とする。$\eps_1+2\eps_2$ の絶対値は $1,3$ をそれぞれ等確率で取るので、四乗平均は $(1+81)/2=41$ である。二乗平均は $5$ である。二乗和に直接対応する指数は $2$ だが、任意の固定 $p$ で両者を定数倍で比較できることが :ref[lem-four-khintchine] の利点である。
:::

### 帯域の符号を変えた部分和

二乗関数の証明に使った符号は、平均を取るためだけのものではない。任意の符号で作用素ノルムが揃っているので、帯域を任意の順で足し戻せる。個々の周波数を一つずつ並べる場合には、そのような一様評価は一般に成り立たない。

:::corollary[滑らかな帯域分解の無条件収束]{#cor-four-lp-unconditional}
$1<p<\infty$ なら、任意の $f\in L^p$ について有限部分集合 $F\subset\Z$ による $\sum_{j\in F}\Delta_j f$ は $F$ が全体へ増える極限で $f$ に $L^p$ 収束する。さらに任意の有界な複素係数 $(c_j)$ に対して $\sum_j c_j\Delta_j f$ は $L^p$ で無条件に収束し、ノルムは $C_{p,n}\sup_j|c_j|\|f\|_p$ 以下である。
:::

:::proof
係数を正規化し $\sup|c_j|\le1$ とする。有限集合 $F$ に対し $\varphi$ の上で $1$ となる補助切断 $\widetilde\Delta_j$ を使うと、Schwartz 関数 $g$ との内積で
$$
\left|\left\langle\sum_{j\in F}c_j\Delta_j f,g\right\rangle\right|
\le \left\|\left(\sum_{j\in F}|\Delta_j f|^2\right)^{1/2}\right\|_p
\left\|\left(\sum_j|\widetilde\Delta_j g|^2\right)^{1/2}\right\|_{p'}
$$
である。恒等式は最初 $L^2\cap L^p$ で Plancherel により成立し、$L^p$ 近似で一般へ移る。補助切断の二乗関数は :ref[thm-four-lp-square] の上側の証明で評価できる。従って有限和のノルムは $C$ と第 1 因子の積以下であり、全和は $\|f\|_p$ で抑えられる。一方、$F$ が固定有限集合 $F_0$ の外にあるときには、第 1 因子を全ての $j\notin F_0$ の二乗和で抑える。$F_0$ を増やすとこの上界は優収束によりゼロへ近づく。従って任意の有界係数について Cauchy 条件が成り立つ。全係数 $1$ の場合は $L^2\cap L^p$ 上で $L^2$ の復元と極限が一致する。稠密性と一様評価で $L^p$ 全体でも極限は $f$ である。
:::

:::example[三つの離れた環の符号を替える]
一次元で $\eta=1$ を $|\xi|\le5/4$ で、$\eta=0$ を $|\xi|\ge7/4$ で満たすように選ぶと、$\varphi_j=1$ となる $2^j$ の近傍がある。実際、$7/8\le2^{-j}|\xi|\le5/4$ なら $\eta(2^{-j}\xi)=1$、$\eta(2^{1-j}\xi)=0$ である。その近傍に台を持つ変換から $u_0,u_3,u_6$ を作る。$f=u_0+u_3+u_6$ に対して $c_0=1,c_3=-1,c_6=1$ とすれば結果は $u_0-u_3+u_6$ となる。二乗関数はどちらも点ごとに $(|u_0|^2+|u_3|^2+|u_6|^2)^{1/2}$ である。$L^2$ ではノルムは厳密に等しく、$L^p$ では :ref[thm-four-lp-square] の定数を通して比較できる。
:::

### Poisson 半群から作る連続な二乗関数

環を離散的に選ぶ代わりに、Poisson 半群を高さ $t$ で微分して帯域を選べる。乗数 $2\pi t|\xi|e^{-2\pi t|\xi|}$ は $t|\xi|$ が小さいときに $0$ へ減り、大きいときに指数的に減る。従って $t^{-1}$ 付近の周波数を連続に測れる。

:::definition[Poisson の $g$ 関数]{#def-four-lp-poisson-g}
$P_t$ を :ref[prop-four-upper-poisson-n] の Poisson 核とする。$f\in\mathcal S$ に対し
$$
g(f)(x)=\left(\int_0^\infty|t\partial_t(P_t*f)(x)|^2\frac{dt}{t}\right)^{1/2}
$$
と置く。$dt/t$ は $t$ を定数倍しても変わらない測度である。
:::

:::lemma[Hilbert 空間値の核に対する弱型評価]{#lem-four-lp-hilbert-valued-cz}
Hilbert 空間 $E,F$ を取り、$U:L^2(\R^n;E)\to L^2(\R^n;F)$ が有界とする。コンパクト台の入力の台の外で $U$ が作用素値核 $K(x,y):E\to F$ による積分で表され、$E\to F$ の作用素ノルムで
$$
\sup_{y\ne z}\int_{|x-z|>2|y-z|}\|K(x,y)-K(x,z)\|\,dx<\infty
$$
を満たすとする。このとき $U$ は $L^1\cap L^2$ 上で弱型 $(1,1)$、$1<p\le2$ で強型 $(p,p)$ である。随伴にも同じ核条件が成り立つなら $2<p<\infty$ でも強型である。
:::

:::proof
$f\in L^1(E)\cap L^2(E)$ を取り、$\|f(x)\|_E$ の平均が $\lambda$ を超える極大二進立方体 $Q_i$ を選ぶ。:ref[thm-four-cz-decomposition] と同じ停止条件により $\sum_i|Q_i|\le\|f\|_1/\lambda$、各平均は $2^n\lambda$ 以下である。良い部分 $g$ は外で $f$、各 $Q_i$ で $f$ の $E$ 値平均とする。三角不等式で $\|g\|_\infty\le C\lambda$、$\|g\|_1\le\|f\|_1$ なので $\|g\|_2^2\le C\lambda\|f\|_1$ である。二乗評価と Chebyshev により $\{\|Ug\|>\lambda/2\}$ の測度は $C\|f\|_1/\lambda$ 以下である。

$b_i=(f-g)1_{Q_i}$ は $E$ 値平均がゼロであり、$\sum_i\|b_i\|_1\le2\|f\|_1$ である。固定倍の拡大 $Q_i^*$ の外では中心 $y_i$ における核を引いて作用させられる。よって Tonelli と仮定により $\int_{(Q_i^*)^c}\|Ub_i\|_F\le C\|b_i\|_1$ となる。拡大の和集合の測度は $C\|f\|_1/\lambda$ 以下で、外ではこの $L^1$ 評価を足して Markov の不等式を用いる。無限和は $\sum b_i$ が $L^2$ で収束することと、拡大の外で出力の和が $L^1$ に収束することで正当化できる。二つの極限は a.e. で一致する。これで弱型が得られた。

$1<p<2$ については各 $\lambda>0$ で $f=f1_{\{\|f\|>\lambda\}}+f1_{\{\|f\|\le\lambda\}}$ と分け、前者には弱型、後者には二乗評価を使う。$\|Uf\|>\lambda$ の測度は $C\lambda^{-1}\int_{\|f\|>\lambda}\|f\|+C\lambda^{-2}\int_{\|f\|\le\lambda}\|f\|^2$ 以下である。$p\lambda^{p-1}$ を掛けて積分すると、一つ目は $\int_0^{\|f(x)\|}\lambda^{p-2}d\lambda$、二つ目は $\int_{\|f(x)\|}^\infty\lambda^{p-3}d\lambda$ となり、どちらも有限で $\|f\|_p^p$ の定数倍を与える。$p>2$ は随伴に $p'$ を適用し、Hölder と $F$ 値 $L^p$ の双対ノルムで得られる。いずれも $L^2$ との共通部分の稠密性で延長する。
:::

:::proposition[連続二乗関数のノルム同値性]{#prop-four-lp-poisson-square}
$1<p<\infty$ に対し $\|g(f)\|_p\asymp\|f\|_p$ である。$p=2$ では $\|g(f)\|_2=\frac12\|f\|_2$ となる。
:::

:::proof
$\mathcal H=L^2((0,\infty),dt/t)$ とし $Uf(x)=(t\partial_tP_t*f(x))_{t>0}$ を考える。Plancherel により
$$
\|Uf\|_{L^2(\mathcal H)}^2
=\int|\widehat f(\xi)|^2
\int_0^\infty (2\pi t|\xi|)^2 e^{-4\pi t|\xi|}\frac{dt}{t}\,d\xi
=\tfrac14\|f\|_2^2
$$
である。内側で $u=2\pi t|\xi|$ に置き、$\int_0^\infty u e^{-2u}du=1/4$ を使った。

$L(x)(t)=t\partial_tP_t(x)$ と置く。Poisson 核の式を微分すると $|L(x)(t)|\le C t(t^2+|x|^2)^{-(n+1)/2}$、$|\nabla_xL(x)(t)|\le C t(t^2+|x|^2)^{-(n+2)/2}$ である。$r=|x|>0$ として $t=ru$ に置き換えると、第 1 式の $\mathcal H$ ノルムの二乗は $C r^{-2n}\int_0^\infty u(1+u^2)^{-(n+1)}du$ 以下である。第 2 式も同様に $C r^{-2n-2}$ 以下となる。差を線分に沿って評価すれば :ref[lem-four-lp-hilbert-valued-cz] の積分条件を満たす。随伴の核は $L(y-x)^*$ であり、ノルムと差に同じ評価が成り立つ。核表示は台の外で $\mathcal H$ ノルムの上界を用いて積分でき、随伴についても台が離れた試験関数の間で Fubini を使えば従う。補題によって $\|Uf\|_{L^p(\mathcal H)}\le C_p\|f\|_p$ が得られる。

逆向きでは $f,h\in\mathcal S$ に上の Plancherel 計算を偏極して $\langle f,h\rangle=4\int_{\R^n}\langle Uf(x),Uh(x)\rangle_{\mathcal H}dx$ を得る。Cauchy–Schwarz、Hölder と $p'$ での上界を使うと、右辺の絶対値は $C\|g(f)\|_p\|h\|_{p'}$ 以下である。$h$ で双対ノルムを取れば逆向きが従う。$L^p$ への延長は稠密性による。各固定 $t$ では微分核が $L^{p'}$ に属し、近似列の畳み込みは一様に収束するので、延長された $Uf$ は同じ畳み込みで a.e. に表される。
:::

:::example[高さを限った二乗エネルギー]
周波数 $\xi\ne0$ の成分が $t\le T$ に与える二乗の重みは
$$
\int_0^{2\pi T|\xi|}u e^{-2u}du
=\frac{1-(1+4\pi T|\xi|)e^{-4\pi T|\xi|}}4
$$
である。$4\pi T|\xi|=1$ なら全体の重み $1/4$ のうち $1-2/e$ の割合がこの範囲に入る。$|\xi|$ を大きくすると同じ $T$ で割合は $1$ へ向かう。高周波の変化が小さい高さで捉えられることを定数まで確認できる。
:::

### 非接錐の面積積分

Poisson の $g$ 関数では、境界の各点の真上で微分を測った。高さに比例する水平範囲にも積分を広げると、非接錐上の二乗関数になる。

:::proposition[Lusin の面積積分]{#prop-four-lusin-area}
$a>0$ とし、$u(z,t)=P_t*f(z)$ と置く。$1<p<\infty$ の $f\in L^p$ に対し
$$
S_a f(x)=\left(\int_0^\infty\int_{|z-x|<at} |t\partial_tu(z,t)|^2
\,\frac{dz\,dt}{t^{n+1}}\right)^{1/2}
$$
と定めると $\|S_a f\|_p\asymp\|f\|_p$ である。$p=2$ では $\|S_a f\|_2^2=\frac14|B(0,a)|\|f\|_2^2$ となる。この積分を :term[Lusin の面積積分]{en="Lusin area integral" yomi="るーじんのめんせきせきぶん"} という。
:::

:::proof
$z=x+tw$ と置き、$\mathcal H_a=L^2(B(0,a)\times(0,\infty),dw\,dt/t)$ に値を取る写像 $Vf(x)(w,t)=(t\partial_tu)(x+tw,t)$ と考える。ここで微分は水平変数を固定して取る。$\|Vf(x)\|_{\mathcal H_a}=S_a f(x)$ である。二乗積分で $x$ を $tw$ だけ移動しても積分値は変わらないので、:ref[prop-four-lp-poisson-square] の二乗計算の $|B(0,a)|$ 倍になる。

核を $L_a(v)(w,t)=(t\partial_tP_t)(v+tw)$ とする。$r=|v|>0$ とし、$0<t\le r/(2(a+1))$ では $|v+tw|\ge r/2$ なので $|L_a|\le Ctr^{-n-1}$、$|\nabla_v L_a|\le Ctr^{-n-2}$ である。残りの高さでは、Poisson 核の微分評価から $|L_a|\le Ct^{-n}$、$|\nabla_vL_a|\le Ct^{-n-1}$ である。これらを二乗して $dw\,dt/t$ で積分すると、$L_a$ とその勾配の $\mathcal H_a$ ノルムはそれぞれ $C_{a,n}r^{-n}$、$C_{a,n}r^{-n-1}$ 以下になる。従って核と随伴核は :ref[lem-four-lp-hilbert-valued-cz] の積分条件を満たす。離れた台の上での核表示は、これらのノルム上界と Fubini から従う。その補題で $\|S_a f\|_p\le C\|f\|_p$ を得る。稠密性で延長した作用は、各固定 $(w,t)$ で導関数核が $L^{p'}$ に入るため同じ畳み込みで表される。

逆向きでは、二乗の等式を偏極して $\langle f,h\rangle=4|B(0,a)|^{-1}\int\langle Vf(x),Vh(x)\rangle_{\mathcal H_a}dx$ とする。右辺に Cauchy–Schwarz と Hölder を用い、$Vh$ を $p'$ で評価すれば $\|f\|_p\le C\|S_a f\|_p$ となる。
:::

:::example[一次元で錐を広げた場合の二乗ノルム]
$n=1$ なら $|B(0,a)|=2a$ なので、開き $a=2$ では $\|S_2 f\|_2=\|f\|_2$ となる。開きが $1$ の場合は $\|S_1 f\|_2=\|f\|_2/\sqrt2$ である。二乗ノルムでは水平方向の積分範囲が倍になることが、そのまま二乗の係数に現れる。
:::

:::proposition[重み付き二乗関数 $g_\lambda^*$]{#prop-four-lp-g-star}
$1<p<\infty$、$\lambda>\max(1,2/p)$ とする。Poisson 拡張 $u(z,t)=P_t*f(z)$ に対し
$$
g_\lambda^* f(x)=\left(\int_0^\infty\int_{\R^n}
(1+|z-x|/t)^{-\lambda n}|t\partial_t u(z,t)|^2\,
\frac{dz\,dt}{t^{n+1}}\right)^{1/2}
$$
と置くと、$\|g_\lambda^*f\|_p\asymp\|f\|_p$ である。
:::

:::proof
開き $1$ の錐内では重みが $2^{-\lambda n}$ 以上なので、:ref[prop-four-lusin-area] で下側を得る。$p=2$ では先に $x$ で積分し、$x=z+tw$ と置くと $\int (1+|w|)^{-\lambda n}dw<\infty$ が二乗積分の定数になる。$p>2$ では非負の試験関数 $h\in L^{(p/2)'}$ を掛けて $x$ で積分する。重みの核 $t^{-n}(1+|\cdot|/t)^{-\lambda n}$ の畳み込みは :ref[prop-four-maximal-radial-envelope] により $CMh$ 以下である。従って $\int |g_\lambda^* f|^2 h\le C\int |g(f)|^2 Mh$ であり、Hölder と $M$ の強型評価、:ref[prop-four-lp-poisson-square] を組み合わせればよい。

$1<p<2$ では開き $a\ge1$ の評価を先に求める。閾値 $s>0$ に対し $E=\{S_1 f>s\}$、$E^*=\{M1_E>1/(2(a+1)^n)\}$ とする。$x\notin E^*$ かつ $|z-x|<at$ なら、$B(z,t)\subset B(x,(a+1)t)$ なので $|B(z,t)\setminus E|\ge|B(z,t)|/2$ である。従って各 $(z,t)$ について
$$
|B(z,at)\setminus E^*|\le 2a^n |B(z,t)\setminus E|.
$$
左辺が正なら $x$ をそこから選んで先の下界を使え、左辺がゼロなら式は自明である。Fubini で $\int_{(E^*)^c}(S_a f)^2\le2a^n\int_{E^c}(S_1 f)^2$ を得る。$|E^*|\le C a^n|E|$ の弱型評価と合わせて
$$
|\{S_a f>s\}|\le C a^n\left(|E|+s^{-2}\int_{E^c}(S_1 f)^2\right)
$$
となる。$ps^{p-1}$ を掛けて積分すると、第 2 項は各点で $v=S_1 f(x)$ として $\int_v^\infty pv^2 s^{p-3}ds=pv^p/(2-p)$ を与える。従って $\|S_a f\|_p^p\le C_p a^n\|S_1 f\|_p^p$ である。

最後に $|z-x|<t$ と $2^{j-1}t\le|z-x|<2^j t$（$j\ge1$）に分ければ、$(g_\lambda^* f)^2\le C_\lambda\sum_{j\ge0}2^{-j\lambda n}(S_{2^j}f)^2$ となる。$p/2<1$ のため和の $p/2$ 乗は各項の $p/2$ 乗の和以下である。上の開きの評価を使うと $\|g_\lambda^*f\|_p^p\le C\sum_{j\ge0}2^{-j(\lambda p/2-1)n}\|S_1 f\|_p^p$ となり、$\lambda>2/p$ で和が収束する。
:::

## Mikhlin–Hörmander の乗算子定理

全周波数で $m$ を微分して一様な上限を求めると、原点の近くの特異性を見落とす。代わりに各環を同じ尺度へ移し、その上で同じ滑らかさを課す。

:::theorem[Mikhlin–Hörmander]{#thm-four-mikhlin-hormander}
$m\in L^\infty$、$s>n/2$ とする。固定した環分割 $\varphi$ に対し

$$
\sup_{j\in\Z}\|m(2^j\cdot)\varphi\|_{H^s(\R^n)}<\infty
$$

なら $T_m$ は弱型 $(1,1)$ であり $1<p<\infty$ で有界である。特に $m$ が原点外で $C^k$、整数 $k>n/2$ で $|\partial^\alpha m(\xi)|\le A|\xi|^{-|\alpha|}$（$|\alpha|\le k$）なら適用できる。
:::

:::proof
$M_j(\xi)=m(2^j\xi)\varphi(\xi)$、$L_j=\check M_j$、$K_j(x)=2^{jn}L_j(2^jx)$ とする。$H^s$ の定義と反射により $\|\langle x\rangle^sL_j\|_2\le C$ である。また $\xi_l M_j$ は $M_j$ に有界台の滑らかな関数を掛けたものとしてよい。周波数の積は $x$ での畳み込みであり、$\langle x\rangle^s\le C_s\langle x-y\rangle^s\langle y\rangle^s$ と Young の不等式から $\|\langle x\rangle^s\nabla L_j\|_2\le C$ を得る。Cauchy–Schwarz で

$$
\|L_j\|_1+\|\nabla L_j\|_1\le C,\qquad
\int_{|x|>R}|L_j(x)|dx\le C R^{-(s-n/2)}\quad(R\ge1)
$$

となる。

有限和の核 $K=\sum K_j$ を考える。$2^j|h|\le1$ の項では差の $L^1$ ノルムは $|h|\|\nabla K_j\|_1\le C2^j|h|$ である。$2^j|h|>1$ では $|x|>2|h|$ の範囲で $|x-h|>|h|$ だから二つの尾の積分は $C(2^j|h|)^{-(s-n/2)}$ 以下である。両方の等比級数を足せば一様な Hörmander 定数を得る。各有限和の核は $L^1$ かつ有界連続なので台の外の核表示が成り立つ。随伴も同じ差の評価を満たす。さらに $L^2$ のノルムは $m$ の有界性と分割の非負性から一様に抑えられる。従って :ref[thm-four-cz] により一様な作用素評価が得られる。

有限和の乗数は $m$ に a.e. 収束し有界なので、$L^2$ では強収束する。$L^2\cap L^p$ に対して a.e. 収束部分列に Fatou を使うと、強型評価は極限へ移る。弱型は閾値を少し下げて同じ部分列に使えば移る。各空間での拡張は稠密性による。最後に微分条件から、固定環上の $\partial^\alpha[m(2^j\cdot)\varphi]$ は連鎖律で一様有界である。:ref[prop-four-sobolev-structure] の整数階の表示で $s=k$ の条件が従う。
:::

:::example[楕円型作用素の乗数を直接調べる]
$m(\xi)=-\xi_j\xi_k/(\xi^{\mathsf T}A\xi)$ で $A$ が正定値なら、$m$ は原点外で滑らかで次数 $0$ の斉次関数である。$\partial^\alpha m$ は次数 $-|\alpha|$ で、単位球面上では有界である。従って定数係数の :ref[thm-four-elliptic-lp] の評価を乗算子定理からも得られる。
:::

:::example[非斉次の平滑化乗数]
$n=1$ で $m(\xi)=(1+\xi^2)^{-1}$ とする。$|\xi m'(\xi)|=2\xi^2/(1+\xi^2)^2\le1/2$ である。最後は $(1+\xi^2)^2\ge4\xi^2$ による。従って一次の条件で $L^p$ 有界性を得る。この例では逆変換を直接計算しても $K(x)=\pi e^{-2\pi|x|}$ となり、$\|K\|_1=1$ なので端点を含めて評価できる。
:::

### 環の切断を選び直す

Hörmander 条件に現れる $\varphi$ は原点を避ける固定環に台を持つ。別の環の形を使っても同じ条件が得られることを確認しておくと、乗数に合わせて計算しやすい切断を選べる。

:::proposition[Hörmander 条件の切断への非依存性]{#prop-four-hormander-cutoff-invariance}
$s\ge0$、$\rho\in C_c^\infty(\R^n\setminus\{0\})$ とする。:ref[def-four-dyadic-frequency] の分割に対し
$$
\sup_j\|\rho\,m(2^j\cdot)\|_{H^s}
\le C_{\rho,s}\sup_j\|\varphi\,m(2^j\cdot)\|_{H^s}.
$$
:::

:::proof
$\rho$ の台の上では $1=\sum_{l\in F}\varphi(2^{-l}\cdot)$ と書ける。ここで $F$ はその台と交わる環を集めた有限集合である。従って
$$
\rho(\xi)m(2^j\xi)=
\sum_{l\in F}\rho(\xi)\bigl[\varphi(\cdot)m(2^{j+l}\cdot)\bigr](2^{-l}\xi).
$$
$s\ge0$ の $H^s$ では、固定された $\rho$ の積と、固定倍率の拡大は有界である。積の方は Fourier 側で $\widehat\rho$ の畳み込みとなり、$\langle x\rangle^s\le C_s\langle x-y\rangle^s\langle y\rangle^s$ から $\|\langle y\rangle^s\widehat\rho(y)\|_1$ と Young の不等式を使える。拡大の方は Fourier 変数を変えれば $\langle 2^l x\rangle^s$ と $\langle x\rangle^s$ の比が固定の $l$ で有界になる。従って各項を右辺の上限で抑え、有限和を取ればよい。
:::

:::example[一つの固定環で導関数を調べる]
$n=3$ では整数階の十分条件に $k=2$ を使える。$m(\xi)=\xi_1^2/|\xi|^2$ は $m(2^j\xi)=m(\xi)$ なので、$\varphi m(2^j\cdot)$ は $j$ に依存しない。固定した環上で $|\xi|\ge1/2$ であるため第 $2$ 階までの連続導関数は有界である。従って $H^2$ ノルムは有限となる。この議論では $\xi=0$ に値を与える必要がない。逆に全空間で $C^2$ であることを要求すると、この代表的な乗数を扱えなくなる。
:::

### 虚数冪の界と高い周波数での振動

Laplace 作用素の虚数冪の乗数は、絶対値が $1$ である。それでも $L^p$ 評価を出すには、位相が各環でどれだけ動くかを測る必要がある。

:::proposition[虚数冪の多項式型評価]{#prop-four-mikhlin-imaginary-powers}
$\tau\in\R$ に対して $m_\tau(\xi)=(4\pi^2|\xi|^2)^{i\tau}$（$\xi\ne0$）とする。整数 $k>n/2$ と $1<p<\infty$ に対し
$$
\|T_{m_\tau}\|_{p\to p}\le C_{n,p,k}(1+|\tau|)^k
$$
である。
:::

:::proof
$\log(4\pi^2|\xi|^2)=L(\xi)$ と書く。$L(r\xi)=L(\xi)+2\log r$ なので、$|\alpha|\ge1$ の微分 $\partial^\alpha L$ は次数 $-|\alpha|$ の斉次関数である。単位球面では滑らかなので、$|\partial^\alpha L(\xi)|\le C_\alpha|\xi|^{-|\alpha|}$ が従う。微分 $\partial^\alpha e^{i\tau L}$ は、$e^{i\tau L}$ と、$(i\tau)^l\prod_{r=1}^l\partial^{\beta_r}L$（$1\le l\le|\alpha|$、各 $|\beta_r|\ge1$、$\sum_r\beta_r=\alpha$）の有限和の積である。これは $\alpha=0$ から始め、各微分が指数因子を微分するか既存の因子を微分するかで帰納できる。このため $|\partial^\alpha m_\tau|\le C_\alpha(1+|\tau|)^{|\alpha|}|\xi|^{-|\alpha|}$ である。:ref[thm-four-mikhlin-hormander] の証明では Hörmander 定数と $L^2$ ノルムに対して線形なノルム上界が得られる。実際、両方の上界で作用素を割ってから弱型・強型評価を使えばよい。従って全ての微分の上限をまとめると所望の評価になる。
:::

:::example[虚数冪の位相は尺度で一定量だけ進む]
$m_\tau(2\xi)=e^{2i\tau\log2}m_\tau(\xi)$ である。例えば $\tau=\pi/(2\log2)$ なら $m_\tau(2\xi)=-m_\tau(\xi)$ となる。各環の入力の $L^2$ ノルムは変わらず、全帯域の $L^2$ ノルムも $1$ 倍である。$L^p$ の上界で $\tau$ が現れるのは、同じ環の内部での位相の変化を制御するためである。
:::

:::example[実数冪の平滑化と端点]
$\alpha>0$ について $m(\xi)=(1+4\pi^2|\xi|^2)^{-\alpha}$ とする。Gamma 積分により
$$
m(\xi)=\frac1{\Gamma(\alpha)}\int_0^\infty t^{\alpha-1}e^{-t}e^{-4\pi^2t|\xi|^2}dt .
$$
最後の指数は質量 $1$ の熱核 $G_t$ の変換である。従って $m$ の逆変換は $t^{\alpha-1}e^{-t}G_t/\Gamma(\alpha)$ を $t$ で積分した非負な $L^1$ 関数で、Tonelli により総質量は $1$ である。Fubini で変換を取れば上の式と一致する。従ってこの乗算子は全ての $L^p$ でノルム $1$ 以下である。例えば $\alpha=1$、$n=1$ では核は $\frac12e^{-|x|}$ となり、積分 $1$ を直接確かめられる。整数階の条件だけから得られる開区間の評価に加え、この乗数の正値な核から端点も評価できた。
:::

## Marcinkiewicz の乗算子定理

一次元では区間への周波数射影が制御できるため、各環で滑らかに変化する必要はない。乗数が跳ぶ場合も、跳躍の絶対値の総和を各尺度で制御すればよい。

:::theorem[一次元 Marcinkiewicz の定理]{#thm-four-marcinkiewicz-multiplier}
$m$ が有界で、正負両方の二進区間 $[2^j,2^{j+1}]$、$[-2^{j+1},-2^j]$ の各内部で全変動が一様に有界なら、$T_m$ は $L^p(\R)$（$1<p<\infty$）で有界である。
:::

:::proof
$P_I$ を区間 $I$ への周波数射影とする。任意の有限個の入力 $g_j$ について、半軸射影 $P_+$ を全ての成分に掛けると

$$
\|(\sum_j|P_+g_j|^2)^{1/2}\|_p\le C_p\|(\sum_j|g_j|^2)^{1/2}\|_p.
$$

これは :ref[lem-four-khintchine] と $P_+=(I+iR_1)/2$ の有界性を符号和に適用して得る。端点が異なる半軸へは成分ごとに異なる変調を前後に掛ける。変調は点ごとの絶対値を変えないため、同じベクトル評価が成り立つ。差を取れば任意の区間列 $I_j$ に対しても成り立つ。

$\widehat{\Delta_j f}$ の台を含む正負の環上で、$m$ を各半軸の二つの隣接二進区間に分ける。各区間で、有界変動の関数は定数と累積変動測度に分解できる。従ってその区間に制限した $m$ は $1_I$ の積分で表せる。具体的に $(a,b)$ 上の右連続代表元は $m(a+)+\mu((a,\xi])$ であり、区間全体の指示関数に定数を掛け、$1_{[t,b)}$ を $d\mu(t)$ で積分する。有限個の境界の値は任意である。各成分のこの表現に現れる複素測度の総変動は共通の $A$ 以下に取れる。

そこで各成分の変動の絶対値を $A$ で割った測度を確率測度に補い、余った確率質量には空区間を割り当てる。極分解で各積分に絶対値 $1$ 以下の因子を含める。この確率を成分ごとに組み合わせ、ベクトル値 Minkowski と先の区間射影の評価を使う。任意の有限部分について

$$
\|(\sum_j|T_m\Delta_j f|^2)^{1/2}\|_p
\le C_p A\|(\sum_j|\Delta_j f|^2)^{1/2}\|_p
$$

を得る。

$f\in\mathcal S$ について $u=T_m f\in L^2$ とし、$\langle u,g\rangle=\sum_j\langle T_m\Delta_j f,\widetilde\Delta_j g\rangle$ を $g\in\mathcal S$ に用いる。これは :ref[thm-four-lp-square] の復元の議論と同じである。上のベクトル評価と滑らかな $\widetilde\Delta_j$ に対する上界により $|\langle u,g\rangle|\le C\|f\|_p\|g\|_{p'}$ となり、双対と稠密性で結論を得る。
:::

:::example[各帯域で不連続な乗数]
正の区間 $[2^j,2^{j+1})$ で $m=(-1)^j$ とし、負の側でゼロとする。区間内部の全変動はゼロで、有界性もあるので定理を適用できる。各境界で跳躍するため、通常の微分条件を全ての $\xi\ne0$ で満たすわけではない。
:::

### 座標ごとに二進区間を選ぶ場合

多変数で球面の環だけを使うと、一つの座標だけで跳躍する乗数を扱いにくい。座標ごとに一次元の分解を行えば、長方形上の射影を使って評価できる。

:::proposition[長方形版 Marcinkiewicz 条件]{#prop-four-rectangular-marcinkiewicz}
各座標を正負の二進区間のいずれかに制限した開長方形を $R=\prod_l(a_l,b_l)$ とする。$m\in L^\infty(\R^n)$ について、各 $R$ 上で $m|_R$ が部分長方形の指示関数の積分
$$
m(\xi)1_R(\xi)=\int 1_{I_1\times\cdots\times I_n}(\xi)\,d\mu_R(I_1,\ldots,I_n)
$$
で表され、$\mu_R$ の総変動が全ての $R$ で $A$ 以下なら、$T_m$ は $1<p<\infty$ で有界である。例えば $m$ が各閉長方形まで滑らかに延長され、各 $E\subset\{1,\ldots,n\}$ について
$$
\int_{\prod_{l\in E}(a_l,b_l)}
\left|\left(\prod_{l\in E}\partial_l\right)m(t_E,a_{E^c})\right|\,dt_E\le A_0
$$
が共通に成り立てば十分である。$E=\varnothing$ の積分は頂点の値の絶対値を意味する。
:::

:::proof
後者から前者を得るには、各変数について微積分の基本定理で $m$ の値を下端での値と導関数の積分に分ける。この操作を $n$ 回繰り返すと、各 $E$ の項は下端から $\xi_l$ までの $\partial_E m$ の積分になる。これを $1_{\{t_l\le\xi_l<b_l\}}$ と残りの座標の全区間の指示関数の積として書けばよい。総変動は $2^nA_0$ 以下である。

証明には一次元の滑らかな切断を座標ごとに使い、$D_{\mathbf j}=\prod_l\Delta^{(l)}_{j_l}$ と置く。$\Delta^{(l)}$ は第 $l$ 座標だけに作用する。一次元の二乗評価は Hilbert 空間値の有限列に対しても成立する。実際、入力成分の独立符号和に :ref[lem-four-khintchine] と一次元の作用素評価を適用すれば上界が得られる。帯域添字についても独立符号を使う際は、まず一方の符号和に対する二乗和を取ってから他方に適用する。ここで Hilbert ノルムに対する符号和評価は、実数・複素数の場合と同じ二乗期待値と四乗期待値で下から抑えられる。上界は $p\ge2$ で Minkowski を $\|\sum_l|X_l|^2\|_{p/2}$ に適用し、$p<2$ なら二乗期待値から従う。これで Fubini を使い各座標へ順に二乗評価を適用でき、
$$
\|(\sum_{\mathbf j}|D_{\mathbf j}f|^2)^{1/2}\|_p\le C_{p,n}\|f\|_p
$$
となる。補助切断にも同じ上界があり、逆方向は :ref[thm-four-lp-square] と同じ双対と復元の計算で従う。座標超平面は周波数側で零測度なので、$L^2$ での復元を妨げない。

長方形への射影は各座標の区間射影の積である。前節の異なる区間のベクトル評価を、座標ごとに Fubini で多変数へ移し、順に適用する。従って帯域ごとに異なる長方形の射影にも $L^p(\ell^2)$ の共通評価がある。各 $D_{\mathbf j}$ の周波数台は高々 $4^n$ 個の二進長方形と交わるので、仮定の表示をそれぞれに使う。前節と同様に総変動を確率測度に補い、ベクトル値 Minkowski で長方形射影の評価を平均すれば
$$
\|(\sum_{\mathbf j}|T_m D_{\mathbf j}f|^2)^{1/2}\|_p
\le C_{p,n}A\|(\sum_{\mathbf j}|D_{\mathbf j}f|^2)^{1/2}\|_p
$$
を得る。補助切断で復元し双対ノルムを取ると $\|T_m f\|_p\le C_{p,n}A\|f\|_p$ となる。
:::

:::example[二次元の一象限]
$m(\xi)=1_{\{\xi_1>0,\xi_2<0\}}$ は各二進長方形内で $0$ または $1$ である。このため定理の総変動の上界は $1$ でよい。直接にも $T_m$ は第 1 座標の正半軸射影と第 2 座標の負半軸射影の積であり、Fubini で $L^p$ 有界となる。
:::

### 円周の整数周波数に移す

一次元の実軸で成立する Marcinkiewicz の条件から、円周上の整数周波数についても同じ結論を導ける。円周の三角多項式に緩やかな包絡を掛けると、実軸上では各整数周波数が狭い帯域に広がる。帯域内で乗数を一定にしておけば、円周と実軸の作用を比較できる。

:::corollary[円周上の Marcinkiewicz 条件]{#cor-four-torus-marcinkiewicz}
有界列 $(a_k)_{k\in\Z}$ について
$$
\sup_{j\ge0}\sum_{k=2^j}^{2^{j+1}-1}
  \bigl(|a_{k+1}-a_k|+|a_{-k-1}-a_{-k}|\bigr)<\infty
$$
とする。$\widehat f(k)\mapsto a_k\widehat f(k)$ は円周上の $L^p$（$1<p<\infty$）に有界な作用素を定める。
:::

:::proof
$m$ を $[k-1/2,k+1/2)$ 上で $a_k$ となる階段関数とする。実軸の二進区間内部の跳躍は、対応する整数区間での隣接差からなる。尺度 $j\ge0$ では正側の $k+1/2$ が $2^j<k+1/2<2^{j+1}$ に入るため、$2^j\le k\le2^{j+1}-1$ の差だけである。負側も同様である。$j<0$ の区間内部では $|k+1/2|\ge1/2$ の跳躍は現れない。従って :ref[thm-four-marcinkiewicz-multiplier] で $T_m$ の $L^p(\R)$ 有界性を得る。

円周の変数を $t=2\pi x$ とし、$P(x)=\sum_{k\in F}c_k e^{2\pi ikx}$ を周期 $1$ の三角多項式とする。$\widehat w\in C_c^\infty((-1,1))$ の非零な $w$ を取り、$w_\eps(x)=w(\eps x)$ とする。$0<\eps<1/2$ なら $w_\eps(x)e^{2\pi ikx}$ の変換は $(k-\eps,k+\eps)$ に含まれるので、$T_m(w_\eps P)=w_\eps P_a$ となる。ここで $P_a(x)=\sum a_kc_k e^{2\pi ikx}$ である。実軸の不等式を $p$ 乗し $\eps$ を掛けると、$P_a,P$ のそれぞれについて $\int_{\R}|w(y)|^p|P(y/\eps)|^pdy$ 型の積分を比較することになる。

任意の連続な周期 $1$ の関数 $v$ と $h\in L^1(\R)$ に対し $\int h(y)v(y/\eps)dy\to\int h\,\int_0^1v$ である。$h$ が有限区間の指示関数のときは区間を周期 $\eps$ の部分に分け、残る長さ $2\eps$ 以下の部分を $\|v\|_\infty$ で評価すればよい。一般の $h$ はその線形結合で $L^1$ 近似し、誤差は $\|v\|_\infty$ で一様に制御する。従って $\eps\downarrow0$ の極限を取って共通の $\|w\|_p^p$ を割ると、周期 $1$ 上の $L^p$ 評価を得る。$t=2\pi x$ に戻せば規約の正規化測度での評価になり、三角多項式の稠密性で延長する。
:::

:::example[整数の二進帯域だけを選ぶ]
$k=0$ では $a_0=0$、$2^{2j}\le |k|<2^{2j+1}$（$j\ge0$）では $a_k=1$、残りの $k$ では $a_k=0$ とする。各正負の二進区間では内部の隣接差はゼロで、次の帯域へ移る一箇所だけが $1$ である。従って上の条件の和は高々 $2$ である。例えば $f(t)=1+2\cos t+3\cos(2t)+4\cos(4t)$ に作用させると、$1,4$ の周波数が残って $2\cos t+4\cos(4t)$ になる。一つの区間射影を多数足して逐次評価すると個数に比例して上界が増えるが、この条件では帯域全体を一度に制御できる。
:::

### 全変動を直接積分できる場合

各帯域の変動に一様な上限を課す定理は、全空間での変動の総和が無限でも使える。一方、全変動が有限なら、半軸射影の積分だけでより直接に評価できる。

:::proposition[全変動が有限な一次元乗数]{#prop-four-marcinkiewicz-global-variation}
一次元で $m$ が有界変動なら、$1<p<\infty$ に対して
$$
\|T_m\|_{p\to p}\le |m(-\infty)|+\operatorname{Var}(m)\|P_+\|_{p\to p}.
$$
:::

:::proof
a.e. の代表元を右連続に取り、有限な複素測度 $d m$ で
$$
m(\xi)=m(-\infty)+\int_{\R}1_{[a,\infty)}(\xi)\,d m(a)
$$
と書ける。$[a,\infty)$ の射影は $M_aP_+M_{-a}$ であり、$a$ に関して強連続である。実際、任意の $f\in L^p$ に対し $M_a f$ は優収束で $L^p$ 連続であり、$P_+$ も有界である。このため $T_m f$ を射影の積分で定義でき、Minkowski で所望の評価を得る。$f\in L^2\cap L^p$ では $L^2$ 内積を試験して Fubini を使い、周波数ごとの積作用素と一致することを確かめられる。
:::

:::example[三つの高さを持つ階段]
$m(\xi)=2$（$\xi<0$）、$m(\xi)=-1$（$0\le\xi<1$）、$m(\xi)=3$（$\xi\ge1$）とすると、$m=2-3\,1_{[0,\infty)}+4\,1_{[1,\infty)}$ である。従って $T_m=2I-3P_++4M_1P_+M_{-1}$ であり、上界は $2+7\|P_+\|_{p\to p}$ となる。$L^2$ では :ref[prop-four-multiplier-l2] から厳密なノルムは $3$ である。一方で正側の各二進区間で交互に符号が変わる先の例では全変動が無限なので、この直接評価には代えられない。
:::

## Besov 空間と Triebel–Lizorkin 空間（概観）

低周波をまとめた非斉次分割を用いる。$D_0=T_\eta$、$D_j=\Delta_j$（$j\ge1$）とすれば、望遠和で $\eta+\sum_{j\ge1}\varphi_j=1$ である。$\mathcal S'$ 上では各 $D_j$ を滑らかな乗数によって定義できる。

:::definition[二つの尺度ノルム]{#def-four-besov-triebel}
$s\in\R$、$1\le p<\infty$、$1\le q\le\infty$ とし、

$$
\|f\|_{B^s_{p,q}}=\| (2^{js}\|D_j f\|_p)_{j\ge0}\|_{\ell^q},\qquad
\|f\|_{F^s_{p,q}}=\left\|\|(2^{js}|D_j f(x)|)_{j\ge0}\|_{\ell^q}\right\|_p
$$

とする。有限となる $\mathcal S'$ の元の空間をそれぞれ :term[Besov 空間]{en="Besov space" yomi="べそふくうかん"}、:term[Triebel–Lizorkin 空間]{en="Triebel–Lizorkin space" yomi="とらいべるりぞるきんくうかん"} という。$q=\infty$ では列の上限を取る。
:::

:::proposition[基本的な一致と包含]{#prop-four-besov-cases}
$B^s_{2,2}=F^s_{2,2}=H^s$ が同値ノルムで成立する。$1<p<\infty$ なら $F^0_{p,2}=L^p$ である。また $1\le p<\infty$ で $B^s_{p,p}=F^s_{p,p}$ が同一ノルムで成立する。
:::

:::proof
二乗の場合と $p=q$ の場合の両ノルムの一致は Tonelli による。二乗の場合、Plancherel で重みは $\sum_{j\ge0}2^{2js}|\varphi_j(\xi)|^2$ となる。ただし第 $0$ 項の $\varphi_0$ は $\eta$ と読む。台の重なりが有界で、各非零項で $2^j$ と $\langle\xi\rangle$ が定数倍の範囲で一致し、切断の和が $1$ なのでこの重みは $\langle\xi\rangle^{2s}$ と比較できる。$\mathcal S'$ で分割ノルムが有限なら、各帯域の変換は局所的に平方可積分な関数である。有限重なりによりそれらの和は局所 $L^2$ で $\widehat f$ に一致し、同じ重みの比較を使える。

$F^0_{p,2}$ では $D_0$ を有限個の環の符号和に加えても、その核は $L^1$ なので :ref[thm-four-lp-square] の上側の証明は保たれる。逆向きでは $D_0$ の台の近くで $1$ の補助切断を加えて同じ再構成を行う。一般の超関数についても、分割の部分和は試験関数との作用で元に戻る。実際、周波数で余りは無限遠に逃げ、試験関数との積は全 Schwartz セミノルムでゼロに近づく。二乗関数が $L^p$ なら双対試験で $f$ 自身が有界な $L^{p'}$ 汎関数を定め、$L^p$ 関数に一致する。
:::

:::example[尺度の指数を読む]
$f$ が一つの帯域 $j$ だけで $D_jf=f$ を満たすなら両ノルムは $2^{js}\|f\|_p$ である。例えば $j=3,s=2$ の重みは $64$ であり、高周波側を低周波側より強く測る。一方 $s<0$ では同じ帯域を軽く測る。
:::

:::example[二つの尺度の和と空間側の積分]
有限な二成分 $u_0,u_1$ が空間で互いに素な台を持ち、$\|u_0\|_p=\|u_1\|_p=1$ とする。列ノルムを先に取ると $\|(|u_0|^q+|u_1|^q)^{1/q}\|_p=2^{1/p}$ だが、積分ノルムを先に取ると $(\|u_0\|_p^q+\|u_1\|_p^q)^{1/q}=2^{1/q}$ である。$p=2,q=4$ ならそれぞれ $\sqrt2,\sqrt[4]2$ になる。これは二つの定義の積分順序の計算例であり、ここでは $u_j$ を同一関数の帯域成分と仮定していない。
:::

### 別の滑らかな切断で測る

滑らかな分割を選び直すと隣の帯域の間で成分が移る。各帯域で重なる成分の個数は有界なので、比較に必要になるのは、近い尺度にある畳み込みの評価である。

:::proposition[分割を替えたときのノルム比較]{#prop-four-lp-besov-cutoffs}
:ref[def-four-dyadic-frequency] の条件を満たす別の $\eta$ から非斉次分割を作っても、$B^s_{p,q}$（$1\le p<\infty$、$1\le q\le\infty$）のノルムは同値である。$F^s_{p,q}$ でも $1<p<\infty$、$1\le q\le\infty$ なら同じことが成り立つ。
:::

:::proof
新しい分割を $E_j$ と書く。$E_j f=\sum_{|j-k|\le2} E_j D_k f$ は超関数として成立する。周波数で $D_k$ の切断の和が $1$ であり、$E_j$ の台に交わるものがその範囲に限られるためである。各 $E_j$ の核は $L^1$ で一様有界で、$2^{js}\le C_s2^{ks}$ がその範囲で成り立つ。Young の不等式で $\|E_j D_k f\|_p\le C\|D_k f\|_p$ を使い、添字を高々二つずつずらした列を足せば $B$ の比較が従う。分割を交換すれば逆向きも得られる。

$F$ の比較では $|E_j g|\le A_j|g|$ と評価する。ここで $A_j$ は $W_j(x)=C2^{jn}(1+2^j|x|)^{-n-1}$ との畳み込みである。各 $W_j$ は偶関数で $L^1$ ノルムが一定であり、:ref[prop-four-maximal-radial-envelope] の球平均の評価から $\sup_j A_j h\le C'Mh$（$h\ge0$）となる。有限列 $(g_j)$ に対し $(A_jg_j)$ の $L^p(\ell^q)$ 評価を示せば、帯域を有限回ずらして加えることで $F$ の比較が従う。

$q=p$ なら Young の不等式を成分ごとに使って $p$ 乗で足す。$1\le q<p$ では Jensen により $(A_j|g_j|)^q\le C_q A_j(|g_j|^q)$ である。非負の $h\in L^{(p/q)'}$ を掛けて積分し、核の偶性を用いると
$$
\int h\sum_j (A_j|g_j|)^q
\le C_q\int (Mh)\sum_j |g_j|^q .
$$
$M$ の $L^{(p/q)'}$ 有界性と双対ノルムから所望の評価を得る。$q=\infty$ なら $\sup_j A_j|g_j|\le C'M(\sup_j|g_j|)$ として同じ結論を得る。残る $p<q<\infty$ では、有限列の $L^p(\ell^q)$ の双対を取り、偶な核の同じ作用素を $L^{p'}(\ell^{q'})$ に適用する。$q'<p'$ なので、すでに示した範囲である。有限列で定数は個数に依存しないため、単調収束で任意の列に移る。二つの分割を交換すれば $F$ のノルムも同値である。
:::

### ノルムの順序による包含

空間積分と尺度の和を入れ替えるとき、どの方向に不等式が得られるかは指数に依存する。これは両方の空間に別々の評価が必要になる理由を説明する。

:::proposition[Besov と Triebel–Lizorkin の基本包含]{#prop-four-lp-besov-triebel-inclusions}
$1\le p<\infty$、$1\le q\le\infty$ に対し
$$
B^s_{p,\min(p,q)}\subset F^s_{p,q}\subset B^s_{p,\max(p,q)}
$$
が連続な包含として成り立つ。
:::

:::proof
$u_j=2^{js}|D_j f|$ と置く。$p\le q<\infty$ のとき、点ごとに $\ell^q$ ノルムは $\ell^p$ ノルム以下である。従って $\|u\|_{L^p(\ell^q)}\le\|u\|_{\ell^p(L^p)}$ を得る。他方、$\ell^{q/p}$ の三角不等式を積分形で $(u_j^p)_j$ に適用すると
$$
\left(\sum_j(\int u_j^p)^{q/p}\right)^{p/q}
\le \int\left(\sum_j(u_j^p)^{q/p}\right)^{p/q}.
$$
$p$ 乗根を取れば $\|u\|_{\ell^q(L^p)}\le\|u\|_{L^p(\ell^q)}$ である。$q=\infty$ でも $\sup_j\int u_j^p\le\int\sup_j u_j^p$ を用いれば同じ方向である。

$q\le p$ では $L^{p/q}$ の三角不等式を $\sum_j u_j^q$ に使い、$\|u\|_{L^p(\ell^q)}\le\|u\|_{\ell^q(L^p)}$ を得る。点ごとの $\ell^p\le\ell^q$ から逆側 $\|u\|_{\ell^p(L^p)}\le\|u\|_{L^p(\ell^q)}$ も得られる。無限列に対しては有限部分で示して単調収束を使う。いずれも $u_j$ の定義を戻せば結論になる。
:::

:::example[二成分が重なると等号が出る場合]
$p=2,q=4$ で $u_0=u_1=1_{[0,1]}$ とすると、$L^2(\ell^4)$ のノルムは $2^{1/4}$ であり、$\ell^4(L^2)$ のノルムも $2^{1/4}$ である。一方 $\ell^2(L^2)$ ノルムは $\sqrt2$ である。台が互いに素な前の例では $L^2(\ell^4)$ と $\ell^2(L^2)$ が一致した。二つの例を比べると、尺度の和の違いが成分の空間上の重なりに依存することが分かる。
:::

### 一階差による $0<s<1$ の Besov 空間

微分が存在しない入力にも、平行移動した差なら意味がある。小さい距離 $t$ での差を $t^s$ と比較すると、帯域の重み $2^{js}$ と同じ尺度を測れる。

:::definition[移動の連続度]{#def-four-lp-modulus}
$f\in L^p$ に対し $\omega_p(f,t)=\sup_{|h|\le t}\|\tau_h f-f\|_p$ とする。$1\le p<\infty$ では移動の $L^p$ 連続性により $t\downarrow0$ でゼロへ向かう。
:::

:::proposition[差分による非整数階 Besov ノルム]{#prop-four-lp-besov-difference}
$0<s<1$、$1\le p<\infty$、$1\le q\le\infty$ とする。$B^s_{p,q}$ の元は $L^p$ 関数を代表元に持ち、同値ノルムは
$$
\|f\|_p+\left\|\bigl(2^{ks}\omega_p(f,2^{-k})\bigr)_{k\ge0}\right\|_{\ell^q}
$$
である。$q<\infty$ では列ノルムを $(\int_0^1(t^{-s}\omega_p(f,t))^q\,dt/t)^{1/q}$ に置き換えてもよい。$q=\infty$ では $0<t\le1$ の上限を取る。
:::

:::proof
$A_j=2^{js}\|D_jf\|_p$ とする。$\|(A_j)\|_{\ell^q}$ が有限なら、$A_j$ は一様有界であり、$\sum_j\|D_jf\|_p\le C_s\|(A_j)\|_{\ell^q}$ である。従って $\sum D_j f$ は $L^p$ で収束し、超関数としては分割の復元により $f$ である。$|h|\le2^{-k}$ に対し、Bernstein の $p=q$ の評価と線分積分を用いて
$$
\|\tau_h D_jf-D_jf\|_p\le C\min(1,2^{j-k})\|D_jf\|_p
$$
となる。$j=0$ でも固定半径の低周波切断で同じ評価が成り立つ。各 $j$ を足して $2^{ks}$ を掛けると
$$
2^{ks}\omega_p(f,2^{-k})
\le C\left(\sum_{j\le k}2^{-(k-j)(1-s)}A_j
+\sum_{j>k}2^{-(j-k)s}A_j\right).
$$
$s$ と $1-s$ が正なので、係数は両側で可和である。列の Young の不等式、すなわち各移動列の $\ell^q$ ノルムに係数を掛けて足す評価により右辺の列ノルムは $C\|(A_j)\|_{\ell^q}$ 以下である。

逆向きでは $j\ge1$ の核の積分がゼロなので
$$
\|D_j f\|_p\le \int |K(y)|\,\omega_p(f,2^{-j}|y|)\,dy .
$$
$h$ を長さ $2^{-j}$ 以下の等しい移動の和に分けると、$\omega_p(f,2^{-j}|y|)\le(1+|y|)\omega_p(f,2^{-j})$ である。$K$ は Schwartz 関数だから、この重みを掛けても積分は有限である。従って $A_j\le C2^{js}\omega_p(f,2^{-j})$ であり、第 $0$ 項は Young の不等式で $C\|f\|_p$ 以下となる。これで同値性が従う。積分版は $2^{-k-1}\le t\le2^{-k}$ で $\omega_p(f,2^{-k-1})\le\omega_p(f,t)\le\omega_p(f,2^{-k})$ と $t^{-s}\asymp2^{ks}$ による。列の最初の項は $2\|f\|_p$ 以下であるから、積分版と比較しても制御できる。
:::

:::example[指示関数における滑らかさの境界]
$f=1_{[0,1]}$ とする。$0\le |h|\le1$ なら二つの区間の対称差の長さは $2|h|$ であり、$\omega_p(f,t)=(2t)^{1/p}$（$0<t\le1$）である。従って $1<p<\infty$、$0<s<1$ の範囲で $f\in B^s_{p,q}$ は $s<1/p$、または $s=1/p$ かつ $q=\infty$ と同値である。例えば $p=2$ では $B^{1/2}_{2,\infty}$ に属するが $B^{1/2}_{2,2}=H^{1/2}$ には属さない。跳躍の位置近くで長さ $t$ の範囲に差が集中するため、点ごとの差は小さくならなくても $L^p$ の差は $t^{1/p}$ で減る。
:::

:::example[低周波入力を一階差で測る]
$f(x)=e^{-\pi x^2}$ に対し $f'(x)=-2\pi x e^{-\pi x^2}$ であり、$\|f'\|_1=2$ である。線分に沿って積分すれば $\omega_1(f,t)\le2t$ となる。従って $s=1/2$ では尺度列は高々 $2\cdot2^{-k/2}$ であり、その $\ell^1$ 和は $2/(1-2^{-1/2})$ 以下である。この計算は Fourier 変換の減衰を使わずに $B^{1/2}_{1,1}$ の有限性を確かめている。
:::

### 連続代表元を得るための尺度和

各帯域で $L^p$ ノルムを一様ノルムへ移すには、Bernstein の因子 $2^{jn/p}$ を払う必要がある。滑らかさの指数がこれより大きいなら等比級数で足せる。臨界指数では尺度方向の絶対和が残る。

:::corollary[Besov 空間から連続関数への埋め込み]{#cor-four-lp-besov-continuity}
$1\le p<\infty$ とする。$s>n/p$ の $B^s_{p,q}$（全ての $q$）、または $B^{n/p}_{p,1}$ の元には無限遠でゼロへ向かう連続な代表元があり、一様ノルムは対応する Besov ノルムで抑えられる。
:::

:::proof
:ref[prop-four-lp-bernstein] で $\|D_j f\|_\infty\le C2^{jn/p}\|D_j f\|_p$ となる。$s>n/p$ なら $\sum_j2^{-j(s-n/p)}$ が有限なので各 $2^{js}\|D_j f\|_p$ の上限を使える。$s=n/p$ では $\ell^1$ ノルムを仮定したため同じ絶対和が有限である。各 $D_j f$ は $L^p$ 関数であり、より広い帯域の Schwartz 核と畳み込んで同じ関数になる。$C_c$ 関数で $L^p$ 近似すると畳み込みは一様収束し、近似した畳み込みは連続で無限遠でゼロへ向かう。このため各 $D_j f$ は同じ性質を持つ。さらに帯域和が一様に収束するため、その和も同じ性質を持つ。超関数としての復元によりこの和は $f$ の代表元である。
:::

:::example[臨界指数で尺度の絶対和が果たす役割]
$n=1$ とし、実偶の非負な $\widehat u\in C_c^\infty$ が $|\xi|=1$ の近傍の二つの小区間に台を持ち、$u(0)=1$ となるよう正規化する。$\eta$ を :ref[cor-four-lp-unconditional] の例のように選べば、$u(2^jx)$ は第 $j$ 帯域だけに属する。$f_N(x)=\sum_{j=1}^N u(2^jx)/j$ と置くと、$2^{j/p}\|D_j f_N\|_p=\|u\|_p/j$ である。従って $q>1$ で $B^{1/p}_{p,q}$ ノルムは $N$ に一様有界だが、$f_N(0)=\sum_{j=1}^N1/j\to\infty$ である。よって臨界指数で $q>1$ にした場合、同じ一様ノルムの評価を得られない。
:::

## 演習問題

::::exercise{level=basic}
$m=2\,1_{[-1,1]}$ の $L^2$ 乗算子ノルムを求めよ。
:::solution
$\|m\|_\infty=2$ であり、:ref[prop-four-multiplier-l2] で等しい。変換が区間内にある非零入力を取れば作用は入力の二倍になる。
:::
::::

::::exercise{level=basic}
$\Delta_j$ の核が $2^{jn}\check\varphi(2^jx)$ であることを確認し、その $L^1$ ノルムを求めよ。
:::solution
逆積分で $\xi=2^j u$ と置くと $2^{jn}$ が出て指数の変数が $2^jx$ になる。絶対値の積分で $y=2^jx$ と置けば係数が相殺され $\|\check\varphi\|_1$ となる。
:::
::::

::::exercise{level=basic}
$\widehat f$ が $12<|\xi|<14$ に台を持つとき、関与し得る帯域番号を求めよ。
:::solution
$[2^{j-1},2^{j+1}]$ が $(12,14)$ と交わるのは $j=3,4$ である。$j\le2$ なら上端は $8$ 以下、$j\ge5$ なら下端は $16$ 以上である。
:::
::::

::::exercise{level=basic}
独立な符号 $\eps_1,\eps_2$ について $\mathbb E|3\eps_1+4i\eps_2|^2$ を求めよ。
:::solution
各符号の二乗は $1$ であり、実部と虚部の二乗は $9$ と $16$ になる。和は常に $25$ なので平均も $25$ である。
:::
::::

::::exercise{level=advanced}
$n=1$ で $m(\xi)=|\xi|^{i\beta}$（$\xi\ne0,\ \beta\in\R$）が $1<p<\infty$ の乗算子であることを示せ。
:::solution
$m=e^{i\beta\log|\xi|}$ なので $|m|=1$、$m'=i\beta m/\xi$ である。$k=1>1/2$ で :ref[thm-four-mikhlin-hormander] の条件を満たす。
:::
::::

::::exercise{level=advanced}
$m(\xi)=1_{[1,2]}(\xi)+2\,1_{[-8,-4]}(\xi)$ が $L^p$ 乗算子であることを二通り示せ。
:::solution
二区間ごとに半軸射影の差とすれば、変調と一次元 Riesz 変換の有界性から従う。また $m$ は各二進区間の内部で定数なので全変動はゼロである。:ref[thm-four-marcinkiewicz-multiplier] も適用できる。
:::
::::

::::exercise{level=advanced}
$p=q$ で二つの非斉次ノルムが一致することを直接示せ。
:::solution
$\|f\|_{F^s_{p,p}}^p=\int\sum_j2^{jsp}|D_j f|^p$ である。Tonelli で順序を交換すると $\sum_j2^{jsp}\int|D_j f|^p=\|f\|_{B^s_{p,p}}^p$ である。無限値も含め等式が成り立つ。
:::
::::

::::exercise{level=advanced}
$f\in L^2$ に対し任意の符号列 $\eps_j$ の $\sum_j\eps_j\Delta_j f$ が $L^2$ で収束することを示せ。
:::solution
各周波数に高々三項しか寄与せず、乗数の絶対値は $\sum_j\varphi_j=1$ 以下である。有限部分和の乗数は $\xi\ne0$ で各点で収束する。差の二乗に $4|\widehat f|^2$ の支配を使い、Plancherel で $L^2$ 収束を得る。
:::
::::

::::exercise{level=basic}
$T=3I-\tau_1-\tau_{-1}$ を一次元で考える。乗数と $L^2$ ノルムを求めよ。
:::solution
移動の変換により $m(\xi)=3-2\cos(2\pi\xi)$ である。値域は $[1,5]$ であり、最大値 $5$ に任意に近い点の集合は正の測度を持つ。従って :ref[prop-four-multiplier-l2] により $\|T\|_{2\to2}=5$ である。
:::
::::

::::exercise{level=basic}
$n=1$ で $\widehat u=1_{[-2,2]}$ とする。$u(0)$、$\|u\|_2$、$u'(0)$ を求めよ。
:::solution
逆積分で $u(0)=4$ である。Plancherel により $\|u\|_2^2=\int_{-2}^2 1\,d\xi=4$ なので $\|u\|_2=2$ である。台が有界なので逆積分を微分でき、$u'(0)=\int_{-2}^2 2\pi i \xi\,d\xi=0$ となる。
:::
::::

::::exercise{level=advanced}
$1<p<\infty$、$0<s<1$ とする。$f=1_{[0,2]}$ に対する差分版の $\ell^1$ 項が有限となる条件と、その値を求めよ。
:::solution
$0<t\le1$ で移動した区間との対称差は $2t$ なので $\omega_p(f,t)=(2t)^{1/p}$ である。従って列の $k$ 項は $2^{1/p}2^{k(s-1/p)}$ である。等比級数は $s<1/p$ で収束し、和は $2^{1/p}/(1-2^{s-1/p})$ である。$s\ge1/p$ では項がゼロに近づかず発散する。
:::
::::

::::exercise{level=advanced}
Poisson の $g$ 関数で、$t\in[T,\infty)$ の部分が周波数 $\xi\ne0$ に与える二乗の重みを求めよ。
:::solution
$v=2\pi t|\xi|$ と置けば重みは $\int_{2\pi T|\xi|}^\infty v e^{-2v}dv$ である。原始関数は $-(2v+1)e^{-2v}/4$ なので、値は $(1+4\pi T|\xi|)e^{-4\pi T|\xi|}/4$ である。$T=0$ では全体の $1/4$ を回収する。
:::
::::

## 確認テスト

::::quiz{answer="25" tol="1e-6"}
ある点で二つの帯域出力が $3,4$、ほかがゼロである。二乗関数の平方を求めよ。
:::explanation
定義から各出力の絶対値二乗の和である。従って $9+16=25$ である。
:::
::::

::::quiz
$L^2$ で有界な可測関数の乗数ノルムを決める量はどれか。

- [x] 本質上限
- [ ] $L^1$ ノルム
- [ ] 原点での値

:::explanation
Plancherel が作用を周波数側の積へ移す。積のノルムは本質上限で抑えられ、本質上限に近い場所を入力の台に取ると同じ値に近づく。原点のような零測度集合の値は関与しない。
:::
::::

::::quiz{answer="16" tol="1e-6"}
$s=2$、$j=2$ に対する非斉次尺度の重み $2^{js}$ を求めよ。
:::explanation
積 $js=4$ なので $2^4=16$ である。
:::
::::

::::quiz
非斉次分割で $H^s$ と同値なノルムを持つのはどれか。

- [x] $B^s_{2,2}$
- [x] $F^s_{2,2}$
- [ ] 任意の $s$ に対する $F^s_{2,2}=L^2$ という等式

:::explanation
二乗積分と二乗和は Tonelli で交換でき、Plancherel で重みは $\langle\xi\rangle^{2s}$ と比較できる。従って両方が $H^s$ を与える。$s>0$ では重みが無限遠で増えるため、全ての $L^2$ 関数の変換がこの重みで二乗可積分になるわけではない。
:::
::::

::::quiz{answer="0.25" tol="1e-6"}
$f\ne0$ かつ $f\in L^2(\R^n)$ とする。Poisson の $g$ 関数について $\|g(f)\|_2^2/\|f\|_2^2$ を求めよ。
:::explanation
$\xi\ne0$ で $v=2\pi t|\xi|$ と置くと、周波数ごとの二乗の重みは $\int_0^\infty v e^{-2v}dv=1/4$ である。Plancherel でこれを $|\widehat f|^2$ に掛けて積分するので、比は $1/4$ となる。
:::
::::

---
title: "局所コンパクト Abel 群上の Fourier 解析"
description: "Haar 測度と双対群から Fourier 反転、Bochner、Pontryagin、Wiener の各定理を証明する。"
notation:
  - symbol: 'G'
    meaning: "局所コンパクト Hausdorff Abel 群"
  - symbol: '\widehat G'
    meaning: "連続指標からなる双対群"
  - symbol: 'dx'
    meaning: "G の Haar 測度"
  - symbol: 'd\chi'
    meaning: "Plancherel の等式で正規化した双対 Haar 測度"
  - symbol: '\lambda(f)'
    meaning: "L² 上の畳み込み作用素"
---

## この章の目標

- 双対群の位相を定め、有限群、円周、整数群、実空間を同じ規約で扱う。
- 畳み込み作用素から双対 Haar 測度の正規化を導く。
- 正定値関数を正測度の変換として表し、二重双対を同定する。
- Fourier 変換の零点と平行移動の稠密性の関係を証明する。

## 双対群

円周上の指数 $e^{ik\theta}$ も、Euclid 空間上の指数 $e^{2\pi i x\cdot\xi}$ も、和を積に移す連続写像である。これを任意の局所コンパクト Abel 群に置き換える。$G$ は加法で書き、Haar 測度 $dx$ を一つ固定する。可換性と :ref[prop-meas-haar-inversion] により $x\mapsto-x$ でも測度は保存される。

:::definition[局所コンパクト群の双対]{#def-four-lca-dual}
連続準同型 $\chi:G\to\{z\in\C:|z|=1\}$ を :term[指標]{en="character" yomi="しひょう"} と呼ぶ。指標全体を $\widehat G$ とし、点ごとの積で群とする。位相は $G$ の各コンパクト集合上一様収束する位相である。この位相を入れた群を :term[双対群]{en="dual group" yomi="そうついぐん"} と呼ぶ。
:::

:::proposition[双対群の局所コンパクト性]{#prop-four-dual-lca}
$\widehat G$ は局所コンパクト Hausdorff Abel 群である。
:::

:::proof
コンパクト上一様収束は群の積と逆元を連続にし、一点集合で値を比べれば Hausdorff である。$G$ の単位元のコンパクト近傍 $K$ を取り、

$$
U=\{\chi:\sup_{x\in K}|\chi(x)-1|\le1\}
$$

とする。これは $\widehat G$ の閉近傍である。整数 $m$ に対し、$x,2x,\ldots,mx\in K$ となるよう $x$ を十分小さい近傍に取る。$\chi\in U$ なら $\chi(x)=e^{it}$ の主偏角 $t$ は $|t|\le\pi/3$ であり、$jt$ の主偏角も $|jt|\le\pi/3$ の範囲に入るよう補正できる。実際 $|t|\le\pi/3$ なので連続する倍数がこの弧から外へ出るとき、円を回り込んで弧に戻ることはできない。帰納的に $|jt|\le\pi/3$ がそのまま成り立ち、$|t|\le\pi/(3m)$ である。よって $U$ は単位元で一様同程度に連続であり、群の法則で各点でも同じ近傍を移動して使える。

$U$ の指標を単位円の $G$ 個の積に入れる。Tychonoff の定理により任意のネットは点ごとに収束する部分ネットを持つ。その極限は準同型で、一様な連続性の評価を引き継ぐので連続である。同じ評価とコンパクト集合の有限被覆を使えば収束はコンパクト上一様である。$K$ 上の不等式も保たれるため極限は $U$ にある。従って $U$ はコンパクトである。
:::

:::example[整数と円周の双対]
$\Z$ は離散群なので、任意の指標は $\chi_t(k)=e^{ikt}$ である。有限集合上一様収束は $k=1$ での収束と同値なので、双対は円周 $\mathbb T$ である。一方、$\mathbb T$ 上の指標を $\R$ に引き戻すと $e^{iax}$ の形で、$2\pi$ 周期性から $a\in\Z$ である。$\R$ 上でこの形を得るには、原点近くで主偏角を連続に取る。十分狭い区間では加法性がそのまま偏角の加法性になり、有理数倍と連続性から偏角は線形である。準同型性で全実軸へ移る。円周の指標 $\chi_k$ は $k\ne0$ なら $\sup_\theta|\chi_k(\theta)-1|=2$ なので、双対の位相は離散である。
:::

:::example[Euclid 空間と積群]
$\R^n$ の指標は座標に分ければ $\chi_\xi(x)=e^{2\pi i x\cdot\xi}$ である。$\xi$ の収束はコンパクト上一様収束を与える。逆に各単位線分上で $\chi_{\xi_j-\xi}\to1$ が一様なら、その座標差がゼロへ近づく。座標差が $1/2$ 以上なら線分上に値 $-1$ が現れ、小さい差なら端点で比較できるためである。有限直積の双対も各因子の双対の直積である。指標は各因子に分解し、コンパクト集合の座標射影はコンパクトであるため位相も一致する。
:::

### 位相の判定とコンパクト開部分群

双対群の位相を調べるとき、まず指標を全て書き出し、次にどのコンパクト集合で値を比べればよいかを考える。離散群側では有限個の値で近傍を記述できる。一方、コンパクト群側では、群全体で十分近い二つの指標は一致する。以下の補題では両者を使い分ける。

:::lemma[離散群とコンパクト群の基本対応]{#lem-four-lca-discrete-compact}
$G$ が離散なら $\widehat G$ はコンパクトであり、$G$ がコンパクトなら $\widehat G$ は離散である。
:::

:::proof
離散群のコンパクト集合は有限集合である。従って $\widehat G$ は積位相を入れた $\prod_{x\in G}\mathbb S^1$ の部分空間であり、$\chi(x+y)=\chi(x)\chi(y)$ と $\chi(0)=1$ は閉条件である。離散群上の写像は常に連続であるから、これらの閉条件だけで $\widehat G$ が得られる。積のコンパクト性により結論が従う。

次に $G$ がコンパクトとする。非自明な指標 $\chi$ には $\chi(x)\ne1$ となる $x$ がある。$0<|\arg\chi(x)|\le\pi/3$ なら、$m=\lfloor\pi/(3|\arg\chi(x)|)\rfloor+1$ として $mx$ で値を取ると主偏角の絶対値が $(\pi/3,2\pi/3]$ に入る。最初から $|\arg\chi(x)|>\pi/3$ なら $x$ 自身でよい。従って非自明な $\chi$ は $\sup_{y\in G}|\chi(y)-1|>1$ を満たす。群全体がコンパクトであるため、$\sup_{y\in G}|\chi(y)-1|<1$ は開近傍を与え、そこには自明な指標しか入らない。他の各点も移動で開集合になる。
:::

:::example[無限個の二値座標]
$G=\prod_{j\ge1}\Z/2\Z$ に積位相を入れる。任意の指標は値が $\pm1$ であり、連続性によりその核は開である。従ってある有限集合 $F$ について、$x_j=0$（$j\in F$）となる $x$ は全て核に入る。$G$ をこの核で割った指標は有限個の座標だけに依存する。有限群で各基底の値を選ぶと、指標はちょうど
$$
\chi_F(x)=(-1)^{\sum_{j\in F}x_j}
$$
の形になる。$F=\{1,3\}$ と $E=\{3,5\}$ なら積は $\chi_{\{1,5\}}$ である。共通座標の寄与が二回現れて消えるためである。双対群は有限台の列の群 $\bigoplus_{j\ge1}\Z/2\Z$ で、位相は離散である。積を直和へ替えると、生成元の像は独立に選べるため、今度は全ての二値列が指標を与える。
:::

:::example[整数進列で与えるコンパクト群]
$p$ を素数とし、$\Z_p$ を $x=(x_j)_{j\ge1}$、$x_j\in\Z/p^j\Z$、$x_{j+1}\equiv x_j\bmod p^j$ を満たす列の群とする。積位相を入れると整合条件で定まる閉部分群なのでコンパクトである。任意の指標は連続性により、ある $j$ に対して $\{x:x_j=0\}$ を $|\chi-1|<1$ の弧へ移す。部分群の像がこの弧に入るなら、上の冪の議論によって像は自明である。従って指標は第 $j$ 座標で因子化し、
$$
\chi(x)=\exp(2\pi i a x_j/p^j),\qquad a\in\Z/p^j\Z
$$
となる。第 $j+1$ 座標で書き直すと添字は $pa$ である。従って双対群は $\Z[1/p]/\Z$ で離散である。例えば $p=3$ で $a/p^j=2/9$、$x_2=4$ なら位相は $8/9$ 周となる。第 $3$ 座標で $x_3=13$ としても $6\cdot13/27$ は $8/9$ と整数しか違わない。
:::

有限群の指標では、積を取って閉じていることと、指標が点を分離することが特に重要であった。上の無限積でも、有限座標の指標の積は有限座標の指標であり、異なる二つの列は必ず一つの座標で区別できる。一般の LCA 群では座標表示がないため、点の分離は次節で畳み込み作用素を使って証明する。

## Fourier 変換と Plancherel の定理

空間上の測度はすでに選んだが、双対側の Haar 測度にも定数倍の自由度がある。この定数を、畳み込み作用素のスペクトルから得る二乗ノルムによって決める。

:::definition[群上の Fourier 変換]{#def-four-lca-transform}
$f\in L^1(G)$ に対し $\widehat f(\chi)=\int_G f(x)\overline{\chi(x)}\,dx$ と定める。$(f*g)(x)=\int f(x-y)g(y)dy$、$f^*(x)=\overline{f(-x)}$ とする。Fubini と Haar 不変性により $\widehat{f*g}=\widehat f\,\widehat g$、$\widehat{f^*}=\overline{\widehat f}$ である。
:::

:::lemma[畳み込み代数の指標]{#lem-four-lca-spectrum}
$L^1(G)$ の非零連続乗法汎関数は $f\mapsto\widehat f(\chi)$ に限られる。その点ごとの収束位相は $\chi$ のコンパクト上一様収束の位相と一致する。
:::

:::proof
乗法汎関数 $\Lambda$ と $\Lambda(h)\ne0$ となる $h$ を取る。$\tau_x f(y)=f(y-x)$ とすると $(\tau_x f)*h=f*(\tau_x h)$ なので $\Lambda(\tau_x f)=b(x)\Lambda(f)$、$b(x)=\Lambda(\tau_x h)/\Lambda(h)$ である。移動の合成から $b$ は準同型である。$L^1$ で移動が連続なので $b$ も連続である。移動連続性は $C_c(G)$ では一様連続性と台の共通なコンパクト性から従い、一般には :ref[prop-meas-radon-lp-density] で近似する。移動でノルムが一定なので $b$ は有界であり、正負の整数倍を取って $|b(x)|=1$ を得る。$f*h=\int f(x)\tau_x h\,dx$ を使って作用させると $\Lambda(f)=\int f b$ である。積分と作用の交換は $x$ のコンパクト範囲で $x\mapsto\tau_x h$ の像を有限近似し、外側を $\|h\|_1\int|f|$ の尾で評価すれば正当化できる。

コンパクト上一様収束なら、$L^1$ の尾を切って汎関数の点ごとの収束を得る。逆に汎関数の収束を仮定し、$\Lambda(h)\ne0$ を固定する。共通ノルム上界 $1$ により、点ごとの収束は $L^1$ のコンパクト集合 $\{\tau_x h:x\in K\}$ 上で一様になる。分母 $\Lambda(h)$ と近似列の分母を使って比を取れば $K$ 上の一様収束が従う。
:::

:::lemma[畳み込み作用素のスペクトル表示]{#lem-four-lca-cstar}
$\lambda(f)$ を $u\mapsto f*u$ の $L^2(G)$ 上の作用素とする。これらの作用素のノルム閉包 $\mathcal C$ は、$\lambda(f)\mapsto\widehat f$ を通じて $C_0(\widehat G)$ と等長な $*$-同型をなす。
:::

:::proof
$\lambda(f)^*=\lambda(f^*)$ は内積の積分を交換すれば従う。可換なノルム閉代数 $\mathcal C$ の各元は正規である。従って Gelfand 変換が等長である。ここで用いる事実は、可換 Banach 代数のスペクトル半径が指標値の絶対値の上限であることと、正規作用素でスペクトル半径が作用素ノルムに等しいことである。後者については、$A^*A$ が自己随伴なので $\|(A^*A)^2\|=\|A^*A\|^2$ である。正規性から $(A^2)^*A^2=(A^*A)^2$ であり、従って $\|A^2\|=\|A\|^2$ となる。倍々に繰り返して半径公式を使えばよい。半径公式は :ref[thm-fa-spectral-radius]、指標によるスペクトルの表示は :ref[thm-fa-maximal-characters] による。非単位的なら単位化の追加指標は $\mathcal C$ 上でゼロなので上限に影響しない。$\mathcal C$ の指標を $\lambda(L^1)$ に制限すると :ref[lem-four-lca-spectrum] の指標になる。

逆に任意の $\chi$ に対し $|\widehat f(\chi)|\le\|\lambda(f)\|$ を示す。非負の $C_c$ 関数 $u$ で $\|u\|_1=1$ を取り、任意のコンパクト $K$ を $x_1+V,\ldots,x_l+V$ で覆う。$V$ は $\|\tau_v u-u\|_1<\eps$ となる単位元近傍である。$\tau_{x_1}^{k_1}\cdots\tau_{x_l}^{k_l}u$ を $0\le k_j<M$ で平均して $v_M$ と置く。各 $\tau_{x_j}$ は平均の両端の面だけを変え、差の $L^1$ ノルムは $2/M$ 以下である。$\tau_v$ は平均と交換するので $\sup_{x\in K}\|\tau_xv_M-v_M\|_1\le\eps+2/M$ である。$w_M=\chi\sqrt{v_M}$ と置くと $\|w_M\|_2=1$ であり、

$$
\|\tau_x w_M-\overline{\chi(x)}w_M\|_2^2
\le\|\tau_x v_M-v_M\|_1 .
$$

$\lambda(f)w_M-\widehat f(\chi)w_M$ は、この差に $f(x)$ を掛けて積分したものである。$K$ 上の上界と $K^c$ の $L^1$ 尾で評価し、$K$、$\eps$、$M$ を選べば任意に小さくできる。これで所望の不等式と指標の延長が従う。指標上の収束は $\lambda(L^1)$ の稠密性と各汎関数のノルム上界により :ref[lem-four-lca-spectrum] の位相に一致する。

一般の Gelfand 変換は指標空間上の $C_0$ 関数である。コンパクト性は単位化の指標空間が双対単位球の閉集合であることから従い、無限遠の追加指標で $0$ を取る。変換の像は等長性により閉で、共役に閉じ、指標を分離し、各指標で全てがゼロになることはない。従って Stone–Weierstrass を一地点コンパクト化に適用すれば像は $C_0(\widehat G)$ 全体である。
:::

この結果から Riemann–Lebesgue の主張 $\widehat f\in C_0$ が従う。さらに $\widehat f=0$ なら $\lambda(f)=0$ である。単位元の各近傍内に台を持つ非負の質量 $1$ の $C_c$ 関数を取り、それらとの畳み込みが $L^1$ で $f$ に近づくことを移動連続性で示せる。各関数は $L^2$ にも属するため畳み込みはゼロであり、従って $f=0$ である。指標が $G$ の点を分離することも分かる。もし $x\ne0$ で全指標が $\chi(x)=1$ を満たせば $\widehat{\tau_x f}=\widehat f$ となるが、$x$ の近傍と単位元の近傍を分けるコンパクト台関数 $f$ を取れば $\tau_x f\ne f$ である。

:::theorem[双対 Haar 測度と Plancherel]{#thm-four-lca-plancherel}
$\widehat G$ 上に一意な Haar 測度 $d\chi$ があり、Fourier 変換は $L^2(G,dx)\to L^2(\widehat G,d\chi)$ のユニタリ写像に拡張される。$f,\widehat f\in L^1$ なら $f(x)=\int_{\widehat G}\widehat f(\chi)\chi(x)d\chi$ が a.e. に成り立つ。
:::

:::proof
$g\in C_c(G)$ に対し、正汎関数 $a\mapsto\langle a g,g\rangle$ を :ref[lem-four-lca-cstar] で $C_0(\widehat G)$ に移す。:ref[thm-meas-rmk] により有限正測度 $\mu_g$ が得られる。$h\in C_c(G)$ とすると $\mu_{\lambda(h)g}=|\widehat h|^2\mu_g$ である。これは試験関数 $a$ を $\lambda(h)^*a\lambda(h)$ に替える計算で分かる。$\lambda(h)g=\lambda(g)h$ だから

$$
|\widehat h|^2\mu_g=|\widehat g|^2\mu_h .
$$

各指標には $\widehat g$ がゼロでないような $g\in C_c(G)$ を選べる。例えば非負非零の $C_c$ 関数にその指標を掛ける。$\widehat g\ne0$ の開集合で $\nu=|\widehat g|^{-2}\mu_g$ と定めれば、上の恒等式で重なり上の値が一致する。そこでコンパクト台の各試験関数を有限個の開集合に従って分割し、局所積分を足す。値は分割によらず正の汎関数になるので、Riesz の表現定理で Radon 測度 $\nu$ を得る。局所的な一致から $\mu_g=|\widehat g|^2\nu$ である。

$g$ を $\chi_0g$ に替えると、$\widehat{\chi_0g}(\chi)=\widehat g(\chi\chi_0^{-1})$ となる。また作用素を変調 $u\mapsto\chi_0u$ で共役すると同じ双対群の移動が現れる。具体的に $M_{\chi_0}^{-1}\lambda(f)M_{\chi_0}=\lambda(\overline{\chi_0} f)$ なので $\mu_{\chi_0g}$ は $\mu_g$ を $\chi\mapsto\chi\chi_0$ で移した測度である。局所表示を比較すれば $\nu$ は移動不変になる。

$\mu_g$ の全質量は $\|g\|_2^2$ である。実際、台が単位元へ縮む質量 $1$ の近似単位 $u_V$ は $g*u_V\to g$ を $L^2$ で満たす。$\widehat u_V\to1$ は双対群のコンパクト集合上で一様である。これには双対群のコンパクト集合の指標が $G$ 上で同程度連続なことを使う。それは $K$ 上の上限が $1$ 未満となる指標の近傍を移動した有限被覆と、:ref[prop-four-dual-lca] の $U$ の同程度連続性に従う。$|\widehat u_V|\le1$ と有限 Radon 測度の正則性で $\int\widehat u_V\,d\mu_g\to\mu_g(\widehat G)$ である。一方、これは $\langle \lambda(u_V)g,g\rangle\to\|g\|_2^2$ に等しい。従って $\nu$ は非零の Haar 測度であり、$\int|\widehat g|^2d\nu=\|g\|_2^2$ である。一意性は Haar 測度の定数倍の一意性（:ref[thm-meas-haar-unique]）から従う。

稠密性により等長写像を拡張し、その像は閉じている。像は $\widehat g$ への $\widehat h$ の乗算に閉じるため、:ref[lem-four-lca-cstar] により任意の $C_0(\widehat G)$ の乗算に閉じる。像に直交する $v\in L^2(\widehat G)$ を取ると、$\widehat g\bar v$ の積分は任意の $C_c$ を掛けてもゼロである。従って $\widehat g\bar v=0$ a.e. である。各コンパクト集合を有限個の $\widehat g\ne0$ の集合で覆えばそこで $v=0$ a.e. となり、$L^2$ 積分の正則性で全体でもゼロである。従って像は全体である。

$f\in L^1\cap L^2$ での一致は $C_c$ による両ノルムでの近似から従う。切断と有限測度上での指示関数の正則近似によりその近似を取れる。$f,\widehat f\in L^1$ とすると、右辺の逆積分 $u$ は有界連続である。$h\in C_c(G)$ に対し $f*h\in L^2$ で、その変換は $\widehat f\,\widehat h$ に等しい。これは $f$ を $L^1$ と $L^2$ の共通部分で $L^1$ 近似し、畳み込みの $L^2$ 評価で極限を取って得る。$\int \widehat f\,|\widehat h|^2d\nu=\langle f*h,h\rangle$ となり、Fubini で $\langle u*h,h\rangle$ も同じ値である。ここで $u$ の畳み込みと内積は $h$ がコンパクト台なので定義できる。$h$ について偏極すれば任意の $h,k\in C_c$ に対して $\langle(f-u)*h,k\rangle=0$ であり、近似単位の局所 $L^1$ 収束で $f=u$ a.e. を得る。局所収束は各コンパクト集合の近傍で $f-u$ を切断してから適用できる。
:::

:::example[測度の対応]
円周の正規化測度 $dt/(2\pi)$ には整数群の数え上げ測度が対応する。有限 $G$ の正規化測度には有限双対群の数え上げ測度が対応する。実空間の Lebesgue 測度には :ref[thm-four-plancherel] により双対の Lebesgue 測度が対応する。一般に $dx$ を $a\,dx$ に替えると $\widehat f$ は $a$ 倍になるので、$d\chi$ は $a^{-1}d\chi$ に変わる。
:::

:::example[整数群の二つの点の二乗エネルギー]
$\Z$ の数え上げ測度で $f(0)=2$、$f(1)=-1$、ほかはゼロとする。変換は $2-e^{-it}$ なので $|\widehat f(t)|^2=5-4\cos t$ である。円周の正規化測度で積分すると $5$ となり、整数上の $2^2+(-1)^2$ に一致する。$dt$ を双対測度に選ぶと余分な $2\pi$ が出る。
:::

:::corollary[群上の Hausdorff–Young 不等式]{#cor-four-lca-hausdorff-young}
双対 Haar 測度を :ref[thm-four-lca-plancherel] に従って選ぶ。$1\le p\le2$ なら $L^1\cap L^p(G)$ の $f$ に対し $\|\widehat f\|_{p'}\le\|f\|_p$ であり、Fourier 変換は $L^p(G)\to L^{p'}(\widehat G)$ に一意に延長される。
:::

:::proof
$p=1$ では指標の絶対値が $1$ なので積分の三角不等式を使う。$p=2$ は Plancherel の等式である。$1<p<2$ では、まず $f$ と双対試験関数をコンパクト台の単関数に取る。両群でその台を含む有限測度集合に制限し、$L^1\to L^\infty$ と $L^2\to L^2$ の間に :ref[thm-meas-riesz-thorin] を適用する。切断後も作用素ノルムは両端で $1$ 以下なので、$1/p=1-\theta/2$ として $\|\widehat f\|_{p'}\le\|f\|_p$ を得る。出力の積分領域をコンパクト集合で試すだけで全ノルムが得られるのは Haar 測度の正則性による。これで両群を可算個のコンパクト集合で覆えるという仮定を置かずに済む。単関数は $L^p$ で稠密であり、$L^1\cap L^p$ では両ノルムで近似できる。$L^1$ による変換の一様収束と $L^{p'}$ 極限は一致するため、延長は積分で定義した変換を保つ。
:::

### コンパクト群で見える正規化

一般論では正汎関数から測度を組み立てた。コンパクト群の場合には指標の直交性だけでも同じ正規化を計算できる。そこで先に質量を $1$ に固定する。

:::proposition[コンパクト群の指標展開]{#prop-four-lca-compact-expansion}
$G$ がコンパクトで $\int_G dx=1$ なら、$\widehat G$ の双対 Haar 測度は数え上げ測度である。任意の $f\in L^2(G)$ は
$$
f=\sum_{\chi\in\widehat G}\widehat f(\chi)\chi,\qquad
\|f\|_2^2=\sum_{\chi\in\widehat G}|\widehat f(\chi)|^2
$$
を満たす。無限和は有限部分集合上の和の極限であり、第 1 式は $L^2$ で成り立つ。
:::

:::proof
:ref[lem-four-lca-discrete-compact] により双対群は離散である。定数関数 $1$ の変換は自明な指標で $1$ となる。他の指標 $\chi$ について $\chi(y)\ne1$ を選べば、Haar 不変性により $\int\overline{\chi(x)}dx=\overline{\chi(y)}\int\overline{\chi(x)}dx$ となるため、変換値は $0$ である。Plancherel を $1$ に使うと自明な指標の質量は $1$ と決まる。Haar 不変性で全ての点の質量が $1$ である。任意の有限部分和の Fourier 変換は対応する部分集合への切断なので、Plancherel により復元と誤差の二乗和が従う。
:::

:::example[三進整数上の剰余類]
$\Z_3$ の正規化 Haar 測度を取り、$f(x)=1_{\{x_2=4\}}$ とする。各剰余類は平行移動で互いに移るため質量は $1/9$ である。従って $\chi_{a/9}$ における変換値は $e^{-2\pi i 4a/9}/9$ である。第 $2$ 座標で因子化しない指標 $\chi$ には $\{x:x_2=0\}$ 内に $\chi(y)\ne1$ となる $y$ があり、$f$ は $y$ の移動で不変なので $\widehat f(\chi)=0$ である。二乗和は $9(1/9)^2=1/9$ となり、$\|f\|_2^2$ と一致する。反転式は、$x_2-4\equiv0\bmod9$ かどうかを九つの根の和で判定する式である。
:::

:::example[二値関数のスペクトル]
$G=\prod_{j\ge1}\Z/2\Z$ 上で $u(x)=1_{\{x_1=x_2=0\}}$ とする。$2u-1$ は値が $\pm1$ なので二乗ノルムは $1$ である。一方、
$$
2u-1=\frac{-1+\chi_{\{1\}}+\chi_{\{2\}}+\chi_{\{1,2\}}}{2}
$$
である。積 $\frac14(1+\chi_{\{1\}})(1+\chi_{\{2\}})$ が $u$ になることから式を得る。四つの変換値の絶対値は全て $1/2$ であり、二乗和は $4/4=1$ となる。第 $3$ 座標を含む項が現れないのは、第 $3$ 座標の移動で $u$ が変わらないためである。
:::

### 測度の変換と一意性

畳み込みの入力には、可積分関数のほかに点質量も使える。点質量はコンパクト群でも離散群でも定義できるが、Haar 測度に関する密度を持つかどうかは群による。そこで変換を最初から測度の積分として定める。

:::definition[Fourier–Stieltjes 変換]{#def-four-fourier-stieltjes}
有限複素 Radon 測度 $\mu$ に対し $\widehat\mu(\chi)=\int_G\overline{\chi(x)}\,d\mu(x)$ と置く。これを :term[Fourier–Stieltjes 変換]{en="Fourier–Stieltjes transform" yomi="ふーりえすてぃるちぇすへんかん"} という。$\mu*\nu$ は加法写像 $(x,y)\mapsto x+y$ による $\mu\otimes\nu$ の像である。
:::

:::proposition[測度変換の積と単射性]{#prop-four-stieltjes-injective}
$\widehat\mu$ は有界連続で $\|\widehat\mu\|_\infty\le|\mu|(G)$ を満たす。$\widehat{\mu*\nu}=\widehat\mu\,\widehat\nu$ であり、$\mu\mapsto\widehat\mu$ は単射である。
:::

:::proof
上界は指標の絶対値が $1$ であることによる。任意の $\eps>0$ に対して $|\mu|(G\setminus K)<\eps$ のコンパクト集合 $K$ を選ぶ。$K$ 内の積分はコンパクト上一様収束で連続であり、外側の差は $2\eps$ 以下である。積の式は $\chi(x+y)=\chi(x)\chi(y)$ と Fubini による。

$\widehat\mu=0$ と仮定し $u\in C_c(G)$ を取る。$\mu*u(x)=\int_G u(x-y)d\mu(y)$ は $L^1$ に属し、Tonelli により $\|\mu*u\|_1\le |\mu|(G)\|u\|_1$ である。その変換は $\widehat\mu\,\widehat u=0$ なので、この節の $L^1$ 変換の単射性で $\mu*u=0$ a.e. となる。近似単位 $u_V\ge0$ で $\int u_V=1$ とし、$h\in C_c(G)$ に対し
$$
\int_G h(x)(\mu*u_V)(x)\,dx
=\int_G\int_G h(y+z)u_V(z)\,dz\,d\mu(y)
$$
を取る。$h$ は全空間で一様連続である。実際、台のコンパクト近傍で有限被覆を使えば移動の差が一様に小さく、近傍の外では小さい移動に対し差がゼロである。従って右辺は $\int h\,d\mu$ に近づく。左辺は常にゼロなので $\mu=0$ である。
:::

:::example[点質量では無限遠へ減衰しない]
$G=\R$、$\mu=2\delta_0+\delta_2$ なら $\widehat\mu(\xi)=2+e^{-4\pi i\xi}$ である。特に $\xi=k/2$ で変換値は $3$ となる。$\mu*\mu=4\delta_0+4\delta_2+\delta_4$ であり、その変換は $4+4e^{-4\pi i\xi}+e^{-8\pi i\xi}$ となる。密度が $L^1$ の場合には Riemann–Lebesgue により変換は無限遠でゼロへ向かうため、この例の点質量を $L^1$ 密度として表せないことが変換からも分かる。
:::

### Plancherel による平行移動の閉包

$L^1$ では一点での変換評価が連続である。$L^2$ の変換では零集合上の値を任意に変えても同じ元である。従って平行移動が生成する部分空間も、台を測度の意味で捉える必要がある。

:::proposition[二乗可積分な平行移動の生成空間]{#prop-four-lca-l2-cyclic}
$f\in L^2(G)$ とし $E=\{\chi:\widehat f(\chi)\ne0\}$ とする。$\tau_x f$（$x\in G$）の線形包の $L^2$ 閉包は、変換が $E^c$ で a.e. ゼロとなる関数全体である。
:::

:::proof
$\widehat{\tau_x f}(\chi)=\overline{\chi(x)}\widehat f(\chi)$ は $L^1\cap L^2$ 上で変数変換により成り立ち、$L^2$ 近似で全体に移る。従って閉包の像は $E$ 上に台を持つ。逆に $g$ が全ての移動に直交すると、Plancherel により
$$
0=\langle\tau_x f,g\rangle
=\int_{\widehat G}\overline{\chi(x)}\widehat f(\chi)\,\overline{\widehat g(\chi)}\,d\chi
$$
である。Hölder で積は $L^1(\widehat G)$ に属する。右辺はこの関数の変換を $\iota(x)$ で評価したものであり、ここで二重双対の全射性を仮定する必要はない。任意の $h\in L^1(G)$ を掛けて $x$ について積分すると、:ref[lem-four-lca-cstar] によって $C_0(\widehat G)$ に稠密な $\widehat h$ を積の前に掛けられる。従って積は a.e. ゼロである。これは $E$ 上で $\widehat g=0$ a.e. と同値であるから、直交補空間を二回取ると結論が従う。
:::

:::example[一点で消える変換の二つの挙動]
$\R$ 上で $f=1_{[-1/2,1/2]}$ とすると $\widehat f(\xi)=\sin(\pi\xi)/(\pi\xi)$ であり、零点は $0$ 以外の整数である。この集合の Lebesgue 測度はゼロなので、平行移動の線形包は $L^2$ で稠密である。一方、$L^1$ での閉包では $\widehat g(1)=0$ という制約が消えない。例えば幅を $1/2$ にした箱の $\xi=1$ での変換値は $\sin(\pi/2)/\pi=1/\pi$ なので、その箱は $L^1$ の閉包に入らない。
:::

## Bochner の定理

指標に非負の重みを掛けて積分すると、平行移動した点の間の行列が正になる。逆にこの行列の正値性だけで、正の周波数測度を復元できる。

:::definition[正定値関数]{#def-four-positive-definite}
連続関数 $b$ が任意の $x_1,\ldots,x_m\in G$、$c_1,\ldots,c_m\in\C$ に対し $\sum_{j,k}c_j\overline{c_k}\,b(x_j-x_k)\ge0$ を満たすとき、:term[正定値関数]{en="positive definite function" yomi="せいていちかんすう"} という。
:::

:::theorem[Bochner]{#thm-four-bochner}
連続関数 $b$ が正定値であることと、$\widehat G$ 上のある有限正 Radon 測度 $\mu$ で $b(x)=\int\chi(x)d\mu(\chi)$ と表せることは同値である。$\mu$ は一意で、全質量は $b(0)$ である。
:::

:::proof
積分表示があれば二次形式は $\int|\sum_jc_j\chi(x_j)|^2 d\mu$ になる。連続性には、双対群のコンパクト集合に属する指標の同程度連続性と有限測度の正則性を使う。

逆に、記号 $e_x$ の有限線形結合に $\langle e_x,e_y\rangle=b(x-y)$ を入れる。正値性で零長さの部分空間を商に取り、完備化して Hilbert 空間を得る。移動 $e_x\mapsto e_{x+t}$ は内積を保存するのでユニタリ作用素 $U_t$ となる。$v=e_0$ とすると $\langle U_t v,v\rangle=b(t)$ である。$b$ の連続性から $\|U_t e_x-e_x\|^2=2(b(0)-\Re b(t))\to0$ であり、稠密性で $U_t$ は強連続である。$\Pi(f)=\int f(x)U_{-x}dx$ は $\|\Pi(f)\|\le\|f\|_1$ の $*$-表現になる。積分は各ベクトルに作用させて定め、コンパクト上の連続性と $L^1$ の尾を使って近似できる。積と随伴は Fubini と $U_{-x}^*=U_x$ から従う。

各 $\Pi(f)$ は正規である。準同型は単位化すると可逆元を可逆元へ移すため、スペクトル半径は $L^1$ の単位化のスペクトル半径以下である。:ref[lem-four-lca-spectrum] で後者は $\|\widehat f\|_\infty$ である。従って $\|\Pi(f)\|\le\|\widehat f\|_\infty$ であり、:ref[lem-four-lca-cstar] により $C_0(\widehat G)$ に表現を延長できる。$a\mapsto\langle\Pi(a)v,v\rangle$ は正である。実際、非負の $a\in C_0$ を $\sqrt a\,\sqrt a$ と書けば作用はノルムの二乗になる。Riesz の表現定理でこれを有限正測度 $\mu$ で表す。$f\in C_c(G)$ を使うと

$$
\int_G f(x)b(-x)dx=\int_{\widehat G}\widehat f(\chi)d\mu(\chi).
$$

Fubini と連続性により $\int\overline{\chi(x)}d\mu(\chi)=b(-x)$ が全点で成り立つ。$x$ を $-x$ に替えれば必要な表示になる。一意性は $\widehat f$ の $C_0$ での稠密性から従い、質量は $x=0$ を代入すれば分かる。
:::

:::example[三角関数と正定値性]
$\mathbb T$ 上で $2+\cos(3t)$ は指標 $0$ に質量 $2$、$\pm3$ に $1/2$ ずつを置いた正測度から来るので正定値である。$1-2\cos t$ は $b(0)=-1$ で、$m=1,c_1=1$ の条件だけで正定値ではないと分かる。
:::

:::example[正定値関数は実部が負でもよい]
$\R$ で $b(x)=e^{2\pi ix}\,e^{-\pi x^2}$ とする。Gauss 変換を $1$ だけ移せば、Bochner 測度は密度 $e^{-\pi(\xi-1)^2}$ を持つ正測度である。従って $b$ は正定値で、質量は $1$ である。ところが $x=1/2$ では $b(x)=-e^{-\pi/4}$ である。正値性を課す対象は各点の実数値ではなく、点の差から作る二次形式である。
:::

### 正値性の行列判定と閉性

正定値関数を探すときは、まず二点だけで成立する必要条件を確かめる。二点で矛盾しない候補については、積分表示や畳み込みから任意の点数の正値性を示す。

:::proposition[正定値関数の基本操作]{#prop-four-positive-definite-operations}
正定値な $b$ について $b(0)\ge0$、$b(-x)=\overline{b(x)}$、$|b(x)|\le b(0)$ である。また
$$
|b(x)-b(y)|^2\le 2b(0)\bigl(b(0)-\Re b(x-y)\bigr).
$$
非負係数による和、点ごとの積、連続な点ごとの極限は正定値性を保つ。
:::

:::proof
一つの点を取ると $b(0)\ge0$ となる。二点の二次形式が全ての複素係数で非負の実数なので、実数と純虚数の係数を代入して非対角成分が共役になることを確かめられる。行列式の非負性は、$b(0)>0$ なら一方の係数を固定し他方で二次形式を最小にして得られる。$b(0)=0$ なら同じ二次形式の交差項は正負いずれにもできるため $b(x)=0$ である。

Bochner 表示を取り、$b(x)-b(y)=\int(\chi(x)-\chi(y))\,d\mu$ に Cauchy–Schwarz を使う。絶対値の二乗の積分は $2\int(1-\Re\chi(x-y))\,d\mu$ となるので所望の評価が従う。和と極限は有限の二次形式で判定できる。積については双対側の二つの正測度の畳み込みの逆変換になることを Fubini で確かめればよい。
:::

:::example[絶対値の上限を通過しても正定値とは限らない]
円周で $b(t)=2+2\cos t-\cos(2t)$ とする。$z=\cos t$ により $b=3+2z-2z^2$ である。例えば $z=1/2$ で $b=7/2>b(0)=3$ なので、二点の判定で正定値でない。二点でこのような違反の出ない例を作るなら、$c(t)=1-\eps(1-\cos t)^2$、$0<\eps\le1/2$ と置く。$-1\le c(t)\le1=c(0)$ なので二点の必要条件は満たすが、$\pm2$ の係数は $-\eps/4$ で負である。Bochner 表示の一意性により、この係数は正測度の質量にはできない。
:::

### 自己相関のスペクトル

確率分布のほかに、信号を平行移動して元との重なりを測る場合にも正定値関数が現れる。$L^2$ の入力なら、積分の値を全ての移動について定められる。

:::proposition[自己相関と二乗スペクトル]{#prop-four-lca-autocorrelation}
$f\in L^2(G)$ に対し $b(x)=\langle f,\tau_x f\rangle$ とする。$b$ は連続な正定値関数であり、その Bochner 測度は $|\widehat f|^2d\chi$ である。逆に $d\chi$ に関して密度 $w\in L^1$ を持つ有限正測度の逆変換は、ある $L^2$ 関数の自己相関になる。
:::

:::proof
内積は第 1 引数線形とする。移動の変換は $\overline{\chi(x)}\widehat f(\chi)$ なので、Plancherel により
$$
b(x)=\int_{\widehat G}|\widehat f(\chi)|^2\chi(x)\,d\chi .
$$
密度は可積分であり、:ref[thm-four-bochner] により結論が従う。逆に $\sqrt w\in L^2$ なので Plancherel の全射性で $\widehat f=\sqrt w$ とできる。位相を掛けても二乗スペクトルは同じであるため、この $f$ は一般には一意でない。
:::

:::example[箱の自己相関]
$\R$ 上で $f=1_{[0,2]}$ を取る。$[0,2]$ と $[x,x+2]$ の重なりの長さは $(2-|x|)_+$ なので、これが自己相関である。変換は $\widehat f(\xi)=e^{-2\pi i\xi}\sin(2\pi\xi)/(\pi\xi)$ であり、Bochner 密度は $\sin^2(2\pi\xi)/(\pi\xi)^2$ である。密度の全積分は $2$ となる。移動した箱でも自己相関は同じであり、この計算は密度だけから箱の原点を決められないことを示す。
:::

:::example[離散スペクトルの相関係数]
円周の正規化測度で $f(t)=2e^{it}+ie^{-2it}$ とする。$\tau_x f(t)=2e^{-ix}e^{it}+ie^{2ix}e^{-2it}$ なので、直交性により $b(x)=4e^{ix}+e^{-2ix}$ である。スペクトルの二つの質量は $4,1$ である。例えば $x=\pi$ で $b(\pi)=-3$ となるが、$|b(\pi)|\le b(0)=5$ は保たれる。
:::

### Euclid 空間での連続性定理

Bochner の表示では、点ごとに指標積分が収束するときに測度も収束するかが問題になる。連続な極限であれば原点近くにおける値を使って、周波数が無限遠へ逃げないことを示せる。ここでは $\R^n$ でこの議論を行う。

:::proposition[確率測度の指標積分と弱収束]{#prop-four-lca-levy-continuity}
$\R^n$ 上の確率測度 $\mu_j$ に対し $b_j(x)=\int e^{2\pi i x\cdot\xi}d\mu_j(\xi)$ とする。各 $x$ で $b_j(x)\to b(x)$ かつ $b$ が $0$ で連続なら、$b$ の Bochner 測度 $\mu$ は確率測度であり、任意の有界連続関数 $h$ に対し $\int h\,d\mu_j\to\int h\,d\mu$ である。
:::

:::proof
有限の二次形式を極限に移せるので $b$ は正定値性を満たし、$b(0)=1$ である。Bochner の定理を使う前にも、その証明で作った $e_x$ の内積に Cauchy–Schwarz を使えば $|b(x)-b(y)|^2\le2(1-\Re b(x-y))$ を得る。従って $b$ は全点で連続になり、:ref[thm-four-bochner] で確率測度 $\mu$ を得る。

$g_\eps(x)=\eps^{-n}e^{-\pi |x|^2/\eps^2}$ と置く。Gauss 変換と Fubini より
$$
\int(1-e^{-\pi\eps^2|\xi|^2})\,d\mu_j(\xi)
=\int g_\eps(x)(1-\Re b_j(x))\,dx .
$$
固定 $\eps$ で右辺は $b$ を代入した値へ収束する。さらに $\eps\downarrow0$ でその値は $0$ に近づく。実際、$x=\eps y$ と変換して $2e^{-\pi|y|^2}$ を優関数に使える。$|\xi|\ge1/\eps$ では左辺の被積分関数が $1-e^{-\pi}$ 以上なので、十分小さな $\eps$ を選べば全ての十分大きな $j$ で尾質量を任意に小さくできる。有限個の残りの確率測度については半径をさらに広げる。

任意の $\phi\in\mathcal S(\R^n)$ に対し、反転と Fubini により $\int\phi\,d\mu_j=\int\widehat\phi(x)b_j(x)\,dx$ となる。優収束で $\mu$ にも同じ極限が得られる。$C_c$ の関数は台を少し広げた $C_c^\infty$ 関数で一様近似できるため、収束は $C_c$ に広がる。最後に有界連続関数 $h$ を大きな球の外で切る。上で得た一様な尾評価と $\mu$ の有限性により切断誤差を全ての $j$ で小さくできる。
:::

:::example[原点の連続性が失われる例]
$\mu_j$ を $[-j,j]$ の一様確率測度とする。$b_j(x)=\sin(2\pi j x)/(2\pi j x)$ なので極限は $0$ で $1$、それ以外で $0$ である。これは $0$ で不連続である。実際、$j>R$ なら $[-R,R]$ に残る質量は $R/j$ であり、任意の固定区間から質量が逃げる。一方、$\nu_j=(\delta_{-1/j}+\delta_{1/j})/2$ なら $b_j(x)=\cos(2\pi x/j)$ は $1$ に収束し、極限測度は $\delta_0$ である。
:::

## Pontryagin の双対定理

$x\in G$ を固定すると、指標の値を $x$ で評価する操作は双対群の指標になる。これで $G$ を二重双対に写せる。全射性は、二重双対側の Plancherel 測度に余分な領域が残らないことから示せる。

:::theorem[Pontryagin の双対定理]{#thm-four-pontryagin}
$\iota(x)(\chi)=\chi(x)$ は $G\to\widehat{\widehat G}$ の位相群同型である。
:::

:::proof
$\iota$ は準同型で、連続性はコンパクトな指標集合の同程度連続性による。単射は :ref[lem-four-lca-cstar] の後で示した指標の分離性による。

まず $\iota$ が固有であり閉像を持つことを示す。$F=\widehat g\,\overline{\widehat h}$（$g,h\in C_c(G)$）は $L^1(\widehat G)$ に属し、その逆積分を評価点 $\iota(x)$ で取ると $(g*h^*)(x)$ になる。これは :ref[thm-four-lca-plancherel] の反転による。後者は連続でコンパクト台に属する。こうした積の線形包は $L^1(\widehat G)$ で稠密である。実際 $L^1$ 関数を絶対値の平方根と位相付きの平方根の積に分け、各 $L^2$ 因子を Plancherel と稠密性で近似すれば Hölder で積も近づく。任意の二重双対の元 $\eta$ に対して、ある $F$ の逆積分 $\int F(\chi)\eta(\chi)d\chi$ が非零である。これらが全てゼロなら稠密性により $L^1$ の全てと $\eta$ の積分がゼロになり、$|\eta|=1$ に矛盾する。

従って $\eta$ の近傍で一つの逆積分が正の絶対値下限を持つようにできる。その近傍の $\iota$ による逆像は、対応する $g*h^*$ の絶対値が下限以上となるコンパクト集合に含まれる。二重双対の任意のコンパクト集合を有限個のこのような近傍で覆えば逆像はコンパクトである。さらに、$\iota(x_\alpha)\to\eta$ なら残りの $x_\alpha$ は一つのコンパクト集合に入るため集積部分ネットを持ち、$\eta$ は像に入る。$\eta=\iota(x)$ の場合は単射性から全ての集積点が $x$ に限られ、$x_\alpha\to x$ である。これで閉像と逆の連続性が従う。

$\rho=\iota_*dx$ は固有性によって二重双対上の Radon 測度である。$F\in L^1\cap L^2(\widehat G)$ に対して

$$
\int_{\widehat{\widehat G}}|\widehat F(\eta)|^2\,d\rho(\eta)
=\int_G\left|\int_{\widehat G}F(\chi)\overline{\chi(x)}\,d\chi\right|^2dx
=\|F\|_2^2 .
$$

最後は $G$ の Plancherel の逆写像である。積分式と逆写像の一致は $C_c(G)$ との内積を取り Fubini で示せる。反射は Haar 測度を保存するので符号はノルムに影響しない。

二重双対の Plancherel Haar 測度を $\rho_0$ とする。どちらの測度でも上のノルム等式が成り立つ。偏極し、$F\in C_c(\widehat G)$ と $h*F$（$h\in L^1(\widehat G)$）の組に適用すると $\int\widehat h|\widehat F|^2d\rho=\int\widehat h|\widehat F|^2d\rho_0$ である。:ref[lem-four-lca-cstar] を $\widehat G$ に適用すれば $\widehat h$ は $C_0$ に稠密なので、重み付き有限測度は一致する。各点で $\widehat F\ne0$ となる $F$ を選べるため、局所的に割って $\rho=\rho_0$ を得る。Haar 測度は空でない開集合に正の質量を持つが、$\rho$ は閉像 $\iota(G)$ に台を持つ。従って閉像の補集合は空である。
:::

:::example[整数群に二回適用する]
$\Z$ を二回双対化すると $\widehat{\mathbb T}=\Z$ であり、$\iota(k)$ は $t\mapsto e^{ikt}$ である。符号が変わる $\mathcal F^2 f(x)=f(-x)$ は負の符号で二回変換した結果である。評価写像の式そのものに負の符号は入らない。
:::

### 閉部分群と消滅化群

指標を閉部分群上で $1$ になるものに制限すると、商群の指標が得られる。この対応は、周波数の集合から元の空間における部分群を復元する手掛かりになる。

:::proposition[消滅化群と商]{#prop-four-lca-annihilator}
閉部分群 $H\le G$ に対し $H^\perp=\{\chi\in\widehat G:\chi(h)=1\ \text{for all }h\in H\}$ は閉部分群であり、商写像 $\pi:G\to G/H$ による合成は $\widehat{G/H}$ から $H^\perp$ への位相群同型である。また :ref[thm-four-pontryagin] の同定の下で $(H^\perp)^\perp=H$ である。
:::

:::proof
$H^\perp$ は各 $h$ の評価連続写像の核の共通部分である。商群には商位相を入れる。$\pi$ は開写像である。実際、開集合 $U$ の像の逆像は $U+H$ で開いている。$H$ が閉であることから商群は Hausdorff であり、コンパクトな単位元近傍の像が近傍を含むので局所コンパクトである。$H$ 上で自明な指標は商位相により連続な指標へ一意に降りる。逆に任意の商の指標は $G$ 上の指標を与える。

$G$ のコンパクト集合の像はコンパクトなので、合成は連続である。商のコンパクト集合 $K$ は $G$ のあるコンパクト集合の像に含まれる。各点の代表元の周囲にコンパクト閉包を持つ開近傍を取り、その像による $K$ の開被覆を有限に減らす。近傍の閉包を有限個合わせれば所望のコンパクト集合である。これにより逆向きのコンパクト上一様収束も移る。最後に $x\notin H$ なら $x+H\ne0$ なので、$G/H$ の分離性により値が $1$ と異なる指標がある。その引き戻しは $H^\perp$ に入り、$x$ を検出する。$H$ の各点は定義により全指標で $1$ になるので等式が従う。
:::

:::example[整数格子の消滅化群]
$G=\R$、$H=a\Z$（$a>0$）なら $\chi_\xi(ak)=e^{2\pi i \xi ak}$ が全ての $k$ で $1$ となる条件は $a\xi\in\Z$ である。従って $H^\perp=a^{-1}\Z$ である。これは :eqref[eq-four-lattice-poisson] の一次元の双対格子と同じである。$a=2$ なら周波数は半整数になる。
:::

### 閉部分群への指標の制限

商群の指標は引き戻しによって得られた。反対に、$G$ 上の指標を部分群に制限すると何が得られるかを調べる。部分群の指標を全群へ延ばせることは、有限 Abel 群では生成元を一つ加える操作で示せた。一般の位相群では延長の連続性も必要である。

:::theorem[閉部分群の指標の延長]{#thm-four-lca-character-extension}
$H$ を $G$ の閉部分群とする。制限写像 $\widehat G\to\widehat H$ は連続な開全射で、核は $H^\perp$ である。従って $\widehat G/H^\perp\simeq\widehat H$ が位相群同型として成り立つ。
:::

:::proof
制限はコンパクト集合上一様収束の定義により連続であり、核は定義から $H^\perp$ である。$D=\widehat G/H^\perp$ と置くと、:ref[prop-four-lca-annihilator] を $\widehat G$ に適用して $\widehat D\simeq(H^\perp)^\perp\simeq H$ を得る。この対応は、$h\in H$ を剰余類 $[\chi]$ 上の指標 $[\chi]\mapsto\chi(h)$ に送る写像である。Pontryagin の定理により $D$ は $\widehat{\widehat D}$ と位相群同型である。双対化は位相群同型を逆向きの位相群同型へ移す。実際、元の写像と逆写像がそれぞれコンパクト集合をコンパクト集合へ送るので、指標の引き戻しとその逆は連続である。そこで $\widehat D\simeq H$ を双対化し $D$ を $\widehat H$ と同定すると、$[\chi]$ に対応する値は $h\mapsto\chi(h)$ となる。従ってこの同型と商写像の合成が制限である。商写像は開であるため制限も開である。
:::

:::example[格子上の指標の延長は一意でない]
$G=\R^2$、$H=2\Z\times\{0\}$ とし、$H$ の指標を $(2k,0)\mapsto e^{ik\pi/3}$ とする。$\chi_\xi(x)=e^{2\pi i x\cdot\xi}$ にこの値を指定すると $4\pi\xi_1\equiv\pi/3\bmod2\pi$ となる。従って全ての延長は $\xi_1=1/12+\ell/2$（$\ell\in\Z$）、$\xi_2\in\R$ である。延長の差は $\frac12\Z\times\R=H^\perp$ であり、これが自由度の全体を表す。
:::

:::example[有限指標を円周の整数モードに延長する]
$\mathbb T$ の部分群 $H=\{2\pi j/6:j=0,\ldots,5\}$ を取る。$H$ 上の指標 $2\pi j/6\mapsto e^{2\pi i 5j/6}$ は $e^{ikt}$ で $k\equiv5\bmod6$ のもの全てに延長される。例えば $k=5$ と $k=-1$ は $H$ 上で同じ値を持つ。消滅化群は $6\Z$ であり、$\widehat H\simeq\Z/6\Z$ となる。
:::

### 商積分と部分群への周期化

閉部分群に沿って積分すると、剰余類の中で代表元を変えても同じ値になる。この関数を商群上で積分して元の積分へ戻すには、商の Haar 測度の正規化が必要である。ここでは $dx$ と $dh$ を先に固定する。

:::proposition[Weil の商積分公式]{#prop-four-lca-quotient-integral}
$H\le G$ を閉部分群とし、$\pi:G\to G/H$ とする。$f\in C_c(G)$ に対し
$$
P_H f(\pi(x))=\int_H f(x+h)\,dh
$$
は $C_c(G/H)$ に属する。商群には一意な Haar 測度 $d\dot x$ を選べて
$$
\int_G f(x)\,dx=\int_{G/H}P_H f(\dot x)\,d\dot x
\label{eq-four-lca-weil-integral}
$$
となる。
:::

:::proof
代表元の独立性は $dh$ の平行移動不変性による。$x$ をコンパクト近傍 $K$ に制限すると、積分に寄与する $h$ は $H\cap(\operatorname{supp}f-K)$ というコンパクト集合内にある。その上で $f$ の平行移動は一様連続なので、$P_Hf\circ\pi$ は連続である。商位相で $P_Hf$ が連続になり、台は $\pi(\operatorname{supp}f)$ に含まれる。

各 $\dot x$ では $P_Hu(\dot x)>0$ となる非負 $u\in C_c(G)$ を選べる。実際、代表元 $x$ の近傍で $u$ を正にすると、$H$ の単位元近傍が正の Haar 測度を持つので積分は正となる。$P_Hu$ の連続性と有限被覆により、任意のコンパクト集合 $C\subset G/H$ に対して $P_Hu>0$ が $C$ の近傍で成り立つよう $u$ を選べる。$\psi\in C_c(G/H)$、$\operatorname{supp}\psi\subset C$ に対し
$$
f(x)=u(x)\frac{\psi(\pi(x))}{P_Hu(\pi(x))}
$$
と置く。分母がゼロの場所では分数を $0$ と定める。$\psi$ の台の近傍では分母が正であり、台の外では分子が局所的にゼロなので、この定義は連続で $f\in C_c(G)$ となる。$P_H f=\psi$ で、$\psi\ge0$ なら $f\ge0$ とできる。

次に $P_Hf=0$ なら $\int f\,dx=0$ を示す。$P_Hu>0$ を $\pi(\operatorname{supp}f)$ の近傍で満たす $u$ を取り、このコンパクト集合上で $1$ になる $\rho\in C_c(G/H)$ を $\{P_Hu>0\}$ 内に台を持つように選ぶ。$v(x)=u(x)(\rho/P_Hu)(\pi(x))$ は $P_Hv=1$ を $\pi(\operatorname{supp}f)$ 上で満たす。そこで
$$
\int_G f(x)P_Hv(\pi(x))\,dx
=\int_G v(y)\left(\int_H f(y-h)\,dh\right)\,dy=0
$$
となる。積分交換は $H$ の変数が二つのコンパクト台の差に制限されるため正当であり、内側の積分では変数を反射して $P_Hf(\pi(y))$ を得る。左辺は $\int f\,dx$ に等しい。

従って $\psi\mapsto\int f\,dx$ は $P_Hf=\psi$ の持ち上げ方によらず、正値の持ち上げを選べるため正汎関数である。非負非零 $\psi$ の持ち上げは非零なので、$dx$ に関する積分も正である。Riesz の表現定理で非零 Radon 測度を得る。商群の移動を $G$ の移動に持ち上げれば、この正汎関数は平行移動で不変である。従って測度は Haar 測度であり、$P_H$ の全射性により正規化も一意である。
:::

:::corollary[消滅化群上の Poisson 公式]{#cor-four-lca-poisson}
:ref[prop-four-lca-quotient-integral] の規約で $H^\perp\simeq\widehat{G/H}$ に $d\dot x$ の双対測度 $d\chi_\perp$ を入れる。$f\in C_c(G)$ かつ $\widehat f|_{H^\perp}\in L^1(H^\perp,d\chi_\perp)$ なら、全ての $x\in G$ で
$$
\int_H f(x+h)\,dh=\int_{H^\perp}\widehat f(\chi)\chi(x)\,d\chi_\perp
$$
が成り立つ。
:::

:::proof
$\chi\in H^\perp$ は $\overline{\chi(x)}$ が $\pi(x)$ だけに依存するため、:eqref[eq-four-lca-weil-integral] を $f\overline\chi$ に適用すれば $\widehat{P_H f}(\dot\chi)=\widehat f(\chi)$ となる。:ref[thm-four-lca-plancherel] の反転公式を商群上で使うと a.e. の等式を得る。左辺は既に連続であり、右辺は可積分の変換を逆積分した連続関数である。空でない開集合の Haar 測度は正なので、二つの連続関数は全点で一致する。
:::

:::example[格子間隔で変わる商の測度]
$G=\R$ で $dx$ は Lebesgue 測度、$H=a\Z$ では数え上げ測度とする。$x$ を長さ $a$ の基本区間に取り、$\int_0^a\sum_k f(x+ka)dx=\int_\R f$ と変数変換すれば、商測度は $[0,a)$ 上の $dx$ に対応する。商の全質量は $a$ なので、双対側では各周波数 $m/a$ に $1/a$ の質量を与える。例えば $a=2$、$f(x)=(1-|x|)_+$ を $x=0$ で周期化すると値は $1$ である。$\widehat f(\xi)=(\sin(\pi\xi)/(\pi\xi))^2$ は箱の畳み込みから得られるので、公式は
$$
1=\frac12+\frac12\sum_{m\text{ odd}}\frac{4}{\pi^2m^2}
$$
となる。正負の奇整数を含めて和を取り、片側だけの和なら $\pi^2/8$ になる。
:::

:::example[有限な輪の点を束ねる]
$G=\Z/12\Z$ では正規化測度、$H=\{0,4,8\}$ では数え上げ測度を用いる。すると $P_Hf$ は四つの剰余類の各々で三点の和を取るため、:eqref[eq-four-lca-weil-integral] を満たす商の一点の重みは $1/12$ である。商の全質量は $1/3$、従って $H^\perp$ の各点の重みは $3$ である。$f=\delta_0+2\delta_5$ なら $P_H f(0)=1$、$P_H f(1)=2$ である。$H^\perp$ は周波数添字 $k=0,3,6,9$ からなる。反転で $3\sum_k \widehat f(k)e^{2\pi i kx/12}$ を計算すると、$3/12$ と四つの根の和が相殺してちょうど $1,2$ が復元される。
:::

### Bohr コンパクト化と概周期関数

双対群に離散位相を入れてから再び双対を取ると、コンパクト群が得られる。その群の連続関数を引き戻すと、指標の有限和で一様近似できる関数がちょうど得られる。

:::proposition[Bohr コンパクト化と移動族]{#prop-four-bohr-compactification}
$\widehat G$ に離散位相を入れた群を $\widehat G_d$ とし、$bG=\widehat{\widehat G_d}$ と置く。評価写像 $\beta:G\to bG$ は連続単射で稠密像を持つ。任意のコンパクト Hausdorff Abel 群 $K$ への連続準同型 $\alpha:G\to K$ は、一意な連続準同型 $\widetilde\alpha:bG\to K$ を通して $\alpha=\widetilde\alpha\circ\beta$ と表せる。この $bG$ を :term[Bohr コンパクト化]{en="Bohr compactification" yomi="ぼーあこんぱくとか"} という。

有界連続関数 $f:G\to\C$ について、次の三条件は同値である。

1. $f=F\circ\beta$ を満たす $F\in C(bG)$ がある。
2. $f$ は $G$ の指標の有限線形結合の一様極限である。
3. $\{\tau_x f:x\in G\}$ の一様ノルムでの閉包はコンパクトである。

これらを満たす $f$ を :term[概周期関数]{en="almost periodic function" yomi="がいしゅうきかんすう"} という。
:::

:::proof
コンパクト性は :ref[lem-four-lca-discrete-compact] による。$\widehat G_d$ のコンパクト集合は有限なので、各指標の連続性から $\beta$ も連続である。指標は $G$ の点を分離するため単射である。$bG$ の指標は :ref[thm-four-pontryagin] により $\widehat G_d$ の元で表される。その指標が $\beta(G)$ の閉包上で $1$ なら、対応する $\chi\in\widehat G_d$ は $G$ の全点で $1$ なので自明である。:ref[prop-four-lca-annihilator] の二重消滅化群の等式から閉包は $bG$ 全体になる。

$\alpha$ の双対写像は $\widehat K$ から $\widehat G_d$ への $\eta\mapsto\eta\circ\alpha$ である。$\widehat K$ は離散なので連続である。この写像を再び双対化し、$\widehat{\widehat K}$ を $K$ と同定すれば $\widetilde\alpha$ を得る。指標上で評価すれば合成の式が分かり、稠密性により $\widetilde\alpha$ は一意である。

指標の有限線形結合は $bG$ 上で共役に閉じ、定数を含み、点を分離する。従って Stone–Weierstrass により $C(bG)$ に一様稠密である。さらに $\beta(G)$ が稠密なので $F\mapsto F\circ\beta$ は一様ノルムの等長写像であり、その像は閉じている。これで条件 1 と 2 が同値である。有限線形結合の移動は、各係数に単位円の位相を掛けるだけなのでコンパクト集合に含まれる。移動は等長だから、一様近似を使えば $f$ の移動全体も全有界となる。有界連続関数の一様ノルムでの完備性により条件 3 が従う。

逆に条件 3 を仮定する。任意の $\eps>0$ に対して移動族の有限 $\eps$-網を取り、その各関数の単位元での連続性を使う。十分小さい $h$ に対し、移動族の任意の $g$ は $|g(h)-g(0)|\le3\eps$ を満たす。網の近い関数との差が二箇所で $\eps$ 以下だからである。$g=\tau_{-y} f$ とすれば $\sup_y|f(y+h)-f(y)|\to0$ を得る。

$d(x,y)=\|\tau_x f-\tau_y f\|_\infty$ は平行移動不変な擬距離である。$N=\{x:d(x,0)=0\}$ で割ると距離となり、条件 3 により全有界である。距離完備化 $K_f$ はコンパクトである。加法は $d(x+y,x'+y')\le d(x,x')+d(y,y')$、反転は等長であるため、両方とも完備化へ延びて $K_f$ は Abel 群となる。自然な準同型 $\alpha_f:G\to K_f$ の連続性は、移動の一様連続性から従う。移動した二関数を $x+y$ で評価すると $|f(x)-f(y)|\le d(x,y)$ を得るので、$f$ は $K_f$ 上の連続関数に延びる。$\alpha_f$ を $\beta$ を通して表せば条件 1 が従う。
:::

:::example[実軸に沿って二つの位相を進める]
$f(x)=e^{2\pi i x}+e^{2\pi i\sqrt2\,x}$ は二つの指標の和である。移動すると二つの係数がそれぞれ単位円上で動くので、移動族は一様ノルムで相対コンパクトである。一方、周期 $T$ があれば二つの異なる周波数の一次独立性から $T,\sqrt2 T$ はともに整数になる。$T=0$ しかないため、概周期性は単一の正の周期を持つことより広い性質である。
:::

## Wiener の Tauber 型定理

平行移動は Fourier 側で指標を掛ける操作である。変換のある周波数がゼロなら、この操作をいくら繰り返してもそこで値を作れない。逆に零点が一つもなければ、局所的な逆元を組み合わせて任意の入力を近似できる。

:::lemma[周波数の局所切断]{#lem-four-lca-local-cutoff}
$\widehat G$ の任意の点 $\chi_0$ とその開近傍 $U$ に対し、$u\in L^1(G)$ で $\widehat u$ が $0\le\widehat u\le1$、$U$ 内のコンパクト台を持ち、$\chi_0$ の近くで $1$ となるものが存在する。
:::

:::proof
単位元の場合、対称なコンパクト近傍 $V$ を $V V V\subset U$ となるように取る。$d\chi$ に関する畳み込みで $b=1_{VV}*1_V/|V|$ と置く。$b$ は $V$ 上で $1$ であり、台は $VVV$ に含まれ、値は $[0,1]$ にある。$b$ は二つのコンパクト台 $L^2$ 関数の畳み込みである。それぞれの逆変換の積は $L^1(G)$ に属し、両因子が有界でもあるので積は $L^2$ にも属する。また双方が周波数側で $L^1$ にも属するため、逆積分の積に Fubini を使えば $b$ の逆積分になる。Plancherel と $L^1\cap L^2$ 上での一致により、この積の変換は $b$ である。一般の中心は変調で移す。
:::

:::theorem[Wiener の Tauber 型定理]{#thm-four-wiener-tauberian}
$f\in L^1(G)$ の平行移動の線形包が $L^1(G)$ で稠密であることと、$\widehat f(\chi)\ne0$ が全ての $\chi\in\widehat G$ で成り立つことは同値である。
:::

:::proof
零点があれば、そこでの変換値を取る非零連続汎関数が全移動を消す。逆のため、$I$ を $f*g$（$g\in L^1$）の閉包とする。この閉包は平行移動の線形包の閉包に等しい。畳み込みを移動した $f$ の積分として、コンパクト上の移動連続性と尾の評価で有限和へ近似できる。一方、$\tau_x f$ は $\tau_x$ を施した近似単位と $f$ の畳み込みで近似できる。

$\chi_0$ を固定し $a=\widehat f(\chi_0)\ne0$ とする。$|\widehat f-a|<|a|/2$ の近傍内で :ref[lem-four-lca-local-cutoff] の $u$ を選ぶ。単位化した代数で $v=f*u+a(1-u)$ を作ると、各指標値は $a+\widehat u(\widehat f-a)$ でありゼロでない。単位化の追加指標でも値は $a$ である。可換 Banach 代数において、全指標値が非零の元は可逆である。実際、非可逆ならその生成するイデアルが極大イデアルに含まれ、商の指標で値がゼロになる（:ref[thm-fa-maximal-characters]）。よって $v$ は可逆である。$\widehat u=1$ の開集合内に $\widehat g$ の台が含まれるなら $(1-u)*g=0$ である。変換が恒等的にゼロになるためである。従って $v*g=f*u*g$ であり、逆元を掛ければ $g\in I$ となる。

$\widehat g$ がコンパクト台なら、台をこのような開集合の有限個で覆う。各開集合内に台を持つ切断 $\widehat{w_j}$ を取り、$\widehat{w_j}=1$ となる小近傍で台を覆う。$g_j=g*w_j*\prod_{l<j}(1-w_l)$ と置けば各 $g_j\in I$ である。変換の側では望遠和により $\sum_j\widehat g_j=\widehat g(1-\prod_j(1-\widehat w_j))=\widehat g$ なので $g\in I$ である。

変換がコンパクト台の $L^1$ 関数は稠密である。$g\in L^1$ を二つの $L^2$ 因子の積に分け、各因子の変換をコンパクト集合へ切る。この近似は Plancherel によって $L^2$ で成り立ち、積は Hölder で $L^1$ 収束する。近似した積の変換は二つのコンパクト台の変換の畳み込みであり、コンパクト台を持つ。これは :ref[lem-four-lca-local-cutoff] と同じ Fubini の計算による。従って $I=L^1(G)$ である。
:::

:::example[Gauss 関数の移動と箱の移動]
$\R$ で $e^{-\pi x^2}$ の変換はどの周波数でも正なので、その移動の線形包は $L^1$ で稠密である。$1_{[-1/2,1/2]}$ は $\widehat f(1)=0$ である。従ってその移動の線形包の閉包には、例えば $e^{-\pi x^2}$ は入らない。
:::

### Wiener の補題との対応

:::corollary[絶対収束級数の逆数]{#cor-four-wiener-inverse-series}
$F(t)=\sum_{k\in\Z}a_ke^{-ikt}$、$\sum_k|a_k|<\infty$ とし、全点で $F(t)\ne0$ とする。ここでは $F$ を $\Z$ 上の $a$ の Fourier 変換として書いているので、円周の規約での Fourier 係数は $\widehat F(k)=a_{-k}$ である。このとき $1/F$ も絶対収束 Fourier 級数を持つ。
:::

:::proof
$a$ を $\Z$ 上の $L^1$ 関数とみなす。:ref[thm-four-wiener-tauberian] により $a*g$ の集合は $\ell^1$ に稠密である。従って $g$ を $\|\delta_0-a*g\|_1<1$ となるよう選べる。Neumann 級数 $\sum_{j\ge0}(\delta_0-a*g)^{*j}$ は絶対収束し、$(a*g)$ の逆元を与える。これに $g$ を畳み込み、$a$ の逆元 $h\in\ell^1$ を得る。Fourier 変換で $F(t)\widehat h(t)=1$ なので $\widehat h$ の絶対収束級数が求める表示である。
:::

:::example[対称な三項式の逆数]
$a_0=3$、$a_{\pm1}=1$、ほかをゼロとすると $F(t)=3+2\cos t$ である。$F\ge1$ であり、直接にも
$$
\frac1{F(t)}=\frac13\sum_{j=0}^\infty\left(-\frac{e^{it}+e^{-it}}3\right)^j
$$
と展開できる。右辺の第 $j$ 項の係数の絶対和は高々 $(1/3)(2/3)^j$ なので、全体は絶対可和である。一般の Wiener の補題では、単一の定数項が残りの係数の絶対和より大きいという仮定を課さずに逆数の絶対可和性を得る。
:::

### 複数の入力で生成する閉イデアル

実際の畳み込みでは、複数の入力のうち一つだけでは零点を消せない場合がある。各周波数で少なくとも一つの変換が値を持てば、その近傍では対応する入力を使って局所逆元を作れる。

:::corollary[平行移動族と共通零点]{#cor-four-lca-wiener-family}
空でない族 $\{f_\alpha\}\subset L^1(G)$ の全ての平行移動の線形包が $L^1(G)$ で稠密であるための必要十分条件は $\bigcap_\alpha\{\widehat f_\alpha=0\}=\varnothing$ である。
:::

:::proof
共通零点があれば、そこでの変換評価が全ての平行移動を消す。逆に共通零点が空とする。有限和 $\sum_{\alpha\in F}f_\alpha*g_\alpha$（$g_\alpha\in L^1$、$F$ 有限）の閉包 $I$ は、:ref[thm-four-wiener-tauberian] の証明と同じ移動近似で全ての平行移動の線形包の閉包に一致する。各 $\chi$ では $\widehat f_\alpha(\chi)\ne0$ になる $\alpha$ を選び、その変換を使って定理の証明と同じ局所切断と逆元を作る。従ってその小近傍内に変換の台を持つ全ての $g$ は $f_\alpha*L^1$ に入る。コンパクトな周波数台については有限被覆と定理の証明の切断の積を使えば $g\in I$ である。こうした $g$ の稠密性により $I=L^1(G)$ となる。
:::

:::example[二種類の差分から一点を作る]
$\Z$ で $f_1=\delta_0+\delta_1$、$f_2=\delta_0-\delta_1$ とする。それぞれの変換は $1+e^{-it}$、$1-e^{-it}$ であり、零点は $t=\pi$ と $t=0$ で異なる。各入力だけでは稠密性が得られないが、$(f_1+f_2)/2=\delta_0$ なので、両者を使えば任意の有限台の列を有限個の移動で正確に表せる。この場合には閉包を取る前から有限台の列が得られる。
:::

:::example[幅が無理数比の二つの箱]
$\R$ で $f_a=1_{[-a/2,a/2]}$ とすると $\widehat f_a(\xi)=\sin(\pi a\xi)/(\pi\xi)$ で、零点は $(1/a)\Z\setminus\{0\}$ である。$a=1$ と $a=\sqrt2$ に共通零点があると整数 $k\ne0$ が $m/\sqrt2$ に等しくなり、$\sqrt2=m/k$ となって矛盾する。従ってこの二種類の箱の移動で $L^1(\R)$ の任意の関数を近似できる。幅を $1$ と $2$ に選ぶと $\xi=1$ が共通零点なので、同じ結論は得られない。
:::

### 一つの平均の極限から別の平均へ

Tauber 型定理は、幅や形を固定した一種類の平均の極限から、別の平均の極限を導くためにも使える。ここでは無限遠へ向かう変数を持つ $\R$ 上で、すべての点で定義できる積分を扱う。

:::proposition[有界関数に対する Wiener の極限移送]{#prop-four-lca-wiener-limit}
$B\in L^\infty(\R)$、$f\in L^1(\R)$ とし、$\widehat f$ に零点がないとする。ある $c\in\C$ に対し
$$
(f*B)(x)\longrightarrow c\int_\R f\qquad(x\to+\infty)
$$
なら、任意の $g\in L^1(\R)$ に対し $(g*B)(x)\to c\int_\R g$ である。
:::

:::proof
$B-c$ に置き換えて $c=0$ を示せばよい。各 $x$ で絶対積分が可能で、$|(g*B)(x)|\le\|g\|_1\|B\|_\infty$ となる。$f$ の移動について $(\tau_a f)*B=\tau_a(f*B)$ であるから、$x\to+\infty$ の極限はゼロである。有限線形結合についてもゼロである。:ref[thm-four-wiener-tauberian] によって任意の $g$ をその有限線形結合 $g_j$ で $L^1$ 近似する。すると $\sup_x|((g-g_j)*B)(x)|\le\|g-g_j\|_1\|B\|_\infty$ となる。まず $j$ を選んでこの上界を小さくし、次に $x$ を大きくして $(g_j*B)(x)$ を小さくすれば結論が従う。
:::

:::example[Gauss 平均から片側平均へ]
$f(x)=e^{-\pi x^2}$ とし、$(f*B)(x)\to c$ を仮定する。この核は積分 $1$ で、変換は全点で正である。従って $g(t)=1_{[0,1]}(t)$ を使えば $\int_{x-1}^x B(y)dy\to c$ となる。$g(t)=2\,1_{[-1/2,0]}(t)$ を使えば $2\int_x^{x+1/2}B(y)dy\to c$ である。二つの平均は区間の向きも長さも違うが、積分が共に $1$ なので同じ極限を持つ。
:::

:::example[零点があると平均を移せない]
$f=1_{[-1/2,1/2]}$、$B(x)=\cos(2\pi x)$ とする。$f*B$ は周期一つ分の積分なので全点で $0$ である。$g=1_{[-1/4,1/4]}$ については、$\cos(2\pi(x-y))$ の加法定理と奇関数の積分を用いると $(g*B)(x)=\pi^{-1}\cos(2\pi x)$ である。従って $x\to+\infty$ の極限を持たない。$f$ の変換が $\pm1$ でゼロになるため、元の平均には $B$ の振動が現れなかった。
:::

### 絶対収束級数と解析関数の合成

逆数は $z\mapsto1/z$ による合成である。同じ代数を使えば、級数の値域を含む複素平面の開集合上で正則な関数についても絶対収束を保てる。級数を項ごとに Taylor 展開する場合には、値域全体が一つの収束円に入るとは限らないので、Cauchy の積分公式を使う。

:::corollary[Wiener–Lévy の合成公式]{#cor-four-lca-wiener-levy}
$F(t)=\sum_{k\in\Z}a_ke^{-ikt}$、$a\in\ell^1(\Z)$ とする。$F(\mathbb T)$ を含む開集合 $U\subset\C$ 上で $\Phi$ が正則なら、$\Phi\circ F$ は絶対収束 Fourier 級数を持つ。
:::

:::proof
$K=F(\mathbb T)$ はコンパクトである。十分細かい正方格子を用い、$K$ に触れる閉正方形とその隣接正方形を有限個集める。格子幅を $K$ から $U^c$ までの距離より十分小さく取れば、その和集合は $U$ に含まれ、$K$ は内部に入る。$U=\C$ の場合は距離の条件は不要である。和集合の各正方形の正向き境界から互いに逆向きに重なる辺を消し、残りの辺を向き付きの線分の和 $\Gamma$ とする。$\Gamma$ は $K$ に交わらない。各正方形が $U$ に含まれるため $\Gamma$ は $U$ で零ホモローグであり、$K$ 上の点に対する巻き数は $1$ である。実際、格子線外の内部の点ではその点を含む正方形だけが寄与する。格子線上であっても外側の境界から離れており、巻き数は連続なので $1$ である。

$z\in\Gamma$ なら $z-F(t)\ne0$ である。:ref[cor-four-wiener-inverse-series] により $z\delta_0-a$ の畳み込み逆元 $r(z)$ は $\ell^1$ に属する。$z$ を少し動かすと
$$
(z+h)\delta_0-a=(z\delta_0-a)*(\delta_0+h r(z))
$$
であり、$|h|\|r(z)\|_1<1$ では最後の因子の逆元が Neumann 級数で連続に変わる。従って $r(z)$ は $\Gamma$ 上連続で、$\ell^1$ の完備性により各線分上の積分を Riemann 和で定義できる。そこで
$$
d=\frac1{2\pi i}\int_\Gamma \Phi(z)r(z)\,dz\in\ell^1
$$
と置く。各周波数 $t$ での変換評価は連続汎関数である。積分に通すと $r(z)$ は $(z-F(t))^{-1}$ へ移るので、:ref[thm-cpx-cauchy-formula] により $\widehat d(t)=\Phi(F(t))$ となる。
:::

:::example[正の Wiener 級数の平方根]
$F(t)=3+2\cos t$ は値域が $[1,5]$ である。右半平面で正の実数上に正の値を取る平方根を選べば、:ref[cor-four-lca-wiener-levy] により $\sqrt{3+2\cos t}$ は絶対収束 Fourier 級数を持つ。この例では $F=3(1+u)$、$u=(e^{it}+e^{-it})/3$ と書けば $\|u\|_{\ell^1}=2/3$ なので、二項級数を代数内で使ってもよい。具体的に第 1 次までの近似は $\sqrt3(1+\cos t/3)$ であり、第 2 次項は $-(\sqrt3/18)\cos^2 t$ である。$F$ が正であっても任意の中心からの二項級数が代数内で絶対収束するとは限らず、その場合にも合成公式は使える。
:::

## 演習問題

::::exercise{level=basic}
$\Z$ の数え上げ測度で $f=2\delta_0+\delta_1$ の変換を求めよ。ここで $\delta_k$ は点 $k$ の指示関数である。
:::solution
指標は $\chi_t(k)=e^{ikt}$ なので $\widehat f(t)=2+e^{-it}$ である。絶対値は $2-1=1$ 以上なので零点を持たない。
:::
::::

::::exercise{level=basic}
円周の Haar 測度を $dt$ に替える。双対測度の一点当たりの重みを求めよ。
:::solution
元の $dt/(2\pi)$ の $2\pi$ 倍なので、:ref[thm-four-lca-plancherel] のノルム式を保つには逆数を掛ける。重みは $1/(2\pi)$ である。
:::
::::

::::exercise{level=basic}
$b(t)=3+2\cos t$ を Bochner の表示にし、測度の質量を求めよ。
:::solution
$2\cos t=e^{it}+e^{-it}$ なので、$\Z$ の $0$ に質量 $3$、$\pm1$ に $1$ ずつを置く。全質量は $5=b(0)$ である。
:::
::::

::::exercise{level=basic}
有限群 $G=\Z/4\Z$ の二重双対で $x=1$ の評価指標の値を全て書け。
:::solution
第 $k$ 指標は $\chi_k(x)=e^{2\pi ikx/4}$ である。$\iota(1)(\chi_k)=\chi_k(1)$ は $k=0,1,2,3$ の順に $(1,i,-1,-i)$ である。
:::
::::

::::exercise{level=advanced}
正定値関数 $b$ について $b(-x)=\overline{b(x)}$ と $|b(x)|\le b(0)$ を証明せよ。
:::solution
:ref[thm-four-bochner] で $b(x)=\int\chi(x)d\mu$ と書く。$\chi(-x)=\overline{\chi(x)}$ なので一つ目の等式が従う。$|\chi(x)|=1$ より三角不等式で $|b(x)|\le\mu(\widehat G)=b(0)$ である。
:::
::::

::::exercise{level=advanced}
$f=\delta_0+\delta_1\in\ell^1(\Z)$ の移動の線形包が稠密か判定せよ。
:::solution
$\widehat f(t)=1+e^{-it}$ は $t=\pi$ でゼロである。全移動の変換もそこでゼロになる。例えば $\delta_0$ の変換はそこで $1$ なので閉包に入らず、稠密でない。
:::
::::

::::exercise{level=advanced}
$G$ が離散群なら $\widehat G$ がコンパクトであることを示せ。逆も示せ。
:::solution
離散群のコンパクト集合は有限なので、双対の位相は円の $G$ 個の積の点ごとの収束である。準同型の条件は閉条件であり、全て連続なので閉部分集合となってコンパクトである。逆にコンパクト群の双対は離散である。非自明指標の像は円の非自明部分群であり、その中に $1$ からの距離が $1$ より大きい点がある。実際、一つの非自明な値の適当な正負の冪を取れば主偏角の絶対値を $\pi/3$ より大きくできる。従って全群上一様に $1$ へ距離 $1$ 以内の指標は自明なものだけである。$\widehat G$ がコンパクトなら、その双対が離散なので :ref[thm-four-pontryagin] により $G$ も離散である。
:::
::::

::::exercise{level=advanced}
$\ell^1(\Z)$ で $f=2\delta_0+\delta_1$ の畳み込み逆元を具体的に求めよ。
:::solution
$\delta_0$ は単位元で $\delta_1^{*k}=\delta_k$ である。$g=(1/2)\sum_{k=0}^\infty(-1/2)^k\delta_k$ は $\ell^1$ で絶対収束する。$2g+\delta_1*g$ を計算すると正添字の係数は全て相殺し、第 $0$ 係数は $1$ である。従って $f*g=\delta_0$ となる。
:::
::::

::::exercise{level=basic}
$\R$ 上の $\mu=\delta_1+2\delta_3$ について $\widehat\mu(1/4)$ と $\mu*\mu$ を求めよ。
:::solution
$\widehat\mu(1/4)=e^{-\pi i/2}+2e^{-3\pi i/2}=-i+2i=i$ である。二つの点質量の畳み込みは添字の和へ移るので、$\mu*\mu=\delta_2+4\delta_4+4\delta_6$ である。これを $1/4$ で変換すると $-1+4-4=-1=i^2$ となり、変換の積とも一致する。
:::
::::

::::exercise{level=basic}
$G=\Z/8\Z$、$H=\{0,2,4,6\}$ とする。$G$ に正規化測度、$H$ に数え上げ測度を入れたとき、商の全質量と $H^\perp$ の一点の双対質量を求めよ。
:::solution
$P_H1=4$ なので商積分公式に $1$ を代入すると、商の全質量は $1/4$ である。$H^\perp$ は $k=0,4$ の指標からなる。正規化測度のときは双対側の一点の質量が $1$ であり、元側の測度を $1/4$ 倍したため逆に双対側は $4$ 倍となる。
:::
::::

::::exercise{level=advanced}
$\Z$ 上の $f=\delta_0-\delta_2$ と $g=\delta_0+\delta_2$ の二つの平行移動の族の閉包が $\ell^1$ 全体か判定し、いずれか一つではどうなるかも答えよ。
:::solution
変換はそれぞれ $1-e^{-2it}$、$1+e^{-2it}$ である。前者は $t=0$ で、後者は $t=\pi/2$ でゼロとなるため、一つだけなら Wiener の定理により稠密でない。両者の和は定数 $2$ なので共通零点はない。また $(f+g)/2=\delta_0$ の平行移動で有限台列を全て生成できるため、両者の場合の閉包は $\ell^1$ 全体である。
:::
::::

::::exercise{level=advanced}
$\R$ の Lebesgue 測度で $f=1_{[0,a]}$、$a>0$ とする。自己相関の値 $b(a/3)$、Bochner 測度の全質量と密度を求めよ。
:::solution
$f$ と $\tau_xf$ の重なりは $(a-|x|)_+$ なので $b(a/3)=2a/3$ である。変換は $\xi\ne0$ で $e^{-\pi i a\xi}\sin(\pi a\xi)/(\pi\xi)$ である。従って :ref[prop-four-lca-autocorrelation] により密度は $\sin^2(\pi a\xi)/(\pi\xi)^2$ で、原点では極限 $a^2$ である。質量は $b(0)=a$ であり、Plancherel の $\|f\|_2^2=a$ とも一致する。
:::
::::

## 確認テスト

::::quiz{answer="4" tol="1e-6"}
$b(t)=2+2\cos t$ の Bochner 測度の全質量を求めよ。
:::explanation
円周上の指標係数は $0$ に $2$、$\pm1$ に $1$ ずつである。従って全質量は $4$ である。$t=0$ の値も同じである。
:::
::::

::::quiz
$G=\Z$ の双対はどれか。

- [ ] $\R$
- [x] $\mathbb T$
- [ ] $\Z$

:::explanation
指標は生成元 $1$ の像によって一意に決まり、その像には円上の任意の点を取れる。コンパクト集合は有限なので、像の収束は各有限集合での一様収束と同値である。
:::
::::

::::quiz{answer="1" tol="1e-6"}
$\Z$ の数え上げ測度に対応する双対 Haar 測度の全質量を求めよ。
:::explanation
入力 $\delta_0$ の二乗ノルムは $1$、変換は定数 $1$ である。Plancherel により定数 $1$ の二乗積分、すなわち双対測度の全質量は $1$ である。
:::
::::

::::quiz
$L^1$ 関数の移動の線形包が稠密であるための条件はどれか。

- [x] 変換が全ての指標で非零である。
- [ ] 変換の零点は一点までである。
- [ ] 変換が全ての指標で非負である。

:::explanation
一点でも零点があれば、そこでの変換評価が全移動を消す。零点がなければ各周波数の近くで局所逆元を作り、コンパクトな周波数台を持つ稠密な集合を近似できる。従って条件は全指標での非消滅である。
:::
::::

::::quiz{answer="3" tol="1e-6"}
$\Z/12\Z$ に正規化 Haar 測度を入れ、$H=\{0,4,8\}$ に数え上げ測度を入れる。商で用いる双対 Haar 測度の一点の重みを求めよ。
:::explanation
$P_H1=3$ なので商積分公式から商の全質量は $1/3$ である。商の Haar 測度が質量 $1$ なら双対は数え上げ測度である。商の測度を $1/3$ 倍したので、Plancherel の正規化では双対の各点の重みは $3$ になる。
:::
::::

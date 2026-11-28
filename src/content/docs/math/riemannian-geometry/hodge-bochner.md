---
title: 'Hodge 理論と Bochner の方法'
description: '微分形式の内積から Hodge 分解を導き、Bochner の公式でコホモロジーと固有値を制約する。'
notation:
  - symbol: '\star'
    meaning: 'Hodge スター作用素'
  - symbol: '\delta'
    meaning: '外微分の形式的随伴'
  - symbol: '\Delta_H'
    meaning: '非負符号の Hodge–Laplacian'
  - symbol: '\mathcal H^p'
    meaning: '調和 p 形式の空間'
  - symbol: '\nabla^*\nabla'
    meaning: '接続 Laplacian'
---

## この章の目標

- スター作用素と余微分の符号を計算する。
- 楕円型方程式の基本定理から Hodge 分解を得る。
- 調和形式を de Rham 類の一意な代表として使う。
- Bochner の公式で消滅定理と固有値の下界を示す。

## Hodge スター作用素と余微分

微分形式には、向きと計量を同時に使って次数を反転する演算がある。これにより外微分の積分に関する随伴を、同じく局所的な微分作用素として書ける。以下では向き付けられた閉 Riemann 多様体、すなわちコンパクトで境界のない多様体を基本とする。外微分と Stokes の定理は [微分形式](/math/manifolds/) を用いる。

:::definition[Hodge スター]{#def-rg-hodge-star}
$p$ 形式の点ごとの内積は、正規直交共標構 $\theta^1,\dots,\theta^n$ の増加添字の外積を正規直交とする。$\alpha\wedge\star\beta=g(\alpha,\beta)\omega_g$ を満たす写像 $\star:\Lambda^pT^*M\to\Lambda^{n-p}T^*M$ を :term[Hodge スター作用素]{en="Hodge star operator" yomi="ほっじすたーさようそ"} という。
:::

内積は $g(\alpha_1\wedge\dots\wedge\alpha_p,\beta_1\wedge\dots\wedge\beta_p)=\det(g(\alpha_i,\beta_j))$ と書けるため標構に依存しない。向き付けた正規直交標構では $\star\theta^I=\epsilon(I,I^c)\theta^{I^c}$ であり、符号は $I,I^c$ をこの順に並べる置換の符号である。従って存在と一意性、等長性、$\star^2=(-1)^{p(n-p)}$ が従う。向きを反転すると $\star$ は $-\star$ になる。

:::definition[余微分]{#def-rg-codifferential}
$p$ 形式上で $\delta=(-1)^{n(p+1)+1}\star d\star$ と定め、これを :term[余微分]{en="codifferential" yomi="よびぶん"} という。次数零では $\delta=0$ とする。
:::

$\alpha$ を $(p-1)$ 形式、$\beta$ を $p$ 形式とする。Stokes の定理により $0=\int d(\alpha\wedge\star\beta)=\int d\alpha\wedge\star\beta+(-1)^{p-1}\int\alpha\wedge d\star\beta$ なので、$\star\delta\beta=(-1)^p d\star\beta$ とすれば随伴条件が成り立つ。$\star^2$ を使うと定義の符号を得る。積分内積を $(\cdot,\cdot)_{L^2}$ と書けば $(d\alpha,\beta)_{L^2}=(\alpha,\delta\beta)_{L^2}$ である。また $d^2=0$ と $\star^2$ の式から $\delta^2=0$ である。

向き付け不可能な場合も、局所的な向きを選んで $\delta$ を定義できる。向きを変えると二つのスターの符号がともに反転して相殺するため、局所式は貼り合わさる。積分内積には計量の体積密度を使う。随伴性は、向きの二重被覆で形式を引き戻して Stokes の定理を使い、両辺を被覆次数 $2$ で割れば従う。従って $\delta,\Delta_H$ の式と Hodge 分解は向きなしで成立する。大域的なスターを用いる双対性には向きを仮定する。

:::example[三次元の符号]
正の共標構 $(dx,dy,dz)$ では $\star dx=dy\wedge dz$、$\star(dy\wedge dz)=dx$ である。$\beta=x\,dx+2y\,dy$ の余微分は $-\star d\star\beta=-3$ である。これは対応するベクトル場の発散の負に等しい。
:::

共形変化 $\widetilde g=e^{2u}g$ では $p$ 形式の内積が $e^{-2pu}$ 倍、体積形式が $e^{nu}$ 倍になるので $\widetilde\star=e^{(n-2p)u}\star$ である。例えば二次元の 1 形式におけるスターは共形変化で不変である。

:::example[異方的な計量のスター]
平面の $g=4dx^2+9dy^2$ で正の共標構は $(2dx,3dy)$ である。$\star(2dx)=3dy$、$\star(3dy)=-2dx$ だから $\star dx=\frac32dy$、$\star dy=-\frac23dx$ となる。$\star^2dx=-dx$ は Euclid 計量と同じだが、各方向の倍率は異なる。また体積形式は $6dx\wedge dy$ である。
:::

### 挿入と外積の随伴

1 形式 $\xi$ と計量双対 $X=\xi^\sharp$ に対し、$\xi\wedge$ の点ごとの随伴は $\iota_X$ である。正規直交余標構の $\xi=\theta^a$ の場合に確認すれば線形性で一般に従う。増加添字の基底形式 $\theta^I$ に $a$ が含まれると外積は零である。含まれないとき、$\theta^a\wedge\theta^I$ を順に並べる符号は、$\iota_{e_a}\theta^{I\cup\{a\}}$ で $a$ を取り除く符号と同じである。従って基底ごとの内積が一致する。この随伴関係から、両空間の作用素に
$$
\iota_X(\xi\wedge\alpha)+\xi\wedge \iota_X\alpha=|\xi|^2\alpha
$$
が成り立つ。積への挿入 $\iota_X(\xi\wedge\alpha)=\xi(X)\alpha-\xi\wedge\iota_X\alpha$ を展開すればよい。

この事実を使うと、余微分の積の公式は
$$
\delta(f\alpha)=f\delta\alpha-\iota_{\grad f}\alpha
$$
と書ける。証明は $(f\alpha,d\beta)=(\alpha,d(f\beta)-df\wedge\beta)$ と積分内積で計算し、両項に随伴性を使うだけである。この等式は試験形式 $\beta$ の台を座標域内に取れるので局所公式でもある。例えば Euclid 平面の $\alpha=x\,dy$ は余閉であり、$\delta(y^2\alpha)=-\iota_{2y\partial_y}(x\,dy)=-2xy$ となる。一方、$d(y^2\alpha)=y^2dx\wedge dy$ である。閉性と余閉性は異なる微分方向を調べている。

### 共形変化した余微分

$p$ 形式への最初のスターには $e^{(n-2p)u}$、その外微分 $(n-p+1)$ 形式への次のスターには $e^{(n-2(n-p+1))u}$ が掛かる。従って指数の和は $-2u$ であり、余微分の積の公式と合わせて
$$
\widetilde\delta\alpha=e^{-2u}\left(\delta\alpha-(n-2p)\iota_{\grad u}\alpha\right)
$$
となる。実際、余微分の式の内側に $e^{(n-2p)u}\alpha$ を代入し、外側の倍率を掛けて、$d(e^{(n-2p)u})=(n-2p)e^{(n-2p)u}du$ を整理した。特に $n=2p$ で閉かつ余閉の条件は共形変化で同じである。積分内積も $e^{(n-2p)u}=1$ だから不変になる。

:::example[同じ 1 形式を二次元と三次元で測る]
Euclid の共標構の中で $\alpha=dx$、$u=x$ とする。元の $\delta dx=0$ で $\iota_{\grad u}dx=1$ である。二次元では $\widetilde\delta dx=0$ だが、三次元では $-e^{-2x}$ となる。実際、三次元での $dx$ の双対ベクトルは $e^{-2x}\partial_x$、密度は $e^{3x}$ なので、発散は $e^{-3x}\partial_x e^x=e^{-2x}$ である。このように、同じ座標関数の微分でも、体積と逆計量の倍率が相殺するかどうかは次元で変わる。
:::

## Hodge–Laplacian

:::definition[Hodge–Laplacian]{#def-rg-hodge-laplacian}
微分形式上の作用素 $\Delta_H=d\delta+\delta d$ を :term[Hodge–Laplacian]{en="Hodge Laplacian" yomi="ほっじらぷらしあん"} という。核の滑らかな元を調和形式と呼ぶ。
:::

随伴性から $(\Delta_H\alpha,\alpha)=\norm{d\alpha}_{L^2}^2+\norm{\delta\alpha}_{L^2}^2$ である。よって調和であることは閉かつ余閉であることと同値である。$d^2=\delta^2=0$ を展開すれば $d,\delta$ は $\Delta_H$ と可換である。スターは閉性と余閉性を交換して調和形式を調和形式に送る。関数上では $\Delta_H=\delta d=-\divg\grad=-\Delta_g$ である。実際 $\star df$ は $\operatorname{grad} f$ を体積形式に挿入した形であり、その外微分は発散に体積形式を掛けたものになる。

### 主部

正規直交共標構で $\varepsilon_i\alpha=\theta^i\wedge\alpha$、$\iota_i=\iota_{e_i}$ とする。挿入の積の規則により $\varepsilon_i\iota_j+\iota_j\varepsilon_i=\delta_{ij}\id$ である。捩率零のため $d=\sum_i\varepsilon_i\nabla_i$ であり、$d\alpha$ の定義の Lie 括弧項を接続で置き換えると各入力の接続項が相殺してこれを得る。$\delta=-\sum_i\iota_i\nabla_i$ は随伴性から従う。点で正規標構を取り、計量接続で部分積分すれば主部がこの式となり、標構の発散と接続係数の補正はその点で零である。従って $\Delta_H$ の二階の主部は $-\sum_i\nabla_i\nabla_i$ で、主シンボルは $|\xi|_g^2\id$ である。これが楕円性を与える。

:::example[平坦トーラス]
$\R^2/(2\pi\Z)^2$ の平行共標構 $dx,dy$ では $\Delta_H$ は各係数に $-\partial_x^2-\partial_y^2$ を作用させる。例えば $\cos(2x)\sin y\,dx$ の固有値は $2^2+1^2=5$ である。
:::

:::example[閉じていても調和でない形式]
周長 $2\pi$ の円の弧長座標で $\alpha=\cos t\,dt$ とする。$d\alpha=0$ は二次の形式がないためだが、$\delta\alpha=-\partial_t(\cos t)=\sin t$ であり余閉ではない。従って調和ではない。実際 $\alpha=d(\sin t)$ で完全成分に属し、$\Delta_H\alpha=\alpha$ である。
:::

### スターとの可換性とエネルギー

随伴性の導出に用いた $\star\delta=(-1)^p d\star$ を $p$ 形式に適用し、$\star^2$ から対応する $\star d=(-1)^{p+1}\delta\star$ を得る。実際、$\delta\star\alpha=(-1)^{n(n-p+1)+1}\star d\star\star\alpha$ の符号指数を整理すると $p+1$ と同じ偶奇になる。従って
$$
\star d\delta\alpha=\delta d\star\alpha,\qquad
\star\delta d\alpha=d\delta\star\alpha
$$
である。最初の式では $\delta\alpha$ の次数は $p-1$ なので、二つの符号はともに $(-1)^p$ で相殺する。二番目は $d\alpha$ の次数が $p+1$ で符号が相殺する。ゆえに $\star\Delta_H=\Delta_H\star$ であり、核だけでなく正の固有空間も次数 $p$ と $n-p$ で対応する。

:::example[体積形式の係数を微分する]
向き付き閉多様体で $\alpha=f\omega_g$ とする。$\star\alpha=f$ なので $\Delta_H(f\omega_g)=(-\Delta_g f)\omega_g$ である。単位二次元球面の高さ $f=\cos\theta$ を使うと $\Delta_H(f\omega_g)=2f\omega_g$ となる。また $df$ は $d$ と $\Delta_H$ の可換性により同じ固有値 $2$ を持つ。$d(f\omega_g)=0$ は最高次であるため自動だが、$f$ が定数でないので調和ではない。
:::

## Hodge の分解定理

:::theorem[Hodge 分解]{#thm-rg-hodge-decomposition}
閉 Riemann 多様体で任意の滑らかな $p$ 形式は
$$
\alpha=h+d\beta+\delta\gamma,\qquad h\in\mathcal H^p,
\label{eq-rg-hodge-decomposition}
$$
と $L^2$ に関して直交分解する。三つの成分は一意であり、$\mathcal H^p$ は有限次元である。
:::

:::proof
解析の入力として、閉多様体上の一様楕円型二階作用素の [楕円型正則性](/math/pde/elliptic-regularity/)、:ref[thm-fa-rellich-kondrachov]{text="Rellich のコンパクト埋め込み"}、:ref[cor-fa-lax-milgram]{text="Lax–Milgram の定理"} を使う。ここでは一次のエネルギー $Q(u)=\norm{du}_{L^2}^2+\norm{\delta u}_{L^2}^2$ に対し $\norm u_{H^1}^2\le C(Q(u)+\norm u_{L^2}^2)$ が成り立つ。理由は、コンパクト台を持つ Euclid 座標で定数係数に凍結すると、外積と挿入の反交換関係から Fourier シンボルの二乗和が $|\xi|^2\id$ となるためである。座標球を縮めて主係数の変動を左辺に吸収し、零次係数と分割関数の微分を $\norm u_{L^2}^2$ で押さえる。有限個の座標球で足せば大域の評価となる。

$Q=0$ の弱解は楕円型正則性で滑らかである。主シンボルが計量の二次形式と恒等行列の積なので、局所標構の成分には同じスカラー楕円作用素が作用する。成分どうしを結ぶのは高々一階の項であり、成分ごとの正則性の評価に移して Sobolev 階数を上げられる。その $L^2$ 単位球は評価で $H^1$ に有界で、Rellich により $L^2$ で相対コンパクトである。無限次元なら正規直交列があって相対コンパクト性に反するので、核は有限次元である。さらに核に直交する $u$ では $Q(u)\ge c\norm u_{L^2}^2$ である。成り立たなければ $\norm{u_j}_{L^2}=1$、$Q(u_j)\to0$ で核に直交する列を選べる。$H^1$ で有界なので弱収束部分列と $L^2$ 強収束を得る。弱下半連続性により極限は $Q=0$、強収束により単位ノルムかつ核に直交しており矛盾する。

$h$ を $\alpha$ の核への $L^2$ 射影とする。上の評価を合わせると、核に直交する $H^1$ 空間で $Q$ の双線形形式は強圧的である。従って Lax–Milgram により $Q(u,v)=(\alpha-h,v)$ を全ての $v$ に対して満たす核に直交した弱解 $u$ がある。核の方向では両辺が零なので試験空間を全 $H^1$ に拡張できる。楕円型正則性により $u$ は滑らかであり、$\alpha=h+\Delta_Hu$ である。$\beta=\delta u,\gamma=du$ と選ぶ。調和成分と残りの直交性は随伴性で、$d\beta$ と $\delta\gamma$ の直交性は $d^2=0$ で従う。直交分解だから各成分は一意である。
:::

:::example[円周の分解]
周長 $2\pi$ の円では $\alpha=(2+\cos t)\,dt$ の調和成分は $2\,dt$ であり、残りは $d(\sin t)$ である。$2$ 形式がないため余完全成分は零である。積分を取ると調和成分の係数は周期平均で求められる。
:::

:::example[トーラスで三成分を確認する]
平坦正方形トーラスで $\alpha=(2+\cos x)dx+\sin x\,dy$ とする。調和成分は $2dx$、完全成分は $d(\sin x)=\cos x\,dx$ である。$\gamma=\cos x\,dx\wedge dy$ とすると $\delta\gamma=-\partial_x(\cos x)\,dy=\sin x\,dy$ なので、これが余完全成分になる。完全成分は $x$ 方向、余完全成分は $y$ 方向で点ごとに直交する。この点ごとの直交は例の特別な性質であり、Hodge 分解が一般に保証するのは積分内積での直交である。
:::

### Green 作用素と射影

:ref[thm-rg-hodge-decomposition] の証明で作った一意な核直交の解を $G_p\alpha$ と書き、調和射影を $H_p$ とする。この $G_p$ を :term[Green 作用素]{en="Green operator" yomi="ぐりーんさようそ"} という。定義から $\Delta_HG_p=I-H_p$ である。逆順についても $G_p\Delta_H\alpha=\alpha-H_p\alpha$ となる。なぜなら $\Delta_H\alpha$ は随伴性で核に直交し、$\alpha-H_p\alpha$ がその核直交解だからである。

$dG_p=G_{p+1}d$ を示す。$\Delta_H dG_p\alpha=d(\alpha-H_p\alpha)=d\alpha$ であり、完全形式 $dG_p\alpha$ は随伴性で調和形式に直交する。従って核直交解の一意性でこの式を得る。$\delta G_p=G_{p-1}\delta$ も同様で、$\delta H_p\alpha=0$ を使う。従って完全成分の射影は $d\delta G_p$、余完全成分の射影は $\delta dG_p$ である。これらは自己共役である。実際、直交分解の成分として他の二成分との内積が零であり、$(P\alpha,\beta)=(P\alpha,P\beta)=(\alpha,P\beta)$ が成立する。ポテンシャル $\beta,\gamma$ を任意に選ぶと核の自由度があるが、ここでは $G_p$ を使って具体的な一組を得た。

:::example[一つの Fourier 周波数の完全成分]
平坦トーラス上で $\alpha=\cos(2x+y)(dx+2dy)$ とする。固有値は $5$ なので $G_1\alpha=\alpha/5$ である。$\delta\alpha=(2+2)\sin(2x+y)=4\sin(2x+y)$ だから
$$
d\delta G_1\alpha=\tfrac45\cos(2x+y)(2dx+dy).
$$
差は $\cos(2x+y)(-3dx+6dy)/5$ であり、係数ベクトル $(-3,6)$ は波数 $(2,1)$ と直交するので余閉である。実際、この差は $\delta(-3\sin(2x+y)dx\wedge dy/5)$ である。$\delta(h\,dx\wedge dy)=h_y\,dx-h_x\,dy$ を用いると両成分が一致する。周波数が非零なので調和成分はない。
:::

:::example[異方的な計量での Hodge 分解]
$x,y$ がともに周期 $2\pi$ の二次元トーラスに、定数の計量行列 $G=\begin{pmatrix}4 & 0\\0 & 9\end{pmatrix}$ を入れる。ここでは $\phi=x$、$\alpha=\cos\phi(c_1\,dx+c_2\,dy)$、$c=\begin{pmatrix}2\\3\end{pmatrix}$ とする。このとき $d\phi$ の成分は $\xi=\begin{pmatrix}1\\0\end{pmatrix}$ であり、コサインの周期境界で値がつながる。

まず $d\phi$ に対応する計量双対ベクトルと二乗長を求める。
$$
\begin{aligned}
\xi^\sharp&=G^{-1}\xi=\begin{pmatrix}\frac{1}{4}\\0\end{pmatrix}\\
\norm\xi_g^2&=\xi^\top G^{-1}\xi=\frac{1}{4}\\
b=c(\xi^\sharp)&=\frac{1}{2}
\end{aligned}
$$
である。定数計量での Levi-Civita 記号は零であり、余微分は $-\sum g^{jk}\partial_j\alpha_k$ になる。従って $\delta\alpha=\frac{1}{2}\sin\phi$ である。一方、$d(\cos\phi)=-\sin\phi\,d\phi$ だから外微分は $-3\sin\phi\,dx\wedge dy$ になる。この二つの微分の係数は、波数と係数ベクトルの組み方によって別々に消える。

Laplacian は定数座標の各成分に $-\sum_{j,k}g^{jk}\partial_j\partial_k$ を作用させる。二階微分で $-\xi_j\xi_k$ を掛けるので $\Delta_H\alpha=\norm\xi_g^2\alpha=\frac{1}{4}\alpha$ である。調和成分は平均によって決まり、ここでは $\cos\phi$ の平均が零なので消える。実際、少なくとも一方の波数が非零整数だから、該当する変数を一周積分すると零になる。よって Green 作用素はこの周波数では $G_1\alpha=\alpha/\frac{1}{4}$ である。

完全成分は $d\delta G_1\alpha$ なので、その係数ベクトルは $c_{\rm ex}=b\xi/\norm\xi_g^2$ である。残りの係数は $c_{\rm co}=c-c_{\rm ex}$ である。指定した数値では
$$
\begin{aligned}
c_{\rm ex}&=\begin{pmatrix}2\\0\end{pmatrix}\\
c_{\rm co}&=\begin{pmatrix}0\\3\end{pmatrix}
\end{aligned}
$$
であり、$c_{\rm co}\cdot \xi^\sharp=0$ だから余閉である。調和成分はすでに除かれているため、残りは Hodge 分解の余完全成分と一致する。このように「波数方向への射影」は余接空間の計量 $G^{-1}$ による射影である。Euclid の成分の積で射影すると、計量に交差項がある場合には結果が異なる。

体積は $(2\pi)^2\sqrt{\det G}=24 \pi^{2}$ であり、$\cos^2\phi$ の平均は $1/2$ である。従って完全成分と余完全成分の積分ノルムの二乗はそれぞれ $12 \pi^{2}$ と $12 \pi^{2}$ となる。両者の和は元の形式の積分ノルムの二乗 $24 \pi^{2}$ である。この例では各係数が $\cos\phi$ に比例するので点ごとの直交性まで成立する。複数の周波数を重ねるときは、積分で周波数の直交性を使う。
:::

:::example[交差項を持つ計量での勾配方向の Hodge 分解]
$x,y$ がともに周期 $2\pi$ の二次元トーラスに、定数の計量行列 $G=\begin{pmatrix}3 & 1\\1 & 2\end{pmatrix}$ を入れる。ここでは $\phi=x + 2 y$、$\alpha=\cos\phi(c_1\,dx+c_2\,dy)$、$c=\begin{pmatrix}4\\-1\end{pmatrix}$ とする。このとき $d\phi$ の成分は $\xi=\begin{pmatrix}1\\2\end{pmatrix}$ であり、コサインの周期境界で値がつながる。

まず $d\phi$ に対応する計量双対ベクトルと二乗長を求める。
$$
\begin{aligned}
\xi^\sharp&=G^{-1}\xi=\begin{pmatrix}0\\1\end{pmatrix}\\
\norm\xi_g^2&=\xi^\top G^{-1}\xi=2\\
b=c(\xi^\sharp)&=-1
\end{aligned}
$$
である。定数計量での Levi-Civita 記号は零であり、余微分は $-\sum g^{jk}\partial_j\alpha_k$ になる。従って $\delta\alpha=-1\sin\phi$ である。一方、$d(\cos\phi)=-\sin\phi\,d\phi$ だから外微分は $9\sin\phi\,dx\wedge dy$ になる。この二つの微分の係数は、波数と係数ベクトルの組み方によって別々に消える。

Laplacian は定数座標の各成分に $-\sum_{j,k}g^{jk}\partial_j\partial_k$ を作用させる。二階微分で $-\xi_j\xi_k$ を掛けるので $\Delta_H\alpha=\norm\xi_g^2\alpha=2\alpha$ である。調和成分は平均によって決まり、ここでは $\cos\phi$ の平均が零なので消える。実際、少なくとも一方の波数が非零整数だから、該当する変数を一周積分すると零になる。よって Green 作用素はこの周波数では $G_1\alpha=\alpha/2$ である。

完全成分は $d\delta G_1\alpha$ なので、その係数ベクトルは $c_{\rm ex}=b\xi/\norm\xi_g^2$ である。残りの係数は $c_{\rm co}=c-c_{\rm ex}$ である。指定した数値では
$$
\begin{aligned}
c_{\rm ex}&=\begin{pmatrix}- \frac{1}{2}\\-1\end{pmatrix}\\
c_{\rm co}&=\begin{pmatrix}\frac{9}{2}\\0\end{pmatrix}
\end{aligned}
$$
であり、$c_{\rm co}\cdot \xi^\sharp=0$ だから余閉である。調和成分はすでに除かれているため、残りは Hodge 分解の余完全成分と一致する。このように「波数方向への射影」は余接空間の計量 $G^{-1}$ による射影である。Euclid の成分の積で射影すると、計量に交差項がある場合には結果が異なる。

体積は $(2\pi)^2\sqrt{\det G}=4 \sqrt{5} \pi^{2}$ であり、$\cos^2\phi$ の平均は $1/2$ である。従って完全成分と余完全成分の積分ノルムの二乗はそれぞれ $\sqrt{5} \pi^{2}$ と $\frac{81 \sqrt{5} \pi^{2}}{5}$ となる。両者の和は元の形式の積分ノルムの二乗 $\frac{86 \sqrt{5} \pi^{2}}{5}$ である。この例では各係数が $\cos\phi$ に比例するので点ごとの直交性まで成立する。複数の周波数を重ねるときは、積分で周波数の直交性を使う。
:::

:::example[一つの周波数が完全になる場合の Hodge 分解]
$x,y$ がともに周期 $2\pi$ の二次元トーラスに、定数の計量行列 $G=\begin{pmatrix}2 & -1\\-1 & 3\end{pmatrix}$ を入れる。ここでは $\phi=2 x + y$、$\alpha=\cos\phi(c_1\,dx+c_2\,dy)$、$c=\begin{pmatrix}6\\3\end{pmatrix}$ とする。このとき $d\phi$ の成分は $\xi=\begin{pmatrix}2\\1\end{pmatrix}$ であり、コサインの周期境界で値がつながる。

まず $d\phi$ に対応する計量双対ベクトルと二乗長を求める。
$$
\begin{aligned}
\xi^\sharp&=G^{-1}\xi=\begin{pmatrix}\frac{7}{5}\\\frac{4}{5}\end{pmatrix}\\
\norm\xi_g^2&=\xi^\top G^{-1}\xi=\frac{18}{5}\\
b=c(\xi^\sharp)&=\frac{54}{5}
\end{aligned}
$$
である。定数計量での Levi-Civita 記号は零であり、余微分は $-\sum g^{jk}\partial_j\alpha_k$ になる。従って $\delta\alpha=\frac{54}{5}\sin\phi$ である。一方、$d(\cos\phi)=-\sin\phi\,d\phi$ だから外微分は 零 になる。この二つの微分の係数は、波数と係数ベクトルの組み方によって別々に消える。

Laplacian は定数座標の各成分に $-\sum_{j,k}g^{jk}\partial_j\partial_k$ を作用させる。二階微分で $-\xi_j\xi_k$ を掛けるので $\Delta_H\alpha=\norm\xi_g^2\alpha=\frac{18}{5}\alpha$ である。調和成分は平均によって決まり、ここでは $\cos\phi$ の平均が零なので消える。実際、少なくとも一方の波数が非零整数だから、該当する変数を一周積分すると零になる。よって Green 作用素はこの周波数では $G_1\alpha=\alpha/\frac{18}{5}$ である。

完全成分は $d\delta G_1\alpha$ なので、その係数ベクトルは $c_{\rm ex}=b\xi/\norm\xi_g^2$ である。残りの係数は $c_{\rm co}=c-c_{\rm ex}$ である。指定した数値では
$$
\begin{aligned}
c_{\rm ex}&=\begin{pmatrix}6\\3\end{pmatrix}\\
c_{\rm co}&=\begin{pmatrix}0\\0\end{pmatrix}
\end{aligned}
$$
であり、$c_{\rm co}\cdot \xi^\sharp=0$ だから余閉である。調和成分はすでに除かれているため、残りは Hodge 分解の余完全成分と一致する。このように「波数方向への射影」は余接空間の計量 $G^{-1}$ による射影である。Euclid の成分の積で射影すると、計量に交差項がある場合には結果が異なる。

体積は $(2\pi)^2\sqrt{\det G}=4 \sqrt{5} \pi^{2}$ であり、$\cos^2\phi$ の平均は $1/2$ である。従って完全成分と余完全成分の積分ノルムの二乗はそれぞれ $\frac{324 \sqrt{5} \pi^{2}}{5}$ と $0$ となる。両者の和は元の形式の積分ノルムの二乗 $\frac{324 \sqrt{5} \pi^{2}}{5}$ である。この例では各係数が $\cos\phi$ に比例するので点ごとの直交性まで成立する。複数の周波数を重ねるときは、積分で周波数の直交性を使う。
:::

## 調和形式と de Rham コホモロジー

:::theorem[調和代表の一意性]{#thm-rg-harmonic-representative}
閉 Riemann 多様体上で各 :ref[def-mfd-de-rham-cohomology]{text="de Rham コホモロジー"} 類は一意の調和代表を持つ。従って $\mathcal H^p\simeq H_{\mathrm{dR}}^p(M)$ である。
:::

:::proof
閉形式 $\alpha=h+d\beta+\delta\gamma$ に $d$ を掛けると $d\delta\gamma=0$ である。随伴性から $\norm{\delta\gamma}_{L^2}^2=(\gamma,d\delta\gamma)=0$ なので $\alpha-h=d\beta$ であり、$h$ は同じ類の代表である。調和形式 $h$ が完全形式 $d\beta$ でもあるなら、$\norm h^2=(\delta h,\beta)=0$ より零である。これで一意性が従う。
:::

向き付けられた閉多様体では $\star$ が $\mathcal H^p$ と $\mathcal H^{n-p}$ を同型で結ぶ。非零の調和形式 $h$ に対して $\int h\wedge\star h=\norm h_{L^2}^2>0$ なので、コホモロジーの積と積分の双対性も非退化である。積分は代表を変えても Stokes により不変であり、これが実係数の Poincaré 双対性の解析的な証明である。

:::example[平坦トーラスの調和 1 形式]
座標 $x,y$ の二つの周期を持つ平坦トーラスでは $\Delta_H$ は係数ごとの負の Laplacian である。調和な係数はコンパクト連結多様体上で一定なので、全調和 1 形式は $a\,dx+b\,dy$ である。従って第一 Betti 数は $2$ である。他の周期の座標を選ぶと、固定する周期積分に応じて係数 $a,b$ が変わる。
:::

:::example[調和代表が計量に依存する]
$S^1=\R/(2\pi\Z)$ に $g=(2+\cos t)^2dt^2$ を入れる。弧長共標構 $\theta=(2+\cos t)dt$ は平行であり、調和 1 形式はその定数倍である。$[dt]$ と同じ類なら一周の積分が $2\pi$ なので、$\int\theta=4\pi$ より調和代表は $\theta/2$ である。実際 $\theta/2-dt=\frac12\cos t\,dt=d(\sin t/2)$ であり、同じ類であることを直接確認できる。
:::

### 調和代表の最小エネルギー性

閉形式 $\alpha$ と同じ類を持つ調和形式を $h$ とすると、任意の他の代表は $h+d\beta$ と書ける。随伴性から $(h,d\beta)=(\delta h,\beta)=0$ なので
$$
\norm{h+d\beta}_{L^2}^2=\norm h_{L^2}^2+\norm{d\beta}_{L^2}^2.
$$
従って $h$ は同じコホモロジー類で積分ノルムを最小にする一意の形式である。逆に閉形式がその類で停留するなら、任意の $\beta$ について $(\alpha,d\beta)=0$ である。$\beta=\delta\alpha$ を選べば $\delta\alpha=0$ となる。従って閉性に加え余閉性が得られ、調和性と同値になる。

:::example[長方形と周期の正規化]
辺の長さが $a,b$ の平坦トーラスで第一、第二の基本円周の周期をそれぞれ $1,0$ にする調和代表は $dx/a$ である。その点ごとの二乗長は $1/a^2$、全体積は $ab$ なので $L^2$ ノルムの二乗は $b/a$ である。例えば $a=2,b=3$ なら $3/2$ である。$dx/a+d(\sin(2\pi x/a))$ は周期を変えない。そのノルムの二乗との差は
$$
\int_0^b\!\int_0^a \left(\frac{2\pi}{a}\right)^2\cos^2(2\pi x/a)\,dx\,dy
=2\pi^2 b/a
$$
である。周期データを保ったまま非定数係数を加えると、この分だけエネルギーが増える。
:::

:::example[余閉な方向が座標に沿わない場合の Hodge 分解]
$x,y$ がともに周期 $2\pi$ の二次元トーラスに、定数の計量行列 $G=\begin{pmatrix}2 & 1\\1 & 2\end{pmatrix}$ を入れる。ここでは $\phi=x$、$\alpha=\cos\phi(c_1\,dx+c_2\,dy)$、$c=\begin{pmatrix}1\\2\end{pmatrix}$ とする。このとき $d\phi$ の成分は $\xi=\begin{pmatrix}1\\0\end{pmatrix}$ であり、コサインの周期境界で値がつながる。

まず $d\phi$ に対応する計量双対ベクトルと二乗長を求める。
$$
\begin{aligned}
\xi^\sharp&=G^{-1}\xi=\begin{pmatrix}\frac{2}{3}\\- \frac{1}{3}\end{pmatrix}\\
\norm\xi_g^2&=\xi^\top G^{-1}\xi=\frac{2}{3}\\
b=c(\xi^\sharp)&=0
\end{aligned}
$$
である。定数計量での Levi-Civita 記号は零であり、余微分は $-\sum g^{jk}\partial_j\alpha_k$ になる。従って $\delta\alpha=0\sin\phi$ である。一方、$d(\cos\phi)=-\sin\phi\,d\phi$ だから外微分は $-2\sin\phi\,dx\wedge dy$ になる。この二つの微分の係数は、波数と係数ベクトルの組み方によって別々に消える。

Laplacian は定数座標の各成分に $-\sum_{j,k}g^{jk}\partial_j\partial_k$ を作用させる。二階微分で $-\xi_j\xi_k$ を掛けるので $\Delta_H\alpha=\norm\xi_g^2\alpha=\frac{2}{3}\alpha$ である。調和成分は平均によって決まり、ここでは $\cos\phi$ の平均が零なので消える。実際、少なくとも一方の波数が非零整数だから、該当する変数を一周積分すると零になる。よって Green 作用素はこの周波数では $G_1\alpha=\alpha/\frac{2}{3}$ である。

完全成分は $d\delta G_1\alpha$ なので、その係数ベクトルは $c_{\rm ex}=b\xi/\norm\xi_g^2$ である。残りの係数は $c_{\rm co}=c-c_{\rm ex}$ である。指定した数値では
$$
\begin{aligned}
c_{\rm ex}&=\begin{pmatrix}0\\0\end{pmatrix}\\
c_{\rm co}&=\begin{pmatrix}1\\2\end{pmatrix}
\end{aligned}
$$
であり、$c_{\rm co}\cdot \xi^\sharp=0$ だから余閉である。調和成分はすでに除かれているため、残りは Hodge 分解の余完全成分と一致する。このように「波数方向への射影」は余接空間の計量 $G^{-1}$ による射影である。Euclid の成分の積で射影すると、計量に交差項がある場合には結果が異なる。

体積は $(2\pi)^2\sqrt{\det G}=4 \sqrt{3} \pi^{2}$ であり、$\cos^2\phi$ の平均は $1/2$ である。従って完全成分と余完全成分の積分ノルムの二乗はそれぞれ $0$ と $4 \sqrt{3} \pi^{2}$ となる。両者の和は元の形式の積分ノルムの二乗 $4 \sqrt{3} \pi^{2}$ である。この例では各係数が $\cos\phi$ に比例するので点ごとの直交性まで成立する。複数の周波数を重ねるときは、積分で周波数の直交性を使う。
:::

## Bochner の公式と Weitzenböck の公式

接続 Laplacian $\nabla^*\nabla=-\tr_g\nabla^2$ と $\Delta_H$ は同じ二階主部を持つが、微分の順序が違うため差に曲率が残る。

:::proposition[Weitzenböck の公式]{#prop-rg-weitzenbock}
$c_i=\varepsilon_i-\iota_i$ とし、$R^\Lambda_{ij}$ を $p$ 形式への誘導接続の曲率作用とすると
$$
\Delta_H=\nabla^*\nabla+\mathcal Q_p,\qquad
\mathcal Q_p=\sum_{i<j}c_ic_j R^\Lambda_{ij}.
\label{eq-rg-weitzenbock}
$$
特に $p=1$ では $\mathcal Q_1$ は Ricci 作用素であり、定断面曲率 $k$ では $\mathcal Q_p=k p(n-p)\id$ である。
:::

:::proof
$d+\delta=\sum_i c_i\nabla_i$、$c_ic_j+c_jc_i=-2\delta_{ij}$ である。点で正規標構を取り二乗すると、対角は $-\sum_i\nabla_i\nabla_i$、異なる添字は $c_ic_j[\nabla_i,\nabla_j]$ である。そこで括弧 $[e_i,e_j]$ は零なので、後者は曲率である。一方 $(d+\delta)^2=\Delta_H$ である。形式上では $R^\Lambda_{ij}$ は各引数に $-R(e_i,e_j)$ を作用させた和である。

1 形式では $R^\Lambda_{ij}$ の後に $c_ic_j$ を掛けたとき、次数 $1$ の成分は $-\varepsilon_i\iota_j+\varepsilon_j\iota_i$ による。その $e^i$ 係数を足すと $-\sum_j(R^\Lambda_{ij}\alpha)(e_j)=\sum_j\alpha(R(e_i,e_j)e_j)$ であり、曲率の対交換と Ricci の定義により $\alpha(\mathrm{Ric}^{\sharp}e_i)$ である。全作用素は次数を保つのでこれが全成分である。

定曲率では $(i,j)$ の回転作用は $e^i$ を $-k e^j$ へ、$e^j$ を $k e^i$ へ移し、他の共標構を零に送る。外積基底 $e^I$ が両方を含むか両方を含まなければ寄与は零である。ちょうど一つを含む場合は $R^\Lambda_{ij}$ がその因子を交換し、$c_ic_j$ が逆に交換するので $k e^I$ が残る。その種類の組は $p(n-p)$ 個である。
:::

:::proposition[Bochner の公式]{#prop-rg-bochner-function}
関数 $f$ に対し
$$
\tfrac12\Delta_g\norm{\grad f}^2
=\norm{\operatorname{Hess}f}^2+
g(\grad f,\grad\Delta_g f)+\mathrm{Ric}(\grad f,\grad f).
\label{eq-rg-bochner-function}
$$
:::

:::proof
任意の 1 形式 $\alpha$ で $\frac12\Delta_g|\alpha|^2=\norm{\nabla\alpha}^2-g(\nabla^*\nabla\alpha,\alpha)$ である。両立接続で二回積の微分をし正規標構で足せば従う。$\alpha=df$ に :eqref[eq-rg-weitzenbock] を使い、$\Delta_H df=d\Delta_H f=-d\Delta_g f$ と $g(df,\cdot)$ の同一視を用いれば主張を得る。
:::

:::example[Euclid 二次関数]
$f=(x^2+y^2)/2$ なら $\norm{\grad f}^2=x^2+y^2$、左辺は $2$ である。右辺で Hessian は単位行列なので二乗ノルム $2$、$\Delta_g f=2$ は定数、曲率は零である。
:::

:::example[球面での曲率項の寄与]
$S^2(a)$ の高さ関数 $f$ について $\Delta_Hdf=(2/a^2)df$ である。$\mathcal Q_1$ は Ricci 作用なので $\mathcal Q_1df=(1/a^2)df$ となり、差から $\nabla^*\nabla df=(1/a^2)df$ が従う。例えば $a=2$ では Hodge 固有値 $1/2$ と接続 Laplacian の固有値 $1/4$ が異なる。
:::

### 定曲率での全次数の曲率項

$g$ の断面曲率が $k$ なら、Weitzenböck の曲率項は $p$ 形式上で $k\,p(n-p)\id$ になる。正規点で $\varepsilon_i$ と $\iota_i$ を使うと $\Delta_H=\nabla^*\nabla-\sum_{i,j}\varepsilon_i\iota_j R^\Lambda(e_i,e_j)$ である。実際、$d\delta+\delta d$ を展開し、$\iota_j\varepsilon_i=\delta_{ij}-\varepsilon_i\iota_j$ を用いると反対称部分が共変微分の交換子となる。$R^\Lambda$ は 共変な引数のそれぞれに $-R$ を作用させるので、定曲率では $R^\Lambda(e_i,e_j)=k(\varepsilon_i\iota_j-\varepsilon_j\iota_i)$ である。これは 1 形式の基底上で $R(e_i,e_j)e_j=ke_i$ と計量両立性から確かめられ、外積への微分の規則で全次数に拡張される。

$i\ne j$ では $(\varepsilon_i\iota_j)^2=0$、$\varepsilon_i\iota_j\varepsilon_j\iota_i$ は、基底形式 $\theta^I$ について $i\in I,j\notin I$ の場合に恒等、それ以外で零となる。最初の操作で $i$ を外して $j$ を入れ、次の操作で逆に戻すため符号は二回同じである。従って $\theta^I$ に寄与する順序対の個数は $p(n-p)$ であり、所望の公式が得られる。$p=1$ では曲率項は $(n-1)k$ で Ricci と一致し、$p=0,n$ では消える。

:::example[勾配の粗い Laplacian と曲率項を分ける]
半径 $2$ の $S^3$ の高さの正規化 $f=x_4/2$ を取る。球面の第二基本形式から $\operatorname{Hess}f=-fg/4$ である。跡を取ると $-\Delta_g f=3f/4$ なので、$d$ との可換性から 1 形式 $\alpha=df$ は $\Delta_H\alpha=3\alpha/4$ を満たす。一方、この球面の $\mathrm{Ric}=g/2$ だから、Weitzenböck の式で粗い Laplacian の項は
$$
\nabla^*\nabla \alpha=(3/4-1/2)\alpha=\alpha/4
$$
となる。符号と尺度を積分でも確かめる。$\norm{\nabla\alpha}^2=\norm{\operatorname{Hess}f}^2=3f^2/16$ であり、関数の部分積分から $\int\norm\alpha^2=(3/4)\int f^2$ である。従って $\int\norm{\nabla\alpha}^2=(1/4)\int\norm\alpha^2$ となる。全固有値 $3/4$ のうち二階の共変微分だけが与える値は $1/4$ で、残る $1/2$ は曲率から来る。
:::

## Bochner の消滅定理

:::theorem[Bochner の消滅定理]{#thm-rg-bochner-vanishing}
閉連結多様体で $\mathrm{Ric}\ge0$ なら全ての調和 1 形式は平行である。この仮定のもとでさらにある一点で $\mathrm{Ric}>0$ なら $H^1_{\mathrm{dR}}(M)=0$ である。一般に $\mathcal Q_p$ が全点で半正定値なら調和 $p$ 形式は平行であり、さらにある一点で正定値なら $H^p_{\mathrm{dR}}(M)=0$ である。
:::

:::proof
調和形式 $\alpha$ について :eqref[eq-rg-weitzenbock] を $\alpha$ に対し積分すると $0=\int(\norm{\nabla\alpha}^2+g(\mathcal Q_p\alpha,\alpha))\mu_g$ である。第一項は接続 Laplacian の部分積分で得る。被積分関数は非負連続なので各点で零となり、$\alpha$ は平行である。正定値の点があればそこで $\alpha=0$ となり、平行移動で連結成分全体に零が広がる。:ref[thm-rg-harmonic-representative] を使えばコホモロジーの結論を得る。
:::

:::example[定正曲率と積]
定曲率 $k>0$ なら $0<p<n$ で $kp(n-p)>0$ となるので、中間次数の実 de Rham コホモロジーが消える。$S^2\times S^1$ は Ricci 曲率が非負だが、円の $dt$ を引き戻した非零平行 1 形式を持つ。円の接方向では Ricci 値が零のため、正定値の仮定が成り立たない。
:::

閉連結の $\mathrm{Ric}\ge0$ の場合、$b_1\le n$ でもある。調和 1 形式の一点での評価は平行性のため単射であり、値域の次元は $n$ だからである。

:::example[正の曲率があっても中間次数が残る積]
向き付けた $S^2(1)\times S^2(1)$ で第 1 因子の体積形式を引き戻した $\alpha$ は、積の接続と体積形式の平行性により平行である。従って $d\alpha=\delta\alpha=0$ で、第 1 球面に積分すると $4\pi$ なので非零のコホモロジー類を与える。Ricci 曲率は正で 1 形式を消すが、2 形式の $\mathcal Q_2$ はこの例で正定値ではない。
:::

### Ricci 非負の場合に残る 1 形式

コンパクトな $M$ で $\mathrm{Ric}\ge0$ なら、調和 1 形式はすべて平行である。積分した Weitzenböck の式で $\int|\nabla\alpha|^2+\int\mathrm{Ric}(\alpha^\sharp,\alpha^\sharp)=0$ となり、どちらも非負だからである。逆に平行なら $d\alpha=\delta\alpha=0$ なので調和である。平行場は一点での値で決まるので $b_1\le n$ である。さらにある一点で Ricci が正定値なら、その点で全ての平行調和 1 形式が零になり、全体で零である。全点で厳密な下界を置く必要はない。

:::example[次数ごとの消滅を比較する]
$S^2(1)\times S^2(1)$ は $\mathrm{Ric}=g$ なので $b_1=0$ である。一方、第 1 因子の面積形式の引き戻し $\omega_1$ は平行である。積の接続が各因子の体積形式を保つことから従う。同様に $\omega_2$ も平行で、従って双方とも調和 2 形式である。対応する球面因子に沿う積分はそれぞれ $(4\pi,0)$、$(0,4\pi)$ となるので、両者は独立な類を表す。Ricci の条件を Weitzenböck の 2 形式の曲率項にそのまま使えないことが、この例で分かる。
:::

### 曲率作用素から全次数を評価する

$X\wedge Y$ を反対称作用素 $Z\mapsto g(Y,Z)X-g(X,Z)Y$ と同一視し、曲率作用素 $\mathcal R$ は $e_i\wedge e_j$ を $R(e_i,e_j)$ に送るものとする。曲率の対交換から $\mathcal R$ は $\Lambda^2$ 上で自己共役である。反対称作用素 $A$ の形式への誘導作用を $\rho(A)$ と書く。これは各引数に $-A$ を入れた和であり、$\rho(e_i\wedge e_j)=\varepsilon_i\iota_j-\varepsilon_j\iota_i$ である。$\rho(A)$ は反対称である。実際、$e^{tA}$ による誘導変換は各次数の内積を保存し、それを微分すると反対称性が従う。

:eqref[eq-rg-weitzenbock] の次数を保つ部分を取ると
$$
\mathcal Q_p=-\sum_{i<j}\rho(e_i\wedge e_j)\rho(\mathcal R(e_i\wedge e_j)).
$$
右辺は $\Lambda^2$ の正規直交基底の選び方によらない。従って正規直交固有基底 $A_a$ と固有値 $\kappa_a$ で書くと $g(\mathcal Q_p\alpha,\alpha)=\sum_a\kappa_a|\rho(A_a)\alpha|^2$ となる。特に $\mathcal R\ge0$ なら $\mathcal Q_p\ge0$ である。さらに $\mathcal R\ge kI$ なら、曲率作用素が $I$ の場合の計算を使って $\mathcal Q_p\ge k p(n-p)I$ を得る。この下界で $0<p<n$ かつ $k>0$ なら中間次数の調和形式は零になる。ここで各固有値を下から押さえるのは曲率作用素への条件であり、各断面の評価は分解可能な二ベクトルに対する二次形式の値を調べる。

## Laplacian の固有値と Lichnerowicz の評価

:::theorem[Lichnerowicz の評価]{#thm-rg-lichnerowicz}
閉連結な $n\ge2$ 次元多様体で $\mathrm{Ric}\ge(n-1)k\,g$、$k>0$ とする。$-\Delta_g$ の任意の正の固有値は $\lambda\ge nk$ を満たす。
:::

:::proof
$-\Delta_g f=\lambda f$ とする。Bochner の公式を積分すると左辺は零であり、$g(\grad f,\grad\Delta_g f)=-\lambda\norm{\grad f}^2$ である。Hessian の固有値の Cauchy–Schwarz により $\norm{\operatorname{Hess}f}^2\ge(\Delta_g f)^2/n$ だから
$$
0\ge \frac{\lambda^2}{n}\int f^2+((n-1)k-\lambda)\int\norm{\grad f}^2.
$$
部分積分により最後の積分は $\lambda\int f^2$ であり、右辺は $(n-1)\lambda(k-\lambda/n)\int f^2$ である。$\lambda>0,f\ne0$ なので $\lambda\ge nk$ が従う。等号なら非負の差が各点で零になり、$\operatorname{Hess}f=-kf\,g$ である。
:::

:::example[球面で等号]
半径 $a$ の球面に周囲の線形関数 $f(x)=v\cdot x$ を制限する。周囲の Hessian は零なので、球面の Hessian は $v$ と第二基本形式の内積である。外向き法線が $x/a$、$\mathrm{II}=-g\,x/a^2$ なので $\operatorname{Hess}f=-fg/a^2$ である。$\lambda=n/a^2$ であり、下界と一致する。
:::

固有値 $\lambda=0$ の固有関数は連結性から定数である。第一正固有値の存在は、平均零の $H^1$ 関数の Rayleigh 商を最小化して得られる。$L^2$ ノルムを $1$ に正規化した最小化列は $H^1$ に有界であり、Rellich による強収束とエネルギーの弱下半連続性から極限が制約を保ち最小値を取る。Euler 方程式は弱い固有値方程式であり、楕円型正則性により固有関数は滑らかである。値が零なら定数になり平均零と単位ノルムに矛盾する。

:::example[Rayleigh 商を実際に計算する]
半径 $2$ の二次元球面で $f=\cos\theta$ とする。$\int f^2\omega_g=4\cdot 2\pi\int_0^\pi\cos^2\theta\sin\theta\,d\theta=16\pi/3$ である。一方 $|\grad f|^2=\sin^2\theta/4$ なのでエネルギーは $2\pi\int_0^\pi\sin^3\theta\,d\theta=8\pi/3$ である。商は $1/2$ であり、Lichnerowicz の $nk=2/4$ を実現する。
:::

### 離散スペクトルと非負性

閉多様体では各次数の $\Delta_H$ に対し滑らかな固有形式からなる $L^2$ 正規直交基底が取れる。構成には $\Delta_H+1$ の逆を使う。:ref[thm-rg-hodge-decomposition] の $H^1$ 評価から、$Q(u,v)+(u,v)_{L^2}$ は全 $H^1$ で強圧的である。従って $h\in L^2$ を右辺とする弱方程式
$$
Q(u,v)+(u,v)_{L^2}=(h,v)_{L^2}
$$
は一意な $H^1$ 解を持ち、$h\mapsto u$ は $L^2\to H^1$ の有界線形写像である。これを $L^2$ 値とみなした作用素 $G$ は Rellich によりコンパクトである。弱方程式の対称性から $(Gh_1,h_2)=(h_1,Gh_2)$、さらに $(Gh,h)=Q(Gh,Gh)+\norm{Gh}_{L^2}^2\ge0$ である。$Gh=0$ なら全 $H^1$ 試験場で $(h,v)=0$ なので $h=0$ である。

:ref[thm-fa-compact-spectral]{text="コンパクト自己共役作用素のスペクトル定理"} により、$G$ の正の固有値 $\nu_j$ に属する固有形式が正規直交基底をなす。固有方程式を弱方程式に代入すると $\Delta_H e_j=(\nu_j^{-1}-1)e_j$ となり、楕円型正則性を反復して $e_j$ は滑らかである。$Q\ge0$ より $\nu_j^{-1}-1\ge0$ であり、固有値 $0$ はちょうど調和形式に対応する。従って Laplacian の有限の固有値は有限重複度であり、無限個あれば無限大にしか集積しない。

正の固有値 $\lambda$ の $p$ 形式 $\alpha$ は、$d\delta\alpha/\lambda$ と $\delta d\alpha/\lambda$ に直交分解する。両成分は $d,\delta$ と $\Delta_H$ の可換性から同じ固有値に属する。余完全な固有形式 $\alpha$ では $\delta\alpha=0$ なので $\|d\alpha\|_{L^2}^2=(\alpha,\delta d\alpha)=\lambda\|\alpha\|_{L^2}^2$ である。従って $d/\sqrt\lambda$ は余完全な $p$ 形式の固有空間から完全な $(p+1)$ 形式の固有空間への等長写像である。逆写像は $\delta/\sqrt\lambda$ である。完全な固有形式 $\beta$ は $d\beta=0$ だから、$d\delta\beta=\lambda\beta$ を使うと全射性と逆の式が従う。例えば関数の正の固有値は、完全 1 形式の固有値と重複度まで一致する。

:::example[長方形トーラスの低い固有値]
$\R^2/(4\pi\Z\times2\pi\Z)$ の平坦計量では、周期に合わせた $e^{i(mx/2+ny)}$ の固有値は $m^2/4+n^2$ である。零を除く最小値は $(m,n)=(\pm1,0)$ の $1/4$ であり、実固有関数の基底は $\cos(x/2),\sin(x/2)$ である。$dx,dy$ は平行であるため、1 形式ではそれぞれの係数に同じ固有値が現れ、第一正固有値の実重複度は $4$ である。
:::

### 等号から Hessian の方程式へ

$n\ge2$、$\mathrm{Ric}\ge(n-1)kg$、$k>0$ とし、$\lambda=nk$ の固有関数 $f\ne0$ があるとする。Lichnerowicz の証明で使った $\norm{\operatorname{Hess}f}^2\ge(\Delta_g f)^2/n$ は、二つの差が Hessian の跡零部分の二乗であるため、積分で等号なら各点で等号になる。従って
$$
\operatorname{Hess}f=-kf g .
$$
この過剰決定方程式が成立すると、多様体が実際に半径 $1/\sqrt k$ の球面になる。証明のため定数倍で $k=1$ とする。$f$ の積分は零なので正の最大値 $a=f(p)$ を取る。単位動径 $\gamma$ に沿い $f''=-f$、$f(0)=a,f'(0)=0$ だから $f=a\cos t$ である。さらに $D_t\grad f=-fT$ と $\grad f(p)=0$ から $\grad f=-a\sin t\,T$ である。

$0<t<\pi$ で角方向を変えた Jacobi 場 $J$ に対し、変分方向に $\grad f=-a\sin t\,T$ を微分すると
$$
-a\cos t\,J=-a\sin t\,D_tJ .
$$
従って $J(t)=\sin t\,E(t)$ で、初期導関数が $E(0)$ の平行場である。ここでは各動径上での恒等式を変分しているので、指数写像が微分同相であることを先に仮定していない。Jacobi 場が非退化だから指数写像は $t<\pi$ で正則である。同じ点への二つの入力の長さが $\pi$ 未満なら、$f$ の値から長さが一致し、$\grad f$ から終端の速度が一致する。逆向きの測地線の一意性で二つは同じである。従ってこの範囲では正規極座標が使え、計量は $dt^2+\sin^2t\,g_{S^{n-1}}$ となる。

時刻 $\pi$ では全角方向の Jacobi 場が零なので、方向の連結性により終点は一点 $q$ である。Bonnet–Myers と Hopf–Rinow により全点は長さ $\pi$ 以下の最短動径の終点だから、$p,q$ を除くところで上の球面計量の表示が成り立つ。曲率は連続なので両点でも曲率 $1$ である。さらに $p$ での単射半径は $\pi$ であり、長さ $\pi$ の動径の中心以外の任意の点は一意の最短入力を持つ。正曲率空間形の記述で被覆球面から商への全射局所等長を取ると、ここで用いた半径 $\pi$ の開球上では北極からの動径の対応により単射である。また対蹠点だけが $q$ に写る。従って全射は一枚であり、球面と等長である。これが :term[Obata の等号剛性]{en="Obata rigidity" yomi="おばたのとうごうごうせい"} である。

:::example[球面とその対蹠商の第一固有値の差]
単位 $S^n$ の高さ関数は $\Delta_g f=-nf$ なので $\lambda=n$ を実現する。対蹠商 $\R P^n$ に滑らかな関数を持ち上げると対蹠写像で偶でなければならず、高さ関数は降りない。剛性により商では $\lambda_1=n$ となることはない。一方、$f^2-1/(n+1)$ は偶で、$\norm{\grad f}^2=1-f^2$ と積の公式から $\Delta_g(f^2)=2-2(n+1)f^2$ となる。従って $\lambda_1$ は $n$ より大きく $2(n+1)$ 以下である。この例では局所 Ricci 曲率は球面と同じだが、全空間で許される関数が商の同一視で変わる。
:::

:::example[二方向で尺度が異なる計量の Hodge 分解]
$x,y$ がともに周期 $2\pi$ の二次元トーラスに、定数の計量行列 $G=\begin{pmatrix}1 & 0\\0 & 16\end{pmatrix}$ を入れる。ここでは $\phi=2 x - 2 y$、$\alpha=\cos\phi(c_1\,dx+c_2\,dy)$、$c=\begin{pmatrix}1\\0\end{pmatrix}$ とする。このとき $d\phi$ の成分は $\xi=\begin{pmatrix}2\\-2\end{pmatrix}$ であり、コサインの周期境界で値がつながる。

まず $d\phi$ に対応する計量双対ベクトルと二乗長を求める。
$$
\begin{aligned}
\xi^\sharp&=G^{-1}\xi=\begin{pmatrix}2\\- \frac{1}{8}\end{pmatrix}\\
\norm\xi_g^2&=\xi^\top G^{-1}\xi=\frac{17}{4}\\
b=c(\xi^\sharp)&=2
\end{aligned}
$$
である。定数計量での Levi-Civita 記号は零であり、余微分は $-\sum g^{jk}\partial_j\alpha_k$ になる。従って $\delta\alpha=2\sin\phi$ である。一方、$d(\cos\phi)=-\sin\phi\,d\phi$ だから外微分は $-2\sin\phi\,dx\wedge dy$ になる。この二つの微分の係数は、波数と係数ベクトルの組み方によって別々に消える。

Laplacian は定数座標の各成分に $-\sum_{j,k}g^{jk}\partial_j\partial_k$ を作用させる。二階微分で $-\xi_j\xi_k$ を掛けるので $\Delta_H\alpha=\norm\xi_g^2\alpha=\frac{17}{4}\alpha$ である。調和成分は平均によって決まり、ここでは $\cos\phi$ の平均が零なので消える。実際、少なくとも一方の波数が非零整数だから、該当する変数を一周積分すると零になる。よって Green 作用素はこの周波数では $G_1\alpha=\alpha/\frac{17}{4}$ である。

完全成分は $d\delta G_1\alpha$ なので、その係数ベクトルは $c_{\rm ex}=b\xi/\norm\xi_g^2$ である。残りの係数は $c_{\rm co}=c-c_{\rm ex}$ である。指定した数値では
$$
\begin{aligned}
c_{\rm ex}&=\begin{pmatrix}\frac{16}{17}\\- \frac{16}{17}\end{pmatrix}\\
c_{\rm co}&=\begin{pmatrix}\frac{1}{17}\\\frac{16}{17}\end{pmatrix}
\end{aligned}
$$
であり、$c_{\rm co}\cdot \xi^\sharp=0$ だから余閉である。調和成分はすでに除かれているため、残りは Hodge 分解の余完全成分と一致する。このように「波数方向への射影」は余接空間の計量 $G^{-1}$ による射影である。Euclid の成分の積で射影すると、計量に交差項がある場合には結果が異なる。

体積は $(2\pi)^2\sqrt{\det G}=16 \pi^{2}$ であり、$\cos^2\phi$ の平均は $1/2$ である。従って完全成分と余完全成分の積分ノルムの二乗はそれぞれ $\frac{128 \pi^{2}}{17}$ と $\frac{8 \pi^{2}}{17}$ となる。両者の和は元の形式の積分ノルムの二乗 $8 \pi^{2}$ である。この例では各係数が $\cos\phi$ に比例するので点ごとの直交性まで成立する。複数の周波数を重ねるときは、積分で周波数の直交性を使う。
:::

:::example[単位球面の次数 3 の高さ多項式]
単位 $S^n$、$n=2$ の最後の座標を高さ $z$ とする。$\operatorname{Hess}z=-zg$、$\norm{\grad z}^2=1-z^2$ である。前者は球面の第二基本形式から得られ、後者は周囲の単位高さ方向から球面法線の成分 $z$ を引いたときの二乗長である。関数 $h(z)$ に Laplacian の連鎖則を使うと
$$
\Delta_g h=(1-z^2)h''-nz h'
$$
となる。最高次が $z^j$ の多項式に作用させると、最高次の係数は $-j(j+n-1)$、それより二つ低い次数の係数は $j(j-1)$ になる。従って固有値の候補は $\lambda=j(j+n-1)$ である。次の次数の係数 $a_{j-2}$ は、$(\Delta_g+\lambda)h$ の $z^{j-2}$ の項を零と置き
$$
a_{j-2}=-\frac{j(j-1)}{\lambda-(j-2)(j+n-3)}
$$
と決める。さらに $z^{j-4}$ 以降も、直前の係数を使って同じ式で消す。分母は $\lambda$ とそれより低い次数の候補固有値の差なので正である。

指定した次元と次数では $\lambda=12$、係数を代入した多項式は
$$
h(z)=z^{3} - \frac{3 z}{5},\qquad \Delta_g h=- 12 z^{3} + \frac{36 z}{5}
$$
である。両辺を項ごとに比べると $-\Delta_g h=\lambda h$ となる。$h$ は周囲の座標関数の多項式を制限したものなので、極でも滑らかである。非零固有値の固有関数は定数と直交する。実際、式を積分すると左辺は発散の積分で零になり、右辺は $\lambda\int h$ なので $\int h=0$ である。

Lichnerowicz の下界はこの球面で $n=2$ であり、一次多項式 $z$ が実現する。ここで求めた多項式は次数が大きいため、より大きい固有値を持つ。Hodge–Laplacian の符号ではこれらの固有値は正となり、$dh$ も同じ固有値を持つ。$d$ が作用素と可換という性質により、関数の固有関数から正の固有値を持つ完全 1 形式を作れる。
:::

:::example[単位球面の次数 3 の高さ多項式]
単位 $S^n$、$n=3$ の最後の座標を高さ $z$ とする。$\operatorname{Hess}z=-zg$、$\norm{\grad z}^2=1-z^2$ である。前者は球面の第二基本形式から得られ、後者は周囲の単位高さ方向から球面法線の成分 $z$ を引いたときの二乗長である。関数 $h(z)$ に Laplacian の連鎖則を使うと
$$
\Delta_g h=(1-z^2)h''-nz h'
$$
となる。最高次が $z^j$ の多項式に作用させると、最高次の係数は $-j(j+n-1)$、それより二つ低い次数の係数は $j(j-1)$ になる。従って固有値の候補は $\lambda=j(j+n-1)$ である。次の次数の係数 $a_{j-2}$ は、$(\Delta_g+\lambda)h$ の $z^{j-2}$ の項を零と置き
$$
a_{j-2}=-\frac{j(j-1)}{\lambda-(j-2)(j+n-3)}
$$
と決める。さらに $z^{j-4}$ 以降も、直前の係数を使って同じ式で消す。分母は $\lambda$ とそれより低い次数の候補固有値の差なので正である。

指定した次元と次数では $\lambda=15$、係数を代入した多項式は
$$
h(z)=z^{3} - \frac{z}{2},\qquad \Delta_g h=- 15 z^{3} + \frac{15 z}{2}
$$
である。両辺を項ごとに比べると $-\Delta_g h=\lambda h$ となる。$h$ は周囲の座標関数の多項式を制限したものなので、極でも滑らかである。非零固有値の固有関数は定数と直交する。実際、式を積分すると左辺は発散の積分で零になり、右辺は $\lambda\int h$ なので $\int h=0$ である。

Lichnerowicz の下界はこの球面で $n=3$ であり、一次多項式 $z$ が実現する。ここで求めた多項式は次数が大きいため、より大きい固有値を持つ。Hodge–Laplacian の符号ではこれらの固有値は正となり、$dh$ も同じ固有値を持つ。$d$ が作用素と可換という性質により、関数の固有関数から正の固有値を持つ完全 1 形式を作れる。
:::

:::example[単位球面の次数 4 の高さ多項式]
単位 $S^n$、$n=4$ の最後の座標を高さ $z$ とする。$\operatorname{Hess}z=-zg$、$\norm{\grad z}^2=1-z^2$ である。前者は球面の第二基本形式から得られ、後者は周囲の単位高さ方向から球面法線の成分 $z$ を引いたときの二乗長である。関数 $h(z)$ に Laplacian の連鎖則を使うと
$$
\Delta_g h=(1-z^2)h''-nz h'
$$
となる。最高次が $z^j$ の多項式に作用させると、最高次の係数は $-j(j+n-1)$、それより二つ低い次数の係数は $j(j-1)$ になる。従って固有値の候補は $\lambda=j(j+n-1)$ である。次の次数の係数 $a_{j-2}$ は、$(\Delta_g+\lambda)h$ の $z^{j-2}$ の項を零と置き
$$
a_{j-2}=-\frac{j(j-1)}{\lambda-(j-2)(j+n-3)}
$$
と決める。さらに $z^{j-4}$ 以降も、直前の係数を使って同じ式で消す。分母は $\lambda$ とそれより低い次数の候補固有値の差なので正である。

指定した次元と次数では $\lambda=28$、係数を代入した多項式は
$$
h(z)=z^{4} - \frac{2 z^{2}}{3} + \frac{1}{21},\qquad \Delta_g h=- 28 z^{4} + \frac{56 z^{2}}{3} - \frac{4}{3}
$$
である。両辺を項ごとに比べると $-\Delta_g h=\lambda h$ となる。$h$ は周囲の座標関数の多項式を制限したものなので、極でも滑らかである。非零固有値の固有関数は定数と直交する。実際、式を積分すると左辺は発散の積分で零になり、右辺は $\lambda\int h$ なので $\int h=0$ である。

Lichnerowicz の下界はこの球面で $n=4$ であり、一次多項式 $z$ が実現する。ここで求めた多項式は次数が大きいため、より大きい固有値を持つ。Hodge–Laplacian の符号ではこれらの固有値は正となり、$dh$ も同じ固有値を持つ。$d$ が作用素と可換という性質により、関数の固有関数から正の固有値を持つ完全 1 形式を作れる。
:::

### 斉次多項式から球面調和関数を作る

Euclid 空間 $\R^{n+1}$ の斉次多項式 $P$ が次数 $j$ で、$\Delta_{\rm E}P=0$ とする。単位球面への制限を $Y$ とすれば $P(ru)=r^jY(u)$ である。極座標の計量 $dr^2+r^2g_{S^n}$ で Laplacian を計算すると
$$
\Delta_{\rm E}=\partial_r^2+\frac nr\partial_r+\frac1{r^2}\Delta_{S^n}
$$
となる。実際、密度の動径因子は $r^n$ なので動径項は $r^{-n}\partial_r(r^n\partial_r)$、角方向の逆計量の倍率が $r^{-2}$ である。$r^jY$ を代入すると
$$
0=\Delta_{\rm E}P=r^{j-2}\bigl(j(j+n-1)Y+\Delta_{S^n}Y\bigr).
$$
従って $Y$ の固有値は $j(j+n-1)$ である。この構成では一つの高さの多項式に限らず、角度や横方向の座標を含む多項式も使える。

この方法で多項式の制限を全て固有関数の線形結合に分解できることを示す。まず次数 $\ell$ の調和斉次多項式 $H$ に対し、$\Delta(r^{2s}H)=2s(2\ell+2s+n-1)r^{2s-2}H$ である。上の極座標の式に $r^{\ell+2s}Y$ を入れ、球面固有値 $\ell(\ell+n-1)$ を引けば係数が得られる。次に次数 $j$ の任意の斉次多項式 $P$ を取る。$\Delta P$ は次数 $j-2$ であるから、次数について帰納的に $\sum_s r^{2s}H_{j-2-2s}$ の形に分解できる。各項の Laplacian の逆像は $r^{2(s+1)}H_{j-2-2s}$ を上の正の係数で割って作れる。それらの和を $P$ から引けば残りは調和である。次数 $0,1$ はそのまま調和なので帰納法が始まる。球面上では $r=1$ となるため、どの多項式の制限も調和斉次多項式の制限の和になる。

この固有関数族は $L^2(S^n)$ で全体を張る。実多項式の制限は定数を含み、各点を座標関数で分離するので、Stone–Weierstrass の定理により連続関数を一様に近似する。連続関数は $L^2$ で稠密である。別の固有値の固有関数があるなら対称性から全てのこの族に直交してしまい、稠密性から零となる。従って球面の関数の固有値は $j(j+n-1)$ で尽くされる。各 $j$ の空間が非零であることは、例えば $\Re(x_1+ix_2)^j$ が $\partial_1^2+\partial_2^2$ で零となるため確かめられる。

:::example[射影空間の第一正固有値]
球面上の斉次次数 $j$ の固有関数は、対蹠変換で $(-1)^j$ 倍になる。前の稠密性の議論を偶関数に限るには、近似多項式 $P$ を $(P(x)+P(-x))/2$ に置き換えればよい。従って $\R P^n$ の固有値は偶数の次数から得られ、第一正固有値は $2(n+1)$ である。例えば $n=3$ では $8$ で、二次調和斉次多項式 $x_1^2-x_2^2$ の制限が実現する。単位球面の第一正固有値は $3$ なので、Ricci 曲率の局所的な下界を変えずに、商の条件だけで第一固有値が増える。
:::

## 演習問題

::::exercise{level=basic}
正の Euclid 三次元共標構で $\star(dx\wedge dz)$ を求めよ。

:::solution
$(dx\wedge dz)\wedge dy=-dx\wedge dy\wedge dz$ なので、元の形式の二乗ノルム $1$ が出るようにすると $\star(dx\wedge dz)=-dy$ である。
:::
::::

::::exercise{level=basic}
$n=6,p=3$ で $\star^2$ の符号を求めよ。

:::solution
$p(n-p)=3\cdot3=9$ なので符号は $(-1)^9=-1$ である。
:::
::::

::::exercise{level=basic}
平坦な周長 $2\pi$ の円で $\cos(3t)\,dt$ の Hodge 固有値を求めよ。

:::solution
平行な $dt$ の係数に $-\partial_t^2$ が作用する。$-\partial_t^2\cos(3t)=9\cos(3t)$ だから固有値は $9$ である。
:::
::::

::::exercise{level=basic}
曲率 $2$ の五次元定曲率模型で $\mathcal Q_2$ の固有値を求めよ。

:::solution
式 :eqref[eq-rg-weitzenbock] の定曲率の場合に $k=2,p=2,n=5$ を入れると $2\cdot2\cdot3=12$ である。
:::
::::

::::exercise{level=advanced}
周長 $2\pi$ の円で $\alpha=(3+2\sin t)\,dt$ を Hodge 分解せよ。

:::solution
係数の平均は $3$ なので調和成分は $3\,dt$ である。残りは $d(-2\cos t)$ であり、$2$ 形式がないため余完全成分は零である。
:::
::::

::::exercise{level=advanced}
閉多様体上の調和関数 $f$ の微分が零であることを随伴性から示せ。

:::solution
$\delta f=0$ なので $(\Delta_Hf,f)=(\delta df,f)=(df,df)$ である。調和性により左辺が零で、右辺が非負連続なノルムの二乗の積分なので $df=0$ である。
:::
::::

::::exercise{level=advanced}
定曲率 $1$ の閉向き付き四次元球面で $b_2$ を求めよ。

:::solution
$\mathcal Q_2=1\cdot2(4-2)\id=4\id$ である。:ref[thm-rg-bochner-vanishing] により調和 $2$ 形式は零なので $b_2=0$ である。
:::
::::

::::exercise{level=advanced}
三次元の閉連結多様体で $\mathrm{Ric}\ge6g$ なら第一正固有値の下界はいくつか。

:::solution
$(n-1)k=6$ から $k=3$ である。:ref[thm-rg-lichnerowicz] により $\lambda_1\ge nk=9$ である。
:::
::::

## 確認テスト

::::quiz
閉多様体上の調和形式を特徴づける条件を選べ。

- [x] 閉かつ余閉
- [ ] 閉であるだけ

:::explanation
積分恒等式 $(\Delta_H\alpha,\alpha)=\norm{d\alpha}^2+\norm{\delta\alpha}^2$ により、核の元は両方零である。逆に両方零なら $\Delta_H\alpha=0$ である。
:::
::::

::::quiz{answer="2" tol="1e-6"}
平坦二次元トーラスの調和 1 形式の次元を求めよ。

:::explanation
平行標構で各係数は調和関数で、連結コンパクト性から定数である。従って $dx,dy$ が基底をなし次元は $2$ である。
:::
::::

::::quiz{answer="4" tol="1e-6"}
単位四次元球面で Lichnerowicz が与える第一正固有値の下界はいくつか。

:::explanation
$\mathrm{Ric}=(4-1)g$ だから $k=1$ であり、$nk=4$ となる。球面の高さ関数はその固有値を実現する。
:::
::::

::::quiz
Hodge 分解の結論を選べ。

- [x] de Rham 類ごとに調和形式が一意に決まる。
- [ ] 調和代表は計量を変えても形式として必ず同じである。

:::explanation
閉形式から完全成分を引けば調和代表になり、完全な調和形式は零である。余微分は計量に依存する。例えば円の非定数の密度を持つ計量では、固定周期の調和 $1$ 形式は弧長形式に比例し、一定密度のものとは異なる。
:::
::::

---
title: 微分形式
description: 交代テンソルの外積から外微分を構成し、Cartan の公式とベクトル解析に結び付ける。
notation:
  - symbol: '\Omega^k(M)'
    meaning: 多様体 M 上の滑らかな k 形式全体
  - symbol: '\mathrm d'
    meaning: 外微分
  - symbol: '\iota_X'
    meaning: ベクトル場 X を最初のスロットに入れる内部積
  - symbol: '\wedge'
    meaning: 交代テンソルの外積
---

## この章の目標

- 交代テンソルを外積の基底で表示し、外積を計算する。
- 外微分を構成して一意性と引き戻しとの交換を証明する。
- 内部積と Lie 微分の関係を Cartan の公式で表す。
- 三次元の勾配・回転・発散を外微分の各次数として理解する。

## 交代テンソルと外積

面積には、二つの辺を交換すると符号が逆になるという性質がある。この性質を任意個の入力に広げると、平行に重なった方向を測ったとき零になる多重線形写像が得られる。

:::definition[交代テンソル]{#def-mfd-alternating-tensor}
実ベクトル空間 $V$ 上の $k$ 重線形写像 $\alpha:V^k\to\R$ が二つの入力を交換するたびに符号を変えるとき :term[交代テンソル]{en="alternating tensor" yomi="こうたいてんそる"} という。全体を $\Lambda^k V^*$ と書き、$k=0$ では $\Lambda^0 V^*=\R$ とする。
:::

同じ入力が二つあれば、その交換で値が負の自分に等しくなるため零になる。入力が一次従属なら、一つを他の線形結合で表し、多重線形性と重複による消失を使うと値は零である。従って $k>\dim V$ なら $\Lambda^k V^*=0$ である。

:::definition[外積]{#def-mfd-wedge-product}
$\alpha\in\Lambda^k V^*$、$\beta\in\Lambda^\ell V^*$ の :term[外積]{en="wedge product" yomi="がいせき"} を
$$
(\alpha\wedge\beta)(v_1,\ldots,v_{k+\ell})
=\frac1{k!\ell!}\sum_{\sigma\in S_{k+\ell}}\sgn(\sigma)
\alpha(v_{\sigma(1)},\ldots,v_{\sigma(k)})
\beta(v_{\sigma(k+1)},\ldots,v_{\sigma(k+\ell)})
\label{eq-mfd-wedge-permutation}
$$
と定める。次数 $0$ の外積はスカラー倍になる。
:::

二入力を入れ替えると、和の置換が符号を反転させながら並べ替わるため、外積は交代的である。各ブロック内を並べ替えた項は、もとの形式の符号変化と置換の符号変化が相殺する。従ってブロック内で添字を昇順にした置換に限って、係数 $1/(k!\ell!)$ を外して足せる。この表示はシャッフル表示と呼ばれる。

:::proposition[外積の結合律と次数付き交換律]{#prop-mfd-wedge-rules}
外積は双線形で結合的である。また $\alpha\wedge\beta=(-1)^{k\ell}\beta\wedge\alpha$ である。
:::

:::proof
双線形性は式 :eqref[eq-mfd-wedge-permutation] の各項が双線形であることから従う。三因子の積をシャッフル表示で展開すると、どちらの括弧でも入力添字を三つの昇順のブロックに分け、分けた順序に対応する置換の符号を付ける。同じ三つのブロックの分け方は一度ずつ現れる。符号は置換の合成で積になるため一致する。交換律は二つのブロックを入れ替える際に、前の $k$ 個が後ろの $\ell$ 個を一つずつ横切るので交換が $k\ell$ 回起きることから従う。
:::

双対基底 $(\varepsilon^i)$ から作る $\varepsilon^{i_1}\wedge\cdots\wedge\varepsilon^{i_k}$（$i_1<\cdots<i_k$）の評価は $\det(\varepsilon^{i_a}(v_b))$ である。これは一形式の外積をシャッフル表示で順に展開すると行列式の置換和になるからである。昇順の基底入力の組に対してこの値は同じ添字集合なら $1$、別なら零である。交代テンソルは昇順の基底入力での値から、符号と多重線形性で全入力に復元できるため、これらの外積は基底をなす。従って次元は $\binom{\dim V}{k}$ である。

:::example[向きが異なる二組の値]
$\R^3$ で $\alpha=\mathrm dx+2\mathrm dy$、$\beta=3\mathrm dy+\mathrm dz$ とする。外積は $3\mathrm dx\wedge\mathrm dy+\mathrm dx\wedge\mathrm dz+2\mathrm dy\wedge\mathrm dz$ である。入力 $(\partial_y,\partial_x)$ では第 1 項が $-3$ を、他は零を返す。この値は $\alpha(\partial_y)\beta(\partial_x)-\alpha(\partial_x)\beta(\partial_y)$ からも得られる。
:::

### 交代化と成分の取り出し

任意の $k$ 重線形写像 $T$ から
$$
\operatorname{Alt}(T)(v_1,\ldots,v_k)
=\frac1{k!}\sum_{\sigma\in S_k}\sgn(\sigma)
T(v_{\sigma(1)},\ldots,v_{\sigma(k)})
$$
を作れる。入力を入れ替えると和の添字を入れ替えた上で符号が反転するため、結果は交代的である。$T$ がすでに交代的なら和の各項が $T$ 自身に等しいので、$\operatorname{Alt}(T)=T$ である。従って交代化を二度適用しても一度と同じになる。外積はこの交代化で
$$
\alpha\wedge\beta=\frac{(k+\ell)!}{k!\ell!}\operatorname{Alt}(\alpha\otimes\beta)
$$
と表せる。係数が必要なのは、例えば $\varepsilon^1\wedge\varepsilon^2$ が $(e_1,e_2)$ で $1$ を返すよう正規化するためである。この場合 $\operatorname{Alt}(\varepsilon^1\otimes\varepsilon^2)$ の値は $1/2$ である。

例えば $T=\mathrm dx\otimes\mathrm dy+2\,\mathrm dy\otimes\mathrm dx$ の交代部は $-\frac12\mathrm dx\wedge\mathrm dy$ となる。実際、$T(\partial_x,\partial_y)=1$、$T(\partial_y,\partial_x)=2$ なので、その差の半分が係数である。これに対し $\mathrm dx\wedge\mathrm dy$ の反対称成分を順序付きの全ての組で足すと、同じ係数を二度数える。一般に成分を昇順 $I$ だけで書くときは係数が $\omega(\partial_{i_1},\ldots,\partial_{i_k})$ に一致する。順序付き添字の全てを使って $\sum_{i_1,\ldots,i_k}\omega_{i_1\ldots i_k}\mathrm dx^{i_1}\wedge\cdots\wedge\mathrm dx^{i_k}$ と書くなら、各 $I$ の並べ替えが $k!$ 回同じ項を作るため、全体を $k!$ で割る必要がある。

### 分解できる形式とできない形式

$x$ ごとに一形式の外積で書けるかどうかは、次数だけでは決まらない。例えば $\R^3$ の二形式 $\omega=a\,\mathrm dx\wedge\mathrm dy+b\,\mathrm dx\wedge\mathrm dz+c\,\mathrm dy\wedge\mathrm dz$ で $a\ne0$ なら
$$
\omega=a(\mathrm dx-(c/a)\mathrm dz)\wedge(\mathrm dy+(b/a)\mathrm dz)
$$
となる。積を展開して $\mathrm dz\wedge\mathrm dy$ の符号を直せば確認できる。$a=0$ で他の係数が非零なら、座標の順序を付け替えて同じことができる。従って三次元では各非零二形式を二つの一形式の外積にできる。一方、四次元の $\omega=\mathrm dx^1\wedge\mathrm dx^2+\mathrm dx^3\wedge\mathrm dx^4$ はこのように分解できない。$\omega=\alpha\wedge\beta$ なら $\omega\wedge\omega$ に $\alpha$ が二回現れ、交換して隣り合わせれば零となるが、この $\omega$ の平方は非零である。次数の基底の表現と一つの外積への分解は異なる。

### 二形式の階数と外積の冪

二形式 $\beta$ に対して線形写像 $v\mapsto\beta(v,\cdot)$ の階数を考える。階数は、どれだけの方向が二形式の評価に関与するかを表す。例えば $\mathrm dx\wedge\mathrm dy$ の核は $x,y$ 以外の座標方向である。

:::proposition[交代二形式の標準形]{#prop-mfd-alternating-normal-form}
有限次元実ベクトル空間 $V$ 上の二形式には、ある双対基底 $\varepsilon^1,\ldots,\varepsilon^n$ で
$$
\beta=\varepsilon^1\wedge\varepsilon^2+\cdots+\varepsilon^{2r-1}\wedge\varepsilon^{2r}
$$
と表せるものがある。ここで $2r$ は二形式の階数である。特に零でない二形式では、$\beta^r$ は非零で $\beta^{r+1}=0$ である。
:::

:::proof
$\beta\ne0$ なら $\beta(v,w)\ne0$ となる組を選び、$w$ を倍率で調整して $\beta(v,w)=1$ とする。任意の $z$ は、$a=\beta(z,w),b=\beta(v,z)$ として $z=av+bw+z_0$ と分解できる。実際、その差 $z_0$ は $v,w$ の両方と評価すると零になる。平面上の制限は非退化なのでこの分解は直和である。残った部分で同じ手順を繰り返し、形式が零になる部分には任意の基底を取る。これで形式の行列が正則な二次ブロックと零ブロックに分かれ、階数は $2r$ となる。外積の冪では同じ二つの添字を使う項が消える。一方、異なる二次因子は交換しても符号が変わらないから $\beta^r=r!\,\varepsilon^1\wedge\cdots\wedge\varepsilon^{2r}$ である。次の冪では全項に重複がある。
:::

:::example[階数を行列の消去でなく外積から計算する]
$\R^4$ で $\beta=(\mathrm dx+\mathrm dy)\wedge\mathrm dz+2\,\mathrm dx\wedge\mathrm dw$ とする。平方では二つの項の交差だけが残り、$\beta^2=4\,\mathrm dy\wedge\mathrm dz\wedge\mathrm dx\wedge\mathrm dw=4\,\mathrm dx\wedge\mathrm dy\wedge\mathrm dz\wedge\mathrm dw$ である。従って階数は $4$ である。これに対し $\mathrm dw$ を $\mathrm dy$ に替えた場合は平方が零であり、形式自体は非零なので階数は $2$ になる。
:::

### 外積を保つ写像の作り方

外積代数 $\Lambda^\bullet V^*=\bigoplus_{k=0}^{\dim V}\Lambda^k V^*$ は、1 形式を生成元とし、各生成元の平方を零とする関係から作れる。正確には、単位元を持つ実結合代数 $A$ と線形写像 $L:V^*\to A$ があり、すべての $\alpha$ について $L(\alpha)^2=0$ なら、$L$ を拡張する単位的代数準同型 $\Lambda^\bullet V^*\to A$ がただ一つある。

実際、$L(\alpha+\beta)^2=0$ を展開すると $L(\alpha)L(\beta)=-L(\beta)L(\alpha)$ である。双対基底の昇順の外積に対し、同じ順の $L(\varepsilon^{i_1})\cdots L(\varepsilon^{i_k})$ を対応させ、空積は $1_A$ に送る。この対応を基底から線形に延ばす。二つの基底外積を掛けると、添字が重なる場合は並べ替えて平方が零になり、重ならない場合は並べ替えの符号が付く。$A$ 側も同じ規則を満たすので積が保たれる。基底の各元は 1 形式の積であるため値は強制され、一意性も従う。

例えば余接空間の線形写像 $B:W^*\to V^*$ は、$B(\alpha)\wedge B(\alpha)=0$ を満たすので、外積代数の準同型に拡張できる。$B$ が微分の双対なら、この拡張が後で定義する形式の引き戻しの各点での作用である。

## 微分形式と外積

:::definition[微分形式]{#def-mfd-differential-form}
$\Lambda^k T^*M$ の滑らかな切断を :term[微分形式]{en="differential form" yomi="びぶんけいしき"} と呼び、$k$ 次の全体を $\Omega^k(M)$ と書く。各座標で
$$
\omega=\sum_{i_1<\cdots<i_k}a_{i_1\ldots i_k}\,\mathrm dx^{i_1}\wedge\cdots\wedge\mathrm dx^{i_k}
$$
と表され、係数が滑らかであることが滑らかな切断の条件である。次数 $0$ では $\Omega^0(M)=C^\infty(M)$ である。
:::

余接束の各ファイバーの変換を交代テンソルに適用すると、基底の外積は別の基底の外積の線形結合になる。係数はヤコビ行列の小行列式で、滑らかな正則変換をなす。従って :ref[thm-mfd-bundle-gluing] により $\Lambda^k T^*M$ はベクトル束である。こうして外積は各ファイバーで取れ、滑らかな係数の有限和と積からできるため滑らかな形式になる。

:::example[平方が零にならない偶数次数の形式]
$\R^4$ の $\omega=\mathrm dx^1\wedge\mathrm dx^2+\mathrm dx^3\wedge\mathrm dx^4$ では
$$
\omega\wedge\omega=2\,\mathrm dx^1\wedge\mathrm dx^2\wedge\mathrm dx^3\wedge\mathrm dx^4.
$$
各項の自己積は添字が重複して零だが、二つの交差項は交換の符号が $(-1)^4=1$ なので加わる。一方、奇数次数の形式の自己積は交換律により自分の負となり零である。
:::

:::example[三次元の基底の変換]
$(u,v,w)=(x+y,y+z,z+x)$ と置く。和を展開すると $\mathrm du\wedge\mathrm dv\wedge\mathrm dw=2\,\mathrm dx\wedge\mathrm dy\wedge\mathrm dz$ である。実際、三つの異なる添字を選ぶ項だけが残り、それらは恒等順序と巡回順序で符号がともに正である。係数 $2$ はヤコビ行列式に一致する。
:::

### 座標の外積基底の小行列式

座標変換のヤコビ行列 $J_{ai}=\partial_i y^a$ に対して、$\mathrm dy^{a_1}\wedge\cdots\wedge\mathrm dy^{a_k}$ を $x$ 座標で展開した $\mathrm dx^{i_1}\wedge\cdots\wedge\mathrm dx^{i_k}$ の係数は $\det(J_{a_s i_t})$ である。実際、形式を $\partial_{x^{i_1}},\ldots,\partial_{x^{i_k}}$ に評価すれば、各 $\mathrm dy^{a_s}$ の値がその行列の成分になり、外積の評価が行列式になる。従って成分の変換では行列の $k$ 行と $k$ 列を選ぶ。$k=1$ では成分そのもの、$k=n$ では全ヤコビ行列式だけを使う。

例えば $(u,v,w)=(x+y,x-y,z+x)$ とすると $\mathrm du\wedge\mathrm dw$ の $(\mathrm dx\wedge\mathrm dy,\mathrm dx\wedge\mathrm dz,\mathrm dy\wedge\mathrm dz)$ の順での係数は $(-1,1,1)$ である。これは $\mathrm du=\mathrm dx+\mathrm dy$、$\mathrm dw=\mathrm dz+\mathrm dx$ を展開しても得られる。最高次の $\mathrm du\wedge\mathrm dv\wedge\mathrm dw$ の係数は $-2$ で、変換後の座標の並びは元の並びと逆の向きを持つ。この符号は積分の変数変換で重要になる。

共変テンソルの引き戻しを形式に適用すると、交代性は保たれる。曲面への制限では、曲面の包含の微分を各入力に適用してから形式を評価する。以下ではこの操作を具体的に計算する。

:::example[曲面に制限すると消える二形式]
$F(u,v)=(u,v,u+2v)$ を取る。$\omega=(\mathrm dz-\mathrm dx-2\mathrm dy)\wedge\mathrm dx$ とすると、周囲では非零の二形式であるが、$F^*\mathrm dz-\mathrm d u-2\mathrm dv=0$ なので引き戻しは零である。曲面を動ける速度に限れば、最初の一形式はどの入力にも零を返す。従って形式の係数を単に周囲の座標に代入するだけでは制限にならず、速度の変換によって各 $\mathrm dx$ も引き戻す必要がある。
:::

:::example[基底の向きが変わる二形式]
平面で $(u,v)=(y,2x+y)$ とする。$\mathrm du\wedge\mathrm dv=\mathrm dy\wedge(2\,\mathrm dx+\mathrm dy)=-2\,\mathrm dx\wedge\mathrm dy$ なので、$\mathrm dx\wedge\mathrm dy=-\tfrac12\,\mathrm du\wedge\mathrm dv$ である。旧座標の入力 $\partial_x,\partial_y$ を新座標に変えると $(0,2),(1,1)$ で、その行列式は $-2$ である。新しい係数と掛けると $1$ となり、旧座標での評価に一致する。
:::

:::example[曲線で消える形式と曲面で残る形式]
座標 $x,y,z$ を持つ $\R^3$ で $\omega=\mathrm dx\wedge\mathrm dz+2\,\mathrm dy\wedge\mathrm dz$ とする。曲線に沿う二つの接ベクトルは比例するので、曲線上に引き戻した二形式は零である。一方、曲面 $(u,v)\mapsto(u,v,uv)$ では $\mathrm dz=v\,\mathrm du+u\,\mathrm dv$ となるので、引き戻しは $(u-2v)\,\mathrm du\wedge\mathrm dv$ である。$(u,v)=(3,1)$ では係数が $1$ である。曲面上でも $u=2v$ に限って零になるが、これは写像の階数が落ちることとは異なる。最初の二成分が恒等なので、曲面の写像の微分は全点で単射である。
:::

## 外微分の存在と一意性

一形式が勾配から来る場合、係数を別の座標方向で微分して比較すると混合微分が一致する。その差を集めて二形式にする操作を、任意の次数に拡張する。

:::theorem[外微分]{#thm-mfd-exterior-derivative}
すべての次数で定義される実線形写像 $\mathrm d:\Omega^k(M)\to\Omega^{k+1}(M)$ で、次を満たすものが一意に存在する。
1. 関数に対して従来の微分 $\mathrm df$ に一致する。
2. $\alpha$ が $k$ 次なら $\mathrm d(\alpha\wedge\beta)=\mathrm d\alpha\wedge\beta+(-1)^k\alpha\wedge\mathrm d\beta$ である。
3. $\mathrm d^2=0$ である。

この写像を :term[外微分]{en="exterior derivative" yomi="がいびぶん"} という。局所座標では
$$
\mathrm d\left(\sum_I a_I\,\mathrm dx^I\right)
=\sum_I\sum_j\partial_j a_I\,\mathrm dx^j\wedge\mathrm dx^I
\label{eq-mfd-exterior-local}
$$
である。ここで $I=(i_1<\cdots<i_k)$、$\mathrm dx^I=\mathrm dx^{i_1}\wedge\cdots\wedge\mathrm dx^{i_k}$ と書く。
:::

:::proof
まず座標内で式 :eqref[eq-mfd-exterior-local] を使って定める。係数の積の微分法則を適用し、後の係数を微分した $\mathrm dx^j$ を前の $k$ 個の因子の後ろに動かすと $(-1)^k$ が出る。これで積公式を得る。また二度適用した式では $\partial_i\partial_j a_I$ が $i,j$ に関して対称であり、$\mathrm dx^i\wedge\mathrm dx^j$ が反対称なので全項が消える。この局所演算は関数の微分と一致する。

別の座標 $y$ で $\alpha=\sum a_J\mathrm dy^J$ と表す。既に作った $x$ での演算は関数の微分に一致するので、$\mathrm dy^i$ に外微分を適用すると $\mathrm d^2y^i=0$ になる。従って積公式で計算すると $\sum_J\mathrm da_J\wedge\mathrm dy^J$ となり、$y$ の式と一致する。ゆえに局所式は貼り合わさる。

一意性のために、公理を満たす任意の演算が局所的であることを確かめる。$\alpha$ が $p$ の近傍で零なら、その近傍に台を持ち $p$ 付近で $1$ の $\chi$ について $\chi\alpha=0$ である。積公式から $p$ で $\mathrm d(\chi\alpha)=\mathrm d\alpha$ となるため、$\mathrm d\alpha$ は $p$ で零である。局所の関数や形式に隆起関数を掛けて大域に延長できるので、座標表示の各項に公理を適用できる。$\mathrm d(\mathrm dx^i)=0$ と積公式により残るのは式 :eqref[eq-mfd-exterior-local] だけである。
:::

:::example[一形式と二形式の計算]
$\R^2$ で $\alpha=x^2y\,\mathrm dx+xy^2\,\mathrm dy$ とすると、
$$
\mathrm d\alpha=(-x^2+y^2)\,\mathrm dx\wedge\mathrm dy.
$$
第 1 係数を $y$ で微分すると $x^2\mathrm dy\wedge\mathrm dx$ となるため負号が付く。第 2 係数では $x$ の微分から $y^2\mathrm dx\wedge\mathrm dy$ が出る。例えば $(x,y)=(2,1)$ で係数は $-3$ である。一方 $\omega=xz\,\mathrm dx\wedge\mathrm dy$ を $\R^3$ で微分すると $x\,\mathrm dx\wedge\mathrm dy\wedge\mathrm dz$ となる。$x$ に関する微分の項は $\mathrm dx$ が重複して消え、$\mathrm dz$ を二つ動かしても符号は正である。
:::

負次数では $\Omega^k(M)=0$ とする。$\mathrm d\alpha=0$ の形式を :term[閉形式]{en="closed form" yomi="へいけいしき"}、$\alpha=\mathrm d\beta$ と表せる形式を :term[完全形式]{en="exact form" yomi="かんぜんけいしき"} という。完全なら閉である。一方、穴の周りの角の 1 形式は局所では角の微分だから閉だが、大域では線積分が $2\pi$ になるため完全ではない。この差を [de Rham コホモロジー](/math/manifolds/de-rham/) で測る。

### ベクトル場への評価で外微分を書く

:::proposition[外微分の座標を使わない評価式]{#prop-mfd-exterior-bracket-formula}
$k$ 形式 $\omega$ とベクトル場 $X_0,\ldots,X_k$ に対し、
$$
\begin{aligned}
\mathrm d\omega(X_0,\ldots,X_k)
={}&\sum_{i=0}^k(-1)^i X_i\bigl(\omega(X_0,\ldots,\widehat X_i,\ldots,X_k)\bigr)\\
&+\sum_{0\le i<j\le k}(-1)^{i+j}
\omega([X_i,X_j],X_0,\ldots,\widehat X_i,\ldots,\widehat X_j,\ldots,X_k)
\end{aligned}
\label{eq-mfd-exterior-bracket-formula}
$$
である。第 2 の和では、括弧を最初のスロットに入れ、残りの場は元の順を保つ。
:::

:::proof
右辺が各スロットで $C^\infty(M)$ 線形であることを確認する。$X_a$ を $fX_a$ に取り替えると、$f$ を外に出した項のほかに、他の場が $f$ を微分する項が現れる。$i\ne a$ の第 1 の和の項では、その余分な項の係数は $(-1)^iX_i(f)$ で、$X_a$ は $X_i$ を省いた残りの順の位置に入る。

$a<i$ なら第 2 の和の $(a,i)$ の項に $[fX_a,X_i]=f[X_a,X_i]-X_i(f)X_a$ が現れる。先頭にある $X_a$ を第 1 の和における位置まで $a$ 回動かすと符号は $(-1)^a$ である。余分な項全体の係数は $-(-1)^{a+i+a}X_i(f)=-(-1)^iX_i(f)$ となり相殺する。$i<a$ なら $(i,a)$ の項に $[X_i,fX_a]=f[X_i,X_a]+X_i(f)X_a$ が現れ、先頭から $a-1$ 回動かすので $(-1)^{i+a+a-1}X_i(f)=-(-1)^iX_i(f)$ となる。従ってこの場合も相殺する。これで右辺は各スロットについてテンソルとして点での値だけから決まる。

従って局所の座標ベクトル場を入力として比較すればよい。座標場の括弧は全部零で、第 2 の和は消える。第 1 の和は :eqref[eq-mfd-exterior-local] を入力に適用して、係数を微分する $\mathrm dx^j$ の入る位置ごとに展開した式である。例えば $\mathrm dx^j$ が $i$ 番目の入力を取る項では、最初のスロットからそこまで $i$ 回移動する符号 $(-1)^i$ が付く。よって両辺は座標入力で等しく、任意の入力でも等しい。
:::

:::example[座標でない枠での外微分]
$x>0$ の平面で $X=\partial_x$、$Y=x\partial_y$ とする。双対枠のうち $Y$ に $1$ を返す形式は $\beta=x^{-1}\mathrm dy$ である。評価 $\beta(X)=0,\beta(Y)=1$ は定数なので、式 :eqref[eq-mfd-exterior-bracket-formula] で $\mathrm d\beta(X,Y)=-\beta([X,Y])$ である。$[X,Y]=\partial_y$ なので値は $-1/x$ となる。座標で直接 $\mathrm d\beta=-x^{-2}\mathrm dx\wedge\mathrm dy$ を計算すると、$(X,Y)$ に対して確かに $-x^{-2}\cdot x=-1/x$ を返す。係数が枠で定数でも、その枠が座標から来なければ外微分は零とは限らない。
:::

### 座標の外積と枠の外積を微分する

任意の滑らかな余接枠 $\theta^1,\ldots,\theta^n$ で、$\omega=\sum_I a_I\theta^{i_1}\wedge\cdots\wedge\theta^{i_k}$ と書ける。座標の外積と違い、一般に $\mathrm d\theta^i$ は零ではない。このため積公式では二種類の項が現れる。第一は $\mathrm da_I\wedge\theta^I$ であり、第二は $a_I$ を固定して各 $\theta^{i_\nu}$ に $\mathrm d$ を作用させ、前の $\nu-1$ 個を通る符号 $(-1)^{\nu-1}$ を掛けた項である。$\theta^{i_\nu}$ が座標関数の微分なら第二の項がすべて消える。一般の枠で係数の微分だけを計算すると、空間で枠自体が変わる分を失う。

例えば極座標の平面 $r>0$ で $\theta^1=\mathrm dr$、$\theta^2=r\,\mathrm d\phi$ と取る。この二つは平面の標準計量での正規直交余接枠だが、$\mathrm d\theta^2=\mathrm dr\wedge\mathrm d\phi=r^{-1}\theta^1\wedge\theta^2$ である。対応する接枠は $E_1=\partial_r$、$E_2=r^{-1}\partial_\phi$ で、その括弧は $[E_1,E_2]=-r^{-1}E_2$ となる。従って :eqref[eq-mfd-exterior-bracket-formula] の括弧項も $\mathrm d\theta^2(E_1,E_2)=r^{-1}$ を返す。例えば $r=2$ でこの値は $1/2$ である。形式の係数をこの余接枠で $1$ と書いていても、外微分が正の値になるのは枠が座標由来でないためである。

一般に接枠 $E_i$ とその双対 $\theta^a$ を取り、$[E_i,E_j]=\sum_a c^a_{ij}E_a$ と書く。$\theta^a(E_j)$ が定数なので
$$
\mathrm d\theta^a=-\sum_{i<j}c^a_{ij}\,\theta^i\wedge\theta^j
$$
である。$E_i,E_j$ に評価すれば左辺は $-\theta^a([E_i,E_j])$、右辺は $-c^a_{ij}$ になる。この公式で枠の括弧の係数を一形式の外微分へ移せる。例えば $\R^3$ で $\theta^1=\mathrm dx,\theta^2=\mathrm dy,\theta^3=\mathrm dz-x\mathrm dy$ とすると、$\mathrm d\theta^3=-\mathrm dx\wedge\mathrm dy$ である。双対の $E_1=\partial_x$、$E_2=\partial_y+x\partial_z$、$E_3=\partial_z$ では $[E_1,E_2]=E_3$ であり、係数 $c^3_{12}=1$ の負が外微分に出ている。

### 外微分の値が一点の近くで決まること

形式の $\mathrm d$ は、その点での係数の値だけには依存しない。例えば $\alpha=x\,\mathrm dy$ と $\beta=0$ は原点で同じ零の値を持つが、$\mathrm d\alpha$ は原点で $\mathrm dx\wedge\mathrm dy$、$\mathrm d\beta=0$ である。一方、値と全ての座標偏微分が一点で同じなら、式 :eqref[eq-mfd-exterior-local] で外微分の値も同じである。座標変更で途中の偏導関数には複数の項が現れるが、外微分に残るのは交代化された一次の情報だけである。

例えば二形式 $\omega=a\,\mathrm dy\wedge\mathrm dz+b\,\mathrm dz\wedge\mathrm dx+c\,\mathrm dx\wedge\mathrm dy$ では、$a,b,c$ の九種類の一階部分微分があっても、外微分の係数は $\partial_x a+\partial_y b+\partial_z c$ だけである。形式に現れる座標を再び外積すると同じ微分が重複するため、例えば $a$ の $y,z$ による偏微分は消える。$a=xy,b=yz,c=zx$ なら外微分の係数は $x+y+z$ となる。座標を固定して計算する際には、各基底外積に含まれない添字から先に調べると必要な項を特定できる。

### 閉性を座標で調べる

$\mathrm d\omega=0$ なら $\omega$ を :term[閉形式]{en="closed form" yomi="へいけいしき"} といい、$\omega=\mathrm d\eta$ と書けるなら :term[完全形式]{en="exact form" yomi="かんぜんけいしき"} という。$\mathrm d^2=0$ から後者は前者である。この逆がどこまで成り立つかは [de Rham コホモロジー](/math/manifolds/de-rham/) で調べる。まず座標では、外微分のすべての昇順係数が零であるかを調べれば閉性が分かる。

:::example[各係数が変化しても二形式が閉じる]
$\R^3$ で $\omega=z\,\mathrm dx\wedge\mathrm dy-x\,\mathrm dy\wedge\mathrm dz$ とする。外微分の係数は $\partial_z z-\partial_x x=1-1=0$ だから閉形式である。この場合 $\eta=xz\,\mathrm dy$ の微分は $(z\,\mathrm dx+x\,\mathrm dz)\wedge\mathrm dy=\omega$ となり、完全でもある。各係数を別々に定数かどうか確認するのでなく、同じ外積に寄与する導関数を符号付きで足す必要がある。
:::

:::example[座標微分でない一形式が閉じない理由]
$(x,y,z)$ 上の $\alpha=\mathrm dz-x\,\mathrm dy$ について $\mathrm d\alpha=-\mathrm dx\wedge\mathrm dy$ である。従ってどの開近傍でも $\alpha=\mathrm df$ と書けない。仮にそのような $f$ があれば $f_x=0$、$f_y=-x$ となるので、$f_{xy}$ と $f_{yx}$ の一方が零、他方が $-1$ となる。外微分はこの矛盾を一つの二形式に記録している。
:::

## 引き戻しと外微分の自然性

:::proposition[引き戻しと外積・外微分]{#prop-mfd-exterior-naturality}
滑らかな $F:M\to N$ について、共変テンソルの引き戻しは交代性を保ち、$F^*(\alpha\wedge\beta)=F^*\alpha\wedge F^*\beta$ および $F^*\mathrm d\alpha=\mathrm d(F^*\alpha)$ が成り立つ。
:::

:::proof
各スロットに $\mathrm dF$ を適用する操作はスロットの入れ替えと交換する。従って交代性を保ち、外積の置換和の各項も引き戻した各因子の積に一致する。外微分については $\alpha=a\,\mathrm dy^{i_1}\wedge\cdots\wedge\mathrm dy^{i_k}$ の一項を扱えばよい。$F^*\mathrm dy^i=\mathrm d(y^i\circ F)$ であり、これを外微分すると $\mathrm d^2=0$ である。従って $\mathrm d(F^*\alpha)$ では係数の $\mathrm d(a\circ F)=F^*\mathrm da$ だけが残り、これは $F^*\mathrm d\alpha$ と等しい。局所の項の和で一般の場合が従う。
:::

:::example[面積形式の引き戻し]
$F:(r,\theta)\mapsto(r\cos\theta,r\sin\theta)$ では $F^*(\mathrm dx\wedge\mathrm dy)=r\,\mathrm dr\wedge\mathrm d\theta$ である。$F^*(x\,\mathrm dy)$ を先に計算して外微分しても同じである。実際 $F^*(x\,\mathrm dy)=r\cos\theta\sin\theta\,\mathrm dr+r^2\cos^2\theta\,\mathrm d\theta$ なので二形式の係数は $2r\cos^2\theta-r(\cos^2\theta-\sin^2\theta)=r$ である。
:::

### 写像の階数と引き戻しの消失

$F:M\to N$ の微分が点 $p$ で階数 $r$ なら、$k>r$ に対して任意の $k$ 形式の引き戻しは $p$ で零になる。実際、$\mathrm dF_p v_1,\ldots,\mathrm dF_p v_k$ は $r$ 次元の部分空間に入って一次従属だからである。逆に $k\le r$ なら $k$ 個の独立な像ベクトルを取り、値域の双対余接ベクトルの外積で値 $1$ を返す形式を作れる。必要なら局所座標でその点の形式を定数係数にして延長し、隆起関数で大域にもできる。従って全ての $k$ 形式の引き戻しがその点で零になるかどうかは、階数の情報を与える。

例えば $F(x,y,z)=(x,x^2,y)$ の階数は $2$ であり、値域の標準の三形式は引き戻すと零である。実際、最初の二成分の微分は $\mathrm dx$ と $2x\,\mathrm dx$ なので外積が零になる。一方、第 1、第 3 成分の微分の外積は $\mathrm dx\wedge\mathrm dy$ が残る。こうすると階数のどの部分が必要な形式を引き戻すのに十分か分かる。

:::example[二つの順番で同じ二形式を得る]
$F(x,y)=(u,v)=(x^2+y,xy)$、$\alpha=u\,\mathrm dv$ とする。$\mathrm d\alpha=\mathrm du\wedge\mathrm dv$ なので、引き戻しは $(2x\,\mathrm dx+\mathrm dy)\wedge(y\,\mathrm dx+x\,\mathrm dy)=(2x^2-y)\,\mathrm dx\wedge\mathrm dy$ である。他方、$F^*\alpha=y(x^2+y)\,\mathrm dx+x(x^2+y)\,\mathrm dy$ の外微分では $3x^2+y-(x^2+2y)=2x^2-y$ が残る。例えば $(1,2)$ で係数は零であり、そこで $F$ の微分の階数が落ちることと一致する。
:::

:::example[引き戻して閉性が生じる]
$\alpha=u\,\mathrm dv$ は平面で閉形式でない。一方、直線への写像 $i(t)=(t,2t)$ では $i^*\alpha=2t\,\mathrm dt=\mathrm d(t^2)$ である。実際 $i^*(\mathrm du\wedge\mathrm dv)=\mathrm dt\wedge2\,\mathrm dt=0$ である。自然性が示すのは、引き戻しと微分を交換できることである。写像の微分の像が小さいと、引き戻しが零でも元の二形式は非零でありうる。
:::

## 内部積と微分形式の Lie 微分

:::definition[内部積]{#def-mfd-interior-product}
ベクトル場 $X$ と $k$ 形式 $\omega$ に対し、$k>0$ なら
$$
(\iota_X\omega)(v_1,\ldots,v_{k-1})=\omega(X,v_1,\ldots,v_{k-1})
$$
と定め、これを :term[内部積]{en="interior product" yomi="ないぶせき"} という。関数の内部積は零とする。
:::

座標の外積基底に適用すると
$$
\iota_X(\mathrm dx^{i_1}\wedge\cdots\wedge\mathrm dx^{i_k})
=\sum_{a=1}^k(-1)^{a-1} X^{i_a}\,\mathrm dx^{i_1}\wedge\cdots\wedge
\widehat{\mathrm dx^{i_a}}\wedge\cdots\wedge\mathrm dx^{i_k}
\label{eq-mfd-interior-coordinates}
$$
である。帽子はその因子を省くことを表す。この式は行列式の最初の列を $X$ として余因子展開した結果である。また外積のシャッフル表示を $X$ が第 1 因子に入る項と第 2 因子に入る項に分けると
$$
\iota_X(\alpha\wedge\beta)=\iota_X\alpha\wedge\beta+(-1)^k\alpha\wedge\iota_X\beta
$$
となる。第 2 因子に入るには $X$ を最初の $k$ スロットの後ろに動かすため、符号 $(-1)^k$ が付く。

:::definition[形式の Lie 微分]{#def-mfd-lie-derivative-form}
$X$ の流れ $\Phi_t$ を用い、形式の :term[Lie 微分]{en="Lie derivative of a form" yomi="りーびぶん"} を
$$
\mathcal L_X\omega=\left.\frac{\mathrm d}{\mathrm dt}\right|_{t=0}\Phi_t^*\omega
$$
と定める。
:::

各点の近くで共通の小さい時間について流れが滑らかなので、微分は滑らかな係数を与える。引き戻しが外積を保つため、$\mathcal L_X$ は積の微分法則を満たす次数 $0$ の演算になる。また :ref[prop-mfd-exterior-naturality] と滑らかな偏微分の交換により $\mathcal L_X\mathrm d=\mathrm d\mathcal L_X$ である。例えば $X=x\partial_x$、$\omega=\mathrm dx$ の場合、$\Phi_t(x)=e^t x$ なので $\Phi_t^*\omega=e^t\mathrm dx$、従って $\mathcal L_X\omega=\mathrm dx$ である。

### 流れと共変テンソルの評価を微分する

形式だけでなく共変テンソルを流れで引き戻して微分できる。共変 $k$ テンソル $T$ とベクトル場 $Y_1,\ldots,Y_k$ に対し
$$
(\mathcal L_X T)(Y_1,\ldots,Y_k)
=X(T(Y_1,\ldots,Y_k))-\sum_{i=1}^k T(Y_1,\ldots,[X,Y_i],\ldots,Y_k)
\label{eq-mfd-lie-derivative-evaluation}
$$
となる。これを導くため、評価値の関数を $\Phi_t$ で引き戻すと、引き戻したテンソルに引き戻した場 $\Phi_t^*Y_i$ を入れたものになることを使う。実際、場を $\Phi_{-t}$ の微分で戻した後、テンソルの引き戻しでは $\Phi_t$ の微分を適用するため、二つの微分が逆同士で消える。時刻零で微分し、右辺でテンソル自体を微分する項と各スロットの場を微分する項に分ける。場の引き戻しの微分は :ref[thm-mfd-commuting-flows] により $[X,Y_i]$ である。これらを反対側に移すと式 :eqref[eq-mfd-lie-derivative-evaluation] を得る。

例えば $X=x\partial_x$、$Y=\partial_x$、$T=\mathrm dx$ とすると $T(Y)=1$ は定数であり、$[X,Y]=-Y$ である。従って式の右辺は $-T(-Y)=1$ となり、$\mathcal L_X\mathrm dx=\mathrm dx$ の値と一致する。$\mathrm dx$ の座標係数は定数だが、入力の場を引き戻したときに変わるため、Lie 微分は零にならない。

### 二つの内部積を続ける

$\iota_X\iota_Y\omega+\iota_Y\iota_X\omega=0$ である。左辺を残りの場に適用すると、$Y,X$ の順で入れた項と $X,Y$ の順で入れた項の和になり、反対称性で消える。$X=Y$ とすれば $\iota_X^2=0$ が従う。例えば $\omega=\mathrm dx\wedge\mathrm dy\wedge\mathrm dz$ で $X=\partial_x,Y=\partial_y$ に対し、$\iota_X\iota_Y\omega=-\mathrm dz$、$\iota_Y\iota_X\omega=\mathrm dz$ である。ここでは記号の右の内部積から適用している。三形式に内部積を二度適用して得る一形式の符号は、最初の二入力をどちらの順に入れるかに依存する。

### 内部積と関係するベクトル場

$F:M\to N$ で場 $X,Y$ が関係するとき、$F^*(\iota_Y\omega)=\iota_X(F^*\omega)$ である。点 $p$ で左辺を評価すると、最初の入力は $Y_{F(p)}$ で、残りは $F$ の微分で移した入力である。右辺では最初の入力を $X_p$ としてから同じ微分を適用する。関係性の式 $\mathrm dF_pX_p=Y_{F(p)}$ により一致する。例えば射影 $F(x,y)=x$ と $X=x\partial_x+y\partial_y$、$Y=x\partial_x$ に 1 形式 $\mathrm dx$ を使うと、どちらの順序でも関数 $x$ を得る。この交換には $F$ が微分同相である必要はない。

流れでも関係性が保たれるので、形式の Lie 微分もこの場合に $F$ の引き戻しと交換する。実際、局所流の合成の式 $F\circ\Phi_t=\Psi_t\circ F$ は $X,Y$ が関係するから成り立つ。両辺で $\omega$ を引き戻して $t=0$ で微分すれば
$$
\mathcal L_X(F^*\omega)=F^*(\mathcal L_Y\omega)
$$
である。あるいは Cartan の公式を両側に使い、内部積と外微分がそれぞれ交換することからも得られる。例えば $F(x,y)=(x,y,x^2)$ の像上に接する場 $Y=\partial_u+2u\partial_w$ と $X=\partial_x$ は関係する。値域の形式 $\mathrm dw$ で $\mathcal L_Y\mathrm dw=\mathrm d(2u)=2\mathrm du$ であり、引き戻して $2\mathrm dx$ を得る。逆に $F^*\mathrm dw=2x\mathrm dx$ を $X$ で Lie 微分しても $2\mathrm dx$ になる。

### 内部積の核と流れによる変化

:::example[三形式から面の向きと倍率を得る]
$\mu=\mathrm dx\wedge\mathrm dy\wedge\mathrm dz$、$X=2\partial_x-\partial_z$ とする。$\iota_X\mu=2\,\mathrm dy\wedge\mathrm dz-\mathrm dx\wedge\mathrm dy$ である。$Y=\partial_y$ をさらに挿すと $\iota_Y\iota_X\mu=2\,\mathrm dz+\mathrm dx$ となる。これを $X$ に適用すると $-2+2=0$ になる。これは同じ $X$ を二度入れれば交代性で零になるという性質の成分表示である。
:::

:::example[係数と基底の両方が流れで変わる]
$X=x\partial_x+2y\partial_y$ の流れは $\Phi_t(x,y)=(e^t x,e^{2t}y)$ である。$\alpha=x\,\mathrm dy-y\,\mathrm dx$ に対し $\Phi_t^*\alpha=e^{3t}\alpha$ なので $\mathcal L_X\alpha=3\alpha$ である。例えば $x\,\mathrm dy$ では係数が $e^t$ 倍になり、基底 $\mathrm dy$ は $e^{2t}$ 倍になる。係数だけを微分すると $x\,\mathrm dy$ しか得られず、実際の $3x\,\mathrm dy$ と一致しない。
:::

## Cartan の公式

:::theorem[Cartan の公式]{#thm-mfd-cartan-formula}
任意のベクトル場 $X$ と形式 $\omega$ に対し
$$
\mathcal L_X\omega=\mathrm d(\iota_X\omega)+\iota_X(\mathrm d\omega)
\label{eq-mfd-cartan}
$$
が成り立つ。
:::

:::proof
右辺を $D=\mathrm d\iota_X+\iota_X\mathrm d$ と置く。二つの積公式を $\alpha\wedge\beta$ に適用すると、$\iota_X\alpha\wedge\mathrm d\beta$ の係数は $(-1)^{k-1}+(-1)^k=0$、$\mathrm d\alpha\wedge\iota_X\beta$ の係数は $(-1)^k+(-1)^{k+1}=0$ になる。残る項は $D\alpha\wedge\beta+\alpha\wedge D\beta$ なので、$D$ も次数 $0$ の導分である。関数 $f$ について $Df=\iota_X\mathrm df=Xf=\mathcal L_X f$ である。また $\mathrm d^2=0$ より $D\mathrm d=\mathrm d\iota_X\mathrm d=\mathrm dD$ である。従って $\mathrm dx^i$ に対しても $D$ と $\mathcal L_X$ はともに $\mathrm d(Xx^i)$ を返す。各座標で関数と $\mathrm dx^i$ は全形式を外積と和によって生成するので、二つの演算が一致する。
:::

:::example[面積の変化を求める]
$X=x\partial_x+2y\partial_y$、$\omega=\mathrm dx\wedge\mathrm dy$ とする。内部積は $x\,\mathrm dy-2y\,\mathrm dx$ であり、その外微分は $3\omega$ である。$\mathrm d\omega=0$ なので $\mathcal L_X\omega=3\omega$ となる。流れで面積が $e^t e^{2t}=e^{3t}$ 倍になることから、時刻 $0$ の変化率は同じ $3$ になる。
:::

一般に $\frac{\mathrm d}{\mathrm dt}\Phi_t^*\omega=\Phi_t^*\mathcal L_X\omega$ である。実際、$\Phi_{t+h}^*=\Phi_t^*\Phi_h^*$ の $h=0$ での微分で得られる。従って流れが形式を保つことと Lie 微分が零であることは同値である。閉形式では Cartan の公式の右辺は $\mathrm d(\iota_X\omega)$ だけになる。このため、閉形式の保存を内部積が閉かどうかで判定できる。

### Lie 微分・外微分・内部積の交換関係

:::proposition[Cartan の交換関係]{#prop-mfd-cartan-commutators}
形式に対して
$$
\begin{aligned}
\mathcal L_X\mathrm d&=\mathrm d\mathcal L_X,\\
\mathcal L_X\iota_Y-\iota_Y\mathcal L_X&=\iota_{[X,Y]},\\
\mathcal L_X\mathcal L_Y-\mathcal L_Y\mathcal L_X&=\mathcal L_{[X,Y]}
\end{aligned}
$$
が成り立つ。
:::

:::proof
最初は外微分の自然性を流れで微分するか、Cartan の公式と $\mathrm d^2=0$ を使えばよい。第 2 式について、引き戻しは各入力に $\mathrm d\Phi_t$ を適用し、$\Phi_t^*Y$ の成分に $\mathrm d\Phi_t$ を適用すれば元の点の $Y$ が戻る。従って $\Phi_t^*(\iota_Y\omega)=\iota_{\Phi_t^*Y}\Phi_t^*\omega$ である。時刻零で微分し、$\frac{\mathrm d}{\mathrm dt}|_0\Phi_t^*Y=[X,Y]$ を使えば式が得られる。第 3 式では $\mathcal L_Y=\mathrm d\iota_Y+\iota_Y\mathrm d$ を入れる。$\mathcal L_X$ を $\mathrm d$ を通して移すと、差は
$$
\mathrm d(\mathcal L_X\iota_Y-\iota_Y\mathcal L_X)
+(\mathcal L_X\iota_Y-\iota_Y\mathcal L_X)\mathrm d
=\mathrm d\iota_{[X,Y]}+\iota_{[X,Y]}\mathrm d
$$
になる。右辺は Cartan の公式で $\mathcal L_{[X,Y]}$ である。
:::

例えば $X=\partial_x,Y=x\partial_y$ と $\omega=\mathrm dy$ について、$\iota_Y\omega=x$、$\mathcal L_X x=1$、$\mathcal L_X\omega=0$ である。従って第 2 式の左辺は $1$ となる。右辺も $[X,Y]=\partial_y$ を入れれば $1$ である。一方 $x$ を掛けた $X$ による Lie 微分は、一般には $x\mathcal L_X$ とは一致しない。Cartan の公式を使うと任意の関数 $f$ と形式に対し
$$
\mathcal L_{fX}\omega=f\mathcal L_X\omega+\mathrm df\wedge\iota_X\omega
$$
である。これは最初の外微分が $f$ にも作用するからであり、$\omega$ の閉性を仮定する必要はない。

### 外微分と内部積で流れを計算する

:::example[同じ例を Cartan の公式で検算する]
$\alpha=x\,\mathrm dy-y\,\mathrm dx$ と $X=x\partial_x+2y\partial_y$ では $\iota_X\alpha=xy$、$\mathrm d\alpha=2\,\mathrm dx\wedge\mathrm dy$ である。従って
$$
\mathrm d\iota_X\alpha+\iota_X\mathrm d\alpha
=(y\,\mathrm dx+x\,\mathrm dy)+(2x\,\mathrm dy-4y\,\mathrm dx)
=3\alpha.
$$
一つの項だけでは正しい倍率が出ない。$\mathrm d\iota_X\alpha$ は評価値の変化を、$\iota_X\mathrm d\alpha$ は二形式から戻した成分を与え、和が流れの微分になる。
:::

:::example[閉形式の流れによる変化は完全になる]
$\omega=\mathrm dx\wedge\mathrm dy$、$X=x\partial_x$ の場合は $\mathrm d\omega=0$ で、$\iota_X\omega=x\,\mathrm dy$ である。Cartan の公式より $\mathcal L_X\omega=\mathrm d(x\,\mathrm dy)=\omega$ となる。実際、流れは $x$ だけを $e^t$ 倍するので、$\Phi_t^*\omega-\omega=(e^t-1)\omega=\mathrm d((e^t-1)x\,\mathrm dy)$ である。一般の閉形式でも $\partial_t\Phi_t^*\omega=\mathrm d(\Phi_t^*\iota_X\omega)$ なので、同じ開集合上で流れが時刻 $0$ から $t$ まで定義されれば、積分して差の原始形式を得られる。外微分と時刻積分を交換できるのは、有限区間上で係数とその偏導関数が滑らかだからである。
:::

## 微分形式と ℝ³ のベクトル解析

標準の座標 $(x,y,z)$ と体積形式 $\mu=\mathrm dx\wedge\mathrm dy\wedge\mathrm dz$ を固定する。$A=(P,Q,R)$ に対し
$$
A^\flat=P\,\mathrm dx+Q\,\mathrm dy+R\,\mathrm dz,\qquad
\iota_A\mu=P\,\mathrm dy\wedge\mathrm dz+Q\,\mathrm dz\wedge\mathrm dx+R\,\mathrm dx\wedge\mathrm dy
$$
と置く。最初の対応は Euclid 計量、二番目は体積形式を使う。それぞれ標準座標で各成分を独立に復元できるため、一形式、二形式との線形同型である。

:::proposition[勾配・回転・発散と外微分]{#prop-mfd-vector-calculus-forms}
Euclid 空間 $\R^3$ で
$$
\mathrm df=(\grad f)^\flat,\qquad
\mathrm d(A^\flat)=\iota_{\rot A}\mu,\qquad
\mathrm d(\iota_A\mu)=(\divg A)\mu .
$$
:::

:::proof
最初の式は関数の微分の座標表示である。二番目では $\mathrm d(P\,\mathrm dx+Q\,\mathrm dy+R\,\mathrm dz)$ を式 :eqref[eq-mfd-exterior-local] で展開する。$\mathrm dy\wedge\mathrm dz$、$\mathrm dz\wedge\mathrm dx$、$\mathrm dx\wedge\mathrm dy$ の係数は順に $\partial_yR-\partial_zQ$、$\partial_zP-\partial_xR$、$\partial_xQ-\partial_yP$ となる。これは回転の成分である。最後の式では各項の係数の微分のうち、すでに微分が含まれる座標方向は外積で消える。そのため $\partial_xP+\partial_yQ+\partial_zR$ が残る。三つの座標順序はいずれも巡回順なので符号は正である。
:::

例えば $A=(-y,x,0)$ の回転は $(0,0,2)$ であり、$\mathrm dA^\flat=2\,\mathrm dx\wedge\mathrm dy$ である。$A=(x,2y,3z)$ では発散は $6$ であり、$\mathrm d(\iota_A\mu)=6\mu$ である。さらに $\mathrm d^2=0$ を関数から始めれば $\rot(\grad f)=0$、一形式から始めれば $\divg(\rot A)=0$ になる。次数を区別しておくと、これらが同じ恒等式の別の表示であることが分かる。

### 計量と向きを使う Hodge の星

Euclid 空間で一形式をベクトル場に対応させたときには内積を使い、二形式を対応させたときには体積形式を使った。次にこれらの操作を、全ての次数で形式から形式への写像としてまとめる。有限次元の正定値内積空間 $V$ が向き付けられているとする。正の正規直交基底の双対を $e^1,\ldots,e^n$ とし、$\Lambda^k V^*$ では昇順の外積 $e^I$ を正規直交基底とする。この内積はどの正規直交基底から作っても同じである。実際、単純な外積同士の内積は、各一形式同士の内積を成分とする $k$ 次行列の行列式である。これを確認するには、基底の外積同士で両方が同じ $I$ なら $1$、異なるなら $0$ になることを見ればよい。双線形性により任意の外積へ拡張される。

体積元 $\mu=e^1\wedge\cdots\wedge e^n$ を使って、:term[Hodge の星作用素]{en="Hodge star operator" yomi="ほっじのほしさようそ"} $\star:\Lambda^kV^*\to\Lambda^{n-k}V^*$ を $\alpha\wedge\star\beta=\inner{\alpha}{\beta}\mu$ で定義できる。存在は $\star e^I=\epsilon(I,I^c)e^{I^c}$ と置けば確認できる。ここで $\epsilon(I,I^c)$ は、並び $(I,I^c)$ を標準の並びに直す符号である。この値なら $\alpha=e^I$ の場合は外積が $\mu$、別の外積基底の場合は添字が重複して零である。各 $\alpha$ との外積の値は $\star\beta$ の全成分を決めるので一意でもある。そして $I$ と $I^c$ の順番を替えると符号が $(-1)^{k(n-k)}$ 倍になるため $\star\star\beta=(-1)^{k(n-k)}\beta$ となる。向きを逆に取れば $\mu$ が負になるため、$\star$ も全次数で負になる。

三次元の標準座標では $\star\mathrm dx=\mathrm dy\wedge\mathrm dz$、$\star\mathrm dy=\mathrm dz\wedge\mathrm dx$、$\star\mathrm dz=\mathrm dx\wedge\mathrm dy$ である。従って $\star(A^\flat)=\iota_A\mu$、$\star\mathrm d(A^\flat)=(\rot A)^\flat$、$\star\mathrm d(\star A^\flat)=\divg A$ と書ける。例えば $A=(0,0,x^2y)$ なら $A^\flat=x^2y\,\mathrm dz$ であり、
$$
\star\mathrm d(A^\flat)=x^2\,\mathrm dx-2xy\,\mathrm dy .
$$
その回転の成分は $(x^2,-2xy,0)$ で、発散は $2x-2x=0$ である。このように $\star$ は形式とベクトル場の対応を揃えるが、$\mathrm d$ とは異なり計量と向きの選択を使っている。

### 外積法則から得るベクトル恒等式

$\R^3$ でスカラー関数 $f$ と場 $A$ に対し
$$
\rot(fA)=\grad f\times A+f\rot A,\qquad
\divg(fA)=\grad f\cdot A+f\divg A
$$
となる。最初は $\mathrm d(fA^\flat)=\mathrm df\wedge A^\flat+f\,\mathrm dA^\flat$ に $\star$ を適用する。$\star(\alpha\wedge\beta)$ の成分は $\alpha,\beta$ の成分の二次小行列式であり、これは $\alpha,\beta$ に対応するベクトルの外積の成分である。二番目は $\mathrm d(f\iota_A\mu)=\mathrm df\wedge\iota_A\mu+f\,\mathrm d(\iota_A\mu)$ であり、最初の項は $(\mathrm df(A))\mu$ になる。基底で確かめるには、$\mathrm df$ の係数と $\iota_A\mu$ の係数のうち同じ添字の積だけが重複を避け、巡回順の正符号で残ることを使う。

十分滑らかな場 $A,B$ について $\iota_{A\times B}\mu=A^\flat\wedge B^\flat$ であることも、同じ二次小行列式から従う。その外微分は $\mathrm dA^\flat\wedge B^\flat-A^\flat\wedge\mathrm dB^\flat$ であるから
$$
\divg(A\times B)=B\cdot\rot A-A\cdot\rot B
$$
となる。例えば $A=(0,x,0)$、$B=(0,0,y)$ なら $A\times B=(xy,0,0)$ で発散は $y$ である。$\rot A=(0,0,1)$、$\rot B=(1,0,0)$ だから、右辺も $y-0=y$ になる。積公式から出る負号は、$A^\flat$ が奇数次数であることに対応する。

### 一つの計量から発散と Laplace 作用素を作る

標準の三次元では $\star\mathrm df$ は $\iota_{\grad f}\mu$ に等しい。その外微分に再び星を適用すると
$$
\star\mathrm d\star\mathrm df=\divg(\grad f)=\partial_x^2 f+\partial_y^2 f+\partial_z^2f .
$$
この符号では通常の Laplace 作用素である。例えば $f=x^2+2y^2-3z^2$ なら、勾配 $(2x,4y,-6z)$ の発散は $2+4-6=0$ である。勾配が零ではないのにこの和が消えることから、$f$ が調和であることと $f$ が定数であることを区別できる。ここでの計算は Euclid 計量と標準の向きを用いている。別の計量では余接ノルムや体積形式が変わるので、Laplace 作用素の係数も変わる。外微分そのものの局所式は変わらないため、どの段階で計量が入るかを区別する。

例えば平面の極座標で $g=\mathrm dr^2+r^2\mathrm d\theta^2$ とすると、正規直交余接枠は $(\mathrm dr,r\mathrm d\theta)$ である。この向きで $\star\mathrm dr=r\mathrm d\theta$、$\star\mathrm d\theta=-\mathrm dr/r$ となる。従って $f=f(r,\theta)$ について
$$
\mathrm d(\star\mathrm df)=
\left(\partial_r(r\partial_r f)+\frac1r\partial_\theta^2 f\right)
\,\mathrm dr\wedge\mathrm d\theta
$$
である。体積形式が $r\mathrm dr\wedge\mathrm d\theta$ なので、最高次形式を $\star$ で関数に直せば
$$
\star\mathrm d\star\mathrm df=\partial_r^2 f+\frac1r\partial_r f+\frac1{r^2}\partial_\theta^2 f.
$$
たとえば $f=r^2\cos(2\theta)$ に代入すると $2\cos(2\theta)+2\cos(2\theta)-4\cos(2\theta)=0$ である。直交座標の $x^2-y^2$ と同じ関数であり、二階微分の和も $2-2=0$ である。この計算に出る $1/r$ の項は原点で座標が使えないことを反映しており、関数そのものがそこで特異になることを意味しない。

### 三次元の曲がった座標で回転を計算する

:::example[円柱座標で高さとともに変わる回転]
$r>0$ の円柱座標で $A=z\partial_\theta$ を取る。Euclid 計量は $\mathrm dr^2+r^2\mathrm d\theta^2+\mathrm dz^2$ なので $A^\flat=zr^2\,\mathrm d\theta$ である。その外微分は
$$
\mathrm dA^\flat=2zr\,\mathrm dr\wedge\mathrm d\theta+r^2\,\mathrm dz\wedge\mathrm d\theta .
$$
正向きの正規直交余接枠 $(\mathrm dr,r\,\mathrm d\theta,\mathrm dz)$ で $\star$ を取ると、$2z\,\mathrm dz-r\,\mathrm dr$ を得る。従って回転は $-r\partial_r+2z\partial_z$ である。直交座標では $A=(-zy,zx,0)$ で、その回転は $(-x,-y,2z)$ になる。二つの結果が一致するのは $r\partial_r=x\partial_x+y\partial_y$ だからである。
:::

:::example[円柱座標で発散の余分な項が生じる理由]
同じ座標で $B=r^2\partial_r$ とする。$\mu=r\,\mathrm dr\wedge\mathrm d\theta\wedge\mathrm dz$ なので $\iota_B\mu=r^3\,\mathrm d\theta\wedge\mathrm dz$、$\mathrm d\iota_B\mu=3r^2\,\mathrm dr\wedge\mathrm d\theta\wedge\mathrm dz$ である。従って $\divg B=3r$ となる。$\partial_r(r^2)=2r$ との差 $r$ は、体積の係数 $r$ が半径とともに増す分である。例えば $r=2$ では発散は $6$ であり、座標成分だけの微分 $4$ ではない。
:::

## 演習問題

::::exercise{level=basic}
$(2\mathrm dx+\mathrm dy)\wedge(\mathrm dx+3\mathrm dy)$ を求めよ。
:::solution
同じ因子の外積は零である。残った交差項は $6\mathrm dx\wedge\mathrm dy+\mathrm dy\wedge\mathrm dx=(6-1)\mathrm dx\wedge\mathrm dy$ である。
:::
::::

::::exercise{level=basic}
$\alpha=x^2\,\mathrm dy+yz\,\mathrm dz$ の外微分を求めよ。
:::solution
$\mathrm d(x^2)\wedge\mathrm dy=2x\,\mathrm dx\wedge\mathrm dy$ であり、$\mathrm d(yz)\wedge\mathrm dz=(z\,\mathrm dy+y\,\mathrm dz)\wedge\mathrm dz=z\,\mathrm dy\wedge\mathrm dz$ である。従ってこれらの和になる。
:::
::::

::::exercise{level=basic}
$F(u,v)=(u+v,u-v)$ による $\mathrm dx\wedge\mathrm dy$ の引き戻しを求めよ。
:::solution
$(\mathrm du+\mathrm dv)\wedge(\mathrm du-\mathrm dv)=-\mathrm du\wedge\mathrm dv+\mathrm dv\wedge\mathrm du=-2\mathrm du\wedge\mathrm dv$ となる。係数はヤコビ行列 $\begin{pmatrix}1&1\\1&-1\end{pmatrix}$ の行列式と一致する。
:::
::::

::::exercise{level=basic}
$X=2\partial_x-\partial_z$ による $\mathrm dx\wedge\mathrm dy\wedge\mathrm dz$ の内部積を求めよ。
:::solution
式 :eqref[eq-mfd-interior-coordinates] より第 1 成分は $2\mathrm dy\wedge\mathrm dz$、第 2 成分は零、第 3 成分は $(-1)^2(-1)\mathrm dx\wedge\mathrm dy$ である。従って $2\mathrm dy\wedge\mathrm dz-\mathrm dx\wedge\mathrm dy$ になる。
:::
::::

::::exercise{level=advanced}
$X=x\partial_y$、$\alpha=y\,\mathrm dx$ の $\mathcal L_X\alpha$ を Cartan の公式と流れの両方で求めよ。
:::solution
$\iota_X\alpha=0$、$\mathrm d\alpha=\mathrm dy\wedge\mathrm dx$ なので内部積は $x\,\mathrm dx$ である。一方流れは $(x,y)\mapsto(x,y+tx)$ で、引き戻しは $(y+tx)\mathrm dx$ になる。その導関数は同じ $x\mathrm dx$ である。
:::
::::

::::exercise{level=advanced}
$A=(yz,zx,xy)$ の回転と発散を、外微分から求めよ。
:::solution
$A^\flat=\mathrm d(xyz)$ なので $\mathrm d(A^\flat)=0$、従って回転は零である。$\iota_A\mu=yz\,\mathrm dy\wedge\mathrm dz+zx\,\mathrm dz\wedge\mathrm dx+xy\,\mathrm dx\wedge\mathrm dy$ の各係数は、対応する微分にすでに現れる変数だけに依存する。従ってその外微分の各項は重複した外積で消え、発散も零になる。
:::
::::

::::exercise{level=advanced}
滑らかな関数 $f$、ベクトル場 $X$ と任意の形式 $\omega$ に対し、$\mathcal L_{fX}\omega=f\mathcal L_X\omega+\mathrm df\wedge\iota_X\omega$ を示せ。
:::solution
$\iota_{fX}\omega=f\iota_X\omega$ である。従って Cartan の公式から $\mathcal L_{fX}\omega=\mathrm d(f\iota_X\omega)+f\iota_X\mathrm d\omega$ となる。積公式により最初の項は $\mathrm df\wedge\iota_X\omega+f\mathrm d(\iota_X\omega)$ である。$f$ を含む二項をまとめると $f\mathcal L_X\omega$ となる。
:::
::::

::::exercise{level=advanced}
$M$ の次元が $n$ のとき、$n$ 形式はすべて閉であることを示せ。さらに $n>0$ のとき、$\R^n$ 上の定数係数の $n$ 形式は完全であることを示せ。
:::solution
外微分の値は $\Omega^{n+1}(M)=0$ に属するので閉である。$\R^n$ 上の定数係数の形式は $c\,\mathrm dx^1\wedge\cdots\wedge\mathrm dx^n$ である。$\beta=cx^1\,\mathrm dx^2\wedge\cdots\wedge\mathrm dx^n$ と置けば式 :eqref[eq-mfd-exterior-local] で $\mathrm d\beta$ が元の形式になる。$n=1$ のとき $\beta=cx^1$ と読む。
:::
::::

::::exercise{level=advanced}
$F(u,v,w)=(u,v,uv+w)$ による $z\,\mathrm dx\wedge\mathrm dy+x\,\mathrm dy\wedge\mathrm dz$ の引き戻しを求め、その外微分を元の外微分の引き戻しと照合せよ。
:::solution
$F^*\mathrm dz=v\,\mathrm du+u\,\mathrm dv+\mathrm dw$ なので、元の二形式 $\omega$ は
$$
F^*\omega=(uv+w)\mathrm du\wedge\mathrm dv+
u\mathrm dv\wedge(v\,\mathrm du+u\,\mathrm dv+\mathrm dw)
=w\,\mathrm du\wedge\mathrm dv+u\,\mathrm dv\wedge\mathrm dw
$$
に写る。最初にあった $uv$ の項は $\mathrm dv\wedge\mathrm du=-\mathrm du\wedge\mathrm dv$ によって消える。この二つの項の外微分を取ると $\mathrm dw\wedge\mathrm du\wedge\mathrm dv+\mathrm du\wedge\mathrm dv\wedge\mathrm dw$ であり、前者は巡回順なので合計 $2\,\mathrm du\wedge\mathrm dv\wedge\mathrm dw$ である。

元の空間では $\mathrm d\omega=\mathrm dz\wedge\mathrm dx\wedge\mathrm dy+\mathrm dx\wedge\mathrm dy\wedge\mathrm dz=2\,\mathrm dx\wedge\mathrm dy\wedge\mathrm dz$ である。これを引き戻すと $2\,\mathrm du\wedge\mathrm dv\wedge(v\,\mathrm du+u\,\mathrm dv+\mathrm dw)$ となり、最初の二項は重複で消える。結果は同じで、:ref[prop-mfd-exterior-naturality] の具体例である。
:::
::::

::::exercise{level=advanced}
$\omega=x\,\mathrm dy\wedge\mathrm dz$ と $X=x\partial_x-y\partial_y$ に対し、流れと Cartan の公式の二通りで $\mathcal L_X\omega$ を計算せよ。
:::solution
流れは $\Phi_t(x,y,z)=(e^tx,e^{-t}y,z)$ である。引き戻すと係数 $x$ に $e^t$、$\mathrm dy$ に $e^{-t}$ が掛かり、$\mathrm dz$ は変わらない。従って $\Phi_t^*\omega=\omega$ であり、Lie 微分は零である。

もう一方の方法では $\iota_X\omega=-xy\,\mathrm dz$ で、その外微分は $-y\,\mathrm dx\wedge\mathrm dz-x\,\mathrm dy\wedge\mathrm dz$ である。一方、$\mathrm d\omega=\mathrm dx\wedge\mathrm dy\wedge\mathrm dz$ の内部積は $x\,\mathrm dy\wedge\mathrm dz+y\,\mathrm dx\wedge\mathrm dz$ である。Cartan の公式で両方を足すと全項が消え、同じ零を得る。$\omega$ 自体は閉ではないので、$\mathrm d(\iota_X\omega)$ の一項だけを計算すると零にならない。
:::
::::

::::exercise{level=advanced}
$X=\partial_x+y\partial_y$、$\alpha=e^{-x}y\,\mathrm dy$ とする。$\mathcal L_X\alpha$、$\mathrm d\mathcal L_X\alpha$、$\mathcal L_X\mathrm d\alpha$ を計算し、後ろ二つが一致することを流れで確認せよ。
:::solution
流れは $\Phi_t(x,y)=(x+t,e^ty)$ である。係数 $e^{-x}y$ は引き戻すと $e^{-(x+t)}e^ty=e^{-x}y$ で変わらないが、$\mathrm dy$ は $e^t\mathrm dy$ になる。従って $\Phi_t^*\alpha=e^t\alpha$ であり、$\mathcal L_X\alpha=\alpha$ である。$\mathrm d\alpha=-e^{-x}y\,\mathrm dx\wedge\mathrm dy$ なので、$\mathrm d\mathcal L_X\alpha$ もこの値である。また $\mathrm dx$ は引き戻しで変わらず、$\mathrm dy$ にだけ $e^t$ が掛かる。係数は先ほどと同じなので $\Phi_t^*\mathrm d\alpha=e^t\mathrm d\alpha$ となり、$\mathcal L_X\mathrm d\alpha=\mathrm d\alpha$ である。

この計算では係数が流れで保存されても一形式は保存されなかった。引き戻しは点での係数を読むことに加え、入力する速度を流れの微分で移す。今回は後者の $e^t$ が残った。外微分との交換は、係数だけでなく微分の因子も正しく引き戻した両方の式で確認できる。
:::
::::

::::exercise{level=advanced}
開集合 $x>0$ に $\alpha=(-y\mathrm dx+x\mathrm dy)/(x^2+y^2)$ を取る。$\mathrm d\alpha$ を直接計算し、その領域での原始関数を求めよ。また同じ式の形式を $\R^2\setminus\{0\}$ 全体へ広げた場合と比較せよ。
:::solution
$P=-y/(x^2+y^2)$、$Q=x/(x^2+y^2)$ と置く。$\partial_xQ=(y^2-x^2)/(x^2+y^2)^2$ であり、$\partial_yP=(y^2-x^2)/(x^2+y^2)^2$ でもある。従って $\mathrm d\alpha=(\partial_xQ-\partial_yP)\mathrm dx\wedge\mathrm dy=0$ となる。$x>0$ では $f=\arctan(y/x)$ が滑らかで
$$
\mathrm df=\frac{1}{1+(y/x)^2}\left(\frac{\mathrm dy}{x}-\frac{y\,\mathrm dx}{x^2}\right)
=\alpha
$$
である。これは角の一つの枝である。

平面の原点を除いた全体でも分母は非零なので外微分の直接計算は成立し、形式は閉である。しかし単位円への引き戻しは $\mathrm dt$ で積分は $2\pi$ になる。原始関数の微分なら閉曲線の積分は零であるから、この全域で原始関数はない。$\arctan(y/x)$ で表せない点があるだけでなく、別の原始関数を選んでも一周の差は除けない。
:::
::::

::::exercise{level=advanced}
$\R^4$ の $\beta=(\mathrm dx+\mathrm dy)\wedge\mathrm dz+\mathrm dy\wedge\mathrm dw$ に対し、$\beta^2$ を求めて階数を判定せよ。さらに $X=\partial_x+2\partial_y$ を挿入した $\iota_X\beta$ と、その外微分を求めよ。

:::solution
二つの二次因子は交換しても符号が変わらず、それぞれの平方は零である。そのため
$$
\beta^2=2(\mathrm dx+\mathrm dy)\wedge\mathrm dz\wedge\mathrm dy\wedge\mathrm dw
=-2\,\mathrm dx\wedge\mathrm dy\wedge\mathrm dz\wedge\mathrm dw.
$$
第 1 の表示から第 2 の表示にするとき、$\mathrm dy$ が重複する項は消え、$\mathrm dz$ と $\mathrm dy$ の交換に負号が付く。平方は非零なので、:ref[prop-mfd-alternating-normal-form] により階数は $4$ である。

$X$ は $\mathrm dz,\mathrm dw$ の両方を消し、$\mathrm dx+\mathrm dy$ と $\mathrm dy$ への評価はそれぞれ $3$ と $2$ である。そのため $\iota_X\beta=3\,\mathrm dz+2\,\mathrm dw$ となる。係数は定数なので $\mathrm d\iota_X\beta=0$ である。また $\beta$ も定数係数で閉じているから、Cartan の公式より $\mathcal L_X\beta=0$ である。流れが平行移動で全ての座標微分を保つことからも同じ結論が得られる。
:::
::::

## 確認テスト

::::quiz{answer="6" tol="1e-6"}
$\R^4$ 上の各点で 2 形式の独立な成分はいくつあるか。
:::explanation
昇順の二添字 $i<j$ に対する $\mathrm dx^i\wedge\mathrm dx^j$ が基底である。従って $\binom42=6$ 成分である。
:::
::::

::::quiz
外微分について正しいものを選べ。

- [x] $\mathrm d^2=0$ である。
- [x] 滑らかな写像による引き戻しと交換する。
- [ ] 二形式に適用すると常に零になる。

:::explanation
混合偏微分の対称性と外積の反対称性により二度の外微分は零になる。引き戻しは各 $\mathrm dy^i$ を $\mathrm d(y^i\circ F)$ に写すため、係数の微分と積公式を適用すると外微分と交換する。一方 $\R^3$ では $\mathrm d(x\,\mathrm dy\wedge\mathrm dz)=\mathrm dx\wedge\mathrm dy\wedge\mathrm dz$ なので二形式の外微分は零とは限らない。
:::
::::

::::quiz{answer="3" tol="1e-6"}
$X=x\partial_x+2y\partial_y$ に沿う $\mathcal L_X(\mathrm dx\wedge\mathrm dy)$ の係数を求めよ。
:::explanation
内部積を取ると $x\,\mathrm dy-2y\,\mathrm dx$ になり、その外微分は $\mathrm dx\wedge\mathrm dy-2\mathrm dy\wedge\mathrm dx=3\mathrm dx\wedge\mathrm dy$ である。元の形式は閉なのでこれは Lie 微分に一致する。
:::
::::

::::quiz
$X$ と閉形式 $\omega$ について $\mathcal L_X\omega=0$ と同値な条件を選べ。

- [ ] $\iota_X\omega=0$
- [x] $\mathrm d(\iota_X\omega)=0$
- [ ] $\omega$ が完全である。

:::explanation
Cartan の公式で $\iota_X(\mathrm d\omega)=0$ となるため、Lie 微分は $\mathrm d(\iota_X\omega)$ と等しい。例えば $\R^2$ で $X=\partial_x$、$\omega=\mathrm dx\wedge\mathrm dy$ なら内部積は非零な $\mathrm dy$ だが閉である。完全性は場の選択に依存しない一方、例えば同じ形式は $x\partial_x$ の流れでは保存されないため、完全であることだけでは十分でない。
:::
::::

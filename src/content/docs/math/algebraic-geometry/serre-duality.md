---
title: 'Serre 双対性'
description: '射影空間の単項式の双対性から、Cohen–Macaulay 射影スキームの Serre 双対性を導く。'
notation:
  - symbol: '\omega_X'
    meaning: '標準層または双対化層'
  - symbol: '\mathcal E xt^j_X(\mathcal F,\mathcal G)'
    meaning: '局所 Ext の層'
  - symbol: '\operatorname{Ext}^j_X(\mathcal F,\mathcal G)'
    meaning: '層の加群の圏の大域 Ext'
---

## この章の目標

- 射影空間の最上次コホモロジーの跡を計算する。
- 局所 Ext と大域 Ext を区別する。
- 双対化層と Serre 双対の対を構成する。
- 曲線と曲面のコホモロジーに双対性を適用する。

## 射影空間上の双対性

$\mathbb P^n$ の $\mathcal O(d)$ の大域切断の単項式に、ちょうど一つの Laurent 単項式を掛けると、$\mathcal O(-n-1)$ の最上次の基底になる。この掛け算が双対性の出発点である。以下、$k$ は体とする。

:::proposition[射影空間の跡と対]{#prop-ag-projective-trace}
$P=\mathbb P^n_k$、$\omega_P=\mathcal O(-n-1)$ とする。$H^n(P,\omega_P)\to k$ を $(X_0\cdots X_n)^{-1}$ の係数を取る写像として定める。積とこの跡で得る $H^i(\mathcal O(d))\times H^{n-i}(\mathcal O(-d)\otimes\omega_P)\to k$ は完全対である。
:::
:::proof
$n>0$ とする。中間次数では両群が零である。零次の非負指数列 $a$ に対応する最上次の指数列は $-1-a$ である。掛け合わせて跡が $1$ となるのは指数の和が全て $-1$ のとき、すなわちこの組だけである。他の組では最上次の商における跡が零なので、二基底の行列は単位行列になる。$i=n$ の場合は二因子を交換すればよい。$n=0$ では空間は一点で、局所基底 $X_0^{-1}$ と通常の積で同じ結論を得る。
:::

### 跡の係数を実際に取り出す

:::example[射影平面の対の行列]
$\mathcal O(1)$ の基底を $X,Y,Z$ とすると、$\mathcal O(-4)$ の最高次基底は $X^{-2}Y^{-1}Z^{-1}$、$X^{-1}Y^{-2}Z^{-1}$、$X^{-1}Y^{-1}Z^{-2}$ である。第一の切断 $X$ を第二の Laurent 単項式と掛けると $X^0Y^{-2}Z^{-1}$ になる。$X$ の指数は非負なので、これは最高次では境界となる。一方第一の Laurent 単項式と掛けると $(XYZ)^{-1}$ なので跡は $1$ である。同様にして対の行列は単位行列となる。
:::

:::example[両側に複数の項がある場合]
$\mathbb P^1$ で $f=2X^2+3XY-Y^2$ と $\gamma=X^{-3}Y^{-1}+4X^{-2}Y^{-2}+5X^{-1}Y^{-3}$ を取る。跡に寄与するのは同じ位置の基底の積だけであり、値は $2+3\cdot4-5=9$ である。他の項には指数が非負の変数があるため、任意の標数でも係数を同じ式で求められる。ここでは $9$ は係数体の元として読む。
:::

### 積の跡と斉次式の線形写像を対応させる

双対対を使うと、Laurent 単項式の商で求めた乗法が、多項式の乗法の転置であると分かる。次数の上がる向きは同じでも、双対を取ると定義域と終域が逆になる。

:::proposition[ねじり層の乗法の随伴]
$P=\mathbb P^n_k$、$n\ge1$ とし、$f$ を次数 $e\ge0$ の斉次式とする。乗法
$$
H^n(P,\mathcal O(-n-1-b-e))
\xrightarrow{f}H^n(P,\mathcal O(-n-1-b))
$$
は、跡による完全対で $H^0(P,\mathcal O(b))\xrightarrow{f}H^0(P,\mathcal O(b+e))$ の双対に対応する。$b\ge0$ かつ $f\ne0$ なら最高次の写像は全射である。
:::

:::proof
$u$ を次数 $b$ の切断、$\gamma$ を次数 $-n-1-b-e$ の最高次の類とする。係数を掛ける写像はすべての局所化における制限と可換なので、Čech 複体の写像を与える。共有する最高次の跡を取れば
$$
\operatorname{tr}\bigl(u(f\gamma)\bigr)=\operatorname{tr}\bigl((fu)\gamma\bigr)
$$
である。この等式は Laurent 多項式の環の結合則であり、境界を掛けても境界になるので類に依存しない。:ref[prop-ag-projective-trace] の二つの完全対の下でこれが双対写像の定義である。$f\ne0$ なら斉次多項式環が整域であるため、切断側の乗法は単射である。有限次元ベクトル空間上でその双対は全射である。
:::

:::example[射影直線で乗法の転置を行列で読む]
標数零で $P=\mathbb P^1$、$f=X+2Y$、$b=1$ とする。切断側の基底をそれぞれ $X,Y$ と $X^2,XY,Y^2$ とすると、$fX=X^2+2XY$、$fY=XY+2Y^2$ なので行列は
$$
\begin{pmatrix}1&0\\2&1\\0&2\end{pmatrix}
$$
である。最高次の基底はこれらに双対な Laurent 指数であるため、対応する行列はその転置になる。例えば入力の係数ベクトル $(u,v,w)$ の像は $(u+2v,v+2w)$ である。核は $v=-2w,u=4w$ を満たすので $(4,-2,1)$ で生成される。最初の二列から行列式 $1$ の小行列を取れるため、全射性はこの計算からも分かる。
:::

:::example[最高次の核は超曲面の正則微分と対応する]
$\mathbb P^2$ の平面四次曲線 $C=V(f)$ を取る。$H^2(\mathcal O(-4))$ の基底は三つあり、$H^2(\mathcal O)=0$ なので $f$ 倍は全てを零に送る。構造層の完全列の連結写像は $H^1(\mathcal O_C)$ をこの三次元空間と同型にする。これに対し、随伴公式による滑らかな $C$ の $\omega_C=\mathcal O_C(1)$ の切断は三つの座標で張られる。この場合に対を計算するには、$H^1(\mathcal O_C)$ の類をまず連結写像で Laurent 基底へ移し、一次座標を掛けて係数を取り出せばよい。後の双対化層の構成は、この計算を高次元の閉埋め込みへ拡張する。
:::

### 射影直線の跡を微分の留数で表す

局所座標 $t$ を使うと、射影直線の最高次の類は重なりの Laurent 微分で表せる。$0,\infty$ 以外の点に極を作らないこの表示では、跡は一つの Laurent 係数を取る操作である。

:::proposition[射影直線での留数の表示]
$P=\mathbb P^1_k$ の有限開を $U_0$、無限遠側を $U_1$、座標を $t=X_1/X_0$ とする。$\omega_P=\Omega^1_{P/k}$ の Čech 一次類は、$k[t,t^{-1}]\,dt$ を二つの開から来る正則微分で割って得る。$dt/t$ の類を $1$ に送る正規化では、跡は $t^{-1}dt$ の係数であり、それは無限遠における留数の負に等しい。
:::

:::proof
$U_0$ の正則微分は $k[t]\,dt$ であり、$U_1$ では座標 $s=1/t$ で $k[s]\,ds$ となる。$ds=-t^{-2}dt$ なので後者は $t^{-2}k[t^{-1}]\,dt$ である。従って重なりでどちらからも来ない指数は $-1$ だけであり、商は $dt/t$ が基底である。$t^m dt=-s^{-m-2}ds$ だから、$s^{-1}ds$ を与えるのも $m=-1$ だけである。そこでの係数は $t^{-1}dt$ の係数の負である。この跡の正規化は、$U_0$ の $X_0^{-2}$ を $dt$ に対応させる $\mathcal O(-2)\simeq\omega_P$ の規約で :ref[prop-ag-projective-trace] と一致する。
:::

:::example[有理微分の複数の項から係数を分ける]
標数零で $\gamma=(4t^2+3t^{-1}-t^{-4})\,dt$ とすると、$4t^2dt$ は有限開から来て、$-t^{-4}dt$ は無限遠側で $s^2ds$ となる。従って類は $3[dt/t]$、跡は $3$ である。無限遠での式は $(-4s^{-4}-3s^{-1}+s^2)\,ds$ なので、その留数は $-3$ である。この二つの計算は重なりにおける類を調べており、一つの正則微分が $P$ 上に延びたと主張しているわけではない。
:::

## Ext 層と局所 Ext

:::definition[局所 Ext と大域 Ext]{#def-ag-sheaf-ext}
$\mathcal H om_{\mathcal O_X}(\mathcal F,-)$ の右導来関手を :term[局所 Ext]{en="sheaf Ext" yomi="きょくしょえくすと"} $\mathcal E xt^j_X$ といい、$\Hom_{\mathcal O_X}(\mathcal F,-)$ の右導来関手を $\operatorname{Ext}^j_X$ と書く。
:::
例えば $\mathcal F$ が有限階数局所自由なら $\mathcal H om(\mathcal F,-)=\mathcal F^\vee\otimes-$ は完全であり、局所 Ext の正次数は零である。大域 Ext はこの場合 $H^j(X,\mathcal F^\vee\otimes\mathcal G)$ であり、一般には
$$
H^p(X,\mathcal E xt^q_X(\mathcal F,\mathcal G))\ \Longrightarrow\
\operatorname{Ext}^{p+q}_X(\mathcal F,\mathcal G)
$$
を持つ。これは :ref[def-ag-sheaf-cohomology] の導来関手と :ref[thm-cat-grothendieck-spectral-sequence]{text="Grothendieck スペクトル系列"} の合成から従う。仮定の確認として、入射 $\mathcal I$ に対し $\mathcal H om(\mathcal F,\mathcal I)$ は脆弱である。開埋め込みからの 零延長 は茎で元の加群または零となる完全な関手であり、$V\subset U$ への制限は $\mathcal F$ の二つの 零延長 の単射に沿う写像の制限に等しい。従って入射性で全射になる。局所自由層の場合は系列に零行だけが残る。 加群層の入射分解は Abel 群の層としても脆弱分解となるので、この系列で取るコホモロジーは前章までの定義と一致する。

例えば滑らかな曲線の点 $p$ の剰余層は局所的に $0\to A\xrightarrow{t}A\to k(p)\to0$ の分解をもつ。双対を取ると $\mathcal E xt^1(k(p),\mathcal O)$ の茎は $k(p)$、$\mathcal H om(k(p),\mathcal O)$ の茎は零になる。

### 局所分解の双対から Ext を計算する

:::example[厚み付き点の Ext]
$A=k[t]_{(t)}$ と $M=A/(t^3)$ に対し、$0\to A\xrightarrow{t^3}A\to M\to0$ が自由分解である。$A$ は整域なので双対した二項複体の零次の核は零、一次の余核は $A/(t^3)$ である。従って $\operatorname{Ext}^1_A(M,A)$ は長さ三で、$1,t,t^2$ が $k$ 基底になる。支持が一点でも Ext には全ての厚みが残る。
:::

:::example[二つのパラメータを使う Koszul 双対]
$A=k[x,y]_{(x,y)}$ の剰余体の分解は $A\xrightarrow{(-y,x)^t}A^2\xrightarrow{(x,y)}A$ である。双対後には列ベクトル $(x,y)^t$ と行ベクトル $(-y,x)$ が現れる。$-ya+xb=0$ は $a=xc,\ b=yc$ と同値なので中間は完全であり、最高次の余核は $k$ である。曲線の一点の場合と比べると、非零 Ext の次数が余次元二へ移る。
:::

### 局所自由分解から Ext 層の茎を計算する

Ext 層と大域 Ext の区別は、まず局所で分解を双対化し、その後に大域切断を取るかどうかにある。連接層の局所表示では、茎に移してから双対を計算しても同じ結果になる。

:::proposition[局所自由分解による Ext 層]
$X$ を局所 Noether スキーム、$\mathcal F$ を連接層、$\mathcal G$ を $\mathcal O_X$ 加群とする。アフィン開で $\mathcal F$ の有限自由層による分解を取ると、双対複体 $\mathcal H om(\mathcal P_\bullet,\mathcal G)$ の第 $j$ コホモロジーは $\mathcal E xt^j_X(\mathcal F,\mathcal G)$ の制限に同型である。特に茎は $\operatorname{Ext}^j_{\mathcal O_{X,x}}(\mathcal F_x,\mathcal G_x)$ である。
:::

:::proof
アフィン上で $\mathcal F$ は有限生成加群に対応し、Noether 性で自由加群からの全射の核は再び有限生成である。従って有限階数の自由加群による分解を各段で作れ、層化は完全である。$\mathcal G$ の入射分解 $\mathcal I^\bullet$ を取り、$\mathcal H om(\mathcal P_i,\mathcal I^j)$ の二重複体を考える。入射層を開へ制限しても入射である。実際、開への制限の左随伴は開からの零延長であり、これは茎で完全なので、随伴により制限は入射性を保つ。各入射層について内部 Hom は第一変数の短完全列に対して完全である。これはさらに小さい開へ制限して Hom を評価しても入射性が残るためである。一方各 $\mathcal P_i$ は有限自由なので第二変数に対して内部 Hom が完全である。二方向に計算すれば、一方で定義による Ext 層、他方で $\mathcal H om(\mathcal P_\bullet,\mathcal G)$ のコホモロジーを得る。全次数ごとの項は有限個なので、二重複体からの比較が適用できる。

自由層の階数が有限であることから、内部 Hom の茎は対応する自由加群と $\mathcal G_x$ の Hom である。茎の完全性は複体のコホモロジーを茎へ移す操作と可換にする。自由分解を茎に移すと $\mathcal F_x$ の自由分解だから、加群の Ext を得る。
:::

:::example[二つの独立な厚みで Ext を求める]
$A=k[x,y]_{(x,y)}$、$B=A/(x^2,y^3)$ とする。$x^2$ は非零因子で、$A/(x^2)$ で $y^3$ は非零因子なので Koszul 分解
$$
0\to A\xrightarrow{(-y^3,x^2)^{\mathsf t}}A^2
\xrightarrow{(x^2,y^3)}A\to B\to0
$$
を使う。$A$ への双対では左の核が零で、中央の条件 $-y^3a+x^2b=0$ は $(a,b)=(x^2c,y^3c)$ となる。実際 $x^2$ と $y^3$ が互いに素であるために $x^2\mid a$ であり、約分すると $b=y^3c$ である。従って非零なのは二次だけで、余核は $B$ になる。$1,y,y^2,x,xy,xy^2$ が $k$ 基底であるから長さは六である。次数二という情報は余次元から、長さ六の情報は二つの厚みから来る。
:::

:::example[ベクトル束の局所 Ext が零でも大域 Ext が残る]
$X=\mathbb P^1_k$ で $\mathcal F=\mathcal O(3)$、$\mathcal G=\mathcal O$ とする。第一変数が局所自由なので $\mathcal E xt^1(\mathcal F,\mathcal G)=0$ である。一方 $\mathcal H om(\mathcal F,\mathcal G)=\mathcal O(-3)$ であり、局所から大域への系列の零行により
$$
\operatorname{Ext}^1_X(\mathcal O(3),\mathcal O)
=H^1(X,\mathcal O(-3))\simeq k^2
$$
である。標準二開での基底は $t^{-1},t^{-2}$ である。例えば第一の類は重なりで $\mathcal O(3)$ の持ち上げに $\mathcal O$ の成分を一つ加えることで拡大を作る。各開で拡大が分裂していても、二つの分裂の差が余境界として消えなければ大域では分裂しない。
:::

### 拡大の Čech 表示と連結写像をつなぐ

局所から大域への系列の低次数では、$\mathcal H=\mathcal H om_X(\mathcal F,\mathcal G)$ とおくと
$$
0\to H^1(X,\mathcal H)\to\operatorname{Ext}_X^1(\mathcal F,\mathcal G)
\to H^0(X,\mathcal E xt_X^1(\mathcal F,\mathcal G))
\xrightarrow{d_2}H^2(X,\mathcal H)
$$
が完全である。全次数一のフィルトレーションの一方の商は $(1,0)$ 項であり、そこへ出入りする二頁以降の微分はない。他方の商は $(0,1)$ 項のうち $(2,0)$ への $d_2$ の核である。ここに入る微分はなく、より高い頁では終域の行が負になる。これが完全性を与える。中央の写像は拡大類を各開に制限する写像であり、その核は局所的に分裂する拡大である。局所分裂の差は $\mathcal H$ の一次余閉鎖となり、分裂の選択を替えると余境界だけ変わる。逆に局所直和 $\mathcal G\oplus\mathcal F$ をこの余閉鎖の非対角成分で貼れば拡大が復元される。

二開で自明な二本の直線束を貼り合わせると、拡大の類は移行行列の非対角成分として表せる。この成分を変えても、各開で持ち上げの基底を変更するだけなら同じ拡大になる。

:::example[自己双対な Euler 列を非自明な拡大として見る]
$\mathbb P^1$ で $t=Y/X$ とし、$\mathcal O(-2)$ の基底を $e_0,e_1=t^{-2}e_0$ とする。商が $\mathcal O$ の拡大の持ち上げを $f_i$ として、$f_1=f_0+t^{-1}e_0$ と貼る。この変換の行列式は $t^{-2}$ なので、中央は階数二の局所自由層である。拡大類は $t^{-1}\in H^1(\mathcal O(-2))$ であり、非零である。実際、定数切断 $1$ の境界は二つの持ち上げの差 $t^{-1}e_0$ なので、もし大域的な商への切断があればこの類が零となり矛盾する。

$\mathbb P^1$ の Euler 列 $0\to\mathcal O(-2)\to\mathcal O(-1)^2\to\mathcal O\to0$ も非自明である。中央の大域切断が零であるため $\mathcal O$ の定数 $1$ は持ち上がらないからである。拡大群が一次元なので、二つの非零拡大は端の束の一方を非零定数倍して同型になる。従ってこの貼り合わせの中央も $\mathcal O(-1)^2$ である。局所で直和に見えるという事実から $\mathcal O(-2)\oplus\mathcal O$ とは判定できない。
:::

## 双対化層

:::definition[双対化層]{#def-ag-dualizing-sheaf}
純次元 $d$ の固有 $k$ スキームで、連接層 $\mathcal F$ に対する $\Hom_X(\mathcal F,\omega)$ と $H^d(X,\mathcal F)^\vee$ が自然に同型になる連接層 $\omega$ を :term[双対化層]{en="dualizing sheaf" yomi="そうついかそう"} という。
:::
$\mathcal F=\omega$ とすると恒等射から $H^d(X,\omega)\to k$ が得られ、これが跡である。自然性により一般の同型は $H^d(\mathcal F)\to H^d(\omega)\to k$ で与えられる。Yoneda の議論と同じく、自然同型を恒等射で評価すると双対化層は跡とともに一意になる。

### 零次元の双対化層を線形代数で作る

:::proposition[有限スキームの双対化加群]
$B$ が体 $k$ 上有限次元の可換代数なら、$\operatorname{Spec}B$ の双対化層は $B^\vee=\Hom_k(B,k)$ に対応し、跡は $1$ での評価である。
:::

:::proof
$\lambda\in B^\vee$ に対し $(b\lambda)(x)=\lambda(bx)$ として $B$ 加群にする。有限 $B$ 加群 $M$ から $B^\vee$ への $B$ 線形写像 $\phi$ を $m\mapsto\phi(m)(1)$ に送る。逆に $k$ 線形形式 $\ell$ から $\phi_\ell(m)(b)=\ell(bm)$ を作ると、$B$ 線形であり、二操作は互いに逆で $M$ に関して自然である。アフィン上の $H^0$ は $M$ なので双対化層の普遍性を満たす。
:::

:::example[双対化加群が自由となる例]
$B=k[t]/(t^3)$ で $\lambda$ を $t^2$ の係数を取る線形形式とする。基底 $1,t,t^2$ に対し $(a,b)\mapsto\lambda(ab)$ の行列は反対角に $1$ が並ぶ。これは可逆なので $B\to B^\vee,\ a\mapsto a\lambda$ は同型である。この同型では大域切断の跡は $t^2$ の係数となり、通常の乗法写像の行列の跡とは異なる。
:::

:::example[双対化加群が自由でない零次元環]
$B=k[x,y]/(x,y)^2$ では極大イデアルに殺される元の空間は $kx+ky$ の二次元である。$B^\vee$ で極大イデアルに殺される形式は $(x,y)$ を消すものに限り一次元である。従って $B$ と $B^\vee$ は $B$ 加群として同型でなく、この局所環の双対化加群は可逆でない。
:::

### 双対化層の一意性と跡の選び方

:::proposition[跡を保つ一意な同型]
二つの対 $(\omega,\operatorname{tr})$、$(\omega',\operatorname{tr}')$ が双対化層の普遍性を満たすなら、$\operatorname{tr}'\circ H^d(u)=\operatorname{tr}$ を満たす同型 $u:\omega\to\omega'$ は一意である。
:::

:::proof
$\omega'$ の普遍性を $\mathcal F=\omega$ に適用して、線形形式 $\operatorname{tr}$ に対応する $u$ を得る。逆に $\omega$ の普遍性を $\omega'$ に適用して $v$ を得る。合成について $\operatorname{tr}\circ H^d(vu)=\operatorname{tr}'\circ H^d(u)=\operatorname{tr}$ なので、普遍性の単射性で $vu=\id_\omega$ である。反対も同じであり、$u$ は一意な同型である。
:::

:::example[跡を定数倍すると同型も変わる]
$\operatorname{tr}'=c\,\operatorname{tr}$（$c\in k^\times$）を同じ $\omega$ 上に置けば、両者とも普遍性を持つ。上の $u$ は $c^{-1}\id_\omega$ である。例えば $c=2$ なら新しい跡のもとで旧来の値を得るには切断を $1/2$ 倍する。跡を明示した双対性と、単に両ベクトル空間の次元が等しいという主張の差が、この定数に現れる。
:::

### Artin 局所環で双対化加群が可逆になる条件

有限な局所環では、双対化加群が自由かどうかを積の対の非退化性で調べられる。その判定には環全体の乗法行列を使う必要はなく、極大イデアルが殺す最も深い部分の次元を見ればよい。

:::definition[Artin 局所環のソクル]
Artin 局所環 $(B,\mathfrak m,K)$ の :term[ソクル]{en="socle" yomi="そくる"} は $\operatorname{Soc}(B)=\{b\in B:\mathfrak m b=0\}$ である。$\mathfrak m$ が作用を殺すので、これは剰余体 $K=B/\mathfrak m$ 上のベクトル空間である。
:::

:::proposition[零次元の Gorenstein 性の判定]
$B$ を体 $k$ 上有限次元の局所代数とする。次の三条件は同値である。
- $B^\vee=\Hom_k(B,k)$ が自由階数一である。
- $\operatorname{Soc}(B)$ が剰余体上一次元である。
- ある $k$ 線形形式 $\lambda:B\to k$ に対して $(a,b)\mapsto\lambda(ab)$ が非退化である。
:::

:::proof
$B^\vee$ のソクルは $\Hom_k(K,k)$ に等しい。実際、$b\in\mathfrak m$ による $\lambda$ の積が零という条件は、$\lambda$ が全ての $bx$ を消すという条件である。$x$ に $1$ を取れば $\lambda(\mathfrak m)=0$ であり、逆も成り立つ。$\Hom_k(K,k)$ は $K$ 加群であり、$k$ 次元が $[K:k]$ なので $K$ 次元は一である。従って $B^\vee\simeq B$ なら $\operatorname{Soc}(B)$ の次元は一である。

逆にソクルが $K$ 上一次元とする。$k$ 線形形式 $\lambda$ を、ソクル上で恒等的に零にならないように選ぶ。例えばソクルの非零元を $k$ 基底の一つへ延長し、その係数を取れる。$0\ne a\in B$ に対し、任意の非零イデアルがソクルに交わることを使う。実際 $\mathfrak m$ は冪零であるため、$Ba$ に $\mathfrak m$ の冪を掛け、非零である最後の段を取ると $\mathfrak m$ に殺される非零部分を得る。ソクルは $K$ 上一次元なので、この $B$ 部分加群はソクルの全体に等しい。従って $Ba$ はソクルを含み、ある $b$ について $\lambda(ab)\ne0$ である。これで対の左核が零となり、有限次元と可換性により非退化である。最後に非退化な対は $a\mapsto(b\mapsto\lambda(ab))$ による同じ次元の $B$ 加群の単射を与えるので、$B^\vee$ は自由階数一である。
:::

:::example[二本の方程式の交叉で双対の基底を作る]
$B=k[x,y]/(x^2,y^3)$ のソクルは $kxy^2$ である。実際、$x$ に殺される元は $x,xy,xy^2$ の張る部分空間で、さらに $y$ に殺される元は最後の一項に限る。$\lambda$ を $xy^2$ の係数とすれば、基底の $1,y,y^2,x,xy,xy^2$ は対によって反対順の $xy^2,xy,x,y^2,y,1$ と組になる。例えば $y$ の双対は $xy$、$x$ の双対は $y^2$ である。積は他の指数の組なら $xy^2$ の係数を与えないので、双対表を直接検算できる。従って双対化加群は $B\lambda$ であり、上で局所 Koszul 分解から得た自由な Ext 加群と整合する。
:::

:::example[非分離な体拡大では通常の体の跡が零でも双対はある]
標数 $p>0$ で $B=k(t)[u]/(u^p-t)$ とする。$t$ は $k(t)$ の $p$ 乗でないので $B$ は次数 $p$ の純非分離拡大である。基礎体を $K=k(t)$ と書く。任意の $b\in B$ に対する $B$ 上の乗法の跡 $\operatorname{Tr}_{B/K}(b)$ は零である。実際、代数閉包上では $B$ は一つの根 $\alpha$ の周囲の $K^{\mathrm{alg}}[\epsilon]/(\epsilon^p)$ になり、$b$ の乗法は $\epsilon$ の冪の基底で対角がすべてその定数項に等しい三角行列になる。対角和は標数 $p$ で零である。一方 $B^\vee$ は体 $B$ 上一次元で、双対化の跡は $B^\vee\to K$ の $1$ での評価である。例えば $u^{p-1}$ の係数を取る $\lambda$ は非零であり、$a\mapsto a\lambda$ は同型である。この同型の下での双対化の跡は $\lambda$ となり、代数の乗法行列の跡とは異なる。
:::

### 有限射に沿って双対化層を移す

有限射では、基底上で加群を忘れても有限性が保たれる。逆に、基底上の線形形式に始域の環の乗法を作用させると、始域の加群に値を取る双対を作れる。零次元の証明で使った操作を、各アフィン開で行って貼り合わせる。

:::proposition[有限射の双対化層]
$f:X\to Y$ を体上の純次元 $d$ の固有スキーム間の有限射とし、$Y$ が双対化層 $\omega_Y$ をもつとする。$f_*\mathcal O_X$ の作用を
$$
(b\lambda)(a)=\lambda(ba)
$$
で与えた層 $\mathcal H om_Y(f_*\mathcal O_X,\omega_Y)$ に対応する $X$ 上の層を $\omega_X$ とすると、これは双対化層である。その跡は
$$
H^d(X,\omega_X)\longrightarrow H^d(Y,\omega_Y)\xrightarrow{\operatorname{tr}_Y}k
$$
であり、第一の写像は基底側の Hom の $1$ での評価から来る。
:::

:::proof
アフィン $V=\operatorname{Spec}A\subset Y$ で $f^{-1}V=\operatorname{Spec}B$、$\omega_Y|_V=\widetilde N$ とする。$B$ は有限 $A$ 加群で、Noether 性から有限表示である。従って $\Hom_A(B,N)$ の局所化は局所化後の Hom に等しい。実際、有限表示から Hom を二つの有限直和の間の写像の核で表せ、局所化はこの核を保つ。従って局所の $B$ 加群 $\Hom_A(B,N)$ は貼り合い、連接層 $\omega_X$ を定める。連接性は $\Hom_A(B,N)$ が有限 $A$ 加群の部分加群であり、$A$ 上有限なら $B$ 上も有限であることから従う。

$B$ 加群 $M$ について、$B$ 線形 $\phi:M\to\Hom_A(B,N)$ を $m\mapsto\phi(m)(1)$ へ送る。逆は $A$ 線形 $\ell:M\to N$ から、$\phi_\ell(m)(b)=\ell(bm)$ を作る写像である。$B$ 線形性は $(b'\phi_\ell(m))(b)=\ell(b' bm)$ で確かめられる。$1$ で評価すれば $\ell$ が戻り、逆の合成も $\phi(bm)(1)=(b\phi(m))(1)=\phi(m)(b)$ から恒等である。これらは制限と可換なので、
$$
\Hom_X(\mathcal F,\omega_X)\simeq\Hom_Y(f_*\mathcal F,\omega_Y)
$$
が得られる。後者は $H^d(Y,f_*\mathcal F)^\vee$ であり、有限射はアフィンなので準連接層の順像のコホモロジーの比較により $H^d(X,\mathcal F)^\vee$ と同型である。全ての写像は $\mathcal F$ に関して自然であるため、双対化層の条件を満たす。$\mathcal F=\omega_X$ の恒等射が基底での評価写像に対応することから跡の式を得る。
:::

:::example[二重被覆の局所双対で分岐点を含めて計算する]
$B=A[u]/(u^2-a)$ とし、$1,u$ を $A$ 基底とする。線形形式 $\lambda$ を $u$ の係数とすれば、対 $\lambda(bc)$ の行列は
$$
\begin{pmatrix}0&1\\1&0\end{pmatrix}
$$
である。実際 $u^2=a$ は定数項なので $\lambda(u^2)=0$ となる。この行列式は $-1$ なので、$a$ が零になっても $\Hom_A(B,A)\simeq B\lambda$ である。これに対し標数が二でない場合、通常の乗法の跡 $\operatorname{Tr}(b)$ は $b=b_0+b_1u$ について $2b_0$ であり、その積の対の行列は $\operatorname{diag}(2,2a)$ である。$a$ が零になる点では後者の対は退化するが、有限射における Hom の記述は変わらない。
:::

## Serre 双対性定理

:::theorem[射影空間上の Serre 双対性]{#thm-ag-serre-ambient}
$P=\mathbb P^n_k$ とする。任意の連接層 $\mathcal F$ と $0\le i\le n$ に対し、Yoneda 積の跡は自然な同型
$$
\operatorname{Ext}_P^{n-i}(\mathcal F,\omega_P)\simeq H^i(P,\mathcal F)^\vee
$$
を与える。
:::
:::proof
$H^i(P,\mathcal F)=\operatorname{Ext}_P^i(\mathcal O_P,\mathcal F)$ なので Yoneda 積は $H^n(P,\omega_P)$ に値を取る。:ref[thm-ag-serre-finiteness] によりコホモロジーは有限次元である。局所 Ext は有限自由層で計算して連接であるため、局所から大域への系列と同じ有限性定理から、大域 Ext も各次数で有限次元である。以下 $n>0$ とする。$n=0$ の場合は一点上の有限次元ベクトル空間の双対性である。有限個の $\mathcal O(d)$ の直和では :ref[prop-ag-projective-trace] と局所自由層の Ext の計算が完全性を与える。

$\mathcal F$ はそのような直和による有限分解をもつ。実際、:ref[lem-ag-extension-twist] から直和 $\mathcal E_0\to\mathcal F$ とその核への直和 $\mathcal E_1$ からの全射を得る。$\mathcal E_1\to\mathcal E_0$ の各成分はねじり層の大域切断、すなわち斉次多項式なので、対応する次数付き自由 $S=k[X_0,\dots,X_n]$ 加群間の写像から誘導される。その余核 $M$ を取ると $\widetilde M=\mathcal F$ である。次数付き加群 $M$ の最小自由分解は、各段階で $S_+$ を法とした空間の斉次基底を持ち上げて作る。持ち上げが生成することは次数の下方有界性と次数に関する帰納法から従う。最小なので微分の成分は $S_+$ に属し、$k$ とテンソルすると全微分が零である。他方、$k$ の Koszul 分解は長さ $n+1$ である。これは変数を順に加え、前段階の商では次の変数が非零因子なので、二項複体とのテンソル写像の錐が次の商の分解になることから分かる。従って $\operatorname{Tor}_j^S(M,k)$ は $j>n+1$ で零であり、最小分解もその次数で零になる。層化は完全だから所望の有限分解を得る。

この分解を $0\to\mathcal K\to\mathcal E\to\mathcal F\to0$ の形の短完全列に分ける。対応する一次の拡大類を $e$ とすると、$u\in H^i(P,\mathcal F)$ の境界は $e\circ u$ である。一方 $v\in\operatorname{Ext}^{n-i-1}_P(\mathcal K,\omega_P)$ の境界は $v\circ e$ である。Yoneda 積の結合則から、両者は $\operatorname{tr}(v\circ e\circ u)$ で対応する。従って境界を含め、反変 Ext の長完全列からコホモロジーの長完全列の線形双対への可換図式が得られる。

負次数のコホモロジーと Ext を零と約束すれば、この比較写像を全ての整数 $i$ に対して定められる。ねじり層では $0\le i\le n$ の外でも両辺は零である。分解の末端から、比較図式の $\mathcal E$ と $\mathcal K$ の全ての項で同型が既知であるとして五項補題を適用すれば、$\mathcal F$ に対する全ての比較写像が同型になる。特に $\operatorname{Ext}_P^j(\mathcal F,\omega_P)$ は $j>n$ で零である。
:::

### 大域 Ext とコホモロジーが一致する次数

:::example[射影平面の点層]
$P=\mathbb P^2$ の有理点 $p$ に対し、$H^0(k(p))=k$、高次群は零である。従って $\operatorname{Ext}_P^2(k(p),\omega_P)=k$、他の次数は零となる。局所計算でも、点の余次元二の Koszul 分解から $\mathcal E xt^2$ が一点支持の層として得られ、その $H^0$ だけが残る。
:::

:::example[負のねじりで双対の側に切断が移る]
$\mathbb P^3$ の $\mathcal F=\mathcal O(-6)$ では $H^3$ の双対は $\Hom(\mathcal O(-6),\mathcal O(-4))=H^0(\mathcal O(2))$ である。四変数二次式には平方四つと相異なる変数の積六つがあるので次元は $10$ である。実際、最高次の負指数の絶対値は四つの正整数で和が $6$ となり、各絶対値から $1$ を引くと同じ二次単項式を数えることになる。
:::

### 点に支持された層で双対性の次数を確かめる

点の剰余層は零次の大域切断しか持たない。それでも双対の側では、周囲に埋め込まれた余次元だけ自由分解を通る必要がある。このとき局所 Ext の最高次と大域 Ext の最高次が一致することを確認できる。

:::example[厚み六の点の双対を射影平面で求める]
$P=\mathbb P^2_k$ の $Z=V(Y^2,Z_0^3)$ を取る。ここで座標を $[X:Y:Z_0]$ とする。台は $[1:0:0]$ の一点であり、$X=1$ の環は $B=k[y,z]/(y^2,z^3)$ である。$H^0(\mathcal O_Z)=B$、高次は零である。周囲の切断で Koszul 分解を作ると
$$
0\to\mathcal O_P(-5)\to\mathcal O_P(-2)\oplus\mathcal O_P(-3)
\to\mathcal O_P\to\mathcal O_Z\to0
$$
となる。$\omega_P=\mathcal O_P(-3)$ に対する局所双対複体では最高次が $\mathcal O_P(2)$ で、そこへの像は $Y^2,Z_0^3$ を掛けて得る部分層である。従って $\mathcal E xt_P^2(\mathcal O_Z,\omega_P)=\mathcal O_Z(2)$ となり、$X$ が単元である近傍で基底を使うと $B$ に同型である。中間の局所 Ext は正則列の双対で消える。局所から大域への系列には二次の行と零次の列の交点だけが残り、$\operatorname{Ext}_P^2\simeq B$ が得られる。双対化層の一意性によって、この $B$ 加群を $B^\vee$ と同一視するときの跡は、基底を調整すれば $yz^2$ の係数と取れる。
:::

### 直和でないベクトル束の双対性を検算する

ねじり層の完全対を用いる証明では、一般の層を短完全列の反復で処理した。Euler 列は、短完全列で現れる境界を計算して双対性を検算できる具体例である。

:::example[射影平面の余接層の一次とその双対]
$P=\mathbb P^2_k$ とする。Euler 列
$$
0\to\Omega^1_{P/k}\to\mathcal O(-1)^3\to\mathcal O\to0
$$
において、中央の全コホモロジーは零であり、右の大域切断は $k$ である。従って $H^1(\Omega^1_{P/k})\simeq k$ で、他の群は零である。この列を $\omega_P$ に対して双対化すると、局所分裂性から
$$
0\to\mathcal O(-3)\to\mathcal O(-2)^3\to
(\Omega^1_{P/k})^\vee\otimes\omega_P\to0
$$
を得る。中央のコホモロジーはすべて零で、左の最高次 $H^2(\mathcal O(-3))$ は $k$ である。従って右の $H^1$ は $k$ となり、ほかの群は零になる。第一の列では定数 $1$ が一次へ上がり、第二の列では一次が跡のある最高次へ上がる。局所の行列式を取ると階数二の束 $\Omega^1$ には $\det\Omega^1=\omega_P$ があるため、外積による $\Omega^1\simeq(\Omega^1)^\vee\otimes\omega_P$ も得られる。これは二つの計算が同じ一次の群に戻る理由である。
:::

:::example[超曲面の切断の次元と大域 Ext を同じ商から得る]
$P=\mathbb P^3_k$ に次数 $e=5$ の超曲面 $X$ を取り、$\mathcal F=i_*\mathcal O_X(1)$ とする。分解は $0\to\mathcal O(-4)\to\mathcal O(1)\to\mathcal F\to0$ である。これを $\omega_P=\mathcal O(-4)$ に対して双対化すると、局所 Ext 一次は $\mathcal O_X$、他の局所 Ext は零になる。従って $\operatorname{Ext}_P^1(\mathcal F,\omega_P)=H^0(X,\mathcal O_X)=k$ である。元の完全列では $\mathcal O(-4)$ の最高次に一つの基底があり、$\mathcal O(1)$ は正次数のコホモロジーを持たないため、$H^2(X,\mathcal O_X(1))\simeq k$ となる。互いに双対になる二群を作るときに、方程式の次数五がねじり一と周囲の標準層の次数 $-4$ を相殺している。
:::

## Cohen–Macaulay 性と双対性

:::theorem[Cohen–Macaulay 射影スキームの双対性]{#thm-ag-serre-duality}
$X\hookrightarrow P=\mathbb P^n_k$ を純次元 $d$ の Cohen–Macaulay 射影スキームとし、$c=n-d$ とする。
$$
i_*\omega_X=\mathcal E xt_P^c(i_*\mathcal O_X,\omega_P)
$$
で定めると $\omega_X$ は双対化層である。連接層 $\mathcal F$ と $0\le j\le d$ に対して $\operatorname{Ext}_X^{d-j}(\mathcal F,\omega_X)\simeq H^j(X,\mathcal F)^\vee$ である。$\mathcal F$ が局所自由なら左辺は $H^{d-j}(X,\mathcal F^\vee\otimes\omega_X)$ である。
:::
:::proof
点で正則局所環 $R=\mathcal O_{P,x}$ と $B=\mathcal O_{X,x}$ を取る。Cohen–Macaulay 性とは $\operatorname{depth}B=\dim B$ である。純次元の仮定と有限型体上の次元公式で $\dim R-\dim B=c$ である。$B$ は射影空間上で得た有限ねじり分解を茎へ移すと有限自由分解をもつ。従って :ref[thm-comm-auslander-buchsbaum]{text="Auslander–Buchsbaum の公式"}から $\operatorname{pd}_R B=c$ であり、$\operatorname{grade}(\ker(R\to B),R)=c$ である。この二つから Ext は $c$ より上と下で零である。$\omega_P$ は局所自由階数一なので同じ消滅が成り立つ。次数 $c$ の Ext は $B$ 加群である。作用は第一変数 $B$ から誘導でき、第二変数からの作用と一致するため、$\ker(R\to B)$ はその Ext を殺す。

入射 $\mathcal I$ に $\mathcal H om_P(i_*\mathcal O_X,\mathcal I)$ を適用した層は $X$ 上の入射層になる。実際、これは完全な $i_*$ に右随伴である。その複体に対して $\Hom_X(\mathcal F,-)$ の導来関手を取ると、合成は $\operatorname{Ext}_P^\bullet(i_*\mathcal F,\omega_P)$ を計算する。一方、その複体のコホモロジーは上の局所消滅により次数 $c$ の $\omega_X$ だけなので、系列は一行に退化して $\operatorname{Ext}_X^{d-j}(\mathcal F,\omega_X)$ を得る。:ref[thm-ag-serre-ambient] と閉埋め込みの順像がコホモロジーを保つことを合わせれば結論となる。$j=d$ は双対化層の定義に一致する。
:::
滑らかな場合、閉埋め込みは局所的に正則列で定まり、Koszul 分解の最高項を双対化すると $\omega_X\simeq\omega_P|_X\otimes\det(\mathcal I/\mathcal I^2)^\vee$ となる。正則列を替えた際の行列式が余法層の行列式の逆であるため貼り合う。余法短完全列の行列式を取るとこれは :ref[def-ag-canonical-sheaf] に一致する。滑らかでない場合にも双対化層は定義できる。例えば正則列一つで定まる任意の超曲面なら二項分解により $\omega_X=\mathcal O_X(e-n-1)$ であり、可逆である。

:::example[埋め込み点で Ext の次数が増える例]
$P=\mathbb P^3$ の平面 $H=V_+(\ell)$ 上に $k$ 有理点 $p$ を取り、$\mathcal I_X=\mathcal I_H\mathcal I_p$ で閉部分スキームを定める。$\mathcal O_X\to\mathcal O_H$ の核は $\mathcal I_H/\mathcal I_H\mathcal I_p\simeq\mathcal O_P(-1)\otimes k(p)$ である。$\mathcal I_H$ は可逆なので、$p$ で基底を選べばこの核を $k(p)$ と同一視できる。$p$ では零次元の部分加群が構造層に残る。局所環の極大イデアルがこの部分加群を殺すため、深さは零であり、次元二と一致しない。

平面の二項自由分解から $\mathcal E xt_P^j(\mathcal O_H,\omega_P)$ は $j=1$ にだけ現れる。点 $p$ では三つの正則パラメータの Koszul 分解があり、その双対は零次と最上次を反転した同じ複体になる。従って $\mathcal E xt_P^j(k(p),\omega_P)$ は $j=3$ にだけ現れ、そこでは一次元の剰余層である。上の短完全列に局所 Ext を適用すると、$\mathcal O_X$ の Ext は次数 $1$ と次数 $3$ の両方で非零である。Cohen–Macaulay の仮定を外すと一行への退化が失われることが、この例で具体的に分かる。
:::

:::proposition[双対化層の深さ]
:ref[thm-ag-serre-duality] の仮定のもとで、各 $x\in X$ について $\operatorname{depth}_{\mathcal O_{X,x}}\omega_{X,x}=\dim\mathcal O_{X,x}$ である。
:::

:::proof
$R=\mathcal O_{P,x}$、$B=\mathcal O_{X,x}$ とし、長さ $c$ の $B$ の最小 $R$ 自由分解を取る。これを $R$ に対して双対化すると、局所 Ext の消滅により $c$ 次だけにコホモロジー $\omega_{X,x}$ をもつ。ここでは $\omega_{P,x}$ の基底を一つ選んだ。最終の自由加群は非零で、そこへ入る写像の成分は極大イデアルに属す。従って Nakayama の補題でその余核 $\omega_{X,x}$ は非零である。双対複体を逆の順に読むと、この余核の長さ $c$ の自由分解であり、同じ理由で最小なので射影次元は $c$ である。Auslander–Buchsbaum の公式により $\omega_{X,x}$ の $R$ 上での深さは $\dim R-c=\dim B$ である。$R$ の作用は $B$ を経由するため、正則列の作用を $B$ で見ても同じであり、$B$ 上の深さもこの値になる。
:::

### 超曲面と正則列で双対化層を計算する

:::example[特異平面三次にも自明な双対化層がある]
尖点三次 $Y^2Z=X^3$ の構造層は $\mathcal O(-3)\to\mathcal O$ の余核である。この分解を $\omega_{\mathbb P^2}=\mathcal O(-3)$ に対して双対化すると、$\mathcal O(-3)\to\mathcal O$ が再び現れ、$\omega_C\simeq\mathcal O_C$ を得る。正規化は射影直線で双対化層が $\mathcal O(-2)$ である。正規化前後の双対化層が同じになるとは限らないことを、この二項計算で確認できる。
:::

:::example[二式の Koszul 分解の最高次数]
$\mathbb P^4$ 内の次数 $2,3$ の完全交叉曲面を取る。分解の最左項は $\mathcal O(-5)$ である。$\omega_{\mathbb P^4}=\mathcal O(-5)$ への双対の最右項は $\mathcal O$ で、そこに入る写像の像は二方程式で生成される。従って $\mathcal E xt^2$ は $\mathcal O_X$ である。これは局所的に正則列であることを使うため、同じ次数の二式でも共通成分があればこの分解は使えない。
:::

### 正則列を取り替えたときの双対化層の移行係数

正則列が一つの開で選ばれていると、最高次の Ext の基底を一つ書ける。基底が開の重なりでどう変わるかを追うと、双対化層がなぜ余法層の行列式の逆を含むかが分かる。

:::proposition[正則埋め込みでの局所 Ext]
$i:X\hookrightarrow P$ のイデアル層 $\mathcal I$ が局所的に長さ $c$ の正則列で生成され、$\mathcal L$ が $P$ 上の可逆層であるとする。このとき $\mathcal E xt_P^j(i_*\mathcal O_X,\mathcal L)$ は $j\ne c$ で零であり、$j=c$ では
$$
i_*\bigl(\mathcal L|_X\otimes\det(\mathcal I/\mathcal I^2)^\vee\bigr)
$$
と同型である。
:::

:::proof
アフィン上で $R$ を周囲の環、$I=(f_1,\dots,f_c)$ とする。$K$ を $R^c=F$ の外積からなる Koszul 複体とし、$\bigwedge^jF$ から $\bigwedge^{j-1}F$ への微分を $f=(f_1,\dots,f_c)$ で縮約する写像にする。正則列であるため $K$ は $R/I$ の分解である。可逆層を自明化し $\Hom_R(K,R)$ を取る。外積による完全対
$$
\bigwedge^jF\otimes\bigwedge^{c-j}F\longrightarrow\det F
$$
を使うと、双対の第 $j$ 項は $\bigwedge^{c-j}F\otimes(\det F)^\vee$ と同一視できる。縮約の規則 $d(\alpha\wedge\beta)=d\alpha\wedge\beta+(-1)^{\deg\alpha}\alpha\wedge d\beta$ から、双対微分はこの同一視で元の微分に符号を掛けたものになる。符号は各次数の基底変更で調整できるため、非零のコホモロジーは元の零次に対応する $j=c$ にだけあり、そこでは $(R/I)\otimes(\det F)^\vee$ となる。

$F\otimes_R R/I\to I/I^2$ が同型であることも Koszul 複体から確認できる。全射性は生成元の定義による。$\sum_i a_i f_i\in I^2$ なら、$I^2$ を $f_i$ で展開して各係数を $I$ に入れられる。差を引くと $f_i$ の関係式になるが、その係数ベクトルは Koszul の二次の像なので全成分が $I$ にある。従って元の $a_i$ もすべて $I$ に入り、単射性を得る。

別の正則列を $g_i=\sum_j a_{ij}f_j$ と表すと、新しい Koszul 複体から古い複体へ、第 $1$ 項で $a_{ij}$ を送る写像があり、外積すれば全次数の複体写像になる。第 $0$ 項は恒等なので、同じ商の恒等射を持ち上げている。第 $c$ 項の行列式は $\det(a_{ij})$ であり、$R/I$ 上では二つの余法基底を結ぶため単元である。最高次の Ext の基底は $\det F$ と双対であるため、局所基底の移行は余法の行列式の逆になる。可逆層 $\mathcal L$ の基底変更も同時に掛けて貼れば、主張の層の同型を得る。
:::

:::example[一次式を入れ替えた余法基底と最高次の符号]
$\mathbb A^2$ の原点の正則列 $(x,y)$ で最高次の余法基底を $[x]\wedge[y]$ とする。順序を $(y,x)$ に替えると余法の向きは $-1$ 倍になる。最高次の Ext も双対の基底を対応させるため $-1$ 倍になる。標数二ではこの二符号は同じである。さらに第一式を $x$ から $3x$ に替えると、標数が三でないとき余法の基底は $3$ 倍、対応する双対基底は $1/3$ 倍となる。局所 Ext の線形空間の次元はどの場合も一だが、周囲から誘導した基底を比較するとこの係数が現れる。
:::

### 最上次の双対化と全次数の双対性を区別する

Cohen–Macaulay の仮定は局所 Ext の二つ以上の行が現れないようにするために使う。最上次のコホモロジーだけを表すとき、余次元より下には最初から Ext がないので、上に別の行があってもその値は変わらない。

:::proposition[純次元の射影スキームの最上次の双対]
$i:X\hookrightarrow\mathbb P^n_k$ を純次元 $d$ の射影スキームとし、$c=n-d$ とする。このとき $\mathcal E xt^c_{\mathbb P^n}(i_*\mathcal O_X,\omega_{\mathbb P^n})$ は $X$ 上の双対化層の順像である。この主張には $X$ の Cohen–Macaulay 性を仮定しない。
:::

:::proof
周囲の正則局所環 $R$ で $B=R/I$ と書く。$X$ が純次元なので $I$ の最小素イデアルの高さは $c$ である。これは体上有限型の次元公式により、各既約成分を周囲のアフィン開に制限して求められる。この高さはその素イデアルを含む点の局所環に移しても同じである。$R$ は正則で、その各素イデアルにおける局所化も正則である。従って [イデアルに含まれる正則列の長さ](/math/commutative-algebra/depth-cohen-macaulay/) の判定から $\operatorname{grade}(I,R)=c$ であり、$\operatorname{Ext}_R^q(B,R)=0$（$q<c$）である。高い次数の消滅はここでは要求しない。

閉埋め込みの順像の右随伴 $\mathcal H om_P(i_*\mathcal O_X,-)$ と $\Hom_X(\mathcal F,-)$ を順に適用した系列の第二頁は
$$
\operatorname{Ext}_X^p\bigl(\mathcal F,\mathcal E xt_P^q(i_*\mathcal O_X,\omega_P)\bigr)
$$
である。ここでは $\mathcal E xt$ が $X$ 上の層になることを用いた。その加群作用は第一変数から来るので、周囲の環の作用のうち $I$ が零になる。全次数 $c$ であり得る項は $(p,q)=(0,c)$ だけである。二頁以降、そこへ入る微分は負の $p$ から来るため零であり、出る微分の終域は $q<c$ なので零である。従って極限でもこの項だけが残り、
$$
\operatorname{Hom}_X(\mathcal F,\omega_X)
\simeq\operatorname{Ext}_P^c(i_*\mathcal F,\omega_P)
\simeq H^d(X,\mathcal F)^\vee
$$
となる。最後は射影空間上の双対性である。自然性から双対化層の普遍性を満たす。
:::

:::example[最上次の双対化層は埋め込み点だけを見分けられない]
先の $\mathcal I_X=\ell\mathcal I_p$ により平面 $H$ に点の厚みを付けた例では、$\mathcal E xt_P^1(\mathcal O_X,\omega_P)\simeq\mathcal E xt_P^1(\mathcal O_H,\omega_P)=\mathcal O_H(-3)$ である。これは $0\to k(p)\to\mathcal O_X\to\mathcal O_H\to0$ の局所 Ext の列で、点の Ext が三次にしかないためである。従って $\omega_X$ は $H$ の双対化層を $X$ 上の層とみなしたものになる。

一方同じ短完全列の大域切断を取ると $H^0(\mathcal O_X)$ は次元二である。実際 $H^0(k(p))=k$、$H^0(\mathcal O_H)=k$ で、点の $H^1$ は零である。また $\omega_X$ の $H^2$ は射影平面上の $\mathcal O(-3)$ の $H^2$ なので次元一である。従って $H^2(\omega_X)$ は $H^0(\mathcal O_X)^\vee$ と同じ次元を持たない。双対化層自体は存在するが、Cohen–Macaulay の場合の全次数の双対性をこの層だけで表すことはできない。
:::

### 負のねじりの消滅と超平面切断の連結性

:::proposition[正規射影多様体の豊富な因子の連結性]{#prop-ag-ample-divisor-connectedness}
代数閉体上の正規整射影多様体 $X$ の次元を $d\ge2$ とする。$\mathcal O_X(D)$ が豊富な空でない有効 Cartier 因子 $D$ の台は連結である。とくに滑らかな既約射影多様体の一般の超平面切断は、次元が二以上なら滑らかかつ既約である。
:::

:::proof
$\mathcal O_X(D)$ の正冪が定める射影埋め込み $i:X\hookrightarrow P=\mathbb P^n$ を選び、この埋め込みのねじりを用いる。$c=n-d$ とする。点 $x$ で $e=\dim\mathcal O_{X,x}$ と置くと、次元公式により $\dim\mathcal O_{P,x}=c+e$ である。:ref[thm-comm-serre-normality] の $(S_2)$ により $\operatorname{depth}\mathcal O_{X,x}\ge\min(2,e)$ だから、Auslander–Buchsbaum の公式で $\mathcal O_{P,x}$ 上の射影次元は
$$
c+e-\operatorname{depth}\mathcal O_{X,x}\le c+\max(0,e-2)\le n-2
$$
となる。射影次元が有限なのは、射影空間上で構造層のねじりによる有限分解を作り茎へ移せるからである。従って $\mathcal E xt_P^{n-1}(i_*\mathcal O_X,\omega_P)=0$ である。

$i_*\mathcal O_X(-m)$ に対する局所 Ext は $i_*\mathcal O_X$ に対する局所 Ext の $m$ ねじりである。これは局所自由分解をねじって双対化すると確認できる。有限個の連接な局所 Ext 層について、十分大きい $m$ で正次数のコホモロジーは Serre の消滅により零になる。従って局所から大域への系列では零次の列だけが残り、$\operatorname{Ext}_P^{n-1}(i_*\mathcal O_X(-m),\omega_P)=0$ である。射影空間上の Serre 双対性から $H^1(X,\mathcal O_X(-m))=0$ となる。

選んだ埋め込みの $\mathcal O_X(1)$ を $\mathcal O_X(aD)$ と書けば、十分大きい $m$ で
$$
0\to\mathcal O_X(-maD)\to\mathcal O_X\to\mathcal O_{maD}\to0
$$
の大域切断は右端へ全射となる。$H^0(X,\mathcal O_X)=k$ であり、$maD$ は空でないので $H^0(maD,\mathcal O_{maD})=k$ である。もし台が非連結なら、二つの開閉部分でそれぞれ $1,0$ と置く切断から、$k$ にない冪等元が得られて矛盾する。従って $D$ は連結である。滑らかな $X$ の一般の超平面切断は平滑性の Bertini 定理により滑らかであり、上の結果で連結である。滑らかなスキームでは局所環が整域なので異なる既約成分は出会えない。従って連結な切断は既約である。
:::

## 例：曲線と曲面における双対性

:::example[平面曲線]
滑らかな平面次数 $e$ 曲線では $H^1(\mathcal O_C)^\vee=H^0(\mathcal O_C(e-3))$ である。$e\ge3$ なら $\mathbb P^2$ の次数 $e-3$ の斉次式が切断を与える。完全列で核の $\mathcal O(-3)$ の $H^0,H^1$ が零だからである。従って次元は $\binom{e-1}{2}$ となる。$e=4$ なら三次元である。
:::

:::example[四次曲面]
滑らかな四次曲面 $X\subset\mathbb P^3$ では $\omega_X=\mathcal O_X$ である。$0\to\mathcal O(-4)\to\mathcal O\to\mathcal O_X\to0$ と :ref[thm-ag-projective-cohom] から $H^0=k,H^1=0$ である。双対性により $H^2=k$ となり $\chi(\mathcal O_X)=2$ である。
:::

### 正則微分の次元から高次群を求める

:::example[平面五次の $\mathcal O_C(1)$]
滑らかな平面五次曲線の $\omega_C$ は $\mathcal O_C(2)$ である。従って $H^1(\mathcal O_C(1))^\vee\simeq H^0(\mathcal O_C(1))\simeq k^3$ である。Euler 標数は超曲面の列から $\chi(\mathcal O_C(1))=\chi(\mathcal O_{\mathbb P^2}(1))-\chi(\mathcal O_{\mathbb P^2}(-4))=3-3=0$ となるので、この二つの次元と一致する。
:::

:::example[五次曲面の幾何種数]
$\mathbb P^3$ 内の滑らかな五次曲面の $\omega_X=\mathcal O_X(1)$ である。$0\to\mathcal O(-4)\to\mathcal O(1)\to\mathcal O_X(1)\to0$ で大域切断を取ると、核の零次と一次が零なので四つの一次座標が基底になる。従って $h^2(\mathcal O_X)=4$ である。また構造層の制限列から $h^1=0,h^0=1$ なので $\chi=5$ となる。
:::

### 有限二重被覆の微分と跡を対応させる

平面曲線の標準層は、周囲の微分から方程式を一つ取り除いて計算できる。射影直線への有限な二重被覆として表した場合には、基底への Hom として計算することもできる。二つの局所表示を比較すると、分岐点を含めてどの微分が正則か分かる。

:::example[奇数次数の二重被覆で標準層を求める]
標数が二でない体で $F(t)$ を重根を持たないモニックな次数 $2g+1$ の多項式とし、$g\ge1$ とする。有限開の被覆を $B_0=k[t,u]/(u^2-F(t))$ とする。無限遠側では $s=1/t$、$v=s^{g+1}u$ と置くと
$$
v^2=sG(s),\qquad G(s)=s^{2g+1}F(1/s),\quad G(0)=1
$$
となる。そこで $B_1=k[s,v]/(v^2-sG(s))$ を使い、重なりを座標変換で貼る。両代数はそれぞれ基底 $1,u$ と $1,v$ を持つ階数二の自由加群なので、有限平坦な射 $f:C\to\mathbb P^1$ を得る。有限開の平滑性は $F$ が重根を持たないことから従う。無限遠の $(s,v)=(0,0)$ では方程式の $s$ 微分が $-G(0)=-1$ であり、ここでも平滑である。無限遠側の他の点は有限開との重なりに入る。

有限開の $B_0$ 上では $\Hom_{k[t]}(B_0,k[t])$ が $u$ の係数を取る形式 $\lambda_0$ を基底とする。無限遠側では対応する基底を $\lambda_1$ とする。重なりで $v=s^{g+1}u$ なので、$\lambda_1(v)=1$、$\lambda_1(1)=0$ から $\lambda_1=s^{-g-1}\lambda_0$ である。基底の標準層の局所基底 $dt,ds$ は $ds=-s^2dt$ を満たす。従って有限射の双対化層の局所基底 $\lambda_0\,dt$ と $\lambda_1\,ds$ は
$$
\lambda_1\,ds=-s^{1-g}\lambda_0\,dt=-t^{g-1}\lambda_0\,dt
$$
の関係をもつ。片方の符号を基底へ吸収すれば $\omega_C\simeq f^*\mathcal O_{\mathbb P^1}(g-1)$ を得る。

微分からも同じ式が出る。$u\ne0$ で $dt/(2u)$ は自由基底であり、$u=0$ で $F'(t)$ が単元なので $du/F'(t)$ として延びる。無限遠側では $ds/(2v)$ が同様に正則で、重なりでは
$$
\frac{ds}{2v}=-s^{1-g}\frac{dt}{2u}
$$
である。従って $\lambda_0\,dt$ と $dt/(2u)$、$\lambda_1\,ds$ と $ds/(2v)$ を対応させると、有限射の Hom の記述と正則微分が一致する。

例えば $g=2$、標数零で $F(t)=t^5-1$ とする。構造層の順像は $1,u$ の二方向で、重なりの第二基底が $s^{g+1}u=t^{-g-1}u$ なので $\mathcal O\oplus\mathcal O(-3)$ である。射影公式と $\omega_C=f^*\mathcal O(1)$ から $f_*\omega_C=\mathcal O(1)\oplus\mathcal O(-2)$ になる。従って $H^0(C,\omega_C)$ は二次元で、微分として $dt/(2u),t\,dt/(2u)$ が基底である。同時に $H^1(C,\mathcal O_C)$ の二次元は、基底上の $\mathcal O(-3)$ の二つの Laurent 類から出る。このように有限射を使うと、一方の群が二次元であることと、双対の微分の具体式を同じ貼り合わせから読める。
:::

### 双対性から Hilbert 多項式の対称性を読む

正と負のねじりでは大域切断の次元が大きく異なる。その差は高次群が担うため、Euler 標数では二つのねじりが多項式の対称性としてつながる。

:::proposition[標準層がねじりで書ける場合の対称性]
$X$ を純次元 $d$ の Cohen–Macaulay 射影スキームとし、$\omega_X\simeq\mathcal O_X(a)$ とする。このとき Hilbert 多項式は
$$
P_X(m)=(-1)^d P_X(a-m)
$$
を満たす。$d\ge1$ で先頭二係数を $b_d,b_{d-1}$ とすると、$b_{d-1}=-da b_d/2$ である。
:::

:::proof
可逆層 $\mathcal O_X(m)$ に Serre 双対性を適用すると、$H^i(\mathcal O_X(m))$ の双対は $H^{d-i}(\mathcal O_X(a-m))$ である。全ての群は有限次元で、次数 $i$ の次元に $(-1)^i$ を掛けて足すと、添字を $d-i$ に取り替える際に全体の符号 $(-1)^d$ が出る。この等式は全ての整数 $m$ で成り立つので多項式として成り立つ。右辺の $m^{d-1}$ の係数は $(-1)^d[(-1)^{d-1}da b_d+(-1)^{d-1}b_{d-1}]=-da b_d-b_{d-1}$ である。左辺の係数 $b_{d-1}$ と比較して公式を得る。
:::

:::example[五次曲面で一次の係数を確認する]
$\mathbb P^3$ の五次曲面 $X$ は $\omega_X=\mathcal O_X(1)$ だから、$P(m)=P(1-m)$ である。次数五より二次の先頭係数は $5/2$ であり、式で $a=1,d=2$ と置くと一次の係数は $-5/2$ となる。$\chi(\mathcal O_X)=5$ なので
$$
P(m)=\frac52 m^2-\frac52m+5
$$
である。実際、$\binom{m+3}{3}-\binom{m-2}{3}$ を展開すると同じ値を得る。例えば $m=2,-1$ は互いに移るため両方で $P=10$ である。$m=-1$ では大域切断が零であるが、双対性で $h^2(\mathcal O_X(-1))=h^0(\mathcal O_X(2))=10$ となり、対称性に寄与する。
:::

:::example[積の曲面で中間次数の双対を測る]
$X=\mathbb P^1\times\mathbb P^1$、$\mathcal L=\mathcal O(2,-3)$ とする。第二因子への射 $f$ の順像は $f_*\mathcal L=\mathcal O(-3)^{\oplus3}$ であり、第一因子の $\mathcal O(2)$ に高次群がないため $R^q f_*\mathcal L=0$（$q>0$）である。従って $h^1(\mathcal L)=3\cdot2=6$ で、零次と二次は零である。標準層は $\mathcal O(-2,-2)$ なので $\mathcal L^\vee\otimes\omega_X=\mathcal O(-4,1)$ である。今度は第一因子への射を使い、第二因子の $\mathcal O(1)$ の切断二つを先に取ると、順像は $\mathcal O(-4)^{\oplus2}$ で高次順像は零である。従って双対の一次は $2\cdot3=6$ になる。一つのねじりが正で一つが負になると、最高次ではなく中間次数にこの対が現れる。
:::

## 演習問題

::::exercise{level=basic}
$\mathbb P^2$ で $X_0^2X_1$ に双対な単項式を求めよ。

:::solution
指数 $(2,1,0)$ と足して $(-1,-1,-1)$ にするため $X_0^{-3}X_1^{-2}X_2^{-1}$ である。
:::
::::

::::exercise{level=basic}
局所自由な $\mathcal L$ の $\mathcal E xt^1_X(\mathcal L,\mathcal G)$ を求めよ。

:::solution
局所で自由加群からの Hom は有限個のコピーを取る完全関手なので一次導来関手は零である。
:::
::::

::::exercise{level=basic}
滑らかな平面五次曲線の $h^1(\mathcal O_C)$ を求めよ。

:::solution
標準層は $\mathcal O_C(2)$ で、二次式六個が全切断の基底である。双対性より $h^1=6$ となる。
:::
::::

::::exercise{level=basic}
$\mathbb P^2$ の $h^2(\mathcal O(-3))$ を双対性で求めよ。

:::solution
双対は $H^0(\mathcal O(3-3))=k$ なので次元は $1$ である。
:::
::::

::::exercise{level=advanced}
$B=k[t]/(t^2)$ に対し $\operatorname{Ext}_{k[t]}^1(B,k[t])$ を求めよ。

:::solution
自由分解 $0\to k[t]\xrightarrow{t^2}k[t]\to B\to0$ を双対化すると同じ掛け算が現れ、その余核は $k[t]/(t^2)$ である。
:::
::::

::::exercise{level=advanced}
滑らかな平面三次曲線で $H^1(\mathcal O_C(1))$ を求めよ。

:::solution
標準層は自明なので双対は $H^0(\mathcal O_C(-1))$ である。超曲面の完全列で $\mathcal O(-1)$ の $H^0$ と $\mathcal O(-4)$ の $H^1$ が零なのでこれは零である。
:::
::::

::::exercise{level=advanced}
滑らかな五次曲面 $X\subset\mathbb P^3$ の $h^2(\mathcal O_X)$ を求めよ。

:::solution
標準層は $\mathcal O_X(1)$ である。四つの一次式は $\mathcal O(-4)$ の $H^0,H^1$ が零なので全切断を与える。双対性で次元は $4$ である。
:::
::::

::::exercise{level=advanced}
滑らかな平面三次曲線で $H^0(\mathcal O_C(-2))$ と $H^1(\mathcal O_C(-2))$ を求めよ。

:::solution
超曲面の完全列で $\mathcal O(-2)$ の $H^0$ と $\mathcal O(-5)$ の $H^1$ が零なので零次は零である。一次の双対は $\mathcal O_C(2)$ の切断であり、$2<3$ のため二次式六個が基底となる。
:::
::::

::::exercise{level=advanced}
$B=k[x,y]/(x^3,y^2)$ とし、$\lambda$ を $x^2y$ の係数を取る線形形式とする。ソクルを決定し、$1,x,x^2,y,xy,x^2y$ に対する双対基底を $\lambda(ab)$ を用いて求めよ。$\mu$ が $xy$ の係数を取るとき、$\mu=b\lambda$ を満たす $b$ も求め、この $\mu$ 自身が双対化加群を生成するか答えよ。

:::solution
六つの単項式が環の基底となる。$x$ 倍で殺される元は $x^2,x^2y$ の張る空間であり、$y$ も掛けて零になるのは $x^2y$ の倍数だけである。従ってソクルは $kx^2y$ であり、一次元である。二単項式の積の $\lambda$ が $1$ になるのは、$x$ の指数の和が $2$、$y$ の指数の和が $1$ に等しいときだけである。従って対応する双対基底は
$$
x^2y,\quad xy,\quad y,\quad x^2,\quad x,\quad 1
$$
である。他の組はその指数に達さないか、指数が上限を超えて零になる。よってこの対の行列は、双対基底を右に置けば単位行列になる。特に $b\mapsto b\lambda$ は $B$ 加群の同型である。

$\mu$ が非零に評価する基底は $xy$ だけなので、対応する双対基底 $x$ を取れば $\mu=x\lambda$ である。実際 $\lambda(x\cdot xy)=1$ で、他の基底には零である。同型の単射性からこれが一意な $b$ である。$x$ は極大イデアルに属し単元でないため、$B\mu=(x)\lambda$ は真の部分加群である。実際、$x^2\ne0$ が $x^2\mu=x^3\lambda=0$ として $\mu$ を殺す。従って $\mu$ は双対化加群を生成せず、$\mu$ が作る積の対は退化する。
:::
::::

::::exercise{level=advanced}
滑らかな六次平面曲線 $C$ に対して $\omega_C$ を求め、$H^0(\mathcal O_C(m))$ と $H^1(\mathcal O_C(m))$ の次元を $m=-1,0,1,2,3,4$ で計算せよ。また Euler 標数の対称性を確かめよ。

:::solution
$\mathbb P^2$ の標準層は $\mathcal O(-3)$ であり、六次の方程式は $\mathcal O(6)$ の切断である。従って随伴公式から $\omega_C=\mathcal O_C(3)$ となる。構造層の列
$$
0\to\mathcal O_{\mathbb P^2}(m-6)\to
\mathcal O_{\mathbb P^2}(m)\to\mathcal O_C(m)\to0
$$
で、周囲の一次群は常に零である。従って $m<0$ の零次は零、$0\le m<6$ の零次は全次数 $m$ の多項式であり、次元は $\binom{m+2}{2}$ である。Serre 双対性は $h^1(\mathcal O_C(m))=h^0(\mathcal O_C(3-m))$ を与える。従って表は
$$
\begin{array}{c|rrrrrr}
m&-1&0&1&2&3&4\\ \hline
h^0&0&1&3&6&10&15\\
h^1&15&10&6&3&1&0\\
\chi&-15&-9&-3&3&9&15
\end{array}
$$
となる。例えば $m=-1$ の $h^1$ は、周囲の $\mathcal O(-7)$ の最高次から出る $\binom{6}{2}=15$ でも確かめられる。

長完全列から $\chi(\mathcal O_C(m))=\binom{m+2}{2}-\binom{m-4}{2}=6m-9$ である。二項係数はここでは負の引数を含めて多項式として評価する。この多項式は $P_C(3-m)=9-6m=-P_C(m)$ を満たし、表でも $m=-1$ と $4$、$0$ と $3$、$1$ と $2$ の値が互いに負になる。次元一なので、双対性で零次と一次を交換する際に Euler 標数の符号が変わる。
:::
::::

::::exercise{level=advanced}
標数が二でない体で $C\to\mathbb P^1$ を $u^2=t^7-1$ の二重被覆の滑らかなモデルとする。ただし標数は七でもないとする。正則微分の基底と $H^1(\mathcal O_C)$ の次元を求めよ。無限遠の局所座標を使って、さらに一つ次数を上げた $t^3dt/(2u)$ が正則でないことを示せ。

:::solution
有限開では $t^7-1$ が重根を持たない。実際導関数 $7t^6$ が零になる点 $t=0$ は根でない。無限遠側で $s=1/t$、$v=s^4u$ とすると $v^2=s-s^8$ である。$s=0$ の点は $(0,0)$ だけであり、$s-s^8$ の導関数はそこで $1$ なので、$v$ が局所パラメータになる。有限射の構造層の順像は $1$ と $u$ の二方向から $\mathcal O\oplus\mathcal O(-4)$ であり、$H^1(\mathcal O_C)$ の次元は射影直線上で $0+3=3$ である。

微分 $\eta=dt/(2u)$ は有限側の分岐点で $du/(7t^6)$ と書け、有限開すべてで基底である。無限遠のパラメータの近くでは $ds/(2v)=dv/(1-8s^7)$ が基底であり、
$$
t^j\eta=-s^{2-j}\,\frac{ds}{2v}
$$
となる。従って $j=0,1,2$ なら正則であり、三つは有限開の係数 $1,t,t^2$ が独立なので一次独立である。双対性から $h^0(\omega_C)=h^1(\mathcal O_C)=3$ であるため、これらが基底である。$j=3$ の場合は単元倍の $dv$ に $s^{-1}$ を掛けるが、$v^2=s(1-s^7)$ かつ $1-s^7$ が単元なので $s$ の位数は二である。従って $j=3$ の微分は無限遠で二位の極をもち、正則でない。
:::
::::

## 確認テスト

::::quiz{answer="3" tol="1e-6"}
滑らかな平面四次曲線の $h^1(\mathcal O_C)$ はいくつか。

:::explanation
標準層の全切断は独立な三つの一次式であり、双対性より次元は $3$ である。
:::
::::

::::quiz
滑らかな $d$ 次元射影多様体で $H^i(\mathcal L)$ の双対を選べ。

- [x] $H^{d-i}(\mathcal L^\vee\otimes\omega_X)$
- [ ] $H^i(\mathcal L^\vee)$
- [ ] $H^d(\mathcal L)$

:::explanation
局所自由層について Serre 双対の Ext は双対層と標準層のテンソル積のコホモロジーになる。従って次数は $d-i$ である。
:::
::::

::::quiz{answer="2" tol="1e-6"}
滑らかな四次曲面の $\chi(\mathcal O_X)$ を求めよ。

:::explanation
零次と二次の次元は $1$、一次は零であるため $1-0+1=2$ となる。
:::
::::

::::quiz
有限階数の局所自由層で常に消えるものを選べ。

- [x] 正次数の局所 Ext
- [ ] 正次数の大域 Ext
- [ ] 層の任意の大域切断

:::explanation
局所 Hom は局所自由の場合に完全だが、大域 Hom にはさらに大域切断が含まれる。従って正次数の局所 Ext が常に零である。
:::
::::

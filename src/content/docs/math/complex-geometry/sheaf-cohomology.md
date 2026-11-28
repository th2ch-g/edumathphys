---
title: 複素多様体上の層係数コホモロジー
description: 連接層と局所分解を導入し、直線束の整数コホモロジー類と双対性、指数公式を説明する
notation:
  - symbol: '\mathcal O(E)'
    meaning: 正則ベクトル束 E の正則切断の層
  - symbol: 'c_1(L)'
    meaning: 直線束 L の第一 Chern 類
  - symbol: '\chi(X,E)'
    meaning: E の正則切断の層コホモロジーの Euler 標数
  - symbol: '\operatorname{Td}(T_X)'
    meaning: 正則接束の Todd 類
---

## この章の目標

- 局所有限生成と関係式の有限生成から連接層を理解する。
- ベクトル束の層係数コホモロジーを Dolbeault 複体に移す。
- 指数完全列から第一 Chern 類を構成する。
- Serre 双対性と Hirzebruch–Riemann–Roch の定理を具体例で使う。

## 連接層の概観

:::definition[連接層]{#def-cg-coherent-sheaf}
$\mathcal O_X$ 加群 $\mathcal F$ が :term[連接層]{en="coherent sheaf" yomi="れんせつそう"} であるとは、局所有限生成で、任意の開部分上の有限階数自由層からの射 $\mathcal O_X^{\oplus r}\to\mathcal F$ の核も局所有限生成であることをいう。
:::

局所有限生成性は各点の近くで有限個の切断が全ての茎を生成することを意味する。したがって連接層は局所的に $\mathcal O_X^{\oplus s}\to\mathcal O_X^{\oplus r}\to\mathcal F\to0$ と表示できる。これを使うと切断の局所条件を正則行列で記録できる。各茎の環が Noether であるというだけでは、関係式が近傍全体で有限生成になることは直ちには従わない。

次に、多変数の割り算を用いて構造層の連接性を証明する。これにより、関係式を一点の茎で選ぶだけでなく、近傍の全ての点で同じ有限個の関係式を使えることが分かる。連接層の間の射の核、像、余核も連接であることは、この節の最後で証明する。

:::example[局所自由層と商]
正則束の切断層は局所的に $\mathcal O_X^r$ なので連接である。滑らかな超曲面 $Y=\{w=0\}$ の関数層を $X$ 上に延ばした層は $\mathcal O_X/(w)$ である。$w$ がその層を消すため $X$ 上の正の階数をもつ自由層にはならないが、正則行列 $(w)$ の余核として連接である。
:::

:::example[厚みを持つ商]
$\mathcal O_{\C,0}/(z^3)$ は $1,z,z^2$ を複素基底にもち、$z^3=0$ である。台は一点だが定数、一次、二次の情報を記録しており、$\mathcal O_{\C,0}/(z)$ とは異なる。
:::

### 関係式を分解で記録する

$X=\C^2$ 上の原点のイデアル層 $\mathcal I_0=(z,w)$ に対し、行列 $(z\ w)$ で $\mathcal O_X^2\to\mathcal I_0$ と写す。関係 $za+wb=0$ の芽について、原点では $w=0$ に制限して $za(z,0)=0$ を得る。従って $a(z,0)=0$ で、Taylor 展開から $a=-wc$ と書ける。元の式に戻せば $b=zc$ である。原点の外では $z,w$ の少なくとも一つが単元なので同じ表示が得られる。従って
$$
0\longrightarrow\mathcal O_X
 \xrightarrow{c\mapsto(-wc,zc)}\mathcal O_X^2
 \xrightarrow{(a,b)\mapsto za+wb}\mathcal O_X
 \longrightarrow\mathcal O_X/\mathcal I_0
 \longrightarrow0
$$
は完全である。このように連接層自体が局所自由でなくても、局所自由層の複体を通じて測れる。例えば $\mathcal O_X/\mathcal I_0$ の台は原点であり、原点以外の茎は零になる。原点の茎は評価 $f\mapsto f(0,0)$ で $\C$ と同型である。

### Oka の連接性定理を割り算から導く

:::theorem[Oka の連接性定理]{#thm-cg-oka-coherence}
複素多様体の構造層 $\mathcal O_X$ は連接である。より具体的には、正則行列による自由層の射の核は、各点の近くで有限個の切断の芽によって生成される。
:::

:::proof
座標近傍で、次元 $n$ に関する帰納法を用いる。$n=0$ では線形代数である。次元 $n$ では先に一行の行列の核を扱う。任意個の行の場合は、第一行の核を生成する列を並べて写像 $G$ とし、残りの行との合成の核を行数に関する帰納法で求めればよい。その生成元を $G$ で送ったものは、元の全ての行を消す。逆に全ての行を消す芽は第一行の核に入るため $G$ で持ち上がり、その持ち上げは合成の核に入る。

一行 $(f_1,\ldots,f_s)$ のうち非零の芽を $f_s$ とする。すべて零なら近傍を縮めて零になり、自由層全体が核である。$f_s$ が単元なら他の成分を任意に置いて最後の成分を解ける。残りの場合、線形座標の変更と :ref[thm-cg-weierstrass-preparation] により、列の単元倍で $f_s$ を $w$ のモニック Weierstrass 多項式 $P$ にできる。$d$ を次数とし、:ref[thm-cg-weierstrass-division] で各 $f_i$ を $P$ で割る。列の加減で各 $f_i$、$i<s$ を剰余 $r_i$ に替えてよい。すなわち核を計算する行は $(r_1,\ldots,r_{s-1},P)$ となる。

$P$ を法とした各 $r_i$ の乗算を $1,w,\ldots,w^{d-1}$ の基底で表す。係数は $z'$ の正則関数であり、これらを横に並べると $d$ 行、$(s-1)d$ 列の正則行列となる。次元 $n-1$ の帰納法でその核を生成する有限個の正則列が基底近傍上にある。各列を $w$ の次数 $d$ 未満の多項式列 $(a_1,\ldots,a_{s-1})$ と読むと、$\sum_i r_i a_i$ は多項式として $P$ で割り切れる。最後の成分 $a_s=-\sum_i r_i a_i/P$ を加えれば、求める行の核に入る。この有限族に、$i<s$ ごとに第 $i$ 成分 $P$、最後の成分 $-r_i$、残りは零の列を加える。

これらの列が近傍の各点で生成することを確かめる。$P\ne0$ の点では後者の列だけで生成できる。$P=0$ の点 $(b,w_0)$ では $z'=b$ の全根を相異なる根ごとに小円板で分ける。基底を縮めれば、各円周は根を避け同じ個数の根を含む。核の芽 $(u_i)$ を $w_0$ の根を囲む円板との積上で取り、他の円板では零に置く。これは互いに離れた積領域上で正則であり、そこでは $\sum_{i<s}r_i u_i+Pu_s=0$ である。各 $u_i$ を、全円周の上での :ref[thm-cg-weierstrass-division] の剰余積分を足して多項式 $R_i$ に移す。この積分の差商の式により $R_i$ は次数 $d$ 未満で、全ての根において $u_i$ と重複度まで一致する。実際、固定した $z'$ について各円板内で $u_i=q_iP+R_i$ となる。これは円周の和での Cauchy 公式を各円板で使えば従う。

従って $\sum_{i<s}r_iR_i$ の $P$ による剰余は零である。$R_i$ の係数列は基底側の行列の核に属するので、そこで選んだ生成元を使って表せる。対応する生成元を $(u_i)$ から引くと、最初の $s-1$ 成分は $P$ の倍数になる。それらは第 $i$ 成分 $P$ の列を引いて消せる。残る最後の成分は $P$ を掛けると零であり、局所環は整域なので零である。従ってどの点の核の芽も有限族から得られ、帰納が完了する。
:::

有限階数自由層が連接であることは上の核の記述と定義から従う。有限表示を持つ層の間で関係と生成元を持ち上げる議論と合わせると、前に述べた連接層の核、像、余核の連接性を使える。

### 関係と台を別々に計算する

有限表示 $\mathcal O_X^m\xrightarrow A\mathcal O_X^r\to\mathcal F\to0$ があるとき、$\mathcal F_p/\mathfrak m_p\mathcal F_p$ は数値行列 $A(p)$ の余核である。実際、茎の写像の余核に剰余体とのテンソル積を取り、テンソル積の右完全性を用いればよい。Nakayama の補題から、この余核が零なら $\mathcal F_p=0$ である。それ以外の点では $\mathcal F_p\ne0$ なので、台は $r$ 次の小行列式がすべて零となる集合である。台だけで局所の層を復元できるわけではない。同じ台の上でも $A$ の成分が高い位数で消えることがあるからである。

:::example[行列で与える台の重複]
$A=\begin{pmatrix}z&1\\0&z\end{pmatrix}$ を円板上で考える。行列式は $z^2$ なので台は原点である。余核で第二列が与える関係は $e_1+ze_2=0$ であり、第一列の関係 $ze_1=0$ は $z^2e_2=0$ となる。従って余核は $\mathcal O/(z^2)$ である。$A(0)$ の余核の次元は $1$ だが、茎を $\C$ ベクトル空間と見たときの次元は $2$ である。
:::

:::example[イデアルの生成数が跳ぶ]
原点のイデアル $\mathcal I=(z,w)$ は原点の外で $\mathcal O$ に等しい。原点では $\mathcal I_0/\mathfrak m_0\mathcal I_0=\mathfrak m_0/\mathfrak m_0^2$ が二次元なので、単一の芽で生成できない。従ってこのイデアル層は連接だが局所自由ではない。局所自由階数は、近傍での最小生成数の一定性を要求する。
:::

### 連接層の局所表示で係数を計算する

:::example[割り算を有限自由加群で再現する]
$P=w^{3} - z^{2}$ を法として $g=w^{6} + w^{4} z + z^{2}$ の余りを求める。$P$ は $w$ のモニック式であるため、係数の変数 $z$ を固定せずに割り算ができる。

各段階で最大の $w$ 次数を消す。$R_0=g$ とすると次の計算になる。

$$
\begin{aligned}
R_{1}&=R_{0}-(w^{3})P=w^{4} z + w^{3} z^{2} + z^{2} \\
R_{2}&=R_{1}-(w z)P=w^{3} z^{2} + w z^{3} + z^{2} \\
R_{3}&=R_{2}-(z^{2})P=w z^{3} + z^{4} + z^{2} \\
q&=w^{3} + w z + z^{2},\qquad r=w z^{3} + z^{4} + z^{2}
\end{aligned}
$$

ここで $g=qP+r$ であり、$r$ は $P$ より低次数である。このような余りが二つあれば、その差が $P$ の倍数かつ低次数となるため零である。

同じ計算を $\C\{z\}[w]/(P)$ の基底 $1,w,\ldots,w^{2}$ で行う。$w$ 倍の各列は、各基底に $w$ を掛けた余りの係数列である。これを $C$ とし、$g(w)$ の $w$ を $C$ に置き換えて最初の基底ベクトルに作用させる。

$$
\begin{aligned}
C&=\left[\begin{matrix}0 & 0 & z^{2}\\1 & 0 & 0\\0 & 1 & 0\end{matrix}\right] \\
g(C)e_1&=\left[\begin{matrix}z^{4} + z^{2}\\z^{3}\\0\end{matrix}\right]
\end{aligned}
$$

各成分を基底に戻すと先ほどの余りと一致する。$z$ を具体値にしても恒等式が保たれる。これに対して $w$ の先頭係数が消える場合は、同じ計算の途中で逆元を使えない。
:::

:::example[割り算を有限自由加群で再現する]
$P=w^{2} + w z^{2} + z$ を法として $g=w^{5} + 2 w z + 3$ の余りを求める。$P$ は $w$ のモニック式であるため、係数の変数 $z$ を固定せずに割り算ができる。

各段階で最大の $w$ 次数を消す。$R_0=g$ とすると次の計算になる。

$$
\begin{aligned}
R_{1}&=R_{0}-(w^{3})P=- w^{4} z^{2} - w^{3} z + 2 w z + 3 \\
R_{2}&=R_{1}-(- w^{2} z^{2})P=w^{3} z^{4} - w^{3} z + w^{2} z^{3} + 2 w z + 3 \\
R_{3}&=R_{2}-(w \left(z^{4} - z\right))P=- w^{2} z^{6} + 2 w^{2} z^{3} - w z^{5} + w z^{2} + 2 w z + 3 \\
R_{4}&=R_{3}-(- z^{6} + 2 z^{3})P=w z^{8} - 3 w z^{5} + w z^{2} + 2 w z + z^{7} - 2 z^{4} + 3 \\
q&=w^{3} - w^{2} z^{2} + w \left(z^{4} - z\right) - z^{6} + 2 z^{3},\qquad r=w \left(z^{8} - 3 z^{5} + z^{2} + 2 z\right) + z^{7} - 2 z^{4} + 3
\end{aligned}
$$

ここで $g=qP+r$ であり、$r$ は $P$ より低次数である。このような余りが二つあれば、その差が $P$ の倍数かつ低次数となるため零である。

同じ計算を $\C\{z\}[w]/(P)$ の基底 $1,w,\ldots,w^{1}$ で行う。$w$ 倍の各列は、各基底に $w$ を掛けた余りの係数列である。これを $C$ とし、$g(w)$ の $w$ を $C$ に置き換えて最初の基底ベクトルに作用させる。

$$
\begin{aligned}
C&=\left[\begin{matrix}0 & - z\\1 & - z^{2}\end{matrix}\right] \\
g(C)e_1&=\left[\begin{matrix}z^{7} - 2 z^{4} + 3\\z^{8} - 3 z^{5} + z^{2} + 2 z\end{matrix}\right]
\end{aligned}
$$

各成分を基底に戻すと先ほどの余りと一致する。$z$ を具体値にしても恒等式が保たれる。これに対して $w$ の先頭係数が消える場合は、同じ計算の途中で逆元を使えない。
:::

### 商層の階数を乗算行列で追う

:::example[割り算を有限自由加群で再現する]
$P=w^{2} + w z - z^{2}$ を法として $g=w^{6} + w + z$ の余りを求める。$P$ は $w$ のモニック式であるため、係数の変数 $z$ を固定せずに割り算ができる。

各段階で最大の $w$ 次数を消す。$R_0=g$ とすると次の計算になる。

$$
\begin{aligned}
R_{1}&=R_{0}-(w^{4})P=- w^{5} z + w^{4} z^{2} + w + z \\
R_{2}&=R_{1}-(- w^{3} z)P=2 w^{4} z^{2} - w^{3} z^{3} + w + z \\
R_{3}&=R_{2}-(2 w^{2} z^{2})P=- 3 w^{3} z^{3} + 2 w^{2} z^{4} + w + z \\
R_{4}&=R_{3}-(- 3 w z^{3})P=5 w^{2} z^{4} - 3 w z^{5} + w + z \\
R_{5}&=R_{4}-(5 z^{4})P=- 8 w z^{5} + w + 5 z^{6} + z \\
q&=w^{4} - w^{3} z + 2 w^{2} z^{2} - 3 w z^{3} + 5 z^{4},\qquad r=w \left(1 - 8 z^{5}\right) + 5 z^{6} + z
\end{aligned}
$$

ここで $g=qP+r$ であり、$r$ は $P$ より低次数である。このような余りが二つあれば、その差が $P$ の倍数かつ低次数となるため零である。

同じ計算を $\C\{z\}[w]/(P)$ の基底 $1,w,\ldots,w^{1}$ で行う。$w$ 倍の各列は、各基底に $w$ を掛けた余りの係数列である。これを $C$ とし、$g(w)$ の $w$ を $C$ に置き換えて最初の基底ベクトルに作用させる。

$$
\begin{aligned}
C&=\left[\begin{matrix}0 & z^{2}\\1 & - z\end{matrix}\right] \\
g(C)e_1&=\left[\begin{matrix}5 z^{6} + z\\1 - 8 z^{5}\end{matrix}\right]
\end{aligned}
$$

各成分を基底に戻すと先ほどの余りと一致する。$z$ を具体値にしても恒等式が保たれる。これに対して $w$ の先頭係数が消える場合は、同じ計算の途中で逆元を使えない。
:::

:::example[割り算を有限自由加群で再現する]
$P=w^{3} + w z^{2} + z^{3}$ を法として $g=w^{7} + w^{2} z$ の余りを求める。$P$ は $w$ のモニック式であるため、係数の変数 $z$ を固定せずに割り算ができる。

各段階で最大の $w$ 次数を消す。$R_0=g$ とすると次の計算になる。

$$
\begin{aligned}
R_{1}&=R_{0}-(w^{4})P=- w^{5} z^{2} - w^{4} z^{3} + w^{2} z \\
R_{2}&=R_{1}-(- w^{2} z^{2})P=- w^{4} z^{3} + w^{3} z^{4} + w^{2} z^{5} + w^{2} z \\
R_{3}&=R_{2}-(- w z^{3})P=w^{3} z^{4} + 2 w^{2} z^{5} + w^{2} z + w z^{6} \\
R_{4}&=R_{3}-(z^{4})P=2 w^{2} z^{5} + w^{2} z - z^{7} \\
q&=w^{4} - w^{2} z^{2} - w z^{3} + z^{4},\qquad r=w^{2} \left(2 z^{5} + z\right) - z^{7}
\end{aligned}
$$

ここで $g=qP+r$ であり、$r$ は $P$ より低次数である。このような余りが二つあれば、その差が $P$ の倍数かつ低次数となるため零である。

同じ計算を $\C\{z\}[w]/(P)$ の基底 $1,w,\ldots,w^{2}$ で行う。$w$ 倍の各列は、各基底に $w$ を掛けた余りの係数列である。これを $C$ とし、$g(w)$ の $w$ を $C$ に置き換えて最初の基底ベクトルに作用させる。

$$
\begin{aligned}
C&=\left[\begin{matrix}0 & 0 & - z^{3}\\1 & 0 & - z^{2}\\0 & 1 & 0\end{matrix}\right] \\
g(C)e_1&=\left[\begin{matrix}- z^{7}\\0\\2 z^{5} + z\end{matrix}\right]
\end{aligned}
$$

各成分を基底に戻すと先ほどの余りと一致する。$z$ を具体値にしても恒等式が保たれる。これに対して $w$ の先頭係数が消える場合は、同じ計算の途中で逆元を使えない。
:::

### 有限表示から核・像・余核の連接性を示す

連接層 $F,G$ の射 $\phi:F\to G$ を取る。点の近傍で自由層 $A\to F$ を全射とする。合成 $A\to G$ の核 $B$ は $G$ の連接性で局所有限生成である。その像は $\ker\phi$ に等しいため、$\ker\phi$ は局所有限生成である。さらに任意の有限自由層から $\ker\phi$ への射の核は、$F$ への合成の核と同じである。$F$ は連接なので核は局所有限生成であり、$\ker\phi$ は連接となる。$\operatorname{im}\phi$ は $F$ の有限個の局所生成元の像で生成される。また $G$ の部分層なので、自由層からの射の核を $G$ への合成で計算できる。従って像も連接である。

最後に $Q=G/\operatorname{im}\phi$ を取る。$Q$ は $G$ の商なので局所有限生成である。任意の有限自由層 $C\to Q$ を与え、点の近くでその基底の像を $G$ へ持ち上げる。像層の有限個の局所生成元を並べて自由層 $D\to\operatorname{im}\phi$ とし、$C\oplus D\to G$ を持ち上げと生成元によって作る。この核は $G$ の連接性から局所有限生成であり、$C$ に射影した像は $C\to Q$ の核である。実際、一方の包含は合成から従い、逆には $C$ の元の持ち上げが像層に入るので $D$ の元で打ち消せばよい。よって $Q$ も連接である。

:::example[余核の射を後から取る]
円板で $\mathcal O/(z^3)\to\mathcal O/(z^2)$ の自然な射を取る。核は $z^2$ の像の定数倍であり、$z$ を掛けると零となるので $\mathcal O/(z)$ に同型である。これに対し $\mathcal O/(z^2)\to\mathcal O/(z^3)$ を $z$ 倍で定めると、核は零、余核が $\mathcal O/(z)$ となる。同じ二層の間でも射を指定してから関係を計算する必要がある。
:::

## Čech コホモロジーと Dolbeault 同型

パラコンパクトな多様体では開被覆の細分を取った Čech コホモロジーが層のコホモロジーと一致する。[層](/math/category-homological/sheaves/) を参照する。ある被覆の全ての空でない有限交差で $\mathcal F$ の正次数コホモロジーが零なら、この被覆だけで計算できる。これが Leray の被覆定理である。証明は :ref[thm-cat-leray-cover] を参照する。そこでは交差上のコホモロジーを先に取る計算と、大域切断を先に取る計算が同じ全複体を与える。

:::theorem[束に値をとる Dolbeault 同型]{#thm-cg-global-dolbeault}
正則ベクトル束 $E$ と正則 $p$ 形式の層 $\Omega_X^p$ に対し
$$
H^q(X,\Omega_X^p\otimes\mathcal O(E))
\cong
\frac{\Ker(\bar\partial_E:A^{p,q}(X,E)\to A^{p,q+1}(X,E))}
{\Img(\bar\partial_E:A^{p,q-1}(X,E)\to A^{p,q}(X,E))}
$$
である。$A^{p,q}(X,E)$ は滑らかな束値微分形式の空間である。
:::

:::proof
正則な局所枠で成分ごとに $\bar\partial$ を取る。遷移行列 $g$ は $\bar\partial g=0$ なので結果は貼り合う。零次の核は正則 $p$ 形式の切断であり、正次数での局所完全性は :ref[thm-cg-local-dolbeault] である。滑らかな形式の層は、$1$ の分割を掛けられるため :ref[def-cg-cech] に続く縮約の計算で非巡回である。従ってこの分解の大域切断を取った複体が層コホモロジーを計算する。
:::

:::example[曲線への制限]
曲線では $A^{0,2}=0$ なので $H^1(X,\mathcal O(E))$ は大域的な $(0,1)$ 形式の余核である。複素二次元では閉性が条件になる。例えば $\bar w\dd\bar z$ の $\bar\partial$ は $\dd\bar w\wedge\dd\bar z$ なので、この形式はそもそも類を表さない。
:::

### 二重複体の符号を固定する

$C^{p,q}=\check C^p(\mathfrak U,\mathcal A^{0,q}(E))$ で $D=\delta+(-1)^p\bar\partial_E$ と置く。制限と微分は可換であるが、$(-1)^p$ が横方向の微分後に反転するので $D^2=0$ である。正則 $1$ コサイクル $c$ に対し $u=hc$ とすると $\delta u=c$ なので、全複体で $c-Du=-\bar\partial_E u$ となる。ここで $c$ を $(1,0)$ の成分、$\bar\partial_Eu$ を $(0,1)$ の成分と見ている。Dolbeault 側を $+\bar\partial_E u$ と定義する直接の比較では、対応する Čech 側の $1$ 次に符号を掛ける規約を採ることになる。どちらでも核と像の判定は同じである。

:::example[局所候補を貼ると現れる形式]
二重交差で $\mathcal O(E)$ 値の $c_{12}$ が与えられ、$1$ の分割が $\rho_1+\rho_2=1$ を満たすとする。$u_1=-\rho_2 c_{12}$、$u_2=\rho_1c_{12}$ とおけば差は $c_{12}$ である。$\bar\partial_Eu_1=-\bar\partial\rho_2\,c_{12}=\bar\partial\rho_1\,c_{12}=\bar\partial_Eu_2$ なので滑らかな形式が貼り合わさる。$\rho_1$ を別の $\rho'_1$ に替えた差は、局所では $(\rho'_1-\rho_1)c_{12}$ で与えられる大域切断の $\bar\partial_E$ となる。
:::

### Dolbeault 代表の平均

:::example[平均と非零 Fourier モードを分離する]
$X=\C/(\Z+6i\Z)$ 上で $\alpha=a\,\mathrm d\bar z$ とし、周期座標を $z=x+iy$ とする。係数と微分の規約は次である。

$$
\begin{aligned}
a(x,y)&=e^{2 i \pi x} + 2 e^{2 i \pi \left(- x + \frac{y}{3}\right)} + 6 \\
\partial_{\bar z}&=\tfrac12(\partial_x+i\partial_y)
\end{aligned}
$$

周期積分では全微分の平均が零になる。したがって定数モードは $\bar\partial$ の方程式の障害である。他方、指数関数 $e^{2\pi i(mx+ny/L)}$ の二つの微分で $2\pi im,2\pi in/L$ が掛かるので、固有値は $\pi(im-n/L)$ になる。

$$
\begin{aligned}
(1,0) &\colon\quad \lambda=i \pi,\quad a=1,\quad a/\lambda=- \frac{i}{\pi} \\
(-1,2) &\colon\quad \lambda=\pi \left(- \frac{1}{3} - i\right),\quad a=2,\quad a/\lambda=\frac{3 \left(-1 + 3 i\right)}{5 \pi} \\
a_{00}&=6
\end{aligned}
$$

定数モードを除いて各係数を固有値で割る。これで得られる関数は平均零で、微分を取り直すと次を満たす。

$$
\begin{aligned}
u&=- \frac{i e^{2 i \pi x}}{\pi} + \frac{9 \left(- \frac{1}{3} + i\right) e^{2 i \pi \left(- x + \frac{y}{3}\right)}}{5 \pi} \\
u_{\bar z}&=a-6
\end{aligned}
$$

全てのモードの周期が元のトーラスと一致するため $u$ は大域的である。平坦直線束の自明な場合に、この平均がコホモロジーの一成分を指定する。
:::

:::example[平均と非零 Fourier モードを分離する]
$X=\C/(\Z+3i\Z)$ 上で $\alpha=a\,\mathrm d\bar z$ とし、周期座標を $z=x+iy$ とする。係数と微分の規約は次である。

$$
\begin{aligned}
a(x,y)&=e^{4 i \pi x} + 2 e^{2 i \pi y} + 3 e^{2 i \pi \left(x - \frac{2 y}{3}\right)} \\
\partial_{\bar z}&=\tfrac12(\partial_x+i\partial_y)
\end{aligned}
$$

周期積分では全微分の平均が零になる。したがって定数モードは $\bar\partial$ の方程式の障害である。他方、指数関数 $e^{2\pi i(mx+ny/L)}$ の二つの微分で $2\pi im,2\pi in/L$ が掛かるので、固有値は $\pi(im-n/L)$ になる。

$$
\begin{aligned}
(2,0) &\colon\quad \lambda=2 i \pi,\quad a=1,\quad a/\lambda=- \frac{i}{2 \pi} \\
(0,3) &\colon\quad \lambda=- \pi,\quad a=2,\quad a/\lambda=- \frac{2}{\pi} \\
(1,-2) &\colon\quad \lambda=\pi \left(\frac{2}{3} + i\right),\quad a=3,\quad a/\lambda=\frac{9 \left(2 - 3 i\right)}{13 \pi} \\
a_{00}&=0
\end{aligned}
$$

定数モードを除いて各係数を固有値で割る。これで得られる関数は平均零で、微分を取り直すと次を満たす。

$$
\begin{aligned}
u&=- \frac{i e^{4 i \pi x}}{2 \pi} - \frac{2 e^{2 i \pi y}}{\pi} + \frac{27 \left(\frac{2}{3} - i\right) e^{2 i \pi \left(x - \frac{2 y}{3}\right)}}{13 \pi} \\
u_{\bar z}&=a-0
\end{aligned}
$$

全てのモードの周期が元のトーラスと一致するため $u$ は大域的である。全てのモードが非零なので、与えた形式は大域的な像に入る。
:::

## 指数完全列と第一 Chern 類

:::proposition[指数完全列]{#prop-cg-exp-sequence}
複素多様体上にはアーベル群の層の完全列
$$
0\longrightarrow \underline{\Z}\longrightarrow\mathcal O_X
 \xrightarrow{f\mapsto\exp(2\pi i f)}\mathcal O_X^\times
 \longrightarrow 1
$$
がある。連結準同型 $H^1(X,\mathcal O_X^\times)\to H^2(X,\underline{\Z})$ は :term[第一 Chern 類]{en="first Chern class" yomi="だいいちちぇるんるい"} $c_1$ である。
:::

:::proof
指数の核は整数値をとる正則関数で、必ず局所定数である。零にならない関数は各点の近くで対数を持つので、末尾は茎で全射である。直線束の遷移 $g_{ij}$ を小さな交差で $g_{ij}=\exp(2\pi i u_{ij})$ と書くと、$u_{ij}+u_{jk}+u_{ki}$ は整数値の $2$ コサイクルである。対数を変更した差は整数値なので、$2$ コサイクルの変化は整数コチェインのコバウンダリである。枠を変えたときは遷移に乗法的コバウンダリが掛かり、足した対数の Čech 微分は二回目で消える。従って同じ類を得る。これは滑らかな非消失関数を使って定義する位相的第一 Chern 類と一致する。実際、滑らかな指数完全列への包含は対数の整数差を同じものに送る。滑らかな関数の層は非巡回なので、その連結写像は滑らかな直線束の分類の同型を与える。
:::

完全列の該当部分は
$$
H^1(X,\Z)\to H^1(X,\mathcal O_X)\to\operatorname{Pic}(X)
 \xrightarrow{c_1}H^2(X,\Z)\to H^2(X,\mathcal O_X)
$$
である。従って $c_1=0$ の直線束は $H^1(\mathcal O_X)$ を整数コホモロジーからの像で割った群に対応する。$c_1$ が零でも正則束として非自明な場合がある。

:::example[射影空間とトーラス]
$\mathbb P^n$ では :ref[thm-cg-projective-line-cohomology] により $\mathcal O$ の一次と二次が零なので、$c_1$ は Picard 群から $H^2(\mathbb P^n,\Z)$ への同型となる。楕円曲線では :ref[thm-cg-abel] により次数零部分が曲線自身と同型であり、同じ $c_1=0$ を持つ多数の非同型束が得られる。
:::

### 整数類の像を正則構造で捉える

滑らかな直線束に正則構造が入るための必要十分条件は、その位相的 $c_1$ が $H^2(\mathcal O_X)$ で零になることである。必要性は指数完全列の完全性による。逆向きでは、その整数類が核に入れば同じクラスをもつ正則直線束が存在する。滑らかな指数完全列の分類により、元の滑らかな束とこの正則束の滑らかな台は同型であり、同型を通して正則構造を移せる。同じ整数類を持つ正則構造の同型類は、一つの代表を選ぶと $\Ker(c_1)$ の元をテンソル積で掛けることにより並べられる。

:::example[トーラスの平坦な遷移]
$\C/(\Z+i\Z)$ で、平行移動 $z\mapsto z+1$ に対してファイバーを $1$ 倍、$z\mapsto z+i$ に対して $-1$ 倍する商直線束を取る。遷移が定数なので平坦な滑らかな接続を持ち、局所対数を加法的に取れるので $c_1=0$ である。もし正則に自明なら、引き戻した非消失整関数 $s$ が $s(z+1)=s(z)$、$s(z+i)=-s(z)$ を満たす。$|s|$ は二周期なので有界であり、Liouville の定理から定数になる。その定数は第二式で零になるため矛盾する。従って位相的に自明でも正則には非自明である。
:::

:::example[整数類でも正則直線束から来ない場合]
$X=(\C/(\Z+i\Z))^2$ とし、$z_j=x_j+iy_j$ とする。閉形式 $\mathrm dx_1\wedge\mathrm dx_2$ は、二つの実円周への射影で積分 $1$ の面積類を引き戻したものなので整数類を表す。定数層の de Rham 分解から構造層の Dolbeault 分解への自然な写像は、形式の $(0,q)$ 成分を取る写像である。実際、$\mathrm d\alpha$ の $(0,q+1)$ 成分は $\bar\partial(\alpha^{0,q})$ で、零次では局所定数関数の包含を延長する。従ってこの整数類の $H^2(X,\mathcal O_X)$ での像は $\frac14\mathrm d\bar z_1\wedge\mathrm d\bar z_2$ で表される。もしこの定数形式が $\bar\partial\beta$ なら、周期的な係数の偏導関数の平均は零なので、その定数係数も零でなければならない。従って像は非零である。指数完全列から、この整数類を第一 Chern 類に持つ正則直線束は存在しない。
:::

### 平坦な遷移と位相的な自明性

格子 $\Lambda\subset\C$ 上の準同型 $\chi:\Lambda\to U(1)$ から、$(z,v)\sim(z+\lambda,\chi(\lambda)v)$ として直線束を作る。格子の実基底それぞれに $\chi$ の値の実対数を選ぶと、準同型 $\ell:\Lambda\to\R$ で $\chi(\lambda)=e^{2\pi i\ell(\lambda)}$ となるものが得られる。局所逆写像の差として与えられる格子値の遷移コサイクルに $\ell$ を合成すれば、指数完全列の連結類は零となる。もし正則に自明なら、零にならない正則切断の持ち上げ $f$ が $f(z+\lambda)=\chi(\lambda)f(z)$ を満たす。その絶対値は周期的で有界なので Liouville の定理で $f$ は定数となり、$\chi$ は自明となる。従って非自明なユニタリ指標は、第一 Chern 類で区別できない正則直線束を与える。

:::example[位数三の正則直線束]
$\Lambda=\Z+i\Z$ で $\chi(1)=e^{2\pi i/3}$、$\chi(i)=1$ とする。$\chi^3=1$ なので三乗束は自明である。$\chi$ と $\chi^2$ は非自明なので、元の束は Picard 群で正確に位数 $3$ を持つ。指数完全列で $c_1=0$ というのは、$H^1(\mathcal O)$ からの像に入ることを意味し、その像が自明という意味にはならない。
:::

### 第一 Chern 類の積分と束の次元

:::example[完全交差の標準束と Chern 数を計算する]
$X\subset\mathbb P^{3}$ を次数 $3$ の横断的な滑らかな完全交差とし、$H=c_1(\mathcal O_X(1))$ とする。法束は $\bigoplus_a\mathcal O_X(d_a)$ であり、次元は $2$ である。

接束の完全列から $c(T_X)c(N_{X/\mathbb P^N})=c(T_{\mathbb P^N}|_X)$ を得る。Euler 列により右辺は $(1+H)^{N+1}$ である。$c(T_X)=1+c_1H+\cdots+c_nH^n$ と書いて係数を比べる。

$$
\begin{aligned}
c(N)&=3 H + 1 \\
c_{1}+(3)&=\binom{4}{1},\qquad c_{1}=1 \\
c_{2}+(3)&=\binom{4}{2},\qquad c_{2}=3 \\
c(T_X)&=3 H^{2} + H + 1
\end{aligned}
$$

法束の行列式を取った随伴公式は $K_X=\mathcal O_X(-1)$ である。これは $c(T_X)$ の一次の係数の符号を反転した値にも等しい。さらに $X$ を一般の $2$ 個の超平面で切り、交点数を積で数えると $\int_XH^2=3$ となる。

$$
\begin{aligned}
\int_X c_n(T_X)&=3\cdot3=9
\end{aligned}
$$

最高 Chern 類は実接束の Euler 類であるため、ここで得た整数は位相的 Euler 標数にも等しい。この三次曲面では $c_1(T_X)=H$ なので、超平面曲線に制限した第一 Chern 類の評価は $3$ である。
:::

:::example[完全交差の標準束と Chern 数を計算する]
$X\subset\mathbb P^{3}$ を次数 $5$ の横断的な滑らかな完全交差とし、$H=c_1(\mathcal O_X(1))$ とする。法束は $\bigoplus_a\mathcal O_X(d_a)$ であり、次元は $2$ である。

接束の完全列から $c(T_X)c(N_{X/\mathbb P^N})=c(T_{\mathbb P^N}|_X)$ を得る。Euler 列により右辺は $(1+H)^{N+1}$ である。$c(T_X)=1+c_1H+\cdots+c_nH^n$ と書いて係数を比べる。

$$
\begin{aligned}
c(N)&=5 H + 1 \\
c_{1}+(5)&=\binom{4}{1},\qquad c_{1}=-1 \\
c_{2}+(-5)&=\binom{4}{2},\qquad c_{2}=11 \\
c(T_X)&=11 H^{2} - H + 1
\end{aligned}
$$

法束の行列式を取った随伴公式は $K_X=\mathcal O_X(1)$ である。これは $c(T_X)$ の一次の係数の符号を反転した値にも等しい。さらに $X$ を一般の $2$ 個の超平面で切り、交点数を積で数えると $\int_XH^2=5$ となる。

$$
\begin{aligned}
\int_X c_n(T_X)&=11\cdot5=55
\end{aligned}
$$

最高 Chern 類は実接束の Euler 類であるため、ここで得た整数は位相的 Euler 標数にも等しい。こちらでは $c_1(T_X)=-H$ であり、同じ評価は $-5$ になる。
:::

## Stein 多様体と Cartan の定理 A・B の主張

:::definition[Stein 多様体]{#def-cg-stein}
複素多様体 $X$ が :term[Stein 多様体]{en="Stein manifold" yomi="しゅたいんたようたい"} であるとは、大域正則関数が点を分離し、各点の局所座標となる関数を大域正則関数から選べ、任意のコンパクト集合 $K$ の正則包
$$
\widehat K=\{x:|f(x)|\le\sup_K|f|\ \text{for all }f\in\mathcal O_X(X)\}
$$
がコンパクトとなることをいう。
:::

:::theorem[Cartan の定理 A と B]{#thm-cg-cartan-ab}
Stein 多様体 $X$ 上の任意の連接解析層 $\mathcal F$ について、定理 A は $\mathcal F_x$ が大域切断の芽で生成されることをいう。定理 B は全ての $q>0$ で $H^q(X,\mathcal F)=0$ であることをいう。
:::

:::proof
主張の範囲には局所自由でない連接層も含まれる。証明では Stein 多様体の消尽と近似を使い、連接層の関係式を保ったまま局所の解を大域化する。これは Dolbeault の局所補題だけでは得られない。一般の場合の証明は :cite[wells-complex] が Stein 理論の文献として挙げる Gunning–Rossi を参照する。
:::

:::example[標準被覆の交差]
$X=\C^a\times(\C^\times)^b$ は座標関数で点を分離し、局所座標も与える。座標関数と非零座標の逆数で正則包の各座標の上界と、非零座標の正の下界が得られる。正則包は閉集合なのでコンパクトである。従って Stein であり、射影空間の標準被覆の計算に定理 B を適用できる。
:::

:::example[商層の切断を持ち上げる]
Stein 多様体上で $0\to\mathcal I\to\mathcal O_X\to\mathcal O_X/\mathcal I\to0$ が連接層の完全列なら、$H^1(\mathcal I)=0$ により末項の大域切断は大域正則関数へ持ち上がる。例えば $\C$ で $\mathcal I=(z^3)$ の場合、$2+3z-z^2$ は原点での各係数 $2,3,-1$ を指定する切断の一つの持ち上げである。点 $p$ のイデアル $\mathcal I_p$ を使うと、$\mathcal I_p\mathcal F$ は連接であり、定理 B から $H^0(X,\mathcal F)\to\mathcal F_p/\mathfrak m_p\mathcal F_p$ は全射となる。剰余体上で基底を持ち上げ、Nakayama の補題を適用すれば定理 A を得る。
:::

### Stein の条件を具体的に検査する

$\C^n$ の閉複素部分多様体 $X$ では座標関数の制限が点を分離し、接空間への制限の適切な $n$ 成分のうち独立なものが局所座標となる。コンパクト集合 $K$ の正則包は、座標関数を用いれば $|z_j|\le\max_K|z_j|$ を全ての $j$ で満たす。正則包は閉じており $X$ が $\C^n$ 内で閉じているので、この制約からコンパクトになる。従ってこの場合は Stein の条件を定義から確認できる。例えば $\C^\times$ は $(z,w)$ 空間の $zw=1$ という閉部分多様体と同型である。

:::example[正則包で穴を埋めるかどうか]
$K=\{|z|=1\}$ を考える。$\C$ での正則包は、最大値原理と座標関数の評価から閉単位円板である。一方 $\C^\times$ では $z$ と $1/z$ がともに正則なので、正則包の点は $|z|\le1$ と $|1/z|\le1$ を満たす。従って正則包は $K$ 自身である。正則包は使える正則関数の環に依存する。
:::

:::example[生成と補間の違い]
$X=\C^2$ の原点のイデアル層で $z,w$ は茎を生成する。これは定理 A の例である。一方 $\mathcal O/\mathfrak m^2$ の切断は原点での値と二つの一次係数であり、$(a,b,c)$ は多項式 $a+bz+cw$ に持ち上がる。一般の連接層で同様の持ち上げの全射性を得るには、完全列の核の $H^1$ を定理 B で消す。
:::

### Stein 上の有限ジェット補間

:::example[点ごとの Taylor 条件をまとめて解く]
多項式 $f$ の次数を $3$ 以下に制限し、以下の Taylor 係数を実現する。階乗で割った導関数を係数として使う。

$$
\begin{aligned}
a=-1,\quad j=0&\colon\quad f^{(0)}(a)/0!= 1 \\
a=-1,\quad j=1&\colon\quad f^{(1)}(a)/1!= 0 \\
a=1,\quad j=0&\colon\quad f^{(0)}(a)/0!= 5 \\
a=1,\quad j=1&\colon\quad f^{(1)}(a)/1!= 1
\end{aligned}
$$

$f=\sum_k c_kz^k$ とし、$z^k$ の $a$ 周りの $j$ 次係数 $\binom kj a^{k-j}$ を並べる。$k<j$ のときは零であり、$a=0$ でも式は多項式の係数として読める。

$$
\begin{aligned}
A&=\left[\begin{matrix}1 & -1 & 1 & -1\\0 & 1 & -2 & 3\\1 & 1 & 1 & 1\\0 & 1 & 2 & 3\end{matrix}\right] \\
Ac&=\left[\begin{matrix}1\\0\\5\\1\end{matrix}\right] \\
\det A&=16 \\
c&=\left[\begin{matrix}\frac{11}{4}\\\frac{11}{4}\\\frac{1}{4}\\- \frac{3}{4}\end{matrix}\right] \\
f(z)&=- \frac{3 z^{3}}{4} + \frac{z^{2}}{4} + \frac{11 z}{4} + \frac{11}{4}
\end{aligned}
$$

この多項式の Taylor 展開は指定の表を再現する。二つの解の差なら $\left(z - 1\right)^{2} \left(z + 1\right)^{2}$ で割り切れ、その次数は $4$ 以上でなければならない。したがって次数の制限内では差は零である。
:::

:::example[点ごとの Taylor 条件をまとめて解く]
多項式 $f$ の次数を $5$ 以下に制限し、以下の Taylor 係数を実現する。階乗で割った導関数を係数として使う。

$$
\begin{aligned}
a=0,\quad j=0&\colon\quad f^{(0)}(a)/0!= 1 \\
a=1,\quad j=0&\colon\quad f^{(0)}(a)/0!= 0 \\
a=1,\quad j=1&\colon\quad f^{(1)}(a)/1!= 0 \\
a=2,\quad j=0&\colon\quad f^{(0)}(a)/0!= 2 \\
a=2,\quad j=1&\colon\quad f^{(1)}(a)/1!= 3 \\
a=2,\quad j=2&\colon\quad f^{(2)}(a)/2!= 1
\end{aligned}
$$

$f=\sum_k c_kz^k$ とし、$z^k$ の $a$ 周りの $j$ 次係数 $\binom kj a^{k-j}$ を並べる。$k<j$ のときは零であり、$a=0$ でも式は多項式の係数として読める。

$$
\begin{aligned}
A&=\left[\begin{matrix}1 & 0 & 0 & 0 & 0 & 0\\1 & 1 & 1 & 1 & 1 & 1\\0 & 1 & 2 & 3 & 4 & 5\\1 & 2 & 4 & 8 & 16 & 32\\0 & 1 & 4 & 12 & 32 & 80\\0 & 0 & 1 & 6 & 24 & 80\end{matrix}\right] \\
Ac&=\left[\begin{matrix}1\\0\\0\\2\\3\\1\end{matrix}\right] \\
\det A&=8 \\
c&=\left[\begin{matrix}1\\\frac{7}{2}\\- \frac{57}{4}\\\frac{119}{8}\\-6\\\frac{7}{8}\end{matrix}\right] \\
f(z)&=\frac{7 z^{5}}{8} - 6 z^{4} + \frac{119 z^{3}}{8} - \frac{57 z^{2}}{4} + \frac{7 z}{2} + 1
\end{aligned}
$$

この多項式の Taylor 展開は指定の表を再現する。二つの解の差なら $z \left(z - 2\right)^{3} \left(z - 1\right)^{2}$ で割り切れ、その次数は $6$ 以上でなければならない。したがって次数の制限内では差は零である。
:::

## Serre 双対性

:::theorem[コンパクト複素多様体の Serre 双対性]{#thm-cg-serre-duality}
複素次元 $n$ のコンパクト複素多様体 $X$ と正則ベクトル束 $E$ に対し、外積と $E\otimes E^*\to\C$ の評価を積分した
$$
H^q(X,\mathcal O(E))\times H^{n-q}(X,\mathcal O(K_X\otimes E^*))
 \longrightarrow\C
$$
は有限次元空間の完全双線形ペアリングである。
:::

:::proof
:ref[thm-cg-global-dolbeault] で両側を $(0,q)$ 形式と $(n,n-q)$ 形式で表す。外積の Leibniz 則により、一方が $\bar\partial$ の像で他方が閉なら積は $\bar\partial$ の像である。積分する最高型では Stokes の定理が使え、ペアリングはコホモロジーに降りる。

滑らかな Hermite 計量を $X,E$ に選ぶ。$\bar\partial_E$ の記号は非零実余接ベクトル $\xi$ について $\xi^{0,1}$ の外積であり、補助ベクトルでの収縮によりこの記号複体は完全である。従って楕円複体の Hodge 理論により各類に一意な調和代表がある。[Hodge 理論と楕円型正則性](/math/riemannian-geometry/hodge-bochner/) を参照する。$L^2$ 内積を $\int_X\alpha\wedge S\beta$ と書くように、外積による点ごとの完全ペアリングを用いて共役線形の束同型 $S:A^{0,q}(E)\to A^{n,n-q}(E^*)$ を作る。この同型は局所正規直交枠で、補集合の添字の形式に係数の共役と向きの符号を掛ける操作である。

任意の試験形式 $\gamma$ に対し積分の部分積分を行うと、$\bar\partial^*\alpha=0$ は $\bar\partial S\alpha=0$ と同値である。逆側から同じ計算をすれば、$\bar\partial\alpha=0$ は $\bar\partial^*S\alpha=0$ と同値である。従って $S$ は調和形式の二空間を入れ替え、非零な調和 $\alpha$ は $S\alpha$ との積分で $\|\alpha\|^2>0$ を与える。これで完全性を得る。
:::

:::example[射影平面の最上位]
$K_{\mathbb P^2}=\mathcal O(-3)$ なので $H^2(\mathbb P^2,\mathcal O(-5))^*\cong H^0(\mathbb P^2,\mathcal O(2))$ である。後者の次元は三変数二次単項式の数 $6$ であり、負指数三つの絶対値の和が $5$ となる Laurent 基底の数とも一致する。
:::

双対性の部分積分では、テスト形式 $\gamma$ を一つ前の次数に取り、複素次元のため最高型の積に対して $\mathrm d$ と $\bar\partial$ が一致することを使う。例えば $\gamma$ が型 $(0,q-1)$ なら
$$
0=\int_X\bar\partial(\gamma\wedge S\alpha)
 =\langle\bar\partial_E\gamma,\alpha\rangle+
   (-1)^{q-1}\int_X\gamma\wedge\bar\partial_{E^*}S\alpha .
$$
任意の $\gamma$ に対して左の内積が零であることは、後ろの積分が全て零であることと同値である。点ごとの完全ペアリングと台を局所化できることから、これは $\bar\partial_{E^*}S\alpha=0$ と同値である。補次数側でも同じ計算を行える。最高外積のペアリングから $(\Omega_X^p)^*\otimes K_X\cong\Omega_X^{n-p}$ なので、$E$ を $\Omega_X^p\otimes E$ に置き換えると、任意の $(p,q)$ 型に対し双対側は $E^*$ 値の $(n-p,n-q)$ 型になる。

### 積分ペアリングと単項式双対

$\mathbb P^n$ の標準被覆で $H^n(\mathcal O(-n-1-k))$ の代表を $Z^{-a-\mathbf 1}$、$a_i\ge0$、$\sum a_i=k$ とする。$H^0(\mathcal O(k))$ の単項式 $Z^b$ を掛けると $Z^{b-a-\mathbf 1}$ になる。総次数は $-n-1$ であるから、全指数が負になるのは $a=b$ の場合に限る。従ってコホモロジーの積の行列は最高次の生成元 $Z_0^{-1}\cdots Z_n^{-1}$ を単位とすれば単位行列になる。積分を用いる場合はこの生成元の跡を非零定数に規格化すれば、同じ双対基底が得られる。

:::example[三つの切断を Laurent 側で識別する]
$\mathbb P^2$ で $k=1$ とする。$Z_1$ に対応する Laurent 代表は $1/(Z_0Z_1^2Z_2)$ である。$Z_1$ を掛けると三成分すべてが指数 $-1$ となるが、$Z_0$ を掛けると $Z_0$ の指数が $0$ になる。その項は最高次のコホモロジーでは零であり、$Z_2$ についても同じである。
:::

### Serre 双対の次数と Euler 標数

:::example[積の直線束でコホモロジーが残る次数]
$X=\prod_j\mathbb P^{n_j}$ とし、外部テンソル積の直線束 $L=\mathcal O(-8)$ を考える。$a$ のべきでコホモロジー次数を記録し、係数を各空間の次元とする。

一因子では符号がすべて揃った Taylor–Laurent 指数だけが残る。指数が非負の場合は全て $0$ 次に、全て負の場合は因子の最高次に属する。負指数を $-b_i-1$ と書けば、弱い合成の個数を二項係数で求められる。

$$
\begin{aligned}
j=1&\colon\quad \mathbb P^{3},\; d=-8,\quad e_i<0,\quad \sum_i e_i=-8,\quad P_j(a)=35 a^{3}
\end{aligned}
$$

積では各因子の Dolbeault Laplacian の非負固有値を加える。零になるには各因子の固有値が零でなければならず、調和代表の次元と次数は多項式の積で合成される。

$$
\begin{aligned}
P_X(a)&=\prod_j P_j(a)=35 a^{3} \\
\dim_\C X&=3 \\
\chi(X,L)&=P_X(-1)=-35
\end{aligned}
$$

随伴の双対束 $K_X\otimes L^*$ は $\mathcal O(4)$ である。$d\ge0$ の因子では双対側が最高次になり、$d\le-n-1$ の因子では逆になる。全消滅の範囲はそのままなので、双対は $q$ 次と $3-q$ 次を対応させる。最高次の双対は四次の多項式を与え、元の次数では Euler 標数に負符号が掛かる。
:::

:::example[積の直線束でコホモロジーが残る次数]
$X=\prod_j\mathbb P^{n_j}$ とし、外部テンソル積の直線束 $L=\mathcal O(-4,3)$ を考える。$a$ のべきでコホモロジー次数を記録し、係数を各空間の次元とする。

一因子では符号がすべて揃った Taylor–Laurent 指数だけが残る。指数が非負の場合は全て $0$ 次に、全て負の場合は因子の最高次に属する。負指数を $-b_i-1$ と書けば、弱い合成の個数を二項係数で求められる。

$$
\begin{aligned}
j=1&\colon\quad \mathbb P^{1},\; d=-4,\quad e_i<0,\quad \sum_i e_i=-4,\quad P_j(a)=3 a \\
j=2&\colon\quad \mathbb P^{2},\; d=3,\quad e_i\ge0,\quad \sum_i e_i=3,\quad P_j(a)=10
\end{aligned}
$$

積では各因子の Dolbeault Laplacian の非負固有値を加える。零になるには各因子の固有値が零でなければならず、調和代表の次元と次数は多項式の積で合成される。

$$
\begin{aligned}
P_X(a)&=\prod_j P_j(a)=30 a \\
\dim_\C X&=3 \\
\chi(X,L)&=P_X(-1)=-30
\end{aligned}
$$

随伴の双対束 $K_X\otimes L^*$ は $\mathcal O(2,-6)$ である。$d\ge0$ の因子では双対側が最高次になり、$d\le-n-1$ の因子では逆になる。全消滅の範囲はそのままなので、双対は $q$ 次と $3-q$ 次を対応させる。積での双対は一次から二次へ移し替える。
:::

:::example[積の直線束でコホモロジーが残る次数]
$X=\prod_j\mathbb P^{n_j}$ とし、外部テンソル積の直線束 $L=\mathcal O(1,-3,-4)$ を考える。$a$ のべきでコホモロジー次数を記録し、係数を各空間の次元とする。

一因子では符号がすべて揃った Taylor–Laurent 指数だけが残る。指数が非負の場合は全て $0$ 次に、全て負の場合は因子の最高次に属する。負指数を $-b_i-1$ と書けば、弱い合成の個数を二項係数で求められる。

$$
\begin{aligned}
j=1&\colon\quad \mathbb P^{1},\; d=1,\quad e_i\ge0,\quad \sum_i e_i=1,\quad P_j(a)=2 \\
j=2&\colon\quad \mathbb P^{2},\; d=-3,\quad e_i<0,\quad \sum_i e_i=-3,\quad P_j(a)=a^{2} \\
j=3&\colon\quad \mathbb P^{1},\; d=-4,\quad e_i<0,\quad \sum_i e_i=-4,\quad P_j(a)=3 a
\end{aligned}
$$

積では各因子の Dolbeault Laplacian の非負固有値を加える。零になるには各因子の固有値が零でなければならず、調和代表の次元と次数は多項式の積で合成される。

$$
\begin{aligned}
P_X(a)&=\prod_j P_j(a)=6 a^{3} \\
\dim_\C X&=4 \\
\chi(X,L)&=P_X(-1)=-6
\end{aligned}
$$

随伴の双対束 $K_X\otimes L^*$ は $\mathcal O(-3,0,2)$ である。$d\ge0$ の因子では双対側が最高次になり、$d\le-n-1$ の因子では逆になる。全消滅の範囲はそのままなので、双対は $q$ 次と $4-q$ 次を対応させる。因子ごとの最高次と零次を合わせると、総次数は最高次より低くなる。
:::

### 双対基底が空になる範囲

:::example[積の直線束でコホモロジーが残る次数]
$X=\prod_j\mathbb P^{n_j}$ とし、外部テンソル積の直線束 $L=\mathcal O(-1,5)$ を考える。$a$ のべきでコホモロジー次数を記録し、係数を各空間の次元とする。

一因子では符号がすべて揃った Taylor–Laurent 指数だけが残る。指数が非負の場合は全て $0$ 次に、全て負の場合は因子の最高次に属する。負指数を $-b_i-1$ と書けば、弱い合成の個数を二項係数で求められる。

$$
\begin{aligned}
j=1&\colon\quad \mathbb P^{2},\; d=-1,\quad -2\le -1<0,\quad P_j(a)=0 \\
j=2&\colon\quad \mathbb P^{2},\; d=5,\quad e_i\ge0,\quad \sum_i e_i=5,\quad P_j(a)=21
\end{aligned}
$$

積では各因子の Dolbeault Laplacian の非負固有値を加える。零になるには各因子の固有値が零でなければならず、調和代表の次元と次数は多項式の積で合成される。

$$
\begin{aligned}
P_X(a)&=\prod_j P_j(a)=0 \\
\dim_\C X&=4 \\
\chi(X,L)&=P_X(-1)=0
\end{aligned}
$$

随伴の双対束 $K_X\otimes L^*$ は $\mathcal O(-2,-8)$ である。$d\ge0$ の因子では双対側が最高次になり、$d\le-n-1$ の因子では逆になる。全消滅の範囲はそのままなので、双対は $q$ 次と $4-q$ 次を対応させる。最初の因子が消滅するため双対側の全次数も消滅する。
:::

:::example[積の直線束でコホモロジーが残る次数]
$X=\prod_j\mathbb P^{n_j}$ とし、外部テンソル積の直線束 $L=\mathcal O(-6,-4)$ を考える。$a$ のべきでコホモロジー次数を記録し、係数を各空間の次元とする。

一因子では符号がすべて揃った Taylor–Laurent 指数だけが残る。指数が非負の場合は全て $0$ 次に、全て負の場合は因子の最高次に属する。負指数を $-b_i-1$ と書けば、弱い合成の個数を二項係数で求められる。

$$
\begin{aligned}
j=1&\colon\quad \mathbb P^{3},\; d=-6,\quad e_i<0,\quad \sum_i e_i=-6,\quad P_j(a)=10 a^{3} \\
j=2&\colon\quad \mathbb P^{1},\; d=-4,\quad e_i<0,\quad \sum_i e_i=-4,\quad P_j(a)=3 a
\end{aligned}
$$

積では各因子の Dolbeault Laplacian の非負固有値を加える。零になるには各因子の固有値が零でなければならず、調和代表の次元と次数は多項式の積で合成される。

$$
\begin{aligned}
P_X(a)&=\prod_j P_j(a)=30 a^{4} \\
\dim_\C X&=4 \\
\chi(X,L)&=P_X(-1)=30
\end{aligned}
$$

随伴の双対束 $K_X\otimes L^*$ は $\mathcal O(2,2)$ である。$d\ge0$ の因子では双対側が最高次になり、$d\le-n-1$ の因子では逆になる。全消滅の範囲はそのままなので、双対は $q$ 次と $4-q$ 次を対応させる。両因子の最高次の積が残るので、双対は総零次の多項式である。
:::

## Hirzebruch–Riemann–Roch の定理の概観

$E$ の形式的 Chern 根 $x_i$ は、基本対称式が Chern 類になる記号である。Chern 指標と Todd 類を
$$
\operatorname{ch}(E)=\sum_i e^{x_i},\qquad
\operatorname{Td}(T_X)=\prod_j\frac{y_j}{1-e^{-y_j}}
$$
と置き、最高次を越える項は切り捨てる。これらは対称式なので実際に束が直線束へ分解する必要はない。Hirzebruch–Riemann–Roch の定理はコンパクト複素多様体と正則ベクトル束について
$$
\chi(X,E)=\sum_{q=0}^n(-1)^q\dim H^q(X,\mathcal O(E))
 =\int_X\operatorname{ch}(E)\operatorname{Td}(T_X)
$$
と述べる。右辺は実次数 $2n$ の部分を積分する。一般の証明は Dolbeault 作用素の解析的指数を位相的指数で計算する指数定理に及ぶため、この章では定理の適用を扱う。証明の枠組みは :cite[wells-complex] を参照する。

$(1-e^{-y})/y=1-y/2+y^2/6+\cdots$ を逆転すると $y/(1-e^{-y})=1+y/2+y^2/12+\cdots$ となる。曲線ではこの一次項までを使うため $\chi(X,L)=\int c_1(L)+c_1(T_X)/2=\deg L+1-g$ が再現される。曲面で直線束 $L$ に対しては
$$
\chi(X,L)=\int_X \tfrac12 c_1(L)^2+
 \tfrac12 c_1(L)c_1(T_X)+\tfrac1{12}\bigl(c_1(T_X)^2+c_2(T_X)\bigr).
$$
最後の項は二つの根 $y_1,y_2$ について $(y_1^2+y_2^2)/12+y_1y_2/4$ を基本対称式に直したものである。

:::example[射影平面]
超平面類を $H$ とすれば接束の完全列から $c_1(T)=3H,c_2(T)=3H^2$ であり、$\int H^2=1$ である。従って $\chi(\mathcal O(k))=(k^2+3k+2)/2$ となる。例えば $k=3$ では $10$、$k=-4$ では $3$ である。前者は $H^0$、後者は $H^2$ の次元として :ref[thm-cg-projective-line-cohomology] と一致する。
:::

### Euler 標数の加法性を完全列で検算する

短完全列 $0\to\mathcal F'\to\mathcal F\to\mathcal F''\to0$ があり全てのコホモロジーが有限次元で有限次数にしか現れなければ、長完全列の交代和から $\chi(\mathcal F)=\chi(\mathcal F')+\chi(\mathcal F'')$ となる。個々の $h^q$ にはこのような加法性はない。連結準同型が隣の次数への移動を記録するためである。例えば滑らかな次数 $d$ の超曲面 $Y\subset\mathbb P^n$ では定義式を掛ける完全列
$$
0\longrightarrow\mathcal O(k-d)\longrightarrow\mathcal O(k)
 \longrightarrow\mathcal O_Y(k)\longrightarrow0
$$
があり、$\chi(\mathcal O_Y(k))=\chi(\mathcal O(k))-\chi(\mathcal O(k-d))$ となる。

:::example[平面三次曲線の値]
$n=2,d=3$ を入れると差は
$$
\frac{(k+1)(k+2)-(k-2)(k-1)}2=3k
$$
になる。曲線の Riemann–Roch でも種数 $1$、超平面束の次数 $3$ なので $\chi=3k$ である。例えば $k=0$ では $\chi=0$ だが、$h^0=h^1=1$ であり、コホモロジー群が両方消えるわけではない。
:::

### 曲面の Riemann–Roch を交点数で使う

曲面で $c_1(T_X)=-c_1(K_X)$ と置けば、直線束の公式は
$$
\chi(L)=\chi(\mathcal O_X)+\tfrac12 c_1(L)\cdot(c_1(L)-c_1(K_X))
$$
となる。$\chi(\mathcal O_X)$ は $L=\mathcal O_X$ を入れて求めた定数項である。この式は Euler 標数だけを与えるので、個々のコホモロジー次元には消滅定理または別の計算が必要である。

:::example[二種類の次数を持つ曲面]
$X=\mathbb P^1\times\mathbb P^1$、各因子からの超平面類を $A,B$ とする。$A^2=B^2=0$、$\int AB=1$、$c_1(T_X)=2A+2B$、$c_2(T_X)=4AB$ である。$L=\mathcal O(a,b)$ なら積の最高項は
$$
ab+(a+b)+1=(a+1)(b+1)
$$
となる。例えば $(a,b)=(2,-3)$ で $\chi=-6$ となる。第一因子に三つの大域切断、第二因子に二つの一次類があり、積の Künneth の公式でも一次に六次元が集中する。
:::

### 完全交差の Chern 係数を指数式に入れる

:::example[完全交差の標準束と Chern 数を計算する]
$X\subset\mathbb P^{3}$ を次数 $4$ の横断的な滑らかな完全交差とし、$H=c_1(\mathcal O_X(1))$ とする。法束は $\bigoplus_a\mathcal O_X(d_a)$ であり、次元は $2$ である。

接束の完全列から $c(T_X)c(N_{X/\mathbb P^N})=c(T_{\mathbb P^N}|_X)$ を得る。Euler 列により右辺は $(1+H)^{N+1}$ である。$c(T_X)=1+c_1H+\cdots+c_nH^n$ と書いて係数を比べる。

$$
\begin{aligned}
c(N)&=4 H + 1 \\
c_{1}+(4)&=\binom{4}{1},\qquad c_{1}=0 \\
c_{2}+(0)&=\binom{4}{2},\qquad c_{2}=6 \\
c(T_X)&=6 H^{2} + 1
\end{aligned}
$$

法束の行列式を取った随伴公式は $K_X=\mathcal O_X(0)$ である。これは $c(T_X)$ の一次の係数の符号を反転した値にも等しい。さらに $X$ を一般の $2$ 個の超平面で切り、交点数を積で数えると $\int_XH^2=4$ となる。

$$
\begin{aligned}
\int_X c_n(T_X)&=6\cdot4=24
\end{aligned}
$$

最高 Chern 類は実接束の Euler 類であるため、ここで得た整数は位相的 Euler 標数にも等しい。この曲面では $\chi(\mathcal O_X)=(c_1^2+c_2)[X]/12=24/12=2$ である。
:::

:::example[完全交差の標準束と Chern 数を計算する]
$X\subset\mathbb P^{4}$ を次数 $2,2$ の横断的な滑らかな完全交差とし、$H=c_1(\mathcal O_X(1))$ とする。法束は $\bigoplus_a\mathcal O_X(d_a)$ であり、次元は $2$ である。

接束の完全列から $c(T_X)c(N_{X/\mathbb P^N})=c(T_{\mathbb P^N}|_X)$ を得る。Euler 列により右辺は $(1+H)^{N+1}$ である。$c(T_X)=1+c_1H+\cdots+c_nH^n$ と書いて係数を比べる。

$$
\begin{aligned}
c(N)&=\left(2 H + 1\right)^{2} \\
c_{1}+(4)&=\binom{5}{1},\qquad c_{1}=1 \\
c_{2}+(8)&=\binom{5}{2},\qquad c_{2}=2 \\
c(T_X)&=2 H^{2} + H + 1
\end{aligned}
$$

法束の行列式を取った随伴公式は $K_X=\mathcal O_X(-1)$ である。これは $c(T_X)$ の一次の係数の符号を反転した値にも等しい。さらに $X$ を一般の $2$ 個の超平面で切り、交点数を積で数えると $\int_XH^2=4$ となる。

$$
\begin{aligned}
\int_X c_n(T_X)&=2\cdot4=8
\end{aligned}
$$

最高 Chern 類は実接束の Euler 類であるため、ここで得た整数は位相的 Euler 標数にも等しい。この曲面では $c_1=H$、$c_2=2H^2$ なので $\chi(\mathcal O_X)=(1+2)\cdot4/12=1$ となる。
:::

:::example[完全交差の標準束と Chern 数を計算する]
$X\subset\mathbb P^{4}$ を次数 $3,3$ の横断的な滑らかな完全交差とし、$H=c_1(\mathcal O_X(1))$ とする。法束は $\bigoplus_a\mathcal O_X(d_a)$ であり、次元は $2$ である。

接束の完全列から $c(T_X)c(N_{X/\mathbb P^N})=c(T_{\mathbb P^N}|_X)$ を得る。Euler 列により右辺は $(1+H)^{N+1}$ である。$c(T_X)=1+c_1H+\cdots+c_nH^n$ と書いて係数を比べる。

$$
\begin{aligned}
c(N)&=\left(3 H + 1\right)^{2} \\
c_{1}+(6)&=\binom{5}{1},\qquad c_{1}=-1 \\
c_{2}+(3)&=\binom{5}{2},\qquad c_{2}=7 \\
c(T_X)&=7 H^{2} - H + 1
\end{aligned}
$$

法束の行列式を取った随伴公式は $K_X=\mathcal O_X(1)$ である。これは $c(T_X)$ の一次の係数の符号を反転した値にも等しい。さらに $X$ を一般の $2$ 個の超平面で切り、交点数を積で数えると $\int_XH^2=9$ となる。

$$
\begin{aligned}
\int_X c_n(T_X)&=7\cdot9=63
\end{aligned}
$$

最高 Chern 類は実接束の Euler 類であるため、ここで得た整数は位相的 Euler 標数にも等しい。この曲面で $c_1=-H$、$c_2=7H^2$ であり、$\chi(\mathcal O_X)=8\cdot9/12=6$ である。
:::

### 最高次だけが残る例の符号

$\mathbb P^3$ 上の $\mathcal O(-4)$ では最高次の Laurent 基底は $1/(Z_0Z_1Z_2Z_3)$ 一つである。従って $\chi(\mathcal O(-4))=(-1)^3=-1$ である。双対側の自明束は $H^0$ にだけ一つの切断を持ち、Euler 標数は $1$ になる。双対性は次元を対応させるが、複素次元が奇数なので両方の Euler 標数の符号は逆になる。 一般にも $h^q(X,E)=h^{n-q}(X,K_X\otimes E^*)$ の両辺に $(-1)^q$ を掛けて和を取ると、$\chi(X,E)=(-1)^n\chi(X,K_X\otimes E^*)$ となる。

## 演習問題

::::exercise{level=basic}
$\mathcal O_{\C,0}/(z^4)$ の複素ベクトル空間の基底を求めよ。
:::solution
Taylor 展開を三次で切って $1,z,z^2,z^3$ で表せる。これらの線形結合の位数は非零なら $3$ 以下であり、$(z^4)$ に入らないので独立である。
:::
::::

::::exercise{level=basic}
直線束 $L,M$ について $c_1(L\otimes M)=c_1(L)+c_1(M)$ を示せ。
:::solution
テンソル積の遷移関数は積である。小交差で各対数を取ると、積の対数として和を使える。三重交差でその差を取る式は加法的なので、整数 $2$ コサイクルとその類も和になる。
:::
::::

::::exercise{level=basic}
$H^2(\mathbb P^2,\mathcal O(-6))$ の次元を Serre 双対性から求めよ。
:::solution
双対は $H^0(\mathbb P^2,\mathcal O(3))$ である。三変数の斉次三次単項式の数は $\binom52=10$ であり、従って求める次元も $10$ である。
:::
::::

::::exercise{level=basic}
$\C^\times$ が正則凸であることを定義から示せ。
:::solution
コンパクト $K$ に対し $M=\max_K|z|$、$m=\min_K|z|>0$ とする。$z$ と $1/z$ を正則包の定義に使うと $\widehat K$ は $m\le |z|\le M$ に入る。正則包は閉なので、このコンパクトな円環の閉部分集合としてコンパクトである。
:::
::::

::::exercise{level=advanced}
複素三次元のトーラス上の自明束の正則 Euler 標数を求めよ。
:::solution
接束は平行移動で自明となるので Todd 類は $1$、自明束の Chern 指標も $1$ である。積の最高次成分は零であり、Riemann–Roch により $\chi=0$ である。
:::
::::

::::exercise{level=advanced}
$H^2(X,\mathcal O_X)=0$ が成り立つなら、どの整数二次コホモロジー類も正則直線束の第一 Chern 類になることを示せ。
:::solution
指数完全列では $H^2(X,\Z)\to H^2(X,\mathcal O_X)$ の核が $c_1$ の像である。終域は零なので核は全体であり、$c_1$ は全射となる。
:::
::::

::::exercise{level=advanced}
コンパクト連結な正次元複素多様体は Stein でないことを示せ。
:::solution
コンパクト連結複素多様体の正則関数は最大値原理により定数である。正次元では任意の座標領域に異なる二点があるが、全正則関数がその二点で同じ値を取る。従って Stein の点の分離条件を満たさない。
:::
::::

::::exercise{level=advanced}
複素曲面で $c_1(T_X)=0$ かつ $\int_X c_2(T_X)=24$ とする。$L$ が $\int_X c_1(L)^2=6$ を満たすとき $\chi(X,L)$ を求めよ。
:::solution
曲面の Riemann–Roch で混合項と $c_1(T_X)^2$ の項は零である。従って $\chi=6/2+24/12=3+2=5$ である。
:::
::::

::::exercise{level=advanced}
$\mathbb P^2$ の滑らかな次数 $5$ の曲線 $Y$ で、$\chi(\mathcal O_Y(2))$ を射影空間の完全列から求めよ。$h^0,h^1$ も求めて交代和を確認せよ。

:::solution
定義式を掛ける列は $0\to\mathcal O_{\mathbb P^2}(-3)\to\mathcal O_{\mathbb P^2}(2)\to\mathcal O_Y(2)\to0$ である。$\mathcal O(2)$ は零次に $6$ 次元だけを持つ。一方 $\mathcal O(-3)$ は最高次に $1$ 次元を持ち、その基底は $1/(Z_0Z_1Z_2)$ である。よって Euler 標数の差は $6-1=5$ となる。長完全列では $H^0(\mathcal O_Y(2))\cong H^0(\mathcal O(2))$、$H^1(\mathcal O_Y(2))\cong H^2(\mathcal O(-3))$ なので $h^0=6,h^1=1$ である。これは随伴で $K_Y=\mathcal O_Y(2)$ となり、曲線の種数 $6$ の標準束について $h^0=6,h^1=1$ を得ることにも一致する。
:::
::::

## 確認テスト

::::quiz{answer="6" tol="1e-6"}
$H^2(\mathbb P^2,\mathcal O(-5))$ の複素次元を求めよ。
:::explanation
Serre 双対性により双対は $H^0(\mathbb P^2,\mathcal O(2))$ である。三変数二次単項式は $\binom42=6$ 個ある。
:::
::::

::::quiz
Cartan の定理 B が主張することを選べ。

- [x] Stein 多様体上の連接層の正次数コホモロジーは零になる。
- [ ] コンパクトな複素多様体上の全ての直線束は自明である。
- [ ] 任意の層の大域切断は各点の茎を生成する。

:::explanation
定理 B には底の Stein 性と層の連接性が必要である。定理 A は同じ仮定の下で大域切断による茎の生成を述べる。
:::
::::

::::quiz{answer="10" tol="1e-6"}
$\chi(\mathbb P^2,\mathcal O(3))$ を求めよ。
:::explanation
$k=3$ を射影平面の公式 $(k^2+3k+2)/2$ に代入すると $10$ である。切断の計算でも、正次数が零で三次斉次多項式の次元が $10$ となる。
:::
::::

::::quiz
指数完全列で遷移関数の局所対数から $c_1$ が得られる理由を選べ。

- [ ] 対数が全ての交差で一意に定まる。
- [x] 三重交差での対数の和が整数となる。
- [ ] 二次コホモロジーは常に零である。

:::explanation
遷移関数の三重交差での積が $1$ なので、$2\pi i$ で規格化した対数の和は整数である。枝を替える差は整数コバウンダリであり、得られる二次コホモロジー類は不変である。
:::
::::

::::quiz{answer="-6" tol="1e-6"}
$\mathbb P^1\times\mathbb P^1$ で $\chi(\mathcal O(2,-3))$ を求めよ。

:::explanation
第一因子 $\mathcal O(2)$ では $H^0$ が三次元、第二因子 $\mathcal O(-3)$ では $H^1$ が二次元である。積のコホモロジーは一次に $3\cdot2=6$ 次元だけ残るので、交代和は $-6$ になる。
:::
::::

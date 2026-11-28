---
title: 入射加群と Gorenstein 環
description: 入射包絡と Matlis 双対を作り、標準加群と局所コホモロジーの関係を調べる。
notation:
  - symbol: 'E_A(k)'
    meaning: 剰余体 k の入射包絡
  - symbol: 'D(M)'
    meaning: Matlis 双対 Hom_A(M,E_A(k))
  - symbol: '\operatorname{id}_A M'
    meaning: 加群 M の入射次元
  - symbol: 'H^i_{\mathfrak m}(M)'
    meaning: 極大イデアルに台をもつ局所コホモロジー
---

## この章の目標

- Baer の判定から入射性を調べる。
- 本質拡大を用いて入射包絡の存在と一意性を証明する。
- 完備局所環で有限加群と Artin 加群の双対を作る。
- Gorenstein 性、標準加群、局所双対性を結び付ける。

## 入射加群と入射包絡

:::definition[入射加群と本質拡大]{#def-comm-injective-hull}
任意の包含 $N\subset M$ に対して $N\to E$ が $M\to E$ に延長できるとき、$E$ を :term[入射加群]{en="injective module" yomi="にゅうしゃかぐん"} という。$M\subset E$ が :term[本質拡大]{en="essential extension" yomi="ほんしつかくだい"} とは、$E$ のどの非零部分加群も $M$ と非零に交わることである。入射加群への本質拡大を :term[入射包絡]{en="injective hull" yomi="にゅうしゃほうらく"} という。
:::

:::theorem[Baer の判定]{#thm-comm-baer}
$E$ が入射であることと、任意のイデアル $I\subset A$ からの写像 $I\to E$ が $A\to E$ に延長できることは同値である。
:::

:::proof
必要性は定義である。逆に $N\subset M$ と $N\to E$ を取り、延長の組を包含で順序付ける。鎖の合併が上界なので、Zorn の補題により極大な延長 $f:L\to E$ がある。$x\in M\setminus L$ とすると、$I=\{a\mid ax\in L\}$ 上 $a\mapsto f(ax)$ が定まる。延長した $A\to E$ の $1$ の像を $e$ とすれば $f(\ell)+ae$ は $\ell+ax$ の表示に依存しない。この写像で $L+Ax$ に延長でき、極大性に反する。従って $L=M$ である。
:::

例えば $\Q/\Z$ は $\Z$ 上入射である。写像 $n\Z\to\Q/\Z$ の値は $n$ の像で決まり、その像を $n$ 倍で持ち上げられるからである。さらに任意の非零 Abel 群から非零写像を受ける。非零巡回部分群からは位数に応じた非零写像を作れ、入射性で延長できる。

:::proposition[入射包絡の存在と一意性]{#prop-comm-injective-envelope}
任意の $A$ 加群には入射包絡があり、元の加群上恒等となる同型を除いて一意である。
:::

:::proof
$J=\Hom_\Z(A,\Q/\Z)$ は $A$ 上入射である。随伴同型 $\Hom_A(N,J)\simeq\Hom_\Z(N,\Q/\Z)$ は $1$ での評価と $f\mapsto(n\mapsto[a\mapsto f(an)])$ で与えられるためである。全ての $M\to J$ を座標とする写像で $M$ は $J$ の直積に埋め込める。非零元を検出できることは $\Q/\Z$ の上の性質から従い、直積の入射性は座標ごとの延長から従う。

この入射加群を $Q$ と書き、$M$ を本質的に含む部分加群のうち極大な $E$ を取る。鎖の合併でも、非零元の作る巡回加群が $M$ と交わるので本質性は保たれる。次に $E$ と零で交わる $Q$ の部分加群のうち極大な $L$ を取る。$Q/L$ では $E$ が本質的である。非零の交わらない部分加群の逆像があれば $L$ を真に増やせるからである。$E\hookrightarrow Q$ を $Q/L\to Q$ に延長すると本質性によって核は零である。本質拡大の推移性により像は $M$ の本質拡大なので極大性で $E$ に等しい。従って $Q=E\oplus L$ であり、$E$ は入射である。二つの入射包絡の間は $M$ の恒等から延長し、本質性により単射である。入射性で分裂し、余分な直和因子は本質性に反するので同型となる。
:::

### 離散付値環の入射包絡

:::example[負のべきで剰余体を包む]
$A=k[[t]]$、$K=k((t))$ とする。このとき $K/A$ は $k$ の入射包絡である。任意の非零 $a\in A$ による乗法は $K/A$ で全射である。代表 $u$ を $K$ で $a$ によって割ればよい。$A$ の非零イデアルは主なので、$(a)\to K/A$ は $a$ の像を $a$ によって割って得た元を $1$ の像に選べば延長できる。従って Baer の判定で入射である。一方 $0:_{K/A}t$ は $kt^{-1}$ であり、任意の非零分数の類は負の最低指数を $-n$ とすると $t^{n-1}$ 倍でこの部分と非零に交わる。従って包含 $k\hookrightarrow K/A$ は本質的である。
:::

:::example[体上の有限次元代数の双対]
$A$ を有限次元の局所 $k$ 代数で剰余体が $k$ であるものとする。$E=\Hom_k(A,k)$ に $(a\lambda)(b)=\lambda(ab)$ で $A$ を作用させる。$\Hom_A(M,E)\simeq\Hom_k(M,k)$ である。順方向は $f\mapsto(m\mapsto f(m)(1))$、逆方向は $\phi\mapsto(m\mapsto(a\mapsto\phi(am)))$ である。従って $E$ は入射である。$\mathfrak m$ で消される双対ベクトルは $\mathfrak m$ 上で消えるものなので、その空間は $(A/\mathfrak m)^*$ に等しく一次元である。任意の非零部分加群は $\mathfrak m$ のべきが零になることから非零のソクルをもち、よってこの一次元空間と交わる。従って $E=E_A(k)$ である。
:::

### 入射加群を素イデアルごとの包絡に分ける

:::theorem[Noether 環上の入射加群の分解]{#thm-comm-injective-decomposition}
Noether 環 $A$ 上の入射加群 $E$ は、$E_A(A/P)$（$P\in\operatorname{Spec}A$）の直和に同型である。$P$ に対応する因子の個数は $\dim_{\kappa(P)}\Hom_{A_P}(\kappa(P),E_P)$ に等しい。
:::

:::proof
入射加群の任意の直和は入射である。実際、イデアル $J$ は有限生成なので、$J$ から直和への準同型は有限個の直和因子にしか非零の像をもたない。各成分を $A$ へ延長できるので、有限個の延長を合わせて Baer の判定を適用できる。

$E$ の中で、$E_A(A/P)$ に同型な部分加群の直和となっている族を極大に選ぶ。鎖の合併で直和性は保たれる。直和の非自明な有限関係があればすべての項が既に鎖の一段に含まれるからである。従って Zorn の補題が使える。この直和は入射なので包含が分裂し、補加群も入射である。補加群が非零なら随伴素を取り、$A/P$ の埋め込みを $E_A(A/P)$ からの写像へ延長する。延長は本質性により単射なので、直和に新たな因子を加えられて矛盾する。

$E_A(A/Q)$ の巡回部分加群の随伴素は $Q$ に限る。実際、埋め込まれた $A/P'$ は $A/Q$ と非零に交わり、その交わりの任意の非零元の零化イデアルは $P'$ と $Q$ の双方でなければならない。従って巡回部分加群の支持は $V(Q)$ に含まれ、$Q\not\subset P$ なら局所化後に全体が消える。$a\notin Q$ は $E_A(A/Q)$ 上単射である。もし核が非零なら $A/Q$ と交わり、$a$ が $A/Q$ 上単射であることに反する。さらに、$e$ を任意に取り $(a)\to E_A(A/Q)$ を $ra\mapsto re$ と定められる。$ra=0$ なら $a(re)=0$ なので $re=0$ である。写像を $A$ へ延長すると $e$ が $a$ 倍で書けるので、$a$ 倍は全射でもある。

$Q\subsetneq P$ なら $a\in P\setminus Q$ は局所化した加群上も可逆であるから、$\mathfrak m_{A_P}$ に消される元はない。$Q=P$ なら $E_A(A/P)$ は $A_P$ 加群であり、その中に $\kappa(P)$ が本質的に入る。もとの本質性により非零な $A_P$ 部分加群は $A/P$ と交わり、従って $\kappa(P)$ とも交わる。このためソクルはちょうど $\kappa(P)$ である。局所化は直和と交換し、ソクルは各座標で判定できるので、$P$ の因子の個数に関する式が従う。
:::

## Matlis 双対性

$(A,\mathfrak m,k)$ を Noether 局所環とし $E=E_A(k)$ とする。$D(M)=\Hom_A(M,E)$ を :term[Matlis 双対]{en="Matlis dual" yomi="まとりすそうつい"} という。これは反変完全で、非零加群を非零に送る。非零巡回部分加群 $A/I$ は $k$ に全射し、その写像を $E$ に取り延長できるからである。

:::theorem[Matlis 双対性]{#thm-comm-matlis}
$A$ が完備なら、$D$ は有限生成加群と Artin 加群の圏の間の反変同値を与える。両圏で評価 $M\to D(D(M))$ は同型である。一般の Noether 局所環でも $\End_A(E)\simeq\widehat A$ である。
:::

:::proof
$0:_E\mathfrak m=k$ である。これが $k$ より大きければ、一次元部分空間で指定の $k$ と交わらないものができる。$A$ 上有限長の $M$ では $D(k)\simeq k$ と完全性から $\ell(D(M))=\ell(M)$ である。評価は任意の非零元を検出する写像を $E$ に延長できるため単射で、有限長では同じ長さなので同型となる。

$E$ の非零巡回加群は、随伴素として $\mathfrak m$ 以外をもてない。$A/P$ が $E$ に埋め込まれると、その中に $k$ が入るが、$A/P$ の全非零元の零化イデアルは $P$ だからである。支持と随伴素の対応から巡回加群の支持は $\{\mathfrak m\}$ であり、その零化イデアルには $\mathfrak m$ のべきが入る。従って $E=\bigcup_n E_n$、$E_n=0:_E\mathfrak m^n=D(A/\mathfrak m^n)$ である。双対で $\End(E_n)$ は $\End(A/\mathfrak m^n)=A/\mathfrak m^n$ と対応し、この同型はスカラー作用で与えられる。自己準同型は $E_n$ を保つので、制限の逆極限により $\End(E)=\widehat A$ を得る。

$E$ の部分加群を双対すると $\widehat A$ の商になる。降鎖する部分加群の消滅イデアル $D(E/N)\subset\widehat A$ は昇鎖をなす。この昇鎖は Noether 性により停止する。厳密な降鎖は完全性と非零性の検出により厳密な昇鎖を与えるため、$E$ は Artin 加群である。有限生成 $M$ の双対は有限個の $E$ の直和に入るので Artin である。一方 Artin 加群の巡回部分は Noether でもあり有限長なので、どの非零部分加群にも $\mathfrak m$ に消される元がある。$0:_M\mathfrak m$ は降鎖条件から有限次元である。これを $E^r$ に入れて $M$ に延長すると本質性で単射となる。従って $D(M)$ は $\widehat A^r$ の商である。

$A$ が完備なら $D(D(A))=\widehat A=A$ であり、有限自由加群でも評価は同型である。有限表示を二重に双対して余核を取ると有限生成加群での同型が従う。Artin 加群では $D(M)$ が有限で、その評価が同型である。評価を双対した写像と $D(M)$ の評価との合成は恒等であるため、$D$ は $M$ の評価を同型に送る。完全性と零の検出により元の評価も同型である。
:::

### 有限商を双対してから極限を取る

:::example[有限次数で作用を追う]
$A=k[[t]]$、$E=K/A$ とする。$D(A/(t^3))$ は $t^{-3},t^{-2},t^{-1}$ の類が張る部分加群である。$t$ は順に $t^{-2},t^{-1},0$ に送る。従って $A/(t^3)\to D(A/(t^3))$、$1\mapsto t^{-3}$ は同型である。自然な商 $A/(t^4)\to A/(t^3)$ を双対すると、$E$ の中の包含となる。上の基底で同定したときは $t^3$ の商から $t^4$ の商への $t$ 倍である。これは双対によって準同型の向きが逆になることも示す。
:::

:::example[双対が無限になる場面]
$D(A)=E$ は $t^{-1},t^{-2},\ldots$ が張るため $k$ 上無限次元である。任意の非零元はある $t^n$ で消えるので、もし有限生成なら全体を消す一つの $t^n$ があるはずである。実際は $t^{-n-1}$ が残るため有限生成でない。それでも各部分加群は $0:_E t^n$ または $E$ 自身に限られる。最低指数 $-n$ をもつ元は $t^{-n}$ の単元倍なので $0:_E t^n$ を生成し、指数が非有界ならその合併が全体になるためである。このように $E$ は Artin 加群になる。
:::

### 完備化の作用を有限長の段から延長する

$E=E_A(k)$ には $\widehat A$ が作用する。$e\in E$ はある $\mathfrak m^n$ で消えるので、$\hat a$ の第 $n$ 段の剰余を代表する $a$ を取って $\hat a e=ae$ と定める。代表を変えても差は $\mathfrak m^n$ に入り、$n$ を大きくしても整合性により値は変わらない。和や積は共通の段で計算すれば $A$ の作用から従う。これは :ref[thm-comm-matlis] で得た $\End_A(E)\simeq\widehat A$ の作用そのものである。同じ方法で、全ての元が極大イデアルのべきで消える任意の加群に $\widehat A$ の作用が入る。その間の $A$ 線形写像は自動的に $\widehat A$ 線形である。各元を消す段で作用を比較すればよい。

$A$ 上 Artin な加群の各巡回部分加群は有限長だから、上の方法で $\widehat A$ 加群となる。逆に、$\widehat A$ 上 Artin な加群も各巡回部分加群が有限長なので、任意の $A$ 部分加群はその作用に対して閉じる。実際、ある元への $\hat a$ の作用は一つの $a\in A$ の作用と同じである。従って $A$ と $\widehat A$ のどちらの部分加群の降鎖も停止する。よって $A$ が完備でなくても Artin $A$ 加群を有限 $\widehat A$ 加群と対応させられる。

ここで使う入射包絡も共通である。$E_{\widehat A}(k)$ を $A$ へ係数制限すると入射である。実際 $\Hom_A(-,E_{\widehat A}(k))$ は $\Hom_{\widehat A}(\widehat A\otimes_A-,E_{\widehat A}(k))$ に同型であり、完備化環の平坦性と入射性から完全である。また $E_{\widehat A}(k)$ は極大イデアル進ねじれなので、$A$ 部分加群も $\widehat A$ 部分加群である。従って $k$ の包含は $A$ 上でも本質的であり、これは $E_A(k)$ である。完備環上の Matlis 双対性を適用すれば、有限 $\widehat A$ 加群と Artin $A$ 加群の反変対応を得る。

:::example[有理関数の局所環と級数環の差]
$A=k[t]_{(t)}$ と $\widehat A=k[[t]]$ を考える。どちらの入射包絡も各有限段で $t^{-n},\ldots,t^{-1}$ に対応する元をもち、その合併である。前者では分数体の商 $k(t)/A$、後者では $k((t))/\widehat A$ である。$k(t)\to k((t))/\widehat A$ の核は $A$ である。分母を $t$ のべきに単元を掛けた形に直すと、正負の位数は両方で同じだからである。また商加群 $k((t))/\widehat A$ の任意の元は有限個の負べきの和であるため写像の像に入る。従って包絡は同型だが、その自己準同型環は $A$ ではなく完備化 $\widehat A$ になる。例えば形式級数による作用は、各負べきの元に対して有限個の係数だけで計算できる。
:::

## 入射次元

$M$ の :term[入射次元]{en="injective dimension" yomi="にゅうしゃじげん"} は入射分解の最短長であり、有限分解がなければ無限である。$\operatorname{id}_A M\le n$ は、全イデアル $I$ で $\operatorname{Ext}_A^{n+1}(A/I,M)=0$ と同値である。必要性は入射分解から従う。逆に入射分解の第 $n$ 余核へ次元移動すると、その加群への $A/I$ の一次 Ext が消える。これは $I$ からの写像が延長できることを意味し、Baer の判定で分解をそこで打ち切れる。

例えば $A=\Z$ では $0\to\Z\to\Q\to\Q/\Z\to0$ が長さ一の入射分解である。$\Q$ と $\Q/\Z$ は割り算できるので Baer の判定で入射である。一方 $\Z$ では $2\Z$ の $2$ を $1$ に送る写像を延長できないため、$\operatorname{id}_{\Z}\Z=1$ である。

### Ext で入射分解の最短長を確認する

:::example[離散付値環の分数体と剰余]
$A=k[[t]]$ の $K$ と $K/A$ は、非零元で割り算できるので Baer の判定により入射である。よって $0\to A\to K\to K/A\to0$ が長さ一の入射分解である。$k$ の自由分解に $\Hom_A(-,A)$ を使うと $\operatorname{Ext}^1_A(k,A)=A/tA=k$ なので、長さを零に短縮できない。一方 $k$ 自身についても $\operatorname{Ext}^1_A(k,k)=k$ である。大域次元が $1$ なので $\operatorname{id}_A k=1$ になる。
:::

$A$ が Noether、$S$ が乗法的集合であるとき、入射 $A$ 加群 $E$ の $S^{-1}E$ は $S^{-1}A$ 上入射である。実際、任意のイデアルはある $J\subset A$ の拡大である。$J$ は有限表示なので、$S^{-1}J\to S^{-1}E$ は Hom の局所化により $J\to E$ の $s^{-1}$ 倍として表される。$J\to E$ を $A$ へ延長して局所化すれば与えられた写像を延長できる。Baer の判定から結論が従う。従って有限長の入射分解も局所化できる。

### 最小入射分解と Bass 数

$M$ をまずその入射包絡 $E^0$ に入れ、余核を $E^1$ に本質的に埋め込む。この操作を反復すると入射分解 $E^\bullet$ が得られる。各 $E^i$ で微分の核が本質的部分加群となるので、これを最小入射分解という。Noether 環では本質拡大 $N\subset L$ は局所化後も本質的である。実際、$L_P$ の非零巡回部分加群はある $A u\subset L$ の局所化である。$N\cap Au$ は $Au$ に本質的であり、両者の随伴素は一致する。$A/Q\subset Au$ の非零の交わりから、$Q$ は $N\cap Au$ でも随伴素になるためである。その結果、有限生成な $Au$ と $N\cap Au$ は同じ支持をもち、$N\cap Au$ の局所化は非零である。従って $L_P$ の非零巡回部分加群は $N_P$ と交わる。

$E^\bullet_P$ は入射分解であり、微分の核が本質的であるため各微分はソクル上で零である。実際ソクルの任意の非零ベクトルが張る単純部分加群は、本質的な核と交わるので核に含まれる。従って
$$
\dim_{\kappa(P)}\operatorname{Ext}_{A_P}^i(\kappa(P),M_P)
=\dim_{\kappa(P)}\Hom_{A_P}(\kappa(P),E_P^i)
$$
である。この数を :term[Bass 数]{en="Bass number" yomi="ばすすう"} と呼ぶ。:ref[thm-comm-injective-decomposition] により、これは最小入射分解の第 $i$ 項に含まれる $E_A(A/P)$ の個数である。例えば $A=k[[t]]$ の分解 $0\to A\to K\to K/A\to0$ は、$A$ が $K$ に本質的に入るため最小である。分数 $u$ の生成する任意の非零部分加群に分母を掛けると $A$ の非零元を得るからである。従って次数零には $(0)$ の因子が一つ、次数一には $(t)$ の因子が一つ現れる。剰余体の Ext だけを調べると後者だけが見える。

### 剰余体だけを用いて入射次元を判定する

:::proposition[有限加群の入射次元の剰余体判定]{#prop-comm-injective-dimension-residue}
Noether 局所環 $(A,\mathfrak m,k)$ 上の有限加群 $M$ について、$\operatorname{id}_A M\le d$ と、全ての $i>d$ で $\operatorname{Ext}_A^i(k,M)=0$ であることは同値である。ただし $d\ge0$ とする。
:::

:::proof
必要性は入射分解から従う。逆向きには、任意の有限加群 $N$ に対する $\operatorname{Ext}_A^i(N,M)$ が $i>d$ で消えることを $\dim N$ に関する帰納法で示す。$N$ が有限長なら組成因子はすべて $k$ で、短完全列の長完全列に仮定を適用して結論を得る。

一般には、有限加群の素イデアル商による有限フィルトレーションを取る。各因子の次元は $\dim N$ 以下なので、今回の次元での証明には $N=A/P$ の場合だけを示せばよい。$\dim A/P>0$ なら $a\in\mathfrak m\setminus P$ を取る。$A/P$ は整域なので $0\to A/P\xrightarrow{a}A/P\to C\to0$ は完全である。$C=A/(P+(a))$ の任意の素鎖の先頭には $P$ を足せるので $\dim C<\dim A/P$ である。帰納法によって $i>d$ で $\operatorname{Ext}_A^{i+1}(C,M)=0$ であり、長完全列により $\operatorname{Ext}_A^i(A/P,M)$ 上の $a$ 倍は全射である。

各項を有限自由とした分解で Ext を計算でき、$M$ は有限なので、この Ext も有限加群である。従って $a\in\mathfrak m$ と中山の補題により零となる。これで帰納が閉じる。特に $N=A/I$ の全てで消滅するので、Baer の判定と次元移動から $\operatorname{id}_A M\le d$ である。
:::

### Bass の公式で入射次元を求める

:::theorem[有限入射次元の Bass の公式]{#thm-comm-bass-formula}
Noether 局所環 $(A,\mathfrak m,k)$ 上の非零有限加群 $M$ が有限入射次元をもつなら、$\operatorname{id}_A M=\operatorname{depth}A$ である。
:::

:::proof
$s=\operatorname{id}_A M$ とする。$s>0$ なら、:ref[prop-comm-injective-dimension-residue] を $s-1$ に適用して $\operatorname{Ext}_A^s(k,M)\ne0$ を得る。実際、それより高い次数は入射次元の定義により零である。$s=0$ でも $\Hom_A(k,M)\ne0$ である。有限非零加群の深さは有限であるから、:ref[thm-comm-depth-ext] によりどこかの次数で Ext が非零であり、$M$ が入射ならその次数は零でなければならない。

$t=\operatorname{depth}A$ とし、長さ $t$ の $A$ 正則列 $\boldsymbol x$ で割った商を $B$ とする。列が極大なので $B$ の深さは零であり、ある単射 $k\hookrightarrow B$ がある。余核を $C$ とすれば、長完全列の一部は
$$
\operatorname{Ext}_A^s(B,M)\longrightarrow\operatorname{Ext}_A^s(k,M)
\longrightarrow\operatorname{Ext}_A^{s+1}(C,M)=0
$$
となる。中央が非零なので $\operatorname{Ext}_A^s(B,M)\ne0$ である。一方 $B$ の Koszul 自由分解は長さ $t$ なので $s\le t$ となる。同じ分解の最後の双対微分の余核は $M/(\boldsymbol x)M$ であり、中山の補題により非零である。従って $\operatorname{Ext}_A^t(B,M)\ne0$ で $t\le s$ となる。ここで列が正則でなければならないのは母環 $A$ 上であり、$M$ 上の正則性は必要ない。
:::

例えば $A=k[[x,y]]$ では大域次元が $2$ なので、全ての非零有限加群の入射次元は有限である。Bass の公式により $A$、$A/(x)$、$k$ の入射次元はどれも $2$ である。一方、それぞれの射影次元は $0,1,2$ である。前二者は長さ零と一の自由分解から、後者は剰余体の Koszul 分解から分かる。入射次元が有限な場合、その値を決めるのは環の深さであり、射影次元の公式に出る加群自身の深さとは異なる。

## Gorenstein 環

局所 Noether 環 $A$ が自身の上有限入射次元をもつとき :term[Gorenstein 環]{en="Gorenstein ring" yomi="ごれんすたいんかん"} と呼ぶ。一般の Noether 環では全ての素局所化でこの条件を課す。正則局所環 $A$ は大域次元が $d=\dim A$ なので入射次元が高々 $d$ である。正則パラメータの Koszul 分解を双対すると、$\operatorname{Ext}_A^d(k,A)=k$ である。最後の余核はパラメータ生成イデアルによる商だからであり、入射次元はちょうど $d$ である。

:::proposition[Artin 局所環の Gorenstein 判定]{#prop-comm-artin-gorenstein}
Artin 局所環 $(A,\mathfrak m,k)$ について、Gorenstein 性、$A$ 自身の入射性、$\dim_k(0:_A\mathfrak m)=1$ は同値である。
:::

:::proof
Artin 局所環は極大イデアルが冪零なので完備である。入射次元が有限なら、有限加群 $D(A)$ の射影次元が有限である。実際、各項が有限自由である $D(A)$ の最小自由分解を取り、その双対を $A=D(D(A))$ の入射分解とする。この時点では分解の長さは仮定しない。$\Hom_A(k,-)$ を適用すると $D(k\otimes_A F_\bullet)$ になるので、十分高次数で $\operatorname{Tor}_i^A(k,D(A))=0$ である。最小自由分解はこの段で停止する。Auslander–Buchsbaum の公式は $\operatorname{depth}A=0$ より $D(A)$ が自由であることを示す。$\ell(D(A))=\ell(A)$ なので階数一であり、$A\simeq D(A)=E$ となる。従って入射である。

入射な $A$ のソクル（極大イデアルに消される部分）を $k^r$ とすれば、包絡は $E^r$ である。実際 $k^r\subset E^r$ は、任意の非零巡回部分加群が有限長なのでソクルと交わるため本質的である。$k^r\subset A$ も本質的で、入射包絡の一意性と長さの比較から $r=1$ である。逆にソクルが一次元なら $A$ は $E$ に埋め込める。同じ長さなので同型となる。
:::

例えば $k[t]/(t^n)$ のソクルは $kt^{n-1}$ なので Gorenstein である。$k[x,y]/(x,y)^2$ のソクルは $kx+ky$ なので Gorenstein でない。有限自己入射次元から Cohen–Macaulay 性を得る一般の議論は、剰余体による入射次元判定と正則元による商を用いて後で証明する。

### 正則元による商への移行

:::proposition[正則元による Ext の移動]{#prop-comm-rees-ext}
$x$ を $A$ の非零因子、$B=A/(x)$ とする。任意の $B$ 加群 $M$ に対して $\operatorname{Ext}_B^i(M,B)\simeq\operatorname{Ext}^{i+1}_A(M,A)$ が成り立つ。従って局所 Gorenstein 環を極大イデアル内の正則元で割った商も Gorenstein である。
:::

:::proof
$A$ の入射分解 $A\to I^\bullet$ を取り、$J^\bullet=\Hom_A(B,I^\bullet)$ と置く。各項は $B$ 上入射である。$B$ 加群 $N$ について $\Hom_B(N,J^i)=\Hom_A(N,I^i)$ が成り立つからである。この複体のコホモロジーは $\operatorname{Ext}^\bullet_A(B,A)$ であり、$x$ による長さ一の自由分解で計算すると次数 $1$ の $B$ に限られる。

$M$ の $B$ 上の自由分解 $F_\bullet$ と $J^\bullet$ に $\Hom_B(F_\bullet,J^\bullet)$ を適用し、全次数を自由分解の添字と入射分解の添字の和とする。入射性により自由分解の方向からコホモロジーを取れば $\Hom_B(M,J^\bullet)=\Hom_A(M,I^\bullet)$ だけが残る。自由性によりもう一方の方向から取れば、次数 $1$ に $\Hom_B(F_\bullet,B)$ だけが残る。各全次数に有限個の項しかないので、二重複体の二つのフィルトレーションにより所望の同型を得る。$A$ の自己入射次元が有限なら、同型により $B$ の自己入射次元も有限である。
:::

特に正則局所環は Gorenstein であるため、正則列の各段でこの命題を使うと局所完全交叉は Gorenstein である。例えば $k[[x,y]]/(x^2,y^3)$ では、正則環の自己入射次元は $2$ であり、二つの正則元で割って自己入射次元 $0$ になる。これはソクルの計算とも一致する。

### ソクルが一次元になる積を具体的に調べる

:::example[乗法と双対の非退化な組]
$A=k[x,y]/(x^2,y^3)$ とし、$\lambda(a)$ を $xy^2$ の係数とする。$1,x,y,xy,y^2,xy^2$ の各単項式 $x^iy^j$ に対して補う単項式 $x^{1-i}y^{2-j}$ を掛けると、$\lambda$ の値は $1$ になる。他の基底元には $0$ になるので、積による組 $(a,b)\mapsto\lambda(ab)$ は非退化である。従って $a\mapsto(b\mapsto\lambda(ab))$ は $A$ を $\Hom_k(A,k)$ に同型に送る。後者が入射なので $A$ は Gorenstein であり、ソクル $kxy^2$ での判定と一致する。
:::

:::example[長さだけから Gorenstein 性は決まらない]
$A=k[t]/(t^3)$ と $B=k[x,y]/(x,y)^2$ の長さは双方 $3$ である。しかし $A$ のソクルは $kt^2$、$B$ のソクルは $kx\oplus ky$ である。従って $A$ は Gorenstein で $B$ はそうでない。$B$ でどの線形汎関数 $\lambda$ を選んでも、$kx\oplus ky$ の中で $\lambda$ に消される非零元 $u$ があり、$b u$ は $b$ の定数項倍の $u$ になる。よって $\lambda(bu)=0$ となって上の組は退化する。
:::

### 正則元による商から Gorenstein 性を戻す

Noether 局所環 $A$ と $\mathfrak m$ 内の非零因子 $x$ に対し、$A/(x)$ が Gorenstein なら $A$ も Gorenstein である。$B=A/(x)$ の自己入射次元を $d$ とする。任意のイデアル $I$ に対し、まず第二変数が $B$ の $\operatorname{Ext}_A^j(A/I,B)$ を考える。$A/I$ の $A$ 上自由分解 $F_\bullet$ を $B$ にテンソルすると、ホモロジーは次数 $0,1$ にしかない。実際、$B$ の $A$ 上の自由分解は $0\to A\xrightarrow{x}A\to B\to0$ である。この $B$ 上自由な複体と $B$ の長さ $d$ の入射分解を Hom で組み合わせる。各入射項への Hom は完全なので、$F_\bullet\otimes_A B$ の方向から先に取れば二つの次数 $0,1$ だけが残り、全次数 $d+1$ より上は消える。逆の方向では各自由項からの Hom は完全なので、全コホモロジーは $\operatorname{Ext}_A^\bullet(A/I,B)$ になる。これで $j>d+1$ における消滅を得る。

次に $0\to A\xrightarrow{x}A\to B\to0$ の Ext 長完全列で $j>d+1$ を取ると、$\operatorname{Ext}_A^j(A/I,A)$ 上の $x$ 倍が全射である。これは有限加群である。Noether 性により $A/I$ の各自由分解項を有限自由に選べ、その Hom 複体の部分商が Ext になるからである。中山の補題でこの Ext は零となる。全イデアルに対する入射次元の判定を使えば $\operatorname{id}_A A\le d+1$ である。下界については、有限入射次元 $d$ が正ならある $B$ イデアル $J$ で $\operatorname{Ext}^d_B(B/J,B)\ne0$ である。$d=0$ なら $B$ 自身からの Hom が非零である。:ref[prop-comm-rees-ext] で次数を一つ移すと $\operatorname{id}_A A\ge d+1$ となる。従って自己入射次元はちょうど一つ増える。

:::example[還元後のソクルから判定する]
$A=k[[x,y]]/(xy)$ は $x+y$ が非零因子で、$A/(x+y)=k[[x]]/(x^2)$ になる。この Artin 環のソクルは $k\bar x$ なので自己入射次元は $0$ である。持ち上げの判定より $A$ の自己入射次元は $1$ となる。$A$ は正則でないが Gorenstein である。
:::

### 有限自己入射次元から Cohen–Macaulay 性を得る

:::theorem[Bass の自己入射次元の定理]{#thm-comm-bass-self-injective}
Noether 局所環 $A$ が有限自己入射次元をもつことと、$A$ が Cohen–Macaulay で $\dim_k\operatorname{Ext}_A^{\dim A}(k,A)=1$ であることは同値である。このとき $\operatorname{id}_A A=\dim A$ である。
:::

:::proof
まず $A$ が Cohen–Macaulay であるとする。$d=\dim A$ 本のパラメータは正則列であり、その商を $B$ とする。:ref[prop-comm-rees-ext] を各段に適用すると
$$
\operatorname{Ext}_A^d(k,A)\simeq\Hom_B(k,B)=0:_B\mathfrak m_B
$$
となる。従って左辺が一次元であることは $B$ が Gorenstein であることと同値である。$B$ の自己入射次元はこのとき零である。正則元による商からの持ち上げにより、$A$ も Gorenstein で、各元を戻すたびに自己入射次元が一つ増える。よって $\operatorname{id}_A A=d$ である。

逆向きには $s=\operatorname{id}_A A<\infty$ とする。$\dim A=0$ の場合は Artin 局所環の判定が結論を与える。$\dim A>0$ で深さが零なら非零の $k\hookrightarrow A$ が存在する。まず $s=0$ なら $A$ は入射なので、この包含は $E_A(k)\to A$ へ延長できる。延長は本質性により単射で、$E_A(k)$ が入射なので分裂する。$A$ 加群 $A$ の直和分解の射影は冪等元による乗法である。局所環の冪等元は $0,1$ しかないため、$A\simeq E_A(k)$ となる。ところが $E_A(k)$ の全ての元は極大イデアルのあるべきで消える。$A$ の $1$ もそうなるので $\mathfrak m$ は冪零となり、次元が正であることに反する。

$s>0$ の場合、$0\to k\to A\to C\to0$ の Ext 長完全列から
$$
\operatorname{Ext}_A^s(A,A)\longrightarrow\operatorname{Ext}_A^s(k,A)
\longrightarrow\operatorname{Ext}_A^{s+1}(C,A)
$$
を取る。左端は $A$ が自由なので零、右端は入射次元から零である。従って $\operatorname{Ext}_A^s(k,A)=0$ であり、次数 $s$ より上も消える。:ref[prop-comm-injective-dimension-residue] で $\operatorname{id}_A A\le s-1$ となるが、これは矛盾である。ゆえに正の次元なら深さは正である。

そこで $x\in\mathfrak m$ を非零因子に選ぶ。有限自己入射次元は $A/xA$ に移る。$x$ は全極小素を避けるため、商の鎖の先頭には元の環で素を一つ足せる。また主イデアルで割った商の次元は $\dim A-1$ 以上である。従って $\dim(A/xA)=\dim A-1$ である。次元による帰納法で商は Cohen–Macaulay である。その深さが $A$ より一つ低いので $A$ も Cohen–Macaulay となる。そこで最初のパラメータ系による同型を使い、Artin 商へ Gorenstein 性を移すと $\operatorname{Ext}_A^d(k,A)$ は一次元と分かる。
:::

:::corollary[Gorenstein 環の最小入射分解]{#cor-comm-gorenstein-bass-numbers}
Gorenstein 局所環 $A$ の最小自己入射分解の第 $i$ 項には、高さ $i$ の素 $P$ の $E_A(A/P)$ がそれぞれ一回現れ、ほかの因子は現れない。
:::

:::proof
入射分解を局所化すると $A_P$ も有限自己入射次元をもつ。従って $e=\dim A_P=\operatorname{ht}P$ とすれば、この局所環の深さも入射次元も $e$ である。剰余体との Ext は深さから $e$ 未満で消え、入射次元から $e$ より大きい次数で消え、次数 $e$ では :ref[thm-comm-bass-self-injective] により一次元になる。Bass 数による入射包絡の個数の公式を使えば主張が従う。
:::

### 完備化と Gorenstein 性

:::proposition[完備化での Gorenstein 判定]{#prop-comm-gorenstein-completion}
Noether 局所環 $A$ が Gorenstein であることと、$\mathfrak m$ 進完備化 $\widehat A$ が Gorenstein であることは同値である。
:::

:::proof
$A$ と $\widehat A$ の次元は Hilbert–Samuel 多項式から等しい。有限自由分解で $\operatorname{Ext}_A^i(k,A)$ を計算して平坦な $\widehat A$ で基底変換すると、$\operatorname{Ext}_{\widehat A}^i(k,\widehat A)$ が得られる。実際各自由項は有限なので Hom がテンソルと交換し、平坦性で核と像も交換する。忠実平坦性から各次数の消滅は同値なので、深さも等しい。従って Cohen–Macaulay 性は一方から他方へ移る。

どちらかが Gorenstein であれば、:ref[thm-comm-bass-self-injective] により双方が Cohen–Macaulay となる。$A$ のパラメータ系 $\boldsymbol x$ の像は $\widehat A$ でもパラメータ系になる。$A/(\boldsymbol x)$ が有限長だから、その完備化 $\widehat A/(\boldsymbol x)\widehat A$ は元の商と同型である。双方で $\boldsymbol x$ は正則列である。従って正則元の商に沿う上昇・下降を繰り返せば、双方の Gorenstein 性は同じ Artin 商の Gorenstein 性に等しい。
:::

## 標準加群

$d$ 次元の Cohen–Macaulay 局所環 $A$ が $n$ 次元正則局所環 $R$ の商として与えられるとき、$\omega_A=\operatorname{Ext}_R^{n-d}(A,R)$ を :term[標準加群]{en="canonical module" yomi="ひょうじゅんかぐん"} と呼ぶ。$A$ が完備でない場合には、$\widehat A\otimes_A\omega$ が $\widehat A$ の標準加群と同型になる有限 $A$ 加群 $\omega$ が存在するとき、それを標準加群と呼ぶ。完備の場合に表示によらない理由は、次節の局所双対性で $\omega_A$ を $H^d_{\mathfrak m}(A)$ の双対として表せることにある。

完全交叉 $A=R/(\boldsymbol f)$ で正則列の長さを $c$ とすると、$d=n-c$ である。Koszul 分解に $\Hom_R(-,R)$ を適用すると、最高次数 $c$ の余核は $R/(\boldsymbol f)$ である。従って $\omega_A\simeq A$ となる。例えば $k[[x,y]]/(x^2,y^3)$ の標準加群は自身であり、ソクルは $kxy^2$ である。

### 標準加群の定義に使う Ext の次数

$A=R/J$ が $d$ 次元 Cohen–Macaulay 局所環で、$R$ が $n$ 次元正則局所環とする。$c=n-d$ と置くと
$$
\operatorname{Ext}_R^j(A,R)=0\qquad(j\ne c)
$$
である。上側の消滅は、$R$ 上の $A$ の射影次元が $c$ であることから従う。実際、$R$ と $A$ の正則列は $A$ への作用を通じて同じ条件を与え、$\operatorname{depth}_R A=d$ であるので、Auslander–Buchsbaum の公式を使える。下側の消滅は、$J$ の grade が $\dim R-\dim A=c$ で、grade が $\operatorname{Ext}_R^j(R/J,R)$ の最初の非零次数に等しいことから従う。ここで $J$ の高さと grade を結ぶのは母環 $R$ の Cohen–Macaulay 性である。従って Ext の集中は局所双対性を使う前にも得られる。

$A$ が完備なら、:ref[thm-comm-cohen-structure] でこの有限表示を選べる。$A$ が完備でなくても $R$ の商として与えられているときは、$\widehat R$ は正則で $\widehat A=\widehat R/J\widehat R$ となる。有限自由分解で Ext を計算して完備化をテンソルすれば
$$
\widehat A\otimes_A\operatorname{Ext}_R^c(A,R)
\simeq
\operatorname{Ext}_{\widehat R}^c(\widehat A,\widehat R)
$$
である。$\operatorname{Ext}_R^c(A,R)$ には $A$ 加群構造がある。$J$ は第一変数 $A$ を消し、そのスカラー作用は第二変数の作用と一致するからである。完備化が平坦で、自由分解の各項が有限であるため、Hom とそのコホモロジーの両方をテンソルと交換できる。従って表示から得た加群は完備化による定義も満たす。

同じ有限 $A$ 加群かどうかは、この状況では完備化で判定できる。有限加群 $U,V$ の完備化が同型なら、その同型は $\Hom_A(U,V)\otimes_A\widehat A$ の元である。$U$ は有限表示なので Hom の基底変換が成り立つためである。$\mathfrak m$ を法とした像が同じになる $A$ 上の準同型 $f:U\to V$ を選べる。実際、テンソルの係数を剰余体で近似すればよい。$f$ は剰余体上の全射と中山の補題で全射である。完備化した $f$ は、もとの同型と合成して有限 $\widehat A$ 加群の全射自己準同型に直せる。Noether 加群の全射自己準同型は単射である。各べきの核の昇鎖が停止する段で、全射性によって最初の核も零になるからである。最後に忠実平坦性によって $\Ker f=0$ となる。よって完備化による標準加群の同型類は一意である。

### 標準加群を行列の転置で求める

:::example[二つの生成元を要する標準加群]
$R=k[[x,y]]$、$A=R/(x^2,xy,y^2)$ とする。自由分解
$$
0\longrightarrow R^2\xrightarrow{\left(\begin{smallmatrix}-y&0\\x&-y\\0&x\end{smallmatrix}\right)}
R^3\xrightarrow{(x^2,xy,y^2)}R\longrightarrow A\longrightarrow0
$$
を使う。実際 $x^2a+xyb+y^2c=0$ なら法 $x$ で $c$ が $x$ の倍数と分かり、$c=xc'$ と置けば $xa+yb+y^2c'=0$ となる。すると $b+yc'=xb'$、$a=-yb'$ と書けるので、二つの列が全関係を生成する。列の独立性は第一、第三座標から従う。標準加群は $\operatorname{Ext}_R^2(A,R)$ なので転置写像の余核である。転置の成分はすべて極大イデアル内にあり、余核を剰余体で割った次元は $2$ である。従って標準加群は $A$ と同型にならない。これはソクル次元 $2$ と一致する。
:::

### 標準加群の最小生成元数をソクルで測る

完備 Cohen–Macaulay 局所環 $A$ とそのパラメータ系 $\boldsymbol x$ に対し、$\omega_A/\boldsymbol x\omega_A\simeq \omega_{A/(\boldsymbol x)}$ である。実際、一元ごとの商で次元が一つ下がり、母環上で Ext の非零な次数が一つ上がる。$0\to A\xrightarrow{x_1}A\to A/x_1A\to0$ に母環上の $\operatorname{Ext}(-,R)$ を適用すると、余核が次の標準加群となる。この操作を繰り返す。Artin 商 $B=A/(\boldsymbol x)$ では $\omega_B\simeq D_B(B)$ である。これを先に確かめるため、母環 $R$ の正則パラメータ系の十分高いべきで生成され、$B$ を消すイデアル $Q$ を選ぶ。$Q$ は正則列で生成されるので :ref[prop-comm-rees-ext] を繰り返すと、$\omega_B=\operatorname{Ext}_R^{\dim R}(B,R)\simeq\Hom_{R/Q}(B,R/Q)$ である。Artin 局所環 $R/Q$ は Gorenstein であり、自身が剰余体の入射包絡になる。従ってこの Hom は $R/Q$ 上の Matlis 双対である。$B$ 上の双対とも一致する。実際、$R/Q\to B$ の核に消される入射包絡の部分加群は、Hom の随伴で $B$ 上入射であり、もとのソクルを本質的に含むからである。双対すると $D_B(\omega_B/\mathfrak m_B\omega_B)$ は $0:_B\mathfrak m_B$ に同型になる。従って
$$
\mu_A(\omega_A)=\dim_k(0:_B\mathfrak m_B)
$$
となる。この数を Cohen–Macaulay 型という。左辺から、どのパラメータ系で割っても同じ数になる。一元ずつの Gorenstein 性の下降・上昇を適用すると、型 $1$ であることと $A$ が Gorenstein であることは同値である。

### 一次元整域でも型が異なる例

:::example[三つの連続した指数で生成される環]
$A=k[[t^3,t^4,t^5]]$ を $k[[t]]$ の定数項と次数 $3$ 以上の項だけをもつ級数の環として考える。$t^3,t^4,t^5$ はこの環を完備な級数環の商として生成する。実際、指数 $3,4,5$ に $3$ を繰り返し足せば全ての次数 $3$ 以上を得る。$k[[t^3]]$ 加群としては $1,t^4,t^5$ が基底である。指数を法 $3$ で分ければ生成性と独立性が分かる。従って有限整拡大であり、$A$ は次元 $1$ の Noether 局所整域である。$t^3$ は非零因子なので深さは $1$ である。

$t^3$ で割った商は $1,\overline{t^4},\overline{t^5}$ を基底とする。積 $t^8,t^9,t^{10}$ はそれぞれ $t^3$ と $t^5,t^6,t^7$ の積なので、二つの非定数基底元の積はすべて零である。従って Artin 商のソクルは二次元であり、型は $2$ になる。この環は Cohen–Macaulay だが Gorenstein ではない。$k[[t^2,t^3]]$ では $k[[t^2]]$ 上の基底が $1,t^3$ であり、$t^2$ で割った商のソクルは $\overline{t^3}$ で張られる。従って後者の型は $1$ である。整域性や次元だけでは二例を区別できない。
:::

## 局所コホモロジー入門と局所双対性

$I$ 進ねじれ $\Gamma_I(M)=\bigcup_n(0:_M I^n)$ の右導来関手を :term[局所コホモロジー]{en="local cohomology" yomi="きょくしょこほもろじー"} $H^i_I(M)$ と呼ぶ。この関手は $\varinjlim_n\operatorname{Ext}_A^i(A/I^n,M)$ で計算できる。入射分解上で $\Gamma_I$ が各 $\Hom(A/I^n,-)$ の合併に等しく、整列した直極限が有限関係を後の段で検査できることから完全であるためである。$I=(x)$ で $x$ が環の非零因子なら、$A/(x^n)$ の長さ一の自由分解を使える。次数一の遷移は $M/x^nM\to M/x^{n+1}M$ の $x$ 倍であり、$m\bmod x^nM\mapsto m/x^n$ により極限は $\operatorname{coker}(M\to M_x)$ である。実際局所化の等号を分母で払えば核は後段で消え、任意の分数を表せる。

:::theorem[完備局所環の局所双対性]{#thm-comm-local-duality}
$R$ を次元 $n$ の完備正則局所環、$A=R/J$、$M$ を有限 $A$ 加群とする。$D_A$ は $A$ の Matlis 双対とする。このとき

$$
D_A(H^i_{\mathfrak m_A}(M))\simeq \operatorname{Ext}_R^{n-i}(M,R).
$$

:::

### Čech 複体で計算する

$I=(x_1,\ldots,x_s)$ の Čech 複体の次数 $j$ は、添字集合 $S$（$|S|=j$）に対する $M_{\prod_{i\in S}x_i}$ の直和である。微分は添字を一つ増やす局所化写像に、増やした位置に応じた交代符号をつける。第二微分では、二つの添字を増やす二通りの符号が逆なので消える。次数零の核は各 $x_i$ のべきで消される元であり、有限個の $i$ に対するべきをまとめれば $\Gamma_I(M)$ となる。

この複体は入射加群で正次数が消える。まず Noether 環上で、入射加群 $E$ の $\Gamma_I(E)$ は再び入射である。実際、その入射包絡は入射性で $E$ に埋め込まれる。包絡内の任意の巡回部分加群の随伴素 $P$ について、対応する $A/P$ は本質性により $I$ 進ねじれ部分と非零に交わる。$A/P$ の非零元は全て零化イデアル $P$ をもつので $I\subset P$ である。支持と随伴素の対応から、その巡回部分加群の支持は $V(I)$ に含まれる。有限生成性から $I$ のべきが巡回部分加群を消す。よって入射包絡は全て $\Gamma_I(E)$ に入り、これ自身が入射である。

$I=(x)$ の場合、$E$ を $\Gamma_I(E)$ と入射な直和補空間 $F$ に分ける。$F$ には $x$ 進ねじれがないため、$ax=0$ なら任意の $f\in F$ について $x(af)=0$ より $af=0$ である。そこで $(x)\to F$、$ax\mapsto af$ はよく定まり、Baer の判定により延長した写像の $1$ の像で $f$ は $x$ 倍に書ける。従って $x$ 倍は $F$ で全単射であり、$E\to E_x$ は全射である。多元の場合は Čech 複体を一元ごとの二項複体で順に作る。最初の二項方向のコホモロジーは、局所化の完全性により次数零の $\Gamma_{(x_1)}(E)$ の局所化に限られる。これは入射加群であるため、残りの元に同じ議論を適用し、正次数は消える。

局所化は完全なので、入射分解を Čech 複体に入れた二重複体では、入射分解の方向を先に計算すると元の加群の Čech 複体を得る。逆方向を先に計算すると入射分解に $\Gamma_I$ を適用した複体を得る。従って Čech 複体のコホモロジーは導来関手で定義したものに一致する。特に商 $R\to A$ と $A$ 加群では、母環のイデアル $I$ と商での拡大 $IA$ に関する複体は同じ分数の列になる。従って二つの環での局所コホモロジーを同一視できる。

:::proof[局所双対性の証明]
$R$ の正則パラメータ系を $x_1,\ldots,x_n$、$J_q=(x_1^q,\ldots,x_n^q)$ とする。$J_q$ は極大イデアルのべきと互いに深いので、先に示した直極限で $\operatorname{Ext}_R^j(R/J_q,R)$ を用いられる。各べきの列は Cohen–Macaulay 環のパラメータ系なので正則列である。Koszul 分解の双対は、最高外積に対する外積のペアリングで元の複体を次数反転したものと同型になる。従ってコホモロジーは最上次数 $n$ の $R/J_q$ に限られる。遷移写像は $\prod_i x_i$ 倍である。実際、$R/J_{q+1}\to R/J_q$ を持ち上げる Koszul 写像では、各外積の該当基底因子に $x_i$ を掛ける。

この遷移は単射である。$J_{q+1}$ を $x_1,\ldots,x_n$ の順に商イデアルで割る操作を考える。ある $x_i$ で割る段階では、他の変数の固定したべきを生成するイデアルで先に割る。その環で $x_i$ は正則である。それらのべきと $x_i$ がパラメータ系なので、任意の順で正則列になるためである。従って $(x_i^{q+1}):x_i=(x_i^q)$ となり、全段で $(J_{q+1}:\prod_i x_i)=J_q$ を得る。

$E=E_R(k)$ とする。$R/J_q$ は正則列による Artin 商なので Gorenstein であり、ソクルは一次元である。$R/J_1=k$ の埋め込みを $E$ に取り、各遷移に沿って入射性で埋め込みを順に延長する。延長は単射である。既に埋め込んだ非零部分は新しい Artin 環加群のソクルと交わり、ソクルが一次元なのでそれを含む。従って核が非零ならソクルを含んでしまう。像は $0:_E J_q=D(R/J_q)$ に入り、両辺の長さが等しいためそこに一致する。$E$ は極大イデアルのべきで消される元の合併なので、直極限は $E$ に等しい。以上により $R$ の Čech 複体のコホモロジーは $E$ が次数 $n$ に集中する形である。

有限自由分解 $F_\bullet\to M$ を $R$ 上で取り、この Čech 複体とテンソルする。Čech の各項が平坦なので一方向のホモロジーは $M$ の Čech 複体を与える。他方から計算すると $F_\bullet\otimes_R E$ が次数 $n$ の行に残る。従って $H^i_{\mathfrak m_R}(M)\simeq H_{n-i}(F_\bullet\otimes_R E)$ である。入射加群 $E$ への Hom は完全で、$\End_R(E)=\widehat R=R$ である。このため双対を取ると右辺は $\operatorname{Ext}_R^{n-i}(M,R)$ となる。$A=R/J$ に対し $0:_E J$ は $A$ 上入射であり $k$ を本質的に含む。従って $E_A(k)\simeq 0:_E J$ であり、$A$ 加群上の双対は $D_A$ と同じになる。Čech 複体の環を換える比較により、求める $A$ 上の同型を得る。
:::

例えば $R=k[[t]]$ では $H^1_{(t)}(R)=k((t))/k[[t]]$ となり、負指数の単項式 $t^{-1},t^{-2},\ldots$ がベクトル空間の基底である。$R/(t^q)$ の元 $t^i$ は極限で $t^{i-q}$ に対応し、$t$ を掛けて指数が $0$ になった項は零になる。双対は $R$ 自身である。これは局所双対性の右辺が $\Hom_R(R,R)$ であることとも一致する。

### 標準加群から Gorenstein 性を判定する

:::proposition[完備 Cohen–Macaulay 環の標準加群]{#prop-comm-canonical-gorenstein}
完備 Cohen–Macaulay 局所環 $A$ で、$A$ が Gorenstein であることと $\omega_A\simeq A$ は同値である。
:::

:::proof
Cohen の構造定理による完備正則局所環 $R$ からの全射を取り、$c=\dim R-\dim A$ とする。$A$ の任意のパラメータ系は正則列であり、そのべきも正則列である。従って Koszul 分解から極限を取ると $H^i_{\mathfrak m}(A)$ は $\dim A$ 以外で零である。局所双対により $\operatorname{Ext}_R^j(A,R)$ も $c$ 以外では零で、次数 $c$ の値が $\omega_A$ である。従って :ref[prop-comm-rees-ext] と同様に、$R$ の入射分解に $\Hom_R(A,-)$ を適用し、自由分解との二重複体で比較すると $\operatorname{Ext}_A^i(M,\omega_A)\simeq\operatorname{Ext}_R^{i+c}(M,R)$ となる。右辺は十分大きい $i$ で零なので、$\omega_A\simeq A$ なら有限自己入射次元をもつ。

逆向きには $A$ が Gorenstein であるとする。次元零では $\omega_A=D_A(A)=E_A(k)$ であり、Artin 環の判定により $A\simeq\omega_A$ である。正の次元ではパラメータ系の先頭元 $x$ を取る。Cohen–Macaulay 性により $x$ は正則である。$A/xA$ は Cohen–Macaulay かつ Gorenstein であり、完備である。上の消滅を $A$ と $A/xA$ に適用すると $R$ 上の長完全列から $0\to\omega_A\xrightarrow{x}\omega_A\to\omega_{A/xA}\to0$ を得る。帰納法で $\omega_{A/xA}\simeq A/xA$ である。生成元を持ち上げると中山の補題によって $A\to\omega_A$ が全射となる。$x$ が $\omega_A$ 上で正則であるため、この全射の核を $x$ で割った加群は $A/xA$ に単射で入る。剰余の写像は同型なので核は中山の補題で零である。
:::

### 二変数の局所コホモロジーを分数で読む

:::example[二つとも負の指数だけが残る]
$R=k[[x,y]]$、$\mathfrak m=(x,y)$ とする。Čech 複体の最上項の余核は $R_{xy}/(R_x+R_y)$ である。分母を $x^a y^b$ として分子の級数を書くと、$x$ の指数が $a$ 以上の部分は $R_y$ に入り、$y$ の指数が $b$ 以上の残りは $R_x$ に入る。従って $x^{-i}y^{-j}$（$i,j\ge1$）の有限和が残る。独立性は、それらの単項式の係数が $R_x$ と $R_y$ のどちらに対しても零になることから従う。例えば $x^{-2}y^{-3}$ は $x^2$ または $y^3$ を掛けると消え、$xy^2$ を掛けるとソクルの $x^{-1}y^{-1}$ になる。局所双対性ではこの加群の双対が $\Hom_R(R,R)=R$ となる。
:::

:::example[有限長加群には零次だけが残る]
$A=k[[x,y]]$、$M=A/(x,y)^2$ とする。$x,y$ は $M$ 上冪零なので $M_x=M_y=M_{xy}=0$ である。従って Čech 複体は零次に $M$ を置いただけになり、$H^0_{\mathfrak m}(M)=M$、$H^i_{\mathfrak m}(M)=0$（$i>0$）である。局所双対性で $D_A(M)\simeq\operatorname{Ext}_A^2(M,A)$ となり、これは前節で計算した余核と同じである。有限長 $M$ の双対は長さを保つので余核の長さは $3$ であるが、最小生成元数は $2$ である。
:::

### Artin 商を逆系から作る

体 $k$ 上の $R=k[[x_1,\ldots,x_r]]$ を取る。:ref[thm-comm-local-duality] の証明と Čech 複体の計算により、$k$ の入射包絡は $e_{\boldsymbol a}$（$\boldsymbol a\in\mathbb N^r$）を $k$ 基底とする加群であり、$x_i e_{\boldsymbol a}=e_{\boldsymbol a-\boldsymbol e_i}$ と作用する。負の成分が出たときは零とする。この基底は Čech 複体における $x_1^{-a_1-1}\cdots x_r^{-a_r-1}$ の類である。例えば一変数では $e_a=t^{-a-1}$ の類で、二変数では両方の負べきが残る。任意の有限和 $F$ は極大イデアルのあるべきで消えるから、べき級数による作用も有限和だけで定義できる。

$A=R/I$ が Artin 環なら、$D_R(A)=0:_E I$ である。実際、準同型は $1$ の像で決まり、その像が $I$ に消されることが必要十分である。$A$ が Gorenstein なら $D_R(A)$ は $A$ と同型なので、一つの $F$ で生成できる。逆に $0:_E I=RF$ なら、その加群の零化イデアルは $I$ である。これには双対による非零性の検出を使う。もし $a\notin I$ なら $A$ 上の $a$ 倍の像は非零であり、その双対も非零なので、$a$ 倍は $D_R(A)$ を消さない。従って $A\to RF$、$1\mapsto F$ は同型であり、$D_R(A)$ は $A$ 上の入射包絡だから $A$ は Gorenstein である。

また、任意の非零 $F\in E$ から $I=\operatorname{Ann}_R(F)$ と置けば、このイデアルは極大イデアルのべきを含む。従って $A=R/I$ は Artin 局所環である。$RF\subset D_R(A)$ であり、$RF\simeq A$ だから両方の長さが等しい。従って $RF=D_R(A)$ となり、$A$ は Gorenstein である。これにより Gorenstein な Artin 商は、一つの負べきの組合せから構成できる。

:::example[一つの基底元による商]
二変数で $F=e_{(2,1)}$ とする。$x^iy^jF$ は $0\le i\le2$、$0\le j\le1$ のとき非零で、異なる組は異なる基底元を与える。それ以外は零である。従って零化イデアルは $(x^3,y^2)$ であり、商の長さは $6$ である。$x^2yF=e_{(0,0)}$ なので、商のソクルは $x^2y$ の類で張られる。単項式の指数を一つずつ減らす操作が、商の乗法の逆向きの計算を与える。
:::

:::example[二つの基底元を足したときの非単項式関係]
$F=e_{(2,0)}+e_{(0,3)}$ とする。$xyF=0$、$x^2F=y^3F=e_{(0,0)}$ なので $J=(xy,x^2-y^3)$ が零化イデアルに入る。$x(x^2-y^3)$ と $xy$ から $x^3$、$y(x^2-y^3)$ と $xy$ から $y^4$ を得る。従って $R/J$ は $1,x,y,y^2,y^3$ で張られる。これらを $F$ に作用させた像は順に $F,e_{(1,0)},e_{(0,2)},e_{(0,1)},e_{(0,0)}$ であり、独立である。よって $I=J$、商の長さは $5$ である。二つの純粋冪が同じソクル元になるため、負べきが二項あってもソクルは一次元になる。
:::

### 局所コホモロジーから深さと次元を読む

:::proposition[局所コホモロジーの最初と最後の非零次数]{#prop-comm-local-cohomology-depth-dimension}
Noether 局所環 $(A,\mathfrak m,k)$ と有限非零加群 $M$ に対し、$H^i_{\mathfrak m}(M)\ne0$ となる最小の $i$ は $\operatorname{depth}_A M$、最大の $i$ は $\dim_A M$ である。
:::

:::proof
完備化は平坦であるため Čech 複体のコホモロジーとテンソルを交換する。完備化は忠実平坦なので、各次数での非零性は変わらない。加群の深さは剰余体との Ext で、次元は Hilbert–Samuel 多項式でそれぞれ保たれる。Ext に使う各自由分解項は有限であり、テンソル後の剰余体や $\mathfrak m$ のべきで割った長さも元と一致するためである。従って $A$ が完備と仮定してよい。

Cohen の構造定理で完備正則局所環 $R$ から $A$ への全射を取り、$n=\dim R$ とする。$M$ の $R$ 上の最小自由分解の長さは $p=n-\operatorname{depth}_R M$ である。これは Auslander–Buchsbaum の公式であり、$\operatorname{depth}_R M=\operatorname{depth}_A M$ は $M$ へのスカラー作用が $A$ を経由するため、正則列が両方で同じ条件となることから従う。$p$ より大きい Ext は消え、$\operatorname{Ext}_R^p(M,R)$ は非零である。最後の Hom の項を剰余体で割ると、直前の微分は最小性により零となり、非零の自由項の剰余が残るからである。局所双対性によって $H^i_{\mathfrak m}(M)$ の最初の非零次数は $n-p=\operatorname{depth}_A M$ となる。

最大次数の上界については、$A/\operatorname{Ann}_A M$ で長さ $d=\dim_A M$ のパラメータ系を取る。この環の極大イデアルとパラメータイデアルの根基は同じなので、両者のねじれ関手は等しい。実際、各イデアルの十分大きいべきが他方に入るためである。局所コホモロジーはこの商環で計算しても変わらず、$d$ 元の Čech 複体には次数 $d$ を超える項がない。従って $i>d$ では消える。

非零性のためには、$R$ 上での $M$ の支持で $\dim R/P=d$ を達成する極小素 $P$ を取る。$R$ は Cohen–Macaulay であり、$P$ 内の極大な $R$ 正則列の長さは
$$
\inf_{Q\supset P}\operatorname{depth}R_Q
=\inf_{Q\supset P}\dim R_Q=\operatorname{ht}P
$$
である。最初の等号は grade の局所的特徴付けと正則性の局所化を使ったものである。この正則列を $R_P$ に局所化すると長さが次元に等しいので、列の生成するイデアルはそこで極大を根基とする。従って $P$ は元の列が生成するイデアルの上の極小素である。その列で割った $R$ は次元 $n-\operatorname{ht}P$ の Cohen–Macaulay 局所環で、非混合性から各極小素の商がその次元をもつ。よって $\dim R/P=n-\operatorname{ht}P=d$ である。

$M_P$ は支持の極小点で非零有限長なので深さは零である。$R_P$ は次元 $n-d$ の正則局所環であり、Auslander–Buchsbaum の公式から $\operatorname{pd}_{R_P} M_P=n-d$ である。その最小分解の最後の双対項は非零の Ext を与える。Ext は有限自由分解で計算するため局所化と交換し、従って $\operatorname{Ext}_R^{n-d}(M,R)\ne0$ である。局所双対性により $H^d_{\mathfrak m}(M)\ne0$ となる。
:::

### Cohen–Macaulay 環の内部で局所双対を表す

$A$ を完備 Cohen–Macaulay 局所環、$d=\dim A$、$\omega_A$ を標準加群とする。任意の有限 $A$ 加群 $M$ について
$$
H^i_{\mathfrak m}(M)\simeq D_A\bigl(\operatorname{Ext}_A^{d-i}(M,\omega_A)\bigr)
$$
が成り立つ。実際 $A=R/J$ となる完備正則局所環 $R$ に対し $c=\dim R-d$ と置く。$A$ の局所コホモロジーは次数 $d$ にしかないので、:ref[thm-comm-local-duality] により $\operatorname{Ext}_R^j(A,R)$ は次数 $c$ の $\omega_A$ に集中する。この集中性から、$R$ の入射分解に $\Hom_R(A,-)$ を施した複体は、各項が $A$ 入射でコホモロジーが次数 $c$ に集中する。$M$ の自由分解とこれを組み合わせ、一方からは $\operatorname{Ext}_A^{d-i}(M,\omega_A)$ を、他方からは $\operatorname{Ext}_R^{c+d-i}(M,R)$ を得る。従ってこの二つが同型である。局所双対性によれば後者は $D_A(H^i_{\mathfrak m}(M))$ である。

有限 $A$ 加群の Ext は有限生成なので、その Matlis 双対は Artin である。ここで $H^i_{\mathfrak m}(M)$ はその双対と同型になる。評価写像が単射であり、$D_A(H^i_{\mathfrak m}(M))$ が有限なので、二重双対側は Artin である。従って $H^i_{\mathfrak m}(M)$ も Artin となり、Matlis 双対性から評価は同型である。これで表示を得る。例えば $M=k$ では $H^0_{\mathfrak m}(k)=k$、高次は零なので、$\operatorname{Ext}^d_A(k,\omega_A)=k$ となる。この式は $\omega_A$ を第二変数にしたものであり、$A$ を第二変数にする場合は Gorenstein 性によって両者を同一視できる。

$M=A$ とすると $H^d_{\mathfrak m}(A)\simeq D_A(\omega_A)$ であり、$M$ を有限直和や商に替えても写像と両立する。同型の選択に依存するのは標準加群を具体的な表示で書くときで、局所コホモロジーから得られる標準加群の同型類は固定される。

### 支持イデアルの根基と Čech 複体

二つの有限生成イデアル $I,J$ の根基が等しいなら、ある $a,b$ で $I^a\subset J$、$J^b\subset I$ が成り立つ。各生成元のあるべきが他方に入り、十分高い全次数の積にはそのいずれかが因子として現れるためである。この両包含から $\Gamma_I(M)=\Gamma_J(M)$ が任意の $M$ で成立する。入射分解に適用すれば全ての局所コホモロジーも同型になる。例えば $A=k[[x,y]]$ の $(x^2,y^3)$ は $(x,y)$ と根基が等しい。Čech 複体でも $M_{x^2}=M_x$、$M_{y^3}=M_y$、$M_{x^2y^3}=M_{xy}$ なので同じ分数が残る。関係式としてのイデアルの違いは有限長商に現れるが、支持のみに依存するこの計算では現れない。

## 演習問題

::::exercise{level=basic}
$\Q/\Z$ が $\Z$ 上入射であることを $3\Z$ の場合に確認せよ。

:::solution
$f(3)=a+\Z$ と書く。$1$ を $a/3+\Z$ に送る写像は $3$ を $a+\Z$ に送るので延長である。
:::
::::

::::exercise{level=basic}
$k[t]/(t^5)$ のソクルを求めよ。

:::solution
$t$ 倍で消える元は $t^4$ の倍数に限られる。従って $kt^4$ である。
:::
::::

::::exercise{level=basic}
$k[x,y]/(x,y)^2$ のソクルの次元を求めよ。

:::solution
定数項が非零なら $x$ を掛けた積は非零である。一方 $x,y$ の線形結合は極大イデアルに消される。従って次元は $2$ である。
:::
::::

::::exercise{level=basic}
$\dim A=3$ の正則局所環の自己入射次元を求めよ。

:::solution
大域次元で上界 $3$ を得る。三元 Koszul 分解を双対した最上次数の余核が $k$ なので下界も $3$ である。
:::
::::

::::exercise{level=advanced}
Noether 局所環 $(A,\mathfrak m,k)$ と $q\ge1$ に対し、$\ell(D(A/\mathfrak m^q))$ を求めよ。

:::solution
$A/\mathfrak m^q$ の組成列を双対すると順序を逆にした列ができ、各単純因子の双対は $k$ である。よって長さは $\ell(A/\mathfrak m^q)$ である。
:::
::::

::::exercise{level=advanced}
$A=k[[t]]$ で $H_{(t)}^1(A)$ の $t^2$ に消される部分の $k$ 基底を求めよ。

:::solution
$k((t))/k[[t]]$ の負指数表示を使うと、$t^2$ を掛けて負指数がなくなるのは $t^{-1},t^{-2}$ の線形結合だけである。
:::
::::

::::exercise{level=advanced}
$k[[x,y]]/(xy)$ の標準加群を求めよ。

:::solution
母環上 $xy$ 倍による長さ一の自由分解があり、その双対の一次余核は母環の $(xy)$ による商である。従って標準加群は自身と同型である。
:::
::::

::::exercise{level=advanced}
$A$ を完備 Noether 局所環とし、$r$ を非負整数とする。$D(A^r)$ と $D(D(A^r))$ を求めよ。

:::solution
有限直和の Hom は $E^r$ であり、再び双対すると $\End_A(E)^r=\widehat A^r=A^r$ である。これは評価と一致する。
:::
::::

::::exercise{level=advanced}
$A=k[x,y]/(x^3,xy,y^2)$ の $\operatorname{Soc}(A)$ と Gorenstein 性を判定せよ。

:::solution
標準単項式は $1,x,x^2,y$ である。$a+bx+cx^2+dy$ に $x$ を掛けると $ax+bx^2$、$y$ を掛けると $ay$ になるので、両方で零となる条件は $a=b=0$ である。従ってソクルは $kx^2\oplus ky$、次元は $2$ である。Artin 局所環の判定より $A$ は Gorenstein でない。
:::
::::

## 確認テスト

::::quiz
Baer の判定で延長する写像の定義域を選べ。

- [ ] 任意の単元の集合。
- [x] 任意のイデアル。
- [ ] 剰余体だけ。

:::explanation
極大な部分加群への延長から一元を加えるとき、その元が既存の部分加群に入る係数がイデアルをなす。そのイデアルからの写像を延長できれば全体へ延長できる。
:::
::::

::::quiz{answer="1" tol="1e-6"}
$k[x,y]/(x^2,y^3)$ のソクルの次元を求めよ。

:::explanation
標準単項式のうち $x$ と $y$ の両方を掛けて零となるのは $xy^2$ だけなので一次元である。
:::
::::

::::quiz
完備局所環で有限生成加群の Matlis 双対は何か。

- [x] Artin 加群。
- [ ] 常に自由加群。
- [ ] 常に有限長加群。

:::explanation
有限自由全射の双対により有限個の $E$ の直和に埋め込める。$E$ は Artin なので部分加群も Artin である。
:::
::::

::::quiz{answer="4" tol="1e-6"}
$A=k[[t]]$ で $0:_{H^1_{(t)}(A)}t^4$ の $k$ 次元を求めよ。

:::explanation
$H^1_{(t)}(A)=A_t/A$ で、$t^4$ に消される基底は $t^{-1},t^{-2},t^{-3},t^{-4}$ である。
:::
::::

::::quiz{answer="6" tol="1e-6"}
$R=k[[x,y]]$ とし、$H^2_{(x,y)}(R)$ のうち $x^3$ と $y^2$ の双方で消される部分の $k$ 次元を答えよ。

:::explanation
負の単項式 $x^{-i}y^{-j}$（$i,j\ge1$）が基底である。$x^3$ 倍で消えるには $1\le i\le3$、$y^2$ 倍で消えるには $1\le j\le2$ が必要十分である。消えない単項式の像は互いに異なるので、線形結合でも追加の打ち消しは起きない。従って次元は $3\cdot2=6$ である。
:::
::::

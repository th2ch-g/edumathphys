---
title: 基数
description: 始順序数による基数の表示、無限基数の算術、共終数、連続体仮説と到達不能基数を扱う。
notation:
  - symbol: '|\!X\!|'
    meaning: 集合 X の基数
  - symbol: '\aleph_\alpha'
    meaning: α 番目の無限基数
  - symbol: '\beth_\alpha'
    meaning: α 段まで冪集合の濃度を反復した値
  - symbol: '\operatorname{cf}(\alpha)'
    meaning: 順序数 α の共終数
---

## この章の目標

- 始順序数で無限基数を表示する。
- 和、積、冪の大小を集合の写像で比較する。
- König の不等式から共終数の制約を導く。
- GCH と到達不能基数の定義を使って基数階層を読む。

## 基数と始順序数

:ref[thm-sets-choice-cycle] により、ZFC では任意の集合を整列できる。各整列に順序型があるので、集合 $X$ と全単射を持つ順序数が存在する。既知の一つの順序型以下から最小のものを取る。それを $|X|$ と書く。自分より小さい順序数と全単射を持たない順序数を :term[始順序数]{en="initial ordinal" yomi="しじゅんじょすう"} と呼び、本章の基数とする。

単射 $X\to Y$ があれば $|X|\le |Y|$ である。逆に順序数の包含を全単射と合成すれば単射が得られる。両方向の単射から全単射を得るには :ref[thm-sets-cantor-bernstein] を使う。無限基数 $\kappa$ より大きい最小の基数を $\kappa^+$ とする。:ref[lem-sets-hartogs] は、任意の集合へ単射を持たない順序数の存在を選択なしに保証する。これにより $\kappa^+$ が存在する。

:::example[順序型と濃度]
$\omega+1$ は最大元を持つので $\omega$ と順序同型でない。一方、最後の点を $0$、各有限点 $n$ を $n+1$ に送ると $\omega$ への全単射を得る。従って $\omega+1$ は始順序数ではなく、基数は $\aleph_0=\omega$ である。
:::

### 初期区間の大小と基数の大小

ある順序数 $\alpha$ の濃度を表す始順序数は、$\alpha+1$ に属する順序数のうち $\alpha$ と全単射を持つ最小元として取れる。この構成で必要なのは $\alpha$ が既に整列していることである。任意の集合 $X$ の濃度を始順序数で表すには、まず $X$ を整列する。

:::example[有限の余りを順序から濃度へ移す]
$\omega+5$ は $\omega$ より大きい順序数である。濃度については末尾の $\omega,\ldots,\omega+4$ を $0,\ldots,4$ に送り、前半の $n$ を $n+5$ に送れば $\omega$ への全単射になる。従ってその基数は $\aleph_0$ である。順序が保たれないのは末尾の五元が像では先頭に置かれるからである。
:::

:::example[集合の大きさから Hartogs 数を得る]
有限集合 $X=\{a,b,c\}$ を取る。順序数 $0,1,2,3$ はそれぞれ $X$ に単射で入る。$4$ は四つの異なる像を要求するため単射がない。この場合の Hartogs 数は $4$ である。一般の集合の Hartogs 数については :ref[lem-sets-hartogs] を参照。一般に $X$ が整列可能なら、Hartogs 数は $|X|$ より大きい次の始順序数である。従って ZFC ではこの対応を任意の集合に適用できる。
:::

### 全射と単射を向きごとに読む

空でない $Y$ に対し単射 $Y\to X$ があれば、$X$ から $Y$ への全射がある。像に属する点には逆写像を使い、像から外れる点には固定した $y_0\in Y$ を返せばよい。逆に全射 $X\to Y$ から単射 $Y\to X$ を取るには、各ファイバーに選択を使える。$X$ が既に整列されていれば各ファイバーの最小元を取れる。従って全射に解を一つずつ選ばせる箇所で、選択の要否を確認する。

### 比較可能性と選択公理

選択公理なしでも、整列可能な二集合は濃度で比較できる。それぞれを順序数へ移せば、二つの順序数は包含で比較できるからである。さらに、任意の二集合の間で少なくとも一方向に単射があるという仮定から、選択公理を逆に導ける。任意の集合 $X$ に Hartogs 数 $h(X)$ を取る。定義により $h(X)$ から $X$ への単射がないため、比較可能性によって $X$ から $h(X)$ への単射がある。その単射で順序数の整列順序を $X$ に戻せるので、任意の集合が整列可能になる。従って :ref[thm-sets-choice-cycle] により選択公理が従う。

ここで Cantor–Bernstein の定理は、両方向に単射が既にあるときに全単射を与える。今の逆向きの証明は、どちらか一方向の単射が存在すること自体を仮定している。例えば $Y$ が整列可能で $X$ が整列可能でないとき、$X\to Y$ の単射は存在しない。逆の単射を考えるときは、$Y$ の大きさに依存する。$Y=\{0\}$ なら $X$ が非空であることだけで単射を取れるが、$Y=h(X)$ ならどちら向きにも単射がない。

:::example[順序数を集合に写した像を整列する]
写像 $f:\alpha\to X$ が全射なら、各 $x\in X$ に対して $f(\xi)=x$ となる最小の $\xi<\alpha$ を割り当てられる。像が異なる点には異なる最小添字が付くため、これで単射 $X\to\alpha$ を得る。従って $X$ は整列可能である。この向きでは、各ファイバーの代表を選ぶのに選択公理を使わずに済む。例えば $\omega$ からの列に繰り返しがあっても、各値の初出の位置だけを残して像を整列できる。
:::

### 無限始順序数に最大元はない

無限順序数 $\alpha$ なら $\omega\subseteq\alpha$ である。$\alpha+1$ から $\alpha$ への単射を、追加点 $\alpha$ は $0$ へ、$n<\omega$ は $n+1$ へ、$\omega\le\xi<\alpha$ は $\xi$ へ送る写像で作れる。各範囲が互いに素であり、どの点にも逆像があるため全単射である。従って後続順序数 $\alpha+1$ は無限始順序数にはなれない。すべての無限基数は順序数としては極限である。

この観察は、後続基数 $\kappa^+$ と順序数の後者 $\kappa+1$ を区別するのにも使える。例えば $\omega_1$ は $\omega$ の後続基数でありながら、自身より小さい順序数の最大元を持たない。無限基数 $\kappa$ について $\kappa+n$ は任意の有限 $n$ で濃度 $\kappa$ のままであるが、$\kappa^+$ はそれらすべてより大きい。$\kappa^+$ の「後続」は、全順序数から始順序数だけに絞ったときの次の位置を表す。

## アレフ階層

:::definition[アレフとベート]{#def-logic-aleph-beth}
$\aleph_0=\beth_0=\omega$ とする。後続段で $\aleph_{\alpha+1}=(\aleph_\alpha)^+$、$\beth_{\alpha+1}=2^{\beth_\alpha}$ とし、極限段では先行する全値の上限とする。ここで $2^\kappa=|\mathcal P(\kappa)|$ である。
:::

極限で無限基数の増加列の上限 $\lambda$ を取ると、$\lambda$ も基数である。もし $\beta<\lambda$ と全単射があれば、列の値の一つ $\kappa$ が $\beta$ より大きい。その $\kappa$ は $\lambda$ に含まれるから $\beta$ へ単射を持つことになり、始順序数であることに反する。従って再帰は基数値を保つ。アレフは後続段で間の基数を飛ばさず、極限で先行値の上限を取る。そのため全無限基数を列挙する。全ての無限基数が現れることは下で、ある基数以上となる最初の添字を取って示す。

:::example[二つの列の差]
$\beth_1$ は自然数集合の部分集合の数であり、Cantor の定理 :ref[thm-sets-cantor] により $\aleph_0$ より大きい。従って $\aleph_1\le\beth_1$ である。等号が成り立つかは後の節で扱う。
:::

### 後続段と極限段を扱う

アレフ列の後続段は一つ上の始順序数である。極限段 $\lambda$ で $\kappa=\sup_{\alpha<\lambda}\aleph_\alpha$ とすれば、$\kappa$ も始順序数である。実際、$\beta<\kappa$ なら $\beta<\aleph_\alpha$ となる $\alpha<\lambda$ があり、もし $\beta$ と $\kappa$ が全単射なら $\aleph_\alpha$ を $\beta$ に単射で移せてしまう。これは $\aleph_\alpha$ の最小性に反する。全無限基数を尽くすのは、後続段で間を飛ばさず、極限段では先行段の上限を直ちに入れるからである。形式的には、ある基数を越える最初の添字を取ると、後続でも極限でもその基数が途中に残ることはない。

:::example[最初の極限添字]
$\aleph_\omega$ は全有限添字の $\aleph_n$ より大きい。$\aleph_n$ と $\aleph_\omega$ の間には $\aleph_{n+1}$ が入るため、どの有限 $n$ の後続でもない。一方 $\aleph_{\omega+1}$ は $\aleph_\omega$ の直後の基数である。添字 $\omega+1$ が後続であることと、$\aleph_{\omega+1}$ が順序数として後続であることは異なり、後者は無限始順序数なので極限順序数である。
:::

:::example[ベート列を基数の添字に換える]
仮に $2^{\aleph_0}=\aleph_2$、$2^{\aleph_2}=\aleph_5$ なら、$\beth_1=\aleph_2$、$\beth_2=\aleph_5$ である。この情報から $2^{\aleph_1}$ が必ず $\aleph_3$ などと決まるわけではない。単調性から直接言えるのは $\aleph_2\le2^{\aleph_1}\le\aleph_5$ であり、一般の冪演算には後続基数の定義とは違う制約が働く。
:::

### 正規な基数列と固定点

添字 $\alpha$ に対する $\aleph_\alpha$ は狭義増加で、極限添字ではそれ以前の上限を取る。従って任意の順序数 $\gamma$ を越える $\aleph$ の固定点がある。具体的に $\alpha_0>\gamma$ を選び、$\alpha_{n+1}=\aleph_{\alpha_n}$ とする。$\aleph_\alpha\ge\alpha$ は超限帰納法で従うので、この列は弱く増加する。途中で等しくなれば固定点であり、増加し続ければ $\lambda=\sup_n\alpha_n$ について $\aleph_\lambda=\sup_n\aleph_{\alpha_n}=\lambda$ となる。$\beth$ 列も狭義増加で極限連続なので同じ構成を使える。固定点が始順序数であることは $\lambda$ が $\aleph_\lambda$ または $\beth_\lambda$ に等しいことから従う。

### アレフの添字を使って全無限基数を拾う

$\alpha\mapsto\aleph_\alpha$ は狭義増加であり、各極限で先行値の上限を取る。狭義増加な順序数列の第 $\alpha$ 項は少なくとも $\alpha$ なので、どの基数 $\kappa$ に対しても、例えば $\aleph_{\kappa+1}$ は $\kappa$ より大きい。そこで $\kappa\le\aleph_\alpha$ となる最小の $\alpha$ を取れる。$\alpha=0$ なら $\kappa=\aleph_0$ である。$\alpha=\beta+1$ なら $\aleph_\beta<\kappa\le(\aleph_\beta)^+$ であり、$\kappa$ は基数だから $\kappa=(\aleph_\beta)^+$ となる。$\alpha$ が極限なら全ての $\beta<\alpha$ で $\aleph_\beta<\kappa$ なので $\aleph_\alpha=\sup_{\beta<\alpha}\aleph_\beta\le\kappa$ となり、やはり等しい。これは極限と後続を分けることで、添字を取り漏らさないことを示す。

$\beth$ 列も狭義増加かつ連続だが、後続段でその時点の冪集合を取るため、間の基数が残ることがある。帰納法により $\aleph_\alpha\le\beth_\alpha$ が全 $\alpha$ で成り立つ。基底は等号である。後続では $\aleph_\alpha\le\beth_\alpha<2^{\beth_\alpha}$ だから、$2^{\beth_\alpha}$ は $\aleph_\alpha$ の次の基数以上である。極限では点ごとの不等式の上限を取る。各列の値が無限基数だから、後続の不等式を「順序数の直後」に置き換えずに扱える。

### 極限添字の共終列を作る

$\lambda$ が非零極限順序数なら、$\aleph_\lambda$ の下の任意の順序数 $x$ は、ある $\alpha<\lambda$ に対して $x<\aleph_\alpha$ を満たす。これは極限段が和集合であることを元ごとに読んだ主張である。例えば $\lambda=\omega\cdot2$ のとき、$\omega+n$ が $\lambda$ に共終なので、$\aleph_{\omega+n}$ は $\aleph_{\omega\cdot2}$ に共終である。$\aleph_{\omega\cdot2}$ の濃度は $\aleph_\omega$ より大きいが、どちらにも可算共終列がある。

同じ理由で、$\lambda$ に共終な任意の添字集合 $A$ に沿って $\aleph_\lambda=\sup_{\alpha\in A}\aleph_\alpha$ となる。もし $A$ が有界なら、ある $\gamma<\lambda$ が全ての $\alpha\in A$ を上から抑える。その場合に得る値は高々 $\aleph_\gamma<\aleph_\lambda$ である。従って極限の階層値を少数の項で近づけるには、添字側の共終性を確認する必要がある。$\beth$ 列でも同じ二方向が成り立つ。

## 基数の和と積

:term[基数の和]{en="cardinal sum" yomi="きすうのわ"} は添字で区別した非交和の濃度であり、積は直積の濃度である。二項の場合、和は $(\{0\}\times X)\cup(\{1\}\times Y)$ で表す。因子を交換する全単射で可換律が得られ、直積の括り直しで結合律が得られる。直積 $X\times(Y\sqcup Z)$ では第二座標の添字で二つの直積へ振り分けられるため、分配律が成立する。

無限基数 $\kappa$ に対して $\kappa^2=\kappa$ という基数の積の結果は :ref[thm-sets-cardinal-square] による。そこでは対を最大座標の順に並べ、各初期区間の基数を小さく抑える。従って $1\le\lambda\le\kappa$ なら $\kappa\lambda=\kappa$ である。実際、一点を固定して $\kappa$ を直積へ入れ、直積を $\kappa\times\kappa$ に入れる。$\kappa+\lambda=\kappa$ も、タグの集合を $\kappa$ に入れて同じ比較を使えば従う。

:::example[可算集合の列を合わせる]
各 $A_n$ が高々可算なら、選択により各 $A_n$ から $\N$ への単射を一つずつ取る。非交和は $\N\times\N$ に入る。普通の和集合には非交和から全射があるので、高々可算である。特に可算個の有限集合の和も高々可算である。
:::

### 添字集合も無限の場合

集合族 $(X_i)_{i\in I}$ の非交和の濃度は、$|I|$ と各 $|X_i|$ の両方で評価する。$\kappa$ が無限基数で $|I|\le\kappa$ かつ $|X_i|\le\kappa$ なら、各 $X_i$ の $\kappa$ への単射をまとめて選び、非交和を $I\times\kappa$ へ単射で移す。:ref[thm-sets-cardinal-square] で上界 $\kappa$ を得る。別々の座標を選ぶ際には選択公理を使う。

:::example[濃度が増加する可算和]
$X_n$ の濃度を $\aleph_n$ とする。非交和の濃度は各 $\aleph_n$ 以上であり、従って $\aleph_\omega$ 以上である。添字集合は可算で各濃度は $\aleph_\omega$ 以下なので、上の評価で高々 $\aleph_\omega$ となる。最大の項がなくても上限が和の濃度になる例である。
:::

:::example[有限点の集合全体を数える]
濃度 $\kappa$ の無限集合 $X$ の有限部分集合は、各 $n$ の順序付き $n$ 組の像として数えられる。例えば三元部分集合なら、相異なる三元を選んだ組から順序を忘れて得る。固定 $n\ge1$ の組の全体は濃度 $\kappa$、$n=0$ では一元である。可算和で上界 $\kappa$、単元集合の族で下界 $\kappa$ が得られる。この計算は有限条件の全体の濃度を求めるときにも使う。
:::

### 可算積と有限積を分ける

無限基数 $\kappa$ に対し、正の有限 $n$ なら $\kappa^n=\kappa$ である。平方の結果を繰り返し適用すればよい。一方 $\kappa^\omega$ は有限列全体の濃度とは限らない。例えば $\kappa=\aleph_0$ のとき、可算列には各二値列を含むため Cantor の定理で非可算である。有限座標の情報をまとめるだけなら $\bigcup_n\kappa^n$ は濃度 $\kappa$ に収まる。従って「有限長が任意」と「一つの無限長」の差が濃度にも現れる。

### 無限の和を添字と上限に分ける

空でない集合だけからなる族 $(X_i)_{i\in I}$ を考える。$\theta=\max(|I|,\sup_{i\in I}|X_i|)$ が無限なら、
$$
\left|\bigsqcup_{i\in I}X_i\right|=\theta
$$
である。各 $X_i$ は非交和に入り、さらに一元ずつ選べば $I$ も非交和へ入るため、濃度は $\theta$ 以上である。逆に各 $X_i$ の $\theta$ への単射を選び、タグと組み合わせることで $I\times\theta$ への単射を作る。この積は $\theta^2=\theta$ 以下である。空の $X_i$ が混ざる場合には、まず非空の座標だけに添字集合を制限する。全てが空なら非交和は空である。式の右辺が有限の場合は、有限の項を足した総数を直接計算する必要がある。

通常の和 $\bigcup_i X_i$ は、タグを忘れることで非交和の全射像になる。従って上界は同様に成り立つが、添字集合の濃度が下界になるとは限らない。例えば $\omega_1$ 個の座標でどれも同じ単元集合を使えば、非交和の濃度は $\aleph_1$、通常の和の濃度は $1$ である。無限集合 $X$ の上で互いに素な非空部分集合が何個あるかを評価するときは、各部分集合から選んだ代表が衝突しないので、添字集合の濃度は $|X|$ 以下となる。

:::example[項の最大がない場合の積と和]
各 $n<\omega$ に濃度 $\aleph_n$ の集合 $X_n$ を取る。非交和の濃度は $\sup_n\aleph_n=\aleph_\omega$ である。直積では各座標から選ぶ関数が一元となる。特にこの直積の元は、非交和のいずれか一つのブロックの元だけでは指定できない。次節では対角法で $\prod_n\aleph_{n+1}>\aleph_\omega$ を導く。有限直積で使う最大値だけでは、この濃度の差を説明できない。
:::

### 有限台の関数を数える

集合 $I$ と、基点 $0$ を含む集合 $A$ を取る。関数 $f:I\to A$ の台を $\{i:f(i)\ne0\}$ とする。$I$ が無限、$|A|\ge2$ であれば、有限台の関数の全体の濃度は $\max(|I|,|A|)$ である。上から抑えるため、台の点とそこでの値を有限列に記録する。$|I\times A|=\max(|I|,|A|)$ であり、この集合の有限列全体も同じ濃度である。逆に非零値を一つ固定して置き場所だけを変えれば $|I|$ 個が得られ、一つの置き場所を固定して値を変えれば $|A|$ 個が得られる。二つの下界を合わせると等号になる。有限部分関数の場合も、値 $0$ の指定と未定義を別記号で区別すれば同様に数えられる。

:::example[有限部分関数と全関数]
$I=\omega_1,A=\{0,1\}$ とする。有限部分関数の濃度は $\aleph_1$ である。各有限条件は $I\times A$ の有限部分集合だから高々 $\aleph_1$ 個であり、各座標で値 $0$ を指定した一座標条件が $\aleph_1$ 個ある。これに対し全二値関数の濃度は $2^{\aleph_1}>\aleph_1$ である。部分関数すべてを $\omega_1$ に延長する規則を固定しても、その像はまだ $\aleph_1$ 個以下なので全二値関数を尽くさない。
:::

## 基数の冪と König の不等式

$\kappa^\lambda$ は $\lambda$ から $\kappa$ への全写像の数である。$2^\lambda$ を部分集合の濃度と同一視するには指示関数を使う。$(\kappa^\lambda)^\mu=\kappa^{\lambda\mu}$ は二引数の関数を、各第一引数から関数を返す写像として読む全単射による。また $\kappa^{\lambda+\mu}=\kappa^\lambda\kappa^\mu$ は非交和上の写像を二つの制限に分けることで従う。

:::theorem[König の不等式]{#thm-logic-konig}
非空の添字集合 $I$ について各 $\kappa_i<\lambda_i$ なら、$\sum_{i\in I}\kappa_i<\prod_{i\in I}\lambda_i$ である。
:::

:::proof
非交和の $i$ 番目の部分 $X_i$ の大きさを $\kappa_i$、直積の $i$ 座標 $Y_i$ の大きさを $\lambda_i$ とする。和から積への任意の写像 $f$ について、$x\in X_i$ が第 $i$ 座標に取る値は最大 $\kappa_i$ 種類なので、$Y_i$ のある値 $y_i$ が使われない。選択で $\vec y$ を作ると、各 $x\in X_i$ の像は少なくとも第 $i$ 座標で異なる。従って全射はない。全基数が比較可能であり、非空な定義域から小さい非空集合への全射は構成できるため、不等式が従う。非交和が空の場合も直積は非空なので厳密である。
:::

:::example[対角の位置]
有限の二座標で $X_0,X_1$ が一元ずつ、$Y_0,Y_1$ が二元ずつなら、和は $2$、積は $4$ である。対角構成は、一つ目の入力が使わない第一座標と二つ目が使わない第二座標を組にする。Cantor の不等式では各 $X_i$ を一元、各 $Y_i$ を二元に取る。
:::

### 指数を小さくすることと値域を小さくすること

$2\le\lambda\le 2^\kappa$ で $\kappa$ が無限なら、$\lambda^\kappa=2^\kappa$ である。値域に二点があることから $2^\kappa\le\lambda^\kappa$ である。逆は $2^\kappa$ への単射を各座標に適用して $\lambda^\kappa\le(2^\kappa)^\kappa=2^{\kappa\cdot\kappa}=2^\kappa$ と計算する。この方法で、例えば $\aleph_0^{\aleph_0}=2^{\aleph_0}$ と、各有限 $m\ge2$ について $m^{\aleph_0}=2^{\aleph_0}$ を同時に得る。

:::example[関数のグラフによる上界]
$\N$ から $\N$ への各関数を $\N\times\N$ の部分集合であるグラフに送ると単射である。$\N\times\N$ が可算なので、関数は高々 $2^{\aleph_0}$ 個である。各部分集合の指示関数は自然数値関数でもあるため逆向きの下界が成り立つ。上の指数法則による計算を、この例では直接の写像で確認できる。
:::

### 指数の書き換えに使う全単射

$C^{A\times B}$ と $(C^B)^A$ の間で、$h(a,b)$ を $a\mapsto(b\mapsto h(a,b))$ に送れば全単射になる。同様に $C^{A\sqcup B}$ は二つの定義域への制限から $C^A\times C^B$ と全単射である。等式 $(\kappa^\lambda)^\mu=\kappa^{\lambda\mu}$ と $\kappa^{\lambda+\mu}=\kappa^\lambda\kappa^\mu$ はこの対応から従う。

:::example[二つの可算指数を比べる]
$\{0,1\}^\omega$ は $\omega^\omega$ に包含される。一方、$\omega^\omega$ の各元をそのグラフに送ると $\mathcal P(\omega\times\omega)$ に単射になる。$\omega\times\omega$ は可算なので両側の濃度は $2^{\aleph_0}$ となる。例えば列 $f(n)=2n$ のグラフは $(0,0),(1,2),(2,4),\ldots$ の集合であり、異なる列はある座標の対が違うため異なる集合に移る。
:::

:::example[König で特異基数の冪が真に増える]
$\kappa=\aleph_\omega$ とすると $\sum_n\aleph_n=\kappa$ である。各 $\aleph_n<\aleph_{n+1}$ に König の不等式を適用し、$\kappa<\prod_n\aleph_{n+1}\le\kappa^{\aleph_0}$ を得る。後の不等式は各座標を $\kappa$ に含めるだけでよい。無限基数を有限個掛けた積がその最大に等しいことと、可算積まで同じ計算ができることは別である。
:::

### 可算指数が基数を保つ条件

$\kappa=2^\mu$ かつ $\mu$ が無限なら、$0<\lambda\le\mu$ に対して $\kappa^\lambda=2^{\mu\lambda}=2^\mu=\kappa$ である。この等式は $\kappa$ の値が何かという情報より、ある冪であるという情報を使う。一般の特異基数の可算指数は異なり得る。既出の $\aleph_\omega$ の例では $\kappa^\omega>\kappa$ となるため、それが $2^\mu$ になるには共終数の制約を満たす必要がある。

:::example[底が二つの場合の比較]
$\kappa=\aleph_\omega$ に対し、$\kappa^\omega>\kappa$ は König から分かる。一方 $\lambda=2^{\aleph_0}$ には $\lambda^\omega=\lambda$ が成立する。両方とも非可算無限基数であるというだけでは可算指数の挙動は決まらない。最初の例では可算共終な増大列、二番目では可算個の可算集合の直積の添字をどのようにまとめられるかが決め手になる。
:::

### 一般の積を定数底の冪で挟む

非空集合族 $(A_i)_{i\in I}$ について、各 $B_i$ が濃度 $\lambda$ を持つなら、写像の組 $(f_i:A_i\to B_i)_{i\in I}$ を、定義域のタグを保った一つの写像にまとめられる。各 $B_i$ を $\lambda$ と同一視すれば
$$
\prod_{i\in I}\lambda^{|A_i|}
   =\lambda^{\sum_{i\in I}|A_i|}
$$
である。逆向きには一つの写像を各タグの部分へ制限する。直積の底の側でも $\bigl(\prod_i\kappa_i\bigr)^\mu=\prod_i\kappa_i^\mu$ となる。こちらでは $\mu$ の各元から組を返す関数を、各座標を返す関数の組へ組み替える。同じ引数の集合を使用しているため、二方向の組み替えが互いに逆である。

特に $\kappa_i\le\theta$ なら $\prod_i\kappa_i\le\theta^{|I|}$ である。逆向きに、全ての座標で $\kappa_i\ge2$ なら、二値を詰めることで $2^{|I|}\le\prod_i\kappa_i$ となる。この二つから、$I$ が無限で全 $i$ に対して $2\le\kappa_i\le2^{|I|}$ なら $\prod_i\kappa_i=2^{|I|}$ が従う。上界には $(2^{|I|})^{|I|}=2^{|I|\cdot|I|}=2^{|I|}$ を使う。各座標の大きさを個別に知っていなくても、一様な挟み撃ちが可能な場合である。

:::example[有限の底が増える直積]
各 $n<\omega$ で $\kappa_n=n+2$ とする。二値列はこの直積の中に入るので下界は $2^{\aleph_0}$ である。上界には各 $\kappa_n$ を $\omega$ に含め、$\prod_n\kappa_n\le\aleph_0^{\aleph_0}=2^{\aleph_0}$ を使える。従って有限の底に共通の有限上界がなくても、直積の濃度は連続体になる。座標を有限個しか取らない積は有限であり、任意の有限積を合わせても可算である。この可算和は無限積そのものとは異なる。
:::

## 共終数と正則基数・特異基数

:::definition[共終数と正則性]{#def-logic-cofinality}
順序数 $\alpha$ の部分集合 $C$ が :term[共終]{en="cofinal" yomi="きょうしゅう"} とは、各 $\beta<\alpha$ に対し $\beta\le\gamma$ となる $\gamma\in C$ があることをいう。その最小濃度を $\operatorname{cf}(\alpha)$ とする。無限基数 $\kappa$ について、$\operatorname{cf}(\kappa)=\kappa$ なら :term[正則基数]{en="regular cardinal" yomi="せいそくきすう"}、小さければ :term[特異基数]{en="singular cardinal" yomi="とくいきすう"} という。
:::

極限順序数 $\alpha$ に共終な集合を最小濃度 $\mu$ で列挙する。その添字の初期区間の像はまだ有界である。そこで順に、前の値より大きく、列挙した値以上の点を取れば、長さ $\mu$ の狭義増加な共終列を得る。$\mu$ が特異なら、その添字列の共終部分だけを取って $\alpha$ により小さい濃度の共終集合ができる。従って $\mu$ は正則である。

:::proposition[後続基数と共終数]{#prop-logic-successor-regular}
無限基数 $\kappa$ の後続基数 $\kappa^+$ は正則である。
:::

:::proof
濃度 $\le\kappa$ の共終族を仮定する。各添字の値は $\kappa^+$ 未満だから、その値の後者も基数高々 $\kappa$ である。それらの初期区間の和が $\kappa^+$ 全体を覆う。しかし和の基数は $\kappa\cdot\kappa=\kappa$ 以下なので矛盾する。
:::

:::example[最初の特異基数]
$\aleph_\omega=\sup_n\aleph_n$ には可算共終列がある。極限なので有限共終集合は最大元の直後を覆えない。従って $\operatorname{cf}(\aleph_\omega)=\aleph_0$ である。
:::

無限基数 $\kappa$ について $\operatorname{cf}(2^\kappa)>\kappa$ が成り立つ。$\lambda=2^\kappa$ に濃度 $\mu\le\kappa$ の共終列があり、その各値を $\lambda_i<\lambda$ の基数で上から抑えられると仮定する。$\lambda$ が後続基数なら直前の命題で既に矛盾する。極限基数なら列の各項より大きい基数を $\lambda$ 未満で取れる。このとき $\sum_{i<\mu}\lambda_i=\lambda$ である。一方 :ref[thm-logic-konig] を各 $\lambda_i<\lambda$ に適用すると $\lambda<\lambda^\mu\le(2^\kappa)^\kappa=2^\kappa$ となり矛盾する。

### 共終列と初期区間の被覆

非零極限順序数 $\alpha$ の共終数は無限である。有限部分集合なら最大元があり、その次の点を覆えない。逆に後続順序数の共終数は $1$、空順序数の共終数は $0$ である。共終列を組み合わせるときは、列の濃度だけでなく、その像が実際に上まで届くかを確認する。

:::example[順序数の濃度と共終数]
$\omega_1+\omega$ は濃度 $\aleph_1$ だが、末尾の $\omega_1+n$ が可算共終なので共終数は $\aleph_0$ である。逆に $\omega+\omega_1=\omega_1$ である。実際各 $\alpha<\omega_1$ の $\omega+\alpha$ は可算で、かつ $\alpha$ 以上なので、その上限は $\omega_1$ に一致する。後者の共終数は $\aleph_1$ である。
:::

:::example[少数の有界集合の和]
$\omega_1$ の可算個の有界部分集合 $A_n$ に対し、各上界 $\beta_n<\omega_1$ を一つ取る。可算個の可算順序数の和は可算なので $\sup_n(\beta_n+1)<\omega_1$ である。この上限が全体の上界となる。$\aleph_\omega$ で $A_n=\aleph_n$ とすると各段は有界だが和は全体になる。同じ添字の濃度でも基数の正則性によって結論が変わる。
:::

### 共終数の二段の計算

$\alpha$ が極限順序数、$\beta>0$ とする。$\beta\alpha=\sup_{\xi<\alpha}\beta\xi$ の共終数は $\operatorname{cf}(\alpha)$ に等しい。上界は $\alpha$ の狭義増加共終列を右乗法で送れば得られる。逆に $\beta\alpha$ の共終集合 $D$ を取る。各 $d\in D$ について $d<\beta\xi_d$ を満たす最小の $\xi_d<\alpha$ を取れる。$\xi_d$ たちが $\alpha$ の有界部分に限られるなら、$D$ もそれに対応する $\beta\xi$ 以下に抑えられる。$\beta\xi<\beta\alpha$ なので矛盾する。従って $\{\xi_d:d\in D\}$ は共終であり、$|D|\ge\operatorname{cf}(\alpha)$ となる。$\xi\mapsto\beta+\xi$ の連続性と狭義増加性から、同じ方法で $\operatorname{cf}(\beta+\alpha)=\operatorname{cf}(\alpha)$ も従う。

### 共終数を添字と値域で読み替える

極限順序数 $\alpha$ の共終数は、狭義増加写像 $f:\delta\to\alpha$ で像が共終になる最小の順序数 $\delta$ としても取れる。もし $\delta$ が始順序数でなければ、$f$ の像は濃度 $|\delta|<\delta$ の共終集合になる。逆に最小濃度の共終集合を列挙し、初期区間の像がまだ有界であることを使って、その列挙を上から追い越す狭義増加写像を作れる。各段では列挙した値、および先行値全部より大きい点を選び、選んだ点が $\alpha$ 未満であるには $\alpha$ が極限であることを用いる。

:::example[添字の長さと共終数が一致しない列]
$\lambda=\omega_1+\omega$ に対し、恒等写像は長さ $\lambda$ の共終列である。しかし最後の部分だけを取り $\omega_1+n$ とすると長さ $\omega$ で共終である。共終数は利用できる共終列のうち最小の長さなので、長い列を一つ示しただけでは決まらない。$\lambda$ は極限を持つため有限の列では共終になれず、従って $\operatorname{cf}(\lambda)=\aleph_0$ になる。
:::

### 強極限基数の指数を共終数で分ける

無限基数 $\kappa$ が :term[強極限基数]{en="strong limit cardinal" yomi="きょうきょくげんきすう"} であるとは、各基数 $\mu<\kappa$ について $2^\mu<\kappa$ となることをいう。共終数より短い列は一つの区間に値域が収まる。共終数以上の長さを使うと、基数全体の部分集合を区間ごとに符号化できる。

:::proposition[強極限基数の二つの指数範囲]{#prop-logic-strong-limit-powers}
$\kappa$ を無限強極限基数とする。$0<\lambda<\operatorname{cf}(\kappa)$ では $\kappa^\lambda=\kappa$ であり、$\operatorname{cf}(\kappa)\le\lambda<\kappa$ では $\kappa^\lambda=2^\kappa$ である。
:::

:::proof
$\lambda<\operatorname{cf}(\kappa)$ なら、任意の関数 $\lambda\to\kappa$ の値域は $\kappa$ のある $\alpha<\kappa$ に含まれる。各 $\alpha$ に対して、そのような関数のグラフは $\lambda\times\alpha$ の部分集合である。この積の濃度は $\kappa$ 未満なので、強極限性によって候補数は $\kappa$ 未満である。すべての $\alpha<\kappa$ について合わせても高々 $\kappa$ 個であり、定数関数から逆の不等式が得られる。

次に $\delta=\operatorname{cf}(\kappa)\le\lambda<\kappa$ とする。$\kappa$ を濃度 $\kappa$ 未満の互いに素な $\delta$ 個の初期区間差に分割する。狭義増加な共終列の各段でまだ入っていない区間を取ればそのような分割になる。各ブロック $B_i$ の冪集合は濃度 $\kappa$ 未満なので、$\mathcal P(B_i)$ を $\kappa$ に入れる単射を一つずつ選ぶ。任意の $A\subseteq\kappa$ は各 $B_i$ との共通部分から復元できる。従ってこの単射で区間別の集合を符号にし、$\mathcal P(\kappa)$ を $\kappa^\delta$ に単射で移せる。これで $2^\kappa\le\kappa^\delta\le\kappa^\lambda$ となる。逆に $\kappa\le2^\kappa$ と $\kappa\lambda=\kappa$ を使えば $\kappa^\lambda\le(2^\kappa)^\lambda=2^\kappa$ である。
:::

:::example[最初のベート極限の可算指数]
$\kappa=\beth_\omega$ は強極限で共終数が $\omega$ である。従って $\kappa^n=\kappa$ が任意の正の有限 $n$ に対して成立する一方、$\kappa^\omega=2^\kappa>\kappa$ になる。後の等式では、$\mathcal P(\kappa)$ の各元を $\beth_n$ 以下の可算個の区間に分け、その情報を一つの $\omega$ 列へ圧縮する。$\kappa$ のすぐ次の基数がこの値になるかどうかは、さらにその点での GCH に依存する。
:::

### 特異基数を共終な積で越える

$\kappa$ を特異基数、$\delta=\operatorname{cf}(\kappa)$ とする。$\kappa$ の下に長さ $\delta$ の共終列がある。特異基数は後続基数ではないので、その列から各順序数を超える基数を $\kappa$ 未満で取り、$(\kappa_i)_{i<\delta}$ が $\kappa$ に共終になるようにできる。非交和の計算で $\sum_{i<\delta}\kappa_i=\kappa$ である。各 $\kappa_i<\kappa_i^+<\kappa$ だから、König の不等式で
$$
\kappa=\sum_{i<\delta}\kappa_i
    <\prod_{i<\delta}\kappa_i^+
    \le\kappa^\delta
$$
となる。従って特異基数に共終数と同じ指数を与えると、値は必ず底より真に増える。強極限性がある場合にはこの不等式を強めて、上の命題のように $\kappa^\delta=2^\kappa$ まで計算できる。強極限性を使わない段階では、各区間の任意の部分集合を一つの $\kappa$ 未満の記号で表せる保証がないため、この等式は同じ議論では従わない。

### 可算でない共終数を計算する

非零極限順序数 $\alpha$ に対し、順序型を指定した共終写像を調べる。$f:\delta\to\alpha$ が狭義増加かつ共終で、$\delta$ に狭義増加な共終写像 $g:\nu\to\delta$ があれば、合成 $f\circ g$ は $\alpha$ に共終である。実際、任意の $x<\alpha$ に対し $x\le f(i)$ となる $i<\delta$ を選び、$i\le g(j)$ となる $j$ を取れば $x\le f(g(j))$ である。逆に、$A\subseteq\alpha$ が共終なら、各 $a\in A$ を $a\le f(i)$ となる最小の $i$ に送る。その像が $\delta$ で有界なら $A$ は $f(i)$ 以下に抑えられ、$\alpha$ が極限であることに反する。従って
$$
\operatorname{cf}(\alpha)=\operatorname{cf}(\delta)
$$
である。ここでは $\delta$ も非零極限である。もし最大添字があれば、$f$ の像が $\alpha$ の一元で上から抑えられるためである。

従って $\operatorname{cf}(\aleph_\lambda)=\operatorname{cf}(\lambda)$、$\operatorname{cf}(\beth_\lambda)=\operatorname{cf}(\lambda)$ が非零極限の $\lambda$ で成り立つ。例えば $\aleph_{\omega_1}$ には $\omega_1$ 長の共終列があり、どの可算列も共終にならない。仮に可算集合が共終なら、それを超えるために必要なアレフの添字だけを集めても可算であり、その上限が $\omega_1$ 未満なので元の集合も有界になる。一方 $\aleph_{\omega_1}$ は $\aleph_{\omega_1+1}$ と違って極限基数である。共終数が非可算という事実だけでは正則性は従わない。

### 関数の値域が有界になる場合を計算する

$\lambda$ を無限基数とし、$\kappa$ は $\lambda<\operatorname{cf}(\kappa)$ を満たす無限基数とする。$f:\lambda\to\kappa$ の値域は $\kappa$ に共終になれない。従って各 $f$ は、ある $\alpha<\kappa$ への写像として数えられる。各 $\alpha$ ごとの個数は $|\alpha|^\lambda$ なので、全体は高々 $\kappa$ 個の集合の和であり、
$$
\kappa^\lambda=\max\left(\kappa,\sup_{\mu<\kappa}\mu^\lambda\right)
$$
となる。右の上限では $\mu$ を基数に限る。上界は和の評価による。下界には定数写像で $\kappa$ を入れ、また各 $\mu$ の写像を $\kappa$ への写像とみなせばよい。

特に $\kappa=\nu^+$ で $0<\lambda\le\nu$、$\nu$ は無限基数とする。$\kappa$ は正則なので同じ有界性が使える。$\alpha<\nu^+$ の濃度は $\nu$ 以下だから
$$
(\nu^+)^\lambda=\max(\nu^+,\nu^\lambda)
$$
となる。この証明は $\lambda$ が正の有限基数の場合にもそのまま通る。例えば $\aleph_1^{\aleph_0}=\max(\aleph_1,\aleph_0^{\aleph_0})=2^{\aleph_0}$ である。この等式は CH を仮定しない。一方 $(\aleph_2)^{\aleph_1}=\max(\aleph_2,(\aleph_1)^{\aleph_1})=2^{\aleph_1}$ であり、後者では $\aleph_1$ の全列を二値の冪で挟んでいる。

### 正則性を関数の有界性に移す

$\kappa$ が無限基数なら、正則であることは「$\lambda<\kappa$ と $f:\lambda\to\kappa$ に対し、その像は $\kappa$ で有界である」と同値である。正則な場合に非有界な像があれば、濃度高々 $\lambda$ の共終集合になる。逆に特異なら、共終列の添字である $\operatorname{cf}(\kappa)$ は $\kappa$ 未満であり、反例の関数をその列から作れる。従って正則な $\kappa$ で複数の小さい範囲を一つの境界で覆うには、各範囲に一つずつ境界を割り当て、その境界の集合自体を同じ考え方で抑えればよい。

:::example[後続基数の正則性とアレフ添字]
$\kappa=\aleph_{\omega+1}$ の直前の基数は $\aleph_\omega$ である。仮に長さ $\omega_1$ の列が $\kappa$ に共終なら、その各項の初期区間は濃度高々 $\aleph_\omega$ であり、$\omega_1$ 個の和も濃度高々 $\aleph_\omega$ になる。この和では $\kappa$ を覆えない。これに対し $\aleph_\omega$ 自身は $\aleph_n$ の可算列で覆える。直前の基数が特異でも、その後続の正則性は失われない。
:::

## 一般連続体仮説

:term[連続体仮説]{en="continuum hypothesis" yomi="れんぞくたいかせつ"} CH は $2^{\aleph_0}=\aleph_1$ である。実数の濃度が $2^{\aleph_0}$ であることは :ref[thm-sets-continuum] による。:term[一般連続体仮説]{en="generalized continuum hypothesis" yomi="いっぱんれんぞくたいかせつ"} GCH は各無限基数 $\kappa$ で $2^\kappa=\kappa^+$ とする主張である。GCH が成り立てば、後続段で同じ値を選び、極限段で同じ上限を取るので $\beth_\alpha=\aleph_\alpha$ が帰納的に従う。逆にこの等号が全添字で成り立てば、$\beth_{\alpha+1}=2^{\aleph_\alpha}=\aleph_{\alpha+1}$ から GCH が従う。相対無矛盾性は [構成可能集合と相対無矛盾性](/math/logic-set-theory/constructible/) と [強制法と独立性証明](/math/logic-set-theory/independence/) で扱う。

### 連続体の候補を制約する

Cantor の定理で $2^{\aleph_0}>\aleph_0$ である。さらに前節の共終数評価から $2^{\aleph_0}\ne\aleph_\omega$ である。実際 $\aleph_\omega$ の共終数は $\aleph_0$ だが連続体の共終数はそれより大きい。この排除には CH を仮定していない。CH は許される候補のうち最初の $\aleph_1$ を指定する等式である。

:::example[CH が答える範囲]
CH の下では $\aleph_0<|X|<|\R|$ を満たす集合は存在しない。右辺が $\aleph_1$ となり、その定義で間に基数がないためである。CH から直ちに $2^{\aleph_1}=\aleph_2$ を計算することはできない。その第二の等式は GCH の別の実例であり、最初の等式からは単調性と Cantor の不等式を適用して下界 $\aleph_2$ を得るところまでである。
:::

:::example[指数法則で読めること]
もし $2^{\aleph_0}=\aleph_2$ なら、$(\aleph_2)^{\aleph_0}=(2^{\aleph_0})^{\aleph_0}=2^{\aleph_0}=\aleph_2$ になる。ここでは二つの可算添字の積が可算であることを使った。一方 $(\aleph_2)^{\aleph_1}=2^{\aleph_1}$ となり、連続体の値だけを代入しても最後の冪は別の基数演算として残る。
:::

### 連続体と順序型を分ける

CH は $\mathcal P(\omega)$ の元を $\omega_1$ で添字付けられるという基数の主張である。通常の実数の順序が $\omega_1$ の整列順序と同型になるわけではない。例えば $\R$ の通常順序には $\{1,1/2,1/3,\ldots\}$ の最小元がないが、$\omega_1$ の各非空部分集合は最小元を持つ。CH の全単射で順序を移して得られる整列順序は、通常の大小とは別の関係となる。

### 冪の値を強制するために必要な基数算術

無限基数 $\kappa\le\lambda$ について、$\kappa$ から二値への各写像は、残りを $0$ で埋めて $\lambda$ への写像に延長できる。この単射により $2^\kappa\le2^\lambda$ である。また $\operatorname{cf}(2^\kappa)>\kappa$ は König の不等式から従う。従ってどの基数でも、冪の値を候補として挙げる時点で、単調性とその値の共終数を確認しなければならない。

例えば $2^{\aleph_0}=\aleph_2$、$2^{\aleph_1}=\aleph_1$ という二つの等式は矛盾する。第二の等式を $\aleph_2$ に直せば単調性には反しないが、それだけで実際に両立することの証明にはならない。また $2^{\aleph_0}=\aleph_{\omega_1}$ は共終数条件には反しない。実際 $\operatorname{cf}(\aleph_{\omega_1})=\aleph_1>\aleph_0$ である。このように候補を算術的に排除する議論と、その候補を持つモデルを構成する議論には違う段階がある。後者の具体例は [強制法と独立性証明](/math/logic-set-theory/independence/) で扱う。正則非可算基数に限れば、この二条件は強制法による実現にも使える。Easton の定理では、GCH を満たす地のモデルで正則非可算基数に基数値を割り当てる定義可能なクラス関数 $F$ が、単調性と $\operatorname{cf}(F(\kappa))>\kappa$ を満たすとする。このとき、基数と共終数を保ち、全ての正則非可算 $\kappa$ で $2^\kappa=F(\kappa)$ とする強制拡大を得られる。全正則点を同時に扱うクラス強制の構成と保存の証明は :cite[jech-set] を参照する。特異基数はこの定理で値を指定する領域に入っていない。

特異基数では、下の基数の冪が上の値を制約することがある。例えば非可算共終数の特異基数より下で GCH が成り立つ場合、その特異基数でも GCH が成り立つ。証明は Fodor の補題を使うので、[無限組合せ論](/math/logic-set-theory/infinite-combinatorics/) で扱う。

### GCH の下で基数の冪を場合分けする

$\kappa,\lambda$ を無限基数とする。GCH の下で $\kappa\le\lambda$ なら $2\le\kappa\le2^\lambda$ なので $\kappa^\lambda=2^\lambda=\lambda^+$ である。次に $\lambda<\kappa$ とする。$\kappa=\nu^+$ が後続基数なら $\lambda\le\nu$ であり、GCH と指数法則により
$$
\kappa^\lambda=(2^\nu)^\lambda=2^{\nu\lambda}=2^\nu=\kappa
$$
となる。この場合は $\operatorname{cf}(\kappa)=\kappa>\lambda$ である。

$\kappa$ が極限基数なら、各無限 $\mu<\kappa$ で $2^\mu=\mu^+<\kappa$ である。有限の $\mu$ でも冪は有限なので $\kappa$ 未満であり、これで $\kappa$ は強極限となる。従って $\lambda<\operatorname{cf}(\kappa)$ のときは $\kappa^\lambda=\kappa$、$\operatorname{cf}(\kappa)\le\lambda<\kappa$ のときは $\kappa^\lambda=2^\kappa=\kappa^+$ である。このように GCH の下では、底と指数の大小および底の共終数から無限基数の冪を計算できる。例えば GCH を仮定すると $\aleph_{\omega_1}^{\aleph_0}=\aleph_{\omega_1}$ だが、$\aleph_{\omega_1}^{\aleph_1}=\aleph_{\omega_1+1}$ になる。底は同じであり、後者の指数が共終数に達したところで値が変わる。

## 到達不能基数

$\kappa$ が非可算正則な極限基数なら弱い意味で到達不能である。さらに全 $\lambda<\kappa$ について $2^\lambda<\kappa$ を要求し、非可算正則と合わせた条件を強い意味での :term[到達不能基数]{en="inaccessible cardinal" yomi="とうたつふのうきすう"} と呼ぶ。後者の冪に関する条件を強極限性という。強極限なら後続基数にはなれない。もし $\kappa=\lambda^+$ なら $2^\lambda\ge\lambda^+$ だからである。GCH の下では逆に、極限基数 $\kappa$ なら各無限 $\lambda<\kappa$ について $2^\lambda=\lambda^+<\kappa$ である。

:::example[仮定ごとの失敗]
$\aleph_1$ は正則だが後続基数なので強極限ではない。$\beth_\omega$ は全 $\lambda<\beth_\omega$ に対し、ある有限 $n$ で $\lambda<\beth_n$ となり $2^\lambda\le\beth_{n+1}<\beth_\omega$ である。しかし $\beth_n$ の列が可算共終なので正則でない。二条件は別々に調べる必要がある。
:::

### 強極限性と正則性を合わせる

強い到達不能基数 $\kappa$ に対し、$1<\lambda<\kappa$ と $\mu<\kappa$ なら $\lambda^\mu<\kappa$ である。グラフへの単射で $\lambda^\mu\le2^{|\lambda\times\mu|}$ となり、積の濃度 $\max(\lambda,\mu)$ 以下（両方有限なら有限）は $\kappa$ 未満である。さらに $\delta<\kappa$ 個の $\kappa$ 未満の基数を取ると、正則性からその上限も $\kappa$ 未満なので、非交和の濃度も $\kappa$ 未満になる。

:::example[可算段のベートをまとめる]
$\beth_0<\kappa$ であり、強極限性で $\beth_n<\kappa$ なら $\beth_{n+1}<\kappa$ となる。さらに $\kappa$ は非可算正則なので、その可算上限 $\beth_\omega$ も $\kappa$ 未満である。強極限性は後続段で使い、正則性は極限段で使っている。
:::

:::example[共終数から正則性を検査する]
強極限基数 $\beth_\omega$ では個々の $\beth_n$ の冪がまだ $\beth_\omega$ 未満になる。それでも濃度 $\beth_n$ の集合を可算個並べた非交和は濃度が $\beth_\omega$ に達する。従って到達不能段の閉性を証明するときに、すべての小濃度の冪が小さいという条件だけで、可算個の小集合の和まで抑えることはできない。
:::

### 到達不能基数の下で反復を閉じる

$\kappa$ が強到達不能であれば、全 $\alpha<\kappa$ に対して $\beth_\alpha<\kappa$ である。基底で $\aleph_0<\kappa$、後続で強極限性、極限で正則性を使う。従って $\beth_\kappa=\sup_{\alpha<\kappa}\beth_\alpha\le\kappa$ である。逆に狭義増加な $\kappa$ 項の列は $\kappa$ 未満の順序数には入らないので、上限は $\kappa$ になる。同様に $\aleph_\alpha$ は $\beth_\alpha$ 以下であり、超限帰納法で後続は Cantor の定理、極限は両列の上限から従う。従って $\aleph_\kappa=\beth_\kappa=\kappa$ である。この固定点性は到達不能性から出るが、固定点だけで共終数が正則になるとは言えない。

### 正則性を和の閉性として使う

無限基数 $\kappa$ の正則性は、濃度が各々 $\kappa$ 未満である $\kappa$ 未満個の集合の和が、いつも濃度 $\kappa$ 未満になることと同値である。正則なら、それらのサイズと添字集合のサイズをひとつの $\kappa$ 未満の基数で抑えられる。その平方はそれ自体に等しいか有限なので、和も同じ範囲に収まる。逆に $\kappa$ が特異なら、共終な $\operatorname{cf}(\kappa)$ 個の初期区間の和が $\kappa$ になる。各区間の濃度は $\kappa$ 未満だから、和の閉性が失敗する。この言い換えを使えば、反復で各段の候補数を抑えた後、まとめた領域の大きさを正則性で評価できる。

### ベート固定点に正則性を足す

非可算基数 $\kappa$ について、強到達不能であることは、正則かつ $\beth_\kappa=\kappa$ であることと同値である。順方向は前の計算による。逆に固定点なら、任意の $\lambda<\kappa$ に対して $\lambda\le\beth_\alpha$ となる $\alpha<\kappa$ がある。$\kappa$ は無限基数だから順序数として極限であり、$\alpha+1<\kappa$ である。従って
$$
2^\lambda\le\beth_{\alpha+1}<\beth_\kappa=\kappa.
$$
よって $\kappa$ は強極限となり、仮定の正則性と合わせて到達不能である。弱到達不能の場合は、条件 $\beth_\kappa=\kappa$ を $\aleph_\kappa=\kappa$ に替える。実際、非可算正則な極限基数の下で $\aleph_\alpha<\kappa$ を帰納するとき、後続では極限基数であること、極限では正則性を使う。逆に $\aleph_\kappa=\kappa$ なら $\kappa$ は極限添字での極限基数であるから、正則性と合わせて弱到達不能である。

:::example[固定点は特異になり得る]
$\alpha_0=\omega$ とし、$\alpha_{n+1}=\beth_{\alpha_n}$ とする。$\beth_\omega>\omega$ であり、$\beth$ は狭義増加なので列は各段で狭義に増える。その上限 $\kappa$ は極限順序数で、$\beth$ の連続性により
$$
\beth_\kappa=\sup_n\beth_{\alpha_n}=\sup_n\alpha_{n+1}=\kappa
$$
である。これは基数の上限であるため基数であり、上の議論で強極限である。一方、列が $\kappa$ に共終で有限共終集合はないため、$\operatorname{cf}(\kappa)=\aleph_0$ である。従ってこの固定点は強到達不能ではない。固定点を得る反復と、その極限値が正則かどうかの判定は別に行う必要がある。
:::

### 小さい集合の反復をどこまで続けられるか

$\kappa$ が強到達不能なら、濃度 $\kappa$ 未満の集合の冪集合を取っても濃度は $\kappa$ 未満に留まる。また長さ $\delta<\kappa$ の列 $(X_i)$ で、各 $X_i$ の濃度が $\kappa$ 未満なら、正則性により濃度の一様な上界 $\mu<\kappa$ がある。非交和は $|\delta|\mu<\kappa$ で抑えられる。直積は各座標を $\mu$ に入れて $\mu^{|\delta|}$ で抑えられ、指数のグラフにより $\kappa$ 未満となる。空の座標があれば直積は空であり、この結論は同様に成立する。従って小さい添字集合で、和・積・冪集合による操作を反復できる。長さが $\kappa$ に達する場合には、例えば $\kappa$ 個の単元集合の非交和だけで濃度が $\kappa$ になるため、長さの制限は外せない。

この閉性と累積階層の各段の濃度の関係は、[集合論のモデルと絶対性](/math/logic-set-theory/models-of-set-theory/) で使う。ここで分かるのは濃度の評価であり、ある階層段が公理を満たすためには、その段に必要な集合が入ることも所属に沿って確認することになる。

### 指数の単調性と極限での跳び

底を上げると、各関数に底の単射を合成できる。指数を上げるときは底に一点を固定し、元の定義域の外ではその値を返すように関数を延長できる。従って底が非空なら、冪は両変数について弱く単調である。底が空の場合、$0^0=1$ だが $0^\lambda=0$ が $\lambda>0$ で成立するので、指数の弱い単調性は成り立たない。さらに $1^\lambda=1$ は全 $\lambda$ で成り立つ。ここで $0^0$ を $1$ とするのは、空集合から空集合への関数が空関数一つだからである。

基数の冪を極限で計算する場合、順序数の冪のように先行値の上限を取ることはできない。例えば各有限 $n>0$ で $\aleph_0^n=\aleph_0$ だが、$\aleph_0^\omega=2^{\aleph_0}>\aleph_0$ である。有限の各制限には可算個の選択肢しかないが、その選択をすべての段で整合的に続ける列は非可算個ある。底を極限に取る場合にも跳びが起こる。例えば $\kappa=\beth_\omega$ なら、固定した $n$ に対して $\beth_n^\omega\le2^{\beth_n}=\beth_{n+1}<\kappa$ であるため、$\sup_n\beth_n^\omega\le\kappa$ である。一方 $\kappa^\omega=2^\kappa>\kappa$ である。ここでは底の極限と指数の計算を入れ替えると、値域がどの先行段にも収まらない関数が漏れる。

:::example[値域に境界がない列を区間ごとに符号化する]
$\kappa=\beth_\omega$ とする。$B_0=\beth_0$、$B_{n+1}=\beth_{n+1}\setminus\beth_n$ とすると、これらは $\kappa$ の分割を与える。$B_n$ の濃度は高々 $\beth_n$ であり、部分集合の候補は高々 $\beth_{n+1}<\kappa$ 個である。その候補を一つの $\kappa$ の元で符号化して $A\subseteq\kappa$ を $(\operatorname{code}_n(A\cap B_n))_{n<\omega}$ に送る。この写像は各区間の情報を復元できるので単射である。たとえ区間ごとに小さな候補集合を使っても、符号の列全体は単一の先行区間に制限されない。このため $\kappa^\omega$ の上界を $\sup_n\beth_n^\omega$ で与えられない。
:::

## 演習問題

::::exercise{level=basic}
二つと三つの元を持つ集合の非交和と直積の濃度を求めよ。
:::solution
非交和では重複しないので $2+3=5$、直積では二つの第一座標ごとに三つの第二座標があり $2\cdot3=6$ となる。
:::
::::

::::exercise{level=basic}
$\N\times\{0,1\}$ を $\N$ へ全単射で移せ。
:::solution
$(n,i)\mapsto2n+i$ は偶奇と商から逆を作れるので全単射である。
:::
::::

::::exercise{level=basic}
$\operatorname{cf}(\omega+7)$ を求めよ。
:::solution
$\omega+7$ は最大元 $\omega+6$ を持つ。これ一元が共終であり空集合は共終でないため共終数は $1$ である。
:::
::::

::::exercise{level=basic}
$\operatorname{cf}(\aleph_\omega)$ を求めよ。
:::solution
$\{\aleph_n:n<\omega\}$ は共終なので高々可算である。極限順序数に有限の共終集合はないため、ちょうど $\aleph_0$ である。
:::
::::

::::exercise{level=advanced}
$\kappa$ が無限なら $(2^\kappa)^\kappa=2^\kappa$ を示せ。
:::solution
関数の組み替えで $(2^\kappa)^\kappa=2^{\kappa\cdot\kappa}$ である。無限基数の積 $\kappa\cdot\kappa=\kappa$ を代入する。
:::
::::

::::exercise{level=advanced}
CH が偽なら $\beth_1$ はどのアレフ以上か。
:::solution
Cantor の定理で $\beth_1\ge\aleph_1$ である。等号が CH なので、偽なら厳密に大きい。次の始順序数の最小性により $\beth_1\ge\aleph_2$ となる。
:::
::::

::::exercise{level=advanced}
非可算正則基数 $\kappa$ に対し、$\kappa$ 未満の個数の $\kappa$ 未満の順序数の上限が $\kappa$ 未満であることを示せ。
:::solution
上限が $\kappa$ なら、その順序数の集合は共終である。共終数の最小性から濃度は $\operatorname{cf}(\kappa)=\kappa$ 以上となり矛盾する。
:::
::::

::::exercise{level=advanced}
$\beth_\omega$ が強極限だが正則でないことを示せ。
:::solution
任意の小さい $\lambda$ に対し $\lambda<\beth_n$ を選ぶと $2^\lambda\le\beth_{n+1}<\beth_\omega$ である。一方、厳密な増加列 $(\beth_n)$ が共終であり、有限共終集合はないため共終数は $\aleph_0$ となる。
:::
::::

::::exercise{level=advanced}
$\kappa$ を非可算強極限基数とし、$\operatorname{cf}(\kappa)=\omega_1$ とする。$\kappa^{\aleph_0}$ と $\kappa^{\aleph_1}$ を比較せよ。
:::solution
第一の指数は共終数未満なので、可算列の値域は一つの $\kappa$ 未満の初期区間に入る。:ref[prop-logic-strong-limit-powers] で $\kappa^{\aleph_0}=\kappa$ である。第二では共終数と同じ長さまで使える。順序数 $\omega_1$ は後続基数であって強極限ではないので、仮定から $\aleph_1<\kappa$ である。従って同じ命題の第二の場合により $\kappa^{\aleph_1}=2^\kappa>\kappa$ となる。
:::
::::

::::exercise{level=advanced}
二つの無限基数 $\lambda,\mu$ に対して $\lambda\le2^\mu$ かつ $2\le\lambda$ とする。$\lambda^\mu=2^\mu$ を示せ。また $\lambda$ または $\mu$ の一方が有限になる場合に使えなくなる等式を指摘せよ。
:::solution
$2\le\lambda$ より $2^\mu\le\lambda^\mu$ である。一方 $\lambda\le2^\mu$ から $\lambda^\mu\le(2^\mu)^\mu=2^{\mu\cdot\mu}$ であり、$\mu$ が無限なので $\mu\cdot\mu=\mu$ である。従って等号になる。この証明では $\lambda$ の無限性を使っておらず、$\lambda$ が有限でも $2\le\lambda\le2^\mu$ なら同じである。しかし有限の $\mu>1$ では $\mu\cdot\mu=\mu$ が失敗する。例えば $\mu=2,\lambda=3$ なら $3^2=9$ と $2^2=4$ は異なる。
:::
::::

::::exercise{level=advanced}
$\kappa$ を無限基数とする。濃度が $\kappa$ 以下の空でない集合を $\kappa^+$ 個足した非交和の濃度を求めよ。各集合が空でもよい場合に残る結論も述べよ。
:::solution
各座標から一点ずつ取れるので $\kappa^+$ が下界である。各集合の濃度は $\kappa$ 以下であり、非交和は $\kappa^+\times\kappa$ へ単射を持つので上界も $\kappa^+$ である。空集合を許す場合は、非空な座標からしか下界が得られない。上界は引き続き $\kappa^+$ であり、例えば全座標が空なら濃度は $0$ になる。
:::
::::

## 確認テスト

::::quiz{answer="8" tol="0"}
三元集合から二元集合への写像は何個か。
:::explanation
三つの入力ごとに二つの値を独立に選ぶ。従って $2^3=8$ 個である。
:::
::::

::::quiz
CH と同値な等式を選べ。

- [ ] $\aleph_1=\aleph_0$
- [x] $\beth_1=\aleph_1$
- [ ] $\beth_1=\aleph_0$

:::explanation
CH は自然数の冪集合の濃度が最小の非可算基数と等しいという主張である。$\beth_1=2^{\aleph_0}$ なので第二の選択肢の等式を選ぶ。
:::
::::

::::quiz{answer="1" tol="0"}
順序数 $7$ の共終数はいくつか。
:::explanation
最大元は $6$ であり、$\{6\}$ が共終である。$7$ は空でないので共終数は一である。
:::
::::

::::quiz
非可算な強到達不能基数の条件を選べ。

- [ ] 後続基数
- [ ] 共終数が可算であること
- [x] 正則かつ強極限

:::explanation
共終数が自身と等しいことが正則性であり、各小さい基数の冪が自身未満であることが強極限性である。この二条件をともに用いる。
:::
::::


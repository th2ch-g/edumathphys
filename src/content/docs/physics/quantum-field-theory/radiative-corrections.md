---
title: "輻射補正入門"
description: "QED の一ループ補正を整理し、異常磁気モーメント、電荷の変化、赤外相殺と Lamb シフトを導く。"
notation:
  - symbol: 'd=4-2\epsilon'
    meaning: 次元正則化の時空次元
  - symbol: 'F_1,F_2'
    meaning: 電磁頂点の形状因子
  - symbol: '\Pi_R'
    meaning: 零運動量で減算した真空偏極
---

## この章の目標

- Feynman パラメータを用いてループ積分を評価する。
- QED の頂点、真空偏極、自己エネルギーを分けて扱う。
- 軟光子の実・仮想補正と束縛準位の先頭対数を求める。

自然単位系 $\hbar=c=1$ を用い、Minkowski 空間の計量は $\eta=\diag(+1,-1,-1,-1)$。Wick 回転後の積分では Euclid の二乗を使うので、同じ外部ベクトルの空間的な移行は $Q^2=-q^2>0$ と書く。

## 1 ループ積分と Feynman パラメータ

分母 $A,B$ の実部が正なら $A^{-1}=\int_0^\infty \dd u\,e^{-uA}$。二本に使い $u=xs,v=(1-x)s$ と置くと Jacobian は $s$。従って
$$
\frac1{A^aB^b}=\frac{\Gamma(a+b)}{\Gamma(a)\Gamma(b)}
\int_0^1\frac{x^{a-1}(1-x)^{b-1}\dd x}{[xA+(1-x)B]^{a+b}}
\label{eq-qft-feynman-parameter}
$$
となる。一般の正の $a,b$ には $u^{a-1}/\Gamma(a)$ を含む Laplace 積分を使う。Minkowski 分母では共通の $i0$ を保って解析接続する。例えば $A=k^2-m^2$、$B=(k+q)^2-m^2$ では $\ell=k+(1-x)q$、分母は $\ell^2-\Delta$、$\Delta=m^2-x(1-x)q^2$。閾値 $q^2=4m^2$ で最初に $\Delta$ が零になる。頂点の三分母には同様に $x+y+z=1$ の単体積分を使う。

### 低運動量での二点積分の展開

同じ質量の二分母の差引き後の核 $B(s)=-\int_0^1\dd x\,\log[1-sx(1-x)/m^2]/(16\pi^2)$ を扱う。$|s|<4m^2$ では対数級数が $x$ に一様に収束する。従って
$$
B(s)=\frac{1}{16\pi^2}\sum_{n=1}^\infty\frac1n
\left(\frac{s}{m^2}\right)^n\frac{(n!)^2}{(2n+1)!}.
$$
積分係数は $\int_0^1 x^n(1-x)^n\dd x$ を $n$ 回部分積分して得る。例えば二次までは $B(s)=(s/(6m^2)+s^2/(60m^4))/(16\pi^2)$。$s=m^2$ で最初の二項の比は $1/10$ である。$s>4m^2$ では級数を使わず、分母の $i0$ から切断をまたぐ値を決める。

### 多分母の単体積分と分岐の始点

$N$ 個の Euclid 分母 $A_i$ に対し、$\Re A_i>0$、$a_i>0$ とする。それぞれのパラメータを $u_i=s x_i$、$\sum_i x_i=1$ に置くと、例えば $x_N=1-\sum_{i<N}x_i$ を使った体積要素は $s^{N-1}\dd s\,\dd x_1\cdots\dd x_{N-1}$ である。実際、Jacobian の最初の列は $(x_i)$、残りの列は $s$ 倍の $e_i-e_N$。最後の行へ他の行を加えれば、最初の列の最後が $1$、残りが零となる。従って
$$
\frac1{\prod_i A_i^{a_i}}
=\frac{\Gamma(\sum_i a_i)}{\prod_i\Gamma(a_i)}
\int_{x_i\ge0}\left(\prod_i\dd x_i\,x_i^{a_i-1}\right)
\frac{\delta(1-\sum_i x_i)}{(\sum_i x_i A_i)^{\sum_i a_i}} .
$$
例えば $N=3$、$a_i=1$ なら三角形の面積は $1/2$ で、$\Gamma(3)=2$ が係数に付く。全て同じ分母のときはこの二つが相殺する。Minkowski の移動後に三つの内部線が同時に殻へ近づく場合には、単体での質量項が零になる点を調べる。途中の分母にそれぞれ別の符号の虚部を入れると、元の時間順序の極には対応しなくなる。

:::example[二種類の質量を持つ泡]
質量 $m_1,m_2>0$ に対して $\Delta=xm_1^2+(1-x)m_2^2-sx(1-x)$。$0<x<1$ で $\Delta=0$ となる $s$ は $m_1^2/(1-x)+m_2^2/x$ である。この式から $(m_1+m_2)^2$ を引くと $(m_1 x-m_2(1-x))^2/[x(1-x)]$。従って最小は $x=m_2/(m_1+m_2)$、閾値は $s=(m_1+m_2)^2$ となる。両端ではどちらかの質量二乗が正なので、先に分母が零になるのは単体の内部である。
:::

### 分岐切断を二体の閾値と比べる

等質量の場合は $B(s)$ に $s+i0$ を入れる。$s>4m^2$ では $1-sx(1-x)/m^2$ が負となる区間の長さが $\sqrt{1-4m^2/s}$。この区間で $\log$ の虚部は $-\pi$ なので
$$
\Im B(s+i0)=\frac{1}{16\pi}\sqrt{1-\frac{4m^2}s}.
$$
上と下の縁の差は $2i$ 倍の虚部である。等質量の二体位相空間を全角積分した値は $\sqrt{1-4m^2/s}/(8\pi)$ なので、切断の差はこの位相空間に $i$ を掛けた量である。実部がどの定数だけずれるかは $B(0)=0$ という減算で決まるが、この切断の差は局所定数に依存しない。例えば $s=16m^2/3$ では速度は $1/2$、虚部は $1/(32\pi)$。

### 真空泡の Schwinger 表示としきい値の制約

二分母の積を $u,v$ の指数にしてから運動量 $k$ を積分する。Euclid で同じ質量 $m$ と外部運動量 $Q$ を使うと、平方完成した定数部は $(u+v)m^2+[uv/(u+v)]Q^2$。Gaussian で $u+v$ の $-d/2$ 乗が出る。$u+v=s$、$x=u/s$ の置換で $Q$ に掛かる係数は $s x(1-x)$。従って紫外は $s\downarrow0$、大距離は $s\to\infty$ に対応する。例えば $d=4$ では測度を含む $s$ の冪は $s^{-1}$。$Q=0$ の同じ積分を引けば、指数の差が $-s x(1-x)Q^2+O(s^2)$ になって紫外端が積分可能になる。一方、$Q^2\to-s_{\rm timelike}-i0$ と続けると指数が負に減衰する条件が $m^2-x(1-x)s_{\rm timelike}>0$。これを破る領域は、元の Euclid の正実パラメータ積分のままでは評価せず解析接続の対数で定める。

:::example[減算後の Frullani 積分]
$a,b>0$ なら $\int_0^\infty(e^{-sa}-e^{-sb})\dd s/s=\log(b/a)$。左辺の $b$ 微分は $1/b$、$b=a$ では零なので式が決まる。積分の収束は原点での指数の差が $O(s)$、無限遠で両指数が減衰するためである。これを $a=m^2+x(1-x)Q^2$、$b=m^2$ に使えば、次元正則化の極を展開せずとも、同じ条件で減算した泡の対数が得られる。二つの正則化の個々の裸積分は異なるが、同じ局所定数を差し引いたこの差は一致する。
:::

### 等質量の泡を閾値まで積分する

$0<s<4m^2$ の減算済みの泡で $v=\sqrt{4m^2/s-1}$ とし、$y=2x-1$ を導入する。すると
$$
1-\frac{s}{m^2}x(1-x)=\frac{y^2+v^2}{1+v^2}
$$
である。偶性で積分を $0<y<1$ へ移す。部分積分により $\int_0^1\log(y^2+v^2)\dd y=\log(1+v^2)-2+2v\arctan(1/v)$。従って $16\pi^2 B(s)=2-2v\arctan(1/v)$。$s\downarrow0$ では $v\to\infty$ なので $\arctan(1/v)=1/v-1/(3v^3)+\cdots$、先頭は $s/(6m^2)$。一方 $s\uparrow4m^2$ では $v\to0$、閾値で $16\pi^2 B=2$ だが、その直前の差は $-\pi v+O(v^2)$。従って値は有限でも $s$ で微分したときに閾値の特異性がある。

$s>4m^2$ では $v\to-i\beta$ と続ける。ただし $\beta=\sqrt{1-4m^2/s}$。同じ結果は $\int_0^1\log[(y^2-\beta^2-i0)/(1-\beta^2)]\dd y$ を $y=\beta$ で分割して求められる。実部の原始関数は
$$
y\log|y^2-\beta^2|-2y+\beta\log\left|\frac{y+\beta}{y-\beta}\right|.
$$
$y=\beta$ で二つの対数の特異性は相殺し、原始関数の左右の値は一致する。下端の値は零、上端では分母の $1-\beta^2$ の対数が相殺する。従って
$$
16\pi^2 B(s+i0)=2-\beta\log\frac{1+\beta}{1-\beta}+i\pi\beta.
$$
例えば $\beta=1/2$ では実部 $2-(\log3)/2$、虚部 $\pi/2$。この虚部は同じパラメータ積分で対数が負実軸を通る領域の長さと一致する。

## Wick 回転と次元正則化・Pauli–Villars 正則化

分母をまとめ、空間的外部運動量で $\Delta>0$ とする。移動後の $\ell^0$ の正エネルギー極は第四象限に、負エネルギー極は第二象限にある。正の実軸から正の虚軸へは第一象限を、負の実軸から負の虚軸へは第三象限を通るため、極を横切らない。時間成分を $\ell^0=i\ell_4$ へ回す :term[Wick 回転]{en="Wick rotation" yomi="うぃっくかいてん"} により、Euclid 分母と Jacobian $i$ を得る。積分路の弧を捨てるには、まず正則化するか十分な収束次数を取る。領域を物理的閾値まで延ばすには $i0$ を付けた解析接続が必要である。

:term[次元正則化]{en="dimensional regularization" yomi="じげんせいそくか"} では Gaussian 積分を解析的に続ける。$\Delta>0$ のとき Schwinger パラメータを導入し $\int\dd^d\ell\,e^{-s\ell^2}=(\pi/s)^{d/2}$ を使えば
$$
\mu^{2\epsilon}\int\frac{\dd^d\ell}{(2\pi)^d}\frac1{(\ell^2+\Delta)^n}
=\frac{\mu^{2\epsilon}}{(4\pi)^{d/2}}
\frac{\Gamma(n-d/2)}{\Gamma(n)}\Delta^{d/2-n}.
\label{eq-qft-loop-master}
$$
$n=2$ では $[1/\bar\epsilon-\log(\Delta/\mu^2)]/(16\pi^2)$ まで得る。$1/\bar\epsilon=1/\epsilon-\gamma_E+\log4\pi$ である。テンソル分子の $\ell_\mu\ell_\nu$ は回転不変性で $\delta_{\mu\nu}\ell^2/d$ へ置く。$d=4$ を先に代入すると極との積で生じる有限項を失う。

Pauli–Villars 正則化では重い質量 $M$ の伝播関数を引く。例えば Euclid で $1/(\ell^2+m^2)-1/(\ell^2+M^2)=(M^2-m^2)/[(\ell^2+m^2)(\ell^2+M^2)]$ となり、二つの逆二乗が相殺して逆四乗になる。必要な差引き回数は発散次数で決まる。対称性を保つには、関連するグラフで正則化を整合させる必要がある。

### 一つの差引きでどこまで収束するか

Euclid 四次元球カットオフで $I(a)=\int_{|\ell|<\Lambda}\dd^4\ell/(2\pi)^4(\ell^2+a)^{-1}$ とする。角面積 $2\pi^2$ と $u=\ell^2$ を使うと $I(a)=[\Lambda^2-a\log(1+\Lambda^2/a)]/(16\pi^2)$。従って一回の Pauli–Villars 差 $I(a)-I(b)$ は二次発散を消すが $(b-a)\log\Lambda^2$ が残る。さらに $a=1,b=4,c=9$ を同じ質量単位の二乗で取ったとき、係数 $(1,-8/5,3/5)$ は $\sum c_i=\sum c_i a_i=0$ を満たす。よってこの組合せは紫外で逆六乗に減衰し、有限部分は $\sum c_i a_i\log a_i/(16\pi^2)$ である。対数の尺度を変えても $\sum c_i a_i=0$ により結果は変わらない。

### Euclid のテンソル積分と全微分

$I_n(\Delta)=\mu^{2\epsilon}\int_\ell(\ell^2+\Delta)^{-n}$ と略記する。積分の収束する $n$ と $d$ で $\partial_{\ell_\mu}[\ell_\mu/(\ell^2+\Delta)^{n-1}]$ の積分は零である。従って
$$
d I_{n-1}=2(n-1)(I_{n-1}-\Delta I_n).
$$
左辺は分子の $\ell_\mu$ の微分、右辺は分母の微分である。:eqref[eq-qft-loop-master] では $\Gamma(n-d/2)=(n-1-d/2)\Gamma(n-1-d/2)$ により同じ関係を得る。解析接続の後にもこの恒等式が保たれるので、真空偏極のように発散積分に現れる全微分を整合して扱える。

:::example[$d$ を先に四へ置いた場合の有限差]
$\int_\ell \ell_\mu\ell_\nu/(\ell^2+\Delta)^3$ で、$1/d=1/4+\epsilon/8+\cdots$ を用いる。残りの積分 $I_2-\Delta I_3$ の極は $1/(16\pi^2\epsilon)$。従って $1/d$ を $1/4$ に置き換えると、$\delta_{\mu\nu}/(128\pi^2)$ の有限部分を失う。両者の極の係数は等しいため、誤りは有限な量を比較して初めて現れる。
:::

:::example[移動を許さない紫外カットオフの境界項]
四次元球 $|k|<\Lambda$ で $\int_k \partial_{k_\mu}[k_\nu/(k^2+\Delta)^2]$ を取る。境界で $k_\nu=\Lambda n_\nu$、外向き法線は $n_\mu$、面積要素は $\Lambda^3\dd\Omega_3$。$n_\mu n_\nu$ の角平均 $\delta_{\mu\nu}/4$ を使うと積分は $\delta_{\mu\nu}\Lambda^4/[32\pi^2(\Lambda^2+\Delta)^2]$。$\Lambda\to\infty$ で有限な $\delta_{\mu\nu}/(32\pi^2)$ が残る。従って発散する項を別々にシフトしてから足す操作は、この正則化では全微分の境界項を変える。電流保存を確認するなら、並進不変な正則化で全ての項を一緒に扱う。
:::

### 大質量の Dirac 場で真空偏極だけを正則化する

外部電磁場の二点を Dirac ループで計算し、同じ電荷で質量だけを $m$ から $M$ に替えた完全なループを引く。すなわち各伝播関数から個別に重い線を引くのではなく、同じ頂点を二個持つループ全体の差にする。一つのループの発散の先頭は質量に依存しないので、対称に積分すれば横方向の対数極は相殺する。横性を縮約の差から示す場合、ベクトル電流の一分母の差を軽い場と重い場で一緒に扱う。必要ならさらに正則化質量と重みを加えて $\sum_i c_i=\sum_i c_i M_i^2=0$ とし、各テンソルの大運動量の級数の先頭二項を消す。分子にも同じ各質量を使うため、$\not q$ を逆核の差に直す関係を変えずに済む。

例えば横性を保つ次元正則化でまず軽い・重い二点の差を作ると、極を引いた横成分は
$$
\Pi_{\rm PV}(-Q^2)=\frac{2\alpha}{\pi}\int_0^1\dd x\,x(1-x)
\log\frac{m^2+x(1-x)Q^2}{M^2+x(1-x)Q^2}
$$
になる。$Q=0$ の値を引くと、分子が軽い場の物理的な差を与え、重い場からは $O(Q^2/M^2)$ しか残らない。従って $M\to\infty$ で同じ有限な補正を再現する。このように電荷の減算を行う前の対数 $\log(m^2/M^2)$ は正則化に依存するが、同じ規格化条件の後では局所定数に吸収される。$M$ を有限の物理的な質量として残す場合は、正則化の差し引きと実在する粒子のループの加算で符号が異なる。

## 電子の頂点補正と異常磁気モーメント

### 磁気モーメントの有限な部分

殻上で $q=p'-p$ とし、頂点から電荷を除いた形を $\Gamma^\mu=\gamma^\mu F_1(q^2)+i\sigma^{\mu\nu}q_\nu F_2(q^2)/(2m)$ と書く。外脚が質量殻上にあるので :eqref[eq-qft-gordon] を適用できる。静的磁場でスピンに結合するのは $F_1+F_2$、電荷を $F_1(0)=1$ と規格化すれば $a=(g-2)/2=F_2(0)$ となる。

ループの光子に赤外調整用質量 $\lambda_\gamma$ を置く。二電子分母に $x,y$、光子に $1-x-y$ を割り当て、$\ell=k-xp'-yp$ を用いると分母の質量項は $m^2(x+y)^2-xyq^2+\lambda_\gamma^2(1-x-y)$。分子 $\gamma^\nu(\not p'-\not k+m)\gamma^\mu(\not p-\not k+m)\gamma_\nu$ で奇の $\ell$ 項は消え、$\ell^2$ 項は $\gamma^\mu$ に比例する。残りを殻方程式と上の恒等式で整理すると、Pauli 項の積分は
$$
F_2(0)=\frac{\alpha}{2\pi}\int_0^1\dd z\,
\frac{2m^2z^2(1-z)}{m^2 z^2+\lambda_\gamma^2(1-z)},\qquad z=x+y.
\label{eq-qft-schwinger-moment}
$$
単体の $\int_0^z\dd x$ が分子の $z$ 一つを供給し、残りの係数は $n=3$ の :eqref[eq-qft-loop-master] から得る。被積分関数は $2(1-z)$ 以下なので、$\lambda_\gamma\to0$ で積分は一になり $a=\alpha/(2\pi)$。$\alpha\simeq1/137$ では約 $1.16\times10^{-3}$ である。

### 頂点分子から Pauli 項を分離する

単体の二変数を $x,y$、$z=x+y$ とする。移動後の二つの分子の定数部分を $A=(1-x)p'-yp$、$B=(1-y)p-xp'$ と置く。Clifford 関係で、一ループ分子の $\ell$ に依らない部分は
$$
\gamma^\nu(\not A+m)\gamma^\mu(\not B+m)\gamma_\nu
=-2\not B\gamma^\mu\not A+4m(A+B)^\mu-2m^2\gamma^\mu
$$
となる。四次元に置いたのは、以下で取り出す Pauli 項の積分が紫外有限だからである。例えば $\bar u'\not p\gamma^\mu\not p'u$ ではまず最初の $\not p$ と $\gamma^\mu$ を反交換し、次に $\not p\not p'$ を入れ替えると
$$
\bar u'\not p\gamma^\mu\not p'u
=\bar u'\bigl[2m(p+p')^\mu-(2p\cdot p'+m^2)\gamma^\mu\bigr]u .
$$
同様に $\not p\gamma^\mu\not p$ と $\not p'\gamma^\mu\not p'$ からはそれぞれ $2m p^\mu-m^2\gamma^\mu$ と $2m p'^\mu-m^2\gamma^\mu$ が出る。$\not p'\gamma^\mu\not p$ は $m^2\gamma^\mu$ だけである。従って $-2\not B\gamma^\mu\not A$ の二つの運動量係数の和は $-4m[2-3z+z^2]$、$4m(A+B)^\mu$ の和は $4m(2-2z)$。二成分の平均が $(p+p')^\mu$ の係数なので $2mz(1-z)$ が残る。差に比例する $q^\mu$ 項は $x\leftrightarrow y$ で奇なので積分で消える。$\ell$ の二次項は回転平均で $\gamma^\mu$ に比例するため、この抽出に寄与しない。

ループ因子は電荷を除いた頂点で $-e^2$、パラメータの $2!$ と Euclid の $I_3=1/(32\pi^2\Delta)$ を合わせると $-e^2/(16\pi^2\Delta)$。殻上の Gordon の式で $(p+p')^\mu=2m\gamma^\mu-i\sigma^{\mu\nu}q_\nu$ の形に移すと、$F_2$ への変換は $-2m$ 倍である。$q^2=0$ では $x$ の範囲の長さが $z$ なので、これらを掛けて :eqref[eq-qft-schwinger-moment] を得る。

### Pauli 形状因子の最初の運動量補正

上で求めた分子は $q^2$ を残しても $F_2$ の射影で同じである。$\lambda_\gamma=0$、$z=x+y$ とし、$x=zt,y=z(1-t)$ と置くと Jacobian は $z$。分母 $z^2[m^2-q^2t(1-t)]$ が分子の $z^2$ を消すため、
$$
F_2(q^2)=\frac{\alpha}{2\pi}\int_0^1\frac{\dd t}{1-q^2t(1-t)/m^2-i0}
$$
となる。$\int_0^1 2(1-z)\dd z=1$ を使った。従って $F_2'(0)=\alpha/(12\pi m^2)$。小さい空間的な運動量移行 $q^2=-m^2/10$ では $F_2(0)$ からの第一補正は $-1/60$ の割合である。電荷形状因子 $F_1$ の微分には赤外に対する別の扱いが必要である。

### 空間的な移行での Pauli 形状因子を積分する

$Q^2=-q^2>0$、$w=Q/(2m)$ とする。前の積分で $y=2t-1$ と置くと $t(1-t)=(1-y^2)/4$。偶関数なので $F_2(-Q^2)/F_2(0)=\int_0^1\dd y/(1+w^2-w^2 y^2)$。ここで $y=w^{-1}\sqrt{1+w^2}\tanh u$ とすれば、分母に現れる $1-\tanh^2u$ と Jacobian が相殺する。従って
$$
\frac{F_2(-Q^2)}{F_2(0)}
=\frac{\operatorname{arctanh}(w/\sqrt{1+w^2})}{w\sqrt{1+w^2}}
=\frac{\operatorname{arsinh}w}{w\sqrt{1+w^2}}.
$$
最後は $\tanh(\operatorname{arsinh}w)=w/\sqrt{1+w^2}$ による。$w\to0$ では $\operatorname{arsinh}w=w-w^3/6+\cdots$、$(1+w^2)^{-1/2}=1-w^2/2+\cdots$ により比は $1-2w^2/3+\cdots=1-Q^2/(6m^2)+\cdots$ となる。$w\gg1$ では比は $\log(2w)/w^2$。$1/Q^2$ による減衰には端点 $t=0,1$ の対数が掛かる。

:::example[四倍の質量二乗に等しい空間的移行]
$Q=2m$ では $w=1$ で比は $\log(1+\sqrt2)/\sqrt2\simeq0.623$。低運動量での一次 Taylor 近似なら $1-Q^2/(6m^2)=1/3$ である。これは小さな運動量に対する勾配の計算が正しくても、$Q/m$ が小さくないと切り詰めた近似の誤差が大きい例である。
:::

### 静磁場との比較で $F_2(0)$ を読む

外部のゆっくり変わる磁場に対する散乱を Breit 枠で考える。移行は $q=(0,\vb q)$、平均空間運動量は零である。Gordon の関係により電流は
$$
\bar u'\left[
\frac{(p'+p)^\mu}{2m}F_1+
\frac{i\sigma^{\mu\nu}q_\nu}{2m}(F_1+F_2)
\right]u .
$$
空間成分では第一項が零である。零運動量近傍で上下成分を静止スピノルに置くと、$\sigma^{ij}$ の上のブロックは $\varepsilon^{ijk}\sigma^k$、外脚規格化を除くと $F_1+F_2$ がスピンと $\vb B$ の結合に掛かる。従って $-\boldsymbol\mu\cdot\vb B$ の $\boldsymbol\mu$ は $q_{\rm el}(F_1+F_2)\boldsymbol S/m$ となる。ここで $q_{\rm el}$ は電子の符号付き電荷、$\boldsymbol S=\boldsymbol\sigma/2$。電荷の規格化を $F_1(0)=1$ として、自由 Dirac 理論の $g=2$ に $2F_2(0)$ が加わる。異常磁気モーメントの符号は $g$ の増減を示し、電子の電荷そのものの符号とは別に定義される。

### 閾値以下の正の時間的移行

$0<q^2<4m^2$、$z=\sqrt{q^2}/(2m)$ とする。Pauli 形状因子の積分で $y=2t-1$ と置くと $F_2(q^2)/F_2(0)=\int_0^1\dd y/[1-z^2+z^2y^2]$。従って
$$
\frac{F_2(q^2)}{F_2(0)}
=\frac{\arctan(z/\sqrt{1-z^2})}{z\sqrt{1-z^2}}
=\frac{\arcsin z}{z\sqrt{1-z^2}} .
$$
例えば $q^2=m^2$ では $\arcsin(1/2)=\pi/6$ から $2\pi/(3\sqrt3)$。$z\uparrow1$ では一次の摂動項が $(1-z^2)^{-1/2}$ の形で増大する。この領域では対を作る運動学が近いため、固定次数の式の適用範囲を別に検討する。

### Dirac 形状因子の減算と赤外発散

Pauli 項と同じ正則化で、電荷に掛かる $F_1(q^2)$ も求められる。二つの電子パラメータを $x,y$、$z=x+y$ とし、$\Delta_q=m^2z^2-xyq^2+(1-z)\lambda_\gamma^2-i0$。Gordon の式を使った後の定数分子で $\gamma^\mu$ に掛かる量は
$$
N_q=m^2(4-4z-2z^2)-2q^2(1-z+xy).
$$
実際、前に整理した $-2\not B\gamma^\mu\not A$ の殻上の $\gamma^\mu$ 成分は $m^2(6-8z+2z^2)-2q^2(1-z+xy)$。ここで $2p\cdot p'=2m^2-q^2$ を使った。さらに $-2m^2\gamma^\mu$ と、運動量型の $2mz(1-z)(p+p')^\mu$ から来る $4m^2z(1-z)\gamma^\mu$ を足す。従って残るのは上の $N_q$ である。$\ell$ の二次の項は回転平均で $(d-2)^2\ell^2\gamma^\mu/d$。この係数は、二か所のガンマ縮約がそれぞれ $(2-d)$ を与えることから従う。Wick 回転で $\ell^2\mapsto-\ell_E^2$ とすると、分子の二次項の積分は $I_2(\Delta_q)-\Delta_q I_3(\Delta_q)$ に比例する。

$F_1(0)=1$ を課すと、$q$ に依存しない対項と外脚の補正は $q=0$ の値と相殺する。$I_2$ の差は $-\log(\Delta_q/\Delta_0)/(16\pi^2)$、$\Delta I_3$ の差は四次元で零。これらに共通のループ係数を掛けて
$$
F_1(q^2)-1=-\frac{\alpha}{4\pi}\int_0^1\dd z\int_0^z\dd x
  \left[\frac{N_q}{\Delta_q}-\frac{N_0}{\Delta_0}
    +2\log\frac{\Delta_q}{\Delta_0}\right]+O(\alpha^2),
\qquad y=z-x .
$$
この差では紫外極が消える。$\lambda_\gamma>0$ を残して $q^2$ で微分してから零に置くと、括弧の微分は $-2(1-z+2xy)/\Delta_0+N_0xy/\Delta_0^2$。$\int_0^z x(z-x)\dd x=z^3/6$ なので
$$
F_1'(0)=-\frac{\alpha}{4\pi}\int_0^1\dd z\,
\left[-\frac{2z(1-z)}{\Delta_0}-\frac{2z^3}{3\Delta_0}
+\frac{m^2(4-4z-2z^2)z^3}{6\Delta_0^2}\right]+O(\alpha^2).
$$
$w=\lambda_\gamma/m$ とすると、$m^2$ を掛けた被積分関数の端点で発散する二項は $-2z/(z^2+w^2)+(2/3)z^3/(z^2+w^2)^2$ と近似できる。差は積分極限を取ることができ、$w=0$ では $4/3-z$。これを確かめるには元の $\Delta_0/m^2=z^2+(1-z)w^2$ と $z^2+w^2$ の差が $-zw^2$ であることを使う。例えば一つ目の分母の置換の誤差は $O(z^2w^2/(z^2+w^2)^2)$ で、$z$ 積分すると $O(w)$。二つ目の分母の平方でも同様に端点で積分可能である。近似した二項の積分は $-(4/3)\log(1/w)-1/3+o(1)$、差の積分は $5/6$。従って同じ調整質量の規約で、
$$
F_1'(0)=\frac{\alpha}{3\pi m^2}
\left[\log\frac{m}{\lambda_\gamma}-\frac38+o(1)\right]+O(\alpha^2).
$$
$F_1(0)$ の電荷の条件は赤外極限でも固定されるが、$F_1'(0)$ には対数が残る。この傾きを自由電子の測定可能な有限半径とみなすことはできない。例えば散乱の測定では軟放射を加え、束縛状態では低エネルギーの束縛準位の寄与を同じ規約へ照合する。

:::example[外部静電場に入る傾きの組]
静電場の Breit 枠で $p'=(E,\vb q/2)$、$p=(E,-\vb q/2)$。同じ二成分スピンに対し $\bar u'\gamma^0u=2m$、$\bar u'u=2E$ である。後者の違いは上下成分の内積の符号から出る。Gordon を使うと $\bar u' i\sigma^{0\nu}q_\nu u/(2m)=2m-2E^2/m=q^2/(2m)$。従って電荷を $2m$ で割った電気的な組は $F_1(q^2)+q^2F_2(q^2)/(4m^2)$。傾きに上の $F_1'(0)$ と $F_2(0)=\alpha/(2\pi)$ を入れると、光子質量での定数 $-3/8$ に対応する項と $F_2(0)/(4m^2)$ が相殺する。外脚を非相対論的規格化にする $m/E$ の因子は、これとは別に木の Darwin 項を作る。
:::

## 真空偏極と電荷のくりこみ

真空偏極の二点ループには閉フェルミオン線の負号がある。運動量で縮約すると $\not q=(\not k+\not q-m)-(\not k-m)$ から二つの一分母積分の差となり、並進を保つ正則化でシフトして相殺する。従ってテンソルは横方向だけである。Euclid の交換に対し荷電の増加を正とするスカラー $\Pi_R(-Q^2)$ を取り、交換振幅の補正を $1+\Pi_R+\cdots$ と定義する。トレースの $q^\mu q^\nu$ の係数は Feynman パラメータ化した分子で $-8x(1-x)$ に比例する。:eqref[eq-qft-loop-master] の $n=2$ を使い、零運動量の値を引くと
$$
\Pi_R(-Q^2)=\frac{2\alpha}{\pi}\int_0^1\dd x\,x(1-x)
\log\left[1+\frac{Q^2}{m^2}x(1-x)\right].
\label{eq-qft-vacuum-polarization}
$$
$Q^2\ll m^2$ では $\int_0^1x^2(1-x)^2\dd x=1/30$ から $\alpha Q^2/(15\pi m^2)$。$Q^2\gg m^2$ では $\int x(1-x)=1/6$、$\int x(1-x)\log[x(1-x)]=-5/18$ から $\alpha[\log(Q^2/m^2)-5/3]/(3\pi)$ となる。例えば $Q/m$ を大きい領域で二倍にすると、この項は $\alpha\log4/(3\pi)$ だけ増す。

### 時間的閾値と吸収部分

:eqref[eq-qft-vacuum-polarization] を $Q^2\mapsto-s-i0$ へ続ける。$s>4m^2$ で $\beta=\sqrt{1-4m^2/s}$ とすると、対数の虚部は $x=(1-\beta)/2$ から $(1+\beta)/2$ まで $-\pi$ である。この区間で $\int x(1-x)\dd x=\beta(3-\beta^2)/12$。従ってこの章の交換補正の符号規約では
$$
\Im\Pi_R(s+i0)=-\frac{\alpha}{3}\beta\left(1+\frac{2m^2}{s}\right).
$$
係数の負号は交換の補正を Euclid で正に取ったことによる。スペクトルの正の重みは $-\Im\Pi_R$ である。例えば $\beta=1/2$ では $-\Im\Pi_R=11\alpha/48$。

### Uehling 補正をスペクトルから組み立てる

一ループ近似において零運動量で引いた交換補正に一回の減算分散関係を使う。複素 $q^2$ の切断が $4m^2$ から始まり、大円で一回引いた被積分関数が消えることは対数のパラメータ積分から従う。上で求めた吸収部分を代入すると
$$
\frac{\Pi_R(-Q^2)}{Q^2}
=\frac{\alpha}{3\pi}\int_{4m^2}^\infty \frac{\dd s}{s(s+Q^2)}
\sqrt{1-\frac{4m^2}s}\left(1+\frac{2m^2}s\right).
$$
例えば $Q^2=0$ の右辺で $y=4m^2/s$ に置換すると、$1/(4m^2)$ と $\int_0^1(1+y/2)\sqrt{1-y}\dd y=4/5$ が出る。これで $\alpha/(15\pi m^2)$ の低運動量係数を再現する。

重い核に対するポテンシャルの補正は $-4\pi Z\alpha\Pi_R(-Q^2)/Q^2$ の空間逆変換である。各 $s$ について $(s+Q^2)^{-1}$ の逆変換は $e^{-\sqrt s r}/(4\pi r)$。従って
$$
\delta V(r)=-\frac{Z\alpha}{r}\frac{\alpha}{3\pi}
\int_{4m^2}^\infty\frac{\dd s}{s}\sqrt{1-\frac{4m^2}s}
\left(1+\frac{2m^2}s\right)e^{-\sqrt{s}r}
\qquad (r>0).
$$
積分の重みは非負なので補正は引力を強める。$r\gg1/m$ では $\sqrt s=2m$ の近傍が支配する。$u=\sqrt s-2m$ と置くと、閾値の重みは $\sqrt{u/m}$ に比例し、指数は $e^{-2mr}e^{-ur}$。従って $\delta V(r)$ は指数減衰 $e^{-2mr}$ にさらに $r^{-5/2}$ を掛けた形となる。短距離の局所近似は、この有限な幅を低運動量の行列要素で接触項へ縮約したものである。

### 横性と光子二点関数の規格化

真空偏極の被積分分子は、四次元のスピノルトレースの規格化を $4$ として
$$
4[k^\mu(k+q)^\nu+k^\nu(k+q)^\mu-\eta^{\mu\nu}(k\cdot(k+q)-m^2)] .
$$
$\ell=k+(1-x)q$ とし、$\ell$ に奇な項を積分で落とす。$q^\mu q^\nu$ の係数は $-8x(1-x)$。残りの $\eta^{\mu\nu}$ の部分は $4[(2/d-1)\ell^2+x(1-x)q^2+m^2]$ である。Wick 回転で $\ell^2\mapsto-\ell_E^2$ とし、全微分の関係で $\int_\ell \ell_E^2/(\ell_E^2+\Delta)^2=d I_1/2$、$(1-d/2) I_1=\Delta I_2$ を使う。従って $\eta^{\mu\nu}$ の積分係数は $4[-\Delta+x(1-x)q^2+m^2]I_2=8x(1-x)q^2 I_2$。これで $q^2\eta^{\mu\nu}-q^\mu q^\nu$ の組が実際に得られる。途中で $d=4$ に固定して二次発散と対数発散を別に扱う手順では、この関係が自動的には成り立たない。

横方向の交換の補正を反復すると、幾何級数で $1/[1-\Pi_R(-Q^2)]$ が現れる。一ループ近似でこの分母を使う場合、厳密に保証される係数は一次 $1+\Pi_R$。例えば $Q\ll m$ では $\Pi_R=O(Q^2/m^2)$ なので、$1/Q^2$ の交換極の留数は減算で固定されたままである。高エネルギー側で $\alpha\log(Q^2/m^2)$ が一に近づけば、一ループの挿入を一回だけ使う近似では足りない。尺度の変化を組織する方法は [くりこみ群](/physics/quantum-field-theory/renormalization-group/) を参照する。

### 大距離での Uehling 補正の係数

真空偏極の分散積分で $w=\sqrt{s}/(2m)\ge1$ とする。すると $\dd s/s=2\dd w/w$、$\sqrt{1-4m^2/s}=\sqrt{w^2-1}/w$ なので
$$
\delta V(r)=-\frac{2Z\alpha^2}{3\pi r}
\int_1^\infty\dd w\, e^{-2mr w}
\left(1+\frac1{2w^2}\right)\frac{\sqrt{w^2-1}}{w^2}.
$$
$r\gg m^{-1}$ では $w=1+t$ と置き、重みの先頭は $3\sqrt{2t}/2$。$t$ が固定値より大きい区間は指数で抑えられ、小さい区間の次の項は $O(t^{3/2})$ になる。従って $\int_0^\infty e^{-2mr t}\sqrt{t}\,\dd t=\Gamma(3/2)/(2mr)^{3/2}$ から
$$
\delta V(r)=-\frac{Z\alpha^2}{4\sqrt{\pi}\,r}\,
\frac{e^{-2mr}}{(mr)^{3/2}}\,[1+O((mr)^{-1})].
$$
係数の質量次元は $1/r$ の一だけである。例えば $mr$ を十分大きい値からもう一つ増やすと、支配的な指数は $e^{-2}$ 倍になり、その上に $(r/(r+1/m))^{5/2}$ が掛かる。電子対の閾値が補正の距離の尺度を決めており、光子の質量が変わったという解釈は必要ない。

### 電荷の空間的な変化率を直接微分する

真空偏極で $X=Q^2/m^2$ と置く。測定される交換の一次補正の勾配は
$$
\frac{\dd\Pi_R(-Q^2)}{\dd\log Q}
=\frac{4\alpha}{\pi}\int_0^1\dd x\,
 x(1-x)\frac{Xx(1-x)}{1+Xx(1-x)} .
$$
被積分関数の分数は $0$ と $1$ の間であり、$X$ とともに単調に増える。従って勾配は非負で $2\alpha/(3\pi)$ を超えない。優収束により $Q/m\to\infty$ では上限に近づく。一方 $Q\ll m$ での先頭は $4\alpha X/(30\pi)=2\alpha X/(15\pi)$。低い運動量では厚みのある仮想対を解像できず、電荷の変化は質量の逆二乗で抑えられる。これは零運動量で電荷を定義する物理的な減算の結論である。[くりこみ群](/physics/quantum-field-theory/renormalization-group/) では、極だけを引いた方式でこの滑らかな変化をどのように照合するかを扱う。

:::example[電荷の大きさとスピンの重みを分ける]
質量 $m_i$、符号付き電荷 $q_i$ を持つ Dirac 場が複数なら、閉ループは各種について別々に足される。共通の電子の $\alpha$ に対して係数は $(q_i/e)^2$ 倍。例えば $m_1\ll Q\ll m_2$、$q_1=e$、$q_2=2e$ なら、第1種の勾配は $2\alpha/(3\pi)$ に近く、第2種は $8\alpha Q^2/(15\pi m_2^2)$ から始まる。両方の質量より上の運動量では勾配は一種の五倍になる。電荷の符号は頂点が二つなので相殺するが、質量の比は分母に残る。
:::

## 電子の自己エネルギーと質量くりこみ

電子の伝播関数の逆を $\not p-m_0-\Sigma(\not p)$ と書く。Feynman ゲージで電子分母に $x$ を割り当てると、移動後の分子は $\gamma^\mu[(1-x)\not p-\not \ell+m]\gamma_\mu=(2-d)[(1-x)\not p-\not \ell]+dm$ である。従って
$$
\Sigma(\not p)=\frac{e^2\mu^{2\epsilon}\Gamma(\epsilon)}{(4\pi)^{d/2}}
\int_0^1\dd x\,
\frac{(2-d)(1-x)\not p+dm}{[xm^2-x(1-x)p^2+(1-x)\lambda_\gamma^2-i0]^\epsilon}
$$
である。極の係数は $\alpha(-\not p+4m)/(4\pi\bar\epsilon)$。殻上では $3\alpha m/(4\pi\bar\epsilon)$ になる。物理質量は逆伝播関数の零点で、$m_{\rm pole}=m_0+\Sigma(\not p=m_{\rm pole})$ を摂動的に解く。対項を用いて質量を $m$ に固定するなら、補正した逆伝播関数が $\not p=m$ で消える条件を課す。留数の規格化はその点で $\not p$ に関する一次係数を調整する。質量のシフトは一ループで $\lambda_\gamma\to0$ が有限だが、一次係数には赤外対数が残り得る。

### 殻上の有限項

前節の自己エネルギーの式を使い、$\lambda_\gamma=0$、$\not p=m$ を外スピノルに作用させる。分母は $x^2m^2$、分子を $m$ で割ると $2+2x-2\epsilon x$ である。:eqref[eq-qft-loop-master] の係数と合わせると
$$
\frac{\Sigma(m)}{m}=\frac{\alpha}{4\pi}
\left[\frac3{\bar\epsilon}+3\log\frac{\mu^2}{m^2}+4\right]+\mathcal O(\epsilon).
$$
有限定数のうち $-1$ は $-2\epsilon x$ と極の積、$+5$ は $-\int_0^1(2+2x)2\log x\,\dd x$ に由来する。例えば $\mu=m$ とした $\overline{\rm MS}$ 質量 $\overline m$ では、ここでの一ループの極質量との差は $m_{\rm pole}/\overline m-1=\alpha/\pi$ である。分子の次元を最初から四に固定すると、有限項が $4$ から $5$ にずれてしまう。

### 逆伝播関数の傾きと赤外端点

自己エネルギーを $\Sigma(\not p)=A(p^2)\not p+B(p^2)m$ と分ける。運動量方向を固定して $\not p$ の正質量枝の固有値を $w=\sqrt{p^2}$ とすれば、殻での微分は $A(m^2)+2m^2(A'(m^2)+B'(m^2))$。Feynman パラメータの分母を微分すると $x(1-x)/\Delta$ が出る。$\lambda_\gamma=0$、$p^2=m^2$ では $\Delta=x^2m^2$ なので、$x\to0$ で $1/(xm^2)$ 型となる。分子はこの端で非零なので積分は対数的に発散する。一方、質量シフトそのものの端点は $\int_0^1\log x\,\dd x$ で有限である。この差によって同じ自己エネルギーから求める質量と留数で赤外極限が異なる。

### 光子に小質量を入れたときの赤外対数

$\lambda_\gamma>0$ を残し、外部運動量の方向を固定して $\not p$ の正の固有値 $w$ で微分する。分子の $w$ に陽に依存する項は一回微分すると $-2(1-x)$ であり、$w=m$ での $\log[x^2m^2+(1-x)\lambda_\gamma^2]$ の積分は $\lambda_\gamma\to0$ で有限である。赤外対数を作るのは分母を微分した部分で、$\alpha/(4\pi)$ を外せば
$$
\int_0^1
\frac{(2+2x)m\cdot 2mx(1-x)}{x^2m^2+(1-x)\lambda_\gamma^2}\,\dd x .
$$
小さい $x$ で分子を $4m^2 x$ に、分母を $x^2m^2+\lambda_\gamma^2$ に置き換えた差は赤外で可積分である。任意の固定した小さい $x_0$ まで近似核を積分すると $2\log[(x_0^2m^2+\lambda_\gamma^2)/\lambda_\gamma^2]$。従って傾きの赤外対数は $\alpha\log(m/\lambda_\gamma)/\pi$。質量の位置にはこの対数がなく、質量殻における場の留数をそのまま赤外安全な測定量として扱うことはできない。

### 共変ゲージを動かしたときの質量の位置

質量ゼロの光子に対する共変ゲージの差は、二点関数の被積分式で $k_\mu k_\nu$ に比例する。正則化を保ち、$Q=\not p-m$、$R=\not p-\not k-m$ とすれば $\not k=Q-R$。従って $\not k\,R^{-1}\not k=QR^{-1}Q-2Q+R$ である。殻上の外スピノルに挟むと $Q$ を含む項が消え、残る $\bar u R u=-\bar u\not k u$ は $k$ について奇。縦の分母 $(k^2+i0)^2$ に対し対称な正則化の積分は零になる。従って一ループの殻上の質量シフトはゲージパラメータに依存しない。一方、逆伝播関数の一次の傾きを取ると $Q$ に比例する項を含むため、同じ議論で留数を固定することはできない。この差は外部の電荷を測る物理的な頂点と、場の規格化で決めた伝播関数の留数を区別する理由の一つである。

## 赤外発散と軟光子放射（Bloch–Nordsieck の定理）

出入りする電荷つき外脚の添字を $i$、出射・入射の符号を $\eta_i=+1,-1$ とする。軟光子の四元運動量 $k$ について、隣の電子分母を $(p_i+\eta_i k)^2-m_i^2\simeq2\eta_i p_i\cdot k$ とし、分子を殻方程式で簡約すると
$$
\mathcal M_{n+\gamma}\simeq J_{\rm soft}\cdot\epsilon^*\,\mathcal M_n,\qquad
J_{\rm soft}^\mu=\sum_i\eta_i q_i\frac{p_i^\mu}{p_i\cdot k}.
$$
電荷保存から $k\cdot J_{\rm soft}=0$。従って実光子の補正は $\int\dd\Pi_k\,[-J_{\rm soft}^2]$ で、$J\propto1/\omega$ と $\dd\Pi_k\propto\omega\dd\omega$ から $\log(\Delta E/\lambda_\gamma)$ を生む。電子質量を有限に保つと角度積分は有限である。

同じ軟領域の仮想補正の実部を考える。先頭軟近似では硬い電流を指定した外部源として光子の各振動子を線形駆動する。駆動の時間順序積は交換子が数なので、正規順序で位相と $\exp[-\frac12\int |J\cdot\epsilon|^2\dd\Pi_k]$ を伴う。これは一振動子で $e^{fa^\dagger-f^*a}=e^{fa^\dagger}e^{-f^*a}e^{-|f|^2/2}$ を適用したものである。零光子確率の一次補正は実放射と逆符号になる。従って分解能 $\Delta E$ 未満の放射を全て足すと $\lambda_\gamma$ の対数が相殺し、分解能より上の有限な対数が残る。これが質量を持つ荷電粒子に対する :term[Bloch–Nordsieck の定理]{en="Bloch–Nordsieck theorem" yomi="ぶろっほのるつじーくのていり"} の軟近似の内容である。

### 速度を変えた電荷からの放射

静止していた電荷が散乱後に速度 $v\vb e_z$ を持つとし、他の放射源の寄与を無視できる外部源模型を考える。軟電流の共通因子 $q/\omega$ を取り出すと $J=(1,0,0,v)/(1-vx)-(1,0,0,0)$、$x=\cos\theta$。従って $-J^2=v^2(1-x^2)/(1-vx)^2$。角度平均は
$$
\frac12\int_{-1}^1\frac{v^2(1-x^2)}{(1-vx)^2}\dd x
=\frac1v\log\frac{1+v}{1-v}-2 .
$$
$u=1-vx$ とすると被積分関数は $-1+2/u+(v^2-1)/u^2$ なので積分できる。$v\ll1$ では $2v^2/3+\order{v^4}$、$v=1/2$ では $2\log3-2$。$v\to1$ の対数は前方に集中するため、固定した質量での軟極限と質量ゼロ極限を区別する。

### 実光子を分解能の内外に分ける

独立な光子モードの変位演算子を $D(f)=\exp(fa^\dagger-f^*a)$ とする。$a|0\rangle=0$ と正規順序化で $\langle n|D(f)|0\rangle=e^{-|f|^2/2}f^n/\sqrt{n!}$。よって放射個数の確率は Poisson 分布であり、全ての $n$ を足すと一である。この議論を軟領域の各モードへ適用すると、各光子のエネルギーが独立に $\Delta E$ 未満なら許す測定では、その領域の確率の和が一になる。$\Delta E$ より上の軟領域で光子を一個も許さない確率は
$$
\exp\left[-\int_{\Delta E<\omega<E_h}\dd\Pi_k\,(-J_{\rm soft}^2)\right].
$$
分解能が各光子のエネルギーでなく全放射エネルギーに課される場合は、光子数を足すときに $\sum_i\omega_i<\Delta E$ を付ける。一次では両者が同じ制約になる。いずれの場合も赤外で無限に小さい光子を区別しないため、最低次の対数の相殺は保たれる。

:::example[角積分と分解能の係数]
速度が $0$ から $v$ へ変わる一電荷の模型では、軟光子の測度の積が $\dd\omega/\omega$ に掛ける係数は $q^2/(4\pi^2)$。実際、$\dd\Pi_k=\omega\dd\omega\,\dd\Omega/[2(2\pi)^3]$ で、$J_{\rm soft}$ の $q/\omega$ を二乗し、立体角の $4\pi$ を掛ける。従って $q^2=4\pi\alpha$ なら、$v=1/2$ で分解能を半分にしたときの先頭の対数補正の差は $-(\alpha/\pi)(2\log3-2)\log2$。入射と出射で速度が同じならベクトル $J_{\rm soft}$ 自体が零になるので、この補正も消える。
:::

### 二つの移動電荷に共通する相対速度

同じ電荷の入射と出射の四速度を $u,u'$ とし、$u^2=u'^2=1$、時間方向が未来向きとする。軟電流の $q$ を除く部分は $u'/(k\cdot u')-u/(k\cdot u)$。軟い対数の係数は Lorentz 不変である。実際、$\dd^3k/(2\omega)$ は不変、$J^2$ はスカラーで、光子の尺度を変えると $J^2$ はその逆二乗。異なる枠で同じエネルギーカットを使う差は、角度依存の有限な倍率を積分の上下端に掛ける効果になる。従って $\log\lambda_\gamma$ の係数は変わらない。$u$ の静止系で $u\cdot u'=w$、相対速度 $v=\sqrt{1-w^{-2}}$ と置くと、前に求めた角係数は
$$
C(w)=\frac{2w\,\operatorname{arcosh}w}{\sqrt{w^2-1}}-2 .
$$
これは $v=\tanh(\operatorname{arcosh}w)$ を入れて得る。弾性の移行 $t=(p'-p)^2<0$ なら $w=1-t/(2m^2)$ だから、$-t\ll m^2$ で $C(w)=-2t/(3m^2)+O(t^2/m^4)$。例えば電荷が電子なら確率の軟対数は $(\alpha/\pi) C(w)$ を含む。木と一ループの頂点振幅の干渉から来る $2F_1'(0)t$ の仮想赤外対数は $(2\alpha/(3\pi m^2))t\log(m/\lambda_\gamma)$。実放射と逆符号であり、一つの電荷が曲げられる場合の小移行での係数が一致する。

### 全放射エネルギーに対する軟光子の和

硬い反応を固定し、角度を積分した先頭の軟放射の重みを $A\,\dd\omega/\omega$ と書く。ここで $A\ge0$ は電荷と外脚速度から決まる。$\lambda_\gamma<\omega<E_h$ の範囲で、仮想の零光子因子を含めた実放射のエネルギー和を $E_{\rm rad}$ とする。$\nu>0$ の Laplace 変換は、光子数を全て足して
$$
\left\langle e^{-\nu E_{\rm rad}}\right\rangle
=\exp\left[A\int_{\lambda_\gamma}^{E_h}
  \frac{\dd\omega}{\omega}(e^{-\nu\omega}-1)\right].
$$
$n$ 光子の積分に $1/n!$ が掛かり、Laplace の重みは各 $\omega$ に分離する。零光子確率は $e^{-A\log(E_h/\lambda_\gamma)}$ なので、指数化した実放射から引く $-1$ が現れる。積分の原点で $(e^{-\nu\omega}-1)/\omega\to-\nu$。従って調整質量を除いた変換は有限であり、一個ごとの上限でなく全エネルギーを測る場合でも相殺を直接示せる。

$E_{\rm rad}<\Delta<E_h$ とする。$\lambda_\gamma\to0$ を取った指数の積分は
$$
-\log(\nu E_h)-\gamma_E-\int_{\nu E_h}^{\infty}\frac{e^{-u}}{u}\dd u.
$$
$\nu E_h$ の微分で両辺の微分が同じになり、原点の定数は $\Gamma(\epsilon)$ の展開 $\Gamma(\epsilon)=1/\epsilon-\gamma_E+\cdots$ から得る。累積確率の Laplace 変換は上の期待値を $\nu$ で割ったもの。最後の積分は $\int_{E_h}^\infty e^{-\nu\omega}\dd\omega/\omega$ なので、これを少なくとも一回含む指数展開の項は、逆変換では $E_{\rm rad}\ge E_h$ にだけ台を持つ。従って $\Delta<E_h$ では $\nu^{-1-A}$ の逆変換だけが寄与し、
$$
P(E_{\rm rad}<\Delta)
=\frac{e^{-\gamma_E A}}{\Gamma(1+A)}
\left(\frac{\Delta}{E_h}\right)^A .
$$
この式は独立な先頭軟放射と上の鋭い上限の模型である。$\log\Gamma(1+A)=-\gamma_E A+O(A^2)$ なので、一次の補正は $A\log(\Delta/E_h)$。各光子に別々の上限を課した式と一次では同じだが、多光子で全エネルギーを制限する効果は $\Gamma(1+A)$ に残る。

## Lamb シフトの計算

### 束縛準位への先頭対数

無限に重い原子核の Coulomb 場 $V=-Z\alpha/r$ にある電子を考える。速度 $Z\alpha\ll1$ とし、低エネルギー側で光子の双極子結合 $-q\vb p\cdot\vb A/m$ に対する二次摂動を取る。中間電子のエネルギーを $E_b$ とすれば分母は $-(\omega+E_b-E_a)$。自由電子の質量補正に対応する $-1/\omega$ を引き、束縛エネルギーより上で展開すると $(E_b-E_a)/\omega^2$ が残る。光子の角平均で横射影は $2\delta_{ij}/3$、$\dd\Pi_k$ の角積分は $\omega\dd\omega/(4\pi^2)$。従って対数の係数は $2\alpha/(3\pi m^2)\sum_b(E_b-E_a)|\mel b{\vb p}a|^2$ である。

完全性により和は $\frac12\sum_i\bra a[p_i,[H,p_i]]\ket a$ となる。$[H,p_i]=i\partial_i V$、$[p_i,i\partial_iV]=\partial_i^2 V$ だから、$\nabla^2V=4\pi Z\alpha\delta^3$ を用いて
$$
\Delta E_{nS}^{\rm LL}=\frac{4\alpha(Z\alpha)^4 m}{3\pi n^3}
\log\frac1{(Z\alpha)^2}
$$
を得る。ここで Coulomb 波動関数の原点値は $|\psi_{nS}(0)|^2=(m Z\alpha)^3/(\pi n^3)$ を使った。低エネルギーでは状態に依存する励起エネルギーが積分の下端を実効的に決め、高エネルギーでは相対論的な頂点と照合して上端を $m$ の程度とする。この上下端の定数倍は非対数項を変える。水素の $2S$ では先頭対数だけでエネルギーは $10^{-6}\,\mathrm{eV}$ の程度で、同じ近似の $2P$ では原点値が零である。精密な差には Bethe 対数、真空偏極、反跳なども必要である。

### 先頭対数と Bethe 対数

非相対論的な低エネルギー積分を $\omega<K$ で止め、束縛エネルギーの尺度 $m(Z\alpha)^2$ と相対論的尺度 $m$ の間に $K$ を取る。双極子近似を一様に使う側ではさらに $K\ll mZ\alpha$ とする。Coulomb の連続スペクトルは無限に広がるため、全ての $\Delta_b$ に対して $\Delta_b\ll K$ と仮定することはできない。質量補正を引いた後の径積分は、$\Delta_b=E_b-E_a>0$ なら
$$
\int_0^K\dd\omega\,\frac{\Delta_b}{\omega+\Delta_b}
=\Delta_b\log\frac{K+\Delta_b}{\Delta_b}.
$$
異なる中間状態で $\log\Delta_b$ が違うため、完全性でまとめられるのは $\log K$ の係数だけである。状態依存の対数重みが Bethe 対数になる。$\Delta_b<0$ の励起準位では極を処方して実部に主値、虚部に放射幅を取る。高エネルギー側の照合による $\log(m/K)$ を加えると任意の分割尺度 $K$ が消える。例えば $\Delta_b$ を全て同じ $\Delta_0$ で近似すると二重交換子から先頭対数が出るが、$2S$ と $2P$ の非対数部分をこの置換で求めることはできない。

### 真空偏極の先頭の局所補正

核の Coulomb ポテンシャルを $\widetilde V(\vb q)=-4\pi Z\alpha/\vb q^2$ とする。$|\vb q|\sim mZ\alpha\ll m$ で :eqref[eq-qft-vacuum-polarization] の先頭項を掛けると $\delta\widetilde V=-4 Z\alpha^2/(15m^2)$。逆 Fourier 変換で $-4Z\alpha^2\delta^3(\vb r)/(15m^2)$ である。従って無限に重い核に対する $nS$ の一次摂動は $-4\alpha(Z\alpha)^4m/(15\pi n^3)$ になる。$nP$ では原点の波動関数が零なのでこの局所項は消える。先頭対数の自己エネルギーとの大きさの比は $1/[5\log(1/(Z\alpha)^2)]$ であり、水素では約 $0.02$ である。符号は短距離で電荷の大きさが増すため引力が強まることと整合する。

### 連続スペクトルを含むエネルギー和

$H=\vb p^2/(2m)+V$ とし、$|a\rangle$ は規格化された束縛固有状態とする。二重交換子を展開すると、
$$
\frac12\langle a|[p_i,[H,p_i]]|a\rangle
=\langle a|p_i(H-E_a)p_i|a\rangle .
$$
従って中間状態の和には束縛固有状態だけでなく、連続スペクトルに対する分解も含まれる。Coulomb 力で $\nabla^2 V=4\pi Z\alpha\delta^3$ となる符号は、$V=-Z\alpha/r$ の径方向微分を小球で積分すれば確認できる。この接触項は $S$ 波で非零である。$P$ 波で消えるのは $\log K$ の係数であり、中間状態のエネルギーが異なるため $\log|\Delta_b|$ で重み付けした和まで消えることは一般に従わない。

:::example[水素における先頭対数の単位]
$Z=1$、$n=2$ とし、$m=0.511\,\mathrm{MeV}$、$\alpha\simeq1/137$ を代入すると、質量因子 $m\alpha^5\simeq1.06\times10^{-5}\,\mathrm{eV}$ に $1/(6\pi)$ と $\log(137^2)\simeq9.84$ が掛かる。従って先頭対数の $2S$ への寄与は約 $5.5\times10^{-6}\,\mathrm{eV}$。周波数に直すには $\Delta E/h$ とし、角周波数 $\Delta E/\hbar$ と区別する。この近似の信頼性は $\alpha$ が小さいことに基づくが、有限な一ループの準位差には状態別の Bethe 対数も要る。
:::

### Bethe 対数を演算子の行列要素として定義する

Coulomb Hamiltonian $H$ の規格化された束縛状態 $|a\rangle=|n\ell m_\ell\rangle$ を取る。無限大の核質量、一電子、非相対論の近似を用いる。$E_B=m(Z\alpha)^2/2$ とし、
$$
\log k_0(n,\ell)=
\frac{n^3}{2m^3(Z\alpha)^4}
\sum_i\left\langle a\left|
p_i(H-E_a)\log\frac{|H-E_a|}{E_B}p_i
\right|a\right\rangle
$$
と定める。対数を含む関数は $x=0$ で $x\log|x|\to0$ と延長し、スペクトル分解で作用させる。これは $K$ に依存しない部分を上の二次摂動から取り出した量である。実際、低エネルギー積分の実部から $\Delta_b\log(K/|\Delta_b|)$ を取り出し、$\log(K/E_B)$ に比例する部分と残りへ分ける。二重交換子により前者の係数は $2m^3(Z\alpha)^4\delta_{\ell0}/n^3$。従って低エネルギー側の対数と状態依存の定数の和は
$$
\Delta E_{\rm low}=
\frac{4\alpha m(Z\alpha)^4}{3\pi n^3}
\left[\delta_{\ell0}\log\frac{K}{E_B}-\log k_0(n,\ell)\right]
$$
となる。右辺は $E_B\ll K\ll mZ\alpha$ の分離における先頭次数を表す。有限な上下端の取り方から出る共通の定数は、同じ正則化で高エネルギー側を計算して決める。$P$ 状態でも $\log k_0(n,1)$ は一般に零でない。完全性で畳み込めるのは $\Delta_b$ で、$\Delta_b\log|\Delta_b|$ を含む和には中間状態の分布が残る。

:::example[分離尺度を二倍にしたとき]
$K$ を $2K$ に変えると $nS$ の低エネルギー項は $4\alpha m(Z\alpha)^4\log2/(3\pi n^3)$ だけ増える。同じ規約での高エネルギー側は $\log(m/K)$ のため同量だけ減る。$\ell>0$ ではこの変化は先頭次数では零である。したがって二つの形の即時的な比較では $\log K$ の係数の打消しまで確かめ、その上で有限な状態依存項を計算する必要がある。
:::

### 放射幅は同じ分母の虚部から出る

励起状態 $|a\rangle$ に $E_b<E_a$ の中間状態がある場合、二次摂動の分母は $E_a-E_b-\omega+i0$。$1/(x+i0)$ の虚部は $-\pi\delta(x)$ なので、準位の虚部は負となり、幅は $-2\Im\Delta E$ である。双極子近似で $\omega_{ab}=E_a-E_b>0$、$\sum_{\lambda}\int\dd\Omega\,|\mel b{\vb p\cdot\boldsymbol\epsilon_\lambda}a|^2=(8\pi/3)|\mel b{\vb p}a|^2$。従って
$$
\Gamma_{a\to b}=\frac{4\alpha}{3m^2}\omega_{ab}|\mel b{\vb p}a|^2
=\frac{4\alpha}{3}\omega_{ab}^3|\mel b{\vb r}a|^2.
$$
最後では $[H,r_i]=-ip_i/m$ を用いた。双極子近似が使えるのは $\omega_{ab}$ と代表半径の積が小さい場合。たとえば $2S\to1S$ で $\mel{1S}{\vb r}{2S}$ は空間反転で零なので、この一光子の寄与は消える。一方、実部の二重交換子は他の中間状態も足すため、$2S$ の先頭対数まで消えるわけではない。エネルギーの実部と幅は同じ処方から得られるが、幅は運動学的に許される特定の中間状態だけを選ぶ。

### Pauli 項が $2P$ へ与える準位補正

静的な Coulomb 外場の電気ポテンシャルエネルギーを $V=-Z\alpha/r$ とする。$F_2(0)=a_e$ の寄与のうちスピン依存部分を低速で取り出す。Gordon の式により、時間成分の Pauli 項は $\bar u'[\gamma^0-(E'+E)/(2m)]u$。外脚の $2m$ を除いた Dirac 密度のスピン項は $i\boldsymbol\sigma\cdot(\vb p'\times\vb p)/(4m^2)$。スカラー密度のスピン項は下成分の符号が逆だからその負である。従って差では $i\boldsymbol\sigma\cdot(\vb p'\times\vb p)/(2m^2)$。Fourier 移行 $\vb q=\vb p'-\vb p$ に対して $i\vb q\,\widetilde V$ を $\nabla V$ に戻すと、作用素は
$$
\delta H_{\rm so}=\frac{a_e}{2m^2}\boldsymbol\sigma\cdot(\nabla V\times\vb p)
=\frac{a_e Z\alpha}{2m^2 r^3}\boldsymbol\sigma\cdot\vb L .
$$
スピンに依らない $q^2$ の成分は接触項なので $P$ 波では消える。水素型の $n=2,\ell=1$ の径関数は $R_{21}(r)=(r/a_B)e^{-r/(2a_B)}/(2\sqrt{6a_B^3})$、$a_B=1/(mZ\alpha)$。これは径方程式に $r e^{-r/(2a_B)}$ を代入し、$\int_0^\infty r^4e^{-r/a_B}\dd r=24a_B^5$ で規格化して得る。従って $\langle r^{-3}\rangle=\int_0^\infty |R_{21}|^2\dd r/r=1/(24a_B^3)$。$\boldsymbol\sigma\cdot\vb L=J^2-L^2-S^2$ の固有値は $j=1/2$ で $-2$、$j=3/2$ で $+1$。従って $2P_{1/2}$ のこの補正は $-\alpha m(Z\alpha)^4/(48\pi)$、$2P_{3/2}$ ではその $-1/2$ 倍である。これに状態依存の Bethe 対数も加わるため、$P$ 波の一ループのエネルギーがこのスピン項だけで決まるわけではない。

### 同じ次元正則化で有限な一ループの部分まで照合する

先頭の $\alpha(Z\alpha)^4m$ では、硬い自己エネルギーの静電場に対する局所演算子は $\nabla^2 V$ とスピン軌道項である。相対論的な振幅と低速の振幅を一光子の外場について比較して係数を決める。複数の外場を含む新しいゲージ不変な項には二つの電場などが必要で、例えば $(\nabla V)^2/m^3$ は $\alpha(Z\alpha)^6m$ から寄与する。運動方程式に比例する時間微分の項は、この次数の殻上での照合から除く。核は無限に重い点電荷とする。

光子を質量ゼロとし、時空の次元を $4-2\epsilon$ とする。共通因子を $C_\epsilon=(4\pi)^\epsilon\Gamma(1+\epsilon)(\mu/m)^{2\epsilon}$ と置く。上の $F_1$ の傾きを再計算する際、$z=0$ だけは次元を四に戻す前に調べる。ループ運動量を零とした分子に左右から外部スピノルを掛けると、$\bar u'\gamma^\nu(\not p'+m)=2p'^\nu\bar u'$ と $(\not p+m)\gamma_\nu u=2p_\nu u$ から $4p\cdot p'\,\bar u'\gamma^\mu u$。従って $z=0$ の $N_q=4m^2-2q^2$ は次元に依らない。$d$ 依存の分子の補正は $O(\epsilon z)$ なので、傾きの端点に極を作らない。

一方、三分母の積分には $\Delta_q^{-1-\epsilon}$ があるため、これを微分すると傾きの $N_0xy$ の部分に $1+\epsilon$ が付く。$x$ を $0$ から $z$ へ積分し $z^{-2\epsilon}$ をくくると、四次元の項は $-4/(3z)+4/3-z$。次元を保ったときの端点には、上の $1+\epsilon$ による $2\epsilon/(3z)$ が加わる。二次のループ分子は傾きを取った後では端点で $O(z)$ なので、その次元依存は有限部分に寄与しない。従って $\int_0^1 z^{-1-2\epsilon}\dd z=-1/(2\epsilon)$ を赤外で収束する側から続けると、$F_1'(0)$ は $-\alpha C_\epsilon/(4\pi m^2)$ に
$$
\left(-\frac43+\frac{2\epsilon}3\right)\left(-\frac1{2\epsilon}\right)
+\int_0^1\left(\frac43-z\right)\dd z+O(\epsilon)
=\frac{2}{3\epsilon}+\frac12+O(\epsilon)
$$
を掛けたものになる。電気的な組の傾き $F_1'(0)+F_2(0)/(4m^2)$ では有限な $-\alpha/(8\pi m^2)$ が打ち消され、$-\alpha C_\epsilon/(6\pi m^2\epsilon)$ が残る。静電場で $q^2=-|\vb q|^2$ なのでこの係数は $\nabla^2 V$ に掛かる。

低エネルギー側は空間次元 $D=3-2\epsilon$ で計算する。質量シフトを引いた双極子積分で、中間状態との差をまず $\Delta=E_b-E_a>0$ とする。一つの $|\mel b{\vb p}a|^2$ に掛かる係数は
$$
\frac{e^2\mu^{2\epsilon}}{m^2}
\frac{D-1}{D}\frac{S_{D-1}}{2(2\pi)^D}
\int_0^\infty \dd k\,k^{D-2}
\left[\frac1k-\frac1{k+\Delta}\right] .
$$
横偏極の角平均が $(D-1)/D$。径積分は $k=\Delta u$ として $\Delta^{1-2\epsilon}\Gamma(1-2\epsilon)\Gamma(2\epsilon)$ であり、$0<\epsilon<1/2$ の領域で Beta 積分として得られる。球面積 $S_{D-1}=2\pi^{D/2}/\Gamma(D/2)$ を入れ、$2\alpha C_\epsilon/(3\pi m^2)$ をくくる。残る角因子は $1+(5/3-2\log2)\epsilon+O(\epsilon^2)$。これは $\Gamma(3/2)=\sqrt\pi/2$ と $\Gamma'(3/2)/\Gamma(3/2)=2-\gamma_E-2\log2$ から出る。後者は Gamma の倍角式を $1/2$ で微分してから漸化式を使えばよい。Gamma 積の極 $1/(2\epsilon)$ とこの因子を掛けると、全体は $\Delta[1/(2\epsilon)+5/6-\log(2\Delta/m)]$ に比例する。励起状態では負の $\Delta$ へ極の処方に従って続け、実部では $\log|\Delta|$ を取る。

二重交換子の恒等式は空間次元を保ったまま $2\sum_b\Delta_b|\mel b{\vb p}a|^2=\langle\nabla^2 V\rangle$ である。従って硬い赤外極と低い側の紫外極は、同じ Coulomb ポテンシャルと同じ次元の行列要素に掛かって相殺する。相殺するまでは両方の係数に $C_\epsilon$ を保つ。次元を三に戻し、Bethe 対数の定義を用いれば一電子の自己エネルギーの実部は
$$
\Delta E_{\rm SE}=
\frac{4\alpha m(Z\alpha)^4}{3\pi n^3}
\left[
\delta_{\ell0}\left(\log\frac1{(Z\alpha)^2}+\frac56\right)
-\log k_0(n,\ell)\right]
+\langle\delta H_{\rm so}\rangle .
$$
ここに真空偏極の接触項を加えると、同じ次数の $nS$ の角括弧の定数は $5/6-1/5=19/30$ となる。$2P_{1/2}$ では $\delta_{\ell0}=0$、スピン軌道項は上で求めた $-\alpha m(Z\alpha)^4/(48\pi)$。有限な準位差を数値にする際は、これらに両状態の Bethe 対数のスペクトル積分を加える。反跳、核の有限半径、高次の輻射補正はこの次数と仮定から外れる。束縛準位の非相対論的な自己エネルギーについては :cite[sakurai-aqm] も参照できる。

## 演習問題

::::exercise{#exr-qft-radiative-corrections-1 level=basic}
等しい二分母 $A=B$ に :eqref[eq-qft-feynman-parameter] を $a=b=1$ で適用せよ。

:::solution
分母は $A^2$ に固定され、$x$ の積分が一なので $1/A^2$ になる。
:::
::::

::::exercise{#exr-qft-radiative-corrections-2 level=basic}
$n=3,d=4$ の Euclid 積分を求めよ。

:::solution
:eqref[eq-qft-loop-master] の比は $\Gamma(1)/\Gamma(3)=1/2$、$\Delta$ の冪は $-1$ なので $1/(32\pi^2\Delta)$ である。
:::
::::

::::exercise{#exr-qft-radiative-corrections-3 level=basic}
$\lambda_\gamma=0$ における :eqref[eq-qft-schwinger-moment] の積分値を求めよ。

:::solution
$z>0$ で $m^2z^2$ が消え、$\int_0^1 2(1-z)\dd z=1$ である。
:::
::::

::::exercise{#exr-qft-radiative-corrections-4 level=basic}
真空偏極の低エネルギー係数に現れる $\int_0^1 x^2(1-x)^2\dd x$ を求めよ。

:::solution
展開して $1/3-2/4+1/5=1/30$ となる。
:::
::::

::::exercise{#exr-qft-radiative-corrections-5 level=advanced}
真空偏極の大運動量展開の定数項を求めよ。

:::solution
対称性から $\int_0^1 x(1-x)\log[x(1-x)]\dd x=2\int_0^1(x-x^2)\log x\,\dd x$。各積分は $-1/4,-1/9$ なので $-5/18$ となる。$\log(Q^2/m^2)$ の係数 $1/6$ の単位では $-5/3$ である。
:::
::::

::::exercise{#exr-qft-radiative-corrections-6 level=advanced}
自己エネルギーの極の分子を積分し殻上で評価せよ。

:::solution
$d\to4$ の極の係数では $\int_0^1[-2(1-x)\not p+4m]\dd x=-\not p+4m$ である。殻上で $\not p=m$ を作用させると $3m$ になる。
:::
::::

::::exercise{#exr-qft-radiative-corrections-7 level=advanced}
分解能 $\Delta E$ と硬い尺度 $E_h$ で軟対数を相殺させよ。

:::solution
共通の角係数を $C$ とすると実放射は $C\log(\Delta E/\lambda_\gamma)$、仮想は $-C\log(E_h/\lambda_\gamma)$。和は $-C\log(E_h/\Delta E)$ である。
:::
::::

::::exercise{#exr-qft-radiative-corrections-8 level=advanced}
先頭対数に限って同じ原子の $1S$ と $2S$ の Lamb シフトの比を求めよ。

:::solution
対数と他の係数は同じで、残りは $n^{-3}$ なので $\Delta E_{1S}/\Delta E_{2S}=2^3=8$ である。
:::
::::

::::exercise{#exr-qft-radiative-corrections-9 level=advanced}
空間的な運動量移行 $Q$ で、既知の一種の電荷に加えて、質量 $M$、電荷が電子の二倍の Dirac 粒子が真空偏極に寄与する。零運動量で差し引いたこの粒子の補正を $Q^4$ まで求めよ。

:::solution
:eqref[eq-qft-vacuum-polarization] で $m\mapsto M$、$\alpha\mapsto4\alpha$ とする。$z=Q^2/M^2$ に対し $\log(1+zx(1-x))=zx(1-x)-z^2x^2(1-x)^2/2+\cdots$。$\int_0^1x^2(1-x)^2\dd x=1/30$、$\int_0^1x^3(1-x)^3\dd x=1/140$ だから
$$
\Delta\Pi_R(-Q^2)=\frac{4\alpha}{15\pi}z-\frac{\alpha}{35\pi}z^2+\order{z^3}.
$$
例えば $Q=M/2$ での第二項と第一項の比は $-3/112$ である。係数が電荷の四乗でなく二乗に比例するのは、二点ループに頂点が二つあるためである。
:::
::::

::::exercise{#exr-qft-radiative-corrections-10 level=advanced}
質量 $m$ の Dirac 場の真空偏極の零運動量減算について、分散表示から空間運動量の展開の $Q^4$ の係数を再導出せよ。

:::solution
分散表示の $1/[s(s+Q^2)]$ を $s^{-2}-Q^2s^{-3}+\cdots$ と展開する。$s\ge4m^2$ なので $Q^2<4m^2$ なら級数が一様に収束する。$Q^4$ の係数は
$$
-\frac{\alpha}{3\pi}\int_{4m^2}^\infty\frac{\dd s}{s^3}
  (1+2m^2/s)\sqrt{1-4m^2/s}
=-\frac{\alpha}{48\pi m^4}\int_0^1 y(1+y/2)\sqrt{1-y}\,\dd y .
$$
ここでは $y=4m^2/s$ を使った。$w=\sqrt{1-y}$ と置くと、最後の積分は $\int_0^1 (3w^2-4w^4+w^6)\dd w=1-4/5+1/7=12/35$。従って $Q^4$ の係数は $-\alpha/(140\pi m^4)$。パラメータ積分で $\log(1+u)$ の二次項を使うと、同じ係数は $-(\alpha/\pi)\int_0^1 x^3(1-x)^3\dd x/m^4$。前者の重みは正の吸収から作られており、展開の次の項が負になる理由も直接分かる。
:::
::::

::::exercise{#exr-qft-radiative-corrections-11 level=advanced}
質量 $m$ の Pauli 形状因子について、空間的な運動量移行が $Q\ll m$ のとき $F_2(-Q^2)/F_2(0)$ の $Q^4$ までの展開を求めよ。一次の減少率と二次の補正の符号を説明せよ。

:::solution
質量ゼロの調整光子を用いた一ループ式は $\int_0^1 [1+(Q^2/m^2)t(1-t)]^{-1}\dd t$。$t(1-t)\le1/4$ なので $Q^2<4m^2$ なら級数を項別に積分できる。一次の係数は $-\int_0^1 t(1-t)\dd t=-1/6$、二次は $\int_0^1 t^2(1-t)^2\dd t=1/30$。従って
$$
\frac{F_2(-Q^2)}{F_2(0)}
=1-\frac{Q^2}{6m^2}+\frac{Q^4}{30m^4}+O(Q^6/m^6).
$$
分母は $Q^2$ とともに増えるので関数そのものは単調に減る。一方その一次の接線との比較では、$1/(1+u)-(1-u)=u^2/(1+u)\ge0$。従って二次の正符号は、減少が一次近似より緩くなることを表す。例えば $Q=m/2$ で二次項は $1/480$、一次の変化は $-1/24$ である。
:::
::::

::::exercise{#exr-qft-radiative-corrections-12 level=advanced}
同じ電荷の入射と出射だけを含む軟放射で、四速度の内積を $w=5/4$ とする。入射静止系での出射速度と、先頭対数の角係数を求めよ。

:::solution
$w=(1-v^2)^{-1/2}$ なので $v=\sqrt{1-16/25}=3/5$。従って角係数は $v^{-1}\log[(1+v)/(1-v)]-2=(5/3)\log4-2$ で、約 $0.31$ である。この係数は速度の向きを指定しなくても、放射を全立体角で積分したときには $w$ だけで決まる。入射と出射が同じ速度なら $w=1$ で角係数は零。実放射と仮想の両方に同じ係数が掛かるので、有限な分解能への変化には残っても調整質量への対数は相殺する。
:::
::::

## 確認テスト

::::quiz
Euclid ループの分子 $\ell_\mu\ell_\nu$ を回転平均した結果を選べ。

- [x] $\delta_{\mu\nu}\ell^2/d$
- [ ] $\delta_{\mu\nu}\ell^2$
- [ ] $\delta_{\mu\nu}\ell^2/(d-1)$

:::explanation
回転で不変な二階テンソルは計量に比例し、両添字を縮約すると $\ell^2$ と $d$ 倍の比例係数が等しくなる。
:::
::::

::::quiz{answer="3" tol="1e-6"}
一ループで殻上の自己エネルギーの極の分子を $m$ の倍数で答えよ。

:::explanation
分子の $x$ 積分が $-\not p+4m$ であり、殻上では $-m+4m=3m$ となる。
:::
::::

::::quiz
軟放射と仮想補正を足すと消えるものを選べ。

- [x] 赤外調整質量への対数依存
- [ ] 検出器の分解能への依存
- [ ] 硬い散乱の振幅そのもの

:::explanation
二つの積分は下端の調整質量と角係数を共有し、確率で符号が逆である。積分の上端が異なるため、分解能を含む有限な対数が残る。
:::
::::

::::quiz{answer="8" tol="1e-6"}
先頭対数で $1S$ の Lamb シフトは $2S$ の何倍か。

:::explanation
$|\psi_{nS}(0)|^2$ が $1/n^3$ に比例するので比は $2^3=8$ となる。
:::
::::

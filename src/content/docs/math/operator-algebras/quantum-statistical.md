---
title: 量子統計力学への応用
description: Weyl 関係、Fock 表現、準局所環と KMS 状態を結び、局所相互作用の伝播を評価する。
notation:
  - symbol: 'W(q,p)'
    meaning: 位置と運動量の移動を表す Weyl ユニタリ
  - symbol: '\mathcal F_-(H)'
    meaning: フェルミ粒子の反対称 Fock 空間
  - symbol: '\mathcal A_\Lambda'
    meaning: 有限領域 Λ の局所作用素環
  - symbol: '\Phi(Z)'
    meaning: 有限部分集合 Z に支えを持つ相互作用のエネルギー
---

## この章の目標

- 正準交換関係を有界な Weyl ユニタリで表す。
- 有限自由度の正則 Weyl 表現の一意性を示す。
- CAR の反対称 Fock 表現とスピン準局所環を構成する。
- Gibbs 状態の KMS 条件と局所的な情報伝播の評価を導く。

## CCR 環と Weyl 関係

位置 $Q$ と運動量 $P$ は一般に非有界で、積の定義域を指定する必要がある。これらの指数なら有界なユニタリであり、積を全空間で扱える。以下、内積は第 1 引数について線形であり、$\hbar$ を明示する。

:::definition[Weyl 関係]{#def-oa-weyl-ccr}
実ベクトル空間 $V=\R^{2d}$ の $v=(q,p)$ と $w=(q',p')$ に対し $\sigma(v,w)=p\cdot q'-q\cdot p'$ とする。:term[Weyl 関係]{en="Weyl relations" yomi="わいるかんけい"} は
$$
W(v)W(w)=e^{-i\sigma(v,w)/(2\hbar)}W(v+w),\qquad W(v)^*=W(-v),\qquad W(0)=1
$$
である。この関係を満たす普遍的ユニタリから作る単位的 C*-環をここでは Weyl CCR 環と呼ぶ。普遍性には強連続性を含めない。$t\mapsto W(tv)$ が各 $v$ で強連続な表現を :term[正則表現]{en="regular representation" yomi="せいそくひょうげん"} という。
:::

普遍環は有限和の *-代数の全表現での上限半ノルムを取り、その核で割って完備化すれば得られる。ノルムの上限は係数の絶対値の和以下である。$q$ は長さ、$p$ は運動量の単位を持ち、$p\cdot q$ と $\hbar$ は同じ作用の単位を持つ。

:::example[Schrödinger 表現]
$L^2(\R^d)$ で
$$
(W(q,p)f)(x)=e^{ip\cdot(x+q/2)/\hbar}f(x+q)
$$
とする。二度作用させると全位相と $W(q+q',p+p')$ の位相の差は $(q\cdot p'-p\cdot q')/(2\hbar)$ である。平行移動と位相乗算はユニタリであり、連続コンパクト台関数上の収束と収縮性から強連続である。一方 $p\ne0$ に対し $\norm{W(0,p)-1}=2$ なのでノルム連続でない。試験関数上で微分すると $Q_jf=x_jf$、$P_jf=-i\hbar\partial_jf$ となり、$[Q_j,P_k]f=i\hbar\delta_{jk}f$ を得る。
:::

### 正則性とノルム連続性を分ける

非退化なシンプレクティック形式で $v\ne0$ とする。Weyl 関係を二度用いると
$$
W(w)W(v)W(w)^*=e^{-i\sigma(w,v)/\hbar}W(v)
$$
になる。$w$ を動かすと係数は円周全体を走る。ユニタリ共役はスペクトルを変えないので、非空な $\sigma(W(v))$ は回転不変で、全円周である。従って $\norm{W(v)-1}=2$ である。正則表現における強連続性と C*-ノルムでの連続性が違う理由は、最大の変化を示すベクトルが $v$ ごとに変わるためである。

:::example[正準対の交換で符号を読む]
$\sigma((q,p),(q',p'))=pq'-qp'$ の規約で $v=(q,0)$、$w=(0,p)$ とする。$W(v)W(w)$ の位相は $e^{iqp/(2\hbar)}$、逆順の位相は $e^{-iqp/(2\hbar)}$ なので、両者の比は $e^{iqp/\hbar}$ である。従って $qp=2\pi\hbar$ なら二つのユニタリは可換であり、$qp=\pi\hbar$ なら反可換になる。$qp$ は作用の次元を持ち、位相に現れるのは無次元の比である。
:::

### 普遍 Weyl 環の単純性

:::proposition[非退化な Weyl CCR 環の単純性]{#prop-oa-weyl-simple}
ここで定義した非退化形式上の単位的 Weyl CCR 環は単純である。従ってその非零表現は忠実である。
:::

:::proof
$V$ をまず離散加法群と見て、$\ell^2(V)$ 上に $W(v)\delta_w=e^{-i\sigma(v,w)/(2\hbar)}\delta_{v+w}$ と置く。各 $W(v)$ は位相付きの基底置換である。形式の双線形性により Weyl 関係が成り立ち、$\delta_0$ への作用で有限 Weyl 和の線形独立性が分かる。

離散群 $V$ の指標の集合 $\widehat{V}$ は $\mathbb T^V$ の閉部分群である。積を保存する条件は座標ごとの閉条件なのでコンパクトである。各指標 $\chi$ に対する $W(v)\mapsto\chi(v)W(v)$ は Weyl 関係を保ち、逆も同じ形である。従って普遍性から自己同型 $\gamma_\chi$ が定まる。有限和上で連続であり等長性から環全体でもノルム連続になる。Haar 平均 $E(x)=\int\gamma_\chi(x)\,\dd\chi$ は定数係数を取る。非零の $v$ には $\chi(v)\ne1$ となる指標がある。例えば $V$ 上の実線形汎関数の指数を取れる。その指標で積分を平行移動すれば $W(v)$ の平均が零になる。有限和の稠密性により $E$ の像はスカラーである。$x\ge0$ で $E(x)=0$ なら、任意の状態 $\varphi$ に対し非負連続関数 $\chi\mapsto\varphi(\gamma_\chi(x))$ の積分が零になる。Haar 測度は非空開集合に正の値を取るため、この関数は恒等指標でも零である。状態が正元を分離するため $x=0$ となる。

他方、実ベクトル $w$ による $\operatorname{Ad}W(w)$ の平均を、$[-R,R]^{2d}$ の規格化 Lebesgue 積分として取る。$W(v)$ に対する平均は、$v\ne0$ なら非退化性で少なくとも一方向に非零周波数を持ち、その方向の指数関数の平均が零へ行く。$v=0$ では係数は $1$ のままである。従ってこの平均は有限和上で $E$ に収束し、縮小性により全環でもノルム収束する。閉両側イデアル内の正元にこの内的共役平均を取れば極限もそのイデアルにある。非零正元なら忠実性により $E(x)$ は正のスカラー倍の単位元なので、非零イデアルは全環となる。
:::

:::example[位置移動と運動量移動を入れ替える]
長さの基準を $L>0$、運動量の基準を $P=\pi\hbar/(2L)$ とする。位置・運動量の組の数字はこれらを単位とする。二つの Weyl ユニタリを掛ける順序を比較する。

係数の組を $(a,b),(c,d)$ と置けば、実際の組は $(aL,bP),(cL,dP)$ である。章の規約で $\sigma(v,w)=(bc-ad)LP$ となる。

$$
\begin{aligned}
(a,b,c,d)&=\left( 1, \  0, \  0, \  1\right)\\
bc-ad&=-1\\
\sigma(v,w)/\hbar&=- \frac{\pi}{2}
\end{aligned}
$$

$W(v)W(w)$ では $-\sigma(v,w)/(2\hbar)$ を指数の位相角に入れる。逆順では $\sigma$ の符号が反転するので位相は複素共役になる。

$$
\begin{aligned}
e^{-i\sigma(v,w)/(2\hbar)}&=e^{\frac{i \pi}{4}}\\
e^{-i\sigma(w,v)/(2\hbar)}&=e^{- \frac{i \pi}{4}}\\
e^{-i\sigma(v,w)/\hbar}&=i
\end{aligned}
$$

両積のユニタリ部分は同じ $W(v+w)$ である。そこで積の比を取ると、残るのは最後の位相だけである。

$$
\begin{aligned}
(a+c,b+d)&=\left( 1, \  1\right)
\end{aligned}
$$

$bc-ad=0$ の場合は位相なしで可換である。零でなくても角が $2\pi$ の整数倍なら可換になる。形式の値が零かどうかと位相が一かどうかを区別する。
:::

:::example[一般の二つの平面方向]
長さの基準を $L>0$、運動量の基準を $P=\pi\hbar/(2L)$ とする。位置・運動量の組の数字はこれらを単位とする。二つの Weyl ユニタリを掛ける順序を比較する。

係数の組を $(a,b),(c,d)$ と置けば、実際の組は $(aL,bP),(cL,dP)$ である。章の規約で $\sigma(v,w)=(bc-ad)LP$ となる。

$$
\begin{aligned}
(a,b,c,d)&=\left( 2, \  1, \  1, \  2\right)\\
bc-ad&=-3\\
\sigma(v,w)/\hbar&=- \frac{3 \pi}{2}
\end{aligned}
$$

$W(v)W(w)$ では $-\sigma(v,w)/(2\hbar)$ を指数の位相角に入れる。逆順では $\sigma$ の符号が反転するので位相は複素共役になる。

$$
\begin{aligned}
e^{-i\sigma(v,w)/(2\hbar)}&=e^{\frac{3 i \pi}{4}}\\
e^{-i\sigma(w,v)/(2\hbar)}&=e^{- \frac{3 i \pi}{4}}\\
e^{-i\sigma(v,w)/\hbar}&=- i
\end{aligned}
$$

両積のユニタリ部分は同じ $W(v+w)$ である。そこで積の比を取ると、残るのは最後の位相だけである。

$$
\begin{aligned}
(a+c,b+d)&=\left( 3, \  3\right)
\end{aligned}
$$

$bc-ad=0$ の場合は位相なしで可換である。零でなくても角が $2\pi$ の整数倍なら可換になる。形式の値が零かどうかと位相が一かどうかを区別する。
:::

:::example[同一直線に並ぶ方向]
長さの基準を $L>0$、運動量の基準を $P=\pi\hbar/(2L)$ とする。位置・運動量の組の数字はこれらを単位とする。二つの Weyl ユニタリを掛ける順序を比較する。

係数の組を $(a,b),(c,d)$ と置けば、実際の組は $(aL,bP),(cL,dP)$ である。章の規約で $\sigma(v,w)=(bc-ad)LP$ となる。

$$
\begin{aligned}
(a,b,c,d)&=\left( 1, \  2, \  2, \  4\right)\\
bc-ad&=0\\
\sigma(v,w)/\hbar&=0
\end{aligned}
$$

$W(v)W(w)$ では $-\sigma(v,w)/(2\hbar)$ を指数の位相角に入れる。逆順では $\sigma$ の符号が反転するので位相は複素共役になる。

$$
\begin{aligned}
e^{-i\sigma(v,w)/(2\hbar)}&=1\\
e^{-i\sigma(w,v)/(2\hbar)}&=1\\
e^{-i\sigma(v,w)/\hbar}&=1
\end{aligned}
$$

両積のユニタリ部分は同じ $W(v+w)$ である。そこで積の比を取ると、残るのは最後の位相だけである。

$$
\begin{aligned}
(a+c,b+d)&=\left( 3, \  6\right)
\end{aligned}
$$

$bc-ad=0$ の場合は位相なしで可換である。零でなくても角が $2\pi$ の整数倍なら可換になる。形式の値が零かどうかと位相が一かどうかを区別する。
:::

:::example[大きさが異なる逆向きの成分]
長さの基準を $L>0$、運動量の基準を $P=\pi\hbar/(2L)$ とする。位置・運動量の組の数字はこれらを単位とする。二つの Weyl ユニタリを掛ける順序を比較する。

係数の組を $(a,b),(c,d)$ と置けば、実際の組は $(aL,bP),(cL,dP)$ である。章の規約で $\sigma(v,w)=(bc-ad)LP$ となる。

$$
\begin{aligned}
(a,b,c,d)&=\left( 1, \  1, \  -2, \  2\right)\\
bc-ad&=-4\\
\sigma(v,w)/\hbar&=- 2 \pi
\end{aligned}
$$

$W(v)W(w)$ では $-\sigma(v,w)/(2\hbar)$ を指数の位相角に入れる。逆順では $\sigma$ の符号が反転するので位相は複素共役になる。

$$
\begin{aligned}
e^{-i\sigma(v,w)/(2\hbar)}&=-1\\
e^{-i\sigma(w,v)/(2\hbar)}&=-1\\
e^{-i\sigma(v,w)/\hbar}&=1
\end{aligned}
$$

両積のユニタリ部分は同じ $W(v+w)$ である。そこで積の比を取ると、残るのは最後の位相だけである。

$$
\begin{aligned}
(a+c,b+d)&=\left( -1, \  3\right)
\end{aligned}
$$

$bc-ad=0$ の場合は位相なしで可換である。零でなくても角が $2\pi$ の整数倍なら可換になる。形式の値が零かどうかと位相が一かどうかを区別する。
:::

:::example[向きの取り方で共役になる位相]
長さの基準を $L>0$、運動量の基準を $P=\pi\hbar/(2L)$ とする。位置・運動量の組の数字はこれらを単位とする。二つの Weyl ユニタリを掛ける順序を比較する。

係数の組を $(a,b),(c,d)$ と置けば、実際の組は $(aL,bP),(cL,dP)$ である。章の規約で $\sigma(v,w)=(bc-ad)LP$ となる。

$$
\begin{aligned}
(a,b,c,d)&=\left( 0, \  1, \  1, \  1\right)\\
bc-ad&=1\\
\sigma(v,w)/\hbar&=\frac{\pi}{2}
\end{aligned}
$$

$W(v)W(w)$ では $-\sigma(v,w)/(2\hbar)$ を指数の位相角に入れる。逆順では $\sigma$ の符号が反転するので位相は複素共役になる。

$$
\begin{aligned}
e^{-i\sigma(v,w)/(2\hbar)}&=e^{- \frac{i \pi}{4}}\\
e^{-i\sigma(w,v)/(2\hbar)}&=e^{\frac{i \pi}{4}}\\
e^{-i\sigma(v,w)/\hbar}&=- i
\end{aligned}
$$

両積のユニタリ部分は同じ $W(v+w)$ である。そこで積の比を取ると、残るのは最後の位相だけである。

$$
\begin{aligned}
(a+c,b+d)&=\left( 1, \  2\right)
\end{aligned}
$$

$bc-ad=0$ の場合は位相なしで可換である。零でなくても角が $2\pi$ の整数倍なら可換になる。形式の値が零かどうかと位相が一かどうかを区別する。
:::

## Stone–von Neumann の一意性定理

:::theorem[有限自由度の正則表現]{#thm-oa-stone-von-neumann}
上記の非退化な $\sigma$ と固定した $\hbar>0$ に対し、Weyl 関係の非零な既約正則表現は Schrödinger 表現とユニタリ同値である。一般の正則表現はその多重度付き表現となる。
:::
:::proof
正則性と Weyl 関係により添字全体に関して強連続である。実際、各座標方向の作用の有限積と連続な位相で $W(v)$ を書ける。$F\in L^1(\R^{2d})$ に対する $T_F=\int F(v)W(v)\,\dd v$ はベクトル上の積分で定まり、$\norm{T_F}\le\norm F_1$ である。積は Weyl 関係を使い、$T_F T_G=T_{F\star G}$、$(F\star G)(v)=\int F(w)G(v-w)e^{-i\sigma(w,v-w)/(2\hbar)}\,\dd w$ となる。対合は $F^*(v)=\overline{F(-v)}$ である。この二つの式は Fubini の定理と変数変換から従い、Schwartz 関数上で閉じる。積の各導関数は、急減少する二関数と座標の多項式で支配されるためである。

Schrödinger 表現では Schwartz 関数 $F$ の積分核は
$$
K_F(x,y)=\int_{\R^d}F(y-x,p)\,e^{ip\cdot(x+y)/(2\hbar)}\,\dd p
$$
である。部分 Fourier 変換と可逆な線形座標変換だから、これは Schwartz 関数空間同士の連続な線形全単射である。特に核による写像は単射であり、階数一核 $\xi(x)\overline{\eta(y)}$（$\xi,\eta$ は Schwartz 関数）の逆像は、それらの積と対合を正確に保つ。これらの核の線形包は Schwartz 位相で稠密である。実際、独立な $x,y$ の切断関数でコンパクト台へ近似し、十分大きな箱で周期化して Fourier 級数を取る。滑らかな周期関数の係数は部分積分で任意の冪より速く減衰するため、部分和は全ての固定階の導関数に関し一様収束する。元の台で $1$ になる積形の切断関数を最後に掛けると、分離した変数の滑らかなコンパクト台関数の有限和が得られる。従って階数一核の逆像の線形包は Schwartz 空間、ひいては $L^1$ で稠密である。

任意の正則表現の積分写像をこの階数一核の代数に制限する。有限個の核の元はそれらのベクトルが張る有限次元空間の行列環に入り、各 *-表現はそこで作用素ノルムに関し縮小的である。そのため表現は $K(L^2(\R^d))$ へ延びる。元の $T_F$ は $L^1$ 有界性によりこの延長と一致する。原点に集中する積分 $1$ の非負な $F$ で $T_F\to1$ が強収束するから非退化である。:ref[thm-oa-compact-representations] により積分表現は多重度付き通常表現となる。$W(v)T_F$ は $F$ の平行移動と Weyl 位相で定まる別の積分である。従って積分の像を絡めるユニタリは稠密な $\spn T_FH$ 上で $W(v)$ も絡める。既約なら多重度は $1$ である。
:::

有限次元の仮定は $V=\R^{2d}$ 上の Haar 測度と Fourier 変換に使った。無限自由度では、次の Fock 表現から正則な非同値表現を構成できる。

:::example[正則性を外した表示]
$d=1$ で、実数を添字とする $\ell^2(\R)$ の正規直交基底 $(\delta_y)_{y\in\R}$ を考える。$W(q,p)\delta_y=e^{ip(y-q/2)/\hbar}\delta_{y-q}$ は基底を位相付きで並べ替えるユニタリである。二回適用すると、同じ添字への移動に対して合成の位相差は $(qp'-pq')/(2\hbar)$ となるので Weyl 関係を満たす。一方 $q\ne0$ なら $\norm{W(q,0)\delta_0-\delta_0}=\sqrt2$ である。従って $q\to0$ で強連続にならず、上の定理の対象には入らない。

この表現は既約である。$W(0,p)$（全 $p$）との可換性から作用素の $(x,y)$ 成分は $x\ne y$ なら零である。実際、対応する対角乗数 $e^{ipx/\hbar}$ と $e^{ipy/\hbar}$ が異なる $p$ を選べる。更に $W(q,0)$ との可換性は全座標で対角値を同じにする。可換子がスカラーなので :ref[prop-oa-schur] が適用できる。従って正則性を外すと、有限自由度であっても非同値な既約表現がある。
:::

### 正則表現の積分像と Weyl ユニタリ

有限自由度の Schrödinger 表現で $W(q,p)$ を平滑でコンパクト台の係数で積分すると、各作用素は平行移動と乗算の積を重ねたものになる。この章の規約では、積分核の変数は $y=x+q$ であり、残った $p$ の積分は Fourier 変換になる。従って核は差 $x-y$ の方向でコンパクト台を持ち、中心 $(x+y)/2$ の方向で急減少する。平方積分可能な核になり、Hilbert–Schmidt 作用素である。実際、核を矩形の指示関数の有限和で $L^2$ 近似すれば、対応する有限階数作用素との作用素ノルム差は核の $L^2$ 差以下になる。各 Weyl ユニタリは無限次元空間の等長全射であり、積分像のコンパクト作用素環に対して乗数として作用する。

### 対称 Fock 空間での正則表現

複素一粒子空間 $K$ の :term[対称 Fock 空間]{en="symmetric Fock space" yomi="たいしょうふぉっくくうかん"} を $\mathcal F_+(K)=\bigoplus_{n\ge0}K^{\otimes_s n}$ とする。指数ベクトル
$$
e(f)=\bigoplus_{n\ge0}\frac{f^{\otimes n}}{\sqrt{n!}}
$$
は $\inner{e(f)}{e(g)}=\exp(\inner f g)$ を満たす。これらの線形包は稠密である。実際、すべての $e(tf)$ と直交すれば、$t$ のべき級数の係数を比較して各粒子数成分がすべての $f^{\otimes n}$ と直交する。対称テンソルはこれらの冪で生成される。例えば $(t_1f_1+\cdots+t_nf_n)^{\otimes n}$ の係数 $t_1\cdots t_n$ を取れば任意の対称化した単純テンソルが得られる。

指数ベクトル上で
$$
D(z)e(f)=\exp\bigl(-\norm z^2/2-\inner f z\bigr)e(f+z)
$$
と置く。指数ベクトルの内積の式を代入すると、この対応は内積を保つ。従って線形包に矛盾なく延び、その閉包で等長になる。二回作用させて係数を比較すると
$$
D(z)D(w)=\exp(i\Im\inner z w)D(z+w)
$$
を得る。$D(-z)$ が逆なのでユニタリである。$\norm{e(f)-e(g)}^2$ は内積の式から $f,g$ に連続であり、稠密な線形包と等長性から $D(z)$ は $z$ について強連続である。有限自由度で長さ尺度 $\ell>0$ を選び、$z=(q/\ell-i\ell p/\hbar)/\sqrt 2$ と置けば $\Im\inner z w=-\sigma(v,w)/(2\hbar)$ となる。従って $D(z)$ はこの章の Weyl 関係を満たす。

:::example[無限個のモードで一意性が破れる]
$K=\ell^2(\N)$ とし、添字 $z$ を有限台の複素列に制限する。Fock 表現 $D(z)$ は既約である。これを示すため、指数ベクトルの内積から得る $\mathcal F_+(K_1\oplus K_2)\cong\mathcal F_+(K_1)\otimes\mathcal F_+(K_2)$ を使う。有限個のモードからなる各因子では、真空を $D(z)$ で移すと全指数ベクトルが得られる。また最初の $N$ モードの正則表現では真空への射影がその von Neumann 環に入る。実際、$D(z)$ の Gaussian 平均
$$
\pi^{-N}\int_{\C^N}e^{-\norm z^2/2}D(z)\,\dd z
$$
は、指数ベクトル間で評価すると真空への射影に等しい。積分する関数の行列係数は $\exp(\inner f g-\norm z^2-\inner f z+\inner z g)$ であり、正規化した複素 Gaussian 積分は $\exp(-\inner f g)$ を与える。これは $e^{-\inner f z}$ と $e^{\inner z g}$ を展開し、各座標で角度が打ち消し合う項だけを積分すれば従う。従って可換子の元は最初の $N$ モードの真空への射影と可換である。これらの射影は全真空への射影に強収束する。可換子の元は全真空をそのスカラー倍に移し、さらに $D(z)\Omega$ の稠密な線形包上で同じスカラーとなる。

有限台の $z$ に $\chi(z)=\exp(i\sum_j\Re z_j)$ を定める。$D'(z)=\chi(z)D(z)$ も同じ Weyl 関係を満たし、既約で、各直線の上で強連続である。ところが最初の $N$ 成分が $\pi/N$ の列 $z_N$ はノルムで零へ行く一方、$\chi(z_N)=-1$ である。従って $D(z_N)\to1$、$D'(z_N)\to-1$ が強く成り立つ。ユニタリ同値であればこの二つの極限は同じになるため、両表現は非同値である。
:::

:::example[標準の Gaussian 波束への射影]
Schrödinger 表現で Gaussian 関数を係数にして Weyl ユニタリを積分し、階数一射影を作る。長さの単位を $L$ とし $r=q/L,\ w=pL/\hbar$、中心座標を $x=(X+Y)/(2L)$ とする。

幅を $kL$ として積分係数 $F(q,p)=(2\pi\hbar)^{-1}\exp[-r^2/(4k^2)-k^2w^2/4]$ を選ぶ。$p$ の積分で変数を $w$ に替えると $\dd p=(\hbar/L)\dd w$ になる。

$$
\begin{aligned}
k&=1\\
2\pi\hbar\,F(q,p)&=e^{- \frac{r^{2}}{4} - \frac{w^{2}}{4}}
\end{aligned}
$$

核 $K_F(X,Y)$ は前節の式により $r=(Y-X)/L$ を代入した後、運動量方向の Fourier 変換で得られる。Gaussian の積分を平方完成すれば中心変数に同じ型の関数が出る。

$$
\begin{aligned}
\frac1{2\pi}\int_{\R} e^{-k^2 w^2/4+iwx}\,\dd w&=\frac{1}{k\sqrt{\pi}}e^{-x^2/k^2}\\
L K_F(X,Y)&=\frac{e^{- \frac{r^{2}}{4} - x^{2}}}{\sqrt{\pi}}
\end{aligned}
$$

$x^2+r^2/4=(X^2+Y^2)/(2L^2)$ なので核は $g(X)\overline{g(Y)}$ と因数分解する。ここで $g(X)=(\pi k^2L^2)^{-1/4}e^{-X^2/(2k^2L^2)}$ であり、二乗の積分は $1$ である。従って積分作用素は $g$ への直交射影になる。

この計算は積分像に有限階数作用素が現れる理由を具体化する。位置幅を広げると係数 $F$ の運動量方向の幅は狭くなり、Fourier 変換後の位置核の幅は広がる。
:::

:::example[幅を二倍に広げた波束]
Schrödinger 表現で Gaussian 関数を係数にして Weyl ユニタリを積分し、階数一射影を作る。長さの単位を $L$ とし $r=q/L,\ w=pL/\hbar$、中心座標を $x=(X+Y)/(2L)$ とする。

幅を $kL$ として積分係数 $F(q,p)=(2\pi\hbar)^{-1}\exp[-r^2/(4k^2)-k^2w^2/4]$ を選ぶ。$p$ の積分で変数を $w$ に替えると $\dd p=(\hbar/L)\dd w$ になる。

$$
\begin{aligned}
k&=2\\
2\pi\hbar\,F(q,p)&=e^{- \frac{r^{2}}{16} - w^{2}}
\end{aligned}
$$

核 $K_F(X,Y)$ は前節の式により $r=(Y-X)/L$ を代入した後、運動量方向の Fourier 変換で得られる。Gaussian の積分を平方完成すれば中心変数に同じ型の関数が出る。

$$
\begin{aligned}
\frac1{2\pi}\int_{\R} e^{-k^2 w^2/4+iwx}\,\dd w&=\frac{1}{k\sqrt{\pi}}e^{-x^2/k^2}\\
L K_F(X,Y)&=\frac{e^{- \frac{r^{2}}{16} - \frac{x^{2}}{4}}}{2 \sqrt{\pi}}
\end{aligned}
$$

$x^2+r^2/4=(X^2+Y^2)/(2L^2)$ なので核は $g(X)\overline{g(Y)}$ と因数分解する。ここで $g(X)=(\pi k^2L^2)^{-1/4}e^{-X^2/(2k^2L^2)}$ であり、二乗の積分は $1$ である。従って積分作用素は $g$ への直交射影になる。

この計算は積分像に有限階数作用素が現れる理由を具体化する。位置幅を広げると係数 $F$ の運動量方向の幅は狭くなり、Fourier 変換後の位置核の幅は広がる。
:::

:::example[広い位置波束と狭い運動量係数]
Schrödinger 表現で Gaussian 関数を係数にして Weyl ユニタリを積分し、階数一射影を作る。長さの単位を $L$ とし $r=q/L,\ w=pL/\hbar$、中心座標を $x=(X+Y)/(2L)$ とする。

幅を $kL$ として積分係数 $F(q,p)=(2\pi\hbar)^{-1}\exp[-r^2/(4k^2)-k^2w^2/4]$ を選ぶ。$p$ の積分で変数を $w$ に替えると $\dd p=(\hbar/L)\dd w$ になる。

$$
\begin{aligned}
k&=3\\
2\pi\hbar\,F(q,p)&=e^{- \frac{r^{2}}{36} - \frac{9 w^{2}}{4}}
\end{aligned}
$$

核 $K_F(X,Y)$ は前節の式により $r=(Y-X)/L$ を代入した後、運動量方向の Fourier 変換で得られる。Gaussian の積分を平方完成すれば中心変数に同じ型の関数が出る。

$$
\begin{aligned}
\frac1{2\pi}\int_{\R} e^{-k^2 w^2/4+iwx}\,\dd w&=\frac{1}{k\sqrt{\pi}}e^{-x^2/k^2}\\
L K_F(X,Y)&=\frac{e^{- \frac{r^{2}}{36} - \frac{x^{2}}{9}}}{3 \sqrt{\pi}}
\end{aligned}
$$

$x^2+r^2/4=(X^2+Y^2)/(2L^2)$ なので核は $g(X)\overline{g(Y)}$ と因数分解する。ここで $g(X)=(\pi k^2L^2)^{-1/4}e^{-X^2/(2k^2L^2)}$ であり、二乗の積分は $1$ である。従って積分作用素は $g$ への直交射影になる。

この計算は積分像に有限階数作用素が現れる理由を具体化する。位置幅を広げると係数 $F$ の運動量方向の幅は狭くなり、Fourier 変換後の位置核の幅は広がる。
:::

:::example[位置を強く局在させる波束]
Schrödinger 表現で Gaussian 関数を係数にして Weyl ユニタリを積分し、階数一射影を作る。長さの単位を $L$ とし $r=q/L,\ w=pL/\hbar$、中心座標を $x=(X+Y)/(2L)$ とする。

幅を $kL$ として積分係数 $F(q,p)=(2\pi\hbar)^{-1}\exp[-r^2/(4k^2)-k^2w^2/4]$ を選ぶ。$p$ の積分で変数を $w$ に替えると $\dd p=(\hbar/L)\dd w$ になる。

$$
\begin{aligned}
k&=\frac{1}{2}\\
2\pi\hbar\,F(q,p)&=e^{- r^{2} - \frac{w^{2}}{16}}
\end{aligned}
$$

核 $K_F(X,Y)$ は前節の式により $r=(Y-X)/L$ を代入した後、運動量方向の Fourier 変換で得られる。Gaussian の積分を平方完成すれば中心変数に同じ型の関数が出る。

$$
\begin{aligned}
\frac1{2\pi}\int_{\R} e^{-k^2 w^2/4+iwx}\,\dd w&=\frac{1}{k\sqrt{\pi}}e^{-x^2/k^2}\\
L K_F(X,Y)&=\frac{2 e^{- r^{2} - 4 x^{2}}}{\sqrt{\pi}}
\end{aligned}
$$

$x^2+r^2/4=(X^2+Y^2)/(2L^2)$ なので核は $g(X)\overline{g(Y)}$ と因数分解する。ここで $g(X)=(\pi k^2L^2)^{-1/4}e^{-X^2/(2k^2L^2)}$ であり、二乗の積分は $1$ である。従って積分作用素は $g$ への直交射影になる。

この計算は積分像に有限階数作用素が現れる理由を具体化する。位置幅を広げると係数 $F$ の運動量方向の幅は狭くなり、Fourier 変換後の位置核の幅は広がる。
:::

:::example[異なる尺度での正規化の確認]
Schrödinger 表現で Gaussian 関数を係数にして Weyl ユニタリを積分し、階数一射影を作る。長さの単位を $L$ とし $r=q/L,\ w=pL/\hbar$、中心座標を $x=(X+Y)/(2L)$ とする。

幅を $kL$ として積分係数 $F(q,p)=(2\pi\hbar)^{-1}\exp[-r^2/(4k^2)-k^2w^2/4]$ を選ぶ。$p$ の積分で変数を $w$ に替えると $\dd p=(\hbar/L)\dd w$ になる。

$$
\begin{aligned}
k&=\frac{3}{2}\\
2\pi\hbar\,F(q,p)&=e^{- \frac{r^{2}}{9} - \frac{9 w^{2}}{16}}
\end{aligned}
$$

核 $K_F(X,Y)$ は前節の式により $r=(Y-X)/L$ を代入した後、運動量方向の Fourier 変換で得られる。Gaussian の積分を平方完成すれば中心変数に同じ型の関数が出る。

$$
\begin{aligned}
\frac1{2\pi}\int_{\R} e^{-k^2 w^2/4+iwx}\,\dd w&=\frac{1}{k\sqrt{\pi}}e^{-x^2/k^2}\\
L K_F(X,Y)&=\frac{2 e^{- \frac{r^{2}}{9} - \frac{4 x^{2}}{9}}}{3 \sqrt{\pi}}
\end{aligned}
$$

$x^2+r^2/4=(X^2+Y^2)/(2L^2)$ なので核は $g(X)\overline{g(Y)}$ と因数分解する。ここで $g(X)=(\pi k^2L^2)^{-1/4}e^{-X^2/(2k^2L^2)}$ であり、二乗の積分は $1$ である。従って積分作用素は $g$ への直交射影になる。

この計算は積分像に有限階数作用素が現れる理由を具体化する。位置幅を広げると係数 $F$ の運動量方向の幅は狭くなり、Fourier 変換後の位置核の幅は広がる。
:::

## CAR 環と Fock 表現

:::definition[反対称 Fock 空間]{#def-oa-fermi-fock}
一粒子 Hilbert 空間 $K$ に対し、$\mathcal F_-(K)=\bigoplus_{n\ge0}\bigwedge^nK$ を :term[反対称 Fock 空間]{en="antisymmetric Fock space" yomi="はんたいしょうふぉっくくうかん"} という。真空 $\Omega$ は $n=0$ の単位ベクトルである。外積の内積を $\inner{f_1\wedge\cdots\wedge f_n}{g_1\wedge\cdots\wedge g_n}=\det(\inner{f_i}{g_j})$ とし、生成作用素を $c^*(f)\eta=f\wedge\eta$、消滅作用素を $c(f)=(c^*(f))^*$ とする。
:::

$c^*(f)$ は $f$ に線形、$c(f)$ は共役線形である。$f/\norm f$ を正規直交基底の先頭に置くと、$c^*(f)$ はその基底元を含まない外積基底に対しその元を先頭に加える等長写像の $\norm f$ 倍であり、既に含む外積基底では零となる。従って有界で $\norm{c^*(f)}=\norm f$ である。その随伴は
$$
c(f)(g_1\wedge\cdots\wedge g_n)=
\sum_{j=1}^n(-1)^{j-1}\inner{g_j}{f}
\,g_1\wedge\cdots\wedge\widehat{g_j}\wedge\cdots\wedge g_n
$$
である。$\widehat{g_j}$ はその因子を除くことを示す。これは行列式を第 1 列で展開して随伴を確認できる。従って $c(f)c^*(g)+c^*(g)c(f)=\inner g f\,1$ であり、生成同士は外積の反対称性、消滅同士はその随伴で反交換する。有限次元 $K$ では占有の $0,1$ の列が基底を与え、前の [AF 環](/math/operator-algebras/af-algebras/) の有限 CAR 行列模型と同じ全行列環となる。無限次元の場合には有限次元部分空間の合併を用いる。$\norm{c^*(f)-c^*(g)}=\norm{f-g}$ により稠密な一粒子ベクトルから全ての生成作用素を近似できる。

:::example[二粒子の符号]
直交する単位ベクトル $e_1,e_2$ について $c^*(e_2)c^*(e_1)\Omega=-e_1\wedge e_2$ である。$c(e_2)(e_1\wedge e_2)=-e_1$ となる。同じ状態 $e_1$ を二度生成すると零になる。
:::

### Jordan–Wigner の行列表示

$Z=\diag(1,-1)$、$c=E_{12}$ とし、
$$
a_j=Z^{\otimes(j-1)}\otimes c\otimes I^{\otimes(n-j)}
$$
と置く。$Zc=-cZ$、$Zc^*=-c^*Z$ なので、$j<k$ の二つの積を交換したときは第 $j$ 因子だけから符号が出る。同じ添字では $c^2=0$ と $cc^*+c^*c=I$ が成り立つ。従ってこの行列族は CAR を満たす。さらに $Z$ は $I-2c^*c$ である。第 1 因子の $c$ と $c^*$ を生成し、それらから $Z$ を作り、次の生成元の前の $Z$ を順に消すと各因子の $c,c^*$ が得られる。これらは全テンソル行列環を生成する。

:::example[二サイトの二粒子を順に消す]
占有数基底 $|n_1,\ldots,n_l\rangle$ 上の Jordan–Wigner 表示で $a_j a_k+a_k a_j$ を試す。消滅は対象の座標を $1$ から $0$ に移すが、それより前の占有数の偶奇も使う。

対象座標が既に $0$ なら結果は零である。それ以外では添字 $j$ より前の $Z$ たちから $(-1)^{\sum_{i<j}n_i}$ が付く。最初の一回の作用を書き出す。

| 入力の占有数 | 最初の作用 | 結果 |
| --- | --- | --- |
| $\left( 1, \  1\right)$ | $a_{1}$ | $1$ 倍の $\left( 0, \  1\right)$ |
| $\left( 1, \  1\right)$ | $a_{2}$ | $-1$ 倍の $\left( 1, \  0\right)$ |

二回目の符号は元の占有数でなく、一回目の後の占有数から計算する。そのため低い添字を先に消す順序と、後で消す順序では符号が逆になる。

| 積 | 結果 |
| --- | --- |
| $a_ja_k$ | $-1$ 倍の $\left( 0, \  0\right)$ |
| $a_ka_j$ | $1$ 倍の $\left( 0, \  0\right)$ |

二サイトが両方占有なら終状態は同じで二つの符号が相殺する。どちらかが空なら両積とも零である。したがってこの基底に限らず全ての占有数ベクトル上で反交換子は零になる。

$$
\begin{aligned}
\{a_j,a_k\}&=0
\end{aligned}
$$

各サイトの単なるテンソル積の消滅作用素では、離れたサイトの元は可換になる。前方の $Z$ の列が反交換の符号を作る。
:::

:::example[離れた二サイトの粒子]
占有数基底 $|n_1,\ldots,n_l\rangle$ 上の Jordan–Wigner 表示で $a_j a_k+a_k a_j$ を試す。消滅は対象の座標を $1$ から $0$ に移すが、それより前の占有数の偶奇も使う。

対象座標が既に $0$ なら結果は零である。それ以外では添字 $j$ より前の $Z$ たちから $(-1)^{\sum_{i<j}n_i}$ が付く。最初の一回の作用を書き出す。

| 入力の占有数 | 最初の作用 | 結果 |
| --- | --- | --- |
| $\left( 1, \  0, \  1\right)$ | $a_{1}$ | $1$ 倍の $\left( 0, \  0, \  1\right)$ |
| $\left( 1, \  0, \  1\right)$ | $a_{3}$ | $-1$ 倍の $\left( 1, \  0, \  0\right)$ |

二回目の符号は元の占有数でなく、一回目の後の占有数から計算する。そのため低い添字を先に消す順序と、後で消す順序では符号が逆になる。

| 積 | 結果 |
| --- | --- |
| $a_ja_k$ | $-1$ 倍の $\left( 0, \  0, \  0\right)$ |
| $a_ka_j$ | $1$ 倍の $\left( 0, \  0, \  0\right)$ |

二サイトが両方占有なら終状態は同じで二つの符号が相殺する。どちらかが空なら両積とも零である。したがってこの基底に限らず全ての占有数ベクトル上で反交換子は零になる。

$$
\begin{aligned}
\{a_j,a_k\}&=0
\end{aligned}
$$

各サイトの単なるテンソル積の消滅作用素では、離れたサイトの元は可換になる。前方の $Z$ の列が反交換の符号を作る。
:::

:::example[前方の占有数が符号を変える]
占有数基底 $|n_1,\ldots,n_l\rangle$ 上の Jordan–Wigner 表示で $a_j a_k+a_k a_j$ を試す。消滅は対象の座標を $1$ から $0$ に移すが、それより前の占有数の偶奇も使う。

対象座標が既に $0$ なら結果は零である。それ以外では添字 $j$ より前の $Z$ たちから $(-1)^{\sum_{i<j}n_i}$ が付く。最初の一回の作用を書き出す。

| 入力の占有数 | 最初の作用 | 結果 |
| --- | --- | --- |
| $\left( 1, \  1, \  1\right)$ | $a_{2}$ | $-1$ 倍の $\left( 1, \  0, \  1\right)$ |
| $\left( 1, \  1, \  1\right)$ | $a_{3}$ | $1$ 倍の $\left( 1, \  1, \  0\right)$ |

二回目の符号は元の占有数でなく、一回目の後の占有数から計算する。そのため低い添字を先に消す順序と、後で消す順序では符号が逆になる。

| 積 | 結果 |
| --- | --- |
| $a_ja_k$ | $-1$ 倍の $\left( 1, \  0, \  0\right)$ |
| $a_ka_j$ | $1$ 倍の $\left( 1, \  0, \  0\right)$ |

二サイトが両方占有なら終状態は同じで二つの符号が相殺する。どちらかが空なら両積とも零である。したがってこの基底に限らず全ての占有数ベクトル上で反交換子は零になる。

$$
\begin{aligned}
\{a_j,a_k\}&=0
\end{aligned}
$$

各サイトの単なるテンソル積の消滅作用素では、離れたサイトの元は可換になる。前方の $Z$ の列が反交換の符号を作る。
:::

:::example[片方のサイトが空の場合]
占有数基底 $|n_1,\ldots,n_l\rangle$ 上の Jordan–Wigner 表示で $a_j a_k+a_k a_j$ を試す。消滅は対象の座標を $1$ から $0$ に移すが、それより前の占有数の偶奇も使う。

対象座標が既に $0$ なら結果は零である。それ以外では添字 $j$ より前の $Z$ たちから $(-1)^{\sum_{i<j}n_i}$ が付く。最初の一回の作用を書き出す。

| 入力の占有数 | 最初の作用 | 結果 |
| --- | --- | --- |
| $\left( 1, \  0, \  1\right)$ | $a_{1}$ | $1$ 倍の $\left( 0, \  0, \  1\right)$ |
| $\left( 1, \  0, \  1\right)$ | $a_{2}$ | $0$ |

二回目の符号は元の占有数でなく、一回目の後の占有数から計算する。そのため低い添字を先に消す順序と、後で消す順序では符号が逆になる。

| 積 | 結果 |
| --- | --- |
| $a_ja_k$ | $0$ |
| $a_ka_j$ | $0$ |

二サイトが両方占有なら終状態は同じで二つの符号が相殺する。どちらかが空なら両積とも零である。したがってこの基底に限らず全ての占有数ベクトル上で反交換子は零になる。

$$
\begin{aligned}
\{a_j,a_k\}&=0
\end{aligned}
$$

各サイトの単なるテンソル積の消滅作用素では、離れたサイトの元は可換になる。前方の $Z$ の列が反交換の符号を作る。
:::

:::example[三つの占有から二つを消す]
占有数基底 $|n_1,\ldots,n_l\rangle$ 上の Jordan–Wigner 表示で $a_j a_k+a_k a_j$ を試す。消滅は対象の座標を $1$ から $0$ に移すが、それより前の占有数の偶奇も使う。

対象座標が既に $0$ なら結果は零である。それ以外では添字 $j$ より前の $Z$ たちから $(-1)^{\sum_{i<j}n_i}$ が付く。最初の一回の作用を書き出す。

| 入力の占有数 | 最初の作用 | 結果 |
| --- | --- | --- |
| $\left( 0, \  1, \  1, \  1\right)$ | $a_{2}$ | $1$ 倍の $\left( 0, \  0, \  1, \  1\right)$ |
| $\left( 0, \  1, \  1, \  1\right)$ | $a_{4}$ | $1$ 倍の $\left( 0, \  1, \  1, \  0\right)$ |

二回目の符号は元の占有数でなく、一回目の後の占有数から計算する。そのため低い添字を先に消す順序と、後で消す順序では符号が逆になる。

| 積 | 結果 |
| --- | --- |
| $a_ja_k$ | $1$ 倍の $\left( 0, \  0, \  1, \  0\right)$ |
| $a_ka_j$ | $-1$ 倍の $\left( 0, \  0, \  1, \  0\right)$ |

二サイトが両方占有なら終状態は同じで二つの符号が相殺する。どちらかが空なら両積とも零である。したがってこの基底に限らず全ての占有数ベクトル上で反交換子は零になる。

$$
\begin{aligned}
\{a_j,a_k\}&=0
\end{aligned}
$$

各サイトの単なるテンソル積の消滅作用素では、離れたサイトの元は可換になる。前方の $Z$ の列が反交換の符号を作る。
:::

## 準局所環と量子スピン系

格子 $\Gamma=\Z^d$ の各点に $m$ 次元 Hilbert 空間を置く。有限 $\Lambda\subset\Gamma$ について $\mathcal A_\Lambda=\bigotimes_{x\in\Lambda}M_m$ とし、$\Lambda\subset\Lambda'$ なら補集合に恒等作用素を掛けて埋め込む。その帰納極限を :term[準局所環]{en="quasi-local algebra" yomi="じゅんきょくしょかん"} $\mathcal A$ と呼ぶ。有限支えの元は稠密である。互いに素な有限領域の元は異なるテンソル因子に作用するため可換である。サイト数 $k$ の局所環の行列サイズは $m^k$ である。

各有限集合 $Z$ に $\Phi(Z)=\Phi(Z)^*\in\mathcal A_Z$ を割り当て、有限領域の Hamiltonian を $H_\Lambda=\sum_{Z\subset\Lambda}\Phi(Z)$ とする。$\Phi(Z)$ はエネルギーの単位、$t$ は時間の単位を持ち、$\alpha_t^\Lambda(a)=e^{itH_\Lambda/\hbar}ae^{-itH_\Lambda/\hbar}$ では指数は無次元である。例えばスピン $1/2$ の一サイトで $H=\epsilon E_{22}$（$\epsilon>0$）とすると $\alpha_t(E_{12})=e^{-it\epsilon/\hbar}E_{12}$ である。

### 支持が交わる場合と交わらない場合

前の定義は、サイトごとに次元 $d_x$ を変えた局所環 $\bigotimes_{x\in\Lambda}M_{d_x}$ にも適用できる。包含では外側の因子に恒等を入れる。二つの作用素の支持が交わらなければ、共通の有限領域に入れたときそれぞれ異なるテンソル因子にしか作用しない。単純テンソル上で積の順序を逆にでき、有限和に広げれば交換子は零である。支持が一点でも重なれば非可換になり得る。例えば二サイトで $X=E_{12}+E_{21}$、$Z=\diag(1,-1)$ とし $a=X\otimes I$、$b=Z\otimes X$ とすると
$$
[a,b]=[X,Z]\otimes X=-2(E_{12}-E_{21})\otimes X
$$
である。二つの因子のノルムの積より交換子ノルムは $2$ になる。

:::example[二つのサイトで密度の重みが異なる]
二つのサイトに独立状態を置く。各状態の密度を規格化し、準局所環に入る二サイトの対角観測の値を計算する。

テンソル積密度 $D_1\otimes D_2$ では、添字順 $(1,1),(1,2),(2,1),(2,2)$ で各確率の積を並べる。全体のトレースは二つのトレースの積なので $1$ である。

$$
\begin{aligned}
D_1&=\begin{pmatrix}\frac{1}{4} & 0\\0 & \frac{3}{4}\end{pmatrix}\\
D_2&=\begin{pmatrix}\frac{2}{3} & 0\\0 & \frac{1}{3}\end{pmatrix}\\
D_1\otimes D_2&=\begin{pmatrix}\frac{1}{6} & 0 & 0 & 0\\0 & \frac{1}{12} & 0 & 0\\0 & 0 & \frac{1}{2} & 0\\0 & 0 & 0 & \frac{1}{4}\end{pmatrix}
\end{aligned}
$$

二サイトの射影 $p$ に対する期待値は、選んだ座標の積確率だけを足す。第一サイトが第一座標である事象と交わる確率を、その全確率で割ると条件付き確率を得る。

$$
\begin{aligned}
\operatorname{diag}(p)&=\left( 1, \  0, \  1, \  0\right)\\
\varphi(p)&=\frac{2}{3}\\
\varphi(E_{11}\otimes1\mid p)&=\frac{1}{4}
\end{aligned}
$$

より大きい領域へ移すときは観測量に単位元をテンソルする。外側の密度のトレースは $1$ なのでこの二サイトの期待値は変わらない。従って有限領域全体で求めた期待値に整合性がある。

条件付き状態にすると二サイト間の相関は一般に変わる。元の状態が積状態であることと、射影の結果を固定した後も積状態であることは別である。
:::

:::example[同じ Gibbs 状態の二回の積]
二つのサイトに独立状態を置く。各状態の密度を規格化し、準局所環に入る二サイトの対角観測の値を計算する。

テンソル積密度 $D_1\otimes D_2$ では、添字順 $(1,1),(1,2),(2,1),(2,2)$ で各確率の積を並べる。全体のトレースは二つのトレースの積なので $1$ である。

$$
\begin{aligned}
D_1&=\begin{pmatrix}\frac{1}{3} & 0\\0 & \frac{2}{3}\end{pmatrix}\\
D_2&=\begin{pmatrix}\frac{1}{3} & 0\\0 & \frac{2}{3}\end{pmatrix}\\
D_1\otimes D_2&=\begin{pmatrix}\frac{1}{9} & 0 & 0 & 0\\0 & \frac{2}{9} & 0 & 0\\0 & 0 & \frac{2}{9} & 0\\0 & 0 & 0 & \frac{4}{9}\end{pmatrix}
\end{aligned}
$$

二サイトの射影 $p$ に対する期待値は、選んだ座標の積確率だけを足す。第一サイトが第一座標である事象と交わる確率を、その全確率で割ると条件付き確率を得る。

$$
\begin{aligned}
\operatorname{diag}(p)&=\left( 1, \  0, \  0, \  1\right)\\
\varphi(p)&=\frac{5}{9}\\
\varphi(E_{11}\otimes1\mid p)&=\frac{1}{5}
\end{aligned}
$$

より大きい領域へ移すときは観測量に単位元をテンソルする。外側の密度のトレースは $1$ なのでこの二サイトの期待値は変わらない。従って有限領域全体で求めた期待値に整合性がある。

条件付き状態にすると二サイト間の相関は一般に変わる。元の状態が積状態であることと、射影の結果を固定した後も積状態であることは別である。
:::

:::example[一方だけが非一様な状態]
二つのサイトに独立状態を置く。各状態の密度を規格化し、準局所環に入る二サイトの対角観測の値を計算する。

テンソル積密度 $D_1\otimes D_2$ では、添字順 $(1,1),(1,2),(2,1),(2,2)$ で各確率の積を並べる。全体のトレースは二つのトレースの積なので $1$ である。

$$
\begin{aligned}
D_1&=\begin{pmatrix}\frac{1}{2} & 0\\0 & \frac{1}{2}\end{pmatrix}\\
D_2&=\begin{pmatrix}\frac{1}{5} & 0\\0 & \frac{4}{5}\end{pmatrix}\\
D_1\otimes D_2&=\begin{pmatrix}\frac{1}{10} & 0 & 0 & 0\\0 & \frac{2}{5} & 0 & 0\\0 & 0 & \frac{1}{10} & 0\\0 & 0 & 0 & \frac{2}{5}\end{pmatrix}
\end{aligned}
$$

二サイトの射影 $p$ に対する期待値は、選んだ座標の積確率だけを足す。第一サイトが第一座標である事象と交わる確率を、その全確率で割ると条件付き確率を得る。

$$
\begin{aligned}
\operatorname{diag}(p)&=\left( 0, \  1, \  1, \  0\right)\\
\varphi(p)&=\frac{1}{2}\\
\varphi(E_{11}\otimes1\mid p)&=\frac{4}{5}
\end{aligned}
$$

より大きい領域へ移すときは観測量に単位元をテンソルする。外側の密度のトレースは $1$ なのでこの二サイトの期待値は変わらない。従って有限領域全体で求めた期待値に整合性がある。

条件付き状態にすると二サイト間の相関は一般に変わる。元の状態が積状態であることと、射影の結果を固定した後も積状態であることは別である。
:::

:::example[第一サイトが高確率の区画を選ぶ]
二つのサイトに独立状態を置く。各状態の密度を規格化し、準局所環に入る二サイトの対角観測の値を計算する。

テンソル積密度 $D_1\otimes D_2$ では、添字順 $(1,1),(1,2),(2,1),(2,2)$ で各確率の積を並べる。全体のトレースは二つのトレースの積なので $1$ である。

$$
\begin{aligned}
D_1&=\begin{pmatrix}\frac{2}{3} & 0\\0 & \frac{1}{3}\end{pmatrix}\\
D_2&=\begin{pmatrix}\frac{3}{4} & 0\\0 & \frac{1}{4}\end{pmatrix}\\
D_1\otimes D_2&=\begin{pmatrix}\frac{1}{2} & 0 & 0 & 0\\0 & \frac{1}{6} & 0 & 0\\0 & 0 & \frac{1}{4} & 0\\0 & 0 & 0 & \frac{1}{12}\end{pmatrix}
\end{aligned}
$$

二サイトの射影 $p$ に対する期待値は、選んだ座標の積確率だけを足す。第一サイトが第一座標である事象と交わる確率を、その全確率で割ると条件付き確率を得る。

$$
\begin{aligned}
\operatorname{diag}(p)&=\left( 1, \  1, \  0, \  0\right)\\
\varphi(p)&=\frac{2}{3}\\
\varphi(E_{11}\otimes1\mid p)&=1
\end{aligned}
$$

より大きい領域へ移すときは観測量に単位元をテンソルする。外側の密度のトレースは $1$ なのでこの二サイトの期待値は変わらない。従って有限領域全体で求めた期待値に整合性がある。

条件付き状態にすると二サイト間の相関は一般に変わる。元の状態が積状態であることと、射影の結果を固定した後も積状態であることは別である。
:::

:::example[一様な二サイト状態]
二つのサイトに独立状態を置く。各状態の密度を規格化し、準局所環に入る二サイトの対角観測の値を計算する。

テンソル積密度 $D_1\otimes D_2$ では、添字順 $(1,1),(1,2),(2,1),(2,2)$ で各確率の積を並べる。全体のトレースは二つのトレースの積なので $1$ である。

$$
\begin{aligned}
D_1&=\begin{pmatrix}\frac{1}{2} & 0\\0 & \frac{1}{2}\end{pmatrix}\\
D_2&=\begin{pmatrix}\frac{1}{2} & 0\\0 & \frac{1}{2}\end{pmatrix}\\
D_1\otimes D_2&=\begin{pmatrix}\frac{1}{4} & 0 & 0 & 0\\0 & \frac{1}{4} & 0 & 0\\0 & 0 & \frac{1}{4} & 0\\0 & 0 & 0 & \frac{1}{4}\end{pmatrix}
\end{aligned}
$$

二サイトの射影 $p$ に対する期待値は、選んだ座標の積確率だけを足す。第一サイトが第一座標である事象と交わる確率を、その全確率で割ると条件付き確率を得る。

$$
\begin{aligned}
\operatorname{diag}(p)&=\left( 0, \  1, \  0, \  1\right)\\
\varphi(p)&=\frac{1}{2}\\
\varphi(E_{11}\otimes1\mid p)&=\frac{1}{2}
\end{aligned}
$$

より大きい領域へ移すときは観測量に単位元をテンソルする。外側の密度のトレースは $1$ なのでこの二サイトの期待値は変わらない。従って有限領域全体で求めた期待値に整合性がある。

条件付き状態にすると二サイト間の相関は一般に変わる。元の状態が積状態であることと、射影の結果を固定した後も積状態であることは別である。
:::

## KMS 状態と平衡状態

物理的な逆温度を $\beta=(k_{\rm B}T)^{-1}$ とすると単位はエネルギーの逆数である。有限次元の Gibbs 状態は $\varphi_\beta(a)=\tr(Da)$、$D=e^{-\beta H}/\tr(e^{-\beta H})$ である。物理時間の $\alpha_t$ に対して :ref[def-oa-kms] の帯の高さは $\hbar\beta$ となる。実際、$F(z)=\tr(Da\,e^{izH/\hbar}be^{-izH/\hbar})$ は帯で有界かつ正則である。$z=t+i\hbar\beta$ では $\alpha_z(b)=D\alpha_t(b)D^{-1}$ なので、トレースの巡回性から上側の境界は $\varphi_\beta(\alpha_t(b)a)$ になる。

全行列環ではこの KMS 状態は一意である。状態の密度を $D'$ とする。二つの帯の関数の差が下側境界で零なら、Schwarz の反射原理で下側に正則に延長でき、恒等定理から内部でも零である。従って任意の KMS 状態の帯の関数は、下側境界を持つ行列指数の式に一致する。$t=0$ での上側境界を全 $a,b$ について比較し、トレースの双線形形式の非退化性を使うと $e^{\beta H}D'a e^{-\beta H}=aD'$ である。$a=1$ から $D'$ は $e^{-\beta H}$ と可換であり、元の式の右に $e^{\beta H}$ を掛けると $e^{\beta H}D'$ が全行列と可換だと分かる。従ってこれはスカラー倍の単位元であり、規格化で $D'=D$ となる。無限系でも、時間発展と状態があれば帯の関係を課せる。

:::example[二準位の比]
一サイトのエネルギーが $0,\epsilon$ の場合、励起の確率は $e^{-\beta\epsilon}/(1+e^{-\beta\epsilon})$ である。$\beta\epsilon=\log3$ なら励起確率 $1/4$、平均エネルギー $\epsilon/4$ となる。
:::

:::example[独立なサイトの無限系]
全サイトに同じ一サイト Hamiltonian $h$ があり、異なるサイト間には相互作用がないとする。有限領域の密度行列は $D_\Lambda=\bigotimes_{x\in\Lambda}D$、$D=e^{-\beta h}/\tr(e^{-\beta h})$ である。追加した因子のトレースは $1$ なので局所状態は整合する。その合併上の正値汎関数はノルム $1$ で、準局所環へ一意に延びる。時間発展も各因子で独立なので、有限支えを保ち、合併上から全体へ延びる。

この積状態は全準局所環で KMS 条件を満たす。局所元 $a,b$ は共通の有限領域に入るので Gibbs の証明を適用できる。一般元を局所元 $a_n,b_n$ でノルム近似する。上下の境界で差の絶対値は $\norm{a_n-a_m}\norm{b_n}+\norm{a_m}\norm{b_n-b_m}$ 以下である。有界正則関数に三線定理を適用すると帯全体で同じ上界を得るため、相関の関数は一様 Cauchy となる。その極限は内部で正則かつ閉じた帯で連続で、両端で所定の値になる。逆に全環の KMS 状態を有限領域に制限すれば、その領域が時間発展で不変なので有限次元の Gibbs 状態に一致する。局所元が稠密であることから積状態は一意である。例えば各サイトの励起確率が $1/4$ のとき、異なる二サイトがともに励起する確率は $1/16$ である。
:::

### KMS 条件から時間不変性を導く

KMS 条件の一方に $a=1$ を入れると、$F$ は帯の上下で同じ境界値を持つ。帯の高さを $b_0=\hbar\beta$ とし、$F$ を虚方向へ $ib_0$ ごとに繰り返して平面全体へ延長する。接合する直線で連続であり、三角形の積分を直線の両側に分けると共通辺の積分が逆向きで相殺する。よって Morera の定理で延長は整関数である。もとの帯で有界なので全平面で有界であり、Liouville の定理から定数である。実軸の値を取れば $\varphi(\alpha_t(b))=\varphi(b)$ が従う。

時間発展が点ごとにノルム連続なら、$b$ を Gaussian で時間平均した
$$
b_n=\sqrt{n/\pi}\int_{\R}e^{-nt^2}\alpha_t(b)\,\dd t
$$
は $b$ にノルム収束する。ここで $n$ は固定した時間単位で測る幅の二乗の逆数である。$\alpha_z(b_n)$ は同じ積分で $e^{-nt^2}$ を $e^{-n(t-z)^2}$ に替えると全平面へ延びる。コンパクトな $z$ の集合で Gaussian の減衰が残り、微分できるためである。こうした :term[解析元]{en="analytic element" yomi="かいせきげん"} 全体が稠密である。$b$ が解析元なら、$F(z)=\varphi(a\alpha_z(b))$ は所定の帯で有界である。実移動は等長で、虚方向の区間がコンパクトだからである。その下側境界は既定であり、反射原理による一意性から KMS 条件は
$$
\varphi(a\alpha_{ib_0}(b))=\varphi(ba)
$$
を含む。逆に全 $a$ とすべての解析元 $b$ でこの等式が成り立てば、$b$ を $\alpha_t(b)$ で置き換えて上側境界を得る。一般の $b$ は解析元でノルム近似し、三線定理を用いて帯全体に移す。従って解析元だけを用いて平衡条件を確かめられる。

### 逆温度の単位と相関の帯

物理時間 $t$ に対し $\alpha_t(a)=e^{itH/\hbar}ae^{-itH/\hbar}$ を用いると、上側の帯の幅は $\hbar\beta$ である。$\beta=1/(k_BT)$ はエネルギーの逆数で、$\hbar$ はエネルギーと時間の積であるため、幅の次元は時間になる。例えば $T=300\,\mathrm K$ で $k_B=1.380649\times10^{-23}\,\mathrm{J\,K^{-1}}$ とすると $k_BT=4.14\times10^{-21}\,\mathrm J$ 程度であり、$\hbar=1.05457\times10^{-34}\,\mathrm{J\,s}$ から帯の幅は $2.55\times10^{-14}\,\mathrm s$ 程度になる。この幅は KMS 相関の複素時間方向の尺度である。

:::example[上側で小さい密度へ切り替わる]
有限系で Gibbs 密度を対角基底で指定し、物理時間の発展 $\alpha_t(x)=e^{itH/\hbar}xe^{-itH/\hbar}$ を使う。

$a=E_{ji}$、$b=E_{ij}$ とする。$ab=E_{jj}$、$ba=E_{ii}$ なので、二つの順序は異なる準位の密度を読む。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{1}{4} & 0\\0 & \frac{3}{4}\end{pmatrix}\\
(i,j)&=\left( 1, \  2\right)\\
\varphi(ab)&=\frac{3}{4}\\
\varphi(ba)&=\frac{1}{4}
\end{aligned}
$$

$d_i/d_j=e^{-\beta(E_i-E_j)}$ より、$\alpha_z(b)=(d_i/d_j)^{-iz/(\hbar\beta)}b$ である。そこで下境界では $F(t)$ は $d_j$ にこの位相を掛けたものになる。上端 $i\hbar\beta$ では位相の指数が $1$ になり、比 $d_i/d_j$ が掛かる。

$$
\begin{aligned}
d_i/d_j&=\frac{1}{3}\\
F(0)&=\frac{3}{4}\\
F(i\hbar\beta)&=\frac{1}{4}
\end{aligned}
$$

時間 $t$ を加えても双方に同じ実時間の位相が掛かるので、$F(t+i\hbar\beta)=\varphi(\alpha_t(b)a)$ になる。有限次元の行列単位は線形基底なので、一般の $a,b$ に対しては各項を足せばよい。

下側と上側の値を比較するには、解析延長の位相だけでなく、境界で演算子を掛ける順序の変化も追う必要がある。
:::

:::example[逆の行列単位で上側の値が増える]
有限系で Gibbs 密度を対角基底で指定し、物理時間の発展 $\alpha_t(x)=e^{itH/\hbar}xe^{-itH/\hbar}$ を使う。

$a=E_{ji}$、$b=E_{ij}$ とする。$ab=E_{jj}$、$ba=E_{ii}$ なので、二つの順序は異なる準位の密度を読む。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{3}{4} & 0\\0 & \frac{1}{4}\end{pmatrix}\\
(i,j)&=\left( 1, \  2\right)\\
\varphi(ab)&=\frac{1}{4}\\
\varphi(ba)&=\frac{3}{4}
\end{aligned}
$$

$d_i/d_j=e^{-\beta(E_i-E_j)}$ より、$\alpha_z(b)=(d_i/d_j)^{-iz/(\hbar\beta)}b$ である。そこで下境界では $F(t)$ は $d_j$ にこの位相を掛けたものになる。上端 $i\hbar\beta$ では位相の指数が $1$ になり、比 $d_i/d_j$ が掛かる。

$$
\begin{aligned}
d_i/d_j&=3\\
F(0)&=\frac{1}{4}\\
F(i\hbar\beta)&=\frac{3}{4}
\end{aligned}
$$

時間 $t$ を加えても双方に同じ実時間の位相が掛かるので、$F(t+i\hbar\beta)=\varphi(\alpha_t(b)a)$ になる。有限次元の行列単位は線形基底なので、一般の $a,b$ に対しては各項を足せばよい。

下側と上側の値を比較するには、解析延長の位相だけでなく、境界で演算子を掛ける順序の変化も追う必要がある。
:::

:::example[三準位間の遷移の相関]
有限系で Gibbs 密度を対角基底で指定し、物理時間の発展 $\alpha_t(x)=e^{itH/\hbar}xe^{-itH/\hbar}$ を使う。

$a=E_{ji}$、$b=E_{ij}$ とする。$ab=E_{jj}$、$ba=E_{ii}$ なので、二つの順序は異なる準位の密度を読む。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{1}{7} & 0 & 0\\0 & \frac{2}{7} & 0\\0 & 0 & \frac{4}{7}\end{pmatrix}\\
(i,j)&=\left( 1, \  3\right)\\
\varphi(ab)&=\frac{4}{7}\\
\varphi(ba)&=\frac{1}{7}
\end{aligned}
$$

$d_i/d_j=e^{-\beta(E_i-E_j)}$ より、$\alpha_z(b)=(d_i/d_j)^{-iz/(\hbar\beta)}b$ である。そこで下境界では $F(t)$ は $d_j$ にこの位相を掛けたものになる。上端 $i\hbar\beta$ では位相の指数が $1$ になり、比 $d_i/d_j$ が掛かる。

$$
\begin{aligned}
d_i/d_j&=\frac{1}{4}\\
F(0)&=\frac{4}{7}\\
F(i\hbar\beta)&=\frac{1}{7}
\end{aligned}
$$

時間 $t$ を加えても双方に同じ実時間の位相が掛かるので、$F(t+i\hbar\beta)=\varphi(\alpha_t(b)a)$ になる。有限次元の行列単位は線形基底なので、一般の $a,b$ に対しては各項を足せばよい。

下側と上側の値を比較するには、解析延長の位相だけでなく、境界で演算子を掛ける順序の変化も追う必要がある。
:::

:::example[同じ重みを持つ二つの準位]
有限系で Gibbs 密度を対角基底で指定し、物理時間の発展 $\alpha_t(x)=e^{itH/\hbar}xe^{-itH/\hbar}$ を使う。

$a=E_{ji}$、$b=E_{ij}$ とする。$ab=E_{jj}$、$ba=E_{ii}$ なので、二つの順序は異なる準位の密度を読む。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{3}{8} & 0 & 0\\0 & \frac{3}{8} & 0\\0 & 0 & \frac{1}{4}\end{pmatrix}\\
(i,j)&=\left( 1, \  2\right)\\
\varphi(ab)&=\frac{3}{8}\\
\varphi(ba)&=\frac{3}{8}
\end{aligned}
$$

$d_i/d_j=e^{-\beta(E_i-E_j)}$ より、$\alpha_z(b)=(d_i/d_j)^{-iz/(\hbar\beta)}b$ である。そこで下境界では $F(t)$ は $d_j$ にこの位相を掛けたものになる。上端 $i\hbar\beta$ では位相の指数が $1$ になり、比 $d_i/d_j$ が掛かる。

$$
\begin{aligned}
d_i/d_j&=1\\
F(0)&=\frac{3}{8}\\
F(i\hbar\beta)&=\frac{3}{8}
\end{aligned}
$$

時間 $t$ を加えても双方に同じ実時間の位相が掛かるので、$F(t+i\hbar\beta)=\varphi(\alpha_t(b)a)$ になる。有限次元の行列単位は線形基底なので、一般の $a,b$ に対しては各項を足せばよい。

下側と上側の値を比較するには、解析延長の位相だけでなく、境界で演算子を掛ける順序の変化も追う必要がある。
:::

:::example[順序を入れ替えた相関]
有限系で Gibbs 密度を対角基底で指定し、物理時間の発展 $\alpha_t(x)=e^{itH/\hbar}xe^{-itH/\hbar}$ を使う。

$a=E_{ji}$、$b=E_{ij}$ とする。$ab=E_{jj}$、$ba=E_{ii}$ なので、二つの順序は異なる準位の密度を読む。

$$
\begin{aligned}
D&=\begin{pmatrix}\frac{1}{3} & 0 & 0\\0 & \frac{1}{6} & 0\\0 & 0 & \frac{1}{2}\end{pmatrix}\\
(i,j)&=\left( 3, \  2\right)\\
\varphi(ab)&=\frac{1}{6}\\
\varphi(ba)&=\frac{1}{2}
\end{aligned}
$$

$d_i/d_j=e^{-\beta(E_i-E_j)}$ より、$\alpha_z(b)=(d_i/d_j)^{-iz/(\hbar\beta)}b$ である。そこで下境界では $F(t)$ は $d_j$ にこの位相を掛けたものになる。上端 $i\hbar\beta$ では位相の指数が $1$ になり、比 $d_i/d_j$ が掛かる。

$$
\begin{aligned}
d_i/d_j&=3\\
F(0)&=\frac{1}{6}\\
F(i\hbar\beta)&=\frac{1}{2}
\end{aligned}
$$

時間 $t$ を加えても双方に同じ実時間の位相が掛かるので、$F(t+i\hbar\beta)=\varphi(\alpha_t(b)a)$ になる。有限次元の行列単位は線形基底なので、一般の $a,b$ に対しては各項を足せばよい。

下側と上側の値を比較するには、解析延長の位相だけでなく、境界で演算子を掛ける順序の変化も追う必要がある。
:::

## 時間発展と Lieb–Robinson 評価（概観）

格子距離を $d$、集合の直径を $\operatorname{diam}Z$ とする。相互作用は有限範囲とし、ある $\mu>0$ に対し $C_\mu=\sup_x\sum_{Z\ni x}\abs Z\norm{\Phi(Z)}e^{\mu\operatorname{diam}Z}<\infty$ とする。$C_\mu$ はエネルギーの単位を持ち、$\mu d$ は無次元である。

:::proposition[局所作用素の交換子評価]{#prop-oa-lieb-robinson}
互いに素な $X,Y\subset\Lambda$ と $a\in\mathcal A_X,b\in\mathcal A_Y$ に対し
$$
\norm{[\alpha_t^\Lambda(a),b]}
\le 2\norm a\norm b\,\abs X\, e^{-\mu d(X,Y)}
\big(e^{2C_\mu\abs t/\hbar}-1\big)
$$
となる。
:::

:::proof
$t\ge0$ とし $D_B(X,t)$ を $a\in\mathcal A_X$ の単位球での左辺の上限とする。$a$ と可換でない Hamiltonian の項は $Z\cap X\ne\varnothing$ に限られる。それらの和を $H_X$ とすると $f(t)=[\alpha_t^\Lambda(a),b]$ の微分は
$$
f'(t)=\frac i\hbar[\alpha_t^\Lambda(H_X),f(t)]
-\frac i\hbar[\alpha_t^\Lambda(a),[\alpha_t^\Lambda(H_X),b]]
$$
となる。第一項は時間依存ユニタリによる共役を生成してノルムを保つため、変分定数法で
$$
D_B(X,t)\le D_B(X,0)
+\frac2\hbar\sum_{Z\cap X\ne\varnothing}\norm{\Phi(Z)}
\int_0^t D_B(Z,s)\,\dd s
$$
を得る。$D_B(X,0)$ は $X\cap Y=\varnothing$ なら零で、一般に $2\norm b$ 以下である。反復で非零項は $X$ から $Y$ まで、隣接する集合が重なる $n$ 個の相互作用集合の鎖になる。その直径の和は $d(X,Y)$ 以上である。各重なりの点を指定して鎖を上から数えると、初めに $\abs X$ 通り、各段で重み $\abs Z\norm{\Phi(Z)}e^{\mu\operatorname{diam}Z}$ の和が $C_\mu$ 以下である。各鎖の積の和は従って $\abs X C_\mu^ne^{-\mu d(X,Y)}$ 以下である。残差は $D_B\le2\norm b$ と同じ重みの評価と反復積分 $t^n/n!$ により零に収束する。$n\ge1$ の指数級数を足し合わせる。$t<0$ は $-\Phi$ の発展に同じ評価を適用する。
:::

$\Lambda\subset\Lambda'$ とし、まず $t\ge0$ を取る。$s\mapsto \alpha_{t-s}^{\Lambda'}(\alpha_s^\Lambda(a))$ を微分すると、二つの Hamiltonian の差との交換子が現れる。$s=0,t$ の間で積分し、$\alpha^{\Lambda'}$ の等長性を使うと二発展の差は
$$
\frac1\hbar\int_0^t\sum_{\substack{Z\subset\Lambda'\\Z\not\subset\Lambda}}
\norm{[\Phi(Z),\alpha_s^\Lambda(a)]}\,\dd s
$$
以下である。負時刻は $-\Phi$ に置き換える。相互作用の範囲を $R$ とし、$\Lambda$ が $X$ の $R$ 近傍を含む段階を取る。和に寄与する $Z$ は $X$ と交わらず、:ref[prop-oa-lieb-robinson] を $\Lambda'$ に零延長した $\Lambda$ 内の相互作用に適用できる。各 $Z$ に一点 $y\in Z\setminus\Lambda$ を取ると $d(X,Z)\ge d(X,y)-R$ である。従って距離の因子を含む和は
$$
\sum_{\substack{Z\subset\Lambda'\\Z\not\subset\Lambda}}\norm{\Phi(Z)}e^{-\mu d(X,Z)}
\le e^{\mu R} C_\mu\sum_{y\notin\Lambda}e^{-\mu d(X,y)}
$$
以下である。最後の級数は $\Z^d$ の球の点数が多項式で抑えられるため収束し、その尾部は $\Lambda$ が広がると零に行く。この評価は時間のコンパクト区間で一様なので、局所元の発展はノルムで極限を持つ。有限体積発展の等長性により全準局所環へ延び、任意の元で時間についてノルム連続である。群法則は、$\alpha_s^\Lambda(a)$ を先に極限 $\alpha_s(a)$ と比較してから $\alpha_t^\Lambda$ を通すと、誤差が等長性で増えないことから従う。従って $\alpha_{-t}$ が逆写像となる。

同じ上界は熱力学極限にも移る。$v_\mu=2C_\mu/(\hbar\mu)$ と置くと、指数の部分は $e^{-\mu(d(X,Y)-v_\mu|t|)}$ で抑えられる。従って距離が $v_\mu|t|$ を大きく超える領域では交換子が指数的に小さい。$v_\mu$ の単位は格子距離を測る長さを時間で割ったものである。

### 有限体積 Gibbs 状態の極限

有限体積の Gibbs 状態に体積外の規格化行列トレースの積を付け、全準局所環の状態 $\varphi_\Lambda$ とする。有限体積の発展は外で恒等であり、共通の有限領域で見れば密度行列はその領域における $e^{-\beta H_\Lambda}$ の規格化である。従って局所元に対して有限次元での KMS の証明が適用でき、ノルム近似で全環にも延びる。

状態空間の弱 * コンパクト性により蓄積点 $\varphi$ がある。前節の Gaussian 平均 $b_n$ を $\alpha^\Lambda$ で作ったものを $b_{n,\Lambda}$ と書く。固定した $n$ では $\alpha^\Lambda$ の実時間でのノルム収束と Gaussian の可積分な上界により、
$$
b_{n,\Lambda}\longrightarrow b_n,\qquad
\alpha^\Lambda_{ib_0}(b_{n,\Lambda})\longrightarrow\alpha_{ib_0}(b_n)
$$
がノルムで成り立つ。二番目にも、虚方向へ移した Gaussian の絶対値は $e^{n b_0^2}e^{-nt^2}$ で抑えられる。従って $\varphi_\Lambda(a\alpha_{ib_0}^\Lambda(b_{n,\Lambda}))=\varphi_\Lambda(b_{n,\Lambda}a)$ の弱 * 極限を取れば、$\varphi$ に対する解析元の関係が得られる。実時間の移動は Gaussian 平均と可換なので、この関係は $\alpha_t(b_n)$ にも適用できる。$b_n\to b$ と三線定理により、$\varphi$ は極限の時間発展に関する KMS 状態となる。この構成で存在が分かり、独立なサイトの例では局所密度の一意性から蓄積点も一つに定まる。

### 交換子の最初の非零項

有限体積では $\alpha_t(a)=a+(it/\hbar)[H,a]+((it/\hbar)^2/2)[H,[H,a]]+\cdots$ がノルム収束する。二サイトで $H=JZ\otimes Z$、$a=X\otimes I$、$b=I\otimes X$ とする。$[a,b]=0$ だが
$$
[[H,a],b]
=J[Z,X]\otimes[Z,X]
=4J(E_{12}-E_{21})\otimes(E_{12}-E_{21})
$$
である。従って短時間の交換子のノルムは $4|J|\,|t|/\hbar+O(t^2)$ となる。初めから作用が伝わっているのではなく、二サイトを同時に含む相互作用が時間展開の一次項で寄与している。$J$ がエネルギーなら係数 $|J|t/\hbar$ は無次元になる。

:::example[短時間に対応する八分の一回転]
二サイトの Ising 相互作用 $H=JZ\otimes Z$ だけを残し、$a=X\otimes I$ と離れたサイトの観測量 $b$ の交換子を近似せずに求める。

$A=Z\otimes Z$ が $A^2=I$ を満たすため、$\theta=Jt/\hbar$ を用い $e^{itH/\hbar}=\cos\theta\,I+i\sin\theta\,A$ と展開できる。$aA=-Aa$ を使うと有限二項に整理できる。

$$
\begin{aligned}
\theta&=\frac{\pi}{8}\\
\cos(2\theta)&=\frac{\sqrt{2}}{2}\\
\sin(2\theta)&=\frac{\sqrt{2}}{2}
\end{aligned}
$$

$Y=\begin{pmatrix}0&-i\\ i&0\end{pmatrix}$ と置くと $ZX=iY$ である。従って $\alpha_t(a)=\cos(2\theta)X\otimes I-\sin(2\theta)Y\otimes Z$ である。第二サイトでテストする観測量の軸を指定し、その交換子を求める。

$$
\begin{aligned}
b&=\begin{pmatrix}0 & 1 & 0 & 0\\1 & 0 & 0 & 0\\0 & 0 & 0 & 1\\0 & 0 & 1 & 0\end{pmatrix}\\
[\alpha_t(a),b]&=\begin{pmatrix}0 & 0 & 0 & \sqrt{2}i\\0 & 0 & -\sqrt{2}i & 0\\0 & -\sqrt{2}i & 0 & 0\\\sqrt{2}i & 0 & 0 & 0\end{pmatrix}
\end{aligned}
$$

この行列の随伴と自身の積はスカラー倍の恒等行列である。そのスカラーの平方根を取ると作用素ノルムになる。

$$
\begin{aligned}
\norm{[\alpha_t(a),b]}&=\sqrt{2}
\end{aligned}
$$

時間の尺度は無次元の $Jt/\hbar$ で決まる。短時間評価での線形項は、この正確な式で $\sin(2\theta)$ を原点の近くで展開して得られる。
:::

:::example[二サイト交換子が最大値に達する]
二サイトの Ising 相互作用 $H=JZ\otimes Z$ だけを残し、$a=X\otimes I$ と離れたサイトの観測量 $b$ の交換子を近似せずに求める。

$A=Z\otimes Z$ が $A^2=I$ を満たすため、$\theta=Jt/\hbar$ を用い $e^{itH/\hbar}=\cos\theta\,I+i\sin\theta\,A$ と展開できる。$aA=-Aa$ を使うと有限二項に整理できる。

$$
\begin{aligned}
\theta&=\frac{\pi}{4}\\
\cos(2\theta)&=0\\
\sin(2\theta)&=1
\end{aligned}
$$

$Y=\begin{pmatrix}0&-i\\ i&0\end{pmatrix}$ と置くと $ZX=iY$ である。従って $\alpha_t(a)=\cos(2\theta)X\otimes I-\sin(2\theta)Y\otimes Z$ である。第二サイトでテストする観測量の軸を指定し、その交換子を求める。

$$
\begin{aligned}
b&=\begin{pmatrix}0 & 1 & 0 & 0\\1 & 0 & 0 & 0\\0 & 0 & 0 & 1\\0 & 0 & 1 & 0\end{pmatrix}\\
[\alpha_t(a),b]&=\begin{pmatrix}0 & 0 & 0 & 2 i\\0 & 0 & - 2 i & 0\\0 & - 2 i & 0 & 0\\2 i & 0 & 0 & 0\end{pmatrix}
\end{aligned}
$$

この行列の随伴と自身の積はスカラー倍の恒等行列である。そのスカラーの平方根を取ると作用素ノルムになる。

$$
\begin{aligned}
\norm{[\alpha_t(a),b]}&=2
\end{aligned}
$$

時間の尺度は無次元の $Jt/\hbar$ で決まる。短時間評価での線形項は、この正確な式で $\sin(2\theta)$ を原点の近くで展開して得られる。
:::

:::example[一時的に交換子が再び零になる]
二サイトの Ising 相互作用 $H=JZ\otimes Z$ だけを残し、$a=X\otimes I$ と離れたサイトの観測量 $b$ の交換子を近似せずに求める。

$A=Z\otimes Z$ が $A^2=I$ を満たすため、$\theta=Jt/\hbar$ を用い $e^{itH/\hbar}=\cos\theta\,I+i\sin\theta\,A$ と展開できる。$aA=-Aa$ を使うと有限二項に整理できる。

$$
\begin{aligned}
\theta&=\frac{\pi}{2}\\
\cos(2\theta)&=-1\\
\sin(2\theta)&=0
\end{aligned}
$$

$Y=\begin{pmatrix}0&-i\\ i&0\end{pmatrix}$ と置くと $ZX=iY$ である。従って $\alpha_t(a)=\cos(2\theta)X\otimes I-\sin(2\theta)Y\otimes Z$ である。第二サイトでテストする観測量の軸を指定し、その交換子を求める。

$$
\begin{aligned}
b&=\begin{pmatrix}0 & 1 & 0 & 0\\1 & 0 & 0 & 0\\0 & 0 & 0 & 1\\0 & 0 & 1 & 0\end{pmatrix}\\
[\alpha_t(a),b]&=\begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{pmatrix}
\end{aligned}
$$

この行列の随伴と自身の積はスカラー倍の恒等行列である。そのスカラーの平方根を取ると作用素ノルムになる。

$$
\begin{aligned}
\norm{[\alpha_t(a),b]}&=0
\end{aligned}
$$

時間の尺度は無次元の $Jt/\hbar$ で決まる。短時間評価での線形項は、この正確な式で $\sin(2\theta)$ を原点の近くで展開して得られる。
:::

:::example[負時刻でも同じ大きさ]
二サイトの Ising 相互作用 $H=JZ\otimes Z$ だけを残し、$a=X\otimes I$ と離れたサイトの観測量 $b$ の交換子を近似せずに求める。

$A=Z\otimes Z$ が $A^2=I$ を満たすため、$\theta=Jt/\hbar$ を用い $e^{itH/\hbar}=\cos\theta\,I+i\sin\theta\,A$ と展開できる。$aA=-Aa$ を使うと有限二項に整理できる。

$$
\begin{aligned}
\theta&=- \frac{\pi}{8}\\
\cos(2\theta)&=\frac{\sqrt{2}}{2}\\
\sin(2\theta)&=- \frac{\sqrt{2}}{2}
\end{aligned}
$$

$Y=\begin{pmatrix}0&-i\\ i&0\end{pmatrix}$ と置くと $ZX=iY$ である。従って $\alpha_t(a)=\cos(2\theta)X\otimes I-\sin(2\theta)Y\otimes Z$ である。第二サイトでテストする観測量の軸を指定し、その交換子を求める。

$$
\begin{aligned}
b&=\begin{pmatrix}0 & 1 & 0 & 0\\1 & 0 & 0 & 0\\0 & 0 & 0 & 1\\0 & 0 & 1 & 0\end{pmatrix}\\
[\alpha_t(a),b]&=\begin{pmatrix}0 & 0 & 0 & -\sqrt{2}i\\0 & 0 & \sqrt{2}i & 0\\0 & \sqrt{2}i & 0 & 0\\-\sqrt{2}i & 0 & 0 & 0\end{pmatrix}
\end{aligned}
$$

この行列の随伴と自身の積はスカラー倍の恒等行列である。そのスカラーの平方根を取ると作用素ノルムになる。

$$
\begin{aligned}
\norm{[\alpha_t(a),b]}&=\sqrt{2}
\end{aligned}
$$

時間の尺度は無次元の $Jt/\hbar$ で決まる。短時間評価での線形項は、この正確な式で $\sin(2\theta)$ を原点の近くで展開して得られる。
:::

:::example[第二の観測量を同じ軸に替える]
二サイトの Ising 相互作用 $H=JZ\otimes Z$ だけを残し、$a=X\otimes I$ と離れたサイトの観測量 $b$ の交換子を近似せずに求める。

$A=Z\otimes Z$ が $A^2=I$ を満たすため、$\theta=Jt/\hbar$ を用い $e^{itH/\hbar}=\cos\theta\,I+i\sin\theta\,A$ と展開できる。$aA=-Aa$ を使うと有限二項に整理できる。

$$
\begin{aligned}
\theta&=\frac{\pi}{6}\\
\cos(2\theta)&=\frac{1}{2}\\
\sin(2\theta)&=\frac{\sqrt{3}}{2}
\end{aligned}
$$

$Y=\begin{pmatrix}0&-i\\ i&0\end{pmatrix}$ と置くと $ZX=iY$ である。従って $\alpha_t(a)=\cos(2\theta)X\otimes I-\sin(2\theta)Y\otimes Z$ である。第二サイトでテストする観測量の軸を指定し、その交換子を求める。

$$
\begin{aligned}
b&=\begin{pmatrix}1 & 0 & 0 & 0\\0 & -1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & -1\end{pmatrix}\\
[\alpha_t(a),b]&=\begin{pmatrix}0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\\0 & 0 & 0 & 0\end{pmatrix}
\end{aligned}
$$

この行列の随伴と自身の積はスカラー倍の恒等行列である。そのスカラーの平方根を取ると作用素ノルムになる。

$$
\begin{aligned}
\norm{[\alpha_t(a),b]}&=0
\end{aligned}
$$

時間の尺度は無次元の $Jt/\hbar$ で決まる。短時間評価での線形項は、この正確な式で $\sin(2\theta)$ を原点の近くで展開して得られる。
:::

## 演習問題

::::exercise{level=basic}
$W(v)W(-v)$ を求めよ。
:::solution
交代形式なので $\sigma(v,-v)=0$ である。Weyl 関係により積は $W(0)=1$ となる。
:::
::::
::::exercise{level=basic}
一粒子の次元が $3$ のとき、二粒子空間と全 Fock 空間の次元を求めよ。
:::solution
二粒子の外積基底は三つの基底から異なる二つを取った組で $3$ 個である。全体には各基底の占有の有無があるので $2^3=8$ 次元である。
:::
::::
::::exercise{level=basic}
$c(e_1)(e_2\wedge e_1)$ を求めよ。ただし $e_1,e_2$ は正規直交とする。
:::solution
第 1 因子の係数は内積が零であり、第 2 因子だけが残る。符号は $(-1)^{2-1}=-1$ なので $-e_2$ である。
:::
::::
::::exercise{level=basic}
二準位系で $e^{-\beta\epsilon}=1/4$ の場合の励起確率を求めよ。
:::solution
二つの重みは $1,1/4$ なので、励起確率は $(1/4)/(1+1/4)=1/5$ である。
:::
::::
::::exercise{level=advanced}
Weyl 関係から $W(v)W(w)W(v)^*W(w)^*=e^{-i\sigma(v,w)/\hbar}1$ を導け。
:::solution
積を逆順にした位相では $\sigma(w,v)=-\sigma(v,w)$ となる。他はどちらも $W(v+w)$ なので、位相の比が $e^{-i\sigma(v,w)/\hbar}$ である。右から逆元を順に掛けると主張を得る。
:::
::::
::::exercise{level=advanced}
二つの局所作用素の支えが交わらなければ交換子が零になることを、共通の有限領域で示せ。
:::solution
共通領域 $X\cup Y$ のテンソル分解では $a=a_X\otimes1$、$b=1\otimes b_Y$ と表せる。両順序の積が $a_X\otimes b_Y$ となるため交換子は零である。
:::
::::
::::exercise{level=advanced}
有限系の Gibbs 状態が時間発展で不変であることを示せ。
:::solution
$D$ は $H$ の関数なので $e^{itH/\hbar}$ と可換である。$\tr(D e^{itH/\hbar}a e^{-itH/\hbar})=\tr(e^{-itH/\hbar}D e^{itH/\hbar}a)=\tr(Da)$ となる。
:::
::::
::::exercise{level=advanced}
サイト間相互作用がなく $\Phi(Z)=0$（$\abs Z\ne1$）のとき、$\alpha_t^\Lambda(a)$ の支えが時間で広がらないことを示せ。
:::solution
Hamiltonian は互いに可換な一サイト項の和なので指数は一サイトユニタリの積となる。$X$ の外のユニタリは $a\in\mathcal A_X$ と可換で共役から消える。残りの共役は $\mathcal A_X$ の中で行われるため支えは $X$ に留まる。
:::
::::

## 確認テスト

::::quiz{answer="16" tol="1e-6"}
スピン $1/2$ の四サイト領域の局所全行列環の行列サイズを答えよ。
:::explanation
各サイトの Hilbert 空間は二次元で、テンソル積の次元は $2^4=16$ である。局所環は $M_{16}$ となる。
:::
::::
::::quiz
Stone–von Neumann の一意性をここで適用する表現の条件はどれか。

- [x] 有限自由度の非退化な Weyl 関係を満たす。
- [x] 正則で既約である。
- [ ] すべての Weyl ユニタリがノルム連続である。

:::explanation
有限次元位相空間での強連続性から積分表現を作り、既約なら多重度が $1$ となる。Schrödinger 表現そのものでも非零の位相乗算と単位元の距離が $2$ なので、ノルム連続性は条件ではない。
:::
::::
::::quiz{answer="0.25" tol="1e-6"}
エネルギー間隔 $\epsilon$ の二準位系で $\beta\epsilon=\log3$ のとき、励起確率を答えよ。
:::explanation
基底と励起の重みは $1$ と $1/3$ なので、励起確率は $(1/3)/(4/3)=1/4$ である。
:::
::::
::::quiz
二つの有限な支えが離れているとき、Lieb–Robinson 評価が与えるものはどれか。

- [x] 固定時刻で支えの距離に対して指数的に小さくなる交換子の上界
- [ ] 全ての非零時刻で交換子が厳密に零になること
- [ ] 有限体積のサイト数に比例する下界

:::explanation
上界には $e^{-\mu d(X,Y)}$ が掛かる。一方、右辺は一般には正なので交換子の厳密な消滅は導かない。相互作用の一様な和で定数を取るため、評価の定数は全体積のサイト数を必要としない。
:::
::::

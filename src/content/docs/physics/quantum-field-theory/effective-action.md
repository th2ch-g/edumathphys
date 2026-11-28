---
title: "有効作用と生成汎関数"
description: "連結生成汎関数の Legendre 変換を使って 1PI 頂点、量子補正、対称性の恒等式を求める。"
notation:
  - symbol: 'W_E[J]'
    meaning: Euclid の連結生成汎関数
  - symbol: '\Phi'
    meaning: 源の下での場の平均
  - symbol: '\Gamma_E[\Phi]'
    meaning: Euclid の 1PI 有効作用
---

## この章の目標

- 連結相関と 1PI 頂点の生成汎関数を区別する。
- Gaussian 揺らぎから一ループの有効作用を求める。
- 対称性を Legendre 変換後の恒等式に移す。

主に正測度を使える Euclid の正則化から始める。通常は $\hbar=c=1$ とし、ループを数えるときだけ $\hbar$ を明示する。実時間に戻す箇所では計量 $\eta=\diag(+1,-1,-1,-1)$ とする。

## 連結 Green 関数の生成汎関数

### 積を除く対数

本章では Euclid 時間を使い、$Z_E[J]=\int\mathcal D\phi\,e^{-S_E[\phi]+J\phi}$、$W_E[J]=\log Z_E[J]$ とする。$J\phi=\int\dd^4x\,J(x)\phi(x)$ の略記である。測度を正に取れる正則化した実場では $W_E^{(1)}(x)=\langle\phi(x)\rangle_J$、$W_E^{(2)}(x,y)=\langle\phi(x)\phi(y)\rangle_J-\langle\phi(x)\rangle_J\langle\phi(y)\rangle_J$ である。実際、規格化された期待値を微分すると分母の微分が第二項を作る。グラフでは連結成分の積が指数をなすので、対数の微分は :term[連結 Green 関数]{en="connected Green function" yomi="れんけつぐりーんかんすう"} を与える。例えば中心化した一変数 $X$ では四次の連結量は $\langle X^4\rangle-3\langle X^2\rangle^2$。係数三は同じ対への三分割から来る。

### 一様源の累積量を有限系で求める

二状態 $X=\pm v$ を等しい確率で取る模型で $Z(J)/Z(0)=\cosh(Jv)$。従って $W(J)-W(0)=(Jv)^2/2-(Jv)^4/12+\order{J^6}$ である。二次累積量は $v^2$、四次は $-2v^4$。モーメントでは $\langle X^4\rangle=v^4$ だから、四次の連結量とは符号も異なる。二つの独立な同じ変数 $X_1,X_2$ を足せば $\cosh^2(Jv)$ となるので累積量は二倍である。一方、四次モーメントは $2v^4+6v^4=8v^4$。対数の操作は独立な成分にまたがる $6v^4$ の寄与を除く。

### 累積量を分割の和から求める

正則化した $Z_E[J]/Z_E[0]$ をある固定した背景の周りで展開する。$W_E$ の $n$ 階微分を $C_n$ と置くと、$Z_E/Z_E[0]=\exp(W_E-W_E[0])$ の微分は外点を集合分割した和になる。実際、$k$ 個の指数因子を選ぶと $1/k!$ が付き、$n$ 個のラベルを $k$ 個のブロックへ割り当てる際のブロック順の $k!$ がそれを消す。よって $n$ 点モーメントは、全ての集合分割 $\mathcal P$ に対し $\prod_{B\in\mathcal P}C_{|B|}(x_B)$ を足したもの。この式で一ブロックの項以外を左へ移せば、連結相関を次数ごとに再帰的に取り出せる。

:::example[奇数点が消える場の六点]
反転対称な一変数で分散を $G$ とする。六点の部分分割は $4+2$ が $\binom62=15$ 通り、$2+2+2$ が $6!/(2^3\,3!)=15$ 通り。よって $\langle X^6\rangle=C_6+15 C_4 G+15G^3$。例えば $X=\pm v$ なら $C_4=-2v^4$、$G=v^2$ なので $C_6=16v^6$。Gaussian なら $C_4=C_6=0$ であり、六点に残る $15G^3$ は独立な二点だけで決まる。
:::

:::example[源を局在させた二点行列]
有限格子で $Z_E[J]\propto e^{J^\mathsf T D J/2}$ とし、$J$ を一つの格子点 $a$ にだけ入れる。すると $\partial_{J_a}W$ は $D_{aa}J_a$、別の点 $b$ の応答は $\partial_{J_b}W=D_{ba}J_a$。$\partial_{J_b}\partial_{J_a}W=D_{ba}$ は源の値に依存しない。一方、$Z_E$ の微分には $\langle\phi_a\rangle_J\langle\phi_b\rangle_J$ が含まれるので、非零源で分散行列を求めるときは対数を使う必要がある。
:::

### 背景の平均値を中心に連結相関を測る

非零の $J_0$ で場の平均を $\Phi_0$ とする。追加の源 $j$ に対して
$$
\log\left\langle e^{j\cdot(\phi-\Phi_0)}\right\rangle_{J_0}
=W_E[J_0+j]-W_E[J_0]-j\cdot\Phi_0 .
$$
従って中心化は一次の連結相関を消すだけで、二次以上の微分は同じ背景 $J_0$ での連結相関である。特に、二次の共分散は $\Phi_0$ をどこで測るかによって変わり得る。二状態の $X=\pm v$ では $\langle X\rangle_J=v\tanh(Jv)$、分散は $v^2\operatorname{sech}^2(Jv)$。三次は二次を $J$ で微分して $-2v^3\tanh(Jv)\operatorname{sech}^2(Jv)$。例えば $Jv=\operatorname{arctanh}(1/2)$ なら分散 $3v^2/4$、三次 $-3v^3/4$。源を加えたことで測度の重みが偏り、もとの反転対称な模型でも非零の三次が生じる。

### 同じ一様源でも連結な体積の次数を区別する

格子の体積を $V$、一様源を $J$、全場の和を $M=\sum_x\phi_x$ とする。$\partial_J W=\langle M\rangle$、$\partial_J^2 W=\operatorname{Var}(M)$。場の平均密度 $m=M/V$ の分散は $\partial_J^2 W/V^2$、源に対する密度の応答は $\partial_J^2 W/V$ であり、規格化に $V$ の一冪の違いがある。翻訳不変な純相では $\operatorname{Var}(M)=V\sum_r G_c(r)$。$G_c$ が十分速く減衰して和が有限なら、密度の分散は $1/V$ に縮む。一方、二相の混合で $M=\pm Vv$ を等確率で取れば分散は $V^2v^2$。この場合は密度の分散が有限に残り、密度応答は体積とともに増える。個々の相における有限な感受率と、混合の重みが源に応答して作る増大は同じものではない。

### 複数の応答を一つの源の微分で取り出す

二変数 $X,Y$ に源 $J_X,J_Y$ を付け、$W=\log\langle e^{J_X X+J_Y Y}\rangle$ とする。二つの順序の微分が交換できる有限の領域では $\partial_{J_X}\partial_{J_Y}W$ は対称で、$X,Y$ の共分散。三次の $\partial_{J_X}^2\partial_{J_Y}W$ は中心化した $\langle(X-\langle X\rangle)^2(Y-\langle Y\rangle)\rangle$。例えば $X$ が分散 $D$ の中心化 Gaussian で $Y=X^2$ なら、零源で二次の交差共分散は零だが三次は $3D^2-D^2=2D^2$。独立な二変数なら全ての混合累積量は零。交差の二次が消えるだけでは、確率分布全体での独立性までは得られない。

この例は生成汎関数を直接評価しても検査できる。$J_Y<1/(2D)$ なら平方完成により
$$
W(J_X,J_Y)=-\tfrac12\log(1-2DJ_Y)
+\frac{D J_X^2}{2(1-2DJ_Y)} .
$$
ここでは $J=0$ の規格化で定数を除いた。$J_X^2J_Y$ の係数は $D^2$ なので、二回の $J_X$ 微分により $2D^2$。$J_Y$ を臨界値へ近づけると Gaussian の二次係数が零へ向かって共分散が増大する。複合演算子に源を入れると測度の安定条件も変わるため、元の線形源の領域をそのまま使えるとは限らない。

## Legendre 変換と 1PI 有効作用

$\Phi=\delta W_E/\delta J$ とする。$W_E^{(2)}$ が可逆な領域で $J=J[\Phi]$ を解き、
$$
\Gamma_E[\Phi]=J\Phi-W_E[J],\qquad
\frac{\delta\Gamma_E}{\delta\Phi}=J,\qquad
\Gamma_E^{(2)}W_E^{(2)}=1
\label{eq-qft-legendre-effective}
$$
と定める。$\delta J$ の係数は $\Phi-W_E^{(1)}=0$ なので第一微分が $J$。その微分と $W_E^{(2)}=\delta\Phi/\delta J$ に連鎖律を適用すれば逆の関係が出る。逆が失われる領域では Legendre–Fenchel 変換の上限で定める。

:term[一粒子既約]{en="one-particle irreducible" yomi="いちりゅうしきやく"} は内部線一本を切っても連結であることをいう。連結グラフを切断可能な線で全て切ると、既約成分を頂点、切った線を辺とする木になる。閉路があればその線は切断可能でないためである。完全二点関数でこの木をつなぐと、例えば微分 $\delta W_E^{(2)}=-W_E^{(2)}(\delta\Gamma_E^{(2)})W_E^{(2)}$ に $\delta/\delta J=W_E^{(2)}\delta/\delta\Phi$ を入れて
$$
W^{(3)}_{abc}=-W^{(2)}_{ai}W^{(2)}_{bj}W^{(2)}_{ck}\Gamma^{(3)}_{ijk}
$$
となる。添字は位置の積分を含む。次の微分では四脚頂点の項と二つの三脚頂点を一本の完全伝播関数でつなぐ三通りが出る。帰納的に、微分は頂点の脚を一つ増やすか辺に三脚頂点を挿入する。これで $\Gamma$ の頂点と完全二点関数が連結木を再生するため、$\Gamma$ の微分が外脚を除いた 1PI 頂点を与える。例えば自由場では $\Gamma_E=S_E+$ 定数である。

### 四点の完全な逆変換

単一変数で $W''=G$、$\Gamma''=1/G$ とすると、$W'''=-G^3\Gamma'''$。さらに $\dd/\dd J=G\dd/\dd\Phi$ であるから
$$
W''''=-G^4\Gamma''''+3G^5(\Gamma''')^2 .
$$
例えば $\Phi=0$ の近くで $W(J)=aJ^2/2+bJ^4/24+\cdots$、$a>0$ とすると $G(0)=a$、奇の微分は零、従って $\Gamma''''(0)=-b/a^4$。これを直接確かめるには $\Phi=aJ+bJ^3/6+\cdots$ を逆に解き $J=\Phi/a-b\Phi^3/(6a^4)+\cdots$ とし、$\Gamma'=J$ を積分すればよい。場では上式の係数三を、外脚を二組へ分ける三通りの添字縮約へ置き換える。

### 正則化された凸性と逆写像の条件

Euclid の正測度に対して $\delta J$ を任意に選ぶと、
$$
\delta J\cdot W_E^{(2)}\cdot \delta J
=\langle(\delta J\cdot\phi)^2\rangle_J
-\langle\delta J\cdot\phi\rangle_J^2\ge0 .
$$
等号が成り立つのは $\delta J\cdot\phi$ が測度の台でほとんど至る所一定のときである。例えば有限格子で密度が全ての実座標に正なら、その方向は $\delta J=0$ に限られる。従って正定値の領域では連鎖律で Hessian の逆を取れる。制約された場やゲージ固定前の冗長な場では、逆を使う前にどの部分空間で応答を定めるかを指定する。上限による Legendre–Fenchel の形は、$J\Phi-W_E(J)$ が $\Phi$ のアフィン関数であり、その上限が凸であることから、微分が存在しない場合にも凸性を保つ。

### 四点関数で既約頂点を取り出す

$G_{ab}=W^{(2)}_{ab}$ と略す。三点の式を $J_d$ で微分し、$G$ の各微分に三点の式を再代入する。四本の外脚の $G$ を共通因子として抜くと、
$$
W^{(4)}_{abcd}
=G_{ai}G_{bj}G_{ck}G_{dl}
\left[-\Gamma^{(4)}_{ijkl}
+\Gamma^{(3)}_{ijm}G_{mn}\Gamma^{(3)}_{nkl}
+\Gamma^{(3)}_{ikm}G_{mn}\Gamma^{(3)}_{njl}
+\Gamma^{(3)}_{ilm}G_{mn}\Gamma^{(3)}_{njk}\right] .
$$
最初の項は三点頂点そのものの微分、残る三項は三本の外脚のそれぞれを微分したものである。外脚を除いただけでは一粒子既約な量にならず、完全な二点線で結んだ三点頂点の積も引く必要がある。例えば反転対称な背景 $\Phi=0$ では $\Gamma^{(3)}=0$ のため、四点で初めてこの積が消え、四点を切断したものに負号を付ければ $\Gamma^{(4)}$ となる。

### Minkowski 規約での符号

実時間で $Z[J]=\int e^{i(S+J\phi)}$ とすれば $W=-i\log Z$、$\Phi=\delta W/\delta J$、$\Gamma=W-J\Phi$ が自然な組である。$\delta\Gamma/\delta\Phi=-J$ なので $\Gamma^{(2)}W^{(2)}=-1$。例えば自由場の $W=-JK^{-1}J/2$ では $\Phi=-K^{-1}J$、逆変換で $\Gamma=\Phi K\Phi/2$ となる。Euclid の定義と符号が異なるのは、実時間の $i$ を外した運動核で作用そのものを木の極限として回復するためである。$\log Z$ のみを共通記号として符号をそのまま移すことはできない。

### 体積当たりの Legendre 変換を曲率で検査する

一様源に対応する $w(J)=W_E(J)/V$ と $v(\varphi)=J\varphi-w(J)$ を取る。$\varphi=\langle M\rangle/V$、$w''(J)=\operatorname{Var}(M)/V$。従って局所的に逆を取れる領域では $v''(\varphi)=V/\operatorname{Var}(M)$。翻訳不変な場では運動量ゼロの完全二点に対応する曲率である。例えば $G_c(p)=[p^2+m^2+\Sigma(p^2)]^{-1}$ なら、$v''=m^2+\Sigma(0)$。粒子質量を定める極の条件は、解析接続した $p^2=-m_{\rm pole}^2$ における逆核の零であり、零運動量の曲率だけでは決まらない。

:::example[運動項の規格化が定数でも極質量と曲率は違う]
小運動量で $\Gamma_E^{(2)}(p)=A+Bp^2$、$A,B>0$ とし、この一次式を極の近くまで使えると仮定する。$v''=A$、極質量二乗は $A/B$、留数は $1/B$。例えば $A=4M^2,B=2$ なら曲率は $4M^2$、極は $2M^2$。場を $B^{1/2}$ 倍に取り直した有効ポテンシャルの曲率は $A/B$ になり、一次の微分近似の範囲で極と一致する。$p^4$ の項が大きければ、運動項だけを規格化しても一般には一致しない。
:::

### 一つの二点線を接続した六点の例

一点の反転対称な背景で $\Gamma'''=0$ とする。$\Gamma'(\Phi)=K\Phi+\lambda\Phi^3/6+g_6\Phi^5/120+\cdots$。逆を $\Phi=J/K+aJ^3+bJ^5+\cdots$ と置けば、
$$
a=-\frac{\lambda}{6K^4},\qquad
b=\frac{\lambda^2}{12K^7}-\frac{g_6}{120K^6}.
$$
第一は $J^3$ を消す条件。第二では $\Phi^3$ の $J^5$ の係数が $3a/K^2$ なので、$K b+(\lambda/6)(3a/K^2)+g_6/(120K^5)=0$ を解いた。従って $W^{(6)}=120 b=10\lambda^2/K^7-g_6/K^6$。係数 $10$ は六外脚を三本ずつの二つへ分ける $\binom63/2$ 通りであり、四点の既約頂点二つを一本の完全二点で結ぶ項を表す。例えば $\Gamma$ に独立した六次項を置かなくても、連結六点は四次項だけから非零になり得る。このように完全二点の外脚を除くだけで得る六点にも、低次の既約頂点の木が残る。

### Legendre 変換で消える定数と変わる源

$x$ の測度を $S(x)$ から $S(x)-a x+c$ に変えると、$W_a(J)=W(J+a)-c$。Legendre 変換の上限で $K=J+a$ と置けば $\Gamma_a(\Phi)=\Gamma(\Phi)-a\Phi+c$。従って線形の改変は $\Gamma$ の傾きだけを平行移動し、一定の $\Phi$ での二点以降の微分を変えない。例えば対称な作用で $a>0$ を入れた真空は $\Gamma'(\Phi)=a$ の点に移るので、真空に評価した二点は元の真空の値と異なり得る。「固定した背景で導関数が変わらない」ことと「平衡点で測る値が変わらない」ことは別である。

有限次元で変数を正則な線形変換 $\phi'=L\phi$ に取り直すと、積 $J\phi$ を保つ源は $J'=L^{-\mathsf T}J$。同じ確率測度の押し出しを使うと $W'(J')=W(L^\mathsf T J')$、$\Gamma'[\Phi']=\Gamma[L^{-1}\Phi']$。Hessian には両側の $L^{-1}$ が付き、連結二点は逆向きに $L$ を持つ。例えば場を二倍にすると分散は四倍、Hessian は四分の一。この変換は線形なので $\Phi'=\langle L\phi\rangle=L\Phi$。非線形な場の取り替えでは平均と関数が一般には交換しないため、同じ合成則を無条件に使えない。

## ループ展開と 1 ループ有効作用

### 背景を固定して揺らぎを積分する

ループの回数を数えるため一時的に $\hbar$ を戻し、$W_E=\hbar\log Z_E$、指数を $(-S_E+J\phi)/\hbar$ とする。停留場 $\phi_c$ は $S_E^{(1)}[\phi_c]=J$ を満たす。$\phi=\phi_c+\sqrt{\hbar}\eta$ と置くと一次は消え、二次は $-\eta S_E^{(2)}[\phi_c]\eta/2$。有限正則化の下で鞍点が支配的で Hessian が正なら、二次揺らぎを積分して $\det S_E^{(2)}[\phi_c]^{-1/2}$ が掛かる。Legendre 変換で $\phi_c=\Phi+\mathcal O(\hbar)$ のずれの一次項は $J-S_E^{(1)}=0$ で消える。従って
$$
\Gamma_E[\Phi]=S_E[\Phi]+\tfrac{\hbar}{2}\Tr\log S_E^{(2)}[\Phi]+\mathcal O(\hbar^2)
\label{eq-qft-one-loop-effective}
$$
である。背景に依存しない行列式は引ける。複素 Grassmann 場の二次積分は行列式を分子に与えるので、$\Gamma_E$ への寄与は $-\hbar\Tr\log K$ となる。相互作用頂点の脚数を $n$ とすると $\hbar^{n/2-1}$ が付き、連結な真空グラフの $W_E$ で $\hbar^{I-V+1}$ となる。この数がループ数である。

### 二点の行列式とループの多重度

有限体積で二つのモードだけを保ち、背景依存の Hessian を $K=\begin{pmatrix}A&b\\b&C\end{pmatrix}$、$A,C>0$、$b^2<AC$ とする。基準 $b=0$ からの一ループの差は $\hbar\log(1-b^2/(AC))/2$。$\Tr\log(K_0+V)=\Tr\log K_0+\Tr\log(1+K_0^{-1}V)$ を小さい $b$ で展開すると、一次項は非対角なので零、二次項のトレースは $2b^2/(AC)$。従って係数は $-\hbar b^2/(2AC)$ となる。例えば $A=2,C=8,b=2$ なら厳密な差は $(\hbar/2)\log(3/4)$、二次の近似は $-\hbar/8$。二種のモードが交互に一周するので、奇数回の挿入は全て消える。

### 二ループで Legendre 変換が橋の線を除く

有限次元の一変数で $S''(\phi_c)=K>0$、$S'''(\phi_c)=B$、$S''''(\phi_c)=C$ とし、$S'(\phi_c)=J$。鞍点の揺らぎを $\sqrt{\hbar}\eta$、$D=1/K$ とする。一ループの $W$ は $w_1=-\frac12\log K$。二ループの $W$ には四次項一つと三次項二つがあり、Gaussian の $\langle\eta^4\rangle=3D^2$、$\langle\eta^6\rangle=15D^3$ を使うと
$$
w_2=-CD^2/8+5B^2D^3/24 .
$$
二つの三次頂点が三本の線でつながる縮約は $3!=6$、各頂点で自己縮約して残り一本でつながる縮約は $3^2=9$ である。前者だけが 1PI である。

$\Phi=\phi_c+\hbar\partial_J w_1+\cdots$。$\partial_J\phi_c=1/K=D$ だから $\partial_J w_1=-BD^2/2$。$J$ を $\Phi$ の関数として逆にすると $J=S'(\Phi)-\hbar K\partial_Jw_1+\cdots$ となる。$J\Phi-W_0(J)$ の二次展開は停留点で一次が消え、$-D(\delta J)^2/2$。$-\hbar w_1(J)$ のずれは $-\hbar(\partial_Jw_1)\delta J$ である。従って二ループでの和は
$$
\Gamma_2=-w_2+\frac{K}{2}(\partial_Jw_1)^2
=\frac18 CD^2-\frac1{12}B^2D^3 .
$$
ここでは全て $\Phi$ で評価する。この相殺で $9/(2(3!)^2)=1/8$ の橋の縮約が消え、$6/(2(3!)^2)=1/12$ の三線の方だけが残る。

:::example[対称な背景では三次頂点が消える]
$S(\phi)=k\phi^2/2+\lambda\phi^4/24$ を $\Phi=0$ で評価すると $K=k,B=0,C=\lambda$。従って $\Gamma_2=\lambda/(8k^2)$。$k=2,\lambda=1$ では $1/32$ である。一方、$\Phi\ne0$ なら $B=\lambda\Phi$ なので三線の縮約も現れる。真空でのループの有無は微分した作用をどの背景で評価するかに依存する。
:::

### 非一様背景で行列式を頂点に展開する

$\phi^4$ 理論で $S_E^{(2)}[\Phi]=K_0+U$、$K_0=-\partial^2+m^2$、$U(x)=\lambda\Phi(x)^2/2$ とする。有限正則化で正の Hessian を仮定し、$\Tr\log$ の差を $\frac12\Tr\log(1+K_0^{-1}U)$ と書く。非可換な $K_0,U$ でも成立するのは、$\Tr\log=\log\det$ と $\det(K_0+U)=\det K_0\det(1+K_0^{-1}U)$ を使うからである。循環積を展開すると
$$
\frac12\Tr(G_0U)-\frac14\Tr(G_0UG_0U)+\frac16\Tr(G_0UG_0UG_0U)+\cdots .
$$
これは $G_0^{1/2}UG_0^{1/2}$ の固有値の絶対値が一より小さいとき通常の級数として収束し、摂動論では形式級数として読む。第一項を背景で二回微分すると $\lambda G_0(x,x)/2$ が二点へ加わる。第二項は $-\lambda^2\int \Phi(x)^2G_0(x,y)^2\Phi(y)^2/16$。一つの外点分割について $x,y$ の選び方が二通り、各点の二つの微分に $2!$ ずつ掛かるため、四点の一チャネルは $-\lambda^2 I(P)/2$。これで行列式が前章の泡の係数を再現する。六回微分する場合は第三項が一ループの三角形を作る。

:::example[定数背景で六次の有限係数を読む]
$m>0$、$|\lambda\varphi^2/(2m^2)|\ll1$ とする。上の第三項は体積当たり $U^3 I_3(m^2)/6$。$I_3(m^2)=1/(32\pi^2m^2)$ なので $\lambda^3\varphi^6/(1536\pi^2m^2)$。これは四次元の局所六次頂点の質量次元 $-2$ を持ち、四点までの新しい紫外対項を増やさない一ループの有限効果である。
:::

### 有効作用を背景付きの積分方程式に戻す

Euclid の定義を $\hbar$ 付きで使う。$Z_E[J]=e^{(J\Phi-\Gamma_E[\Phi])/\hbar}$ に $\phi=\Phi+\eta$ を入れ、$J=\Gamma_E^{(1)}[\Phi]$ を用いると
$$
e^{-\Gamma_E[\Phi]/\hbar}
=\int\mathcal D\eta\,
  \exp\left[-\frac{S_E[\Phi+\eta]
  -\int\Gamma_E^{(1)}[\Phi]\eta}{\hbar}\right]
$$
を得る。これは元の積分と同じ測度での厳密な式であり、右辺の $\eta$ の平均は $\langle\phi\rangle_J-\Phi=0$。従って量子の背景展開では、古典の $S^{(1)}$ を一回だけ引くのではなく、一定の平均場を保つ線形の打消しを各次数で調整する。木では $\Gamma_E^{(1)}=S_E^{(1)}$ なので通常の一次の相殺になるが、一ループでは $\tfrac12\Tr(G\,\delta S^{(2)})$ が線形の打消しに加わる。

有限次元の添字を使って $G_{ij}=(S_E^{(2)})^{-1}_{ij}$、$B_{ijk}=S_E^{(3)}$、$C_{ijkl}=S_E^{(4)}$ とする。二ループの既約な項は
$$
\Gamma_{E,2}=\frac18 C_{ijkl}G_{ij}G_{kl}
-\frac1{12}B_{ijk} B_{abc}G_{ia}G_{jb}G_{kc}.
$$
四脚一つに三種の対を作る $3/4!=1/8$ と、三脚二つで全線が頂点をまたぐ $3!/[2(3!)^2]=1/12$ が係数を決める。符号は指数の $-S_E$ と作用への負の対数から出る。一線だけを頂点間に結ぶパターンは $9/[2(3!)^2]=1/8$。各端で二線を自己縮約した形なので、一次の打消しの係数 $\Gamma_{E,1}^{(1)}{}_i=B_{ijk}G_{jk}/2$ を用いた線形項がちょうど同じ分を消す。これは一変数で見た相殺を一般の添字に書いたものである。場では添字の和を座標積分に読み替え、正則化した $G(x,y)$ を使う。

:::example[四次スカラーの非一様な二ループ作用]
$S_E^{(3)}$ の局所係数は $\lambda\Phi(x)$、$S_E^{(4)}$ は $\lambda$。従って二ループの先頭は
$$
\frac{\lambda}{8}\int\dd^4x\,G(x,x)^2
-\frac{\lambda^2}{12}\int\dd^4x\,\dd^4y\,\Phi(x)G(x,y)^3\Phi(y).
$$
真空 $\Phi=0$ に評価するだけなら第二項は零である。背景で二回微分してから零にすると、第二項はサンセットの $-\lambda^2/6$ を与える。第一項は $G$ 自身も背景に依存するので、背景微分には $G$ の微分が入る。固定背景の値を読む操作と、既約頂点を求める微分は順序が異なる。
:::

### 重い場と軽い場の揺らぎを同時に対角化する

一定背景で二自由度の質量二乗の Hessian を $H=\begin{pmatrix}A&B\\B&C\end{pmatrix}$ とし、$A>0,AC>B^2$ とする。各運動量で核は $p^2 1+H$。共通の直交変換で $H$ を対角化すれば、固有値は
$$
M_\pm^2=\frac{A+C\pm\sqrt{(A-C)^2+4B^2}}2 .
$$
従って一ループのポテンシャルは $M_+^2,M_-^2$ の対数の和である。別の計算として、重い方とみなす第二の揺らぎを先に積分する。平方完成で $\det(p^2+C)^{-1/2}$ と、第一の揺らぎの逆核 $p^2+A-B^2/(p^2+C)$ が残る。積分を続けると二つの行列式の積が、もとの $\det(p^2 1+H)^{-1/2}$。従って二つの積分順は同じ結果を与える。前者は完全な有効作用、後者で第二の揺らぎだけを消す段階は第一の揺らぎをまだ積分していない低エネルギー記述に対応する。

:::example[運動量を保った核と零運動量近似]
$C\gg A$、$|p^2|\ll C$ とする。第一の逆核は $A-B^2/C+(1+B^2/C^2)p^2+\cdots$。従って二次の質量の照合には $-B^2/C$ だけでよいが、運動項の係数も同じ展開で $B^2/C^2$ だけ変わる。例えば $A=M^2/4,C=M^2,B=M^2/4$ では零運動量の逆核は $3M^2/16$、低運動量の $p^2$ の係数は $17/16$。二次の核だけを低い尺度へ展開するときにも、一定場の行列式を取り出すだけでは運動項が分からないことが確認できる。
:::

### Dirac 場の行列式を実スカラーと比べる

一定背景 $\varphi$ に Yukawa 質量 $m_F(\varphi)=m+g\varphi$ を持つ一つの Dirac 場を取る。Euclid 運動量 $p$ に対する核は $i\gamma_E\cdot p+m_F$。$p\ne0$ で $\gamma_E\cdot p$ の固有値は $\pm|p|$ が二重ずつである。二乗が $p^2$ を与え、任意の直交するガンマとの反交換により固有空間の次元が等しいからである。従って行列式は $(p^2+m_F^2)^2$。Grassmann 積分を作用へ戻すと $\Gamma_{1,F}=-2\int_p\log(p^2+m_F^2)$、一実スカラーの $+\tfrac12$ の $-4$ 倍である。従って規約をそろえれば、
$$
V_{1,F}
=-\frac{m_F^4(\varphi)}{16\pi^2}
\left[\log\frac{m_F^2(\varphi)}{\mu^2}-\frac32\right].
$$
例えば $m=0$ の Yukawa 項なら $\varphi^4\log\varphi^2$ の係数は $-g^4/(16\pi^2)$。スカラーや光子による係数にこれを足す場合、高い次数の近似を必要としない小さい対数の範囲で符号を比較する。フェルミオンの数を二重に数える必要はなく、粒子と反粒子を含む四つの成分は既に行列式に入っている。

:::example[一様な背景で補正が打ち消し合う条件]
質量二乗が同じ $M^2(\varphi)$ の実スカラーを四成分、Dirac 場を一種類持つ二次の模型を考える。同じ $\overline{\rm MS}$ の基準での一ループのポテンシャルは成分数 $4-4=0$ により消える。一方、二種の揺らぎの演算子が非一様な背景でも同じスペクトルを持つとは限らないため、微分を含む作用まで消えるとはこの例だけからは言えない。一様な比較では質量の固有値だけを使い、背景が変化する場合は導関数を含む核を調べる。
:::

## 有効ポテンシャルと Coleman–Weinberg 機構

一様な $\Phi=\varphi$ では Hessian が $-\partial_E^2+M^2(\varphi)$、$M^2=V''(\varphi)$ となる。$\frac12\int_p\log(p^2+M^2)$ を $M^2$ で微分すると :eqref[eq-qft-loop-master] の $I_1/2$。極を引いて積分し直すと
$$
V_1(\varphi)=\frac{M^4(\varphi)}{64\pi^2}
\left[\log\frac{M^2(\varphi)}{\mu^2}-\frac32\right]
\label{eq-qft-scalar-one-loop-potential}
$$
となる。$M^2<0$ では一様背景の Gaussian 自体が安定でなく、同じ実積分でこの式を使えない。厳密な Legendre 有効ポテンシャルは $W^{(2)}$ の共分散の正性と Legendre 変換から凸である。摂動的な一様背景の式にはこの凸性を点ごとに要求できない。

質量ゼロのスカラー QED で $\varphi=\sqrt{2}|\phi|$ とすると光子質量は $e\varphi$ である。Landau ゲージで Euclid 運動量に横なベクトル揺らぎの行列式が $\varphi^4\log\varphi^2$ に $B=3e^4/(64\pi^2)$ を与える。$\lambda$ を $e^4/(16\pi^2)$ と同程度に取り、スカラーの $\lambda^2$ のループを高次として落とすと、有限な四次係数とくりこみ尺度を組み合わせて
$$
V(\varphi)=B\varphi^4\left(\log\frac{\varphi^2}{v^2}-\tfrac12\right)+C
$$
と書ける。定数 $-1/2$ は $\partial_\varphi V(v)=0$ で決まる。$V'(\varphi)=4B\varphi^3\log(\varphi^2/v^2)$、$V''(v)=8Bv^2>0$。無次元の古典結合から量子補正で尺度が現れるこの例が :term[Coleman–Weinberg 機構]{en="Coleman–Weinberg mechanism" yomi="こーるまんわいんばーぐきこう"} である。大きな対数の領域は固定次数の式で扱えず、背景の位置自体は一般にゲージ依存である。

### 一様場の二つの固有値を別々に積分する

二成分の $V=\lambda(\varphi^2-v^2)^2/4$ を取る。一様な半径 $R$ で、半径と直交する単位ベクトルに Hessian を作用させると、それぞれ $M_\parallel^2=\lambda(3R^2-v^2)$、$M_\perp^2=\lambda(R^2-v^2)$ になる。$R>v$ なら両者が正なので、一ループの有効ポテンシャルは :eqref[eq-qft-scalar-one-loop-potential] を二質量について足す。例えば $R^2=2v^2$ では質量二乗の比は $5:1$ で、対数の前に掛かる質量四乗の比は $25:1$ である。同じ相互作用でも、揺らぎの固有値を一つとみなすとこの重みを失う。$R\downarrow v$ では $M_\perp^4\log M_\perp^2\to0$ だが、その二階以上の背景微分では赤外極限を個別に調べる必要がある。

### 凸なポテンシャルと二相の混合

正規化体積 $\mathcal V$ の模型で相が $\varphi=\pm v$ の二つだけを取り、密度源を $J$ とする。$w_{\mathcal V}(J)=\mathcal V^{-1}\log[2\cosh(\mathcal V Jv)]$ だから $\partial_Jw_{\mathcal V}=v\tanh(\mathcal V Jv)$。有限体積では零源で平均は零である。一方、先に $\mathcal V\to\infty$ とすると $w_\infty=v|J|$。Legendre–Fenchel 変換 $\sup_J[J\Phi-v|J|]$ は $|\Phi|\le v$ で零、範囲外で発散する。この模型で中間の $\Phi$ は二相の重みを変えて得られる。実際 $\Phi=v(2p-1)$ は $+v$ の相の確率を $p$ とした平均である。場の理論でも相分離の界面費用が体積あたりで消える極限では、このような平坦部が現れる。

### ベクトル行列式の有限定数

スカラー QED で $\varphi$ を一定に取り、$m_A^2=e^2\varphi^2$ とする。二次作用のゲージ・位相場の横縦分解では、横射影に $p^2+m_A^2$、縦に $p^2/\xi+m_A^2$ が現れる。混合項は $\partial\cdot A$ に比例し、Landau ゲージ $\xi\downarrow0$ で横成分に寄与しない。有限の非零 $p$ で縦の規格化を背景なしで割った後の比は $1+\xi m_A^2/p^2$。従って Landau ゲージで背景依存のベクトル項は $d-1$ 成分の和である。位相と半径のスカラー Hessian による質量二乗は $\lambda\varphi^2$ に比例するので、$\lambda$ を $e^4/(16\pi^2)$ の大きさに取る近似ではその行列式はここでの先頭に含めない。

一実スカラーの未減算の一ループの項は $M^4[-1/\bar\epsilon+\log(M^2/\mu^2)-3/2]/(64\pi^2)$。$d-1=3-2\epsilon$ を掛けると、極との積から有限な $+2$ が残る。従って $\overline{\rm MS}$ のベクトル項は
$$
V_{1,A}=\frac{3m_A^4}{64\pi^2}
\left(\log\frac{m_A^2}{\mu^2}-\frac56\right).
$$
$3$ は四次元での横部分空間の次元、$5/6=3/2-2/3$ は次元を留めたまま極を引いた結果である。古典ポテンシャルを $\lambda\varphi^4/24$ とすると、停留条件は
$$
\frac{\lambda}{6}+4B\left[\log(e^2v^2/\mu^2)-\frac13\right]=0,\qquad
B=\frac{3e^4}{64\pi^2}.
$$
例えば $\mu=ev$ を選ぶと $\lambda=8B$。この条件を消去して、前の $B\varphi^4(\log(\varphi^2/v^2)-1/2)$ の形に戻る。ここでは停留点が摂動的に制御できるように $\mu$ を選んだ。有限の係数を別の方式で定義すれば停留条件中の定数も変化する。

### 正の Hessian を持つ一様背景のくりこみ

$m^2>0$、$V_0=m^2\varphi^2/2+\lambda\varphi^4/4!$ とし、$M^2=m^2+\lambda\varphi^2/2>0$ の範囲を考える。未減算の $V_1$ の極は $-M^4/(64\pi^2\bar\epsilon)$。よって定数、$\varphi^2$、$\varphi^4$ の極はそれぞれ $-m^4/(64\pi^2\bar\epsilon)$、$-\lambda m^2\varphi^2/(64\pi^2\bar\epsilon)$、$-\lambda^2\varphi^4/(256\pi^2\bar\epsilon)$。四次の対項を $\delta\lambda\,\varphi^4/24$ とすると、必要な $\delta\lambda=3\lambda^2/(32\pi^2\bar\epsilon)$ で、前章の三チャネルの泡と一致する。二点でも $\delta m^2/2$ の係数が一致する。定数は真空エネルギーを固定する場合にのみ独立の条件が要る。背景一様の行列式からは、運動項のくりこみ条件を決められない。

例えば小さい $\lambda\varphi^2/m^2$ で :eqref[eq-qft-scalar-one-loop-potential] を展開する。$f(M^2)=(M^2)^2(\log(M^2/\mu^2)-3/2)$ と置けば、$f'(m^2)=2m^2[\log(m^2/\mu^2)-1]$、$f''(m^2)=2\log(m^2/\mu^2)$。従って零背景の二点の補正は $\lambda m^2[\log(m^2/\mu^2)-1]/(32\pi^2)$、四点の補正は $3\lambda^2\log(m^2/\mu^2)/(32\pi^2)$。$V''$ と $V''''$ はいずれも零外部運動量の 1PI 頂点であり、有限運動量の減算方式では条件をそのまま使わない。

### Abel ゲージで位相揺らぎと縦揺らぎを追う

実背景 $\varphi>0$ の周りで $\phi=(\varphi+h+i\chi)/\sqrt2$ とする。$\chi$ は位相方向の実揺らぎである。Euclid の $D=\partial+iea$ では二次混合が $m_A a_\mu\partial_\mu\chi$、$m_A=e\varphi$。ゲージ条件を $F=\partial\cdot a+\xi m_A\chi$ とし、$F^2/(2\xi)$ を加えると交差項 $m_A\chi\partial\cdot a$ が積分で混合を消す。局所変換 $\delta a_\mu=\partial_\mu\omega$、$\delta\chi=-e(\varphi+h)\omega$ で $\delta F=(\partial^2-\xi m_A^2-\xi e m_Ah)\omega$。従って背景依存の正規化行列式があり、二次では $p^2+\xi m_A^2$ である。線形な Maxwell のゲージ条件では行列式が場に依らなかったが、この条件は物質の位相を含むので結果が変わる。

古典の位相の Hessian を $m_G^2=V_0'(\varphi)/\varphi$ とすると、縦偏極は半個分の $\log(p^2+\xi m_A^2)$、位相は半個分の $\log(p^2+m_G^2+\xi m_A^2)$、行列式は一個分の $-\log(p^2+\xi m_A^2)$ を与える。縦の $1/\xi$ は背景と独立なため規格化で引いた。半径揺らぎと横偏極の二次核は $\xi$ に依存しない。従って正の二次核を持つ範囲で同じ正則化を用いると
$$
\partial_\xi V_1=
-\frac{m_A^2}{2}\, m_G^2
\int_p\frac1{(p^2+\xi m_A^2)(p^2+m_G^2+\xi m_A^2)} .
$$
右辺に $m_G^2=V_0'/\varphi$ を入れれば、一ループのゲージ依存が古典の方程式に比例する形になる。通常のループ展開で $V_0'=0$ の点を使うと、一ループの極値の値はこの変化から守られる。他方、停留点の位置は $V'$ の条件で決まり、$\xi$ を変えると動き得る。Coleman–Weinberg のように木の結合もループ次数に合わせて数える場合は、同じ次数に寄与する項をその数え方でまとめてから極値を比較する。

### 負の曲率での解析接続を安定な積分と区別する

一様な背景の Hessian が $p^2+M^2$ で、$M^2=-\kappa^2<0$ とする。このとき $|p|<\kappa$ のモードは負の二次係数を持つため、実の揺らぎの Gaussian は発散する。対数を $M^2=-\kappa^2+i0$ へ形式的に続ければ $\log(p^2-\kappa^2+i0)$ の負の区間で虚部 $\pi$。従って体積当たりの一ループ項の虚部は
$$
\frac{\pi}{2}\int_{|p|<\kappa}\frac{\dd^4p}{(2\pi)^4}
=\frac{\kappa^4}{64\pi}.
$$
四次元球の体積 $\pi^2\kappa^4/2$ を使った。これは :eqref[eq-qft-scalar-one-loop-potential] の対数の形式的な虚部にも一致する。逆側の $-i0$ に続ければ符号は逆である。実の正測度に対する厳密な $\Gamma_E$ は凸な実汎関数なので、この継続を厳密な Legendre 変換と解釈しない。一様な鞍点に負モードがあるという事実を、近似の計算から読み取る。一方、準安定相の崩壊率を求めるには、どの背景の境界条件で積分するかと、非一様な鞍点を別に扱う必要がある。

### 多相の凸包を上限から作る

規格化した体積当たりの $w(J)$ が二つの安定な一様相から $w(J)=\max(J\varphi_- -f_-,J\varphi_+-f_+)$ と近似される場合を考える。ここで $\varphi_-<\varphi_+$、$f_\pm$ は二相の自由エネルギー密度である。二つの直線が交わる源は $J_c=(f_+-f_-)/(\varphi_+-\varphi_-)$。$\varphi_-\le\Phi\le\varphi_+$ について $J\Phi-w(J)$ の傾きは $J<J_c$ で $\Phi-\varphi_-\ge0$、$J>J_c$ で $\Phi-\varphi_+\le0$。従って上限は $J_c$ で取られ、
$$
v(\Phi)=f_-+J_c(\Phi-\varphi_-)
$$
となる。二相の自由エネルギーを結ぶ直線である。例えば $f_+=f_-$ なら区間は平坦。等しくない場合も、共存する源は零とは限らない。同じ傾きの直線を作る二相の値は、一つの背景で Hessian を計算するだけでは決まらない。有限系で二相からの指数を和にしてから対数を取ると折れ曲がりは丸まり、体積を先に無限大に取る操作が区別を作る。

### 有限温度では零周波数のモードが別の冪を作る

温度 $T$ の一実スカラー揺らぎで、一定背景の Hessian 質量を $M\ge0$ とする。真空部分を引いた一ループの熱ポテンシャルは、周期の時間核の行列式から
$$
V_{1,T}(M)=T\int\frac{\dd^3 k}{(2\pi)^3}\log(1-e^{-\sqrt{k^2+M^2}/T})
$$
である。$M^2$ を微分すると $\int\dd^3k/(2\pi)^3\, n_B(\omega_k)/(2\omega_k)$。$M=0$ では $k=T y$ とし、$\int_0^\infty y/(e^y-1)\dd y=\sum_{n\ge1}1/n^2=\pi^2/6$ から $T^2/24$。従って熱ポテンシャルの最初の背景依存は $M^2T^2/24$ である。例えば $\lambda\varphi^4/24$ の揺らぎで $M^2=\lambda\varphi^2/2$ とすれば、熱質量の二乗の先頭は $\lambda T^2/24$。この符号は正の結合なら一様な原点を安定にする向きである。

次の非解析な項はボースの零周波数から出る。差分の質量微分を $M>0$ と零の間で取り、零周波数だけについて
$$
\frac{T}{2}\int\frac{\dd^3 k}{(2\pi)^3}
  \left[\frac1{k^2+M^2}-\frac1{k^2}\right]
=-\frac{TM}{8\pi}
$$
を得る。球積分で $\int_0^\infty\dd k/(k^2+M^2)=\pi/(2M)$。$M^2$ に関して積分すると $-T M^3/(12\pi)$ である。零以外の Matsubara モードには $2\pi T$ 以上の間隔があり、小 $M$ での解析な展開からこの $M^3$ は出ない。真空での $M^4\log M^2$ の非解析性とも次数が異なる。フェルミオンの反周期モードは零を持たないため、同じ機構での $M^3$ は生じない。質量が零の点の近くでは、熱質量を無視した固定ループ展開が赤外で不安定になり得る。ここでは一ループ行列式のどの部分が非解析な冪を供給するかを計算した。

## 対称性と Ward–Takahashi 恒等式の汎関数的導出

線形な変換 $\delta\phi^a=\epsilon(T\phi)^a$ が作用と測度を保つとする。源つき積分で変数変換して $0=J_aT^a{}_b W_E^{(1)b}$ を得る。従って :eqref[eq-qft-legendre-effective] を用い
$$
\int\dd^4x\,\frac{\delta\Gamma_E}{\delta\Phi^a(x)}(T\Phi)^a(x)=0
\label{eq-qft-effective-ward}
$$
となる。局所的な $\epsilon(x)$ では作用の変化がカレントと $\partial_\mu\epsilon$ の積になるため、源の接触項を含む恒等式を得る。例えば $O(2)$ の一様な二成分なら $\varphi_1\partial_2 V-\varphi_2\partial_1 V=0$ なので、半径一定の円上で $V$ が一定である。変換が非線形なときは $\langle\delta\phi\rangle$ を $\delta\Phi$ へ置けるとは限らず、必要な合成演算子へ別の源を導入する。測度の Jacobian が非自明な場合はその寄与も残る。

### 局所変換から接触項を展開する

保存する内部電流だけを持つ Euclid の線形対称性に対し、$\delta S_E=\int j^\mu\partial_\mu\epsilon$ とする規約を採る。変数変換で測度が不変なら
$$
0=\int\dd^4x\,\epsilon(x)
 \left(\langle\partial_\mu j^\mu(x)\rangle_J+J_a(x)(T\Phi)^a(x)\right).
$$
境界で $\epsilon=0$ として作用変分を部分積分した結果である。これを $J_b(y)$ で微分して $J=0$ にすると
$$
\partial_\mu\langle j^\mu(x)\phi^b(y)\rangle_c
=-\delta^4(x-y)\langle(T\phi)^b(y)\rangle.
$$
さらにもう一回微分すれば、接触項は二つの外点のどちらにも現れる。例えば $O(2)$ で $T\phi=(-\phi_2,\phi_1)$、$\langle\phi_1\rangle=v$ の状態を選ぶと、$\phi_2$ の挿入には $-v\delta^4(x-y)$ が残る。例えば純相を無限体積極限で選択してから源を外した場合、対称性で動く場が接触項を持つことが分かる。この型の相関を [対称性の自発的破れ](/physics/quantum-field-theory/symmetry-breaking/) で使う。

### 停留背景上の Hessian の零方向

:eqref[eq-qft-effective-ward] をもう一度 $\Phi$ で微分すると
$$
\int\dd^4x\,\Gamma^{(2)}_{ba}(y,x)(T\Phi)^a(x)
+\Gamma^{(1)}_a(y)T^a{}_b=0
$$
となる。停留背景で第二項は零なので、$T\Phi$ が非零なら Hessian の零方向である。一様背景ではこれは運動量ゼロの逆伝播関数に対する条件となる。例えば $V=\frac{\lambda}{4}(\Phi^2-v^2)^2$ の背景 $(v,0)$ では $T\Phi=(0,v)$ を Hessian $\diag(2\lambda v^2,0)$ が消す。有限体積の凸な有効作用と、一つの純相を源で選択してから用いる相関関数では極限の順序を指定する。

### 外部接続を源として局所恒等式を保つ

内部対称性に対し、$B_\mu$ を外部の接続として、物質作用を $\phi$ と $B_\mu$ の同時局所変換で不変にする。例えば実二成分を直交行列で回し、$D_\mu=\partial_\mu+B_\mu T$ とする。$\delta\phi=\epsilon T\phi$ に対して $\delta B_\mu=-\partial_\mu\epsilon$ とすれば、微分項の変化は相殺する。対称な正則化の下で物質場を変数変換し、物質源についてのみ Legendre 変換すると、
$$
\partial_\mu\frac{\delta\Gamma_E[\Phi,B]}{\delta B_\mu(x)}
+\frac{\delta\Gamma_E[\Phi,B]}{\delta\Phi^a(x)}(T\Phi)^a(x)=0 .
$$
第一項は $\delta B_\mu$ の積分を部分積分したもの。$\Phi$ の微分を取ると、同じ式の第二項は外脚の変換と $\delta(x-y)$ を作る。従って外部接続を導入した場合、保存電流の挿入と変換の接触項が一つの作用の変分で得られる。$B=0$ に置いてから局所変換を行うと、第一項を作る変分を失う。

:::example[一様な対称背景の制約]
$O(2)$ の不変な有効ポテンシャルを $V=U(r)$、$r=\Phi_1^2+\Phi_2^2$ と書く。$\partial_a V=2\Phi_a U'$、$V_{ab}=2\delta_{ab}U'+4\Phi_a\Phi_b U''$。半径 $v>0$ の停留点では $U'(v^2)=0$ なので横の固有値は零、半径方向は $4v^2 U''(v^2)$。半径方向の正値性には $U''(v^2)>0$ が要る。Ward 恒等式から直接決まるのは横方向だけである。
:::

### 二つの電流微分から接触項を取り出す

実二成分の $O(2)$ 模型で $D_\mu\phi=\partial_\mu\phi+B_\mu T\phi$。物質作用の外部接続による一回微分は $j_\mu=(D_\mu\phi)^\mathsf T T\phi$。同じ場を固定してさらに微分すると、$(T\phi)^2\delta_{\mu\nu}\delta(x-y)$。外部接続のみによる正規化した $W_E[B]=\log\int e^{-S_E[\phi,B]}$ を二回微分すれば、
$$
\frac{\delta^2W_E}{\delta B_\mu(x)\delta B_\nu(y)}
=\langle j_\mu(x)j_\nu(y)\rangle_c
-\delta_{\mu\nu}\delta(x-y)\langle (T\phi(x))^2\rangle .
$$
従って保存電流の応答の横性は、電流二点の連結相関だけでなく作用の二階微分も含めて成り立つ。実際、$W_E[B-\partial\epsilon]=W_E[B]$ を $\epsilon$ と $B$ について微分すると、上の応答の発散は零となる。この接触項はスカラー QED の二光子頂点と同じく、共変微分を二乗したことで現れる。電流を $B=0$ での式に固定してから微分すると、作用の二階項を落とし、保存の条件が成立しなくなる。

### 非線形な変換では平均値だけで恒等式が閉じない

有限正則化で場の変数変換 $\delta\phi^a=\epsilon R^a(\phi)$ を考え、作用と測度の積が不変とする。源付きでは $0=J_a\langle R^a(\phi)\rangle_J$。Legendre 変換後は $\Gamma^{(1)}_a\langle R^a(\phi)\rangle_{J[\Phi]}=0$ であり、一般には $\langle R^a(\phi)\rangle_J=R^a(\Phi)$ ではない。例えば二次の変換 $R^a(\phi)=R^a_{bc}\phi^b\phi^c$ なら、平均は $R^a_{bc}[\Phi^b\Phi^c+W^{(2)}_{bc}]$。高次の積にはさらに高次の連結量が現れる。従って相関が有限でない局所積を含む場合、変換に対応する合成演算子の源と正則化を指定する必要がある。

別の例として、$X$ が分散 $D$ の Gaussian で、源を $J X^2$ に付けると、$\partial_J W$ は $X^2$ の平均であって $\langle X\rangle_J^2$ ではない。反転対称性で後者は零、前者は $D/(1-2DJ)$。この差は量子の揺らぎの最も簡単な例である。場の変換が線形なら共分散を追加せずに $\Gamma^{(1)} T\Phi=0$ が閉じる。非線形な場合に同じ簡潔な形を使うには、拡張した生成汎関数で合成源の微分を用いる。

## 背景場の方法

### 背景を変換しながら量子場を積分する

$\phi=\bar\phi+\eta$ と分ける。線形対称性なら $\bar\phi,\eta$ を一緒に変換し、$S^{(2)}[\bar\phi]$ を対応する共役変換へ移せる。有限次元で直交変換なら Hessian の行列式は同じであり、対称な正則化では :eqref[eq-qft-one-loop-effective] も不変になる。ゲージ理論では背景接続 $\bar A$ を接続として、量子揺らぎ $a=A-\bar A$ を随伴テンソルとして変換する。背景共変な条件 $\bar D_\mu a^\mu=0$ を使えばゲージ固定も背景変換を保つ。この区別により二次揺らぎの行列式から背景の $F_{\mu\nu}^2$ の係数を読める。[非可換ゲージ理論の量子化](/physics/quantum-field-theory/gauge-quantization/) で計算する。

例えば二成分の $V=\frac{\lambda}{4}(\Phi^2-v^2)^2$ では一様背景の Hessian は $\lambda(\Phi^2-v^2)1+2\lambda\Phi\Phi^{\mathsf T}$ である。半径方向と直交方向は固有ベクトルであり、背景を回すと固有値を保つ。行列式で計算した一ループ補正が角度によらないことを直接確認できる。

### 背景 Hessian から運動項の補正を求める

実揺らぎの二次演算子を $K=-\partial^2+M^2+U(x)$、$M>0$ とする。$\Tr\log$ の二次項は $-\Tr(G_0UG_0U)/4$ である。Fourier 変換すると $-\frac14\int_p U(p)U(-p)I(p)$、$I(p)=\int_k[(k^2+M^2)((k+p)^2+M^2)]^{-1}$。小さい $p$ で $I(p)-I(0)=-p^2/(96\pi^2M^2)+\cdots$ だから、二微分の有限な寄与は $\int(\partial U)^2/(384\pi^2M^2)$ である。例えば $U=h\Phi$ なら、標準の係数 $1/2$ をくくった波動関数の補正は $h^2/(192\pi^2M^2)$。背景を一様に置いた計算はこの項を捉えないため、有効ポテンシャルからだけでは波動関数の係数を求められない。

### 背景場を変えて同じ全場を記述する

普通の線形分解 $\phi=B+\eta$ では、背景 $B$ を固定して量子場 $\eta$ に源 $J$ を加えると
$$
Z_B[J]=\int\mathcal D\eta\,e^{-S[B+\eta]+J\eta}
=e^{-J B}Z[J]
$$
となる。並進不変な規格化測度で変数を $\phi$ に移した。従って $\langle\eta\rangle_{B,J}=\Phi-B$。量子場の平均で Legendre 変換すれば、$J(\Phi-B)-[W(J)-JB]=\Gamma_E[\Phi]$ となる。背景に対する恒等式 $\delta_B W_B=-J$ も同じ式に含まれる。この結論の仮定は、分解の外に背景 $B$ を陽に含む項がないことである。背景ゲージ固定を行う場合は $\bar D_\mu a^\mu$ に背景が陽に入るため、単純な変数移動だけで $Z_B=e^{-JB}Z$ になるとは限らない。

### 緩やかな背景で高次の微分を見積もる

実揺らぎの二次項で、$x(1-x) p^2/M^2$ に対する対数級数をもう一つ保つ。$\int_0^1 x^2(1-x)^2\dd x=1/30$ なので
$$
I(p)-I(0)=-\frac{p^2}{96\pi^2M^2}
+\frac{p^4}{960\pi^2M^4}+\cdots .
$$
$-\Tr(G_0UG_0U)/4$ に代入すると四微分の作用は $-\int(\partial^2 U)^2/(3840\pi^2M^4)$。正負の交代は元の $I$ の対数級数に由来する。Fourier 波数を $|p|=M/2$ とした場合、この四微分項は二微分項の $-1/40$。一方、元の泡の解析接続は時間的な $p^2=-4M^2$ に切断を持つため、任意に大きい外部運動量にこの微分展開を延長できない。

### 背景の共変な二次演算子を直接作る

行列値の接続を $A=\bar A+a$ とする。共変微分が $D'_\mu=U D_\mu U^{-1}$ になるように $\bar A$ と全ての $A$ に同じ接続の変換則を選べば、その差は $a'=UaU^{-1}$。$\bar D_\mu a_\nu=\partial_\mu a_\nu-ig[\bar A_\mu,a_\nu]$ であり、$\partial_\mu U$ を含む項は接続の変分で消えるため $\bar D a$ も随伴に変換する。曲率を二次まで展開すると
$$
F_{\mu\nu}(A)=\bar F_{\mu\nu}+\bar D_\mu a_\nu-\bar D_\nu a_\mu-ig[a_\mu,a_\nu].
$$
これは $[D_\mu,D_\nu]=-igF_{\mu\nu}$ を展開したもの。Euclid の作用 $\tfrac12\int\tr F_{\mu\nu}F_{\mu\nu}$ の二次には、線形揺らぎの二乗と $-ig\,\tr(\bar F_{\mu\nu}[a_\mu,a_\nu])$ が現れる。背景共変なゲージ固定は $\xi^{-1}\int\tr(\bar D_\mu a_\mu)^2$ で、これも随伴変換のトレースだから不変。従って量子場の二次演算子は背景変換で共役され、有限な正則化が同じ共変性を保つなら行列式も不変である。

例えば背景が $\bar F=0$ で $\bar D_\mu=U\partial_\mu U^{-1}$ の形に取れる局所領域では、境界条件まで同じ変換で移せるとき、二次演算子を自由背景のものへ共役できる。固有値が等しいので背景を変えない差を引いた行列式は零である。局所的に曲率が零でも、大域的な境界条件を別のものへ変える変換なら、同じスペクトルを仮定できない。非可換の量子場に必要なゲージ条件と行列式の詳しい計算は [非可換ゲージ理論の量子化](/physics/quantum-field-theory/gauge-quantization/) に続く。

### 微分展開から非局所応答へ戻る

背景 $U(x)$ の二次で厳密な一ループの泡の差は、外部運動量 $p$ に対する $I(p)-I(0)=-\int_0^1\log[1+x(1-x)p^2/M^2]\dd x/(16\pi^2)$。微分展開の係数は、$p^2$ の $n$ 次について
$$
\frac{(-1)^n}{16\pi^2n M^{2n}}
\,\frac{(n!)^2}{(2n+1)!}.
$$
$n$ 次の係数の積分を $n$ 回部分積分すると Beta の値を得る。$|p^2|<4M^2$ なら対数級数は一様に積分でき、数学上のこの Taylor 展開には収束円がある。一方、低エネルギーの有効作用を有限回の微分へ打ち切る場合は、$p^2/(4M^2)\ll1$ のとき精度がよい。例えば $p^2=M^2$ では第二項の第一項に対する比は絶対値 $1/10$、第三項と第二項の比は $1/7$。一様な有効ポテンシャルはこの差に $p=0$ を入れた値しか含まないため、曲がりのある背景の応答には微分項が必要である。

:::example[微分項から空間的な変調の重みを測る]
周期箱で $U(x)=u\cos(kx_1)$ を取り、$k\ne0$ を箱に合う波数とする。$\int U^2/V=u^2/2$、$\int(\partial U)^2/V=k^2u^2/2$。従って二微分の一ループのエネルギー密度差は $k^2u^2/(768\pi^2M^2)$。四微分項は $-k^4u^2/(7680\pi^2M^4)$。同じ局所分布の振幅 $u$ でも変調の波長を半分にすると先頭の微分補正は四倍になり、高次の補正比も四倍になる。従って同じ一様ポテンシャルの式を、空間の変化が大きい背景に適用できるとは限らない。
:::

### 一定磁場の行列式を Landau 準位で求める

複素スカラー一種類が質量 $M>0$、符号付き電荷 $q$ を持ち、一定の Euclid 磁場 $F_{12}=B$ に結合するとする。ゲージを $A_2=Bx_1$、他を零に取る。$x_2$ の Fourier 運動量 $k_2$ を固定すれば
$$
-D^2+M^2=-\partial_1^2+(k_2+qBx_1)^2+p_3^2+p_4^2+M^2 .
$$
$x_1$ の核は中心 $-k_2/(qB)$ の調和振動子で、固有値は $(2n+1)|qB|$、$n\ge0$。大きな矩形 $L_1L_2$ では $k_2$ が $2\pi/L_2$ ごとに違い、中心が箱に入る $k_2$ の幅が $|qB|L_1$ なので、単位面積あたりの縮退は $|qB|/(2\pi)$。残る二方向の Gaussian の積分 $1/(4\pi s)$ と準位の幾何級数 $1/[2\sinh(|qB|s)]$ を掛けると
$$
\frac1V\Tr e^{-s(-D^2+M^2)}
=\frac{e^{-sM^2}}{(4\pi s)^2}\,
\frac{|qB|s}{\sinh(|qB|s)} .
$$
$B\to0$ では最後の因子が一になり、四次元の自由熱核を回復する。

複素ボースの積分では $\Gamma_1=\Tr\log(-D^2+M^2)$。Schwinger 表示で基準 $B=0$ の値と二次の電荷項を引けば、有限な残りは
$$
\frac{\Delta\Gamma_1}{V}
=-\frac1{16\pi^2}\int_0^\infty
\frac{e^{-sM^2}}{s^3}\!
\left[\frac{|qB|s}{\sinh(|qB|s)}-1+\frac{(qBs)^2}{6}\right]\dd s .
$$
対数の差を指数の差へ移すときの負号を残した。$|qB|\ll M^2$ では、括弧の先頭が $7(qBs)^4/360$。積分で $\int_0^\infty s e^{-sM^2}\dd s=1/M^4$ より
$$
\frac{\Delta\Gamma_1}{V}=-\frac{7(qB)^4}{5760\pi^2 M^4}
+O((qB)^6/M^8).
$$
例えば $B$ を小さいまま二倍にするとこの補正は十六倍となる。二次は電荷の規格化に移し、四次は局所な非線形応答として残る。一定磁場では Hessian に負固有値がないため実の行列式を使える。時間的な電場へ式を続ける場合は、指数の核の極の処方を別に指定する。

## 演習問題

::::exercise{#exr-qft-effective-action-1 level=basic}
一変数の $W(J)=J^2/(2K)$、$K>0$ の Legendre 変換を求めよ。

:::solution
$\Phi=J/K$ より $J=K\Phi$。$J\Phi-W=K\Phi^2/2$ となる。
:::
::::

::::exercise{#exr-qft-effective-action-2 level=basic}
分散 $v$ の中心化 Gaussian の四次連結量を求めよ。

:::solution
四次モーメントは $3v^2$、二点ずつの分割が三つあるので差は零である。
:::
::::

::::exercise{#exr-qft-effective-action-3 level=basic}
二個の等しい固有値 $M^2$ だけを持つ有限の Hessian で一ループ行列式の項を求めよ。

:::solution
$\frac12\Tr\log H=\frac12(2\log M^2)=\log M^2$ である。基準を選べば対数の引数を無次元にできる。
:::
::::

::::exercise{#exr-qft-effective-action-4 level=basic}
$\phi^4$ のポテンシャル $\lambda\varphi^4/24$ の Hessian の質量項を求めよ。

:::solution
二回微分で $\lambda(4\cdot3)\varphi^2/24=\lambda\varphi^2/2$ である。
:::
::::

::::exercise{#exr-qft-effective-action-5 level=advanced}
$W^{(3)}$ の公式を一変数で導け。

:::solution
$W''=1/\Gamma''(\Phi)$ かつ $\partial_J\Phi=W''$ である。連鎖律で $W'''=-\Gamma'''(\partial_J\Phi)/(\Gamma'')^2=-\Gamma'''(W'')^3$ となる。
:::
::::

::::exercise{#exr-qft-effective-action-6 level=advanced}
Coleman–Weinberg 形で $V(v)-V(0)$ を求めよ。

:::solution
$\varphi^4\log\varphi^2\to0$ なので $V(0)=C$。$V(v)=C-Bv^4/2$ から差は $-Bv^4/2$ である。
:::
::::

::::exercise{#exr-qft-effective-action-7 level=advanced}
$O(2)$ Ward 恒等式から非零停留点で角方向 Hessian 固有値が零であることを示せ。

:::solution
$\partial_i V(T\Phi)^i=0$ を $\Phi^j$ で微分する。停留点で $\partial_i V=0$ だから $V_{ji}(T\Phi)^i=0$。非零の $\Phi$ では $T\Phi$ も非零である。
:::
::::

::::exercise{#exr-qft-effective-action-8 level=advanced}
一様背景の例の半径方向 Hessian を求めよ。

:::solution
$R=|\Phi|$ とすれば $\Phi\Phi^{\mathsf T}\Phi=R^2\Phi$ なので半径方向固有値は $\lambda(R^2-v^2)+2\lambda R^2=\lambda(3R^2-v^2)$ である。
:::
::::

::::exercise{#exr-qft-effective-action-9 level=advanced}
二相だけを持つ体積 $\mathcal V$ の模型で $w(J)=\log[\cosh(\mathcal V vJ)]/\mathcal V$ とする。有限体積の Legendre 変換を $|\Phi|<v$ で求め、$\Phi=0$ での曲率を計算せよ。

:::solution
$\Phi=v\tanh(\mathcal V vJ)$ なので $r=\Phi/v$ と置くと $J=\operatorname{arctanh}(r)/(\mathcal V v)$。また $\cosh(\operatorname{arctanh}r)=1/\sqrt{1-r^2}$。これを $\Gamma=J\Phi-w$ に入れると
$$
\Gamma(\Phi)=\frac1{\mathcal V}\left[r\operatorname{arctanh}r+\tfrac12\log(1-r^2)\right].
$$
微分すると角括弧の微分は $\operatorname{arctanh}r$、二階は $1/(1-r^2)$。従って $\Gamma''(0)=1/(\mathcal V v^2)>0$。$w''(0)=\mathcal V v^2$ の逆にもなっている。有限体積では凸であり、$\mathcal V\to\infty$ で区間内の全ての曲率が点ごとに零へ向かう。体積と源の極限を区別すれば、有限体積の解析的な応答と平坦な無限体積極限がつながる。
:::
::::

::::exercise{#exr-qft-effective-action-10 level=advanced}
一変数の反転対称な既約作用で $\Gamma''(0)=K>0$、$\Gamma^{(4)}(0)=\lambda$、$\Gamma^{(6)}(0)=0$ とする。六次の連結相関が零かどうか、理由とともに示せ。

:::solution
$\Gamma'(\Phi)=K\Phi+\lambda\Phi^3/6+O(\Phi^7)$ を逆に解く。$\Phi=J/K-\lambda J^3/(6K^4)+bJ^5+\cdots$ を代入すると、$J^5$ の係数の条件は $Kb-\lambda^2/(12K^6)=0$。従って $b=\lambda^2/(12K^7)$、$\Phi=W'$ なので $W^{(6)}(0)=5! b=10\lambda^2/K^7$。既約な六点は零でも、二つの四点が一本の完全二点でつながる連結な木がある。外部六本と内部一本を合わせた $K^{-7}$、三本ずつへ割り付ける十通りという規格化もこの答えと一致する。
:::
::::

::::exercise{#exr-qft-effective-action-11 level=advanced}
背景付きの質量核 $H=\begin{pmatrix}4&1\\1&3\end{pmatrix}M^2$ を持つ二つの実スカラーで、二成分を同時に対角化する計算と重い成分を先に積分する計算が同じ零運動量の行列式を与えることを確かめよ。

:::solution
同時に対角化すれば固有値は $M^2(7\pm\sqrt5)/2$。積は $11M^4$。第二成分を先に積分すると行列式の因子は $3M^2$、第一成分の残りの核は $4M^2-M^4/(3M^2)=11M^2/3$。積はやはり $11M^4$ で、一ループ作用の寄与 $\tfrac12\log(11M^4)$ に一致する。基準の行列式で割れば対数の引数を無次元にできる。先に消した成分の行列式も残す必要があり、残った第一成分の核だけを用いると二自由度を積分した結果にはならない。
:::
::::

::::exercise{#exr-qft-effective-action-12 level=advanced}
三次元の空間と周期の Euclid 時間の実スカラーで、有限温度の零周波数モードが一ループに与える $M$ の三次の係数を再計算せよ。質量二乗 $M^2$ について零で減算してよい。

:::solution
$M^2$ の微分は $T/2$ と $\int\dd^3 k/(2\pi)^3[(k^2+M^2)^{-1}-(k^2)^{-1}]$ の積。球積分にすると $-TM^2/(4\pi^2)\int_0^\infty\dd k/(k^2+M^2)=-TM/(8\pi)$。$M^2$ を零から $M^2$ まで積分すれば $-TM^3/(12\pi)$ である。$M^2$ の Taylor 展開を先に取ると $k\to0$ で個々の項が発散するため、この三次の冪を求められない。零周波数だけを分けて積分してから小質量へ進む順が必要である。
:::
::::

## 確認テスト

::::quiz
連結相関を生成するものを選べ。

- [x] $\log Z_E$
- [ ] $1/Z_E$
- [ ] $\det Z_E$

:::explanation
独立な連結成分の積を全て足す組合せ論は指数形を与える。その対数を取ると成分一個の寄与が残る。
:::
::::

::::quiz{answer="8" tol="1e-6"}
本文の Coleman–Weinberg 形で $V''(v)/(Bv^2)$ を答えよ。

:::explanation
$V'=4B\varphi^3\log(\varphi^2/v^2)$ を微分する。$\varphi=v$ で対数部分が零になり、対数の微分 $2/v$ から $8Bv^2$ を得る。
:::
::::

::::quiz
Euclid の複素 Grassmann Gaussian が有効作用に与える項を選べ。

- [ ] $+\frac12\Tr\log K$
- [x] $-\Tr\log K$
- [ ] $+\Tr K^{-1}$

:::explanation
積分は $\det K$ であり、有効作用では指数の $e^{-S_E}$ に合わせて積分結果の負の対数を取る。
:::
::::

::::quiz{answer="2" tol="1e-6"}
$W(J)=J^2/4$ のとき $\Gamma''$ を答えよ。

:::explanation
$W''=1/2$ である。Legendre 変換の二階微分はその逆なので二になる。
:::
::::

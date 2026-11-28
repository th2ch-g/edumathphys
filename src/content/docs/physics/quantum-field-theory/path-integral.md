---
title: "経路積分による定式化"
description: "時間分割を起点に場の生成汎関数を構成し、フェルミオン、電磁場、Euclid 時間へ拡張する。"
notation:
  - symbol: 'Z[J]'
    meaning: 源つき真空汎関数
  - symbol: '\eta,\bar\eta'
    meaning: Grassmann 奇の源
  - symbol: 'S_E'
    meaning: Euclid 作用
---

## この章の目標

- 短時間核を合成して経路積分を作る。
- 自由場の Gaussian 汎関数から相関と摂動展開を求める。
- フェルミオンと電磁場の積分および場の変数変換を扱う。

自然単位系 $\hbar=c=1$ を用いる。実時間の場の計量は $\eta=\diag(+1,-1,-1,-1)$、Euclid 時間を用いる箇所では $\delta_{\mu\nu}$ とする。

## 量子力学の経路積分の復習

### 時間を分割して完全性を挿入する

$H=p^2/(2m)+V(q)$、$T=N\delta t$ とする。積公式で短時間の運動項とポテンシャル項に分け、座標と運動量の完全系を交互に挿入する。$\braket{q}{p}=e^{ipq}/\sqrt{2\pi}$ を使うと各区間の核は
$$
\int\frac{\dd p}{2\pi}\exp\left[ip(q_{j+1}-q_j)-i\delta t
\bigl(p^2/(2m)+V(q_j)\bigr)\right].
$$
$p$ の平方完成で $(m/(2\pi i\delta t))^{1/2}\exp[i m(q_{j+1}-q_j)^2/(2\delta t)-i\delta t V(q_j)]$ を得る。微小な負の虚時間を付けて Gaussian の枝を決める。中間の $q_j$ を全て積分し $N\to\infty$ とするものが :term[経路積分]{en="path integral" yomi="けいろせきぶん"} である。自由粒子なら Fourier 積の位相が $e^{-iTp^2/(2m)}$ に合成され、核の同じ式で $\delta t$ を $T$ に替えたものになる。一般の運動量依存相互作用では、座標と運動量の順序を指定してから時間分割する必要がある。短時間核の合成は :ref[def-aqm-propagator] の伝播核の考え方である。

### 調和振動子の Euclid 核

$S_E[q]=\int_0^T m(\dot q^2+\omega^2q^2)/2\,\dd\tau$ とする。端点 $q_i,q_f$ をつなぐ古典解は
$$
q_c(\tau)=\frac{q_i\sinh[\omega(T-\tau)]+q_f\sinh(\omega\tau)}{\sinh(\omega T)} .
$$
$q=q_c+\eta$、$\eta(0)=\eta(T)=0$ と置くと交差積は運動方程式で消え、揺らぎの積分は端点に依らない。$S_E[q_c]=m[q_c\dot q_c]_0^T/2$ なので指数は
$$
-\frac{m\omega}{2\sinh(\omega T)}[(q_i^2+q_f^2)\cosh(\omega T)-2q_iq_f].
$$
揺らぎの係数は時間合成の Gaussian 積分と $T\downarrow0$ の自由核規格化で $\sqrt{m\omega/(2\pi\sinh(\omega T))}$ となる。実際、係数 $m\omega[\coth(\omega T_1)+\coth(\omega T_2)]/2$ の Gaussian を中間端点で積分すると、$\sinh(\omega(T_1+T_2))$ の加法公式でこの規格化が保たれる。時間を一周させて端点を積分すると $1/[2\sinh(\omega T/2)]$、エネルギー $\omega(n+1/2)$ の Boltzmann 重みを足したものと一致する。

### 時間片の核を畳み込む

自由粒子の実時間核は $K_T(x,y)=\sqrt{m/(2\pi iT)}\exp[im(x-y)^2/(2T)]$ とする。根の位相は $T$ に小さい負の虚部を入れて定める。$K_{T_2}(x,z)K_{T_1}(z,y)$ の指数を $z$ について平方完成すると中心は $(T_1x+T_2y)/(T_1+T_2)$、二次の係数は $m(T_1+T_2)/(2T_1T_2)$ である。Gauss 積分の因子が二つの規格化を一つにまとめ、定数部分は $m(x-y)^2/[2(T_1+T_2)]$ になる。従って畳み込みは $K_{T_1+T_2}$ である。例えば $T_2=2T_1$ では中心は $(x+2y)/3$ であり、古典直線経路を時刻 $T_1$ で評価した値に一致する。

### 時間分割の順序と適用条件

実時間の座標間の核は通常の確率密度と異なり、振動する位相を持つ。例えば $\delta t>0$ で運動項だけを積分するときは、$\delta t$ を $\delta t-i\epsilon$ へ置いた収束する運動量積分から始め、最後に $\epsilon\downarrow0$ を取る。この処方で得る平方根は正の実時間に対して $1/\sqrt i=e^{-i\pi/4}$ を含む。$\delta t\to0$ での極限は、滑らかな波束に核を作用させて単位作用素へ強収束するという意味である。

$T=p^2/(2m)$ と $V(q)$ が可換でないため、一片の $e^{-i\delta t T}e^{-i\delta t V}$ は $e^{-i\delta t(T+V)}$ と一般に一致しない。有限次元の切断の下では両者を二次まで展開でき、その差は $-(\delta t)^2[T,V]/2$。有限時間を $N$ 片に分ければ差のノルムは $N O((\delta t)^2)\to0$ になる。実際、積の差を一片ずつ取り替えた和にすると各因子はユニタリでノルム一である。連続系で同じ結論を使うには自己共役性および積公式の条件が要る。例えば有界な実 $V$ と自己共役な自由 $T$ の場合、$T+V$ は $T$ の定義域で自己共役で、積公式が強収束する。一般の特異ポテンシャルで項別のノルム展開を用いることはできない。作用の時間格子は、どの演算子を時間片のどこへ置いたかも記録する。

:::example[位置に線形な力の短時間積]
$V(q)=-f q$ とし、密な共通定義域の上で $[T,V]=ifp/m$ を得る。実際、$[p,q]=-i$ から $[p^2,q]=-2ip$。時間片の順序を反転した積では二次の交換子が逆符号になる。従って一片の誤差を抑えたい場合は $e^{-i\delta t V/2}e^{-i\delta t T}e^{-i\delta t V/2}$ を使える。指数を二次まで掛けると交差項が $-(\delta t)^2(TV+VT)/2$ になり、元の時間発展の二次と一致する。無界演算子のこの展開は共通の適切な波束上で読む。
:::

### 位相と揺らぎの規格化を独立に読む

調和振動子の Euclid 核について $\omega T\gg1$ では、$\sinh(\omega T)\simeq e^{\omega T}/2$、$\coth(\omega T)\simeq1$ となる。従って核は $(m\omega/\pi)^{1/2}e^{-\omega T/2}e^{-m\omega(q_f^2+q_i^2)/2}$。指数の $q$ 依存は基底状態の波動関数の積であり、$e^{-\omega T/2}$ は揺らぎの規格化から出ている。一方 $\omega\to0$ で $T$ を固定すると、$\sinh(\omega T)\sim\omega T$ から $\sqrt{m/(2\pi T)}\,e^{-m(q_f-q_i)^2/(2T)}$。この自由核を一周させて座標の全実線を積分すると空間体積で発散するので、振動子の有限なトレースの $\omega\downarrow0$ にも発散がある。

### 運動量を含む Hamiltonian の離散化

自由運動項以外にも運動量が現れる場合、古典の $H(q,p)$ だけでは各短時間核を決められない。演算子の順序を指定し、短時間核を $\langle q'|H|q\rangle$ から作る。例えば $H=(f(q)p+p f(q))/2$ なら
$$
\langle q'|H|q\rangle
=\frac{f(q')+f(q)}2
\int\frac{\dd p}{2\pi}\,p\,e^{ip(q'-q)}
$$
である。$f$ が滑らかなら、$(f(q')+f(q))/2=f((q'+q)/2)+O((q'-q)^2)$。運動量の積分はデルタの一階微分なので、二次の差を掛けても分布として零になる。従ってこの対称順序は中点で $f$ を評価する短時間核と対応する。一方 $f(q)p$ を選ぶと $f(q')$ が掛かり、$p f(q)$ では $f(q)$ が掛かる。二つの差は $[f(q),p]=if'(q)$ に対応する。場の作用を時間連続で書くときにも、時間格子上の同時刻の積はこの順序に依存する。例えば $f(q)=q$ の対称化では $(qp+pq)/2=qp-i/2$ であり、順序の変更が定数を生む。

### 揺らぎの行列式を格子の漸化式で求める

調和振動子の Euclid 時間 $T=N\delta\tau$ を $N$ 片に分け、両端を固定した揺らぎを $\eta_0=\eta_N=0$ とする。離散作用の二次核は、共通の $m/\delta\tau$ を除いて、対角が $2+(\omega\delta\tau)^2$、隣接要素が $-1$ の $(N-1)$ 次行列である。各時刻に $\eta_j^2$ が隣接する二つの時間差から出るため、対角の運動項は $2$ となる。この行列の $j$ 次の先頭小行列の行列式を $D_j$ とすると、端の行の余因子展開から
$$
D_j=[2+(\omega\delta\tau)^2]D_{j-1}-D_{j-2},\qquad
D_0=1,\quad D_1=2+(\omega\delta\tau)^2 .
$$
$2\cosh u=2+(\omega\delta\tau)^2$ とすれば、解は $D_j=\sinh[(j+1)u]/\sinh u$。初期値と $\sinh$ の加法公式で確認できる。自由核の規格化の由来となる $\omega=0$ の行列式は $D_{N-1}=N$。従って振動子の揺らぎ積分と自由な揺らぎ積分の比は $\sqrt{N\sinh u/\sinh(Nu)}$。$\delta\tau\to0$ で $u=\omega\delta\tau+O((\delta\tau)^3)$ なので比は $\sqrt{\omega T/\sinh(\omega T)}$ となる。自由粒子の係数 $\sqrt{m/(2\pi T)}$ を掛けて、先に得た調和振動子の規格化を回復する。古典解の作用は端点に依存するが、この行列式は Dirichlet の揺らぎだけを数えるため端点の値に依存しない。

## スカラー場の経路積分と生成汎関数

箱と格子でスカラー場を有限個の座標とみなす。各時刻の $q_j$ を格子点ごとの $\phi_j(\vb x)$ で置き換え、時空の各点に源 $J(x)\phi(x)$ を加える。時間の両端で真空を射影し、源なしの値で割ると
$$
Z[J]=\frac{\int\mathcal D\phi\,e^{iS[\phi]+i\int J\phi}}
{\int\mathcal D\phi\,e^{iS[\phi]}},\qquad
G_n=\left.\frac1{i^n}\frac{\delta^n Z}{\delta J(x_1)\cdots\delta J(x_n)}\right|_{J=0}.
\label{eq-qft-path-generator}
$$
挿入は時間分割された時間発展の中に置かれるため $G_n$ は時間順序積である。有限の時空格子や虚時間で積分の意味を定め、連続実時間の式には $i0$ を保つ。例えば自由真空の一点関数は $\phi\mapsto-\phi$ で奇なので零である。

### 有限格子で生成汎関数を定義する

実 Euclid 格子の点数を $N$ とし、正定値の自由二次形式を $K$ とする。格子間隔の因子を $K,J$ に含めて $Z[J]=\int_{\R^N}\dd^N\phi\,e^{-\phi^\mathsf T K\phi/2+J^\mathsf T\phi}$ と書く。$\phi\mapsto\phi+K^{-1}J$ で平方完成すれば $Z[J]/Z[0]=e^{J^\mathsf T K^{-1}J/2}$。二つの点で $K=\begin{pmatrix}3&-1\\-1&3\end{pmatrix}$ なら $K^{-1}=\begin{pmatrix}3&1\\1&3\end{pmatrix}/8$ である。例えば源 $J=(h,0)$ は平均 $\langle\phi_1\rangle=3h/8$ だけでなく $\langle\phi_2\rangle=h/8$ も作る。逆行列の非対角成分が、源のない格子点への応答を表す。

### 空間格子から作用と真空境界を得る

間隔 $a$ の周期空間格子で、位置 $j$ の場を $\phi_j$、連続場の運動量密度に $a^3$ を掛けた共役を $P_j$ とする。離散 Hamiltonian は
$$
H=\sum_j\left[
 \frac{P_j^2}{2a^3}+
 \frac{a^3}{2}\left(m^2\phi_j^2+\sum_i\frac{(\phi_{j+\widehat i}-\phi_j)^2}{a^2}\right)
 +a^3 V_I(\phi_j)\right].
$$
時間片ごとに $P_j$ を積分すると、各点のガウス指数は $ia^3(\phi_{j,t+\delta t}-\phi_{j,t})^2/(2\delta t)$。従って離散的な時間微分の二乗と空間微分の二乗は Minkowski 作用では逆符号で現れる。空間の差分の Fourier 固有値は $\sum_i 4\sin^2(k_i a/2)/a^2$。例えば一次元の $k a=\pi$ のモードでは固有値 $4/a^2$。低い波数では $k^2$ へ近づくが、格子端では連続表示の $\pi^2/a^2$ をそのまま使えない。

有限体積で $H$ に孤立した基底状態 $|\Omega\rangle$ があり、端点に置く試験状態 $|B\rangle$ が $\langle\Omega|B\rangle\ne0$ を満たすとする。挿入を含む有限時間区間の両側へ $e^{-T H}$ を掛け、同じ源なしの行列要素で割る。スペクトル展開では $e^{-TH}|B\rangle=e^{-TE_0}(|\Omega\rangle\langle\Omega|B\rangle+O(e^{-T\Delta}))$。従って両端の規格化因子と $e^{-TE_0}$ が比で消え、$T\to\infty$ で真空期待値を得る。ここで $\Delta>0$ は有限系でのギャップである。実時間では偏角をわずかに負へ傾けた時間発展に対応し、端点をまず無限遠へ送る。源を置いた有限区間は、この射影を行う区間の外側と区別する。

:::example[有限格子で零モードを検出する]
上の自由系で $m=0$ とすると、全格子点が同じ値を取る一様な場の空間エネルギーは零。Euclid で時間も周期的なら一様な時空定数は二次形式の零固有ベクトルになる。規格化の分母を単に無限の Gaussian として扱うと共分散が定まらない。例えば空間平均を零に制約すればこの方向を積分から除ける。別に小質量 $m>0$ を残すと、規格化された一様モードの分散は $1/m^2$ に比例し、全場の質量ゼロ極限を無条件に取れないことが分かる。
:::

### 時間順序と因果応答が異なる境界値問題

古典源を作用に加えると、真空から真空への振幅を源で微分した量は、時間順序 Green 関数に従う。一方、ある初期状態の中で後の時刻の期待値を測るなら、時間発展の両側で源が作用する。例えば自由な一振動子 $H_0=(p^2+E^2 q^2)/2$、$H_I(t)=-J(t)q(t)$ とする。一次の状態の変化 $U_I=1+i\int\dd t'\,J(t')q(t')+\cdots$ を $U_I^\dagger q(t)U_I$ へ入れると
$$
\delta\langle q(t)\rangle
=i\int_{t_0}^t\dd t'\,J(t')\langle[q(t),q(t')]\rangle .
$$
$q(t)=(a e^{-iEt}+a^\dagger e^{iEt})/\sqrt{2E}$ だから、$i[q(t),q(t')]=\sin(E(t-t'))/E$。従って応答は遅延核を持つ。例えば $J(t)=j\delta(t)$ で $t_0<0$ なら $\delta\langle q(t)\rangle=j\theta(t)\sin(Et)/E$。これに対して、真空間の $Z[J]$ の一微分は $ij e^{-iE|t|}/(2E)$。同じ $J$ でも、前者は初期状態の時間発展、後者は未来の真空への重なりで規格化した行列要素である。物理的な因果応答を求めるときは、一方を他方へ置き換えられない。

正則化した作用でも二つの計算は同じ時間片から作れる。初期状態の行列要素では前向きの $U_I$ と後向きの $U_I^\dagger$ の間に観測量を置くので、前向きと後向きの経路をそれぞれ積分する。外部源を二本に別々につけると、順序の異なる相関を区別できる。同じ源を両側で使い、観測量が単位なら $U_I^\dagger U_I=1$ になる。これは真空間の $Z[J]/Z[0]$ が同じ源に対して位相や真空遷移を持つのと異なる規格化である。

## 自由場の生成汎関数と摂動展開

部分積分して $S_0=\frac12\phi K\phi$、$K=-(\Box+m^2)+i0$ と書く。積分変数を $\phi+K^{-1}J$ へずらすと $\frac12\phi K\phi+J\phi=\frac12(\phi+K^{-1}J)K(\phi+K^{-1}J)-\frac12 JK^{-1}J$。$\Delta_F=iK^{-1}$ だから
$$
Z_0[J]=e^{-\frac12 J\Delta_F J}.
\label{eq-qft-free-generating-functional}
$$
二微分で $\Delta_F$、四微分で三つのペアリングが出る。$\exp(iS_I[(1/i)\delta/\delta J])$ を $Z_0[J]$ に作用させれば相互作用の級数になる。各微分が $e^{i\int J\phi}$ から場を取り出すためである。最後に $J=0$ での値で割れば真空泡が相殺する。例えば四次の一次で四つの外点がすべて一頂点へ結ばれる項は :eqref[eq-qft-phi4-position] と一致する。

### 分母を展開して真空因子を除く

単一実変数で $Z_g(J)=\int\dd x\,e^{-x^2/(2D)-gx^4/4!+Jx}$ とする。$g\downarrow0$ の一次係数は Gauss の可積分なモーメントで微分して求められる。$\langle x^2\rangle_g$ の分子と分母を一次で割ると
$$
\langle x^2\rangle_g=D-\frac{g}{24}\left(\langle x^6\rangle_0-\langle x^2\rangle_0\langle x^4\rangle_0\right)+\order{g^2}
=D-\frac{g D^3}{2}+\order{g^2}.
$$
$\langle x^{2n}\rangle_0=(2n-1)D\langle x^{2n-2}\rangle_0$ は部分積分で従う。$15D^3-3D^3$ の差の $12$ は、二外点を同じ頂点の異なる脚につなぐ選択数 $4\cdot3$ である。従って有限積分の規格化が場の真空泡の除去を再現する。

### 正規化された連結四点の最低次

Euclid の自由共分散を $D$ とし、$\lambda\int\phi^4/4!$ を一次で加える。四つの外点を $x_1,\ldots,x_4$ とした完全四点の一次の分子には $-\lambda\int\dd^4z\,\langle\prod_i\phi(x_i)\phi(z)^4\rangle_0/4!$ が現れる。外部点が全て $z$ の異なる脚へつながる選択数は $4!$。従って連結四点は
$$
-\lambda\int\dd^4z\, D(x_1-z)D(x_2-z)D(x_3-z)D(x_4-z)+O(\lambda^2).
$$
二つの外点だけを $z$ へ結び、残る外部二点を直接結ぶ項は、補正された二点関数と自由な二点関数の積である。四つの $z$ の脚を互いに縮約する項は全外点と切れており規格化の分母が消す。従って分母での真空泡の除去と、四点から二点の積を引く連結化は区別する。例えば対称な一変数の $D=1$ の模型では一次の二点は $1-\lambda/2$、完全四点は $3-4\lambda$ であり、連結四点は $(3-4\lambda)-3(1-\lambda)= -\lambda+O(\lambda^2)$ である。

### 古典源で自由場の応答を計算する

$Z_0[J]=\exp[-J\Delta_F J/2]$ を一回微分し、$Z_0[J]$ で割ると $\langle\phi\rangle_J=i\Delta_F J=-K^{-1}J$。従って $(\Box+m^2)\langle\phi\rangle_J=J$ になる。ただしこの式の境界条件は Feynman の処方である。例えば固定空間運動量、一時刻の源 $J(t)=j\delta(t)$ なら平均は $ij e^{-iE|t|}/(2E)$。$t<0$ でも零ではなく、後で測った状態で規格化する真空間の行列要素であることが現れる。因果的な初期値応答では遅延 Green 関数 $\theta(t)\sin(Et)/E$ が現れる。微分方程式の右辺が同じでも核は境界条件で区別される。

### 二点関数の二次補正から独立な組合せを分ける

一変数の正則化された $\phi^4$ 積分で共分散を $D$ とし、モーメントを $M_{2n}=(2n-1)!!D^n$ と書く。$e^{-\lambda\phi^4/24}$ を二次まで展開した分子と分母は
$$
\begin{aligned}
N&=M_2-\lambda M_6/24+\lambda^2 M_{10}/(2\cdot24^2),\\
Z&=1-\lambda M_4/24+\lambda^2 M_8/(2\cdot24^2).
\end{aligned}
$$
従って $N/Z$ の二次の係数は $[M_{10}-M_2M_8-2M_6M_4+2M_2M_4^2]/(2\cdot24^2)$。数値を入れると $(945-105-90+18)D^5/(2\cdot24^2)=2D^5/3$。一回補正した逆核が $\Sigma_1=\lambda D/2$ なら、伝播関数の逆を展開した $D\Sigma_1D\Sigma_1 D$ は $\lambda^2D^5/4$。従って残る二次の逆核は $\Sigma_2=-(2/3-1/4)\lambda^2 D^3=-5\lambda^2 D^3/12$。この例には運動量はないが、発生する Wick 縮約の数と逆核へ移す符号は場のグラフと同じ手順で分けられる。

### 自由な核を二次の相互作用で照合する

正定値の有限核 $K$ に対し、小さい対称行列 $U$ を加える。全体 $K+U$ が正定値なら、相互作用として $-\phi^\mathsf T U\phi/2$ を指数に展開した積分は、厳密な共分散 $(K+U)^{-1}$ を再現する。実際、$D=K^{-1}$ とすれば
$$
(K+U)^{-1}=D-DUD+DUDUD-\cdots .
$$
各項は内部で $U$ を一回ずつ挿入した連結鎖である。源と関係しない閉じた鎖は $\det(1+DU)^{-1/2}$ の展開を作り、同じ積分の $J=0$ の値で割ると消える。例えば $K=\diag(2,2)$、$U=\begin{pmatrix}0&v\\v&0\end{pmatrix}$、$|v|<2$ なら $(K+U)^{-1}_{12}=-v/(4-v^2)=-v/4-v^3/16+\cdots$。非対角の応答には奇数回の挿入が必要で、対角の補正には偶数回が必要である。外部源の微分からでも、一つの内部線を逆核に挿入する方法からでも同じ選択則が得られる。

:::example[場の一点にだけ二次の欠陥を入れる]
$U=u\,e_a e_a^\mathsf T$、$u$ は実数とする。$(K+U)X=1$ を $D$ で左から掛ければ $X=D-u D e_a e_a^\mathsf T X$。$e_a^\mathsf T$ を掛けて解くと $e_a^\mathsf T X=e_a^\mathsf T D/(1+uD_{aa})$。従って $D'_{ij}=D_{ij}-uD_{ia}D_{aj}/(1+uD_{aa})$。$u>0$ なら $D'_{aa}=D_{aa}/(1+uD_{aa})$ となってその点の揺らぎが小さくなる。他点の分散も $D_{ia}$ が非零なら下がる。これは相互作用が一点にあっても、場の伝播が別の点へ応答を運ぶ例である。
:::

## Grassmann 数とフェルミオンの経路積分

:term[Grassmann 数]{en="Grassmann number" yomi="ぐらすまんすう"} は奇変数 $\theta_i$ が $\theta_i\theta_j=-\theta_j\theta_i$ を満たす代数の元である。特に $\theta_i^2=0$。積分は $\int\dd\theta\,1=0$、$\int\dd\theta\,\theta=1$ と定める。この定義では定数の平行移動が積分値を変えない。線形変換 $\theta=M\chi$ の最上次数の係数には $\det M$ が掛かるので、積分測度の変換は逆行列式を与える。

独立な $\psi_i,\bar\psi_i$ の順序を $\int\mathcal D\bar\psi\mathcal D\psi\,e^{-\bar\psi K\psi}=\det K$ となるように固定する。実際、指数の最上次数では各 $\bar\psi_i,\psi_j$ が一度ずつ現れ、奇変数を基準順へ戻す符号が行列式の順列符号になる。Euclid 源を $+\bar\eta\psi+\bar\psi\eta$ とすれば、平方完成で残る因子は $e^{+\bar\eta K^{-1}\eta}$ である。左・右微分の順序を保って $\psi\bar\psi$ の二点を取ると $K^{-1}$ になる。ボース Gaussian の逆平方根との符号・冪の違いは、一ループ有効作用で現れる。

### Grassmann Gaussian と行列式

$N$ 組の独立な奇変数 $\psi_i,\bar\psi_i$ を用い、積分の順序を $\int e^{-\bar\psi\psi}\dd\bar\psi\,\dd\psi=1$ が成り立つように固定する。$e^{-\bar\psi M\psi}$ の最高次数は各 $\bar\psi_i$ を一度ずつ選び、全 $\psi_j$ も一度ずつ選んだ項だけである。並べ替えの符号は列の置換符号で、和が $\det M$ となる。源 $\bar\eta,\eta$ を加えると
$$
-\bar\psi M\psi+\bar\eta\psi+\bar\psi\eta
=-(\bar\psi-\bar\eta M^{-1})M(\psi-M^{-1}\eta)+\bar\eta M^{-1}\eta .
$$
最高次数の係数は奇変数の平行移動で不変だから $Z[\bar\eta,\eta]=\det M\,e^{\bar\eta M^{-1}\eta}$。例えば $M=\begin{pmatrix}2&1\\1&3\end{pmatrix}$ では $\det M=5$、収縮 $M^{-1}_{12}=-1/5$ である。

### 場の二点関数を奇の源から読む

正規化した Euclid の $Z[\bar\eta,\eta]=\exp(\bar\eta K^{-1}\eta)$ を用いる。左微分 $\partial_{\bar\eta_i}$ を先に作用させると $(K^{-1})_{ij}\eta_j$。続いて $\partial_{\eta_j}$ を左から作用させると $(K^{-1})_{ij}$。同じ結果は指数を直接微分しても得られる。最初に $\partial_{\bar\eta_i}$ を作用させると $\psi_i$、次に $\partial_{\eta_j}$ がこの奇因子を通過して $-\psi_i\partial_{\eta_j} e^{\cdots}$ となるが、$\partial_{\eta_j}(\bar\psi_k\eta_k)=-\bar\psi_j$ なので符号が二つ相殺する。逆順の二微分は $-K^{-1}_{ij}$ で、$\bar\psi_j\psi_i$ の順に対応する。このように源の順序を決めると、四点以上でも途中で規約を変えずに計算できる。

:::example[二組の奇変数の交差縮約]
$D=K^{-1}$ とし、$\langle\psi_1\bar\psi_1\psi_2\bar\psi_2\rangle=D_{11}D_{22}-D_{12}D_{21}$。二番目は $\bar\psi_1$ と $\bar\psi_2$ のラベルを入れ替えるため逆符号になる。例えば $K=\begin{pmatrix}2&1\\1&3\end{pmatrix}$ なら $D=\begin{pmatrix}3&-1\\-1&2\end{pmatrix}/5$ なので四点は $(6-1)/25=1/5$。$\det D=1/\det K$ と一致する。個々の縮約が負の値を取っても、組合せの交換符号は縮約の数値の符号とは独立である。
:::

### フェルミオンの一モードのトレース

規格化 $\langle\bar\theta|\theta'\rangle=1+\bar\theta\theta'$ のコヒーレント状態を使い、測度は $\int\dd\bar\theta\dd\theta\,e^{-\bar\theta\theta}=1$ とする。一モードの対角演算子 $O=A|0\rangle\langle0|+B|1\rangle\langle1|$ には $\langle-\bar\theta|O|\theta\rangle=A-B\bar\theta\theta$。重み $1-\bar\theta\theta$ を掛けると最上次数は $-(A+B)\bar\theta\theta$ なので積分は $A+B$。非対角演算子の行列要素は奇変数について一次で、積分は零。よって全線形演算子に対しこの規則がトレースを与える。例えば $O=e^{-\beta\omega a^\dagger a}$ なら $A=1$、$B=e^{-\beta\omega}$ である。端の負号を入れない同じ式では $A-B$ となり、粒子数の偶奇を重みとして入れた別のトレースになる。

### Dirac 演算子の Euclid 化

計量は Euclid で $\delta_{\mu\nu}$ とし、$\gamma_E^4=\gamma^0$、$\gamma_E^j=-i\gamma^j$ と選ぶ。すると全てが Hermite で $\{\gamma_E^\mu,\gamma_E^\nu\}=2\delta^{\mu\nu}$。$t=-i\tau$ なので Dirac の指数から $S_E=\int\dd^4x_E\,\bar\psi(\gamma_E^\mu\partial_\mu+m)\psi$ が得られる。$\bar\psi$ は独立な積分変数として保つ。Fourier 規約 $\partial_\mu\mapsto ip_\mu$ では $(i\gamma_E\cdot p+m)(-i\gamma_E\cdot p+m)=p_E^2+m^2$。従って Euclid の逆核は $(-i\gamma_E\cdot p+m)/(p_E^2+m^2)$ であり、ボース場とは分子が異なる。例えば $\vb p=0$ なら正の Euclid 時間で核は $(1+\gamma_E^4)e^{-m\tau}/2$。$p_4$ の極 $im$ における分子が $m(1+\gamma_E^4)$ を与えるためである。

### コヒーレント状態の時間片から一次の作用を作る

一モードで $a|\theta\rangle=\theta|\theta\rangle$、$\langle\bar\theta|a^\dagger=\langle\bar\theta|\bar\theta$ となる奇変数の状態を選ぶ。$\langle\bar\theta|\theta'\rangle=e^{\bar\theta\theta'}$。規格化された測度で $\int\dd\bar\theta\dd\theta\,e^{-\bar\theta\theta}|\theta\rangle\langle\bar\theta|=1$ は、粒子数零と一の各成分を積分すれば従う。例えば真空成分は $\int e^{-\bar\theta\theta}=1$、粒子数一成分では奇係数を場の基底ベクトルの外に移す符号により、最上次数は $-\bar\theta\theta$ になる。時間片にそれを挿入すると、$H=E a^\dagger a$ では
$$
\langle\bar\theta_{j+1}|e^{-\delta\tau H}|\theta_j\rangle
=1+e^{-\delta\tau E}\bar\theta_{j+1}\theta_j .
$$
実際、真空と一粒子の重みがそれぞれ $1,e^{-\delta\tau E}$。測度の指数を各片の指数に集めると、閉じた時間格子で
$$
S_E=\sum_j\bar\theta_{j+1}(\theta_{j+1}-e^{-\delta\tau E}\theta_j).
$$
端のトレースには反周期の符号が入る。$\theta_{j+1}-\theta_j\simeq\delta\tau\,\partial_\tau\theta$ と $1-e^{-\delta\tau E}\simeq E\delta\tau$ から連続の二次作用は $\bar\theta(\partial_\tau+E)\theta$。フェルミオンの時間運動項が一次になるのは、各時間片の重なりが $\bar\theta_{j+1}\theta_j$ を一回ずつ含むためである。空間の有限モードを Dirac Hamiltonian に従って結合すれば、時間差分の項に空間と質量の行列が加わる。

:::example[反周期符号を行列式から検査する]
$N$ 片で $r=e^{-\beta E/N}$ とする。時間を一つ進める反周期の移動行列 $P$ は $P^N=-1$。行列式は $\det(1-rP)=1+r^N$。実際、$P$ の固有値は $z^N=-1$ の $N$ 個の根なので、$\prod_z(1-rz)$ は $r^N[(1/r)^N+1]$。従って各分割で得た積分は $1+e^{-\beta E}$ となり、トレースの二状態の和と一致する。もし周期境界にすると $P^N=1$ になって $1-e^{-\beta E}$ に変わる。時間上の一つの境界符号が熱的なトレースを変えることが分かる。
:::

### Dirac 行列式を有限なループ展開にする

有限正則化の Dirac 演算子を $K=K_0+V$、自由逆核を $S_0=K_0^{-1}$ とする。場を積分して自由の規格化を除くと $\det(1+S_0V)$ が掛かる。$\|S_0V\|<1$ で対数を展開すれば
$$
\log\frac{\det K}{\det K_0}
=\sum_{n\ge1}\frac{(-1)^{n+1}}n\Tr[(S_0V)^n].
$$
規格化したボースの重み $e^{-S_E}$ の中に取り込むなら、作用にはこの対数の負が加わる。$1/n$ は閉線を巡る基準点の $n$ 通りに対する因子であり、外点を微分してそのラベルを固定すると、同じ一周を重複して数えない。例えば $K_0=m1$、$V$ が二モードで非対角成分 $v$ だけを持つ場合、行列式の比は $1-v^2/m^2$。対数の二次項は $-v^2/m^2$。実ボースの場合に同じ二次核を積分すると、行列式の冪が $-1/2$ になるので、重みの対数は $+v^2/(2m^2)$。同じ演算子の固有値を使っても、Grassmann とボースの積分ではループの符号と多重度が変わる。

### 実条件のある奇変数を積分する

Dirac 場の独立な $\bar\psi,\psi$ と違い、実条件を持つフェルミオンは一組の奇変数 $\chi_1,\ldots,\chi_{2n}$ で書ける。反対称行列 $A=-A^\mathsf T$ に対し、測度を $\int\dd\chi_{2n}\cdots\dd\chi_1\,\chi_1\cdots\chi_{2n}=1$ とする。指数 $e^{\chi^\mathsf T A\chi/2}$ の最高次数の係数は
$$
\operatorname{Pf} A
=\frac{1}{2^n n!}\sum_{\pi\in S_{2n}} \sgn(\pi)
 \prod_{j=1}^{n} A_{\pi(2j-1),\pi(2j)} .
$$
これが :term[Pfaffian]{en="Pfaffian" yomi="ぱふぃあん"} である。指数の $n$ 次だけが最高次数に寄与し、奇変数を標準順へ戻す符号が順列の符号になる。行列変換 $A\mapsto M^\mathsf T A M$ では、最高次数の外積が $\det M$ 倍なので $\operatorname{Pf}(M^\mathsf T A M)=\det M\,\operatorname{Pf} A$。非退化な反対称形式に対して、$A(u,v)\ne0$ の二ベクトルを一対選んで $A(u,v)=1$ に規格化し、残りをその対の直交補空間へ射影する操作を繰り返す。これで対角ブロック $\begin{pmatrix}0&1\\-1&0\end{pmatrix}$ の直和へ合同変換できる。各ブロックの行列式も Pfaffian も一である。変換則により $(\operatorname{Pf} A)^2=\det A$、退化する場合も多項式の連続性で同じである。

:::example[四つの奇変数の符号]
$n=2$ なら $\operatorname{Pf} A=A_{12}A_{34}-A_{13}A_{24}+A_{14}A_{23}$。三つの完全ペアの Wick 符号と一致する。例えば $A_{12}=2,A_{34}=3$ が唯一の非零の上三角要素なら Pfaffian は $6$、行列式は $36$。$A_{12}$ を $-2$ に替えると Pfaffian は $-6$ だが、行列式は同じである。従って Grassmann 積分を単に非負の平方根 $\sqrt{\det A}$ とすることはできず、測度の順序と連続な枝を保つ必要がある。
:::

## 電磁場の経路積分量子化

電磁場の二次演算子には縦方向の零固有値がある。ゲージ条件 $F[A]=\partial\cdot A=f$ を選ぶと $F[A+\partial\omega]=F[A]+\Box\omega$。従って $\int\mathcal D\omega\,\delta(F[A+\partial\omega]-f)|\det\Box|=1$ が成り立つ。ここでは境界条件で $\Box$ の零モードを除く。ゲージ変換で被積分関数と測度が不変なら、変数を移し $\int\mathcal D\omega$ を規格化で除ける。$f$ を $e^{-i\int f^2/(2\xi)}$ で平均すると $-(\partial\cdot A)^2/(2\xi)$ が加わり :eqref[eq-qft-photon-propagator] の逆を得る。$|\det\Box|$ は場に依存しないので規格化で消える。

保存された源 $J$ では $\int J\cdot\partial\omega=0$ であり、Gaussian の縦成分は $k\cdot J=0$ で消える。従って例えば保存源の二次の真空振幅は $\xi$ に依存しない。

### 保存源でゲージパラメータを除く

Euclid の非零モードで $K=k^2P_T+(k^2/\xi)P_L$ とする。有限モードのボース積分は $(\det K)^{-1/2}\exp(J^\mathsf T K^{-1}J/2)$ を与える。$Z[J]/Z[0]$ で行列式が相殺し、$k\cdot J=0$ なら $P_L J=0$ なので指数も $\xi$ に依らない。例えば $k=(2,0,0,0)$、$J=(0,3,4,0)$ なら指数は $25/8$。$J=(1,3,4,0)$ に変えると $\xi/8$ が加わり、非保存源の生成汎関数はゲージに依存する。箱の $k=0$ 成分はこの逆行列の計算から除き、境界条件と残留ゲージを別に扱う。

### 射影演算子で伝播関数を反転する

非零 Euclid 運動量に対して $P_L^{\mu\nu}=k^\mu k^\nu/k^2$、$P_T=1-P_L$ とする。作用の Maxwell 項は $k^2P_T$ だけを持ち、ゲージ固定を加えると $\xi>0$ で
$$
K^{-1}(k)=\frac{P_T}{k^2}+\xi\frac{P_L}{k^2}.
$$
積 $K K^{-1}$ では $P_TP_L=0$、$P_T^2=P_T$、$P_L^2=P_L$ により単位行列を得る。Landau ゲージの $\xi\downarrow0$ はこの計算の後で取れ、伝播関数は横射影に近づく。Minkowski へ戻すときは極の処方を自由光子の伝播関数と合わせる。

:::example[場の強さの相関が $\xi$ に依存しない]
$F_{\mu\nu}(k)=i(k_\mu A_\nu-k_\nu A_\mu)$ と書く。縦の共分散をどちらか一端で縮約すると $(k_\mu k_\nu-k_\nu k_\mu)=0$。例えば $k=(0,0,0,K)$ のとき $F_{41}=iKA_1$ なので $\langle F_{41}(k)F_{41}(-k)\rangle=K^2K^{-2}=1$。一方 $A_4$ の共分散は $\xi/K^2$。前者はゲージ不変な観測量、後者はゲージ固定後の場そのものなので、異なる振る舞いになる。運動量保存の共通のデルタはここで省いた。
:::

### 一つの有限なゲージ軌道で規格化を検査する

二変数の Euclid 積分で $A=(x,y)$、ゲージ変換を $(x,y)\mapsto(x+\omega,y)$ とし、作用を $S=\kappa y^2/2$、$\kappa>0$ とする。ゲージ条件を $F[A]=a x+b y=f$、$a\ne0$ とすると $\int\dd\omega\,\delta(a(x+\omega)+by-f)|a|=1$。この一は全ての点の軌道を同じ重みで数える。デルタを $x$ について積分すれば、Jacobian の $1/|a|$ が行列式を消し、$y$ だけの Gaussian が残る。条件の $b$ を変えても $\langle y^2\rangle=1/\kappa$。一方、$x=(f-by)/a$ なので、ゲージに依存する量 $\langle x^2\rangle$ は $b$ に依存する。例えば $f=0$ なら $b^2/(a^2\kappa)$ である。$f$ を分散 $\xi$ の Gaussian で平均すると $\xi/a^2$ が加わる。これが縦の共分散の $\xi$ 依存を一自由度で表した例である。

### 電磁場を消した保存電流の非局所項

Euclid の Abel 場で $K=k^2(P_T+\xi^{-1}P_L)$、結合を $i q\int A_\mu j_\mu$ とした作用を考える。ここで $j_\mu$ は実の外部電流、または独立な Grassmann 場の偶の双一次として扱う。二次の $A$ を積分すると、規格化後の重みには
$$
\exp\left[-\frac{q^2}{2}\int j_\mu K^{-1}_{\mu\nu}j_\nu\right]
$$
が掛かる。平方完成の源は $-iq j$ なので二乗で負号が付く。積分自体は実 $A$ 上の Gaussian の Fourier 変換であり、変数を虚方向へ移す説明は有限正則化の下での解析接続として理解できる。実の保存源については $P_L j=0$ で、指数は $\xi$ に依存しない。例えば時間に依存しない $j^0(\vb x)$ だけの源では、時間体積を除くと $q^2/2$ と $\int_{\vb{k}} |\widetilde j^0|^2/|\vb k|^2$ の積が作用に現れ、二源の交差項は瞬時の Coulomb 相互作用に一致する。

物質場も積分する QED で $j$ を場の双一次にした場合、各積分配置で $j$ が保存していると仮定して縦項を最初から落とすことはできない。保存則は物質の運動方程式や Ward 恒等式を通じた相関の条件であり、オフシェルの道筋ごとに成立する条件とは違う。従ってゲージ不変な挿入の値で $\xi$ の打消しを示す場合は、ゲージ変換の変数変換や、Ward 恒等式の接触項を保つ。

## Schwinger–Dyson 方程式

有限個の積分変数の平行移動は測度を変えない。境界で被積分関数が減衰すると仮定し、$\delta/\delta\phi(x)$ の全微分を積分する。積分の挿入を $F[\phi]$ とすると
$$
\left\langle \frac{\delta F}{\delta\phi(x)}\right\rangle_J
+i\left\langle F\left(\frac{\delta S}{\delta\phi(x)}+J(x)\right)\right\rangle_J=0
\label{eq-qft-schwinger-dyson}
$$
である。これを :term[Schwinger–Dyson 方程式]{en="Schwinger–Dyson equation" yomi="しゅうぃんがーだいそんほうていしき"} という。$F=\phi(y)$ なら第一項は $\delta^4(x-y)$。自由場では $(\Box+m^2)\Delta_F=-i\delta^4$ を再現する。四次相互作用ではこの関係に $\lambda\langle\phi^3(x)\phi(y)\rangle/6$ が加わり、高次の相関へ接続する。点が重なる積は別に正則化し、古典の運動方程式を相関の中で無条件に零と置かない。

### 有限積分の方程式が示す階層

一変数で $S_E(x)=m^2x^2/2+\lambda x^4/24$、$\lambda>0$ とし、$Z_E(J)=\int_\R\dd x\,e^{-S_E+Jx}$ を考える。全微分を積分して得る微分方程式は
$$
\left(m^2\partial_J+\frac{\lambda}{6}\partial_J^3-J\right)Z_E(J)=0 .
$$
境界の指数は四次項で減衰する。$J$ の一微分後に零へ置くと $m^2\langle x^2\rangle+\lambda\langle x^4\rangle/6=1$。例えば $\lambda=0$ なら分散は $1/m^2$、$\lambda$ 一次で $\langle x^4\rangle=3/m^4+\mathcal O(\lambda)$ を代入すれば $\langle x^2\rangle=1/m^2-\lambda/(2m^6)+\cdots$ である。この点では四次相互作用が分散を減らす。ただし同じ階層の式だけで全ての解を特定するには、積分が満たす境界条件や対称性も指定する必要がある。

### 接触項をデルタとして残す

Euclid で $\delta/\delta\phi(x)$ を $e^{-S[\phi]}\prod_j\phi(y_j)$ に作用させて積分する。測度が並進不変で境界項が消える正則化では
$$
\left\langle\frac{\delta S}{\delta\phi(x)}\prod_j\phi(y_j)\right\rangle
=\sum_j\delta(x-y_j)\left\langle\prod_{i\ne j}\phi(y_i)\right\rangle .
$$
例えば自由場で二点なら $(-\partial^2+m^2)G(x-y)=\delta(x-y)$。四次相互作用なら左辺に $\lambda\langle\phi(x)^3\phi(y)\rangle/6$ が現れる。自由な先頭の値で $\langle\phi(x)^3\phi(y)\rangle_0=3G(0)G(x-y)$ だから、二点関数の逆核の補正は $\lambda G(0)/2$ となる。$G(0)$ は正則化を残したまま計算し、くりこみ条件と合わせる。

### 源の汎関数微分として階層を組み立てる

Euclid で $S_E=\int[(\partial\phi)^2/2+m^2\phi^2/2+\lambda\phi^4/4!]$ を採る。正則化された生成汎関数に対して
$$
\left[
(-\partial_x^2+m^2)\frac{\delta}{\delta J(x)}
+\frac{\lambda}{6}\frac{\delta^3}{\delta J(x)^3}-J(x)
\right]Z_E[J]=0.
$$
全微分を取った後、指数から $\phi$ を引き出す微分へ置き換えた式である。$J(y)$ で一回微分すると右辺の源の微分から $\delta(x-y)$ を得る。$J$ を二回微分した場合、$J=0$ での右辺は二つのデルタとそれぞれに対応する一点関数になる。反転対称な真空では一点関数は零なので、この方程式は奇数点の関数が零という性質と整合する。三回微分すれば通常の四点と接触位置の二点の関係になる。

:::example[質量と四次結合のくりこみ前の符号]
一変数で共分散を $D$、相互作用を $\lambda x^4/24$ とする。$x^3$ を挿入した全微分の平均から $3\langle x^2\rangle=D^{-1}\langle x^4\rangle+\lambda\langle x^6\rangle/6$。$\langle x^2\rangle=D-\lambda D^3/2+O(\lambda^2)$ と自由の $\langle x^6\rangle=15D^3$ を代入すると、$\langle x^4\rangle=3D^2-4\lambda D^4+O(\lambda^2)$。同じ組合せ論を生成汎関数の摂動展開に適用した場合と一致する。方程式は正則化した測度による関係であり、同一点の積を連続極限へ持ち込むには演算子の規格化も要る。
:::

### Grassmann の全微分と電子の逆核

Grassmann の左微分の積分は零である。実際、ある変数で一回微分するとその変数がなくなり、最高次数の係数を抽出する積分では寄与しない。$S_E=\bar\psi K\psi+S_I$ とし、$\partial_{\bar\psi_i}$ を $\bar\psi_j e^{-S_E}$ に作用させる。$\bar\psi_j$ は奇変数なので Leibniz 則の第二項に負号が付き、指数の微分にさらに負号が付く。場を順に戻すと、
$$
\left\langle\frac{\partial S_E}{\partial\bar\psi_i}\bar\psi_j\right\rangle=\delta_{ij}.
$$
自由の場合は $K_{ik}\langle\psi_k\bar\psi_j\rangle=\delta_{ij}$。順序を逆にした相関では右辺の符号も逆になる。この形は場を有限個にしたときの正確な積分恒等式である。例えば $S_E=\bar\psi (K_0+g A)\psi+S_E[A]$ とすれば $K_0\langle\psi\bar\psi\rangle$ に $g\langle A\psi\bar\psi\rangle$ が加わり、二点だけでは閉じない。位相因子や行列の向きは、奇の微分をどちらから作用させたかを固定して扱う。

### Schwinger–Dyson の階層と境界条件を同時に使う

一変数の Euclid の六次まで数えるため、$M_n=\langle x^n\rangle$、$S_E=x^2/(2D)+\lambda x^4/24$ とする。全微分を $x^n e^{-S_E}$ に作用させると、
$$
D^{-1}M_{n+1}+\frac{\lambda}{6}M_{n+3}=nM_{n-1}
\quad(n\ge1).
$$
$M_0=1$、奇数次は反転対称性で零。ただし $\lambda>0$ を固定した非摂動の漸化式では、偶数次の $M_2$ を与えないと残りが定まらない。例えば $n=1$ から $M_4=(6/\lambda)(1-M_2/D)$ なので $M_4\ge M_2^2$ は $M_2$ に上限を課すが、一つの不等式だけでは値を決めない。もとの実積分では、正測度、対称性、全モーメントの有限性が共に成り立つ。これらを切り離して形式的な常微分方程式だけを解くと、別の積分経路や正測度を持たない関数まで含める。

一方 $\lambda$ の形式級数を $\lambda=0$ の Gaussian から始める場合、次数を一つ上げるために必要な高次モーメントは前の次数で既知になる。例えば $n=1$ は $M_2=D-\lambda D M_4/6$、$n=3$ は $M_4=3D M_2-\lambda D M_6/6$。$M_4$ の一次を $M_2$ の二次へ入れれば $M_2=D-\lambda D^3/2+2\lambda^2D^5/3+\cdots$ を得る。これは摂動係数の再帰を与えるが、級数が収束することまでは示さない。実際、この積分で指数の級数を絶対値にして足すと $e^{|\lambda|x^4/24}$ となり、Gaussian の測度に対して可積分ではない。従って固定した有限の次数での Taylor 係数と、無限級数を積分の内側で和にする操作は区別する。

## Euclid 化と統計力学との対応

$t=-i\tau$ を取ると $iS=-S_E$、自由スカラーで $S_E=\int\dd^4x_E[(\partial_E\phi)^2+m^2\phi^2]/2$ である。正の質量、または零モードを除いた系で Euclid の Gaussian は収束し、二点関数は $1/(p_E^2+m^2)$ になる。時間方向の Fourier 積分の極 $\pm iE_{\vb p}$ から固定空間運動量で $e^{-E_{\vb p}|\tau|}/(2E_{\vb p})$ を得る。相互作用する場合も点を $\tau>0$ に離して完全状態を挿入すれば $|\langle0|\phi|n\rangle|^2e^{-(E_n-E_0)\tau}$ の和になるため、最も低い結合するギャップを長時間減衰から読める。

逆温度 $\beta$ では $\operatorname{Tr}e^{-\beta H}$ を時間分割する。ボース座標では行列要素の最初と最後が同じになり $\phi(\beta)=\phi(0)$。フェルミオンのコヒーレント状態のトレースでは $\operatorname{Tr}O=\int\dd\bar\theta\dd\theta\,e^{-\bar\theta\theta}\langle-\theta|O|\theta\rangle$ から反周期境界になる。この負号は一モードで $O=1$ と $O=a^\dagger a$ の係数を最上次数まで展開して確認でき、一般の積モードへ広がる。従って周波数はボースで $2\pi n/\beta$、フェルミオンで $(2n+1)\pi/\beta$。Euclid の重みは古典統計の Boltzmann 重みに対応するが、実時間の相関を取り出すには解析接続の条件が要る。

### 周期化した二点関数

自由ボース場の空間運動量を固定し、$E=\sqrt{\vb p^2+m^2}>0$ とする。温度 $1/\beta$ の周期応答は真空の $e^{-E|\tau|}/(2E)$ を整数シフトで足すと得られる。$0\le\tau\le\beta$ なら、正と負のシフトを幾何級数に分けて
$$
G_\beta(\tau)=\frac{e^{-E\tau}+e^{-E(\beta-\tau)}}{2E(1-e^{-\beta E})}.
$$
この関数は区間内で $(-\partial_\tau^2+E^2)G_\beta=0$ を満たし、接合点で導関数の跳びが $-1$ なので周期デルタの Green 関数になる。例えば $\tau=0$ では $G_\beta(0)=(1+2n_B(E))/(2E)$、ここで $n_B(E)=1/(e^{\beta E}-1)$ である。$\beta\to\infty$ で真空値 $1/(2E)$ に戻り、高温では $1/(\beta E^2)$ へ近づく。

### 有限温度でボースとフェルミの周波数が違う理由

熱的トレースで演算子を周回させると、時間順序化されたフェルミオン二点関数には二つの奇演算子を交換する負号が出る。従って Euclid 時間長 $\beta$ でボースは周期的、フェルミは反周期的である。Fourier 基底 $e^{-i\omega_n\tau}$ に課すと $e^{-i\omega_n\beta}=\pm1$ なので $\omega_n^{\rm B}=2\pi n/\beta$、$\omega_n^{\rm F}=(2n+1)\pi/\beta$。例えば温度 $T=1/\beta$ でフェルミオンの最小非零周波数は $\pi T$、ボースは零モードを持つ。そのため高温の長距離極限でボースの零モードだけが低いスケールに残る場合がある。

### 反射正値性を自由な格子極限で確認する

固定した空間運動量で $m^2+|\vb k|^2=E^2>0$ とする。自由 Euclid 核は $D(\tau)=e^{-E|\tau|}/(2E)$。正時間 $\tau_i>0$ に置いた線形挿入の係数を $z_i$ とし、一方の時刻を反射する形式を作ると
$$
\sum_{i,j}z_i^*D(-\tau_i-\tau_j)z_j
=\frac1{2E}\left|\sum_j e^{-E\tau_j}z_j\right|^2\ge0 .
$$
これは :term[反射正値性]{en="reflection positivity" yomi="はんしゃせいちせい"} の自由な二点での形である。正時間の二つの挿入が同じ一モードの一粒子状態へ写るので、この行列は階数一となる。複数の質量を正のスペクトルで足せば、各質量について同じ正の形を積分するため不等式が保たれる。一般の非 Gaussian 理論の全挿入に対する条件は二点だけでは検査できない。

:::example[二つの極の残留係数で正値性を検査する]
$D(\tau)=e^{-|\tau|}-c e^{-2|\tau|}$、$c>0$ とし、$\tau_1=1,\tau_2=2$ で $z=(e^{-2},-e^{-1})$ を取る。正の第一項の線形和は $e^{-3}-e^{-3}=0$。第二項の線形和は $e^{-4}-e^{-5}\ne0$ なので反射形式は負になる。従って運動量空間の核が正の実数になるだけでは正定値な Hilbert 空間を保証しない。例えば小さな $c$ ではこの例の Fourier 核 $2/(p^2+1)-4c/(p^2+4)$ が全実数 $p$ で正にもなり得る。
:::

:::example[温度零と時間無限遠の順序]
有限温度の自由ボース核は $0\le\tau\le\beta$ で $D_\beta(\tau)=[e^{-E\tau}+e^{-E(\beta-\tau)}]/[2E(1-e^{-E\beta})]$。$\tau$ を固定して $\beta\to\infty$ とすれば反対側から回り込む項が消え、温度零の核になる。一方 $\tau=\beta-s$ とし $s>0$ を固定した極限では回り込む項 $e^{-Es}/(2E)$ が残る。周期化の近似を有限温度で使うときは、考えている時間差が周期の両端からどれほど離れているかを調べる。
:::

### 熱的境界条件と伝播方向を比べる

一フェルミオンモードの Hamiltonian を $H=Ea^\dagger a$ とし、$E>0$ とする。熱平均の分母は $1+e^{-\beta E}$、$a(\tau)=e^{-E\tau}a$。従って $0<\tau<\beta$ の時間順序相関は $\langle a(\tau)a^\dagger(0)\rangle=e^{-E\tau}/(1+e^{-\beta E})$。$-\beta<\tau<0$ では順序が逆になり、$-\langle a^\dagger(0)a(\tau)\rangle=-e^{-E(\tau+\beta)}/(1+e^{-\beta E})$。これで反周期性を直接確かめられる。$\tau\downarrow0$ と $\tau\uparrow0$ の差は一であり、正準反交換子の大きさに等しい。従って $G$ は $(\partial_\tau+E)G=\delta$ を反周期な分布として満たす。Fourier 基底 $e^{-i\omega_n\tau}$ では係数 $1/(E-i\omega_n)$ を得る。

:::example[零温度での片側の時間発展]
上の一モードで $\beta\to\infty$ とすると、正の $\tau$ では $e^{-E\tau}$、負の $\tau$ では零になる。負の時間側の熱的な項には $e^{-\beta E}$ が含まれるためである。一方、実スカラーの固定空間運動量の時間順序相関は $e^{-E|\tau|}/(2E)$ で、両側に時間発展がある。Dirac 場には粒子・反粒子の独立なモードがあるので、一つの $a$ だけを取ったこの例を Dirac 全体の核と同一視しない。
:::

### Matsubara 和から自由エネルギーを得る

一つのボースのモード $E>0$ について、円周 $\beta$ での二次作用の固有値は $(2\pi n/\beta)^2+E^2$。全行列式は定数を含むため、$E$ による微分を先に取る。実ボースでは
$$
\frac{\partial \log Z}{\partial E}
=-E\sum_{n\in\Z}\frac1{(2\pi n/\beta)^2+E^2}
=-\beta E\,G_\beta(0).
$$
二番目では周期の Fourier 基底の規格化 $1/\beta$ を使った。既に得た $G_\beta(0)$ により $-\beta(1+2n_B(E))/2$。$E$ に関して積分すると $\log Z=-\beta E/2-\log(1-e^{-\beta E})+$ 定数である。定数は同じ振動子の Hilbert 空間で状態数を数える規格化に合わせる。従って自由エネルギー $F=E/2+\beta^{-1}\log(1-e^{-\beta E})$。第一項が真空エネルギー、第二項が熱励起である。場の箱でモードを足すと第一項は紫外正則化を要するが、$E\to\infty$ の第二項は指数で減衰する。例えば $E\beta\gg1$ で熱項は $-\beta^{-1}e^{-\beta E}$ から始まる。周期方向の行列式の比だけで、温度零との差を計算できる。

### 反射正値性を Wick の全ての多項式へ広げる

正定値の自由 Euclid Gaussian とし、時間反射 $\theta:\tau\mapsto-\tau$ が共分散を保つとする。正時間に台を持つ試験関数 $f$ で線形挿入 $\phi(f)$ を作り、$\Theta$ が時刻を反転し複素係数を共役するよう定める。正時間の線形挿入の交差共分散 $(f,h)_\theta=\langle(\Theta\phi(f))\phi(h)\rangle$ は、先の $e^{-E(\tau+\tau')}$ の式に空間運動量を積分した半正定値形式である。正時間だけの多項式を Wick の正規積の和へ展開する。これは最外次の恒等項と低次への自己縮約からなる三角な基底変換なので、有限次数の全ての多項式にできる。反射した負時間側にも同じ正規順序を使う。反射が共分散を保つので同じ側の縮約は相殺する。

従って $n$ 次と $m$ 次の Wick 積の反射形式では、縮約は正時間側と負時間側を結ぶものだけである。$n\ne m$ なら全ての因子を結べず零。$n=m$ なら
$$
\left\langle\Theta:\phi(f_1)\cdots\phi(f_n):
\,:\phi(h_1)\cdots\phi(h_n):\right\rangle
=\sum_{\pi\in S_n}\prod_i(f_i,h_{\pi(i)})_\theta .
$$
交差共分散の零ベクトルを除いて一粒子の内積空間を作ると、右辺は対称テンソル積の内積に $n!$ を掛けたもの。実際、テンソルを $n!$ 通りの置換で平均してから内積を取れば、二つの順列の相対置換だけが残る。よって種類や次数の異なる正規積の線形結合に対しても $\langle(\Theta F)F\rangle\ge0$。これで自由 Gaussian において二点の反射正値性から多項式の反射正値性が従う理由が分かる。相互作用する非 Gaussian では、同じ Wick 分解はないので二点だけの判定をこの証明で置き換えることはできない。

### 有限時間のトレースと真空極限の自由エネルギー

質量ゼロの自由な一実スカラー場を箱で考え、一様な零モードを除く。長さを無限大にして温度 $T$、体積 $V$ あたりの熱自由エネルギーを求める。各モードの真空項を引いた熱項を積分すると
$$
\frac{F_{\rm therm}}{V}=\frac{T}{2\pi^2}\int_0^\infty k^2\log(1-e^{-k/T})\dd k
=-\frac{\pi^2}{90}T^4 .
$$
最後では $y=k/T$ と置き、$-\log(1-e^{-y})=\sum_{n\ge1}e^{-ny}/n$ を使う。正項の和なので積分と交換でき、各 $\int_0^\infty y^2e^{-ny}\dd y=2/n^3$。従って積分は $-2\sum_n n^{-4}=-\pi^4/45$ である。例えば Euler のこの和は $-\pi<x<\pi$ での $x^2$ の Fourier 係数 $4(-1)^n/n^2$ と Parseval から得る。平均の二乗を引いた $x^2$ の分散は $\pi^4(1/5-1/9)$、Fourier 余弦二乗の平均が $1/2$ なので $8\sum_n n^{-4}$ に等しい。自由エネルギーから圧力は $\pi^2 T^4/90$、エネルギー密度は $F/V-T\partial_T(F/V)=\pi^2 T^4/30$。熱項の低 $k$ での対数は $k^2\log k$ となって積分可能であり、箱の一様モードの問題を除いた積分極限として得た量である。

:::example[反周期モードと熱自由エネルギーの符号]
一フェルミオンモードのトレースは $1+e^{-E/T}$。粒子・反粒子を合わせて同じエネルギーの独立なモードが $g_F$ 個あれば、質量ゼロの熱自由エネルギーは $-g_F T/(2\pi^2)\int_0^\infty k^2\log(1+e^{-k/T})\dd k$。級数の重みはボースの正項の和に $(-1)^{n-1}$ を掛けたもの。偶数項を二回引くので $1-2\cdot2^{-4}=7/8$ 倍になる。従って $F_{\rm therm}/V=-g_F(7/8)\pi^2T^4/90$。$g_F$ は実スカラー一成分に対する独立な励起の数であり、四成分 Dirac 場なら粒子の二スピンと反粒子の二スピンから $g_F=4$ である。
:::

## 演習問題

::::exercise{#exr-qft-path-integral-1 level=basic}
自由粒子の短時間核で $p$ の平方完成の中心を求めよ。

:::solution
指数を $-i\delta t[p-m(q_{j+1}-q_j)/\delta t]^2/(2m)+i m(q_{j+1}-q_j)^2/(2\delta t)$ と書けるので、中心は $m(q_{j+1}-q_j)/\delta t$ である。
:::
::::

::::exercise{#exr-qft-path-integral-2 level=basic}
$\int\dd\theta\,(a+b\theta)$ を求めよ。

:::solution
定数の積分は零、$\theta$ の積分は一なので $b$ である。
:::
::::

::::exercise{#exr-qft-path-integral-3 level=basic}
:eqref[eq-qft-free-generating-functional] から自由二点関数を求めよ。

:::solution
指数の二次項 $-\frac12 J\Delta_F J$ を二微分すると $-\Delta_F$。:eqref[eq-qft-path-generator] の $i^{-2}=-1$ を掛けて $\Delta_F$ となる。
:::
::::

::::exercise{#exr-qft-path-integral-4 level=basic}
Euclid 時間で $\beta=2$ の最小正フェルミオン周波数を求めよ。

:::solution
反周期条件は $e^{i\omega\beta}=-1$ なので、最小正値は $\pi/\beta=\pi/2$ である。
:::
::::

::::exercise{#exr-qft-path-integral-5 level=advanced}
一変数の Euclid Gaussian で $\langle\phi^4\rangle$ を計算せよ。

:::solution
$Z_E[J]=e^{J^2/(2K)}$ の $J^4$ の係数は $1/(8K^2)$。四微分して $4!/(8K^2)=3/K^2$ になる。
:::
::::

::::exercise{#exr-qft-path-integral-6 level=advanced}
:eqref[eq-qft-schwinger-dyson] を自由場、$F=\phi(y)$、$J=0$ で確認せよ。

:::solution
$\delta S/\delta\phi=-(\Box+m^2)\phi$ より $\delta^4-i(\Box+m^2)\Delta_F=0$。これは $(\Box+m^2)\Delta_F=-i\delta^4$ と同じである。
:::
::::

::::exercise{#exr-qft-path-integral-7 level=advanced}
Abel ゲージの行列式が規格化で消える理由を示せ。

:::solution
$F[A+\partial\omega]-F[A]=\Box\omega$ の微分は全ての $A$ で $\Box$。従って分子・分母の積分の外に同じ数として出せる。
:::
::::

::::exercise{#exr-qft-path-integral-8 level=advanced}
$\theta'=a\theta$、$a\ne0$ の一変数 Berezin 測度を求めよ。

:::solution
$\theta'=a\theta$ を積分して一になるには $\dd\theta'=a^{-1}\dd\theta$ が必要であり、定数の積分はどちらでも零なので全多項式で成立する。
:::
::::

::::exercise{#exr-qft-path-integral-9 level=advanced}
二つの実格子変数に対し $Z(J)=\int\dd x_1\dd x_2\,\exp(-x^\mathsf TKx/2+J^\mathsf Tx)$、$K=\begin{pmatrix}3&-1\\-1&3\end{pmatrix}$ とする。零源で $x_1^2x_2^2$ の期待値を求め、$x_1^2$ と $x_2^2$ の共分散を求めよ。

:::solution
共分散行列は $D=K^{-1}=\begin{pmatrix}3&1\\1&3\end{pmatrix}/8$。生成汎関数を $J_1$ と $J_2$ で二度ずつ微分すると、$x_1$ 対と $x_2$ 対を独立につなぐ一通りと、二組とも点をまたぐ二通りが出る。従って $\langle x_1^2x_2^2\rangle=D_{11}D_{22}+2D_{12}^2=9/64+2/64=11/64$。共分散は積の平均から平均の積 $D_{11}D_{22}$ を引いて $1/32$ である。$D_{12}$ が負でもこの共分散は同じ値になる。二次の観測量は交差縮約を二つ必要とするためである。
:::
::::

::::exercise{#exr-qft-path-integral-10 level=advanced}
Euclid 時間を三片に分けた調和振動子で、端点の固定された揺らぎの二次核を求め、自由粒子に対する揺らぎ積分の比を求めよ。時間片を $\delta\tau$ とする。

:::solution
内部時刻は二つ。共通の $m/\delta\tau$ を除くと対角は $a=2+(\omega\delta\tau)^2$、隣接成分は $-1$ である。行列式は $a^2-1$。$\omega=0$ では $a=2$ なので行列式は $3$。実の Gaussian 積分は行列式の逆平方根に比例するから、求める比は $\sqrt{3/(a^2-1)}$。例えば $\omega\delta\tau=1$ なら $a=3$ で比は $\sqrt{3/8}$。各時間片を短くしたときの連続極限では、時間片の数も $T/\delta\tau$ に増やしてから漸化式を使う。この二自由度の式で $\delta\tau\to0$ とするなら全時間も同時に零になる。
:::
::::

::::exercise{#exr-qft-path-integral-11 level=advanced}
二点の自由 Euclid 格子で $K=\begin{pmatrix}3&-1\\-1&3\end{pmatrix}$ とする。第1点に $u=2$ の二次の欠陥を加えた後の二点の共分散行列を求め、第2点の分散がどう変わるか述べよ。

:::solution
二次核は $K'=\begin{pmatrix}5&-1\\-1&3\end{pmatrix}$。逆は $D'=\begin{pmatrix}3&1\\1&5\end{pmatrix}/14$。欠陥を入れる前の第2点の分散は $3/8$。差は $5/14-3/8=-1/56$ なので、直接の二次項を変えなかった場所でも分散が下がる。これは $D'_{22}-D_{22}=-u D_{21}^2/(1+u D_{11})$ に $D_{21}=1/8$、$D_{11}=3/8$ を代入しても得られる。一方、もし元のオフ対角の共分散が零なら、この公式では第2点の分散は変わらない。
:::
::::

::::exercise{#exr-qft-path-integral-12 level=advanced}
温度 $T$、質量ゼロの自由スカラー一成分の熱自由エネルギーから、エントロピー密度 $s_{\rm th}$ とエネルギー密度 $u_{\rm th}$ の比を求めよ。零点エネルギーは除く。

:::solution
本文の $f=F_{\rm therm}/V=-\pi^2 T^4/90$ から $s_{\rm th}=-\partial_T f=2\pi^2 T^3/45$。従って $u_{\rm th}=f+Ts_{\rm th}=\pi^2 T^4/30$。比は $u_{\rm th}/(Ts_{\rm th})=3/4$、また $u_{\rm th}=3P$ である。粒子の化学ポテンシャルはここでは零なので、自由エネルギーと内部エネルギーの差は $TS$ だけである。成分数や統計によって $f$ の前の係数を替えても、同じ $T^4$ の冪を持つ限りこの比は保たれる。
:::
::::

::::exercise{#exr-qft-path-integral-13 level=advanced}
固定した空間運動量の自由 Euclid 場で、$0<\tau_1<\tau_2$ に支持を持つ二点の線形挿入から、反射内積が零になる非零の組を作れ。単粒子エネルギーを $E>0$ とする。

:::solution
反射形式は $|\sum_j e^{-E\tau_j}z_j|^2/(2E)$。従って $z_1=e^{-E\tau_2}$、$z_2=-e^{-E\tau_1}$ なら和が零である。この組は形式的な二成分ベクトルとしては非零だが、反射内積では状態を作らない。両挿入が一つの固定運動量モードへ写るため、二つの係数のうち一方向だけを測っている。例えば $\tau_2=2\tau_1$ では $z_2/z_1=-e^{E\tau_1}$。同じ時間配置で別のエネルギーのモードを足せば、そのエネルギーに対する零方向は一般に異なり、反射形式が正になることがある。
:::
::::

## 確認テスト

::::quiz
自由実スカラー場の奇数点関数はどうなるか。

- [x] 零になる
- [ ] 常に一になる
- [ ] 正の値を持つ

:::explanation
$\phi\mapsto-\phi$ が自由作用と測度を保ち、奇数個の積だけに負号を付けるため、積分は自身の負に等しい。
:::
::::

::::quiz{answer="3" tol="1e-6"}
ガウス変数 $X$ の分散が一で平均が零である。$\langle X^4\rangle$ を答えよ。

:::explanation
生成関数 $e^{J^2/2}$ の四次係数は $1/8$ なので四微分の値は $24/8=3$ である。
:::
::::

::::quiz
有限温度フェルミオンに課す境界条件を選べ。

- [ ] 周期
- [x] 反周期
- [ ] 両端で零

:::explanation
フェルミオンのコヒーレント状態でトレースを取ると端点の奇変数が負号つきで同一視される。従って時間方向に一周すると符号が変わる。
:::
::::

::::quiz{answer="2" tol="1e-6"}
$\beta=\pi$ のボース場で最小正の周波数を答えよ。

:::explanation
周期条件は $\omega\beta=2\pi n$ で、最小正値は $2\pi/\pi=2$ である。
:::
::::

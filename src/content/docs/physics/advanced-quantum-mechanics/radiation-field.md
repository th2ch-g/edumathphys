---
title: 輻射場の量子化
description: 電磁場を横波モードに分けて量子化し、原子の放出、散乱、放射補正を導く。
notation:
  - symbol: 'a_{\vb k\lambda}'
    meaning: 光子モードの消滅演算子
  - symbol: '\epsilon_0'
    meaning: 真空の誘電率
  - symbol: '\vb d'
    meaning: 原子の電気双極子
  - symbol: 'A_{eg}'
    meaning: 自然放出率
---

## この章の目標

- 横波モードを正準量子化し、光子のエネルギーと偏光を求める。
- 数状態とコヒーレント状態の光検出統計を比較する。
- 真空場の寄与から Casimir 力と自然放出を導く。
- 双極子選択則を使い、散乱と Lamb シフトの起源を説明する。

## 自由電磁場の正準量子化と光子

真空に電荷と電流がないとき、Coulomb ゲージ $\divg\vb A=0$ と $\phi=0$ を取る。境界は体積 $\mathcal V$ の周期箱とし、周期箱の $\vb k=0$ の自由度を除き、$\omega_k>0$ の輻射モードだけを扱う。$\vb E=-\dot{\vb A}$、$\vb B=\rot\vb A$ で、Lagrange 関数は $\epsilon_0\int(\dot{\vb A}^2-c^2(\rot\vb A)^2)\dd^3r/2$ である。共役運動量は $\epsilon_0\dot{\vb A}=-\epsilon_0\vb E$。実の正規直交横波基底に展開すると、各座標 $q$ は $\epsilon_0(\dot q^2-\omega^2 q^2)/2$ の調和振動子になる。平面進行波へ戻して

$$
\vb A(\vb r)=\sum_{\vb k,\lambda}
\sqrt{\frac{\hbar}{2\epsilon_0\mathcal V\omega_k}}
(\vb e_{\vb k\lambda}\,a_{\vb k\lambda}e^{i\vb k\cdot\vb r}
 +\vb e_{\vb k\lambda}^*\,a_{\vb k\lambda}^\dagger e^{-i\vb k\cdot\vb r}),
\quad \omega_k=c|\vb k|.
\label{eq-aqm-vector-potential}
$$

$\vb e_{\vb k\lambda}$ は $\vb k$ と直交する二つの単位偏光である。実定在波の二モードから $\pm\vb k$ の二進行波へユニタリに移るので、二重に数えていない。:eqref[eq-aqm-boson-actions] のボソンの代数を各モードに課す。実モードごとに共役な $q,p_q=\epsilon_0\dot q$ を取り、$a=\sqrt{\epsilon_0\omega/(2\hbar)}\,q+i p_q/\sqrt{2\epsilon_0\hbar\omega}$ とすれば、$[q,p_q]=i\hbar$ は $[a,a^\dagger]=1$ になる。逆に $q,p_q$ を $a,a^\dagger$ で表して $p_q^2/(2\epsilon_0)+\epsilon_0\omega^2 q^2/2$ に代入すると、$aa$ と $a^\dagger a^\dagger$ が打ち消し合う。進行波への基底変換は同じ $\omega$ の中だけなので、
$H_{\mathrm{rad}}=\sum_{\vb k\lambda}\hbar\omega_k(a^\dagger a+1/2)$ になる。$\vb E$ の消滅項は $i\sqrt{\hbar\omega_k/(2\epsilon_0\mathcal V)}\,\vb e_{\vb k\lambda}a e^{i\vb k\cdot\vb r}$ である。随伴の生成項には逆の $i$ の符号が付く。

$\sum_\lambda e_i e_j^*=\delta_{ij}-\hat k_i\hat k_j$ は横波平面への射影である。従って $\comm{A_i(\vb r)}{\Pi_j(\vb r')}=i\hbar\,\delta^\perp_{ij}(\vb r-\vb r')$ で、$\delta^\perp$ はこの射影の Fourier 和である。縦成分を独立な振動子に加えると Gauss の制約に反する。:term[光子]{en="photon" yomi="こうし"} はモードの一量子で、エネルギー $\hbar\omega$、運動量 $\hbar\vb k$ を持つ。運動量は $\vb r$ の並進でモード振幅が $e^{i\vb k\cdot\vb r}$ の位相を得ることから決まる。

:::example[可視光の一量子]
波長 $500\,\mathrm{nm}$、$hc\simeq1240\,\mathrm{eV\,nm}$ とすれば $E=hc/\lambda\simeq2.48\,\mathrm{eV}$ である。角周波数は $2\pi c/\lambda\simeq3.77\times10^{15}\,\mathrm{s^{-1}}$ となる。
:::

:::example[円偏光]
$\vb k$ を $z$ 方向に取り、$\vb e_\pm=(\hat{\vb x}\pm i\hat{\vb y})/\sqrt2$ を選ぶ。能動 $z$ 回転で $\vb e_+\mapsto e^{-i\theta}\vb e_+$、$\vb e_-\mapsto e^{i\theta}\vb e_-$ なので、この向きの角運動量の固有値は $\pm\hbar$ である。$\vb e_z$ は横波条件を満たさない。
:::

### 一光子の規格化と検出確率の密度

$\ket{1_f}=\sum_{\vb k\lambda}f_{\vb k\lambda}a^\dagger_{\vb k\lambda}\ket0$ とする。真空の正準交換関係により $\norm{1_f}^2=\sum|f|^2$ である。正規化すれば平均エネルギーは零点を除き $\sum\hbar\omega_k|f|^2$ であり、エネルギー分散は同じ重みでの $\hbar\omega_k$ の分散である。検出器が局所の固定偏光 $\vb d$ に結合すると、一光子を吸収する振幅は $\bra0\vb d\cdot\vb E^{(+)}(\vb r,t)\ket{1_f}$ に比例する。これは $f$ を単純に逆 Fourier 変換したものに加えて、モードごとの $\sqrt{\omega_k}$ と横偏光の射影を含む。異なる光子数の成分と同様、モードの重ね合わせでは吸収した先の場が同じ真空なので振幅を足す。

:::example[二つの進行波から定在する検出分布を作る]
同じ周波数と共通の横偏光の $\vb k,-\vb k$ を取り、$\ket{1_f}=(\ket{1_{\vb k}}+e^{i\varphi}\ket{1_{-\vb k}})/\sqrt2$ とする。$\vb k$ を $z$ 方向に選べば吸収振幅の位置の因子は $(e^{ikz}+e^{i\varphi}e^{-ikz})/\sqrt2$ である。強度は $1+\cos(2kz-\varphi)$ となる。$\varphi=\pi$ では原点で零、$z=\lambda/4$ で最大 2 になる。運動量を測ると二結果は位相に依存せず等確率である。
:::

### 横射影を場の等時刻交換子で確かめる

正周波数の部分を $\vb A^{(+)}$、その随伴を $\vb A^{(-)}$ とし、$\boldsymbol\Pi=\epsilon_0\dot{\vb A}$ を使う。正周波数因子の微分は $-i\omega$、負周波数は $+i\omega$ である。混合する二つの交換子はそれぞれ $+i\hbar/(2\mathcal V)$ に横偏光の成分と位相を掛けたものになる。$\vb k$ の和で一方を $-\vb k$ に替えると、偏光射影 $\delta_{ij}-\hat k_i\hat k_j$ が変わらず位相も揃うので二つが足される。同種の周波数の交換子は零である。よって
$$
[A_i(\vb r),\Pi_j(\vb r')]=\frac{i\hbar}{\mathcal V}
\sum_{\vb k\ne0}(\delta_{ij}-\hat k_i\hat k_j)e^{i\vb k\cdot(\vb r-\vb r')}
$$
となる。周期箱の零波数を除いた横試験関数に積分すれば、射影が各 Fourier 成分を変えないので通常の $i\hbar$ 倍の恒等作用になる。縦試験関数には射影が零であり、残りの縦部分を物理振動子の交換関係に加えない。

:::example[斜め方向の横偏光を足す]
$\hat{\vb k}=(1,1,0)/\sqrt2$ なら偏光基底を $(1,-1,0)/\sqrt2$ と $(0,0,1)$ に選べる。外積を足した行列は $\begin{pmatrix}1/2&-1/2&0\\-1/2&1/2&0\\0&0&1\end{pmatrix}$ である。試験ベクトルを $\vb e_x$ に取ると横成分のノルム二乗は $1/2$、$\vb e_z$ なら 1 になる。直線偏光の吸収率が軸によって変わる計算でも同じ射影が現れる。
:::

### 横波の運動量を場の積から求める

周期境界で古典の場の運動量は $\vb P=\epsilon_0\int\vb E\times\vb B\,\dd^3x$ である。成分で $(\vb E\times(\grad\times\vb A))_i=\sum_j E_j\partial_i A_j-E_j\partial_j A_i$ と書く。第二項を部分積分すると、境界と $\grad\cdot\vb E$ が零で消える。従って $P_i=-\int \sum_j\Pi_j\partial_i A_j\dd^3x$ である。量子では Hermite に整えるため積を対称化する。:eqref[eq-aqm-vector-potential] を代入すると、座標積分が $a_{\vb k}^\dagger a_{\vb k'}$ に対し $\vb k=\vb k'$、二つの生成または二つの消滅には $\vb k'=-\vb k$ を課す。後者は波数の符号が逆の項が足されて消える。前者と逆順の積は各モードで $\hbar k_i(a^\dagger a+aa^\dagger)/2$ になる。真空定数 $\sum_{\vb k\lambda}\hbar\vb k/2$ は $\pm\vb k$ を対にする切断で零なので、場の運動量は $\sum_{\vb k\lambda}\hbar\vb k N_{\vb k\lambda}$ である。エネルギーの零点項と異なり、運動量では反対方向が相殺する。

:::example[同じ周波数の二方向を重ねたときの運動量分散]
$\vb k,-\vb k$ の一光子状態を確率重み $3/4,1/4$ の純粋な重ね合わせにする。$k$ の方向への成分 $P_\parallel$ は $\pm\hbar k$ なので平均は $\hbar k/2$、分散は $\hbar^2k^2-(\hbar k/2)^2=3\hbar^2k^2/4$ である。一方エネルギーは $\hbar c k$ に固定され、分散は零である。二つの向きの相対位相は運動量とエネルギーの分布には出ないが、局所の電場を測る吸収には干渉として出る。
:::

## 光子の状態（Fock 状態・コヒーレント状態）

一モードの数状態 $\ket n$ では $\langle a\rangle=0$ だが $\langle a^\dagger a\rangle=n$ である。従って平均電場は零でもエネルギーと光検出率は有限である。一方 :ref[def-aqm-coherent]{full} の状態なら $\langle a\rangle=\alpha$ なので、電場の平均は複素振幅 $\alpha$ を持つ古典モードになる。古典振動の位相は数基底の異なる成分の相対位相にある。

$X=(a+a^\dagger)/\sqrt2$、$P=(a-a^\dagger)/(i\sqrt2)$ とすると $\comm{X}{P}=iI$ である。数状態では $\langle X^2\rangle=\langle P^2\rangle=n+1/2$ となる。各二次式で $a^2,a^{\dagger2}$ の期待値が零だからである。コヒーレント状態では $b=a-\alpha$ が状態を消し、$\comm{b}{b^\dagger}=I$ なので $X,P$ の分散はいずれも $1/2$ である。

光子を短時間のうちに二個数える場合には、最初の吸収が $a$ を、次がさらに $a$ を掛ける。二個同時の理想化した規格化相関は $g^{(2)}(0)=\langle a^{\dagger2}a^2\rangle/\langle a^\dagger a\rangle^2$ である。これが :term[二次コヒーレンス]{en="second-order coherence" yomi="にじこひーれんす"} である。平均数が正の場合、数状態なら $1-1/n$、コヒーレント状態なら $1$ となる。真空では平均数が零なのでこの規格化相関は定義しない。熱状態で $p_n=(1-z)z^n$、$0<z<1$ とすると母関数は $(1-z)/(1-zt)$ で、二階微分は $2[z/(1-z)]^2$ となるため $g^{(2)}(0)=2$ である。

:::example[単一光子]
$\ket1$ に対し一回の吸収は真空を作り、二回目は零である。従って $g^{(2)}(0)=0$ である。理想的な二分岐の検出器では同時クリックがない。
:::

:::example[二光子数状態と同じ強度のコヒーレント光]
$\ket2$ と $\ket{\alpha=\sqrt2}$ はともに平均数 2 である。前者の二光子階乗モーメントは 2、後者は 4 なので、$g^{(2)}(0)$ はそれぞれ $1/2,1$ である。
:::

### 偏光測定とモードの部分情報

同じ波数の二直交偏光を $H,V$ とする。一光子空間で $a_H^\dagger\ket0,a_V^\dagger\ket0$ を基底に選ぶ。角 $\theta$ の偏光子は透過側の生成演算子を $\cos\theta\,a_H^\dagger+\sin\theta\,a_V^\dagger$ と定める。初期係数 $f_H,f_V$ なら透過確率は $|\cos\theta f_H+\sin\theta f_V|^2$ である。二直交出力を足すと規格化の 1 に戻る。場のコヒーレント状態について、線形なモード回転を用いると $a_H\ket{\alpha_H,\alpha_V}=\alpha_H\ket{\alpha_H,\alpha_V}$ 等から、出力の消滅演算子にも固有値を持つ。従って出力も各モードのコヒーレント状態の直積になり、透過側の平均数は同じ振幅を二乗した値になる。

:::example[損失で二次相関が変わるか]
入力 $\ket n$、他方の入力真空、透過率 $\eta$ の理想ビーム分岐を考える。出力で $(\sqrt\eta\,a^\dagger+\sqrt{1-\eta}\,b^\dagger)^n/\sqrt{n!}$ を展開すると、透過数 $m$ は $\binom nm\eta^m(1-\eta)^{n-m}$ である。平均は $n\eta$、二次階乗モーメントは $n(n-1)\eta^2$ で、$\eta>0$ かつ $n\ge1$ なら $g^{(2)}(0)=1-1/n$ のままである。$n=2,\eta=1/2$ の透過数の確率は $(1/4,1/2,1/4)$ となる。二個の同時クリック数は減るが、平均数の二乗で割ると損失率が消える。
:::

### 光子数と電場の周波数成分を区別する

一モードで $\vb E^{(+)}(t)=i\mathcal E_0\,\vb e\,a e^{-i\omega t}$ とし、局所偏光の係数を除く。この自由な単一モードでは、状態の定常性を仮定せずに一次相関 $G^{(1)}(t,t+\tau)=\langle E^{(-)}(t)E^{(+)}(t+\tau)\rangle$ は $\mathcal E_0^2\langle N\rangle e^{-i\omega\tau}$ である。同じ平均を持つ数、熱、コヒーレントの各準備では、この一モード・無損失の一次相関だけで区別できない。二つの異なる周波数の場が独立で、$\langle a_1^\dagger a_2\rangle=0$ なら $G^{(1)}$ は二つの指数の和になる。相関の絶対値が下がるのはこの例ではモードごとの位相差であり、各モード自体が減衰するという仮定はしていない。

:::example[二周波数が等強度で一次干渉が消える時刻]
$G^{(1)}$ を $\tau=0$ で規格化し、二つの寄与を各 $1/2$ とする。絶対値は $|\cos((\omega_1-\omega_2)\tau/2)|$ である。初めの零は $\pi/|\omega_1-\omega_2|$、次の完全な回復は $2\pi/|\omega_1-\omega_2|$ である。理想的な二つの狭いモードでは長時間で単調に消えるとは限らない。連続な周波数帯なら包絡は帯域の Fourier 変換になり、線形スペクトルと干渉のデータを相互に検算できる。
:::

### 偏光の Stokes 成分で位相を復元する

二つの偏光モードの双一次 $S_0=N_H+N_V$、$S_1=N_H-N_V$、$S_2=a_H^\dagger a_V+a_V^\dagger a_H$、$S_3=-i(a_H^\dagger a_V-a_V^\dagger a_H)$ を取る。一光子空間ではそれぞれ $I,Z,X,Y$ になる。例えば $a_H^\dagger a_V$ は $V$ を $H$ に移す一つの行列単位である。従って偏光の一次密度は $\rho_{\rm pol}=(I+\sum_{j=1}^3 \langle S_j\rangle\sigma_j^{\rm pol})/2$ であり、$\sigma_j^{\rm pol}$ はこの順序で $Z,X,Y$ を表す。二つの直線偏光だけを測ると $S_1$ しか分からない。$H\pm V$ の基底で $S_2$、$H\pm iV$ の基底で $S_3$ を得る。周波数や進行方向を別に持つ一光子状態で、各モードに共通の偏光ラベルを指定した場合は、モードごとの $S_j$ を足せば、それらを跡消去した偏光密度が得られる。進行方向が異なる横波には、その方向ごとの二つの横偏光の基底を指定する必要がある。

:::example[同じ二直線偏光の度数で位相が逆になる]
$\ket f=(\sqrt3\ket H+i\ket V)/2$ は $(\langle S_1\rangle,\langle S_2\rangle,\langle S_3\rangle)=(1/2,0,\sqrt3/2)$ である。成分を共役すると $S_3$ だけが逆符号になる。$H,V$ の測定では両者とも $3/4,1/4$、$H\pm V$ ではともに $1/2,1/2$ である。$H+iV$ の方向に射影する確率は最初の状態で $(1+\sqrt3/2)/2$、共役で $(1-\sqrt3/2)/2$ となる。平均電場は一光子の各状態で零でも、偏光の相対位相はこの一次の偏光測定で得られる。
:::

## 真空のゆらぎと Casimir 効果

各モードの真空は消滅演算子に消されるが、電場の二乗には $aa^\dagger$ が寄与する。例えば実の偏光で一モードの偏光方向成分の分散は $\hbar\omega/(2\epsilon_0\mathcal V)$ である。モード全体を一点で足した量は発散し得るので、観測する場の平均範囲と境界による差を指定する。

面積が十分大きい二理想導体平板を距離 $d$ に置く。接線電場の零条件により垂直波数は $n\pi/d$ となる。$n\ge1$ で二偏光、$n=0$ で一偏光である。平行波数の大きさを $k$ と置き、滑らかな高波数 cutoff $C_\Lambda(u)$ を使う。零点エネルギーの面積密度は

$$
\sum_{n=0}^\infty{}' F_\Lambda(n),\qquad
F_\Lambda(x)=\frac{\hbar c}{2\pi}
   \int_{\pi x/d}^{\infty} u^2 C_\Lambda(u)\dd u .
$$

プライムは最初の半重みで、$u=\sqrt{k^2+(\pi x/d)^2}$ を用いた。平板間の自由空間の寄与は $\int_0^\infty F_\Lambda(x)\dd x$ である。板を動かすと外側の自由空間の体積も変わるので、この体積に比例する項は力の計算で相殺する。板ごとの位置に依存しない自己エネルギーを除き、板間隔に依存する相互作用をこの差から求める。Euler–Maclaurin の半重み付き和から積分を引くと $-B_2 F_\Lambda'(0)/2!-B_4 F_\Lambda'''(0)/4!+\cdots$ になる。符号は下端での導関数なので負である。$C_\Lambda(0)=1$、$F_\Lambda'(0)=0$、$F_\Lambda'''(0)=-(\hbar c/\pi)(\pi/d)^3$ である。例えば $C_\Lambda=e^{-u/\Lambda}$ なら四次より高い導関数と Euler–Maclaurin の剰余は $\Lambda$ の負冪で消える。剰余は六階導関数の全積分で評価でき、$F_\Lambda(x)=\Lambda^3 F_1(x/\Lambda)$ と尺度を替えると $\Lambda^{-2}$ になる。$B_4=-1/30$ を入れて

$$
\frac{E_{\mathrm C}}{\mathcal A}=-\frac{\pi^2\hbar c}{720d^3},\qquad
P_{\mathrm C}=-\frac{\pi^2\hbar c}{240d^4}.
\label{eq-aqm-casimir}
$$

ここで圧力は面積密度を $d$ で微分して負号を取ったものである。:term[Casimir 効果]{en="Casimir effect" yomi="かしみーるこうか"} のこの式は零温度での理想平板の結果である。

:::example[距離依存性]
間隔を $d$ から $d/2$ に縮めると引力の圧力の大きさは $16$ 倍になる。$\hbar c\simeq3.16\times10^{-26}\,\mathrm{J\,m}$ を用いれば $d=1\,\mathrm{\mu m}$ で約 $1.30\times10^{-3}\,\mathrm{Pa}$ である。
:::

:::example[真空電場のモード体積]
同じ周波数で量子化体積を八分の一にする理想モードでは、真空電場の標準偏差が $\sqrt8$ 倍になる。共振器のモードを使う場合には、一様な箱の代わりに規格化したモード関数の原子位置での値を使う。
:::

### 指数 cutoff の和を展開する

$x=\pi/(d\Lambda)$ とする。上の $F_\Lambda(n)$ の積分を部分積分すると
$$
F_\Lambda(n)=\frac{\hbar c\Lambda^3}{2\pi}e^{-nx}[(nx)^2+2nx+2]
$$
となる。$g(x)=\sum_{n=0}^\infty{}'e^{-nx}=1/2+(e^x-1)^{-1}$ とすると、和の無次元部は $x^2g''-2xg'+2g$ である。一方、積分では $g$ の代わりに $1/x$ を使える。$g(x)=x^{-1}+x/12-x^3/720+\order{x^5}$ なので、線形項の係数は $0-2+2=0$、三次項は $6-6+2=2$ となる。従って差の最初は $-x^3/360$ で、$\Lambda^3x^3=(\pi/d)^3$ を掛けて :eqref[eq-aqm-casimir] の有限項が得られる。これで、各モードを先に無限大へ送った二つの発散和を差し引く操作を避けられる。

:::example[面積から力を求める]
$d=1\,\mathrm{\mu m}$ で十分広い平板の対向面積を $10^{-6}\,\mathrm{m^2}$ とすれば、理想模型の力の絶対値は約 $1.30\times10^{-9}\,\mathrm N$ である。$d=0.5\,\mathrm{\mu m}$ なら約 $2.08\times10^{-8}\,\mathrm N$ になる。同じ面積でも反射率、温度、端の寄与が無視できない条件では、この理想間隔依存からずれる。
:::

### 二平板のモード数を境界条件で追う

平行波数を $x$ 方向に取り、$z=0,d$ が導体面とする。横電気型では電場を $y$ だけに取り、$E_y$ は二面で零なので $\sin(n\pi z/d)$、$n\ge1$ である。横磁気型を $B_y$ の形で書けば、自由 Maxwell 方程式の $\grad\times\vb B$ の $x$ 成分から接線電場 $E_x$ は $\partial_zB_y$ に比例する。従って $B_y$ は余弦となり $n=0$ も許される。この零成分では電場は法線方向で、$k_\parallel\ne0$ なら伝播する一偏光になる。$n=0,k_\parallel=0$ の静的な成分を同じ振動子の一種として足す必要はない。こうして $n\ge1$ に二偏光、零に一偏光という半重みが生じる。横面積を連続極限にするときだけ $k_\parallel=0$ の単一点の測度は零になる。

:::example[最初の垂直励起の周波数]
$d=1\,\mathrm{\mu m}$ の理想平板で $n=1$、$k_\parallel=0$ の角周波数は $\pi c/d\simeq9.42\times10^{14}\,\mathrm{s^{-1}}$ であり、エネルギーは約 $0.620\,\mathrm{eV}$ である。$n=2$ の同じ横波数では二倍になる。一方 $n=0$ のモードは $\omega=ck_\parallel$ で、垂直量子化によるこのギャップを持たない。全モードの最低周波数を一律に $\pi c/d$ と置くと、差分の和の半重みも失う。
:::

### 高温の理想平板で零周波数が残る

$\beta=1/(k_BT)$ とし、各モードの自由エネルギーを $k_BT\log[2\sinh(\beta\hbar\omega/2)]$ とする。積の恒等式 $\sinh x/x=\prod_{n=1}^\infty(1+x^2/(\pi n)^2)$ を用いれば、境界差の自由エネルギーは虚周波数 $\xi_n=2\pi n k_BT/\hbar$ ごとの対数行列式へ分けられる。この恒等式は有限分割した $-\partial_\tau^2+\omega^2$ の周期固有値の積からも得られ、自由な分母で割って零周波数を分けたものである。ここではまず各偏光の固有値に $\xi_n^2$ を足した一方向の作用素を考える。固定横波数 $k$ の Dirichlet の境界行列式は $K=\sqrt{k^2+\xi_n^2/c^2}$ として $\sinh(Kd)/K$ に比例する。理由は $y''=K^2 y$、$y(0)=0,y'(0)=1$ の右端値による。格子で $-\partial_z^2+K^2$ を作ると主小行列式が $D_m=(2+\varepsilon^2 K^2)D_{m-1}-D_{m-2}$ を満たし、同じ漸化式が零エネルギーの境界解にも成り立つ。Neumann では同じ正の垂直波数に $K^2$ の一つが加わる。これは $d$ に依存しない。$d$ に線形な自由空間項と独立な表面の定数を引くと、残る対数は各偏光で $\log(1-e^{-2Kd})$ である。零周波数の $\log(\beta\hbar\omega)$ は $\log(\omega^2)$ の半分を含み、積の正周波数項はそれぞれ $\log(\omega^2+\xi_n^2)$ を一つ含む。従って二偏光を足した有限温度の相互作用は
$$
\frac{F_{\rm int}}{\mathcal A}
=2k_BT\sum_{n=0}^{\infty}{}'\int\frac{\dd^2 k}{(2\pi)^2}
\log(1-e^{-2d\sqrt{k^2+\xi_n^2/c^2}}).
$$
温度を零に近づけると $k_BT\sum_n{}'$ は $(\hbar/2\pi)\int_0^\infty\dd\xi$ になる。$\xi/c$ と二つの平行波数で半空間の球座標を取ると、自由エネルギーの面積密度は $(\hbar c/(2\pi^2))\int_0^\infty K^2\log(1-e^{-2dK})\dd K$ になる。対数を展開し $\int_0^\infty K^2e^{-2mdK}\dd K=1/(4m^3d^3)$ と $\zeta(4)=\pi^4/90$ を使うと、零温度のエネルギーに一致する。

$k_BT\gg\hbar c/(2\pi d)$ では非零の $\xi_n$ が $e^{-2\xi_n d/c}$ を含むので、理想完全反射の模型の先頭は二偏光の零項である。零項は Gaussian の自由エネルギーの $1/2$ を持つので、合計は
$$
\frac{F_{\rm int}}{\mathcal A}\simeq
k_BT\int\frac{\dd^2k}{(2\pi)^2}\log(1-e^{-2kd})
=-\frac{k_BT\,\zeta(3)}{8\pi d^2}.
$$
最後は $-\log(1-e^{-2kd})=\sum_{m\ge1}e^{-2mkd}/m$ と $\int_0^\infty k e^{-2mkd}\dd k=1/(4m^2d^2)$ による。温度一定で微分すると圧力は $-k_BT\zeta(3)/(4\pi d^3)$ である。これは理想な境界が零周波数でも二偏光に反射する場合であり、実金属の低周波応答をこの仮定だけから決めない。

:::example[零温度と高温で間隔の冪を比べる]
高温の式で $d$ を二倍にすると力の面積密度の絶対値は $1/8$ になる。零温度の式では $1/16$ である。二つの極限を同じ条件で同時に適用はできないので、間隔と温度を指定して $\xi_1d/c$ の大きさを見る。例えば $k_BT d/(\hbar c)=1$ の理想模型では $\xi_1d/c=2\pi$ で、非零の周波数の指数に $e^{-4\pi}$ 以下の因子が出る。
:::

## 原子と輻射場の相互作用

固定原子核の近くの荷電粒子に対し、SI 単位の最小結合は $(\vb p-Q\vb A)^2/(2m)+Q\phi$ である。線形項は $-Q(\vb p\cdot\vb A+\vb A\cdot\vb p)/(2m)$ で、Coulomb ゲージでは二項の差 $\comm{p_i}{A_i}$ の和が零になる。二次項 $Q^2 A^2/(2m)$ も同じ規約から出る。原子の大きさ $a$ に比べて波長が長い $ka\ll1$ なら、$\vb A$ の内部座標依存を先頭で落とせる。電気双極子 $\vb d=\sum_i Q_i\vb r_i$ は単位 $\mathrm{C\,m}$ を持つ。

古典場の双極子近似では、波動関数を $e^{-i\vb d\cdot\vb A(t)/\hbar}$ で変換すると、この共役作用で $\vb p_i\mapsto \vb p_i+Q_i\vb A$ となる。時間微分から $+\vb d\cdot\dot{\vb A}=-\vb d\cdot\vb E$ が加わる。量子場での線形結合も $\comm{H_{\mathrm{at}}}{\vb r}=-i\hbar\vb p/m$ により共鳴の行列要素で等価である。例えば放出では $\mel g{\vb p/m}e=-i\omega_0\mel g{\vb r}e$ であり、生成側の電場が $-i\omega \vb A^{(-)}$ なので $\omega=\omega_0$ で一致する。量子場の完全な表示変換では二次の自己エネルギーなども変換するので、線形項だけを任意の非共鳴計算で置き換えない。

:::example[一光子放出行列要素]
$H_I=-\vb d\cdot\vb E$ でモードに最初 $n$ 光子があると、放出行列要素の絶対値は $\sqrt{\hbar\omega/(2\epsilon_0\mathcal V)}\,|\vb d_{ge}\cdot\vb e^*|\sqrt{n+1}$ である。$n=0$ でも 1 の項が残る。
:::

:::example[回転波近似]
二準位原子で $H_I$ の相互作用表示には $a\sigma_+$ の位相 $e^{i(\omega_0-\omega)t}$ と $a^\dagger\sigma_+$ の位相 $e^{i(\omega_0+\omega)t}$ が現れる。結合が弱く離調が小さいとき、後者の時間積分は高速振動で抑えられる。共鳴の一モードで $\hbar g(a\sigma_++a^\dagger\sigma_-)$ を残すと $\ket{e,n},\ket{g,n+1}$ 間の非対角係数は $\hbar g\sqrt{n+1}$ で、遷移確率は $\sin^2(g\sqrt{n+1}t)$ になる。
:::

### 離調した単一モード

回転波近似で得る :term[Jaynes–Cummings 模型]{en="Jaynes–Cummings model" yomi="じぇいんずかみんぐすもけい"} の閉じた対 $\ket{e,n},\ket{g,n+1}$ で平均エネルギーを除くと Hamiltonian は $\hbar(\delta\sigma_z/2+g\sqrt{n+1}\sigma_x)$、$\delta=\omega_0-\omega$ である。二乗が $\hbar^2(\delta^2/4+g^2(n+1))I$ となるから、$\Omega_n^2=\delta^2+4g^2(n+1)$ として
$$
p_{e\to g}(t)=\frac{4g^2(n+1)}{\Omega_n^2}\sin^2(\Omega_n t/2).
$$
この対は励起数 $a^\dagger a+\ketbra e e$ が保存されるため閉じている。捨てた項はこの数を二つ変え、高速位相を持つ。その抑制には $g\sqrt{n+1}\ll\omega+\omega_0$ が必要で、離調や近傍の他準位への結合も考慮する。

:::example[離調で振幅が抑えられる場合]
$n=0,\delta=2g$ にすると $\Omega_0=2\sqrt2 g$ で確率の上限は $1/2$ である。離調が零なら上限は 1 になる。固定 $\delta=2g$ のまま $n=3$ に増やすと上限は $4/5$ になる。誘導放出を含む $\sqrt{n+1}$ の結合が離調に対して大きくなるためである。
:::

### 二次結合と和則を静的なゲージで照合する

一荷電粒子に、全空間で一定で時間に依存しない $\vb A=A\hat{\vb x}$ を入れる。全空間の束縛問題で、この $\vb A$ は $e^{iqAx/\hbar}$ により運動量を平行移動しただけなので束縛の固有値は不変である。摂動で見ると $A^2$ の一次補正 $q^2A^2/(2m)$ に、線形項 $-qAp_x/m$ の二次補正
$$
-\frac{q^2A^2}{m^2}
\sum_{n\ne i}\frac{|\mel n{p_x}i|^2}{E_n-E_i}
$$
が加わる。ここでは孤立した非縮退準位を取り、連続状態があれば和に積分を含める。交換子 $[H,x]=-i\hbar p_x/m$ によって分子を $m^2(E_n-E_i)^2|\mel n x i|^2/\hbar^2$ に直す。完全性で $\sum_n(E_n-E_i)|\mel n x i|^2=\frac12\mel i{[x,[H,x]]}i=\hbar^2/(2m)$ である。二次補正は $-q^2A^2/(2m)$ となり、二つは相殺する。有限な数の中間状態への切断ではこの和則を満たすか確かめる必要がある。単連結という条件を外すと一定に見える接線の $\vb A$ に磁束情報が残る場合がある。

:::example[振動子で和則が上下二準位だけに閉じる]
$H=p_x^2/(2m)+m\Omega^2x^2/2$ の $\ket n$ では $x=\sqrt{\hbar/(2m\Omega)}(b+b^\dagger)$ なので、上への重み付き項は $(n+1)\hbar^2/(2m)$、下への項は $-n\hbar^2/(2m)$ である。$n=2$ ではそれぞれ係数 $3,-2$ であり、合計で 1 になる。二準位に切って上側のみを残すと係数 3 となる。線形の共鳴過程に適した準位の切断も、静的な二次補正にそのまま適切とは限らない。
:::

### 遠い共鳴の二次シフト

回転波近似の一対で $H_{\rm pair}=\hbar\begin{pmatrix}\delta/2&g\sqrt{n+1}\\g\sqrt{n+1}&-\delta/2\end{pmatrix}$ とする。$\delta>0$、$g\sqrt{n+1}\ll\delta$ なら固有値の平方根を展開して
$$
E_\pm=\pm\hbar\left[\delta/2+\frac{g^2(n+1)}{\delta}
-\frac{g^4(n+1)^2}{\delta^3}+\cdots\right].
$$
正の枝はほぼ $\ket{e,n}$、負の枝はほぼ $\ket{g,n+1}$ である。固有ベクトルの反対成分の先頭は $g\sqrt{n+1}/\delta$ だから、遷移が強く抑えられても位相の二次シフトは残る。実際、時間発展に対して先頭のエネルギー差を使えば、観測時刻が長いと小さなシフトでも位相が大きくなる。非共鳴の上向きと下向きでは分母の符号が逆になるので、二枝のシフトは逆になる。

:::example[光子数を変えたときの分散的なシフト]
$\delta=20g$、$n=3$ なら混合係数の先頭は $1/10$、上の準位の二次シフトは $\hbar g/5$ である。共鳴からの距離に対しては $1/100$ の割合だが、時間 $5/g$ では位相の変化は約 1 rad になる。$n=0$ でもシフトが $\hbar g/20$ 残るのは、上の状態から $\ket{g,1}$ へ仮想遷移できるためである。この模型は一つの共振器モードに対するもので、自由空間の全モードによる Lamb シフトをそのまま表す数ではない。
:::

## 自然放出と Einstein 係数

原子の二準位間隔を $\hbar\omega_0$ とし、真空を初期場に取る。:ref[law-qm-golden] の率 $2\pi|\mel f{H_I}i|^2\delta(E_f-E_i)/\hbar$ を連続モードに足す。エネルギーのデルタ関数は $\delta(\omega-\omega_0)/\hbar$、和は $\mathcal V\sum_\lambda\int k^2\dd k\,\dd\Omega/(2\pi)^3$ に替える。偏光和が $|\vb d|^2-|\hat{\vb k}\cdot\vb d|^2$ で、角平均は $2|\vb d|^2/3$ である。これは $\int\hat k_i\hat k_j\dd\Omega=(4\pi/3)\delta_{ij}$ から分かる。従って

$$
A_{eg}=\frac{\omega_0^3|\vb d_{ge}|^2}{3\pi\epsilon_0\hbar c^3}.
\label{eq-aqm-spontaneous-rate}
$$

これは :term[自然放出]{en="spontaneous emission" yomi="しぜんほうしゅつ"} の率で単位は $\mathrm{s^{-1}}$ である。弱結合で浴の相関時間が寿命よりずっと短ければ、率方程式 $\dot p_e=-A_{eg}p_e$ により寿命 $1/A_{eg}$ を得る。

非縮退二準位で等方的な光の単位角周波数当たりのエネルギー密度を $u_\omega$ とする。モード密度は $\omega^2/(\pi^2c^3)$ なので、平均占有 $\bar n$ に対し $u_\omega=\hbar\omega^3\bar n/(\pi^2c^3)$ である。誘導放出と吸収の行列要素は随伴であり、率は $A_{eg}\bar n$ となる。従って Einstein 係数は $B_{eg}=B_{ge}=A_{eg}\pi^2 c^3/(\hbar\omega_0^3)$ である。縮退度 $g_g,g_e$ の各内部準位を一様に占有するなら、行列要素の二乗和が上下の全組について同じなので $g_gB_{ge}=g_eB_{eg}$ となる。

:::example[双極子の放射パターン]
実の $\vb d$ に対して方向 $\theta$ は双極子と波数のなす角とすると $\dd A/\dd\Omega=(3A/8\pi)\sin^2\theta$ である。軸上に放射せず、垂直方向で最大となる。$\int\sin^2\theta\dd\Omega=8\pi/3$ で規格化を確認できる。
:::

:::example[熱平衡]
上向き率は $A\bar n$、下向き率は $A(\bar n+1)$ である。平衡確率比を $p_e/p_g=e^{-\beta\hbar\omega_0}$ として両流を等しくすると $\bar n/(1+\bar n)=e^{-\beta\hbar\omega_0}$、よって $\bar n=1/(e^{\beta\hbar\omega_0}-1)$ となる。
:::

### 自然幅の Lorentz 線

励起の振幅が Markov 近似で $c_e(t)=e^{-i\omega_0t-A t/2}$、$t\ge0$ とする。これは励起の確率減衰率が $A$ であることと整合する。放出モードの振幅には $\int_0^\infty e^{i(\omega-\omega_0)t-A t/2}\dd t=1/[A/2-i(\omega-\omega_0)]$ が掛かる。密度と結合を線幅の範囲で一定にすれば強度は Lorentz 型となり、角周波数の半値全幅は $A$、通常の周波数の半値全幅は $A/(2\pi)$ である。この近似は $\omega_0\gg A$ とし、負周波数への延長を規格化の便宜に使える。極端に短い時刻や連続体の端では指数近似から外れる。

:::example[寿命から線幅を予測する]
$A=6.3\times10^8\,\mathrm{s^{-1}}$ なら寿命は約 $1.59\,\mathrm{ns}$、通常の周波数の半値全幅は約 $100\,\mathrm{MHz}$ である。中心から $\Delta\omega=A$ だけ離れると、ピークに対する強度比は $(A/2)^2/[(A/2)^2+A^2]=1/5$ になる。線幅に含める緩和が自然放出だけという条件での比較である。
:::

### 一光子吸収の断面積と電離のしきい値

入射偏光を $\vb e$、原子を $\ket i$ とし、入射一光子が吸収されるとする。行列要素の二乗は $\hbar\omega|\mel f{\vb d\cdot\vb e}i|^2/(2\epsilon_0\mathcal V)$ である。黄金律の率を光子の流束 $c/\mathcal V$ で割ると
$$
\sigma_{\rm abs}(\omega)=\frac{\pi\omega}{\epsilon_0 c}
\sum_f|\mel f{\vb d\cdot\vb e}i|^2\delta(E_f-E_i-\hbar\omega)
$$
となる。電離では $\sum_f$ の連続状態を積分し、反跳なしで連続スペクトルの下端を零に取ると $\hbar\omega\ge -E_i$ が必要である。例えば基底の水素では電離エネルギーが約 $13.6\,\mathrm{eV}$ なので、$15.0\,\mathrm{eV}$ の光子を一個吸収した電子の運動エネルギーは約 $1.4\,\mathrm{eV}$ である。断面積には選択則と連続状態の行列要素も入る。

### 熱光と多重度を含む占有方程式

下と上の多重項をそれぞれ内部で一様とする。上の平均自然率を $A$、同じエネルギーのモード平均を $\bar n$ とする。上下の全行列要素の二乗和は同じなので下多重項の一準位当たりの平均吸収率は $(g_e/g_g)A\bar n$、上からの全放出率は $A(\bar n+1)$ である。全占有 $p_e$ は $\dot p_e=(g_e/g_g)A\bar n(1-p_e)-A(\bar n+1)p_e$ を満たす。定常比は $(g_e/g_g)\bar n/(\bar n+1)$ である。熱平衡 $\bar n=(e^{\beta\hbar\omega}-1)^{-1}$ を代入すると Gibbs 分布の多重度と一致する。偏極した光で内部占有を一様に維持できない場合は、個別の磁気準位を率方程式に残す。

:::example[個々の上準位と上多重項の占有を区別する]
$g_g=1,g_e=3$、$\bar n=1$ なら上向きは $3A$、下向きは $2A$ である。定常の上多重項は $3/5$、上の各内部準位は一様なら $1/5$、下は $2/5$ になる。下から始めると $p_e(t)=3(1-e^{-5At})/5$ である。各上準位の占有が下より小さくても、上の総占有は大きくなることがある。分光強度の釣り合いを議論するときはこの二種類を混同しない。
:::

### 浴の記憶を積分して指数を導く

一励起に制限し、相互作用表示で励起原子の振幅を $c_e(t)$、基底原子と一光子 $\lambda$ の振幅を $c_\lambda(t)$ とする。結合周波数を $g_\lambda$、離調を $\Delta_\lambda=\omega_\lambda-\omega_0$ とすると
$$
\dot c_e=-i\sum_\lambda \bar g_\lambda e^{-i\Delta_\lambda t}c_\lambda,\qquad
\dot c_\lambda=-i g_\lambda e^{i\Delta_\lambda t}c_e .
$$
初期 $c_\lambda(0)=0$ の第二式を積分して戻せば $\dot c_e(t)=-\int_0^t F(\tau)c_e(t-\tau)\dd\tau$、$F(\tau)=\sum_\lambda|g_\lambda|^2 e^{-i\Delta_\lambda\tau}$ である。連続体で $F$ が短時間に減り、$c_e$ の変化がそれより遅ければ、$c_e(t-\tau)\simeq c_e(t)$ として上限を無限へ伸ばす。$\int_0^\infty e^{-i\Delta\tau-\eta\tau}\dd\tau\to\pi\delta(\Delta)-i\operatorname{PV}(1/\Delta)$ であり、実部が $A/2$、虚部が準位の位相を変える。振幅は $e^{-A t/2}$ の絶対値、占有は $e^{-At}$ を得る。有限の箱でモードを離散のまま保つとこの速い相関減衰の仮定は一般に成り立たない。

:::example[一つのモードだけでは指数でない]
一モードで $\Delta_\lambda=0$ なら $F=|g|^2$ と一定である。積分式をもう一回微分すると $\ddot c_e=-|g|^2c_e$、$\dot c_e(0)=0$ なので $c_e=\cos(|g|t)$ になる。$t=\pi/(2|g|)$ で全励起が光へ移り、$t=\pi/|g|$ で原子側へ戻る。連続なモードの自然放出と同じ一励起・同じ結合から出発しても、環境のスペクトルによって長時間の解が変わる。
:::

## 双極子選択則と寿命

電気双極子は回転階数 1、パリティ奇である。以下で $q=-1,0,1$ は双極子演算子の球成分の添字である。:ref[thm-aqm-wigner-eckart] と :ref[prop-aqm-parity-selection] から $|J_f-J_i|\le1\le J_f+J_i$、$M_f-M_i=q$、$\pi_f=-\pi_i$ が必要である。特に $0\to0$ は消える。一電子の軌道で $\pi=(-1)^\ell$ なら三角条件と合わせて $\Delta\ell=\pm1$ である。スピンに作用しない双極子では、軌道とスピンが分離できる場合にスピン状態の内積がさらに課される。

複数の放出先に対し光子と原子の終状態が直交すれば確率率を足す。従って $\tau_i=(\sum_f A_{if})^{-1}$、分岐比 $A_{if}/\sum_{f'}A_{if'}$ になる。電気双極子で禁制でも、磁気双極子、電気四重極、二光子過程などの振幅が残り得る。

:::example[水素の $2p\to1s$]
Bohr 半径を $a_0$ とする。規格化動径関数は $R_{10}=2a_0^{-3/2}e^{-r/a_0}$、$R_{21}=(r/a_0)e^{-r/(2a_0)}/(2\sqrt6\,a_0^{3/2})$ である。動径積分 $\int r^3 R_{10}R_{21}\dd r=(256/(81\sqrt6))a_0$ と、$m=0$ での $\int\overline{Y_{00}}\cos\theta\,Y_{10}\dd\Omega=1/\sqrt3$ を掛ける。$|\vb d|=256 e a_0/(81\sqrt{18})$、$\hbar\omega_0\simeq10.2\,\mathrm{eV}$ で :eqref[eq-aqm-spontaneous-rate] は約 $6.3\times10^8\,\mathrm{s^{-1}}$ となる。微細構造を無視した寿命は約 $1.6\,\mathrm{ns}$ である。
:::

:::example[競合する二つの線]
同じ双極子強度の二経路で $\omega_2=2\omega_1$ なら $A_2=8A_1$ であり、分岐比は $1/9,8/9$ になる。全寿命は $1/(9A_1)$ である。
:::

### 磁気成分を分けて角度を変える

実の $z$ 双極子のパターンは $\sin^2\theta$ である。円偏光型の双極子 $\vb d=d(\hat{\vb x}+i\hat{\vb y})/\sqrt2$ なら $\abs{\hat{\vb k}\cdot\vb d}^2=|d|^2\sin^2\theta/2$ であり、和を取った放出は $(1+\cos^2\theta)/2$ に比例する。どちらも角積分は $8\pi/3$ なので、同じ $|d|$ と $\omega$ なら寿命は同じである。固定方向では横射影が異なるため、一方向の光量をそのまま全寿命に変換するには角分布を使う。

:::example[軸を中心とする円錐の割合]
$\Theta=\pi/3$ の円錐で $x=\cos\theta$ とすれば、直線双極子は $\frac34\int_{1/2}^1(1-x^2)\dd x=5/32$ である。円型では $\frac38\int_{1/2}^1(1+x^2)\dd x=19/64$ となる。二つの規格化した全放出率は同じでも、この検出器で見える円型の割合は直線型の $19/10$ 倍である。
:::

### 長波長展開で電気四重極を識別する

局所場の $E_i(\vb r)$ を原子中心の周りで展開する。電気的な結合を、原子中の座標 $\vb r_a$ に沿って $\int_0^1 E_i(s\vb r_a)\dd s$ と掛けた多極表示で考えると、第一項は $-\sum_a q_a r_{ai}E_i(0)$、次の項は $-\tfrac12\sum_a q_a r_{ai} r_{aj}\partial_j E_i(0)$ になる。$\int_0^1 s\dd s=1/2$ が係数である。真空の横波で $\divg\vb E=0$ なら $r_i r_j$ のトレース成分は縮約から消え、対称で跡なしの階数 2 が残る。座標を二つ掛けるからパリティは偶である。磁気双極子も長波長展開の次の次数に現れ得る。例えば局所的に一定な $\vb B$ に対し、対称ゲージ $\vb A=\vb B\times\vb r/2$ では最小結合の線形項が $-Q\vb L\cdot\vb B/(2m)$ となる。これは $\vb A\cdot\vb p=\vb B\cdot(\vb r\times\vb p)/2$ による。粒子にスピンの磁気モーメントがあればそれを軌道のモーメントに加え、全体を $\vb\mu$ と書くと結合は $-\vb\mu\cdot\vb B$ になる。$\vb\mu$ は階数 1 の軸性ベクトルなので、電気双極子と同じ三角条件で、パリティが同じ準位を結ぶ。光の磁場の一量子振幅は電場の $1/c$ で、横偏光を九十度回した方向にある。従って同じ角積分で、磁気双極子の自然率は $\mu_0\omega_0^3|\vb\mu_{ge}|^2/(3\pi\hbar c^3)$ となる。四重極の代表振幅は $k a$ 倍の尺度で、行列要素の大小は個別に確認する。

:::example[電気双極子が禁制でも四重極は候補になる]
同じパリティの $J_i=0$ から $J_f=2$ を考える。階数 1 は三角条件とパリティの両方に反し、階数 2 は両者を満たす。$M_i=0$ なら各 $q=-2,\ldots,2$ に $M_f=q$ が必要である。すべてが観測角で同じ強度になるという意味ではなく、横電場勾配のどの成分を照らすかが決める。例えば $a=0.1\,\mathrm{nm}$、$\lambda=500\,\mathrm{nm}$ では $ka\simeq1.26\times10^{-3}$ で、同程度の内部尺度で振幅の二乗の比は $1.6\times10^{-6}$ 程度になる。実際の寿命には磁気双極子や二光子過程との比較も入る。
:::

:::example[水素の 21 cm 線]
水素の $1s$ では、電子と陽子のスピンの三重項から一重項への遷移で軌道のパリティが変わらない。電気双極子は禁制だが、磁気双極子は許される。例えば $M=0$ の三重項 $(\ket{\uparrow\downarrow}+\ket{\downarrow\uparrow})/\sqrt2$ に電子の $\sigma_z$ を掛けると一重項になる。陽子の磁気モーメントを無視し、電子のモーメントを $-\mu_B\boldsymbol\sigma$ とすれば、磁気双極子の行列要素の大きさは $\mu_B=e\hbar/(2m)$ である。$\omega_0/(2\pi)\simeq1.42\,\mathrm{GHz}$ では自然放出率は約 $2.9\times10^{-15}\,\mathrm{s}^{-1}$、寿命は約 $1.1\times10^7$ 年である。磁気モーメントの小ささに加え、光学遷移より低い周波数が三乗で効く。
:::

### 振動子強度の和で遷移の規格化を検査する

一荷電粒子で $\Delta E_{ni}=E_n-E_i$ とし、方向 $x$ の振動子強度を $f_{ni}^{(x)}=(2m/\hbar^2)\Delta E_{ni}|\mel n x i|^2$ と定める。符号は $\Delta E_{ni}$ をそのまま使う。二重交換子 $[x,[H,x]]=\hbar^2/m$ と完全性から $\sum_n f_{ni}^{(x)}=1$ である。連続スペクトルがあれば和に積分も含む。基底からの吸収では各強度は非負で、分離した線の部分和が 1 を超えるなら規格化またはモデルの成分に再検査が必要である。励起状態からの和には負の下向き項があるので、上向きの部分和が 1 を超えることはこの和則と矛盾しない。電気双極子の $x$ 成分の放出率には $|\mel i{q x}n|^2$ が掛かるので、同じ強度を持つ遷移でも率は $\omega^2$ に比例して変わる。逆に $x$ 偏光の吸収では、上の断面積を一本の線の周りで角周波数積分すると $\pi Q^2 f_{ni}^{(x)}/(2\epsilon_0 mc)$ になる。デルタ関数を $\omega$ で積分すると $1/\hbar$ が出て、$\omega$ は $\Delta E_{ni}/\hbar$ で評価するからである。

:::example[等しい振動子強度の二本の線の寿命]
二本で比較する $|\Delta E|$ の比を 2 とし、同じ方向の振動子強度の絶対値が等しいとする。双極子の二乗は高い方で半分になる。自然放出率は $\omega^3 |d|^2$ なので、高い方の率は $2^3/2=4$ 倍である。同じ $d$ と仮定した場合は 8 倍であり、何を同じと置いたかで予測が違う。複数の偏光成分がある場合は $\sum_i|d_i|^2$ の換算にも同じ注意が必要である。
:::

## 光の散乱（Rayleigh 散乱・Thomson 散乱・Raman 散乱）

弱い単色電場の複素振幅を $\vb E_0$ とする。束縛電子の線形モデル $m(\ddot{\vb r}+\omega_b^2\vb r)=Q\vb E_0e^{-i\omega t}$ を解くと双極子振幅は $\vb d_0=Q^2\vb E_0/[m(\omega_b^2-\omega^2)]$ である。遠方放射場の横成分は $\omega^2[\vb d_0-\hat{\vb n}(\hat{\vb n}\cdot\vb d_0)]/(4\pi\epsilon_0 c^2r)$ に比例する。これは遅延ベクトルポテンシャルの双極子項 $\vb A\simeq\mu_0\dot{\vb d}(t-r/c)/(4\pi r)$ の時間微分と横射影から出る。入射と出射の平均 Poynting 流の比を取ると偏光 $\vb e,\vb e'$ に対し

$$
\frac{\dd\sigma}{\dd\Omega}
=\frac{\omega^4|\alpha_{\mathrm{pol}}(\omega)|^2}{16\pi^2\epsilon_0^2c^4}
 |\vb e'{}^*\cdot\vb e|^2,\qquad
\alpha_{\mathrm{pol}}=\frac{Q^2}{m(\omega_b^2-\omega^2)} .
$$

$\omega\ll\omega_b$ なら $\omega^4$ の :term[Rayleigh 散乱]{en="Rayleigh scattering" yomi="れいりーさんらん"}、$\omega_b=0$ なら周波数に依存しない :term[Thomson 散乱]{en="Thomson scattering" yomi="とむそんさんらん"} になる。自由電子では係数は $r_e^2$、$r_e=e^2/(4\pi\epsilon_0 mc^2)$ である。反跳 $\hbar\omega/(mc^2)$ を無視できる範囲を前提とする。

原子が別の状態 $\ket f$ で終わる :term[Raman 散乱]{en="Raman scattering" yomi="らまんさんらん"} では $\hbar\omega'=\hbar\omega+E_i-E_f>0$ が必要である。二次摂動 $H_I(E_{\rm in}-H_0+i0)^{-1}H_I$ に中間状態を挿入すると、偏光を含む原子部分は

$$
\sum_n\left[
\frac{\mel f{\vb d\cdot\vb e'{}^*}n\mel n{\vb d\cdot\vb e}i}{E_i+\hbar\omega-E_n+i0}
+\frac{\mel f{\vb d\cdot\vb e}n\mel n{\vb d\cdot\vb e'{}^*}i}{E_i-\hbar\omega'-E_n+i0}
\right].
$$

前者は吸収が先で光子が一時的に零個、後者は放出が先で二個になる。両項に各電場の因子 $\sqrt{\hbar\omega/(2\epsilon_0\mathcal V)}$ と $\sqrt{\hbar\omega'/(2\epsilon_0\mathcal V)}$ が掛かる。中間状態はエネルギー殻上にある必要がなく、両順序は振幅として加える。

:::example[Thomson 全断面積]
入射偏光が直線なら出射偏光和は $1-|\hat{\vb n}\cdot\vb e|^2$ である。積分は $8\pi/3$ なので $\sigma_T=8\pi r_e^2/3$。$r_e\simeq2.82\,\mathrm{fm}$ なら $\sigma_T\simeq6.66\times10^{-29}\,\mathrm{m^2}$ である。
:::

:::example[Stokes と反 Stokes の出射エネルギー]
物質の振動状態を $0.10\,\mathrm{eV}$ 励起する場合、$2.00\,\mathrm{eV}$ の光は $1.90\,\mathrm{eV}$ に下がる。逆に既に励起した振動を同じ差だけ下げれば $2.10\,\mathrm{eV}$ になる。二経路の強度には初期振動準位の存在確率も掛かる。
:::

### 一光子の断面積への換算

前の二順序の和を $\mathcal M_{fi}$ と置く。箱規格化で入射一光子から出射一光子への行列要素の絶対値二乗は $\hbar^2\omega\omega'|\mathcal M_{fi}|^2/(4\epsilon_0^2\mathcal V^2)$ である。一つの出射偏光を固定してモードを積分すると $\mathcal V k'^2\dd k'\dd\Omega'/(2\pi)^3$ が掛かる。Fermi の黄金律のデルタ関数で $\dd k'$ を積分すれば $1/(\hbar c)$ である。最後に一光子の入射流束密度 $c/\mathcal V$ で割ると

$$
\frac{\dd\sigma_{fi}}{\dd\Omega'}
=\frac{\omega\omega'^3}{16\pi^2\epsilon_0^2c^4}|\mathcal M_{fi}|^2 .
$$

$\mathcal M_{fi}$ の単位は双極子二乗をエネルギーで割った値であり、右辺は面積となる。上式は光子数の流束から求めた断面積である。光の出射エネルギー流束と入射エネルギー流束を比較するならさらに $\omega'/\omega$ が掛かる。

基底 $\ket i$ に戻る散乱で中間準位が一つ、遷移双極子が実である場合、偏光因子を除いた和の分母は $1/(-\Delta+\hbar\omega)+1/(-\Delta-\hbar\omega)=-2\Delta/(\Delta^2-\hbar^2\omega^2)$ となる。$\omega\ll\Delta/\hbar$ では有限な定数に収束し、断面積の $\omega^4$ が現れる。共鳴近くで分母が小さくなると中間準位の幅と高次の相互作用を含める必要がある。

### 最小結合での Thomson 項

自由な電子の低エネルギー散乱を最小結合で計算する。軌道運動の反跳をエネルギー保存では無視し、行列要素には入射と出射の波数の差だけを保つ。$Q^2\vb A^2/(2m)$ を :eqref[eq-aqm-vector-potential] で展開すると、入射光子を消し出射光子を作る項が二通りの順序で現れる。二つのモードが異なれば演算子は可換で、同じ偏光内積を生む。従って行列要素は電子の運動量保存因子を除いて
$$
\frac{Q^2\hbar}{2m\epsilon_0\mathcal V\sqrt{\omega\omega'}}\,
 \vb e'{}^*\cdot\vb e
$$
である。初期電子を静止とする。このとき $\vb p\cdot\vb A$ の第一頂点は零である。実際、横偏光なので $\vb p\cdot\vb A=\vb A\cdot\vb p$ であり、右側の運動量が初期平面波を消す。従ってこの近似では線形頂点二回の寄与を要さず、二次項が散乱を担う。電子のスピンの結合と有限の反跳は $\hbar\omega/(mc^2)$ による補正を生じる。

上の光子流束の計算にこの行列要素を入れ、$\omega'=\omega$ とすれば
$$
\frac{\dd\sigma}{\dd\Omega'}=
\left(\frac{Q^2}{4\pi\epsilon_0 mc^2}\right)^2
|\vb e'{}^*\cdot\vb e|^2
$$
となる。古典的な自由電子の結果と一致する。束縛電子では運動量の線形頂点が一般に非零なので、同じゲージと近似次数で両方の寄与を扱う。

### 振動 Raman 線の温度依存

振動を角周波数 $\Omega$ の調和振動子とし、偏極率を核座標 $Q_{\rm vib}$ で一次展開する。定数項は振動量子数を変えない。一次項は $Q_{\rm vib}=\sqrt{\hbar/(2M\Omega)}(b+b^\dagger)$ なので Stokes では $n+1$、反 Stokes では $n$ が強度に掛かる。非共鳴で同じ偏極率微分を用い、熱平衡の幾何分布 $p_n=(1-e^{-\beta\hbar\Omega})e^{-n\beta\hbar\Omega}$ で平均すれば $\bar n/(\bar n+1)=e^{-\beta\hbar\Omega}$ となる。光子数断面積の $\omega\omega'^3$ も保持すると両線の比は
$$
\frac{\sigma_{\rm aS}}{\sigma_{\rm S}}=
\left(\frac{\omega+\Omega}{\omega-\Omega}\right)^3e^{-\beta\hbar\Omega}.
$$
エネルギー強度の比ではもう一乗が掛かる。

:::example[二種類の検出量の比]
$\hbar\omega=2\,\mathrm{eV},\hbar\Omega=0.1\,\mathrm{eV}$、$e^{-\beta\hbar\Omega}=1/4$ とする。光子数の強度比は $(21/19)^3/4\simeq0.338$、エネルギー強度では $(21/19)^4/4\simeq0.373$ となる。励起準位の熱占有だけなら $1/4$ だが、出射光の状態密度と電場の振幅が異なる。
:::

### 自由電子の偏光平均と反跳の条件

入射波数の方向を $z$ とし、非偏極の入射で $\hat{\vb x},\hat{\vb y}$ を等確率とする。出射の偏光和 $1-|\hat{\vb n}\cdot\vb e|^2$ を平均すると $1-\sin^2\theta/2=(1+\cos^2\theta)/2$ である。従って $\dd\sigma/\dd\Omega=r_e^2(1+\cos^2\theta)/2$ となる。有限反跳があるときに光の周波数を求めるだけなら、静止電子の四元運動量 $p$、入出射光子 $k,k'$ を使う。$(p+k-k')^2=p^2$ から $p\cdot(k-k')=k\cdot k'$ である。エネルギーと角度を入れると
$$
\omega'=\frac{\omega}{1+(\hbar\omega/(mc^2))(1-\cos\theta)}
$$
である。Thomson の角分布はこの反跳が小さい極限のもので、有限な比での偏光と断面積の完全な補正は相対論的な相互作用を使う。

:::example[可視光の角度で反跳の大きさを計算する]
$E_\gamma=2.5\,\mathrm{eV}$、電子の静止エネルギーを $511\,\mathrm{keV}$ とする。後方での $\hbar\omega'/(\hbar\omega)$ は $[1+5/(5.11\times10^5)]^{-1}$ であり、失う比は約 $9.8\times10^{-6}$ である。直角なら先頭はその半分となる。断面積の Thomson 近似とこの波長のずれを、同じ小パラメータで評価できる。
:::

### 偏光で Raman の散乱経路が選ばれる

非共鳴で入出射の偏光を明示した有効な偏極率テンソル $\alpha_{ij}(Q)$ を考える。振動の一量子遷移の振幅は $\sum_{ij}\bar e'_i\,(\partial_Q\alpha_{ij})e_j\langle n+1|Q|n\rangle$ に比例する。偏極率そのものと、その振動座標に対する微分は別の行列である。例えば振動の対称性によって微分が対角の差の形になると、同じ方向の入出射と直交した方向で強度の違いが出る。全偏光和で測る場合は横波の二方向だけを足し、出射の進行方向と平行な成分を独立に数えない。

:::example[四十五度の偏光で交差成分が現れる]
入出射がともに $z$ 軸に沿う配置を考え、$xy$ ブロックの $\partial_Q\alpha$ を $\diag(A,-A)$ とする。入射が $\hat x$ なら同じ偏光に $A$、直交 $\hat y$ には零である。入射を $(\hat x+\hat y)/\sqrt2$ にすれば、テンソルを掛けたベクトルは $A(\hat x-\hat y)/\sqrt2$ となる。従って同偏光は零、直交偏光にすべての振幅が入る。別の振動でテンソルが $AI$ に比例するなら、入射をどの角度にしても逆の結果になる。強度の違いは振動の占有だけでなく、どのテンソル成分と偏光を組んだかを示す。
:::

## Lamb シフトの概観

:term[Lamb シフト]{en="Lamb shift" yomi="らむしふと"} は、原子と量子化された電磁場の相互作用による準位のずれである。例えば非相対論的な自己エネルギーでは真空から一光子を作る中間状態が寄与する。最小結合の $-(Q/m)\vb p\cdot\vb A$ で双極子近似を使うと、二次補正は角度と偏光を積分して

$$
\Delta E_n=
-\frac{Q^2}{6\pi^2\epsilon_0m^2 c^3}
\sum_m|\vb p_{mn}|^2\,
  \operatorname{PV}\int_0^\Lambda\frac{\omega\,\dd\omega}{\omega+(E_m-E_n)/\hbar}.
$$

ここで $\Lambda$ は角周波数 cutoff である。係数にはベクトルポテンシャルの $\omega^{-1}$、モード密度の $\omega^2$、分母の $\hbar^{-1}$ が入る。積分の線形発散に当たる自由粒子の質量補正を差し引くと、$\omega/(\omega+\Omega)-1=-\Omega/(\omega+\Omega)$ から

$$
\Delta E_n^{\rm sub}=
\frac{Q^2}{6\pi^2\epsilon_0m^2c^3\hbar}
\sum_m|\vb p_{mn}|^2(E_m-E_n)\log\left|\frac{\Lambda+(E_m-E_n)/\hbar}{(E_m-E_n)/\hbar}\right|.
$$

これは cutoff 付きの低エネルギー模型で、励起状態の実部には主値を使う。零エネルギー差の項は極限で零である。対数を共通の代表値 $L$ で近似すると、加重和は $\frac12\sum_i\langle[p_i,[H_{\rm at},p_i]]\rangle=\frac{\hbar^2}{2}\langle\nabla^2 V\rangle$ となる。二重交換子は完全性を挟むと $2(E_m-E_n)$ を掛け、座標表示では $[V,p_i]=i\hbar\partial_i V$ だからである。$V=-Ze^2/(4\pi\epsilon_0r)$ の Laplacian は $(Ze^2/\epsilon_0)\delta^{(3)}(\vb r)$ となり、$s$ 状態に大きな正のずれが出る。相対論的な高周波側との接続、真空偏極と実際の対数平均を含む精密値の計算は再正規化を要するため :cite[bjorken-drell] に委ねる。

:::example[準位による違い]
Coulomb 固有関数では $R_{n\ell}(r)$ が原点で $r^\ell$ に比例する。従って接触項の期待値は $2p$ で零、$2s$ では正である。後者は $|\psi_{2s}(0)|^2=1/(8\pi a_0^3)$ で、$1s$ の $1/8$ となる。
:::

:::example[エネルギー尺度]
周波数差が $1\,\mathrm{GHz}$ なら $h\nu\simeq4.14\times10^{-6}\,\mathrm{eV}$ である。水素の原子内遷移の eV 尺度に比べ小さいが、[Dirac の Coulomb 問題](/physics/advanced-quantum-mechanics/dirac-applications/) で同じエネルギーとなる $2s_{1/2}$ と $2p_{1/2}$ を区別できる。
:::

### 準位のずれを寿命と同じ分母から読む

励起準位の自己エネルギーに $(E_n-E_m-\hbar\omega+i0)^{-1}$ が現れる。分布の恒等式を使うと、実部は主値、虚部は $-\pi\delta(E_n-E_m-\hbar\omega)$ である。放出先があれば、虚部は正の確率減衰率 $A$ に対して $-\hbar A/2$ になる。実際、黄金律では同じ二乗行列要素とデルタに $2\pi/\hbar$ が掛かるからである。このため $E_n+\Delta E_n-i\hbar A/2$ を振幅の指数に入れると位相のずれと指数減衰が同時に出る。基底では放出先のデルタがなく、仮想状態から実部が出てもこの次数の寿命幅は生じない。

:::example[接触項と原子番号の関係]
非相対論的な水素型イオンでは $\psi_{ns}(0)$ の二乗は $Z^3/(\pi a_0^3 n^3)$ である。Coulomb の Laplacian の係数は $Z$ に比例するので、同じ代表対数を仮定した接触近似では $Z^4/n^3$ に比例する。同じ $2s$ を $Z=1,2$ で比べると近似の倍率は 16 であり、$1s$ と $2s$ の比は 8 である。実際の対数や質量、有限核や相対論の補正はこの冪則だけでは決まらない。
:::

### 高周波の質量項を引く前と後の積分

$\Omega=(E_m-E_n)/\hbar$ とする。cutoff の積分は $\int_0^\Lambda\omega/(\omega+\Omega)\dd\omega=\Lambda-\Omega\log|(\Lambda+\Omega)/\Omega|$ である。$\Omega<0$ なら内部の特異点の両側を同じ幅で切って主値を取る。差し引く線形項は $\sum_m|\mel m{\vb p}n|^2=\langle n|\vb p^2|n\rangle$ を伴う。従ってこれは一粒子の $\vb p^2$ の係数のずれ、すなわち質量の再定義に対応する。一方、残る項では重み付きの準位差が必要で、全体を振動子数のような単純な変数に置き換えられない。周波数 cutoff の変化が残るのは、低周波だけを使った計算なので高周波側の補償をまだ行っていないためである。

:::example[振動子で上下二つの分母を残す]
軌道が一方向の振動子なら、$p$ は上・下の二準位だけを結び、$\Omega$ は $+\omega_b,-\omega_b$ となる。$\Lambda>\omega_b$ とする。加重和の対数部分は共通因子 $\hbar^2m\omega_b^2/2$ を除き
$$
(n+1)\log\frac{\Lambda+\omega_b}{\omega_b}
-n\log\frac{\Lambda-\omega_b}{\omega_b}.
$$
$\Lambda/\omega_b$ が大きければ先頭は $\log(\Lambda/\omega_b)$、次の項は $(2n+1)\omega_b/\Lambda$ である。例えば $n=1,\Lambda/\omega_b=10$ では $2\log11-\log9$ となり、先頭の $\log10$ に正の補正が加わる。上と下の符号を共に正と扱うと、先頭に誤って $2n+1$ が掛かる。
:::

### 微細構造と同じ次元で輻射の尺度を読む

$V=-Z e^2/(4\pi\epsilon_0r)$ とし、接触近似の代表対数を $L$ とする。Laplacian と水素型 $s$ の原点値を先の式に入れ、$\alpha=e^2/(4\pi\epsilon_0\hbar c)$、$a_0=\hbar/(mc\alpha)$ を使うと、$ns$ の接触シフトは
$$
\Delta E_{\rm contact}\simeq \frac{4}{3\pi}\alpha(Z\alpha)^4\frac{mc^2}{n^3}L
$$
となる。実際、$q^2/(6\pi^2\epsilon_0 m^2c^3\hbar)$ と $\hbar^2/2$ の積に $Ze^2/\epsilon_0$ および $Z^3/(\pi a_0^3 n^3)$ を掛ける。Dirac の微細構造の代表尺度 $(Z\alpha)^4mc^2$ と比べると、さらに $\alpha$ と対数・係数が掛かる。$L$ は実際には全中間準位の重みや高エネルギーの補正で決まり、一つの経験的な定数をすべての原子に使えるとは限らない。

:::example[対数を固定した尺度の推定]
水素で $n=2$ とし、$L=8$ を尺度の例として入れると係数の $L/n^3$ は 1 である。$mc^2\simeq511\,\mathrm{keV}$、$\alpha\simeq1/137$ を使うと約 $4.5\,\mathrm{\mu eV}$ になる。$\Delta E/h$ は約 $1.1\,\mathrm{GHz}$ の尺度である。これは接触のみの値であり、$2p$ の非接触の寄与、有限核の大きさ、真空偏極、高次の項まで含めた分光値としては使わない。
:::

## 演習問題

::::exercise{level=basic}
角周波数が同じ二モードに光子を合計 2 個入れる。零点を除くエネルギーと状態数を求めよ。
:::solution
一光子当たり $\hbar\omega$ なので全励起エネルギーは $2\hbar\omega$ である。占有は $(2,0),(1,1),(0,2)$ の 3 状態である。
:::
::::

::::exercise{level=basic}
光子数状態 $\ket3$ の $g^{(2)}(0)$ を求めよ。
:::solution
$aa\ket3=\sqrt6\ket1$ なので分子は 6、分母は $3^2=9$ である。比は $2/3$ となる。
:::
::::

::::exercise{level=basic}
双極子が $z$ 方向の場合、放出の半頂角 $\pi/2$ の円錐に入る割合を求めよ。
:::solution
角分布は $\sin^2\theta$ であり、$\theta$ と $\pi-\theta$ の重みは測度 $\sin\theta\dd\theta$ も含めて等しい。よって上半球への割合は $1/2$ である。
:::
::::

::::exercise{level=basic}
$\omega_0$ を一定として双極子行列要素を 2 倍にする。自然放出寿命は何倍になるか。
:::solution
率は $|\vb d|^2$ に比例し 4 倍になる。自然放出が唯一の崩壊なら寿命は率の逆数なので $1/4$ 倍である。
:::
::::

::::exercise{level=advanced}
等方的な放射場で $\bar n=2$ の二準位原子の定常励起確率を求めよ。
:::solution
吸収率は $2A$、放出率は $3A$ である。$2A(1-p_e)=3Ap_e$ を解くと $p_e=2/5$ を得る。
:::
::::

::::exercise{level=advanced}
直線偏光での自由電子散乱について、出射方向と入射偏光の角度が $\pi/6$ のとき、最大方向との強度比を求めよ。
:::solution
偏光和は $\sin^2\theta$ なので指定方向の強度は $\sin^2(\pi/6)=1/4$ に比例する。最大方向は $\pi/2$ で重み 1 のため比は $1/4$ である。
:::
::::

::::exercise{level=advanced}
電気双極子で $J_i=1,M_i=1$ から $J_f=1$ に至る磁気量子数を各 $q$ について調べよ。パリティは反対とする。
:::solution
$q=1$ は $M_f=2$ となり範囲外である。$q=0$ は $M_f=1$、$q=-1$ は $M_f=0$ である。$1\otimes1$ の $J=1$ は交換反対称で、最高重み $(\ket{1,0}-\ket{0,1})/\sqrt2$ を降ろすと $(\ket{1,-1}-\ket{-1,1})/\sqrt2$ になる。従って二つの許される成分の強度はともに換算行列要素の二乗の $1/6$ である。
:::
::::

::::exercise{level=advanced}
$H=p^2/(2m)+m\Omega^2x^2/2$ の全固有状態について $\sum_m (E_m-E_n)|\mel m p n|^2$ を求めよ。
:::solution
二重交換子 $[p,[H,p]]=\hbar^2 m\Omega^2 I$ である。完全性を挿入した加重和はこの期待値の半分なので $\hbar^2m\Omega^2/2$ となり $n$ に依存しない。梯子では上側の寄与が $n+1$、下側の寄与が $-n$ に比例する。
:::
::::

::::exercise{level=advanced}
入射が $z$ 方向、入射偏光が $\hat{\vb x}$ の Thomson 散乱で、出射を $z$ 軸の周りの半頂角 $\Theta=\pi/3$ の円錐で集める。偏光を読まないときの全断面積に対する割合を求めよ。方位角の平均が、$z$ 向きの振動双極子の放出と同じかも調べよ。

:::solution
出射の単位方向を $\hat{\vb n}=(\sin\theta\cos\phi,\sin\theta\sin\phi,\cos\theta)$ とする。偏光和では射影 $I-\hat{\vb n}\hat{\vb n}^{\mathsf T}$ により $r_e^2(1-\sin^2\theta\cos^2\phi)$ となる。方位角を $0$ から $2\pi$ まで足すと $2\pi r_e^2(1-\sin^2\theta/2)=\pi r_e^2(1+\cos^2\theta)$ である。従って $x=\cos\theta$ に替えた円錐の断面積は
$$
\pi r_e^2\int_{1/2}^1(1+x^2)\dd x
=\pi r_e^2\left[\frac12+\frac{1-1/8}{3}\right]
=\frac{19\pi r_e^2}{24}.
$$
全角では $\pi r_e^2(2+2/3)=8\pi r_e^2/3$ なので、割合は $19/64$ である。$z$ 向きの振動双極子では方位角に依らず $\sin^2\theta$ となり、同じ円錐の割合は $\tfrac34\int_{1/2}^1(1-x^2)\dd x=5/32$ になる。入射の伝播は $z$ でも、Thomson の加速度は入射電場の $\hat{\vb x}$ 方向であるため、円錐の軸は振動方向と一致しない。
:::
::::

## 確認テスト

::::quiz
自由電磁場の非零波数における独立な偏光を選べ。

- [ ] 縦偏光だけ
- [x] 二つの横偏光
- [ ] 三つの空間成分が独立

:::explanation
Coulomb ゲージの $\vb k\cdot\vb A_{\vb k}=0$ により振幅は $\vb k$ に垂直な二次元平面にある。縦成分を独立に量子化しない。
:::
::::

::::quiz{answer="16" tol="1e-9"}
理想平板間隔を半分にしたとき、Casimir 圧力の絶対値の倍率を求めよ。
:::explanation
$P_{\mathrm C}\propto-d^{-4}$ である。したがって圧力の大きさは $(1/2)^{-4}=16$ 倍になる。
:::
::::

::::quiz
原子のエネルギーが増加する Raman 遷移の出射光を選べ。

- [x] 入射より低い周波数
- [ ] 入射より高い周波数
- [ ] 必ず同じ周波数

:::explanation
エネルギー保存は $\hbar\omega'=\hbar\omega-(E_f-E_i)$ である。$E_f>E_i$ の分だけ出射光のエネルギーが小さくなる。
:::
::::

::::quiz{answer="8" tol="1e-9"}
双極子行列要素が同じ二遷移で自然放出周波数の比が 2 である。高い周波数の率を低い周波数の率で割れ。
:::explanation
偏光和とモード密度を含めた率は $\omega_0^3|\vb d|^2$ に比例する。双極子強度が共通なので比は $2^3=8$ である。
:::
::::

---
title: "因果構造と特異点定理"
description: "因果関係と Cauchy 面を定義し、収束方程式から Penrose と Hawking の不完備性定理を導く。"
notation:
  - symbol: "I^+(p),J^+(p)"
    meaning: "点 p の時間的未来と因果的未来"
  - symbol: "\\theta"
    meaning: "測地線束の膨張率"
  - symbol: "\\sigma_{ab},\\omega_{ab}"
    meaning: "せん断と渦度"
  - symbol: "\\lambda"
    meaning: "光的測地線のアフィン径数"
---

## この章の目標

- 因果関係と大域的双曲性が初期値問題に与える意味を説明する。
- 測地線束の膨張、せん断、渦度を分解する。
- 焦点形成の時間上限を求める。
- Penrose と Hawking の定理が要求する仮定と結論を区別する。

## 時間的曲線・因果的曲線と因果関係

時間向きが選べる滑らかな Lorentz 時空を考える。区分的に滑らかで、各区間の非零接線が未来時間的である曲線を時間的曲線と呼ぶ。未来時間的または光的な接線をほとんど至る所で持つ局所 Lipschitz 曲線を :term[因果的曲線]{en="causal curve" yomi="いんがてききょくせん"} という。止まる区間は除き、補助弧長などで非零接線がほとんど至る所にあるように径数を取る。Lipschitz 曲線を含めるのは、因果曲線の列を収束させたり、端点で曲線を延長したりするためである。$p\ll q$ は $p$ から $q$ に時間的曲線があること、$p\le q$ は因果的曲線があるか $p=q$ であることを表す。これらの終点集合を $I^+(p),J^+(p)$ と書く。時間を逆にすると $I^-,J^-$ になる。曲線をつなぐと両関係の推移性が従う。時間的な終点は初期点・終点を少し変えても時間的に結べるので $I^+(p)$ は開集合である。局所の凸正規近傍でこれを確かめるには曲線の最後の時間的区間を取り、その端点の指数座標が時間円錐の内部であることを使えばよい。

Minkowski 平面で $p=(0,0)$、$q=(T,X)$ なら、光速の不等式を積分して $T\ge\abs X$ が必要である。逆に直線を取れば十分であり、$T>\abs X$ が時間的未来になる。例えば $(2,2)$ は $J^+(p)\setminus I^+(p)$、$(2,1)$ は $I^+(p)$ にある。一般の時空では局所円錐だけから遠方の $J^+$ は決まらない。集合内のどの二点も $p\ll q$ を満たさない場合、その集合を :term[非時間的関連集合]{en="achronal set" yomi="ひじかんてきかんれんしゅうごう"} と呼ぶ。ここでは $p=q$ も禁止する。さらに相異なる二点が因果的にも結ばれない場合、その集合を :term[非因果的関連集合]{en="acausal set" yomi="ひいんがてきかんれんしゅうごう"} と呼ぶ。

### 因果的に接続可能かを座標で判定する

Minkowski 空間で各未来因果曲線の空間速度は $|\vb v|\le1$ である。$p$ と $q$ を結ぶ曲線の座標時刻差を $\Delta t$ とすれば、三角不等式から $|\Delta\vb x|\le\int|\vb v|\dd t\le\Delta t$ である。逆にこの不等式と $\Delta t>0$ が成立すれば線分の速度は $|\Delta\vb x|/\Delta t\le1$ なので因果曲線を与える。時間的な場合は不等号が厳しい。従って $J^+(p)$ は閉円錐、$I^+(p)$ はその内部である。一般の時空で $J^+$ が閉になるかは、曲線が途中で時空の欠けた点に近づけるかにも依存する。

:::example[一点を除いた平面での因果関係]
原点を除いた Minkowski 平面で $p=(-1,-1),q=(1,1)$ を取る。光的に離れた二点を結ぶには、速度の積分が上の不等式で等号を取るので $x=t$ の直線に限られる。この直線は除いた原点を通るため $q\notin J^+(p)$ である。一方 $q_\epsilon=(1,1-\epsilon)$ に直線で進むと、その直線は原点を外れて時間的なので $q_\epsilon\in I^+(p)$ である。$\epsilon\to0^+$ で $q$ に近づくため $J^+(p)$ は閉でない。
:::

### 時間的な区間を含む因果経路

凸正規近傍において、近傍内で結べる点の因果関係は指数座標の円錐で判定できる。光的に離れた二点間の因果経路は一本の光的測地線に限られる。これは Lorentz 版の Gauss の補題を使い、円錐の内部では放射方向の固有時が曲線長を上から押さえ、円錐表面では因果的接線がその光的生成子と平行でなければならないことから分かる。例えば Minkowski の光円錐 $t=|\vb x|$ 上の接線 $(\dot t,\dot{\vb x})$ は $\dot t=\hat{\vb x}\cdot\dot{\vb x}$ を満たすので、因果的であるには Cauchy–Schwarz の等号が必要である。

この局所事実から $p\ll q\le r$ なら $p\ll r$ を示せる。$q$ から $r$ の曲線を有限個の凸正規近傍に収まる区間に分ける。最初の区間の始点を少しだけ近傍内の時間的過去に押すと、$I^+(p)$ の開性により押した点も $I^+(p)$ にある。区間に対し、元の始点が終点の光円錐表面にある場合も、時間的過去に押した点を終点の過去の円錐の内部に選べる。従って終点も $I^+(p)$ にある。同様に各区間の終点を順に含めればよい。時間を逆にすれば $p\le q\ll r$ にも同じ帰結がある。さらに因果曲線が折れ曲がった二つの異なる光的方向を含めば、折れ点の近くの両端は凸正規近傍で時間的に結べる。従って因果境界をたどる曲線には、時間的区間も光的な折れ点も現れない。

:::example[二つの光信号が張る経路]
平坦な二次元で $(0,0)$ から $(1,1)$、次に $(2,0)$ へ光で進む。それぞれの区間は光的だが、全体の端点の時刻差は $2$、位置差は $0$ なので、その二点の間には固有時 $2$ の静止した時間的経路がある。従って $(2,0)$ は $(0,0)$ の因果的未来の中でも内部にある。この例では因果経路の固有時が零でも、端点が光的に離れているとは限らない。
:::

## 大域的双曲性と Cauchy 面

:term[Cauchy 面]{en="Cauchy surface" yomi="こーしーめん"} はすべての延長不能な因果曲線とちょうど一度交わる超曲面である。強因果性とは、各点の任意に小さい近傍の中に、因果曲線が離れた後で戻らない近傍を選べることである。強因果性と全ての $p,q$ の因果ダイヤモンド $J^+(p)\cap J^-(q)$ のコンパクト性を満たす時空を :term[大域的双曲性]{en="global hyperbolicity" yomi="たいいきてきそうきょくせい"} があるという。Minkowski 時空では $t$ は未来曲線上で増加し、曲線の空間速度は $1$ 以下である。$t$ が有限で止まると空間位置も極限を持ち延長できる。従って各平面 $t=t_0$ は Cauchy 面である。一方、原点を抜いた時空で原点の前後を結ぶダイヤモンドには、原点に収束して時空内に極限を持たない列がある。

以下で使う大域的基本事実を述べる。大域的双曲時空は Cauchy 面を持ち、因果的に結べる二点間で固有時を最大にする因果測地線が存在する。Cauchy 面 $\Sigma$ があるとき、点 $p$ から未来側の $\Sigma$ へ到達できれば $J^+(p)\cap\Sigma$ はコンパクトで、最長曲線がある。閉集合としての因果境界の生成子やこれらの最大化定理には、因果曲線列の極限定理と固有時の上半連続性が必要なので、証明は :cite[hawking-ellis] を参照する。以下の特異点定理では、この大域的な存在結果と次節の局所的な収束計算を組み合わせる。

### 初期値の情報が届く範囲

空間的な初期面の閉部分 $S$ について、その点を通るすべての過去延長不能な因果曲線が $S$ と交わる点の集合を :term[未来依存領域]{en="future domain of dependence" yomi="みらいいぞんりょういき"} $D^+(S)$ とする。時間を反転したものを $D^-(S)$ とし、その和を $D(S)$ と書く。$S$ が非時間的関連集合であるとき、$H^+(S)=\overline{D^+(S)}\setminus I^-(D^+(S))$ を :term[未来 Cauchy 地平線]{en="future Cauchy horizon" yomi="みらいこーしーちへいせん"} とする。単に $D^+(S)$ の全境界を取ると、初期面の内部まで含む。これを除くため、依存領域の点に向かって時間的に進める点を取り除く。Hawking の不完備性定理で用いる Cauchy 面は、依存領域が時空全体になる場合である。

:::example[有限の初期区間]
Minkowski 平面で $S$ を $t=0$ の閉区間 $-L\le x\le L$ とする。未来の点 $(t,x)$ から過去に進み、$t=0$ で交わり得る位置は $[x-t,x+t]$ である。全てが $S$ に入る条件は $t\ge0$ と $\abs x+t\le L$ であり、上側の境界は $\abs x+t=L$ の光線になる。例えば $L=3$、$(t,x)=(2,1)$ は境界上、$(2,2)$ は依存領域外である。後者の過去向き光線の一つは $t=0$ で $x=4$ に到達し、初期区間と交わらない。
:::

### 時刻面が Cauchy 面になるか

完全な Riemann 空間 $(S,h)$ に対して $g=-\dd t^2+h$ の積を考える。未来因果曲線を $t$ でパラメータ化すれば空間速度の $h$ ノルムは $1$ 以下である。もし延長不能な未来曲線の $t$ の上端が有限なら、時刻が上端へ近づくとき空間成分は Cauchy 列となる。完全性により点へ収束し、そこから未来時間的な線分をつないで延長できる。従って上端は無限であり、過去側も同様である。各切片には単調性により一度だけ交わるため、Cauchy 面である。

:::example[空間の端から情報が入る場合]
空間を開区間 $0<x<1$ に替えると、$x=t$、$0<t<1$ の光線はどちらの端でも時空に点がないため延長できない。この曲線は $t=2$ の面に交わらない。$g$ の成分が内点で正則でも、空間の完全性を取り除くと上の証明の収束先を失う。
:::

### 因果階層を時間関数で検査する

滑らかな関数 $t$ の勾配が過去時間的なら、未来因果曲線の接線 $k$ に対して $\dd t(k)=g(\nabla t,k)>0$ となる。これは局所 Lorentz 枠で時間成分の積が空間成分の積より厳密に大きいことから従う。従って同一点に戻る閉じた未来因果曲線は存在しない。また各点で勾配が時間的という条件は計量を小さく変えても保たれる。各点で少し幅を広げた光円錐を持つ滑らかな計量を選び、その計量でも同じ関数が増加する場合、元の時空を安定因果的という。時間が増える条件だけでは、その時間面にすべての延長不能曲線が交わるとは結論できない。曲線が有限の $t$ で空間境界に出る可能性も調べる。

:::example[円筒に空間を巻く場合と時間を巻く場合]
$\dd s^2=-\dd t^2+\dd x^2$ で $x$ を長さ $L$ だけ同一視すれば、$t$ は滑らかな実数値の時間関数であり、各時刻面は完全な円である。完全な空間と時刻の直積の議論から Cauchy 面になる。逆に $t$ を周期 $L$ で同一視して $x$ をそのままにすると、$x=\mathrm{const}$ の一周が閉じた未来時間的曲線になる。どちらも局所の接続と曲率は零であり、因果階層を分けるのは局所計量の式でなく、大域的な同一視である。
:::

### 因果ダイヤモンドと最大固有時

ダイヤモンドがコンパクトであることは、点の列だけでなく曲線の列を扱う際に使う。二点 $p\ll q$ に対し、$p$ から $q$ の未来因果曲線はすべて $K=J^+(p)\cap J^-(q)$ に入る。強因果性を用いて $K$ を因果曲線が再入しない座標近傍で有限被覆する。各近傍は有界な座標箱に収め、座標時刻に対する空間速度を有界にできる大きさに取る。この速度評価はコンパクトな閉包上で円錐が座標時刻面に接しないことから従う。補助 Riemann 計量を一つ固定すれば、各近傍での曲線長は座標時刻の幅の定数倍で抑えられる。再入しないので総曲線長にも一様な上限がある。

長さの上限に近づく曲線の列を補助弧長に比例する径数で共通の閉区間に移すと、同じ Lipschitz 定数を持つ。そのため Arzelà–Ascoli により一様収束する部分列を取れる。極限曲線が因果的であることは、十分小さい凸正規近傍で近傍内の因果関係が閉であることを使えば局所的に調べられる。固有時の上半連続性を使うと、この極限が長さの最大値を取る。上半連続性の要点は、小区間に分け、各区間で端点を結ぶ放射測地線の固有時を上界にすることにある。分割を細かくすると、極限曲線の固有時はこれらの上界の和の下限になる。この最後の等式を局所 Lipschitz 曲線に対して正当化する因果曲線の長さ理論は :cite[hawking-ellis] を参照する。

最大曲線が各凸正規近傍で最大を取らなければ、その部分をより長い曲線で置き換えられる。従って各小区間は放射測地線であり、折れ点も作れない。端点 $p,q$ が時間的に結べるなら正の長さを持つので、最大曲線は時間的測地線になる。始点を空間的 Cauchy 面 $\Sigma$ 上で動かす場合は、可能な始点が $J^-(q)\cap\Sigma$ に収まる。点から Cauchy 面へのこの集合のコンパクト性を使えば、最大化する列の始点にも収束部分列を取れる。点の間の議論と合わせると、終点を固定し始点の位置を自由にした最大化が意味を持つ。

:::example[曲線長が任意に増える静的パッチ]
反 de Sitter の中心を通る動径面を $r=L\tan\rho$、$t=L\tau$ で書くと、動径計量は $L^2\sec^2\rho\,(-\dd\tau^2+\dd\rho^2)$ であり、外端は $\rho=\pi/2$ にある。中心の二点の時刻差を $\Delta\tau>\pi$ とし、外向きと内向きの光路を $\rho_R<\pi/2$ で折り返す。途中で静止できる時間は $\Delta\tau-2\rho_R>\Delta\tau-\pi>0$ である。その区間の固有時は $L\sec\rho_R(\Delta\tau-2\rho_R)$ だから $\rho_R\uparrow\pi/2$ で発散する。この経路は因果的である。初めと終わりを少し時間的な速度に変えても、$\Delta\tau-\pi$ の余裕があれば発散は残る。従ってこの時空では中心の二点間ですら有限の最大固有時を得られない。曲線が境界へ逃げるため、先のコンパクト性の議論を適用できない。
:::

### Cauchy 面への時間方向の射影

滑らかな未来時間的ベクトル場 $X$ は局所の時間的ベクトルを分割の一で貼り合わせて作れる。未来時間的円錐が凸なので貼り合わせ後も時間的である。完全な補助 Riemann 計量に対し、例えば $X/\sqrt{1+|X|^2}$ と再規格化すればノルムが $1$ 以下になる。流線が有限パラメータで終わるなら補助計量で Cauchy 列となり、滑らかな常微分方程式を到達点から続けられるので矛盾する。従って再規格化した場は完全である。

滑らかな空間的 Cauchy 面 $\Sigma$ があれば、各流線はその面を横切って一度だけ通る。実際、流線は時間的であり、時空内の一点へ漸近して終わることは、非零の場のフローボックスで座標に直して否定できる。このため流れ $\R\times\Sigma\to M$ は全単射である。また初期面で流れと接空間が横断的なので微分が非退化である。よって微分同相となり、$M$ が連結なら $\Sigma$ も連結である。以下の Penrose の証明では、$X$ の流線に沿うこの射影を用いる。$\Sigma$ は滑らかな空間的 Cauchy 面を選んでよいという大域的双曲性の標準結果に基づく。

## 測地線の束と Raychaudhuri 方程式

時間的な単位測地線束の接線を $u$ とし、$u^\perp$ で $B_{ab}=h_a{}^ch_b{}^d\nabla_d u_c$ を定める。$\theta=\tr B$、$\sigma_{ab}=B_{(ab)}-\theta h_{ab}/3$、$\omega_{ab}=B_{[ab]}$ がそれぞれ :term[膨張率]{en="expansion scalar" yomi="ぼうちょうりつ"}、せん断、渦度である。測地線束の分離に $\nabla_u\xi=B\xi$ が成り立つ。平行輸送した枠で二度微分すると $(\dot B+B^2)\xi$ なので、:eqref[eq-gr-jacobi] と比較して $\dot B^i{}_j=-B^i{}_k B^k{}_j-R^i{}_{0j0}$ である。跡を取ると

$$
\dot\theta=-\frac{\theta^2}{3}-\sigma_{ab}\sigma^{ab}+\omega_{ab}\omega^{ab}-R_{ab}u^a u^b .
\label{eq-gr-raychaudhuri-timelike}
$$

反対称部分の行列二乗の跡は負になるため渦度項が正である。光的測地線束では、$k^\perp$ で $k$ 自身の方向を除いた二次元の横断空間を用いる。接線 $k$ がアフィンであれば同じ行列計算の跡が

$$
\dv{\theta}{\lambda}=-\frac{\theta^2}{2}-\sigma_{ab}\sigma^{ab}+\omega_{ab}\omega^{ab}-R_{ab}k^a k^b
\label{eq-gr-raychaudhuri-null}
$$

を与える。これらが :term[Raychaudhuri 方程式]{en="Raychaudhuri equation" yomi="らいちゃうどぅりーほうていしき"} である。横方向の Jacobi 行列を $A$ とすれば $B=\dot A A^{-1}$、$\theta=(\log\det A)^{\cdot}$ となるので、面積または体積の局所的な拡大率を表す。例えば平坦時空の外向き球面光束は $A\propto r I$、$\dd r/\dd\lambda=1$ であり、$\theta=2/r$ は $\theta'=-\theta^2/2$ を満たす。

### 光的束の横断計量と規格化

光的接線 $k$ の直交空間では、$g(k,\cdot)=0$ なので $v$ と $v+ak$ は全ての横断内積で同じである。補助光的ベクトル $\ell$ を $k\cdot\ell=-1$ と取り、$v+(v\cdot\ell)k$ を代表に選ぶと $k,\ell$ の両方に直交する。従って商 $k^\perp/\R k$ の計量は二次元で正定値である。$k$ がアフィン光的測地線の接線なら、$\nabla_k k=0$ と $\nabla_v(k^2)=0$ から $[v]\mapsto[\nabla_vk]$ がこの商の作用素になる。補助 $\ell$ は代表を選ぶだけなので、跡とせん断のノルムは商の上だけで定義できる。

地平線の Killing 時刻など、アフィンでない径数を使う場合は $\ell_1=bk$、$b>0$ と置く。$\nabla_{\ell_1}\ell_1=k(b)\ell_1$、横射影では $\nabla_v(bk)$ の $v(b)k$ が消えるので $\theta_1=b\theta$ である。そこで $\ell_1(\theta_1)=b^2 k(\theta)+k(b)\theta_1$ となり、面積の式にも同じ径数変換が入る。例えば Minkowski の発散する球面光束でアフィン径数 $r>0$ の代わりに $s=\log(r/R)$ を使うと $\ell_1=rk$、$\theta_1=2$ で一定である。非アフィン方程式の $-\theta_1^2/2=-2$ と $k(b)\theta_1=2$ が相殺する。この例のように、径数の変更によって膨張率の数値が変わっても、面積変化の符号は保たれる。

### 加速する束と収束

単位時間的場 $u$ の加速度を $a^a=u^b\nabla_b u^a$ とする。$\theta=\nabla_a u^a$ を速度に沿って微分し、共変微分の順を交換すると
$$
u(\theta)=\nabla_a a^a-(\nabla_a u^b)(\nabla_b u^a)-R_{ab}u^au^b .
$$
$\nabla_b u^a=B^a{}_b-a^au_b$ であり、$B$ は両側で $u$ と直交し $a\cdot u=0$ なので積の跡は $\tr B^2$ である。従って測地線束と比べて発散 $\nabla_a a^a$ が追加される。例えば Rindler の静止束では $u=(1+ax)^{-1}\partial_t$、$a^x=a/(1+ax)$ である。体積密度 $1+ax$ を使うと $\nabla_a a^a=(1+ax)^{-1}\partial_x[(1+ax)a^x]=0$ であり、$\theta=0$ と平坦性に一致する。

:::example[宇宙の共動観測者で膨張率と Ricci 収束を照合する]
平坦な共動切片の $a(t)>0$ の計量で $u=\partial_t$ とする。座標速度が零でも空間の距離は時間とともに変わるため、速度勾配と固有加速度は別々に計算する。 座標の順序を $(t,x,y,z)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & a^{2}{\left(t \right)} & 0 & 0\\0 & 0 & a^{2}{\left(t \right)} & 0\\0 & 0 & 0 & a^{2}{\left(t \right)}\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & \frac{1}{a^{2}{\left(t \right)}} & 0 & 0\\0 & 0 & \frac{1}{a^{2}{\left(t \right)}} & 0\\0 & 0 & 0 & \frac{1}{a^{2}{\left(t \right)}}\end{matrix}\right]\\
\det g&=- a^{6}{\left(t \right)}
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{t}_{x x}$ | $- \frac{1}{2}\left(-\partial_{t}\left(a^{2}{\left(t \right)}\right)\right)$ | $a{\left(t \right)} a'{\left(t \right)}$ |
| $\Gamma^{t}_{y y}$ | $- \frac{1}{2}\left(-\partial_{t}\left(a^{2}{\left(t \right)}\right)\right)$ | $a{\left(t \right)} a'{\left(t \right)}$ |
| $\Gamma^{t}_{z z}$ | $- \frac{1}{2}\left(-\partial_{t}\left(a^{2}{\left(t \right)}\right)\right)$ | $a{\left(t \right)} a'{\left(t \right)}$ |
| $\Gamma^{x}_{t x}$ | $\frac{1}{2 a^{2}{\left(t \right)}}\left(\partial_{t}\left(a^{2}{\left(t \right)}\right)\right)$ | $\frac{a'{\left(t \right)}}{a{\left(t \right)}}$ |
| $\Gamma^{y}_{t y}$ | $\frac{1}{2 a^{2}{\left(t \right)}}\left(\partial_{t}\left(a^{2}{\left(t \right)}\right)\right)$ | $\frac{a'{\left(t \right)}}{a{\left(t \right)}}$ |
| $\Gamma^{z}_{t z}$ | $\frac{1}{2 a^{2}{\left(t \right)}}\left(\partial_{t}\left(a^{2}{\left(t \right)}\right)\right)$ | $\frac{a'{\left(t \right)}}{a{\left(t \right)}}$ |

Ricci テンソルは第一・第三添字を縮約して得る。曲率が残る場合も相殺が起こるので、各項を落とさずに加える。

| 成分 | 縮約の項 | 代入 | 結果 |
|---|---|---|---|
| $R_{t t}$ | $R^{x}_{t x t}+R^{y}_{t y t}+R^{z}_{t z t}$ | $- \frac{a''{\left(t \right)}}{a{\left(t \right)}} - \frac{a''{\left(t \right)}}{a{\left(t \right)}} - \frac{a''{\left(t \right)}}{a{\left(t \right)}}$ | $- \frac{3 a''{\left(t \right)}}{a{\left(t \right)}}$ |
| $R_{x x}$ | $R^{t}_{x t x}+R^{y}_{x y x}+R^{z}_{x z x}$ | $a{\left(t \right)} a''{\left(t \right)} + \left(a'{\left(t \right)}\right)^{2} + \left(a'{\left(t \right)}\right)^{2}$ | $a{\left(t \right)} a''{\left(t \right)} + 2 \left(a'{\left(t \right)}\right)^{2}$ |
| $R_{y y}$ | $R^{t}_{y t y}+R^{x}_{y x y}+R^{z}_{y z y}$ | $a{\left(t \right)} a''{\left(t \right)} + \left(a'{\left(t \right)}\right)^{2} + \left(a'{\left(t \right)}\right)^{2}$ | $a{\left(t \right)} a''{\left(t \right)} + 2 \left(a'{\left(t \right)}\right)^{2}$ |
| $R_{z z}$ | $R^{t}_{z t z}+R^{x}_{z x z}+R^{y}_{z y z}$ | $a{\left(t \right)} a''{\left(t \right)} + \left(a'{\left(t \right)}\right)^{2} + \left(a'{\left(t \right)}\right)^{2}$ | $a{\left(t \right)} a''{\left(t \right)} + 2 \left(a'{\left(t \right)}\right)^{2}$ |

逆計量と縮約すると Ricci スカラーとなる。さらに $G_{ab}=R_{ab}-g_{ab}R/2$ を成分ごとに計算する。

$$
\begin{aligned}
R&=\sum_{a,b}g^{ab}R_{ab}\\
&=\frac{a{\left(t \right)} a''{\left(t \right)} + 2 \left(a'{\left(t \right)}\right)^{2}}{a^{2}{\left(t \right)}} + \frac{a{\left(t \right)} a''{\left(t \right)} + 2 \left(a'{\left(t \right)}\right)^{2}}{a^{2}{\left(t \right)}} + \frac{a{\left(t \right)} a''{\left(t \right)} + 2 \left(a'{\left(t \right)}\right)^{2}}{a^{2}{\left(t \right)}} + \frac{3 a''{\left(t \right)}}{a{\left(t \right)}}\\
&=\frac{6 \left(a{\left(t \right)} a''{\left(t \right)} + \left(a'{\left(t \right)}\right)^{2}\right)}{a^{2}{\left(t \right)}}
\end{aligned}
$$

非零の Einstein テンソル成分は次のとおりである。

| 成分 | Ricci 項 | 跡の項 | 差 |
|---|---|---|---|
| $G_{t t}$ | $- \frac{3 a''{\left(t \right)}}{a{\left(t \right)}}$ | $- \frac{3 \left(a{\left(t \right)} a''{\left(t \right)} + \left(a'{\left(t \right)}\right)^{2}\right)}{a^{2}{\left(t \right)}}$ | $\frac{3 \left(a'{\left(t \right)}\right)^{2}}{a^{2}{\left(t \right)}}$ |
| $G_{x x}$ | $a{\left(t \right)} a''{\left(t \right)} + 2 \left(a'{\left(t \right)}\right)^{2}$ | $3 a{\left(t \right)} a''{\left(t \right)} + 3 \left(a'{\left(t \right)}\right)^{2}$ | $- 2 a{\left(t \right)} a''{\left(t \right)} - \left(a'{\left(t \right)}\right)^{2}$ |
| $G_{y y}$ | $a{\left(t \right)} a''{\left(t \right)} + 2 \left(a'{\left(t \right)}\right)^{2}$ | $3 a{\left(t \right)} a''{\left(t \right)} + 3 \left(a'{\left(t \right)}\right)^{2}$ | $- 2 a{\left(t \right)} a''{\left(t \right)} - \left(a'{\left(t \right)}\right)^{2}$ |
| $G_{z z}$ | $a{\left(t \right)} a''{\left(t \right)} + 2 \left(a'{\left(t \right)}\right)^{2}$ | $3 a{\left(t \right)} a''{\left(t \right)} + 3 \left(a'{\left(t \right)}\right)^{2}$ | $- 2 a{\left(t \right)} a''{\left(t \right)} - \left(a'{\left(t \right)}\right)^{2}$ |

アフィン径数の微分を点で表し、$w_a$ を座標 $a$ の速度成分とする。ここでは添字を下げる操作を行っていない。測地線の右辺 $-\Gamma^a{}_{bc}w_bw_c$ では $b\ne c$ の項を二度数える。

$$
\begin{aligned}
\ddot{t}&=- \left(w_{x}^{2} + w_{y}^{2} + w_{z}^{2}\right) a{\left(t \right)} a'{\left(t \right)}\\
\ddot{x}&=- \frac{2 w_{t} w_{x} a'{\left(t \right)}}{a{\left(t \right)}}\\
\ddot{y}&=- \frac{2 w_{t} w_{y} a'{\left(t \right)}}{a{\left(t \right)}}\\
\ddot{z}&=- \frac{2 w_{t} w_{z} a'{\left(t \right)}}{a{\left(t \right)}}
\end{aligned}
$$

また計量ノルムは

$$
g(\dot x,\dot x)=- w_{t}^{2} + w_{x}^{2} a^{2}{\left(t \right)} + w_{y}^{2} a^{2}{\left(t \right)} + w_{z}^{2} a^{2}{\left(t \right)}
$$
である。これを径数で微分して上の加速度を代入すると、係数の微分と速度の微分の和が零になる。従って初期のノルムが保存される。

各座標を固定しても $\Gamma^a{}_{tt}=0$ なので $u$ は測地線である。$\nabla_j u^i=H\delta^i_j$ で、せん断と渦度は零、$\theta=3H$ となる。従って $-\theta^2/3-R_{tt}=-3H^2+3a''/a=3H'$ で、直接微分と一致する。$a(t)=e^{Ht}$ では $H$ が一定であり、$R_{tt}<0$ が膨張率二乗の項を正確に打ち消す。
:::

## 共役点と焦点

測地線上の二点で零になる非零の Jacobi 場がある場合、後の点は初めの点の :term[共役点]{en="conjugate point" yomi="きょうやくてん"} である。面から直交して出る測地線族について法線指数写像の微分が退化する点を :term[焦点]{en="focal point" yomi="しょうてん"} という。横の Jacobi 行列が $\det A=0$ となることに対応する。直交法線束では渦度が零である。実際、局所で法線一形式はレベル面関数の微分の倍であり、その外微分の横射影は零になる。Ricci の収束条件があると十分滑らかな区間で $\dot\theta\le-\theta^2/d$、$d=3$ または $2$ である。$\theta_0<0$ のとき $1/\theta$ の微分は $-\dot\theta/\theta^2\ge1/d$ である。負のまま続くなら $\lambda-\lambda_0=d/\abs{\theta_0}$ までに矛盾する。測地線がその先まで正則なら、その前またはその点で Jacobi 行列が退化する。

焦点を通過した時間的測地線が最大でない理由は第二変分で分かる。直交変分の指標形式は $I(V,V)=\int[g(D_\tau V,D_\tau V)-g(R(V,u)u,V)]\dd\tau$ と始点面に由来する境界項の和で、固有時の第二変分はその負である。Jacobi 場 $J$ を最初の焦点まで取り、その後に零で延長すると、部分積分と始点の直交変分条件により $I(J,J)=0$ となるが、$I(J,W)=g(D_\tau J,W)$ の焦点境界項が任意の試験場 $W$ に対して残る。焦点で $D_\tau J\ne0$ なので $W$ を選んでこれを負にできる。$V=J+\eps W$ では小さい正の $\eps$ で $I(V,V)<0$ となる。滑らかな変分へ近似しても符号は保たれる。光束でも以下で横断指標形式を使い、焦点後の生成子が因果境界から外れることを示す。定曲率 $K=-1/L^2$ の時間的偏差では点から出る場が $J=L\sin(\tau/L)J'(0)$ なので、$\tau=\pi L$ に共役点がある。

### 第二変分の境界項

時間的測地線の固有時の二階変分を、始点を空間的な面に置き、終点を固定して計算する。基準曲線の単位接線を $u$、時刻と変分の座標ベクトルを $T,V$ とし、基準曲線では $g(u,V)=0$ にする。$\sqrt{-g(T,T)}$ の微分と $\nabla_VT=\nabla_TV$ により

$$
\delta^2\tau=-\int\bigl[g(D_\tau V,D_\tau V)+g(\nabla_V\nabla_TV,u)\bigr]\dd\tau
$$

となる。平方根からもう一つ $-g(D_\tau V,u)^2$ が出るが、$g(V,u)=0$ と $\nabla_u u=0$ により零である。交換子の定義を用い $\nabla_V\nabla_TV=\nabla_T\nabla_VV+R(V,T)V$ と置く。第一項は端点での $-[g(\nabla_VV,u)]$、第二項は曲率の反対称性により $+g(R(V,u)u,V)$ を与える。始点が面上にあり $u$ が法線なので $g(\nabla_VV,u)=-K(V,V)$ である。ここで $K$ は始点で速度 $u$ と同じ向きの単位法線を使う。従って固有時の二階変分は前に用いた指標形式の負であり、指標形式の始点項は $+K(V,V)$ である。始点も一点に固定する場合は $V=0$ でこの項が消える。

### 光的焦点を越えた区間を時間的に変形する

空間的二面 $S$ から直交して出るアフィン光的測地線 $\gamma$ を $[0,L]$ で取る。横断場 $V$ を $g(V,k)=0$、$V(0)\in T S$、$V(L)=0$ とし、また開始点の変分は $S$ 内に保つ。曲線の変分で速度を $T$、変分方向を $V$ とすると、基準曲線で $\partial_s g(T,T)=2g(D_\lambda V,k)=0$ である。$\tfrac12\partial_s^2g(T,T)$ は
$$
Q(\lambda)=g(D_\lambda V,D_\lambda V)-g(R(V,k)k,V)
+\frac{\dd}{\dd\lambda}g(\nabla_VV,k).
$$
最後の項を積分すると、終点は固定なので零、始点は $g(\nabla_VV,k)=-g(V,\nabla_V k)$ から符号を反転して $+g(V,\nabla_V k)$ を与える。よって $\int Q$ は光的な横断指標形式である。

$L$ より手前に焦点があるとき、横断 Jacobi 場を焦点まで取りその後を零でつなぐと、始点の法線変分条件と Jacobi 方程式によりその指標形式は零になる。一方、焦点でその横断微分は零でない。値と微分がともに零なら横断 Jacobi 方程式の一意性から全てが零になるからである。焦点の近くだけで値を持つ $W$ を取り、部分積分で現れる $g(D_\lambda J,W)$ を負に選ぶと $J+\epsilon W$ の指標形式は、小さな正の $\epsilon$ で負になる。区分的な場を滑らかに近似しても負号を保てる。こうして $\bar Q=L^{-1}\int_0^LQ<0$ となる変分がある。

さらに曲線変分の二次の加速度を変えると、$g(\nabla_VV,k)$ に任意の滑らかな関数 $F(\lambda)$ を足せる。ただし端点を保つため $F(0)=F(L)=0$ とする。実際 $k$ に対し $g(Z,k)=1$ となる補助場 $Z$ を選び、曲線を二次の量 $s^2 FZ/2$ だけ動かせばよい。$F(\lambda)=\int_0^\lambda(\bar Q-Q(\lambda'))\dd\lambda'$ とすれば新しい $Q$ は一定の負値 $\bar Q$ となる。コンパクト区間で十分小さな $s\ne0$ を選ぶと変形した曲線のノルムは全域で負であり、未来向きを保つ。従って焦点を過ぎた終点は $I^+(S)$ に入る。この議論は光的曲線の長さの符号だけを比較する場合と異なり、変形後の接線が各点で時間的であることまで示す。

### 平坦時空でも焦点はできる

Minkowski 時空で $t=0,r=R$ の球から未来内向きの光を出す。角方向の分離は $r(\lambda)=R-\lambda$ に比例するので、横断の Jacobi 行列は初期枠で $A=(1-\lambda/R)I_2$ である。膨張率は $\theta=-2/(R-\lambda)$ で、零曲率でも $\theta'=-\theta^2/2$ に従い $\lambda=R$ に発散する。その点でも一束ごとの光線は正則に延長できるため、発散は光線束による座標の退化である。異なる初期角の光線が同一点で会うことと、測地線が時空内で延長不能になることを区別する。

:::example[共役点のある一定負曲率背景]
半径 $L>0$ の反 de Sitter の静的計量を用いる。時間的測地線に対する偏差で $g(U,U)=-1$ の符号を保って計算する。 座標の順序を $(t,r,\theta,\phi)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}-1 - \frac{r^{2}}{L^{2}} & 0 & 0 & 0\\0 & \frac{1}{1 + \frac{r^{2}}{L^{2}}} & 0 & 0\\0 & 0 & r^{2} & 0\\0 & 0 & 0 & r^{2} \sin^{2}{\left(\theta \right)}\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}- \frac{L^{2}}{L^{2} + r^{2}} & 0 & 0 & 0\\0 & \frac{L^{2} + r^{2}}{L^{2}} & 0 & 0\\0 & 0 & \frac{1}{r^{2}} & 0\\0 & 0 & 0 & \frac{1}{r^{2} \sin^{2}{\left(\theta \right)}}\end{matrix}\right]\\
\det g&=- r^{4} \sin^{2}{\left(\theta \right)}
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{t}_{t r}$ | $- \frac{L^{2}}{2 \left(L^{2} + r^{2}\right)}\left(\partial_{r}\left(-1 - \frac{r^{2}}{L^{2}}\right)\right)$ | $\frac{r}{L^{2} + r^{2}}$ |
| $\Gamma^{r}_{t t}$ | $\frac{L^{2} + r^{2}}{2 L^{2}}\left(-\partial_{r}\left(-1 - \frac{r^{2}}{L^{2}}\right)\right)$ | $\frac{r \left(L^{2} + r^{2}\right)}{L^{4}}$ |
| $\Gamma^{r}_{r r}$ | $\frac{L^{2} + r^{2}}{2 L^{2}}\left(\partial_{r}\left(\frac{1}{1 + \frac{r^{2}}{L^{2}}}\right)+\partial_{r}\left(\frac{1}{1 + \frac{r^{2}}{L^{2}}}\right)-\partial_{r}\left(\frac{1}{1 + \frac{r^{2}}{L^{2}}}\right)\right)$ | $- \frac{r}{L^{2} + r^{2}}$ |
| $\Gamma^{r}_{\theta \theta}$ | $\frac{L^{2} + r^{2}}{2 L^{2}}\left(-\partial_{r}\left(r^{2}\right)\right)$ | $- \frac{r \left(L^{2} + r^{2}\right)}{L^{2}}$ |
| $\Gamma^{r}_{\phi \phi}$ | $\frac{L^{2} + r^{2}}{2 L^{2}}\left(-\partial_{r}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $- \frac{r \left(L^{2} + r^{2}\right) \sin^{2}{\left(\theta \right)}}{L^{2}}$ |
| $\Gamma^{\theta}_{r \theta}$ | $\frac{1}{2 r^{2}}\left(\partial_{r}\left(r^{2}\right)\right)$ | $\frac{1}{r}$ |
| $\Gamma^{\theta}_{\phi \phi}$ | $\frac{1}{2 r^{2}}\left(-\partial_{\theta}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $- \frac{\sin{\left(2 \theta \right)}}{2}$ |
| $\Gamma^{\phi}_{r \phi}$ | $\frac{1}{2 r^{2} \sin^{2}{\left(\theta \right)}}\left(\partial_{r}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{r}$ |
| $\Gamma^{\phi}_{\theta \phi}$ | $\frac{1}{2 r^{2} \sin^{2}{\left(\theta \right)}}\left(\partial_{\theta}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{\tan{\left(\theta \right)}}$ |

曲率の計算では :eqref[eq-gr-curvature-components] の微分差と積差を分ける。第二列は $\partial_c\Gamma^a{}_{db}$、第三列は $-\partial_d\Gamma^a{}_{cb}$ であり、第四列は接続の積を符号付きで足したものである。三列を足して最後の列を得る。

| 成分 | 第一の微分 | 第二の微分 | 接続の積 | 結果 |
|---|---|---|---|---|
| $R^{t}_{r t r}$ | $0$ | $\frac{2 r^{2}}{\left(L^{2} + r^{2}\right)^{2}} - \frac{1}{L^{2} + r^{2}}$ | $- \frac{r^{2}}{\left(L^{2} + r^{2}\right)^{2}} - \frac{r^{2}}{\left(L^{2} + r^{2}\right)^{2}}$ | $- \frac{1}{L^{2} + r^{2}}$ |
| $R^{t}_{\theta t \theta}$ | $0$ | $0$ | $- \frac{r^{2}}{L^{2}}$ | $- \frac{r^{2}}{L^{2}}$ |
| $R^{t}_{\phi t \phi}$ | $0$ | $0$ | $- \frac{r^{2} \sin^{2}{\left(\theta \right)}}{L^{2}}$ | $- \frac{r^{2} \sin^{2}{\left(\theta \right)}}{L^{2}}$ |
| $R^{r}_{\theta r \theta}$ | $- \frac{2 r^{2}}{L^{2}} - \frac{L^{2} + r^{2}}{L^{2}}$ | $0$ | $\frac{r^{2}}{L^{2}} + \frac{L^{2} + r^{2}}{L^{2}}$ | $- \frac{r^{2}}{L^{2}}$ |
| $R^{r}_{\phi r \phi}$ | $- \frac{2 r^{2} \sin^{2}{\left(\theta \right)}}{L^{2}} - \frac{\left(L^{2} + r^{2}\right) \sin^{2}{\left(\theta \right)}}{L^{2}}$ | $0$ | $\frac{r^{2} \sin^{2}{\left(\theta \right)}}{L^{2}} + \frac{\left(L^{2} + r^{2}\right) \sin^{2}{\left(\theta \right)}}{L^{2}}$ | $- \frac{r^{2} \sin^{2}{\left(\theta \right)}}{L^{2}}$ |
| $R^{\theta}_{\phi \theta \phi}$ | $- \cos{\left(2 \theta \right)}$ | $0$ | $\frac{\sin{\left(2 \theta \right)}}{2 \tan{\left(\theta \right)}} - \frac{\left(L^{2} + r^{2}\right) \sin^{2}{\left(\theta \right)}}{L^{2}}$ | $- \frac{r^{2} \sin^{2}{\left(\theta \right)}}{L^{2}}$ |

第一添字を計量で下ろすと $R_{abcd}=g_{ak}R^k{}_{bcd}$ となる。対の反対称性と交換対称性を用いて、非零の独立な組を下表に示す。

| 成分 | 結果 |
|---|---|
| $R_{t r t r}$ | $\frac{1}{L^{2}}$ |
| $R_{t \theta t \theta}$ | $\frac{r^{2} \left(L^{2} + r^{2}\right)}{L^{4}}$ |
| $R_{t \phi t \phi}$ | $\frac{r^{2} \left(L^{2} + r^{2}\right) \sin^{2}{\left(\theta \right)}}{L^{4}}$ |
| $R_{r \theta r \theta}$ | $- \frac{r^{2}}{L^{2} + r^{2}}$ |
| $R_{r \phi r \phi}$ | $- \frac{r^{2} \sin^{2}{\left(\theta \right)}}{L^{2} + r^{2}}$ |
| $R_{\theta \phi \theta \phi}$ | $- \frac{r^{4} \sin^{2}{\left(\theta \right)}}{L^{2}}$ |

Ricci テンソルは第一・第三添字を縮約して得る。曲率が残る場合も相殺が起こるので、各項を落とさずに加える。

| 成分 | 縮約の項 | 代入 | 結果 |
|---|---|---|---|
| $R_{t t}$ | $R^{r}_{t r t}+R^{\theta}_{t \theta t}+R^{\phi}_{t \phi t}$ | $\frac{L^{2} + r^{2}}{L^{4}} + \frac{L^{2} + r^{2}}{L^{4}} + \frac{L^{2} + r^{2}}{L^{4}}$ | $\frac{3 \left(L^{2} + r^{2}\right)}{L^{4}}$ |
| $R_{r r}$ | $R^{t}_{r t r}+R^{\theta}_{r \theta r}+R^{\phi}_{r \phi r}$ | $- \frac{1}{L^{2} + r^{2}} - \frac{1}{L^{2} + r^{2}} - \frac{1}{L^{2} + r^{2}}$ | $- \frac{3}{L^{2} + r^{2}}$ |
| $R_{\theta \theta}$ | $R^{t}_{\theta t \theta}+R^{r}_{\theta r \theta}+R^{\phi}_{\theta \phi \theta}$ | $- \frac{r^{2}}{L^{2}} - \frac{r^{2}}{L^{2}} - \frac{r^{2}}{L^{2}}$ | $- \frac{3 r^{2}}{L^{2}}$ |
| $R_{\phi \phi}$ | $R^{t}_{\phi t \phi}+R^{r}_{\phi r \phi}+R^{\theta}_{\phi \theta \phi}$ | $- \frac{r^{2} \sin^{2}{\left(\theta \right)}}{L^{2}} - \frac{r^{2} \sin^{2}{\left(\theta \right)}}{L^{2}} - \frac{r^{2} \sin^{2}{\left(\theta \right)}}{L^{2}}$ | $- \frac{3 r^{2} \sin^{2}{\left(\theta \right)}}{L^{2}}$ |

逆計量と縮約すると Ricci スカラーとなる。さらに $G_{ab}=R_{ab}-g_{ab}R/2$ を成分ごとに計算する。

$$
\begin{aligned}
R&=\sum_{a,b}g^{ab}R_{ab}\\
&=- \frac{3}{L^{2}} - \frac{3}{L^{2}} - \frac{3}{L^{2}} - \frac{3}{L^{2}}\\
&=- \frac{12}{L^{2}}
\end{aligned}
$$

非零の Einstein テンソル成分は次のとおりである。

| 成分 | Ricci 項 | 跡の項 | 差 |
|---|---|---|---|
| $G_{t t}$ | $\frac{3 \left(L^{2} + r^{2}\right)}{L^{4}}$ | $- \frac{6 \left(-1 - \frac{r^{2}}{L^{2}}\right)}{L^{2}}$ | $- \frac{3 \left(L^{2} + r^{2}\right)}{L^{4}}$ |
| $G_{r r}$ | $- \frac{3}{L^{2} + r^{2}}$ | $- \frac{6}{L^{2} \left(1 + \frac{r^{2}}{L^{2}}\right)}$ | $\frac{3}{L^{2} + r^{2}}$ |
| $G_{\theta \theta}$ | $- \frac{3 r^{2}}{L^{2}}$ | $- \frac{6 r^{2}}{L^{2}}$ | $\frac{3 r^{2}}{L^{2}}$ |
| $G_{\phi \phi}$ | $- \frac{3 r^{2} \sin^{2}{\left(\theta \right)}}{L^{2}}$ | $- \frac{6 r^{2} \sin^{2}{\left(\theta \right)}}{L^{2}}$ | $\frac{3 r^{2} \sin^{2}{\left(\theta \right)}}{L^{2}}$ |

成分表では $R_{abcd}=-(g_{ac}g_{bd}-g_{ad}g_{bc})/L^2$ が成り立つ。単位時間的速度 $U$ に直交する Jacobi 場では $R(U,J)U=-J/L^2$ なので、始点で零の場は $J=L\sin(\tau/L)W$ である。$W$ は平行輸送した空間方向である。$\tau=\pi L$ で再び零になるが、計量の曲率が発散するわけではない。例えば中心の時間的測地線では $t$ が固有時であり、曲線はその時刻を越えて正則に続く。
:::

## Penrose の特異点定理

:::theorem[Penrose の特異点定理]{#thm-gr-penrose-singularity}
連結で時間向きのある大域的双曲時空に、非コンパクトな Cauchy 面とコンパクトで境界のない滑らかな未来捕捉面 $S$ があるとする。$S$ はコンパクトな空間的二次元面で、両方の未来光的法線の膨張率が負である。全光的 $k$ に $R_{ab}k^ak^b\ge0$ が成立するなら、未来光的測地線が少なくとも一本不完備である。
:::

:::proof
未来の光的測地線が全て完備と仮定する。滑らかな未来時間的単位場 $u$ を使い $-k\cdot u=1$ で $S$ の法線を規格化する。法線の組は $S$ 上のコンパクトな二重被覆であり、膨張率は連続で負なので $\theta_0\le-\delta<0$ の一様な定数 $\delta$ がある。未来境界 $E^+(S)=J^+(S)\setminus I^+(S)$ は $S$ 自身との交点を除き、$S$ から直交して出る光的測地線で、途中に焦点のないものの上にある。これは因果境界の生成子定理であり、$J^+(S)$ の閉性と合わせて :cite[hawking-ellis] を用いる。:eqref[eq-gr-raychaudhuri-null] による焦点はアフィン径数 $2/\delta$ 以下なので、$E^+(S)$ はコンパクトな初期法線と閉区間 $[0,2/\delta]$ の指数写像像の閉部分集合となる。従ってコンパクトである。

一方、$E^+(S)=\partial I^+(S)$ は非空である。実際、連続な Cauchy 時間関数を取り、コンパクトな $S$ 上での最小値より過去には $I^+(S)$ がなく、$S$ の点から時間的未来へ進めば $I^+(S)$ に入る。$J^+(S)$ の閉性と、各因果的未来点から少し未来へ進むと $I^+(S)$ に入ることから $\overline{I^+(S)}=J^+(S)$ である。境界上の点のごく近くの同じ時間的流線で、未来側は $I^+(S)$ の中、過去側は閉包の外にある。上下の点を少しずらしてもこの関係が保たれるので、局所のフローボックスで近くの各流線が境界と一度ずつ交わる。その交点を挟む上下の時刻を任意に近づけられるため、境界は連続なグラフである。

完全な時間的ベクトル場の流線に沿い Cauchy 面へ射影すると、非時間的関連性により単射になる。グラフの局所性から像は開であり、コンパクト性から閉でもある。Cauchy 面は時空の連結性により連結であり、非空の像がその全体になる。しかし面は非コンパクトと仮定した。矛盾により未来光的完備性は成立しない。
:::

例えば球対称な崩壊時空で未来内部の面積半径を両方向の未来光が減らすと、球の膨張率は $\theta=2k(r)/r<0$ となる。定理の結論は有限アフィン径数で延長できない測地線の存在であり、どの曲率不変量がどう発散するかは別の問題である。

:::example[球対称の捕捉面を座標で判定する]
$\dd s^2=-f\dd v^2+2\dd v\,\dd r+r^2\dd\Omega^2$、$f=1-2m(v)/r$ の入射計量を用いる。定数 $r,v$ の球面の二つの未来向き光的法線は $k_+=\partial_v+(f/2)\partial_r$、$k_-=-\partial_r$ と取れる。$g(k_+,k_-)=-1$ である。 座標の順序を $(v,r,\theta,\phi)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}-1 + \frac{2 m{\left(v \right)}}{r} & 1 & 0 & 0\\1 & 0 & 0 & 0\\0 & 0 & r^{2} & 0\\0 & 0 & 0 & r^{2} \sin^{2}{\left(\theta \right)}\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}0 & 1 & 0 & 0\\1 & \frac{r - 2 m{\left(v \right)}}{r} & 0 & 0\\0 & 0 & \frac{1}{r^{2}} & 0\\0 & 0 & 0 & \frac{1}{r^{2} \sin^{2}{\left(\theta \right)}}\end{matrix}\right]\\
\det g&=- r^{4} \sin^{2}{\left(\theta \right)}
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{v}_{v v}$ | $\frac{1}{2}\left(-\partial_{r}\left(-1 + \frac{2 m{\left(v \right)}}{r}\right)\right)$ | $\frac{m{\left(v \right)}}{r^{2}}$ |
| $\Gamma^{v}_{\theta \theta}$ | $\frac{1}{2}\left(-\partial_{r}\left(r^{2}\right)\right)$ | $- r$ |
| $\Gamma^{v}_{\phi \phi}$ | $\frac{1}{2}\left(-\partial_{r}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $- r \sin^{2}{\left(\theta \right)}$ |
| $\Gamma^{r}_{v v}$ | $\frac{1}{2}\left(\partial_{v}\left(-1 + \frac{2 m{\left(v \right)}}{r}\right)+\partial_{v}\left(-1 + \frac{2 m{\left(v \right)}}{r}\right)-\partial_{v}\left(-1 + \frac{2 m{\left(v \right)}}{r}\right)\right)+\frac{r - 2 m{\left(v \right)}}{2 r}\left(-\partial_{r}\left(-1 + \frac{2 m{\left(v \right)}}{r}\right)\right)$ | $\frac{r^{2} m'{\left(v \right)} + r m{\left(v \right)} - 2 m^{2}{\left(v \right)}}{r^{3}}$ |
| $\Gamma^{r}_{v r}$ | $\frac{1}{2}\left(\partial_{r}\left(-1 + \frac{2 m{\left(v \right)}}{r}\right)\right)$ | $- \frac{m{\left(v \right)}}{r^{2}}$ |
| $\Gamma^{r}_{\theta \theta}$ | $\frac{r - 2 m{\left(v \right)}}{2 r}\left(-\partial_{r}\left(r^{2}\right)\right)$ | $- r + 2 m{\left(v \right)}$ |
| $\Gamma^{r}_{\phi \phi}$ | $\frac{r - 2 m{\left(v \right)}}{2 r}\left(-\partial_{r}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $- \left(r - 2 m{\left(v \right)}\right) \sin^{2}{\left(\theta \right)}$ |
| $\Gamma^{\theta}_{r \theta}$ | $\frac{1}{2 r^{2}}\left(\partial_{r}\left(r^{2}\right)\right)$ | $\frac{1}{r}$ |
| $\Gamma^{\theta}_{\phi \phi}$ | $\frac{1}{2 r^{2}}\left(-\partial_{\theta}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $- \frac{\sin{\left(2 \theta \right)}}{2}$ |
| $\Gamma^{\phi}_{r \phi}$ | $\frac{1}{2 r^{2} \sin^{2}{\left(\theta \right)}}\left(\partial_{r}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{r}$ |
| $\Gamma^{\phi}_{\theta \phi}$ | $\frac{1}{2 r^{2} \sin^{2}{\left(\theta \right)}}\left(\partial_{\theta}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{\tan{\left(\theta \right)}}$ |

曲率の計算では :eqref[eq-gr-curvature-components] の微分差と積差を分ける。第二列は $\partial_c\Gamma^a{}_{db}$、第三列は $-\partial_d\Gamma^a{}_{cb}$ であり、第四列は接続の積を符号付きで足したものである。三列を足して最後の列を得る。

| 成分 | 第一の微分 | 第二の微分 | 接続の積 | 結果 |
|---|---|---|---|---|
| $R^{v}_{\theta v \theta}$ | $0$ | $0$ | $- \frac{m{\left(v \right)}}{r}$ | $- \frac{m{\left(v \right)}}{r}$ |
| $R^{v}_{\phi v \phi}$ | $0$ | $0$ | $- \frac{m{\left(v \right)} \sin^{2}{\left(\theta \right)}}{r}$ | $- \frac{m{\left(v \right)} \sin^{2}{\left(\theta \right)}}{r}$ |
| $R^{r}_{r v r}$ | $0$ | $- \frac{2 m{\left(v \right)}}{r^{3}}$ | $0$ | $- \frac{2 m{\left(v \right)}}{r^{3}}$ |
| $R^{r}_{\theta v \theta}$ | $2 m'{\left(v \right)}$ | $0$ | $- \frac{\left(- r + 2 m{\left(v \right)}\right) m{\left(v \right)}}{r^{2}} - \frac{r^{2} m'{\left(v \right)} + r m{\left(v \right)} - 2 m^{2}{\left(v \right)}}{r^{2}}$ | $m'{\left(v \right)}$ |
| $R^{r}_{\theta r \theta}$ | $-1$ | $0$ | $\frac{r - 2 m{\left(v \right)}}{r} + \frac{m{\left(v \right)}}{r}$ | $- \frac{m{\left(v \right)}}{r}$ |
| $R^{r}_{\phi v \phi}$ | $2 \sin^{2}{\left(\theta \right)} m'{\left(v \right)}$ | $0$ | $\frac{\left(r - 2 m{\left(v \right)}\right) m{\left(v \right)} \sin^{2}{\left(\theta \right)}}{r^{2}} - \frac{\left(r^{2} m'{\left(v \right)} + r m{\left(v \right)} - 2 m^{2}{\left(v \right)}\right) \sin^{2}{\left(\theta \right)}}{r^{2}}$ | $\sin^{2}{\left(\theta \right)} m'{\left(v \right)}$ |
| $R^{r}_{\phi r \phi}$ | $- \sin^{2}{\left(\theta \right)}$ | $0$ | $\frac{\left(r - 2 m{\left(v \right)}\right) \sin^{2}{\left(\theta \right)}}{r} + \frac{m{\left(v \right)} \sin^{2}{\left(\theta \right)}}{r}$ | $- \frac{m{\left(v \right)} \sin^{2}{\left(\theta \right)}}{r}$ |
| $R^{\theta}_{\phi \theta \phi}$ | $- \cos{\left(2 \theta \right)}$ | $0$ | $\frac{\sin{\left(2 \theta \right)}}{2 \tan{\left(\theta \right)}} - \frac{\left(r - 2 m{\left(v \right)}\right) \sin^{2}{\left(\theta \right)}}{r}$ | $\frac{2 m{\left(v \right)} \sin^{2}{\left(\theta \right)}}{r}$ |

第一添字を計量で下ろすと $R_{abcd}=g_{ak}R^k{}_{bcd}$ となる。対の反対称性と交換対称性を用いて、非零の独立な組を下表に示す。

| 成分 | 結果 |
|---|---|
| $R_{v r v r}$ | $- \frac{2 m{\left(v \right)}}{r^{3}}$ |
| $R_{v \theta v \theta}$ | $\frac{r^{2} m'{\left(v \right)} + r m{\left(v \right)} - 2 m^{2}{\left(v \right)}}{r^{2}}$ |
| $R_{v \theta r \theta}$ | $- \frac{m{\left(v \right)}}{r}$ |
| $R_{v \phi v \phi}$ | $\frac{\left(r^{2} m'{\left(v \right)} + r m{\left(v \right)} - 2 m^{2}{\left(v \right)}\right) \sin^{2}{\left(\theta \right)}}{r^{2}}$ |
| $R_{v \phi r \phi}$ | $- \frac{m{\left(v \right)} \sin^{2}{\left(\theta \right)}}{r}$ |
| $R_{\theta \phi \theta \phi}$ | $2 r m{\left(v \right)} \sin^{2}{\left(\theta \right)}$ |

Ricci テンソルは第一・第三添字を縮約して得る。曲率が残る場合も相殺が起こるので、各項を落とさずに加える。

| 成分 | 縮約の項 | 代入 | 結果 |
|---|---|---|---|
| $R_{v v}$ | $R^{r}_{v r v}+R^{\theta}_{v \theta v}+R^{\phi}_{v \phi v}$ | $- \frac{2 \left(r - 2 m{\left(v \right)}\right) m{\left(v \right)}}{r^{4}} + \frac{r^{2} m'{\left(v \right)} + r m{\left(v \right)} - 2 m^{2}{\left(v \right)}}{r^{4}} + \frac{r^{2} m'{\left(v \right)} + r m{\left(v \right)} - 2 m^{2}{\left(v \right)}}{r^{4}}$ | $\frac{2 m'{\left(v \right)}}{r^{2}}$ |
| $R_{v r}$ | $R^{v}_{v v r}+R^{\theta}_{v \theta r}+R^{\phi}_{v \phi r}$ | $- \frac{m{\left(v \right)}}{r^{3}} - \frac{m{\left(v \right)}}{r^{3}} + \frac{2 m{\left(v \right)}}{r^{3}}$ | $0$ |
| $R_{\theta \theta}$ | $R^{v}_{\theta v \theta}+R^{r}_{\theta r \theta}+R^{\phi}_{\theta \phi \theta}$ | $- \frac{m{\left(v \right)}}{r} - \frac{m{\left(v \right)}}{r} + \frac{2 m{\left(v \right)}}{r}$ | $0$ |
| $R_{\phi \phi}$ | $R^{v}_{\phi v \phi}+R^{r}_{\phi r \phi}+R^{\theta}_{\phi \theta \phi}$ | $- \frac{m{\left(v \right)} \sin^{2}{\left(\theta \right)}}{r} - \frac{m{\left(v \right)} \sin^{2}{\left(\theta \right)}}{r} + \frac{2 m{\left(v \right)} \sin^{2}{\left(\theta \right)}}{r}$ | $0$ |

逆計量と縮約すると Ricci スカラーとなる。さらに $G_{ab}=R_{ab}-g_{ab}R/2$ を成分ごとに計算する。

$$
\begin{aligned}
R&=\sum_{a,b}g^{ab}R_{ab}\\
&=0\\
&=0
\end{aligned}
$$

非零の Einstein テンソル成分は次のとおりである。

| 成分 | Ricci 項 | 跡の項 | 差 |
|---|---|---|---|
| $G_{v v}$ | $\frac{2 m'{\left(v \right)}}{r^{2}}$ | $0$ | $\frac{2 m'{\left(v \right)}}{r^{2}}$ |

面積要素は $r^2\sin\theta\,\dd\theta\,\dd\phi$ だから、各光束の膨張率は $2k_\pm(r)/r$ である。従って $\theta_+=f/r,\theta_-=-2/r$ となり、$0<r<2m(v)$ の球は未来捕捉面である。$m'\ge0$ なら Ricci による光的収束条件は任意の光的 $k$ に対して $2m'(k^v)^2/r^2\ge0$ である。定理を適用するには、さらに非コンパクトな Cauchy 面の存在を全時空で調べる。
:::

### 未来ではなく過去に捕捉される面

球対称な平坦 Robertson–Walker の膨張枝で $r=a(t)\chi$ を面積半径とする。共動の直交枠で未来光的方向を $k_\pm=\partial_t\pm a^{-1}\partial_\chi$ と取ると $k_\pm(r)=Hr\pm1$ であり、
$$
\theta_\pm=2(Hr\pm1)/r
$$
となる。$H>0,r>1/H$ では未来の両膨張率が正なので、過去向きに両法線を反転すると両方負になる。従ってその球は過去捕捉面である。零収束条件と非コンパクトな Cauchy 面の仮定も成り立てば、Penrose の定理を時間反転して過去光的不完備性を得る。例えば $a=t^{2/3}$、$t>0$ の塵模型では $H=2/(3t)$ で、半径 $r=3t$ なら $Hr=2$ である。この幾何では動径零測地線の $\dd\lambda/\dd t$ は $a(t)$ に比例するので、$t\downarrow0$ までの積分 $\int_0^{t_0} t^{2/3}\dd t=3t_0^{5/3}/5$ が有限であることも直接確かめられる。空間切片に正の曲率を持つ閉じた模型なら非コンパクト性が失われるため、同じ面の符号だけでこの形の定理を適用しない。

## Hawking の特異点定理

:::theorem[Hawking の宇宙論的特異点定理の一形]{#thm-gr-hawking-singularity}
大域的双曲時空に滑らかな空間的 Cauchy 面 $\Sigma$ があり、その未来法線の膨張率が $\theta\ge\theta_*>0$ とする。全時間的ベクトル $u$ に $R_{ab}u^au^b\ge0$ が成立すれば、時空は過去時間的に不完備である。
:::

:::proof
過去向きの法線では初期膨張率が $-\theta\le-\theta_*$ となる。:eqref[eq-gr-raychaudhuri-timelike] により、その法線測地線が十分続けば $3/\theta_*$ 以内で焦点を生じる。$p$ を $\Sigma$ の時間的過去の点とする。Cauchy 面への最大化定理により $\Sigma$ と $p$ を結ぶ最長時間的曲線があり、自由な面上の端点の第一変分により面に直交する。途中に焦点があると長さの最大性を失うので、この長さは $3/\theta_*$ 以下である。従って $\Sigma$ から出る任意の過去向き時間的曲線長も同じ上限を持つ。完備なら面からの単位時間的測地線をその固有時より長く延長でき、矛盾する。
:::

例えば平坦な物質宇宙の $t=t_0$ で $\theta=3H=2/t_0$ なので上限は $3t_0/2$ である。実際の共動線の過去までの固有時は $t_0$ であり、上限が厳密な寿命を与えるわけではない。正の $\Lambda$ は $R_{ab}u^au^b$ に負に寄与し、de Sitter 時空にこの収束仮定を課せない。

### 膨張面から過去への上限を具体化する

平坦な塵の膨張モデル $a(t)\propto t^{2/3}$、$t>0$ では、$t=t_0$ の共動法線の膨張率は $3H=2/t_0$ である。加速度は零で $\ddot a/a=-2/(9t^2)$ なので $R_{00}=2/(3t^2)>0$ である。完全流体の跡反転で $R_{ab}u^au^b=8\pi G\rho[(u\cdot n)^2-1/2]>0$ が任意の単位時間的 $u$ に成り立つ。従って定理の上限は $3t_0/2$ である。実際に共動法線の固有時は $t$ そのものだから過去端まで $t_0$ であり、上限より短い。Raychaudhuri の不等式で Ricci の非負項を落としているため、一般に等号は要求されない。

### 正の宇宙項が定理の仮定を外す例

:::example[指数膨張の時間面で収束を検算する]
平坦パッチの de Sitter で $a=e^{Ht}$、$H>0$ とする。共動する未来法線の膨張率は $3H$ で一定だが、$R_{ab}=3H^2g_{ab}$ により $R_{ab}u^au^b=-3H^2$ である。Raychaudhuri の式では $-\theta^2/3=-3H^2$ と $-R_{ab}u^au^b=+3H^2$ が相殺して $\dot\theta=0$ となる。従って膨張率に正の下限があっても、時間的収束の条件がなければ Hawking の固有時上限は従わない。
:::

### 収束条件は物質と宇宙項のどこに現れるか

四次元で完全流体 $T_{ab}=(\rho+p)n_an_b+pg_{ab}$ を取り、未来単位時間的 $u$ について $\gamma=-u\cdot n\ge1$ とする。跡は $T=-\rho+3p$ なので、:eqref[eq-gr-trace-einstein] の縮約は
$$
R_{uu}=8\pi G\left[(\rho+p)(\gamma^2-1)+\frac{\rho+3p}{2}\right]-\Lambda
$$
である。従って全ての $u$ に時間的収束を課すには $\rho+p\ge0$ と $4\pi G(\rho+3p)\ge\Lambda$ が必要十分である。$\rho+p<0$ なら $\gamma$ を増やすと負になる。$\rho+p\ge0$ なら最小値は $\gamma=1$ にある。光的な場合は $R_{kk}=8\pi G(\rho+p)(k\cdot n)^2$ なので、宇宙項と圧力の跡は消える。例えば $p=-\rho$、$\rho>0$ では、$\Lambda=0$ でも光的収束は等号で成立し、時間的収束は成立しない。

Hawking の証明で $\Sigma$ から $p$ への最大曲線を過去向きに取ると、端点の第一変分が $g(u,V(0))$ に比例する。$V(0)$ は $\Sigma$ の接方向を自由に動けるので、曲線は直交する。もし長さが $3/\theta_*$ より大きければ、その手前で焦点があるため第二変分の議論と矛盾する。従って初期面からの最大固有時は一様に有界である。途中で曲がる観測者もこの上限を超えられないが、任意の観測者が同じ固有時で終わるという主張は従わない。

:::example[空間がコンパクトな膨張型の模型]
$t>0$ で $a(t)=t^{1/2}$ とし、平坦空間を三次元トーラスに同一視する。$t$ を使う各面の間では $|\dd\vb x/\dd t|\le1/a(t)$ であり、$0<t_1<t_2$ の内部で因果曲線が時間を止めれば位置も収束して延長できる。従って全ての時刻切片が Cauchy 面になる。$\Lambda=0$ の放射流体では $p=\rho/3$、$\rho=3/(32\pi Gt^2)$ なので時間的収束が成立する。$t=t_0$ で未来法線の膨張率は $3/(2t_0)$ であり、過去への上限は $2t_0$ となる。共動の測地線は $t=0$ まで固有時 $t_0$ であり、反対に未来の固有時は無限に続く。空間の非コンパクト性を欠くこの場合にも、時間的な過去不完備性の定理は適用できる。
:::

### Hawking–Penrose の定理との仮定の違い

ここまでの二つの定理では、大域的双曲性を仮定している。:term[Hawking–Penrose の特異点定理]{en="Hawking–Penrose singularity theorem" yomi="ほーきんぐぺんろーずのとくいてんていり"} の一形は、閉じた時間的曲線がなく、すべての因果的 $k$ に $R_{kk}\ge0$ が成り立ち、すべての延長不能な因果測地線が曲率の一般性条件を満たすときに使える。最後の条件は、各測地線上の少なくとも一点で $k_{[a}R_{b]cd[e}k_{f]}k^c k^d\ne0$ となることである。ここで角括弧は $(a,b)$ と $(e,f)$ をそれぞれ反対称化する。

さらに、コンパクトで縁のない非時間的関連集合がある、閉じた捕捉面がある、または一点から出る一方の時間向きのすべての光線がそれぞれ途中で負の膨張率を得る、のいずれかを仮定する。このとき不完備な因果測地線が存在する。ここでは Penrose の形のように不完備なものを未来光的と指定したり、Hawking の形のように面からの固有時に一様上限を与えたりはしない。証明は因果測地線全体の共役点の存在と、完全であれば仮定から得られる最大化測地線との矛盾を使う。Cauchy 面への射影では扱えない因果境界の議論が必要なので :cite[hawking-ellis] を参照する。

例えば Minkowski 時空の空間三方向を周期同一視すれば、コンパクトで縁のない時刻面があり、閉じた時間的曲線はない。また曲率零なので収束条件は等号で成立する。しかし全ての因果測地線で曲率の一般性条件が破れ、測地線は平坦な被覆で直線として完全に延長できる。従って一般性条件を単に省略することはできない。

## 宇宙検閲仮説

:term[宇宙検閲仮説]{en="cosmic censorship conjecture" yomi="うちゅうけんえつかせつ"} には異なる形がある。弱い形は適切な一般的な漸近平坦初期データの崩壊で、遠方の観測者に特異な領域が直接むき出しにならないという主張である。強い形は最大大域双曲発展が、指定した正則性の時空として一般には延長できないという主張である。「一般的」と延長の正則性の指定が定式化の一部になる。例えば理想的な Reissner–Nordström 解は内側に Cauchy 地平線をもち解析拡張できるが、特殊な定常解なので、それだけで一般的な初期データに対する主張を否定できない。Penrose の定理は不完備性を証明するが、大域的にどの観測者が何を受信できるかまでは結論しない。内部地平線の不安定性を証明するには非線形発展と延長の正則性を解析する必要がある。

### 完備性・拡張可能性・可視性

時間的測地線 $\gamma(s)$ が有限な上端 $s_*$ までしか定義されないとする。もしある滑らかな等長拡張において $\gamma$ と速度が正則な点に近づくなら、測地線方程式の初期値問題をその点で解いて続けられる。例えば $t>0$ の Minkowski 部分だけを取ると、静止測地線は過去に有限の固有時で終わるが、$t=0$ を越えて Minkowski 時空へ延長できる。これは不完備性の主張だけから曲率の発散や延長不能性が従わないことを示す。強い宇宙検閲で拡張の正則性を指定する理由は、連続な計量として続けられるかと、古典的な場の方程式を満たす二回微分可能な計量として続けられるかで結論が異なり得るためである。

### 特異領域が見えるかを光で判定する

漸近平坦な時空で未来零無限遠 $\mathscr I^+$ が定義できる場合、外まで信号を届けられる領域は $J^-(\mathscr I^+)$ の物理時空内の部分である。事象地平線はその物理時空内の境界になる。弱い宇宙検閲の典型的な定式化では、適切な孤立系の初期データに対し、未来零無限遠の完全性やその過去での予測可能性を要求する。ここで使う漸近平坦性と共形境界の定義は [漸近的平坦性と質量](/physics/general-relativity/asymptotic-flatness/) を参照する。

例えば形式的に負質量 $M<0$ の Schwarzschild 解を考える。$f=1+2G|M|/r>0$ は全ての $r>0$ で正である。動径の光的測地線では定常性より $E=f\dot t$ が一定で、零条件により $\dot r=E>0$ の出射枝がある。$r=0$ から有限アフィン径数しか離れていない点から、$r\to\infty$ までこの枝を延長できる。出射の遅れ時刻は $u=t-\int\dd r/f$ で一定であり、$r=0$ 付近のこの積分は収束する。従って有限の遅れ時刻で原点に近い領域から無限遠へ光が届く。一方、外部の曲率不変量 $48(GM)^2/r^6$ は $r\downarrow0$ で発散する。初期切片は原点で完全でなく、そこから正則な崩壊が生じた模型ではない。従ってこれは弱い宇宙検閲で可視性を区別する例であり、正則な初期値を要求した主張の反例ではない。

### 初期データの最大発展と Cauchy 地平線

強い宇宙検閲でいう最大発展は、与えた初期面を Cauchy 面として含む発展のうちの最大のものである。例えば $t>0$、$|\vb x|<t$ の Minkowski 未来円錐で、$\tau=\sqrt{t^2-|\vb x|^2}$ は増加する時間関数である。座標を $t=\tau\cosh\chi$、$r=\tau\sinh\chi$ とすれば $g=-\dd\tau^2+\tau^2h_{\mathbb H^3}$ であり、$\eta=\log\tau$ を使うと $g=\tau^2(-\dd\eta^2+h_{\mathbb H^3})$ となる。双曲空間が完全で $\eta$ が全実数を動くため、$\tau=\tau_0$ は円錐内で Cauchy 面になる。しかし Minkowski 全体では円錐の外を通る延長不能な光線がこの面と交わらない。延長した計量が平坦でも、円錐内の Cauchy 面の予測範囲が Minkowski 全体へ広がるわけではない。この例は、拡張を考える際に初期面の Cauchy 性を課したかどうかを区別するためのものである。

## 演習問題

::::exercise{level=basic}
Minkowski 平面で $(0,0)$ から $(3,4)$ に因果曲線で行けるか。

:::solution
空間距離 $4$ が時間差 $3$ を超え、光速の不等式に反するので行けない。
:::
::::

::::exercise{level=basic}
半径 $r$ の球に沿う光束で $k(r)=-1$ の膨張率を求めよ。

:::solution
面積 $4\pi r^2$ の対数微分は $2k(r)/r=-2/r$ である。
:::
::::

::::exercise{level=basic}
初期膨張率 $-4$ のねじれない光束に Ricci 収束条件を課したとき、焦点形成の上限を求めよ。

:::solution
光束の横方向は二次元なので $2/\abs{-4}=1/2$ である。ここでは測地線がそこまで正則に延長するとする。
:::
::::

::::exercise{level=basic}
$K=-1/L^2$ の例で最初の時間的共役点を求めよ。

:::solution
偏差の式は $J''=-J/L^2$ であり、$J(0)=0$ なら $J=L\sin(\tau/L)J'(0)$ なので最初の正の零は $\pi L$ である。
:::
::::

::::exercise{level=advanced}
$\dot\theta=-\theta^2/3$、$\theta(0)=-3/T$ を解け。

:::solution
$(1/\theta)^{\cdot}=1/3$ なので $1/\theta=(t-T)/3$、従って $\theta=-3/(T-t)$ である。$t=T$ に負へ発散する。
:::
::::

::::exercise{level=advanced}
$\Lambda$ が零収束条件の左右を変えない理由を説明せよ。

:::solution
:eqref[eq-gr-trace-einstein] を $k^ak^b$ と縮約する。$g(k,k)=0$ なので跡項と $\Lambda$ 項が消え、$R_{ab}k^ak^b=8\pi G T_{ab}k^ak^b$ となる。
:::
::::

::::exercise{level=advanced}
時間向きのある時空において閉じた時間的曲線の点で強因果性が成り立たないことを示せ。

:::solution
その曲線は点の十分小さい近傍から一度離れ、同じ点へ戻る。従ってどれだけ小さな因果曲線の再入を禁じる近傍を選んでも、その閉曲線が再入する。
:::
::::

::::exercise{level=advanced}
de Sitter 時空の $\Lambda>0$ の真空で、Hawking の定理の時間的曲率条件を判定せよ。

:::solution
$R_{ab}=\Lambda g_{ab}$ なので単位時間的 $u$ に対し $R_{uu}=-\Lambda<0$ となる。時間的収束条件は成立しない。
:::
::::

::::exercise{level=advanced}
平坦時空で三方向の Jacobi 行列を $A(\tau)=\diag(1-b\tau,1,1)$、$b>0$ とする。$\theta$、$\sigma_{ij}\sigma^{ij}$ と焦点時刻を求め、Raychaudhuri 方程式と比較せよ。

:::solution
$B=\dot A A^{-1}=\diag(-b/(1-b\tau),0,0)$ なので $\theta=-b/(1-b\tau)$ である。$B$ は対称で、せん断は $\diag(2\theta/3,-\theta/3,-\theta/3)$ だから二乗は $2\theta^2/3$ となる。渦度と曲率は零なので右辺は $-\theta^2/3-2\theta^2/3=-\theta^2$ であり、$\dot\theta=-b^2/(1-b\tau)^2$ と一致する。焦点は $A$ の第一成分が消える $\tau=1/b$ にある。せん断を落とした不等式から得る上限 $3/b$ より早い。曲率が零でも、方向ごとに異なる初期速度勾配がせん断を作る。
:::
::::

## 確認テスト

::::quiz
Penrose の定理で結論される命題を選べ。

- [x] 未来光的に不完備な測地線が存在する。
- [ ] 全ての未来時間的測地線が有限の固有時で終わる。
- [ ] 同じ場所で全ての曲率成分が発散する。

:::explanation
捕捉面による焦点形成と非コンパクト Cauchy 面を組み合わせると未来の全光的完備性が矛盾する。不完備な対象は少なくとも一本の光的測地線であり、曲率発散は仮定にも結論にも含めていない。
:::
::::

::::quiz{answer="2" tol="1e-9"}
光的 Raychaudhuri 方程式の膨張率二乗項を $-\theta^2/d$ としたとき、四次元の $d$ を答えよ。

:::explanation
光的方向の直交空間は三次元だが接線自身を含む。それを商に取ると横断空間は二次元なので $d=2$ である。
:::
::::

::::quiz
Hawking の過去不完備性の上の形で必要な条件を選べ。

- [x] 初期面の未来向き膨張率に一様な正の下限がある。
- [x] 時間的収束条件が成り立つ。
- [ ] Cauchy 面が非コンパクトであることが必須である。

:::explanation
過去向き法線の初期膨張率は一様に負になる。時間的収束条件と合わせて最長曲線の長さを制限する。ここでは Cauchy 面の非コンパクト性を使わない。
:::
::::

::::quiz{answer="3" tol="1e-9"}
未来向きの初期膨張率の一様下限 $\theta_*=1$ に対する Hawking の定理の固有時上限を答えよ。

:::explanation
過去向きの膨張率は $-1$ 以下であり、時間的束の横次元が三なので $3/\theta_*=3$ である。
:::
::::

::::quiz{answer="-2" tol="1e-7"}
平坦時空で半径 $R=3$ の球から内向き零測地線束を出す。$r=R-\lambda$ で $\lambda=2$ における横断膨張率を求めよ。

:::explanation
面積は $r^2$ に比例し、$\dd r/\dd\lambda=-1$ なので $\theta=-2/r$ である。$\lambda=2$ では $r=1$ で、$\theta=-2$ である。
:::
::::

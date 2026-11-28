---
title: "Einstein 方程式"
description: "物質の保存則と弱重力極限を使って Einstein 方程式を構成し、エネルギー条件と自由度を調べる。"
notation:
  - symbol: "T_{ab}"
    meaning: "エネルギー運動量テンソル"
  - symbol: "\\rho,p"
    meaning: "流体静止系のエネルギー密度と圧力"
  - symbol: "\\Lambda"
    meaning: "宇宙定数"
  - symbol: "\\bar h_{ab}"
    meaning: "跡反転した計量摂動"
---

## この章の目標

- 粒子、完全流体、電磁場で応力と流束を表す。
- 局所的な物理法則を共変化し、その仮定を確認する。
- 曲率と物質を結ぶ方程式の係数を Newton 極限で定める。
- エネルギー条件と計量の二つの伝播自由度を区別する。

## エネルギー運動量テンソル（塵・完全流体・電磁場）

:term[エネルギー運動量テンソル]{en="stress-energy tensor" yomi="えねるぎーうんどうりょうてんそる"} $T^{ab}$ は規格化した観測者 $u$ に対しエネルギー密度 $T_{ab}u^au^b$、空間方向へのエネルギー流束 $-h^a{}_b T^{bc}u_c$ を与える。ただし $h^a{}_b=\delta^a_b+u^au_b$ である。両添字を空間射影したものが応力である。共通の速度 $u$ を持ち衝突を無視できる粒子の集団を塵と呼ぶ。静止系では $T^{00}=\rho$、他成分は零なので $T^{ab}=\rho u^a u^b$ である。

:ref[def-sr-fluid]{full} で導入した完全流体では、静止系の応力が回転で不変なので $p\delta^{ij}$ であり、流束は零である。その定義を局所枠に適用すると

$$
T^{ab}=(\rho+p)u^a u^b+p g^{ab}.
\label{eq-gr-perfect-fluid}
$$

保存則 $\nabla_a T^{ab}=0$ を $u_b$ と $h^c{}_b$ で射影する。$u_b\nabla_a u^b=0$ によって

$$
u(\rho)+(\rho+p)\nabla_a u^a=0,\qquad
(\rho+p)a^c=-h^{cb}\nabla_b p
\label{eq-gr-fluid-conservation}
$$

を得る。塵の密度が非零なら $a=0$ で、$\nabla_a(\rho u^a)=0$ となる。例えば均一な流体要素の体積が二倍になれば塵の $\rho$ は半分となり、$p=\rho/3$ なら $\rho V^{4/3}$ が一定なので $2^{-4/3}$ 倍になる。

電磁場の二形式は局所的に $F=\dd A$ と書く。電磁気の単位はこの章で有理化し $\mu_0=1$ とする。:ref[law-em-covariant]{full} と同じく $\nabla_a F^{ba}=J^b$、力密度を $F^b{}_cJ^c$ と取る。従って $\nabla_aF^{ab}=-J^b$ である。対称テンソルは

$$
T^{ab}_{\rm em}=F^{ac}F^b{}_c-\tfrac14g^{ab}F_{cd}F^{cd}.
$$

$\nabla_{[a}F_{bc]}=0$ に $F^{ac}$ を掛けて $F^{ac}\nabla_a F_{bc}=\tfrac14\nabla_b(F_{cd}F^{cd})$ を得る。よって $\nabla_aT^{ab}_{\rm em}=-F^b{}_cJ^c$ で、荷電体の力密度と和を取ると零である。純電場 $F_{01}=-E$ ではエネルギー密度 $E^2/2$、応力は $\diag(-E^2,E^2,E^2)/2$ になる。電場に平行な引張りと直交方向の圧力を区別できる。

### 運動する観測者が測るエネルギー流

完全流体の $T^{ab}=(\rho+p)u^au^b+pg^{ab}$ を、流体に対して速度 $v$ で動く枠で読む。例えば新しい枠で流体が $x$ 方向へ速度 $v$ を持つとすると $u^{\hat a}=(\gamma,\gamma v,0,0)$ であり、
$$
T^{\hat0\hat0}=\gamma^2(\rho+p)-p,\quad
T^{\hat0\hat1}=\gamma^2v(\rho+p),\quad
T^{\hat1\hat1}=\gamma^2v^2(\rho+p)+p
$$
である。横二方向の成分は $p$ のままである。従って静止系で完全流体の圧力が等方的でも、運動する枠では流束と余分な縦の応力が現れる。跡は $-T^{\hat0\hat0}+T^{\hat1\hat1}+2p=-\rho+3p$ で、$\gamma^2(1-v^2)=1$ により保たれる。

:::example[放射流体の縦の応力を検算する]
$p=\rho/3,v=3/5$ なら $\gamma=5/4$ である。従って $T^{\hat0\hat0}=7\rho/4$、$T^{\hat0\hat1}=5\rho/4$、$T^{\hat1\hat1}=13\rho/12$ である。横二方向を足した跡は $-21\rho/12+13\rho/12+8\rho/12=0$ になる。運動に伴って観測密度が増しても、放射の無跡性は保たれる。
:::

### 観測者ごとの分解と復元

任意の対称二階テンソル $T^{ab}$ に対し、単位時間ベクトル $u^a$ で $\rho=T(u,u)$、$j^a=-h^a{}_c T^{cb}u_b$、$S^{ab}=h^a{}_c h^b{}_d T^{cd}$ を定める。恒等写像を $h^a{}_b-u^au_b$ と分解して二つの添字に作用させると
$$
T^{ab}=\rho u^au^b+u^aj^b+j^au^b+S^{ab}
$$
となる。従って $\rho$ 一成分、$j$ 三成分、$S$ 六成分で十成分を復元できる。等方圧は $p=h_{ab}S^{ab}/3$ で、残り $S^{ab}-ph^{ab}$ は跡が零である。完全流体はある局所枠で $j=0$ かつ応力が等方的な模型である。一般の物質では、エネルギー流が零になる枠があっても応力が等方的になるとは限らない。

:::example[逆向きの二つの光束]
Minkowski 枠で $k_\pm=(1,\pm1,0,0)$ とし、各々のエネルギー密度の係数を同じ $a>0$ とする。$T^{ab}=a k_+^a k_+^b+a k_-^a k_-^b$ を計算すると $\diag(2a,2a,0,0)$ である。静止枠で流束は互いに打ち消すが、$x$ 方向だけに圧力がある。従って $j=0$ だけから完全流体の形は導けない。片方だけにすると $J^a=-T^a{}_bu^b$ はどの未来時間的な $u$ でも光的であり、ある未来単位 $u$ に平行な非零の $J$ を作れない。この光束にはエネルギーの静止枠がない。
:::

### 点粒子のテンソルを分布として求める

端点を考える領域の外に持つ時間的世界線 $z(\tau)$ に対し、共変デルタを $\int \delta_g(x,z)f(x)\sqrt{-g}\,\dd^4x=f(z)$ で定める。点粒子の作用 $S=-m\int\dd\tau$ の逆計量による変分は、世界線を固定して $\delta g_{ab}=-g_{ac}g_{bd}\delta g^{cd}$ を使うと
$$
\delta S=-\frac m2\int u_a u_b\delta g^{ab}(z)\dd\tau
$$
となる。従って $T^{ab}=m\int u^au^b\delta_g(x,z(\tau))\dd\tau$ である。コンパクトな台を持つ試験一形式 $\psi$ を用い、発散を分布の意味で評価すると
$$
\int \psi_b\nabla_a T^{ab}\,\epsilon_g
=-m\int u^a u^b\nabla_a\psi_b\,\dd\tau
=m\int \psi_b a^b\,\dd\tau.
$$
最後は積の法則で $\dv{}{\tau}(u^b\psi_b)$ を取り出し、台の外で端点項を消した。従って、相互作用のない単一の粒子のテンソルが保存することは、その世界線が測地線であることと結び付く。同じ領域に端点があるなら、物質がそこで加わったり失われたりする項が必要である。

### 電磁応力の三次元表示

観測者の正規直交枠で $F_{0i}=-E_i$、$F_{ij}=\varepsilon_{ijk}B_k$ とする。時刻の添字を上げた $F^{0i}=E_i$ と $\sum_j\varepsilon_{ijk}\varepsilon_{j\ell m}$ の縮約を使うと
$$
T^{00}=\frac{E^2+B^2}{2},\qquad
T^{0i}=(\boldsymbol E\times\boldsymbol B)^i,\qquad
T^{ij}=-E_iE_j-B_iB_j+\delta^{ij}\frac{E^2+B^2}{2}
$$
となる。実際、$F_{ab}F^{ab}=2(B^2-E^2)$ であり、空間積の磁場部分は $\delta^{ij}B^2-B_iB_j$ である。$\rho^2-\norm{\boldsymbol E\times\boldsymbol B}^2=(E^2-B^2)^2/4+(\boldsymbol E\cdot\boldsymbol B)^2$ なので、エネルギーの流れはどの観測者に対しても未来因果的または零となる。等号は $E^2=B^2$ かつ $\boldsymbol E\cdot\boldsymbol B=0$ の場合である。

:::example[一方向の電磁波の応力]
$\boldsymbol E=(A,0,0)$、$\boldsymbol B=(0,A,0)$ とすると、エネルギー密度は $A^2$、流束は $(0,0,A^2)$、空間応力は $\diag(0,0,A^2)$ になる。このため $T^{ab}=A^2 k^ak^b$、$k=(1,0,0,1)$ と書ける。一方、$\boldsymbol E=(A,0,0)$、$\boldsymbol B=0$ なら流束は零で、$x$ 方向の応力は $-A^2/2$、残りは $A^2/2$ である。同じ場の振幅でも、流束と応力を記述するには電場と磁場の両方が必要になる。
:::

## 曲がった時空での物理法則と最小結合

:term[最小結合]{en="minimal coupling" yomi="さいしょうけつごう"} は平坦時空の作用にある $\eta$ を $g$ に、直交座標の微分を共変微分に、積分密度を $\sqrt{-\det g}$ に替え、曲率を陽に含む項を加えない処方である。例えばスカラーの作用密度は $-\sqrt{-\det g}(g^{ab}\partial_a\phi\partial_b\phi+m^2\phi^2)/2$ である。$\delta\phi$ の微分を部分積分すると $\frac1{\sqrt{-\det g}}\partial_a(\sqrt{-\det g}g^{ab}\partial_b\phi)-m^2\phi=0$ を得る。$\Gamma^a{}_{ab}=\partial_b\log\sqrt{-\det g}$ は :eqref[eq-gr-christoffel] の縮約と行列式の微分で従うので、この式は $(\nabla_a\nabla^a-m^2)\phi=0$ である。

電磁場の $F_{ab}=2\partial_{[a}A_{b]}=2\nabla_{[a}A_{b]}$ は接続の対称性により同じ形を保つ。一方、平坦時空で偏微分を交換して作った二つの式は、ベクトル場の共変微分に替えると曲率項だけ異なる。このため最小結合は作用の形を指定して行う。共変性自体は、例えば $R\phi^2$ の追加を禁じない。

### 反対称場の発散に必要な接続

反対称テンソル $F^{ab}$ の発散を展開すると
$$
\nabla_aF^{ab}=\partial_aF^{ab}+\Gamma^a{}_{ac}F^{cb}+\Gamma^b{}_{ac}F^{ac}.
$$
最後の項は接続の下添字対称性で消え、中央の項は体積密度の対数微分なので、$\nabla_aF^{ab}=(-g)^{-1/2}\partial_a(\sqrt{-g}F^{ab})$ である。一方、電磁場 $F=\dd A$ では $\dd F=\dd^2 A=0$ であり、偏微分を共変微分に替えても完全反対称化したときに接続が消える。このため接続の成分を全部書き下さずに曲がった背景で Maxwell 方程式を解ける。

:::example[球面を貫く電場の保存]
静的球対称で $g_{tt}g_{rr}=-1$、$r$ は面積半径とする。静止枠の外向き電場を $E(r)$ とすると、$F_{tr}=-E(r)$、$F^{tr}=E(r)$、$\sqrt{-g}=r^2\sin\theta$ なので、真空の $t$ 成分は $(r^2E)'=0$ である。$E=Q/(4\pi r^2)$ と規格化すると、半径 $r$ の球で電束 $4\pi r^2E=Q$ を得る。$\dd F=0$ も $\dd r\wedge\dd t\wedge\dd r=0$ から直ちに成立する。
:::

### 電荷保存と台の条件

曲がった時空でも反対称場の方程式は $\partial_a(\sqrt{-g}F^{ba})=\sqrt{-g}J^b$ と書ける。さらに $\partial_b$ を作用させると、混合偏微分が交換し $F$ の反対称性と打ち消す。従って $\partial_b(\sqrt{-g}J^b)=0$ であり、$\nabla_bJ^b=0$ を得る。保存式に計量の微分が明示的に出ないようにできるのは、$J$ を体積密度と掛けて扱ったからである。

例えば静的な球の中に静止する電荷密度があるとする。正規直交枠の電荷密度を $\rho_{\rm e}$ とすれば $J^t=\rho_{\rm e}/\sqrt{-g_{tt}}$ である。対角な $g_{tt}=-N^2,g_{rr}=A^2$ では外向き電場に対し $F_{tr}=-NAE$、$F^{tr}=E/(NA)$、$\sqrt{-g}=NAr^2\sin\theta$ である。従って
$$
\partial_r(r^2E)=Ar^2\rho_{\rm e}
$$
となる。半径 $r$ の球を貫く電束 $4\pi r^2 E$ の変化は、空間切片の固有体積 $4\pi Ar^2\dd r$ に乗る電荷に等しい。中心まで正則なら $r^2E$ の積分定数は零となる。真空の殻では積分定数が残り、殻の内側にある総電荷で決まる。

### 共変化では区別が必要な微分の順序

ベクトル $V$ について Minkowski 背景では $\partial_a\partial_bV^a=\partial_b\partial_a V^a$ である。Levi-Civita 接続で同じ二通りを書けば
$$
\nabla_a\nabla_bV^a-\nabla_b\nabla_a V^a
=R_{cb}V^c
$$
となる。これは上添字への曲率作用の縮約であり、Ricci が非零の背景では別の方程式を与える。従って平坦時空の成分式を出発点にする場合、平坦時空での恒等変形によって異なる曲率結合が生じないかを確認する必要がある。

:::example[曲率と結合するスカラー]
スカラーの作用に $-\tfrac12\int \xi R\phi^2\sqrt{-g}\,\dd^4x$ を加える。背景計量を固定し、$\phi$ にコンパクトな台をもつ変分を取る。新しい項の変分は $-\xi R\phi\,\delta\phi$ であり、場の式は $(\Box-m^2-\xi R)\phi=0$ になる。定曲率 $K$ の四次元背景では $R=12K$ で、係数は $m^2+12\xi K$ となる。平坦背景の運動方程式は同一でも、この曲率依存性は異なる。これは最小結合を使うかどうかで物理的に選択する項である。
:::

### 共形変換に合わせる結合を決める

四次元で $\tilde g=e^{2f}g$ とし、スカラー場も $\tilde\phi=e^{-f}\phi$ に移す。スカラーの波動演算子を体積密度の式で書けば $\sqrt{-\tilde g}=e^{4f}\sqrt{-g}$、$\tilde g^{ab}=e^{-2f}g^{ab}$ なので
$$
\tilde\Box\psi=e^{-2f}[\Box\psi+2\nabla^a f\,\partial_a\psi]
$$
となる。$\psi=e^{-f}\phi$ を代入すると $\partial\phi$ の交差項が消えて
$$
\tilde\Box\tilde\phi=e^{-3f}[\Box\phi-(\Box f+(\nabla f)^2)\phi]
$$
を得る。一方、Schouten の共形変換を跡に取れば $\tilde R=e^{-2f}[R-6\Box f-6(\nabla f)^2]$ である。従って質量を零にした $(\Box-\xi R)\phi=0$ は、$\xi=1/6$ のときこの変換に共変になる。曲率結合の値は座標共変性だけでは決まらないが、共形変換とともに場を所定の重みで移すという要求を加えると、この値を選べる。

:::example[共形に平坦な背景での平面波]
$\tilde g=a(\eta)^2\eta_{ab}$、$a>0$ とする。質量零で結合が $1/6$ なら $\tilde\phi=a^{-1}\exp[i(\boldsymbol k\cdot\boldsymbol x-|\boldsymbol k|\eta)]$ は方程式を満たす。実際、平坦な $\Box$ を指数関数に作用させると $|\boldsymbol k|^2-\boldsymbol k^2=0$ であり、上の変換でその零が保たれる。一方、最小結合 $\xi=0$ では $a''\ne0$ のとき曲率の項が相殺しない。どの作用から出発したかは、曲がった背景での解に現れる。
:::

## Einstein 方程式の導出

物質側が対称二階テンソルで発散零なので、同じ型の曲率から発散零の組を探す。曲率に線形でその微分を含まず、計量だけで縮約するという仮定では候補が $aR_{ab}+bRg_{ab}+d g_{ab}$ である。四階曲率から二添字を残す縮約は Ricci またはその符号違いで、全縮約には $R$ だけが残るからである。:eqref[eq-gr-einstein-div] により発散は $(a/2+b)\partial_cR$（$c$ は自由添字） である。任意の計量で零にするには $b=-a/2$ とする。$a\ne0$ を仮定して全体を割り、物質との結合定数を Newton 極限で定めると

$$
G_{ab}+\Lambda g_{ab}=8\pi G T_{ab}.
\label{eq-gr-einstein-field}
$$

これを :term[Einstein 方程式]{en="Einstein field equations" yomi="あいんしゅたいんほうていしき"} という。ここで $8\pi G$ は後の Newton 極限で決まる。この議論は採用する候補を絞ったものであり、一般共変性だけから高階微分や追加の場が排除されるわけではない。作用で同じ式が得られることは [作用原理と Hamilton 形式](/physics/general-relativity/action-formulation/) で示す。

四次元で跡を取ると $-R+4\Lambda=8\pi G T$ なので

$$
R_{ab}=8\pi G(T_{ab}-\tfrac12Tg_{ab})+\Lambda g_{ab}.
\label{eq-gr-trace-einstein}
$$

電磁場の跡は $F_{ab}F^{ab}-F_{ab}F^{ab}=0$ である。従って電磁真空では $\Lambda=0$ なら $R=0$ であるが、$R_{ab}$ は電磁場に比例して残る。

### スカラー場を物質源にする

最小結合のスカラーについて $V(\phi)$ をポテンシャルとし、物質作用を $S=-\int\sqrt{-g}[(\nabla\phi)^2/2+V]\,\dd^4x$ とする。$g^{ab}$ を変分する。微分の縮約から $\partial_a\phi\partial_b\phi\,\delta g^{ab}/2$、体積密度から $-\sqrt{-g}g_{ab}\delta g^{ab}/2$ が出る。従って $\delta S=-\tfrac12\int\sqrt{-g}T_{ab}\delta g^{ab}$ で読むテンソルは

$$
T_{ab}=\partial_a\phi\,\partial_b\phi-g_{ab}[(\nabla\phi)^2/2+V].
$$

発散を計算すると第一項の一部 $\partial^a\phi\,\nabla_a\partial_b\phi$ と $(\nabla\phi)^2/2$ の微分が消え、$\nabla^aT_{ab}=(\Box\phi-V'(\phi))\partial_b\phi$ である。場の方程式が成り立つと $T$ の発散も零になる。

:::example[時間にだけ依存するスカラー]
局所枠で $\phi=\phi(t)$ とすると、静止観測者では $\rho=\dot\phi^2/2+V$、$p=\dot\phi^2/2-V$ である。例えば $\dot\phi^2/2=2V$ で $V>0$ なら $p/\rho=1/3$ である。反対に $\dot\phi^2/2\ll V$ なら $p/\rho$ は $-1$ に近く、同じ場でも運動とポテンシャルの比で重力源の性質が変わる。
:::

### 次元と跡反転

$n$ 次元で同じ形の場の式 $R_{ab}-Rg_{ab}/2+\Lambda g_{ab}=\varkappa T_{ab}$ を仮定する。跡を取ると $(1-n/2)R+n\Lambda=\varkappa T$ である。$n>2$ ならこれを解いて
$$
R_{ab}=\varkappa\left(T_{ab}-\frac{T}{n-2}g_{ab}\right)
+\frac{2\Lambda}{n-2}g_{ab}.
$$
$n=4$ で :eqref[eq-gr-trace-einstein] を得る。$n=2$ では跡から $R$ を解けない。実際、曲率の独立成分が一つの二次元では $R_{ab}=Rg_{ab}/2$ であり、Einstein テンソルは恒等的に零となる。従って次元を変えたときの跡操作には例外がある。

:::example[三次元真空と四次元真空の違い]
$n=3,\Lambda=0$ の真空では $R_{ab}=0$ である。三次元では Ricci から Riemann が復元できるので曲率全体が零である。一方、四次元で Ricci が零でも Weyl が残り得る。従って場の式の跡反転が次元ごとに与える情報を、自由な曲率成分が何個あるかと合わせて読む必要がある。
:::

:::example[空間の膨張が要求する密度と圧力]
$a(t)>0$ の平坦切片の計量を考え、宇宙項を零とする。この幾何を生む完全流体の応力を、静止枠の Einstein テンソルから読み取る。 座標の順序を $(t,x,y,z)$ と固定する。成分行列と逆行列は次の組になる。

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

曲率の計算では :eqref[eq-gr-curvature-components] の微分差と積差を分ける。第二列は $\partial_c\Gamma^a{}_{db}$、第三列は $-\partial_d\Gamma^a{}_{cb}$ であり、第四列は接続の積を符号付きで足したものである。三列を足して最後の列を得る。

| 成分 | 第一の微分 | 第二の微分 | 接続の積 | 結果 |
|---|---|---|---|---|
| $R^{t}_{x t x}$ | $a{\left(t \right)} a''{\left(t \right)} + \left(a'{\left(t \right)}\right)^{2}$ | $0$ | $- \left(a'{\left(t \right)}\right)^{2}$ | $a{\left(t \right)} a''{\left(t \right)}$ |
| $R^{t}_{y t y}$ | $a{\left(t \right)} a''{\left(t \right)} + \left(a'{\left(t \right)}\right)^{2}$ | $0$ | $- \left(a'{\left(t \right)}\right)^{2}$ | $a{\left(t \right)} a''{\left(t \right)}$ |
| $R^{t}_{z t z}$ | $a{\left(t \right)} a''{\left(t \right)} + \left(a'{\left(t \right)}\right)^{2}$ | $0$ | $- \left(a'{\left(t \right)}\right)^{2}$ | $a{\left(t \right)} a''{\left(t \right)}$ |
| $R^{x}_{y x y}$ | $0$ | $0$ | $\left(a'{\left(t \right)}\right)^{2}$ | $\left(a'{\left(t \right)}\right)^{2}$ |
| $R^{x}_{z x z}$ | $0$ | $0$ | $\left(a'{\left(t \right)}\right)^{2}$ | $\left(a'{\left(t \right)}\right)^{2}$ |
| $R^{y}_{z y z}$ | $0$ | $0$ | $\left(a'{\left(t \right)}\right)^{2}$ | $\left(a'{\left(t \right)}\right)^{2}$ |

第一添字を計量で下ろすと $R_{abcd}=g_{ak}R^k{}_{bcd}$ となる。対の反対称性と交換対称性を用いて、非零の独立な組を下表に示す。

| 成分 | 結果 |
|---|---|
| $R_{t x t x}$ | $- a{\left(t \right)} a''{\left(t \right)}$ |
| $R_{t y t y}$ | $- a{\left(t \right)} a''{\left(t \right)}$ |
| $R_{t z t z}$ | $- a{\left(t \right)} a''{\left(t \right)}$ |
| $R_{x y x y}$ | $a^{2}{\left(t \right)} \left(a'{\left(t \right)}\right)^{2}$ |
| $R_{x z x z}$ | $a^{2}{\left(t \right)} \left(a'{\left(t \right)}\right)^{2}$ |
| $R_{y z y z}$ | $a^{2}{\left(t \right)} \left(a'{\left(t \right)}\right)^{2}$ |

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

$a(t)=(t/t_0)^{2/3}$ とすると $p=0$ まで簡約される。$a=(t/t_0)^{1/2}$ とすると $p=\rho/3$ となる。この差は $R_{tt}$ に入る $\ddot a$ に関係する。いずれも共動観測者には熱流がなく、密度と圧力はそれぞれ $G_{tt}/(8\pi G)$ と $G_{ii}/(8\pi G a^2)$ である。
:::

### 物質作用の座標不変性と発散

スカラーなどの物質場を $\psi$ とし、作用 $S_m[g,\psi]$ の計量変分を $\delta S_m=-\tfrac12\int T_{ab}\delta g^{ab}\epsilon_g$ と規格化する。台がコンパクトなベクトル $X$ の流れで計量も物質も同時に引き戻す。座標不変な作用では積分の値は変わらない。物質が自分の場の方程式を満たす場合、物質変分の項は零である。計量の逆の変分は $\mathcal L_Xg^{ab}=-2\nabla^{(a}X^{b)}$ だから、
$$
0=\int T_{ab}\nabla^aX^b\,\epsilon_g
=-\int(\nabla^aT_{ab})X^b\,\epsilon_g
$$
である。任意の $X$ を取れるため $\nabla^aT_{ab}=0$ となる。もし別の外部場を固定した作用なら、その外部場を一緒に変換したときの変分が追加され、考えている物質部分の発散が単独で零になるとは限らない。例えば電磁場と荷電体の間では力密度が二つの発散の中で相殺する。

### Einstein テンソルが必要になる弱い場の例

曲率テンソルの候補として $R_{ab}$ だけを取ると、$\nabla^aR_{ab}=\partial_bR/2$ なので $T_{ab}$ の発散が零という条件の下では $R$ が一定でなければならない。跡方程式と合わせると $T$ も一定になる。例えば空間内で密度が変わる静的な低圧物質では $T\simeq-\rho$ は一定ではない。従ってこの候補では一般の源を扱えない。$R_{ab}-Rg_{ab}/2$ の両項は、その発散が計量の恒等式として打ち消し合うように組み合わせている。

:::example[跡を取った真空方程式の意味]
四次元で $\Lambda=0$ の場合、真空方程式は $R_{ab}=0$ と同値である。例えば曲率成分がある点で $R_{0101}=-2b$、$R_{0202}=R_{0303}=b$、$R_{2323}=2b$、$R_{1212}=R_{1313}=-b$ の直交枠を取る。Ricci の時間成分は $-2b+b+b=0$、第一空間成分は $2b-b-b=0$ であり、残りも相殺する。従って $b\ne0$ なら真空の代数条件を満たしても全曲率は消えない。:ref[def-gr-riemann]{full} の意味での曲率と、場の方程式が点ごとに決める Ricci の違いが分かる。
:::

## 宇宙項

定数 $\Lambda$ は長さの逆二乗の単位を持つ。右辺へ移すと $T^\Lambda_{ab}=-\Lambda g_{ab}/(8\pi G)$ となる。完全流体の形と比べると $\rho_\Lambda=\Lambda/(8\pi G)$、$p_\Lambda=-\rho_\Lambda$ で、任意の観測者に対して同一である。例えば $\Lambda=3/L^2$ の真空では $R=12/L^2$ であり、$L$ を倍にすれば曲率は $1/4$ になる。$R_{ab}=\Lambda g_{ab}$ は Weyl 曲率を制限しないため、これだけで定曲率とは結論できない。$C=0$ を加えると前章の分解で $K=\Lambda/3$ が従う。

### 真空エネルギーによる密度と応力

局所 Lorentz 変換に不変な対称二階テンソルを考える。空間回転によって静止枠で $\diag(\rho,p,p,p)$ となる。$x$ 方向の boost を施すと $tx$ 成分は速度に比例して $(\rho+p)$ を含むため、どの速度でも混合成分が零なら $p=-\rho$ である。この形は $T^{\rm vac}_{ab}=-\rho_{\rm vac}g_{ab}$ であり、任意の局所 Lorentz 変換に不変である。さらに発散が零なら $\partial_b\rho_{\rm vac}=0$ である。場の式の右辺から左辺へ移すと観測に現れる宇宙項は $\Lambda+8\pi G\rho_{\rm vac}$ となる。

:::example[観測した密度の組み合わせ]
裸の宇宙項を $\Lambda=-L^{-2}$、真空密度を $\rho_{\rm vac}=1/(4\pi G L^2)$ とする。和は $-L^{-2}+2L^{-2}=L^{-2}$ で、他の物質を含まない解では $R_{ab}=L^{-2}g_{ab}$ となる。同じ有効宇宙項を与える別の分解もあり、古典背景の測定だけでは二つの寄与を区別できない。
:::

:::example[宇宙項と質量を同じ曲率で区別する]
$f=1-2m/r-H^2r^2$ が正である静的パッチを使う。質量の曲率と空間全体に入る宇宙項を同じ式で求める。 座標の順序を $(t,r,\theta,\phi)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}H^{2} r^{2} + \frac{2 m}{r} - 1 & 0 & 0 & 0\\0 & \frac{1}{- H^{2} r^{2} - \frac{2 m}{r} + 1} & 0 & 0\\0 & 0 & r^{2} & 0\\0 & 0 & 0 & r^{2} \sin^{2}{\left(\theta \right)}\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}\frac{r}{H^{2} r^{3} + 2 m - r} & 0 & 0 & 0\\0 & - \frac{H^{2} r^{3} + 2 m - r}{r} & 0 & 0\\0 & 0 & \frac{1}{r^{2}} & 0\\0 & 0 & 0 & \frac{1}{r^{2} \sin^{2}{\left(\theta \right)}}\end{matrix}\right]\\
\det g&=- r^{4} \sin^{2}{\left(\theta \right)}
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{t}_{t r}$ | $\frac{r}{2 \left(H^{2} r^{3} + 2 m - r\right)}\left(\partial_{r}\left(H^{2} r^{2} + \frac{2 m}{r} - 1\right)\right)$ | $\frac{H^{2} r^{3} - m}{r \left(H^{2} r^{3} + 2 m - r\right)}$ |
| $\Gamma^{r}_{t t}$ | $- \frac{H^{2} r^{3} + 2 m - r}{2 r}\left(-\partial_{r}\left(H^{2} r^{2} + \frac{2 m}{r} - 1\right)\right)$ | $\frac{\left(H^{2} r^{3} - m\right) \left(H^{2} r^{3} + 2 m - r\right)}{r^{3}}$ |
| $\Gamma^{r}_{r r}$ | $- \frac{H^{2} r^{3} + 2 m - r}{2 r}\left(\partial_{r}\left(\frac{1}{- H^{2} r^{2} - \frac{2 m}{r} + 1}\right)+\partial_{r}\left(\frac{1}{- H^{2} r^{2} - \frac{2 m}{r} + 1}\right)-\partial_{r}\left(\frac{1}{- H^{2} r^{2} - \frac{2 m}{r} + 1}\right)\right)$ | $- \frac{H^{2} r^{3} - m}{r \left(H^{2} r^{3} + 2 m - r\right)}$ |
| $\Gamma^{r}_{\theta \theta}$ | $- \frac{H^{2} r^{3} + 2 m - r}{2 r}\left(-\partial_{r}\left(r^{2}\right)\right)$ | $H^{2} r^{3} + 2 m - r$ |
| $\Gamma^{r}_{\phi \phi}$ | $- \frac{H^{2} r^{3} + 2 m - r}{2 r}\left(-\partial_{r}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\left(H^{2} r^{3} + 2 m - r\right) \sin^{2}{\left(\theta \right)}$ |
| $\Gamma^{\theta}_{r \theta}$ | $\frac{1}{2 r^{2}}\left(\partial_{r}\left(r^{2}\right)\right)$ | $\frac{1}{r}$ |
| $\Gamma^{\theta}_{\phi \phi}$ | $\frac{1}{2 r^{2}}\left(-\partial_{\theta}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $- \frac{\sin{\left(2 \theta \right)}}{2}$ |
| $\Gamma^{\phi}_{r \phi}$ | $\frac{1}{2 r^{2} \sin^{2}{\left(\theta \right)}}\left(\partial_{r}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{r}$ |
| $\Gamma^{\phi}_{\theta \phi}$ | $\frac{1}{2 r^{2} \sin^{2}{\left(\theta \right)}}\left(\partial_{\theta}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{\tan{\left(\theta \right)}}$ |

曲率の計算では :eqref[eq-gr-curvature-components] の微分差と積差を分ける。第二列は $\partial_c\Gamma^a{}_{db}$、第三列は $-\partial_d\Gamma^a{}_{cb}$ であり、第四列は接続の積を符号付きで足したものである。三列を足して最後の列を得る。

| 成分 | 第一の微分 | 第二の微分 | 接続の積 | 結果 |
|---|---|---|---|---|
| $R^{t}_{r t r}$ | $0$ | $- \frac{3 H^{2} r}{H^{2} r^{3} + 2 m - r} - \frac{\left(- 3 H^{2} r^{2} + 1\right) \left(H^{2} r^{3} - m\right)}{r \left(H^{2} r^{3} + 2 m - r\right)^{2}} + \frac{H^{2} r^{3} - m}{r^{2} \left(H^{2} r^{3} + 2 m - r\right)}$ | $- \frac{\left(H^{2} r^{3} - m\right)^{2}}{r^{2} \left(H^{2} r^{3} + 2 m - r\right)^{2}} - \frac{\left(H^{2} r^{3} - m\right)^{2}}{r^{2} \left(H^{2} r^{3} + 2 m - r\right)^{2}}$ | $- \frac{H^{2} r^{3} + 2 m}{r^{2} \left(H^{2} r^{3} + 2 m - r\right)}$ |
| $R^{t}_{\theta t \theta}$ | $0$ | $0$ | $\frac{H^{2} r^{3} - m}{r}$ | $\frac{H^{2} r^{3} - m}{r}$ |
| $R^{t}_{\phi t \phi}$ | $0$ | $0$ | $\frac{\left(H^{2} r^{3} - m\right) \sin^{2}{\left(\theta \right)}}{r}$ | $\frac{\left(H^{2} r^{3} - m\right) \sin^{2}{\left(\theta \right)}}{r}$ |
| $R^{r}_{\theta r \theta}$ | $3 H^{2} r^{2} - 1$ | $0$ | $- \frac{H^{2} r^{3} - m}{r} + \frac{- H^{2} r^{3} - 2 m + r}{r}$ | $\frac{H^{2} r^{3} - m}{r}$ |
| $R^{r}_{\phi r \phi}$ | $\left(3 H^{2} r^{2} - 1\right) \sin^{2}{\left(\theta \right)}$ | $0$ | $- \frac{\left(H^{2} r^{3} - m\right) \sin^{2}{\left(\theta \right)}}{r} - \frac{\left(H^{2} r^{3} + 2 m - r\right) \sin^{2}{\left(\theta \right)}}{r}$ | $\frac{\left(H^{2} r^{3} - m\right) \sin^{2}{\left(\theta \right)}}{r}$ |
| $R^{\theta}_{\phi \theta \phi}$ | $- \cos{\left(2 \theta \right)}$ | $0$ | $\frac{\sin{\left(2 \theta \right)}}{2 \tan{\left(\theta \right)}} + \frac{\left(H^{2} r^{3} + 2 m - r\right) \sin^{2}{\left(\theta \right)}}{r}$ | $\frac{\left(H^{2} r^{3} + 2 m\right) \sin^{2}{\left(\theta \right)}}{r}$ |

第一添字を計量で下ろすと $R_{abcd}=g_{ak}R^k{}_{bcd}$ となる。対の反対称性と交換対称性を用いて、非零の独立な組を下表に示す。

| 成分 | 結果 |
|---|---|
| $R_{t r t r}$ | $- \frac{H^{2} r^{3} + 2 m}{r^{3}}$ |
| $R_{t \theta t \theta}$ | $\frac{\left(H^{2} r^{3} - m\right) \left(H^{2} r^{3} + 2 m - r\right)}{r^{2}}$ |
| $R_{t \phi t \phi}$ | $\frac{\left(H^{2} r^{3} - m\right) \left(H^{2} r^{3} + 2 m - r\right) \sin^{2}{\left(\theta \right)}}{r^{2}}$ |
| $R_{r \theta r \theta}$ | $- \frac{H^{2} r^{3} - m}{H^{2} r^{3} + 2 m - r}$ |
| $R_{r \phi r \phi}$ | $- \frac{\left(H^{2} r^{3} - m\right) \sin^{2}{\left(\theta \right)}}{H^{2} r^{3} + 2 m - r}$ |
| $R_{\theta \phi \theta \phi}$ | $r \left(H^{2} r^{3} + 2 m\right) \sin^{2}{\left(\theta \right)}$ |

Ricci テンソルは第一・第三添字を縮約して得る。曲率が残る場合も相殺が起こるので、各項を落とさずに加える。

| 成分 | 縮約の項 | 代入 | 結果 |
|---|---|---|---|
| $R_{t t}$ | $R^{r}_{t r t}+R^{\theta}_{t \theta t}+R^{\phi}_{t \phi t}$ | $\frac{\left(H^{2} r^{3} - m\right) \left(H^{2} r^{3} + 2 m - r\right)}{r^{4}} + \frac{\left(H^{2} r^{3} - m\right) \left(H^{2} r^{3} + 2 m - r\right)}{r^{4}} + \frac{\left(H^{2} r^{3} + 2 m\right) \left(H^{2} r^{3} + 2 m - r\right)}{r^{4}}$ | $\frac{3 H^{2} \left(H^{2} r^{3} + 2 m - r\right)}{r}$ |
| $R_{r r}$ | $R^{t}_{r t r}+R^{\theta}_{r \theta r}+R^{\phi}_{r \phi r}$ | $- \frac{H^{2} r^{3} - m}{r^{2} \left(H^{2} r^{3} + 2 m - r\right)} - \frac{H^{2} r^{3} - m}{r^{2} \left(H^{2} r^{3} + 2 m - r\right)} - \frac{H^{2} r^{3} + 2 m}{r^{2} \left(H^{2} r^{3} + 2 m - r\right)}$ | $- \frac{3 H^{2} r}{H^{2} r^{3} + 2 m - r}$ |
| $R_{\theta \theta}$ | $R^{t}_{\theta t \theta}+R^{r}_{\theta r \theta}+R^{\phi}_{\theta \phi \theta}$ | $\frac{H^{2} r^{3} - m}{r} + \frac{H^{2} r^{3} - m}{r} + \frac{H^{2} r^{3} + 2 m}{r}$ | $3 H^{2} r^{2}$ |
| $R_{\phi \phi}$ | $R^{t}_{\phi t \phi}+R^{r}_{\phi r \phi}+R^{\theta}_{\phi \theta \phi}$ | $\frac{\left(H^{2} r^{3} - m\right) \sin^{2}{\left(\theta \right)}}{r} + \frac{\left(H^{2} r^{3} - m\right) \sin^{2}{\left(\theta \right)}}{r} + \frac{\left(H^{2} r^{3} + 2 m\right) \sin^{2}{\left(\theta \right)}}{r}$ | $3 H^{2} r^{2} \sin^{2}{\left(\theta \right)}$ |

逆計量と縮約すると Ricci スカラーとなる。さらに $G_{ab}=R_{ab}-g_{ab}R/2$ を成分ごとに計算する。

$$
\begin{aligned}
R&=\sum_{a,b}g^{ab}R_{ab}\\
&=3 H^{2} + 3 H^{2} + 3 H^{2} + 3 H^{2}\\
&=12 H^{2}
\end{aligned}
$$

非零の Einstein テンソル成分は次のとおりである。

| 成分 | Ricci 項 | 跡の項 | 差 |
|---|---|---|---|
| $G_{t t}$ | $\frac{3 H^{2} \left(H^{2} r^{3} + 2 m - r\right)}{r}$ | $6 H^{2} \left(H^{2} r^{2} + \frac{2 m}{r} - 1\right)$ | $- \frac{3 H^{2} \left(H^{2} r^{3} + 2 m - r\right)}{r}$ |
| $G_{r r}$ | $- \frac{3 H^{2} r}{H^{2} r^{3} + 2 m - r}$ | $\frac{6 H^{2}}{- H^{2} r^{2} - \frac{2 m}{r} + 1}$ | $\frac{3 H^{2} r}{H^{2} r^{3} + 2 m - r}$ |
| $G_{\theta \theta}$ | $3 H^{2} r^{2}$ | $6 H^{2} r^{2}$ | $- 3 H^{2} r^{2}$ |
| $G_{\phi \phi}$ | $3 H^{2} r^{2} \sin^{2}{\left(\theta \right)}$ | $6 H^{2} r^{2} \sin^{2}{\left(\theta \right)}$ | $- 3 H^{2} r^{2} \sin^{2}{\left(\theta \right)}$ |

Ricci は $3H^2g_{ab}$ である。従って真空方程式の宇宙項は $\Lambda=3H^2$ であり、$m$ はこの縮約で相殺する。曲率二乗では $48m^2/r^6+24H^4$ が残る。断面成分に $m$ の項があったので、それが Ricci で消えても全曲率からは消えない。地平線の根は $m$ と $H$ の両方に依存する。
:::

### 正と負の宇宙項を局所で識別する

最大対称な四次元真空では $K=\Lambda/3$ である。時間的測地線に沿う平行枠で初期の相対速度を零にすると、偏差方程式は $J''=KJ$ となる。従って $K>0$ なら $J=L_0\cosh(\sqrt K\,\tau)$、$K<0$ なら $J=L_0\cos(\sqrt{-K}\,\tau)$ である。短い時間ではどちらも $J=L_0(1+\Lambda\tau^2/6)+\order{\tau^4}$ となり、$\Lambda$ の符号に応じて伸びと縮みが入れ替わる。この計算は一定曲率の真空であることを使う。一般の宇宙項付き真空には Weyl の寄与が加わり、三軸の効果は同じとは限らない。

:::example[局所のエネルギー密度と次元]
曲率半径を $L$、$\Lambda=\pm3/L^2$ とする。SI に戻すと宇宙項の有効エネルギー密度は $\rho_{\Lambda,E}=\Lambda c^4/(8\pi G)$ で、その質量密度換算は $\rho_{\Lambda,E}/c^2=\Lambda c^2/(8\pi G)$ である。例えば $L$ を三倍にすると両密度は $1/9$ になり、偏差方程式での時間尺度は三倍になる。SI の固有時を秒で表すなら、$K$ の長さの逆二乗を時間の逆二乗へ換える際に $c^2$ を掛ける。
:::

:::example[同じ代数条件を満たす二つの真空]
$\Lambda>0$ のとき $R_{ab}=\Lambda g_{ab}$ から全ての観測者に対して $R_{ab}u^au^b=-\Lambda$ となる。先の球対称な例で $m=0$ とすれば Weyl が零で、$m\ne0$ とすれば $R_{abcd}R^{abcd}$ の $48m^2/r^6$ が残る。従って Ricci の測定だけでは両者を区別できないが、各方向の潮汐の差には質量の寄与を読むことができる。宇宙項の測定と、自由な曲率成分の測定が独立の情報を与える例である。
:::

## Newton 極限

計量を静的で $g_{00}=-(1+2\Phi)$、$g_{0i}=0$、$g_{ij}=\delta_{ij}+\order{\Phi}$ と仮定する。$\abs{\Phi}\ll1$、$v^2\ll1$、$p\ll\rho$ とし、$\Lambda=0$ とする。測地線の空間式の主要項は $\ddot x^i=-\Gamma^i{}_{00}=-\partial_i\Phi$ である。従って $\Phi$ は Newton ポテンシャルと同定される。曲率の主要項は $R_{00}=\partial_i\Gamma^i{}_{00}=\Delta\Phi$ である。$T_{00}\simeq\rho,T\simeq-\rho$ を :eqref[eq-gr-trace-einstein] に入れると右辺は $4\pi G\rho$ となり :eqref[eq-gr-newton-poisson] と一致する。未定係数を $k$ とした場合は右辺が $k\rho/2$ なので $k=8\pi G$ が決まる。SI で $T$ をエネルギー密度で書くと係数は $8\pi G/c^4$ となる。

圧力も数えると、流体静止枠で $R_{\hat0\hat0}=4\pi G(\rho+3p)-\Lambda$ である。ただし高圧で静的な計量の全成分を $\Phi$ 一つで決めることはできない。$\Lambda>0$ の弱重力真空では $\Delta(-\Lambda r^2/6)=-\Lambda$ であり、$\vb a=\Lambda\vb r/3$ の効果が現れる。

### 圧力が残る弱い静的極限

静的な弱い場で $R_{00}\simeq\Delta\Phi$ とする。物質に静止する時間軸が計量の静的な時間軸に一致し、その正規直交枠で $T_{\hat a\hat b}=\diag(\rho,p_1,p_2,p_3)$ なら $T=-\rho+p_1+p_2+p_3$ である。:eqref[eq-gr-trace-einstein] により
$$
\Delta\Phi\simeq4\pi G(\rho+p_1+p_2+p_3)-\Lambda
$$
となる。これは点での関係であり、圧力が大きい有限容器の遠方重力を求めるときは容器の応力も源に加える。平坦背景の静的な保存源が空間でコンパクトな台を持つなら $\partial_i T^{ij}=0$ であり、$\partial_i(x^kT^{ij})=T^{kj}$ を全空間で積分すると $\int T^{kj}\dd^3x=0$ になる。内部の正の圧力は容器の応力と積分で相殺し得る。相殺を使えるのは系全体の保存が成り立つ場合である。

:::example[宇宙項の静的ポテンシャル]
球対称の真空で $\Delta\Phi=-\Lambda$ とすると $(r^2\Phi')'=-\Lambda r^2$ である。積分して $\Phi=-GM/r-\Lambda r^2/6+C$ となる。粒子の動径加速度は $-GM/r^2+\Lambda r/3$ である。正の $\Lambda$ と $M$ で両項が等しい半径は $(3GM/\Lambda)^{1/3}$ となる。この半径まで弱い場が保たれるかは $GM/r$ と $\abs{\Lambda}r^2$ の双方で確認する。
:::

### 空間計量の摂動を同時に解く

静的に線形化し、空間座標を等方的に選べる場合に $h_{00}=-2\Phi,h_{ij}=-2\Psi\delta_{ij},h_{0i}=0$ とする。Christoffel 記号の一次項は
$$
\Gamma^0{}_{0i}=\partial_i\Phi,\quad
\Gamma^i{}_{00}=\partial_i\Phi,\quad
\Gamma^i{}_{jk}=-\delta^i_k\partial_j\Psi-\delta^i_j\partial_k\Psi+\delta_{jk}\partial^i\Psi.
$$
接続の積を落として $R_{ab}=\partial_c\Gamma^c{}_{ab}-\partial_b\Gamma^c{}_{ac}$ を使うと
$$
\begin{aligned}
R_{00}&=\Delta\Phi,\\
R_{ij}&=\partial_i\partial_j(\Psi-\Phi)+\delta_{ij}\Delta\Psi,\\
R&=4\Delta\Psi-2\Delta\Phi.
\end{aligned}
$$
従って $G_{00}=2\Delta\Psi$、$G_{ij}=(\partial_i\partial_j-\delta_{ij}\Delta)(\Psi-\Phi)$ である。$\Lambda=0$ で応力を先頭で無視するなら $G_{ij}=0$ である。跡から $\Delta(\Psi-\Phi)=0$、その結果 Hessian も零になる。連結領域では差はアフィン関数であり、全空間の孤立系で両方が遠方で零なら差も零になる。この結論には静的な境界条件と空間計量の形を使っている。

:::example[平坦な空間切片だけを固定してはいけない理由]
上の線形化で $\Psi=0$ と仮定すれば、$\Phi$ をどう選んでも $G_{00}=0$ である。従って $T_{00}=\rho\ne0$ の源と両立しない。例えば $\Phi=a(x^2+y^2+z^2)$ で $\Psi=0$ なら $R_{00}=6a$ だが、$G_{00}$ は零で $G_{ij}=4a\delta_{ij}$ である。$R_{00}$ だけを Newton の Poisson 方程式に合わせても、場の方程式の残る成分を満たすとは限らない。$\Psi=\Phi$ を用いれば $G_{ij}=0$ となり、$G_{00}=12a=8\pi G\rho$ を同時に満たせる。
:::

### SI に戻して時刻と密度の尺度をそろえる

座標を $x^0=ct$、ポテンシャルを $\Phi_{\rm SI}$ とすると、$g_{00}=-1-2\Phi_{\rm SI}/c^2$ である。低速の式を秒で書くには $\dd^2x^i/\dd t^2=-c^2\Gamma^i{}_{00}=-\partial_i\Phi_{\rm SI}$ とする。$R_{00}\simeq\Delta\Phi_{\rm SI}/c^2$ に対し、静止質量密度 $\rho_{\rm m}$ のエネルギー密度は $\rho_{\rm m}c^2$ である。従って跡反転した源は $\rho_{\rm m}c^2/2$ が主要項であり、結合定数 $8\pi G/c^4$ を掛けて $4\pi G\rho_{\rm m}/c^2$ になる。両辺の $c^{-2}$ を消すと通常の Poisson 方程式が得られる。圧力を Pa で測ると、その単位は $\mathrm{J\,m^{-3}}$ と同じなので、活性質量密度には $3p/c^2$ が加わる。

## エネルギー条件

:term[エネルギー条件]{en="energy conditions" yomi="えねるぎーじょうけん"} は観測者による測定や測地線の収束に必要な代数的不等式である。零条件は全光的 $k$ について $T(k,k)\ge0$、弱条件は全時間的 $v$ について $T(v,v)\ge0$ である。優勢条件は全未来時間的 $v$ に対し $-T^a{}_b v^b$ が未来因果的または零であること、強条件は全時間的 $v$ に対し $(T_{ab}-Tg_{ab}/2)v^av^b\ge0$ である。

静止枠で異方的応力を対角化でき、$T_{AB}=\diag(\rho,p_1,p_2,p_3)$ とする。$v=v^0(1,\vb w)$、$\abs{\vb w}<1$ なら $T(v,v)/(v^0)^2=\rho+\sum_i p_i(w^i)^2$ である。$\vb w=0$ と各軸に沿う光極限から、弱条件の必要条件は $\rho\ge0$ と全ての $\rho+p_i\ge0$ である。式を $\rho(1-\abs{\vb w}^2)+\sum_i(\rho+p_i)(w^i)^2$ と書けば十分性も分かる。光極限で零条件には後者だけが残る。流束は $v^0(\rho,-p_iw^i)$ なので優勢条件は $\rho\ge\abs{p_i}$ と同値である。跡反転したテンソルの時間成分は $(\rho+\sum_i p_i)/2$、零評価は元と同じなので、強条件は $\rho+\sum_i p_i\ge0$ かつ $\rho+p_i\ge0$ である。正の $\rho_\Lambda$ の流体は弱・零・優勢条件を満たすが $\rho_\Lambda+3p_\Lambda=-2\rho_\Lambda$ なので強条件を破る。

### 主圧力でエネルギー条件を判定する

直交枠で $T_{ab}=\diag(\rho,p_1,p_2,p_3)$ とする。未来時間的な $u=\gamma(1,v^i)$ に対し $T_{ab}u^au^b=\gamma^2(\rho+\sum_i p_iv_i^2)$ である。従って弱エネルギー条件は $\rho\ge0$ かつ $\rho+p_i\ge0$ と同値である。必要性には $v=0$ と軸方向の $\abs v\to1$ を用い、十分性には $\rho(1-\sum v_i^2)+\sum(\rho+p_i)v_i^2$ と書き換える。光的な場合は $\sum v_i^2=1$ であり、零エネルギー条件には $\rho\ge0$ を別に加える必要はない。優勢条件は $-T^a{}_b u^b=\gamma(\rho,-p_iv^i)$ が未来因果的であることから $\rho\ge\abs{p_i}$ と同値である。強条件は $T_{ab}-Tg_{ab}/2$ に弱条件の計算を適用すると $\rho+\sum_i p_i\ge0$ と $\rho+p_i\ge0$ になる。

:::example[一方向の張力がある場合]
$(\rho,p_1,p_2,p_3)=(2,-2,1,1)$ とする。零条件の三つの和は $(0,3,3)$、主圧力の絶対値は $(2,1,1)$、$\rho+\sum_i p_i=2$ である。従って四条件をすべて満たす。$p_1=-3$ に替えると $x^1$ 方向の光線に対し $T_{ab}k^ak^b=-1$ になり、零条件が破れる。圧力を平均するとこの失敗を見落とす。
:::

### 零流体と観測者を極限へ動かす検査

光的なベクトル $k$ と $\mu\ge0$ に対して $T_{ab}=\mu k_a k_b$ を考える。これは静止系で対角化する方法をそのまま適用できない。任意の実ベクトル $v$ に対して $T(v,v)=\mu(k\cdot v)^2\ge0$ であり、跡は零なので零・弱・強条件を満たす。$k$ が未来向きなら未来時間的 $v$ に対し $-T^a{}_bv^b=-\mu(k\cdot v)k^a$ が未来光的または零なので優勢条件も成り立つ。例えば光束と同じ方向に速度を光速へ近づけた観測者の測る密度は小さくなるが、有限の速度では非零の光束を静止させられない。

### エネルギー条件の含意と反例

弱条件から零条件が従うのは、任意の光的ベクトルを同じ時間向きの時間的ベクトルで近似し、二次形式の連続性を使えるからである。優勢条件から弱条件が従うのは、未来向きの二本の因果的なベクトルの内積が非正だからである。強条件から零条件も同様の極限で従い、極限では計量の跡項が消える。一方、強条件と弱条件のどちらか一方から他方が従うとは限らない。例えば $\rho=-1$、$p_i=2$ では $\rho+p_i=1$、$\rho+\sum p_i=5$ で強条件を満たすが、静止観測者の密度は負である。正の宇宙項を有効応力に入れた場合は、弱条件を満たしながら強条件を破る。

:::example[物質の条件と曲率の条件]
光的ベクトルでは :eqref[eq-gr-trace-einstein] により $R_{ab}k^a k^b=8\pi G\,T_{ab}k^a k^b$ となり、宇宙項は寄与しない。時間的単位ベクトル $u$ なら $R_{ab}u^a u^b=8\pi G(T_{ab}-Tg_{ab}/2)u^a u^b-\Lambda$ である。従って物質が強条件を満たしていても、$\Lambda>0$ があると時間的な Ricci の縮約が負になる場合がある。例えば物質を含まない一定曲率真空では $R_{ab}u^au^b=-\Lambda$ を得る。
:::

## Einstein 方程式の自由度とゲージ

計量の十成分は座標の取り方に依存する。線形近似 $g=\eta+h$ で $h\mapsto h-2\partial_{(a}\xi_{b)}$ は微小な座標変換を表す。$\bar h_{ab}=h_{ab}-\eta_{ab}h/2$ とし、$\partial^a\bar h_{ab}=0$ を課すと、真空の方程式は $\Box\bar h_{ab}=0$ となる。以下で線形化を計算する。波の伝播と偏光の物理は [線形化重力と重力波](/physics/general-relativity/gravitational-waves/) で扱う。零ノルムを持つ非零の波数共変ベクトル $k$ のモードでは $k^a\bar h_{ab}=0$ が独立な四条件である。残る座標変換の $\Box\xi=0$ も四つの振幅を持ち、モードへの作用は単射である。実際 $\delta\bar h=0$ なら跡反転して $k_a\xi_b+k_b\xi_a=0$ であり、$k$ と縮約が非零になるベクトルを一つ入れると $\xi=0$ が従う。従って平坦背景上の伝播振幅は $10-4-4=2$ 成分である。非線形での初期データと拘束は [作用原理と Hamilton 形式](/physics/general-relativity/action-formulation/) で扱う。

### 調和座標を選ぶ方程式

各座標をスカラーとして波動演算子に入れると $\Box_g x^\mu=-g^{ab}\Gamma^\mu{}_{ab}$ である。従って $\Box_g x^\mu=0$ を課すと接続の縮約を消せる。Ricci の二階微分の主要項を展開すると
$$
R_{\mu\nu}=-\tfrac12g^{ab}\partial_a\partial_b g_{\mu\nu}
+\tfrac12\bigl(g_{\nu\sigma}\partial_\mu \Gamma^\sigma
+g_{\mu\sigma}\partial_\nu\Gamma^\sigma\bigr)
+\text{一階微分以下の項}
$$
となる。ここで $\Gamma^\sigma=g^{ab}\Gamma^\sigma{}_{ab}$ である。二階微分のうち微分添字が $a,b$ にまとまる項が第一項であり、$\partial_\mu,\partial_\nu$ を含む項を接続の縮約の微分にまとめた。従って調和座標で真空方程式の主要部は十成分それぞれに波動演算子を持つ。座標を指定することと、初期面で拘束を満たすことは別の手続きである。

:::example[同じ Minkowski 時空で調和性を失う座標]
Cartesian 座標で $t,x$ は調和的だが $X=e^x$ と置くと $\Box X=e^x$ である。この座標では空間計量は $\dd X^2/X^2$ になり、$\Gamma^X{}_{XX}=-1/X$、$g^{XX}=X^2$ なので $-g^{XX}\Gamma^X{}_{XX}=X$ を得る。曲率が零であっても、調和座標であるかは座標の選択に依存する。
:::

### 線形化で波動方程式が生じる箇所

背景を Minkowski とし、添字を $\eta$ で動かす。$\Gamma^c{}_{ab}$ の一次式を $\partial_c\Gamma^c{}_{ab}-\partial_b\Gamma^c{}_{ac}$ へ入れると
$$
2R^{(1)}_{ab}
=\partial_a\partial^c h_{bc}+\partial_b\partial^c h_{ac}
-\Box h_{ab}-\partial_a\partial_b h,
\qquad
R^{(1)}=\partial_a\partial_bh^{ab}-\Box h.
$$
$C_b=\partial^a\bar h_{ab}$ と置き $R^{(1)}_{ab}-\eta_{ab}R^{(1)}/2$ を整理すれば
$$
G^{(1)}_{ab}=-\tfrac12\Box\bar h_{ab}
+\partial_{(a} C_{b)}-\tfrac12\eta_{ab}\partial^cC_c .
$$
微小座標変換で $\delta h_{ab}=-2\partial_{(a}\xi_{b)}$ とする規約では、$\delta \bar h_{ab}=-\partial_a\xi_b-\partial_b\xi_a+\eta_{ab}\partial_c\xi^c$ である。発散を取ると $\delta C_b=-\Box\xi_b$ が残る。従って $\Box\xi_b=C_b$ を満たす $\xi$ を解けば調和ゲージにでき、すでに調和ゲージなら残留変換は $\Box\xi_b=0$ を満たす。ゲージの選択に用いる波動方程式と、物理的な曲率の伝播を表す波動方程式は、このように役割が異なる。

:::example[非零モードの二つの振幅を具体的に取り出す]
波を $+z$ に進む $\exp[i\omega(z-t)]$、$\omega\ne0$ とする。$h_{ab}$ の複素振幅に対して残留変換を使い、$h_{00}$ を $\xi_0$、$h_{0i}$ を $\xi_i$ で順に零にできる。なぜなら $k_0=-\omega\ne0$ で、$\delta h_{00}=-2ik_0\xi_0$、$\delta h_{0i}=-i(k_0\xi_i+k_i\xi_0)$ だからである。その後、$k^a\bar h_{a0}=0$ は $\bar h_{00}=h/2=0$ を課し、残る $k^a\bar h_{ai}=0$ は $h_{zi}=0$ を課す。こうして残るのは $h_{xx}=-h_{yy}$ と $h_{xy}=h_{yx}$ である。例えば $h_{xx}=A\cos[\omega(z-t)]$ の場合、線形曲率は $R^{(1)}_{0x0x}=-\partial_t^2h_{xx}/2=\omega^2 h_{xx}/2$ となる。従って残った振幅には測地線偏差で読める曲率が対応する。
:::

### 調和条件を初期値から保つ

調和座標を用いて Ricci の波動方程式を解くときは、座標条件が途中でも成り立つことを確認する必要がある。座標の縮約 $\Gamma^\mu=g^{ab}\Gamma^\mu{}_{ab}$ をこの座標で成分とする場を $C^\mu$ とする。これは指定した座標での場であり、別の座標に移すときに $\Gamma^\mu$ をそのまま変換するという意味ではない。跡反転した方程式の左辺と右辺の差を $Q_{ab}=R_{ab}-\Lambda g_{ab}-8\pi G(T_{ab}-Tg_{ab}/2)$ とする。波動型の縮約方程式を $Q_{ab}-\nabla_{(a}C_{b)}=0$ と書けば、$\partial C$ の項は Ricci の主要部にある座標条件の微分を消す。従って $C=0$ を示せば元の方程式に戻る。

$\nabla^a T_{ab}=0$ が物質の方程式から成り立つと仮定する。Bianchi 恒等式から $\nabla^a(Q_{ab}-Qg_{ab}/2)=0$ なので、縮約方程式を代入して
$$
\begin{aligned}
0&=2\nabla^a\nabla_{(a}C_{b)}-\nabla_b\nabla_a C^a\\
&=\Box_g C_b+R_b{}^aC_a
\end{aligned}
$$
となる。最後の項は共変微分の縮約された交換子である。この連立波動方程式の初期値問題が一意な領域では、空間的初期面で $C=0$ と法線微分が零なら、その領域で $C=0$ が保たれる。初期面で $C=0$ だけを課しても、法線微分は一般には零とは限らない。初期面における場の方程式の拘束とゲージ条件が、残る初期値を定める。この構成を非線形な初期値の幾何として扱う方法は [作用原理と Hamilton 形式](/physics/general-relativity/action-formulation/) に進むと得られる。

:::example[平坦背景の線形拘束を切片で確認する]
線形化した真空の縮約方程式を $\Box\bar h_{ab}=0$ とする。従って $G^{(1)}_{ab}=\partial_{(a}C_{b)}-\eta_{ab}\partial^cC_c/2$ である。$t=0$ の切片全体で $C_a=0$ を課すと接線微分 $\partial_i C_a$ はすべて零なので、
$$
G^{(1)}_{00}\big|_{t=0}=\tfrac12\partial_t C_0\big|_{t=0},
\qquad
G^{(1)}_{0i}\big|_{t=0}=\tfrac12\partial_t C_i\big|_{t=0}.
$$
従って $G^{(1)}_{0a}=0$ を同じ切片で要求すれば、$C$ の法線微分も消える。実際、$\Box C_a=\partial^b\Box\bar h_{ab}=0$ であるから、零の初期値を持つ解は $C_a=0$ である。逆に $\bar h_{00}=t^2+x^2$、他成分を零に取ると $\Box\bar h_{00}=-2+2=0$ で、$C_0=-2t$ である。$t=0$ で $C=0$ でも $\partial_t C_0=-2$ なので、元の方程式の $00$ 成分は満たさない。ゲージ条件の初期値と場の拘束を別々に使う理由が分かる。
:::

### 局所保存と座標エネルギーの差

対称な応力テンソルの発散を $\nabla_aT^a{}_b=0$ と書く。密度を用いると
$$
\partial_a(\sqrt{-g}\,T^a{}_b)
=\sqrt{-g}\,\Gamma^c{}_{ab}T^a{}_c
=\tfrac12\sqrt{-g}\,T^{ac}\partial_b g_{ac}
$$
である。最後の等号では $a,c$ の交換で二項が消える。$b=t$ で右辺が消えるのは $(\partial_t g)_{ac}=0$ などの条件がある場合である。一般の時空では局所共変的な保存が成り立っても、座標時間に関する物質エネルギーを空間の座標範囲を固定したまま切片ごとに積分すると、値は一般には時間に依存する。例えば計量 $-\dd t^2+a(t)^2\delta_{ij}\dd x^i\dd x^j$ の共動流体では、座標体積 $V_0$ の領域について $\partial_t(\rho a^3 V_0)=-p\,\partial_t(a^3 V_0)$ となる。右辺は圧力と体積の変化の積であり、流体の保存式を積分しても同じ結果を得る。

## 演習問題

::::exercise{level=basic}
完全流体の跡を求めよ。

:::solution
$u^au_a=-1$、$g^{ab}g_{ab}=4$ だから $T=-(\rho+p)+4p=-\rho+3p$ である。
:::
::::

::::exercise{level=basic}
$\Lambda=2/L^2$ の真空の $R$ を求めよ。

:::solution
$T=0$ の跡方程式から $R=4\Lambda=8/L^2$ である。
:::
::::

::::exercise{level=basic}
$\rho=4,p=2$ の完全流体が満たすエネルギー条件を判定せよ。

:::solution
$\rho>0$、$\rho+p=6$、$\rho-\abs p=2$、$\rho+3p=10$ が非負であり、四条件すべてを満たす。
:::
::::

::::exercise{level=basic}
初期体積 $V_0$ の塵の体積が $3V_0$ になったときの密度を求めよ。

:::solution
:eqref[eq-gr-fluid-conservation] で $p=0$、$\nabla_a u^a=\dot V/V$ なので $(\rho V)^{\cdot}=0$ である。従って $\rho=\rho_0/3$ となる。
:::
::::

::::exercise{level=advanced}
$\Lambda$ を時空関数としたとき、同じ式で $\nabla_a T^{ab}=0$ を仮定して何が従うか示せ。

:::solution
Bianchi により $\nabla_a G^{ab}=0$、計量整合性より宇宙項の発散は $\partial^b\Lambda$ である。したがって $\partial^b\Lambda=0$ となり、$\Lambda$ は各連結成分で一定である。
:::
::::

::::exercise{level=advanced}
静的なポテンシャル $\Phi=A(x^2+y^2+z^2)$ を作る Newton 質量密度を求めよ。

:::solution
$\Delta\Phi=6A$ であり、Poisson 方程式から $\rho=3A/(2\pi G)$ である。$A>0$ なら力 $-2A\vb x$ は内向きである。
:::
::::

::::exercise{level=advanced}
$\rho=1$、$(p_1,p_2,p_3)=(-1,2,0)$ の対角応力について四条件を判定せよ。

:::solution
$\rho+p_i$ は $0,3,1$、$\rho+\sum p_i=2$ なので零・弱・強条件を満たす。だが $p_2=2>\rho$ なので優勢条件を破る。
:::
::::

::::exercise{level=advanced}
電磁場で光的 $k$ に対し零条件を示せ。

:::solution
計量項が消え $T_{ab}k^a k^b=w_cw^c$、$w_c=k^aF_{ac}$ である。反対称性から $w\cdot k=0$ であり、正規枠で $k$ を $(1,1,0,0)$ の倍とすれば $w^0=w^1$、従って $w_cw^c=(w^2)^2+(w^3)^2\ge0$ である。
:::
::::

::::exercise{level=advanced}
四次元の対角応力を $(\rho,p_1,p_2,p_3)=(\rho,-\rho,-\rho,\rho)$、$\rho>0$ とする。四つのエネルギー条件を判定し、静的弱場の活性重力密度を求めよ。さらに符号規約を $G_{ab}=8\pi G T_{ab}$ として $R_{\hat0\hat0}$ を求めよ。

:::solution
零条件の三つの和は $(0,0,2\rho)$ で非負であり、$\rho>0$ なので弱条件も成り立つ。$\abs{p_i}=\rho$ なので優勢条件を満たす。$\rho+p_1+p_2+p_3=0$ であるため強条件も境界で成り立つ。活性密度の点での組はこの最後の和で零である。跡は $-\rho+p_1+p_2+p_3=-2\rho$ なので
$$
R_{\hat0\hat0}=8\pi G\left(\rho-\tfrac12(-1)(-2\rho)\right)=0.
$$
この一成分が零でも $R_{\hat3\hat3}=16\pi G\rho$ は非零であるため、曲率全体が零とはいえない。
:::
::::

## 確認テスト

::::quiz
$\nabla_aG^{ab}=0$ の根拠を選べ。

- [x] 縮約した Bianchi 恒等式
- [ ] 物質が塵であること
- [ ] 時空が真空であること

:::explanation
第二 Bianchi 恒等式を二度縮約すると $\nabla_aR^{ab}=\partial^bR/2$ である。従って計量整合性から $G$ の発散は零となり、この計算は物質の内容によらない。
:::
::::

::::quiz{answer="2" tol="1e-9"}
$\rho=1,p=1$ の完全流体の跡を求めよ。

:::explanation
静止枠で時間成分は $-1$、空間成分の和は $3$ を与えるので $T=-1+3=2$ である。
:::
::::

::::quiz
正の宇宙定数を完全流体として見たときに満たす条件をすべて選べ。

- [x] 零条件
- [x] 優勢条件
- [ ] 強条件

:::explanation
$p_\Lambda=-\rho_\Lambda$ なので $\rho+p=0$ で零条件、$\rho=\abs p$ で優勢条件が成立する。一方 $\rho+3p=-2\rho_\Lambda<0$ で強条件は破れる。
:::
::::

::::quiz{answer="2" tol="1e-9"}
平坦背景の非零 Fourier モードで、真空の計量の伝播自由度を答えよ。

:::explanation
対称計量の十成分に対して調和ゲージが四条件を課し、残留ゲージでさらに四成分を除く。$10-4-4=2$ である。
:::
::::

::::quiz{answer="6" tol="1e-7"}
主圧力 $(p_1,p_2,p_3)=(-2,1,3)$ の点で $\rho=4$ とする。弱い静的な極限の $\rho+\sum p_i$ を求めよ。

:::explanation
圧力の跡は $-2+1+3=2$ である。従って $00$ 成分の跡反転に現れる組は密度を加えた $6$ である。
:::
::::

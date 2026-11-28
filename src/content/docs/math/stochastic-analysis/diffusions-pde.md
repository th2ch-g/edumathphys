---
title: 拡散過程と偏微分方程式
description: 拡散過程から半群と生成作用素を作り、期待値、密度、退出分布が満たす方程式を導く。
notation:
  - symbol: 'P_t'
    meaning: 拡散の推移半群
  - symbol: '\mathcal L^*'
    meaning: 生成作用素の形式随伴
  - symbol: 'I_{0T}'
    meaning: 小雑音極限の作用汎関数
---

## この章の目標

- 強解から Markov 半群を構成する。
- 後ろ向き方程式と前向き方程式を、それぞれの未知量から導く。
- Feynman–Kac 公式を終端値問題と境界値問題に用いる。
- 不変分布への収束と小雑音極限を計算する。

## Markov 過程と推移半群

この章ではまず時間に依存しない大域 Lipschitz 係数の SDE を用いる。$X^x$ は初期値 $x$ から始めた解とし、$P_tf(x)=\mathbb E f(X_t^x)$ とする。現在の状態が分かれば、その後の解は独立な将来増分を使って作り直せる。

:::proposition[Markov 性と半群]{#prop-sa-diffusion-semigroup}
有界 Borel $f$ について $\mathbb E[f(X_{s+t})\mid\mathcal F_s]=P_tf(X_s)$ が成り立つ。従って $P_{s+t}=P_sP_t$、$P_0=I$ である。同じ条件付き式は有限な停止時刻 $s$ にも成り立つ。
:::

:::proof
時刻 $s$ に解を再出発し、未来増分で逐次近似を作る。初期値が有限個の値なら、各初期値での解と未来増分が過去と独立なことから式が成り立つ。一般の初期値に対しては初期値の連続性と差の積率評価で近似し、まず有界連続 $f$ について結論を得る。$P_tf$ も初期値の連続性で連続である。生成する有界連続関数から単調類で有界 Borel へ延ばす。停止時刻 $\tau$ については上から可算格子で近似し、格子での式に任意の $A\in\mathcal F_\tau$ の指示関数を掛ける。連続な $f$ と $P_tf$、路の連続性で極限を取れる。最後に再び単調類を使う。
:::

$P_t$ は正を保ち $P_t1=1$、最大ノルムを増やさない。例えば Brown 運動なら $P_t$ は Gauss 密度との畳み込みであり、半群の積は独立正規増分の分散の加法性に対応する。加法的 OU では :eqref[eq-sa-ou] の平均と分散を持つ正規核になる。

### 推移確率を核として扱う

集合 $A$ に対して $p_t(x,A)=P_x(X_t\in A)$ と置く。固定した $x$ では確率測度であり、固定した Borel 集合 $A$ では $x$ の Borel 関数である。後者は、連続初期値への依存から有界連続関数の期待値が連続となることと、単調類によって従う。この二つの性質を持つ関数を :term[推移確率核]{en="transition probability kernel" yomi="すいいかくりつかく"} という。半群の等式を指示関数に適用すると Chapman–Kolmogorov の関係
$$
p_{s+t}(x,A)=\int p_s(x,\dd y)p_t(y,A)
$$
を得る。連続時間の有限次元分布も、初期測度にこの核を時間の順で掛ければ求まる。

:::example[三時刻の条件付き平均]
$\dd X=2\,\dd t+\dd W$、$X_0=0$ とする。$X_1=y$ を条件に $X_3$ は $N(y+4,2)$ である。$X_1$ 自体は $N(2,1)$ なので、$X_3$ は二つの独立な正規変数を足して $N(6,3)$ となる。$X_2$ を間に挟むと、条件付き分布の分散は区間ごとに $1+1$、平均の増分は $2+2$ になる。同じ核の畳み合わせをどこで区切っても $p_3$ が変わらないことを数値で確認できる。
:::

### 消える関数上の Feller 半群

最大ノルムだけでは、有界連続関数で各点の連続性と時間のノルム連続性を同時に得られない。空間の無限遠で零に近づく連続関数の空間 $C_0(\R^d)$ は一様連続関数からなり、そこで自然な定義域を選べる。

:::proposition[有界係数の Feller 性]{#prop-sa-pde-feller}
$b,\sigma$ が有界大域 Lipschitz なら、$P_t$ は $C_0(\R^d)$ を保ち、最大ノルムで強連続な正の縮小半群である。
:::

:::proof
$|b|\le K$、$\norm{\sigma}\le K$ とする。等長性で
$$
\mathbb E|X_t^x-x|^2\le 2K^2(t^2+t)
$$
が $x$ に一様に成り立つ。$f\in C_0$ は有界一様連続であり、$\omega_f(\eta)=\sup_{|x-y|\le\eta}|f(x)-f(y)|$ とすると
$$
\sup_x|P_tf(x)-f(x)|
\le\omega_f(\eta)+4\norm f_\infty K^2(t^2+t)/\eta^2.
$$
先に $t\downarrow0$、次に $\eta\downarrow0$ で強連続性が従う。$P_tf$ の空間連続性は :ref[prop-sa-diffusion-semigroup] の証明と同じ初期値の評価を使う。また $|y|\ge R$ で $|f(y)|\le\delta$ となる $R$ を取れば、$|x|>R$ について
$$
|P_tf(x)|\le\delta+\norm f_\infty P_x(|X_t^x-x|\ge|x|-R)
\le\delta+\frac{2\norm f_\infty K^2(t^2+t)}{(|x|-R)^2}.
$$
$|x|\to\infty$ と $\delta\downarrow0$ で $P_tf\in C_0$ が従う。正性と縮小性は期待値の性質であり、強連続性は半群の性質で任意時刻へ移せる。
:::

:::example[OU の場合の無限遠]
一次元の $X_t^x=e^{-t}x+Z_t$、$Z_t\sim N(0,(1-e^{-2t})/2)$ も Feller 半群を与える。$t$ を固定すると $e^{-t}|x|\to\infty$ なので、$f\in C_0$ に対して $P_tf(x)\to0$ は Gauss 裾で従う。$0\le t\le1$ では $e^{-t}\ge e^{-1}$ なのでこの尾部の評価は $t$ に一様である。コンパクトな $x$ では $|(e^{-t}-1)x|$ と $Z_t$ の分散が一様に零へ行く。初期値を無限遠まで動かすとドリフトが有界でないが、領域を分ければ強連続性が得られる。
:::

### 時刻依存係数と二パラメータの族

$b(t,x),\sigma(t,x)$ が時刻に依存すると、$P_{s,t}f(x)=\mathbb E[f(X^{s,x}_t)]$ を用いる。$P_{s,t}=P_{s,u}P_{u,t}$ は中間時刻での条件付き期待値から従い、一般には $t-s$ だけで決まらない。例えば $\dd X=t\,\dd t+\dd W$ の $(s,t)$ 区間の平均増分は $(t^2-s^2)/2$ である。長さ $1$ の区間でも、$(0,1)$ で $1/2$、$(1,2)$ で $3/2$ になる。時刻も状態に加えた $(s+t,X_{s+t})$ は時間同次の過程である。将来の進み方が現在の時刻と空間位置の組から決まるためである。

:::example[境界で殺したときの推移半群の質量]
$X_t=1+B_t$ を零に達した時点で殺し、$P_t^D f(x)=\mathbb E_x[f(X_t)1_{\{\tau_0>t\}}]$ とする。定数 $1$ を入れると生存確率 $P_t^D1(1)=2\Phi(1/\sqrt t)-1$ である。これは $B$ の下向き最大値に反射原理を使って求められ、例えば $t=4$ では $2\Phi(1/2)-1$ になる。半群性は、中間時刻に生存している路だけに残りの生存条件を付けると Markov 性から従う。殺された状態を墓場点として加えれば、その点の質量を含めた全体の推移確率が $1$ となる。
:::

## 生成作用素と Dynkin の公式

:::proposition[Dynkin の公式]{#prop-sa-dynkin}
$f\in C_c^2(\R^d)$ と停止時刻 $\tau\le T$ に対し
$$
\mathbb E_x f(X_\tau)-f(x)=\mathbb E_x\int_0^\tau\mathcal L f(X_s)\,\dd s,\qquad
\mathcal L f=b\cdot\grad f+\tfrac12\tr(a\nabla^2 f).
\label{eq-sa-dynkin}
$$
特に $\lim_{t\downarrow0}(P_tf(x)-f(x))/t=\mathcal L f(x)$ である。
:::

:::proof
伊藤の公式の確率積分は被積分関数が有界なので有限区間で二乗可積分である。有界停止でその平均を零にすると公式が得られる。$\tau=t$ とし、$\mathcal Lf$ の有界連続性と $X_s\to x$ を使って時間平均の極限を取る。
:::

これは :term[生成作用素]{en="infinitesimal generator" yomi="せいせいさようそ"} の点ごとの表示である。Banach 空間上の強生成作用素を定めるにはノルムと定義域を別に選ぶ。例えば Brown 半群は $C_0$ 上で強連続である。実際、$C_0$ 関数は一様連続で、Gauss 増分が小さい領域では連続度で、大きい領域では裾確率で一様に誤差を抑えられる。一般の強連続半群と生成作用素の性質は :ref[def-fa-semigroup] と :ref[lem-fa-generator] を参照する。

### 強生成作用素と積分の微分

$C_0$ 上で $P_t$ が強連続な場合、$t^{-1}(P_tf-f)$ が最大ノルムで収束する $f$ を強生成作用素 $A$ の定義域 $D(A)$ に入れる。この極限は点ごとの形式微分より強い条件を課す。:ref[prop-sa-pde-feller] の仮定の下で $f\in C_c^2$ なら、$\mathcal Lf$ は $C_0$ に属するので
$$
\frac{P_tf-f}{t}-\mathcal Lf
=\frac1t\int_0^t(P_s\mathcal Lf-\mathcal Lf)\,\dd s
$$
の最大ノルムは零へ収束する。従って $f\in D(A)$ かつ $Af=\mathcal Lf$ である。停止で得た Dynkin 公式を半群のノルムに移した形である。

:::example[点ごとの生成作用素と強い定義域の違い]
一次元 Brown 運動で $f(x)=|x|^\alpha\chi(x)$、$1<\alpha<2$ とし、$\chi=1$ が原点の近くで成り立つ非負の滑らかな切断を使う。$P_tf(0)/t=t^{\alpha/2-1}\mathbb E[|Z|^\alpha\chi(\sqrt{t}Z)]$、$Z\sim N(0,1)$ である。期待値は正の $\mathbb E|Z|^\alpha$ へ収束するので、原点での差商は無限大へ発散する。$f\in C_0$ で一階微分が連続でも、生成作用素の定義域に入るとは限らない。
:::

### 割引した時間積分とレゾルベント

:::proposition[レゾルベントの確率表現]{#prop-sa-pde-resolvent}
強連続な $C_0$ 縮小半群 $P_t$ に対し $\lambda>0$、$f\in C_0$ とする。
$$
R_\lambda f(x)=\int_0^\infty e^{-\lambda t}P_tf(x)\,\dd t
             =\mathbb E_x\int_0^\infty e^{-\lambda t}f(X_t)\,\dd t
$$
は $D(A)$ に属し、$(\lambda-A)R_\lambda f=f$、$\norm{R_\lambda}\le1/\lambda$ である。また $u\in D(A)$ について $R_\lambda(\lambda-A)u=u$ である。
:::

:::proof
強連続性と縮小性から最初の積分は $C_0$ に収束し、ノルムは $\norm f/\lambda$ 以下である。Fubini で期待値に移せる。$h>0$ について半群の等式と積分の置換を使うと
$$
P_hR_\lambda f=e^{\lambda h}\left(R_\lambda f-\int_0^h e^{-\lambda r}P_rf\,\dd r\right).
$$
差商を最大ノルムで取ると $AR_\lambda f=\lambda R_\lambda f-f$ を得る。逆側は $u\in D(A)$ について $\dd P_tu/\dd t=P_tAu$ の半群微分を使う。これは $\frac{P_{t+h}u-P_tu}h=P_t\frac{P_hu-u}{h}$ の収束と $P_tAu$ の強連続性から従う。$\frac{\dd}{\dd t}(e^{-\lambda t}P_tu)=-e^{-\lambda t}P_t(\lambda-A)u$ を積分し、無限遠でのノルムが零になることを用いる。
:::

:::example[Brown 運動のレゾルベント核]
一次元 Brown 運動に対し $k_\lambda(z)=e^{-\sqrt{2\lambda}|z|}/\sqrt{2\lambda}$ と置く。$z\ne0$ で $\lambda k_\lambda=k_\lambda''/2$ であり、原点で導関数の右極限と左極限の差は $-2$ である。従って分布微分では $(\lambda-\partial_{zz}/2)k_\lambda=\delta_0$ となる。$f$ を滑らかなコンパクト台とし、$u=k_\lambda*f$ と置けば $(\lambda-\partial_{xx}/2)u=f$、$u,u',u''\in C_0$ である。時間積分した Dynkin 公式はこの $u$ でも、$u'$ の有界性と等長性によって使える。$\lambda u-u''/2=f$ を代入して指数割引した伊藤公式の時間上限を無限にすると、$u=R_\lambda f$ と分かる。$\lambda=2$ では核は $e^{-2|z|}/2$、全質量は $1/2$ になる。
:::

### 局所の計算が停止の期待値になるとき

停止時刻の有界性を外す際は、確率積分の可積分性か、停止点の関数とドリフト積分の優収束を確認する。例えば $X_t=2W_t$ を $|X|=3$ まで止めると、$f(x)=x^2$ に対して $\mathcal L f=4$。有界時刻を併用して $4\mathbb E(\tau\wedge t)=\mathbb E X_{\tau\wedge t}^2\le9$ となるため、まず $\mathbb E\tau\le9/4$ を得る。有限性を確認した後は終端値が $9$ に収束するので $\mathbb E\tau=9/4$ である。$f$ がコンパクト台でなくても停止領域で $f$ とその導関数が有界なら、停止した伊藤公式の正当化は同じである。

:::example[Dynkin の公式で停止前の二乗位置を積分する]
$B_0=0$ とし、$\tau$ は $(-2,2)$ の退出時刻とする。$f(x)=x^4$ では $\mathcal L f=6x^2$ である。$f$ を $[-2,2]$ で同じ値を持つコンパクト台の滑らかな関数に延長し、Dynkin の公式を $t\wedge\tau$ に適用すると $\mathbb E B_{t\wedge\tau}^4=6\mathbb E\int_0^{t\wedge\tau} B_s^2\dd s$ となる。左辺は $16$ 以下で、$\tau<\infty$ a.s. により $16$ に収束する。右辺は単調収束するため、退出までの位置の二乗積分の平均は $16/6=8/3$ となる。位置の終端モーメントから、途中の状態に重みを置く占有量が得られた。
:::

## Kolmogorov の後ろ向き方程式

終端関数 $g$ を固定して初期点 $x$ を変えると、$P_tg(x)$ が後ろ向き方程式の未知関数になる。一方、時刻を終端までの残り時間で表すか現在時刻で表すかによって符号が変わる。

:::proposition[後ろ向き方程式の検証]{#prop-sa-backward}
$u$ が $[0,T]\times\R^d$ 上有界かつ $C^{1,2}$ で $u_t+\mathcal L u=0$、$u(T,x)=g(x)$ を満たすなら
$$
u(t,x)=P_{T-t}g(x).
$$
:::

:::proof
時刻 $t$ から開始した解に伊藤の公式を使うと、$u(s,X_s)$ の時間項は零である。有界な局所マルチンゲールなので終端の期待値が初期値に等しい。従って表現も一意性も得られる。
:::

存在を半群から示す場合には正則性の仮定を要する。例えば Brown 運動では $v(t,x)=P_tg(x)$ は $t>0$ で Gauss 核の微分を積分に移せる。核 $p_t(z)$ は $\partial_t p_t=\Delta p_t/2$ を満たすため $v_t=\Delta v/2$ である。有界連続な $g$ に対し $v(t,x)\to g(x)$ は近似恒等作用素の性質で得られる。:ref[thm-pde-heat-cauchy] に対応する。例として $g(x)=\cos(2x)$ なら $v(t,x)=e^{-2t}\cos(2x)$ で、時刻 $1$ における $x=0$ の値は $e^{-2}$ である。

### 半群の微分順序

半群の定義域の議論では、$g\in D(A)$ なら $P_tg\in D(A)$ かつ
$$
\partial_t P_tg=A P_tg=P_tAg
$$
である。実際、$P_hP_t=P_tP_h$ から $\frac1h(P_h-I)P_tg=P_t\frac1h(P_h-I)g$ をノルムで極限すればよい。この意味で後ろ向き方程式はまず抽象的な半群の微分方程式として解ける。$A$ が $P_tg$ に対して $\mathcal L$ と一致することを別に確かめれば、偏微分の形が得られる。試験関数が十分滑らかである場合でも、係数や境界が特異なら自動的に空間二階微分が存在するわけではない。

:::proposition[滑らかな係数での古典的後ろ向き解]{#prop-sa-pde-backward-smooth}
$b,\sigma$ は時間に依存せず、値と三階までの導関数が有界連続な $C^3$ 関数とする。$g$ は値と三階までの導関数が有界連続な $C^3$ 関数とする。すると $u(t,x)=P_{T-t}g(x)$ は $C^{1,2}$ で $u_t+\mathcal L u=0$、$u(T,x)=g(x)$ を満たす。ここでは $\sigma\sigma^\top$ の正定値性を仮定しない。
:::

:::proof
初期値方向 $v$ の一次変分 $J^v$ は :ref[prop-sa-sde-initial-derivative] の線形方程式を満たす。その証明の差商評価を $L^p$ 最大ノルムに移すこともできる。実際、係数差は一様有界で確率収束し、支配関数には一次変分の $p$ 次最大値を取る。一次変分の全有限積率は、線形方程式と :ref[prop-sa-sde-moment] の証明で有限である。

二つの初期値方向 $v,w$ に対し、二次変分 $K^{v,w}$ を初期値零の線形方程式で定める。ドリフトは $Db(X)K^{v,w}+D^2b(X)[J^v,J^w]$、$k$ 番目のノイズ係数は $D\sigma_k(X)K^{v,w}+D^2\sigma_k(X)[J^v,J^w]$ である。右辺の非斉次項は一次変分の四次積率から二乗可積分である。一次変分方程式を二つの初期点で引いて $h$ で割ると、$Db$ の差商は
$$
\int_0^1 D^2b(X^x+\theta(X^{x+hw}-X^x))
      \left[\frac{X^{x+hw}-X^x}{h},\,\cdot\right]\,\dd\theta .
$$
これは $L^4$ で $D^2b(X^x)[J^w,\cdot]$ へ収束する。差商の $L^4$ 収束と $D^2b$ の Lipschitz 性を用い、残る積は Hölder と高次積率で抑えればよい。拡散項も同じであり、$L^2$ 最大ノルムの安定性評価で一次変分の差商が $K^{v,w}$ へ収束する。

これにより期待値を二回微分できる。具体的に $D_vP_tg=\mathbb E[Dg(X_t)J_t^v]$、$D_{v,w}^2P_tg=\mathbb E[D^2g(X_t)[J_t^v,J_t^w]+Dg(X_t)K_t^{v,w}]$ である。空間点を変えると $X,J,K$ は最大ノルムの各有限 $L^p$ で連続である。$K$ に対しても、方程式の係数および一次変分の連続性と、Cauchy–Schwarz で非斉次項の連続性を得てから同じ安定性評価を使う。時刻の連続性には各路の連続性と最大積率を使う。従って導関数は時空で連続有界である。

固定した $t$ で $P_tg$ に Dynkin 公式を適用すると $\lim_{h\downarrow0}(P_hP_tg-P_tg)/h=\mathcal LP_tg$ が各点で得られる。一方 $g$ に対する有界係数の Dynkin 公式を半群に入れれば、その差商の極限は $P_t\mathcal Lg$ でもあり、これは時間に連続である。積分の形 $P_tg=g+\int_0^tP_s\mathcal Lg\,\dd s$ に戻すと二側の時間微分が得られる。残り時間 $T-t$ を代入して符号を変えれば主張が従う。
:::

:::example[多項式の閉じた後ろ向き方程式]
一次元 OU $\dd X=-2X\,\dd t+2\,\dd W$ で $g(x)=x^2$ とする。$\mathcal Lg=4-4g$ なので、$P_tg$ も $1$ と $x^2$ の線形結合で求められる。実際、$v(t,x)=1+(x^2-1)e^{-4t}$ は $v_t=2v_{xx}-2xv_x$、$v(0,x)=x^2$ を満たす。:eqref[eq-sa-ou] の平均二乗と条件付き分散を足しても同じ式になる。検証で有界解を仮定した命題を直接使うことはできないが、OU の四次最大値が有限なので、$v(T-s,X_s)$ の確率項を二乗可積分にして検証できる。
:::

:::example[実部を取って三角関数の期待値を解く]
定数ドリフト $b=3$、拡散係数 $\sigma=2$ で $\mathcal L=3\partial_x+2\partial_{xx}$ とする。$g(x)=\cos x$ に対して複素指数 $e^{ix}$ は固有値 $3i-2$ を持つ。従って $P_tg(x)=e^{-2t}\cos(x+3t)$。初期位置 $x=0$、$t=\pi/3$ では $-e^{-2\pi/3}$ である。減衰は分散 $4t$ に、位相の移動は平均移動 $3t$ に対応する。
:::

### 時刻依存の終端条件を展開する

時間依存の作用素 $\mathcal L_s$ では、固定した終端時刻 $T$ に対し $u(s,x)=P_{s,T}g(x)$ と置く。古典解が存在し、$u$ が有界なら、時刻 $s$ から出発した $u(r,X_r)$ に伊藤公式を適用して $u_s+\mathcal L_su=0$ の検証ができる。実際、ドリフトは零で全体が有界な局所マルチンゲールとなる。存在を確かめる場合は短時間の関係 $P_{s,T}=P_{s,s+h}P_{s+h,T}$ を用い、$P_{s+h,T}g$ の空間導関数と $s$ の微分を正当化した上で $h\downarrow0$ とする。$P_{s,T}$ を $P_{T-s}$ と置き換えると、時刻 $s$ で使う係数が失われる。

:::example[分散が時刻とともに増える場合]
$\dd X=\sqrt{2t}\,\dd W_t$ とし、終端値を $g(x)=\cos(2x)$ とする。時刻 $s$ から $T$ への増分の分散は $\int_s^T2r\,\dd r=T^2-s^2$ だから
$$
u(s,x)=e^{-2(T^2-s^2)}\cos(2x).
$$
$\partial_su=4su$、$\mathcal L_su=s u_{xx}=-4su$ で打ち消し合う。例えば $T=2$ のとき $s=1,x=0$ では $e^{-6}$ になる。長さ $1$ の区間でも分散は開始時刻で異なり、開始時刻 $s$ の偏微分を取る理由が見える。
:::

## Fokker–Planck 方程式

:::proposition[前向き方程式の弱形式]{#prop-sa-fokker-planck}
$X_t$ の分布を $\mu_t$ とする。任意の $f\in C_c^2$ について
$$
\int f\,\dd\mu_t-\int f\,\dd\mu_0=\int_0^t\int\mathcal L f\,\dd\mu_s\,\dd s.
$$
密度 $p(t,x)$ がある場合、分布の意味で
$$
\partial_t p=\mathcal L^*p
=-\sum_i\partial_i(b_i p)+\tfrac12\sum_{i,j}\partial_{ij}(a_{ij}p)
\label{eq-sa-fokker-planck}
$$
である。
:::

:::proof
前半は :eqref[eq-sa-dynkin] と Fubini で得られる。時間にもコンパクト台の試験関数を掛けて積分し、時間の一階微分を分布微分に移す。空間の一次項を一回、二次項を二回部分積分すれば形式随伴が得られる。境界項は試験関数のコンパクト台により消える。
:::

例えば $b=-\kappa x,a=\nu^2$ の OU で $p_*(x)\propto e^{-\kappa x^2/\nu^2}$ とすると $(\nu^2/2)p_*'=bp_*$ で流束が零である。従って $\mathcal L^*p_*=0$ である。密度の方程式では終端位置の変数を微分し、後ろ向き方程式では初期位置の変数を微分する。ドリフト項の微分が積の外にあるか内にあるかはこの違いに伴う。

### 前向き方程式の一意性を双対性で示す

典型的な一意性の仮定を明確にする。$b,\sigma$ は有界であり、任意の $g\in C_c^\infty$ について $u(s,x)=P_{t-s}g(x)$ が $0\le s\le t$ で古典的な後ろ向き方程式を満たすとする。さらに $u$ とその時間一階・空間二階までの導関数は有界とする。この条件は定数係数拡散の Gauss 核では、導関数を $g$ に作用させて確かめられる。

:::proposition[双対性による前向き方程式の一意性]{#prop-sa-forward-uniqueness}
上の仮定の下で、同じ初期確率測度を持つ弱連続な確率測度の族で、:ref[prop-sa-fokker-planck] の弱形式を満たすものは一意である。
:::

:::proof
まず時空で連続な $\partial_s v,\partial_i v,\partial_{ij}v$ を持ち、全時刻で共通の空間コンパクト台を持つ試験関数 $v(s,x)$ に弱形式を拡張する。時間を分割し、各小区間で $v$ の変化と測度の変化に分ける。前者は $\partial_s v$ の積分であり、後者には固定した左端の試験関数で弱形式を使う。微分の一様連続性と測度の弱連続性によって、和の極限は
$$
\mu_t(v(t,\cdot))-\mu_0(v(0,\cdot))
=\int_0^t\mu_s(\partial_s v+\mathcal L v)\,\dd s
$$
となる。$v(s,x)=\chi(x/R)u(s,x)$ とする。$\chi$ は単位球上 $1$、半径 $2$ の外で零の滑らかな切断関数である。$u_s+\mathcal L u=0$ なので右辺には $\mathcal L$ が $\chi$ に当たる項と、その一階微分と $u$ の一階微分の積だけが残る。係数と導関数の有界性により、それらは一様に $O(R^{-1})$ である。$R\to\infty$ で有界収束を使うと、$\mu_t(g)=\mu_0(P_tg)$ である。$C_c^\infty$ は有限 Borel 測度を分離するので一意である。
:::

例えば $\dd X=c\,\dd t+\nu\,\dd W$ で初期分布が $N(m_0,v_0)$ なら、正規増分との畳み込みで $N(m_0+ct,v_0+\nu^2t)$ が前向き方程式を解く。この双対性により、初期値を合わせた別の確率密度の族を同じ弱形式の下で作ることはできない。

### 確率流束と境界での保存

一変数で密度 $p$ と係数が滑らかなら、:eqref[eq-sa-fokker-planck] は $\partial_t p=-\partial_x J$、$J=bp-\tfrac12\partial_x(ap)$ と書ける。区間 $(l,r)$ の全質量の変化は $J(l)-J(r)$ である。反射境界で流出のないモデルなら流束を両端で零とし、殺すモデルなら境界値を零にする例が現れる。同じ内部の式でも、境界条件が異なると質量の時間変化が異なる。

:::example[反射と吸収の質量を比べる]
$(0,\pi)$ で Brown 運動を考える。反射の前向き方程式には Neumann 条件 $p_x(0)=p_x(\pi)=0$ を付ける。$p_0(x)=(1+\cos x/2)/\pi$ なら
$$
p(t,x)=\frac1\pi\big(1+\tfrac12e^{-t/2}\cos x\big)
$$
は非負で質量 $1$ を保つ。一方、吸収では $p_0(x)=\sin x/2$ とすれば $p(t,x)=e^{-t/2}\sin x/2$ は端点で零、質量は $e^{-t/2}$ である。$J=-p_x/2$ より吸収の場合の流束の差は $-e^{-t/2}/2$ であり、全質量の微分と一致する。吸収後も墓場状態を付けて確率を数えるなら、失われた質量は墓場状態の原子になる。
:::

### 雑音が一成分でも密度は二次元に広がる

:::example[積分 Brown 運動の前向き方程式]
$V_t=W_t$、$Q_t=\int_0^tV_s\,\dd s$ をともに零から始める。状態 $(q,v)$ で生成作用素は $v\partial_q+\tfrac12\partial_{vv}$ である。Fubini で $Q_t=\int_0^t(t-s)\,\dd W_s$ と書けるため、$(Q_t,V_t)$ は共分散
$$
\begin{pmatrix}
t^3/3&t^2/2\\
t^2/2&t
\end{pmatrix}
$$
の正規分布である。行列式は $t^4/12>0$ であり、$t>0$ では二変数の密度がある。固定した時刻の拡散行列では $q$ 方向に二階項がないが、ドリフトが過去の $V$ を時間積分するので独立な変動方向が加わる。前向き方程式は $p_t=-v p_q+p_{vv}/2$ であり、滑らかな正規密度は弱形式からこの式を古典的に満たす。
:::

:::example[決定論的な成分が残る場合]
同じ一次元ノイズでも、$\dd Q=2\,\dd t$、$\dd V=\dd W$ とすると $(Q_0,V_0)=(0,0)$ からの分布は $\delta_{2t}(\dd q)\otimes N(0,t)(\dd v)$ である。分布としては $\partial_t\mu=-2\partial_q\mu+\tfrac12\partial_{vv}\mu$ を解くが、二次元 Lebesgue 測度に対する密度は存在しない。マルチンゲール問題や前向き方程式の弱形式が密度の存在を仮定しない理由である。
:::

### 逆向きの双対性と時空密度

時刻に依存する弱形式が成り立ち、積分の交換に必要な有界性がある場合、$u(s,x)$ を試験関数に使うことで、$\partial_s u+\mathcal L_su=0$ なら $\int u(s,x)\mu_s(\dd x)$ は一定である。初期値 $x$ と終端変数 $y$ を持つ滑らかな核 $p(s,x;t,y)$ が存在するとき、この事実は次の二つの微分方程式と対応する。
$$
\partial_s p+\mathcal L_{s,x}p=0,\qquad
\partial_t p-\mathcal L^*_{t,y}p=0.
$$
第二式は初期測度を $\delta_x$ として導く。第一式は終端に任意の滑らかな $g(y)$ を掛けて積分すると後ろ向き方程式になり、微分と積分の交換が可能な正則性の下で試験関数の任意性から得る。例えば定数 $b$、$\sigma$ の正規核では $y-x-b(t-s)$ が核の変数であり、$\partial_x$ と $\partial_y$ の符号は反対である。分散も $t-s$ に比例するため、$s$ と $t$ の符号が逆になる点まで直接検証できる。

### 半直線で失われる質量を計算する

標準 Brown 運動を $x>0$ から始めて $0$ で殺す。$q_t(z)=(2\pi t)^{-1/2}e^{-z^2/(2t)}$ に対し、殺した核は
$$
p^+(t,x,y)=q_t(y-x)-q_t(y+x),\qquad y>0
$$
である。実際、$0$ に到達してから $y$ に終わる路を最初の到達後に符号反転すると、$-y$ に終わる路になる。強 Markov 性と Brown 増分の対称性により二つの事象の終点密度は等しく、元の密度からこの分を引けばよい。$p^+(t,x,0)=0$ であり、全質量の時間微分は $-\tfrac12\partial_y p^+(t,x,0)$ である。熱方程式を $y>0$ で積分し、無限遠の導関数が零になることを使った。微分すると
$$
\frac12\partial_y p^+(t,x,0)
=\frac{x}{\sqrt{2\pi t^3}}e^{-x^2/(2t)}
$$
なので、これが原点の初到達時刻の密度である。例えば $x=2,t=4$ では流出率は $e^{-1/2}/(4\sqrt{2\pi})$ となる。この計算では境界値が零であっても境界への流束は零にならない。

## Feynman–Kac の公式

:::theorem[終端値問題の Feynman–Kac 公式]{#thm-sa-feynman-kac}
有界な連続 $c,f$ と有界古典解 $u$ が
$$
u_t+\mathcal L u-cu+f=0,\qquad u(T,x)=g(x)
$$
を満たすとする。$t$ で $x$ から始めた拡散を用いると
$$
u(t,x)=\mathbb E_{t,x}\left[
e^{-\int_t^T c(X_r)\,\dd r}g(X_T)+
\int_t^T e^{-\int_t^s c(X_r)\,\dd r} f(X_s)\,\dd s
\right].
\label{eq-sa-feynman-kac}
$$
:::

:::proof
$D_s=e^{-\int_t^s c(X_r)\,\dd r}$ とする。積の公式で $D_su(s,X_s)$ のドリフトは $-D_s f(X_s)$ である。そこに $D_sf$ の積分を足すと局所マルチンゲールになり、有限時間のすべての項が有界なので真のマルチンゲールである。終端と初期の期待値を等置する。
:::

例えば Brown 運動で $c=3,f=2,g=1$ なら $u=e^{-3(T-t)}+\frac23(1-e^{-3(T-t)})$ である。これは空間一様の常微分方程式 $u_t-3u+2=0$ の解でもある。正の $c$ は遅い時刻の寄与を小さくする。

### 割引した作用素の合成と Duhamel 公式

時間に依らない有界連続な $c$ に対し、$Q_tg(x)=\mathbb E_x[D_tg(X_t)]$、$D_t=\exp(-\int_0^t c(X_s)\dd s)$ と定める。割引は区間を二つに分けると積になる。Markov 性で時刻 $s$ から再出発すれば $Q_{s+t}=Q_sQ_t$ である。$c\ge0$ の場合、$Q_t1\le1$ となり、保存的な半群とは異なり終端までの質量が減る。例えば独立な平均 $1$ の指数変数 $E$ を追加し、累積 $\int_0^t c(X_s)\dd s$ が $E$ を超えた時刻で殺せば、$X$ を条件とした時刻 $t$ の生存確率は $D_t$ である。

:::proposition[積分形の Feynman–Kac 方程式]{#prop-sa-pde-duhamel-kac}
$g,c$ を有界 Borel 関数、$V_t=Q_tg$ とする。すると
$$
V_t=P_tg-\int_0^t P_s(cV_{t-s})\,\dd s
$$
である。この式を満たす有界可測な族は有限時間上で一意である。$c\ge0$ なら $Q_t$ は正を保ち最大ノルムで縮小する。
:::

:::proof
固定時刻 $t$ で関数 $s\mapsto\exp(-\int_s^t c(X_r)\dd r)$ は絶対連続である。その導関数は $c(X_s)\exp(-\int_s^t c(X_r)\dd r)$ なので
$$
D_t=1-\int_0^t c(X_s)e^{-\int_s^t c(X_r)\dd r}\,\dd s.
$$
両辺に $g(X_t)$ を掛け、Fubini と時刻 $s$ の条件付き期待値により積分式を得る。二つの解の差を $D(t,x)$ とすると、$|D(t,x)|\le \|c\|_\infty\int_0^t P_s(|D(t-s,\cdot)|)(x)\dd s$。有限区間での一様上界を $M$ とし、この不等式に $M$ を入れる操作を $n$ 回反復すると $|D(t,x)|\le M\|c\|_\infty^n t^n/n!$ である。$n\to\infty$ で零になり一意である。正性と縮小性は $0<D_t\le1$ から従う。
:::

同じ仮定で時間に依らない有界な源 $f$ が入る場合、残り時間を $t$ として $V_t=Q_tg+\int_0^t Q_sf\dd s$ となる。Markov な割引の伝播と Fubini により、これは $V_t=P_tg+\int_0^tP_s(f-c V_{t-s})\,\dd s$ を満たす。積分式は空間微分を要求しない。滑らかな係数で古典解を求めるときは、先に $V_t$ の空間二階微分を正当化し、そこで初めて微分形へ移る。

:::example[三角関数の源と定数割引]
Brown 運動、$g(x)=\cos(2x)$、$f(x)=3\cos(2x)$、$c=1$ とする。$P_t\cos(2x)=e^{-2t}\cos(2x)$ なので、$Q_t$ の減衰は $e^{-3t}$ である。従って
$$
V_t(x)=e^{-3t}\cos(2x)+3\int_0^t e^{-3s}\cos(2x)\,\dd s=\cos(2x)
$$
となる。$-cV$ の減衰と熱作用素の減衰を源が合わせて補う。実際、$\tfrac12(\cos2x)''-\cos2x+3\cos2x=0$ である。
:::

:::proposition[滑らかなデータからの Feynman–Kac 古典解]{#prop-sa-pde-kac-smooth}
$b,\sigma,c,f,g$ を時刻に依存しない $C^3$ 関数とし、値と三階までの導関数を有界とする。$v(t,x)=Q_tg(x)+\int_0^t Q_s f(x)\,\dd s$ は有限時間区間上で $C^{1,2}$ であり、$v_t=\mathcal Lv-cv+f$、$v(0,x)=g(x)$ を満たす。
:::

:::proof
状態に $A_t=\int_0^t c(X_s)\,\dd s$ を加える。拡大した拡散のドリフトは $(b,c)$、拡散行列は $\sigma$ に零の行を加えたものであり、:ref[prop-sa-pde-backward-smooth] の仮定を満たす。有限時間 $T$ で $A_0=0$ とすれば $|A_t|\le T\|c\|_\infty$ である。この範囲で $e^{-a}g(x)$ と一致し、その外では $a$ について滑らかな切断を掛けた終端関数を使う。拡大した半群をこの関数に施し $a=0$ と置けば $Q_tg(x)$ が得られる。従って $Q_tg,Q_tf$ の空間二階までの微分は時空で連続有界であり、時間積分でも微分と積分を交換できる。

固定した $t$ に対し、半群の合成から $v(t+h)=Q_hv(t)+\int_0^hQ_s f\,\dd s$ である。割引の積に伊藤公式を使えば、任意の有界 $C^2$ 関数 $w$ で導関数が有界なとき、
$$
Q_h w(x)-w(x)=\int_0^h Q_s((\mathcal L-c)w)(x)\,\dd s
$$
となる。よって右差商は $(\mathcal L-c)v(t,x)+f(x)$ へ収束する。左差商では同じ式を $w=v(t-h)$ に適用する。$(\mathcal L-c)v$ が連続有界で $X_s^x\to x$ だから、その積分を $h$ で割っても同じ極限になる。したがって時間微分が存在し、右辺の時空連続性から古典解となる。
:::

### 有界でない非負ポテンシャル

$c$ の非負性があれば、経路に沿った積分が有限という仮定の下、割引そのものは有界である。

:::example[二次ポテンシャルでの割引]
Brown 運動に $c(x)=\lambda x^2$（$\lambda>0$）を取る。$V(t,x)=\mathbb E_x\exp(-\lambda\int_0^t X_s^2\dd s)$ を求めるには、$V_t=\tfrac12V_{xx}-\lambda x^2V$、$V(0,x)=1$ に対して $V=A(t)e^{-B(t)x^2/2}$ と置く。候補の偏微分から
$$
B'=2\lambda-B^2,\qquad A'/A=-B/2,\qquad B(0)=0,\quad A(0)=1
$$
を得る。$\omega=\sqrt{2\lambda}$ とすると $B=\omega\tanh(\omega t)$、$A=\cosh(\omega t)^{-1/2}$。候補は $0<V\le1$ である。確認のため、コンパクト退出まで止めて $D_sV(t-s,X_s)$ に伊藤公式を用いる。ドリフトが消え、全体が有界局所マルチンゲールなので、停止を外して上の期待値と一致する。非負の $c$ によって局所化後の一様有界性を確保している。例えば $\lambda=2,x=0,t=(\log 2)/2$ では $\cosh(2t)=5/4$ なので $V=2/\sqrt5$ である。
:::

:::example[時刻だけで変わる割引率]
終端 $T=2$、終端値 $1$、源項零、殺し率 $V(s,x)=2s$ とする。Feynman–Kac の期待値の被積分確率変数は $\exp(-\int_t^2 2s\dd s)=e^{t^2-4}$ に等しく、状態に依存しない。従って $u(t,x)=e^{t^2-4}$ である。$\partial_t u=2tu$、$\mathcal L u=0$ なので $(\partial_t+\mathcal L-V)u=0$、$u(2,x)=1$ を直接確かめられる。時刻 $0$ からの係数は $e^{-4}$ であり、殺し率の終端値 $4$ を全期間に用いた $e^{-8}$ と異なる。必要なのは各時刻の率の積分である。
:::

## 退出時刻と境界値問題

:::proposition[退出で止めた表現]{#prop-sa-exit-pde}
$D$ を有界開集合、$\tau$ をその退出時刻とする。$\tau<\infty$ a.s.、$\mathbb E_x\tau<\infty$ とする。$c,f$ は有界連続で $c\ge0$ とする。$u\in C^2(D)\cap C(\overline D)$ が $\mathcal L u-cu=-f$、$u|_{\partial D}=g$ を満たすなら
$$
u(x)=\mathbb E_x\left[D_\tau g(X_\tau)+\int_0^\tau D_s f(X_s)\,\dd s\right],
\qquad D_s=e^{-\int_0^s c(X_r)\,\dd r}.
$$
:::

:::proof
コンパクトな部分領域の退出と時刻 $n$ の小さい方で止めて、:ref[thm-sa-feynman-kac] の積の計算を使う。連続路が $D$ に留まる任意のコンパクトな時間範囲は、あるコンパクト部分領域に含まれる。従って部分領域を増やすと停止時刻は $\tau$ に近づく。$D_s\le1$ と $u$ の有界性で終端に有界収束、積分項には $\norm f_\infty\tau$ による優収束を使う。
:::

例えば $(0,4)$ で Brown 運動の平均退出時間は $\frac12u''=-1$、$u(0)=u(4)=0$ の解 $x(4-x)$ である。先に $\tau\wedge n$ でこの非負解を用いれば $\mathbb E(\tau\wedge n)\le x(4-x)$ となり、平均退出時刻有限の仮定を確認できる。

### 割引退出時刻の変換

:::example[有限区間からの退出時刻の Laplace 変換]
標準 Brown 運動が $(-a,a)$ から出る時刻を $\tau$ とする。$\lambda>0$ に対し $u(x)=\mathbb E_x e^{-\lambda\tau}$ は $u''/2=\lambda u$、$u(-a)=u(a)=1$ に対応する。$k=\sqrt{2\lambda}$ とすれば
$$
u(x)=\frac{\cosh(kx)}{\cosh(ka)}.
$$
候補は有界なので、指数付き伊藤公式と $\tau<\infty$ で検証できる。$\lambda\downarrow0$ の展開は $1-\lambda(a^2-x^2)+O(\lambda^2)$ である。$\tau$ の平均が $a^2-x^2$ という Dynkin の結果を零での微分からも確認できる。例えば $a=2,x=0,\lambda=1/2$ なら $\mathbb E_0e^{-\tau/2}=1/\cosh 2$ である。
:::

:::example[退出時刻と出口を同時に重み付けする]
$(0,a)$ では上端の出口だけに重みを付けて $u(x)=\mathbb E_x[e^{-\lambda\tau};X_\tau=a]$ とする。境界値は $u(0)=0,u(a)=1$ であり、$u(x)=\sinh(kx)/\sinh(ka)$ となる。$\lambda\downarrow0$ では $x/a$ に戻る。$a=2$、$x=1$ なら反射対称性から上端の確率は $1/2$ で、条件付き変換は $2\sinh(k)/\sinh(2k)=1/\cosh(k)$。中央出発の場合、退出方向と退出時間が独立であることもこの対称性から分かる。実際、位置の反転で退出時間を保存し上側と下側の事象を交換する。
:::

### 殺した半群と固有関数

$P_t^Df(x)=\mathbb E_x[f(X_t);t<\tau_D]$ とすると、Markov 性で $P_{s+t}^D=P_s^DP_t^D$ である。$P_t^D1(x)$ は生存確率であり $1$ 以下である。$u\in C^2(D)\cap C(\overline D)$ が $\mathcal Lu=-\lambda u$、$u|_{\partial D}=0$ を満たすなら、止めた $e^{\lambda t}u(X_t)$ に伊藤公式を適用して $P_t^Du=e^{-\lambda t}u$ を得る。有限の $t$ で $u$ は有界で、退出した路では端の値が零だからである。

:::example[第一固有関数を初期重みに選ぶ]
Brown 運動の $(0,\pi)$ では $u(x)=\sin x$ が $\lambda=1/2$ の固有関数である。始点 $x=\pi/2$ から時刻 $t$ まで生き残った路の重み $\mathbb E_x[\sin X_t;t<\tau]$ は $e^{-t/2}$ である。$0<\sin x\le1$ を用いれば $P_x(t<\tau)\ge e^{-t/2}$。時間が十分長いときの指数の大きさをこの一個の固有関数から読み取れる。開始分布を密度 $\sin x/2$ にしたときは、:eqref[eq-sa-fokker-planck] の吸収例と対応して、条件付きの生存分布は同じ密度のままで生存確率だけが $e^{-t/2}$ になる。
:::

### 反射境界の Neumann 条件

半直線では、退出して境界値を読む代わりに、反射の押し戻しを境界の積分に使える。$R=x+W+K$ を $x\ge0$ からの反射 Brown 運動とし、$K$ は $R=0$ の時刻だけで増加する。$c>0$、有界連続な $f$、定数 $h$ に対し、$u$ が有界で、内部で $C^2$、原点まで二階導関数が連続に延び、
$$
\tfrac12 u''-cu=-f \quad (x>0),\qquad u'(0)=h
$$
を満たすとする。$u'(0)$ は内部向きの導関数である。$e^{-ct}u(R_t)$ の積の公式では、実時間のドリフトは $-e^{-ct}f(R_t)\,\dd t$、境界項は $h e^{-ct}\,\dd K_t$ である。有限時間にコンパクト領域の退出を併用して期待値を取れば、停止した確率項が消える。停止を外すときは $u,f$ の有界性と $\mathbb E K_t<\infty$ を使える。反射問題の明示式より $\mathbb E_x K_t\le\sqrt{2t/\pi}$ であり、部分積分で $\mathbb E_x\int_0^\infty e^{-cs}\dd K_s<\infty$ となる。時刻を無限へ送ると
$$
u(x)=\mathbb E_x\int_0^\infty e^{-cs}f(R_s)\,\dd s
     - h\,\mathbb E_x\int_0^\infty e^{-cs}\,\dd K_s
$$
を得る。したがって $h=0$ の同次 Neumann 条件では内部の源だけが積分に残る。この表現は上のクラスでの一意性も与える。

例えば $f=0$ のとき $u(x)=-h e^{-\sqrt{2c}x}/\sqrt{2c}$ は微分方程式と境界条件を満たす。$h=1$ を入れると、割引した押し戻しの期待値は $e^{-\sqrt{2c}x}/\sqrt{2c}$ である。原点で始めれば $c=1/2$ で値は $1$ となる。右局所時間で表す場合、:ref[prop-sa-reflecting-transition] の反射表示と局所時間の規格化により $L_t^0(R)=2K_t$ なので、境界積分の係数を半分にする。

### 非負性からの比較原理

同じ作用素と $c\ge0$ を用いる二つの境界値解 $u_1,u_2$ について、境界値と源がともに順序付けられていれば :ref[prop-sa-exit-pde] の被積分関数も同じ順序になる。従って $u_1\le u_2$ が従い、データが同じなら一意である。例えば境界値が $0$、源が $1$ の解は、非負割引 $c$ を追加すると小さくなる。さらに $c=\lambda>0$ なら
$$
0\le\mathbb E_x\int_0^\tau e^{-\lambda s}\dd s
=\frac{1-\mathbb E_x e^{-\lambda\tau}}{\lambda}\le1/\lambda.
$$
$(-2,2)$ の中央、$\lambda=1/2$ では値は $2(1-1/\cosh2)$ であり、無割引の平均退出時間 $4$ より小さい。右辺に退出時刻の変換を代入すれば、連立の境界値問題を解き直さずに占有時間の割引を計算できる。

## 不変測度とエルゴード性

確率測度 $\pi$ が $\int P_tf\,\dd\pi=\int f\,\dd\pi$ を全ての $t$ と有界 Borel $f$ で満たすとき :term[不変測度]{en="invariant measure" yomi="ふへんそくど"} という。$\mathcal L^*p=0$ は定常性の候補を与える。候補が本当に不変かどうかは推移半群または前向き方程式の一意性で確認する。

:::proposition[OU の不変分布と時間平均]{#prop-sa-ou-ergodic}
$\kappa>0,\nu>0$ の OU の不変確率測度は $\pi=N(m,\nu^2/(2\kappa))$ だけである。$X_0\sim\pi$ のとき、有界 Lipschitz $f$ の時間平均は $\pi(f)$ に $L^2$ と a.s. で収束する。
:::

:::proof
:eqref[eq-sa-ou] で正規初期値と独立ノイズの和を取れば平均と分散が保たれる。任意の固定初期値から正規分布が $\pi$ へ弱収束するので、もし別の不変測度があれば、この収束を有界連続試験関数で積分し、優収束により $\pi$ に一致する。

同じ Brown 運動で $X^x$ と独立な $\pi$ 初期値 $Y_0$ の OU を結合すると、$\abs{X_t^x-Y_t}=e^{-\kappa t}|x-Y_0|$ である。従って $|P_tf(x)-\pi(f)|\le\mathrm{Lip}(f)e^{-\kappa t}\mathbb E|x-Y_0|$。定常初期値で Markov 性を使うと $|\operatorname{Cov}(f(X_s),f(X_t))|\le Ce^{-\kappa(t-s)}$（$s\le t$）である。二重積分で時間平均の分散は $2C/(\kappa T)$ 以下となる。$T=n^2$ に対して Chebyshev と Borel–Cantelli を使えば概収束する。有界な $f$ では $n^2$ から $(n+1)^2$ の間での時間平均の差も一様に零となる。
:::

### 時間平均した測度から不変測度を作る

長時間分布の極限が存在しなくても、時刻を平均した測度が緊密なら不変測度を構成できる。これは定常状態の候補を作る手続きであり、その一意性や経路の時間平均の結論は別に示す。

:::proposition[時間平均の緊密性による不変測度の存在]{#prop-sa-pde-krylov-bogoliubov}
$P_t$ が有界連続関数を有界連続関数に移すとする。ある $x$ に対し $\nu_T=T^{-1}\int_0^T p_t(x,\cdot)\dd t$ が $T\ge1$ で緊密なら、$T_n\to\infty$ として取った弱収束部分列の極限は不変測度である。
:::

:::proof
緊密性と :ref[thm-prob-prokhorov]{full} により $\nu_{T_n}\Rightarrow\pi$ となる $T_n\uparrow\infty$ の部分列を選べる。有界連続 $f$ と $s\ge0$ に対し
$$
\nu_T(P_sf)-\nu_T(f)
=\frac1T\left(\int_T^{T+s}P_tf(x)\,\dd t-\int_0^sP_tf(x)\,\dd t\right).
$$
絶対値は $2s\|f\|_\infty/T$ 以下である。$P_sf$ が有界連続なので両辺を弱収束で極限に移すと $\pi(P_sf)=\pi(f)$。有界連続関数は確率測度を分離するため、測度として不変であり有界 Borel にも等式が従う。
:::

:::example[Lyapunov 評価で尾を抑える]
$\dd X=-X^3\,\dd t+\sqrt2\,\dd W$ を考える。$V(x)=1+x^2$ に対し $\mathcal LV=2-2x^4\le4-V(x)$ である。実際 $2-2x^4+V=3+x^2-2x^4\le25/8<4$ である。指数付きの停止した伊藤公式から $\mathbb EV(X_t)\le e^{-t}V(x)+4(1-e^{-t})$ が従う。非爆発性は :ref[prop-sa-sde-lyapunov] で確認できる。時刻平均の測度でもこの上界を使えば、$|y|>R$ の質量は $\max(V(x),4)/(1+R^2)$ 以下となる。推移関数の空間連続性は解をコンパクトまで局所化して大域 Lipschitz の初期値連続性を使い、退出確率は $V$ によってコンパクトな初期値上で一様に抑えられる。従って命題の仮定が成り立ち、不変測度を少なくとも一つ得る。
:::

### 詳細釣り合いとエネルギー

$\mathcal L=\tfrac12\Delta-\tfrac12\grad U\cdot\grad$ とし、$e^{-U}$ が可積分とする。コンパクト台の試験関数 $f,g$ に対し
$$
\int g\,\mathcal L f\,e^{-U}\,\dd x=-\frac12\int \grad g\cdot\grad f\,e^{-U}\,\dd x.
$$
実際 $\mathcal Lf=e^U\operatorname{div}(e^{-U}\grad f)/2$ なので部分積分で出る。この対称性は、定常分布候補 $e^{-U}/\int e^{-U}$ の下で生成作用素のエネルギーが非正であることを示す。また $\int\mathcal Lf\,e^{-U}\dd x=0$ が従うため候補は定常前向き方程式を満たす。ただしこの微分の議論から半群の不変性へ進むには、弱前向き方程式の一意性などを確認する。例えば有界滑らかな係数と有界導関数の条件では :ref[prop-sa-forward-uniqueness] と :ref[prop-sa-pde-backward-smooth] が使える。

:::example[有界な復元力が持つ定常密度]
$U(x)=2\log\cosh x$、$\dd X=-\tanh X\,\dd t+\dd W$ とする。ドリフトと三階までの導関数は有界である。候補 $\pi(\dd x)=\tfrac12\operatorname{sech}^2 x\,\dd x$ は $\int\operatorname{sech}^2 x\,\dd x=\tanh x$ により正規化されている。$\pi$ を初期分布とする定数の測度族は定常前向き方程式を解き、双対性による一意性からこの拡散の不変測度である。確率 $P_\pi(-\log 2\le X\le \log 2)$ は $\tanh(\log 2)=(4-1)/(4+1)=3/5$ である。定常分布を正規と仮定する必要はない。
:::

### 収縮から定常分布の一意性へ

共通 Brown 運動で二つの解を結合する。各 $x,y$ で
$$
2(x-y)\cdot(b(x)-b(y))+\norm{\sigma(x)-\sigma(y)}_{\mathrm F}^2
\le -\kappa|x-y|^2
$$
とする。差の二乗の伊藤公式に $e^{\kappa t}$ を掛け、局所化したマルチンゲール項の平均を消す。非負性と Fatou から $\mathbb E|X_t^x-X_t^y|^2\le e^{-\kappa t}|x-y|^2$ である。従って有界 Lipschitz $f$ に対して
$$
|P_tf(x)-P_tf(y)|\le\operatorname{Lip}(f)e^{-\kappa t/2}|x-y|.
$$
二つの不変分布を独立に初期値に取れば、二つの期待値の差は上式右辺と $2\|f\|_\infty$ の小さい方で抑えられる。固定した初期点対でこの上界は零へ行き、有界収束で二つの測度は有界 Lipschitz 関数上で一致する。閉集合の指示関数を Lipschitz な距離関数から近似すれば測度として一致する。初期分布の二次モーメントを仮定せずに一意性が得られる。

:::example[非定常点から定常時間平均へ]
$X$ を $x$ からの OU、$Y$ を定常 OU とし、同じ Brown 運動で動かすと $X_t-Y_t=e^{-\kappa t}(x-Y_0)$ である。$f$ が Lipschitz なら
$$
\left|\frac1T\int_0^T(f(X_t)-f(Y_t))\,\dd t\right|
\le\operatorname{Lip}(f)|x-Y_0|\frac{1-e^{-\kappa T}}{\kappa T}.
$$
右辺は a.s. 零へ行くため、:ref[prop-sa-ou-ergodic] の定常時間平均は任意の $x$ からも同じ極限を持つ。$\kappa=2$ で $|x-Y_0|=3$ の標本なら、時間平均差の上界は $\operatorname{Lip}(f)\,3/(2T)$ である。初期分布の影響が短時間に集中する理由はこの結合にある。
:::

### Poisson 方程式で時間積分を分解する

$\mathcal Lu=f-\pi(f)$ を満たす $u$ があり、伊藤公式の確率積分が定義できるとする。すると
$$
\int_0^T(f(X_t)-\pi(f))\,\dd t
=u(X_T)-u(X_0)-\int_0^T\grad u(X_t)^\top\sigma(X_t)\,\dd W_t.
$$
$u(X_T)-u(X_0)$ が時間に一様な二次モーメントを持てば、時間平均ではこの部分は $L^2$ で $T^{-1}$ の大きさになる。積分の二次変分密度 $\norm{\grad u^\top\sigma}^2$ の定常平均が有限なら、確率項は $L^2$ で $\sqrt T$ の大きさであり、平均を取れば $T^{-1/2}$ の大きさになる。この方法は時間相関を二重積分する方法とは異なる計算で同じ尺度を与える。

:::example[OU の積分を確率項へ移す]
定常 OU $\dd X=-\kappa X\,\dd t+\nu\,\dd W$ で $f(x)=x$、$u(x)=-x/\kappa$ を取ると $\mathcal Lu=x$ である。従って
$$
T^{-1/2}\int_0^T X_s\,\dd s
=\frac{X_0-X_T}{\kappa\sqrt T}+\frac{\nu W_T}{\kappa\sqrt T}.
$$
第一項は定常な二次積率から $L^2$ で零に近づき、第二項は各 $T$ で $N(0,\nu^2/\kappa^2)$ に従う。よって全体はその正規分布へ弱収束する。例えば $\kappa=2,\nu=3$ なら極限分散は $9/4$ である。二つの項の独立性は必要ない。第一項が確率収束で零になることと、第二項の分布が一定であることを使った。
:::

## 小雑音極限と Freidlin–Wentzell の大偏差

$\dd X^\eps=b(X^\eps)\,\dd t+\sqrt{\eps}\,\sigma(X^\eps)\,\dd W$ とする。典型的な路は $\dot x=b(x)$ に近い。:ref[thm-sa-bdg] と Lipschitz 評価を組み合わせると $\mathbb E\sup_{t\le T}\norm{X^\eps_t-x(t)}^2\le C_T\eps$ が得られる。実際、ドリフト差の二乗最大値を差の時間積分で抑え、確率項を $\eps$ 倍の二乗積分で抑える。線形成長の積率評価は $0<\eps\le1$ で一様なので Grönwall を適用できる。

加法的な定数 $\sigma$ の場合、大偏差原理は :ref[thm-prob-contraction]{full} と Brown 路の Schilder の定理から導ける。Schilder の定理では $\sqrt\eps W$ の一様ノルムでの良い率関数は、絶対連続な $h$ で $\frac12\int_0^T\norm{\dot h}^2\,\dd t$、その他で $\infty$ である。格子近似から路全体へ移す証明を後述する。入力 $h$ に対し $x_h=x_0+\int b(x_h)\,\dd t+\sigma h$ を解く写像は Grönwall により一様ノルムで連続である。よって縮約により率関数は
$$
I_{0T}(\phi)=\inf\left\{\tfrac12\int_0^T\norm{v(t)}^2\,\dd t:
\dot\phi=b(\phi)+\sigma v,\quad \phi(0)=x_0\right\}
\label{eq-sa-action}
$$
となる。非絶対連続な $\phi$ や空の制約は $\infty$ とする。状態依存 $\sigma$ を含む有界大域 Lipschitz 係数でも、同じ制御表示で大偏差原理が成り立つ。以下では時間格子で係数を固定した近似を使って証明する。

:::example[時刻を固定した最小作用]
$b=0,\sigma=1,T=2$ で $0$ から $3$ へ向かうとき、Cauchy–Schwarz により $\int_0^2\dot\phi^2\ge3^2/2$ である。等号は一定速度 $3/2$ で達成され、最小作用は $9/4$ である。これは終端の $N(0,2\eps)$ の指数部と一致する。
:::

### Brown 路の Schilder 定理を格子から導く

終端の正規分布の指数を知るだけでは、途中の広がりを含む路全体の確率を求められない。そこで格子の独立増分から多時刻の指数を得た後、格子間の振れを指数の尺度で制御する。以下の証明は加法ノイズの拡散に必要な路空間での主張を与える。

:::theorem[一様位相での Schilder の定理]{#thm-sa-pde-schilder}
$H^\eps=\sqrt\eps W$ を $C_0([0,T];\R^m)$、すなわち時刻零で零の連続路の空間に入れ、一様ノルムを使う。$h$ が絶対連続で $\dot h\in L^2$ なら $J(h)=\tfrac12\int_0^T|\dot h|^2\,\dd t$、それ以外では $J(h)=\infty$ とする。任意の閉集合 $F$、開集合 $G$ に対し
$$
\begin{aligned}
\limsup_{\eps\downarrow0}\eps\log P(H^\eps\in F)&\le-\inf_{h\in F}J(h),\\
\liminf_{\eps\downarrow0}\eps\log P(H^\eps\in G)&\ge-\inf_{h\in G}J(h).
\end{aligned}
$$
さらに $\{J\le K\}$ は各有限 $K$ でコンパクトである。
:::

:::proof
**格子の分布。** $t_i=iT/n$ とし、$\Pi_n h$ をこの格子での折れ線補間とする。$\Pi_n H^\eps$ の分布を増分 $z_i=h(t_i)-h(t_{i-1})$ で表すと、密度は
$$
(2\pi \eps T/n)^{-mn/2}
\exp\left(-\frac{1}{2\eps}\sum_{i=1}^n\frac{|z_i|^2}{T/n}\right)
$$
である。指数の分子 $J_n(z)$ は折れ線の $J$ に等しい。固定した $n$ で、閉集合上の上限は $-\inf J_n$、開集合上の下限は $-\inf J_n$ になる。上限はまずコンパクトな部分で密度の上界と体積を使う。外側の尾は、例えば $\mathbb E\exp(J_n/(2\eps))=2^{mn/2}$ を用いて、$\{J_n>M\}$ の確率を $2^{mn/2}e^{-M/(2\eps)}$ で抑える。固定した $M$ で密度の正規化因子は $\eps\log$ の極限で消え、次に $M\to\infty$ とする。下限は開集合の任意の点の周りの小さい球で、$J_n$ の連続性と正の体積を使う。増分の座標と折れ線の座標は互いに連続な線形写像で対応する。

**補間誤差。** 各格子間で $h-\Pi_nh$ のノルムは左端からの最大増分ノルムの二倍以下である。一成分で反射原理と正規裾の評価を用い、$m$ 成分と $n$ 区間の和を取ると、
$$
P\big(\norm{H^\eps-\Pi_nH^\eps}_\infty>\delta\big)
\le 4mn\exp\left(-\frac{n\delta^2}{8mT\eps}\right).
$$
実際、ベクトルノルムが $\delta/(2\sqrt\eps)$ を超えるなら、ある成分で絶対値が $\delta/(2\sqrt{m\eps})$ を超える。各区間長は $T/n$ である。

**率関数のコンパクト性。** $\{J\le K\}$ の路は
$$
|h(t)-h(s)|\le\sqrt{2K}|t-s|^{1/2},\qquad |h(t)|\le\sqrt{2KT}
$$
を満たす。:ref[thm-anl-ascoli-interval]{full} を有限個の座標に用いれば、一様収束する部分列がある。さらに導関数を $L^2$ の弱コンパクト性で部分列に取り、弱極限を $v$ とする。各 $t$ で $\int_0^t\dot h_n$ は指示関数との内積なので $\int_0^t v$ へ収束する。よって一様極限は $h(t)=\int_0^t v$ である。弱下半連続性により $J(h)\le K$ となり、集合はコンパクトである。

**上限と下限。** 閉集合 $F$ に対し $\delta$-近傍の閉包を $F^\delta=\{h:\operatorname{dist}(h,F)\le\delta\}$ と書く。$H^\eps\in F$ のとき、補間誤差が $\delta$ 以下なら補間路は $F^\delta$ に属する。従って極限上限は $-\inf_{F^\delta}J$ と $-n\delta^2/(8mT)$ の大きい方以下である。格子の率関数は折れ線での $J$ なので、折れ線に制限した下限を全路の下限へ弱めてよい。先に $n\to\infty$、次に $\delta\downarrow0$ とする。水準集合のコンパクト性により $\inf_{F^\delta}J\uparrow\inf_F J$ である。もし有限値までしか増えないなら、ほぼ最小化する路の収束部分列が $F$ に入り、下半連続性でその値が達成されるからである。

開集合 $G$ では $J(h)<\infty$ の一つの $h\in G$ を固定し、半径 $3\delta$ の球が $G$ に含まれるように取る。$\Pi_nh$ は一様に $h$ へ収束し、各区間の平均速度への Cauchy–Schwarz によって $J(\Pi_nh)\le J(h)$ である。十分大きい固定 $n$ で $\Pi_n H^\eps$ が $\Pi_nh$ の $\delta$ 球に入り、誤差が $\delta$ 以下なら $H^\eps\in G$ である。前の事象の確率の指数下限は $-J(h)$ 以上で、後の誤差確率の指数上限は $-n\delta^2/(8mT)$ である。この値が $-J(h)-1$ より小さくなるように $n$ を選ぶと、差し引いた確率の指数下限も $-J(h)$ 以上である。$h$ について下限を取ると主張を得る。
:::

加法ノイズで $P_\eps$ を $X^\eps$ の路の法則とする。$H^\eps$ から解への写像は :eqref[eq-sa-action] の直前で説明した Grönwall 評価で連続なので、:ref[thm-prob-contraction] と :ref[thm-sa-pde-schilder] により $P_\eps$ の率関数は :eqref[eq-sa-action] になる。$F$ と $G$ の下限が一致する場合は、対応する事象の確率が指数の尺度で明確に決まる。境界条件にちょうど接する路では閉集合と内部の下限が異なり得るので、いきなり指数の等号を書けない。

### 状態依存の雑音を含む Freidlin–Wentzell の原理

:::theorem[有界 Lipschitz 係数の経路大偏差]{#thm-sa-pde-freidlin-wentzell}
$b,\sigma$ を有界大域 Lipschitz とし、初期点 $x_0$ を固定する。$\dd X^\eps=b(X^\eps)\,\dd t+\sqrt\eps\,\sigma(X^\eps)\,\dd W$ の $[0,T]$ 上の路の法則は、一様位相において速度 $\eps^{-1}$ の大偏差原理を満たす。良い率関数は :eqref[eq-sa-action] の $I_{0T}$ である。拡散行列に正定値性を要求しない。
:::

:::proof
**格子で係数を固定する。** 格子幅を $h=T/n$、時刻 $s$ の左端を $\eta_h(s)$ とする。$Y^{\eps,h}$ は各区間で $b,\sigma$ を $Y^{\eps,h}_{\eta_h(s)}$ に固定した Euler 過程である。同じ $W$ で $X^\eps$ と結合する。$Y^{\eps,h}=F_h(\sqrt\eps W)$ と書ける。$F_h$ は各区間で左端の値にドリフトの直線と入力路の増分を足す写像なので、有限個の区間について帰納すれば一様ノルムで連続である。

任意の $r>0$ を取り、$Y^{\eps,h}$ が格子の左端から距離 $r$ に達する時刻で止める。各区間の右端も、そこで左端を更新する前の値を調べる。$h$ が十分小さければ格子区間のドリフト増分は $r/2$ 以下である。係数は有界なので、雑音の最大値を標準 Brown 運動の成分の最大値で抑え、反射原理と和の上界を使う。時刻 $T$ までに止まる確率は $C n\exp(-r^2/(C\eps h))$ 以下である。$C$ は $r,h,\eps$ によらない。

**二つの近似の誤差。** $D=X^\eps-Y^{\eps,h}$ とする。上の停止前では $D$ のドリフト $\beta$ と、$\sqrt\eps$ を除いたノイズ係数 $\Gamma$ は $|\beta|+\|\Gamma\|_{\mathrm F}\le L(|D|+r)$ を満たす。$F(z)=\log(1+|z|^2/r^2)$ に伊藤公式を使う。$|\grad F|\le 2|z|/(r^2+|z|^2)$、$\|\nabla^2F\|\le6/(r^2+|z|^2)$ なので、ドリフトは上から $C$、マルチンゲールの二次変分密度は $C\eps$ で抑えられる。ここで $0<\eps\le1$ とする。$D_0=0$ だから、格子増分による停止前に $|D|$ が $\delta$ に達するには、このマルチンゲールの上限が $\log(1+\delta^2/r^2)-CT$ 以上にならなければならない。停止した指数マルチンゲールより、二次変分が $v$ 以下のマルチンゲールが高さ $a>0$ に達する確率は $e^{-a^2/(2v)}$ 以下である。これを $v=C\eps T$ に使い、格子増分による停止確率と足し合わせる。まず $h\downarrow0$、次に $r\downarrow0$ とすれば
$$
\lim_{h\downarrow0}\limsup_{\eps\downarrow0}
\eps\log P(\|X^\eps-Y^{\eps,h}\|_\infty>\delta)=-\infty
\qquad(\delta>0)
$$
を得る。

**制御方程式の極限。** $J(w)<\infty$ なら $v=\dot w$ として、$\dot\phi=b(\phi)+\sigma(\phi)v$、$\phi(0)=x_0$ の解を $F_0(w)$ と書く。この方程式は可積分な Lipschitz 係数 $L(1+|v|)$ を持つので、積分方程式の反復と Grönwall で一意な絶対連続解を持つ。$\{J\le M\}$ 上で $F_h(w)$ の格子内増分は有界な係数と Cauchy–Schwarz により $C(h+\sqrt{2Mh})$ 以下である。二つの制御路の差は、この値と差そのものを $L(1+|v|)$ で重み付けした時間積分で抑えられる。$\int_0^T|v|\le\sqrt{2TM}$ なので、Grönwall で $F_h\to F_0$ が $\{J\le M\}$ 上で一様に成り立つ。従って $F_0$ は各水準集合上で連続である。この集合のコンパクト性と $J$ の下半連続性により、$I(\phi)=\inf_{F_0(w)=\phi}J(w)$ の有限な下限は達成され、$\{I\le M\}=F_0(\{J\le M\})$ はコンパクトである。

**確率の上限と下限。** 閉集合 $F$ と $0\le M<\inf_F I$ を取る。コンパクト集合 $F_0(\{J\le M\})$ と $F$ は交わらないので、それらの距離の半分より小さい $\delta>0$ を選べる。$h$ が十分小さいと $F_h(\{J\le M\})$ は $F$ の閉 $\delta$ 近傍に当たらない。Schilder の上限を $F_h$ によるこの近傍の逆像に使うと、$Y^{\eps,h}$ がそこへ入る確率の指数上限は $-M$ 以下である。さらに $h$ を小さく取り、誤差が $\delta$ を超える確率の指数上限も $-M$ 以下にする。和の上界を取り、$M\uparrow\inf_F I$ とすれば閉集合の上限を得る。$\inf_F I=0$ では確率の上界 $1$ を使う。

開集合 $G$ について $F_0(w)\in G$ かつ $J(w)<\infty$ を満たす $w$ を選ぶ。$F_0(w)$ の半径 $3\delta$ の球を $G$ に含める。小さい固定 $h$ なら $F_h(w)$ は $F_0(w)$ から距離 $\delta$ 以内である。$Y^{\eps,h}$ が $F_h(w)$ の開 $\delta$ 球へ入る確率は、Schilder の下限と $F_h$ の連続性から指数下限 $-J(w)$ 以上を持つ。誤差が $\delta$ を超える確率の指数上限を $-J(w)-1$ 未満にしておけば、差し引いてもこの下限は保たれ、その事象では $X^\eps\in G$ である。すべての $w$ について下限を取ると開集合下限を得る。最後に $v=\dot w$ と書けば $I=I_{0T}$ である。
:::

### 線形系の最小エネルギー

:::proposition[終端制約と可制御 Gram 行列]{#prop-sa-pde-linear-action}
$\dot\phi=A\phi+\Sigma v$、$\phi(0)=x$ を考える。$Q_T=\int_0^T e^{A(T-s)}\Sigma\Sigma^\top e^{A^\top(T-s)}\,\dd s$ と置く。終端 $\phi(T)=y$ を指定すると、$z=y-e^{AT}x$ が $\operatorname{Im}Q_T$ にある場合の最小作用は $\tfrac12 z^\top Q_T^+z$ であり、それ以外では無限である。$Q_T^+$ は Moore–Penrose 逆行列である。
:::

:::proof
$Lv=\int_0^T e^{A(T-s)}\Sigma v(s)\dd s$ と置く。終端制約は $Lv=z$ であり、$L^*w$ は $s\mapsto\Sigma^\top e^{A^\top(T-s)}w$ である。よって $LL^*=Q_T$。$w^\top Q_Tw=\|L^*w\|_{L^2}^2$ なので $\ker Q_T=\ker L^*$ である。値域が有限次元で閉じているため、その直交補を取ると $\operatorname{Im}L=\operatorname{Im}Q_T$ になる。到達可能なら $v_*=L^*Q_T^+z$ が制約を満たす。他の解との差は $\ker L$ に入り、$v_*$ と直交するので $\|v\|^2=\|v_*\|^2+\|v-v_*\|^2$。これで最小値とその制御を得る。
:::

:::example[小雑音 OU の終端]
$\dot\phi=-2\phi+3v$、$x=0$ とする。$Q_T=9(1-e^{-4T})/4$ なので、終端 $y=3$ への最小作用は $2/(1-e^{-4T})$ である。ノイズを加えた $X_T^\eps$ の分散は $\eps Q_T$ なので、正規密度の指数と一致する。$T\to\infty$ で作用は $2$ に下がる。一方 $T\downarrow0$ では $1/(2T)$ の大きさで発散する。雑音を使う時間が短いほど大きい制御エネルギーが要る。
:::

:::example[雑音が入らない成分への最小作用]
$\dot q=v$、$\dot z=q$、初期値 $(q,z)=(0,0)$ とする。$T=1$ の終端 Gram 行列は座標順を $(q,z)$ として
$$
Q_1=\begin{pmatrix}1&1/2\\1/2&1/3\end{pmatrix},\qquad
Q_1^{-1}=\begin{pmatrix}4&-6\\-6&12\end{pmatrix}.
$$
従って終端 $(0,1)$ の最小作用は $6$ である。制御は $v_*(s)=-6+12(1-s)$ であり、$\int_0^1 v_*\,\dd s=0$、$\int_0^1(1-s)v_*(s)\dd s=1$ を満たす。二つ目の成分に雑音を直接入れなくても、ドリフトを介して有限コストで到達できる。
:::

### 勾配系の準ポテンシャル

固定した初期点と終端点の間で、最小作用を時間にもわたって下限したものを :term[準ポテンシャル]{en="quasipotential" yomi="じゅんぽてんしゃる"} という。$b=-\grad U$、$\sigma=I$ の場合には任意の絶対連続路に対し
$$
\frac12\int_0^T|\dot\phi+\grad U(\phi)|^2\,\dd t
=\frac12\int_0^T|\dot\phi-\grad U(\phi)|^2\,\dd t
+2(U(\phi(T))-U(\phi(0))).
$$
これは二つの二乗の差 $4\dot\phi\cdot\grad U$ を時間積分した式である。$x_*$ が臨界点で、$y$ からの勾配降下が $x_*$ に近づくとする。$\grad U$ は $x_*$ の近くで Lipschitz と仮定する。すると準ポテンシャルは $2(U(y)-U(x_*))$ である。下界は上の恒等式から出る。上界を作るには、$y$ から降下時間 $R$ 後の点 $z_R$ までの路を逆向きにたどる。この区間は $\dot\phi=\grad U$ なので二乗の第一項が零である。$x_*$ から $z_R$ への最初の区間を時間 $1$ の直線で補う。その作用は Lipschitz 性と $\grad U(x_*)=0$ により $C|z_R-x_*|^2$ 以下で、$R\to\infty$ で消える。従って下界に近づける。

例えば $U(x)=x^2$、$\sigma=1$、$x_*=0$ では $y=3$ までの準ポテンシャルは $18$ である。線形系の $Q_T=(1-e^{-4T})/4$ で最小作用を取ると $18/(1-e^{-4T})$ で、時間について下限すれば同じ $18$ となる。出口の最尤点や平均退出時間の指数をこの準ポテンシャルから得る定理では、固定した有限時間の LDP に加え、繰返し試行、安定領域への再帰、初期点に一様な評価が必要となる。

### 一次元で退出時間の指数を直接求める

偶関数 $U\in C^2([-a,a])$ が $[0,a]$ で狭義増加とする。$0$ から始まる $\dd X=-U'(X)\,\dd t+\sqrt\eps\,\dd W$ の $(-a,a)$ からの退出を $\tau$ とする。平均退出時間の境界値問題は $\eps u''/2-U'u'=-1$、$u(\pm a)=0$ である。対称解は $u'(0)=0$ なので、$(e^{-2U/\eps}u')'=-(2/\eps)e^{-2U/\eps}$ を二回積分して
$$
u(0)=\frac2\eps\int_0^a\int_0^y e^{2(U(y)-U(z))/\eps}\,\dd z\,\dd y
$$
を得る。候補 $u$ は区間内で非負である。時刻 $t$ と退出で止めた伊藤公式から $\mathbb E(\tau\wedge t)\le u(0)$ が従い、先に平均が有限と分かる。その後に終端値の有界収束を使うと $u(0)=\mathbb E\tau$ となる。

指数の最大値は $2(U(a)-U(0))$ である。上界にはこの最大値と積分領域の面積を使う。下界には $(y,z)=(a,0)$ 近くの正面積の長方形を選ぶと、連続性で指数を最大値に任意に近づけられる。前因子 $2/\eps$ は $\eps\log$ の極限で消える。従って
$$
\lim_{\eps\downarrow0}\eps\log\mathbb E_0\tau=2(U(a)-U(0)).
$$
例えば $U(x)=x^2$、$a=3$ では指数は $18$ である。一般の安定領域での退出位置や平均退出時間には、吸引域と境界に関する追加条件が必要である。一般の領域の退出定理の証明は、初期点に一様な出口への到達確率と吸引域へ戻る時間を制御する必要があるため、:cite[dembo-zeitouni] を参照する。

## 演習問題

::::exercise{level=basic}
$\dd X=2\,\dd W$ の生成作用素を求めよ。
:::solution
$a=4$、$b=0$ を :eqref[eq-sa-dynkin] に代入すれば $\mathcal L f=2f''$ である。
:::
::::

::::exercise{level=basic}
Brown 半群を $g(x)=\sin(3x)$ に施せ。
:::solution
$P_t\sin(3x)=\Im(e^{3ix}\mathbb Ee^{3iB_t})=e^{-9t/2}\sin(3x)$ である。
:::
::::

::::exercise{level=basic}
定数 $c=2,f=0,g=3$ の終端値問題を解け。
:::solution
:eqref[eq-sa-feynman-kac] により $u(t,x)=3e^{-2(T-t)}$ である。直接微分でも $u_t=2u$、$u(T,x)=3$ を確認できる。
:::
::::

::::exercise{level=basic}
標準 Brown 運動の $(-2,3)$ からの平均退出時間を $x=1$ で求めよ。
:::solution
$u(x)=(x+2)(3-x)$ は二階微分が $-2$ で両端が零である。止めた公式で $\mathbb E_1\tau=u(1)=6$ である。
:::
::::

::::exercise{level=advanced}
OU の不変分布が $N(0,2)$ となるよう、$\kappa=1$ のとき $\nu>0$ を決めよ。
:::solution
分散の不変条件は $\nu^2/(2\kappa)=2$ であるから $\nu=2$ である。分散 $2$ が時刻 $t$ で $2e^{-2t}$ と $2(1-e^{-2t})$ に分かれ、和が $2$ になる。
:::
::::

::::exercise{level=advanced}
$\mathcal L f=f'+f''$ に対する前向き方程式を、定数ドリフトと拡散の法則で検証せよ。
:::solution
解の密度は初期値 $0$ のとき $p(t,x)=(4\pi t)^{-1/2}e^{-(x-t)^2/(4t)}$ である。$q(t,y)$ を同じ分散の中心正規密度とすれば $q_t=q_{yy}$、$p(t,x)=q(t,x-t)$ である。連鎖律で $p_t=p_{xx}-p_x$ を得る。
:::
::::

::::exercise{level=advanced}
$b=1,\sigma=1,T=2$ で $0$ から $4$ への最小作用を求めよ。
:::solution
$\int_0^2(\dot\phi-1)\,\dd t=4-2=2$。Cauchy–Schwarz で作用は $2^2/(2\cdot2)=1$ 以上で、$\phi(t)=2t$ が等号を与える。
:::
::::

::::exercise{level=advanced}
$P_t$ に対し二つの確率測度の積分が全有界連続関数上で一致するとき、二測度が同じである理由を示せ。
:::solution
閉集合 $F$ の指示関数は $f_n(x)=\max(0,1-nd(x,F))$ の減少極限である。従って有界収束から閉集合上で測度が一致する。閉集合は π 系で Borel σ-加法族を生成するため、π–λ 定理により全体で一致する。
:::
::::

::::exercise{level=advanced}
$\dd X=-X\,\dd t+\sqrt\eps\,\dd W$、$X_0=0$ とする。時刻 $T$ に $X_T\ge2$ となる確率について $\lim_{\eps\downarrow0}\eps\log P(X_T\ge2)$ を求めよ。
:::solution
終端 Gram 行列は $Q_T=(1-e^{-2T})/2$ なので、終端値 $y$ の率関数は $y^2/(1-e^{-2T})$ である。これは $y\ge2$ で増加し、閉集合 $[2,\infty)$ と開集合 $(2,\infty)$ で下限は等しく $4/(1-e^{-2T})$ である。従って極限は $-4/(1-e^{-2T})$ である。開集合からの下限と閉集合からの上限が同じなので、境界点に確率が集まるかどうかについて別の仮定は要らない。
:::
::::

::::exercise{level=advanced}
Brown 運動を $x=1$ から始め、$0$ で殺す。時刻 $t=1$ まで生存する確率を標準正規分布関数 $\Phi$ で表せ。
:::solution
半直線の殺した核を積分する。$\int_0^\infty q_1(y-1)\dd y=\Phi(1)$、$\int_0^\infty q_1(y+1)\dd y=1-\Phi(1)$ なので、生存確率は $2\Phi(1)-1$ である。全質量と流束の関係から、失われる側の確率は原点到達時刻の分布関数になる。
:::
::::

## 確認テスト

::::quiz
後ろ向き方程式で空間微分を行う変数を選べ。

- [x] 初期位置
- [ ] 未来の位置の密度だけ
- [ ] Brown 運動の標本番号

:::explanation
$P_tg(x)=\mathbb E_xg(X_t)$ は終端関数を固定し初期点 $x$ を動かした関数である。前向き方程式では初期法則を固定して現在位置の分布を動かし、形式随伴を用いる。
:::
::::

::::quiz{answer="6" tol="1e-9"}
Brown 運動が $0$ から $(-2,3)$ を出る平均時間を求めよ。

:::explanation
$u(x)=(x+2)(3-x)$ は $\mathcal Lu=-1$、境界値零を満たす。退出表現で $u(0)=6$ が平均時間となる。
:::
::::

::::quiz
率関数の制御 $v$ が零であるときの路を選べ。

- [x] $\dot\phi=b(\phi)$ を満たす路
- [ ] 常に任意の定数路
- [ ] 常に速さ $1$ の直線

:::explanation
制約 $\dot\phi=b(\phi)+\sigma v$ で $v=0$ とすれば決定論的なドリフト方程式だけが残る。作用は被積分関数が零なので零になる。
:::
::::

::::quiz{answer="2" tol="1e-9"}
$b=0,\sigma=1,T=1$ で $0$ から $2$ へ向かう最小作用を求めよ。

:::explanation
$\int_0^1\dot\phi\,\dd t=2$ なので Cauchy–Schwarz から作用は $2^2/2=2$ 以上である。$\phi(t)=2t$ が等号に達する。
:::
::::

::::quiz{answer="0.5" tol="1e-9"}
Brown 運動のレゾルベントで $\lambda=2$ とするとき、その積分核の全質量はいくらか。

:::explanation
熱核は各時刻に全質量 $1$ を持つ。従って割引核の全質量は $\int_0^\infty e^{-2t}\dd t=1/2$ である。明示核 $e^{-2|z|}/2$ を $z$ で積分しても同じ値になる。
:::
::::

::::quiz
密度を持たない初期測度について、前向き方程式を扱う方法を選べ。

- [x] 試験関数で測度を積分した弱形式を使う
- [ ] 初期時刻の測度を常に Lebesgue 密度とみなす
- [ ] 前向き方程式の符号を反転させる

:::explanation
:ref[prop-sa-fokker-planck] の弱形式は $\mu_t$ を測度として試験関数と積分する。例えば $\mu_0=\delta_x$ でも $\int f\,\dd\mu_0=f(x)$ と定まる。後から密度が存在すれば、部分積分で密度に対する分布微分の方程式へ移れる。
:::
::::

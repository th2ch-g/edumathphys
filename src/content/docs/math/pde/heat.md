---
title: 熱方程式
description: 熱核による解の表示と一意性を証明し、無限遠での増大条件と時間を逆にたどる問題を調べる。
notation:
  - symbol: '\kappa'
    meaning: 拡散係数
  - symbol: 'K_t'
    meaning: 全空間の熱核
  - symbol: 'E_r'
    meaning: 基本解の等位集合で定めた熱球
---

## この章の目標

- 熱核の形を求め、正規化と畳み込み性を確かめる。
- 全空間の初期値問題を解き、平均値公式と最大値原理を導く。
- 増大条件を外すと一意性が失われる理由を説明する。
- 固有関数展開とエネルギー法で有界領域の解を調べる。

## 基本解と熱核

温度 $u$ と位置 $x$ を観測すると、Fourier の熱伝導則は熱流束を $-\lambda Du$ とする。密度 $\rho$ と比熱 $c_p$ と熱伝導率 $\lambda$ が一定で熱源がないとすると、内部エネルギーの収支は $\rho c_p u_t=\lambda\Delta u$ である。以下の熱方程式はこの定数係数の拡散を表す。$\kappa=\lambda/(\rho c_p)$ の単位は $\mathrm{m^2\,s^{-1}}$ である。例えば $\kappa=10^{-6}\,\mathrm{m^2\,s^{-1}}$ なら $1\,\mathrm s$ での一方向の標準偏差 $\sqrt{2\kappa t}$ はおよそ $1.4\,\mathrm{mm}$ である。

:::definition[熱核]{#def-pde-heat-kernel}
$\kappa>0$ とし、$u_t-\kappa\Delta u=0$ を :term[熱方程式]{en="heat equation" yomi="ねつほうていしき"} と呼ぶ。$t>0$ で

$$
K_t(x)=(4\pi\kappa t)^{-n/2}\exp(-|x|^2/(4\kappa t))
\label{eq-pde-heat-kernel}
$$

を :term[熱核]{en="heat kernel" yomi="ねつかく"} と呼ぶ。
:::

:::proposition[熱核の規格化と半群性]{#prop-pde-heat-kernel}
$K_t$ は熱方程式を満たし、$\int K_t=1$、$K_t*K_s=K_{t+s}$ である。また時空間で $K_t=0$（$t<0$）と延長すると $(\partial_t-\kappa\Delta)K=\delta_{(0,0)}$ となる。
:::

:::proof
総量を時間によらず一定にする相似形 $t^{-n/2}F(x/\sqrt t)$ を代入すると $\divg(\kappa DF+yF/2)=0$ になる。放射対称で原点まで滑らかなら $r^{n-1}(\kappa F'(r)+rF(r)/2)$ は定数 $0$ である。解は $F=C e^{-|y|^2/(4\kappa)}$ で、一次元 Gauss 積分を各座標に使えば $C=(4\pi\kappa)^{-n/2}$ である。畳み込みでは

$$
\frac{|x-y|^2}{t}+\frac{|y|^2}{s}
=\frac{t+s}{ts}\left|y-\frac{s x}{t+s}\right|^2+\frac{|x|^2}{t+s}
$$

となるため、$y$ を積分すれば $K_{t+s}$ を得る。基本解の性質は $t=\eps$ から部分積分して境界項 $\int K_\eps\phi(\eps,x)\,\dd x$ を取り、$x=\sqrt\eps\,y$ の置換で極限が $\phi(0,0)$ となることから従う。
:::

:::example[広がった熱分布]
一次元で $\kappa=1$ とする。$K_1(x)=e^{-x^2/4}/\sqrt{4\pi}$ で、$K_1(2)/K_1(0)=e^{-1}$ である。$\int x^2K_t(x)\,\dd x=2t$ は $xK_t=-2t K_t'$ を積分して得られる。
:::

:::example[連続する二段階の拡散]
時刻 $0$ に $K_2$ の分布から出発すると、時間 $3$ 後の分布は $K_3*K_2=K_5$ である。一次元では初期の分散 $4\kappa$ に追加の分散 $6\kappa$ が足される。
:::

### 微分した熱核の大きさ

$K_t(x)=t^{-n/2}K_1(x/\sqrt t)$ と書く。$x$ で一回微分すると $t^{-1/2}$ が出る。時間で一回微分すると、前の係数の微分と変数 $x/\sqrt t$ の微分の両方に $t^{-1}$ が出る。帰納的に任意の $m\ge0$、多重指数 $\alpha$ について
$$
\partial_t^m D^\alpha K_t(x)
=t^{-n/2-m-|\alpha|/2}P_{m,\alpha}(x/\sqrt t)
e^{-|x|^2/(4\kappa t)}
$$
となる。$P_{m,\alpha}$ は多項式で、時間微分では $y\cdot D_y$ が加わり、空間微分では $D_y$ が加わる。多項式に $e^{-|y|^2/(8\kappa)}$ を掛けると有界なので、絶対値は $C t^{-n/2-m-|\alpha|/2}e^{-|x|^2/(8\kappa t)}$ 以下である。変数変換により $L^1$ ノルムは $C t^{-m-|\alpha|/2}$ 以下になる。これを畳み込みに使えば有界な初期値から $t>0$ で任意階の微分が作れる。

:::example[偶関数の最大勾配が現れる場所]
一次元で $|\partial_xK_t|=|x|K_t/(2\kappa t)$ である。$x>0$ における $x e^{-x^2/(4\kappa t)}$ の微分は $e^{-x^2/(4\kappa t)}(1-x^2/(2\kappa t))$ なので最大は $x=\sqrt{2\kappa t}$ である。$\kappa=2,\ t=1$ なら $x=2$ で最大値 $e^{-1/2}/(2\sqrt{8\pi})$ である。中心の値は高いが、対称性により最も急な場所は中心から外れる。
:::

:::example[Gaussian の質量と幅を別々に追う]
$n=2$、$\kappa=1$、$g(x)=e^{-|x|^2/4}$ の熱の広がりを調べる。中心の高さは変わるが全量は変わらない。

中心からの距離を $r$、空間次元を $n$ とする。Gauss 型の解を保つには $K=C(t)e^{-r^2/(4d(t))}$ と置く。熱方程式に代入すると $r^2$ の係数から $d^{\prime}=\kappa$、定数項から $C^{\prime}/C=-n\kappa/(2d)$ が得られる。初期データの幅と高さに合わせると、次の解になる。

$$
\begin{aligned}
d(t)&=t + 1\\
u(t,r)&=\frac{e^{- \frac{r^{2}}{4 t + 4}}}{t + 1}\\
u_t&=\frac{\partial}{\partial t} \frac{e^{- \frac{r^{2}}{4 t + 4}}}{t + 1}=\frac{\left(r^{2} - 4 t - 4\right) e^{- \frac{r^{2}}{4 t + 4}}}{4 \left(t + 1\right)^{3}}\\
u_r&=\frac{\partial}{\partial r} \frac{e^{- \frac{r^{2}}{4 t + 4}}}{t + 1}=- \frac{r e^{- \frac{r^{2}}{4 t + 4}}}{2 \left(t + 1\right)^{2}}\\
u_{rr}&=\frac{\partial^{2}}{\partial r^{2}} \frac{e^{- \frac{r^{2}}{4 t + 4}}}{t + 1}=\frac{\left(r^{2} - 2 t - 2\right) e^{- \frac{r^{2}}{4 t + 4}}}{4 \left(t + 1\right)^{3}}\\
\Delta u&=\frac{\left(r^{2} - 4 t - 4\right) e^{- \frac{r^{2}}{4 t + 4}}}{4 \left(t + 1\right)^{3}}
\end{aligned}
$$

放射状の Laplacian は $u_{rr}+(n-1)u_r/r$ である。$r=0$ の値は $u_r/r$ の極限で読む。次に、中心の値と半値となる半径を比べ、指定された時刻と位置に代入する。

$$
\begin{aligned}
u_t-\kappa\Delta u&=\frac{4 r^{2} e^{- \frac{r^{2}}{4 t + 4}}}{\left(t + 1\right) \left(4 t + 4\right)^{2}} - \frac{e^{- \frac{r^{2}}{4 t + 4}}}{\left(t + 1\right)^{2}} - \frac{\left(r^{2} - 4 t - 4\right) e^{- \frac{r^{2}}{4 t + 4}}}{4 \left(t + 1\right)^{3}}=0\\
u(t,0)&=\frac{1}{t + 1},\qquad r=0\\
r_{1/2}(t)&=2 \sqrt{t + 1} \sqrt{\log{\left(2 \right)}}\\
u(t_*,r_*)&=\frac{1}{2 e^{\frac{1}{2}}},\qquad t=1,r=2\\
u_r(t_*,r_*)&=- \frac{1}{4 e^{\frac{1}{2}}},\qquad t=1,r=2
\end{aligned}
$$

全量を求めるには直交座標で積に戻す。一方向の因子を $g$ とし、その中心を $x_0$ とする。他の方向でも各中心座標に関して同じ Gauss 積分が得られる。積分が分離されるので、全量は一方向の積分の $n$ 乗、正規化した一方向の分散は二次の積分を零次の積分で割った値になる。

$$
\begin{aligned}
g(t,x)&=\frac{e^{- \frac{x^{2}}{4 t + 4}}}{\sqrt{t + 1}}\\
M_1&=\int\limits_{-\infty}^{\infty} \frac{e^{- \frac{x^{2}}{4 t + 4}}}{\sqrt{t + 1}}\, dx=\frac{\sqrt{\pi} \sqrt{4 t + 4}}{\sqrt{t + 1}}\\
M_2&=\int\limits_{-\infty}^{\infty} \frac{x^{2} e^{- \frac{x^{2}}{4 t + 4}}}{\sqrt{t + 1}}\, dx=\frac{\sqrt{\pi} \left(4 t + 4\right)^{\frac{3}{2}}}{2 \sqrt{t + 1}}\\
\text{全量}&=\frac{\pi \left(4 t + 4\right)}{t + 1}=4 \pi\\
\text{一方向の分散}&=2 t + 2=2 \left(t + 1\right)
\end{aligned}
$$

時刻 $1$ で分散は一方向当たり $4$、中心値は $1/2$ である。熱核と違い初期幅が正なので、時刻零でも尖っていない。初期幅 $A$ を零にする極限で熱核を得るには、同時に初期の全量を $1$ に保つ必要がある。
:::

### 分散を和として求める別の計算

一次元で平均 $a,b$、分散 $2\kappa t,2\kappa s$ の二つの正規化された Gauss 密度を畳み込むと、平均は $a+b$ となる。実際 $x=y+z$ として二重積分の一次項を分ける。中心からの差の二乗は $(y-a)^2+2(y-a)(z-b)+(z-b)^2$ であり、交差項の積分は両側の中心化した一階モーメントが零なので消える。従って分散は $2\kappa(t+s)$ である。例えば $\kappa=1,t=1,s=2,a=2,b=-1$ なら平均は $1$、分散は $6$ である。平方完成で得た半群性は、拡散の長さを二乗して足すことと整合する。

## 初期値問題の解の表示

:::theorem[全空間の初期値問題]{#thm-pde-heat-cauchy}
$g$ は有界な一様連続関数とする。$u(t,x)=K_t*g(x)$ は $t>0$ で滑らかで、$u_t=\kappa\Delta u$ を満たし、$t\downarrow0$ で $g$ に一様収束する。右辺 $f$ が時空間で滑らかで空間の導関数も有界なら、非斉次方程式には

$$
u(t)=K_t*g+\int_0^t K_{t-s}*f(s)\,\dd s
\label{eq-pde-heat-duhamel}
$$

を用いる。
:::

:::proof
正の時間から離れた区間では熱核の導関数は多項式と Gaussian の積であり、局所的に共通の可積分な上界を持つ。従って積分の下で方程式が成り立つ。$u-g=\int K_t(y)(g(x-y)-g(x))\,\dd y$ であり、小さい $|y|$ では一様連続性、大きい $|y|$ では Gaussian の裾が $0$ に収束することを使う。非斉次項を $t$ で微分すると、端点で $K_0*f(t)=f(t)$ が現れる。内部積分で $\partial_t K$ を $\kappa\Delta K$ と書き、空間の二階微分を $f$ に移すと $s\uparrow t$ でも有界に積分できる。これが $\kappa\Delta$ を非斉次項に作用させたものなので結論を得る。
:::

初期値が $L^p$（$1\le p<\infty$）なら、$\norm{K_t*g-g}_p\le\int K_t(y)\norm{g(\cdot-y)-g}_p\,\dd y$ と $L^p$ 平行移動の連続性で強収束する。滑らかな台付き関数の密度から平行移動の連続性が従う。核の非負性と質量 $1$ により $\norm{K_t*g}_p\le\norm g_p$ である。$g\in L^1$ なら Fubini で積分値も保存される。

:::example[周期的な初期温度]
$g(x)=4+\cos(2x)$、$\kappa=3$ とする。$\cos(2x)$ は $\partial_{xx}$ の固有値 $-4$ を持つので $u=4+e^{-12t}\cos(2x)$ が解である。核表示でこの係数は $\int K_t(y)\cos(2y)\,\dd y=e^{-12t}$ に対応する。
:::

:::example[温度の段差]
$g=\mathbf1_{(0,\infty)}$ のとき $u(t,x)=\frac1{\sqrt\pi}\int_{-\infty}^{x/\sqrt{4\kappa t}}e^{-z^2}\,\dd z$ であり、$x=0$ では $1/2$ である。任意の正の時刻で全空間に正の値を持つため、熱方程式は波動方程式と異なる伝播を示す。
:::

### 初期値が測度の場合

有限符号付き測度 $\mu$ に対し $u(t,x)=\int K_t(x-y)\dd\mu(y)$ と置く。前節の熱核の微分評価は $t$ に正の下限があれば $y$ に一様なので、積分の下で微分して熱方程式を得る。$\phi\in C_c^\infty$ に対して $\int u(t,x)\phi(x)\dd x=\int (K_t*\phi)(y)\dd\mu(y)$ であり、右辺は $K_t*\phi\to\phi$ の一様収束から $\int\phi\dd\mu$ へ収束する。従って初期値は測度として実現される。符号付きの $\mu$ に対しても Fubini により $\int u(t,x)\dd x=\mu(\R^n)$ であり、符号があっても $\norm{u(t)}_1\le|\mu|(\R^n)$ である。非負で零でない測度なら任意の正時刻で解は全空間で正になる。

:::example[二点の正負の熱源]
一次元で $\mu=2\delta_1-\delta_{-1}$ とすれば $u=2K_t(x-1)-K_t(x+1)$ で全質量は $1$ である。符号を変える点は $2\exp(x/(\kappa t))=1$、すなわち $x=-\kappa t\log2$ である。$x$ がこれより大きければ正、より小さければ負である。中心 $x=0$ では二点からの核の大きさが同じなので $u=K_t(1)>0$。総量の符号だけから点ごとの符号は決まらない。
:::

:::example[位置のずれた初期分布を点ごとに読む]
一次元で $\kappa=2$、$g(x)=e^{-(x-3)^2/8}$ とする。中心からの距離 $r=|x-3|$ で値を記録する。

中心からの距離を $r$、空間次元を $n$ とする。Gauss 型の解を保つには $K=C(t)e^{-r^2/(4d(t))}$ と置く。熱方程式に代入すると $r^2$ の係数から $d^{\prime}=\kappa$、定数項から $C^{\prime}/C=-n\kappa/(2d)$ が得られる。初期データの幅と高さに合わせると、次の解になる。

$$
\begin{aligned}
d(t)&=2 t + 2\\
u(t,r)&=\frac{\sqrt{2} e^{- \frac{r^{2}}{8 t + 8}}}{\sqrt{2 t + 2}}\\
u_t&=\frac{\partial}{\partial t} \frac{\sqrt{2} e^{- \frac{r^{2}}{8 t + 8}}}{\sqrt{2 t + 2}}=\frac{\left(r^{2} - 4 t - 4\right) e^{- \frac{r^{2}}{8 t + 8}}}{8 \left(t + 1\right)^{\frac{5}{2}}}\\
u_r&=\frac{\partial}{\partial r} \frac{\sqrt{2} e^{- \frac{r^{2}}{8 t + 8}}}{\sqrt{2 t + 2}}=- \frac{r e^{- \frac{r^{2}}{8 t + 8}}}{4 \left(t + 1\right)^{\frac{3}{2}}}\\
u_{rr}&=\frac{\partial^{2}}{\partial r^{2}} \frac{\sqrt{2} e^{- \frac{r^{2}}{8 t + 8}}}{\sqrt{2 t + 2}}=\frac{\left(r^{2} - 4 t - 4\right) e^{- \frac{r^{2}}{8 t + 8}}}{16 \left(t + 1\right)^{\frac{5}{2}}}\\
\Delta u&=\frac{\left(r^{2} - 4 t - 4\right) e^{- \frac{r^{2}}{8 t + 8}}}{16 \left(t + 1\right)^{\frac{5}{2}}}
\end{aligned}
$$

放射状の Laplacian は $u_{rr}+(n-1)u_r/r$ である。$r=0$ の値は $u_r/r$ の極限で読む。次に、中心の値と半値となる半径を比べ、指定された時刻と位置に代入する。

$$
\begin{aligned}
u_t-\kappa\Delta u&=\frac{8 \sqrt{2} r^{2} e^{- \frac{r^{2}}{8 t + 8}}}{\sqrt{2 t + 2} \left(8 t + 8\right)^{2}} - \frac{\sqrt{2} e^{- \frac{r^{2}}{8 t + 8}}}{\left(2 t + 2\right)^{\frac{3}{2}}} - \frac{\left(r^{2} - 4 t - 4\right) e^{- \frac{r^{2}}{8 t + 8}}}{8 \left(t + 1\right)^{\frac{5}{2}}}=0\\
u(t,0)&=\frac{\sqrt{2}}{\sqrt{2 t + 2}}=\frac{1}{\sqrt{t + 1}},\qquad r=0\\
r_{1/2}(t)&=2 \sqrt{2 t + 2} \sqrt{\log{\left(2 \right)}}\\
u(t_*,r_*)&=\frac{\sqrt{2}}{2 e^{\frac{1}{4}}},\qquad t=1,r=2\\
u_r(t_*,r_*)&=- \frac{\sqrt{2}}{8 e^{\frac{1}{4}}},\qquad t=1,r=2
\end{aligned}
$$

全量を求めるには直交座標で積に戻す。一方向の因子を $g$ とし、その中心を $x_0$ とする。他の方向でも各中心座標に関して同じ Gauss 積分が得られる。積分が分離されるので、全量は一方向の積分の $n$ 乗、正規化した一方向の分散は二次の積分を零次の積分で割った値になる。

$$
\begin{aligned}
g(t,x)&=\frac{\sqrt{2} e^{- \frac{\left(x - 3\right)^{2}}{8 t + 8}}}{\sqrt{2 t + 2}}\\
M_1&=\int\limits_{-\infty}^{\infty} \frac{\sqrt{2} e^{- \frac{\left(x - 3\right)^{2}}{8 t + 8}}}{\sqrt{2 t + 2}}\, dx=\frac{\sqrt{2} \sqrt{\pi} \sqrt{8 t + 8}}{\sqrt{2 t + 2}}\\
M_2&=\int\limits_{-\infty}^{\infty} \frac{\sqrt{2} \left(x - 3\right)^{2} e^{- \frac{\left(x - 3\right)^{2}}{8 t + 8}}}{\sqrt{2 t + 2}}\, dx=\frac{\sqrt{2} \sqrt{\pi} \left(8 t + 8\right)^{\frac{3}{2}}}{2 \sqrt{2 t + 2}}\\
\text{全量}&=\frac{\sqrt{2} \sqrt{\pi} \sqrt{8 t + 8}}{\sqrt{2 t + 2}}=2 \sqrt{2} \sqrt{\pi}\\
\text{一方向の分散}&=4 t + 4=4 \left(t + 1\right)
\end{aligned}
$$

中心は $3$ のままであり、同じ $r$ を持つ $x=1,5$ で温度は等しい。両点の $x$ 微分は符号が反対である。すなわち、熱流が中心から両側へ出ることと、全体の質量の中心が動かないことは同時に成り立つ。
:::

### 差の平滑化と初期値の復元

$g\in L^p(\R^n)$、$1\le p\le q\le\infty$ とする。畳み込みの Young の不等式で $1+1/q=1/p+1/r$ の $r$ を選べば
$$
\norm{D^\alpha(K_t*g)}_q
\le \norm{D^\alpha K_t}_r\norm g_p
\le C t^{-|\alpha|/2-n(1/p-1/q)/2}\norm g_p .
$$
最後のべきは $x=\sqrt t y$ によって $L^r$ ノルムに現れる $t^{n/(2r)}$ から出る。$p<\infty$ での初期値の強収束を示すには、差を $\int K_t(y)(g(\cdot-y)-g)\dd y$ と書く。$L^p$ の平行移動は連続であり、核の質量は $|y|<\delta$ に集中するので、差は零へ行く。$p=\infty$ では平行移動の強連続性が一般にはない。有界一様連続な $g$ なら一様ノルムで収束する。

:::example[段差で一様な初期収束が壊れる]
一次元で $g=\mathbf1_{(0,\infty)}$ とし、$u(t,x)=K_t*g$ を考える。対称性から $u(t,0)=1/2$ である。$x>0$ では $u\in(1/2,1)$、$x<0$ では $u\in(0,1/2)$ であり、$x$ を零へ近づけるとそれぞれの誤差は $1/2$ に近づく。従って任意の $t>0$ で $\norm{u(t)-g}_\infty=1/2$ である。一方固定した $x\ne0$ で $t\downarrow0$ とすれば $g(x)$ に収束する。$x=0$ の値だけを変更しても本質上限や畳み込みは変化しない。
:::

## 平均値公式

時刻 $t_0$ の値を過去の値から表す。空間の球だけを使うと未知関数が満たす時間微分が反映されない。そこで熱核の等位面を境界に選ぶ。

:::theorem[熱球の平均値公式]{#thm-pde-heat-mean}
$\Phi(y,s)=K_s(y)$、$E_r=\{(y,s):s>0,\ \Phi(y,s)>r^{-n}\}$ とする。$(x_0-y,t_0-s)$ で写した $\overline E_r$ の近傍に古典解 $u$ があれば

$$
u(t_0,x_0)=\frac{1}{r^n}\int_{E_r}u(t_0-s,x_0-y)\frac{|y|^2}{4\kappa s^2}\,\dd y\,\dd s .
\label{eq-pde-heat-mean}
$$

:::

:::proof
$w(y,s)=u(t_0-s,x_0-y)$ と置くと $w_s=-\kappa\Delta w$ である。$q\ge r^{-n}$ に対して

$$
\partial_s(w(\Phi-q))=\kappa\divg\big(wD_y\Phi-(\Phi-q)D_yw\big)
$$

である。従って時空間のベクトル場 $(-\kappa(wD_y\Phi-(\Phi-q)D_yw),\,w(\Phi-q))$ の発散は零である。$\{\Phi>q,s>\eps\}$ 上で発散定理を使う。外側の等位面では $\Phi-q=0$ で、外法線は $-D_{y,s}\Phi/|D_{y,s}\Phi|$ である。従ってその面を通る流束は $\kappa w|D_y\Phi|^2/|D_{y,s}\Phi|$ である。底面流束の負号を移すと $\int_{\{\Phi(y,\eps)>q\}}w(y,\eps)(\Phi(y,\eps)-q)\dd y$ と等しい。この底面は $y=\sqrt\eps z$ によって半径が無限大に向かい、元の座標での体積は $0$ に向かう。従って近似恒等式により極限は $w(0,0)$ である。等位面は $D_{y,s}\Phi\ne0$ であり、$q^{-2}\dd q$ を掛けて $r^{-n}$ から $\infty$ まで積分し、余面積公式で体積積分に変える。左辺は $r^n w(0,0)$ であり、右辺の係数は $\kappa|D_y\Phi|^2/\Phi^2=|y|^2/(4\kappa s^2)$ になる。定数 $w=1$ では非負積分の交換で規格化が従う。これで重みの積分が有限だと分かり、一般の $w$ に対してもその絶対値の上界を使って積分を交換できる。
:::

:::example[平均する時刻の範囲]
$K_s(0)=(4\pi\kappa s)^{-n/2}$ が $r^{-n}$ を超える条件は $s<r^2/(4\pi\kappa)$ である。従って $r=2$、$\kappa=1$ のとき平均には中心時刻より $1/\pi$ より古い時刻は現れない。
:::

:::example[熱多項式の平均]
一次元では $u(t,x)=x^2+2\kappa t$ が熱方程式を満たす。$(0,t_0)$ の熱球では $w=y^2+2\kappa(t_0-s)$ であるから、:eqref[eq-pde-heat-mean] と定数解の規格化を合わせると、同じ重みでの $y^2$ の平均は $2\kappa$ と $s$ の平均の積に等しい。
:::

### 熱球の逆向きの判定

$E_r$ は $(y,s)=(rz,r^2\tau)$ で $E_1$ を移した集合である。重みに体積要素を合わせた $r^{-n}|y|^2\dd y\dd s/(4\kappa s^2)$ はこの変換で不変なので、中心 $(t_0,x_0)$ での Taylor 展開を $E_1$ で積分できる。反転と回転の対称性により一次空間項と Hessian の非対角項は零である。$A=\int_{E_1}\tau\,\dd m$、$B=\int_{E_1}z_1^2\dd m$ と書く。ここで $\dd m=|z|^2\dd z\dd\tau/(4\kappa\tau^2)$ は全質量 $1$ である。熱方程式を満たす多項式 $u(t,x)=x_1^2+2\kappa t$ に平均値公式を用いると $B=2\kappa A$ である。従って一般の $C^{2,1}$ 関数では平均と中心値の差が
$$
r^2 A(-u_t+\kappa\Delta u)(t_0,x_0)+o(r^2)
$$
になる。余項が一様に $o(r^2)$ なのは $E_1$ が有界で、空間と時間の導関数が連続だからである。$A>0$ なので全ての十分小さい熱球で平均が値に等しければ PDE が従う。

上解と劣解の平均不等式も同じ流束から得られる。$u_t-\kappa\Delta u\ge0$ なら、反転した $w$ では $w_s+\kappa\Delta w\le0$ である。発散定理に使ったベクトル場の発散は $(\Phi-q)(w_s+\kappa\Delta w)\le0$ となる。従って等位面を通る流束は、底面を通る流束の負号を取った値以下になる。等位面を積分すれば熱球平均は中心値以下である。$-u$ に適用すると、劣解では平均が中心値以上になる。この向きと Taylor 展開を合わせると、滑らかな上解と劣解を小さい熱球の平均不等式から判定できる。

:::example[候補の二次係数を平均だけで決める]
$u=x_1^2+ct$ を考える。中心 $(0,0)$ での値は $0$、熱球平均は $r^2(2\kappa-c)A$ である。$c=2\kappa$ でだけ平均と値は一致する。例えば $\kappa=3$ で $c=5$ なら平均の方が大きく、$u_t-\kappa\Delta u=-1$ の符号と一致する。
:::

:::example[熱球で平均される第四次の解]
一次元で $\kappa=1$、$g=x^4+x^2$、$P=\partial_{xx}$ とする。小さい熱球で二次モーメントだけでは足りない項を調べる。

時間一階の方程式なので $g,Pg,P^2g,\ldots$ を順に作る。空間の多項式に二階の作用素を掛けるとこの場合は次数が少なくとも二つ下がる。有限回の後に零になることを確かめると、級数は多項式に打ち切れる。

$$
\begin{aligned}
P^{0}g&=x^{4} + x^{2}\\
P^{1}g&=2 \left(6 x^{2} + 1\right)\\
P^{2}g&=24\\
P^{3}g&=0
\end{aligned}
$$

$P^j g$ に $t^j/j!$ を掛けて和を作る。時間微分で階乗が一つ減ることにより、空間作用素を一回作用させた列と一致する。さらに時刻を零に戻すと初期データが得られる。

$$
\begin{aligned}
u(t,x)&=12 t^{2} + 2 t \left(6 x^{2} + 1\right) + x^{4} + x^{2}\\
u_t&=\frac{\partial}{\partial t} \left(12 t^{2} + 2 t \left(6 x^{2} + 1\right) + x^{4} + x^{2}\right)=2 \left(12 t + 6 x^{2} + 1\right)\\
Pu&=2 \left(12 t + 6 x^{2} + 1\right)\\
u(0,x)&=x^{4} + x^{2},\qquad t=0
\end{aligned}
$$

空間の導関数を測ると、時間の多項式の各係数が局所でどう効くかが分かる。初期関数との差は単に平行移動で置き換えず、以下では導関数と時刻を代入して計算する。

$$
\begin{aligned}
\partial_{x} u&=\frac{\partial}{\partial x} \left(12 t^{2} + 2 t \left(6 x^{2} + 1\right) + x^{4} + x^{2}\right)=2 x \left(12 t + 2 x^{2} + 1\right)\\
g&=0,\qquad x=0\\
u&=14,\qquad t=1,x=0\\
u_t&=26,\qquad t=1,x=0
\end{aligned}
$$

中心を $(0,t_0)$ とし、熱球の座標で $u(t_0-s,y)$ を展開する。$y^2-2s$ の平均は零であるため、係数が $t_0$ に比例する部分は相殺する。さらに $y^4-12sy^2+12s^2$ の平均も零になる。これは $u(t,x)=x^4+12tx^2+12t^2$ に平均値公式を適用した帰結であり、熱球の高次モーメントは時間と空間に独立には分かれない。
:::

### 一次元熱球の時間モーメントを積分する

$n=1$ で $S=r^2/(4\pi\kappa)$ とおく。$K_s(y)>r^{-1}$ の条件は $0<s<S$ かつ $|y|<a(s)$、$a(s)=(2\kappa s\log(S/s))^{1/2}$ である。時刻の平均分布を求めるため、重み $r^{-1}y^2/(4\kappa s^2)$ を $-a(s)<y<a(s)$ で積分する。得られる密度は
$$
m(s)=\frac{a(s)^3}{6r\kappa s^2}
=\frac{(2\kappa)^{3/2}}{6r\kappa}s^{-1/2}\log(S/s)^{3/2}.
$$
$t=\log(S/s)$ と置けば $s=Se^{-t}$、$\dd s=-Se^{-t}\dd t$ である。従って $\int_0^S s^j m(s)\dd s$ は定数倍の $S^{j+1/2}\int_0^\infty e^{-(j+1/2)t}t^{3/2}\dd t$ になる。$j=0$ の積分で割り、前節の規格化 $\int m=1$ を用いると
$$
\int_0^S s^j m(s)\dd s=S^j(2j+1)^{-5/2}\qquad(j\ge0).
$$
例えば $r=2\sqrt\pi$、$\kappa=1$ なら $S=1$ で、平均時間差は $3^{-5/2}$ である。時刻の範囲は $[0,1]$ だが、重みは古い時刻に一様ではない。この計算で得た $j=1$ の平均に $2\kappa$ を掛けると位置の二次モーメントになることは、熱多項式 $x^2+2\kappa t$ の平均から確かめられる。

## 最大値原理と一意性

:::theorem[放物境界での最大値原理]{#thm-pde-heat-maximum}
$\Omega$ は有界で、$u$ は $\overline\Omega\times[0,T]$ で連続、内部で $C^{2,1}$ とする。$u_t-\kappa\Delta u\le0$ なら $\sup u$ は初期面または側面で達成される。
:::

:::proof
$u_\eps=u-\eps t$ は厳密な不等式を満たす。まず $0<S<T$ に時刻を切る。$[0,S]$ 上で初期面と側面を除いた最大点では空間 Hessian は非正、時間微分は非負である。時刻 $S$ でも左からの時間微分は非負なので、厳密な不等式と矛盾する。従って $u_\eps$ は初期面か側面で最大である。$\eps\downarrow0$ とした後で $S\uparrow T$ とすれば、連続性により終時刻まで結論が成り立つ。
:::

下からは $-u$ に適用できるため、同じ初期値と側面値の解は一意である。全空間で有界な解に対しては $b=1+|x|^2+(2n\kappa+1)t$ を使うと $b_t-\kappa\Delta b=1$ である。差 $w$ の初期値が $0$ なら、任意の $\eps>0$ について十分大きい球の側面で $w-\eps b\le0$ となるので全空間で $w\le\eps b$ である。$-\!w$ にも適用して一意性が従う。

指数増大 $|w|\le C e^{A|x|^2}$ まで許す場合には、短時間で $b=(1-4\kappa Bt)^{-n/2}\exp(B|x|^2/(1-4\kappa Bt))$ を用いる。$B>A$ とし $t<1/(4\kappa B)$ に限定すると、この $b$ は熱方程式を満たし、球の側面で $w$ より速く増大する。前の比較を繰り返せば短区間で一意であり、同じ幅の区間を順に適用すれば有限時間まで一意である。

:::example[外力に対する差の上界]
同じ初期値と側面値の二解で、外力の差が絶対値 $3$ 以下とする。差を $w$ とすると $w-3t$ と $-w-3t$ はともに劣解で、放物境界上で非正である。従って $|w|\le3t$ である。
:::

:::example[終時刻に最大値を置くこととの違い]
$u(t,x)=e^{-t}\sin x$ を $0<x<\pi$ で考える。正の最大値 $1$ は $t=0$ にあり、終時刻で観測される空間最大値は $e^{-T}$ である。側面と終時刻だけでは時間全体の最大値を押さえられない。
:::

### 強最大値原理と過去への伝播

:::proposition[熱方程式の強最大値原理]{#prop-pde-heat-strong-maximum}
$\Omega$ は連結開集合で、$u$ は $(0,T)\times\Omega$ で $C^{2,1}$ とし、$u_t-\kappa\Delta u\le0$ とする。$u\le M$ が全域で成立し、$0<t_0<T$ の点 $(t_0,x_0)$ で $u=M$ なら、任意の $0<t\le t_0$ と $x\in\Omega$ で $u(t,x)=M$ である。
:::

:::proof
上解 $M-u$ に熱球の平均不等式を適用する。中心値が零であるため、非負な積分値も零である。空間中心の線を除いて重みは正なので、連続性により熱球の内部全体で $u=M$ である。任意の $x\in\Omega$ を $x_0$ と折れ線で結ぶ。折れ線はコンパクトなので領域の外部まで正の距離がある。時間差 $t_0-t>0$ を $N$ 等分し、一段の時間差を $h=(t_0-t)/N$ とする。定数 $C_0>\sqrt{4\pi\kappa}$ を固定して熱球の尺度を $r=C_0\sqrt h$ とすれば、時間差 $h$ の断面は半径 $C_1\sqrt h$ の球になる。実際 $K_h(y)>r^{-n}$ を解くと $C_1=\sqrt{4\kappa n\log(C_0/\sqrt{4\pi\kappa})}>0$ が得られる。熱球全体も尺度変換で空間幅 $O(\sqrt h)$、時間幅 $O(h)$ である。折れ線を同じ長さの $N$ 区間に分けると一段の空間差は $O(1/N)$ であり、十分大きな $N$ では $C_1\sqrt h$ より小さい。さらに各熱球は領域内かつ時刻 $0$ より後に入る。従って熱球の等号を $N$ 回使い、$(t,x)$ に到達できる。$t=t_0$ での等号は連続性から従う。
:::

例えば初期値 $g\ge0$ が恒等的に $0$ でない場合、全空間の有界熱解は任意の正の時刻に各点で正である。もし途中で $u=0$ なら $-u$ が最大値 $0$ を取るので、初期時刻に近い全ての時刻でも $u=0$ となり、初期値への収束と矛盾する。

### 全空間での時間をまたぐ Harnack 評価

$g\ge0$ とし、熱核との畳み込みが有限とする。$0<s<t$ では
$$
u(s,x)\le
\left(\frac{t}{s}\right)^{n/2}
\exp\left(\frac{|x-y|^2}{4\kappa(t-s)}\right)u(t,y).
$$
実際、同じ初期点 $z$ に対する二つの熱核の比を取り、指数を平方完成すると
$$
\frac{|y-z|^2}{t}-\frac{|x-z|^2}{s}
=\frac{|x-y|^2}{t-s}
-\frac{t-s}{st}\left|z-\frac{tx-sy}{t-s}\right|^2 .
$$
第2項が非正なので、$K_s(x-z)$ は $K_t(y-z)$ の上記定数倍以下である。これに $g(z)$ を掛けて積分する。例えば $x=y$ では $u(s,x)\le(t/s)^{n/2}u(t,x)$ である。この向きは、早い時刻の正の値が後の時刻にも下界を与えることに対応する。

### 吸収項と有界な外力

$u_t-\kappa\Delta u+c(t,x)u=f$、$c\ge0$ とする。放物境界で $|u|\le M$、$\sup_x|f(t,x)|\le F(t)$ で、$F$ は連続と仮定する。$W(t)=M+\int_0^t F(s)\dd s$ は $W_t+cW\ge|f|$ を満たすので、$u-W$ と $-u-W$ は劣解で放物境界で非正である。最大値原理の証明で内部の正の最大点を考えると $c(u-W)\ge0$ が加わるだけなので、$|u|\le W$ を得る。$c\ge c_0>0$ が分かるなら、$W'+c_0W=F$ として $W=e^{-c_0t}M+\int_0^t e^{-c_0(t-s)}F(s)\dd s$ に改善できる。ただしこの改善には時間ごとの側面境界値も $W(t)$ 以下という条件が要る。

:::example[定数吸収がある強制の上界]
$u_t-\Delta u+2u=6$、初期値 $0$、斉次 Dirichlet 境界を考える。下解 $0$ と上解 $W=3(1-e^{-2t})$ で挟める。例えば $t=(\log2)/2$ で $0\le u\le3/2$ である。外力の積分だけの上界はこの時刻で $3\log 2$ であり、吸収を使った方が小さい。領域に境界がなく初期値が空間一定なら $W$ 自体が解になり、上界に達する。
:::

:::example[最大値は下がっても局所温度が上がる]
一次元で $\kappa=1$、$g=e^{-x^2/4}$ を取り、観測点 $x=4$ と中心を比べる。

中心からの距離を $r$、空間次元を $n$ とする。Gauss 型の解を保つには $K=C(t)e^{-r^2/(4d(t))}$ と置く。熱方程式に代入すると $r^2$ の係数から $d^{\prime}=\kappa$、定数項から $C^{\prime}/C=-n\kappa/(2d)$ が得られる。初期データの幅と高さに合わせると、次の解になる。

$$
\begin{aligned}
d(t)&=t + 1\\
u(t,r)&=\frac{e^{- \frac{r^{2}}{4 t + 4}}}{\sqrt{t + 1}}\\
u_t&=\frac{\partial}{\partial t} \frac{e^{- \frac{r^{2}}{4 t + 4}}}{\sqrt{t + 1}}=\frac{\left(r^{2} - 2 t - 2\right) e^{- \frac{r^{2}}{4 t + 4}}}{4 \left(t + 1\right)^{\frac{5}{2}}}\\
u_r&=\frac{\partial}{\partial r} \frac{e^{- \frac{r^{2}}{4 t + 4}}}{\sqrt{t + 1}}=- \frac{r e^{- \frac{r^{2}}{4 t + 4}}}{2 \left(t + 1\right)^{\frac{3}{2}}}\\
u_{rr}&=\frac{\partial^{2}}{\partial r^{2}} \frac{e^{- \frac{r^{2}}{4 t + 4}}}{\sqrt{t + 1}}=\frac{\left(r^{2} - 2 t - 2\right) e^{- \frac{r^{2}}{4 t + 4}}}{4 \left(t + 1\right)^{\frac{5}{2}}}\\
\Delta u&=\frac{\left(r^{2} - 2 t - 2\right) e^{- \frac{r^{2}}{4 t + 4}}}{4 \left(t + 1\right)^{\frac{5}{2}}}
\end{aligned}
$$

放射状の Laplacian は $u_{rr}+(n-1)u_r/r$ である。$r=0$ の値は $u_r/r$ の極限で読む。次に、中心の値と半値となる半径を比べ、指定された時刻と位置に代入する。

$$
\begin{aligned}
u_t-\kappa\Delta u&=\frac{4 r^{2} e^{- \frac{r^{2}}{4 t + 4}}}{\sqrt{t + 1} \left(4 t + 4\right)^{2}} - \frac{e^{- \frac{r^{2}}{4 t + 4}}}{2 \left(t + 1\right)^{\frac{3}{2}}} - \frac{\left(r^{2} - 2 t - 2\right) e^{- \frac{r^{2}}{4 t + 4}}}{4 \left(t + 1\right)^{\frac{5}{2}}}=0\\
u(t,0)&=\frac{1}{\sqrt{t + 1}},\qquad r=0\\
r_{1/2}(t)&=2 \sqrt{t + 1} \sqrt{\log{\left(2 \right)}}\\
u(t_*,r_*)&=\frac{\sqrt{2}}{2 e^{2}},\qquad t=1,r=4\\
u_r(t_*,r_*)&=- \frac{\sqrt{2}}{2 e^{2}},\qquad t=1,r=4
\end{aligned}
$$

全量を求めるには直交座標で積に戻す。一方向の因子を $g$ とし、その中心を $x_0$ とする。他の方向でも各中心座標に関して同じ Gauss 積分が得られる。積分が分離されるので、全量は一方向の積分の $n$ 乗、正規化した一方向の分散は二次の積分を零次の積分で割った値になる。

$$
\begin{aligned}
g(t,x)&=\frac{e^{- \frac{x^{2}}{4 t + 4}}}{\sqrt{t + 1}}\\
M_1&=\int\limits_{-\infty}^{\infty} \frac{e^{- \frac{x^{2}}{4 t + 4}}}{\sqrt{t + 1}}\, dx=\frac{\sqrt{\pi} \sqrt{4 t + 4}}{\sqrt{t + 1}}\\
M_2&=\int\limits_{-\infty}^{\infty} \frac{x^{2} e^{- \frac{x^{2}}{4 t + 4}}}{\sqrt{t + 1}}\, dx=\frac{\sqrt{\pi} \left(4 t + 4\right)^{\frac{3}{2}}}{2 \sqrt{t + 1}}\\
\text{全量}&=\frac{\sqrt{\pi} \sqrt{4 t + 4}}{\sqrt{t + 1}}=2 \sqrt{\pi}\\
\text{一方向の分散}&=2 t + 2=2 \left(t + 1\right)
\end{aligned}
$$

時間微分の符号は $r^2-2d(t)$ の符号である。従って $r=4$ では $d<8$ の間は温度が上がり、$t=7$ で最大になる。一方、中心での最大値は最初から減少する。最大値原理が比較するのは全空間での最大値であり、各点の時間単調性を主張しない。
:::

## Tychonoff の非一意性の例

:::proposition[初期値がゼロの非自明解]{#prop-pde-tychonoff}
一次元全空間で $t>0$ に

$$
u(t,x)=\sum_{k=0}^\infty\frac{g^{(k)}(t)}{(2k)!\,\kappa^k}x^{2k},
\qquad g(t)=e^{-1/t^2}
\label{eq-pde-tychonoff}
$$

と置くと、$u_t=\kappa u_{xx}$ であり、$t=0$ に $0$ と連続に延長できる。しかし $u(t,0)=g(t)>0$ である。
:::

:::proof
複素円 $|z-t|=t/4$ は $|\arg z|\le\arctan(1/3)$ にあり、$|z|\le5t/4$ である。従って $\Re(z^{-2})\ge c/t^2$ を満たす $c>0$ がある。Cauchy の微分公式で $|g^{(k)}(t)|\le k!(4/t)^k e^{-c/t^2}$ である。$(2k)!\ge(k!)^2$ より絶対和は $\exp(-c/t^2+4x^2/(\kappa t))$ 以下である。特に有界な $x$ の範囲で初期値は $0$ に収束する。$t$ が正の下限を持つコンパクト集合上では $g^{(k+1)}$ の評価により、微分した項の和も一様収束する。空間で二階微分して添字を一つずらすと、時間微分と項ごとに一致する。
:::

:::example[最初の二つの係数]
$g'=2t^{-3}e^{-1/t^2}$ なので $u=e^{-1/t^2}+t^{-3}e^{-1/t^2}x^2/\kappa+O(x^4)$ である。原点では $\kappa u_{xx}=2t^{-3}e^{-1/t^2}=u_t$ となる。
:::

:::example[無限遠での増大の帰結]
どの $T>0$ に対しても、この解は $[0,T]\times\R$ 上で $C e^{A|x|^2}$ の形の一様評価を満たさない。満たすとすれば前節の短時間の比較で $u(t,0)=0$ が正の短時間に従い、明示値と矛盾する。
:::

### ゼロ時刻まで全階の滑らかさ

Tychonoff 級数を時間で $m$ 回、空間で $\ell$ 回微分する。空間微分した項の分母は $(2k-\ell)!$、分子には $g^{(k+m)}(t)$ が現れる。$k$ が $\ell$ より大きいとき
$$
\frac{(k+m)!}{(2k-\ell)!}
\le C_{m,\ell}\frac{(1+k)^{m+\ell}}{k!}
$$
である。実際 $(k+m)!/k!\le C_m(1+k)^m$、$(2k)!/(2k-\ell)!\le(2k)^\ell$、$(2k)!\ge(k!)^2$ を用いる。$|x|\le R$ を $R\ge1$ に広げ、$k$ の多項式因子を $C2^k$ で抑えると、微分した級数の絶対和は $C t^{-m}\exp(-c/t^2+C_R/t)$ 以下となる。有限個の小さい $k$ の項も同じ形で抑えられる。従って全階の導関数は有界な空間範囲で $t\downarrow0$ に一様に零へ行く。$t<0$ で零に延長し、時間の差商に平均値の定理を使えば、各階で継ぎ目の微分が零になることを帰納的に確かめられる。

:::example[零の Taylor データと正の温度]
$\kappa=1$、$t=1$ で $u(1,0)=e^{-1}$、$u_{xx}(1,0)=2e^{-1}$ である。これに対し $t=0$ では両者は零であり、任意次数の時間微分も零である。この解は時空間で滑らかだが、$t=0$ の原点の近傍で時間に解析的ではない。解析的なら Taylor 級数が恒等的に零となり、正の小時刻での値と矛盾する。
:::

:::example[有限に止まる級数と零の Taylor データの比較]
一次元で $\kappa=1$、$g=3x^6-x^2$ とし、有限多項式解を構成する。

時間一階の方程式なので $g,Pg,P^2g,\ldots$ を順に作る。空間の多項式に二階の作用素を掛けるとこの場合は次数が少なくとも二つ下がる。有限回の後に零になることを確かめると、級数は多項式に打ち切れる。

$$
\begin{aligned}
P^{0}g&=3 x^{6} - x^{2}\\
P^{1}g&=2 \left(45 x^{4} - 1\right)\\
P^{2}g&=1080 x^{2}\\
P^{3}g&=2160\\
P^{4}g&=0
\end{aligned}
$$

$P^j g$ に $t^j/j!$ を掛けて和を作る。時間微分で階乗が一つ減ることにより、空間作用素を一回作用させた列と一致する。さらに時刻を零に戻すと初期データが得られる。

$$
\begin{aligned}
u(t,x)&=360 t^{3} + 540 t^{2} x^{2} + 2 t \left(45 x^{4} - 1\right) + 3 x^{6} - x^{2}\\
u_t&=\frac{\partial}{\partial t} \left(360 t^{3} + 540 t^{2} x^{2} + 2 t \left(45 x^{4} - 1\right) + 3 x^{6} - x^{2}\right)=2 \left(540 t^{2} + 540 t x^{2} + 45 x^{4} - 1\right)\\
Pu&=2 \left(540 t^{2} + 540 t x^{2} + 45 x^{4} - 1\right)\\
u(0,x)&=3 x^{6} - x^{2},\qquad t=0
\end{aligned}
$$

空間の導関数を測ると、時間の多項式の各係数が局所でどう効くかが分かる。初期関数との差は単に平行移動で置き換えず、以下では導関数と時刻を代入して計算する。

$$
\begin{aligned}
\partial_{x} u&=\frac{\partial}{\partial x} \left(360 t^{3} + 540 t^{2} x^{2} + 2 t \left(45 x^{4} - 1\right) + 3 x^{6} - x^{2}\right)=2 x \left(540 t^{2} + 180 t x^{2} + 9 x^{4} - 1\right)\\
g&=0,\qquad x=0\\
u&=44,\qquad t=\frac{1}{2},x=0\\
u_t&=268,\qquad t=\frac{1}{2},x=0
\end{aligned}
$$

ここでは時刻零での各係数が空間の初期値から順に決まり、しかも列が有限で止まる。Tychonoff の例では逆に時刻零で全階の係数が零なのに正時刻で値が現れる。滑らかさだけでは Taylor 表示との一致を保証できず、解析的であることを別に使う必要がある。
:::

## 有界領域における固有関数展開

斉次 Dirichlet 条件を持つ $-\Delta$ の $L^2$ 正規直交固有系を $\phi_j$、固有値を $\lambda_j>0$ とする。この固有系の存在は [楕円型方程式の弱解](/math/pde/weak-solutions/) で説明する。$\int u\phi_j$ を時間微分し、境界項が消えることを使うと係数は $\dot u_j=-\kappa\lambda_j u_j$ を満たす。従って

$$
u(t)=\sum_j e^{-\kappa\lambda_jt}g_j\phi_j,\qquad
g_j=\int_\Omega g\phi_j\,\dd x .
\label{eq-pde-heat-eigen}
$$

任意の $g\in L^2$ について Parseval により定義でき、優収束で $u(t)\to g$ が $L^2$ で成り立つ。$t\ge\tau>0$ では $\sup_{\lambda>0}\lambda^m e^{-\kappa\lambda\tau}<\infty$ なので、時間微分や $(-\Delta)^m$ も $L^2$ で項別に適用できる。最低固有値が正なので $\norm{u(t)}_2\le e^{-\kappa\lambda_1t}\norm g_2$ である。境界が断熱なら Neumann の定数固有モードが $\lambda=0$ に現れ、その係数は減衰しない。

非斉次項 $f\in C([0,T];L^2(\Omega))$ に対しては、$S(t)$ をこの固有展開の解作用素とすると
$$
u(t)=S(t)g+\int_0^t S(t-s)f(s)\,\dd s
$$
となる。積分は $L^2$ の値を取る積分である。$S(0)=I$ として各係数に積分因子を使えばこの式が得られ、$\norm{S(t)}\le1$ と $L^2$ での強連続性により右辺は時刻に連続である。有限個のモードに切った式は係数方程式を満たす。空間試験関数を Dirichlet Laplacian の定義域に取って、作用素を試験関数側へ移せば、$L^2$ 収束で時間の弱形式に移れる。これにより元の初期値と外力に対する解が得られる。また、この弱形式で各 $\phi_j$ を試験に使えば係数は同じ積分因子で定まり、一意である。

:::example[区間の二つのモード]
$(0,\pi)$ の斉次 Dirichlet 条件で $g=3\sin x-2\sin3x$ なら $u=3e^{-\kappa t}\sin x-2e^{-9\kappa t}\sin3x$ である。第 3 モードと第 1 モードの係数の比は $-(2/3)e^{-8\kappa t}$ である。
:::

:::example[一定の初期値と角での不整合]
$(0,\pi)$ で $g=1$ の係数は、規格化していない正弦基底 $\sin(jx)$ では $2(1-(-1)^j)/(\pi j)$ である。$t>0$ では級数が端点まで一様に収束し、各項が端点で $0$ なので $u(t,x)\to0$（$x\downarrow0$）となる。一方 $g=1$ であるから、各 $t>0$ で $\sup_{0<x<\pi}|u(t,x)-1|\ge1$ であり、初期値には一様収束しない。$L^2$ 収束では端点の値を指定しないので式は有効である。
:::

### 半直線と区間の熱核

半直線の零境界に対して初期値を奇延長し、全空間の熱核で解く。$K_t$ は偶なので奇性を保ち、$x>0$ での解の核は

$$
K_t^D(x,y)=K_t(x-y)-K_t(x+y),\qquad y>0
\label{eq-pde-heat-reflection}
$$

である。Neumann 条件なら偶延長し、二項の差を和に替える。Dirichlet 核は $x,y>0$ で $|x-y|<x+y$ なので正である。$y$ による積分は $1$ 未満となる。実際、前の項の積分は $\int_{-\infty}^{x}K_t(z)\dd z$、後の項の積分は $\int_x^\infty K_t(z)\dd z$ なので、差は $\int_{-x}^{x}K_t(z)\dd z$ である。Neumann 核では二項の和が $1$ になり、温度が一定なら保たれる。

長さ $L$ の区間では、奇延長を両端で繰り返すと $2L$ 周期になる。従って Dirichlet 核は

$$
\sum_{m\in\Z}\big(K_t(x-y+2mL)-K_t(x+y+2mL)\big)
$$

である。これは全直線の積分を周期ごとに分け、奇な反射による符号を付けたものである。正の時刻のコンパクト区間では Gaussian の減衰で級数とその導関数が一様収束する。初期値への収束と境界条件は延長と全空間表示から、固有展開との一致は一意性から従う。境界に達する熱の効果は、項ごとの反射としても正弦モードの減衰としても表せる。

:::example[境界近くで減る熱量]
半直線で初期値 $1$ を置く。Dirichlet の解は $\int_{-x}^x K_t(z)\,\dd z=\operatorname{erf}(x/\sqrt{4\kappa t})$、Neumann の解は定数 $1$ である。前者は $x>0$ を固定して $t$ を大きくすると $0$ に近づく。ここで $\operatorname{erf}(w)=2\pi^{-1/2}\int_0^w e^{-s^2}\,\dd s$ である。
:::

### Robin 条件と熱量の散逸

境界で $\partial_\nu u+\beta u=0$（$\beta\ge0$）とする。このとき $-\Delta$ の形式は $H^1(\Omega)$ で
$$
a(u,v)=\int Du\cdot Dv+\beta\int_{\partial\Omega}uv .
$$
有界 Lipschitz 領域ならトレース評価により有界であり、$a+\inner{\cdot}{\cdot}_{L^2}$ は強圧的である。コンパクトな埋め込みと対称性を使うと、Dirichlet の場合と同様に固有関数展開が得られる。零固有値があると $a(u,u)=0$ であるから各連結成分で $u$ は定数である。$\beta>0$ なら境界項も零にするため定数は零である。$\beta=0$ では定数に対応する零モードが残る。滑らかな解の総熱量の微分は $\kappa\int_{\partial\Omega}\partial_\nu u=-\kappa\beta\int_{\partial\Omega}u$ なので、非負解では Robin 境界から熱が失われる。

:::example[区間の偶モード]
区間 $(-1,1)$ で $X(x)=\cos(kx)$ を取る。両端条件は $-k\sin k+\beta\cos k=0$ になり、$k\tan k=\beta$ である。$\beta=\pi/4$ なら $k=\pi/4$ がこれを満たし、$\lambda=\pi^2/16$ である。$u=e^{-\kappa\pi^2t/16}\cos(\pi x/4)$ の総熱量は $e^{-\kappa\pi^2t/16}2\sin k/k$ で、その微分は $-2\kappa k\sin k\,e^{-\kappa k^2t}$、境界での損失は $-2\kappa\beta\cos k\,e^{-\kappa k^2t}$ なので一致する。
:::

:::example[二つの熱モードと一次の源]
$(0,\pi)$、$\kappa=1$、零境界で、初期値 $\sin x+\sin2x$、外力 $2e^{-t}\sin x$ とする。第1モードは外力と自然な減衰の指数が同じなので、時間に比例する係数が現れる。

まず初期値を固有関数ごとに分け、対応する時間方程式を解く。外力がある場合も、同じ固有関数に属する外力だけを時間方程式に入れる。以下では得られた表示の微分を取って、式の各項の関係を確かめる。

$$
\begin{aligned}
u&=\left(2 t + 1\right) e^{- t} \sin{\left(x \right)} + e^{- 4 t} \sin{\left(2 x \right)}\\
u_t&=\frac{\partial}{\partial t} \left(\left(2 t + 1\right) e^{- t} \sin{\left(x \right)} + e^{- 4 t} \sin{\left(2 x \right)}\right)=- \left(2 t e^{3 t} \sin{\left(x \right)} - e^{3 t} \sin{\left(x \right)} + 4 \sin{\left(2 x \right)}\right) e^{- 4 t}\\
u_{xx}&=\frac{\partial^{2}}{\partial x^{2}} \left(\left(2 t + 1\right) e^{- t} \sin{\left(x \right)} + e^{- 4 t} \sin{\left(2 x \right)}\right)=- \left(2 t e^{3 t} \sin{\left(x \right)} + e^{3 t} \sin{\left(x \right)} + 4 \sin{\left(2 x \right)}\right) e^{- 4 t}\\
f&=- \left(2 t + 1\right) e^{- t} \sin{\left(x \right)} + \left(\left(2 t + 1\right) \sin{\left(x \right)} + 4 e^{- 3 t} \sin{\left(2 x \right)}\right) e^{- t} + 2 e^{- t} \sin{\left(x \right)} - 4 e^{- 4 t} \sin{\left(2 x \right)}=2 e^{- t} \sin{\left(x \right)}
\end{aligned}
$$

境界条件を各モードで確かめてから、時刻零のデータに戻る。端点値がすべて零でも時間の係数は任意に選べない。時刻を零に戻して初めて各係数の初期条件が確認できる。

$$
\begin{aligned}
u(0,x)&=\sin{\left(x \right)} + \sin{\left(2 x \right)},\qquad t=0\\
u(t,0)&=0,\qquad x=0\\
u(t,L)&=0,\qquad x=\pi
\end{aligned}
$$

エネルギーを空間積分で測る。異なるモードの積の積分は零であり、一つずつ積分して足し合わせられる。その時間微分を取り、解の式を使って積分する方法と比較する。

$$
\begin{aligned}
e(t,x)&=\frac{\left(\left(2 t + 1\right) e^{- t} \sin{\left(x \right)} + e^{- 4 t} \sin{\left(2 x \right)}\right)^{2}}{2}\\
\mathcal E(t)&=\frac{\pi \left(4 t^{2} e^{6 t} + 4 t e^{6 t} + e^{6 t} + 1\right) e^{- 8 t}}{4}\\
\mathcal E^{\prime}&=\frac{d}{d t} \frac{\pi \left(4 t^{2} e^{6 t} + 4 t e^{6 t} + e^{6 t} + 1\right) e^{- 8 t}}{4}=- \frac{\pi \left(4 t^{2} e^{6 t} - e^{6 t} + 4\right) e^{- 8 t}}{2}\\
\text{入力と散逸の和}&=\int\limits_{0}^{\pi} \left(2 \left(\left(2 t + 1\right) e^{- t} \sin{\left(x \right)} + e^{- 4 t} \sin{\left(2 x \right)}\right) e^{- t} \sin{\left(x \right)} - \left(\left(2 t + 1\right) e^{- t} \cos{\left(x \right)} + 2 e^{- 4 t} \cos{\left(2 x \right)}\right)^{2}\right)\, dx=- \frac{\pi \left(4 t^{2} e^{6 t} - e^{6 t} + 4\right) e^{- 8 t}}{2}
\end{aligned}
$$

最後に一つの観測点で変位または温度と勾配を計算する。点で値が相殺されることと、区間全体のエネルギーが零になることは区別する。

$$
\begin{aligned}
u&=\frac{1}{2} + \log{\left(2 \right)}=\frac{1 + 2 \log{\left(2 \right)}}{2},\qquad t=\log{\left(2 \right)},x=\frac{\pi}{2}\\
u_x&=- \frac{1}{8},\qquad t=\log{\left(2 \right)},x=\frac{\pi}{2}\\
\mathcal E&=\frac{\pi \left(65 + 256 \log{\left(2 \right)}^{2} + 256 \log{\left(2 \right)}\right)}{1024},\qquad t=\log{\left(2 \right)}
\end{aligned}
$$

第1モードの係数は $(e^t u_1)^{\prime}=2$ から、第2モードは $(e^{4t}u_2)^{\prime}=0$ から得られる。観測点 $x=\pi/2$ で第2モードの値は零だが、その導関数は零でない。従って一点で温度だけを観測すると、初期データの第2モードは見えない。
:::

:::example[同じ初期値から吸収と加熱を比較する]
長さ $\pi$ の区間で $\kappa=2$、反応項 $+u$、初期値 $\sin2x$、源 $18\sin2x$ とする。第2モードの時間応答を求める。

まず $u(t,x)=a(t)\sin(mx)$ を代入する。空間二階の微分は $-m^2$ 倍であり、時間微分は $a^{\prime}$ だけに作用する。従って時間依存係数を含む場合でも、他の空間モードを結合させなければ一次の方程式に帰着する。ここで $\mu(t)$ を空間作用素のこのモードでの係数とする。

$$
\begin{aligned}
m&=2\\
\mu(t)&=9\\
f_m&=18\\
M(t) &= \int_0^t \mu(s)\,\mathrm{d}s=\int\limits_{0}^{t} 9\, dz=9 t\\
I(t)&=e^{M(t)}=e^{9 t}
\end{aligned}
$$

積分因子を掛けると $(Ia)^{\prime}=If_m$ である。零時刻から積分し、初期係数を加えてから因子で割る。外力を単に時間積分するのと違い、積分の中で指数の重みが掛かる。

$$
\begin{aligned}
If_m&=18 e^{9 t}\\
\int_0^t I(s)f_m(s)\,\mathrm{d}s&=\int\limits_{0}^{t} 18 e^{9 z}\, dz=2 e^{9 t} - 2\\
a(t)&=\left(2 e^{9 t} - 1\right) e^{- 9 t}\\
a^{\prime}&=\frac{d}{d t} \left(2 e^{9 t} - 1\right) e^{- 9 t}=9 e^{- 9 t}\\
a^{\prime}+\mu a&=18\\
a(0)&=1,\qquad t=0
\end{aligned}
$$

正弦の二乗積分は $\pi/2$ である。解の二乗ノルムの半分を時間微分し、散逸および入力から直接求めた変化率と比べる。係数が時間に依存しても、ここで微分する二乗ノルム自体には時間依存の作用素が入らない。

$$
\begin{aligned}
u(t,x)&=\left(2 e^{9 t} - 1\right) e^{- 9 t} \sin{\left(2 x \right)}\\
E(t)&=\frac{\pi \left(2 e^{9 t} - 1\right)^{2} e^{- 18 t}}{4}\\
E^{\prime}&=\frac{d}{d t} \frac{\pi \left(2 e^{9 t} - 1\right)^{2} e^{- 18 t}}{4}=\frac{9 \pi \left(2 e^{9 t} - 1\right) e^{- 18 t}}{2}\\
\frac{\pi}2 a(f_m-\mu a)&=\frac{9 \pi \left(2 e^{9 t} - 1\right) e^{- 18 t}}{2}\\
u(t_*,x_*)&=\frac{3}{2},\qquad t=\frac{\log{\left(2 \right)}}{9},x=\frac{\pi}{4}\\
E(t_*)&=\frac{9 \pi}{16},\qquad t=\frac{\log{\left(2 \right)}}{9}
\end{aligned}
$$

このモードでは拡散が $8$、吸収が $1$ なので、減衰率は $9$ である。定常係数は $18/9=2$、観測時刻は差を半分に減らす時刻なので係数は $3/2$ である。初期に正負があるため、これは一点の値を読んだ結果であり、全区間で温度が非負になるわけではない。
:::

## エネルギー法と後退一意性

:::proposition[エネルギー散逸と後退一意性]{#prop-pde-heat-backward}
斉次 Dirichlet 境界の古典解に対し

$$
\frac{\dd}{\dd t}\frac12\norm{u(t)}_2^2=-\kappa\norm{Du(t)}_2^2 .
\label{eq-pde-heat-energy}
$$

また $L^2$ 固有関数展開で定まる解が時刻 $T>0$ で $0$ なら、初期値も $0$ である。ただし逆写像は $L^2$ ノルムで連続でない。
:::

:::proof
$u$ を掛けて積分し、Green の公式で境界項 $u\partial_\nu u$ を消せばエネルギー式を得る。後退一意性は各固有係数で $e^{-\kappa\lambda_jT}g_j=0$ が $g_j=0$ を意味することから従う。一方、$g=\phi_j$ はノルム $1$ であり、$u(T)$ のノルムは $e^{-\kappa\lambda_jT}\to0$ であるため、逆写像は連続でない。
:::

別の証明で時間の凸性が見える。$A=-\kappa\Delta$、$F=(u,u)>0$ とすると、対称性から $F'=-2(Au,u)$、$F''=4(Au,Au)$ である。従って $F''F-(F')^2\ge0$ であり、$\log F$ は凸である。$t_1<t<t_2$ で $F(t)\le F(t_1)^{(t_2-t)/(t_2-t_1)}F(t_2)^{(t-t_1)/(t_2-t_1)}$ を得る。初期に $F>0$ ならこの不等式の右端を最初の零点へ近づけると、途中の固定時刻の $F$ が $0$ と矛盾する。

:::example[エネルギーの減衰率]
$(0,\pi)$ で $u=e^{-4t}\sin2x$ とする。$\norm{u}_2^2=\pi e^{-8t}/2$、$\norm{u_x}_2^2=2\pi e^{-8t}$ なので $\frac12(\norm{u}_2^2)'=-2\pi e^{-8t}$ と一致する。
:::

:::example[逆向きの誤差の増幅]
$(0,\pi)$ で $\kappa=1$、$T=1$ の終値に $e^{-100}\sin10x$ の誤差があると、それを逆にたどった初期誤差は $\sin10x$ になる。周波数がさらに高くなると終値の誤差を任意に小さくできるため、終値からの復元は安定しない。
:::

### 周波数を切って後退問題を安定化する

時刻 $T$ の観測値が $y=S(T)g+\eta$、$\norm\eta_2\le\delta$ とする。$A=-\kappa\Delta$ の固有値を $\mu_j$ とし、$\mu_j\le\Lambda$ のモードだけで $g_\Lambda=\sum_{\mu_j\le\Lambda}e^{T\mu_j}(y,\phi_j)\phi_j$ と置く。低周波の誤差ノルムは Parseval により $e^{T\Lambda}\delta$ 以下で、高周波の切捨ては $\norm{A^s g}_2\le M$（$s>0$）なら $M\Lambda^{-s}$ 以下である。二つは異なるモードにあるため二乗ノルムは和になる。とくに $0<\delta<M$ に対して $\Lambda=(2T)^{-1}\log(M/\delta)$ とすれば低周波の誤差は $\sqrt{M\delta}$、高周波は対数の負のべきで減る。従って条件付きの安定化は可能だが、全周波を同時に逆にすることはできない。

:::example[二つのモードから一つだけ戻す]
$(0,\pi)$ で $\kappa=T=1$、正規化した基底を $\phi_j$ とし、$g=\phi_1+\phi_3/9$ とする。$\Lambda=4$ では第 $1,2$ モードだけ戻す。誤差が $\eta=\delta\phi_2$ なら復元値は $\phi_1+e^4\delta\phi_2$ であり、誤差は $\sqrt{e^8\delta^2+1/81}$ である。$\Lambda=9$ に上げればこの例の切捨ては消えるが、第 $3$ モードに誤差があれば $e^9$ 倍される。切捨ての選択は信号と誤差の両方に依存する。
:::

:::example[熱エネルギーが一つの指数にならない場合]
$(0,\pi)$ で $\kappa=2$、外力零、初期値 $\sin x+3\sin(2x)$ をとる。

まず初期値を固有関数ごとに分け、対応する時間方程式を解く。外力がある場合も、同じ固有関数に属する外力だけを時間方程式に入れる。以下では得られた表示の微分を取って、式の各項の関係を確かめる。

$$
\begin{aligned}
u&=e^{- 2 t} \sin{\left(x \right)} + 3 e^{- 8 t} \sin{\left(2 x \right)}\\
u_t&=\frac{\partial}{\partial t} \left(e^{- 2 t} \sin{\left(x \right)} + 3 e^{- 8 t} \sin{\left(2 x \right)}\right)=- 2 \left(e^{6 t} + 24 \cos{\left(x \right)}\right) e^{- 8 t} \sin{\left(x \right)}\\
u_{xx}&=\frac{\partial^{2}}{\partial x^{2}} \left(e^{- 2 t} \sin{\left(x \right)} + 3 e^{- 8 t} \sin{\left(2 x \right)}\right)=- \left(e^{6 t} + 24 \cos{\left(x \right)}\right) e^{- 8 t} \sin{\left(x \right)}\\
f&=2 \left(\sin{\left(x \right)} + 12 e^{- 6 t} \sin{\left(2 x \right)}\right) e^{- 2 t} - 2 e^{- 2 t} \sin{\left(x \right)} - 24 e^{- 8 t} \sin{\left(2 x \right)}=0
\end{aligned}
$$

境界条件を各モードで確かめてから、時刻零のデータに戻る。端点値がすべて零でも時間の係数は任意に選べない。時刻を零に戻して初めて各係数の初期条件が確認できる。

$$
\begin{aligned}
u(0,x)&=\sin{\left(x \right)} + 3 \sin{\left(2 x \right)}=\left(6 \cos{\left(x \right)} + 1\right) \sin{\left(x \right)},\qquad t=0\\
u(t,0)&=0,\qquad x=0\\
u(t,L)&=0,\qquad x=\pi
\end{aligned}
$$

エネルギーを空間積分で測る。異なるモードの積の積分は零であり、一つずつ積分して足し合わせられる。その時間微分を取り、解の式を使って積分する方法と比較する。

$$
\begin{aligned}
e(t,x)&=\frac{\left(e^{- 2 t} \sin{\left(x \right)} + 3 e^{- 8 t} \sin{\left(2 x \right)}\right)^{2}}{2}\\
\mathcal E(t)&=\frac{\pi \left(e^{12 t} + 9\right) e^{- 16 t}}{4}\\
\mathcal E^{\prime}&=\frac{d}{d t} \frac{\pi \left(e^{12 t} + 9\right) e^{- 16 t}}{4}=- \pi \left(e^{12 t} + 36\right) e^{- 16 t}\\
\text{入力と散逸の和}&=\int\limits_{0}^{\pi} \left(- 2 \left(e^{- 2 t} \cos{\left(x \right)} + 6 e^{- 8 t} \cos{\left(2 x \right)}\right)^{2}\right)\, dx=- \pi \left(e^{12 t} + 36\right) e^{- 16 t}
\end{aligned}
$$

最後に一つの観測点で変位または温度と勾配を計算する。点で値が相殺されることと、区間全体のエネルギーが零になることは区別する。

$$
\begin{aligned}
u&=\frac{11 \sqrt{3}}{32},\qquad t=\frac{\log{\left(2 \right)}}{2},x=\frac{\pi}{3}\\
u_x&=\frac{1}{16},\qquad t=\frac{\log{\left(2 \right)}}{2},x=\frac{\pi}{3}\\
\mathcal E&=\frac{73 \pi}{1024},\qquad t=\frac{\log{\left(2 \right)}}{2}
\end{aligned}
$$

初期には第2モードのエネルギーが第1の9倍である。一方、正時刻での比は $9e^{-12t}$ となり、$t>(\log9)/12$ では低いモードが優勢になる。対数エネルギーの傾きが一定でないのは、異なる減衰率の重みの比が変わるからである。
:::

:::example[負の外力がエネルギーをさらに減らす]
$(0,\pi)$ で $\kappa=1$、初期値 $3\sin x$、源 $-2e^{-t}\sin x$ とする。$t=1$ でのエネルギー収支を測る。

まず $u(t,x)=a(t)\sin(mx)$ を代入する。空間二階の微分は $-m^2$ 倍であり、時間微分は $a^{\prime}$ だけに作用する。従って時間依存係数を含む場合でも、他の空間モードを結合させなければ一次の方程式に帰着する。ここで $\mu(t)$ を空間作用素のこのモードでの係数とする。

$$
\begin{aligned}
m&=1\\
\mu(t)&=1\\
f_m&=- 2 e^{- t}\\
M(t) &= \int_0^t \mu(s)\,\mathrm{d}s=\int\limits_{0}^{t} 1\, dz=t\\
I(t)&=e^{M(t)}=e^{t}
\end{aligned}
$$

積分因子を掛けると $(Ia)^{\prime}=If_m$ である。零時刻から積分し、初期係数を加えてから因子で割る。外力を単に時間積分するのと違い、積分の中で指数の重みが掛かる。

$$
\begin{aligned}
If_m&=-2\\
\int_0^t I(s)f_m(s)\,\mathrm{d}s&=\int\limits_{0}^{t} \left(-2\right)\, dz=- 2 t\\
a(t)&=- \left(2 t - 3\right) e^{- t}\\
a^{\prime}&=\frac{d}{d t} \left(- \left(2 t - 3\right) e^{- t}\right)=\left(2 t - 5\right) e^{- t}\\
a^{\prime}+\mu a&=- 2 e^{- t}\\
a(0)&=3,\qquad t=0
\end{aligned}
$$

正弦の二乗積分は $\pi/2$ である。解の二乗ノルムの半分を時間微分し、散逸および入力から直接求めた変化率と比べる。係数が時間に依存しても、ここで微分する二乗ノルム自体には時間依存の作用素が入らない。

$$
\begin{aligned}
u(t,x)&=- \left(2 t - 3\right) e^{- t} \sin{\left(x \right)}\\
E(t)&=\frac{\pi \left(2 t - 3\right)^{2} e^{- 2 t}}{4}\\
E^{\prime}&=\frac{d}{d t} \frac{\pi \left(2 t - 3\right)^{2} e^{- 2 t}}{4}=- \frac{\pi \left(2 t - 5\right) \left(2 t - 3\right) e^{- 2 t}}{2}\\
\frac{\pi}2 a(f_m-\mu a)&=- \frac{\pi \left(2 t - 5\right) \left(2 t - 3\right) e^{- 2 t}}{2}\\
u(t_*,x_*)&=e^{-1},\qquad t=1,x=\frac{\pi}{2}\\
E(t_*)&=\frac{\pi}{4 e^{2}},\qquad t=1
\end{aligned}
$$

$0\le t<3/2$ では係数が正、源の係数が負なので、拡散による散逸と源による仕事は共に負である。$t=3/2$ に解は全区間で零になるが、外力はその後も続くので負の係数が現れる。後退一意性は同次式の性質であり、既知の外力がない問題と混同できない。
:::

### 全空間でも Fourier 乗数は零にならない

$g\in L^2(\R^n)$ とし、Fourier 変換を :ref[def-four-euclid-transform] の規約で取る。$u(T)=K_T*g$ なら $\widehat u(T,\xi)=e^{-4\pi^2\kappa T|\xi|^2}\widehat g(\xi)$ である。乗数はどの有限周波数でも正なので、終値が零なら $\widehat g=0$ a.e.、Plancherel から $g=0$ である。それでも逆の有界性はない。実際、各 $R$ に対し $R\le|\xi|\le R+1$ に台を持ち $L^2$ ノルムが $1$ の Fourier データを選べば
$$
\norm{K_T*g_R}_2\le e^{-4\pi^2\kappa T R^2}
$$
である。従って固有値が離散でない全空間でも、高周波を選ぶだけで同じ結論が出る。一方、周波数を $|\xi|\le R_0$ に制限したクラス上では、乗数は $e^{-4\pi^2\kappa T R_0^2}$ 以上であり、逆の作用素ノルムはその逆数以下になる。制限された帯域での安定性は、全空間での後退問題の適切性とは異なる。

## 演習問題

::::exercise{level=basic}
一次元で $\kappa=2$ の熱核の分散を時刻 $3$ に求めよ。

:::solution
平均は対称性から $0$ である。$xK_t=-2\kappa t\,\partial_xK_t$ を用いて $\int x^2K_t=2\kappa t$ なので、分散は $12$ である。
:::
::::

::::exercise{level=basic}
$(0,\pi)$ の斉次 Dirichlet 問題で $\kappa=1$、$g=2\sin3x$ の解を求めよ。

:::solution
$-\partial_{xx}\sin3x=9\sin3x$ なので $u=2e^{-9t}\sin3x$ である。初期条件と両端の条件も満たす。
:::
::::

::::exercise{level=basic}
$\kappa=1$ で空間に一定の条件 $g=2$、$f(t,x)=3t$ を課す。解を求めよ。

:::solution
熱核は定数を保つので :eqref[eq-pde-heat-duhamel] は $u=2+\int_0^t3s\,\dd s=2+3t^2/2$ を与える。空間 Laplace 作用素は $0$ であり $u_t=3t$ である。
:::
::::

::::exercise{level=basic}
$E_r$ に現れる過去への時間差の上限を $r=4$、$\kappa=2$ で求めよ。

:::solution
$K_s(y)\le K_s(0)$ なので条件は $s<r^2/(4\pi\kappa)$ である。値を入れて上限は $2/\pi$ になる。
:::
::::

::::exercise{level=advanced}
$g(x)=x^2$ の熱核との畳み込みを一次元で求めよ。

:::solution
$\int K_t(y)(x-y)^2\,\dd y=x^2+2\kappa t$ である。交差項は奇関数の積分で消え、$y^2$ の積分は核の分散である。$u_t=2\kappa=\kappa u_{xx}$ である。
:::
::::

::::exercise{level=advanced}
区間 $(0,\pi)$ の斉次 Dirichlet 条件で $u_t=u_{xx}+\sin x$、初期値 $0$ を解け。

:::solution
右辺の空間モードは $\sin x$ のみなので、係数は $a'+a=1$、$a(0)=0$ を満たす。よって $u=(1-e^{-t})\sin x$ である。Duhamel 表示でも $\int_0^t e^{-(t-s)}\,\dd s=1-e^{-t}$ になる。
:::
::::

::::exercise{level=advanced}
熱核の $\norm{K_t}_{L^2(\R^n)}^2$ を畳み込みで求めよ。

:::solution
偶性により $\int K_t(y)^2\,\dd y=(K_t*K_t)(0)=K_{2t}(0)=(8\pi\kappa t)^{-n/2}$ である。時間が増えるにつれて質量は変わらず $L^2$ ノルムは減少する。
:::
::::

::::exercise{level=advanced}
$t>0$ の Tychonoff の例で $x^4$ の係数を計算し、前の係数との熱方程式の関係を確かめよ。

:::solution
$g''=e^{-1/t^2}(4t^{-6}-6t^{-4})$ なので係数は $g''/(24\kappa^2)$ である。この項を $x$ で二階微分し $\kappa$ を掛けると $x^2$ の係数は $g''/(2\kappa)$ である。これは $x^2$ の元の係数 $g'/(2\kappa)$ の時間微分と一致する。
:::
::::

## 確認テスト

::::quiz{answer="8" tol="1e-6"}
一次元で $\kappa=2$、$t=2$ の熱核の分散を求めよ。

:::explanation
分散は $2\kappa t$ なので $2\cdot2\cdot2=8$ である。
:::
::::

::::quiz
有界領域の熱方程式の最大値原理で使う境界はどれか。

- [x] 初期面と側面
- [ ] 終時刻と側面
- [ ] 側面だけ

:::explanation
内部最大点では $\Delta u\le0$ で、過去からの時間微分は非負である。終時刻の内部でも左微分を使える。従って $u-\eps t$ の最大は初期面と側面に限られる。
:::
::::

::::quiz{answer="9" tol="1e-6"}
$\kappa=1$ の区間 $(0,\pi)$ で $\sin3x$ モードの減衰係数を $e^{-\alpha t}$ と書く。$\alpha$ を求めよ。

:::explanation
$\partial_{xx}\sin3x=-9\sin3x$ であるため係数は $e^{-9t}$ であり、$\alpha=9$ である。
:::
::::

::::quiz
終値から初期値を求める熱方程式について成立するものはどれか。

- [x] 固有関数展開の解では終値が等しければ初期値も等しい
- [ ] $L^2$ での観測誤差は逆向きでも増幅しない
- [ ] 高周波ほど前向きの減衰が小さい

:::explanation
各モードの係数 $e^{-\kappa\lambda_jT}$ は正なので一意性がある。しかし $\lambda_j\to\infty$ でこの係数は $0$ に近づく。逆向きに割り算をすると高周波成分の誤差が増幅する。
:::
::::

::::quiz{answer="14" tol="1e-6"}
一次元の $u_t=3u_{xx}$ で初期値 $u(0,x)=x^2$ とする。二次の多項式解について $u(5/3,2)$ を求めよ。

:::explanation
$u=x^2+at$ と置けば $u_{xx}=2$ なので $a=6$ である。この式に $t=5/3,x=2$ を代入すると $4+6(5/3)=14$ である。多項式の係数比較では空間変化の二階が定数を作る。
:::
::::


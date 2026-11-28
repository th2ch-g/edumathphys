---
title: 解の表示のための諸方法
description: 境界条件、対称性、変換を使って偏微分方程式を簡約し、短波長と細かな周期構造の極限を考える。
notation:
  - symbol: '\widehat u'
    meaning: 空間変数についての Fourier 変換
  - symbol: '\widetilde u'
    meaning: 時間についての Laplace 変換
  - symbol: '\bar a'
    meaning: 一次元の周期係数に対応する有効係数
---

## この章の目標

- 境界条件から分離定数と展開係数を求める。
- 相似形と進行波を代入して常微分方程式を作る。
- Fourier、Laplace、Cole–Hopf、Legendre の各変換で解を求める。
- 幾何光学の位相と振幅の式を導き、周期媒質の有効量を計算する。

## 変数分離法

:term[変数分離法]{en="separation of variables" yomi="へんすうぶんりほう"} では、未知関数を独立変数ごとの積として探す。この形は解全体を表す仮定ではない。線形方程式なら個々の積の解を重ね合わせ、初期値や境界値に合わせられる。

区間 $(0,L)$ の熱方程式に $u=T(t)X(x)$ を代入すると、非自明な積に対し $T'=-\kappa\lambda T$、$-X''=\lambda X$ である。Dirichlet 両端条件から $\lambda=(m\pi/L)^2$、$X=\sin(m\pi x/L)$、$m\ge1$ となる。実際、$\lambda\le0$ なら $\int|X'|^2=\lambda\int |X|^2$ で $X=0$ となり、正の $\lambda$ では $X(0)=0$ から正弦のみが残る。係数は $g$ の正弦 Fourier 係数 $2L^{-1}\int_0^Lg(x)\sin(m\pi x/L)\,\dd x$ である。$L^2$ 完全性と収束は [Fourier 級数の応用](/math/fourier-analysis/l2-fourier-series/) と :eqref[eq-pde-heat-eigen] を用いる。

:::example[四辺のうち一辺から与える値]
矩形 $(0,\pi)\times(0,a)$ の上辺で $g(x)=2\sin3x$、他の辺で $0$ を与える。$\Delta u=0$ を $X(x)Y(y)$ で分離すると $X=\sin3x$、$Y''=9Y$ である。底辺で $Y=0$、上辺で $Y=2$ なので $u=2\sinh(3y)\sin(3x)/\sinh(3a)$ である。上辺が一般の $g$ なら各 $m$ について上端の係数を正弦係数に替える。高さ $y<a$ での比 $\sinh(my)/\sinh(ma)$ は高周波を指数的に減衰させる。
:::

一般の連続な $g$ について閉矩形までの連続性を要求するなら $g(0)=g(\pi)=0$ とする。この条件の下では奇かつ $2\pi$ 周期に延長した $g$ が連続であり、:ref[thm-four-fejer] により正弦多項式で一様近似できる。各多項式の調和拡張は有限和であり、:ref[cor-pde-laplace-dirichlet] によって解の差の一様ノルムは境界データの差以下である。従って閉矩形上の一様極限が存在し、平均値公式から内部で調和になる。任意の内点では指数減衰により通常の分離級数も絶対収束し、Fejér 平均の係数を極限へ送ると同じ値を得る。

:::example[極座標での分離]
円盤の $\Delta u=u_{rr}+r^{-1}u_r+r^{-2}u_{\theta\theta}$ に $R(r)\Theta(\theta)$ を入れる。$2\pi$ 周期性から $\Theta''=-m^2\Theta$、$m\in\{0,1,\dots\}$ である。$m>0$ に対する半径の解は $r^m,r^{-m}$ で、原点まで有界にするには前者を選ぶ。$m=0$ では $1,\log r$ から定数を選ぶ。従って単位円盤で $g=1+3\cos2\theta$ の拡張は $1+3r^2\cos2\theta$ である。
:::

### 同次化してから固有関数を使う

不均一な境界値を直ちに正弦級数で展開すると、各項の境界値が零なので、境界まで一様に収束する場合には和の境界値も零になる。例えば区間 $(0,L)$ で $u(t,0)=a(t)$、$u(t,L)=b(t)$ なら、$G(t,x)=a(t)(1-x/L)+b(t)x/L$ を引く。$v=u-G$ は零境界となり、$v_t-\kappa v_{xx}=f-G_t$ である。空間について $G_{xx}=0$ なので追加された外力の係数は
$$
\frac2L\int_0^L G_t(t,x)\sin(m\pi x/L)\dd x
=\frac2{m\pi}(a'(t)-(-1)^m b'(t)).
$$
これは定数項と $x$ の項を一回ずつ部分積分すれば得られる。従って同次化により境界入力を各モードの強制項へ移せる。例えば $L=\pi,\kappa=1$、$a=t$、$b=0$、$f=0$、$u(0)=0$ なら $v_m=-2(1-e^{-m^2t})/(\pi m^3)$ である。時間積分は $v_m'+m^2v_m=-2/(m\pi)$ を解いた結果であり、各境界値は加え戻した $G$ によって実現する。

:::example[長さの異なる区間で分離定数を規格化する]
$(0,2\pi)$ で拡散率 $4$、初期値 $2\sin(x/2)+\sin(3x/2)$ とする。区間の長さの効果を分離定数と時間係数で確かめる。

まず初期値を固有関数ごとに分け、対応する時間方程式を解く。外力がある場合も、同じ固有関数に属する外力だけを時間方程式に入れる。以下では得られた表示の微分を取って、式の各項の関係を確かめる。

$$
\begin{aligned}
u&=2 e^{- t} \sin{\left(\frac{x}{2} \right)} + e^{- 9 t} \sin{\left(\frac{3 x}{2} \right)}\\
u_t&=\frac{\partial}{\partial t} \left(2 e^{- t} \sin{\left(\frac{x}{2} \right)} + e^{- 9 t} \sin{\left(\frac{3 x}{2} \right)}\right)=- \left(2 e^{8 t} \sin{\left(\frac{x}{2} \right)} + 9 \sin{\left(\frac{3 x}{2} \right)}\right) e^{- 9 t}\\
u_{xx}&=\frac{\partial^{2}}{\partial x^{2}} \left(2 e^{- t} \sin{\left(\frac{x}{2} \right)} + e^{- 9 t} \sin{\left(\frac{3 x}{2} \right)}\right)=- \frac{\left(2 e^{8 t} \sin{\left(\frac{x}{2} \right)} + 9 \sin{\left(\frac{3 x}{2} \right)}\right) e^{- 9 t}}{4}\\
f&=\left(2 \sin{\left(\frac{x}{2} \right)} + 9 e^{- 8 t} \sin{\left(\frac{3 x}{2} \right)}\right) e^{- t} - 2 e^{- t} \sin{\left(\frac{x}{2} \right)} - 9 e^{- 9 t} \sin{\left(\frac{3 x}{2} \right)}=0
\end{aligned}
$$

境界条件を各モードで確かめてから、時刻零のデータに戻る。端点値がすべて零でも時間の係数は任意に選べない。時刻を零に戻して初めて各係数の初期条件が確認できる。

$$
\begin{aligned}
u(0,x)&=2 \sin{\left(\frac{x}{2} \right)} + \sin{\left(\frac{3 x}{2} \right)},\qquad t=0\\
u(t,0)&=0,\qquad x=0\\
u(t,L)&=0,\qquad x=2 \pi
\end{aligned}
$$

エネルギーを空間積分で測る。異なるモードの積の積分は零であり、一つずつ積分して足し合わせられる。その時間微分を取り、解の式を使って積分する方法と比較する。

$$
\begin{aligned}
e(t,x)&=\frac{\left(2 e^{- t} \sin{\left(\frac{x}{2} \right)} + e^{- 9 t} \sin{\left(\frac{3 x}{2} \right)}\right)^{2}}{2}\\
\mathcal E(t)&=\frac{\pi \left(2 e^{8 t} - 2 e^{4 t} + 1\right) \left(2 e^{8 t} + 2 e^{4 t} + 1\right) e^{- 18 t}}{2}\\
\mathcal E^{\prime}&=\frac{d}{d t} \frac{\pi \left(2 e^{8 t} - 2 e^{4 t} + 1\right) \left(2 e^{8 t} + 2 e^{4 t} + 1\right) e^{- 18 t}}{2}=- \pi \left(4 e^{16 t} + 9\right) e^{- 18 t}\\
\text{入力と散逸の和}&=\int\limits_{0}^{2 \pi} \left(- 4 \left(e^{- t} \cos{\left(\frac{x}{2} \right)} + \frac{3 e^{- 9 t} \cos{\left(\frac{3 x}{2} \right)}}{2}\right)^{2}\right)\, dx=- \pi \left(4 e^{16 t} + 9\right) e^{- 18 t}
\end{aligned}
$$

最後に一つの観測点で変位または温度と勾配を計算する。点で値が相殺されることと、区間全体のエネルギーが零になることは区別する。

$$
\begin{aligned}
u&=\frac{511}{512},\qquad t=\log{\left(2 \right)},x=\pi\\
u_x&=0,\qquad t=\log{\left(2 \right)},x=\pi\\
\mathcal E&=\frac{262145 \pi}{524288},\qquad t=\log{\left(2 \right)}
\end{aligned}
$$

分離定数は $1/4$ と $9/4$ であり、拡散率 $4$ を掛けると減衰率が $1,9$ になる。各正弦の二乗積分は $\pi$ である。区間を倍にすると固有値は4分の1になり、拡散率を4倍にするとこの違いが相殺される。
:::

:::example[円板の外側で成長する枝を除く]
単位円周で $g=3\cos2\theta+\sin5\theta$、無限遠で零になる解を求める。分離定数は内側と同じだが半径方向の条件が違う。

円周データを定数と三角関数に分ける。角振動数 $j$ について半径の方程式は $R^{\prime\prime}+r^{-1}R^{\prime}-j^2r^{-2}R=0$ である。元の条件に合う根を選び、境界半径で係数を規格化する。この例での関数は次のとおりである。

$$
\begin{aligned}
u(r,\theta)&=\frac{3 \cos{\left(2 \theta \right)}}{r^{2}} + \frac{\sin{\left(5 \theta \right)}}{r^{5}}\\
g&=\sin{\left(5 \theta \right)} + 3 \cos{\left(2 \theta \right)},\qquad r=1\\
u_r&=\frac{\partial}{\partial r} \left(\frac{3 \cos{\left(2 \theta \right)}}{r^{2}} + \frac{\sin{\left(5 \theta \right)}}{r^{5}}\right)=- \frac{6 r^{3} \cos{\left(2 \theta \right)} + 5 \sin{\left(5 \theta \right)}}{r^{6}}\\
u_\theta&=\frac{\partial}{\partial \theta} \left(\frac{3 \cos{\left(2 \theta \right)}}{r^{2}} + \frac{\sin{\left(5 \theta \right)}}{r^{5}}\right)=- \frac{6 r^{3} \sin{\left(2 \theta \right)} - 5 \cos{\left(5 \theta \right)}}{r^{5}}
\end{aligned}
$$

二階微分を別々に計算する。半径方向の二階微分だけを角度方向の二階微分と足しても相殺しない。極座標の Laplacian の一階項も入れると各振動数で相殺する。

$$
\begin{aligned}
u_{rr}&=\frac{\partial^{2}}{\partial r^{2}} \left(\frac{3 \cos{\left(2 \theta \right)}}{r^{2}} + \frac{\sin{\left(5 \theta \right)}}{r^{5}}\right)=\frac{6 \left(3 r^{3} \cos{\left(2 \theta \right)} + 5 \sin{\left(5 \theta \right)}\right)}{r^{7}}\\
u_{\theta\theta}&=\frac{\partial^{2}}{\partial \theta^{2}} \left(\frac{3 \cos{\left(2 \theta \right)}}{r^{2}} + \frac{\sin{\left(5 \theta \right)}}{r^{5}}\right)=- \frac{12 r^{3} \cos{\left(2 \theta \right)} + 25 \sin{\left(5 \theta \right)}}{r^{5}}\\
\Delta u&=\frac{- \frac{6 \cos{\left(2 \theta \right)}}{r^{3}} - \frac{5 \sin{\left(5 \theta \right)}}{r^{6}}}{r} + \frac{6 \left(3 \cos{\left(2 \theta \right)} + \frac{5 \sin{\left(5 \theta \right)}}{r^{3}}\right)}{r^{4}} - \frac{12 \cos{\left(2 \theta \right)} + \frac{25 \sin{\left(5 \theta \right)}}{r^{3}}}{r^{4}}=0
\end{aligned}
$$

境界における値の平均と法線微分を読む。$j\ne0$ の角方向の平均は零であり、全体の平均は定数項である。またエネルギーは体積積分の方法と Green の境界積分の方法の両方で計算できる。

$$
\begin{aligned}
\text{円周の平均}&=\int\limits_{- \pi}^{\pi} \frac{\sin{\left(5 \theta \right)} + 3 \cos{\left(2 \theta \right)}}{2 \pi}\, d\theta=0\\
\partial_\nu u&=5 \sin{\left(5 \theta \right)} + 6 \cos{\left(2 \theta \right)}\\
\int_{-\pi}^{\pi}|Du|^2\,\mathrm{d}\theta&=\frac{72 \pi}{r^{6}} + \frac{50 \pi}{r^{12}}=\frac{2 \pi \left(36 r^{6} + 25\right)}{r^{12}}\\
\int |Du|^2&=\int\limits_{1}^{\infty} r \left(\frac{72 \pi}{r^{6}} + \frac{50 \pi}{r^{12}}\right)\, dr=23 \pi\\
\int_{\partial\Omega} u\partial_\nu u&=\int\limits_{- \pi}^{\pi} \left(\sin{\left(5 \theta \right)} + 3 \cos{\left(2 \theta \right)}\right) \left(5 \sin{\left(5 \theta \right)} + 6 \cos{\left(2 \theta \right)}\right)\, d\theta=23 \pi
\end{aligned}
$$

次に同じ振動数の半径因子を比較し、角度を固定した一点で値と半径方向の傾きを調べる。角度の取り方によって、振幅のあるモードでも点での値が零になることがある。

$$
\begin{aligned}
u&=- \frac{\sqrt{2}}{64},\qquad r=2,\theta=\frac{\pi}{4}\\
u_r&=\frac{5 \sqrt{2}}{128},\qquad r=2,\theta=\frac{\pi}{4}\\
u_\theta&=- \frac{3}{2} - \frac{5 \sqrt{2}}{64}=- \frac{5 \sqrt{2} + 96}{64},\qquad r=2,\theta=\frac{\pi}{4}
\end{aligned}
$$

外部では負の指数の半径因子を使う。例えば半径 $2$ では二次モードは $3/4$、五次は $1/32$ の振幅になる。外部領域の境界での法線は原点に向くため、エネルギー式に入れる微分の符号も内側の場合と逆である。
:::

### 混合境界条件が作る半整数の固有値

$(0,L)$ で左 Dirichlet $X(0)=0$、右 Neumann $X'(L)=0$ とする。$-X''=\lambda X$ に $X$ を掛けて積分すると境界項は零で、$\lambda\int X^2=\int(X')^2$ である。従って非零解に対し $\lambda>0$。$X=A\sin(\sqrt\lambda x)$ なので、右端は $\cos(\sqrt\lambda L)=0$ を要求し、$\sqrt\lambda=(j+1/2)\pi/L$（$j\ge0$）である。直交性は異なる二つの ODE に相手の関数を掛け、積分して引くと従う。完全性は左に奇、右に偶に延長し $4L$ 周期とした Fourier 級数から得られる。このとき $2L$ の平行移動で符号が反転するので、周期 $4L$ における奇数の正弦モードだけが残る。

例えば $L=\pi$ で $g=3\sin(x/2)-\sin(3x/2)$ とすれば、熱方程式では $\kappa=4$ のとき $u=3e^{-t}\sin(x/2)-e^{-9t}\sin(3x/2)$。時刻 $\log2$ の右端では二つの正弦の値が $1,-1$ なので $u=3/2+1/512$ である。第二モードも右端の導関数は零であり、値が零かどうかを調べて固有値を選ぶわけではない。

## 相似解と進行波解

:term[相似解]{en="self-similar solution" yomi="そうじかい"} は、尺度を変えた変数で同じ形を保つ解である。拡散では $x\mapsto\lambda x$ と $t\mapsto\lambda^2 t$ が対応するため $\xi=x/\sqrt t$ が自然である。一次元の熱方程式に $u(t,x)=F(\xi)$ を代入すると $\kappa F''+\xi F'/2=0$ となる。段差データ $u(0,x)=0$（$x<0$）、$1$（$x>0$）に対応するには $F(-\infty)=0$、$F(+\infty)=1$ とする。従って $F'=(4\pi\kappa)^{-1/2}e^{-\xi^2/(4\kappa)}$ となり、熱方程式の章で得た熱核表示と一致する。初期質量を規格化した点源なら $t^{-1/2}$ も必要になる。

:::example[非線形拡散の有限伝播]
一次元で $u_t=(u^2)_{xx}$、$u\ge0$ とする。移動する台をもつ $u=t^{-\alpha}F(xt^{-\beta})$ の質量を一定にすると $\alpha=\beta$ であり、時間のべきの一致は $\alpha+1=2\alpha+2\beta$ である。従って $\alpha=\beta=1/3$ となる。ODE は $(F^2)''+\frac13(\xi F)'=0$ である。原点対称に取れば定数は $0$ となり、$F>0$ で $2F'=-\xi/3$ だから $F=(A-\xi^2/12)_+$ である。台の端で $F^2$ とその一階微分が $0$ であり、時刻微分の界面項も $u=0$ によって消える。従って弱形式が全空間で成り立つ。$A=1$ のとき半幅は $\sqrt{12}\,t^{1/3}$ である。
:::

:term[進行波]{en="traveling wave" yomi="しんこうは"} では $u(t,x)=U(x-ct)$ とおく。$c$ は形を変えずに進む速度である。粘性 Burgers 方程式 $u_t+uu_x=\nu u_{xx}$ では積分して $\nu U'=-cU+U^2/2+C$ となる。

:::example[衝撃波を滑らかな層にする]
左 $a$、右 $b$、$a>b$ の両極限を持つ定常形では $c=(a+b)/2$ が必要で、$\nu U'=(U-a)(U-b)/2$ になる。分離積分して

$$
U(\xi)=b+\frac{a-b}{1+\exp((a-b)\xi/(2\nu))}
\label{eq-pde-viscous-burgers-wave}
$$

を得る。$a=3,b=1$ では $c=2$、中心値は $2$、中心勾配は $-1/(2\nu)$ である。従って $\nu$ が小さいほど層は鋭くなり、跳躍の速度は :eqref[eq-pde-rankine-hugoniot] と一致する。
:::

### 一般の多孔質媒質の相似指数

$u_t=\Delta(u^m)$、$m>1$ を $\R^n$ で考え、$u=t^{-\alpha}F(|x|t^{-\beta})$ とする。質量保存には $\alpha=n\beta$、両辺の時間の指数が合うには $\alpha+1=m\alpha+2\beta$ が必要である。従って $\beta=(n(m-1)+2)^{-1}$。$r=|x|t^{-\beta}$ として PDE に代入し、$r^{n-1}$ を掛けると
$$
(r^{n-1}(F^m)')'+\beta(r^n F)'=0
$$
である。原点で $(F^m)^{\prime}(0)=0$ となる滑らかな放射対称プロファイルを選ぶと、積分定数は零になり、$F>0$ では $mF^{m-2}F'=-\beta r$ である。二度目の積分により $F=(A-\beta(m-1)r^2/(2m))_+^{1/(m-1)}$ を得る。台の端で $F$ と $(F^m)'=-\beta rF$ は零なので、領域を分けて部分積分しても分布の界面項は生じない。$t\downarrow0$ で台の半径は $O(t^\beta)$ へ縮まり質量は一定なので初期値は点源になる。

:::example[二次元二乗拡散の規格化]
$n=m=2$ なら $\alpha=1/2$、$\beta=1/4$、$F(r)=(A-r^2/16)_+$ である。質量は $2\pi\int_0^{4\sqrt A}(A-r^2/16)r\dd r=8\pi A^2$ なので、質量を $2\pi$ とすれば $A=1/2$ である。台半径は $2\sqrt2 t^{1/4}$、中心値は $1/(2\sqrt t)$ である。線形熱核では正時刻に全空間へ広がるため、この有限の台は拡散率が $u=0$ で退化することを反映する。
:::

:::example[相似変数に初期幅を残す]
三次元で $\kappa=1$、$g(x)=e^{-|x|^2/12}$ を考える。相似解の $t$ をずらすと、時間零でも有限な値を持つ初期分布になる。

中心からの距離を $r$、空間次元を $n$ とする。Gauss 型の解を保つには $K=C(t)e^{-r^2/(4d(t))}$ と置く。熱方程式に代入すると $r^2$ の係数から $d^{\prime}=\kappa$、定数項から $C^{\prime}/C=-n\kappa/(2d)$ が得られる。初期データの幅と高さに合わせると、次の解になる。

$$
\begin{aligned}
d(t)&=t + 3\\
u(t,r)&=\frac{3 \sqrt{3} e^{- \frac{r^{2}}{4 t + 12}}}{\left(t + 3\right)^{\frac{3}{2}}}\\
u_t&=\frac{\partial}{\partial t} \frac{3 \sqrt{3} e^{- \frac{r^{2}}{4 t + 12}}}{\left(t + 3\right)^{\frac{3}{2}}}=\frac{3 \sqrt{3} \left(r^{2} - 6 t - 18\right) e^{- \frac{r^{2}}{4 t + 12}}}{4 \left(t + 3\right)^{\frac{7}{2}}}\\
u_r&=\frac{\partial}{\partial r} \frac{3 \sqrt{3} e^{- \frac{r^{2}}{4 t + 12}}}{\left(t + 3\right)^{\frac{3}{2}}}=- \frac{3 \sqrt{3} r e^{- \frac{r^{2}}{4 t + 12}}}{2 \left(t + 3\right)^{\frac{5}{2}}}\\
u_{rr}&=\frac{\partial^{2}}{\partial r^{2}} \frac{3 \sqrt{3} e^{- \frac{r^{2}}{4 t + 12}}}{\left(t + 3\right)^{\frac{3}{2}}}=\frac{3 \sqrt{3} \left(r^{2} - 2 t - 6\right) e^{- \frac{r^{2}}{4 t + 12}}}{4 \left(t + 3\right)^{\frac{7}{2}}}\\
\Delta u&=\frac{3 \sqrt{3} \left(r^{2} - 6 t - 18\right) e^{- \frac{r^{2}}{4 t + 12}}}{4 \left(t + 3\right)^{\frac{7}{2}}}
\end{aligned}
$$

放射状の Laplacian は $u_{rr}+(n-1)u_r/r$ である。$r=0$ の値は $u_r/r$ の極限で読む。次に、中心の値と半値となる半径を比べ、指定された時刻と位置に代入する。

$$
\begin{aligned}
u_t-\kappa\Delta u&=\frac{12 \sqrt{3} r^{2} e^{- \frac{r^{2}}{4 t + 12}}}{\left(t + 3\right)^{\frac{3}{2}} \left(4 t + 12\right)^{2}} - \frac{9 \sqrt{3} e^{- \frac{r^{2}}{4 t + 12}}}{2 \left(t + 3\right)^{\frac{5}{2}}} - \frac{3 \sqrt{3} \left(r^{2} - 6 t - 18\right) e^{- \frac{r^{2}}{4 t + 12}}}{4 \left(t + 3\right)^{\frac{7}{2}}}=0\\
u(t,0)&=\frac{3 \sqrt{3}}{\left(t + 3\right)^{\frac{3}{2}}},\qquad r=0\\
r_{1/2}(t)&=2 \sqrt{t + 3} \sqrt{\log{\left(2 \right)}}\\
u(t_*,r_*)&=\frac{3 \sqrt{3}}{8 e},\qquad t=1,r=4\\
u_r(t_*,r_*)&=- \frac{3 \sqrt{3}}{16 e},\qquad t=1,r=4
\end{aligned}
$$

全量を求めるには直交座標で積に戻す。一方向の因子を $g$ とし、その中心を $x_0$ とする。他の方向でも各中心座標に関して同じ Gauss 積分が得られる。積分が分離されるので、全量は一方向の積分の $n$ 乗、正規化した一方向の分散は二次の積分を零次の積分で割った値になる。

$$
\begin{aligned}
g(t,x)&=\frac{\sqrt{3} e^{- \frac{x^{2}}{4 t + 12}}}{\sqrt{t + 3}}\\
M_1&=\int\limits_{-\infty}^{\infty} \frac{\sqrt{3} e^{- \frac{x^{2}}{4 t + 12}}}{\sqrt{t + 3}}\, dx=\frac{\sqrt{3} \sqrt{\pi} \sqrt{4 t + 12}}{\sqrt{t + 3}}\\
M_2&=\int\limits_{-\infty}^{\infty} \frac{\sqrt{3} x^{2} e^{- \frac{x^{2}}{4 t + 12}}}{\sqrt{t + 3}}\, dx=\frac{\sqrt{3} \sqrt{\pi} \left(4 t + 12\right)^{\frac{3}{2}}}{2 \sqrt{t + 3}}\\
\text{全量}&=\frac{3 \sqrt{3} \pi^{\frac{3}{2}} \left(4 t + 12\right)^{\frac{3}{2}}}{\left(t + 3\right)^{\frac{3}{2}}}=24 \sqrt{3} \pi^{\frac{3}{2}}\\
\text{一方向の分散}&=2 t + 6=2 \left(t + 3\right)
\end{aligned}
$$

中心からの距離を $\sqrt{3+t}$ で割った値に固定すると、指数因子は一定である。中心値を掛け合わせて全量の変化を消すには三次元の Jacobian $(3+t)^{3/2}$ を使う。$t$ だけで相似変数を作る点源解と異なり、初期分布の二次モーメントが既に正である。
:::

## Fourier 変換と Laplace 変換による解法

本章は [Fourier 解析](/math/fourier-analysis/fourier-transform/) の規約 $\widehat f(\xi)=\int_{\R^n}f(x)e^{-2\pi i x\cdot\xi}\,\dd x$ に従う。Schwartz 関数では部分積分により $\widehat{\partial_j f}=2\pi i \xi_j\widehat f$ である。従って $u_t-\kappa\Delta u=f$ の空間変換は

$$
\partial_t\widehat u+4\pi^2\kappa|\xi|^2\widehat u=\widehat f
$$

である。積分因子で解いて逆変換すると :eqref[eq-pde-heat-duhamel] を得る。波動方程式では $u_{tt}$ から二階の ODE となり、変位の乗数は $\cos(2\pi c|\xi|t)$、速度の乗数は $\sin(2\pi c|\xi|t)/(2\pi c|\xi|)$ である。$\xi=0$ の除算は極限 $t$ で読む。

時間変換 $\widetilde u(s,x)=\int_0^\infty e^{-st}u(t,x)\,\dd t$ を :term[Laplace 変換]{en="Laplace transform" yomi="らぷらすへんかん"} という。時間について $C^2$ とし、$u,u_t,u_{tt}$ が $Ce^{\gamma t}$ 以下の大きさを持つとする。$\Re s>\gamma$ なら部分積分で $\widetilde{u_t}=s\widetilde u-u(0)$、$\widetilde{u_{tt}}=s^2\widetilde u-su(0)-u_t(0)$ となる。空間微分も変換する場合は、空間導関数について局所的に共通な指数上界を使って積分と交換する。境界に入力があれば、空間 ODE の境界条件にもその変換を課す。

:::example[半直線への階段状入力]
$x>0$、$\kappa=1$、初期値 $0$、境界値 $u(t,0)=1$ とする。変換すると $\widetilde u_{xx}=s\widetilde u$ であり、無限遠で有界なら $\widetilde u=s^{-1}e^{-x\sqrt{s}}$ である。ここで $\Re s>0$ では実部が正の平方根を選ぶ。時間側では $u=\operatorname{erfc}(x/(2\sqrt t))$ が候補になる。ここで $\operatorname{erfc}(z)=2\pi^{-1/2}\int_z^\infty e^{-r^2}\,\dd r$ である。直接微分で $u_t=u_{xx}$、境界値 $1$、$x>0$ の初期値 $0$ を確認できる。従って候補の変換も同じ一意な空間 ODE を満たす。このデータは角 $(0,0)$ で整合しないので、そこでの連続性は要求しない。
:::

:::example[零周波数の意味]
$f=0$ の熱方程式で $g\in L^1$ なら $\widehat u(t,0)=\widehat g(0)$ であり、総熱量が保存される。全空間波では、例えば $g=0$、$h$ が滑らかでコンパクトな台を持つなら、有限伝播により各時刻で $u$ も可積分である。従って $\widehat u(t,0)=t\widehat h(0)$ は総変位が時間に比例することを表す。
:::

### 半空間の Poisson 核を Fourier 乗数で求める

上半空間の調和関数 $u(x',h)$ について、接変数 $x'\in\R^{n-1}$ のみを Fourier 変換する。境界値 $g$ を Schwartz 関数とすれば、各周波数で $\partial_{hh}\widehat u=4\pi^2|\xi|^2\widehat u$、$\widehat u(\xi,0)=\widehat g(\xi)$ となる。高さについて有界な枝を選び、
$$
\widehat u(\xi,h)=e^{-2\pi|\xi|h}\widehat g(\xi)
$$
を得る。$\xi=0$ では ODE の解は高さのアフィン関数であり、有界性により定数 $\widehat g(0)$ となる。逆変換では $\widehat g$ が急減少するため、内部で二階微分を積分に入れられる。また $h\downarrow0$ では優収束により元の境界値に一様収束する。従って Laplace 方程式の章で得た有界な Poisson 積分に一致する。

特に二次元では、乗数の逆変換を正負の周波数に分けると
$$
\int_{\R}e^{-2\pi h|\xi|}e^{2\pi i x\xi}\,\dd\xi
=\frac{1}{2\pi(h-ix)}+\frac{1}{2\pi(h+ix)}
=\frac{h}{\pi(x^2+h^2)}
$$
となる。この例では、境界からの距離 $h$ に応じて高い振動数が強く減衰することが乗数から直接読める。

### 時間 Laplace 変換で得られる境界核

半直線の熱方程式で初期値零、境界データ $b(t)$ とする。変換後の解は $\widetilde b(s)e^{-x\sqrt{s/\kappa}}$ である。階段応答は $\operatorname{erfc}(x/(2\sqrt{\kappa t}))$ なので、その時間微分
$$
P(t,x)=\frac{x}{2\sqrt{\pi\kappa}t^{3/2}}e^{-x^2/(4\kappa t)}
$$
の変換は $e^{-x\sqrt{s/\kappa}}$ である。$x>0$ で階段応答の初期値が $0$ なので、部分積分で端点項は出ない。従って時間畳み込み $u(t,x)=\int_0^t P(t-s,x)b(s)\dd s$ は境界の入力を与える。$P(t,x)\ge0$、$\int_0^\infty P(t,x)\dd t=1$ であり、$x\downarrow0$ で任意の正の時間より後の尾部質量は零へ行く。従って有限時刻で連続な $b$ なら、この畳み込みは $x\downarrow0$ で $b(t)$ に収束する。$x>0$ で微分すると $P_t=\kappa P_{xx}$ であり、$P(0,x)=0$ なので畳み込みも PDE を満たす。

:::example[Fourier 乗数が Gaussian の幅を足す]
一次元で $g(x)=e^{-(x+1)^2/4}$、$\kappa=3$ とする。Fourier 規約は指数が $-2\pi i x\xi$ のものである。

中心からの距離を $r$、空間次元を $n$ とする。Gauss 型の解を保つには $K=C(t)e^{-r^2/(4d(t))}$ と置く。熱方程式に代入すると $r^2$ の係数から $d^{\prime}=\kappa$、定数項から $C^{\prime}/C=-n\kappa/(2d)$ が得られる。初期データの幅と高さに合わせると、次の解になる。

$$
\begin{aligned}
d(t)&=3 t + 1\\
u(t,r)&=\frac{e^{- \frac{r^{2}}{12 t + 4}}}{\sqrt{3 t + 1}}\\
u_t&=\frac{\partial}{\partial t} \frac{e^{- \frac{r^{2}}{12 t + 4}}}{\sqrt{3 t + 1}}=\frac{3 \left(r^{2} - 6 t - 2\right) e^{- \frac{r^{2}}{12 t + 4}}}{4 \left(3 t + 1\right)^{\frac{5}{2}}}\\
u_r&=\frac{\partial}{\partial r} \frac{e^{- \frac{r^{2}}{12 t + 4}}}{\sqrt{3 t + 1}}=- \frac{r e^{- \frac{r^{2}}{12 t + 4}}}{2 \left(3 t + 1\right)^{\frac{3}{2}}}\\
u_{rr}&=\frac{\partial^{2}}{\partial r^{2}} \frac{e^{- \frac{r^{2}}{12 t + 4}}}{\sqrt{3 t + 1}}=\frac{\left(r^{2} - 6 t - 2\right) e^{- \frac{r^{2}}{12 t + 4}}}{4 \left(3 t + 1\right)^{\frac{5}{2}}}\\
\Delta u&=\frac{\left(r^{2} - 6 t - 2\right) e^{- \frac{r^{2}}{12 t + 4}}}{4 \left(3 t + 1\right)^{\frac{5}{2}}}
\end{aligned}
$$

放射状の Laplacian は $u_{rr}+(n-1)u_r/r$ である。$r=0$ の値は $u_r/r$ の極限で読む。次に、中心の値と半値となる半径を比べ、指定された時刻と位置に代入する。

$$
\begin{aligned}
u_t-\kappa\Delta u&=\frac{12 r^{2} e^{- \frac{r^{2}}{12 t + 4}}}{\sqrt{3 t + 1} \left(12 t + 4\right)^{2}} - \frac{3 e^{- \frac{r^{2}}{12 t + 4}}}{2 \left(3 t + 1\right)^{\frac{3}{2}}} - \frac{3 \left(r^{2} - 6 t - 2\right) e^{- \frac{r^{2}}{12 t + 4}}}{4 \left(3 t + 1\right)^{\frac{5}{2}}}=0\\
u(t,0)&=\frac{1}{\sqrt{3 t + 1}},\qquad r=0\\
r_{1/2}(t)&=2 \sqrt{3 t + 1} \sqrt{\log{\left(2 \right)}}\\
u(t_*,r_*)&=\frac{1}{2 e^{\frac{1}{4}}},\qquad t=1,r=2\\
u_r(t_*,r_*)&=- \frac{1}{8 e^{\frac{1}{4}}},\qquad t=1,r=2
\end{aligned}
$$

全量を求めるには直交座標で積に戻す。一方向の因子を $g$ とし、その中心を $x_0$ とする。他の方向でも各中心座標に関して同じ Gauss 積分が得られる。積分が分離されるので、全量は一方向の積分の $n$ 乗、正規化した一方向の分散は二次の積分を零次の積分で割った値になる。

$$
\begin{aligned}
g(t,x)&=\frac{e^{- \frac{\left(x + 1\right)^{2}}{12 t + 4}}}{\sqrt{3 t + 1}}\\
M_1&=\int\limits_{-\infty}^{\infty} \frac{e^{- \frac{\left(x + 1\right)^{2}}{12 t + 4}}}{\sqrt{3 t + 1}}\, dx=\frac{\sqrt{\pi} \sqrt{12 t + 4}}{\sqrt{3 t + 1}}\\
M_2&=\int\limits_{-\infty}^{\infty} \frac{\left(x + 1\right)^{2} e^{- \frac{\left(x + 1\right)^{2}}{12 t + 4}}}{\sqrt{3 t + 1}}\, dx=\frac{\sqrt{\pi} \left(12 t + 4\right)^{\frac{3}{2}}}{2 \sqrt{3 t + 1}}\\
\text{全量}&=\frac{\sqrt{\pi} \sqrt{12 t + 4}}{\sqrt{3 t + 1}}=2 \sqrt{\pi}\\
\text{一方向の分散}&=6 t + 2=2 \left(3 t + 1\right)
\end{aligned}
$$

初期値の Fourier 変換は $\sqrt{4\pi}e^{-4\pi^2\xi^2+2\pi i\xi}$ である。熱の乗数 $e^{-12\pi^2t\xi^2}$ を掛けると二次指数に幅 $1+3t$ が現れ、逆変換で上の式を得る。平行移動は位相だけを変え、減衰する速さは平行移動に依存しない。
:::

:::example[時間 Laplace 変換を簡単な固有成分で逆にする]
$(0,\pi)$ の熱方程式で第2正弦モードの初期係数 $2$、源 $3e^{-t}$ を使う。時間積分の表示と部分分数分解を比べる。

まず $u(t,x)=a(t)\sin(mx)$ を代入する。空間二階の微分は $-m^2$ 倍であり、時間微分は $a^{\prime}$ だけに作用する。従って時間依存係数を含む場合でも、他の空間モードを結合させなければ一次の方程式に帰着する。ここで $\mu(t)$ を空間作用素のこのモードでの係数とする。

$$
\begin{aligned}
m&=2\\
\mu(t)&=4\\
f_m&=3 e^{- t}\\
M(t) &= \int_0^t \mu(s)\,\mathrm{d}s=\int\limits_{0}^{t} 4\, dz=4 t\\
I(t)&=e^{M(t)}=e^{4 t}
\end{aligned}
$$

積分因子を掛けると $(Ia)^{\prime}=If_m$ である。零時刻から積分し、初期係数を加えてから因子で割る。外力を単に時間積分するのと違い、積分の中で指数の重みが掛かる。

$$
\begin{aligned}
If_m&=3 e^{3 t}\\
\int_0^t I(s)f_m(s)\,\mathrm{d}s&=\int\limits_{0}^{t} 3 e^{3 z}\, dz=e^{3 t} - 1\\
a(t)&=\left(e^{t} + 1\right) \left(e^{2 t} - e^{t} + 1\right) e^{- 4 t}\\
a^{\prime}&=\frac{d}{d t} \left(e^{t} + 1\right) \left(e^{2 t} - e^{t} + 1\right) e^{- 4 t}=- \left(e^{3 t} + 4\right) e^{- 4 t}\\
a^{\prime}+\mu a&=\left(e^{t} + 1\right) \left(2 e^{2 t} - e^{t}\right) e^{- 4 t} + \left(e^{2 t} - e^{t} + 1\right) e^{- 3 t}=3 e^{- t}\\
a(0)&=2,\qquad t=0
\end{aligned}
$$

正弦の二乗積分は $\pi/2$ である。解の二乗ノルムの半分を時間微分し、散逸および入力から直接求めた変化率と比べる。係数が時間に依存しても、ここで微分する二乗ノルム自体には時間依存の作用素が入らない。

$$
\begin{aligned}
u(t,x)&=\left(e^{t} + 1\right) \left(e^{2 t} - e^{t} + 1\right) e^{- 4 t} \sin{\left(2 x \right)}\\
E(t)&=\frac{\pi \left(e^{t} + 1\right)^{2} \left(e^{2 t} - e^{t} + 1\right)^{2} e^{- 8 t}}{4}\\
E^{\prime}&=\frac{d}{d t} \frac{\pi \left(e^{t} + 1\right)^{2} \left(e^{2 t} - e^{t} + 1\right)^{2} e^{- 8 t}}{4}=- \frac{\pi \left(e^{t} + 1\right) \left(e^{3 t} + 4\right) \left(e^{2 t} - e^{t} + 1\right) e^{- 8 t}}{2}\\
\frac{\pi}2 a(f_m-\mu a)&=- \frac{\pi \left(e^{t} + 1\right) \left(e^{3 t} + 4\right) \left(e^{2 t} - e^{t} + 1\right) e^{- 8 t}}{2}\\
u(t_*,x_*)&=\frac{9}{16},\qquad t=\log{\left(2 \right)},x=\frac{\pi}{4}\\
E(t_*)&=\frac{81 \pi}{1024},\qquad t=\log{\left(2 \right)}
\end{aligned}
$$

時間変換で $(s+4)\widetilde a=2+3/(s+1)$ となる。部分分数は $3/((s+4)(s+1))=1/(s+1)-1/(s+4)$ であり、変換を戻すと $e^{-t}+e^{-4t}$ になる。初期項を変換時に落とすと、$t=0$ で零となる外力の寄与だけが残ってしまう。
:::

:::example[時間の分母と指数の逆変換]
零境界を持つ拡散の第2モードに、源 $f_m=4+3e^{-t}$、初期値 $2$ を加える。時間の減衰率は $4$ とする。

時間変換の変数を $s$ とし、$\mathcal L[f](s)=\int_0^\infty e^{-st}f(t)\dd t$ とする。大きい $s$ の範囲で積分が収束する。時間微分を部分積分すると、$t=0$ の値が分子の項になる。まず外力の変換、微分から来る分母、初期データから来る分子を分ける。

$$
\begin{aligned}
f(t)&=4 + 3 e^{- t}\\
\widetilde f(s)&=\frac{7 s + 4}{s \left(s + 1\right)}\\
\text{微分の多項式}&=s + 4\\
\text{初期項}&=2\\
\widetilde a(s)&=\frac{2 s + 1}{s^{2} + s}
\end{aligned}
$$

次に分母を因数分解する。一つの極が重なるときは単純な指数の係数に加え、時間のべきが掛かる項を求める。例えば $(s+\lambda)^{-2}$ は $te^{-\lambda t}$、分母の二次式 $(s+\lambda)^2+\omega^2$ は減衰する正弦または余弦に対応する。

$$
\begin{aligned}
\widetilde a(s)&=\frac{1}{s + 1} + \frac{1}{s}\\
a(t)&=\left(e^{t} + 1\right) e^{- t}\\
a^{\prime}&=\frac{d}{d t} \left(e^{t} + 1\right) e^{- t}=- e^{- t}
\end{aligned}
$$

最後に時刻を戻して初期条件を調べ、次に得た時間関数を元の左辺に代入する。変換領域の式を満たしても、逆変換で定数や初期項を落とすとこの二つの確認のどちらかに失敗する。

$$
\begin{aligned}
a(0)&=2,\qquad t=0\\
\text{時間方程式の左辺}&=3 \left(e^{t} + 1\right) e^{- t} + 1=4 + 3 e^{- t}\\
a(t_*)&=\frac{3}{2},\qquad t=\log{\left(2 \right)}\\
a^{\prime}(t_*)&=- \frac{1}{2},\qquad t=\log{\left(2 \right)}
\end{aligned}
$$

定数の源は係数を $1$ に引き寄せ、指数の源は $e^{-t}$ の寄与を残す。この初期値の選択では $e^{-4t}$ の係数が相殺する。減衰率 $4$ が式の分母に現れても、初期項と強制項の和の分子によってその極が消える場合がある。
:::

## Cole–Hopf 変換と Legendre 変換

:::proposition[Cole–Hopf 変換]{#prop-pde-cole-hopf}
$w>0$ が $w_t=\nu w_{xx}$ を満たせば $U=-2\nu\log w$ は $U_t+\tfrac12U_x^2=\nu U_{xx}$ を満たし、$u=U_x$ は粘性 Burgers 方程式を満たす。逆に滑らかな $U$ がこの式を満たせば $e^{-U/(2\nu)}$ は熱方程式を満たす。
:::

:::proof
$U_t=-2\nu w_t/w$、$U_{xx}=-2\nu(w_{xx}/w-(w_x/w)^2)$、$\tfrac12U_x^2=2\nu^2(w_x/w)^2$ を代入すると

$$
U_t+\tfrac12U_x^2-\nu U_{xx}=-2\nu(w_t-\nu w_{xx})/w
\label{eq-pde-cole-hopf}
$$

である。$x$ で微分すれば $u$ の式を得る。$u$ の式から $U$ を逆に作るとき、空間積分の定数は時刻に依存する。$U_t+\frac12U_x^2-\nu U_{xx}$ は空間一定なので、この余分な時間関数を $U$ の加法項で消せる。
:::

:::example[正の熱解から波の層を作る]
$w(t,x)=1+e^{x+t}$、$\nu=1$ とすれば熱方程式を満たす。従って $u=-2e^{x+t}/(1+e^{x+t})$ は粘性 Burgers 解である。左状態は $0$、右状態は $-2$、速度は $-1$ であり、:eqref[eq-pde-viscous-burgers-wave] と一致する。
:::

凸関数 $U(t,\cdot)$ の :term[Legendre 変換]{en="Legendre transform" yomi="るじゃんどるへんかん"} は $V(t,p)=\sup_x(p\cdot x-U(t,x))$ である。$D_x^2U>0$ で $U$ が超線形なら最適点 $x$ は一つで $D_xU=p$ になる。時間微分では $p\cdot x_t-D_xU\cdot x_t=0$ なので $V_t=-U_t$ である。従って $U_t+H(D_xU)=0$ は $V_t=H(p)$ に変わる。$V(t)=V(0)+tH$ の凸共役を取ると $U$ を復元できる。凸性が失われる場合は逐点での勾配の逆写像が使えないため、[Hopf–Lax の公式](/math/pde/nonlinear-first-order/) のような最適化の形が必要になる。

:::example[二次関数の変換]
$U(0,x)=x^2$、$H(p)=p^2/2$ なら $V(0,p)=p^2/4$ である。$V(t,p)=(1+2t)p^2/4$ を $p$ について変換すると $U=x^2/(1+2t)$ を得る。勾配の逆写像が $x=(1+2t)p/2$ であり、特性から求めた式と同じになる。
:::

### Cole–Hopf に残る時間だけの規格化

$u=U_x$ とし、粘性 Burgers 方程式を $x$ について積分すると $U_t+U_x^2/2-\nu U_{xx}=C(t)$ となる。$U$ から $\int_0^t C(s)\dd s$ を引いても $u$ は変わらない。$w=e^{-U/(2\nu)}$ とすると
$$
w_t-\nu w_{xx}=-\frac{w}{2\nu}\left(U_t+\frac{U_x^2}{2}-\nu U_{xx}\right)
$$
であり、時間だけの規格化を済ませれば熱方程式になる。$w_0>0$ の熱核畳み込みが有限なら $w>0$ なので $u=-2\nu w_x/w$ で逆変換できる。複数の $w$ の正の和は熱方程式を満たすが、対応する $u$ は単純な和ではなく、その位置での $w$ の比を重みにした平均になる。

:::example[二つの指数モードが作る移動層]
$w=e^{k_1x+\nu k_1^2t}+e^{k_2x+\nu k_2^2t}$、$k_1=-3/(2\nu)$、$k_2=-1/(2\nu)$ とすると $u=(3w_1+w_2)/(w_1+w_2)$ である。$w_2/w_1=e^{(x-2t)/\nu}$ なので中心は $x=2t$、値は $2$。左側では $w_1$ が優勢になって $u\to3$、右側では $u\to1$ である。速度 $2$ は両熱モードの指数が等しくなる直線の傾きとして得られる。
:::

:::example[Legendre 変換で勾配を反転する]
$F(p)=3p^2/2-2p+1$ の凸共役を求めるには $F(a)-av$ を最小化して負号を付ける。$v=4$ で二つの試行値を比べる。

まず候補の選択を一変数 $a$ の問題にする。停留条件は $J_a=0$ であり、境界や無限遠も許される場合は停留点と分けて扱う。次の式では一階微分を取って零になる点を求め、二階微分を残しておく。

$$
\begin{aligned}
J(a)&=\frac{3 a^{2}}{2} - a v - 2 a + 1\\
J_a&=\frac{\partial}{\partial a} \left(\frac{3 a^{2}}{2} - a v - 2 a + 1\right)=3 a - v - 2\\
J_{aa}&=\frac{\partial^{2}}{\partial a^{2}} \left(\frac{3 a^{2}}{2} - a v - 2 a + 1\right)=3\\
a_*&=\frac{v}{3} + \frac{2}{3}\\
J_a(a_*)&=0,\qquad a=\frac{v}{3} + \frac{2}{3}
\end{aligned}
$$

停留値と、その点から $h$ だけ動かした差を計算する。一次項が相殺しているため、展開の先頭から局所の振る舞いが分かる。一方、停留点が全体で最小か最大かを結論するには、許容範囲全体でこの差の符号を調べる必要がある。

$$
\begin{aligned}
J(a_*)&=- v \left(\frac{v}{3} + \frac{2}{3}\right) - \frac{2 v}{3} + \frac{3 \left(\frac{v}{3} + \frac{2}{3}\right)^{2}}{2} - \frac{1}{3}=- \frac{v^{2} + 4 v - 2}{6},\qquad a=\frac{v}{3} + \frac{2}{3}\\
J_{aa}(a_*)&=3\\
J(a_*+h)-J(a_*)&=- 2 h - v \left(h + \frac{v}{3} + \frac{2}{3}\right) - \frac{2 v}{3} + \frac{3 \left(h + \frac{v}{3} + \frac{2}{3}\right)^{2}}{2} + \frac{v^{2} + 4 v - 2}{6} - \frac{1}{3}=\frac{3 h^{2}}{2}\\
\text{二次係数}&=\frac{3}{2}
\end{aligned}
$$

指定した値で目的関数を比較する。$a$ 自体の値と、そこでの目的関数の値は別の量である。解の族を作る場合は、選んだ $a$ を元の式に代入して得られる $J(a_*;x)$ が包絡線または値関数になる。

$$
\begin{aligned}
a_*\text{ の具体値}&=2\\
J_*\text{ の具体値}&=-5\\
J&=1,\qquad a=0\\
J&=-5,\qquad a=2\\
\text{二つの値の差}&=-6\\
\text{停留点まわりの差}&=\frac{3 h^{2}}{2}
\end{aligned}
$$

凸共役は $F^*(v)=(v+2)^2/6-1$ である。$(F^*)^{\prime}(v)=(v+2)/3$ は $F^{\prime}(p)=3p-2$ の逆関数であり、もう一度凸共役を取ると元へ戻る。最後の等号では係数 $3$ を逆数にし、一次項の平行移動も保持する必要がある。
:::

### Legendre 座標で二階の係数を反転する

滑らかな厳密凸解 $u(t,x)$ の勾配 $x\mapsto D_xu$ が局所的に可逆であり、$D_x^2u>0$ とする。$p=D_xu(t,x)$ に対し $v(t,p)=p\cdot x-u(t,x)$ と置く。$p$ を固定したとき $x_t$ の項は $p-D_xu=0$ で打ち消され、$v_t=-u_t$ になる。$p$ で微分しても同じ打消しで $D_pv=x$、さらに $D_p^2v=(D_x^2u)^{-1}$ となる。従って $u_t+H(D_xu)=0$ は $v_t=H(p)$ に変わり、$v(t,p)=g^*(p)+tH(p)$ である。一方、二階の PDE は Hessian の逆を含むことがある。従っていつも線形化されるのではなく、元の方程式で勾配と Hessian がどこに入るかに依存する。

:::example[初期値の二次係数を反転して時間発展させる]
$g(x)=x^2+2x+3$、$H(p)=p^2$ とする。$g'(x)=2x+2$ なので $x=(p-2)/2$、そのため $g^*(p)=(p-2)^2/4-3$。$v=g^*+tp^2$ は $t\ge0$ で二階微分 $1/2+2t>0$ を持つ。逆変換の停留条件は $x=(p-2)/2+2tp$、よって $p=2(x+1)/(1+4t)$。計算すると
$$
u(t,x)=2+\frac{(x+1)^2}{1+4t}.
$$
$t=0$ で $g$ に一致し、$u_t=-4(x+1)^2/(1+4t)^2$ と $u_x^2$ は相殺する。例えば $t=1/4,x=1$ で $u=4$ になる。時刻を負に延ばすと $t=-1/4$ で Legendre 座標の凸性と逆変換が失われる。
:::

### Cole–Hopf の初期値に必要な大きさ

有界な初期速度 $g$ から原始関数 $U_0(x)=\int_0^xg$ を作ると $|U_0(x)|\le\norm g_\infty|x|$ である。従って $w_0=e^{-U_0/(2\nu)}$ は指数増大し得るが、任意の $t>0$ で Gauss 核との畳み込みは収束する。実際、線形の増大 $C|y|$ は $|x-y|^2/(8\nu t)$ と定数で抑えられ、残りの Gauss 因子は可積分である。$w_0$ を定数倍しても $-2\nu(\log w)_x$ は同じなので、原始関数の積分定数に依存しない。例えば $g=b$ が定数なら $w_0=e^{-bx/(2\nu)}$、熱方程式の解は $w=e^{-bx/(2\nu)+b^2t/(4\nu)}$ であり、対数微分で元の定数 $b$ に戻る。初期速度が有界なのに $w_0$ が有界でない場合もあり、熱核表示で許される増大と初期速度の増大は区別する。

### 小さい粘性の極限と最小作用

$U_0$ を大域 Lipschitz な初期値とする。Cole–Hopf の熱核表示は、$t>0$ に対し
$$
U^\nu(t,x)=-2\nu\log\left(\frac{1}{\sqrt{4\pi\nu t}}
\int_\R e^{-J(y)/(2\nu)}\,\dd y\right),\qquad
J(y)=U_0(y)+\frac{(x-y)^2}{2t}
$$
である。$J$ は無限遠で二次以上の速さで増大するので最小値 $J_*$ を取る。任意の $\delta>0$ について、最小点の近くに $J\le J_*+\delta$ となる正の長さの区間がある。その区間で積分を下から抑えれば $\limsup_{\nu\downarrow0}U^\nu(t,x)\le J_*+\delta$ になる。

逆向きには、十分大きい $R>0$ で $J-J_*\ge |y|^2/(8t)$（$|y|>R$）とできる。$e^{J_*/(2\nu)}$ を掛けた積分は $|y|\le R$ で $2R$ 以下、外側で全直線上の Gaussian 積分以下である。また熱核の係数の対数は $\nu$ を掛ければ零に収束する。よって $\liminf_{\nu\downarrow0}U^\nu(t,x)\ge J_*$ であり、$\delta\downarrow0$ で極限が $J_*$ になる。この最小作用は :eqref[eq-pde-hopf-lax] の $H(p)=p^2/2$ の場合である。例えば $U_0=-|x|$ では $x=0$ に二つの最小点があり、極限には角が残る。

## 漸近解析と幾何光学近似

大きな振動数 $1/\eps$ に対し、滑らかな位相 $\phi$ と振幅 $a$ を用いて $u^\eps=a e^{i\phi/\eps}$ と置く。これを :term[幾何光学近似]{en="geometric optics approximation" yomi="きかこうがくきんじ"} の出発点とする。例えば $\eps^2\Delta u+n(x)^2u=0$ に代入して指数因子を除くと

$$
(n^2-|D\phi|^2)a+i\eps(2D\phi\cdot Da+(\Delta\phi)a)+\eps^2\Delta a
\label{eq-pde-wkb-residual}
$$

となる。各次数を順に消すと $|D\phi|^2=n^2$ と $2D\phi\cdot Da+(\Delta\phi)a=0$ を得る。前者を :term[アイコナール方程式]{en="eikonal equation" yomi="あいこなーるほうていしき"} と呼ぶ。実振幅に対し後者に $a$ を掛けると $\divg(a^2D\phi)=0$ になり、光線束に沿った振幅の変化を表す。

:::example[平面波と球面波]
$n=2$ が一定の屈折率を表すとする。$\phi=2x_1$ と一定振幅なら残差も $0$ である。三次元で $\phi=2r$ なら $\Delta\phi=4/r$ であり、振幅が動径だけに依存すると輸送式から $a=C/r$ になる。$r>0$ で $\Delta(1/r)=0$ なので $e^{2ir/\eps}/r$ は原点の外で厳密解である。原点では位相と振幅の仮定が崩れる。
:::

:::example[一変数の屈折率]
$n(x)=1+x$ を $x>-1$ で考える。正の枝を使うと $\phi=x+x^2/2$、$a=(1+x)^{-1/2}$ である。規格化した残差は $\eps^2 a''=3\eps^2(1+x)^{-5/2}/4$ である。従って $x=-1$ から離れた範囲では残差が二次で小さいが、その境界まで一様な評価はこの計算から得られない。残差の小ささを真の解との差に直すには、境界条件を持つ作用素の逆の評価も必要である。
:::

### 波動方程式の位相と時間を含む輸送

定数 $c$ の波動作用素 $P=\partial_{tt}-c^2\Delta$ で $u^\eps=a e^{i\phi/\eps}$ と置く。指数を除いた $\eps^2 Pu^\eps$ の零次項は $(-\phi_t^2+c^2|D\phi|^2)a$、一次項は

$$
i\eps\big(2\phi_t a_t-2c^2D\phi\cdot Da+(P\phi)a\big)
$$

であり、二次項は $\eps^2Pa$ である。時間と空間を別々に二度微分するとこの二次式と交差項が得られる。$\phi_t\ne0$ の枝では $dX/dt=-c^2D\phi/\phi_t$ の軌道に沿って、振幅は $d a/dt=-(P\phi)a/(2\phi_t)$ の ODE に従う。$a$ を掛けると

$$
\partial_t(a^2\phi_t)-c^2\divg(a^2D\phi)=0
\label{eq-pde-wave-ray-conservation}
$$

である。位相の速度と輸送の向きは一次の符号を含めて決める。

例えば $\phi=x_1-ct$ では $\phi_t=-c$、$D\phi=e_1$、$P\phi=0$ なので $a_t+c a_{x_1}=0$ となる。$a=A(x_1-ct)$ とすれば $Pa=0$ でもあり、合成波は厳密解である。別の例として三次元の $\phi=r-ct$ では $P\phi=-2c^2/r$ なので $a_t+c a_r+(c/r)a=0$ である。従って $a=A(r-ct)/r$ となる。$r u$ が一次元波動方程式を満たすことから、原点を除けば近似が厳密である。$r=0$ では球座標の速度場が定義できず、この形から正則な解まで主張できない。

残差の評価を実際の誤差に移す場合、例えば固定した有界領域で $u^\eps$ が真の解と同じ Dirichlet データ、同じ初期変位と初速度を持つとする。$e=u-u^\eps$ の外力は $-P u^\eps$、初期エネルギーは零である。実部と虚部のエネルギーを合わせ、境界で $e_t=0$ を使うと $\sqrt{2E_e(t)}\le\int_0^t\norm{P u^\eps(s)}_2\,\dd s$ となる。残差を $\eps^2 P$ について計算した場合はここで $\eps^{-2}$ を戻す必要がある。零次と一次だけを消した場合、エネルギーの絶対誤差評価は一般に $O(1)$ である。そこで振幅に $\eps a_1$ を加え、一次の輸送作用素を $L_\phi=2\phi_t\partial_t-2c^2D\phi\cdot D+P\phi$ とおく。$iL_\phi a_1=-Pa$ と選ぶと $\eps^2 P$ の残差から二次項も消え、同じ時間区間で $P$ の残差を一階小さくできる。この評価には、補正後の近似が初期条件と境界条件を満たすことも必要である。

### 行路の広がりを振幅で補う

定数波速 $c$ に対し $u=a(r)e^{i(r-ct)/\eps}$ とすると、$r>0$ で $\Delta u=u_{rr}+(n-1)u_r/r$ である。位相の項は $c^2\phi_r^2=\phi_t^2$ で消え、一階の項は $2a'+(n-1)a/r=0$ を要求する。従って $a(r)=C r^{-(n-1)/2}$ であり、$r^{n-1}|a|^2$ が一定となる。これは球面の断面積の増加に応じた振幅である。残差は $-c^2(a''+(n-1)a'/r)e^{i\phi/\eps}$ となる。$a=r^{-p}$、$p=(n-1)/2$ を代入すれば $a''+(n-1)a'/r=(n-1)(3-n)r^{-p-2}/4$。三次元では零なのでこの形は原点を除いて厳密解である。二次元では $r^{-5/2}/4$ が残り、原点に近づくと評価は悪化する。従って幾何光学近似の残差は波長だけでなく振幅が変化する距離にも依存する。

### Laplace の方法で指数積分の主項を出す

連続な $a$ が有界区間に台を持ち、滑らかな実関数 $\phi$ がその区間で唯一の最小点 $x_0$ を持つとする。$x_0$ は内点、$\phi''(x_0)=d>0$ とする。$I(\eps)=\int a(x)e^{-\phi(x)/\eps}\dd x$ では、$x_0$ の固定した小近傍を除く部分は $|I|\le C e^{-(\phi(x_0)+\delta)/\eps}$ である。これはコンパクト性から $\phi-\phi(x_0)$ が正の下限 $\delta$ を持つためである。小近傍で $x=x_0+\sqrt\eps y$ と置けば
$$
\frac{e^{\phi(x_0)/\eps}I(\eps)}{\sqrt\eps}
\longrightarrow a(x_0)\int_{\R}e^{-dy^2/2}\dd y
=a(x_0)\sqrt{2\pi/d}
$$
となる。実際 Taylor 公式から固定した $y$ で指数は $-dy^2/2$ へ収束し、近傍を十分小さく取れば $\phi(x)-\phi(x_0)\ge d(x-x_0)^2/4$ なので $C e^{-dy^2/4}$ が可積分な支配関数になる。従って主要項の幅は $\sqrt\eps$ である。最小点が複数なら近傍ごとの寄与を加え、二階微分が零ならこの尺度は変わる。

:::example[正の費用の比から選ばれる状態]
$I(\eps)=\int_{-1}^1 (1+x^2)e^{-2x^2/\eps}\dd x$ では $x_0=0,d=4,a(0)=1$ なので $I/\sqrt\eps\to\sqrt{\pi/2}$ である。置換 $x=\sqrt\eps y$ を直接使うと、奇数項はなく、次の寄与は $\eps^{3/2}\int_{\R}y^2e^{-2y^2}\dd y=\eps^{3/2}\sqrt{\pi/2}/4$ になる。有限区間を全直線へ広げた差は、$|y|\ge\eps^{-1/2}$ における Gauss の尾なのでどのべきよりも小さい。例えば固定点 $x=1/2$ の被積分関数は $5e^{-1/(2\eps)}/4$ であり、中心近傍の寄与と比べて指数的に小さい。
:::

### 停留位相法

$a\in C_c^\infty(\R)$ とし、実の滑らかな $\phi$ は $a$ の台上でただ一つの停留点 $x_0$ を持ち、$d=\phi''(x_0)\ne0$ とする。:term[停留位相法]{en="stationary phase method" yomi="ていりゅういそうほう"} では、指数が振動する積分の主要項を
$$
\int_\R a(x)e^{i\phi(x)/\eps}\dd x
=e^{i\phi(x_0)/\eps+i\operatorname{sgn}(d)\pi/4}
a(x_0)\sqrt{\frac{2\pi\eps}{|d|}}+O(\eps)
$$
と評価する。

停留点の近傍に滑らかな切断を置く。その外では $\phi'$ が零にならず、$e^{i\phi/\eps}=(\eps/(i\phi'))\partial_x e^{i\phi/\eps}$ を使って一回部分積分すれば $O(\eps)$ である。近傍では
$$
y=(x-x_0)\left(\frac{2(\phi(x)-\phi(x_0))}{d(x-x_0)^2}\right)^{1/2}
$$
とおく。根号内は Taylor の積分余項により滑らかに $1$ へ延びる。したがって $y'(x_0)=1$ であり、この座標では $\phi=\phi(x_0)+dy^2/2$ である。変換後の振幅を $b$ とし、原点近くで $1$ の切断 $\chi$ を使って $b=b(0)\chi+yb_1$ と書く。$b_1$ は滑らかで台がコンパクトなので、$yb_1$ の積分は $\partial_y e^{idy^2/(2\eps)}=(idy/\eps)e^{idy^2/(2\eps)}$ で部分積分して $O(\eps)$ となる。

残った項で $\chi$ を $1$ に替えた差も、原点から離れた所で同じ部分積分を使って $O(\eps)$ である。このとき無限遠で $1/y$ が消え、$( (1-\chi)/y)'$ が可積分である。全直線の振動積分は、$e^{-\delta y^2}$ を掛けて Gaussian 積分 $\sqrt{\pi/(\delta-id/(2\eps))}$ を計算し、$\delta\downarrow0$ として求められる。減衰付きと減衰なしの積分が同じ極限を持つのは、$|y|>R$ の部分を分母 $2(\delta-id/(2\eps))y$ で部分積分すれば、尾が $\delta\ge0$ に一様に $O(\eps/R)$ になるためである。$b(0)=a(x_0)$ を戻すと上の式を得る。例えば $\phi(x)=3x^2$、$a(0)=1$ なら主要項は $e^{i\pi/4}\sqrt{\pi\eps/3}$ である。実指数の積分では Gaussian の減衰が寄与を局所化し、振動積分では停留点の外の相殺が同じ役割を持つ。

### 位相と振幅の誤差をそれぞれ計算する

一次元の波動作用素 $P=\partial_t^2-c^2\partial_x^2$ で $c$ を一定とし、$S=x-ct$、$u=a(t,x)e^{iS/\eps}$ とする。積の微分から
$$
e^{-iS/\eps}P u
=P a-\frac{2ic}{\eps}(a_t+ca_x)
$$
となる。位相 $S$ 自身の二階微分は零であり、二次の位相項も $S_t^2-c^2S_x^2=0$ で消える。振幅を $a=A(x-ct)$ と取れば輸送式も $Pa$ も零になる。従ってこの一方向の平面波は任意の滑らかな $A$ で厳密解である。逆に $a=x^2$ と固定してしまうと、$a_t+ca_x=2cx$ なので残差は $-2c^2-4ic^2x/\eps$ になる。位相だけを選んでも残差が小さくならないことが、この場合には微分の次数で確かめられる。

:::example[放射状の厳密解と原点の条件]
定常方程式 $\eps^2(u_{rr}+2u_r/r)+u=0$（$r>0$）で位相を $S=r$ とする。$u=e^{ir/\eps}a(r)$ を代入すると残差は $e^{ir/\eps}[\eps^2(a''+2a'/r)+2i\eps(a'+a/r)]$ である。例えば $a=2/r$ とすれば $a'+a/r=0$、$a''+2a'/r=0$ で完全に消える。内側の原点で有界にしたければ、正負の位相の解を組み合わせる。$u=\sin(r/\eps)/r$ は $r\downarrow0$ で $1/\eps$ へ収束する。出て行く波だけの振幅 $1/r$ は原点までの正則性を満たさないため、境界条件によって許される位相の組合せが変わる。
:::

## 均質化法（概観）

:::theorem[一次元周期媒質の均質化]{#thm-pde-homogenization-one-d}
連続な $1$ 周期関数 $a$ が $0<m\le a\le M$ を満たすとする。$F(x)=\int_0^x f$、$f\in C([0,1])$ とし、$-(a(x/\eps)u_\eps')'=f$、$u_\eps(0)=u_\eps(1)=0$ を考える。$u_\eps$ は一様に $u_0$ に収束し、$-\bar a\,u_0''=f$、ここで

$$
\bar a=\left(\int_0^1\frac{\dd y}{a(y)}\right)^{-1}.
\label{eq-pde-homogenized-coefficient}
$$

:::

:::proof
一度積分すると $a(x/\eps)u_\eps'=C_\eps-F(x)$ であり、端点の差が $0$ なので $C_\eps=\int_0^1F/a(x/\eps)\big/\int_0^1 1/a(x/\eps)$ となる。$b=1/a$、$\bar b=\int_0^1b$ とし、$B(y)=\int_0^y(b-\bar b)$ と置くと $B$ は有界かつ周期的である。任意の $C^1$ の $\psi$ に対して

$$
\int_0^x (b(s/\eps)-\bar b)\psi(s)\,\dd s
=\eps[B(s/\eps)\psi(s)]_0^x-\eps\int_0^x B(s/\eps)\psi'(s)\,\dd s.
$$

右辺は $x$ について一様に $O(\eps)$ である。$\psi=1,F$ で用いると $C_\eps\to\int_0^1F$ となる。さらに $u_\eps(x)=\int_0^x b(s/\eps)(C_\eps-F(s))\,\dd s$ に入れれば、一様極限の導関数が $\bar b(\int_0^1F-F(x))$ であると分かる。
:::

### 周期セル問題と有効行列

高次元でもある一つの平均を係数に直接取ればよいわけではない。一定の勾配 $e_i$ を与えたとき、セル内で勾配を微調整して流束の発散を消す必要がある。ここでセルとは向かい合う面を同一視した $Y=(0,1)^n$ である。

:::theorem[周期な対称係数の均質化]{#thm-pde-periodic-homogenization}
$A(y)$ は有界可測、$Y$ 周期、対称一様楕円とする。有界領域 $\Omega$ と $f\in L^2(\Omega)$ について $u_\eps\in H^1_0(\Omega)$ が $-\divg(A(x/\eps)Du_\eps)=f$ を満たすとする。すると $u_\eps\rightharpoonup u$ が $H_0^1$ で、$u_\eps\to u$ が $L^2$ で成り立ち、$-\divg(\bar A Du)=f$ である。有効行列は

$$
\bar A e_i=\int_Y A(y)(e_i+D\chi_i(y))\,\dd y,\qquad
\divg_y(A(y)(e_i+D\chi_i))=0,\quad \int_Y\chi_i=0
\label{eq-pde-periodic-cell}
$$

で与えられる。
:::

:::proof[セル問題と周期関数の平均]
平均零の周期 $H^1$ 空間で $\int_Y A D\chi_i\cdot Dv=-\int_Y A e_i\cdot Dv$ を解く。周期 Fourier 級数では零モードがなく、$\int_Y|v|^2\le (2\pi)^{-2}\int_Y|Dv|^2$ なので Lax–Milgram の定理が適用できる。全方向 $\xi,\eta$ について $\chi_\xi=\sum_i\xi_i\chi_i$ と書けば

$$
\eta^\top\bar A\xi
=\int_Y(\eta+D\chi_\eta)^\top A(\xi+D\chi_\xi)
$$

である。追加した $D\chi_\eta$ の項はセル式により零である。従って $\bar A$ は対称で、$\xi=\eta$ として楕円性と周期勾配の平均零を使うと $\xi^\top\bar A\xi\ge\lambda|\xi|^2$ である。

周期的な $q\in L^2(Y)$ は $q(\cdot/\eps)\rightharpoonup \int_Y q$ が局所 $L^2$ で成り立つ。滑らかな台付き試験関数との積分を各 $\eps$ セルに分け、試験関数をセルの一点の値で置換すると、一様連続性とセル上の $\int_Y|q|$ によって誤差は零に収束する。残りは $\int_Yq$ を掛けた Riemann 和である。コンパクト集合を覆うセル数は $O(\eps^{-n})$ なので局所 $L^2$ ノルムは一様有界であり、稠密性から弱収束が従う。
:::

:::proof[極限の流束の同定]
$u_\eps$ を試験関数にすると $\norm{Du_\eps}_2$ が一様有界である。Rellich のコンパクト性により部分列で $u_\eps\rightharpoonup u$ が $H_0^1$ で、$u_\eps\to u$ が $L^2$ で成り立つ。流束 $F_\eps=A(x/\eps)Du_\eps$ は $L^2$ で有界なので、弱極限を $F$ とおけば $-\divg F=f$ である。

$w_i^\eps(x)=x_i+\eps\chi_i(x/\eps)$ とし、$J_i^\eps=A(x/\eps)Dw_i^\eps$ とする。$w_i^\eps\to x_i$ は局所 $L^2$ で強収束し、$J_i^\eps\rightharpoonup\bar A e_i$ は局所 $L^2$ で弱収束し、$\divg J_i^\eps=0$ である。$\phi\in C_c^\infty(\Omega)$ に対し対称性から

$$
\begin{aligned}
\int\phi\, F_\eps\cdot Dw_i^\eps
&=\int f\phi w_i^\eps-\int w_i^\eps F_\eps\cdot D\phi,\\
\int\phi\, F_\eps\cdot Dw_i^\eps
&=\int\phi\, J_i^\eps\cdot Du_\eps=-\int u_\eps J_i^\eps\cdot D\phi
\end{aligned}
$$

である。最初の等式は弱式に $\phi w_i^\eps$ を入れ、次は $\divg J_i^\eps=0$ に $\phi u_\eps$ を入れた。各右辺の極限は片方が強収束するので取れる。$-\divg F=f$ を再び使うと最初の極限は $\int\phi F_i$、二番目の極限は $\int\phi\bar A e_i\cdot Du$ である。よって $F=\bar A Du$ である。極限問題は強圧性により一意なので、元の列の全てが同じ極限を持つ。
:::

この議論が直接与えるのは、変位の強収束と勾配や流束の弱収束である。例えば層に垂直な向きの $1+\chi_i'(y)$ は $y$ に依存するので、勾配の振動は変位の極限だけには残らない。勾配を強いノルムで近似するには、セル補正を含め、境界付近にその補正をどう接続するかを調べる。

:::example[層状の係数]
二次元の $A=a(y_1)I$ では $e_1$ 方向の流束を一定にするため $1+\chi_1'=\bar a/a$、$e_2$ では $\chi_2=0$ がセル問題を満たす。従って $\bar A$ の対角成分は調和平均と算術平均である。二つの層で $a=1,4$ が同じ厚さなら、積分を区分ごとに読んでそれぞれ $8/5,5/2$ になる。流れが層を横切るか層に沿うかで有効係数が異なる。
:::

:::example[二つの尺度で見える勾配]
一次元で $f=0$、$u_\eps(0)=0$、$u_\eps(1)=1$ とし、各周期の前半で $a=1$、後半で $a=4$ とする。$\eps=1/m$ なら区間に周期が整数個入るので、流束は正確に $8/5$ である。各層の勾配は $8/5$ と $2/5$ で、平均は $1$ である。従って $u_\eps\to x$ は一様だが、$\int_0^1|u_\eps'-1|^2=(3/5)^2$ はどの $m$ でも変わらない。変位の強収束から勾配の強収束は従わない。
:::

### 層状媒質の二つの有効係数

$A(y)=a(y_1)I$、$a$ は周期 $1$、正の上下界を持つとする。第 $1$ 方向で $\chi_1$ を $y_1$ だけの関数とすればセル方程式は $(a(1+\chi_1'))'=0$ である。従って流束は定数 $J$。平均を取ると $1=J\int_0^1 a^{-1}$ なので $J=(\int_0^1 a^{-1})^{-1}$ である。一方 $j>1$ では $\divg(a e_j)=\partial_j a=0$ なので $\chi_j=0$ が解である。このとき平均流束は $(\int_0^1 a)e_j$ となる。セル解の一意性は周期的零平均の差を試験し、勾配が零であることから従う。

:::example[同じ材料から異方性が出る]
周期の前半で $a=1$、後半で $a=9$ とすれば、層に直交する方向は $J=(1/2+1/18)^{-1}=9/5$、平行方向は $5$ である。媒質の各点で行列はスカラー倍だが、セルの解は直交方向だけ勾配を再配分する。例えば平均勾配を $1$ としたとき、各層での局所勾配は $9/5$ と $1/5$ で、両者の平均は $1$、流束はどちらも $9/5$ である。
:::

### 有効係数の変分表示と上下界

周期セル $Y=(0,1)^n$ で対称な $A$ が $\lambda I\le A\le\Lambda I$ を満たすとする。平均零の周期関数 $w$ に対し
$$
J_\xi(w)=\int_Y (\xi+Dw)\cdot A(\xi+Dw)
$$
を最小化する。セル解 $\chi_\xi$ の弱式は $\int_Y A(\xi+D\chi_\xi)\cdot D\phi=0$ である。従って $w=\chi_\xi+\phi$ として展開すると、差は $\int_Y D\phi\cdot A D\phi\ge0$ であり、一意な最小点になる。ここで試験関数 $\chi_\xi$ を使って最小値から勾配の一項を消せば $\min J_\xi=\xi\cdot A_{\rm hom}\xi$ である。線形性により $\chi_\xi$ は方向ごとのセル解の線形結合であり、二方向 $\xi,\eta$ の内積を
$$
\eta\cdot A_{\rm hom}\xi=\int_Y (\eta+D\chi_\eta)\cdot A(\xi+D\chi_\xi)
$$
と書ける。従って $A_{\rm hom}$ は対称である。

上界は $w=0$ を入れて $\xi\cdot A_{\rm hom}\xi\le\xi\cdot(\int_Y A)\xi$。下界は任意の定ベクトル $p$ について $z\cdot Az\ge2p\cdot z-p\cdot A^{-1}p$ とする。これは $(z-A^{-1}p)\cdot A(z-A^{-1}p)\ge0$ の展開である。$z=\xi+Dw$ を入れて周期性で $\int Dw=0$ を使い、$p$ について二次式を最大化すると
$$
\xi\cdot\left(\int_Y A^{-1}\right)^{-1}\xi
\le \xi\cdot A_{\rm hom}\xi .
$$
層状のスカラー係数では直交方向で下界、平行方向で上界を実現する。多次元の一般セルでは $A^{-1}p$ の振動部分が勾配として表せるとは限らないため、下界の等号をそのまま仮定できない。

### 一次元の周期係数で収束の速さを計算する

$-(a(x/\eps)u_\eps')'=f(x)$、$u_\eps(0)=u_\eps(L)=0$ とし、$a$ は周期 $1$、$0<\lambda\le a\le\Lambda$、$f\in L^1(0,L)$ とする。$b=1/a$、平均 $\bar b$、$F(x)=\int_0^x f$ と置く。流束は $C_\eps-F(x)$ であり、
$$
C_\eps=\frac{\int_0^L b(x/\eps)F(x)\dd x}{\int_0^L b(x/\eps)\dd x}
$$
となる。$B(y)=\int_0^y(b-\bar b)$ は有界周期関数なので $b(x/\eps)-\bar b=(\eps B(x/\eps))'$ である。従って任意の $G\in W^{1,1}(0,L)$ に対し
$$
\left|\int_0^z (b(x/\eps)-\bar b)G(x)\dd x\right|
\le \eps\norm B_\infty (|G(z)|+|G(0)|+\norm{G'}_1)
$$
が $0\le z\le L$ に一様に成り立つ。分母は $L/\Lambda$ 以上で、分子には $G=F$、分母には $G=1$ を用いれば $|C_\eps-C_0|\le C\eps$、$C_0=L^{-1}\int_0^L F$ である。さらに解の式 $u_\eps(z)=\int_0^z b(x/\eps)(C_\eps-F(x))\dd x$ に同じ評価を用い、$\norm{u_\eps-u_0}_\infty\le C\eps$ を得る。ここで $u_0'=\bar b(C_0-F)$ なので有効係数は $\bar b^{-1}$ である。

これに対し勾配の差には $(b(x/\eps)-\bar b)(C_0-F)$ が残るため、一般には同じ $L^2$ 強収束は出ない。流束は $C_\eps-F$ そのものなので強収束する。変位、勾配、流束を区別して収束を述べるのは、振動する係数の影響が三者で違うためである。

### 正弦状に変わる逆係数を積分する

区間 $(0,1)$ で $a(x/\eps)^{-1}=2+\sin(2\pi x/\eps)$、$f=0$、$u(0)=0,u(1)=1$ とする。流束は定数 $C_\eps$ で
$$
u_\eps(x)=C_\eps\left(2x+\frac{\eps}{2\pi}(1-\cos(2\pi x/\eps))\right).
$$
$\eps=1/N$（$N\in\N$）に沿えば端点の余弦は $1$ なので $C_\eps=1/2$。従って $u_\eps-x=(\eps/(4\pi))(1-\cos(2\pi Nx))$ であり、一様ノルムは厳密に $\eps/(2\pi)$ である。一方 $u_\eps'-1=\tfrac12\sin(2\pi Nx)$ の $L^2$ ノルムは $1/(2\sqrt2)$ で、零へ行かない。流束は各 $N$ で既に有効係数 $1/2$ と等しい。従って三つの量の収束の違いを一つの解から比較できる。

同じ周期係数を二次元の層状係数 $A=a(y_1)I$ に使うと、平行方向の有効係数は $\int_0^1(2+\sin(2\pi s))^{-1}\dd s=1/\sqrt3$ で、直交方向の $1/2$ と異なる。積分では角をずらして $\sin$ を $\cos$ にし、$z=\tan(\theta/2)$ で一周期を実直線へ写す。半角公式から $\int_{-\pi}^\pi\dd\theta/(2+\cos\theta)=\int_{-\infty}^\infty2\dd z/(3+z^2)=2\pi/\sqrt3$ である。従って有効な異方性は不連続な材料の境界を用いなくても現れる。

## 演習問題

::::exercise{level=basic}
$(0,2)$ の斉次 Dirichlet 条件で、空間作用素 $-\partial_{xx}$ の第 2 固有値を求めよ。

:::solution
端点条件により $\sqrt\lambda L=m\pi$ である。$L=2,m=2$ を入れると $\lambda=\pi^2$ である。
:::
::::

::::exercise{level=basic}
矩形 $(0,\pi)\times(0,1)$ で上辺のみ値 $\sin x$ を持つ調和関数を求めよ。

:::solution
分離して $Y''=Y$、$Y(0)=0$、$Y(1)=1$ となるため $u=\sin x\sinh y/\sinh1$ である。左右の境界値 $0$ と $\Delta u=0$ を満たす。
:::
::::

::::exercise{level=basic}
$\kappa=2$ の熱方程式で $\xi=3/(2\pi)$ の Fourier モードの減衰指数を求めよ。

:::solution
$\widehat u_t=-4\pi^2\kappa\xi^2\widehat u$ なので指数は $-18t$ である。
:::
::::

::::exercise{level=basic}
$a=1,9$ の二つの層が同じ厚さのとき、層を横切る有効係数を求めよ。

:::solution
逆係数の平均は $(1+1/9)/2=5/9$ なのでその逆数 $9/5$ が有効係数である。
:::
::::

::::exercise{level=advanced}
式 :eqref[eq-pde-viscous-burgers-wave] で $a=4,b=0,\nu=1$ の中心勾配と速度を求めよ。

:::solution
$U=4/(1+e^{2\xi})$ なので $U'=-8e^{2\xi}/(1+e^{2\xi})^2$ であり、$U'(0)=-2$ である。速度は $(4+0)/2=2$ である。
:::
::::

::::exercise{level=advanced}
$U(0,x)=2x^2$、$H=p^2/2$ を Legendre 変換で解け。

:::solution
$V(0,p)=p^2/8$ なので $V(t,p)=(1+4t)p^2/8$ である。逆変換では $p=4x/(1+4t)$ を使い、$U=2x^2/(1+4t)$ を得る。
:::
::::

::::exercise{level=advanced}
$n(x)=e^x$ の一次元幾何光学で $\phi'>0$ を使い、振幅と残差を求めよ。

:::solution
$\phi'=e^x$ だから $\phi=e^x$ とできる。$(a^2\phi')'=0$ なので $a=e^{-x/2}$ を選べる。式 :eqref[eq-pde-wkb-residual] は $\eps^2 e^{-x/2}/4$ に等しい。
:::
::::

::::exercise{level=advanced}
一次元で $a$ を正の上下界を持つ可測な $1$ 周期係数とする。$-(a(x/\eps)u_\eps')'=0$、$u_\eps(0)=0$、$u_\eps(1)=1$ を解き、一様極限を求めよ。

:::solution
流束は一定だから $u_\eps(x)=\int_0^x a(s/\eps)^{-1}\dd s/\int_0^1 a(s/\eps)^{-1}\dd s$ である。分子は $x\bar b+O(\eps)$、分母は $\bar b+O(\eps)$ で $\bar b>0$ なので $u_\eps\to x$ が一様に従う。
:::
::::

## 確認テスト

::::quiz
分離解を重ね合わせられる直接の理由はどれか。

- [x] 方程式が線形である
- [ ] 固有値が全て等しい
- [ ] 分離定数が存在しない

:::explanation
線形作用素では $L(u+v)=Lu+Lv$ であり、斉次方程式の二つの解の和も解である。境界条件と級数の収束は別に確認する。
:::
::::

::::quiz{answer="2" tol="1e-6"}
粘性 Burgers の左状態 $4$、右状態 $0$ の波の速度を求めよ。

:::explanation
端の状態では $\nu U'=0$ なので保存型の流束差を状態差で割る。結果は $(4+0)/2=2$ である。
:::
::::

::::quiz
幾何光学で最大次数の項を消して得るものはどれか。

- [x] 位相のアイコナール方程式
- [ ] 振幅の二階微分方程式
- [ ] 変分問題の境界条件

:::explanation
指数関数を二度微分して生じる最高次数の項は $-|D\phi|^2a$ に比例する。屈折率項と合わせて $|D\phi|^2=n^2$ を得る。
:::
::::

::::quiz{answer="1.6" tol="1e-6"}
係数 $1$ と $4$ が同じ厚さで直列に並ぶ媒質の有効係数を求めよ。

:::explanation
平均逆係数は $(1+1/4)/2=5/8$ なので、その逆数 $8/5=1.6$ である。
:::
::::

::::quiz{answer="1.8" tol="1e-6"}
一次元の周期係数 $a$ が一周期の前半で $1$、後半で $9$ のとき、有効係数を求めよ。

:::explanation
セルの流束 $a(1+\chi^{\prime})$ を定数 $J$ とおく。$\chi^{\prime}$ の周期平均が零なので、$1=J(\tfrac12/1+\tfrac12/9)=5J/9$ である。従って $J=9/5$。有効係数は勾配 $1$ に対する流束なので $9/5$ である。
:::
::::


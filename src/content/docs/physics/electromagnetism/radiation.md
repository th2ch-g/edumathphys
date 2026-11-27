---
title: 電磁放射
description: 点電荷の遅延場から加速による放射を取り出し、局在電流の多重極放射、アンテナの指向性、自由電子と小粒子の散乱を導く。
notation:
  - symbol: 't_*'
    meaning: 点電荷の遅延時刻
  - symbol: '\kappa'
    meaning: 遅延時刻を微分するときの係数
  - symbol: '\vb{n}'
    meaning: 源から観測点への単位ベクトル
  - symbol: 'Q_{ij}'
    meaning: 跡を除いた電気四重極モーメント
  - symbol: '\sigma_s'
    meaning: 全散乱断面積
---

## この章の目標

- 運動する点源で遅延時刻の Jacobian を計算する。
- 近接場と放射場を距離依存性で区別する。
- 電気双極子・磁気双極子・電気四重極の放射電力を導く。
- アンテナの放射抵抗と粒子の散乱断面積を強度から求める。

## Liénard–Wiechert ポテンシャル

### 動いている源のどの位置を見るか

電荷 $q$ が軌道 $\vb{z}(t')$ を速度 $\vb{v}$ で進む。
観測点 $\vb{r}$ に時刻 $t$ に届く光は

$$
t-t_*=\frac{R(t_*)}{c},\qquad R(t')=|\vb{r}-\vb{z}(t')|,\qquad \vb{n}(t')=(\vb{r}-\vb{z}(t'))/R(t')
\label{eq-em-retarded-time}
$$

を満たす時刻 $t_*$ に出た。
以後点源の式の右辺はすべて $t_*$ で評価する。
$\vb{\beta}=\vb{v}/c$、

$$
\kappa=1-\vb{n}\cdot\vb{\beta}
$$

とすると $\dd(t'+R/c)/\dd t'=\kappa>0$ である。

従って亜光速で過去光円錐と交わる軌道では交点は一つである。

自由空間の遅延核 :eqref[eq-em-retarded-green] に $\rho=q\delta^3(\vb{r}'-\vb{z}(t'))$、$\vb{J}=\rho\vb{v}(t')$ を入れる。
空間を先に積分すると $\delta(t-t'-R(t')/c)/R(t')$ の時間積分が残る。

変数変換 $\delta(f)=\delta(t'-t_*)/|f'(t_*)|$ が分母に $\kappa$ を与える。
従って :term[Liénard–Wiechert ポテンシャル]{en="Liénard–Wiechert potentials" yomi="りえなーるうぃーへるとぽてんしゃる"} は

$$
\phi=\frac{q}{4\pi\varepsilon_0\kappa R},\qquad \vb{A}=\frac{\vb{v}}{c^2}\phi
\label{eq-em-lw}
$$

となる。
静止電荷の Coulomb 電位に遅延位置を代入するだけでは観測時間と放出時間の対応が漏れる。

### 速度場と加速度場

遅延条件を $t$ と $\vb{r}$ で暗黙微分すると $\partial_t t_*=1/\kappa$、$\nabla t_*=-\vb{n}/(c\kappa)$ である。

$$
D=\kappa R=R-(\vb{r}-\vb{z})\cdot\vb{\beta}
$$

と置く。

源時刻を固定した空間勾配は $\vb{n}-\vb{\beta}$ であり、観測位置を固定した源時刻微分は

$$
D'=-c\,\vb{n}\cdot\vb{\beta}+c\beta^2-R\,\vb{n}\cdot\dot{\vb{\beta}}.
$$

従って全勾配 $\nabla D=\vb{n}-\vb{\beta}-D'\vb{n}/(c\kappa)$、$\partial_tD=D'/\kappa$ である。
$\vb{E}=-\nabla\phi-\partial_t(\vb{\beta}\phi/c)$ に代入すれば、$k_eq/D^2$ をくくった括弧は $(\vb{n}-\vb{\beta})(1-D'/(c\kappa))-D\dot{\vb{\beta}}/(c\kappa)$ となる。

速度と加速度に分けて三重積を戻すと

$$
\vb{E}=k_eq\left[\frac{(1-\beta^2)(\vb{n}-\vb{\beta})}{\kappa^3R^2}+\frac{\vb{n}\times[(\vb{n}-\vb{\beta})\times\dot{\vb{\beta}}]}{c\kappa^3R}\right],\qquad \vb{B}=\frac1c\vb{n}\times\vb{E}.
\label{eq-em-lw-field}
$$

磁場は $\nabla\times(\vb{\beta}\phi/c)=[\nabla\phi\times\vb{\beta}+\phi\,\nabla t_*\times\dot{\vb{\beta}}]/c$ に上の微分を入れれば得られる。
$\dot{\vb{\beta}}=\vb{a}/c$ は源時刻についての微分である。

速度場は $R^{-2}$、加速度場は $R^{-1}$ である。
遠方の球面を横切る電力を積分すると、面積の $R^2$ に対して有限に残るのは後者の二乗だけである。

:::example[等速直線運動]
瞬時位置からの変位を

$$
\vb{s}=\vb{r}-\vb{z}(t)=R(\vb{n}-\vb{\beta})
$$

とすると、

$$
R^2\kappa^2=s^2-|\vb{\beta}\times\vb{s}|^2
$$

である。
従って等速場は

$$
\vb{E}=k_eq(1-\beta^2)\vb{s}/[s^2-|\vb{\beta}\times\vb{s}|^2]^{3/2}
$$

となり、瞬時位置から放射状に向く。

横方向で静止 Coulomb 場の $\gamma$ 倍、縦方向で $1/\gamma^2$ 倍である。
加速度が零なので放射電力も零である。
:::

:::example[パルスの到着間隔]
源が観測者に向かって $0.80c$ で動く瞬間には $\kappa=0.20$ である。
源の近接する二時刻の間隔が $1.00\,\mathrm{ns}$ なら、到着間隔は一次近似で $0.200\,\mathrm{ns}$ になる。
放出単位時間あたりの電力を求めるときにも同じ係数が必要になる。
:::

### 遅延位置を方程式から求める例

一定速度 $v$ で $x$ 軸を進む電荷が $t=0$ に原点を通るとする。
観測点を $(0,b,0)$ とし、最接近の実験室時刻 $t=0$ で測る。
遅延条件は $-ct_*=\sqrt{v^2t_*^2+b^2}$ である。
両辺を二乗すると $t_*=\pm\gamma b/c$ と出るが、光が過去から届く条件によって負根を選ぶ。

従って $R=\gamma b$、$\vb{n}=(\beta,1/\gamma,0)$、$\kappa=1-\beta^2$、$\kappa R=b/\gamma$ である。
電位は $k_eq\gamma/b$、電場は

$$
\vb{E}=k_eq\gamma\hat{\vb{y}}/b^2
$$

となる。

遅延位置は左にずれているにもかかわらず、等速場が最接近点で純粋に横を向くことを直接確かめられる。

### 放射パルスを観測時刻で積分する

加速度が短い源時刻区間だけ非零の粒子を十分遠くで観測する。観測方向 $\vb n$ は区間を通じて一定とし、$1/R$ の係数に現れる $R$ も一定とみなすが、遅延時刻の微分 $\dd t=(1-\vb n\cdot\boldsymbol\beta)\dd t_*$ は保つ。このとき電場の時間積分は粒子が区間の途中でどう加速したかによらず、両端の速度だけで表せる。$\kappa=1-\vb n\cdot\boldsymbol\beta$ と置いて商を微分すると

$$
\dv{}{t_*}\frac{\vb n\times(\vb n\times\boldsymbol\beta)}{\kappa}
=\frac{\vb n\times[(\vb n-\boldsymbol\beta)\times\dot{\boldsymbol\beta}]}{\kappa^2}.
$$

実際、分子は $\vb n(\vb n\cdot\boldsymbol\beta)-\boldsymbol\beta$ であり、その微分に $\kappa$ を掛けてから、分母の微分に由来する項を加えれば右辺の三重積になる。Liénard–Wiechert の加速度場に $\dd t$ を掛けると分母が $\kappa^3$ から $\kappa^2$ になるので、

$$
\int \vb E_{\rm rad}\dd t=\frac{k_eq}{cR}
\left[\frac{\vb n\times(\vb n\times\boldsymbol\beta)}{1-\vb n\cdot\boldsymbol\beta}\right]_{\rm before}^{\rm after}.
$$

受信アンテナがこの時間積分に比例する信号を返すなら、始めと終わりが同じ速度のパルスでは正負の波形が相殺する。しかし放射エネルギーの積分には $|\vb E|^2$ が入るため、同じ条件でも正の値を取りうる。極めて低い周波数の成分を推定するとき、電場の面積と電力の面積を区別する必要がある。


:::example[観測角度とパルスの時間圧縮]
速度と観測方向のなす角を $\theta$ とすれば、短い源時刻間隔 $\Delta t_*$ の到着間隔は $\Delta t=(1-\beta\cos\theta)\Delta t_*$ である。有限区間にわたる平均をこの微分式で近似するには、区間内で角度と速度の変化が小さい必要がある。

$x=\theta$ をラジアンで測り、$\beta=0.90$、$Y_*=\Delta t_*=1.00\,\mathrm{ns}$、$Y=\Delta t$ とする。

$$
f(x)=- \frac{9 \cos{\left(x \right)} - 10}{10},\qquad Y=Y_* f(x).
$$

$x=0$ のまわりで先頭の二項を残した近似は

$$
f_{\rm a}(x)=\frac{9 x^{2}}{20} + \frac{1}{10}.
$$

次の表の差は $Y_*[f(x)-f_{\rm a}(x)]$ である。

| $x$ | $f(x)$ | $Y$ $[\mathrm{ns}]$ | 厳密値と近似値の差 $[\mathrm{ns}]$ |
|---:|---:|---:|---:|
| $0$ | $0.1$ | $0.1$ | $0$ |
| $\frac{1}{10}$ | $0.1045$ | $0.1045$ | $-3.749\times10^{-6}$ |
| $\frac{3}{10}$ | $0.1402$ | $0.1402$ | $-3.028\times10^{-4}$ |
| $\frac{3}{5}$ | $0.2572$ | $0.2572$ | $-0.004802$ |

前方の時間圧縮は有限であり、$1-\beta$ が小さいときほど小さな角度の変化が相対的に強く効く。領域を角度だけで指定せず、$\theta^2$ と $1-\beta$ の大きさを比べると近似の意味がわかる。
:::

## 双極子放射

### 局在源で二つの近似を順番に行う

源の半径を $a$、観測距離を $r$、$k=\omega/c$ とする。
遠方の $1/r$ 項を取るには $r\gg a$ で振幅を $1/r$ とし、位相の曲率を無視するには $ka^2/r\ll1$ とする。
さらに放射が近接場より優勢な波動域 $kr\gg1$ を考える。

$\vb{n}=\vb{r}/r$、$t_0=t-r/c$ とおくと遅延ポテンシャルは

$$
\vb{A}_{\rm far}=\frac{\mu_0}{4\pi r}\int\vb{J}(\vb{r}',t_0+\vb{n}\cdot\vb{r}'/c)\,\dd^3r'.
$$

ここではまだ $ka\ll1$ を仮定していない。
$1/r$ 項の空間微分は遅延時刻の微分から来るので $\vb{B}_{\rm rad}=-\vb{n}\times\partial_t\vb{A}_{\rm far}/c$、

$$
\vb{E}_{\rm rad}=c\vb{B}_{\rm rad}\times\vb{n}=-\partial_t\vb{A}_{{\rm far},\perp}
$$

である。

ここで

$$
\vb{X}_\perp=\vb{X}-\vb{n}(\vb{n}\cdot\vb{X})
$$

とする。

さらに $ka\ll1$ なら源内の時間差で展開できる。
連続の式に $x_i'$ を掛けて積分し、表面項を消すと $\dot p_i=\int J_i\,\dd^3r'$ である。
よって先頭の :term[電気双極子放射]{en="electric dipole radiation" yomi="でんきそうきょくしほうしゃ"} は

$$
\vb{E}_{p}=-\frac{\ddot{\vb{p}}_\perp(t_0)}{4\pi\varepsilon_0 c^2 r},\qquad \frac{\dd P_p}{\dd\Omega}=\frac{|\vb{n}\times\ddot{\vb{p}}|^2}{16\pi^2\varepsilon_0 c^3},\qquad P_p=\frac{|\ddot{\vb{p}}|^2}{6\pi\varepsilon_0 c^3}.
\label{eq-em-dipole-radiation}
$$

角度分布は $r^2\varepsilon_0 c E_p^2$ であり、双極子軸を極軸として $\int\sin^2\theta\,\dd\Omega=8\pi/3$ を使った。
単色のピーク振幅で時間平均をとる場合はさらに $1/2$ を掛ける。

### 次の次数にある二つの放射

時間差の一次項を求める。
$M_{ij}=\int x'_j J_i\,\dd^3r'$ とし、$D_{ij}=\int\rho x_i'x_j'\,\dd^3r'$、$Q_{ij}=3D_{ij}-\delta_{ij}D_{kk}$、$\vb{m}=\tfrac12\int\vb{r}'\times\vb{J}\,\dd^3r'$ と定める。

連続の式を部分積分すれば $\dot D_{ij}=M_{ij}+M_{ji}$、反対称部分は

$$
M_{ij}-M_{ji}=-2\epsilon_{ijk}m_k
$$

となる。
従って

$$
M_{ij}n_j=\dot D_{ij}n_j/2+(\vb{m}\times\vb{n})_i
$$

である。

横射影で跡の項は消えるので、放射場はこの次数まで

$$
\vb{E}_{\rm rad}=-\frac{\mu_0}{4\pi r}\left[\ddot{\vb{p}}_\perp+\frac{\ddot{\vb{m}}\times\vb{n}}{c}+\frac{(\dddot{\vb{Q}}\vb{n})_\perp}{6c}\right].
\label{eq-em-higher-radiation}
$$

すべての量は $t_0$ で評価する。
磁気双極子単独の全電力は電気双極子で $\ddot p$ を $\ddot m/c$ へ置き換えた値である。
四重極の積分には、任意の実対称・無跡な $X_{ij}$ について

$$
\int |(\vb{X}\vb{n})_\perp|^2\,\dd\Omega=\int[|\vb{X}\vb{n}|^2-(\vb{n}\cdot\vb{X}\vb{n})^2]\,\dd\Omega=\frac{4\pi}{5}\sum_{ij}X_{ij}^2
$$

を使う。
実際、等方性から $\int n_i n_j\,\dd\Omega=4\pi\delta_{ij}/3$、四次平均は三つの $\delta\delta$ の組の和の $4\pi/15$ 倍である。
係数は添字を縮約して $\vb{n}^2=1$ とすれば決まる。複素振幅を用いるときは一方の因子を複素共役とし、右辺の和を $\sum_{ij}|X_{ij}|^2$ とする。

よって

$$
P_m=\frac{|\ddot{\vb{m}}|^2}{6\pi\varepsilon_0 c^5},\qquad P_Q=\frac{\sum_{ij}(\dddot Q_{ij})^2}{720\pi\varepsilon_0 c^5}.
\label{eq-em-multipole-power}
$$

電気双極子と他の二項との干渉は角度を反転すれば積分で消える。
磁気双極子と四重極の干渉は $\epsilon_{ijk}$ と対称テンソルの縮約で消える。
角度別の強度には干渉が残りうる。
任意の源を同じ精度で展開する際には、先頭の双極子とさらに高次の振幅との干渉も次数を確認する。

:::example[磁気ループと軸対称な四重極]
面積 $A$ の十分小さな電流ループを $I_0\cos\omega t$ で駆動すれば $m_z=AI_0\cos\omega t$ であり、平均電力は $\omega^4 A^2 I_0^2/(12\pi\varepsilon_0c^5)$ である。

軸対称な四重極では $Q_{xx}=Q_{yy}=-Q_{zz}/2$ なので

$$
\sum_{ij}(\dddot Q_{ij})^2=3(\dddot Q_{zz})^2/2
$$

である。
軸上にも赤道面にも四重極の放射はない。
いずれも $\vb{Q}\vb{n}$ が $\vb{n}$ に平行になるからである。
:::

### 線状アンテナと放射抵抗

長さ $\ell\ll\lambda$ の細線にほぼ一様な電流 $I=I_0\cos\omega t$ が流れる理想電流要素を考える。
$\dot p=I\ell$ より平均放射電力は $\mu_0\omega^2 I_0^2\ell^2/(12\pi c)$ である。

$\langle P\rangle=I_0^2R_{\rm rad}/2$ で定めた :term[放射抵抗]{en="radiation resistance" yomi="ほうしゃていこう"} は

$$
R_{\rm rad}=Z_0(k\ell)^2/(6\pi)\simeq80\pi^2(\ell/\lambda)^2\,\Omega
$$

となる。

これは熱損失と区別される。
中央給電の短い開放線で電流が端に向けて三角形状に減るなら、$\int I(z)\dd z/I_0=\ell/2$ なので抵抗は一様電流の四分の一である。

半波長ダイポールでは長波長展開を使えない。
細線上で $I(z)=I_0\cos kz,\ |z|<\lambda/4$ と近似し、遠方積分の位相 $e^{-ikz\cos\theta}$ を残す。
余弦を指数の和にして積分すると

$$
\int I(z)e^{-ikz\cos\theta}\dd z=2I_0\cos[(\pi/2)\cos\theta]/(k\sin^2\theta)
$$

である。

横射影の $\sin\theta$ を掛けた遠方電場の振幅は

$$
|E_\theta|=\frac{Z_0 I_0}{2\pi r}\left|\frac{\cos[(\pi/2)\cos\theta]}{\sin\theta}\right|.
$$

平均強度を積分すれば

$$
R_{\rm rad}=(Z_0/2\pi)\int_0^\pi\cos^2[(\pi/2)\cos\theta]\,\dd\theta/\sin\theta\simeq73.1\,\Omega
$$

となる。
:term[指向性]{en="directivity" yomi="しこうせい"} $D=4\pi(\dd P/\dd\Omega)_{\max}/P$ は短い双極子で $1.5$、この半波長模型で約 $1.64$ となる。

### 振動する点双極子の近接場から遠方場まで

局在双極子の点源模型は $\rho=-\vb{p}(t)\cdot\nabla\delta^3(\vb{r})$、

$$
\vb{J}=\dot{\vb{p}}(t)\delta^3(\vb{r})
$$

である。
この組は連続の式を満たす。
遅延積分で空間微分を部分積分すると、電位には $1/r$ の微分と遅延時刻の微分が現れ、

$$
\phi=k_e\left[\frac{\vb{n}\cdot\vb{p}}{r^2}+\frac{\vb{n}\cdot\dot{\vb{p}}}{cr}\right],\qquad \vb{A}=\frac{\mu_0\dot{\vb{p}}}{4\pi r}
$$

となる。
右辺の各モーメントは $t_0=t-r/c$ にとる。
$\partial_i n_j=(\delta_{ij}-n_i n_j)/r$、$\partial_i t_0=-n_i/c$ を使って微分し、線形演算 $\mathcal D(\vb{b})=3\vb{n}(\vb{n}\cdot\vb{b})-\vb{b}$ を使うと、$r>0$ における全場は

$$
\begin{aligned}
\vb{E}&=k_e\left[\frac{\mathcal D(\vb{p})}{r^3}+\frac{\mathcal D(\dot{\vb{p}})}{cr^2}-\frac{\ddot{\vb{p}}_\perp}{c^2r}\right],\\
\vb{B}&=\frac{\mu_0}{4\pi}\left[\frac{\dot{\vb{p}}\times\vb{n}}{r^2}+\frac{\ddot{\vb{p}}\times\vb{n}}{cr}\right].
\end{aligned}
$$

単色源なら電場の三項の大きさは一般の観測角で $1:kr:(kr)^2$ の比である。
近接域の $r^{-3}$ 成分はほぼ静電場であり、$kr\sim1$ では三項の位相も考える必要がある。
$kr\gg1$ で横向きの放射項が残る。
双極子軸上では $\ddot{\vb{p}}_\perp=0$ なので遠くへ行ってもこの項が復活することはない。

:::example[同じ全電流でも給電位置と位相が指向性を変える]
二つの同じ短い双極子を $z$ 軸に平行に置き、中心を $x=\pm d/2$ に取る。
同相の駆動なら、遠方での和は一個の放射振幅の $2\cos(kd\,n_x/2)$ 倍である。
$d=\lambda/2$ のとき $x$ 軸上は二電流が打ち消し合い、$y$ 軸上では振幅が倍で強度が四倍になる。

給電位相を一方だけ $\pi$ 変えれば因子は $2i\sin(kd\,n_x/2)$ に変わり、この二方向の極大と零点が入れ替わる。
全放射電力を求めるには角度ごとの干渉も含めて積分する。
:::

### 大きさが波長に近い源の球面多重極

$ka\ll1$ でなくても、源を囲む球の外で角度を球面調和関数に分解できる。
時間因子を $e^{-i\omega t}$ とし、$\ell\ge1$ について、球面の勾配 $\nabla_\Omega$ を用いて接線ベクトル

$$
\vb{X}_{\ell m}=\frac{-i\,\vb{n}\times\nabla_\Omega Y_{\ell m}}{\sqrt{\ell(\ell+1)}}
$$

を定める。
内積の角度積分で部分積分と $-\nabla_\Omega^2Y_{\ell m}=\ell(\ell+1)Y_{\ell m}$ を使えば、$\vb{X}_{\ell m}$ は正規直交する。
$\vb{n}\times\vb{X}_{\ell m}$ も正規直交し、二族の交差内積は球面上の回転をもつ場と勾配の内積であるため部分積分で零になる。

動径関数 $z_\ell(kr)$ を球 Bessel 方程式

$$
\frac1{r^2}\partial_r(r^2\partial_r z_\ell)+\left(k^2-\frac{\ell(\ell+1)}{r^2}\right)z_\ell=0
$$

の解とする。
$\vb{M}_{\ell m}=z_\ell\vb{X}_{\ell m}$、

$$
\vb{N}_{\ell m}=k^{-1}\nabla\times\vb{M}_{\ell m}
$$

と置く。
$\vb{M}$ は $z_\ell Y_{\ell m}$ に角運動量演算子 $-i\vb{r}\times\nabla$ を掛けた定数倍である。

この演算子は Laplace 演算子と可換なので $(\nabla^2+k^2)\vb{M}=0$ である。
また $\vb{r}\times\nabla f$ の発散は零なので $\nabla\times\vb{N}=k\vb{M}$ となる。

従って

$$
(\vb{E},\vb{H})=(\vb{M},-i\vb{N}/Z_0),\qquad (\vb{N},-i\vb{M}/Z_0)
$$

はそれぞれ無源の Maxwell 方程式を満たす。
前者は $E_r=0$ の磁気型、後者は $H_r=0$ の電気型である。
$\ell=0$ では $\nabla_\Omega Y_{00}=0$ なので放射モードは存在しない。
原点での正則性を課す内部解には球 Bessel 関数 $j_\ell$、外向き波を課す外部解には球 Hankel 関数 $h_\ell^{(1)}$ を使う。

後者は大きな $kr$ で $(-i)^{\ell+1}e^{ikr}/(kr)$ に近づく。

球面で測った複素遠方電場を $\vb{E}=(e^{ikr}/r)\vb{f}(\Omega)+O(r^{-2})$ と書くと、$\vb{f}$ は接線場なので上の二族へ直交射影できる。
例えば磁気型係数は、動径関数の漸近的な定数を振幅へ吸収して $a_{\ell m}=\int\vb{X}_{\ell m}^*\cdot\vb{f}\,\dd\Omega$ と定められる。

電気型は $\vb{n}\times\vb{X}_{\ell m}$ への同じ射影で $b_{\ell m}$ を得る。
球面上で滑らかな接線場を勾配とその $90^\circ$ 回転に分けるには、発散と面内回転をそれぞれ球面調和関数に展開して $\ell(\ell+1)$ で割ればよい。

球面上にはこの二条件がともに零の非零接線場はない。
実際に回転が零なら球が単連結なので勾配と書け、発散も零ならそのポテンシャルは極値原理により定数である。
従って展開は完全である。
直交性を用いると全平均電力は

$$
\langle P\rangle=\frac1{2Z_0}\int|\vb{f}|^2\dd\Omega=\frac1{2Z_0}\sum_{\ell m}(|a_{\ell m}|^2+|b_{\ell m}|^2).
$$

小さな源について最初に求めた二つの双極子はそれぞれの $\ell=1$、電気四重極は電気型の $\ell=2$ に対応する。

### 半波長アンテナの電流分布と放射抵抗

細い直線状の中心給電アンテナを理想化し、全長 $2l$、軸方向のピーク電流を $I(z)=I_f\sin[k(l-|z|)]/\sin(kl)$ とする。端で電流が消え、中心で給電電流 $I_f$ になる定在波近似である。導体の有限半径、給電隙間、導電損失による分布の補正は無視する。軸を $z$ とすると遠方場は、電流の横射影から因子 $\sin\theta$、位置による遅延から $e^{-ikz\cos\theta}$ を含む。$u=\cos\theta$ とし、正負の $z$ をまとめて二つの正弦・余弦の積を積分すれば

$$
\int_{-l}^l I(z)e^{-ikuz}\dd z
=\frac{2I_f}{k\sin kl}\frac{\cos(klu)-\cos(kl)}{1-u^2}.
$$

ここでは $\int_0^l\sin[k(l-z)]\cos(kuz)\dd z$ を部分積分しても同じ結果を得る。$2l=\lambda/2$ で $kl=\pi/2$ だから遠方電場の大きさは $Z_0 |I_f|\,|\cos(\pi u/2)|/(2\pi r\sqrt{1-u^2})$ である。軸上で分母も零になるが、分子は $1-|u|$ に一次で消えるため場は零となる。全放射電力を $|I_f|^2 R_r/2$ と書くと

$$
R_r=\frac{Z_0}{2\pi}\int_{-1}^1\frac{\cos^2(\pi u/2)}{1-u^2}\dd u
=\frac{Z_0}{4\pi}\int_0^{2\pi}\frac{1-\cos t}{t}\dd t.
$$

最後は部分分数 $1/(1-u^2)=[1/(1-u)+1/(1+u)]/2$ と偶性を使い、$t=\pi(1-u)$ と置いた。被積分関数は $t=0$ で有限に零へ近づく。$Z_0\simeq376.73\,\Omega$ を入れれば $R_r=73.1\,\Omega$ である。最大強度は赤道上であり、指向性は $Z_0/(\pi R_r)\simeq1.64$ となる。放射抵抗は平均電力を給電点の電流で換算した値で、線材の発熱抵抗とは別に求める。

### 多重極が干渉するときの後方放射

電気双極子の遠方電場と磁気双極子の遠方電場を同時に測る場合、位相を失わずにベクトルを足す必要がある。単色ピーク振幅を用い、基準の電気双極子モーメントを $p_*$ とする。$\vb d=\vb p/p_*$、$\vb h=\vb m/(cp_*)$、対称無跡行列 $Q$ の無次元化を $W=-ikQ/(6p_*)$ と書けば、共通の遠方係数を取り除いた電場は

$$
\vb f(\vb n)=\vb d-\vb n(\vb n\cdot\vb d)
+\vb h\times\vb n+W\vb n-\vb n(\vb n\cdot W\vb n).
$$

符号は時間因子 $e^{-i\omega t}$ に対応する。角分布は $\dd\bar P/\dd\Omega=P_*|\vb f|^2$、$P_*=\omega^4 p_*^2/(32\pi^2\varepsilon_0c^3)$ である。全立体角では異なる種類の多重極の干渉は消えるが、片側の半球だけを集めると消えない場合がある。以下では $p_*=1.00\times10^{-12}\,\mathrm{C\,m}$、周波数 $100\,\mathrm{MHz}$ を用いる。立体角に比例する係数を測るには、方位角平均と極角の積分を区別する。

:::example[回転する電気双極子]
同じ強さの二つの直交電気双極子を四分の一周期ずらして駆動する。実空間でベクトルの先端が円を描くため、一個の直線双極子のような軸上の暗点が残るかを調べる。無次元の素モーメントを選ぶと

$$
\begin{aligned}
\vb d&=\left(\begin{matrix}1 & i & 0\end{matrix}\right),\quad \vb h=\left(\begin{matrix}0 & 0 & 0\end{matrix}\right),\\
W&=\left(\begin{matrix}0 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{matrix}\right).
\end{aligned}
$$

極軸での場をまず調べる。$\vb n=\pm\hat{\vb z}$ では横射影は第三成分を落とす演算であり、電気双極子の項は方向反転で変わらず、残りの二項は符号を変える。

$$
\begin{aligned}
\vb f(+\hat{\vb z})&=\left(\begin{matrix}1 & i & 0\end{matrix}\right),\\
\vb f(-\hat{\vb z})&=\left(\begin{matrix}1 & i & 0\end{matrix}\right).
\end{aligned}
$$

半球全体を見るため、$u=\cos\theta$ と置き $X=\sqrt{1-u^2}\cos\varphi$、$Y=\sqrt{1-u^2}\sin\varphi$ を代入する。$X,Y$ の奇数乗を一つでも含む項は方位角平均で零となる。残る項をまとめて $F(u)=\langle|\vb f|^2\rangle_\varphi$ と書くと

$$
F(u)=u^{2} + 1.
$$

例えば二次の方位角平均は $\langle X^2\rangle=\langle Y^2\rangle=(1-u^2)/2$ である。立体角は $\dd\varphi\,\dd u$ なので、正負の半球を別々に積分すると

$$
\begin{aligned}
P_+/(2\pi P_*)&=\int_0^1 F(u)\dd u=\frac{4}{3},\\
P_-/(2\pi P_*)&=\int_{-1}^0 F(u)\dd u=\frac{4}{3},\\
\langle \cos\theta\rangle&=\frac{\int_{-1}^1 uF(u)\dd u}{\frac{8}{3}}=0.
\end{aligned}
$$

角度分解した検出器の較正値も示す。以下は各極角における方位角平均であり、極方向の開口幅を一定としただけの測定ではさらに $\sin\theta$ の立体角重みが必要である。

| 極角 | $F$ | 単位立体角あたりの平均放射電力 $[\mathrm{\mu W/sr}]$ |
|---:|---:|---:|
| $0^\circ$ | $2$ | $4.137$ |
| $45^\circ$ | $\frac{3}{2}$ | $3.103$ |
| $90^\circ$ | $1$ | $2.069$ |
| $135^\circ$ | $\frac{3}{2}$ | $3.103$ |
| $180^\circ$ | $2$ | $4.137$ |

全電力は $34.66\,\mathrm{\mu W}$ である。観測方向の極軸では横方向に回転する場がそのまま見える。二つの直線双極子を非干渉で足しても強度分布は同じになりうるが、偏光の測定で区別できる。
:::

:::example[電気・磁気双極子による一方向の暗点]
同位相の電気双極子と磁気双極子を直交させる。両者の規格化振幅を一致させれば、極軸の片側では二本が同じ方向、他方では逆方向を向く。無次元の素モーメントを選ぶと

$$
\begin{aligned}
\vb d&=\left(\begin{matrix}1 & 0 & 0\end{matrix}\right),\quad \vb h=\left(\begin{matrix}0 & 1 & 0\end{matrix}\right),\\
W&=\left(\begin{matrix}0 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{matrix}\right).
\end{aligned}
$$

極軸での場をまず調べる。$\vb n=\pm\hat{\vb z}$ では横射影は第三成分を落とす演算であり、電気双極子の項は方向反転で変わらず、残りの二項は符号を変える。

$$
\begin{aligned}
\vb f(+\hat{\vb z})&=\left(\begin{matrix}2 & 0 & 0\end{matrix}\right),\\
\vb f(-\hat{\vb z})&=\left(\begin{matrix}0 & 0 & 0\end{matrix}\right).
\end{aligned}
$$

半球全体を見るため、$u=\cos\theta$ と置き $X=\sqrt{1-u^2}\cos\varphi$、$Y=\sqrt{1-u^2}\sin\varphi$ を代入する。$X,Y$ の奇数乗を一つでも含む項は方位角平均で零となる。残る項をまとめて $F(u)=\langle|\vb f|^2\rangle_\varphi$ と書くと

$$
F(u)=\left(u + 1\right)^{2}.
$$

例えば二次の方位角平均は $\langle X^2\rangle=\langle Y^2\rangle=(1-u^2)/2$ である。立体角は $\dd\varphi\,\dd u$ なので、正負の半球を別々に積分すると

$$
\begin{aligned}
P_+/(2\pi P_*)&=\int_0^1 F(u)\dd u=\frac{7}{3},\\
P_-/(2\pi P_*)&=\int_{-1}^0 F(u)\dd u=\frac{1}{3},\\
\langle \cos\theta\rangle&=\frac{\int_{-1}^1 uF(u)\dd u}{\frac{8}{3}}=\frac{1}{2}.
\end{aligned}
$$

各極角で方位角平均した強度を物理単位に戻すと次の値になる。

| 極角 | $F$ | 単位立体角あたりの平均放射電力 $[\mathrm{\mu W/sr}]$ |
|---:|---:|---:|
| $0^\circ$ | $4$ | $8.274$ |
| $45^\circ$ | $\frac{\left(\sqrt{2} + 2\right)^{2}}{4}$ | $6.028$ |
| $90^\circ$ | $1$ | $2.069$ |
| $135^\circ$ | $\frac{\left(-2 + \sqrt{2}\right)^{2}}{4}$ | $0.1775$ |
| $180^\circ$ | $0$ | $0$ |

全電力は $34.66\,\mathrm{\mu W}$ である。後方軸上の一点で放射が消えても、後半球全体の放射が消えるわけではない。有限な開口の検出器に届く量は角度積分から求める。
:::

:::example[磁気双極子の駆動位相を反転する]
前例の磁気双極子だけを半周期反転する。入力の各モーメントの絶対値を変えずに放射の向きと全電力を別々に制御できるかを調べる。無次元の素モーメントを選ぶと

$$
\begin{aligned}
\vb d&=\left(\begin{matrix}1 & 0 & 0\end{matrix}\right),\quad \vb h=\left(\begin{matrix}0 & -1 & 0\end{matrix}\right),\\
W&=\left(\begin{matrix}0 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{matrix}\right).
\end{aligned}
$$

極軸での場をまず調べる。$\vb n=\pm\hat{\vb z}$ では横射影は第三成分を落とす演算であり、電気双極子の項は方向反転で変わらず、残りの二項は符号を変える。

$$
\begin{aligned}
\vb f(+\hat{\vb z})&=\left(\begin{matrix}0 & 0 & 0\end{matrix}\right),\\
\vb f(-\hat{\vb z})&=\left(\begin{matrix}2 & 0 & 0\end{matrix}\right).
\end{aligned}
$$

半球全体を見るため、$u=\cos\theta$ と置き $X=\sqrt{1-u^2}\cos\varphi$、$Y=\sqrt{1-u^2}\sin\varphi$ を代入する。$X,Y$ の奇数乗を一つでも含む項は方位角平均で零となる。残る項をまとめて $F(u)=\langle|\vb f|^2\rangle_\varphi$ と書くと

$$
F(u)=\left(u - 1\right)^{2}.
$$

例えば二次の方位角平均は $\langle X^2\rangle=\langle Y^2\rangle=(1-u^2)/2$ である。立体角は $\dd\varphi\,\dd u$ なので、正負の半球を別々に積分すると

$$
\begin{aligned}
P_+/(2\pi P_*)&=\int_0^1 F(u)\dd u=\frac{1}{3},\\
P_-/(2\pi P_*)&=\int_{-1}^0 F(u)\dd u=\frac{7}{3},\\
\langle \cos\theta\rangle&=\frac{\int_{-1}^1 uF(u)\dd u}{\frac{8}{3}}=- \frac{1}{2}.
\end{aligned}
$$

各極角で方位角平均した強度を物理単位に戻すと次の値になる。

| 極角 | $F$ | 単位立体角あたりの平均放射電力 $[\mathrm{\mu W/sr}]$ |
|---:|---:|---:|
| $0^\circ$ | $0$ | $0$ |
| $45^\circ$ | $\frac{\left(-2 + \sqrt{2}\right)^{2}}{4}$ | $0.1775$ |
| $90^\circ$ | $1$ | $2.069$ |
| $135^\circ$ | $\frac{\left(\sqrt{2} + 2\right)^{2}}{4}$ | $6.028$ |
| $180^\circ$ | $4$ | $8.274$ |

全電力は $34.66\,\mathrm{\mu W}$ である。全電力の値が同じでも縦方向へ運ぶ運動量は符号を変える。遠方パターンを積分する際に方向余弦を掛け忘れるとこの反転を検出できない。
:::

:::example[四重極の軸と赤道の暗線]
対角な電気四重極だけを駆動する。単に行列を方向ベクトルに掛けた結果には縦成分が含まれる。放射電場に寄与する横射影を取る前後で強度分布がどう変わるかを調べる。無次元の素モーメントを選ぶと

$$
\begin{aligned}
\vb d&=\left(\begin{matrix}0 & 0 & 0\end{matrix}\right),\quad \vb h=\left(\begin{matrix}0 & 0 & 0\end{matrix}\right),\\
W&=\left(\begin{matrix}- \frac{1}{2} & 0 & 0\\0 & - \frac{1}{2} & 0\\0 & 0 & 1\end{matrix}\right).
\end{aligned}
$$

極軸での場をまず調べる。$\vb n=\pm\hat{\vb z}$ では横射影は第三成分を落とす演算であり、電気双極子の項は方向反転で変わらず、残りの二項は符号を変える。

$$
\begin{aligned}
\vb f(+\hat{\vb z})&=\left(\begin{matrix}0 & 0 & 0\end{matrix}\right),\\
\vb f(-\hat{\vb z})&=\left(\begin{matrix}0 & 0 & 0\end{matrix}\right).
\end{aligned}
$$

半球全体を見るため、$u=\cos\theta$ と置き $X=\sqrt{1-u^2}\cos\varphi$、$Y=\sqrt{1-u^2}\sin\varphi$ を代入する。$X,Y$ の奇数乗を一つでも含む項は方位角平均で零となる。残る項をまとめて $F(u)=\langle|\vb f|^2\rangle_\varphi$ と書くと

$$
F(u)=- \frac{9 u^{2} \left(u - 1\right) \left(u + 1\right)}{4}.
$$

例えば二次の方位角平均は $\langle X^2\rangle=\langle Y^2\rangle=(1-u^2)/2$ である。立体角は $\dd\varphi\,\dd u$ なので、正負の半球を別々に積分すると

$$
\begin{aligned}
P_+/(2\pi P_*)&=\int_0^1 F(u)\dd u=\frac{3}{10},\\
P_-/(2\pi P_*)&=\int_{-1}^0 F(u)\dd u=\frac{3}{10},\\
\langle \cos\theta\rangle&=\frac{\int_{-1}^1 uF(u)\dd u}{\frac{3}{5}}=0.
\end{aligned}
$$

各極角で方位角平均した強度を物理単位に戻すと次の値になる。

| 極角 | $F$ | 単位立体角あたりの平均放射電力 $[\mathrm{\mu W/sr}]$ |
|---:|---:|---:|
| $0^\circ$ | $0$ | $0$ |
| $45^\circ$ | $- \frac{9 \left(-2 + \sqrt{2}\right) \left(\sqrt{2} + 2\right)}{32}$ | $1.164$ |
| $90^\circ$ | $0$ | $0$ |
| $135^\circ$ | $- \frac{9 \left(-2 + \sqrt{2}\right) \left(\sqrt{2} + 2\right)}{32}$ | $1.164$ |
| $180^\circ$ | $0$ | $0$ |

全電力は $7.798\,\mathrm{\mu W}$ である。節の間には環状の極大が現れる。軸上と赤道上の読みが両方零でも、全放射電力は有限である。
:::

:::example[電気双極子と四重極の干渉]
軸対称な二種類の電気多重極を同時に駆動する。物理的な四重極に含まれる余分な時間微分の位相を補償して、遠方振幅の規格化係数を実数にそろえる。無次元の素モーメントを選ぶと

$$
\begin{aligned}
\vb d&=\left(\begin{matrix}0 & 0 & 1\end{matrix}\right),\quad \vb h=\left(\begin{matrix}0 & 0 & 0\end{matrix}\right),\\
W&=\left(\begin{matrix}- \frac{1}{2} & 0 & 0\\0 & - \frac{1}{2} & 0\\0 & 0 & 1\end{matrix}\right).
\end{aligned}
$$

極軸での場をまず調べる。$\vb n=\pm\hat{\vb z}$ では横射影は第三成分を落とす演算であり、電気双極子の項は方向反転で変わらず、残りの二項は符号を変える。

$$
\begin{aligned}
\vb f(+\hat{\vb z})&=\left(\begin{matrix}0 & 0 & 0\end{matrix}\right),\\
\vb f(-\hat{\vb z})&=\left(\begin{matrix}0 & 0 & 0\end{matrix}\right).
\end{aligned}
$$

半球全体を見るため、$u=\cos\theta$ と置き $X=\sqrt{1-u^2}\cos\varphi$、$Y=\sqrt{1-u^2}\sin\varphi$ を代入する。$X,Y$ の奇数乗を一つでも含む項は方位角平均で零となる。残る項をまとめて $F(u)=\langle|\vb f|^2\rangle_\varphi$ と書くと

$$
F(u)=- \frac{\left(u - 1\right) \left(u + 1\right) \left(3 u + 2\right)^{2}}{4}.
$$

例えば二次の方位角平均は $\langle X^2\rangle=\langle Y^2\rangle=(1-u^2)/2$ である。立体角は $\dd\varphi\,\dd u$ なので、正負の半球を別々に積分すると

$$
\begin{aligned}
P_+/(2\pi P_*)&=\int_0^1 F(u)\dd u=\frac{103}{60},\\
P_-/(2\pi P_*)&=\int_{-1}^0 F(u)\dd u=\frac{13}{60},\\
\langle \cos\theta\rangle&=\frac{\int_{-1}^1 uF(u)\dd u}{\frac{29}{15}}=\frac{12}{29}.
\end{aligned}
$$

各極角で方位角平均した強度を物理単位に戻すと次の値になる。

| 極角 | $F$ | 単位立体角あたりの平均放射電力 $[\mathrm{\mu W/sr}]$ |
|---:|---:|---:|
| $0^\circ$ | $0$ | $0$ |
| $45^\circ$ | $- \frac{\left(-2 + \sqrt{2}\right) \left(4 + 3 \sqrt{2}\right)^{2} \left(\sqrt{2} + 2\right)}{64}$ | $4.392$ |
| $90^\circ$ | $1$ | $2.069$ |
| $135^\circ$ | $- \frac{\left(-4 + 3 \sqrt{2}\right)^{2} \left(-2 + \sqrt{2}\right) \left(\sqrt{2} + 2\right)}{64}$ | $0.003806$ |
| $180^\circ$ | $0$ | $0$ |

全電力は $25.13\,\mathrm{\mu W}$ である。奇関数である交差項は全球で消える一方、両半球の差では二倍になって残る。全球の電力だけを測ると二つの源の相対位相の情報を失う。
:::

:::example[電気と磁気の位相が直角の場合]
二つの双極子が同じ軸配置にあっても相対位相を直角にする。強度に現れる干渉は振幅積の実部なので、同位相のときの方向性が残るかを調べる。無次元の素モーメントを選ぶと

$$
\begin{aligned}
\vb d&=\left(\begin{matrix}1 & 0 & 0\end{matrix}\right),\quad \vb h=\left(\begin{matrix}0 & i & 0\end{matrix}\right),\\
W&=\left(\begin{matrix}0 & 0 & 0\\0 & 0 & 0\\0 & 0 & 0\end{matrix}\right).
\end{aligned}
$$

極軸での場をまず調べる。$\vb n=\pm\hat{\vb z}$ では横射影は第三成分を落とす演算であり、電気双極子の項は方向反転で変わらず、残りの二項は符号を変える。

$$
\begin{aligned}
\vb f(+\hat{\vb z})&=\left(\begin{matrix}1 + i & 0 & 0\end{matrix}\right),\\
\vb f(-\hat{\vb z})&=\left(\begin{matrix}1 - i & 0 & 0\end{matrix}\right).
\end{aligned}
$$

半球全体を見るため、$u=\cos\theta$ と置き $X=\sqrt{1-u^2}\cos\varphi$、$Y=\sqrt{1-u^2}\sin\varphi$ を代入する。$X,Y$ の奇数乗を一つでも含む項は方位角平均で零となる。残る項をまとめて $F(u)=\langle|\vb f|^2\rangle_\varphi$ と書くと

$$
F(u)=u^{2} + 1.
$$

例えば二次の方位角平均は $\langle X^2\rangle=\langle Y^2\rangle=(1-u^2)/2$ である。立体角は $\dd\varphi\,\dd u$ なので、正負の半球を別々に積分すると

$$
\begin{aligned}
P_+/(2\pi P_*)&=\int_0^1 F(u)\dd u=\frac{4}{3},\\
P_-/(2\pi P_*)&=\int_{-1}^0 F(u)\dd u=\frac{4}{3},\\
\langle \cos\theta\rangle&=\frac{\int_{-1}^1 uF(u)\dd u}{\frac{8}{3}}=0.
\end{aligned}
$$

各極角で方位角平均した強度を物理単位に戻すと次の値になる。

| 極角 | $F$ | 単位立体角あたりの平均放射電力 $[\mathrm{\mu W/sr}]$ |
|---:|---:|---:|
| $0^\circ$ | $2$ | $4.137$ |
| $45^\circ$ | $\frac{3}{2}$ | $3.103$ |
| $90^\circ$ | $1$ | $2.069$ |
| $135^\circ$ | $\frac{3}{2}$ | $3.103$ |
| $180^\circ$ | $2$ | $4.137$ |

全電力は $34.66\,\mathrm{\mu W}$ である。この場合、二種類の双極子は強度では干渉しない。双方が存在することは偏光を分解した振幅や位相から調べることができる。
:::


同位相の双極子対で得た角分布を極座標に描くと、後方の暗点の周囲にも弱い放射が広がることが見える。図の半径は電場振幅でなく単位立体角あたりの電力に比例する。

```tikz caption="電気・磁気双極子が同位相で干渉したときの角分布。前方を右に取り、半径を最大強度で規格化した。"
\begin{tikzpicture}
\draw[->] (-.8,0)--(3.4,0) node[right] {$+z$};
\draw[dashed] (0,-1.4)--(0,1.4);
\draw[thick,domain=0:360,samples=81,variable=\t] plot ({.7*(1+cos(\t))*(1+cos(\t))*cos(\t)},{.7*(1+cos(\t))*(1+cos(\t))*sin(\t)});
\fill (0,0) circle (1.5pt);
\end{tikzpicture}
```

## Larmor の公式

### 遅延場から粒子の放射電力へ

単一点電荷が非相対論的に運動すると $\ddot{\vb{p}}=q\vb{a}$ である。
従って :eqref[eq-em-dipole-radiation] は :term[Larmor の公式]{en="Larmor formula" yomi="らーもあのこうしき"}

$$
P=\frac{q^2a^2}{6\pi\varepsilon_0 c^3}
\label{eq-em-larmor}
$$

となる。
力が速度に垂直でも加速度があれば放射する。
仕事と放射の瞬時関係には、粒子付近の場に蓄えられるエネルギーも関わるため、機械エネルギーの変化を常に $-P$ と等置することはできない。

相対論的な粒子では :eqref[eq-em-lw-field] の加速度場から遠方の強度を求める。
観測時間あたりの電力 $R^2\varepsilon_0cE_{\rm rad}^2$ を放出時間あたりに直すため $\dd t/\dd t_*=\kappa$ を掛ける。

従って

$$
\frac{\dd P}{\dd\Omega}=\frac{q^2}{16\pi^2\varepsilon_0c}\frac{|\vb{n}\times[(\vb{n}-\vb{\beta})\times\dot{\vb{\beta}}]|^2}{(1-\vb{n}\cdot\vb{\beta})^5}.
$$

全電力を求めるには、速度を $z$ 方向に取り、$u=\cos\theta$、$\kappa=1-\beta u$ と置く。加速度を速度に平行な成分と垂直な成分に分け、垂直成分の向きを $x$ 軸に取る。上式の分子の交差項は $\cos\varphi$ に比例するので、方位角について平均すると消える。二重外積を展開し、$\langle n_x^2\rangle_\varphi=(1-u^2)/2$ を使うと、平行成分の係数は $1-u^2$、垂直成分の係数は $\kappa^2-(1-\beta^2)(1-u^2)/2$ となる。

従って必要な積分は二つである。$w=1-\beta u$ と変数を換え、$w=1-\beta$ から $1+\beta$ までのべき積分を計算すると

$$
\int_{-1}^1\frac{\dd u}{\kappa^3}=\frac{2}{(1-\beta^2)^2},\qquad
\int_{-1}^1\frac{1-u^2}{\kappa^5}\,\dd u=\frac{4}{3(1-\beta^2)^3}.
$$

$\beta=0$ でも連続に同じ値を取る。方位角積分の $2\pi$ と $\dot{\vb\beta}=\vb a/c$ を戻すと、$a_\parallel^2$ の係数には $\gamma^6$、$a_\perp^2$ の係数には $\gamma^4$ が掛かる。ここで $\gamma=(1-\beta^2)^{-1/2}$ である。

従って Liénard の公式は

$$
P=\frac{q^2}{6\pi\varepsilon_0c^3}\gamma^6(a^2-|\vb{\beta}\times\vb{a}|^2).
\label{eq-em-lienard-power}
$$

等速円運動では $\vb{a}\perp\vb{v}$ なので係数は $\gamma^4$、縦加速なら $\gamma^6$ である。
$\gamma\gg1$ で $\kappa\simeq(\gamma^{-2}+\theta^2)/2$ となる方向では放射が小角度に集中する。

実際のシンクロトロン放射の周波数スペクトルには、観測する円弧全体の位相積分が必要となる。
:cite[jackson-em]

:::example[一電子の放射]
低速電子が局所的に $a=1.00\times10^{16}\,\mathrm{m/s^2}$ で加速すると、$e=1.6022\times10^{-19}\,\mathrm C$ を代入して

$$
P=5.71\times10^{-22}\,\mathrm W
$$

となる。

二個の電子の放射電力を足すときは、互いの位相がそろうかどうかを確かめる。
完全に同相の同じ振幅なら電場が倍となり電力は四倍である。
:::

:::example[垂直加速の円軌道とエネルギー収支]
一様磁場 $B$ に電子が垂直に入射する。
磁気力の運動方程式から瞬間的に $a=e v B/(\gamma m_e)$ である。
これを横加速の Liénard の式に入れると

$$
P=\frac{e^4B^2}{6\pi\varepsilon_0m_e^2c}\gamma^2\beta^2.
$$

$\gamma=10,\ B=0.100\,\mathrm T$ なら $P=1.57\times10^{-14}\,\mathrm W$ 程度である。
磁気力自身は仕事をしないので、長い時間にわたって同じエネルギーで周回させるには、この電力を別の電場から補う必要がある。

短い区間の軌道を放射なしで計算できるかは、一周で失うエネルギーを運動エネルギーと比較して判定する。
:::

### 角度分布を積分して縦加速の倍率を確かめる

速度と加速度をともに $z$ 軸へ向けると三重積の二乗は $\dot\beta^2(1-u^2)$、$u=\cos\theta$ になる。
放出時間基準の角度積分で必要なのは

$$
\int_{-1}^1\frac{1-u^2}{(1-\beta u)^5}\dd u=\frac{4}{3(1-\beta^2)^3}.
$$

これを得るには $s=1-\beta u$ と置き、$\beta^{-3}\int_{1-\beta}^{1+\beta}[(\beta^2-1)s^{-5}+2s^{-4}-s^{-3}]\dd s$ を積分する。
両端を代入して整理すれば右辺である。

方位角の $2\pi$ を掛けて $a=c\dot\beta$ を使うと、Larmor の値の $\gamma^6$ 倍が得られる。
観測時刻と放出時刻を取り違えて分母を六乗のまま積分するとこの一致が失われる。

### 放射反作用と近接場のエネルギー

非相対論的な点粒子の形式的な Abraham–Lorentz 力を

$$
\vb{F}_{\rm self}=q^2\dot{\vb{a}}/(6\pi\varepsilon_0c^3)
$$

とする。
これは有限サイズの電荷の自己場を極限に取り、電磁的な慣性の発散を質量に吸収して得る式である。

点極限の正当化や物理的な解の選択には自己力の理論が必要なため、ここではこの式の導出まで立ち入らない。
:cite[jackson-em]

この力の仕事が放射と矛盾しないことは計算できる。
$\vb{v}\cdot\dot{\vb{a}}=\dd(\vb{v}\cdot\vb{a})/\dd t-a^2$ なので、$C_q=q^2/(6\pi\varepsilon_0c^3)$ とおけば $P_{\rm ext}=\dot K+P_{\rm rad}+\dot E_{\rm near}$、

$$
E_{\rm near}=-C_q\vb{v}\cdot\vb{a}
$$

となる。

加速のない二時刻の間、または一周期の積分では近接場の全微分が消え、外力の仕事と放射・運動エネルギーの増分とが一致する。

運動方程式を $m(\vb a-\tau_q\dot{\vb a})=\vb F_{\rm ext}$ と書くと、$\tau_q=C_q/m$ は時間の次元を持つ。外力が零でも $\vb a=\vb a_0 e^{t/\tau_q}$ が形式解である。この暴走解があるため、任意の初期加速度を指定して点粒子の式を厳密に時間発展させることには問題がある。外力が時間尺度 $T\gg\tau_q$ で滑らかに変わるとき、$\tau_q$ の一次までなら自己力の中で $\dot{\vb a}$ を最低次の $\dot{\vb F}_{\rm ext}/m$ に置き換えられる。その結果 $m\vb a=\vb F_{\rm ext}+\tau_q\dot{\vb F}_{\rm ext}$ を得る。ここで外力の全微分は最低次の軌道に沿って評価する。

### 力の向きと放射電力

角度積分から得た :eqref[eq-em-lienard-power] を加速度の成分に分けて書くと

$$
P=\frac{q^2}{6\pi\varepsilon_0 c^3}
(\gamma^6a_\parallel^2+\gamma^4 a_\perp^2).
$$

一定の力で比較するため、質量 $m$ の粒子に対し運動量 $\vb p=m\gamma\vb v$ と力 $\vb F=\dot{\vb p}$ を用いる。$\gamma$ を微分すると $\dot\gamma=\gamma^3\vb v\cdot\vb a/c^2$ である。$\vb F=\dd(\gamma m\vb v)/\dd t$ から $F_\parallel=m\gamma^3a_\parallel$、$F_\perp=m\gamma a_\perp$ なので、$P=q^2(F_\parallel^2+\gamma^2F_\perp^2)/(6\pi\varepsilon_0m^2 c^3)$ である。同じ加速度で比べる実験と同じ力で比べる実験を取り違えると、速度依存の指数を誤る。


:::example[円運動の電子が一周で失うエネルギー]
半径 $b$ の円運動では $a=v^2/b$ が速度に垂直である。Liénard 電力に周回時間 $2\pi b/v$ を掛けると一周の損失は $q^2\gamma^4\beta^3/(3\varepsilon_0 b)$ である。軌道と速度が一周内でほぼ一定であるとする。

$x=\beta$、電子の $b=0.100\,\mathrm m$ に対して $Y_*=e^2/(3\varepsilon_0 b)=6.03\times10^{-8}\,\mathrm{eV}$ である。$Y$ は一周の放射エネルギーである。

$$
f(x)=\frac{x^{3}}{\left(x - 1\right)^{2} \left(x + 1\right)^{2}},\qquad Y=Y_* f(x).
$$

変数を少し動かしたときの読みの変化は、接線の傾きから計算できる。基準値 $x_0$ の近くで

$$
\begin{aligned}
\frac{\Delta Y}{Y_*}= f\!\left(x_0+\Delta x\right)-f(x_0)= f'(x_0)\Delta x+O((\Delta x)^2),\\
f'(x)=- \frac{x^{2} \left(x^{2} + 3\right)}{\left(x - 1\right)^{3} \left(x + 1\right)^{3}}.
\end{aligned}
$$

同じ尺度で値と傾きを並べる。

| $x$ | $f(x)$ | $Y$ $[\mathrm{eV}]$ | $\dd Y/\dd x$ $[\mathrm{eV}]$ |
|---:|---:|---:|---:|
| $\frac{1}{5}$ | $0.008681$ | $5.236\times10^{-10}$ | $8.290\times10^{-9}$ |
| $\frac{3}{5}$ | $0.5273$ | $3.181\times10^{-8}$ | $2.783\times10^{-7}$ |
| $\frac{9}{10}$ | $20.19$ | $1.218\times10^{-6}$ | $2.714\times10^{-5}$ |
| $\frac{99}{100}$ | $2450$ | $1.478\times10^{-4}$ | $0.02986$ |

同じ半径を保つために必要な磁場は運動量とともに増える。この表は磁場を一定にして速度を変えた場合とは条件が異なる。実際には放射損失を無視できない周回数で加速または軌道修正が必要となる。
:::

## 散乱

### 散乱断面積を測る

入射平面波の平均強度を $I_{\rm in}$、散乱波が単位立体角へ運ぶ平均電力を $\dd P_s/\dd\Omega$ とすると、:term[微分散乱断面積]{en="differential scattering cross section" yomi="びぶんさんらんだんめんせき"} は $(\dd P_s/\dd\Omega)/I_{\rm in}$ である。

単位は $\mathrm{m^2/sr}$ であり、全立体角で積分すると全散乱断面積 $\sigma_s$ になる。

自由な低速電子では $m_e\vb{a}=-e\vb{E}_{\rm in}$ と近似する。
入射場をピーク $E_0$、偏光単位ベクトルを実の $\vb{e}$ とし、散乱後の偏光を区別しない。
双極子の角度分布を入射強度 $\varepsilon_0 c E_0^2/2$ で割ると :term[Thomson 散乱]{en="Thomson scattering" yomi="とむそんさんらん"} の式

$$
\frac{\dd\sigma_T}{\dd\Omega}=r_e^2[1-(\vb{n}\cdot\vb{e})^2],\quad r_e=\frac{e^2}{4\pi\varepsilon_0m_ec^2},\quad \sigma_T=\frac{8\pi}{3}r_e^2
\label{eq-em-thomson}
$$

が得られる。
非偏光入射では横二偏光について平均するため、入射方向との角を $\theta$ として

$$
\dd\sigma_T/\dd\Omega=r_e^2(1+\cos^2\theta)/2
$$

となる。
$r_e=2.818\times10^{-15}\,\mathrm m$、

$$
\sigma_T=6.65\times10^{-29}\,\mathrm{m^2}
$$

である。

電子の反跳や量子効果が重要になる光子エネルギーでは、この古典的導出の仮定が失われる。

### 束縛電子と小さな球

一様な入射電場で誘起される双極子を

$$
\vb{p}_0=\alpha(\omega)\vb{E}_0
$$

とする。
双極子放射の時間平均を強度で割れば

$$
\sigma_s=\frac{k^4}{6\pi\varepsilon_0^2}|\alpha|^2,\qquad k=\omega/c.
\label{eq-em-rayleigh}
$$

Lorentz 振動子の $\alpha=(e^2/m)/(\omega_0^2-\omega^2-i\gamma\omega)$ を用いると、自由電子の $\sigma_T$ に $\omega^4/[(\omega_0^2-\omega^2)^2+\gamma^2\omega^2]$ を掛けた形になる。

$\omega\ll\omega_0$ では $\sigma_s\propto\lambda^{-4}$ であり、これを :term[Rayleigh 散乱]{en="Rayleigh scattering" yomi="れいりーさんらん"} と呼ぶ。

逆に束縛と摩擦が無視できる高周波側では Thomson 値に近づく。

半径 $a$、比誘電率 $\varepsilon_r$ の非磁性球で、外部・内部の波長がともに球より十分長ければ、[誘電体](/physics/electromagnetism/dielectrics/) の静電解 $\alpha=4\pi\varepsilon_0a^3(\varepsilon_r-1)/(\varepsilon_r+2)$ を使える。

よって

$$
\sigma_s=(8\pi/3)k^4a^6|(\varepsilon_r-1)/(\varepsilon_r+2)|^2
$$

である。
波長に比べて大きい粒子の解には球面の全多重極で Maxwell 境界条件を接続する必要があるため、この近似の適用外である。
:cite[jackson-em]

:::example[二種類の波長]
分極率の波長変化が小さく、どちらの波長にも Rayleigh 近似が成り立つとする。
$450\,\mathrm{nm}$ と $650\,\mathrm{nm}$ の散乱断面積の比は $(650/450)^4=4.35$ である。
入射強度や検出器の感度が異なる場合は、観測した色の強度比を直接この比とみなせない。
:::

:::example[散乱と吸収のエネルギー収支]
振動双極子に入射場がする平均仕事は $\langle\dot{\vb{p}}\cdot\vb{E}_{\rm in}\rangle=\omega\Im\alpha\,E_0^2/2$ なので、消衰断面積は

$$
\sigma_{\rm ext}=k\Im\alpha/\varepsilon_0
$$

となる。

吸収は

$$
\sigma_{\rm abs}=\sigma_{\rm ext}-\sigma_s
$$

である。
損失のない球に実の静的 $\alpha$ を使うと、散乱は計算できても同じ精度で消衰を求められない。
放射によって必要な虚部は

$$
\Im\alpha=k^3|\alpha|^2/(6\pi\varepsilon_0)
$$

であり、$1/\alpha=1/\alpha_{\rm stat}-ik^3/(6\pi\varepsilon_0)$ とすればこの収支を満たす。
:::

### 偏光を選んで観測する散乱

散乱方向に直交する単位ベクトル $\vb{e}_s$ の偏光だけを計測する場合、横射影された双極子場をさらに $\vb{e}_s$ に射影する。
真空中の自由電子ではその微分断面積は $r_e^2|\vb{e}_s^*\cdot\vb{e}|^2$ である。

散乱方向に直交する二つの基底で足せば $r_e^2[1-|\vb{n}\cdot\vb{e}|^2]$ に戻る。
これは偏光検出器の向きによって同じ散乱方向でも受信電力が変わることを表す。

:::example[直角方向で生じる直線偏光]
非偏光が $+z$ へ進み、$+x$ 方向の散乱を観測するとする。
入射 $x$ 偏光では電子加速度が観測方向と平行で、この方向には放射がない。
一方、入射 $y$ 偏光では $y$ 偏光の散乱が生じる。
従って直角方向の散乱光は完全な直線偏光である。

一般角では、散乱面に垂直な成分と平行な成分に対する強度が $1:\cos^2\theta$ となる。
これらから直線偏光度 $(1-\cos^2\theta)/(1+\cos^2\theta)$ が得られる。
:::

### 散乱断面積と観測した減衰率を対応させる

独立に散乱する同種の粒子が数密度 $n_s$ で入った薄い層を考える。
面積 $A$、厚さ $\dd z$ の層の粒子数は $n_sA\dd z$ である。
入射方向のビームから各粒子が奪う電力は $\sigma_{\rm ext}I$ なので、層を通した後の直進強度は

$$
\dd I=-n_s\sigma_{\rm ext}I\dd z
$$

に従う。

従って

$$
I(z)=I(0)e^{-n_s\sigma_{\rm ext}z}.
$$

例えば $n_s=10^{18}\,\mathrm{m^{-3}}$、$\sigma_{\rm ext}=2.0\times10^{-16}\,\mathrm{m^2}$ なら強度の減衰長は $5.0\,\mathrm{mm}$ である。
ここで数えているのは入射方向のコヒーレントなビームから失われる電力であり、失われた電力のすべてが熱になるわけではない。

周囲から散乱光が再び検出器に入る場合には、この直進成分だけで総検出強度を記述できない。

### 前方散乱だけで消衰を測れる理由

直進光が失う電力は、全方向に散乱した電力と物体内で吸収した電力の和である。一方、その和は前方の散乱振幅を一つ測るだけで分かる。前方だけが特別になる理由を、入射場が物体にする仕事から導く。

真空中の局在した物体に、単位偏光 $\vb{e}$、進行方向 $\hat{\vb{z}}$ の平面波 $\vb{E}_{\rm in}=E_0\vb{e}e^{ikz}$ を入れる。時間因子は $e^{-i\omega t}$、$E_0$ は正のピーク振幅とする。入射によって物体に生じた全電流を $\vb{J}$ とし、束縛電荷の動きもこれに含める。散乱遠方場を $\vb{E}_s=E_0 e^{ikr}\vb{f}(\vb{n})/r+O(r^{-2})$ と書けば、遅延ベクトルポテンシャルの遠方展開と電場の横射影から

$$
\vb{f}(\vb{n})=\frac{i\omega\mu_0}{4\pi E_0}
\left[\int \vb{J}(\vb{r}')e^{-ik\vb{n}\cdot\vb{r}'}\,\dd^3r'\right]_\perp
$$

となる。$\vb{f}$ は長さの単位をもち、$\dd\sigma_s/\dd\Omega=|\vb{f}|^2$ である。この式は源の大きさを波長より小さいとは仮定していない。

全場が電流にする仕事は入射場の仕事と散乱場の仕事に分かれる。周期平均の Poynting の定理を散乱場だけに適用すると、散乱場のする仕事は外へ出る散乱電力の負値である。一方、定常状態で全場が物体に与えた正味の仕事は吸収電力に等しい。従って入射場のする仕事 $P_{\rm ext}$ が散乱電力と吸収電力の和になる。

$$
P_{\rm ext}=\tfrac12\Re\int\vb{E}_{\rm in}^*\cdot\vb{J}\,\dd^3r'
=\frac{E_0}{2}\Re\int\vb{e}^*\cdot\vb{J}\,e^{-ikz'}\,\dd^3r'.
$$

最後の積分に現れる位相は、前方 $\vb{n}=\hat{\vb{z}}$ の散乱振幅に現れる位相と等しい。$\vb{e}$ はこの方向に垂直なので横射影も内積に影響しない。入射強度 $I=E_0^2/(2Z_0)$ と $Z_0=\omega\mu_0/k$ を使えば、:term[光学定理]{en="optical theorem" yomi="こうがくていり"} は

$$
\sigma_{\rm ext}=\frac{P_{\rm ext}}{I}
=\frac{4\pi}{k}\Im\left[\vb{e}^*\cdot\vb{f}(\hat{\vb{z}})\right]
$$

となる。前方の強度 $|\vb{f}|^2$ だけでは虚部を決められない。実験でこの関係を利用できるのは、散乱場と入射場の干渉が相対位相を反映するためである。

:::example[一つの電気双極子が散乱できる断面積の上限]
等方的な電気双極子だけが励起される場合、$u=k^3\alpha/(6\pi\varepsilon_0)$ と置く。双極子放射と入射場の仕事から $\sigma_s=(6\pi/k^2)|u|^2$、$\sigma_{\rm ext}=(6\pi/k^2)\Im u$ である。受動的な散乱体では吸収が非負なので、$u=X+iY$ の取りうる範囲は $X^2+(Y-1/2)^2\le1/4$ となる。この円板で $|u|^2$ が最大になるのは $u=i$ であり、最大散乱断面積は $6\pi/k^2$ である。一方、吸収に比例する $Y-|u|^2=1/4-X^2-(Y-1/2)^2$ は $u=i/2$ で最大になる。

真空波長 $600\,\mathrm{nm}$ の場合、最大散乱断面積は $1.72\times10^{-13}\,\mathrm{m^2}$、最大吸収断面積はその四分の一の $4.30\times10^{-14}\,\mathrm{m^2}$ である。前者で吸収は零、後者では散乱と吸収が等しい。別の磁気応答や高次の多重極も励起されると、全断面積にはそれらの寄与も加わる。
:::

### 散乱による運動量移行と放射圧断面積

入射方向を $z$ とし、弾性散乱と吸収を許す静止した小粒子を考える。方向 $\theta$ に散乱された光が持ち去る縦運動量は、同じエネルギーの入射光の $\cos\theta$ 倍である。運動量の収支を全方向でとると、縦の力 $F_z=(I_i/c)\sigma_{\rm pr}$ に対する断面積は

$$
\sigma_{\rm pr}=\sigma_{\rm abs}+\int(1-\cos\theta)
\frac{\dd\sigma}{\dd\Omega}\dd\Omega
=\sigma_{\rm ext}-g\sigma_{\rm sca},\qquad
g=\frac{1}{\sigma_{\rm sca}}\int\cos\theta\frac{\dd\sigma}{\dd\Omega}\dd\Omega .
$$

$g$ は方向余弦の散乱重み付き平均なので $-1\le g\le1$ である。Thomson 散乱のように前後対称なら $g=0$、前方に鋭く集中する場合は $g$ が一に近い。従って入射光束から除かれたエネルギー量が等しい二粒子でも、前方散乱の強い粒子のほうが受ける力は小さい。吸収後の再放射が方向性を持つ場合は、その反動も別途追加する。


:::example[有限開口で集める Thomson 散乱]
非偏光光に対する $\dd\sigma/\dd\Omega=r_e^2(1+u^2)/2$ を、前方の円錐内だけで積分する。円錐の半頂角 $\Theta$ に対して下限 $x=\cos\Theta$、上限一とし、全球での積分で割ると $\int_x^1(1+u^2)\dd u/(8/3)$ となる。

$Y$ は全散乱光に対する収集率、$Y_*=100$ とし百分率で表す。電子一個の全断面積 $6.65\times10^{-29}\,\mathrm{m^2}$ にこの分率を掛ければ有効な収集断面積になる。

$$
f(x)=- \frac{\left(x - 1\right) \left(x^{2} + x + 4\right)}{8},\qquad Y=Y_* f(x).
$$

二点の読みの差を求めるには、基準を同じに保ち、振幅の換算を最後に一度だけ行う。最初と最後の設定で

$$
\Delta Y/Y_* = f\!\left(1\right)-f\!\left(0\right)=-0.5.
$$

第四列は初期設定からの変化である。

| $x$ | $f(x)$ | $Y$ $[\%]$ | $Y-Y(x_0)$ $[\%]$ |
|---:|---:|---:|---:|
| $0$ | $0.5$ | $50$ | $0$ |
| $\frac{1}{4}$ | $0.4043$ | $40.43$ | $-9.57$ |
| $\frac{1}{2}$ | $0.2969$ | $29.69$ | $-20.31$ |
| $\frac{9}{10}$ | $0.07137$ | $7.138$ | $-42.86$ |
| $1$ | $0$ | $0$ | $-50$ |

$90^\circ$ の円錐で半分が集まる。実際の検出器の面積と距離が決まっている場合、面積を距離の二乗で割る近似は小開口でのみ正しく、広い開口ではこのように向きごとの立体角を積分する。
:::

### 一様な球による Mie 散乱と Rayleigh 極限

波長に比べて粒子が大きくなると、内部を一様に分極した球の分極率だけでは散乱を予測できない。真空中に均一な非磁性球を置き、その半径を $a$、屈折率を $m$、外側の波数を $k$ とする。球面に接する場を、電場が半径成分を持たない TE 成分と、磁場が半径成分を持たない TM 成分に分ける。各 $\ell\ge1$ と方位角次数は独立である。放射条件には $e^{-i\omega t}$ に対して外向きの球 Hankel 関数 $h_\ell^{(1)}$ を選ぶ。

各偏光で接線場を作る動径関数を $u$ とする。$\psi_\ell(s)=s j_\ell(s)$、$\xi_\ell(s)=s h_\ell^{(1)}(s)$ は $u''+[1-\ell(\ell+1)/s^2]u=0$ のそれぞれ正則、外向きの解である。方位角次数を $\nu$ と書く。例えば TE 成分では、接線電場を $u(r)\vb X_{\ell\nu}/r$ と書くと、回転から得る接線磁場は $(\partial_r u)/(i\omega\mu r)$ に比例する。よって $u,\mu^{-1}\partial_r u$ を接続する。TM では電場と磁場を交換して $u,\epsilon^{-1}\partial_r u$ を接続する。この二つを区別するため、外側で正則な入射成分の振幅を一とし、外部解を $\psi_\ell(kr)-s_\ell\xi_\ell(kr)$、内部解を $c_\ell\psi_\ell(mkr)$ と表す。

$x=ka$ と置く。電位の鏡像とは異なり、球面で場を零にする条件ではない。両側の接線場を等置すると、偏光に応じて $h=m$（TE）または $h=1/m$（TM）として

$$
\begin{aligned}
\psi_\ell(x)-s_\ell\xi_\ell(x)&=c_\ell\psi_\ell(mx),\\
\psi_\ell'(x)-s_\ell\xi_\ell'(x)&=hc_\ell\psi_\ell'(mx).
\end{aligned}
$$

第二式の $h$ は引数を半径で微分するときの $m$ と構成定数の比から来る。$c_\ell$ を消去すれば

$$
s_\ell=\frac{\psi_\ell(mx)\psi_\ell'(x)-h\psi_\ell'(mx)\psi_\ell(x)}
{\psi_\ell(mx)\xi_\ell'(x)-h\psi_\ell'(mx)\xi_\ell(x)}.
\label{eq-em-mie-channel}
$$

通常、TM の係数を $a_\ell$、TE の係数を $b_\ell$ と呼ぶ。$m=1$ で分子は消え、すべての次数で散乱は起きない。吸収しない実数 $m$ の場合、分母は実部が分子 $N$ に等しく、$N+iD$ と書ける。このため $s_\ell=N/(N+iD)$ は $\Re s_\ell=|s_\ell|^2$ を満たす。この規約で各成分の消衰と散乱が等しいことは、入射の正則解を内向きと外向きに半分ずつ分けるとすぐ確かめられる。外向きの振幅比は $1-2s_\ell$ で、上の等式によりその絶対値は一である。

$\ell=1$ では $\psi_1(s)=s^2/3+O(s^4)$、$\xi_1(s)=-i/s+O(s)$ である。$x\ll1,|m|x\ll1$ とし、TM の $h=1/m$ で代入すると

$$
a_1=-\frac{2i}{3}\frac{m^2-1}{m^2+2}x^3+O(x^5).
$$

これに対し TE では分子の三次項が消え $b_1=O(x^5)$ となる。従って通常の小さな非磁性球では電気双極子が優勢であり、散乱電場の振幅に静電分極率と同じ物性因子が現れる。粒子の外側で長波長でも、高屈折率の球では $|m|x$ が小さいという条件が先に破れる場合がある。

### 平面波に対する各多重極の散乱断面積

:eqref[eq-em-mie-channel] の各係数から、球が平面波から奪う電力を求める。次数による重みは、入射偏光を球面上のどのモードが受け取るかで決まる。入射波をピーク振幅一、方向 $\hat{\vb z}$、単位偏光 $\vb e$ とする。大きな球面上の位相 $e^{ikr\cos\theta}$ を任意の滑らかな角度関数 $g$ に掛けて積分する。$u=\cos\theta$ とし、$g$ を方位角で平均したものを $\bar g(u)$ と書く。部分積分で

$$
2\pi\int_{-1}^1e^{ikru}\bar g(u)\,\dd u
=\frac{2\pi}{ikr}\{e^{ikr}g(\hat{\vb z})-e^{-ikr}g(-\hat{\vb z})\}+O((kr)^{-2})
$$

を得る。従って外向き部分を規格直交な接線関数 $\vb U_{\ell\nu}$ に射影した係数は、$e^{ikr}/r$ を除いて $2\pi\vb U_{\ell\nu}^*(\hat{\vb z})\cdot\vb e/(ik)$ である。この計算では滑らかな関数を通じて射影している。平面波の振幅そのものが球面上の各点で $1/r$ に減衰するという意味ではない。

ここで $\vb U$ は電気型と磁気型のいずれかの基底である。球によって外向き振幅が $1-2s_\ell$ 倍になるので、入射分を引いた散乱振幅の係数はこの射影の $-2s_\ell$ 倍である。各族について、同じ方向での和 $\sum_\nu\vb U_{\ell\nu}\vb U_{\ell\nu}^\dagger$ は接平面への射影に比例する。実際、球面の回転で同じ次数の基底は互いに混ざるだけで、さらに $\nu$ と $-\nu$ を足すと行列は実対称になる。跡は方向によらず、その球面積分は正規直交性により $2\ell+1$ である。接平面が二次元であるため単位偏光に対する和は $(2\ell+1)/(8\pi)$ になる。

直交する各成分の散乱電力を足す場合は係数の絶対値を二乗する。一方、前方で光学定理を使う場合は振幅のまま足す。前方偏光振幅が $i\sum_{\ell}(2\ell+1)(a_\ell+b_\ell)/(2k)$ となるので、それぞれ

$$
\begin{aligned}
\sigma_{\rm sca}&=\frac{2\pi}{k^2}\sum_{\ell=1}^\infty(2\ell+1)(|a_\ell|^2+|b_\ell|^2),\\
\sigma_{\rm ext}&=\frac{2\pi}{k^2}\sum_{\ell=1}^\infty(2\ell+1)\Re(a_\ell+b_\ell),\qquad
\sigma_{\rm abs}=\sigma_{\rm ext}-\sigma_{\rm sca}.
\end{aligned}
\label{eq-em-mie-cross-sections}
$$

無損失なら各次数で二つの和の項が等しくなる。小球の $a_1$ の先頭項を第一式に入れると $\sigma_{\rm sca}=(8\pi/3)k^4a^6|(m^2-1)/(m^2+2)|^2$ となり、静電分極率から求めた Rayleigh 散乱と一致する。実数 $m$ に対するこの $a_1$ の先頭項は純虚数である。消衰を同じ精度で求めるには $\Re a_1=|a_1|^2$ により $x^6$ 次の実部も残す必要がある。例えば $\lambda=600\,\mathrm{nm}$ で一つの電気双極子成分が $a_1=1$ に達すると、断面積は $6\pi/k^2=1.72\times10^{-13}\,\mathrm{m^2}$ となる。係数が実数 $a_1=1/2$ の受動的な粒子では、同じ成分の消衰電力が散乱と吸収に等分される。

## 演習問題

::::exercise{level=basic}
遅延時刻で $\vb{n}\cdot\vb{\beta}=0.60$ のとき、源時刻が $2.0\,\mathrm{ns}$ 進む間の到着時刻の増分を一次近似で求めよ。

:::solution

$$
\dd t=\kappa\dd t_*=(1-0.60)\dd t_*
$$

より $0.80\,\mathrm{ns}$ である。
:::
::::

::::exercise{level=basic}
電気双極子の加速度方向から $30^\circ$ と $90^\circ$ の方向で、単位立体角あたりの電力比を求めよ。

:::solution
角度分布は $\sin^2\theta$ に比例する。
従って

$$
\sin^230^\circ/\sin^290^\circ=1/4
$$

である。
:::
::::

::::exercise{level=basic}
一定の双極子モーメント振幅で周波数を二倍にする。
平均放射電力の倍率を求めよ。

:::solution
二階時間微分の振幅は $\omega^2p_0$ である。
その二乗を平均するため電力は $\omega^4$ に比例し、$2^4=16$ 倍になる。
:::
::::

::::exercise{level=basic}
$|\vb{v}|=0.80c$ の粒子を速度に直交して加速する。
同じ実験室加速度に対する Larmor 値との比を求めよ。

:::solution
横加速なので比は $\gamma^4$ である。
$\gamma=1/\sqrt{1-0.8^2}=5/3$ より、比は $625/81=7.72$ となる。
:::
::::

::::exercise{level=advanced}
同じ放射周波数で $m_0$ の磁気双極子と $p_0$ の電気双極子が単独で同じ平均電力をもつ条件を求めよ。

:::solution
電気双極子の電力には $\omega^4p_0^2/c^3$、磁気双極子には $\omega^4m_0^2/c^5$ が現れ、残る係数は等しい。
従って $m_0=cp_0$ である。
:::
::::

::::exercise{level=advanced}
$Q_{zz}=Q_0\cos\omega t$、$Q_{xx}=Q_{yy}=-Q_{zz}/2$、他の成分が零の電気四重極の平均全放射電力を求めよ。

:::solution
三階微分の二乗平均は $zz$ 成分で $\omega^6Q_0^2/2$、$xx,yy$ はその四分の一ずつである。
和は $3\omega^6Q_0^2/4$ なので

$$
\langle P_Q\rangle=\omega^6Q_0^2/(960\pi\varepsilon_0c^5)
$$

である。
:::
::::

::::exercise{level=advanced}
長さ $\lambda/50$ の一様電流要素の放射抵抗と指向性を求めよ。

:::solution

$$
R_{\rm rad}=Z_0(2\pi/50)^2/(6\pi)=0.316\,\Omega
$$

である。
指向性はピークの $\sin^2\theta=1$ を全方向平均 $2/3$ で割り、$D=1.5$ となる。
:::
::::

::::exercise{level=advanced}
非偏光の Thomson 散乱で前半球と後半球の散乱断面積を求めよ。

:::solution
分布は $u=\cos\theta$ について偶である。
前半球積分は $2\pi r_e^2\int_0^1(1+u^2)\dd u/2=4\pi r_e^2/3$、後半球も同じである。
従ってそれぞれ全断面積の半分を占める。
:::
::::

::::exercise{level=advanced}
二個の等しい電荷 $q$ が原点を挟んで $\vb{z}_\pm(t)=\pm b\cos(\omega t)\hat{\vb{z}}$ と運動する。
速度は十分小さい。
電気双極子、磁気双極子、電気四重極のどれが最初に放射するかを調べ、その平均電力を求めよ。

:::solution
二電荷の位置和は零なので $\vb{p}=0$ である。
各速度は原点からの位置と平行で $\vb{m}=0$ でもある。
全電荷 $2q$ は時間変化しない。
残る四重極は $Q_{zz}=4qb^2\cos^2\omega t=2qb^2(1+\cos2\omega t)$、$Q_{xx}=Q_{yy}=-Q_{zz}/2$ である。

従って放射周波数は $2\omega$ である。
変動振幅 $Q_0=2qb^2$ を軸対称四重極の式に入れて

$$
\langle P\rangle=\frac{(2\omega)^6(2qb^2)^2}{960\pi\varepsilon_0c^5}=\frac{4q^2b^4\omega^6}{15\pi\varepsilon_0c^5}.
$$

この式は $2\omega b/c\ll1$ での先頭項である。
それぞれの電荷に Larmor の公式を適用して単純に足せば、打ち消される双極子場の電力まで数えてしまう。
:::
::::

::::exercise{level=advanced}
真空中の点状散乱体が純粋に電気双極子として応答する。
線形・受動の条件だけから、この一つの散乱チャネルがもちうる最大の散乱断面積を波長 $\lambda$ で表せ。

:::solution
無次元量 $z=k^3\alpha/(6\pi\varepsilon_0)$ を導入する。
本文の全散乱・消衰の式から $\sigma_s=(6\pi/k^2)|z|^2$、

$$
\sigma_{\rm ext}=(6\pi/k^2)\Im z
$$

である。

吸収が負にならない条件は

$$
|z|^2\le\Im z,\qquad (\Re z)^2+(\Im z-1/2)^2\le1/4.
$$

従って許される $z$ は複素平面で中心 $i/2$、半径 $1/2$ の円板内にある。
原点から最遠の点は $z=i$ であり、そこで $\sigma_s=6\pi/k^2=3\lambda^2/(2\pi)$、$\sigma_{\rm abs}=0$ となる。

この上限は一つの等方電気双極子チャネルについてのもので、より大きな散乱体には別の多重極も加わる。
:::
::::

::::exercise{level=advanced}
非相対論的な荷電粒子が短い時間だけ加速し、始めと終わりは等速運動であった。
遠方の放射電場の時間積分を、速度の変化で表せ。
観測方向と距離は加速区間内で一定と近似する。

:::solution
非相対論的な加速度場は

$$
\vb{E}_{\rm rad}=-q\vb{a}_\perp(t-r/c)/(4\pi\varepsilon_0c^2r)
$$

である。
遅延時刻で積分しても $\dd t_*=\dd t$ がこの近似の精度で成り立つので、

$$
\int_{-\infty}^\infty\vb{E}_{\rm rad}\dd t=-\frac{q(\vb{v}_f-\vb{v}_i)_\perp}{4\pi\varepsilon_0c^2r}.
$$

同じ速度に戻る運動では電場の時間積分が零になるが、電場の二乗で決まる放射エネルギーは正になりうる。
例えば一度加速し、時間を置いて逆向きに減速する二つのパルスは、電場の面積が相殺してもエネルギーは足し合わされる。
:::
::::

::::exercise{level=advanced}
同位相で規格化した電気・磁気双極子が $\vb d=\hat{\vb x},\vb h=\hat{\vb y}$ を持ち、四重極が零である。前半球で集めたエネルギーは全体の何割か。また縦の反動を全放射電力で表せ。

:::solution
$u=\cos\theta$ とすると角度分布は $(1+u)^2$ に比例する。$\int_0^1(1+u)^2\dd u=7/3$、$\int_{-1}^1(1+u)^2\dd u=8/3$ だから前半球の割合は $7/8$ である。縦の運動量にはさらに $u/c$ を掛ける。$\int_{-1}^1u(1+u)^2\dd u=4/3$ なので光の縦運動量は単位時間あたり $P/(2c)$、他の外力がなければ源の反動は $-P\hat{\vb z}/(2c)$ である。
:::
::::

## 確認テスト

::::quiz{answer="16" tol="1e-8"}
Rayleigh 領域で分極率を一定として波長を半分にすると断面積は何倍になるか。

:::explanation
誘起双極子の平均放射は $k^4|\alpha|^2$ に比例する。
波長を半分にすれば波数は倍になるので $16$ 倍になる。
:::
::::

::::quiz
放射場の距離依存性を選べ。

- [ ] $R^{-3}$
- [x] $R^{-1}$
- [ ] $R^{-2}$

:::explanation
強度は場の二乗で、球面積は $R^2$ に比例する。
遠方の有限な電力を担う項は場で $R^{-1}$ である。
:::
::::

::::quiz{answer="1.5" tol="1e-8"}
理想的な短い電気双極子アンテナの指向性を答えよ。

:::explanation
指向性は最大強度と方向平均の比である。
$\sin^2\theta$ の平均は $2/3$、最大は $1$ なので $3/2$ となる。
:::
::::

::::quiz
同じ $a$ で Liénard の公式が縦加速に与える倍率を、非相対論的な値に対して選べ。

- [x] $\gamma^6$
- [ ] $\gamma^4$
- [ ] $\gamma^2$

:::explanation
縦加速では $\vb{\beta}\times\vb{a}=0$ のため、Liénard の公式の括弧は $a^2$ となり $\gamma^6$ がそのまま残る。
:::
::::

::::quiz
電気四重極テンソルを $Q_{ij}=\int \rho(3x_ix_j-r^2\delta_{ij})\dd^3r$ と定義する。
全四重極放射電力に現れるものを選べ。

- [ ] $Q_{ij}$ 自身の二乗和
- [ ] 二階時間微分の二乗和
- [x] 三階時間微分の二乗和

:::explanation
電流の空間一次モーメントの対称部分に、電荷の二次モーメントの一階微分が現れる。
源内の遅延の展開と放射電場の時間微分でさらに二階が加わり、放射電場が三階微分に比例する。
その二乗を角度積分したものが電力である。
:::
::::

::::quiz{answer="0.5" tol="1e-8"}
無偏光の Thomson 散乱で、全散乱電力のうち後半球に向かう割合を答えよ。

:::explanation
角度分布は $1+\cos^2\theta$ に比例して前後対称である。
同じ面積の前後半球の積分が等しいため、割合は $1/2$ となる。
:::
::::

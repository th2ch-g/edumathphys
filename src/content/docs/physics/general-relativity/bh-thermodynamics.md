---
title: "ブラックホール熱力学"
description: "定常ブラックホールの地平線力学を導き、表面重力、面積、熱力学と半古典蒸発の関係を説明する。"
notation:
  - symbol: "\\kappa"
    meaning: "無限遠で規格化した表面重力"
  - symbol: "A_H"
    meaning: "地平線の断面積"
  - symbol: "\\Omega_H,\\Phi_H"
    meaning: "地平線の角速度と電位"
  - symbol: "S_{\\rm BH},T_H"
    meaning: "ブラックホールのエントロピーと Hawking 温度"
---

## この章の目標

- Killing 地平線の表面重力を計算する。
- Raychaudhuri 方程式から面積が減らない条件を説明する。
- Kerr–Newman 族のパラメータを変分して力学の第一法則を得る。
- 半古典的な温度、エントロピー、蒸発の帰結を調べる。

## Killing 地平線と表面重力

:term[Killing 地平線]{en="Killing horizon" yomi="きりんぐちへいせん"} は Killing 場 $\chi$ が非零な部分で法線となる零超曲面である。Killing 方程式から $\chi^b\nabla_b\chi_a=-\tfrac12\nabla_a(\chi^2)$ であり、面上の $\chi^2=0$ からその勾配は面に垂直である。従って

$$
\chi^b\nabla_b\chi^a=\kappa\chi^a,\qquad \nabla_a(\chi^2)=-2\kappa\chi_a
\label{eq-gr-surface-gravity}
$$

を満たす係数 $\kappa$ がある。これが :term[表面重力]{en="surface gravity" yomi="ひょうめんじゅうりょく"} である。$\chi$ を定数倍すれば $\kappa$ も同じ倍になる。Schwarzschild では無限遠での単位時間に合わせて $\chi=\partial_v$ と取る。:eqref[eq-gr-ingoing-ef] の地平線上で $\chi_a=(\dd r)_a$、$\chi^2=-f$ なので $\kappa=f'(2m)/2=1/(4m)$ である。静止観測者の固有加速度は $f'/(2\sqrt f)$ だが、これに赤方偏移係数 $\sqrt f$ を掛けて地平線に近づけると $\kappa$ になる。

Kerr–Newman では $m=GM$、$a=J/M$ とし、$q$ を幾何単位の電荷とする。$\Delta=r^2-2mr+a^2+q^2$ の二根を $r_+\ge r_-$、$\Sigma=r^2+a^2\cos^2\theta$、$\Sigma_+=\Sigma(r_+,\theta)$、$\mathcal D=r_+^2+a^2$ とする。地平線を生成する Killing 場を $\chi=\partial_t+\Omega_H\partial_\phi$、$\Omega_H=a/\mathcal D$ と取る。ここで無限遠で規格化するのは時間並進成分 $\partial_t$ であり、回転を含む $\chi$ 自身は遠方で一般に空間的になる。

動径の交差項が現れる正則な入射座標では $\chi_a=(\Sigma_+/\mathcal D)(\dd r)_a$ である。一方、$\chi^2$ の動径微分は地平線で $-\Delta'(r_+)\Sigma_+/\mathcal D^2$ となる。下で両式を計算する。:eqref[eq-gr-surface-gravity] に代入すると $\kappa=(r_+-r_-)/(2\mathcal D)$ であり、二根が一致する極限では零である。

### 入射座標で回転地平線の法線を求める

Kerr–Newman 計量の時間と方位角を含む部分を
$$
-\frac{\Delta}{\Sigma}(\dd t-a\sin^2\theta\,\dd\phi)^2
+\frac{\sin^2\theta}{\Sigma}[(r^2+a^2)\dd\phi-a\dd t]^2
$$
と書く。$\dd v=\dd t+(r^2+a^2)\dd r/\Delta$、$\dd\tilde\phi=\dd\phi+a\dd r/\Delta$ と置くと、第二の括弧の動径項は相殺する。第一の括弧は $\dd v-a\sin^2\theta\,\dd\tilde\phi-\Sigma\dd r/\Delta$ となる。二乗から出る $-\Sigma\dd r^2/\Delta$ が元の動径計量を消し、$2(\dd v-a\sin^2\theta\,\dd\tilde\phi)\dd r$ が残る。$\chi=\partial_v+\Omega_H\partial_{\tilde\phi}$ を地平線で各括弧に入れると、第二の括弧は零であり、交差項だけが $\chi^\flat$ に寄与する。係数は $1-a\Omega_H\sin^2\theta=\Sigma_+/\mathcal D$ である。

一般の $r$ での $\chi^2$ は第一の括弧から $-\Delta(1-a\Omega_H\sin^2\theta)^2/\Sigma$、第二から $\sin^2\theta[(r^2+a^2)\Omega_H-a]^2/\Sigma$ を得る。後者の括弧は $r_+$ で零なので第一微分にも寄与しない。従って $\partial_r\chi^2|_{r_+}=-\Delta'(r_+)\Sigma_+/\mathcal D^2$ で、角微分は地平線全体で $\chi^2=0$ のため零である。これで表面重力の比から緯度依存が消える理由が分かる。例えば $a=0$ なら $\mathcal D=r_+^2$ で、入射計量の $f'(r_+)/2$ と一致する。

### Killing パラメータとアフィン径数

地平線上で $\chi=\dd/\dd v$ とし、$\kappa$ が生成子に沿って一定である場合を考える。$\chi$ 流は計量と $\chi$ を保つので、:eqref[eq-gr-surface-gravity] の両辺を流れに沿って運ぶと係数 $\kappa$ も保たれる。アフィン接線を $k=(\dd v/\dd\lambda)\chi$ とする。その加速度が零であるには、$\dd^2v/\dd\lambda^2+\kappa(\dd v/\dd\lambda)^2=0$ が必要である。従って $\dd\lambda/\dd v$ は $e^{\kappa v}$ に比例し、$\kappa\ne0$ では $\lambda-\lambda_0=B e^{\kappa v}$ となる。$\kappa=0$ なら両パラメータはアフィン変換で結ばれる。例えば $\kappa>0$ の二分岐地平線では $v\to-\infty$ の極限が有限のアフィン径数 $\lambda_0$ に対応する。パラメータの無限大だけで、曲線がそこで幾何学的に不完備とは判定できない。

:::example[Rindler 地平線における符号と規格化]
Minkowski 平面の零座標 $U=t-x,V=t+x$ に対し、boost 場を $\chi=a(V\partial_V-U\partial_U)$、$a>0$ とする。未来枝 $U=0,V>0$ では $\nabla_\chi\chi=a^2V\partial_V=a\chi$ だから $\kappa=a$ である。過去枝 $V=0,U<0$ では $\nabla_\chi\chi=a^2U\partial_U=-a\chi$ である。$\chi$ が零になる交差点ではこの式だけから $\kappa$ を割って求められないが、各枝の極限を取れる。右楔の $\rho=\sqrt{-UV}$ では $\chi^2=-a^2\rho^2$ で、$\rho=1/a$ 上で $\chi$ のノルムが $-1$ になる。この例には無限遠の時間並進による規格化がなく、どの軌道を基準時計にするかが係数を決める。
:::

:::example[宇宙論的 Killing 地平線との比較]
de Sitter の静的領域 $f=1-H^2r^2>0$ を用い、中心で $\partial_t$ を単位に取る。無限遠での規格化を持つブラックホール地平線と比べるため、局所加速度を求める。 座標の順序を $(t,r,\theta,\phi)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}H^{2} r^{2} - 1 & 0 & 0 & 0\\0 & \frac{1}{- H^{2} r^{2} + 1} & 0 & 0\\0 & 0 & r^{2} & 0\\0 & 0 & 0 & r^{2} \sin^{2}{\left(\theta \right)}\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}\frac{1}{\left(H r - 1\right) \left(H r + 1\right)} & 0 & 0 & 0\\0 & - \left(H r - 1\right) \left(H r + 1\right) & 0 & 0\\0 & 0 & \frac{1}{r^{2}} & 0\\0 & 0 & 0 & \frac{1}{r^{2} \sin^{2}{\left(\theta \right)}}\end{matrix}\right]\\
\det g&=- r^{4} \sin^{2}{\left(\theta \right)}
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{t}_{t r}$ | $\frac{1}{2 \left(H r - 1\right) \left(H r + 1\right)}\left(\partial_{r}\left(H^{2} r^{2} - 1\right)\right)$ | $\frac{H^{2} r}{\left(H r - 1\right) \left(H r + 1\right)}$ |
| $\Gamma^{r}_{t t}$ | $- \frac{\left(H r - 1\right) \left(H r + 1\right)}{2}\left(-\partial_{r}\left(H^{2} r^{2} - 1\right)\right)$ | $H^{2} r \left(H r - 1\right) \left(H r + 1\right)$ |
| $\Gamma^{r}_{r r}$ | $- \frac{\left(H r - 1\right) \left(H r + 1\right)}{2}\left(\partial_{r}\left(\frac{1}{- H^{2} r^{2} + 1}\right)+\partial_{r}\left(\frac{1}{- H^{2} r^{2} + 1}\right)-\partial_{r}\left(\frac{1}{- H^{2} r^{2} + 1}\right)\right)$ | $- \frac{H^{2} r}{\left(H r - 1\right) \left(H r + 1\right)}$ |
| $\Gamma^{r}_{\theta \theta}$ | $- \frac{\left(H r - 1\right) \left(H r + 1\right)}{2}\left(-\partial_{r}\left(r^{2}\right)\right)$ | $r \left(H r - 1\right) \left(H r + 1\right)$ |
| $\Gamma^{r}_{\phi \phi}$ | $- \frac{\left(H r - 1\right) \left(H r + 1\right)}{2}\left(-\partial_{r}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $r \left(H r - 1\right) \left(H r + 1\right) \sin^{2}{\left(\theta \right)}$ |
| $\Gamma^{\theta}_{r \theta}$ | $\frac{1}{2 r^{2}}\left(\partial_{r}\left(r^{2}\right)\right)$ | $\frac{1}{r}$ |
| $\Gamma^{\theta}_{\phi \phi}$ | $\frac{1}{2 r^{2}}\left(-\partial_{\theta}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $- \frac{\sin{\left(2 \theta \right)}}{2}$ |
| $\Gamma^{\phi}_{r \phi}$ | $\frac{1}{2 r^{2} \sin^{2}{\left(\theta \right)}}\left(\partial_{r}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{r}$ |
| $\Gamma^{\phi}_{\theta \phi}$ | $\frac{1}{2 r^{2} \sin^{2}{\left(\theta \right)}}\left(\partial_{\theta}\left(r^{2} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{\tan{\left(\theta \right)}}$ |

アフィン径数の微分を点で表し、$w_a$ を座標 $a$ の速度成分とする。ここでは添字を下げる操作を行っていない。測地線の右辺 $-\Gamma^a{}_{bc}w_bw_c$ では $b\ne c$ の項を二度数える。

$$
\begin{aligned}
\ddot{t}&=- \frac{2 H^{2} r w_{r} w_{t}}{\left(H r - 1\right) \left(H r + 1\right)}\\
\ddot{r}&=- \frac{r \left(H^{6} r^{4} w_{t}^{2} + H^{4} r^{4} w_{\phi}^{2} \sin^{2}{\left(\theta \right)} + H^{4} r^{4} w_{\theta}^{2} - 2 H^{4} r^{2} w_{t}^{2} - 2 H^{2} r^{2} w_{\phi}^{2} \sin^{2}{\left(\theta \right)} - 2 H^{2} r^{2} w_{\theta}^{2} - H^{2} w_{r}^{2} + H^{2} w_{t}^{2} + w_{\phi}^{2} \sin^{2}{\left(\theta \right)} + w_{\theta}^{2}\right)}{\left(H r - 1\right) \left(H r + 1\right)}\\
\ddot{\theta}&=\frac{r w_{\phi}^{2} \sin{\left(2 \theta \right)} - 4 w_{\theta} w_{r}}{2 r}\\
\ddot{\phi}&=- \frac{2 w_{\phi} \left(r w_{\theta} + w_{r} \tan{\left(\theta \right)}\right)}{r \tan{\left(\theta \right)}}
\end{aligned}
$$

また計量ノルムは

$$
g(\dot x,\dot x)=r^{2} w_{\phi}^{2} \sin^{2}{\left(\theta \right)} + r^{2} w_{\theta}^{2} + \frac{w_{r}^{2}}{- H^{2} r^{2} + 1} + w_{t}^{2} \left(H^{2} r^{2} - 1\right)
$$
である。これを径数で微分して上の加速度を代入すると、係数の微分と速度の微分の和が零になる。従って初期のノルムが保存される。

静止観測者を保つ加速度の動径成分は $f'/2=-H^2r$ である。従って固有加速度の大きさに $\sqrt f$ を掛けると $H^2r$ であり、地平線 $r=H^{-1}$ に近づくにつれて $H$ へ行く。零座標で表面重力の符号を付ける際には、選んだ未来地平線への拡張と Killing ベクトルの向きを固定する。温度で使うのはこの規格化における大きさである。
:::

### 非退化の極限と極限解の固有距離

球対称な荷電解を $f=(r-r_+)(r-r_-)/r^2$ とし、$\delta=r_+-r_->0$ とする。静的時刻面で $r=r_++x$、$x>0$ なら
$$
\dd\ell=\frac{r_++x}{\sqrt{x(x+\delta)}}\,\dd x .
$$
$x\ll r_+$ で分子を $r_+$ とすれば、地平線から $x$ までの距離は $2r_+\operatorname{arsinh}\sqrt{x/\delta}$ である。$x\ll\delta$ なら $2r_+\sqrt{x/\delta}$、一方、固定した小さな $x$ に対して $\delta\to0$ を取ると $r_+\log(4x/\delta)$ が先頭となり発散する。$\delta=0$ の厳密な極限解で直接計算しても $\dd\ell=(r_++x)\dd x/x$ なので、外の任意の点から地平線への距離は無限である。従って非退化の近傍を Euclid 平面の原点として扱う手続きでは $x\ll\delta$ が重要であり、その近傍は極限に向かうにつれて狭くなる。

## Hawking の面積定理

滑らかな事象地平線の生成子が未来完備で地平線にとどまる場合をまず考える。この仮定は、初期膨張率から得る焦点までのアフィン径数を確保するために使う。光的超曲面の生成子には渦度がなく、零収束条件 $R_{ab}k^ak^b\ge0$ の下で :eqref[eq-gr-raychaudhuri-null] は $\theta'\le-\theta^2/2$ を与える。もし $\theta<0$ なら有限のアフィン径数で焦点を作り、その後の生成子は非時間的関連な因果境界にとどまれない。従って $\theta\ge0$ である。$\dd(\log\dd A)/\dd\lambda=\theta$ なので既存の生成子が運ぶ面積要素は減らない。

この形の :term[Hawking の面積定理]{en="Hawking area theorem" yomi="ほーきんぐのめんせきていり"} は、同じ生成子を一度ずつ横切る滑らかな二つの断面のうち、未来側の面積が以上になると述べる。より一般の事象地平線では新たな生成子が過去端点から入り、断面に寄与する。その場合に面積を比較するには、地平線の非滑らかな集合と断面間の写像も制御する。未来完備性を直接仮定しない古典的な証明では、外部の強い漸近予測可能性の下で因果境界が保たれることを使う。この大域的な扱いは :cite[hawking-ellis] を参照する。

例えば二つの非回転穴の初期質量がそれぞれ $M$、合体後も非回転なら $16\pi G^2M_f^2\ge32\pi G^2M^2$ なので $M_f\ge\sqrt2 M$ である。初期の総エネルギーを $2M$ と近似すれば、放射可能な割合のこの上界は $1-1/\sqrt2\simeq0.293$ となる。これは面積だけから得る上限であり、実際の合体の出力を求めるには時間発展が必要である。

### 収束による面積増加を球で検算する

入射する球対称の正エネルギー流で $\dd s^2=-[1-2m(v)/r]\dd v^2+2\dd v\dd r+r^2\dd\Omega^2$ とする。外向き生成子の傾きは $\dd r/\dd v=(1-2m(v)/r)/2$ である。最終的に $m(v)=m_f$ となるなら、未来事象地平線を追うには終端条件 $r=2m_f$ からこの方程式を過去へ積分する。一方、$\dd r/\dd v=0$ の球は各時刻の $r=2m(v)$ である。これは膨張が零になる局所条件であり、同じ面が常に地平線の生成子であるためには $m'=0$ が必要である。$\dot m\ge0$ で最終定常な模型では、差 $d=r_H-2m$ は
$$
d'-\frac{d}{2r_H}=-2m'
$$
を満たす。未来で $d=0$ の条件から後ろ向きに積分すれば $d\ge0$ で、$r'_H=d/(2r_H)\ge0$ である。従って $4\pi r_H^2$ が非減少であることを、球対称の方程式から直接得る。

### 断面の面積を生成子に沿って比較する

初めの断面で角座標 $x^A$ を取り、生成子に沿ってその値を一定に保つ。$\lambda$ が生成子ごとのアフィン径数であるとする。二つの断面が $\lambda_1(x),\lambda_2(x)$ で与えられるとき、各断面の接線は $\partial_A+(\partial_A\lambda_i)k$ である。地平線上では $k^2=0$、$g(k,\partial_A)=0$ なので、誘導計量には $\partial_A\lambda_i$ が現れない。従って同じ角座標の面積密度の比は
$$
\frac{\sqrt{\det q(\lambda_2(x),x)}}{\sqrt{\det q(\lambda_1(x),x)}}
=\exp\left(\int_{\lambda_1(x)}^{\lambda_2(x)}\theta(\lambda,x)\dd\lambda\right)\ge1
$$
である。ここでは $\lambda_2(x)\ge\lambda_1(x)$ とし、両断面の間で生成子のラベルが有効である。従って面積の比較は同時刻の切断だけに限られない。等号の場合、滑らかな部分では連続な非負の $\theta$ の積分が零なので $\theta$ は零である。Raychaudhuri 方程式に戻ると、その区間でせん断と $R_{kk}$ も零である。物質に零エネルギー条件を課した Einstein 方程式なら $T_{kk}=0$ も従う。例えば定常な Killing 地平線では $\mathcal L_\chi q=0$ なので、地平線上で面積密度が変わらない。

### 殻が来る前の事象地平線

薄い入射殻を $v=0$ に置き、$v<0$ は平坦、$v>0$ は質量パラメータ $m$ の Schwarzschild とする。殻の外側の未来地平線は $r=2m$ である。先の外向き光路の式を平坦領域に続けると $r_H=2m+v/2$ となり、$v=-4m$ で中心から出る。面積は $-4m<v<0$ で $4\pi(2m+v/2)^2$ であり、殻の到着以前に増える。一方、この平坦部分ではどの正の半径の球も外向き膨張率が正である。将来まで逃げられるかで事象地平線を選ぶことと、各球の膨張率から局所的に限界面を探すことは異なる操作である。無限に薄い殻は理想化であり、滑らかな入射分布では $r_H$ を同じ終端条件から求められる。

### 収束不等式の時間尺度

:::example[負の膨張率で進む光束を地平線に置けるか]
光束に $\theta(0)=-1/L$ を与え、Ricci 項とせん断が零とする。$\theta'=-\theta^2/2$ の解は $\theta=-1/[L-\lambda/2]$ である。面積要素は $\dd\log A/\dd\lambda=\theta$ を積分して $A(\lambda)=A(0)(1-\lambda/(2L))^2$ となり、$\lambda=2L$ で焦点を作る。生成子を無限の未来まで持つ滑らかな事象地平線の一部にこの区間を置けば、焦点の後に非時間的関連な境界であり続ける条件と矛盾する。初期の膨張率の符号を代入するだけで、この局所計算から大域条件に必要な矛盾を作れる。
:::

## ブラックホール力学の 4 法則

### 零法則と第一法則

零法則は定常平衡の地平線で $\kappa$ が一定となる条件を述べる。滑らかな二分岐 Killing 地平線では、交差断面で $\chi=0$、$\nabla_a\chi_b$ は断面接方向で零で反対称だから法線二平面の boost 生成子になる。そこで $\kappa^2=-(\nabla_a\chi_b)(\nabla^a\chi^b)/2$ である。Killing 恒等式により $\nabla_c\nabla_a\chi_b$ は曲率と $\chi$ の積だから断面で零であり、$\kappa^2$ の接方向微分は零である。Killing 流でも一定なので連結な地平線に広がる。一般の定常な地平線の零法則には物質条件を用いる形もあり、:cite[wald-gr] を参照する。

第一法則を Kerr–Newman 族で確認する。$m=GM$、$J=Ma$、$Q=4\pi Q_0$、$q^2=GQ^2/(4\pi)$、$\Phi_H=Q_0r_+/\mathcal D$ とする。根の式 $r_+^2-2mr_++a^2+q^2=0$ を微分すると

$$
(r_+-m)\delta r_+=r_+\delta m-a\delta a-q\delta q .
$$

$A_H=4\pi\mathcal D$ なので $\kappa\delta A_H/(8\pi)=(r_+-m)(r_+\delta r_++a\delta a)/\mathcal D$ である。ここへ $\Omega_H\delta(ma)=a(a\delta m+m\delta a)/\mathcal D$ を加えると $\delta a$ が消え、$\delta m-q r_+\delta q/\mathcal D$ となる。全体を $G$ で割り、電荷の規格化を入れると

$$
\delta M=\frac{\kappa}{8\pi G}\delta A_H+\Omega_H\delta J+\Phi_H\delta Q .
\label{eq-gr-bh-first-law}
$$

長さを一様に倍する操作で $M,Q$ は長さの一次、$A_H,J$ は二次なので、同じ式にこの変分を使えば $M=\kappa A_H/(4\pi G)+2\Omega_H J+\Phi_H Q$ を得る。これは Smarr 式である。

### 地平線の Komar 積分から Smarr 式を導く

定常軸対称の真空で、地平線から無限遠までの領域に他の境界がないとする。地平線の生成子 $\chi=\xi+\Omega_H\varphi$ に対し、Komar の質量の規約で無限遠の積分は $M-2\Omega_HJ$ である。ここで遠方の $M,J$ は時間並進と回転それぞれの Komar 電荷であり、この定常真空の設定で ADM の値と一致する。真空では $*\dd\chi^\flat$ が閉形式なので、同じ積分を地平線上の断面へ移せる。

そこで別の光的法線 $\ell$ を $\chi\cdot\ell=-1$ に選ぶ。Killing 方程式と $\nabla_\chi\chi=\kappa\chi$ から
$$
(\dd\chi^\flat)_{ab}\chi^a\ell^b=2\ell_b\nabla_\chi\chi^b=-2\kappa .
$$
時間と外向きの向きを無限遠の球に合わせると、断面への Hodge 双対の引き戻しは $-2\kappa\,\dd A$ である。従って $M-2\Omega_H J=\kappa A_H/(4\pi G)$ となる。この証明は面積の係数を法線平面から得る。電磁場がある場合には $\dd *\dd\chi^\flat$ の Ricci 源が零でなく、同じ真空の Stokes の議論では $\Phi_H Q$ を落としてしまう。

実際、非回転の荷電穴なら半径 $r$ の Komar 質量は $(m-q^2/r)/G$ である。$r=r_+$ で $\kappa A_H/(4\pi G)=(r_+-r_-)/(2G)=(m-q^2/r_+)/G$ となる。差 $q^2/(Gr_+)=\Phi_H Q$ は外部電場を含めた体積積分そのものである。面積だけを読む地平線の積分と、無限遠で読む質量を同一視できるかどうかは外側の物質によって変わる。

### Kerr–Newman の地平線電位をそろえる

電磁場の作用を $-\tfrac14\int\sqrt{-g}F_{ab}F^{ab}$ とする有理化規約で、静的遠方の電荷を $Q=4\pi Q_0$ と読む。定常ゲージを無限遠で $A(\xi)\to0$ と選べば、$A=-Q_0 r(\dd t-a\sin^2\theta\dd\phi)/\Sigma$ である。地平線で
$$
-A(\chi)=\frac{Q_0 r_+}{\Sigma_+}(1-a\Omega_H\sin^2\theta)
=\frac{Q_0r_+}{\mathcal D}
$$
となり緯度に依存しない。第一法則では無限遠との差を $\Phi_H$ に使う。例えば時間に比例するゲージ関数で $A_t$ を定数だけ動かしても、無限遠と地平線で $A(\chi)$ は同じだけ動くので差は保たれる。なお $q^2=GQ^2/(4\pi)$ では $q\,\delta q=GQ\,\delta Q/(4\pi)$ が成り立つ。従って幾何学的根を微分して生じた $qr_+\delta q/(G\mathcal D)$ は $\Phi_H\delta Q$ である。

### 微小な吸収過程での面積変化

平衡な非退化地平線に小さな物質流を入れ、十分遅い時刻に再び平衡になるとする。Killing 時刻 $v$ を生成子 $\chi$ のパラメータに取ると $\nabla_\chi\chi=\kappa\chi$ なので、Raychaudhuri 方程式には非アフィン項 $+\kappa\theta$ が加わる。これはアフィン接線 $k$ から $\chi=bk$ へ替えたとき $\theta_\chi=b\theta_k$、$\chi(b)/b=\kappa$ を用いて得る。平衡の $\theta$ とせん断が零なので、一次で

$$
(\partial_v-\kappa)\delta\theta=-8\pi G T_{ab}\chi^a\chi^b
$$

となる。未来の平衡条件 $\delta\theta\to0$ を課すと、$\delta\theta(v)=8\pi G\int_v^\infty e^{-\kappa(v'-v)}T_{\chi\chi}(v')\dd v'$ である。初期・終端の両方で膨張が消える場合、時刻を積分して $\delta A_H=(8\pi G/\kappa)\int T_{\chi\chi}\dd v\,\dd A$ となる。この右辺の流束は生成子 $\chi$ に対応する流入エネルギーである。中性物質で定常な端状態なら $\delta M-\Omega_H\delta J$ に対応し、:eqref[eq-gr-bh-first-law] の面積係数が時空の収束からも得られる。摂動の二次ではせん断による重力波の寄与も入る。

:::example[短い吸収パルスに対する面積変化]
背景の断面で $\int T_{\chi\chi}\dd A$ が $F_0>0$ に等しい区間を $0<v<L$、それ以外を零とする。一次の断面積の微分は、前節の解を角積分して
$$
\frac{\dd\,\delta A_H}{\dd v}=\frac{8\pi G F_0}{\kappa}
\begin{cases}
e^{\kappa v}(1-e^{-\kappa L}),&v<0,\\
1-e^{-\kappa(L-v)},&0\le v\le L,\\
0,&v>L.
\end{cases}
$$
前の区間に非零の変化が現れるのは、事象地平線を将来の平衡条件で選んだためである。二つの非零区間を積分すると、過去側の $(1-e^{-\kappa L})/\kappa$ とパルス中のその負の補正が打ち消し合い、$\delta A_H=8\pi G F_0L/\kappa$ となる。量 $F_0L$ は流入エネルギーの摂動である。有限の振幅や極限地平線については、二次の膨張率やせん断を捨て、この一次式を $1/\kappa$ のまま利用できない。
:::

### 第二法則と第三法則

第二法則の古典版は前節の面積非減少である。第三法則の到達不能性は、非退化な状態から有限の物理過程で極限状態を作れるかを問う。古典的な到達不能性の定理には、見かけの地平線を追跡できること、近傍の応力が有界で弱いエネルギー条件を満たすことなどが要る。任意の時間依存計量で定常時の Killing 表面重力 $\kappa$ を各時刻に定義する議論ではない。内部の捕捉面と時間発展の制御を伴う定理については :cite[poisson-toolkit] を参照する。定常解の族では例えば $Q=0$ の極限 Kerr に $a=m$、$A_H=8\pi m^2>0$ の点がある。従って絶対零度への到達不能性と、零温度の極限でエントロピーが状態に依存しないという Nernst 型の主張は異なる。

### 非可逆な変化と回転エネルギー

Kerr の面積を $A_H=8\pi m r_+$ とする。面積から $m_{\rm irr}=\sqrt{A_H/(16\pi)}$ を定義すると $m_{\rm irr}^2=mr_+/2$ である。根の関係 $r_+^2-2mr_++a^2=0$ と $GJ=ma$ を使い、
$$
m^2=m_{\rm irr}^2+\frac{(GJ)^2}{4m_{\rm irr}^2}
$$
を得る。実際、右辺第二項は $ma^2/(2r_+)=m(2m-r_+)/2$ である。面積が減らない過程で $J$ を減らしても、最終的な $m$ は初期の $m_{\rm irr}$ 以上である。極限回転の初期状態なら $r_+=m$ であり、静止質量エネルギーの抽出割合の上限は $1-1/\sqrt2$ となる。これは古典的面積法則を使ったエネルギーの上限であり、ある具体的な過程で等号を実現するという結論は含まない。

### 全微分から回転と電荷の係数を導く

Kerr–Newman の地平線で $\mathcal D=r_+^2+a^2$、$a=J/M$ とする。$\mathcal D-2mr_++q^2=0$ より $2mr_+=\mathcal D+q^2$ である。$r_+^2=\mathcal D-a^2$ を使い二乗して整理すると
$$
m^2=\frac{(\mathcal D+q^2)^2}{4\mathcal D}
+\frac{(GJ)^2}{\mathcal D}.
$$
従って幾何学的面積は $\mathcal D$ で、$J$ と $q$ を独立変数として質量を微分できる。まず $\mathcal D,q$ 一定で $\partial M/\partial J=GJ/(m\mathcal D)=a/\mathcal D=\Omega_H$ である。また $\partial m/\partial q=q(\mathcal D+q^2)/(2m\mathcal D)=q r_+/\mathcal D$ で、$q^2=GQ^2/(4\pi)$ を入れると $\partial M/\partial Q=Qr_+/(4\pi\mathcal D)$ となる。定常族の面積項・回転項・電荷項は、同じ質量関数の三つの偏微分である。このため、異なる変数を固定した熱容量や感受率を比較する際は、何を変えたかも記録する。

:::example[電荷を保って回転をなくす場合の質量下限]
面積非減少と電荷保存だけを使うと、回転が零の最終状態の $m_f$ は $(\mathcal D_i+q^2)/(2\sqrt{\mathcal D_i})$ 以上である。ただし初期の外側地平線では $\mathcal D_i\ge q^2$ である。実際、関数 $(D+q^2)/(2\sqrt D)$ の導関数は $(D-q^2)/(4D^{3/2})$ で、$D\ge D_i$ で非負である。$q=0$ なら $\sqrt{\mathcal D_i}/2=m_{\rm irr}$ に戻る。面積の下限だけで境界を与えたので、放電がある過程にこの固定電荷の下限を直接適用しない。
:::

### 極限に近い穴へ電荷を加える試験粒子

Reissner–Nordström で $m>|q|>0$ とし、$q>0$ の場合を取る。静的な電位は $\Phi(r)=Q_0/r$、同符号の電荷 $e>0$ を持つ粒子の正準 Killing エネルギーは $E=-g(\xi,p)-e A(\xi)$ である。地平線上では $\xi$ が未来光的、粒子の運動量 $p$ が未来時間的なので $-g(\xi,p)>0$ である。従って地平線を横切る粒子は $E>e\Phi_H$ を満たす。この不等式は第一法則の $\delta M-\Phi_H\delta Q>0$ と同じ向きである。極限背景 $r_+=q$ の接線変分なら $G\,e\Phi_H=\delta q$ なので $\delta m>\delta q$ となる。

一方、非退化な背景では $G e\Phi_H=(q/r_+)\delta q<\delta q$ である。$m-q$ が小さいと、背景を保った試験粒子の条件と形式的な最終値 $m+GE<q+\delta q$ が同時に満たされる領域が生じ得る。しかし極限からの距離が二次の摂動と同程度なら、粒子の場と放射、地平線の変位を落とした一次計算では最終値を判定できない。この例の一次不等式は、極限への到達不能性の一般証明を与えるものではなく、第一法則を有限の操作へ外挿するとどの次数が問題になるかを示す。

## Bekenstein–Hawking エントロピー

粒子を落とすと外から見える物質の状態が失われ得るので、外部の物質エントロピーだけで第二法則を表せない。温度を次節の $T_H=\hbar\kappa/(2\pi k_B)$ とし、:eqref[eq-gr-bh-first-law] の面積項を $T_H\delta S$ と同定すると、非退化な族に沿って $\delta S=k_B\delta A_H/(4G\hbar)$ である。半古典的重力で加法定数を零に取った

$$
S_{\rm BH}=\frac{k_B A_H}{4G\hbar}
\label{eq-gr-bh-entropy}
$$

を :term[Bekenstein–Hawking エントロピー]{en="Bekenstein–Hawking entropy" yomi="べけんしゅたいんほーきんぐえんとろぴー"} という。SI では分子に $c^3$ を掛ける。Schwarzschild で $A_H=16\pi G^2M^2$ なので $S_{\rm BH}/k_B=4\pi GM^2/\hbar$ である。太陽質量でおよそ $10^{77}$ となる。外部場のエントロピー $S_{\rm out}$ を加えた $S_{\rm gen}=S_{\rm BH}+S_{\rm out}$ の非減少が一般化第二法則の主張である。量子論の $S_{\rm out}$ は短距離で発散するので、この和は重力側の定数も同じ規約でくりこむ必要がある。

### エントロピーと合体で失えるエネルギー

十分離れている二つの静止 Schwarzschild 穴が合体し、合体後も回転しないとする。初期エネルギーを $M_1+M_2$ と近似する。古典面積定理から $16\pi G^2 M_f^2\ge16\pi G^2(M_1^2+M_2^2)$ である。従って放射へ出る質量換算分は $M_1+M_2-\sqrt{M_1^2+M_2^2}$ 以下である。等質量では全質量のうち $1-1/\sqrt2$ が上限である。角運動量が残る場合、同じ $M_f$ での Kerr 面積は Schwarzschild 面積以下なので、この上限はさらに締まる。面積定理が適用される古典過程だけを仮定したので、波形や実際の放射効率の値を得るには別の運動方程式が必要である。

### Euclid 作用から Schwarzschild のエントロピーを求める

Euclid 計量を $f\dd\tau^2+f^{-1}\dd r^2+r^2\dd\Omega^2$、$f=1-2m/r$ とする。地平線を滑らかな原点として含めるには $\tau$ に周期 $\beta=8\pi m$ を与える。作用を $\hbar=1$ で書き、外側の境界を $r=R$ に置く。Ricci スカラーが零なので、基準背景を差し引く Gibbons–Hawking–York 項だけが残る。原点には境界を置かない。外向き法線は $n=\sqrt f\,\partial_r$ であり、四次元の体積密度は $r^2\sin\theta$ だから
$$
K=\frac1{R^2}\frac{\dd}{\dd R}(R^2\sqrt f)
=\frac{2\sqrt f}R+\frac{f'}{2\sqrt f}.
$$
同じ境界計量を持つ平坦な基準計量では円の固有周期が $\beta\sqrt f$、球面半径が $R$ で、外的曲率は $K_0=2/R$ である。従って
$$
I_E=-\frac{\beta}{2G}\left[2R(f-\sqrt f)+\frac{R^2 f'}2\right].
$$
$R\to\infty$ で $2R(f-\sqrt f)\to-2m$、$R^2f'/2=m$ であるため $I_E\to\beta m/(2G)=\beta^2/(16\pi G)$ となる。これはこの正則な Euclid 解の作用である。これを鞍点近似 $\log Z\simeq-I_E$ に入れると $E=-\partial_\beta\log Z=\beta/(8\pi G)=m/G$、さらに $S/k_B=\log Z+\beta E=\beta^2/(16\pi G)=A_H/(4G)$ となる。$\hbar$ を戻すと $A_H/(4G\hbar)$ である。負の熱容量があるため、無限遠で熱浴と安定に平衡すると主張するためにこの鞍点近似だけを使うことはできない。ここでは滑らかさの周期と境界項が面積則の係数にどう入るかを検算した。

### 地平線近傍で物体を静かに下ろす模型

非退化の静的地平線の近くでは、地平線からの固有距離 $\ell$ に対して lapse は $N\simeq\kappa\ell$ である。局所エネルギー $E_{\rm loc}$ の物体を半径の尺度 $R$ に対して $\ell\simeq R$ までゆっくり下ろす。物体が小さく、自己重力と支持装置の寄与を無視できる理想化では、放した後に穴へ入る無限遠エネルギーは $E_\infty\simeq \kappa R E_{\rm loc}$ である。従って第一法則と $T_H$ による穴のエントロピー増分は
$$
\Delta S_{\rm BH}\simeq\frac{E_\infty}{T_H}
\simeq\frac{2\pi k_B R E_{\rm loc}}{\hbar}
$$
となる。$c$ を戻すと右辺の分母は $\hbar c$ である。この尺度は、物体が外部から消えるときに一般化エントロピーの増加を考える動機になる。実際の支持や放射、量子場による浮力を考慮した一般の第二法則の証明ではなく、$R$ より近くで落とす設定では物体全体が地平線の外にある条件も必要になる。

:::example[同じエネルギーでサイズを変える]
同じ $E_{\rm loc}$ の物体について、サイズを $R$ から $2R$ に変え、両方とも $\kappa R\ll1$ とする。放す位置の最小の距離を各サイズに合わせると、無限遠のエネルギーと穴のエントロピー増分はともに二倍になる。一方、同じ物体と高さを固定して $\kappa$ だけを変えると $E_\infty$ と $T_H$ がともに $\kappa$ に比例し、比の先頭は変わらない。局所でのエネルギーと、穴の質量を変える無限遠でのエネルギーを分けることでこの相殺が見える。
:::

### 円錐の曲率から面積項を取り出す

非退化 Schwarzschild の Euclid 計量で、周期 $\beta$ を平滑条件 $\beta_H=8\pi m$ から独立に変える。以下では $\hbar=k_B=1$ とする。$\alpha=\beta/\beta_H$ とすれば、地平線に直交する小円の周長は $2\pi\alpha\ell$ である。先端を半径 $\ell$ の滑らかな帽子で置き換え、二次元の Gauss–Bonnet 公式を用いる。周上の測地曲率は $\ell^{-1}$、その積分が $2\pi\alpha$ なので、帽子に集中する Gauss 曲率の積分は $2\pi(1-\alpha)$ になる。スカラー曲率は Gauss 曲率の二倍であり、四次元で残りの地平線断面を積分すると $4\pi(1-\alpha)A_H$ が残る。角面積の $\ell$ による変化は高次なので帽子の極限で消える。

Euclid 作用の符号を用いると、この寄与は $-(1-\alpha)A_H/(4G)$ である。一方、外の境界項は周期に比例して $\beta M/2$ である。$\alpha A_H/(4G)=\beta M/2$ を代入すると、外から固定した $M$ と周期 $\beta$ に対する作用は
$$
I_E(\beta,M)=\beta M-A_H(M)/(4G)
$$
となる。従って $(\beta\partial_\beta-1)I_E$ を $M$ 固定で計算すると $A_H/(4G)$ が残る。また $\partial_M I_E=\beta-8\pi G M$ であり、$M$ に関する停留条件は円錐の欠損がない条件に一致する。ただし $\partial_M^2 I_E=-8\pi G<0$ であり、無限遠で温度を固定したまま $M$ を変えるこの一変数の模型では安定な最小値にならない。

:::example[周期を半分にした幾何の作用]
$m$ を固定して $\beta=\beta_H/2$ を取る。このとき外の境界の寄与は $\beta_H M/4=A_H/(8G)$ で、円錐の寄与は $-A_H/(8G)$ である。和は零になる。正則な解で作用が零になるという意味ではない。正則な周期へ戻すと円錐の寄与が消え、外側の寄与だけから $I_E=A_H/(4G)$ を得る。周期を動かす計算で先端の曲率を忘れると、面積項が出ない。
:::

### 一般化エントロピーを熱浴モデルで調べる

温度 $T$ の大きい熱浴へ小さい系を落とし、浴の温度変化を一次では無視できる模型を考える。系の Hamiltonian $H$ の Gibbs 状態を $\sigma=e^{-H/(k_BT)}/Z$、系の初期状態を $\rho$ とする。無次元の von Neumann エントロピーを $s(\rho)=-\tr(\rho\log\rho)$、エネルギーを $E(\rho)=\tr(\rho H)$ とすれば、
$$
D(\rho\|\sigma)=\tr[\rho(\log\rho-\log\sigma)]
=-s(\rho)+\frac{E(\rho)}{k_BT}+\log Z .
$$
$\rho$ の系を浴に吸収し、$\sigma$ の系を外に戻す理想化で、浴のエネルギー変化は $E(\rho)-E(\sigma)$、外部のエントロピー変化は $k_B[s(\sigma)-s(\rho)]$ である。従って両者のエントロピー変化の和は $k_B D(\rho\|\sigma)\ge0$ になる。ここで浴をブラックホールと考えると、浴の一次のエントロピー変化に第一法則を用いたことに相当する。

有限次元で正値な $\sigma$ について、この不等式を確かめる。$\rho,\sigma$ の固有値を $p_i,q_j$、固有基底の重なりの絶対値二乗を $w_{ij}$ とする。各行と列の和は $1$ である。従って $D=\sum_{ij}p_i w_{ij}\log(p_i/q_j)$ である。$p_i>0$ の項で $-\log x\ge1-x$ を使うと $D\ge\sum_{i:p_i>0,j}(p_i-q_j)w_{ij}\ge0$ となる。最後の不等式は、$p_i$ の和が $1$ で $q_j$ の寄与の和が $1$ 以下だからである。例えばエネルギーがともに零の二準位で $\rho$ が純粋、$\sigma=I/2$ なら $\Delta S_{\rm gen}=k_B\log2$ である。一方、初めから $\rho=\sigma$ なら各変化は零である。この模型は有限系と平衡浴での収支を示しており、時空中の一般の分割で $S_{\rm out}$ をくりこんで定義する手続きは別に要る。

### 面積エントロピーが作用の規格化を測る

Einstein–Hilbert ラグランジアンのスカラー部を $L=R/(16\pi G)$ とし、計量を固定して曲率について微分する。$R$ は曲率の二重縮約なので、曲率の反対称性を入れると
$$
P^{abcd}=\frac{\partial L}{\partial R_{abcd}}
=\frac{1}{32\pi G}(g^{ac}g^{bd}-g^{ad}g^{bc})
$$
である。地平線断面の法線二平面の体積二形式を $\epsilon_{ab}$ とし、$\epsilon_{ab}\epsilon^{ab}=-2$ に規格化する。$\epsilon$ の縮約では括弧の二項が同符号になり、$P^{abcd}\epsilon_{ab}\epsilon_{cd}=-1/(8\pi G)$ となる。従って曲率にのみ依存する共変作用で使う地平線エントロピーの積分 $-2\pi k_B\hbar^{-1}\int P^{abcd}\epsilon_{ab}\epsilon_{cd}\dd A$ は、この作用について面積則を再現する。ここでの計算は係数の照合である。一般の作用でこの積分が第一法則に現れることを示すには、微分同相対称性の Noether 電荷と境界条件の議論が必要であり、上の円錐計算は Einstein–Hilbert 項についてそれを直接確認する。

例えば同じスカラー部を $(R+bR^2)/(16\pi G)$ とするなら $P^{abcd}$ は上の形の $(1+2bR)$ 倍である。従って $R$ が断面上で一定のとき、この積分は $(1+2bR) k_B A_H/(4G\hbar)$ となる。Ricci 平坦な Schwarzschild では補正は零である。この比較では背景解や境界条件が作用を変えても同じと仮定する必要はない。積分の式に指定した解の地平線での $R$ を入れて評価する。

## Hawking 温度と熱力学との対応

非退化な静的地平線では $f\simeq2\kappa(r-r_+)$ である。固有距離 $\ell=\sqrt{2(r-r_+)/\kappa}$ で近傍の動径計量は $-\kappa^2\ell^2\dd t^2+\dd\ell^2$ になる。$t=-it_E$ とすると $\ell=0$ は Euclid 平面の原点であり、角度 $\kappa t_E$ の周期が $2\pi$ のとき滑らかになる。熱状態の虚時間周期が $\hbar/(k_BT)$ であるため、平衡に対応する温度は

$$
T_H=\frac{\hbar\kappa}{2\pi k_B},\qquad T_{\rm Schwarz}=\frac{\hbar}{8\pi G k_B M}.
\label{eq-gr-hawking-temperature}
$$

この周期の議論は熱平衡の温度を指定する。崩壊からの粒子生成の計算は [曲がった時空の場の量子論入門](/physics/general-relativity/qft-curved/) で導く。SI で Schwarzschild の分子は $\hbar c^3$ であり、太陽質量では約 $6.2\times10^{-8}\,\mathrm K$ である。$T_H\propto M^{-1}$ だから $c=1$ で熱容量 $\dd M/\dd T=-M/T$ は負になる。例えば質量を二倍にすると温度は半分、面積とエントロピーは四倍である。

### 放射の次数と蒸発の時刻尺度

Schwarzschild の黒体近似で温度は $T\propto M^{-1}$、面積は $A\propto M^2$ である。平衡黒体の単位面積当たりの出力が $T^4$ に比例するのは、各モードのエネルギー $\hbar\omega$ と三次元の状態密度 $\omega^2\dd\omega$ に熱占有を掛け、$x=\hbar\omega/(k_BT)$ に変数を替えるためである。従って質量減少を $\dot M=-B/M^2$ と近似すれば $M(t)^3=M_0^3-3Bt$ となる。初期質量を二倍にすると、この模型の時刻尺度は八倍である。係数 $B$ には粒子種と角モードごとの透過率が入り、質量閾値を越えると一定とは限らない。また小さな質量まで半古典近似を外挿することは蒸発の終点の理論にはならない。

### 電荷一定の熱容量が符号を変える点

Reissner–Nordström で $q$ を固定し、$r=r_+$ を変数に取る。根の式によって $m=(r+q^2/r)/2$、$r_-=q^2/r$ である。$\hbar=k_B=1$ の温度は $T=(1-q^2/r^2)/(4\pi r)$、質量は $M=m/G$ である。各導関数は
$$
\frac{\dd M}{\dd r}=\frac{1-q^2/r^2}{2G},\qquad
\frac{\dd T}{\dd r}=\frac{3q^2-r^2}{4\pi r^4}.
$$
従って固定電荷で熱容量は $C_Q=2\pi r^2(r^2-q^2)/[G(3q^2-r^2)]$ となる。非退化で $q<r<\sqrt3q$ なら正、$r>\sqrt3q$ なら負である。$r=\sqrt3q$ では温度が極大を取るため $\dd T/\dd r=0$ だが、$\dd M/\dd r\ne0$ であるため比は発散する。ここでの符号は電荷を固定した変分のものであり、電位を固定する場合の熱力学的安定性をこの一条件だけでは決められない。

:::example[同じ電荷で温度が等しくなる二つの半径]
$T$ を $0<T<T_{\max}$ に固定すると、$T(r)$ は $r=q$ で零、$r=\sqrt3q$ で極大、無限遠で再び零に近づくため、内側と外側で一つずつ解がある。例えば $r=2q$ における温度は $3/(32\pi q)$ である。この温度の他方の半径を $xq$ とすると $3x^3-8x^2+8=0$ で、既知の根 $x=2$ を除くと $3x^2-2x-4=0$ となる。正の根は $(1+\sqrt{13})/3$ で、およそ $1.54$ である。二つの半径は同じ温度を持つが熱容量の符号が違う。
:::

### 局所の温度と時刻の赤方偏移

静的領域の lapse を $N=\sqrt{-\xi^2}$ とする。静止観測者の固有時は $\dd\tau=N\dd t$ である。時刻 $t$ に関する熱平衡状態の虚時間周期を $\beta_t=\hbar/(k_BT_\infty)$ とすれば、観測者の固有時については $\beta_{\rm loc}=N\beta_t$ である。従って
$$
T_{\rm loc}N=T_\infty .
$$
これを :term[Tolman の関係]{en="Tolman relation" yomi="とるまんのかんけい"} という。同じ結論は Boltzmann 因子でも得られる。静止枠でエネルギー $E_{\rm loc}=-p\cdot\xi/N$ を持つ粒子について、保存エネルギーは $E_\infty=N E_{\rm loc}$ である。従って $e^{-E_\infty/(k_BT_\infty)}=e^{-E_{\rm loc}/(k_B(T_\infty/N))}$ となる。

非退化な地平線近くでは $N\simeq\kappa\ell$ なので、$T_\infty=T_H$ を入れると $T_{\rm loc}\simeq\hbar/(2\pi k_B\ell)$ である。同じ静止観測者の固有加速度の大きさは $|D\log N|\simeq1/\ell$ となる。例えば Schwarzschild の $r=4m$ では $N=1/\sqrt2$ なので温度は無限遠の $\sqrt2$ 倍である。地平線へ近づく静止系の温度が発散することから、自由落下系の局所応力の発散を結論してはならない。静止するための加速度も発散し、比較する運動と状態が変わる。

### 有限半径に熱浴を置いたときの安定性

壁の面積半径 $R$ を固定し、壁の固有温度を外部から指定する。$\hbar=k_B=1$ とし、$0<2m<R$ とする。Tolman の関係から $T_R=[8\pi m\sqrt{1-2m/R}]^{-1}$ である。対数を微分すると
$$
\frac{\dd\log T_R}{\dd m}
=-\frac1m+\frac{1}{R(1-2m/R)}
=\frac{3m-R}{m(R-2m)} .
$$
壁で測る熱力学的エネルギーを導くには、正則な Euclid 解の作用を先の有限半径の値で評価する。$\beta_R=\beta\sqrt f$、$\beta=8\pi m$ を使って整理すると
$$
I_E=\beta_R E_R-\frac{A_H}{4G},\qquad
E_R=\frac{R}{G}(1-\sqrt{1-2m/R}).
$$
実際、先の作用の $2R(f-\sqrt f)$ の項は $\beta_R E_R$ になり、残る $-\beta m/(2G)$ が $-A_H/(4G)$ である。$\dd E_R/\dd m=1/(G\sqrt{1-2m/R})$ で正なので壁での熱容量の符号は $3m-R$ と同じである。従って $m<R/3$ の枝では負、$R/3<m<R/2$ では正である。例えば一定の小さな $m$ に対し $R\to\infty$ にすれば負の枝に入る。壁のある場合に正の枝が現れるのは、$m$ が増えると壁に対する赤方偏移も強くなり、それが無限遠の温度減少を上回るためである。これは球対称な質量変分についての局所的安定性であり、その枝が全状態中で自由エネルギーを最小にするかは別途比較する。

### 黒体近似の定数を戻す

光子二偏光の理想黒体を考える。平衡のエネルギー密度は $2\int\dd^3k/(2\pi)^3\,\hbar c k/(e^{\hbar c k/(k_BT)}-1)$ である。被積分が方向に依存しないので、外向き半球で速度の法線成分 $c\cos\theta$ を平均すると、全方向平均は $c/4$ になる。従って流束は
$$
F=\frac{k_B^4T^4}{4\pi^2 c^2\hbar^3}\int_0^\infty\frac{x^3}{e^x-1}\dd x
=\frac{\pi^2 k_B^4}{60c^2\hbar^3}T^4 .
$$
最後の積分は $e^{-nx}$ の和に展開して $\sum 6/n^4=\pi^4/15$ を使った。面積を $A_H=16\pi G^2 M^2/c^4$ と同一視し、$T_H=\hbar c^3/(8\pi G k_B M)$ を入れると $P=A_HF=\hbar c^6/(15360\pi G^2M^2)$ である。質量減少は $c^2\dot M=-P$ だから、この模型を $M=0$ まで外挿した時刻は
$$
t_{\rm model}=\frac{5120\pi G^2M_0^3}{\hbar c^4}.
$$
これは面積を放射面として用いた黒体近似である。実際には遠方へ伝播するモードの透過率と放射する粒子種が係数を変える。例えば同じ光子だけの模型で $M_0$ を $1/10$ にすると、温度は十倍になり、出力は百倍、時刻尺度は千分の一になる。この比も、粒子種の質量閾値を越えて係数が変わればそのまま使えない。

:::example[温度と質量の変化を同時に追う]
$\dot M=-B/M^2$ が一定の $B>0$ で成り立つと仮定する。初期質量の半分に達する時刻は $(M_0^3-(M_0/2)^3)/(3B)=7M_0^3/(24B)$ である。模型の全時刻の $7/8$ が既に経過しており、その時点で温度は初期の二倍、面積は四分の一である。従って残りの時刻は全時刻の $1/8$ であり、質量の残りの割合とは違う。初期の質量が十分大きく、この半分まで半古典近似が成立する間なら、模型の終点を観測したとみなす必要もない。
:::

## 情報喪失問題の概観

:term[情報喪失問題]{en="black hole information problem" yomi="じょうほうそうしつもんだい"} は、崩壊の純粋状態が完全な蒸発の後にどう記述されるかを問う。半古典近似の出射モードと穴の内側の相方を二つの系 $A,B$ とし、単一周波数の模型を $\ket\Psi=\sqrt{1-z}\sum_{n\ge0}z^{n/2}\ket n_A\ket n_B$、$z=e^{-\hbar\omega/(k_BT_H)}$ とする。内側の跡を取ると $\rho_A=(1-z)\sum z^n\ket n\bra n$ であり、純度は $(1-z)^2/(1-z^2)=(1-z)/(1+z)<1$ である。全体は純粋でも放射部分だけは混合である。もし蒸発後に相方や相関が一切残らないとすれば、閉じた系のユニタリ時間発展と整合しない。

放射全体が純粋になるユニタリな過程では、放射の細粒度エントロピーは終盤に減少する必要がある。これを表す曲線を Page 曲線と呼ぶ。局所的な熱スペクトルは弱い相関の有無を決めないため、熱分布から最終状態の純度は推定できない。量子重力でどの自由度に相関が保持され、どの観測量に現れるかは近似の枠を越える。ここでの二モード計算は蒸発の終点の記述を与えない。

### 二モードの細粒度エントロピー

外側の対角重みを $p_n=(1-z)z^n$ とすると、平均数は $\sum n p_n=z/(1-z)$ である。従って von Neumann エントロピーは

$$
-\tr(\rho_A\log\rho_A)=-\log(1-z)-\frac{z}{1-z}\log z .
$$

これは $k_B$ を単位として表した値である。例えば $z=1/2$ なら平均数は $1$、エントロピーは $2\log2$ である。一方、同じ状態の $A,B$ を合わせた密度演算子は階数一なのでエントロピーは零である。$A$ のエントロピーは熱状態の重みだけから計算できるが、純粋な全体系のどこに相方があるかを決めるには相関も測る必要がある。

### 局所の熱分布から全体の純度は決まらない

二つの二準位系で $\ket\psi=\sqrt{1-p}\ket{00}+\sqrt p\ket{11}$ を取る。密度演算子 $\rho=\ket\psi\bra\psi$ の部分跡では $\ket{00}\bra{11}$ とその共役が消え、各側では $\diag(1-p,p)$ となる。一方、混合状態 $(1-p)\ket{00}\bra{00}+p\ket{11}\bra{11}$ も同じ部分跡を持つ。前者の全体の純度は $1$、後者は $(1-p)^2+p^2$ であり、例えば $p=1/4$ では後者は $5/8$ である。この二つを区別するには、例えば $\ket{00}\bra{11}+\ket{11}\bra{00}$ の期待値を測る。前者では $2\sqrt{p(1-p)}$、後者では零である。個別に測った熱的な占有数が同じでも、全体に相関が残っている場合がある。

### Page 曲線が現れる有限系の計算

全体系を $A\otimes B$、次元を $d_A,d_B$ とし、全体は純粋状態とする。Schmidt 係数を $\sqrt{p_i}$ とすれば、部分系 $A$ の固有値は $p_i$、非零な個数 $r$ は $\min(d_A,d_B)$ 以下である。全体を基底で行列に並べたときの特異値分解が Schmidt 分解であり、部分跡で直交する相方の積だけが残る。従って $S_A=S_B=-\sum_{i=1}^r p_i\log p_i$ となる。この和は $\log r$ 以下である。実際、$\log$ の凹性から $\sum p_i\log(1/p_i)\le\log\sum p_i(1/p_i)=\log r$ である。等号は正の $p_i$ が等しいときである。例えば $d_A=8,d_B=2$ の時点で、放射を $A$ とすると純粋な全体系での放射のエントロピーは高々 $\log 2$ となる。初期に $d_A=1$、終期に $d_B=1$ とする有限模型では、両端で $S_A=0$ であり、途中では小さい方の次元で上限が決まる。この計算は Page 曲線の上限を示す。実際の蒸発でその上限付近の値を取るかは、時間発展や相関を与える特定の理論で判定する。

### 放射の粗視化では捨てる相関を計算する

:::example[純粋な二体系の相互情報量]
二つの二準位系で $\ket\psi=(\ket{00}+\ket{11})/\sqrt2$ を取る。全体系の細粒度エントロピーは零、部分系は双方とも $\diag(1/2,1/2)$ なので各エントロピーは $\log2$ である。相互情報量 $I(A:B)=S_A+S_B-S_{AB}$ は $2\log2$ になる。混合の $\tfrac12\ket{00}\bra{00}+\tfrac12\ket{11}\bra{11}$ に替えると、各部分系の値は同じだが全体系のエントロピーが $\log2$ なので相互情報量は $\log2$ になる。占有数の分布だけでは二つの状態を区別できず、放射全体の細粒度の議論には異なるモード間の相関が要る。
:::

### 無作為な純粋状態で粗視化と細粒度を比べる

$A,B$ の次元を $a,b$ とし、積の空間上で一様な単位ベクトル $\ket\psi$ を選ぶ。$\rho=\ket\psi\bra\psi$ の二コピーの平均は $I$ と交換作用素 $F_{AB}$ を用いて
$$
\mathbb E(\rho\otimes\rho)=\frac{I+F_{AB}}{ab(ab+1)}
$$
となる。これを成分で確かめる。$d=ab$ として $\psi_i$ の独立な位相の不変性を使うと、$\mathbb E\,\psi_i\psi_j\overline{\psi_k\psi_l}$ は $(i,j)$ と $(k,l)$ の添字の組が同じ場合だけ非零である。置換不変性により $i\ne j$ の組の係数を $v=\mathbb E|\psi_i|^2|\psi_j|^2$、同じ添字の四乗を $u=\mathbb E|\psi_i|^4$ と書ける。最初の二成分を $(\psi_1+\psi_2)/\sqrt2$ に変えた四乗平均も $u$ であり、位相の平均で残るのは $(2u+4v)/4$ だから $u=2v$ である。規格化の二乗から $1=d u+d(d-1)v=d(d+1)v$ となり、この式を得る。$d=1$ でも式を直接確認できる。

$A$ の二コピーだけを交換する作用素を $F_A$ と置くと、$\tr(\rho_A^2)=\tr[(\rho\otimes\rho)F_A]$ である。実際、積基底で $A$ の添字を入れ替えて跡を取ると、$B$ の二つの添字はそれぞれ独立に足し合わされ、$A$ の二つの部分跡の積になる。$\tr F_A=ab^2$、$\tr(F_{AB}F_A)=\tr F_B=a^2b$ だから
$$
\mathbb E\,\tr(\rho_A^2)=\frac{a+b}{ab+1}.
$$
従って $a\ll b$ では純度は平均でほぼ $1/a$ となり、逆の $b\ll a$ ではほぼ $1/b$ となる。これは純度から求めた二次 Rényi エントロピー $s_2=-\log\tr\rho_A^2$ の尺度を与える。Jensen の不等式で $\mathbb E s_2\ge-\log[(a+b)/(ab+1)]$、また各状態で $s_2\le s_A\le\log\min(a,b)$ である。$s_2\le s_A$ は固有値 $p_i$ を重みとして $-\log$ の凸性を適用すれば従う。例えば $a=8,b=2$ なら平均純度は $10/17$ であり、平均の $s_A$ は $\log(17/10)$ 以上 $\log2$ 以下である。無作為な状態は時間発展についての仮定であり、半古典的な各モードの計算から正当化されるものではない。

### 放射の独立な対生成を繰り返したとき

早い放射を $R$、新しい出射モードを $B$、内部の相方を $C$ とする。新しい対 $BC$ が純粋で、既存の放射から独立だと仮定する。既存の内部自由度を跡で消すと $\rho_{RB}=\rho_R\otimes\rho_B$ である。固有値は各側の固有値の積だから、$s(RB)-s(R)=s(B)\ge0$ となる。新しい対に非自明な絡み合いがあればこの差は正であり、完全蒸発をユニタリな過程とみなしたときに終盤で必要となる減少は起こらない。

対の独立性を課さない場合、相互情報量 $I(R:B)=s(R)+s(B)-s(RB)$ から $s(RB)-s(R)=s(B)-I(R:B)$ になる。従って減少には $I(R:B)>s(B)$ が必要である。例えば $R,B$ が二準位の最大絡み合い対なら、$B$ 単独の固有値は $(1/2,1/2)$ のまま、$s(RB)=0$ で $s(RB)-s(R)=-\log2$ となる。この模型では $B$ が既に $R$ と純粋な対を作るので、独立な別の $C$ と同時に同じ純粋な対を作れない。実際、$RB$ の状態が階数一なら、任意の $RBC$ 状態はその一次元部分空間に台を持ち、$RB$ と $C$ のテンソル積になる。観測した局所スペクトルだけで独立な対生成を仮定すれば、放射のエントロピー変化を誤る。

:::example[初期状態の違いを放射の位相に残す]
二つの直交する初期状態がそれぞれ $\ket{\psi_\pm}=(\ket{00}\pm\ket{11})/\sqrt2$ という最終放射へ写る有限模型を取る。この写像は内積を保つので、この部分空間の情報を失わない。各出射モード単独の状態はどちらも $I/2$ であり、占有数も同じである。一方、$X\otimes X$（$X\ket0=\ket1$、$X\ket1=\ket0$）の期待値はそれぞれ $+1,-1$ となる。初期情報を区別するためにどの相関が必要かが、この例では明示できる。
:::

## 演習問題

::::exercise{level=basic}
Schwarzschild の質量が三倍なら $\kappa$ は何倍か。

:::solution
$\kappa=1/(4GM)$ だから $1/3$ 倍である。
:::
::::

::::exercise{level=basic}
$M$ が二倍の Schwarzschild のエントロピーの比を求めよ。

:::solution
$A_H=16\pi G^2M^2$ と :eqref[eq-gr-bh-entropy] により四倍である。
:::
::::

::::exercise{level=basic}
二モード模型で $z=1/2$ の外側の純度を求めよ。

:::solution
$\tr\rho_A^2=(1-z)/(1+z)=(1/2)/(3/2)=1/3$ である。
:::
::::

::::exercise{level=basic}
Reissner–Nordström の $r_+=3L,r_-=L$ の $\kappa$ を求めよ。

:::solution
$a=0$ なので $\kappa=(r_+-r_-)/(2r_+^2)=1/(9L)$ である。
:::
::::

::::exercise{level=advanced}
Schwarzschild について第一法則の面積の係数を直接確かめよ。

:::solution
$\dd A_H=32\pi G^2M\,\dd M$、$\kappa=1/(4GM)$ なので $(\kappa/(8\pi G))\dd A_H=\dd M$ である。
:::
::::

::::exercise{level=advanced}
定常 Schwarzschild の不可逆質量を $M_{\rm irr}=\sqrt{A_H/(16\pi G^2)}$ で定める。同質量の極限 Kerr と比較せよ。

:::solution
Schwarzschild では $M_{\rm irr}=M$、極限 Kerr では $A_H=8\pi G^2M^2$ なので $M_{\rm irr}=M/\sqrt2$ である。
:::
::::

::::exercise{level=advanced}
遠方で時計の規格化を変えて $\chi$ を $b\chi$、$b>0$ とすると $\kappa$ がどう変わるか示せ。

:::solution
左辺は $b^2\nabla_\chi\chi=b^2\kappa\chi$ である。右辺を新しい生成子で表せば $b\kappa(b\chi)$ なので新係数は $b\kappa$ である。
:::
::::

::::exercise{level=advanced}
光子の黒体放射を $P\propto A_HT_H^4$ と近似する。Schwarzschild の蒸発時間の質量依存を求めよ。

:::solution
$A_H\propto M^2$、$T_H^4\propto M^{-4}$ なので $\dot M=-C/M^2$ となる。積分すると $M^3=M_0^3-3Ct$ であり、外挿の寿命は $M_0^3/(3C)$ である。
:::
::::

::::exercise{level=advanced}
質量比 $3:4$ の二つの非回転ブラックホールが十分遠くから合体する。古典面積定理で制約される放射のエネルギーの上限を、初期総質量エネルギーに対する割合で求めよ。最終穴も回転しないとする。

:::solution
質量を $3M_0,4M_0$ と置く。面積非減少から $M_f^2\ge9M_0^2+16M_0^2=25M_0^2$、従って $M_f\ge5M_0$ である。初期質量は結合エネルギーを無視して $7M_0$ であるから、放射は高々 $2M_0$、総エネルギーで割って $2/7$ である。この比は実際の運動の効率ではなく、終状態の面積だけで得る上限である。最終穴が回転する場合は同じ質量で面積が小さくなるため、この質量下限は依然必要である。
:::
::::

## 確認テスト

::::quiz
古典的な面積非減少の証明で使う条件を選べ。

- [x] 零収束条件と未来生成子の大域条件
- [ ] すべての静止観測者の固有加速度が有限
- [ ] 放射の熱スペクトル

:::explanation
零収束条件で負の膨張は有限径数の焦点を生じる。未来の因果境界を生成し続けるという大域条件に反するため、膨張は非負になる。
:::
::::

::::quiz{answer="0.25" tol="1e-9"}
Schwarzschild の $m\kappa$ を答えよ。

:::explanation
入射座標で $\kappa=f'(2m)/2$、$f'=2m/r^2$ なので、$\kappa=1/(4m)$ である。従って $m\kappa=1/4$ となる。
:::
::::

::::quiz
二モード模型で、外側だけの状態が混合になる操作を選べ。

- [x] 内側の状態について部分跡を取る。
- [ ] 全体系に同じユニタリを施す。
- [ ] 時間並進の Killing ベクトルを規格化する。

:::explanation
$\ket n_B$ の直交性により部分跡は $n$ が等しい項だけを残す。複数の $n$ の重みが非零なら外側は混合状態となる。全体の純度の変化を意味しない。
:::
::::

::::quiz{answer="0.5" tol="1e-9"}
Schwarzschild の質量が二倍になったときの Hawking 温度の比を求めよ。

:::explanation
:eqref[eq-gr-hawking-temperature] により $T_H\propto M^{-1}$ なので、新しい温度は半分である。
:::
::::

::::quiz{answer="0.707106781187" tol="1e-7"}
等質量で初期の回転が零の二つの穴が合体した。面積定理から全初期質量に対する最終質量比の下限を数値で求めよ。

:::explanation
面積の非減少は $M_f^2\ge M_1^2+M_2^2=2M_1^2$ を与える。合計 $2M_1$ で割れば最終質量比は少なくとも $1/\sqrt2$ である。
:::
::::

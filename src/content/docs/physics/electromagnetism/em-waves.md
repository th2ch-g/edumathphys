---
title: 電磁波
description: 真空波の偏光から出発し、物質の応答が屈折・吸収・反射を決める仕組みと、金属壁で許される伝播・共振モードを導く。
notation:
  - symbol: '\vb{k}'
    meaning: 波数ベクトル
  - symbol: 'Z'
    meaning: 波動インピーダンス
  - symbol: '\tilde n'
    meaning: 複素屈折率
  - symbol: '\beta'
    meaning: 導波管の軸方向波数
  - symbol: 'k_c'
    meaning: 遮断波数
---

## この章の目標

- 横波の場、強度、偏光を複素振幅から復元する。
- Lorentz 振動子から分散と吸収を導く。
- 偏光別の反射率・透過率、Brewster 角、全反射を計算する。
- 導波管と空洞で境界条件から固有周波数を求める。

## 真空中の平面波

### 場を一組の振幅で表す

真空で電荷も電流もない領域を考える。
実在する場を $\vb{E}=\Re[\vb{E}_0 e^{i(\vb{k}\cdot\vb{r}-\omega t)}]$ と書き、磁場も同じ位相因子で表す。
以下で $\omega>0$ である。
微分は $\nabla\to i\vb{k}$、$\partial_t\to-i\omega$ に置き換わるので :eqref[eq-em-maxwell] から

$$
\vb{k}\cdot\vb{E}_0=\vb{k}\cdot\vb{B}_0=0,\qquad \vb{B}_0=\frac{\vb{k}\times\vb{E}_0}{\omega},\qquad \vb{k}\times\vb{B}_0=-\frac{\omega}{c^2}\vb{E}_0.
\label{eq-em-plane-wave}
$$

後二式から $\vb{k}\times(\vb{k}\times\vb{E}_0)=-k^2\vb{E}_0$ を用いて $\omega=ck$ を得る。
電場と磁場の大きさは $E=cB$ であり、$\vb{E},\vb{B},\vb{k}$ は右手系をつくる。

:term[波動インピーダンス]{en="wave impedance" yomi="はどういんぴーだんす"} は一方向進行波の $E/H$ である。
真空では

$$
Z_0=\sqrt{\mu_0/\varepsilon_0}\simeq376.73\,\Omega
$$

である。

単色波どうしの積の時間平均は

$$
\langle\Re(Ae^{-i\omega t})\Re(Be^{-i\omega t})\rangle=\Re(AB^*)/2
$$

である。
従って

$$
\langle\vb{S}\rangle=\tfrac12\Re(\vb{E}_0\times\vb{H}_0^*)=\frac{\vb{E}_0\cdot\vb{E}_0^*}{2Z_0}\hat{\vb{k}}.
\label{eq-em-wave-intensity}
$$

実効値を使うなら $1/2$ は不要である。

:::example[強度から電場を求める]
強度 $1.00\,\mathrm{kW/m^2}$ の直線偏光はピーク電場 $\sqrt{2Z_0 I}=868\,\mathrm{V/m}$、ピーク磁場 $E_0/c=2.90\,\mathrm{\mu T}$ をもつ。
エネルギー密度の平均は

$$
I/c=3.34\times10^{-6}\,\mathrm{J/m^3}
$$

である。
:::

:::example[定在波の節]
互いに逆向きに進む等振幅の波を重ねると、電場は $2E_0\cos kz\cos\omega t\,\hat{\vb{x}}$、磁場は $2(E_0/c)\sin kz\sin\omega t\,\hat{\vb{y}}$ になる。
電場の節は磁場の腹である。

瞬時にはエネルギーが行き来するが、時間平均の流れは零である。
:::

### 有限時間のパルスが運ぶ量

平面波に有限の持続時間をもたせると、単位面積を通過する全エネルギーを計算できる。
$E_x(t,z)=E_p\exp[-(t-z/c)^2/(2\tau^2)]$、$B_y=E_x/c$ とする。
一点を波が通過するとき、時間積分した強度は

$$
W/A=\frac{E_p^2}{Z_0}\int_{-\infty}^{\infty} e^{-(t-z/c)^2/\tau^2}\dd t=\frac{\sqrt\pi\tau E_p^2}{Z_0}.
$$

この式にさらに $1/2$ を掛ける必要はない。
ここで与えたのは搬送波の包絡線ではなく、瞬時の実電場そのものだからである。
例えば $E_p=100\,\mathrm{V/m}$、$\tau=1.0\,\mathrm{ns}$ なら

$$
W/A=4.70\times10^{-8}\,\mathrm{J/m^2}
$$

である。

波束全体の運動量を面積で割った値は $W/(Ac)$ となる。

### 有限幅の光束を平面波の重ね合わせで表す

真空中の単色場に時間因子 $e^{-i\omega t}$ を選ぶ。$z=0$ での横方向分布を Fourier 分解すると、各成分は $k_z^2=k_0^2-k_x^2-k_y^2$ に従う。$k_0=\omega/c$ である。正の $z$ へ伝わる解では $k_z$ の実部を正に、純虚数のときは虚部を正に取る。従って横波数が $k_0$ を超える細かい模様は、真空でも界面から離れるにつれて指数的に減衰する。これは物質による吸収を仮定した結果ではない。縦電場は $k_z\ne0$ の成分ごとに $E_z=-(k_x E_x+k_y E_y)/k_z$ で定まり、単に横成分をそのまま全空間へ運んだだけでは Maxwell の拘束式を満たさない。

$k_\perp/k_0\ll1$ のとき $k_z=k_0-k_\perp^2/(2k_0)+O(k_\perp^4/k_0^3)$ である。電場から $e^{ik_0z}$ をくくった包絡線を $u$ とすると、Fourier 成分の微分 $\partial_z\tilde u=-ik_\perp^2\tilde u/(2k_0)$ を実空間へ戻して $2ik_0\partial_z u+\nabla_\perp^2u=0$ を得る。例えば $u(x,0)=u_0\exp(-x^2/w_0^2)$ という一方向の Gaussian 分布なら、Fourier 空間での二次の指数をまとめて逆変換すれば

$$
u(x,z)=\frac{u_0}{\sqrt{1+iz/z_R}}
\exp\left[-\frac{x^2}{w_0^2(1+iz/z_R)}\right],\qquad
z_R=k_0 w_0^2/2.
$$

二乗振幅の幅は $w(z)=w_0\sqrt{1+(z/z_R)^2}$ である。光束が細いほど横波数の幅が広くなるため、短い距離で広がる。$w_0=0.500\,\mathrm{mm}$、真空波長 $500\,\mathrm{nm}$ なら $z_R=1.57\,\mathrm m$ である。この例は一方向だけで集束する光を扱う。両方向を同じ幅で集束させれば二つの独立な Gaussian 積分が掛かり、軸上振幅の分母は平方根から一次に変わる。

## 偏光

### 位相差を測る

$+z$ に進む波の横成分を $(E_x,E_y)=(A\cos\psi,B\cos(\psi+\delta))$、$\psi=kz-\omega t$ とする。
$\sin\psi$ を消去すると

$$
(E_x/A)^2+(E_y/B)^2-2(E_x/A)(E_y/B)\cos\delta=\sin^2\delta.
$$

従って固定点で電場の端は楕円を描く。
$\delta=0,\pi$ なら直線、$A=B$ かつ $\delta=\pm\pi/2$ なら円である。
回転の左右という名称は観察方向で入れ替わるので、ここでは振幅と位相で指定する。

完全に偏光した単色波は :term[Jones ベクトル]{en="Jones vector" yomi="じょーんずべくとる"} $(A,B e^{i\delta})^\top$ で表せる。
偏光素子を直交二成分の線形変換とみなせば、単位透過軸 $\vb{e}$ の理想直線偏光子は $\vb{e}\vb{e}^\top$ である。

角 $\alpha$ の直線偏光に対し振幅比は $\cos\alpha$ なので、強度は $I_0\cos^2\alpha$ となる。

:::example[四分の一波長板]
二つの主軸に沿う屈折率を $n_x,n_y$ とすると、厚さ $d$ で生じる位相差は

$$
\delta=2\pi(n_y-n_x)d/\lambda_0
$$

である。
$\lambda_0=600\,\mathrm{nm}$、$|n_y-n_x|=0.010$ なら最小厚は $15.0\,\mathrm{\mu m}$ となる。

主軸に $45^\circ$ の入力は等振幅なので円偏光に変わる。
主軸に平行な入力では偏光は変わらない。
:::

### 部分偏光と Stokes パラメータ

互いに位相が揺らぐ波は一つの Jones ベクトルでは記述できない。
横成分の複素振幅の二次平均を $J_{ij}=\langle E_i E_j^*\rangle$ とし、$S_0=J_{xx}+J_{yy}$、$S_1=J_{xx}-J_{yy}$、$S_2=2\Re J_{yx}$、$S_3=2\Im J_{yx}$ と定める。

これは本章での $S_3$ の向きの規約である。
$J$ は非負な Hermite 行列なので

$$
4\det J=S_0^2-S_1^2-S_2^2-S_3^2\ge0
$$

である。
完全偏光は $J$ の階数が $1$ で等号、非偏光は $J$ が単位行列の定数倍で $S_1=S_2=S_3=0$ となる。

:::example[偏光子を回しても区別できない二つの波]
円偏光と非偏光はどちらも理想直線偏光子で強度の半分を通す。
円偏光に適切な四分の一波長板を先に通せば直線偏光となり、後段の偏光子を回すと消光点が現れる。
非偏光の $J$ はユニタリ変換で不変なので、同じ操作でも消光点は生じない。
:::

### 偏光状態を四回の強度測定から求める

理想偏光子を通した強度は、共通係数 $1/(2Z)$ を除けば $\vb{e}^\dagger J\vb{e}$ である。
直交する二軸で測った強度をこの係数で割り $I_x,I_y$ とすれば $S_0=I_x+I_y,\ S_1=I_x-I_y$ である。

$45^\circ$ の直線偏光子なら $I_{45}=(S_0+S_2)/2$ なので $S_2=2I_{45}-S_0$ を得る。
最後に $\vb{e}=(1,i)^\top/\sqrt2$ だけを透過させる円偏光分析器を使えば $I_+=(S_0+S_3)/2$ である。

これで $S_3$ の符号を観測する規約まで定まる。
例えば $(I_x,I_y,I_{45},I_+)=(3,1,2,3)$ なら $(S_0,S_1,S_2,S_3)=(4,2,0,2)$、偏光度は

$$
\sqrt{S_1^2+S_2^2+S_3^2}/S_0=1/\sqrt2
$$

となる。

### 波長板を回して偏光の楕円を調べる

吸収しない複屈折板の主軸を $1,2$ とする。垂直入射で界面反射を無視し、板を通過する間の相対位相を $\delta=2\pi(n_2-n_1)d/\lambda_0$ とする。主軸座標での Jones 行列は全体位相を除いて $\operatorname{diag}(1,e^{i\delta})$ である。第一軸が $x$ 軸に対して $\alpha$ をなすなら、入射 $x$ 偏光の主軸成分は $(\cos\alpha,-\sin\alpha)$ である。第二軸への位相を加え、実験室座標へ戻すと

$$
e_x=\cos^2\alpha+e^{i\delta}\sin^2\alpha,\qquad
e_y=(1-e^{i\delta})\sin\alpha\cos\alpha .
$$

これにより交差した検光子を通る強度は $I_y/I_i=\sin^2(2\alpha)\sin^2(\delta/2)$ となる。半波長板では出射偏光が直線のまま $2\alpha$ 回転する。四分の一波長板で $\alpha=45^\circ$ なら、主軸成分が等振幅で位相差が直角になるため円偏光である。強度 $I_y$ が同じでも円偏光と直線偏光を区別できない場合があるので、回転方向を判定するには別の分析器または位相測定が必要となる。


:::example[回転検光子で直線偏光の向きを測る]
$S_1/S_0=3/5,S_2/S_0=4/5,S_3=0$ の光を考える。偏光度は一で、二成分の相対位相は零だから直線偏光である。線形検光子の向きを $\alpha$ とすれば射影強度は $[S_0+S_1\cos2\alpha+S_2\sin2\alpha]/2$ である。

$x=\alpha/\pi$、$Y_*=10.0\,\mathrm{W/m^2}$ を入射強度とし、$Y$ を透過強度とする。

$$
f(x)=\frac{4 \sin{\left(2 \pi x \right)} + 3 \cos{\left(2 \pi x \right)} + 5}{10},\qquad Y=Y_* f(x).
$$

変数を少し動かしたときの読みの変化は、接線の傾きから計算できる。基準値 $x_0$ の近くで

$$
\begin{aligned}
\frac{\Delta Y}{Y_*}= f\!\left(x_0+\Delta x\right)-f(x_0)= f'(x_0)\Delta x+O((\Delta x)^2),\\
f'(x)=\frac{\pi \left(- 3 \sin{\left(2 \pi x \right)} + 4 \cos{\left(2 \pi x \right)}\right)}{5}.
\end{aligned}
$$

同じ尺度で値と傾きを並べる。

| $x$ | $f(x)$ | $Y$ $[\mathrm{W/m^2}]$ | $\dd Y/\dd x$ $[\mathrm{W/m^2}]$ |
|---:|---:|---:|---:|
| $0$ | $0.8$ | $8$ | $25.13$ |
| $\frac{1}{8}$ | $0.995$ | $9.95$ | $4.443$ |
| $\frac{1}{4}$ | $0.9$ | $9$ | $-18.85$ |
| $\frac{1}{2}$ | $0.2$ | $2$ | $-25.13$ |

最大透過の向きは $2\alpha=\arctan(4/3)$ であり、強度軸の回転角は Stokes 成分から得た角の半分である。分析器を半回転すれば元の読みへ戻る。
:::

## 物質中の電磁波と分散

### 物質中の分散関係

均一・等方・線形な媒質で

$$
\vb{D}=\varepsilon(\omega)\vb{E},\ \vb{B}=\mu(\omega)\vb{H}
$$

とする。
自由電荷・電流がなければ、真空の導出で $\varepsilon_0,\mu_0$ を置き換えて $k^2=\omega^2\mu\varepsilon$、$Z=\omega\mu/k$ となる。

非磁性媒質では $k=(\omega/c)\tilde n$、$\tilde n^2=\varepsilon_r$ である。
受動な媒質へ $+z$ に進む枝は $\Im k\ge0$ を選ぶ。
振幅は $e^{-\Im k\,z}$、強度は $e^{-2\Im k\,z}$ で減衰する。

### Lorentz 振動子

束縛された電子の変位 $x$ に線形復元力と摩擦を仮定する。
磁気力や隣の分子の局所場を無視し、

$$
m(\ddot x+\gamma\dot x+\omega_0^2x)=-eE
$$

とする。
定常応答 $x=x_0 e^{-i\omega t}$ と $P=-Nex$ を代入すれば

$$
\varepsilon_r(\omega)=1+\frac{\omega_p^2}{\omega_0^2-\omega^2-i\gamma\omega},\qquad \omega_p^2=\frac{Ne^2}{m\varepsilon_0}.
\label{eq-em-lorentz-dispersion}
$$

実部は $\omega_0^2-\omega^2$ とともに符号を変える寄与をもち、虚部は $\omega_p^2\gamma\omega/[(\omega_0^2-\omega^2)^2+\gamma^2\omega^2]$ で正である。
平均吸収は

$$
\langle E\dot P\rangle=\omega\varepsilon_0\Im(\varepsilon_r)|E_0|^2/2
$$

となる。

振動子の式に $N\dot x$ を掛ければ、電場の仕事は運動・ばねエネルギーの増加と $Nm\gamma\dot x^2$ の散逸に分かれる。
これで虚部の符号を検算できる。

$\omega_0=0$ は自由電子の Drude 模型であり、

$$
\sigma(\omega)=Ne^2/[m(\gamma-i\omega)]
$$

となる。
無衝突では $\varepsilon_r=1-\omega_p^2/\omega^2$ なので、$\omega<\omega_p$ は進行波を許さない。

Ohm 則 $\vb{J}=\sigma\vb{E}$ を直接使う場合は $\varepsilon_{\rm eff}=\varepsilon+i\sigma/\omega$ にまとめられる。
良導体 $\sigma\gg\omega\varepsilon$ では $k\simeq(1+i)\sqrt{\mu\sigma\omega/2}$ となり、振幅の :term[表皮深さ]{en="skin depth" yomi="ひょうひふかさ"} は

$$
\delta_s=\sqrt{2/(\mu\sigma\omega)}
$$

である。

:::example[銅の表皮深さ]
$\sigma=5.8\times10^7\,\mathrm{S/m}$、$\mu\simeq\mu_0$ を一定とし、$f=1.00\,\mathrm{MHz}$ とすると $\delta_s=66.1\,\mathrm{\mu m}$ となる。

厚さが $\delta_s$ の層を通ると場の振幅が $1/e$、強度が $1/e^2$ になる。
:::

### 波束の速度

透明域で狭い波数帯の重ね合わせをとる。
$\omega(k)=\omega_*+(k-k_*)\omega'(k_*)+\cdots$ を Fourier 積分の位相に入れると、搬送波に掛かる包絡線は $z-\omega'(k_*)t$ の関数になる。

従って位相速度は $v_p=\omega/k$、群速度は

$$
v_g=\dd\omega/\dd k=c/[n+\omega(\dd n/\dd\omega)]
$$

である。
二次項はパルスを広げる。
強い吸収や異常分散では波形が大きく変わるため、包絡線のピーク速度から情報の先端速度は求められない。

吸収と屈折を独立に選べるかどうかは、因果性にもかかわる。
以下では応答が十分速く減衰する場合にこの関係を導く。
:cite[jackson-em]

:::example[プラズマ中の透明帯]
無衝突の自由電子なら $\omega^2=\omega_p^2+c^2k^2$ である。
$v_p=\omega/k>c$ だが $v_g=c^2k/\omega<c$ で、積は $c^2$ である。
$\omega=2\omega_p$ では $v_g/c=\sqrt3/2$ となる。
:::

### 因果性が分散と吸収を結ぶ

因果的な線形応答を $P(t)=\varepsilon_0\int_0^\infty h(s)E(t-s)\dd s$ と書く。
時間因子 $e^{-i\omega t}$ に対し $\chi(z)=\int_0^\infty h(s)e^{izs}\dd s$ は $\Im z>0$ で解析的である。

実の $h$ に対して $\chi(-\omega)=\chi(\omega)^*$ である。
ここでは $h$ が十分速く減衰し、$\chi(z)$ は上半円で無限遠に向かって一様に零へ近づくと仮定する。
瞬時応答の定数項があればそれを引いて議論する。

上半平面の $z$ に対して Cauchy の積分公式を使い、半円の寄与を消すと

$$
\chi(z)=(2\pi i)^{-1}\int_{-\infty}^\infty\chi(s)/(s-z)\dd s
$$

である。
$z=\omega+i0$ として $1/(s-\omega-i0)=\operatorname{PV}(1/(s-\omega))+i\pi\delta(s-\omega)$ を使う。

後者は幅 $\eta$ の核 $\eta/[(s-\omega)^2+\eta^2]$ の面積が $\pi$ で、実部が対称主値に収束することから従う。
左辺の半分を移項し、実部と虚部を分けると

$$
\begin{aligned}
\Re\chi(\omega)&=\frac1\pi\operatorname{PV}\int_{-\infty}^\infty\frac{\Im\chi(s)}{s-\omega}\dd s,\\
\Im\chi(\omega)&=-\frac1\pi\operatorname{PV}\int_{-\infty}^\infty\frac{\Re\chi(s)}{s-\omega}\dd s.
\end{aligned}
\label{eq-em-kramers-kronig}
$$

これが Kramers–Kronig 関係である。
例えば静的な感受率は

$$
\Re\chi(0)=(2/\pi)\int_0^\infty \Im\chi(s)\dd s/s
$$

であり、広い周波数域の吸収スペクトルから得られる。
単一周波数の吸収だけでは決められない。

:::example[緩和型媒質で静的感受率を復元する]
$h(s)=(\chi_s/\tau)e^{-s/\tau}$ とすると

$$
\chi(\omega)=\chi_s/(1-i\omega\tau)
$$

である。
損失部は $\chi_s\omega\tau/(1+\omega^2\tau^2)$ となる。

静的な関係に入れると

$$
(2\chi_s/\pi)\int_0^\infty\dd u/(1+u^2)=\chi_s
$$

である。
Lorentz 振動子とは異なる応答でも同じ因果性の条件を満たす。
:::

### 透明だが分散する媒質に蓄えられるエネルギー

非磁性で $\gamma=0$ の Lorentz 振動子一種を考え、共鳴点を除く。
真空の電場エネルギーに束縛電子の運動エネルギーとばねエネルギーを足して、物質中での仕事の行き先を追う。
ピーク振幅 $x_0=-e E_0/[m(\omega_0^2-\omega^2)]$ を用いると、その時間平均は

$$
\overline w_{\rm osc}=\frac{Nm}{4}(\omega^2+\omega_0^2)|x_0|^2.
$$

一方、$\varepsilon=\varepsilon_0[1+\omega_p^2/(\omega_0^2-\omega^2)]$ を微分すると

$$
\dd(\omega\varepsilon)/\dd\omega=\varepsilon_0[1+\omega_p^2(\omega_0^2+\omega^2)/(\omega_0^2-\omega^2)^2]
$$

である。

真空電場と振動子の和をこの形にまとめれば

$$
\overline w=\tfrac14\frac{\dd(\omega\varepsilon)}{\dd\omega}|E_0|^2+\tfrac14\mu_0|H_0|^2.
$$

$\varepsilon$ だけを使う非分散媒質の式では、共鳴近くに蓄えられる振動エネルギーを取り落とす。
透明な進行帯で $\varepsilon=\varepsilon_0n^2$、$H_0=nE_0/Z_0$ を代入すると

$$
\overline w=\varepsilon_0 n[n+\omega\,\dd n/\dd\omega]|E_0|^2/2
$$

となる。

時間平均の流れは $\overline S=\varepsilon_0 c n|E_0|^2/2$ なので、

$$
\overline S/\overline w=c/(n+\omega\,\dd n/\dd\omega)=v_g
$$

である。
損失の大きい媒体ではこの無減衰模型のエネルギー式を使わない。

### 吸収スペクトルの重みと振動子強度

因果的な感受率は虚周波数 $z=i\xi$、$\xi>0$ でも評価できる。Kramers–Kronig 関係の積分表示で正負の実周波数を組み合わせると、$\Im\chi$ の奇性から

$$
\chi(i\xi)=\frac{2}{\pi}\int_0^\infty\frac{s\Im\chi(s)}{s^2+\xi^2}\dd s.
$$

具体的には上半面の Cauchy 表示から、下半面で解析的な $\chi(z^*)^*$ の実軸積分（上半面の点に関して下側に閉じると零）を引く。二つの境界値の差が $2i\Im\chi$ であるため、虚部だけの核が得られる。受動な振動子なら被積分関数は非負であり、$\xi=0$ の値から高周波に向かって単調に減少する。

Lorentz 模型では直接代入して $\chi(i\xi)=\omega_p^2/(\omega_0^2+\xi^2+\gamma\xi)$ となる。両辺に $\xi^2$ を掛けて $\xi\to\infty$ とすれば、積分中の $\xi^2/(s^2+\xi^2)$ は単調に一へ近づく。よって

$$
\int_0^\infty s\Im\chi(s)\dd s=\frac{\pi}2\omega_p^2.
$$

減衰率 $\gamma>0$ を変えると吸収線の幅と高さが変わるが、この重み付き面積は変わらない。一方、静的感受率を決めるのは $\Im\chi(s)/s$ の積分であり、低周波側に重みがある。共鳴を二つに分ける場合には、それぞれの振動子強度を足した値が上の和則の右辺に現れる。


:::example[有限の測定帯域から静的感受率を再構成する]
緩和時間 $\tau$ を持つ媒質の損失は $\chi_s\omega\tau/[1+(\omega\tau)^2]$ である。静的な分散関係を周波数 $\Omega$ までしか積分しないとき、得られる割合は $2\int_0^{\Omega\tau}(1+u^2)^{-1}\dd u/\pi$ となる。

$x=\Omega\tau$、$Y$ を復元した割合の百分率、$Y_*=100$ とする。$\tau=1.00\,\mathrm{ns}$ なら $x=1$ は $159\,\mathrm{MHz}$ に相当する。

$$
f(x)=\frac{2 \operatorname{atan}{\left(x \right)}}{\pi},\qquad Y=Y_* f(x).
$$

変数を少し動かしたときの読みの変化は、接線の傾きから計算できる。基準値 $x_0$ の近くで

$$
\begin{aligned}
\frac{\Delta Y}{Y_*}= f\!\left(x_0+\Delta x\right)-f(x_0)= f'(x_0)\Delta x+O((\Delta x)^2),\\
f'(x)=\frac{2}{\pi \left(x^{2} + 1\right)}.
\end{aligned}
$$


| $x$ | $f(x)$ | $Y$ $[\%]$ | $\dd Y/\dd x$ $[\%]$ |
|---:|---:|---:|---:|
| $\frac{1}{10}$ | $0.06345$ | $6.345$ | $63.03$ |
| $1$ | $0.5$ | $50$ | $31.83$ |
| $3$ | $0.7952$ | $79.52$ | $6.366$ |
| $10$ | $0.9365$ | $93.65$ | $0.6303$ |

損失ピークまでを覆っても全体の半分しか回収できない。少数の代表周波数だけで静的応答を推定する場合には、未測定の低周波端と高周波端の寄与を評価する必要がある。
:::

## 境界での反射と屈折（Fresnel の式）

### 接線位相と振幅の接続

$z=0$ を境に屈折率 $n_1,n_2>0$、インピーダンス $Z_1,Z_2>0$ の損失のない等方媒質を置く。
入射面は $xz$、角は法線から測る。
界面のすべての $x,t$ で波を足せるには、周波数と接線波数が共通である。
従って反射角は入射角 $\theta_i$ に等しく、

$$
n_1\sin\theta_i=n_2\sin\theta_t
$$

が成立する。

$s$ 偏光は $\hat{\vb{s}}=\hat{\vb{y}}$、$p$ 偏光はそれぞれの波ごとに $\hat{\vb{p}}=\hat{\vb{s}}\times\hat{\vb{k}}$ と定める。
この選択では反射波の $\hat{\vb{p}}$ の $x$ 成分は入射波と逆向きになる。

$r,t$ はこれらの基底に沿う電場の振幅を入射振幅で割った値とする。
自由面電荷・面電流がないので接線 $E,H$ が連続する。

$$
C_i=\cos\theta_i,\ C_t=\cos\theta_t
$$

とおくと、$s$ 波では $1+r_s=t_s$、$C_i(1-r_s)/Z_1=C_t t_s/Z_2$ である。

$p$ 波では $C_i(1-r_p)=C_t t_p$、$(1+r_p)/Z_1=t_p/Z_2$ である。
連立して :term[Fresnel の式]{en="Fresnel equations" yomi="ふれねるのしき"} を得る。

$$
\begin{aligned}
r_s&=\frac{Z_2C_i-Z_1C_t}{Z_2C_i+Z_1C_t},&t_s&=\frac{2Z_2C_i}{Z_2C_i+Z_1C_t},\\
r_p&=\frac{Z_1C_i-Z_2C_t}{Z_1C_i+Z_2C_t},&t_p&=\frac{2Z_2C_i}{Z_1C_i+Z_2C_t}.
\end{aligned}
\label{eq-em-fresnel}
$$

入射波が伝播波で、両側が正の実数 $\varepsilon,\mu$ をもつ半無限媒質であるとする。時間平均した法線方向の電力で定義すると、$R=|r|^2$、$T=[Z_1\operatorname{Re}(C_t)/(Z_2 C_i)]|t|^2$ である。これは $\operatorname{Re}(\vb E\times\vb H^*)/2$ を界面の法線に射影して得られる。全反射では、透過側で減衰する根を選ぶと $C_t$ は純虚数である。従って界面の向こうに場が存在しても $T=0$ となる。以下の実数による検算では透過波も伝播波とする。
例えば $s$ 波で共通分母の二乗から反射の分子の二乗を引くと $4 Z_1Z_2 C_i C_t$ が得られ、$R+T=1$ となる。

$p$ 波でも同じ積が残る。

:::example[空気からガラスへの垂直入射]
非磁性媒質で $n_1=1,\ n_2=1.5$ とすると $Z_2/Z_1=2/3$ である。
$r_s=-0.20,\ r_p=+0.20$ となる。
反射の $p$ 基底の反転を含めれば電場の位相変化は両者とも同じである。
反射率は $0.040$、透過率は $0.960$ となる。
:::

### 吸収媒質の界面に入る電力

透過側で誘電率や透磁率が複素数になるときは、接線電場を両側で同じ座標に沿って測ると扱いやすい。時間因子を $e^{-i\omega t}$ とする。受動な半無限媒質では $k_z^2=\omega^2\mu\varepsilon-k_x^2$ のうち、内部へ進むほど場が減衰する根を選ぶ。損失を零にする極限で二根が実数になる場合は、外向きに電力を運ぶ根を取る。

$s$ 波では $\vb k\times\vb E=\omega\mu\vb H$ から $-H_x/E_y=k_z/(\omega\mu)=Y_s$ を得る。$p$ 波では $\vb k\times\vb H=-\omega\varepsilon\vb E$ から $H_y/E_x=\omega\varepsilon/k_z=Y_p$ を得る。いずれも界面を通る平均電力は $\operatorname{Re}(Y)|E_{\rm tan}|^2/2$ となる。入射側には損失がなく、入射波は伝播するものとする。接線電場の振幅比を $r_{\rm tan},t_{\rm tan}$ とすると、反射で磁場と電場の比の符号が反転するため

$$
1+r_{\rm tan}=t_{\rm tan},\qquad Y_1(1-r_{\rm tan})=Y_2t_{\rm tan}.
$$

従って $r_{\rm tan}=(Y_1-Y_2)/(Y_1+Y_2)$、媒質に入る電力の割合は $(\operatorname{Re}Y_2/Y_1)|t_{\rm tan}|^2=1-|r_{\rm tan}|^2$ である。先の偏光基底では $r_{\rm tan}=r_s$ または $r_{\rm tan}=-r_p$ なので、どちらでも反射率は一致する。例えば真空から複素屈折率 $2+i$ の非磁性媒質へ垂直に入射すれば、$Y_2/Y_1=2+i$ であり、反射率は $|(-1-i)/(3+i)|^2=1/5$ である。残る $4/5$ は界面を通過し、その後に媒質内で吸収される。

### Brewster 角と全反射

非磁性媒質で $Z_i=Z_0/n_i$ を使えば $r_p=0$ の条件は

$$
n_2\cos\theta_i=n_1\cos\theta_t
$$

である。
Snell の式と掛け合わせると異なる屈折率では $\theta_i+\theta_t=\pi/2$、従って $\tan\theta_B=n_2/n_1$ を得る。

空気から $n_2=1.5$ のガラスなら $\theta_B=56.3^\circ$ で、反射光には $s$ 成分だけが残る。

逆向きで $n_1>n_2$、$\sin\theta_i>n_2/n_1$ なら透過側では法線波数が純虚数になる。
界面から遠ざかる向きで $k_{2z}=i\kappa$、$\kappa=(\omega/c)\sqrt{n_1^2\sin^2\theta_i-n_2^2}$ と選ぶ。

すると $C_t=ic\kappa/(n_2\omega)$ なので :eqref[eq-em-fresnel] の反射係数の大きさは両方 $1$ になる。
透過場は $e^{-\kappa z}$ で減衰し、単一界面の損失なしの条件では法線方向の平均電力は零である。

界面に沿う電力は残りうる。

:::example[臨界角を越えた光]
$n_1=1.5$ のガラスから空気への臨界角は $\sin^{-1}(1/1.5)=41.8^\circ$ である。
$\lambda_0=600\,\mathrm{nm}$、$\theta_i=45^\circ$ では振幅の侵入長

$$
1/\kappa=\lambda_0/[2\pi\sqrt{1.5^2/2-1}]=270\,\mathrm{nm}
$$

となる。

近くに別の高屈折率物質を置けば二界面条件を解く必要があり、透過電力が生じうる。
:::

### 膜の二つの界面での多重反射

非磁性の三層に垂直入射し、中央層の厚さを $d$、片道位相を $\delta=n_2\omega d/c$ とする。
同じ実験室の接線軸で電場を測ると、$r_{ij}=(n_i-n_j)/(n_i+n_j)$、$t_{ij}=1+r_{ij}$ である。

第一面で反射した光のほかに、第二面まで一往復して戻る $t_{12}r_{23}t_{21}e^{2i\delta}$ がある。
以後の各往復には $r_{21}r_{23}e^{2i\delta}$ が掛かる。
等比級数を足し、$r_{21}=-r_{12}$ と $t_{12}t_{21}=1-r_{12}^2$ を入れれば

$$
r=\frac{r_{12}+r_{23}e^{2i\delta}}{1+r_{12}r_{23}e^{2i\delta}},\qquad t=\frac{t_{12}t_{23}e^{i\delta}}{1+r_{12}r_{23}e^{2i\delta}}.
\label{eq-em-slab-amplitudes}
$$

層が対称 $n_1=n_3$ なら、どちらの面が一枚では反射しても $\delta=m\pi$ では $r=0$ となる。
板厚を変化させて透過を測ると、媒質内の半波長ごとに極大が並ぶ。

### 全反射した有限幅の光をどこで観測するか

全反射の反射率は $1$ だが、反射振幅の位相は角度によって変わる。その変化を測るには、ただ一つの平面波の代わりに、狭い角度幅をもつビームを界面に当てる。非磁性で損失のない二媒質とし、$z=0$ の界面に沿う波数を $K$ と書く。$k_j=n_j\omega/c$、$q=\sqrt{k_1^2-K^2}>0$、$\kappa=\sqrt{K^2-k_2^2}>0$ とすると、$s$ 偏光の反射振幅は

$$
r_s=\frac{q-i\kappa}{q+i\kappa}=e^{i\phi_s(K)},\qquad
\phi_s=-2\arctan(\kappa/q)
$$

である。実数の $\kappa$ は透過側の振幅が減衰する条件から選ばれている。入射ビームの界面上の電場を $\int A(K)e^{iKx}\dd K$ とし、その中心波数を $K_0$ とする。反射場では被積分関数に $r_s(K)$ が掛かる。$\phi_s$ を $K_0$ のまわりで一次まで展開すれば、元の包絡線の引数 $x$ は $x+\phi'_s(K_0)$ に変わる。従って反射ビームの中心は幾何光学の反射点から

$$
\Delta x=-\frac{\dd\phi_s}{\dd K}\bigg|_{K_0}
=\frac{2K_0}{q\kappa}
$$

だけ界面に沿って進む。最後の等号には、$q'=-K/q$、$\kappa'=K/\kappa$ を $\arctan(\kappa/q)$ の微分に代入した。これが :term[Goos–Hänchen シフト]{en="Goos–Hänchen shift" yomi="ぐーすへんしぇんしふと"} である。臨界角の直上で式が大きくなると、ビームに含まれる角度の一部が全反射の範囲から外れたり、位相の二次以降の変化を無視できなくなったりする。その場合は元の角度積分を評価する。

:::example[同じ侵入長をビームのずれと透過率から測る]
先のガラスと空気の例では $\theta_i=45^\circ$ なので $K=q$ である。十分太い $s$ 偏光ビームのずれは $\Delta x=2/\kappa=540\,\mathrm{nm}$ となる。これは界面に沿う長さであり、界面に垂直に測る侵入長の二倍である。

別の測り方として、空気層を厚さ $d$ にし、その向こうに最初と同じガラスを置く。空気層では電場の振幅が減る解と増す解の両方を許さなければ、二つ目の界面条件を満たせない。$s$ 偏光について $E_y$ と $\partial_z E_y$ は連続である。第二面を位相の基準にして透過振幅を $t$ とすれば、その直前の二量は $(t,iqt)$ である。厚さ $d$ だけ逆向きに進めて第一面に移すと、第一面での条件は

$$
1+r=t[\cosh(\kappa d)-i(q/\kappa)\sinh(\kappa d)],\qquad
iq(1-r)=t[-\kappa\sinh(\kappa d)+iq\cosh(\kappa d)].
$$

二式から $r$ を消し、両側の媒質が同じなので $T=|t|^2$ と置けば

$$
T=\left[1+\frac{(q^2+\kappa^2)^2}{4q^2\kappa^2}\sinh^2(\kappa d)\right]^{-1}.
$$

$d=1/\kappa$、同じ $45^\circ$ では $q/\kappa=3$ で、$T=0.207$ である。空気層を薄くすると $T\to1$、厚くすると $T\to0$ となる。単一界面では法線方向に電力が流れないことと、最後の界面まで届く二つの解が有限の電力を運ぶことは、異なる境界条件の帰結である。
:::

### 複数の平行界面に接線成分を引き継ぐ

薄膜の反射係数を一枚ずつ掛けても、多重反射の位相は再現できない。そこで電場そのものと、磁場の接線成分を同じ界面で引き継ぐ。非磁性、垂直入射とし、$q=Z_0H_y$、$e=E_x$ と書く。屈折率 $n$ の一様な膜中を進む順方向・逆方向の振幅を $A,B$ とすれば、膜の左端で $e=A+B,q=n(A-B)$、右端で $e=Ae^{i\delta}+Be^{-i\delta}$、$q=n(Ae^{i\delta}-Be^{-i\delta})$ である。$\delta=n k_0d$ として右端から左端へ戻す行列は

$$
M(n,\delta)=\begin{pmatrix}\cos\delta&-i\sin\delta/n\\-in\sin\delta&\cos\delta\end{pmatrix},
\qquad \det M=1.
\label{eq-em-film-transfer}
$$

左から順に番号を付けた膜の積は $M_1M_2\cdots M_N$ となる。基板には戻り波がないとし、その屈折率を $n_s$ とすれば右端の状態は $t(1,n_s)^{\mathsf T}$ である。積をこの列に掛けた結果を $(a,b)^{\mathsf T}$ と書く。入射側の屈折率を $n_0$、入射振幅を一とすれば $1+r=ta,n_0(1-r)=tb$ なので、

$$
t=\frac{2n_0}{n_0a+b},\qquad r=\frac{n_0a-b}{n_0a+b}.
$$

透明な膜では $R=|r|^2$ と $T=(n_s/n_0)|t|^2$ の和が一である。各界面でのエネルギー流は $\Re(eq^*)/(2Z_0)$ に等しく、$e,q$ が連続だからである。以下の各例では真空波長 $\lambda_0=600\,\mathrm{nm}$ を共通に用い、位相厚から実際の膜厚も求める。

:::example[単層の無反射膜]
基板の手前に、光学厚を四分の一波長にした膜を置く。二つの反射が相殺する条件を単に位相だけでなく、磁場を含めて確かめる。入射側と出射側は $n_0=1$、$n_s=\frac{9}{4}$ とする。

| 層（入射側から） | 屈折率 | 位相厚 $\delta/\pi$ | 厚さ $[\mathrm{nm}]$ |
|---:|---:|---:|---:|
| 1 | $\frac{3}{2}$ | $\frac{1}{2}$ | $100$ |

膜厚は $d=\lambda_0\delta/(2\pi n)$ から求めた。各層を基板側からたどれば、全体の行列と右向きの出射状態は

$$
\begin{aligned}
M_{\rm all}&=\left(\begin{matrix}0 & - \frac{2 i}{3}\\- \frac{3 i}{2} & 0\end{matrix}\right),\\
(a,b)^{\mathsf T}&=M_{\rm all}(1,n_s)^{\mathsf T}=\left(\begin{matrix}- \frac{3 i}{2}\\- \frac{3 i}{2}\end{matrix}\right).
\end{aligned}
$$

入口で入射波と戻り波の電場は加算、磁場は減算される。したがって

$$
\begin{aligned}
r&=0,\quad t=\frac{2 i}{3},\\
R&=0,\quad T=1.
\end{aligned}
$$

実際に各界面のすぐ両側で連続する場を復元できる。表では入射電場の振幅を $1\,\mathrm{V/m}$ とし、$e$ および $q=Z_0H_y$ を同じ単位で書く。最後の列は入射電力に対する正味の流れの比である。

| 面 | $e$ | $q$ | $\Re(eq^*)/n_0$ |
|---|---:|---:|---:|
| 0 | $1$ | $1$ | $1$ |
| 1 | $\frac{2 i}{3}$ | $\frac{3 i}{2}$ | $1$ |

電場だけの比を二乗する方法では、媒質ごとの波のインピーダンスを取り落とす。ここでは接線磁場も追ったため、すべての面で同じ透過電力を回収できる。この整合には振幅条件も必要である。膜の屈折率が両側の幾何平均となるため、境界で二つの反射の大きさがそろう。
:::

:::example[半波長の膜を通過した位相]
周囲と屈折率の違う平行平板を真空で挟む。界面は一枚ずつなら光を反射するが、設計波長では二つの界面が一組となって透過する。入射側と出射側は $n_0=1$、$n_s=1$ とする。

| 層（入射側から） | 屈折率 | 位相厚 $\delta/\pi$ | 厚さ $[\mathrm{nm}]$ |
|---:|---:|---:|---:|
| 1 | $2$ | $1$ | $150$ |

膜厚は $d=\lambda_0\delta/(2\pi n)$ から求めた。各層を基板側からたどれば、全体の行列と右向きの出射状態は

$$
\begin{aligned}
M_{\rm all}&=\left(\begin{matrix}-1 & 0\\0 & -1\end{matrix}\right),\\
(a,b)^{\mathsf T}&=M_{\rm all}(1,n_s)^{\mathsf T}=\left(\begin{matrix}-1\\-1\end{matrix}\right).
\end{aligned}
$$

入口で入射波と戻り波の電場は加算、磁場は減算される。したがって

$$
\begin{aligned}
r&=0,\quad t=-1,\\
R&=0,\quad T=1.
\end{aligned}
$$

実際に各界面のすぐ両側で連続する場を復元できる。表では入射電場の振幅を $1\,\mathrm{V/m}$ とし、$e$ および $q=Z_0H_y$ を同じ単位で書く。最後の列は入射電力に対する正味の流れの比である。

| 面 | $e$ | $q$ | $\Re(eq^*)/n_0$ |
|---|---:|---:|---:|
| 0 | $1$ | $1$ | $1$ |
| 1 | $-1$ | $-1$ | $1$ |

強度だけなら膜のない場合と等しくても透過振幅の符号が反転する。参照光との干渉によってこの位相を検出できる。
:::

:::example[四分の一波長層の一対]
異なる屈折率の二層を一周期だけ積む。それぞれの膜厚を調整して位相厚を等しくすると、二層の行列積が対角になる。入射側と出射側は $n_0=1$、$n_s=1$ とする。

| 層（入射側から） | 屈折率 | 位相厚 $\delta/\pi$ | 厚さ $[\mathrm{nm}]$ |
|---:|---:|---:|---:|
| 1 | $2$ | $\frac{1}{2}$ | $75$ |
| 2 | $\frac{3}{2}$ | $\frac{1}{2}$ | $100$ |

膜厚は $d=\lambda_0\delta/(2\pi n)$ から求めた。各層を基板側からたどれば、全体の行列と右向きの出射状態は

$$
\begin{aligned}
M_{\rm all}&=\left(\begin{matrix}- \frac{3}{4} & 0\\0 & - \frac{4}{3}\end{matrix}\right),\\
(a,b)^{\mathsf T}&=M_{\rm all}(1,n_s)^{\mathsf T}=\left(\begin{matrix}- \frac{3}{4}\\- \frac{4}{3}\end{matrix}\right).
\end{aligned}
$$

入口で入射波と戻り波の電場は加算、磁場は減算される。したがって

$$
\begin{aligned}
r&=- \frac{7}{25},\quad t=- \frac{24}{25},\\
R&=\frac{49}{625},\quad T=\frac{576}{625}.
\end{aligned}
$$

実際に各界面のすぐ両側で連続する場を復元できる。表では入射電場の振幅を $1\,\mathrm{V/m}$ とし、$e$ および $q=Z_0H_y$ を同じ単位で書く。最後の列は入射電力に対する正味の流れの比である。

| 面 | $e$ | $q$ | $\Re(eq^*)/n_0$ |
|---|---:|---:|---:|
| 0 | $\frac{18}{25}$ | $\frac{32}{25}$ | $\frac{576}{625}$ |
| 1 | $\frac{16 i}{25}$ | $\frac{36 i}{25}$ | $\frac{576}{625}$ |
| 2 | $- \frac{24}{25}$ | $- \frac{24}{25}$ | $\frac{576}{625}$ |

電場と磁場の成分が逆の倍率で変わるため、エネルギーを失わずに入口で見える波のアドミタンスが変わる。
:::

:::example[二周期の誘電体反射鏡]
同じ一対を二回重ねる。各面での反射が位相を保つとき、反射率を一周期の値の二倍にしてよいかを確かめる。入射側と出射側は $n_0=1$、$n_s=1$ とする。

| 層（入射側から） | 屈折率 | 位相厚 $\delta/\pi$ | 厚さ $[\mathrm{nm}]$ |
|---:|---:|---:|---:|
| 1 | $2$ | $\frac{1}{2}$ | $75$ |
| 2 | $\frac{3}{2}$ | $\frac{1}{2}$ | $100$ |
| 3 | $2$ | $\frac{1}{2}$ | $75$ |
| 4 | $\frac{3}{2}$ | $\frac{1}{2}$ | $100$ |

膜厚は $d=\lambda_0\delta/(2\pi n)$ から求めた。各層を基板側からたどれば、全体の行列と右向きの出射状態は

$$
\begin{aligned}
M_{\rm all}&=\left(\begin{matrix}\frac{9}{16} & 0\\0 & \frac{16}{9}\end{matrix}\right),\\
(a,b)^{\mathsf T}&=M_{\rm all}(1,n_s)^{\mathsf T}=\left(\begin{matrix}\frac{9}{16}\\\frac{16}{9}\end{matrix}\right).
\end{aligned}
$$

入口で入射波と戻り波の電場は加算、磁場は減算される。したがって

$$
\begin{aligned}
r&=- \frac{175}{337},\quad t=\frac{288}{337},\\
R&=\frac{30625}{113569},\quad T=\frac{82944}{113569}.
\end{aligned}
$$

実際に各界面のすぐ両側で連続する場を復元できる。表では入射電場の振幅を $1\,\mathrm{V/m}$ とし、$e$ および $q=Z_0H_y$ を同じ単位で書く。最後の列は入射電力に対する正味の流れの比である。

| 面 | $e$ | $q$ | $\Re(eq^*)/n_0$ |
|---|---:|---:|---:|
| 0 | $\frac{162}{337}$ | $\frac{512}{337}$ | $\frac{82944}{113569}$ |
| 1 | $\frac{256 i}{337}$ | $\frac{324 i}{337}$ | $\frac{82944}{113569}$ |
| 2 | $- \frac{216}{337}$ | $- \frac{384}{337}$ | $\frac{82944}{113569}$ |
| 3 | $- \frac{192 i}{337}$ | $- \frac{432 i}{337}$ | $\frac{82944}{113569}$ |
| 4 | $\frac{288}{337}$ | $\frac{288}{337}$ | $\frac{82944}{113569}$ |

対角行列の各成分は周期数に対してべきで変わる。これを入口での振幅比に戻すため、エネルギー反射率は単純には加算されない。
:::

:::example[整合点を外れた単層膜]
整合用と同じ材料を設計厚の半分にする。位相が整合条件を外れたとき、残る反射には位相差も生じる。入射側と出射側は $n_0=1$、$n_s=\frac{9}{4}$ とする。

| 層（入射側から） | 屈折率 | 位相厚 $\delta/\pi$ | 厚さ $[\mathrm{nm}]$ |
|---:|---:|---:|---:|
| 1 | $\frac{3}{2}$ | $\frac{1}{4}$ | $50$ |

膜厚は $d=\lambda_0\delta/(2\pi n)$ から求めた。各層を基板側からたどれば、全体の行列と右向きの出射状態は

$$
\begin{aligned}
M_{\rm all}&=\left(\begin{matrix}\frac{\sqrt{2}}{2} & - \frac{\sqrt{2} i}{3}\\- \frac{3 \sqrt{2} i}{4} & \frac{\sqrt{2}}{2}\end{matrix}\right),\\
(a,b)^{\mathsf T}&=M_{\rm all}(1,n_s)^{\mathsf T}=\left(\begin{matrix}\frac{\sqrt{2} \left(2 - 3 i\right)}{4}\\\frac{\sqrt{2} \left(9 - 6 i\right)}{8}\end{matrix}\right).
\end{aligned}
$$

入口で入射波と戻り波の電場は加算、磁場は減算される。したがって

$$
\begin{aligned}
r&=- \frac{65}{313} - \frac{60 i}{313},\quad t=\frac{8 \sqrt{2} \left(13 + 12 i\right)}{313},\\
R&=\frac{25}{313},\quad T=\frac{288}{313}.
\end{aligned}
$$

実際に各界面のすぐ両側で連続する場を復元できる。表では入射電場の振幅を $1\,\mathrm{V/m}$ とし、$e$ および $q=Z_0H_y$ を同じ単位で書く。最後の列は入射電力に対する正味の流れの比である。

| 面 | $e$ | $q$ | $\Re(eq^*)/n_0$ |
|---|---:|---:|---:|
| 0 | $\frac{248}{313} - \frac{60 i}{313}$ | $\frac{378}{313} + \frac{60 i}{313}$ | $\frac{288}{313}$ |
| 1 | $\frac{8 \sqrt{2} \left(13 + 12 i\right)}{313}$ | $\frac{18 \sqrt{2} \left(13 + 12 i\right)}{313}$ | $\frac{288}{313}$ |

膜を置いた表面を基準に振幅の実部と虚部を比較する。複素反射振幅が生じても、材料そのものに損失があることを意味しない。
:::

:::example[層の順序を反転させる]
高屈折率の基板へ低屈折率側の膜から入射する。入射側媒質と基板が違うと、積層の順序は反射の大きさにも関わる。入射側と出射側は $n_0=1$、$n_s=\frac{9}{4}$ とする。

| 層（入射側から） | 屈折率 | 位相厚 $\delta/\pi$ | 厚さ $[\mathrm{nm}]$ |
|---:|---:|---:|---:|
| 1 | $\frac{3}{2}$ | $\frac{1}{2}$ | $100$ |
| 2 | $2$ | $\frac{1}{2}$ | $75$ |

膜厚は $d=\lambda_0\delta/(2\pi n)$ から求めた。各層を基板側からたどれば、全体の行列と右向きの出射状態は

$$
\begin{aligned}
M_{\rm all}&=\left(\begin{matrix}- \frac{4}{3} & 0\\0 & - \frac{3}{4}\end{matrix}\right),\\
(a,b)^{\mathsf T}&=M_{\rm all}(1,n_s)^{\mathsf T}=\left(\begin{matrix}- \frac{4}{3}\\- \frac{27}{16}\end{matrix}\right).
\end{aligned}
$$

入口で入射波と戻り波の電場は加算、磁場は減算される。したがって

$$
\begin{aligned}
r&=- \frac{17}{145},\quad t=- \frac{96}{145},\\
R&=\frac{289}{21025},\quad T=\frac{20736}{21025}.
\end{aligned}
$$

実際に各界面のすぐ両側で連続する場を復元できる。表では入射電場の振幅を $1\,\mathrm{V/m}$ とし、$e$ および $q=Z_0H_y$ を同じ単位で書く。最後の列は入射電力に対する正味の流れの比である。

| 面 | $e$ | $q$ | $\Re(eq^*)/n_0$ |
|---|---:|---:|---:|
| 0 | $\frac{128}{145}$ | $\frac{162}{145}$ | $\frac{20736}{21025}$ |
| 1 | $\frac{108 i}{145}$ | $\frac{192 i}{145}$ | $\frac{20736}{21025}$ |
| 2 | $- \frac{96}{145}$ | $- \frac{216}{145}$ | $\frac{20736}{21025}$ |

逆側からの入射で相反性を検査するには、膜の順序に加えて入射媒質と基板も交換し、電力流で振幅を規格化する。
:::


:::example[全反射で外側に残る場の深さ]
屈折率 $n_1=1.50$ の媒質から空気 $n_2=1.00$ へ光を入れる。接線波数の連続から $x=n_1\sin\theta_i$ を用いると外部の法線波数は $ik_0\sqrt{x^2-1}$ となる。$x>1$ での強度の侵入長は振幅の侵入長の半分である。

真空波長 $600\,\mathrm{nm}$、$Y$ は強度が $1/e$ まで下がる長さ、$Y_*=\lambda_0/(4\pi)=47.75\,\mathrm{nm}$ とする。

$$
f(x)=\frac{1}{\sqrt{x - 1} \sqrt{x + 1}},\qquad Y=Y_* f(x).
$$

二点の読みの差を求めるには、基準を同じに保ち、振幅の換算を最後に一度だけ行う。最初と最後の設定で

$$
\Delta Y/Y_* = f\!\left(\frac{7}{5}\right)-f\!\left(\frac{101}{100}\right)=-6.033.
$$

第四列は初期設定からの変化である。

| $x$ | $f(x)$ | $Y$ $[\mathrm{nm}]$ | $Y-Y(x_0)$ $[\mathrm{nm}]$ |
|---:|---:|---:|---:|
| $\frac{101}{100}$ | $7.053$ | $336.8$ | $0$ |
| $\frac{11}{10}$ | $2.182$ | $104.2$ | $-232.6$ |
| $\frac{5}{4}$ | $1.333$ | $63.66$ | $-273.1$ |
| $\frac{7}{5}$ | $1.021$ | $48.73$ | $-288$ |

臨界角の直上では外側に広い場ができる。無吸収の半空間へ法線方向の平均電力が流れなくても場は存在し、近くに第二の高屈折率媒質を置くと減衰場を介してエネルギーが渡る。
:::

## 導波管

### 金属壁による横方向の量子化

完全導体の筒の内部を無源で一様な媒質とし、$z,t$ 依存性を $e^{i\beta z-i\omega t}$ とする。
媒質の波速 $v=1/\sqrt{\mu\varepsilon}$、$k=\omega/v$、$k_c^2=k^2-\beta^2$ と置く。

各縦成分 $F=E_z,H_z$ は $(\nabla_t^2+k_c^2)F=0$ を満たす。
添字 $t$ は横断面成分である。
二つの回転方程式の横成分から横場を消去すると、$k_c\ne0$ では

$$
\vb{E}_t=\frac{i}{k_c^2}(\beta\nabla_t E_z-\omega\mu\hat{\vb{z}}\times\nabla_t H_z),\qquad \vb{H}_t=\frac{i}{k_c^2}(\beta\nabla_t H_z+\omega\varepsilon\hat{\vb{z}}\times\nabla_t E_z).
\label{eq-em-guide-fields}
$$

例えば Faraday 則の横成分は

$$
i\omega\mu\vb{H}_t=\hat{\vb{z}}\times(i\beta\vb{E}_t-\nabla_tE_z)
$$

である。
Ampère 則の対応式をこれに代入し、$\hat{\vb{z}}\times(\hat{\vb{z}}\times\vb{E}_t)=-\vb{E}_t$ を使えば分母に $k^2-\beta^2$ が現れる。

:term[TM モード]{en="transverse magnetic mode" yomi="てぃーえむもーど"} では $H_z=0$ であり、壁に接する $E_z$ は零となる。
:term[TE モード]{en="transverse electric mode" yomi="てぃーいーもーど"} では $E_z=0$ であり、壁に沿う横電場を零にするには $\partial_n H_z=0$ が必要である。
長方形 $0<x<a,\ 0<y<b$ の分離解は TM では $\sin(m\pi x/a)\sin(n\pi y/b)$、TE では余弦どうしの積になる。

$$
k_c^2=(m\pi/a)^2+(n\pi/b)^2,\qquad f_c=\frac{v k_c}{2\pi},\quad \beta^2=\omega^2/v^2-k_c^2.
\label{eq-em-guide-cutoff}
$$

TM は $m,n\ge1$、TE は $m,n\ge0$ で両方零を除く。
TE の定数 $H_z$ は時間変動する軸磁束をもつが、壁に沿った電場の周回積分が零であるため許されない。
$f<f_c$ なら伝播せず指数減衰となる。
伝播帯の位相・群速度は

$$
v_p=\omega/\beta,\ v_g=v^2\beta/\omega
$$

で、積は $v^2$ である。

:term[TEM モード]{en="transverse electromagnetic mode" yomi="てぃーいーえむもーど"} では $E_z=H_z=0$ で $k_c=0$ である。
この場合は :eqref[eq-em-guide-fields] で割らず、横断面の $\nabla_t\times\vb{E}_t=0$、$\nabla_t\cdot\vb{E}_t=0$ から調和電位を解く。

境界が連結した中空管では壁が等電位なので :ref[thm-em-uniqueness] により横電場も零である。
同軸線には二つの電極があり、電位差を指定できるため遮断のない TEM 波が存在する。

:::example[長方形導波管]
$a=22.86\,\mathrm{mm}$、$b=10.16\,\mathrm{mm}$ の真空管で最低次は $\mathrm{TE}_{10}$ であり $f_c=c/(2a)=6.56\,\mathrm{GHz}$ である。
$10.0\,\mathrm{GHz}$ では $v_g/c=\sqrt{1-(f_c/f)^2}=0.755$、管内波長は

$$
\lambda_g=(c/f)/\sqrt{1-(f_c/f)^2}=39.7\,\mathrm{mm}
$$

になる。
:::

### 円形導波管と空洞共振器

半径 $a$ の円形断面では正則な固有関数は $J_m(k_cs)e^{im\varphi}$ である。
TM では $J_m(k_ca)=0$、TE では $J'_m(k_ca)=0$ の正根を採る。
最低次は $k_ca\simeq1.841$ の $\mathrm{TE}_{11}$ であり、角度依存の余弦・正弦の二通りが縮退する。

TM の最小根は $J_0$ の最初の零点 $2.405$ である。

長方形管の両端 $z=0,d$ を完全導体で閉じると、進行波と逆行波が重なる。
TE では横電場の係数を $\sin(p\pi z/d)$ とし $p\ge1$、TM では $E_z\propto\cos(p\pi z/d)$ としてその微分すなわち接線電場を端で消し $p\ge0$ となる。

従って

$$
f_{mnp}=\frac v2\sqrt{(m/a)^2+(n/b)^2+(p/d)^2}.
\label{eq-em-cavity}
$$

例えば真空の $a=d=3.00\,\mathrm{cm}$ の空洞の $\mathrm{TE}_{101}$ は $7.07\,\mathrm{GHz}$ である。
壁に微小損失があれば空洞の蓄積エネルギー $U$ と平均損失 $P_{\rm loss}$ から $Q=\omega U/P_{\rm loss}$ を定める。

狭帯域で損失率が一定なら無駆動で $\dot U=-\omega U/Q$ となり、エネルギーの減衰時定数は $Q/\omega$ である。

:::example[横成分から導波管の搬送電力を求める]
矩形管の $\mathrm{TE}_{10}$ で

$$
H_z=H_a\cos(\pi x/a)e^{i\beta z}
$$

とする。
$k_c=\pi/a$ と置き :eqref[eq-em-guide-fields] に代入すると

$$
E_y=i\frac{\omega\mu}{k_c}H_a\sin(k_cx)e^{i\beta z},\qquad H_x=-i\frac{\beta}{k_c}H_a\sin(k_cx)e^{i\beta z}.
$$

残りの横成分は零である。
$z$ 方向の平均流は $-\Re(E_yH_x^*)/2$ なので、断面積分は

$$
P=ab\,\omega\mu\beta|H_a|^2/(4k_c^2)
$$

となる。
単位長さあたりの平均エネルギーは $\int[\varepsilon|E_y|^2+\mu(|H_z|^2+|H_x|^2)]\,\dd A/4$ であり、$k^2=k_c^2+\beta^2$ を使うと $U'=ab\mu k^2|H_a|^2/(4k_c^2)$ である。

従って $P/U'=\omega\beta/k^2=v_g$ となる。
無損失導波管では群速度がエネルギーの輸送速度にも等しい。
:::

### 壁損失を空洞の寿命に結びつける

表皮深さ $\delta_s$ が壁の厚さと曲率半径より十分小さい良導体では、壁内の接線電流は $e^{-s/\delta_s}$ の振幅で減る。
表面抵抗を $R_s=1/(\sigma\delta_s)$ とする。
実際、平面壁内の波動インピーダンスは $Z=\omega\mu/k=(1-i)R_s$ なので、平均電力の流入は $\Re Z\,|H_t|^2/2$ である。

弱い損失なら完全導体として求めた固有場を使い、

$$
P_{\rm wall}=\frac{R_s}{2}\oint|H_t|^2\dd S,\qquad Q_{\rm wall}=\omega\,\frac{\int(\varepsilon|E|^2+\mu|H|^2)\dd V/4}{P_{\rm wall}}
$$

とできる。
誘電体損失や結合孔への流出も小さければ各電力を足し、$1/Q_{\rm total}=\sum_j1/Q_j$ を得る。
例えば壁の $Q$ が $10^4$、外部結合の $Q$ も $10^4$ なら、測定した共振幅は $Q_{\rm total}=5000$ に対応する。

### 遮断のすぐ上で導波管にパルスを送る

横波数 $k_c$ が決まると、透明な一様媒質での軸方向波数は $\beta(\omega)=v^{-1}\sqrt{\omega^2-\omega_c^2}$、$\omega_c=vk_c$ である。
一方向へ長さ $L$ だけ進むスペクトル成分が受ける位相は $\beta L$ である。

周波数微分から

$$
t_g=\frac{\dd(\beta L)}{\dd\omega}=\frac{L\omega}{v\sqrt{\omega^2-\omega_c^2}},\qquad \beta''=-\frac{\omega_c^2}{v(\omega^2-\omega_c^2)^{3/2}} .
$$

遮断の近くでは群遅延が長く、しかも二階微分の絶対値が大きい。
このため入力の周波数帯域にわたって群遅延が変わり、時間波形が広がる。
例えば上記 $\mathrm{TE}_{10}$ の $10.0\,\mathrm{GHz}$ の狭帯域パルスを $1.00\,\mathrm m$ 伝送すると、群遅延は約 $4.42\,\mathrm{ns}$ である。

同じ長さを真空中で進む光の $3.34\,\mathrm{ns}$ より長い。

### 共振器の自由減衰と線幅

一つの空洞モードが弱く減衰し、放出がなくなった後の平均エネルギーが $U(0)\exp(-\omega_0t/Q)$ で減るとする。
振幅はその平方根なので、複素振幅の自由応答は $a(t)=a(0)e^{-i\omega_0t-\kappa t/2}$、$\kappa=\omega_0/Q$ である。

単色源 $se^{-i\omega t}$ が同じモードを駆動するとき、線形な一モード方程式 $\dot a+(i\omega_0+\kappa/2)a=se^{-i\omega t}$ から

$$
|a_\omega|^2=\frac{|s|^2}{(\omega-\omega_0)^2+(\kappa/2)^2}
$$

となる。
半値点の全幅は $\kappa$ である。
例えば $f_0=7.00\,\mathrm{GHz}$、$Q=7000$ なら線幅は $1.00\,\mathrm{MHz}$、エネルギー減衰時間は $Q/(2\pi f_0)=159\,\mathrm{ns}$ である。

振幅の減衰時間はその二倍なので、時間領域でどちらを測定したかを区別する。

### 誘電体スラブの導波条件を境界から導く

厚さ $2a$、屈折率 $n_1$ の平板の両側を $n_2<n_1$ で囲み、すべて非磁性かつ透明とする。$z$ 方向へ $e^{i\beta z}$ で伝わる TE 波の $E_y=f(x)$ は $f''+(n^2 k_0^2-\beta^2)f=0$ を満たす。閉じ込めには $n_2 k_0<\beta<n_1 k_0$ が必要である。$h^2=n_1^2k_0^2-\beta^2$、$q^2=\beta^2-n_2^2k_0^2$ と置く。偶モードを内部 $A\cos(hx)$、右外部 $Be^{-q(x-a)}$ とすれば、接線 $E_y,H_z$ の連続性は $A\cos(ha)=B$ と $hA\sin(ha)=qB$ を与える。奇モードでは内部を $A\sin(hx)$ として同様に解く。従って

$$
u=ha,\quad w=qa,\quad u^2+w^2=V^2=(ak_0)^2(n_1^2-n_2^2),\qquad
w=u\tan u\quad(\text{偶}),\quad w=-u\cot u\quad(\text{奇}).
$$

第一の偶枝は $0<u<\pi/2$ で零から無限大へ単調に増す。$w=\sqrt{V^2-u^2}$ は $0<u<V$ で減少するため、すべての $V>0$ に交点があり基本 TE モードには有限の遮断厚がない。第一の奇枝は $u=\pi/2$ で初めて零から立ち上がるため、$V>\pi/2$ が必要である。金属の矩形導波管における遮断と混同しないためには、外部の減衰場を境界条件に含めることが重要である。TM 波では $H_y$ と $\epsilon^{-1}\partial_xH_y$ が連続するので、偶枝の式は $w=(n_2^2/n_1^2)u\tan u$ に変わる。


次の図では、界面に平行な波数 $\beta$ を固定している。コア内での $h$ と外での $q$ は、この同じ $\beta$ から決まる。

```tikz caption="偶 TE モードの横電場の断面。中央に腹を持ち、誘電体の外側にも減衰しながら続く。"
\begin{tikzpicture}
\draw[->] (-3.2,0)--(3.3,0) node[right] {$x$};
\draw[dashed] (-1.3,-.3)--(-1.3,1.6);
\draw[dashed] (1.3,-.3)--(1.3,1.6);
\draw[thick] (-3,.13) .. controls (-2.2,.21) and (-1.65,.41) .. (-1.3,.65) .. controls (-.45,1.25) and (.45,1.25) .. (1.3,.65) .. controls (1.65,.41) and (2.2,.21) .. (3,.13);
\node[below] at (-1.3,0) {$-a$};\node[below] at (1.3,0) {$a$};
\node[above] at (0,1.25) {$n_1$};\node[above] at (2.4,.5) {$n_2$};
\end{tikzpicture}
```


:::example[矩形導波管の遮断から離れた群速度]
真空で満たした導波管の一つの伝播モードは $\beta^2=(\omega^2-\omega_c^2)/c^2$ を満たす。従って $v_g=\dd\omega/\dd\beta=c^2\beta/\omega$ は平面波の光速より小さい。

$x=\omega_c/\omega$、$Y=v_g$、$Y_*=c=299.792\,\mathrm{Mm/s}$ とする。$1\,\mathrm{Mm}=10^6\,\mathrm m$ で、$x$ は零と一の間である。

$$
f(x)=\sqrt{1 - x} \sqrt{x + 1},\qquad Y=Y_* f(x).
$$

$x=0$ のまわりで先頭の二項を残した近似は

$$
f_{\rm a}(x)=1 - \frac{x^{2}}{2}.
$$

次の表の差は $Y_*[f(x)-f_{\rm a}(x)]$ である。

| $x$ | $f(x)$ | $Y$ $[\mathrm{Mm/s}]$ | 厳密値と近似値の差 $[\mathrm{Mm/s}]$ |
|---:|---:|---:|---:|
| $\frac{1}{20}$ | $0.9987$ | $299.4$ | $-2.345\times10^{-4}$ |
| $\frac{1}{4}$ | $0.9682$ | $290.3$ | $-0.1511$ |
| $\frac{1}{2}$ | $0.866$ | $259.6$ | $-2.691$ |
| $\frac{9}{10}$ | $0.4359$ | $130.7$ | $-47.7$ |

遮断に近いほど小さい $x$ の展開は精度を失う。導波管長を群速度で割って波束の到着を推定するには、波束内の周波数成分がすべて十分に遮断から離れている必要がある。
:::

:::example[共振器から外へ取り出せる蓄積エネルギー]
励振を止めた単一モードが壁で熱になる速さを $\kappa_i U$、結合口から出る速さを $\kappa_e U$ とする。$\dot U=-(\kappa_i+\kappa_e)U$ を零時刻から積分すると、外部へ渡るエネルギーは $U_0\kappa_e/(\kappa_i+\kappa_e)$ となる。

$x=\kappa_e/\kappa_i$、$Y$ を取り出された全エネルギー、$Y_*=U_0=1.00\,\mathrm{mJ}$ とする。

$$
f(x)=\frac{x}{x + 1},\qquad Y=Y_* f(x).
$$

二点の読みの差を求めるには、基準を同じに保ち、振幅の換算を最後に一度だけ行う。最初と最後の設定で

$$
\Delta Y/Y_* = f\!\left(10\right)-f\!\left(\frac{1}{10}\right)=0.8182.
$$

第四列は初期設定からの変化である。

| $x$ | $f(x)$ | $Y$ $[\mathrm{mJ}]$ | $Y-Y(x_0)$ $[\mathrm{mJ}]$ |
|---:|---:|---:|---:|
| $\frac{1}{10}$ | $0.09091$ | $0.09091$ | $0$ |
| $\frac{1}{2}$ | $0.3333$ | $0.3333$ | $0.2424$ |
| $1$ | $0.5$ | $0.5$ | $0.4091$ |
| $3$ | $0.75$ | $0.75$ | $0.6591$ |
| $10$ | $0.9091$ | $0.9091$ | $0.8182$ |

外部結合を強めると取り出せる割合は増すが、同時にエネルギーが減衰する時定数は短くなる。等しい二つの減衰率を選ぶと、自由減衰における全エネルギーの半分ずつが二つの出口へ向かう。
:::

### 板を閉じた空洞の電場と磁場

完全導体の直方体 $0<x<a,0<y<b,0<z<d$ に真空を入れ、$y$ 方向に一様なモードを考える。時間因子を $e^{-i\omega t}$ とすると、

$$
E_y=E_0\sin(\pi x/a)\sin(\pi z/d)
$$

は $x,z$ の四壁で接線電場を消し、$y$ の二壁では電場が法線を向く。内部で発散は零であり、波動方程式から $\omega^2=c^2\pi^2(a^{-2}+d^{-2})$ を得る。Faraday の法則 $\nabla\times\vb E=i\omega\vb B$ を用いれば $B_x=-\partial_z E_y/(i\omega)$、$B_z=\partial_x E_y/(i\omega)$ である。電場の腹では磁場の両成分が零になり、電場と磁場の振幅の位相差は直角である。閉じた空洞の静止した場の模様と一方向の進行波を区別する手掛かりとなる。

積分 $\int_0^a\sin^2(\pi x/a)\dd x=a/2$ などを用いれば、電気エネルギーの周期平均は $\varepsilon_0|E_0|^2abd/16$ である。磁気エネルギーについては二つの微分の二乗が $\pi^2(a^{-2}+d^{-2})/\omega^2=1/c^2$ を与え、同じ値になる。電気双極子プローブを中央へ置けばこのモードに強く結合し、磁気ループを同じ位置に置いても一次では結合しない。壁の近くで磁場を検出する場合は、ループの法線を局所の磁束に合わせる必要がある。

## 演習問題

::::exercise{level=basic}
$100\,\mathrm{MHz}$ の真空波の波長と周期を求めよ。

:::solution
$\lambda=c/f=2.998\,\mathrm m$、$T=1/f=10.0\,\mathrm{ns}$ である。
磁場も同じ空間・時間周期で変化する。
:::
::::

::::exercise{level=basic}
直線偏光が透過軸と $60^\circ$ をなす。
損失のない理想偏光子を通過する強度比を求めよ。

:::solution
透過振幅は元の振幅の $\cos60^\circ=1/2$ であり、強度比はその二乗の $1/4$ となる。
:::
::::

::::exercise{level=basic}
非磁性媒質 $n_1=1$ から $n_2=2$ への垂直入射で $R,T$ を求めよ。

:::solution
接線電場の反射係数は $(1-2)/(1+2)=-1/3$、透過係数は $2/3$ である。
従って $R=1/9$、$T=(n_2/n_1)(2/3)^2=8/9$ となる。
:::
::::

::::exercise{level=basic}
真空で $a=3.00\,\mathrm{cm}$、$b=1.00\,\mathrm{cm}$ の導波管の $\mathrm{TE}_{10}$ と $\mathrm{TE}_{20}$ の遮断周波数を求めよ。

:::solution
$n=0$ では $f_c=mc/(2a)$ である。
$m=1$ で $5.00\,\mathrm{GHz}$、$m=2$ で $9.99\,\mathrm{GHz}$ となる。
:::
::::

::::exercise{level=advanced}
減衰なしの Lorentz 模型で負の誘電率となる正周波数の範囲を求めよ。

:::solution

$$
\varepsilon_r=(\omega_0^2+\omega_p^2-\omega^2)/(\omega_0^2-\omega^2)
$$

である。
分子が正で分母が負となる $\omega_0<\omega<\sqrt{\omega_0^2+\omega_p^2}$ だけで負になる。

非磁性ならこの領域では波数は純虚数である。
:::
::::

::::exercise{level=advanced}
非磁性の単層反射防止膜を垂直入射で設計せよ。
両外側の屈折率を $n_1,n_3$ とする。

:::solution
膜の屈折率を $n_2$ とし、膜内片道位相を $\delta$ とする。
境界条件または多重反射の等比和から $r=(r_{12}+r_{23}e^{2i\delta})/(1+r_{12}r_{23}e^{2i\delta})$、$r_{ij}=(n_i-n_j)/(n_i+n_j)$ である。

$\delta=\pi/2$ なら必要条件は $r_{12}=r_{23}$、すなわち $n_2^2=n_1n_3$ である。
最小厚は $d=\lambda_0/(4n_2)$ となる。
:::
::::

::::exercise{level=advanced}
同軸線の TEM 波の特性インピーダンス $V/I$ を求めよ。
内外半径は $a,b$ とする。

:::solution
横電場 $E_s=V/[s\ln(b/a)]$ に対し $H_\varphi=E_s/Z$、$Z=\sqrt{\mu/\varepsilon}$ である。
内導体電流は $I=2\pi sH_\varphi$ なので $V/I=Z\ln(b/a)/(2\pi)$ となる。

媒質のインピーダンスだけでなく半径比にも依存する。
:::
::::

::::exercise{level=advanced}
円筒空洞（半径 $a$、高さ $d$）の $\mathrm{TM}_{010}$ 共振周波数を求め、高さへの依存を説明せよ。

:::solution
$E_z$ は方位角・高さに一定で、半径方向は $J_0(k_cs)$ である。
側壁で $J_0(k_ca)=0$、第一根は $2.4048$ である。
従って $f=2.4048v/(2\pi a)$ となる。
端板では場は法線方向で、横電場は零なので高さによる位相条件は追加されない。
:::
::::

::::exercise{level=advanced}
単軸結晶で主誘電率が $(\varepsilon_\perp,\varepsilon_\perp,\varepsilon_\parallel)$、透磁率が $\mu_0$ である。
波数を $\vb{k}=(k_x,0,k_z)$ として、二種類の平面波の分散関係を求めよ。

損失・周波数分散は無視する。

:::solution
二つの回転方程式から

$$
\vb{k}\times(\vb{k}\times\vb{E}_0)=-\omega^2\mu_0\boldsymbol{\varepsilon}\vb{E}_0
$$

である。
$E_y$ だけが非零の場合、$k_x^2+k_z^2=\omega^2\mu_0\varepsilon_\perp$ が得られ、位相速度は方向によらない。

$xz$ 面内の場合、二成分は

$$
\begin{pmatrix}
\omega^2\mu_0\varepsilon_\perp-k_z^2&k_x k_z\\
k_x k_z&\omega^2\mu_0\varepsilon_\parallel-k_x^2
\end{pmatrix}\begin{pmatrix}
E_x\\
E_z
\end{pmatrix}=0.
$$

非自明な振幅があるためには行列式が零であり、

$$
k_x^2/\varepsilon_\parallel+k_z^2/\varepsilon_\perp=\omega^2\mu_0
$$

となる。
前者が常光線、後者が異常光線に対応する。
後者の群速度はこの楕円の法線方向を向くので、一般に $\vb{k}$ と平行でない。
:::
::::

::::exercise{level=advanced}
非磁性のガラス $n_1=1.5$ から空気 $n_2=1$ に $\theta_i=45^\circ$ で入射する。
全反射しても偏光が変わる例として、$s,p$ 各反射係数の位相を本章の偏光基底で求めよ。

:::solution
$\sin\theta_t=1.5/\sqrt2$ より $\cos\theta_t=i/\sqrt8$ である。
$s$ 波では $r_s=(n_1\cos\theta_i-n_2\cos\theta_t)/(n_1\cos\theta_i+n_2\cos\theta_t)$ なので、分子・分母をそろえて $(3-i)/(3+i)$ となる。

従って

$$
\arg r_s=-2\arctan(1/3)=-36.87^\circ
$$

である。

$p$ 波では

$$
r_p=(n_2\cos\theta_i-n_1\cos\theta_t)/(n_2\cos\theta_i+n_1\cos\theta_t)=(4-3i)/(4+3i)
$$

であり、

$$
\arg r_p=-2\arctan(3/4)=-73.74^\circ
$$

となる。

両方とも絶対値は $1$ だが位相差は $-36.87^\circ$ である。
両成分を含む直線偏光が反射すると一般には楕円偏光になる。
:::
::::

::::exercise{level=advanced}
真空の一辺 $a$ の立方体空洞について、すべての異なるモードを数えた高周波の累積個数と周波数密度を導け。
壁の表面積に比例する低次の補正は無視する。

:::solution
三つの軸方向の波数間隔は各 $\pi/a$ である。
固有モードを正の八分空間の格子点に対応させると、三つの指数が正の格子点には TE と TM の二つの偏光がある。
一つの指数が零の場合の違いは境界上の格子点で、体積に比例する先頭項には寄与しない。

$k=2\pi f/c$ 以下の個数は

$$
N(f)\simeq 2\frac18\frac{4\pi k^3/3}{(\pi/a)^3}=\frac{8\pi a^3}{3c^3}f^3.
$$

従って単位体積・単位周波数あたりのモード密度は

$$
a^{-3}\dd N/\dd f=8\pi f^2/c^3
$$

である。
この連続近似は波長が空洞より短く、多数の離散準位を含む周波数窓で平滑化した場合に用いる。
:::
::::

::::exercise{level=advanced}
$a=30.0\,\mathrm{mm}$、$b=10.0\,\mathrm{mm}$、$d=40.0\,\mathrm{mm}$ の完全導体壁で囲んだ真空空洞を考える。$E_y=E_0\sin(\pi x/a)\sin(\pi z/d)$、$E_x=E_z=0$ の定在波の周波数と、ピーク電場 $E_0=100\,\mathrm{V/m}$ に対する蓄積エネルギーを求めよ。

:::solution
$f=c\sqrt{a^{-2}+d^{-2}}/2=6.25\,\mathrm{GHz}$ である。電気と磁気の二つの時間平均は等しいので、和は $\varepsilon_0 E_0^2abd/8=1.33\times10^{-13}\,\mathrm J$ となる。理想空洞では漏れも損失もなく、この平均の和は任意の時刻の全エネルギーに等しい。実際、二種類の瞬間エネルギーは $\cos^2\omega t$ と $\sin^2\omega t$ に比例して交互に最大となる。
:::
::::

## 確認テスト

::::quiz{answer="0.04" tol="1e-6"}
空気から屈折率 $1.5$ の非磁性ガラスへの垂直入射の反射率を小数で答えよ。

:::explanation
振幅比の大きさは $(1.5-1)/(1.5+1)=0.2$ なので、反射率は $0.2^2=0.04$ である。
:::
::::

::::quiz
全反射で、損失のない半無限の透過媒質に生じるものを選べ。

- [ ] すべての点で電場が零になる。
- [x] 法線方向に指数減衰する場がある。
- [ ] 法線波数が実数で正になる。

:::explanation
接線波数を保存すると

$$
k_{2z}^2=k_2^2-k_\parallel^2<0
$$

になる。
正の虚部を選ぶと $e^{ik_{2z}z}$ は減衰し、法線平均電力は零である。
:::
::::

::::quiz{answer="15" tol="0.01"}
波長 $600\,\mathrm{nm}$、複屈折の差 $0.010$ の四分の一波長板の最小厚を $\mathrm{\mu m}$ で答えよ。

:::explanation
光路差 $|\Delta n|d$ を $\lambda_0/4$ とするので

$$
d=600\,\mathrm{nm}/0.040=15\,\mathrm{\mu m}
$$

である。
:::
::::

::::quiz
中空の矩形完全導体管で伝播できるものを選べ。
$a>b$、周波数は最低遮断のすぐ上とする。

- [ ] $\mathrm{TEM}$
- [ ] $\mathrm{TM}_{00}$
- [x] $\mathrm{TE}_{10}$

:::explanation
単一導体断面には TEM の電位差が存在しない。
TM は両指数が正であり、TE の許される最小の横波数は $\pi/a$ である。
:::
::::

::::quiz
無損失の非磁性媒質の界面で、Brewster 角で零になる量を選べ。
二つの屈折率は異なるとする。

- [ ] $s$ 偏光の透過振幅
- [x] $p$ 偏光の反射振幅
- [ ] 二つの偏光の反射振幅

:::explanation
Snell の式と接線場の連続性から $\tan\theta_B=n_2/n_1$ のとき

$$
n_2\cos\theta_i=n_1\cos\theta_t
$$

となる。
この組合せが分子に現れる $r_p$ だけが消える。
:::
::::

::::quiz{answer="1" tol="0.001"}
共振周波数 $7.00\,\mathrm{GHz}$、品質係数 $7000$ の空洞の電力半値全幅を $\mathrm{MHz}$ で答えよ。

:::explanation
狭い単一共鳴では $Q=f_0/\Delta f$ である。
従って

$$
\Delta f=7.00\times10^9/7000=1.00\times10^6\,\mathrm{Hz}
$$

となる。
:::
::::

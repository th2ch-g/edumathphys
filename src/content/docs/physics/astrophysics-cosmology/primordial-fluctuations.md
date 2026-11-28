---
title: インフレーションによる原始ゆらぎ
description: 量子場のモードを曲率と重力波のスペクトルへ変換し、ガウス性を検査する量を導く。
notation:
  - symbol: 'v'
    meaning: 正準スカラーモード
  - symbol: 'z_s'
    meaning: 曲率から正準場への変換係数
  - symbol: '\mathcal P_\zeta'
    meaning: 無次元曲率スペクトル
  - symbol: 'r'
    meaning: テンソル・スカラー比
  - symbol: 'f_{\rm NL}'
    meaning: 局所型非ガウス性の係数
---

## この章の目標

- de Sitter 上で短波長真空を設定する。
- Mukhanov–Sasaki 方程式から振幅と傾きを求める。
- テンソルの二偏光と局所型二次場の観測量を計算する。

## de Sitter 時空中の量子場

$c=\hbar=k_B=1$ と計量符号 $(-,+,+,+)$ を用いる。$\eta$ は共形時刻、プライムはその微分、上の点は宇宙時刻の微分とする。Fourier 展開は $e^{i\vb k\cdot\vb x}$ を用い、$k=|\vb k|>0$ とする。平坦な de Sitter パッチを $\eta<0$、$a=-1/(H\eta)$、$H>0$ とする。質量 $m$ の最小結合実場 $\varphi$ の二次作用は $\frac12\int a^2[(\varphi')^2-(\grad\varphi)^2-a^2m^2\varphi^2]\dd\eta\dd^3x$ である。$u=a\varphi$ を代入して交差項を部分積分すると、質量項が $a^2m^2-a''/a$ となる。Fourier モードは $u_k''+[k^2+a^2m^2-a''/a]u_k=0$ である。実場を消滅演算子と生成演算子で展開し、等時刻正準交換関係から $u_k u_k^{*'}-u_k^*u_k'=i$ を得る。平坦時空における正準交換関係は :eqref[eq-qft-scalar-ccr] を参照する。

:term[Bunch–Davies 真空]{en="Bunch–Davies vacuum" yomi="ばんちでいびすしんくう"} は $-k\eta\gg1$ で $u_k\to e^{-ik\eta}/\sqrt{2k}$ とする初期条件である。$m=0$ のモード方程式には

$$
u_k=\frac{e^{-ik\eta}}{\sqrt{2k}}\left(1-\frac{i}{k\eta}\right) \label{eq-cosmo-desitter-mode}
$$

が成り立つ。$a''/a=2/\eta^2$ を入れた微分方程式を直接満たし、短波長で規格化が成立する。$|\eta|\to0$ で $u_k$ は増えるが $\varphi_k=u_k/a$ は有限に凍結する。厳密に質量零の場の赤外領域には初期時刻や空間体積による扱いが必要で、ここでは各非零 $k$ の振幅を求める。

:::example[凍結値]
$|\varphi_k|^2=H^2(1+k^2\eta^2)/(2k^3)$ なので $-k\eta=1$ での値は遅い極限の二倍になる。
:::


:::example[波数の比]
凍結した $|\varphi_k|^2$ は波数二倍で八分の一になるが、$\log k$ 幅当たりには後述の $k^3$ が掛かる。
:::

### 質量がある場の長波長指数

最小結合の自由場で de Sitter 背景 $a=(-H\eta)^{-1}$ とする。正準変数 $u=a\varphi$ に対する有効質量は $a^2m^2-a''/a=(m^2/H^2-2)/\eta^2$ である。$x=-k\eta$、$u=x^{1/2}w(x)$ と置くと $w''+w'/x+(1-\nu^2/x^2)w=0$、$\nu^2=9/4-m^2/H^2$ になる。$0<\nu\le3/2$ の範囲で Bunch–Davies 解は
$$
u=\frac{\sqrt{-\pi\eta}}2 e^{i\pi(2\nu+1)/4}H_\nu^{(1)}(-k\eta).
$$
Hankel 関数の大きな $x$ の漸近 $\sqrt{2/(\pi x)}e^{i(x-\pi\nu/2-\pi/4)}$ で $e^{-ik\eta}/\sqrt{2k}$ に戻る。小さい $x$ で $H_\nu^{(1)}\sim-i\Gamma(\nu)(2/x)^\nu/\pi$ を使うと
$$
\mathcal P_\varphi(k,\eta)
\simeq \left(\frac H{2\pi}\right)^2
 2^{2\nu-3}\left|\frac{\Gamma(\nu)}{\Gamma(3/2)}\right|^2
 (-k\eta)^{3-2\nu}.
$$
この形は実の正の $\nu$ に限る。正の質量では時間に伴って長波長振幅が減る。曲率摂動の保存と、任意の自由場の振幅が一定になることは別である。

:::example[指数が二になる場]
$m^2=2H^2$ では $\nu=1/2$、有効質量項は零になり、$u$ は全時刻で平面波である。従って $\varphi$ のパワーは $(H/(2\pi))^2(k\eta)^2$ と減衰する。
:::

### 正準交換関係と二点関数

$u$ の作用から $\pi_u=u'$ を正準運動量に選ぶ。これは前の部分積分の後で選んだ作用に対応する。場を
$$
\hat u(\eta,\vb x)=\int\frac{\dd^3k}{(2\pi)^3}
 [u_k(\eta)\hat a_{\vb k}e^{i\vb k\cdot\vb x}
  +u_k^*(\eta)\hat a^\dagger_{\vb k}e^{-i\vb k\cdot\vb x}]
$$
と展開する。$[\hat a_{\vb k},\hat a^\dagger_{\vb p}]=(2\pi)^3\delta^3(\vb k-\vb p)$ とおけば $[\hat u(\vb x),\hat\pi_u(\vb y)]$ の積分核は $u_k u_k^{*'}-u_k^*u_k'$ である。これを $i$ に定めると $i\delta^3(\vb x-\vb y)$ になる。モード方程式の周波数二乗は実なので、この Wronskian の時間微分は $u_k u_k^{*''}-u_k^*u_k''=0$ である。従って初期真空での規格化は後まで保たれる。

真空を $\hat a_{\vb k}|0\rangle=0$ とすると、二点関数には消滅演算子と生成演算子の順の積だけが残る。このため等時刻の $\varphi$ の二点は
$$
\langle\hat\varphi(\vb x)\hat\varphi(\vb y)\rangle
=\int\frac{\dd^3k}{(2\pi)^3}\frac{|u_k|^2}{a^2}
e^{i\vb k\cdot(\vb x-\vb y)}
$$
である。同一点で角度を積分すると $k^3|u_k|^2/(2\pi^2a^2)$ が $\log k$ 当たりの分散を与える。厳密な質量零で全ての長波長を含めると $\int_0 \dd k/k$ が発散する。有限体積や初期の赤外尺度を指定すれば、有限な波数帯の比較ができる。逆に同一点の紫外発散は、有限の波数帯を指定する問題とは別に、局所演算子の繰り込みを必要とする。

:::example[有限な波数帯の分散]
全モードが十分に長波長で、$k_{\max}/k_{\min}=e^4$ の帯域を取る。質量零なら帯の分散は $(H/2\pi)^2\log(k_{\max}/k_{\min})=4(H/2\pi)^2$ である。波数の上限と下限をともに同じ倍率で動かしても、理想的な平坦スペクトルならこの値は変わらない。
:::

### 非最小結合と初期状態の選択

作用に $-\xi R\varphi^2/2$ を加えた場合、de Sitter の $R=12H^2$ より上の $m^2$ を $m^2+12\xi H^2$ に置き換える。特に $m=0,\xi=1/6$ では $a^2\xi R=a''/a$ であり、正準場のモード式は全波数で $u_k''+k^2u_k=0$ となる。この場の物理振幅は $a^{-1}$ に減る。質量零でも曲率との結合によって凍結の有無が異なる。

短波長の基準モードを $u_k^{\rm BD}$ とする。別の純粋な Gaussian 状態はモードを $\alpha_k u_k^{\rm BD}+\beta_k u_k^{\rm BD*}$ として与えられる。Wronskian の交差項は消え、複素共役の Wronskian の符号は逆なので、規格化条件は $|\alpha_k|^2-|\beta_k|^2=1$ である。:eqref[eq-cosmo-desitter-mode] の凍結極限では先頭成分が純虚数なので、その場のパワーに掛かる係数は $|\alpha_k-\beta_k|^2$ である。位相差がパワーに残る。任意の $\beta_k$ を全ての波数に与えると励起エネルギー $\int k^3|\beta_k|^2\dd k/a^4$ が大きくなり得るため、背景を固定した計算には紫外の減衰と十分小さい励起が必要である。

:::example[同じ励起数で異なるパワー]
実数 $b>0$ に対して $\alpha=\cosh b$、$\beta=\pm\sinh b$ とする。基準真空からの励起数は両方で $\sinh^2 b$ だが、長波長パワーはそれぞれ $e^{\mp2b}$ 倍である。粒子数だけで二点関数の位相情報は決まらない。
:::

### 複素モードと実確率変数の関係

場は Hermite なので、Fourier 成分について $\hat\varphi_{\vb k}^\dagger=\hat\varphi_{-\vb k}$ である。この二成分を独立な複素乱数として扱うと自由度を二度数える。有限の周期箱で $\vb k$ と $-\vb k$ の組から一方を選び、実部と虚部に分ければそれぞれが等分散の実変数になる。Gaussian 真空の等時刻では場どうしが可換なので、その確率分布で場の等時刻の積を再現できる。後の線形進化で比較に用いるのはこの等時刻共分散であり、異なる時刻の量子交換関係とは区別する。

正準化前の場の二時刻交換関係はモードごとに $\varphi_k(\eta)\varphi_k^*(u)-\varphi_k^*(\eta)\varphi_k(u)$ を係数に持つ。$m=0$ の de Sitter で全体位相を定数だけ回して $\varphi_k=H(1+ik\eta)e^{-ik\eta}/\sqrt{2k^3}$ とする。短い $k\eta$ の展開は
$$
(1+ik\eta)e^{-ik\eta}=1+\tfrac12(k\eta)^2-\tfrac{i}{3}(k\eta)^3+\order{(k\eta)^4}
$$
である。長波長では交換関係を担う虚部が定数成分に比べて小さい。従って後の古典的な線形発展へ、ほぼ同じ実位相を持つランダムな初期振幅を渡せる。時間に依らない確率分布で全ての演算子を同時に記述できることまでは、この近似で仮定していない。

:::example[赤外と紫外で異なる窓の効果]
正則化した二点の距離依存を比べると、角度平均に $\sin(kr)/(kr)$ が掛かる。$k_{\rm min}r\ll1$ で長波長の足し合わせは距離が変わってもほぼ同じである。一方、差 $\langle[\varphi(\vb x)-\varphi(\vb y)]^2\rangle$ では $2[1-\sin(kr)/(kr)]=(kr)^2/3+\order{(kr)^4}$ が現れ、質量零の平坦なスペクトルでも長波長からの寄与は積分可能になる。この計算は紫外側の平滑化を保った比較である。
:::

### 重い場の長波長と境界の対数

$m^2/H^2$ が $9/4$ より大きい場合、有効指数は $\nu=i\lambda$、$\lambda=\sqrt{m^2/H^2-9/4}$ である。長波長のモード方程式では $k^2$ が小さく、$u\propto(-\eta)^\alpha$ を入れると $\alpha(\alpha-1)+m^2/H^2-2=0$ となる。従って $\alpha=1/2\pm i\lambda$ であり、物理場 $u/a$ は $(-\eta)^{3/2}$ の包絡で $\log(-\eta)$ に対し振動する。宇宙時刻で $-\eta\propto e^{-Ht}$ なので、これは減衰する固有時間の振動でもある。$m^2/H^2=9/4$ では二根が重なり、第二解は $(-\eta)^{1/2}\log(-\eta)$ になる。二根の差が零へ近づく際の $\{(-\eta)^{i\lambda}-(-\eta)^{-i\lambda}\}/(2i\lambda)$ の極限からもこの対数を得る。

同じ結果は $k=0$ の固有時方程式 $\ddot\varphi+3H\dot\varphi+m^2\varphi=0$ で見ることができる。$e^{st}$ と置くと $s=-3H/2\pm\sqrt{9H^2/4-m^2}$ である。$m\ll H$ なら遅い根は $-m^2/(3H)+\order{m^4/H^3}$ であり、Hubble 時間に比べてゆっくり減衰する。この二つの表示で時間依存を照合すると、正準場 $u$ の増減と測るべき物理場 $\varphi$ の増減を区別できる。

:::example[重い自由場の振動]
$m^2=13H^2/4$ では $\lambda=1$ である。長波長の場は $e^{-3Ht/2}$ を包絡として、固有時で角周波数 $H$ で振動する。一周期の長さ $2\pi/H$ に対し、包絡の振幅は $e^{-3\pi}$ 倍になる。短波長で設定する真空が同じであっても、長波長の時間発展は質量零と異なる。
:::

## Mukhanov–Sasaki 方程式

単一正準場が平坦な背景を支配し $\dot\phi\ne0$ とする。換算 Planck 質量は $M_{\rm P}^2=(8\pi G)^{-1}$、膨張率は $H=\dot a/a$、$\epsilon_H=-\dot H/H^2$ である。共動ゲージ $\delta\phi=0$ の空間計量を $a^2e^{2\mathcal R}\delta_{ij}$ と置く。正則な断熱長波長解では共動密度が勾配の二次に抑えられ、$\mathcal R$ と $\zeta$ は先頭次数で一致する。時間発展やその近似の一様性は :eqref[eq-cosmo-zeta-conservation] と運動方程式で確認する。lapse を $1+\alpha$、shift の共変空間成分を $\partial_i\beta$ とする。ADM の運動量制約の縦成分は $\partial_i(H\alpha-\dot{\mathcal R})=0$、エネルギー制約は $\nabla^2(\beta+\mathcal R/H)=a^2\epsilon_H\dot{\mathcal R}$ となる。均一な残余は除く。これらを代入した作用は $M_{\rm P}^2\int a^3\epsilon_H[(\dot{\mathcal R})^2-a^{-2}(\grad\mathcal R)^2]\dd t\dd^3x$ である。制約を作用へ代入する過程は下の小節で導く。

$z_s^2=2a^2M_{\rm P}^2\epsilon_H$、$v=z_s\mathcal R$ と置く。$(v'-(z_s'/z_s)v)^2$ の交差項を部分積分し $\frac12\int[(v')^2-(\grad v)^2+(z_s''/z_s)v^2]\dd\eta\dd^3x$ を得る。これを変分した

$$
v_k''+(k^2-z_s''/z_s)v_k=0 \label{eq-cosmo-mukhanov-sasaki}
$$

を :term[Mukhanov–Sasaki 方程式]{en="Mukhanov–Sasaki equation" yomi="むはのふささきほうていしき"} と呼ぶ。

:::example[長波長]
$k=0$ の式を $\mathcal R$ に戻すと $(z_s^2\mathcal R')'=0$ であり、解は定数と $\int \dd\eta/z_s^2$ の和になる。
:::


:::example[減衰モード]
$z_s\propto a\propto(-\eta)^{-1}$ なら積分の時間依存部は $\eta^3$ で、膨張中に消える。
:::

### 制約を消去する二次作用

前の式の導出を補う。空間計量を $h_{ij}=a^2e^{2\mathcal R}\delta_{ij}$、lapse を $N=1+\alpha$ とする。$D_i$ を $h_{ij}$ の共変微分、$N_i=\partial_i\beta$ とし、$E_{ij}=(\dot h_{ij}-D_i N_j-D_j N_i)/2$ と置く。部分積分で時空曲率から境界項を除いた ADM 作用は
$$
S=\frac12\int\sqrt h\left\{M_{\rm P}^2\left[N R^{(3)}+N^{-1}(E_{ij}E^{ij}-E^2)\right]+N^{-1}\dot\phi^2-2NV\right\}\dd^3x\,\dd t .
$$
共動ゲージでは $\phi$ が時刻だけの関数なので、空間勾配は零である。時間微分の運動項が $N^{-1}\dot\phi^2$ として残る。$R^{(3)}=-a^{-2}e^{-2\mathcal R}[4\nabla^2\mathcal R+2(\nabla\mathcal R)^2]$ であり、$E_{ij}E^{ij}-E^2$ の背景を含む項は $-6(H+\dot{\mathcal R})^2+4(H+\dot{\mathcal R})D^2\beta$ となる。残る二次の shift 項は $(\partial_i\partial_j\beta)^2-(\nabla^2\beta)^2$ に比例し、空間積分で相殺する。これは二回の部分積分で両者が $\int\beta\nabla^4\beta$ になるためである。

$N^{-1}=1-\alpha+\alpha^2+\order{\alpha^3}$ と $\sqrt h=a^3e^{3\mathcal R}$ を展開する。背景式と $\dd(a^3H)/\dd t=a^3 H^2(3-\epsilon_H)$ を使い、$\mathcal R\dot{\mathcal R}$ を時間について部分積分すると、$\mathcal R^2$ の項が消える。二次作用を $M_{\rm P}^2/2$ で割った被積分関数は
$$
\begin{aligned}
&a^3[-6\dot{\mathcal R}^2+12H\alpha\dot{\mathcal R}-2(3-\epsilon_H)H^2\alpha^2]\\
&\quad +a[2(\nabla\mathcal R)^2-4\alpha\nabla^2\mathcal R+4(\dot{\mathcal R}-H\alpha)\nabla^2\beta] .
\end{aligned}
$$
$\beta$ の変分で $\alpha=\dot{\mathcal R}/H$ となる。$\alpha$ の変分から $\nabla^2(\beta+\mathcal R/H)=a^2\epsilon_H\dot{\mathcal R}$ が出る。制約を代入すると $\dot{\mathcal R}^2$ の係数は $2a^3\epsilon_H$ である。勾配の交差項 $-4(a/H)\dot{\mathcal R}\nabla^2\mathcal R$ は、空間で部分積分すると $2(a/H)\dd[(\nabla\mathcal R)^2]/\dd t$ となる。時間積分で $\dd(a/H)/\dd t=a(1+\epsilon_H)$ を使うと、元の $2a(\nabla\mathcal R)^2$ と合わせた係数は $-2a\epsilon_H$ となる。従って前の二次作用と :eqref[eq-cosmo-mukhanov-sasaki] が得られる。

### 長波長解が成長する場合

ほぼ一定のポテンシャルで $\ddot\phi+3H\dot\phi\simeq0$、$H$ はほぼ一定とする。$\dot\phi\propto a^{-3}$ から $z_s=a\dot\phi/H\propto a^{-2}$ である。従って長波長積分 $\int\dd\eta/z_s^2=\int\dd t/(az_s^2)$ は、加算定数を除いて $a^3$ に比例する。通常の吸引解で減衰すると見なした第二モードが、ここでは増える。この近似は平坦な領域を通る有限区間に適用し、$\dot\phi=0$ の厳密な de Sitter に共動ゲージをそのまま延ばさない。

:::example[有限区間での増幅]
第二モードが支配し、この背景が $a_f/a_i=2$ まで続けば曲率の振幅は八倍、二点パワーは六十四倍になる。初期の定数モードが優勢なら、この比を全振幅に直接掛けることはできない。
:::

### 平坦面の場と共動面の曲率

空間の平坦ゲージでは $C_{\rm flat}=0$ とし、そこでの場の摂動を $Q$ とする。受動変換で $\delta\phi\mapsto\delta\phi-\bar\phi'T$ であり、$T=Q/\bar\phi'$ を選ぶと共動面へ移る。曲率の定義 $\mathcal R=-C_{\rm com}$ より $\mathcal R=-\mathcal H Q/\bar\phi'=-HQ/\dot{\bar\phi}$ である。一般ゲージで $Q=\delta\phi+(\bar\phi'/\mathcal H)C$ と書くと時間変換で相殺する。$z_s=a\dot{\bar\phi}/H$ を符号付きに選べば $v=z_s\mathcal R=-aQ$ である。モードの全体符号は二点関数に影響しないが、別の場との相関を計算するときは固定する。$Q$ は $\dot{\bar\phi}$ で割らないため、共動座標が退化する時刻で場の摂動を記述する際にも用いられる。

### 曲率モードを滑らかにつなぐ

曲率の作用から正準運動量は $z_s^2\mathcal R'$、運動方程式は $(z_s^2\mathcal R')'=-k^2z_s^2\mathcal R$ となる。背景が有限で、$z_s^2$ が零から離れ、表面の局在源がない転換を考える。曲率に跳びがあると $\int z_s^2(\mathcal R')^2\dd\eta$ が発散するため、有限作用の解では $\mathcal R$ は連続である。さらに運動方程式を転換の前後の狭い区間で積分すると $z_s^2\mathcal R'$ も連続である。背景変化を段差に近似する場合、$v=z_s\mathcal R$ をそのまま両側で連続にすれば良いとは限らない。

:::example[運動量が固定された転換]
有限の $\mathcal R'$ を持つモードが転換を通り、$z_s^2$ がすばやく四倍になるとする。勾配で加わる積分が小さければ、転換直後の $\mathcal R'$ は直前の四分の一である。曲率そのものは連続なので、時間微分を変更した時点で振幅まで四分の一にはならない。
:::

勾配の補正は積分方程式でも追える。$\eta_i$ で $s_i=z_s(\eta_i)^2\mathcal R'(\eta_i)$ と置けば二回積分して
$$
\mathcal R(\eta)=\mathcal R_i+s_i\int_{\eta_i}^{\eta}\frac{\dd u}{z_s(u)^2}
-k^2\int_{\eta_i}^{\eta}\frac{\dd u}{z_s(u)^2}
      \int_{\eta_i}^{u}z_s(w)^2\mathcal R(w)\dd w
$$
となる。初期に $s_i=0$ であっても、右端が有限の二次補正を作る。$z_s$ が急減する区間では外側の積分が大きくなるため、単に $k/(aH)$ の瞬間値を調べるだけでは曲率の一定近似を保証できない。

### 背景変化を含む有効ポテンシャル

$\epsilon_1=-\dd\log H/\dd\log a$、$\epsilon_2=\dd\log\epsilon_1/\dd\log a$、$\epsilon_2\epsilon_3=\dd\epsilon_2/\dd\log a$ とする。$z_s^2=2M_{\rm P}^2a^2\epsilon_1$ と $\mathcal H'=\mathcal H^2(1-\epsilon_1)$ を用いると
$$
\frac{z_s'}{z_s}=\mathcal H(1+\epsilon_2/2),\qquad
\frac{z_s''}{z_s}=\mathcal H^2
\left[2-\epsilon_1+\left(\frac32-\frac12\epsilon_1\right)\epsilon_2
 +\frac{\epsilon_2^2}4+\frac{\epsilon_2\epsilon_3}2\right].
$$
二式目は $(z_s'/z_s)'+(z_s'/z_s)^2$ である。真空初期条件を置くには $k^2\gg|z_s''/z_s|$ が必要である。$\epsilon_2$ が大きく速く変わる場合は、$\mathcal H$ だけで初期値や凍結時刻を判断しない。有限の $\eta_i$ に初期条件を置く実装では、初期値に $\omega_k^2=k^2-z_s''/z_s>0$ と $|\omega_k'|/\omega_k^2\ll1$ を確認し、開始時刻を変えて振幅が収束するか調べる。

:::example[スローロールと同じ有効ポテンシャルを持つ区間]
ほぼ定数 $H$ の区間で $\dot\phi\propto a^{-3}$ なら $\epsilon_2\simeq-6$、$\epsilon_1\simeq0$ である。有効ポテンシャルの括弧は $2-9+9=2$ となるため、通常のスローロールと同様の先頭の $v$ 方程式を持つ。しかし $z_s\propto a^{-2}$ なので $v/z_s$ の成長は全く異なる。モードの形だけで曲率の保存を判定できない。
:::

### 運動項が一般の単一場と音速

場の Lagrangian 密度を $P(X,\phi)$、$X=-g^{\mu\nu}\partial_\mu\phi\partial_\nu\phi/2$ とする。背景で圧力 $p=P$、エネルギー密度 $\rho=2XP_X-P$ であり、$\rho+p=2XP_X=2M_{\rm P}^2H^2\epsilon_H$ である。最後の等式は平坦な場優勢の Friedmann 式から出る。$P_X>0$、$P_X+2XP_{XX}>0$ として勾配項と時間微分項の符号が正の領域に限る。共動面 $\delta\phi=0$ では $X=\dot{\bar\phi}^{\,2}/(2N^2)$ であり、$\delta p=P_X\delta X$、$\delta\rho=(P_X+2XP_{XX})\delta X$ なので、静止枠音速は
$$
c_s^2=\frac{P_X}{P_X+2XP_{XX}}
$$
である。背景の $\bar p'/\bar\rho'$ は一般にこの値と異なる。

共動ゲージの ADM 作用で物質側は $N\sqrt h\,P(\dot{\bar\phi}^{\,2}/(2N^2),\bar\phi)$ となる。$N$ 微分一回は $P-2XP_X=-\rho$、二回は $2X(P_X+2XP_{XX})/N$ である。従って二次作用で変更される lapse の二次係数は、正準の場合の $M_{\rm P}^2a^3 H^2\epsilon_H\alpha^2$ から $M_{\rm P}^2a^3 H^2\epsilon_H\alpha^2/c_s^2$ への移行で得られる。場に空間勾配がなく物質側に shift は現れないので、shift の制約 $\alpha=\dot{\mathcal R}/H$ は同じである。従って最終作用は
$$
S_2=M_{\rm P}^2\int a^3\frac{\epsilon_H}{c_s^2}
\left[\dot{\mathcal R}^{\,2}-\frac{c_s^2}{a^2}(\nabla\mathcal R)^2\right]\dd t\dd^3x .
$$
ここから $z_s^2=2a^2M_{\rm P}^2\epsilon_H/c_s^2$、$v_k''+(c_s^2k^2-z_s''/z_s)v_k=0$ を得る。短波長で定数音速なら $v_k\sim e^{-ic_sk\eta}/\sqrt{2c_sk}$ である。$c_s k$ が一定の有効波数になり、物理的な出入りの境界は $c_sk=aH$ で与えられる。時間で強く変わる音速を同じ定数の解で扱うことはできない。

:::example[簡単な運動項での音速]
$P=X+\lambda X^2-V(\phi)$ なら $P_X=1+2\lambda X$、$P_X+2XP_{XX}=1+6\lambda X$ である。$\lambda X=1/2$ では $c_s^2=2/4=1/2$ である。$X$ の係数が同じ標準場でも、運動エネルギーが有限なら伝播速度が変わる。
:::

## スカラー型スペクトルとスペクトル指数

等方な二点関数を $\langle\zeta_{\vb k}\zeta_{\vb k'}\rangle=(2\pi)^3\delta^3(\vb k+\vb k')P_\zeta(k)$ と定めると、分散の $\dd\log k$ 当たりの寄与は $\mathcal P_\zeta=k^3P_\zeta/(2\pi^2)$ である。ゆっくり変わる $H,\epsilon_H$ のもとで :eqref[eq-cosmo-desitter-mode] を $v$ に用い、$v/z_s$ を長波長で評価する。

$$
\mathcal P_\zeta\simeq\left.\frac{H^2}{8\pi^2M_{\rm P}^2\epsilon_H}\right|_{k=aH},\quad n_s-1\defeq\dv{\log\mathcal P_\zeta}{\log k}\simeq-2\epsilon_H-\epsilon_2 . \label{eq-cosmo-primordial-scalar}
$$

ここで $\epsilon_2=\dd\log\epsilon_H/\dd\log a$ である。$\dd\log k=(1-\epsilon_H)\dd\log a$ と $\dd\log H/\dd\log a=-\epsilon_H$ を使った。:eqref[eq-cosmo-inflaton] とスローロール条件により $\epsilon_2\simeq4\epsilon_V-2\eta_V$、従って $n_s-1\simeq-6\epsilon_V+2\eta_V$ となる。

:::example[弱い赤い傾き]
$\epsilon_H=0.005,\epsilon_2=0.025$ なら $n_s\simeq0.965$ である。
:::


:::example[振幅の縮退]
同じ $H$ で $\epsilon_H$ を二倍にすると $\mathcal P_\zeta$ は半分となる。
:::

### 一定 $\epsilon_H$ のスペクトル指数

$\epsilon_H=\epsilon\in(0,1)$ が一定なら、$H\propto a^{-\epsilon}$、$aH=-[(1-\epsilon)\eta]^{-1}$ である。従って $a,z_s\propto(-\eta)^{-p}$、$p=(1-\epsilon)^{-1}$ であり $z_s''/z_s=p(p+1)/\eta^2$ となる。上の Bessel 方程式に対応する $\nu$ は $p+1/2$ である。長波長で $v/z_s$ は時間に依存しない先頭項を持ち、その波数べきは $\mathcal P_\zeta\propto k^{3-2\nu}$ である。よって
$$
n_s-1=-\frac{2\epsilon}{1-\epsilon}.
$$
前に求めたスローロール式で $\epsilon_2=0$ とした結果は、この厳密な式の一次展開である。

:::example[有限の遅い変化の補正]
$\epsilon=1/20$ なら $n_s-1=-2/19\simeq-0.1053$ である。一次近似の $-0.10$ との差 $-1/190$ は二次以上の効果にあたる。
:::

### 交差時刻とスペクトルの走り

波数を比較するため基準 $k_*$ を選び、$\mathcal P_\zeta(k_*)=A_s$ と書く。$\epsilon_1=\epsilon_H$、$\epsilon_{j+1}=\dd\log|\epsilon_j|/\dd\log a$ を、それぞれ分母が零でない区間で定義する。一次の $n_s-1=-2\epsilon_1-\epsilon_2$ に $\dd/\dd\log k=(1-\epsilon_1)^{-1}\dd/\dd\log a$ を適用すると、先頭の二次の走りは
$$
\alpha_s\defeq\dv{n_s}{\log k}\simeq-2\epsilon_1\epsilon_2-\epsilon_2\epsilon_3 .
$$
$\epsilon_2=0$ などでは積 $\epsilon_2\epsilon_3$ を $\dd\epsilon_2/\dd\log a$ として評価できる。$L=\log(k/k_*)$ と置けば $\log\mathcal P_\zeta=\log A_s+(n_{s,*}-1)L+\alpha_{s,*}L^2/2+\order{L^3}$ である。広い波数域で背景が変わる場合、定数の傾きや短い Taylor 展開の外挿ではなくモード式を解く必要がある。

:::example[二次ポテンシャルの走り]
$N$ を膨張終了までの対数膨張量とすれば、$n_s-1\simeq-4/(2N+1)$ である。$\dd/\dd\log k\simeq-\dd/\dd N$ を使うと $\alpha_s\simeq-8/(2N+1)^2$ となる。$N=50$ では約 $-7.84\times10^{-4}$ である。残りの $N$ は時間とともに減るため、時刻の向きを誤ると走りの符号が逆になる。
:::

### Hankel 解の振幅と交差時刻の規約

一定 $\epsilon_1=\epsilon$ の場合の厳密解を長波長で規格化する。$\nu=(1-\epsilon)^{-1}+1/2$、$x=-k\eta$、$aH=[(1-\epsilon)|\eta|]^{-1}$ とする。Hankel 解から
$$
\mathcal P_\zeta(k,\eta)\simeq
\frac{H(\eta)^2}{8\pi^2 M_{\rm P}^2\epsilon}
\left[\frac{2^{\nu-3/2}\Gamma(\nu)}{\Gamma(3/2)}\right]^2
 (1-\epsilon)^2 x^{3-2\nu}.
$$
$H^2\propto(-\eta)^{2\epsilon/(1-\epsilon)}$ と $x^{3-2\nu}\propto(-\eta)^{-2\epsilon/(1-\epsilon)}$ が相殺する。従って右辺の長波長漸近式の係数を $H_*$ で表すために、交差規約 $k=a_*H_*$ から $x_*=1/(1-\epsilon)$ を代入できる。結果は
$$
\mathcal P_\zeta(k)\simeq
\frac{H_*^2}{8\pi^2M_{\rm P}^2\epsilon}
\left[\frac{2^{\nu-3/2}\Gamma(\nu)}{\Gamma(3/2)}\right]^2
(1-\epsilon)^{2\nu-1}
$$
である。これは凍結振幅を交差時刻の背景量で表した式であり、交差の瞬間の $|v/z_s|^2$ の値とは異なる。

### スローロール振幅の一次補正

小さい $\epsilon_1,\epsilon_2$ が交差付近でほぼ一定なら $\mathcal H\simeq-(1+\epsilon_1)/\eta$ であり、$z_s''/z_s\simeq(2+3\epsilon_1+3\epsilon_2/2)/\eta^2$ となる。$\nu=3/2+\epsilon_1+\epsilon_2/2+\order{\epsilon_i^2}$ を使う。$H^2/\epsilon_1$ の時間変化は $x^{3-2\nu}$ と先頭で相殺し、$x_*=(1-\epsilon_1)^{-1}$ なので $x_*^{3-2\nu}=1+\order{\epsilon_i^2}$ である。Hankel 関数の係数を $\nu=3/2$ の周りで対数微分すると $2[\log2+\Gamma'(3/2)/\Gamma(3/2)]$ であり、$\Gamma'(3/2)/\Gamma(3/2)=2-\gamma_{\rm E}-2\log2$ を用いると $-2C_0$、$C_0=\log2+\gamma_{\rm E}-2$ である。これは $\Gamma(z+1)=z\Gamma(z)$ と倍角公式を微分すれば得られる。さらに $(aH\eta)^{-2}\simeq1-2\epsilon_1$ なので、交差時の先頭振幅に
$$
1-2(C_0+1)\epsilon_1-C_0\epsilon_2
$$
を掛ければよい。これは振幅の相対一次補正である。傾きを一次まで計算するだけなら、この因子の交差時刻での微分は二次なので落とせる。運動が吸引解を外れる有限区間にこの展開を適用するには条件を改めて確かめる。

### 音速の変化とテンソルとの比較

$H,\epsilon_H,c_s$ が音響交差付近でゆっくり変わると仮定する。短波長真空から長波長へ移る近似では :eqref[eq-cosmo-desitter-mode] の $k$ を $c_sk$ に置き換え、$z_s^2\propto a^2\epsilon_H/c_s^2$ で割る。すると三乗の $c_sk$ と $z_s^{-2}$ の $c_s^2$ が組み合わさり
$$
\mathcal P_\zeta\simeq\left.\frac{H^2}{8\pi^2M_{\rm P}^2\epsilon_H c_s}\right|_{c_sk=aH}.
$$
$s=\dd\log c_s/\dd\log a$ とおくと、交差時刻で $\dd\log k=(1-\epsilon_H-s)\dd\log a$ である。従って傾きの先頭は $-2\epsilon_H-\epsilon_2-s$ となる。テンソルの伝播は Einstein 重力では $c=1$ なので、その交差時刻は $k=aH$ である。同じ $k$ の比を取るときは二つの交差時刻で $H$ が異なる。$|\epsilon_H\log c_s|\ll1$ など、間の背景変化を無視できる場合に限り $r\simeq16\epsilon_H c_s$、$r\simeq-8c_s n_t$ と書ける。$c_s$ が小さく交差時刻が遠く離れる場合は $H_t^2/H_s^2$ を戻す必要がある。

### モデルを変えると動く $n_s$ と $r$

正の $\phi$ で $V=A\phi^p$、$A,p>0$ のスローロールを比較する。背景のポテンシャルパラメータは $\epsilon_V=M_{\rm P}^2p^2/(2\phi^2)$、$\eta_V=M_{\rm P}^2p(p-1)/\phi^2$ である。$\epsilon_V=1$ で終了すると近似し、:eqref[eq-cosmo-inflation-efolds] から $\phi_*^2=M_{\rm P}^2(2p N+p^2/2)$ となる。従って
$$
\epsilon_{V,*}=\frac{p}{4N+p},\quad
\eta_{V,*}=\frac{2(p-1)}{4N+p},\quad
n_s-1\simeq-\frac{2(p+2)}{4N+p},\quad
r\simeq\frac{16p}{4N+p}.
$$
$p=2$ なら $n_s-1=-4/(2N+1)$ と $r=16/(2N+1)$ である。$A$ を変えて背景の高さを変えても、$N$ が同じなら傾きと比は先頭で変わらない。逆に同じポテンシャルでも、再加熱の背景を変えて観測する $k$ に対する $N$ が変われば予測点が動く。$\epsilon_V=1$ の終点そのものは一次の近似であるため、精密な比較には背景とモードを積分する。

### 真空から凍結までモードを積分する

一定の $\epsilon_H$ を持つべき膨張を数値例に取る。$p=(1-\epsilon_H)^{-1}$、$x=-k\eta$、$w(x)=\sqrt{2k}\,v_k(\eta)$ なら
$$
w_{xx}+\left(1-\frac{p(p+1)}{x^2}\right)w=0.
$$
真空は大きい $x$ で $w=e^{ix}$、$w_x=iw$ である。時間は $x$ が減る向きに進む。係数が実なので $\operatorname{Im}(w^*w_x)$ は一定で、真空初期値なら $1$ になる。長波長で $w\propto x^{-p}$、$z_s\propto x^{-p}$ なので $x^{2p}|w|^2$ が有限な定数へ近づく。Hankel 解の小引数展開を使えば、その値は $2^{2p}\Gamma(p+1/2)^2/\pi$ である。従って Wronskian の保存と振幅の一致を別々に検査できる。

次のコードは大きい $x$ の有限時刻から小さい $x$ まで進める。$x$ が小さい側では $\omega^2$ が $x^{-2}$ で変化するため、一定の絶対刻みだけでなく相対刻み $|\Delta x|/x$ に上限を設ける。初期真空を平面波に置いた誤差と、積分法の誤差はそれぞれ初期時刻、刻みを変えて検査する。

~~~python
from cmath import exp
from math import gamma, pi

def solve_mode(power, initial_x, final_x=0.025, max_step=0.04):
    x = initial_x
    mode = exp(1j * x)
    derivative = 1j * mode

    def force(position, value):
        return (power * (power + 1) / position ** 2 - 1) * value

    while x > final_x:
        step = -min(max_step, 0.025 * x, x - final_x)
        midpoint = x + step / 2
        first_u, first_v = derivative, force(x, mode)
        second_u = derivative + step * first_v / 2
        second_v = force(midpoint, mode + step * first_u / 2)
        third_u = derivative + step * second_v / 2
        third_v = force(midpoint, mode + step * second_u / 2)
        fourth_u = derivative + step * third_v
        fourth_v = force(x + step, mode + step * third_u)
        mode += step * (first_u + 2 * second_u + 2 * third_u + fourth_u) / 6
        derivative += step * (first_v + 2 * second_v + 2 * third_v + fourth_v) / 6
        x += step

    power_at_end = final_x ** (2 * power) * abs(mode) ** 2
    wronskian = (mode.conjugate() * derivative).imag
    return power_at_end, wronskian

for epsilon in (0.0, 0.05, 0.15):
    power = 1 / (1 - epsilon)
    limit = 2 ** (2 * power) * gamma(power + 0.5) ** 2 / pi
    estimate, conserved = solve_mode(power, initial_x=240)
    print(epsilon, f"{estimate:.6f}", f"{limit:.6f}", f"{conserved:.8f}")
~~~

例えば $\epsilon_H=0,0.05,0.15$ で末端の振幅はそれぞれ $1.000611,1.083186,1.330655$、解析上の凍結値は $1,1.082587,1.330051$ である。Wronskian は全行で $0.99999967$ となる。$\epsilon_H=0$ の行は正準モードの極限を調べる検算である。厳密な de Sitter で $\mathcal R$ を $\dot\phi$ で割って定義する計算ではない。初期値を大きくしても終点 $x$ が有限なら、長波長の振幅に $\order{x^2}$ の補正が残る。例えば $p=1$ の厳密解では $x^2|w|^2=1+x^2$ なので、終点 $x=0.025$ の補正だけで $0.000625$ である。計算の精度を真空の早さと終点の遅さで別々に評価すると、数値誤差を物理的な凍結の時間変化と取り違えずに済む。

### 傾きの二次補正と背景量の階層

スローロール振幅の一次補正を対数に移すと、交差時刻で
$$
\log\mathcal P_\zeta
= 2\log H-\log\epsilon_1-\log(8\pi^2M_{\rm P}^2)
-2(C_0+1)\epsilon_1-C_0\epsilon_2+\order{\epsilon_i^2}
$$
である。階層 $\dd\epsilon_j/\dd\log a=\epsilon_j\epsilon_{j+1}$ に従い高次が一段ずつ小さくなると仮定し、$\dd/\dd\log k=(1-\epsilon_1)^{-1}\dd/\dd\log a$ を使う。先頭の二項の微分は $-2\epsilon_1-\epsilon_2$ であり、交差の Jacobian からさらに $-2\epsilon_1^2-\epsilon_1\epsilon_2$ が出る。振幅の補正を微分した項は $-2(C_0+1)\epsilon_1\epsilon_2-C_0\epsilon_2\epsilon_3$ である。従って
$$
n_s-1=-2\epsilon_1-\epsilon_2-2\epsilon_1^2
-(2C_0+3)\epsilon_1\epsilon_2-C_0\epsilon_2\epsilon_3
+\order{\epsilon_i^3}.
$$
一定 $\epsilon_1$ では $-2\epsilon_1-2\epsilon_1^2$ が残り、厳密な $-2\epsilon_1/(1-\epsilon_1)$ の二次までに一致する。これは交差時の振幅補正だけでなく、波数と時刻の関係にも補正があることを示す。例えば $\epsilon_1=0.05$、$\epsilon_2=0$ では一次が $-0.10$、二次が $-0.105$、厳密値が $-0.10526\ldots$ である。高次の背景量が急に増える場合は、この階層で誤差を見積もれない。

## 原始重力波とテンソル・スカラー比

横・跡なしの $h_{ij}$ の自由作用は $(M_{\rm P}^2/8)\int a^2[(h_{ij}')^2-(\partial_l h_{ij})^2]\dd\eta\dd^3x$ である。:eqref[eq-cosmo-tensor-wave] のテンソル方程式を変分で再現し、係数は Einstein–Hilbert 作用の正規化から決まる。偏光テンソルを $e^\lambda_{ij}e^{\lambda'}_{ij}=\delta_{\lambda\lambda'}$ と規格化し $v_\lambda=aM_{\rm P}h_\lambda/2$ とすれば、二つの正準場が :eqref[eq-cosmo-desitter-mode] と同じ力学を持つ。二偏光の無次元パワーを足し、

$$
\mathcal P_t\simeq 2H_*^2/(\pi^2M_{\rm P}^2),\qquad r=\mathcal P_t/\mathcal P_\zeta\simeq16\epsilon_{H,*},\qquad n_t\simeq-2\epsilon_{H,*} \label{eq-cosmo-primordial-tensor}
$$

を得る。$r$ は :term[テンソル・スカラー比]{en="tensor-to-scalar ratio" yomi="てんそるすからーひ"} であり、$*$ は $k=aH$ での値である。一場・正準・ゆっくりした変化という仮定の下で $r\simeq-8n_t$ が成立する。

:::example[重力波の振幅]
$\epsilon_H=0.001$ なら $r\simeq0.016$ である。
:::


:::example[重力波の傾き]
$r=0.08$ に対応する先頭近似の $n_t$ は $-0.01$ である。
:::

### パワー上限をエネルギー尺度に直す

:eqref[eq-cosmo-primordial-scalar] と :eqref[eq-cosmo-primordial-tensor] の比を定義通り使い、基準波数での $A_s=\mathcal P_\zeta(k_*)$ をとる。重力波のパワーは $rA_s=2H_*^2/(\pi^2M_{\rm P}^2)$ なので $H_*/M_{\rm P}=\pi\sqrt{rA_s/2}$ となる。スローロールで $V\simeq3M_{\rm P}^2H^2$ なら
$$
V_*^{1/4}\simeq M_{\rm P}(3\pi^2 rA_s/2)^{1/4}
$$
である。比較する $r,A_s$ は同じ基準波数で定義しなければならない。別の過程が原始テンソルを作る場合、同じ換算にその寄与を混ぜられない。

:::example[仮の比を入れる]
$A_s=2.1\times10^{-9}$、$r=0.01$ とすれば $H_*/M_{\rm P}\simeq1.02\times10^{-5}$、$V_*^{1/4}/M_{\rm P}\simeq4.20\times10^{-3}$ となる。テンソル比を十六倍にすると、後者は二倍になる。
:::

### 二偏光の作用を ADM 形式から得る

スカラー型を外し、空間計量を $h^{(3)}_{ij}=a^2(e^h)_{ij}$ とする。$h$ は横・跡なしなので $\det e^h=e^{\tr h}=1$ であり、体積要素に一次・二次の $h$ は出ない。テンソルはスカラー・ベクトル制約を一次では源に持たないため、二次作用では $N=1$、$N_i=0$ を代入できる。二次の lapse と shift の項に掛かる係数は背景制約で零である。

$E^i{}_j=H\delta^i_j+(e^{-h}\partial_t e^h)^i{}_j/2$ なので、その跡は厳密に $3H$ である。従って $\tr E^2-(\tr E)^2=-6H^2+\dot h_{ij}\dot h_{ij}/4+\order{h^3}$ となる。空間曲率の積分は、一階の接続を二つ掛けた形へ部分積分できる。直交座標で一次接続は $\Gamma^i{}_{jk}=(\partial_jh_{ik}+\partial_kh_{ij}-\partial_i h_{jk})/2$ である。二次の $\sqrt{\gamma}R^{(3)}$ の空間積分は、境界項を除き $a\int[\delta^{ij}\Gamma^k{}_{il}\Gamma^l{}_{jk}-\delta^{ij}\Gamma^k{}_{ij}\Gamma^l{}_{kl}]\dd^3x$ となる。跡の項は $\Gamma^l{}_{kl}=\partial_k h_{ll}/2=0$ で消える。第一項を展開すると微分添字の交差項は部分積分と $\partial_i h_{ij}=0$ で消え、$-a(\partial_l h_{ij})^2/4$ が残る。Einstein–Hilbert 作用の係数 $M_{\rm P}^2/2$ を掛け、共形時刻へ戻せば前の作用となる。

### 一定の減速率でのテンソル整合関係

$0<\epsilon_H<1$ が一定なら $z_s/a=\sqrt{2\epsilon_H}M_{\rm P}$ も一定であり、$v_k$ と $v_{\lambda,k}$ は同じ有効ポテンシャルと同じ真空初期値を持つ。この場合は近似を使わず $r=16\epsilon_H$ であり、$n_t=n_s-1=-2\epsilon_H/(1-\epsilon_H)$ となる。$\epsilon_H=-n_t/(2-n_t)$ なので $r=-16n_t/(2-n_t)$ である。通常の $-8n_t$ との違いは二次以上である。

### 宇宙を伝わった重力波のエネルギー

原始テンソルの各偏光の振幅が $h_\lambda(k,\eta)=T_h(k,\eta)h_{\lambda,\rm prim}(k)$ と伝わるとする。応力源を無視すれば $T_h$ は :eqref[eq-cosmo-tensor-wave] の自由方程式を満たす。地平線の外で正則な定数モードへの規格化を $T_h\to1$ とする。地平線内で波長より大きく背景変化の尺度より小さい領域にわたりエネルギーを平均する。二次作用の運動項と勾配項が周期平均で同じなので $\rho_{\rm GW}=M_{\rm P}^2\langle h_{ij}'h_{ij}'\rangle/(4a^2)$ である。$\rho_{\rm crit}=3M_{\rm P}^2H^2$ と短波長での $|T_h'|^2\simeq k^2|T_h|^2$ の周期平均を使うと
$$
\frac1{\rho_{\rm crit}}\dv{\rho_{\rm GW}}{\log k}
\simeq\frac1{12}\left(\frac{k}{aH}\right)^2
\mathcal P_t(k)\,\overline{|T_h|^2}.
$$
非零の $h_{ij}$ があっても、地平線の外では局所的な放射エネルギーをこの式で定義できない。

:::example[放射期を伝わる平坦スペクトル]
放射期の自由解 $T_h=\sin(k\eta)/(k\eta)$ を使うと、深い地平線内で $\overline{|T_h|^2}\simeq1/[2(k\eta)^2]$、$aH=1/\eta$ である。従って放射期のこの量は $\mathcal P_t(k)/24$ である。物質期へ進むと波のエネルギーは放射として減り、臨界密度は物質のように減るため、両者の比も低下する。自由流粒子の応力がある場合は、同じ初期パワーでも $T_h$ が変わる。
:::

### 二つの場の軌道が曲がる場合

平坦な場空間を持つ二つの正準スカラーを取り、背景速度の長さを $\dot\sigma=(\dot\phi_1^2+\dot\phi_2^2)^{1/2}>0$、接線の単位ベクトルを $\vb e_\sigma$ とする。法線 $\vb e_s$ を $\dot{\vb e}_\sigma=\dot\vartheta\,\vb e_s$ で規約する。背景の運動式を法線に射影すると $V_s=-\dot\sigma\dot\vartheta$ である。共動面では接線の場摂動は零だが、法線の $Q_s$ は残る。正準な場では $\delta p-\delta\rho=-2\delta V=-2V_s Q_s$ である。共動面の密度は曲率の制約で勾配の二次であり、長波長の :eqref[eq-cosmo-bardeen] に続く共動曲率の発展に入れると
$$
\dot{\mathcal R}\simeq-\frac{2H\dot\vartheta}{\dot\sigma}Q_s
$$
になる。符号は空間計量を $a^2 e^{2\mathcal R}$ とし、$\mathcal R=-C_{\rm com}$ とした規約による。$\dot\vartheta=0$ ならこの線形転換はない。角度の変化があれば、地平線外でも曲率は時間変化し得る。

ある期間の線形転換を $\mathcal R_f=\mathcal R_*+T_{\mathcal R s}\mathcal S_*$、$\mathcal S=HQ_s/\dot\sigma$ と表す。初期に二つのパワーが同じで相関が零なら $P_{\mathcal R,f}=(1+T_{\mathcal R s}^2)P_{\mathcal R,*}$ である。テンソルがその間に追加で生成されない場合、テンソル比は初期値からこの因子で割られる。初期に相関があれば $2T_{\mathcal R s}P_{\mathcal R\mathcal S,*}$ も足す。終時刻の $r$ だけから初期の $\epsilon_H$ を戻すには、一場の保存条件が必要である。

### テンソル波が背景の式を出入りする条件

膨張中の Hubble 半径の出入りを、$k=aH$ という瞬間の規約で指定する。長さの比較で使う物理波長は $2\pi a/k$ だが、この交差はモード方程式の $k^2$ と $\mathcal H^2$ の尺度に合わせるので $2\pi$ を入れない。例えば一定 $w>-1/3$ の後期背景なら $aH\propto a^{-(1+3w)/2}$ なので、高波数のテンソルが先に入る。$k$ が入る時点の尺度因子は $a_k\propto k^{-2/(1+3w)}$ である。地平線の外では一定、深く中へ入ってから $h\propto a^{-1}$ と近似すれば、同じ観測時刻で $|T_h|^2\propto a_k^2$、$\Omega_{\rm GW}(k)\propto k^2 a_k^2\mathcal P_t(k)$ となる。このため背景が一定 $w$ に支配された区間で入ったモードの傾きは
$$
n_t+2-\frac4{1+3w}=n_t+\frac{2(3w-1)}{1+3w}
$$
である。放射では $n_t$、物質では $n_t-2$ である。転換に近い波数や、自由流による応力を受ける波数ではこの単純なべき則からずれる。原始の傾きと観測時刻の傾きを同じものとして扱うには、入った後の背景の伝播を考慮する必要がある。

## 非 Gauss 性の概観

三点関数の Fourier 係数を $\langle\zeta_1\zeta_2\zeta_3\rangle=(2\pi)^3\delta^3(\sum_i\vb k_i)B_\zeta(k_1,k_2,k_3)$ とする。平均零の Gaussian な場では奇数モーメントは符号反転で零となる。実空間で $\zeta=g+(3/5)f_{\rm NL}(g^2-\langle g^2\rangle)$ とする局所型では、二次項を一脚から選んだ四つの $g$ の Wick 収縮のうち、平均値の引き算で内部収縮は消える。残る二通りの外部収縮から、$f_{\rm NL}$ の一次で $B_\zeta=(6/5)f_{\rm NL}[P_g(k_1)P_g(k_2)+P_g(k_2)P_g(k_3)+P_g(k_3)P_g(k_1)]$ となる。三角形の形依存で絞られた局所型、正三角形付近を重く見る型などを区別できる。単一場の長波長が局所の座標の拡大縮小だけを与える場合、短波長パワーの $\log k$ 微分が絞られた三点を決める。高次作用を使う精密な相関関数は :cite[mukhanov-cosmology] の摂動論を基礎とする。

:::example[正三角形]
同じ $P_g=P$ の三脚では $B_\zeta/P^2=18 f_{\rm NL}/5$ である。
:::


:::example[絞られた三角形]
長波長のパワーが短波長より大きいと、長–短の積を含む二項が主になる。
:::

### 長波長と短波長の相関

単一の断熱的な時計があり、十分長い $\zeta_L$ が短波長の局所力学を変えず、空間の拡大縮小だけを与える場合を考える。局所座標 $\tilde{\vb x}=e^{\zeta_L}\vb x$ で二点関数が同じなら、座標 $\vb x$ での短波長場は $g(e^{\zeta_L}\vb x)$ になる。二つの Fourier 変換には $e^{-6\zeta_L}$ が掛かり、運動量デルタ関数の変換で $e^{3\zeta_L}$ が戻る。従って
$$
P_S(k\mid\zeta_L)=e^{-3\zeta_L}P_S(k e^{-\zeta_L}).
$$
線形化すると $\delta P_S=-(n_s-1)\zeta_LP_S$ である。これと長波長の振幅を相関させると、絞られた極限の三点関数は $(1-n_s)P_LP_S$ になる。局所型の定義では支配する二項が $(12/5)f_{\rm NL}P_LP_S$ なので、対応値は $5(1-n_s)/12$ である。この導出は、長波長が局所の物理的な初期状態を変えないという仮定を使う。

:::example[小さな赤い傾き]
$n_s=0.96$ の場合、この極限の局所型に換算した値は $f_{\rm NL}=1/60$ である。短波長二本の波数に比べ長波長が十分小さくない三角形では、この極限式から三点関数全体を推定できない。
:::

### 局所膨張の差から曲率の非線形項を作る

波長が Hubble 半径より十分長く、勾配の寄与を捨てられる各領域を局所的な一様宇宙として扱う。初期の空間平坦面から、終点の一様密度面までの積分膨張を $N(\phi_i)$ とする。終点までの局所の尺度因子の比は $e^N$ なので、背景と比べた終点の空間計量の内在的な倍率は $e^{2(N-\bar N)}$ である。従ってこの近似で $\zeta=\delta N$ である。初期の場が一つの Gaussian 摂動 $\delta\phi$ に支配されるとし、速度などの初期情報が同じ場の値から決まる吸引解を仮定する。Taylor 展開で
$$
\zeta=N_{,\phi}\delta\phi+\frac12N_{,\phi\phi}[(\delta\phi)^2-\langle(\delta\phi)^2\rangle]+\cdots .
$$
$g=N_{,\phi}\delta\phi$ とすると、局所型の定義と比べて $f_{\rm NL}^{(\delta N)}=(5/6)N_{,\phi\phi}/N_{,\phi}^2$ である。これは初期場自身の三点関数を無視した分である。実際の単一場の整合関係を得るには、初期場での相互作用や波数が地平線を出る時刻の寄与も合わせる。

:::example[二次ポテンシャルでの局所変換分]
スローロールの $N=(\phi^2-\phi_{\rm end}^2)/(4M_{\rm P}^2)$ なら $N_{,\phi}=\phi/(2M_{\rm P}^2)$、$N_{,\phi\phi}=1/(2M_{\rm P}^2)$ である。従って局所変換分は $(5/3)M_{\rm P}^2/\phi^2$ となり、$\phi=10M_{\rm P}$ では $1/60$ である。これはこの展開で得た一寄与であり、絞られた三点関数の全振幅を置き換えない。
:::

### 局所型の三点と四点

複数の正準軽場の初期揺らぎ $Q_i$ が独立で同じ Gaussian パワー $P_Q$ を持つとする。吸引解のもとで $\zeta=N_i Q_i+\tfrac12N_{ij}(Q_iQ_j-\langle Q_iQ_j\rangle)$ と展開する。二点の先頭は $P_\zeta=N_l N_l P_Q$ である。三点では二次項の二つの $Q_i$ を残りの線形な脚へ割り当てる二通りがあり、$1/2$ が消える。従って $(6/5)f_{\rm NL}=N_iN_jN_{ij}/(N_lN_l)^2$ となる。分子は初期揺らぎが変換の二次の方向とどの程度そろうかを表す。

四点の連結部では、二次頂点二つの間に内部の Wick 収縮がある。各頂点の外側の相手が一つずつ固定されると、係数は $(N_{ij}N_j)(N_{ik}N_k)$ となる。これを外側三本の $P_\zeta$ の規格化で割った量が $\tau_{\rm NL}=(N_{ij}N_j)^2/(N_lN_l)^3$ である。$a_i=N_i$、$b_i=N_{ij}N_j$ と置けば $(a_ib_i)^2\le (a_i a_i)(b_jb_j)$ なので $\tau_{\rm NL}\ge(6f_{\rm NL}/5)^2$ を得る。等号は $b$ が $a$ に平行な場合である。ここでは初期場の非 Gaussian 性と高次の寄与を無視した。

:::example[三点が零でも残る四点]
$N_1=1,N_2=0,N_{11}=0,N_{12}=N_{21}=1$ とすると $b=(0,1)$ であり、$f_{\rm NL}=0$、$\tau_{\rm NL}=1$ である。三点の振幅だけでは、線形の曲率と直交する方向から来る非線形変換を除けない。
:::

### 相互作用が作る三点の時間積分

初期を Gaussian 真空とし、終時刻での積 $Q$ を考える。相互作用表示の時間発展 $U=1-i\int H_I\dd\eta+\cdots$ を $\langle U^\dagger Q U\rangle$ に代入すると、一次の変化は
$$
-i\int\langle[Q,H_I(\eta)]\rangle_0\dd\eta
=2\,\operatorname{Im}\int\langle Q H_I(\eta)\rangle_0\dd\eta
$$
である。右の等式では、実場の Fourier 三点に対応する係数について複素共役を用いる。これは確率分布へ局所な二次変換を施す計算と違い、相互作用が起きた時刻を積分する。

具体例として de Sitter 上の質量零の試験場と、摂動的な相互作用 Hamiltonian $H_I=\lambda\int a(\varphi')^3\dd^3x/3!$ を仮定する。$\lambda$ は定数とする。この相互作用を背景を無視できる小さい結合の有効模型として扱う。外側の三脚は終時刻 $\eta\to0^-$ で $\varphi_k(0)=H/\sqrt{2k^3}$ である。相互作用の各脚は $\varphi_k^{*'}(\eta)=H k^2\eta e^{ik\eta}/\sqrt{2k^3}$ である。Wick 収縮の $3!$ が分母を打ち消す。運動量のデルタ関数を外した integrand は $K=k_1+k_2+k_3$ を用いて
$$
-\frac{\lambda H^5}{8k_1k_2k_3}\eta^2e^{iK\eta}.
$$
初期に $e^{\varepsilon\eta}$、$\varepsilon>0$ を掛けた積分は $2/(\varepsilon+iK)^3\to 2i/K^3$ なので
$$
B_\varphi(k_1,k_2,k_3)=-\frac{\lambda H^5}{2k_1k_2k_3 K^3}.
$$
毎脚に時間微分があるため、$k_L\ll k_S$ で長波長に $\lambda H^5/(k_L k_S^5)$ の尺度しか現れない。局所型の $P_\varphi(k_L) P_\varphi(k_S)\propto1/(k_L^3 k_S^3)$ で割った比には二次の $(k_L/k_S)^2$ が付く。この例は場の三点を作る独立の相互作用の計算であり、正準単一場の曲率の全三点には重力制約などに由来する別の頂点も入る。

### 原始パワーに関する仮定と形の比較

二つの波数で同じ局所型の振幅を比較するためには、まず $B$ に使う $P_\zeta$ の規格化を定める。例えば尺度不変な $P_\zeta=A/k^3$ で $k_S=k_2=k_3$、$k_L=k_1$ とすれば局所型の $B=(6/5)f_{\rm NL}A^2(2/k_L^3k_S^3+1/k_S^6)$ である。$k_L/k_S$ が $1/10$ なら最後の項は最初の二項の和の $1/2000$ にすぎない。$k_L=k_S$ に近づくと全項を数える必要がある。$f_{\rm NL}$ を異なる形のテンプレートで推定した数値は、同じ三角形の同じ振幅を規格化していなければ直接比較できない。観測の制限には CMB と大規模構造への転送も入るため、原始の二点・三点だけから観測空間の振幅を等しいとみなさない。

### 崩壊する第二の場が作る局所型

インフレーション中のエネルギーは主に別の場に担わせ、軽い試験場 $\sigma$ の初期揺らぎを Gaussian とする。後に二次ポテンシャルで速く振動する区間では、固有時間で一周期を平均すると運動エネルギーとポテンシャルエネルギーが等しく、$\rho_\sigma$ は圧力なしの成分として減る。局所振幅を $\bar\sigma+\delta\sigma$ とし、共通の初期平坦面で $\rho_\sigma/\bar\rho_\sigma=(1+g)^2$、$g=\delta\sigma/\bar\sigma$ とする。物質の局所保存から曲率は $\zeta_\sigma=(1/3)\log(\rho_\sigma/\bar\rho_\sigma)=(2/3)g-(1/3)g^2+\cdots$ である。放射の曲率は零とし、長波長で崩壊が全密度一定の面に一斉に起こる模型を取る。

崩壊直前の背景分率を $\Omega=\bar\rho_\sigma/(\bar\rho_\sigma+\bar\rho_r)$ とする。局所が背景より $\zeta$ だけ多く膨張した崩壊面では、放射と場の密度比はそれぞれ $e^{-4\zeta}$、$e^{3(\zeta_\sigma-\zeta)}$ である。全密度を一定にする条件は
$$
(1-\Omega)e^{-4\zeta}+\Omega e^{3(\zeta_\sigma-\zeta)}=1.
$$
$\zeta=r_D\zeta_\sigma+b\zeta_\sigma^2+\cdots$ とすると、一次で $r_D=3\Omega/(4-\Omega)$ を得る。二次では $(4-\Omega)b=[16(1-\Omega)r_D^2+9\Omega(1-r_D)^2]/2$ だから $b=r_D(1-r_D)(3+r_D)/2$ である。これを $g$ の式へ戻す。線形曲率は $g_\zeta=(2r_D/3)g$、二次の係数は $-r_D/3+4b/9$ となる。平均値を差し引き、この係数を $(2r_D/3)^2$ で割って $(3/5)f_{\rm NL}$ に等置すれば
$$
f_{\rm NL}=\frac{5}{4r_D}-\frac53-\frac{5r_D}6.
$$
この値は突然の崩壊、二次ポテンシャル、放射の原始曲率を無視できるという仮定に依存する。$r_D\to0$ で係数が大きくなるのは、同じ $g$ で二次変換が線形曲率の二乗ほどには速く減らないためである。摂動の小ささは $g$ や崩壊時の密度揺らぎで別に確かめる。

:::example[崩壊分率が変換を変える]
$r_D=1$ で $f_{\rm NL}=-5/4$ となる。これは場が全密度を担えば $\zeta=\zeta_\sigma$ となり、対数の二次係数だけが残る結果と一致する。$r_D=1/10$ なら $f_{\rm NL}=10.75$ であり、同じ線形パワーであっても三点が変わり得る。
:::

### 局所型が二点へ与える高次の補正

小さい係数 $\lambda=3f_{\rm NL}/5$ による局所変換 $\zeta=g+\lambda(g^2-\langle g^2\rangle)$ を有限の波数帯で考える。二点関数の $\lambda$ 一次の寄与は Gaussian 場の三次なので零である。非零の $\vb k$ について二次の項同士を Wick 展開すると、一つの頂点の二本を内部でつなぐ収縮は零波数に限られ、残りの二通りが同じ積分を与える。従ってこの模型で
$$
P_\zeta(k)=P_g(k)+2\lambda^2\int\frac{\dd^3p}{(2\pi)^3}
P_g(p)P_g(|\vb k-\vb p|)
$$
となる。初期場に三次や四次の独立な統計を仮定する場合は、この式の外にそれらの寄与も現れる。

平坦な無次元パワーを想定すると、$p\ll k$ では被積分項は $P_g(k)P_g(p)$ に近い。角度積分後の小波数の分散を $\sigma_L^2=\int_{k_{\min}}^{k_L}\mathcal P_g(p)\dd\log p$ とすれば、$|\vb k-\vb p|\ll k$ の別の領域からも同じ寄与がある。従ってこの二領域による補正は $4\lambda^2 P_g(k)\sigma_L^2$ である。帯域を指定せず無限の波長まで平坦とすれば $\sigma_L^2$ は対数で発散する。有限体積で長波長の実現を比較するという観測の条件と、普遍的な熱平均を取るという理論の条件をそろえる必要がある。$f_{\rm NL}$ の大きさが小さくても、$\lambda^2\sigma_L^2$ を別に確かめることが展開の条件である。

### スカラーとテンソルを基準波数で比べる

固定した基準 $k_0$ で $\mathcal P_\zeta=A_s(k/k_0)^{n_s-1}$、$\mathcal P_t=A_t(k/k_0)^{n_t}$ とする。別の波数での比は
$$
r(k)=r(k_0)(k/k_0)^{n_t-n_s+1}.
$$
例えば $n_s=0.965$、$n_t=-0.01$、$r(k_0)=0.08$ なら、$k=10 k_0$ で比は $0.08\times10^{0.025}\simeq0.0847$ になる。$r$ の値を比較する際は、基準波数だけでなくテンソルの傾きを観測から自由に求めたか、整合関係を仮定したかも区別する。後者を仮定すれば $r$ を動かすと $n_t$ も変わるため、同じ観測データでも別の模型を検定していることになる。

## 演習問題

::::exercise{#exr-cosmo-primordial-fluctuations-1 level=basic}
de Sitter の $a''/a$ を求めよ。
:::solution
$a=-1/(H\eta)$ を二回微分して割ると $2/\eta^2$ である。
:::
::::

::::exercise{#exr-cosmo-primordial-fluctuations-2 level=advanced}
式 :eqref[eq-cosmo-desitter-mode] をモード方程式へ代入せよ。
:::solution
$u_k=e^{-ik\eta}b/\sqrt{2k}$ と置くと左辺は同じ指数因子に $b''-2ikb'-2b/\eta^2$ を掛けたものになる。$b=1-i/(k\eta)$ ではこれが零である。
:::
::::

::::exercise{#exr-cosmo-primordial-fluctuations-3 level=basic}
波数二倍で一定 $\mathcal P_\zeta$ の次元付き $P_\zeta$ は何倍か。
:::solution
$P_\zeta$ は $k^{-3}$ に比例し八分の一である。
:::
::::

::::exercise{#exr-cosmo-primordial-fluctuations-4 level=advanced}
$k=0$ で $v=z_s$ が解になる理由を示し、独立解を与えよ。
:::solution
代入で $z_s''-(z_s''/z_s)z_s=0$。$v=z_s u$ に変換すると $(z_s^2u')'=0$ なので $v=z_s\int\dd\eta/z_s^2$ がもう一つの解である。
:::
::::

::::exercise{#exr-cosmo-primordial-fluctuations-5 level=basic}
$\epsilon_H=0.01,\epsilon_2=0.03$ のスカラー指数を求めよ。
:::solution
$1-2\epsilon_H-\epsilon_2=0.95$ である。
:::
::::

::::exercise{#exr-cosmo-primordial-fluctuations-6 level=advanced}
二次ポテンシャルで $N=50$ の $n_s$ を先頭で求めよ。
:::solution
$\epsilon_V=\eta_V=1/(2N+1)$ なので $n_s=1-4/101\simeq0.9604$ である。
:::
::::

::::exercise{#exr-cosmo-primordial-fluctuations-7 level=basic}
$r=0.064$ の模型の $\epsilon_H$ を求めよ。
:::solution
一場の式で $r/16=0.004$ である。
:::
::::

::::exercise{#exr-cosmo-primordial-fluctuations-8 level=advanced}
局所型変換の一地点三次モーメントを先頭次数で求めよ。
:::solution
$g$ の分散を $\sigma^2$ とすると $3(3/5)f_{\rm NL}(\langle g^4\rangle-\sigma^4)=(18/5)f_{\rm NL}\sigma^4$ である。Gaussian の四次モーメントは三つの Wick 対で $3\sigma^4$ となる。
:::
::::

::::exercise{#exr-cosmo-primordial-fluctuations-9 level=advanced}
質量零で $\xi=1/6$ の自由場を de Sitter に置き、一定の共形時刻で無次元パワーの $k$ 依存性を求めよ。
:::solution
$a^2\xi R=a''/a$ なので正準モードは $u_k=e^{-ik\eta}/\sqrt{2k}$ である。$\mathcal P_\varphi=k^3|u_k|^2/(2\pi^2a^2)=k^2/(4\pi^2a^2)$ となる。一定時刻で $k^2$ に比例し、一定波数では $a^{-2}$ に減衰する。
:::
::::

::::exercise{#exr-cosmo-primordial-fluctuations-10 level=advanced}
$V=A\phi^4$ の正準模型で $N=60$ の $n_s$ と $r$ を先頭近似で求めよ。同じ $N$ の二次ポテンシャルと比べよ。
:::solution
$p=4$ で $\epsilon_V=4/244=1/61$、$\eta_V=6/244=3/122$ である。従って $n_s=1-12/244=58/61\simeq0.951$、$r=64/244=16/61\simeq0.262$ となる。$p=2$ では $n_s=1-8/242=117/121\simeq0.967$、$r=32/242=16/121\simeq0.132$ である。両者の差はポテンシャルの形から出る。ここでは終了を $\epsilon_V=1$ と近似し、二つの背景で $N$ を同じ値に固定した。
:::
::::

::::exercise{#exr-cosmo-primordial-fluctuations-11 level=advanced}
一定 $\epsilon_H$ の一場正準模型で $n_t=-1/19$ とする。厳密な背景のべき則から $r$ を求め、一次の整合関係が与える比と比較せよ。
:::solution
$n_t=-2\epsilon_H/(1-\epsilon_H)$ なので $\epsilon_H=1/39$ である。スカラーとテンソルの正準モードは先頭だけでなく同じ方程式に従うので $r=16/39\simeq0.4103$ である。一次の $-8n_t$ は $8/19\simeq0.4211$ で、差は $8/741$ である。$n_t$ の値を同じに固定した比較のため、$\epsilon_H$ と $n_t$ の間の高次補正が差を作る。
:::
::::

::::exercise{#exr-cosmo-primordial-fluctuations-12 level=advanced}
二つの軽い正準場から曲率とエントロピーが同じパワーで生まれ、初期相関は零とする。地平線外の転換が $\mathcal R_f=\mathcal R_*+2\mathcal S_*$ で、テンソルの振幅は変わらない。終時刻の曲率とテンソル比を初期値に対する倍率で表せ。また $\mathcal R_f$ と $\mathcal S_*$ の相関係数を求めよ。
:::solution
共通の初期パワーを $P$ とすれば、無相関なので $P_{\mathcal R,f}=P+4P=5P$ である。テンソルの初期と終時刻のパワーが同じなので $r_f=r_*/5$ となる。交差パワーは $P_{\mathcal R_f\mathcal S_*}=2P$、二つの自己パワーの幾何平均は $\sqrt5P$ である。従って相関係数は $2/\sqrt5$ になる。地平線外で線形転換があれば、最初は独立な二つの揺らぎから終時刻の相関が作られる。
:::
::::

::::exercise{#exr-cosmo-primordial-fluctuations-13 level=advanced}
$m^2=2H^2$ の最小結合の試験場と、$m=0$ かつ $\xi=1/6$ の場を、同じ de Sitter 背景に置く。二点パワーが同じになる理由と、異なる背景ではどの条件が変わるかを述べよ。
:::solution
de Sitter では $R=12H^2$ なので、後者の $\xi R$ は $2H^2$ である。両場で $a^2(m^2+\xi R)-a''/a=0$ となり、同じ真空初期値なら正準モードと $\mathcal P_\varphi=k^2/(4\pi^2a^2)$ が一致する。一方、平坦 FLRW の $R=6a''/a^3$ を用いると、後者の共形結合では背景の時間発展が変わっても $a^2R/6=a''/a$ が恒等的に成り立つ。前者の一定質量で同じ相殺が成り立つには $m^2=a''/a^3$ が区間全体で一定でなければならない。従って一つの背景でパワーが一致しても、作用の結合まで同じとは言えない。
:::
::::

## 確認テスト

::::quiz
Bunch–Davies 初期条件を課す領域はどれか。

- [x] $k\gg aH$
- [ ] $k=0$
- [ ] $k\ll aH$

:::explanation
$k^2$ が曲率の項より大きく正周波数平面波が定義できる短波長側である。
:::
::::

::::quiz{answer="0.032" tol="1e-6"}
$\epsilon_H=0.002$ に対する先頭近似の $r$ は何か。
:::explanation
比は $16\epsilon_H$ なので $0.032$ である。
:::
::::

::::quiz{answer="18" tol="1e-6"}
局所型で $f_{\rm NL}=5$ の正三角形における $B_\zeta/P^2$ を求めよ。
:::explanation
三通りの積が全て $P^2$ なので $(6/5)5\cdot3=18$ である。
:::
::::

::::quiz
長波長曲率の一定モードが優勢になる条件を選べ。

- [x] 一定モードの係数が非零で、第二解の寄与が減衰する。
- [ ] 第二解が急増する。
- [ ] 常に一定モードが零である。

:::explanation
第二解 $\int\dd\eta/z_s^2$ の時間依存部が減衰する背景では、非零の一定モードを持つ解はその一定値に近づく。一定モードの係数が零なら、減衰する解だけが残る。
:::
::::

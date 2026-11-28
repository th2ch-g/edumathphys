---
title: 光の電磁波理論
description: 一様媒質中の電磁波を Maxwell 方程式から導き、波の位相、エネルギー流、運動量および光子の測定を結び付ける。
notation:
  - symbol: '\omega'
    meaning: 角周波数
  - symbol: '\lambda_0,\lambda'
    meaning: 真空中および媒質中の波長
  - symbol: '\eta'
    meaning: 媒質の波動インピーダンス
  - symbol: '\vb S'
    meaning: Poynting ベクトル
---

## この章の目標

- 一様媒質で電磁場の波動方程式を導き、発散の条件を同時に課す。
- 実場と複素振幅を使い分け、平面波と球面波の位相を読む。
- 照度、光速、放射圧を別々の測定量として計算する。
- 光子エネルギーと検出率を結び付け、干渉と単一光子検出を説明する。

## Maxwell 方程式から波動方程式へ

真空や空間的に一様なガラス中では、波を追う範囲内に自由電荷や自由電流がなくても、電場と磁場は互いの時間変化を通して伝わる。物質による応答を構成関係に入れた上で、この結合から波の速さを求めよう。

### 一様な線形等方媒質

実数の正定数 $\eps,\mu$ に対して $\vb D=\eps\vb E$、$\vb B=\mu\vb H$ とする。この式は、応答が場に比例し、向きに依存せず、観測する帯域で遅延を無視できる近似である。自由電荷密度と自由電流密度をゼロにすると、SI 系の Maxwell 方程式は
$$
\begin{aligned}
\divg\vb D&=0,&\divg\vb B&=0,\\
\rot\vb E&=-\partial_t\vb B,&
\rot\vb H&=\partial_t\vb D .
\end{aligned}
$$

:::proposition[電磁場の波動方程式]
自由電荷と自由電流がなく、$\eps,\mu$ が正の定数である領域で、Maxwell 方程式の解は次を満たす。
$$
\nabla^2\vb F-\mu\eps\partial_t^2\vb F=0,\qquad
\vb F\in\{\vb E,\vb H\},\qquad v=(\mu\eps)^{-1/2}.
\label{eq-optics-em-wave}
$$
:::

:::proof
Faraday の式の両辺に回転を取る。成分で $\epsilon_{ijk}\epsilon_{klm}=\delta_{il}\delta_{jm}-\delta_{im}\delta_{jl}$ を使うと、
$$
(\rot\rot\vb E)_i=\partial_i\partial_jE_j-\partial_j\partial_j E_i .
$$
$\eps$ が一定なので $\divg\vb E=0$ である。したがって左辺は $-\nabla^2\vb E$、右辺は $-\mu\partial_t(\rot\vb H)=-\mu\eps\partial_t^2\vb E$ となる。磁場については Ampère–Maxwell の式から $\rot\rot\vb H=\eps\partial_t(\rot\vb E)=-\eps\mu\partial_t^2\vb H$ となる。$\divg\vb H=0$ なので左辺は $-\nabla^2\vb H$ である。両式の時間二階微分の係数は $v^{-2}$ である。
:::

この導出で捨てた発散条件も、解を選ぶ際に必要である。例えば $\vb E=\hat{\vb z}f(z-vt)$ も波動方程式を満たすが、$f'\ne0$ なら $\divg\vb E=f'$ となる。このため、自由電荷ゼロという仮定とは両立しない。回転の二式から $\partial_t\divg\vb D=\partial_t\divg\vb B=0$ が従うので、初期時刻で発散条件を満たす場はその後も満たす。

:::example[進行波の相棒]
$\vb E=\hat{\vb x}f(z-vt)$ とすると $\rot\vb E=\hat{\vb y}f'$ なので、$\vb B=\hat{\vb y}f/v$ を組に選べば $-\partial_t\vb B=\hat{\vb y}f'$ となる。さらに $\rot\vb H=-\hat{\vb x}f'/(\mu v)=\eps\partial_t\vb E$ であり、両発散もゼロである。$f$ を任意の滑らかなパルスに選べるのは媒質の分散を無視したためである。
:::

:::example[静電場を加える]
一定の電場を進行波に加えても Maxwell 方程式は変わらない。定数の時間微分と空間微分がゼロだからである。単色表示で $\omega>0$ の成分を選ぶときは、この静的な成分を別に扱う。
:::

### 電流がある領域では何が変わるか
真空に電荷密度 $\rho$ と電流密度 $\vb J$ があるときは $\divg\vb E=\rho/\eps_0$ である。上と同じ回転の計算から
$$
\begin{aligned}
(\nabla^2-c^{-2}\partial_t^2)\vb E
 &=\eps_0^{-1}\grad\rho+\mu_0\partial_t\vb J,\\
(\nabla^2-c^{-2}\partial_t^2)\vb B&=-\mu_0\rot\vb J .
\end{aligned}
$$
第一式では回転の積の $\grad\divg\vb E$ が残り、第二式では磁場の発散がゼロなので電流の回転だけが残る。電荷と電流は独立に与えられない。Ampère–Maxwell の式の発散を取ると $\partial_t\rho+\divg\vb J=0$ になるためである。例えば局在電流を正弦駆動して光を出すと、電流のない離れた位置では同次方程式を解くが、波の振幅は源と境界で定まる。

:::example[電荷密度のない一様電流]
無限空間の理想化で $\vb J=J_0\cos\omega t\,\hat{\vb x}$、$\rho=0$ とする。$\vb B=0$、$\vb E=-(J_0/\eps_0\omega)\sin\omega t\,\hat{\vb x}$ を取ると Ampère–Maxwell の式の $\vb J+\eps_0\partial_t\vb E$ がゼロになり、回転と発散の残りの式も満たす。時間変化する場があっても、その項を空間伝播する波に限って考えることはできない。
:::

### 初期値から二つの場を復元する
一様な無源媒質で初期値 $\vb E(\vb r,0)=\vb E_i(\vb r)$、$\vb H(\vb r,0)=\vb H_i(\vb r)$ が十分滑らかに減衰するとする。両者が発散ゼロの場合、初期時間微分を $\partial_t\vb E(\vb r,0)=\eps^{-1}\rot\vb H_i$ として電場の波動方程式を解く。次に $\vb H(t)=\vb H_i-\mu^{-1}\int_0^t\rot\vb E(t')\dd t'$ と置く。電場の発散は波動方程式を満たし、初期値と初期微分がともにゼロなので消える。実際、空間 Fourier 変換すると、波数ごとの調和振動子が初期値ゼロで始まる。磁場の発散は定義からゼロである。残る Ampère–Maxwell の式について、残差 $\vb R=\rot\vb H-\eps\partial_t\vb E$ を時間微分すると $\partial_t\vb R=-\mu^{-1}\rot\rot\vb E-\eps\partial_t^2\vb E=0$ になる。初期時刻で $\vb R=0$ なので、全時刻で四式を満たす。従って電場だけの波動方程式を使うときは、電場の初期微分にも磁場との整合を課せばよい。

:::example[任意の横初期値を二方向へ分ける]
一次元で $\vb E=E_x(z,t)\hat{\vb x}$、$\vb H=H_y(z,t)\hat{\vb y}$ とする。$\eta=\sqrt{\mu/\eps}$、$E_\pm(z)=[E_i(z)\pm\eta H_i(z)]/2$ と置くと解は $E_x=E_+(z-vt)+E_-(z+vt)$、$H_y=\eta^{-1}[E_+(z-vt)-E_-(z+vt)]$。初期値を代入すると元の $E_i,H_i$ に一致し、各項は対応する向きの Maxwell 方程式を満たす。例えば $H_i=0$ の電場分布からは、二つの等振幅の逆向きパルスが出る。$H_i=E_i/\eta$ なら逆向きパルスは出ない。
:::

### 導電性媒質の初期値と電荷緩和

絶縁体の波動方程式に減衰項を加えると、振幅の減少を扱える。ただし、導電性媒質では電荷密度も変化するため、横波と縦方向の電荷の運動を分ける。以下では $\sigma$、$\eps$、$\mu$ を空間的に一様な正定数とし、$\vb J=\sigma\vb E$ と置く。電荷保存式と Gauss の法則をつなぐと
$$
\partial_t\rho=-\divg\vb J=-\frac{\sigma}{\eps}\rho,\qquad
\rho(\vb r,t)=\rho(\vb r,0)e^{-t/\tau_q},\qquad \tau_q=\frac{\eps}{\sigma}
$$
となる。$\tau_q$ は :term[電荷緩和時間]{en="charge relaxation time" yomi="でんかかんわじかん"} である。この証明は領域内部で成り立つ。媒質の端で表面電荷が蓄積するときは、端に流れ込む電流を別に数える。

$\rho=0$ の初期値から始めれば、その後も $\divg\vb E=0$ である。Faraday の式を回転させ、$\rot\vb H=\sigma\vb E+\eps\partial_t\vb E$ を入れると
$$
\nabla^2\vb E-\mu\sigma\partial_t\vb E-\mu\eps\partial_t^2\vb E=0
$$
を得る。波数 $k$ の横 Fourier 成分の振幅を $A(t)$ とすると $A''+2\Gamma A'+\omega_k^2 A=0$ であり、$\Gamma=\sigma/(2\eps)$、$\omega_k=k/\sqrt{\mu\eps}$ である。$A=e^{-\Gamma t}F$ と置けば $F''+(\omega_k^2-\Gamma^2)F=0$。従って $\omega_k>\Gamma$ のとき減衰しながら振動し、等号では $F$ が時間の一次式、逆の不等号では $F$ が指数関数の和になる。最後の場合も $A$ の二つの指数率 $-\Gamma\pm\sqrt{\Gamma^2-\omega_k^2}$ は $k>0$ で負である。これは実波数を固定した時間発展の分類である。外部から実周波数で駆動する問題では、空間減衰を複素波数に含める。

:::example[低い導電率で駆動を止める]
弱く導電する誘電体で、初期電荷と横モードが同時に残ったとする。両者の時間尺度を比較する。

入力は誘電率 $\epsilon\simeq 2\times 10^{-11}\,\mathrm{F/m}$、導電率 $\sigma\simeq 1\times 10^{-5}\,\mathrm{S/m}$、無損失なら得られる角周波数 $\omega_k=30000000\,\mathrm{s^{-1}}$、初期時刻からの時間 $t\simeq 2\times 10^{-6}\,\mathrm{s}$ である。

電荷密度が $1/e$ に減るまでの時間を $\eps$ と $\sigma$ から求める。

$$
\begin{aligned}
\tau_q&=\frac{\epsilon}{\sigma},\\
\tau_q/(\mathrm{s})&\simeq \frac{\left(2\times 10^{-11}\right)}{\left(1\times 10^{-5}\right)}\simeq 2\times 10^{-6}.
\end{aligned}
$$

横モードの振幅因子を取り出すと、指数率の絶対値は電荷緩和率の半分である。

$$
\begin{aligned}
\Gamma&=\frac{1}{2 \tau_q},\\
\Gamma/(\mathrm{s^{-1}})&\simeq \frac{1}{2 \cdot \left(2\times 10^{-6}\right)}=250000.
\end{aligned}
$$

真空に近い波動を仮定する前に、減衰率と無損失の角周波数を比較する。

$$
\begin{aligned}
r&=\frac{\Gamma}{\omega_k},\\
r&= \frac{\left(250000\right)}{\left(30000000\right)}=\frac{1}{120}.
\end{aligned}
$$

初期電荷密度に緩和の指数を掛ける。

$$
\begin{aligned}
q_{\mathrm{rem}}&=e^{- \frac{t}{\tau_q}},\\
q_{\mathrm{rem}}&\simeq e^{- \frac{\left(2\times 10^{-6}\right)}{\left(2\times 10^{-6}\right)}}\simeq 0.36788.
\end{aligned}
$$

不足減衰の場合はこの因子が横モードの振動の包絡線を与える。過減衰の場合は、これに非振動の $F$ が掛かる。

$$
\begin{aligned}
a_{\mathrm{env}}&=e^{- \Gamma t},\\
a_{\mathrm{env}}&\simeq e^{- \left(250000\right) \cdot \left(2\times 10^{-6}\right)}\simeq 0.60653.
\end{aligned}
$$

$r\ll1$ なので、横モードは振動を保つ。電荷密度と横モードの振幅を同じ割合で減らすモデルでは、この差を再現できない。
:::

:::example[臨界減衰で初期値を選ぶ]
波数を調整して横モードを臨界減衰に置く。初期微分を $\dot A(0)=-\Gamma A(0)$ に選んだときの残留振幅を考える。

入力は誘電率 $\epsilon\simeq 4\times 10^{-11}\,\mathrm{F/m}$、導電率 $\sigma=\frac{1}{1250}\,\mathrm{S/m}$、無損失なら得られる角周波数 $\omega_k=10000000\,\mathrm{s^{-1}}$、初期時刻からの時間 $t\simeq 1\times 10^{-7}\,\mathrm{s}$ である。

電荷密度が $1/e$ に減るまでの時間を $\eps$ と $\sigma$ から求める。

$$
\begin{aligned}
\tau_q&=\frac{\epsilon}{\sigma},\\
\tau_q/(\mathrm{s})&\simeq \frac{\left(4\times 10^{-11}\right)}{\left(\frac{1}{1250}\right)}\simeq 5\times 10^{-8}.
\end{aligned}
$$

横モードの振幅因子を取り出すと、指数率の絶対値は電荷緩和率の半分である。

$$
\begin{aligned}
\Gamma&=\frac{1}{2 \tau_q},\\
\Gamma/(\mathrm{s^{-1}})&\simeq \frac{1}{2 \cdot \left(5\times 10^{-8}\right)}=10000000.
\end{aligned}
$$

真空に近い波動を仮定する前に、減衰率と無損失の角周波数を比較する。

$$
\begin{aligned}
r&=\frac{\Gamma}{\omega_k},\\
r&= \frac{\left(10000000\right)}{\left(10000000\right)}=1.
\end{aligned}
$$

初期電荷密度に緩和の指数を掛ける。

$$
\begin{aligned}
q_{\mathrm{rem}}&=e^{- \frac{t}{\tau_q}},\\
q_{\mathrm{rem}}&\simeq e^{- \frac{\left(1\times 10^{-7}\right)}{\left(5\times 10^{-8}\right)}}\simeq 0.13534.
\end{aligned}
$$

不足減衰の場合はこの因子が横モードの振動の包絡線を与える。過減衰の場合は、これに非振動の $F$ が掛かる。

$$
\begin{aligned}
a_{\mathrm{env}}&=e^{- \Gamma t},\\
a_{\mathrm{env}}&\simeq e^{- \left(10000000\right) \cdot \left(1\times 10^{-7}\right)}\simeq 0.36788.
\end{aligned}
$$

$r=1$ で $F''=0$ となり、選んだ初期微分では $F'(0)=0$ である。そのため、ここでは $a_{\mathrm{env}}$ がそのまま振幅比になる。
:::

:::example[高導電率で小さな波数を与える]
電荷を速く緩和させる媒質で、小さな $k$ の横モードを同時に励起する。

入力は誘電率 $\epsilon\simeq 1\times 10^{-11}\,\mathrm{F/m}$、導電率 $\sigma=\frac{1}{250}\,\mathrm{S/m}$、無損失なら得られる角周波数 $\omega_k=20000000\,\mathrm{s^{-1}}$、初期時刻からの時間 $t\simeq 1\times 10^{-8}\,\mathrm{s}$ である。

電荷密度が $1/e$ に減るまでの時間を $\eps$ と $\sigma$ から求める。

$$
\begin{aligned}
\tau_q&=\frac{\epsilon}{\sigma},\\
\tau_q/(\mathrm{s})&\simeq \frac{\left(1\times 10^{-11}\right)}{\left(\frac{1}{250}\right)}\simeq 2.5\times 10^{-9}.
\end{aligned}
$$

横モードの振幅因子を取り出すと、指数率の絶対値は電荷緩和率の半分である。

$$
\begin{aligned}
\Gamma&=\frac{1}{2 \tau_q},\\
\Gamma/(\mathrm{s^{-1}})&\simeq \frac{1}{2 \cdot \left(2.5\times 10^{-9}\right)}=200000000.
\end{aligned}
$$

真空に近い波動を仮定する前に、減衰率と無損失の角周波数を比較する。

$$
\begin{aligned}
r&=\frac{\Gamma}{\omega_k},\\
r&= \frac{\left(200000000\right)}{\left(20000000\right)}=10.
\end{aligned}
$$

初期電荷密度に緩和の指数を掛ける。

$$
\begin{aligned}
q_{\mathrm{rem}}&=e^{- \frac{t}{\tau_q}},\\
q_{\mathrm{rem}}&\simeq e^{- \frac{\left(1\times 10^{-8}\right)}{\left(2.5\times 10^{-9}\right)}}\simeq 0.018316.
\end{aligned}
$$

不足減衰の場合はこの因子が横モードの振動の包絡線を与える。過減衰の場合は固有解ごとに減衰率が異なる。ここでは遅い固有解だけを励起した場合を計算する。

$$
\begin{aligned}
\frac{A_{\rm slow}(t)}{A_{\rm slow}(0)}&=\exp[-(\Gamma-\sqrt{\Gamma^2-\omega_k^2})t],\\
\frac{A_{\rm slow}(10\,\mathrm{ns})}{A_{\rm slow}(0)}&=\exp[-2(1-\sqrt{0.99})]\simeq 0.99002.
\end{aligned}
$$

この時間では電荷密度の約 $1.8\%$ が残るのに対し、遅い横モードの振幅は約 $99\%$ が残る。一般の初期値では、より速く減衰する横モードもこの解に加わる。
:::

## 平面波と球面波の複素表示

### 位相と横波条件

本サイトでは単色場を $\vb E(\vb r,t)=\Re[\vb{\mathcal E}(\vb r)e^{-i\omega t}]$ と書く。周波数を $\nu$、角周波数を $\omega=2\pi\nu>0$ とする。この規約では時間微分は $-i\omega$ 倍である。:eqref[eq-optics-em-wave] に代入すると $\nabla^2\vb{\mathcal E}+k^2\vb{\mathcal E}=0$、$k=\omega\sqrt{\eps\mu}$ を得る。これを :term[Helmholtz 方程式]{en="Helmholtz equation" yomi="へるむほるつほうていしき"} と呼ぶ。

:term[平面波]{en="plane wave" yomi="へいめんは"} では $\vb{\mathcal E}=\vb E_0 e^{i\vb k\cdot\vb r}$ とする。空間微分は $i\vb k$ 倍になり、四つの Maxwell 方程式は
$$
\abs{\vb k}^2=\omega^2\eps\mu,\qquad
\vb k\cdot\vb E_0=0,\qquad
\vb H_0=\frac{\vb k\times\vb E_0}{\omega\mu}
=\frac{\hat{\vb k}\times\vb E_0}{\eta},\qquad
\eta=\sqrt{\mu/\eps}
\label{eq-optics-plane}
$$
と整理できる。実際、Faraday の式が三番目の等式を与え、この $\vb H_0$ を Ampère–Maxwell の式へ代入すると $\vb k\times(\vb k\times\vb E_0)=-k^2\vb E_0$ となる。電場と磁場はいずれも進行方向と垂直である。$\eta$ を :term[波動インピーダンス]{en="wave impedance" yomi="はどういんぴーだんす"} と呼び、単位は $\Omega$ である。真空では約 $377\,\Omega$ となる。この平面波の議論は :eqref[eq-em-plane-wave] と同じである。

:::example[位相面の移動]
$\cos(kz-\omega t+\phi)$ の山は $kz-\omega t+\phi=2\pi m$ を満たす。$m$ を固定して微分すると $\dd z/\dd t=\omega/k>0$ である。反対に $\cos(-kz-\omega t+\phi)$ では山が $-z$ 方向へ進む。$e^{-i\omega t}$ を採用したとき、空間指数の符号で向きを区別できる。
:::

### 球対称なスカラー解

波面の曲率を残したいときには :term[球面波]{en="spherical wave" yomi="きゅうめんは"} を使う。$\vb r$ のノルムを $r$ とし、動径だけに依存する複素振幅を $U(r)$ とする。$\partial_i r=x_i/r$ から $\nabla^2 U=U''+2U'/r=(rU)''/r$ である。したがって $r>0$ で
$$
U(r)=\frac{A e^{ikr}+B e^{-ikr}}{r}.
\label{eq-optics-spherical}
$$
$A$ の項は外向き、$B$ の項は内向きである。外向き波の遠方では、放射成分の $1/r$ の振幅が逆二乗の照度に対応する。源の外側で $\vb E=\hat{\vb x}U(r)$ と機械的に置くと $\divg\vb E=\partial_xU$ が一般にゼロにならない。電磁場を作るには横波条件を課す必要がある。例えば任意の定ベクトル $\vb a$ を用いて $\vb{\mathcal E}=\rot(\vb aU)$ と置けば、回転と Laplace 演算子が可換なので Helmholtz 方程式を満たし、$\divg\vb{\mathcal E}=0$ である。対応する磁場は $\vb{\mathcal H}=(i\omega\mu)^{-1}\rot\vb{\mathcal E}$ で求まる。

:::example[遠方の局所平面波]
外向き解の空間微分は $U'(r)=(ik-1/r)U(r)$ である。$kr\gg1$ なら $1/r$ の項は $k$ に比べて小さい。このとき $\vb{\mathcal E}\simeq ik(\hat{\vb r}\times\vb a)U$ であり、方向を固定した小領域では :eqref[eq-optics-plane] を進行方向 $\hat{\vb r}$ に適用できる。
:::

### 球面波を小さな平面で観測する
軸からの横距離を $\rho_\perp$、源から観測面までを $z$ とする。$r=\sqrt{z^2+\rho_\perp^2}=z+\rho_\perp^2/(2z)+O(\rho_\perp^4/z^3)$ なので、源の横位置が軸上なら平面に届く球面波の位相は先頭の補正で $k\rho_\perp^2/(2z)$ である。振幅の差が相対的に小さいことと、位相の差が 1 ラジアン未満であることは別々の条件である。後者には $k\rho_\perp^2/(2z)\ll1$ が要る。この区別が [回折理論](/physics/optics/diffraction/) で現れる。

:::example[振幅は一様でも位相は曲がる]
$\lambda=500\,\mathrm{nm}$、$z=1\,\mathrm m$、観測面の半径 $1\,\mathrm{mm}$ とすると、外縁での振幅の相対差は約 $5\times10^{-7}$ である。一方の位相差は $2\pi$ 程度であるから、面全体を同相の平面波で置き換えると干渉の計算を誤る。
:::

### 外向き条件、点源と有限の波面

点源の解を選ぶときは、無限遠から入る波を指定しないという条件がいる。スカラー場 $U$ に対する :term[Sommerfeld の放射条件]{en="Sommerfeld radiation condition" yomi="ぞんまーふぇるとのほうしゃじょうけん"} は $r(\partial_rU-ikU)\to0$ である。:eqref[eq-optics-spherical] の外向き成分では左辺が $-Ae^{ikr}/r$ になり、内向き成分では $-2ikBe^{-ikr}-Be^{-ikr}/r$ となる。従って球対称解でこの条件が成り立つのは $B=0$ の場合である。

点源を分布として規格化するには $G=e^{ikr}/(4\pi r)$ の微分を原点のまわりで積分する。半径 $a$ の球で $4\pi a^2\partial_r G(a)=e^{ika}(ika-1)\to-1$、$k^2G$ の体積積分は $O(a^2)$ である。より正確には、試験関数 $\varphi$ を掛けて原点を除く部分で Green の公式を使うと、原点側の境界で領域の外向き法線は $-\hat{\vb r}$ である。従って $\int G\nabla^2\varphi\,\dd V$ の境界項は $\int_{r=a}(-G\partial_r\varphi+\varphi\partial_rG)\,\dd A$ となり、前項は $O(a)$、後項は $-\varphi(0)$ に収束する。従って $(\nabla^2+k^2)G=-\delta$ である。この負号が [回折理論](/physics/optics/diffraction/) の境界積分の符号を定める。

球面波を平面上に制限する場合は振幅と位相を別々に比較する。$u=\rho_\perp/z$ と置くと、二次位相の余りは $kz[\sqrt{1+u^2}-1-u^2/2]$ である。平方根の有理化を二度使えば、この式は $-kz u^4/[2(\sqrt{1+u^2}+1)^2]$ と書ける。そのため、二次位相モデルの絶対位相誤差は $kz u^4/8$ 以下である。近軸応用では、小角であることに加えてこの誤差を確かめる。

:::example[小さな受光面]
小さな検出器に外向き球面波を入れ、正面からの一様な平面波で置き換えられるか判定する。

入力は媒質中の波長 $\lambda\simeq 6\times 10^{-7}\,\mathrm m$、軸上の距離 $z=2\,\mathrm m$、評価面の半径 $a\simeq 1\times 10^{-4}\,\mathrm m$ である。

幾何学的な距離は $zg$ である。外縁の距離が軸上とどれだけ違うかを測る。

$$
\begin{aligned}
g&=\sqrt{\frac{a^{2}}{z^{2}} + 1},\\
g&\simeq \sqrt{\frac{\left(\left(1\times 10^{-4}\right)\right)^{2}}{\left(2\right)^{2}} + 1}\simeq 1.
\end{aligned}
$$

同じ球面波の二点で $1/r$ の比を取る。

$$
\begin{aligned}
A_{\mathrm{edge}}/A_{\mathrm{axis}}&=\frac{1}{g},\\
A_{\mathrm{edge}}/A_{\mathrm{axis}}&\simeq \frac{1}{\left(1.000000001\right)}\simeq 1.
\end{aligned}
$$

位相の二次近似を用いて、面全体を同相と見なせるかを調べる。

$$
\begin{aligned}
\phi_2&=\frac{\pi a^{2}}{\lambda z},\\
\phi_2/(\mathrm{rad})&\simeq \frac{\pi \cdot \left(\left(1\times 10^{-4}\right)\right)^{2}}{\left(6\times 10^{-7}\right) \cdot \left(2\right)}\simeq 0.02618.
\end{aligned}
$$

二次項までの位相を使うなら、残った四次以上の誤差をこの値で見積もる。

$$
\begin{aligned}
\lvert\delta\phi\rvert&=\frac{\pi a^{4}}{\lambda z^{3} \left(g + 1\right)^{2}},\\
\lvert\delta\phi\rvert/(\mathrm{rad})&\simeq \frac{\pi \cdot \left(\left(1\times 10^{-4}\right)\right)^{4}}{\left(6\times 10^{-7}\right) \cdot \left(2\right)^{3} \cdot \left(\left(1.000000001\right) + 1\right)^{2}}\simeq 1.6362\times 10^{-11}.
\end{aligned}
$$

軸上から外縁までの二次位相を周期の数で表し、干渉に効く位相の変化を読む。

$$
\begin{aligned}
N&=\frac{\phi_2}{2 \pi},\\
N&\simeq \frac{\left(0.02617993878\right)}{2 \cdot \pi}=\frac{1}{240}.
\end{aligned}
$$

振幅差と位相差がともに小さい。従って正面付近の複素振幅を一定とする近似が使える。
:::

:::example[位相だけを保存する面]
開口内をほぼ一様に照らしながら、波面の曲率をレンズで測ることを考える。

入力は媒質中の波長 $\lambda\simeq 6\times 10^{-7}\,\mathrm m$、軸上の距離 $z=1\,\mathrm m$、評価面の半径 $a=\frac{3}{1000}\,\mathrm m$ である。

幾何学的な距離は $zg$ である。外縁の距離が軸上とどれだけ違うかを測る。

$$
\begin{aligned}
g&=\sqrt{\frac{a^{2}}{z^{2}} + 1},\\
g&= \sqrt{\frac{\left(\frac{3}{1000}\right)^{2}}{\left(1\right)^{2}} + 1}\simeq 1.0000045.
\end{aligned}
$$

同じ球面波の二点で $1/r$ の比を取る。

$$
\begin{aligned}
A_{\mathrm{edge}}/A_{\mathrm{axis}}&=\frac{1}{g},\\
A_{\mathrm{edge}}/A_{\mathrm{axis}}&\simeq \frac{1}{\left(1.0000045\right)}\simeq 0.9999955.
\end{aligned}
$$

位相の二次近似を用いて、面全体を同相と見なせるかを調べる。

$$
\begin{aligned}
\phi_2&=\frac{\pi a^{2}}{\lambda z},\\
\phi_2/(\mathrm{rad})&\simeq \frac{\pi \cdot \left(\frac{3}{1000}\right)^{2}}{\left(6\times 10^{-7}\right) \cdot \left(1\right)}\simeq 47.124.
\end{aligned}
$$

二次項までの位相を使うなら、残った四次以上の誤差をこの値で見積もる。

$$
\begin{aligned}
\lvert\delta\phi\rvert&=\frac{\pi a^{4}}{\lambda z^{3} \left(g + 1\right)^{2}},\\
\lvert\delta\phi\rvert/(\mathrm{rad})&\simeq \frac{\pi \cdot \left(\frac{3}{1000}\right)^{4}}{\left(6\times 10^{-7}\right) \cdot \left(1\right)^{3} \cdot \left(\left(1.0000045\right) + 1\right)^{2}}\simeq 1.0603\times 10^{-4}.
\end{aligned}
$$

軸上から外縁までの二次位相を周期の数で表し、干渉に効く位相の変化を読む。

$$
\begin{aligned}
N&=\frac{\phi_2}{2 \pi},\\
N&\simeq \frac{\left(47.1238898\right)}{2 \cdot \pi}=\frac{15}{2}.
\end{aligned}
$$

振幅比はほぼ $1$ だが、位相は数周期を越えて変化する。四次以上の位相誤差は小さいため、二次位相を残したモデルが使える。
:::

:::example[近軸でも位相の四次項が必要な面]
同じ程度の開口角で、距離と波長の比を大きく取る。振幅近似と位相近似を別々に判定する。

入力は媒質中の波長 $\lambda\simeq 5\times 10^{-7}\,\mathrm m$、軸上の距離 $z=4\,\mathrm m$、評価面の半径 $a=\frac{1}{10}\,\mathrm m$ である。

幾何学的な距離は $zg$ である。外縁の距離が軸上とどれだけ違うかを測る。

$$
\begin{aligned}
g&=\sqrt{\frac{a^{2}}{z^{2}} + 1},\\
g&= \sqrt{\frac{\left(\frac{1}{10}\right)^{2}}{\left(4\right)^{2}} + 1}\simeq 1.0003.
\end{aligned}
$$

同じ球面波の二点で $1/r$ の比を取る。

$$
\begin{aligned}
A_{\mathrm{edge}}/A_{\mathrm{axis}}&=\frac{1}{g},\\
A_{\mathrm{edge}}/A_{\mathrm{axis}}&\simeq \frac{1}{\left(1.000312451\right)}\simeq 0.99969.
\end{aligned}
$$

位相の二次近似を用いて、面全体を同相と見なせるかを調べる。

$$
\begin{aligned}
\phi_2&=\frac{\pi a^{2}}{\lambda z},\\
\phi_2/(\mathrm{rad})&\simeq \frac{\pi \cdot \left(\frac{1}{10}\right)^{2}}{\left(5\times 10^{-7}\right) \cdot \left(4\right)}\simeq 15708.
\end{aligned}
$$

二次項までの位相を使うなら、残った四次以上の誤差をこの値で見積もる。

$$
\begin{aligned}
\lvert\delta\phi\rvert&=\frac{\pi a^{4}}{\lambda z^{3} \left(g + 1\right)^{2}},\\
\lvert\delta\phi\rvert/(\mathrm{rad})&\simeq \frac{\pi \cdot \left(\frac{1}{10}\right)^{4}}{\left(5\times 10^{-7}\right) \cdot \left(4\right)^{3} \cdot \left(\left(1.000312451\right) + 1\right)^{2}}\simeq 2.4536.
\end{aligned}
$$

軸上から外縁までの二次位相を周期の数で表し、干渉に効く位相の変化を読む。

$$
\begin{aligned}
N&=\frac{\phi_2}{2 \pi},\\
N&\simeq \frac{\left(15707.96327\right)}{2 \cdot \pi}=2500.
\end{aligned}
$$

開口角は約 $0.025\,\mathrm{rad}$ であるが、位相の許容誤差を $0.1\,\mathrm{rad}$ とすると四次以上を捨てられない。振幅比だけを判定に使うと位相を扱い損ねる。
:::

## Poynting ベクトルと放射照度

### 場のエネルギーの連続の式

照度を電場だけで表せるのは進行平面波などに限られる。反射波や定在波も含めるため、まずエネルギーの流れを場から定義する。ベクトル恒等式と Maxwell 方程式より
$$
\begin{aligned}
\divg(\vb E\times\vb H)
&=\vb H\cdot(\rot\vb E)-\vb E\cdot(\rot\vb H)\\
&=-\vb H\cdot\partial_t\vb B-\vb E\cdot\partial_t\vb D-\vb E\cdot\vb J_{\mathrm f}.
\end{aligned}
$$
時間独立な $\eps,\mu$ を用いると
$$
\partial_t u+\divg\vb S=-\vb J_{\mathrm f}\cdot\vb E,\qquad
u=\tfrac12(\eps\abs{\vb E}^2+\mu\abs{\vb H}^2),\qquad
\vb S=\vb E\times\vb H.
\label{eq-optics-poynting}
$$
$\vb S$ が :term[Poynting ベクトル]{en="Poynting vector" yomi="ぽいんてぃんぐべくとる"} であり、単位は $\mathrm{W/m^2}$ である。領域で積分すると、エネルギーの変化率は境界を通る流出と電流が受ける仕事率の和の負値に等しい。吸収に分散が伴う媒質では、この瞬時の $u$ をそのまま使えない。

### 時間平均と受光面

$\Re a=(a+a^*)/2$ を使って積を開くと、一周期平均では $e^{\pm2i\omega t}$ を含む項だけが消える。したがって
$$
\langle\vb S\rangle=\tfrac12\Re(\vb{\mathcal E}\times\vb{\mathcal H}^{*}),\qquad
I=\langle\vb S\rangle\cdot\hat{\vb n}.
\label{eq-optics-irradiance}
$$
$\hat{\vb n}$ は受光面へ入る向きの法線である。この単位面積当たりの入射パワーを :term[放射照度]{en="irradiance" yomi="ほうしゃしょうど"} という。正面から入る平面波では :eqref[eq-optics-plane] によって $I=\abs{\vb E_0}^2/(2\eta)$ となる（:eqref[eq-em-wave-intensity]）。$\abs{\vb E_0}^2=\vb E_0\cdot\vb E_0^*$ は複素振幅のノルムの二乗である。直線偏光なら $\abs{\vb E_0}$ は実電場の最大値に等しい。円偏光では実電場の大きさは一定で、$\abs{\vb E_0}/\sqrt 2$ である。どの偏光でも実効値 $E_{\mathrm{rms}}=\sqrt{\langle|\vb E|^2\rangle}$ を使えば $I=E_{\mathrm{rms}}^2/\eta$ となる。

:::example[検出器の傾き]
ビームと法線が $60^\circ$ をなす面の面積を $2\,\mathrm{cm^2}$、ビームに垂直な面での照度を $100\,\mathrm{W/m^2}$ とする。面全体が一様に照らされるなら入射パワーは $100(2\times10^{-4})\cos60^\circ=0.010\,\mathrm W$ である。
:::

:::example[電場の腹で流れはゼロ]
反対向きに進む波を同じ実振幅 $E_a$ で重ねると、$\mathcal E_x=2E_a\cos kz$、$\mathcal H_y=2iE_a\sin kz/\eta$ になる。積 $\mathcal E_x\mathcal H_y^*$ は純虚数である。このため電場の腹でも時間平均の流れはゼロであり、$|\mathcal E_x|^2/(2\eta)$ をそこに適用できない。
:::

### 局所のエネルギーと定在波
前述の定在波では、実場を用いて $E_x=2E_a\cos kz\cos\omega t$、$H_y=(2E_a/\eta)\sin kz\sin\omega t$ と表せる。時間平均の電気エネルギー密度は $\eps E_a^2\cos^2 kz$、磁気エネルギー密度は $\eps E_a^2\sin^2 kz$ なので、和は $\eps E_a^2$ で位置によらない。ただし瞬時には $\vb S$ が振動する。瞬時の $\vb S$ は $S_z=(E_a^2/\eta)\sin(2kz)\sin(2\omega t)$ となり、:eqref[eq-optics-poynting] でエネルギー密度の時間変化と釣り合う。

:::example[有限面で受けるパワー]
ビームを $I(r)=I_0 e^{-2r^2/w^2}$ とする。半径 $a$ の円形検出器が同軸にあるとき、受光パワーは $2\pi I_0\int_0^a r e^{-2r^2/w^2}\dd r=(\pi I_0w^2/2)(1-e^{-2a^2/w^2})$ である。$a=w$ なら全パワーの約 $86.5$ パーセントを受ける。
:::

### 複素 Poynting ベクトルの虚部

時間平均では消える交換を定量化するため、単色場の :term[複素 Poynting ベクトル]{en="complex Poynting vector" yomi="ふくそぽいんてぃんぐべくとる"} を $\vb S_c=\vb{\mathcal E}\times\vb{\mathcal H}^{*}/2$ と置く。実定数 $\eps,\mu$ の媒質で電流の複素振幅を $\vb{\mathcal J}$ とすると、
$$
\begin{aligned}
\divg\vb S_c&=\tfrac12[\vb{\mathcal H}^{*}\cdot\rot\vb{\mathcal E}-\vb{\mathcal E}\cdot\rot\vb{\mathcal H}^{*}]\\
&=-\tfrac12\vb{\mathcal E}\cdot\vb{\mathcal J}^{*}
+\frac{i\omega}{2}\bigl(\mu|\vb{\mathcal H}|^2-\eps|\vb{\mathcal E}|^2\bigr).
\end{aligned}
$$
時間因子を $e^{-i\omega t}$ としたので、ここで $\rot\vb{\mathcal H}^*=\vb{\mathcal J}^*+i\omega\eps\vb{\mathcal E}^*$ となる。自由電流がなければ実部は発散ゼロ、虚部の発散は磁気と電気の平均エネルギー密度の差の $2\omega$ 倍である。従って虚部をそのまま検出器に流入する平均パワーと解釈すると誤る。

前向きと後ろ向きの振幅の絶対値が異なる場合を調べる。$\mathcal E_x=E_+e^{ikz}+E_-e^{-ikz}$、$\eta\mathcal H_y=E_+e^{ikz}-E_-e^{-ikz}$ と置くと
$$
2\eta S_{c,z}=|E_+|^2-|E_-|^2-2i\Im(E_+E_-^*e^{2ikz})
$$
である。交差項が純虚数なので平均のエネルギー流に寄与しない。一方、$|\mathcal E_x|^2$ の交差項は実数である。電場に応答する小さな散乱体で空間的に走査するときは、平均の正味エネルギー流と違う曲線が得られる。

:::example[腹があってもパワーは流れない]
真空近似で反射して戻る波を重ねる。腹に相当する位置を取る。

入力は前向き電場の振幅 $E_+=30\,\mathrm{V/m}$、後ろ向き電場の振幅 $E_-=30\,\mathrm{V/m}$、波動インピーダンス $\eta=377\,\Omega$、二波の相対位相 $\varphi=0\,\mathrm{rad}$ である。

$$
\begin{aligned}
\lvert\mathcal E_x\rvert^2&=E_+^{2} + 2 E_+ E_- \cos{\left(\varphi \right)} + E_-^{2},\\
\lvert\mathcal E_x\rvert^2/(\mathrm{V^2/m^2})&= \left(30\right)^{2} + 2 \cdot \left(30\right) \cdot \left(30\right) \cdot \cos{\left(\left(0\right) \right)} + \left(30\right)^{2}=3600.
\end{aligned}
$$

逆進波は磁場に負符号で寄与するため、電場で足した交差項を引く。

$$
\begin{aligned}
\lvert\mathcal H_y\rvert^2&=\frac{E_+^{2} - 2 E_+ E_- \cos{\left(\varphi \right)} + E_-^{2}}{\eta^{2}},\\
\lvert\mathcal H_y\rvert^2/(\mathrm{A^2/m^2})&= \frac{\left(30\right)^{2} - 2 \cdot \left(30\right) \cdot \left(30\right) \cdot \cos{\left(\left(0\right) \right)} + \left(30\right)^{2}}{\left(377\right)^{2}}=0.
\end{aligned}
$$

両方向の進行波が運ぶパワーを、向きに応じて足し引きする。

$$
\begin{aligned}
\Re S_c&=\frac{E_+^{2} - E_-^{2}}{2 \eta},\\
\Re S_c/(\mathrm{W/m^2})&= \frac{\left(30\right)^{2} - \left(30\right)^{2}}{2 \cdot \left(377\right)}=0.
\end{aligned}
$$

交換を表す虚部は、二波の相対位相の正弦で定まる。

$$
\begin{aligned}
\Im S_c&=- \frac{E_+ E_- \sin{\left(\varphi \right)}}{\eta},\\
\Im S_c/(\mathrm{W/m^2})&= - \frac{\left(30\right) \cdot \left(30\right) \cdot \sin{\left(\left(0\right) \right)}}{\left(377\right)}=0.
\end{aligned}
$$

電場計で測った振幅から進行波の公式を使った場合の値も比較する。

$$
\begin{aligned}
I_{\mathrm{est}}&=\frac{\lvert\mathcal E_x\rvert^2}{2 \eta},\\
I_{\mathrm{est}}/(\mathrm{W/m^2})&= \frac{\left(3600\right)}{2 \cdot \left(377\right)}=\frac{1800}{377}.
\end{aligned}
$$

ここでは電場だけから求めた $I_{\mathrm{est}}$ が非零である。定在波の平均流はゼロである。
:::

:::example[部分反射がある位置]
進行波より小さな逆進波がある。二波の相対位相が四分の一周期の位置で流れを評価する。

入力は前向き電場の振幅 $E_+=30\,\mathrm{V/m}$、後ろ向き電場の振幅 $E_-=15\,\mathrm{V/m}$、波動インピーダンス $\eta=377\,\Omega$、二波の相対位相 $\varphi\simeq 1.5708\,\mathrm{rad}$ である。

$$
\begin{aligned}
\lvert\mathcal E_x\rvert^2&=E_+^{2} + 2 E_+ E_- \cos{\left(\varphi \right)} + E_-^{2},\\
\lvert\mathcal E_x\rvert^2/(\mathrm{V^2/m^2})&\simeq \left(30\right)^{2} + 2 \cdot \left(30\right) \cdot \left(15\right) \cdot \cos{\left(\left(1.570796327\right) \right)} + \left(15\right)^{2}=1125.
\end{aligned}
$$

逆進波は磁場に負符号で寄与するため、電場で足した交差項を引く。

$$
\begin{aligned}
\lvert\mathcal H_y\rvert^2&=\frac{E_+^{2} - 2 E_+ E_- \cos{\left(\varphi \right)} + E_-^{2}}{\eta^{2}},\\
\lvert\mathcal H_y\rvert^2/(\mathrm{A^2/m^2})&\simeq \frac{\left(30\right)^{2} - 2 \cdot \left(30\right) \cdot \left(15\right) \cdot \cos{\left(\left(1.570796327\right) \right)} + \left(15\right)^{2}}{\left(377\right)^{2}}\simeq 0.0079153.
\end{aligned}
$$

両方向の進行波が運ぶパワーを、向きに応じて足し引きする。

$$
\begin{aligned}
\Re S_c&=\frac{E_+^{2} - E_-^{2}}{2 \eta},\\
\Re S_c/(\mathrm{W/m^2})&= \frac{\left(30\right)^{2} - \left(15\right)^{2}}{2 \cdot \left(377\right)}=\frac{675}{754}.
\end{aligned}
$$

交換を表す虚部は、二波の相対位相の正弦で定まる。

$$
\begin{aligned}
\Im S_c&=- \frac{E_+ E_- \sin{\left(\varphi \right)}}{\eta},\\
\Im S_c/(\mathrm{W/m^2})&\simeq - \frac{\left(30\right) \cdot \left(15\right) \cdot \sin{\left(\left(1.570796327\right) \right)}}{\left(377\right)}=- \frac{450}{377}.
\end{aligned}
$$

電場計で測った振幅から進行波の公式を使った場合の値も比較する。

$$
\begin{aligned}
I_{\mathrm{est}}&=\frac{\lvert\mathcal E_x\rvert^2}{2 \eta},\\
I_{\mathrm{est}}/(\mathrm{W/m^2})&= \frac{\left(1125\right)}{2 \cdot \left(377\right)}=\frac{1125}{754}.
\end{aligned}
$$

$\Re S_c$ と $\Im S_c$ がともに残る。前後の波のエネルギーを区別すると、虚部の値を吸収に対応させる必要がないことが分かる。
:::

:::example[電場振幅の小さな位置]
低インピーダンスの媒質で入射と反射の電場が打ち消し合う位置を選ぶ。

入力は前向き電場の振幅 $E_+=20\,\mathrm{V/m}$、後ろ向き電場の振幅 $E_-=16\,\mathrm{V/m}$、波動インピーダンス $\eta=250\,\Omega$、二波の相対位相 $\varphi=\pi\,\mathrm{rad}$ である。

$$
\begin{aligned}
\lvert\mathcal E_x\rvert^2&=E_+^{2} + 2 E_+ E_- \cos{\left(\varphi \right)} + E_-^{2},\\
\lvert\mathcal E_x\rvert^2/(\mathrm{V^2/m^2})&\simeq \left(20\right)^{2} + 2 \cdot \left(20\right) \cdot \left(16\right) \cdot \cos{\left(\left(3.141592654\right) \right)} + \left(16\right)^{2}=16.
\end{aligned}
$$

逆進波は磁場に負符号で寄与するため、電場で足した交差項を引く。

$$
\begin{aligned}
\lvert\mathcal H_y\rvert^2&=\frac{E_+^{2} - 2 E_+ E_- \cos{\left(\varphi \right)} + E_-^{2}}{\eta^{2}},\\
\lvert\mathcal H_y\rvert^2/(\mathrm{A^2/m^2})&\simeq \frac{\left(20\right)^{2} - 2 \cdot \left(20\right) \cdot \left(16\right) \cdot \cos{\left(\left(3.141592654\right) \right)} + \left(16\right)^{2}}{\left(250\right)^{2}}\simeq 0.020736.
\end{aligned}
$$

両方向の進行波が運ぶパワーを、向きに応じて足し引きする。

$$
\begin{aligned}
\Re S_c&=\frac{E_+^{2} - E_-^{2}}{2 \eta},\\
\Re S_c/(\mathrm{W/m^2})&= \frac{\left(20\right)^{2} - \left(16\right)^{2}}{2 \cdot \left(250\right)}=\frac{36}{125}.
\end{aligned}
$$

交換を表す虚部は、二波の相対位相の正弦で定まる。

$$
\begin{aligned}
\Im S_c&=- \frac{E_+ E_- \sin{\left(\varphi \right)}}{\eta},\\
\Im S_c/(\mathrm{W/m^2})&=-\frac{20\cdot16\sin{\pi}}{250}=0.
\end{aligned}
$$

電場計で測った振幅から進行波の公式を使った場合の値も比較する。

$$
\begin{aligned}
I_{\mathrm{est}}&=\frac{\lvert\mathcal E_x\rvert^2}{2 \eta},\\
I_{\mathrm{est}}/(\mathrm{W/m^2})&= \frac{\left(16\right)}{2 \cdot \left(250\right)}=\frac{4}{125}.
\end{aligned}
$$

この場合の $I_{\mathrm{est}}$ は平均流より小さい。場の定在波成分による誤差は受光位置で符号が変わる。
:::

## 屈折率と物質中の光速

真空中の定数を $\eps_0,\mu_0$、速度を $c=(\eps_0\mu_0)^{-1/2}$ とする。境界の両側で単色波の $\omega$ が同じでも、波数 $k$ は変わる。位相を追う速さ $v_{\mathrm{ph}}$ と真空中の値を比較して
$$
n=\frac{c}{v_{\mathrm{ph}}}=\frac{ck}{\omega}=\sqrt{\eps_r\mu_r},\qquad
\lambda=\frac{2\pi}{k}=\frac{\lambda_0}{n},\qquad
\lambda_0=\frac{2\pi c}{\omega}
\label{eq-optics-index}
$$
と定義する。$n$ が :term[屈折率]{en="refractive index" yomi="くっせつりつ"} で、$\eps_r=\eps/\eps_0$、$\mu_r=\mu/\mu_0$ である。光周波数で多くの透明誘電体は $\mu_r\simeq1$ と扱えるが、磁気応答が必要な場合は $\eta=\eta_0\sqrt{\mu_r/\eps_r}$ と $n$ の両方を独立に追う。

### 位相と包絡線

狭い帯域で $k(\omega)$ が滑らかな実関数なら、中心 $\omega_0$ の近くで $k(\omega_0+\Omega)=k_c+k'_c\Omega+O(\Omega^2)$ と置ける。二次の位相誤差 $|k''_c z\Omega^2/2|$ が小さければ、周波数成分を重ねて
$$
\int a(\Omega)e^{i[k(\omega_0+\Omega)z-(\omega_0+\Omega)t]}\dd\Omega
\simeq e^{i(k_cz-\omega_0t)}\int a(\Omega)e^{-i\Omega(t-k'_c z)}\dd\Omega
$$
となる。包絡線は時間 $k'_c z$ だけ遅れ、:term[群速度]{en="group velocity" yomi="ぐんそくど"} は $v_{\mathrm g}=1/k'_c=c/(n+\omega_0 n')$ である。帯域内で吸収が大きく異なると包絡線は変形するので、ピークの移動だけで情報の速度を推定できない。媒質応答の遅延は [物質中の光と分散](/physics/optics/dispersion-absorption/) で扱う。

:::example[水中の波長]
$\lambda_0=600\,\mathrm{nm}$、$n=4/3$ なら $\lambda=450\,\mathrm{nm}$ である。真空と水中で周波数はともに $c/\lambda_0$ で約 $5.0\times10^{14}\,\mathrm{Hz}$ となる。媒質中の波長を使うときは、媒質中の位相速度 $c/n$ を $\lambda$ で割って同じ周波数を得る。
:::

:::example[群屈折率]
中心で $n=1.50$、$\omega_0 n'=0.02$ なら $v_{\mathrm g}=c/1.52$ である。真空中を同じ距離 $L$ だけ進む場合に比べて、包絡線の遅れは $0.52 L/c$ となる。
:::

### 分散式を測る
媒質を通過する距離 $L$ を固定し、真空の同じ距離と干渉させる。反射を無視すると位相差は $\varphi=(n(\omega)-1)\omega L/c$ である。従って $\varphi/\omega$ から位相屈折率、$\dd\varphi/\dd\omega$ から真空に対する群遅延を得られる。位相は干渉から $2\pi$ を法としてしか読めないので、帯域を連続走査してつながる位相を選ぶ必要がある。絶対の整数次数を決められない場合でも、位相差の変化から群遅延を求められる。

:::example[分散のある板]
厚さ $30\,\mathrm{mm}$、位相屈折率 1.50、群屈折率 1.52 の板では、真空に対する群遅延は $0.52(0.030\,\mathrm m)/c\simeq52\,\mathrm{ps}$ である。屈折率を 1.50 で一定とみなすと $50\,\mathrm{ps}$ になり、差は $2\,\mathrm{ps}$ である。
:::

### 位相屈折率から群屈折率を測り直す

真空波長で材料を記述すると $\omega\lambda_0=2\pi c$ なので、$\omega\dd n/\dd\omega=-\lambda_0\dd n/\dd\lambda_0$ である。例えば透明帯域内の経験式 $n=A+B/\lambda_0^2$ では
$$
n_{\mathrm g}=A+\frac{3B}{\lambda_0^2},\qquad
\beta_2=\frac{\dd^2 k}{\dd\omega^2}
=-\frac{\lambda_0^2}{2\pi c^2}\frac{\dd n_{\mathrm g}}{\dd\lambda_0}
=\frac{3B}{\pi c^2\lambda_0}.
$$
ここで $B$ は面積の単位を持ち、$\beta_2$ は $\mathrm{s^2/m}$ である。最後の式は経験式が有効な帯域の局所微分に限る。包絡線を単なる時間移動で近似するには、観測帯域の端 $\Omega$ における残留位相 $|\beta_2 L\Omega^2/2|$ が十分小さくなければならない。

同じ位相屈折率の板でも、$B$ が違えば包絡線の到着時刻が変わる。二つの波長で $n_j$ を測ったとき、$u_j=\lambda_j^{-2}$ と置けば $B=(n_2-n_1)/(u_2-u_1)$、$A=(n_1u_2-n_2u_1)/(u_2-u_1)$ である。波長間隔が狭いと、分母 $u_2-u_1$ が小さくなる。例えば二点の屈折率誤差が独立で標準偏差が同じ $\sigma_n$ なら、$B$ の誤差の標準偏差は $\sqrt2\sigma_n/|u_2-u_1|$ である。一方で間隔を広げすぎると、この二項だけで屈折率を表す誤差が大きくなる。差分で遅延を推定するときは両方を確認する。

:::example[狭帯域で遅延だけを測る]
透明板を一方の腕だけに挿入する。光源の帯域内で分散パラメータを一定とみなす。

入力は経験式の定数項 $A=\frac{37}{25}$、経験式の分散係数 $B\simeq 4\times 10^{-15}\,\mathrm{m^2}$、真空波長 $\lambda_0\simeq 6\times 10^{-7}\,\mathrm m$、板の厚さ $L=\frac{1}{50}\,\mathrm m$、真空の光速 $c=299792458\,\mathrm{m/s}$、中心から測った角周波数の幅 $\Omega\simeq 2\times 10^{11}\,\mathrm{rad/s}$ である。

まず単色波の波長を決める位相屈折率を求める。

$$
\begin{aligned}
n&=A + \frac{B}{\lambda_0^{2}},\\
n&\simeq \left(\frac{37}{25}\right) + \frac{\left(4\times 10^{-15}\right)}{\left(\left(6\times 10^{-7}\right)\right)^{2}}=\frac{671}{450}.
\end{aligned}
$$

包絡線を追うため、波長微分の補正を加える。

$$
\begin{aligned}
n_{\mathrm g}&=A + \frac{3 B}{\lambda_0^{2}},\\
n_{\mathrm g}&\simeq \left(\frac{37}{25}\right) + \frac{3 \cdot \left(4\times 10^{-15}\right)}{\left(\left(6\times 10^{-7}\right)\right)^{2}}=\frac{227}{150}.
\end{aligned}
$$

真空での移動時間 $L/c$ を、媒質での時間 $n_{\mathrm g}L/c$ から引く。

$$
\begin{aligned}
\Delta t&=\frac{L \left(n_{\mathrm g} - 1\right)}{c},\\
\Delta t/(\mathrm s)&= \frac{\left(\frac{1}{50}\right) \cdot \left(\left(\frac{227}{150}\right) - 1\right)}{\left(299792458\right)}\simeq 3.4246\times 10^{-11}.
\end{aligned}
$$

時間移動から外れる量は、さらに群屈折率を周波数で微分して求める。

$$
\begin{aligned}
\beta_2&=\frac{3 B}{\pi c^{2} \lambda_0},\\
\beta_2/(\mathrm{s^2/m})&\simeq \frac{3 \cdot \left(4\times 10^{-15}\right)}{\pi \cdot \left(299792458\right)^{2} \cdot \left(6\times 10^{-7}\right)}\simeq 7.0834\times 10^{-26}.
\end{aligned}
$$

包絡線の二次位相を、スペクトルの端で評価する。

$$
\begin{aligned}
\delta\phi_2&=\frac{L \beta_2 \Omega^{2}}{2},\\
\delta\phi_2/(\mathrm{rad})&\simeq \frac{\left(\frac{1}{50}\right) \cdot \left(7.083350254\times 10^{-26}\right) \cdot \left(\left(2\times 10^{11}\right)\right)^{2}}{2}\simeq 2.8333\times 10^{-5}.
\end{aligned}
$$

帯域端の二次位相が十分小さい。ピークの時間差を群遅延として扱える。
:::

:::example[波長を短くした板]
短波長側で、前例と同じ経験式が有効とする。厚さを増やしたときに、位相屈折率で代用できるかを調べる。

入力は経験式の定数項 $A=\frac{37}{25}$、経験式の分散係数 $B\simeq 4\times 10^{-15}\,\mathrm{m^2}$、真空波長 $\lambda_0\simeq 4\times 10^{-7}\,\mathrm m$、板の厚さ $L=\frac{3}{50}\,\mathrm m$、真空の光速 $c=299792458\,\mathrm{m/s}$、中心から測った角周波数の幅 $\Omega\simeq 2\times 10^{11}\,\mathrm{rad/s}$ である。

まず単色波の波長を決める位相屈折率を求める。

$$
\begin{aligned}
n&=A + \frac{B}{\lambda_0^{2}},\\
n&\simeq \left(\frac{37}{25}\right) + \frac{\left(4\times 10^{-15}\right)}{\left(\left(4\times 10^{-7}\right)\right)^{2}}=\frac{301}{200}.
\end{aligned}
$$

包絡線を追うため、波長微分の補正を加える。

$$
\begin{aligned}
n_{\mathrm g}&=A + \frac{3 B}{\lambda_0^{2}},\\
n_{\mathrm g}&\simeq \left(\frac{37}{25}\right) + \frac{3 \cdot \left(4\times 10^{-15}\right)}{\left(\left(4\times 10^{-7}\right)\right)^{2}}=\frac{311}{200}.
\end{aligned}
$$

真空での移動時間 $L/c$ を、媒質での時間 $n_{\mathrm g}L/c$ から引く。

$$
\begin{aligned}
\Delta t&=\frac{L \left(n_{\mathrm g} - 1\right)}{c},\\
\Delta t/(\mathrm s)&= \frac{\left(\frac{3}{50}\right) \cdot \left(\left(\frac{311}{200}\right) - 1\right)}{\left(299792458\right)}\simeq 1.1108\times 10^{-10}.
\end{aligned}
$$

時間移動から外れる量は、さらに群屈折率を周波数で微分して求める。

$$
\begin{aligned}
\beta_2&=\frac{3 B}{\pi c^{2} \lambda_0},\\
\beta_2/(\mathrm{s^2/m})&\simeq \frac{3 \cdot \left(4\times 10^{-15}\right)}{\pi \cdot \left(299792458\right)^{2} \cdot \left(4\times 10^{-7}\right)}\simeq 1.0625\times 10^{-25}.
\end{aligned}
$$

包絡線の二次位相を、スペクトルの端で評価する。

$$
\begin{aligned}
\delta\phi_2&=\frac{L \beta_2 \Omega^{2}}{2},\\
\delta\phi_2/(\mathrm{rad})&\simeq \frac{\left(\frac{3}{50}\right) \cdot \left(1.062502538\times 10^{-25}\right) \cdot \left(\left(2\times 10^{11}\right)\right)^{2}}{2}\simeq 1.275\times 10^{-4}.
\end{aligned}
$$

群屈折率と位相屈折率の差は $2B/\lambda_0^2$ なので、長い板では同じ帯域でも遅延が区別しやすい。
:::

:::example[広帯域で二次位相が効く場合]
同じ透明帯域の近傍で短いパルスを使い、板による変形の有無を見積もる。

入力は経験式の定数項 $A=\frac{37}{25}$、経験式の分散係数 $B\simeq 4\times 10^{-15}\,\mathrm{m^2}$、真空波長 $\lambda_0\simeq 6\times 10^{-7}\,\mathrm m$、板の厚さ $L=\frac{1}{5}\,\mathrm m$、真空の光速 $c=299792458\,\mathrm{m/s}$、中心から測った角周波数の幅 $\Omega\simeq 2\times 10^{13}\,\mathrm{rad/s}$ である。

まず単色波の波長を決める位相屈折率を求める。

$$
\begin{aligned}
n&=A + \frac{B}{\lambda_0^{2}},\\
n&\simeq \left(\frac{37}{25}\right) + \frac{\left(4\times 10^{-15}\right)}{\left(\left(6\times 10^{-7}\right)\right)^{2}}=\frac{671}{450}.
\end{aligned}
$$

包絡線を追うため、波長微分の補正を加える。

$$
\begin{aligned}
n_{\mathrm g}&=A + \frac{3 B}{\lambda_0^{2}},\\
n_{\mathrm g}&\simeq \left(\frac{37}{25}\right) + \frac{3 \cdot \left(4\times 10^{-15}\right)}{\left(\left(6\times 10^{-7}\right)\right)^{2}}=\frac{227}{150}.
\end{aligned}
$$

真空での移動時間 $L/c$ を、媒質での時間 $n_{\mathrm g}L/c$ から引く。

$$
\begin{aligned}
\Delta t&=\frac{L \left(n_{\mathrm g} - 1\right)}{c},\\
\Delta t/(\mathrm s)&= \frac{\left(\frac{1}{5}\right) \cdot \left(\left(\frac{227}{150}\right) - 1\right)}{\left(299792458\right)}\simeq 3.4246\times 10^{-10}.
\end{aligned}
$$

時間移動から外れる量は、さらに群屈折率を周波数で微分して求める。

$$
\begin{aligned}
\beta_2&=\frac{3 B}{\pi c^{2} \lambda_0},\\
\beta_2/(\mathrm{s^2/m})&\simeq \frac{3 \cdot \left(4\times 10^{-15}\right)}{\pi \cdot \left(299792458\right)^{2} \cdot \left(6\times 10^{-7}\right)}\simeq 7.0834\times 10^{-26}.
\end{aligned}
$$

包絡線の二次位相を、スペクトルの端で評価する。

$$
\begin{aligned}
\delta\phi_2&=\frac{L \beta_2 \Omega^{2}}{2},\\
\delta\phi_2/(\mathrm{rad})&\simeq \frac{\left(\frac{1}{5}\right) \cdot \left(7.083350254\times 10^{-26}\right) \cdot \left(\left(2\times 10^{13}\right)\right)^{2}}{2}\simeq 2.8333.
\end{aligned}
$$

二次位相が無視できない。時間差だけで波形全体を記述できず、周波数成分ごとの位相を残す。
:::

## 光の運動量と放射圧

### 真空中の運動量保存

真空での場の運動量密度を $\vb g=\eps_0\vb E\times\vb B=\vb S/c^2$ とする。Lorentz 力密度は $\vb f=\rho\vb E+\vb J\times\vb B$ である。Maxwell 方程式で $\rho=\eps_0\divg\vb E$、$\vb J=(\rot\vb B)/\mu_0-\eps_0\partial_t\vb E$、$\partial_t\vb B=-\rot\vb E$ を代入すると
$$
\vb f+\partial_t\vb g
=\eps_0[(\divg\vb E)\vb E+(\rot\vb E)\times\vb E]
+\mu_0^{-1}(\rot\vb B)\times\vb B.
$$
成分で $[(\rot\vb E)\times\vb E]_i=E_j\partial_j E_i-\tfrac12\partial_i|\vb E|^2$ であり、磁場でも同じ恒等式が成り立つ。$\divg\vb B=0$ を合わせると、$f_i+\partial_tg_i=\partial_j T_{ij}$、ただし
$$
T_{ij}=\eps_0(E_iE_j-\tfrac12\delta_{ij}|\vb E|^2)
+\mu_0^{-1}(B_iB_j-\tfrac12\delta_{ij}|\vb B|^2)
$$
となる。これは :term[Maxwell 応力テンソル]{en="Maxwell stress tensor" yomi="まくすうぇるおうりょくてんそる"} である。時間平均が定常なら、物体を囲む境界で $\langle T_{ij}\rangle n_j$ を積分すれば平均力が得られる。

### 吸収と反射

進行平面波では $\langle u\rangle=I/c$、$\langle\vb g\rangle=I\hat{\vb k}/c^2$ である。したがって真空中で光エネルギー $U$ が一方向へ運ばれるとき、運動量は $U\hat{\vb k}/c$ である。正面にパワー $P$ を受ける面が、エネルギーの割合 $R$ を後方へ、$T$ を同方向の真空へ返すとする。出ていく法線運動量の率は $(T-R)P/c$ なので、面が受け取る力は
$$
F=(1+R-T)P/c .
\label{eq-optics-pressure}
$$
面積で割ったものを :term[放射圧]{en="radiation pressure" yomi="ほうしゃあつ"} という。完全吸収は $I/c$、完全反射は $2I/c$ である。完全鏡に角度 $\theta$ で入射すると、面積当たりの入射パワーは $I\cos\theta$、反転する運動量成分は入射運動量の $\cos\theta$ 倍である。従って圧力は $2I\cos^2\theta/c$ となり、接線力はゼロである。誘電体内では場だけの運動量と物質の力学的運動量を分ける必要があるため、ここでは真空中の入出射だけを比較した。

:::example[レーザーによる力]
$30\,\mathrm{mW}$ が完全鏡に正面から当たるなら $F=2(0.030)/c\simeq2.0\times10^{-10}\,\mathrm N$ である。同じ光を完全吸収面に当てると力はその半分になる。
:::

:::example[運動量の接線成分]
斜入射で完全に吸収する面では、入射光が持つ接線運動量も物体へ移る。従って面積当たりの力の大きさは $I\cos\theta/c$、うち接線成分は $I\cos\theta\sin\theta/c$ である。
:::

### 斜入射で反射と吸収が混ざる面
反射は鏡面反射、透過は同じ真空に直進し、入射角を $\theta$ とする。面積当たりの入射パワーを $I\cos\theta$ とすれば、入射と透過は同符号の法線運動量、反射は逆符号である。従って法線の圧力は $(1+R-T)I\cos^2\theta/c$。接線成分は反射でも保存されるので、面への接線力密度は $(1-R-T)I\sin\theta\cos\theta/c$ になる。光を拡散反射する面では出る光の全方向にわたって運動量を積分する必要があり、$R$ だけでは圧力を決められない。

:::example[反射で運動量が倍になる成分]
$\theta=45^\circ$、$R=1/2,T=0$ の面への法線圧力は $3I/(4c)$、接線力密度は $I/(4c)$ である。両成分の比は 3 で、合力の方向は入射光の方向と一致しない。
:::

### 定在波の応力から鏡の力を求める
理想金属の平面鏡に真空側から垂直入射する。表面で接線電場がゼロになるように、反射電場は入射の負になる。一方、反射によって進行方向も反転するため表面の磁場は $2E_a/\eta_0$ の振幅で足し合わさる。真空側の $T_{zz}$ は $-B_y^2/(2\mu_0)$ で、平均の大きさは $(2\mu_0 E_a/\eta_0)^2/(4\mu_0)=\eps_0 E_a^2$ である。入射照度 $I=E_a^2/(2\eta_0)$ に対して $\eps_0\eta_0=1/c$ なので、圧力の大きさは $2I/c$ となる。鏡の直前で平均のエネルギー流がゼロであっても、応力がゼロになるわけではない。

:::example[光のエネルギーを蓄える箱]
真空と理想鏡だけからなる箱に、軸と平行な光がエネルギー $U$ を持って往復しているとする。鏡間隔を $L$ とすれば一枚の鏡に当たる間隔は $2L/c$、到達ごとに渡す運動量は $2U/c$。したがって平均力は $U/L$ である。共振モードの光子数を一定にして鏡をゆっくり動かす場合にも、$k=m\pi/L$ なので $U\propto L^{-1}$、$-\dd U/\dd L=U/L$ となる。
:::

### 光が吸収と鏡面反射に分かれる場合

真空中の平らな面に、照度 $I$ の光が面の法線と $\theta$ の角度で入射する。法線は面へ向かう向き、接線は入射光の接線成分の向きに取る。反射は鏡面反射とし、透過も平行な真空へ出るとする。割合を $R,T$、吸収割合を $A=1-R-T$ とすると、面積当たりの入射パワーは $I\cos\theta$ である。出射光は接線方向に $(R+T)I\cos\theta\sin\theta/c$、法線方向に $(T-R)I\cos^2\theta/c$ の運動量を運び去る。入射の運動量流から引くと
$$
p_\perp=(2R+A)\frac{I\cos^2\theta}{c},\qquad
p_\parallel=A\frac{I\cos\theta\sin\theta}{c}
$$
となる。鏡面反射は接線成分を変えず、吸収は入射した成分をすべて面に渡す。散乱がある場合は、出射パワーを各方向の単位ベクトルで重み付けして積分する。

この式は物体の平均運動量の収支である。例えば光を吸収する物体が直後に熱放射でエネルギーを放出しても、放出が物体の静止系で等方なら時間平均の放射反跳はゼロである。放出が片側に偏る場合はその運動量流を引く。電磁場の力を取り出すという目的で十分大きな閉曲面を取れば、内部の吸収層の細かい構造を知らなくても、出入りの運動量を測ることで総力を求められる。

:::example[吸収膜を傾ける]
反射を抑えた吸収膜に斜めから一様な光を当てる。膜の支持部が受ける二方向の力を求める。

入力はビームに垂直な面での照度 $I=5000\,\mathrm{W/m^2}$、試料面積 $S\simeq 2\times 10^{-5}\,\mathrm{m^2}$、入射角 $\theta\simeq 0.7854\,\mathrm{rad}$、反射率 $R=0$、透過率 $T=0$、真空の光速 $c=299792458\,\mathrm{m/s}$ である。

反射と透過で出ていかない割合をエネルギー保存から求める。

$$
\begin{aligned}
A&=- R - T + 1,\\
A&= - \left(0\right) - \left(0\right) + 1=1.
\end{aligned}
$$

ビームの進行方向から見た投影面積を使う。

$$
\begin{aligned}
P&=I S \cos{\left(\theta \right)},\\
P/(\mathrm W)&\simeq \left(5000\right) \cdot \left(2\times 10^{-5}\right) \cdot \cos{\left(\left(0.7853981634\right) \right)}\simeq 0.070711.
\end{aligned}
$$

法線成分は、鏡面反射で反転する分と吸収で失われる分の和である。

$$
\begin{aligned}
F_\perp&=\frac{P \left(A + 2 R\right) \cos{\left(\theta \right)}}{c},\\
F_\perp/(\mathrm N)&\simeq \frac{\left(0.07071067812\right) \cdot \left(\left(1\right) + 2 \cdot \left(0\right)\right) \cdot \cos{\left(\left(0.7853981634\right) \right)}}{\left(299792458\right)}\simeq 1.6678\times 10^{-10}.
\end{aligned}
$$

接線方向の変化は吸収されたパワーだけが担う。

$$
\begin{aligned}
F_\parallel&=\frac{A P \sin{\left(\theta \right)}}{c},\\
F_\parallel/(\mathrm N)&\simeq \frac{\left(1\right) \cdot \left(0.07071067812\right) \cdot \sin{\left(\left(0.7853981634\right) \right)}}{\left(299792458\right)}\simeq 1.6678\times 10^{-10}.
\end{aligned}
$$

熱などに変換される光パワーも、接線力と同じ吸収率から求める。

$$
\begin{aligned}
P_{\mathrm{abs}}&=A P,\\
P_{\mathrm{abs}}/(\mathrm W)&\simeq \left(1\right) \cdot \left(0.07071067812\right)\simeq 0.070711.
\end{aligned}
$$

法線成分と接線成分は等しい。光を完全に吸収するため、力は入射方向を向く。
:::

:::example[半透鏡で熱の寄与を分ける]
吸収の小さい平行平板の前後を真空とし、透過光を入射と平行に取り出す。

入力はビームに垂直な面での照度 $I=5000\,\mathrm{W/m^2}$、試料面積 $S\simeq 2\times 10^{-5}\,\mathrm{m^2}$、入射角 $\theta\simeq 0.5236\,\mathrm{rad}$、反射率 $R=\frac{9}{20}$、透過率 $T=\frac{1}{2}$、真空の光速 $c=299792458\,\mathrm{m/s}$ である。

反射と透過で出ていかない割合をエネルギー保存から求める。

$$
\begin{aligned}
A&=- R - T + 1,\\
A&= - \left(\frac{9}{20}\right) - \left(\frac{1}{2}\right) + 1=\frac{1}{20}.
\end{aligned}
$$

ビームの進行方向から見た投影面積を使う。

$$
\begin{aligned}
P&=I S \cos{\left(\theta \right)},\\
P/(\mathrm W)&\simeq \left(5000\right) \cdot \left(2\times 10^{-5}\right) \cdot \cos{\left(\left(0.5235987756\right) \right)}\simeq 0.086603.
\end{aligned}
$$

法線成分は、鏡面反射で反転する分と吸収で失われる分の和である。

$$
\begin{aligned}
F_\perp&=\frac{P \left(A + 2 R\right) \cos{\left(\theta \right)}}{c},\\
F_\perp/(\mathrm N)&\simeq \frac{\left(0.08660254038\right) \cdot \left(\left(\frac{1}{20}\right) + 2 \cdot \left(\frac{9}{20}\right)\right) \cdot \cos{\left(\left(0.5235987756\right) \right)}}{\left(299792458\right)}\simeq 2.3766\times 10^{-10}.
\end{aligned}
$$

接線方向の変化は吸収されたパワーだけが担う。

$$
\begin{aligned}
F_\parallel&=\frac{A P \sin{\left(\theta \right)}}{c},\\
F_\parallel/(\mathrm N)&\simeq \frac{\left(\frac{1}{20}\right) \cdot \left(0.08660254038\right) \cdot \sin{\left(\left(0.5235987756\right) \right)}}{\left(299792458\right)}\simeq 7.2219\times 10^{-12}.
\end{aligned}
$$

熱などに変換される光パワーも、接線力と同じ吸収率から求める。

$$
\begin{aligned}
P_{\mathrm{abs}}&=A P,\\
P_{\mathrm{abs}}/(\mathrm W)&\simeq \left(\frac{1}{20}\right) \cdot \left(0.08660254038\right)\simeq 0.0043301.
\end{aligned}
$$

法線力には反射が大きく効くが、接線力と熱は吸収の割合に従って小さくなる。
:::

:::example[同じ光を反射膜で受ける]
吸収を無視できる鏡で、角度を大きく取ったときの力を見積もる。

入力はビームに垂直な面での照度 $I=5000\,\mathrm{W/m^2}$、試料面積 $S\simeq 2\times 10^{-5}\,\mathrm{m^2}$、入射角 $\theta\simeq 1.0472\,\mathrm{rad}$、反射率 $R=1$、透過率 $T=0$、真空の光速 $c=299792458\,\mathrm{m/s}$ である。

反射と透過で出ていかない割合をエネルギー保存から求める。

$$
\begin{aligned}
A&=- R - T + 1,\\
A&= - \left(1\right) - \left(0\right) + 1=0.
\end{aligned}
$$

ビームの進行方向から見た投影面積を使う。

$$
\begin{aligned}
P&=I S \cos{\left(\theta \right)},\\
P/(\mathrm W)&\simeq \left(5000\right) \cdot \left(2\times 10^{-5}\right) \cdot \cos{\left(\left(1.047197551\right) \right)}=\frac{1}{20}.
\end{aligned}
$$

法線成分は、鏡面反射で反転する分と吸収で失われる分の和である。

$$
\begin{aligned}
F_\perp&=\frac{P \left(A + 2 R\right) \cos{\left(\theta \right)}}{c},\\
F_\perp/(\mathrm N)&\simeq \frac{\left(\frac{1}{20}\right) \cdot \left(\left(0\right) + 2 \cdot \left(1\right)\right) \cdot \cos{\left(\left(1.047197551\right) \right)}}{\left(299792458\right)}\simeq 1.6678\times 10^{-10}.
\end{aligned}
$$

接線方向の変化は吸収されたパワーだけが担う。

$$
\begin{aligned}
F_\parallel&=\frac{A P \sin{\left(\theta \right)}}{c},\\
F_\parallel/(\mathrm N)&\simeq \frac{\left(0\right) \cdot \left(\frac{1}{20}\right) \cdot \sin{\left(\left(1.047197551\right) \right)}}{\left(299792458\right)}=0.
\end{aligned}
$$

熱などに変換される光パワーも、接線力と同じ吸収率から求める。

$$
\begin{aligned}
P_{\mathrm{abs}}&=A P,\\
P_{\mathrm{abs}}/(\mathrm W)&= \left(0\right) \cdot \left(\frac{1}{20}\right)=0.
\end{aligned}
$$

熱と接線力はゼロである。面に垂直入射したときに比べて法線力は $\cos^2\theta=1/4$ 倍になる。
:::

## 光子と光の二重性

### エネルギーの一量子と検出率

放射を弱くしても、十分な感度の検出器では出力が連続的に小さくなるのでなく、個々の検出事象の頻度が下がる。真空中の一つの電磁場モードを量子化すると、場の Hamiltonian は調和振動子になり、隣り合う準位の差は $\hbar\omega$ である。この量子を :term[光子]{en="photon" yomi="こうし"} と呼ぶ。エネルギー $E_\gamma$ と真空中の運動量 $\vb p_\gamma$ は
$$
E_\gamma=\hbar\omega=h\nu=hc/\lambda_0,\qquad
\vb p_\gamma=\hbar\vb k .
\label{eq-optics-photon}
$$
後者は平面波モードの空間並進の位相から読める。モード関数 $e^{i\vb k\cdot\vb r}$ に並進の生成子 $-i\hbar\grad$ を作用させると固有値は $\hbar\vb k$ であり、そのモードの一量子は同じ運動量を持つ。真空では大きさが $E_\gamma/c$ に一致する。熱平衡でのモードの占有数は :ref[def-stm-photons-phonons] で記述する。

単色光の平均パワーが $P$、一光子に対する検出効率が $\eta_{\mathrm d}$ なら、飽和や暗計数を無視した期待検出率は $\eta_{\mathrm d}P/E_\gamma$ である。この $\eta_{\mathrm d}$ は無次元で、波動インピーダンス $\eta$ とは異なる。

:::example[光電子の阻止電圧]
一個の光子で金属から電子を取り出すと仮定する。仕事関数を $\Phi$、最大運動エネルギーを $K_{\max}$ とすると、エネルギー収支は $K_{\max}=h\nu-\Phi$ である。$h\nu=3\,\mathrm{eV}$、$\Phi=2\,\mathrm{eV}$ なら、電子を阻止する電位差の大きさは $K_{\max}/e=1\,\mathrm V$ である。周波数を固定して照度を上げると、弱い光の範囲では電子の数が増える。
:::

### 経路振幅と検出確率

二つの経路が同じ検出結果を与え、経路を識別する記録が残らないとする。量子状態の線形性により単一光子の検出振幅は $a_1+a_2$ であり、検出確率は $|a_1|^2+|a_2|^2+2\Re(a_1 a_2^*)$ となる。経路標識がそれぞれ正規化された状態 $\ket{m_1},\ket{m_2}$ に残る場合は、ベクトル $a_1\ket{m_1}+a_2\ket{m_2}$ の二乗ノルムを取るので、交差項は $2\Re(a_1a_2^*\braket{m_2}{m_1})$ となる。直交する標識なら縞は消える。検出事象は点に局在するが、多数回の試行ではこの交差項を含む分布が見える。

### 一モードの古典エネルギーと量子化

真空の周期境界内で、実の横モード $\vb f(\vb r)$ を $\int|\vb f|^2\dd V=1$ と規格化する。$\rot\rot\vb f=(\omega^2/c^2)\vb f$ とし、ベクトルポテンシャルを $q(t)\vb f/\sqrt{\eps_0}$ と置く。スカラーポテンシャルがゼロの横波では $\vb E=-\dot q\,\vb f/\sqrt{\eps_0}$、$\vb B=q\,\rot\vb f/\sqrt{\eps_0}$ である。周期境界で回転を部分積分すると $\int|\rot\vb f|^2\dd V=(\omega/c)^2$ となる。これらを電気・磁気エネルギーに代入すると、一モードの全エネルギーは $(\dot q^2+\omega^2 q^2)/2$ である。

正準運動量を $p=\dot q$ とし、量子化の手順は [調和振動子](/physics/quantum-mechanics/harmonic-oscillator/) に従う。交換関係 $[\hat q,\hat p]=i\hbar$ のもとで $\hat a=(\omega\hat q+i\hat p)/\sqrt{2\hbar\omega}$ と置けば、$[\hat a,\hat a^\dagger]=1$ であり、積を展開すると $\hat H=\hbar\omega(\hat a^\dagger\hat a+1/2)$ が得られる。従って $\hat a^\dagger$ を一回作用させるごとにエネルギーは $\hbar\omega$ 増える。ここで用いた実モードは定在モードである。運動量を指定するには、同じ周波数の正弦と余弦のモードを組み合わせて空間並進の固有モードを選ぶ。

:::example[一モードで数を変える]
真空波長 $600\,\mathrm{nm}$ のモードに光子を二個加えると、真空状態に対するエネルギーの増分は $2hc/\lambda_0\simeq6.62\times10^{-19}\,\mathrm J$ である。零点エネルギーは差を取ると消える。
:::

:::example[二つの出力で保存される確率]
対称な二経路の位相差を $\phi$ とし、出口の振幅を $(1+e^{i\phi})/2$ と $(1-e^{i\phi})/2$ と置く。確率は $\cos^2(\phi/2)$ と $\sin^2(\phi/2)$ で、和は 1 である。$\phi=\pi$ で暗くなる出口の確率は、もう一つの出口へ移る。
:::

### 光電検出と統計量
一定のパワーから求めた計数率は長い時間の平均値である。独立の低確率事象が多数ある理想化では、平均 $\bar N$ の Poisson 分布を得る。導出は、$M$ 区間に分けて各区間の検出確率を $\bar N/M$ とした二項分布で $M\to\infty$ とする。$N$ を固定した確率は $\bar N^N e^{-\bar N}/N!$ に近づく。この分布で $\langle N(N-1)\rangle=\bar N^2$ は和の添字を二つずらすと得られ、分散は $\bar N$ となる。実在の光では光源の揺らぎが加わり、また単一光子源のように独立事象を仮定できない場合もある。[コヒーレンス](/physics/optics/coherence/) では同時計数に現れる違いを扱う。

:::example[単色光の計数]
真空波長 $600\,\mathrm{nm}$、パワー $1\,\mathrm{pW}$、検出効率 0.50 なら $hc/\lambda_0\simeq3.31\times10^{-19}\,\mathrm J$ なので期待検出率は約 $1.51\times10^6\,\mathrm{s^{-1}}$。Poisson 近似で $1\,\mathrm{ms}$ の測定なら平均は約 1510、標準偏差は約 39 である。
:::

### 状態を一光子ずつ更新する
理想的な光子計数で完全に検出される一光子の入力を、損失のない二経路素子へ入れる。出力の状態が $a\ket{1,0}+b\ket{0,1}$、$|a|^2+|b|^2=1$ なら、二検出器の数 $N_1,N_2$ は一試行で $(1,0)$ または $(0,1)$ になる。各平均は $|a|^2,|b|^2$ であるが、同じ試行で $N_1N_2=0$。場を連続波の照度だけで記述すると平均パワーの比は再現できるものの、この同時計数の条件は別の仮定なしには出ない。独立の Poisson 計数を仮定する場合とも区別する。

:::example[対称な分割で試行を繰り返す]
一光子を毎回用意し、$|a|^2=|b|^2=1/2$ の分割器に 200 回入れる。第 1 検出器の計数は二項分布で平均 100、分散 50。第 2 の計数は第 1 の計数を 200 から引いた値なので、二つの計数の共分散は $-50$ になる。同じ平均を持つ二つの独立な Poisson 計数では共分散はゼロになる。
:::

### 光子数の統計と検出の損失

光子のエネルギーが一定でも、ある時間窓に現れる個数は光源で異なる。平均 $m$ の Poisson 分布では $p_j=e^{-m}m^j/j!$ が $j$ 個の確率である。生成関数 $G(s)=\sum_{j\ge0}p_js^j$ は指数の級数から $e^{m(s-1)}$。$G'(1)=m$、$G''(1)=m^2$ なので、分散は $G''(1)+G'(1)-[G'(1)]^2=m$ となる。従って標準偏差を平均で割った値は $m^{-1/2}$ である。これは平均強度が一定の理想的なコヒーレント光での計数モデルであり、光源のゆらぎや検出器の暗計数は別に加える。

各光子を独立に確率 $\eta$ で記録する検出器を考える。もとの個数を $j$ に固定すると検出数の生成関数は二項定理により $(1-\eta+\eta s)^j$ になる。従って検出後は $G(1-\eta+\eta s)$ で表される。Poisson 入射なら $e^{\eta m(s-1)}$ であり、平均と分散の両方が $\eta m$ になる。一方、常に $j$ 個の入射なら検出数は平均 $\eta j$、分散 $\eta(1-\eta)j$ である。平均だけでは両入射を区別できないため、光統計の測定では多数の窓で分散も調べる。この関係を [コヒーレンス](/physics/optics/coherence/) では二検出器の計数相関へ拡張する。

:::example[一個未満の窓で光を数える]
平均が小さな入射を短い窓で繰り返し測る。暗計数と検出器の不感時間を無視する。

入力は平均入射パワー $P\simeq 1\times 10^{-13}\,\mathrm W$、計数窓の長さ $\Delta t\simeq 1\times 10^{-6}\,\mathrm s$、真空波長 $\lambda_0\simeq 6\times 10^{-7}\,\mathrm m$、検出効率 $\eta=\frac{1}{2}$、光量子の換算定数 $hc\simeq 1.9864\times 10^{-25}\,\mathrm{J\,m}$ である。

パワーを光子数へ換えるため、まず一個のエネルギーを求める。

$$
\begin{aligned}
E_\gamma&=\frac{hc}{\lambda_0},\\
E_\gamma/(\mathrm J)&\simeq \frac{\left(1.986445857\times 10^{-25}\right)}{\left(6\times 10^{-7}\right)}\simeq 3.3107\times 10^{-19}.
\end{aligned}
$$

窓の中の平均エネルギーを一光子のエネルギーで割る。

$$
\begin{aligned}
m&=\frac{P \Delta t}{E_\gamma},\\
m&\simeq \frac{\left(1\times 10^{-13}\right) \cdot \left(1\times 10^{-6}\right)}{\left(3.310743095\times 10^{-19}\right)}\simeq 0.30205.
\end{aligned}
$$

独立な検出損失を入れて実際の検出数の尺度を得る。

$$
\begin{aligned}
\mu&=\eta m,\\
\mu&\simeq \left(\frac{1}{2}\right) \cdot \left(0.3020469941\right)\simeq 0.15102.
\end{aligned}
$$

Poisson 計数で検出数をゼロに固定する。

$$
\begin{aligned}
p_0&=e^{- \mu},\\
p_0&\simeq e^{- \left(0.151023497\right)}\simeq 0.85983.
\end{aligned}
$$

Poisson 分布の分散が平均に等しいことから、相対的な統計誤差を求める。

$$
\begin{aligned}
\sigma_N/\mu&=\frac{1}{\sqrt{\mu}},\\
\sigma_N/\mu&\simeq \frac{1}{\sqrt{\left(0.151023497\right)}}\simeq 2.5732.
\end{aligned}
$$

一回の窓では検出なしが多い。平均の推定には独立な窓を足し合わせる。
:::

:::example[窓を延ばして相対誤差を下げる]
前例と同じパワーを入れ、窓だけを $100$ 倍にする。

入力は平均入射パワー $P\simeq 1\times 10^{-13}\,\mathrm W$、計数窓の長さ $\Delta t\simeq 1\times 10^{-4}\,\mathrm s$、真空波長 $\lambda_0\simeq 6\times 10^{-7}\,\mathrm m$、検出効率 $\eta=\frac{1}{2}$、光量子の換算定数 $hc\simeq 1.9864\times 10^{-25}\,\mathrm{J\,m}$ である。

パワーを光子数へ換えるため、まず一個のエネルギーを求める。

$$
\begin{aligned}
E_\gamma&=\frac{hc}{\lambda_0},\\
E_\gamma/(\mathrm J)&\simeq \frac{\left(1.986445857\times 10^{-25}\right)}{\left(6\times 10^{-7}\right)}\simeq 3.3107\times 10^{-19}.
\end{aligned}
$$

窓の中の平均エネルギーを一光子のエネルギーで割る。

$$
\begin{aligned}
m&=\frac{P \Delta t}{E_\gamma},\\
m&\simeq \frac{\left(1\times 10^{-13}\right) \cdot \left(1\times 10^{-4}\right)}{\left(3.310743095\times 10^{-19}\right)}\simeq 30.205.
\end{aligned}
$$

独立な検出損失を入れて実際の検出数の尺度を得る。

$$
\begin{aligned}
\mu&=\eta m,\\
\mu&\simeq \left(\frac{1}{2}\right) \cdot \left(30.20469941\right)\simeq 15.102.
\end{aligned}
$$

Poisson 計数で検出数をゼロに固定する。

$$
\begin{aligned}
p_0&=e^{- \mu},\\
p_0&\simeq e^{- \left(15.1023497\right)}\simeq 2.7614\times 10^{-7}.
\end{aligned}
$$

Poisson 分布の分散が平均に等しいことから、相対的な統計誤差を求める。

$$
\begin{aligned}
\sigma_N/\mu&=\frac{1}{\sqrt{\mu}},\\
\sigma_N/\mu&\simeq \frac{1}{\sqrt{\left(15.1023497\right)}}\simeq 0.25732.
\end{aligned}
$$

平均は $100$ 倍になり、相対標準偏差は $1/10$ になる。測りたい時間変動より窓が長い場合は、時間分解能との交換になる。
:::

:::example[同じパワーで波長を変える]
検出器の効率が波長で変わる場合を考える。波長を前例の二倍にする一方で効率を半分にする。

入力は平均入射パワー $P\simeq 1\times 10^{-13}\,\mathrm W$、計数窓の長さ $\Delta t\simeq 1\times 10^{-4}\,\mathrm s$、真空波長 $\lambda_0\simeq 1.2\times 10^{-6}\,\mathrm m$、検出効率 $\eta=\frac{1}{4}$、光量子の換算定数 $hc\simeq 1.9864\times 10^{-25}\,\mathrm{J\,m}$ である。

パワーを光子数へ換えるため、まず一個のエネルギーを求める。

$$
\begin{aligned}
E_\gamma&=\frac{hc}{\lambda_0},\\
E_\gamma/(\mathrm J)&\simeq \frac{\left(1.986445857\times 10^{-25}\right)}{\left(1.2\times 10^{-6}\right)}\simeq 1.6554\times 10^{-19}.
\end{aligned}
$$

窓の中の平均エネルギーを一光子のエネルギーで割る。

$$
\begin{aligned}
m&=\frac{P \Delta t}{E_\gamma},\\
m&\simeq \frac{\left(1\times 10^{-13}\right) \cdot \left(1\times 10^{-4}\right)}{\left(1.655371547\times 10^{-19}\right)}\simeq 60.409.
\end{aligned}
$$

独立な検出損失を入れて実際の検出数の尺度を得る。

$$
\begin{aligned}
\mu&=\eta m,\\
\mu&\simeq \left(\frac{1}{4}\right) \cdot \left(60.40939882\right)\simeq 15.102.
\end{aligned}
$$

Poisson 計数で検出数をゼロに固定する。

$$
\begin{aligned}
p_0&=e^{- \mu},\\
p_0&\simeq e^{- \left(15.1023497\right)}\simeq 2.7614\times 10^{-7}.
\end{aligned}
$$

Poisson 分布の分散が平均に等しいことから、相対的な統計誤差を求める。

$$
\begin{aligned}
\sigma_N/\mu&=\frac{1}{\sqrt{\mu}},\\
\sigma_N/\mu&\simeq \frac{1}{\sqrt{\left(15.1023497\right)}}\simeq 0.25732.
\end{aligned}
$$

入射光子数は二倍になるが、平均検出数は前例に一致する。分光測定で強度を比べるには、この波長依存性を補正する必要がある。
:::

## 演習問題

::::exercise{level=basic}
$\eps=4\eps_0,\mu=\mu_0$ の媒質で $v$、$n$、$\eta/\eta_0$ を求めよ。

:::solution
:eqref[eq-optics-em-wave] から $v=c/\sqrt{4}=c/2$、:eqref[eq-optics-index] から $n=2$ となる。:eqref[eq-optics-plane] のインピーダンス比は $\sqrt{1/4}=1/2$ である。
:::
::::

::::exercise{level=basic}
真空で $\lambda_0=500\,\mathrm{nm}$ の平面波を $+z$ 方向へ進める。$z=125\,\mathrm{nm}$ と原点の複素振幅比を求めよ。

:::solution
進行因子は $e^{i2\pi z/\lambda_0}$ である。$z/\lambda_0=1/4$ より比は $e^{i\pi/2}=i$ である。
:::
::::

::::exercise{level=basic}
真空中で $\eta_0=377\,\Omega$ とする。平面波の電場の実効値が $377\,\mathrm{V/m}$ のとき、正面の照度を求めよ。

:::solution
実効値は二乗の時間平均の平方根である。従って照度は $377^2/377=377\,\mathrm{W/m^2}$ となる。
:::
::::

::::exercise{level=basic}
吸収のない媒質で外向きの球面放射を遠方で観測する。放射源から $2\,\mathrm m$ の位置で照度が $9\,\mathrm{W/m^2}$ である。同じ放射方向で距離 $6\,\mathrm m$ における照度を求めよ。

:::solution
固定した微小立体角の中で放射パワーは保存され、その断面積は $r^2$ に比例する。従って逆二乗則から $9(2/6)^2=1\,\mathrm{W/m^2}$ となる。
:::
::::

::::exercise{level=advanced}
進行平面波に対して電気エネルギー密度と磁気エネルギー密度が等しいことを、実場を使って示せ。

:::solution
:eqref[eq-optics-plane] の比例係数 $\eta^{-1}$ は実数なので、実場でも $\abs{\vb H}^2=\abs{\vb E}^2/\eta^2$ である。$\mu/\eta^2=\eps$ を代入すると、両密度は瞬時に $\eps\abs{\vb E}^2/2$ となる。
:::
::::

::::exercise{level=advanced}
吸収率 $1/4$、パワー反射率 $1/2$ の面がある。残りは真空へ直進するとき、垂直入射で $Fc/P$ を求めよ。

:::solution
透過率は $1-1/4-1/2=1/4$ である。:eqref[eq-optics-pressure] から $Fc/P=1+1/2-1/4=5/4$ である。これは反射光による寄与 $2R=1$ と吸収の寄与 $1/4$ の和でもある。
:::
::::

::::exercise{level=advanced}
$\vb E=\hat{\vb z}\cos(kz-\omega t)$ が $k>0$ のとき電荷のない真空で進行平面波を表せるか、理由を示せ。

:::solution
空間発散は $\partial_z E_z=-k\sin(kz-\omega t)$ である。これは恒等的にゼロでない。真空で $\rho=0$ なら $\eps_0\divg\vb E=0$ が必要なので、磁場の選び方で補うことはできない。
:::
::::

::::exercise{level=advanced}
同じパワーの $400\,\mathrm{nm}$ と $800\,\mathrm{nm}$ の真空光を、同じ量子効率の検出器に順に入れる。後者の期待計数率は前者の何倍か。

:::solution
光子一個のエネルギーは $hc/\lambda_0$ である。このため率は固定パワーでは波長に比例し、比は $800/400=2$ になる。
:::
::::

::::exercise{level=advanced}
理想の計数器で、入射の平均光子数が $3$ の Poisson 光と、各窓に必ず $3$ 個の光子が入る光を測る。検出効率を $1/3$ とし、各光子が独立に検出されるとする。各光の平均、分散、検出ゼロの確率を求めよ。

:::solution
入射が Poisson なら検出後の平均と分散はともに $(1/3)3=1$ で、ゼロの確率は $e^{-1}$。各窓に $3$ 個入る場合は、各光子の検出を成功確率 $1/3$ の独立な Bernoulli 試行として足す。平均は $3(1/3)=1$、分散は $3(1/3)(2/3)=2/3$ である。ゼロの確率は三つとも検出されない確率 $(2/3)^3=8/27$ となる。平均が等しくても独立な窓での個数分布は異なる。
:::
::::

## 確認テスト

::::quiz
$e^{-i\omega t}$ 規約で $+z$ に進む成分を選べ。ただし $k>0$ とする。

- [x] $e^{ikz}$
- [ ] $e^{-ikz}$
- [ ] $e^{-kz}$

:::explanation
一定位相 $kz-\omega t=\mathrm{const.}$ を微分すると、位相面は $\omega/k$ の速度で $+z$ に進む。二番目は逆向きで、最後の式は位相の空間変化を表さない。
:::
::::

::::quiz{answer="400" tol="1e-6"}
$\lambda_0=600\,\mathrm{nm}$、$n=1.5$ のとき、媒質中の波長を nm 単位で答えよ。

:::explanation
境界で周波数が共通なので、波長は位相速度に比例する。比が $1/n$ であるため $\lambda=600/1.5=400\,\mathrm{nm}$ である。
:::
::::

::::quiz
同振幅の平面波二つが逆向きに進む定在波で、時間平均 Poynting ベクトルはどうなるか。

- [ ] 電場の腹で最大となる。
- [x] どこでもゼロとなる。
- [ ] 磁場の腹で最大となる。

:::explanation
足し合わせた電場は $2E_a\cos kz$、磁場は $2iE_a\sin kz/\eta$ である。その共役積の実部が各位置でゼロなので、平均エネルギー流はない。
:::
::::

::::quiz{answer="2" tol="1e-6"}
完全鏡に正面から照度 $I$ の光を入れるとき、圧力を $I/c$ の倍数で答えよ。

:::explanation
入射光の運動量率が面積当たり $I/c$ で、反射後は $-I/c$ である。鏡への移行は両者の差なので圧力は $2I/c$ となる。
:::
::::

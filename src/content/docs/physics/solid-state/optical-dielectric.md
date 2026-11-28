---
title: 誘電性と光学的性質
description: 局所電場、因果的な誘電応答、光と電子・イオンの結合から物質の光学スペクトルを導く。
notation:
  - symbol: '\epsilon_r(\omega)'
    meaning: 複素比誘電率
  - symbol: '\tilde n=n+i\kappa'
    meaning: 複素屈折率
  - symbol: '\omega_p'
    meaning: 自由キャリアのプラズマ角振動数
  - symbol: '\omega_{\rm TO},\omega_{\rm LO}'
    meaning: 長波長の横および縦光学角振動数
  - symbol: 'E_B,a_X'
    meaning: 励起子の束縛エネルギーと Bohr 半径
---

## この章の目標

- 微視的な分極から巨視的な誘電率を求める。
- 因果律と散逸を複素屈折率に結び付ける。
- 電子プラズマと極性光学フォノンの集団モードを区別する。
- バンド間吸収端、励起子、強誘電体のソフト化を計算する。

SI 単位を用いる。時間因子は $e^{-i\omega t}$ とし、透磁率は特に断らなければ $\mu_0$ とする。分極と束縛電荷の定義は :ref[def-em-polarization] および式 :eqref[eq-em-bound-charge] を使う。

波長が単位胞より十分長い範囲では、一つの胞の中で電場がどの原子にどう作用するかと、試料全体を通る光がどう変わるかを順に考えるとよい。本章では応答が線形で、特に断らなければ空間分散を無視できる等方体を扱う。この時間規約では、受動媒質の相対誘電率の虚部は正周波数で非負である。外部から励起した媒質で光が増幅される場合は、虚部が負になり得る。

## 誘電関数と局所場（Clausius–Mossotti の式）

### 分極する粒子が感じる場

分極 $\vb P$ は単位体積の双極子モーメントで、巨視的には $\vb D=\epsilon_0\vb E+\vb P$ である。分子が感じる局所場は平均場と一般に異なる。等方、または立方対称な媒質で小さい球を原子を中心にくり抜いて考える。外側の分極が一様なら空洞表面の束縛面電荷は $-\vb P\cdot\hat{\vb r}$ である。中心に向く電場の向きの負符号と合わせ、$\int \hat r_i\hat r_j\dd\Omega=4\pi\delta_{ij}/3$ を使えば空洞由来の場は $\vb P/(3\epsilon_0)$ となる。近傍の離散双極子を球の内部に戻したとき、立方対称な同一双極子の和は零である。電場カーネルが跡零の二階テンソルで、立方対称なら恒等行列に比例するためである。従って :term[局所場]{en="local field" yomi="きょくしょば"} は $\vb E_{\rm loc}=\vb E+\vb P/(3\epsilon_0)$ となる。

密度 $N$、分極率 $\alpha$ の分子に $\vb P=N\alpha\vb E_{\rm loc}$ を代入し、$\vb P=\epsilon_0(\epsilon_r-1)\vb E$ と置けば
$$
 \frac{\epsilon_r-1}{\epsilon_r+2}=\frac{N\alpha}{3\epsilon_0}
 \label{eq-ss-clausius-mossotti}
$$
を得る。$\alpha$ の単位は $\mathrm{C\,m^2V^{-1}}$ である。この等式を強い異方性や近接相関がある系に適用する場合、離散場の項が消えるとは限らない。

:::example[希薄と高密度の比較]
$x=N\alpha/(3\epsilon_0)$ とすると $\epsilon_r=(1+2x)/(1-x)$ である。$x=0.1$ なら $\epsilon_r=4/3$ であり、局所場補正を入れない一次値 $1+3x=1.3$ より大きい。$x\to1$ で発散するが、線形な分子模型がその密度まで保たれるかは別の問題である。
:::

:::example[振動子の分極率]
電子の変位に $m\ddot x+m\gamma\dot x+m\omega_0^2x=-eE_{\rm loc}$ を仮定する。双極子は $-ex$ なので $\alpha(\omega)=e^2/[m(\omega_0^2-\omega^2-i\gamma\omega)]$ となる。静的極限は正、共鳴より十分上では実部が負となる。
:::

### 試料の境界と空洞の境界

外部電源や電極が作る場を $\vb E_{\rm a}$ とする。試料に一様な分極が生じると外表面には $\sigma_b=\vb P\cdot\hat{\vb n}$ が現れる。例えば分極が広い平板に垂直なら、二枚の面電荷が平板内に $-\vb P/\epsilon_0$ の場を作る。厚さに比べて横幅が大きい試料の中央に限れば、この場に $\vb E_{\rm a}$ を加えた量が巨視的な内部場 $\vb E$ である。球形試料の中心では外表面の場は $-\vb P/(3\epsilon_0)$ である。実際、表面の $\sigma_b=P\cos\theta$ から中心への場の $z$ 成分を積分すると $-P(4\pi\epsilon_0)^{-1}\int\cos^2\theta\dd\Omega=-P/(3\epsilon_0)$ となる。ここで極軸を分極の向きに取った。

一方、内部の球形空洞では、媒質から空洞へ向く法線が $-\hat{\vb r}$ である。表面電荷の符号が試料外表面の場合と逆なので、中心の場は $+\vb P/(3\epsilon_0)$ となる。これが局所場補正であり、試料外表面の電荷が作る反電場とは区別する。空洞の半径は格子間隔より大きく、巨視的な場の変化の尺度より小さく選ぶ。空洞内に戻す双極子の近接場の係数は、点双極子からの位置を $\vb R$ として $(3R_i R_j-R^2\delta_{ij})/(4\pi\epsilon_0R^5)$ である。各項の跡が零であり、立方対称な集合で和をとると軸の交換で三対角項が等しく、軸の回転で非対角項が零になる。従って近接双極子が同じ分極をもつなら和は零である。

:::example[同じ外場で形状を変える]
一様分極が実現する模型で $\epsilon_r=4$ とし、平板に垂直な外場を $\vb E_{\rm a}$ とする。$\vb E=\vb E_{\rm a}-\vb P/\epsilon_0$ と $\vb P=3\epsilon_0\vb E$ から $\vb E=\vb E_{\rm a}/4$ となる。球の中心では $\vb E=\vb E_{\rm a}-\vb P/(3\epsilon_0)$ なので $\vb E=\vb E_{\rm a}/2$ である。どちらも材料の比誘電率は同じだが、外場に対して測る分極の大きさは内部場が違うために変わる。
:::

### 周波数依存する局所場

波長の長い領域で局所場の式を各周波数に用いる。$x(\omega)=N\alpha(\omega)/(3\epsilon_0)$ と置いて式 :eqref[eq-ss-clausius-mossotti] を解くと $\epsilon_r=1+3x/(1-x)$ である。複素共役を掛けて分母を実にすれば

$$
\Im\epsilon_r=\frac{3\Im x}{|1-x|^2}
$$

を得る。従って受動的な分子の $\Im\alpha(\omega)$ が正なら、この模型でも巨視的応答の散逸は正である。同じ分子でも分母の大きさによって吸収は変わる。非磁性で透明な周波数なら $\epsilon_r=n^2$ であり、式 :eqref[eq-ss-clausius-mossotti] は屈折率から分極率を求める Lorentz–Lorenz の式になる。

単一振動子の場合、$P=N\alpha (E+P/(3\epsilon_0))$ の分母を払うと

$$
\epsilon_r(\omega)-1
=\frac{Ne^2/(\epsilon_0m)}
 {\omega_0^2-Ne^2/(3\epsilon_0m)-\omega^2-i\gamma\omega}
\label{eq-ss-lorentz-local-frequency}
$$

である。誘起双極子の場により復元係数が下がるため、分子単独の周波数より低い側に巨視的な共鳴がある。$\omega_0^2-Ne^2/(3\epsilon_0m)>0$ はこの線形模型の安定条件である。零まで近づくと静的応答が増すが、実際の相転移を判定するには原子間相互作用、結晶の対称性、非調和項も必要である。

:::example[透明領域で屈折率から分極率を出す]
$n=1.50$ なら $(n^2-1)/(n^2+2)=(2.25-1)/(2.25+2)=5/17$ である。従って $\alpha=(15/17)\epsilon_0/N$ となる。静的な比誘電率に遅いイオン変位が寄与する材料でも、光学周波数で同じイオンが追随するとは限らない。比較する際は $\alpha$ と $n$ の周波数をそろえる。
:::

### 量子状態からの分極率

一つの分子の非縮退な基底状態を $\ket{0}$、励起状態を $\ket{j}$ とし、双極子の一成分を $d$、遷移角振動数を $\omega_{j0}=(E_j-E_0)/\hbar>0$ とする。静電場による摂動は $-dE_{\rm loc}$ である。反転対称な基底状態では一次のエネルギー変化が零なので、二次摂動は $-E_{\rm loc}^2\sum_{j\ne0}|\mel{j}{d}{0}|^2/(E_j-E_0)$ である。このエネルギーを $E_{\rm loc}$ について微分して符号を反転すると、双極子の一次項が得られる。従って

$$
\alpha(0)=2\sum_{j\ne0}\frac{|\mel{j}{d}{0}|^2}{E_j-E_0}
$$

を得る。双極子が結ぶ励起が低いほど分極しやすい。時間依存の場に対しても、一次応答の核は摂動 $-dE_{\rm loc}$ による状態の変化から $(i/\hbar)\theta(t)\langle[d(t),d(0)]\rangle$ となる。実際、相互作用表示で状態の一次補正は $(i/\hbar)\int^t d(t')E_{\rm loc}(t')\ket{0}\dd t'$ であり、その共役と合わせると応答は交換子になる。完全系を挿入すると交換子は $\sum_j|\mel{j}{d}{0}|^2(e^{-i\omega_{j0}t}-e^{i\omega_{j0}t})$ である。従って $\Im z>0$ について指数積分を行えば

$$
\alpha(z)=\sum_{j\ne0}\frac{2\omega_{j0}|\mel{j}{d}{0}|^2}
 {\hbar(\omega_{j0}^2-z^2)}
\label{eq-ss-molecular-polarizability}
$$

となる。放射や他の粒子との衝突を無視した結果なので、実軸には $z=\omega+i0^+$ として近づく。$\omega=0$ では静的摂動の式に戻る。

:::example[二つの遷移の静的寄与]
双極子の絶対値が同じ $d_0$ の二遷移で、励起エネルギー比が $1:3$ とする。静的寄与の比は $3:1$ なので、低い遷移だけで全体の $3/4$ を占める。一方、式 :eqref[eq-ss-molecular-polarizability] の分子の重みは $\omega_{j0}$ に比例するため $1:3$ になる。静的分極の寄与と、周波数の二乗を分母にもつ振動子としての重みは異なる。
:::

### 振動子強度を拘束する交換関係

一電子の Hamiltonian を $H=\vb p^2/(2m)+V(\vb r)$、双極子の一成分を $d=-ex$ とする。局所ポテンシャルは $x$ と交換するので、$[H,d]=ie\hbar p_x/m$、$[d,[H,d]]=e^2\hbar^2/m$ である。他方、二重交換子を展開して基底状態で平均すると

$$
\langle0|[d,[H,d]]|0\rangle
=2\sum_j(E_j-E_0)|\mel{j}{d}{0}|^2
$$

である。中間に固有状態の完全系を入れ、端に作用する $H$ を $E_0$ で置き換えた。離散準位のほかに連続準位があればその積分も含める。無次元の振動子強度を $f_{j0}=2m\omega_{j0}|\mel{j}{d}{0}|^2/(e^2\hbar)$ とすると、従って $\sum_jf_{j0}=1$ である。そのため式 :eqref[eq-ss-molecular-polarizability] の高周波側の先頭は $-e^2/(mz^2)$ になる。例えば虚軸 $z=i\xi$ では $\xi^2/(\omega_{j0}^2+\xi^2)$ が $1$ へ増加するので、高周波極限と和を交換できる。多電子の場合、座標だけの電子間相互作用は全双極子と交換するため、同じ計算で全電子の寄与を数える。

:::example[調和振動子が強度を使い切る場合]
一電子の調和振動子では $x$ は生成と消滅の演算子の和であり、基底状態からは第一励起状態だけを結ぶ。$|\mel{1}{d}{0}|^2=e^2\hbar/(2m\omega_0)$ なので $f_{10}=1$ である。式 :eqref[eq-ss-molecular-polarizability] に代入すると $e^2/[m(\omega_0^2-z^2)]$ となり、無損失の古典振動子の分極率と一致する。
:::

## 光学定数と Kramers–Kronig の関係

### 反射、吸収、因果律

一様な局所媒質で Maxwell 方程式に平面波を代入すると $k^2=\epsilon_r\omega^2/c^2$ である。受動な半無限媒質では、奥に進む伝播波または奥へ向かって減衰する波を選び、$\tilde n^2=\epsilon_r$、$k=\tilde n\omega/c$ とする。伝播波の強度吸収係数は $\alpha_{\rm abs}=2\kappa\omega/c$ である。真空からの垂直入射で電場と磁場の接線を連続にすると振幅反射係数は $(1-\tilde n)/(1+\tilde n)$ であり、その絶対値二乗が反射率 $R$ である。例えば $\tilde n=2$、損失零では $R=1/9$ である。

線形応答の核 $\chi(t)$ を $t<0$ で零とすれば $\chi(z)=\int_0^\infty\chi(t)e^{izt}\dd t$ は $\Im z>0$ で解析的である。ここで核は実で、必要な収束と上半平面での十分な減衰を仮定する。高周波定数を引いた $\chi=\epsilon_r-\epsilon_\infty$ に Cauchy 積分公式を用い、実軸上の極を小さい上半円で避けると
$$
 \Re\chi(\omega)=\frac1\pi\mathop{\rm P}\int_{-\infty}^{\infty}
 \frac{\Im\chi(\omega')}{\omega'-\omega}\dd\omega'
 =\frac2\pi\mathop{\rm P}\int_0^\infty
 \frac{\omega'\Im\chi(\omega')}{\omega'^2-\omega^2}\dd\omega'.
$$
くぼみの円弧が $-i\pi\chi(\omega)$ を与え、実軸の主値積分が $i\pi\chi(\omega)$ になることを用いた。虚部の式は実部と虚部を交換し符号を反転した全実軸の式である。正周波数だけにする際は $\chi(-\omega)=\chi(\omega)^*$ を用いた。これが :term[Kramers–Kronig の関係]{en="Kramers–Kronig relations" yomi="くらまーすくろーにっひのかんけい"} である。式 :eqref[eq-em-kramers-kronig] と同じ因果性を物質の応答に適用している。直流導電成分が $\omega=0$ に極を持つ場合は、極の項を明示してからこの式を使う。

:::example[静的応答を吸収から読む]
絶縁体で $\Im\epsilon_r(\omega')\ge0$ とし、式の $\omega=0$ を取ると $\epsilon_r(0)-\epsilon_\infty=(2/\pi)\int_0^\infty\Im\epsilon_r(\omega')\dd\omega'/\omega'$ である。同じ面積の狭い吸収線なら中心周波数を半分にすると静的寄与はおよそ二倍になる。
:::

:::example[侵入長の定義]
$\kappa=0.50$、真空波長 $\lambda_0=500\,\mathrm{nm}$ では強度の $1/e$ 長さ $1/\alpha_{\rm abs}=\lambda_0/(4\pi\kappa)\simeq79.6\,\mathrm{nm}$ である。振幅の $1/e$ 長さはその二倍である。
:::

### 誘電関数から反射と吸収を計算する

$\epsilon_r=\epsilon'+i\epsilon''$、$\tilde n=n+i\kappa$ とする。二乗をとれば $n^2-\kappa^2=\epsilon'$、$2n\kappa=\epsilon''$ であり、$n^2+\kappa^2=|\epsilon_r|$ である。受動媒質に入る正周波数の波に対し $n,\kappa\ge0$ を選ぶので

$$
n=\sqrt{\frac{|\epsilon_r|+\epsilon'}{2}},\qquad
\kappa=\sqrt{\frac{|\epsilon_r|-\epsilon'}{2}} .
$$

無損失で $\epsilon'<0$ の場合は $n=0$、$\kappa>0$ となる。この場合にも振幅は内部で指数減衰するが、時間平均の熱吸収は零である。散逸による減衰と、伝播できない周波数でのエバネッセントな減衰は $n,\kappa$ を共に見て区別する。

実電場を $\Re[\vb E e^{-i\omega t}]$ とすると、分極電流の複素振幅は $\vb J_P=-i\omega\epsilon_0(\epsilon_r-1)\vb E$ である。従って時間平均の吸収電力密度は $\frac12\Re(\vb E^*\cdot\vb J_P)=\frac12\omega\epsilon_0\epsilon''|\vb E|^2$ である。受動媒質で $\omega>0$ なら $\epsilon''\ge0$ となる。非磁性の進行波の Poynting 流束は $I(z)=\frac12\epsilon_0 cn|\vb E(0)|^2e^{-2\kappa\omega z/c}$ である。$-\partial_z I=\omega\epsilon_0 n\kappa|\vb E(z)|^2$ なので、$2n\kappa=\epsilon''$ により両方のエネルギー収支が一致する。

真空から垂直に入る場合、境界上の電場と磁場から $1+r=t$、$1-r=\tilde n t$ である。従って $r=(1-\tilde n)/(1+\tilde n)$、$t=2/(1+\tilde n)$ になる。内部の $z=0$ での電力流束と入射流束の比には屈折率の実部 $n$ が掛かり、

$$
 R=\frac{(1-n)^2+\kappa^2}{(1+n)^2+\kappa^2},\qquad
 \frac{I(0^+)}{I_{\rm in}}=n|t|^2
 =\frac{4n}{(1+n)^2+\kappa^2}=1-R
 \label{eq-ss-optical-reflectance}
$$

となる。これは式 :eqref[eq-em-fresnel] の垂直入射の場合に対応する。半無限の損失媒質では奥へ入った電力は全て熱に変わる。有限の薄膜では裏面からも透過するため、この $1-R$ を試料全体の吸収率とみなせない。

:::example[振幅透過係数と電力の違い]
$\tilde n=3+i$ なら反射率は $5/17$、境界直下へ入る電力の比は $12/17$ である。振幅透過の絶対値二乗は $4/17$ なので、これをそのまま電力の比にすると因子 $3$ を落とす。媒体中の電場と磁場の振幅比が真空と違うために係数が必要である。
:::

### 因果的な復元と直流極

実数の時間核から $\chi(-\omega)=\conj{\chi(\omega)}$ が成り立つので、上の実部の式と対をなす正周波数での虚部の式は

$$
\Im\chi(\omega)
=-\frac{2\omega}{\pi}\mathop{\rm P}\int_0^\infty
\frac{\Re\chi(\omega')}{\omega'^2-\omega^2}\dd\omega'
\label{eq-ss-optical-kk-imag}
$$

である。全実軸の積分で、実部が偶関数であるために $\omega'>0$ と $\omega'<0$ の二つの分母の差が $2\omega/(\omega'^2-\omega^2)$ になる。正しい高周波背景を引くのは、Cauchy 積分の大半円を零にするためである。材料の共鳴域だけを測ったときに用いる $\epsilon_\infty$ は、さらに高エネルギーの吸収線の寄与を定数で近似した背景である。その背後の全ての電子の運動も含めて無限周波数まで追う場合は背景が変わり得る。

有限の直流伝導度があると $\chi$ に $i\sigma_{\rm dc}/(\epsilon_0\omega)$ が現れる。通常の実周波数積分でこのまま静的実部を出そうとすると、低周波側の積分が発散する。$\sigma(z)=\sigma_{\rm dc}+\order{z}$ で零点のまわりに正則な模型では、この項を $i\sigma_{\rm dc}/(\epsilon_0z)$ として分け、残余に分散関係を使う。因果的な境界値 $i/(\omega+i0^+)=i\,\operatorname{P}(1/\omega)+\pi\delta(\omega)$ には、虚部の主値と実部のデルタ関数の両方が入る。例えば Drude 模型の残余の静的実部は $\epsilon_\infty-\omega_p^2/\gamma^2$ であり、直流伝導がある系全体の有限な静的誘電率を意味する値ではない。

:::example[帯域の端にある吸収線]
狭い線の $\int\epsilon''\dd\omega=A$ が同じ二試料を考える。一方の中心が $\omega_0$、他方が $2\omega_0$ なら、静的増分は近似的に $(2A/\pi)(1/\omega_0)$ と $(2A/\pi)(1/(2\omega_0))$ である。線の重みを無視して可視光域だけから静的分極率を推定すると、測定域より低い周波数の線の寄与を落とすことがある。
:::

### 振動子の強度と積分則

一つの Lorentz 振動子を $\chi=S/(\omega_0^2-\omega^2-i\gamma\omega)$、$S>0$ とする。無損失線の極限 $\gamma\downarrow0$ で $\omega=\omega_0+\delta$ に近いとき、虚部は $[S/(2\omega_0)](\gamma/2)/[\delta^2+(\gamma/2)^2]$ に近づく。この Lorentz 関数は面積 $\pi$ のデルタ関数へ収束するので、正周波数の線面積は $\pi S/(2\omega_0)$ である。Kramers–Kronig の静的な式に入れれば $\chi(0)=S/\omega_0^2$ を得る。発熱が消えた極限でも、共鳴位置には分布として吸収の重みが残る。

因果性の積分を $z=i\xi$、$\xi>0$ で評価すると $\chi(i\xi)=(2/\pi)\int_0^\infty\omega'\Im\chi(\omega')/(\omega'^2+\xi^2)\dd\omega'$ となる。例えば減衰した Lorentz 振動子では $\xi^2\chi(i\xi)\to S$ である。積分側の係数 $\xi^2/(\omega'^2+\xi^2)$ が $1$ へ単調に増えるので、受動媒質で積分がこの有限値に収束する場合に

$$
\int_0^\infty \omega\Im\chi(\omega)\dd\omega=\frac{\pi}{2}S
\label{eq-ss-dielectric-spectral-weight}
$$

を得る。複数の独立振動子なら $S$ をその和に置き換えればよい。低周波応答が強い振動子でも、面積にさらに $\omega$ を掛けた積分は静的増分と違う重みで線を数える。

## 金属の光学応答とプラズマ振動

### Drude 誘電関数

電流は $\vb J=\partial_t\vb P=-i\omega\vb P$ と表せるので、背景誘電率を $\epsilon_\infty$ とすれば $\epsilon_r=\epsilon_\infty+i\sigma/(\epsilon_0\omega)$ である。式 :eqref[eq-ss-drude] で緩和率を $\gamma=1/\tau$ と置くと
$$
 \epsilon_r(\omega)=\epsilon_\infty-\frac{\omega_p^2}{\omega(\omega+i\gamma)},\qquad
 \omega_p^2=\frac{ne^2}{m\epsilon_0}
 \label{eq-ss-drude-permittivity}
$$
を得る。無損失で $\omega<\omega_p/\sqrt{\epsilon_\infty}$ なら $\epsilon_r<0$ で、横波の $k$ は虚数となる。界面から金属の奥に向かって場は減衰し、半無限媒質の反射率は $1$ となる。上側では伝播できるが、その周波数で背景が一定という近似は実際のバンド間吸収の影響を受ける。

外部自由電荷のない縦波は $\vb k\cdot\vb D=0$ を満たす。縦電場が非零なら $\epsilon_r(\omega)=0$ であり、無損失では $\omega_p/\sqrt{\epsilon_\infty}$ に振動する。この値は電子流体を微小に動かした模型でも得られる。変位 $\vb u$ に伴う縦の分極は $-ne\vb u$ なので $\vb D=0$ から $\vb E=ne\vb u/(\epsilon_0\epsilon_\infty)$、電子の力は $-ne^2\vb u/(\epsilon_0\epsilon_\infty)$ である。

:::example[プラズマ端の変化]
$m$ と背景が一定で密度を四倍にするとプラズマ端の角振動数は二倍になる。逆に $\epsilon_\infty=4$ の背景は同じ $\omega_p$ で端を半分にする。
:::

:::example[損失と低周波極限]
$\omega\ll\gamma$ では $\epsilon_r\simeq\epsilon_\infty+i\sigma_0/(\epsilon_0\omega)$ の主要虚部が現れる。良導体なら $\tilde n\simeq(1+i)\sqrt{\sigma_0/(2\epsilon_0\omega)}$ である。振幅減衰長は $c/(\kappa\omega)=\sqrt{2/(\mu_0\sigma_0\omega)}$ となる。
:::

### 緩和があるときのプラズマ端

式 :eqref[eq-ss-drude-permittivity] の分母の複素共役を掛けると

$$
\epsilon'=\epsilon_\infty-\frac{\omega_p^2}{\omega^2+\gamma^2},\qquad
\epsilon''=\frac{\omega_p^2\gamma}{\omega(\omega^2+\gamma^2)}
$$

となる。実部が零になる正周波数は、$\Omega=\omega_p/\sqrt{\epsilon_\infty}>\gamma$ のとき $\sqrt{\Omega^2-\gamma^2}$ である。この実部の零と、縦モードの複素角振動数は区別する。後者は $\epsilon_r(z)=0$ を解いた $z_\pm=-i\gamma/2\pm\sqrt{\Omega^2-\gamma^2/4}$ である。等価に、時間領域の運動方程式は $\ddot u+\gamma\dot u+\Omega^2u=0$ である。従って $\gamma<2\Omega$ では減衰振動し、$\gamma>2\Omega$ では二つの減衰率を持つ緩和となる。

外部電荷を与える場合、縦電場は $\epsilon_r$ の逆数を含む。実際、Fourier 成分 $\rho_{\rm ext}$ に対し $i\vb k\cdot\vb D=\rho_{\rm ext}$、$\vb E=-i\vb k\phi$ なので $\phi=\rho_{\rm ext}/(\epsilon_0\epsilon_r k^2)$ である。従って縦の損失を測るときは $\epsilon''$ とともに :term[損失関数]{en="energy loss function" yomi="そんしつかんすう"} $\Im[-1/\epsilon_r]$ が現れる。Drude 模型で逆数をとると

$$
\Im\left[-\frac1{\epsilon_r}\right]
=\frac1{\epsilon_\infty}
\frac{\Omega^2\gamma\omega}{(\omega^2-\Omega^2)^2+\gamma^2\omega^2}
\label{eq-ss-drude-loss}
$$

となり、$\gamma\ll\Omega$ では $\Omega$ 付近に狭いピークを持つ。これは自由電子の衝突による $\epsilon''$ の低周波での増大とは別の特徴である。

:::example[三種類の周波数を比較する]
$\Omega=10\gamma$ とすると、$\epsilon'=0$ の角振動数は $\sqrt{99}\gamma\simeq9.950\gamma$ である。複素縦モードの実部は $\sqrt{99.75}\gamma\simeq9.987\gamma$、減衰率は $\gamma/2$ である。無損失の周波数 $10\gamma$ に近いが、減衰を含む観測では零やピークの定義によって数値が変わる。
:::

### 縦波の波数依存と静的な遮蔽

$k$ 依存を扱うには、密度勾配も復元力を与えることを考える。一様な背景中の電子を縦方向へ変位 $\vb u$ だけ動かすと、連続の式の線形化より $\delta n=-n\grad\cdot\vb u$ である。圧力変化を $\delta P=(\partial P/\partial n)\delta n$ と近似する。ここで微分は仮定した流体の状態方程式に沿う。圧力勾配による一電子当たりの力は $-\grad\delta P/n$ である。$\beta^2=(\partial P/\partial n)/m$ とし、平面波に対して $m(\partial_t^2+\gamma\partial_t)\vb u=-e\vb E-\grad\delta P/n$ を書けば、縦方向圧力の力は $-m\beta^2k^2\vb u$ となる。$\vb P_e=-ne\vb u$ を使うと

$$
\epsilon_L(k,\omega)=\epsilon_\infty-
 \frac{\omega_p^2}{\omega^2+i\gamma\omega-\beta^2k^2}
\label{eq-ss-longitudinal-plasma}
$$

である。無損失の零は $\omega^2=\Omega^2+\beta^2k^2$ であり、$k=0$ で局所模型に戻る。逆に固定 $k>0$ で $\omega=0$ とすれば $\epsilon_L(k,0)=\epsilon_\infty[1+\Omega^2/(\beta^2k^2)]$ となる。したがって $\phi(k)$ の分母は $\epsilon_0\epsilon_\infty[k^2+\Omega^2/\beta^2]$ で、静電場には遮蔽の長さ $\beta/\Omega$ が現れる。流体模型では個々の電子正孔励起への減衰を扱っていない。衝突なしの量子的な $k$ 依存を求めるには [バンド計算の方法と電子間相互作用](/physics/solid-state/band-calculation/) の密度応答を用いる。

:::example[長波長側の分散]
$\beta k=0.20\Omega$ では無損失の縦モードは $\sqrt{1.04}\Omega\simeq1.020\Omega$ であり、$k=0$ から約 $2.0\%$ 上がる。同じ $k$ の静的比誘電率は $\epsilon_\infty(1+25)=26\epsilon_\infty$ になる。有限周波数の集団振動と静的遮蔽は同じ模型でも異なる極限である。
:::

### 表面に沿って進む電磁場

金属と誘電体を平面 $z=0$ で接触させ、界面に沿う波数 $K$ を持つ横磁気偏光を考える。無損失とし、誘電体 $z>0$ の比誘電率を $\epsilon_d>0$、金属 $z<0$ を $\epsilon_m(\omega)$ とする。磁場の $y$ 成分が界面から離れると減る条件は、減衰定数 $\kappa_j=\sqrt{K^2-\epsilon_j\omega^2/c^2}>0$ を取り、$H_y$ を誘電体で $e^{-\kappa_d z}$、金属で $e^{+\kappa_m z}$ に比例させることである。共通の因子 $e^{iKx-i\omega t}$ は省いた。

接線の $H_y$ は連続である。また $\rot\vb H=-i\omega\epsilon_0\epsilon_j\vb E$ から $E_x=(i\omega\epsilon_0\epsilon_j)^{-1}\partial_zH_y$ であり、接線の $E_x$ の連続条件は $\kappa_d/\epsilon_d+\kappa_m/\epsilon_m=0$ となる。よって両者の誘電率は逆符号を持つ。両辺を二乗し、$\kappa_j^2$ を代入して

$$
K^2=\frac{\omega^2}{c^2}
 \frac{\epsilon_d\epsilon_m(\omega)}{\epsilon_d+\epsilon_m(\omega)}
\label{eq-ss-surface-plasmon-dispersion}
$$

を得る。元の符号条件を合わせると、無損失の束縛モードには $\epsilon_m<-\epsilon_d$ が必要である。この :term[表面プラズモンポラリトン]{en="surface plasmon polariton" yomi="ひょうめんぷらずもんぽらりとん"} は、金属の電荷の振動と界面に束縛された光の混成である。式の $K$ は誘電体を進む同じ周波数の光の波数より大きい。平らな界面で平面波の接線波数を保存したまま入射するだけでは励起できないため、例えば周期構造により接線波数を補う方法がある。

:::example[表面とバルクの周波数]
無損失 Drude 金属に $\epsilon_\infty=1$ を用いる。式 :eqref[eq-ss-surface-plasmon-dispersion] で $K$ が大きくなる極限の周波数は $\epsilon_m+\epsilon_d=0$、すなわち $\omega_p/\sqrt{1+\epsilon_d}$ である。真空なら $\omega_p/\sqrt2$、$\epsilon_d=3$ の媒質に接すれば $\omega_p/2$ となる。同じ金属のバルクの縦モードはこの模型では $\omega_p$ に留まる。光の波長より短い尺度まで $K$ を大きくとると、媒質の空間分散を無視した近似も見直す必要がある。
:::

## フォノンポラリトンと Lyddane–Sachs–Teller の関係

### イオン変位と電場が連成する

極性のある一つの光学分枝について相対変位 $\vb u$、換算質量 $M$、有効電荷 $Z^*e$、密度 $N$ を取る。巨視的な電場を零に保ったときの復元角振動数を $\omega_{\rm TO}$ とし、減衰を省くと $M(\omega_{\rm TO}^2-\omega^2)\vb u=Z^*e\vb E$ である。$\vb P_{\rm ion}=NZ^*e\vb u$ を背景 $\epsilon_\infty$ に足すと
$$
 \epsilon_r(\omega)=\epsilon_\infty+
 \frac{N(Z^*e)^2}{\epsilon_0 M(\omega_{\rm TO}^2-\omega^2)}.
$$
縦モード条件 $\epsilon_r(\omega_{\rm LO})=0$ を解けば $\omega_{\rm LO}^2-\omega_{\rm TO}^2=N(Z^*e)^2/(\epsilon_0 M\epsilon_\infty)$ である。これを誘電関数に戻して
$$
 \epsilon_r(\omega)=\epsilon_\infty\frac{\omega_{\rm LO}^2-\omega^2}{\omega_{\rm TO}^2-\omega^2},
 \qquad \frac{\epsilon_r(0)}{\epsilon_\infty}
 =\frac{\omega_{\rm LO}^2}{\omega_{\rm TO}^2}
 \label{eq-ss-lst}
$$
を得る。一極一零の因数分解を $\omega=0$ で評価したものが :term[Lyddane–Sachs–Teller の関係]{en="Lyddane–Sachs–Teller relation" yomi="りだんざっくすてらーのかんけい"} である。複数の無損失極を持つスカラー応答でも、分母と分子の次数が同じなら比の積に一般化する。

横モードでは $c^2k^2=\omega^2\epsilon_r(\omega)$ を代入すると、$\omega^2$ に関する二次式になる。一つの根は $k=0$ で零、もう一つは $\omega_{\rm LO}^2$ である。$k$ が大きい側では下枝は $\omega_{\rm TO}$ に下から近づき、上枝は背景中の光の分散に近づく。この混成を :term[フォノンポラリトン]{en="phonon polariton" yomi="ふぉのんぽらりとん"} という。$\omega_{\rm TO}<\omega<\omega_{\rm LO}$ では実 $k$ がなく、半無限の無損失媒質は光を全反射する。この高反射の領域を :term[残留線帯]{en="Reststrahlen band" yomi="ざんりゅうせんたい"} という。

:::example[静的誘電率]
$\epsilon_\infty=4$ と $\omega_{\rm LO}/\omega_{\rm TO}=1.5$ を取ると、$\epsilon_r(0)=9$ になる。光学分枝の縦横差が電子背景に加わる静的分極を測っている。
:::

:::example[低波数分枝]
下枝の小さい $\omega$ では $k=\sqrt{\epsilon_r(0)}\,\omega/c$ であり、音速でなく光速を静的な屈折率で割った傾きになる。先の $\epsilon_r(0)=9$ ならその傾きは $c/3$ である。
:::

### 長距離 Coulomb 力による縦横差

光学モードの変位が縦方向であれば、一様な極限に近づいても $\rho_b=-i\vb k\cdot\vb P_{\rm ion}$ が生じる。その Coulomb ポテンシャルには $1/k^2$ があり、一階微分して得る電場は $\vb E_L=-\vb P_{{\rm ion},L}/(\epsilon_0\epsilon_\infty)$ となる。実際、外部電荷がなければ $\epsilon_0\epsilon_\infty\vb E_L+\vb P_{{\rm ion},L}=0$ である。$\vb E_L$ がイオンに加える力は $-(NZ^{*2}e^2/(\epsilon_0\epsilon_\infty))\vb u_L$ であり、バネを硬くする。横変位では $\vb k\cdot\vb P_{\rm ion}=0$ で、この静電的な復元力は生じない。従って $k\to0$ の極限でも、縦と横を区別して近づくと周波数差が残る。

有効電荷 $Z^*e$ は分極変化と変位を結ぶ量であり、イオンを孤立原子から作るときの形式電荷と違い得る。一般の胞では $\partial P_i/\partial u_{s,j}$ に体積を掛けた行列を Born 有効電荷とする。全ての原子を同じ変位だけ並進させても中性胞の双極子は変わらないため、原子についてこの行列を足すと零である。光学分枝が赤外活性となるには、固有ベクトルに沿ってこの電荷を掛けて足した有効双極子が非零でなければならない。

:::example[同位体置換と静的値]
一分枝の模型で復元係数 $K=M\omega_{\rm TO}^2$ と有効電荷が変わらず、換算質量だけを $2M$ とする。$\omega_{\rm TO}^2$ と $\omega_{\rm LO}^2-\omega_{\rm TO}^2$ はともに半分になるため、$\omega_{\rm LO}/\omega_{\rm TO}$ は同じである。従って式 :eqref[eq-ss-lst] の静的比には質量が入らない。現実の同位体置換で構造や非調和効果も変わる場合は、単に質量を変えた計算からずれ得る。
:::

### 複数の光学振動と積の公式

等方で無損失のモデルで、赤外活性を持つ相異なる横振動の二乗を $t_1<\cdots<t_s$、振動子の強度を $A_j>0$ とする。$x=\omega^2$ と置けば $\epsilon_r(x)=\epsilon_\infty+\sum_{j=1}^s A_j/(t_j-x)$ である。極を除いた各実区間で導関数は $\sum_j A_j/(t_j-x)^2>0$ である。$t_j$ のすぐ上で $-\infty$、次の $t_{j+1}$ のすぐ下で $+\infty$ なので、その間に一つだけ零 $\ell_j$ がある。最後の区間では右端の極限が $\epsilon_\infty>0$ なので最後の零 $\ell_s$ がある。最初の極より下では全項が正である。従って

$$
t_1<\ell_1<t_2<\ell_2<\cdots<t_s<\ell_s,\qquad
\frac{\epsilon_r(x)}{\epsilon_\infty}
=\prod_{j=1}^s\frac{\ell_j-x}{t_j-x}
\label{eq-ss-lst-multiple}
$$

である。後者は共通分母へまとめて分子の零点と最高次の係数を比較して得る。$x=0$ とすれば $\epsilon_r(0)/\epsilon_\infty=\prod_j\ell_j/t_j$ になる。非活性な分枝は $A_j=0$ なので、この積には含めない。異方体では測る偏光や波数方向で縦の条件も変わるため、このスカラー式を任意の方向のモードにそのまま適用できるとは限らない。

:::example[二つの赤外活性分枝]
周波数二乗をある単位で $t_1=1,t_2=4,\ell_1=2,\ell_2=8$ とし、背景を $\epsilon_\infty=3$ とする。静的比誘電率は $3(2/1)(8/4)=12$ である。部分分数に直すと $\epsilon_r=3+7/(1-x)+8/(4-x)$ であり、確かに両極の強度が正で、$x=0$ で $3+7+2=12$ となる。
:::

### ポラリトンの分散を解く

$L=\omega_{\rm LO}^2$、$T_0=\omega_{\rm TO}^2$、$C=c^2k^2/\epsilon_\infty$ とすると、横波の方程式は $x(L-x)=C(T_0-x)$ である。従って

$$
x_\pm=\frac{L+C\pm\sqrt{(L-C)^2+4C(L-T_0)}}2,\qquad \omega_\pm=\sqrt{x_\pm}.
\label{eq-ss-phonon-polariton-branches}
$$

$L>T_0>0$、$C>0$ なら、二次多項式 $x^2-(L+C)x+CT_0$ は $x=0$ で正、$T_0$ で負、$L$ でも負であり、$x\to\infty$ で正になる。従って $0<x_-<T_0<L<x_+$ である。$C\ll L$ では $x_-\simeq CT_0/L$ となり、式 :eqref[eq-ss-lst] の静的屈折率が下枝の速度を決める。$C\to\infty$ では積 $x_-x_+=CT_0$ と $x_+/C\to1$ から $x_-\to T_0$ である。共鳴から上の枝は背景の光に近くなり、下の枝はイオンの振動に近くなる。

:::example[背景の光が縦振動の位置を通る波数]
$L=4,T_0=1$ を同じ周波数二乗単位で与え、$C=4$ の波数を見る。式 :eqref[eq-ss-phonon-polariton-branches] では $x_\pm=4\pm2\sqrt3$ であり、下枝は $x<1$、上枝は $x>4$ にある。背景中の光だけなら $x=C$ を通るところで、混成のため二つの解に分かれる。
:::

### 自由キャリアを含む極性結晶の縦モード

ドープした極性結晶では同じ縦電場が電子とイオンを動かす。長波長で双方の減衰を無視し、イオンに一つの光学分枝、自由キャリアに Drude 模型を使う。$x=\omega^2$、$T_0=\omega_{\rm TO}^2$、$L=\omega_{\rm LO}^2$、$\Omega^2=\omega_p^2/\epsilon_\infty$ とする。このとき誘電関数は $\epsilon_r/\epsilon_\infty=(L-x)/(T_0-x)-\Omega^2/x$ である。縦モードではこれが零なので、分母を払うと $x^2-(L+\Omega^2)x+T_0\Omega^2=0$ を得る。$L>T_0>0$ と $\Omega>0$ なら、この多項式は $x=0$ で正、$x=T_0$ で $T_0(T_0-L)<0$、$x=L$ で $\Omega^2(T_0-L)<0$ となる。従って一方の根は $0<x_-<T_0$、他方は $x_+>L$ にある。

希薄な極限 $\Omega^2\ll L$ では、根の和と積から $x_+\to L$、$x_-\simeq T_0\Omega^2/L=\omega_p^2/\epsilon_r(0)|_{n=0}$ である。従って低い分枝で自由キャリアがゆっくり動くときは、イオンの静的分極も遮蔽に加わる。$\Omega^2\gg L$ では $x_+\simeq\Omega^2$、$x_-\to T_0$ であり、速い分枝を遮蔽するのは主に電子背景である。

:::example[電子とイオンが同じ周波数域にある場合]
$T_0=1$、$L=\Omega^2=4$ を同じ周波数二乗単位で与える。このとき零点は $x_\pm=4\pm2\sqrt3$ となる。電子だけの縦モードとイオンだけの縦モードはいずれも $x=4$ であるが、両者を合わせると二つの異なる周波数に分かれる。
:::

## バンド間遷移と吸収端

### 結合状態密度

光の波数が小さく、電気双極子行列要素がしきい値付近で有限とする。吸収は Fermi の黄金律により $[f_v(\vb k)-f_c(\vb k)]|\mel{c\vb k}{H'}{v\vb k}|^2\delta(E_c(\vb k)-E_v(\vb k)-\hbar\omega)$ の積分に比例する。吸収の率から誘導放出の率を引くと、占有因子は $f_v(1-f_c)-f_c(1-f_v)=f_v-f_c$ になる。三次元の直接ギャップで両端を放物線とすれば遷移エネルギーは $E_g+\hbar^2k^2/(2m_r)$、$m_r^{-1}=m_c^{-1}+m_v^{-1}$ である。式 :eqref[eq-ss-semiconductor-joint-dos] の球殻積分により、滑らかな周波数前因子を除き吸収端は $(\hbar\omega-E_g)^{1/2}$ となる。帯端で禁制で、最初の非零の行列要素が波数に一次なら、行列要素の二乗から $k^2$ が加わり、べき指数は $3/2$ になる。

間接ギャップで、必要なフォノンの振動数と行列要素を端近くで一定とみなすと、電子と正孔の余分な運動エネルギーを独立に積分する。余剰 $\Delta$ に対して
$$
 \int_0^\Delta\sqrt{E}\sqrt{\Delta-E}\dd E=\frac{\pi}{8}\Delta^2 .
$$
$E=\Delta\sin^2\theta$ の置換で右辺を得る。$\Delta=\hbar\omega-E_g\pm\hbar\omega_{\rm ph}$ に対して許容側の二乗の立ち上がりとなり、フォノンを吸収する正号には $n_{\rm ph}$、生成する負号には $n_{\rm ph}+1$ が掛かる。

:::example[吸収端の傾き]
余剰光エネルギーを四倍にすると許容直接遷移の吸収は約二倍、同じ条件の間接遷移では約十六倍になる。ここでは背景、占有、行列要素を一定とした。
:::

:::example[帯の充填による青方移動]
n 型で伝導帯の $|\vb k|<k_F$ が満たされ、価電子帯が満たされているとする。直接遷移は Pauli 禁制により $k\ge k_F$ に限られ、吸収しきい値は $E_g+\hbar^2 k_F^2/(2m_r)$ に上がる。これは相互作用によるギャップの変化とは別の占有効果である。
:::

### 電流の行列要素から吸収を求める

独立電子の局所周期ポテンシャルを仮定し、電荷 $-e$、裸の質量 $m$ の電子を考える。光の空間変化を原始胞の尺度で無視すると、ベクトルポテンシャル $\vb A(t)$ による最低次の摂動は $e\vb A(t)\cdot\vb p/m$ である。これは運動エネルギー $(\vb p+e\vb A)^2/(2m)$ の一次項である。非局所ポテンシャルを使う計算では一般に速度演算子もその Hamiltonian に合わせる必要がある。

$\vb E(t)=\Re[E_0\hat{\vb e}e^{-i\omega t}]$ とし、$\vb E=-\partial_t\vb A$ を用いると、吸収を駆動する $e^{-i\omega t}$ の摂動の振幅は $eE_0\hat{\vb e}\cdot\vb p/(2im\omega)$ である。空間的に一様なこの演算子は格子並進と交換するので、前後で同じ結晶波数の Bloch 状態だけを結ぶ。実際、並進を行列要素の右と左から作用させると、$\vb k$ と $\vb k'$ の位相が全ての格子ベクトルに対して等しくなければ、その行列要素は零である。光の波数を保つと波数差はその値に逆格子ベクトルを加えたものになる。

黄金律の率を光子のエネルギー $\hbar\omega$ と掛け、単位体積の吸収電力を $\omega\epsilon_0\epsilon''|E_0|^2/2$ と比較する。誘導放出を引き、スピンを状態和に含めると

$$
\epsilon''_{\hat{\vb e}}(\omega)
=\frac{\pi e^2}{\epsilon_0m^2\omega^2 V}
\sum_{v,c,\vb k}
(f_{v\vb k}-f_{c\vb k})
\left|\mel{c\vb k}{\hat{\vb e}\cdot\vb p}{v\vb k}\right|^2
\delta(E_{c\vb k}-E_{v\vb k}-\hbar\omega)
\label{eq-ss-interband-epsilon}
$$

となる。外場の振幅と時間因子の規約が合っていれば、黄金律での $1/2$ と散逸電力での $1/2$ が上の係数に含まれる。光がなくても状態が崩壊する自発放出は、この古典場による線形吸収式の誘導放出とは別である。強い外部注入により $f_c>f_v$ になる領域では、その遷移からの寄与が負になり増幅に変わる。

:::example[反転対称性による選択則]
反転不変な波数で、価電子帯と伝導帯の状態に反転固有値 $s_v,s_c=\pm1$ があるとする。反転は運動量を $-\vb p$ へ送るため、行列要素 $M$ は $M=-s_cs_v M$ を満たす。従って両状態の固有値が同じなら $M=0$ である。異なる場合は反転だけから零とは判定できず、偏光によって別の対称性が制限を与えることがある。例えば鏡映に関して両状態が偶なら、鏡面に垂直な運動量だけが奇となり、その偏光に対する遷移は消える。
:::

### 次元、偏光、吸収端の読み取り

直接放物線帯の遷移エネルギー $\hbar\omega=E_g+\hbar^2 k^2/(2m_r)$ について、一つの谷の三次元の結合状態密度は、内部縮退度を $g_s$ として

$$
g_{\rm joint}(\hbar\omega)
=\frac{g_s}{4\pi^2}\left(\frac{2m_r}{\hbar^2}\right)^{3/2}
\sqrt{\hbar\omega-E_g}
$$

である。$\hbar\omega>E_g$ について、$g_s(2\pi)^{-3}4\pi k^2\dd k$ をエネルギー幅 $\hbar^2k\dd k/m_r$ で割った。二次元の帯では $2\pi k\dd k$ に換えるため、単位面積当たりの結合状態密度が一定 $g_s m_r/(2\pi\hbar^2)$ となる。従って理想的な同じ種類の許容遷移でも、三次元では平方根で、量子井戸の一つの部分バンド対では階段状に立ち上がる。線幅と励起子はこの立ち上がりを変える。

弱い吸収では $2n\kappa=\epsilon''$ から $\alpha_{\rm abs}=\omega\epsilon''/(cn)$ である。式 :eqref[eq-ss-interband-epsilon] を使い $n$ と行列要素をほぼ一定と見なすと、$(\alpha_{\rm abs}\omega)^2$ は直接許容端で $\hbar\omega-E_g$ に比例する。その範囲の直線を吸収零へ外挿すれば $E_g$ を見積もれる。間接遷移の二乗の立ち上がりに同じ縦軸を使うと直線にならない。励起子の線や欠陥に起因する低エネルギー側の吸収がある場合も、全域を一つの直線で合わせる根拠はない。

:::example[有効質量で吸収端が動く場合]
電子の質量 $m_c=0.10m$、正孔の質量 $m_v=0.40m$ なら換算質量は $m_r=0.080m$ である。充填された電子の準位が伝導帯底から $40\,\mathrm{meV}$ までなら、同じ波数でのしきい値の移動は $(m_c/m_r)40=50\,\mathrm{meV}$ になる。このうち $10\,\mathrm{meV}$ は、その波数の価電子状態が価電子帯の頂上より低いことによる。
:::

### 間接遷移の時間順序

前後の谷の波数差を光だけで埋められない場合に、光と格子変位をそれぞれ一回使う。光の吸収の演算子を $V_\gamma$、エネルギー $s\hbar\omega_{\rm ph}$ を電子へ移す演算子を $V_{\rm ph}^{(s)}$ とする。$s=+1$ はフォノン吸収、$s=-1$ は放出である。中間電子状態を $\ket a$ とすると、二次の行列要素は

$$
\begin{aligned}
M^{(s)}_{fi}=\sum_a\bigg[
&\frac{\mel f{V_{\rm ph}^{(s)}}a\mel a{V_\gamma}i}
 {E_i+\hbar\omega-E_a+i0^+}\\
&+\frac{\mel f{V_\gamma}a\mel a{V_{\rm ph}^{(s)}}i}
 {E_i+s\hbar\omega_{\rm ph}-E_a+i0^+}\bigg].
\end{aligned}
$$

例えば光を先に吸収する項では、相互作用表示の早い時刻 $t_1$ に位相 $e^{i[(E_a-E_i)/\hbar-\omega]t_1}$ が付く。相互作用を遠い過去にゆっくり立ち上げ、$t_1$ を $-\infty$ から後の時刻 $t_2$ まで積分すると、最初の分母に比例する因子が出る。逆の順序では最初に受け取るエネルギーが $s\hbar\omega_{\rm ph}$ に変わる。残りの時間積分の長時間極限はエネルギー保存 $E_f-E_i=\hbar\omega+s\hbar\omega_{\rm ph}$ を課す。二項は同じ最終状態を作るため、率を計算する前に振幅を足す。共鳴する中間状態がある場合は寿命を無視した分母をそのまま用いると発散し、逐次過程との区別が必要になる。

:::example[低温で残る間接吸収]
$E_g=1.20\,\mathrm{eV}$、$\hbar\omega_{\rm ph}=0.040\,\mathrm{eV}$ とする。フォノンを吸収する端は $1.16\,\mathrm{eV}$、作る端は $1.24\,\mathrm{eV}$ である。光子が $1.22\,\mathrm{eV}$ では前者だけがエネルギー条件を満たすが、温度零では元のフォノンがいないので吸収が消える。光子が $1.26\,\mathrm{eV}$ なら放出側も開き、$n_{\rm ph}+1$ の $1$ により温度零でも寄与が残る。
:::

## 励起子

### 電子と正孔の相対運動

:term[励起子]{en="exciton" yomi="れいきし"} は、電子と正孔が束縛した中性励起である。三次元の放物線バンドで遮蔽を定数 $\epsilon$ とする。$\vb R=(m_e^*\vb r_e+m_h^*\vb r_h)/(m_e^*+m_h^*)$ と $\vb r=\vb r_e-\vb r_h$ に変数を替えると、勾配の交差項は相殺する。全質量 $M=m_e^*+m_h^*$ と換算質量 $m_r$ を使い、エネルギーは $E_g+\hbar^2K^2/(2M)$ に相対 Hamiltonian $-\hbar^2\nabla_r^2/(2m_r)-e^2/(4\pi\epsilon r)$ の固有値を足したものとなる。

相対運動は [量子力学](/physics/quantum-mechanics/) の Coulomb 問題と同じで、質量と電荷結合を置き換える。従って $E_\nu=E_g-E_B/\nu^2+\hbar^2K^2/(2M)$、$E_B=m_re^4/[2(4\pi\epsilon)^2\hbar^2]$ である。最低状態は半径 $a_X=4\pi\epsilon\hbar^2/(m_re^2)$ の指数関数を相対方程式に代入しても確かめられる。光は同じ位置で電子正孔対を作る局所行列要素を持つため、包絡近似の遷移は相対波動関数の $\varphi(0)$ に比例する。この近似では原点で非零の $s$ 状態が主要な線を与える。

:::example[Wannier 型の尺度]
$m_r=0.10m_e$、$\epsilon=10\epsilon_0$ とすると $a_X=100a_B\simeq5.29\,\mathrm{nm}$、$E_B\simeq13.6\,\mathrm{meV}$ である。第一線はギャップ下 $13.6\,\mathrm{meV}$、第二線は $3.4\,\mathrm{meV}$ にある。
:::

:::example[局在励起の限界]
同一サイトの励起エネルギーが $E_0$、隣接サイトへの励起移動が $-J$ の一次元 Frenkel 型模型では、Bloch 和から $E=E_0-2J\cos Ka$ である。励起の重心は動くが、相対運動を連続空間の水素型で近似する根拠はない。
:::

### 束縛状態の規格化と有効質量

電子と正孔の有効質量を $m_e^*,m_h^*$ とすると、$\grad_e=(m_e^*/M)\grad_{\vb R}+\grad_{\vb r}$、$\grad_h=(m_h^*/M)\grad_{\vb R}-\grad_{\vb r}$ である。これを二つの運動エネルギーに代入すれば、$\grad_{\vb R}\cdot\grad_{\vb r}$ の係数が相殺し、$\grad_{\vb R}^2$ と $\grad_{\vb r}^2$ の係数はそれぞれ $1/M$ と $1/m_r$ になる。従って重心の波数 $\vb K$ は良い量子数であり、各内部状態の分散に共通して $\hbar^2K^2/(2M)$ が加わる。

三次元で $\varphi_{1s}=A e^{-r/a_X}$ と置くと、$\nabla^2 e^{-r/a_X}=(a_X^{-2}-2/(a_X r))e^{-r/a_X}$ である。Coulomb 項と係数 $1/r$ が消える条件は $a_X=4\pi\epsilon\hbar^2/(m_re^2)$ であり、固有値は $-\hbar^2/(2m_ra_X^2)=-E_B$ である。規格化 $\int_0^\infty4\pi r^2|A|^2e^{-2r/a_X}\dd r=1$ から $A=1/\sqrt{\pi a_X^3}$ となる。強い空間変化をもつ格子の周期部分はここに含まれないため、$a_X$ が格子間隔より十分大きいことが Wannier 型包絡近似の目安になる。

水素型の正規化した $s$ 状態の動径波動関数は

$$
R_{\nu0}(r)=\left(\frac2{\nu a_X}\right)^{3/2}
\sqrt{\frac{(\nu-1)!}{2\nu\,\nu!}}\,
e^{-r/(\nu a_X)} L_{\nu-1}^{1}\left(\frac{2r}{\nu a_X}\right)
$$

である。ここでは式 :eqref[eq-qm-hydrogen-radial] の Coulomb 固有状態の形を、$a_X$ へ置き換えて用いた。Laguerre 多項式は原点で $L_{\nu-1}^{1}(0)=\nu$ なので、角度関数 $Y_{00}=1/\sqrt{4\pi}$ を掛けると

$$
|\varphi_{\nu s}(0)|^2=\frac{1}{\pi a_X^3\nu^3}
\label{eq-ss-exciton-envelope-origin}
$$

となる。角運動量が非零の正則な動径解は原点付近で $r^\ell$ に比例するため、その原点値は零である。

### 光が選ぶ重心波数と内部状態

相対波動関数を箱で $\varphi(\vb r)=V^{-1/2}\sum_{\vb k}A_{\vb k}e^{i\vb k\cdot\vb r}$ と展開し、$\sum_{\vb k}|A_{\vb k}|^2=1$ とする。光の摂動が作る直接遷移の電子と正孔の波数は、光波数を零とすれば互いに反対であり重心波数は零である。バンド間行列要素を包絡が占める狭い波数域で定数 $M_0$ とすれば、励起子への行列要素は $M_0\sum_{\vb k}A_{\vb k}^*=M_0\sqrt V\,\conj{\varphi(0)}$ である。この和が原点値を測る理由は、局所的な双極子が電子と正孔を同じ位置に作るためである。光波数を保てば対応する $\vb K$ にずれるが、重心の運動エネルギーは $\hbar^2K^2/(2M)$ で見積もれる。

この近似では、単位体積当たりの許容な束縛線の重みは式 :eqref[eq-ss-exciton-envelope-origin] に従い、緩やかに変わる光学係数を除けば $\nu^{-3}$ である。ただしバンド間の双極子が対称性で零なら、原点値が非零でもこの項による線は現れない。ギャップを超えた連続状態の吸収にも Coulomb 引力が作用し、自由な電子正孔の平方根のしきい値から強度が変わる。連続状態の寄与は、束縛状態とは別に散乱状態の原点値を使って数える必要がある。

:::example[エネルギー間隔と強度を比較する]
$E_B=16\,\mathrm{meV}$ の許容直接励起子なら、$1s$ と $2s$ の線の差は $16(1-1/4)=12\,\mathrm{meV}$ である。$2s$ の線強度は包絡近似では $1s$ の $1/8$ である。線の差が $E_B$ の $3/4$ であるのに対し、強度の比はエネルギーの差から直接は決まらない。
:::

### 連続準位にも残る Coulomb 引力

ギャップより上で相対運動の波数を $k$、$\eta=1/(ka_X)>0$ とする。引力の散乱波には原点付近に集まりやすい性質がある。[量子散乱の Coulomb 問題](/physics/advanced-quantum-mechanics/scattering/) の解で、同章が用いる符号付き Coulomb パラメータを $-\eta$ とすれば、入射を $z$ 方向とした散乱波は

$$
\varphi_{\vb k}^{(+)}(\vb r)
=V^{-1/2}e^{\pi\eta/2}\Gamma(1-i\eta)
 e^{ikz}{}_1F_1(i\eta;1;ik(r-z))
$$

である。この式の係数は、遠方の Coulomb 位相を含む入射波の絶対値を自由な入射波と同じ $V^{-1/2}$ にしたものである。合流型超幾何関数の原点値は級数の定数項 $1$ であるから、吸収の行列要素に効く原点確率は自由平面波の値に対して

$$
\mathcal C(k)=V|\varphi_{\vb k}^{(+)}(0)|^2
=e^{\pi\eta}|\Gamma(1-i\eta)|^2
=\frac{2\pi\eta}{1-e^{-2\pi\eta}}
\label{eq-ss-exciton-continuum}
$$

倍となる。最後に Gamma 関数の漸化式と反射公式を使った。具体的には $\Gamma(1+i\eta)=i\eta\Gamma(i\eta)$ と $\Gamma(i\eta)\Gamma(1-i\eta)=\pi/\sin(\pi i\eta)$ から積は $\pi\eta/\sinh(\pi\eta)$ になる。$k a_X\gg 1$ なら $\mathcal C\to1$ であり、$k\downarrow0$ では $\mathcal C\simeq2\pi/(ka_X)$ なので自由な結合状態密度の $k$ に比例する因子を打ち消す。

一つの光学チャネルについて共通の滑らかな光学係数を外すと、束縛線の重みは $\sum_\nu(\pi a_X^3\nu^3)^{-1}\delta(\hbar\omega-E_g+E_B/\nu^2)$、連続側の重みは $\mathcal C(k)m_r k/(2\pi^2\hbar^2)$ である。後者は相対運動の球殻を数えてから原点値の係数を掛けた。束縛線の系列と連続側を合わせると、直接許容端での電子正孔の相互作用の効果が見える。ここでは有効質量、定数遮蔽、希薄励起、無限の寿命を仮定した。現実の線幅を再現するには散乱や不均一性による広がりを加える。

:::example[束縛線が密になる極限]
大きな $\nu$ では隣接線の間隔は $E_B(\nu^{-2}-(\nu+1)^{-2})\simeq2E_B/\nu^3$ である。線の重みをこの間隔で割ると $1/(2\pi a_X^3E_B)=m_r/(\pi\hbar^2a_X)$ へ近づく。これは上の連続側の式で $k\downarrow0$ とした値と一致する。従って線が十分密で個々を分解できないときには、ギャップの両側で平均したスペクトルの重みが滑らかにつながる。
:::

### 二次元での束縛と局在極限

電子と正孔を同じ平面に閉じ込め、相互作用は三次元の Coulomb 法則の平面内での値 $-e^2/(4\pi\epsilon r)$ を使う理想模型を考える。二次元で $\nabla^2 e^{-\beta r}=(\beta^2-\beta/r)e^{-\beta r}$ なので、$1/r$ を相殺するには $\beta=2/a_X$ が必要である。平面での規格化から $|\varphi(0)|^2=2\beta^2/\pi=8/(\pi a_X^2)$ で、束縛エネルギーは $4E_B$ となる。

これが最低固有値であることは、任意の十分滑らかで減衰する規格化関数 $\psi$ について

$$
\frac{2m_r}{\hbar^2}\langle H_{\rm rel}\rangle+\beta^2
=\int_{\R^2}|\grad\psi+\beta\hat{\vb r}\psi|^2\dd^2r\ge0
$$

から分かる。右辺を展開し、交差項を $\divg\hat{\vb r}=1/r$ で部分積分した。原点を囲む小円からの境界項は半径に比例して消える。$\psi$ が上の指数関数なら右辺が零になる。有限厚みや周囲と異なる誘電率がある薄膜では相互作用がこの $1/r$ 形から離れ、とくに束縛の倍率をそのまま使えない。

:::example[局在励起のうち光と結び付く状態]
同じ向きの遷移双極子を持つ $N$ サイトの環で、サイト $j$ に励起を局在させた状態を $\ket j$ とする。波数 $K$ の状態は $\ket K=N^{-1/2}\sum_j e^{iKja}\ket j$ である。一様光に対して各サイトの振幅が同じ $d_0$ なら全振幅は $(d_0/\sqrt N)\sum_j e^{-iKja}$ である。公比が $1$ の $K=0$ では振幅は $\sqrt N d_0$、他の許容波数では幾何級数の和が消える。従って相対運動が水素型でない Frenkel 励起でも、並進位相により光への結合が選ばれる。
:::

## 強誘電体とソフトモード

### 分極を秩序変数として扱う

外部電場零で対称性に関連した二つ以上の分極状態をもち、電場で切り替えられる物質を :term[強誘電体]{en="ferroelectric" yomi="きょうゆうでんたい"} という。一軸の均一分極 $P$ に対し、反転対称な相からの Landau 展開を
$$
 f(P)=\tfrac12a(T-T_c)P^2+\tfrac14bP^4-EP,\quad a,b>0
$$
とする。停留条件は $E=a(T-T_c)P+bP^3$ で、$E=0$ で $T>T_c$ なら最小は $P=0$、$T<T_c$ なら $P=\pm\sqrt{a(T_c-T)/b}$ である。曲率は上側で $a(T-T_c)$、下側の極小で $2a(T_c-T)$ なので、小信号の $\partial P/\partial E$ はその逆数である。この温度依存を Curie–Weiss 型と呼ぶ。一次転移の場合には四次係数が負で六次を正にするなど熱力学が変わる。

変位型で分極に比例する座標を有効慣性 $M_P$ で運動させれば $M_P\ddot{\delta P}=-f''(P)\delta P$ である。従って均一な横光学モードの $\omega^2=f''(P)/M_P$ は連続転移に近づくと零に近づく。式 :eqref[eq-ss-lst] で $\omega_{\rm LO}$ が有限なら、横モードの軟化が静的誘電率の増大に対応する。秩序無秩序型では、複数の局所極性状態間の移動による緩和が重要になることがある。

:::example[転移温度を挟む応答]
$T=T_c+\Delta$ と $T=T_c-\Delta$ を比べると、分極率 $\partial P/\partial E$ は下側で上側の半分になる。下側では一つの局所最小の近傍の線形応答を取り、分域の移動や分極反転は含めない。
:::

:::example[安定性の確認]
$T<T_c$ の $P=0$ は停留解だが $f''<0$ である。微小擾乱が成長するので平衡状態とみなせない。計算上の虚数フォノンが構造の不安定性を示すのは、この曲率を動力学で読んだものである。
:::

### 短絡と開放で変わる分極の安定性

単一のソフトな分極成分を $P$、速い背景分極を $\epsilon_0(\epsilon_b-1)E$ とする。外部の自由電荷で変位電束 $D$ を固定した平板を考えると $D=\epsilon_0\epsilon_b E+P$ であり、電場エネルギーは $(D-P)^2/(2\epsilon_0\epsilon_b)$ である。電荷固定の自由エネルギーは、外場項を省いた物質の $f_0(P)$ にこれを加えたものである。これを微分すると $f_0'(P)=(D-P)/(\epsilon_0\epsilon_b)=E$ となり、物質の構成式は同じだが $E$ 自体が $P$ で変わる。電極を短絡して $E=0$ に保つ場合と比べると、電極を開放し $D=0$ にした場合の曲率には $1/(\epsilon_0\epsilon_b)$ が加わる。従って均一な法線分極が作る反電場は、その状態を不利にする。試料は分極の向きの異なる分域に分かれて反電場を減らす場合がある。

強誘電体の小信号の絶対感受率を $\partial P/\partial E$ とすれば、比誘電率のソフト成分はこれを $\epsilon_0$ で割った値である。従って短絡に近い電場制御の条件で $T>T_c$ の均一応答は $\epsilon_r(0)=\epsilon_b+1/[\epsilon_0 a(T-T_c)]$ となる。実験で誘電率の温度依存を比較する際は、分域の移動や電極との不完全な遮蔽の寄与も考える。

:::example[反電場で常誘電相の範囲が広がる模型]
$1/(\epsilon_0\epsilon_b a)=20\,\mathrm K$ とする。均一で電極を開放した $D=0$ の曲率は $a(T-T_c+20\,\mathrm K)$ であり、短絡での $T_c$ より $20\,\mathrm K$ 低いところで零になる。例えば $T=T_c-10\,\mathrm K$ では短絡模型の均一な $P=0$ は不安定だが、開放模型で均一に限れば曲率は正である。空間不均一な分域が許される場合の転移温度は、この均一模型からは決まらない。
:::

### 六次の自由エネルギーと一次転移

一軸で反転対称な自由エネルギーを $f_0(P)=a_2P^2/2+bP^4/4+c_6P^6/6$ とし、$b<0,c_6>0$ を仮定する。$a_2$ は温度とともに増すとする。$P\ne0$ の停留条件は $x=P^2$ に対して $a_2+bx+c_6 x^2=0$ である。この条件から $a_2$ を消して、分極状態と $P=0$ の自由エネルギー差を計算すると $-b x^2/4-c_6x^3/3$ になる。これが零となる正の根と、そのときの二次係数は

$$
x_*=-\frac{3b}{4c_6},\qquad a_{2,*}=\frac{3b^2}{16c_6}
\label{eq-ss-ferroelectric-first-order}
$$

である。$P=0$ の曲率は $a_{2,*}>0$、分極した解の曲率は停留条件を使って $2x(b+2c_6x)$ であり、$x=x_*$ で正である。従って安定な二つの種類の状態が等しい自由エネルギーを持つ。分極は零から有限の値へ跳ぶので一次転移である。加熱側では分極した局所最小が $a_2=b^2/(4c_6)$ まで存在し、冷却側では $P=0$ が $a_2=0$ まで局所最小となる。これらの点は各局所最小の消滅であり、平衡転移の点とは区別する。

:::example[曲率が零になる前に転移する]
自由エネルギー密度と分極を適切な単位で割って $b=-4,c_6=3$ とする。転移では $x_* =1$、$a_{2,*}=1$ であり、分極は $0$ から $\pm1$ のいずれかへ跳ぶ。$P=0$ 側の曲率は $1$、分極側の曲率は $2(1)(-4+6)=4$ であり、どちらの均一モードもこの点では軟化し切っていない。誘電率の強い温度依存と、横モード周波数が厳密に零になることは別に確かめる必要がある。
:::

### 連続転移の分域壁と電場反転

四次模型で $T<T_c$ とし $A=a(T_c-T)>0$、$P_s=\sqrt{A/b}$ と置く。壁の法線を $x$ 方向、分極を壁に平行な方向に取り、反電場を無視できる中性の $180^\circ$ 壁を考える。傾きのエネルギー密度を $g(\partial_xP)^2/2$、$g>0$ とする。バルクの値を引いたエネルギー密度は $g(P')^2/2+b(P^2-P_s^2)^2/4$ である。変分を取れば $gP''=bP(P^2-P_s^2)$ である。$P'$ を掛けて積分し、両端で $P'(\pm\infty)=0$、$P(\pm\infty)=\pm P_s$ とすれば $g(P')^2/2=b(P^2-P_s^2)^2/4$ になる。

正の傾きの解を $P(0)=0$ から積分すると

$$
P(x)=P_s\tanh(x/\lambda),\qquad
\lambda=\sqrt{2g/A},\qquad
\sigma_{\rm wall}=\frac{4gP_s^2}{3\lambda}
\label{eq-ss-ferroelectric-wall}
$$

である。最後の面積当たりのエネルギーは二つの密度項が等しいので $\int g(P')^2\dd x$ であり、$\int_{-\infty}^\infty\operatorname{sech}^4u\dd u=4/3$ を使った。この積分は $v=\tanh u$ で $\int_{-1}^1(1-v^2)\dd v$ に換わる。従って $T\uparrow T_c$ では壁が広がる一方で壁エネルギーは低下する。

一様な分極反転を強制するなら、停留式 $E=-AP+bP^3$ と曲率 $-A+3bP^2=0$ が同時に成り立つ点で準安定枝が消える。分極 $P>0$ の枝に逆向きの場を掛けると $P=\sqrt{A/(3b)}$ で消え、この場の大きさは $2A^{3/2}/(3\sqrt{3b})$ である。これは均一な四次模型の反転限界である。実際には核生成や分域壁の移動でそれより小さな場でも反転が起きる場合があり、測定した履歴曲線は欠陥や駆動速度にも左右される。

:::example[壁の太さとエネルギー]
同じ $a,b,g$ のもとで温度差 $T_c-T$ を四倍にすると、$P_s$ は二倍、$\lambda$ は半分、$\sigma_{\rm wall}$ は八倍になる。壁を形成する熱力学的な費用と、壁を既に持つ試料での反転の容易さは同じ量で測れない。
:::

### 緩和する双極子と振動する分極

分極が原子の連続的な変位を表す模型では、電場を零に保った曲率を $K_P=f_0''(P_0)$、減衰係数を $\Gamma_P$ とすると

$$
M_P\ddot{\delta P}+\Gamma_P\dot{\delta P}+K_P\delta P=\delta E
$$

である。従って $\delta P/\delta E=(K_P-M_P\omega^2-i\Gamma_P\omega)^{-1}$ であり、弱減衰の共鳴は $\sqrt{K_P/M_P}$ にある。連続転移では曲率が零に近づき、静的分極率の増大とソフトモードの周波数低下が同じ係数から出る。

一方、各場所に双極子 $+d$ と $-d$ の二つの配置があり、場に応じて配置間を移る模型を考える。密度を $N$、一つの双極子が感じる平均場を $E+\lambda_dP$ とする。二状態の Boltzmann 重みから局所平衡の分極は $P_{\rm loc}=Nd\tanh[d(E+\lambda_d P)/(k_BT)]$ である。小さい $E,P$ では $P_{\rm loc}\simeq C(T)(E+\lambda_dP)$、$C(T)=Nd^2/(k_BT)$ となる。$\lambda_d>0$ として $T_c=\lambda_d Nd^2/k_B$ を定めると、$T>T_c$ の静的応答は $P/E=C(T)/(1-T_c/T)$ であり、これも Curie–Weiss 型である。

有限時間では $\dot P=-(P-P_{\rm loc})/\tau_d$ と仮定する。すると動的応答は

$$
\frac{P}{E}=\frac{C(T)}{1-T_c/T-i\omega\tau_d}
\label{eq-ss-ferroelectric-relaxational}
$$

であり、共鳴の代わりに緩和の角振動数 $(1-T_c/T)/\tau_d$ が現れる。静的な温度依存だけでは二つの極限を区別しづらい。測定では虚部のスペクトルと実部の分散を比較する。実際の材料では振動と配置の移動の両方が寄与することもある。

:::example[静的な増大と緩和の遅れ]
$T=1.25T_c$ なら $1-T_c/T=0.20$ である。配置間の時間 $\tau_d$ を固定したとき、巨視的な緩和時間は $5\tau_d$ となる。$\omega\tau_d=0.20$ では応答の虚部と実部が等しく、どちらも静的応答の半分である。これは一階緩和の特徴であり、慣性を持つ振動子の共鳴と同じ形ではない。
:::

### ひずみを許したときの応答

一軸のひずみを $\eta$、弾性定数を $C>0$ とする。反転で $P$ は符号を変えるが $\eta$ は変わらないため、反転対称な自由エネルギーで許される最低次の結合は $-Q\eta P^2$ である。弾性エネルギー $C\eta^2/2$ を加え、応力が零なら $\eta$ について極小にした値は $\eta=QP^2/C$ である。これを代入すると分極の四次項に $-Q^2P^4/(2C)$ が加わる。従って四次係数は固定ひずみの $b$ から $b-2Q^2/C$ に変わる。負になれば安定な自由エネルギーのために六次などを残す必要があり、上の一次転移の模型が現れる。

この二乗のひずみ応答を :term[電歪]{en="electrostriction" yomi="でんわい"} という。自発分極 $P_s$ のまわりで機械的には自由なら、電場へのひずみの一次応答は

$$
\left.\pdv{\eta}{E}\right|_{E=0}
=\frac{2Q P_s}{C}\left.\pdv{P}{E}\right|_{E=0}
$$

である。従って同じ分域で小さい電場を掛ければひずみも一次で変わり、分域の向きを逆にするとその符号が反転する。高温の反転対称相で $P_s=0$ ならこの一次の効果は消えるが、$P^2$ による二次応答は残る。

:::example[固定された試料と自由な試料]
各係数をそろえた単位で $b=4$、$2Q^2/C=5$ とする。固定ひずみなら四次係数が正なので、他の不安定性を無視すると均一な四次模型の連続転移が可能である。一方、ひずみを緩和させると係数は $-1$ になり、正の六次項とともに一次転移が起こり得る。機械的な境界条件が分極の平衡を変える例である。
:::

## 演習問題

::::exercise{level=basic}
$N\alpha/(3\epsilon_0)=1/4$ で比誘電率を求めよ。
:::solution
式 :eqref[eq-ss-clausius-mossotti] を $\epsilon_r$ について解いて $(1+2/4)/(1-1/4)=2$ となる。
:::
::::

::::exercise{level=basic}
真空から損失のない屈折率 $n=3$ の媒質に垂直入射する。反射率を求めよ。
:::solution
振幅は $(1-n)/(1+n)=-1/2$ なので、強度比は $1/4$ である。
:::
::::

::::exercise{level=basic}
$\epsilon_\infty=3$ と $\omega_{\rm LO}/\omega_{\rm TO}=2$ の誘電体で、$\epsilon_r(0)$ を求めよ。
:::solution
式 :eqref[eq-ss-lst] の比は周波数比の二乗であり、$3\cdot2^2=12$ となる。
:::
::::

::::exercise{level=basic}
三次元の許容直接遷移で、しきい値からの余分なエネルギーを九倍にする。滑らかな前因子を一定とみなした吸収の倍率を求めよ。
:::solution
結合状態密度は余剰の平方根だから $\sqrt9=3$ 倍である。
:::
::::

::::exercise{level=advanced}
無損失 Drude 媒質に対し、分散を $k$ の式として解け。
:::solution
$c^2k^2=\omega^2\epsilon_\infty-\omega_p^2$ なので、正周波数は $\sqrt{(c^2k^2+\omega_p^2)/\epsilon_\infty}$ となる。$k=0$ が横波の下限で、$\vb D=0$ の縦モードの周波数とも一致する。
:::
::::

::::exercise{level=advanced}
間接吸収の積分 $\int_0^\Delta\sqrt{E(\Delta-E)}\dd E$ を幾何で求めよ。
:::solution
$u=E-\Delta/2$ で $\sqrt{(\Delta/2)^2-u^2}$ になる。積分は半径 $\Delta/2$ の半円面積だから $\pi\Delta^2/8$ である。
:::
::::

::::exercise{level=advanced}
励起子で $m_r$ を三倍、$\epsilon$ を二倍とした場合の $a_X,E_B$ の比を求めよ。
:::solution
$a_X\propto\epsilon/m_r$ なので $2/3$ 倍、$E_B\propto m_r/\epsilon^2$ なので $3/4$ 倍である。
:::
::::

::::exercise{level=advanced}
強誘電体の四次模型で $T<T_c$、$E=0$ とする。$P=0$ から平衡分極までのエネルギー密度の減少を求めよ。
:::solution
$P^2=a(T_c-T)/b$ を代入すると、二次項は $-a^2(T_c-T)^2/(2b)$、四次項は $a^2(T_c-T)^2/(4b)$ である。よって減少の大きさは $a^2(T_c-T)^2/(4b)$ である。
:::
::::

::::exercise{level=basic}
式 :eqref[eq-ss-optical-reflectance] で $n=\kappa=1$ の半無限媒質の反射率と、媒質へ入る電力の比を求めよ。
:::solution
反射率の分子は $(1-1)^2+1^2=1$、分母は $(1+1)^2+1^2=5$ なので $R=1/5$ である。境界に入る電力の比は $4n/5=4/5$ であり、両者の和は $1$ である。媒質中に損失があるので、十分深い場所ではこの入射分は全て減衰する。
:::
::::

::::exercise{level=advanced}
式 :eqref[eq-ss-lst-multiple] で $t_1=1,t_2=3,\ell_1=2,\ell_2=6$ とする。どのような $\epsilon_\infty>0$ に対しても、静的な比 $\epsilon_r(0)/\epsilon_\infty$ と二つの極の強度 $A_1,A_2$ を求めよ。
:::solution
原点で積を取ると比は $(2/1)(6/3)=4$ である。$A_1$ は $\epsilon_r(x)$ に $1-x$ を掛けて $x=1$ とした値で、$\epsilon_\infty(2-1)(6-1)/(3-1)=5\epsilon_\infty/2$ となる。同様に $A_2=\epsilon_\infty(2-3)(6-3)/(1-3)=3\epsilon_\infty/2$ である。いずれも正であり、$\epsilon_\infty+A_1/1+A_2/3=4\epsilon_\infty$ になる。
:::
::::

::::exercise{level=advanced}
光学チャネルが一つの三次元励起子について、独立電子正孔の場合のギャップ上の重み $m_r k/(2\pi^2\hbar^2)$ に式 :eqref[eq-ss-exciton-continuum] の係数を掛けよ。結果がしきい値で有限になることを示せ。
:::solution
$\eta=1/(ka_X)$ と置き、結合状態密度に $2\pi\eta/[1-e^{-2\pi\eta}]$ を掛ける。積は $m_r/[\pi\hbar^2a_X(1-e^{-2\pi/(ka_X)})]$ である。$k\downarrow0$ では指数の項は零に近づくので、有限で正の $m_r/(\pi\hbar^2a_X)$ に収束する。Coulomb 結合がないときの係数 $1$ とは、端の極限が異なる。
:::
::::

::::exercise{level=advanced}
式 :eqref[eq-ss-ferroelectric-wall] の壁の面積当たりエネルギーを、$A,b,g$ のみで表せ。温度を固定して $g$ を四倍にしたときの壁幅とエネルギーの比を求めよ。
:::solution
$P_s^2=A/b$ と $\lambda=\sqrt{2g/A}$ を代入すると $\sigma_{\rm wall}=(2\sqrt2/3b)\sqrt g\,A^{3/2}$ である。壁幅もエネルギーも $\sqrt g$ に比例するので、いずれも二倍になる。傾きの費用が増えると壁は幅を広げ、勾配を弱めて対応する。
:::
::::

## 確認テスト

::::quiz
一定背景を引いた因果的な誘電応答の実部を広い周波数域で復元する情報はどれか。

- [x] 虚部の周波数分布
- [ ] 単一周波数の反射強度のみ
- [ ] 物質の原始胞体積のみ

:::explanation
線形因果応答は上半平面で解析的であり、Kramers–Kronig の関係で各周波数の実部は虚部の主値積分から求まる。積分には周波数域全体の情報と背景の扱いが必要になる。
:::
::::

::::quiz{answer="2" tol="1e-6"}
背景と質量が一定の Drude 媒質で、電子密度が四倍のときプラズマ角振動数は何倍か。
:::explanation
式 :eqref[eq-ss-drude-permittivity] で $\omega_p^2$ は密度に比例するため、$\omega_p$ は二倍である。
:::
::::

::::quiz
無損失かつ一つの赤外活性モードをもつ媒質の Reststrahlen 帯に対応する周波数領域はどれか。

- [x] $\omega_{\rm TO}<\omega<\omega_{\rm LO}$
- [ ] $\omega>\omega_{\rm LO}$ 全体
- [ ] $0<\omega<\omega_{\rm TO}$ 全体

:::explanation
式 :eqref[eq-ss-lst] ではその区間で分子が正、分母が負である。したがって $\epsilon_r<0$ で横波が減衰し、半無限の無損失媒質では完全反射する。
:::
::::

::::quiz{answer="0.25" tol="1e-6"}
水素型励起子の第二準位の束縛エネルギーを、第一準位の値で割れ。
:::explanation
束縛は主量子数 $\nu$ に対し $E_B/\nu^2$ であるため、比は $1/2^2=1/4$ となる。
:::
::::

::::quiz
近赤外の範囲で損失が零で、金属模型の比誘電率の実部が負である。この周波数の横波について正しいものはどれか。

- [x] 半無限の金属で場は奥に向かって減衰し、時間平均の熱吸収は零である。
- [ ] 屈折率の実部が負で、普通の進行波が奥まで進む。
- [ ] 電場の振幅が減る分は必ず同じ割合で熱に変わる。

:::explanation
波数は純虚数となり、奥で発散しない指数減衰の解を選ぶ。比誘電率の虚部が零なので、吸収電力密度 $\omega\epsilon_0\epsilon''|\vb E|^2/2$ は零である。入口での電力流束も $\tilde n$ の実部に比例するため零であり、反射率は $1$ である。
:::
::::

::::quiz{answer="8" tol="1e-6"}
水素型の包絡近似で、$1s$ の束縛吸収線の重みは $2s$ の何倍か。共通の光学係数は同じとする。
:::explanation
線の重みは式 :eqref[eq-ss-exciton-envelope-origin] の原点確率に比例する。主量子数が $1$ と $2$ なので比は $2^3=8$ である。
:::
::::

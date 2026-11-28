---
title: "くりこみ群"
description: "尺度を変えた有効作用と相関関数を結び、固定点、演算子混合、重い粒子の効果を計算する。"
notation:
  - symbol: '\beta_g'
    meaning: 尺度に対する結合 g の微分
  - symbol: '\gamma_\phi'
    meaning: 場の異常次元
  - symbol: 'y_a'
    meaning: 赤外に向けた線形化の固有指数
---

## この章の目標

- 運動量殻の積分と尺度の再設定を区別する。
- 裸量の不変性から Callan–Symanzik 方程式を導く。
- 固定点の固有方向と有効理論の係数の発展を求める。

自然単位系 $\hbar=c=1$ を用いる。殻積分は Euclid の運動量で行い、散乱の不変量に戻す箇所では $\eta=\diag(+1,-1,-1,-1)$ を採る。長距離へ流す Wilson の変数と、微分を $\log\mu$ で取る方式では尺度の向きが逆である。

## Wilson のくりこみ群と有効作用

Euclid の紫外カットオフ $\Lambda$ を持つ場を低モード $\phi_< $ と $\Lambda/b<|p|<\Lambda$ の高モード $\phi_>$ に分ける。高モードを消した有限尺度作用を
$$
e^{-S_{\Lambda/b}[\phi_<]}=\int\mathcal D\phi_>\,e^{-S_\Lambda[\phi_<+\phi_>]}
$$
で定める。異なる二つの殻を積分する操作は順次積分と同じなので合成できる。この操作を尺度の再設定と合わせるものを :term[Wilson のくりこみ群]{en="Wilsonian renormalization group" yomi="うぃるそんのくりこみぐん"} という。高モードは積分されるため、一般に逆写像は不要である。

四次スカラー相互作用 $\lambda\phi^4/24$ では、高モードの連結累積量を使い $S_{\rm eff}=S_<+\langle S_I\rangle_>-\frac12\langle S_I^2\rangle_{>,c}+\cdots$。一頂点で高脚二本を縮約する六通りから $\delta m^2=(\lambda/2)\int_{\rm shell}\dd^4k/[(2\pi)^4(k^2+m^2)]$。二頂点の四点では三チャネルと二線の対称因子から、外部運動量が殻より十分小さいとき $\delta\lambda=-(3\lambda^2/2)\int_{\rm shell}\dd^4k/[(2\pi)^4(k^2+m^2)^2]$。$m\ll\Lambda$ なら四次元の球面積 $2\pi^2$ を用いて殻の後者は $\log b/(8\pi^2)$。符号が負なので殻を下ろすと結合は小さくなる。

### 殻を有限幅に取ったときの質量依存

四次元で $0<k_1<k_2$ の殻に対し
$$
I_1=\frac1{16\pi^2}\left[k_2^2-k_1^2-m^2\log\frac{k_2^2+m^2}{k_1^2+m^2}\right]
$$
である。角面積 $2\pi^2$ を掛け、$u=k^2$ として $u/(u+m^2)=1-m^2/(u+m^2)$ を積分した。$I_2=-\partial_{m^2}I_1$ も同じ殻に制限される。例えば $m=0$、$k_2=2k_1$ では $I_1=3k_1^2/(16\pi^2)$、$I_2=\log2/(8\pi^2)$。$m\gg k_2$ なら $I_2\simeq(k_2^4-k_1^4)/(32\pi^2m^4)$。Wilson の殻積分では、重い場が低い殻で伝播しにくいことを積分から直接読める。

### 流れを連続にして Wilson 作用の微分方程式を作る

有限格子で低モードの共分散を正定値行列 $C_s$ とし、積分済みの高モードの共分散を $A_s$ とする。$C_s+A_s=C_{\rm all}$、$s$ を大きくすると積分するモードを増やす。相互作用の重みを $F_s(\phi)=\int\dd\mu_{A_s}(\eta)\,e^{-S_{I,0}(\phi+\eta)}$、$S_{I,s}=-\log F_s$ と置く。微小な共分散 $\delta A=\dot A_s\delta s$ の独立 Gaussian を足すと、Taylor の二次が $\delta F=\delta A_{ij}\partial_i\partial_jF/2$、奇の項は零、高次は $\order{(\delta s)^2}$ となる。従って
$$
\partial_s S_{I,s}=
\frac12\dot A_{s,ij}
\left(\partial_i\partial_jS_{I,s}-\partial_iS_{I,s}\partial_jS_{I,s}\right).
\label{eq-qft-wilson-gaussian-flow}
$$
これがこの規約での Wilson 作用の連続な殻積分である。全作用には残りの共分散 $C_s$ による自由二次項を別に持たせる。座標の再スケーリングを併せると固定点を無次元の場と結合で定義できる。連続場では添字を運動量と成分へ拡張し、行列の縮約を運動量積分に替える。

:::example[四次相互作用の一次での流れ]
一自由度の $S_{I,s}=\lambda\phi^4/24+\delta m^2(s)\phi^2/2+\cdots$ を考える。$\lambda$ の一次では勾配の二乗を落とす。二階微分は $\lambda\phi^2/2+\delta m^2$ なので、$\partial_s\delta m^2=\lambda\dot A_s/2$。これは $\langle(\phi+\eta)^4\rangle=\phi^4+6A_s\phi^2+3A_s^2$ の二次の係数を直接読む方法と一致する。例えば $A_s$ を $0$ から $2$ に増やすと $\delta m^2=\lambda$。最後の $3A_s^2$ は真空項であり、質量には寄与しない。
:::

:::example[二次相互作用の厳密な流れ]
$S_{I,0}=b\phi^2/2$、一変数の $A_s=s\ge0$、$b>0$ とする。平方完成から $F_s=(1+bs)^{-1/2}\exp[-b\phi^2/(2(1+bs))]$。従って $b_s=b/(1+bs)$ は $\dot b_s=-b_s^2$ を満たし、定数は $\log(1+bs)/2$。:eqref[eq-qft-wilson-gaussian-flow] の $\partial^2S=b_s$ が定数項の微分、$(\partial S)^2=b_s^2\phi^2$ が逆符号の二次係数の微分を与える。例えば $b=2,s=1/2$ では $b_s=1$ である。
:::

### 殻を積分した後で紫外尺度を戻す

$b>1$ として殻を消した後、$p'=bp$、従って $x'=x/b$ と置く。自由運動項を $1/2$ に保つには $\phi(x)=b^{-(d-2)/2}\phi'(x')$。実際、$\dd^d x=b^d\dd^d x'$、$\partial_x=b^{-1}\partial_{x'}$ なので二つの指数が相殺する。この再設定により $\tfrac12 m^2\phi^2$ の係数は $b^2m^2$、$\lambda\phi^4/4!$ の係数は $b^{4-d}\lambda$。殻の量子補正はこの再設定の前の次元付き係数に加わる。

四次元で $b=e^{\delta \ell}$、$\delta\ell\ll1$ とし、カットオフは再設定の後も $\Lambda$ に戻す。$r=m^2/\Lambda^2$ とすれば殻積分 $I_1\simeq\Lambda^2\delta\ell/[8\pi^2(1+r)]$ なので
$$
\frac{\dd r}{\dd\ell}=2r+\frac{\lambda}{16\pi^2(1+r)}+\cdots,\qquad
\frac{\dd\lambda}{\dd\ell}=-\frac{3\lambda^2}{16\pi^2(1+r)^2}+\cdots .
$$
これは局所ポテンシャルの定数・二点・四点までを、一ループで射影した式である。四点の泡に外部運動量を入れずに係数だけを求めたので、急なカットオフをまたぐ運動量で同じ近似を使うことはできない。一高線の交換も二頂点から現れ得るが、一つの高線に流れるのは外部低運動量の和であり、一定の背景の四次係数には寄与しない。滑らかなカットオフを採れば有限な係数の形は変わる。

:::example[臨界面の座標は方式に依存する]
四次元の弱結合域で $\lambda$ の変化を一時的に高次とし、$r$ の式を一次まで使うと $r(\ell)=-\lambda/(32\pi^2)+[r(0)+\lambda/(32\pi^2)]e^{2\ell}$。初期値を $-\lambda/(32\pi^2)$ へ調整すれば、増大する質量方向の係数をこの近似の範囲で消せる。一方、質量なしのタドポールを次元正則化で計算すると尺度がなく零であり、臨界質量はその方式の座標で零に置ける。この二つの値が異なるのは尺度 $\Lambda$ に比例する局所質量の定義が違うためである。
:::

### 低モードの源に対する観測量を保つ

高低分割の Gaussian 測度が直積になり、全体の有限積分が収束するとする。低モードだけに源 $J_<\phi_<$ を付けるなら、先に高モードを積分しても分子と分母の比は同じになる。従って低モードの相関は全て $S_{\Lambda/b}$ から再現される。この等式は低運動量の有効作用を局所項で打ち切る前には厳密である。一方、高モードを含む演算子 $\mathcal O[\phi_<,\phi_>]$ も観測したいなら、低作用だけではその情報を取り出せない。その演算子を入れた重みを高モードで積分し、同じ低作用で割った関数を、低モードの挿入として用いる。

:::example[観測量の挿入も殻で変わる]
自由な二モードを、独立な低成分 $x$ と高成分 $y$ に分け、分散を $D_x,D_y$ とする。$\mathcal O=(x+y)^2$ に高モードの平均を取ると $\mathcal O_<=x^2+D_y$。従って零源の平均は $D_x+D_y$ であり、単に $y=0$ と書いた場合には $D_y$ を失う。一方、低成分の $x$ だけを観測する場合はもとの演算子をそのまま使える。殻を消す操作が作用だけの変更で済むかどうかは、どの観測量を測るかで決まる。
:::

### 一様な背景でポテンシャル全体の殻を計算する

Euclid の一成分で作用を $\int[(\partial\phi)^2/2+V(\phi)]$ に切り詰め、一様な低背景 $\varphi$ の上で薄い殻の揺らぎを二次まで残す。高モードは零運動量を持たないので、線形項 $V'(\varphi)\int\phi_>$ は零。二次核は $k^2+V''(\varphi)$。選んだ殻でこれが正なら、Gaussian の行列式から体積当たり
$$
\delta V(\varphi)=\frac12\int_{\rm shell}\frac{\dd^dk}{(2\pi)^d}
\log\left[1+\frac{V''(\varphi)}{k^2}\right]
$$
を得る。場と無関係な自由の行列式は引いた。例えば $V=m^2\varphi^2/2+\lambda\varphi^4/24$ では分子の $V''=m^2+\lambda\varphi^2/2$。対数を $\lambda\varphi^2$ で展開すると、二次は $\lambda\varphi^2/[4(k^2+m^2)]$、四次は $-\lambda^2\varphi^4/[16(k^2+m^2)^2]$。二次と四次の作用の規格化に戻せば前の $\delta m^2$ と $\delta\lambda$ である。六次は $+\lambda^3\varphi^6/[48(k^2+m^2)^3]$。最初に作用へ六次を置かなくても、殻の積分により対称性の許す項が生じる。

$\varphi=\Lambda^{(d-2)/2}\widehat\varphi$、$v(\widehat\varphi)=\Lambda^{-d}V(\varphi)$ とし、自由の運動項の規格化で尺度を戻すと
$$
\partial_\ell v=dv-\frac{d-2}{2}\widehat\varphi\,v'
+\frac{S_{d-1}}{2(2\pi)^d}\log(1+v'')
$$
の一ループの式を得る。$S_{d-1}$ は単位球面の面積である。対数の係数は殻の体積 $S_{d-1}\Lambda^d\delta\ell$ と揺らぎの $1/2$。この式は外部運動量に依存する項、場の規格化の流れ、高次のループを取り除いて一様ポテンシャルへ射影したものである。例えば $v''<-1$ の背景では高モードの二次核が不安定であり、対数を通常の実数としてこの式へ入れられない。低モードを残す作用の流れと、全てのモードを積分した Legendre の有効ポテンシャルは異なる定義を持つ。

### Gaussian の殻で混ざる次数を整理する

微小な高モードの共分散を $a$ とする。一点の多項式を高モードで平均すると、微分演算子 $e^{(a/2)\partial_\phi^2}$ が作用する。実際、Taylor 展開の偶数項の $\langle y^{2n}\rangle=(2n-1)!!a^n$ と $(2n)!$ の比は $(a/2)^n/n!$。従って $\phi^6$ は $\phi^6+15a\phi^4+45a^2\phi^2+15a^3$ に変わる。多項式をそのまま基底に取ると、たとえ自由固定点の近くの一次の流れでも、殻の自己縮約によってより低い次数へ混ざる。

低と高を合わせた一点の分散を $D$ とし、Wick 多項式の生成関数を $e^{t\phi-Dt^2/2}$ と定める。高モード $y$ の分散が $a$ なら、$\phi=x+y$ の高モード平均は $e^{tx-(D-a)t^2/2}$。従って $n$ 次の各 Wick 多項式は、同じ $n$ のまま低モードの分散に対応するものへ写る。例えば四次では $\phi^4-6D\phi^2+3D^2\mapsto x^4-6(D-a)x^2+3(D-a)^2$。このように自由固定点の固有方向を求めるときは、下位の演算子と混ざる部分をあらかじめ基底に含めておける。発散を持つ連続場では $D$ に正則化が必要だが、この有限の計算で何を引いているかが見える。

## Callan–Symanzik 方程式

裸場 $\phi_0=Z_\phi^{1/2}\phi$、裸相関 $G_0^{(n)}=Z_\phi^{n/2}G_R^{(n)}$ とする。固定した裸パラメータで $\mu\partial_\mu G_0^{(n)}=0$。$\beta_g=\mu\,\dd g/\dd\mu$、$\beta_{m^2}=\mu\,\dd m^2/\dd\mu$、$\gamma_\phi=\frac12\mu\,\dd\log Z_\phi/\dd\mu$ と置き積を微分すると
$$
\left(\mu\partial_\mu+\beta_g\partial_g+\beta_{m^2}\partial_{m^2}
+n\gamma_\phi\right)G_R^{(n)}=0
\label{eq-qft-callan-symanzik}
$$
を得る。これを :term[Callan–Symanzik 方程式]{en="Callan–Symanzik equation" yomi="からんしまんじくほうていしき"} と呼ぶ。$\dot g=\beta_g,\dot m^2=\beta_{m^2}$ の特性曲線に沿って $G_R(t)=e^{-n\int_0^t\gamma_\phi(t')\dd t'}G_R(0)$。固定点の質量ゼロ二点では次元解析と合わせて場の次元は $(d-2)/2+\gamma_\phi^*$ になる。

### 一つの物理尺度で対数を減らす

質量ゼロで単一尺度 $|p|>0$ を持つ振幅を考え、その正準次元を担う $|p|$ の冪を除いた関数を $F(p/\mu,g)$ とし、外脚の異常次元の係数までまとめて $\gamma_F$ と書くと $(\partial_{\log\mu}+\beta\partial_g+\gamma_F)F=0$ である。$t=\log(|p|/\mu)$ に沿って
$$
F(p/\mu,g)=\exp\left(\int_0^t\gamma_F(g(t'))\dd t'\right)F(1,g(t)).
$$
右辺の $F$ を基準点で計算すると、同じスケールの大きい対数が減る。例えば $\beta=bg^2$、$\gamma_F=ag$ とすれば指数因子は $(1-bgt)^{-a/b}$。展開の一次は $1+agt$ で、二次では $a(a+b)g^2t^2/2$ が出る。複数の独立な物理尺度が非常に離れていると、単一の $\mu$ の選択だけでは全ての対数を小さくできない。

### 境界値の選択と先頭対数の再帰

同じ物理振幅が、固定した $|p|$ と無次元結合 $g$ について $F=F(g,L)$、$L=\log(\mu/|p|)$ と表せるとする。異常次元がない場合、:eqref[eq-qft-callan-symanzik] は $\partial_L F=-\beta(g)\partial_g F$。$\beta=bg^2$、$F(g,0)=g$ なら $F=\sum_{n\ge0}a_n g^{n+1}L^n$ として $(n+1)a_{n+1}=-(n+1)b a_n$ を得る。従って $a_n=(-b)^n$。この先頭対数の和は $g/(1+bgL)$ となり、特性曲線で $\mu$ から $|p|$ に移した :eqref[eq-qft-log-running] と一致する。全ての $g^{n+1}L^n$ が重要でも、$g$ 自体が十分小さければ境界値での余分なループは順序を保って抑制できる。

:::example[固定次数と対数の再和の比較]
$bgL=1/2$ なら、二次の先頭対数までの部分和を $g$ で割った値は $1-1/2+1/4=3/4$。無限和は $2/3$ で、違いは $1/12$。同じ展開係数でも $L$ が大きく $bgL$ が一に近いと部分和は遅く収束する。対数をまとめた式はこの近似の $\beta$ だけを使うため、未知の高次の $\beta$ と境界値の有限な係数は別に評価する。
:::

### 特性方程式から二点関数の冪を読む

固定点で $\mu$ しか参照尺度がない質量ゼロの理論を考える。運動量保存のデルタを外したスカラー二点を $G(p)=|p|^{-2}F(|p|/\mu)$ と書く。:eqref[eq-qft-callan-symanzik] に $\beta=0$、$\gamma_\phi=\gamma_\phi^*$ を入れると、$z=|p|/\mu$ に対して $zF'(z)=2\gamma_\phi^*F(z)$。従って $G(p)\propto |p|^{-2+2\gamma_\phi^*}$。逆 Fourier 変換で $p$ を $b^{-1}p$ に変えると位置の二点は $|x|^{-(d-2+2\gamma_\phi^*)}$ の冪になる。ただし接触項は別に扱う。従って場の次元は $(d-2)/2+\gamma_\phi^*$。例えば $\gamma_\phi^*>0$ なら、同じ正準次元の自由場に比べて長距離の二点は速く減衰する。この結論は質量の関連方向を調整した固定点でのものである。

### 有限な場の再定義と異常次元

場を $\phi'=f(g)\phi$ と取り直す。裸場が同じなら $Z_\phi'=Z_\phi/f(g)^2$。従って
$$
\gamma_\phi'=\gamma_\phi-\beta_g\,\partial_g\log f(g).
$$
$f$ が固定点の近くで滑らかで零を持たなければ、$\beta(g_*)=0$ から $\gamma_\phi'(g_*')=\gamma_\phi(g_*)$。例えば $f=1+c g+\cdots$、$\beta=b g^2+\cdots$ の場合、二次の $\gamma_\phi$ の係数は $-bc$ だけ変わり得る。一方、固定点における二点の冪はこの有限な場の規格化に依存しない。同じ記号 $\gamma$ に対して数値を比較する前に、尺度微分の符号と場の再定義の向きを確認する。

### 既約頂点と相関関数で異常次元の符号を比べる

裸の連結相関の $n$ 本の外脚には $Z_\phi^{n/2}$ が掛かる。一方、作用の値を裸場と同じにする定義では $\Gamma_R[\Phi]=\Gamma_0[Z_\phi^{1/2}\Phi]$。従って $n$ 階の既約頂点は $\Gamma_R^{(n)}=Z_\phi^{n/2}\Gamma_0^{(n)}$ であり、裸量で尺度微分を取ると $(\mu\partial_\mu+\beta_g\partial_g+\beta_{m^2}\partial_{m^2}-n\gamma_\phi)\Gamma_R^{(n)}=0$。符号が連結相関と逆なのは、Legendre 変換で行列として逆関係にある二点の各外脚で確認できる。例えば自由でない二点では $\Gamma_R^{(2)}G_R^{(2)}=1$ なので、片方の微分に $2\gamma_\phi$、他方にその負が出て積の尺度微分が消える。

### 質量微分から挿入のある Callan–Symanzik 方程式へ移る

質量付き二点以上の連結相関で接触項を同じ規約に調整し、$G_R(p;m^2,g,\mu)$ とする。外のデルタを除いた工学次元を $d_G$ とすれば、一斉に物理運動量と $m,\mu$ を拡大する操作は
$$
\left(\sum_i p_i\cdot\partial_{p_i}+2m^2\partial_{m^2}
+\mu\partial_\mu-d_G\right)G_R=0
$$
と書ける。これは無次元比だけが残るという次元解析である。:eqref[eq-qft-callan-symanzik] から $\mu\partial_\mu$ を消すと
$$
\left(\sum_i p_i\cdot\partial_{p_i}-d_G-\beta_g\partial_g
-n\gamma_\phi\right)G_R
=-(2m^2-\beta_{m^2})\partial_{m^2}G_R.
$$
右辺は、質量項で作用を微分した挿入の相関として書き換えられる。例えば先頭の $\partial_{m^2}G_R^{(2)}=-1/(p^2+m^2)^2$ は、質量の局所二点を二本の自由線に挟んだもの。この式の質量項を小さいと見なせるのは、外部運動量を大きくしたときにその挿入が抑制される場合である。外部角度が特異な配置だと、尺度が大きくても一つの移行は小さくなり得る。

:::example[自由二点で質量微分の係数を確かめる]
$G=(p^2+m^2)^{-1}$ では $d_G=-2$、$\beta_g=\beta_{m^2}=\gamma_\phi=0$。従って $(p\cdot\partial_p+2)G=-2p^2/(p^2+m^2)^2+2/(p^2+m^2)=2m^2/(p^2+m^2)^2$。これは右辺 $-2m^2\partial_{m^2}G$ と一致する。例えば $p^2\gg m^2$ では尺度の破れは $G$ 自身に対して $O(m^2/p^2)$ である。臨界点を外れたまま $p\to0$ とすると逆に無視できない。
:::

### 先頭の対数に一つ高い係数を加える

$\beta(g)=b_1 g^2+b_2 g^3+\cdots$、$b_1\ne0$ とし、$t=\log(\mu/\mu_0)$、$g(0)=g_0$。$w=1-b_1 g_0t$ を定義し、先頭を $g^{(0)}=g_0/w$ とする。$b_2$ の寄与を一回だけ $\delta g$ に入れると
$$
\frac{\dd\delta g}{\dd t}-2b_1\frac{g_0}{w}\delta g
=b_2\frac{g_0^3}{w^3},\qquad \delta g(0)=0 .
$$
積分因子は $w^2$。従って $w^2\delta g=b_2g_0^3\int_0^t\dd s/(1-b_1g_0s)=-(b_2/b_1)g_0^2\log w$。これにより $b_1g_0t$ を固定して $g_0\to0$ とする数え方では、
$$
g(t)=\frac{g_0}{w}-\frac{b_2}{b_1}\frac{g_0^2\log w}{w^2}
+O(g_0^3)
$$
となる。$w$ が正で零から離れた領域を仮定した。この式を短い $t$ で展開すれば $g_0+b_1g_0^2t+(b_1^2t^2+b_2t)g_0^3+\cdots$。$b_1^2$ の二次対数だけなら先頭の解に既に含まれ、$b_2$ は同じ $g_0^3$ での一つの対数を供給する。境界の振幅自身に有限な $O(g_0^2)$ 項がある場合は、それにも先頭の流れを適用して同じ精度で合わせる。

## β 関数と異常次元

:term[ベータ関数]{en="beta function" yomi="べーたかんすう"} は尺度に対する係数の微分である。:eqref[eq-qft-phi4-counterterms] を用い $a=3/(32\pi^2)$、$\lambda_0=\mu^{2\epsilon}(\lambda+a\lambda^2/\bar\epsilon+\cdots)$ とする。固定裸量で微分すると $0=2\epsilon(\lambda+a\lambda^2/\bar\epsilon)+\beta_\lambda(1+2a\lambda/\bar\epsilon)+\cdots$。$\beta_\lambda=-2\epsilon\lambda+b\lambda^2+\cdots$ を代入すると $\epsilon\to0$ で $b=2a$。従って
$$
\beta_\lambda=-2\epsilon\lambda+\frac{3\lambda^2}{16\pi^2}+\cdots,\qquad
\gamma_\phi=\mathcal O(\lambda^2).
\label{eq-qft-phi4-beta}
$$
一ループの二点が運動量に依存しないことが後者の理由である。:term[異常次元]{en="anomalous dimension" yomi="いじょうじげん"} は場や演算子の次元の量子補正である。単一の無次元結合で $g'=g+cg^2+\cdots$ に変えると $\beta'=(\partial g'/\partial g)\beta$。例えば $\beta$ が $bg^2$ から始まるなら先頭係数 $b$ は変わらない。

### 結合の座標を変えて普遍な係数を確認する

$\beta(g)=b_1g^2+b_2g^3+\cdots$、$g'=g+ag^2$ とする。逆関数は $g=g'-ag'^2+\cdots$。$\beta'=(1+2ag)\beta$ を $g$ で展開すると三次係数は $b_2+2ab_1$ だが、$b_1g^2=b_1g'^2-2ab_1g'^3+\cdots$ により $2ab_1$ が相殺する。従ってこの単一結合の解析的な再定義では最初の二係数が同じである。例えば $b_1=2,b_2=-3,a=4$ では $g$ で書いた三次係数は $13$、$g'$ で書けば $-3$ に戻る。複数結合ではこの簡単な二係数の議論を各成分へそのまま適用できない。

### 多成分の四次結合で係数を数える

$N$ 成分のスカラー $\phi_i$ に $V=\lambda(\phi_i\phi_i)^2/4!$ を用いる。四階微分した木のテンソルは $\lambda T_{ijkl}$、ここで
$$
T_{ijkl}=\tfrac13(\delta_{ij}\delta_{kl}+\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk}).
$$
例えば $N=1$ なら $T=1$ である。泡の一チャネルの二頂点を縮約すると
$$
\sum_{a,b} T_{ijab}T_{abkl}
=\tfrac19[(N+4)\delta_{ij}\delta_{kl}+2\delta_{ik}\delta_{jl}
+2\delta_{il}\delta_{jk}].
$$
第一の $N$ は二つの頂点とも内部を $\delta_{ab}$ で縮約した項、$4$ は一方だけがそうした項の和である。残りの四項は二通りの外部の縮約を二つずつ作る。チャネルを全て足すと係数は $N+8$ になり、$\lambda^2/(32\pi^2\bar\epsilon)$ との積から
$$
\delta\lambda=\frac{(N+8)\lambda^2}{96\pi^2\bar\epsilon},\qquad
\beta_\lambda=-\varepsilon\lambda+\frac{(N+8)\lambda^2}{48\pi^2}+\cdots
\quad(d=4-\varepsilon).
$$
質量での二脚縮約は $\sum_a T_{ijaa}=(N+2)\delta_{ij}/3$ なので、$\beta_{m^2}=(N+2)\lambda m^2/(48\pi^2)+\cdots$。従って正の固定点 $\lambda_*=48\pi^2\varepsilon/(N+8)+O(\varepsilon^2)$ での質量方向の指数は $y_t=2-(N+2)\varepsilon/(N+8)+O(\varepsilon^2)$。$N=1$ では既出の結果と一致し、例えば $N=2$ では一次の補正係数は $2/5$ である。

### 有限変換で変わる係数を三次より先まで調べる

$\beta(g)=b_1g^2+b_2g^3+b_3g^4+\cdots$ とし、$g'=g+a g^2+c g^3+\cdots$ と置く。逆は $g=g'-ag'^2+(2a^2-c)g'^3+\cdots$。$\beta'(g')=(1+2ag+3c g^2)\beta(g)+\cdots$ に代入すると、$g'^4$ の係数は
$$
b_3'=b_3-a b_2+(c-a^2)b_1.
$$
例えば元の $\beta$ が $b_1g^2$ だけで、$c=0$ なら第四次の係数は $-a^2b_1$ になる。この係数が出ても二つの方式は同じ軌道を異なる結合で表している。これを $\beta=0$ の解の比較に使うなら、有限の再定義 $g_*'=h(g_*)$ も同じ次数で用いる。異なる次数で切り詰めた $\beta$ の零をそのまま比べると、近似誤差まで方式の不一致と解釈することになる。

### $O(N)$ の二ループの場の指数

質量ゼロで非零の Euclid 運動量をとり、サンセットに二つのテンソル $T_{iabc},T_{jabc}$ を入れる。縮約は
$$
\sum_{a,b,c} T_{iabc}T_{jabc}=\frac{N+2}{3}\delta_{ij}.
$$
実際、一頂点には外脚の $i$ を三つの内部添字の一つとつなぐ三項がある。相手の $j$ も同じ内部線を選ぶ三通りでは、残り二内部線の色が自由で係数 $N\delta_{ij}$。異なる線を選ぶ六通りでは、三線の色が外脚で固定され $\delta_{ij}$。頂点の $1/3$ を二つ掛けて $(3N+6)/9$。一成分の :eqref[eq-qft-two-loop-field-counterterm] をこの係数で置き換えれば
$$
\gamma_\phi=\frac{N+2}{36}\frac{\lambda^2}{(16\pi^2)^2}+O(\lambda^3).
$$
非零固定点 $\lambda_* /(16\pi^2)=3\varepsilon/(N+8)+O(\varepsilon^2)$ を入れると
$$
\eta=\frac{N+2}{2(N+8)^2}\varepsilon^2+O(\varepsilon^3).
$$
$N=1$ で $\eta=\varepsilon^2/54$、$N=2$ で $\varepsilon^2/50$。同じ固定点でも $\nu$ は $\varepsilon$ の一次で変化し、$\eta$ は二次から始まる。この差は質量タドポールには一ループで極がある一方、場の二微分の極は二ループから出ることと対応する。

### 二実場の四次テンソルから複数のベータ関数を作る

二つの独立な反転対称性を持つ $V=\lambda_1\phi^4/24+\lambda_2\chi^4/24+h\phi^2\chi^2/4$ を考える。四次の一ループの極は Hessian $H$ の $-\tr H^2/(64\pi^2\bar\epsilon)$。これは :eqref[eq-qft-scalar-one-loop-potential] の減算前の極を、Hessian の固有値について足すと得られる。質量を零にして四次だけを選べば、
$$
H_{11}=\frac{\lambda_1}2\phi^2+\frac h2\chi^2,\quad
H_{22}=\frac{\lambda_2}2\chi^2+\frac h2\phi^2,\quad
H_{12}=h\phi\chi .
$$
従って $\tr H^2$ の $\phi^4,\chi^4$ の係数はそれぞれ $(\lambda_1^2+h^2)/4,(\lambda_2^2+h^2)/4$、$\phi^2\chi^2$ は $h(\lambda_1+\lambda_2)/2+2h^2$。これらの局所極を引き、$\lambda_i,h$ の裸の係数を固定する。四次の対項は結合の二次式なので、$d=4-2\epsilon$ で $\beta=-2\epsilon g$ をその微分に入れると、一ループの有限係数は極の留数の二倍になる。よって $d=4-\varepsilon$ の式は
$$
\begin{aligned}
\beta_{\lambda_i}&=-\varepsilon\lambda_i+
 \frac{3(\lambda_i^2+h^2)}{16\pi^2},\\
\beta_h&=-\varepsilon h+\frac{h(\lambda_1+\lambda_2)+4h^2}{16\pi^2}.
\end{aligned}
$$
$h=0$ なら二場の流れは独立。$\lambda_1=\lambda_2=\lambda$、$h=\lambda/3$ なら $O(2)$ のポテンシャルになり $\beta_h=\beta_\lambda/3$。従って反転だけを仮定した式でも、回転対称な部分空間は流れの中で保たれる。

## ランニング結合定数と固定点

四次元で $\beta_g=b g^2$ を積分すると
$$
g(\mu)=\frac{g(\mu_0)}{1-bg(\mu_0)\log(\mu/\mu_0)} .
\label{eq-qft-log-running}
$$
$b>0$ なら高尺度で増え、分母の零は有限次数の摂動計算の破綻を示す。例えば $bg(\mu_0)\log(\mu/\mu_0)=1/2$ ならこの解で $g$ は二倍になる。

:term[固定点]{en="fixed point" yomi="こていてん"} は無次元結合の $\beta$ が全て零の点である。摂動が $\delta g$ なら一次の発展は $\dd \delta g_i/\dd\log\mu=B_{ij}\delta g_j$、$B=\partial\beta|_*$。$d=4-\varepsilon$ の四次スカラーでは :eqref[eq-qft-phi4-beta] の $2\epsilon$ を $\varepsilon$ と置き、$\lambda_*=16\pi^2\varepsilon/3$。その点で $\beta'=\varepsilon$ だから $\delta\lambda\propto\mu^\varepsilon$ であり、$\varepsilon>0$ では質量を臨界値に調整した軌道が赤外で近づく。

### $4-\varepsilon$ 次元の臨界指数

:eqref[eq-qft-two-loop-field-counterterm] を固定裸量で微分する。$\beta_\lambda=-2\epsilon\lambda+\cdots$ なので $\frac12\dd\log Z_\phi/\dd\log\mu=\lambda^2/[12(16\pi^2)^2]+\cdots$。$d=4-\varepsilon$ で $\lambda_*=16\pi^2\varepsilon/3+\cdots$ を入れると $\eta=2\gamma_\phi^*=\varepsilon^2/54+\mathcal O(\varepsilon^3)$。同じ固定点の無次元質量は $\beta_r=(-2+\varepsilon/3)r+\cdots$ なので、関連指数は $y_t=2-\varepsilon/3+\cdots$。相関長の指数は $\nu=1/y_t=1/2+\varepsilon/12+\cdots$ である。$\varepsilon=1$ に直接代入する値は低次数の外挿であり、展開の制御は $\varepsilon\ll1$ にある。

### 一ループで固定点間の流れを解く

$\beta_\lambda=-\varepsilon\lambda+b\lambda^2$ を $u=1/\lambda$ に直すと $\dd u/\dd t=\varepsilon u-b$、$t=\log(\mu/\mu_0)$。従って
$$
\lambda(\mu)^{-1}=\frac{b}{\varepsilon}+
\left(\lambda(\mu_0)^{-1}-\frac b\varepsilon\right)
\left(\frac\mu{\mu_0}\right)^\varepsilon .
$$
$\lambda(\mu_0)>0$ で有限なら、$\mu<\mu_0$ の範囲で固定点 $\lambda_*=\varepsilon/b$ に近づく。例えば $\lambda(\mu_0)=\lambda_*/2$、$(\mu/\mu_0)^\varepsilon=1/4$ なら $\lambda(\mu)=4\lambda_*/5$。$\varepsilon\to0$ の極限は、二項を別々に取ると発散するが和の指数を展開すると $1/\lambda(\mu_0)-b\log(\mu/\mu_0)$ に戻る。

### 軌道の線形化と座標の取り替え

結合ベクトルを $g$、正則で可逆な座標変換を $g'=h(g)$ とする。$\beta'=Dh\,\beta$ を固定点で微分すると、二階微分 $D^2h$ に掛かる $\beta$ が零なので
$$
B'=Dh(g_*)\,B\,Dh(g_*)^{-1}.
$$
従って固有値は有限の結合の組の取り方に依存しない。赤外へ $\mu$ を下げると、固有値 $\omega>0$ の線形なずれは $\mu^\omega$ として減衰する。一方 $\omega<0$ では増大する。固定点での長さの再設定 $b$ での指数 $y$ は $\mu\mapsto\mu/b$ から $-\omega$ である。縮退した固有値の行列が対角化できなければ、例えば $B=\omega 1+E_{12}$ で $e^{tB}=e^{\omega t}(1+tE_{12})$。このように冪に対数が掛かる場合があるので、固有値だけから全ての成分の比が定数になるとは言えない。

### 次元を作る積分定数と摂動式の適用範囲

例えば $\beta(g)=-b g^3$、$b>0$ の小結合理論を考える。$\dd(g^{-2})/\dd\log\mu=2b$ なので、固定した基準点のデータから
$$
\Lambda=\mu_0\exp[-1/(2bg(\mu_0)^2)],\qquad
g(\mu)^2=\frac1{2b\log(\mu/\Lambda)}
$$
と書ける。$\Lambda$ の尺度微分を取ると、明示の $\log\mu_0$ の微分と $\beta$ を通じた $g^{-2}$ の微分が相殺する。基準相互作用は無次元だが、走る軌道の境界条件が次元一の量を与える。$\mu$ が $\Lambda$ に近づくと流れの小結合の仮定が崩れるため、この解だけから物理的な極の存在を結論できない。例えば $2bg(\mu_0)^2=1/4$ なら $\Lambda/\mu_0=e^{-4}$。この指数的な尺度は、$g^2$ の有限次の Taylor 展開では得られない。具体的なゲージ理論で $b$ がどう決まるかは [非可換ゲージ理論の量子化](/physics/quantum-field-theory/gauge-quantization/) に続く。

### 一ループの QED で交換の有効電荷を比較する

質量を無視できる範囲の一 Dirac 種で、:eqref[eq-qft-wt] と電荷の関係から $\beta_e=e^3/(12\pi^2)$。$\alpha=e^2/(4\pi)$ なので $\dd\alpha/\dd\log\mu=2\alpha^2/(3\pi)$。これを積分すると
$$
\alpha(Q)=\frac{\alpha(Q_0)}
 {1-\frac{2\alpha(Q_0)}{3\pi}\log(Q/Q_0)} .
$$
例えば $Q,Q_0\gg m$ で :eqref[eq-qft-vacuum-polarization] の差を取ると、$\Pi_R(-Q^2)-\Pi_R(-Q_0^2)=2\alpha\log(Q/Q_0)/(3\pi)$。従って交換振幅を $\alpha[1+\Pi_R+\cdots]$ で比べた対数は、ランニングを一次だけ展開した値と同じである。全体の規格化を $Q_0$ で同じにしてから比べたので、物理電荷の減算に出た有限な $-5/3$ はここで消える。大きい対数を繰り返すときも基準点の有限な照合が要り、電荷の符号の規約や光子二点の符号を先に揃える。

:::example[どこまで一ループの走りを使うか]
$\alpha(Q_0)>0$ で分母の対数項が $1/2$ ならこの式の電荷は二倍。先頭の対数を一度だけ足すと $3/2$ 倍なので、対数の冪がまだ小さくない。この例の分母が零になる形式的な尺度は $Q_0e^{3\pi/(2\alpha(Q_0))}$。その近傍では $\alpha(Q)$ 自身が小さいという仮定が失われるため、実際の光子の伝播関数にこの式どおりの物理的な極があるとは判定できない。質量閾値を跨ぐ場合は係数自体も切り替わる。
:::

### 二結合の固定点で混合した方向を判定する

前の二場の模型で $\lambda_1=\lambda_2=\lambda$ とし、$\varepsilon>0$。$x=\lambda/(16\pi^2\varepsilon)$、$y=h/(16\pi^2\varepsilon)$ と置く。流れを $\varepsilon$ で割れば $-x+3(x^2+y^2)$ と $-y+2xy+4y^2$。$y=0$ では $x=0,1/3$。$y\ne0$ では $y=(1-2x)/4$ を第一式へ入れ、$60x^2-28x+3=0$ から $(x,y)=(3/10,1/10),(1/6,1/6)$ を得る。このうち $(3/10,1/10)$ は回転対称の $y=x/3$。安定性の行列は
$$
B=\varepsilon
\begin{pmatrix}-1+6x&6y\\2y&-1+2x+8y\end{pmatrix}.
$$
回転対称点では跡が $6\varepsilon/5$、行列式が $\varepsilon^2/5$ で、固有値は $\varepsilon,\varepsilon/5$。この二結合の面では両方が赤外で減衰する。対して $(1/3,0)$ では $\varepsilon,-\varepsilon/3$ なので、両場を結ぶ相互作用は赤外で増える。どちらの点でも質量の関連方向は別に調整する必要がある。四次の面での安定性だけで全パラメータが固定点へ近づくと判定してはならない。

:::example[座標を回すと同じ固定点を別の数で表せる]
$(x,y)=(1/6,1/6)$ では元の場に混合項がある。$\phi=(u+v)/\sqrt2,\chi=(u-v)/\sqrt2$ を代入すると、$\phi^4+\chi^4=(u^4+v^4)/2+3u^2v^2$、$\phi^2\chi^2=(u^4+v^4-2u^2v^2)/4$。$h=\lambda$ のとき混合の係数は消え、四次の自己結合は $2\lambda$。従って回転した二場の座標では $(1/3,0)$ の解に対応する。$B$ の固有値も $\varepsilon,-\varepsilon/3$ で同じである。固定点の座標だけを見て、独立な普遍性の類が増えたと解釈しない。
:::

## 合成演算子のくりこみと演算子混合

同一点の場の積は外脚の波動関数のくりこみだけで有限になるとは限らない。例えば $\phi^2$ の挿入の二点頂点に四次結合を一つ加えると二線泡が付き、:eqref[eq-qft-loop-master] の $I_2$ が現れる。源 $J_O$ を局所的に入れ、$J_O$ に一次で許される対項を取り出して挿入演算子を定める。外脚と対称性が同じ演算子同士は発散の多項式を共有でき、基底を選び $O_{0i}=Z_{ij}O_{Rj}$ と書く。すると
$$
\mu\frac{\dd O_R}{\dd\mu}=-\gamma O_R,\qquad
\gamma=Z^{-1}\mu\frac{\dd Z}{\dd\mu}.
\label{eq-qft-operator-mixing}
$$
係数を行ベクトル $C$ として物理量 $C O_R$ を保つと $\mu\,\dd C/\dd\mu=C\gamma$ である。例えば $\gamma=\begin{pmatrix}0&a\\0&0\end{pmatrix}$ が一定なら $O_{R2}$ は一定、$O_{R1}(\mu)=O_{R1}(\mu_0)-a\log(\mu/\mu_0)O_{R2}(\mu_0)$。この混合を無視すると $C O_R$ に尺度依存性を残す。

### 尺度ごとの混合行列が可換でない場合

:eqref[eq-qft-operator-mixing] で二つの連続区間にわたり $\gamma=\gamma_1,\gamma_2$ が一定なら、先に $\gamma_1$ が作用する演算子の発展は $O(t_1+t_2)=e^{-t_2\gamma_2}e^{-t_1\gamma_1}O(0)$。例えば $\gamma_1=E_{12}$、$\gamma_2=E_{21}$、$t_1=t_2=1$ とすると各行列は冪零二次なので積は $1-E_{12}-E_{21}+E_{22}$。順序を逆にすると最後が $E_{11}$ になる。従って一般の尺度依存 $\gamma$ に必要なのは順序付き指数であり、単に $\gamma$ を積分して指数を取る方法は非可換な混合を失う。係数は反対向きの逆行列で発展するので $CO$ は同じである。

### 質量源としての $\phi^2$ 挿入

Euclid 作用の $m^2$ を $m^2+h(x)$ に替えると、$\delta/\delta h$ は $\phi^2/2$ の挿入を作る。規格化相関では
$$
\partial_{m^2}\langle X\rangle
=-\int\dd^4x\,\bigl(\langle X\,\phi^2(x)/2\rangle-\langle X\rangle\langle\phi^2(x)/2\rangle\bigr)
+\langle\partial_{m^2}X\rangle .
$$
従って質量の対項を同時に微分すれば、積分された挿入の局所対項が得られる。裸の質量が $m_0^2=Z_{m^2}m^2$ で、$Z_{m^2}$ が $m$ に依存しない方式では、$O_R=\partial S_0/\partial m^2$ は $Z_{m^2} O_0$ であり、$O_0=\phi_0^2/2$。これは外脚との接触項や真空項の規格化を別に整えた挿入を指す。前節の規約なら $O_0=Z_{m^2}^{-1}O_R$ だから、$\gamma_{\phi^2}=-\mu\dd\log Z_{m^2}/\dd\mu=\beta_{m^2}/m^2$。Wilson–Fisher 固定点の $\phi^2$ の次元は $d-2+\gamma_{\phi^2}^*$、従って $y_t=2-\gamma_{\phi^2}^*$ になる。

例えば自由格子で $O_0=\phi_a^2/2$ を一点に置くと真空期待値は $D_{aa}/2$。$O_R=O_0-D_{aa}1/2$ と規格化すれば、この定数は消えるが、離れた二点との連結挿入は変わらない。曲がった背景や異なる質量の場がある場合は、対称性が同じ低次元の演算子を質量や曲率の冪で補った項も混合の候補になる。連結相関で定数が消えるからといって、局所演算子の全ての規格化から定数を常に省けるわけではない。

### 演算子基底を取り替えたときの混合

可逆な有限行列 $R(g)$ で $O'_R=R(g)O_R$ とする。:eqref[eq-qft-operator-mixing] を微分し、右に $R^{-1}$ を掛けると
$$
\gamma'=R\gamma R^{-1}-(\mu\dd R/\dd\mu)R^{-1}.
$$
固定点で行列が滑らかであり、陽に尺度を含まないなら第二項は $\beta$ が零なので消える。従って固定点での混合の固有値は基底に依らない。逆に流れの途中では $\gamma$ を各尺度ごとに対角化しても、その基底が回る分の微分が残る。例えば $\gamma=0$ で $R=\begin{pmatrix}1&a g\\0&1\end{pmatrix}$ とすると $\gamma'=-a\beta_g E_{12}$。この非対角成分は基底の取り直しから生じる。対応する係数は $C'=CR^{-1}$ なので、$\gamma'$ を含めて流しても積 $C'O'_R$ は同じである。

### 二種類の質量挿入で混合を計算する

二実場 $\phi,\chi$ の相互作用を $\lambda_\phi\phi^4/4!+\lambda_\chi\chi^4/4!+h\phi^2\chi^2/4$ とする。$\phi$ の二点に付く $\chi$ のタドポールは $h I_1(m_\chi^2)/2$。これは頂点の二つの $\phi$ を外脚に結ぶ $2!$ と、残る $\chi$ を縮約する一通りから出る。同じ方法で $\chi$ の二点には $h I_1(m_\phi^2)/2$。$I_1$ の極と裸質量の微分から、質量二乗の列 $r=(m_\phi^2,m_\chi^2)^\mathsf T$ は一ループで
$$
\frac{\dd r}{\dd\log\mu}
=\frac1{16\pi^2}
\begin{pmatrix}\lambda_\phi &h\\ h&\lambda_\chi\end{pmatrix}r
$$
を満たす。局所挿入を質量微分で $O=(\phi^2/2,\chi^2/2)^\mathsf T$ と定義する。裸質量 $r_0=Z_m r$ から、連鎖律で $O_R=Z_m^\mathsf T O_0$。従って異常次元の行列は質量の発展行列の転置である。この例の行列は対称なので両者は同じだが、一般の基底で転置を落とすと、質量と挿入の積に正しい尺度依存が残らない。真空からの混合と、挿入点の接触項は別に規格化している。

:::example[対称な二つの質量の和と差]
$\lambda_\phi=\lambda_\chi=\lambda$ なら和と差の固有方向で、次元付き質量二乗の固有係数は $(\lambda\pm h)/(16\pi^2)$。例えば $h=\lambda/2>0$ なら和の係数が差の三倍。無次元質量で書くとそれぞれから $2$ を引く。弱結合では両方とも赤外へ増大するが、二つの増え方が同じとは限らない。局所挿入の和と差にも同じ固有異常次元が現れ、それぞれを測る源は異なる。
:::

## 質量パラメータの発展と重い粒子の decoupling

対項 :eqref[eq-qft-phi4-counterterms] と $\delta Z=0$ を使い、$m_0^2=m^2[1+\lambda/(32\pi^2\bar\epsilon)]$ の固定微分を取る。先頭の $\beta_\lambda=-2\epsilon\lambda$ を掛けると $\beta_{m^2}=\lambda m^2/(16\pi^2)+\cdots$。従って無次元 $r=m^2/\mu^2$ では $\beta_r=(-2+\lambda/(16\pi^2))r+\cdots$。固定点近傍の $r$ は赤外で増える。

重い粒子が質量 $M$ を持つと、外部運動量が $M$ よりはるかに小さい振幅では重い内部線による非局所効果を局所展開できる。例えば :eqref[eq-qft-vacuum-polarization] の重い Dirac 場の寄与は物理的な零運動量減算後に $\alpha Q^2/(15\pi M^2)$ から始まる。しかし質量と独立な極だけを引く $\overline{\rm MS}$ の $\beta$ には重い粒子が低尺度でも寄与する。そこで $\mu\simeq M$ で測定振幅を一致させ、重い場を除いた理論へ接続する。これが :term[デカップリング]{en="decoupling" yomi="でかっぷりんぐ"} の実装である。重い粒子の結合が $M$ とともに増える場合、単純に $1/M$ の抑制を仮定できない。

### 有限な閾値照合の模型

低エネルギーの逆電荷を $a^{-1}$ とし、重い場のある理論での一ループ係数を $b_{\rm H}$、ない理論での係数を $b_{\rm L}$ とする。各範囲で $\dd a^{-1}/\dd\log\mu=b_i$ を使い、閾値 $M$ で連続とする近似では
$$
a^{-1}(\mu)=a^{-1}(\mu_0)+b_{\rm H}\log(M/\mu_0)+b_{\rm L}\log(\mu/M)
$$
となる。例えば $b_{\rm H}=4,b_{\rm L}=1$、$\mu_0=10M,\mu=M/10$ では差は $-5\log10$。重い場を低い尺度まで使うと $-8\log10$ となり、閾値で理論を接続した値より $3\log10$ 小さい。正確な照合では同じ低エネルギー振幅を両理論で計算し、くりこみ方式に依存する有限項もこの境界条件に入れる。

### QED の閾値を一ループで照合する

重い電荷一種類の質量を $M$ とし、全理論の $\overline{\rm MS}$ 結合を $\alpha_{\rm H}(\mu)$ とする。:eqref[eq-qft-vacuum-polarization] の零運動量減算の代わりに極だけを引いた修正は、この章で使う交換振幅の規約で
$$
\Pi_{\rm H}(-Q^2)=\frac{2\alpha_{\rm H}}{\pi}
\int_0^1\dd x\,x(1-x)
\log\frac{M^2+x(1-x)Q^2}{\mu^2}.
$$
$Q\ll M$ で残る定数は $\alpha_{\rm H}\log(M^2/\mu^2)/(3\pi)$。軽い場からの寄与を双方で共通に取れば、$1/Q^2$ の交換の係数を一致させる条件は
$$
\alpha_{\rm L}(\mu)=\alpha_{\rm H}(\mu)
\left[1+\frac{\alpha_{\rm H}(\mu)}{3\pi}\log\frac{M^2}{\mu^2}\right]
+O(\alpha_{\rm H}^3).
$$
$Q^2/M^2$ の部分は別の局所演算子に対応する。この式を $\log\mu$ で微分すると、明示対数は $-2\alpha_{\rm H}^2/(3\pi)$ を与え、全理論の $\beta_\alpha$ に含まれる重い粒子の $+2\alpha_{\rm H}^2/(3\pi)$ を消す。したがって照合の境界条件を任意の近い尺度で課しても、予測の差は高次となる。$\mu=M$ で大きい対数を避けることができ、この次数の特定の方式では定数項も消える。一般の重い相互作用や高次では、同じ選択で有限な照合項まで零になるとは限らない。

### 重いスカラーの四次結合と軽い質量を照合する

Euclid の二実場に $h X^2\phi^2/4$ を入れ、$X$ の質量 $M$ を軽い外部運動量より十分大きく取る。$U=h\phi^2/2$ と書けば重い場の行列式の零運動量での展開は $U I_1(M^2)/2-U^2 I_2(M^2)/4+\cdots$。既出のループ積分の極だけを引いた有限部分から、同じ低エネルギーの二点と四点を再現する条件は
$$
\begin{aligned}
m_{\rm L}^2(\mu)&=m_{\rm H}^2(\mu)
  +\frac{h M^2}{32\pi^2}\left(\log\frac{M^2}{\mu^2}-1\right)+O(h^2),\\
\lambda_{\rm L}(\mu)&=\lambda_{\rm H}(\mu)
  +\frac{3h^2}{32\pi^2}\log\frac{M^2}{\mu^2}+\cdots .
\end{aligned}
$$
二点は $m^2\phi^2/2$ と比較し、四点は $\lambda\phi^4/4!$ と比較した。両理論の軽い場だけのループは共通に取った。例えば $\mu=M$ では四次のこの次数の差は零だが、軽い質量には $-hM^2/(32\pi^2)$ の有限な差が残る。尺度微分を取ると対数の係数から、全理論の $\beta_{m^2}$ にある $h M^2/(16\pi^2)$ と $\beta_\lambda$ にある $3h^2/(16\pi^2)$ が低理論で消える。デカップリングは高次元の非局所効果が抑制されるという主張であり、関連な質量パラメータまで重い質量の逆冪で変わるという意味にはならない。

### 閾値をまたいだ重軽の寄与の違いを追う

二つの重い Dirac 場が同じ電磁場へ結合し、質量 $M_2>M_1$、電荷 $q_i$ を持つとする。$\mu>M_2$ から $M_1$ 未満まで下げる場合、$\mu=M_2$ で二番目だけを外し、次に $\mu=M_1$ で一番目を外す。各 $\beta_\alpha$ の係数を $b_l$ と $b_i=2(q_i/e)^2/(3\pi)$ の和で書くと、逆結合の微分はその負である。$\mu=M_i$ で一ループの $\overline{\rm MS}$ の対数が零になる条件を使えば、
$$
\alpha^{-1}(Q)=\alpha^{-1}(Q_0)
-b_l\log\frac{Q}{Q_0}
-b_2\log\frac{M_2}{Q_0}
-b_1\log\frac{M_1}{Q_0},
\qquad Q<M_1<M_2<Q_0 .
$$
$Q$ の微分には $b_l$ しか残らないが、低エネルギーで測る電荷の境界値そのものには $M_1,M_2$ の対数が残る。全理論の結合を同じ $Q_0$ で固定して重い質量を変えると、低エネルギーでの逆電荷は変わり得る。逆に $Q$ の電荷を実験で固定した上で $Q\ll M_i$ の変化を見るなら、重い場の非局所効果は $Q^2/M_i^2$ から始まる。境界値をどの尺度で固定したかで、比較している物理量が違う。

### 質量と四次結合の連立の走りを解く

四次元の一成分で $\beta_\lambda=3\lambda^2/(16\pi^2)$、$\beta_{m^2}=\lambda m^2/(16\pi^2)$ の一ループを使う。$\lambda>0$ の範囲で二式を割ると $\dd\log m^2/\dd\lambda=1/(3\lambda)$。よって $m^2(\mu)=m^2(\mu_0)[\lambda(\mu)/\lambda(\mu_0)]^{1/3}$ である。ここでは $m^2>0$ を仮定したが、同じ線形方程式は負の $m^2$ の符号も保存する。$m^2=0$ は独立な解である。この方式の質量二乗と $\Lambda$ の二乗を直接比較するのでなく、無次元比 $r=m^2/\mu^2$ を使うと、
$$
\frac{r(\mu)}{r(\mu_0)}
=\left(\frac{\mu_0}{\mu}\right)^2
\left[1-\frac{3\lambda(\mu_0)}{16\pi^2}\log\frac{\mu}{\mu_0}\right]^{-1/3}.
$$
例えば $\mu<\mu_0$ では $\lambda$ と $m^2$ は減るが、$r$ は小結合域で増え得る。これは次元付き質量が減ることと、その尺度の測定において質量効果が目立つことが両立する例である。$r$ が一に近づくと質量を無視した赤外の近似は崩れ、同じ無質量の相関の式でさらに外挿できない。

$4-\varepsilon$ 次元で非零固定点への流れを追う場合にも、質量は積分で求まる。$\beta_\lambda=-\varepsilon\lambda+b\lambda^2$、$\beta_{m^2}=c\lambda m^2$ とすると、固定点にいない軌道では $\dd\log|m^2|/\dd\lambda=c/(b\lambda-\varepsilon)$。従って
$$
m^2(\mu)=m^2(\mu_0)
\left|\frac{b\lambda(\mu)-\varepsilon}{b\lambda(\mu_0)-\varepsilon}\right|^{c/b}.
$$
固定点にちょうど置いた軌道は別に $\lambda=\varepsilon/b$ を直接使い、$m^2(\mu)=m^2(\mu_0)(\mu/\mu_0)^{c\varepsilon/b}$。軌道が固定点へ近づく場合、$\lambda-\lambda_*$ の大きさが $\mu^\varepsilon$ で減るため、こちらも同じ先頭の冪を再現する。$r$ に戻せば、その固定点での関連指数が $2-c\varepsilon/b$ になる。

## 関連・周辺・無関連演算子

固定点で次元 $\Delta$ の局所演算子 $\mathcal O$ に $u\int\dd^d x\,\mathcal O$ を加える。長さを $b>1$ 倍再設定すると測度は $b^d$、演算子は $b^{-\Delta}$ なので $u'=b^{d-\Delta}u$。この固有指数を $y=d-\Delta$ とし、$y>0,=0,<0$ をそれぞれ :term[関連]{en="relevant" yomi="かんれん"}、:term[周辺]{en="marginal" yomi="しゅうへん"}、:term[無関連]{en="irrelevant" yomi="むかんれん"} と呼ぶ。量子補正や混合があれば演算子の固有方向で評価する。$y=0$ では $\beta$ の高次の符号が必要で、四次元の小さい正の $\lambda\phi^4$ は赤外で弱くなる。Gaussian 四次元で $\phi^2,\phi^4,\phi^6$ の $y$ は $2,0,-2$ である。相関長が一つの関連方向 $r$ のみで決まると仮定すると、$|r|b^y\sim1$ で停止する尺度から $\xi\propto |r|^{-1/y}$ になる。

### 周辺方向と対数補正

長距離への変数を $\ell=\log b$ とし、周辺の結合 $u$ が $\dd u/\dd\ell=-a u^2$、関連方向 $r$ が $\dd\log r/\dd\ell=y+c u$ を満たすとする。$a>0,u_0>0$ なら $u=u_0/(1+a u_0\ell)$、従って $r=r_0e^{y\ell}(1+a u_0\ell)^{c/a}$。$r$ が $1$ に達するまで流した尺度から $\xi\propto |r_0|^{-1/y}[\log(1/|r_0|)]^{-c/(ay)}$ が従う。実際 $r=1$ の式を解いて $e^\ell=|r_0|^{-1/y}(1+a u_0\ell)^{-c/(ay)}$ とし、$\ell\simeq\log(1/|r_0|)/y$ を対数因子に使った。例えば $a=c=1,y=2$ なら対数に $-1/2$ の冪が付く。

### 無関連方向が残す最初の補正

固定点の近くに一つの関連方向 $r$ と一つの無関連方向 $u$ があり、線形化で $r'=b^{y_t}r$、$u'=b^{-\omega}u$、$y_t,\omega>0$ とする。相関長が同次式 $\xi(r,u)=b\,\xi(b^{y_t}r,b^{-\omega}u)$ を満たすなら、$b=|r|^{-1/y_t}$ として
$$
\xi(r,u)=|r|^{-1/y_t}\xi(\operatorname{sgn} r,u|r|^{\omega/y_t}).
$$
第二引数が零で滑らかなら、その一次は $u|r|^{\omega/y_t}$ の相対補正になる。例えば $y_t=2$、$\omega=1$ なら、$|r|$ を $1/100$ にしたとき無関連方向の係数は初期値の $1/10$ に減る。観測量が第二引数の零で特異ならこの Taylor 展開は使えない。例えば Landau の一様ポテンシャル $r\varphi^2/2+u\varphi^4/4!$ は、$r<0$ で $\varphi^2=-6r/u$ に極小を持つ。Gaussian 固定点の四次項が無関連となる $d>4$ でも、その極小の式では $u$ を零にできない。これは無関連な係数が観測量の極限を特異にする場合である。

### 臨界面からの距離と有限体積

系の長さを $L$、一つの関連方向を $r$、場の固定点の次元を $\Delta_\phi$ とする。空間の二点相関を $C(x;r,L)$ と書き、同じ形の境界条件で尺度変換できるとする。このとき $C(x;r,L)=b^{-2\Delta_\phi}C(x/b;b^{y_t}r,L/b)$。$b=L$ と選ぶと $C(x;r,L)=L^{-2\Delta_\phi}\mathcal C(x/L,rL^{y_t})$。$rL^{y_t}$ が小さければ、実際の相関長が箱の外にあるため有限サイズの関数で臨界の近傍を比べる。例えば感受率 $\chi=\int\dd^dx\, C(x;r,L)$ が長距離の積分で支配されるときは、臨界で $\chi\propto L^{d-2\Delta_\phi}=L^{2-\eta}$。積分の下端の短距離に非普遍な接触項があれば、その項は別に扱う。$rL^{y_t}$ が大きい領域では $b\sim|r|^{-1/y_t}$ を選んだ形が適する。この二つはどの尺度で流れを止めるかの違いである。

### 固定点で感受率と相関長の指数を結ぶ

一つの関連な温度方向を $r$ とし、相関長は $\xi\propto |r|^{-\nu}$。固定点の場の次元が $\Delta_\phi=(d-2+\eta)/2$ なら、短距離を除いた連結二点は $C(x;r)=|x|^{-(d-2+\eta)} f_\pm(|x|/\xi)$ と書けると仮定する。ここで $\pm$ は $r$ の符号、$f_\pm$ は大きい引数で十分速く減衰する関数。感受率の長距離成分は
$$
\chi\sim\int_0^\infty \dd R\, R^{1-\eta} f_\pm(R/\xi)
\propto \xi^{2-\eta}
$$
であり、$\eta<2$ かつ小さい引数で $f_\pm$ が有限なら積分の下端も有限。従って感受率の臨界指数は $\gamma_{\rm sus}=\nu(2-\eta)$。例えば一成分の $4-\varepsilon$ 次元では $\nu=1/2+\varepsilon/12+\cdots$、$\eta=O(\varepsilon^2)$ から $\gamma_{\rm sus}=1+\varepsilon/6+\cdots$。$\gamma_{\rm sus}$ は相関を空間全体で積分した量の指数であり、場の $\gamma_\phi$ と区別する。危険な無関連演算子が $f_\pm$ の極限に必要なら、この単純な尺度形の仮定を再検討する。

### Gaussian 固定点で四次が危険な無関連変数になる例

$d>4$ の一様ポテンシャル $V=r\varphi^2/2+u\varphi^4/4!$ で $u>0$ とする。Gaussian の $y_r=2$、$y_u=4-d<0$。$r<0$ の極小は $\varphi_0^2=-6r/u$、有効ポテンシャルの極小の差は $V(\varphi_0)-V(0)=-3r^2/(2u)$。実際、二次と四次に $\varphi_0^2$ を代入するとそれぞれ $-3r^2/u,+3r^2/(2u)$。従って $r$ への指数は磁化に相当する $\varphi_0$ では $1/2$ だが、自由の $y_r$ と場の次元だけで $\Delta_\phi/y_r=(d-2)/4$ と見積もると一致しない。流した後に $u$ を零と置くと、極小を決める分母を失うためである。

尺度 $b=|r|^{-1/2}$ まで流すと $u'=b^{4-d}u$。$\varphi_0=b^{-(d-2)/2}\sqrt{6/u'}$ と書いて $u'$ を戻せば $b^{-1}\sqrt{6/u}=\sqrt{6|r|/u}$。従って尺度変換自体が破れたのではなく、尺度を止めた後の関数が無関連な引数 $u'\downarrow0$ に対して特異である。この違いを確かめれば、相関長の指数と秩序変数の指数に単純な関係を使える仮定がどこにあるか分かる。

## 演習問題

::::exercise{#exr-qft-renormalization-group-1 level=basic}
四次元の $m=0$ の薄い殻で $\int_{\rm shell}\dd^4k/[(2\pi)^4 k^4]$ を求めよ。

:::solution
角積分をして $2\pi^2/(2\pi)^4$ を掛け、径方向の $\int_{\Lambda/b}^{\Lambda}\dd k/k=\log b$ より $\log b/(8\pi^2)$ である。
:::
::::

::::exercise{#exr-qft-renormalization-group-2 level=basic}
Gaussian 四次元で $\phi^6$ の $y$ を求めよ。

:::solution
$[\phi]=1$ から $\Delta=6$ なので $y=4-6=-2$ である。
:::
::::

::::exercise{#exr-qft-renormalization-group-3 level=basic}
:eqref[eq-qft-log-running] で $bg(\mu_0)\log(\mu/\mu_0)=1/4$ とする。$g(\mu)/g(\mu_0)$ を求めよ。

:::solution
分母が $1-1/4=3/4$ なので比は $4/3$ である。
:::
::::

::::exercise{#exr-qft-renormalization-group-4 level=basic}
一ループで $\gamma_\phi=0$ となる理由を述べよ。

:::solution
二点の一ループは同一点で閉じるタドポールで外部運動量が流れない。二点の $p^2$ の発散がないため場のくりこみはこの次数では不要である。
:::
::::

::::exercise{#exr-qft-renormalization-group-5 level=advanced}
$\gamma=\diag(\gamma_1,\gamma_2)$ が一定のときの $C_i(\mu)$ を求めよ。

:::solution
$\dd C_i/\dd\log\mu=C_i\gamma_i$ を積分すると $C_i(\mu_0)(\mu/\mu_0)^{\gamma_i}$ になる。演算子は逆の冪で発展するため積は一定である。
:::
::::

::::exercise{#exr-qft-renormalization-group-6 level=advanced}
四次スカラーの非零固定点での $\beta_\lambda'$ を求めよ。

:::solution
$\beta=-\varepsilon\lambda+b\lambda^2$ で $\lambda_*=\varepsilon/b$。微分は $-\varepsilon+2b\lambda_*=\varepsilon$ である。
:::
::::

::::exercise{#exr-qft-renormalization-group-7 level=advanced}
無次元質量の Gaussian 近傍で $\mu$ を半分にする。$r$ は何倍か。

:::solution
$\beta_r=-2r$ なので $r(\mu)=r(\mu_0)(\mu/\mu_0)^{-2}$。比は $(1/2)^{-2}=4$ である。
:::
::::

::::exercise{#exr-qft-renormalization-group-8 level=advanced}
$g'=g+cg^2$、$\beta_g=bg^2+\mathcal O(g^3)$ のとき $\beta_{g'}$ の先頭係数を求めよ。

:::solution
連鎖律で $(1+2cg)bg^2+\mathcal O(g^3)=bg^2+\mathcal O(g^3)$。逆変換の先頭は $g=g'+\mathcal O(g'^2)$ なので $b$ がそのまま残る。
:::
::::

::::exercise{#exr-qft-renormalization-group-9 level=advanced}
$\dd g/\dd t=-bg^3$ と $\dd O/\dd t=-ag^2 O$ が一定の $a,b>0$、$t=\log(\mu/\mu_0)$ で成り立つとする。$g(0)=g_0>0$ から $g(t),O(t)$ を求め、$a=2b$ のとき $g(t)=g_0/2$ となる尺度で $O$ の比を求めよ。

:::solution
$\dd(g^{-2})/\dd t=2b$ なので $g(t)^2=g_0^2/(1+2b g_0^2t)$。次に $\int_0^t g(s)^2\dd s=(1/(2b))\log(1+2b g_0^2t)$ を代入すると
$$
O(t)=O(0)(1+2b g_0^2t)^{-a/(2b)} .
$$
指定の点では $1+2b g_0^2t=4$、$a/(2b)=1$ だから比は $1/4$。係数 $C$ を反対の符号で発展させれば $C(t)=4C(0)$ で積は一定である。解を赤外へ延長するなら分母が正で、結合が摂動領域にある範囲に限る。
:::
::::

::::exercise{#exr-qft-renormalization-group-10 level=advanced}
一ループの二場模型で $x=3/10,y=1/10$ の回転対称な固定点における安定性行列を対角化せよ。$x,y$ の定義は本文のものを使う。

:::solution
行列を $\varepsilon$ で割ると $\begin{pmatrix}4/5&3/5\\1/5&2/5\end{pmatrix}$。固有値一のベクトルは $(3,1)^\mathsf T$。これは $y=x/3$ に沿って回転対称性を保つ。固有値 $1/5$ のベクトルは $(1,-1)^\mathsf T$。これは回転対称な比を破るが、赤外でそのずれは $\mu^{\varepsilon/5}$ として減衰する。二方向とも固有値が正なのは $\log\mu$ で微分する規約であり、長さを拡げる指数では負になる。
:::
::::

::::exercise{#exr-qft-renormalization-group-11 level=advanced}
一ループの質量を保った $\phi^4$ 理論で、基準尺度から下げて四次結合が $1/8$ になったとする。四次元で両尺度とも同じ質量非依存の方式を使う範囲で、$m^2$ と $m$ の比を求めよ。

:::solution
連立式を割ると $\dd\log m^2=(1/3)\dd\log\lambda$ だから $m^2(\mu)/m^2(\mu_0)=(1/8)^{1/3}=1/2$。正の質量そのものの比は $1/\sqrt2$ である。無次元比 $r$ についてはさらに $\mu_0^2/\mu^2$ が掛かる。この大きな尺度変更で実際に低エネルギーの観測を記述する場合、質量より十分低い領域の照合も別に考える。
:::
::::

::::exercise{#exr-qft-renormalization-group-12 level=advanced}
一成分の運動量殻の局所ポテンシャルで、四次の弱い結合 $\lambda$ しかない出発点から、殻積分の後に初めて現れる $\phi^6/6!$ の係数を求めよ。一様背景の一ループを使う。

:::solution
背景依存質量は $m^2+\lambda\varphi^2/2$。対数の $(\lambda\varphi^2/2)^3$ の項は $+\lambda^3\varphi^6/[3\cdot8(k^2+m^2)^3]$、実場の行列式にさらに $1/2$ が掛かる。従って $6!$ で規格化した係数は $+\lambda^3(6!/48)\int_{\rm shell}\dd^d k/[(2\pi)^d(k^2+m^2)^3]=15\lambda^3 I_3^{\rm shell}$。これは短い殻で新たに許される局所項であり、場の再スケーリングによる正準次元の項はまだ加えていない。四次元でこの係数の次元は $-2$ である。
:::
::::

::::exercise{#exr-qft-renormalization-group-13 level=advanced}
二場の質量の発展で、対角結合が同じ $\lambda$、混合結合が $h=\lambda/2$ とする。質量二乗の和と差の一ループの発展率を比べよ。

:::solution
質量行列を $(1,1)^\mathsf T$ と $(1,-1)^\mathsf T$ に作用させると固有値は $(\lambda+h)/(16\pi^2)$ と $(\lambda-h)/(16\pi^2)$。従って和の尺度微分に掛かる係数は $3\lambda/(32\pi^2)$、差では $\lambda/(32\pi^2)$ で、前者が三倍。この比は次元付きの質量二乗の発展率である。無次元化すると両方の係数から二を引くので、無次元の関連指数どうしには同じ三倍の関係はない。
:::
::::

## 確認テスト

::::quiz
Wilson の作用を低尺度へ移すときに最初にする操作を選べ。

- [x] 高運動量モードを積分する
- [ ] 低運動量モードを無条件に零にする
- [ ] 外部粒子の全状態を和にする

:::explanation
低モードを引数として固定し、高モードを積分する。得られる重みが低モードに対する新しい有効作用を定める。
:::
::::

::::quiz{answer="2" tol="1e-6"}
四次元 Gaussian 固定点における $\phi^2$ の $y$ を答えよ。

:::explanation
$\phi$ の次元は一なので $\Delta=2$。$y=4-2=2$ である。
:::
::::

::::quiz
質量付き重い場を低エネルギー理論から除く際に必要な操作を選べ。

- [x] 重い尺度近傍で係数を照合する
- [ ] 低エネルギーの全ての局所演算子を零にする
- [ ] 外部運動量を重い極へ置く

:::explanation
重い場の非局所性は低運動量で局所項に置き換わる。両理論の振幅が一致する条件でその係数を決める。
:::
::::

::::quiz{answer="4" tol="1e-6"}
固定点で $y=2$ の摂動を長さの倍率 $b=2$ で再設定する。係数の倍率を答えよ。

:::explanation
線形化では $u'=b^yu$ なので倍率は $2^2=4$ である。
:::
::::

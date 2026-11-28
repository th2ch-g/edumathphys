---
title: "線形化重力と重力波"
description: "Einstein 方程式の一次展開を行い、重力波の偏光、潮汐応答、二次の平均エネルギーを導く。"
notation:
  - symbol: "h_{ab}"
    meaning: "計量摂動"
  - symbol: "\\bar h_{ab}"
    meaning: "跡反転摂動"
  - symbol: "\\Box"
    meaning: "平坦背景の波動演算子"
  - symbol: "h_+,h_\\times"
    meaning: "二つの偏光振幅"
---

## この章の目標

- 計量摂動の一次 Einstein テンソルを計算する。
- 座標の自由度を用いて二つの偏光を分離する。
- 自由粒子の間隔の変化を測地線偏差から導く。
- 短波長での平均応力と真空曲率上の波動方程式を得る。

## 線形化 Einstein 方程式

$g_{ab}=\eta_{ab}+h_{ab}$、$\abs{h_{ab}}\ll1$ とし、添字の操作に $\eta$ を用いる。逆は $g^{ab}=\eta^{ab}-h^{ab}+\order{h^2}$ である。:eqref[eq-gr-christoffel] から $\Gamma^{(1)a}{}_{bc}=(\partial_b h^a{}_c+\partial_c h^a{}_b-\partial^a h_{bc})/2$ を得る。接続の積を捨て、$h=h^a{}_a$ とすると

$$
R^{(1)}_{ab}=\tfrac12(\partial_a\partial^c h_{cb}+\partial_b\partial^c h_{ca}-\Box h_{ab}-\partial_a\partial_b h),\quad
R^{(1)}=\partial^a\partial^b h_{ab}-\Box h.
$$

:term[跡反転]{en="trace reversal" yomi="あとはんてん"} $\bar h_{ab}=h_{ab}-\eta_{ab}h/2$ は四次元では二度施すと元に戻る。$C_b=\partial^a\bar h_{ab}$ とおけば

$$
G^{(1)}_{ab}=-\tfrac12\Box\bar h_{ab}+\partial_{(a}C_{b)}-\tfrac12\eta_{ab}\partial^c C_c .
\label{eq-gr-linear-einstein}
$$

偏微分が交換するので発散は零である。従って $\Lambda=0$ の一次方程式 $G^{(1)}=8\pi G T$ の源は背景上で $\partial^aT_{ab}=0$ を満たす。例えば定常な低速源では $\bar h_{00}=-4\Phi$、$\bar h_{ij}$ は低次で無視でき、$\Delta\Phi=4\pi G\rho$ に戻る。

### 線形化の次数と跡反転

$g_{ab}=\eta_{ab}+\epsilon h_{ab}$ とし、$\epsilon$ の次数を帳簿として扱う。逆行列を $g^{ab}=\eta^{ab}+\epsilon k^{ab}+\cdots$ と置いて積を取ると、一次の係数は $h^a{}_b+k^a{}_b$ である。従って $k^{ab}=-h^{ab}$ であり、接続の一次項では $g^{-1}$ の摂動は計量の微分と掛かって二次になる。四次元で $\bar h_{ab}=h_{ab}-\eta_{ab}h/2$ とすると $\bar h=-h$、従って $h_{ab}=\bar h_{ab}-\eta_{ab}\bar h/2$ である。これを曲率の一次式に入れ、$D_b=\partial^a\bar h_{ab}$ と置くと
$$
2G_{ab}^{(1)}=-\Box\bar h_{ab}
+\partial_aD_b+\partial_bD_a-\eta_{ab}\partial_cD^c .
$$
実際、$R_{ab}^{(1)}$ の各項のうち $\partial_a\partial_b\bar h$ は跡反転から生じた二項と $-\partial_a\partial_bh$ の間で消える。跡の補正を加えると $\eta_{ab}\Box\bar h$ も消える。さらに発散を取れば $-\Box D_b+\Box D_b+\partial_b\partial_cD^c-\partial_b\partial_cD^c=0$ である。従って外から与える線形源は $\partial_aT^{ab}=0$ を満たす必要がある。

:::example[静的な線形源で時間成分を戻す]
低速の質量源で $T_{00}=\rho$ が先頭、他成分をこの次数で零とする。遅延のない静的解では $\bar h_{00}=4G\int\rho(\vb y)/|\vb x-\vb y|\dd^3y=-4\Phi$ である。跡は $\bar h=4\Phi$ だから $h_{00}=-2\Phi$、$h_{ij}=-2\Phi\delta_{ij}$ となる。時間成分は加速度を、空間成分は光の距離を変えるため、光の伝播では同じ次数で両方を残す。
:::

### 遅延 Green 関数を原点で規格化する

Minkowski 背景の三次元波動演算子に対する遅延核を $G_{\rm ret}(t,\vb x)=\delta(t-r)/(4\pi r)$、$r=|\vb x|$ とする。$r>0$ で球対称場に $r$ を掛ければ $\partial_t^2-\partial_r^2$ が働くので、$\delta(t-r)$ は外向きの同次解である。原点の係数を確かめるには、滑らかな試験関数を掛けて $r\ge\epsilon$ で二回部分積分する。時間の境界は台が有界なら消え、球面の法線方向の項を評価する。$G_{\rm ret}$ の $1/r$ の微分は $-\delta(t-r)/(4\pi r^2)$ を含み、内側境界の向きと面積 $4\pi\epsilon^2$ により、極限で試験関数の $(0,0)$ での値を返す。$\delta'(t-r)/r$ と試験関数の空間微分の項には、面積を掛けた後で $\epsilon$ が残るため消える。従って $(\partial_t^2-\Delta)G_{\rm ret}=\delta(t)\delta^{(3)}(\vb x)$ である。この核を :eqref[eq-gr-wave-equation] に畳み込むと、源の過去光円錐上での値を積分する遅延解を得る。

:::example[一点で有限時間だけ動くスカラー源]
波動演算子 $\partial_t^2-\Delta$ の右辺を $S(t)\delta^{(3)}(\vb x)$ とし、$S$ の台を $0<t<T$ とする。解は $S(t-r)/(4\pi r)$ であるから、固定 $r$ で信号の台は $r<t<r+T$ で、後には残らない。この例はスカラーとして波動核だけを確認するものだ。重力の $\bar h_{ab}$ では同じ核を各成分に使うが、任意に時間変化する質量の点源をそのまま置くと $\partial_aT^{ab}=0$ を破る。波動方程式の因果性と、許される源の保存則を両方課す必要がある。
:::

### 背景に対する摂動の割り振り

平坦背景の Einstein テンソルは零なので、$g(\epsilon)=\eta+\epsilon h+\order{\epsilon^2}$ に対して一次の源も $T(\epsilon)=\epsilon T^{(1)}+\order{\epsilon^2}$ とする。有限の密度がある場を固定して、背景を平坦だと仮定したまま $\epsilon$ を零にすると、零次の方程式が成り立たない。平坦背景の線形源とは、$GT$ が摂動と同じ次数になるような弱重力近似の源である。物質の運動が保持されるための応力も $T$ に含める。例えば質量 $M(t)$ を持つ静止点源だけを $T^{00}=M(t)\delta^{(3)}$ と置けば、保存則は $\dot M=0$ を要求する。質量を変える過程を記述するには、エネルギーを運ぶ $T^{i0}$ とその運動量の収支が要る。

平坦の一点からの典型的な大きさを $d$、摂動の変化時間を $T$ とする。振幅 $h$ が小さくて $d/T$ が大きい場合、$h$ の時間微分が空間微分より重要になる。逆に非相対論的にゆっくり動く源では $d/T\ll1$ であるから、近接領域では空間微分が先頭になる。源の運動が低速であっても遠方の外向き波では $\partial_r h\simeq-\partial_t h$ となる。低速という条件は源のサイズと波長の比較であり、遠方の波の伝播速度を低くする条件ではない。

### 静的な球の内外で跡反転摂動を解く

総質量 $M$、半径 $R$ の一様球を弱重力 $GM/R\ll1$ で扱う。維持に必要な圧力は重力による Newton 平衡では $\rho GM/R$ の程度であるため、源の先頭を $T_{00}=\rho$ だけとする。球対称の Poisson 方程式
$$
\frac1{r^2}\dv{}{r}(r^2\bar h_{00}')=-16\pi G\rho
$$
を積分する。中心で $r^2\bar h_{00}'=0$ とし、外部で定数項が零となる条件を使う。$\rho=3M/(4\pi R^3)$ を入れると
$$
\bar h_{00}(r)=
\begin{cases}
2GM(3R^2-r^2)/R^3,&0\le r\le R,\\
4GM/r,&r\ge R .
\end{cases}
$$
値と一次微分が表面で連続なのは、体積の密度が跳んでいても、表面にデルタ関数の質量を置かなかったためである。$\bar h$ を戻すと $h_{00}=\bar h_{00}/2$、$h_{ij}=(\bar h_{00}/2)\delta_{ij}$ であり、中心の Newton ポテンシャルは $-3GM/(2R)$ になる。

:::example[内部では曲率が一定になる先頭近似]
内部で $h_{00}=GM(3R^2-r^2)/R^3$ である。静的な一次曲率は $R_{i0j0}=-\partial_i\partial_j h_{00}/2=(GM/R^3)\delta_{ij}$ となる。従って等方的に $\ddot\xi^i=-GM\xi^i/R^3$ の偏差が生じる。外部の静的な球の潮汐は動径伸長と横収縮の組で、真空で跡は零である。内部での跡 $3GM/R^3=4\pi G\rho$ は、Ricci の源として密度があることを表す。
:::

## ゲージ変換と Lorenz ゲージ

微小座標変換 $x'^a=x^a+\xi^a$ では $h'_{ab}=h_{ab}-\partial_a\xi_b-\partial_b\xi_a$ となる。跡反転すると $C'_b=C_b-\Box\xi_b$ である。局所の波動方程式 $\Box\xi_b=C_b$ を解いて :term[Lorenz ゲージ]{en="Lorenz gauge" yomi="ろーれんつげーじ"} $C_b=0$ を取れる。すると

$$
\Box\bar h_{ab}=-16\pi G T_{ab}.
\label{eq-gr-wave-equation}
$$

残留変換は $\Box\xi_b=0$ を満たす。例えば任意の $F(t-z)$ について $\xi_x=F(t-z)$ はこの残留自由度である。座標変換 $x'^a=x^a+\xi^a$ による $\delta h_{tx}=-F'$ と $\delta h_{zx}=F'$ は曲率を生まない。実際、$\Gamma^{(1)}$ の変化は $-\partial_b\partial_c\xi^a$ で、曲率を作る偏微分の差で消える。源がある領域で $C=0$ を取ることと、次節の真空の条件を全領域に課せるかどうかは区別する。

### 座標の変化と純ゲージの曲率

計量を引き戻す規約で $x^a\mapsto x^a+\epsilon\xi^a(x)$ とすると、平坦背景の変化は $h_{ab}=\partial_a\xi_b+\partial_b\xi_a$ である。具体的に $\xi_x=Ax^2/2,\xi_t=\xi_y=\xi_z=0$ なら $h_{xx}=2Ax$ のみ残る。一方で曲率の一次式はこの $h$ の二階微分の組なので零である。より一般の $\xi$ でも、各曲率成分の三階偏微分は符号が反対の項と相殺する。従って線形化曲率は平坦背景でのゲージ変換に不変である。曲がった背景では同じ引き戻しにより $\delta R$ が $\mathcal L_\xi\bar R$ だけ変わるため、背景の曲率がある場合には個々の摂動成分をそのままゲージ不変とは呼べない。

:::example[残った調和ゲージで縦成分を作る]
$\xi_a=B_a e^{ik\cdot x}$、$k^2=0$ とすると $\Box\xi_a=0$ である。座標変換 $x'^a=x^a+\xi^a$ による跡反転の変化は $-i(k_aB_b+k_bB_a-\eta_{ab}k\cdot B)e^{ik\cdot x}$ となる。$k^a$ と縮約すると $-ik^2B_b e^{ik\cdot x}$ だけが残り零である。従って Lorenz 条件を保ったまま縦成分が変わり得る。波の曲率はこの自由度では変わらないため、偏光の計数には条件を課した後に残るゲージも含める。
:::

### ゲージ条件を初期面から伝播させる

Lorenz ゲージで解くときは、十個の波動方程式の解に四つの拘束を保つ必要がある。保存源に対して :eqref[eq-gr-wave-equation] の発散を取ると $\Box C_b=-16\pi G\partial^aT_{ab}=0$ である。従って初期面 $t=t_0$ で $C_b=0,\ \partial_t C_b=0$ を課せば、同次波動方程式の解の一意性により全時刻で $C_b=0$ となる。コンパクトな初期データの一意性は、零の初期データを持つ二解の差 $C$ に対して $\int[(\partial_t C)^2+|\grad C|^2]\dd^3x/2$ を微分し、部分積分すれば示せる。保存される非負の積分が零なので、$C$ は定数であり、初期値からその定数が零になる。初期データが全空間で有界な範囲に限られない場合にも、各点の過去光円錐に限って同じ局所一意性を使える。

初期の $C_b$ は $\partial_t\bar h_{0b}=\partial_i\bar h_{ib}$ である。一方、$\partial_t C_b=0$ には波動方程式から二階の時刻微分を消して
$$
-\Delta\bar h_{0b}+\partial_i\partial_t\bar h_{ib}-16\pi G T_{0b}=0
$$
が要る。これは単に $\bar h_{0b}$ の初速を選ぶ条件と同じではない。既に一次 Einstein 方程式を満たす摂動から出発するなら、$\Box\xi_b=C_b$ を解いたあとに、この二条件は変換された場で満たされる。遅延解を採る場合には、源が過去に存在しないか、十分な減衰で部分積分ができれば、$C_b$ は保存則から直接零となる。

:::example[初期に $C_b$ だけを零にしても不十分]
真空で $\bar h_{00}(t,z)=A\cos(kz)\cos(kt)$、他成分を零とする。これは個別の波動方程式を解き、$t=0$ に $C_b=0$ を満たす。しかし $C_0=-\partial_t\bar h_{00}$ であるため、$\partial_t C_0(0,z)=Ak^2\cos(kz)$ である。$A\ne0,k\ne0$ なら初期面から外れると $C_0$ が現れ、元の Einstein 方程式は一般には満たされない。波動演算子を各成分に適用できることと、ゲージ条件が保たれることを分けて初期条件を用意する。
:::

### 時刻座標のずれから見かけの重力場を作る

$\xi^0=\varepsilon(\vb x)$、$\xi^i=0$ という時間に依らない微小座標変換では $\xi_0=-\varepsilon$ であり、$h'_{0i}=\partial_i\varepsilon$、他成分は零になる。$R^{(1)}_{ab}$ や線形化 Riemann に二階微分で入るが、微分の交換ですべて消える。静的な交差項を角速度のように読み取ると誤るのは、勾配部分の $\rot\vb h_0$ が零であるためである。反対に、定常で $\rot\vb h_0\ne0$ の項はこの時刻変換では消去できない。動く観測者の物理量を求めるときには、座標場とともに世界線や局所枠も同じ変換で移す。

:::example[純ゲージの縦波と実際の横波]
同じ遅延変数 $u=t-z$ で、$\xi_x=A\sin(\omega u)/\omega$ が作る真空摂動は $h_{0x}=-A\cos(\omega u)$、$h_{zx}=A\cos(\omega u)$ である。$R_{x0z0}^{(1)}$ に入る $\partial_z\partial_t h_{x0}$ と $-\partial_t^2h_{xz}$ は相殺する。これに $h_{xx}=-h_{yy}=B\cos(\omega u)$ を加えると、$R_{x0x0}^{(1)}=B\omega^2\cos(\omega u)/2$ が残る。$A$ をどの値に選んでも、この局所潮汐測定は変わらない。
:::

### 線形場の初期値にある四つの拘束

計量摂動を $h_{00}=-2\alpha$、$h_{0i}=\beta_i$、$h_{ij}=q_{ij}$ とし、$q=q_{ii}$ とする。これは時刻方向、時空の交差項、切片の空間計量の変化を分けた表記である。一次の Ricci と跡を $G^{(1)}_{00}=R^{(1)}_{00}+R^{(1)}/2$ に入れると、$\alpha$ と $\beta_i$、および二階の時間微分は消えて
$$
2G^{(1)}_{00}=\partial_i\partial_j q_{ij}-\Delta q
$$
となる。混合成分には跡の項がなく、
$$
2G^{(1)}_{0i}=\partial_j\dot q_{ji}-\partial_i\dot q
+\partial_i\partial_j\beta_j-\Delta\beta_i
$$
である。ここでも $q_{ij}$ の二階時間微分は出ない。従って初期のエネルギー密度と運動量密度は、空間の計量とその時間変化に制約を与える。これらの式を消して、全ての初期成分を独立な波として指定することはできない。

真空で空間の非零波数を $z$ 軸へ合わせれば、最初の式は $q_{xx}+q_{yy}=0$ を与える。$\beta_i=0$ の座標を取ると混合成分は $\dot q_{zi}=0$（$i=x,y$）、$\dot q_{xx}+\dot q_{yy}=0$ を与える。時刻に依らない $q_{zi}$ は、同じ波数をもつ空間ゲージ $\xi_i$ を用いて取り除ける。残る縦成分の扱いには残留する時刻の変換と空間成分の方程式も必要である。伝播する真空平面波に絞ると、次節の Lorenz 条件と残留変換がこれをまとめて行う。$\vb k=0$ の成分にはこの波数での除算が使えない。

:::example[束縛された密度のある静的な初期データ]
静的な $q_{ij}=-2\Phi\delta_{ij}$ なら $q=-6\Phi$ なので $2G^{(1)}_{00}=4\Delta\Phi$ である。従って $G^{(1)}_{00}=8\pi G\rho$ は $\Delta\Phi=4\pi G\rho$ と同じ式になる。この静的な空間摂動を $\rho$ を保ったまま零にすると、初期面だけで既に $00$ の方程式と矛盾する。源のある場所に真空の TT ゲージの条件を課せない理由を、この拘束からも見られる。
:::

## TT ゲージと重力波の偏光

非零周波数の真空平面波 $\bar h_{ab}=\Re(B_{ab}e^{i\omega(z-t)})$ を考える。以下では複素振幅で計算し、最後に実部を取る。波数は $k^2=0$、$k^a B_{ab}=0$ である。跡反転を元に戻した振幅を $A_{ab}$ とする。残留変換の振幅 $\epsilon_a$ で $\delta A_{00}=-2ik_0\epsilon_0$、$\delta A_{0i}=-i(k_0\epsilon_i+k_i\epsilon_0)$ なので、$k_0\ne0$ から四成分を消せる。$A_{0a}=0$ のもとで $k^a(A_{a0}-\eta_{a0}A/2)=0$ は $A=0$、空間成分の式は $A_{3i}=0$ を与える。これを :term[横波・無跡ゲージ]{en="transverse-traceless gauge" yomi="おうはむせきげーじ"}、略して TT ゲージという。残る二成分は

$$
h_{ij}^{\rm TT}=
\begin{pmatrix}
h_+&h_\times&0\\
h_\times&-h_+&0\\
0&0&0
\end{pmatrix}.
$$

横面の基底を $\psi$ 回すと $(h'_+,h'_\times)=(h_+\cos2\psi+h_\times\sin2\psi,-h_+\sin2\psi+h_\times\cos2\psi)$ である。これは両添字を回転させた結果である。例えば $\psi=\pi/4$ では純粋なプラス偏光は符号付きのクロス偏光として見える。

### 任意方向の波を横面へ射影する

進行方向を単位ベクトル $n$ とすれば $P_{ij}=\delta_{ij}-n_i n_j$ は横面への正射影である。任意の対称行列 $A$ に対し

$$
A^{\rm TT}_{ij}=P_i{}^kP_j{}^l A_{kl}-\tfrac12P_{ij}P^{kl}A_{kl}
$$

とすると $n^i P_{ij}=0$ で横方向に射影され、$P^i{}_i=2$ から跡が零になる。逆に横かつ無跡の $A$ には射影を施しても変わらない。従ってこの操作は横かつ無跡の部分空間への射影である。さらに対称行列の内積 $A:B=A_{ij}B_{ij}$ では、横かつ無跡の $B$ に対し $A^{\rm TT}:B=A:B$ である。この性質を使うと遠方の放射の角平均を計算できる。

:::example[斜め方向から見る]
$n=(1,1,0)/\sqrt2$、$A=\diag(2,0,0)$ とする。横枠 $e_1=(1,-1,0)/\sqrt2$、$e_2=(0,0,1)$ で $PAP$ は $\diag(1,0)$ であり、跡を半分ずつ引けば $\diag(1/2,-1/2)$ となる。従ってこの枠ではプラス偏光だけが残る。一方 $n=(1,0,0)$ なら横面に原行列の成分がないので TT 部分は零である。
:::

### 偏光基底を回す

横断面の基底を $(p,q)$ とし、$e_+=p\otimes p-q\otimes q$、$e_\times=p\otimes q+q\otimes p$ とする。$p'=p\cos\psi+q\sin\psi$、$q'=-p\sin\psi+q\cos\psi$ を代入すると
$$
e_+'=\cos2\psi\,e_++\sin2\psi\,e_\times,\qquad
e_\times'=-\sin2\psi\,e_++\cos2\psi\,e_\times
$$
である。テンソル自体は固定なので係数は逆向きに回り、$(h'_+,h'_\times)=(h_+\cos2\psi+h_\times\sin2\psi,-h_+\sin2\psi+h_\times\cos2\psi)$ である。従って角度 $\pi$ の基底回転で係数は元に戻り、$\pi/4$ の回転で純粋な $+$ は純粋な $\times$ になる。$\abs{h_+}^2+\abs{h_\times}^2$ は変わらないので、振幅の大きさを偏光軸から独立に比べられる。

:::example[横断面で直交する二つの方向の伸び]
$h_+=A,h_\times=0$ の瞬間に基線の方向を $n=p\cos\alpha+q\sin\alpha$ とすると $h^{\rm TT}_{ij}n^in^j=A\cos2\alpha$ である。$p$ 方向が最大の正、$q$ 方向が同じ大きさで負であり、$\alpha=\pi/4$ の基線では一次の伸びが零になる。これは波が消えたのではなく、基線が二つの主軸に等しい成分を持つためである。
:::

### 波数ごとの射影と空間全体の条件

空間 Fourier 変換で各 $\vb k\ne0$ に対し $n_i=k_i/|\vb k|$ とおく。空間摂動は二つのスカラー、二成分の横ベクトル、二成分の TT テンソルに分解できる。実際、$h_{ij}$ を $n$ と横面に分けると、$n^in^jh_{ij}$ が縦のスカラー、$P^{ij}h_{ij}$ が横の跡、$P_i{}^kn^jh_{kj}$ が横ベクトルである。残りは $P_i{}^kP_j{}^\ell h_{k\ell}$ の無跡部分になる。各成分を足すと
$$
h_{ij}=H_{\parallel} n_i n_j+\frac{H_\perp}2 P_{ij}
+n_i V_j+n_j V_i+h_{ij}^{\rm TT}.
$$
この式の各成分は互いに直交するため分解は一意である。座標変換の空間部分は Fourier 空間では $-i(k_i\xi_j+k_j\xi_i)$ である。各項に $k$ が入るので空間 TT 射影はこの変換を消す。平坦背景では $\xi_0$ は $h_{ij}$ の変化に直接入らないため、空間 TT 部分はこの意味でゲージ不変である。

逆 Fourier 変換すれば $P_{ij}=\delta_{ij}-\partial_i\partial_j\Delta^{-1}$ に相当する。逆 Laplace 演算子を含むため、空間全体の TT 部分を求める操作は局所的ではない。境界条件を指定せず一定成分 $\vb k=0$ を同じ式に入れることもできない。源の近くでは、TT 以外の成分にも密度や応力が作る測定可能な曲率がある。一方、外向き波の $1/r$ の先頭では、各方向について $\vb n$ を伝播方向として代数的な横射影を使える。静的な場や近接場を計算している場合と、遠方の放射での方向ごとの射影を区別する。

:::example[空間の TT 成分だけで静的な潮汐は作れない]
$h_{ij}=-2\Phi(\vb x)\delta_{ij}$ なら任意の非零波数に対し $P_i{}^kP_j{}^\ell h_{k\ell}=-2\widetilde\Phi P_{ij}$ で、横の跡を引くと零になる。ところが $h_{00}=-2\Phi$ があり $R_{i0j0}=\partial_i\partial_j\Phi$ は一般に非零である。従って「TT 部分が零」から一般の空間領域で物理効果が零と結論できない。
:::

### 二つのヘリシティを回転から読む

$e_+,e_\times$ を横の線形偏光の基底とする。複素結合 $e_{\rm R}=(e_++i e_\times)/\sqrt2$ と $e_{\rm L}=(e_+-i e_\times)/\sqrt2$ を取れば、横基底を $\psi$ 回したとき $e'_{\rm R}=e^{-2i\psi}e_{\rm R}$、$e'_{\rm L}=e^{2i\psi}e_{\rm L}$ になる。これは各添字が基底の回転に従うためで、位相の角度が二倍になる。波の向きを反転する場合は右手系の横枠も選び直すので、右・左の呼び方には伝播方向と時間位相の規約を添える。

:::example[楕円偏光の主軸と一周期の変形]
$h_+=A\cos\omega u$、$h_\times=B\sin\omega u$ とする。各瞬間の潮汐主軸は $2\alpha=\arg(h_++i h_\times)$ で表せる。$A=B>0$ なら $2\alpha=\omega u$ なので、主伸長方向は一周期に角 $\pi$ だけ回り、無向の軸として元に戻る。$B=0$ のときは主軸が固定され、$h_+$ が符号を変えると伸びと縮みが入れ替わる。一般に $A,B$ の比が半周期中の主軸の変化を定める。
:::

## 重力波が自由粒子に及ぼす効果

:ref[def-gr-riemann]{full} の符号で $R^{(1)}_{i0j0}=-\ddot h^{\rm TT}_{ij}/2$ である。TT ゲージでは $h_{0a}=0$ なので $\Gamma^i{}_{00}=0$ となり、最初に静止する自由粒子は座標上で静止したままである。しかし粒子間の固有間隔は変わる。長さが波長より短い配置で、自由落下の観測枠の分離を $\xi^i$ とすると :eqref[eq-gr-jacobi] が $\ddot\xi^i=\ddot h_{ij}^{\rm TT}\xi^j/2$ を与える。一次では右辺の $\xi$ に未摂動値 $\xi_0$ を使う。波の前に $h=\dot h=0$ かつ相対速度零なら $\xi^i=(\delta^i_j+h^i{}_j/2)\xi_0^j$ となる。腕の単位方向 $\vb n$ と長さ $L$ には $\delta L/L=h_{ij}n^in^j/2$ である。$L=4\,\mathrm{km}$、$h_+=10^{-21}$、$\vb n=\hat x$ なら $\delta L=2\times10^{-18}\,\mathrm m$ である。$\hat y$ の腕では同じ大きさで符号が逆になる。

### 曲率を二回積分して応答を得る

平行枠で $R_{i0j0}^{(1)}=-\ddot h_{ij}^{\rm TT}/2$ である。短い基線の初期分離を $\xi_0^i$ とすると、一次の偏差は $\delta\ddot\xi^i=\ddot h^i{}_j \xi_0^j/2$ に従う。積分定数も書くと
$$
\delta\xi^i(t)=\delta\xi^i(t_0)+(t-t_0)\delta\dot\xi^i(t_0)
+\tfrac12[h^i{}_j(t)-h^i{}_j(t_0)-(t-t_0)\dot h^i{}_j(t_0)]\xi_0^j.
$$
従って波の前に $h=\dot h=0$、装置も相対速度零であった場合に $\delta\xi^i=h^i{}_j\xi_0^j/2$ が得られる。$h$ が有限な別の定数へ移ると、曲率が再び零になっても分離の変位が残る。初期の定数を決めずに曲率を積分すると、これを初期の物体配置と取り違える。

:::example[有限時間に立ち上がるひずみ]
$h_+$ が $0$ から $A$ まで滑らかに変わり、その前後で $\dot h_+=0$ とする。$x$ 方向の分離 $L$ は終状態で $AL/2$ 増えるが、相対速度の変化は $[\dot h_+]L/2=0$ である。$y$ 方向では同じ量だけ縮む。曲率の時間積分が零でも、二回積分した変位は零とは限らない。
:::

### TT 座標の距離を落下観測者の枠に変換する

波が $+z$ に進み、中心の粒子は TT 座標で原点に静止するとする。原点に沿って $\tau=t$ であり、空間の正規直交枠は一次まで $e_{\hat i}{}^j=\delta_i^j-h_i{}^j/2$ である。$\Gamma^i{}_{0j}=\dot h^i{}_j/2$ なので時間微分と接続項が相殺し、この枠は平行輸送される。原点に直交して出る短い空間測地線上の近接点に対し、固有枠での分離は $X^{\hat i}=(\delta^i_j+h^i{}_j/2)\Delta x^j+\order{L^2\partial h}$ である。$\Delta x^j$ が一定なら、固有枠での $X^{\hat i}$ を二回微分して $\ddot X^{\hat i}=\ddot h^i{}_j\Delta x^j/2$ が得られる。

次の次数では、空間測地線の時刻成分も変わる。初速を原点の枠 $X^{\hat i}e_{\hat i}$ に取ると、終点の座標は $x^a=x^a_{\rm origin}+X^{\hat i}e_{\hat i}{}^a-\Gamma^a{}_{ij} X^{\hat i}X^{\hat j}/2+\cdots$ である。$\Gamma^0{}_{ij}=\dot h_{ij}/2$ だから $t=\tau-(\dot h_{ij}/4)X^{\hat i}X^{\hat j}+\cdots$ となり、原点の Fermi 同時面と $t$ の同時面は一般にずれる。長さと波長の比を小さくしたときにしか、原点の一つの曲率を配置全体に使えない。

:::example[横向きの輪を同時に観察する]
波のない初期状態で $X_0^2+Y_0^2=L^2$ とし、$h_+=A,h_\times=0$ の瞬間を考える。$(X,Y)=((1+A/2)X_0,(1-A/2)Y_0)$ だから半軸は $L(1\pm A/2)$ である。面積比は $1+\order{A^2}$ となり、一次では変わらない。純粋な $h_\times=A$ は横枠を $\pi/4$ 回せば同じ半軸の組になり、伸びの主方向だけが回る。この輪で $\ddot X,\ddot Y$ を調べる場合、ひずみの符号だけでなく二階微分も必要である。
:::

### 相対速度があるときの記憶

$h_{ij}(t)$ の前後で $\ddot h=0$ であっても、前後で $\dot h$ が異なるなら運動の変化が残る。初期の分離を固定する先頭近似で、$\Delta\dot\xi^i=[\dot h^i{}_j]_{\rm before}^{\rm after}\xi_0^j/2$ である。終状態の $\dot h$ が一定かつ非零なら、この表現を無限に長い時間へ外挿すると $|h|\ll1$ を破る。有限時間の波を考えるときは、局所で曲率が消えた後の粒子の相対運動を、そこでの慣性枠と初期値から延長する。

:::example[二つの潮汐パルスを間隔を空けて入れる]
横方向の一成分について、曲率を $\delta$ 関数の極限で $\ddot h(t)=B[\delta(t)-\delta(t-T)]$ とする。先に滑らかで短いパルスと考え、その幅を $T$ より十分短くする。第一のパルスで $\delta\dot\xi=B\xi_0/2$ となり、第二のパルスで一次の相対速度は元に戻る。それまでに $\delta\xi=BT\xi_0/2$ の変位が蓄積する。$|BT|\ll1$ の範囲で、曲率の積分は零でも位置の記憶を持つ例である。
:::

### 潮汐の時間成分と空間成分を同じ波から作る

局所直交枠で線形 Riemann の成分を取り、$E_{ij}=R^{(1)}_{i0j0}$、$B_{ij}=\epsilon_i{}^{kl}R^{(1)}_{klj0}/2$ と定義する。ここで $\epsilon_{xyz}=1$ とする。真空の制約で $E$ の跡は零であり、$B$ の跡は循環恒等式により零である。波が $+z$ に進むとき、$h_{0a}=0$ の座標で
$$
R^{(1)}_{klj0}=\tfrac12\partial_t(\partial_k h_{lj}-\partial_l h_{kj}),\qquad
E_{ij}=-\tfrac12\partial_t^2h_{ij}
$$
である。$\partial_k h_{ij}=-n_k\partial_t h_{ij}$ を入れれば
$$
B_{ij}=\epsilon_i{}^{kl}n_k E_{lj}
$$
となる。従って $B$ は $E$ の第一添字に横面の九十度回転を施したものである。横面にある $E$ は対称無跡なので、回転後も対称無跡であり、ノルムが保存される。伝播方向を反転すれば同じ $E$ に対する $B$ の符号が反転する。時刻成分の潮汐だけでは、ある点で同じ $E$ を作る逆向きの波を区別できない。

:::example[同じ伸縮でも伝播方向を見分ける]
$E_{xx}=C$、$E_{yy}=-C$、$E_{xy}=0$ の瞬間に、$+z$ へ進む波では $B_{yx}=B_{xy}=C$ である。$-z$ へ進む波に同じ $E$ を指定すると $B_{xy}=-C$ である。この違いは $R_{klj0}$ の混合成分にある。二つの相対する波を重ねれば、場所と時刻によって $E$ と $B$ の相殺の仕方も異なる。
:::

## 重力波のエネルギー運動量（Isaacson テンソル）

:term[Isaacson テンソル]{en="Isaacson tensor" yomi="あいざくそんてんそる"} は、波長 $\lambda$ より大きく背景の変化尺度 $\mathcal L$ より小さい領域で、二次の Einstein テンソルを平均して定義する。$\lambda\ll\mathcal L$ とし、平均で全微分の寄与を境界補正として捨てる。背景を局所 Minkowski とし、先頭の TT 波を用いる。二次 Ricci の接続微分は全微分であり、$\Gamma^{(1)c}{}_{cd}=0$ なので

$$
\langle R_{ab}^{(2)}\rangle=-\langle\Gamma^{(1)c}{}_{bd}\Gamma^{(1)d}{}_{ac}\rangle .
$$

二つの接続を展開する。微分添字が別の因子の摂動添字と縮約する項は部分積分して $\partial^ch_{cd}=0$ を使うと消える。微分二本が縮約する項は部分積分で $\Box h=0$ を使うと消える。残る項は $\langle\partial_a h_{cd}\partial_b h^{cd}\rangle/4$ である。この跡も同じ波動方程式で零となるため

$$
t^{\rm GW}_{ab}=-\frac{\langle G^{(2)}_{ab}\rangle}{8\pi G}
=\frac1{32\pi G}\langle\nabla_a h^{\rm TT}_{cd}\nabla_b h_{\rm TT}^{cd}\rangle .
\label{eq-gr-isaacson}
$$

背景での式は短波長の先頭次数である。例えば $h_+=A\cos\omega(t-z)$、$h_\times=0$ なら $t_{00}=t^{0z}=A^2\omega^2/(32\pi G)$ である。$h_{xx}$ と $h_{yy}$ の両方が寄与し、正弦二乗の平均 $1/2$ と相殺する。平均を取らない局所の二次項には座標依存の全微分が含まれる。

### エネルギー流束を SI で見積もる

同じ方向へ進む二偏光を $h_+=A\cos(\omega t-kz)$、$h_\times=B\sin(\omega t-kz)$ とする。二つの偏光テンソルは縮約が零で、各二乗は $2$ であるから $\langle\dot h_{ij}\dot h^{ij}\rangle=\omega^2(A^2+B^2)$ となる。SI では時間微分を秒で取り、$c^3/(32\pi G)$ を掛ければ $\mathrm{W\,m^{-2}}$ での流束を得る。$A=10^{-21},B=0,f=\omega/(2\pi)=100\,\mathrm{Hz}$ とすると約 $1.6\times10^{-3}\,\mathrm{W\,m^{-2}}$ である。同じ振幅でも周波数を二倍にすれば流束は四倍となる。長い平均区間での値であり、単発の短い波形なら $\dot h^2$ を時間積分して単位面積当たりのエネルギーを求める。

:::example[平面波の一次の振幅と二次の収束を分ける]
横方向の尺度を $A(u),B(u)>0$ とする Rosen 形の計量を使う。$u$ は進行方向の零変数である。横断面の面積要素は $AB$ であり、真空の式が各尺度とその積に課す条件を読む。 座標の順序を $(u,v,x,y)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}0 & -1 & 0 & 0\\-1 & 0 & 0 & 0\\0 & 0 & A^{2}{\left(u \right)} & 0\\0 & 0 & 0 & B^{2}{\left(u \right)}\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}0 & -1 & 0 & 0\\-1 & 0 & 0 & 0\\0 & 0 & \frac{1}{A^{2}{\left(u \right)}} & 0\\0 & 0 & 0 & \frac{1}{B^{2}{\left(u \right)}}\end{matrix}\right]\\
\det g&=- A^{2}{\left(u \right)} B^{2}{\left(u \right)}
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{v}_{x x}$ | $- \frac{1}{2}\left(-\partial_{u}\left(A^{2}{\left(u \right)}\right)\right)$ | $A{\left(u \right)} A'{\left(u \right)}$ |
| $\Gamma^{v}_{y y}$ | $- \frac{1}{2}\left(-\partial_{u}\left(B^{2}{\left(u \right)}\right)\right)$ | $B{\left(u \right)} B'{\left(u \right)}$ |
| $\Gamma^{x}_{u x}$ | $\frac{1}{2 A^{2}{\left(u \right)}}\left(\partial_{u}\left(A^{2}{\left(u \right)}\right)\right)$ | $\frac{A'{\left(u \right)}}{A{\left(u \right)}}$ |
| $\Gamma^{y}_{u y}$ | $\frac{1}{2 B^{2}{\left(u \right)}}\left(\partial_{u}\left(B^{2}{\left(u \right)}\right)\right)$ | $\frac{B'{\left(u \right)}}{B{\left(u \right)}}$ |

曲率の計算では :eqref[eq-gr-curvature-components] の微分差と積差を分ける。第二列は $\partial_c\Gamma^a{}_{db}$、第三列は $-\partial_d\Gamma^a{}_{cb}$ であり、第四列は接続の積を符号付きで足したものである。三列を足して最後の列を得る。

| 成分 | 第一の微分 | 第二の微分 | 接続の積 | 結果 |
|---|---|---|---|---|
| $R^{v}_{x u x}$ | $A{\left(u \right)} A''{\left(u \right)} + \left(A'{\left(u \right)}\right)^{2}$ | $0$ | $- \left(A'{\left(u \right)}\right)^{2}$ | $A{\left(u \right)} A''{\left(u \right)}$ |
| $R^{v}_{y u y}$ | $B{\left(u \right)} B''{\left(u \right)} + \left(B'{\left(u \right)}\right)^{2}$ | $0$ | $- \left(B'{\left(u \right)}\right)^{2}$ | $B{\left(u \right)} B''{\left(u \right)}$ |

第一添字を計量で下ろすと $R_{abcd}=g_{ak}R^k{}_{bcd}$ となる。対の反対称性と交換対称性を用いて、非零の独立な組を下表に示す。

| 成分 | 結果 |
|---|---|
| $R_{u x u x}$ | $- A{\left(u \right)} A''{\left(u \right)}$ |
| $R_{u y u y}$ | $- B{\left(u \right)} B''{\left(u \right)}$ |

Ricci テンソルは第一・第三添字を縮約して得る。曲率が残る場合も相殺が起こるので、各項を落とさずに加える。

| 成分 | 縮約の項 | 代入 | 結果 |
|---|---|---|---|
| $R_{u u}$ | $R^{x}_{u x u}+R^{y}_{u y u}$ | $- \frac{B''{\left(u \right)}}{B{\left(u \right)}} - \frac{A''{\left(u \right)}}{A{\left(u \right)}}$ | $- \frac{A{\left(u \right)} B''{\left(u \right)} + B{\left(u \right)} A''{\left(u \right)}}{A{\left(u \right)} B{\left(u \right)}}$ |

逆計量と縮約すると Ricci スカラーとなる。さらに $G_{ab}=R_{ab}-g_{ab}R/2$ を成分ごとに計算する。

$$
\begin{aligned}
R&=\sum_{a,b}g^{ab}R_{ab}\\
&=0\\
&=0
\end{aligned}
$$

非零の Einstein テンソル成分は次のとおりである。

| 成分 | Ricci 項 | 跡の項 | 差 |
|---|---|---|---|
| $G_{u u}$ | $- \frac{A{\left(u \right)} B''{\left(u \right)} + B{\left(u \right)} A''{\left(u \right)}}{A{\left(u \right)} B{\left(u \right)}}$ | $0$ | $- \frac{A{\left(u \right)} B''{\left(u \right)} + B{\left(u \right)} A''{\left(u \right)}}{A{\left(u \right)} B{\left(u \right)}}$ |

アフィン径数の微分を点で表し、$w_a$ を座標 $a$ の速度成分とする。ここでは添字を下げる操作を行っていない。測地線の右辺 $-\Gamma^a{}_{bc}w_bw_c$ では $b\ne c$ の項を二度数える。

$$
\begin{aligned}
\ddot{u}&=0\\
\ddot{v}&=- w_{x}^{2} A{\left(u \right)} A'{\left(u \right)} - w_{y}^{2} B{\left(u \right)} B'{\left(u \right)}\\
\ddot{x}&=- \frac{2 w_{u} w_{x} A'{\left(u \right)}}{A{\left(u \right)}}\\
\ddot{y}&=- \frac{2 w_{u} w_{y} B'{\left(u \right)}}{B{\left(u \right)}}
\end{aligned}
$$

また計量ノルムは

$$
g(\dot x,\dot x)=- 2 w_{u} w_{v} + w_{x}^{2} A^{2}{\left(u \right)} + w_{y}^{2} B^{2}{\left(u \right)}
$$
である。これを径数で微分して上の加速度を代入すると、係数の微分と速度の微分の和が零になる。従って初期のノルムが保存される。

$R_{uu}=-(A''/A+B''/B)$ である。$A=1+\epsilon f, B=1-\epsilon f$ とすると一次では $R_{uu}=0$ になるが、二次では $2\epsilon^2 f f''$ が残る。$A=1+\epsilon f+\epsilon^2p,B=1-\epsilon f+\epsilon^2q$ として二次まで真空条件を課せば $(p+q)''=2 f f''$ である。横面積の二次係数 $p+q-f^2$ は、従って二回微分すると $-2(f')^2$ になる。一次で反対符号の横伸縮が加わって面積変化を消しても、二次の収束が残ることが分かる。
:::

### 短波長での振幅・曲率・出力の尺度

波を $h_{ab}=\epsilon A_{ab}\cos(\varphi/\epsilon)$ と書き、背景と振幅・位相の変化の代表長さを $L$ とする。$h$ 自体は $\epsilon$ の次数だが、位相を微分するごとに $\epsilon^{-1}$ が出る。従って波の曲率の一次項は $\order{1/(\epsilon L^2)}$、二つの一階微分の積は $\order{1/L^2}$ である。後者は同じ次数の背景曲率に作用できる。一次の曲率が大きく振動しても、一周期で符号を変えるため粗視化すると平均の主項は二次から来る。例えば振幅 $A$ の単一調和成分では $h\partial_t^2 h=-\omega^2 A^2\cos^2\omega t$、$(\partial_t h)^2=\omega^2A^2\sin^2\omega t$ である。和は $\partial_t(h\partial_t h)$ で、整数周期で積分すれば零である。従って短波長の平均では、この二種類の二次項が逆符号で同じ流束の構成に現れる。有限な平均区間では端点の境界項が残るので、波長より十分長く背景尺度より短い窓を選ぶ。

### 二次の平均応力が満たす保存則

局所平坦背景で波動方程式を使い、
$$
\partial^a(\partial_a h_{cd}\,\partial_b h^{cd})
=(\Box h_{cd})\partial_b h^{cd}+\tfrac12\partial_b(\partial_a h_{cd}\,\partial^a h^{cd})
$$
と書く。第一項は真空で零である。第二項は、$\partial_a h_{cd}\partial^a h^{cd}$ 自体が真空中で $\partial^a(h_{cd}\partial_a h^{cd})$ であることから、短波長平均の先頭で消える。よって $\partial^at^{\rm GW}_{ab}=0$ が平均の次数で成り立つ。曲がった背景でも、位相の先頭で $t_{ab}^{\rm GW}$ は $k_ak_b$ に振幅の二乗を掛けた形になる。$k$ の測地線条件と束の輸送則を用いれば、対応する背景共変発散も先頭で零になる。

微小ゲージ変換を速い位相の振幅で考える。$h_{ab}$ の振幅に $\delta A_{ab}=k_a B_b+k_b B_a$ が加わるとき、真空波の横断条件を保つには跡反転した変化に $k^2=0$ が効く。残留変換のうち元の無跡性も保つものでは $k\cdot B=0$ である。このとき $A^{ab}\delta A_{ab}=0$、$\delta A_{ab}\delta A^{*ab}=0$ は $k^aA_{ab}=0$ と $k^2=k\cdot B=0$ からそれぞれ従う。従って位相勾配の積を掛けた平均応力の先頭は変わらない。無跡性が保たれない変換では TT 成分に射影し直す。同じ方向の TT 射影は $k_a B_b+k_b B_a$ の空間成分を消すので、物理的な二つの偏光は変わらない。この議論は高周波の順序で背景の分割を保つゲージ変換に関するものである。

:::example[逆向きの二つの波を平均する]
平坦背景で $h_+=A\cos\omega(t-z)+A\cos\omega(t+z)$ とする。各進行波の平均応力は同じエネルギー密度 $A^2\omega^2/(32\pi G)$ を持つ。時間と空間の両方を波長以上で平均すれば、逆向きの二位相の積は零になる。従って $t_{00}=t_{zz}=A^2\omega^2/(16\pi G)$、$t^{0z}=0$ である。流束は打ち消し合うがエネルギー密度は加わる。時間だけを平均する場合には、例えば $\sin\omega(t-z)\sin\omega(t+z)$ の平均に $\cos2\omega z/2$ が残るため、粗視化する範囲の指定が要る。
:::

### 短波長の波束を検出面で比較する

局所で単一方向に進むほぼ調和的な波を $h_+=A(t-z)\cos[\omega(t-z)]$ とし、包絡線の変化時間を $\Delta t$、$\omega\Delta t\gg1$ とする。時間微分では $A'\cos\omega(t-z)$ と $-\omega A\sin\omega(t-z)$ がある。短波長の先頭で後者が支配し、波束の単位面積当たりのエネルギーは
$$
\mathcal F=\frac{c^3\omega^2}{32\pi G}\int |A(t)|^2\dd t
$$
となる。球面状に広がり $A\propto1/R$ なら、球面一枚の積分では $R^2$ が相殺する。一方向での $\mathcal F$ を球面全体へ換算するには、偏光振幅の方向依存も必要である。

:::example[Gaussian 包絡線のフルエンス]
$A=A_0e^{-t^2/(2\sigma^2)}$ なら $\int |A|^2\dd t=A_0^2\sqrt\pi\sigma$ である。従って同じ周波数とピークひずみで $\sigma$ を二倍にすると、単位面積当たりのエネルギーも二倍になる。$\sigma$ が波周期に近づくと包絡線の微分も同じ次数になり、全波形の時間微分を積分する必要がある。振幅と周波数だけでは有限波束の総エネルギーは決まらない。
:::

### TT 二次作用で運動量の規格化を検算する

真空の平坦背景で $S=(16\pi G)^{-1}\int \sqrt{-g}R\,\dd^4x$ を使う。逆計量の変分に対し $\delta\sqrt{-g}=-\sqrt{-g}g_{ab}\delta g^{ab}/2$ である。また $g^{ab}\delta R_{ab}=g^{ab}(\nabla_c\delta\Gamma^c{}_{ab}-\nabla_b\delta\Gamma^c{}_{ac})$ は計量整合性により発散であり、変分の台を有界に取れば積分で消える。従って $\delta S=(16\pi G)^{-1}\int\sqrt{-g}G_{ab}\delta g^{ab}\,\dd^4x$ である。一般の境界条件は [作用原理と Hamilton 形式](/physics/general-relativity/action-formulation/) で扱う。$g=\eta+\epsilon h$ として $\partial_\epsilon g^{ab}=-h^{ab}+\order{\epsilon}$、$G_{ab}(\epsilon)=\epsilon G^{(1)}_{ab}+\order{\epsilon^2}$ だから
$$
\dv{S}{\epsilon}=-\frac{\epsilon}{16\pi G}
\int h^{ab}G^{(1)}_{ab}[h]\,\dd^4x+\order{\epsilon^2}.
$$
従って $\epsilon^2$ の係数は $-(32\pi G)^{-1}\int h^{ab}G^{(1)}_{ab}\dd^4x$ である。TT 条件を満たす自由場では $G^{(1)}_{ij}=-\Box h_{ij}/2$ であり、部分積分して
$$
S^{(2)}_{\rm TT}=\frac1{64\pi G}\int
[\dot h_{ij}\dot h_{ij}-(\partial_kh_{ij})(\partial_kh_{ij})]\,\dd^4x
$$
となる。ここでは真空の物理的な波の変数だけを用いた。源がある領域や境界条件を指定していないゲージでは、この作用だけから全部の Einstein 方程式を復元できない。

時間並進の保存エネルギー密度は $\mathcal E=[\dot h_{ij}^2+(\partial_kh_{ij})^2]/(64\pi G)$、流束は $\mathcal S^k=-\dot h_{ij}\partial_k h_{ij}/(32\pi G)$ である。実際、時刻で微分し $\ddot h=\Delta h$ を使うと $\partial_t\mathcal E=-\partial_k\mathcal S^k$ になる。一方向の進行波で $\partial_kh=-n_k\dot h$ を使えば $\mathcal S^k=n_k\mathcal E$ となり、平均後に :eqref[eq-gr-isaacson] と一致する。この一致は伝播する物理自由度の規格化を確認する。一般座標で局所の二次のエネルギーを定義する際には、部分積分の項やゲージの選び方に依存する。

:::example[同じ周波数でも直交する二偏光の応力]
$h_+=A\cos\omega(t-z)$、$h_\times=B\cos[\omega(t-z)+\delta]$ とする。空間添字の二重縮約は二偏光間で直交するため交差項が消え、$\mathcal E$ の周期平均は $\omega^2(A^2+B^2)/(32\pi G)$ である。従って $\delta=0$ の直線偏光も、$A=B,\delta=\pi/2$ の円偏光も、同じ $A^2+B^2$ なら平均の出力は同じである。偏光を決めるには平均エネルギーのほかに相対位相を測る。
:::

## 曲がった背景時空上の摂動

Ricci 平坦な背景 $\bar g$ で、$\Gamma^{(1)}$ の偏微分を背景共変微分に替えて $\delta\Gamma$ を得る。$\delta R_{ab}=\bar\nabla_c\delta\Gamma^c{}_{ab}-\bar\nabla_b\delta\Gamma^c{}_{ac}$ であり、微分の順を入れ替えると背景曲率が残る。例えば $\bar\nabla_c\bar\nabla_a h^c{}_b=\bar\nabla_a\bar\nabla_c h^c{}_b-\bar R^d{}_{bca}h^c{}_d$ である。Ricci の縮約は零で、残る二項は曲率の対称性により $-2\bar R_{acbd}h^{cd}$ になる。調和ゲージ $\bar\nabla^a\bar h_{ab}=0$ のもとで無源の式は

$$
\bar\Box\bar h_{ab}+2\bar R_{acbd}\bar h^{cd}=0.
$$

ここで $\bar h$ は背景で跡反転した量で、Ricci 平坦性から曲率との縮約で跡反転の計量項は消える。短波長で $\bar h=Ae^{i\varphi/\eps}$ と置くと $\eps^{-2}$ の次数が $k^ak_a=0$、$k_a=\partial_a\varphi$ を与える。$k$ は勾配なので $k^b\bar\nabla_b k_a=\bar\nabla_a(k^2)/2=0$ となる。次の次数は $2\bar\nabla_k A+(\bar\nabla_a k^a)A=0$ である。振幅のノルム二乗について $\bar\nabla_a(\abs A^2 k^a)=0$ が従う。球面波では $\nabla_a k^a=2/r$、$k^r=1$ なので $A\propto1/r$ となる。

:::example[正確な平面波計量と横方向の測地線]
線形近似と比較するため $\dd s^2=-2\dd u\,\dd v+\dd x^2+\dd y^2+a(x^2-y^2)\dd u^2$ を用いる。$a$ は逆長さの二乗を単位とする定数である。 座標の順序を $(u,v,x,y)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}a \left(x^{2} - y^{2}\right) & -1 & 0 & 0\\-1 & 0 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}0 & -1 & 0 & 0\\-1 & - a \left(x - y\right) \left(x + y\right) & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1\end{matrix}\right]\\
\det g&=-1
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{v}_{u x}$ | $- \frac{1}{2}\left(\partial_{x}\left(a \left(x^{2} - y^{2}\right)\right)\right)$ | $- a x$ |
| $\Gamma^{v}_{u y}$ | $- \frac{1}{2}\left(\partial_{y}\left(a \left(x^{2} - y^{2}\right)\right)\right)$ | $a y$ |
| $\Gamma^{x}_{u u}$ | $\frac{1}{2}\left(-\partial_{x}\left(a \left(x^{2} - y^{2}\right)\right)\right)$ | $- a x$ |
| $\Gamma^{y}_{u u}$ | $\frac{1}{2}\left(-\partial_{y}\left(a \left(x^{2} - y^{2}\right)\right)\right)$ | $a y$ |

曲率の計算では :eqref[eq-gr-curvature-components] の微分差と積差を分ける。第二列は $\partial_c\Gamma^a{}_{db}$、第三列は $-\partial_d\Gamma^a{}_{cb}$ であり、第四列は接続の積を符号付きで足したものである。三列を足して最後の列を得る。

| 成分 | 第一の微分 | 第二の微分 | 接続の積 | 結果 |
|---|---|---|---|---|
| $R^{v}_{x u x}$ | $0$ | $a$ | $0$ | $a$ |
| $R^{v}_{y u y}$ | $0$ | $- a$ | $0$ | $- a$ |

第一添字を計量で下ろすと $R_{abcd}=g_{ak}R^k{}_{bcd}$ となる。対の反対称性と交換対称性を用いて、非零の独立な組を下表に示す。

| 成分 | 結果 |
|---|---|
| $R_{u x u x}$ | $- a$ |
| $R_{u y u y}$ | $a$ |

Ricci テンソルは第一・第三添字を縮約して得る。曲率が残る場合も相殺が起こるので、各項を落とさずに加える。

| 成分 | 縮約の項 | 代入 | 結果 |
|---|---|---|---|
| $R_{u u}$ | $R^{x}_{u x u}+R^{y}_{u y u}$ | $- a + a$ | $0$ |

逆計量と縮約すると Ricci スカラーとなる。さらに $G_{ab}=R_{ab}-g_{ab}R/2$ を成分ごとに計算する。

$$
\begin{aligned}
R&=\sum_{a,b}g^{ab}R_{ab}\\
&=0\\
&=0
\end{aligned}
$$

Einstein テンソルはすべての成分で零となる。

アフィン径数の微分を点で表し、$w_a$ を座標 $a$ の速度成分とする。ここでは添字を下げる操作を行っていない。測地線の右辺 $-\Gamma^a{}_{bc}w_bw_c$ では $b\ne c$ の項を二度数える。

$$
\begin{aligned}
\ddot{u}&=0\\
\ddot{v}&=2 a w_{u} \left(w_{x} x - w_{y} y\right)\\
\ddot{x}&=a w_{u}^{2} x\\
\ddot{y}&=- a w_{u}^{2} y
\end{aligned}
$$

また計量ノルムは

$$
g(\dot x,\dot x)=- w_{u} w_{v} + w_{u} \left(a w_{u} \left(x^{2} - y^{2}\right) - w_{v}\right) + w_{x}^{2} + w_{y}^{2}
$$
である。これを径数で微分して上の加速度を代入すると、係数の微分と速度の微分の和が零になる。従って初期のノルムが保存される。

$\ddot u=0$ なので $\dot u\ne0$ なら $u$ をアフィン径数に取れる。このとき $x''=ax,y''=-ay$ となる。$a>0$ で二本の測地線を同じ横初速 $x'=y'=0$ で出すと、$x$ 側の初期分離は $\cosh(\sqrt a\,u)$、$y$ 側は $\cos(\sqrt a\,u)$ を掛けた値に変わる。$R_{uu}$ で相殺する曲率が、偏差では方向ごとに現れる。
:::

### 短波長近似で偏光を平行輸送する

波の振幅テンソルを $A_{ab}=\mathcal A e_{ab}$ とし、複素振幅のノルムが非零の横断偏光を $e_{ab}e^{*ab}=1$ に規格化する。上の輸送式と共役を縮約すれば、実の振幅 $\mathcal A$ は $k(\mathcal A)=-(\nabla\cdot k)\mathcal A/2$ に従う。元の輸送式へ戻すと $\nabla_k e_{ab}=0$ になる。横断条件は $\nabla_k k^a=0$ と合わせて伝播で保たれる。例えば光束の断面積がある区間で四倍になると、$\nabla\cdot k=k(\log A_{\rm screen})$ により $\mathcal A$ は半分となる。この近似では偏光テンソルを固定した座標軸で一定にする必要はなく、光線に沿って平行輸送した横断枠で成分が一定である。焦点で断面積が零になるとこの輸送式の振幅は発散するため、波長を無視する近似をそこで続けられない。

### 物質のある背景と真空背景の差

背景共変微分で $C_b=\bar\nabla^a(h_{ab}-\bar g_{ab}h/2)$ と置く。能動的な変化を $h_{ab}\mapsto h_{ab}+\bar\nabla_a\xi_b+\bar\nabla_b\xi_a$ として発散を取ると
$$
\delta C_b=\bar\Box\xi_b+[\bar\nabla^a,\bar\nabla_b]\xi_a
=\bar\Box\xi_b+\bar R_b{}^a\xi_a .
$$
第二項の符号は、共変添字に交換子を作用させたあと Ricci を縮約して正になる。従って Ricci 平坦な背景では残留条件は $\bar\Box\xi_b=0$ である。物質や宇宙項のある背景で同じ座標条件を使うときは Ricci 項を落とせない。また背景が物質を持つ場合、ゲージ変換で $\delta T$ に $\mathcal L_\xi\bar T$ が加わる。計量側だけを変換すると、線形化 Einstein 方程式の両辺が対応しない。

### 相対論的宇宙の膨張中を伝わるテンソル波

共形時刻 $\eta$ と平坦な空間切片で $\dd s^2=a^2(\eta)[-\dd\eta^2+(\delta_{ij}+\gamma_{ij})\dd x^i\dd x^j]$ を考える。ここでは背景の Einstein 方程式を再導出せず、等方圧を持つ一様な流体との解で $a(\eta)>0$ が指定されているとする。$\partial_i\gamma_{ij}=\gamma_{ii}=0$、流体の混合空間応力の TT 摂動がない場合を取る。背景との接続差を求めると $\delta\Gamma^0{}_{ij}=\mathcal H\gamma_{ij}+\gamma'_{ij}/2$、$\delta\Gamma^i{}_{0j}=\gamma'^i{}_j/2$、空間成分は平坦摂動の接続である。ここで $\mathcal H=a'/a$、プライムは共形時刻微分である。これらを Ricci に代入すると
$$
\delta R_{ij}=(\mathcal H'+2\mathcal H^2)\gamma_{ij}
+\mathcal H\gamma'_{ij}+\tfrac12(\gamma''_{ij}-\Delta\gamma_{ij}).
$$
混合成分では背景 $R_{ij}=(\mathcal H'+2\mathcal H^2)\delta_{ij}$ に逆計量の変分が掛かり、第一項を消す。スカラー曲率の跡の項も TT には寄与しない。従って $\gamma''_{ij}+2\mathcal H\gamma'_{ij}-\Delta\gamma_{ij}=0$ であり、波数 $k$ のモードで $\mu=a\gamma$ とすると $\mu''+(k^2-a''/a)\mu=0$ となる。ここでは背景を Ricci 平坦と仮定した式を流用せず、混合成分にして流体の応力の条件を使った。

:::example[放射優勢の尺度因子で伝わる波]
$a=A\eta$、$A>0,\eta>0$ とする。すると $a''=0$ で $\mu$ は $\cos(k\eta)$ と $\sin(k\eta)$ の線形結合になる。$\gamma$ はこれを $a$ で割った量である。$k\eta\gg1$ では振幅が $a^{-1}$、物理波長が $a$ に比例し、固有時の時間微分の先頭は $a^{-2}$ に比例する。従って局所で測る平均エネルギー密度は $a^{-4}$ に薄まる。この結果は波長が背景の曲率半径より十分短いときのエネルギーの解釈である。背景の膨張の由来は [相対論的宇宙論](/physics/general-relativity/cosmology/) で扱う。
:::

### 曲率がある領域で幾何光学の次の次数を読む

$\bar h_{ab}=e^{i\varphi/\epsilon}(A_{ab}+\epsilon B_{ab}+\cdots)$ とする。曲がった背景の波動演算子を展開すると、任意の振幅 $A$ に対し
$$
e^{-i\varphi/\epsilon}\bar\Box(e^{i\varphi/\epsilon}A)
=-\epsilon^{-2}k^2 A
+i\epsilon^{-1}[2\bar\nabla_k A+(\bar\nabla\cdot k)A]+\bar\Box A.
$$
$k^2=0$ を選んで先頭の輸送を解くと、次の $B$ には $i[2\bar\nabla_k B+(\bar\nabla\cdot k)B]=-\bar\Box A-2\bar R_{acbd}A^{cd}$ が入る。テンソル波の偏光や波束の形を、光線に沿う輸送だけで扱えるのはこの右辺の補正が小さい範囲である。ゲージ条件も同時に展開し、先頭で $k^a A_{ab}=0$、次で $\bar\nabla^a A_{ab}+ik^aB_{ab}=0$ を課す。真空で背景の曲率が非零のとき、位相速度の先頭が光速だからといって、有限波長で曲率の影響を全部捨てることはできない。

例えば球面状の波束で出発点からの面積が $A_{\rm screen}(s)$ で変わるとき、$\mathcal A(s)\sqrt{A_{\rm screen}(s)}$ は先頭で一定である。ただし小さい角幅の束を曲率が集めて断面が零に近づける場合、横方向の振幅微分が大きくなる。そこで $\bar\Box A$ を小さな補正として扱う順序は破れる。これは座標成分だけの特異性とは限らず、光線で一点に集中させた近似の特異性である。有限波長での場を得るには、異なる光線の位相も含む波動問題として扱う必要がある。

### 超長波長で宇宙膨張のモードを積分する

異方応力がない平坦な膨張背景で、前のテンソル方程式を固有時に直す。$\dd\eta=\dd t/a$ なので $\gamma'=a\dot\gamma$、$\gamma''=a^2(\ddot\gamma+H\dot\gamma)$ である。従って $\ddot\gamma_k+3H\dot\gamma_k+(k^2/a^2)\gamma_k=0$ になる。$k/a\ll H$ で勾配項を落とせば $\dd(a^3\dot\gamma_k)/\dd t=0$ であり、
$$
\gamma_k(t)=C_1+C_2\int^t\frac{\dd t'}{a(t')^3}.
$$
例えば $a\propto t^{2/3}$ の膨張では第二項は定数成分を引けば $t^{-1}$、$a\propto e^{Ht}$ なら $e^{-3Ht}$ で減衰する。一定のモードを、短波長で得られた $\gamma\propto1/a$ に置き換えてはいけない。短波長に入る時刻の近傍では両項を含む元の式を使う。

:::example[平坦な de Sitter のテンソル波の二つの極限]
$\eta<0$、$a=-1/(H\eta)$ なら $a''/a=2/\eta^2$ であり、$\mu''+[k^2-2/\eta^2]\mu=0$ である。$\mu=(1-i/(k\eta))e^{-ik\eta}$ を微分すると方程式を満たし、共役な解も独立である。$|k\eta|\gg1$ では $\mu$ はほぼ平面波で $\gamma=\mu/a$ の振幅は $1/a$ で減る。$|k\eta|\ll1$ では $\mu$ の先頭が $\eta^{-1}$ であり、$\gamma$ は定数に近づく。同じ解で両方を確認できる。後者は波の周期より短い時空尺度を平均する Isaacson の仮定を満たさない。
:::

## 演習問題

::::exercise{level=basic}
進行波 $e^{i(kz-\omega t)}$ が平坦真空波動方程式を満たす条件を求めよ。

:::solution
$\Box=-\partial_t^2+\partial_z^2$ を作用させると $(\omega^2-k^2)$ 倍になる。従って $\omega^2=k^2$ である。
:::
::::

::::exercise{level=basic}
$\bar h_{00}=4A$ だけが非零のとき $h_{00}$ と $h_{ij}$ を求めよ。

:::solution
跡は $\bar h=-4A$ である。$h_{ab}=\bar h_{ab}-\eta_{ab}\bar h/2$ から $h_{00}=2A$、$h_{ij}=2A\delta_{ij}$ となる。
:::
::::

::::exercise{level=basic}
$h_+=0,h_\times=A$ で $\vb n=(1,1,0)/\sqrt2$ の長さの相対変化を求めよ。

:::solution
$n^ih_{ij}n^j=2n^xn^y A=A$ である。従って $\delta L/L=A/2$ である。
:::
::::

::::exercise{level=basic}
プラス偏光 $h_+=A\cos\omega t$ の $R_{x0x0}$ を求めよ。

:::solution
$-\ddot h_{xx}/2=A\omega^2\cos(\omega t)/2$ である。
:::
::::

::::exercise{level=advanced}
静的な点質量の Lorenz ゲージの $\bar h_{00}$ を求めよ。

:::solution
:eqref[eq-gr-wave-equation] が $\Delta\bar h_{00}=-16\pi GM\delta^{(3)}$ となる。$\Delta(1/r)=-4\pi\delta^{(3)}$ と遠方での零条件より $\bar h_{00}=4GM/r$ である。
:::
::::

::::exercise{level=advanced}
同一方向の円偏光 $h_+=A\cos(\omega(t-z))$、$h_\times=A\sin(\omega(t-z))$ の $t_{00}$ を求めよ。

:::solution
横行列の二乗縮約は $2(\dot h_+^2+\dot h_\times^2)=2A^2\omega^2$ なので、:eqref[eq-gr-isaacson] から $A^2\omega^2/(16\pi G)$ である。
:::
::::

::::exercise{level=advanced}
源がない波で $\xi$ の残留変換が波動方程式とゲージ条件を保つことを示せ。

:::solution
$\delta C_b=-\Box\xi_b=0$ であり、$\Box$ を $\delta\bar h_{ab}$ に作用させると偏微分との交換で $\Box\xi$ の微分だけになるので零である。
:::
::::

::::exercise{level=advanced}
波が二倍遠い位置まで球面状に伝わるとき、波長が一定なら強度の比を求めよ。

:::solution
振幅輸送で $A\propto1/r$ である。周波数一定なら $t^{0r}\propto A^2$ なので強度は $1/4$ となり、球面積四倍と相殺する。
:::
::::

::::exercise{level=advanced}
$x,y$ 軸に沿う二本の基線と、$x$ 軸から $\pi/4$ 回した一つの基線で、角度を固定して同じ平面 TT 波を測る。一次の長さの相対変化を $s_0,s_{\pi/2},s_{\pi/4}$ とする。二つの偏光を復元し、三測定が整合する条件を求めよ。

:::solution
基線の単位ベクトル $n=(\cos\alpha,\sin\alpha,0)$ を使うと
$$
s_\alpha=\tfrac12 n^in^jh_{ij}^{\rm TT}
=\tfrac12(h_+\cos2\alpha+h_\times\sin2\alpha).
$$
従って $h_+=2s_0=-2s_{\pi/2}$、$h_\times=2s_{\pi/4}$ である。三測定の整合条件は $s_0+s_{\pi/2}=0$ で、これは横断二成分の跡が零であることに対応する。例えば $s_0=2\times10^{-21},s_{\pi/4}=-10^{-21}$ なら $(h_+,h_\times)=(4,-2)\times10^{-21}$ である。有限な光路の遅延がある場合は、比を取る前に基線に応じた時刻応答を入れる。
:::
::::

## 確認テスト

::::quiz
平坦背景の調和ゲージの残留変換に必要な条件を選べ。

- [x] $\Box\xi_a=0$
- [ ] $\xi_a=0$ が必須である。
- [ ] 真空の曲率が全て零である。

:::explanation
ゲージ関数は $C_a\mapsto C_a-\Box\xi_a$ と変わるので、もとの $C_a=0$ を保つには $\Box\xi_a=0$ でよい。
:::
::::

::::quiz{answer="2" tol="1e-9"}
非零周波数の単一方向の真空重力波の偏光の数を答えよ。

:::explanation
TT ゲージの横面には対称二階行列の三成分があり、無跡条件で一成分を引くため $2$ である。
:::
::::

::::quiz
TT ゲージで波によって変化し得る量を選べ。

- [x] 自由粒子の物理的な距離
- [x] 直交枠での潮汐曲率
- [ ] 初めに静止した自由粒子の空間座標が必ず変わる。

:::explanation
$h_{ij}$ は空間計量と $R_{i0j0}=-\ddot h_{ij}/2$ に現れる。一方、$\Gamma^i{}_{00}=0$ なので初めに静止した粒子の座標は一次で一定である。
:::
::::

::::quiz{answer="0.5" tol="1e-9"}
プラス偏光を $x$ 軸の短い腕で測るとき、$(\delta L/L)/h_+$ を求めよ。

:::explanation
腕方向では $n^in^jh_{ij}=h_+$ である。長さは計量による長さの二乗の平方根なので、一次の変化率は $h_+/2$ である。
:::
::::

::::quiz{answer="-0.25" tol="1e-7"}
純粋なプラス偏光を $\alpha=\pi/3$ で傾けた短い基線で測る。一次の $(\delta L/L)/h_+$ を求めよ。

:::explanation
横断の縮約は $h_+\cos2\alpha$ で、長さには計量の平方根による $1/2$ が掛かる。比は $\cos(2\pi/3)/2=-1/4$ である。
:::
::::

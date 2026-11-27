---
title: 微小振動
description: 平衡点での二次形式を質量で規格化し、連成振動の周波数と各質点の振幅を求める。
notation:
  - symbol: 'M,K'
    meaning: 平衡点での質量行列と剛性行列
  - symbol: '\eta'
    meaning: 平衡位置からの変位の列
  - symbol: 'a_s,Q_s'
    meaning: モード形と基準座標
  - symbol: '\omega_s'
    meaning: 基準角振動数
---

## この章の目標

- 線形化で残す項を、変位と速度の次数から選ぶ。
- 二つの二次形式を同時に対角化し、初期条件から各振幅を求める。
- 振動子、振り子、直線三原子分子でモードを物理的に読む。
- 零固有値、縮退、不安定方向を区別する。

## 平衡点まわりの展開

二つの振り子を弱いばねで結ぶと、一方を放した後でもう一方が揺れ始める。
各振り子は単独の周期を保たないが、二つを一緒に動かす特別な形には固有の周期がある。
その形と周期を同時に求める。

### 拘束を先に解く

:ref[def-am-generalized-coordinates] により定常なホロノミック拘束を解き、独立な座標 $q$ で $L=T-V$ と書く。速度依存力、駆動、散逸は当面含めない。静止した状態 $q_0$ で $\partial_i V(q_0)=0$ がすべての $i$ について成り立つとき、そこは :term[平衡点]{en="equilibrium point" yomi="へいこうてん"} である。直交座標で拘束を残して微分する場合、拘束面に垂直な力は零でなくてもよい。

質点の位置を $\vb r_b(q)$ とすれば、:eqref[eq-am-coordinate-velocity] から

$$
T=\frac12\sum_{i,j=1}^n M_{ij}(q)\dot q^i\dot q^j,\qquad
M_{ij}(q)=\sum_b m_b\frac{\partial\vb r_b}{\partial q^i}\cdot\frac{\partial\vb r_b}{\partial q^j}.
$$

$M$ は対称で、任意の列 $v$ に対して $v^\top Mv=\sum_bm_b|\sum_i v_i\partial_i\vb r_b|^2$ である。
正則な座標と正の質量なら、$v\ne0$ に対し正である。

### 二次までの Lagrangian

$\eta=q-q_0$ とし、$\eta,\dot\eta$ をともに微小な振幅パラメータの一次と数える。
$M(q)-M(q_0)$ の寄与は $\eta\dot\eta^2$ の次数だから三次である。
一方、$V$ の一次項は平衡条件で消える。
定数を落とした二次の Lagrangian は

$$
L_2=\tfrac12\dot\eta^\top M\dot\eta-\tfrac12\eta^\top K\eta,\qquad
K_{ij}=\partial_i\partial_jV(q_0),\quad M=M(q_0)
\label{eq-am-quadratic-lagrangian}
$$

である。
混合偏微分が交換できる $C^2$ 級の $V$ なら $K$ も対称である。
Lagrange 方程式は $\sum_j M_{ij}\ddot\eta^j=-\sum_jK_{ij}\eta^j$。
直交変位を m で測れば $M_{ij}$ は kg、$K_{ij}$ は $\mathrm{N/m}$、$M^{-1}K$ は $\mathrm{s^{-2}}$ を単位とする。
角度を含む座標では、行列の各成分が同じ単位とは限らない。

:::example[下向きと上向きの振り子]
長さ $l$ の振り子で下向きから $\theta$ を測ると $V=-mgl\cos\theta$。
平衡は $0,\pi$ である。
$M=ml^2$、$K=mgl\cos\theta_0$ だから、下では $\ddot\eta=-(g/l)\eta$、上では $\ddot\eta=(g/l)\eta$ である。
$l=0.25\,\mathrm m,\ g=9.8\,\mathrm{m/s^2}$ の下向き振動は $\omega\simeq6.26\,\mathrm{s^{-1}}$。
角度 $0.1$ rad なら $\sin\theta$ の線形近似の相対誤差は約 $\theta^2/6=0.0017$ である。
:::

:::example[二次の判定ができない平衡]
$V=\alpha x^4$ では $K=0$ で線形方程式は $\ddot x=0$。
$\alpha>0$ なら有限振幅で復元力 $-4\alpha x^3$ があり、$\alpha<0$ なら外向きである。
零固有値の物理的意味は対称性や高次項から決める必要がある。
:::

### 重力で移動した平衡点を選ぶ

鉛直ばねの自然長から下向きに測った変位を $z$ とする。
伸びが $z$ なら $V=kz^2/2-mgz$ である。
$z=0$ のまわりで展開すると一次項が残り、線形方程式には定数外力が現れる。
一方、$z_0=mg/k$ として $z=z_0+\eta$ と置けば

$$
V=\tfrac12 k\eta^2-\frac{m^2g^2}{2k},\qquad m\ddot\eta+k\eta=0.
$$

$k=98\,\mathrm{N/m}$、$m=1\,\mathrm{kg}$ なら平衡は自然長より $0.10\,\mathrm m$ 下で、角振動数は $\sqrt{98}\,\mathrm{s^{-1}}$。
重力はここでは平衡を移し、ばねの曲率を変えない。
剛性が伸びに依存するばねでは移動先で曲率を取り直すため、重力で周波数も変わりうる。

### 曲面上の曲率を正しく測る

振り子の重力ポテンシャル $V=mgz$ の空間 Hessian は零である。
それでも円弧に沿えば曲率が生じる。
拘束面の座標 $q$ で微分すれば

$$
\frac{\dd^2V(\vb r(q))}{\dd q^2}
=\vb r'^\top(\nabla^2 V)\vb r'+\nabla V\cdot\vb r'' .
$$

第二項は経路自体が曲がる効果である。
$z=-l\cos q$ なら下端で $mgz''=mgl$ を得る。
したがって空間 Hessian を接平面に単純に射影するだけではばね定数は求まらない。
平衡で拘束力が支えていた方向の勾配も必要である。

:::example[安定な平衡をエネルギーで囲む]
$M>0,\ K>0$ の自然系では、十分小さい領域で $V(q)-V(q_0)\ge c|\eta|^2$ となる $c>0$ が取れる。
これは Taylor の二次項の最小値が単位球面上で正だからである。
半径 $r$ の境界で $V$ は少なくとも $cr^2$ 高い。
初期の余剰エネルギーがそれより小さければ境界に達せず、運動エネルギーの正値性から速度も小さく保たれる。
これで元の非線形方程式についても平衡の安定性が分かる。
ただし摩擦がなければ平衡へ戻って静止するとは限らない。
:::

### 二次形式が退化する平衡と分岐

二次の位置エネルギーが正なら、十分小さいエネルギーをもつ軌道は極小のまわりに閉じ込められる。実際、平衡を囲む小球の表面で $V-V_*$ の最小値を取れば、正定値な Hessian と Taylor の定理によりこの値は正である。それより小さい全エネルギーでは球面に到達できない。質量行列も近傍で一様に正定値なので速度も小さい。これが静止状態の Lyapunov 安定性を与える。この論法には非線形の運動を線形の解で無期限に近似できるという仮定は要らない。

零固有値がある場合には、高次の項を調べる前に、拘束や硬い方向への変位を正しく処理する必要がある。二変数の静的ポテンシャル $V(x,y)$ で $V_{yy}>0$ とする。$V_y=0$ を満たす枝を $y=y_*(x)$ と書けば $y_*'=-V_{yx}/V_{yy}$ であり、枝に沿う曲率は $V_{xx}-V_{xy}^2/V_{yy}$。硬い方向を固定した測定と、力を零に保って緩和させた測定とでは柔らかさが違う。この消去は静的な判定であり、振動数を出すには枝の運動に伴う慣性も考えなければならない。

:::example[曲がった谷での四次の安定性]
適切な長さとエネルギーを単位にとり $V=(y-x^2)^2/2+\alpha x^4$ とする。原点の Hessian は $\diag(0,1)$ で、$\alpha$ の符号を含まない。$y=0$ の断面だけを調べると $(1/2+\alpha)x^4$ となるが、谷底 $y=x^2$ では $\alpha x^4$。従って $\alpha>0$ なら厳密な極小、$\alpha=0$ では連続な極小、$\alpha<0$ なら極小でない。例えば $\alpha=-1/4$ では両座標軸上の値は正でも、谷底に沿って原点より低くなる。高次まで調べるときも一本の直線上の検査だけで安定性を結論してはならない。
:::

### 対称な二本のばねの座屈

同じばねを $(\pm a,0)$ に固定し、他端を $y$ 軸上を動く質点に集める。自然長を $l$、各ばねの剛性を $k$ とする。左右の長さは等しいので $V=k(\sqrt{a^2+y^2}-l)^2$、力は $-2ky(1-l/\sqrt{a^2+y^2})$。常に $y=0$ が平衡であり、$l>a$ でのみ別の二つの平衡 $y=\pm\sqrt{l^2-a^2}$ がある。中心で $K_0=2k(1-l/a)$、両側の平衡で $K_s=2k(1-a^2/l^2)$ である。臨界長で二次の復元力は消えるが $V=ky^4/(4a^2)+O(y^6)$ なので、振幅に依存する非線形振動は残る。

:::example[自然長が 0.80 m の場合]
$a=1\,\mathrm m,k=10\,\mathrm{N/m},m=2\,\mathrm{kg}$ とする。中心を基準にしてポテンシャルを展開すれば、SI の数値で

$$
V(y)-V(0)=2 y^{2} \left(y^{2} + 1\right)+O(y^6).
$$

中心での曲率は $4\,\mathrm{N/m}$。線形角振動数は $\sqrt{2}\,\mathrm{s^{-1}}$ となる。二本は初めから張っている。横変位の一次では長さが変わらなくても、既存の張力の方向が変わるため復元力は一次になる。
:::

:::example[自然長が 1 m の場合]
$a=1\,\mathrm m,k=10\,\mathrm{N/m},m=2\,\mathrm{kg}$ とする。中心を基準にしてポテンシャルを展開すれば、SI の数値で

$$
V(y)-V(0)=\frac{5 y^{4}}{2}+O(y^6).
$$

中心での曲率は $0\,\mathrm{N/m}$。二次近似では周期を決められない。この場合に角振動数零という結果を入れて一定速度の運動を予測すると、有限振幅で蓄えられる四次の弾性エネルギーを取り落とす。
:::

:::example[自然長が 1.25 m の場合]
$a=1\,\mathrm m,k=10\,\mathrm{N/m},m=2\,\mathrm{kg}$ とする。中心を基準にしてポテンシャルを展開すれば、SI の数値で

$$
V(y)-V(0)=\frac{5 y^{2} \left(5 y^{2} - 4\right)}{8}+O(y^6).
$$

中心での曲率は $-5\,\mathrm{N/m}$。線形成長率は $\frac{\sqrt{10}}{2}\,\mathrm{s^{-1}}$ となる。側方の二点 $y=\pm \frac{3}{4}\,\mathrm m$ に移れば正の角振動数 $\frac{3 \sqrt{10}}{5}\,\mathrm{s^{-1}}$ を得る。中心で小さな横変位を与えると増大する。端の位置が固定されていても、対称な装置が片側の平衡を選ぶ。
:::

## 基準振動と基準座標

### 特別な振れ方から一般解へ

全座標が同じ時間関数に比例する解 $\eta=a\cos(\omega t+\delta)$ を :term[基準振動]{en="normal mode" yomi="きじゅんしんどう"} という。
代入すると

$$
Ka=\omega^2Ma,\qquad \det(K-\omega^2M)=0.
\label{eq-am-normal-pencil}
$$

後者は非零の $a$ が存在する条件である。
質量が異なると $M^{-1}K$ は通常の内積では対称でない。
振幅を比較するには質量で重みを付ける。

:::theorem[質量行列と剛性行列の同時対角化]{#thm-am-normal-modes}
実対称正定値 $M$ と実対称 $K$ に対して実の可逆行列 $A=(a_1,\ldots,a_n)$ が存在し、

$$
A^\top MA=I,\qquad A^\top KA=\diag(\lambda_1,\ldots,\lambda_n),\qquad Ka_s=\lambda_sMa_s
$$

を満たす。
$K$ が正定値ならすべての $\lambda_s$ は正である。
:::

:::proof
まず任意の実対称行列 $D$ を考える。
単位球面上で連続関数 $x^\top Dx$ は最大値をもつ。
その点 $u$ で任意の $w\perp u$ への方向微分は $2w^\top Du=0$ だから $Du$ は $u$ に平行である。
さらに $u^\top Dw=(Du)^\top w=0$ なので $u^\perp$ は不変である。
この空間上で同じ議論を次元が零になるまで繰り返せば、正規直交固有基底を得る。

これを $M=U\diag(m_s)U^\top$ に用いると、単位固有ベクトル上の二次形式の値 $m_s$ は正である。
そこで $C=U\diag(m_s^{-1/2})U^\top$ とすると $CMC=I$。
$CKC$ も実対称だから直交行列 $O$ により $O^\top CKCO=\diag(\lambda_s)$ となる。
$A=CO$ は所望の二式を満たす。
$A^\top KA=A^\top MA\diag(\lambda_s)$ の左から $(A^\top)^{-1}$ を掛けると固有方程式も得られる。
最後に $\lambda_s=a_s^\top Ka_s$ より正値性が従う。
:::

### 初期条件を射影する

:term[基準座標]{en="normal coordinate" yomi="きじゅんざひょう"} を $Q=A^{-1}\eta=A^\top M\eta$ と定義すると

$$
L_2=\tfrac12\sum_{s=1}^n(\dot Q_s^2-\lambda_sQ_s^2),\qquad \ddot Q_s=-\lambda_sQ_s .
$$

$\lambda_s>0$ なら $\omega_s=\sqrt{\lambda_s}$ とし、

$$
Q_s(t)= (a_s^\top M\eta(0))\cos\omega_st+
(a_s^\top M\dot\eta(0))\frac{\sin\omega_st}{\omega_s}.
$$

零固有値なら $Q_s(0)+t\dot Q_s(0)$、負なら $\cos,\sin/\omega$ を $\cosh(\sqrt{-\lambda_s}t),\sinh(\sqrt{-\lambda_s}t)/\sqrt{-\lambda_s}$ に替える。
縮退した固有空間では質量内積に関して任意の正規直交基底を選べる。
個別のモード形の選択は一意でないが、その空間への射影は一意である。

:::example[質量での規格化]
$M=\diag(1,4)\,\mathrm{kg}$ で $a=(2,1)^\top$ と同じ向きの振幅を使いたい。
$a^\top Ma=8\,\mathrm{kg}$ なので $a/\sqrt{8\,\mathrm{kg}}$ とする。
これと質量直交するベクトルは $2b_1+4b_2=0$ を満たす。
例えば $b=(-2,1)^\top$ の通常の内積は $a^\top b=-3$ だが、運動エネルギーに交差項は生じない。
:::

### 直交変換と合同変換の違い

実対称行列一つの対角化は :ref[thm-linalg-real-spectral] でも示されている。
二つの行列に同じ直交変換を使って両方を対角化できるなら、その二行列は互いに可換でなければならない。
実際、同じ基底で対角なら積の順序を交換でき、基底を戻しても可換である。
振動問題ではこの可換性を仮定しない。
変位の変換を二次形式の両側に代入する**合同変換**を使い、運動エネルギーを先に単位行列に直す。

:::example[重みを外して固有値を求める]
SI の数値で

$$
M=\begin{pmatrix}4&0\\0&1\end{pmatrix},\qquad
K=\begin{pmatrix}8&-2\\-2&2\end{pmatrix}
$$

とする。
$M$ と $K$ の積は非対角成分の左右の値が掛ける順序で入れ替わるので可換でない。
$C=\diag(1/2,1)$ により

$$
CKC=\begin{pmatrix}2&-1\\-1&2\end{pmatrix}
$$

となり、直交する $(1,1)/\sqrt2,(1,-1)/\sqrt2$ の固有値は $1,3$。
元の振幅は $C$ を掛けて $(1/2,1)/\sqrt2,(1/2,-1)/\sqrt2$ となる。
軽い質点は重い質点の二倍の振幅を持ち、$A^\top MA=I$ である。
単に $K$ の固有ベクトルを使うと、この振幅比は得られない。
:::

### 一つの質点にだけ力を加えた後の応答

$t=0$ の短い時間だけ外力を加え、その力積を列 $j$ にまとめる。
衝撃の間に変位がほとんど変わらないとすると、運動方程式を時間積分して $M[\dot\eta(0^+)-\dot\eta(0^-)]=j$ となる。
静止した平衡から押したとき

$$
Q_s(0^+)=0,\qquad \dot Q_s(0^+)=a_s^\top j .
$$

従って運動量の強制は、変位の初期条件とは異なる重みでモードに投影される。
例えば等質量の二振動子に逆向きで等しい力積 $(j,-j)$ を与えると、同位相形との内積は零である。
二つの質点が動き出しても同位相振動にはエネルギーが入らない。

### 質量で規格化したモードから周波数応答を復元する

動かす点と測定点が異なると、ある振動数では外力に対する変位が消える。この零点を :term[反共振]{en="antiresonance" yomi="はんきょうしん"} という。粘性のない系に調和力 $\vb f\cos\Omega t$ を加え、同じ時間依存をもつ特解 $\vb d\cos\Omega t$ を探す。$z=\Omega^2$ とすると $(K-zM)\vb d=\vb f$。質量直交規格化した列 $u_j$ の完全性 $\sum_j u_j u_j^\top M=1$ により

$$
(K-zM)^{-1}=\sum_j\frac{u_j u_j^\top}{\omega_j^2-z}.
\label{eq-am-mode-resolvent}
$$

実際、右辺に左から $K-zM$ を掛けると $\sum_j M u_j u_j^\top=1$ になる。右辺は固有振動数を除いて成立する。重根では、その固有空間に属する項をまとめれば基底の取り方に依存しない。励振方向 $\vb f$ との内積または測定方向との内積が零のモードは測定した成分に現れない。初期条件による自由振動はこの特解とは別に重ねる。従って実験で特解を観測するには、自由振動の寄与を分離する方法も指定する必要がある。

:::example[動吸振器を主質量に取り付ける]
変位を m、外力を N とし、質量行列と剛性行列の成分にそれぞれ kg、$\mathrm{N/m}$ を用いる。装置と入力を

$$
M=\left(\begin{matrix}2 & 0\\0 & 1\end{matrix}\right),\quad K=\left(\begin{matrix}12 & -4\\-4 & 4\end{matrix}\right),\quad \vb f=\left(\begin{matrix}1\\0\end{matrix}\right).
$$

とする。$z$ は $\mathrm{s^{-2}}$ で測った角振動数の二乗である。直接解くために左辺の動剛性とその行列式を作ると

$$
D=K-zM=\left(\begin{matrix}12 - 2 z & -4\\-4 & 4 - z\end{matrix}\right),\qquad \det D=2 \left(z - 8\right) \left(z - 2\right).
$$

余因子行列に $\vb f$ を掛けて $\det D$ で割り、成分ごとに約分すれば

$$
\vb d=\left(\begin{matrix}- \frac{z - 4}{2 \left(z - 8\right) \left(z - 2\right)}\\\frac{2}{\left(z - 8\right) \left(z - 2\right)}\end{matrix}\right),\qquad \vb a^\top\vb d=- \frac{z - 4}{2 \left(z - 8\right) \left(z - 2\right)},\quad \vb a=\left(\begin{matrix}1\\0\end{matrix}\right).
$$

例えば $z=4\,\mathrm{s^{-2}}$ で $\vb d=\left(\begin{matrix}0 & - \frac{1}{4}\end{matrix}\right)^\top\,\mathrm m$ となる。変位が外力に対して負の成分をもつ場合、その成分は逆位相で振れる。検算には $K\vb d$ だけを外力と比較せず、慣性項 $-zM\vb d$ も含める。主質量が静止しても付加質量は動く。そのばねの張力が外力を打ち消すので、外力の作用点で応答が零になる。
:::

:::example[反対向きに押して重心の応答を測る]
変位を m、外力を N とし、質量行列と剛性行列の成分にそれぞれ kg、$\mathrm{N/m}$ を用いる。装置と入力を

$$
M=\left(\begin{matrix}1 & 0\\0 & 1\end{matrix}\right),\quad K=\left(\begin{matrix}5 & -2\\-2 & 5\end{matrix}\right),\quad \vb f=\left(\begin{matrix}1\\-1\end{matrix}\right).
$$

とする。$z$ は $\mathrm{s^{-2}}$ で測った角振動数の二乗である。直接解くために左辺の動剛性とその行列式を作ると

$$
D=K-zM=\left(\begin{matrix}5 - z & -2\\-2 & 5 - z\end{matrix}\right),\qquad \det D=\left(z - 7\right) \left(z - 3\right).
$$

余因子行列に $\vb f$ を掛けて $\det D$ で割り、成分ごとに約分すれば

$$
\vb d=\left(\begin{matrix}- \frac{1}{z - 7}\\\frac{1}{z - 7}\end{matrix}\right),\qquad \vb a^\top\vb d=0,\quad \vb a=\left(\begin{matrix}1\\1\end{matrix}\right).
$$

例えば $z=1\,\mathrm{s^{-2}}$ で $\vb d=\left(\begin{matrix}\frac{1}{6} & - \frac{1}{6}\end{matrix}\right)^\top\,\mathrm m$ となる。観測量は二変位の和である。反対称な特解だけが励起され、対称な観測には任意の非共鳴周波数で信号がない。
:::

:::example[片側から内部のばねを通して測定する]
変位を m、外力を N とし、質量行列と剛性行列の成分にそれぞれ kg、$\mathrm{N/m}$ を用いる。装置と入力を

$$
M=\left(\begin{matrix}1 & 0\\0 & 2\end{matrix}\right),\quad K=\left(\begin{matrix}5 & -2\\-2 & 8\end{matrix}\right),\quad \vb f=\left(\begin{matrix}1\\0\end{matrix}\right).
$$

とする。$z$ は $\mathrm{s^{-2}}$ で測った角振動数の二乗である。直接解くために左辺の動剛性とその行列式を作ると

$$
D=K-zM=\left(\begin{matrix}5 - z & -2\\-2 & 8 - 2 z\end{matrix}\right),\qquad \det D=2 \left(z - 6\right) \left(z - 3\right).
$$

余因子行列に $\vb f$ を掛けて $\det D$ で割り、成分ごとに約分すれば

$$
\vb d=\left(\begin{matrix}- \frac{z - 4}{\left(z - 6\right) \left(z - 3\right)}\\\frac{1}{\left(z - 6\right) \left(z - 3\right)}\end{matrix}\right),\qquad \vb a^\top\vb d=\frac{1}{\left(z - 6\right) \left(z - 3\right)},\quad \vb a=\left(\begin{matrix}0\\1\end{matrix}\right).
$$

例えば $z=2\,\mathrm{s^{-2}}$ で $\vb d=\left(\begin{matrix}\frac{1}{2} & \frac{1}{4}\end{matrix}\right)^\top\,\mathrm m$ となる。低周波で結合ばねは外力を遠い質点へ渡す。特解の符号まで測れば、二つの共鳴をまたぐときの位相反転も識別できる。
:::

:::example[固定端三質点の中央を押す]
変位を m、外力を N とし、質量行列と剛性行列の成分にそれぞれ kg、$\mathrm{N/m}$ を用いる。装置と入力を

$$
M=\left(\begin{matrix}1 & 0 & 0\\0 & 1 & 0\\0 & 0 & 1\end{matrix}\right),\quad K=\left(\begin{matrix}2 & -1 & 0\\-1 & 2 & -1\\0 & -1 & 2\end{matrix}\right),\quad \vb f=\left(\begin{matrix}0\\1\\0\end{matrix}\right).
$$

とする。$z$ は $\mathrm{s^{-2}}$ で測った角振動数の二乗である。直接解くために左辺の動剛性とその行列式を作ると

$$
D=K-zM=\left(\begin{matrix}2 - z & -1 & 0\\-1 & 2 - z & -1\\0 & -1 & 2 - z\end{matrix}\right),\qquad \det D=- \left(z - 2\right) \left(z^{2} - 4 z + 2\right).
$$

余因子行列に $\vb f$ を掛けて $\det D$ で割り、成分ごとに約分すれば

$$
\vb d=\left(\begin{matrix}\frac{1}{z^{2} - 4 z + 2}\\- \frac{z - 2}{z^{2} - 4 z + 2}\\\frac{1}{z^{2} - 4 z + 2}\end{matrix}\right),\qquad \vb a^\top\vb d=- \frac{z - 2}{z^{2} - 4 z + 2},\quad \vb a=\left(\begin{matrix}0\\1\\0\end{matrix}\right).
$$

例えば $z=3\,\mathrm{s^{-2}}$ で $\vb d=\left(\begin{matrix}-1 & 1 & -1\end{matrix}\right)^\top\,\mathrm m$ となる。中央に節をもつ反対称モードは中央の力を受けない。その固有振動数が全体の行列式に見えても中央の伝達関数の極にはならない。
:::

## 連成振動の例

### 同じ振動子を二つつなぐ

同じ二質点 $m$ をそれぞれ壁にばね $k$ でつなぎ、互いをばね $k_c$ で結ぶ。
伸び $x_1,x_2$ に対し

$$
V=\tfrac{k}{2}(x_1^2+x_2^2)+\tfrac{k_c}{2}(x_2-x_1)^2,\qquad
K=\begin{pmatrix}k+k_c&-k_c\\-k_c&k+k_c\end{pmatrix}.
$$

同位相 $(1,1)$ では中央ばねは伸びず $\omega_+^2=k/m$。
逆位相 $(1,-1)$ では $\omega_-^2=(k+2k_c)/m$。
$x_1(0)=d,x_2(0)=0$、初速零なら二つを等振幅で足して

$$
x_1=\tfrac d2(\cos\omega_+t+\cos\omega_-t),\qquad
x_2=\tfrac d2(\cos\omega_+t-\cos\omega_-t)
$$

となる。
$k_c\ll k$ では $\omega_-$ と $\omega_+$ が近く、速い振動の包絡線が差周波数でゆっくり変わる。
これを :term[うなり]{en="beats" yomi="うなり"} という。
$m=1\,\mathrm{kg},k=100\,\mathrm{N/m},k_c=10.5\,\mathrm{N/m}$ では $10,11\,\mathrm{s^{-1}}$、最初に $x_1$ の包絡線が零になる時刻は $\pi\,\mathrm s$ である。

### 直線三原子分子の縦振動

自由な両端の質量を $m$、中央を $m_c$ とし、等しい結合ばね $k$ で一直線につなぐ。
平衡からの軸方向変位を $(x_1,x_2,x_3)$ とすれば

$$
M=\diag(m,m_c,m),\quad
V=\tfrac k2[(x_2-x_1)^2+(x_3-x_2)^2],\quad
K=k\begin{pmatrix}1&-1&0\\-1&2&-1\\0&-1&1\end{pmatrix}.
$$

両端が逆向きの $(1,0,-1)$ なら中央への二つの力は打ち消し合う。
両端が同じ変位の $(1,h,1)$ なら、非零周波数の運動では重心が振動しないので $2m+m_ch=0$。
一様並進を加えて三つの独立な解がそろう。

| 変位比 | $\omega^2$ | 物理的な動き |
| --- | --- | --- |
| $(1,1,1)$ | $0$ | 並進し、両方の結合長が一定 |
| $(1,0,-1)$ | $k/m$ | 両方の結合が同時に伸縮する対称伸縮 |
| $(1,-2m/m_c,1)$ | $k/m+2k/m_c$ | 一方が伸びると他方が縮む逆対称伸縮 |

最後の固有値は第一行 $k(1-h)=\omega^2m$ から得られ、中央の行 $2k(h-1)=\omega^2m_ch$ も満たす。
質量規格化の分母は順に $\sqrt{2m+m_c},\sqrt{2m},\sqrt{2m+4m^2/m_c}$。
例えば $m=2.0\times10^{-26}\,\mathrm{kg},m_c=m,k=200\,\mathrm{N/m}$ なら非零角振動数は $10^{14},\sqrt3\,10^{14}\,\mathrm{s^{-1}}$ の桁である。
ここでは縦変位だけを許した。
三次元の直線分子の二方向の曲げ振動には、結合角の変化に対する復元力も必要である。

```tikz caption="等質量の直線三原子。上段は対称伸縮、下段は逆対称伸縮で、矢印は同じ瞬間の変位を表す。"
\draw (0,1.6) -- (3.6,1.6);
\draw (0,0) -- (3.6,0);
\foreach \x in {0,1.8,3.6}{
  \fill (\x,1.6) circle (2pt);
  \fill (\x,0) circle (2pt);
}
\draw[->,thick] (0,1.85) -- (-.6,1.85);
\draw[->,thick] (3.6,1.85) -- (4.2,1.85);
\node at (1.8,1.94) {$0$};
\draw[->,thick] (0,.25) -- (.45,.25);
\draw[->,thick] (1.8,.25) -- (.9,.25);
\draw[->,thick] (3.6,.25) -- (4.05,.25);
\node[below] at (1.8,1.55) {$\omega^2=k/m$};
\node[below] at (1.8,-.05) {$\omega^2=3k/m$};
```

### 三次元の直線分子に残る二つの曲げ振動

縦振動だけの三自由度模型を、空間内の三原子全体に拡張する。
平衡位置を $x=-d,0,d$ として、二つの結合の長さは $d$ とする。
全変位は九成分である。
並進を三つ、実際に配位を変える回転を二つ取り除くと、独立な分子振動は $9-3-2=4$ 個になる。
分子軸のまわりの回転は、軸上にある点核を動かさないから三番目の回転変位には数えない。

横の一方向の変位を $y_1,y_2,y_3$ とする。
二つの結合の傾きの差は一次まで

$$
\delta\beta_y=\frac{y_3-y_2}{d}-\frac{y_2-y_1}{d}
=\frac{y_1-2y_2+y_3}{d}
$$

である。
$\delta\beta_y$ は折れ曲がる角の符号付き小変位で、もう一つの方向 $z$ にも同じ式がある。
直線のまわりで等方的な曲げエネルギーを

$$
V_{\rm bend}=\frac{\kappa}{2}\big[(\delta\beta_y)^2+(\delta\beta_z)^2\big]
$$

とする。
角は無次元なので $\kappa$ の単位は J である。
各横方向の剛性は、$b=(1,-2,1)^\top$ として $K_\perp=(\kappa/d^2)bb^\top$。
階数が一なので、逆行列を作らなくても固有モードをすべて求められる。

| 一つの横方向の変位 | 動き | $\omega^2$ |
| --- | --- | --- |
| $(1,1,1)$ | 横への並進 | $0$ |
| $(-1,0,1)$ | 分子全体の小回転 | $0$ |
| $(1,-2m/m_c,1)$ | 結合角が変わる曲げ | $(\kappa/d^2)(2/m+4/m_c)$ |

表の最初の二つの変位は $b$ と内積が零なので復元力を受けない。
第三の形は $M^{-1}b$ に比例し、$K_\perp M^{-1}b=(\kappa/d^2)(b^\top M^{-1}b)b$ より表の固有値を得る。
$y,z$ について同じ曲げ振動数が一つずつあるため、二つの縦伸縮に二重縮退の曲げを加えて四つとなる。
軸方向のばねだけにして $\kappa=0$ と置けば、追加の二方向も二次近似では零になり、実在分子の曲げ周波数を予言できない。

例えば $d=1.0\times10^{-10}\,\mathrm m$、$m_c=m=2.0\times10^{-26}\,\mathrm{kg}$、
$\kappa=3.0\times10^{-20}\,\mathrm J$ なら $\omega_{\rm bend}=3.0\times10^{13}\,\mathrm{s^{-1}}$。
縦のばね定数が同じでも、角を変えるエネルギーが違えば曲げの周波数は異なる。

### 二重振り子

:eqref[eq-am-double-pendulum] で二つの質量を $m$、長さを $l$ とし、ともに下向きから角を測る。二番目の質点の水平速度は一次まで $l(\dot\theta_1+\dot\theta_2)$ だから、

$$
M=ml^2\begin{pmatrix}2&1\\1&1\end{pmatrix},\quad
K=mgl\begin{pmatrix}2&0\\0&1\end{pmatrix}.
$$

$s=\omega^2 l/g$ とすると永年方程式は $2(1-s)^2-s^2=0$、根は $2\mp\sqrt2$ である。
第二行 $-s a_1+(1-s)a_2=0$ より振幅比は低周波で $a_2/a_1=\sqrt2$、高周波で $-\sqrt2$。
質量行列にも交差項がある系の例である。

### 三原子分子の同位体置換を見分ける

同じ核配置で中央の質量だけを置換し、平衡位置と二つの $k$ がほぼ等しいままとみなす。
中央が静止する対称伸縮では $\omega_s^2=k/m$ が変わらない。
一方、逆対称伸縮では中央も動くので $\omega_a^2=k/m+2k/m_c$ が変わる。
例えば $m_c=m$ から $2m$ へ変えると $\omega_a$ の比は $\sqrt{2/3}$、$\omega_s$ の比は $1$ となる。
どの核が大きく動くかを調べれば、周波数変化の向きだけでなく変化しない線も予測できる。

ただし、この結論で対称・逆対称と呼ぶのは**結合長の変化**である。
右端と左端では外向き変位の正方向が逆なので、同時に二本が伸びると直交座標の変位成分は $( -d,0,d)$ となる。
変位列の数字が左右同じかどうかだけで名前を付けると取り違える。

### 自由鎖から固定端の長い鎖へ

同質量 $m$ の $N$ 個の質点を同じばね $k$ で直列につなぎ、さらに両端をばねで壁につなぐ。
仮想的な端 $x_0=x_{N+1}=0$ を追加すれば、全質点に同じ式

$$
m\ddot x_j=k(x_{j+1}-2x_j+x_{j-1}),\qquad 1\le j\le N
$$

が使える。
形 $a_j=\sin(j\alpha)$ を代入すると、隣り合う二つの和は $2\cos\alpha\,a_j$。
右端の条件は $(N+1)\alpha=s\pi$ だから

$$
\omega_s=2\sqrt{k/m}\sin\frac{s\pi}{2(N+1)},\qquad s=1,\ldots,N .
$$

$\alpha=0,\pi$ では全変位が零になりモードではない。
得られた $N$ 個の固有値は互いに異なるためモードは独立で、すべてがそろった。
$N=3$ の中央のモードは $(1,0,-1)$ で $\omega^2=2k/m$。
分子と似た形でも壁へのばねがあるので周波数は異なる。
間隔を $a$ とすると長波長で $\omega_s\simeq a\sqrt{k/m}\,s\pi/[(N+1)a]$ となり、[弦の振動](/physics/analytical-mechanics/field-mechanics/)につながる。

### 二種類の質量が交互に並ぶ格子

同じばね定数 $k$ の一次元鎖で、単位胞 $j$ の二つの変位を $u_j,v_j$ とする。質量をそれぞれ $m,M$、単位胞の間隔を $a$ とする。ばねの長さの変化を隣り合う変位の差で表すと

$$
m\ddot u_j=k(v_j+v_{j-1}-2u_j),\qquad
M\ddot v_j=k(u_j+u_{j+1}-2v_j).
$$

周期境界を課して $(u_j,v_j)=(U,V)\exp[i(j\kappa a-\omega t)]$ の実部を取る。$\kappa$ に $2\pi/a$ を足しても同じ位相になるので $|\kappa|\le\pi/a$ を調べればよい。隣接胞の係数は $e^{\pm i\kappa a}$。二元の連立方程式で行列式を零とおくと

$$
mM\omega^4-2k(m+M)\omega^2+4k^2\sin^2(\kappa a/2)=0 .
$$

よって二つの枝は
$$
\omega_\pm^2=k(1/m+1/M)\ \pm\
k\sqrt{(1/m+1/M)^2-4\sin^2(\kappa a/2)/(mM)}.
\label{eq-am-diatomic-chain}
$$

小波数では下の根を直接引くと小さい差になる。根の積を上の根で割れば $\omega_-^2\simeq ka^2\kappa^2/[2(m+M)]$ が得られる。長波長の極限で二質点が同方向に進む :term[音響分枝]{en="acoustic branch" yomi="おんきょうぶんし"} である。上の根では $\kappa=0$ でも有限の $\omega$ があり、$mU+MV=0$ を満たして逆向きに振動する。これを :term[光学分枝]{en="optical branch" yomi="こうがくぶんし"} と呼ぶ。名称は、異符号の電荷が逆向きに動く結晶では電場と結合できることに由来する。ここで扱うばね系に電荷を仮定する必要はない。

:::example[長波長での二枝の分離]
$m=1\,\mathrm{kg},M=3\,\mathrm{kg},k=12\,\mathrm{N/m},a=0.20\,\mathrm m$ とし、$\kappa a=0$ を選ぶ。左右のばねを足す際、隣の胞の位相を入れると、変位振幅に作用する剛性行列と質量行列は

$$
K(\kappa)=\left(\begin{matrix}24 & -24\\-24 & 24\end{matrix}\right),\qquad M_{\rm cell}=\left(\begin{matrix}1 & 0\\0 & 3\end{matrix}\right).
$$

$z=\omega^2$ とおいた特性式は $3 z \left(z - 32\right)=0$。各根を元の二つの力の式に戻して零空間を求めると

$$
\begin{aligned}\omega^2=0,\quad (U,V)^\top\ \propto\ \left(\begin{matrix}1\\1\end{matrix}\right)\\ \omega^2=32,\quad (U,V)^\top\ \propto\ \left(\begin{matrix}-3\\1\end{matrix}\right).\end{aligned}
$$

となる。周波数の二乗の単位は $\mathrm{s^{-2}}$ である。振動するモードでは、$U,V$ に共通の位相因子 $e^{i\delta}$ を掛けると時間原点が $\delta/\omega$ だけ移る。したがって振幅比は同じである。下の根は全体の並進で、復元力は零である。上の根は一つの胞の重心を動かさないため、質量が大きい方の振幅が小さい。
:::

:::example[胞の間で位相が四分の一周ずれる波]
$m=1\,\mathrm{kg},M=3\,\mathrm{kg},k=12\,\mathrm{N/m},a=0.20\,\mathrm m$ とし、$\kappa a=\frac{\pi}{2}$ を選ぶ。左右のばねを足す際、隣の胞の位相を入れると、変位振幅に作用する剛性行列と質量行列は

$$
K(\kappa)=\left(\begin{matrix}24 & -12 + 12 i\\-12 - 12 i & 24\end{matrix}\right),\qquad M_{\rm cell}=\left(\begin{matrix}1 & 0\\0 & 3\end{matrix}\right).
$$

$z=\omega^2$ とおいた特性式は $3 \left(z^{2} - 32 z + 96\right)=0$。各根を元の二つの力の式に戻して零空間を求めると

$$
\begin{aligned}\omega^2=- 4 \left(-4 + \sqrt{10}\right),\quad (U,V)^\top\ \propto\ \left(\begin{matrix}-1 + \frac{\sqrt{10}}{2} - \frac{\sqrt{10} i}{2} + i\\1\end{matrix}\right)\\ \omega^2=4 \left(\sqrt{10} + 4\right),\quad (U,V)^\top\ \propto\ \left(\begin{matrix}- \frac{\sqrt{10}}{2} - 1 + i + \frac{\sqrt{10} i}{2}\\1\end{matrix}\right).\end{aligned}
$$

となる。周波数の二乗の単位は $\mathrm{s^{-2}}$ である。二つの副格子の複素振幅比には位相差がある。実部だけを最初から未知量にするとこの解を落とす。
:::

:::example[Brillouin ゾーンの端]
$m=1\,\mathrm{kg},M=3\,\mathrm{kg},k=12\,\mathrm{N/m},a=0.20\,\mathrm m$ とし、$\kappa a=\pi$ を選ぶ。左右のばねを足す際、隣の胞の位相を入れると、変位振幅に作用する剛性行列と質量行列は

$$
K(\kappa)=\left(\begin{matrix}24 & 0\\0 & 24\end{matrix}\right),\qquad M_{\rm cell}=\left(\begin{matrix}1 & 0\\0 & 3\end{matrix}\right).
$$

$z=\omega^2$ とおいた特性式は $3 \left(z - 24\right) \left(z - 8\right)=0$。各根を元の二つの力の式に戻して零空間を求めると

$$
\begin{aligned}\omega^2=8,\quad (U,V)^\top\ \propto\ \left(\begin{matrix}0\\1\end{matrix}\right)\\ \omega^2=24,\quad (U,V)^\top\ \propto\ \left(\begin{matrix}1\\0\end{matrix}\right).\end{aligned}
$$

となる。周波数の二乗の単位は $\mathrm{s^{-2}}$ である。隣の胞の符号が逆になって両側の力が相殺し、一方の副格子だけが振れる。二つの端の周波数の差が禁制帯を作る。
:::

:::example[同じ質量を二個ずつ数えた場合]
$m=2\,\mathrm{kg},M=2\,\mathrm{kg},k=12\,\mathrm{N/m},a=0.20\,\mathrm m$ とし、$\kappa a=\pi$ を選ぶ。左右のばねを足す際、隣の胞の位相を入れると、変位振幅に作用する剛性行列と質量行列は

$$
K(\kappa)=\left(\begin{matrix}24 & 0\\0 & 24\end{matrix}\right),\qquad M_{\rm cell}=\left(\begin{matrix}2 & 0\\0 & 2\end{matrix}\right).
$$

$z=\omega^2$ とおいた特性式は $4 \left(z - 12\right)^{2}=0$。各根を元の二つの力の式に戻して零空間を求めると

$$
\begin{aligned}\omega^2=12,\quad (U,V)^\top\ \propto\ \left(\begin{matrix}1\\0\end{matrix}\right)\\ \omega^2=12,\quad (U,V)^\top\ \propto\ \left(\begin{matrix}0\\1\end{matrix}\right).\end{aligned}
$$

となる。周波数の二乗の単位は $\mathrm{s^{-2}}$ である。縮退は単位胞を物理的な最小周期の二倍に選んだ結果である。質量差を戻せば境界で縮退が解ける。
:::

単位胞を二つ並べると、波数零という条件の意味が明瞭になる。図のどちらの段でも隣の胞は同じ動きを繰り返す。ただし同一の胞の内部では、光学分枝の二粒子が異なる向きに動く。重い方を軽い方の三倍の質量とすると、その変位の大きさは三分の一である。図の矢印はある瞬間の変位の比であり、波が伝わる速さを表してはいない。

```tikz caption="波数零の二原子鎖。上段が音響分枝の零モード、下段が光学分枝で、丸の大きさは質量の違いを表す。"
\foreach \y in {0,1.6}{
 \draw (0,\y)--(5.4,\y);
 \foreach \x in {0.6,3.2}{\filldraw[fill=white] (\x,\y) circle (.12);}
 \foreach \x in {1.9,4.5}{\filldraw[fill=white] (\x,\y) circle (.23);}
}
\foreach \x in {0.6,1.9,3.2,4.5}{\draw[->,thick] (\x,1.98)--++(.53,0);}
\foreach \x in {0.6,3.2}{\draw[->,thick] (\x,.38)--++(.75,0);}
\foreach \x in {1.9,4.5}{\draw[->,thick] (\x,.38)--++(-.25,0);}
\node[left] at (0,1.6) {$\omega_-$};
\node[left] at (0,0) {$\omega_+$};
```

## 一般固有値問題としての定式化

### Rayleigh 商と周波数の上下限

試しに一つの変位形 $v$ を指定し、その形を保って動かすときの角振動数の二乗は

$$
\mathcal R(v)=\frac{v^\top Kv}{v^\top Mv}.
$$

これを :term[Rayleigh 商]{en="Rayleigh quotient" yomi="れいりーしょう"} という。
$v=\sum_s c_sa_s$ と展開すれば $\mathcal R=\sum_s\lambda_sc_s^2/\sum_sc_s^2$ である。
従って最小固有値と最大固有値の間にあり、最小モードの試算値は真の最小値の上限になる。
既知の低次モードと質量直交する試行空間に絞れば、残りの最小値に同じ議論が使える。

:::example[近似形の精度]
先の二重振り子で $v=(1,1)$ を試すと $\mathcal R=(3/5)g/l$ である。
正確な最低値は $(2-\sqrt2)g/l\simeq0.586g/l$。
試行形は両腕が同じ角で動くという制限を加えており、その分だけ固有値を高く見積もる。
:::

### 二番目以降の周波数も上下から挟む

最低モードより上の周波数を求めるとき、すべての低次モードを厳密に知っているとは限らない。
その場合は試行ベクトル一つの代わりに、一緒に試せる変位の空間を使う。
$\lambda_1\le\cdots\le\lambda_n$ と番号を振ると

$$
\lambda_j=\min_{\dim S=j}\ \max_{v\in S,\ v\ne0}\mathcal R(v)
$$

となる。
証明には質量正規直交基底を使えば足りる。
任意の $j$ 次元部分空間 $S$ の中で $a_1,\ldots,a_{j-1}$ に質量直交する条件は $j-1$ 個なので、これを満たす非零の $v$ がある。
その $v$ には第 $j$ 以降の成分しかなく $\mathcal R(v)\ge\lambda_j$。
従って右辺は $\lambda_j$ 以上である。
一方 $S=\spn(a_1,\ldots,a_j)$ を選べばその最大値は $\lambda_j$。
二つの不等式を合わせて式を得る。

測定や近似計算で少数の変位しか許さない場合、試行基底を列に持つ行列 $B$ に対して $\eta=Bz$ と置く。
こうすると小さい行列 $B^\top KB,B^\top MB$ を使って同じ一般固有値問題が解ける。
上述の最小化の範囲を狭くしたので、試行空間の第 $j$ 固有値は真の $\lambda_j$ 以上である。
基底を一つ追加して空間を広げれば上限は増えない。
これが Rayleigh–Ritz 法の変位法である。

:::example[三原子の試行形と零モード]
自由な等質量三原子で $(1,0,0)$ を試すと Rayleigh 商は $k/m$。
最低振動数が $k/m$ だと読んでしまうと一様並進の零モードを失う。
重心を固定するなら試行形も $x_1+x_2+x_3=0$ に制限する必要がある。
この平面で $(1,-1,0)$ を試すと

$$
\mathcal R=\frac{k[(1+1)^2+(0+1)^2]}{m(1+1)}
=\frac{5k}{2m}.
$$

最低の正固有値 $k/m$ に対する上限である。
さらに $(0,1,-1)$ を加えると重心固定平面全体を張り、その中の二固有値 $k/m,3k/m$ が正確に戻る。
:::

### 行列式とモードの確認

大きな系では $\det(K-\lambda M)$ は次数 $n$ の多項式になる。
根を得た後、各ベクトルについて残差 $Ka_s-\lambda_sMa_s$、全体について $A^\top MA$ を確認する。
異なる固有値に属するモードの質量直交性は、対称性から

$$
a_t^\top Ka_s=\lambda_s a_t^\top Ma_s
=\lambda_t a_t^\top Ma_s
$$

としても示せる。
縮退時に固有ベクトルの本数が足りる理由は :ref[thm-am-normal-modes] にある。
一般行列の対角化可能性の判定は :ref[thm-linalg-diagonalization-criterion] を参照されたい。

例えば $M=\diag(1,2)$、$K=\begin{pmatrix}2&-1\\-1&3\end{pmatrix}$ を SI の数値で与えると、行列式は $2\lambda^2-7\lambda+5$、
根は $1,5/2$。
固有ベクトル $(1,1),(2,-1)$ の質量内積は $2-2=0$ であり、質量ノルムの二乗は $3,6$ である。

### 一つのばねや質量を変える感度

固有値 $\lambda_s$ が単純で、行列 $K,M$ が制御パラメータ $h$ で滑らかに変化するとする。
$Ka_s=\lambda_sMa_s$ を微分し、左から $a_s^\top$ を掛ける。
$a_s^\top K=\lambda_s a_s^\top M$ を使うと $a_s'$ を含む両項が消え、質量規格化のもとで

$$
\lambda_s'=a_s^\top(K'-\lambda_s M')a_s
$$

を得る。
例えば質量行列が対角で一つの質量 $m_j$ を増やし、ばねを変えなければ $\partial\lambda_s/\partial m_j=-\lambda_s(a_s)_j^2$。
節に当たる質点では一次変化が零になる。
縮退点ではどのモードを追うかが一意にならないので、摂動を縮退部分空間に制限して対角化する。

:::example[ばね一本が測る変位差]
二振動子の結合ばねの定数 $k_c$ を変える。
$K'=\begin{pmatrix}1&-1\\-1&1\end{pmatrix}$ なので補正は $[(a_s)_1-(a_s)_2]^2$。
同位相では零、逆位相では質量規格化した $(1,-1)/\sqrt{2m}$ に対して $2/m$ である。
厳密式 $(k+2k_c)/m$ を微分しても同じ値が得られる。
:::

### 調和外力と共鳴

小さな力 $f\cos\nu t$ を加え、減衰を無視する。
一つずつ基準方程式に投影すると

$$
\ddot Q_s+\omega_s^2Q_s=f_s\cos\nu t,\qquad f_s=a_s^\top f .
$$

初期変位と初速がともに零で $\nu\ne\omega_s$ なら

$$
Q_s(t)=\frac{f_s}{\omega_s^2-\nu^2}(\cos\nu t-\cos\omega_st).
$$

$\nu=\omega_s>0$ では極限をとって $Q_s=f_st\sin\omega_st/(2\omega_s)$ となる。零モードに一定力が作用する $\nu=\omega_s=0$ の場合は、直接 $\ddot Q_s=f_s$ を積分して $Q_s=f_st^2/2$ を得る。
長時間には振幅が成長するため、線形化が有効な時間もそれに応じて限られる。
外力の周波数が一致しても $f_s=0$ ならそのモードは励起されない。
例えば二質点を同じ方向に同じ力で駆動すれば逆位相モードにエネルギーは入らない。

:::example[弱い散逸とモードの独立性]
抵抗の行列を $D$ として $M\ddot\eta+D\dot\eta+K\eta=0$ とする。
$D=\gamma M$ なら $A^\top DA=\gamma I$ なので各モードに同じ $\gamma\dot Q_s$ が加わる。
$\gamma<2\omega_s$ での時間依存は $e^{-\gamma t/2}\cos(\sqrt{\omega_s^2-\gamma^2/4}\,t+\delta)$ である。
一方、二質点の片方だけにダンパーを付けると、同位相・逆位相座標に散逸の交差項が現れる。
保存系でのモードがそのまま減衰系の独立なモードになるとは限らない。
:::

### 減衰があるときの振幅と吸収仕事率

質量で規格化した一つのモードに対し、$\ddot Q+\gamma\dot Q+\omega_0^2Q=f\cos\nu t$、$\gamma>0$、$\omega_0>0$ とする。過渡解は時間とともに減衰する。残る周期解を $Q=\Re(Ae^{i\nu t})$ と置くと

$$
A=\frac{f}{\omega_0^2-\nu^2+i\gamma\nu},\qquad
|A|^2=\frac{f^2}{(\omega_0^2-\nu^2)^2+\gamma^2\nu^2}.
$$

一定時間に平均した入力仕事率は $\langle f\cos\nu t\,\dot Q\rangle=(f/2)\Re(i\nu A)$ である。$A$ を代入すると $\gamma\nu^2|A|^2/2$ となり、$\gamma\langle\dot Q^2\rangle$ に等しい。平均のエネルギー収支も、$\dd(\dot Q^2/2+\omega_0^2Q^2/2)/\dd t=f\cos\nu t\,\dot Q-\gamma\dot Q^2$ を一周期積分すれば確かめられる。

$z=\nu^2$ とすれば、変位振幅の分母は $(\omega_0^2-z)^2+\gamma^2 z$ である。その微分は $2(z-\omega_0^2)+\gamma^2$ だから、最大振幅は $z=\omega_0^2-\gamma^2/2>0$ で生じる。この根が正でなければ最大は静的な入力で生じる。一方、平均仕事率の周波数依存部分は $z/[(\omega_0^2-z)^2+\gamma^2z]$ であり、その微分の符号は $\omega_0^4-z^2$ の符号に一致する。従って平均仕事率は $\nu=\omega_0$ で最大となる。例えば $\omega_0=5\,\mathrm{s^{-1}}$、$\gamma=2\,\mathrm{s^{-1}}$ なら、変位が最大になる周波数は $\sqrt{23}\,\mathrm{s^{-1}}$ であり、仕事率の最大は $5\,\mathrm{s^{-1}}$ にある。

### 実験で観測する振幅とモードエネルギー

基準座標の規格化によって力の加わり方と観測の仕方を一つの式で比較できる。
安定な非縮退の系に角振動数 $\nu$ の力の複素振幅 $f$ を加え、共鳴から離れているとする。
$\eta(t)=\Re(\hat\eta e^{i\nu t})$ として強制振動だけを取り出すと

$$
(K-\nu^2M)\hat\eta=f .
$$

$A^\top M A=I$ を使えば、解は

$$
\hat\eta=\sum_s a_s\frac{a_s^\top f}{\omega_s^2-\nu^2}.
$$

実際、右辺に $K-\nu^2M$ を掛けると $\sum_s Ma_sa_s^\top f=f$。
最後の等式は $AA^\top=M^{-1}$ から従う。
測定器が変位の線形結合 $b^\top\eta$ を読むなら、第 $s$ の共鳴への寄与には $(b^\top a_s)(a_s^\top f)$ が掛かる。
励起する位置と観測する位置のどちらに節があっても、その測定での共鳴ピークは消える。

:::example[二振動子の反共鳴]
同じ二振動子のうち左だけに力 $F\cos\nu t$ を加える。
$D=k+k_c-m\nu^2$ とおくと、特解の振幅は

$$
\begin{pmatrix}D&-k_c\\-k_c&D\end{pmatrix}
\begin{pmatrix}\hat x_1\\\hat x_2\end{pmatrix}
=\begin{pmatrix}F\\0\end{pmatrix}.
$$

行列を逆にすると

$$
\hat x_1=\frac{DF}{D^2-k_c^2},\qquad
\hat x_2=\frac{k_cF}{D^2-k_c^2}.
$$

$k_c>0$ で $\nu^2=(k+k_c)/m$ のとき、分母は $-k_c^2$ で有限だが、加振点の変位は零になる。
隣の質点には $\hat x_2=-F/k_c$ の運動が残り、そのばね力が外力を打ち消す。
この周波数は二つの共鳴の間にある。
無減衰系では初期条件で自由振動も重なるので、ここで零といったのは駆動周波数成分の振幅である。
:::

### 運動の近似をどの範囲まで信頼するか

線形振動子を長時間追うとき、瞬間の力の誤差が小さいことと位相が正しいことは別である。
単振り子を小振幅 $A$ で放すと、エネルギーから厳密な四分の一周期は

$$
\frac T4=\sqrt{\frac{l}{2g}}\int_0^A\frac{\dd\theta}{\sqrt{\cos\theta-\cos A}}.
$$

$\sin(\theta/2)=\sin(A/2)\sin s$ と置けば

$$
\frac T4=\sqrt{\frac lg}\int_0^{\pi/2}\frac{\dd s}
{\sqrt{1-\sin^2(A/2)\sin^2s}} .
$$

分母を二次まで展開し $\int_0^{\pi/2}\sin^2s\,\dd s=\pi/4$ を使うと $T=2\pi\sqrt{l/g}[1+A^2/16+O(A^4)]$。
例えば $A=0.20$ rad なら一周期の相対誤差は約 $0.0025$ だが、百周期では約四分の一周期のずれに積み上がる。
精密に周期を測るときには、線形化の誤差を何周期にわたって評価するかも指定する。

### 初期エネルギーから実験で現れるモードを予測する

測定器が周波数だけを読んだ場合、系がもつ固有振動数の一部しか見えていないことがある。
自由振動で $Q_s(t)=C_s\cos\omega_st+D_s\sin\omega_st$ と書けば、第 $s$ モードのエネルギーは

$$
E_s=\frac{\omega_s^2}{2}(C_s^2+D_s^2)
=\frac12[(a_s^\top M\dot\eta_0)^2+
\omega_s^2(a_s^\top M\eta_0)^2].
$$

$\omega_s>0$ では初期変位と初速の両方がそのモードに対して質量直交するときに限り $E_s=0$ となる。零モードでは $E_s=\dot Q_s^2/2$ であり、静止した変位があってもエネルギーは零である。
一点で変位の振幅が零でも、その一点がモードの節であれば系全体の $E_s$ は零とは限らない。

:::example[二つのエネルギーの割合]
等質量二振動子で一方だけを $d$ ずらして放す。
正規化した二固有ベクトルは $(1,\pm1)/\sqrt{2m}$ なので、

$$
Q_+(0)=Q_-(0)=\sqrt{m/2}\,d .
$$

二つの初期振幅は等しいがエネルギーは

$$
E_+=\frac{k d^2}{4},\qquad E_-=\frac{(k+2k_c)d^2}{4}.
$$

和は $(k+k_c)d^2/2$ で、実際に変形させた左のばねと中央ばねの初期エネルギーの和に一致する。
$k_c=k/2$ なら高周波モードが全エネルギーの $2/3$ を担う。
自由振動の間に左右の質点の揺れの大きさは変化しても、この二つのモードエネルギーは個別に保存する。
:::

### 速度に比例する力があるときの固有値問題

回転系や磁場では速度の一次に比例する力が仕事をせず、散逸を表す力とは異なる形で振動数を変える。平衡点の近くで $M\ddot q+(D+G)\dot q+Kq=0$ とし、$M,K,D$ は対称、$G^\top=-G$ とする。変位速度を左から掛けると $\dot E=-\dot q^\top D\dot q$。実際 $\dot q^\top G\dot q$ は転置で自分自身の負に等しいので零である。$D$ が半正定値ならエネルギーは増えないが、$D$ と $K$ が同じモードで対角化されるとは限らない。解を $v e^{st}$ とすると

$$
[s^2M+s(D+G)+K]v=0 .
$$

これは行列係数をもつ二次の固有値問題である。保存系で用いた $Kv=\omega^2Mv$ は速度項が零の場合に対応する。特に $D=\alpha M+\beta K,\ G=0$ なら質量座標で各モードに $s^2+(\alpha+\beta\omega_j^2)s+\omega_j^2=0$ を得る。この :term[比例減衰]{en="proportional damping" yomi="ひれいげんすい"} ではモードの形を変えずに減衰率を別々に計算できる。一つのダンパーを任意の位置につけた系では、一般にこの単純化はできない。

:::example[磁場で分裂する等方な振動]
両軸の質量を $1\,\mathrm{kg}$ とし、$K=\left(\begin{matrix}4 & 0\\0 & 4\end{matrix}\right)\,\mathrm{N/m}$、$D+G=\left(\begin{matrix}0 & -3\\3 & 0\end{matrix}\right)\,\mathrm{kg/s}$ とする。初期位置における外力は釣り合うものとする。$e^{st}$ の係数とその行列式は

$$
P(s)=\left(\begin{matrix}s^{2} + 4 & - 3 s\\3 s & s^{2} + 4\end{matrix}\right),\qquad \det P(s)=\left(s^{2} + 1\right) \left(s^{2} + 16\right).
$$

二つの実方程式を第一式に第二式の $i$ 倍を加えてまとめる。$x+iy$ の特性式は $\left(s - i\right) \left(s + 4 i\right)=0$、$x-iy$ ではその共役になる。従って $Z=x+iy=Ae^{it}+Be^{-4it}$（$t$ は秒で測った数値）である。$x=\Re Z$、$y=\Im Z$ として実変位が得られる。二つの円偏光は角振動数 $1,4\,\mathrm{s^{-1}}$ で逆向きに回る。例えば初速度 $(1,2)\,\mathrm{m/s}$ の瞬間にも $\dot E=-\dot q^\top D\dot q=0\,\mathrm W$ である。単に剛性行列を対角化すると、同じ角振動数のままという誤った予測になる。
:::

:::example[復元力が消えた磁場中の粒子]
両軸の質量を $1\,\mathrm{kg}$ とし、$K=\left(\begin{matrix}0 & 0\\0 & 0\end{matrix}\right)\,\mathrm{N/m}$、$D+G=\left(\begin{matrix}0 & -2\\2 & 0\end{matrix}\right)\,\mathrm{kg/s}$ とする。初期位置における外力は釣り合うものとする。$e^{st}$ の係数とその行列式は

$$
P(s)=\left(\begin{matrix}s^{2} & - 2 s\\2 s & s^{2}\end{matrix}\right),\qquad \det P(s)=s^{2} \left(s^{2} + 4\right).
$$

二つの実方程式を第一式に第二式の $i$ 倍を加えてまとめる。$x+iy$ の特性式は $s \left(s + 2 i\right)=0$、$x-iy$ ではその共役になる。従って $x+iy=A+Be^{-2it}$（$t$ は秒で測った数値）であり、時間指数は $0,-2i\,\mathrm{s^{-1}}$ である。任意の初期位置と速度に対して $A,B$ が決まる。例えば初速度 $(1,2)\,\mathrm{m/s}$ の瞬間には $\dot E=-\dot q^\top D\dot q=0\,\mathrm W$。回転する解のほかに、任意の位置に静止する解がある。零根を有限周期の振動として数えることはできない。
:::

:::example[臨界減衰に調整した二方向の台]
両軸の質量を $1\,\mathrm{kg}$ とし、$K=\left(\begin{matrix}4 & 0\\0 & 4\end{matrix}\right)\,\mathrm{N/m}$、$D+G=\left(\begin{matrix}4 & 0\\0 & 4\end{matrix}\right)\,\mathrm{kg/s}$ とする。初期位置における外力は釣り合うものとする。$e^{st}$ の係数とその行列式は

$$
P(s)=\left(\begin{matrix}s^{2} + 4 s + 4 & 0\\0 & s^{2} + 4 s + 4\end{matrix}\right),\qquad \det P(s)=\left(s + 2\right)^{4}.
$$

二つの実方程式を第一式に第二式の $i$ 倍を加えてまとめる。$x+iy$ の特性式は $\left(s + 2\right)^{2}=0$、$x-iy$ ではその共役になる。従って各方向の変位は $(A+Bt)e^{-2t}$（$t$ は秒で測った数値）である。二方向それぞれで初期位置と速度を指定して $A,B$ を決める。例えば初速度 $(1,2)\,\mathrm{m/s}$ の瞬間には $\dot E=-\dot q^\top D\dot q=-20\,\mathrm W$。二方向で同じ臨界減衰が起きる。重根では指数関数だけでは初期条件をすべて表せず、時間を掛けた解が必要になる。
:::

:::example[粘性と横向きの力が共存する台]
両軸の質量を $1\,\mathrm{kg}$ とし、$K=\left(\begin{matrix}2 & 0\\0 & 2\end{matrix}\right)\,\mathrm{N/m}$、$D+G=\left(\begin{matrix}1 & -1\\1 & 1\end{matrix}\right)\,\mathrm{kg/s}$ とする。初期位置における外力は釣り合うものとする。$e^{st}$ の係数とその行列式は

$$
P(s)=\left(\begin{matrix}s^{2} + s + 2 & - s\\s & s^{2} + s + 2\end{matrix}\right),\qquad \det P(s)=s^{4} + 2 s^{3} + 6 s^{2} + 4 s + 4.
$$

二つの実方程式を第一式に第二式の $i$ 倍を加えてまとめる。$x+iy$ の特性式は $s^{2} + s \left(1 + i\right) + 2=0$、$x-iy$ ではその共役になる。従って $x+iy$ の時間指数は $[-(1+i)\pm\sqrt{-8+2i}]/2$（$\mathrm{s^{-1}}$）である。平方根の実部の絶対値は $\sqrt{\sqrt{17}-4}<1$ だから両方とも減衰する。初期条件から二つの指数関数の係数を求め、その和の実部と虚部を取ればよい。例えば初速度 $(1,2)\,\mathrm{m/s}$ の瞬間には $\dot E=-\dot q^\top D\dot q=-5\,\mathrm W$。エネルギーを失う速さは粘性だけが決めるが、粒子の軌跡は横向きの力によって曲がる。
:::

### ばね定数を周期的に変えるパラメータ励振

外力を加えずに支持点を揺らしても振幅が増すことがある。平衡位置の近くで $\ddot x+\omega_0^2[1+h\cos(2\omega t)]x=0$ とする。$|h|\ll1$、$\omega-\omega_0$ は $h\omega_0$ の桁とする。位相を外部の半周波数に合わせて $x=a(t)\cos\omega t+b(t)\sin\omega t$ と置く。$a,b$ の一周期あたりの相対変化は小さいと仮定し、二階の時間微分を一次の近似で落とす。$\cos2s\cos s=(\cos3s+\cos s)/2,\ \cos2s\sin s=(\sin3s-\sin s)/2$ だから、基本調波の係数を零とおくと
$$
\dot a=\frac{\omega_0^2-\omega^2-h\omega_0^2/2}{2\omega}\,b,\qquad
\dot b=-\frac{\omega_0^2-\omega^2+h\omega_0^2/2}{2\omega}\,a .
$$
三倍周波数の強制項は固有周波数から離れており、振幅 $O(h)$ の有界な波形補正を与える。この補正を再び基本調波へ戻す項は二次なので、一次の増大率を変えない。上の定数行列の固有値を $\pm s$ とすると
$$
s^2=\frac{(h\omega_0^2/2)^2-(\omega_0^2-\omega^2)^2}{4\omega^2}.
$$
これは基本の不安定帯の一次近似である。中心 $\omega=\omega_0$ では $s=|h|\omega_0/4$、近くでの半幅は $|\omega-\omega_0|<|h|\omega_0/4$ となる。励振角振動数は $2\omega$ なので、実際のモーターの周波数をそのまま $\omega$ に入れてはいけない。厳密に静止した初期値 $x=\dot x=0$ は解のままで、増大は初期のずれのうち増大固有方向の成分に起こる。

例えば $\omega_0=10\,\mathrm{s^{-1}},h=0.04$ とし、剛性を $20\,\mathrm{s^{-1}}$ で変える。振幅の $e$ 倍時間は約 $10\,\mathrm s$。弱い粘性による項 $2\gamma\dot x$ も加えれば、基本調波での二式に $-\gamma a,-\gamma b$ が加わる。従って中心で $\gamma=0.06\,\mathrm{s^{-1}}$ なら正味の増大率は約 $0.04\,\mathrm{s^{-1}}$、$\gamma=0.12\,\mathrm{s^{-1}}$ なら二枝とも減衰する。瞬間エネルギー $m\dot x^2/2+k(t)x^2/2$ の変化は $\dot kx^2/2-2m\gamma\dot x^2$。エネルギー供給の符号が変位の位相に依存するため、同じ剛性変化でも増大する位相と減衰する位相がある。

### 非線形な復元力で共鳴周波数が移る

線形の共鳴解が成長すると、復元力の高次項が駆動周波数と自由振動数の差に匹敵することがある。弱い減衰と三次の復元力をもつ振動子
$$
\ddot x+\gamma\dot x+\omega_0^2x+\alpha x^3=F\cos\nu t
$$
を考える。$\gamma/\omega_0$、$|\nu-\omega_0|/\omega_0$、$|\alpha|A^2/\omega_0^2$ を同じ小さい次数とし、$|F|/(\omega_0^2A)$ もその次数とする。$x=A\cos(\nu t-\delta)$ を先頭の近似とすると、$\cos^3 s=(3\cos s+\cos 3s)/4$ から基本調波の復元力に $3\alpha A^3/4$ が加わる。三倍調波は共鳴から離れており、$x$ に対して相対的に一次の波形補正を与える。その補正が基本調波に戻るのは二次以降である。基本調波の余弦と正弦をそれぞれ比較して
$$
F\cos\delta=A(\omega_0^2-\nu^2+3\alpha A^2/4),\qquad
F\sin\delta=\gamma\nu A
$$
を得る。二乗和を取れば、振幅は $F^2=A^2[(\omega_0^2-\nu^2+3\alpha A^2/4)^2+\gamma^2\nu^2]$ を満たす。$\alpha>0$ なら振幅が大きくなるほど、同じ位相遅れを得るための周波数は高くなる。$\alpha<0$ なら逆になる。$A^2$ に関するこの三次式には複数の正根があり得る。どの周期応答が安定かを知るには、初期条件を変えたときの振幅と位相の時間発展も調べる必要がある。

駆動と減衰を零にした場合は $\omega^2=\omega_0^2+3\alpha A^2/4+O(A^4)$ である。角振動数の補正は $3\alpha A^2/(8\omega_0)$。例えば $\omega_0=10\,\mathrm{s^{-1}}$、$\alpha A^2=4\,\mathrm{s^{-2}}$ なら、硬化による一次の増分は $0.15\,\mathrm{s^{-1}}$ となる。

### 速い駆動が平衡点を変える

高周波の力は平均が零でも、変位と力の相関を通じてゆっくりした運動を変える。$m\ddot x=-V_0'(x)+f(x)\cos\Omega t$ を考え、$x=X+\xi$ とおく。$X$ が駆動一周期の間ほぼ一定なら、零平均の速い成分は $m\ddot\xi=f(X)\cos\Omega t$ から $\xi=-f(X)\cos\Omega t/(m\Omega^2)$ と近似できる。速い変位が $f$ の変化尺度より小さければ $f(X+\xi)=f(X)+f'(X)\xi+\cdots$ と展開できる。平均では後者の寄与が $-f(X)f'(X)/(2m\Omega^2)$ となる。従って主要な平均運動を $V_{\rm eff}(X)=V_0(X)+f(X)^2/(4m\Omega^2)$ で記述できる。

長さ $l$ の振り子の支点を鉛直に $a\cos\Omega t$ だけ動かす。下向きからの角を $\theta$ とすると、支持点の加速度を重力に加えた式は $ml^2\ddot\theta=-ml(g-a\Omega^2\cos\Omega t)\sin\theta$ である。上の式の質量を $ml^2$、周期力の振幅を $mla\Omega^2\sin\theta$ に置き換えれば
$$
V_{\rm eff}=-mgl\cos\theta+\frac{ma^2\Omega^2}{4}\sin^2\theta
$$
となる。$a/l\ll1,\ \Omega\gg\sqrt{g/l}$ の範囲で、$a^2\Omega^2>2gl$ なら上向き $\theta=\pi$ での平均ポテンシャルの曲率 $m(a^2\Omega^2/2-gl)$ は正になる。例えば $l=0.20\,\mathrm m$、$g=9.8\,\mathrm{m/s^2}$、$a=0.01\,\mathrm m$、$\Omega=300\,\mathrm{s^{-1}}$ なら $a^2\Omega^2=9\,\mathrm{m^2/s^2}$ で、しきい値 $2gl=3.92\,\mathrm{m^2/s^2}$ を上回る。しきい値のごく近くを判定するには高次の平均化も要る。

### 試行モードの誤差を力の残差で測る

近似した変位パターン $v$ を与えられたとき、Rayleigh 商がもっともらしい値でも固有応答が再現されるとは限らない。$v^\top Mv=1$ と規格化し、$\lambda=v^\top Kv$、残った力を $r=Kv-\lambda Mv$ とする。質量直交規格化した厳密な基底に $v=\sum_j a_j u_j$ と展開すれば、$\sum_j a_j^2=1$ であり
$$
r^\top M^{-1}r=\sum_j(\lambda_j-\lambda)^2a_j^2 .
$$
従って少なくとも一つの固有値は $\lambda$ から $\sqrt{r^\top M^{-1}r}$ 以内にある。また狙った固有空間以外の固有値がすべて $\lambda$ から距離 $\delta>0$ 以上なら、その空間の外の重みの和は $r^\top M^{-1}r/\delta^2$ 以下である。近い二つの固有値を分離するには、その間隔に比べて十分に小さな残差を要する。

:::example[検出器で読んだ変位比を点検する]
質量がどちらも $1\,\mathrm{kg}$、剛性行列が対角で $4,9\,\mathrm{N/m}$ の二つの振動子を使う。測定方向を規格化して $v=(\sqrt{99},1)/10$ とすれば、第一軸が主で第二軸が少し混ざる。SI の数値では $\lambda=(99\times4+9)/100=4.05$、残差の二乗は $0.99(4-4.05)^2+0.01(9-4.05)^2=0.2475$。少なくとも一つの根が $4.05\pm0.498\,\mathrm{s^{-2}}$ にある。第二固有値が $9$ と分かっていれば、外部の重みは $0.2475/4.95^2\simeq0.0101$ 以下であり、実際の $0.0100$ と照合できる。周波数だけの誤差を振幅比の誤差と同じとみなすと、この混入量を正しく評価できない。
:::

## 演習問題

::::exercise{level=basic}
$V=k(x-a)^2/2+Cx$、$k>0$ の平衡位置と質量 $m$ の微小振動数を求めよ。

:::solution
$V'=k(x-a)+C=0$ より $x_0=a-C/k$。
曲率 $V''=k$ に変化がないので角振動数は $\sqrt{k/m}$ である。
:::
::::

::::exercise{level=basic}
$m=0.5\,\mathrm{kg}, k=8\,\mathrm{N/m},k_c=5\,\mathrm{N/m}$ の連成二振動子の二つの角振動数を求めよ。

:::solution
同位相で $\sqrt{8/0.5}=4\,\mathrm{s^{-1}}$、逆位相で $\sqrt{(8+10)/0.5}=6\,\mathrm{s^{-1}}$ である。
:::
::::

::::exercise{level=basic}
等質量二振動子の初期変位 $(d,-d)$、速度 $(u,u)$ を時間発展させよ。

:::solution
変位は逆位相、速度は同位相にだけ射影される。
従って $x_1=d\cos\omega_-t+(u/\omega_+)\sin\omega_+t$、
$x_2=-d\cos\omega_-t+(u/\omega_+)\sin\omega_+t$ であり、微分して初期速度も一致する。
:::
::::

::::exercise{level=basic}
直線三原子で $m_c=2m$ の非零角振動数の比と高周波の振幅比を求めよ。

:::solution
高い方の二乗は $k/m+2k/(2m)=2k/m$、低い方は $k/m$。
周波数比は $\sqrt2$、高い方の形は $(1,-1,1)$ である。
:::
::::

::::exercise{level=advanced}
$M$ 正定値のもと、二つの固有値が等しい場合にも質量直交化で固有方程式が保たれることを示せ。

:::solution
同じ $\lambda$ の解 $a,b$ の任意の線形結合は $Kv=\lambda Mv$ を満たす。
$b'=b-a(a^\top Mb)/(a^\top Ma)$ とすれば $a^\top Mb'=0$。
独立な $b$ を選べば $b'\ne0$、正定値性により質量ノルムで割れる。
これを反復すればよい。
:::
::::

::::exercise{level=advanced}
平衡付近で $K$ に負方向があると指数的に増える解があることを示せ。

:::solution
$v^\top Kv<0$ の $v$ を質量正規直交モードで展開すると $\sum_s\lambda_s c_s^2<0$。
従ってある $\lambda_j<0$。
そのモードで $Q_j=\exp(\sqrt{-\lambda_j}t)$ と選べば運動方程式を満たす。
:::
::::

::::exercise{level=advanced}
壁のない三原子系を初速零で $(0,d,0)$ から放す。
$m_c=2m$ のときの位置を求めよ。

:::solution
$(0,d,0)=(d/2)(1,1,1)-(d/2)(1,-1,1)$。
前者は静止した並進モード、後者は $\omega=\sqrt{2k/m}$。
従って両端は $d(1-\cos\omega t)/2$、中央は $d(1+\cos\omega t)/2$ である。
:::
::::

::::exercise{level=advanced}
座標を $\eta=Bz$（定数可逆行列）と変えても角振動数が変わらないことを証明せよ。

:::solution
新しい行列は $B^\top KB,B^\top MB$。
行列式の積の公式より $\det[B^\top(K-\lambda M)B]=(\det B)^2\det(K-\lambda M)$ である。
係数が零でないので根と重複度が一致する。
:::
::::

::::exercise{level=advanced}
$M=\diag(1,2)$、$K=\begin{pmatrix}2&-1\\-1&3\end{pmatrix}$ を SI の数値で用いる。
第一質点だけを $d$ ずらし、静かに放した後の各点の変位を求めよ。

:::solution
$(1,1)$ と $(2,-1)$ の形の固有値はそれぞれ $1,5/2$。
初期変位を $(d,0)=a(1,1)+b(2,-1)$ とすると、第二成分より $a=b$、第一成分より $a=b=d/3$ である。
各初速が零なので

$$
\begin{aligned}
x_1(t)&=\tfrac d3[\cos t+2\cos(\sqrt{5/2}\,t)],\\
x_2(t)&=\tfrac d3[\cos t-\cos(\sqrt{5/2}\,t)]
\end{aligned}
$$

となる（$t$ は秒単位の数値）。
$t=0$ の加速度は $(-2d,d/2)$ となり、$-M^{-1}K(d,0)^\top$ に一致する。
:::
::::

::::exercise{level=advanced}
物体の向きによらない二次元の等方振動子で、$M=mI_2$、$K=kI_2+\epsilon\begin{pmatrix}0&h\\h&0\end{pmatrix}$ とする。
小摂動でどの軸が選ばれ、周波数がどれだけ分かれるか。

:::solution
$\epsilon=0$ では全方向で $\omega_0=\sqrt{k/m}$。
非対角項は $(1,1),(1,-1)$ を固有軸に選び、二つの $\omega_\pm^2=(k\pm\epsilon h)/m$ を与える。
$|\epsilon h|<k$ なら両方向とも安定で、

$$
\omega_\pm=\omega_0\left(1\pm\frac{\epsilon h}{2k}\right)+O(\epsilon^2).
$$

従って分裂幅は一次まで $\omega_0|\epsilon h|/k$。
もとの $x,y$ 軸で対角成分だけを読めばどちらも補正零となり、この分裂を見落とす。
:::
::::

::::exercise{level=advanced}
二つの直交方向に同じ質量 $m$ と復元係数 $k$ をもつ振動子を磁場に入れる。速度力を $m\omega_c(\dot y,-\dot x)$ とする。二つの正の固有振動数の積と差を求め、$m=1\,\mathrm{kg},k=4\,\mathrm{N/m},\omega_c=3\,\mathrm{s^{-1}}$ を代入せよ。

:::solution
$Z=x+iy$ を用いると $\ddot Z+i\omega_c\dot Z+(k/m)Z=0$。$Z=e^{-i\omega t}$ に対する根は $\omega=(\omega_c\pm\sqrt{\omega_c^2+4k/m})/2$。二根は異符号で逆向きの回転を表すから、絶対値をとった二周波数の積は $k/m$、大きい方から小さい方を引けば $|\omega_c|$ となる。数値は $4$ と $1\,\mathrm{s^{-1}}$。磁場を反転しても大小の値は変わらず、回転の向きが交換される。
:::
::::

## 確認テスト

::::quiz{answer="5" tol="1e-6"}
$m=2\,\mathrm{kg},k=50\,\mathrm{N/m}$ のばね振動子の角振動数を毎秒で答えよ。

:::explanation
$\omega^2=k/m=25\,\mathrm{s^{-2}}$ の正の平方根は $5\,\mathrm{s^{-1}}$ である。
:::
::::

::::quiz
異なる固有値をもつ基準振動が必ず直交する内積はどれか。

- [ ] $a^\top b$ が常に適切
- [x] $a^\top Mb$
- [ ] $a^\top M^{-1}b$ が常に適切

:::explanation
運動エネルギーの交差項は $a^\top Mb$ に比例するため、質量行列で重みを付ける。
:::
::::

::::quiz{answer="1" tol="1e-6"}
軸方向だけを許した自由な直線三原子の零モード数はいくつか。
結合ばねはともに正とする。

:::explanation
弾性エネルギーが零なのは $x_1=x_2=x_3$ に限られるので、核は一次元である。
:::
::::

::::quiz
Rayleigh 商の最小値は何に等しいか。

- [ ] 最大角振動数
- [x] 一般固有値の最小値
- [ ] 質量の最小値

:::explanation
質量規格化した基底では Rayleigh 商は固有値の重み付き平均であり、最小固有値の固有ベクトルで下限に達する。
:::
::::

::::quiz
自由な直線三原子分子で両端が同質量である。
中央原子の同位体だけを重くし、二次の力の定数は固定する。
どの角振動数が変わらないか。

- [ ] 逆対称伸縮
- [x] 対称伸縮
- [ ] 二方向の曲げ

:::explanation
対称伸縮の変位比は $(1,0,-1)$ で中央の運動エネルギーが零になる。
固有値は $k/m$ であり中央質量を含まない。
他の二種類では中央原子も動く。
:::
::::

::::quiz{answer="2" tol="1e-6"}
二質点模型で $\omega_+=4,\omega_-=8$（毎秒）である。
二つの質量規格化した基準座標の初速は零で、低周波振幅は高周波振幅の何倍なら両エネルギーが等しいか。
振幅の絶対値で答えよ。

:::explanation
静止して放せば $E_s=\omega_s^2 Q_s(0)^2/2$。
従って $|Q_+(0)|/|Q_-(0)|=\omega_-/\omega_+=2$ である。
:::
::::

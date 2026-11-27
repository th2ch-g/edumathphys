---
title: 剛体の運動
description: 姿勢と角速度を区別し、慣性主軸に沿う運動方程式から自由回転と重い対称こまを解く。
notation:
  - symbol: 'R'
    meaning: 物体系から空間系への回転行列
  - symbol: '\boldsymbol\omega'
    meaning: 角速度
  - symbol: 'I_1,I_2,I_3'
    meaning: 主慣性モーメント
  - symbol: '\phi,\theta,\psi'
    meaning: 節線を用いる三つの Euler 角
---

## この章の目標

- 剛体の速度と運動エネルギーを並進と回転に分ける。
- 慣性主軸を求め、角運動量と角速度の向きを比べる。
- Euler 方程式を導き、三つの主軸回転の安定性を判定する。
- 対称こまの章動、歳差、定常歳差を保存量から区別する。

## 剛体の自由度と運動の分解

### 距離から姿勢を定める

互いの距離が変わらない質点の集まりを :term[剛体]{en="rigid body" yomi="ごうたい"} と呼ぶ。
ここでは変形や内部振動の時間尺度を分離し、質量分布を一定とみなす。
一点を定め、その点から別の一点への向きを選ぶと二つの角が必要になる。
一直線上にない第三点を指定すると、その直線のまわりの回転も決まる。
鏡映を除く一つの連結成分では、重心の三座標と姿勢の三自由度で配置が決まる。
質点が一直線上だけに並ぶ模型では直線のまわりの回転を位置から識別できない。

空間に固定した右手系と剛体に固定した右手系を選び、$R$ を物体系の成分から空間系の成分への行列とする。
長さと向きを保つ条件は $R^\top R=I,\det R=1$ である。
物体系の位置 $\rho_b$ を重心から測ると、

$$
\vb r_b=\vb X+R\rho_b,\qquad \sum_bm_b\rho_b=0.
$$

$RR^\top=I$ を微分すると $\dot R R^\top$ は反対称である。
三次元の反対称行列は三つの独立成分を持ち、ただ一つの $\boldsymbol\omega$ に対して

$$
\dot R R^\top=
\begin{pmatrix}0&-\omega_z&\omega_y\\\omega_z&0&-\omega_x\\-\omega_y&\omega_x&0\end{pmatrix}
$$

と書ける。
この行列の作用は $\boldsymbol\omega\times$ である。
こうして :term[角速度]{en="angular velocity" yomi="かくそくど"} を姿勢の微分から定義できる。

### 輸送公式とエネルギー

物体成分 $a(t)$ をもつ任意のベクトルを微分すると、

$$
\frac{\dd}{\dd t}(Ra)=R\dot a+\boldsymbol\omega\times(Ra).
\label{eq-am-transport}
$$

物体系に対して動く点の加速度も、同じ公式を二回使えば得られる。原点を $\vb X$、相対位置の物体成分を $a(t)$ とすると
$$
\ddot{\vb r}=\ddot{\vb X}+R\ddot a+
2\boldsymbol\omega\times(R\dot a)+
\dot{\boldsymbol\omega}\times(Ra)+
\boldsymbol\omega\times(\boldsymbol\omega\times Ra).
$$
一つの $\boldsymbol\omega\times(R\dot a)$ は $R\dot a$ の微分から、もう一つは $\boldsymbol\omega\times Ra$ の微分から出る。その和が Coriolis 加速度である。例えば一定の $\boldsymbol\omega=\Omega\vb e_z$ の回転台で、中心から半径方向に速さ $v_r$ で動く点では、順回転方向の加速度に $2\Omega v_r$ が寄与する。物体に固定した点では $a$ は一定なので、この項は消える。

特に各点の速度は $\dot{\vb X}+\boldsymbol\omega\times R\rho_b$。
二乗して和をとると交差項は $\sum_bm_bR\rho_b=0$ により消え、

$$
T=\tfrac12 M|\dot{\vb X}|^2+\tfrac12\sum_bm_b|\boldsymbol\omega\times R\rho_b|^2 .
$$

重心以外を基準点に選ぶ場合も、空間で静止する支点ならその点の速度が零なので回転の式だけを使える。

:::example[転がる円板]
半径 $a$ の円板が水平面上をすべらず転がる。
重心速度 $v$ と角速度の大きさは $v=a\omega$ である。
最下点の並進速度と回転速度は相殺され、最上点では加わる。
$v=2\,\mathrm{m/s}$ ならそれぞれ $0$ と $4\,\mathrm{m/s}$。
瞬間的に静止する接触点は支点が固定されたヒンジと異なり、時刻とともに別の物体点に移る。
:::

### 瞬間の回転軸と有限回転

ある時刻に基準点の速度を $\vb v_O$ とする。
$\boldsymbol\omega\ne0$ なら、基準点から

$$
\vb d=\frac{\boldsymbol\omega\times\vb v_O}{|\boldsymbol\omega|^2}
$$

だけ離れた点の速度は $\vb v_O+\boldsymbol\omega\times\vb d$ である。
三重積の公式を用いると、これは $\boldsymbol\omega(\boldsymbol\omega\cdot\vb v_O)/|\boldsymbol\omega|^2$、
すなわち角速度に平行な成分だけになる。
$\vb d$ に角速度方向の任意の変位を加えても同じである。
この直線を軸とした回転と、軸方向の並進で瞬間の運動を表せる。
転がる円板では内積が零なので、軸上の点の速度も零である。

有限な姿勢の差も、適切な一本の軸を選べば回転角一つで表せる。これが :term[Euler の回転定理]{en="Euler's rotation theorem" yomi="おいらーのかいてんていり"} である。実際、二つの姿勢の差を表す回転行列 $U$ の固有値は絶対値が一である。三次元の実行列では非実固有値は共役対になり、その積は一である。従って $\det U=1$ なら、残る実固有値は $+1$ となる。全固有値が実でも $-1$ の個数は偶数なので、少なくとも一本の不動軸がある。軸に垂直な平面は直交性により不変で、そこへの制限は行列式一の二次元回転である。恒等変換以外では固定軸は一本に定まる。

一定な単位軸 $\vb n$ と回転角 $\alpha$ に対して、任意のベクトルを $\vb n$ に平行な部分と垂直な部分に分ける。
平行部分は変わらず、垂直部分は平面内で $\alpha$ 回るから

$$
R_{\vb n}(\alpha)\vb v
=\vb v\cos\alpha+(\vb n\times\vb v)\sin\alpha+
\vb n(\vb n\cdot\vb v)(1-\cos\alpha).
$$

この Rodrigues の公式を $\alpha=0$ で微分すれば瞬間の回転 $\vb n\times\vb v$ が戻る。
異なる軸の有限回転では順序が重要である。
例えば $\vb e_z$ をまず $x$ 軸まわりに $90^\circ$ 回し、次に $z$ 軸まわりに $90^\circ$ 回すと $\vb e_x$。
逆順では $-\vb e_y$ となる。

### 瞬間ねじ軸と軸上の速度

回転と並進を分ける基準点は任意に選べるが、ある点で回転の寄与が消えても、その点が空間に固定されているとは限らない。ある瞬間に物体上の基準点 $O$ の速度が $\vb v_O$ であれば、そこから $\vb r$ だけ離れた点の速度は $\vb v_O+\boldsymbol\omega\times\vb r$ である。$\omega\ne0$ のとき、軸に沿う速度は全点で共通である。軸に垂直な速度だけが消える点は
$$
\vb r_c=\frac{\boldsymbol\omega\times\vb v_O}{\omega^2},\qquad
\vb r=\vb r_c+\lambda\boldsymbol\omega .
\label{eq-am-instant-screw}
$$
実際、三重積を展開すると $\boldsymbol\omega\times\vb r_c=-\vb v_O+\boldsymbol\omega(\boldsymbol\omega\cdot\vb v_O)/\omega^2$。軸上の速度は $h\boldsymbol\omega$、$h=(\boldsymbol\omega\cdot\vb v_O)/\omega^2$ となる。$h$ は一ラジアンあたりの軸方向変位であり、瞬間の速度場が一定の間だけ有限変位にそのまま使える。$\omega=0$ の純並進はこの表示の例外である。

:::example[接触点が静止する転がり]
$\boldsymbol\omega=\left(\begin{matrix}0 & 2 & 0\end{matrix}\right)^\top\,\mathrm{s^{-1}}$、$\vb v_O=\left(\begin{matrix}3 & 0 & 0\end{matrix}\right)^\top\,\mathrm{m/s}$ とする。軸の位置を決めるために、速度を回転方向に平行な成分と垂直な成分に分解する。射影の長さを角速度の大きさで割れば $h=0\,\mathrm m$。外積の式から

$$
\vb r_c=\left(\begin{matrix}0 & 0 & - \frac{3}{2}\end{matrix}\right)^\top\,\mathrm m,\quad \boldsymbol\omega\times\vb r_c=\left(\begin{matrix}-3 & 0 & 0\end{matrix}\right)^\top\,\mathrm{m/s},\quad \vb v_c=\left(\begin{matrix}0 & 0 & 0\end{matrix}\right)^\top\,\mathrm{m/s}.
$$

軸から垂直距離 $d$ の点では、軸上の速度とこれに垂直な速さ $|\omega|d$ が加わる。従って $d=0.50\,\mathrm m$ の位置での速さは $1\,\mathrm{m/s}$。軸は基準点から下へずれる。基準点を車軸上の中心とするなら、接触点が静止する車輪の半径は $1.5\,\mathrm m$ である。直前の $0.50\,\mathrm m$ は、ねじ軸から観測点までの距離である。
:::

:::example[軸に平行な送りがある回転]
$\boldsymbol\omega=\left(\begin{matrix}0 & 0 & 4\end{matrix}\right)^\top\,\mathrm{s^{-1}}$、$\vb v_O=\left(\begin{matrix}0 & 2 & 1\end{matrix}\right)^\top\,\mathrm{m/s}$ とする。軸の位置を決めるために、速度を回転方向に平行な成分と垂直な成分に分解する。射影の長さを角速度の大きさで割れば $h=\frac{1}{4}\,\mathrm m$。外積の式から

$$
\vb r_c=\left(\begin{matrix}- \frac{1}{2} & 0 & 0\end{matrix}\right)^\top\,\mathrm m,\quad \boldsymbol\omega\times\vb r_c=\left(\begin{matrix}0 & -2 & 0\end{matrix}\right)^\top\,\mathrm{m/s},\quad \vb v_c=\left(\begin{matrix}0 & 0 & 1\end{matrix}\right)^\top\,\mathrm{m/s}.
$$

軸から垂直距離 $d$ の点では、軸上の速度とこれに垂直な速さ $|\omega|d$ が加わる。従って $d=0.50\,\mathrm m$ の位置での速さは $\sqrt{5}\,\mathrm{m/s}$。軸上の点も毎秒送り方向に動くため、空間に静止する点は存在しない。ねじのピッチを一周あたりで指定する場合は $2\pi h$ に直す。
:::

:::example[座標軸に平行でない回転]
$\boldsymbol\omega=\left(\begin{matrix}1 & 2 & 2\end{matrix}\right)^\top\,\mathrm{s^{-1}}$、$\vb v_O=\left(\begin{matrix}2 & -1 & 3\end{matrix}\right)^\top\,\mathrm{m/s}$ とする。軸の位置を決めるために、速度を回転方向に平行な成分と垂直な成分に分解する。射影の長さを角速度の大きさで割れば $h=\frac{2}{3}\,\mathrm m$。外積の式から

$$
\vb r_c=\left(\begin{matrix}\frac{8}{9} & \frac{1}{9} & - \frac{5}{9}\end{matrix}\right)^\top\,\mathrm m,\quad \boldsymbol\omega\times\vb r_c=\left(\begin{matrix}- \frac{4}{3} & \frac{7}{3} & - \frac{5}{3}\end{matrix}\right)^\top\,\mathrm{m/s},\quad \vb v_c=\left(\begin{matrix}\frac{2}{3} & \frac{4}{3} & \frac{4}{3}\end{matrix}\right)^\top\,\mathrm{m/s}.
$$

軸から垂直距離 $d$ の点では、軸上の速度とこれに垂直な速さ $|\omega|d$ が加わる。従って $d=0.50\,\mathrm m$ の位置での速さは $\frac{5}{2}\,\mathrm{m/s}$。三軸すべての速度が観測されても、速度の大きさが最小になる点は直線を作る。その最小値は共通の軸方向速度の大きさである。
:::

## 慣性テンソルと主軸

### 軸と角運動量

以下では支点または重心から測った物体成分を $\rho_b$、物体系での角速度を $\omega=R^\top\boldsymbol\omega$ とする。
三重積 $\rho\times(\omega\times\rho)=|\rho|^2\omega-(\rho\cdot\omega)\rho$ から、回転の角運動量は $\vb L_{\rm body}=I\omega$ となる。
ただし

$$
I_{ij}=\sum_bm_b\bigl(|\rho_b|^2\delta_{ij}-\rho_{bi}\rho_{bj}\bigr),\qquad T_{\rm rot}=\tfrac12\omega^\top I\omega.
\label{eq-am-inertia}
$$

この対称行列を :term[慣性テンソル]{en="inertia tensor" yomi="かんせいてんそる"} という。
単位は $\mathrm{kg\,m^2}$ で、連続体では質量和を $\int \dd m$ に置き換える。
任意の軸 $\vb n$ に関して $\vb n^\top I\vb n=\sum_bm_b|\vb n\times\rho_b|^2\ge0$。
零になるのは全質量が基準点を通るその軸に集中する場合に限る。従ってそのような軸が存在しなければ正定値である。

:ref[thm-am-normal-modes] の証明で示した実対称行列の直交対角化により、直交固有軸を選べる。これを :term[慣性主軸]{en="principal axis of inertia" yomi="かんせいしゅじく"} といい、固有値 $I_i$ を主慣性モーメントという。主軸回転では $\vb L$ と角速度が平行になる。縮退した固有平面ではその中のどの単位ベクトルも主軸である。

### 支点を移す

重心から $\vb d$ だけ移した点について $\rho_b-\vb d$ を :eqref[eq-am-inertia] に入れる。
一次の項は質量重み付きの $\rho_b$ の和が零なので消え、

$$
I_{\rm new}=I_{\rm cm}+M(|\vb d|^2 I_{\rm unit}-\vb d\vb d^\top)
$$

となる。
これが平行軸の定理である。
重心を通る軸と、それに平行な移動先の軸との差は $M[|\vb d|^2-(\vb n\cdot\vb d)^2]$ である。ここで $\vb n$ は軸の単位ベクトルであり、角括弧は両軸間の垂直距離の二乗である。重心からずれた二軸を比較するには、それぞれの重心軸からの距離の二乗の差を取る。

:::example[細い棒と円板]
一様な細い棒の中央を通る垂直軸では $I=(M/l)\int_{-l/2}^{l/2}x^2\dd x=Ml^2/12$。
端を通る平行軸では $Ml^2/12+Ml^2/4=Ml^2/3$ となる。
一様円板の対称軸では質量素片 $2Mr\,\dd r/a^2$ を用い $I_3=\int_0^a2Mr^3\dd r/a^2=Ma^2/2$。
平面薄板では $I_x+I_y=I_z$ が被積分関数で成り立つので、円板の直径については各 $Ma^2/4$ である。
:::

:::example[座標軸が主軸でないとき]
質量 $m$ の二点を $\pm(d,d,0)$ に置く。
面内ブロックは $2md^2\begin{pmatrix}1&-1\\-1&1\end{pmatrix}$ である。
$(1,1,0)$ 方向の固有値は零、$(1,-1,0)$ 方向は $4md^2$。
細い直線に沿って質量が並ぶ極限に対応する。
:::

### 対称性を使って積分の前に主軸を選ぶ

質量分布が $x\mapsto-x$ に関して対称なら、$\int xy\,\dd m$ と $\int xz\,\dd m$ はともに零である。
従ってその鏡面に垂直な $x$ 軸は主軸となる。
ただし鏡面内の二つの軸がどちら向きでも主軸になるわけではない。
さらに $z$ 軸まわりの任意角の回転が対称なら、面内で方向を変えても慣性モーメントが同じなので $I_x=I_y$ が従う。

:::example[一様直方体の主軸]
重心を原点とし、長さ $a,b,c$ の三辺を $x,y,z$ 軸に平行に置く。
三つの鏡面により非対角成分は零である。
密度一定なら

$$
\int y^2\,\dd m=\frac Mb\int_{-b/2}^{b/2}y^2\dd y=\frac{Mb^2}{12}.
$$

同様に $I_x=M(b^2+c^2)/12$、$I_y=M(c^2+a^2)/12$、$I_z=M(a^2+b^2)/12$。
例えば $a:b:c=3:2:1$ なら $I_x:I_y:I_z=5:10:13$。
一番長い辺に沿う軸の慣性モーメントが最小である。
:::

:::example[一様な球]
任意の軸を主軸に選べて三つの固有値は等しい。
テンソルの跡は $\tr I=2\int |\rho|^2\dd m$ である。
一様な半径 $a$ の球では $\dd m=3Mr^2\dd r/a^3$ だから

$$
3I_1=2\frac{3M}{a^3}\int_0^a r^4\dd r=\frac{6Ma^2}{5},\qquad I_1=\frac{2Ma^2}{5}.
$$

表面だけに質量がある薄い球殻なら $\int r^2\dd m=Ma^2$ であり、同質量・同半径でも $I_1=2Ma^2/3$ と大きくなる。
:::

### 軸を変えて測った慣性から質量分布を読む

ねじり振り子では一本の軸まわりの慣性モーメントしか直接には測れない。同じ固定点を通る複数の軸について $I(\vb n)=\vb n^\top I\vb n$ を測れば、全テンソルを復元できる。座標軸についての三測定で対角成分を求め、さらに $\vb n=(\vb e_i+\vb e_j)/\sqrt2$ について $I_{ij}=I(\vb n)-(I_{ii}+I_{jj})/2$ を使う。質量を正とする物体では $S=\sum_a m_a\vb r_a\vb r_a^\top$ が半正定値で、$I=(\tr S)1-S$。従って $S=(\tr I)1/2-I$ も半正定値でなければならない。正定値な行列を任意に測定結果として与えてよいわけではない。例えば主値 $(1,1,3)$ は全方向の回転エネルギーを正にするが、第三主値が他二つの和を超え、正の質量では実現しない。

:::example[平面上の三個のおもり]
組 $(m_a,\vb r_a^\top)$ を kg と m の数値で $(1,\left(\begin{matrix}1 & 0 & 0\end{matrix}\right)),\ (2,\left(\begin{matrix}0 & 1 & 0\end{matrix}\right)),\ (1,\left(\begin{matrix}1 & 1 & 0\end{matrix}\right))$ とする。すべての軸を座標原点に通す。位置成分の積を質量で重み付けして足し、そこから慣性を求めると

$$
S=\left(\begin{matrix}2 & 1 & 0\\1 & 3 & 0\\0 & 0 & 0\end{matrix}\right),\qquad I=(\tr S)1-S=\left(\begin{matrix}3 & -1 & 0\\-1 & 2 & 0\\0 & 0 & 5\end{matrix}\right).
$$

対角方向 $xy,xz,yz$ の三軸で測る値はそれぞれ $\frac{3}{2},4,\frac{7}{2}\,\mathrm{kg\,m^2}$ である。例えば $xy$ の測定値から $(I_{xx}+I_{yy})/2=\frac{5}{2}$ を引けば $I_{xy}=-1$ に戻る。主値は特性式 $\det(\lambda 1-I)=\left(\lambda - 5\right) \left(\lambda^{2} - 5 \lambda + 5\right)=0$ の根である。位置の積と慣性乗積は非対角で符号が逆になるので、実測値の再現はこの符号を点検する方法になる。全点が測定原点を通る平面上にあるため、一つの主慣性モーメントは他の二つの和に等しくなる。
:::

:::example[対角方向に質量のある枠]
組 $(m_a,\vb r_a^\top)$ を kg と m の数値で $(2,\left(\begin{matrix}1 & 1 & 0\end{matrix}\right)),\ (2,\left(\begin{matrix}-1 & -1 & 0\end{matrix}\right)),\ (1,\left(\begin{matrix}0 & 0 & 2\end{matrix}\right)),\ (1,\left(\begin{matrix}0 & 0 & -2\end{matrix}\right))$ とする。すべての軸を座標原点に通す。位置成分の積を質量で重み付けして足し、そこから慣性を求めると

$$
S=\left(\begin{matrix}4 & 4 & 0\\4 & 4 & 0\\0 & 0 & 8\end{matrix}\right),\qquad I=(\tr S)1-S=\left(\begin{matrix}12 & -4 & 0\\-4 & 12 & 0\\0 & 0 & 8\end{matrix}\right).
$$

対角方向 $xy,xz,yz$ の三軸で測る値はそれぞれ $8,10,10\,\mathrm{kg\,m^2}$ である。例えば $xy$ の測定値から $(I_{xx}+I_{yy})/2=12$ を引けば $I_{xy}=-4$ に戻る。主値は特性式 $\det(\lambda 1-I)=\left(\lambda - 16\right) \left(\lambda - 8\right)^{2}=0$ の根である。位置の積と慣性乗積は非対角で符号が逆になるので、実測値の再現はこの符号を点検する方法になる。原点は重心である。対称な二点の対は一次のモーメントを消すが、慣性乗積には同符号で寄与する。
:::

:::example[三方向に厚みのある質点系]
組 $(m_a,\vb r_a^\top)$ を kg と m の数値で $(1,\left(\begin{matrix}1 & 0 & 1\end{matrix}\right)),\ (1,\left(\begin{matrix}0 & 1 & 0\end{matrix}\right)),\ (1,\left(\begin{matrix}-1 & 0 & 1\end{matrix}\right)),\ (2,\left(\begin{matrix}0 & 0 & -1\end{matrix}\right))$ とする。すべての軸を座標原点に通す。位置成分の積を質量で重み付けして足し、そこから慣性を求めると

$$
S=\left(\begin{matrix}2 & 0 & 0\\0 & 1 & 0\\0 & 0 & 4\end{matrix}\right),\qquad I=(\tr S)1-S=\left(\begin{matrix}5 & 0 & 0\\0 & 6 & 0\\0 & 0 & 3\end{matrix}\right).
$$

対角方向 $xy,xz,yz$ の三軸で測る値はそれぞれ $\frac{11}{2},4,\frac{9}{2}\,\mathrm{kg\,m^2}$ である。例えば $xy$ の測定値から $(I_{xx}+I_{yy})/2=\frac{11}{2}$ を引けば $I_{xy}=0$ に戻る。主値は特性式 $\det(\lambda 1-I)=\left(\lambda - 6\right) \left(\lambda - 5\right) \left(\lambda - 3\right)=0$ の根である。位置の積と慣性乗積は非対角で符号が逆になるので、実測値の再現はこの符号を点検する方法になる。座標軸が主軸でも主慣性モーメントは等しくない。厚みをもつため三角不等式はすべて厳密になる。
:::

## Euler 角

### 回転の順序を固定する

姿勢を :term[Euler 角]{en="Euler angles" yomi="おいらーかく"} で示す。
右手系の能動回転で

$$
R=R_z(\phi)R_x(\theta)R_z(\psi)
$$

とする。
順に、空間鉛直軸まわりの方位、回転後の $x$ 軸（節線）まわりの傾き、物体の第三軸まわりの自転である。
$0<\theta<\pi$ なら節線は鉛直軸と第三軸の双方に直交する。
各因子を一つずつ微分すると、空間角速度は三つの瞬間回転軸に沿うベクトルの和となる。
物体系に戻すと、鉛直軸の成分は $(\sin\theta\sin\psi,\sin\theta\cos\psi,\cos\theta)$、節線の成分は $(\cos\psi,-\sin\psi,0)$。
従って

$$
\begin{aligned}
\omega_1&=\dot\phi\sin\theta\sin\psi+\dot\theta\cos\psi,\\
\omega_2&=\dot\phi\sin\theta\cos\psi-\dot\theta\sin\psi,\\
\omega_3&=\dot\phi\cos\theta+\dot\psi .
\end{aligned}
\label{eq-am-euler-angles}
$$

角速度は単位 $\mathrm{rad/s}$ のベクトルで、三つの角度の時間微分そのものを成分とするベクトルではない。

:::example[座標の極で失われる情報]
$\theta=0$ で $R=R_z(\phi+\psi)$ となるため、$\phi$ と $\psi$ の個別の値は定まらない。
$\dot\theta=0,\dot\phi=2,\dot\psi=3$（毎秒）なら $\omega_3=5\,\mathrm{s^{-1}}$ である。
同じ姿勢変化を $\dot\phi=0,\dot\psi=5$ と書いてもよい。
剛体の運動が特異になったわけではない。
:::

:::example[一定の傾きでの歳差]
$\theta=\pi/3,\psi=0,\dot\theta=0,\dot\phi=2,\dot\psi=3$ とする。
その瞬間の物体成分は $(0,\sqrt3,4)\,\mathrm{s^{-1}}$ である。
第三成分には、鉛直まわりの回転の射影も含まれる。
:::

### 角速度から観測角を復元する

横成分の式を $\cos\psi,-\sin\psi$ の組と $\sin\psi,\cos\psi$ の組に射影すると

$$
\begin{aligned}
\dot\theta&=\omega_1\cos\psi-\omega_2\sin\psi,\\
\dot\phi&=(\omega_1\sin\psi+\omega_2\cos\psi)/\sin\theta,\\
\dot\psi&=\omega_3-\dot\phi\cos\theta .
\end{aligned}
$$

計測した $\omega$ を姿勢に直すときには、現在の $\theta,\psi$ も必要である。
例えば $\theta=\pi/6,\psi=\pi/2$、$\omega=(1,2,3)\,\mathrm{s^{-1}}$ なら $\dot\theta=-2,\ \dot\phi=2,\ \dot\psi=3-\sqrt3$（毎秒）となる。
$\theta$ が小さいと、有限な横成分でも $\dot\phi$ と $\dot\psi$ が個別には大きくなりうる。
この座標の悪条件を避けたいときは、代わりに $\dot R=R[\omega]_\times$ を直接積分し、$R^\top R=I$ を確認する。

### 共役運動量を軸への射影として読む

一般の剛体で $T=\sum_i I_i\omega_i^2/2$ と書く。
Euler 角の任意の一つを $\alpha$ とすると

$$
p_\alpha=\sum_{i=1}^3 I_i\omega_i\frac{\partial\omega_i}{\partial\dot\alpha}.
$$

$\partial\omega/\partial\dot\alpha$ は、$\alpha$ を動かしたときの瞬間回転軸を物体系で表した単位ベクトルである。
従って $p_\phi$ は鉛直への射影、$p_\psi$ は第三軸への射影、$p_\theta$ は節線への射影となる。
どの射影が保存するかはトルクと対称性で別途決まり、角度を導入しただけで三つが保存するわけではない。

### Euler パラメータで有限回転を合成する

傾きが零を通り抜ける計算では Euler 角の二つの大きな速度が打ち消し合う。そこで単位 :term[四元数]{en="quaternion" yomi="しげんすう"} $a=(a_0,\vb a)$ で姿勢を記す。積を
$$
(a_0,\vb a)(b_0,\vb b)=(a_0b_0-\vb a\cdot\vb b,\
a_0\vb b+b_0\vb a+\vb a\times\vb b)
$$
と定義する。結合則を確かめるため、各四元数を複素行列
$$
M(a)=\begin{pmatrix}a_0-ia_3&-a_2-ia_1\\a_2-ia_1&a_0+ia_3\end{pmatrix}
$$
に対応させる。純ベクトルの基底に対応する行列を順に $E_1,E_2,E_3$ とおく。直接掛けると $E_j^2=-I_2$、$E_1E_2=E_3$、$E_2E_3=E_1$、$E_3E_1=E_2$ となり、逆順の積は符号が反転する。よって積を成分ごとに展開すれば $M(ab)=M(a)M(b)$ である。行列の成分から四つの実数が復元できるので、行列積の結合則は四元数にも成り立つ。共役を $a^*=(a_0,-\vb a)$ とすれば $aa^*=(a_0^2+\vb a^2,\vb0)$。行列式もこの長さの二乗である。単位四元数では $M(a)^*M(a)=I_2$ かつ $\det M(a)=1$ となるので、二つの複素成分で姿勢を表す :term[Cayley–Klein パラメータ]{en="Cayley–Klein parameters" yomi="けーりーくらいんぱらめーた"} が得られる。任意の $SU(2)$ 行列の第一列を $(\alpha,\beta)^\top$ とすると、直交性と行列式一から第二列は $(-\bar\beta,\bar\alpha)^\top$。従ってこれは上の $M(a)$ の形であり、単位四元数全体と一対一に対応する。

ベクトルを $(0,\vb x)$ として挟む積 $a(0,\vb x)a^*$ は純ベクトルになり、その成分は
$$
R(a)\vb x=(a_0^2-|\vb a|^2)\vb x+2\vb a(\vb a\cdot\vb x)+2a_0\vb a\times\vb x .
\label{eq-am-quaternion-rotation}
$$
$a=(\cos(\theta/2),\vb n\sin(\theta/2))$ を代入すると軸方向を固定し、垂直方向を有向角 $\theta$ だけ回す Rodrigues の公式になる。$a,-a$ は同じ姿勢である。右から純ベクトルを掛けて $\dot a=a(0,\boldsymbol\omega_{\rm body})/2$ とすれば、微小な物体軸回転の合成になる。空間軸角速度を使う場合は左から掛ける。どちらも $a_0\dot a_0+\vb a\cdot\dot{\vb a}=0$ を満たし、単位長さを保つ。

:::example[異なる空間軸の四分の一回転]
初期姿勢を単位行列とし、空間固定の単位軸 $\left(\begin{matrix}1 & 0 & 0\end{matrix}\right)^\top$ のまわりに $\frac{\pi}{2}$、続いて $\left(\begin{matrix}0 & 1 & 0\end{matrix}\right)^\top$ のまわりに $\frac{\pi}{2}$ だけ能動回転する。各回転をスカラー部が先の四成分の列で表すと

$$
a=\left(\begin{matrix}\frac{\sqrt{2}}{2}\\\frac{\sqrt{2}}{2}\\0\\0\end{matrix}\right),\qquad b=\left(\begin{matrix}\frac{\sqrt{2}}{2}\\0\\\frac{\sqrt{2}}{2}\\0\end{matrix}\right),\qquad ba=\left(\begin{matrix}\frac{1}{2}\\\frac{1}{2}\\\frac{1}{2}\\- \frac{1}{2}\end{matrix}\right).
$$

ベクトルを挟む積を二回行うので最後の回転が左に来る。この積から得られる回転行列と、物体につけた印の移る先は

$$
R(ba)=\left(\begin{matrix}0 & 1 & 0\\0 & 0 & -1\\-1 & 0 & 0\end{matrix}\right),\qquad R(ba)\left(\begin{matrix}0\\0\\1\end{matrix}\right)=\left(\begin{matrix}0\\-1\\0\end{matrix}\right).
$$

順序を逆にした列は $ab=\left(\begin{matrix}\frac{1}{2} & \frac{1}{2} & \frac{1}{2} & \frac{1}{2}\end{matrix}\right)^\top$ である。水平に倒した印をどちらへ送るかは二番目の軸で決まる。逆順の積とスカラー部が同じでも回転軸は異なる。
:::

:::example[直交する半回転を続ける]
初期姿勢を単位行列とし、空間固定の単位軸 $\left(\begin{matrix}1 & 0 & 0\end{matrix}\right)^\top$ のまわりに $\pi$、続いて $\left(\begin{matrix}0 & 1 & 0\end{matrix}\right)^\top$ のまわりに $\pi$ だけ能動回転する。各回転をスカラー部が先の四成分の列で表すと

$$
a=\left(\begin{matrix}0\\1\\0\\0\end{matrix}\right),\qquad b=\left(\begin{matrix}0\\0\\1\\0\end{matrix}\right),\qquad ba=\left(\begin{matrix}0\\0\\0\\-1\end{matrix}\right).
$$

ベクトルを挟む積を二回行うので最後の回転が左に来る。この積から得られる回転行列と、物体につけた印の移る先は

$$
R(ba)=\left(\begin{matrix}-1 & 0 & 0\\0 & -1 & 0\\0 & 0 & 1\end{matrix}\right),\qquad R(ba)\left(\begin{matrix}1\\1\\1\end{matrix}\right)=\left(\begin{matrix}-1\\-1\\1\end{matrix}\right).
$$

順序を逆にした列は $ab=\left(\begin{matrix}0 & 0 & 0 & 1\end{matrix}\right)^\top$ である。逆順の四元数は符号が違うが姿勢は等しい。三次元の有限回転がいつでも可換であるという意味ではない。
:::

:::example[一軸の回転をちょうど打ち消す]
初期姿勢を単位行列とし、空間固定の単位軸 $\left(\begin{matrix}0 & 0 & 1\end{matrix}\right)^\top$ のまわりに $\frac{\pi}{3}$、続いて $\left(\begin{matrix}0 & 0 & 1\end{matrix}\right)^\top$ のまわりに $- \frac{\pi}{3}$ だけ能動回転する。各回転をスカラー部が先の四成分の列で表すと

$$
a=\left(\begin{matrix}\frac{\sqrt{3}}{2}\\0\\0\\\frac{1}{2}\end{matrix}\right),\qquad b=\left(\begin{matrix}\frac{\sqrt{3}}{2}\\0\\0\\- \frac{1}{2}\end{matrix}\right),\qquad ba=\left(\begin{matrix}1\\0\\0\\0\end{matrix}\right).
$$

ベクトルを挟む積を二回行うので最後の回転が左に来る。この積から得られる回転行列と、物体につけた印の移る先は

$$
R(ba)=\left(\begin{matrix}1 & 0 & 0\\0 & 1 & 0\\0 & 0 & 1\end{matrix}\right),\qquad R(ba)\left(\begin{matrix}1\\0\\1\end{matrix}\right)=\left(\begin{matrix}1\\0\\1\end{matrix}\right).
$$

順序を逆にした列は $ab=\left(\begin{matrix}1 & 0 & 0 & 0\end{matrix}\right)^\top$ である。同じ軸に限れば角度の加減算が積の公式から出る。単位四元数への復帰は角速度積分の符号を点検する例になる。
:::

## Euler の運動方程式

### 回転する観測軸から見たトルク

重心または空間固定点まわりで $\dd\vb L_{\rm space}/\dd t=\vb N_{\rm space}$ である。
物体系の成分を用い :eqref[eq-am-transport] により

$$
I\dot\omega+\omega\times(I\omega)=N .
$$

主軸で成分に分けた式が :term[Euler の運動方程式]{en="Euler equations" yomi="おいらーのうんどうほうていしき"} である。

$$
\begin{aligned}
I_1\dot\omega_1+(I_3-I_2)\omega_2\omega_3&=N_1,\\
I_2\dot\omega_2+(I_1-I_3)\omega_3\omega_1&=N_2,\\
I_3\dot\omega_3+(I_2-I_1)\omega_1\omega_2&=N_3 .
\end{aligned}
\label{eq-am-euler-body}
$$

$\omega$ との内積をとると $\dot T_{\rm rot}=N\cdot\omega$。
また $I\omega$ との内積から $\dd |I\omega|^2/\dd t=2(I\omega)\cdot N$。
自由回転では物体系で $\vb L$ が動いても、これら二つは一定である。

### 中間軸での不安定性

$I_1<I_2<I_3$、外力モーメント零とする。
第三軸まわりの回転 $(0,0,\Omega)$ を少し乱すと、

$$
I_1\delta\dot\omega_1=(I_2-I_3)\Omega\,\delta\omega_2,\qquad
I_2\delta\dot\omega_2=(I_3-I_1)\Omega\,\delta\omega_1 .
$$

よって $\delta\ddot\omega_1=-\Omega^2(I_3-I_2)(I_3-I_1)\delta\omega_1/(I_1I_2)$。
振動するため最大軸は線形安定であり、番号を巡回すると最小軸も安定となる。
中間軸では

$$
\delta\ddot\omega_1=\Omega^2\frac{(I_2-I_1)(I_3-I_2)}{I_1I_3}\delta\omega_1
$$

となり成長解がある。

端の二軸で安定性が非線形にも保たれる理由は、保存する $L^2-2I_3T=\sum_i I_i(I_i-I_3)\omega_i^2$ を見ると分かる。
第三軸以外の二項はともに負で、初期に小さい横成分は永久に小さい。
$T$ の保存が縦成分の変化も抑える。
最小軸では符号を逆にして同じ議論が使える。
中間軸では二項の符号が異なり、この拘束は働かない。

:::example[回転の乱れが育つ時間]
$I_1:I_2:I_3=2:3:4$、$\Omega=4\,\mathrm{s^{-1}}$ なら中間軸の成長率は $4/\sqrt8=\sqrt2\,\mathrm{s^{-1}}$。
小さな初期ずれが $e$ 倍になる時間は約 $0.71\,\mathrm s$ である。
大振幅になった後は線形解を外挿できない。
:::

### 保存面の交線と分離軌道

自由剛体で物体系の角運動量を $L_i=I_i\omega_i$ とすると、運動に沿って

$$
L_1^2+L_2^2+L_3^2=L^2,\qquad
\frac{L_1^2}{I_1}+\frac{L_2^2}{I_2}+\frac{L_3^2}{I_3}=2T
$$

である。
従って運動は球面と楕円体の交線上で起きる。
球面上でエネルギーを比較すれば $L^2/(2I_3)\le T\le L^2/(2I_1)$ であり、極値に近い交線は第三または第一軸を囲む閉曲線になる。
$T=L^2/(2I_2)$ の境界では

$$
(1/I_1-1/I_2)L_1^2=(1/I_2-1/I_3)L_3^2
$$

となり、二つの平面と球の交線が中間軸で交わる。
中間軸回転の線形不安定性をエネルギー面の形からも読み取れる。

さらに $u=L_2$ を選ぶと $L_1^2+L_3^2=L^2-u^2$ と $L_1^2/I_1+L_3^2/I_3=2T-u^2/I_2$ から

$$
L_1^2=\frac{2T-L^2/I_3-(1/I_2-1/I_3)u^2}{1/I_1-1/I_3},
\qquad L_3^2=L^2-u^2-L_1^2 .
$$

Euler 方程式で $\dot u=(1/I_1-1/I_3)L_3L_1$ なので、$\dd t=\dd u/[(1/I_1-1/I_3)L_3L_1]$ の一変数積分に帰着する。
符号は初期値に合わせ、転回ごとに替える。
一般には平方根内が四次式になり、逆関数は楕円関数を用いて表せる。
次に保存面の幾何を調べ、その後でこの積分を反転して自由な非対称こまの周期を求める。

球面の北半分を $L_1,L_3$ 平面へ射影すると、$L_2=\sqrt{L^2-L_1^2-L_3^2}$ で値を戻せる。図では $I_1:I_2:I_3=2:3:4$ としたため、分離平面の式は $L_3=\pm\sqrt2 L_1$ である。一本の軌道が横切って交差する図ではない。交点に厳密に初期値を置けば、物体成分はその点に静止する。

```tikz caption="角運動量球と二つの分離平面を、中間軸の正の側から見た図。中間軸回転は中央で交わる分離曲線の共通点である。"
\draw (0,0) circle (1.8);
\draw[->] (-2.2,0)--(2.4,0) node[right] {$L_1$};
\draw[->] (0,-2.1)--(0,2.25) node[above] {$L_3$};
\draw[thick] (-1.04,-1.47)--(1.04,1.47);
\draw[thick] (-1.04,1.47)--(1.04,-1.47);
\fill (0,0) circle (2pt);
\node at (.90,-.30) {$L_2=L$};
```

### 物体が回る速さを二つの空間で比べる

支点に対する角運動量の空間表示は一定でも、物体上に印を付けて観察すると主軸の向きが変わる。
非対称な自由剛体では、物体系と空間系における $\omega$ の先端の軌跡は別の曲線になる。
運動エネルギーを $T>0$ として、空間固定の $\vb L$ を法線とする平面

$$
\boldsymbol\omega\cdot\vb L=2T
$$

を考える。
この等式は $T=\boldsymbol\omega\cdot\vb L/2$ にほかならない。
平面までの距離は $2T/L$ で一定である。

物体系では $\omega^\top I\omega=2T$ が楕円体になる。
その点での法線は $I\omega$ に平行であり、空間に戻せば $\vb L$ と平行なので、楕円体は常に上の平面に接する。
楕円体を剛体に固定した補助的な図形とみなすと、接点の位置は $\boldsymbol\omega$ 方向にあるから、瞬間の回転速度による接点の速度は $\boldsymbol\omega\times\boldsymbol\omega=0$。
従って楕円体はこの固定平面上をすべらず転がるように見える。
これが Poinsot の幾何学的構成である。
図形を位置空間に描くには角速度の単位を長さに換算する一定係数を掛けるが、接線と回転の幾何は変わらない。

:::example[球形こまの極限]
$I_1=I_2=I_3=C$ の場合は $\vb L=C\boldsymbol\omega$。
Euler 方程式は $C\dot\omega_i=0$ となり、物体でも空間でも角速度方向が変わらない。
補助楕円体は球になり、接点は一点にとどまる。
これは「主軸まわりなら定常」という結果の極限で、すべての軸が主軸となる例である。
:::

### 非対称こまの周期と分離軌道の厳密な時刻表示

回転の一周期をエネルギーと角運動量から計算するには、二つの保存面の交線を時間でたどる。$a=1/I_1>b=1/I_2>c=1/I_3$、$l=|\vb L|$、$d=2T/l^2$ と書く。$u=L_2/l$ を選べば両保存則から
$$
L_1^2/l^2=\frac{d-c-(b-c)u^2}{a-c},\qquad
L_3^2/l^2=\frac{a-d-(a-b)u^2}{a-c}.
$$
$\dot L_2=(a-c)L_1L_3$ を二乗して
$$
\dot u^2=l^2(a-b)(b-c)(A-u^2)(B-u^2),\qquad
A=\frac{d-c}{b-c},\quad B=\frac{a-d}{a-b}.
$$
$c<d<a$、$d\ne b$ なら小さい方の根 $D=\min(A,B)$ と大きい方 $F=\max(A,B)$ はともに正で異なる。$u=\sqrt D\sin\alpha$ と置くと $\dot\alpha^2=\nu^2(1-k^2\sin^2\alpha)$、$\nu=l\sqrt{(a-b)(b-c)F},\ k^2=D/F$ である。初め $u=0,\dot u>0$ と選び、$\nu t=\int_0^\alpha(1-k^2\sin^2 v)^{-1/2}\dd v$ を反転して $u=\sqrt D\,\operatorname{sn}(\nu t,k)$ を得る。これが Jacobi 楕円関数のこの問題で必要な定義である。$L_1,L_3$ は上の平方根で復元し、零を通る成分を滑らかに符号反転する。物体成分の周期は $4K(k)/\nu$、$K(k)$ は積分上限を $\pi/2$ とした値である。物体成分が一周しても空間内の姿勢まで一周したとは限らない。

境界 $d=b$ では $A=B=1$、$\dot u=\nu(1-u^2)$ の上昇枝を積分すると $u=\tanh[\nu(t-t_0)]$。従って残り二成分は $\operatorname{sech}[\nu(t-t_0)]$ に比例する。中間軸へ到達するには無限時間を要し、両側から近づく周期は発散する。$k\to1$ で $K(k)$ の被積分関数が上限近くで $1/\cos v$ に近づくことからもこの発散が分かる。

:::example[最大慣性軸を囲む回転]
$I_1=2,I_2=3,I_3=4$（$\mathrm{kg\,m^2}$）、$l=6\,\mathrm{kg\,m^2/s}$ とし、$T=\frac{21}{4}\,\mathrm J$ を指定する。$u=0$ で $L_1,L_3$ をともに正にとると

$$
L_1(0)=l\sqrt{(d-c)/(a-c)}=\sqrt{6},\quad L_3(0)=\sqrt{30}.
$$

初期の $L_2$ の増加は無トルクでも起きる。これは空間では一定の角運動量を物体の軸に射影しているためである。動径に相当する変数 $u$ の二つの根は $A=\frac{1}{2}, B=\frac{5}{4}$。従って $k^2=\frac{2}{5},\ \nu=\frac{\sqrt{10}}{4}\,\mathrm{s^{-1}}$。$K$ の積分を評価すると周期は約 $8.994\,\mathrm s$。上の解では第三成分の符号が変わらず、第一成分が折返しのたびに符号を変える。
:::

:::example[中間軸へ近づく運動]
$I_1=2,I_2=3,I_3=4$（$\mathrm{kg\,m^2}$）、$l=6\,\mathrm{kg\,m^2/s}$ とし、$T=6\,\mathrm J$ を指定する。$u=0$ で $L_1,L_3$ をともに正にとると

$$
L_1(0)=l\sqrt{(d-c)/(a-c)}=2 \sqrt{3},\quad L_3(0)=2 \sqrt{6}.
$$

初期の $L_2$ の増加は無トルクでも起きる。これは空間では一定の角運動量を物体の軸に射影しているためである。動径に相当する変数 $u$ の二つの根は $A=1, B=1$。従って $k^2=1,\ \nu=\frac{\sqrt{2}}{2}\,\mathrm{s^{-1}}$。$u$ が零から $1/2$ に達する時間は $\log 3/(2\nu)=0.777\,\mathrm s$ である。周期解の公式に有限の数を代入せず、双曲線関数の枝を取る。この上昇枝は $t\to-\infty$ で $u\to-1$、$t\to+\infty$ で $u\to1$ となり、二つの中間軸回転をつなぐ。$L_1,L_3$ の符号を互いに逆に選べば下降枝となる。
:::

:::example[最小慣性軸を囲む回転]
$I_1=2,I_2=3,I_3=4$（$\mathrm{kg\,m^2}$）、$l=6\,\mathrm{kg\,m^2/s}$ とし、$T=\frac{15}{2}\,\mathrm J$ を指定する。$u=0$ で $L_1,L_3$ をともに正にとると

$$
L_1(0)=l\sqrt{(d-c)/(a-c)}=2 \sqrt{6},\quad L_3(0)=2 \sqrt{3}.
$$

初期の $L_2$ の増加は無トルクでも起きる。これは空間では一定の角運動量を物体の軸に射影しているためである。動径に相当する変数 $u$ の二つの根は $A=2, B=\frac{1}{2}$。従って $k^2=\frac{1}{4},\ \nu=1\,\mathrm{s^{-1}}$。$K$ の積分を評価すると周期は約 $6.743\,\mathrm s$。この枝では第一成分の符号が変わらず、第三成分が符号を変える。楕円関数の母数を同じ範囲に取るため、根の大小を交換する。
:::

### 物体成分の解から空間での姿勢を戻す

非対称こまの軸端を空間で追うには、楕円関数で求めた物体系の角運動量に一つの角の積分を加える。空間の鉛直を不変な $\vb L$ の方向に選ぶ。物体第三軸の傾きは $\cos\theta=L_3/L$ で決まるが、第三軸の方位 $\phi$ は依然未知である。物体成分で見た鉛直単位ベクトルは $(\sin\theta\sin\psi,\sin\theta\cos\psi,\cos\theta)$。従って $\sin\theta\ne0$ で $\sin\psi=L_1/(L\sin\theta),\ \cos\psi=L_2/(L\sin\theta)$。:eqref[eq-am-euler-angles] の第一、第二式の和を適切にとれば
$$
\dot\phi=\frac{L(L_1^2/I_1+L_2^2/I_2)}{L_1^2+L_2^2}.
$$
この式の右辺は各 $L_i(t)$ を知ればすべて求まる。一周期 $T_b$ でこれを積分した $\Delta\phi$ は一般に $2\pi$ の整数倍でない。従って物体の同じ点に塗った印は、一周期後にも空間で異なる位置に見える。

例として $(I_1,I_2,I_3)=(2,3,4)\,\mathrm{kg\,m^2}$、ある時刻の物体成分を $(L_1,L_2,L_3)=(2,0,2)\,\mathrm{kg\,m^2/s}$ とする。$L=2\sqrt2\,\mathrm{kg\,m^2/s}$ から $\theta=\pi/4,\psi=\pi/2$、$\dot\phi=\sqrt2\,\mathrm{s^{-1}}$。それでも $\omega_3=L_3/I_3=1/2\,\mathrm{s^{-1}}$ なので $\dot\psi=\omega_3-\dot\phi\cos\theta=-1/2\,\mathrm{s^{-1}}$。物体の自転角の変化率と第三軸への角速度の射影は異なる。周期全体の方位を楕円積分の標準形へ還元する場合は第三種の積分も現れる。自由こまの解のこの表示については :cite[landau-mechanics] を参照できる。

## 自由な対称こまと重いこま

### 二種類の歳差

$I_1=I_2=I,\ I_3=C$ と置く。
トルクのない場合 $\dot\omega_3=0$。
$w=\omega_1+i\omega_2$ とすれば :eqref[eq-am-euler-body] の最初の二式は $\dot w=i(C-I)\omega_3w/I$ だから、

$$
w(t)=w(0)e^{i(C-I)\omega_3t/I}.
$$

これは物体から見た角速度の回転である。
空間から見るには $\vb L=I\boldsymbol\omega+(C-I)\omega_3\vb e_3$ と解き直す。
$\dot{\vb e}_3=\boldsymbol\omega\times\vb e_3=(\vb L/I)\times\vb e_3$ より、対称軸は固定した $\vb L$ のまわりに角速度 $L/I$ で円錐を描く。
この運動を :term[歳差]{en="precession" yomi="さいさ"} という。
$\vb L$ との傾きは一定である。

### 重力で傾きも変わる場合

固定支点から対称軸の正方向に距離 $l$ の位置に重心がある。
$I,C$ は支点まわりの主慣性モーメントである。
$\theta$ を上向き鉛直から測り $h=Mgl>0$ とすれば $V=h\cos\theta$。
:eqref[eq-am-euler-angles] を使って

$$
L=\tfrac I2(\dot\theta^2+\dot\phi^2\sin^2\theta)
+\tfrac C2(\dot\psi+\dot\phi\cos\theta)^2-h\cos\theta .
$$

$\phi,\psi$ は循環座標である。
$b=p_\psi=C\omega_3$ は第三軸への角運動量射影、$a=p_\phi=I\sin^2\theta\,\dot\phi+b\cos\theta$ は鉛直軸への射影である。
これらをエネルギーに入れて

$$
\begin{aligned}
E'&=E-b^2/(2C)=\tfrac I2\dot\theta^2+U(\theta),\\
U(\theta)&=\frac{(a-b\cos\theta)^2}{2I\sin^2\theta}+h\cos\theta,\\
\dot\phi&=\frac{a-b\cos\theta}{I\sin^2\theta},\qquad
\dot\psi=b/C-\dot\phi\cos\theta .
\end{aligned}
\label{eq-am-top-potential}
$$

$u=\cos\theta$ と置くと $I^2\dot u^2=f(u)=2I(E'-hu)(1-u^2)-(a-bu)^2$。
物理域 $-1\le u\le1$ で $f\ge0$ の区間を選ぶ。
隣り合う単純根の間で傾きが往復する運動を :term[章動]{en="nutation" yomi="しょうどう"} という。
周期は $2I\int_{u_-}^{u_+}\dd u/\sqrt{f(u)}$。
この積分と $\dot\phi,\dot\psi$ の積分で姿勢を再構成できる。
$a-bu$ が区間内で符号を変えると、歳差は途中で逆行する。
極 $u=\pm1$ に達するには $a=\pm b$ が必要である。

```tikz caption="固定支点の対称こま。傾角は上向き鉛直から測り、重心には下向きの重力が働く。"
\coordinate (O) at (0,0);
\coordinate (C) at (1.35,1.8);
\draw[->] (O) -- (0,2.8) node[above] {$z$};
\draw[thick,->] (O) -- (1.95,2.6) node[above] {$\mathbf e_3$};
\draw[->,thick] (C) -- (1.35,.65) node[right] {$M\mathbf g$};
\draw (0,.7) arc[start angle=90,end angle=53.13,radius=.7];
\node at (.3,.95) {$\theta$};
\node[right] at (.73,.85) {$l$};
\fill (C) circle (2pt) node[right] {$C$};
\fill (O) circle (1.7pt) node[below] {$O$};
```

### 定常歳差と直立回転

$\theta$ が一定で極以外なら、$\theta$ の方程式を $\sin\theta$ で割ると

$$
I\cos\theta\,\Omega_p^2-b\Omega_p+h=0,\qquad \Omega_p=\dot\phi .
\label{eq-am-steady-precession}
$$

$\cos\theta>0$ では $b^2\ge4Ih\cos\theta$ が必要である。
高速自転で $b>0$ なら、遅い根は $\Omega_p=2h/(b+\sqrt{b^2-4Ih\cos\theta})\simeq h/b$、速い根はおよそ $b/(I\cos\theta)$。
例えば $I=0.02\,\mathrm{kg\,m^2}$、$h=0.1\,\mathrm J$、$b=0.4\,\mathrm{kg\,m^2/s}$、$\cos\theta=1/2$ なら遅い根は約 $0.252\,\mathrm{s^{-1}}$ である。

直立 $\theta=0$ では節線を使えないが、$a=b$ と置いた有効ポテンシャルの極限はとれる。
$(1-\cos\theta)^2/\sin^2\theta=\tan^2(\theta/2)$ より

$$
U=h+\left(\frac{b^2}{8I}-\frac h2\right)\theta^2+O(\theta^4).
$$

従ってこの保存量面での直立は $b^2>4Ih$ で局所極小、$b^2<4Ih$ で極大になる。
境界では二次の判定は使えない。極で正則な座標を使うと、直立の線形応答を確認できる。対称軸の単位ベクトルを $\vb n$ とすれば $\vb L=b\vb n+I\vb n\times\dot{\vb n}$。その微分を重力トルク $h\vb e_z\times\vb n$ と等置する。$\vb n=\vb e_z+(x,y,0)+O(x^2+y^2)$ と置き、両辺に $-\vb e_z\times$ を作用させると $Z=x+iy$ は

$$
I\ddot Z-ib\dot Z-hZ=0
$$

に従う。時間指数は $i[b\pm\sqrt{b^2-4Ih}]/(2I)$ である。$b^2<4Ih$ なら成長指数があり、$b^2>4Ih$ なら異なる二つの振動成分がある。等号では重根となり、一般の線形解は $t e^{ibt/(2I)}$ を含む。この場合は有限の摂動に対する安定性を線形解から決められない。実際、$a=b$ の面に制限した有効ポテンシャルは、後の演習で求めるように四次の極小をもつ。

下向きの回転 $\vb n=-\vb e_z$ の近くでは、同じ操作により $I\ddot Z+ib\dot Z+h Z=0$ となる。判別式は $-b^2-4Ih<0$ で、どの $b$ でも指数は相異なる純虚数である。$a=-b$ の有効ポテンシャルも $\theta=\pi$ で $U=-h+[b^2/(8I)+h/2](\pi-\theta)^2+O((\pi-\theta)^4)$ と展開でき、極小である。

### 章動するこまを初期値から再構成する

転回点を図で読む代わりに、根の位置を符号だけで判定できる数値を使ってみよう。
以下は $I,C$ を $\mathrm{kg\,m^2}$、$h,E'$ を J、運動量を $\mathrm{kg\,m^2/s}$ で測った値とする。

$$
I=C=1,\quad h=1,\quad a=0,\quad b=1,\quad E'=1.
$$

例えば $\theta=\pi/2,\ \dot\phi=0,\ \dot\psi=1,\ \dot\theta=-\sqrt2$ がこれを実現する。
$f(u)$ は

$$
f(u)=2u^3-3u^2-2u+2
$$

という三次式である。
転回点はこの式そのものの符号で調べればよい。
$f(-1)=-1,\ f(0)=2,\ f(1)=-1$ かつ十分大きな正の $u$ で $f>0$ なので、$(-1,0)$、$(0,1)$、$(1,\infty)$ に根が一つずつある。
最初の二つだけが物理域の転回点である。
初めは $\dot u=-\sin\theta\,\dot\theta>0$。
正の転回点まで進んで折り返し、次に負の転回点まで進む。

この途中で $\dot\phi=-u/(1-u^2)$ は $u=0$ を横切るたびに符号を変える。
方位が逆行するのはエネルギーの符号が変わるためではない。
固定した鉛直運動量 $a=0$ と、自転運動量の射影 $bu$ がつり合うように歳差の寄与が変化するからである。
方位の正味の進みは、一章動周期に対して

$$
\Delta\phi=2\int_{u_-}^{u_+}
\frac{a-bu}{(1-u^2)\sqrt{f(u)}}\,\dd u
$$

となる。
根における特異性は単純根なら積分可能である。

### 高速歳差の近似をトルクから点検する

大きな自転に比べて軸の動きが遅いときには、$\vb L\simeq b\vb e_3$ となる。
重心に働く重力の支点まわりのトルクは $\vb N=h\vb e_z\times\vb e_3$ であり、鉛直まわりの定常歳差なら $\dot{\vb L}\simeq b\Omega_p\vb e_z\times\vb e_3$。
従って $\Omega_p\simeq h/b$。
二次式の速い根を捨ててよいのは、$\vb L$ がほぼ自転だけで決まる初期条件を与えた場合である。
角速度を大きくしただけでは初期の横運動量が消える保証はないため、一般にはこの歳差に章動が重なる。

### 軸を傾けて自転だけを与えたこま

$0<\theta_0<\pi/2$ で $\dot\phi(0)=\dot\theta(0)=0$ とし、$\dot\psi(0)=b/C$ だけを与える。
この初期状態は定常歳差ではない。
保存量は

$$
u_0=\cos\theta_0,\qquad a=bu_0,\qquad E'=hu_0 .
$$

これを三次式に代入すると、初期の根を陽にくくれて

$$
f(u)=(u_0-u)\,[2Ih(1-u^2)-b^2(u_0-u)] .
$$

$u_0$ で角速度が零でも $\theta$ 方程式は $I\ddot\theta(0)=h\sin\theta_0>0$ を与えるため、まず倒れる向きに傾いて $u$ は減少する。
$u$ が少し減ると $\dot\phi=b(u_0-u)/[I(1-u^2)]$ が生じる。
すなわち歳差は初期から強制されたのではなく、重力のトルクから生じる。

$b^2\gg Ih$ なら章動の幅は小さい。
$\delta=u_0-u$ と置いて第二因子を一次まで展開すると

$$
2Ih(1-u_0^2)+(4Ihu_0-b^2)\delta+O(Ih\delta^2)
$$

となる。
従って次の転回点までの幅は主要項で $\delta_{\max}\simeq2Ih\sin^2\theta_0/b^2$、角度にして $\Delta\theta\simeq2Ih\sin\theta_0/b^2$。
高速自転にすると自転の逆数の二乗で章動幅が小さくなる。

転回点近くの一周期をさらに追うには、$f(u)\simeq b^2\delta(\delta_{\max}-\delta)$ と近似する。
初期で $\delta=0$ に合う解は

$$
\delta(t)\simeq\frac{\delta_{\max}}2[1-\cos(|b|t/I)].
$$

この式を $\dot\phi$ に代入して一周期平均すれば $\langle\dot\phi\rangle\simeq h/b$。
高速こまの遅い歳差は、初期条件によっては微小だが速い章動を平均して初めて見える。

:::example[定常歳差を用意する初期条件]
$I=0.10\,\mathrm{kg\,m^2},\ C=0.05\,\mathrm{kg\,m^2},\ h=0.20\,\mathrm J$、
$\theta=\pi/3,\ b=0.30\,\mathrm{kg\,m^2/s}$ とする。
定常歳差条件は SI の数値で

$$
0.05\Omega_p^2-0.30\Omega_p+0.20=0 .
$$

従って $\Omega_p=3\pm\sqrt5$（毎秒）。
どちらも同じ傾角・自転成分 $\omega_3=6\,\mathrm{s^{-1}}$ を持つが、必要な鉛直角運動量は

$$
a=I\sin^2\theta\,\Omega_p+b\cos\theta=0.075\Omega_p+0.15
$$

と異なる。
用意すべき Euler 角の初速も $\dot\psi=6-\Omega_p/2$ で区別される。
単に $\dot\psi=6$ として傾けて放すと、いずれの定常解の初期条件にも一致しない。
:::

### 固定支点と重心のどちらで力を数えるか

重いこまの支点は静止しているため、その反力は支点まわりのトルクを持たず、仕事もしない。
支点まわりでは重力トルク一つで回転運動を解ける。
一方、重心まわりでは重力の腕は零で、回転を変えるトルクは支点反力から来る。
二つの表示が一致することは

$$
\vb L_O=\vb L_C+\vb r_C\times M\dot{\vb r}_C
$$

を微分すれば分かる。
$\dot{\vb L}_C=(-\vb r_C)\times\vb R$、$M\ddot{\vb r}_C=M\vb g+\vb R$ を使うと $\dot{\vb L}_O=\vb r_C\times M\vb g$ で、
反力 $\vb R$ は消える。
慣性モーメントだけを重心のものに替えて重力トルクは支点のものを使うと、この相殺を壊す。

### 支点にかかる荷重を歳差から求める

定常歳差では重心は水平半径 $l\sin\theta$、高さ $l\cos\theta$ の円を描く。
鉛直軸から外へ向く水平単位ベクトルを $\vb e_r$ とすれば、加速度は $-l\sin\theta\,\Omega_p^2\vb e_r$ である。
従ってこまが支点から受ける反力は

$$
\vb R=Mg\vb e_z-Ml\sin\theta\,\Omega_p^2\vb e_r .
$$

鉛直成分は重力とつり合い、水平成分が重心の向心力を与える。
支点がこまから受ける力はこの逆向きである。
遅い歳差と速い歳差で同じ傾角が実現しても、水平荷重は $\Omega_p^2$ に比例するので大きく違う。

:::example[自由な対称こまの姿勢を Euler 角で戻す]
トルクがない場合に空間の $z$ 軸を保存する $\vb L$ に合わせる。
対称軸を $\vb L$ と平行でない向きに選ぶと $\dot\theta=0,\ \dot\phi=L/I$。
また $L\cos\theta=C\omega_3$ なので、

$$
\dot\psi=\omega_3-\dot\phi\cos\theta=(1-C/I)\omega_3 .
$$

従って $\dot\psi$ は物体から見た角速度の横成分の位相速度 $(C/I-1)\omega_3$ と符号が逆になる。
節線を基準に測る自転角と、物体軸から測る横成分の位相では観測している相手が逆だからである。
$C=2I,\omega_3=1\,\mathrm{s^{-1}}$ なら $\dot\psi=-1\,\mathrm{s^{-1}}$ だが、第三軸に沿う物理的な角速度は正であり、
自転が反転したわけではない。
:::

### 撃力の作用点と支点の衝撃

こまの軸受を含め、固定支点は力積を伝達できる。短時間の撃力なら、衝突中の姿勢変化と重力の力積を無視して積分形の運動方程式を使う。平面剛体に固定支点 $O$ があり、重心は支点から距離 $a$ にある。支点と重心を結ぶ線上で距離 $b$ の点を、この線に垂直に力積 $J$ で打つ。支点の反力積を同じ向きを正として $S$ とすれば
$$
I_O\Delta\omega=bJ,\qquad M a\Delta\omega=J+S,\qquad
S=J(Mab/I_O-1).
$$
従って $b_*=I_O/(Ma)$ に当てれば支点の反力積が零になる。この点を :term[撃心]{en="center of percussion" yomi="げきしん"} という。$a>0$、$I_C>0$ の場合は $b_*=a+I_C/(Ma)>a$ である。この条件は「撃った直後に支点が静止する」だけからは得られない。支点はどこを撃っても固定されており、その仕事が零でも一般に力積は零でない。

物理振り子の微小振動方程式 $I_O\ddot\theta=-Mga\theta$ の換算長も同じ $b_*$ になる。重心から見て反対側の距離 $a'=b_*-a=I_C/(Ma)$ に新しい支点を取ると $a'+I_C/(Ma')=a'+a=b_*$。従って二つの支点で小振動周期が等しい。この交換性には同じ剛体を用いること、二軸が平行なこと、小振幅であることが必要である。

:::example[端で支えた細棒を重心で打つ]
支点軸に関するデータを $M=3\,\mathrm{kg},a=1\,\mathrm m,I_C=1\,\mathrm{kg\,m^2}$ とする。支点から $b=1\,\mathrm m$ の点に $J=4\,\mathrm{N\,s}$ を与える。平行軸の定理で $I_O=4\,\mathrm{kg\,m^2}$ だから、角速度の跳びと支点の反力積は

$$
\Delta\omega=\frac{bJ}{I_O}=1\,\mathrm{s^{-1}},\qquad S=Ma\Delta\omega-J=-1\,\mathrm{N\,s}.
$$

力積の符号は打撃方向を正としている。撃心までの距離は $\frac{4}{3}\,\mathrm m$。初め静止していれば、衝突直後の運動エネルギーは $I_O(\Delta\omega)^2/2=2\,\mathrm J$ である。
:::

:::example[撃心を打つ薄い円板]
支点軸に関するデータを $M=2\,\mathrm{kg},a=1\,\mathrm m,I_C=1\,\mathrm{kg\,m^2}$ とする。支点から $b=\frac{3}{2}\,\mathrm m$ の点に $J=3\,\mathrm{N\,s}$ を与える。平行軸の定理で $I_O=3\,\mathrm{kg\,m^2}$ だから、角速度の跳びと支点の反力積は

$$
\Delta\omega=\frac{bJ}{I_O}=\frac{3}{2}\,\mathrm{s^{-1}},\qquad S=Ma\Delta\omega-J=0\,\mathrm{N\,s}.
$$

撃心までの距離は $\frac{3}{2}\,\mathrm m$。初め静止していれば、衝突直後の運動エネルギーは $I_O(\Delta\omega)^2/2=\frac{27}{8}\,\mathrm J$ である。
:::

:::example[重心より支点寄りに当てる振り子]
支点軸に関するデータを $M=4\,\mathrm{kg},a=\frac{1}{2}\,\mathrm m,I_C=1\,\mathrm{kg\,m^2}$ とする。支点から $b=\frac{1}{4}\,\mathrm m$ の点に $J=2\,\mathrm{N\,s}$ を与える。平行軸の定理で $I_O=2\,\mathrm{kg\,m^2}$ だから、角速度の跳びと支点の反力積は

$$
\Delta\omega=\frac{bJ}{I_O}=\frac{1}{4}\,\mathrm{s^{-1}},\qquad S=Ma\Delta\omega-J=- \frac{3}{2}\,\mathrm{N\,s}.
$$

撃心までの距離は $1\,\mathrm m$。初め静止していれば、衝突直後の運動エネルギーは $I_O(\Delta\omega)^2/2=\frac{1}{16}\,\mathrm J$ である。
:::

### 一章動で進む方位と軸端の軌跡

軸の先端に光点を置いて球面に軌跡を記録すると、同じ二つの緯度の間を往復しても、方位が単調に進む線、輪を描いて逆行する線、尖点をもつ線が現れる。$b\ne0$ とし、方位速度が零になる候補 $u_0=a/b$ をまず求める。$u_0$ が許容区間の外なら方位は単調、内部なら折返す。単純転回点と一致すれば、そこで $\dot\theta$ と $\dot\phi$ がともに零になる。根の近くで $u-u_0$ は経過時間の二乗、$\dot\phi$ も二乗に比例するため、方位のずれは三乗に比例する。球面の極以外では緯度のずれの三乗と方位のずれの二乗が同程度になり、尖点を作る。

数値積分では二つの物理根を $A<B$、第三根を $D>B$ として $u=A+(B-A)\sin^2 v$ と置く。すると正の平方根を取る半周期で
$$
\frac{\dd u}{\sqrt{f(u)}}=\frac{2\,\dd v}{\sqrt{2Ih(D-u)}},\qquad 0\le v\le\pi/2 .
$$
根を始点にしたまま $1/\sqrt{f}$ を評価する必要がなくなる。往復時間と方位の増分にはそれぞれ $2I$ と $2(a-bu)/(1-u^2)$ を掛けて積分する。方位の増分を章動周期で割った値は平均歳差速度であり、ある時点で観測した角速度とは異なる。

:::example[単調に進む方位]
数値は $I=C=1\,\mathrm{kg\,m^2},h=1\,\mathrm J$ とし、$a=\frac{3}{2},b=1$（$\mathrm{kg\,m^2/s}$）、$E'=2\,\mathrm J$ とする。この条件での三次式は $f(u)=\frac{8 u^{3} - 20 u^{2} + 4 u + 7}{4}$。根を小さい順に求めると $-0.46540,\ 0.91852,\ 2.04688$ となる。最初の二つの間では $f\ge0$ であり、第三根は球面上の緯度を表さない。ここで $u_0=\frac{3}{2}$ を二根と比較する。方位速度の零点が物理範囲の外にあり、緯度が折り返しても方位の向きは変わらない。上の置換を使って時間と方位を別々に積分すれば、一章動の周期は約 $3.3900\,\mathrm s$、一章動の正味の方位は約 $6.7233$ ラジアンである。したがって平均歳差速度は約 $1.9832\,\mathrm{s^{-1}}$。単に二つの転回点での角速度を算術平均すると、各緯度に滞在する時間の違いを落とす。
:::

:::example[途中で逆行する方位]
数値は $I=C=1\,\mathrm{kg\,m^2},h=1\,\mathrm J$ とし、$a=\frac{1}{3},b=1$（$\mathrm{kg\,m^2/s}$）、$E'=2\,\mathrm J$ とする。この条件での三次式は $f(u)=\frac{18 u^{3} - 45 u^{2} - 12 u + 35}{9}$。根を小さい順に求めると $-0.86556,\ 0.91773,\ 2.44784$ となる。最初の二つの間では $f\ge0$ であり、第三根は球面上の緯度を表さない。ここで $u_0=\frac{1}{3}$ を二根と比較する。積分する向きを動径の枝だけで決めると、負の方位の寄与を正に数えてしまう。軸端の輪を作る区間では符号を保って足す。上の置換を使って時間と方位を別々に積分すれば、一章動の周期は約 $2.9335\,\mathrm s$、一章動の正味の方位は約 $0.6143$ ラジアンである。したがって平均歳差速度は約 $0.2094\,\mathrm{s^{-1}}$。
:::

:::example[ゆっくり動き出す尖点]
数値は $I=C=1\,\mathrm{kg\,m^2},h=1\,\mathrm J$ とし、$a=\frac{1}{2},b=1$（$\mathrm{kg\,m^2/s}$）、$E'=\frac{1}{2}\,\mathrm J$ とする。この条件での三次式は $f(u)=\frac{\left(2 u - 1\right) \left(4 u^{2} - 2 u - 3\right)}{4}$。根を小さい順に求めると $-0.65139,\ 0.50000,\ 1.15139$ となる。最初の二つの間では $f\ge0$ であり、第三根は球面上の緯度を表さない。ここで $u_0=\frac{1}{2}$ を二根と比較する。上の転回点では二つの速度が同時に消える。重力のトルクは零でないため、軸端はそこで止まり続けない。上の置換を使って時間と方位を別々に積分すれば、一章動の周期は約 $4.1999\,\mathrm s$、一章動の正味の方位は約 $2.7524$ ラジアンである。したがって平均歳差速度は約 $0.6554\,\mathrm{s^{-1}}$。
:::

## 演習問題

::::exercise{level=basic}
一様な細棒 $M=3\,\mathrm{kg},l=2\,\mathrm m$ の中央と端で棒に垂直な軸の慣性モーメントを求めよ。

:::solution
中央で $Ml^2/12=1\,\mathrm{kg\,m^2}$、端で平行軸項 $M(l/2)^2=3\,\mathrm{kg\,m^2}$ を加え $4\,\mathrm{kg\,m^2}$ である。
:::
::::

::::exercise{level=basic}
円板がすべらず速さ $v$ で転がるとき、全運動エネルギーを求めよ。

:::solution
$T=Mv^2/2+(Ma^2/2)(v/a)^2/2=3Mv^2/4$。
回転エネルギーが全体の三分の一である。
:::
::::

::::exercise{level=basic}
$I=\diag(2,3,4)\,\mathrm{kg\,m^2}$、$\omega=(1,0,2)\,\mathrm{s^{-1}}$ で $\vb L,T$ を求めよ。

:::solution
$\vb L=I\omega=(2,0,8)\,\mathrm{kg\,m^2/s}$ である。
$T=\omega\cdot\vb L/2=(2+16)/2=9\,\mathrm J$。
:::
::::

::::exercise{level=basic}
$\theta=\pi/2,\psi=0$ の瞬間の物体角速度を Euler 角の速度で表せ。

:::solution
鉛直軸は物体第二軸と一致し、節線は第一軸である。
従って $(\omega_1,\omega_2,\omega_3)=(\dot\theta,\dot\phi,\dot\psi)$ となる。
:::
::::

::::exercise{level=advanced}
自由対称こまで $C=3I/2,\omega_3=8\,\mathrm{s^{-1}}$ のとき横成分の物体系での周期を求めよ。

:::solution
$w$ の位相速度は $(C-I)\omega_3/I=4\,\mathrm{s^{-1}}$。
横成分が零でなければ周期は $2\pi/4=\pi/2\,\mathrm s$ である。
:::
::::

::::exercise{level=advanced}
三つの主慣性モーメントは $I_1+I_2\ge I_3$ を満たすことを示し、等号の条件を述べよ。

:::solution
物体主軸で $I_1+I_2-I_3=\int[(y^2+z^2)+(z^2+x^2)-(x^2+y^2)]\dd m=2\int z^2\dd m\ge0$。
質量が支点を通る $z=0$ 平面に集中するときに限り等号である。
:::
::::

::::exercise{level=advanced}
重い対称こまを $\theta=\pi/2$ で定常歳差させる条件を求めよ。
$h>0$ とする。

:::solution
二次式の二次項が消え、$-b\Omega_p+h=0$ となる。
従って $b\ne0,\ \Omega_p=h/b$。
さらに $a=I\Omega_p$ を初期の鉛直角運動量に選べば一定傾角が実現する。
:::
::::

::::exercise{level=advanced}
上向き直立で $a=b,\ b^2=4Ih$ のとき、有効ポテンシャルの最初の非定数項を求めよ。

:::solution
$U=2h\tan^2(\theta/2)+h\cos\theta$。
$\tan^2(\theta/2)=\theta^2/4+\theta^4/24+O(\theta^6)$ と $\cos\theta=1-\theta^2/2+\theta^4/24+O(\theta^6)$ を加えると $U=h+h\theta^4/8+O(\theta^6)$。
この保存量面で平坦な極小となる。
:::
::::

::::exercise{level=advanced}
細い棒二本を原点で直角に接合した。
一本は $0\le x\le l$、他方は $0\le y\le l$ にあり、各棒の質量は $m$ である。
原点まわりと重心まわりの主軸を求めよ。

:::solution
原点まわりでは各点で $xy=0$ なので面内の非対角成分は零。
$I_{xx}=I_{yy}=ml^2/3,\ I_{zz}=2ml^2/3$ であり、面内の任意の軸が主軸である。
重心は $(l/4,l/4,0)$。
平行軸項を引けば重心まわりの面内ブロックは

$$
ml^2\begin{pmatrix}5/24&1/8\\1/8&5/24\end{pmatrix}.
$$

対角線 $(1,1)$ が固有値 $ml^2/3$、$(1,-1)$ が $ml^2/12$ である。
面外は薄板の関係から $5ml^2/12$。
支点を移すと同じ剛体でも主軸の縮退が解ける。
:::
::::

::::exercise{level=advanced}
重い対称こまの定常歳差の傾き $\theta_*$ を少し変える。
ただし保存量 $a,b$ を固定し、$0<\theta_*<\pi$ とする。
微小章動の角振動数を定常歳差速度 $\Omega_p$ で表せ。

:::solution
$u=\cos\theta$、$\Omega(u)=(a-bu)/[I(1-u^2)]$ とする。
有効ポテンシャルの微分は $U_u=h-b\Omega+Iu\Omega^2$ である。
平衡では $U_u=0$。
$\Omega'=(2Iu\Omega-b)/[I(1-u^2)]$ より

$$
U_{uu}=I\Omega^2+\frac{(2Iu\Omega-b)^2}{I(1-u^2)} .
$$

$U_{\theta\theta}=(1-u^2)U_{uu}-uU_u$ なので、$\theta_*$ で

$$
\omega_{\rm nut}^2=\frac{U_{\theta\theta}}I
=(1-u_*^2)\Omega_p^2+(2u_*\Omega_p-b/I)^2 .
$$

$h>0$ の定常歳差では $\Omega_p\ne0$ だから正である。
これは $a,b$ を固定した傾角の安定性であり、歳差速度を外部から強制して一定に保つ問題とは条件が異なる。
:::
::::

::::exercise{level=advanced}
一定の物体角速度 $\omega_3=\Omega$ だけをもつ剛体が、初め単位姿勢から動く。四元数を時間について解き、一周後の符号と観測した姿勢を比較せよ。

:::solution
$\dot a=a(0,\Omega\vb e_3)/2$ より $\dot a_0=-\Omega a_3/2,\dot a_3=\Omega a_0/2$、他の二成分は零のままである。従って $a=(\cos(\Omega t/2),0,0,\sin(\Omega t/2))$。$\Omega t=2\pi$ で $a=(-1,0,0,0)$ だが、挟む積はすべてのベクトルを元に戻す。姿勢の時系列を滑らかに補間するとき、各時刻でスカラー部を強制的に正にすると半回転を通る際に不連続を入れてしまう。
:::
::::

## 確認テスト

::::quiz{answer="6" tol="1e-6"}
非共線な点を含む自由剛体の配位の自由度はいくつか。

:::explanation
重心に3、姿勢に3の独立座標があるから合計6である。
:::
::::

::::quiz
自由な非対称剛体で不安定な定常回転軸はどれか。

- [ ] 最大主慣性モーメントの軸
- [x] 中間主慣性モーメントの軸
- [ ] 最小主慣性モーメントの軸

:::explanation
中間軸での摂動には $\delta\ddot\omega_1=c\,\delta\omega_1$、$c>0$ の成長解がある。
:::
::::

::::quiz{answer="3" tol="1e-6"}
$L=6\,\mathrm{kg\,m^2/s}$、$I=2\,\mathrm{kg\,m^2}$ の自由対称こまの第三軸が角運動量に傾いている。
空間歳差の角速度を毎秒で答えよ。

:::explanation
$\dot{\vb e}_3=(\vb L/I)\times\vb e_3$ より $L/I=3\,\mathrm{s^{-1}}$ である。
:::
::::

::::quiz
重い対称こまで支点まわりのどの角運動量成分が保存するか。
すべて選べ。

- [x] 鉛直への射影
- [x] 対称軸への射影
- [ ] 任意の空間水平軸への射影

:::explanation
$\phi,\psi$ が循環座標であり、それぞれの共役運動量は鉛直射影 $a$ と対称軸射影 $b$ である。
重力トルクの水平成分は一般には零でない。
:::
::::

::::quiz
剛体が慣性主軸まわりに自由回転するとき、常に成り立つものを選べ。
慣性モーメントは正とする。

- [x] 角速度と角運動量が平行になる
- [ ] 中間軸の回転でも小さなずれは必ず振動する
- [ ] どの主軸でも同じ角速度ならエネルギーが等しい

:::explanation
主軸では $I\omega=I_i\omega$ だから二つのベクトルは平行になる。
三軸の慣性が異なる場合、中間軸は不安定であり、同じ $\omega$ でもエネルギー $I_i\omega^2/2$ は軸に依存する。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
傾角一定の重いこまの重心が半径 $0.25\,\mathrm m$ の水平円を角速度 $2\,\mathrm{s^{-1}}$ で回る。
質量 $0.5\,\mathrm{kg}$ のこまが支点から受ける水平力の大きさを N で答えよ。

:::explanation
向心加速度は $0.25\times2^2=1\,\mathrm{m/s^2}$ である。
重力には水平成分がないので反力の水平成分が $0.5\times1=0.5\,\mathrm N$ を担う。
:::
::::

---
title: 光学系と収差
description: 厚いレンズの基準面、開口と視野、光学機器、単色・色収差と点像の評価を扱う。
notation:
  - symbol: '\Phi'
    meaning: 光学系の換算屈折力
  - symbol: "H,H'"
    meaning: 物体側および像側の主点
  - symbol: 'W'
    meaning: 波面収差
  - symbol: 'S'
    meaning: Strehl 比
---

## この章の目標

- 厚いレンズの行列から基準面と焦点を求める。
- 絞りを光学系の両側から見て、明るさと視野を区別する。
- 回転対称系の収差を瞳と視野の次数で分類する。
- 波面誤差から焦点強度の低下を見積もる。

## 厚いレンズと主点・節点

### 主面への分解
厚さ $d$、内部屈折率 $n_\ell$、二面の屈折力 $\Phi_1,\Phi_2$ のレンズを考える。$e=d/n_\ell$ とすると :eqref[eq-optics-ray-matrix] の積は
$$
L(\Phi_2)T(e)L(\Phi_1)=
\begin{pmatrix}1-e\Phi_1&e\\-\Phi&1-e\Phi_2\end{pmatrix},
\qquad\Phi=\Phi_1+\Phi_2-e\Phi_1\Phi_2 .
\label{eq-optics-thick-power}
$$
一般に $\det M=1,\ C=-\Phi\ne0$ なら $M=T(e_2)L(\Phi)T(e_1)$ と分解できる。実際、積の対角成分を比較して $e_1=(1-D)/\Phi$、$e_2=(1-A)/\Phi$ を得る。残る右上成分は行列式の条件で一致する。入力頂点から $n_1e_1$ 進んだ光軸上の点を $H$、出力頂点から $n_2 e_2$ 戻った点を $H'$ とする。両点を通る横断面が :term[主面]{en="principal plane" yomi="しゅめん"} である。主面間の写像 $L(\Phi)$ は高さを保つ。従って、全系は主面を基準に $n_1/a+n_2/b=\Phi$ で結像する。

### 焦点と節点
入射が平行なら像側焦点は $H'$ の右 $n_2/\Phi$ にある。出射が平行なら物体点は $H$ の左 $n_1/\Phi$ にある。:term[節点]{en="nodal point" yomi="せってん"} は、前側の点へ向けた光線が後側の点から同じ角度で出る一対の点である。前側節点を $H$ の右 $\delta$ に置くと、角 $\theta$ の光線の主面での高さは $-\delta\theta$。出射換算角は $n_1\theta+\Phi\delta\theta$ なので、角を保つには $\delta=(n_2-n_1)/\Phi$ である。像側の角も $\theta$ だから、高さを保つ主面間の写像の後で距離 $\delta$ 進めば高さは $-\delta\theta+\delta\theta=0$ になる。従って後側節点も $H'$ の右 $\delta$ にある。同じ両外側媒質では主点と節点が一致する。

:::example[焦点距離と後側作動距離]
$\Phi_1=\Phi_2=5\,\mathrm{m^{-1}},e=0.02\,\mathrm m$、外側は空気とする。$\Phi=9.5\,\mathrm{m^{-1}}$、両 $e_i=0.10/9.5\,\mathrm m$。像側焦点と出力頂点の距離は $(1-0.10)/9.5\simeq94.7\,\mathrm{mm}$ である。
:::
:::example[屈折力がゼロの系]
$C=0$ なら平行光は平行に出る。$A$ が有限でも焦点距離 $n_2/\Phi$ は有限に取れず、上の主面分解は使えない。このときは行列自体の $A,D,B$ で倍率と位置を追う。
:::

:::example[平凸レンズを裏返す]
外部は空気、$n_\ell=1.50,d=6\,\mathrm{mm},|R|=50\,\mathrm{mm}$ とする。曲面一つの力は $10\,\mathrm{m^{-1}}$、平面の力はゼロ、$e=4\,\mathrm{mm}$。曲面を先にすれば $A=1-0.004(10)=0.96$ で後頂点からの焦点距離は $A/\Phi=96\,\mathrm{mm}$。平面を先にすれば $A=1$ で $100\,\mathrm{mm}$ になる。主点からの像側焦点距離は両者とも $100\,\mathrm{mm}$ である。
:::

### 頂点基準で像と倍率を求める

二主点の位置が測れなくても、頂点での光線行列 $M$ を知れば結像を予測できる。以下では外側を空気とし、物面を入射頂点の左に距離 $a$、像面を出射頂点の右に符号付き距離 $b$ に置く。全行列は $T(b)MT(a)$。その右上は $Aa+B+b(Ca+D)$ である。これがゼロなら物面の高さを固定して出発角を変えても像高は変わらない。従って $b=-(Aa+B)/(Ca+D)$。左上成分 $A+bC$ にこの $b$ を代入すると $(AD-BC)/(Ca+D)=1/(Ca+D)$ である。分母がゼロなら物面は前焦点面にあり、有限の像面はない。有限共役が焦点近くにある場合には、主面の小さな位置誤差も像距離に大きく効く。

:::example[厚い対称レンズ]
二面が同じ正の屈折力を持つレンズに、有限距離の物体を置く。

入力は第一面の屈折力 $\Phi_1=5\,\mathrm{m^{-1}}$、第二面の屈折力 $\Phi_2=5\,\mathrm{m^{-1}}$、換算厚さ $e=\frac{1}{50}\,\mathrm m$、物面と入射頂点の距離 $a=\frac{1}{4}\,\mathrm m$ である。

第一面での角変化を内部伝播の高さへ換算する。

$$
\begin{aligned}
A&=- e \Phi_1 + 1,\\
A&= - \left(\frac{1}{50}\right) \cdot \left(5\right) + 1=\frac{9}{10}.
\end{aligned}
$$

高さゼロから入れると第二面への高さは $e p_{\rm ray}$ である。

$$
\begin{aligned}
D&=- e \Phi_2 + 1,\\
D&= - \left(\frac{1}{50}\right) \cdot \left(5\right) + 1=\frac{9}{10}.
\end{aligned}
$$

二面の屈折力に距離の補正を加える。

$$
\begin{aligned}
\Phi&=- e \Phi_1 \Phi_2 + \Phi_1 + \Phi_2,\\
\Phi/(\mathrm{m^{-1}})&= - \left(\frac{1}{50}\right) \cdot \left(5\right) \cdot \left(5\right) + \left(5\right) + \left(5\right)=\frac{19}{2}.
\end{aligned}
$$

入射が平行な場合に出力から軸へ達する距離を取る。

$$
\begin{aligned}
b_{\infty}&=\frac{A}{\Phi},\\
b_{\infty}/(\mathrm m)&= \frac{\left(\frac{9}{10}\right)}{\left(\frac{19}{2}\right)}=\frac{9}{95}.
\end{aligned}
$$

出射角に物面の出発角の情報が残っても、像面の高さでは消えるように伝播距離を選ぶ。

$$
\begin{aligned}
b&=\frac{- A a - e}{D - \Phi a},\\
b/(\mathrm m)&= \frac{- \left(\frac{9}{10}\right) \cdot \left(\frac{1}{4}\right) - \left(\frac{1}{50}\right)}{\left(\frac{9}{10}\right) - \left(\frac{19}{2}\right) \cdot \left(\frac{1}{4}\right)}=\frac{49}{295}.
\end{aligned}
$$

行列式が一なので同じ分母の逆数が高さの倍率になる。

$$
\begin{aligned}
M&=\frac{1}{D - \Phi a},\\
M&= \frac{1}{\left(\frac{9}{10}\right) - \left(\frac{19}{2}\right) \cdot \left(\frac{1}{4}\right)}=- \frac{40}{59}.
\end{aligned}
$$

後焦点と像面を区別すると、無限遠で合焦した位置からセンサーをどれだけ動かすか分かる。
:::

:::example[後面が平面のレンズ]
曲面が最初に来る平凸レンズを使う。

入力は第一面の屈折力 $\Phi_1=10\,\mathrm{m^{-1}}$、第二面の屈折力 $\Phi_2=0\,\mathrm{m^{-1}}$、換算厚さ $e=\frac{1}{100}\,\mathrm m$、物面と入射頂点の距離 $a=\frac{1}{5}\,\mathrm m$ である。

第一面での角変化を内部伝播の高さへ換算する。

$$
\begin{aligned}
A&=- e \Phi_1 + 1,\\
A&= - \left(\frac{1}{100}\right) \cdot \left(10\right) + 1=\frac{9}{10}.
\end{aligned}
$$

高さゼロから入れると第二面への高さは $e p_{\rm ray}$ である。

$$
\begin{aligned}
D&=- e \Phi_2 + 1,\\
D&= - \left(\frac{1}{100}\right) \cdot \left(0\right) + 1=1.
\end{aligned}
$$

二面の屈折力に距離の補正を加える。

$$
\begin{aligned}
\Phi&=- e \Phi_1 \Phi_2 + \Phi_1 + \Phi_2,\\
\Phi/(\mathrm{m^{-1}})&= - \left(\frac{1}{100}\right) \cdot \left(10\right) \cdot \left(0\right) + \left(10\right) + \left(0\right)=10.
\end{aligned}
$$

入射が平行な場合に出力から軸へ達する距離を取る。

$$
\begin{aligned}
b_{\infty}&=\frac{A}{\Phi},\\
b_{\infty}/(\mathrm m)&= \frac{\left(\frac{9}{10}\right)}{\left(10\right)}=\frac{9}{100}.
\end{aligned}
$$

出射角に物面の出発角の情報が残っても、像面の高さでは消えるように伝播距離を選ぶ。

$$
\begin{aligned}
b&=\frac{- A a - e}{D - \Phi a},\\
b/(\mathrm m)&= \frac{- \left(\frac{9}{10}\right) \cdot \left(\frac{1}{5}\right) - \left(\frac{1}{100}\right)}{\left(1\right) - \left(10\right) \cdot \left(\frac{1}{5}\right)}=\frac{19}{100}.
\end{aligned}
$$

行列式が一なので同じ分母の逆数が高さの倍率になる。

$$
\begin{aligned}
M&=\frac{1}{D - \Phi a},\\
M&= \frac{1}{\left(1\right) - \left(10\right) \cdot \left(\frac{1}{5}\right)}=-1.
\end{aligned}
$$

平面では光軸上の換算角は変えないが、内部伝播が高さを変えるため頂点からの距離に差が出る。
:::

:::example[前面が平面のレンズ]
同じレンズを反転し、物体距離を同じ頂点基準で測る。

入力は第一面の屈折力 $\Phi_1=0\,\mathrm{m^{-1}}$、第二面の屈折力 $\Phi_2=10\,\mathrm{m^{-1}}$、換算厚さ $e=\frac{1}{100}\,\mathrm m$、物面と入射頂点の距離 $a=\frac{1}{5}\,\mathrm m$ である。

第一面での角変化を内部伝播の高さへ換算する。

$$
\begin{aligned}
A&=- e \Phi_1 + 1,\\
A&= - \left(\frac{1}{100}\right) \cdot \left(0\right) + 1=1.
\end{aligned}
$$

高さゼロから入れると第二面への高さは $e p_{\rm ray}$ である。

$$
\begin{aligned}
D&=- e \Phi_2 + 1,\\
D&= - \left(\frac{1}{100}\right) \cdot \left(10\right) + 1=\frac{9}{10}.
\end{aligned}
$$

二面の屈折力に距離の補正を加える。

$$
\begin{aligned}
\Phi&=- e \Phi_1 \Phi_2 + \Phi_1 + \Phi_2,\\
\Phi/(\mathrm{m^{-1}})&= - \left(\frac{1}{100}\right) \cdot \left(0\right) \cdot \left(10\right) + \left(0\right) + \left(10\right)=10.
\end{aligned}
$$

入射が平行な場合に出力から軸へ達する距離を取る。

$$
\begin{aligned}
b_{\infty}&=\frac{A}{\Phi},\\
b_{\infty}/(\mathrm m)&= \frac{\left(1\right)}{\left(10\right)}=\frac{1}{10}.
\end{aligned}
$$

出射角に物面の出発角の情報が残っても、像面の高さでは消えるように伝播距離を選ぶ。

$$
\begin{aligned}
b&=\frac{- A a - e}{D - \Phi a},\\
b/(\mathrm m)&= \frac{- \left(1\right) \cdot \left(\frac{1}{5}\right) - \left(\frac{1}{100}\right)}{\left(\frac{9}{10}\right) - \left(10\right) \cdot \left(\frac{1}{5}\right)}=\frac{21}{110}.
\end{aligned}
$$

行列式が一なので同じ分母の逆数が高さの倍率になる。

$$
\begin{aligned}
M&=\frac{1}{D - \Phi a},\\
M&= \frac{1}{\left(\frac{9}{10}\right) - \left(10\right) \cdot \left(\frac{1}{5}\right)}=- \frac{10}{11}.
\end{aligned}
$$

合成屈折力は変わらないが、頂点基準の後焦点も倍率も前例と異なる。
:::

## 絞りと瞳・視野

### 開口絞りと入射瞳・射出瞳
軸上物点から出た近軸光線のうち、通過できる角を最も強く制限する部品を :term[開口絞り]{en="aperture stop" yomi="かいこうしぼり"} という。各部品までの転送行列の右上を $B_j$、半径を $a_j$ とすれば、その面の高さは軸上物点に対して $B_j p$ である。制約は $|p|\le a_j/|B_j|$ で、$B_j=0$ の面は軸上点を制限しない。開口絞りを物体側の前群で逆向きに結像した像が :term[入射瞳]{en="entrance pupil" yomi="にゅうしゃひとみ"}、後群で像側に結像した像が :term[射出瞳]{en="exit pupil" yomi="しゃしゅつひとみ"} である。瞳は虚像や無限遠の像にもなり得る。物体側から見た瞳の半角を $\alpha$ とすると:term[開口数]{en="numerical aperture" yomi="かいこうすう"} は $\mathrm{NA}=n\sin\alpha$ である。

物点から開口絞りの中心を通る光線を :term[主光線]{en="chief ray" yomi="しゅこうせん"} と呼ぶ。物面や像面と共役な面に置き、見える物体の範囲を画定する絞りを :term[視野絞り]{en="field stop" yomi="しやしぼり"} という。共役面から離れた部品は、まず斜めの光束を部分的に遮り、視野端を徐々に暗くする。主光線が遮られる角と光束全体が遮られる角は一般に異なる。軸外で光束の一部だけが遮られる場合を口径食という。物点の高さ $y_o$ を含む一般の制約は $|A_j y_o+B_j p|\le a_j$ であり、全部品の許容区間の共通部分を取れば軸外の通過光束を求められる。

:::example[絞りの選定]
二面で $B_1=0.10\,\mathrm m,\ B_2=0.20\,\mathrm m$、半径を $4,6\,\mathrm{mm}$ とする。換算角の上限は $0.04,0.03$ であり、半径が大きい第 2 面が開口絞りになる。
:::
:::example[口径食の区間]
軸外点からの光束に対して二面の角制約が $[-0.04,0.04]$ と $[0,0.06]$ なら、通るのは $[0,0.04]$ である。主光線が通過していても光束全体は欠け得る。
:::

### 射出瞳が無限遠にある系
絞りを後群の前焦点面へ置くと、絞り中心からの主光線は後群で光軸に平行になる。これは :term[像側テレセントリック]{en="image-space telecentricity" yomi="ぞうがわてれせんとりっく"} の条件である。像面が軸方向に少しずれても各主光線の高さが変わらないため、ぼけ中心の倍率は一次では変わらない。物体側でも同様に入射瞳を無限遠へ置ける。物体側主光線が平行なら、厚さのある部品を測る際に、奥行による透視の倍率差を抑えられる。光軸に対称な開口を仮定しており、口径食で点像の重心が変わる場合には別に評価する。

:::example[絞りの像の倍率]
絞りを焦点距離 $20\,\mathrm{mm}$ の後群の前方 $30\,\mathrm{mm}$ に置く。薄レンズでは絞りの共役像は後方 $60\,\mathrm{mm}$、横倍率は $-2$。半径 $1\,\mathrm{mm}$ の絞りの射出瞳は半径 $2\,\mathrm{mm}$ になる。
:::

### 明るさと透過の上限
近軸の円形光束で、物体面の微小面積を $\dd A$、受ける立体角を $\dd\Omega$ とする。二横方向について換算角の小面積は $n^2\dd\Omega$ である。各方向で行列式 1 の光線行列を用いると、四次元の光線の体積要素 $n^2\dd A\dd\Omega$ が保存する。これは近軸でのエタンデュの保存である。同じ媒質で開口絞りの大きさを固定したなら、面積を狭めて集光するには角度の広がりが要る。絞りを閉じると、透過する光線だけで数えたエタンデュは減る。

### 放射輝度と画面周辺の照度
放射輝度を $\mathcal L=\dd P/(\cos\theta\,\dd A\,\dd\Omega)$ と定義する。ここで $\theta$ は光線と面法線の角である。無損失で近軸の光学系では、細い光線管のパワー $\dd P$ と $n^2\dd A\,\dd\Omega$ がそれぞれ保存する。従って $\mathcal L/n^2$ が光線管に沿って一定である。

像側で、射出瞳を面積 $A_{\rm p}$ の小さな平面とし、像面とは平行で軸方向に距離 $b$ 離れているとする。像点から瞳中心への角を $\theta$ とすると、斜距離は $b/\cos\theta$ で、見込む立体角は $A_{\rm p}\cos^3\theta/b^2$ になる。像面の投影からもう一つ $\cos\theta$ が掛かり、放射輝度と透過率が一定なら像面の放射照度は $\cos^4\theta$ に比例する。瞳の口径食や形の変化はここに含めていない。

:::example[口径食がなくても生じる周辺減光]
この平面瞳のモデルで像側角が $30^\circ$ の点は、同じ放射輝度の軸上点に対して $\cos^4 30^\circ=9/16$ の照度になる。像側テレセントリックな系では主光線の $\theta$ がゼロになるので、この配置に由来する角度因子は生じない。
:::

### 一次元の瞳で口径食を計算する

物点からの光線を $p=n\theta$ で数える。二面に対して $y_j=A_j y_o+B_jp$ とし、ここでは $B_j>0$、半口径を $a_j>0$ とする。第 $j$ 面の角範囲は $[-a_j/B_j-A_jy_o/B_j,\ a_j/B_j-A_jy_o/B_j]$ である。二面とも通過するには下限の最大と上限の最小の間に入る必要があり、幅は差の正の部分になる。角当たりの入射パワーを一定とした一次元モデルなら、軸上での幅で割ると相対透過量が得られる。二次元の円瞳では範囲の面積を使う。縦横が分離した矩形絞りなら二方向の幅を掛けてよい。幅が消える前にも、共通区間の中点が移って像の非対称なけられを生じる。

:::example[軸上の開口絞りを選ぶ]
半径が異なる二つの開口を通す。

入力は第一面の半口径 $a_{1}=\frac{1}{250}\,\mathrm m$、第二面の半口径 $a_{2}=\frac{3}{500}\,\mathrm m$、第一面の角から高さへの係数 $B_{1}=\frac{1}{10}\,\mathrm m$、第二面の角から高さへの係数 $B_{2}=\frac{1}{5}\,\mathrm m$、物点の高さ $y_o=0\,\mathrm m$ である。

以下では第一面の像高係数を $1$ とする。

$$
\begin{aligned}
u_1&=\frac{a_{1}}{B_{1}},\\
u_1&= \frac{\left(\frac{1}{250}\right)}{\left(\frac{1}{10}\right)}=\frac{1}{25}.
\end{aligned}
$$

第二面の像高係数も $1$ とする。

$$
\begin{aligned}
u_2&=\frac{a_{2}}{B_{2}},\\
u_2&= \frac{\left(\frac{3}{500}\right)}{\left(\frac{1}{5}\right)}=\frac{3}{100}.
\end{aligned}
$$

すべての下限を満たす最小の換算角を取る。

$$
\begin{aligned}
p_-&=\max\left(- u_1 - \frac{y_o}{B_{1}}, - u_2 - \frac{y_o}{B_{2}}\right),\\
p_-&= \max\left(- \left(\frac{1}{25}\right) - \frac{\left(0\right)}{\left(\frac{1}{10}\right)}, - \left(\frac{3}{100}\right) - \frac{\left(0\right)}{\left(\frac{1}{5}\right)}\right)=- \frac{3}{100}.
\end{aligned}
$$

一つでも上限を越える光線は除く。

$$
\begin{aligned}
p_+&=\min\left(u_1 - \frac{y_o}{B_{1}}, u_2 - \frac{y_o}{B_{2}}\right),\\
p_+&= \min\left(\left(\frac{1}{25}\right) - \frac{\left(0\right)}{\left(\frac{1}{10}\right)}, \left(\frac{3}{100}\right) - \frac{\left(0\right)}{\left(\frac{1}{5}\right)}\right)=\frac{3}{100}.
\end{aligned}
$$

角の共通区間がなければ幅はゼロである。

$$
\begin{aligned}
\Delta p&=\max\left(0, - p_- + p_+\right),\\
\Delta p&= \max\left(0, - \left(- \frac{3}{100}\right) + \left(\frac{3}{100}\right)\right)=\frac{3}{50}.
\end{aligned}
$$

一定の角密度の仮定で、幅を軸上の幅で規格化する。

$$
\begin{aligned}
g&=\frac{\Delta p}{2 \min\left(u_1, u_2\right)},\\
g&= \frac{\left(\frac{3}{50}\right)}{2 \cdot \min\left(\left(\frac{1}{25}\right), \left(\frac{3}{100}\right)\right)}=1.
\end{aligned}
$$

角半幅が小さいのは第二面である。半径だけで開口絞りを選べない。
:::

:::example[視野を外へ広げる]
前例と同じ絞りで物点を光軸から離す。

入力は第一面の半口径 $a_{1}=\frac{1}{250}\,\mathrm m$、第二面の半口径 $a_{2}=\frac{3}{500}\,\mathrm m$、第一面の角から高さへの係数 $B_{1}=\frac{1}{10}\,\mathrm m$、第二面の角から高さへの係数 $B_{2}=\frac{1}{5}\,\mathrm m$、物点の高さ $y_o=\frac{1}{125}\,\mathrm m$ である。

以下では第一面の像高係数を $1$ とする。

$$
\begin{aligned}
u_1&=\frac{a_{1}}{B_{1}},\\
u_1&= \frac{\left(\frac{1}{250}\right)}{\left(\frac{1}{10}\right)}=\frac{1}{25}.
\end{aligned}
$$

第二面の像高係数も $1$ とする。

$$
\begin{aligned}
u_2&=\frac{a_{2}}{B_{2}},\\
u_2&= \frac{\left(\frac{3}{500}\right)}{\left(\frac{1}{5}\right)}=\frac{3}{100}.
\end{aligned}
$$

すべての下限を満たす最小の換算角を取る。

$$
\begin{aligned}
p_-&=\max\left(- u_1 - \frac{y_o}{B_{1}}, - u_2 - \frac{y_o}{B_{2}}\right),\\
p_-&= \max\left(- \left(\frac{1}{25}\right) - \frac{\left(\frac{1}{125}\right)}{\left(\frac{1}{10}\right)}, - \left(\frac{3}{100}\right) - \frac{\left(\frac{1}{125}\right)}{\left(\frac{1}{5}\right)}\right)=- \frac{7}{100}.
\end{aligned}
$$

一つでも上限を越える光線は除く。

$$
\begin{aligned}
p_+&=\min\left(u_1 - \frac{y_o}{B_{1}}, u_2 - \frac{y_o}{B_{2}}\right),\\
p_+&= \min\left(\left(\frac{1}{25}\right) - \frac{\left(\frac{1}{125}\right)}{\left(\frac{1}{10}\right)}, \left(\frac{3}{100}\right) - \frac{\left(\frac{1}{125}\right)}{\left(\frac{1}{5}\right)}\right)=- \frac{1}{25}.
\end{aligned}
$$

角の共通区間がなければ幅はゼロである。

$$
\begin{aligned}
\Delta p&=\max\left(0, - p_- + p_+\right),\\
\Delta p&= \max\left(0, - \left(- \frac{7}{100}\right) + \left(- \frac{1}{25}\right)\right)=\frac{3}{100}.
\end{aligned}
$$

一定の角密度の仮定で、幅を軸上の幅で規格化する。

$$
\begin{aligned}
g&=\frac{\Delta p}{2 \min\left(u_1, u_2\right)},\\
g&= \frac{\left(\frac{3}{100}\right)}{2 \cdot \min\left(\left(\frac{1}{25}\right), \left(\frac{3}{100}\right)\right)}=\frac{1}{2}.
\end{aligned}
$$

どちらかの面だけを通る角が増え、共通範囲は狭くなる。
:::

:::example[両面を通れない物点]
さらに物点を外へ動かす。

入力は第一面の半口径 $a_{1}=\frac{1}{250}\,\mathrm m$、第二面の半口径 $a_{2}=\frac{3}{500}\,\mathrm m$、第一面の角から高さへの係数 $B_{1}=\frac{1}{10}\,\mathrm m$、第二面の角から高さへの係数 $B_{2}=\frac{1}{5}\,\mathrm m$、物点の高さ $y_o=\frac{2}{125}\,\mathrm m$ である。

以下では第一面の像高係数を $1$ とする。

$$
\begin{aligned}
u_1&=\frac{a_{1}}{B_{1}},\\
u_1&= \frac{\left(\frac{1}{250}\right)}{\left(\frac{1}{10}\right)}=\frac{1}{25}.
\end{aligned}
$$

第二面の像高係数も $1$ とする。

$$
\begin{aligned}
u_2&=\frac{a_{2}}{B_{2}},\\
u_2&= \frac{\left(\frac{3}{500}\right)}{\left(\frac{1}{5}\right)}=\frac{3}{100}.
\end{aligned}
$$

すべての下限を満たす最小の換算角を取る。

$$
\begin{aligned}
p_-&=\max\left(- u_1 - \frac{y_o}{B_{1}}, - u_2 - \frac{y_o}{B_{2}}\right),\\
p_-&= \max\left(- \left(\frac{1}{25}\right) - \frac{\left(\frac{2}{125}\right)}{\left(\frac{1}{10}\right)}, - \left(\frac{3}{100}\right) - \frac{\left(\frac{2}{125}\right)}{\left(\frac{1}{5}\right)}\right)=- \frac{11}{100}.
\end{aligned}
$$

一つでも上限を越える光線は除く。

$$
\begin{aligned}
p_+&=\min\left(u_1 - \frac{y_o}{B_{1}}, u_2 - \frac{y_o}{B_{2}}\right),\\
p_+&= \min\left(\left(\frac{1}{25}\right) - \frac{\left(\frac{2}{125}\right)}{\left(\frac{1}{10}\right)}, \left(\frac{3}{100}\right) - \frac{\left(\frac{2}{125}\right)}{\left(\frac{1}{5}\right)}\right)=- \frac{3}{25}.
\end{aligned}
$$

角の共通区間がなければ幅はゼロである。

$$
\begin{aligned}
\Delta p&=\max\left(0, - p_- + p_+\right),\\
\Delta p&= \max\left(0, - \left(- \frac{11}{100}\right) + \left(- \frac{3}{25}\right)\right)=0.
\end{aligned}
$$

一定の角密度の仮定で、幅を軸上の幅で規格化する。

$$
\begin{aligned}
g&=\frac{\Delta p}{2 \min\left(u_1, u_2\right)},\\
g&= \frac{\left(0\right)}{2 \cdot \min\left(\left(\frac{1}{25}\right), \left(\frac{3}{100}\right)\right)}=0.
\end{aligned}
$$

個々の面を通る方向は存在するが、二面を順に通る方向はない。
:::

## 光学機器（眼・カメラ・顕微鏡・望遠鏡）

### 眼とカメラ
眼の網膜までの距離はほぼ固定であり、光学系の屈折力を変えて合焦する。像距離を $b$ とすれば必要屈折力は単純化した空気中モデルで $1/a+1/b$ である。遠方から距離 $a$ へ合わせるには $1/a$ を加える。カメラでは薄レンズの位置を動かしてセンサーと像面を一致させる。無限遠で焦点距離 $f$、入射瞳直径 $D$ のとき :term[F 値]{en="f-number" yomi="えふち"} は $N=f/D$。円形瞳の角面積は近軸で $\pi(D/2f)^2$ なので、一様被写体に対する軸上照度は透過率一定で $1/N^2$ に比例する。

### 顕微鏡と望遠鏡
明視距離を $d_0=0.25\,\mathrm m$ とする。接眼レンズの前焦点面に高さ $h$ の中間像を置くと平行出射の角は大きさ $h/f_e$ である。裸眼で $d_0$ から見る角 $h/d_0$ との比は $d_0/f_e$。顕微鏡ではこれに対物レンズの横倍率 $|b/a|$ を掛ける。Kepler 型望遠鏡では対物焦点面と接眼前焦点面を共有させる。遠方物体の角 $\alpha$ に対し中間像高さは $f_o\alpha$、出射角は $-f_o\alpha/f_e$ である。入射瞳の径 $D$ に対する射出瞳径は $D f_e/f_o$ になる。これは同じ媒質中で $C=0$、$AD=1$ のためである。

:::example[顕微鏡の倍率]
対物倍率が $-20$、接眼焦点距離が $25\,\mathrm{mm}$ なら、角倍率の大きさは $20(250/25)=200$ である。倍率を増しても開口数が同じなら回折で失われた細部は戻らない。
:::
:::example[カメラの絞り]
$N=2$ から $N=4$ では軸上の照度が $1/4$ になる。同じ露光量を得るには露光時間を 4 倍にする。
:::

### 解像と倍率を別に評価する
円形の無収差瞳で独立な二物点を写す場合、円孔回折の第一零点を尺度にした物体側の距離は $0.61\lambda_0/\mathrm{NA}$ である。Airy パターンの第一零点は [回折理論](/physics/optics/diffraction/) で導く。近軸の Lagrange 不変量から $\mathrm{NA}_{\mathrm{object}}\simeq |M|\mathrm{NA}_{\mathrm{image}}$ なので、像側の第一零点の半径を $|M|$ で割っても同じ物体側距離を得る。これは近軸で、かつ光路が回転対称で大きな収差がない場合の基準である。倍率を増やす役割は既に伝達された形を大きくすることにある。

:::example[解像の目安]
$\lambda_0=550\,\mathrm{nm}$、対物の開口数 $0.55$ なら Rayleigh の物体側距離は約 $0.61\,\mathrm{\mu m}$ である。接眼倍率だけを倍にしても、この尺度は変わらない。
:::

### 虚像を作る接眼と鏡筒長
Kepler 型は対物と接眼がともに正で、焦点を共有するためレンズ間隔は $f_o+f_e$。接眼を負レンズ $f_e<0$ にすると、対物が作るはずの実像の手前に接眼を置き、接眼から見て虚物とする。間隔は同じ式 $f_o+f_e$ だが負の分だけ短い。式 :eqref[eq-optics-ray-matrix] の二レンズ積でも $C=0$ で、角倍率は $-f_o/f_e>0$。これが Galilei 型で、正立する。ただし軸外の光束が有限の接眼で遮られやすく、視野と鏡筒長の双方を決める必要がある。

:::example[間隔と倍率の比較]
$f_o=200\,\mathrm{mm}$、$f_e=\pm40\,\mathrm{mm}$ なら鏡筒長は $240\,\mathrm{mm}$ と $160\,\mathrm{mm}$、角倍率は $-5$ と $+5$ である。同じ倍率の大きさでも像の向きと機械長は異なる。
:::

### 虚像を近点に置く拡大鏡
正の薄レンズを目のすぐ前に置き、物体を前焦点の内側に置く。明視距離 $d_0$ に虚像を作るなら $b=-d_0$ である。薄レンズ式を解くと $1/a=1/f+1/d_0$ で、物体の高さを $h$ とすれば虚像の高さは $d_0h/a$ となる。目に入る角は近軸で $h/a$ なので、裸眼で明視距離から見る場合の角 $h/d_0$ に対して角倍率は $1+d_0/f$ である。目を遠方に合わせる場合は $a=f$ とし、角倍率は $d_0/f$ になる。

:::example[目の調節による拡大率の差]
$f=50\,\mathrm{mm}$、$d_0=250\,\mathrm{mm}$ の拡大鏡では、近点に虚像を置く場合は $a=125/3\,\mathrm{mm}$ で 6 倍、平行に出す場合は $a=50\,\mathrm{mm}$ で 5 倍になる。
:::

### 近接撮影で実効 F 値が増える理由

薄レンズの位置に直径 $D$ の絞りを置き、物面と像面を平行とする。像面で瞳を見込む半角は近軸で $D/(2b)$、従って一様な放射輝度の軸上像で照度はこの二乗に比例する。無限遠での像距離 $f$ を有限の $b$ で置き換えると、実効 F 値は $N_{\rm eff}=b/D$ となる。薄レンズ式から $b=f(1+b/a)$ であるため $N_{\rm eff}=N(1+|M|)$。焦点を近くに調整したときの露光量を保つなら、露光時間を $(1+|M|)^2$ 倍にする。これは瞳倍率一の近軸モデルであり、実際の光学系で射出瞳が異なる場合の半角は射出瞳から測る。

:::example[遠くの平面を撮る]
薄レンズで十分遠くにある一様な物面を撮る。

入力は焦点距離 $f=\frac{1}{20}\,\mathrm m$、レンズ面の開口径 $D=\frac{1}{50}\,\mathrm m$、物体距離 $a=10\,\mathrm m$、真空波長 $\lambda_0\simeq 5.5\times 10^{-7}\,\mathrm m$ である。

撮影距離から像が合う位置を求める。

$$
\begin{aligned}
b&=\frac{a f}{a - f},\\
b/(\mathrm m)&= \frac{\left(10\right) \cdot \left(\frac{1}{20}\right)}{\left(10\right) - \left(\frac{1}{20}\right)}=\frac{10}{199}.
\end{aligned}
$$

センサー上の縮小または拡大の程度を求める。

$$
\begin{aligned}
M_a&=\frac{b}{a},\\
M_a&= \frac{\left(\frac{10}{199}\right)}{\left(10\right)}=\frac{1}{199}.
\end{aligned}
$$

無限遠へ合わせた場合の開口を基準に取る。

$$
\begin{aligned}
N&=\frac{f}{D},\\
N&= \frac{\left(\frac{1}{20}\right)}{\left(\frac{1}{50}\right)}=\frac{5}{2}.
\end{aligned}
$$

有限の像距離で見込む円錐角に換算する。

$$
\begin{aligned}
N_{\rm eff}&=\frac{b}{D},\\
N_{\rm eff}&= \frac{\left(\frac{10}{199}\right)}{\left(\frac{1}{50}\right)}=\frac{500}{199}.
\end{aligned}
$$

像側の立体角の比の逆数を取る。

$$
\begin{aligned}
g&=\frac{N_{\rm eff}^{2}}{N^{2}},\\
g&= \frac{\left(\frac{500}{199}\right)^{2}}{\left(\frac{5}{2}\right)^{2}}\simeq 1.0101.
\end{aligned}
$$

一様な円形瞳を仮定し、同じ角度の縮小が回折に与える効果も求める。

$$
\begin{aligned}
r_{\rm A}&=\frac{61 N_{\rm eff} \lambda_0}{50},\\
r_{\rm A}/(\mathrm m)&\simeq \frac{61 \cdot \left(\frac{500}{199}\right) \cdot \left(5.5\times 10^{-7}\right)}{50}\simeq 1.6859\times 10^{-6}.
\end{aligned}
$$

倍率が小さく、公称値と実効値は近い。
:::

:::example[等倍撮影]
物体を焦点距離の二倍の位置まで近づける。

入力は焦点距離 $f=\frac{1}{20}\,\mathrm m$、レンズ面の開口径 $D=\frac{1}{50}\,\mathrm m$、物体距離 $a=\frac{1}{10}\,\mathrm m$、真空波長 $\lambda_0\simeq 5.5\times 10^{-7}\,\mathrm m$ である。

撮影距離から像が合う位置を求める。

$$
\begin{aligned}
b&=\frac{a f}{a - f},\\
b/(\mathrm m)&= \frac{\left(\frac{1}{10}\right) \cdot \left(\frac{1}{20}\right)}{\left(\frac{1}{10}\right) - \left(\frac{1}{20}\right)}=\frac{1}{10}.
\end{aligned}
$$

センサー上の縮小または拡大の程度を求める。

$$
\begin{aligned}
M_a&=\frac{b}{a},\\
M_a&= \frac{\left(\frac{1}{10}\right)}{\left(\frac{1}{10}\right)}=1.
\end{aligned}
$$

無限遠へ合わせた場合の開口を基準に取る。

$$
\begin{aligned}
N&=\frac{f}{D},\\
N&= \frac{\left(\frac{1}{20}\right)}{\left(\frac{1}{50}\right)}=\frac{5}{2}.
\end{aligned}
$$

有限の像距離で見込む円錐角に換算する。

$$
\begin{aligned}
N_{\rm eff}&=\frac{b}{D},\\
N_{\rm eff}&= \frac{\left(\frac{1}{10}\right)}{\left(\frac{1}{50}\right)}=5.
\end{aligned}
$$

像側の立体角の比の逆数を取る。

$$
\begin{aligned}
g&=\frac{N_{\rm eff}^{2}}{N^{2}},\\
g&= \frac{\left(5\right)^{2}}{\left(\frac{5}{2}\right)^{2}}=4.
\end{aligned}
$$

一様な円形瞳を仮定し、同じ角度の縮小が回折に与える効果も求める。

$$
\begin{aligned}
r_{\rm A}&=\frac{61 N_{\rm eff} \lambda_0}{50},\\
r_{\rm A}/(\mathrm m)&\simeq \frac{61 \cdot \left(5\right) \cdot \left(5.5\times 10^{-7}\right)}{50}\simeq 3.355\times 10^{-6}.
\end{aligned}
$$

実効 F 値は二倍、露光時間は四倍になる。像面の回折半径も大きくなる。
:::

:::example[近接で絞りを閉じる]
等倍に保ったまま開口径を小さくした場合を比較する。

入力は焦点距離 $f=\frac{1}{20}\,\mathrm m$、レンズ面の開口径 $D=\frac{1}{100}\,\mathrm m$、物体距離 $a=\frac{1}{10}\,\mathrm m$、真空波長 $\lambda_0\simeq 5.5\times 10^{-7}\,\mathrm m$ である。

撮影距離から像が合う位置を求める。

$$
\begin{aligned}
b&=\frac{a f}{a - f},\\
b/(\mathrm m)&= \frac{\left(\frac{1}{10}\right) \cdot \left(\frac{1}{20}\right)}{\left(\frac{1}{10}\right) - \left(\frac{1}{20}\right)}=\frac{1}{10}.
\end{aligned}
$$

センサー上の縮小または拡大の程度を求める。

$$
\begin{aligned}
M_a&=\frac{b}{a},\\
M_a&= \frac{\left(\frac{1}{10}\right)}{\left(\frac{1}{10}\right)}=1.
\end{aligned}
$$

無限遠へ合わせた場合の開口を基準に取る。

$$
\begin{aligned}
N&=\frac{f}{D},\\
N&= \frac{\left(\frac{1}{20}\right)}{\left(\frac{1}{100}\right)}=5.
\end{aligned}
$$

有限の像距離で見込む円錐角に換算する。

$$
\begin{aligned}
N_{\rm eff}&=\frac{b}{D},\\
N_{\rm eff}&= \frac{\left(\frac{1}{10}\right)}{\left(\frac{1}{100}\right)}=10.
\end{aligned}
$$

像側の立体角の比の逆数を取る。

$$
\begin{aligned}
g&=\frac{N_{\rm eff}^{2}}{N^{2}},\\
g&= \frac{\left(10\right)^{2}}{\left(5\right)^{2}}=4.
\end{aligned}
$$

一様な円形瞳を仮定し、同じ角度の縮小が回折に与える効果も求める。

$$
\begin{aligned}
r_{\rm A}&=\frac{61 N_{\rm eff} \lambda_0}{50},\\
r_{\rm A}/(\mathrm m)&\simeq \frac{61 \cdot \left(10\right) \cdot \left(5.5\times 10^{-7}\right)}{50}\simeq 6.71\times 10^{-6}.
\end{aligned}
$$

倍率は同じであるが、絞りによる追加の照度低下と回折半径の増加を生じる。
:::

## 単色収差と Seidel の 5 収差

### 回転対称系での四次展開
理想的な参照球面からの光路長の差を :term[波面収差]{en="wave aberration" yomi="はめんしゅうさ"} $W$ とする。瞳座標を二次元ベクトル $\vb u$、視野を $\vb h$ で表し、回転と鏡映に対して対称な光学系を考える。反射対称性を含め、二つのベクトルの相対位置は $u^2,h^2,\vb u\cdot\vb h$ で決まる。四次の不変多項式はこれら二次式の二次積で張られる。実際、二次元の回転を複素座標で表すと、不変な四次単項式は複素座標とその共役を二つずつ含む。これらを内積と向き付き面積の積で書けるが、鏡映は向き付き面積の符号を反転させ、その平方は $u^2 h^2-(\vb u\cdot\vb h)^2$ に等しい。瞳によらない $h^4$ を除くと
$$
W=A u^4+B u^2(\vb u\cdot\vb h)+C(\vb u\cdot\vb h)^2
+D h^2u^2+E h^2(\vb u\cdot\vb h).
\label{eq-optics-seidel}
$$
順に球面収差、コマ収差、非点収差、瞳方向に等方な像面湾曲、歪曲の項である。ただし、この書き方で等方な項の係数 $D$ には非点収差に結び付いた寄与も含まれる。Petzval 湾曲だけを表す組合せは後で導く。瞳方向の位相勾配が光線角を変えるので、横光線誤差は幾何学的な比例係数を除いて $\grad_{\vb u}W$ に従う。例えば $\vb h$ を $x$ 方向に取ると、$C$ 項は $u_x^2$ だけに作用し、$D$ 項は両方向の焦点を同時に動かす。$E$ 項の勾配は瞳位置によらないため、同一物点からの像を横にずらし、格子の形を曲げる。

:::example[軸上で開口を絞る]
$W=A u^4$ なら誤差角は $4A u^3$ に比例する。瞳半径を半分にすると縁での光線誤差は $1/8$ になる。ただし絞ると回折像は広がる。
:::
:::example[視野とコマ]
$B u^2hu_x$ の勾配は $Bh(3u_x^2+u_y^2,2u_xu_y)$ である。$h$ の符号を反転すると、像の非対称な伸びの向きも反転する。
:::

### コマを抑える Abbe の正弦条件
回転対称系で、軸上の物点と像点が有限開口の全光線について同じ光路長で結ばれるとする。次に物点と像点を子午面内でそれぞれ横に $h$ と $Mh$ ずらす。固定した瞳点で前後の経路を停留させると、光路長の一次変化は :eqref[eq-optics-optical-variation] の端点項だけになり、$h(n_2 M\sin u_2-n_1\sin u_1)$ である。$u_1,u_2$ は軸上光線の物側・像側の符号付き角であり、同じ軸で測る。

この一次変化が光線によらなければ、軸上結像とともに小視野での結像も一次まで保たれる。軸上の中心光線では両角がゼロなので、必要十分な一次条件は全開口で $n_1\sin u_1=n_2 M\sin u_2$ となることである。これを :term[Abbe の正弦条件]{en="Abbe sine condition" yomi="あっべのせいげんじょうけん"} という。近軸で $\sin u\simeq u$ とすると Lagrange の光学不変量と一致する。有限開口で成立すれば、視野 $h$ に一次の波面収差が消え、Seidel のコマ項も消える。

:::example[有限開口で角を換算する]
両側が空気で $M=-2$、軸上のある光線の物側角が $30^\circ$ であるとする。正弦条件が成り立つ系では、像側角は $\arcsin(-1/4)\simeq-14.48^\circ$ となる。
:::

### 焦点の曲がりと非点収差を分ける
視野を $x$ 方向に取ると、:eqref[eq-optics-seidel] のうち瞳に二次の項は $h^2[(C+D)u_x^2+D u_y^2]$ となる。像側が屈折率 $n_2$ の一様媒質で、物理座標 $\vb u$ で表した射出瞳から参照像面までの距離を $\ell$ とする。参照像点を光の進行側へ $\Delta z$ 移すと、瞳から軸近くの像点までの距離は
$$
\sqrt{(\ell+\Delta z)^2+u^2}-\sqrt{\ell^2+u^2}
=\Delta z-\frac{\Delta z\,u^2}{2\ell^2}+\cdots
$$
と変わる。ここで視野二次の焦点移動を求めるため、$\ell$ と投影係数には軸上値を使う。ピストンを除き光路長換算すると両方向に $-n_2\Delta z\,u^2/(2\ell^2)$ が加わる。従って $x$ 方向の焦点ずれは $2\ell^2h^2(C+D)/n_2$、$y$ 方向は $2\ell^2h^2D/n_2$。二焦点の差は $C$ で決まり、$C=0$ でも $D$ が残れば最良像面は曲がる。各焦点ずれから、非点収差由来の湾曲と各面の屈折力で決まる湾曲を分離できる。その組合せを多面系の光路長から求める。

### 軸上の収差を実際の曲面から計算する
曲率半径 $\rho$ の凹面鏡の頂点から入射側へ座標 $s$ を取り、$s(r)=\rho-\sqrt{\rho^2-r^2}$ とする。入射平面 $s=S$ から鏡を経て $s=f$ の軸上点へ向かう幾何学的な長さは $S-s(r)+\sqrt{(f-s(r))^2+r^2}$。軸上の長さ $S+f$ を引き、$s(r)=r^2/(2\rho)+r^4/(8\rho^3)+O(r^6)$ を代入する。距離の平方根は $f-s+r^2/[2(f-s)]-r^4/[8(f-s)^3]+O(r^6)$ なので
$$
W(r)=\left(\frac1{2f}-\frac1\rho\right)r^2+
\left(\frac{1}{4\rho f^2}-\frac1{8f^3}-\frac1{4\rho^3}\right)r^4+O(r^6).
$$
$f=\rho/2$ では二次が消えるが、四次係数は $-1/(4\rho^3)$ として残る。各瞳点から参照像点までの光路差で定義した符号であり、収差を打ち消す補正板には反対符号を与える。この例では、近軸の計算に四次項を加えると球面収差が現れる。後で述べる光路長の展開では、多面系の波面収差を求める。反射系も含めて個々の面の Seidel 係数を標準形で一覧するには、曲面ごとに座標と符号を定める必要があるため、:cite[born-wolf] の幾何学的収差の扱いを参照する。

:::example[瞳端の波面差]
$\rho=0.25\,\mathrm m$、瞳端 $r=10\,\mathrm{mm}$ なら四次近似で $W(r)=-0.16\,\mathrm{\mu m}$。真空波長 $500\,\mathrm{nm}$ の約 $0.32$ 波に相当するため、小さい角度でも位相差を無視できない。
:::

### 斜入射の球面が作る二焦線
一枚の曲面の主光線が通る位置を原点とし、そこでの面を $z=(x^2+y^2)/(2R)+O((x^2+y^2)^2)$ と置く。入射面を $xz$ とし、物点は $-a(\sin i,0,\cos i)$、像候補は $b(\sin t,0,\cos t)$ にあるとする。$a,b$ は主光線沿いの正距離、$i,t$ は Snell 則を満たす。曲面点 $\vb r=(x,y,z)$ までの入射光路の展開は
$$
n_1|\vb r+a\hat{\vb s}_1|=
n_1\left[a+x\sin i+z\cos i+\frac{x^2\cos^2 i+y^2}{2a}\right]+O(|(x,y)|^3).
$$
像側には、一次の符号を反転して $b,t,n_2$ を入れる。接線の一次項は Snell 則で消える。曲面高さの二次項を含めて、$x^2$ と $y^2$ の係数をそれぞれゼロにする条件は
$$
\begin{aligned}
n_1\cos^2 i/a+n_2\cos^2 t/b_T&=(n_2\cos t-n_1\cos i)/R,\\
n_1/a+n_2/b_S&=(n_2\cos t-n_1\cos i)/R .
\end{aligned}
$$
$b_T$ が入射面内、$b_S$ が直交方向に変化する近傍光線の結像距離である。曲面に対し両方向の二次係数が同時に消えるとは限らず、二つの焦線の間に非点収差が見える。垂直入射なら $\cos i=\cos t=1$ となり、二式が共に球面の近軸公式へ戻る。

:::example[斜めの平行入射]
$n_1=1,n_2=\sqrt2,i=45^\circ$ では $t=30^\circ$。光源を無限遠に置くと $b_T/b_S=\cos^2 t=3/4$ である。球面が集光力を持つ場合、入射面内の近傍光線の焦線がより手前にできる。
:::

### 多面系の四次波面収差を計算する手順
透過だけを持つ回転対称な系を考える。面頂点の距離を $d_j>0$、第 $j$ 区間の屈折率を $n_j$、各面での横座標を $\vb r_j$ とする。終点は近軸像点、途中の一面の座標を指定した瞳点とし、残りの面では光路長を停留させる。面の高さを $s_j=s_j^{(2)}+s_j^{(4)}+O(r_j^6)$ と分ける。球面なら $s_j^{(2)}=r_j^2/(2R_j)$、$s_j^{(4)}=r_j^4/(8R_j^3)$ であり、外部の物面・像面を平面とすれば対応する $s_j$ はゼロである。

$Q_j=|\vb r_{j+1}-\vb r_j|^2$、$\delta s_j=s_{j+1}-s_j$ と置くと、区間長は $\sqrt{(d_j+\delta s_j)^2+Q_j}$。横座標の四次までの光路長は
$$
\begin{aligned}
L^{(2)}&=\sum_j n_j\left[\delta s_j^{(2)}+\frac{Q_j}{2d_j}\right],\\
L^{(4)}&=\sum_j n_j\left[\delta s_j^{(4)}
-\frac{Q_j\delta s_j^{(2)}}{2d_j^2}-\frac{Q_j^2}{8d_j^3}\right]
\end{aligned}
$$
である。これを導くには、距離を $(d_j+\delta s_j)\sqrt{1+Q_j/(d_j+\delta s_j)^2}$ と書き、平方根の $Q_j^2$ の項まで残せばよい。$d_j+\delta s_j>0$ は十分小さい開口で成立する。

まず $L^{(2)}$ を変分して近軸光線の交点を求める。例えば $\partial_{\vb r_j}Q_{j-1}=2(\vb r_j-\vb r_{j-1})$ により、区間の傾きに $n$ を掛けた量の差が面の勾配に比例し、曲面の屈折式 :eqref[eq-optics-surface-power] に一致する。これらを高さの並び $\vb r_j^{(1)}$ とする。次に真の停留経路に対して、物点高さと瞳高さを一律に $\epsilon$ 倍したとき $\vb r_j(\epsilon)=\epsilon \vb r_j^{(1)}+O(\epsilon^3)$ を使う。一次解が一意で変分の Hessian が非退化なら、陰関数定理でこの展開ができ、左右の回転対称性により偶数の項が消える。これを $L^{(2)}$ に代入して生じる四次の補正は、その勾配と三次の座標補正の内積である。しかし指定した瞳・物・像の座標には補正がなく、自由な面の勾配は近軸解でゼロ。従って波面の四次項は $L^{(4)}$ を近軸交点で評価するだけで得られる。瞳によらない値を引けば :eqref[eq-optics-seidel] の係数を比較できる。光軸上で指定点を通る経路に退化がある場合は、停留経路を区間ごとに取り直す必要がある。

### 多面系の Petzval 和を光路長から求める
式 :eqref[eq-optics-seidel] の $D-C/2$ を直接取り出す。各面の近軸交点は $\vb r_j=a_j\vb u+b_j\vb h$ と書ける。$\vb u,\vb h$ のどちらか一方をゼロにして得る二本の子午光線から、係数 $a_j,b_j$ を読めばよい。四次多項式に対して、$u^2h^2$ の係数から $(\vb u\cdot\vb h)^2$ の係数の半分を引く操作を $\mathcal P$ とする。

任意の実数 $a,b$ について $\mathcal P(|a\vb u+b\vb h|^4)=2a^2b^2-(4a^2b^2)/2=0$ である。また $\mathcal P(|a\vb u+b\vb h|^2|c\vb u+d\vb h|^2)=(ad-bc)^2$ となる。従って前節の $L^{(4)}$ のうち、$\delta s_j^{(4)}$ と $Q_j^2$ は寄与しない。$Q_j\delta s_j^{(2)}$ だけが残り、$Q_j$ 中の係数 $(a_{j+1}-a_j,b_{j+1}-b_j)$ と、$s_j^{(2)}$ または $s_{j+1}^{(2)}$ 中の係数の行列式は、二乗するとどちらも $(a_jb_{j+1}-b_ja_{j+1})^2$ である。

二本の光線に対する Lagrange の不変量は $H=(n_j/d_j)(a_jb_{j+1}-b_ja_{j+1})$ である。これは伝播での二本の換算角が $n_j(a_{j+1}-a_j)/d_j$ と $n_j(b_{j+1}-b_j)/d_j$ になることから従う。従って
$$
D-\frac{C}{2}
=-\frac{H^2}{4}\sum_j \frac{1}{n_j}\left(\frac{1}{R_{j+1}}-\frac{1}{R_j}\right)
=-\frac{H^2}{4}\sum_{\rm surfaces}\frac{\Phi_j}{n_{\rm in,j}n_{\rm out,j}} .
\label{eq-optics-petzval}
$$
ここで $\Phi_j=(n_{\rm out,j}-n_{\rm in,j})/R_j$ であり、物面と参照像面は平面とした。二つ目の等号は、同じ曲面の前後の区間の寄与をまとめた結果である。従って面の和は絞りの位置に依存しない。これが :term[Petzval 和]{en="Petzval sum" yomi="ぺっつばーるわ"} である（:cite[born-wolf]）。

像側の屈折率を $n_2$、射出瞳から像面までの距離を $\ell$、横倍率を $M$ とする。瞳座標 $\vb u$ が物理的な射出瞳の高さなら、軸上光線の換算角は $-n_2u/\ell$。像面でのもう一方の光線の高さが $Mh$ なので、$H^2=(n_2 M/\ell)^2$ となる。像高 $Y=Mh$ に対して Petzval 像面は $\Delta z=Y^2/(2R_P)$ と書ける。前節のデフォーカスの式に $D-C/2$ を使うと、$1/R_P=-n_2\sum_j\Phi_j/(n_{\rm in,j}n_{\rm out,j})$ である。この符号では正の空気中の薄レンズで像面の周囲がレンズ側に曲がる。$C\ne0$ のときは、子午・球欠の二像面のいずれも一般には Petzval 像面と一致しない。

:::example[屈折力を残して Petzval 和を打ち消す]
空気中の薄レンズで材質の屈折率を $n$ とすると、二面の Petzval 和は全屈折力を $n$ で割った値になる。式 :eqref[eq-optics-petzval] を接触薄レンズ近似に使う。$n_1=2$、$\Phi_1=8\,\mathrm{m^{-1}}$ のレンズと、屈折率 $1.5$、屈折力 $-6\,\mathrm{m^{-1}}$ のレンズを組み合わせると、面の和は $8/2-6/1.5=0$。一方、合成屈折力は $2\,\mathrm{m^{-1}}$ で、空気中の焦点距離は $0.50\,\mathrm m$ となる。
:::

:::example[平行板による軸上四次項]
光が一様板の厚さ $d$ を一定横波数 $k_0p$ で通過する。板の屈折率を $n$ とすると平面波の位相の光路長換算は $d\sqrt{n^2-p^2}$。同じ厚さの空気との差は $d[(n-1)+(1-1/n)p^2/2+(1-1/n^3)p^4/8+\cdots]$。中心の定数位相を引き、空気側の観測面を軸方向へ $d(1-1/n)$ 進めると、その追加位相の二次項が板の二次項を消す。残る四次の係数は $d(1/n-1/n^3)/8$。$d=10\,\mathrm{mm}$、$n=1.5$ なら係数は約 $0.463\,\mathrm{mm}$ である。従って $p=0.10$ での補正後の波面差は約 $46\,\mathrm{nm}$。焦点位置を動かしても有限の開口で同相には戻らない。
:::

### コマの横光線誤差と瞳での平均

軸外の一点について、参照する像点へ向かう波の位相を引いた差を $k_0W$ とする。像側が空気で、瞳が像から距離 $L$ にある近軸配置では、瞳の物理座標で測った $\grad W$ が角度の誤差を与え、基準面まで進むと横の誤差は $L\grad W$ となる。ここで $W=B x(x^2+y^2)$ のコマだけを残し、視野に比例する係数を $B$ に含める。半径 $a$ の一様円瞳で $\langle x^2\rangle=\langle y^2\rangle=a^2/4$、$\langle xy\rangle=0$ である。これは $x=r\cos\varphi,y=r\sin\varphi$ として円平均を取ると従う。従って光線の $x$ ずれは平均 $LBa^2$、$y$ ずれの平均はゼロである。中心光線の位置を測っても、瞳内の光線を等しく測った重心には一般に一致しない。

端の点 $(a,0)$ でのずれは $3LBa^2$、$(0,a)$ では $LBa^2$。一方、原点ではゼロである。瞳の一部を開いた場合は平均の測度が変わり、重心は開口の対称性にも依存する。

:::example[開口を保った小視野]
視野が小さく、コマ係数も小さい場合を評価する。

入力はコマの係数 $B=\frac{1}{10}\,\mathrm{m^{-2}}$、瞳と像の距離 $L=\frac{1}{10}\,\mathrm m$、瞳半径 $a=\frac{1}{250}\,\mathrm m$、真空波長 $\lambda_0\simeq 5\times 10^{-7}\,\mathrm m$ である。

視野方向の瞳端で固定視野での三次の瞳依存を評価する。

$$
\begin{aligned}
W_{\rm edge}&=B a^{3},\\
W_{\rm edge}/(\mathrm m)&= \left(\frac{1}{10}\right) \cdot \left(\frac{1}{250}\right)^{3}\simeq 6.4\times 10^{-9}.
\end{aligned}
$$

位相を一波長当たりの光路差に換算する。

$$
\begin{aligned}
N_W&=\frac{W_{\rm edge}}{\lambda_0},\\
N_W&\simeq \frac{\left(6.4\times 10^{-9}\right)}{\left(5\times 10^{-7}\right)}=\frac{8}{625}.
\end{aligned}
$$

横端では勾配の三倍の成分が現れる。

$$
\begin{aligned}
\Delta x_{\rm edge}&=3 B L a^{2},\\
\Delta x_{\rm edge}/(\mathrm m)&= 3 \cdot \left(\frac{1}{10}\right) \cdot \left(\frac{1}{10}\right) \cdot \left(\frac{1}{250}\right)^{2}\simeq 4.8\times 10^{-7}.
\end{aligned}
$$

一様円瞳で二次の勾配を面積平均する。

$$
\begin{aligned}
\bar x&=B L a^{2},\\
\bar x/(\mathrm m)&= \left(\frac{1}{10}\right) \cdot \left(\frac{1}{10}\right) \cdot \left(\frac{1}{250}\right)^{2}\simeq 1.6\times 10^{-7}.
\end{aligned}
$$

幾何学的誤差を比較する尺度として無収差の回折半径も求める。

$$
\begin{aligned}
r_{\rm A}&=\frac{61 L \lambda_0}{100 a},\\
r_{\rm A}/(\mathrm m)&\simeq \frac{61 \cdot \left(\frac{1}{10}\right) \cdot \left(5\times 10^{-7}\right)}{100 \cdot \left(\frac{1}{250}\right)}\simeq 7.625\times 10^{-6}.
\end{aligned}
$$

横光線のずれは回折半径より小さい。幾何学的な重心は像の広がりそのものを与えない。
:::

:::example[同じ開口で視野を増す]
係数が視野に比例する範囲で、コマの大きい視野点を選ぶ。

入力はコマの係数 $B=1\,\mathrm{m^{-2}}$、瞳と像の距離 $L=\frac{1}{10}\,\mathrm m$、瞳半径 $a=\frac{1}{250}\,\mathrm m$、真空波長 $\lambda_0\simeq 5\times 10^{-7}\,\mathrm m$ である。

視野方向の瞳端で固定視野での三次の瞳依存を評価する。

$$
\begin{aligned}
W_{\rm edge}&=B a^{3},\\
W_{\rm edge}/(\mathrm m)&= \left(1\right) \cdot \left(\frac{1}{250}\right)^{3}\simeq 6.4\times 10^{-8}.
\end{aligned}
$$

位相を一波長当たりの光路差に換算する。

$$
\begin{aligned}
N_W&=\frac{W_{\rm edge}}{\lambda_0},\\
N_W&\simeq \frac{\left(6.4\times 10^{-8}\right)}{\left(5\times 10^{-7}\right)}=\frac{16}{125}.
\end{aligned}
$$

横端では勾配の三倍の成分が現れる。

$$
\begin{aligned}
\Delta x_{\rm edge}&=3 B L a^{2},\\
\Delta x_{\rm edge}/(\mathrm m)&= 3 \cdot \left(1\right) \cdot \left(\frac{1}{10}\right) \cdot \left(\frac{1}{250}\right)^{2}\simeq 4.8\times 10^{-6}.
\end{aligned}
$$

一様円瞳で二次の勾配を面積平均する。

$$
\begin{aligned}
\bar x&=B L a^{2},\\
\bar x/(\mathrm m)&= \left(1\right) \cdot \left(\frac{1}{10}\right) \cdot \left(\frac{1}{250}\right)^{2}\simeq 1.6\times 10^{-6}.
\end{aligned}
$$

幾何学的誤差を比較する尺度として無収差の回折半径も求める。

$$
\begin{aligned}
r_{\rm A}&=\frac{61 L \lambda_0}{100 a},\\
r_{\rm A}/(\mathrm m)&\simeq \frac{61 \cdot \left(\frac{1}{10}\right) \cdot \left(5\times 10^{-7}\right)}{100 \cdot \left(\frac{1}{250}\right)}\simeq 7.625\times 10^{-6}.
\end{aligned}
$$

光線誤差と回折半径の比が増す。視野による像質の変化を軸上の点だけで予測できない。
:::

:::example[開口を絞って比較する]
前例の視野点を維持し、瞳半径を半分にする。

入力はコマの係数 $B=1\,\mathrm{m^{-2}}$、瞳と像の距離 $L=\frac{1}{10}\,\mathrm m$、瞳半径 $a=\frac{1}{500}\,\mathrm m$、真空波長 $\lambda_0\simeq 5\times 10^{-7}\,\mathrm m$ である。

視野方向の瞳端で固定視野での三次の瞳依存を評価する。

$$
\begin{aligned}
W_{\rm edge}&=B a^{3},\\
W_{\rm edge}/(\mathrm m)&= \left(1\right) \cdot \left(\frac{1}{500}\right)^{3}\simeq 8\times 10^{-9}.
\end{aligned}
$$

位相を一波長当たりの光路差に換算する。

$$
\begin{aligned}
N_W&=\frac{W_{\rm edge}}{\lambda_0},\\
N_W&\simeq \frac{\left(8\times 10^{-9}\right)}{\left(5\times 10^{-7}\right)}=\frac{2}{125}.
\end{aligned}
$$

横端では勾配の三倍の成分が現れる。

$$
\begin{aligned}
\Delta x_{\rm edge}&=3 B L a^{2},\\
\Delta x_{\rm edge}/(\mathrm m)&= 3 \cdot \left(1\right) \cdot \left(\frac{1}{10}\right) \cdot \left(\frac{1}{500}\right)^{2}\simeq 1.2\times 10^{-6}.
\end{aligned}
$$

一様円瞳で二次の勾配を面積平均する。

$$
\begin{aligned}
\bar x&=B L a^{2},\\
\bar x/(\mathrm m)&= \left(1\right) \cdot \left(\frac{1}{10}\right) \cdot \left(\frac{1}{500}\right)^{2}\simeq 4\times 10^{-7}.
\end{aligned}
$$

幾何学的誤差を比較する尺度として無収差の回折半径も求める。

$$
\begin{aligned}
r_{\rm A}&=\frac{61 L \lambda_0}{100 a},\\
r_{\rm A}/(\mathrm m)&\simeq \frac{61 \cdot \left(\frac{1}{10}\right) \cdot \left(5\times 10^{-7}\right)}{100 \cdot \left(\frac{1}{500}\right)}\simeq 1.525\times 10^{-5}.
\end{aligned}
$$

光線誤差は四分の一、回折半径は二倍になる。絞りを変えると両者は逆方向に変化する。
:::

## 色収差とアクロマート

### 分散の相殺
薄レンズの屈折力は波長に依存する。空気中で形状係数を $K_i=1/R_{i1}-1/R_{i2}$ とすれば $\Phi_i(\lambda)=K_i[n_i(\lambda)-1]$ である。基準のスペクトル線を $d$、青と赤の線を $F,C$ とし、
$$
V_i=(n_{id}-1)/(n_{iF}-n_{iC})
$$
を :term[Abbe 数]{en="Abbe number" yomi="あっべすう"} と呼ぶ。接触した薄レンズ二枚では $\Phi_d=\Phi_{1d}+\Phi_{2d}$ で、青と赤の屈折力の差は $\Phi_{1d}/V_1+\Phi_{2d}/V_2$ である。差をゼロにする条件と所要の合成屈折力を連立して
$$
\Phi_{1d}=\frac{V_1}{V_1-V_2}\Phi_d,\qquad
\Phi_{2d}=-\frac{V_2}{V_1-V_2}\Phi_d .
$$
$V_1>V_2>0,\Phi_d>0$ であれば、分散の小さい正レンズと分散の大きい負レンズを組み合わせる。二波長で一致しても中間の波長が一致するとは限らず、残差を二次スペクトルという。

:::example[色消しの配分]
$V_1=60,V_2=30,\Phi_d=5\,\mathrm{m^{-1}}$ なら $\Phi_{1d}=10\,\mathrm{m^{-1}}$、$\Phi_{2d}=-5\,\mathrm{m^{-1}}$ である。青と赤の差は $10/60-5/30=0$ になる。
:::
:::example[横色収差の区別]
物体距離と像面を固定し、系に複数の間隔があると、行列の倍率 $A(\lambda)$ が波長で変わり得る。軸上焦点の一致だけでは、画面周辺の異なる色の位置を一致させる条件は得られない。
:::

### 距離を持つ二枚の一次色消し
空気中で間隔 $d$ の二薄レンズは $\Phi(\lambda)=\Phi_1+\Phi_2-d\Phi_1\Phi_2$。設計波長で微分すると $\Phi'=(1-d\Phi_2)\Phi_1'+(1-d\Phi_1)\Phi_2'$ を得る。これは合成焦点距離の色変化を消す条件を与える。しかし後頂点からの焦点は $(1-d\Phi_1)/\Phi$ なので、分子も波長で変わる。焦点距離を一致させることと、固定した検出器の位置で色を合焦させることは、頂点から離れた複合系では別の設計条件である。

:::example[二色の焦点差を測る]
空気中で固定形状の正薄レンズが $n_d=1.50$ で $f_d=100\,\mathrm{mm}$、青で $n_F=1.51$ とする。屈折力は $n-1$ に比例するので青の焦点は $f_F=100(0.50/0.51)\,\mathrm{mm}\simeq98.0\,\mathrm{mm}$。同じ物距離が十分大きいとき、青の像は基準色より約 $2.0\,\mathrm{mm}$ 手前になる。
:::

### 部分分散から二次スペクトルを計算する

接触薄レンズで第 $i$ 材料の短波長側の :term[相対部分分散]{en="relative partial dispersion" yomi="そうたいぶぶんぶんさん"} を $P_i=(n_{iF}-n_{id})/(n_{iF}-n_{iC})$ とする。形状係数は波長に依存しないので、$F$ 線と $d$ 線での力の差は $\sum_i P_i\Phi_{id}/V_i$ となる。$F,C$ が一致するように二枚の力を決めると、共有の力は $\Phi_F=\Phi_d[1+(P_1-P_2)/(V_1-V_2)]$。従って基準焦点との差は $1/\Phi_F-1/\Phi_d$ である。二枚の材料の $P_i$ が一致すればこのモデルでは三波長が一致する。$V_1-V_2$ が小さい場合には、必要な個々のレンズの力が大きいうえに、同じ部分分散差でも残りの色ずれが増える。

:::example[二色はそろうが基準色が外れる]
部分分散が異なる二材料を接触させる。

入力は第一材料の Abbe 数 $V_{1}=60$、第二材料の Abbe 数 $V_{2}=30$、第一材料の部分分散 $P_{1}=\frac{17}{25}$、第二材料の部分分散 $P_{2}=\frac{7}{10}$、基準の合成屈折力 $\Phi_d=10\,\mathrm{m^{-1}}$ である。

色消しの連立条件から正側の力を決める。

$$
\begin{aligned}
\Phi_{1d}&=\frac{\Phi_d V_{1}}{V_{1} - V_{2}},\\
\Phi_{1d}/(\mathrm{m^{-1}})&= \frac{\left(10\right) \cdot \left(60\right)}{\left(60\right) - \left(30\right)}=20.
\end{aligned}
$$

合成屈折力を保って反対符号の力を選ぶ。

$$
\begin{aligned}
\Phi_{2d}&=- \frac{\Phi_d V_{2}}{V_{1} - V_{2}},\\
\Phi_{2d}/(\mathrm{m^{-1}})&= - \frac{\left(10\right) \cdot \left(30\right)}{\left(60\right) - \left(30\right)}=-10.
\end{aligned}
$$

短波長側の部分分散だけを重み付けして足す。

$$
\begin{aligned}
\Delta\Phi_F&=\frac{P_{1} \Phi_{1d}}{V_{1}} + \frac{P_{2} \Phi_{2d}}{V_{2}},\\
\Delta\Phi_F/(\mathrm{m^{-1}})&= \frac{\left(\frac{17}{25}\right) \cdot \left(20\right)}{\left(60\right)} + \frac{\left(\frac{7}{10}\right) \cdot \left(-10\right)}{\left(30\right)}=- \frac{1}{150}.
\end{aligned}
$$

焦点の位置は合成の力の逆数になる。

$$
\begin{aligned}
f_F&=\frac{1}{\Phi_d + \Delta\Phi_F},\\
f_F/(\mathrm m)&= \frac{1}{\left(10\right) + \left(- \frac{1}{150}\right)}=\frac{150}{1499}.
\end{aligned}
$$

基準の $d$ 線と同一の軸上座標で差を取る。

$$
\begin{aligned}
\Delta f&=f_F - \frac{1}{\Phi_d},\\
\Delta f/(\mathrm m)&= \left(\frac{150}{1499}\right) - \frac{1}{\left(10\right)}\simeq 6.6711\times 10^{-5}.
\end{aligned}
$$

力が小さくなる側に共通の色が動くので、その焦点が基準色より後ろに来る。
:::

:::example[部分分散が一致する理想例]
Abbe 数は前例のまま、部分分散だけを等しくする。

入力は第一材料の Abbe 数 $V_{1}=60$、第二材料の Abbe 数 $V_{2}=30$、第一材料の部分分散 $P_{1}=\frac{69}{100}$、第二材料の部分分散 $P_{2}=\frac{69}{100}$、基準の合成屈折力 $\Phi_d=10\,\mathrm{m^{-1}}$ である。

色消しの連立条件から正側の力を決める。

$$
\begin{aligned}
\Phi_{1d}&=\frac{\Phi_d V_{1}}{V_{1} - V_{2}},\\
\Phi_{1d}/(\mathrm{m^{-1}})&= \frac{\left(10\right) \cdot \left(60\right)}{\left(60\right) - \left(30\right)}=20.
\end{aligned}
$$

合成屈折力を保って反対符号の力を選ぶ。

$$
\begin{aligned}
\Phi_{2d}&=- \frac{\Phi_d V_{2}}{V_{1} - V_{2}},\\
\Phi_{2d}/(\mathrm{m^{-1}})&= - \frac{\left(10\right) \cdot \left(30\right)}{\left(60\right) - \left(30\right)}=-10.
\end{aligned}
$$

短波長側の部分分散だけを重み付けして足す。

$$
\begin{aligned}
\Delta\Phi_F&=\frac{P_{1} \Phi_{1d}}{V_{1}} + \frac{P_{2} \Phi_{2d}}{V_{2}},\\
\Delta\Phi_F/(\mathrm{m^{-1}})&= \frac{\left(\frac{69}{100}\right) \cdot \left(20\right)}{\left(60\right)} + \frac{\left(\frac{69}{100}\right) \cdot \left(-10\right)}{\left(30\right)}=0.
\end{aligned}
$$

焦点の位置は合成の力の逆数になる。

$$
\begin{aligned}
f_F&=\frac{1}{\Phi_d + \Delta\Phi_F},\\
f_F/(\mathrm m)&= \frac{1}{\left(10\right) + \left(0\right)}=\frac{1}{10}.
\end{aligned}
$$

基準の $d$ 線と同一の軸上座標で差を取る。

$$
\begin{aligned}
\Delta f&=f_F - \frac{1}{\Phi_d},\\
\Delta f/(\mathrm m)&= \left(\frac{1}{10}\right) - \frac{1}{\left(10\right)}=0.
\end{aligned}
$$

この三波長は一致する。波長の全範囲で分散が消えることまでは保証しない。
:::

:::example[Abbe 数が近い材料対]
二枚の力を大きくしないと全体を色消しできない組合せを考える。

入力は第一材料の Abbe 数 $V_{1}=45$、第二材料の Abbe 数 $V_{2}=40$、第一材料の部分分散 $P_{1}=\frac{17}{25}$、第二材料の部分分散 $P_{2}=\frac{7}{10}$、基準の合成屈折力 $\Phi_d=10\,\mathrm{m^{-1}}$ である。

色消しの連立条件から正側の力を決める。

$$
\begin{aligned}
\Phi_{1d}&=\frac{\Phi_d V_{1}}{V_{1} - V_{2}},\\
\Phi_{1d}/(\mathrm{m^{-1}})&= \frac{\left(10\right) \cdot \left(45\right)}{\left(45\right) - \left(40\right)}=90.
\end{aligned}
$$

合成屈折力を保って反対符号の力を選ぶ。

$$
\begin{aligned}
\Phi_{2d}&=- \frac{\Phi_d V_{2}}{V_{1} - V_{2}},\\
\Phi_{2d}/(\mathrm{m^{-1}})&= - \frac{\left(10\right) \cdot \left(40\right)}{\left(45\right) - \left(40\right)}=-80.
\end{aligned}
$$

短波長側の部分分散だけを重み付けして足す。

$$
\begin{aligned}
\Delta\Phi_F&=\frac{P_{1} \Phi_{1d}}{V_{1}} + \frac{P_{2} \Phi_{2d}}{V_{2}},\\
\Delta\Phi_F/(\mathrm{m^{-1}})&= \frac{\left(\frac{17}{25}\right) \cdot \left(90\right)}{\left(45\right)} + \frac{\left(\frac{7}{10}\right) \cdot \left(-80\right)}{\left(40\right)}=- \frac{1}{25}.
\end{aligned}
$$

焦点の位置は合成の力の逆数になる。

$$
\begin{aligned}
f_F&=\frac{1}{\Phi_d + \Delta\Phi_F},\\
f_F/(\mathrm m)&= \frac{1}{\left(10\right) + \left(- \frac{1}{25}\right)}=\frac{25}{249}.
\end{aligned}
$$

基準の $d$ 線と同一の軸上座標で差を取る。

$$
\begin{aligned}
\Delta f&=f_F - \frac{1}{\Phi_d},\\
\Delta f/(\mathrm m)&= \left(\frac{25}{249}\right) - \frac{1}{\left(10\right)}=\frac{1}{2490}.
\end{aligned}
$$

同じ部分分散差でも焦点ずれが大きい。個々の屈折力が増すため、単色収差への影響も別に調べる。
:::

## 波面収差と Strehl 比

理想の像点で全瞳からの波を足すとき、波面収差をもつ瞳の各点は位相 $k_0W$ を加える。瞳の実非負振幅を $P$ とし、基準像点が実際のピーク位置になるように傾きや焦点を調整したとすると、同じ振幅瞳での理想ピークに対する :term[Strehl 比]{en="Strehl ratio" yomi="しゅとれーるひ"} は
$$
S=\frac{|\int P e^{ik_0W}\dd^2u|^2}{|\int P\dd^2u|^2}.
\label{eq-optics-strehl}
$$
$P/\int P$ で重みづけた平均を $\langle\cdot\rangle$ と書く。$W$ に定数を加えても絶対値は不変なので $\langle W\rangle=0$ と置ける。$\max|k_0W|\ll1$ なら指数関数を展開して $\langle e^{ik_0W}\rangle=1-k_0^2\langle W^2\rangle/2+O(k_0^3\langle|W|^3\rangle)$ を得る。三次の項は純虚数なので、絶対値二乗には四次から寄与が現れる。従って $S=1-k_0^2\sigma_W^2+O(k_0^4\langle W^4\rangle)\simeq e^{-k_0^2\sigma_W^2}$。後者は二次まで同じ展開を持つ表式であり、一般の大収差への拡張ではない。

:::example[平坦な位相のずれ]
瞳の半分ずつの $W$ が $+w,-w$ なら平均複素振幅は $\cos k_0w$ である。$w=\lambda_0/8$ のとき基準点の強度比は $1/2$ になる。
:::
:::example[RMS 波面誤差]
$\sigma_W=\lambda_0/20$ なら指数近似は $e^{-(2\pi/20)^2}\simeq0.906$ である。ここで用いる長さはピーク間差ではなく、平均を引いた RMS である。
:::

### 一様円瞳のデフォーカス
瞳半径で規格化した座標を $0\le u\le1$ とし $W=w u^2$ と置く。円瞳の面積平均の測度は $2u\dd u$ なので、$t=u^2$ は $[0,1]$ 上の一様分布である。基準像点の振幅は $\int_0^1 e^{ik_0wt}\dd t=e^{ik_0w/2}\operatorname{sinc}(k_0w/2)$。この強度比は $\operatorname{sinc}^2(k_0w/2)$ で、定数を引いた波面誤差の分散は $w^2/12$。小さい $w$ で二乗の展開が $1-k_0^2w^2/12$ となることは :eqref[eq-optics-strehl] の二次展開と一致する。焦点面を動かせるなら、測定面での低下と最良焦点での低下は区別して記す。

:::example[半波のデフォーカス]
$w=\lambda_0/2$ では基準点の強度比が $4/\pi^2\simeq0.405$ となる。二次近似は $1-\pi^2/12\simeq0.178$ で、瞳端の位相が $\pi$ も異なる場合には精度が足りない。
:::

### 収差が混ざると最良の焦点面が動く
一様円瞳で、規格化半径の二乗を $t$ と置く。軸上球面収差 $W_4 t^2$ にデフォーカス $W_2t$ を足す。$t$ は $[0,1]$ の一様分布なので $\operatorname{Var}t=1/12$、$\operatorname{Cov}(t^2,t)=1/4-1/6=1/12$。波面分散を $W_2$ に関して最小にするには $W_2=-W_4$。最小分散は $W_4^2[\operatorname{Var}(t^2)-(\operatorname{Cov}(t^2,t))^2/\operatorname{Var}t]=W_4^2/180$ となる。小収差で Strehl 比を大きくする最良焦点は、二次項を幾何学的にゼロとした近軸焦点から移動する。

### 正弦状の波面誤差で Strehl 比を積分する

振幅の一様な矩形瞳で $W=w\cos(2\pi x/\Lambda)$、瞳幅が $\Lambda$ の整数倍とする。中心の振幅比は周期平均 $\frac1{2\pi}\int_0^{2\pi}e^{iv\cos t}\dd t$、$v=k_0w$ である。指数展開はコンパクトな区間で一様収束するから項ごとに平均できる。$\cos t=(e^{it}+e^{-it})/2$ を二項展開すると偶数べきで等しい指数の対だけが残り、$\langle \cos^{2j}t\rangle=(2j)!/[2^{2j}(j!)^2]$。従って平均振幅は $\sum_{j\ge0}(-1)^j(v/2)^{2j}/(j!)^2=J_0(v)$ である。中心の強度比は $J_0(v)^2$、一方で分散は $w^2/2$ である。小さい $v$ では両者から同じ二次評価が得られる。大きな誤差があると光が中心から横の回折次数へ移るため、中心で求めた比を最大強度の比と呼ぶには中心が最大であることを確認する。

:::example[小さな周期的加工痕]
波面の周期的な誤差が波長よりずっと小さい場合を評価する。

入力は誤差の振幅 $w\simeq 1\times 10^{-8}\,\mathrm m$、真空波長 $\lambda_0\simeq 5\times 10^{-7}\,\mathrm m$ である。

光路誤差を位相へ換算する。

$$
\begin{aligned}
v&=\frac{2 \pi w}{\lambda_0},\\
v/(\mathrm{rad})&\simeq \frac{2 \cdot \pi \cdot \left(1\times 10^{-8}\right)}{\left(5\times 10^{-7}\right)}\simeq 0.12566.
\end{aligned}
$$

正弦の二乗の周期平均は $1/2$ である。

$$
\begin{aligned}
\sigma_W&=\frac{\sqrt{2} w}{2},\\
\sigma_W/(\mathrm m)&\simeq \frac{\sqrt{2} \cdot \left(1\times 10^{-8}\right)}{2}\simeq 7.0711\times 10^{-9}.
\end{aligned}
$$

瞳の位相の周期平均を積分する。

$$
\begin{aligned}
A&=J_{0}\left(v\right),\\
A&\simeq J_{0}\left(\left(0.1256637061\right)\right)\simeq 0.99606.
\end{aligned}
$$

位相誤差を含めた複素振幅の二乗を取る。

$$
\begin{aligned}
S&=A^{2},\\
S&\simeq \left(0.9960560529\right)^{2}\simeq 0.99213.
\end{aligned}
$$

同じ RMS の小収差評価を求める。

$$
\begin{aligned}
S_{\rm M}&=e^{- \frac{v^{2}}{2}},\\
S_{\rm M}&\simeq e^{- \frac{\left(0.1256637061\right)^{2}}{2}}\simeq 0.99214.
\end{aligned}
$$

正確な周期平均と二次近似の違いを比較する。

$$
\begin{aligned}
\Delta S&=S - S_{\rm M},\\
\Delta S&\simeq \left(0.9921276605\right) - \left(0.9921354055\right)\simeq -7.745\times 10^{-6}.
\end{aligned}
$$

二つの強度比はよく一致する。誤差の形より RMS が効く領域である。
:::

:::example[十分には小さくない位相]
同じ瞳で振幅を六倍にする。

入力は誤差の振幅 $w\simeq 6\times 10^{-8}\,\mathrm m$、真空波長 $\lambda_0\simeq 5\times 10^{-7}\,\mathrm m$ である。

光路誤差を位相へ換算する。

$$
\begin{aligned}
v&=\frac{2 \pi w}{\lambda_0},\\
v/(\mathrm{rad})&\simeq \frac{2 \cdot \pi \cdot \left(6\times 10^{-8}\right)}{\left(5\times 10^{-7}\right)}\simeq 0.75398.
\end{aligned}
$$

正弦の二乗の周期平均は $1/2$ である。

$$
\begin{aligned}
\sigma_W&=\frac{\sqrt{2} w}{2},\\
\sigma_W/(\mathrm m)&\simeq \frac{\sqrt{2} \cdot \left(6\times 10^{-8}\right)}{2}\simeq 4.2426\times 10^{-8}.
\end{aligned}
$$

瞳の位相の周期平均を積分する。

$$
\begin{aligned}
A&=J_{0}\left(v\right),\\
A&\simeq J_{0}\left(\left(0.7539822369\right)\right)\simeq 0.86285.
\end{aligned}
$$

位相誤差を含めた複素振幅の二乗を取る。

$$
\begin{aligned}
S&=A^{2},\\
S&\simeq \left(0.8628483467\right)^{2}\simeq 0.74451.
\end{aligned}
$$

同じ RMS の小収差評価を求める。

$$
\begin{aligned}
S_{\rm M}&=e^{- \frac{v^{2}}{2}},\\
S_{\rm M}&\simeq e^{- \frac{\left(0.7539822369\right)^{2}}{2}}\simeq 0.75258.
\end{aligned}
$$

正確な周期平均と二次近似の違いを比較する。

$$
\begin{aligned}
\Delta S&=S - S_{\rm M},\\
\Delta S&\simeq \left(0.7445072695\right) - \left(0.7525825354\right)\simeq -0.0080753.
\end{aligned}
$$

指数近似との差が現れる。仕様が数パーセントの差を問う場合には、瞳全体の複素平均を使う。
:::

:::example[中心が抑制される位相]
誤差の振幅をさらに増し、中心振幅が零になる領域へ近づける。

入力は誤差の振幅 $w\simeq 1.9\times 10^{-7}\,\mathrm m$、真空波長 $\lambda_0\simeq 5\times 10^{-7}\,\mathrm m$ である。

光路誤差を位相へ換算する。

$$
\begin{aligned}
v&=\frac{2 \pi w}{\lambda_0},\\
v/(\mathrm{rad})&\simeq \frac{2 \cdot \pi \cdot \left(1.9\times 10^{-7}\right)}{\left(5\times 10^{-7}\right)}\simeq 2.3876.
\end{aligned}
$$

正弦の二乗の周期平均は $1/2$ である。

$$
\begin{aligned}
\sigma_W&=\frac{\sqrt{2} w}{2},\\
\sigma_W/(\mathrm m)&\simeq \frac{\sqrt{2} \cdot \left(1.9\times 10^{-7}\right)}{2}\simeq 1.3435\times 10^{-7}.
\end{aligned}
$$

瞳の位相の周期平均を積分する。

$$
\begin{aligned}
A&=J_{0}\left(v\right),\\
A&\simeq J_{0}\left(\left(2.387610417\right)\right)\simeq 0.0089689.
\end{aligned}
$$

位相誤差を含めた複素振幅の二乗を取る。

$$
\begin{aligned}
S&=A^{2},\\
S&\simeq \left(0.008968896645\right)^{2}\simeq 8.0441\times 10^{-5}.
\end{aligned}
$$

同じ RMS の小収差評価を求める。

$$
\begin{aligned}
S_{\rm M}&=e^{- \frac{v^{2}}{2}},\\
S_{\rm M}&\simeq e^{- \frac{\left(2.387610417\right)^{2}}{2}}\simeq 0.057825.
\end{aligned}
$$

正確な周期平均と二次近似の違いを比較する。

$$
\begin{aligned}
\Delta S&=S - S_{\rm M},\\
\Delta S&\simeq \left(8.044110703\times 10^{-5}\right) - \left(0.0578245559\right)\simeq -0.057744.
\end{aligned}
$$

二次までの情報だけでは中心強度を予測できない。全像の最大を探すには瞳の Fourier 変換が要る。
:::

## 演習問題

::::exercise{level=basic}
厚いレンズの二面が平面のとき行列を求めよ。

:::solution
$\Phi_1=\Phi_2=0$ なので $L$ は両方とも単位行列であり、行列は $T(d/n_\ell)$ になる。
:::
::::

::::exercise{level=basic}
同一媒質中で主点と節点が一致する理由を示せ。

:::solution
主面間では高さゼロの光線の $p$ が変わらない。両側の $n$ が同じなら角も変わらず、これが節点の条件である。
:::
::::

::::exercise{level=basic}
$f=50\,\mathrm{mm}$、入射瞳径 $25\,\mathrm{mm}$ のとき F 値を求めよ。

:::solution
$N=f/D=50/25=2$ である。
:::
::::

::::exercise{level=basic}
$f_o=1\,\mathrm m,\ f_e=0.05\,\mathrm m$ の Kepler 型望遠鏡の角倍率を求めよ。

:::solution
共有する焦点面で高さが $f_o\alpha$、出射角が $-f_o\alpha/f_e$ なので倍率は $-20$ である。
:::
::::

::::exercise{level=advanced}
$\Phi_1=6,\Phi_2=4\,\mathrm{m^{-1}},e=0.05\,\mathrm m$ とする。$\Phi$ と $e_1,e_2$ を求めよ。

:::solution
$\Phi=10-0.05(24)=8.8\,\mathrm{m^{-1}}$。$1-D=e\Phi_2=0.20$、$1-A=e\Phi_1=0.30$ なので $e_1=0.20/8.8\,\mathrm m$、$e_2=0.30/8.8\,\mathrm m$ である。
:::
::::

::::exercise{level=advanced}
$V_1=50,V_2=25,\Phi_d=8\,\mathrm{m^{-1}}$ を満たす接合色消しレンズの二屈折力を求めよ。

:::solution
色消し条件は $\Phi_{1d}/50+\Phi_{2d}/25=0$。従って $\Phi_{1d}=-2\Phi_{2d}$、和の条件から $\Phi_{2d}=-8\,\mathrm{m^{-1}}$、$\Phi_{1d}=16\,\mathrm{m^{-1}}$ となる。
:::
::::

::::exercise{level=advanced}
$W=E h^3 u_x$ の収差は同じ物点の光線をどのように変えるか。

:::solution
$\grad_{\vb u} W=(Eh^3,0)$ である。すべての瞳位置で一定なので、像点がまとめて動く。変位が $h$ に比例せず $h^3$ に比例するため、視野全体では倍率が位置に依存する。
:::
::::

::::exercise{level=advanced}
基準点での式 :eqref[eq-optics-strehl] が $1$ を超えないことを示せ。

:::solution
$P\ge0$ と $|e^{ik_0W}|=1$ より、積分の三角不等式で $|\int P e^{ik_0W}|\le\int P$。両辺を二乗して正の分母で割ると $S\le1$ である。
:::
::::

::::exercise{level=advanced}
一様円瞳で、瞳の半径を規格化した座標の二乗を $t$ とする。波面差が $W=s t^2+d t$ のとき、小収差で最良の焦点を選ぶ $d$ と、残る波面差の RMS を求めよ。ピストンを除いた RMS とする。

:::solution
円瞳の測度は $t\in[0,1]$ で一様で、$\langle t^j\rangle=1/(j+1)$。よって $\operatorname{Var}(t)=1/12$、$\operatorname{Cov}(t,t^2)=1/12$、$\operatorname{Var}(t^2)=1/5-1/9=4/45$。分散を平方完成すると $(d+s)^2/12+s^2/180$。小収差では Strehl 比を上げるためこの分散を最小にする。従って $d=-s$、残る RMS は $|s|/\sqrt{180}$ である。
:::
::::

## 確認テスト

::::quiz
軸上の通過角を制限する部品を選べ。

- [x] 開口絞り
- [ ] 視野絞り
- [ ] 主面

:::explanation
各面の半径を $a_j$、行列の右上を $B_j$ とすれば、軸上で $|p|\le a_j/|B_j|$ を最も厳しく課す面が開口絞りである。
:::
::::

::::quiz{answer="4" tol="1e-6"}
$f=80\,\mathrm{mm}$、入射瞳径 $20\,\mathrm{mm}$ のとき F 値を答えよ。

:::explanation
焦点距離を入射瞳径で割ると $80/20=4$ になる。
:::
::::

::::quiz
瞳方向の勾配が瞳位置によらない Seidel 収差を選べ。

- [x] 歪曲
- [ ] コマ収差
- [ ] 非点収差

:::explanation
式 :eqref[eq-optics-seidel] の最後の項が瞳に線形である。その勾配は同じ視野点の全光線で一定になる。
:::
::::

::::quiz{answer="1" tol="1e-6"}
波面が瞳全体で一定なら Strehl 比を答えよ。

:::explanation
式 :eqref[eq-optics-strehl] の分子の積分から定数位相がくくり出せる。絶対値では消えるため分子と分母が一致する。
:::
::::

---
title: "超伝導の微視的理論"
description: "Cooper 対の形成から BCS ギャップ、熱力学、Nambu 表示と磁場の排除を導く。"
notation:
  - symbol: 'g,\omega_D'
    meaning: 対形成の引力とエネルギー切断
  - symbol: '\Delta,E_k'
    meaning: ギャップと準粒子エネルギー
  - symbol: '\nu_0'
    meaning: 一スピンの Fermi 面状態密度
  - symbol: '\mathcal G_N'
    meaning: Nambu 行列の Green 関数
---

## この章の目標

- Cooper 問題と有限密度の BCS 変分の違いを説明する。
- Bogoliubov 変換を実行し、ギャップと転移温度を求める。
- 超伝導の自由エネルギーと比熱を計算する。
- Nambu 表示から対振幅と電磁応答の関係を示す。

## Cooper 不安定性

### Fermi 海の上の二粒子

[固体の超伝導](/physics/solid-state/superconductivity/) で扱う零抵抗と磁束排除を微視的に説明するため、まず電子の対ができる条件を調べる。埋まった Fermi 海を固定し、その上の全運動量零の異スピン二電子だけを変分する。$0<\xi_k<\omega_D$ に限り、行列要素 $-g/V$ の引力を仮定する。束縛エネルギーを $B>0$ とすれば二体方程式は $(2\xi_k+B)\varphi_k=(g/V)\sum_p\varphi_p$ である。一スピンあたり単位体積の状態密度 $\nu_0$ を一定とすると、非零解の条件は
$$
1=g\nu_0\int_0^{\omega_D}\frac{\dd\xi}{2\xi+B}
=\frac{g\nu_0}{2}\log\frac{2\omega_D+B}{B},\qquad
B=\frac{2\omega_D}{e^{2/(g\nu_0)}-1}.
\label{eq-mb-cooper-binding}
$$
$B\to0^+$ で積分が発散するので任意の正の $g$ に対し解がある。この対形成の発散を :term[Cooper 不安定性]{en="Cooper instability" yomi="くーぱーふあんていせい"} と呼ぶ。対の波動関数は $\varphi_k\propto(2\xi_k+B)^{-1}$ であり、殻の下端に集中する。Fermi 海自体の再配置を含めた多体系のギャップは、二粒子の $B$ と同じ量ではない。

:::example[弱引力の指数]
$g\nu_0=1/5$、$\omega_D=20\,\mathrm{meV}$ なら $B=40/(e^{10}-1)\,\mathrm{meV}\simeq0.00182\,\mathrm{meV}$ である。有限箱で準位間隔がこの値より大きいと、連続積分を用いる条件が崩れる。
:::

:::example[束縛波動関数の殻内での比]
$B=\omega_D/2$ の Cooper 対なら、$\xi=0$ の振幅に対し殻の端 $\xi=\omega_D$ の振幅は $B/(2\omega_D+B)=1/5$ である。同じ数の一粒子状態を含む小区間どうしなら、端の重みは下端の $1/25$ の桁になる。
:::

### 束縛波動関数の規格化

対の運動量を零に固定し、Fermi 海の上で $0<\xi<\omega_D$ の状態だけを残す。行列要素を $-g/V$ とすると、固有方程式は $(2\xi+B)\varphi(\xi)=(g/V)\sum_{\vb k'}\varphi(\xi_{\vb k'})$ であり、束縛エネルギー $B>0$ で $\varphi(\xi)=C/(2\xi+B)$ となる。片スピンの状態密度 $V\nu_0$ を一定とした規格化は
$$
1=V\nu_0|C|^2\int_0^{\omega_D}\frac{\dd\xi}{(2\xi+B)^2}
=\frac{V\nu_0|C|^2}{2}\left(\frac1B-\frac1{2\omega_D+B}\right).
$$
従って弱結合の $B\ll\omega_D$ では $|C|^2\simeq2B/(V\nu_0)$ であり、確率の多くは Fermi 海のすぐ外の幅 $O(B)$ に集まる。この対の指数尺度を、海の占有も変分する BCS ギャップの指数と区別する。

:::example[上の半分にある確率]
$B=2\omega_D$ とし、$\xi>\omega_D/2$ の確率を上の原始関数で割り算すると
$(1/(3\omega_D)-1/(4\omega_D))/(1/(2\omega_D)-1/(4\omega_D))=1/3$ である。$B\ll\omega_D$ では同じ比は小さくなり、Fermi 海の境界への集中を確かめられる。
:::

### 有限箱の対感受率と熱力学極限

Cooper 問題でエネルギーを $2\mu$ から測る。離散スペクトルでその直上に厳密な零の軌道がないとすれば、
$$
C_L(B)=\frac1V\sum_{0<\xi_k<\omega_D}\frac1{2\xi_k+B}
$$
は $B\ge0$ で有限であり、$B>0$ の束縛解がある条件は $g C_L(0)>1$ である。これは $C_L$ が正の $B$ で狭義単調に減り、無限大で零に行くことから従う。有限箱の十分弱い引力でも最低の自由対準位 $2\xi_{\min}$ に対しては下がる固有値がある。実際 $B$ を負の領域 $-2\xi_{\min}<B<0$ まで許すと、最小準位の項が端で発散する。この根は $2\mu$ より下とは限らない。どの基準より下がることを束縛と呼ぶかを明記すれば、有限箱の見かけの閾値と矛盾しない。

等間隔 $\xi_j=j\delta$、$j=1,\dots,M$ とし、$1/V=\nu_0\delta$ になる一組の状態の縮約模型では $C_L(0)=(\nu_0/2)\sum_{j=1}^M 1/j$ である。$M=\omega_D/\delta$ を大きくすると和は $\log M+\gamma_E+O(M^{-1})$ になり、$g C_L(0)>1$ の閾値は零に近づく。無限体積でどれほど弱い引力でも対の解があるのは、有限箱の閾値を越えた極限として理解できる。

:::example[離散準位二つの閾値]
$\xi_j=j\delta$、$M=2$ なら $C_L(0)=3\nu_0/4$ である。$g\nu_0=1$ では $B>0$ の解がなく、$g\nu_0=2$ ではある。後者では $1=2\delta[(2\delta+B)^{-1}+(4\delta+B)^{-1}]$ を解くと $B/\delta=\sqrt5-1$ となる。分母が正の枝を取った。
:::

### 正常相の熱分布で対の発散を確認する

二粒子の束縛問題に対し、正常な多体状態が微弱な対外場へ応答する量は、Fermi 面の両側を含む。$b_{\vb k}=c_{-\vb k\downarrow}c_{\vb k\uparrow}$ とし、自由な平均で対を作って戻す虚時間相関を計算する。$0<\tau<\beta$ では $\langle b(\tau)b^\dagger(0)\rangle=e^{-2\xi\tau}[1-f(\xi)]^2$ であり、$\tau$ を積分すると
$$
[1-f(\xi)]^2\frac{1-e^{-2\beta\xi}}{2\xi}=\frac{1-2f(\xi)}{2\xi}.
$$
$\xi=0$ は連続に $1/(4T)$ と読む。分子は空の二軌道に対を足す効果と、満ちた二軌道から取り除く効果の差である。殻全体を足せば $\nu_0\int_0^{\omega_D}\tanh(\xi/(2T))\dd\xi/\xi$ になる。対外場への一次応答と相互作用が作る自己無撞着な外場を合わせると、この感受率に $g$ を掛けたものが $1$ に達した点で正常解が不安定になる。従って Fermi 海を固定した Cooper 問題の $\nu_0/2$ と、海の両側を変分する多体の対応答の $\nu_0$ の差は、状態の数え方の差から出る。

:::example[符号が反対でも同じ対応答]
$\xi=E$ と $\xi=-E$ の二組では $\tanh(\xi/(2T))/(2\xi)$ は同じ正の値である。$\beta E=\log3$ ならそれぞれ $1/(4E)$ で、両方を残せば $1/(2E)$ となる。
:::

## BCS 波動関数と変分計算

### 対の重ね合わせと変分

縮約 BCS Hamiltonian を $K=\sum_{k\sigma}\xi_k n_{k\sigma}-(g/V)\sum_{kp}b_k^\dagger b_p$、$b_k=c_{-k\downarrow}c_{k\uparrow}$ とする。対を作る範囲は $|\xi|<\omega_D$ である。各 $(k\uparrow,-k\downarrow)$ が一組であり、全ての $k$ を走っても各スピン軌道は一度だけ使う。:term[BCS 波動関数]{en="BCS wave function" yomi="びーしーえすはどうかんすう"} は
$$
\ket{\Psi}=\prod_{k}(u_k+v_k b_k^\dagger)\ket0,\quad
|u_k|^2+|v_k|^2=1,\qquad
\Delta=\frac gV\sum_k u_k^* v_k .
\label{eq-mb-bcs-state}
$$
対のない殻の外では自由占有を用いる。$p\ne k$ の $\langle b_k^\dagger b_p\rangle=(u_k^*v_k)^*u_p^*v_p$ が積に分かれ、対角でのずれは $O(V^0)$ の全エネルギーなので、熱力学極限の変分密度から除く。すると $E=2\sum\xi |v|^2-V|\Delta|^2/g$ となる。共通の位相を選び $u=\cos(\theta/2),v=\sin(\theta/2)$、$\Delta>0$ とする。$\theta_k$ を微分すると $\xi_k\sin\theta_k-\Delta\cos\theta_k=0$ である。各組の有効エネルギー $\xi_k(1-\cos\theta_k)-\Delta\sin\theta_k$ を下げる枝を選び、
$$
E_k=\sqrt{\xi_k^2+|\Delta|^2},\qquad
u_k^2=\frac{1+\xi_k/E_k}{2},\quad |v_k|^2=\frac{1-\xi_k/E_k}{2},
\quad u_k^*v_k=\frac{\Delta}{2E_k}
\label{eq-mb-bcs-coherence}
$$
を得る。粒子数の期待値は $2\sum |v_k|^2$ であり、指定密度では $\mu$ も調節する。一組の粒子数は $0,2$ の二値を取るため、分散は $4u_k^2|v_k|^2$ である。

:::example[殻の中心と端]
$\xi=0$ では $u^2=|v|^2=1/2$、一組の分散は $1$ である。$\xi=3\,\mathrm{meV}$、$\Delta=4\,\mathrm{meV}$ では $E=5\,\mathrm{meV}$、$u^2=4/5$、$|v|^2=1/5$ である。
:::

:::example[二組の粒子数分布]
二組で $|v_1|^2=|v_2|^2=1/4$ とする。電子数 $0,2,4$ の確率は順に $9/16,6/16,1/16$ であり、平均は $1$ である。分散は一組の $4(1/4)(3/4)=3/4$ を足して $3/2$ になる。
:::

### 固定粒子数の成分を取り出す

$P_k^\dagger=c_{k\uparrow}^\dagger c_{-k\downarrow}^\dagger$ とし、別の $k$ の対演算子は四回の交換で可換、各 $(P_k^\dagger)^2=0$ である。$u_k\ne0$ なら
$$
\prod_k(u_k+v_kP_k^\dagger)\ket0
=\left(\prod_k u_k\right)\exp\left(\sum_k\frac{v_k}{u_k}P_k^\dagger\right)\ket0 .
$$
$2M$ 粒子への射影は指数の $M$ 次だけを取る。従って全体の粒子数を固定しても、どの $M$ 本の対軌道が占有されるかの振幅は $v_k/u_k$ の積で重みづけされる。固定粒子数の状態では $\langle c c\rangle=0$ になるが、$\langle P_k^\dagger P_{k'}\rangle$ は粒子数を保つので対の干渉を測れる。

:::example[二対軌道で一対を射影する]
両方で $u=v=1/\sqrt2$ とすると粒子数 $0,2,4$ の確率は $1/4,1/2,1/4$ である。二粒子へ射影して規格化すると $(P_1^\dagger+P_2^\dagger)\ket0/\sqrt2$ であり、$\langle P_1^\dagger P_2\rangle=1/2$ となる。どちらか一方を確率 $1/2$ で占有する混合状態ではこの非対角成分が消える。
:::

### 対の位相をそろえる変分の条件

殻内の $v_k$ に位相 $\phi_k$ を許す。運動エネルギーは $|v_k|$ だけで決まる。相互作用は $-g/V$ に $\left|\sum_k |u_kv_k|e^{i\phi_k}\right|^2$ を掛けたものなので、三角不等式により $u_kv_k\ne0$ の組の位相がすべて同じときに最小になる。すべての位相を同時に動かしてもエネルギーは変わらない。この自由度が対の $U(1)$ 位相である。一方、引力が一定でなく符号の違う行列要素を持つ場合には、三角不等式だけから位相をそろえる結論は得られない。角度依存する対は [超伝導の発展的話題](/physics/many-body/superconductivity-advanced/) で扱う。

組の偶数占有部門を二準位として扱うと、正準演算子を取り替える意味が明確になる。$s^+=b^\dagger$、$s^-=b$、$s^z=(n_\uparrow+n_\downarrow-1)/2$ とすると、空状態と対状態で $s^z$ は $-1/2$ と $1/2$ である。$[s^+,s^-]=2s^z$ は $b^\dagger b=n_\uparrow n_\downarrow$ と $b b^\dagger=(1-n_\uparrow)(1-n_\downarrow)$ の差から従う。$\Delta$ を実に選べば組の平均場は $\xi+2\xi s^z-2\Delta s^x$ であり、最低状態の擬スピンは有効場 $(-\Delta,0,\xi)$ の逆向きである。したがって $\langle s^x\rangle=\Delta/(2E)$、$\langle s^z\rangle=-\xi/(2E)$ となる。この :term[Anderson 擬スピン]{en="Anderson pseudospin" yomi="あんだーそんぎすぴん"} は軌道の空と対を表し、一電子の本来のスピンを測るものではない。

:::example[位相の打ち消しによるエネルギー差]
二組とも $|u_k v_k|=1/2$ とする。位相差が零なら $\left|\sum_k u_k^*v_k\right|^2=1$、位相差が $\pi$ なら零である。縮約 Hamiltonian の二組の対角期待値はどちらの位相でも同じであり、エネルギー差は非対角の対移動による $g/V$ である。
:::

### 二粒子密度行列に残る巨視的固有値

位相を選んだ BCS の積状態で $P_{kp}=\langle b_k^\dagger b_p\rangle$ とする。非対角は $\kappa_k^*\kappa_p$、$\kappa_k=u_k^*v_k$ である。対角は $|v_k|^2$、階数一の表式との差は $|v_k|^2-|\kappa_k|^2=|v_k|^4$ である。対角補正の作用素ノルムは $1$ 以下なので、最大固有値は $\sum_k|\kappa_k|^2$ と高々 $1$ だけ異なる。弱結合でこの和は
$$
\frac1V\sum_k|\kappa_k|^2=\frac{\nu_0\Delta}{2}\arctan\frac{\omega_D}{\Delta}
 \simeq\frac{\pi\nu_0\Delta}{4}
$$
である。これは対角補正に比べ体積に比例して大きい。粒子数を変えない観測量だけでも、多数の対軌道にまたがるコヒーレンスを示せる。この計算で $\sum |v|^2$ は深い占有側も数えるので、相関行列の巨視的固有値とは別の量である。

:::example[殻の中心だけを残す模型]
$\xi_k=0$ の $M$ 組で $u=v=1/\sqrt2$ とする。このとき $P$ のすべての非対角要素は $1/4$、対角は $1/2$ である。一様なベクトルの固有値は $(M+1)/4$、それに直交する空間の固有値は $1/4$ である。一方、位相を各組で独立に平均した混合状態では非対角が零で、最大固有値は $1/2$ のままである。
:::

## Bogoliubov–de Gennes 変換

### Nambu スピノルと正準変換

対の殻での平均場は $K_{\rm MF}=\sum_k[\xi_k(n_{k\uparrow}+n_{-k\downarrow})-\Delta b_k^\dagger-\Delta^*b_k]+V|\Delta|^2/g$ である。$\Psi_k=(c_{k\uparrow},c_{-k\downarrow}^\dagger)^T$ とすると、この二次部分は $\Psi_k^\dagger h_k\Psi_k+\xi_k$ で、
$$
h_k=\begin{pmatrix}\xi_k&-\Delta\\-\Delta^*&-\xi_k\end{pmatrix},
\quad h_k^2=E_k^2 I.
\label{eq-mb-bdg}
$$
Nambu の第二成分の正規順序から定数 $\xi_k$ が出る。実正の $\Delta$ では列 $(u,-v)^T,(v,u)^T$ が固有値 $E,-E$ の正規化固有ベクトルである。直交変換を $\Psi_k$ に施し $(\gamma_{k\uparrow},\gamma_{-k\downarrow}^\dagger)^T$ を得る。変換がユニタリなので正準反交換関係を保つ。負の固有値の演算子を正規順序に戻すと
$$
K_{\rm MF}=V|\Delta|^2/g+\sum_k(\xi_k-E_k)
+\sum_{k,\sigma}E_k\gamma_{k\sigma}^\dagger\gamma_{k\sigma}.
$$
殻の外の自由電子の項は以降の正常相との差で相殺する。準粒子真空は各組で $u+v b^\dagger$ になり、変分状態と一致する。非一様な場合は位置の一体核と対ポテンシャルをブロックに持つ固有値問題になり、これを :term[Bogoliubov–de Gennes 方程式]{en="Bogoliubov–de Gennes equation" yomi="ぼごりゅーぼふどじぇんぬほうていしき"} と呼ぶ。

:::example[粒子と正孔の混合]
$\xi=0,\Delta>0$ では固有ベクトルが $(1,-1)^T/\sqrt2,(1,1)^T/\sqrt2$ である。正の励起は電子だけを追加する操作と異なり、電子の追加と正孔の追加が同じ振幅で入る。
:::

:::example[非零の分散での固有ベクトル]
$\xi=5,\Delta=12$ なら $E=13$、$u^2=9/13,|v|^2=4/13$ となる。正の固有ベクトルは $(3,-2)^T/\sqrt{13}$ である。実際に $h(3,-2)^T=(39,-26)^T$ となり、固有値 $13$ を確認できる。
:::

### 行列を射影してコヒーレンス因子を読む

一対の Hermite な BdG 行列を $M=\begin{pmatrix}\xi&-\Delta\\-\bar\Delta&-\xi\end{pmatrix}$ とする。二乗すると混成項が相殺して $M^2=E^2 I$、$E=\sqrt{\xi^2+|\Delta|^2}$ となる。$E>0$ なら $P_\pm=(I\pm M/E)/2$ は互いに直交する射影であり、$MP_\pm=\pm EP_\pm$ である。従って固有ベクトルの位相を任意に選ばなくても、対角の粒子重みは $(1\pm\xi/E)/2$、非対角は $\mp\Delta/(2E)$ と読める。$E=0$ の縮退点では個別の射影は分離せず、二つの合計 $I$ が定義できる。

:::example[粒子優勢な枝]
$\xi=3\,\mathrm{meV}$、$\Delta=4\,\mathrm{meV}$ とすると $E=5\,\mathrm{meV}$、正の枝の粒子重みは $4/5$、負の枝では $1/5$ である。異常成分は正負の枝で $-2/5,2/5$ となる。全周波数で異常成分の重みを足すと零になるのは、二つの消滅演算子の等時刻の反交換子が零だからである。
:::

### 複素ギャップと粒子–正孔の相手

$u=\sqrt{(1+\xi/E)/2}\ge0$、$v=(\Delta/|\Delta|)\sqrt{(1-\xi/E)/2}$ とし、$\Delta\ne0$ を取る。実のギャップに使った固有ベクトルは、複素の場合には $(u,-v^*)^T$ と $(v,u)^T$ になる。従って
$$
c_{k\uparrow}=u\gamma_{k\uparrow}+v\gamma_{-k\downarrow}^\dagger,\qquad
c_{-k\downarrow}^\dagger=-v^*\gamma_{k\uparrow}+u\gamma_{-k\downarrow}^\dagger .
$$
正準性には列の直交に加えて、もとの二つの成分の反交換関係を使う。逆変換の $\gamma_{k\uparrow}=u c_{k\uparrow}-v c_{-k\downarrow}^\dagger$ を $u|0\rangle+v b_k^\dagger|0\rangle$ に作用させると二項は相殺する。相手の $\gamma_{-k\downarrow}$ に対しても、$c_{-k\downarrow}b_k^\dagger|0\rangle=-c_{k\uparrow}^\dagger|0\rangle$ から同じ結果を得る。

完全なスピン付き Nambu 空間では作成演算子と消滅演算子を重複して表すので、二次形式を全空間に書く際に $1/2$ を入れ、正規順序の定数も戻す。一方、本文の二成分 $(c_{k\uparrow},c_{-k\downarrow}^\dagger)$ を全 $k$ で用いれば各物理軌道は一つの組にだけ現れる。第二成分を正規順序に直せば正の励起の数は二スピン分である。粒子–正孔対称性で現れる負エネルギーを、さらに独立な熱励起として数えると自由エネルギーの係数を誤る。

:::example[ギャップの位相が純虚数の場合]
$\xi=0$、$\Delta=iD$、$D>0$ なら正固有ベクトルは $(1,i)^T/\sqrt2$ である。第 1 成分へ行列を作用させると $(-iD)i/\sqrt2=D/\sqrt2$ となる。位相は異常成分の符号へ現れるが、粒子を足すスペクトルの重みは実ギャップと同じ $1/2$ である。
:::

### 非一様な対と基底の取り方

連続の一体作用素を $h_0=-\nabla^2/(2m)+U-\mu$ とし、スピンに依らない実の外部ポテンシャルと一重項の局所対ポテンシャルを取る。空間の Nambu Hamiltonian は
$$
\mathcal H_{\rm BdG}=
\begin{pmatrix}
 h_0&-\Delta(\vb r)\\
 -\Delta^*(\vb r)&-h_0^*
\end{pmatrix}
$$
となる。位置の複素共役を含むのは、第二成分が正孔に対応し一体作用素の転置が入るためである。磁場を含む場合、一体作用素に入るベクトルポテンシャルの複素共役で電荷の符号が反転する。有限箱で同じ自己共役境界条件を持つ二成分に作用させ、固有ベクトルを $\int\dd^3r(|u|^2+|v|^2)=1$ で規格化する。全スピン空間なら粒子–正孔共役はエネルギーの符号を反転する。上の縮約一重項ブロックで、$h_0$ がスピンに依らず局所 $\Delta$ が転置対称なら、$(u,v)^T$ に対応する負の解を $(-v^*,u^*)^T$ と取れる。これは行列を直接作用させて確かめられる。場を解の基底に展開した後、$\Delta(\vb r)=g\langle\psi_\downarrow(\vb r)\psi_\uparrow(\vb r)\rangle$ を満たすよう対ポテンシャルを更新する。局所相互作用を全波数で使う場合は、有限の殻か散乱長による正則化を別に指定する。

:::example[局所的な外部ポテンシャルの二つの符号]
$\Delta=0$ とし、一体作用素に実数 $U_0$ を一様に加える。空間一様な準粒子の固有値は $\pm(\epsilon_k+U_0-\mu)$ である。対角二成分に同じ $+U_0$ を入れると粒子–正孔の写像が破れ、単に化学ポテンシャルを動かした結果と一致しない。
:::

## ギャップ方程式と転移温度

### ギャップ方程式

準粒子を Fermi 分布で占有すると $\langle b_k\rangle=\Delta[1-2f(E_k)]/(2E_k)$ である。変換を代入すると二つの余分な準粒子占有が対振幅をそれぞれ減らす。このため
$$
\frac1{g\nu_0}=\int_0^{\omega_D}\frac{\tanh(E/2T)}{E}\dd\xi,\qquad E=\sqrt{\xi^2+|\Delta|^2}.
\label{eq-mb-bcs-gap}
$$
$T=0$ の積分は $\operatorname{arsinh}(\omega_D/|\Delta_0|)$ であり、$\Delta_0=\omega_D/\sinh[1/(g\nu_0)]\simeq2\omega_De^{-1/(g\nu_0)}$ である。$T_c$ では $\Delta=0$ と置く。$X=\omega_D/T_c\gg1$ なら
$$
\int_0^X\frac{\tanh(x/2)}x\dd x
=\log X+\log(2e^{\gamma_E}/\pi)+o(1).
$$
定数を得るには部分積分し、$\log X-\int_0^\infty \log x\,(\partial_x\tanh(x/2))\dd x$ を評価する。積分のモーメントは $2\Gamma(s+1)(1-2^{1-s})\zeta(s)$ である。まず $\Re s>1$ で $\partial_x\tanh(x/2)=2\sum_{n\ge1}(-1)^{n-1}ne^{-nx}$ を項別積分し、収束するモーメントの領域 $\Re s>-1$ に延長する。$s=0$ で微分し、$\Gamma'(1)=-\gamma_E$、$\zeta(0)=-1/2$、$\zeta'(0)=-\tfrac12\log(2\pi)$ を用いると上の定数になる。従って
$$
T_c\simeq\frac{2e^{\gamma_E}}{\pi}\omega_De^{-1/(g\nu_0)},\qquad
\frac{2\Delta_0}{T_c}\simeq2\pi e^{-\gamma_E}\simeq3.53 .
$$
$T$ を K で表す場合は分母を $k_BT_c$ と読む。$g\nu_0=0.2$、$\omega_D=20\,\mathrm{meV}$ なら $\Delta_0\simeq0.270\,\mathrm{meV}$、$T_c\simeq1.77\,\mathrm K$ である。ギャップの存在温度は、エネルギーに依存する結合や状態密度ではこの数値と異なる。

### 粒子数の自己無撞着な決定

有限温度の Bogoliubov 変換で、上スピンの占有は $u^2 f(E)+|v|^2[1-f(E)]$ である。下スピンも同じなので、殻内の全粒子数は $N_{\rm shell}=\sum_k[1-\xi_k\tanh(E_k/(2T))/E_k]$ になる。殻外の占有を加え、指定した $N$ に等しくすることで $\mu$ を決める。$\nu_0$ が一定で殻が $\xi=0$ の周りに対称なら、奇関数の項は積分で消える。これが弱結合の対称殻でギャップ方程式を粒子数の式から独立に評価できる理由である。

:::example[一組の有限温度占有]
$\xi/E=3/5$、$f(E)=1/4$ なら一スピンの占有は $(4/5)(1/4)+(1/5)(3/4)=7/20$ であり、一組で $7/10$ である。零温度では各スピンで $|v|^2=1/5$ に戻る。
:::

:::example[ギャップの切断依存]
同じ $g\nu_0$ で有効な切断 $\omega_D$ を二倍にした場合、零温度の厳密な殻積分の解でも $\Delta_0$ は二倍になる。弱結合の転移温度も同じ倍率なので、比 $2\Delta_0/T_c$ は変わらない。状態密度の変化を同時に入れた比較にはこの比例を使えない。
:::

### 対角化後のギャップ方程式の一意性

一定のカットオフと引力を取る。正の温度でギャップ積分の被積分関数は $F(E)=\tanh(E/(2T))/(2E)$、$E=\sqrt{\xi^2+|\Delta|^2}$ である。$x=E/(2T)$ とすると $F'(E)$ の符号は $x\operatorname{sech}^2 x-\tanh x$ の符号である。後者は $x=0$ で零、その微分は $-2x\operatorname{sech}^2x\tanh x<0$ なので $E>0$ で負になる。従って各 $T>0$ でギャップ積分は $|\Delta|$ に関して厳密に減少し、無限大では零である。$|\Delta|=0$ での値は $T$ を上げると減少するので、線形化した方程式が等号になる温度の下にだけ正の解が一つある。対の位相はこの方程式から決まらない。

:::example[臨界での極限が有限である理由]
$\xi=0$ では $\tanh(E/(2T))/(2E)$ は $E\downarrow0$ で $1/(4T)$ へ近づく。例えば $T=2\,\mathrm{meV}$ なら $1/8\,\mathrm{meV^{-1}}$ である。零温度では同じ点の極限が発散する。従って $T$ の微小な値を先に零へ置くと、線形化した積分の赤外切断を見失う。
:::

### 散乱長で正則化した連続模型

殻に限定された金属とは別に、三次元連続空間の二成分ガスで短距離相互作用を使う場合を考える。相互作用を $U_0<0$ と書く。二体の相対運動の $T$ 行列を $T^{-1}(z)=U_0^{-1}-\int_{\vb k}(z-2\epsilon_k)^{-1}$ とすると、$T(0)=4\pi a_s/m$ の定義から
$$
U_0^{-1}=\frac{m}{4\pi a_s}-\int_{\vb k}^{\Lambda}\frac1{2\epsilon_k}.
$$
ここで $\epsilon_k=k^2/(2m)$、$a_s$ は散乱長、$\int_{\vb k}=\int\dd^3 k/(2\pi)^3$ である。同じ切断を用いて $-1/U_0=\int_{\vb k}^{\Lambda}(1-2f(E_k))/(2E_k)$ に代入し、$\Lambda$ を外すと
$$
\frac m{4\pi a_s}=\int_{\vb k}\left[\frac1{2\epsilon_k}-\frac{1-2f(E_k)}{2E_k}\right],
\qquad
n=\int_{\vb k}\left[1-\frac{\epsilon_k-\mu}{E_k}(1-2f(E_k))\right].
$$
大きな $k$ で二項の先頭は相殺し、粒子数も占有が $O(k^{-4})$ で収束する。この模型ではギャップと化学ポテンシャルを同時に求めるため、一定状態密度の対称殻の式は用いない。

零密度の二体極は $a_s>0$ で $B=1/(m a_s^2)$ に現れる。実際、$T(-B)^{-1}=m/(4\pi a_s)-\int_{\vb k}[(2\epsilon_k)^{-1}-(2\epsilon_k+B)^{-1}]$ であり、半径方向の積分は $m\sqrt{m B}/(4\pi)$ である。分母が零のとき $1/a_s=\sqrt{mB}>0$ となる。零温度のギャップ方程式で密度と $\Delta$ を小さくし、$\mu<0$ を仮定すると $E_k\to\epsilon_k+|\mu|$ である。上の二体条件と同じ積分になり、$2|\mu|=B$ へ近づく。粒子あたりの化学ポテンシャルが束縛対のエネルギーの半分になることを表す。この平均場の補間は密な強相互作用での相関まで正確に与えるものではない。

:::example[散乱長を二倍にした二体極]
正の $a_s$ の浅い二体状態では $B\propto a_s^{-2}$ なので、散乱長が二倍なら束縛エネルギーは四分の一になる。一方、負の $a_s$ ではこの二体極は存在しない。有限密度の Fermi 面があれば、負の $a_s$ の弱い引力でも Cooper の対数が現れる。
:::

### 同位体効果に入る仮定

格子振動のエネルギー切断を $\omega_D\propto M^{-1/2}$ とし、無次元引力 $g\nu_0$ は同位体で一定と仮定する。このとき弱結合の転移温度から $-\partial\log T_c/\partial\log M=1/2$ である。調和近似でばね定数を固定すると、フォノンによる静的な結合の $|g_{\rm ep}|^2/\omega_D\propto 1/(M\omega_D^2)$ が一定になるので、この仮定の一例が得られる。Coulomb 擬ポテンシャル $\mu^*(M)=\mu/[1+\mu\log(W/\omega_D)]$ を引く場合、引力 $\lambda$ が質量に依らない二殻模型では $T_c\propto\omega_D e^{-1/(\lambda-\mu^*)}$ となる。$\partial_{\log\omega_D}\mu^*=(\mu^*)^2$ を使うと同位体指数は
$$
\alpha_{\rm iso}=\frac12\left[1-\frac{(\mu^*)^2}{(\lambda-\mu^*)^2}\right].
$$
$\lambda>\mu^*$ かつ弱結合で質量くりこみを単純化した場合の式である。従って指数が $1/2$ から外れることだけで、格子の寄与の有無は判定できない。

:::example[反発を入れた指数]
$\lambda=3/10$、$\mu^*=1/10$ なら、上の弱結合模型で $\alpha_{\rm iso}=3/8$ である。反発を零にすれば $1/2$ に戻る。これは基準にした電子帯域と裸の反発、ばね定数を固定した比較である。
:::

### 対の有限運動量で正常相の対数を切る

正常でギャップのない低温の金属で、対の中心の全運動量を $\vb q$ とする。$q\ll p_F$、殻が帯域より十分狭く、粒子–正孔について対称なものとする。二つの一体エネルギーは $\xi_+=\xi+\delta$、$\xi_-=\xi-\delta$、$\delta=v_F\hat k\cdot\vb q/2$ である。二軌道の虚時間の積を前節と同じように積分すると、線形の対の感受率の一軌道の因子は
$$
\frac{[1-f(\xi_+)][1-f(\xi_-)](1-e^{-\beta(\xi_++\xi_-)})}{\xi_++\xi_-}
=\frac{1-f(\xi_+)-f(\xi_-)}{2\xi}.
$$
$T=0$ では $|\xi|<|\delta|$ で一方だけが占有されており、分子は零である。従って $\xi$ の対数の切断が $|\delta|$ に移る。角度を平均しても、三次元等方で $\int_0^1\log(1/x)\dd x=1$ が有限なので、$q\ne0$ の対数は有限である。$\delta_q=v_Fq/2\ll\omega_D$ とすると、殻での感受率は $\nu_0[\log(\omega_D/\delta_q)+1]$ である。動径の殻を線形化し、両軌道が同じ対称な切断に入るかによる端の誤差を無視した値である。有限温度なら $T$ による切断も競合する。位相を持たない一定の引力では、平衡の対がまず $\vb q=0$ を選ぶ理由の一つである。

Zeeman のエネルギーを上下に $\mp h$ とすると、同じ式の分子に入るずれが $\delta-h$ になる。$q$ を有限にしたとき、一部の方向ではこのずれを相殺できる。全角度でそれができるわけではなく、有限運動量で最も強い対の不安定性があるかは角度平均で判定する。後の Pauli 限界の比較は、こうした有限運動量の対を変分に入れない二つの一様相の比較である。

### 低温でのギャップの補正を解く

$\Delta=\Delta_0-\delta$、$\delta\ll\Delta_0$ とし、弱結合なので $\omega_D\gg\Delta_0$ とする。ギャップ方程式を零温度のものから引くと
$$
\operatorname{arsinh}(\omega_D/\Delta)-\operatorname{arsinh}(\omega_D/\Delta_0)
=2\int_0^{\omega_D}\frac{f(\sqrt{\xi^2+\Delta^2})}{\sqrt{\xi^2+\Delta^2}}\dd\xi.
$$
左辺は $\delta/\Delta_0$ が先頭である。右辺は Boltzmann 因子が小さいので、先頭の式では $\Delta_0$ を用いる。小さい $\xi$ で $E\simeq\Delta_0+\xi^2/(2\Delta_0)$ とし、上端を無限へ延ばすと積分は $\sqrt{2\pi T/\Delta_0}\,e^{-\Delta_0/T}$ である。したがって
$$
\Delta(T)=\Delta_0-\sqrt{2\pi\Delta_0 T}\,e^{-\Delta_0/T}
  [1+O(T/\Delta_0)]
$$
となる。最後の誤差の形は $\omega_D/\Delta_0$ を十分大きくした弱結合極限である。有限の切断では左辺の傾きに $\omega_D/\sqrt{\omega_D^2+\Delta_0^2}$ が掛かる。指数的な補正が $\delta/T$ でも小さくなるため、右辺に $\Delta$ ではなく $\Delta_0$ を入れる操作は自己整合する。

:::example[低温のギャップを固定してよい精度]
$T=\Delta_0/10$ とすると、この先頭の相対補正は $\sqrt{\pi/5}\,e^{-10}\simeq3.6\times10^{-5}$ である。ギャップの差を比較する測定でも、この値より小さい精度を求めるときは指数的な温度依存を保持する。
:::

## 超伝導状態の熱力学

### 自由エネルギーと比熱差

独立な準粒子をトレースすると
$$
\Omega_s=V|\Delta|^2/g+\sum_k(\xi_k-E_k)-2T\sum_k\log(1+e^{-E_k/T}).
\label{eq-mb-bcs-omega}
$$
$\Delta=0$ を引けば正常相との差になり、熱力学極限で $\partial_{\Delta^*}\Omega_s=0$ は :eqref[eq-mb-bcs-gap] と一致する。零温度では単位体積あたりの差は $\Delta^2/g-2\nu_0\int_0^{\omega_D}(\sqrt{\xi^2+\Delta^2}-\xi)\dd\xi$ である。原始関数 $\frac12[\xi\sqrt{\xi^2+\Delta^2}+\Delta^2\operatorname{arsinh}(\xi/\Delta)-\xi^2]$ を使い、ギャップ方程式で対数を消すと $-\nu_0\omega_D(\sqrt{\omega_D^2+\Delta^2}-\omega_D)\simeq-\nu_0\Delta_0^2/2$ となる。

$T_c$ 近傍で $\delta\Omega/V=a|\Delta|^2+b|\Delta|^4/2$ と書く。式 :eqref[eq-mb-pair-sum] からギャップ積分に和を入れ、$\Delta^2$ で微分して $b=\nu_0 T_c\sum_n\int_\R\dd\xi/(\omega_n^2+\xi^2)^2=7\zeta(3)\nu_0/(8\pi^2T_c^2)$ を得る。積分は $\pi/(2|\omega_n|^3)$ であり、正の奇数の三乗逆数の和が $7\zeta(3)/8$ である。二次係数は正常ギャップ積分の差から $a\simeq\nu_0(T-T_c)/T_c$ となる。最小値 $-a^2/(2b)$ を温度で二度微分すると $\Delta C/V=T_c(a')^2/b$ であり、正常相の $C_n/V=2\pi^2\nu_0T_c/3$ で割ると $\Delta C/C_n=12/[7\zeta(3)]\simeq1.43$ になる。低温の準粒子分布は $e^{-\Delta/T}$ を持ち、ノードのないギャップでは励起の比熱は指数的に減る。

:::example[凝縮エネルギーの尺度]
$\nu_0=2\,\mathrm{eV^{-1}}$ を一単位胞あたりとし $\Delta_0=1\,\mathrm{meV}$ を入れると、差は約 $-1\,\mathrm{\mu eV}$ 毎単位胞である。比熱の跳びが有限でも、凝縮エネルギーは電子の帯域幅より小さい。
:::

:::example[転移直下のギャップ]
自由エネルギーを最小にすると $\Delta(T)^2\simeq8\pi^2T_c(T_c-T)/[7\zeta(3)]$ である。例えば $(T_c-T)/T_c=0.01$ では $\Delta/T_c\simeq0.306$ となる。この近似は転移直下の極限であり、零温度の値の計算には使わない。
:::

### 零温度の凝縮エネルギーを端点まで積分する

単位体積の片スピン状態密度 $\nu_0$ を一定とし、両側のカットオフを $W$ とする。零温度で正常状態との差を単位体積で表すと $\Delta^2/g+2\nu_0\int_0^W(\xi-\sqrt{\xi^2+\Delta^2})\dd\xi$ である。$\Delta>0$ のとき平方根の原始関数は $[\xi\sqrt{\xi^2+\Delta^2}+\Delta^2\operatorname{arsinh}(\xi/\Delta)]/2$ である。ギャップ方程式 $1/g=\nu_0\operatorname{arsinh}(W/\Delta)$ を用いると、逆双曲線関数が相殺され、
$$
\frac{\Omega_s-\Omega_n}V=\nu_0 W(W-\sqrt{W^2+\Delta^2})
=-\frac{\nu_0 W\Delta^2}{W+\sqrt{W^2+\Delta^2}}
$$
を得る。弱結合では $-\nu_0\Delta^2/2$ となる。積分を先に展開して $\sqrt{\xi^2+\Delta^2}\simeq\xi+\Delta^2/(2\xi)$ とすると $\xi=0$ で発散するので、閉じた式で相殺を行ってから展開する。

:::example[有限帯域での差]
$W=3\,\mathrm{meV}$、$\Delta=4\,\mathrm{meV}$ なら差は $-6\nu_0\,\mathrm{meV^2}$ となる。弱結合式の $-8\nu_0\,\mathrm{meV^2}$ と比べると、$\Delta/W$ が小さくないときに端点を無限に延ばす誤差が確認できる。
:::

### 低温比熱の前因子

三次元の等方的なギャップが全 Fermi 面上にあり、片スピンの正常状態密度 $\nu_0$ が滑らかとする。$T\ll\Delta_0$ では $f(E)\simeq e^{-E/T}$ であり、準粒子のエントロピーは一軌道あたり $(1+E/T)e^{-E/T}$ が先頭である。これは $-f\log f-(1-f)\log(1-f)$ を小さい $f$ で展開して得られる。$E=\Delta_0+\xi^2/(2\Delta_0)+O(\xi^4/\Delta_0^3)$ とすると、熱で効く殻は $|\xi|\sim\sqrt{\Delta_0T}$ である。二スピンの $\xi$ 積分を Gauss 積分で行い、$E/T$ の先頭に $\Delta_0/T$ を入れれば
$$
\frac{S_s}V\simeq2\nu_0 \sqrt{2\pi\Delta_0 T}\,
\frac{\Delta_0}{T}e^{-\Delta_0/T},\qquad
\frac{C_s}V\simeq 2\nu_0\sqrt{2\pi\Delta_0 T}\,
\frac{\Delta_0^2}{T^2}e^{-\Delta_0/T}.
$$
$C=T\partial_T S$ の最大の二次因子は指数の微分から来る。有限温度の $\Delta(T)$ のずれは、この範囲ではギャップ方程式に現れる熱占有と同じ指数で抑えられている。これを既に指数で小さいエントロピーへ戻す効果は高次である。励起の縮退や節がある場合は小さいエネルギーでの状態数を数え直す。

### Zeeman エネルギーと常伝導への転移

軌道効果を無視し、上と下のスピンをエネルギー $-h,+h$ だけ移す。均一な一重項では Bogoliubov の二枝は $E_k-h$ と $E_k+h$ になる。$T=0$、$h<\Delta_0$ の範囲ではまだ熱力学的な占有がなく、一重項の自由エネルギーはこの近似で変化しない。正常相の磁化の変分は $\delta n_\uparrow-\delta n_\downarrow=2\nu_0 h$ であり、$-\partial\Omega_n/\partial h$ に等しい。積分して $\Omega_n(h)/V=\Omega_n(0)/V-\nu_0 h^2$ である。従って零磁場の弱結合の凝縮エネルギー $-\nu_0\Delta_0^2/2$ と比較すると
$$
\frac{\Omega_s(h)-\Omega_n(h)}V\simeq-\frac{\nu_0\Delta_0^2}2+\nu_0h^2
$$
となり、二つの一様相の自由エネルギーは $h=\Delta_0/\sqrt2$ で交差する。この :term[Pauli 限界]{en="Pauli paramagnetic limit" yomi="ぱうりげんかい"} は、均一な対が粒子の枝が閉じる $h=\Delta_0$ より前に自由エネルギーで負けることを表す。有限運動量の対を許す極低温の清浄系や、強いスピン–軌道相互作用がある場合の相境界はこの比較だけでは決まらない。

:::example[磁場で二枝がまだ開いている場合]
$g$ 因子を $2$ とし、$h=\mu_B B$、$\Delta_0=1\,\mathrm{meV}$ とすると一様相の交差は約 $12.2\,\mathrm T$ である。そのとき小さい枝に残る最小エネルギーは約 $0.293\,\mathrm{meV}$ である。磁場をかけて得た自由エネルギーの交差は、スペクトルが先に零になることを前提にしない。
:::

### 自由エネルギーの勾配から安定な枝を選ぶ

自由エネルギー差を $x=|\Delta|^2\ge0$ の関数として書くと、
$$
\frac1V\frac{\partial(\Omega_s-\Omega_n)}{\partial x}
=\frac1g-\nu_0\int_0^{\omega_D}
\frac{\tanh(\sqrt{\xi^2+x}/(2T))}{\sqrt{\xi^2+x}}\dd\xi .
$$
右側の積分は $x$ で厳密に減るので、この導関数は増加する。$T<T_c$ では $x=0$ で負、$x\to\infty$ では正となるため、非零のギャップの枝は大域的な最小である。$T\ge T_c$ では端 $x=0$ が最小である。位相の回転による縮退は $x$ の最小化には含めていない。固定した一様な殻と一定の状態密度を使うこの模型では、ギャップの方程式に非零の根が見つかったとき、それが自由エネルギーの最大であるという別の解釈はこの単調性で排除できる。スピンを分裂させたり複数のバンドの自由度を許したりすると、同じ単一変数の単調性は前提から導き直す必要がある。

### スピン磁化率の低温則

一様な弱い Zeeman 場 $h$ に対し、上の準粒子は $E-h$、下は $E+h$ である。スピン密度の差は対の生成・消滅と可換であり、一重項の真空は磁化を持たない。従って $\delta m=n_\uparrow-n_\downarrow$ の変化は $V^{-1}\sum_k[f(E_k-h)-f(E_k+h)]\simeq2h\nu_0\int_\R[-f'(E)]\dd\xi$ となる。片スピンの状態密度はここでも $\nu_0$ であり、正常の磁化率は $2\mu_m^2\nu_0$ である。比は
$$
Y(T)=\int_\R[-f'(\sqrt{\xi^2+\Delta^2})]\dd\xi.
$$
これを :term[Yosida 関数]{en="Yosida function" yomi="よしだかんすう"} と呼ぶ。$\Delta=0$ なら $-f'(|\xi|)$ は $-f'(\xi)$ と同じ偶関数であり、全積分は $1$ である。低温では $-f'(E)\simeq e^{-E/T}/T$ を用い、殻の Gauss 積分から $Y(T)\simeq\sqrt{2\pi\Delta_0/T}\,e^{-\Delta_0/T}$ となる。これはスピン回転を保つ相のスピン応答であり、軌道磁性やスピン–軌道結合のある試料の総磁化には別の寄与も入る。

:::example[温度で測る磁化率の比]
$T=\Delta_0/10$ では $Y(T)\simeq\sqrt{20\pi}e^{-10}\simeq3.60\times10^{-4}$ である。正常のパウリ磁化率を基準にした比である。低温で粒子のトンネルのギャップとスピン応答の両方に同じ $\Delta_0$ が現れるが、前因子は測る演算子により異なる。
:::

## Nambu–Gor'kov Green 関数

### 行列逆と異常相関

有限温度の :term[Nambu–Gor'kov Green 関数]{en="Nambu–Gor'kov Green function" yomi="なんぶごるこふぐりーんかんすう"} は $\mathcal G_N(\tau)=-\langle T_\tau \Psi_k(\tau)\Psi_k^\dagger(0)\rangle$ である。運動方程式または固有ベクトル分解で
$$
\mathcal G_N(k,z)=(zI-h_k)^{-1}
=\frac{1}{z^2-E_k^2}\begin{pmatrix}z+\xi_k&-\Delta\\-\Delta^*&z-\xi_k\end{pmatrix}.
\label{eq-mb-nambu-green}
$$
$11$ 成分の留数は $+E_k$ で $u^2$、$-E_k$ で $|v|^2$ である。$12$ 成分は粒子数を二つ変える異常振幅を持つ。$z=i\omega_n$ で成分を和に取り :eqref[eq-mb-pair-sum] を適用すれば $\langle b_k\rangle$ が得られる。固定粒子数の有限系では異常平均は零だが、対の二点相関は有限であり得る。ここでは対の位相を選んだ平均場を用いる。

:::example[正常極限]
$\xi>0$ で $\Delta\to0$ とすると $u^2\to1$ なので、$11$ 成分の負エネルギー極は消える。異常成分も零になり、$11$ 成分は $1/(z-\xi)$ に戻る。
:::

### トンネル分光で現れる状態密度

正常側の一スピン状態密度を $\nu_0$ とし、一定で両側に十分な帯域があると仮定する。正のエネルギー $\omega$ のトンネルで見えるスペクトルの重みは $\nu_0\int \dd\xi\,u_\xi^2\delta(\omega-\sqrt{\xi^2+\Delta^2})$ である。$\omega>\Delta$ なら $\xi=\pm\sqrt{\omega^2-\Delta^2}$ の二根があり、同じ Jacobian $\omega/\sqrt{\omega^2-\Delta^2}$ を持つ。両根の $u_\xi^2$ を足すと $1$ である。負のエネルギーでは同じ計算で $|v_\xi|^2$ が入り、両根の和はやはり $1$ である。従って $\nu_s(\omega)/\nu_0=|\omega|/\sqrt{\omega^2-\Delta^2}$ が帯域内のギャップの外に得られ、内側の状態密度は零である。温度を上げると測定の微分伝導度はこの関数と $-\partial_\omega f$ の畳み込みになる。

:::example[端の重みと高いエネルギー]
$\Delta=4\,\mathrm{meV}$、$\omega=5\,\mathrm{meV}$ なら $\nu_s/\nu_0=5/3$ となる。$\omega$ がギャップより十分大きい領域ではこの比は $1$ に近づく。ギャップ直上の発散は積分可能であり、元の $\xi$ の区間を $\omega$ の区間に移す Jacobian がその由来である。
:::

### 異常伝播の虚時間和と対振幅

式 :eqref[eq-mb-nambu-green] の非対角成分は $-\Delta/(z^2-E^2)$ である。式 :eqref[eq-mb-pair-sum] を使えば $\beta^{-1}\sum_n\mathcal G_{12}(i\omega_n)=\Delta\tanh(\beta E/2)/(2E)$ となる。この和は $O(\omega_n^{-2})$ なので等時刻左極限を取れる。$12$ 成分の左極限では時間順序の交換により $\mathcal G_{12}(0^-)=\langle c_{-k\downarrow}c_{k\uparrow}\rangle$ である。この対振幅を $\Delta=(g/V)\sum_k\langle b_k\rangle$ へ代入するとギャップ方程式に戻る。

:::example[対振幅の温度依存]
$E=5\,\mathrm{meV}$、$\Delta=4\,\mathrm{meV}$ を固定して計算する。$T=0$ で振幅は $2/5$、$\beta E=\log3$ では $1/5$ である。実際の自己無撞着計算では $\Delta$ 自体も変わるため、この比は同じ BdG 行列に対する占有の効果を示す。
:::

### 密度頂点と対頂点の違い

$\Psi^\dagger\tau_3\Psi=n_{k\uparrow}+n_{-k\downarrow}-1$ なので、電荷の密度頂点は Nambu 空間で $\tau_3$ である。$\tau_i$ は Pauli 行列で、実の $\Delta$ では $h_k=\xi_k\tau_3-\Delta\tau_1$ である。行列要素を正負の固有ベクトルで挟むと、正の励起の対角密度は $u^2-v^2=\xi/E$ である。他方、正負の枝をつなぐ密度頂点の行列要素は $2uv=\Delta/E$ である。従って分光の単一粒子の重み $u^2$ と、密度の遷移に現れる振幅は異なる。一般の波数では例えば $\langle +,\vb k'|\tau_3|+,\vb k\rangle=u_{k'}u_k-v_{k'}v_k$ となり、ギャップの位相や軌道に応じた符号の干渉がある。

一様な密度に対し固定した平均場だけを用いると、枝間の $\Delta/E$ は周波数 $2E$ に見かけの吸収を作る。しかし元の数を保存する Hamiltonian では、厳密な全密度の交換子は零なので、非零周波数の一様密度応答があってはならない。対ポテンシャルの位相が外場に応じて動く頂点を含めると、この矛盾は解消する。すなわち密度の頂点を平均場の行列と整合させる必要がある。次節でその項を Ward 恒等式に入れる。

:::example[粒子–正孔混合の中心]
$\xi=0$ では $u^2=v^2=1/2$、通常の密度の正枝への対角行列要素は零、正負枝間の行列要素の絶対値は $1$ である。この二つは同じ行列の異なる成分であり、正枝の電子付加スペクトル重みの $1/2$ とも区別される。
:::

### 正常金属とのトンネル電流をスペクトルから求める

金属の状態密度と界面の一電子の行列要素が、ギャップと測定する電圧の窓で一定とする。電圧の規約を $\mu_N-\mu_S=eV$ とし、超伝導側の周波数をその化学ポテンシャルから測る。一次の遷移の黄金律で、左から右の占有因子 $f(\omega-eV)[1-f(\omega)]$ から逆の因子を引くと $f(\omega-eV)-f(\omega)$ が残る。粒子の伝播のスペクトルを用いて $\nu_s(\omega)$ が重みになる。従って正常状態の伝導度 $G_N$ による規格化で
$$
I(V)=\frac{G_N}{e}\int_\R
 \frac{\nu_s(\omega)}{\nu_0}[f(\omega-eV)-f(\omega)]\,\dd\omega .
$$
$eV$ の微分により $G(V)/G_N$ は $[-f'(\omega-eV)]\nu_s(\omega)/\nu_0$ の積分である。零温度で $V>0$ のとき積分の窓は $[0,eV]$ である。$eV>\Delta$ なら $\omega/\sqrt{\omega^2-\Delta^2}$ の原始関数から $I(V)=(G_N/e)\sqrt{(eV)^2-\Delta^2}$ であり、それより下では零である。背景の状態密度と粒子–正孔の対称性を保つ場合、$V<0$ 側は電流の符号だけ反転する。透過が大きく対をまとめて運ぶ寄与がある場合には、ここで残した一電子の最低次だけでは不十分である。

:::example[状態密度の発散と電流の連続性]
$eV=\Delta+\epsilon$、$0<\epsilon\ll\Delta$ とすると電流は $(G_N/e)\sqrt{2\Delta\epsilon}$ に近づく。微分伝導度は $\sqrt{\Delta/(2\epsilon)}$ の桁で大きくなるが、電流自身は閾値で連続である。熱で幅がある場合には状態密度を微分した値ではなく、状態密度と熱窓の畳み込みを使う。
:::

## 電磁応答と Meissner 効果

### 位相と磁場の応答

粒子電荷を $-e$ とし、対の位相を $\theta$ とする。連続な静的系の位相のエネルギー密度は $n_s(\hbar\nabla\theta+2e\vb A)^2/(8m)$ である。対の速度が $\vb v_s=(\hbar\nabla\theta+2e\vb A)/(2m)$ で、単位体積の運動費用が $mn_sv_s^2/2$ であることから係数が得られる。清浄で Galilei 不変の零温度では全体の一様ブーストがこの運動費用を与えるため $n_s=n$ である。有限温度では励起エネルギーに $\vb k\cdot\vb v_s$ が加わり、熱自由エネルギーの二次微分が $2\int\dd^3k/(2\pi)^3\, k^2f'(E)/(3m)$ を $n_s$ に足す。ここで $\hbar=1$ の運動量を使い、$f'<0$ なので熱励起が剛性を減らす。

横方向の $\vb A$ に対して位相を一定に取ると $\vb j=-n_se^2\vb A/m$ になる。SI の静的 Maxwell 方程式 $\nabla\times\vb B=\mu_0\vb j$ から $\nabla^2\vb B=\vb B/\lambda_L^2$、$\lambda_L=\sqrt{m/(\mu_0n_se^2)}$ を得る。半無限体では表面から $e^{-x/\lambda_L}$ で減衰する。この :term[Meissner 効果]{en="Meissner effect" yomi="まいすなーこうか"} は静的横応答であり、正常金属の一様な動的 Drude 重みとは極限の順序が異なる。縦応答で純ゲージの $\vb A$ を消すには同時に位相を変分する。Nambu の粒子数頂点にこの位相応答が入ることで、連続の式と整合する。

:::example[侵入長]
$n_s=10^{28}\,\mathrm{m^{-3}}$、電子質量なら $\lambda_L\simeq53\,\mathrm{nm}$ である。剛性が四分の一になると侵入長は二倍になる。
:::

:::example[磁場の減衰を測る]
$\lambda_L=50\,\mathrm{nm}$ の試料で表面から $100\,\mathrm{nm}$ の点では、磁場は表面値の $e^{-2}\simeq0.135$ になる。磁気エネルギーの密度は磁場の二乗に比例するので $e^{-4}\simeq0.0183$ 倍になる。
:::

### 位相剛性から侵入長を求める

電子の電荷を $-e$、対の位相を $\theta$ とする。対の質量を $2m$ とすれば流速は $\vb v_s=(\hbar\grad\theta+2e\vb A)/(2m)$ であり、電子数で測った超流動密度を $n_s$ とすれば運動エネルギー密度は $n_sm|\vb v_s|^2/2$ である。$\vb A$ の微分の負を取ると $\vb j_s=-en_s\vb v_s$ となる。渦のない横応答では位相の勾配の横成分は零なので $\vb j_s=-e^2n_s\vb A/m$ である。静磁場で $\rot\vb B=\mu_0\vb j_s$、$\divg\vb B=0$ を使ってもう一度回転を取ると $\nabla^2\vb B=\vb B/\lambda_L^2$、$\lambda_L=\sqrt{m/(\mu_0 e^2 n_s)}$ となる。半空間 $x>0$ に一定方向の場を入れ、奥で発散しない解を取れば $B(x)=B(0)e^{-x/\lambda_L}$ である。

:::example[超流動密度の低下]
$n_s$ が四分の一になると侵入長は二倍になる。例えば $\lambda_L=50\,\mathrm{nm}$ の状態で深さ $100\,\mathrm{nm}$ の比は $e^{-2}$ であり、密度低下後は $e^{-1}$ である。磁場が静的で横方向のため、この薄れ方は電場に対する直流抵抗から直接計算できない。
:::

### Nambu Ward 恒等式と位相の応答

外場が無いときの逆行列を $G_N^{-1}(K)=i\omega_nI-\xi_k\tau_3+\Delta\tau_1$ とする。粒子数を回転する生成子が $\tau_3$ であることから、電荷を除いた完全頂点は :eqref[eq-mb-ward] を行列に拡張した
$$
q_0\Gamma^0-\vb q\cdot\vb\Gamma
=G_N^{-1}(K+Q)\tau_3-\tau_3 G_N^{-1}(K)
$$
を満たす。実の $\Delta$ の固定した平均場に右辺を評価すると、$q_0=0,\vb q=0$ のための極限に $[\Delta\tau_1,\tau_3]=-2i\Delta\tau_2$ が残る。裸頂点の左辺は $Q$ に比例して消えるので、このままでは一致しない。対の位相を回すと非対角の $\tau_2$ 成分が一次で現れる。従って位相を外場で変分した頂点が必要になる。以下の位相の作用は、低周波数で電荷を保存する頂点を構成する一つの方法である。

電子電荷 $-e$ で $\vb A\mapsto\vb A+\nabla\Lambda$、$\Phi\mapsto\Phi-\partial_t\Lambda$ と取ると、位相は $\theta\mapsto\theta-2e\Lambda/\hbar$ と変わる。従って空間の組み合わせ $\hbar\nabla\theta+2e\vb A$ と時間の組み合わせ $\hbar\partial_t\theta-2e\Phi$ は不変である。密度の変位 $\delta n$ を位相と共役とし、低速でのエネルギーの二次費用を $\delta n^2/(2\chi_n)$、$\chi_n=\partial n/\partial\mu>0$ と置く。位相の自由度だけの時間項は、Lagrangian の
$$
-\frac12\delta n(\hbar\partial_t\theta-2e\Phi)
-\frac{\delta n^2}{2\chi_n}
$$
を $\delta n$ について平方完成して得る。結果は $\chi_n(\hbar\partial_t\theta-2e\Phi)^2/8$ である。空間項は自由エネルギーから負符号で入り、$-n_s(\hbar\nabla\theta+2e\vb A)^2/(8m)$ である。$e=0$ なら変分方程式は速度 $\sqrt{n_s/(m\chi_n)}$ の波動方程式である。三次元の自由状態密度を用いた弱結合、零温度の中性系なら $n_s=n$、$\chi_n=2\nu_0$ なので速度は $v_F/\sqrt3$ になる。温度や周波数によって励起の減衰が無視できない場合には、局所な二次作用だけでは不十分である。

### Coulomb 力と静的横応答

縦密度の運動を求めるときは、密度の二次費用の逆圧縮率へ Coulomb 核 $v_C(q)=e^2/(\epsilon_{\rm vac}q^2)$ を加える。長波長でエネルギーが二次で、背景誘電率を真空値とした場合である。位相と密度の二つの変分式を用いると、縦波は
$$
\omega^2=\frac{n_s}{m}\,[\chi_n^{-1}+v_C(q)]q^2
=\frac{n_s e^2}{m\epsilon_{\rm vac}}+\frac{n_s}{m\chi_n}q^2
$$
を満たす。純粋な位相作用で得る係数は、位相と密度だけの低周波数の模型の係数である。通常の金属でプラズマ周波数が $2\Delta/\hbar$ より大きい場合、実際の高周波数の分散をこの低周波数の式で外挿しない。清浄な Galilei 不変の零温度で $q=0$ の重みを固定するのは全粒子数と質量の総和則である。一方、静的横ベクトルポテンシャルには位相の勾配が結合しない。従って静的な磁場の侵入長は同じ $n_s/m$ を含むが、縦の集団モードと同じ極限を測るものではない。

:::example[侵入長とプラズマの尺度の比較]
上の理想化した局所の式では $\omega_p\lambda_L=(\mu_0\epsilon_{\rm vac})^{-1/2}=c$ である。$\lambda_L=50\,\mathrm{nm}$ とすると対応する $\hbar\omega_p$ の桁は $3.95\,\mathrm{eV}$ である。ギャップが $\mathrm{meV}$ の系に対しては低周波数の位相作用をその周波数まで使えないことが、この桁の比較から分かる。
:::

### 対波動関数の空間の大きさ

一様で弱結合の金属の対振幅は運動量表示で $\kappa_k=\Delta/(2\sqrt{\xi_k^2+\Delta^2})$ である。硬い殻の不連続は位置の長距離に人為的な裾を作るので、位置の二次モーメントを計算する際は殻の端を滑らかにする。$\omega_D\gg\Delta$ で低エネルギーの部分が支配的と仮定し、片スピン状態密度と $v_F$ を定数とする。Parseval の関係より対の平均二乗距離は $\int_{\vb k}|\nabla_{\vb k}\kappa_k|^2/\int_{\vb k}|\kappa_k|^2$ である。$\hbar$ を戻し、$\nabla_k\xi=\hbar\vb v_F$ を用いると
$$
\langle r^2\rangle \simeq\hbar^2 v_F^2
\frac{\int_\R \xi^2(\xi^2+\Delta^2)^{-3}\dd\xi}
{\int_\R(\xi^2+\Delta^2)^{-1}\dd\xi}
=\frac{\hbar^2 v_F^2}{8\Delta^2}.
$$
分子は $\pi/(8\Delta^3)$、分母は $\pi/\Delta$ であり、$\xi=\Delta\tan\varphi$ で評価できる。磁場応答の非局所性や温度による GL 回復長を定義する際の長さと、この二次モーメントにはそれぞれ定義がある。いずれも弱結合の低温の領域では微視的な距離の尺度として $\hbar v_F/\Delta$ が現れる。

:::example[対が多数の粒子間隔にまたがる場合]
$v_F=10^6\,\mathrm{m\,s^{-1}}$、$\Delta=1\,\mathrm{meV}$ では、この定義の対の大きさは約 $233\,\mathrm{nm}$ である。粒子間隔が $0.5\,\mathrm{nm}$ なら数百倍である。弱結合で多くの対が重なり、固い二粒子分子の気体として数えられないことと整合する。
:::

### 静的な横電流の微視的な小波数展開

清浄で等方的な三次元の放物線帯、$T=0$、$\Delta>0$ の一重項を取る。磁場と電流が横方向で $\vb q$ を $z$ 軸、ベクトルポテンシャルを $x$ 軸に取る。$\hbar=1$ とする。粒子と正孔ではベクトルポテンシャルで運動量の引数が逆にずれ、さらに正孔のエネルギーの符号が逆である。この二つが重なるため、線形の電流の頂点は $e v_x I$ の大きさになる。すでに見たスカラー密度の $\tau_3$ とは異なる。帯の二次の微分には反磁性の $ne^2/m$ が入り、横応答では位相を固定してよい。

$\xi_\pm=\xi_{\vb k\pm\vb q/2}$、$E_\pm=\sqrt{\xi_\pm^2+\Delta^2}$ とする。枝の射影子を $P_\pm=(I\pm h/E)/2$ と書き、トレースを取ると異なる符号の枝をつなぐ二つの組の重みの和は $1-(\xi_+\xi_-+\Delta^2)/(E_+E_-)$ となる。零温度でその遷移の分母は $E_++E_-$ であり、応答によるエネルギーの二次変化は仮想励起により負である。従って $\delta j_x=-K_T(q) A_x$ の核は
$$
K_T(q)=\frac{ne^2}{m}
-e^2\int_{\vb k} v_x^2
\frac{1-(\xi_+\xi_-+\Delta^2)/(E_+E_-)}{E_++E_-}.
$$
$\int_{\vb k}$ は片スピンの軌道を走る。二成分の組の遷移を二つ数えたため、ここでさらにスピンの係数を掛けない。$q=0$ で分子は零となり、応答は反磁性の値である。正常の静的横極限では Fermi 面を横切る粒子–正孔励起があるため、$q=0$ をこのギャップのある式へ先に代入した結果を使えない。

$\Delta\ll E_F$、$v_Fq\ll\Delta$ で殻を線形化し、$\delta=v_z q/2$ とする。単位ベクトル $(\xi,-\Delta)/E$ の $\xi$ 微分の二乗は $\Delta^2/E^4$ であり、$\xi\pm\delta$ での二つの単位ベクトルの内積は $1-2\delta^2\Delta^2/E^4+O(\delta^4)$ になる。従って分数は $\delta^2\Delta^2/E^5$ が先頭である。$E_++E_-=2E+O(q^2)$ を使った。$\xi$ 積分 $\int_\R\Delta^2(\xi^2+\Delta^2)^{-5/2}\dd\xi=4/(3\Delta^2)$ と、球面の $\langle v_x^2v_z^2\rangle=v_F^4/15$ を用いる。後者は $\langle\hat k_z^4\rangle=1/5$ と $\hat k_x^2+\hat k_y^2=1-\hat k_z^2$ から出る。従って
$$
\frac{K_T(q)}{ne^2/m}=1-\frac{v_F^2q^2}{30\Delta^2}+\cdots .
$$
$n/m=2\nu_0 v_F^2/3$ を用いた。SI では $v_Fq$ に $\hbar$ を掛ける。この式は長波長の滑らかな磁場に対する局所 London 核の最初の補正である。磁場が $\hbar v_F/\Delta$ と同じ尺度で変化する場合は二次だけでは不足し、非局所の全核や表面の境界条件を調べる。

:::example[局所近似の小波数の誤差]
$v_F q=\Delta/2$ とし、弱結合の殻の式をそのまま評価すれば、$q^2$ 項は局所値の $1/120$ である。$\Delta$ が四分の一になると同じ $q$ で $v_F q/\Delta=2$ となる。その場合、二次の小波数補正だけで誤差を見積もる前提が崩れる。
:::

## 演習問題

::::exercise{level=basic}
$\Delta\ne0$ とし、$\xi=0$ の対の $u^2$ と $|v|^2$ を求めよ。

:::solution
$u$ の式に $\xi/E=0$ を代入すると、双方とも $1/2$ になる。対の占有の平均は $1/2$、電子数の平均は $1$ である。
:::
::::

::::exercise{level=basic}
$\xi=3$、$\Delta=4$ の準粒子エネルギーを求めよ。

:::solution
$h^2=(\xi^2+\Delta^2)I$ であり、正の固有値は $\sqrt{9+16}=5$ である。
:::
::::

::::exercise{level=basic}
$g\nu_0=2/\log 3$ の二粒子問題で $B/\omega_D$ を求めよ。

:::solution
式 :eqref[eq-mb-cooper-binding] から $B/\omega_D=2/(3-1)=1$ である。
:::
::::

::::exercise{level=basic}
同じ $m,e$ の超伝導体で $n_s$ が四倍になったとき、侵入長の倍率を求めよ。

:::solution
$\lambda_L\propto n_s^{-1/2}$ なので $1/2$ 倍である。
:::
::::

::::exercise{level=advanced}
実の $\Delta$ で $h^2=E^2 I$ と Nambu Green 関数を示せ。

:::solution
非対角の積は $-\xi\Delta+\xi\Delta=0$ であり、対角の積は $\xi^2+\Delta^2$ である。従って $(zI-h)(zI+h)=(z^2-E^2)I$ であり、逆は $(zI+h)/(z^2-E^2)$ である。
:::
::::

::::exercise{level=advanced}
Nambu 関数の $11$ 成分を部分分数にして全スペクトル重みを求めよ。

:::solution
$z=E$ の留数は $(E+\xi)/(2E)$、$z=-E$ では $(-E+\xi)/(-2E)$ である。和は $1$ になり、正準反交換関係の総和則を再現する。
:::
::::

::::exercise{level=advanced}
零温度のギャップ方程式の積分を評価し、弱結合極限の指数を導け。

:::solution
$\xi=\Delta\sinh u$ と置けば $\dd\xi/E=\dd u$ なので、積分は $\operatorname{arsinh}(\omega_D/\Delta)$ である。それを $1/(g\nu_0)$ に等置し、$\sinh y\simeq e^y/2$ を使うと $2\omega_D e^{-1/(g\nu_0)}$ になる。
:::
::::

::::exercise{level=advanced}
自由エネルギー差が $a'(T-T_c)x+bx^2/2$、$x=|\Delta|^2$、$b>0$ のとき、$a'>0$ として比熱の跳びを導け。

:::solution
$T<T_c$ の停留点は $x=a'(T_c-T)/b$ であり、値は $-a'^2(T-T_c)^2/(2b)$ となる。$T>T_c$ では $x=0$ が最小で差は零である。$-T\partial_T^2$ の差を $T_c$ に近づければ $T_c a'^2/b$ になる。
:::
::::

::::exercise{level=advanced}
弱結合の零温度で $h=\Delta_0/2$ の軌道効果を無視した Zeeman 場を加える。均一な一重項と常伝導相の大ポテンシャル差と、準粒子の小さい方のギャップを求めよ。

:::solution
小さい枝は $E-h$ なので最小は $\Delta_0/2>0$ である。そのため対の状態ではこの範囲で粒子の占有を変えない。一方、正常相は磁化により単位体積あたり $\nu_0h^2$ だけ大ポテンシャルが下がる。従って差は $-\nu_0\Delta_0^2/2+\nu_0\Delta_0^2/4=-\nu_0\Delta_0^2/4$ で、この二候補の比較では均一な一重項が低い。
:::
::::

::::exercise{level=advanced}
一電子トンネルの零温度の式で $eV=2\Delta>0$ とする。$I/(G_N\Delta/e)$ と $G(V)/G_N$ を求めよ。

:::solution
電流は積分の上端までの $\sqrt{(eV)^2-\Delta^2}$ から $I/(G_N\Delta/e)=\sqrt3$ となる。上端で微分すると、$G(V)/G_N=(eV)/\sqrt{(eV)^2-\Delta^2}=2/\sqrt3$ である。いずれも同じ正常状態の伝導度で規格化した。
:::
::::

## 確認テスト

::::quiz
引力が弱くても三次元の Cooper 問題で束縛解が現れる理由を選べ。

- [x] Fermi 面に近い対状態の積分が対数発散する。
- [ ] 有効質量が必ず無限大になる。
- [ ] 状態密度が必ず無限大になる。

:::explanation
積分 $\int_0^{\omega_D}\dd\xi/(2\xi+B)$ は一定の有限状態密度でも $B\to0$ で発散する。そのため任意の正の $g$ で束縛条件を満たせる。
:::
::::

::::quiz{answer="5" tol="1e-6"}
$\xi=3$、$\Delta=4$ の準粒子の正のエネルギーを答えよ。

:::explanation
二つの成分の二乗和の平方根で $\sqrt{3^2+4^2}=5$ となる。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
$\Delta\ne0$、$\xi=0$ における粒子成分の正の極の重みを答えよ。

:::explanation
$11$ 成分の留数は $(1+\xi/E)/2$ であり、$\xi=0$ なら $1/2$ である。
:::
::::

::::quiz
均一な完全ギャップの弱結合 BCS 相について正しいものを選べ。

- [x] 準粒子の熱励起は低温で指数的に減る。
- [x] 静的横方向の磁場は侵入長で減衰する。
- [ ] 単一粒子ギャップと Cooper 二体束縛エネルギーが等しい。

:::explanation
最低の準粒子エネルギーが $\Delta_0>0$ なので Boltzmann 因子が現れる。横応答と Maxwell 方程式から磁場の指数減衰を得る。一方、二体と多体では弱引力で指数に入る係数がそれぞれ $2$ と $1$ である。
:::
::::

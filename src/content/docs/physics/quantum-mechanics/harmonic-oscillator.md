---
title: 調和振動子
description: 放物線型の井戸を交換関係と微分方程式の両方で解き、古典軌道を追う最小波束を構成する。
notation:
  - symbol: 'a^\dagger'
    meaning: 生成演算子
  - symbol: '\ell'
    meaning: 振動子の長さ
  - symbol: '\alpha'
    meaning: コヒーレント状態の複素振幅
---

## この章の目標

- Hamilton 演算子を因数分解して準位と規格化を決める。
- Hermite 多項式を用いて座標波動関数を求める。
- コヒーレント状態の確率分布、幅、時間発展を計算する。

## 生成・消滅演算子

安定平衡点 $x_0$ の周りで $V(x)=V(x_0)+V''(x_0)(x-x_0)^2/2+\cdots$ である。
変位が小さく三次以上を無視できれば、原点とエネルギーの基準を取り直して $H=P^2/(2m)+m\omega^2X^2/2$ となる。
ここで $\omega=\sqrt{V''(x_0)/m}>0$ は角振動数である。

長さ $\ell=\sqrt{\hbar/(m\omega)}$ によって位置と運動量を無次元化し、:term[消滅演算子]{en="annihilation operator" yomi="しょうめつえんざんし"} とその随伴である :term[生成演算子]{en="creation operator" yomi="せいせいえんざんし"} を

$$
\begin{aligned}
a&=\frac1{\sqrt2}\left(\frac X\ell+i\frac{\ell P}{\hbar}\right), \\
a^\dagger&=\frac1{\sqrt2}\left(\frac X\ell-i\frac{\ell P}{\hbar}\right)
\end{aligned}
$$

と置く。
$[X,P]=i\hbar I$ より $[a,a^\dagger]=I$。
積 $a^\dagger a$ は古典的な平方の和に $i[X,P]/(2\hbar)=-1/2$ を加えたものである。
従って:term[個数演算子]{en="number operator" yomi="こすうえんざんし"} $N=a^\dagger a$ で

$$
\begin{aligned}
H&=\hbar\omega(N+1/2), \\
X&=\ell(a+a^\dagger)/\sqrt2, \\
P&=\hbar(a-a^\dagger)/(i\sqrt2\ell).\label{eq-qm-oscillator}
\end{aligned}
$$

これらの積は例えば急減少する滑らかな関数の共通領域で計算する。
$\braket\psi{N\psi}=\norm{a\psi}^2\ge0$ なので下限は $\hbar\omega/2$ である。

:::example[赤外振動の尺度]
$\nu=\omega/(2\pi)=30.0\,\mathrm{THz}$ とすれば $\hbar\omega=h\nu=0.124\,\mathrm{eV}$、:term[零点エネルギー]{en="zero-point energy" yomi="れいてんえねるぎー"}は $0.0620\,\mathrm{eV}$ である。
角振動数に $h$ を掛けると $2\pi$ 倍のずれが生じる。
:::

### 不確定性から見た底のエネルギー

古典的には $X=P=0$ とすれば底で静止できる。
量子力学では中心を原点に置いても二つの幅が残る。
$\Delta X=s$ と固定すれば、運動エネルギーは $\expval{P^2}/(2m)\ge\hbar^2/(8ms^2)$、ポテンシャルの期待値は少なくとも $m\omega^2s^2/2$ である。
従って

$$
\expval H\ge\frac{\hbar^2}{8ms^2}+\frac{m\omega^2s^2}{2}\ge\frac{\hbar\omega}{2}.
$$

最後は相加相乗平均であり、等号の条件は $s^2=\ell^2/2$ である。最初の不等式でも等号にするには、位置と運動量の平均がともに零で、その幅の積が $\hbar/2$ でなければならない。座標表示の基底 Gaussian がこれらを同時に満たす。
運動量と位置が釣り合うこの長さが、因数分解に現れた $\ell$ の意味である。

:::example[一定力でずれた井戸]
$H=P^2/(2m)+m\omega^2X^2/2-FX$ なら、平方完成で平衡点は $X_0=F/(m\omega^2)$。
$X$ を $X-X_0$ に替えて同じ因数分解ができ、全準位が $-F^2/(2m\omega^2)$ だけ下がるが間隔は等しいままである。
:::

### 平衡点を決めてから振動数を定める

二原子分子の伸縮では原子核の位置をそれぞれ $x_1,x_2$ とし、ポテンシャルを $V(x_1-x_2)$ と近似する。重心 $X=(m_1x_1+m_2x_2)/(m_1+m_2)$、相対変位 $r=x_1-x_2-r_0$ に変えると、運動エネルギーは $P_X^2/[2(m_1+m_2)]+p_r^2/(2\mu)$、$\mu=m_1m_2/(m_1+m_2)$ となる。平衡位置で $V'(r_0)=0$ だから、二次の係数 $k=V''(r_0)>0$ から $\omega=\sqrt{k/\mu}$ が決まる。従って一方の原子を同位体に替えたときには、電子状態に由来する $k$ を固定できる近似の下で、換算質量を変更する必要がある。重い方だけを固定端とみなして軽い原子の質量を直接入れると、両者が同程度の質量の分子で周波数を過小評価する。

### 因数分解は定数だけずれた二つの問題を結ぶ

無次元座標 $\xi=x/\ell$ で $a=(\xi+\partial_\xi)/\sqrt2$、$a^\dagger=(\xi-\partial_\xi)/\sqrt2$ とする。関数 $f$ にまとまりを保って作用させると
$$
2a^\dagger a f=\xi^2 f+\xi f'-(\xi f)'-f''=(-\partial_\xi^2+\xi^2-1)f.
$$
順序を逆にすると最後の定数が $+1$ に変わる。底状態の $a u_0=0$ は $u_0'=-\xi u_0$ なので $u_0\propto e^{-\xi^2/2}$。一方 $a^\dagger f=0$ の解 $e^{+\xi^2/2}$ は規格化できない。最低状態は存在しても最高状態はないという非対称性が、二つの一次方程式の遠方での振る舞いに現れる。


### 同位体置換した分子から結合の曲率を測る

伸縮の基本線からポテンシャルの底の曲率を推定する。二つの核はともに運動し、電子状態は置換の前後でほぼ同じとして調和近似を適用する。

相対座標の質量は $m_1m_2/(m_1+m_2)$。振動の角周波数が曲率と換算質量の比の平方根であることから通常周波数へ換算する。

**順方向の予測。** 以下では丸める前に同じ定数を使い、最後に表示する値だけを丸める。

| 入力 | 記号 | 設定値 |
|---|---|---|
| 軽い原子の質量 | $m_1$ | $1.6735\times10^{-27}\,\mathrm{kg}$ |
| 重い原子の質量 | $m_2$ | $5.81\times10^{-26}\,\mathrm{kg}$ |
| 平衡点での曲率 | $K$ | $480\,\mathrm{N\,m^{-1}}$ |



$$
\nu=\frac{\sqrt{K} \sqrt{m_1 + m_2}}{2 \pi \sqrt{m_1} \sqrt{m_2}}\simeq 8.6455807007\times10^{13}\,\mathrm{Hz}.
$$

**測定値からの逆算。** 次に平衡点での曲率だけを未知とする。振動の通常周波数を測り、同じ関係を未知量について解く。ここで物理的な範囲にある解を選ぶ。

$$
K=\frac{4 \pi^{2} m_1 m_2 \nu^{2}}{m_1 + m_2}\simeq 480\,\mathrm{N\,m^{-1}}.
$$

この二つの計算は同じ理想関係を逆にたどる。実際の記録では、観測値や別に校正した入力にも誤差がある。各入力の相対標準偏差を下表のように仮定し、誤差は小さく互いに独立とする。対数を微分した係数 $c_i=(u_i/v)(\partial v/\partial u_i)$ は、微小な相対変化の比を符号つきで表す。近い測定値の差が分母にある場合には、その差に比べても小さな変化に限ってこの線形近似を使う。

| 揺れる量 | 対数微分 | 係数の値 | 入力の幅（%） | 推定値への寄与（%） |
|---|---|---|---|---|
| 振動の通常周波数 | $2$ | $2$ | $0.2$ | $0.4$ |
| 軽い原子の質量 | $\frac{m_2}{m_1 + m_2}$ | $0.972$ | $0.01$ | $0.00972$ |
| 重い原子の質量 | $\frac{m_1}{m_1 + m_2}$ | $0.027997$ | $0.01$ | $2.7997\times10^{-4}$ |

分散を加算してから平方根を取ると、独立な誤差について

$$
\frac{\sigma_v}{|v|}\simeq\sqrt{\sum_i c_i^2\frac{\sigma_i^2}{u_i^2}}\simeq 0.0040012,\qquad \sigma_v\simeq 1.9206\,\mathrm{N\,m^{-1}}.
$$

軽い核をほぼ二倍の質量の同位体にすると、曲率を保っても吸収線が移動する。

| 軽い原子の質量 | 振動の通常周波数 |
|---|---|
| $1.6735\times10^{-27}\,\mathrm{kg}$ | $8.6456\times10^{13}\,\mathrm{Hz}$ |
| $3.3445\times10^{-27}\,\mathrm{kg}$ | $6.2005\times10^{13}\,\mathrm{Hz}$ |

置換で井戸の非調和性や平均核間距離も目に見えて変わる場合、この推定は各同位体で調和近似した有効曲率として解釈する。


### Gaussian 幅が最適でなくても底状態を近似する

$x=X/\ell$ とする。正確な基底波形を知らず、やや幅の狭い Gaussian $e^{-x^2}$ から計算を始めた場合を考える。偶の多項式を順に掛け、波形を広げる自由度を足す。

無次元座標 $x$ で積分し、エネルギーの一単位を $b=1.9878\times10^{-20}\,\mathrm J$ とする。試行関数の線形結合 $\psi=\sum_j c_j v_j$ の係数をすべて自由にする。今回の候補は

$$
(v_0,\ldots,v_{2})=\left(e^{- x^{2}},\quad x^{2} e^{- x^{2}},\quad x^{4} e^{- x^{2}}\right),\quad -\infty<x<\infty.
$$

である。分母の規格化と分子の二つのエネルギーを別々に積分する。実関数の基底を選んだので、以下の行列は実対称である。

$$
G=\begin{pmatrix}\frac{\sqrt{2} \sqrt{\pi}}{2} & \frac{\sqrt{2} \sqrt{\pi}}{8} & \frac{3 \sqrt{2} \sqrt{\pi}}{32}\\\frac{\sqrt{2} \sqrt{\pi}}{8} & \frac{3 \sqrt{2} \sqrt{\pi}}{32} & \frac{15 \sqrt{2} \sqrt{\pi}}{128}\\\frac{3 \sqrt{2} \sqrt{\pi}}{32} & \frac{15 \sqrt{2} \sqrt{\pi}}{128} & \frac{105 \sqrt{2} \sqrt{\pi}}{512}\end{pmatrix},\qquad K=\begin{pmatrix}\frac{\sqrt{2} \sqrt{\pi}}{4} & - \frac{\sqrt{2} \sqrt{\pi}}{16} & - \frac{9 \sqrt{2} \sqrt{\pi}}{64}\\- \frac{\sqrt{2} \sqrt{\pi}}{16} & \frac{7 \sqrt{2} \sqrt{\pi}}{64} & \frac{21 \sqrt{2} \sqrt{\pi}}{256}\\- \frac{9 \sqrt{2} \sqrt{\pi}}{64} & \frac{21 \sqrt{2} \sqrt{\pi}}{256} & \frac{225 \sqrt{2} \sqrt{\pi}}{1024}\end{pmatrix},\qquad V=\begin{pmatrix}\frac{\sqrt{2} \sqrt{\pi}}{16} & \frac{3 \sqrt{2} \sqrt{\pi}}{64} & \frac{15 \sqrt{2} \sqrt{\pi}}{256}\\\frac{3 \sqrt{2} \sqrt{\pi}}{64} & \frac{15 \sqrt{2} \sqrt{\pi}}{256} & \frac{105 \sqrt{2} \sqrt{\pi}}{1024}\\\frac{15 \sqrt{2} \sqrt{\pi}}{256} & \frac{105 \sqrt{2} \sqrt{\pi}}{1024} & \frac{945 \sqrt{2} \sqrt{\pi}}{4096}\end{pmatrix}.
$$

例えば左上の成分は、順に次の三つの積分である。ここで運動エネルギーは部分積分後の二次形式を用い、境界条件を満たす候補を選んでいる。

$$
\begin{aligned}G_{00}&=\int_{-\infty}^{\infty}e^{- 2 x^{2}}\,\dd x=\frac{\sqrt{2} \sqrt{\pi}}{2},\\K_{00}&=\int_{-\infty}^{\infty}2 x^{2} e^{- 2 x^{2}}\,\dd x=\frac{\sqrt{2} \sqrt{\pi}}{4},\\V_{00}&=\int_{-\infty}^{\infty}\frac{x^{2} e^{- 2 x^{2}}}{2}\,\dd x=\frac{\sqrt{2} \sqrt{\pi}}{16}.\end{aligned}
$$

エネルギーを $b\epsilon$ と書き、$c^\dagger(H-\epsilon G)c$ を $c_j^*$ で微分すると $(H-\epsilon G)c=0$。$H=K+V$ である。候補の数を一つずつ増やしたときに行列式がどう変わるかを追う。

| 候補の数 | $\det(H-zG)$ | 最低根 $z$ | 最低変分値 |
|---|---|---|---|
| 1 | $- \frac{\sqrt{2} \sqrt{\pi} \left(8 z - 5\right)}{16}$ | $0.625$ | $1.2424\times10^{-20}\,\mathrm J$ |
| 2 | $\frac{\pi \left(64 z^{2} - 240 z + 107\right)}{1024}$ | $0.51715$ | $1.028\times10^{-20}\,\mathrm J$ |
| 3 | $- \frac{3 \sqrt{2} \pi^{\frac{3}{2}} \left(512 z^{3} - 4800 z^{2} + 10792 z - 4275\right)}{524288}$ | $0.50236$ | $9.9859\times10^{-21}\,\mathrm J$ |

最低根を選ぶ理由は、$G$ が正定値なので $G^{-1/2}HG^{-1/2}$ の通常のエルミート固有値問題に直せるからである。試行空間に新しい関数を加えた後も、その係数を零にすれば前の状態を再現できる。従って表の最低値は増加しない。正確な値は $\hbar\omega/2$ である。この表では幅を固定しても多項式の係数を動かすことで、狭すぎる裾が徐々に補われる。


## エネルギー固有値と固有状態

:eqref[eq-qm-oscillator] に現れる $N$ に対し、交換子の積の公式から $[N,a]=-a$、$[N,a^\dagger]=a^\dagger$。
従って $N\ket\lambda=\lambda\ket\lambda$ なら、$a\ket\lambda$ は消えるか固有値 $\lambda-1$ をもつ。
規格化した固有状態で $\norm{a\ket\lambda}^2=\lambda$、$\norm{a^\dagger\ket\lambda}^2=\lambda+1$ である。
$\lambda\ge0$ である。整数でない場合は $\lfloor\lambda\rfloor$ 回下げた後も固有値が $(0,1)$ にあるため、さらに一度下げたベクトルは非零で負の固有値を持つ。これは $N$ の非負性に反する。この議論は梯子演算を繰り返せる共通領域で行う。後の座標解で構成する完全正規直交系により、ほかのスペクトルがないことを確認する。
そこで底では $a\ket0=0$、上へは任意回進める。

$$
\begin{aligned}
\ket n&=\frac{(a^\dagger)^n}{\sqrt{n!}}\ket0, \\
a\ket n&=\sqrt n\ket{n-1}, \\
a^\dagger\ket n&=\sqrt{n+1}\ket{n+1}.
\end{aligned}
$$

異なる $n$ はエルミート演算子の異なる固有値で直交する。
エネルギーは $E_n=\hbar\omega(n+1/2)$ であり、絶対零度でも振動のエネルギーが残る。
底状態の存在と全空間に対する完全性は座標解で確認する。

:::example[位置の遷移行列要素]
$\mel m X n=\ell(\sqrt n\,\delta_{m,n-1}+\sqrt{n+1}\,\delta_{m,n+1})/\sqrt2$。
例えば $\mel2 X1=\ell$ だが $\mel3 X1=0$ である。
直線的な外力が直接結合する準位が読み取れる。
:::

:::example[運動と位置の揺らぎ]
$\ket n$ では $\expval X=\expval P=0$、$a^2,(a^\dagger)^2$ の対角要素も零である。
従って $\expval{X^2}=\ell^2(n+1/2)$、$\expval{P^2}=(\hbar/\ell)^2(n+1/2)$。
運動とポテンシャルに平均エネルギーが半分ずつ入る。
:::

### 時間相関と一量子の吸収

Hamilton 演算子との交換子は $[H,a]=-\hbar\omega a$。
従って Heisenberg 描像で $a(t)=a(0)e^{-i\omega t}$ である。
底状態の平均位置は零だが、二時刻の積の平均は

$$
\mel0{X(t)X(0)}0=\frac{\ell^2}{2}\mel0{(ae^{-i\omega t}+a^\dagger e^{i\omega t})(a+a^\dagger)}0=\frac{\ell^2}{2}e^{-i\omega t}.
$$

$aa^\dagger$ の項だけが残る。
順を逆にすれば複素共役となる。
従って時刻の異なる位置の交換子は $-i\ell^2\sin\omega t$。
基底状態でも位置のゆらぎに特徴的な周波数があることが分かる。

:::example[二次のプローブ]
$X^2$ を掛けたときは準位番号が $0,\pm2$ だけ変わりうる。
特に $X^2\ket0=(\ell^2/2)(\ket0+\sqrt2\ket2)$ なので $\mel2{X^2}0=\ell^2/\sqrt2$。
線形・二次の異なる結合を測れば、同じ井戸でも別の準位差が現れる。
:::

### 有限行列の上端で交換関係が崩れる

数値計算で $\ket0,\ldots,\ket{M-1}$ だけを残す。下げる操作は空間内で閉じるが、最高状態を上げようとすると次の状態が欠けている。射影した $a_M,a_M^\dagger$ の交換子を各基底に作用させると
$$
[a_M,a_M^\dagger]=I-M\ketbra{M-1}{M-1}.
$$
実際、下の $M-1$ 本には $(n+1)-n=1$ が掛かる。最高状態では第一項が零、第二項が $M-1$ なので $-(M-1)$ が掛かる。これで交換子の跡が零になることとも一致する。従って $a$ と $a^\dagger$ の積を使う計算では、一段途中で切れた経路が最終結果に戻って来ないかを確かめる。例えば $\langle M-1|X^2|M-1\rangle$ で、最初に $\ket M$ へ上がってから戻る寄与まで切ってしまうと誤る。演算子の積を作ってから射影する操作と、射影した演算子を掛ける操作は同じとは限らない。

### 二つのモードのエネルギーと縮退

独立な二方向の振動数が $\omega_x=\omega_y=\omega$ なら $H=\hbar\omega(N_x+N_y+1)$。全励起数 $N$ に対して $(n_x,n_y)=(0,N),(1,N-1),\ldots,(N,0)$ の $N+1$ 状態がある。三方向なら最初の $n_x$ を固定し残り二方向を数えて $\sum_{k=0}^N(k+1)=(N+1)(N+2)/2$ である。この縮退は各一次元問題の非縮退性と矛盾しない。振動数を異方的にして比を無理数にすれば、$\omega_x(n_x-n'_x)+\omega_y(n_y-n'_y)=0$ を満たす異なる整数対がなくなり、二次元でも縮退は消える。有理数比なら異なる方向の励起を交換して同じエネルギーを保てる場合が残る。


### 二次の位置結合で一段の遷移が消える理由

分極率の位置依存性が二次から始まる場合には $X^2$ の行列要素が必要になる。第二励起状態からの分岐を、途中の状態を省かずに求める。

$n=2$ を選び、質量 $1.66\times10^{-27}\,\mathrm{kg}$、通常周波数 $3.0\times10^{13}\,\mathrm{Hz}$ とする。

初めの規格化したケットを $\left(1\right)\ket{2}$ とする。以下の表では共通な次元を最後にまとめ、途中は無次元の係数で追う。

**最初の位置演算。** 無次元化した $a+a^\dagger$ は下に $\sqrt n$、上に $\sqrt{n+1}$ の係数を与える。


| 始状態 | 入ってくる係数 | 終状態 | この段の行列要素 | 終状態への寄与 |
|---|---|---|---|---|
| $2$ | $1$ | $1$ | $\sqrt{2}$ | $\sqrt{2}$ |
| $2$ | $1$ | $3$ | $\sqrt{3}$ | $\sqrt{3}$ |

同じ終状態への行を振幅の段階で合算すると、この段の後の列は

$$
\ket{w_1}=\left(\sqrt{2}\right)\ket{1}+\left(\sqrt{3}\right)\ket{3}.
$$



**二つ目の位置演算。** 二段目で元の準位に戻る二通りの経路は同じ終状態に達するため、符号を保って合算する。


| 始状態 | 入ってくる係数 | 終状態 | この段の行列要素 | 終状態への寄与 |
|---|---|---|---|---|
| $1$ | $\sqrt{2}$ | $0$ | $1$ | $\sqrt{2}$ |
| $1$ | $\sqrt{2}$ | $2$ | $\sqrt{2}$ | $2$ |
| $3$ | $\sqrt{3}$ | $2$ | $\sqrt{3}$ | $3$ |
| $3$ | $\sqrt{3}$ | $4$ | $2$ | $2 \sqrt{3}$ |

同じ終状態への行を振幅の段階で合算すると、この段の後の列は

$$
\ket{w_2}=\left(\sqrt{2}\right)\ket{0}+\left(5\right)\ket{2}+\left(2 \sqrt{3}\right)\ket{4}.
$$

となる。ここまでの積を $D$ と呼び、実際の演算子は $bD$、$b=1.6851\times10^{-22}\,\mathrm{m^2}$ とする。元のケットへの重なりと、準備から到達できる全終状態への強度の和は異なる。

$$
\mel\psi D\psi=5,\qquad \mel\psi{D^\dagger D}\psi=39.
$$

対角の重なりは $\langle X^2\rangle$ の計算に使い、二乗和は完全性により $\langle X^4\rangle$ に対応する。 行列要素の二乗が全強度に占める割合を知るには、全強度が非零なら各成分の絶対値二乗を全強度で割る。

| 最終量子数 | 絶対値二乗（$b^2$ 単位） | 全強度に占める割合 |
|---|---|---|
| $0$ | $2$ | $\frac{2}{39}$ |
| $2$ | $25$ | $\frac{25}{39}$ |
| $4$ | $12$ | $\frac{4}{13}$ |

位置を二回掛けるとパリティが二回反転する。最終表に偶数の量子数だけが残るのはこの対称性の結果である。


### 分子の赤外振動が室温で励起される割合

光の電磁モードと同じエネルギー間隔が分子の伸縮にも現れる。分子の回転や並進を除き、一つの調和モードが熱平衡にあると仮定して励起の内訳を調べる。

熱浴と平衡にある独立な振動モードで、$E_n=\hbar\omega(n+1/2)$ とする。$q=e^{-\hbar\omega/(k_{\rm B}T)}$ と書けば隣り合う Boltzmann 重みの比は $q$。従って正規化と一次・二次の和は

$$
\begin{aligned}p_n&=(1-q)q^n,\qquad\sum_{n=0}^\infty p_n=1,\\ \bar n&=(1-q)q\frac{\dd}{\dd q}(1-q)^{-1}=\frac q{1-q},\\ \langle n(n-1)\rangle&=(1-q)q^2\frac{\dd^2}{\dd q^2}(1-q)^{-1}=\frac{2q^2}{(1-q)^2}.\end{aligned}
$$

このため $\operatorname{var}(n)=\bar n(1+\bar n)$ である。零点エネルギーは分布の幅に寄与しないが、エネルギーの平均には加わる。周波数 $\nu=3\times10^{13}\,\mathrm{Hz}$、温度 $T=300\,\mathrm K$ に対し、$h\nu/(k_{\rm B}T)=4.7992$、隣接確率の比は $0.008236$ となる。

| 量子数 $n$ | 確率 $p_n$ | 準位エネルギー $[\mathrm J]$ | 平均エネルギーへの寄与 $[\mathrm J]$ |
|---|---|---|---|
| 0 | $0.99176$ | $9.9391\times10^{-21}$ | $9.8572\times10^{-21}$ |
| 1 | $0.0081681$ | $2.9817\times10^{-20}$ | $2.4355\times10^{-22}$ |
| 2 | $6.7273\times10^{-5}$ | $4.9696\times10^{-20}$ | $3.3432\times10^{-24}$ |
| 3 | $5.5406\times10^{-7}$ | $6.9574\times10^{-20}$ | $3.8548\times10^{-26}$ |
| 4 | $4.5632\times10^{-9}$ | $8.9452\times10^{-20}$ | $4.0819\times10^{-28}$ |

表に載せていない高い準位の確率の和は幾何級数の残り $q^5=3.7894\times10^{-11}$ である。無限和を保った最終値は

$$
\bar E=\hbar\omega(\bar n+1/2)\simeq 1.0104\times10^{-20}\,\mathrm J,\qquad \Delta E=\hbar\omega\sqrt{\bar n(1+\bar n)}\simeq 1.819\times10^{-21}\,\mathrm J.
$$

$k_{\rm B}T$ と $h\nu$ の大小を比べずに振動一自由度に古典の等分配を適用すれば、振動励起に貯まるエネルギーを大きく見積もる。


### トラップの励起側帯と脱励起側帯から温度を決める

閉じ込めた単一モードに弱い線形の位置結合を加え、内部準位の遷移に伴って振動量子を一個加える線と一個引く線を別々に測る。二本の検出効率と内部行列要素を同じにそろえ、熱平衡を仮定する。



| 設定する量 | 記号 | 計算に用いる値 |
|---|---|---|
| トラップの振動数 | $\nu$ | $2.2\times10^{6}\,\mathrm{Hz}$ |
| 一量子を減らす線の計数 | $N_-$ | $240$ |
| 一量子を増やす線の計数 | $N_+$ | $1000$ |

**二線の面積比。** 生成では行列要素二乗が n+1、消滅では n なので熱平均の比は平均励起数とそれに一を加えた数の比になる。


$$
\begin{aligned}q&=\frac{N_-}{N_+},\\ q&\simeq \frac{240}{1000}\simeq 0.24.\end{aligned}
$$

**平均振動量子数。** 測った比を n の平均について解く。


$$
\begin{aligned}\bar n&=\frac{q}{1 - q},\\ \bar n&\simeq \frac{0.24}{1 - 0.24}\simeq 0.31579.\end{aligned}
$$

**振動温度。** 熱平衡なら隣り合う占有の比が q である。対数を取ると温度を求められる。


$$
\begin{aligned}T&=- \frac{h \nu}{k_{\rm B} \ln{\left(q \right)}},\\ \frac{T}{\mathrm{K}}&\simeq - \frac{(6.6260701\times10^{-34}) (2.2\times10^{6})}{(1.380649\times10^{-23}) \ln{\left(0.24 \right)}}\simeq 7.3984\times10^{-5}.\end{aligned}
$$

**基底状態の割合。** 幾何分布の規格化条件から、最低準位の占有は一から隣接比を引いた値である。


$$
\begin{aligned}p_0&=1 - q,\\ p_0&\simeq 1 - 0.24\simeq 0.76.\end{aligned}
$$

**零点を含む平均エネルギー。** 励起数の平均に二分の一を加え、一量子のエネルギーを掛ける。


$$
\begin{aligned}U&=h \nu \left(\bar n + \frac{1}{2}\right),\\ \frac{U}{\mathrm{J}}&\simeq (6.6260701\times10^{-34}) (2.2\times10^{6}) \left(0.31578947 + \frac{1}{2}\right)\simeq 1.1892\times10^{-27}.\end{aligned}
$$

増加側の計数を一定に取り、減少側の計数を変える。同じ面積比の変化でも高温側ほど推定温度に大きく効く。


| $N_-$ | $\bar n$ | $T\;[\mathrm{K}]$ | $p_0$ |
|---|---|---|---|
| $50$ | $0.052632$ | $3.5245\times10^{-5}$ | $0.95$ |
| $120$ | $0.13636$ | $4.9797\times10^{-5}$ | $0.88$ |
| $240$ | $0.31579$ | $7.3984\times10^{-5}$ | $0.76$ |
| $500$ | $1$ | $1.5232\times10^{-4}$ | $0.5$ |
| $800$ | $4$ | $4.7316\times10^{-4}$ | $0.2$ |

赤側帯が弱いのは最低準位からさらに下へ移れないためである。熱平衡が成立しなければ面積比から平均励起数はなお得られるが、基底占有率や温度までは一意に決まらない。


## 座標表示と Hermite 多項式

### 固有方程式と漸化式

$\xi=x/\ell$、$\epsilon=E/(\hbar\omega)$ とすると $u''+(2\epsilon-\xi^2)u=0$ である。
大きな $|\xi|$ で減衰する要因を $e^{-\xi^2/2}$ と見て $u=e^{-\xi^2/2}h(\xi)$ と置くと

$$
h''-2\xi h'+(2\epsilon-1)h=0,\qquad c_{j+2}=\frac{2j+1-2\epsilon}{(j+2)(j+1)}c_j
$$

となる。
$h=\sum_j c_j\xi^j$ の偶系列と奇系列は独立である。
$2\epsilon-1=2n$ のとき次数 $n$ の系列が打ち切れ、残りの独立解は両端で規格化できない（一次元束縛準位の非縮退性）。
逆に、打ち切れない級数が規格化の条件を満たさないことを微分方程式から確かめよう。
仮に $u$ が二乗可積分なら、$[u(\xi)\pm u(-\xi)]/2$ も二乗可積分なので、偶解と奇解を別々に調べれば十分である。
パリティを $s=0,1$ として $d_k=c_{2k+s}$ とおく。
上の漸化式から

$$
\frac{d_{k+1}}{d_k}
=\frac{4k+2s+1-2\epsilon}{(2k+s+2)(2k+s+1)},\qquad
\lim_{k\to\infty}(k+1)\frac{d_{k+1}}{d_k}=1
$$

である。無限に続く非零の系列では末尾の符号がそろう（初期係数の位相は実に選べる）。
$1/2<b<1$ を固定すれば、十分大きい $k$ 以降は $|d_{k+1}|\ge b|d_k|/(k+1)$ となる。
従って $\xi>0$ で末尾の絶対値を足すと、$\xi^s e^{b\xi^2}$ の正の定数倍から有限多項式を引いたもの以上になる。
最初の有限個の項はこの成長を相殺できない。
$u=e^{-\xi^2/2}h$ に戻しても $b>1/2$ だから無限遠で発散する。
よって打ち切り条件 $\epsilon=n+1/2$ は必要でもあり、代数的解法と独立に同じスペクトルが得られる。

座標表示で $a=(\xi+\partial_\xi)/\sqrt2$、$a^\dagger=(\xi-\partial_\xi)/\sqrt2$ である。
$a u_0=0$ を解いて Gaussian 積分を用いると $u_0=\pi^{-1/4}\ell^{-1/2}e^{-\xi^2/2}$。
昇る演算子を繰り返し用いて

$$
\begin{aligned}
u_n(x)&=\frac{H_n(\xi)e^{-\xi^2/2}}{\pi^{1/4}\sqrt{2^n n!\ell}}, \\
H_n(\xi)&=(-1)^ne^{\xi^2}\frac{\dd^n}{\dd\xi^n}e^{-\xi^2}.
\end{aligned}
$$

実際 $(\xi-\partial_\xi)f=-e^{\xi^2/2}\partial_\xi(e^{-\xi^2/2}f)$ を反復すればこの Rodrigues 公式を得る。
$e^{2\xi t-t^2}=e^{\xi^2}e^{-(\xi-t)^2}$ を $t$ で Taylor 展開すると係数が $H_n/n!$ となる。
両辺を $t,\xi$ で微分して比較すれば $H_{n+1}=2\xi H_n-2nH_{n-1}$、$H_n'=2nH_{n-1}$ である。

### 完全性

すべての $u_n$ と直交する $f\in L^2$ があれば、先頭係数が $2^n$ なので全単項式について $\int f(\ell\xi)e^{-\xi^2/2}\xi^n\dd\xi=0$ である。
$F(z)=\int f(\ell\xi)e^{-\xi^2/2}e^{z\xi}\dd\xi$ は整関数である。
実際 Cauchy–Schwarz により任意の有界な $z$ の範囲で Gaussian が全次数の微分を支配する。
原点での全導関数が零なので $F=0$。
虚軸上で Fourier 変換の一意性を使えば $f=0$ である。
従って得た列が状態空間全体を張る。有限個の $u_n$ の線形結合上で、微分演算子 $H$ は $E_n$ を掛ける対角演算子に対応する。その閉包の定義域は $\sum_n E_n^2|c_n|^2<\infty$ を満たす係数列である。実際、この条件を満たす列は部分和で状態と $H$ の作用をともにノルム近似できる。逆に両方が収束すれば成分ごとの極限からこの条件が従う。随伴の定義を各 $u_n$ に適用しても同じ条件と作用を得るため、この閉包は自己共役である。$E_n$ は有限の集積点を持たず、これらの値以外では $(E_n-z)^{-1}$ が有界なので、スペクトルは得られた固有値だけからなる。

:::example[節とパリティ]
$H_0=1,H_1=2\xi,H_2=4\xi^2-2$。
第二励起状態の節は $x=\pm\ell/\sqrt2$ であり、反転に対する符号は $(-1)^2=+1$ である。
多項式の次数と内部の節の数が対応する。
:::

### 密度を古典の転回点と比較する

基底エネルギー $E_0=\hbar\omega/2$ に対応する古典的転回点は $\pm\ell$。
規格化した $u_0$ を使うと、その外側にある確率は

$$
2\int_\ell^\infty |u_0(x)|^2\dd x=\frac{2}{\sqrt\pi}\int_1^\infty e^{-\xi^2}\dd\xi=\operatorname{erfc}(1)\simeq0.1573 .
$$

この領域でも局所密度があることと、そこへ向かう定常な流れがあることは別である。
$u_0$ を実に選べるため、底状態の流束は全域で零である。

:::example[規格化を積分で追う]
$H_2=4\xi^2-2$ を使うと、$\int e^{-\xi^2}H_2^2\dd\xi=16(3\sqrt\pi/4)-16(\sqrt\pi/2)+4\sqrt\pi=8\sqrt\pi$。
波動関数の規格化因子に入る $2^2 2!\sqrt\pi$ と一致する。
次数が高い場合は同じ結果を生成関数からまとめて得られる。
:::

実際、二つの生成関数の積を Gaussian 重みで積分すると、平方完成により $\int e^{-\xi^2}e^{2t\xi-t^2}e^{2s\xi-s^2}\dd\xi=\sqrt\pi e^{2st}$。
$s^mt^n$ の係数を左右で比較すれば $\int e^{-\xi^2}H_mH_n\dd\xi=\sqrt\pi\,2^n n!\delta_{mn}$ となる。

### 母関数から規格化と直交性を同時に得る

Hermite 多項式の母関数を $G(\xi,t)=\exp(2\xi t-t^2)=\sum_{n\ge0}H_n(\xi)t^n/n!$ とする。$e^{-\xi^2}G(\xi,t)=e^{-(\xi-t)^2}$ なので、$t$ の $n$ 階微分を零で取れば $(-1)^n\partial_\xi^n e^{-\xi^2}$ になる。これで Rodrigues の式と一致する。二つの変数 $s,t$ を用意して積を積分すると
$$
\int_{\mathbb R}e^{-\xi^2}G(\xi,s)G(\xi,t)\,\dd\xi
=\int_{\mathbb R}e^{-(\xi-s-t)^2+2st}\,\dd\xi
=\sqrt\pi\,e^{2st}.
$$
右辺には $s^n t^n$ だけが含まれる。左右で係数を比較すれば $\int e^{-\xi^2}H_m H_n\,\dd\xi=\sqrt\pi\,2^nn!\delta_{mn}$。展開と積分を交換できるのは、$s,t$ を有界な範囲に限ったとき、各有限階微分が Gaussian に一次指数と多項式を掛けた可積分関数で押さえられるためである。

### 運動量表示での節と位相

空間単位 $\ell$ に対応する運動量単位は $\hbar/\ell$ である。$q=p\ell/\hbar$ とおくと正規化した無次元 Fourier 変換 $\mathcal F f(q)=(2\pi)^{-1/2}\int e^{-iq\xi}f(\xi)\dd\xi$ が得られる。基底 Gaussian は同じ Gaussian に移る。$\mathcal F[\xi f]=i\partial_q\mathcal Ff$、$\mathcal F[\partial_\xi f]=iq\mathcal F f$ であるから、一回上げるごとに $\mathcal F a^\dagger=-i a^\dagger_q\mathcal F$。帰納的に $\mathcal F u_n=(-i)^n u_n$ となる。従って固有状態の位置密度と運動量密度は尺度を合わせれば同じ形である。一方、たとえば $\ket0+\ket1$ の重ね合わせでは相対位相が変わるので、この一致を各成分の和にそのまま移すことはできない。


### 位置分布の特性関数から高次モーメントを読む

同じ平均と標準偏差をもつ分布でも、干渉計へ入れたときの応答は同じとは限らない。空間周波数 $k$ の位相板を通した経路と基準経路との重なりは $\chi_n(k)=\langle n|e^{ikX}|n\rangle$ である。$x_0=\sqrt{\hbar/(2m\omega)}$、$\eta=kx_0$ と置き、座標表示で $\chi_n(k)$ を計算する。二つの Hermite 母関数の積に位相因子を加え、平方完成する。

$$
\int_{\R} e^{-\xi^2} G(\xi,s)G(\xi,t)e^{i\sqrt2\eta\xi}\,\dd\xi
=\sqrt{\pi}\,e^{-\eta^2/2}e^{2st}e^{i\sqrt2\eta s}e^{i\sqrt2\eta t}.
$$

左辺の $s^nt^n$ の係数に $n!/(2^n\sqrt\pi)$ を掛けると規格化された $\chi_n(k)$ が得られる。右辺で $e^{2st}$ から次数 $n-r$ を選ぶと、残りの二つの指数からはどちらも次数 $r$ を選ぶ。従って寄与する項は $0\le r\le n$ に限られる。

$$
\chi_n(k)=e^{-\eta^2/2}\sum_{r=0}^n\frac{n!(-\eta^2)^r}{(n-r)!(r!)^2}
=e^{-\eta^2/2}L_n(\eta^2).
$$

ここで最後の等号を Laguerre 多項式 $L_n$ の定義に使ってもよい。特性関数は偶関数なので奇数次のモーメントはすべて零。二次まで読むと $\langle X^2\rangle=x_0^2(2n+1)$ であり、四次まででは $\langle X^4\rangle=3x_0^4(2n^2+2n+1)$ となる。後者は生成消滅演算子を四回作用させて $n$ へ戻る経路の和と一致する。基底状態 $n=0$ だけでは $\langle X^4\rangle=3\langle X^2\rangle^2$、すなわち全位置分布が Gaussian であることを再現する。

$n=1$ の場合は $\chi_1=e^{-\eta^2/2}(1-\eta^2)$。位置密度は常に非負でも、異なる位置で付いた位相を平均するため $\eta=1$ で干渉が消え、それを越すと干渉縞の位相が $\pi$ だけ変わる。位相板を付けた側の粒子が吸収されたという解釈なら、二経路を単独で通した透過率まで下がるはずである。この例では各経路の透過率は一であり、相対位相の分布が可視度を減らしている。


### 微分方程式と強度表で Hermite 解の低準位を検査する

調和振動子の方程式から Gaussian を外した演算子は $D=-\partial_x^2+2x\partial_x$ であり、$D H_n=2nH_n$。母関数の係数を実際に展開して比較する。

ここで $x$ は無次元変数である。波動関数に共通な因子を外し、その絶対値二乗と座標の測度をまとめて $w(x)=e^{- x^{2}}$ とする。残りの多項式 $p_n$ に作用する演算を $D$ と書く。微分方程式に代入して、上の次数の項から順に係数を照合していく。

| 次数 | $p_n(x)$ | 固有値 $D p_n/p_n$ | $\int w p_n^2\,\dd x$ |
|---|---|---|---|
| 0 | $1$ | $0$ | $\sqrt{\pi}$ |
| 1 | $2 x$ | $2$ | $2 \sqrt{\pi}$ |
| 2 | $2 \left(2 x^{2} - 1\right)$ | $4$ | $8 \sqrt{\pi}$ |
| 3 | $4 x \left(2 x^{2} - 3\right)$ | $6$ | $48 \sqrt{\pi}$ |
| 4 | $4 \left(4 x^{4} - 12 x^{2} + 3\right)$ | $8$ | $384 \sqrt{\pi}$ |

表の最後の行は、多項式をそのまま微分して再現できる。各次数の係数を追うために、一階と二階の導関数を分けて書く。

$$
\begin{aligned}p_{4}^\prime(x)&=32 x \left(2 x^{2} - 3\right),\\p_{4}^{\prime\prime}(x)&=96 \left(2 x^{2} - 1\right),\\D p_{4}&=32 \left(4 x^{4} - 12 x^{2} + 3\right)=8\,p_{4}.\end{aligned}
$$

固有値が異なる解の直交性を、実際の積分でも確かめられる。例えば最初の行と最後の行の積に重みを掛け、各単項式ごとに積分する。

| 次数 $k$ | 積の係数 $a_k$ | $\int a_k x^k w(x)\,\dd x$ |
|---|---|---|
| 0 | $12$ | $12 \sqrt{\pi}$ |
| 2 | $-48$ | $- 24 \sqrt{\pi}$ |
| 4 | $16$ | $12 \sqrt{\pi}$ |

列の和は零になる。積分を打ち切った区間ではこの相殺が成立する保証はない。

さらに $n=1$ に座標 $x$ を掛ける。得られる状態の分布を求めるには、積の多項式をもう一度各規格化関数に射影する。

| 終状態 $m$ | $\langle m\vert x\vert n\rangle$ | 行列要素の二乗 | 終状態の物理量 |
|---|---|---|---|
| 0 | $\frac{\sqrt{2}}{2}$ | $\frac{1}{2}$ | $9.9391\times10^{-21}\,\mathrm J$ |
| 1 | $0$ | $0$ | $2.9817\times10^{-20}\,\mathrm J$ |
| 2 | $1$ | $1$ | $4.9696\times10^{-20}\,\mathrm J$ |
| 3 | $0$ | $0$ | $6.9574\times10^{-20}\,\mathrm J$ |
| 4 | $0$ | $0$ | $8.9452\times10^{-20}\,\mathrm J$ |

全ノルムは $\int w\,|x p_n|^2\,\dd x/\|p_n\|^2=\frac{3}{2}$、表の二乗和は $\frac{3}{2}$ である。一致することから、この入力に $x$ を一回掛けた列はすでに表の空間の中にある。ここでは周波数 $3.0\times10^{13}\,\mathrm{Hz}$ を用いた。


### 励起状態の位置を検出する窓

第一励起状態では中央に節がある。検出器を片側だけに置くと、左右対称な密度の条件付き平均を読み取ることになる。変数は $x=X/\ell$ とする。

この数値例の質量は $1.66\times10^{-27}\,\mathrm{kg}$、周波数は $3.0\times10^{13}\,\mathrm{Hz}$ である。

以下では無次元の読みを $x$、その一単位に対応する実際の目盛を $b=1.8358\times10^{-11}\,\mathrm m$ とする。まず全領域を積分したときに確率が一になるよう係数を決める。形だけを表す関数を $w$ とすると

$$
\begin{aligned}w(x)&=x^{2} e^{- x^{2}},\qquad -\infty<x<\infty,\\ Z&=\int_{-\infty}^{\infty}w(x)\,\dd x=\frac{\sqrt{\pi}}{2},\\ f(x)&=\frac{w(x)}{Z}=\frac{2 x^{2} e^{- x^{2}}}{\sqrt{\pi}}.\end{aligned}
$$

目盛を有限幅に区切ると、区間に入る確率は曲線下の面積で決まる。実際の幅に対する密度は $f(x)/b$ である。したがって幅を変換してから積分しても、無次元の区間の両端で原始関数を引いても同じ値を得る。ここで使う原始関数は

$$
F(x)=- \frac{\left(2 x - \sqrt{\pi} e^{x^{2}} \operatorname{erf}{\left(x \right)}\right) e^{- x^{2}}}{2 \sqrt{\pi}},\qquad F^\prime(x)=f(x).
$$

各区間に入る割合を求める。境界ちょうどにある点は連続分布の確率を持たないため、隣り合う区間に二重に数える心配はない。

| 無次元の区間 | $F(右端)-F(左端)$ | 確率の概数 |
|---|---|---|
| $-\infty\ \text{–}\ -1$ | $- \frac{- e \sqrt{\pi} - 2 + e \sqrt{\pi} \operatorname{erf}{\left(1 \right)}}{2 e \sqrt{\pi}}$ | 0.2862 |
| $-1\ \text{–}\ 0$ | $\frac{-2 + e \sqrt{\pi} \operatorname{erf}{\left(1 \right)}}{2 e \sqrt{\pi}}$ | 0.2138 |
| $0\ \text{–}\ 1$ | $\frac{-2 + e \sqrt{\pi} \operatorname{erf}{\left(1 \right)}}{2 e \sqrt{\pi}}$ | 0.2138 |
| $1\ \text{–}\ \infty$ | $- \frac{- e \sqrt{\pi} - 2 + e \sqrt{\pi} \operatorname{erf}{\left(1 \right)}}{2 e \sqrt{\pi}}$ | 0.2862 |

次に、このうち $-1<x<0$ だけを受ける検出器を置く。ゲートを通過した事象を分母にして平均すると

$$
\begin{aligned}p_{\rm gate}&=\frac{-2 + e \sqrt{\pi} \operatorname{erf}{\left(1 \right)}}{2 e \sqrt{\pi}},\\\int_{-1}^{0}x f(x)\,\dd x&=- \frac{-2 + e}{e \sqrt{\pi}},\\ \overline{x}_{\rm gate}&=- \frac{2 \left(-2 + e\right) \left(2 + e \sqrt{\pi} \operatorname{erf}{\left(1 \right)}\right)}{-4 + \pi e^{2} \operatorname{erf}^{2}{\left(1 \right)}}\simeq -0.69731.\end{aligned}
$$

となる。これを元の単位に戻した条件付き平均は $-1.2801\times10^{-11}\,\mathrm m$ である。一方、入射した集団全体では二つの積分

$$
\begin{aligned}\langle x\rangle&=\int_{-\infty}^{\infty}x f(x)\,\dd x=0,\\\langle x^2\rangle&=\int_{-\infty}^{\infty}x^2 f(x)\,\dd x=\frac{3}{2}\end{aligned}
$$

を使える。差 $\langle x^2\rangle-\langle x\rangle^2=\frac{3}{2}$ の平方根に $b$ を掛ければ、位置の標準偏差は $2.2484\times10^{-11}\,\mathrm m$ となる。

左右をすべて集計すれば平均は零だが、節の近くでの検出は抑えられる。平均位置から最も見つかりやすい位置を推定できない例である。


### 波束の重心が静止しても幅が振動する場合

偶数パリティの三つの準位だけを重ねる。位置そのものの平均は常に零だが、位置の二乗は二段離れた準位を結ぶ。

基底は $n=0,2,4$ の順とし、共通の零点エネルギーを引く。検出量 $A$ はこの三準位へ射影した $X^2/x_0^2$、$x_0=\sqrt{\hbar/(2m\omega)}$ である。この測定では外の準位を除くので、表示する分散も射影した装置に対する値である。


$$
H=\epsilon M,\quad M=\begin{pmatrix}0 & 0 & 0\\0 & 2 & 0\\0 & 0 & 4\end{pmatrix},\quad A=\begin{pmatrix}1 & \sqrt{2} & 0\\\sqrt{2} & 5 & 2 \sqrt{3}\\0 & 2 \sqrt{3} & 9\end{pmatrix},\quad \ket{\psi(0)}=\begin{pmatrix}\frac{\sqrt{3}}{3}\\\frac{\sqrt{3} i}{3}\\\frac{\sqrt{3}}{3}\end{pmatrix}.
$$

この式で $A$ は示した目盛で測った無次元の読みで、エネルギー単位は $\epsilon=6.6261\times10^{-28}\,\mathrm J$ である。$M$ の各固有値を $u$ として $w_u=P_u\ket{\psi(0)}$ と置く。すると時刻 $t$ の列は $w_u$ のそれぞれに $e^{-iu\epsilon t/\hbar}$ を掛けた和になる。射影の段階では各列を規格化し直さない。

| 無次元エネルギー $u$ | 成分 $w_u$ | このエネルギーを測る確率 |
|---|---|---|
| $0$ | $\begin{pmatrix}\frac{\sqrt{3}}{3}\\0\\0\end{pmatrix}$ | $\frac{1}{3}$ |
| $2$ | $\begin{pmatrix}0\\\frac{\sqrt{3} i}{3}\\0\end{pmatrix}$ | $\frac{1}{3}$ |
| $4$ | $\begin{pmatrix}0\\0\\\frac{\sqrt{3}}{3}\end{pmatrix}$ | $\frac{1}{3}$ |

検出量の平均と分散を順に求める。まず $A$ の平均について、エネルギーの異なる二列を結ぶ係数を $a_{uv}=w_u^\dagger A w_v$ とする。ケットの位相とブラの共役位相が逆向きになるので、角周波数にはエネルギー差だけが入る。

$$
\langle A\rangle_t=\sum_{u,v} a_{uv}e^{i(u-v)\epsilon t/\hbar},\qquad a_{vu}=a_{uv}^*.
$$



$$
(a_{uv})=\begin{pmatrix}\frac{1}{3} & \frac{\sqrt{2} i}{3} & 0\\- \frac{\sqrt{2} i}{3} & \frac{5}{3} & - \frac{2 \sqrt{3} i}{3}\\0 & \frac{2 \sqrt{3} i}{3} & 3\end{pmatrix}.
$$



| 正の差 $v-u$ | 通常の振動数 | 余弦の係数 | 正弦の係数 |
|---|---|---|---|
| $2$ | $2\times10^{6}\,\mathrm{Hz}$ | $0$ | $\frac{2 \left(- 2 \sqrt{3} + \sqrt{2}\right)}{3}$ |
| $4$ | $4\times10^{6}\,\mathrm{Hz}$ | $0$ | $0$ |

同じ差を与える準位対の寄与は振幅として足した。振動数が許されていても、行列要素や初期位相によって係数が消える場合がある。$s=\epsilon t/\hbar$ と書くと今回の全信号は

$$
\langle A\rangle_t=\frac{2 \left(- 2 \sqrt{3} + \sqrt{2}\right) \sin{\left(2 s \right)}}{3} + 5,\qquad \overline{\langle A\rangle}=5.
$$

上線は観測窓を無限に長くした時間平均である。異なるエネルギー間の指数関数を積分して時間で割ると零に近づく。一方、同じエネルギー内の重ね合わせは射影 $P_u$ の内部に残る。全ての測定結果を捨ててエネルギーだけ先に測った状態は $\rho_{\rm d}=\sum_u w_u w_u^\dagger$ なので、この状況での一定の平均も同じ値になる。

$$
\rho_{\rm d}=\begin{pmatrix}\frac{1}{3} & 0 & 0\\0 & \frac{1}{3} & 0\\0 & 0 & \frac{1}{3}\end{pmatrix},\qquad \tr(\rho_{\rm d}A)=5.
$$

個別の信号に常にこの定数が見えるかは別に検査する。各時刻で分散も求めるには、同じ列に $A^2$ を作用させる。平均が零でも単発の読みが揃うとは限らない。

| $s$ | 経過時間 $[\mathrm s]$ | 平均 | 分散 |
|---|---|---|---|
| $0$ | $0$ | $5$ | $\frac{4 \left(\sqrt{6} + 15\right)}{3}$ |
| $\frac{\pi}{8}$ | $6.25\times10^{-8}$ | $\frac{17-2\sqrt{6}}{3}$ | $\frac{128-16\sqrt{6}}{9}$ |
| $\frac{\pi}{4}$ | $1.25\times10^{-7}$ | $\frac{- 4 \sqrt{3} + 2 \sqrt{2} + 15}{3}$ | $- \frac{4 \left(-31 - \sqrt{6} + 6 \sqrt{2} + 12 \sqrt{3}\right)}{9}$ |
| $\frac{\pi}{2}$ | $2.5\times10^{-7}$ | $5$ | $\frac{4 \left(\sqrt{6} + 15\right)}{3}$ |
| $\pi$ | $5\times10^{-7}$ | $5$ | $\frac{4 \left(\sqrt{6} + 15\right)}{3}$ |

周波数が二倍になるのは、位置二乗の非対角要素が二段離れた準位にだけ現れるためである。


## コヒーレント状態

単一のエネルギー固有状態は密度が時間変化しない。
振動しながら形を保つ波束を作るには、多数の準位に位相をそろえて振幅を割り振る。

:::definition[コヒーレント状態]{#def-qm-coherent}
$a\ket\alpha=\alpha\ket\alpha$、$\alpha\in\C$ を満たす規格化された状態を :term[コヒーレント状態]{en="coherent state" yomi="こひーれんとじょうたい"} という。
:::

$\ket\alpha=\sum c_n\ket n$ を代入すれば $\sqrt{n+1}c_{n+1}=\alpha c_n$。
$\sum|\alpha|^{2n}/n!=e^{|\alpha|^2}$ で規格化すると

$$
\begin{aligned}
\ket\alpha&=e^{-|\alpha|^2/2}\sum_{n=0}^\infty\frac{\alpha^n}{\sqrt{n!}}\ket n, \\
p_n&=e^{-|\alpha|^2}\frac{|\alpha|^{2n}}{n!}.
\end{aligned}
$$

平均と分散はともに $|\alpha|^2$ である。
実際に階乗を一つ、二つとずらして和をとると $\expval N=|\alpha|^2$、$\expval{N(N-1)}=|\alpha|^4$ が得られる。
固有値 $\alpha$ は複素でもよい。
$a$ はエルミート演算子ではないからである。

:ref[def-qm-coherent] の状態の各係数に $e^{-i(n+1/2)\omega t}$ を掛けると $\ket{\psi(t)}=e^{-i\omega t/2}\ket{\alpha e^{-i\omega t}}$。
$\alpha_t=\alpha e^{-i\omega t}$ と書けば

$$
\begin{aligned}
x_c&=\sqrt2\ell\,\Re\alpha_t, \\
p_c&=\sqrt2(\hbar/\ell)\Im\alpha_t, \\
\Delta X&=\ell/\sqrt2, \\
\Delta P&=\hbar/(\sqrt2\ell).
\end{aligned}
$$

分散は中心化した $a-\alpha_t$ が状態を消し、その随伴との交換子が一であることから底状態と同じ値になる。
中心は古典振動し、波束は常に :ref[thm-qm-uncertainty] の最小値をとる。
座標方程式 $(\xi+\partial_\xi)u=\sqrt2\alpha u$ を積分すれば $u\propto e^{-\xi^2/2+\sqrt2\alpha\xi}$ となり、密度が平行移動した Gaussian であることも分かる。

:::example[三つの長さの変位]
$\alpha=3/\sqrt2$ を実数に選ぶと $x_c(0)=3\ell,p_c(0)=0$。
四分の一周期後には $x_c=0,p_c=-3\hbar/\ell$ である。
一方、準位番号の平均は $9/2$ で時刻に依存しない。
:::

:::example[異なる波束の重なり]
級数で内積を計算すると $\braket\beta\alpha=\exp[-(|\alpha|^2+|\beta|^2)/2+\beta^*\alpha]$。
$\beta=0,\alpha=1$ では重なり確率は $e^{-1}$ である。
異なる固有値でも直交しないのは $a$ が正規演算子でないためである。
:::

### 変位演算子と直交基底の違い

$d=\sqrt2\ell\Re\alpha$、$p_d=\sqrt2\hbar\Im\alpha/\ell$ と置く。$K=\alpha a^\dagger-\alpha^*a=i(p_dX-dP)/\hbar$ を生成子とする変位を、座標表示で

$$
(D(s\alpha)f)(x)=e^{is p_d(x-sd/2)/\hbar}f(x-sd)
$$

と定める。移動と絶対値一の因子からユニタリであり、二つの実数 $s,s'$ の作用を合成すると $D((s+s')\alpha)$ になる。微分すると $\partial_sD(s\alpha)f=K D(s\alpha)f$ が急減少する滑らかな $f$ で成り立つので、$D(\alpha)=e^K$ という記法を使える。

座標の掛け算と微分にこの式を適用すると、$D^\dagger XD=X+d$、$D^\dagger PD=P+p_d$ を得る。従って $D^\dagger aD=a+\alpha I$ であり、底状態に $D$ を作用させると $a$ の固有値 $\alpha$ の規格化状態になる。位相の一致も $\braket{0}{D(\alpha)0}=e^{-|\alpha|^2/2}>0$ を Gaussian 積分で確認できる。二つの変位を順に作用させて位相を整理すると

$$
D(\alpha)D(\beta)=e^{(\alpha\beta^*-\alpha^*\beta)/2}D(\alpha+\beta).
$$

実際、二番目の変位 $\beta$ の位相を $x$ から $x-d_\alpha$ へ移すため、一つの変位 $D(\alpha+\beta)$ に比べて位相の指数に $i(p_\alpha d_\beta-p_\beta d_\alpha)/(2\hbar)=(\alpha\beta^*-\alpha^*\beta)/2$ が残る。
順序の違いは位相にのみ残る。
コヒーレント状態は互いに直交しないが、演習で示す恒等分解により任意の状態の展開に使える。

:::example[古典振幅と量子揺らぎを比較する]
実数 $\alpha>0$ なら最大振幅 $x_{\max}=\sqrt2\ell\alpha$ に対する位置幅の比は $\Delta X/x_{\max}=1/(2\alpha)$。
$\alpha=10$ なら $5\%$ である。
同時に $N$ の相対幅は $1/\sqrt{\expval N}=1/10$ となる。
平均励起数が大きいほど、揺らぎの絶対値を零とせず古典的な運動を観測できる。
:::

### 外力で真空からコヒーレント状態を作る

$H(t)=H_0-F(t)X$ とし、底状態から有限の時間だけ力を加える。$X=\ell(a+a^\dagger)/\sqrt2$ を使うと、Heisenberg 方程式は $\dot a=-i\omega a+i\ell F(t)/(\sqrt2\hbar)$ になる。両辺に $e^{i\omega t}$ を掛けて積を微分すれば
$$
a_H(t)=a e^{-i\omega t}+\alpha(t),\qquad
\alpha(t)=\frac{i\ell}{\sqrt2\hbar}\int_0^t e^{-i\omega(t-s)}F(s)\,\dd s .
$$
$U^\dagger aU\ket0=\alpha(t)\ket0$ だから $a\,U\ket0=\alpha(t)U\ket0$。従って時刻 $t$ の状態は全体位相を除いて $\ket{\alpha(t)}$ である。外力が与えた平均励起エネルギーは $\hbar\omega|\alpha|^2$ に等しい。持続時間 $\tau$ の一定力 $F$ なら $\alpha(\tau)=F\ell(1-e^{-i\omega\tau})/(\sqrt2\hbar\omega)$。一周期の整数倍で力を止めると励起は消え、半周期で止めると位相が加算される。力の時間積分だけから励起エネルギーを決められない理由である。

### スクイーズされた真空と周波数の急変

実数 $r$ に対し $b=a\cosh r+a^\dagger\sinh r$ とする。$[b,b^\dagger]=\cosh^2r-\sinh^2r=1$。$b\ket{0_r}=0$ の座標方程式は $(e^r\xi+e^{-r}\partial_\xi)u=0$ なので $u\propto e^{-e^{2r}\xi^2/2}$。従って $\Delta X=e^{-r}\ell/\sqrt2$、$\Delta P=e^r\hbar/(\sqrt2\ell)$、積は $\hbar/2$ のままである。逆変換 $a=b\cosh r-b^\dagger\sinh r$ を使えば $\langle a^\dagger a\rangle=\sinh^2 r$。一方の幅を小さくするとエネルギーが増える。この :term[スクイーズ状態]{en="squeezed state" yomi="すくいーずじょうたい"} は、一定の中心変位だけを加えたコヒーレント状態と区別される。


### 短い力積を加えた後の励起数

一定力を短時間だけ加える。力を加え続けている間にも波束は回転するため、単純な力積近似が正しいかを角度で判定する。

本文の一定力の解で $F\ell/(\sqrt2\hbar\omega)=1$ と選ぶと、励起数は $|1-e^{-ix}|^2=2(1-\cos x)$、瞬間的な力積なら $x^2$。

入力の目盛を $b=5.3052\times10^{-15}\,\mathrm s$、出力の目盛を $d=1$ と定める。入力を $b$ で割った値を $x$ として、平均励起数の値を $d f(x)$ と書くと

$$
f(x)=- 2 \left(\cos{\left(x \right)} - 1\right),\qquad f_{\rm app}(x)=x^{2}.
$$

となる。ここで $f_{\rm app}$ は比較に使う近似式である。例えば入力が $5.3052\times10^{-15}\,\mathrm s$ のとき、$x=1$ なので

$$
\begin{aligned}\bar n&=d f(1)=d\left(- 2 \left(-1 + \cos{\left(1 \right)}\right)\right)\simeq 0.9194,\\ \bar n_{\rm app}&=d f_{\rm app}(1)=d\left(1\right)\simeq 1.\end{aligned}
$$

同じ近似の精度をほかの目盛でも調べる。下表の差は近似値から元の式の値を引いているので、符号から過大評価か過小評価かを判定できる。

| $x$ | 入力 $[\mathrm s]$ | 元の式 $[1]$ | 近似 $[1]$ | 差 $[1]$ |
|---|---|---|---|---|
| $\frac{1}{10}$ | $5.3052\times10^{-16}$ | $0.0099917$ | $0.01$ | $8.3306\times10^{-6}$ |
| $\frac{1}{2}$ | $2.6526\times10^{-15}$ | $0.24483$ | $0.25$ | $0.0051651$ |
| $1$ | $5.3052\times10^{-15}$ | $0.9194$ | $1$ | $0.080605$ |
| $2$ | $1.061\times10^{-14}$ | $2.8323$ | $4$ | $1.1677$ |
| $3$ | $1.5915\times10^{-14}$ | $3.98$ | $9$ | $5.02$ |

入力に小さな幅があるときには、表の行を選ぶだけでは出力の幅を決められない。連鎖律により、実際の入力 $u=bx$ に対する傾きは

$$
\frac{\dd \bar n}{\dd u}=\frac{d}{b} f^\prime(x),\qquad f^\prime(x)=2 \sin{\left(x \right)}.
$$

である。入力が基準値から $1\%$ 上にずれた場合を直接計算した差と、接線を使った推定を照合する。

$$
\begin{aligned}\delta \bar n_{\rm lin}&=d f^\prime(1)\,\frac{1}{100}\simeq 0.016829,\\ \delta \bar n_{\rm exact}&=d\{f(\frac{101}{100})-f(1)\}\simeq 0.016883.\end{aligned}
$$

力せき近似の励起は単調に増えるが、厳密な状態は調和運動の一周期後に底へ戻る。この差は摩擦で振動が減衰した結果ではない。


### ばねの平衡位置を急にずらした後のエネルギー

イオンのトラップの曲率を保ったまま、底の位置だけを距離 $d$ だけ瞬時に動かす。切替え中に波動関数は変わらない。新しい底を原点にすると、もとの底状態は位置だけずれた Gaussian である。



| 設定する量 | 記号 | 計算に用いる値 |
|---|---|---|
| イオンの質量 | $m$ | $6.64\times10^{-26}\,\mathrm{kg}$ |
| トラップの通常周波数 | $\nu$ | $1\times10^{6}\,\mathrm{Hz}$ |
| 位置の急変量 | $d$ | $3\times10^{-8}\,\mathrm{m}$ |

**角振動数。** 通常の振動数から一周期の位相を掛けて変換する。


$$
\begin{aligned}\omega&=2 \pi \nu,\\ \frac{\omega}{\mathrm{s^{-1}}}&\simeq 2 \pi (1\times10^{6})\simeq 6.2832\times10^{6}.\end{aligned}
$$

**振動子長。** 運動項とポテンシャル項の尺度を合わせる。


$$
\begin{aligned}\ell&=\sqrt{\frac{\hbar}{m \omega}},\\ \frac{\ell}{\mathrm{m}}&\simeq \sqrt{\frac{(1.0545718\times10^{-34})}{(6.64\times10^{-26}) (6.2831853\times10^{6})}}\simeq 1.5899\times10^{-8}.\end{aligned}
$$

**新しい井戸での励起数。** 新しい中心から見た波束の中心は負の変位を持つ。コヒーレント振幅の二乗で個数を求める。


$$
\begin{aligned}\bar n&=\frac{d^{2}}{2 \ell^{2}},\\ \bar n&\simeq \frac{\left((3\times10^{-8})\right)^{2}}{2 \left((1.5898791\times10^{-8})\right)^{2}}\simeq 1.7803.\end{aligned}
$$

**底状態を測る確率。** コヒーレント状態の Poisson 分布の零番目を選ぶ。


$$
\begin{aligned}p_0&=e^{- \bar n},\\ p_0&\simeq e^{- 1.7802636}\simeq 0.16859.\end{aligned}
$$

**加えられた平均仕事。** 振動数は保たれ零点エネルギーは相殺するので、余分な準位間隔の平均だけが残る。


$$
\begin{aligned}W&=\bar n \hbar \omega,\\ \frac{W}{\mathrm{J}}&\simeq 1.7802636 (1.0545718\times10^{-34}) (6.2831853\times10^{6})\simeq 1.1796\times10^{-27}.\end{aligned}
$$

変位を段階的に増すと、基底状態への残存率は変位の二乗の指数関数で下がる。


| $d\;[\mathrm{m}]$ | $\bar n$ | $p_0$ | $W\;[\mathrm{J}]$ |
|---|---|---|---|
| $1\times10^{-8}$ | $0.19781$ | $0.82053$ | $1.3107\times10^{-28}$ |
| $2\times10^{-8}$ | $0.79123$ | $0.45329$ | $5.2427\times10^{-28}$ |
| $3\times10^{-8}$ | $1.7803$ | $0.16859$ | $1.1796\times10^{-27}$ |

位置の符号を逆にしても個数分布は変わらない。しかし運動の四分の一周期後に測る平均運動量は逆向きになる。


### 周波数を二倍にした振動子で出る偶数の励起

トラップ中の電子の有効質量を $m_e$ とし、中心を動かさずに角周波数を $\omega_0=2\pi\times10^6\,\mathrm{s^{-1}}$ から $\omega_1=2\omega_0$ に急変させる。操作時間を両周期より短くして波形を凍結する。旧基底は新しい数演算子のどの成分に移るだろうか。エネルギー平均だけでなく離散的な検出分布を求める。

$$
a_j=\sqrt{\frac{m\omega_j}{2\hbar}}X+\frac{iP}{\sqrt{2m\hbar\omega_j}},\quad r=\frac12\ln(\omega_1/\omega_0),\quad a_0=a_1\cosh r-a_1^\dagger\sinh r.
$$

新しい消滅・生成演算子で $X,P$ を書いて $a_0$ に代入すると、係数は $(e^r+e^{-r})/2$ と $-(e^r-e^{-r})/2$ になる。旧真空の展開を $\ket{0_0}=\sum c_n\ket{n_1}$ として、$a_0\ket{0_0}=0$ の各成分を比べる。

$$
\sqrt{n+1}c_{n+1}\cosh r=\sqrt n\,c_{n-1}\sinh r,\qquad c_{-1}=0.
$$

$n=0$ から $c_1=0$。漸化式を二段ずつ使うと全奇数成分が消え、偶数成分は $c_{2k}=c_0\sqrt{(2k)!}(\tanh r)^k/(2^k k!)$。規格化には二項級数 $\sum_k\binom{2k}{k}y^k/4^k=(1-y)^{-1/2}$ を $y=\tanh^2 r$ に適用する。従って $|c_0|^2=1/\cosh r$ である。

| 新しい励起数 | 確率 | 累積確率 | 測定エネルギー $[\mathrm J]$ |
|---|---|---|---|
| 0 | $0.942809$ | $0.942809$ | $6.6261\times10^{-28}$ |
| 2 | $0.0523783$ | $0.9951873$ | $3.313\times10^{-27}$ |
| 4 | $0.00436486$ | $0.9995522$ | $5.9635\times10^{-27}$ |
| 6 | $4.04153\times10^{-4}$ | $0.9999563$ | $8.6139\times10^{-27}$ |
| 8 | $3.92927\times10^{-5}$ | $0.9999956$ | $1.1264\times10^{-26}$ |
| 10 | $3.92927\times10^{-6}$ | $0.9999996$ | $1.3915\times10^{-26}$ |

確率母関数 $G(z)=\sum p_n z^n$ は、先ほどの和で $y=z^2\tanh^2 r$ と置けばよい。微分の際に $G^{\prime\prime}(1)$ は $\langle n(n-1)\rangle$ を返すから、分散を作るときに一次微分を戻す。

$$
G(z)=\frac{1}{\cosh r\sqrt{1-z^2\tanh^2r}},\quad \langle n\rangle=\sinh^2r,\quad (\Delta n)^2=2\langle n\rangle(1+\langle n\rangle).
$$

ここでは $\tanh r=1/3$、$\langle n\rangle=1/8$、$(\Delta n)^2=9/32$。新しい Hamilton 演算子の平均は $\hbar\omega_1(1/2+1/8)=8.2826\times10^{-28}\,\mathrm J$。初期の位置・運動量幅を新しいポテンシャルのエネルギーへ代入しても同じになる。

$$
\frac{\langle P^2\rangle_0}{2m}+\frac{m\omega_1^2\langle X^2\rangle_0}{2}=\frac{\hbar\omega_0}{4}+\frac{\hbar\omega_1^2}{4\omega_0}.
$$

旧真空は新しい固有状態一つにはならないので、待ってから測ると位置分布の幅は振動する。$X(t)=X\cos\omega_1t+[P/(m\omega_1)]\sin\omega_1t$ と初期共分散零から、幅の二乗は $\hbar[\cos^2\omega_1t+(\omega_0/\omega_1)^2\sin^2\omega_1t]/(2m\omega_0)$。最大幅 $3.0352\times10^{-6}\,\mathrm m$ の半分まで縮み、一方で数分布は時間一定である。


### 平衡位置と平均運動量を閉じた順路で動かす

位置と運動量の平均を元に戻す四つの操作を考える。振動子の一方の干渉経路だけでこれを行うと、両経路の波束の中心が一致していても再合成後の明暗が動く。操作が短く、その間の自由回転を無視するか補償したものとする。

$D(\alpha)=\exp(\alpha a^\dagger-\alpha^*a)$ とし、$A=\alpha a^\dagger-\alpha^*a,\ B=\beta a^\dagger-\beta^*a$ とおく。$[A,B]=\alpha\beta^*-\alpha^*\beta$ は数なので $e^{sA}Be^{-sA}=B+s[A,B]$。$F(s)=e^{sA}e^{sB}$ を微分すると $F^\prime=(A+B+s[A,B])F$、初期値は $I$。右辺の生成子は異なる時点で可換であるから積分し、$s=1$ を入れると


$$
D(\alpha)D(\beta)=e^{(\alpha\beta^*-\alpha^*\beta)/2}D(\alpha+\beta),\quad D(\alpha)D(\beta)D(-\alpha)D(-\beta)=e^{\alpha\beta^*-\alpha^*\beta}I.\label{eq-qm-displacement-phase}
$$

順序を逆にすれば指数の符号が変わる。単独の振動子で位置分布だけを見るとこの位相は共通で消える。一方、二経路の片側で四操作を行い、もう片側で位相を保つと、周回が作る位相は経路間の相対位相になる。

```tikz caption="無次元振幅の複素平面。右端の演算子から作用させると四つの変位で閉路ができる。"
\begin{tikzpicture}[>=Stealth,x=1cm,y=1cm]
\draw[->] (-2.0,0)--(0.9,0) node[right] {$\mathrm{Re}\,\alpha$};\draw[->] (0,-2.2)--(0,0.8) node[above] {$\mathrm{Im}\,\alpha$};
\draw[thick,->] (0,0)--(0,-1.6) node[midway,right] {$-\beta$};
\draw[thick,->] (0,-1.6)--(-1.4,-1.6) node[midway,below] {$-\alpha$};
\draw[thick,->] (-1.4,-1.6)--(-1.4,0) node[midway,left] {$\beta$};
\draw[thick,->] (-1.4,0)--(0,0) node[midway,above] {$\alpha$};
\end{tikzpicture}
```


## 演習問題

::::exercise{level=basic}
$n=2$ から一回下げ、再び一回上げた結果を求めよ。

:::solution
$a\ket2=\sqrt2\ket1$、続いて $a^\dagger\sqrt2\ket1=2\ket2$。
これは $N\ket2=2\ket2$ に等しい。
:::
::::

::::exercise{level=basic}
$\ket1$ の $\Delta X\Delta P$ を求めよ。

:::solution
両分散に $n+1/2=3/2$ が掛かるので積は $3\hbar/2$ である。
:::
::::

::::exercise{level=basic}
$H_3$ を求め、三つの節を示せ。

:::solution
$H_3=2\xi(4\xi^2-2)-4(2\xi)=8\xi^3-12\xi$。
従って $x=0,\pm\sqrt{3/2}\ell$ が節である。
:::
::::

::::exercise{level=basic}
$\alpha=2i$ の $p_0,p_1$ と初期平均運動量を求めよ。

:::solution
$|\alpha|^2=4$ だから $p_0=e^{-4},p_1=4e^{-4}$。
$\Im\alpha=2$ を用い $p_c=2\sqrt2\hbar/\ell$ である。
:::
::::

::::exercise{level=advanced}
$\ket n$ で $\expval{X^4}$ を求めよ。

:::solution
$(a+a^\dagger)^2\ket n$ の $n-2,n,n+2$ 成分は $\sqrt{n(n-1)},2n+1,\sqrt{(n+1)(n+2)}$。
ノルム二乗を足し $\ell^4/4$ を掛けると $3\ell^4(2n^2+2n+1)/4$ となる。
存在しない下の準位の項は零とする。
:::
::::

::::exercise{level=advanced}
温度 $T$ の熱浴と平衡にある振動子について、$\beta=1/(k_{\rm B}T)$、$p_n=e^{-\beta E_n}/Z$ として平均エネルギーを求めよ。

:::solution
$q=e^{-\beta\hbar\omega}$ として $Z=e^{-\beta\hbar\omega/2}/(1-q)$。
$-\partial_\beta\ln Z$ を計算すると $\hbar\omega[1/2+1/(e^{\beta\hbar\omega}-1)]$ で、温度依存部分が黒体放射に寄与する。
:::
::::

::::exercise{level=advanced}
$D(\alpha)=\exp(\alpha a^\dagger-\alpha^*a)$ が底状態をコヒーレント状態に移すことを示せ。

:::solution
座標表示の変位は $f(x-d)$ に位相 $e^{ip_d(x-d/2)/\hbar}$ を掛ける。従って位置には $d$、運動量には位相を微分して得る $p_d$ が加わる。:eqref[eq-qm-oscillator] の定義に戻すと $D^\dagger aD=a+\alpha I$ であり、$aD\ket0=\alpha D\ket0$ となる。座標表示の変位はノルムを保つ。
:::
::::

::::exercise{level=advanced}
コヒーレント状態による恒等分解 $\int\ketbra\alpha\alpha\,\dd^2\alpha/\pi=I$ を示せ。

:::solution
$\alpha=re^{i\theta}$ と置くと $\ketbra n m$ の係数の角積分は $2\pi\delta_{nm}$。
残る積分は $2\int_0^\infty e^{-r^2}r^{2n+1}\dd r/n!=1$ だからすべての対角要素が一となる。
:::
::::

::::exercise{level=advanced}
振動子を角振動数 $\omega_i$ の底状態に準備し、周波数を瞬時に $\omega_f$ に変える。
直後に新しい底状態が得られる確率と平均エネルギーを求めよ。

:::solution
$\ell_j=\sqrt{\hbar/(m\omega_j)}$ と書く。
直後の波形は初期の Gaussian のままである。
二つの底状態の内積は Gaussian 積分より

$$
\begin{aligned}
\braket{0_f}{0_i}&=\sqrt{\frac{2\ell_i\ell_f}{\ell_i^2+\ell_f^2}}, \\
P_{0_f}&=\frac{2\sqrt{\omega_i\omega_f}}{\omega_i+\omega_f}.
\end{aligned}
$$

初期の二乗平均 $\expval{P^2}=m\hbar\omega_i/2$、$\expval{X^2}=\hbar/(2m\omega_i)$ を新しい $H$ に入れると $\expval{H_f}=\hbar(\omega_i+\omega_f^2/\omega_i)/4$。
新しい底との差は $\hbar(\omega_i-\omega_f)^2/(4\omega_i)\ge0$ である。
波形のパリティは変わらないため、励起確率があるのは偶数 $n$ のみである。
:::
::::

::::exercise{level=advanced}
$\alpha$ が実のとき、二つのコヒーレント状態の重ね合わせ $\ket\alpha\pm\ket{-\alpha}$ の規格化と準位番号のパリティを求めよ。

:::solution
$\braket\alpha{-\alpha}=e^{-2\alpha^2}$ なので規格化因子は $[2(1\pm e^{-2\alpha^2})]^{-1/2}$。
負号側には $\alpha\ne0$ が必要である。
個数基底を展開すると $\ket n$ の係数には $\alpha^n[1\pm(-1)^n]$ が現れる。
従って正号側は偶数準位だけ、負号側は奇数準位だけからなる。

それぞれ一周期の半分だけ発展させると $\alpha\mapsto-\alpha$ となる。
二つの波束が入れ替わるので、重ね合わせ全体の射線はその時点でも最初と同じになる。
これは単一のコヒーレント状態が半周期で反対側へ移ることと整合する。
:::
::::

::::exercise{level=advanced}
次の計算には :eqref[eq-qm-displacement-phase] を使える。$\alpha=1/2,\ \beta=i\pi/4$ とし、上式の順で一方の経路に作用させる。等分した二経路を同位相で再合成すると明側の確率はいくらか。逆向きの閉路と確率だけで見分けられるか。

:::solution
$\alpha\beta^*-\alpha^*\beta=-i\pi/4$ なので、この経路だけに $e^{-i\pi/4}$ が掛かる。明側の振幅は $(1+e^{-i\pi/4})/2$、確率は $(1+1/\sqrt2)/2$。逆回りは位相が共役になり同じ確率である。基準側に既知の $\pi/2$ の位相を加えれば、余弦の代わりに正弦を読み、向きを判別できる。
:::
::::


## 確認テスト

::::quiz{answer="0.5" tol="1e-6"}
零点エネルギーを $\hbar\omega$ 単位で答えよ。

:::explanation
$a\ket0=0$ と因数分解から $E_0/(\hbar\omega)=1/2$ である。
:::
::::

::::quiz{answer="2" tol="1e-6"}
$a^\dagger\ket3$ の規格化された $\ket4$ に対する係数はいくつか。

:::explanation
生成演算子の係数は $\sqrt{n+1}$ なので $\sqrt4=2$ である。
:::
::::

::::quiz
底状態の座標波動関数の形はどれか。

- [ ] 正弦波
- [x] Gaussian 関数
- [ ] 両側の単純指数関数

:::explanation
一階方程式 $u_0'=-\xi u_0$ を積分すると $e^{-\xi^2/2}$ に比例する。
:::
::::

::::quiz
コヒーレント状態で時間がたつと変化するものはどれか。

- [ ] エネルギー分布
- [ ] 波束の位置幅
- [x] 零でない振幅の偏角

:::explanation
固有状態の位相が準位番号に比例してずれるため $\alpha_t=\alpha e^{-i\omega t}$。
絶対値と位置幅は一定である。
:::
::::

::::quiz{answer="10" tol="0.0001"}
三次元等方振動子で全励起数 $N=3$ をもつ空間状態はいくつか。

:::explanation
$n_x=0,1,2,3$ に対し $n_y+n_z=3-n_x$ の解数は $4,3,2,1$ なので和は $10$。スピンの自由度は含まない。
:::
::::

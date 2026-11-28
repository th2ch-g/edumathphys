---
title: 物質中の光と分散
description: 電荷の運動と因果応答から複素屈折率を求め、パルスの分散と粒子による散乱を計算する。
notation:
  - symbol: '\chi'
    meaning: 電気感受率
  - symbol: '\tilde n=n+i\kappa'
    meaning: 複素屈折率
  - symbol: '\alpha'
    meaning: 強度の吸収係数
  - symbol: '\beta_j'
    meaning: 伝播定数の角周波数に関する第 j 微分
  - symbol: '\sigma_{\mathrm{sca}}'
    meaning: 散乱断面積
---

## この章の目標

- 束縛電荷と自由電荷の応答から複素誘電率を導く。
- 因果律が実部と虚部を結び付けることを示す。
- 群遅延とパルス幅の変化を計算する。
- 小粒子の散乱と球の多重極散乱の関係を説明する。

## Lorentz 振動子モデル

透明に見える物質でも、電場をかけた瞬間に分極が最大になるとは限らない。束縛電荷の慣性と摩擦を残すと、同じ機構から位相のずれと吸収を得られる。

### 一つの共鳴とその符号
独立な電荷 $q$、質量 $m$、数密度 $N$ を考え、変位 $x$ は一次元、ばね定数は $m\omega_0^2$、摩擦は $m\gamma\dot x$ とする。弱い光で
$$
\ddot x+\gamma\dot x+\omega_0^2 x=(q/m)E(t).
$$
共鳴角周波数 $\omega_0$ と減衰率 $\gamma$ は $\mathrm{s^{-1}}$ を単位とする。$P=Nqx$ とし、$\eps_0\chi E$ と同一視する。採用した時間因子のもとで $\dot x\mapsto-i\omega x$ なので
$$
\chi(\omega)=\frac{F}{\omega_0^2-\omega^2-i\gamma\omega},\qquad F=\frac{Nq^2}{\eps_0m}.
\label{eq-optics-lorentz}
$$
これが :term[Lorentz 振動子]{en="Lorentz oscillator" yomi="ろーれんつしんどうし"} の感受率である。分母に共役を掛けると $\chi''=F\gamma\omega/[(\omega_0^2-\omega^2)^2+\gamma^2\omega^2]$ となる。$F,\gamma,\omega>0$ なら正である。摩擦が受け取る単位体積当たりの平均仕事率は $Nm\gamma\omega^2|x|^2/2=(\omega\eps_0/2)\chi''|E|^2$ であり、電場から物質へ正のエネルギーが渡る。独立な複数の振動子では分極の加法性から $\chi$ を加算する。

### 時間応答を復元する
$\gamma<2\omega_0$ とし $\omega_d=\sqrt{\omega_0^2-\gamma^2/4}$ と置く。微小時刻まで静止した電荷を短い電場で蹴ると、速度は電場の時間積分に $(q/m)$ を掛けた量だけ飛ぶ。この初期値で同次式を解くと、分極の畳み込み核は $h(t)=F e^{-\gamma t/2}\sin(\omega_d t)/\omega_d$（$t>0$）である。$\int_0^\infty h(t)e^{i\omega t}\dd t$ で正弦を二指数に分けて積分すれば :eqref[eq-optics-lorentz] になる。応答は過去の電場だけで決まる。

:::example[零周波数と共鳴]
$\omega=0$ では $\chi=F/\omega_0^2$ が実数で、静的変位は電場と同相である。$\omega=\omega_0$ では $\chi=iF/(\gamma\omega_0)$ であり、実場の分極の山は電場の山より四分の一周期遅れる。
:::
:::example[電子とイオン]
同じ $\omega_0,\gamma$ と電荷の大きさを仮定し、一方の質量を二倍にすると $F$ は半分になる。符号を反転しただけの電荷では $F$ は変わらない。実材料の共鳴を比べるには質量だけでなく復元力と減衰も調べる必要がある。
:::

### 強度でなく平均仕事率を積分する
自由な減衰振動のエネルギー $\mathcal E=(m/2)(\dot x^2+\omega_0^2 x^2)$ を微分すると $\dot{\mathcal E}=qE\dot x-m\gamma\dot x^2$ である。定常正弦応答で周期平均の左辺はゼロとなり、入力された平均仕事が摩擦に等しい。駆動を止めると不足減衰の振幅は $e^{-\gamma t/2}$ を持ち、振幅寿命と平均エネルギーの寿命が異なる。位相規約を $e^{-i\omega t}$ から逆に変えれば $\chi$ の虚部も反転するが、仕事率は変わらない。

### 密な媒質での局所電場
独立な振動子模型では、各電荷を駆動する場を巨視的な $\vb E$ とした。等方な密な媒質の球形空洞近似では、分極が $\vb P$ で一様な周囲を除去したとき、球の内面にできる表面電荷は $-\vb P\cdot\hat{\vb r}$ である。その中心で各面積素片の場を積分すると $\vb E_{\rm loc}-\vb E=(4\pi\eps_0)^{-1}\int \hat{\vb r}(\vb P\cdot\hat{\vb r})\dd\Omega=\vb P/(3\eps_0)$。最後の等号は $\int \hat r_i\hat r_j\dd\Omega=(4\pi/3)\delta_{ij}$ による。球内の近傍粒子が作る余分な場を平均でゼロと仮定し、単一粒子の分極率を $\alpha_a$ とすれば $P=N\alpha_a(E+P/(3\eps_0))$。この連立式を解いて
$$
\frac{\eps_r-1}{\eps_r+2}=\frac{N\alpha_a}{3\eps_0}.
$$
希薄極限では $\chi\simeq N\alpha_a/\eps_0$ へ戻る。特定の格子方向や短距離相関が支配する媒質では、球内の近接場を別に求める必要がある。この関係を :term[Clausius–Mossotti の式]{en="Clausius–Mossotti relation" yomi="くらうじうすもそってぃのしき"} という。球対称の幾何が与える係数 $1/3$ は、この章の後半で扱う小球の静電分極にも現れる。

:::example[時間の減衰尺度を読む]
不足減衰で $\gamma=2.0\times10^8\,\mathrm{s^{-1}}$ なら、自由振動の包絡線が $1/e$ に下がる時間は $2/\gamma=10\,\mathrm{ns}$。振幅二乗の包絡線は $1/\gamma=5\,\mathrm{ns}$ で $1/e$ になる。線幅を寿命へ換算するときは、どちらの量が緩和しているかを確かめる。
:::

### 共鳴の極と自由減衰を結び付ける

無源の変位の指数を $e^{-i\omega t}$ に選ぶと、Lorentz 分母の零点は $\omega=-i\gamma/2\pm\sqrt{\omega_0^2-\gamma^2/4}$。$\gamma>0$ のとき両方の極は下半面にある。過減衰で平方根が純虚数でも、その絶対値は $\gamma/2$ より小さいからである。従って遅い成分を含めて自由変位は減衰する。一方、変位の強制振動の絶対値を最大にするには分母の二乗 $(\omega_0^2-\omega^2)^2+\gamma^2\omega^2$ を最小にする。$\omega^2$ で微分して $\omega_{\rm max}^2=\omega_0^2-\gamma^2/2$。これは正のときだけ有限周波数のピークになる。自由振動の周波数と強制変位のピークを、吸収パワーのピークからも区別する。外部電場の振幅が一定なら平均吸収は $\omega\chi^{\prime\prime}$ に比例する。この量は $\omega>0$ で $F\gamma/[(\omega_0^2/\omega-\omega)^2+\gamma^2]$ だから、$\gamma>0$ なら $\omega=\omega_0$ で最大になる。

:::example[細い線を共鳴中心で駆動]
減衰の小さいモデル振動子を無減衰の角周波数で駆動する。

入力は無減衰の角周波数 $\omega_0\simeq 1\times 10^{15}\,\mathrm{s^{-1}}$、減衰定数 $\gamma\simeq 1\times 10^{13}\,\mathrm{s^{-1}}$、駆動角周波数 $\omega\simeq 1\times 10^{15}\,\mathrm{s^{-1}}$ である。

二次方程式の根で減衰と自由振動をまとめる。

$$
\begin{aligned}
\omega_+&=- \frac{i \gamma}{2} + \sqrt{- \frac{\gamma^{2}}{4} + \omega_0^{2}},\\
\omega_+/(\mathrm{s^{-1}})&\simeq - \frac{i \cdot \left(1\times 10^{13}\right)}{2} + \sqrt{- \frac{\left(\left(1\times 10^{13}\right)\right)^{2}}{4} + \left(\left(1\times 10^{15}\right)\right)^{2}}\simeq 9.9999\times 10^{14}-5\times 10^{12}i.
\end{aligned}
$$

変位応答の分母が最小になる周波数を求める。

$$
\begin{aligned}
\omega_{\rm max}&=\sqrt{\max\left(0, - \frac{\gamma^{2}}{2} + \omega_0^{2}\right)},\\
\omega_{\rm max}/(\mathrm{s^{-1}})&\simeq \sqrt{\max\left(0, - \frac{\left(\left(1\times 10^{13}\right)\right)^{2}}{2} + \left(\left(1\times 10^{15}\right)\right)^{2}\right)}\simeq 9.9997\times 10^{14}.
\end{aligned}
$$

Lorentz 分母の有理化で同相応答を出す。

$$
\begin{aligned}
\chi'/\chi(0)&=\frac{\omega_0^{2} \left(\omega_0^{2} - \omega^{2}\right)}{\gamma^{2} \omega^{2} + \left(\omega_0^{2} - \omega^{2}\right)^{2}},\\
\chi'/\chi(0)&\simeq \frac{\left(\left(1\times 10^{15}\right)\right)^{2} \cdot \left(\left(\left(1\times 10^{15}\right)\right)^{2} - \left(\left(1\times 10^{15}\right)\right)^{2}\right)}{\left(\left(1\times 10^{13}\right)\right)^{2} \cdot \left(\left(1\times 10^{15}\right)\right)^{2} + \left(\left(\left(1\times 10^{15}\right)\right)^{2} - \left(\left(1\times 10^{15}\right)\right)^{2}\right)^{2}}=0.
\end{aligned}
$$

遅れの成分は分母の共役によって正となる。

$$
\begin{aligned}
\chi''/\chi(0)&=\frac{\gamma \omega_0^{2} \omega}{\gamma^{2} \omega^{2} + \left(\omega_0^{2} - \omega^{2}\right)^{2}},\\
\chi''/\chi(0)&\simeq \frac{\left(1\times 10^{13}\right) \cdot \left(\left(1\times 10^{15}\right)\right)^{2} \cdot \left(1\times 10^{15}\right)}{\left(\left(1\times 10^{13}\right)\right)^{2} \cdot \left(\left(1\times 10^{15}\right)\right)^{2} + \left(\left(\left(1\times 10^{15}\right)\right)^{2} - \left(\left(1\times 10^{15}\right)\right)^{2}\right)^{2}}=100.
\end{aligned}
$$

同相成分と直交位相成分を合成する。

$$
\begin{aligned}
A/A(0)&=\sqrt{(\chi'/\chi(0))^{2} + (\chi''/\chi(0))^{2}},\\
A/A(0)&= \sqrt{\left(0\right)^{2} + \left(100\right)^{2}}=100.
\end{aligned}
$$

強制変位は駆動と直交位相で、自由振動とピークのずれは小さい。
:::

:::example[広い線の低周波側を観測]
減衰を増やし、無減衰の共鳴より下の位置を調べる。

入力は無減衰の角周波数 $\omega_0\simeq 1\times 10^{15}\,\mathrm{s^{-1}}$、減衰定数 $\gamma\simeq 7\times 10^{14}\,\mathrm{s^{-1}}$、駆動角周波数 $\omega\simeq 8\times 10^{14}\,\mathrm{s^{-1}}$ である。

二次方程式の根で減衰と自由振動をまとめる。

$$
\begin{aligned}
\omega_+&=- \frac{i \gamma}{2} + \sqrt{- \frac{\gamma^{2}}{4} + \omega_0^{2}},\\
\omega_+/(\mathrm{s^{-1}})&\simeq - \frac{i \cdot \left(7\times 10^{14}\right)}{2} + \sqrt{- \frac{\left(\left(7\times 10^{14}\right)\right)^{2}}{4} + \left(\left(1\times 10^{15}\right)\right)^{2}}\simeq 9.3675\times 10^{14}-3.5\times 10^{14}i.
\end{aligned}
$$

変位応答の分母が最小になる周波数を求める。

$$
\begin{aligned}
\omega_{\rm max}&=\sqrt{\max\left(0, - \frac{\gamma^{2}}{2} + \omega_0^{2}\right)},\\
\omega_{\rm max}/(\mathrm{s^{-1}})&\simeq \sqrt{\max\left(0, - \frac{\left(\left(7\times 10^{14}\right)\right)^{2}}{2} + \left(\left(1\times 10^{15}\right)\right)^{2}\right)}\simeq 8.6891\times 10^{14}.
\end{aligned}
$$

Lorentz 分母の有理化で同相応答を出す。

$$
\begin{aligned}
\chi'/\chi(0)&=\frac{\omega_0^{2} \left(\omega_0^{2} - \omega^{2}\right)}{\gamma^{2} \omega^{2} + \left(\omega_0^{2} - \omega^{2}\right)^{2}},\\
\chi'/\chi(0)&\simeq \frac{\left(\left(1\times 10^{15}\right)\right)^{2} \cdot \left(\left(\left(1\times 10^{15}\right)\right)^{2} - \left(\left(8\times 10^{14}\right)\right)^{2}\right)}{\left(\left(7\times 10^{14}\right)\right)^{2} \cdot \left(\left(8\times 10^{14}\right)\right)^{2} + \left(\left(\left(1\times 10^{15}\right)\right)^{2} - \left(\left(8\times 10^{14}\right)\right)^{2}\right)^{2}}=\frac{225}{277}.
\end{aligned}
$$

遅れの成分は分母の共役によって正となる。

$$
\begin{aligned}
\chi''/\chi(0)&=\frac{\gamma \omega_0^{2} \omega}{\gamma^{2} \omega^{2} + \left(\omega_0^{2} - \omega^{2}\right)^{2}},\\
\chi''/\chi(0)&\simeq \frac{\left(7\times 10^{14}\right) \cdot \left(\left(1\times 10^{15}\right)\right)^{2} \cdot \left(8\times 10^{14}\right)}{\left(\left(7\times 10^{14}\right)\right)^{2} \cdot \left(\left(8\times 10^{14}\right)\right)^{2} + \left(\left(\left(1\times 10^{15}\right)\right)^{2} - \left(\left(8\times 10^{14}\right)\right)^{2}\right)^{2}}=\frac{350}{277}.
\end{aligned}
$$

同相成分と直交位相成分を合成する。

$$
\begin{aligned}
A/A(0)&=\sqrt{(\chi'/\chi(0))^{2} + (\chi''/\chi(0))^{2}},\\
A/A(0)&= \sqrt{\left(\frac{225}{277}\right)^{2} + \left(\frac{350}{277}\right)^{2}}\simeq 1.5021.
\end{aligned}
$$

同相成分と吸収成分がともに残る。自由振動の実部と強制変位のピークも区別できる。
:::

:::example[過減衰の光学応答]
減衰が特に大きい模型の低周波側を評価する。

入力は無減衰の角周波数 $\omega_0\simeq 1\times 10^{15}\,\mathrm{s^{-1}}$、減衰定数 $\gamma\simeq 3\times 10^{15}\,\mathrm{s^{-1}}$、駆動角周波数 $\omega\simeq 2\times 10^{14}\,\mathrm{s^{-1}}$ である。

二次方程式の根で減衰と自由振動をまとめる。

$$
\begin{aligned}
\omega_+&=- \frac{i \gamma}{2} + \sqrt{- \frac{\gamma^{2}}{4} + \omega_0^{2}},\\
\omega_+/(\mathrm{s^{-1}})&\simeq - \frac{i \cdot \left(3\times 10^{15}\right)}{2} + \sqrt{- \frac{\left(\left(3\times 10^{15}\right)\right)^{2}}{4} + \left(\left(1\times 10^{15}\right)\right)^{2}}\simeq -3.8197\times 10^{14}i.
\end{aligned}
$$

変位応答の分母が最小になる周波数を求める。

$$
\begin{aligned}
\omega_{\rm max}&=\sqrt{\max\left(0, - \frac{\gamma^{2}}{2} + \omega_0^{2}\right)},\\
\omega_{\rm max}/(\mathrm{s^{-1}})&\simeq \sqrt{\max\left(0, - \frac{\left(\left(3\times 10^{15}\right)\right)^{2}}{2} + \left(\left(1\times 10^{15}\right)\right)^{2}\right)}=0.
\end{aligned}
$$

Lorentz 分母の有理化で同相応答を出す。

$$
\begin{aligned}
\chi'/\chi(0)&=\frac{\omega_0^{2} \left(\omega_0^{2} - \omega^{2}\right)}{\gamma^{2} \omega^{2} + \left(\omega_0^{2} - \omega^{2}\right)^{2}},\\
\chi'/\chi(0)&\simeq \frac{\left(\left(1\times 10^{15}\right)\right)^{2} \cdot \left(\left(\left(1\times 10^{15}\right)\right)^{2} - \left(\left(2\times 10^{14}\right)\right)^{2}\right)}{\left(\left(3\times 10^{15}\right)\right)^{2} \cdot \left(\left(2\times 10^{14}\right)\right)^{2} + \left(\left(\left(1\times 10^{15}\right)\right)^{2} - \left(\left(2\times 10^{14}\right)\right)^{2}\right)^{2}}=\frac{200}{267}.
\end{aligned}
$$

遅れの成分は分母の共役によって正となる。

$$
\begin{aligned}
\chi''/\chi(0)&=\frac{\gamma \omega_0^{2} \omega}{\gamma^{2} \omega^{2} + \left(\omega_0^{2} - \omega^{2}\right)^{2}},\\
\chi''/\chi(0)&\simeq \frac{\left(3\times 10^{15}\right) \cdot \left(\left(1\times 10^{15}\right)\right)^{2} \cdot \left(2\times 10^{14}\right)}{\left(\left(3\times 10^{15}\right)\right)^{2} \cdot \left(\left(2\times 10^{14}\right)\right)^{2} + \left(\left(\left(1\times 10^{15}\right)\right)^{2} - \left(\left(2\times 10^{14}\right)\right)^{2}\right)^{2}}=\frac{125}{267}.
\end{aligned}
$$

同相成分と直交位相成分を合成する。

$$
\begin{aligned}
A/A(0)&=\sqrt{(\chi'/\chi(0))^{2} + (\chi''/\chi(0))^{2}},\\
A/A(0)&= \sqrt{\left(\frac{200}{267}\right)^{2} + \left(\frac{125}{267}\right)^{2}}\simeq 0.88333.
\end{aligned}
$$

極は虚軸にあり、変位の最大は静的な駆動で得られる。
:::

## 複素屈折率と吸収係数

### 屈折率の枝と二つの減衰距離
非磁性媒質で波面を $e^{ik_0\tilde n z}$ と書くと $\tilde n^2=1+\chi$ である。外向きで減衰する解を取ると $\tilde n=n+i\kappa$、$\kappa\ge0$ である。$\eps'=1+\chi'$、$\eps''=\chi''>0$ とすると $n^2-\kappa^2=\eps'$ と $2n\kappa=\eps''$ が成り立つ。さらに $n^2+\kappa^2=\sqrt{(\eps')^2+(\eps'')^2}$ なので両成分を個別に復元できる。$\eps''\ll\eps'$ かつ $\eps'>0$ なら $\kappa\simeq\eps''/(2\sqrt{\eps'})$ である。

磁場と電場の比は一様媒質内で一定なので、法線への平均エネルギー流は電場振幅の二乗と同じ率で減衰する。:term[吸収係数]{en="absorption coefficient" yomi="きゅうしゅうけいすう"} $\alpha$ は
$$
I(z)=I(0)e^{-\alpha z},\qquad\alpha=2k_0\kappa=\frac{4\pi\kappa}{\lambda_0}.
\label{eq-optics-absorption}
$$
$\alpha$ の単位は $\mathrm{m^{-1}}$ であり、強度の $1/e$ 長は $\alpha^{-1}$、振幅の $1/e$ 長は $2\alpha^{-1}$ である。前進平面波では $I=n|E|^2/(2\eta_0)$ である。時間平均した Poynting の定理で分極電流 $-i\omega\eps_0\chi E$ の仕事を取ると、$-\dd I/\dd z=(\omega\eps_0/2)\chi''|E|^2$ である。式 :eqref[eq-optics-absorption] へ $2n\kappa=\chi''$ を入れても同じ仕事率を得る。

:::example[弱い吸収のスケール]
$\lambda_0=500\,\mathrm{nm},\kappa=10^{-5}$ なら $\alpha=80\pi\,\mathrm{m^{-1}}$ で、強度の $1/e$ 長は約 $4.0\,\mathrm{mm}$ である。
:::
:::example[透過率と厚さ]
界面反射を取り除いた試料で厚さ $2d$ の透過率が $1/4$ なら、厚さ $d$ では $1/2$ になる。反射があり単一通過で近似する場合、厚さに依存しない界面透過率が前に掛かる。ただし、この近似を実際の厚さゼロまで延長すると、二面の振幅が干渉する効果を落とす。
:::

:::example[複素誘電率から場の減衰へ]
$\eps_r=3+4i$ なら $|\eps_r|=5$ であり、$n=\sqrt{(5+3)/2}=2,\kappa=\sqrt{(5-3)/2}=1$。真空波長 $600\,\mathrm{nm}$ での強度の $1/e$ 長は $600/(4\pi)\,\mathrm{nm}\simeq47.7\,\mathrm{nm}$。空気側から入ったパワー全体のうち、界面に入る割合を知るには :eqref[eq-optics-metal] も用いる。
:::

### 背景光を含む吸収測定
吸収を光学濃度 $A_{10}=-\log_{10}(I/I_0)$ で記録すると、界面損失を除いた透過の傾きは $A_{10}=\alpha L/\ln 10$。反射のある同じ材質の試料を厚さだけ変え、界面透過が厚さによらず、内部多重反射が無視できるとすれば、$I(L_2)/I(L_1)=e^{-\alpha(L_2-L_1)}$ で界面の定数因子が消える。試料内に定在波が生じる薄膜では、この単純な消去ができない。層の振幅を重ねて全透過を求める必要がある。

:::example[厚さを変えて吸収係数を読む]
同じ界面を持つ試料で厚さを $1\,\mathrm{mm}$ から $3\,\mathrm{mm}$ に変えたとき、透過が半分になったとする。この近似で $\alpha=\ln2/(2\,\mathrm{mm})\simeq347\,\mathrm{m^{-1}}$ である。
:::

### 表面反射が吸収測定に入る場合

反射率 $R$ の二つの平行面を持つ吸収板で、内部の一方向の透過率を $q=e^{-\alpha L}$ とする。界面自体の吸収を無視でき、往復干渉が帯域平均で消える場合、各透過パワーに対し最初の寄与は $(1-R)^2q$、一往復ごとの比は $R^2q^2$。等比級数より $T=(1-R)^2q/(1-R^2q^2)$ である。このモデルは両側の面で同じパワー分割を仮定する。具体的な複素屈折率の強吸収媒質では、界面外と内部の接線場を使って薄膜の複素振幅を求めてから平均する。位相がよく保たれる板では、強度の等比級数では干渉縞を再現できない。

:::example[透明な板]
内部吸収の無視できる板を透過で評価する。

入力は各面の反射率 $R=\frac{1}{25}$、吸収係数 $\alpha=0\,\mathrm{m^{-1}}$、厚さ $L=\frac{1}{1000}\,\mathrm m$ である。

表面での反射を除き、内部のパワー減衰を求める。

$$
\begin{aligned}
q&=e^{- L \alpha},\\
q&= e^{- \left(\frac{1}{1000}\right) \cdot \left(0\right)}=1.
\end{aligned}
$$

往復時に二回ずつ界面反射と内部通過がある。

$$
\begin{aligned}
g&=R^{2} q^{2},\\
g&= \left(\frac{1}{25}\right)^{2} \cdot \left(1\right)^{2}=\frac{1}{625}.
\end{aligned}
$$

前面と背面を一回ずつ透過する。

$$
\begin{aligned}
T_0&=q \left(1 - R\right)^{2},\\
T_0&= \left(1\right) \cdot \left(1 - \left(\frac{1}{25}\right)\right)^{2}=\frac{576}{625}.
\end{aligned}
$$

時間的に分離または位相平均された透過を足す。

$$
\begin{aligned}
T&=\frac{T_0}{1 - g},\\
T&= \frac{\left(\frac{576}{625}\right)}{1 - \left(\frac{1}{625}\right)}=\frac{12}{13}.
\end{aligned}
$$

測定値を内部透過だけに割り当てた場合の見かけの値を求める。

$$
\begin{aligned}
\alpha_{\rm app}&=- \frac{\log{\left(T \right)}}{L},\\
\alpha_{\rm app}/(\mathrm{m^{-1}})&= - \frac{\log{\left(\left(\frac{12}{13}\right) \right)}}{\left(\frac{1}{1000}\right)}\simeq 80.043.
\end{aligned}
$$

透過が一未満でも、板の内部で吸収されたことにはならない。
:::

:::example[厚い吸収板]
吸収が厚さにわたって効く場合を評価する。

入力は各面の反射率 $R=\frac{1}{25}$、吸収係数 $\alpha=500\,\mathrm{m^{-1}}$、厚さ $L=\frac{1}{100}\,\mathrm m$ である。

表面での反射を除き、内部のパワー減衰を求める。

$$
\begin{aligned}
q&=e^{- L \alpha},\\
q&= e^{- \left(\frac{1}{100}\right) \cdot \left(500\right)}\simeq 0.0067379.
\end{aligned}
$$

往復時に二回ずつ界面反射と内部通過がある。

$$
\begin{aligned}
g&=R^{2} q^{2},\\
g&\simeq \left(\frac{1}{25}\right)^{2} \cdot \left(0.006737946999\right)^{2}\simeq 7.264\times 10^{-8}.
\end{aligned}
$$

前面と背面を一回ずつ透過する。

$$
\begin{aligned}
T_0&=q \left(1 - R\right)^{2},\\
T_0&\simeq \left(0.006737946999\right) \cdot \left(1 - \left(\frac{1}{25}\right)\right)^{2}\simeq 0.0062097.
\end{aligned}
$$

時間的に分離または位相平均された透過を足す。

$$
\begin{aligned}
T&=\frac{T_0}{1 - g},\\
T&\simeq \frac{\left(0.006209691954\right)}{1 - \left(7.263988762\times 10^{-8}\right)}\simeq 0.0062097.
\end{aligned}
$$

測定値を内部透過だけに割り当てた場合の見かけの値を求める。

$$
\begin{aligned}
\alpha_{\rm app}&=- \frac{\log{\left(T \right)}}{L},\\
\alpha_{\rm app}/(\mathrm{m^{-1}})&\simeq - \frac{\log{\left(\left(0.006209692405\right) \right)}}{\left(\frac{1}{100}\right)}\simeq 508.16.
\end{aligned}
$$

往復透過は小さい。透過指数の補正には主に最初の二面透過を使えばよい。
:::

:::example[反射が大きい薄板]
界面に部分反射面を持つ弱吸収板を考える。

入力は各面の反射率 $R=\frac{2}{5}$、吸収係数 $\alpha=30\,\mathrm{m^{-1}}$、厚さ $L=\frac{1}{500}\,\mathrm m$ である。

表面での反射を除き、内部のパワー減衰を求める。

$$
\begin{aligned}
q&=e^{- L \alpha},\\
q&= e^{- \left(\frac{1}{500}\right) \cdot \left(30\right)}\simeq 0.94176.
\end{aligned}
$$

往復時に二回ずつ界面反射と内部通過がある。

$$
\begin{aligned}
g&=R^{2} q^{2},\\
g&\simeq \left(\frac{2}{5}\right)^{2} \cdot \left(0.9417645336\right)^{2}\simeq 0.14191.
\end{aligned}
$$

前面と背面を一回ずつ透過する。

$$
\begin{aligned}
T_0&=q \left(1 - R\right)^{2},\\
T_0&\simeq \left(0.9417645336\right) \cdot \left(1 - \left(\frac{2}{5}\right)\right)^{2}\simeq 0.33904.
\end{aligned}
$$

時間的に分離または位相平均された透過を足す。

$$
\begin{aligned}
T&=\frac{T_0}{1 - g},\\
T&\simeq \frac{\left(0.3390352321\right)}{1 - \left(0.1419072699\right)}\simeq 0.3951.
\end{aligned}
$$

測定値を内部透過だけに割り当てた場合の見かけの値を求める。

$$
\begin{aligned}
\alpha_{\rm app}&=- \frac{\log{\left(T \right)}}{L},\\
\alpha_{\rm app}/(\mathrm{m^{-1}})&\simeq - \frac{\log{\left(\left(0.3951032565\right) \right)}}{\left(\frac{1}{500}\right)}\simeq 464.3.
\end{aligned}
$$

往復の寄与が見える。最初の透過だけでなく、和を取ってから指数の減衰と分ける。
:::

## 正常分散と異常分散

### 位相の傾き
$n$ が複素屈折率の実部であるとき、$\dd n/\dd\omega>0$ の挙動を :term[正常分散]{en="normal dispersion" yomi="せいじょうぶんさん"}、逆向きを :term[異常分散]{en="anomalous dispersion" yomi="いじょうぶんさん"} と呼ぶ。$\lambda_0=2\pi c/\omega$ で描けば傾きの符号は反転する。振動子の共鳴から離れた弱吸収域で $\chi'\simeq F/(\omega_0^2-\omega^2)$ なので $\dd\chi'/\dd\omega\simeq 2F\omega/(\omega_0^2-\omega^2)^2>0$ である。

共鳴を横切るとこの近似は使えない。$D=\omega_0^2-\omega^2$ とおき、$\chi'=FD/(D^2+\gamma^2\omega^2)$ を微分して
$$
\frac{\dd\chi'}{\dd\omega}=
\frac{2F\omega[(\omega_0^2-\omega^2)^2-\gamma^2\omega_0^2]}{[(\omega_0^2-\omega^2)^2+\gamma^2\omega^2]^2}
$$
を得る。希薄媒質で $n\simeq1+\chi'/2$ とできれば、$|D|<\gamma\omega_0$ の域で $n$ は下がる。強い吸収の域では $n$ の微分に虚部も寄与するので、実誘電率の傾きと同一視できない。

:::example[共鳴の中心]
$\omega=\omega_0$ では $\dd\chi'/\dd\omega=-2F/(\gamma^2\omega_0)<0$ である。減衰を小さくすると下り坂が急になり、同時に共鳴点での吸収も大きくなる。
:::
:::example[分散式の係数]
$\omega\ll\omega_0$ で減衰を無視すれば $\chi'\simeq F/\omega_0^2+F\omega^2/\omega_0^4$。弱応答で $n\simeq1+\chi'/2$ なら、真空波長の逆二乗に比例する最初の補正を得る。
:::

### 透明域で波長の式に直す
各共鳴の減衰を無視できる透明域で、$n^2-1=\sum_j F_j/(\omega_j^2-\omega^2)$ とする。$\lambda_j=2\pi c/\omega_j$、$B_j=F_j/\omega_j^2$ を使うと
$$
n^2=1+\sum_j B_j\frac{\lambda_0^2}{\lambda_0^2-\lambda_j^2}.
$$
これが Lorentz 模型から得られる Sellmeier 形である。実数の係数だけの式は共鳴付近の吸収を表せないため、透明域で使う。$\lambda_0\gg\lambda_j$ で各項を展開すれば、屈折率の波長依存の先頭に $\lambda_0^{-2}$、その次に $\lambda_0^{-4}$ が現れる。

:::example[一共鳴の傾き]
$B_1=1,\ \lambda_0=2\lambda_1$ なら $n^2=1+4/3=7/3$。波長を増やすと第二項は 1 へ下がり、屈折率は $\sqrt2$ に近づく。
:::

### Cauchy 型の分散で角周波数微分を計算する

透明域で真空波長を $l$ とし $n=A+B/l^2+C/l^4$ を用いる。$\dd l/\dd\omega=-l/\omega$ なので $n_g=n-l\,\dd n/\dd l=A+3B/l^2+5C/l^4$。角周波数による波数の二階微分は $\beta_2=-l^2(\dd n_g/\dd l)/(2\pi c^2)$。単色屈折率に対して小さい係数でも、高次の波長依存性は群屈折率や二階微分に大きく現れる。分散の符号を報告するときは $\dd n/\dd l$ の符号と $\beta_2$ の符号を区別する。

:::example[長波長側の透明域]
二項の分散係数が正の材料を長波長側で評価する。

入力は定数項 $A=\frac{29}{20}$、二次の分散係数 $B\simeq 3\times 10^{-15}\,\mathrm{m^2}$、四次の分散係数 $C\simeq 1\times 10^{-28}\,\mathrm{m^4}$、真空波長 $\lambda_0\simeq 8\times 10^{-7}\,\mathrm m$、板厚 $L=\frac{1}{50}\,\mathrm m$ である。

経験式をその波長で評価する。

$$
\begin{aligned}
n&=A + \frac{B}{\lambda_0^{2}} + \frac{C}{\lambda_0^{4}},\\
n&\simeq \left(\frac{29}{20}\right) + \frac{\left(3\times 10^{-15}\right)}{\left(\left(8\times 10^{-7}\right)\right)^{2}} + \frac{\left(1\times 10^{-28}\right)}{\left(\left(8\times 10^{-7}\right)\right)^{4}}\simeq 1.4549.
\end{aligned}
$$

真空波長について微分する。

$$
\begin{aligned}
\dd n/\dd\lambda_0&=- \frac{2 B}{\lambda_0^{3}} - \frac{4 C}{\lambda_0^{5}},\\
\dd n/\dd\lambda_0/(\mathrm{m^{-1}})&\simeq - \frac{2 \cdot \left(3\times 10^{-15}\right)}{\left(\left(8\times 10^{-7}\right)\right)^{3}} - \frac{4 \cdot \left(1\times 10^{-28}\right)}{\left(\left(8\times 10^{-7}\right)\right)^{5}}\simeq -12939.
\end{aligned}
$$

位相だけの項に周波数微分を足す。

$$
\begin{aligned}
n_g&=- \lambda_0 \dd n/\dd\lambda_0 + n,\\
n_g&\simeq - \left(8\times 10^{-7}\right) \cdot \left(-12939.45312\right) + \left(1.454931641\right)\simeq 1.4653.
\end{aligned}
$$

同じ板の二種類の時間を比較する。

$$
\begin{aligned}
\Delta t&=\frac{L \left(n_g - n\right)}{299792458},\\
\Delta t/(\mathrm s)&\simeq \frac{\left(\frac{1}{50}\right) \cdot \left(\left(1.465283203\right) - \left(1.454931641\right)\right)}{299792458}\simeq 6.9058\times 10^{-13}.
\end{aligned}
$$

群屈折率をもう一度周波数で微分する。

$$
\begin{aligned}
\beta_2&=\frac{\frac{3 B}{\lambda_0} + \frac{10 C}{\lambda_0^{3}}}{89875517873681764 \pi},\\
\beta_2/(\mathrm{s^2/m})&\simeq \frac{\frac{3 \cdot \left(3\times 10^{-15}\right)}{\left(8\times 10^{-7}\right)} + \frac{10 \cdot \left(1\times 10^{-28}\right)}{\left(\left(8\times 10^{-7}\right)\right)^{3}}}{89875517873681764 \cdot \pi}\simeq 4.6761\times 10^{-26}.
\end{aligned}
$$

主に $B$ が分散を決める。
:::

:::example[短波長側で高次項を残す]
同じ式を使用できる短波長側で群遅延を比較する。

入力は定数項 $A=\frac{29}{20}$、二次の分散係数 $B\simeq 3\times 10^{-15}\,\mathrm{m^2}$、四次の分散係数 $C\simeq 1\times 10^{-28}\,\mathrm{m^4}$、真空波長 $\lambda_0\simeq 4\times 10^{-7}\,\mathrm m$、板厚 $L=\frac{1}{50}\,\mathrm m$ である。

経験式をその波長で評価する。

$$
\begin{aligned}
n&=A + \frac{B}{\lambda_0^{2}} + \frac{C}{\lambda_0^{4}},\\
n&\simeq \left(\frac{29}{20}\right) + \frac{\left(3\times 10^{-15}\right)}{\left(\left(4\times 10^{-7}\right)\right)^{2}} + \frac{\left(1\times 10^{-28}\right)}{\left(\left(4\times 10^{-7}\right)\right)^{4}}=\frac{377}{256}.
\end{aligned}
$$

真空波長について微分する。

$$
\begin{aligned}
\dd n/\dd\lambda_0&=- \frac{2 B}{\lambda_0^{3}} - \frac{4 C}{\lambda_0^{5}},\\
\dd n/\dd\lambda_0/(\mathrm{m^{-1}})&\simeq - \frac{2 \cdot \left(3\times 10^{-15}\right)}{\left(\left(4\times 10^{-7}\right)\right)^{3}} - \frac{4 \cdot \left(1\times 10^{-28}\right)}{\left(\left(4\times 10^{-7}\right)\right)^{5}}\simeq -1.3281\times 10^{5}.
\end{aligned}
$$

位相だけの項に周波数微分を足す。

$$
\begin{aligned}
n_g&=- \lambda_0 \dd n/\dd\lambda_0 + n,\\
n_g&\simeq - \left(4\times 10^{-7}\right) \cdot \left(-1.328125\times 10^{5}\right) + \left(\frac{377}{256}\right)=\frac{1953}{1280}.
\end{aligned}
$$

同じ板の二種類の時間を比較する。

$$
\begin{aligned}
\Delta t&=\frac{L \left(n_g - n\right)}{299792458},\\
\Delta t/(\mathrm s)&= \frac{\left(\frac{1}{50}\right) \cdot \left(\left(\frac{1953}{1280}\right) - \left(\frac{377}{256}\right)\right)}{299792458}\simeq 3.5441\times 10^{-12}.
\end{aligned}
$$

群屈折率をもう一度周波数で微分する。

$$
\begin{aligned}
\beta_2&=\frac{\frac{3 B}{\lambda_0} + \frac{10 C}{\lambda_0^{3}}}{89875517873681764 \pi},\\
\beta_2/(\mathrm{s^2/m})&\simeq \frac{\frac{3 \cdot \left(3\times 10^{-15}\right)}{\left(4\times 10^{-7}\right)} + \frac{10 \cdot \left(1\times 10^{-28}\right)}{\left(\left(4\times 10^{-7}\right)\right)^{3}}}{89875517873681764 \cdot \pi}\simeq 1.3503\times 10^{-25}.
\end{aligned}
$$

$C$ の項の比が増え、波数の二階微分ではさらに重く掛かる。
:::

:::example[二階微分の零点付近]
透明域の局所近似として二係数が逆符号の式を使う。

入力は定数項 $A=\frac{29}{20}$、二次の分散係数 $B\simeq 3\times 10^{-15}\,\mathrm{m^2}$、四次の分散係数 $C\simeq -3.24\times 10^{-28}\,\mathrm{m^4}$、真空波長 $\lambda_0\simeq 6\times 10^{-7}\,\mathrm m$、板厚 $L=\frac{1}{50}\,\mathrm m$ である。

経験式をその波長で評価する。

$$
\begin{aligned}
n&=A + \frac{B}{\lambda_0^{2}} + \frac{C}{\lambda_0^{4}},\\
n&\simeq \left(\frac{29}{20}\right) + \frac{\left(3\times 10^{-15}\right)}{\left(\left(6\times 10^{-7}\right)\right)^{2}} + \frac{\left(-3.24\times 10^{-28}\right)}{\left(\left(6\times 10^{-7}\right)\right)^{4}}=\frac{1747}{1200}.
\end{aligned}
$$

真空波長について微分する。

$$
\begin{aligned}
\dd n/\dd\lambda_0&=- \frac{2 B}{\lambda_0^{3}} - \frac{4 C}{\lambda_0^{5}},\\
\dd n/\dd\lambda_0/(\mathrm{m^{-1}})&\simeq - \frac{2 \cdot \left(3\times 10^{-15}\right)}{\left(\left(6\times 10^{-7}\right)\right)^{3}} - \frac{4 \cdot \left(-3.24\times 10^{-28}\right)}{\left(\left(6\times 10^{-7}\right)\right)^{5}}\simeq -11111.
\end{aligned}
$$

位相だけの項に周波数微分を足す。

$$
\begin{aligned}
n_g&=- \lambda_0 \dd n/\dd\lambda_0 + n,\\
n_g&\simeq - \left(6\times 10^{-7}\right) \cdot \left(-11111.11111\right) + \left(\frac{1747}{1200}\right)=\frac{117}{80}.
\end{aligned}
$$

同じ板の二種類の時間を比較する。

$$
\begin{aligned}
\Delta t&=\frac{L \left(n_g - n\right)}{299792458},\\
\Delta t/(\mathrm s)&= \frac{\left(\frac{1}{50}\right) \cdot \left(\left(\frac{117}{80}\right) - \left(\frac{1747}{1200}\right)\right)}{299792458}\simeq 4.4475\times 10^{-13}.
\end{aligned}
$$

群屈折率をもう一度周波数で微分する。

$$
\begin{aligned}
\beta_2&=\frac{\frac{3 B}{\lambda_0} + \frac{10 C}{\lambda_0^{3}}}{89875517873681764 \pi},\\
\beta_2/(\mathrm{s^2/m})&\simeq \frac{\frac{3 \cdot \left(3\times 10^{-15}\right)}{\left(6\times 10^{-7}\right)} + \frac{10 \cdot \left(-3.24\times 10^{-28}\right)}{\left(\left(6\times 10^{-7}\right)\right)^{3}}}{89875517873681764 \cdot \pi}=0.
\end{aligned}
$$

この選択では $\beta_2$ がゼロになる。屈折率の一次微分までゼロになる必要はない。
:::

## 金属の光学と Drude モデル

### 自由電荷とプラズマ端
復元力をなくして $\omega_0=0$ とする。これが :term[Drude モデル]{en="Drude model" yomi="どるーでもでる"} である。分極から束縛電子などの背景 $\eps_\infty$ を分けると、
$$
\eps_r(\omega)=\eps_\infty-\frac{\omega_p^2}{\omega(\omega+i\gamma)},\qquad
\sigma(\omega)=\frac{\eps_0\omega_p^2}{\gamma-i\omega},\qquad
\omega_p^2=F.
\label{eq-optics-drude}
$$
導電率は $\vb J=-i\omega\vb P$ から得られ、単位は $\mathrm{S/m}$ である。$\gamma=0$、$\eps_\infty>0$ なら $\eps_r<0$ となる $\omega<\omega_p/\sqrt{\eps_\infty}$ では波数が純虚数になる。半無限媒質の垂直入射では :eqref[eq-optics-metal] から $R=1$。内部の減衰場の存在だけでは発熱を意味しない。有限の $\gamma$ では $\Re\sigma>0$ となって発熱する。

:::example[閾値]
$\eps_\infty=4$、無減衰で $\omega=\omega_p/4$ なら $\eps_r=4-16=-12$。法線波数は $i\sqrt{12}k_0$ である。$\omega=\omega_p$ なら $\eps_r=3>0$ で伝播できる。
:::
:::example[低周波と高周波]
$\omega=\gamma$ では $\sigma=\eps_0\omega_p^2(1+i)/(2\gamma)$ である。実部は直流値の半分で、虚部と等しい。直流値をそのまま光周波数へ延長する近似は $\omega/\gamma\ll1$ に限られる。
:::

### 無損失での境界付近
$\eps_\infty=1,\ \gamma=0$ では $k^2=(\omega^2-\omega_p^2)/c^2$。伝播域で $\omega^2=\omega_p^2+c^2k^2$ を微分すると $v_g=c^2k/\omega$、位相速度は $\omega/k$ で、積は $c^2$ となる。プラズマ端の直上で $v_{\rm ph}$ が大きくなっても $v_g$ はゼロに近づく。端の直下では実の波数を用いたこの群速度の式を適用できない。

### 界面に束縛された金属の波
誘電率 $\eps_d>0$ の誘電体と $\eps_m<0$ の無損失金属を接する。界面に沿う波数を $K$ とし、両側に減衰する TM 波を考える。各側の減衰定数は $\kappa_j=\sqrt{K^2-\omega^2\mu_0\eps_j}>0$。接線磁場を同じ境界値とすると、接線電場の連続は $\kappa_d/\eps_d+\kappa_m/\eps_m=0$ になる。これは磁場の法線微分が両側で逆符号であり、TM の接線電場が $\eps^{-1}$ 倍の法線微分に比例するためである。二乗して整理すると
$$
K^2=\omega^2\mu_0\frac{\eps_d\eps_m}{\eps_d+\eps_m}.
$$
元の符号条件に戻すと $\eps_m<-\eps_d$ で両側の減衰が成立する。このとき $K$ は誘電体の一様平面波の波数より大きい。平面波を単純に界面へ当てても接線波数が足りないので、格子やプリズムで付加の接線成分を与える。有限の金属損失では $K$ も複素数になり、界面沿いにも減衰する。

:::example[無損失 Drude 媒質を端の上で測る]
$\eps_\infty=1$、$\omega=2\omega_p$ なら $k=(\sqrt3/2)\omega/c$ で、群速度は $(\sqrt3/2)c$、位相速度は $(2/\sqrt3)c$。位相速度だけが $c$ を超え、狭いパルスの包絡線の速さは下がっている。
:::

### 導電率の実部と虚部の極限

自由電子の電荷を $-e$、速度を $v$ とする。$m(\partial_t+\gamma)v=-eE$ と $\mathcal J=-Ne v$ をつなぐと $\sigma(\omega)=\eps_0\omega_p^2/(\gamma-i\omega)$。ここで $\omega_p^2=Ne^2/(m\eps_0)$ である。実部は散逸を、虚部は電流の直交位相成分を表す。$\eps_r=1+i\sigma/(\eps_0\omega)$ へ代入すれば Drude 誘電率を得る。$\omega\ll\gamma$ では実導電率がほぼ一定、$\omega\gg\gamma$ では虚部優勢である。後者を直流導電率の式で置き換えると、プラズマ端での反射率の変化を再現できない。

:::example[衝突が速い帯域]
衝突率より低い駆動周波数から見る。

入力はプラズマ角周波数 $\omega_p\simeq 1\times 10^{16}\,\mathrm{s^{-1}}$、衝突率 $\gamma\simeq 1\times 10^{14}\,\mathrm{s^{-1}}$、駆動角周波数 $\omega\simeq 1\times 10^{13}\,\mathrm{s^{-1}}$ である。

電流の応答が静的に近いかを見る。

$$
\begin{aligned}
\omega/\gamma&=\frac{\omega}{\gamma},\\
\omega/\gamma&\simeq \frac{\left(1\times 10^{13}\right)}{\left(1\times 10^{14}\right)}=\frac{1}{10}.
\end{aligned}
$$

自由電子電流の直交位相から実部を得る。

$$
\begin{aligned}
\epsilon_r'&=- \frac{\omega_p^{2}}{\gamma^{2} + \omega^{2}} + 1,\\
\epsilon_r'&\simeq - \frac{\left(\left(1\times 10^{16}\right)\right)^{2}}{\left(\left(1\times 10^{14}\right)\right)^{2} + \left(\left(1\times 10^{13}\right)\right)^{2}} + 1\simeq -9900.
\end{aligned}
$$

散逸する電流から正の虚部を得る。

$$
\begin{aligned}
\epsilon_r''&=\frac{\gamma \omega_p^{2}}{\omega \left(\gamma^{2} + \omega^{2}\right)},\\
\epsilon_r''&\simeq \frac{\left(1\times 10^{14}\right) \cdot \left(\left(1\times 10^{16}\right)\right)^{2}}{\left(1\times 10^{13}\right) \cdot \left(\left(\left(1\times 10^{14}\right)\right)^{2} + \left(\left(1\times 10^{13}\right)\right)^{2}\right)}\simeq 99010.
\end{aligned}
$$

非磁性として、正の虚部を持つ平方根を選ぶ。

$$
\begin{aligned}
\tilde n&=\sqrt{\epsilon_r' + i \epsilon_r''},\\
\tilde n&\simeq \sqrt{\left(-9899.990099\right) + i \cdot \left(99009.90099\right)}\simeq 211.66+233.88i.
\end{aligned}
$$

垂直入射の Fresnel 係数を二乗する。

$$
\begin{aligned}
R&=\left|{\frac{\tilde n - 1}{\tilde n + 1}}\right|^{2},\\
R&\simeq \left|{\frac{\left(211.6643917+233.8841696i\right) - 1}{\left(211.6643917+233.8841696i\right) + 1}}\right|^{2}\simeq 0.99153.
\end{aligned}
$$

誘電率の虚部が実部の大きさを上回る。直流近似はこの領域で始まる。
:::

:::example[プラズマ端の近く]
駆動をプラズマ角周波数へ近づける。

入力はプラズマ角周波数 $\omega_p\simeq 1\times 10^{16}\,\mathrm{s^{-1}}$、衝突率 $\gamma\simeq 1\times 10^{14}\,\mathrm{s^{-1}}$、駆動角周波数 $\omega\simeq 1\times 10^{16}\,\mathrm{s^{-1}}$ である。

電流の応答が静的に近いかを見る。

$$
\begin{aligned}
\omega/\gamma&=\frac{\omega}{\gamma},\\
\omega/\gamma&\simeq \frac{\left(1\times 10^{16}\right)}{\left(1\times 10^{14}\right)}=100.
\end{aligned}
$$

自由電子電流の直交位相から実部を得る。

$$
\begin{aligned}
\epsilon_r'&=- \frac{\omega_p^{2}}{\gamma^{2} + \omega^{2}} + 1,\\
\epsilon_r'&\simeq - \frac{\left(\left(1\times 10^{16}\right)\right)^{2}}{\left(\left(1\times 10^{14}\right)\right)^{2} + \left(\left(1\times 10^{16}\right)\right)^{2}} + 1\simeq 9.999\times 10^{-5}.
\end{aligned}
$$

散逸する電流から正の虚部を得る。

$$
\begin{aligned}
\epsilon_r''&=\frac{\gamma \omega_p^{2}}{\omega \left(\gamma^{2} + \omega^{2}\right)},\\
\epsilon_r''&\simeq \frac{\left(1\times 10^{14}\right) \cdot \left(\left(1\times 10^{16}\right)\right)^{2}}{\left(1\times 10^{16}\right) \cdot \left(\left(\left(1\times 10^{14}\right)\right)^{2} + \left(\left(1\times 10^{16}\right)\right)^{2}\right)}\simeq 0.009999.
\end{aligned}
$$

非磁性として、正の虚部を持つ平方根を選ぶ。

$$
\begin{aligned}
\tilde n&=\sqrt{\epsilon_r' + i \epsilon_r''},\\
\tilde n&\simeq \sqrt{\left(9.9990001\times 10^{-5}\right) + i \cdot \left(0.0099990001\right)}\simeq 0.071062+0.070354i.
\end{aligned}
$$

垂直入射の Fresnel 係数を二乗する。

$$
\begin{aligned}
R&=\left|{\frac{\tilde n - 1}{\tilde n + 1}}\right|^{2},\\
R&\simeq \left|{\frac{\left(0.07106155796+0.07035449537i\right) - 1}{\left(0.07106155796+0.07035449537i\right) + 1}}\right|^{2}\simeq 0.75328.
\end{aligned}
$$

位相屈折率が一より小さく、吸収も残る。反射の変化は急になる。
:::

:::example[プラズマ端より上]
衝突は同じまま駆動をさらに増す。

入力はプラズマ角周波数 $\omega_p\simeq 1\times 10^{16}\,\mathrm{s^{-1}}$、衝突率 $\gamma\simeq 1\times 10^{14}\,\mathrm{s^{-1}}$、駆動角周波数 $\omega\simeq 2\times 10^{16}\,\mathrm{s^{-1}}$ である。

電流の応答が静的に近いかを見る。

$$
\begin{aligned}
\omega/\gamma&=\frac{\omega}{\gamma},\\
\omega/\gamma&\simeq \frac{\left(2\times 10^{16}\right)}{\left(1\times 10^{14}\right)}=200.
\end{aligned}
$$

自由電子電流の直交位相から実部を得る。

$$
\begin{aligned}
\epsilon_r'&=- \frac{\omega_p^{2}}{\gamma^{2} + \omega^{2}} + 1,\\
\epsilon_r'&\simeq - \frac{\left(\left(1\times 10^{16}\right)\right)^{2}}{\left(\left(1\times 10^{14}\right)\right)^{2} + \left(\left(2\times 10^{16}\right)\right)^{2}} + 1\simeq 0.75001.
\end{aligned}
$$

散逸する電流から正の虚部を得る。

$$
\begin{aligned}
\epsilon_r''&=\frac{\gamma \omega_p^{2}}{\omega \left(\gamma^{2} + \omega^{2}\right)},\\
\epsilon_r''&\simeq \frac{\left(1\times 10^{14}\right) \cdot \left(\left(1\times 10^{16}\right)\right)^{2}}{\left(2\times 10^{16}\right) \cdot \left(\left(\left(1\times 10^{14}\right)\right)^{2} + \left(\left(2\times 10^{16}\right)\right)^{2}\right)}\simeq 0.00125.
\end{aligned}
$$

非磁性として、正の虚部を持つ平方根を選ぶ。

$$
\begin{aligned}
\tilde n&=\sqrt{\epsilon_r' + i \epsilon_r''},\\
\tilde n&\simeq \sqrt{\left(0.7500062498\right) + i \cdot \left(0.001249968751\right)}\simeq 0.86603+7.2167\times 10^{-4}i.
\end{aligned}
$$

垂直入射の Fresnel 係数を二乗する。

$$
\begin{aligned}
R&=\left|{\frac{\tilde n - 1}{\tilde n + 1}}\right|^{2},\\
R&\simeq \left|{\frac{\left(0.8660293128+7.216665373\times 10^{-4}i\right) - 1}{\left(0.8660293128+7.216665373\times 10^{-4}i\right) + 1}}\right|^{2}\simeq 0.0051546.
\end{aligned}
$$

実部が正で、虚部が小さい。低周波側での反射を高周波へ延長できない。
:::

## Kramers–Kronig の関係

### 因果律と積分路
実で可積分な応答核 $h(t)$ が $t<0$ でゼロとする。さらに $\chi(z)=\int_0^\infty h(t)e^{izt}\dd t$ が上半平面で正則、閉半円上で一様にゼロへ減衰し、実軸で滑らかとする。振動子のように極が下半平面にある有理関数はこの仮定を満たす。定数の瞬時応答を含む場合は先に引く。

$\chi(z)/(z-\omega)$ を上半平面で積分する。実軸上の $z=\omega$ を上側の微小半円で避けると、半円は時計向きなので寄与は $-i\pi\chi(\omega)$。大半円の寄与はゼロで、囲む極もない。従って $\mathrm{PV}\int_{-\infty}^\infty\chi(\xi)/(\xi-\omega)\dd\xi=i\pi\chi(\omega)$ である。実部と虚部を分けると :term[Kramers–Kronig の関係]{en="Kramers–Kronig relations" yomi="くらまーすくろーにひのかんけい"} が得られる（:eqref[eq-em-kramers-kronig]）。
$$
\chi'(\omega)=\frac1\pi\mathrm{PV}\int_{-\infty}^\infty\frac{\chi''(\xi)}{\xi-\omega}\dd\xi,\qquad
\chi''(\omega)=-\frac1\pi\mathrm{PV}\int_{-\infty}^\infty\frac{\chi'(\xi)}{\xi-\omega}\dd\xi .
\label{eq-optics-kk}
$$
$\mathrm{PV}$ は特異点の両側を同じ幅だけ除く主値である。$h$ が実なので $\chi(-\omega)=\chi(\omega)^*$。$\chi^\prime$ は偶関数、$\chi^{\prime\prime}$ は奇関数である。負の周波数を $\xi\mapsto-\xi$ で折り返すと、第一式には $1/(\xi-\omega)+1/(\xi+\omega)$、第二式には $1/(\xi-\omega)-1/(\xi+\omega)$ が現れる。従って $\omega>0$ では
$$
\chi^\prime(\omega)=\frac{2}{\pi}\mathrm{PV}\int_0^\infty\frac{\xi\chi^{\prime\prime}(\xi)}{\xi^2-\omega^2}\dd\xi,\qquad
\chi^{\prime\prime}(\omega)=-\frac{2\omega}{\pi}\mathrm{PV}\int_0^\infty\frac{\chi^\prime(\xi)}{\xi^2-\omega^2}\dd\xi .
$$
測定が有限帯域なら、その外側にも寄与がある。

:::example[静的感受率]
$\chi''(\xi)\ge0$ の模型では $\chi'(0)=(2/\pi)\int_0^\infty\chi''(\xi)/\xi\,\dd\xi\ge0$。この式を Drude の導電性極へ直接使うと積分が発散するため、零周波数で有限という仮定を見直す必要がある。
:::
:::example[測定帯域より上の吸収]
損失ピークが $\xi>\omega$ にあり他に寄与がなければ、正周波数式の分母が正の範囲から $\chi'(\omega)$ に正の寄与がある。測定波長で吸収が見えなくても屈折率は真空値と異なり得る。
:::

### 高周波の減衰と積分強度
単一 Lorentz 振動子で $\omega>0$ の吸収全体を測る。上半平面で $z\chi(z)$ に極がなく、半円上で $z\chi(z)\sim-F/z$ である。半円を反時計回りに積分すると寄与は $-i\pi F$ に近づく。実軸上の $x\chi(x)$ の実部は奇関数なので対称に積分すると消え、虚部は偶関数である。従って
$$
\int_0^\infty \omega\chi''(\omega)\dd\omega=\frac{\pi F}{2}.
$$
$\gamma$ を変えて吸収の幅と高さが変わっても、この重み付き面積は一定である。高周波で定数となる背景応答を含む場合には、半円を評価する前に背景を引く。

:::example[吸収の積分から振動子強度を戻す]
単一共鳴の重み付き面積 $\int_0^\infty\omega\chi''(\omega)\dd\omega$ が $\pi(10^{28}\,\mathrm{s^{-2}})$ なら $F=2\times10^{28}\,\mathrm{s^{-2}}$。同じ $F$ で減衰率だけを二倍にしたデータでは共鳴ピークは半分になるが、この積分値は保たれる。
:::

### 導電性媒質の零周波数の極
直流導電率 $\sigma_0>0$ を持つとき、分極に自由電流を含めた感受率には $i\sigma_0/(\eps_0\omega)$ が現れる。そこで $\chi_\infty$ とこの極を引いた
$$
\chi_{\rm reg}(\omega)=\chi(\omega)-\chi_\infty-\frac{i\sigma_0}{\eps_0\omega}
$$
を考える。残りが原点で有限かつ上半面で正則であり、前節の減衰条件を満たせば、式 :eqref[eq-optics-kk] を $\chi_{\rm reg}$ に適用できる。実周波数 $\omega\ne0$ で逆に虚部を戻す場合は、積分で得た $\Im\chi_{\rm reg}$ に $\sigma_0/(\eps_0\omega)$ を加える。正周波数の積分から実部を求める場合にも、吸収の虚部から同じ項を先に引いておく。

:::example[Drude の極を分ける]
$\gamma>0$ の Drude 模型で $\sigma_0=\eps_0\omega_p^2/\gamma$ である。部分分数分解すると
$$
-\frac{\omega_p^2}{\omega(\omega+i\gamma)}
=\frac{i\omega_p^2}{\gamma\omega}-\frac{\omega_p^2}{\gamma(\gamma-i\omega)}.
$$
第二項が正則な残りで、その静的実部は $-\omega_p^2/\gamma^2$ である。実際、第二項の虚部を $u$ で割ると $-\omega_p^2/[\gamma(\gamma^2+u^2)]$ となる。これを正周波数で積分し $2/\pi$ 倍しても同じ静的実部になる。
:::

### 単一緩和で因果性と主値積分を確かめる

応答核 $h(t)=(\chi_s/\tau)e^{-t/\tau}\Theta(t)$ を用いると $\chi(\omega)=\int h(t)e^{i\omega t}\dd t=\chi_s/(1-i\omega\tau)$。上半面に極がなく、静的な Kramers–Kronig 関係は
$$
\frac{2}{\pi}\int_0^\infty\frac{\chi''(u)}{u}\dd u
=\frac{2\chi_s}{\pi}\int_0^\infty\frac{\tau\,\dd u}{1+\tau^2u^2}=\chi_s
$$
になる。積分上限を $\Omega$ に切った推定値は真値の $2\arctan(\Omega\tau)/\pi$ 倍。正の吸収しかない受動モデルでは、未測定の高周波分をゼロにするとこの静的値を過小評価する。定数の高周波背景を持つ応答では、関係式を $\chi-\chi_\infty$ に適用する。

:::example[緩和の肩まで測る]
帯域上限が緩和率に一致する測定を考える。

入力は静的感受率 $\chi_s=\frac{1}{2}$、緩和時間 $\tau\simeq 1\times 10^{-12}\,\mathrm s$、測定帯域上限 $\Omega\simeq 1\times 10^{12}\,\mathrm{s^{-1}}$ である。

緩和の速さと帯域上限を同じ無次元量にする。

$$
\begin{aligned}
x&=\tau \Omega,\\
x&\simeq \left(1\times 10^{-12}\right) \cdot \left(1\times 10^{12}\right)=1.
\end{aligned}
$$

単一緩和の同相成分を計算する。

$$
\begin{aligned}
\chi'(\Omega)&=\frac{\chi_s}{x^{2} + 1},\\
\chi'(\Omega)&= \frac{\left(\frac{1}{2}\right)}{\left(1\right)^{2} + 1}=\frac{1}{4}.
\end{aligned}
$$

積分の被積分関数に入る吸収成分も求める。

$$
\begin{aligned}
\chi''(\Omega)&=\frac{\chi_s x}{x^{2} + 1},\\
\chi''(\Omega)&= \frac{\left(\frac{1}{2}\right) \cdot \left(1\right)}{\left(1\right)^{2} + 1}=\frac{1}{4}.
\end{aligned}
$$

帯域の正負を対称に含めて積分の割合を求める。

$$
\begin{aligned}
f&=\frac{2 \operatorname{atan}{\left(x \right)}}{\pi},\\
f&= \frac{2 \cdot \operatorname{atan}{\left(\left(1\right) \right)}}{\pi}=\frac{1}{2}.
\end{aligned}
$$

全周波数の積分との差を取る。

$$
\begin{aligned}
\Delta\chi_s&=\chi_s \left(1 - f\right),\\
\Delta\chi_s&= \left(\frac{1}{2}\right) \cdot \left(1 - \left(\frac{1}{2}\right)\right)=\frac{1}{4}.
\end{aligned}
$$

虚部はそこで最大になるが、静的感受率の積分の半分は外に残る。
:::

:::example[帯域を十倍にする]
同じ試料で測定端を上げる。

入力は静的感受率 $\chi_s=\frac{1}{2}$、緩和時間 $\tau\simeq 1\times 10^{-12}\,\mathrm s$、測定帯域上限 $\Omega\simeq 1\times 10^{13}\,\mathrm{s^{-1}}$ である。

緩和の速さと帯域上限を同じ無次元量にする。

$$
\begin{aligned}
x&=\tau \Omega,\\
x&\simeq \left(1\times 10^{-12}\right) \cdot \left(1\times 10^{13}\right)=10.
\end{aligned}
$$

単一緩和の同相成分を計算する。

$$
\begin{aligned}
\chi'(\Omega)&=\frac{\chi_s}{x^{2} + 1},\\
\chi'(\Omega)&= \frac{\left(\frac{1}{2}\right)}{\left(10\right)^{2} + 1}=\frac{1}{202}.
\end{aligned}
$$

積分の被積分関数に入る吸収成分も求める。

$$
\begin{aligned}
\chi''(\Omega)&=\frac{\chi_s x}{x^{2} + 1},\\
\chi''(\Omega)&= \frac{\left(\frac{1}{2}\right) \cdot \left(10\right)}{\left(10\right)^{2} + 1}=\frac{5}{101}.
\end{aligned}
$$

帯域の正負を対称に含めて積分の割合を求める。

$$
\begin{aligned}
f&=\frac{2 \operatorname{atan}{\left(x \right)}}{\pi},\\
f&= \frac{2 \cdot \operatorname{atan}{\left(\left(10\right) \right)}}{\pi}\simeq 0.93655.
\end{aligned}
$$

全周波数の積分との差を取る。

$$
\begin{aligned}
\Delta\chi_s&=\chi_s \left(1 - f\right),\\
\Delta\chi_s&\simeq \left(\frac{1}{2}\right) \cdot \left(1 - \left(0.9365489651\right)\right)\simeq 0.031726.
\end{aligned}
$$

積分の大半を回収できる。虚部の裾が積分に残る点に注意する。
:::

:::example[緩和が速い試料]
同じ測定端で緩和時間の短い試料へ替える。

入力は静的感受率 $\chi_s=\frac{1}{2}$、緩和時間 $\tau\simeq 1\times 10^{-14}\,\mathrm s$、測定帯域上限 $\Omega\simeq 1\times 10^{13}\,\mathrm{s^{-1}}$ である。

緩和の速さと帯域上限を同じ無次元量にする。

$$
\begin{aligned}
x&=\tau \Omega,\\
x&\simeq \left(1\times 10^{-14}\right) \cdot \left(1\times 10^{13}\right)=\frac{1}{10}.
\end{aligned}
$$

単一緩和の同相成分を計算する。

$$
\begin{aligned}
\chi'(\Omega)&=\frac{\chi_s}{x^{2} + 1},\\
\chi'(\Omega)&= \frac{\left(\frac{1}{2}\right)}{\left(\frac{1}{10}\right)^{2} + 1}=\frac{50}{101}.
\end{aligned}
$$

積分の被積分関数に入る吸収成分も求める。

$$
\begin{aligned}
\chi''(\Omega)&=\frac{\chi_s x}{x^{2} + 1},\\
\chi''(\Omega)&= \frac{\left(\frac{1}{2}\right) \cdot \left(\frac{1}{10}\right)}{\left(\frac{1}{10}\right)^{2} + 1}=\frac{5}{101}.
\end{aligned}
$$

帯域の正負を対称に含めて積分の割合を求める。

$$
\begin{aligned}
f&=\frac{2 \operatorname{atan}{\left(x \right)}}{\pi},\\
f&= \frac{2 \cdot \operatorname{atan}{\left(\left(\frac{1}{10}\right) \right)}}{\pi}\simeq 0.063451.
\end{aligned}
$$

全周波数の積分との差を取る。

$$
\begin{aligned}
\Delta\chi_s&=\chi_s \left(1 - f\right),\\
\Delta\chi_s&\simeq \left(\frac{1}{2}\right) \cdot \left(1 - \left(0.06345103486\right)\right)\simeq 0.46827.
\end{aligned}
$$

測定は緩和の中心まで届かず、静的感受率の推定に大きな補正が必要になる。
:::

## 光パルスの伝播と群速度分散

パルスの中心を $\omega_c$、離調を $\Omega$、実伝播定数を $\beta(\omega)$ とする。帯域が狭く吸収を無視できる媒質で、包絡線は各周波数に $e^{i[\beta(\omega_c+\Omega)-\beta(\omega_c)]z}$ を掛けて求める。$\beta_j=\beta^{(j)}(\omega_c)$ を伝播定数の第 $j$ 微分の中心での値とする。$\beta_1 z$ を遅延として引いた時間 $T=t-\beta_1z$ で二次まで残せば
$$
A(z,T)=\int a(\Omega)e^{i\beta_2 z\Omega^2/2-i\Omega T}\frac{\dd\Omega}{2\pi},\qquad
\partial_z A=-\frac{i\beta_2}{2}\partial_T^2A .
\label{eq-optics-gvd}
$$
$\beta_2$ が :term[群速度分散]{en="group velocity dispersion" yomi="ぐんそくどぶんさん"} で、単位は $\mathrm{s^2/m}$。群速度は $1/\beta_1$、群速度分散の符号は $\dd n/\dd\omega$ の符号だけでは決まらない。

### ガウスパルスを積分する
$A(0,T)=A_0 e^{-T^2/(2T_0^2)}$ のとき $a$ は $e^{-T_0^2\Omega^2/2}$ に比例する。Gauss 積分で $T_0^2$ を $T_0^2-i\beta_2z$ に置き換えると
$$
A(z,T)=\frac{A_0}{\sqrt{1-i b}}\exp\left[-\frac{T^2}{2T_0^2(1-i b)}\right],\qquad b=\beta_2z/T_0^2 .
$$
平方根は $z=0$ で正で連続な枝を取る。$|A|^2$ の $1/e$ 半幅は $T_0\sqrt{1+b^2}$ となり、エネルギー積分は一定である。位相の $T^2$ の項から瞬時周波数偏移 $-\partial_T\arg A$ は $bT/[T_0^2(1+b^2)]$ になる。$\beta_2>0$ の場合は高い周波数が遅れて現れる。

:::example[分散長]
$T_0=100\,\mathrm{fs}$、$\beta_2=20\,\mathrm{fs^2/mm}$ なら $T_0^2/|\beta_2|=500\,\mathrm{mm}$。この距離で幅は $\sqrt2$ 倍となる。
:::
:::example[帯域の変更]
パルス幅を半分にすると分散長は四分の一になる。さらに $|\beta_3z|/T_0^3$ も増えるので、極端に短いパルスに二次近似を使う際は三次位相を評価する。
:::

### 予め二次位相があるパルス
二次分散媒質へ入れる前からスペクトルに $e^{i\Phi_2\Omega^2/2}$ を掛けたとする。一次遅延を引くと出力の二次位相係数は $\Phi_2+\beta_2L$ になる。ガウススペクトルの絶対値が同じなら、幅は $T_0\sqrt{1+[(\Phi_2+\beta_2L)/T_0^2]^2}$ である。したがって $\Phi_2\beta_2<0$ なら適当な伝播長でパルスを再圧縮できる。強度のピークを増す作用は、エネルギーの増幅とは区別する。

:::example[逆符号の分散]
$T_0=50\,\mathrm{fs}$、$\Phi_2=-1000\,\mathrm{fs^2}$、$\beta_2=25\,\mathrm{fs^2/mm}$ なら $L=40\,\mathrm{mm}$ で二次位相が消える。入射の幅は約 $53.9\,\mathrm{fs}$、出射の幅は $50\,\mathrm{fs}$ である。
:::

### 波長で測った群遅延
$\beta=n(\omega)\omega/c$ なので $\beta_1=(n+\omega n')/c$、$\beta_2=(2n'+\omega n'')/c$。真空波長を独立変数とするなら $n_g=n-\lambda_0\dd n/\dd\lambda_0$。もう一度波長で微分すると $\dd n_g/\dd\lambda_0=-\lambda_0\dd^2n/\dd\lambda_0^2$。連鎖律 $\dd\lambda_0/\dd\omega=-\lambda_0^2/(2\pi c)$ と組み合わせて
$$
\beta_2=\frac{\lambda_0^3}{2\pi c^2}\frac{\dd^2 n}{\dd\lambda_0^2}.
$$
従って正常分散 $\dd n/\dd\omega>0$ と $\beta_2>0$ は定義から異なる。材料の表が波長の屈折率を与えるなら、群遅延は一階微分、群速度分散は二階微分を取って求める。

### 二次位相が時間帯域幅積を変える

真空中での搬送周波数からのずれを $\Omega$、包絡線のスペクトルを $\exp(-T_0^2\Omega^2/2)$ とする。二次位相 $\exp(iB\Omega^2/2)$ を加えた Fourier 積分は $\exp[-t^2/(2(T_0^2-iB))]/\sqrt{T_0^2-iB}$ に比例する。平方完成は $T_0^2>0$ なら複素係数でも成立する。強度の $1/e$ 半幅は $T_0\sqrt{1+(B/T_0^2)^2}$。強度スペクトルは二次位相では変わらないので、半値全幅の積は $\Delta t_{\rm FWHM}\Delta\nu_{\rm FWHM}=(2\ln 2/\pi)\sqrt{1+(B/T_0^2)^2}$ になる。実数の $B$ の符号は時間の強度だけでは分からない。スペクトル位相または追加の分散を使って区別する。

:::example[分散を加えてパルスを伸ばす]
位相のないガウススペクトルを分散媒質に入れる。

入力は無チャープ強度の半幅 $T_0\simeq 1\times 10^{-13}\,\mathrm s$、入射時の二次位相係数 $B_{\rm in}=0\,\mathrm{s^2}$、媒質の波数の二階微分 $\beta_2\simeq 2\times 10^{-26}\,\mathrm{s^2/m}$、伝播距離 $L=\frac{1}{2}\,\mathrm m$ である。

媒質の二次位相を入射位相に足す。

$$
\begin{aligned}
B_{\rm out}&=B_{\rm in} + L \beta_2,\\
B_{\rm out}/(\mathrm{s^2})&\simeq \left(0\right) + \left(\frac{1}{2}\right) \cdot \left(2\times 10^{-26}\right)\simeq 1\times 10^{-26}.
\end{aligned}
$$

複素係数の逆数の実部から時間幅を求める。

$$
\begin{aligned}
g&=\sqrt{\frac{B_{\rm out}^{2}}{T_0^{4}} + 1},\\
g&\simeq \sqrt{\frac{\left(\left(1\times 10^{-26}\right)\right)^{2}}{\left(\left(1\times 10^{-13}\right)\right)^{4}} + 1}\simeq 1.4142.
\end{aligned}
$$

ガウスの強度が半分になる二点を求める。

$$
\begin{aligned}
\Delta t&=2 T_0 g \sqrt{\log{\left(2 \right)}},\\
\Delta t/(\mathrm s)&\simeq 2 \cdot \left(1\times 10^{-13}\right) \cdot \left(1.414213562\right) \cdot \sqrt{\log{\left(2 \right)}}\simeq 2.3548\times 10^{-13}.
\end{aligned}
$$

角周波数の幅を $2\pi$ で割る。

$$
\begin{aligned}
\Delta\nu&=\frac{\sqrt{\log{\left(2 \right)}}}{\pi T_0},\\
\Delta\nu/(\mathrm{Hz})&\simeq \frac{\sqrt{\log{\left(2 \right)}}}{\pi \cdot \left(1\times 10^{-13}\right)}\simeq 2.6501\times 10^{12}.
\end{aligned}
$$

強度の二つの半値全幅の積を取る。

$$
\begin{aligned}
P&=\Delta t \Delta\nu,\\
P&\simeq \left(2.354820045\times 10^{-13}\right) \cdot \left(2.650103635\times 10^{12}\right)\simeq 0.62405.
\end{aligned}
$$

二次位相の無次元値は一で、時間幅が $\sqrt2$ 倍になる。
:::

:::example[逆分散で補償する]
先に与えた二次位相を媒質で打ち消す。

入力は無チャープ強度の半幅 $T_0\simeq 1\times 10^{-13}\,\mathrm s$、入射時の二次位相係数 $B_{\rm in}\simeq -1\times 10^{-26}\,\mathrm{s^2}$、媒質の波数の二階微分 $\beta_2\simeq 2\times 10^{-26}\,\mathrm{s^2/m}$、伝播距離 $L=\frac{1}{2}\,\mathrm m$ である。

媒質の二次位相を入射位相に足す。

$$
\begin{aligned}
B_{\rm out}&=B_{\rm in} + L \beta_2,\\
B_{\rm out}/(\mathrm{s^2})&\simeq \left(-1\times 10^{-26}\right) + \left(\frac{1}{2}\right) \cdot \left(2\times 10^{-26}\right)=0.
\end{aligned}
$$

複素係数の逆数の実部から時間幅を求める。

$$
\begin{aligned}
g&=\sqrt{\frac{B_{\rm out}^{2}}{T_0^{4}} + 1},\\
g&\simeq \sqrt{\frac{\left(0\right)^{2}}{\left(\left(1\times 10^{-13}\right)\right)^{4}} + 1}=1.
\end{aligned}
$$

ガウスの強度が半分になる二点を求める。

$$
\begin{aligned}
\Delta t&=2 T_0 g \sqrt{\log{\left(2 \right)}},\\
\Delta t/(\mathrm s)&\simeq 2 \cdot \left(1\times 10^{-13}\right) \cdot \left(1\right) \cdot \sqrt{\log{\left(2 \right)}}\simeq 1.6651\times 10^{-13}.
\end{aligned}
$$

角周波数の幅を $2\pi$ で割る。

$$
\begin{aligned}
\Delta\nu&=\frac{\sqrt{\log{\left(2 \right)}}}{\pi T_0},\\
\Delta\nu/(\mathrm{Hz})&\simeq \frac{\sqrt{\log{\left(2 \right)}}}{\pi \cdot \left(1\times 10^{-13}\right)}\simeq 2.6501\times 10^{12}.
\end{aligned}
$$

強度の二つの半値全幅の積を取る。

$$
\begin{aligned}
P&=\Delta t \Delta\nu,\\
P&\simeq \left(1.665109222\times 10^{-13}\right) \cdot \left(2.650103635\times 10^{12}\right)\simeq 0.44127.
\end{aligned}
$$

時間帯域幅積はもとの最小値へ戻る。
:::

:::example[過補償した時間波形]
入射時の二次位相の大きさを増し、補償を越えた符号の例を作る。

入力は無チャープ強度の半幅 $T_0\simeq 1\times 10^{-13}\,\mathrm s$、入射時の二次位相係数 $B_{\rm in}\simeq -2\times 10^{-26}\,\mathrm{s^2}$、媒質の波数の二階微分 $\beta_2\simeq 2\times 10^{-26}\,\mathrm{s^2/m}$、伝播距離 $L=\frac{1}{2}\,\mathrm m$ である。

媒質の二次位相を入射位相に足す。

$$
\begin{aligned}
B_{\rm out}&=B_{\rm in} + L \beta_2,\\
B_{\rm out}/(\mathrm{s^2})&\simeq \left(-2\times 10^{-26}\right) + \left(\frac{1}{2}\right) \cdot \left(2\times 10^{-26}\right)\simeq -1\times 10^{-26}.
\end{aligned}
$$

複素係数の逆数の実部から時間幅を求める。

$$
\begin{aligned}
g&=\sqrt{\frac{B_{\rm out}^{2}}{T_0^{4}} + 1},\\
g&\simeq \sqrt{\frac{\left(\left(-1\times 10^{-26}\right)\right)^{2}}{\left(\left(1\times 10^{-13}\right)\right)^{4}} + 1}\simeq 1.4142.
\end{aligned}
$$

ガウスの強度が半分になる二点を求める。

$$
\begin{aligned}
\Delta t&=2 T_0 g \sqrt{\log{\left(2 \right)}},\\
\Delta t/(\mathrm s)&\simeq 2 \cdot \left(1\times 10^{-13}\right) \cdot \left(1.414213562\right) \cdot \sqrt{\log{\left(2 \right)}}\simeq 2.3548\times 10^{-13}.
\end{aligned}
$$

角周波数の幅を $2\pi$ で割る。

$$
\begin{aligned}
\Delta\nu&=\frac{\sqrt{\log{\left(2 \right)}}}{\pi T_0},\\
\Delta\nu/(\mathrm{Hz})&\simeq \frac{\sqrt{\log{\left(2 \right)}}}{\pi \cdot \left(1\times 10^{-13}\right)}\simeq 2.6501\times 10^{12}.
\end{aligned}
$$

強度の二つの半値全幅の積を取る。

$$
\begin{aligned}
P&=\Delta t \Delta\nu,\\
P&\simeq \left(2.354820045\times 10^{-13}\right) \cdot \left(2.650103635\times 10^{12}\right)\simeq 0.62405.
\end{aligned}
$$

最初の例と同じ時間幅だが、位相の符号は逆である。
:::

## Rayleigh 散乱と Mie 散乱

### 双極子が支配する小球
真空中の半径 $a$ の非磁性球で、相対誘電率を $\eps_r$ とする。$k_0a\ll1$ かつ $|\sqrt{\eps_r}|k_0a\ll1$ とし、静電分極の共鳴から離れているとする。場が球内でほぼ一様なら静電ポテンシャルで局所分極を求められる。入射場を $E_0\hat{\vb z}$ とし外部を $-E_0r\cos\theta+B\cos\theta/r^2$、内部を $-C r\cos\theta$ と置く。$r=a$ でポテンシャルと法線電束を連続にすると $C=E_0-B/a^3$ と $\eps_r C=E_0+2B/a^3$。従って双極子は $p=4\pi\eps_0 B$ で
$$
p=4\pi\eps_0a^3\frac{\eps_r-1}{\eps_r+2}E_0.
$$
遠方では双極子の電場振幅は $k_0^2p\sin\theta/(4\pi\eps_0 r)$ となる（:eqref[eq-em-dipole-radiation]）。照度比に $r^2$ を掛けると微分散乱断面積となる。$\int\sin^2\theta\,\dd\Omega=8\pi/3$ より、:term[Rayleigh 散乱]{en="Rayleigh scattering" yomi="れいりーさんらん"} の全断面積は
$$
\sigma_{\mathrm{sca}}=\frac{8\pi}{3}k_0^4a^6\left|\frac{\eps_r-1}{\eps_r+2}\right|^2 .
\label{eq-optics-rayleigh-scatter}
$$
吸収が小さく材料分散を無視できる範囲では $a^6/\lambda_0^4$ に比例する。$\eps_r\simeq-2$ では静電近似の分母が小さく、散逸や放射反作用の補正が必要である。

### 有限の大きさの球
粒径が波長の程度になると、球の中の位相が変化して双極子一つでは足りない。:term[Mie 散乱]{en="Mie scattering" yomi="みーさんらん"} では、各角運動量 $\ell\ge1$ の電気型と磁気型の球面波を使う。球の外側では入射に正則な球 Bessel 関数 $j_\ell$、散乱に外向きの球 Hankel 関数 $h_\ell^{(1)}$、内側には原点で正則な $j_\ell$ を取る。$x=k_0a$、$m=\sqrt{\eps_r}$、$\psi_\ell(u)=u j_\ell(u)$、$\xi_\ell(u)=u h_\ell^{(1)}(u)$ とすると、非磁性球での標準的な散乱係数は
$$
\begin{aligned}
a_\ell&=\frac{m\psi_\ell(mx)\psi_\ell'(x)-\psi_\ell(x)\psi_\ell'(mx)}
{m\psi_\ell(mx)\xi_\ell'(x)-\xi_\ell(x)\psi_\ell'(mx)},\\
b_\ell&=\frac{\psi_\ell(mx)\psi_\ell'(x)-m\psi_\ell(x)\psi_\ell'(mx)}
{\psi_\ell(mx)\xi_\ell'(x)-m\xi_\ell(x)\psi_\ell'(mx)} .
\end{aligned}
$$
電気型では $m\psi_\ell(mx)(\psi_\ell'-a_\ell\xi_\ell')=\psi_\ell'(mx)(\psi_\ell-a_\ell\xi_\ell)$、磁気型では $m$ の位置を入れ替える。これらは接線場の二条件から内側の振幅を消去した式である（:eqref[eq-em-mie-channel]）。直交性で異なる成分の積分が消え、全散乱断面積は $(2\pi/k_0^2)\sum_{\ell=1}^\infty(2\ell+1)(|a_\ell|^2+|b_\ell|^2)$ となる。入射平面波の展開とその次数ごとの重みの導出は :eqref[eq-em-mie-cross-sections] にある。$x,|m|x\ll1$ で $\psi_1(x)\sim x^2/3,\ \xi_1(x)\sim-i/x$ を入れると $a_1\sim-2i x^3(m^2-1)/[3(m^2+2)]$ となり、断面積は式 :eqref[eq-optics-rayleigh-scatter] に一致する。

:::example[波長の比較]
無分散の同じ小球では $400\,\mathrm{nm}$ の散乱断面積は $800\,\mathrm{nm}$ の 16 倍である。粒径が大きくなって小球条件が破れると、この比は使えない。
:::
:::example[コントラストのない球]
$m=1$ のとき両係数の分子は各 $\ell$ でゼロとなる。球の半径が大きくても、内外の物性が同じなら散乱しない。
:::

### 吸収と散乱を分けた尺度
静電近似の双極子係数を $\alpha_d=p/E_0$ とする。双極子近似で粒子が入射場から取り出す平均パワーは $\langle \dot p E\rangle=(\omega/2)\Im\alpha_d |E_0|^2$ である。照度 $|E_0|^2/(2\eta_0)$ で割ると消衰断面積 $\sigma_{\rm ext}=k_0\Im\alpha_d/\eps_0$ を得る。消衰は散乱と内部吸収の和である。無損失の静電近似では $\alpha_d$ が実となって消衰ゼロだが、有限波数の散乱式は正であり、異なる次数の近似を混ぜたために保存則が合わない。実の静電係数 $\alpha_{\rm es}$ に対し $\alpha_d=\alpha_{\rm es}/[1-i k_0^3\alpha_{\rm es}/(6\pi\eps_0)]$ と放射反作用の項を入れると、$\Im\alpha_d=k_0^3|\alpha_d|^2/(6\pi\eps_0)$ となり、双極子の散乱断面積と消衰断面積が等しくなる。ここで補正したのは放射の虚部であり、有限サイズに由来する実部の補正まで求めたわけではない。球についての散乱と消衰の区別は :eqref[eq-em-mie-cross-sections] も参照する。

:::example[体積と断面積の比]
固定の小球材質と波長で半径を半分にすると双極子モーメントは $1/8$、散乱は $1/64$ になる。内部吸収が支配的で静電係数の虚部を使えるなら、吸収の最初の項は体積に比例して $1/8$ になる。
:::

### Rayleigh 領域で透過と散乱断面積を分ける

希薄な同一粒子の数密度を $N$ とする。独立散乱で光のうち検出器の直接光から失われる分を $-\dd I=N\sigma_{\rm ext}I\dd z$ と書くと、光学的厚さは $NL\sigma_{\rm ext}$。無吸収の小球では $\sigma_{\rm ext}=\sigma_{\rm sca}$ である。透明な非磁性背景に置く場合も球内外の静電境界条件は同じ形となり、$\eps_r$ を粒子と背景の誘電率比 $m^2$、$k_0$ を背景の波数 $2\pi/\lambda$ に置き換えて式 :eqref[eq-optics-rayleigh-scatter] を使う。遠方の散乱照度と入射照度で同じ背景のインピーダンスが現れ、比を取ると打ち消される。無吸収なのに透過が減るのは散乱方向へエネルギーが移るためである。有限の開口の検出器が散乱光の一部を拾う場合、直接光の測定と全パワーの測定を区別する。濃度を増した媒体で多重散乱を無視できないときは、散乱後に別の粒子で進行方向が変わる過程も必要になる。

:::example[希薄な試料を長波長で見る]
無吸収の微粒子を小さなサイズパラメータで評価する。

入力は粒子半径 $a\simeq 1\times 10^{-8}\,\mathrm m$、背景媒質中での波長 $\lambda\simeq 6\times 10^{-7}\,\mathrm m$、背景に対する粒子屈折率比 $m=\frac{3}{2}$、数密度 $N\simeq 1\times 10^{18}\,\mathrm{m^{-3}}$、光路長 $L=\frac{1}{100}\,\mathrm m$ である。

小球近似の適用条件を調べる。

$$
\begin{aligned}
x&=\frac{2 \pi a}{\lambda},\\
x&\simeq \frac{2 \cdot \pi \cdot \left(1\times 10^{-8}\right)}{\left(6\times 10^{-7}\right)}\simeq 0.10472.
\end{aligned}
$$

静電境界条件に由来する球の分極率の係数を入れる。

$$
\begin{aligned}
f&=\frac{m^{2} - 1}{m^{2} + 2},\\
f&= \frac{\left(\frac{3}{2}\right)^{2} - 1}{\left(\frac{3}{2}\right)^{2} + 2}=\frac{5}{17}.
\end{aligned}
$$

幾何面積に無次元の散乱効率を掛ける。

$$
\begin{aligned}
\sigma_{\rm sca}&=\frac{8 \pi a^{2} f^{2} x^{4}}{3},\\
\sigma_{\rm sca}/(\mathrm{m^2})&\simeq \frac{8 \cdot \pi \cdot \left(\left(1\times 10^{-8}\right)\right)^{2} \cdot \left(\frac{5}{17}\right)^{2} \cdot \left(0.1047197551\right)^{4}}{3}\simeq 8.7152\times 10^{-21}.
\end{aligned}
$$

単位長さ当たりの消散を粒子数で数える。

$$
\begin{aligned}
\tau&=L N \sigma_{\rm sca},\\
\tau&\simeq \left(\frac{1}{100}\right) \cdot \left(1\times 10^{18}\right) \cdot \left(8.715157554\times 10^{-21}\right)\simeq 8.7152\times 10^{-5}.
\end{aligned}
$$

一次方程式を積分し、前方の直接成分を求める。

$$
\begin{aligned}
T&=e^{- \tau},\\
T&\simeq e^{- \left(8.715157554\times 10^{-5}\right)}\simeq 0.99991.
\end{aligned}
$$

直接光の減少が小さく、光学的厚さをその減少率の近似に使える。
:::

:::example[波長を半分にする]
材質の分散を無視して波長の効果だけを測る。

入力は粒子半径 $a\simeq 1\times 10^{-8}\,\mathrm m$、背景媒質中での波長 $\lambda\simeq 3\times 10^{-7}\,\mathrm m$、背景に対する粒子屈折率比 $m=\frac{3}{2}$、数密度 $N\simeq 1\times 10^{18}\,\mathrm{m^{-3}}$、光路長 $L=\frac{1}{100}\,\mathrm m$ である。

小球近似の適用条件を調べる。

$$
\begin{aligned}
x&=\frac{2 \pi a}{\lambda},\\
x&\simeq \frac{2 \cdot \pi \cdot \left(1\times 10^{-8}\right)}{\left(3\times 10^{-7}\right)}\simeq 0.20944.
\end{aligned}
$$

静電境界条件に由来する球の分極率の係数を入れる。

$$
\begin{aligned}
f&=\frac{m^{2} - 1}{m^{2} + 2},\\
f&= \frac{\left(\frac{3}{2}\right)^{2} - 1}{\left(\frac{3}{2}\right)^{2} + 2}=\frac{5}{17}.
\end{aligned}
$$

幾何面積に無次元の散乱効率を掛ける。

$$
\begin{aligned}
\sigma_{\rm sca}&=\frac{8 \pi a^{2} f^{2} x^{4}}{3},\\
\sigma_{\rm sca}/(\mathrm{m^2})&\simeq \frac{8 \cdot \pi \cdot \left(\left(1\times 10^{-8}\right)\right)^{2} \cdot \left(\frac{5}{17}\right)^{2} \cdot \left(0.2094395102\right)^{4}}{3}\simeq 1.3944\times 10^{-19}.
\end{aligned}
$$

単位長さ当たりの消散を粒子数で数える。

$$
\begin{aligned}
\tau&=L N \sigma_{\rm sca},\\
\tau&\simeq \left(\frac{1}{100}\right) \cdot \left(1\times 10^{18}\right) \cdot \left(1.394425209\times 10^{-19}\right)\simeq 0.0013944.
\end{aligned}
$$

一次方程式を積分し、前方の直接成分を求める。

$$
\begin{aligned}
T&=e^{- \tau},\\
T&\simeq e^{- \left(0.001394425209\right)}\simeq 0.99861.
\end{aligned}
$$

断面積は十六倍になる。短波長ほど散乱が増える条件をこの模型で切り出せる。
:::

:::example[屈折率を背景に近づける]
短波長のまま粒子の屈折率比を一へ近づける。

入力は粒子半径 $a\simeq 1\times 10^{-8}\,\mathrm m$、背景媒質中での波長 $\lambda\simeq 3\times 10^{-7}\,\mathrm m$、背景に対する粒子屈折率比 $m=\frac{21}{20}$、数密度 $N\simeq 1\times 10^{18}\,\mathrm{m^{-3}}$、光路長 $L=\frac{1}{100}\,\mathrm m$ である。

小球近似の適用条件を調べる。

$$
\begin{aligned}
x&=\frac{2 \pi a}{\lambda},\\
x&\simeq \frac{2 \cdot \pi \cdot \left(1\times 10^{-8}\right)}{\left(3\times 10^{-7}\right)}\simeq 0.20944.
\end{aligned}
$$

静電境界条件に由来する球の分極率の係数を入れる。

$$
\begin{aligned}
f&=\frac{m^{2} - 1}{m^{2} + 2},\\
f&= \frac{\left(\frac{21}{20}\right)^{2} - 1}{\left(\frac{21}{20}\right)^{2} + 2}=\frac{41}{1241}.
\end{aligned}
$$

幾何面積に無次元の散乱効率を掛ける。

$$
\begin{aligned}
\sigma_{\rm sca}&=\frac{8 \pi a^{2} f^{2} x^{4}}{3},\\
\sigma_{\rm sca}/(\mathrm{m^2})&\simeq \frac{8 \cdot \pi \cdot \left(\left(1\times 10^{-8}\right)\right)^{2} \cdot \left(\frac{41}{1241}\right)^{2} \cdot \left(0.2094395102\right)^{4}}{3}\simeq 1.7595\times 10^{-21}.
\end{aligned}
$$

単位長さ当たりの消散を粒子数で数える。

$$
\begin{aligned}
\tau&=L N \sigma_{\rm sca},\\
\tau&\simeq \left(\frac{1}{100}\right) \cdot \left(1\times 10^{18}\right) \cdot \left(1.759451136\times 10^{-21}\right)\simeq 1.7595\times 10^{-5}.
\end{aligned}
$$

一次方程式を積分し、前方の直接成分を求める。

$$
\begin{aligned}
T&=e^{- \tau},\\
T&\simeq e^{- \left(1.759451136\times 10^{-5}\right)}\simeq 0.99998.
\end{aligned}
$$

波長が短くても、分極コントラストを減らせば散乱を抑えられる。
:::

## 演習問題

::::exercise{level=basic}
$\omega=\omega_0$ での Lorentz 感受率を求めよ。

:::solution
分母は $-i\gamma\omega_0$ となるので $\chi=iF/(\gamma\omega_0)$ である。虚部が正なので仕事率も正になる。
:::
::::

::::exercise{level=basic}
$\lambda_0=600\,\mathrm{nm}$、$\kappa=10^{-4}$ での吸収係数を求めよ。

:::solution
$\alpha=4\pi(10^{-4})/(6\times10^{-7})=(2000\pi/3)\,\mathrm{m^{-1}}$ である。
:::
::::

::::exercise{level=basic}
強度が厚さ $d$ で $1/8$ になるとき $\alpha d$ を求めよ。

:::solution
表面での損失を除けば $e^{-\alpha d}=1/8$。対数を取り $\alpha d=\log 8=3\log 2$ である。
:::
::::

::::exercise{level=basic}
$\eps_\infty=1,\gamma=0,\omega=2\omega_p$ での Drude 誘電率を求めよ。

:::solution
:eqref[eq-optics-drude] より $1-1/4=3/4$。実の正の屈折率は $\sqrt3/2$ である。
:::
::::

::::exercise{level=advanced}
Lorentz 模型の極が下半平面にあることを $\gamma<2\omega_0$ の場合に示せ。

:::solution
二次方程式 $\omega^2+i\gamma\omega-\omega_0^2=0$ の根は $\pm\sqrt{\omega_0^2-\gamma^2/4}-i\gamma/2$。$\gamma>0$ なら両根の虚部が負となる。
:::
::::

::::exercise{level=advanced}
二次分散だけのガウスパルスについてピーク強度と幅の積が保存されることを示せ。

:::solution
絶対値の二乗の前置因子は $(1+b^2)^{-1/2}$、半幅は $T_0(1+b^2)^{1/2}$ なので積は初期値と等しい。
:::
::::

::::exercise{level=advanced}
小球の半径を二倍、真空波長を二倍にしたとき無分散の散乱断面積はいくらになるか。

:::solution
両方で $k_0a$ が同じで小球近似は維持される。$a^6/\lambda_0^4$ は $2^6/2^4=4$ 倍になる。
:::
::::

::::exercise{level=advanced}
$\beta(\omega)=n_c\omega/c+B\omega^3$ とする。中心 $\omega_c$ の群遅延と群速度分散を求めよ。

:::solution
一階微分より距離 $L$ の群遅延は $L(n_c/c+3B\omega_c^2)$。二階微分より $\beta_2=6B\omega_c$ である。
:::
::::

::::exercise{level=advanced}
正の二次スペクトル位相 $B_{\rm in}=3T_0^2$ を持つガウスパルスを、$\beta_2<0$ の一様媒質に入れる。二次位相がゼロになる長さと、その位置の時間半値全幅を入力半値全幅で割った値を求めよ。

:::solution
二次位相は $B_{\rm out}=3T_0^2+\beta_2 L$ なので、消える長さは $3T_0^2/|\beta_2|$。時間半値全幅の共通因子は $2\sqrt{\ln2}\,T_0$ である。入力には幅の倍率 $\sqrt{1+3^2}=\sqrt{10}$ があるが、補償位置にはその倍率がない。従って比は $1/\sqrt{10}$ になる。
:::
::::

## 確認テスト

::::quiz
受動媒質を $e^{ik_0(n+i\kappa)z-i\omega t}$ で表す。減衰に必要な符号を選べ。

- [x] $\kappa>0$
- [ ] $\kappa<0$
- [ ] $\kappa=0$

:::explanation
空間振幅の絶対値は $e^{-k_0\kappa z}$。$k_0>0$ で $z$ とともに減るには $\kappa>0$ が必要である。
:::
::::

::::quiz{answer="5" tol="1e-6"}
$\alpha=200\,\mathrm{m^{-1}}$ の強度の $1/e$ 長を mm 単位で答えよ。

:::explanation
$\alpha z=1$ が条件で、$z=1/200\,\mathrm m=5\,\mathrm{mm}$ となる。
:::
::::

::::quiz
分散長で変換限界のガウスパルスはどうなるか。

- [x] 幅は $\sqrt2$ 倍になる。
- [ ] 幅は半分になる。
- [ ] 幅は減り、エネルギーが増える。

:::explanation
$\abs{\beta_2}z/T_0^2=1$ なので $T_0\sqrt{1+b^2}=\sqrt2 T_0$ となる。吸収を無視した式ではエネルギーは保存する。
:::
::::

::::quiz{answer="16" tol="1e-6"}
真空波長を半分にすると Rayleigh 散乱断面積は何倍か。無分散として答えよ。

:::explanation
同じ小球に対して $\sigma_{\mathrm{sca}}\propto\lambda_0^{-4}$ なので $(1/2)^{-4}=16$ 倍である。
:::
::::

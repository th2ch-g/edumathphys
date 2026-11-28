---
title: 恒星構造の基礎方程式
description: 球対称な恒星の力学・熱力学・輸送の方程式を導き、対流と輻射の役割を比較する。
notation:
  - symbol: 'm(r)'
    meaning: 半径内の質量
  - symbol: 'P'
    meaning: 全圧力
  - symbol: 'L'
    meaning: 光度
  - symbol: '\kappa_R'
    meaning: Rosseland 平均不透明度
  - symbol: '\nabla'
    meaning: 温度の圧力に対する対数勾配
---

## この章の目標

- 静水圧平衡とビリアル定理を導く。
- 物質と放射の状態方程式を使い分ける。
- 輻射拡散と対流を判定し、混合距離理論で輸送量を見積もる。

## 静水圧平衡とビリアル定理

### 球殻に働く力

SI 単位を用いる。非回転、球対称でゆっくり変化する星では、厚さ $\dd r$ の球殻の単位体積当たりの力は $-\dd P/\dd r-\rho Gm/r^2$ である。運動の時間尺度が $\sqrt{R^3/(GM)}$ より長ければ慣性項を無視できる。こうして :term[静水圧平衡]{en="hydrostatic equilibrium" yomi="せいすいあつへいこう"} の式を得る。

$$
P'=-Gm\rho/r^2,\qquad m'=4\pi r^2\rho . \label{eq-cosmo-hydrostatic}
$$

近い半径で $g=Gm/r^2$ と $\rho$ が定数なら、圧力が変化する長さは $H_P=P/(\rho g)$ である。球殻を内側から組み立てる仕事は $\Omega=-\int_0^M Gm\,\dd m/r$ である。:eqref[eq-cosmo-hydrostatic] に $4\pi r^3$ を掛けて積分し、中心で $r^3P\to0$ を使うと

$$
3\int P\,\dd V+\Omega=4\pi R^3P_s . \label{eq-cosmo-virial-star}
$$

単原子理想気体なら内部エネルギー $U=(3/2)\int P\,\dd V$ である。この積分関係を :term[ビリアル定理]{en="virial theorem" yomi="びりあるていり"} と呼ぶ。外圧 $P_s=0$ なら全エネルギーは $\Omega/2=-U$ なので、収縮に伴って光を放出しつつ温度が上がる。放射が優勢なら内部エネルギーは $3\int P_{\rm rad}\dd V$ となり、同じ係数 $2$ のビリアル関係は使えない。

:::example[一様密度球]
$m=M(r/R)^3$ を積分すると $P=P_s+3GM^2(1-r^2/R^2)/(8\pi R^4)$ となる。$r=R/2$ では過剰圧力は中心の $3/4$ である。
:::


:::example[太陽の時間尺度]
$G=6.67\times10^{-11}\,\mathrm{m^3\,kg^{-1}\,s^{-2}}$、$R=7.0\times10^8\,\mathrm m$、$M=2.0\times10^{30}\,\mathrm{kg}$ なら $\sqrt{R^3/(GM)}\simeq1.6\times10^3\,\mathrm s$ である。
:::

### 非定常なビリアル式と圧力境界

固定質量の流体球を考え、物質面での外圧を一様な $P_s$ とする。慣性モーメント $I=\int r^2\dd m$ を流体粒子に沿って二回微分すると
$$
\frac12\ddot I=2K+\int \vb r\cdot\dot{\vb v}\dd m,\qquad
K=\frac12\int |\vb v|^2\dd m .
$$
運動方程式の圧力項は $-\int\vb r\cdot\grad P\dd V$ である。$\divg(P\vb r)=3P+\vb r\cdot\grad P$ を使い、表面で $\int\vb r\cdot\dd\vb A=3V$ とすると $3\int P\dd V-3P_sV$ になる。重力項は殻に対して $-\int_0^M Gm\,\dd m/r(m)=\Omega$ である。従って
$$
\frac12\ddot I=2K+3\int P\dd V-3P_sV+\Omega .
$$
表面圧と内部圧は差で効く。定常で流れが小さく、単原子理想気体であれば $3\int P\dd V=2U$ となる。そのとき $P_s\simeq0$ なら $E=U+\Omega=\Omega/2$ である。収縮で $\Omega$ が下がると $U$ は上がり、減った重力エネルギーの半分が放射できる。この議論は一連の準静的平衡を仮定する。

:::example[密度一定の圧力支持]
$m(r)=M(r/R)^3$ を静水圧式へ代入して $P(r)-P_s=3GM^2(1-r^2/R^2)/(8\pi R^4)$ を得る。この超過圧を積分すると $3\int(P-P_s)\dd V=3GM^2/(5R)$ となる。一方、殻の結合を積分して $\Omega=-3GM^2/(5R)$ となるので、定常式は境界圧の値にかかわらず閉じている。
:::

### 動径変位の安定性を試す

全質量を保存する相似変位 $r\to\lambda r$ を考える。断熱指数 $\Gamma$ が一定の熱量方程式 $u=P/[(\Gamma-1)\rho]$ を使い、各質量殻のエントロピーを固定する。体積が $\lambda^3$ 倍になるので内部エネルギーは $U_0\lambda^{-s}$、$s=3(\Gamma-1)$ となる。表面圧を無視し、重力エネルギーを $\Omega_0/\lambda$ とすると
$$
E'(\lambda)=-sU_0\lambda^{-s-1}-\Omega_0\lambda^{-2}.
$$
平衡では $\Omega_0=-sU_0$ である。二階微分は $E''(1)=s(s-1)U_0$ となる。$\Gamma<4/3$ ではこの試験変位だけで局所極小が失われる。$\Gamma>4/3$ での正符号は相似変位に対する判定であり、一般の変位を扱うには運動方程式の線形化が必要である。

:::example[圧力の硬さの比較]
$\Gamma=5/3$ では $s=2$ だから $E''(1)=2U_0$ である。$\Gamma=6/5$ では $s=3/5$ だから $E''(1)=-6U_0/25$ である。どちらも同じ質量で平衡を構成しても、圧縮に対する回復力の符号が異なる。
:::

### 質量殻の動径固有振動

ビリアル式の相似変位は星全体を同じ比率で動かす。星の一部だけが圧縮される場合にも安定かを調べるには、動径に依存する変位を許す。球対称な非回転星を Newton 重力で扱い、熱輸送や組成変化は一周期の間に効かないとする。静止した殻の半径を $r$ とし、質量座標 $m$ を固定した変位を $\xi(r,t)=r\eta(r)e^{-i\omega t}$ と書く。振幅は複素数で扱い、物理変位は実部をとる。固定質量での変化を $\Delta$ で表す。

薄い球殻の体積は $4\pi r^2\dd r$ だから、一次の相対変化は $2\xi/r+\partial_r\xi$ である。従って $\Delta\rho/\rho=-(3\eta+r\eta')e^{-i\omega t}$ となり、局所断熱応答で $\Delta P=-\Gamma_1P(3\eta+r\eta')e^{-i\omega t}$ である。ここで $\Gamma_1=(\partial\log P/\partial\log\rho)_{s,\mathrm{composition}}$ であり、星の背景内で $r$ とともに変化してよい。

質量殻の運動方程式を $\ddot r=-4\pi r^2\partial_mP-Gm/r^2$ と書く。固定 $m$ で一次化すると
$$
\ddot\xi=-8\pi r\xi\,\partial_mP
              -4\pi r^2\partial_m\Delta P+\frac{2Gm}{r^3}\xi .
$$
背景静水圧式 $\partial_mP=-Gm/(4\pi r^4)$ と $\partial_m=(4\pi r^2\rho)^{-1}\partial_r$ を使う。圧力の面積変化と重力の半径変化の項は足して $4Gm\xi/r^3$ である。時間因子を除いて $\Delta P$ を代入し、$r^3$ を掛けると
$$
\dv{}{r}\!\left(\Gamma_1P r^4\dv{\eta}{r}\right)
+r^3\dv{[(3\Gamma_1-4)P]}{r}\eta+\omega^2\rho r^4\eta=0 .
\label{eq-cosmo-radial-pulsation}
$$
中心の正則な解では $\eta$ は有限であり、滑らかな球対称変位なら $\eta'(0)=0$ である。真空を外部とする自由表面では、変位した表面の圧力が外圧と一致するので $\Delta P(R)=0$ である。$\eta$ が有界で有限なエネルギーを持ち、$P=0$ の正則な自由表面を考えると、部分積分の境界項 $\Gamma_1P r^4\eta^*\eta'$ は零である。$P_s>0$ の場合は $\Delta P(R)=0$ を用いて別の表面項が残る。

:::proposition[断熱動径モードのエネルギー積分]
正則な自由表面を持ち、$\Gamma_1P>0$、$\rho>0$ が内部で成り立つとする。上の境界条件を満たす非零の固有モードについて
$$
\omega^2=
\frac{\int_0^R \Gamma_1P r^4|\eta'|^2\dd r
 -\int_0^R r^3[(3\Gamma_1-4)P]'|\eta|^2\dd r}
{\int_0^R\rho r^4|\eta|^2\dd r}.
$$
:::

:::proof
:eqref[eq-cosmo-radial-pulsation] に $\eta^*$ を掛けて積分する。第一項を部分積分し、両端で零となる項を落とせば分子と分母が出る。分母は正で分子は実数だから $\omega^2$ は実である。従って負の固有値なら $\omega$ は純虚数で、時間的に成長する解を持つ。
:::

$\Gamma_1$ 一定かつ $\Gamma_1>4/3$ なら、内部で $P'<0$ なので分子の第二項も正である。この場合、全動径固有モードの $\omega^2$ は正になる。$\Gamma_1$ が場所で変わる場合、単にその中心値を $4/3$ と比べるだけではこの判定は使えない。

:::example[一様密度の厳密な相似振動]
一様密度星で $\Gamma_1$ を一定とし、$\eta(r)$ を非零の定数にする。$P'=-\rho GM r/R^3$ なので :eqref[eq-cosmo-radial-pulsation] は全半径で $\omega^2=(3\Gamma_1-4)GM/R^3$ に簡約する。$\Gamma_1=5/3$ なら $\omega=\sqrt{GM/R^3}$ であり、$\Gamma_1=4/3$ なら零になる。一様な平衡密度を仮定していても、摂動の圧縮応答は有限の $\Gamma_1$ によって定めている。
:::

:::example[相似モードの臨界指数近く]
同じ星で $\Gamma_1=4/3+1/30$ なら $\omega^2=GM/(10R^3)$ で、周期は $2\pi\sqrt{10R^3/(GM)}$ となる。$\Gamma_1=4/3-1/30$ では同じ式の値が負になり、成長時間は $\sqrt{10R^3/(GM)}$ である。安定側の周期と不安定側の指数時間は、同じ固有値の絶対値から得るが $2\pi$ の差がある。
:::

### 動径モードの直交性と単位

異なる固有値 $\omega_1^2\ne\omega_2^2$ の二つの正則解を $\eta_1,\eta_2$ とする。:eqref[eq-cosmo-radial-pulsation] の第一の式に $\eta_2^*$、第二の共役式に $\eta_1$ を掛けて引く。圧力勾配の項は同じなので消え、微分項の差は
$$
\dv{}{r}[\Gamma_1 P r^4(\eta_2^*\eta_1'-\eta_1\eta_2^{*\prime})]
$$
となる。中心と自由表面で境界値が零なら、全積分は $(\omega_1^2-\omega_2^2)\int_0^R\rho r^4\eta_2^*\eta_1\dd r=0$ を与える。従って重み $\rho r^4$ の積分で直交する。この重みは、物理変位が $r\eta$ であり、運動エネルギーの体積積分にさらに $4\pi r^2$ が掛かることと一致する。$\omega^2$ を与えるエネルギー積分で、分子と分母の次元比は $\mathrm{s^{-2}}$ である。

:::example[表面変位とモードの運動エネルギー]
一様密度球の相似振動で $\eta$ は定数とする。速度振幅は $\omega r|\eta|$ なので、速度が最大の位相での運動エネルギーは $\frac12\omega^2|\eta|^2\int r^2\dd m=3M\omega^2R^2|\eta|^2/10$ である。表面の速度振幅だけを全質量に割り当てる評価の $3/5$ になる。星の内部では表面より速度が小さいためである。
:::

## 質量保存とエネルギー保存

### 質量座標での熱収支

質量が保存される各球殻を内包質量の質量座標 $m$ で識別すると $\pdv{r}{m}=(4\pi r^2\rho)^{-1}$、$\pdv{P}{m}=-Gm/(4\pi r^4)$ となる。単位質量当たりの内部エネルギーを $u$、核反応による静止質量の減少率をエネルギーに換算して $\epsilon_{\rm nuc}$、ニュートリノ損失を $\epsilon_\nu$ と置く。球殻へ入る光度と出る光度の差を熱力学第一法則に代入すると

$$
 \pdv{L}{m}=\epsilon_{\rm nuc}-\epsilon_\nu-\left(\pdv{u}{t}\right)_m-P\left(\pdv{(1/\rho)}{t}\right)_m . \label{eq-cosmo-stellar-energy}
$$

固定組成では最後の二項は $-T(\partial s/\partial t)_m$ に等しい。組成変化を内部エネルギーと核発熱の双方に重複して計上しない。ここでの $L$ は、$\epsilon_\nu$ に数えたニュートリノを除く輸送光度である。中心で $r=L=0$、光球で圧力と温度を大気解に接続する。

:::example[一様な加熱]
熱平衡で $\epsilon_{\rm nuc}=2\times10^{-4}\,\mathrm{W\,kg^{-1}}$ と一定なら $m=10^{30}\,\mathrm{kg}$ で $L=2\times10^{26}\,\mathrm W$ である。
:::


:::example[収縮光度]
外圧零の一様単原子気体球で $\Omega=-3GM^2/(5R)$ なら、核加熱なしで $L=-(3GM^2/10R^2)\dot R$ である。収縮率の絶対値を二倍にすれば同じ半径での光度は二倍になる。
:::

### 光度の中心条件と発熱殻

球対称の質量座標で、中心からの光度を $L(m)$ と定める。有効発熱率を $\epsilon_{\rm eff}$ と書けば、中心で有限な $\epsilon_{\rm eff}$ のもとで
$$
L(m)=\int_0^m\epsilon_{\rm eff}(m')\dd m'.
$$
$\epsilon_{\rm eff}=\epsilon_c+\order{r^2}$ なら $L=\epsilon_c m+\order{m^{5/3}}$ である。滑らかな球対称場は中心で勾配が零であり、$m\propto r^3$ だからである。従って $F=L/(4\pi r^2)\propto r$ であり、中心で光束は有限で零となる。数値積分では $r=0$ の除算を避けてこの級数を境界値に使う。

:::example[厚さのある燃焼殻]
$m_1<M$、$m_2>m_1$ として、$\epsilon_{\rm eff}=\epsilon_0$ が $m_1<m<m_2$ のみにあるとする。$m<m_1$ で $L=0$、殻内で $L=\epsilon_0(m-m_1)$、外側で $L=\epsilon_0(m_2-m_1)$ となる。$m_2-m_1=0.02M$ なら表面光度は $0.02M\epsilon_0$ である。外層の持つ質量が大きくても、熱を蓄えず発熱しなければ光度勾配は零である。
:::

### 局所時間尺度の区別

音速 $c_s$ が各層での圧力変化を伝え、熱輸送は別の時間でエントロピーを変える。星全体の基準値は $t_{\rm dyn}\sim\sqrt{R^3/(GM)}$、$t_{\rm KH}\sim GM^2/(RL)$、$t_{\rm nuc}\sim qfMc^2/L$ である。最後の $f$ は燃料として使える質量分率、$q$ はその質量をエネルギーに換算する効率である。局所には、質量区間 $\Delta m$ の熱容量を使った $t_{\rm th}\sim c_PT\Delta m/L$ が使える。光度の変化が短時間なら、その外側の熱容量の小さな層だけが熱的に応答できる。

:::example[小さな外層の熱時間]
外層の温度が $10^5\,\mathrm K$、定圧比熱が $3\times10^4\,\mathrm{J\,kg^{-1}K^{-1}}$、質量が $10^{21}\,\mathrm{kg}$ とする。$L=3\times10^{26}\,\mathrm W$ なら熱時間は $10^4\,\mathrm s$ である。外層を通る光度の変動と星全体の進化時間を比較するときは、この局所熱容量を用いる。
:::

### 収縮と組成変化を含む熱力学の帳尻

質量殻の一質量当たりの内部エネルギーを $u(s,v,Y_1,\ldots,Y_N)$ とし、$v=1/\rho$、$Y_i$ は一質量当たりの粒子数とする。独立な組成座標に関して $\mu_i$ を定めれば第一法則は
$$
\dd u=T\dd s-P\dd v+\sum_i\mu_i\dd Y_i
$$
である。粒子数に制約があれば $\sum_i\mu_i\dd Y_i$ は許された組成変化に沿った微分として扱う。$\mu_i$ が粒子一個当たりのエネルギーなら、和は $\mathrm{J\,kg^{-1}}$ の次元を持つ。内部エネルギーにどの静止質量を含めたかで、反応による見かけの発熱項が変わる。

静止質量を除いた $u$ で熱を追い、核反応で減った静止質量から持ち込む率を $\epsilon_{\rm nuc}$ と定める。圧縮仕事は一質量当たり $-P\,\partial_t v$ である。熱と放射の流出 $-\partial_mL$ とニュートリノ損失 $-\epsilon_\nu$ を加えると
$$
\partial_m L=\epsilon_{\rm nuc}-\epsilon_\nu
       -(\partial_t u)_m-P(\partial_t v)_m .
$$
シェルを囲う二つの物質面の間で境界が動くが、質量座標を固定しているため、熱輸送を $L$ で数えれば別の移流項を加えない。組成が固定なら最後の二項は $-T\partial_t s$ である。組成が変わる場合は、採用した $u$ に対応する $\sum_i\mu_i\partial_t Y_i$ を含めて変換する。一方、静止質量を含む $u$ を使うと、反応の静止質量減少はすでに $(\partial_t u)_m$ に入るので同じ $\epsilon_{\rm nuc}$ をさらに足すと二重計上になる。

:::example[断熱圧縮と発光の区別]
発熱も放射もない理想気体の殻を断熱的に圧縮する。単原子気体なら $u\propto\rho^{2/3}$ で、密度を八倍にすると比内部エネルギーは四倍である。一方、$\partial_t u=-P\partial_t v$ なので、エネルギー式で $-\partial_t u-P\partial_t v=0$ となる。圧縮で温度が上がることだけから光度が増えるとはいえない。
:::

:::example[収縮で表面へ出せるエネルギー]
単原子理想気体が支える一様密度球を、準静的な相似模型として $R_0$ から $R_0/2$ へ縮める。外圧と発熱は無視する。ビリアル式より $E=-3GM^2/(10R)$ なので、放出エネルギーは $E(R_0)-E(R_0/2)=3GM^2/(10R_0)$ である。重力エネルギーだけの減少量はその二倍で、差は内部エネルギーの増加になる。
:::

### 中心境界の級数を計算する

中心付近で $\rho=\rho_c+\rho_2r^2+\order{r^4}$、$\epsilon_{\rm eff}=\epsilon_c+\epsilon_2r^2+\order{r^4}$ とする。質量と光度の式を積分すると
$$
\begin{aligned}
m(r)&=4\pi\left(\frac{\rho_c r^3}{3}+\frac{\rho_2r^5}{5}\right)+\order{r^7},\\
L(r)&=4\pi\left(\frac{\rho_c\epsilon_c r^3}{3}
 +\frac{(\rho_c\epsilon_2+\rho_2\epsilon_c)r^5}{5}\right)+\order{r^7}.
\end{aligned}
$$
静水圧式の最初の非零次は $P'=-4\pi G\rho_c^2 r/3$ なので $P(r)=P_c-2\pi G\rho_c^2r^2/3+\order{r^4}$ である。輻射拡散を使える中心で $F=L/(4\pi r^2)\simeq\rho_c\epsilon_c r/3$ を $F=-16\sigma T^3 T'/(3\kappa\rho)$ と比べると $T(r)=T_c-\kappa_c\rho_c^2\epsilon_c r^2/(32\sigma T_c^3)+\order{r^4}$ である。$\sigma$ は Stefan–Boltzmann 定数とする。温度勾配の係数の比をとれば、中心で $m$ と $L$ がともに零でも $L/m\to\epsilon_c$ と分かる。

:::example[数値の開始点と倍距離]
十分中心に近い関数を半径 $r_0$ と $2r_0$ で比べる。$L\propto r^3$ なので比は先頭で八倍であるが、フラックス $L/(4\pi r^2)$ の比は二倍、中心からの温度差の比は四倍である。この三つの比は中心境界を実装したとき、数値の桁を独立に調べるために使える。
:::

### 時間発展で静水圧近似を使う区間

静水圧式では圧力と重力を釣り合わせ、運動の慣性項を落とす。半径が時間 $t_{\rm ev}$ で自身と同程度だけ変わる滑らかな変化なら、慣性加速度の尺度は $R/t_{\rm ev}^2$、重力加速度は $GM/R^2$ である。その比は $(t_{\rm dyn}/t_{\rm ev})^2$、$t_{\rm dyn}=\sqrt{R^3/(GM)}$ である。この評価には、急な振動や衝撃が別に重なっていないことが必要である。圧力の応答が局所でも速いかを確認するため、層の厚さに対する音の横断時間も比べる。

強い熱源がある場合、熱的な定常性も独立に点検する。星の表面光度を $L_{\rm surf}$、核燃焼とニュートリノ損失の積分を $L_{\rm nuc},L_\nu$ とすると、組成を固定した近似で
$$
L_{\rm surf}=L_{\rm nuc}-L_\nu-\int_0^M T(\partial_t s)_m\dd m .
$$
ここで $\partial_t s<0$ の層は正の光度寄与を持つ。中心が静水圧に近くてもエントロピーが変わるなら、外面の光度と核発熱の差は残る。温度上昇とエントロピー上昇を混同しない。圧縮されながら熱を失う場合は、温度が上がってもエントロピーが下がることがある。

:::example[準静的収縮で発熱源を見積もる]
ある時刻に $L_{\rm nuc}=0.3L_{\rm surf}$、$L_\nu=0.1L_{\rm surf}$ とする。上の式では $-\int T\partial_t s\dd m=0.8L_{\rm surf}$ である。仮に $t_{\rm ev}/t_{\rm dyn}=10^3$ なら慣性の尺度は重力の $10^{-6}$ にすぎない。したがって力学は静水圧で扱えても、熱収支でエントロピー項を落とすと光度の大半を取り逃す。
:::

## 恒星内部の状態方程式（理想気体・輻射・縮退）

### 粒子数と圧力の分担

完全電離した理想気体で核種 $i$ の質量分率を $X_i$、質量数を $A_i$、電荷数を $Z_i$ とする。核数密度が $\rho X_i/(A_i m_u)$ であり、一核当たりに自由電子が $Z_i$ 個あるので、$P_g=\rho k_BT/(\mu m_u)$、$\mu^{-1}=\sum_i(1+Z_i)X_i/A_i$ となる。電離が不完全なら各イオン段階の占有率が必要である。

:ref[def-cosmo-planck-function] の等方放射では $P_{\rm rad}=a_RT^4/3$、$a_R=4\sigma/c$ である。非相対論的な電子が強く縮退した場合、状態数から $n_e=p_F^3/(3\pi^2\hbar^3)$ であり、運動量フラックス $P=2\int \dd^3p\,pv/(3(2\pi\hbar)^3)$ に $v=p/m_e$ を入れると $P=n_ep_F^2/(5m_e)$ となる。$v\simeq c$ では同じ積分から $P=n_ep_Fc/4$ を得る。条件は $k_BT\ll E_F$ である。有限温度での補正は :eqref[eq-stm-fermi-gas] から求められる。ここで $E_F$ は静止質量を除く Fermi エネルギーである。

:::example[電離した H–He 混合]
$X_{\rm H}=3/4$、$X_{\rm He}=1/4$ なら $\mu^{-1}=2(3/4)+(3/4)(1/4)=27/16$、$\mu=16/27$ となる。
:::


:::example[圧縮の比較]
密度を八倍にすると非相対論的縮退圧は $8^{5/3}=32$ 倍、超相対論的縮退圧は $8^{4/3}=16$ 倍になる。
:::

### 気体と放射を一緒に圧縮する

完全電離、一定組成、非縮退の単原子気体と平衡放射を混ぜる。$\mathcal R_g=k_B/(\mu m_u)$、$\beta=P_g/P$ とすると、$P_g=\mathcal R_g\rho T$、$P_r=a_RT^4/3$ より $\chi_\rho=(\partial\log P/\partial\log\rho)_T=\beta$、$\chi_T=(\partial\log P/\partial\log T)_\rho=4-3\beta$ となる。単位質量内エネルギーは $u=(3/2)\mathcal R_gT+a_RT^4/\rho$ である。$c_V=(\partial u/\partial T)_\rho$ を評価すると $\rho Tc_V/P=(3/2)\beta+12(1-\beta)$ となる。

第一法則 $T\dd s=\dd u+P\dd(1/\rho)$ を $\dd s=0$ で用いる。密度に対する $u$ の微分と仕事項から、$\rho T c_V\dd\log T=(P_g+4P_r)\dd\log\rho$ となる。従って
$$
\Gamma_3-1=\left(\pdv{\log T}{\log\rho}\right)_s
=\frac{4-3\beta}{12-21\beta/2},\qquad
\Gamma_1=\left(\pdv{\log P}{\log\rho}\right)_s
=\beta+(4-3\beta)(\Gamma_3-1).
$$
:term[断熱温度勾配]{en="adiabatic temperature gradient" yomi="だんねつおんどこうばい"} は $\nabla_{\rm ad}=(\Gamma_3-1)/\Gamma_1$ である。$\beta=1$ からは $\Gamma_1=5/3,\nabla_{\rm ad}=2/5$、$\beta=0$ からは $\Gamma_1=4/3,\nabla_{\rm ad}=1/4$ を再現する。放射圧の寄与がある大質量星では $2/5$ を固定して対流を判断すると、半径ごとの安定性を誤ることがある。部分電離や縮退がある域では $u$ と $P$ の両方をその状態方程式で微分する。

:::example[気体と放射の圧力が等しい層]
$\beta=1/2$ なら $\Gamma_3-1=(5/2)/(27/4)=10/27$、$\Gamma_1=77/54$、従って $\nabla_{\rm ad}=20/77\simeq0.260$ になる。
:::

### 熱力学微分と一般の断熱勾配

質量組成を固定し、比体積 $v=1/\rho$ をとる。Gibbs の関係 $\dd h=T\dd s+v\dd P$ で $h$ は単位質量当たりのエンタルピーである。Gibbs 自由エネルギー $g=h-Ts$ の混合偏微分を交換すると $(\partial s/\partial P)_T=-(\partial v/\partial T)_P$ を得る。$c_P=T(\partial s/\partial T)_P$ および $\delta=-(\partial\log\rho/\partial\log T)_P$ を用いると
$$
\dd s=c_P\,\dd\log T-\frac{P\delta}{\rho T}\,\dd\log P,
\qquad \nabla_{\rm ad}=\left(\pdv{\log T}{\log P}\right)_s=\frac{P\delta}{\rho Tc_P}.
$$
温度が上がっても電離にエネルギーが使われる領域では $c_P$ が増える。このとき同じ圧力変化に対して断熱的な温度変化は小さくなり、対流の判定へ影響する。

:::example[理想単原子気体との比較]
$P/\rho=R_sT$、$\delta=1$、$c_P=5R_s/2$ から $\nabla_{\rm ad}=2/5$ となる。模型として $P/\rho$ と $\delta$ を固定したまま $c_P$ が二倍になれば $\nabla_{\rm ad}=1/5$ になる。実際の部分電離層では $\delta$ も変化するので、両方の微分を同じ状態方程式から評価する。
:::

### 電離を含む平均分子量

水素だけからなる非縮退気体で、核の数密度を $n_H$、電離度を $x=n_p/n_H=n_e/n_H$ とする。重い核の質量を $m_H$ とすれば、粒子数密度は $(1+x)n_H$ なので $\mu=1/(1+x)$ である。電離エネルギーを $\chi$ とし、励起を無視すると比内部エネルギーは
$$
u=\frac32(1+x)\frac{k_BT}{m_H}+\frac{x\chi}{m_H}.
$$
密度を固定した微分に $(\partial x/\partial T)_\rho$ が入る。例えば Saha 形 $x^2/(1-x)=A(T)/n_H$、$A\propto T^{3/2}e^{-\chi/(k_BT)}$ の両辺を対数微分すると
$$
\left(\pdv{x}{\log T}\right)_\rho=\frac{x(1-x)}{2-x}\left(\frac32+\frac{\chi}{k_BT}\right).
$$
$x=0$ と $x=1$ の極限で電離度変化の寄与が消え、中間で熱容量が大きくなる。

:::example[半電離での比熱]
$x=1/2$、$\chi/(k_BT)=10$ では $\partial x/\partial\log T=23/12$ である。電離エネルギー項だけでも $c_V/(k_B/m_H)$ に $10\times23/12\simeq19.2$ を加える。一方、電離度を固定した並進項は $9/4$ である。
:::

### 三種類の断熱応答を状態方程式から得る

組成固定の $P(\rho,T)$ で
$$
\chi_\rho=\left(\pdv{\log P}{\log\rho}\right)_T,\qquad
\chi_T=\left(\pdv{\log P}{\log T}\right)_\rho
$$
とする。Helmholtz 自由エネルギーの微分 $\dd f=-s\dd T-P\dd v$ から $(\partial s/\partial v)_T=(\partial P/\partial T)_v$ である。$c_V=T(\partial s/\partial T)_v$ と $\dd v=-\dd\log\rho/\rho$ を使うと
$$
\dd s=c_V\dd\log T-\frac{P\chi_T}{\rho T}\dd\log\rho.
$$
断熱線では $\Gamma_3-1=(\partial\log T/\partial\log\rho)_s=P\chi_T/(\rho T c_V)$ である。$P$ の全微分から $\Gamma_1=\chi_\rho+\chi_T(\Gamma_3-1)$ を得る。$\nabla_{\rm ad}=(\Gamma_3-1)/\Gamma_1$ であり、別の記法 $\Gamma_2$ では $(\Gamma_2-1)/\Gamma_2=\nabla_{\rm ad}$ と定める。断熱音速は $c_s^2=(\partial P/\partial\rho)_s=\Gamma_1P/\rho$ である。従って $\Gamma_1$ は力学、$\nabla_{\rm ad}$ は温度勾配の判定に直接現れる。

非縮退単原子理想気体と放射が熱平衡にあり $P=P_g+P_r$、$\beta=P_g/P$ とする。$P_g=R_s\rho T$ と $P_r=a_{\rm R}T^4/3$ から $\chi_\rho=\beta$、$\chi_T=4-3\beta$ である。気体の比内部エネルギーは $(3/2)R_sT$、放射の比内部エネルギーは $a_{\rm R}T^4/\rho$ だから
$$
c_V=\frac{P}{\rho T}\left[\frac32\beta+12(1-\beta)\right].
$$
代入すると
$$
\Gamma_3-1=\frac{4-3\beta}{12-(21/2)\beta},\qquad
\Gamma_1=\beta+\frac{(4-3\beta)^2}{12-(21/2)\beta}.
$$
純気体極限では $\Gamma_3-1=2/3,\Gamma_1=5/3$、放射優勢極限では $\Gamma_3-1=1/3,\Gamma_1=4/3$ である。これらは $\Gamma_1$ を圧力で単純に重み付けして平均すると一般には得られない。二成分の温度がともに変わり、両方の熱容量が断熱条件に入るためである。

:::example[圧力を等しく分担する気体と放射]
$\beta=1/2$ なら $\Gamma_3-1=10/27$、$\Gamma_1=77/54$、$\nabla_{\rm ad}=20/77$ である。同じ $P/\rho$ で単原子理想気体と比べると、音速二乗の比は $(77/54)/(5/3)=77/90$ である。光速で動く成分が圧力を担っても、Newton 的なモデルで混合物の断熱応答はこの比だけ変わる。
:::

### 非理想性が小さいという判定

一種類のイオンの数密度を $n_i$、平均間隔の尺度を $a_i=(3/(4\pi n_i))^{1/3}$ とする。電荷 $Ze$ のイオンについて $a_i$ における Coulomb エネルギーと熱エネルギーの比
$$
\Gamma_{\rm C}=\frac{Z^2 e^2}{4\pi\epsilon_0a_i k_BT}
$$
を作る。弱結合 $\Gamma_{\rm C}\ll1$ なら典型的な相互作用エネルギーは熱より小さい。電子 Fermi 運動量 $p_F=\hbar(3\pi^2n_e)^{1/3}$ から電子の Fermi 運動エネルギーを $E_F=\sqrt{p_F^2c^2+m_e^2c^4}-m_ec^2$ と作ると、$k_BT/E_F\ll1$ は電子の熱励起が Fermi 面近傍に限られる判定である。イオンの強結合と電子の縮退は、別々の比で調べる。強結合の定量的な自由エネルギーや混合物の結晶化は多体相関を必要とするため、理想気体の式の補正だけからは決定できない。

:::example[密度を八倍にして圧縮する]
温度、電離、組成を固定する。$a_i$ は半分になり、$\Gamma_{\rm C}$ は二倍になる。非相対論的な電子の $E_F\propto n_e^{2/3}$ は四倍なので、$k_BT/E_F$ は四分の一になる。圧縮するとイオンの理想気体性が弱まる一方、電子の縮退近似はよくなることがある。
:::

### 化学組成と状態変数の整合性

質量分率 $X_i$、質量数 $A_i$、イオンの平均電荷 $Z_i^*$ を用いる。核の総質量を $A_i m_u$ と近似すると、原子核の数密度は $n_i=\rho X_i/(A_i m_u)$ である。電荷中性により $n_e=\sum_i Z_i^*n_i$ なので $1/\mu_e=\sum_i X_i Z_i^*/A_i$ である。非縮退でイオン・電子を理想気体として足すなら $1/\mu=\sum_i X_i(1+Z_i^*)/A_i$ である。部分電離の $Z_i^*$ は $P,T$ に依存し、この依存を固定した微分と平衡に沿う微分は一般に異なる。

分子ができる領域では、たとえば水素核二個が一分子として一つの並進粒子を作るため、核の数を並進粒子数と同一視できない。各分子を一種類の粒子として力学的圧力に数え、核の数は反応の保存則に用いる。気体が十分希薄なら分子の内部自由度は並進圧力を直接増やさないが、熱容量を大きくして断熱指数を変える。励起や解離の時間が振動より長ければ、固有振動の $\Gamma_1$ に平衡の化学変化を含めない。

:::example[核の組成一定で電離だけが変わる層]
質量分率が水素 $3/4$、He $1/4$ で両方が完全電離しているとき、$1/\mu=2(3/4)+(3/4)(1/4)=27/16$、$1/\mu_e=3/4+(1/2)(1/4)=7/8$ である。双方が中性原子なら $1/\mu=3/4+(1/4)(1/4)=13/16$ となる。同じ密度と温度で完全電離時の理想気体圧は中性時の $27/13$ 倍である。この比は放射圧や分子化を含まない。
:::

## 輻射によるエネルギー輸送と不透明度

### 光の拡散と平均

平均自由行程が温度変化の長さより短い局所熱平衡を仮定する。:eqref[eq-cosmo-transfer] を一反復すると $I_\nu=B_\nu-(\rho\kappa_\nu)^{-1}\mu\,\partial_rB_\nu$ となる。角度平均の $\langle\mu^2\rangle=1/3$ から $F_\nu=-4\pi\partial_rB_\nu/(3\rho\kappa_\nu)$ を得る。ゆえに、:term[Rosseland 平均]{en="Rosseland mean opacity" yomi="ろっせらんどへいきん"} は $\partial_TB_\nu$ で重み付けした調和平均
$$
\kappa_R^{-1}=\frac{\int\kappa_\nu^{-1}\partial_TB_\nu\,\dd\nu}{\int\partial_TB_\nu\,\dd\nu}
$$
である。散乱があるときの拡散には輸送不透明度を用いる。$L=4\pi r^2F$ と $\int B_\nu\dd\nu=\sigma T^4/\pi$ により

$$
\dv{T}{r}=-\frac{3\kappa_R\rho L}{16\pi a_R c r^2T^3},\qquad \nabla_{\rm rad}\defeq\left(\dv{\log T}{\log P}\right)_{\rm rad}=\frac{3\kappa_R LP}{16\pi a_R c Gm T^4}. \label{eq-cosmo-radiative-gradient}
$$

完全電離では電子散乱の $\kappa_{\rm es}=n_e\sigma_T/\rho$ が密度に依存しない。自由–自由吸収の近似 $\kappa_\nu\propto\rho T^{-1/2}\nu^{-3}(1-e^{-h\nu/k_BT})$ に $x=h\nu/k_BT$ を適用すると $\kappa_R\propto\rho T^{-7/2}$ となる。束縛–束縛や束縛–自由遷移が強ければ単一のべき則は崩れる。精密な不透明度表の作成には多数の原子準位を含む量子計算が必要であり :cite[kippenhahn-stellar] を参照する。

:::example[透明な帯域]
等重みの二帯域で $\kappa_\nu=1,9\,\mathrm{m^2\,kg^{-1}}$ なら $\kappa_R=1.8\,\mathrm{m^2\,kg^{-1}}$ である。小さい不透明度の帯域が輸送を担う。
:::


:::example[勾配の変化]
$m,P,T$ を固定して $\kappa_R$ と $L$ をそれぞれ二倍にすると、必要な $\nabla_{\rm rad}$ は四倍になる。
:::

### 拡散係数と境界で破れる条件

比強度を $I_\nu=B_\nu+\delta I_\nu$ とし、平均自由行程 $\ell_\nu=(\rho\kappa_\nu)^{-1}$ が温度の変化尺度より小さいとする。静止媒質の吸収・放射による輸送式は $\hat{\vb n}\cdot\grad I_\nu=(B_\nu-I_\nu)/\ell_\nu$ である。一次では $\delta I_\nu=-\ell_\nu\hat{\vb n}\cdot\grad B_\nu$ となる。角度の積分 $\int n_i n_j\dd\Omega=4\pi\delta_{ij}/3$ から、拡散フラックス
$$
\vb F_\nu=-\frac{4\pi\ell_\nu}{3}\pdv{B_\nu}{T}\grad T
$$
を得る。$u_\nu=4\pi B_\nu/c$ を用いれば、エネルギー密度に対する拡散係数は $D_\nu=c\ell_\nu/3$ である。散乱主体の場合も熱化が成り立つ深さで用い、強い前方散乱では角再分布を反映した輸送不透明度を使う。

:::example[二つの周波数窓]
Rosseland 平均で低不透明度帯の重みを $1/4$、高不透明度帯を $3/4$ に変える。$\kappa_1=1$、$\kappa_2=9$ を同じ単位で与えると、$\kappa_R=[1/4+(3/4)/9]^{-1}=3$ である。低不透明度帯は重みが四分の一でも、全フラックスの四分の三を運ぶ。
:::

### 不透明度が変わるときの勾配

質量殻における輻射勾配は $\nabla_{\rm rad}\propto \kappa LP/(mT^4)$ である。状態を固定した比較では
$$
\Delta\log\nabla_{\rm rad}=\Delta\log\kappa+\Delta\log L+\Delta\log P-\Delta\log m-4\Delta\log T
$$
が厳密に成り立つ。構造が緩和した後には $L,P,T$ も変わるため、ここでの比較は局所応答を表す。

:::example[不透明度の山と対流]
ある層の $\nabla_{\rm rad}=0.30$、$\nabla_{\rm ad}=0.40$ とする。組成により $\kappa$ だけが $1.5$ 倍に変われば $\nabla_{\rm rad}=0.45$ となる。これで輻射単独の温度勾配は不安定になるので、新しい構造を求める際には対流も含める。
:::

### 伝導と輻射を同じ温度差に結び付ける

散乱または吸収を何度も受ける粒子の平均自由行程 $\ell$ が勾配の長さより短ければ、等方な定常拡散で熱流を $\vb F=-K\grad T$ と書く。粒子速度 $v$ と体積熱容量 $C$ を用いた尺度は $K\simeq C v\ell/3$ である。これは温度分布 $T(\vb r-\ell\hat{\vb n})$ を一次まで展開し、$\int n_i n_j\dd\Omega/(4\pi)=\delta_{ij}/3$ で平均すると得られる。平均自由行程が粒子エネルギーに依存すれば、$C$ に対応する微分重みで $v\ell$ を積分する。縮退電子では Fermi 面近くの励起が熱を運ぶ。圧力を担う電子がすべて同じ割合で熱伝導へ寄与するとは限らない。

放射拡散を $K_{\rm rad}=16\sigma T^3/(3\rho\kappa_R)$ と書く。同じ $\grad T$ によって流れる伝導束を $K_{\rm cond}$ とすると、和の係数は $K_{\rm eff}=K_{\rm rad}+K_{\rm cond}$ である。$K_{\rm cond}$ を $16\sigma T^3/(3\rho\kappa_{\rm cond})$ と形式的に表せば、$\kappa_{\rm eff}^{-1}=\kappa_R^{-1}+\kappa_{\rm cond}^{-1}$ となる。伝導の $\kappa_{\rm cond}$ は光子の吸収係数を意味しない。並列の輸送を共通の恒星構造方程式に入れるための換算である。

:::example[必要な勾配の変化]
ある層で $K_{\rm cond}=4K_{\rm rad}$ とする。全フラックスを固定すると必要な温度勾配は輻射だけの値の $1/5$ になる。伝導は全フラックスの $4/5$ を担う。固定温度勾配で比較した場合は、逆に同じ層が運べるフラックスが五倍になる。
:::

### 灰色の境界を適用するための条件

平面近似で高さを $z$、外向きの重力加速度の大きさを $g$、内向き光学深さを $\tau(z)=\int_z^\infty\kappa\rho\dd z'$ とする。静水圧式との比で $\dd P/\dd\tau=g/\kappa$ である。$g,\kappa$ がほぼ一定なら $P(\tau)-P(0)=g\tau/\kappa$ となる。灰色の温度関係 :eqref[eq-cosmo-grey-atmosphere] と合わせて、$\tau=2/3$ で $P,T$ の境界値を与えられる。一方、$\kappa$ が圧力に依存するなら $\tau$ から単純比例で $P$ を得られない。

例えば $g$ 一定、固定した温度近似で $\kappa=\kappa_0 P^\alpha$ とし、$\alpha>-1$ とする。$P(0)=0$ から積分すると
$$
P(\tau)^{1+\alpha}=\frac{(1+\alpha)g}{\kappa_0}\tau.
$$
$\alpha>0$ なら表面では不透明度が零に近づくため、$\dd P/\dd\tau$ は大きくなるが $P(\tau)$ 自体は有限である。実際に温度も変わる場合には、$\kappa(P,T)$ と温度境界を一緒に積分する。

:::example[光学深さと圧力の非比例]
$\alpha=1$ の境界模型では、$\tau$ が $2/3$ から $8/3$ に増えると $P$ は二倍である。$\kappa$ を一定にした境界では四倍である。恒星内部と大気の接続深さを変えるとき、二つの境界式は同じ圧力を与えない。
:::

### Planck 平均と拡散の平均を使い分ける

薄く、熱源を黒体で近似できる物質が真の吸収をもつとする。角度積分した放射率は $4\pi\rho\int\kappa_\nu B_\nu(T)\dd\nu$ であり、これを $4\pi\rho\kappa_P\int B_\nu\dd\nu$ と等置すると、
$$
\kappa_P=\frac{\int_0^\infty\kappa_\nu B_\nu\dd\nu}{\int_0^\infty B_\nu\dd\nu}
$$
となる。これは :term[Planck 平均]{en="Planck mean opacity" yomi="ぷらんくへいきん"} である。放射を作るときは大きい吸収係数を持つ帯域ほど熱放射する。一方、厚い媒質の Rosseland 平均では自由行程の長い帯域の寄与が大きい。薄い媒質が別の天体によって加熱されるときには、吸収率を $\int\kappa_\nu J_\nu\dd\nu$ と評価する必要がある。吸収側の $J_\nu$ が物質の $B_\nu(T)$ でなければ、同じ Planck 平均を加熱と冷却に使えない。

:::example[光が出やすい帯域と逃げやすい帯域]
ある二帯域で $\kappa_1=2$、$\kappa_2=18$ とし、$B_\nu$ の帯域重みも $\partial_T B_\nu$ の帯域重みも便宜上それぞれ等しいとする。$\kappa_P=10$、$\kappa_R=2/(1/2+1/18)=3.6$ である。薄い状態から厚い状態に変わる際に、同じ平均を使うとこの二帯域でフラックスを大きく見誤る。
:::

## 対流と Schwarzschild の判定条件

### 浮力による局所安定性

温度勾配を $\nabla=\dd\log T/\dd\log P$ とする。流体塊を $\xi>0$ だけ持ち上げ、周囲と圧力だけを平衡させる。熱交換しない塊と周囲の温度差は $\Delta\log T=-(\nabla_{\rm ad}-\nabla)\xi/H_P$ である。固定組成かつ定圧で $\delta=-(\partial\log\rho/\partial\log T)_{P,\mu}>0$ なら密度差は $\Delta\log\rho=-\delta\Delta\log T$ となる。浮力の加速度は $-g\Delta\log\rho$ なので $\ddot\xi=-N^2\xi$、$N^2=(g\delta/H_P)(\nabla_{\rm ad}-\nabla)$ である。

伝導を無視し、全光度を輻射で運ぶ状態の判定に $\nabla=\nabla_{\rm rad}$ を用いる。$\nabla_{\rm rad}>\nabla_{\rm ad}$ が :term[Schwarzschild の判定条件]{en="Schwarzschild criterion" yomi="しゅわるつしるとのはんていじょうけん"} である。組成勾配を許すと、塊は上昇前の $\mu$ を保つので $\Delta\log\mu=\nabla_\mu\xi/H_P$ が加わる。従って Ledoux 条件は $\nabla_{\rm rad}>\nabla_{\rm ad}+(\varphi/\delta)\nabla_\mu$、$\varphi=(\partial\log\rho/\partial\log\mu)_{P,T}$ となる。

:::example[理想気体]
単原子気体は $P\rho^{-5/3}=$一定だから $\nabla_{\rm ad}=2/5$ である。$\nabla_{\rm rad}=0.5$ なら組成一様で不安定になる。
:::


:::example[組成成層]
同じ勾配でも $\varphi/\delta=1$、$\nabla_\mu=0.2$ なら判定の閾値が $0.6$ に上がり、変位は安定になる。
:::

### 浮力振動数と組成の凍結

上向きの変位を $\xi$ とし、背景の圧力勾配を $\partial_r\log P=-1/H_P$ とする。塊が局所の圧力と均衡しながら断熱かつ同じ組成を保つと、周囲との密度差は
$$
\Delta\log\rho=\frac{\xi}{H_P}
   [\delta(\nabla_{\rm ad}-\nabla)+\varphi\nabla_\mu]
$$
である。ここで $\varphi=(\partial\log\rho/\partial\log\mu)_{P,T}$、$\nabla_\mu=\dd\log\mu/\dd\log P$ である。圧力平衡のため圧力の寄与は塊と背景で消える。一方、温度の差は $-(\nabla_{\rm ad}-\nabla)\xi/H_P$ で、組成の差は $\nabla_\mu\xi/H_P$ である。加速度は $-g\Delta\log\rho$ なので
$$
\ddot\xi+N^2\xi=0,\qquad
N^2=\frac{g}{H_P}[\delta(\nabla_{\rm ad}-\nabla)+\varphi\nabla_\mu].
$$
$N^2>0$ なら :term[浮力振動数]{en="buoyancy frequency" yomi="ふりょくしんどうすう"} $N$ の振動、$N^2<0$ なら成長率 $\sqrt{-N^2}$ の不安定性となる。組成を固定する条件は粒子拡散が変位より遅いことを意味する。

:::example[組成勾配による安定化]
$\delta=\varphi=1$、$\nabla-\nabla_{\rm ad}=0.03$、$\nabla_\mu=0.05$ とする。$N^2=0.02g/H_P>0$ である。例えば $g/H_P=2\times10^{-6}\,\mathrm{s^{-2}}$ なら $N=2\times10^{-4}\,\mathrm{s^{-1}}$、周期は約 $3.14\times10^4\,\mathrm s$ である。組成変化を無視した判定ではこの安定化を取り逃す。
:::

### 内部重力波の局所分散式

安定層の小領域で密度の背景変化を浮力項だけに残す Boussinesq 近似をとる。流体の速度摂動を $\vb u$、浮力を $b$、圧力摂動を基準密度で割ったものを $\pi$ とすると
$$
\partial_t\vb u=-\grad\pi+b\hat{\vb z},\quad
\partial_tb=-N^2 u_z,\quad
\divg\vb u=0.
$$
Fourier 成分に対し運動式を波数 $\vb k$ に垂直に射影する。$\pi$ が消え、鉛直成分は $-i\omega u_z=b(1-k_z^2/k^2)$ となる。浮力式を代入して
$$
\omega^2=N^2\frac{k_\perp^2}{k_\perp^2+k_z^2}
$$
を得る。浮力が戻す運動に水平発散が必要なのは、全体で速度の発散が零だからである。

:::example[伝わる周波数の上限]
$k_z=\sqrt3 k_\perp$ では $\omega=N/2$ である。$N$ より大きい周波数に対しては、この近似式から実数の $k_z$ を得られず、鉛直伝搬解がない。
:::

### 熱と組成の拡散時間が違う層

線形の局所模型で、温度から来る浮力 $b_T$、組成から来る浮力 $b_\mu$、鉛直速度 $w$ をとる。温度と組成の分子拡散係数をそれぞれ $\chi_T,\chi_\mu$、動粘性を $\nu$ とし、波数 $\vb k$ の発展を $e^{st}$ で追う。この小節の $\chi_T$ は熱拡散係数であり、状態方程式の対数微分とは別の量である。非圧縮条件で圧力を消すと
$$
(s+\nu k^2)w=\frac{k_\perp^2}{k^2}(b_T+b_\mu).
$$
安定化する符号を正として断熱浮力を $N_T^2,N_\mu^2$ に分ける。温度と組成は背景勾配を $w$ で運びつつ拡散するので
$$
(s+\chi_T k^2)b_T=-N_T^2w,\qquad
(s+\chi_\mu k^2)b_\mu=-N_\mu^2w .
$$
例えば $N_T^2>0$、$N_\mu^2<0$ のとき、断熱的には $N_T^2+N_\mu^2>0$ なら安定である。それでも移動する塊の温度差が組成差より早く拡散すると、安定化する熱の寄与だけが小さくなりうる。上の三式を連立すると
$$
\begin{aligned}
0={}&(s+\nu k^2)(s+\chi_Tk^2)(s+\chi_\mu k^2)\\
 &+\frac{k_\perp^2}{k^2}
 [N_T^2(s+\chi_\mu k^2)+N_\mu^2(s+\chi_Tk^2)].
\end{aligned}
$$
$k_\perp>0$ で定数項が負なら、左辺は $s=0$ で負、$s\to+\infty$ で正なので正の実根を少なくとも一つ持つ。これは二重拡散による不安定性を判定する十分条件である。この局所模型は一定の背景勾配と弱い摂動に限り、非線形の混合率までは決めない。

:::example[断熱的には安定でも実根が成長する場合]
無次元化して $k^2=k_\perp^2=1$、$\nu=0$、$\chi_T=10$、$\chi_\mu=1$、$N_T^2=2$、$N_\mu^2=-1$ とする。断熱の和は $1$ である。一方、分散多項式は $s^3+11s^2+11s-8$ であり、$s=0$ で $-8$、$s=1$ で $15$ である。この区間で導関数は正なので、成長する実根が一個存在する。
:::

### 圧力平衡に必要な速さ

対流の塊で周囲と同じ圧力を仮定するには、塊を横切る音の時間 $t_s\sim\ell/c_s$ が、移動 $t_{\rm adv}\sim\ell/v$ より短い必要がある。比は Mach 数 $\mathcal M=v/c_s$ である。$v\sim c_s$ になると圧縮性を捨てた局所の浮力計算は信頼しにくい。一方、熱拡散時間は $t_{\rm diff}\sim\ell^2/\chi$ である。この比は $t_{\rm diff}/t_{\rm adv}=\ell v/\chi$ であり、これを :term[Péclet 数]{en="Peclet number" yomi="ぺくれすう"} と呼ぶ。音が速くても熱は遅い場合があり、圧力平衡の条件と断熱移動の条件は独立である。

:::example[圧力は同じでも温度差は残らない塊]
$\ell=10^6\,\mathrm m$、$c_s=10^5\,\mathrm{m\,s^{-1}}$、$v=10\,\mathrm{m\,s^{-1}}$、$\chi=10^8\,\mathrm{m^2\,s^{-1}}$ とする。$t_s=10\,\mathrm s$、$t_{\rm adv}=10^5\,\mathrm s$、$t_{\rm diff}=10^4\,\mathrm s$ である。圧力平衡はよいが、移動中に拡散で熱差が弱まるので、断熱の混合距離近似には補正が要る。
:::

### 組成を消さずに中立条件を求める

熱交換と組成交換が遅い塊で $N^2=0$ の条件は、$H_P>0,\delta>0$ のとき $\nabla=\nabla_{\rm ad}+(\varphi/\delta)\nabla_\mu$ である。輻射だけで全光度を運ぶ状態について、この条件より大きい $\nabla_{\rm rad}$ が必要なら、微小変位の浮力が同じ向きに働く。温度よりも圧力が密度を決める縮退領域では $\delta$ が小さくなるため、両辺を $\delta$ で割った書き方は不安定評価に向かない。$N^2$ に入る組み合わせ $\delta(\nabla_{\rm ad}-\nabla)+\varphi\nabla_\mu$ をそのまま比較すれば、各項が密度へ与える応答から符号を判断できる。

:::example[平均分子量勾配が逆向きの場合]
$\nabla_{\rm ad}-\nabla=0.03$ とし、$\delta=0.2$、$\varphi=1$、$\nabla_\mu=-0.01$ を与える。括弧は $0.006-0.01=-0.004$ である。熱的な差だけなら復元力があるが、上に進む塊が周囲より軽くなる組成差の寄与が勝る。安定化する熱勾配があっても、組成差まで含めて判断する必要がある。
:::

## 混合距離理論

### 有限距離を進む塊

:term[混合距離]{en="mixing length" yomi="こんごうきょり"} $\ell=\alpha_{\rm MLT}H_P$ を、塊が周囲と混ざるまでに進む代表距離とする。効率のよい対流では塊の勾配 $\nabla_e\simeq\nabla_{\rm ad}$ とみなす。上昇距離を $\ell/2$ と置く規約では $\Delta T/T\simeq(\ell/2H_P)(\nabla-\nabla_e)$ となる。密度差に対応する浮力が距離とともに線形に増えると仮定し、仕事を積分することで $v^2\simeq g\delta\,\ell^2(\nabla-\nabla_e)/(4H_P)$ を得る。幾何係数は混合の定義で変わる。

上昇流と下降流のエンタルピー差を平均したフラックスの尺度は $F_{\rm conv}\sim\rho c_P v\Delta T$ であり、$F_{\rm conv}\propto\ell^2(\nabla-\nabla_e)^{3/2}$ となる。$F_{\rm rad}+F_{\rm conv}=L/(4\pi r^2)$ を解いて環境の $\nabla$ を求める。放射損失が速いと塊温度が周囲へ緩和され $\nabla_e\to\nabla$ となるので、表面付近では断熱近似を緩める必要がある。温度の拡散率を $\chi$ とすれば $v\ell/\chi$ が大きいときに断熱近似が適切である。

:::example[勾配の感度]
ほかを固定して $\nabla-\nabla_e$ を四倍にすると速度は二倍、フラックスは八倍となる。
:::


:::example[二つの時間尺度]
$\ell=10^6\,\mathrm m$、$v=100\,\mathrm{m\,s^{-1}}$、$\chi=10^5\,\mathrm{m^2\,s^{-1}}$ なら移動時間は $10^4\,\mathrm s$、拡散時間は $10^7\,\mathrm s$ で、塊はほぼ断熱である。
:::

### 塊からの熱拡散を含める

一定の上昇速度 $v$ を仮定し、局所的に周囲と塊の温度差を $\Delta T$ とする。上昇が作る差は $\dd\Delta T/\dd t=Tv(\nabla-\nabla_{\rm ad})/H_P$ である。塊の外へ熱が拡散する時間を $t_{\rm th}$ とし、熱損失率を $\Delta T/t_{\rm th}$ で近似すると
$$
\dv{\Delta T}{t}+\frac{\Delta T}{t_{\rm th}}=\frac{Tv}{H_P}(\nabla-\nabla_{\rm ad}).
$$
$\Delta T(0)=0$ から時間 $\ell/(2v)$ まで積分する。断熱で予想される $\Delta T_{\rm ad}$ との比は
$$
\frac{\Delta T}{\Delta T_{\rm ad}}=\frac{1-e^{-X}}{X},\qquad X=\frac{\ell}{2v t_{\rm th}} .
$$
$X\ll1$ では比が $1-X/2+\order{X^2}$ となり、$X\gg1$ では $X^{-1}$ まで低下する。$t_{\rm th}\sim\ell^2/\chi$ であれば $X^{-1}$ は移流が拡散より効く指標であり、密度低下に伴う輻射拡散の効率の上昇が表層対流へ影響する。ここでは速度は一定と仮定したため、実際の浮力の仕事と整合させるには、塊の速度もこの温度差から更新する必要がある。

:::example[移動と拡散が競合する]
$X=1$ なら温度差は断熱推定の $1-e^{-1}\simeq0.632$ 倍となる。同じ $v$ なら熱フラックスの比も同じである。
:::

### 混合距離の閉包を単調な方程式として解く

高効率の極限で塊の勾配を $\nabla_{\rm ad}$ と近似し、局所フラックスの輻射成分を $F_{\rm rad}=A\nabla$、対流成分を $F_{\rm conv}=B(\nabla-\nabla_{\rm ad})^{3/2}$ とおく。$A,B>0$ は局所係数である。指数 $3/2$ は、ある距離を動いた塊の温度差が $\nabla-\nabla_{\rm ad}$ に比例し、浮力で得る速度がその平方根に比例するためである。

$F>A\nabla_{\rm ad}$ のとき、$x=\sqrt{\nabla-\nabla_{\rm ad}}$ と置くと $Bx^3+Ax^2=F-A\nabla_{\rm ad}$ である。左辺は $x\ge0$ で連続かつ狭義増加なので、正の解は一意である。上限 $\sqrt{(F-A\nabla_{\rm ad})/A}$ と零で根を挟める。$F\le A\nabla_{\rm ad}$ では対流束を零として $\nabla=F/A$ となる。

:::example[輻射と対流でフラックスを分担する]
共通のフラックス単位で $A=1$、$B=2$、$\nabla_{\rm ad}=0.4$、$F=0.9$ とする。$x=0.5$ なら $2x^3+x^2=0.5$ なので $\nabla=0.65$ となる。輻射が $0.65$、対流が $0.25$ を担い、両者の和は $0.9$ である。輻射だけで全フラックスを運ぶ勾配は $0.9$ であり、対流で実際の勾配は小さくなる。
:::

### 熱輸送と化学混合の時間差

乱流を長さ $\ell$ と速度 $v$ の移動で近似する。三次元で独立な飛程を $n$ 回繰り返すと二乗変位平均は $n\ell^2$、時間は $n\ell/v$ である。従って等方拡散 $\langle r^2\rangle=6Dt$ に対応する尺度は $D\sim v\ell/6$ である。係数の数値は飛程長の分布や相関に依存する。質量分率 $X(r,t)$ に対する混合の保存形は
$$
\left(\pdv{X}{t}\right)_m=\frac{1}{\rho r^2}\partial_r(r^2\rho D\partial_r X)
$$
である。質量殻に沿う微分を使うため、物質境界での流束を零にすると、両辺に $4\pi\rho r^2$ を掛けて積分して全成分質量が保存されることが分かる。

:::example[局所混合の尺度]
$v=100\,\mathrm{m\,s^{-1}}$、$\ell=6\times10^6\,\mathrm m$ を等方飛程模型へ入れると $D=10^8\,\mathrm{m^2\,s^{-1}}$ となる。距離 $10^8\,\mathrm m$ の混合時間の尺度 $L^2/D$ は $10^8\,\mathrm s$ である。混合を瞬時とみなすには、核反応と境界移動の時間がこの値より十分長いことを確認する。
:::

### 熱漏れを含む局所閉包の一つの導出

周囲より高温の塊が圧力平衡で上昇し、長さ $\ell$ の区間で背景勾配を一定とみなす。温度相対差を $q=\Delta T/T$、塊の速さを $v$ とし、熱拡散の緩和時間 $t_{\rm diff}=\ell^2/\chi$ を用いる。移動に伴う断熱差の生成と拡散損失を残した模型は
$$
\dot q=\frac{v}{H_P}(\nabla-\nabla_{\rm ad})-\frac{q}{t_{\rm diff}} .
$$
一定の $v$ で $q(0)=0$ とすると、移動時間 $t_{\rm adv}=\ell/v$ 後の差は
$$
q=\frac{v t_{\rm diff}}{H_P}(\nabla-\nabla_{\rm ad})
 (1-e^{-t_{\rm adv}/t_{\rm diff}}).
$$
断熱差 $q_{\rm ad}=(\ell/H_P)(\nabla-\nabla_{\rm ad})$ と比べ、温度差は $\mathcal P(1-e^{-1/\mathcal P})$ 倍である。$\mathcal P=t_{\rm diff}/t_{\rm adv}$ を大きくすると Taylor 展開で比は $1-1/(2\mathcal P)+\order{\mathcal P^{-2}}$、小さくするとほぼ $\mathcal P$ となる。したがって高速移動では断熱近似、遅い移動では形成した温度差を最終値ごと小さく見積もる必要がある。

この線形緩和模型の $v$ は指定された速さである。自己無撞着に閉じるには浮力の仕事から $v^2\sim g\delta q\ell$ と置き、対流フラックスを $F_{\rm conv}\sim\rho c_PT\,v q$ とする。$\mathcal P\gg1$ なら $q\simeq q_{\rm ad}$ から $F_{\rm conv}\propto (\nabla-\nabla_{\rm ad})^{3/2}$ を得る。$\mathcal P\ll1$ では $q\simeq v t_{\rm diff}(\nabla-\nabla_{\rm ad})/H_P$ であり、非零速度の枝では $v\sim g\delta\ell t_{\rm diff}(\nabla-\nabla_{\rm ad})/H_P$ となる。従って $F_{\rm conv}\propto(\nabla-\nabla_{\rm ad})^3$ の尺度が出る。これは固定形状の局所モデルの指数であり、粘性によるしきい値や幾何係数は含まない。

:::example[拡散と移動時間が等しい場合]
$\mathcal P=1$ なら、同じ速度と背景の過断熱差で $q/q_{\rm ad}=1-e^{-1}\simeq0.632$ である。$\mathcal P=0.1$ なら比は $0.1(1-e^{-10})\simeq0.100$ である。温度差を直接用いる診断では、この比をフラックスの中の $q$ に入れる。自己無撞着な速度を求める場合は $v$ も同時に変わる。
:::

### 局所安定境界を越える塊の運動

中立境界の外側で $N^2>0$ が一定である層を考える。境界からの距離 $\xi$ を微小とし、熱漏れと抗力を無視すると、$\ddot\xi+N^2\xi=0$ である。境界から速度 $v_0$ で出た塊は $\xi=v_0\sin(Nt)/N$ となる。最初の折り返しは $t=\pi/(2N)$、距離は $v_0/N$ である。より一般に $N^2$ が距離で変わり、浮力が局所変位に積分される弱い勾配の模型では、軌道に沿って浮力の仕事と運動エネルギーを比べて転回点を求める。境界の局所安定判定だけでは、有限の慣性を持つ塊の到達距離は決まらない。

:::example[境界外へ入る距離]
$v_0=30\,\mathrm{m\,s^{-1}}$、$N=10^{-3}\,\mathrm{s^{-1}}$ とすると折り返し距離は $3\times10^4\,\mathrm m$ である。この距離が背景の勾配尺度と同程度なら、$N$ 一定の仮定を見直す。塊が折り返せることだけから、その領域すべてが一様に化学混合すると結論することもできない。
:::

### 保存形で化学混合を離散化する

混合を記述するのに質量座標 $m$ を使えば、$\partial_r=(4\pi r^2\rho)\partial_m$ から
$$
(\partial_tX)_m=\partial_m[(4\pi r^2\rho)^2D\,\partial_mX].
$$
$K_m=(4\pi r^2\rho)^2D$ の単位は $\mathrm{kg^2\,s^{-1}}$ であり、$J_m=-K_m\partial_mX$ の単位は $\mathrm{kg\,s^{-1}}$ である。質量幅 $\Delta m_i$ のセルの平均値を $X_i$ とすると
$$
\Delta m_i\dot X_i=J_{i-1/2}-J_{i+1/2}.
$$
この形では共有面の流束は一方で流出、隣で流入であり、離散和でも相殺する。境界で流束を零とすれば、時間更新で $\sum_i\Delta m_i X_i$ が一定に保たれる。したがって差分でセルの勾配を評価する方法を変えても、各面の両側で同じ流束を使う限り全成分質量の収支が閉じる。

:::example[二つの質量殻間での混合]
同じ質量幅 $\Delta m$ の二セルを取る。セル中心の距離も $\Delta m$ と置き、共有面で $K_m$ を一定とする。すると $\dot X_1=K_m(X_2-X_1)/(\Delta m)^2$、$\dot X_2=-\dot X_1$ である。平均 $(X_1+X_2)/2$ は保存し、差は $\exp[-2K_m t/(\Delta m)^2]$ 倍に減る。初期値 $X_1=1,X_2=0$ なら差が $1/2$ になる時刻に $X_1=3/4,X_2=1/4$ である。陽解法で同じ式を進める場合、差を見れば一歩が長すぎて符号が反転するかも調べられる。
:::

## 演習問題

::::exercise{#exr-cosmo-stellar-structure-1 level=basic}
一様密度星で $r=R/\sqrt2$ における $P-P_s$ の中心値に対する比を求めよ。
:::solution
静水圧平衡を積分した結果 $P-P_s\propto1-r^2/R^2$ に代入し $1/2$ を得る。
:::
::::

::::exercise{#exr-cosmo-stellar-structure-2 level=advanced}
一様球の重力エネルギーを求めよ。
:::solution
$r=R(m/M)^{1/3}$ により $\Omega=-(G M^{1/3}/R)\int_0^M m^{2/3}\dd m=-3GM^2/(5R)$ となる。
:::
::::

::::exercise{#exr-cosmo-stellar-structure-3 level=basic}
核発熱が $\epsilon(m)=\epsilon_0(1-m/M)$ の定常星の全光度を求めよ。
:::solution
:term[質量座標]{en="mass coordinate" yomi="しつりょうざひょう"} で損失なしの式を積分すると $L(M)=\epsilon_0[M-M/2]=\epsilon_0M/2$ である。
:::
::::

::::exercise{#exr-cosmo-stellar-structure-4 level=basic}
全電離 He の $\mu$ と電子一個当たりの原子質量単位数 $\mu_e$ を求めよ。
:::solution
質量 $4m_u$ にイオン一個、電子二個があるので $\mu=4/3$、$\mu_e=2$ となる。
:::
::::

::::exercise{#exr-cosmo-stellar-structure-5 level=advanced}
輻射勾配を求める際、なぜ $L/r^2$ の因子が消えるかを示せ。
:::solution
$T'/T$ を $P'/P$ で割ると、輸送式の $L/r^2$ と静水圧式の $Gm/r^2$ の比になる。共通の $r^{-2}$ と $\rho$ が消え、残る幾何に起因する因子は $L/(Gm)$ である。
:::
::::

::::exercise{#exr-cosmo-stellar-structure-6 level=advanced}
理想気体で $\nabla=0.3$、$\gamma=5/3$、$g/H_P=10^{-6}\,\mathrm{s^{-2}}$ の浮力周期を求めよ。
:::solution
$N^2=10^{-6}(0.4-0.3)=10^{-7}\,\mathrm{s^{-2}}$ なので $2\pi/N\simeq1.99\times10^4\,\mathrm s$ である。
:::
::::

::::exercise{#exr-cosmo-stellar-structure-7 level=advanced}
混合距離を二倍にして対流フラックスを保つ。必要な過断熱差の比を求めよ。
:::solution
効率のよい近似で $4(\Delta\nabla_{\rm new}/\Delta\nabla_{\rm old})^{3/2}=1$ なので $2^{-4/3}$ になる。
:::
::::

::::exercise{#exr-cosmo-stellar-structure-8 level=advanced}
表面圧力のある理想気体球で $U$ を $\Omega$ と外圧で表せ。
:::solution
ビリアル関係で $3\int P\dd V=2U$ と置くと、$U=(4\pi R^3 P_s-\Omega)/2$ となる。対流の有無はこの積分関係に影響しない。
:::
::::

::::exercise{#exr-cosmo-stellar-structure-9 level=advanced}
一様密度で真空を外側とする Newton 恒星を考える。$\Gamma_1$ は一定で $5/3$ とする。固有振動において半径全体を同一比率で動かすモードの振動数を求め、この結果を相似変位でのエネルギーの二階微分から確かめよ。

:::solution
相似変位なら $\eta$ は定数である。一様密度の静水圧式は $P'=-\rho GM r/R^3$ なので、:eqref[eq-cosmo-radial-pulsation] から $\omega^2=(3\Gamma_1-4)GM/R^3=GM/R^3$ となる。

次に半径の比を $\lambda$ とすれば、運動エネルギーは $\dot\lambda^2\int r^2\dd m/2=3MR^2\dot\lambda^2/10$ である。単原子理想気体で平衡時の熱エネルギーを $U_0$ とすると、ビリアル式から $U_0=-\Omega_0/2=3GM^2/(10R)$ である。相似変位の $E(\lambda)=U_0\lambda^{-2}+\Omega_0\lambda^{-1}$ は $E''(1)=2U_0$ を満たす。従って運動係数 $3MR^2/5$ で割った小振動の $\omega^2$ は $GM/R^3$ となり、固有方程式の値と一致する。運動エネルギーの重みまで積分することで、ビリアル式の符号判定に加えて周期も得られる。
:::
::::

## 確認テスト

::::quiz
球対称な恒星で局所の重力に寄与する質量はどれか。

- [x] 半径内の質量
- [ ] 全質量
- [ ] 外側の質量

:::explanation
球殻定理により外の球殻の重力が打ち消し合うので、内側の質量だけが寄与する。
:::
::::

::::quiz{answer="0.5" tol="1e-6"}
完全電離した純 H の $\mu$ を求めよ。
:::explanation
陽子一個につき電子一個あるので、粒子数は核数の二倍で $\mu=1/2$ である。
:::
::::

::::quiz{answer="0.4" tol="1e-6"}
単原子理想気体の断熱勾配を求めよ。
:::explanation
断熱で $T\propto\rho^{2/3}$、$P\propto\rho^{5/3}$ なので勾配は $2/5$ である。
:::
::::

::::quiz
放射拡散について正しいものを選べ。

- [x] Rosseland 平均は調和平均である。
- [ ] 平均自由行程が長いほど拡散近似がよい。
- [ ] 中心の光度は一般に有限の正値である。

:::explanation
一次の角度非対称成分が輸送を担い、その係数が $\kappa_\nu^{-1}$ に比例するので、$\partial_TB_\nu$ で重み付けした調和平均を取る。
:::
::::

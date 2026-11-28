---
title: Dirac 方程式の応用と 1 粒子理論の限界
description: 中心力の動径方程式から水素準位を求め、離散対称性と負エネルギー枝の意味を調べる。
notation:
  - symbol: '\kappa'
    meaning: Dirac の中心力量子数
  - symbol: '\mathcal Y_{\kappa m}'
    meaning: スピノル球面調和関数
  - symbol: 'G,F'
    meaning: 大成分と小成分の動径関数
  - symbol: '\psi^{\mathrm c}'
    meaning: 荷電共役したスピノル
---

## この章の目標

- 中心力で Dirac スピノルを二本の動径関数に分離する。
- Coulomb 束縛エネルギーを求め、微細構造の縮退を説明する。
- 粒子と反粒子、空間反転と時間反転の作用を区別する。
- Zitterbewegung、Weyl 場、Majorana 条件から一粒子表示の限界を調べる。

## 中心力場中の Dirac 方程式

:ref[def-aqm-clifford] と同じ計量と Dirac 表現で、$H=c\boldsymbol\alpha\cdot\vb p+\beta mc^2+V(r)I$ を考える。$V$ は静電ポテンシャルエネルギーであり、Dirac 行列の質量項に加える Lorentz スカラーとは結合が異なる。スピンを $\hbar\boldsymbol\Sigma/2$ とすると $\vb J=\vb L+\hbar\boldsymbol\Sigma/2$ は $H$ と可換である。例えば $[L_i,p_j]=i\hbar\epsilon_{ijk}p_k$ と $[\hbar\Sigma_i/2,\alpha_j]=i\hbar\epsilon_{ijk}\alpha_k$ を使うと、$\sum_j\alpha_jp_j$ に対する二つの交換子が反対符号で消える。

二成分角関数を、軌道角運動量 $\ell$ とスピン $1/2$ を結合して作る。最高重みの係数は :eqref[eq-aqm-cg-top] の漸化式で求め、他の $m$ の係数は全角運動量の降下演算子で得る。$\boldsymbol\sigma\cdot\vb L/\hbar=(\vb J^2-\vb L^2-3\hbar^2/4)/\hbar^2$ の固有値は $j=\ell+1/2$ で $\ell$、$j=\ell-1/2$ で $-\ell-1$ である。そこで :term[Dirac の中心力量子数]{en="Dirac angular quantum number" yomi="でぃらっくのちゅうしんりょくりょうしすう"} を $\kappa=-(\ell+1)$ または $+\ell$ とし、角関数 $\mathcal Y_{\kappa m}$ に $(\boldsymbol\sigma\cdot\vb L/\hbar+1)\mathcal Y_{\kappa m}=-\kappa\mathcal Y_{\kappa m}$ を課す。$\boldsymbol\sigma\cdot\hat{\vb r}$ は回転スカラーであり、この角演算子と反交換する。実際、Pauli 積と $[L_i,\hat r_j]=i\hbar\epsilon_{ijk}\hat r_k$ を用いると、$\hat{\vb r}\cdot\vb L=0$ と合わせて $\{\boldsymbol\sigma\cdot\hat{\vb r},\boldsymbol\sigma\cdot\vb L\}=-2\hbar\boldsymbol\sigma\cdot\hat{\vb r}$ となる。従って $\boldsymbol\sigma\cdot\hat{\vb r}$ は同じ $j,m$ の $\kappa$ を $-\kappa$ に替える。二乗が $I$ なので位相を $\boldsymbol\sigma\cdot\hat{\vb r}\mathcal Y_{\kappa m}=-\mathcal Y_{-\kappa m}$ と揃えられる。

定常状態を $\psi=(G\mathcal Y_{\kappa m},iF\mathcal Y_{-\kappa m})^\mathsf T/r$ と置く。$\hat{\vb r}\times\vb p=\vb L/r$ を Pauli 積に用いると

$$
\boldsymbol\sigma\cdot\vb p=-i\hbar\boldsymbol\sigma\cdot\hat{\vb r}
\left(\partial_r-\frac{\boldsymbol\sigma\cdot\vb L}{\hbar r}\right).
$$

これを上成分に作用させて $i\hbar(G'+\kappa G/r)\mathcal Y_{-\kappa m}/r$ を得る。下成分も同様に取り、$M=mc^2$ とすると

$$
G'+\frac{\kappa}{r}G=\frac{M+E-V}{\hbar c}F,\qquad
F'-\frac{\kappa}{r}F=\frac{M-E+V}{\hbar c}G.
\label{eq-aqm-dirac-radial}
$$

規格化は $\int_0^\infty(|G|^2+|F|^2)\dd r=1$ である。四成分のパリティは $\psi(\vb r)\mapsto\beta\psi(-\vb r)$ であり、上下で軌道量子数の偶奇が異なるため、全体では上成分の $(-1)^\ell$ となる。

:::example[同じ $j$ の異なるパリティ]
$j=1/2$ では、$s_{1/2}$ に $\ell=0,\kappa=-1$、$p_{1/2}$ に $\ell=1,\kappa=1$ が対応する。スピノルの下成分の軌道角運動量はそれぞれ $1,0$ となる。
:::

:::example[大成分の低速極限]
$E\simeq M$、$|V|\ll M$ なら $F\simeq\hbar(G'+\kappa G/r)/(2mc)$ となる。第 2 式の左に入れると、$\kappa(\kappa+1)=\ell(\ell+1)$ より $G''-\ell(\ell+1)G/r^2$ が出て、Schrödinger の動径方程式が回復する。
:::

### Dirac の中心力演算子と射影

$\mathcal K=-\beta(\boldsymbol\Sigma\cdot\vb L/\hbar+I)$ とする。上成分の $\mathcal Y_{\kappa m}$ で固有値は $\kappa$、下成分の $\mathcal Y_{-\kappa m}$ でも $\beta$ の負号があるため固有値は $\kappa$ である。$\mathcal K$ は $V,\beta$ と可換で、運動項とは二つの反交換を使う。実際、$\boldsymbol\sigma\cdot\vb L/\hbar+1$ は $\boldsymbol\sigma\cdot\hat{\vb r}$ と反交換し、$\partial_r-\boldsymbol\sigma\cdot\vb L/(\hbar r)$ と可換なので $\boldsymbol\sigma\cdot\vb p$ と反交換する。さらに $\beta$ と $\boldsymbol\alpha\cdot\vb p$ は反交換するため、積の $\mathcal K$ は可換になる。

:::example[同じ全角運動量の二つのチャンネル]
$j=3/2$ では $\ell=1,\kappa=-2$ と $\ell=2,\kappa=+2$ がある。上成分が $p$ の枝では下成分が $d$ であり、上成分が $d$ の枝では下成分が $p$ である。全体のパリティは前者が $-1$、後者が $+1$ となる。ここでの中心力 $V(r)I$ では二チャンネルの行列要素が零である。パリティと $\mathcal K$ のどちらを使っても、同じ $j,m$ の間でも相互作用が混ぜないことを確認できる。
:::

### 動径表示での境界形式と自由な例

二成分の列 $(G,F)$ に対する動径 Hamiltonian は
$$
H_\kappa=\begin{pmatrix}
M+V&\hbar c(-\partial_r+\kappa/r)\\
\hbar c(\partial_r+\kappa/r)&-M+V
\end{pmatrix}.
$$
二つの滑らかな動径関数 $u=(G_u,F_u)$、$v=(G_v,F_v)$ で $\langle u,H_\kappa v\rangle-\langle H_\kappa u,v\rangle$ を計算する。ポテンシャル部分と $\kappa/r$ の部分は消え、部分積分で $\hbar c[\bar F_u G_v-\bar G_u F_v]_0^\infty$ になる。原点と無限遠でこの境界項を消す条件を揃えれば異なるエネルギーの束縛解は直交する。原点で二つの独立解がともに二乗可積分になる特異な力では、単に可積分性だけでこの境界形式が零にはならない。ここで選ぶ境界条件は、正則な有限核または弱結合につながる特定の枝を指す。

:::example[自由な $s_{1/2}$ を境界へ合わせる]
$V=0$、$\kappa=-1$、$E>M$ とする。二本の動径式を消去すると $G''+k^2G=0$、$k=\sqrt{E^2-M^2}/(\hbar c)$ である。原点正則な $G=rj_0(kr)$ は $\sin kr/k$ に等しい。第一式で $G'-G/r=krj_0'(kr)=-krj_1(kr)$ となるので、$F=-(\hbar ck/(M+E))rj_1(kr)$ である。上の $s$ 波に対し下は $p$ 波であり、原点では $G\sim r$、$F\sim r^2$ となる。$E=5M/4$ なら係数の絶対値は $(3/4)/(9/4)=1/3$ で、遠方で両部分波の振幅の比の尺度となる。位相がずれているため、点ごとの $F/G$ を同じ一定値とみなさない。
:::

### 中心力散乱の正枝の外向き波

自由 Hamiltonian を $H_0$、$E>M$ とする。外向き resolvent は
$$
G_0^{(+)}=(E+H_0)[(E+i0)^2-H_0^2]^{-1}
$$
と境界値で書ける。厳密に虚部を残す段では分子も $E+i\varepsilon+H_0$ とし、最後に極限を取る。$H_0^2=M^2+c^2p^2$ なので、逆の二次作用素の位置核は $-e^{ikr}/(4\pi\hbar^2c^2 r)$ である。ここで $k=p/\hbar$、$cp=\sqrt{E^2-M^2}$ である。遠方で $H_0$ をこの外向き位相に作用させると、先頭の $1/r$ には $H_0(\vb p')$ が掛かり、$\vb p'=p\hat{\vb r}$ である。従って分子は $E+H_0(\vb p')=2E P_+(\vb p')$ である。負枝の成分には同じ正のエネルギーの外向き波の留数がないことを確認できる。

静的でスピンに依存しない $V(\vb r)$ に対し、第一 Born では積分の内部のスピノルを $u_s(\vb p)e^{i\vb p\cdot\vb r/\hbar}$ で置き換える。$u^\dagger u=1$ の規約で、出射スピン $s'$ の振幅は
$$
f_{s's}(\vb p',\vb p)
=-\frac{E}{2\pi\hbar^2c^2}\widetilde V((\vb p'-\vb p)/\hbar)\,u_{s'}(\vb p')^\dagger u_s(\vb p).
$$
同じエネルギーで両群速度の大きさが $c^2p/E$ なので、断面積の速度の比は 1 である。入射スピンを二状態で平均し、出射の二状態を足すと、スピン部分は $\frac12\tr[P_+(\vb p')P_+(\vb p)]$ になる。$P_+=(I+H_0/E)/2$ を展開すると $H_0$ 自身のトレースは零、$\tr[H_0(\vb p')H_0(\vb p)]=4(M^2+c^2\vb p'\cdot\vb p)$ である。従ってこの平均は $1-\beta_v^2\sin^2(\theta/2)$、$\beta_v=cp/E$ となる。

:::example[後方散乱のスピン係数を読む]
$V=g/r$ を遮蔽の有限範囲から非前方の極限へ近づける。$\widetilde V(q)=4\pi g/q^2$ を入れると、スピンを除く強度は $g^2/[4p^2v^2\sin^4(\theta/2)]$、$v=c^2p/E$ である。これは $pv$ により非相対論的な $2E_{\rm kin}$ へ近づく。第一 Born の適用には $|g|/(\hbar v)$ が小さいことを要する。例えば $E=2M$、$\theta=\pi$ なら $\beta_v^2=3/4$、スピン平均の係数は $1/4$ である。直角では $\sin^2(\theta/2)=1/2$ なので $5/8$ である。スピンなしの運動学の係数に掛かる補正が角度で異なるため、一つの規格化定数だけで全角分布を直せない。
:::

## 水素原子の厳密解と微細構造

点核の静電力 $V=-\alpha_Z\hbar c/r$ を仮定し、反跳と輻射補正を無視する。ここで $\alpha_Z=Ze^2/(4\pi\epsilon_0\hbar c)$ である。$0<\alpha_Z<1$ とし、原点では弱結合につながる動径関数の枝 $r^\gamma$ を採用し、$-M<E<M$ の束縛域を考える。:eqref[eq-aqm-dirac-radial] の $1/r$ 係数行列は

$$
B=\begin{pmatrix}-\kappa&\alpha_Z\\-\alpha_Z&\kappa\end{pmatrix}.
$$

$B^2=(\kappa^2-\alpha_Z^2)I$ なので採用する動径関数の指数は $\gamma=\sqrt{\kappa^2-\alpha_Z^2}$ である。四成分の波動関数には $1/r$ があるため、動径関数が原点で消えても波動関数が有界とは限らない。無限遠の定数係数行列の固有値は $\pm\lambda$、$\lambda=\sqrt{M^2-E^2}/(\hbar c)$ である。$s=\sqrt{(M-E)/(M+E)}$ を使い $G=w_-+w_+$、$F=s(-w_-+w_+)$ と置くと、減衰・増大の定数部分が対角化される。

### 級数の終端

$\rho=2\lambda r$、$\nu=E\alpha_Z/(\hbar c\lambda)$、$\xi=M\alpha_Z/(\hbar c\lambda)$ とする。変換された動径式は

$$
\begin{aligned}
w_-'&=\left(-\frac12+\frac{\nu}{\rho}\right)w_-+\frac{-\kappa+\xi}{\rho}w_+,\\
w_+'&=\left(\frac12-\frac{\nu}{\rho}\right)w_++\frac{-\kappa-\xi}{\rho}w_- ,
\end{aligned}
$$

ここでプライムは $\rho$ 微分である。例えば変換後の対角成分は $\pm\alpha_Z(1/s-s)/2=\pm\nu$ である。$w_-=\rho^\gamma e^{-\rho/2}u$ とし、第 1 式から $w_+=\rho^\gamma e^{-\rho/2}[\rho u'+(\gamma-\nu)u]/(-\kappa+\xi)$ を得る。これを第 2 式に代入し、$\nu^2-\xi^2=-\alpha_Z^2$ を使うと

$$
\rho u''+(2\gamma+1-\rho)u'+(\nu-\gamma)u=0.
$$

級数の係数比は $(n+\gamma-\nu)/[(n+1)(n+2\gamma+1)]$ である。[Klein–Gordon 方程式](/physics/advanced-quantum-mechanics/klein-gordon/) と同様、終端しない正則解には増大成分があるため $\nu=\gamma+n_r$、$n_r\ge0$ であり、この解は $E>0$ を持つ。

$n_r=0$ なら $u$ は定数、$w_+=0$ が必要なので第 2 式の混合係数 $-\kappa-\xi$ が零でなければならない。従って $\kappa<0$ だけを許す。消去の分母が零の例は $\xi=\kappa>0$ で $\nu=\pm\gamma$ となる。この場合 $w_-\propto e^{-\rho/2}\rho^\nu$ であり、$\nu<0$ なら原点の条件でその係数を零にするので増大する $w_+$ も捨てる。$\nu>0$ なら、$w_+$ の積分因子 $e^{-\rho/2}\rho^\nu$ を掛けた方程式の源は定数倍の $e^{-\rho}\rho^{2\nu-1}$ である。その全積分は非零なので、両端で積分因子付きの $w_+$ を零にできず束縛解は存在しない。$\xi=\kappa$ の場合を除けば、正則な多項式 $u$ と上式で復元した $w_+$ は両方とも原点で $
ho^gamma$ を因子に持ち、無限遠で指数減衰する。従って終端条件は、ここで選んだ境界条件のもとで十分でもある。

従ってエネルギーは

$$
E_{n_r\kappa}
=mc^2\left[1+\frac{\alpha_Z^2}{(n_r+\sqrt{\kappa^2-\alpha_Z^2})^2}\right]^{-1/2},
\qquad n=n_r+|\kappa|.
\label{eq-aqm-dirac-coulomb}
$$

$n$ が同じとき、$|\kappa|=j+1/2$ の値だけでエネルギーが決まる。異なる $m$ は回転対称性で縮退するが、$\kappa$ の符号間の縮退は Coulomb の動径方程式の形に依存する。

:::example[微細構造の展開]
$\gamma=|\kappa|-\alpha_Z^2/(2|\kappa|)+\order{\alpha_Z^4}$ を代入すると

$$
\frac{E}{mc^2}=1-\frac{\alpha_Z^2}{2n^2}
-\frac{\alpha_Z^4}{2n^3}\left(\frac1{j+1/2}-\frac{3}{4n}\right)
+\order{\alpha_Z^6}.
$$

$n=2$ で $j=3/2$ と $1/2$ の差は $mc^2\alpha_Z^4/32$ である。水素で $mc^2\simeq511\,\mathrm{keV}$、$\alpha_Z\simeq1/137$ とすれば約 $4.5\times10^{-5}\,\mathrm{eV}$ である。
:::

:::example[最低状態]
$n_r=0,\kappa=-1$ なら $E=M\sqrt{1-\alpha_Z^2}$、$\lambda=mc\alpha_Z/\hbar$ となる。$w_+=0$ なので $F/G=-s=-\alpha_Z/(1+\sqrt{1-\alpha_Z^2})$ で一定である。$G=N r^\gamma e^{-\lambda r}$ の規格化には $(1+s^2)N^2\Gamma(2\gamma+1)/(2\lambda)^{2\gamma+1}=1$ を用いる。
:::

### 低速展開からの微細構造

厳密な準位と独立に、:eqref[eq-aqm-fw] の三補正が同じ四次項を出すことを確認する。$g=\alpha_Z\hbar c$、$a=\hbar^2/(mg)$、$L=\ell+1/2$、非相対論のエネルギーを $\mathcal E=-g/(2a n^2)$ とする。期待値 $\langle r^{-1}\rangle=1/(a n^2)$ は $\partial_g\mathcal E=-\langle r^{-1}\rangle$ から得る。一般に規格化固有ベクトルの固有値をパラメータで微分すると、状態の微分の寄与は規格化の微分に比例して消え、$\partial\mathcal E=\langle\partial H\rangle$ である。

動径式の $\ell(\ell+1)$ を連続変数として扱う。正則境界を保つなら多項式の終端により $n=n_r+\ell+1$ は同じ式で連続に変わる。$\partial_\ell H=\hbar^2(2\ell+1)/(2m r^2)$ と $\partial_\ell\mathcal E=g/(a n^3)$ を比べ、$\langle r^{-2}\rangle=1/(a^2 n^3 L)$ を得る。$\ell>0$ のときは動径方程式を $u'$ で掛けて積分する。運動項は $-\hbar^2[u'^2]/(4m)=0$ となり、残りを部分積分して $\langle V_{\mathrm{eff}}'\rangle=0$ である。$\ell>0$ では $u\propto r^{\ell+1}$ なので境界項が原点でも消える。$V_{\mathrm{eff}}=-g/r+\hbar^2\ell(\ell+1)/(2mr^2)$ より

$$
\langle r^{-3}\rangle=\frac{1}{a^3n^3\ell(\ell+1)L}\qquad(\ell>0).
$$

$\vb p^2\psi=2m(\mathcal E+g/r)\psi$ のノルム二乗から $\langle p^4\rangle=4m^2(\mathcal E^2+2g\mathcal E\langle r^{-1}\rangle+g^2\langle r^{-2}\rangle)$ である。従って運動補正は $M\alpha_Z^4[3/(8n^4)-1/(2n^3L)]$ になる。$\ell>0$ では Darwin 項は零で、スピン軌道は $M\alpha_Z^4\zeta/[4n^3\ell(\ell+1)L]$ を足す。ここで $\zeta=2\vb L\cdot\vb S/\hbar^2$ の固有値は $j=\ell+1/2$ で $\ell$、$j=\ell-1/2$ で $-\ell-1$ である。どちらを加えても全体の $L^{-1}$ は $(j+1/2)^{-1}$ に置き換わる。

$\ell=0$ では軌道角運動量が零でスピン軌道項も零である。非相対論的 Coulomb 波動関数の $|\psi_{n0}(0)|^2=1/(\pi a^3n^3)$ を使い、Darwin 項の $\nabla^2V=4\pi g\delta^{(3)}(\vb r)$ から $+M\alpha_Z^4/(2n^3)$ を得る。ここでも運動補正と合わせると $j=1/2$ の四次項に一致する。$s$ 波で $\langle r^{-3}\rangle$ を先に計算してスピン軌道の零と掛けると不定な形になるので、$\vb L\psi=0$ を先に使う。

### 弱い磁場で微細構造の各準位を分ける

磁場が十分弱く、Zeeman エネルギーが異なる $j$ の微細構造の差より小さいとする。固定した $n,\ell,j$ で $m_j$ の縮退だけを扱う。点核 Coulomb 模型で同じ $j$ の別の $\ell$ が縮退していても、磁場の一次の演算子はパリティ偶なので両者を混ぜない。低速の電子では磁場の一次の項は $\mu_B(\vb L+2\vb S)\cdot\vb B/\hbar$ である。同じ $j$ の空間に制限すると $\vb S$ は階数 1 のベクトル演算子で、Wigner–Eckart の定理によりその行列要素は $\vb J$ に比例する。比例係数は $\vb J$ との内積を取り、$2\vb S\cdot\vb J=\vb J^2+\vb S^2-\vb L^2$ を使って決まる。異なる動径準位間では、この先頭の演算子が動径座標に依存しないため、同じ $\ell$ での動径直交性が効く。従って
$$
\Delta E=\mu_B B\, g_j m_j,\qquad
g_j=1+\frac{j(j+1)+s(s+1)-\ell(\ell+1)}{2j(j+1)},\qquad s=\frac12
$$
となる。これが :term[Landé の g 因子]{en="Landé g factor" yomi="らんでのじーいんし"} の低速の値である。強い磁場では異なる $j$ も混ざるため、固定した $j$ への射影は使えない。

:::example[$p$ 軌道の二つの Zeeman 間隔]
$\ell=1$ では $j=1/2$ の $g_j$ は $2/3$、$j=3/2$ では $4/3$ である。隣り合う $m_j$ の差は 1 なので、二つの多重項の等間隔はそれぞれ $(2/3)\mu_B B$、$(4/3)\mu_B B$ となる。ここには点核の相対論的な磁気能率の補正や束縛電子の輻射補正を入れていない。
:::

### 最低準位の動径確率

$\kappa=-1$、$n_r=0$ では $E=M\gamma$、$\gamma=\sqrt{1-\alpha_Z^2}$、$\lambda=mc\alpha_Z/\hbar$ である。$G=Cr^\gamma e^{-\lambda r}$ を第一動径式に入れると、$F/G=-(\alpha_Z/(1+\gamma))$ が定数である。実際、$G'+\kappa G/r=[(\gamma-1)/r-\lambda]G$ であり、係数を $M+E+\alpha_Z\hbar c/r$ で割ると $r$ に依存しない。この比を二乗して $\int\dd r (G^2+F^2)=1$ にすると、動径確率密度は
$$
P(r)=\frac{(2\lambda)^{2\gamma+1}}{\Gamma(2\gamma+1)}r^{2\gamma}e^{-2\lambda r}.
$$
これは角積分後の確率密度であり、上成分の頂点だけを見て求めるものではない。$\Gamma$ の積分から $\langle r\rangle=(2\gamma+1)/(2\lambda)$、微分して最大位置は $\gamma/\lambda$ である。

:::example[二つの代表半径を比べる]
弱結合では $\gamma\to1$ なので最大位置は $1/\lambda$、平均は $3/(2\lambda)$ に近づく。例えば $\alpha_Z=3/5$ の点核模型なら $\gamma=4/5$、両者を $1/\lambda$ で割ると $4/5,13/10$ である。小成分の確率は $(1-\gamma)/2=1/10$ となる。同じ模型でも、動径確率と成分の割合は異なる観測量を表す。
:::

### Laguerre 多項式から規格化をまとめる

$n_r\ge1$ とし $a=2\gamma$ と置く。前の微分方程式の多項式を $u=L_{n_r}^{(a)}(\rho)$ と選ぶ。係数表示 $L_n^{(a)}(x)=\sum_{j=0}^n(-1)^j\binom{n+a}{n-j}x^j/j!$ から
$$
\rho u'-n_ru=-(n_r+2\gamma)L_{n_r-1}^{(a)}(\rho)
$$
である。実際、左辺の各係数は $(j-n_r)$ を掛け、$(n_r-j)\binom{n_r+a}{n_r-j}=(n_r+a)\binom{n_r+a-1}{n_r-j-1}$ を使えば右辺になる。したがって
$$
\begin{aligned}
w_-&=C\rho^\gamma e^{-\rho/2}L_{n_r}^{(2\gamma)}(\rho),\\
w_+&=C\rho^\gamma e^{-\rho/2} B\,L_{n_r-1}^{(2\gamma)}(\rho),
\quad B=\frac{n_r+2\gamma}{\kappa-\xi}.
\end{aligned}
$$
積分は $\dd r=\dd\rho/(2\lambda)$ である。$G=w_-+w_+$、$F=s(-w_-+w_+)$ を入れると交差項は Laguerre の重み $\rho^{2\gamma}e^{-\rho}$ で直交する。したがって
$$
C^{-2}=\frac{1+s^2}{2\lambda}
\left[\frac{\Gamma(n_r+2\gamma+1)}{n_r!}
+B^2\frac{\Gamma(n_r+2\gamma)}{(n_r-1)!}\right].
$$
直交性の証明は [Klein–Gordon 方程式の Coulomb 固有状態](/physics/advanced-quantum-mechanics/klein-gordon/) と同じ Rodrigues の部分積分を使う。$n_r=0$ は $w_+=0$ の場合で、最初の項だけが残る。この式では $C$ に $r^\gamma$ ではなく $\rho^\gamma$ の係数を含めているため、最低状態の別表示では $(2\lambda)^\gamma$ を吸収する。

:::example[動径の節が上と下でずれる]
$n_r=1$ なら $L_1^{(2\gamma)}=2\gamma+1-\rho$、$L_0^{(2\gamma)}=1$ である。従って $G$ の多項式の零は $\rho=2\gamma+1+B$、$F$ の零は $\rho=2\gamma+1-B$ で、正の根だけが原点の外の節になる。$\kappa=-1,\alpha_Z\ll1$ なら $\gamma\simeq1,\nu\simeq2,\xi\simeq2$ なので $B\simeq-1$、$G$ の節は $\rho\simeq2$、$F$ は $\rho\simeq4$ である。このように $F$ は $G$ より小さくても波形全体が定数倍とは限らない。$\kappa$ の符号を替えて同じエネルギーにしても、多項式の組み合わせは変わる。
:::

### 質量微分から成分の割合を検査する

点核の $V=-\alpha_Z\hbar c/r$ を質量に依存しないものとして固定する。規格化固有状態で $H(m)\ket{\psi(m)}=E(m)\ket{\psi(m)}$ を微分し左から $\bra\psi$ を掛ければ、$H\ket{\psi'}$ の寄与は右辺の $E\ket{\psi'}$ と相殺する。従って $E'(m)=c^2\langle\beta\rangle$ である。エネルギーが :eqref[eq-aqm-dirac-coulomb] で $m$ に比例するなら $\langle\beta\rangle=E/M$ となり、
$$
\int_0^\infty |F|^2\dd r=\frac{1-E/M}{2}
$$
になる。Laguerre の計算でも $w_-$ と $w_+$ の交差項が積分で消え、下の比は $s^2/(1+s^2)$ になるため一致する。有限核の半径など別の固定長さを入れると $E/M$ が質量に依存し得るので、微分と比の等式はそのまま使えない。

:::example[同じエネルギーの二つの状態を成分で比べる]
$n=2,|\kappa|=1$ の二状態では、$E$ が同じなので下成分の全ノルムは等しい。弱結合では $1-E/M\simeq\alpha_Z^2/8$ だから約 $\alpha_Z^2/16$ である。しかし $s_{1/2}$ の下成分の軌道は $p$、$p_{1/2}$ の下成分は $s$ である。動径積分の一つの割合が等しくても、角度分布や原点近くの重みは同じではない。外のポテンシャルと重ねると符号の縮退を破れるのはこのためである。
:::

### 有限な核の半径が $s$ 準位を動かす

球対称な一様電荷の核半径を $R$ とし、電子のポテンシャルを $-g/r$ に外側で合わせる。核内で $V'(r)=gr/R^3$ なのは、半径 $r$ の球内に含まれる電荷が全体の $(r/R)^3$ だからである。$V(R)=-g/R$ から積分すると、内部は $-g(3-r^2/R^2)/(2R)$ となる。点核との差は
$$
\delta V(r)=g\left(\frac1r-\frac3{2R}+\frac{r^2}{2R^3}\right)\mathbf1_{r<R}
$$
である。非相対論的な水素型 $s$ 状態、$R\ll a_Z$ なら、内部で $|\psi(r)|^2\simeq|\psi(0)|^2$ とする一次摂動が使える。動径積分は
$$
4\pi g|\psi(0)|^2
\left[\frac{R^2}{2}-\frac{R^2}{2}+\frac{R^2}{10}\right]
=\frac{2\pi gR^2}{5}|\psi(0)|^2
$$
である。点核の引力を内部で弱めたので、シフトは正になる。非相対論的な $\ell>0$ では $|\psi|^2$ の先頭が $r^{2\ell}$ のため、積分は $R^{2\ell+2}$ の尺度でさらに小さい。Dirac の高い $Z$ では点核の上・下の成分が $r^{\gamma-1}$ を含むので、$|\psi(0)|^2$ を単純な有限値とする近似を外れる。ここでは弱結合の摂動を検算するための模型である。

:::example[$2s$ と $2p$ の点核での縮退を持ち上げる]
$|\psi_{ns}(0)|^2=1/(\pi a_Z^3 n^3)$ なので、$2s$ の有限核の先頭は $1s$ の $1/8$ である。$R/a_Z=10^{-4}$ とし、エネルギーを $g/a_Z$ で割ると $2s$ のシフトは $(2/5)10^{-8}/8=5\times10^{-10}$ である。同じ $n$ の非相対論的 $2p$ はさらに $(R/a_Z)^2$ の尺度で抑制される。回転対称な有限核でも $\kappa$ の符号間の縮退は保つ必要がない。輻射の Lamb シフトと同じ向きの差を生む場合があるが、原因はポテンシャルの短距離の形が変わることにある。
:::

### 特異な原点で何を境界条件にしたかを追う

原点での係数行列は $B^2=(\kappa^2-\alpha_Z^2)I$ を満たす。$\gamma>0$ のとき二つの固有方向に沿った局所解の主な冪は $r^\gamma,r^{-\gamma}$ である。動径ノルムでは $r^2$ をあらためて掛けず $\dd r$ で積分するので、特異な枝 $r^{-\gamma}$ が平方可積分な条件は $\gamma<1/2$ である。$\gamma=1/2$ では $\int\dd r/r$ が発散する。従って例えば $|\kappa|=1$、$\sqrt3/2<\alpha_Z<1$ では両方が積分でき、選んだ境界条件は「可積分な解」とだけ言って決められない。この章の厳密式は $r^\gamma$ につながる枝を指定した結果である。$\gamma=0$ では二つの固有値が合流し、$B^2=0$ でも $B$ が零でないため、$e^{B\log r}=I+B\log r$ が現れる。$\gamma$ が純虚数となる範囲では $r^{\pm i|\gamma|}$ が原点に向かって位相を回す。有限核を解く場合には原点までこの外部の点核の力を延長しない。自己共役拡張を一般に分類するには非有界な微分演算子の領域論が要るため、ここでは :eqref[eq-aqm-dirac-radial] の原点の二つの冪と境界形式までを使う。一般の作用素の定義域については [非有界作用素](/math/functional-analysis/unbounded-operators/) を参照せよ。

:::example[両方の積分が有限な点核の例]
$|\kappa|=1,\alpha_Z=2\sqrt6/5$ なら $\gamma=1/5$ である。区間 $(0,R)$ の $r^{-2\gamma}$ の積分は $(5/3)R^{3/5}$、正則枝の $r^{2\gamma}$ は $(5/7)R^{7/5}$ であり、どちらも有限である。このためノルムのみを比べても係数の線形結合を決められない。境界の固定による固有値の式を使うときは、指数の正負と核の模型を明記する。
:::

## 空孔理論と陽電子

自由 Dirac Hamiltonian の負枝は大きな運動量でいくらでも低い。電子一個が負枝へ落ちることを防ぐ歴史的な模型では、負枝をすべてフェルミオンで占有する。反交換関係により同じ状態に二個入れられない。体積と運動量を有限にして有限個の準位として考え、満たした基準からの差を取ると、準位の電子を一個除く操作はエネルギー $+E$、運動量 $-\vb p$、電荷 $+e$ を加える。これが :term[空孔]{en="hole" yomi="くうこう"} による陽電子の描像である。

場の演算子では正周波数係数に $a$、負周波数係数に $d^\dagger$ を割り当てる。負枝のエネルギー項 $-E d d^\dagger$ は反交換により $E d^\dagger d-E$ になる。電荷も $-e(a^\dagger a+d d^\dagger)$ から真空定数を引いて $-e a^\dagger a+e d^\dagger d$ となる。従って場の真空から陽電子を作ると正エネルギーを持つ。この書き方では占有された無限個の電子を直接数えず、規格化した真空と励起を使う。

:::example[対の静止エネルギー]
電子・陽電子を重心系で静止させるなら必要な静止エネルギーは $2mc^2\simeq1.022\,\mathrm{MeV}$ である。真空中の単一光子の四元運動量の二乗は零だが、二粒子の合計は少なくとも $(2mc)^2$ なので自由空間の一光子だけでは生成できない。
:::

:::example[自由対の電荷]
$a^\dagger d^\dagger\ket0$ の全電荷は $-e+e=0$ である。一方、Hamiltonian の二つの正項によりエネルギーは $E_a+E_d$ である。電荷零をエネルギー零と取り違えない。
:::

### 負枝のモードを粒子と空孔に読み替える

自由場を箱に入れ、運動量とスピンの有限個のモードだけで代数を確かめる。正枝の消滅を $a_s$、負枝を占める粒子の消滅を $b_s$ としたとき、素朴なエネルギーは $\sum_s E_s(a_s^\dagger a_s-b_s^\dagger b_s)$ である。負枝をすべて埋めた参照状態では、$b_s$ は空孔を作るので $d_s^\dagger=b_s$、$d_s=b_s^\dagger$ と置ける。反交換関係により $-E_s b_s^\dagger b_s=E_s d_s^\dagger d_s-E_s$ となる。参照エネルギーを引くと $\sum_s E_s(a_s^\dagger a_s+d_s^\dagger d_s)$ になる。電荷も $q(a_s^\dagger a_s+b_s^\dagger b_s)$ から参照値を引けば $q(a_s^\dagger a_s-d_s^\dagger d_s)$ となり、空孔に逆符号の電荷が付く。

:::example[一つの負枝を空けて正枝を占める]
同じ正エネルギーの絶対値 $E$ の二モードを取る。負枝だけ占めた参照状態のエネルギーは $-E$、電荷は $q$ である。負枝を空けて正枝を占めた状態のエネルギーは $+E$、電荷は $q$ である。差はエネルギー $2E$、電荷零であり、粒子一個と反粒子一個の励起と読める。連続無限個のモードではこの参照定数をそのまま有限な和として扱わず、生成・消滅演算子の正規順序を定義する。
:::

### 正エネルギーの二体状態でもしきい値の拘束が残る

空孔が正のエネルギーを持つように数え直しても、粒子と反粒子を作る相互作用には全電荷、エネルギー、運動量の保存が必要である。たとえば二つの質量がともに $m$ の励起の全運動量を $\vb P$ とする。固定 $\vb P$ において最小の全エネルギーは、二つの粒子が同じ速度を持つ場合である。実際、$E(p)=\sqrt{M^2+c^2p^2}$ は厳密凸であり、$\tfrac12[E(\vb p_1)+E(\vb p_2)]\ge E((\vb p_1+\vb p_2)/2)$ である。Hessian の直交方向の固有値は $c^2/E$、平行方向は $m^2c^6/E^3$ なので両方正である。従って $E_{\min}(\vb P)=\sqrt{(2M)^2+c^2P^2}$ になる。

:::example[全運動量を持つ対のエネルギー]
$P=2mc$ の対では最小エネルギーは $2\sqrt2\,M$ であり、各粒子に $mc$ を分ける。この時の運動エネルギーの合計は $2(\sqrt2-1)M$ である。一方、一方に全 $2mc$ を渡して片方を静止させると $(\sqrt5+1)M$ となり、最小値より大きい。生成される二粒子のエネルギーを常に $2M$ と置けるのは全運動量が零の系での閾値に限る。
:::

### 真空からの励起と外場による混合を分ける

有限のモードで正枝の粒子を $a$、負枝の空孔を $d$ とする。二次の一粒子演算子の行列要素が元の正負を結ぶ項は、$a^\dagger b+b^\dagger a$ の形である。空孔表示 $b=d^\dagger$ に替えれば $a^\dagger d^\dagger+da$ になる。これは電荷を変えず粒子と反粒子を一つずつ作る演算子であり、元の粒子の一個の移動と同じ代数から出る。係数が静的か時間依存か、初期と終状態が自由か外場の固有状態かにより、どの準備を真空と呼ぶかが変わる。自由基底でこの項が非零なのは結合が存在することを示すが、エネルギー保存を満たさずに定常な真空が粒子を出し続けることを意味しない。

:::example[外場に対応する二モードの最小模型]
$H=E(a^\dagger a+d^\dagger d)+g(a^\dagger d^\dagger+da)$、$g$ は実とする。偶の $(\ket0,a^\dagger d^\dagger\ket0)$ のブロックは $\begin{pmatrix}0&g\\g&2E\end{pmatrix}$ である。$E>0$ で $|g|\ll E$ ならその最低固有値は $E-\sqrt{E^2+g^2}\simeq-g^2/(2E)$、対の確率は $[1-E/\sqrt{E^2+g^2}]/2\simeq g^2/(4E^2)$ である。静的な真空をこの最低状態として準備すれば、その対成分は時間で増えない。自由の裸の $\ket0$ から結合を急に入れた場合は別で、この二次元の回転で振動的な移行が生じる。外場を戻した後の自由粒子の数を調べるなら、切り替えの時間依存も問題の一部になる。
:::

## 荷電共役・パリティ・時間反転

Dirac 表現の波動関数で次の規約を取る。$K$ は成分の複素共役である。荷電共役は $\psi^{\rm c}=B\psi^*$、$B=i\gamma^2$。Pauli 行列の実虚に注意すると $B\gamma^{\mu*}B^{-1}=-\gamma^\mu$ である。電荷 $q$ の式 $(i\hbar\gamma^\mu\partial_\mu-q\gamma^\mu A_\mu-mc)\psi=0$ を共役し $B$ を掛ければ、微分項の二つの負号が相殺し、電荷が $-q$ の式を得る。この :term[荷電共役]{en="charge conjugation" yomi="かでんきょうやく"} は一粒子波動関数では反線形である。場の Fock 空間で粒子演算子と反粒子演算子を交換する操作とは区別する。

パリティでは $\psi_P(t,\vb r)=\beta\psi(t,-\vb r)$ を取る。$\beta\gamma^0\beta=\gamma^0$、$\beta\gamma^i\beta=-\gamma^i$ が空間微分の符号を補う。外場も $\Phi_P(t,\vb r)=\Phi(t,-\vb r)$、$\vb A_P(t,\vb r)=-\vb A(t,-\vb r)$ に移す。時間反転は $\psi_T(t,\vb r)=R\psi(-t,\vb r)^*$、$R=-i\Sigma_y$ と置く。$R\beta^*R^{-1}=\beta$、$R\alpha_i^*R^{-1}=-\alpha_i$ であり、$\vb p^*=-\vb p$ なので自由 $H$ は保たれる。時間微分の符号は $i$ の共役で補われる。電磁場は $\Phi_T(t,\vb r)=\Phi(-t,\vb r)$、$\vb A_T(t,\vb r)=-\vb A(-t,\vb r)$ である。この規約で $R$ は実で $R^2=-I$ なので :ref[thm-aqm-kramers] の条件に合う。

:::example[平面波を荷電共役する]
$u e^{-i(Et-\vb p\cdot\vb r)/\hbar}$ は $B u^* e^{i(Et-\vb p\cdot\vb r)/\hbar}$ になる。自由 Hamiltonian の負枝として見れば四元波数の符号が逆になる。場の演算子ではこのモードが正のエネルギーの反粒子を作る。
:::

:::example[一定外場での対称性]
一様な $\vb B=B\hat{\vb z}$ は軸性ベクトルなので空間反転で保たれるが、時間反転では $\vb B\mapsto-\vb B$ となる。外場を固定した一粒子 Hamiltonian に時間反転を適用すると、一般に同じ Hamiltonian に戻らない。
:::

### 電荷共役と電流の符号

一粒子の $\psi^{\rm c}=B\psi^*$ では、$B^\dagger B=I$ なので密度 $\psi^\dagger\psi$ は同じである。$\alpha_i=\gamma^0\gamma^i$ に対して $B\alpha_i^*B^{-1}=\alpha_i$ であることから、空間成分の速度密度も同じになる。物理的電荷カレントでは同時に $q\mapsto-q$ なので全成分の符号が替わる。一方スカラー $\bar\psi\psi$ では $B\beta^*B^{-1}=-\beta$ のため、可換な一粒子スピノルの値は符号が逆になる。場の積ではフェルミオン演算子の順序の符号も現れるため、この可換成分の計算をそのまま場の双一次形式へ移さない。

:::example[静止スピノルの三変換]
Dirac 表現で $u=(1,0,0,0)^\mathsf T$ とする。$\beta u=u$ であり、$R u^*=(0,1,0,0)^\mathsf T$ は同じ枝の逆スピンである。$i\gamma^2u^*=(0,0,0,1)^\mathsf T$ は $H=\beta M$ の負枝である。いずれも密度は 1 だが、荷電共役の前後のみ $\bar uu$ が $+1$ から $-1$ へ替わる。時間反転を二度行うと $-u$、荷電共役を二度行うと $u$ となる。
:::

### 離散変換を行列と引数に分ける

ここでの一粒子規約では $P=\beta$ と空間引数の反転、$T=R K$ と時間引数の反転、$C=B K$ を使う。$R=-i\diag(\sigma_y,\sigma_y)$、$B=i\gamma^2$ は Dirac 表現で実である。行列 $B$ は $\beta$ と反交換、$R$ は $\beta$ と可換である。そのため位相を固定したこの規約では $CP=-PC$、$TP=PT$ となる。$R$ と $B$ はブロック形で $\sigma_y$ だけを含んで可換なので $CT=TC$ である。パリティの演算子に位相を選び直すと反線形演算子との交換で位相の共役が現れるため、これらの符号をあらゆる規約に独立なものと扱わない。$C$ と $T$ はそれぞれ反線形なので合成には二回の共役があり、$CPT$ をこの一粒子の古典波動関数に適用した場合は線形になる。場の量子論における CPT の主張には局所性など別の仮定があり、この行列の計算で一般の相互作用場の定理を証明したことにはならない。

:::example[反粒子を作る波の内在パリティ]
$p=0$ で $u$ を正枝の固有ベクトルとする。$\beta u=+u$ であり、$\beta(Bu^*)=-Bu^*$ なので共役波のパリティ固有値は逆である。同じ箱で静止して座標反転の位相がない場合でも、上下の成分で符号が替わることを確かめられる。$R K$ では $\beta(Ru^*)=+Ru^*$ なのでエネルギー枝もパリティも保つ。空間部が異なる場合は軌道の $(-1)^\ell$ もこの内在成分に掛ける。
:::

### 外場と観測量の変換を別々に照合する

パリティで四元カレントは $(j^0,\vb j)(t,-\vb r)$ の時間成分を保ち空間成分を反転する。同じように時間反転では転置でなく複素共役が行列の成分に掛かり、$R\alpha_i^*R^\dagger=-\alpha_i$ により $\vb j$ を反転する。電磁場の Lorentz 力 $q(\vb E+\vb v\times\vb B)$ と比べると、空間反転では $\vb E$ が奇、$\vb B$ が偶になり、時間反転では $\vb E$ が偶、$\vb B$ が奇になる。これは $\vb E=-\grad\Phi-\dot{\vb A}$、$\vb B=\grad\times\vb A$ に前の外場変換を代入しても確かめられる。$\Phi$ の変換を固定したまま磁場だけの符号を議論すると、時間依存の誘導電場に対する整合性を落とすことがある。

:::example[交差した電場と磁場の対称性]
固定時刻で空間一様な $\vb E=E\hat x$、$\vb B=B\hat z$ を考える。空間反転だけでは $\vb E$ が逆、時間反転だけでは $\vb B$ が逆なので、両方非零ならどちらも単独の対称性ではない。一粒子方程式で $C$ を行うと $q\mapsto-q$ の外場を保った式に移る。同じ電荷の式として比較するには外場を同時に $A_\mu\mapsto-A_\mu$ とみなす。すると $P,T$ の積で $\vb E,\vb B$ がともに逆になるものを、$C$ に対応する外場の反転で元に戻せる。座標引数も反転するため、非一様な場ではこの成分の符号だけで結論は決まらない。
:::

## Zitterbewegung

自由の速度演算子 $c\alpha_i$ は $H$ と必ずしも可換でない。:eqref[eq-aqm-dirac-alpha] から $\{H,\alpha_i\}=2cp_i$ なので $D_i=\alpha_i-cp_iH^{-1}$ は $H$ と反交換する。$e^{iHt/\hbar}D_i=D_i e^{-iHt/\hbar}$ が指数級数から従うため $\alpha_i(t)=cp_iH^{-1}+D_i e^{-2iHt/\hbar}$ である。$m>0$ なら $H^{-1}$ が有界なので、Heisenberg 方程式 $\dot x_i=c\alpha_i(t)$ を積分して

$$
x_i(t)=x_i(0)+c^2p_iH^{-1}t
+\frac{i\hbar c}{2}D_iH^{-1}(e^{-2iHt/\hbar}-I).
\label{eq-aqm-zitter}
$$

この振動を :term[Zitterbewegung]{en="Zitterbewegung" yomi="つぃったーべゔぇーぐんぐ"} と呼ぶ。$D_i$ は正負を入れ替えるため、振動部分の期待値は両方の枝を含む状態で生じる。例えば $H\ket u=E\ket u$ なら $H D_i\ket u=-E D_i\ket u$ で、正枝に射影した $D_i$ は零である。

:::example[静止時の尺度]
$p=0$ で $\omega_{\rm Z}=2mc^2/\hbar$、座標の振動係数の尺度は $\hbar/(2mc)$ である。電子ならおよそ $1.55\times10^{21}\,\mathrm{s^{-1}}$ と $1.93\times10^{-13}\,\mathrm m$ である。この周波数を一粒子が放射する古典的軌道の周波数としては解釈しない。
:::

:::example[枝を混ぜた状態]
$\vb p=0$、$z$ 上向きの二成分 $\chi$ で $\psi(0)=(\chi,\chi)^\mathsf T/\sqrt2$ と取る。上と下の位相は $e^{\mp iMt/\hbar}$ なので $\langle c\alpha_z(t)\rangle=c\cos(2Mt/\hbar)$ となる。時間積分は $\langle x_z(t)-x_z(0)\rangle=(\hbar c/(2M))\sin(2Mt/\hbar)$ を与える。これは箱規格化した内部スピノルの速度から得る変位の式である。
:::

### 波束の干渉と運動量の分布

固定運動量で $\psi(p)=u_+(p)+u_-(p)$ と枝を分ける。速度の振動部分の期待値は
$$
2c\,\Re\int\dd^3p\, e^{-2iE_pt/\hbar}
u_-(p)^\dagger D_i(p)u_+(p)
$$
である。$D_i$ は枝を入れ替え、速度演算子は運動量を変えないのでこの形になる。有限な幅の波束では $E_p$ が異なる成分が異なる角速度で回る。狭い分布を非零の中心 $p_0$ に置き、一次元で係数が分布内でほぼ一定とする。$E_p\simeq E_0+v_g(p-p_0)$、運動量の確率密度が標準偏差 $\Delta p$ の Gaussian なら、振動の包絡は $e^{-2(v_g\Delta p\,t/\hbar)^2}$ となる。これは Gaussian の特性関数であり、二次のエネルギー補正が小さい時刻に使う。

:::example[二つの近いエネルギーのうなり]
運動量が異なる二モードで振動速度の寄与がそれぞれ $A\cos(2E_1t/\hbar)$、$A\cos(2E_2t/\hbar)$ とする。和は $2A\cos((E_1+E_2)t/\hbar)\cos((E_1-E_2)t/\hbar)$ になる。最初に包絡が零になる時刻は $\pi\hbar/(2|E_1-E_2|)$ である。各モードでの正負の干渉が消えたわけではなく、平均に寄与する二つの振動が相殺した時刻である。
:::

### 位置演算子の表現と振動のブロック

自由の Foldy–Wouthuysen 写像 $U(\vb p)$ で $H$ は $\beta E$ になる。運動量表示で元の位置 $x_i=i\hbar\partial_{p_i}$ をこの表示へ移すと
$$
U x_iU^\dagger=i\hbar\partial_{p_i}+i\hbar U(\partial_{p_i}U^\dagger)
$$
である。第二項は $U U^\dagger=I$ の微分により Hermite である。運動量微分だけの時間発展は、$e^{i\beta Et/\hbar}i\hbar\partial_{p_i}e^{-i\beta Et/\hbar}$ から $i\hbar\partial_{p_i}+\beta(c^2p_i/E)t$ となる。一方、第二項のうち $\beta$ に奇なブロックは共役で $e^{\pm2iEt/\hbar}$ が掛かる。従ってどの位置演算子を測るかまで指定し、状態だけを新しい表示に移して演算子を据え置かない。正枝のみに制限すれば奇ブロックは消える。

:::example[低運動量で振動係数の長さを確かめる]
$p=0$ 近くで $U=I+\beta\alpha_i p_i/(2mc)+\cdots$ と展開すると $i\hbar U\partial_{p_i}U^\dagger=-i\hbar\beta\alpha_i/(2mc)$ になる。$\beta\alpha_i$ は上と下を結ぶのでこの長さは Compton 波長の半分である。$\beta$ に奇な行列を固定したスピンの正負二枝に制限すれば、回転する非対角の角周波数は $2M/\hbar$ である。正枝だけの表示で明示的な高速項がないことと、元の位置演算子で二枝の干渉があることは両立する。
:::

### 静止近傍の波束では包絡が冪で弱まる

前の非零の中心運動量の展開と違い、$\vb p=0$ を中心とする狭い波束では $v_g=0$ であり二次から位相が変わる。$\psi$ の正負干渉の係数が運動量分布内でほぼ一定と仮定し、$E_p=M+p^2/(2m)+\cdots$ を使う。三軸の運動量がそれぞれ分散 $s_p^2$ の独立 Gaussian なら
$$
\int\frac{\dd^3p}{(2\pi s_p^2)^{3/2}}
 e^{-p^2/(2s_p^2)}e^{-ip^2t/(m\hbar)}
=\left(1+\frac{2is_p^2t}{m\hbar}\right)^{-3/2}.
$$
これは各軸の複素 Gaussian の積で、$t=0$ で 1 に揃えた枝を使う。従って $e^{-2iMt/\hbar}$ を除いた絶対値は $[1+(2s_p^2t/(m\hbar))^2]^{-3/4}$ である。二次相での計算なので、典型波数 $p\sim s_p$ の四次の残りが $t$ と掛かって小さい範囲でこの式を使う。中心が非零なら一次の幅による Gaussian 型の包絡が先に現れ得る。

:::example[運動量幅の二乗が時間尺度を決める]
$t=m\hbar/(2s_p^2)$ では絶対値は $2^{-3/4}\simeq0.595$ である。$s_p$ を半分にして同じ時刻を比べれば括弧の二乗変数が $1/16$ になり、絶対値は $(17/16)^{-3/4}\simeq0.956$ となる。同じ包絡の低下を得るには時間を四倍にする必要がある。これは環境へ情報を逃がす緩和ではなく、異なる運動量の位相を足した期待値の平均化である。
:::

## Weyl 方程式と Majorana 方程式

Weyl 表現で $\gamma^0=\begin{pmatrix}0&I\\I&0\end{pmatrix}$、$\gamma^i=\begin{pmatrix}0&\sigma_i\\-\sigma_i&0\end{pmatrix}$ とする。$\gamma^5=\diag(-I,I)$ なので $\psi=(\psi_L,\psi_R)^\mathsf T$ と書ける。Hamiltonian の運動項は $c\diag(-\boldsymbol\sigma\cdot\vb p,+\boldsymbol\sigma\cdot\vb p)$、質量項は $mc^2$ を上下の非対角に持つ。$m=0$ では独立な二つの :term[Weyl 方程式]{en="Weyl equation" yomi="わいるほうていしき"} となる。正のエネルギーでは $\psi_R$ が正、$\psi_L$ が負のヘリシティを持つ。空間反転は $\gamma^0$ で左右を交換する。

中性場について $\psi^{\rm c}=\psi$ を課すと、この表現では $\psi_R=-i\sigma_y\psi_L^*$ になる。上の Dirac 方程式に入れると

$$
i\hbar\partial_t\psi_L=-c\boldsymbol\sigma\cdot\vb p\,\psi_L-i mc^2\sigma_y\psi_L^* .
$$

これは :term[Majorana 方程式]{en="Majorana equation" yomi="まよらなほうていしき"} の二成分表示である。下の式はこの複素共役から従う。実際、$-i\sigma_y\sigma_i^*(i\sigma_y)=-\sigma_i$ と $\vb p^*=-\vb p$ がある。複素位相 $\psi\mapsto e^{i\theta}\psi$ に対し $\psi^{\rm c}\mapsto e^{-i\theta}\psi^{\rm c}$ なので、同じ実条件を保つ位相は $e^{i\theta}=\pm1$ に限る。通常の連続電荷対称性を持つ荷電粒子にはこの条件を課せない。

:::example[質量零の右成分]
$\vb p=p\hat{\vb z}$、$p>0$ では $H_R=cp\sigma_z$ なので正エネルギーの Weyl モードは $(1,0)^\mathsf T$、負は $(0,1)^\mathsf T$ である。左成分では正負の割り当てが反対になる。
:::

:::example[自己共役条件の成分]
$\psi_L=(1,i)^\mathsf T/\sqrt2$ なら $-i\sigma_y=\begin{pmatrix}0&-1\\1&0\end{pmatrix}$ なので $\psi_R=(i,1)^\mathsf T/\sqrt2$ となる。全スピノルを規格化するならさらに全体に $1/\sqrt2$ を掛ける。
:::

### 質量があるとカイラリティとヘリシティが分かれる

Weyl 表現で $p>0$、$\boldsymbol\sigma\cdot\hat{\vb p}$ の固有値を $s=\pm1$ とする。同じスピン固有ベクトルの左右成分の Hamiltonian は
$$
H_s=\begin{pmatrix}-scp&M\\M&scp\end{pmatrix}
$$
である。$H_s^2=E^2I$ を使うと正枝の射影は $(I+H_s/E)/2$ になる。その対角成分は正枝での左・右成分の重みなので $w_L=(1-scp/E)/2$、$w_R=(1+scp/E)/2$ である。質量零で $E=cp$ なら片方が消える。有限質量ではスピンの向きは一定でも、左右の両方が残る。

:::example[右ヘリシティに含まれる左成分]
$s=+1$、$p=3mc/4$ なら $E=5M/4$ である。従って $w_R=4/5$、$w_L=1/5$ となり、カイラリティの期待値 $w_R-w_L$ は $3/5$ である。速度 $3c/5$ で共に動く系へ boost すると運動量は零になる。そこでヘリシティの軸は定義できないが、左右のノルムは等しくなる。
:::

### Weyl 成分の保存の破れを方程式から読む

Weyl 基底で $i\hbar\dot\psi_L=-c\boldsymbol\sigma\cdot\vb p\,\psi_L+M\psi_R$、$i\hbar\dot\psi_R=c\boldsymbol\sigma\cdot\vb p\,\psi_R+M\psi_L$ である。空間積分で運動量項に境界流束がないとすると
$$
\dv{}{t}\|\psi_L\|^2=\frac{2M}{\hbar}\Im\inner{\psi_L}{\psi_R},\qquad
\dv{}{t}\|\psi_R\|^2=-\frac{2M}{\hbar}\Im\inner{\psi_L}{\psi_R}.
$$
質量が零なら左右それぞれが保存される。有限質量では総ノルムだけがこの式から保存され、同じヘリシティ内で左右が移り得る。固定した正エネルギーの一固有状態では左右の比は一定であるため、時間に依らないカイラリティ期待値が出ても、演算子が Hamiltonian と可換であることを意味しない。

:::example[左だけを準備した自由モード]
静止で初期を $(\chi,0)$ として $\|\chi\|=1$ とする。質量項は左右を交換し二乗が $M^2I$ なので、発展後は $(\cos(Mt/\hbar)\chi,-i\sin(Mt/\hbar)\chi)$ になる。時間 $\pi\hbar/(2M)$ で右だけになる。この準備は正枝の固有状態ではなく、静止の左右を等しく含む正・負枝の両方を使っている。実際、両枝の位相差 $2Mt/\hbar$ が左右の重みの変化を作る。
:::

### Majorana の左右の関係が Lorentz 変換で保たれる理由

Weyl 表現の回転では左右とも $U=e^{-i\theta\boldsymbol\sigma\cdot\hat{\vb n}/2}$、同じ方向の boost では $S_L=e^{-\zeta\boldsymbol\sigma\cdot\hat{\vb n}/2}$ と $S_R=e^{+\zeta\boldsymbol\sigma\cdot\hat{\vb n}/2}$ である。$\epsilon=-i\sigma_y$ とすると $\epsilon\sigma_i^*\epsilon^{-1}=-\sigma_i$ である。従って回転では複素共役の $i$ の符号と Pauli の負号が相殺し、$\epsilon U^*\epsilon^{-1}=U$ である。boost では指数の係数が実なので Pauli の負号だけが残り、$\epsilon S_L^*\epsilon^{-1}=S_R$ である。変換の積にも同じ関係が伝わる。よって $\psi_R=\epsilon\psi_L^*$ を一つの系で課せば、固有順時の Lorentz 変換後も同じ式が成り立つ。右・左を独立に持つ Dirac 場と異なり、この条件では左のデータから右が決まる。

:::example[一方向の boost で二つの成分を追う]
左のスピンを $(1,0)^\mathsf T$、右は $\epsilon(1,0)^\mathsf T=(0,1)^\mathsf T$ とする。$\hat z$ の rapidity $\zeta$ では左の上は $e^{-\zeta/2}$、右の下も $e^{-\zeta/2}$ を得る。右の下に $e^{+\zeta/2}$ を掛けると左との共役条件を失う。Weyl の二成分への boost は通常のノルムを保つ回転ではなく、二つの表現で符号が異なることがこの照合で見える。
:::

## 1 粒子理論の限界と場の理論の必要性

一粒子理論で正枝だけを選ぶと、自由な時間発展について閉じる。しかし空間依存の摂動 $V(\vb r)$ は運動量を変え、運動量に依存する正負射影と一般に可換でない。:eqref[eq-aqm-dirac-u] の上・下成分の関係を比較すれば、静止の正枝スピノルを運動量だけ変えたものに、新しい運動量の負枝成分が含まれることが分かる。時間変化する外場では、場の表示でこのような混合を粒子・反粒子の生成消滅を含めて記述する。

自由空間の正エネルギーに射影したデルタ局在は、射影 $P_+(\vb p)$ が $\vb p$ に依存するため位置空間でデルタのままではない。例えば $P_+(\vb p)$ は $1/\sqrt{p^2c^2+m^2c^4}$ を含むので有限階の局所微分作用素ではない。場の理論では局所場の因果性は正負周波数の両方を含む交換子または反交換子で実装する。電磁場との相互作用や対生成の遷移確率を一貫して扱うには、Fock 空間の異なる粒子数を結ぶ Hamiltonian が必要である。相対論的な場の因果性とスピン・統計の証明は、局所場の公理とスペクトル条件を要するので [場の量子論](/physics/quantum-field-theory/) に委ねる。

:::example[正枝への射影の失敗]
初期スピノルを $(\chi,0)^\mathsf T$ とし、空間位相 $e^{ipz/\hbar}$ を掛ける。$p\ne0$ の正射影の重みは $\mel{\chi;0}{(I+H/E)/2}{\chi;0}=(1+M/E)/2$ である。$E=2M$ では正枝 $3/4$、負枝 $1/4$ となる。この計算は位相を掛ける操作の枝成分を測ったもので、外場による対の生成率は別に力学を指定して求める。
:::

:::example[外場が反跳を受け持つ場合]
重い核の近くで光子が対を作るときは、核も四元運動量の保存に入る。核質量を $M_N$ とすると静止核への光子の全不変質量二乗は $M_N^2c^4+2M_Nc^2 E_\gamma$ であり、閾値で $(M_N+2m)^2c^4$ になる。従って $E_\gamma=2mc^2(1+m/M_N)$ で、重い核では静止エネルギー二個分に近い。
:::

### 自由な一光子と二光子の閾値

一光子の四元運動量の二乗は零である。正エネルギーの質量 $m$ の二粒子について、片方の静止系を取ると $p_1\cdot p_2=mE_2\ge m^2c^2$ なので $(p_1+p_2)^2\ge4m^2c^2$ となる。従って真空中の単一光子から、運動量を受け取る別の物体なしにこの二粒子だけを作ることはできない。二光子のエネルギーを $E_1,E_2$、なす角を $\theta$ とすると、全四元運動量の二乗に $c^2$ を掛けた値は $2E_1E_2(1-\cos\theta)$ である。この値が $(2mc^2)^2$ 以上なら重心系で許容される。

:::example[交差角で閾値が変わる]
一つ目の光子が $E_1=2mc^2$ である。正面衝突 $\theta=\pi$ では $E_2\ge mc^2/2$、直角では $E_2\ge mc^2$ である。同方向では全不変質量が零なので有限な $E_2$ で閾値に届かない。これは四元運動量の保存による条件であり、実際の率には偏光と場の相互作用の行列要素が必要である。
:::

### Klein 領域で Dirac 段差の流束を合わせる

一次元で同じスピンの二成分を選び、$\alpha=\sigma_x,\beta=\sigma_z$ として $H=c\sigma_xp+\sigma_z M+V$ とする。左で $E>M$ とし、入射と反射のスピノルは $(1,\eta)^\mathsf T$ と $(1,-\eta)^\mathsf T$、$\eta=c p/(E+M)>0$ である。右が $E'=E-V<-M$、群速度が右向きなら運動量は $-P$、$P>0$ である。下成分の比は $\eta'=(-cP)/(E'+M)>0$ なので、右流束は $2c\eta'|t|^2>0$ になる。一階の空間微分方程式でデルタ源が無ければ波動関数が連続であり、$1+r=t$、$\eta(1-r)=\eta't$ を得る。従って
$$
r=\frac{\eta-\eta'}{\eta+\eta'},\qquad
\mathcal R=|r|^2,\quad
\mathcal T=\frac{\eta'}{\eta}|t|^2=\frac{4\eta\eta'}{(\eta+\eta')^2}.
$$
この単粒子の Dirac 流束は正で $\mathcal R+\mathcal T=1$ である。負枝を右へ取る境界条件の解は、同じ電荷の正のエネルギー粒子が右の遠方まで伝播したという意味ではない。場の記述では左右における粒子・反粒子モードを再定義し、何を準備したかを指定する必要がある。これはスカラー Klein–Gordon の符号付き流束との相違でもある。

:::example[左右の運動量の大きさが等しい場合]
左 $E=5M/4$、右 $E'=-5M/4$ とする。$cp=cP=3M/4$ だから $\eta=1/3,\eta'=3$ である。$r=-4/5,t=1/5$ となり、$\mathcal R=16/25,\mathcal T=9/25$ で和は 1 である。四成分の問題でも運動と同じ軸のスピンを固定すれば等価な二成分に落とせる。スピンを逆に取ると下の二比の符号がともに反転し、流束を $\alpha$ の符号も含めて評価すれば同じ比を得る。
:::

### 正負の射影に共通する運動量の移行

自由状態を $\ket{u(p)}e^{ipz/\hbar}$ とし、位置に依存する摂動の一 Fourier 成分が $p$ を $p'$ に動かすとする。スピノル成分にはこの摂動が単位行列として働く。$z$ 軸の同じスピンを選び、$\cos\theta_p=M/E_p,\sin\theta_p=cp/E_p$ とすれば、正枝の二成分は $(\cos(\theta_p/2),\sin(\theta_p/2))$、負枝は $(-\sin(\theta_p/2),\cos(\theta_p/2))$ にできる。同じ列に $p'$ の射影を適用すると、負枝の振幅は $\sin((\theta_p-\theta_{p'})/2)$ である。同じ運動量なら零だが、$p$ と $p'$ が違うと新しい枝への射影が混ざる。これは空間的な行列要素の要因である。定常な摂動ではさらにエネルギー保存、時変な摂動ではその周波数の時間積分が掛かる。

:::example[運動量と枝の混合を別々に制御する]
$p=3mc/4$ から $p'=-3mc/4$ へ移す Fourier 成分を考える。$E_{p'}=E_p=5M/4$、$\theta_{p'}=-\theta_p$ なので負枝への投影の重みは $\sin^2\theta_p=9/25$ である。正枝から負枝へ本当に遷移するには、単粒子の自由エネルギーの差 $-2E_p$ に対応する時間依存が必要である。静的な弾性散乱で同じ $p'$ が出ても、負の枝が同じエネルギーの出射状態であるわけではない。場の表示では負枝の係数を反粒子として読み替えて、外部の系がどれだけエネルギーを供給したかを追う。
:::

## 演習問題

::::exercise{level=basic}
$d_{3/2}$ と $d_{5/2}$ の $\kappa$ を求めよ。
:::solution
$d$ は $\ell=2$ である。$j=\ell-1/2$ では $\kappa=\ell=2$、$j=\ell+1/2$ では $\kappa=-(\ell+1)=-3$ なので、順に $2,-3$ である。
:::
::::

::::exercise{level=basic}
$n=2,\kappa=2$ が許されるかを判断せよ。
:::solution
$n_r=n-|\kappa|=0$ で、$\kappa>0$ である。終端次数零では $-\kappa-\xi=0$ が必要だが $\xi>0$ なので満たせない。従ってこの束縛状態は存在しない。
:::
::::

::::exercise{level=basic}
Dirac のパリティの下で $j^0$ と $j^i$ の符号を求めよ。
:::solution
$j^0=c\psi^\dagger\psi$ では $\beta^\dagger\beta=I$ なので偶である。$j^i=c\psi^\dagger\alpha_i\psi$ では $\beta\alpha_i\beta=-\alpha_i$ なので奇である。比較する座標は $\vb r$ と $-\vb r$ である。
:::
::::

::::exercise{level=basic}
$\vb p=0$ の Zitterbewegung の周期を求めよ。
:::solution
時間因子の角周波数は $2mc^2/\hbar$ である。従って周期は $2\pi$ を角周波数で割った $\pi\hbar/(mc^2)$ である。
:::
::::

::::exercise{level=advanced}
水素の最低状態で全ノルムのうち小成分の割合を求めよ。
:::solution
$F=-sG$ なので割合は $s^2/(1+s^2)$ である。$s^2=(M-E)/(M+E)$ を入れると $(1-E/M)/2=(1-\sqrt{1-\alpha_Z^2})/2$ となる。弱結合では $\alpha_Z^2/4+\order{\alpha_Z^4}$ である。
:::
::::

::::exercise{level=advanced}
$E_{2,3/2}-E_{2,1/2}$ の先頭項を弱結合展開から求めよ。
:::solution
展開の共通部分を引くと、$j+1/2=2,1$ の逆数の差 $1/2-1=-1/2$ だけが残る。係数は $-mc^2\alpha_Z^4/(2\cdot2^3)$ なので積は $+mc^2\alpha_Z^4/32$ である。
:::
::::

::::exercise{level=advanced}
$D_i=\alpha_i-cp_iH^{-1}$ が異なる枝しか結ばないことを射影で示せ。
:::solution
$\{D_i,H\}=0$ で $E>0$ は運動量のスカラーなので $P_\pm D_i=D_iP_\mp$ である。従って $P_\pm D_iP_\pm=D_iP_\mp P_\pm=0$ である。どちらか一枝の状態だけでは :eqref[eq-aqm-zitter] の振動期待値は零である。
:::
::::

::::exercise{level=advanced}
Majorana 条件 $\psi^{\rm c}=\psi$ を満たす非零場に定数位相を掛けて同じ条件が保たれる場合を求めよ。
:::solution
反線形なので $(e^{i\theta}\psi)^{\rm c}=e^{-i\theta}\psi$ である。元と等しいには $e^{2i\theta}=1$、つまり $\theta\in\pi\Z$ が必要十分である。
:::
::::

::::exercise{level=advanced}
スピンに陽に依存しない弱い遮蔽 Coulomb ポテンシャルの第一 Born 近似で、遮蔽を外した非前方の強度を使う。$E=2mc^2$ で $\theta=\pi$ と $\pi/2$ の非偏極断面積の比を求めよ。さらに、同じ $E$ のまま運動量を同じ軸で反転させるとき、正枝の上向きスピンの重なりを直接求めて後方の係数を照合せよ。

:::solution
運動量の大きさは $p=\sqrt3\,mc$、群速度は $v=c\sqrt3/2$ である。スピン係数は $1-(3/4)\sin^2(\theta/2)$ であり、後方で $1/4$、直角で $5/8$ となる。遮蔽を外した Coulomb の Fourier 強度にはさらに $\sin^{-4}(\theta/2)$ が掛かる。二角度で共通の $g^2/(4p^2v^2)$ を除けば、後方は $1/4$、直角は $(5/8)\cdot 4=5/2$ なので、求める比は $1/10$ である。スピンのない Rutherford の角度因子だけなら $1/4$ なので、運動量の向きによる重なりの違いで後方がさらに弱くなる。

正枝の $z$ 上向きについて $p_z=p$ の上下の振幅は $\sqrt3/2,1/2$ である。$p_z=-p$ の正枝では下の成分の符号だけが反転する。従ってスピノルの重なりは $3/4-1/4=1/2$、絶対値二乗は $1/4$ である。反転先の逆向きスピンの上・下二成分はどちらも初期スピノルのスピンと直交するので、この角度ではスピンを反転した遷移振幅は零である。入口をもう一方のスピンにしても後方の係数は同じになる。これが非偏極の跡で得た値と一致する。
:::
::::

::::exercise{level=advanced}
自由な Weyl 基底の静止モードで初期の左と右を $(\chi,i\chi)/\sqrt2$、$\|\chi\|=1$ とする。左右の確率、時間零での確率の増減、正エネルギーの重みを求めよ。

:::solution
静止の質量 Hamiltonian は左右を入れ替える $M\sigma_x$ である。$u=Mt/\hbar$ と書けば、時間発展の二次元行列は $I\cos u-i\sigma_x\sin u$ なので、左は $(\cos u+\sin u)\chi/\sqrt2$、右は $i(\cos u-\sin u)\chi/\sqrt2$ である。それぞれの確率は $(1+\sin2u)/2$、$(1-\sin2u)/2$ になる。初期では左が $M/\hbar$ の率で増え、右が同じ率で減る。これは $2M\Im\inner{\psi_L}{\psi_R}/\hbar$ に $\inner{\psi_L}{\psi_R}=i/2$ を代入した値にも一致する。

正枝の射影は $(I+\sigma_x)/2$ であり、初期の $\sigma_x$ の平均は二つの係数の積の実部の二倍で零である。従って正枝の重みは $1/2$ で、時間発展で保たれる。$u=\pi/4$ で左が 1 となってもエネルギー枝がすべて正になるわけではない。左右の確率の変化は質量により二つの枝の相対位相が回ることから出る。
:::
::::

## 確認テスト

::::quiz
点核 Coulomb 模型の $2s_{1/2}$ と $2p_{1/2}$ のエネルギーを比較せよ。

- [x] 等しい
- [ ] $2s_{1/2}$ が必ず高い
- [ ] $2p_{1/2}$ が許されない

:::explanation
両者は $n=2$、$|\kappa|=1$、$n_r=1$ であり、符号が違ってもエネルギー式は同じである。$n_r=0$ だけにある正の $\kappa$ の制限にも触れない。
:::
::::

::::quiz{answer="-1" tol="1e-9"}
$p_{1/2}$ のスピノルのパリティ固有値を求めよ。
:::explanation
上成分の軌道量子数が $\ell=1$ なので $(-1)^\ell=-1$ である。下成分は $\ell'=0$ だが、空間反転時の $\beta$ の下成分の符号が負なので同じ固有値になる。
:::
::::

::::quiz
正のエネルギーの質量零粒子で右カイラリティに対応するものを選べ。

- [x] 正のヘリシティ
- [ ] 負のヘリシティ
- [ ] ヘリシティ零

:::explanation
右成分の Hamiltonian は $c\boldsymbol\sigma\cdot\vb p$ である。正エネルギーなら $\boldsymbol\sigma\cdot\hat{\vb p}$ の固有値が $+1$ で、ヘリシティは $+1/2$ になる。
:::
::::

::::quiz{answer="2" tol="1e-9"}
静止時の Zitterbewegung の角周波数を $mc^2/\hbar$ の単位で答えよ。
:::explanation
二つの枝は $+mc^2,-mc^2$ で、干渉の角周波数は差の大きさを $\hbar$ で割った $2mc^2/\hbar$ である。
:::
::::

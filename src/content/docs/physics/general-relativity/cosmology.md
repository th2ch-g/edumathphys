---
title: "相対論的宇宙論"
description: "一様等方な時空の力学を導き、宇宙論的赤方偏移、距離、因果的地平線を計算する。"
notation:
  - symbol: "a(t)"
    meaning: "宇宙の尺度因子"
  - symbol: "H"
    meaning: "Hubble パラメータ"
  - symbol: "k"
    meaning: "空間曲率の符号"
  - symbol: "\\eta"
    meaning: "共形時刻"
  - symbol: "D_L,D_A"
    meaning: "光度距離と角径距離"
---

## この章の目標

- 宇宙原理から Robertson–Walker 計量の形を得る。
- Friedmann 方程式を導き、放射・物質・宇宙項の解を求める。
- 赤方偏移と距離の関係を光の伝播から計算する。
- 加速膨張と地平線の関係を積分で判定する。

## 宇宙原理と Robertson–Walker 計量

:term[宇宙原理]{en="cosmological principle" yomi="うちゅうげんり"} は共動観測者が測る一定の宇宙時刻の空間に、一様性と等方性を課す仮定である。等方性から空間内の方向に区別がないため、各点の二平面の断面曲率が一致する。一様性によりその値は位置によらない。等方性は時間空間の混合成分を零にし、共動時刻を固有時に取れる。同じ回転と空間移動を異なる時刻に用いれば、空間計量の時間変化は全方向に共通の倍率 $a(t)^2$ になる。:term[Robertson–Walker 計量]{en="Robertson–Walker metric" yomi="ろばーとそんうぉーかーけいりょう"} は

$$
\dd s^2=-\dd t^2+a(t)^2[\dd\chi^2+S_k(\chi)^2\dd\Omega^2],\qquad
S_{+1}(\chi)=\sin\chi,\quad S_0(\chi)=\chi,\quad S_{-1}(\chi)=\sinh\chi.
\label{eq-gr-rw}
$$

球対称な空間の動径 Jacobi 場の長さを $S$ とすると $S''=-kS$、原点で $S(0)=0,S'(0)=1$ であり、上の三関数になる。これは :eqref[eq-gr-jacobi] の空間版である。空間の実際の曲率は $k/a^2$ である。例えば $k=1$ では $\chi=\pi/2$ の球の面積が $4\pi a^2$ で最大になり、$\chi=\pi$ で反対側の極に達する。局所計量だけでは許される大域的同一視まで決まらない。

### 固有距離と面積半径で空間曲率を読む

一定時刻の計量 $\dd\ell^2+a^2 S_k(\ell/a)^2\dd\Omega^2$ で、中心からの固有距離を $\ell=a\chi$ とする。小さな $\ell$ に対し $S_k(\chi)=\chi-k\chi^3/6+\order{\chi^5}$ なので、球の面積は
$$
4\pi\ell^2[1-k\ell^2/(3a^2)+\order{(\ell/a)^4}]
$$
となる。固有体積はこれを $\ell$ で積分するので $4\pi\ell^3[1-k\ell^2/(5a^2)+\order{(\ell/a)^4}]/3$ である。例えば $k=+1$ で $\ell=\pi a/2$ までの体積を厳密に積分すると $4\pi a^3\int_0^{\pi/2}\sin^2\chi\dd\chi=\pi^2a^3$ である。閉じた単位三球の全体積はこの二倍であり、球面の面積の極大を越えても空間は続く。

:::example[動径の形から Robertson–Walker の曲率を作る]
$a(t)>0,S(\chi)>0$ のパッチを取り、最初は動径関数を任意にする。等方な空間の条件は、空間曲率を計算してから $S''=-kS$ と $(S')^2=1-kS^2$ として入れる。 座標の順序を $(t,\chi,\theta,\phi)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & a^{2}{\left(t \right)} & 0 & 0\\0 & 0 & S^{2}{\left(\chi \right)} a^{2}{\left(t \right)} & 0\\0 & 0 & 0 & S^{2}{\left(\chi \right)} a^{2}{\left(t \right)} \sin^{2}{\left(\theta \right)}\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & \frac{1}{a^{2}{\left(t \right)}} & 0 & 0\\0 & 0 & \frac{1}{S^{2}{\left(\chi \right)} a^{2}{\left(t \right)}} & 0\\0 & 0 & 0 & \frac{1}{S^{2}{\left(\chi \right)} a^{2}{\left(t \right)} \sin^{2}{\left(\theta \right)}}\end{matrix}\right]\\
\det g&=- S^{4}{\left(\chi \right)} a^{6}{\left(t \right)} \sin^{2}{\left(\theta \right)}
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{t}_{\chi \chi}$ | $- \frac{1}{2}\left(-\partial_{t}\left(a^{2}{\left(t \right)}\right)\right)$ | $a{\left(t \right)} a'{\left(t \right)}$ |
| $\Gamma^{t}_{\theta \theta}$ | $- \frac{1}{2}\left(-\partial_{t}\left(S^{2}{\left(\chi \right)} a^{2}{\left(t \right)}\right)\right)$ | $S^{2}{\left(\chi \right)} a{\left(t \right)} a'{\left(t \right)}$ |
| $\Gamma^{t}_{\phi \phi}$ | $- \frac{1}{2}\left(-\partial_{t}\left(S^{2}{\left(\chi \right)} a^{2}{\left(t \right)} \sin^{2}{\left(\theta \right)}\right)\right)$ | $S^{2}{\left(\chi \right)} a{\left(t \right)} \sin^{2}{\left(\theta \right)} a'{\left(t \right)}$ |
| $\Gamma^{\chi}_{t \chi}$ | $\frac{1}{2 a^{2}{\left(t \right)}}\left(\partial_{t}\left(a^{2}{\left(t \right)}\right)\right)$ | $\frac{a'{\left(t \right)}}{a{\left(t \right)}}$ |
| $\Gamma^{\chi}_{\theta \theta}$ | $\frac{1}{2 a^{2}{\left(t \right)}}\left(-\partial_{\chi}\left(S^{2}{\left(\chi \right)} a^{2}{\left(t \right)}\right)\right)$ | $- S{\left(\chi \right)} S'{\left(\chi \right)}$ |
| $\Gamma^{\chi}_{\phi \phi}$ | $\frac{1}{2 a^{2}{\left(t \right)}}\left(-\partial_{\chi}\left(S^{2}{\left(\chi \right)} a^{2}{\left(t \right)} \sin^{2}{\left(\theta \right)}\right)\right)$ | $- S{\left(\chi \right)} \sin^{2}{\left(\theta \right)} S'{\left(\chi \right)}$ |
| $\Gamma^{\theta}_{t \theta}$ | $\frac{1}{2 S^{2}{\left(\chi \right)} a^{2}{\left(t \right)}}\left(\partial_{t}\left(S^{2}{\left(\chi \right)} a^{2}{\left(t \right)}\right)\right)$ | $\frac{a'{\left(t \right)}}{a{\left(t \right)}}$ |
| $\Gamma^{\theta}_{\chi \theta}$ | $\frac{1}{2 S^{2}{\left(\chi \right)} a^{2}{\left(t \right)}}\left(\partial_{\chi}\left(S^{2}{\left(\chi \right)} a^{2}{\left(t \right)}\right)\right)$ | $\frac{S'{\left(\chi \right)}}{S{\left(\chi \right)}}$ |
| $\Gamma^{\theta}_{\phi \phi}$ | $\frac{1}{2 S^{2}{\left(\chi \right)} a^{2}{\left(t \right)}}\left(-\partial_{\theta}\left(S^{2}{\left(\chi \right)} a^{2}{\left(t \right)} \sin^{2}{\left(\theta \right)}\right)\right)$ | $- \frac{\sin{\left(2 \theta \right)}}{2}$ |
| $\Gamma^{\phi}_{t \phi}$ | $\frac{1}{2 S^{2}{\left(\chi \right)} a^{2}{\left(t \right)} \sin^{2}{\left(\theta \right)}}\left(\partial_{t}\left(S^{2}{\left(\chi \right)} a^{2}{\left(t \right)} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{a'{\left(t \right)}}{a{\left(t \right)}}$ |
| $\Gamma^{\phi}_{\chi \phi}$ | $\frac{1}{2 S^{2}{\left(\chi \right)} a^{2}{\left(t \right)} \sin^{2}{\left(\theta \right)}}\left(\partial_{\chi}\left(S^{2}{\left(\chi \right)} a^{2}{\left(t \right)} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{S'{\left(\chi \right)}}{S{\left(\chi \right)}}$ |
| $\Gamma^{\phi}_{\theta \phi}$ | $\frac{1}{2 S^{2}{\left(\chi \right)} a^{2}{\left(t \right)} \sin^{2}{\left(\theta \right)}}\left(\partial_{\theta}\left(S^{2}{\left(\chi \right)} a^{2}{\left(t \right)} \sin^{2}{\left(\theta \right)}\right)\right)$ | $\frac{1}{\tan{\left(\theta \right)}}$ |

曲率の計算では :eqref[eq-gr-curvature-components] の微分差と積差を分ける。第二列は $\partial_c\Gamma^a{}_{db}$、第三列は $-\partial_d\Gamma^a{}_{cb}$ であり、第四列は接続の積を符号付きで足したものである。三列を足して最後の列を得る。

| 成分 | 第一の微分 | 第二の微分 | 接続の積 | 結果 |
|---|---|---|---|---|
| $R^{t}_{\chi t \chi}$ | $a{\left(t \right)} a''{\left(t \right)} + \left(a'{\left(t \right)}\right)^{2}$ | $0$ | $- \left(a'{\left(t \right)}\right)^{2}$ | $a{\left(t \right)} a''{\left(t \right)}$ |
| $R^{t}_{\theta t \theta}$ | $S^{2}{\left(\chi \right)} a{\left(t \right)} a''{\left(t \right)} + S^{2}{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2}$ | $0$ | $- S^{2}{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2}$ | $S^{2}{\left(\chi \right)} a{\left(t \right)} a''{\left(t \right)}$ |
| $R^{t}_{\phi t \phi}$ | $S^{2}{\left(\chi \right)} a{\left(t \right)} \sin^{2}{\left(\theta \right)} a''{\left(t \right)} + S^{2}{\left(\chi \right)} \sin^{2}{\left(\theta \right)} \left(a'{\left(t \right)}\right)^{2}$ | $0$ | $- S^{2}{\left(\chi \right)} \sin^{2}{\left(\theta \right)} \left(a'{\left(t \right)}\right)^{2}$ | $S^{2}{\left(\chi \right)} a{\left(t \right)} \sin^{2}{\left(\theta \right)} a''{\left(t \right)}$ |
| $R^{\chi}_{\theta \chi \theta}$ | $- S{\left(\chi \right)} S''{\left(\chi \right)} - \left(S'{\left(\chi \right)}\right)^{2}$ | $0$ | $S^{2}{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + \left(S'{\left(\chi \right)}\right)^{2}$ | $- \left(- S{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + S''{\left(\chi \right)}\right) S{\left(\chi \right)}$ |
| $R^{\chi}_{\phi \chi \phi}$ | $- S{\left(\chi \right)} \sin^{2}{\left(\theta \right)} S''{\left(\chi \right)} - \sin^{2}{\left(\theta \right)} \left(S'{\left(\chi \right)}\right)^{2}$ | $0$ | $S^{2}{\left(\chi \right)} \sin^{2}{\left(\theta \right)} \left(a'{\left(t \right)}\right)^{2} + \sin^{2}{\left(\theta \right)} \left(S'{\left(\chi \right)}\right)^{2}$ | $- \left(- S{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + S''{\left(\chi \right)}\right) S{\left(\chi \right)} \sin^{2}{\left(\theta \right)}$ |
| $R^{\theta}_{\phi \theta \phi}$ | $- \cos{\left(2 \theta \right)}$ | $0$ | $S^{2}{\left(\chi \right)} \sin^{2}{\left(\theta \right)} \left(a'{\left(t \right)}\right)^{2} - \sin^{2}{\left(\theta \right)} \left(S'{\left(\chi \right)}\right)^{2} + \frac{\sin{\left(2 \theta \right)}}{2 \tan{\left(\theta \right)}}$ | $- \left(- S^{2}{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + \left(S'{\left(\chi \right)}\right)^{2} - 1\right) \sin^{2}{\left(\theta \right)}$ |

第一添字を計量で下ろすと $R_{abcd}=g_{ak}R^k{}_{bcd}$ となる。対の反対称性と交換対称性を用いて、非零の独立な組を下表に示す。

| 成分 | 結果 |
|---|---|
| $R_{t \chi t \chi}$ | $- a{\left(t \right)} a''{\left(t \right)}$ |
| $R_{t \theta t \theta}$ | $- S^{2}{\left(\chi \right)} a{\left(t \right)} a''{\left(t \right)}$ |
| $R_{t \phi t \phi}$ | $- S^{2}{\left(\chi \right)} a{\left(t \right)} \sin^{2}{\left(\theta \right)} a''{\left(t \right)}$ |
| $R_{\chi \theta \chi \theta}$ | $- \left(- S{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + S''{\left(\chi \right)}\right) S{\left(\chi \right)} a^{2}{\left(t \right)}$ |
| $R_{\chi \phi \chi \phi}$ | $- \left(- S{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + S''{\left(\chi \right)}\right) S{\left(\chi \right)} a^{2}{\left(t \right)} \sin^{2}{\left(\theta \right)}$ |
| $R_{\theta \phi \theta \phi}$ | $- \left(- S^{2}{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + \left(S'{\left(\chi \right)}\right)^{2} - 1\right) S^{2}{\left(\chi \right)} a^{2}{\left(t \right)} \sin^{2}{\left(\theta \right)}$ |

Ricci テンソルは第一・第三添字を縮約して得る。曲率が残る場合も相殺が起こるので、各項を落とさずに加える。

| 成分 | 縮約の項 | 代入 | 結果 |
|---|---|---|---|
| $R_{t t}$ | $R^{\chi}_{t \chi t}+R^{\theta}_{t \theta t}+R^{\phi}_{t \phi t}$ | $- \frac{a''{\left(t \right)}}{a{\left(t \right)}} - \frac{a''{\left(t \right)}}{a{\left(t \right)}} - \frac{a''{\left(t \right)}}{a{\left(t \right)}}$ | $- \frac{3 a''{\left(t \right)}}{a{\left(t \right)}}$ |
| $R_{\chi \chi}$ | $R^{t}_{\chi t \chi}+R^{\theta}_{\chi \theta \chi}+R^{\phi}_{\chi \phi \chi}$ | $- \frac{- S{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + S''{\left(\chi \right)}}{S{\left(\chi \right)}} - \frac{- S{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + S''{\left(\chi \right)}}{S{\left(\chi \right)}} + a{\left(t \right)} a''{\left(t \right)}$ | $- \frac{- S{\left(\chi \right)} a{\left(t \right)} a''{\left(t \right)} - 2 S{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + 2 S''{\left(\chi \right)}}{S{\left(\chi \right)}}$ |
| $R_{\theta \theta}$ | $R^{t}_{\theta t \theta}+R^{\chi}_{\theta \chi \theta}+R^{\phi}_{\theta \phi \theta}$ | $- \left(- S{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + S''{\left(\chi \right)}\right) S{\left(\chi \right)} + \left(S^{2}{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} - \left(S'{\left(\chi \right)}\right)^{2} + 1\right) + S^{2}{\left(\chi \right)} a{\left(t \right)} a''{\left(t \right)}$ | $S^{2}{\left(\chi \right)} a{\left(t \right)} a''{\left(t \right)} + 2 S^{2}{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} - S{\left(\chi \right)} S''{\left(\chi \right)} - \left(S'{\left(\chi \right)}\right)^{2} + 1$ |
| $R_{\phi \phi}$ | $R^{t}_{\phi t \phi}+R^{\chi}_{\phi \chi \phi}+R^{\theta}_{\phi \theta \phi}$ | $- \left(- S{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + S''{\left(\chi \right)}\right) S{\left(\chi \right)} \sin^{2}{\left(\theta \right)} - \left(- S^{2}{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + \left(S'{\left(\chi \right)}\right)^{2} - 1\right) \sin^{2}{\left(\theta \right)} + S^{2}{\left(\chi \right)} a{\left(t \right)} \sin^{2}{\left(\theta \right)} a''{\left(t \right)}$ | $- \left(- S^{2}{\left(\chi \right)} a{\left(t \right)} a''{\left(t \right)} - 2 S^{2}{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + S{\left(\chi \right)} S''{\left(\chi \right)} + \left(S'{\left(\chi \right)}\right)^{2} - 1\right) \sin^{2}{\left(\theta \right)}$ |

逆計量と縮約すると Ricci スカラーとなる。さらに $G_{ab}=R_{ab}-g_{ab}R/2$ を成分ごとに計算する。

$$
\begin{aligned}
R&=\sum_{a,b}g^{ab}R_{ab}\\
&=- \frac{- S{\left(\chi \right)} a{\left(t \right)} a''{\left(t \right)} - 2 S{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + 2 S''{\left(\chi \right)}}{S{\left(\chi \right)} a^{2}{\left(t \right)}} - \frac{- S^{2}{\left(\chi \right)} a{\left(t \right)} a''{\left(t \right)} - 2 S^{2}{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + S{\left(\chi \right)} S''{\left(\chi \right)} + \left(S'{\left(\chi \right)}\right)^{2} - 1}{S^{2}{\left(\chi \right)} a^{2}{\left(t \right)}} + \frac{S^{2}{\left(\chi \right)} a{\left(t \right)} a''{\left(t \right)} + 2 S^{2}{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} - S{\left(\chi \right)} S''{\left(\chi \right)} - \left(S'{\left(\chi \right)}\right)^{2} + 1}{S^{2}{\left(\chi \right)} a^{2}{\left(t \right)}} + \frac{3 a''{\left(t \right)}}{a{\left(t \right)}}\\
&=- \frac{2 \left(- 3 S^{2}{\left(\chi \right)} a{\left(t \right)} a''{\left(t \right)} - 3 S^{2}{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + 2 S{\left(\chi \right)} S''{\left(\chi \right)} + \left(S'{\left(\chi \right)}\right)^{2} - 1\right)}{S^{2}{\left(\chi \right)} a^{2}{\left(t \right)}}
\end{aligned}
$$

非零の Einstein テンソル成分は次のとおりである。

| 成分 | Ricci 項 | 跡の項 | 差 |
|---|---|---|---|
| $G_{t t}$ | $- \frac{3 a''{\left(t \right)}}{a{\left(t \right)}}$ | $\frac{- 3 S^{2}{\left(\chi \right)} a{\left(t \right)} a''{\left(t \right)} - 3 S^{2}{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + 2 S{\left(\chi \right)} S''{\left(\chi \right)} + \left(S'{\left(\chi \right)}\right)^{2} - 1}{S^{2}{\left(\chi \right)} a^{2}{\left(t \right)}}$ | $- \frac{- 3 S^{2}{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + 2 S{\left(\chi \right)} S''{\left(\chi \right)} + \left(S'{\left(\chi \right)}\right)^{2} - 1}{S^{2}{\left(\chi \right)} a^{2}{\left(t \right)}}$ |
| $G_{\chi \chi}$ | $- \frac{- S{\left(\chi \right)} a{\left(t \right)} a''{\left(t \right)} - 2 S{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + 2 S''{\left(\chi \right)}}{S{\left(\chi \right)}}$ | $- \frac{- 3 S^{2}{\left(\chi \right)} a{\left(t \right)} a''{\left(t \right)} - 3 S^{2}{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + 2 S{\left(\chi \right)} S''{\left(\chi \right)} + \left(S'{\left(\chi \right)}\right)^{2} - 1}{S^{2}{\left(\chi \right)}}$ | $\frac{- 2 S^{2}{\left(\chi \right)} a{\left(t \right)} a''{\left(t \right)} - S^{2}{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + \left(S'{\left(\chi \right)}\right)^{2} - 1}{S^{2}{\left(\chi \right)}}$ |
| $G_{\theta \theta}$ | $S^{2}{\left(\chi \right)} a{\left(t \right)} a''{\left(t \right)} + 2 S^{2}{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} - S{\left(\chi \right)} S''{\left(\chi \right)} - \left(S'{\left(\chi \right)}\right)^{2} + 1$ | $3 S^{2}{\left(\chi \right)} a{\left(t \right)} a''{\left(t \right)} + 3 S^{2}{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} - 2 S{\left(\chi \right)} S''{\left(\chi \right)} - \left(S'{\left(\chi \right)}\right)^{2} + 1$ | $\left(- 2 S{\left(\chi \right)} a{\left(t \right)} a''{\left(t \right)} - S{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + S''{\left(\chi \right)}\right) S{\left(\chi \right)}$ |
| $G_{\phi \phi}$ | $- \left(- S^{2}{\left(\chi \right)} a{\left(t \right)} a''{\left(t \right)} - 2 S^{2}{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + S{\left(\chi \right)} S''{\left(\chi \right)} + \left(S'{\left(\chi \right)}\right)^{2} - 1\right) \sin^{2}{\left(\theta \right)}$ | $- \left(- 3 S^{2}{\left(\chi \right)} a{\left(t \right)} a''{\left(t \right)} - 3 S^{2}{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + 2 S{\left(\chi \right)} S''{\left(\chi \right)} + \left(S'{\left(\chi \right)}\right)^{2} - 1\right) \sin^{2}{\left(\theta \right)}$ | $\left(- 2 S{\left(\chi \right)} a{\left(t \right)} a''{\left(t \right)} - S{\left(\chi \right)} \left(a'{\left(t \right)}\right)^{2} + S''{\left(\chi \right)}\right) S{\left(\chi \right)} \sin^{2}{\left(\theta \right)}$ |

$R_{\chi\theta\chi\theta}$ を $g_{\chi\chi}g_{\theta\theta}$ で割ると、$a$ の時間変化に加え $-S''/(a^2S)$ が現れる。角二方向の同じ比には $(1-(S')^2)/(a^2S^2)$ が現れる。両者が同じ $k/a^2$ になるのは $S$ に上の二条件を課したときである。$\sin\chi,\chi,\sinh\chi$ は微分で条件を満たし、中心で $S(0)=0,S'(0)=1$ なので選択できる。
:::

### 一様等方性が空間の曲率を制限する

共動する観測者が時刻を共有する切片があるとする。観測者に直交する空間計量の時刻変化を局所直交枠で測ると、等方性のもとで対称二階テンソルの無跡部分は零になる。実際、どの単位ベクトルに沿って縮約しても同じ値なら、互いに直交する三方向が同じ固有値を持つ。その時間変化は一点で空間計量そのものに比例する。一様性により比例係数は時刻だけの関数なので、$\partial_t g_{ij}=2H(t)g_{ij}$ を積分して $g_{ij}=a^2(t)\gamma_{ij}$ と書ける。

同じ時刻に等方性を課せば一点でどの二次元断面も等しい曲率を持つ。また一様性が点の移動を許すため、曲率値は切片上で一定である。三次元の Riemann テンソルは断面曲率から偏極によって復元できる。具体的に固定した一点で二形式を三次元ベクトルと対応させると、曲率は二形式の空間上の対称二次形式である。任意の二形式が単純二形式である三次元では、全断面に同じ曲率 $k$ が成り立つと、その二次形式は $k$ 倍の恒等形式になる。従って $\gamma$ の曲率は $R_{ijkl}[\gamma]=k(\gamma_{ik}\gamma_{jl}-\gamma_{il}\gamma_{jk})$ となる。非零の $k$ の大きさを $\gamma$ の定数倍で吸収できるが、その際は $a$ も逆向きに変わる。実際の空間曲率は任意の規格化で $k/a^2$ である。

:::example[同じ Hubble 率で異なる空間曲率]
ある時刻の $H$ が同じ二つの模型を考える。$k=0$ の模型と $k=1,a=2/H$ の模型では、同じ動径固有距離 $d=a\chi$ の球の面積はそれぞれ $4\pi d^2$ と $4\pi a^2\sin^2(d/a)$ である。$d\ll a$ なら後者は前者より相対的に $d^2/(3a^2)$ 小さい。局所の膨張率を合わせても、同時刻の空間全体の曲率までは定まらない。
:::

## Friedmann 方程式

背景空間計量を $\gamma_{ij}$、$H=\dot a/a$ とする。時間を含む非零接続は $\Gamma^0{}_{ij}=a\dot a\,\gamma_{ij}$、$\Gamma^i{}_{0j}=H\delta^i_j$ で、純空間は $\gamma$ の接続である。曲率では $R^0{}_{i0j}=a\ddot a\gamma_{ij}$、$R^i{}_{jkl}=(k+\dot a^2)(\delta^i_k\gamma_{jl}-\delta^i_l\gamma_{jk})$ となる。前者では $(a\dot a)^{\cdot}-a\dot aH=a\ddot a$、後者では空間曲率 $k$ と時間添字の接続積 $\dot a^2$ が加わる。縮約すると $R_{00}=-3\ddot a/a$、$R_{ij}=[a\ddot a+2(\dot a^2+k)]\gamma_{ij}$ である。完全流体と :eqref[eq-gr-einstein-field] を使うと :term[Friedmann 方程式]{en="Friedmann equations" yomi="ふりーどまんほうていしき"} は

$$
H^2=\frac{8\pi G}3\rho-\frac{k}{a^2}+\frac{\Lambda}3,\qquad
\frac{\ddot a}a=-\frac{4\pi G}3(\rho+3p)+\frac{\Lambda}3 .
\label{eq-gr-friedmann}
$$

$\nabla_a u^a=(a^3)^{\cdot}/a^3=3H$ を :eqref[eq-gr-fluid-conservation] に入れて $\dot\rho+3H(\rho+p)=0$ を得る。$H\ne0$ では第一式を微分して保存則を用いれば第二式になる。一方、反転点を扱うには $H$ で割った式を根拠にせず、空間の Einstein 方程式を使う。例えば $k=0,\Lambda=0$ の臨界密度は $\rho=3H^2/(8\pi G)$ である。$H=70\,\mathrm{km\,s^{-1}\,Mpc^{-1}}$ を例に SI の質量密度へ直せば約 $9.2\times10^{-27}\,\mathrm{kg\,m^{-3}}$ である。

:::example[平坦切片の計量から時空の曲率を求める]
$a(t)>0$ を任意とし、$k=0$ の共動 Cartesian 座標で曲率を計算する。切片内の接続が零でも、計量の時間微分から時空の接続が生じる。 座標の順序を $(t,x,y,z)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & a^{2}{\left(t \right)} & 0 & 0\\0 & 0 & a^{2}{\left(t \right)} & 0\\0 & 0 & 0 & a^{2}{\left(t \right)}\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & \frac{1}{a^{2}{\left(t \right)}} & 0 & 0\\0 & 0 & \frac{1}{a^{2}{\left(t \right)}} & 0\\0 & 0 & 0 & \frac{1}{a^{2}{\left(t \right)}}\end{matrix}\right]\\
\det g&=- a^{6}{\left(t \right)}
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{t}_{x x}$ | $- \frac{1}{2}\left(-\partial_{t}\left(a^{2}{\left(t \right)}\right)\right)$ | $a{\left(t \right)} a'{\left(t \right)}$ |
| $\Gamma^{t}_{y y}$ | $- \frac{1}{2}\left(-\partial_{t}\left(a^{2}{\left(t \right)}\right)\right)$ | $a{\left(t \right)} a'{\left(t \right)}$ |
| $\Gamma^{t}_{z z}$ | $- \frac{1}{2}\left(-\partial_{t}\left(a^{2}{\left(t \right)}\right)\right)$ | $a{\left(t \right)} a'{\left(t \right)}$ |
| $\Gamma^{x}_{t x}$ | $\frac{1}{2 a^{2}{\left(t \right)}}\left(\partial_{t}\left(a^{2}{\left(t \right)}\right)\right)$ | $\frac{a'{\left(t \right)}}{a{\left(t \right)}}$ |
| $\Gamma^{y}_{t y}$ | $\frac{1}{2 a^{2}{\left(t \right)}}\left(\partial_{t}\left(a^{2}{\left(t \right)}\right)\right)$ | $\frac{a'{\left(t \right)}}{a{\left(t \right)}}$ |
| $\Gamma^{z}_{t z}$ | $\frac{1}{2 a^{2}{\left(t \right)}}\left(\partial_{t}\left(a^{2}{\left(t \right)}\right)\right)$ | $\frac{a'{\left(t \right)}}{a{\left(t \right)}}$ |

曲率の計算では :eqref[eq-gr-curvature-components] の微分差と積差を分ける。第二列は $\partial_c\Gamma^a{}_{db}$、第三列は $-\partial_d\Gamma^a{}_{cb}$ であり、第四列は接続の積を符号付きで足したものである。三列を足して最後の列を得る。

| 成分 | 第一の微分 | 第二の微分 | 接続の積 | 結果 |
|---|---|---|---|---|
| $R^{t}_{x t x}$ | $a{\left(t \right)} a''{\left(t \right)} + \left(a'{\left(t \right)}\right)^{2}$ | $0$ | $- \left(a'{\left(t \right)}\right)^{2}$ | $a{\left(t \right)} a''{\left(t \right)}$ |
| $R^{t}_{y t y}$ | $a{\left(t \right)} a''{\left(t \right)} + \left(a'{\left(t \right)}\right)^{2}$ | $0$ | $- \left(a'{\left(t \right)}\right)^{2}$ | $a{\left(t \right)} a''{\left(t \right)}$ |
| $R^{t}_{z t z}$ | $a{\left(t \right)} a''{\left(t \right)} + \left(a'{\left(t \right)}\right)^{2}$ | $0$ | $- \left(a'{\left(t \right)}\right)^{2}$ | $a{\left(t \right)} a''{\left(t \right)}$ |
| $R^{x}_{y x y}$ | $0$ | $0$ | $\left(a'{\left(t \right)}\right)^{2}$ | $\left(a'{\left(t \right)}\right)^{2}$ |
| $R^{x}_{z x z}$ | $0$ | $0$ | $\left(a'{\left(t \right)}\right)^{2}$ | $\left(a'{\left(t \right)}\right)^{2}$ |
| $R^{y}_{z y z}$ | $0$ | $0$ | $\left(a'{\left(t \right)}\right)^{2}$ | $\left(a'{\left(t \right)}\right)^{2}$ |

第一添字を計量で下ろすと $R_{abcd}=g_{ak}R^k{}_{bcd}$ となる。対の反対称性と交換対称性を用いて、非零の独立な組を下表に示す。

| 成分 | 結果 |
|---|---|
| $R_{t x t x}$ | $- a{\left(t \right)} a''{\left(t \right)}$ |
| $R_{t y t y}$ | $- a{\left(t \right)} a''{\left(t \right)}$ |
| $R_{t z t z}$ | $- a{\left(t \right)} a''{\left(t \right)}$ |
| $R_{x y x y}$ | $a^{2}{\left(t \right)} \left(a'{\left(t \right)}\right)^{2}$ |
| $R_{x z x z}$ | $a^{2}{\left(t \right)} \left(a'{\left(t \right)}\right)^{2}$ |
| $R_{y z y z}$ | $a^{2}{\left(t \right)} \left(a'{\left(t \right)}\right)^{2}$ |

Ricci テンソルは第一・第三添字を縮約して得る。曲率が残る場合も相殺が起こるので、各項を落とさずに加える。

| 成分 | 縮約の項 | 代入 | 結果 |
|---|---|---|---|
| $R_{t t}$ | $R^{x}_{t x t}+R^{y}_{t y t}+R^{z}_{t z t}$ | $- \frac{a''{\left(t \right)}}{a{\left(t \right)}} - \frac{a''{\left(t \right)}}{a{\left(t \right)}} - \frac{a''{\left(t \right)}}{a{\left(t \right)}}$ | $- \frac{3 a''{\left(t \right)}}{a{\left(t \right)}}$ |
| $R_{x x}$ | $R^{t}_{x t x}+R^{y}_{x y x}+R^{z}_{x z x}$ | $a{\left(t \right)} a''{\left(t \right)} + \left(a'{\left(t \right)}\right)^{2} + \left(a'{\left(t \right)}\right)^{2}$ | $a{\left(t \right)} a''{\left(t \right)} + 2 \left(a'{\left(t \right)}\right)^{2}$ |
| $R_{y y}$ | $R^{t}_{y t y}+R^{x}_{y x y}+R^{z}_{y z y}$ | $a{\left(t \right)} a''{\left(t \right)} + \left(a'{\left(t \right)}\right)^{2} + \left(a'{\left(t \right)}\right)^{2}$ | $a{\left(t \right)} a''{\left(t \right)} + 2 \left(a'{\left(t \right)}\right)^{2}$ |
| $R_{z z}$ | $R^{t}_{z t z}+R^{x}_{z x z}+R^{y}_{z y z}$ | $a{\left(t \right)} a''{\left(t \right)} + \left(a'{\left(t \right)}\right)^{2} + \left(a'{\left(t \right)}\right)^{2}$ | $a{\left(t \right)} a''{\left(t \right)} + 2 \left(a'{\left(t \right)}\right)^{2}$ |

逆計量と縮約すると Ricci スカラーとなる。さらに $G_{ab}=R_{ab}-g_{ab}R/2$ を成分ごとに計算する。

$$
\begin{aligned}
R&=\sum_{a,b}g^{ab}R_{ab}\\
&=\frac{a{\left(t \right)} a''{\left(t \right)} + 2 \left(a'{\left(t \right)}\right)^{2}}{a^{2}{\left(t \right)}} + \frac{a{\left(t \right)} a''{\left(t \right)} + 2 \left(a'{\left(t \right)}\right)^{2}}{a^{2}{\left(t \right)}} + \frac{a{\left(t \right)} a''{\left(t \right)} + 2 \left(a'{\left(t \right)}\right)^{2}}{a^{2}{\left(t \right)}} + \frac{3 a''{\left(t \right)}}{a{\left(t \right)}}\\
&=\frac{6 \left(a{\left(t \right)} a''{\left(t \right)} + \left(a'{\left(t \right)}\right)^{2}\right)}{a^{2}{\left(t \right)}}
\end{aligned}
$$

非零の Einstein テンソル成分は次のとおりである。

| 成分 | Ricci 項 | 跡の項 | 差 |
|---|---|---|---|
| $G_{t t}$ | $- \frac{3 a''{\left(t \right)}}{a{\left(t \right)}}$ | $- \frac{3 \left(a{\left(t \right)} a''{\left(t \right)} + \left(a'{\left(t \right)}\right)^{2}\right)}{a^{2}{\left(t \right)}}$ | $\frac{3 \left(a'{\left(t \right)}\right)^{2}}{a^{2}{\left(t \right)}}$ |
| $G_{x x}$ | $a{\left(t \right)} a''{\left(t \right)} + 2 \left(a'{\left(t \right)}\right)^{2}$ | $3 a{\left(t \right)} a''{\left(t \right)} + 3 \left(a'{\left(t \right)}\right)^{2}$ | $- 2 a{\left(t \right)} a''{\left(t \right)} - \left(a'{\left(t \right)}\right)^{2}$ |
| $G_{y y}$ | $a{\left(t \right)} a''{\left(t \right)} + 2 \left(a'{\left(t \right)}\right)^{2}$ | $3 a{\left(t \right)} a''{\left(t \right)} + 3 \left(a'{\left(t \right)}\right)^{2}$ | $- 2 a{\left(t \right)} a''{\left(t \right)} - \left(a'{\left(t \right)}\right)^{2}$ |
| $G_{z z}$ | $a{\left(t \right)} a''{\left(t \right)} + 2 \left(a'{\left(t \right)}\right)^{2}$ | $3 a{\left(t \right)} a''{\left(t \right)} + 3 \left(a'{\left(t \right)}\right)^{2}$ | $- 2 a{\left(t \right)} a''{\left(t \right)} - \left(a'{\left(t \right)}\right)^{2}$ |

$G_{tt}=3(\dot a/a)^2$ である。空間成分を $a^2$ で割ると $G^i{}_i=-(2\ddot a/a+(\dot a/a)^2)$ なので、圧力方程式には $\ddot a$ が残る。$a$ が一定なら時空の曲率まで零であるが、切片の平坦性は $a$ の時間変化を制約しない。
:::

## 宇宙論的赤方偏移と Hubble の法則

共動光源と観測者の間の動径零線では $\dd\chi=\pm\dd t/a$ である。同じ共動距離を後続の波頂にも適用すると、$\delta t_{\rm o}/a_{\rm o}=\delta t_{\rm e}/a_{\rm e}$ となる。従って :term[宇宙論的赤方偏移]{en="cosmological redshift" yomi="うちゅうろんてきせきほうへんい"} は

$$
1+z=\frac{\nu_{\rm e}}{\nu_{\rm o}}=\frac{a_{\rm o}}{a_{\rm e}},\qquad
\chi=\frac1{a_{\rm o}}\int_0^z\frac{\dd z'}{H(z')}
\label{eq-gr-cosmo-redshift}
$$

である。距離式では $\dd z/\dd t=-(1+z)H$ を使った。$z\ll1$ なら現在の動径固有距離 $D=a_{\rm o}\chi$ に対して $z\simeq H_{\rm o}D$ であり、SI では $cz\simeq H_{\rm o}D$ となる。例えば小さい赤方偏移で $z=0.01$、$H_{\rm o}=70\,\mathrm{km\,s^{-1}\,Mpc^{-1}}$ なら約 $43\,\mathrm{Mpc}$ である。遠方では $cz$ を共動天体の相対速度とみなすと、距離や同時面の選択を見失う。

### 赤方偏移の時刻ドリフト

同じ共動源の異なる放射時刻を比べる。観測者の時刻 $t_o$ と放射者の時刻 $t_e$ に対して $\chi=\int_{t_e}^{t_o}\dd t/a(t)$ は一定なので、微分して $\dd t_e/\dd t_o=a(t_e)/a(t_o)=1/(1+z)$ となる。$1+z=a(t_o)/a(t_e)$ の微分から
$$
\frac{\dd z}{\dd t_o}=(1+z)H_o-H_e
$$
を得る。$H_e$ を $H_o$ で置き換えると膨張の時間変化を捨てることになる。例えば平坦な物質のみの模型では $H_e=H_o(1+z)^{3/2}$ なので、$z>0$ でドリフトは負である。一方、指数膨張で両端が同じ定数 $H$ を持てば $Hz>0$ である。いずれも共動源の比較であり、固有運動をする源では運動の時間変化が加わる。

### 低赤方偏移の観測から減速率を求める

現在を $t_0$、尺度因子の現在値を $a_0=a(t_0)$ とする。曲率を $k=\pm1$ に規格化した場合も使えるよう、比 $a(t)/a_0$ を展開する。Taylor 展開を
$$
\frac{a(t_0-\Delta t)}{a_0}=1-H_0\Delta t-\tfrac12 q_0 H_0^2\Delta t^2+\order{\Delta t^3},
\qquad q_0=-\frac{\ddot a(t_0)}{a_0 H_0^2}
$$
と書く。逆数から $z=H_0\Delta t+(1+q_0/2)H_0^2\Delta t^2+\cdots$ である。これを解くと $H_0\Delta t=z-(1+q_0/2)z^2+\cdots$ となる。一方で現在の動径固有距離 $D_C$ は $a_0\int_{t_0-\Delta t}^{t_0}\dd t/a(t)=\Delta t+H_0\Delta t^2/2+\cdots$ なので $H_0D_C=z-(1+q_0)z^2/2+\cdots$ である。曲率は横断距離を三次から変えるため、光度距離は
$$
H_0D_L=z+\tfrac12(1-q_0)z^2+\order{z^3}
$$
となる。従って一次の Hubble の法則は $H_0$ を、二次は加減速の $q_0$ を含む。SI では左辺を $H_0D_L/c$ とする。

:::example[低赤方偏移で二つの模型を比較する]
平坦な塵では $q_0=1/2$、指数膨張では $q_0=-1$ である。例えば $z=0.05$ で $H_0D_L$ の二次項は塵で $0.000625$、指数膨張で $0.0025$ となる。同じ $H_0$ を固定してもこの次数の距離は異なる。厳密な指数膨張の平坦模型では $D_C=z/H_0$ なので $D_L=z(1+z)/H_0$ がそのまま成り立ち、ここでの二次形と一致する。
:::

## 宇宙論における距離（光度距離・角径距離）

物体の横方向の固有サイズ $s$ と角度 $\delta\theta$ から :term[角径距離]{en="angular diameter distance" yomi="かくけいきょり"} を $D_A=s/\delta\theta=a_{\rm e}S_k(\chi)$ と定める。ここでは共役点や多重像に達する前の光線束を使い、$k=1$ では $0<\chi<\pi$ とする。固有光度 $L$ を等方放射する光源から受ける流束は、光のエネルギーと到達率が各 $1/(1+z)$ 倍になり、光線を横切る現在の面積が $4\pi a_{\rm o}^2 S_k(\chi)^2$ になるため

$$
F=\frac{L}{4\pi a_{\rm o}^2S_k(\chi)^2(1+z)^2},\qquad
D_L=(1+z)a_{\rm o} S_k(\chi)=(1+z)^2D_A .
\label{eq-gr-distance-duality}
$$

ここで :term[光度距離]{en="luminosity distance" yomi="こうどきょり"} $D_L$ を $F=L/(4\pi D_L^2)$ で定義し、光子の吸収を無視した。例えば $z=2$ では同じ光源の $D_L$ は $D_A$ の九倍である。平坦な物質優勢モデルでは $H(z)=H_0(1+z)^{3/2}$ なので $a_0\chi=2H_0^{-1}[1-(1+z)^{-1/2}]$ である。$D_A$ を微分すると最大は $1+z=9/4$、すなわち $z=5/4$ であり、赤方偏移に対して角の大きさは単調でない。

### 光度距離と角径距離をモデル内で比較する

平坦な物質のみの模型で $H(z)=H_0(1+z)^{3/2}$ とすると、共動距離は $D_C=2[1-(1+z)^{-1/2}]/H_0$ である。光度距離は $(1+z)D_C$、角径距離は $D_C/(1+z)$ となる。$s=1+z$ とすれば $D_A=2(s^{-1}-s^{-3/2})/H_0$ であり、その微分の零点は $-s^{-2}+(3/2)s^{-5/2}=0$、従って $s=9/4$ である。そこで $D_A=8/(27H_0)$ を取る。$z$ がこれより大きくなると、同じ横方向の固有サイズを持つ物体の角度は再び大きくなる。この模型では $D_L$ は同じ範囲で単調増加であるため、距離の種類を指定することが必要である。

:::example[同じ赤方偏移に割り当てる四つの距離と時間]
平坦な物質だけの宇宙で、$H(z)=H_0(1+z)^{3/2}$ とする。共動距離は $\int_0^z H^{-1}\dd z'$、光度距離と角径距離はそれぞれこれの $(1+z)$ 倍とその逆数倍である。観測時刻から放射時刻までの経過は $\int_0^z[(1+z')H(z')]^{-1}\dd z'$ となる。どれも $c=1$ で $H_0$ を掛けて無次元にそろえる。

ここでは$z$ は共動光源の赤方偏移であるとする。計算に用いる式を先に並べる。

$$
\begin{aligned}
H_0D_C&=2 - \frac{2}{\sqrt{z + 1}}\\
H_0D_L&=\left(1 - \frac{1}{\sqrt{z + 1}}\right) \left(2 z + 2\right)\\
H_0D_A&=\frac{2 - \frac{2}{\sqrt{z + 1}}}{z + 1}\\
H_0(t_0-t_e)&=\frac{2}{3} - \frac{2}{3 \left(z + 1\right)^{\frac{3}{2}}}
\end{aligned}
$$

次の条件を比較する。入力を区別してから計算し、定数の規格化をそろえた値を表に示す。

| $z$ | $H_0D_C$ | $H_0D_L$ | $H_0D_A$ | $H_0(t_0-t_e)$ |
|---|---|---|---|---|
| $0$ | $0$ | $0$ | $0$ | $0$ |
| $\frac{5}{4}$ | $\frac{2}{3}$ | $\frac{3}{2}$ | $\frac{8}{27}$ | $\frac{38}{81}$ |
| $3$ | $1$ | $4$ | $\frac{1}{4}$ | $\frac{7}{12}$ |
| $8$ | $\frac{4}{3}$ | $12$ | $\frac{4}{27}$ | $\frac{52}{81}$ |
| $15$ | $\frac{3}{2}$ | $24$ | $\frac{3}{32}$ | $\frac{21}{32}$ |

$z=0$ では四量とも零である。最大角径距離は $z=5/4$ に現れ、その後は小さくなる。他方、光度距離は増え、時間差はモデルの年齢 $2/(3H_0)$ に近づく。宇宙の距離を語るときは、光の流束・角度・経過時間のどれから定義したかを指定する。
:::

### 距離双対性に必要な光の保存

共動する放射源が固有時 $\Delta t_e$ の間に光子 $N_\gamma$ 個を等方に放つとする。一光子のエネルギーは到着時に $1/(1+z)$ 倍、到着間隔は $(1+z)\Delta t_e$ である。伝播中に光子が失われなければ、観測者を通る球で受ける単位面積の流束は
$$
F=\frac{L_e}{4\pi(a_o S_k(\chi))^2(1+z)^2}
$$
である。光度距離の定義 $F=L_e/(4\pi D_L^2)$ と、放射面で測る $D_A=a_e S_k(\chi)$ を比べれば :eqref[eq-gr-distance-duality] になる。途中の透明度を $P$ とすると $F$ に $P$ が掛かり、流束から推定した光度距離は $D_L/\sqrt P$ である。例えば光子の半分だけが到達する $P=1/2$ では、固有の光度を既知とした推定値が $\sqrt2$ 倍になる。角径距離の幾何の式にその因子は入らない。

### 光線束の曲率と角径距離

同じ観測点から出る近接光線の分離は、アフィン径数 $\lambda$ に沿う横方向の測地線偏差で決まる。Robertson–Walker の対称性から、動径光に沿う二つの横方向には同じ潮汐が働き、円形断面はせん断を持たない。横枠 $e_A$ で偏差式を縮約すると
$$
\frac{\dd^2 D_A}{\dd\lambda^2}=-\frac12R_{ab}k^ak^b\,D_A
=-4\pi G(\rho+p)(u\cdot k)^2D_A
$$
となる。ここで $k^a=\dd x^a/\dd\lambda$、$u$ は共動四元速度である。第一式では零方向と補助の零方向を取り、逆計量を横の二方向と二つの零方向の組に分けた。Riemann の反対称性により $R_{ab}k^ak^b$ の縮約に零方向の組は寄与しないため、横方向の潮汐の和がちょうど $R_{ab}k^ak^b$ になる。第二式では Einstein 方程式を零ベクトルで縮約した。従って $\Lambda$ はこの局所収束の式に直接は入らないが、$a(t)$ やアフィン径数と赤方偏移の関係を変える。

観測点で $|u\cdot k|=1$ となるように尺度を選び、過去向きに $\lambda$ を増やす。面積半径は観測点で零で、$D_A'(0)=1$ である。$k^t=-a_0/a$、$k^\chi=a_0/a^2$ を $D_A=aS_k(\chi)$ の微分へ入れると
$$
\frac{\dd^2 (aS_k)}{\dd\lambda^2}
=\left(\frac{a_0}{a}\right)^2
\left(\dot H-\frac{k}{a^2}\right)aS_k.
$$
途中で $S_k''=-kS_k$ を使い、$S_k'$ を持つ交差項は打ち消し合う。Friedmann 方程式から $\dot H-k/a^2=-4\pi G(\rho+p)$ なので、上の光線束の式が再現する。例えば平坦な de Sitter の膨張枝では $\rho+p=0$ であり、$D_A=\lambda$ となる。実際、$H$ 一定なら $\lambda=\int_{t_e}^{t_o}a(t)\dd t/a_0=z/[H(1+z)]$ であり、光度距離から求めた $D_A$ と一致する。

## 物質・放射・宇宙項が支配する宇宙モデル

一定の状態方程式 $p=w\rho$ なら保存則から $\rho\propto a^{-3(1+w)}$ となる。$k=\Lambda=0$ の単一成分の膨張解では $\dot a\propto a^{-(1+3w)/2}$ なので $w>-1$ なら $a\propto t^{2/[3(1+w)]}$ である。放射 $w=1/3$ は $a\propto t^{1/2}$、物質 $w=0$ は $a\propto t^{2/3}$、正の宇宙項だけなら $a\propto e^{\sqrt{\Lambda/3}t}$ である。膨張とともに放射の密度は $a^{-4}$、塵の密度は $a^{-3}$ で減る。従って放射の塵に対する割合は $a^{-1}$ で減る。例えば $a$ を二倍にすると $\rho_r/\rho_m$ は半分になる。平坦宇宙の物質と正の宇宙項では $\ddot a>0$ の条件は $\rho_m<2\rho_\Lambda$ である。曲率項が支配する $k=-1,\rho=\Lambda=0$ では $a=t$ であり、座標変換 $T=t\cosh\chi,R=t\sinh\chi$ によって Minkowski 時空の光円錐内の表示になる。

### 密度パラメータから観測量へ

現在値に添字 $0$ を付け、$\Omega_{m0}=8\pi G\rho_{m0}/(3H_0^2)$、$\Omega_{r0}=8\pi G\rho_{r0}/(3H_0^2)$、$\Omega_{\Lambda0}=\Lambda/(3H_0^2)$、$\Omega_{k0}=-k/(a_0^2H_0^2)$ と定める。:eqref[eq-gr-friedmann] に各成分の保存式を代入すると

$$
\frac{H(z)^2}{H_0^2}=\Omega_{r0}(1+z)^4+\Omega_{m0}(1+z)^3+\Omega_{k0}(1+z)^2+\Omega_{\Lambda0}.
$$

右辺の係数の和は $1$ である。曲率の項は $(1+z)^2$ で変わるが、流体の圧力として物質の式に代入する量ではない。膨張枝ではこの正根を :eqref[eq-gr-cosmo-redshift] に使う。また同じ $\dd z/\dd t$ の式から、放射と観測の宇宙時刻の差は $\int_0^z\dd z'/[(1+z')H(z')]$ になる。

:::example[過去を見る時間と距離の違い]
物質だけの平坦な宇宙で $z=3$ の光を見る。源が光を出した宇宙時刻は $t_{\rm e}=t_0(1+z)^{-3/2}=t_0/8$ なので、経過時間は $7t_0/8$ である。一方、現在の動径距離は $2H_0^{-1}(1-1/2)=H_0^{-1}=3t_0/2$ である。経過した固有時間に光速を掛けた距離とは異なる。光の通過中に共動点間の固有距離が変わったためである。
:::

### 塵と宇宙項が共存する平坦模型

現在 $a_0=1$ と規格化し、$H_\Lambda=\sqrt{\Lambda/3}>0$、$C=8\pi G\rho_{m0}/3>0$ とする。第一 Friedmann 式は $\dot a^2/a^2=C/a^3+H_\Lambda^2$ である。$y=a^{3/2}$ と置くと膨張枝で $\dot y=(3/2)\sqrt{C+H_\Lambda^2 y^2}$ になる。$a(0)=0$ として分離積分すれば
$$
a(t)=\left(\frac{\sqrt C}{H_\Lambda}\sinh\frac{3H_\Lambda t}{2}\right)^{2/3}.
$$
初期には $a\simeq(3\sqrt C t/2)^{2/3}$、後期には定数倍した $e^{H_\Lambda t}$ である。加速度は $\ddot a/a=H_\Lambda^2-C/(2a^3)$ なので、加速の開始は $a^3=C/(2H_\Lambda^2)$ である。密度が等しいのは $a^3=C/H_\Lambda^2$ であり、宇宙項は圧力ももつため二つの時刻は違う。例えば $C/H_\Lambda^2=1/2$ なら、加速開始の赤方偏移は $4^{1/3}-1$、密度が等しい赤方偏移は $2^{1/3}-1$ となる。

### 閉じた塵宇宙の折返しを解く

宇宙項がなく曲率 $k=1$ の塵宇宙を考える。$\rho a^3$ が一定なので $\dot a^2=C/a-1$、$C=(8\pi G/3)\rho a^3$ となる。右辺が零になる最大尺度は $a=C$ である。そこで $a=C(1-\cos\eta)/2$、$\dd t=a\dd\eta$ と選ぶと
$$
t=\frac C2(\eta-\sin\eta),\qquad \dot a=\frac{\sin\eta}{1-\cos\eta}.
$$
後者の二乗は $(1+\cos\eta)/(1-\cos\eta)=C/a-1$ である。$0<\eta<\pi$ で膨張、$\pi<\eta<2\pi$ で収縮であり、折返しで $H=0$ になる。微分した式を $\dot a$ で割ると折返しで使えないが、空間の Einstein 方程式は $\ddot a=-C/(2a^2)$ を直接与える。そのため最大尺度で $\ddot a=-1/(2C)<0$ となり、時間を越えて収縮側へ滑らかに続く。

:::example[最大膨張までの時間と面積半径]
$\eta=\pi$ で $a=C,t=\pi C/2$ である。共動角 $\chi=\pi/3$ の球の面積半径は $\sqrt3C/2$、その球までの固有距離は $\pi C/3$ である。後者の方が大きいのは正曲率の球の面積半径が固有距離より小さくなるためである。寿命は $t(2\pi)=\pi C$ で、最大膨張の時刻はそのちょうど半分である。
:::

### 放射と塵が共存する解を共形時刻で積分する

平坦な模型で宇宙項を零とし、尺度因子を無次元に取る。密度の保存から $H^2=A/a^4+B/a^3$、$A,B>0$ と書く。宇宙時刻では $\dot a$ の右辺に $a^{-1}$ が残るので、$\dd t=a\dd\eta$ として
$$
(a')^2=A+Ba
$$
を得る。膨張側で $a(0)=0$ とすれば $\sqrt{A+Ba}=\sqrt A+B\eta/2$、従って
$$
a(\eta)=\sqrt A\,\eta+\frac B4\eta^2,\qquad
t(\eta)=\frac{\sqrt A}{2}\eta^2+\frac{B}{12}\eta^3 .
$$
$\eta$ が小さいと $a\propto\eta$、$t\propto\eta^2$ で $a\propto t^{1/2}$ となる。後半では $a\propto\eta^2$、$t\propto\eta^3$ から $a\propto t^{2/3}$ である。従って二つの極限を別の初期条件で貼る必要がなく、同じ解で滑らかにつながる。支配が交代する条件は密度が等しい $a=A/B$ であり、尺度因子の二つの項が等しい条件とは異なる。

:::example[二つの成分が等しい時刻]
両密度が等しい点では $A+Ba=2A$ なので $\eta_{\rm eq}=2(\sqrt2-1)\sqrt A/B$ である。共動粒子地平線は $\eta_{\rm eq}$ であり、固有距離は $a_{\rm eq}\eta_{\rm eq}=2(\sqrt2-1)A^{3/2}/B^2$ となる。一方、その時刻の Hubble 半径は $a_{\rm eq}^2/\sqrt{2A}=A^{3/2}/(\sqrt2 B^2)$ である。比は $4-2\sqrt2$ で、全過去の膨張とその時刻の膨張率が関わる量は一致しない。$A=B=t_s^{-2}$ と選べば $a_{\rm eq}=1$、$\eta_{\rm eq}=2(\sqrt2-1)t_s$ なので両方を時間の単位 $t_s$ で比べられる。
:::

## de Sitter 時空と反 de Sitter 時空

真空の定曲率 $K=\pm1/L^2$ の例を埋め込みで構成する。五次元で一つの時間方向がある平坦時空の二次曲面 $X\cdot X=L^2$ は :term[de Sitter 時空]{en="de Sitter spacetime" yomi="どじったーじくう"} である。二つの時間方向をもつ平坦時空の $X\cdot X=-L^2$ は :term[反 de Sitter 時空]{en="anti-de Sitter spacetime" yomi="はんどじったーじくう"} の局所模型である。法線 $n=X/L$ はノルム $\pm1$、平坦接続で $\nabla_v n=v/L$ である。接続を接成分と法線成分に分け、二度微分の接成分を反対称化すると曲率は $\pm(g_{ac}g_{bd}-g_{ad}g_{bc})/L^2$ になる。平坦空間の線形等長変換のうち二次形式を保つものが十次元作用を与え、両方とも局所的に最大対称である。

de Sitter の大域座標 $X^0=L\sinh(t/L)$、$\vb X=L\cosh(t/L)\hat n$ は $k=1$、$a=L\cosh(t/L)$ を与える。AdS では二つの時間座標を $(\sqrt{L^2+r^2}\cos(t/L),\sqrt{L^2+r^2}\sin(t/L))$、空間座標を $r\hat n$ とすると $\dd s^2=-(1+r^2/L^2)\dd t^2+(1+r^2/L^2)^{-1}\dd r^2+r^2\dd\Omega^2$ である。時間を実数へ展開して閉じた時間的曲線を除く。$r=0$ から無限遠に向かう光の座標時間は $\int_0^\infty\dd r/(1+r^2/L^2)=\pi L/2$ であり、共形境界は時間的になる。

### 反 de Sitter の共形境界へ届く光

四次元の静的形を
$$
\dd s^2=-(1+r^2/L^2)\dd t^2+(1+r^2/L^2)^{-1}\dd r^2+r^2\dd\Omega^2
$$
とする。$r=L\tan\chi,t=LT$ と置くと、$0\le\chi<\pi/2$ で $\dd s^2=L^2\sec^2\chi(-\dd T^2+\dd\chi^2+\sin^2\chi\dd\Omega^2)$ になる。角度で閉じる埋込み時刻ではなく $T\in\R$ の普遍被覆を取る。物理的な動径光の座標速度は $\dd\chi/\dd T=\pm1$ なので、中心から境界までの座標時間は $\pi L/2$ である。境界の誘導計量は $-\dd T^2+\dd\Omega^2$ で時間的である。このため境界が完全に未来にしか関わらない面とは異なり、場の波動方程式を立てる際は境界から入り得る自由度を指定する必要がある。物理計量のアフィン径数では共形因子の二乗が効くため、有限の $T$ で境界に達することは有限のアフィン径数を意味しない。

:::example[Milne と de Sitter で膨張の意味を比べる]
座標変換 $T=t\cosh\chi$、$r=t\sinh\chi$ により Minkowski の未来光円錐内の計量を $t>0,\chi>0$ で書く。空間形状は $k=-1$、尺度因子は $a(t)=t$ である。 座標の順序を $(t,\chi,\theta,\phi)$ と固定する。成分行列と逆行列は次の組になる。

$$
\begin{aligned}
(g_{ab})&=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & t^{2} & 0 & 0\\0 & 0 & t^{2} \sinh^{2}{\left(\chi \right)} & 0\\0 & 0 & 0 & t^{2} \sin^{2}{\left(\theta \right)} \sinh^{2}{\left(\chi \right)}\end{matrix}\right]\\
(g^{ab})&=\left[\begin{matrix}-1 & 0 & 0 & 0\\0 & \frac{1}{t^{2}} & 0 & 0\\0 & 0 & \frac{1}{t^{2} \sinh^{2}{\left(\chi \right)}} & 0\\0 & 0 & 0 & \frac{1}{t^{2} \sin^{2}{\left(\theta \right)} \sinh^{2}{\left(\chi \right)}}\end{matrix}\right]\\
\det g&=- t^{6} \sin^{2}{\left(\theta \right)} \sinh^{4}{\left(\chi \right)}
\end{aligned}
$$

逆行列を掛けて単位行列を得ることを確かめてから、:eqref[eq-gr-christoffel] を使う。同じ下添字の順序を二度計算せず、$b\le c$ の組を選ぶ。中央の列には、座標に依存する計量成分の微分だけを残した。最後の列は微分を実行した値であり、表に出ない成分は添字交換で得られるものを除いて零である。

| 成分 | 微分への代入 | 結果 |
|---|---|---|
| $\Gamma^{t}_{\chi \chi}$ | $- \frac{1}{2}\left(-\partial_{t}\left(t^{2}\right)\right)$ | $t$ |
| $\Gamma^{t}_{\theta \theta}$ | $- \frac{1}{2}\left(-\partial_{t}\left(t^{2} \sinh^{2}{\left(\chi \right)}\right)\right)$ | $t \sinh^{2}{\left(\chi \right)}$ |
| $\Gamma^{t}_{\phi \phi}$ | $- \frac{1}{2}\left(-\partial_{t}\left(t^{2} \sin^{2}{\left(\theta \right)} \sinh^{2}{\left(\chi \right)}\right)\right)$ | $t \sin^{2}{\left(\theta \right)} \sinh^{2}{\left(\chi \right)}$ |
| $\Gamma^{\chi}_{t \chi}$ | $\frac{1}{2 t^{2}}\left(\partial_{t}\left(t^{2}\right)\right)$ | $\frac{1}{t}$ |
| $\Gamma^{\chi}_{\theta \theta}$ | $\frac{1}{2 t^{2}}\left(-\partial_{\chi}\left(t^{2} \sinh^{2}{\left(\chi \right)}\right)\right)$ | $- \sinh{\left(\chi \right)} \cosh{\left(\chi \right)}$ |
| $\Gamma^{\chi}_{\phi \phi}$ | $\frac{1}{2 t^{2}}\left(-\partial_{\chi}\left(t^{2} \sin^{2}{\left(\theta \right)} \sinh^{2}{\left(\chi \right)}\right)\right)$ | $- \sin^{2}{\left(\theta \right)} \sinh{\left(\chi \right)} \cosh{\left(\chi \right)}$ |
| $\Gamma^{\theta}_{t \theta}$ | $\frac{1}{2 t^{2} \sinh^{2}{\left(\chi \right)}}\left(\partial_{t}\left(t^{2} \sinh^{2}{\left(\chi \right)}\right)\right)$ | $\frac{1}{t}$ |
| $\Gamma^{\theta}_{\chi \theta}$ | $\frac{1}{2 t^{2} \sinh^{2}{\left(\chi \right)}}\left(\partial_{\chi}\left(t^{2} \sinh^{2}{\left(\chi \right)}\right)\right)$ | $\frac{\cosh{\left(\chi \right)}}{\sinh{\left(\chi \right)}}$ |
| $\Gamma^{\theta}_{\phi \phi}$ | $\frac{1}{2 t^{2} \sinh^{2}{\left(\chi \right)}}\left(-\partial_{\theta}\left(t^{2} \sin^{2}{\left(\theta \right)} \sinh^{2}{\left(\chi \right)}\right)\right)$ | $- \frac{\sin{\left(2 \theta \right)}}{2}$ |
| $\Gamma^{\phi}_{t \phi}$ | $\frac{1}{2 t^{2} \sin^{2}{\left(\theta \right)} \sinh^{2}{\left(\chi \right)}}\left(\partial_{t}\left(t^{2} \sin^{2}{\left(\theta \right)} \sinh^{2}{\left(\chi \right)}\right)\right)$ | $\frac{1}{t}$ |
| $\Gamma^{\phi}_{\chi \phi}$ | $\frac{1}{2 t^{2} \sin^{2}{\left(\theta \right)} \sinh^{2}{\left(\chi \right)}}\left(\partial_{\chi}\left(t^{2} \sin^{2}{\left(\theta \right)} \sinh^{2}{\left(\chi \right)}\right)\right)$ | $\frac{\cosh{\left(\chi \right)}}{\sinh{\left(\chi \right)}}$ |
| $\Gamma^{\phi}_{\theta \phi}$ | $\frac{1}{2 t^{2} \sin^{2}{\left(\theta \right)} \sinh^{2}{\left(\chi \right)}}\left(\partial_{\theta}\left(t^{2} \sin^{2}{\left(\theta \right)} \sinh^{2}{\left(\chi \right)}\right)\right)$ | $\frac{1}{\tan{\left(\theta \right)}}$ |

曲率の計算では :eqref[eq-gr-curvature-components] の微分差と積差を分ける。第二列は $\partial_c\Gamma^a{}_{db}$、第三列は $-\partial_d\Gamma^a{}_{cb}$ であり、第四列は接続の積を符号付きで足したものである。三列を足して最後の列を得る。

| 成分 | 第一の微分 | 第二の微分 | 接続の積 | 結果 |
|---|---|---|---|---|
| $R^{\theta}_{\phi \theta \phi}$ | $- \cos{\left(2 \theta \right)}$ | $0$ | $\sin^{2}{\left(\theta \right)} \sinh^{2}{\left(\chi \right)} - \sin^{2}{\left(\theta \right)} \cosh^{2}{\left(\chi \right)} + \frac{\sin{\left(2 \theta \right)}}{2 \tan{\left(\theta \right)}}$ | $- \left(- \sinh^{2}{\left(\chi \right)} + \cosh^{2}{\left(\chi \right)} - 1\right) \sin^{2}{\left(\theta \right)}$ |

第一添字を計量で下ろすと $R_{abcd}=g_{ak}R^k{}_{bcd}$ となる。対の反対称性と交換対称性を用いて、非零の独立な組を下表に示す。

| 成分 | 結果 |
|---|---|
| $R_{\theta \phi \theta \phi}$ | $- t^{2} \left(- \sinh^{2}{\left(\chi \right)} + \cosh^{2}{\left(\chi \right)} - 1\right) \sin^{2}{\left(\theta \right)} \sinh^{2}{\left(\chi \right)}$ |

Ricci テンソルは第一・第三添字を縮約して得る。曲率が残る場合も相殺が起こるので、各項を落とさずに加える。

| 成分 | 縮約の項 | 代入 | 結果 |
|---|---|---|---|
| $R_{\theta \theta}$ | $R^{\phi}_{\theta \phi \theta}$ | $\sinh^{2}{\left(\chi \right)} - \cosh^{2}{\left(\chi \right)} + 1$ | $\sinh^{2}{\left(\chi \right)} - \cosh^{2}{\left(\chi \right)} + 1$ |
| $R_{\phi \phi}$ | $R^{\theta}_{\phi \theta \phi}$ | $- \left(- \sinh^{2}{\left(\chi \right)} + \cosh^{2}{\left(\chi \right)} - 1\right) \sin^{2}{\left(\theta \right)}$ | $- \left(- \sinh^{2}{\left(\chi \right)} + \cosh^{2}{\left(\chi \right)} - 1\right) \sin^{2}{\left(\theta \right)}$ |

逆計量と縮約すると Ricci スカラーとなる。さらに $G_{ab}=R_{ab}-g_{ab}R/2$ を成分ごとに計算する。

$$
\begin{aligned}
R&=\sum_{a,b}g^{ab}R_{ab}\\
&=- \frac{- \sinh^{2}{\left(\chi \right)} + \cosh^{2}{\left(\chi \right)} - 1}{t^{2} \sinh^{2}{\left(\chi \right)}} + \frac{\sinh^{2}{\left(\chi \right)} - \cosh^{2}{\left(\chi \right)} + 1}{t^{2} \sinh^{2}{\left(\chi \right)}}\\
&=- \frac{2 \left(- \sinh^{2}{\left(\chi \right)} + \cosh^{2}{\left(\chi \right)} - 1\right)}{t^{2} \sinh^{2}{\left(\chi \right)}}
\end{aligned}
$$

非零の Einstein テンソル成分は次のとおりである。

| 成分 | Ricci 項 | 跡の項 | 差 |
|---|---|---|---|
| $G_{t t}$ | $0$ | $\frac{- \sinh^{2}{\left(\chi \right)} + \cosh^{2}{\left(\chi \right)} - 1}{t^{2} \sinh^{2}{\left(\chi \right)}}$ | $- \frac{- \sinh^{2}{\left(\chi \right)} + \cosh^{2}{\left(\chi \right)} - 1}{t^{2} \sinh^{2}{\left(\chi \right)}}$ |
| $G_{\chi \chi}$ | $0$ | $- \frac{- \sinh^{2}{\left(\chi \right)} + \cosh^{2}{\left(\chi \right)} - 1}{\sinh^{2}{\left(\chi \right)}}$ | $\frac{- \sinh^{2}{\left(\chi \right)} + \cosh^{2}{\left(\chi \right)} - 1}{\sinh^{2}{\left(\chi \right)}}$ |

アフィン径数の微分を点で表し、$w_a$ を座標 $a$ の速度成分とする。ここでは添字を下げる操作を行っていない。測地線の右辺 $-\Gamma^a{}_{bc}w_bw_c$ では $b\ne c$ の項を二度数える。

$$
\begin{aligned}
\ddot{t}&=- t \left(w_{\chi}^{2} + w_{\phi}^{2} \sin^{2}{\left(\theta \right)} \sinh^{2}{\left(\chi \right)} + w_{\theta}^{2} \sinh^{2}{\left(\chi \right)}\right)\\
\ddot{\chi}&=\frac{t w_{\phi}^{2} \sin^{2}{\left(\theta \right)} \sinh{\left(\chi \right)} \cosh{\left(\chi \right)} + t w_{\theta}^{2} \sinh{\left(\chi \right)} \cosh{\left(\chi \right)} - 2 w_{\chi} w_{t}}{t}\\
\ddot{\theta}&=- \frac{4 t w_{\chi} w_{\theta} \cosh{\left(\chi \right)} - t w_{\phi}^{2} \sin{\left(2 \theta \right)} \sinh{\left(\chi \right)} + 4 w_{\theta} w_{t} \sinh{\left(\chi \right)}}{2 t \sinh{\left(\chi \right)}}\\
\ddot{\phi}&=- \frac{2 w_{\phi} \left(t w_{\chi} \tan{\left(\theta \right)} \cosh{\left(\chi \right)} + t w_{\theta} \sinh{\left(\chi \right)} + w_{t} \tan{\left(\theta \right)} \sinh{\left(\chi \right)}\right)}{t \tan{\left(\theta \right)} \sinh{\left(\chi \right)}}
\end{aligned}
$$

また計量ノルムは

$$
g(\dot x,\dot x)=t^{2} w_{\chi}^{2} + t^{2} w_{\phi}^{2} \sin^{2}{\left(\theta \right)} \sinh^{2}{\left(\chi \right)} + t^{2} w_{\theta}^{2} \sinh^{2}{\left(\chi \right)} - w_{t}^{2}
$$
である。これを径数で微分して上の加速度を代入すると、係数の微分と速度の微分の和が零になる。従って初期のノルムが保存される。

曲率成分は相殺してすべて零であり、$a''=0$、$(a')^2+k=0$ と一致する。$H=1/t$ は非零であっても、ここでは切片を定める観測者の束が広がることを表す。宇宙項に支配される de Sitter では $a\propto e^{Ht}$ と取れるうえに Ricci も非零なので、膨張の数値だけでは二つの幾何を同定できない。
:::

### de Sitter の平坦座標と静的座標

五次元の埋め込みで平坦な空間座標 $\vb x$ と時刻 $t$ を使い、
$$
X^0=L\sinh(t/L)+\frac{|\vb x|^2}{2L}e^{t/L},\quad
X^4=L\cosh(t/L)-\frac{|\vb x|^2}{2L}e^{t/L},\quad
X^i=e^{t/L}x^i
$$
と置く。$X^0+X^4=Le^{t/L}>0$ であり、$(X^4-X^0)(X^4+X^0)+|\vb X|^2=L^2$ である。微分して $-\dd X^{0\,2}+\dd X^{4\,2}=\dd(X^4+X^0)\dd(X^4-X^0)$ とすると、空間との交差項と $|\vb x|^2$ に比例する時間項が $\sum_i\dd X^{i\,2}$ の項と打ち消し合う。従って誘導計量は $-\dd t^2+e^{2t/L}\dd\vb x^2$ である。この座標に $X^0+X^4\le0$ の点は含まれず、$a\to0$ の過去極限でも曲率は一定である。

静的座標では $X^i=r\hat n^i$、$(X^0,X^4)=\sqrt{L^2-r^2}(\sinh(T/L),\cosh(T/L))$ と置く。この座標は $0\le r<L$、$X^4>|X^0|$ を覆い、誘導計量は
$$
\dd s^2=-(1-r^2/L^2)\dd T^2+(1-r^2/L^2)^{-1}\dd r^2+r^2\dd\Omega^2
$$
である。動径光では $\dd T=\pm\dd r/(1-r^2/L^2)$ なので、静的時刻で $r=L$ に近づくと対数発散になる。地平線を越える際に時空の曲率が発散するわけではないことは、埋め込みで確かめられる。$r=0$ の観測者に信号を送れる領域の未来境界は $X^0=X^4$ であり、$X^0+X^4>0$ の中では $r=L$ に当たる。平坦な膨張座標でも物理距離 $a(t)\int_t^\infty\dd t'/a(t')=L$ になる。

:::example[同じ定曲率でも二つの観測者で地平線が変わる]
de Sitter の大域座標の二つの共動観測者を $S^3$ の対蹠点に置く。$a=L\cosh(t/L)$ とし、$\sin\eta=\tanh(t/L)$ と置けば $\dd\eta=\dd t/a$、$-\pi/2<\eta<\pi/2$ である。ある観測者の世界線の未来端へ向かう光の共動距離は $\pi/2-\eta$ が上限であるため、二者の最短距離 $\pi$ には有限時刻の間に達しない。一方、各人に対して近傍には信号が届く。地平線は曲率の大小から決まるのではなく、どの世界線に向けた因果的過去を取るかで決まる。
:::

## 粒子の地平線と事象の地平線

共形時刻 $\dd\eta=\dd t/a(t)$ を使うと動径光は $\dd\chi/\dd\eta=\pm1$ である。下端 $t_{\rm i}$ と上端 $t_{\rm f}$ に対して

$$
\chi_{\rm p}(t)=\int_{t_{\rm i}}^t\frac{\dd t'}{a(t')},\qquad
\chi_{\rm e}(t)=\int_t^{t_{\rm f}}\frac{\dd t'}{a(t')}
$$

がそれぞれ :term[粒子の地平線]{en="particle horizon" yomi="りゅうしのちへいせん"} と宇宙論的事象の地平線の共動範囲を与える。ここでは非コンパクトな空間または最短動径距離を考える。前者はすでに受信できた範囲、後者はこれから発した信号が受信側に届き得る範囲である。$a=A t^p$、$0<p<1$ の模型では前者が有限で物理距離 $t/(1-p)$、後者は未来までこの膨張が続くと発散する。$a=e^{Ht}$ の未来では後者の物理距離は $a(t)\int_t^\infty e^{-Ht'}\dd t'=H^{-1}$ となる。Hubble 半径 $H^{-1}$ は局所的な膨張率から定まり、積分で定まる地平線とは一般に一致しない。

### Hubble 半径と二種類の地平線

冪乗 $a(t)=a_0(t/t_0)^p$、$p>0,t>0$ とする。粒子地平線の固有半径は $a(t)\int_0^t\dd t'/a(t')$ である。$p<1$ なら積分して $t/(1-p)$、$p\ge1$ なら下端が発散する。事象の地平線には $a(t)\int_t^\infty\dd t'/a(t')$ を使い、$p>1$ なら $t/(p-1)$、$p\le1$ なら上端が発散する。$p=1$ はいずれも対数発散である。これらは光が伝播した時刻の範囲に依存する。一方、$H^{-1}=t/p$ は時刻 $t$ での膨張率だけから作った長さである。

:::example[Hubble 半径を横切る内向き光]
平坦な切片で共動半径 $\chi$ の光の固有半径を $D=a\chi$ とすると、内向き光は $\dot D=HD-1$ に従う。例えば $p=2/3$ では粒子地平線は $3t$、Hubble 半径は $3t/2$ である。光が $D>H^{-1}$ から出発しても、拡大する Hubble 半径の内側へ入れば $\dot D<0$ に変わる。事象の地平線はこの模型では存在しないため、開始時の $\dot D>0$ だけで永遠に届かないとは言えない。
:::

### 共形時刻で光線の到達範囲を計算する

平坦な模型で $\dd\eta=\dd t/a$ とすれば光は $\chi-\chi_0=\pm(\eta-\eta_0)$ である。時刻 $\eta_1$ の共動観測者がそれまでに受信した信号を考えると、最初の時刻を $\eta_{\min}$ として半径の上限は $\eta_1-\eta_{\min}$ である。同じ観測者が将来放った信号がどこまで達するかには $\eta_{\max}-\eta_1$ を用いる。どちらも共動半径なので、特定の時刻の固有距離へ直すときはその切片の $a$ を掛ける。$\eta$ の両端が有限なら、同じ模型で二種類の地平線が成り立つ場合もある。

:::example[前半が減速で後半が加速する模型]
$t=t_*$ の前後で $a(t_*)=1$ とし、$0<t<t_*$ では $a=(t/t_*)^{2/3}$、$t>t_*$ では $a=(t/t_*)^2$ とする。$a$ は連続だが $\dot a$ が跳ぶため、これは光の積分を比較するための区分模型である。前半の共形時間の長さは $\int_0^{t_*}(t/t_*)^{-2/3}\dd t=3t_*$ である。後半の長さは $\int_{t_*}^\infty(t/t_*)^{-2}\dd t=t_*$ であり、合計は $4t_*$ になる。例えば $t=2t_*$ のとき、後半をすでに進んだ共形時間は $t_*/2$、残りも $t_*/2$ である。従って共動での粒子地平線と事象地平線はそれぞれ $7t_*/2$ と $t_*/2$ であり、その切片の $a=4$ を掛けた固有半径は $14t_*$ と $2t_*$ である。一方、$t\ge t_*$ の膨張率は $2/t$ なので同じ時刻の Hubble 半径は $t_*$ である。
:::

:::example[受信の有無は積分で判定する]
同じ模型で、$t=t_*$ のとき観測者と共動距離 $3t_*/4$ の点から光を内向きに放つ。$t$ までに進む共動距離は $t_*[1-t_*/t]$ である。この値が $3t_*/4$ になるには $t=4t_*$ が必要で、有限時刻に到達する。放射時の固有距離は $3t_*/4$ で、同時刻の Hubble 半径 $t_*/2$ より大きい。従って放射直後には固有半径が一度増えるが、やがて半径が減って到達する。これに対し初期の共動距離が $5t_*/4$ なら、未来の共形時間の全長 $t_*$ を超えるので到達しない。
:::

### 有限の共形時刻とコンパクトな空間

$k=1$ の単連結な球面切片で最短共動距離の上限は $\pi$ である。例えば閉じた塵模型では $0<\eta<2\pi$ である。$\eta=\pi/2$ の共動観測者が既に受け取れる信号の源は最短距離 $\chi<\pi/2$ に限られるが、未来に向けては $3\pi/2$ の共形時間が残るため、対蹠点にいる共動観測者にも有限時刻で送れる。逆に $\eta=3\pi/2$ ではすべての共動点から過去の信号が届き得るが、これからの信号の範囲は $\chi<\pi/2$ である。共形時間の積分が有限でも、その値が空間の直径より大きければ、最短経路に関する地平線を持たない時刻がある。ちょうど上限に一致する場合は、端の時刻を時空に含まないため対蹠点への到達は別に判定する。

## 演習問題

::::exercise{level=basic}
平坦で放射だけの宇宙において $a$ が三倍になるとエネルギー密度は何倍か。

:::solution
$w=1/3$ の保存則で $\rho a^4$ が一定なので $3^{-4}=1/81$ 倍である。
:::
::::

::::exercise{level=basic}
$a_{\rm e}=a_{\rm o}/4$ の共動光源の $z$ を求めよ。

:::solution
:eqref[eq-gr-cosmo-redshift] で $1+z=4$ なので $z=3$ である。
:::
::::

::::exercise{level=basic}
$z=1$、$D_A=1\,\mathrm{Gpc}$ の天体の $D_L$ を求めよ。

:::solution
:eqref[eq-gr-distance-duality] から $D_L=(1+1)^2D_A=4\,\mathrm{Gpc}$ である。
:::
::::

::::exercise{level=basic}
平坦な物質のみのモデルで $t=t_0$ の $H$ を求めよ。

:::solution
$a=A t^{2/3}$ なので $H=\dot a/a=2/(3t_0)$ である。
:::
::::

::::exercise{level=advanced}
$a=A t^{2/3}$ の粒子の地平線の物理距離を求めよ。

:::solution
共動距離は $\int_0^t\dd t'/(A t'^{2/3})=3t^{1/3}/A$ である。$a(t)$ を掛けると $3t$ となる。
:::
::::

::::exercise{level=advanced}
$a=A t^2$ が未来まで続くモデルで事象の地平線と Hubble 半径を比較せよ。

:::solution
$a(t)\int_t^\infty\dd t'/(At'^2)=t$ である。一方 $H=2/t$、$H^{-1}=t/2$ なので事象の地平線はその二倍である。
:::
::::

::::exercise{level=advanced}
平坦な物質と宇宙項で現在 $\rho_{m0}/\rho_\Lambda=1/3$ とする。加速の開始赤方偏移を求めよ。

:::solution
開始は $\rho_m=2\rho_\Lambda$ なので $(1+z)^3/3=2$ から $z=6^{1/3}-1$ となる。
:::
::::

::::exercise{level=advanced}
$a=L\cosh(t/L),k=1$ が真空の Friedmann 方程式を満たすことを確認せよ。

:::solution
$H=L^{-1}\tanh(t/L)$ で $H^2+a^{-2}=L^{-2}(\tanh^2+\operatorname{sech}^2)=L^{-2}$ である。また $\ddot a/a=L^{-2}$ である。従って $\Lambda=3/L^2$ として両式を満たす。
:::
::::

::::exercise{level=advanced}
平坦な物質のみの宇宙で、放射時刻が現在の $1/8$ の光源を考える。赤方偏移、共動距離、光度距離、角径距離を求めよ。現在の尺度因子を $1$ とする。

:::solution
$a(t)=(t/t_0)^{2/3}$ なので $a_e=(1/8)^{2/3}=1/4$、$1+z=4$ である。共動距離は
$$
D_C=\int_{t_e}^{t_0} (t/t_0)^{-2/3}\dd t
=3t_0[1-(t_e/t_0)^{1/3}]=3t_0/2.
$$
$H_0=2/(3t_0)$ だから $D_C=H_0^{-1}$ である。平坦模型では $D_L=(1+z)D_C=4H_0^{-1}$、$D_A=D_C/(1+z)=(4H_0)^{-1}$ となる。両者の比は $16=(1+z)^2$ で、独立に計算して距離双対性を検算できる。
:::
::::

## 確認テスト

::::quiz
光度距離と角径距離の差に寄与する要素を選べ。

- [x] 受光一個当たりのエネルギーの変化
- [x] 到着率の変化
- [ ] 光子が一般に静止質量を持つこと

:::explanation
光子数保存と零伝播を用い、エネルギーと到着率にそれぞれ $(1+z)^{-1}$ が掛かる。さらに放射後の面積半径は放射時より $1+z$ 倍になるので $D_L=(1+z)^2D_A$ を得る。
:::
::::

::::quiz{answer="3" tol="1e-9"}
$a_{\rm o}/a_{\rm e}=4$ のとき赤方偏移を答えよ。

:::explanation
$1+z=a_{\rm o}/a_{\rm e}=4$ なので $z=3$ である。
:::
::::

::::quiz
事象の地平線の判定に直接必要な情報を選べ。

- [x] 未来にわたる $\int\dd t/a(t)$
- [ ] 現在の $a(t)$ だけ
- [ ] 現在の物質密度だけ

:::explanation
光の共動移動距離は $\int\dd t/a(t)$ である。非コンパクトな空間では、未来の積分が有限なら、その値以上の共動距離を隔てた点の間で今から発した光は有限の未来時刻に届かない。コンパクトな空間では最短距離の上限とも比較する。
:::
::::

::::quiz{answer="3" tol="1e-9"}
$k=0$ の純粋な de Sitter 膨張で $\Lambda/H^2$ を求めよ。

:::explanation
物質密度が零なので :eqref[eq-gr-friedmann] から $H^2=\Lambda/3$ となる。従って比は $3$ である。
:::
::::

::::quiz{answer="0.25" tol="1e-7"}
平坦で物質だけの宇宙において $z=3$ の $H_0D_A$ を答えよ。

:::explanation
積分で $H_0D_C=2[1-(1+z)^{-1/2}]=1$ である。角径距離はこれを $1+z=4$ で割るので $H_0D_A=1/4$ となる。
:::
::::

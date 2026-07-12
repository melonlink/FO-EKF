# 心率激励下分数阶 ECG 动力学的可辨识性与鲁棒观测：理论规格 v0.1

- **冻结日期**：2026-07-11
- **工作分支**：`dev-codex`
- **研究阶段**：P0 模型与理论主线冻结；尚未进入性能实现
- **适用对象**：表面 ECG 动力学代理（ECG dynamical surrogate / digital shadow）
- **不适用称谓**：当前不称 cardiac digital twin，不把代理状态解释为真实心肌内部状态

> **S0--S3 后审计修正（2026-07-11）**：本文中的三点“精确复圆盘可行集”只对差商代数关系精确，不是强制正实公共阻尼、全部心率/谐波共享物理参数的完整可行集；三心率也尚未证明为最小实验。完整修订路线及两心率物理可辨识推导见 `theory_strengthening_decision_2026-07-11.md`。在联合前向圆盘集合和有限记录到响应圆盘的桥定理完成前，T0 不视为关闭。

## 1. 决策结论

原计划的“分数阶 McSharry + FO-EKF + 一般 Mittag-Leffler 稳定性”不能作为论文创新主线：分数阶 ECG、状态/参数/阶次联合估计、完整历史 Kalman 滤波、有限记忆误差观测器以及分数阶 LMI 观测器均已有直接先行工作。

本项目改为以下理论问题：

> 在未知 ECG 形态幅值、未知阻尼、未知公共导联增益和公共相位规范的条件下，能否利用同一受试者在多个稳定心率下的同一形态谐波，全局唯一地辨识固定分数阶阶次，并把有限驻留、形态漂移、有限记忆、测量污染和学习残差逐项传播到可计算的观测误差管？

总判定为 **CONDITIONAL GO**。可守的完整贡献链是：

1. 用已有的非周期性定理排除“把完整 ECG 极限环直接 Caputo 分数阶化”的错误架构；
2. 固定 McSharry 相位规范并删除表面 ECG 无法观测的径向伪状态；
3. 推导有限下限 Caputo 形态通道的渐近周期响应和有限驻留误差界；
4. 构造三心率 gauge-invariant 逆频响差商，证明其对阶次全局单射；
5. 给出噪声、形态漂移和心率间隔下的显式阶次误差界与更新/拒绝门；
6. 设计 approximation-aware 分数阶观测器，给出 Mittag-Leffler ISS 误差管；
7. 用 full-memory 合成真值先证伪，再用 Fantasia 做探索性受试者内多心率检验，用 PTB-XL/Challenge 做预测与鲁棒性辅助验证。

这里的创新是这条 **端到端、可证伪的组合链**，不是其中任意一个标准工具。

## 2. 文献碰撞后保留与删除的主张

| 候选主张 | 审查 | 处理 |
|---|---|---|
| 首个分数阶 ECG / 分数阶 McSharry | Takha et al. 2024 已只对第三个 PQRST 方程分数阶化并做 MIT-BIH 拟合 | 删除 |
| 首个 FO-EKF 或状态—参数—阶次联合估计 | Sierociuk 2006、Sierociuk & Macias 2021 等已覆盖 | 删除 |
| 首个完整跨时协方差 FKF | Sadeghian & Salarieh 2013 已用 cumulative vector 构造一般线性 FKF | 删除 |
| 首个有限记忆截断误差观测/区间界 | Rauh & Lahme 2025 已覆盖点值、集合值和区间观测 | 删除 |
| 首个 Caputo observer + LMI + UUB | Naifar & Ben Makhlouf 2021 等已有 practical Mittag-Leffler 稳定性 | 删除 |
| 首个分数阶系统频域可辨识性 | Nazarian et al. 2010、Alavi et al. 2015 等已有一般理论 | 删除 |
| ECG 相位 gauge 与 McSharry 径向不可观的形式化商空间分析 | 未发现针对该模型的完整定理；Sameni 2007 已实际采用 phase-only 降阶 | 作为支撑贡献，不声称首次降阶 |
| 三心率逆频响差商对 \(\alpha\) 的全局单射、恢复式和模型拒绝检验 | 未发现 ECG 直接先例；单独数学引理偏薄 | 与有限驻留、漂移界和 ISS 观测器合并为核心 |
| 近似误差逐项进入最终观测误差管 | 通用 ISS/LMI 已有，ECG 混合阶结构下的端到端映射仍可守 | 核心工程理论贡献 |

正式投稿前仍需做一次按标题、公式和引用网络展开的系统检索；当前只允许使用窄限定语，不使用未经证实的 `first`。

## 3. 建议论文题目与贡献措辞

建议题目：

> **Gauge-Invariant Multi-Rate Identifiability and Robust Observation of Mixed Integer-Fractional ECG Morphology Dynamics**

更偏应用的题目：

> **Heart-Rate-Excited Identification and Approximation-Aware Observation of Fractional ECG Dynamics**

可守贡献措辞：

> We derive a gauge-invariant multi-rate certificate that globally identifies the fixed fractional order of a stable ECG morphology channel despite unknown damping, morphology amplitude, common lead gain, and common phase gauge. Finite dwell, morphology drift, measurement error, and memory approximation are then propagated into explicit order and Mittag-Leffler observer-error bounds, with rejection gates that prevent unsupported online order updates.

不再把 `FO-EKF` 放在题目或首条贡献中。FO-EKF 只可作为概率层、比较基线或后续扩展。

## 4. 模型定义与量纲

### 4.1 标量可辨识性模型

固定参考时间 \(\tau_0>0\)，定义无量纲时间 \(s=t/\tau_0\)。对第 \(r\) 个近似恒定心率区间，采用

\[
\frac{\mathrm d\phi}{\mathrm ds}=\omega_r,\qquad
{}^C D_{0+}^{\alpha}z(s)=-\lambda z(s)+u(\phi(s))+d(s),
\]

\[
y(s)=g z(s)+b+v(s),
\]

其中：

- \(0<\alpha\le 1\) 为校准和观测区间内固定的常阶；
- \(\lambda>0\) 为无量纲稳定系数；
- \(u(\phi+2\pi)=u(\phi)\) 为 ECG 形态驱动；
- \(g\ne0\) 是固定导联增益，\(b\) 是基线；
- \(d,v\) 分别是有界模型失配和测量污染；
- \(z\) 是表面 ECG 形态代理状态，不是心肌电位、传导速度或兴奋性真值。

对固定谐波 \(m\ge1\)，令 \(q_m\ne0\) 为 \(u\) 的复 Fourier 系数，并定义该谐波的无量纲角频率

\[
\nu_r=m\tau_0\Omega_r=m\omega_r>0.
\]

以下省略下标 \(m\)，写成 \(q\) 与 \(\nu_r\)。固定 \(\tau_0\) 是必要的；若把参考时间也作为未知量，则只能恢复时间缩放组合，不能单独解释 \(\lambda\) 和 \(q\)。

### 4.2 向量观测器扩展

观测器使用更一般的形态状态：

\[
{}^C D^{\alpha}z
=Az+B_\phi g_\phi(\phi)+B_r r^\star(\chi)+B_d d,
\]

\[
y=Cz+c_\phi(\phi)+v.
\]

该扩展用于误差证书；三心率闭式辨识定理首先作用于已固定相位规范的标量可观测商通道。多状态/多导联的结构可辨识性必须另做秩检验，不能由标量定理自动推出。

## 5. 假设清单

后续每个定理必须显式引用下列假设，不能把条件藏在实现中。

- **A1 固定阶次**：\(\alpha\in I=[\alpha_-,\alpha_+]\subset(0,1]\) 在一次校准/测试段内不变；本论文不证明变阶系统。
- **A2 固定导数定义**：连续理论使用固定下限的 Caputo 导数；严格稳态频响用等价的无限历史/Weyl 响应描述，并显式保留有限下限瞬态。
- **A3 稳定形态通道**：标量模型 \(\lambda>0\)；向量模型具有共同的 LMI 证书。
- **A4 公共相位规范**：所有区间均用 R 峰事件固定 \(\phi(t_R)=0\)，且预处理群延迟已经统一校准。
- **A5 同一形态实验**：用于同一三元组的 \(q,\lambda,\alpha,g\) 在三种心率下相同，或其变化有独立上界。不得对各区间分别做幅值归一化。
- **A6 充分驻留**：每种心率保持足够时间，有限下限瞬态和前一心率历史污染进入已量化的 \(\varepsilon_r(T_r)\)。
- **A7 复谐波可测**：使用带共同相位参考的复 Fourier 系数；只有功率或幅值不足以直接使用主定理。
- **A8 非退化激励**：\(0<\nu_1<\nu_2<\nu_3\)，且逆频响差 \(|W_2-W_1|\) 高于噪声阈值。
- **A9 有界输入语义**：ISS 定理中的噪声、残差和截断项是有界分段连续函数，不是理想白噪声。
- **A10 同一信息边界**：所有预测基线只能使用预测时刻以前的数据，真实数据不得用未来零相位滤波或测试段全局归一化。

## 6. 架构约束：为什么不能把整个极限环直接分数阶化

### 命题 0（已有定理的 ECG 架构推论，不列为原创）

在相关正则性和有限下限条件下，自治、时不变、非整数阶 Caputo 系统不存在非恒定精确周期解。因此，把 McSharry 的相位极限环整体写成自治 Caputo 系统，不能同时保留非恒定精确周期轨道。

这不是本论文的新定理。Tavazoei & Haeri 2009 已证明一般非周期性；Yazdani & Salarieh 2011 说明在下限趋于负无穷的稳态意义下可检测周期响应。Takha et al. 2024 已在数值上保留前两个整数阶相位方程、只分数阶化第三个形态方程。

本项目据此采用三角结构：

\[
\text{整数阶相位/心率}\longrightarrow\text{稳定 Caputo 形态通道}.
\]

即使采用这一结构，有限下限 Caputo 解也不是从 \(s=0\) 起严格周期；正确表述是“趋向无限历史定义的周期稳态响应”。

## 7. 相位 gauge 与最小可观测商空间

### 命题 1（相位 gauge）

对任意 \(\delta\in S^1\)，变换

\[
\phi'(s)=\phi(s)+\delta,\qquad
q_m'=q_m e^{-im\delta}
\]

不改变 \(u(\phi)\)、\(z\) 和全部输出。McSharry 的 Gaussian 参数形式中对应

\[
\phi'=\phi+\delta,\qquad \theta_i'=\theta_i+\delta.
\]

所以未固定相位原点时，初始相位与全部波形位置参数不具结构可辨识性，灵敏度矩阵至少有一维零空间。

**规范固定**：

\[
\theta_R=0,qquad \phi(t_R)=0,qquad \omega>0,
\]

并固定 P-Q-R-S-T 标签顺序、\(b_i>0\)、R 峰局部唯一性和无混叠心率范围。

### 命题 2（McSharry 径向状态不可观）

经典前两个状态写成极坐标后为

\[
\dot r=(1-r)r,qquad \dot\phi=\omega,
\]

而形态方程和表面 ECG 输出只依赖 \((\phi,z)\)，不依赖 \(r\)。相同 \(\phi(0),z(0)\) 下，任意不同 \(r(0)>0\) 产生同一表面 ECG，故径向状态结构不可观。

因此首篇论文直接令 \(r=1\)，使用 phase-only 模型；不在 FO-EKF 中估计笛卡尔 \(x,y\) 的径向自由度。Sameni 2007 已实际采用 \((\phi,z)\) 降阶，本项目的贡献是把该选择写成明确的可观测商空间条件，而不是声称首次降阶。

## 8. 有限下限 Caputo 响应与有限驻留误差

考虑复数化标量通道

\[
{}^C D_{0+}^{\alpha}z(s)=-\lambda z(s)+u(s),\qquad z(0)=z_0,
\]

其中 \(0<\alpha\le1,\lambda>0\)，\(u\) 有界。定义 Mittag-Leffler 核

\[
k_\alpha(\sigma)=\sigma^{\alpha-1}E_{\alpha,\alpha}(-\lambda\sigma^\alpha).
\]

对 \(0<\alpha\le1\)，\(k_\alpha\ge0\)，且

\[
\int_T^\infty k_\alpha(\sigma)\,\mathrm d\sigma
=\frac{E_\alpha(-\lambda T^\alpha)}{\lambda}.
\]

### 定理 1（渐近周期响应）

若 \(u_p\) 是有界 \(T_p\)-周期函数，定义无限历史稳态

\[
z_{\mathrm{ss}}(s)=\int_0^\infty k_\alpha(\sigma)u_p(s-\sigma)\,\mathrm d\sigma.
\]

则 \(z_{\mathrm{ss}}\) 是 \(T_p\)-周期的；有限下限解

\[
z(s)=z_0E_\alpha(-\lambda s^\alpha)
+\int_0^s k_\alpha(\sigma)u_p(s-\sigma)\,\mathrm d\sigma
\]

满足

\[
|z(s)-z_{\mathrm{ss}}(s)|
\le
\left(|z_0|+\frac{\|u_p\|_\infty}{\lambda}\right)
E_\alpha(-\lambda s^\alpha).
\]

当 \(0<\alpha<1\) 时，右侧为 \(O(s^{-\alpha})\)，不是指数衰减。

**证明**：周期性由卷积平移直接得到。有限解与无限历史解之差等于初值项减去 \([s,\infty)\) 的核尾；用核非负性和上述尾积分恒等式即可。

### 推论 1.1（心率切换后的历史污染界）

设在时刻 \(S\) 切换到新周期输入 \(u_{\mathrm{new}}\)，在此之前任意历史输入满足 \(\|u_{\mathrm{pre}}\|_\infty\le U_{\mathrm{pre}}\)，新输入满足 \(\|u_{\mathrm{new}}\|_\infty\le U_{\mathrm{new}}\)。在切换后驻留 \(T\) 时刻，固定下限 Caputo 解与“新输入从负无穷已存在”的周期稳态之间满足

\[
|z(S+T)-z_{\mathrm{ss,new}}(S+T)|
\le
\left(|z_0|+\frac{U_{\mathrm{pre}}+U_{\mathrm{new}}}{\lambda}\right)
E_\alpha(-\lambda T^\alpha).
\]

证明把 \([0,S]\) 的旧输入卷积和新稳态的 \([T,\infty)\) 核尾分别界定。该界保守但显式说明：小 \(\alpha\) 时驻留要求很长，不能把短 ECG 片段直接当稳态。

若

\[
u_p(s)=\sum_m q_m e^{im\omega s},
\]

则稳态 Fourier 系数为

\[
Z_m=\frac{q_m}{\lambda+(im\omega)^\alpha},
\]

主值支路取

\[
(i\nu)^\alpha=\nu^\alpha e^{i\pi\alpha/2},\qquad \nu>0.
\]

有限驻留估计必须写成

\[
\widehat Z_r(T_r)=Z_r+\varepsilon_r(T_r),
\]

其中 \(\varepsilon_r\) 同时包含上述 Mittag-Leffler 瞬态、有限窗谱泄漏、形态漂移、预处理和测量误差。

## 9. 核心定理：三心率商不变量全局辨识阶次

### 定理 2（全局单射与闭式恢复）

取同一受试者、同一导联、同一形态谐波在三个稳态心率下的复系数：

\[
Z_r=\frac{q}{\lambda+(i\nu_r)^\alpha},
\qquad
0<\nu_1<\nu_2<\nu_3,
\]

其中 \(q\ne0,\lambda>0,\alpha>0\) 跨三个区间不变。令

\[
W_r=\frac1{Z_r},\qquad
R=\frac{W_3-W_1}{W_2-W_1}.
\]

则

\[
R=F(\alpha)
=\frac{\nu_3^\alpha-\nu_1^\alpha}
       {\nu_2^\alpha-\nu_1^\alpha},
\]

且 \(F\) 在 \(\alpha>0\) 上严格递增，因此 \(R\) 全局唯一确定 \(\alpha\)。

**证明**：令

\[
a=\log(\nu_2/\nu_1),\qquad
b=\log(\nu_3/\nu_1),\qquad0<a<b.
\]

则

\[
F(\alpha)=\frac{e^{b\alpha}-1}{e^{a\alpha}-1},
\]

且

\[
\frac{\mathrm d}{\mathrm d\alpha}\log F(\alpha)
=\frac{b}{1-e^{-b\alpha}}-
  \frac{a}{1-e^{-a\alpha}}.
\]

对固定 \(\alpha>0\)，函数

\[
h(x)=\frac{x}{1-e^{-\alpha x}}
\]

严格递增，因为

\[
h'(x)=
\frac{1-(1+\alpha x)e^{-\alpha x}}
     {(1-e^{-\alpha x})^2}>0,
\]

最后一步来自 \(e^{\alpha x}>1+\alpha x\)。所以 \(F'(\alpha)>0\)。证毕。

若 \(0<\alpha\le1\)，理想数据还必须满足

\[
\frac{\log(\nu_3/\nu_1)}{\log(\nu_2/\nu_1)}
<R\le
\frac{\nu_3-\nu_1}{\nu_2-\nu_1},
\]

且 \(R\) 必须为正实数。

得到 \(\alpha\) 后，令

\[
B=\frac{W_2-W_1}{\nu_2^\alpha-\nu_1^\alpha},
\qquad
A=W_1-B\nu_1^\alpha,
\]

即可恢复

\[
q=\frac{e^{i\pi\alpha/2}}{B},
\qquad
\lambda=Aq.
\]

未知但跨三个区间恒定的导联增益 \(g\) 只会把 \(q\) 替换为 \(gq\)。因此 \(\alpha,\lambda\) 仍可恢复，而物理 \(q\) 与 \(g\) 不能分开。

### 为什么该量是商空间不变量

\(W_r\) 可写为

\[
W_r=\frac{\lambda}{q}+
\frac{e^{i\pi\alpha/2}}{q}\nu_r^\alpha=A+B\nu_r^\alpha.
\]

差分先消除平移 \(A=\lambda/q\)，差分比再消除公共复缩放 \(B\)。因此它同时对未知形态幅值、公共导联增益和公共相位旋转不变。

## 10. 有限噪声阶次界与更新门

### 定理 3（逆频响扰动传播）

设 \(|\widehat Z_r-Z_r|\le\epsilon_{Z,r}<|Z_r|\)。则

\[
|\widehat W_r-W_r|
\le
\frac{\epsilon_{Z,r}}
{|Z_r|(|Z_r|-\epsilon_{Z,r})}
=:\epsilon_{W,r}.
\]

令 \(\epsilon_W=\max_r\epsilon_{W,r}\)、\(D=W_2-W_1\)、\(N=W_3-W_1\)。若 \(2\epsilon_W<|D|\)，则

\[
|\widehat R-R|
\le
\frac{2\epsilon_W(|D|+|N|)}
{|D|(|D|-2\epsilon_W)}
=:\epsilon_R.
\]

在紧区间 \(I=[\alpha_-,\alpha_+]\subset(0,1]\) 上令

\[
m_I=\inf_{\alpha\in I}F'(\alpha)>0.
\]

投影反演 \(\widehat\alpha=F^{-1}(\Pi_{F(I)}\operatorname{Re}\widehat R)\) 满足

\[
|\widehat\alpha-\alpha|
\le\frac{\epsilon_R}{m_I}.
\]

**证明**：第一式由倒数差的精确恒等式得到；第二式由复商扰动恒等式和 \(|D+\delta D|\ge|D|-2\epsilon_W\) 得到；第三式使用均值定理。

### S3 执行后补充：可由观测量直接计算的后验闭式界

上式适合真值审计。实际算法只有观测复系数 \(\widetilde Z_r\) 和误差半径 \(\epsilon_r\)。若

\[
|\widetilde Z_r|>\epsilon_r,
\]

则

\[
|\widetilde W_r-W_r|
\le b_r=
\frac{\epsilon_r}
{|\widetilde Z_r|(|\widetilde Z_r|-\epsilon_r)}.
\]

令

\[
\widetilde B=\widetilde W_2-\widetilde W_1,
\qquad
\widetilde R=\frac{\widetilde W_3-\widetilde W_1}{\widetilde B}.
\]

若 \(|\widetilde B|>b_1+b_2\)，则共享的 \(W_1\) 误差可保留其正确系数，得到

\[
|\widetilde R-R|
\le
\epsilon_R^{\mathrm{chain}}
=
\frac{
b_3+|\widetilde R|b_2+|1-\widetilde R|b_1
}{
|\widetilde B|-b_1-b_2
}.
\]

若响应圆盘包含零或差分分母条件失败，该三元组不是“误差较大但仍可用”，而是没有有限倒数/差商证书，必须拒绝。

进一步令

\[
a=\log(\nu_2/\nu_1),\qquad
b=\log(\nu_3/\nu_1),\qquad0<a<b,
\]

并写成

\[
F'(\alpha)=F(\alpha)g(\alpha),
\qquad
g(\alpha)=
\frac b{1-e^{-b\alpha}}-
\frac a{1-e^{-a\alpha}}.
\]

因为

\[
g'(\alpha)=
\left[\frac{a}{2\sinh(a\alpha/2)}\right]^2-
\left[\frac{b}{2\sinh(b\alpha/2)}\right]^2>0,
\]

所以 \(F''(\alpha)=F(\alpha)[g(\alpha)^2+g'(\alpha)]>0\)。因此可严格取

\[
m_I=F'(\alpha_-),
\]

无需通过数值网格猜测 \(F'\) 的最小点。

### S3 执行后补充：三点差商的复圆盘外松弛

倒数映射可精确处理。若

\[
Z_r\in\mathcal D(\widetilde Z_r,\epsilon_r),
\qquad |\widetilde Z_r|>\epsilon_r,
\]

则

\[
W_r\in\mathcal D(c_r,s_r),
\]

\[
c_r=
\frac{\overline{\widetilde Z_r}}
{|\widetilde Z_r|^2-\epsilon_r^2},
\qquad
s_r=
\frac{\epsilon_r}
{|\widetilde Z_r|^2-\epsilon_r^2}.
\]

对候选实数 \(R>1\)，三点圆盘满足线性差商关系的闭包条件为

\[
|c_3-c_1-R(c_2-c_1)|
\le s_3+(R-1)s_1+Rs_2.
\]

代入 \(R=F(\alpha)\) 得到一维代数相容集合：

\[
\mathcal A_{\mathrm{disk}}
=\left\{\alpha\in I:
|c_3-c_1-F(\alpha)(c_2-c_1)|
\le s_3+[F(\alpha)-1]s_1+F(\alpha)s_2
\right\}.
\]

该集合没有强制由同一正实 \(\lambda\) 和公共 \(q\) 产生，也没有强制多个三元组/谐波共享同一组潜变量，因此只是完整物理阶次集合的外松弛。实现必须保留所有连通分量；固定网格只能给近似结果，不能称完备算法。倒数圆盘含零或差分无法分离时只能判定该快速证书不可用，不能据此判定完整物理模型无解。

### 形态和阻尼漂移如何进入 \(\epsilon_Z\)

对 \(d_r=\lambda+(i\nu_r)^\alpha\)，若 \(q_r=q+\delta q_r\)、\(\lambda_r=\lambda+\delta\lambda_r\) 且 \(|\delta\lambda_r|<|d_r|\)，则

\[
\left|
\frac{q+\delta q_r}{d_r+\delta\lambda_r}
-\frac{q}{d_r}
\right|
\le
\frac{|\delta q_r||d_r|+|q||\delta\lambda_r|}
{|d_r|(|d_r|-|\delta\lambda_r|)}.
\]

该项必须与有限驻留、谱估计和测量误差相加，成为可验证的 \(\epsilon_{Z,r}\)。若没有 \(\delta q_r,\delta\lambda_r\) 的可信上界，就不能在线更新 \(\alpha\)。

### 必须同时通过的阶次更新门

1. **G0 协议门**：公共相位锚、共同 \(q,\lambda\) 或其独立漂移上界、禁止逐段归一化、延迟已校准、驻留预算已计算；
2. **驻留门**：定理 1 的瞬态上界进入总 \(\epsilon_r\) 后仍可形成有限阶次区间；
3. **倒数门**：所有 \(|\widetilde Z_r|>\epsilon_r\)；
4. **心率分离门**：\(|\widetilde W_2-\widetilde W_1|>b_1+b_2\)；
5. **复相位门**：\(|\operatorname{Im}\widetilde R|\le\epsilon_R^{\mathrm{chain}}+\tau_{\mathrm{num}}\)；
6. **范围门**：\(\operatorname{dist}(\operatorname{Re}\widetilde R,F(I))\le\epsilon_R^{\mathrm{chain}}+\tau_{\mathrm{num}}\)；
7. **条件数门**：主论文要求 \(\epsilon_R^{\mathrm{chain}}/F'(\alpha_-)\le0.02\)，探索性上限为 0.05；
8. **多三元组门**：五个以上心率的不同三元组给出相容阶次；
9. **多谐波门**：至少两个非零形态谐波给出相容阶次；
10. **漂移门**：允许 \(q_r,\lambda_r\) 漂移的竞争模型没有显著优于共同参数模型；
11. **预测门**：冻结参数后能改善留一心率的更新前预测。

门控结果使用三类语义：`PASS`、`REJECT_NUMERIC` 和 `NOT_CERTIFIABLE`。`REJECT_NUMERIC` 只用于差商外松弛已经与误差预算不相容的情况；G0 失败、倒数圆盘含零、差分无法分离或非空阶次集合过宽都属于 `NOT_CERTIFIABLE`。后者不能据此宣称完整代理模型无解。任一关键门失败时冻结 \(\alpha\)；不是调低阈值直到得到可用数值。

## 11. 反例与证伪条件

### 反例 1：每个心率独立形态

若每个区间允许任意 \(q_r\)，则对任意候选 \(\alpha,\lambda\) 都可取

\[
q_r=Z_r[\lambda+(i\nu_r)^\alpha],
\]

所以 \(\alpha\) 完全不可辨识。

### 反例 2：独立相位旋转

共同相位旋转会并入公共 \(q\) 并被差商消掉；每个区间独立的 \(\delta_r\) 会令

\[
\widetilde Z_r=e^{-i\delta_r}Z_r,
\]

从而破坏不变量。R 峰相位对齐和预处理延迟校准是定理条件，不是普通预处理细节。

### 反例 3：幅值逐段归一化

对每个心率区间单独做 min-max 或方差归一化等价于引入独立增益 \(g_r\)，会破坏公共复缩放。Takha 2024 的逐搏归一化流程不能直接用于本定理验证。

### 反例 4：短驻留与错误离散频响

有限下限 Caputo 解在短区间内包含初始化和前一心率历史；有限记忆 GL 的离散频响也不是精确 \((i\nu)^\alpha\)，而是

\[
h^{-\alpha}(1-e^{-i\nu h})^\alpha
=(i\nu)^\alpha e^{-i\alpha\nu h/2}
\left[\frac{2\sin(\nu h/2)}{\nu h}\right]^\alpha.
\]

因此闭式辨识必须先在 full-memory/连续稳态层验证，采样与截断误差另行计入 \(\epsilon_Z\)。

### 反例 5：接近整数阶时的实际病态

Takha 2024 报告的优化阶次约为 0.96–0.99。该区间中常用心率组合对 \(\alpha\) 的差商变化可能非常小；结构可辨识不代表实际可辨识。必须报告 \(F'\)、阶次区间和误差上界，不能只报告优化器点估计。

## 12. Approximation-aware Mittag-Leffler ISS 观测器

观测器定义为

\[
{}^C D^\alpha\widehat z
=A\widehat z+B_\phi g_\phi(\widehat\phi)
+B_r r_\psi(\widehat\chi)
+L[y-C\widehat z-c_\phi(\widehat\phi)].
\]

若相位由 R 峰直接给定，\(\widehat\phi=\phi\)，误差 \(e=z-\widehat z\) 为

\[
{}^C D^\alpha e=(A-LC)e+w,
\]

\[
w=B_r\varepsilon_\psi+B_dd+B_mr_{\mathrm{mem}}-Lv,
\qquad
\varepsilon_\psi=r^\star-r_\psi.
\]

### 定理 4（显式 Mittag-Leffler ISS 误差管）

给定 \(\rho>0\)。若存在 \(P=P^\top\succ0,Y\) 满足

\[
A^\top P+PA-C^\top Y^\top-YC\preceq-2\rho P,
\]

取 \(L=P^{-1}Y\)，且 \(\|w(s)\|_P\le\overline w_P\)，则令

\[
m_\alpha(s)=E_\alpha(-\rho s^\alpha),
\]

有

\[
\boxed{
\|e(s)\|_P
\le
m_\alpha(s)\|e(0)\|_P
+\frac{\overline w_P}{\rho}[1-m_\alpha(s)]
}
\]

以及

\[
\boxed{
\limsup_{s\to\infty}\|e(s)\|_P
\le\frac{\overline w_P}{\rho}.
}
\]

欧氏范数界为

\[
\|e(s)\|_2
\le
\sqrt{\kappa(P)}m_\alpha(s)\|e(0)\|_2
+\frac{\overline w_P}{\rho\sqrt{\lambda_{\min}(P)}}
[1-m_\alpha(s)].
\]

扰动预算可逐项写为

\[
\overline w_P
\le b_r\overline\varepsilon_\psi
+b_d\overline d
+b_m\overline r_{\mathrm{mem}}
+b_v\overline v,
\]

\[
b_i=\sqrt{\lambda_{\max}(B_i^\top PB_i)},
\qquad
b_v=\sqrt{\lambda_{\max}(L^\top PL)}.
\]

**证明骨架**：LMI 先给出经典半群的 \(P\)-收缩

\[
\|e^{(A-LC)\sigma}\|_P\le e^{-\rho\sigma}.
\]

利用 Caputo 解算子的正 subordination 表示，得到

\[
\|E_\alpha((A-LC)s^\alpha)\|_P
\le E_\alpha(-\rho s^\alpha),
\]

\[
\|s^{\alpha-1}E_{\alpha,\alpha}((A-LC)s^\alpha)\|_P
\le s^{\alpha-1}E_{\alpha,\alpha}(-\rho s^\alpha).
\]

代入变常数公式并使用

\[
\int_0^s\sigma^{\alpha-1}E_{\alpha,\alpha}(-\rho\sigma^\alpha)\,\mathrm d\sigma
=\frac{1-E_\alpha(-\rho s^\alpha)}{\rho}
\]

即可。\(\alpha=1\) 退化为指数 ISS；\(0<\alpha<1\) 只能声称 Mittag-Leffler/代数衰减。

该 LMI 是充分非必要条件，且本身不是创新。\(\rho P\) 使联合最大化为双线性问题，实现时固定 \(\rho\) 二分求可行性。

### 推论 4.1（在线学习残差的小增益条件）

必须拆分

\[
r^\star(x)-r_\psi(\widehat x)
=[r^\star(x)-r_\psi(x)]
+[r_\psi(x)-r_\psi(\widehat x)].
\]

第一项可作为经验证的加性误差。若第二项满足

\[
\|B_\psi[r_\psi(x)-r_\psi(\widehat x)]\|_P
\le\ell_{\mathrm{eff}}\|e\|_P,
\qquad \ell_{\mathrm{eff}}<\rho,
\]

则把定理 4 中的 \(\rho\) 替换为

\[
\rho_{\mathrm{eff}}=\rho-\ell_{\mathrm{eff}}.
\]

测试集平均 RMSE 不能证明 \(\overline\varepsilon_\psi\) 或 \(\ell_{\mathrm{eff}}\)。若大模型没有区间/Lipschitz 证书，它只能作为经验预测基线，不能进入稳定性定理。

### 推论 4.2（相位估计误差）

若

\[
\dot e_\phi=-k_\phi e_\phi+d_\phi,
\]

则

\[
|e_\phi(s)|\le e^{-k_\phi s}|e_\phi(0)|+rac{\overline d_\phi}{k_\phi}.
\]

若相位通道对形态误差的增益为 \(b_\phi\ell_\phi\)，则

\[
\limsup\|e(s)\|_P
\le
\frac{\overline w_{0,P}
+b_\phi\ell_\phi\overline d_\phi/k_\phi}{\rho}.
\]

若 R 峰造成跳变/重置，还需证明跳映射在 \(P\) 范数下非扩张；上述连续定理不自动覆盖混杂跳变。

## 13. 有限记忆的显式精度律

对 \(0<\alpha<1\)，定义 GL 历史系数

\[
a_j(\alpha)=(-1)^{j+1}{\alpha\choose j}>0.
\]

截断到 \(L\) 个历史点后的尾质量为

\[
\tau_L(\alpha)
=\sum_{j=L+1}^\infty a_j(\alpha)
=\frac{\Gamma(L+1-\alpha)}
{\Gamma(1-\alpha)\Gamma(L+1)}
\sim\frac{L^{-\alpha}}{\Gamma(1-\alpha)}.
\]

令无量纲采样步长 \(h_s=h/\tau_0\)。若 \(\|z_k-z_0\|\le M_z\)，把截断项作为右端扰动时有

\[
\overline r_{\mathrm{mem}}
\le h_s^{-\alpha}M_z\tau_L(\alpha).
\]

所以记忆导致的最终观测半径满足

\[
R_\infty^{\mathrm{mem}}(L)
\le
\frac{b_mh_s^{-\alpha}M_z}{\rho}
\frac{\Gamma(L+1-\alpha)}
{\Gamma(1-\alpha)\Gamma(L+1)}
=O(L^{-\alpha}).
\]

该式直接揭示小 \(\alpha\) 下短记忆的风险。它是已有 GL 尾界与定理 4 的组合，不声称 GL 尾式本身是新结果。

若采用 positive-SOE/diffusive 近似并能证明

\[
\|K_\alpha-K_q\|_{L^1(0,T)}\le\epsilon_q,
\]

且驱动有界 \(\|f\|\le M_f\)，则

\[
\overline r_{\mathrm{mem}}\le M_f\epsilon_q,
\qquad
R_\infty^{\mathrm{mem}}\le b_mM_f\epsilon_q/\rho.
\]

仅假设 \(\dot z\) 有界，不能推出连续短记忆在无限时域的一致误差界。

## 14. 与 FO-EKF 的严格边界

定理 4 是确定性、固定阶、连续时间 observer 证书。它**不证明**：

- EKF 协方差最优或一致；
- Gaussian 白噪声下的均方稳定；
- NEES/NIS 校准；
- 变阶 \(\alpha(t)\)；
- 采样/保持和 R 峰跳变稳定；
- 非线性参数联合估计的全局收敛。

若保留 \(L\) 个历史状态，正确的概率实现需增广

\[
\xi_k=[x_k^\top,x_{k-1}^\top,\ldots,x_{k-L+1}^\top]^\top
\]

并传播完整协方差。当前状态预测协方差包含

\[
P^-_{00}=\sum_{i=1}^L\sum_{j=1}^L A_iP_{ij}A_j^\top+Q.
\]

忽略 \(i\ne j\) 的跨历史项一般既不保证保守，也不保证乐观。Sadeghian & Salarieh 2013 已给出随时间增长的完整线性 FKF，因此完整协方差不能单独作为创新。

首篇建议：

1. 用定理 4 的观测器作为有严格保证的安全主方法；
2. 把相关协方差一致的有限记忆 EKF/UKF 作为概率比较层；
3. 把截断偏差保留为确定性误差管，不随意 Gaussian 化；
4. 只有另立随机定理后才声称均方稳定或协方差一致。

## 15. 与定理逐项对应的验证协议

### 15.1 合成理论闭环（先做，失败即停）

**S0 解析恢复**

- \(\alpha\in\{0.5,0.7,0.9,0.97,1.0\}\)，多组 \(q,\lambda\)；
- 至少五个充分分离的心率、至少两个非零谐波；
- 用不同三元组恢复 \(\alpha\)，验证公共增益和公共相位旋转不变性；
- noise-free 恢复误差只允许数值舍入量级。

**S1 finite-terminal/full-memory 真值**

- 高精度 full-memory Caputo 求解器生成，不使用滤波器自身离散器；
- 多驻留时间、前一心率历史、不同初值；
- 逐点检查定理 1 的瞬态上界，不允许挑选平均意义通过。

**S2 反例与拒绝门**

- 独立 \(q_r\)、独立 \(\lambda_r\)、独立相位旋转、逐段幅值归一化；
- 心率间隔、SNR、窗长、固定延迟和相位误差扫描；
- 无效实验应被 \(\operatorname{Im}R\)、范围、条件数或多三元组一致性门拒绝。

**S3 有限噪声阶次界**

- 对每次试验保存实际 \(|\widehat\alpha-\alpha|\) 和理论上界；
- 所有满足假设的有界扰动试验均不得突破上界；
- 单独报告上界松弛比，不用平均 RMSE 掩盖越界。

**S4 observer 误差管**

- 有界测量噪声、模型扰动、相位误差、学习残差、有限记忆误差分别和组合注入；
- 多 seed、多初值、多 \(\alpha\)；
- 检查 \(\|e(s)\|_P\) 是否始终落在定理 4 的管内；
- \(\alpha=1\) 验证指数退化，\(\alpha<1\) 验证代数尾。

**S5 记忆—精度律**

- 扫描 \(L\) 或 SOE 模态数；
- 报告实测最坏误差、理论 \(O(L^{-\alpha})\) 包络、RAM 和延迟；
- 若所需 SOE 模态数与直接历史长度同量级，则计算创新失败。

### 15.2 Fantasia 多心率可行性预审

本地只读审计了 `../PUBLIC/fantasia` 的 40 位受试者 R 峰注释：

- 38/40 的 120 秒窗口平均心率跨度至少 5 bpm；
- 在窗口内 HR 变异系数 \(\le5\%\)、相邻三档平均 HR 至少相差 3 bpm 的条件下，15/40 具有三个候选窗口；
- 若相邻三档至少相差 5 bpm，仅 4/40 满足；
- 这些数字只说明“可做探索性筛查”，尚未证明形态 \(q\) 和阻尼 \(\lambda\) 跨窗口稳定。

因此 Fantasia 协议应为：

1. 受试者内选择至少三个、优选五个低变异心率窗口；
2. 同一 ECG 通道、同一物理增益、统一因果预处理；
3. R 峰相位锚定，禁止逐窗口幅值归一化；
4. 至少两个形态谐波、所有可用心率三元组；
5. 报告 \(\operatorname{Im}R\)、恢复 \(\lambda\) 的实正一致性、三元组/谐波阶次离散度；
6. 与允许 \(q_r\) 漂移的竞争模型比较；
7. 用留一心率预测而不是全数据后验拟合评估；
8. 受试者是统计单位，使用 subject-bootstrap 区间。

Fantasia 不能提供真实 \(\alpha\) 标签。它只能检验模型是否自洽、是否可证伪以及是否改善未来预测，不能证明 \(\alpha\) 是某个心肌生理常数。

### 15.3 PTB-XL 与 Challenge 的角色

- PTB-XL 每位记录通常只有 10 秒，不能承担三心率稳态辨识；用于冻结 \(\alpha\) 后的短时形态、噪声鲁棒性和五诊断超类分层。
- Challenge 2021 排除 PTB/PTB-XL 来源后用于锁模外部预测验证。
- 两者都不得通过同一患者的相邻窗口随机拆分夸大样本量。

### 15.4 数据/大模型的正确角色

大模型不用于“证明定理”。可使用容量匹配的 TCN/GRU/状态空间模型承担：

- 纯数据预测基线；
- 学习 residual 的候选，但只有独立误差和 Lipschitz 证书后才能进入 ISS 管；
- 检查分数阶结构是否在更高容量基线下仍带来留一心率预测收益。

所有模型使用相同的历史、相同的患者划分和调参预算。若大模型显著胜出而分数阶模型未提供更好的校准、效率或可解释的拒绝能力，应诚实把论文降为边界/负结果。

## 16. Go/No-Go 门

| 门 | 通过标准 | 失败动作 |
|---|---|---|
| T0 理论自洽 | 定理 1–4 的假设、支路、量纲、有限下限和证明无缺口 | 停止编码，修正理论 |
| T1 结构可辨识 | noise-free 多三元组恢复；反例明确失败；相位/增益 gauge 处理闭合 | 固定 \(\alpha\)，删除在线辨识 |
| T2 实际可辨识 | 预设噪声和形态漂移下的阶次区间有限，实际误差不突破理论界 | 仅保留模型拒绝结果 |
| T3 观测证书 | 存在共同 \(P,Y,\rho>0\)；所有有界扰动合成轨迹在误差管内 | 删除全局稳定主张或缩小运行域 |
| T4 真实数据自洽 | Fantasia 多三元组/多谐波一致性和留一心率预测同时通过 | 论文转为假设失效/负结果 |
| T5 方法收益 | 相对容量匹配整数阶和统计长记忆基线的主指标改善有受试者级置信区间 | 不声称分数阶优势 |
| T6 称谓 | 只有代理预测结论；没有解剖/内部状态真值就不升级 digital twin | 保持 ECG surrogate/digital shadow |

## 17. 目标期刊的现实排序

1. **ISA Transactions**：需完整有限驻留/漂移鲁棒理论和受控多心率验证；
2. **Biomedical Signal Processing and Control**：理论完整、真实数据主要为公开 ECG 时最现实；
3. **Physiological Measurement**：若重点是测量条件、模型拒绝和受试者内一致性；
4. **IEEE TBME**：只有取得受控起搏/运动平台稳定阶段数据、独立生理验证和更强多导联结果后再冲刺；
5. 当前理论不足以投 Automatica 或 IEEE TAC。

## 18. 接下来的执行顺序

1. 把定理 1–4 转成独立 LaTeX 推导，补齐引用和符号表；
2. 编写无数据依赖的解析恢复、单调性、扰动界和反例单元测试；
3. 用 full-memory Caputo 合成器完成 S0–S3；
4. 解 LMI 并完成 S4–S5；
5. 只有 T0–T3 通过后，才读取 Fantasia 波形做真实数据实验；
6. Fantasia 只作探索性验证；若“共同形态跨心率”失败，优先形成理论边界/负结果，不通过调参隐藏；
7. 若需要额外受控多心率数据库，先核实授权和体积；任何超过 100 MB 的下载必须再次取得用户批准。

## 19. 关键一手来源

- Tavazoei & Haeri, non-existence of periodic solutions, [Automatica 2009](https://doi.org/10.1016/j.automatica.2009.04.001)
- Yazdani & Salarieh, steady-state periodic response distinction, [Automatica 2011](https://doi.org/10.1016/j.automatica.2011.04.013)
- Zhang & Zhou, nonexistence and asymptotically periodic solutions, [CNSNS 2013](https://doi.org/10.1016/j.cnsns.2012.07.004)
- Takha, Talbi & Ravier, fractional McSharry model, [Medical Engineering & Physics 2024](https://doi.org/10.1016/j.medengphy.2024.104237)
- Guermah, Djennoune & Bettayeb, discrete fractional observability, [IJAMCS 2008](https://doi.org/10.2478/v10006-008-0019-6)
- Djennoune, Bettayeb & Al-Saggaf, output-memory observability, [IJAMCS 2019](https://doi.org/10.2478/amcs-2019-0014)
- Nazarian, Haeri & Tavazoei, frequency-domain identifiability limitations, [ISA Transactions 2010](https://doi.org/10.1016/j.isatra.2009.11.007)
- Alavi et al., structural identifiability of fractional models, [arXiv 1511.01402](https://doi.org/10.48550/arXiv.1511.01402)
- Sadeghian et al., general linear fractional Kalman filter, [Mechatronics 2013](https://doi.org/10.1016/j.mechatronics.2013.02.006)
- Rauh & Lahme, finite-memory observer error correction, [Annual Reviews in Control 2025](https://doi.org/10.1016/j.arcontrol.2025.101018)
- Naifar & Ben Makhlouf, polytopic fractional observer, [Mathematical Problems in Engineering 2021](https://doi.org/10.1155/2021/6699756)
- Reif et al., discrete-time EKF stochastic stability, [IEEE TAC 1999](https://doi.org/10.1109/9.754809)

## 20. 最终可行性判断

**理论上可行，但只在严格条件下成立。**

最强的新意不是再造一个 FO-EKF，而是把真实 ECG 中通常被当作 nuisance 的心率变化转化为可检验激励，并构造一个对相位、增益和形态尺度规范不敏感的阶次证书；随后让所有近似误差进入可计算的观测误差管。该理论具有明确反例和拒绝条件，因而可被数据证伪。

最大风险是不同心率下 ECG 形态和阻尼并不保持不变。现有 Fantasia 足以做探索性验证，但不足以单独支撑强生理结论。若真实数据拒绝共同形态假设，仍可形成一篇有价值的“分数阶 ECG 阶次不可辨识边界”论文；不能退回到逐搏优化 MSE 来掩盖理论失败。

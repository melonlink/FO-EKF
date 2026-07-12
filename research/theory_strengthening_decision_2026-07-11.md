# FO-ECG 论文理论强化决策：联合标量谐波可行集与证书驱动的心率设计

- **决策日期**：2026-07-11
- **工作分支**：`dev-codex`
- **状态**：理论优先；暂停新的真实数据与大模型验证
- **结论**：继续完善当前理论，只加入与同一辨识问题不可分割的联合集合证书和心率/驻留设计；不叠加 FO-EKF、UKF、神经网络或另一套 observer 作为首篇主方法

## 1. 最终选择

首篇论文只回答一个问题：

> 在有限 ECG 记录、公共相位/增益 gauge、心率相关形态漂移、有限驻留和有界测量误差同时存在时，分数阶阶次最多能被确定到多窄的集合；应选择哪些心率窗口、顺序和驻留时间，才能最小化最坏不可区分阶次宽度？

采用一个统一的复合方法，而不是多个算法拼接：

1. **声明的标量谐波代理的联合 set-membership 证书**：同时使用全部心率和全部谐波，强制公共实正阻尼、公共阶次以及每个谐波跨心率共享的形态参数或其有界 restitution 漂移；
2. **certificate-aware 心率/驻留设计**：直接最小化上述联合可行阶次集合的最坏直径，而不是再使用普通 Fisher 信息、局部 MSE 或通用最优输入设计。

observer、FO-EKF、有限记忆实现和学习残差移到第二篇。原因不是它们没有价值，而是未知阶次 observer 会新增算子失配、完整历史、采样、切换和 R 峰重置稳定性问题，无法由当前的固定阶次 Mittag-Leffler ISS 定理自动覆盖。

## 2. 文献查重后的创新边界

以下宽泛主张均不可再作为创新：

- 分数阶模型的有界频域 set-membership 辨识：Malti et al. 2010 和 Khemane et al. 2012 已直接覆盖；
- 有界噪声下同时估计分数阶阶次和系数：Zhang et al. 2021 已覆盖；
- 为分数阶参数/阶次选择最优频率：Abrashov et al. 2017 和 Malti et al. 2022 已覆盖；
- 用 LMI/Fisher 信息设计分数阶辨识输入谱：Jakowluk and Swiercz 2025 已覆盖。

仍可守的窄创新必须是下列结构的联合结果：

1. ECG 公共相位和公共复增益商空间下的物理可辨识性；
2. 有限下限 Caputo 历史、有限驻留和 ECG restitution 漂移共同进入确定性数据圆盘；
3. 对全部心率、全部谐波和公共实正阻尼的**联合标量代理**可行阶次集合；
4. 以最坏可行集合直径为目标的被动心率窗口/驻留选择；未来有受控运动或起搏数据后，才称主动实验设计。

真实 ECG 的形态随心率变化且存在 QT/RR hysteresis，因此跨心率固定 $q,\lambda$ 只能是零漂移特例，不能作为一般生理事实。

针对本轮两心率公式做了按模型结构、公式片段和引用网络的定向检索：尚未发现“未知公共复增益 + 正实阻尼 + 两个完整复响应 + 严格最小频点数”的直接同构定理；但分数阶模型的两/多频辨识、Nyquist 余切几何和多正弦最小二乘均有相邻先例。因此当前只称内部新引理，正式系统综述完成前不使用 `first`。

## 3. 对现有 S0--S3 结论的必要修正

现有 S0--S3 仍证明了以下内容：

- 理想三心率差商公式、$F'(\alpha)>0$ 和 $F''(\alpha)>0$ 的实现正确；
- 有限下限 Caputo 驻留上界在所测合成案例中无越界；
- $Z\to W\to R\to\alpha$ 的链式扰动界在声明的合成扰动模型下有覆盖；
- G0 的 `NOT_CERTIFIABLE` 语义确实不可由纯数据检验替代。

但它们没有关闭完整理论 T0，原因是：

1. 当前所谓“精确复圆盘可行集”只精确描述三点代数关系
   \[
   W_3-W_1=F(\alpha)(W_2-W_1),
   \]
   没有强制恢复的 $\lambda$ 为正实数，也没有强制不同三元组和谐波共享同一组代理参数；因此它只是声明标量代理完整集合的外松弛。
2. 固定 4097 网格不是完备集合算法，可能漏掉窄分量或切触根。
3. 逆响应圆盘包含零或两圆盘重叠，只表示倒数/差商精度无法认证，不等价于声明标量代理无解。
4. 当前 $\epsilon_Z$ 是人工给定的合成半径，尚未从有限 ECG 窗、HRV、采样、因果滤波延迟和测量噪声推出。
5. 三心率不是已证明的最小实验；利用 $\lambda>0$ 的物理相位约束，两心率已经可以在理想模型下全局辨识阶次。

所以当前状态应写成：**S0--S3 在声明的合成模型下通过；标量代理的完整 T0 尚未通过。**

## 4. 新的理论核心 I：两心率全局物理可辨识

考虑同一非零谐波在两个不同心率下的理想复响应

\[
Z_r=\frac{q}{\lambda+(i\nu_r)^\alpha},
\qquad
0<\nu_1<\nu_2,
\quad \alpha\in I=[\alpha_-,1],\quad \alpha_->0,
\quad \lambda>0,
\quad q\ne0.
\]

取主值支路

\[
(i\nu)^\alpha=\nu^\alpha e^{i\pi\alpha/2}.
\]

令

\[
W_r=\frac1{Z_r},\qquad
D=W_2-W_1,\qquad
C=\frac{W_1}{D},\qquad
r=\frac{\nu_2}{\nu_1}>1.
\]

### 定理 A（两心率唯一辨识，待系统查重后编号）

理想模型产生的数据满足 $D\ne0$ 和 $\operatorname{Im}C<0$。真实阶次 $\alpha$ 是方程

\[
H_C(\beta)=1,
\]

在 $0<\beta\le1$ 内的唯一解，其中

\[
H_C(\beta)
=(r^\beta-1)
\left[
\operatorname{Re}C+
\operatorname{Im}C\cot\left(\frac{\pi\beta}{2}\right)
\right].
\]

得到 $\alpha$ 后，测量单位中的有效形态系数 $q$ 和 $\lambda$ 唯一恢复为

\[
q_\alpha=
\frac{e^{i\pi\alpha/2}(\nu_2^\alpha-\nu_1^\alpha)}{D},
\]

\[
\lambda_\alpha=
\nu_1^\alpha e^{i\pi\alpha/2}
\left[C(r^\alpha-1)-1\right]
\in\mathbb R_{>0}.
\]

### 证明

理想数据给出

\[
C=
\frac{\lambda+e^{i\pi\alpha/2}\nu_1^\alpha}
{e^{i\pi\alpha/2}(\nu_2^\alpha-\nu_1^\alpha)}
=
\frac{1+(\lambda/\nu_1^\alpha)e^{-i\pi\alpha/2}}
{r^\alpha-1},
\]

故 $\operatorname{Im}C<0$。对任意候选 $\beta$，由两点差分只能得到上述唯一候选 $q_\beta$，再代回第一点得到 $\lambda_\beta$。令

\[
t=r^\beta-1,\qquad
\theta=\frac{\pi\beta}{2},qquad
C=u+iv.
\]

则

\[
\operatorname{Im}\lambda_\beta=0
\Longleftrightarrow
(tu-1)\sin\theta+tv\cos\theta=0
\Longleftrightarrow
t(u+v\cot\theta)=1.
\]

即 $H_C(\beta)=1$。为直接证明全局唯一性，令

\[
A=\frac{\pi\alpha}{2},\qquad
B=\frac{\pi\beta}{2},\qquad
\kappa=\frac{\lambda}{\nu_1^\alpha}>0.
\]

把理想数据产生的 $C$ 代入候选方程，并使用恒等式

\[
\cos A-\sin A\cot B
=\frac{\sin(B-A)}{\sin B},
\]

可得

\[
r^\beta-r^\alpha
+\kappa(r^\beta-1)
\frac{\sin(B-A)}{\sin B}
=0.
\]

若 $\beta>\alpha$，上式两项均严格为正；若 $\beta<\alpha$，两项均严格为负。因此只有 $\beta=\alpha$ 能成立，真实阶次是唯一根。该符号论证直接覆盖 $\alpha=1$ 或 $\beta=1$ 的端点。

在根处又有 $tu-1=-tv\cot\theta$，于是

\[
\operatorname{Re}\lambda_\beta
=
-\nu_1^\beta\frac{tv}{\sin\theta}>0.
\]

所以唯一实数候选自动满足正阻尼条件，且 $D\ne0$ 保证 $q_\beta\ne0$。此外，

\[
\lim_{\beta\to0^+}H_C(\beta)
=\frac{2\operatorname{Im}C\log r}{\pi}<0.
\]

$\alpha=0$ 时频率项退化为常数且 $D=0$，所以必须排除；当 $\alpha_-\to0$ 或 $r\to1$ 时，虽然每个固定参数仍结构唯一，但不存在统一的鲁棒辨识裕度。后续设计必须预注册 $\alpha_->0$ 和最小心率分离。证毕。

### 直接推论

- 一个心率下，对任意候选 $\alpha,\lambda>0$ 都可选择 $q=Z[\lambda+(i\nu)^\alpha]$，故阶次不可辨识；在“两个已知不同正频率、可测完整复响应、跨率公共 $q$、公共 $\lambda\in\mathbb R_{>0}$、固定主值支路且无逐段相位/延迟”的声明模型内，两个心率是最小实验数。
- $C$ 对全部 $Z_r\mapsto gZ_r$、$g\in\mathbb C\setminus\{0\}$ 不变，因此公共导联增益和公共相位旋转被消除；$\alpha,\lambda$ 不变，而恢复的 $q$ 是测量单位中的有效系数 $gq_{\rm physical}$。
- 三心率差商不再是“最小辨识定理”，而是无需单独求解 $\lambda$ 的闭式过定约束，可用于模型拒绝、提高鲁棒性和检验跨心率参数共享。
- 只有幅值、逐段独立相位、逐段独立增益或自由 $q_r$ 时，上述定理失效，必须返回 `NOT_CERTIFIABLE`。

对非理想数据，$\operatorname{Im}C<0$ 只是必要条件而不是充分条件。无根、仅数值边界根、候选阻尼非正或根区间过宽必须分别进入模型拒绝或不可认证语义；小 $\beta$ 的实现需使用 `expm1` 和稳定的余切展开。

该定理目前是新的内部推导；正式论文中在完成按公式和引用网络的系统查重前，不使用 `first`。

## 5. 新的理论核心 II：声明标量谐波代理的联合前向圆盘集合

不再优先对 $Z$ 求倒数。设有限记录给出

\[
|\widetilde Z_{r,m}-Z_{r,m}|\le\epsilon_{r,m}.
\]

对每个预注册的入选谐波，还要求测量单位中的有效形态系数满足独立校准的紧幅值先验

\[
0<q_{\min,m}\le |q_m|\le q_{\max,m}<\infty,
\qquad
\mathcal E_m
=\{q\in\mathbb C:q_{\min,m}\le|q|\le q_{\max,m}\}.
\]

该闭环带同时排除零谐波退化并保证参数集紧致；上下界都不能由同一批辨识残差反向调节。

对候选 $\beta\in I$、$\ell\in\Lambda=[\lambda_-,\lambda_+]\subset\mathbb R_{>0}$，定义

\[
d_{r,m}(\beta,\ell)
=\ell+(i\nu_{r,m})^\beta
\]

和复平面圆盘

\[
\mathcal Q_{r,m}(\beta,\ell)
=
\mathcal D\left(
d_{r,m}\widetilde Z_{r,m},
|d_{r,m}|\epsilon_{r,m}
\right).
\]

### 定理 B（无漂移联合代理可行集）

固定 $\beta,\ell$ 后，存在一个跨全部心率共享的非零形态系数 $q_m$，使

\[
\left|
\widetilde Z_{r,m}-
\frac{q_m}{d_{r,m}(\beta,\ell)}
\right|
\le\epsilon_{r,m},
\qquad \forall r,
\]

当且仅当

\[
\left(\bigcap_{r=1}^{K}
\mathcal Q_{r,m}(\beta,\ell)\right)
\cap\mathcal E_m
\ne\varnothing.
\]

因此完整的联合阶次集合恰为

\[
\mathcal A_{\rm phys}
=
\operatorname{proj}_{\beta}
\left{
(\beta,\ell)\in I\times\Lambda:
\left(\bigcap_r\mathcal Q_{r,m}(\beta,\ell)\right)
\cap\mathcal E_m
\ne\varnothing,
\ \forall m\in\mathcal M
\right}.
\]

**证明**：原测量约束与

\[
|d_{r,m}\widetilde Z_{r,m}-q_m|
\le |d_{r,m}|\epsilon_{r,m}
\]

等价；右式恰表示 $q_m\in\mathcal Q_{r,m}$。对所有 $r$ 同时成立且满足幅值先验，等价于圆盘交与紧环带 $\mathcal E_m$ 相交。对所有谐波同时成立时，$\beta,\ell$ 共享，而每个 $q_m$ 是对应闭交集中的一点。证毕。

该构造同时解决四个问题：

1. $\ell$ 从定义上就是公共实正阻尼；
2. 所有心率共享同一个 $q_m$，不是逐三元组事后求交；
3. 所有谐波共享同一 $\beta,\ell$；
4. 不做倒数，所以观测圆盘包含零时仍有合法的物理可行性语义。

固定 $\beta,\ell$ 时，普通圆盘交是二维凸可行问题，可先用 SOCP 判断；再与环带 $\mathcal E_m$ 相交需要二维精确圆几何或有保证的最小/最大模计算，不能由一次普通 SOCP 代替。由于 $I\times\Lambda$ 以及全部参数/误差集合均取紧集，完整投影仍必须使用带区间包络的分支定界；固定网格只能做近似可视化，不能称 complete/certified。

### 有界形态漂移的精确扩展

若

\[
q_{r,m}=\bar q_m+\delta q_{r,m},
\qquad
|\delta q_{r,m}|\le\rho^q_{r,m},
\qquad q_{\min,m}\le|\bar q_m|\le q_{\max,m},
\]

则只需把圆盘半径精确扩大为

\[
|d_{r,m}|\epsilon_{r,m}+\rho^q_{r,m}.
\]

这是两个复圆盘的 Minkowski 和，不是启发式误差相加，但“精确”严格限于每个 $(r,m)$ 的 $\delta q_{r,m}$ 是相互独立的复圆盘变量、且幅值下界只施加在名义 $\bar q_m$ 上。若 restitution 是跨谐波共享增益/相位、仅实幅值漂移、相关不确定集，或还要求每个 $q_{r,m}$ 非零，则必须保留共享变量和原约束；简单扩圆盘只能作为外包络。

公共未知增益 $g$ 下，$\widetilde Z,\epsilon,\rho^q,q_{\min},q_{\max}$ 必须分别按 $g,|g|,|g|,|g|,|g|$ 在同一测量单位中缩放。未独立标定该尺度时只能恢复有效形态系数，不能解释物理 $q_m$。

### 有界阻尼漂移的精确扩展

若

\[
\lambda_r=\bar\lambda+\delta\lambda_r,
\qquad
|\delta\lambda_r|\le\rho^\lambda_r,
\]

取紧先验 $\bar\lambda\in\bar\Lambda\subset\mathbb R_{>0}$。把 $\delta\lambda_r$ 保留为每个心率的实变量，并要求它在全部谐波中共享。联合集合为

\[
\operatorname{proj}_{\beta}
\left\{
(\beta,\bar\lambda,\delta\lambda_1,\ldots,\delta\lambda_K)
\in I\times\bar\Lambda\times\mathbb R^K:
\begin{array}{l}
|\delta\lambda_r|\le\rho^\lambda_r,\\
\bar\lambda+\delta\lambda_r\in\Lambda_r\subset\mathbb R_{>0},
\quad \forall r,\\
\left(\bigcap_r\mathcal Q_{r,m}
(\beta,\bar\lambda+\delta\lambda_r)\right)
\cap\mathcal E_m
\ne\varnothing,
\ \forall m
\end{array}
\right}.
\]

若需要把 $\bar\lambda$ 解释成唯一的名义阻尼，还必须增加锚定规范（例如 $\sum_r w_r\delta\lambda_r=0$ 或固定一个 $\delta\lambda_r=0$）；否则 $\bar\lambda$ 与漂移存在平移冗余，但这不妨碍对 $\beta$ 做集合投影。

当 $q$ 与 $\lambda$ 漂移同时存在时，必须在同一联合投影中令

\[
d_{r,m}
=\bar\lambda+\delta\lambda_r+(i\nu_{r,m})^\beta,
\qquad
\operatorname{rad}\mathcal Q_{r,m}
=|d_{r,m}|\epsilon_{r,m}+\rho^q_{r,m},
\]

并保留跨谐波共享的同一个 $\delta\lambda_r$。为加速计算可以把阻尼漂移外包成更大圆盘，但必须明确标成 outer enclosure；不能再称声明代理的精确集。

本节得到的是**所声明标量 harmonic surrogate 的联合集合**。它尚未施加真实波形的共轭对称、跨谐波生理耦合或相关估计误差；若这些被列入物理模型，当前集合相对于更强模型仍是外松弛。

## 6. 必须补齐的数据入口定理

上述定理只有在 $\epsilon_{r,m}$ 有协议级覆盖保证时才能用于真实 ECG。下一步必须对明确的相位域解调器或有限窗复 Fourier/最小二乘估计器证明

\[
Z_{r,m}\in
\mathcal D(\widetilde Z_{r,m},\epsilon_{r,m}),
\]

其中

\[
\epsilon_{r,m}
\le
\epsilon^{\rm dwell}_{r,m}
+\epsilon^{\rm history}_{r,m}
+\epsilon^{\rm HRV}_{r,m}
+\epsilon^{\rm window}_{r,m}
+\epsilon^{\rm sample}_{r,m}
+\epsilon^{\rm delay}_{r,m}
+\epsilon^{\rm noise}_{r,m}.
\]

要求如下：

- 窗口跨整数个 R--R 周期，或显式保留非整数窗泄漏；
- 使用统一 R 峰相位锚和因果预处理；
- HRV 由相位速度偏差上界进入，而不是用“窗口平均心率”替代；
- Caputo 前史和切换瞬态使用 Mittag-Leffler 核尾进入平均复系数误差；
- 采样和数值积分误差有确定性包络；
- 测量噪声采用与定理一致的有界语义；随机置信区间不能直接冒充确定性圆盘。

在该桥定理完成以前，不开始 Fantasia 波形上的阶次更新，也不把人工相对误差半径当作真实覆盖证书。

## 7. 理论核心 III：最坏不可区分宽度与心率/驻留设计

令设计

\[
\mathcal D=
(\{\nu_r\},\{T_r\},\sigma,\mathcal M)
\]

包含心率、驻留时间、访问顺序和谐波集合。定义允许不确定性下阶次 $\alpha$ 的全部可能观测集合 $\mathcal Y_{\mathcal D}(\alpha)$，以及一致性集合

\[
\mathcal A_{\mathcal D}(y)
=\{\alpha:y\in\mathcal Y_{\mathcal D}(\alpha)\}.
\]

最坏不可区分宽度定义为

\[
\Delta(\mathcal D)
=
\sup\left\{
|\alpha-\alpha'|:
\mathcal Y_{\mathcal D}(\alpha)
\cap
\mathcal Y_{\mathcal D}(\alpha')
\ne\varnothing
\right\}.
\]

需证明

\[
\sup_y\operatorname{diam}\mathcal A_{\mathcal D}(y)
=\Delta(\mathcal D),
\]

并在生理心率区间、最小心率间隔、总时间预算和谐波质量约束下求

\[
\mathcal D^\star
\in\arg\min_{\mathcal D\in\mathfrak D}
\Delta(\mathcal D).
\]

首篇至少要完成：

1. 最优设计存在性；
2. 添加有效观测只能使 $\Delta$ 非增，但未必严格改善；
3. 两心率物理证书和三心率差商证书在 close-rate 极限下的不同条件数阶；
4. 有界 $q,\lambda$ restitution 和 dwell/history 预算如何改变最优窗口、顺序和驻留；
5. 一个有全局上下界或可证明近似比的求解方法，而不是局部黑箱优化。

Fantasia 中心率不能主动施加，因此只称 **certificate-aware retrospective rate-window selection**。未来受控运动、药理或起搏实验才使用 optimal experiment design。

## 8. 是否加入低阶 restitution 新模型

并行评估了一个更强但也更冒险的扩展。令 $x=\log\nu$，在逆频响域使用

\[
W_m(x)=A_m(x)+e^{\alpha x}B_m(x)+r_m(x),
\qquad
A_m,B_m\in\mathcal P_p,
\qquad |r_m(x_k)|\le\varepsilon_{m,k}.
\]

它把随心率变化的 nuisance 限制在固定维数空间，不再允许每个心率任意漂移。扩展 Chebyshev 系给出的正确门槛是：

| 目标 | 最低不同心率数 |
|---|---:|
| $\alpha$ 已知时恢复 $A_m,B_m$ | $2p+2$ |
| 正则局部联合辨识 | $2p+3$ |
| 无附加系数约束下的最坏全局唯一 | $3p+3$ |

原先猜测“$2p+3$ 足以全局唯一”是错误的；$p=1$ 时五个心率存在两个不同阶次产生完全相同五点数据的精确反例，六个心率才给无约束最坏全局保证。$p=0$ 时局部和全局门槛都退化为三心率差商。

该路线的优点是有尖锐样本复杂度、非正则分层和“模型余项偏差--阶次可辨识性”权衡，理论上明显强于任意有界漂移。缺点同样实质：

1. 独立多项式 $A_m,B_m$ 放松了原物理模型中的 $A_m/B_m=\lambda e^{-i\pi\alpha/2}$ 耦合；
2. $p=1$ 的全局结果至少需要六个充分分离且有驻留保证的心率，当前 Fantasia 预审很可能不足；
3. 若真实 $B_m$ 的次数低于声明的 $p$，阶次落在 Jacobian 非正则分层，噪声会被 nuisance 多项式吸收；
4. ECT 零点计数和指数多项式恢复本身是已有数学工具，创新仍须来自 ECG 物理推导、余项证书和最坏歧义设计。

**决策**：暂不把该模型并入第一版主定理链。先关闭保留完整物理耦合的 R0--R4；只有独立校准能支持低阶 restitution 余项界、且数据协议能提供至少六个有效心率时，再把 $p=1$ 作为增强定理或第二篇的一般化。这样既保留新方法方向，也不为追求公式数量牺牲物理可验证性。

相关数学背景包括指数多项式的 Chebyshev 系和广义 Prony 恢复；正式采用时需单独完成定向查重。

## 9. 第一篇论文的最小定理链

1. **P1 商空间与协议命题**：公共相位/增益 gauge、可辨识对象、因果预处理边界；
2. **T1 有限记录到复圆盘包络**：dwell、history、HRV、window、sampling、delay、noise 全部进入 $\epsilon_Z$；
3. **T2 最小结构可辨识性**：一心率不可辨、两心率全局唯一、三/多心率为过定模型检验；
4. **P2 不可能性命题**：自由 $q_r$、独立相位/增益、只有幅值、无漂移上界时不可辨识；
5. **T3 联合标量代理可行集**：全部 $K\times M$ 约束、公共 $\lambda>0$、有界 restitution、完备投影算法；
6. **C3 覆盖和拒绝语义**：空集为模型/误差预算不相容，非空宽集为 `NOT_CERTIFIABLE`，窄集才为 `PASS`；
7. **T4 实际可辨识性与设计**：最坏不可区分宽度、条件数渐近律、心率/驻留/顺序设计。

现有三心率闭式公式降为 T2 的推论；有限驻留 ML 界降为 T1 的误差项；不再把多个标准模块并列为独立创新。

## 10. 执行门与停止条件

在继续任何大规模数值或真实数据工作前，依次关闭：

| 门 | 通过条件 | 失败动作 |
|---|---|---|
| R0 两心率定理 | 完成严格证明、符号审计、边界/反例与定向查重 | 回退到三心率，但删除“最小”措辞 |
| R1 联合集合定理 | 前向圆盘等价、共享参数、漂移语义和完备算法证明闭合 | 只保留结构辨识，不声称 robust exact set |
| R2 数据到圆盘 | 每个误差项有可观测/协议上界和覆盖证明 | 禁止进入真实阶次认证 |
| R3 minimax 设计 | 存在性、下界、条件数律和有保证算法闭合 | 设计仅作经验窗口筛选，不列理论贡献 |
| R4 系统查重 | 按公式、标题和引用网络未发现同构定理 | 使用窄措辞，不声称 `first` |

只有 R0--R4 通过后，才重新设计合成验证；只有 R2 通过后，才读取 Fantasia 波形做探索性联合证书。任何额外数据或软件预计超过 100 MB，必须先取得用户批准。

## 11. 期刊与成果判断

- **ISA Transactions**：R0--R4 完整、联合证书有真实 ECG case study 后有现实可行性；
- **Automatica**：还需把结果抽象为一般“未知公共复 nuisance 下的有限记录多实验分数阶辨识”，ECG 只作实例，并给出必要充分联合可辨识和 minimax 设计；
- **IEEE TAC**：当前标量路线仍偏窄，通常还需 MIMO/多状态或更一般的信息下界，暂不作为首篇目标。

建议题目：

> **Gauge-Invariant Finite-Record Set Identification of Fractional ECG Dynamics under Rate-Dependent Restitution**

若 T4 足够强，可用：

> **Joint Set Identifiability and Certificate-Aware Rate-Window Selection for Fractional ECG Dynamics**

最终判定为 **CONDITIONAL GO**。相比继续加入新滤波器，先把“两心率最小可辨识 + 联合前向圆盘标量代理集合 + finite-record 误差入口 + 最坏集合直径设计”证明完整，理论强度和创新可守性都更高。

## 12. 本轮新增的一手文献边界

- Malti et al., bounded frequency-domain set membership, [CNSNS 2010](https://doi.org/10.1016/j.cnsns.2009.05.005)
- Khemane et al., robust fractional set-membership estimation, [Signal Processing 2012](https://doi.org/10.1016/j.sigpro.2011.12.008)
- Nazarian et al., frequency-content identifiability limits, [ISA Transactions 2010](https://doi.org/10.1016/j.isatra.2009.11.007)
- Nazarian and Haeri, order-distribution frequency-domain identification, [Signal Processing 2010](https://doi.org/10.1016/j.sigpro.2010.02.008)
- Zhang et al., operational-matrix set membership, [Journal of the Franklin Institute 2021](https://doi.org/10.1016/j.jfranklin.2021.10.020)
- Abrashov et al., robust fractional experiment design, [IEEE TAC 2017](https://doi.org/10.1109/TAC.2016.2614910)
- Malti et al., one/multiple optimal frequencies for fractional models, [CNSNS 2022](https://doi.org/10.1016/j.cnsns.2022.106337)
- Jakowluk and Swiercz, LMI-based optimal input spectrum, [Applied Sciences 2025](https://doi.org/10.3390/app152312665)
- Ramirez et al., T-wave morphology restitution, [JAHA 2017](https://doi.org/10.1161/JAHA.116.005310)
- Andrsova et al., rate-change and QT/RR hysteresis, [Frontiers in Physiology 2022](https://doi.org/10.3389/fphys.2021.814542)

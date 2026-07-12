# FO-EKF R0–R4 理论阶段执行报告

日期：2026-07-11

分支：`dev-codex`
阶段边界：只做理论闭环、文献边界和小规模合成算法验证；不进入真实 ECG 阶次结论。

## 1. 总结论

本方向值得继续，但当前只能给出：

> **理论核心有条件可行；真实数据验证暂不放行。**

两心率完整复响应在“公共复形态/导联增益 + 公共正实阻尼 + 固定分数阶支路”的声明模型下，可以严格、全局、唯一地恢复分数阶阶次；一心率严格不足。该结果已经提升为公共复增益商空间上的全局单射/微分同胚结构，并给出了边界和病态性。

真正阻止当前论文进入真实 ECG 的不是优化器，也不是深度模型，而是 R2：现有数据尚不能给出具有确定性覆盖保证的复响应圆盘。R 峰时序误差、有限 Caputo 前史、谐波尾、抗混叠/采样误差、滤波相对群延迟、确定性噪声和标量代理模型余项都缺少独立上界。

因此下一步应继续补强理论与实验协议，不应立即叠加新模型或大模型。机器学习后续只能用于验证或提高中心估计精度，不能替代覆盖性误差界。

## 2. R0–R4 门控

| 门 | 判定 | 已完成 | 尚未闭合 |
|---|---|---|---|
| R0 最小全局辨识 | **PASS** | 一率不可辨；两率必要且充分；公共复 gauge 完全不变量；全局唯一性；端点与条件数 | 仅在预注册假设成立时有效 |
| R1 联合集合反演 | **解析理论 PASS；机器实现 OPEN** | 固定参数圆盘交-环带等价；独立复形态漂移；共享实阻尼漂移；多率多谐波联合盒原型 | 定向舍入、完整投影、rate-ID 对齐、阻尼漂移 gauge |
| R2 有限记录到圆盘 | **条件定理 PASS；协议 FAIL** | WLS 有限样本界；Caputo history/dwell；HRV 核传播；泄漏、采样、延迟、噪声；最小驻留/SNR 推论 | 缺少所有必要的独立数值上界 |
| R3 最坏歧义设计 | **THEOREM PASS** | 最小可靠集合；最坏直径恒等式；圆盘分离阈值 2；最优设计存在性；近率条件数 | 全局认证数值优化器 |
| R4 创新边界 | **NARROW PASS** | 完成最接近公式与方法的碰撞审计 | 不能把通用分数阶辨识、圆盘集合或实验设计单独宣称为创新 |

阶段总门：

- 理论论文路线：**CONDITIONAL GO**
- 当前真实 ECG 阶次认证：**NO-GO / NOT_CERTIFIABLE**
- 深度学习或大模型验证：**暂缓**

## 3. R0：理论核心

固定参考时间 \(\tau_\star\)，令频率无量纲化。对一个正谐波：

\[
Z_j=\frac{q}{\lambda+(i\nu_j)^\alpha},\qquad
0<\nu_1<\nu_2,
\]

其中 \(q\neq0\) 是未知公共复系数，\(\lambda>0\)，
\(0<\alpha\le1\)。

令

\[
W_j=Z_j^{-1},\quad D=W_2-W_1,\quad
C=\frac{W_1}{D},\quad r=\frac{\nu_2}{\nu_1}>1.
\]

公共复增益 \(Z_j\mapsto gZ_j\) 被 \(C\) 完全消去。进一步令
\(\kappa=\lambda/\nu_1^\alpha\)，得到商空间映射

\[
(\alpha,\kappa)\mapsto
C=\frac{1+\kappa e^{-i\pi\alpha/2}}{r^\alpha-1}.
\]

其像为

\[
\left\{u+iv:u\ge\frac1{r-1},\ v<0\right\},
\]

在 \(0<\alpha<1,\kappa>0\) 的内部是微分同胚，Jacobian 严格为正：

\[
\det\frac{\partial(\Re C,\Im C)}
{\partial(\alpha,\kappa)}
=
\frac{r^\alpha\log r\sin(\pi\alpha/2)}
{(r^\alpha-1)^3}
+\frac{\kappa\pi}{2(r^\alpha-1)^2}>0.
\]

阶次是下式在 \((0,1]\) 中的唯一根：

\[
H_C(\beta)
=(r^\beta-1)
\left[
\Re C+\Im C\cot\left(\frac{\pi\beta}{2}\right)
\right]=1.
\]

关键边界：

- \(\alpha=1\) 是合法模型边界，不能因为根在边界就拒绝；
- \(\alpha\to0^+\) 或 \(r\to1^+\) 时，结构唯一仍可存在，但不存在统一抗噪裕度；
- \(H_C\) 对任意非理想数据不保证全区间单调，圆盘数据必须保留全部可行分支；
- 两率只完成恰定辨识，至少三率才能形成独立的精确模型拒绝自由度。

## 4. R1：联合集合算法及审计修正

固定 \((\alpha,\bar\lambda,\delta\lambda_r)\) 时，

\[
|\widetilde Z_{r,m}-Z_{r,m}|\le\epsilon_{r,m}
\]

与公共形态系数 \(q_m\) 落入复圆盘交严格等价：

\[
q_m\in
\bigcap_r
\mathcal D\left(
d_{r,m}\widetilde Z_{r,m},
|d_{r,m}|\epsilon_{r,m}+\rho^q_{r,m}
\right).
\]

代码原型已经保留了：

- 各率共享同一 \(q_m\)；
- 所有谐波共享 \(\alpha,\bar\lambda\)；
- \(\delta\lambda_r\) 是跨谐波共享的实变量，未错误膨胀成独立复误差；
- outer-union 圆盘、uniform inner core 和有限预算分支定界；
- `FEASIBLE / INFEASIBLE / UNKNOWN` 三值语义。

第二轮审计发现，原浮点容差曾允许一个低于环带下界
\(0.5\times10^{-10}\) 的点被提升为 `INNER_INCLUDED`。本阶段已：

1. 将最终可行见证改为对未膨胀约束复核；
2. 对 inner inclusion 增加严格向内余量；
3. outer/inner 冲突时强制降级为 `UNKNOWN`；
4. 增加该反例和声明参数域回归测试；
5. 明确代码只是浮点原型，不能充当论文级机器证书。

outer exclusion 仍没有定向舍入误差证明，所以 R1 不能写成“已完成认证求解器”。

## 5. R2：为什么现在不能跑真实数据结论

已建立条件性半径分解：

\[
\epsilon_{r,m}
=\epsilon^{\rm hist}
+\epsilon^{\rm dwell}
+\epsilon^{\rm HRV}
+\epsilon^{\rm omitted}
+\epsilon^{\rm sample}
+\epsilon^{\rm delay}
+\epsilon^{\rm noise}
+\epsilon^{\rm model}
+\epsilon^{\rm filt-init}.
\]

理论已覆盖：

- Mittag–Leffler 核下的有限记录前史与切换驻留；
- HRV/相位变形先在输入侧成界，再通过分数阶核传播；
- 非整数窗和未建模谐波泄漏；
- 时间戳、ADC、插值及梯形求积误差；
- 频响幅值与相对群延迟；
- 确定性能量噪声；
- 最小证书驻留和最小确定性 SNR。

但是以下数值必须来自协议、硬件校准或独立数据，不能由同一拟合残差反推：

1. R 峰绝对时序误差和潜在搏内相位误差；
2. 抗混叠/带宽或导数上界；
3. 确定性噪声能量或幅值界；
4. 前端真实频响、率间相对群延迟和滤波初始化尾；
5. Caputo 记录前历史、前段输入和状态幅值界；
6. 未建模谐波尾的绝对可和界；
7. 表面 ECG 相对所声明标量代理的模型余项界。

缺少任一项时，正确输出是 `NOT_CERTIFIABLE`，不是缩小半径，也不是拒绝模型。

## 6. R3：最坏情况设计

对设计 \(\mathcal D\)，令 \(\mathcal Y_{\mathcal D}(\alpha)\) 为阶次
\(\alpha\) 的全部有界误差输出，\(\mathcal A_{\mathcal D}(y)\) 为
一致阶次集。已经证明：

\[
\inf_{S\ {m reliable}}\sup_y\operatorname{diam}S(y)
=
\sup_y\operatorname{diam}\mathcal A_{\mathcal D}(y)
=
\max\{|\alpha-\beta|:
\mathcal Y_{\mathcal D}(\alpha)\cap
\mathcal Y_{\mathcal D}(\beta)\ne\varnothing\}.
\]

对乘积复圆盘，定义归一化最小分离 \(\mu_{\mathcal D}(\delta)\)，则

\[
\Delta(\mathcal D)<\delta
\quad\Longleftrightarrow\quad
\mu_{\mathcal D}(\delta)>2.
\]

紧设计域上最优设计存在。该阈值把“选什么心率、多少搏、哪些谐波”
转化为可审计的最坏集合直径问题。

近心率分析还表明：

- 物理两率恢复映射的局部最坏 Jacobian 为 \(\Theta(h^{-1})\)；
- 三率差商恢复映射为 \(\Theta(h^{-2})\)；
- 后者只是该差商算法的条件数，不是所有三率估计器的信息下界。

因此心率窗口不能过近，数值设计必须与 R2 的误差半径共同优化。

## 7. R4：可保留的创新性

不能宣称为创新的内容：

- 一般分数阶频域辨识；
- 有界复数/圆盘集合参数估计；
- 一般鲁棒或集合成员实验设计；
- 单纯将分数阶参数用于 ECG 特征。

最强近邻包括：

- Malti 等对同构模型 \(K/[1+(\tau s)^\nu]\) 做了有界复频域集合反演：[CNSNS 2010](https://doi.org/10.1016/j.cnsns.2009.05.005)；
- Khemane 等已研究复矩形、极坐标与圆盘集合估计：[Signal Processing 2012](https://doi.org/10.1016/j.sigpro.2011.12.008)；
- 分数阶鲁棒实验设计已有工作：[IEEE TAC 2017](https://doi.org/10.1109/TAC.2016.2614910)；
- 一般集合成员最坏实验设计已有系统结果：[Automatica 2020](https://doi.org/10.1016/j.automatica.2020.109036)；
- ECG 分数阶频域参数已被用于 PVC 分类特征：[Biomedical Engineering 2021](https://doi.org/10.1515/bmt-2020-0170)。

当前可辩护的窄创新组合是：

> 在固定无量纲化和公共相位协议下，对未知公共复形态/导联 gauge、
> 公共正实阻尼的 ECG 分数阶形态通道，证明两心率必要且充分的全局
> 商空间辨识，并把它嵌入有限记录有界误差、多率多谐波一致性证书和
> 最坏集合直径实验设计。

论文不得使用宽泛的优先权表述，也不得把圆盘、集合辨识或实验设计本身列为新方法。

## 8. 下一阶段建议

按理论优先顺序执行：

1. **R1-CERT**：用定向舍入区间复几何替换当前浮点外排；加入 rate-ID、归一化尺度和阻尼漂移 gauge。
2. **R2-PROTOCOL**：形成可填写的硬件/数据协议表，逐项获得上述七类独立上界。
3. **R3-SOLVER**：实现 \(\max_{\mathcal D}\mu_{\mathcal D}(\delta)\) 的全局分支定界，并只使用 R2 产生的半径。
4. 只有 R1 与 R2 同时闭合后，才进入公开 ECG 的预注册验证；深度模型放在此后作为比较或中心估计器。

## 9. 本阶段验证

- 全量单元测试：`33 passed`；
- Ruff 静态检查：通过；
- 新增 Python 文件格式检查：通过；
- 500 组随机真值压力检查：
  - 真值固定参数误判不可行：`0`；
  - 包含真值的小参数盒误外排：`0`；
- 理论稿 LaTeX 连续两次编译：通过，无未解析引用、Overfull 或 Underfull；
- 生成稿：10 页，约 294 KiB，并抽查首页、中间公式页、门控/参考文献页。

这些数值测试只验证原型一致性，不能替代 R1 所缺的定向舍入证明。

## 10. 本阶段资源与边界

- 未下载任何数据、模型、软件或新文献；
- 未触发任何超过 100 MB 或大小未知的下载；
- 未复制上级 `PUBLIC` 数据；
- `main` 未修改；阶段成果只进入 `dev-codex`。

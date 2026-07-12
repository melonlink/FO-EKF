# FO-EKF R1-CERT 与 R2 协议闭环阶段报告

日期：2026-07-11

分支：`dev-codex`

阶段边界：继续理论闭环和数据前协议验证；不读取真实 ECG，不生成患者阶次结论，不引入大模型。

## 1. 阶段结论

本阶段选择了“保留现有分数阶谐波代理 + 增加严格集合认证与数据到圆盘协议”的复合路线，没有另换滤波器或叠加黑盒模型。

结论是：

> **论文路线有条件可行，理论创新已形成可投稿的核心结构；R1 认证核心通过，R2 工具链通过，但真实 ECG 的独立数值误差界仍未闭合。**

可辩护的创新不是“首次使用 FO-EKF”，也不是一般区间算法，而是以下组合：

1. 公共复形态/导联 gauge 下的两心率全局阶次辨识；
2. 带显式 rate-ID、共享阻尼漂移和非凸复环带的多率多谐波集合反演；
3. 具有 `INNER_INCLUDED / OUTER_EXCLUDED / UNKNOWN` 三值语义的阶次投影夹逼；
4. 从有限 ECG 记录到复响应圆盘的因果、可预注册、确定性/概率性分离协议。

该组合比继续增加一个新滤波方法更符合控制科学与控制工程论文的理论主线。

## 2. 并行工作流

| 工作流 | 内容 | 结果 |
|---|---|---|
| W1 理论证明与完备性边界 | gauge、参数盒 hull/core、q-box 排除、投影夹逼、UNKNOWN 终止语义 | 完成 |
| W2 R1 严格实现与独立审计 | Python-FLINT Arb/ACB、共享 rate-ID、反例、并发和分支覆盖 | 认证核心通过 |
| W3 R2 协议与对抗审计 | 冻结哈希、11 项半径传播、泄漏/概率/单位/因果门 | 工具链通过 |
| W4 论文与复现检查 | 理论稿、协议规范、模板、环境和 PDF 检查 | 完成 |

## 3. R1-CERT：已完成的理论与实现

### 3.1 认证坐标

- 分数幂只按主支实公式计算，不调用通用复数幂；
- 所有谐波共享不可变的 rate-ID 和同一个基础频率向量；
- 同一心率的阻尼偏移在所有谐波中是同一个变量；
- 选定参考心率并令其偏移精确为零，从而消除
  \(\bar\lambda+\delta\lambda_r\) 的平移 gauge；
- 参数盒内所有物理阻尼必须严格为正。

锚定坐标与物理阻尼向量一一对应，因此不改变模型分母、可行集或阶次投影。

### 3.2 参数盒外包与内核

对每个速率--谐波圆盘，Arb/ACB 计算：

- 覆盖参数盒内全部圆盘的 outer disk；
- 在半径严格非负时，包含于参数盒内每个圆盘的 uniform inner core。

普通浮点几何只能提出形态系数候选。`INNER_INCLUDED` 只有在 Arb/ACB 再次证明同一候选位于所有 inner core 和形态环带内时才能签发。

### 3.3 多圆盘 q-box 空交证书

快速外排包括：

1. 两个 outer disk 严格分离；
2. 一个 outer disk 严格位于环带内孔；
3. 一个 outer disk 严格位于 \(q_{\max}\) 外。

若三项均不能决定，则对
\([-q_{\max},q_{\max}]^2\) 做精确有理二分。每个 q-box 只有在严格证明位于环带外或至少一个 outer disk 外时才删除。全部 q-box 被删除才签发空交证书。

这已覆盖“三圆两两相交但总交为空”的经典反例。相切、零裕度、深度或节点预算耗尽均保留 `UNKNOWN`。

### 3.4 阶次投影夹逼

参数分支树叶节点分为 \(I,O,U\)：

\[
\mathcal A_{\rm in}
=\bigcup_{B\in I}\operatorname{proj}_\alpha B,\qquad
\mathcal A_{\rm out}
=\bigcup_{B\in I\cup U}\operatorname{proj}_\alpha B,
\]

\[
\mathcal A_{\rm in}\subseteq
\operatorname{proj}_\alpha F
\subseteq\mathcal A_{\rm out}.
\]

分支路径会从初始盒重新播放，并检查 prefix-free 与精确 Kraft 和。仅做“叶盒体积相加等于父盒”不再被当作无缺口证明。

只有 \(U=\varnothing\) 时，一次运行的投影才能称为 exact。有限资源下保留 UNKNOWN 是正确性的一部分，不是算法失败。

### 3.5 审计后修正

独立审计实际复现并修复了两条高风险路径：

- 一旦同一精度出现 outer/inner 内部矛盾，结果会粘性保持 UNKNOWN；升精度和继续分支都不能重新签发证书；
- 覆盖审计改为完整二叉路径重放，能够拒绝“等体积重叠 + 等体积缺口”的伪覆盖。

## 4. R2：从冻结证据到圆盘的 fail-closed 合同

### 4.1 协议与结果分离

- `protocol.config_sha256` 只哈希预注册部分，明确排除运行后 `window_result`；
- 每个结果保存所引用的 `protocol_sha256`；
- 每个结果另有排除自引用字段后的 `result_sha256`；
- 修改运行结果不会改变协议哈希，但会破坏结果哈希；
- 协议仍需由 Git commit/tag 或外部冻结记录提供时间锚，校验器本身不能证明一个哈希是在何时发布的。

### 4.2 11 项单向半径传播

对窗口--谐波 \(w\)，使用认证上界

\[
\kappa_w=\frac{\overline a_w}{\underline H_w},\qquad
\overline a_w\ge\sqrt{(G^{-1})_{mm}}.
\]

history、dwell、HRV、phase-anchor、sampling、measurement-noise 和 model-residual 的残差 W 范数先乘 \(\kappa_w\)；window leakage、quadrature、delay 和 filter initialization 使用已在系数层给出的半径。

所有传播、逐项和、总半径和 union-bound 概率都只允许向上舍入。报告值小于重算下界即为 `NOT_CERTIFIABLE`，不使用对称 `isclose` 签发覆盖结论。

测量噪声、输出端模型余项、采样误差和相对延迟还从原始校准量重新检查 `computed_*` 下界，因而不能把非零原始误差隐藏成零。

### 4.3 信息边界

校验器现在 fail closed 地检查：

- 固定 delayed-causal WLS、恒频名义相位和 R 峰左右锚；
- 无逐窗幅值归一化、无逐窗自由相位、无 acausal 预处理；
- 同一 R1 联合组的 subject、lead、physical gain 一致；
- 校准、辨识和验证 split 不同，各 bound 不得引用辨识/验证 split 作为校准；
- 来源类型使用固定 allowlist，不能在配置中临时加入“目标拟合残差”；
- Gram 最小特征值为外向认证下界、条件数为外向认证上界，并保存证书哈希；
- 全局 `computed_*` 明确覆盖全部满足冻结窗口设计的窗口；
- 概率方法不得从 union bound 偷换成无证据的 direct joint bound；
- 概率覆盖范围必须同时覆盖全部进入 R1 的 rate--harmonic 圆盘；
- R2 永远不能单独返回 `REJECT_MODEL`。

未填写模板被实测返回 `NOT_CERTIFIABLE`，退出码为 2。

## 5. 验证证据

### 5.1 自动化

- 全仓：107 passed；
- R1-CERT 专项：21 passed；
- R2 协议与对抗专项：53 passed；
- Ruff 全仓 lint：通过；
- 本阶段变更 Python 文件的 Ruff format check：通过；全仓格式检查仍提示本阶段未改动的历史文件 `src/fo_ekf/certificate.py`；
- Python-FLINT 0.9.0 官方自测：227 tests + 3,092 doctests 通过；
- 理论稿：`pdflatex` 两次编译，12 页，无未解析引用、Warning、Overfull 或 Underfull。

### 5.2 R1 独立压力验证

固定随机种子 20260711：

- 240 个含解析真点的随机参数盒，false OUTER = 0；
- 173 个 INNER 证书，抽取 11,072 个参数点并检查 66,432 次实际前向约束，违例 = 0；
- 120 次含已知可行 q 的 q-box 检查，误排 = 0；
- 两圆外切和 \(q_{\max}\) 相切保持 UNKNOWN，false OUTER = 0；
- 16 线程共 320 次混合精度签发，关系与精度不一致 = 0，`flint.ctx.prec` 无泄漏；
- 40 次真实 branch、140 个叶节点，路径重放/prefix/Kraft 全通过，1,200 个随机点全部被覆盖。

这些压力结果是实现回归证据，不替代理论中的全称证明。

### 5.3 R2 对抗回归

下列原有假 PASS 路径均被拒绝：

- 来源上界很大而结果半径人为缩小；
- 11 项半径和总半径全部填零；
- calibration split 指向 identification/validation split；
- 不同受试者、导联或增益窗口进入同一个联合组；
- 任意估计器、自由相位、非 R 峰锚点或固定秒窗；
- 把冻结的 union bound 偷换成无证据的 direct joint bound；
- 总半径或联合失败概率利用浮点容差向下漏算；
- 非零噪声、模型、采样或延迟原始界被隐藏为零 `computed_*`。

## 6. 当前 gate

| Gate | 当前判定 | 未闭合项 |
|---|---|---|
| R0 全局结构辨识 | PASS | 仅在声明模型和紧先验内成立 |
| R1 集合认证 | CERTIFIED CORE PASS | UNKNOWN 叶、相切、有限预算；尚无参数化 interval-Newton 内投影 |
| R2 协议工具 | TOOLING PASS | 真实数据的独立数值界和实际冻结记录 |
| R2 真实 ECG 圆盘 | NOT_CERTIFIABLE | history/HRV/phase/tail/sampling/delay/noise/model 等校准 |
| R3 最坏歧义设计 | THEOREM PASS | 全局认证设计优化器 |
| R4 创新边界 | NARROW PASS | 不宣称通用 FO 辨识或区间算法优先权 |

## 7. 资源和数据边界

- 本阶段没有读取或复制真实 ECG；
- 没有下载数据或模型；
- 新增的 Python-FLINT 0.9.0 wheel 为 9.15 MiB；
- 版本发现阶段下载的 0.8.0 wheel 为 10.4 MB；
- 两者均低于 100 MB 审批阈值；
- 后续任何单个或计划批量超过 100 MB、或大小未知且可能超过 100 MB 的下载仍须先审批。

## 8. 论文可行性与下一理论阶段

当前材料足以形成一篇理论导向论文的骨架，但还不能形成真实 ECG 结论。论文应以“有限记录下、带共享速率参数和三值证书语义的分数阶 ECG 集合辨识”为题眼，而不是以通用 FO-EKF 或数字孪生为创新标题。

下一阶段仍应优先理论，不应立刻进入大规模数值或深度模型：

1. 构造参数化 interval Newton/Krawczyk tube，证明
   \(\forall\alpha\in A,\exists(\lambda(\alpha),q(\alpha))\)，解决 uniform fixed-q inner 过强的问题；
2. 为 R1 证书保存输入摘要、严格裕量、q-box 覆盖树和独立 replay checker；
3. 实现 R2 的严格 Gram、Mittag--Leffler history/dwell、采样和延迟 bound engine；
4. 上述门通过并冻结真实协议后，才读取校准 split；最终 identification split 只作盲验证；
5. 深度模型以后只能作为预测基线或在独立校准集上构造明确的概率界，不能反向设定确定性半径。

因此下一步推荐执行：**T1 参数化 interval-Newton/Krawczyk 内投影证书 + 可重放证书格式**。

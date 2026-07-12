# FO-EKF T1 解析选择器可行管与重放证书阶段报告

日期：2026-07-12

分支：dev-codex

阶段边界：继续理论证明和无数据认证；不读取真实 ECG，不训练模型，不下载软件或数据。

## 1. 阶段结论

T1 已把上一阶段过强的统一见证

\[
\exists q\ \forall\alpha\in A
\]

修正并强化为论文实际需要的量词：

\[
\boxed{
\forall\alpha\in A\ \exists
(\bar\lambda(\alpha),\delta\lambda(\alpha),q(\alpha))
}.
\]

每个签发结果均由 Arb 对原始圆盘、形态环带和物理阻尼约束重新验证。

文献复核表明，参数依赖严格解包络、存在量词投影和 validated continuation 均已有先例，不能把通用 interval Newton/Krawczyk 当作创新。当前模型具有更强的特殊结构：消去固定坐标并锚定 damping gauge 后，中心残差对 nuisance 变量是仿射的，因此可得到常正定 Gram 矩阵和显式解析选择器。Krawczyk 在本模型中退化为精确可核验的线性恒等式。

> **T1 理论链与机器证书闭环通过。论文的可辩护创新空间得到加强，但仍是“FO-ECG 特化理论 + 内外投影夹逼 + 有来源的有限记录误差链”的组合创新，不是通用 Krawczyk、通用分数阶辨识或 FO-EKF 首创。**

## 2. 先验工作与路线修正

- [Rump, Acta Numerica 2010](https://doi.org/10.1017/S096249291000005X)：参数区间内逐参数存在唯一根的 interval verification 语义；
- [Neumaier, JMAA 1989](https://doi.org/10.1016/0022-247X(89)90357-0)：参数依赖方程的严格灵敏度和仿射解包络；
- [Goldsztejn and Jaulin, CP 2006](https://doi.org/10.1007/11889205_16)：存在量词约束的内外投影；
- [Ponleitner and Schichl, JGO 2021](https://doi.org/10.1007/s10898-021-01082-3)：参数盒内可行区域和解曲线包络；
- [Duff and Lee, ISSAC 2024](https://doi.org/10.1145/3666000.3669699)：parametric Krawczyk homotopy tracking。

因此本阶段采用“模型结构优先”的复合路线：

1. 用 FO-ECG 特有的 gauge-reduced 仿射结构构造解析 selector；
2. 用 exact-rational Gram 逆证明 selector 方程；
3. 用 Arb 独立验证原始不等式；
4. 与既有 outer disk / q-box cover 组合成更强的投影 sandwich。

## 3. 新增理论结果

### 3.1 固定坐标的精确消元

所有零宽 damping 和 offset 坐标先固定并代入中心项。若剩余 \(p\) 个自由阻尼坐标和 \(M\) 个复形态系数，则 reduced nuisance 向量维数为

\[
n=p+2M.
\]

固定坐标消元是精确仿射限制，不改变可行纤维或阶次投影。满秩只需在 reduced \(B_J\) 上检查；必要行数条件变为

\[
2RM\ge p+2M.
\]

### 3.2 全自由阻尼时的显式秩判据

当全部 \(R\) 个物理速率阻尼均自由、且所有响应中心非零时，中心残差矩阵失去满列秩，当且仅当存在与谐波无关的非零实比例

\[
\widetilde Z_{r,m}
=\kappa_r\widetilde Z_{r_0,m},
\qquad
\kappa_r\in\mathbb R\setminus\{0\}.
\]

该判据把 rate separation、跨谐波相位信息和 Gram 可逆性连接起来，而不是只在数值上调用一次矩阵求逆。

### 3.3 解析 least-squares selector

中心残差的实表示为

\[
r(\alpha,y)=By-g(\alpha).
\]

冻结正权重 \(W\)，若 reduced \(B\) 满列秩，则

\[
H=B^\top WB\succ0,
\qquad
y_\star(\alpha)=H^{-1}B^\top Wg(\alpha)
\]

是唯一的 least-squares selector。这里的唯一性只属于冻结的 selector 方程，不代表物理 nuisance 唯一，更不代表 \(\alpha\) 唯一。

主值分数幂满足

\[
\frac{d^k}{d\alpha^k}(i\nu)^\alpha
=(i\nu)^\alpha(\log\nu+i\pi/2)^k,
\qquad k=1,2,
\]

因此 \(y_\star\) 是候选阶次域某个开邻域上实解析映射的限制。

### 3.4 条件数和二阶 predictor 界

令

\[
\mu=\lambda_{\min}(H),\qquad
\sigma=\|B^\top W\|_2,\qquad
L_{r,m}=\sqrt{(\log\nu_{r,m})^2+\pi^2/4}.
\]

则

\[
\|y_\star^{(k)}(\alpha)\|_2
\le
\frac{\sigma}{\mu}
\left(
\sum_{r,m}
|\widetilde Z_{r,m}|^2
\nu_{r,m}^{2\alpha}L_{r,m}^{2k}
\right)^{1/2},
\qquad k=1,2.
\]

响应阵列逼近实比例可分离形态时，\(\mu\to0\)，selector 变化会被放大，宽阶次区间更难认证。

对中心 \(\alpha_0\) 和半宽 \(h\)，二阶 predictor 满足

\[
y_{\star,j}(\alpha)
\in
y_{\star,j}(\alpha_0)
+y'_{\star,j}(\alpha_0)(\alpha-\alpha_0)
+[-M_jh^2/2,M_jh^2/2].
\]

### 3.5 Exact-Gram Krawczyk 特化

输入浮点数按精确二进制有理数解释，使用精确有理矩阵 \(C=H^{-1}\)。对任意中心 \(y^0\) 和搜索盒 \(X\)，参数 Krawczyk 像为

\[
\begin{aligned}
K
&=y^0-C(Hy^0-\boldsymbol h(A))
+(I-CH)(X-y^0)\\
&=C\boldsymbol h(A),
\end{aligned}
\]

因为 \(CH=I\) 且 \(I-CH=0\)。state 方向固定点映射的收缩因子精确为零；这不表示参数 \(\alpha\) 方向发生收缩。重放时重新构造 exact-rational \(H^{-1}\)，不信任近似逆。

### 3.6 原始不等式的第二层证明

selector 方程本身不是生理约束。只有继续证明

\[
\sup_{\alpha\in A}
|q_{\star,m}(\alpha)
-d^\star_{r,m}(\alpha)\widetilde Z_{r,m}|
\le
\inf_{\alpha\in A}
\left(
|d^\star_{r,m}(\alpha)|\epsilon_{r,m}
+\rho^q_{r,m}
\right),
\]

\[
q_{\min,m}
\le
\inf_{\alpha\in A}|q_{\star,m}(\alpha)|
\le
\sup_{\alpha\in A}|q_{\star,m}(\alpha)|
\le q_{\max,m},
\]

并证明完整物理 damping path 位于先验域，才能签发内投影。残差直接按

\[
r_\star(\alpha)
=(BH^{-1}B^\top W-I)g(\alpha)
\]

计算，以保留 selector 与 \(\alpha\) 的相关性。

## 4. 结果语义

T1 只允许三种结果：

- **ROBUST_INNER**：全部坐标、圆盘和环带约束具有严格正裕量；
- **CLOSED_INNER**：闭不等式已被严格证明，但至少一项裕量为零；
- **UNKNOWN**：失秩、越界、无法定号、相切、精度或分支预算耗尽。

T1 不存在 OUTER_EXCLUDED 输出。selector 不可行并不排除另一个 nuisance witness 可能可行。

## 5. 重放证书与 alpha-only 分支

单区间证书保存完整输入 manifest 和 SHA-256、后端版本、exact-rational variable box 和预条件器、contraction norm、analytic image、原始约束裕量以及精度记录。

SHA-256 只提供内容完整性，不提供签名、作者身份或冻结时间证明。Git commit/tag 或外部注册仍负责时间锚。

alpha-only 分支具有以下约束：

- 只二分 \(\alpha\)，不先分裂全部 nuisance 坐标；
- 每个叶节点保存独立且绑定该叶输入的 JSON 证书；
- 路径必须 prefix-free，exact Kraft 和必须等于 1；
- 从根盒重放每条二叉路径，验证无缺口；
- 只有数值不可判的 UNKNOWN 继续细分；
- 结构性或资源性 UNKNOWN 立即终止；
- 深度或叶数预算耗尽时保留 UNKNOWN。

## 6. 验证证据

- T1 专项：17 passed；
- 全仓：124 passed；
- Ruff lint：通过；
- Ruff format check：通过；
- git diff check：通过；
- 理论证明独立审计：PASS；
- 理论稿双遍编译：16 页，无未解析引用、Overfull 或 Underfull；
- PDF 视觉检查：通过。

关键对抗反例：

1. 行数大于变量数但响应实比例可分离，Gram 仍精确奇异，返回 UNKNOWN；
2. 两个盘的中心均不属于另一个盘，但 least-squares 中点位于透镜交集；
3. 端点圆盘严格分离，不存在统一固定 \(q\)，旧 uniform-inner 为 UNKNOWN，而解析 \(q(\alpha)\) 为 ROBUST_INNER；
4. 零半径精确 residual 只签 CLOSED_INNER；
5. 数学边界无法数值定号时保持 UNKNOWN；
6. \(10^{-9},1,10^9\) 共同测量尺度均保持关系；
7. 父区间不可判但两个子区间均可认证，覆盖端点无缝且逐叶可重放；
8. 修改预条件器、输入、关系或精度，即使重新计算摘要也无法绕过 replay；
9. 结构性失秩不会被无意义地细分至叶预算。

## 7. 创新边界

可以主张：

1. gauge-reduced、多心率、多谐波、主值支路固定的 FO-ECG 可行集特化；
2. 全自由阻尼响应阵列的实比例可分离秩判据；
3. 解析 nuisance selector、conditioning bound 与原始圆盘裕量的两层证明；
4. analytic inner tube 与 q-box outer cover 的阶次投影夹逼；
5. finite-record data-to-disk 来源链与 R1 三值语义的组合。

不能主张：

- 首次提出 interval Newton、Krawczyk 或 validated continuation；
- 首次进行存在量词参数投影；
- 首次进行分数阶 set-membership 辨识；
- selector 唯一等于物理参数唯一；
- tube 失败等于模型被排除；
- 仍有 UNKNOWN 时投影已经 exact；
- 当前结果已经是患者级 cardiac digital twin。

## 8. 当前 gate

| Gate | 当前状态 | 未闭合项 |
|---|---|---|
| R0 结构辨识 | PASS | 限于声明模型和先验 |
| R1 outer / uniform inner | CERTIFIED CORE PASS | 相切与有限预算保留 UNKNOWN |
| T1 analytic inner tube | PASS | selector 是充分构造，不保证覆盖全部可行纤维 |
| T1 replay / alpha cover | PASS | 自摘要不提供外部时间锚 |
| R2 协议工具 | TOOLING PASS | 独立物理/数值界尚未校准 |
| R2 真实 ECG | NOT_CERTIFIABLE | 未进入真实数据阶段 |
| 创新边界 | NARROW PASS | 不宣称通用区间算法优先权 |

## 9. 资源边界

- 没有读取、复制或处理真实 ECG；
- 没有下载论文 PDF、数据、模型或软件；
- 仅向文献 manifest 增加 metadata-only 条目；
- 没有触发 100 MB 审批阈值；
- 后续任何单项或计划批量超过 100 MB、或大小不明且可能超过 100 MB 的下载，仍必须先申请审批。

## 10. 下一理论阶段建议

下一阶段仍应先做理论，而不是直接启动大数据或深度模型：

1. 实现可重放的 \(\lambda_{\min}(H)\) 下界和二阶 predictor tube，使 conditioning theorem 进入机器证书；
2. 推导“R2 每项圆盘半径预算到 T1 robust margin”的单调鲁棒性定理，给出数据质量必须达到的显式阈值；
3. 完成 Mittag--Leffler history/dwell、采样、延迟和 Gram 的独立 bound engine；
4. 只有当解析 selector 在预注册无数据压力域中覆盖率不足时，再进入 nonlinear max-margin/KKT selector；
5. T2 通过后才读取 calibration split，最终 identification/validation split 继续保持盲法。

这条路线保持论文以控制理论和严格数值证明为主，数据只承担后续验证理论的角色。

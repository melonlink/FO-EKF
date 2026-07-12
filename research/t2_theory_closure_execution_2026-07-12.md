# T2 理论闭环与机器证书执行报告

日期：2026-07-12

分支：`dev-codex`

阶段边界：纯理论、确定性数值证书与对抗性单元测试；不读取真实 ECG，不训练模型，不下载数据或软件。

## 1. 阶段结论

本阶段不更换 FO-EKF 主方法，而是在 T1 的“解析 selector 内管 + 外盘投影”上增加三个相互闭合的控制理论部件：

1. 用精确有理 Gram 谱裕量和二阶解析预测器，提高 T1 内管的可验证条件数与区间紧致度；
2. 用冻结旧 selector 路径的单调裕量定理，把已独立认证的 R2 数据盘变化条件地传播到 T1 内证书；
3. 用 Mittag–Leffler 双曲包络、Hermitian Gershgorin 球界和向外舍入，建立可重放的 R2 组件界引擎。

这一路线比立即叠加新的滤波器或大模型更合适。原因是当前首要缺口不是表达能力，而是从有限记录原始假设到 T1 几何证书之间仍有未机器闭合的量词和误差传播。T2 已关闭其中可由数学和确定性算术关闭的部分，同时没有把尚未校准的真实 ECG 假设伪装成通过。

当前可支持的论文核心主张是：

> 对声明的分数阶 ECG 谐波代理，在公共 gauge、锚定阻尼坐标和独立有限记录误差合同下，可构造精确 Gram 解析 witness；可用只沿分数阶坐标的验证管证明整段阶次均有可行 nuisance witness；并可把有限记录组件界和后续响应盘扰动以可重放、失败闭合的证书链传入该 witness。

这不是“新 Krawczyk 方法”或“新 Mittag–Leffler 不等式”的优先权主张。可辩护创新仍是模型特定的 gauge 消元、秩阈值、解析 witness、内外投影语义、有限记录误差和试验设计的一体化链条。

## 2. T2-A：Gram 谱裕量与二阶解析预测器

### 2.1 保证为正的精确谱下界

实现中的无权 Gram 矩阵为

\[
H=B^\top B,\qquad C=H^{-1}.
\]

所有元素均由输入二进制浮点数的精确有理表示构造。当秩检验通过时，\(H\succ0\)。令

\[
s_C=\|C\|_\infty,\qquad \mu_0=s_C^{-1}.
\]

因为 \(C\succ0\)，

\[
\lambda_{\min}(H)^{-1}
=\lambda_{\max}(C)
\le \|C\|_\infty,
\]

故 \(\mu_0\le\lambda_{\min}(H)\) 是永远为正、无需浮点特征值求解器的基线证书。

为提高紧致度，对有理候选 \(\widehat\mu\) 精确分解

\[
H-\widehat\mu I=LDL^\top.
\]

若 \(L\) 为单位下三角且全部有理 pivot \(d_i>0\)，则

\[
x^\top(H-\widehat\mu I)x=(L^\top x)^\top D(L^\top x)>0
\]

对全部 \(x\ne0\) 成立，因而 \(\widehat\mu<\lambda_{\min}(H)\)。有限次有理二分只负责锐化；资源耗尽退回 \(\mu_0\)，不因锐化失败产生 `UNKNOWN`。

同时保存

\[
h_\infty=\|H\|_\infty,\qquad
\frac{\sigma^2}{\lambda_{\min}(H)^2}
\le\frac{h_\infty}{\widehat\mu^2},
\quad \sigma=\|B^\top\|_2,
\]

从而机器证书不仅说明“满秩”，还给出 selector 导数增益的单向上界。

### 2.2 二阶 predictor

对

\[
y_\star(\alpha)=P g(\alpha),\qquad P=H^{-1}B^\top,
\]

固定区间中点 \(\alpha_0\) 和半宽 \(h\)。由

\[
\frac{d^k}{d\alpha^k}(i\nu)^\alpha
=(i\nu)^\alpha(\log\nu+i\pi/2)^k,\quad k=1,2,
\]

Arb 可分别包住 \(y_\star(\alpha_0)\)、\(y_\star'(\alpha_0)\) 和整区间上的 \(y_\star''\)。若

\[
U_j\ge\sup_{\alpha\in A}|(y_\star''(\alpha))_j|,
\]

则

\[
(y_\star(\alpha))_j
\in y_{\star,j}(\alpha_0)
+y'_{\star,j}(\alpha_0)[-h,h]
+[-U_jh^2/2,U_jh^2/2].
\]

同一构造也用于保持相关性的残差路径

\[
r_\star(\alpha)=(BP-I)g(\alpha).
\]

最终 enclosure 是 direct interval 与 Taylor2 enclosure 的交；真实解析路径同时属于两者，因此交集仍然可靠且不会变宽。replay 重新计算两个来源、交集和全部原始 disk/annulus 不等式。

## 3. T2-B：R2 到 T1 的冻结路径裕量定理

对旧 T1 witness，记

\[
\mu_{r,m}(\alpha)=
|d^\star_{r,m}|\epsilon_{r,m}+\rho^q_{r,m}
-|q^\star_m-d^\star_{r,m}\widetilde Z_{r,m}|
\ge\underline\mu_{r,m}.
\]

新盘满足

\[
\widetilde Z'=\widetilde Z+\xi,\qquad
\epsilon'=\epsilon+a,\qquad
\rho^{q\prime}=\rho^q+b.
\]

若

\[
0<d^-\le|d^\star|\le d^+,\quad
|\xi|\le\zeta,\quad a\ge\underline a,\quad b\ge\underline b,
\]

令 \(c=\zeta-\underline a\)，并定义

\[
L^q=
\begin{cases}
d^+c-\underline b,&c\ge0,\\
d^-c-\underline b,&c<0.
\end{cases}
\]

则三角不等式直接给出

\[
\mu'_{r,m}(\alpha)\ge\underline\mu_{r,m}-L^q_{r,m}.
\]

全部新剩余裕量严格为正时输出 `PRESERVED_ROBUST`；全部非负时输出 `PRESERVED_CLOSED`；任何不确定或负预算只输出 `UNKNOWN`。该结论固定旧 selector 路径，量词为

\[
\exists y^\star(\cdot)\ \forall\omega\ \forall\alpha\in A,
\]

比“对每个扰动重新找一条路径”更强，也使 replay 不会暗中切换 witness。

必须特别记录方向性：

- 中心不变而响应半径或独立复 morphology 半径增大，会保留 inner feasibility；
- 但更大的盘削弱模型排除能力并扩大 outer projection；
- 因此该定理是中心漂移、半径收缩和先验收紧的容许预算，不是孤立的数据质量指标；
- 新中心会改变最小二乘矩阵，不能用新中心重建 selector 后冒充旧证书转移。

## 4. T2-C：R2 组件界引擎

### 4.1 Mittag–Leffler 的可计算统一包络

Pollard 证明了 \(E_\alpha(-x)\) 在 \(0<\alpha\le1\) 下完全单调。Simon 的 Theorem 4 证明，对 \(0<\alpha<1,\ x\ge0\)，

\[
\frac{1}{1+\Gamma(1-\alpha)x}
\le E_\alpha(-x)
\le\frac{1}{1+\Gamma(1+\alpha)^{-1}x}.
\]

由于 \(\Gamma(1+\alpha)\le1\)，可取更简单但完全可靠的

\[
E_\alpha(-x)\le\frac1{1+x}.
\]

对 \(I=[\alpha_-,\alpha_+]\)、\(\lambda\ge\lambda_-\)，定义

\[
p_I(t)=
\begin{cases}
0,&t=0,\\
t^{\alpha_+},&0<t\le1,\\
t^{\alpha_-},&t\ge1,
\end{cases}
\qquad
R_I(t)=\frac1{1+\lambda_-p_I(t)}.
\]

则在完整 \((\alpha,\lambda)\) 盒上

\[
0<E_\alpha(-\lambda t^\alpha)\le R_I(t).
\]

由此得到

\[
B_{\rm hist}(s)
\le
\left(|z_0|+\frac{U_0}{\lambda_-}+C_{\rm pre}\right)R_I(S+s),
\]

\[
B_{\rm dwell}(s)
\le
\frac{U_{\rm pre}+U_0}{\lambda_-}R_I(s).
\]

第二式舍弃了一个非负的相减项，较保守但不需要不受控的 Mittag–Leffler 数值积分。

主要文献：

- H. Pollard, 1948, [DOI](https://doi.org/10.1090/S0002-9904-1948-09132-7)；
- T. Simon, 2014, [DOI](https://doi.org/10.1214/EJP.v19-3058)。

### 4.2 WLS Gram 与其余组件

对 Hermitian Gram 的球包络，令

\[
g_-=\min_j\left(\underline G_{jj}
-\sum_{k\ne j}\overline{|G_{jk}|}\right),
\quad
g_+=\max_j\left(\overline G_{jj}
+\sum_{k\ne j}\overline{|G_{jk}|}\right).
\]

若 \(g_->0\)，Hermitian Gershgorin 定理给出

\[
\lambda_{\min}(G)\ge g_->0,\qquad
\kappa_2(G)\le g_+/g_-,\qquad
\sqrt{(G^{-1})_{mm}}\le1/\sqrt{g_-}.
\]

若该充分条件失败，只能说明当前包络不够紧，不能宣布 Gram 奇异。

同一引擎使用球算术重算：

- history/dwell 的样本加权范数；
- ADC、时间抖动、插值和抗混叠残差的 sampling 范数；
- 离散 WLS 的零求积项，或有独立 \(M_2\) 时的复合梯形余项；
- 相对幅值、通道群延迟和相位偏差的 delay 项；
- Gram 最小特征值、最大特征值、条件数和 WLS amplification。

这里必须区分两个不能混用的算子合同：

- `discrete_wls` 使用样本 Gram、\(G^{-1}\) 和 WLS row norm，求积项严格为零；
- `trapezoidal_continuous_lockin` 只在完整周期连续正交基和精确 composite-trapezoid 权下使用 \(G=I\)，再独立加入 panelwise \(h_j^3M_{2,j}/(12T|\widehat H_m|)\)。

history/dwell 也不能只靠目标谐波的 \(|\widehat H_m|\) 传播任意前端瞬态。机器合同因而只允许：

- 确切 identity 前端；或
- 独立提供每个样本处已经包含 \((|h|*B)(s_n)\) 的 preconvolved output envelope。

前端幅频或群延迟的单频校准不能替代这个时域卷积合同。所有 evidence 还必须绑定与请求一致的 signal unit、秒制 time unit 和适用于该组件的来源类别。

输入 manifest 把每个原始数值、来源 hash、校准 split、单位、作用域、谐波与样本顺序绑定到 SHA-256。replay 不信任派生数字，而是在证书精度下重建全部结果。该组件证书不等同于完整 `PASS_DETERMINISTIC`。

## 5. 已关闭与仍开放的理论边界

| 项目 | T2 状态 | 结论 |
|---|---:|---|
| T1 Gram 满秩后的定量谱裕量 | 已关闭 | exact inverse-inf 基线，可选 exact shifted-LDL 锐化 |
| selector/residual 的二阶路径包络 | 已关闭 | \(O(h^2)\) remainder，可重放 |
| 已认证 R2 disk 变化到旧 T1 witness | 条件关闭 | robust/closed/unknown 三值、冻结路径；bundle 绑定待 T3 |
| history/dwell 数值传播 | 条件关闭 | primitive bound 成立时可签发组件上界 |
| sampling/quadrature/delay/Gram | 条件关闭 | primitive calibration 成立时可签发组件上界 |
| 真实 ECG 的 R2 总盘 | 未关闭 | 仍缺独立生理与仪器校准 |
| 全局 R1 投影无 UNKNOWN 叶 | 未关闭 | 仍是有限资源 completeness 问题 |
| 内部心脏状态或 observer mismatch 稳定性 | 不在当前主张 | 需要另一个含全历史算子的理论问题 |

真实 ECG 前仍需外部或独立校准的项目至少包括：

- \(z_0,U_{\rm pre},U_0,C_{\rm pre}\) 的适用性与独立性；
- R 峰误差、HRV 操作相位和任何潜在生理相位误差；
- 抗混叠、\(M_1,M_2\)、时间戳抖动、插值和 ADC 界；
- 前端幅值、ECG 与 R-peak 通道对齐、相对延迟和传递函数下界；
- harmonic tail、filter initialization、measurement noise、model residual。

缺一项时仍为 `NOT_CERTIFIABLE`，不是模型被拒绝。

## 6. 验证记录

代码合并后的实际验证：

- 聚焦证书测试：`62 passed`
  - `tests/test_tube_certificate.py`：25；
  - `tests/test_margin_transfer.py`：14；
  - `tests/test_r2_bound_certificate.py`：23。
- 全仓测试：`python -m pytest -q`，`169 passed`。
- 静态检查：`python -m ruff check src tests`，通过。
- 本阶段六个 Python 文件格式检查：`ruff format --check`，通过。
- 补丁空白检查：`git diff --check`，通过。
- 论文：`latexmk -pdf -interaction=nonstopmode -halt-on-error` 两轮编译成功，19 页；无 undefined reference/citation，只有 gate 表格的一处非阻断 `Underfull \hbox`。
- PDF 视觉 QA：按共享文档运行时用 `D:\TOOLS\AI\anaconda3\python.exe`、`pdf2image` 和指定 Poppler 在 110 dpi 栅格化全部 19 页；检查总览以及谱裕量、Taylor2、R2 组件定理、margin transfer、gate 表和参考文献页面，未发现截断、重叠、越界或空白异常。

交叉审计额外发现并关闭：

1. 非 identity 前端不能仅凭目标谐波 \(|\widehat H_m|\) 传播 history/dwell；v2 现在要求 identity 或独立 preconvolved output envelopes。
2. continuous trapezoid lock-in 不能复用任意离散 WLS Gram；v2 已拆成完整周期 \(G=I\) 合同并强制 exact composite 权。
3. 深层嵌套 JSON 的 `RecursionError` 已在三类 replay 中统一失败闭合。
4. alpha branch replay 现在校验叶深度预算，并重新签发 UNKNOWN 叶来核对停止原因，不能再靠重签 JSON 伪造 resource/structural 分类。
5. margin transfer 的真实三速率、自由 damping 端到端测试验证：\(c\ge0\) 选 \(d^+\)，\(c<0\) 选 \(d^-\)，且实际使用 Taylor2 收紧后的 `selector_image`。

本阶段没有访问 `../PUBLIC`，没有训练模型，没有下载数据或软件，也没有触发超过 100 MB 的审批事项。

## 7. 下一阶段建议

进入 T3 前仍应保持“理论先行”：

1. 将 R2 组件证书与冻结协议 TOML 做显式 bundle 绑定，使 `bound_engine_required=true` 时旧的纯标量证据失败闭合；
2. 构造合成窗口的端到端 `R2 component certificate -> T1 margin transfer -> alpha tube replay`；
3. continuous lock-in 模式由样本时刻在引擎内部直接生成 composite-trapezoid 权，避免不可二进制精确表示的合法网格产生不必要的 `UNKNOWN`；
4. 对 Gram Gershgorin 失败但实际正定的样例增加 verified LDL/Cholesky 分支，降低不必要的 `UNKNOWN`；
5. 只有上述闭环通过后，再请求使用上一级 `PUBLIC` 中的小规模 ECG 清单；任何新增下载若超过 100 MB，必须先取得用户批准。

# T3 端到端证书链与数值完备性执行报告

日期：2026-07-12

分支：`dev-codex`

阶段边界：理论、机器证书、协议约束和合成对抗测试；本阶段不把真实 ECG 的未知工程界填成通过，也不把条件性证书解释成心脏内部状态或数字孪生证明。

## 1. 阶段结论

T3 已关闭 T2 留下的两个机器闭环缺口：

1. R2 冻结协议、逐窗口组件证书、结果盘、旧 T1 analytic-selector tube 和 margin transfer 现在由一个可重放 bundle 显式绑定；
2. continuous lock-in 权重与离散 WLS Gram 的两处不必要 `UNKNOWN` 已分别由“内部生成 composite-trapezoid 权”和“验证型 Hermitian LDL* 备用证明”关闭。
3. 离散路径进一步由 R2 v4 冻结整数 sample indices/window endpoints；独立
   Arb/ACB exact-integer-payload certificate 从数字样值直接包络 phase、\(G,b\)、
   joint solve 和前端除法，关闭 stored binary64 centre 的完整数值差。

因此当前实现已经能机械证明下列条件命题：

> 若冻结的 R2 协议为确定性通过，且每一个 T1 rate--harmonic 盘均有与协议、记录、估计器、单位、split、整数窗口几何、独立 primitive evidence、可重放 WLS execution 和 exact-functional interval 一致的 R2 v4 组件证书，则可从 R2 结果重建完整新响应盘族，并验证旧 T1 robust witness 在该盘族下仍被保留。

任何缺证、hash 不一致、单位或 split 泄漏、Gram/R2 半径虚报、窗口缺失、先验不一致或 margin 不足都只返回 `NOT_CERTIFIABLE`。bundle 没有 outer 或 reject 关系，因而证据链断裂不能被误读为模型排除。

这仍不是“真实 ECG 已通过 R2”。真实记录中的 history、谐波尾、R 峰误差、抗混叠、时钟、前端、噪声和 model residual 必须来自独立校准或合格概率合同；目标窗残差不能回填这些量。

## 2. T3-A：sealed R2--T1 bundle

新增 `src/fo_ekf/r2_t1_bundle.py`。签发前依次执行：

1. 重放完整 R2 protocol，并只接受 `PASS_DETERMINISTIC`；
2. 重放旧 T1 tube，并只接受 `ROBUST_INNER`；
3. 核对冻结协议中绑定的 T1 certificate hash 和 input hash；
4. 要求每个 canonical `(harmonic, rate)` 恰有一个 R2 component link；
5. 核对 protocol/component/window 的 protocol hash、record manifest、estimator manifest、signal unit、校准与辨识 split、alpha/lambda 先验、harmonic basis、频率、R 峰锚点、完整 RR 数和样本数；
6. 重放每个 R2 v4 certificate，并确认 window result 没有高报 Gram 下界、低报条件数/row norm 或低报已证组件半径；
7. 重放 hash-bound WLS execution 与 exact integer-payload interval，严格核对
   sample indices、window start/end、digital/time/weight hashes、target 和 centre；
8. 从 `z_tilde` 和 `radius_total_deterministic` 构造完整新响应盘族；
9. 对中心移动作精确二进制有理上取整，再调用 frozen-selector margin transfer；
10. 把旧 T1、全部 R2/WLS/interval components 和 margin-transfer 子证书嵌入最终 bundle，重新计算时逐个重放。

bundle 只表达：

- `PRESERVED_ROBUST`；
- `PRESERVED_CLOSED`；
- `NOT_CERTIFIABLE`。

它不检查新盘的 outer projection，也不输出 `REJECT_MODEL`。若将来需要排除结论，仍须由完整 R1 joint outer cover 在无 `UNKNOWN` 叶时单独签发。

初版 declared-disk composition 已被后续 T5 正式工件取代。当前生产链绑定原始
数字 payload hash、整数相位算子、完整 complex WLS \(G,b\)、求解结果和独立
exact-functional interval，因而合成路径可准确称为 samples-to-disk-to-frozen-
witness composition。普通 solve-only execution record 单独仍保持
`NOT_CERTIFIABLE`；只有独立 interval replay 成功后才关闭完整数值 functional。
真实 Fantasia 包没有把该合成闭环冒充 acquisition/model coverage，仍为
`NOT_CERTIFIABLE`。

## 3. T3-B：continuous lock-in 单一权重合同

R2 bound schema 当前为 v4。对 `trapezoidal_continuous_lockin`，调用者必须令 `raw_weights=None`；引擎从输入时间戳的精确二进制有理表示内部生成

\[
 w_0=\frac{t_1-t_0}{2T},\qquad
 w_j=\frac{t_{j+1}-t_{j-1}}{2T},\qquad
 w_N=\frac{t_N-t_{N-1}}{2T}.
\]

这样，同一组时间戳只有一个授权的 composite-trapezoid 算子。合法的 0.1、0.35 等非 dyadic 浮点时间不再因调用方用另一条浮点运算路径重算权重而产生假 `UNKNOWN`；proof 中保存派生权及 `weight_contract`，replay 再次从时间戳生成。

离散 WLS 仍要求调用者声明正权并由引擎精确归一化，两个算子合同不混用。

## 4. T3-C：验证型 Hermitian LDL* Gram 备用证明

Gershgorin 给出快速且可靠的

\[
 \lambda_{\min}(G)\ge g_-,\qquad
 \lambda_{\max}(G)\le g_+,
\]

但 `g_-<=0` 不代表真实 Gram 非正定。v4 保留在 v3 引入的 verified fallback：当 Gershgorin 不足以达到冻结的最小特征值和条件数阈值时，定义

\[
 \mu_{\rm req}
 =\max\left\{\mu_{\min},\frac{g_+}{\kappa_{\max}}\right\}
\]

并对 ACB 球矩阵 `G-mu_req I` 递归计算无主元 Hermitian

\[
 G-\mu_{\rm req}I=LDL^*.
\]

每个真实 pivot 都由 ball arithmetic 包含；若全部 pivot 的认证下端严格为正，则 Sylvester 准则证明

\[
 \lambda_{\min}(G)>\mu_{\rm req},\qquad
 \kappa_2(G)\le g_+/\mu_{\rm req}\le\kappa_{\max}.
\]

如果任一 pivot 的符号未决，只返回 `gram_thresholds_not_certified`。实现没有从普通浮点特征值“猜测”正定性。

对抗样例使用一个实际正定、但 Gershgorin 下端为负的非均匀四点 Gram；LDL* 路径签发并可重放。完全 alias 的奇异样例仍失败闭合。

## 5. T3-D：split 与时间泄漏协议加固

协议模板和 validator 现在强制声明：

- calibration、identification、validation 三者两两不相交；
- subject-specific identification 在 validation 之前；
- 最小时间保护间隔已给出；
- split manifest 与 record-interval manifest 均有 SHA-256；
- 每个 window result 绑定 record manifest 和 estimator manifest SHA-256；
- `bound_engine_required=true`。

这里的“不相交”分两层解释：calibration 受试者与评估受试者做 patient-level 隔离；同一评估受试者的 identification/validation 则共享 subject identity、但样本区间不重叠并保持严格时间顺序。协议用 `split_contract=calibration_subject_disjoint_then_within_subject_time` 固定该含义，避免把 subject-specific 验证与患者隔离写成相互矛盾的要求。

这些字段只证明配置内部的声明与 hash 一致，不替代外部数据审计；但它们阻断了“同一目标窗既校准误差盘、又做辨识和验证”的无痕循环。

## 6. 验证记录

合并 exact-index/WLS interval/T5/T6 代码后的实际质量门：

- 全仓 pytest：`254 passed`；
- T5/T6/Fantasia provenance 聚焦测试：`19 passed`；
- two-rate sweep + multirate 聚焦测试：`16 passed`；
- R2 protocol/bound/bundle/WLS 聚焦路径包含在全仓门中；
- Ruff lint：通过；
- 本阶段文件 Ruff format：通过；
- `git diff --check`：通过。

主要对抗测试覆盖：缺盘 link、component/record hash 不一致、evidence 语义不一致、T1 baseline 或 q prior 不一致、中心大幅移动、重签 tampered proof、重复 JSON key、非 dyadic composite 权、Gershgorin 失败但 LDL* 成功、真实 alias Gram 失败闭合，以及 record/estimator manifest 缺失。

## 7. 本地数据执行与最终边界

只读清单确认：

- Fantasia 共有 40 条长时记录；本研究选择的四条均为 250 Hz（数据库另有
  一条 `f2y02` 为 333 Hz），每条约 2 小时并带 R 峰标注，适合
  subject-specific 多心率与未来时段验证；
- 首轮固定四条 Fantasia 记录输入约 28.98 MiB，连同元数据与 hash 读取仍低于 100 MiB；
- PTB-XL fold 10 的重复患者 100 Hz 子集约 12.25 MiB，适合作为之后的跨日期个体形态迁移验证；
- 全 Fantasia、完整 PTB-XL 或 INCART 全库均超过 100 MiB，不在未审批时执行全量读取或新增下载。

后续任务均已执行：reference-pair/三率一致性与 no-go 定理已进入匿名补充材料；
T5 完整合成链在 clean production commit 上签发并经篡改对照；T6 对 1260 个
两率全局逆、共同复增益和三率一致性样例完成确定性 sweep；Fantasia 四记录包从
官方 v1.0.0 文件重新生成并实际重放 1859/1859。最终边界不变：只有独立 primitive
校准充分时，真实数据才允许进入 `PASS_DETERMINISTIC/PASS_PROBABILISTIC`；当前
四记录结果保持 `NOT_SUPPORTED`（经验模型比较）与 `NOT_CERTIFIABLE`（R2 coverage）。

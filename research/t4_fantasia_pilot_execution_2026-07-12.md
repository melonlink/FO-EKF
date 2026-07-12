# T4 Fantasia 锁定小样本与 WLS execution 执行报告

日期：2026-07-12

分支：dev-codex

## 1. 结论先行

本阶段证明了时域工程管线可以从外部 WFDB 记录稳定地重建同一组整窗复 WLS
系数，但没有得到分数阶特异性证据：

- 四条记录的自由分数阶拟合全部落在 \(\alpha=1\) 上边界；
- 相对设计锁定的嵌套 \(\alpha=1\) 模型，分数阶模型没有一个实质胜场；
- 两条锁定工程验证记录中，只有一条优于常数 Fourier 基线，也只有一条优于
  rate-shuffle 负对照；
- 因而 distinct_fractional_order_evidence 与
  confirmatory_rate_response_evidence 均为 NOT_SUPPORTED；
- 该结果只支持 segmentation、固定相位、联合 WLS、哈希绑定与重放管线可行，
  不支持心脏数字孪生、内部心脏状态或生理分数阶机制。

这是一项应保留的阴性结果。它与理论中的“有限频率只能辨识声明模型族内的有效
阶次”以及“速率特异复形态会破坏阶次可辨识性”一致，阻止项目把普通的速率相关
形态变化误报为分数阶记忆。

## 2. 冻结设计

只读输入位于仓库外部的 Fantasia 目录。使用四条预先指定记录：

| 角色 | 年轻 | 年长 |
|---|---|---|
| 工程校准 | f1y01 | f1o01 |
| 锁定工程验证 | f2y01 | f2o01 |

唯一源文件共 30,386,171 B，低于用户规定的 100 MB 审批阈值。仓库没有复制原始
波形。每个 subject 内使用前 1800 s 做 identification，2100 s 后做 validation，
中间留出 300 s 时间保护；本阶段不把四条记录称为独立临床验证 cohort。

数据源为 PhysioNet Fantasia Database v1.0.0，DOI
`10.13026/C2RG61`，文件许可为 Open Data Commons Attribution License
v1.0。发布与引用必须同时保留 Iyengar et al. (1996) 原研究和 Goldberger et al.
(2000) PhysioNet 标准引用。结果包没有波形样值，但 execution records 保存了从
官方 beat annotations 派生的选中 R 峰样本索引，因此数据署名与许可边界仍适用；
详见仓库 `DATA_ATTRIBUTION.md`。

每个窗口固定为 16 个相邻 RR 间隔。相位、基、权和系数约定为

\[
 \phi_n=2\pi N\frac{n-a}{b-a},\qquad
 \mathcal K_{\rm est}=\{-3,-2,-1,0,1,2,3\},
\]
\[
 G=\Psi^*W\Psi,\qquad b=\Psi^*Wy,\qquad
 \widetilde C=G^{-1}b,\qquad
 \widetilde Z_m=\widetilde C_m/\widehat H_m .
\]

这里 \(W\) 是一次归一化的 identity raw weights，\(\widehat H_m=1\)，且实 ECG
重建约定为

\[
 y=Z_0+2\operatorname{Re}\sum_{m>0}Z_m e^{im\phi}.
\]

模型使用 \(\tau_\star=1\,{\rm s}\)，所以
\(\nu=2\pi m f\tau_\star\) 和 \(\lambda\) 无量纲，\(q_m\) 与 ECG 同为 mV。
不得把一侧实正弦幅值或逐搏重新定相混入该约定。

当完整 16-RR 窗超过 4096 点时，使用与样值无关的冻结选择

\[
 n_j=a+\left\lfloor \frac{j(b-a)}{4096}\right\rfloor,
 \qquad j=0,\ldots,4095 .
\]

1859 个窗口中有 233 个使用该路径。数字样值、时间、权、相位和基均使用相同索引；
完整窗口与选中 payload 分别哈希。该抽样尚无抗混叠界，因此不能签发真实 R2。

## 3. 可重放证据

结果包位于 research/results/fantasia_pilot_2026-07-12/，总大小
4,727,899 B，不包含原始样值。主要对象为：

- design_lock.json：冻结约定和角色；
- source_inventory.json：外部文件大小、SHA-256 和逐记录组合哈希；
- window_manifest.csv.gz：1859 个窗口的几何与 execution 哈希；
- wls_coefficients.csv.gz：5577 个正谐波系数；
- wls_execution_records.jsonl.gz：压缩后的完整执行记录；
- model_fits.json 与 heldout_metrics.csv：四种模型/对照的时序留出结果；
- r2_status.json：逐记录缺失 primitive；
- SHA256SUMS：结果包校验。

从外部记录重新读取数字 payload 后，重放器同时检查：

1. 外部 .hea/.dat/.ecg 组合哈希与 source inventory；
2. execution 的 record/protocol/estimator/window 绑定；
3. manifest wrapper 的 record、role、split、window index、窗口起止和样本数；
4. 选中 payload、时间、权、基、\(G\)、\(b\)、\(\widetilde C\)、
   \(\widetilde Z\) 和 execution SHA-256。

最终结果为 1859/1859 MATCH。`design_lock.json` 与 `summary.json` 绑定
production implementation SHA-256
`8bda0793202aa7c57e6adb0888741efef4b482115e6d4ed2a4d736cfb7b2cea7`
以及 Python/NumPy/SciPy/WFDB/python-flint/FLINT 版本。相对于记录下来的
binary64 \(G,b,\widehat H\)，
Arb 对联合求解和除法给出 1859/1859 outward dyadic enclosure，最大绝对上界为
3.616468081368432e-15。确定性重跑后的 SHA256SUMS 文件哈希为
c6d4018f84e90a2abe7130b384e9b0ca8c5ac7859bc04bf6b19d98964d80518a。
所有未压缩的 sealed JSON/CSV/SHA 文本均以显式 LF 字节写出；因此 Git
规范化不会改变其内容，候选包与 canonical 包的 10 个文件逐字节一致。

该数值包络目前不含 \(G,b\) 的 binary64 累积舍入；它也不含 R 峰定时、采集噪声、
抗混叠、前史、有限驻留、模型余项或前端校准误差。以上缺口在证据对象中被显式列出，
而不是由目标窗拟合残差回填。

## 4. 留出结果

四种设计锁定比较为：

1. fractional_common：共享 \(\alpha,\lambda\)，谐波特异 \(q_m\)；
2. alpha_one_nested：同一模型固定 \(\alpha=1\)；
3. constant_fourier：与速率无关的复谐波均值；
4. rate_shuffle_negative_control：打乱 identification 的 rate--response 配对。

validation 汇总 RMSE（mV）如下：

| record | fractional | \(\alpha=1\) | constant | rate shuffle |
|---|---:|---:|---:|---:|
| f1y01 | 0.056631 | 0.056631 | 0.056638 | 0.056688 |
| f1o01 | 0.020851 | 0.020851 | 0.020974 | 0.020920 |
| f2y01 | 0.019044 | 0.019044 | 0.019296 | 0.019282 |
| f2o01 | 0.009714 | 0.009714 | 0.009458 | 0.009457 |

实质改善阈值冻结为相对 RMSE 0.001。自由分数阶相对嵌套
\(\alpha=1\) 的比值在四条记录上均数值等于 1，且四个自由阶次拟合返回值均为 1。
这些是冻结网格与多起点局部搜索的返回解，不是经认证的全局最优解。
锁定的 f2y01 对常数/打乱对照有约 1.30%/1.23% 改善，而 f2o01
反而差约 2.71%/2.72%。两条锁定记录方向不一致，不能构成确认性结论。

## 5. R2 状态

所有真实记录均保持 NOT_CERTIFIABLE。至少缺少：

- annotation timing error bound；
- HRV 与 whole-window phase-warp bound；
- measurement noise envelope；
- model discrepancy envelope；
- finite-dwell transient bound；
- prehistory/memory-tail bound；
- sampling/quadrature 与 subsampling anti-alias bound；
- omitted-harmonic tail/window-leakage bound；
- ADC gain 与 analog front-end calibration bound；
- \(G,b\) 累积及完整 exact-sample WLS interval enclosure。

因此本包不得改名为真实 ECG 的 PASS_DETERMINISTIC，也不得用于 R1 模型排除。
下一阶段只在独立 primitive 校准齐备时升级 R2；否则保留阴性/未知状态。

## 6. 复现命令

    .\.venv\Scripts\python.exe scripts\run_fantasia_pilot.py
    .\.venv\Scripts\python.exe scripts\replay_fantasia_wls.py --bundle-dir artifacts\fantasia_pilot_candidate
    .\.venv\Scripts\python.exe scripts\replay_fantasia_wls.py --bundle-dir research\results\fantasia_pilot_2026-07-12
    .\.venv\Scripts\python.exe -m pytest tests\test_wls_execution.py tests\test_fantasia_pilot.py tests\test_paths.py -q

生成器连续两次运行的 10 个文件逐字节一致；候选目录与 canonical 目录逐文件
SHA-256 一致。增强后的 replayer 先检查九项 `SHA256SUMS` allowlist、官方数据
哈希、config/implementation/runtime provenance，再读取外部记录；最终 canonical
重放为 1859/1859 MATCH、0 failures。最终全仓质量门见仓库根目录
`REPRODUCIBILITY.md` 与阶段提交记录。

# R2-PROTOCOL：有限 ECG 记录到复响应圆盘的预注册协议

- 日期：2026-07-11
- 适用分支：dev-codex
- 对应理论门：R2 data-to-disk
- 配置模板：config/r2_protocol_template.toml
- 当前状态：协议和模板草案已形成；机器传播审计与数值上界均未闭合，R2 仍为 NOT_CLOSED

## 1. 目的与允许结论

本协议只回答：对预注册的受试者、心率窗口和谐波，有限 ECG 记录能否给出

\[
Z_{r,m}\in\mathcal D(\widetilde Z_{r,m},\epsilon_{r,m})
\]

的可复核覆盖保证。\(Z_{r,m}\) 是声明的标量分数阶谐波代理在公共 R 峰相位 gauge 下的复响应，不是心肌内部状态或无条件的生理常数。

允许的 R2 状态只有：

- **PASS_DETERMINISTIC**：全部必需项有确定性上界，总半径有限且通过预注册宽度门；
- **PASS_PROBABILISTIC**：得到明确联合覆盖概率的概率圆盘，不得冒充确定性圆盘；
- **NOT_CERTIFIABLE**：任一必需上界缺失、来源不独立、量纲不闭合、非有限或总圆盘过宽；
- **EXCLUDE_PROTOCOL**：窗口违反信息边界，例如非因果滤波、逐窗幅值归一化或自由相位对齐。

REJECT_MODEL 不属于 R2 状态。只有 R2 圆盘通过后，R1 联合物理可行集为空，才允许写成“声明代理与预注册误差集不相容”。不得为获得 PASS 而在看过目标拟合残差后调整上界。

## 2. 量纲与目标模型

在任何参数区间和误差预算前固定参考时间 \(\tau_\star>0\)：

\[
\xi=t/\tau_\star,\qquad \Omega_r=\omega_r\tau_\star .
\]

全部分数阶动力学相对于无量纲时间 \(\xi\) 书写：

\[
{}^C D_\xi^\alpha z
=-\lambda z+\sum_{k\in\mathbb Z}q_k e^{ik\phi(\xi)},\qquad
Z_{r,k}=\frac{q_k}{\lambda+(ik\Omega_r)^\alpha}.
\]

| 量 | 单位 |
|---|---|
| \(t,T,\Delta\tau\) | s |
| \(\xi,\Omega_r,\lambda\) | 1 |
| \(\omega_r\) | rad/s |
| 相位及相位误差 | rad |
| \(z,q_k,Z_{r,k},\epsilon_{r,k}\) | 同一 ECG 测量单位 |
| 确定性能量噪声界 \(N_2\) | ECG-unit \(\sqrt{\mathrm s}\) |
| 导数界 \(M_1,M_2\) | ECG-unit/s、ECG-unit/s\(^2\) |
| 前端频响 \(H\)、WLS Gram 矩阵 | 1 |

若 ADC counts 未通过一个跨全部窗口共同的标定换成 mV，可以在 counts 中认证；但 \(q,Z,\epsilon\) 和全部幅值先验必须同时使用 counts。

## 3. 冻结的估计器与公共 gauge

### 3.1 延迟因果窗口

窗口左、右端 \(a,b\) 由 R 峰锚定并跨越 \(N\) 个完整 R--R 周期：

\[
T=b-a,\qquad
\widehat\omega_r=\frac{2\pi N}{T},\qquad
\widehat\phi(t)=\widehat\omega_r(t-a).
\]

估计值仅在右端 R 峰到达后输出；这是 delayed-causal，不是零延迟估计。不得使用 \(b\) 之后的 ECG。

一个对所有窗口共同的弧度相位偏置可吸收到有效 \(q_m\)。以下操作必须 EXCLUDE_PROTOCOL：

- 每窗单独最大化复相关或旋转到最佳相位；
- 每窗单独做 min--max、方差、峰值或能量归一化；
- 把共同秒延迟当成共同复相位；秒延迟的旋转 \(e^{-im\omega_r\tau}\) 随心率变化；
- 使用 filtfilt、中心平滑或其他依赖未来样本的预处理，却仍声称因果。

### 3.2 加权最小二乘

对预注册谐波集合 \(\mathcal M\)，令

\[
\Psi_{n,k}=e^{ik\widehat\phi(t_n)},\quad
W=\operatorname{diag}(w_n),\quad
\sum_nw_n=1,\quad G=\Psi^*W\Psi.
\]

必须验证 \(G\succ0\) 及预注册的最小特征值/条件数门。因果前端名义频响为 \(\widehat H_{r,m}\ne0\) 时，先估计前端输出系数，再在系数层校正：

\[
\widetilde{\boldsymbol C}=G^{-1}\Psi^*W\boldsymbol y,\qquad
\widetilde Z_{r,m}=\frac{\widetilde C_{r,m}}{\widehat H_{r,m}}.
\]

identity 前端取 \(\widehat H_{r,m}=1\)。公共未知导联增益可保留，但必须跨窗口相同，且 \(q,Z,\epsilon\) 使用同一有效单位。

实 ECG 必须预注册双边 Fourier 集（含 DC、正负频率并施加共轭对称），或使用经单独证明的延迟因果解析信号实现；不能默认非因果理想 Hilbert 变换。

## 4. 证据等级和独立性

每项必须记录 source_kind、source_reference、校准方法、单位、数值、适用范围和独立性。

### 4.1 确定性证据

- 传感器/ADC 的硬件规格上界；
- 与目标辨识记录独立的有界校准；
- 经证明的物理或工程包络；
- 不依赖目标拟合残差的窗口几何量，如 \(T,N,G,h_{mk}\) 和 R--R 相位偏移；
- 由上述上界经解析不等式或有保证区间运算传播的量。

### 4.2 概率性证据

Gaussian、bootstrap、经验分位数、置信区间和预测区间只能形成概率圆盘。覆盖范围必须预注册为 R1 将使用的全部窗口、心率和谐波圆盘同时覆盖。每个误差来源及其适用的窗口--心率--谐波组合记录失败概率；令 \(F_{j,r,m,w}\) 表示对应失败事件。没有直接联合校准时，只能对全部这些事件使用预注册 union bound

\[
\Pr\!\left(\bigcup_{j,r,m,w}F_{j,r,m,w}\right)
\le\sum_{j,r,m,w}\delta_{j,r,m,w}
\le\delta_{\rm joint}<1,
\]

最终必须保留

\[
\Pr\left\{
\bigcap_{(r,m,w)\in\mathcal S}
|\widetilde Z_{r,m,w}-Z_{r,m,w}|
\le\epsilon^{\rm prob}_{r,m,w}
\right\}\ge1-\delta_{\rm joint},
\]

其中 \(\mathcal S\) 是冻结的全部选定圆盘集合。不能删除概率符号，也不能只给逐圆盘概率后把 R1 联合空集写成确定性拒绝。确定性项的失败概率为 0。

### 4.3 禁止反向定界

以下任一情况均为 NOT_CERTIFIABLE：

- 用同一窗口、同一代理的拟合残差最大值设定模型余项，再认证同一次拟合；
- 看到阶次区间过宽后调小噪声、前史、漂移或相位上界；
- 先删除残差大的窗口，再称剩余样本满足预注册覆盖；
- 将交叉验证 RMSE、标准差或 SNR 点估计直接当确定性上界；
- 校准与验证数据未按冻结清单隔离，或拆分存在受试者泄漏。

## 5. 总 data-to-disk 预算

令

\[
a_m=\sqrt{(G^{-1})_{mm}},\qquad
\|\boldsymbol r_j\|_W\le\bar r_j.
\]

确定性系数半径取

\[
\begin{aligned}
\epsilon_{r,m}^{\rm det}
={}&\frac{a_m}{|\widehat H_{r,m}|}
\left(
\bar r_{\rm hist}+\bar r_{\rm dwell}
+\bar r_{\rm HRV}+\bar r_{\rm phase}
+\bar r_{\rm sample}+\bar r_{\rm model}
+\bar r_{\rm measurement}
\right)\\
&+\epsilon^{\rm omitted}_{r,m}
+\epsilon^{\rm delay}_{r,m}
+\epsilon^{\rm filt-init}_{r,m}
+\epsilon^{\rm quad}_{r,m}.
\end{aligned}
\]

总和是保证覆盖的外圆盘，不是相关误差下的最小圆盘。每个残差项的总预算必须满足 \(\bar r\ge\sum_j\bar r_j\)，不能写反。

### 5.1 机器传播合同

模板中的 `computed_*` 数值必须是在冻结设计所允许的全部窗口上统一有效的最坏上界，不是某一目标窗口上的事后估计。`evidence_policy.bound_numeric_scope` 因而固定为
`uniform_all_windows_satisfying_frozen_design`。

每个窗口--谐波结果必须保存 Gram 最小特征值的认证下界、条件数的认证上界、`gram_bounds_outward_certified=true` 和证书哈希。令

\[
\overline a_m\ge\sqrt{(G^{-1})_{mm}},\qquad
\underline H_m\le|\widehat H_m|,\qquad
\kappa_m=\frac{\overline a_m}{\underline H_m}.
\]

只知道认证下界 \(\underline\lambda_{\min}(G)>0\) 时，可保守取
\(\overline a_m=1/\sqrt{\underline\lambda_{\min}(G)}\)。验证器按下表单向传播；报告值只能大于或等于向上舍入后的重算下界。

| 结果字段 | 冻结证据到半径的最低传播 |
|---|---|
| `radius_history` | \(\kappa_m\,\bar r_{\rm hist}\) |
| `radius_dwell` | \(\kappa_m\,\bar r_{\rm dwell}\) |
| `radius_hrv` | \(\kappa_m\,\bar r_{\rm HRV}\) |
| `radius_phase_anchor` | \(\kappa_m\,\bar r_{\rm phase}\) |
| `radius_sampling` | \(\kappa_m\,\bar r_{\rm sample}\) |
| `radius_measurement_noise` | \(\kappa_m\,\bar r_{\rm measurement}\) |
| `radius_model_residual` | \(\kappa_m\,\bar r_{\rm model}\) |
| `radius_window_leakage` | `computed_coefficient_radius` |
| `radius_quadrature` | `quadrature_radius` |
| `radius_delay` | `computed_coefficient_radius` |
| `radius_filter_initialization` | `computed_coefficient_radius` |

测量噪声、输出端模型余项、采样误差和相对延迟还要从各自的原始校准字段重新得到 `computed_*` 的下界，防止把非零原始误差隐藏为零。逐项和、总半径与 union-bound 失败概率均使用单向向上舍入；对称的 `isclose` 不能签发覆盖证书。

## 6. 各误差项的预注册合同

### 6.1 finite-lower / pre-record history

记录起点为 0，当前心率段从无量纲时刻 \(S_0\) 开始。令

\[
M_{\alpha,\lambda}(s)=E_\alpha(-\lambda s^\alpha).
\]

新名义输入满足 \(\|u_0\|_\infty\le U_0\) 时，

\[
B_{\rm hist}(s)=
\left(|z_0|+\frac{U_0}{\lambda}\right)
M_{\alpha,\lambda}(S_0+s).
\]

若真实过程还有未由有限下限初值表示的记录前 history，必须另给 \(C_{\rm pre}\) 并加入 \(C_{\rm pre}M_{\alpha,\lambda}(S_0+s)\)。缺少适用的 \(z_0,U_0,C_{\rm pre}\) 独立上界时为 NOT_CERTIFIABLE。

必填：记录起点、心率段起点、\(|z_0|\)、\(U_0\)、\(C_{\rm pre}\)、\((\alpha,\lambda)\) 区间、来源及最终 \(\bar r_{\rm hist}\)。

### 6.2 previous-rate memory / finite dwell

前段输入满足 \(\|u_{\rm pre}\|_\infty\le U_{\rm pre}\) 时，

\[
B_{\rm dwell}(s)=
\frac{U_{\rm pre}+U_0}{\lambda}
\left[M_{\alpha,\lambda}(s)
-M_{\alpha,\lambda}(S_0+s)\right].
\]

对 \(j\in\{\mathrm{hist},\mathrm{dwell}\}\)，

\[
\bar r_j\ge
\sup_{(\alpha,\lambda)\in I\times\Lambda}
\left(\sum_nw_nB_j(s_n)^2\right)^{1/2}.
\]

窗口必须给出从段起点到估计窗左端的 burn-in。不能因“看起来平稳”把 dwell 设为零。

### 6.3 HRV 与操作性相位扭曲

WLS 基始终使用 3.1 节由窗口两端确定的名义恒频相位。另用逐搏 R 峰线性插值得到操作性相位，只把它相对名义恒频相位的偏差作为 HRV 扰动传播；不声称不可观测的细胞振荡器相位。该偏差由时间戳计算，但 R 峰检测误差仍需独立校准。若

\[
|R_j-\widehat R_j|\le\Delta_R,
\]

则第 \(j\) 搏内

\[
\Delta_{\phi,R,j}\le
\frac{2\pi\Delta_R}{\widehat T_j-2\Delta_R},
\qquad \widehat T_j>2\Delta_R.
\]

若仍声称潜在生理相位，必须另给 \(\Delta_{\rm latent}\)；R 峰本身不能证明该界。定义

\[
d(x)=2\sin\left(\frac{\min\{x,\pi\}}2\right),\qquad
\Delta U(s)=\sum_{k\in\mathbb Z}\overline Q_k
d(|k|\Delta_\phi(s)).
\]

HRV 经 Caputo 动力学传播后的包络为

\[
B_{\rm HRV}(s)=
\int_0^s
\sigma^{\alpha-1}E_{\alpha,\alpha}(-\lambda\sigma^\alpha)
\Delta U(s-\sigma)\,d\sigma.
\]

若 \(\Delta U\le\overline{\Delta U}\)，则

\[
B_{\rm HRV}(s)\le
\frac{\overline{\Delta U}}{\lambda}
[1-E_\alpha(-\lambda s^\alpha)].
\]

延长 dwell 会减小前史项，却可能增大单恒频近似的 HRV 项；不得预设总半径随窗长单调下降。

### 6.4 R 峰 anchor 的剩余相位误差

若估计基与公共 gauge 间仍有非共同偏差 \(\Delta_{\rm anchor}(t)\)，令 \(B_{\rm phase}\) 表示该误差在固定前端输出端的独立包络。identity 前端下，一个充分界为

\[
B_{\rm phase}(t)\le
\sum_{k\in\mathcal M}\overline Z_k
d(|k|\Delta_{\rm anchor}(t)).
\]

非 identity 前端必须把相应的传递幅值包络纳入 \(B_{\rm phase}\)。最终登记

\[
\bar r_{\rm phase}\ge
\left(\sum_nw_nB_{\rm phase}(t_n)^2\right)^{1/2}.
\]

它与 6.3 的动力学 HRV 项必须定义为不重叠来源；若重复覆盖，应明确是保守外包络。

### 6.5 非整数窗、其他谐波与谐波尾

一般 WLS 的 alias/leakage 因子为

\[
h_{mk}=e_m^*G^{-1}\Psi^*W\boldsymbol\psi_k,
\]

\[
\epsilon^{\rm omitted}_{r,m}\le
\sum_{k\notin\mathcal M}
\frac{|h_{mk}|\,\overline{|H_{r,k}|}}
{|\widehat H_{r,m}|}\overline Z_{r,k}.
\]

必须预注册无限谐波尾的绝对可和包络：

\[
\sum_{k\in\mathbb Z}
\overline{|H_{r,k}|}\overline Z_{r,k}<\infty.
\]

对单谐波连续矩形 lock-in，\(L=\widehat\omega_rT\) 时

\[
h_{mk}=e^{i(k-m)L/2}
\operatorname{sinc}\left(\frac{(k-m)L}{2}\right).
\]

\(L=2\pi N\) 只消除名义恒频谐波的窗泄漏，不消除 HRV、历史、噪声或模型误差。无谐波尾上界时为 NOT_CERTIFIABLE。

### 6.6 采样、反混叠、时间戳和求积

精确采样的有限谐波模型由离散 WLS 直接恢复，不另产生求积误差。ADC 量化、采样抖动、R 峰端点插值、丢样填补和抗混叠余项必须作为样本扰动：

\[
\bar r_{\rm sample}\ge
\left\{\sum_nw_n
\left[M_1\Delta t_n+\frac{\Delta_{\rm ADC}}2
+\epsilon_{{\rm interp},n}\right]^2
\right\}^{1/2}.
\]

若实现采用梯形公式逼近连续 lock-in，还需

\[
\epsilon^{\rm quad}_{r,m}\le
\frac{1}{12T|\widehat H_{r,m}|}
\sum_j\Delta t_j^3M_{2,j},
\]

\[
M_{2,j}\ge
\sup_{[t_j,t_{j+1}]}
|x''-2im\widehat\omega_rx'-(m\widehat\omega_r)^2x|.
\]

普通样本有限差分不是 \(M_1,M_2\) 的确定性证明。没有硬件抗混叠规格、带宽/导数上界或区间重建证书时为 NOT_CERTIFIABLE。

### 6.7 前端幅值和相对群延迟

令

\[
\frac{H_{r,m}}{\widehat H_{r,m}}
=(1+\delta_a)e^{-i\delta_\psi},\quad
|\delta_a|\le a_H,
\]

\[
|\delta_\psi|\le
|m|\widehat\omega_r\Delta_\tau+\Delta_{\psi,0}.
\]

则

\[
\epsilon^{\rm delay}_{r,m}\le
\overline Z_{r,m}\left[
a_H+(1+a_H)
d(|m|\widehat\omega_r\Delta_\tau+\Delta_{\psi,0})
\right].
\]

必须校准 ECG、R 峰标注通道和数字滤波链之间的相对延迟，不只是软件滤波器的理论群延迟。

### 6.8 因果滤波器初始化

若稳定因果滤波器 \(h\) 在 \(t_0\) 以零状态启动，此前输入满足 \(|y|\le Y_{\max}\)，

\[
\epsilon^{\rm filt-init}_{r,m}\le
\frac{Y_{\max}}{T|\widehat H_{r,m}|}
\int_a^b\int_{t-t_0}^{\infty}|h(s)|\,ds\,dt.
\]

若恢复滤波器状态，必须记录来源和完整性；不能默认 IIR 状态在每窗开始时稳态。

### 6.9 测量噪声

\[
\|\boldsymbol v\|_W\le\bar v_W,\qquad
\epsilon^{\rm noise}_{r,m}\le
\frac{a_m\bar v_W}{|\widehat H_{r,m}|}.
\]

当 \(w_n\) 是归一化时间权且独立协议给出

\[
\|v\|_{L_2(a,b)}\le N_2,
\]

可取 \(\bar v_W\le N_2/\sqrt T\)。只有点态界时，最坏系数误差通常不随 \(T\) 下降。概率噪声只能进入 PASS_PROBABILISTIC。

### 6.10 声明代理的模型余项

非标量 ECG 成分、基线变化、未建模非平稳形态和前端非线性记为 \(r_{\rm model}\)。直接输出余项必须界定 \(\|r_{\rm model}\|_W\)；只有动力学输入扰动 \(d\) 才允许经核传播，例如

\[
|d|\le\bar d
\Longrightarrow
|r_{\rm model}(s)|\le
\frac{\bar d}{\lambda}[1-E_\alpha(-\lambda s^\alpha)].
\]

模型余项必须来自独立校准集、工程规格或冻结的外部包络，不得由本窗口同一代理的拟合残差反推。

## 7. 最小驻留和最小 SNR 门

### 7.1 证书意义下的最小 burn-in

\[
C_H(\alpha,\lambda)=
\max\left\{
|z_0|+\frac{U_0}{\lambda},
\frac{U_{\rm pre}+U_0}{\lambda}
\right\}.
\]

对固定窗口设计，除 history/dwell 外的半径为 \(\epsilon_{\rm rest}\)，总半径门为 \(\epsilon_\star\)。定义

\[
B_{\min}=\inf\left\{
B\ge0:
\frac{a_m}{|\widehat H_{r,m}|}
\sup_{(\alpha,\lambda)\in I\times\Lambda}
C_H(\alpha,\lambda)E_\alpha(-\lambda B^\alpha)
\le\epsilon_\star-\epsilon_{\rm rest}
\right\}.
\]

只有 \(\epsilon_\star>\epsilon_{\rm rest}\) 且 \(B_{\min}<\infty\) 才可通过。物理段驻留必须满足

\[
T_{\rm dwell}\ge
\tau_\star B_{\min}+\frac{2\pi N}{\widehat\omega_r}.
\]

这是外包络证书的最小驻留，不是心脏达到稳态的真实最短时间。

### 7.2 最小确定性能量 SNR

独立先验给出 \(|Z_{r,m}|\ge\underline Z_{r,m}>0\) 时，定义

\[
\mathrm{SNR}_{E,r,m}=
\frac{|\widehat H_{r,m}|\underline Z_{r,m}\sqrt T}{N_2}.
\]

目标相对半径为 \(\eta\)，且
\(\eta\underline Z_{r,m}>\epsilon_{\rm rest}\) 时，充分条件为

\[
\mathrm{SNR}_{E,r,m}\ge
\frac{a_m}
{\eta-\epsilon_{\rm rest}/\underline Z_{r,m}}.
\]

没有独立的 \(\underline Z_{r,m}\) 或确定性 \(N_2\) 时，不得报告确定性最低 SNR 已通过。

## 8. 决策顺序

每个窗口、谐波按顺序判定：

1. 协议在读取目标 ECG 或拟合结果前冻结，否则 NOT_CERTIFIABLE。
2. 相同导联/增益、共同 gauge、因果前端、无逐窗归一化；明确违反为 EXCLUDE_PROTOCOL。
3. \(\tau_\star\)、信号单位及全部项量纲闭合，否则 NOT_CERTIFIABLE。
4. \(G\succ0\)、条件数合格、\(\widehat H_{r,m}\ne0\)、窗口和权重可复现，否则 NOT_CERTIFIABLE。
5. 所有适用项 provided=true、来源可追踪且未使用目标拟合残差，否则 NOT_CERTIFIABLE。
6. 确定性请求中只要一个必需项是概率性，即 NOT_CERTIFIABLE；可另报 PASS_PROBABILISTIC。
7. 谐波尾、前史、抗混叠/导数、延迟和总半径必须有限，否则 NOT_CERTIFIABLE。
8. 总半径、相对半径、最小驻留和最低 SNR 同时满足冻结门，才为 PASS_DETERMINISTIC。
9. 概率协议还必须给出全部选定窗口--心率--谐波圆盘的同时覆盖概率，才为 PASS_PROBABILISTIC。
10. 只有通过的圆盘进入 R1；R1 联合集为空才可能 REJECT_MODEL。概率圆盘只能给出置信度限定的拒绝；非空但过宽仍是 NOT_CERTIFIABLE。

## 9. 预注册、冻结和审计产物

冻结前必须填写：

- 协议 ID、版本、冻结时间、Git commit 和协议 SHA-256；该哈希只覆盖预注册部分并明确排除运行后的 `window_result`；
- 数据角色与受试者级拆分 ID，不把原始 ECG 复制进仓库；
- \(\tau_\star,I,\Lambda,\mathcal M\)、信号单位和信号表示；
- R 峰算法/版本及时间误差校准；
- 前端因果性、系数、状态策略、频响和延迟校准；
- 每项的来源、单位、证据语义、独立性、数值和适用范围；
- 总/相对半径、Gram、最小驻留和最低 SNR 门；
- 概率圆盘的逐项 \(\delta_j\) 与联合覆盖方法。

运行后只允许填写：窗口/心率 ID、受试者、导联与物理增益 ID，R 峰端点、\(N,T,\widehat\omega_r\)、外向认证的 Gram 下/上界及证书哈希、\(\overline a_m\)、\(\widetilde Z\)、冻结公式传播出的逐项/总半径、状态和原因代码。同一 R1 联合组必须保持受试者、导联和物理增益一致；结果键 \((\texttt{window\_id},m)\) 必须唯一。每一结果行必须保存所引用的 `protocol_sha256`，并对排除自身 `result_sha256` 字段后的结果内容生成 `result_sha256`。修改运行结果不得改变冻结的协议哈希。

运行后不得改上界来源、公式、阈值、谐波集合或排除规则。修订必须生成新版本并保留旧结果。

原始 ECG 继续只读地从 ECG_DATA_DIR 解析；派生数据写 ECG_PROCESSED_DIR，结果写 ECG_OUTPUT_DIR。Git 只保存协议、配置、公式、代码和小型汇总。任何单个或计划批量下载超过 100 MB、或大小未知且可能超过 100 MB 的数据/软件，必须先获用户批准。

## 10. R2 关闭条件

文档和模板完成不等于 R2 通过。只有以下条件全部满足才可关闭：

1. 每个必需项有独立、量纲一致、可审计的数值上界；
2. 覆盖性由解析不等式或有保证区间证书证明；有限合成压力试验只作实现回归，不替代对全部有界扰动的全称证明；
3. R 峰时序、抗混叠/采样、群延迟、前史、谐波尾、噪声和模型余项全部关闭；
4. 最小驻留和最低 SNR 门在冻结协议下可满足；
5. 没有用同一拟合残差反向设定半径；
6. 确定性与概率性结论在文件名、表格和论文措辞中明确分离。

在此之前，真实 ECG 不能进入确定性阶次认证。

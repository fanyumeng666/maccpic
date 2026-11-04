# ccp1d 图册与数据索引（学习用）

所有图都由 `tools/` 下的脚本从代码的 ASCII 输出画出，数据文件在 `results/`，重画命令写在每节末尾。按"看什么—怎么解释—对应面试题"组织；面试题编号对应 `program_flow.md` 第 5 节的速答卡。

## 0. 建议的学习顺序

1. `doc/fig/program_flow.png` + `doc/program_flow.md`：先把一步循环（push → mcc → deposit → poisson → diag）和每个模块的算法讲清楚。
2. 第 2 节初始条件：装载、噪声、尺度（λ_D/dx、ωpe·dt、ν_max·dt）——面试最常问的数值约束都在这里。
3. 第 3 节无碰撞鞘层：最简单的物理检验，能说明为什么需要碰撞。
4. 第 4、4b 节 Case 1：基准对比、剖面、能谱、x–t 图、功率去向、三把能量尺；收敛图说明"误差有多大"。
5. 第 5–7 节扩展：Case 2 的压强依赖、双频、电子群对 Boltzmann 方程——展示代码不只是复现一个数。
6. 第 8 节测试表：每个数字都能在 `PROGRESS.md` 里找到出处。

## 1. 程序结构

**`doc/fig/program_flow.png`**（`tools/draw_flow.py`）：启动 → 每步 push / mcc / deposit / poisson / diag → 每周期 history 与续算 → 结束输出与分析工具；绿色是已完成并通过测试的部分。配 `doc/program_flow.md` 逐模块讲解。

## 2. 初始条件（`doc/fig/init_*.png`，`tools/plot_initial.py`，数据 `results/init_case1/`）

| 图 | 看什么 | 要点 |
| --- | --- | --- |
| `init_phase_space.png` | t = 0 的 x–v_x 相空间 | 均匀随机位置、Maxwell 速度；电子与离子热速度差 855 倍 |
| `init_velocities.png` | v_x 直方图（Poisson 误差棒）、电子能量分布、折叠的 ± 尾巴 | 3 v_th 以外每侧只有 88 ± 9 个粒子，两侧看起来"不对称"是统计涨落；T_e = 2.59 eV 时超过电离阈值 24.59 eV 的电子只有约 3×10⁻⁴ |
| `init_grid.png` | t = 0 的 n_e、n_i、ρ、φ、E | ρ 是 5 % 的散粒噪声，φ 却有 ±14 V 的长波结构：泊松方程把长波放大 1/k² |
| `init_drive.png` | 驱动电压与所有尺度 | λ_D/dx = 1.43，ωpe·dt = 0.166，离子渡越 58 个周期 |
| `init_xsec.png` | 6 张截面、碰撞频率、ν_max | ν_max = 8.65×10⁷ s⁻¹ 在 8.2 eV，P = 1.6 %/步（题 1、4） |

讲解：`doc/initial_conditions.md`。重画：`cd runs/init_case1 && mpirun -np 2 ../../ccp1d && python3 ../../tools/plot_initial.py .`（输入 `inputs/init_case1.in`）。

## 3. 无碰撞鞘层（`doc/fig/sheath_collisionless.png`，`tools/check_sheath.py`，数据 `results/sheath/`）

Case 1 的等离子体、无气体、两壁接地、10 µs。能量守恒到 5×10⁻⁴，Γ_e = Γ_i 到 1.6 %。最值得记的现象：一维无碰撞时电子的 T_x 从 2.57 eV 塌到 0.13 eV 而 T_⊥ 不变，浮动电势只有 0.19 V——没有碰撞把 x 和 ⊥ 方向的能量混合，鞘层电势撑不住；这正是 Case 1 需要碰撞和射频加热的原因（题 16）。

## 4. Case 1 基准（`results/case1/`，4 核 170 s）

| 图 | 脚本 | 看什么 |
| --- | --- | --- |
| `compare_turner.png` | `compare_turner.py` | 时间平均 n_i(x) 与 Turner 五个代码的参考曲线叠画 + 相对误差：中心 +2.1 %，129 点 rms 1.5 %，电极处 3.62e13 vs 3.615e13（题 12） |
| `profiles.png` | `plot_profiles.py` | n_e/n_i、⟨φ⟩ = 223 V 的平台、⟨E⟩、电子能量三条线（⟨ε⟩、T_x、T_eff）、电离率、加热剖面（精确做功 vs 网格 J·E）、离子通量、离子能量 |
| `dist.png` | `plot_dist.py` | IEDF 的电荷交换多峰结构（最大 240 eV ≈ V₀/2）、壁面电子能谱、EEPF 双温（5.2 与 11.0 eV，题 7、8） |
| `xt.png` | `plot_xt.py` | 8 张 x–t 相位图：鞘层边缘 0.26–2.72 cm 摆动，电子加热紧贴扩张中的鞘层边缘（随机加热，题 5、6），电离跟在后面 |
| `diag_checks.png` | `check_diag.py` | 鞘层边缘随相位（左 0.26–2.72 cm，电极最负时最宽）、总电流 J_e + J_i + ε₀∂E/∂t 沿 x 均匀到 1.7 %（振幅 9.2 A/m²，Turner 设计四个算例电流幅值都约 10 A/m²）、EEPF 两段斜率 5.2 / 11.0 eV、IEDF 平均 78.5 eV 对比 e⟨φ_c⟩ = 223 eV |
| `power.png` | `plot_power.py` | 125 W/m² 的去向：离子撞壁 28 %、电子撞壁 8 %、激发+电离 19 %、离子→气体 45 %（题 13） |
| `history.png` | `plot_history.py` | 粒子数先降后平（稳态电子只剩初始的 19 %、离子 30 %：HANDOVER 里"涨到 3–5 倍"的估计是错的，PLAN §1）、密度、通量、能量、功率、守恒量 |
| `convergence.png` | `convergence.py` | 基线、换种子、dt/2、nx 256、nppc 256/1024 的剖面与偏差；换种子噪声地板 max 1.4 % / rms 0.6 %（题 11） |

数据：`profiles.dat`（14 列时间平均剖面，列名在文件头）、`history.dat`（31 列每周期）、`iedf_*.dat`、`eedf_*.dat`、`eepf.dat`、`info.txt`（参数、导出量、能量平衡三把尺）。x–t 矩阵（`xt_*.dat`）不进 git，重跑 `inputs/case1.in` 三分钟即得。收敛算例的小文件在 `results/conv/`。

数值要点（题 4、11）：稳态 ωpe·dt = 0.123、λ_D/dx = 3.59（Turner 表 III：0.121、3.72）；网格乘积 ⟨J_e·E⟩ 只有电子精确做功的 54 %，dt 减半后 77 %，dx 减半不变——是电子一步走 v·dt 的中心化误差，所以加热剖面用精确做功沉积。

## 4b. 与 Turner 表 III 的对照（`tools/compare_table3.py`，数据 `data/turner2013_table3.dat`）

Turner 等 2013 的表 III 给了四个算例的中面离子密度、电子温度、电子/离子吸收功率、电极离子电流和数值品质因数，是 Case 1 以外唯一的定量参考。`python3 tools/compare_table3.py <run_dir>` 自动识别算例并打印对照（T18 判据：n_i ±5 %，S_e、S_i ±10 %，J_i ±5 %）。Case 1：

| 量 | ccp1d | Turner | 偏差 |
| --- | --- | --- | --- |
| n_i 中面 [m⁻³] | 1.432e14 | 1.40e14 | +2.3 % |
| kT_e 中面 [eV]：2/3⟨ε⟩ / T_x | 8.88 / 9.79 | 9.36 | −5.1 % / +4.6 % |
| S_e 电子功率 [W/m²]（精确做功） | 34.05 | 34.3 | −0.7 % |
| S_i 离子功率 [W/m²] | 90.73 | 90.6 | +0.1 % |
| J_i 电极离子电流 [A/m²] | 0.219 | 0.219 | −0.1 % |
| ωpe·dt / λ_D/dx | 0.123 / 3.59 | 0.121 / 3.72 | +1.3 % / −3.5 % |
| ν_e,max·dt | 0.0158 | 0.0158 | +0.2 % |
| N_D 每德拜长度粒子数 | 998 | 1042 | −4 % |
| N_P 总粒子数（电子+离子） | 3.16e4 | 3.19e4 | −0.9 % |

要点：积分量（功率、电流）在 1 % 内，密度 +2 %；Turner 的 S_e = 34.3 W/m² 对应我们的**精确做功**而不是网格乘积 ⟨J_e·E⟩（18.3），说明他们的 J·E 是时间中心化的——面试里被问"J·E 为什么对不上"时这是最好的证据（题 4）。电子温度相差 5 %，但论文没有给出 kT_e 的定义（2/3 平均能量还是 x 方向温度），我们的两个定义正好把 9.36 eV 夹在中间，其他积分量都吻合，所以更可能是定义差异。

## 5. Case 2 试跑：压强依赖（`results/case2/`）

（运行中，结果见下文更新。）

## 6. 双频扩展（`results/dual/`）

| 图 | 看什么 |
| --- | --- |
| `results/dual/dual_frequency.png`（`dual_frequency.py`） | V_lf = 0/50/100/150 V 的 IEDF、最大能量 240 → 298 eV（≈ 0.5(V₀+V_lf)）、宽度、密度 −9/−26/−45 %、功率 |
| `results/dual/dual_100_lf/xt_lf_resolved.png`（`plot_xt.py`，以低频为主周期） | 一个低频周期里 10 次高频鞘层振荡骑在低频包络上，电极最正时鞘层塌到电极 |

物理结论（题 10）：低频抬高最大离子能量，但在 6.7 cm、30 mTorr 的几何里鞘层按 s ∝ V^{3/4} 变厚，200 V 时两侧鞘层相遇、放电熄灭（`dual_200/`、`dual_300/` 的 history 是证据）；工业双频反应器间隙小、密度高、鞘层毫米级，解耦才成立。

## 7. 电子群与两项 Boltzmann 方程（`results/swarm/`）

**`swarm.png`**（`compare_swarm.py`）：PIC 电子群（点）与两项 Boltzmann 求解器（线，`boltzmann2term.py`，同一套截面）的平均能量、漂移速度、α/N 随 E/N，以及 7 个 E/N 的能量分布。20–200 Td 差 0.1–1 %。两个教训：两项方程的动量转移截面要包含所有各向同性过程；高 E/N 下氦的高能电子逃逸，要用脉冲 Townsend 模型（题 20）。数据：`swarm_pic.dat`、`swarm_boltzmann.dat`（含 `hi/` 的 300、500 Td）、`swarm_eedf_*Td.dat`、`swarm_b2t_eedf_*Td.dat`。重跑：`make check-swarm`（8 分钟）。

## 8. 验证与测试结果（数字备查）

| 测试 | 结果 |
| --- | --- |
| T1 泊松 | 解析解 1e-14，二阶收敛比 4.0004，离散高斯定律 1e-13 |
| T2 沉积 | 权重、电荷守恒、端点半格、矩、1 与 2 rank 一致 |
| T3 随机数 | 矩、状态存取逐位；`check-rng`：2×10⁷ 粒子的尾部计数对称到 ±2.6σ 内，1000 个种子的均值/偏度/尾部 z 值 rms = 1.00 |
| T4 截面 | 查表 = 原表（电子 1e-15，离子 1e-7）；ν_max 8.654e7 / 2.793e7 |
| T5–T7 MCC | 碰撞计数 = n_gσv × 0.992，\|z\| ≤ 1.7；弹性公式 1e-15、各向同性 χ² = 12/20/18（19 自由度）、新离子 T_g；10 eV 离子 20 µs 热化到 300 K 偏差 0.05 %，电子弹性冷却与 ODE 差 0.02 % |
| T8 续算 | 20 周期 vs 10+10 逐位相同（含碰撞） |
| T10 Langmuir | 频率与蛙跳理论差 3e-4 |
| T11 鞘层 | 能量 5e-4、Γ_e = Γ_i 1.6 %、T_⊥ 不变、对称 0.14 % |
| T12 电子群 | 20–200 Td 平均能量、漂移速度差 0.1–1 %，α/N 2–5 % |
| T13/T14 Case 1 | 13 项自检全过；中心 +2.1 %，rms 1.5 % |
| T15 收敛 | 加密算例与基线差 ≤ 2 倍换种子噪声 |
| T16 自洽 | 23 项全过；总电流沿 x 均匀 1.7 % |
| T17 双频 | E_max 单调上升 ≈ 0.5(V₀+V_lf)，宽度净增 4.5σ |
| T18 表 III | Case 1 的 S_e、S_i、J_i 与 Turner 差 ≤ 0.7 %，n_i +2.3 % |

命令：`make check`（2 分钟）、`make check-rng`、`make check-swarm`；每个算例目录 `python3 tools/check_run.py <dir> --full` 和 `python3 tools/check_diag.py <dir>`。

# ccp1d 学习计划（面试准备版）

目标：两周内能不看稿讲清楚这个代码做什么、每一步用什么方法、怎么证明它是对的、结果说明了什么物理；并能回答追问。每天 2–3 小时；也可以压缩成 5 个整天。每一阶段都有"读、跑、算、答"四件事，最后的"自测"过了再进下一阶段。

## 0. 文档地图（先知道有什么）

| 文档 | 内容 | 什么时候读 |
| --- | --- | --- |
| `README.md`（ccp1d） | 现状、编译、测试、运行、文件表 | 第 1 天 |
| `doc/program_flow.md` + `doc/fig/program_flow.png` | 逐模块讲解：做什么、算法、为什么、追问；完成状态；20 题速答卡 | 第 1 天通读，之后反复 |
| `doc/code_guide.md` | 逐文件、逐例程的说明；输出格式；测试与工具清单；改什么去哪里 | 第 2–5 天对照源码读 |
| `doc/initial_conditions.md` + `doc/fig/init_*.png` | 初始条件讲解：装载、噪声、泊松放大长波、尺度、截面与零碰撞、无碰撞鞘层 | 第 2–3 天 |
| `doc/figures.md` | 图册与数据索引：每张图看什么、关键数字、重画命令、测试数字表 | 第 5–7 天 |
| `doc/theory.md` + `results/case*/theory.txt` | 与经典公式的对照和解释 | 第 7 天 |
| `PLAN.md` | 实施计划、对交接文件的 15 处修正、测试矩阵 T1–T18 | 第 6 天 |
| `PROGRESS.md` | 十次会话的日志：每一步做了什么、遇到什么问题、怎么解决 | 随时查数字出处 |
| `../HANDOVER.md` | 原始规格 | 想知道"任务是什么"时 |

## 1. 之前讲解过的内容在哪里

| 对话里讲过的 | 写进了 |
| --- | --- |
| 修改计划与测试方案 | `PLAN.md` |
| 初始条件五张图的讲解 | `doc/initial_conditions.md` §1–4 |
| "+v 高速粒子更多"是统计涨落，怎么量化 | `doc/initial_conditions.md`"尾巴看起来不对称怎么办"；`tests/random/test_tails.f90`、`test_seeds.f90`；`make check-rng` |
| 功能流程图与每个模块的作用算法 | `doc/program_flow.md`（全文）、`doc/fig/program_flow.png` |
| 网格 ⟨J_e·E⟩ 只有做功的 54 %：是网格太粗吗？PIC 的空间/时间步长限制？质量比用多少？ | `doc/program_flow.md` diag 一节 + 速答卡 17、18；`PROGRESS.md` 2026-11-03（dt/2、dx/2 的定量验证） |
| M2 碰撞：截面、零碰撞法、电子/离子碰后速度 | `doc/program_flow.md` mcc 一节；`doc/code_guide.md` §7–8；`doc/initial_conditions.md` 截面图 |
| M5 收敛：噪声地板、dt/2、nx 256、nppc | `doc/figures.md` §4 convergence；`PROGRESS.md` 2026-11-03 |
| M6 诊断：x–t 图、总电流检验、功率去向、三把能量尺 | `doc/program_flow.md` diag 一节；`doc/figures.md` §4 |
| M7 双频：为什么解耦在这里不完全成立 | `doc/figures.md` §6；速答卡 19 |
| M3 电子群 vs 两项 Boltzmann：两个坑 | `doc/figures.md` §7；速答卡 20；`tools/boltzmann2term.py` 文件头 |
| Case 2 试跑与压强依赖 | `doc/figures.md` §5；`results/case_comparison.png` |
| Turner 表 III 对照 | `doc/figures.md` §4b；`tools/compare_table3.py` |
| 数据能和理论合上吗 | `doc/theory.md` |
| 随机数是自己写的还是库函数 | `doc/program_flow.md` random 一节；`doc/code_guide.md` §3 |

## 2. 预备知识自查（第 0 天，半天）

能默写并说出物理意义：Debye 长度 λ_D = sqrt(ε₀kT_e/ne²)、电子等离子体频率 ω_pe、Bohm 速度 u_B = sqrt(eT_e/M)、浮动电势、鞘层与预鞘、Child 定律 J ∝ V^{3/2}/s²、Maxwell 分布与 EEPF 的半对数直线、碰撞频率 ν = n_gσv 与自由程 1/(n_gσ)、Townsend 系数 α。数值方面：中心差分、三对角方程、蛙跳格式的二阶与辛性质、蒙特卡洛抽样（反函数法、接受-拒绝）、随机噪声 ∝ 1/√N。不会的查 Lieberman & Lichtenberg（§2、§3、§6）和 Birdsall & Langdon（§2–4）。

## 3. 分阶段计划

### 第 1 天：全貌

- 读：`README.md`；`program_flow.md` §0–2（不求全懂）；`figures.md` §0–1。
- 跑：`./setup.sh`（新容器）→ `cd ccp1d && make && make check`（2 分钟，应 ALL TESTS PASSED）→ 测试模式 Case 1：`mkdir -p runs/t && cp inputs/case1_test.in runs/t/ccp1d.in && cd runs/t && mpirun -np 2 ../../ccp1d`（8 秒）→ `python3 ../../tools/plot_history.py . && python3 ../../tools/check_run.py .`。
- 算：用 info.txt 里的 n0、Te0 手算 ω_pe、λ_D、ωpe·dt、λ_D/dx，和 info.txt 对上。
- 答：一步循环的四件事是什么？粒子和网格各自存什么？为什么每个进程都持有整个网格？

### 第 2 天：场与沉积（poisson、deposit）

- 读：`code_guide.md` §4–5；`poisson.f90`、`deposit.f90` 全文（166 行）；`tests/poisson`、`tests/deposit`；`initial_conditions.md` 关于 ρ、φ 的段落。
- 跑：`./tests/poisson/test_poisson`；把 `test_poisson.f90` 里的 nx 改一改看二阶收敛比。
- 算：在纸上写出 Thomas 算法的前扫和回代；证明面场 E_{j+1/2} 满足离散高斯定律；证明 CIC 沉积的电荷守恒 Σρ_j dx = Σqw。
- 答：为什么 CIC 而不是 NGP？沉积和插值用同一形函数有什么好处（无自力、动量守恒）？端点为什么乘 2？泊松方程为什么把噪声的长波放大（1/k²）？

### 第 3 天：推进、装载、无碰撞鞘层（push、init、random）

- 读：`code_guide.md` §3、§6、§10；`push.f90`、`init.f90`、`random.f90`；`initial_conditions.md` 全文；`tests/langmuir`。
- 跑：`./tests/langmuir/langmuir`（看频率与 (2/dt)asin(ω_pe dt/2) 的差）；无碰撞鞘层 `inputs/sheath.in`（2 核 21 秒）+ `tools/check_sheath.py`；初始条件图 `inputs/init_case1.in` + `tools/plot_initial.py`。
- 算：推导做功恒等式 ½m(v⁺² − v⁻²) = qE_p dt·½(v⁺ + v⁻)；推导蛙跳的数值频率；算一下 3 v_th 以外应有多少粒子（N/2·erfc(3/√2)）。
- 答：蛙跳为什么要回推半步？被吸收的粒子怎么删（O(1) 交换）？为什么 T_⊥ 在无碰撞时必须不变、而 T_x 会塌？随机数为什么自己写、怎么做到续算逐位一致？

### 第 4 天：碰撞（xsec、mcc）

- 读：`code_guide.md` §7–8；`xsec.f90`、`mcc.f90`；`program_flow.md` mcc 一节；`tests/mcc/test_mcc.f90` 的头注释；`data/SOURCES.md`。
- 跑：`./tests/mcc/test_xsec`、`./tests/mcc/test_mcc`；看 `doc/fig/init_xsec.png`。
- 算：从截面表手算 8.2 eV 处的 ν_max 和 P = 1 − e^{−ν_max dt}；推导弹性碰撞的能量损失公式 ε' = ε[1 − (2m/M)(1 − cosχ)]；证明背向散射 v' = v_cm − ½g = v_n。
- 答：零碰撞法为什么严格等价？一步最多一次碰撞带来多大偏差（0.8 %）、测试里怎么修正？电离后能量为什么平分？离子截面为什么用质心系能量？ν_max 为什么要能量上限？

### 第 5 天：诊断与能量（diag）

- 读：`code_guide.md` §9、§13；`diag.f90`（678 行，按例程表的顺序读）；`program_flow.md` diag 一节；`figures.md` §4。
- 跑：全尺寸 Case 1（4 核 3 分钟）`inputs/case1.in`；`plot_profiles`、`plot_dist`、`plot_xt`、`check_diag`、`plot_power`。
- 算：从 info.txt 验证 closure_residual 的公式；用 history 验证 KE − KE(0) = Work_cum − E_wall_cum；用 profiles 对 P_work_e 做梯形积分和 P_part_e 比。
- 答：三把能量尺各查什么？为什么网格 ⟨J_e·E⟩ 只有 54 %、为什么是一阶 dt 误差、为什么不影响密度？x–t 图上怎么看出随机加热？总电流为什么必须处处相同？续算文件为什么分成粒子类和网格类？

### 第 6 天：基准、收敛、表 III

- 读：`PLAN.md` §1（对交接文件的修正）、§4（测试矩阵）；`figures.md` §4、§4b；`tools/compare_turner.py`、`convergence.py`、`compare_table3.py`。
- 跑：`compare_turner.py runs/case1`；`compare_table3.py runs/case1`；有时间再跑一个换种子（seed = 2）看噪声地板。
- 算：中心密度的统计噪声估计（每格 280 粒子、12800 步平均）；λ_D/dx 稳态值用 info 的 n_e、T_eff 重算。
- 答：怎么证明代码是对的（分层）？收敛判据是什么？为什么 nppc 256 略差而 Turner 选 512？Turner 的 S_e 为什么对应精确做功？

### 第 7 天：理论对照与压强依赖

- 读：`theory.md` 全文；`results/case1/theory.txt`、`results/case2/theory.txt`；`figures.md` §5；`results/case_comparison.png`。
- 跑：`python3 tools/theory_check.py runs/case1`；Case 2 太长（68 分钟）可以直接用 `results/case2/`。
- 算：用 info.txt 的 J_i 和 ⟨φ_c⟩ 手算射频 Child 定律的 s_m；用 EEPF 算 ν_m/ω。
- 答：哪些公式两个算例都合上、哪些只在 100 mTorr 合上、为什么（局域 vs 非局域）？T_e 是"哪个 T_e"？

### 第 8 天：扩展

- 读：`figures.md` §6–7；`tools/dual_frequency.py`；`tools/boltzmann2term.py` 文件头和 `tests/swarm/swarm.f90` 文件头；`code_guide.md` §15。
- 跑：`make check-swarm`（8 分钟）；双频 `inputs/dual_100.in`（4 核几分钟）+ `dual_frequency.py`。
- 答：双频为什么能/不能解耦？两项近似的动量转移截面要包含什么？高 E/N 下为什么要脉冲 Townsend 模型？MPI 怎么分、扩展性如何？

### 第 9–10 天：面试演练

- 把 `program_flow.md` §5 的 20 题各讲一遍，录音回听，每题 1 分钟内。
- 准备一个 10 分钟的介绍：问题（Turner 基准是什么）→ 方法（四步循环 + 零碰撞）→ 验证分层（单元 → 守恒 → 标定 → 基准 → 收敛）→ 结果（+2.1 %、表 III ≤ 1.2 %、理论对照）→ 扩展（Case 2、双频、swarm）→ 局限（kT_e 定义、网格 J·E、Case 3–4 需要集群）。
- 最可能的追问（准备答案）：
  1. 你自己写的部分和复用的部分各是什么？（em2 只借了 ran2 的算法和项目结构；一维 PIC、泊松、碰撞、诊断、工具全部新写。）
  2. 跑多久、多大规模？（Case 1 4 核 3 分钟、Case 2 68 分钟；粒子 10⁴–10⁵；MPI 粒子分解。）
  3. 最难的 bug 是什么？（离子截面网格原点错位 8 %；网格 J·E 的中心化；两项求解器的奇异漂移项和漏掉的非弹性动量转移；固定粒子数电子群在高 E/N 没有稳态。）
  4. 如果要做二维或加磁场怎么办？（§15；Boris 推进；二维泊松用多重网格或 FFT+三对角。）
  5. 和商业/开源代码比呢？（和 Turner 的五个代码、WarpX 的例子同一基准；差别在诊断精确记账和理论对照。）

## 4. 自测清单（全部能做到再去面试）

- 不看稿在白板上画出一步循环和时间层（x、ρ、φ、E 在整数步，v 在半整数步）。
- 推导：Thomas 算法、CIC 电荷守恒、做功恒等式、零碰撞法的等价性、弹性能量损失、背向散射 = 电荷交换。
- 说出每个数：ωpe·dt 0.12、λ_D/dx 3.6、P_coll 1.6 %、中心密度 +2.1 %、rms 1.5 %、噪声地板 1.4 %、三把尺 P_part/P_JE/P_VI = 125/108/125 W/m²、格式能量误差 0.5 %、Case 2 −0.4 %。
- 说出每个误差来源与大小：统计噪声、dt 一阶中心化（只影响网格 J·E）、有限网格、每格粒子数、截面数据本身（几 %）。
- 对每张图能讲"看什么、为什么是这样"。

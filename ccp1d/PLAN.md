# ccp1d 实施计划与测试方案

日期 2026-11-01。依据 `HANDOVER.md`（规格）以及本次对仓库、em2 源码、截面表和参考剖面的核对写成。
回答两个问题：代码怎么写（§2–§3），怎么证明它是对的（§4）。§1 是对交接文件的修正，**先读**。

## 0. 一页摘要

- 总体按 HANDOVER §6 的 M0 → M7 顺序做，目标和三档验收不变。
- 先定下 §1 的 15 条修正。影响代码结构的有五条：`dran2` 改成可存取状态的模块；测试做成独立的小程序链接同一套模块；功率用三种独立量测（P_part、P_JE、P_VI）；截面查表用线性插值；ν_max 加运行时守卫。
- 测试分五层：L0 构建、L1 单元测试（`make check`，约一分钟）、L2 物理标定（Langmuir 振荡、无碰撞鞘层、电子群）、L3 集成验收（测试模式自检 → Case 1 全尺寸 → 收敛）、L4 诊断自洽与双频。每项都有数字判据、命令和预计用时。
- 第一次编码会话做 M0 + M1，清单在 §5。

## 1. 对交接文件的修正与补充

| # | 交接文件的说法 | 核对结果 | 处理 |
| --- | --- | --- | --- |
| 1 | 稳态粒子数约为初始的 3–5 倍，maxnp 取初始的 8 倍 | 参考剖面的空间平均是 **0.29 n0**（中心 1.40e14 < n0 = 2.56e14，端点 3.6e13），稳态宏粒子数应**降到初始的约 30%**（每种约 1.9 万，中心每格约 280 个） | `history.dat` 里粒子数应先降后平，升高即是 bug；maxnp = 2× 初始已够（留 8× 也无害，内存 < 20 MB）；噪声按每格 150–280 个估计 |
| 2 | 电极节点的密度没有说明 | CIC 在节点 0 和 nx 只覆盖半个格体积，按 1/dx 归一后密度低一半 | `profiles.dat` 两端节点的 n_e、n_i 乘 2 并在表头注明；`compare_turner.py` 同时报告 129 点和内点 1..127 的偏差；端点若仍系统性偏离，只以内点为准 |
| 3 | `dran2.f90` 原样复制；续算"存随机数状态或接受不同" | `dran2` 的状态是函数内 `SAVE` 变量，外部存不了 | 同一算法放进模块 `random`（状态改为模块变量），加 `rng_get_state` / `rng_set_state`；Box–Muller 不缓存第二个数。代价 20 行，收益是续算逐位可复现，T8 因此成为最强的续算测试 |
| 4 | 截面查表 `k = int(ε/Δε)` 取整 | 热离子 ε_rel ≈ 0.04 eV，而 σ_iso 在 0.01 → 0.04 eV 间从 7.6e-19 降到 3.8e-19 m²（2 倍）；电子在阈值附近同理 | 在 Δε = 0.01 eV 均匀网格上做**线性插值**（多一次乘加）；T4 用原表直接插值做基准 |
| 5 | ν_max 的能量上限 1000 eV 写死 | Case 1 中离子 ε_rel ≤ 约 V0/2 = 225 eV、电子 < 450 eV，上限够用；改参数后未必 | namelist 参数 `numax_emax_e`、`numax_emax_i`（默认 1000 eV）；运行时若某粒子 ν_tot > ν_max 则计数，写进 `history.dat` 和 `info.txt` |
| 6 | 能量平衡只用 ⟨J·E⟩ | 一把尺分不清是物理错了还是诊断错了 | 三把尺：**P_part** = 推进时逐粒子做功 Σ q E_p·½(v⁻+v⁺)·dt（对蛙跳是恒等式：½m(v⁺²−v⁻²) = qE_p dt·½(v⁺+v⁻)）；**P_JE** = ∫(⟨J_e E⟩+⟨J_i E⟩)dx；**P_VI** = ⟨V(t)·I_tot(t)⟩，I_tot = 传导 + 位移电流，1D 中与 x 无关，周期平均后严格等于 ∫⟨J·E⟩dx。P_part 的收支（壁面动能通量 + 非弹性 + 弹性 + 储能变化）必须闭合到 < 0.1%（记账检查）；P_JE、P_VI 与 P_part 的差才是离散化误差（预期几 %）。面试报 P_VI |
| 7 | M3 要"无壁、周期、关电离产物"，暗示主程序加测试开关 | 开关越多主循环越难信 | 核心例程做成模块（random、xsec、mcc、poisson、deposit、push），`tests/*/` 各是一个小 program 链接这些模块、自带循环；主程序不出现测试分支；Makefile 加 `make check` |
| 8 | 没提一步跨多格 | 450 eV 电子一步走 4.4 dx | 吸收判据只看 x < 0 或 x > L，不做"一步一格"的假设；静电码的 CIC 沉积对此无要求 |
| 9 | `history.dat` 的列 | 分段计划需要每周期墙钟 | 加累计墙钟秒数一列、ν_max 超限计数一列 |
| 10 | 出错退出方式没说 | em2 卡死就是一个 rank `stop` 了、另一个在 `MPI_ALLREDUCE` 里等 | 一律 `MPI_ABORT`：np > maxnp、文件读失败、NaN（每周期检查 Σφ² 有限）都走这条路 |
| 11 | 云端"2 核 7 GB" | 本次容器 4 核 15 GB，每次不同 | 每次 `nproc`、`free -g` 后再定 `-np`；新容器都要 `./setup.sh`；numpy 2.x 没有 `np.trapz`，脚本用 `np.trapezoid` |
| 12 | M1 判据"中心相对壁面电势 ≈ 3.5 T_e，跑 2 µs" | 3.53 T_e 只在离子通量达到 Bohm 值后成立，离子从中心走到壁要 L/2/c_s = 4.3 µs，2 µs 时电势会偏高 | 跑 10 µs；与瞬时 T_e（= 2⟨ε_e⟩/3）算的 3.53 T_e 比，容差 20%；更硬的判据改为能量守恒和整体高斯定律（T11） |
| 13 | 测试模式"只看形状和量级" | nx = 32 时 λ_D/dx ≈ 0.4–0.5（T_e ≈ 2.5 eV）、ωpe·dt = 0.25–0.33，数值加热明显 | 只用它查崩溃、符号、对称性、量级，不用它调物理 |
| 14 | 双频参数没定 | — | f_lf = f/10 = 1.356 MHz（整除），V_lf = 0 / 100 / 200 V；n_avg 和 restart_every 取 10 的倍数（40、160） |
| 15 | M2 的期望值 n_g σ_k v | 一步最多一次碰撞使测得频率低一个因子 (1−e^(−ν_max dt))/(ν_max dt) = 0.992（Case 1 的 dt） | 比较时乘此因子，或测试用 dt/10 |

核对无误、可以直接用的数字：电子 ν_max = 8.654e7 s⁻¹（8.24 eV），P_e = 1.58%；离子 ν_max = 2.793e7 s⁻¹（上限处），P_i = 0.51%；w = 2.617e8 m⁻²；初始每种 65536 个；dx = 0.5234 mm；初始 ωpe·dt = 0.166、λ_D/dx = 1.43。
另：Turner 表 III 的 λ_D/dx = 3.72 配合 n = 1.35e14 反推 T_eff ≈ 9 eV，说明他用的是中心电子平均能量；`info.txt` 同时输出 ⟨ε_e⟩ 和按 T_eff = 2⟨ε_e⟩/3 算的 λ_D/dx，便于对照。

## 2. 代码结构

### 2.1 目录

```
ccp1d/
  Makefile          make（-O3 → ccp1d）/ make debug（→ ccp1d_dbg）/ make check / make clean
  *.f90             源文件（2.2），可执行文件生成在本目录
  inputs/           namelist 模板：case1.in  case1_test.in  sheath.in  dual_100.in  dual_200.in  conv_*.in
  tests/            单元与标定测试，每个子目录一个 program 或 run.sh（§4）
  tools/            Python 分析脚本（HANDOVER §5）+ check_run.py + common.py
  runs/             运行目录（.gitignore 已排除）
  results/          关键算例的小文件（进 git）
  data/             已有
  PLAN.md  PROGRESS.md  README.md
```

### 2.2 模块与依赖（也是编译顺序）

| 文件 | 内容 | 公开接口 | 来自 em2 |
| --- | --- | --- | --- |
| `input.f90` | namelist 变量、`read_input`、参数检查、不依赖截面的导出量（dx、dt、w、nsteps、maxnp）、`print_input(unit)` | — | input 的组织方式 |
| `global.f90` | CODATA 常数；`q(2)、m(2)`；网格数组 `phi, ex, rho, rho_s(0:nx,2), jx_s(0:nx,2), eps_s(0:nx), ioniz(0:nx)`；粒子数组 `x, vx, vy, vz (maxnp,2)`、`np(2)`（integer(8)）；计数器、壁面与碰撞账本；MPI 变量 | `add_particle(isp,x,vx,vy,vz)`、`remove_particle(isp,i)`（末尾搬入）、`abort_run(msg)`（MPI_ABORT） | global |
| `random.f90` | ran2 算法（状态为模块变量）、Box–Muller | `rng_init(seed)`、`rand_u()`、`rand_gauss()`、`rng_get_state(s)`、`rng_set_state(s)` | dran2（算法原样） |
| `xsec.f90` | 读 6 个表（电荷交换不读），插值到 Δε = 0.01 eV 的均匀网格，阈值，ν_max | `xsec_load(dir)`、`xsec_e(eps, sig(4))`、`xsec_i(eps_rel, sig(2))`（线性插值）、`compute_numax()` | 新 |
| `poisson.f90` | 常系数三对角的 Thomas 前扫一次预存，每步只做回代；E 的中心差分 + 两端二阶单侧差分 | `poisson_init()`、`poisson_solve(rho, v_left, v_right, phi, ex)` | 新 |
| `deposit.f90` | CIC 沉积 ρ（每步）与 ρ_e、ρ_i、J_e、J_i、能量密度（测量窗口内）；`MPI_ALLREDUCE` | `deposit_rho()`、`deposit_moments()` | setrho 的权重、curnt 的归约 |
| `push.f90` | 插值 E、推 v 和 x、逐粒子做功累加、电极吸收（能量直方图、计数、电荷） | `push_species(isp)` | move 的结构（先把 E 乘 (q/m)dt 再插值） |
| `mcc.f90` | 零碰撞法；电子弹性/激发×2/电离；离子各向同性/背向；电离产生新粒子；能量账本；电离位置记录 | `mcc_electrons()`、`mcc_ions()`、`scatter_isotropic(v)` | 新 |
| `diag.f90` | 每步：窗口累加、I_tot、能量；每周期：history 行、NaN 检查、续算文件；结束：profiles、iedf/eedf、eepf、xt_*、info | `diag_step()`、`diag_period()`、`diag_final()`、`restart_write()`、`restart_read()` | output + cont |
| `init.f90` | 分配；装粒子（位置均匀随机，速度 `rand_gauss`）；E = 0；蛙跳回推 −dt/2（通用写法） | `init_fresh()`、`init_restart()` | creator + setv 的职责 |
| `main.f90` | MPI、读入、新算/续算分支、主循环、计时、`max_wall_minutes` 优雅退出 | — | main |

用 `use mpi` 而不是 `include 'mpif.h'`（有类型检查）。所有文件 `implicit none`。

### 2.3 每步顺序与时间层

第 n 步开始：x 在 t_n，v 在 t_{n−1/2}，φ、E 在 t_n（电压 V(t_n)）。

1. `push_species(1)`、`push_species(2)`：E_p 插值，v_x += (q/m)E_p dt，x += v_x dt；做功 q E_p ½(v⁻+v⁺) dt 累加到 P_part；动能历史用时间中心形式 ½m v⁻·v⁺；x < 0 或 x > L 的粒子记账后删除。
2. `mcc_electrons()`、`mcc_ions()`：新电子、新离子放在父粒子位置，速度与父粒子同一时间层。
3. `deposit_rho()`（窗口内再 `deposit_moments()`），归约。
4. `poisson_solve` 用 V(t_{n+1})。
5. `diag_step()`；每周期末 `diag_period()`。

### 2.4 namelist（HANDOVER §4.1 之外新增）

```fortran
  numax_emax_e = 1000d0, numax_emax_i = 1000d0   ! ν_max 搜索的能量上限 [eV]
  maxnp_factor = 4                               ! maxnp = maxnp_factor × 初始每种粒子数 / NPROCS
  max_wall_minutes = 0                           ! > 0：墙钟到时跑完当前周期、写续算文件、正常退出
```

`max_wall_minutes` 让每个后台段都是同一条命令（续算文件在就接着跑），不必手算分段。

### 2.5 续算文件

- `restart_grid.bin`（rank 0 写，所有 rank 读）：格式版本号、nx、NPROCS、it、time、周期号、φ(0:nx)、窗口是否已开始、窗口累加量全部、壁面与碰撞账本。
- `restart_rankNNN.bin`：np(2)、各物种 x、vx、vy、vz、随机数状态。
- 读入时检查 NPROCS、nx 一致；若 n_periods 被改大而窗口已开始，重置窗口并在 info 里警告。
- 每 `restart_every` 周期和结束时写；文件名固定、覆盖式。

### 2.6 Makefile 要点

- `FCFLAGS_OPT = -O3 -march=native`；`FCFLAGS_DBG = -g -O0 -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -Wall -Wextra -finit-real=snan -fimplicit-none`（snan 比 `-finit-local-zero` 更能暴露未初始化的读）。
- 目标文件分两个目录（`obj/`、`obj_dbg/`），两种构建互不覆盖。
- 模块依赖写全：`main.o: global.o input.o random.o ... diag.o`。
- `make check`：用 debug 标志编译 `tests/` 下所有 program 并依次运行，任何一个 FAIL 则返回非零。
- 编译时把 `git rev-parse --short HEAD` 写进 `version.f90`，`info.txt` 里带 commit。

## 3. 里程碑

| # | 做什么 | 验收（§4 编号） | 预计 |
| --- | --- | --- | --- |
| M0 | Makefile、input、global、random、poisson、deposit、push、diag（info + history + 最简 profiles）、init、main、续算 | T1–T3、T8（n_g = 0 时续算逐位一致）、T10 | 1 个会话 |
| M1 | `inputs/sheath.in`（n_g = 0、V0 = 0、10 µs）+ `plot_history.py` | T11 | 同一会话或下一个 |
| M2 | xsec、mcc、`tests/mcc/` | T4–T7 | 1 个会话 |
| M3 | `tests/swarm/`（独立 program：固定 E、周期、电离只记能量不产生粒子） | T12 | 半个会话 |
| M4 | Case 1 测试模式 → `check_run.py` → 全尺寸后台跑 → `compare_turner.py` → `results/case1/` | T13、T14、T9 | 1 个会话 + 后台 |
| M5 | 收敛：换种子、dt/2、nx = 256、nppc 256 和 1024 | T15 | 后台 |
| M6 | 四类诊断齐全、能量平衡三把尺、全部作图脚本 | T16 | 1 个会话 |
| M7 | 双频 | T17 | 1 个会话 |
| 文档 | `ccp1d/README.md`、面试材料 | 别人照 README 能跑通 | 半个会话 |

每完成一个 T 就 commit 一次，提交信息里写判据和实测数字。

## 4. 测试方案

### 4.0 原则

- 新例程先用 `ccp1d_dbg`（`-fcheck=all -ffpe-trap`）跑过测试，再用 `-O3`。
- 每个测试打印 `PASS`/`FAIL` 和实测数字，非零退出码表示失败；`make check` 跑完 L1 全部不超过 2 分钟。
- 统计性判据用 3σ；"在噪声内"的参照物是**同参数换种子重跑**，不是拍脑袋的百分数。
- 单个测试在 2 核 10 分钟内跑完，更长的放后台（`max_wall_minutes`）。
- 物理量和轴都带单位。

### 4.1 L1 单元测试（`make check`，全部 Fortran program，debug 编译）

| 编号 | 目录 | 做法 | 判据 |
| --- | --- | --- | --- |
| T1 poisson | `tests/poisson/` | (a) ρ = 0、φ(0) = V0：φ 线性、E = −V0/L 常数；(b) ρ = 常数：抛物线；(c) ρ = sin(πx/L)：nx = 64/128/256 的误差比；(d) 随机 ρ：离散方程残差；(e) 整体高斯定律 ε0(E_nx − E_0) 与 Σρ_j dx | (a)(b) 最大误差 < 1e-10（相对）；(c) 误差比 ≈ 4（二阶）；(d) 残差 < 1e-12；(e) 偏差 O(dx²)，nx 加倍缩小 4 倍 |
| T2 deposit | `tests/deposit/` | (a) 单粒子在 x：两节点权重 (1−s)、s，Σρ_j dx = q w；(b) 粒子恰在节点上；(c) N 个均匀随机粒子：内点 ρ ≈ n0 q，端点 ≈ 一半；(d) J 同上乘 v_x；(e) 1 rank 与 2 rank（各一半粒子）结果比较 | (a)(b) 1e-14；(c) 内点偏差 < 3σ = 3/√(每格粒子数)，端点为内点的 0.5 ± 3σ；(e) 相对 1e-13 |
| T3 random | `tests/random/` | (a) 均匀分布均值、方差；(b) 高斯均值、方差、峰度；(c) 存/取状态后序列逐位相同；(d) 不同种子序列不同 | (a)(b) 3σ 内（N = 1e7）；(c) 逐位；(d) 前 1000 个无一相同 |
| T4 xsec | `tests/xsec/` | (a) 表格点上查表等于表值；(b) 表格点之间等于原表线性插值；(c) 749.99 eV 以上取末值；(d) 阈值以下为 0；(e) ν_max | (a) 1e-12；(b) 相对 < 1e-3（离子 0.01–0.1 eV 段单独检查）；(e) 电子 8.654e7 ± 0.5%（≈ 8.2 eV），离子 2.793e7 ± 0.5% |
| T5 mcc 频率 | `tests/mcc/` | 不解场、E = 0、无壁。电子全给 ε0 = 10、30、100 eV（每次碰后把能量重置回 ε0，只计数），N = 1e5、2000 步；离子给 10 eV、T_g = 0（期望值精确）；统计各过程事件数/(N·步数·dt)，并统计空碰撞份额 | 每个过程 = n_g σ_k v × (1−e^(−ν_max dt))/(ν_max dt)，3σ Poisson 内；Σν_k < ν_max；空碰撞份额 = 1 − Σν_k/ν_max |
| T6 mcc 运动学 | `tests/mcc/` | 逐次碰撞前后比较：弹性 ε' = ε[1 − (2m_e/M)(1−cosχ)]（χ 用实际新旧速度夹角算）；激发 ε' = ε − ε_th；电离 ε1 = ε2 = (ε − 24.59)/2；离子各向同性：|g'| = |g|、v_cm 不变；背向：v_i' = v_n；cosχ 直方图均匀、φ 均匀；新离子三个分量的方差 = k T_g/M | 能量关系相对 1e-12；直方图 χ² 检验 p > 0.01（20 格、1e6 样本）；方差 3σ 内 |
| T7 mcc 弛豫 | `tests/mcc/` | (a) 10 eV 离子在 300 K He 中，dt = 1e-9 s 跑 20 µs：平均能量 → 3/2 k T_g = 0.0388 eV，各分量方差 → k T_g/M；(b) 10 eV 电子只开弹性（低于阈值），初始冷却斜率 dε/dt = −(2m_e/M) ν_el(ε) ε | (a) 终态 ± 5%；(b) 斜率 ± 5%（N = 1e5，取前 2 µs） |
| T8 restart | `tests/restart/run.sh` | 测试模式参数：20 周期一口气跑 vs 10 周期 + 续算 10 周期 | `history.dat`、`profiles.dat` 逐位相同（`cmp`）；若 MPI 归约顺序造成末位差异，放宽到相对 1e-12 |
| T9 MPI | `tests/mpi/run.sh` | 同一输入 `-np 1` 与 `-np 2` | 中心密度、粒子数、P_part 的差在换种子噪声内；2 核加速比 > 1.6 |

### 4.2 L2 物理标定

| 编号 | 做法 | 判据 | 用时 |
| --- | --- | --- | --- |
| T10 Langmuir 振荡 | `tests/langmuir/`：独立 program，冷电子、离子不动、两壁接地 V = 0，电子位移 δx = A sin(πx/L)（A = 0.1 dx），Case 1 的 n0、dx、dt；记录场能随时间 | 场能振荡频率 = 2ωpe，ωpe 偏差 < 1%（蛙跳的理论偏差 (ωpe dt)²/24 = 0.1%）；振幅 50 个周期内不增长 > 5% | 秒级 |
| T11 无碰撞鞘层（M1） | 主程序，`inputs/sheath.in`：n_g = 0、V0 = 0、Case 1 其余参数，10 µs（≈ 2.3 个离子渡越时间），每 0.1 µs 写一行 history | (i) 无 NaN，粒子数单调减；(ii) 能量守恒：动能 + 场能 + 已送到壁面的动能 = 常数，漂移 < 1% 初始动能；(iii) 整体高斯定律每步成立到相对 1e-10；(iv) φ 剖面左右对称（< 噪声）；(v) 渡越时间之后 φ_center − φ_wall = 3.53 T_e(t) ± 20%，T_e = 2⟨ε_e⟩/3；(vi) 之后 Γ_e ≈ Γ_i（每壁，± 10%）；(vii) 壁面 IEDF 峰 ≈ e(φ_center − φ_wall)（顺带检验 IEDF 诊断） | 2–5 分钟 |
| T12 电子群（M3） | `tests/swarm/`：独立 program，固定均匀 E、周期边界、只有电子、电离只扣能量不产粒子；E/N = 50、100、200、500 Td；输出平均能量、漂移速度、α/N | 内部精确判据：稳态 e·E·v_d = 每电子单位时间碰撞能量损失（± 2%）；结果表存 `results/swarm/`，以后与 BOLSIG+/LXCat Biagi 比（预期 ± 10%，因为截面一致、各向同性） | 5 分钟 |

### 4.3 L3 集成与基准

**T13 Case 1 测试模式自检**（`tools/check_run.py <run_dir>`，M4 起每次改动后都跑，也是回归基线）：

| 检查 | 判据（测试模式 / 全尺寸） |
| --- | --- |
| 所有输出无 NaN/Inf | 必须 |
| 粒子数先降后平；最后 20% 周期的线性漂移 | < 5% / < 2% |
| 0 < (N_i − N_e)/N_i（鞘层正电荷） | < 15%（两个 5 mm 鞘层的物理值约 5–10%） |
| φ(0, t) = V0 sin(2πft)（从 xt_phi 第 0 列） | 相对 1e-6 |
| 中心 n_i 与 1.405e14 | 2 倍以内 / < 5%（目标 2%） |
| 剖面左右对称 max\|n(x) − n(L−x)\|/n_c（参考值本身 0.1%） | < 5% / < 2% |
| Γ_i 左 ≈ Γ_i 右；每壁 Γ_e ≈ Γ_i（稳态每周期净电荷为零） | < 10% / < 5% |
| 离子撞壁能量：最大 < V0 + 50 eV，平均在 0.2–0.6 V0 | 必须 |
| 中心 ⟨ε_e⟩ | 3–15 eV（Turner 表 III 反推 T_eff ≈ 9 eV） |
| 能量记账：P_part 收支闭合 | < 0.5% / < 0.1% |
| P_JE/P_part、P_VI/P_part | 10% 内 / 5% 内 |
| ν_max 超限计数 | = 0 |
| 每周期墙钟 | 记录，用来排全尺寸 |

**T14 Case 1 全尺寸**（`compare_turner.py`）：中心密度偏差 < 5%（目标 2%）；内点最大相对偏差 < 10%，均方相对偏差 < 5%；T13 全尺寸一列全部通过。通过后把 info、history、profiles、iedf、eepf 复制到 `results/case1/`，并把 `compare_turner.png` 放进 README。

**T15 收敛**（`convergence.py`，5 个后台算例）：基线、换种子、steps_per_period 800、nx 256、nppc 256、nppc 1024。判据：任意两个的剖面偏差（5 点滑动平均后）≤ 2 × "基线 vs 换种子"的偏差；中心密度两两 < 2%。换种子那一个就是噪声地板。

### 4.4 L4 诊断自洽与扩展

**T16 诊断自洽**（M6，`check_run.py` 增加）：

- ∫IEDF dε = history 里的离子通量（同一批计数，相对 1e-6）；IEDF 双峰，峰间距与 Benoit-Cattin–Bernard 公式同量级（He 离子轻，渡越时间 ≈ 射频周期，分裂明显）。
- ∫F(ε)dε = 中心 n_e（± 1%）；EEPF 半对数图低能段直线，斜率给 T_e，与 ⟨ε_e⟩ 一致到 30%（双温时取低能段）。
- xt_ioniz 对 x 和相位积分 = history 里每周期电离数（± 1%）；xt_PeE 积分 = P_JE 的电子部分（± 1%）。
- xt_phi 第 0 列 = V0 sin（检验相位分格）。
- 定性：鞘层边缘随相位振荡；J_e·E 极大值在鞘层边缘（随机加热）；电离峰在鞘层边缘附近。

**T17 双频**（M7）：`dual_100.in`、`dual_200.in`（f_lf = 1.356 MHz，n_avg = 40）。判据：IEDF 的 10–90% 能量宽度随 V_lf 单调增大；最大离子能量 ≈ V0 + V_lf；T13 的全部检查仍通过；密度变化如实报告，不设通过线。

### 4.5 回归

每次改动后：`make check`（L1）+ 测试模式 + `check_run.py`，并与 `results/case1_test/` 里存的基线比（中心密度、粒子数、P_part 在噪声内）。改了物理（碰撞、边界）的提交必须重跑 T5–T7 和测试模式。

## 5. 第一次编码会话清单（M0 + M1）

1. `./setup.sh`；`nproc`。
2. Makefile（两套标志、两套目标目录、`check` 目标）。
3. 按 §2.2 顺序写 input → global → random → poisson → deposit → push → diag → init → main；每写完一个模块就 `make debug` 一次。
4. T1、T2、T3 通过 → commit "M0: skeleton, poisson/deposit/random tests pass"。
5. T10 Langmuir → commit。
6. 主程序跑通测试模式参数（n_g = 0、V0 = 0 时就是鞘层问题的短版），T8 续算逐位一致 → commit。
7. `inputs/sheath.in` 跑 10 µs，`plot_history.py` 画能量、粒子数、电势，T11 的 7 条逐条记入提交信息 → commit。
8. 更新 `PROGRESS.md`（做了什么、数字、下一步、未解决），push。

## 6. 风险与预案

- **每步耗时未知。** 全尺寸前先跑 10 周期看 `history.dat` 的墙钟列；粗估稳态每步约 4 万粒子、1–2 ms，512000 步约 15–40 分钟，2 小时后台上限之内；`max_wall_minutes` 兜底。
- **中心密度偏差达不到 5%。** 排查顺序：电极节点因子 2（§1 第 2 条）；泊松用的是 V(t_{n+1}) 而不是 V(t_n)；E 两端的单侧差分符号；电离能量平分和新离子取 T_g；离子碰撞用 ε_rel（μ = M/2）而不是实验室能量；截面线性插值；ν_max 超限计数；换种子看噪声。
- **数值加热。** 全尺寸 λ_D/dx = 1.4–1.9（初始到稳态）本来就不宽裕，T11 的能量漂移数字是基线；测试模式的加热不必修。
- **MPI 归约顺序破坏逐位复现。** T8 放宽到 1e-12，info 里注明。
- **云端网络只通 github。** LXCat/BOLSIG+ 对照在本地做，云端只存表（T12）。

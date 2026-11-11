# ccp1d 代码说明：每个文件、每个例程做什么、用什么方法

读法：先看 `program_flow.md`（为什么这样设计、算法推导、面试追问），再用本文逐个文件对照源码。本文按"文件 → 数据 → 例程"的顺序写，例程表里给出调用者、输入输出和方法；最后是输出文件格式、测试和工具的清单，以及"想改什么去哪里改"的索引。源码总共约 1 900 行 Fortran（12 个文件）加 2 000 行 Python 工具。

## 0. 文件一览与依赖

| 文件 | 行数 | 角色 | 依赖（use） |
| --- | --- | --- | --- |
| `input.f90` | 150 | namelist 参数、合法性检查、导出量 dx/dt/weight/maxnp、驱动电压 V(t) | — |
| `global.f90` | 87 | 物理常数、物种、MPI 变量、网格数组、粒子数组、步数计数器、增删粒子 | input, mpi |
| `random.f90` | 75 | ran2 随机数（状态可存取）、Box–Muller 高斯 | — |
| `poisson.f90` | 80 | 一维泊松方程 Thomas 解法、电场 | —（独立，便于单元测试） |
| `deposit.f90` | 86 | CIC 沉积：电荷密度、数密度、物种矩 | global |
| `push.f90` | 135 | 蛙跳推进、电极吸收、做功/动能/Ramo 电流累加 | global |
| `xsec.f90` | 197 | 截面表读入、均匀能量网格、ν_max | global |
| `mcc.f90` | 176 | 零碰撞法蒙特卡洛碰撞（电子 4 过程、离子 2 过程）、碰撞账本 | global, random, xsec |
| `diag.f90` | 678 | history、测量窗口累加、剖面/能谱/x–t 输出、能量平衡、续算文件 | 以上全部 + version |
| `init.f90` | 55 | 新算例装载与蛙跳起步 | global, random, poisson, deposit, diag |
| `main.f90` | 109 | 程序入口与时间循环 | 全部 |
| `version.f90` | 4 | 由 `make` 生成，记录 git 提交号写进 info.txt | — |

编译顺序即上表顺序（Makefile 的 `DEPS` 规则）。所有物理量用 SI 单位，只有截面表和输出文件里的能量用 eV。

## 1. input.f90 —— 参数

**namelist `&ccp1d`**（默认值是 Turner Case 1）：

| 变量 | 含义 | Case 1 默认 |
| --- | --- | --- |
| `L`, `nx` | 电极间距 [m]、格数 | 0.067, 128 |
| `freq`, `V0` | 驱动频率 [Hz]、幅值 [V]，φ(0,t) = V₀ sin(2πft) | 13.56e6, 450 |
| `freq_lf`, `V0_lf` | 第二（低）频率项，加在 φ(0,t) 上 | 0, 0 |
| `n_gas`, `T_gas` | 氦气密度 [m⁻³]（0 = 无碰撞）、温度 [K] | 9.64e20, 300 |
| `n0`, `Te0`, `Ti0` | 初始密度 [m⁻³]、初始温度 [K] | 2.56e14, 30000, 300 |
| `nppc` | 每格每物种初始宏粒子数 | 512 |
| `steps_per_period`, `n_periods` | 每周期步数、总周期数（续算时是累计目标） | 400, 1280 |
| `periods_this_run` | > 0 时本次只跑这么多周期然后写续算文件 | 0 |
| `n_avg`, `n_phase` | 测量窗口 = 最后 n_avg 个周期；x–t 图每周期的相位格数（须整除 steps_per_period） | 32, 100 |
| `seed` | 随机种子（每个 rank 用 seed + rank） | 12345 |
| `xsec_dir` | 截面表目录 | `../../data/he` |
| `restart`, `restart_every` | 是否续算、每多少周期写续算文件 | F, 160 |
| `iedf_emax/iedf_de`, `eepf_emax/eepf_de` | 壁面能谱、体内 EEPF 的直方图范围与分辨率 [eV] | 600/1, 100/0.5 |
| `numax_emax_e/i` | 搜索 ν_max 的能量上限 [eV] | 1000, 1000 |
| `maxnp_factor` | 粒子槽数 = 该因子 × 每 rank 初始粒子数 | 4 |
| `max_wall_minutes` | > 0 时到墙钟后写续算文件干净退出 | 0 |
| `dump_initial` | 输出 t = 0 的粒子和网格（画初始条件图用） | F |

**导出量**：`dx = L/nx`，`period = 1/freq`，`dt = period/steps_per_period`，`np_init_total = nx·nppc`（每物种、所有 rank 合计），`weight = n0·L/np_init_total`（一个宏粒子代表单位电极面积上多少个真实粒子，单位 m⁻²），`maxnp = maxnp_factor × ceil(np_init_total/nranks)`，直方图格数 `nbin_iedf`、`nbin_eepf`。

| 例程 | 调用者 | 做什么 |
| --- | --- | --- |
| `read_input(nranks)` | main | 打开第一个命令行参数指定的文件（默认 `ccp1d.in`），读 namelist，然后 `derive_parameters` |
| `derive_parameters(nranks)` | read_input、各测试 | 11 项合法性检查（nx ≥ 4、steps_per_period 是 n_phase 的倍数、n_avg ≤ n_periods……），算导出量。测试程序直接给 namelist 变量赋值再调它 |
| `v_drive(t)` | main、init、diag | 驱动电极电势 V₀ sin(2πft) + V_lf sin(2πf_lf t)（pure 函数） |
| `print_input(u)` | diag.write_info_header | 把所有参数和导出量按 `key = value  # unit` 写进 info.txt |

## 2. global.f90 —— 数据与常数

- **常数**：CODATA 2018 的 e、m_e、ε₀、k_B；He⁺ 质量 6.67×10⁻²⁷ kg 按 Turner 原文。
- **物种**：`nsp = 2`，`q = (−e, +e)`，`m = (m_e, M)`，`qm = q/m`，名字 `e`、`He+`。整个代码里物种下标 1 = 电子、2 = 离子。
- **网格**：节点 j = 0..nx 在 x_j = j·dx；数组 `phi`、`ex`、`ex_old`（上一步的场，给 J·E 和位移电流做时间中心化）、`rho`（所有 rank 求和后的电荷密度）。
- **粒子**："数组的结构"：`x(maxnp,nsp)`、`vx`、`vy`、`vz`，`np(isp)` 是在用的槽数。每个 rank 只持有自己的粒子，但持有整个网格。
- **计数**：`it`（已走步数，integer(8)）、`iperiod`、`time = it·dt`、`mcc_on`。

| 例程 | 做什么 |
| --- | --- |
| `allocate_arrays` | 分配并清零网格和粒子数组 |
| `add_particle(isp, x, vx, vy, vz)` | 追加到末尾；槽用完时 `abort_run`（提示增大 maxnp_factor） |
| `remove_particle(isp, i)` | 把最后一个粒子搬进槽 i，`np = np − 1`，O(1)；PIC 不关心粒子顺序 |
| `abort_run(msg)` | 打印后 `MPI_ABORT`：任一 rank 出错都让整个作业退出，避免其他 rank 在 ALLREDUCE 里死等 |

## 3. random.f90 —— 随机数

Numerical Recipes 的 `ran2`：两个线性同余生成器（乘子 40014 模 2147483563，乘子 40692 模 2147483399，Schrage 法避免溢出）相减，再经 32 格 Bays–Durham 洗牌表，周期约 2×10¹⁸。em2 的 `dran2` 是同一算法；这里改成模块，状态是 35 个整数（idum、idum2、iy、iv(32)）。

| 例程 | 做什么 |
| --- | --- |
| `rng_init(seed)` | 8 次预热后填洗牌表；main 用 `seed + myrank`，各 rank 序列不同 |
| `rand_u()` | (0, 1) 内均匀数，上限 1 − 1.2×10⁻⁷，永远不返回 0 或 1 |
| `rand_gauss()` | Box–Muller：sqrt(−2 ln u₁)·cos(2πu₂)，不缓存第二个数，保证状态完整 |
| `rng_get_state / rng_set_state(s(35))` | 续算文件存取；T8 验证续算逐位一致 |

用到随机数的地方：装载（init，每粒子 7 个数）、碰撞判定与过程选择、散射方向、靶原子速度、新离子速度（mcc）。

## 4. poisson.f90 —— 场方程

方程 (φ_{j−1} − 2φ_j + φ_{j+1})/dx² = −ρ_j/ε₀，j = 1..nx−1，边界 φ₀ = V(t)，φ_nx = 0。三对角系数 a = 1、b = −2、c = 1 每步相同，所以 Thomas 算法的前扫系数只算一次。模块私有、自带 ε₀，不依赖 global，便于 `tests/poisson` 单独测。

| 例程 | 做什么 |
| --- | --- |
| `poisson_init(nx, dx)` | 预存 `rb(j) = 1/(b − a·cp(j−1))` 和 `cp(j) = c·rb(j)` |
| `poisson_solve(rho, v_left, v_right, phi, ex)` | 右端 d_j = −dx²ρ_j/ε₀，边界值移到 d₁、d_{nx−1}；前扫 d_j ← (d_j − d_{j−1})·rb_j；回代 φ_j = d_j − cp_j φ_{j+1}。电场：内部中心差分 E_j = −(φ_{j+1} − φ_{j−1})/2dx，两端二阶单侧差分。每步 O(nx)，每个 rank 都解一遍（ρ 已经 ALLREDUCE，结果相同，省一次广播） |
| `efield_faces(phi, ef)` | 面上的场 E_{j+1/2} = −(φ_{j+1} − φ_j)/dx：离散高斯定律 ε₀(E_{nx−1/2} − E_{1/2}) = Σ_{j=1}^{nx−1} ρ_j dx 对它精确成立，场能 ½ε₀ΣE²dx 也用它 |

## 5. deposit.f90 —— 粒子到网格

CIC/线性权重：j = ⌊x/dx⌋，s = x/dx − j，节点 j 得 (1−s)/dx，节点 j+1 得 s/dx，所以节点值直接是密度。`min(int(xd), nx−1)` 防止舍入把 x 刚好落在 L 的情况推出数组（推进时 x ≥ L 的粒子已被吸收）。

| 例程 | 归约 | 做什么 |
| --- | --- | --- |
| `deposit_rho()` | MPI_ALLREDUCE | 两物种电荷密度 [C/m³]，所有 rank 得到完整 ρ，每步一次 |
| `deposit_number_density(isp, dens)` | MPI_ALLREDUCE | 一个物种的数密度（history 的中心密度、初始状态输出用） |
| `deposit_species_moments(isp, dens, flux, ekd, ekxd)` | 不归约 | 单位权重的密度、n·v_x、动能密度、x 方向动能密度，rank 本地；diag 在窗口内每步累加，窗口结束时归约一次 |

## 6. push.f90 —— 粒子推进

进入时 x 在 t_n、v_x 在 t_{n−1/2}；离开时 x 在 t_{n+1}、v_x 在 t_{n+1/2}。

| 例程 | 做什么 |
| --- | --- |
| `push_init` | 分配节点数组 `a(0:nx)`、做功沉积数组 `work_dep(0:nx,nsp)`、吸收缓冲区（4096，可增长） |
| `push_species(isp)` | 见下 |
| `record_absorbed(side, isp, e)` | 把被吸收粒子的（哪一侧、物种、动能）放进缓冲区，diag 在本步结束时做直方图 |
| `kinetic_sums(isp, ekx, ekperp)` | 本 rank 粒子的 x 方向与垂直方向动能和（history 的 T_x、T_⊥） |
| `kinetic_energy(isp)` | 本 rank 粒子总动能（窗口起止的储能） |

`push_species` 的循环：

1. 先算节点量 `a_j = (q/m)·dt·E_j`，粒子循环里只剩一次插值 `ep = a_j(1−s) + a_{j+1}s`。
2. `vnew = vold + ep`；`x ← x + vnew·dt`。
3. 累加 `wsum += ep·(vold + vnew)`（乘 m/2 就是做功 Σ qE_p dt·½(v⁻+v⁺)，它**恒等于**动能增量，所以能量记账能闭合到舍入误差）；`vsum += vnew`（Shockley–Ramo 电流 (w/L)Σqv_x）。窗口内（`tally_work`）再把每个粒子的做功按**旧位置**的 CIC 权重沉积到 `work_dep`——这就是加热剖面 `P_work`，积分精确等于总做功。
4. 若 `x < 0` 或 `x ≥ L`：记壁面计数和动能（分侧、分物种），`record_absorbed`，`remove_particle`——注意循环用 `do while`，删除后不递增 i，因为最后一个粒子被搬进了槽 i，下一轮会处理它。否则累加留在域内粒子的动能。
5. 出口：`work_step`、`ekin_step`、`vxsum_step`、`nwall_step(2,nsp)`、`ewall_step(2,nsp)`，供 diag 读取。

没有磁场，v_y、v_z 推进时不变（无碰撞时 T_⊥ 必须严格不变，T11 用它检验）。

## 7. xsec.f90 —— 截面

| 例程 | 做什么 |
| --- | --- |
| `xsec_load()` | 读 4 张电子表（弹性、激发 1、激发 2、电离）和 2 张离子表（各向同性、背向）；阈值 = 各阈值过程表的首个能量；网格原点 e0 = 各表最小首能量（电子 0 eV，离子 1×10⁻⁴ eV），步长 Δε 电子 0.01、离子 0.001 eV，格数 = (numax_emax − e0)/Δε；`to_grid` 把每张表线性插值到网格（原点对齐表的节点，所以插值精确复现原表）；最后在网格上找 ν_max = max n_g σ_tot v（电子 v = sqrt(2ε/m_e)；离子用相对能量，g = sqrt(2ε_rel/μ)，μ = M/2），乘 (1 + 10⁻⁶) 的余量 |
| `xsec_read_table(fname, e, s, n)` | 两列表，跳过 # 行，检查能量单调递增 |
| `to_grid(e, s, n, e0, de, ng, sg, threshold)` | 表以下取首值（阈值过程取 0），表以上取末值，中间线性插值 |
| `xsec_e(eps, sig(4))` | O(1) 查表：k = ⌊(ε − e0)/Δε⌋，两点线性插值；超出网格取末值 |
| `xsec_i(eps_rel, sig(2))` | 同上，自变量是质心系相对能量 |

公开变量：`eth_e(4)`（阈值 0、19.82、20.61、24.59 eV）、`numax_e/i` 与对应能量 `enumax_e/i`、`xsec_loaded`。

## 8. mcc.f90 —— 蒙特卡洛碰撞

| 例程 | 做什么 |
| --- | --- |
| `mcc_init` | P_coll = 1 − exp(−ν_max dt)（电子 1.6 %、离子 0.5 %，Case 1）；靶原子热速度 σ_v = sqrt(kT_g/M)；μ = M/2；电离位置缓冲区 |
| `mcc_step` | 清零本步账本，先电子后离子 |
| `isotropic_direction(ux,uy,uz)` | 球面均匀单位矢量：cosθ = 1 − 2R，φ = 2πR' |
| `mcc_electrons` | 对进入时已有的 np(1) 个电子（本步新生的不再测）：(1) `rand_u() ≥ P_coll` 则跳过；(2) ε = ½m_e v²，`xsec_e`，ν_k = n_g σ_k v，累计 c₁ < c₂ < c₃ < ν_tot；ν_tot > ν_max 时计数（守卫）；(3) r = R·ν_max，r ≥ ν_tot 为空碰撞；(4) 抽新方向；弹性：ε' = ε[1 − (2m_e/M)(1 − cosχ)]，χ 是新旧方向夹角，损失记入 `eloss_el_step`；激发：ε' = ε − 阈值，阈值记入 `einel_step`；电离：剩余能量平分，第二个电子独立方向、`add_particle(1)`，新 He⁺ 在同一位置取气体 Maxwell 速度、`add_particle(2)`，热能记入 `enewion_step`，位置记入 `record_ionization` |
| `mcc_ions` | 对每个离子：P 测试；抽靶原子速度 v_n（三个高斯）；g = v_i − v_n，ε_rel = ½μg²，`xsec_i`；ν₁（各向同性）、ν₂（背向）；空碰撞判定；各向同性：v' = v_cm + ½|g|n̂；背向：v' = v_n（即电荷交换）；离子动能变化记入 `egas_i_step` |
| `record_ionization(x)` | 本步电离位置，diag 用 CIC 做电离率剖面和 x–t 图 |

账本（单位权重，能量 J）：`n_test_step`、`n_el_step`、`n_exc_step(2)`、`n_ion_step`、`n_icoll_step(2)`、`numax_exceed_step`、`eloss_el_step`、`einel_step`、`enewion_step`、`egas_i_step`。

## 9. diag.f90 —— 诊断、输出、续算

**两类累加器**：粒子类（前缀 `w_`、`p_`、`c_`）是本 rank 的"单位权重"和，窗口结束或每周期用 `reduce0` 求到 rank 0；网格类（φ、E、x–t 的 φ、计数）在每个 rank 上相同，不归约。续算文件按这个划分存：粒子类进 `restart_rankNNN.bin`，网格类进 `restart_grid.bin`。

| 例程 | 何时 | 做什么 |
| --- | --- | --- |
| `diag_init` | 启动 | 分配窗口数组（剖面 (0:nx,nsp)、x–t (0:nx, 0:n_phase−1, 7)、壁面直方图 (nbin, 2, nsp)、EEPF），清零，记墙钟起点 |
| `diag_window_check` | 每周期开头 | iperiod ≥ n_periods − n_avg 时打开窗口：清零、记窗口起点的动能/场能/时间，`tally_work = .true.`；若 n_periods 被改大则重置 |
| `diag_step` | 每步解场之后 | (a) Ramo 电流 I = (w/L)Σqv_x + (ε₀/L)dV/dt，与 V(t_{n+1/2}) 相乘累加成 P_VI；(b) 本周期累加做功、壁面能量、壁面计数、电离/激发数、ν_max 超限数；(c) 离散高斯定律残差（取周期内最大值）；(d) 窗口内：相位格 p = (it mod spp)·n_phase/spp；每物种 `deposit_species_moments` 累加密度、通量、能量密度；J·E 用 E_mid = ½(E_old + E_new)；x–t 七个量；φ、E 平均；P_VI 的传导/位移分量；碰撞账本；电离位置 CIC；壁面能谱直方图（从 push 的缓冲区）；中心 10 % 的 EEPF 直方图；最后清空吸收缓冲区 |
| `history_open(append)` | 启动 | 新建（写 31 列表头）或追加 history.dat |
| `diag_period` | 每周期末 | 中心 5 个节点的密度、动能和；25 个数 `reduce0`；rank 0 写一行（通量 = 计数·w/T，撞击能 = 能量/计数，功率 = 做功·w/T，T_x = 2⟨ε_x⟩，T_⊥ = ⟨ε_⊥⟩）；按进度打印；清零周期量 |
| `write_info_header` | 新算例启动 | info.txt 前半：git 提交、日期、参数、初始等离子体参数（ωpe·dt、λ_D/dx、v_th·dt/dx、Bohm 速度、离子渡越时间）、碰撞参数（ν_max、P_coll、阈值、气压） |
| `dump_initial_state` | dump_initial | `init_particles_rankNNN.dat`（物种、x、v）和 `init_grid.dat`（x、ρ、φ、E、n_e、n_i） |
| `diag_final` | iperiod == n_periods | 窗口终点的动能/场能；归约全部粒子类累加器；写 profiles.dat、iedf/eedf、eepf.dat、8 张 x–t 图、info.txt 结果段（见 §13 的格式与能量平衡定义） |
| `write_xt(fn, what, a, factor)` | diag_final | 一张 x–t 矩阵：行 = 相位格，列 = 节点，除以该相位格的步数 |
| `restart_write` | 每 restart_every 周期、本次结束、墙钟到 | 每 rank：魔数 `CCP1DR03`、nx、nprocs、it、iperiod、粒子、随机数状态、累计量、窗口标志与粒子类窗口累加器；rank 0 另写网格文件：φ、E、time、网格类窗口量 |
| `restart_read` | restart = .true. | 反向读入，检查魔数、nx、nprocs、maxnp；恢复 `tally_work`；在 info.txt 追加 `restart_at_period` |
| `trapz`, `wfield`, `reduce0`, `wall_clock`, `window_reset` | 内部 | 梯形积分；场能 ½ε₀ΣE_{j+1/2}²dx；MPI_REDUCE 到 rank 0；墙钟；清零 |

## 10. init.f90 与 main.f90

**`init_fresh`**：每 rank 装 nr = N/nprocs 个（余数给前几个 rank）；每粒子按固定顺序抽 7 个数：x = L·U，v_x、v_y、v_z = σ_v·高斯（分开写成独立语句，因为 Fortran 不规定同一表达式里函数引用的求值顺序，而抽样顺序必须可复现）；沉积、解 t = 0 的场（V(0) = 0）；可选输出初始状态；最后 v_x ← v_x − ½(q/m)E_p dt 把速度回推半步。

**`main`**：MPI 初始化 → `read_input` → `rng_init(seed + rank)` → 分配 → `poisson_init`、`push_init`、`diag_init` → n_gas > 0 时 `xsec_load`、`mcc_init` → 续算则 `restart_read` + 追加 history，否则 `init_fresh` + info 表头 + history 第 0 行。时间循环按周期：`diag_window_check`；每步 `push_species(1)`、`push_species(2)`、`mcc_step`、`deposit_rho`、`ex_old = ex`、`it += 1`、`poisson_solve(rho, V(t_{n+1}), 0)`、`diag_step`；周期末 `diag_period`；rank 0 判断墙钟并广播；按 `restart_every`、本次终点或墙钟写续算文件；只有 iperiod == n_periods 才 `diag_final`。停止条件三选一：周期到、`periods_this_run` 到、`max_wall_minutes` 到。

## 11. Makefile

`make`（-O3 -march=native，目标 `./ccp1d`，对象在 obj/）、`make debug`（-O0 -fcheck=all -ffpe-trap=invalid,zero,overflow -finit-real=snan -Wall，目标 `./ccp1d_dbg`，对象在 obj_dbg/）、`make check`（T1–T8、T10，用调试对象链接）、`make check-rng`（三个统计测试，优化编译）、`make check-swarm`（T12：编译 swarm，跑 PIC 电子群，解两项 Boltzmann，比较）、`make clean`。`version.f90` 每次 make 由 git 提交号生成，只在变化时重写（避免无谓重编译）。

## 12. 测试（tests/）

| 测试 | 文件 | 检验什么 | 方法 |
| --- | --- | --- | --- |
| T1 | poisson/test_poisson | 泊松求解器 | 三个解析解（ρ = 0 线性电势、均匀 ρ 抛物线、正弦 ρ）、nx 加倍误差比 4（二阶）、离散残差、离散高斯定律 |
| T2 | deposit/test_deposit | CIC 沉积 | 单粒子权重、粒子在节点上、均匀装载的内部密度与端点半格、物种矩、1 rank 与 2 rank 结果一致（同一批粒子分给两个 rank） |
| T3 | random/test_random | 生成器 | 均匀/高斯的各阶矩、状态存取逐位、不同种子不同序列 |
| T3+ | random/test_tails, test_seeds, test_corr | 装载统计 | 2×10⁷ 粒子的尾部计数对称性；1000 个种子复现 Case 1 一个 rank 的装载并收集 6 个统计量的 z 值（rms 应为 1）；滞后 1–16 的序列相关、同一粒子三个分量的相关 |
| T4 | mcc/test_xsec | 截面查表 | 与原表直接插值比较、表外行为、阈值、ν_max |
| T5–T7 | mcc/test_mcc | 碰撞模块（无场无壁） | T5 单能电子（10/30/100 eV）与离子的各过程计数 = n_gσv·(1−e^{−ν_max dt})/(ν_max dt)；T6 运动学：弹性能量损失公式、激发/电离能量、新方向各向同性 χ²、新离子热速度、离子碰撞的质心系不变量；T7 弛豫：10 eV 离子热化到气温、阈值以下电子按 2m/M·ν_el·ε 冷却 |
| T8 | restart/run.sh | 续算 | 6 周期直跑 vs 3 + 3 续算，history（除墙钟列）、剖面、能谱、x–t 逐位相同，2 rank，含碰撞 |
| T10 | langmuir/langmuir | 推进 + 场 | 冷 Langmuir 振荡：静止离子格点、电子正弦位移，模幅值以 (2/dt)·asin(ω_pe dt/2) 振荡且幅值不变 |
| T11 | `inputs/sheath.in` + tools/check_sheath | 无碰撞鞘层 | 能量守恒、记账恒等式、高斯定律、对称、Γ_e = Γ_i、T_⊥ 不变 |
| T12 | swarm/swarm + tools/boltzmann2term, compare_swarm | 碰撞的宏观效果 | 均匀电场、周期边界、脉冲 Townsend 电子群 vs 两项 Boltzmann 方程 |
| T13–T18 | tools/check_run, compare_turner, convergence, check_diag, dual_frequency, compare_table3 | 算例层面 | 见 §14 |

## 13. 输出文件格式

- **history.dat**（每周期一行，31 列）：1 period，2 t，3 N_e，4 N_i，5 ⟨n_e⟩，6 ⟨n_i⟩，7 n_e 中心，8 n_i 中心（中心 5 节点均值，周期末），9–12 Γ_i,L Γ_i,R Γ_e,L Γ_e,R，13–14 E_i,L E_i,R（平均撞击能），15 P_e，16 P_i（场对物种的精确做功），17 P_VI，18 电离数，19 激发数，20 KE_e，21 KE_i，22 W_field，23 累计壁面能量，24 累计做功，25 高斯残差最大值，26 墙钟，27 T_x,e，28 T_⊥,e，29 T_x,i，30 T_⊥,i，31 ν_max 超限数。恒等式：KE_e + KE_i − KE(0) = Work_cum − E_wall_cum。
- **profiles.dat**（时间平均，14 列）：x，n_e，n_i，φ，E，⟨ε_e⟩，电离率，⟨J_e·E⟩，⟨J_i·E⟩，⟨n_i v_i⟩，⟨ε_x,e⟩，⟨ε_i⟩，P_work,e，P_work,i。电极节点的密度、速率、J·E、通量、P_work 乘 2（半格体积）。
- **iedf_left/right.dat、eedf_left/right.dat**：能量格中心 [eV]、通量每 eV [m⁻²s⁻¹eV⁻¹]；头两行给总通量。left = x = 0 驱动电极。
- **eepf.dat**：中心 10 % 的 F(ε)（∫F dε = n_e）和 f_p = F/√ε。
- **xt_*.dat**（phi、ne、ni、Je、Ji、PeE、ioniz、Pe_work）：行 = 相位格，列 = 节点。
- **info.txt**：`key = value  # unit`，三段：参数与导出量、初始等离子体与碰撞参数、结果（中心密度、T_eff、T_x、ωpe·dt、λ_D/dx、四个通量、平均撞击能、能量平衡）。能量平衡的定义：P_part（精确做功）、P_wall（壁面动能通量）、P_dKE（储能变化）、P_inelastic（激发 + 电离阈值）、P_elastic_e、P_gas_i（离子交给气体）、P_newion（新离子热能）；`closure_residual = P_part + P_newion − walls − dKE − inelastic − elastic − P_gas_i` 应为舍入误差；P_JE（网格 J·E 积分）、P_VI（电源）、`P_numerical = P_VI − dW/dt − P_part`（格式能量误差）。
- **restart_rankNNN.bin、restart_grid.bin**：流式二进制，魔数 CCP1DR03。
- **init_particles_rankNNN.dat、init_grid.dat**：dump_initial 的 t = 0 状态。

## 14. 工具（tools/，只依赖 numpy、matplotlib）

| 工具 | 用法 | 做什么 |
| --- | --- | --- |
| `common.py` | 被其他工具 import | `load_dat`、`load_info`、`history_columns`、`profile_columns`、`check`、`check_range`、`expected_phase_drive`（按步精确平均的相位格驱动电压）、`turner_table3`（按 n_gas、V0 识别算例） |
| `plot_history.py` | `<run>` | 粒子数、密度、通量、能量、功率随周期；打印稳态漂移 |
| `plot_profiles.py` | `<run>` | 8 张时间平均剖面 |
| `plot_dist.py` | `<run>` | IEDF、壁面电子能谱、EEPF |
| `plot_xt.py` | `<run>` | 8 张 x–t 图叠加鞘层边缘（n_e/n_i = 0.5） |
| `plot_power.py` | `<run>` | 功率去向柱状图 |
| `plot_initial.py` | `<run> [out]` | 初始条件 5 张图（需 dump_initial） |
| `check_sheath.py` | `<run>` | T11 |
| `check_run.py` | `<run> [--full]` | T13：13 项自检（有限值、稳态漂移、相位格驱动电压、中心密度对 Turner、对称、通量平衡、离子能量、能量闭合、三把尺、ν_max、电离 = 损失、墙钟） |
| `compare_turner.py` | `<run> [...]` | T14：Case 1 剖面对 129 点参考曲线 |
| `convergence.py` | `<base> <dirs...>` | T15：加密算例与基线的差对比换种子噪声 |
| `check_diag.py` | `<run>` | T16：23 项诊断自洽 + diag_checks.png |
| `dual_frequency.py` | `<run0> <runs...>` | T17 |
| `compare_table3.py` | `<run> [--case N]` | T18：对 Turner 表 III |
| `compare_cases.py` | `<runs...> [--out]` | 多算例按 Turner 图 3–10 的单位并排 |
| `theory_check.py` | `<run>` | 经典公式对照（`doc/theory.md`） |
| `boltzmann2term.py` | `[dir] [--validate] [--en ...]` | 两项 Boltzmann 求解器（Scharfetter–Gummel，F = n/√ε，时间增长模型） |
| `compare_swarm.py` | `<dir>` | T12 比较 |
| `draw_flow.py` | `[out_dir]` | 功能流程图 |

## 15. 想改什么去哪里改

| 想做的事 | 改哪里 |
| --- | --- |
| 换气体 / 加一个碰撞过程 | `data/<gas>/` 放表；`xsec.f90` 的 `nproc_e`、`fname_e`、阈值逻辑；`mcc_electrons` 的分支与账本；`tests/mcc` 的期望值 |
| 各向异性散射 | `mcc.f90`：把 `isotropic_direction` 换成按微分截面抽 χ 并相对入射方向旋转 |
| 二次电子发射 | `push_species` 吸收分支：按 γ 概率 `add_particle(1, 壁面位置内侧, 发射速度)`，账本加发射能量 |
| 改边界电压波形（如非正弦、直流偏置、外电路） | `input.v_drive`；外电路需要在 `main` 的解场前解电路方程给 `v_left` |
| 新的剖面诊断 | `diag_step` 窗口内累加 + `diag_final` 的 profiles 列 + `tools/common.profile_columns` + 续算文件的读写列表（两处都要加，否则 T8 会失败） |
| 新的 history 列 | `diag_period` 的 buf 与写格式 + `history_open` 表头 + `tools/common.history_columns` |
| 改成二维 | 网格数组和沉积/插值全部要改；泊松改成二维求解器（em2 的共轭梯度思路）；本代码按一维设计，不建议硬改 |
| 性能 | 粒子循环已是 SoA、无分支的插值；进一步可做粒子按格排序（缓存）、OpenMP（需要每线程独立的随机数状态） |

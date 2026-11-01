# ccp1d 交接文件：在云端做出能用的一维 PIC-MCC

写给下一个 Claude 云端会话，也写给我自己。日期 2026-11-01。读完本文件再读 `ccp1d/PROGRESS.md`（进度），然后开工。

## 0. 目标

在 `ccp1d/` 下新写一个**一维有界静电 PIC-MCC** 程序（Fortran 90 + MPI，代码风格沿用 `pic/` 里的 em2），复现 Turner 等 2013 氦气射频容性放电基准的 Case 1，**输出的数据能直接用 Python 做分析和画图**。这是一个求职项目：代码要能跑、结果要能验证、图要能拿进面试。

验收分三档：

| 档 | 内容 | 判据 |
| --- | --- | --- |
| 一（必须） | 无碰撞鞘层；MCC 单元测试；Case 1 测试模式和全尺寸；收敛测试 | 时间平均离子密度剖面与 `ccp1d/data/turner2013_case1_ni.dat` 吻合：中心密度偏差 < 5%（目标 2%），形状一致；dt、dx、每格粒子数各变一倍结果不变 |
| 二 | 四类诊断：IEDF、EEPF、x–t 相位图、功率分解 | 每个量有图、单位正确、能用物理解释；能量平衡闭合 |
| 三 | 一个扩展算例：双频放电 | 低频电压升高时 IEDF 向高能展宽 |

em2（`pic/`）**不要动**，它是另一个项目；唯一的已知 bug 记在第 8 节。

## 1. 每次会话的开头和结尾

开头：

```sh
git clone <仓库地址> mccpic && cd mccpic     # 已有则 git pull
./setup.sh                                   # 装 gfortran + Open MPI，约 1 分钟；容器每次都是新的
cat HANDOVER.md ccp1d/PROGRESS.md
```

如果 GitHub 没接上，代码以 `mccpic-<日期>.tar.gz` 的形式附在会话里：`tar xzf` 解开即可，里面带 `.git` 历史。

结尾（不做等于白干——容器一段时间不活动就被回收，没 push 的东西全丢）：

```sh
# 更新 ccp1d/PROGRESS.md：做了什么、验证结果、下一步、未解决的问题
git add -A && git commit -m "..." && git push
```

云端约束（2026-11-01 实测）：Ubuntu 24.04，**2 核、7 GB 内存**；前台命令最长 10 分钟，后台最长 2 小时（跑完有通知）；`github.com` 可访问（clone、raw 文件），其他网站多数不通；apt、pip 可用。

## 2. 仓库现状

```
mccpic/
  README.md                    仓库总览
  HANDOVER.md                  本文件
  setup.sh                     装编译器
  pic/                         em2：二维电磁相对论 PIC（双 Harris 电流片磁重联），已整理、已在云端编译运行
  ccp1d/
    PROGRESS.md                进度（每次会话更新）
    data/SOURCES.md            数据来源和许可
    data/he/*.dat              氦截面表：电子弹性/激发×2/电离，He⁺ 各向同性/背向（energy[eV]  sigma[m^2]）
    data/turner2013_case1_ni.dat   Case 1 参考离子密度剖面，129 点，x[m] n_i[m^-3]
```

ccp1d 还没有一行代码。设计思路的长版在 Claude 项目文档 `claude/ccp1d_MCC改造思路.md`，em2 的整理记录在 `claude/pic代码整理与错误报告.md`；读不到也没关系，本文第 3–5 节自足。

em2 里可以直接搬的东西：`pic/dran2.f90`（随机数，原样复制）、`pic/move.f90` 的结构（先把 E 乘以 (q/m)·dt 再按粒子插值）、`pic/setrho.f90` 的线性权重、`pic/curnt.f90` 末尾的 `MPI_ALLREDUCE` 写法、`pic/output.f90`/`cont.f90` 的续算机制、`pic/input.f90` + `global.f90` 的模块组织、`pic/Makefile`。不能搬的：场求解（em2 推 Maxwell 方程，受 c·dt < dx 限制，Case 1 的 dt 比这个限制大 100 倍）、周期边界、固定粒子数的槽位布局、ωpe = c = 1 的单位。

## 3. 物理模型与参数

**几何与边界。** x ∈ [0, L]，L = 0.067 m。x = 0 的电极加电压 φ(0, t) = V0·sin(2π f t)，x = L 的电极接地 φ(L) = 0。粒子到达 x < 0 或 x > L 即被吸收删除（记录能量和种类），无二次电子、无反射。

**粒子。** 电子（m_e，−e）和 He⁺（M = 6.67×10⁻²⁷ kg，+e），1d3v：位置一维，速度三分量。背景氦气固定：密度 n_g，温度 T_g = 300 K，Maxwell 分布，不随时间变。所有宏粒子权重相同（单位电极面积上代表的真实粒子数 w，单位 m⁻²）。

**初始条件。** n_e = n_i = n0 均匀；电子 Maxwell 分布 T_e0 = 30000 K（2.585 eV），离子 T_i0 = 300 K；t = 0 时 V = 0 且准中性，E = 0。

**碰撞（基准规定）。** 电子–He：弹性、激发 19.82 eV、激发 20.61 eV、电离 24.59 eV；**所有电子过程各向同性散射**；电离后剩余能量**两个电子严格各分一半**，新离子取气体温度的 Maxwell 速度。He⁺–He：各向同性弹性 + 背向散射，都在质心系定义，靶原子速度每次从 T_g 的 Maxwell 分布抽；不另算电荷交换（背向散射就是它）。截面线性插值，超出表范围取最后一个值。常数用 CODATA（基准用 2006 版，差别可忽略）。

**四个算例：**

| | Case 1 | Case 2 | Case 3 | Case 4 |
| --- | --- | --- | --- | --- |
| V0 [V] | 450 | 200 | 150 | 120 |
| n_g [m⁻³] | 9.64e20 | 3.21e21 | 9.64e21 | 3.21e22 |
| n0 [m⁻³] | 2.56e14 | 5.12e14 | 5.12e14 | 3.84e14 |
| 网格数 nx | 128 | 256 | 512 | 512 |
| 每周期步数 | 400 | 800 | 1600 | 3200 |
| 初始每格粒子数（每种） | 512 | 256 | 128 | 64 |
| 总周期数 | 1280 | 5120 | 5120 | 15360 |
| 平均窗口 | 最后 32 周期 | 同 | 同 | 同 |

f = 13.56 MHz，T_g = 300 K，T_e0 = 30000 K 四例共用。云端只做 Case 1（2 核约 1 小时）；Case 2–4 要十几倍时间，留给集群。

**Case 1 的导出量（已算好，用来自检）：**

- dx = 0.5234 mm，dt = 1/(400 f) = 1.8437×10⁻¹⁰ s，总步数 512000
- 权重 w = n0·L/(nx·512) = 2.617×10⁸ m⁻²；初始每种 65536 个宏粒子；稳态时粒子数约为初始的 3–5 倍（中心密度升到 1.4×10¹⁴），maxnp 取初始的 8 倍
- 气压 3.99 Pa = 30 mTorr；电子平均自由程约 2 cm（和 6.7 cm 间距可比，所以是"低气压"算例，随机加热重要）
- 零碰撞：电子 ν_max = n_g·max[σ_tot(ε)·v] = 8.65×10⁷ s⁻¹（在 8.2 eV 取到），每步碰撞概率 P = 1 − e^(−ν_max dt) = **1.6%**；离子在 ε_rel ≤ 1000 eV 内 ν_max = 2.79×10⁷ s⁻¹，P = **0.5%**（离子的 σ·g 随能量单调增，ν_max 取决于能量上限，1000 eV 够用）
- 初始 ωpe·dt = 0.166，λ_D/dx = 1.43；Turner 表 III 用稳态值：ωpe·dt = 0.121，λ_D/dx = 3.72
- 参考结果：中心 n_i = 1.405×10¹⁴ m⁻³，电极处 3.6×10¹³ m⁻³；五个基准代码之间差别 < 1%

**测试模式**（WarpX 例子的做法，几分钟跑完，调试用）：nx = 32，每格 256 粒子，每周期 200 步（dt 加倍），320 周期。只看形状和量级。

## 4. 代码规格

### 4.1 单位与输入

内部全部 SI；能量只在查截面表和输出时换成 eV。参数用 **namelist** 文件 `ccp1d.in` 读入（不要像 em2 那样编译期常数：四个算例和收敛测试要换参数，不该重新编译），数组按 nx 和 maxnp 分配。建议的 namelist：

```fortran
&ccp1d
  L = 0.067d0, nx = 128
  freq = 13.56d6, V0 = 450d0
  freq_lf = 0d0, V0_lf = 0d0            ! 双频扩展：V(t) = V0 sin(2π f t) + V0_lf sin(2π f_lf t)
  n_gas = 9.64d20, T_gas = 300d0
  n0 = 2.56d14, Te0 = 30000d0, Ti0 = 300d0
  nppc = 512                            ! 初始每格每种宏粒子数
  steps_per_period = 400
  n_periods = 1280                      ! 本段要跑到的周期数（续算时是累计值）
  n_avg = 32                            ! 测量窗口：最后 n_avg 个周期
  n_phase = 100                         ! x–t 图每周期的相位格数
  seed = 12345
  xsec_dir = '../../data/he'
  restart = .false., restart_every = 160   ! 每 160 周期写一次续算文件
  iedf_emax = 600d0, iedf_de = 1d0         ! IEDF 直方图范围和分辨率 [eV]
  eepf_emax = 100d0, eepf_de = 0.5d0
/
```

一次运行放在自己的目录里（`ccp1d/runs/case1/`），`ccp1d.in` 放在该目录，可执行文件在上一级：`cd runs/case1 && mpirun --allow-run-as-root -np 2 ../../ccp1d`。

### 4.2 文件与职责

| 文件 | 职责 | 对应 em2 |
| --- | --- | --- |
| `input.f90` | namelist 读入、参数检查、导出量（dx、dt、w、ν_max、ωpe·dt、λ_D/dx）；rank 0 打印 | input |
| `global.f90` | 常数、网格数组 `phi, ex, rho, rho_e, rho_i, jx_e, jx_i (0:nx)`、粒子数组、计数器 | global |
| `main.f90` | MPI、新算/续算分支、主循环、计时 | main |
| `init.f90` | 装粒子（位置均匀随机，速度 Box–Muller 抽 Maxwell）；装截面；E = 0；蛙跳起步 | creator+smear+setv |
| `move.f90` | 插值 E、推速度和位置、电极吸收 | move + curnt 的位置部分 |
| `deposit.f90` | ρ（测量期间还有 ρ_e、ρ_i、J_e、J_i）的线性权重沉积 + `MPI_ALLREDUCE` | curnt 的 ρ 部分、setrho |
| `poisson.f90` | 三对角求解 φ，E = −dφ/dx | 全新（em2 的 cgm 是周期二维修正，用不上） |
| `xsec.f90` | 读 7 个表（电荷交换不用），插值到均匀能量网格，算 ν_max | 全新 |
| `mcc.f90` | 零碰撞法、各过程的碰后速度、电离产生新粒子、能量账本 | 全新 |
| `diag.f90` | 时间平均、直方图、x–t 累加、history、info、续算文件（只 rank 0 写） | output + cont |
| `dran2.f90` | 随机数 | 原样复制 |
| `Makefile` | `make`（-O3）、`make debug`（见 7）、依赖写全 | Makefile |

### 4.3 数据结构

每个进程只存自己的粒子：`real(8), allocatable :: x(:,:), vx(:,:), vy(:,:), vz(:,:)`，维度 `(maxnp, nsp)`，`integer(8) :: np(nsp)`。位置用米。删除粒子：把最后一个搬到被删的位置，`np = np − 1`（顺序无关）；新增：追加到末尾；`np > maxnp` 时报错退出并提示增大 maxnp。所有进程都持有完整网格（nx+1 个节点，很小），不做区域分解。

### 4.4 算法

**时间推进（蛙跳）。** 第 n 步开始时有 φ、E 在 t_n（电压 V(t_n)），速度在 t_{n−1/2}，位置在 t_n：

1. `move`：E_p = 线性插值 E 到粒子位置；v_x += (q/m)·E_p·dt；x += v_x·dt；v_y、v_z 不变。x < 0 或 x > L：记录 ½m(v_x²+v_y²+v_z²)（IEDF/EEDF 直方图）、计数（通量）、电荷（电极电流），删除。
2. `mcc`：电子、离子各做一遍零碰撞判断，改速度，电离产生新粒子（新电子、新离子放在父粒子位置）。
3. `deposit`：ρ(t_{n+1})，各进程求和。
4. `poisson`：φ(0) = V(t_{n+1})，φ(L) = 0，解 φ 和 E。
5. `diag`：进入测量窗口后累加。

初始 E = 0，所以 −dt/2 的速度回推是零，但仍按通用写法写一次（以后改初值不用动）。

**沉积（CIC）。** j = floor(x/dx)，s = x/dx − j：`rho(j) += q·w·(1−s)/dx`，`rho(j+1) += q·w·s/dx`（单位 C/m³）。电极节点 0 和 nx 的 ρ 不参与泊松方程（电势固定），不用处理半格体积。J 同样沉积 q·w·v_x·权重/dx。

**泊松。** 节点 j = 0..nx，x_j = j·dx。φ_{j−1} − 2φ_j + φ_{j+1} = −ρ_j dx²/ε0，j = 1..nx−1，Thomas 算法（系数每步相同，只换右端）。E_j = −(φ_{j+1} − φ_{j−1})/(2dx)；两端二阶单侧差分 E_0 = −(−3φ_0 + 4φ_1 − φ_2)/(2dx)，E_nx 对称。

**零碰撞法（每种粒子）。** 启动时在 0–1000 eV 的能量网格上算 ν_max = n_g·max[σ_tot(ε)·v(ε)]（电子：v = √(2ε/m_e)；离子：用 ε_rel 和 g = √(2ε_rel/μ)，μ = M/2）。每步每个粒子：抽 R1，R1 < P = 1 − e^(−ν_max dt) 才进入碰撞；算当前能量的各过程频率 ν_k = n_g σ_k v，抽 R2·ν_max，落在累积区间 [0, ν_1)、[ν_1, ν_1+ν_2)… 的哪一段就发生哪个过程，落在 Σν_k 之外是空碰撞。截面查表：启动时把表插值到 Δε = 0.01 eV 的均匀网格，查表 `k = int(ε/Δε)`，O(1)。

**电子碰后速度**（各向同性：直接在球面均匀抽新方向 cosχ = 1 − 2R，φ = 2πR'，不需要相对入射方向旋转）：

- 弹性：ε' = ε·[1 − (2m_e/M)(1 − cosχ)]，χ 为新旧方向夹角（Vahedi–Surendra；基准没写但这是 19.8 eV 以下电子唯一的能量损失渠道，对剖面影响很小）
- 激发：ε' = ε − ε_th
- 电离：ε' = (ε − 24.59 eV)/2 给原电子和新电子各一份，各自独立抽方向；新离子速度从 T_g 的 Maxwell 分布抽（每分量 σ_v = √(k_B T_g/M) = 788 m/s）

**离子碰后速度。** 抽靶原子速度 v_n（三个高斯分量，σ_v = 788 m/s）；g = v_i − v_n；ε_rel = ½μg²；查 σ_iso(ε_rel)、σ_back(ε_rel)，ν_k = n_g σ_k |g|。v_cm = (v_i + v_n)/2。各向同性：质心系均匀抽新方向 n̂，v_i' = v_cm + (|g|/2)·n̂。背向：n̂ = −ĝ，代入得 **v_i' = v_n**——离子换成中性原子的速度，一行代码。

**随机数。** `dran2`，每个进程种子 `idum = −(seed + MYRANK)`，只在启动时初始化一次；续算文件里存随机数状态（或接受续算后序列不同，在 info 里注明）。`dran2` 有 `SAVE` 状态，不能多线程调用：只用 MPI，不开 OpenMP。

**MPI。** 粒子在初始化时按进程均分（每个进程各装 N/NPROCS 个，位置都在整个 [0, L] 上均匀分布）；每步一次 `MPI_ALLREDUCE` ρ（nx+1 个数，测量期间再加 ρ_e、ρ_i、J_e、J_i）；直方图和计数器只在输出时归约。电离和壁面损失都是随机的，负载自然均衡。结果只在进程数固定时可重复。

**续算。** 每 `restart_every` 个周期和运行结束时写 `restart_rankNNN.bin`（np、x、v、随机数状态）和 `restart_grid.bin`（步数、时间、φ、测量窗口的累加量）；`restart = .true.` 时读入继续。Case 1 全尺寸在云端要分段跑（见 6），续算是必需品。

## 5. 输出数据规格（重点：要能直接分析）

原则：**全部 ASCII，`#` 开头为注释头，`numpy.loadtxt` 直接能读，列名和单位写在头里**；只有 rank 0 写文件；时间平均量只在测量窗口（最后 n_avg 个周期）累加；文件名固定、覆盖式；每次运行一个目录，`ccp1d.in` 和输出放在一起。x–t 矩阵和续算文件较大，不进 git；关键算例的小文件（info、profiles、iedf、eepf、history）复制到 `ccp1d/results/<名字>/` 进 git。

| 文件 | 何时写 | 内容 | 用途 |
| --- | --- | --- | --- |
| `info.txt` | 开始时写参数，结束时追加结果 | `key = value  # unit` 格式。输入参数全部；导出量 dx、dt、w、maxnp、ν_max(e)、ν_max(i)、P_e、P_i、初始 ωpe·dt 和 λ_D/dx；结束时：稳态中心 n_e、ωpe·dt 和 λ_D/dx（用稳态值）、最终粒子数、每电极的平均离子通量和平均撞击能量、各种粒子吸收功率、能量平衡各项（见下）、每周期电离数、墙钟时间、git commit、日期 | 一眼看出参数对不对；面试题 4、11 的数字 |
| `history.dat` | 每周期一行 | 周期号、t[s]、N_e、N_i（宏粒子）、⟨n_e⟩、⟨n_i⟩[m⁻³]、中心 n_i、离子通量 左/右[m⁻²s⁻¹]、电子通量 左/右、离子平均撞击能量 左/右[eV]、P_e、P_i[W/m²]、本周期电离数、激发数 | 看是否到稳态（Case 1 要 1280 周期才算），续算衔接是否平滑 |
| `profiles.dat` | 结束 | x[m]、n_e、n_i[m⁻³]、φ[V]、E[V/m]、电子平均能量[eV]、电离率[m⁻³s⁻¹]、⟨J_e·E⟩、⟨J_i·E⟩[W/m³]；时间平均；nx+1 行 | **和 Turner 对比**（n_i 列）；鞘层厚度；加热在哪发生 |
| `iedf_left.dat`、`iedf_right.dat` | 结束 | ε_bin 中心[eV]、通量能量分布[m⁻²s⁻¹eV⁻¹]（撞墙宏粒子计数 × w /(Δε·T_avg)） | 双峰结构、峰间距（面试题 8） |
| `eedf_left.dat`、`eedf_right.dat` | 结束 | 同上，电子 | 可选 |
| `eepf.dat` | 结束 | ε[eV]、f_p(ε) = F(ε)/√ε [eV⁻³ᐟ² m⁻³]，取中心 10% 区域，时间平均；归一化 ∫F dε = n_e | 半对数图上 Maxwell 是直线、斜率 −1/T_e；双温结构（面试题 7） |
| `xt_phi.dat`、`xt_ne.dat`、`xt_ni.dat`、`xt_Je.dat`、`xt_Ji.dat`、`xt_PeE.dat`、`xt_ioniz.dat` | 结束 | 矩阵 n_phase 行 × (nx+1) 列，头里写 nx、n_phase、dx、周期；每行是一个相位格在测量窗口内的平均 | 鞘层振荡、随机加热发生在鞘层边缘、电离在哪发生（面试题 6） |
| `restart_*.bin` | 每 restart_every 周期 | 二进制，格式写在 info.txt 和 README | 分段跑 |

**能量平衡（最硬的自检，写进 info.txt 末尾）。** 测量窗口内每单位电极面积、每单位时间：吸收功率 P_abs = ∫(⟨J_e·E⟩ + ⟨J_i·E⟩)dx 应等于 粒子带到两电极的动能通量 + 非弹性损失（激发数×19.82/20.61 eV + 电离数×24.59 eV）+ 电子弹性损失（Σ Δε）+ 储能变化（稳态 ≈ 0）。各项分别列出，再给闭合误差百分比。闭合到几个百分点以内才算对。

**分析脚本（`ccp1d/tools/`，Python，只依赖 numpy、matplotlib）：**

| 脚本 | 读 | 出 |
| --- | --- | --- |
| `compare_turner.py <run_dir>` | profiles.dat + data/turner2013_case1_ni.dat | n_i(x) 两条曲线叠画 + 相对误差子图；打印中心密度、最大相对偏差、均方偏差。网格不同于 129 点时把参考值插到本网格（照 WarpX 的 analysis_1d.py 做法） |
| `plot_history.py <run_dir>` | history.dat | 粒子数、中心密度、功率 vs 周期；判断稳态 |
| `plot_profiles.py <run_dir>` | profiles.dat | n_e/n_i/φ/E/⟨ε_e⟩/电离率/功率 的多子图 |
| `plot_dist.py <run_dir>` | iedf_*.dat、eepf.dat | IEDF（线性）和 EEPF（半对数） |
| `plot_xt.py <run_dir>` | xt_*.dat | 每个量一张 x–t 伪彩图，横轴 x[cm]、纵轴相位 t/T |
| `convergence.py <dir1> <dir2> ...` | 多个 profiles.dat | 剖面叠画 + 两两偏差表 |

每个脚本输出 PNG 到 run 目录，并把关键数字打印到终端。第一次能画出 `compare_turner.py` 的图，就是第一档的里程碑。

## 6. 任务顺序与验收

按顺序做，每步有便宜的验收，做完一步 commit 一次（提交信息里写验收结果）。

| # | 任务 | 验收 | 预计 |
| --- | --- | --- | --- |
| M0 | 骨架：Makefile、input、global、main、init、deposit、poisson、move、diag（只写 info 和 history）、续算 | `tests/poisson/`：ρ = 常数 → φ 为抛物线（与解析解比到 1e-10）；ρ = 0 → 直线；V(t) 正确 | 1 个会话 |
| M1 | 无碰撞鞘层：n_g = 0、V0 = 0，均匀等离子体两壁吸收，跑 2 µs | 鞘层形成；中心相对壁面电势 ≈ (T_e/2)·ln(M/(2π m_e)) ≈ 3.5 T_e（氦）；电子先快速损失后准稳态；无 NaN；能量账对得上 | 同一会话 |
| M2 | xsec + mcc + 单元测试 `tests/mcc/`（不解场，E = 0） | (a) 电子全给固定能量 ε0 = 10、30、100 eV，统计各过程的碰撞次数/时间，与 n_g σ_k(ε0) v 比，统计误差内一致，且 Σν_k < ν_max；(b) 每种过程碰前后能量守恒（电离：ε = ε1 + ε2 + 24.59）；(c) 静止离子在 300 K 气体里应热化到 300 K，10 eV 离子应减速 | 1 个会话 |
| M3 | 电子群测试：均匀 E、无壁（周期或再注入）、只算电子、关电离产物 | E/N = 50、100、200、500 Td 下的平均能量和漂移速度；对照值用 LXCat Biagi 数据库的 swarm 表或 Mac 上的 BOLSIG+（云端拿不到就先存结果） | 半个会话 |
| M4 | Case 1 测试模式 → 全尺寸 | 测试模式：形状对、量级对（几分钟）。全尺寸：`compare_turner.py` 中心密度偏差 < 5%，形状吻合；1280 周期分 8 段（每段 160 周期，预计 8 分钟左右，先跑 10 个周期按墙钟时间定），后台跑，每段结束 commit history/profiles | 1 个会话 + 后台 |
| M5 | 收敛：dt/2、nx = 256、nppc 256 和 1024 各一次 | `convergence.py` 表：剖面差别在统计噪声内 | 后台 |
| M6 | 诊断四件套 + 能量平衡 + 图 | 四类图齐全，能量闭合 < 5% | 1 个会话 |
| M7 | 双频扩展 | V_lf 升高 IEDF 展宽、密度变化相对小 | 1 个会话 |
| 文档 | `ccp1d/README.md`：一条命令复现、Case 1 对比图放顶部、输出文件说明、与 em2 的关系 | 别人按 README 能跑通 | 半个会话 |

时间估算（2 核）：Case 1 全尺寸每步约 30–40 万粒子 × （推进 + MCC + 沉积），大约 1 小时；测试模式 1–3 分钟。第一次全尺寸前先用 `history.dat` 看 10 个周期的墙钟时间再排分段。

## 7. 工程约定与已知坑

- 编译：`make` 用 `-O3 -march=native`；`make debug` 用 `-g -O0 -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -Wall -finit-local-zero`。**新写的每个例程先在 debug 模式下跑过测试再开 -O3**（em2 的越界 bug 就是 -O3 下静默跑了多年，debug 模式一秒定位，见第 8 节）。Makefile 把模块依赖写全（`main.o: global.o input.o` 之类），em2 的没写全。
- 运行：`mpirun --allow-run-as-root -np 2 ../../ccp1d`；进程数 > 核数要加 `--oversubscribe`。
- gfortran 的 stdout 在管道里是全缓冲的：进度写到 `run.log` 并 `flush`，或设 `GFORTRAN_UNBUFFERED_ALL=1`。
- 长算例用后台运行（≤ 1.5 小时一段），每段结束必须有续算文件；不要指望一个会话跑完 1280 周期。
- 粒子计数和总数用 `integer(8)`。
- 不做任何滤波（基准代码都不平滑）。
- 每个物理量输出带单位；每张图的轴带单位。
- git：一个功能一个 commit；`runs/` 整个目录在 `.gitignore`，关键结果拷到 `results/`。
- 不确定基准的某条规定时，查 Turner 2013 原文（arXiv:1211.5246）和 WarpX 的 `capacitive_discharge` 例子（它已对照过 Turner），不要猜。

## 8. em2 的一个待修 bug（不影响 ccp1d，顺手记下）

2026-11-01 在云端用缩小参数（nx = 128、ny = 256、npcel = 4、nproc = 2）跑 em2，rank 0 卡死不退出，文件名里出现 `*****`。`-fcheck=all` 立刻定位：`setv.f90:83` 下标 259 超出上界 258。根源在 `pic/creator.f90` 装背景格点的循环：

```fortran
do i=1,npback/4
  x(MYRANK*maxnp+i,isp)=mod((i-0.5d0)*ddx,xn/NPROCS) + xn*MYRANK/NPROCS
  y(MYRANK*maxnp+i,isp)=ddy*(i/nnx+0.5d0)      ! 没有检查 y 是否超过 yn/2
enddo
```

`nnx = int(xn/ddx/NPROCS)` 向下取整，每行粒子变少、行数变多，格点顶端超出 `yn/2`；上半区复制时再加 `yn/2`，粒子就落到盒子外。`setv`/`move` 越界读，`curnt` 越界写，把 rank 0 的全局变量 `it`、`time` 写坏（文件名 `*****`），循环永远到不了 `it == nt`，在 `MPI_ALLREDUCE` 里空转。原始参数（1024/0.2/4 = 1280 恰为整数）只超出 0.1 格，落在保护格里被周期折回吃掉，所以多年没暴露。修法：背景循环里加一行

```fortran
if(y(MYRANK*maxnp+i,isp) >= 0.5d0*yn) y(MYRANK*maxnp+i,isp) = y(MYRANK*maxnp+i,isp) - 0.5d0*yn
```

（电流片那两段已有同样的保护。）改完要在 `pic/README.md` 和项目文档的错误报告里补一条。

## 9. 面试材料（做完 M6 后整理）

| 图 / 数字 | 来自 | 回答的问题 |
| --- | --- | --- |
| n_i(x) 与 Turner 叠画 + 误差 | compare_turner.py | 怎么证明代码是对的（题 12） |
| 收敛表 | convergence.py | 数值条件：ωpe·dt、λ_D/dx、每格粒子数（题 11） |
| ν_max、P_e = 1.6%、P_i = 0.5% | info.txt | 零碰撞法、时间步对碰撞的限制（题 1、4） |
| IEDF 双峰 | plot_dist.py | 峰间距与射频周期、离子渡越时间（题 8） |
| EEPF 半对数 | plot_dist.py | 低气压双温结构（题 7） |
| x–t 的 φ 和 J_e·E | plot_xt.py | 鞘层振荡、随机加热 vs 欧姆加热（题 5、6） |
| 能量平衡表 | info.txt | 代码自检 |
| 双频 IEDF 对比 | M7 | 为什么工业上用双频（题 10） |
| `mcc.f90` 里弹性、激发、电离、离子碰撞四段代码 | — | 题 2、3 |

仓库 README 顶部放第一张图；5–8 页幻灯片按这张表的顺序讲。

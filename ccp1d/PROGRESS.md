# ccp1d 进度

每次会话结束前在最上面加一段：日期、做了什么、验证结果、下一步、未解决的问题。

## 2026-11-01（第三次会话：M0 + M1 + 初始条件讲解）

- **M0 完成**：`ccp1d/` 下 10 个源文件（input/global/random/poisson/deposit/push/diag/init/main + 生成的 version）、Makefile（`make`/`make debug`/`make check`）、namelist 模板 `inputs/`。collisionless 全流程能跑：装载、CIC 沉积、三对角泊松、蛙跳推进、电极吸收、history/profiles/IEDF/EEDF/EEPF/x–t 输出、续算文件、`periods_this_run` 和 `max_wall_minutes` 分段。
- **单元测试全过**（`make check`，debug 编译）：T1 泊松（解析解 1e-14、二阶收敛比 4.0004、高斯定律 1e-13）；T2 沉积（权重、电荷守恒、端点半格、矩、1 与 2 rank 一致）；T3 随机数；T8 续算 20 周期 vs 10+10 **逐位相同**（history、profiles、直方图、x–t 全部 cmp 通过）；T10 Langmuir 振荡频率与蛙跳理论值差 3e-4、与 ωpe 差 8.6e-4、振幅 1500 步不增长。
- **M1 无碰撞鞘层**（`inputs/sheath.in`，10 µs，2 核 21 s）：总能量守恒到 5e-4（80 % 初始能量进壁）、记账恒等式 1e-8、高斯定律 2e-13、两半离子数差 0.14 %、每壁 Γ_e = Γ_i 到 1.6 %、T_⊥,e 不变（0.5 %）。物理发现：一维无碰撞等离子体的电子 T_x 从 2.57 eV 塌到 0.13 eV 而 T_⊥ 不变，浮动电势只有 0.19 V——教科书 3.5 T_e 不适用，`check_sheath.py` 据此改为用 T_x 的宽松判据。结果在 `results/sheath/`。
- **能量三把尺已实现**：P_part（逐粒子做功，恒等式闭合 1e-16）、P_JE（网格 J·E）、P_VI（Shockley–Ramo 电流 I = (w/L)Σqv_x 的 ⟨VI⟩）；`P_numerical = P_VI − dW/dt − P_part` 是格式能量误差（烟雾测试里 0.9 %）。
- **初始条件讲解**：`doc/initial_conditions.md` + `doc/fig/init_*.png`（相空间、速度/能量分布、t=0 网格量、驱动与尺度、截面与 ν_max）+ 鞘层图；`tools/plot_initial.py` 从 `dump_initial = .true.` 的输出画图。实测随机装载的初始 φ 噪声 rms 5.7 V（长波被泊松放大），E rms 630 V/m。
- 与 PLAN 的差异：history 列增加到 31（加 P_VI、T_x/T_⊥ 四列）；profiles 增加 ⟨ε_x,e⟩ 和 ⟨ε_i⟩ 两列；x–t 的 7 个文件、EEPF、IEDF 已在 M0 实现（M6 只剩画图和自洽检查）；`deposit` 的矩例程返回 rank 局部和，窗口结束时一次归约。
- 性能：sheath 算例（2 × 65536 粒子起，54400 步）2 核 21 s；init_case1 一个周期 1.7 s → Case 1 全尺寸 1280 周期估计 **20–40 分钟**（无碰撞时），MCC 会再加一些。
- 下一步：M2（`xsec.f90`、`mcc.f90`、`tests/mcc/` T4–T7），然后 Case 1 测试模式 → 全尺寸。
- 未解决：em2 的 creator 越界 bug 未修（HANDOVER §8）；Langmuir 测试的冷等离子体有网格加热（能量只报告不判定，见 tests/langmuir）。

## 2026-11-01（第二次会话：计划）

- 核对了 HANDOVER、em2 源码、7 个截面表和 Turner 参考剖面，写出 `ccp1d/PLAN.md`：实施计划（模块、接口、每步顺序、namelist 新增项、续算格式、Makefile）和五层测试方案（T1–T17，每项有数字判据和命令）。
- 对交接文件的修正（详见 PLAN §1）：稳态宏粒子数应降到初始的约 30%（参考剖面平均 0.29 n0），不是升高 3–5 倍；电极节点密度要乘 2；dran2 改成可存取状态的模块以便续算逐位复现；截面查表改线性插值（热离子段 σ 变化 2 倍）；功率用 P_part / P_JE / P_VI 三把尺；M1 要跑 10 µs（离子渡越 4.3 µs）；MCC 频率测试的期望值乘 0.992。
- 重算 ν_max：电子 8.654e7 s⁻¹（8.24 eV）、离子 2.793e7 s⁻¹，与 HANDOVER 一致。
- 本次容器 4 核 15 GB（与交接写的 2 核 7 GB 不同，每次要看）；gfortran 未装（需 `./setup.sh`）；numpy 2.x 无 `np.trapz`。
- ccp1d 代码：仍未开始。下一步按 PLAN §5 做 M0 + M1。
- 未解决：em2 的 creator 越界 bug 未修（HANDOVER §8）。

## 2026-11-01

- 仓库建好：`pic/`（em2，已整理）、`ccp1d/data/`（氦截面 7 个表 + Turner Case 1 参考剖面）、`setup.sh`、`HANDOVER.md`。
- 云端环境实测：gfortran 13.3 + Open MPI 4.1.6 装得上，em2 编译运行正常；2 核 7 GB。
- 从截面表算出 Case 1 的 ν_max：电子 8.65e7 s⁻¹（P = 1.6%/步），离子 2.79e7 s⁻¹（P = 0.5%/步），见 HANDOVER §3。
- ccp1d 代码：**尚未开始**。下一步 M0（骨架 + 泊松测试），见 HANDOVER §6。
- 远端仓库 https://github.com/fanyumeng666/maccpic 已建，全部提交已推送（2026-11-01）；Claude GitHub App 已对该仓库授权，云端会话可直接 clone/push。
- 未解决：em2 的 creator 越界 bug 未修（HANDOVER §8）。

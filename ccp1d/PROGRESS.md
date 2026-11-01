# ccp1d 进度

每次会话结束前在最上面加一段：日期、做了什么、验证结果、下一步、未解决的问题。

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

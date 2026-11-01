# ccp1d 进度

每次会话结束前在最上面加一段：日期、做了什么、验证结果、下一步、未解决的问题。

## 2026-11-01

- 仓库建好：`pic/`（em2，已整理）、`ccp1d/data/`（氦截面 7 个表 + Turner Case 1 参考剖面）、`setup.sh`、`HANDOVER.md`。
- 云端环境实测：gfortran 13.3 + Open MPI 4.1.6 装得上，em2 编译运行正常；2 核 7 GB。
- 从截面表算出 Case 1 的 ν_max：电子 8.65e7 s⁻¹（P = 1.6%/步），离子 2.79e7 s⁻¹（P = 0.5%/步），见 HANDOVER §3。
- ccp1d 代码：**尚未开始**。下一步 M0（骨架 + 泊松测试），见 HANDOVER §6。
- 远端仓库 https://github.com/fanyumeng666/maccpic 已建；本地历史已接在其初始提交之后。
- 未解决：em2 的 creator 越界 bug 未修（HANDOVER §8）。

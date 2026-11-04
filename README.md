# mccpic

Particle-in-cell plasma simulation codes (Fortran 90, MPI).

| Directory | What it is | Status |
| --- | --- | --- |
| `pic/` | **em2** — two-dimensional, relativistic, electromagnetic PIC code for magnetic reconnection in a double Harris current sheet. Original author S. Matsukiyo (2012); cleaned up and corrected in October 2026. | working |
| `ccp1d/` | one-dimensional bounded electrostatic PIC/MCC code for a radio-frequency capacitive discharge (Fortran 90 + MPI, namelist input, ASCII output). Reproduces Case 1 of the Turner et al. (2013) helium benchmark to +2.1 % at the centre and 1.5 % rms over the gap, with convergence, diagnostic self-consistency, electron-swarm and dual-frequency studies. Entry points: [`ccp1d/README.md`](ccp1d/README.md), the figure guide [`ccp1d/doc/figures.md`](ccp1d/doc/figures.md), the module guide [`ccp1d/doc/program_flow.md`](ccp1d/doc/program_flow.md); specification in [`HANDOVER.md`](HANDOVER.md), plan and test matrix in [`ccp1d/PLAN.md`](ccp1d/PLAN.md), log in [`ccp1d/PROGRESS.md`](ccp1d/PROGRESS.md) | working, milestones M0-M7 done |

Documentation for em2 (Chinese / English / Japanese): [`pic/README.md`](pic/README.md) (parameters, build, run, output files) and [`pic/RECONNECTION_MODEL.md`](pic/RECONNECTION_MODEL.md) (how the reconnection initial state is set up).

Authorship and license: `ccp1d/` and the repository-level files are by Yumeng Fan (MIT, see `LICENSE`). `pic/` (em2) is by S. Matsukiyo and is not covered by the MIT license; it is included as the basis of `ccp1d/`. Cross-section data: CDLA-Permissive-2.0, see `ccp1d/data/SOURCES.md`.

Fresh container: `./setup.sh` installs gfortran, Open MPI and the Python plotting libraries.

Build, test and run ccp1d:

```sh
cd ccp1d && make && make check                      # -O3 executable + unit tests (about 2 min)
mkdir -p runs/case1 && cp inputs/case1.in runs/case1/ccp1d.in
cd runs/case1 && mpirun --allow-run-as-root -np 4 ../../ccp1d && python3 ../../tools/compare_turner.py .
```

Build and run em2:

```sh
cd pic
make clean && make          # mpif90; edit input.f90 first (nx, ny, nproc, ...)
mpirun -np 4 ./em2          # -np must equal nproc in input.f90
```

---

等离子体 PIC 模拟代码（Fortran 90，MPI 并行）。

- `pic/`：**em2**，二维相对论性电磁 PIC 代码，模拟双 Harris 电流片的磁重联。原作者 S. Matsukiyo（2012），2026 年 10 月整理、修正。说明见 `pic/README.md`（参数、编译、运行、输出）和 `pic/RECONNECTION_MODEL.md`（重联初始状态的实现），中、英、日三语。
- `ccp1d/`（已完成 M0–M7）：一维有界静电 PIC/MCC 代码，模拟射频容性放电。复现 Turner 等（2013）氦气基准 Case 1：中心离子密度偏差 +2.1 %，129 点均方偏差 1.5 %；另有收敛测试、诊断自洽检查、电子群对比两项 Boltzmann 解、双频扩展、Case 2 试算。入口：`ccp1d/README.md`（编译、测试、运行）、`ccp1d/doc/figures.md`（图册与数据索引）、`ccp1d/doc/program_flow.md`（逐模块讲解）、`HANDOVER.md`（规格）、`ccp1d/PLAN.md`（计划与测试矩阵）、`ccp1d/PROGRESS.md`（日志）。新容器先运行 `./setup.sh` 装编译器。
- 作者与许可：`ccp1d/` 和仓库级文件由 Fan Yumeng 编写（MIT，见 `LICENSE`）；`pic/`（em2）作者为 S. Matsukiyo，不在 MIT 范围内，作为 `ccp1d/` 的基础收录。截面数据为 CDLA-Permissive-2.0，见 `ccp1d/data/SOURCES.md`。

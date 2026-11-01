# mccpic

Particle-in-cell plasma simulation codes (Fortran 90, MPI).

| Directory | What it is | Status |
| --- | --- | --- |
| `pic/` | **em2** — two-dimensional, relativistic, electromagnetic PIC code for magnetic reconnection in a double Harris current sheet. Original author S. Matsukiyo (2012); cleaned up and corrected in October 2026. | working |
| `ccp1d/` | one-dimensional bounded electrostatic PIC/MCC code for a radio-frequency capacitive discharge, built on the particle routines of `pic/`; first target is the Turner et al. (2013) helium benchmark. Cross-section tables and the Case 1 reference profile are in `ccp1d/data/`; the specification and work plan are in [`HANDOVER.md`](HANDOVER.md) | in preparation |

Documentation for em2 (Chinese / English / Japanese): [`pic/README.md`](pic/README.md) (parameters, build, run, output files) and [`pic/RECONNECTION_MODEL.md`](pic/RECONNECTION_MODEL.md) (how the reconnection initial state is set up).

Authorship and license: `ccp1d/` and the repository-level files are by Yumeng Fan (MIT, see `LICENSE`). `pic/` (em2) is by S. Matsukiyo and is not covered by the MIT license; it is included as the basis of `ccp1d/`. Cross-section data: CDLA-Permissive-2.0, see `ccp1d/data/SOURCES.md`.

Fresh container: `./setup.sh` installs gfortran, Open MPI and the Python plotting libraries.

Build and run em2:

```sh
cd pic
make clean && make          # mpif90; edit input.f90 first (nx, ny, nproc, ...)
mpirun -np 4 ./em2          # -np must equal nproc in input.f90
```

---

等离子体 PIC 模拟代码（Fortran 90，MPI 并行）。

- `pic/`：**em2**，二维相对论性电磁 PIC 代码，模拟双 Harris 电流片的磁重联。原作者 S. Matsukiyo（2012），2026 年 10 月整理、修正。说明见 `pic/README.md`（参数、编译、运行、输出）和 `pic/RECONNECTION_MODEL.md`（重联初始状态的实现），中、英、日三语。
- `ccp1d/`（准备中）：一维有界静电 PIC/MCC 代码，模拟射频容性放电，复用 `pic/` 的粒子推进等例程；第一个目标是复现 Turner 等（2013）的氦气基准算例。截面数据和参考曲线在 `ccp1d/data/`，设计规格、输出格式和任务顺序见 `HANDOVER.md`，进度见 `ccp1d/PROGRESS.md`。新容器先运行 `./setup.sh` 装编译器。
- 作者与许可：`ccp1d/` 和仓库级文件由 Fan Yumeng 编写（MIT，见 `LICENSE`）；`pic/`（em2）作者为 S. Matsukiyo，不在 MIT 范围内，作为 `ccp1d/` 的基础收录。截面数据为 CDLA-Permissive-2.0，见 `ccp1d/data/SOURCES.md`。

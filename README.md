# mccpic

Particle-in-cell plasma simulation codes (Fortran 90, MPI).

| Directory | What it is | Status |
| --- | --- | --- |
| `pic/` | **em2** — two-dimensional, relativistic, electromagnetic PIC code for magnetic reconnection in a double Harris current sheet. Original author S. Matsukiyo (2012); cleaned up and corrected in October 2026. | working |
| `ccp1d/` | one-dimensional bounded electrostatic PIC/MCC code for a radio-frequency capacitive discharge, built on the particle routines of `pic/`; first target is the Turner et al. (2013) helium benchmark | planned |

Documentation for em2 (Chinese / English / Japanese): [`pic/README.md`](pic/README.md) (parameters, build, run, output files) and [`pic/RECONNECTION_MODEL.md`](pic/RECONNECTION_MODEL.md) (how the reconnection initial state is set up).

Build and run em2:

```sh
cd pic
make clean && make          # mpif90; edit input.f90 first (nx, ny, nproc, ...)
mpirun -np 4 ./em2          # -np must equal nproc in input.f90
```

---

等离子体 PIC 模拟代码（Fortran 90，MPI 并行）。

- `pic/`：**em2**，二维相对论性电磁 PIC 代码，模拟双 Harris 电流片的磁重联。原作者 S. Matsukiyo（2012），2026 年 10 月整理、修正。说明见 `pic/README.md`（参数、编译、运行、输出）和 `pic/RECONNECTION_MODEL.md`（重联初始状态的实现），中、英、日三语。
- `ccp1d/`（计划中）：一维有界静电 PIC/MCC 代码，模拟射频容性放电，复用 `pic/` 的粒子推进等例程；第一个目标是复现 Turner 等（2013）的氦气基准算例。

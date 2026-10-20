# em2 first step memo (2026/10/20)

[中文](#中文) | [English](#english) | [日本語](#日本語)

Adapted from S. Matsukiyo, "em1d first step memo". Only the parts that differ for em2 have been changed.
The reconnection setup itself is explained in [RECONNECTION_MODEL.md](RECONNECTION_MODEL.md).

---

## 中文

### em2 入门说明（2026/10/20）

本文简要说明二维周期边界全粒子（PIC，Particle-In-Cell）模拟代码 em2 的用法，根据 S. Matsukiyo 老师的《em1d first step memo》改写，只改了与 em2 不同的部分。PIC 模拟被视为无碰撞等离子体的第一性原理模拟：同时求解大量超粒子（电子、离子）的运动方程和麦克斯韦方程组，自洽地再现系统的时空演化。基本算法见 Birdsall and Langdon (1985) 等。em2 是 2D3V 代码（空间二维 x、y，速度三个分量，电子和离子都按相对论处理）。当前版本设置的是双 Harris 电流片的磁重联，模型说明见 RECONNECTION_MODEL.md。下面说明如何编译、运行 em2 以及如何查看输出。代码中的归一化如下：

```
t → ωpe·t        v → v/c        x → x/(c/ωpe)
```

其中 ωpe 用背景密度 n0 定义。密度以 n0 为单位，磁场以 m_e·ωpe/e 为单位（因此 B0 = |wce0|），温度 temp(j) 以 m_j·c² 为单位。

### Step 1. 参数设置

使用 input.f90。下面以当前的双电流片磁重联算例为例。

```
nx = 1024, ny = 512    : x、y 方向的网格数
npcel = 50             : 每个网格每种粒子的超粒子数（背景等离子体；应取 2k²，见模型说明）
nsp = 2                : 粒子种类数（1：电子，2：离子）
nproc = 4              : MPI 进程数（必须与 mpirun 的进程数一致）
nthreads = 4           : OpenMP 线程数（Makefile 没有打开 OpenMP，目前不起作用）
irun = 1               : 运行编号（1：新运行；2、3…：续跑）
nt0 = 20480            : 每次运行的时间步数
isav = 20              : 能量和场时间序列的输出间隔（每 isav 步）；也是 E、B 的滤波间隔
iphs = 1000            : 二维场快照的输出间隔（每 iphs 步）
idata = 1024           : 每个 ex***.hst 类文件中的采样次数（nt0/isav/idata 为整数时每个文件正好写满）
dx = dy = 0.1118       : 网格间距，单位 c/ωpe（要求 dx = dy）
dt = dx*0.7            : 时间步长
wce0 = -1.498147003616 : 远离电流片处的磁场 B0 对应的电子回旋频率，单位 ωpe（|wce0| = B0）
theta = 0              : 背景磁场的方向（本版本不起作用）
rmass = 100            : 离子与电子的质量比
rtemp = 4              : Ti/Te，用来分配电流片的漂移速度（必须等于 temp(2)*rmass/temp(1)）
csw = 2.23606797749979 : 电流片半厚度，单位 c/ωpe（= 20 dx）
csden = 5              : 电流片峰值密度与背景密度之比
trace, parentwave, mno1, b01, v01, v02, freq1 --- 本版本不使用
temp(1) = 4.4888888888D-02        : 电子温度，单位 m_e c²
temp(2) = 1.7955555555D-01/rmass  : 离子温度，单位 m_i c²（Ti = 4 Te）
aist(1), aist(2) = 0   : 温度各向异性，T_x = temp*(1+aist)
va(j), ve(j) = 0       : 本版本不使用（电流片的漂移速度由 wce0、csw、csden、rtemp 决定）
dn(j) = 1              : 相对密度（在 com_den 中）
```

temp(j) 等中的 j 表示粒子种类，j = 1 为电子，j = 2 为离子，当前版本只有这两种。系统大小 Lx = nx·dx = 114.5、Ly = ny·dy = 57.2（c/ωpe），x、y 方向都是周期边界。二维时时间步长要满足 Courant 条件 c·dt < dx/√2 ≈ 0.707 dx，这里 dt = 0.7 dx，已接近上限，不要再加大。数值稳定性要求 dx ≲ λD（德拜长度），λD = (temp(1)/n)^0.5，n 为当地密度：背景区（n = 1）dx ≈ 0.53 λD，电流片中心（n = 1 + csden = 6）dx ≈ 1.3 λD。另外，wce0、csden 和温度要满足压力平衡 B0²/2 = csden·(Te + Ti)（现在 1.122 = 5 × 0.2244），改其中任何一个都要一起调整。

改网格数（nx、ny）时还要注意：
- output.f90 按格式 (1024e13.5) 写场，nx > 1024 时一行会被拆成几行；
- fields.f90 只对 iy = 88..168 这些行求场能量（按 ny = 512 写死）；
- 粒子快照每 4000 步输出一次，写死在 output.f90 中（并且只在 iphs 的整数倍步输出）；
- maxnp（global.f90）会自动计算，但粒子是一行一行加载的，见模型说明。

### Step 2. 编译和运行

用 make 命令编译（选项见 Makefile：编译器 mpif90，选项 -O3）。改过 input.f90 之后，make 之前一定要先 make clean（Makefile 没有写模块依赖）。

```
$ make clean
$ make
```

make 之后生成可执行文件 em2。用 qsub 命令提交作业（PBS），运行脚本是 em2d.sh。

```
$ qsub em2d.sh
```

em2d.sh 申请 1 个节点、4 个进程（#PBS -l nodes=1:ppn=4），用 mpirun 以 4 个进程运行 em2；ppn 必须等于 input.f90 中的 nproc。脚本中的 ulimit -s unlimited（解除栈大小限制）请保留；OMP_NUM_THREADS 目前不起作用。

用 qstat 查看作业状态：S 栏为 R 表示运行中，Q 表示排队中；作业结束后不再显示（或显示为 C）。标准输出写入 em2d.sh.o<作业号>，错误信息写入 em2d.sh.e<作业号>。要取消排队中的作业或强制结束运行中的作业，用 qdel 指定作业号。

```
$ qstat
$ qdel <作业号>
```

不用作业系统时，可以直接运行：

```
$ ulimit -s unlimited
$ mpirun -np 4 ./em2 > em2.log
```

运行时每个进程每 10 步打印一行 'it = ...'，Poisson 修正（cgm）每步打印一行 'n= ... rr0= ... rr1= ...'。如果出现 'MAXITS exceeded in cgm'，说明 Poisson 修正没有收敛。

续跑：把 input.f90 中的 irun 改为 2（3、…），make clean; make，然后在同一目录下提交。cont.f90 读入上次运行结束时写出的 fld.fin、rhoj.fin 和 *_rank###.fin。第一次运行到第 nt0 − isav 步结束，第 irun 次运行到第 irun·nt0 − isav 步结束；输出文件名中的编号（ene002.hst 等）随 irun 增加。cf-p1 文件的编号和帧的位置由步数决定，续跑时接着上一次运行的帧往下写，不会覆盖。

按当前参数，每个进程约需 1.2 GB 内存；一次运行的输出约 50 GB，大部分是粒子快照（每次约 6 GB）。

### Step 3. 输出数据与作图

运行结束后，确认生成了以下文件（在 output.f90 和 init.f90 中指定）：

```
1) ene***.hst                                         能量时间序列
2) xvelc***.hst, xvion***.hst                         粒子轨道
3) ex***.hst, ey***.hst, ez***.hst,
   bx***.hst, by***.hst, bz***.hst                    沿 x = Lx/2 的电磁场 y–t 数据
4) ex_TTTTT.dat … bz_TTTTT.dat                        电磁场的二维快照
5) xp*_TTTTT_rank###.datap                            粒子快照
6) cf-p1_###.bin                                      每 100 步的 Bz（二进制）
7) bx0, by0                                           初始的 Bx、By
8) fld.fin, rhoj.fin, *_rank###.fin                   最终状态（续跑用）
```

（fort.999 是第 1 步 Jz 的调试输出，可以删除。）

与 em1d 不同，除 6) 和 8) 外都是用空格分隔的文本数据，可以直接用 gnuplot 读入作图。em1d 的解析程序（ana.dir 下的 combifld、ptc1d、wkspc、distfunc）不包含在内，一维版本也不能用于 em2。

1) 是能量的时间序列，每 isav 步一行（第一行是第 1 步），共 10 列，从左到右为：

```
时刻 | 电场 3 分量   | 磁场 3 分量   | 电子   离子  | 合计
 t   | Ex² Ey² Ez²   | Bx² By² Bz²   | K(1)   K(2)  | total
```

场能量是 iy = 88..168 行（下电流片 ±2 csw 的范围）内 ½E²、½B² 的平均值（再乘 81/80）；动能 K(1)、K(2) 是各种粒子的总动能除以全部电子的静止质量，即每个电子平均的值（单位 m_e c²）。两者的归一化不同，所以 total 列不是守恒的总能量。第一行（t = dt）的场能量为 0（这时还没有计算）。用 gnuplot 以各能量为纵轴、时间为横轴作图：

```
gnuplot> set logscale y
gnuplot> plot 'ene001.hst' using 1:2 w l, '' using 1:3 w l, '' using 1:4 w l,\
> '' using 1:5 w l, '' using 1:6 w l, '' using 1:7 w l
```

磁重联进行时 Bx²（第 5 列）减小，By²（第 6 列）增大。

2) 是每 isav 步一行的粒子轨道：xvelc 为进程 1 上的前 100 个电子，xvion 为进程 2 上的前 100 个离子（需要 nproc ≥ 3；进程数更少时会写到 fort.70 或 fort.71）。每行依次为 t 和 100 个粒子的 (x, y, vx, vy, vz)，共 501 列（x、y 以 c/ωpe 为单位）：第 k 个粒子的 x、y、vx、vy、vz 在第 5k − 3 … 5k + 1 列。例：第 1 个粒子的 x–y 轨道

```
gnuplot> plot 'xvelc001.hst' using 2:3 w l
```

（2026/10/20 以前的版本用定宽格式 (f10.3,500d12.5) 写这两个文件，负数前面没有空格，不能直接这样读。）这些粒子是从 y ≈ 0 附近出发的背景粒子；要跟踪电流片中的粒子，需修改 output.f90 中的粒子序号范围。粒子数据（2) 和 5)）是在推动粒子之后输出的：标为第 it 步的数据，位置是 (it + 1)·dt 时的，速度是 (it + ½)·dt 时的。

3) 是 x = Lx/2（ix = nx/2）处电磁场沿 y 的分布（iy = 1..ny），每 isav 步一行，每个文件 idata 行。可以直接画 y–t 图：

```
gnuplot> set pm3d map
gnuplot> splot 'bz001.hst' matrix
```

matrix 的序号从 0 开始：横轴为 iy − 1，乘以 dy 即为 Ex、Bx 的 y（其余分量再加 dy/2）；纵轴为记录序号 n，n = 0 是第 1 步，n ≥ 1 是第 n·isav 步（t ≈ n·isav·dt）。x = Lx/2 正好穿过下电流片的 O 点和上电流片的 X 点（见模型说明）。对这些数据做二维傅里叶变换可得到 ω–k 图，例如用 Python：

```python
import numpy as np
isav, dt, dy = 20, 0.07826, 0.1118
d = np.loadtxt('ez001.hst')                                   # [时间, y]
F = np.fft.fftshift(np.abs(np.fft.fft2(d))**2)
w = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(d.shape[0], isav*dt))   # ω
k = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(d.shape[1], dy))        # k_y
```

4) 是每 iphs 步输出一次的二维快照，文件名中的 TTTTT 是时刻 t（单位 1/ωpe）的整数部分（例如第 1000 步为 bz_00078.dat）。每个文件 ny 行 × nx 列（第 iy 行第 ix 列）；6 个分量依次由进程 mod(3, nproc)、mod(4, nproc)、…、mod(8, nproc) 写出。

```
gnuplot> set pm3d map
gnuplot> splot 'bz_00078.dat' matrix
```

横轴为 ix − 1，纵轴为 iy − 1（matrix 的序号从 0 开始）。网格是交错的：记 (j, i) = (ix − 1, iy − 1)，Ex、Bx 位于 (j + ½, i)，Ey、By 位于 (j, i + ½)，Ez、Bz 位于 (j + ½, i + ½)（单位为网格，乘以 dx 即为 c/ωpe）。

5) 是粒子快照，每 4000 步输出一次（写死在 output.f90 中，并且只在 iphs 的整数倍步输出），xp1 为电子、xp2 为离子，每个进程一个文件。每行为 x、y（c/ωpe）和 vx、vy、vz（c），包含该进程的全部粒子（当前参数下每个文件约 1170 万行、0.77 GB）。和 2) 一样，位置比文件名中的时刻晚一步、速度晚半步，所以 xp*_00000 也已经推动过一步。例：离子的 y–vz 相空间散点图（每 100 个粒子画一个）

```
gnuplot> plot 'xp2_00313_rank000.datap' using 2:5 every 100 w d
```

速度分布函数可以对某一区域内粒子的速度做直方图得到（例如 Python 的 numpy.histogram）。

6) 是进程 0 每 100 步写出的 Bz（包括保护格，(nx+3)×(ny+3) 个双精度数，x 方向变化最快），每个文件 30 帧，流式（stream）二进制。用 Python 读入：

```python
import numpy as np
nx, ny = 1024, 512
bz = np.fromfile('cf-p1_001.bin', dtype='<f8').reshape(-1, ny+3, nx+3)   # [帧, iy, ix]
```

第 n 个文件的第 k 帧（k 从 0 开始）对应第 100·(30·(n − 1) + k) 步，续跑时也一样。只有调用 output 的步才会写帧，所以 isav 要能整除 100（现在是 20）。

7) 是初始的 Bx、By，按 iy = 0..ny+1 的顺序，每个 iy 写 nx + 2 个数（ix = 0..nx+1）；格式是 (1024e13.5)，所以 nx = 1024 时每个 iy 占两行（1024 个 + 2 个）。每个进程都写同样的内容。

8) 是续跑用的最终状态（Fortran 无格式顺序文件），续跑前不要删除。

### 附：源文件

```
main.f90      主程序（时间循环）               input.f90    参数
global.f90    共用变量，粒子数 maxnp            init.f90     初始磁场，粒子数组
creator.f90   粒子加载（位置、速度、漂移）      smear.f90    打乱速度的顺序
setv.f90      把初速度推回 −dt/2               move.f90     推动速度（Boris 法）
curnt.f90     推动位置，分配电荷和电流          fields.f90   场的推进、滤波、Poisson 修正
smooth.f90    滤波（smooth2）                  output.f90   输出
cont.f90      续跑时读入                        dran2.f90    随机数
setrho.f90    （未使用）                        em2d.sh      PBS 运行脚本
temp1.data, temp2.data   以前的后处理输出，与运行无关
```

---

## English

### em2 first step memo (2026/10/20)

This is a brief manual explaining how to use the two-dimensional periodic full particle-in-cell (PIC) code em2. It is adapted from S. Matsukiyo's "em1d first step memo"; only the parts that differ for em2 have been changed. PIC simulation is regarded as a first principle simulation of a collisionless plasma. The equations of motion of a number of super-particles (electrons and ions) and Maxwell's equations are simultaneously solved. The basic algorithm is explained by many authors (e.g., Birdsall and Langdon, 1985). em2 is a 2D3V code (two spatial dimensions x, y and three velocity components; relativistic electrons and ions). The present version sets up magnetic reconnection in a double Harris current sheet; the model is described in RECONNECTION_MODEL.md. In the following we explain how to compile em2, execute it, and look at the output. The normalization used in the code is as follows.

```
t → ωpe·t        v → v/c        x → x/(c/ωpe)
```

Here ωpe is defined with the background density n0. Density is normalized by n0, the magnetic field by m_e·ωpe/e (so that B0 = |wce0|), and the temperature temp(j) by m_j·c².

### Step 1. Parameter settings

The 'input.f90' is used. The present double current sheet reconnection run is taken as the example.

```
nx = 1024, ny = 512    : number of grid cells in x and y
npcel = 50             : number of super-particles per cell for each species (background; use 2k², see the model document)
nsp = 2                : number of species (1: electrons, 2: ions)
nproc = 4              : number of MPI processes (must equal the number given to mpirun)
nthreads = 4           : number of OpenMP threads (no effect: OpenMP is not enabled in the Makefile)
irun = 1               : index of run (1: new run; 2, 3, ...: continuation runs)
nt0 = 20480            : iteration of time steps in each run
isav = 20              : rate of output of energies and field histories (every 'isav' time steps); also the interval of E, B filtering
iphs = 1000            : rate of output of 2D field snapshots (every 'iphs' time steps)
idata = 1024           : number of samplings in each ex***.hst-type file (every file is full when nt0/isav/idata is an integer)
dx = dy = 0.1118       : grid spacing in unit of c/ωpe (dx = dy is assumed)
dt = dx*0.7            : time step
wce0 = -1.498147003616 : electron gyro frequency of the asymptotic field B0 in unit of ωpe (|wce0| = B0)
theta = 0              : direction of a uniform background field (no effect in this version)
rmass = 100            : ion-to-electron mass ratio
rtemp = 4              : Ti/Te used to share the current-sheet drift (must equal temp(2)*rmass/temp(1))
csw = 2.23606797749979 : half-thickness of the current sheets in unit of c/ωpe (= 20 dx)
csden = 5              : peak density of each current sheet relative to the background
trace, parentwave, mno1, b01, v01, v02, freq1 --- not used in this version
temp(1) = 4.4888888888D-02        : temperature of electrons in unit of me c²
temp(2) = 1.7955555555D-01/rmass  : temperature of ions in unit of mi c² (Ti = 4 Te)
aist(1), aist(2) = 0   : temperature anisotropy, T_x = temp*(1+aist)
va(j), ve(j) = 0       : not used (the current-sheet drift is set from wce0, csw, csden and rtemp)
dn(j) = 1              : relative density (in com_den)
```

The 'j' in temp(j), etc., represents particle species: electrons (j=1) and ions (j=2); this version has only these two. The system size is Lx = nx·dx = 114.5 and Ly = ny·dy = 57.2 (c/ωpe), periodic in both x and y. In 2D the time step must satisfy the Courant condition c·dt < dx/√2 ≈ 0.707 dx; here dt = 0.7 dx, close to the limit, so do not increase it. For numerical stability, dx ≲ λD (Debye length) is necessary, with λD = (temp(1)/n)^0.5 for the local density n: dx ≈ 0.53 λD in the background (n = 1) and dx ≈ 1.3 λD at the sheet centre (n = 1 + csden = 6). wce0, csden and the temperatures must also satisfy the pressure balance B0²/2 = csden·(Te + Ti) (now 1.122 = 5 × 0.2244); change them together.

When changing the grid (nx, ny), note also:
- output.f90 writes the fields with the format (1024e13.5), so for nx > 1024 one row is split over several lines;
- fields.f90 sums the field energy over the rows iy = 88..168 only (fixed for ny = 512);
- particle snapshots are written every 4000 steps, fixed in output.f90 (and only at multiples of iphs);
- maxnp (global.f90) is computed automatically, but the particles are loaded row by row; see the model document.

### Step 2. Compilation and execution

Use the make command (see the 'Makefile' for options: mpif90, -O3). After editing input.f90, always run make clean before make (the Makefile has no module dependencies).

```
$ make clean
$ make
```

After the make command, an executable file 'em2' is produced. Use the qsub command for job submission (PBS). The shell script for execution is 'em2d.sh'.

```
$ qsub em2d.sh
```

em2d.sh requests one node with 4 processes (#PBS -l nodes=1:ppn=4) and runs em2 with mpirun on 4 processes; ppn must equal nproc in input.f90. Keep the line ulimit -s unlimited (removes the stack size limit); OMP_NUM_THREADS has no effect at present.

To see the job status, the qstat command is used: R in the S column means running and Q queued; a finished job disappears (or is shown as C). Standard output goes to em2d.sh.o<job-ID> and error messages to em2d.sh.e<job-ID>. If one forces the job to quit, the qdel command is used.

```
$ qstat
$ qdel <job-ID>
```

Without a job system, run it directly:

```
$ ulimit -s unlimited
$ mpirun -np 4 ./em2 > em2.log
```

During the run each process prints 'it = ...' every 10 steps, and the Poisson correction (cgm) prints 'n= ... rr0= ... rr1= ...' every step. 'MAXITS exceeded in cgm' means that the Poisson correction did not converge.

Continuation run: set irun = 2 (3, ...) in input.f90, make clean; make, and submit in the same directory. cont.f90 reads fld.fin, rhoj.fin and *_rank###.fin written at the end of the previous run. The first run stops at step nt0 − isav, and run number irun at step irun·nt0 − isav; the number in the output file names (ene002.hst, etc.) increases with irun. The cf-p1 file number and the position of a frame follow from the step number, so a continuation run writes on after the frames of the previous run without overwriting them.

With the present parameters each process needs about 1.2 GB of memory, and one run writes about 50 GB, mostly particle snapshots (about 6 GB each time).

### Step 3. Data analysis and visualization

The following files are output (specified in 'output.f90' and 'init.f90').

```
1) ene***.hst                                         Energy time history
2) xvelc***.hst, xvion***.hst                         Particle trajectories
3) ex***.hst, ey***.hst, ez***.hst,
   bx***.hst, by***.hst, bz***.hst                    y-t data of EM fields along x = Lx/2
4) ex_TTTTT.dat ... bz_TTTTT.dat                      2D snapshots of EM fields
5) xp*_TTTTT_rank###.datap                            Snapshot data for particles
6) cf-p1_###.bin                                      Bz every 100 steps (binary)
7) bx0, by0                                           Initial Bx and By
8) fld.fin, rhoj.fin, *_rank###.fin                   Data of final state quantities (for continuation)
```

(fort.999 is a debug output of Jz at step 1 and can be deleted.)

Unlike em1d, all files except 6) and 8) are whitespace-separated text and can be plotted directly with gnuplot. The em1d analysis programs (combifld, ptc1d, wkspc and distfunc in ana.dir) are not included, and the 1D versions do not apply to em2.

1) contains the time series of energies, one line every 'isav' steps (the first line is step 1), in 10 columns:

```
Time | E-field        | B-field        | Elec.  Ion   | Total
 t   | Ex² Ey² Ez²    | Bx² By² Bz²    | K(1)   K(2)  | total
```

The field energies are averages of ½E² and ½B² over the rows iy = 88..168 (the lower current sheet ± 2 csw), multiplied by 81/80. The kinetic energies K(1) and K(2) are the total kinetic energy of each species divided by the total rest mass of the electrons, i.e., per electron in unit of me c². Because the two are normalized differently, the 'total' column is not the conserved total energy. The field energies in the first line (t = dt) are zero (not yet computed). An example of a plot by using gnuplot:

```
gnuplot> set logscale y
gnuplot> plot 'ene001.hst' using 1:2 w l, '' using 1:3 w l, '' using 1:4 w l,\
> '' using 1:5 w l, '' using 1:6 w l, '' using 1:7 w l
```

As reconnection proceeds, Bx² (column 5) decreases and By² (column 6) increases.

2) contains trajectories every 'isav' steps: xvelc the first 100 electrons on process 1, xvion the first 100 ions on process 2 (this needs nproc ≥ 3; with fewer processes the data go to fort.70 or fort.71). Each line has t and (x, y, vx, vy, vz) of the 100 particles, 501 columns in total (x, y in c/ωpe): x, y, vx, vy, vz of particle k are in columns 5k − 3 ... 5k + 1. Example: the x-y orbit of the first particle,

```
gnuplot> plot 'xvelc001.hst' using 2:3 w l
```

(Versions before 2026/10/20 wrote these two files in the fixed-width format (f10.3,500d12.5), with no blank before negative numbers, which cannot be read this way.) These are background particles starting near y ≈ 0; to follow particles in a current sheet, change the slot range in output.f90. Particle data (2) and 5)) are written after the particle push: in the data labelled step it, positions are at (it + 1)·dt and velocities at (it + ½)·dt.

3) contains the profiles of the field components along y (iy = 1..ny) at x = Lx/2 (ix = nx/2), one line every 'isav' steps and 'idata' lines per file. A y-t diagram is plotted directly:

```
gnuplot> set pm3d map
gnuplot> splot 'bz001.hst' matrix
```

The matrix indices start at 0: the horizontal axis is iy − 1, which times dy is y for Ex and Bx (add dy/2 for the other components); the vertical axis is the record number n, where n = 0 is step 1 and n ≥ 1 is step n·isav (t ≈ n·isav·dt). The line x = Lx/2 passes through the O point of the lower sheet and the X point of the upper sheet (see the model document). A 2D Fourier transform gives the ω-k diagram, e.g., in Python:

```python
import numpy as np
isav, dt, dy = 20, 0.07826, 0.1118
d = np.loadtxt('ez001.hst')                                   # [time, y]
F = np.fft.fftshift(np.abs(np.fft.fft2(d))**2)
w = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(d.shape[0], isav*dt))   # omega
k = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(d.shape[1], dy))        # k_y
```

4) are 2D snapshots written every 'iphs' steps. TTTTT in the file name is the integer part of the time t in 1/ωpe (e.g., bz_00078.dat at step 1000). A file has ny lines × nx columns (line iy, column ix); the six components are written by processes mod(3, nproc), mod(4, nproc), ..., mod(8, nproc) in turn.

```
gnuplot> set pm3d map
gnuplot> splot 'bz_00078.dat' matrix
```

The horizontal axis is ix − 1 and the vertical axis iy − 1 (matrix indices start at 0). The grid is staggered: with (j, i) = (ix − 1, iy − 1), Ex and Bx are at (j + ½, i), Ey and By at (j, i + ½), and Ez and Bz at (j + ½, i + ½), in units of cells (multiply by dx for c/ωpe).

5) is the snapshot data of particles, written every 4000 steps (fixed in output.f90, and only at multiples of 'iphs'): xp1 for electrons and xp2 for ions, one file per process. Each line has x, y (c/ωpe) and vx, vy, vz (c), for all particles of the process (with the present parameters about 11.7 million lines, 0.77 GB per file). As in 2), positions are one step and velocities half a step later than the time in the file name, so even xp*_00000 has been pushed once. Example: a scatter plot of ions in y-vz phase space (every 100th particle),

```
gnuplot> plot 'xp2_00313_rank000.datap' using 2:5 every 100 w d
```

Velocity distribution functions are obtained by histogramming the velocities of the particles in a region (e.g., numpy.histogram in Python).

6) is Bz written by process 0 every 100 steps, including the guard cells ((nx+3)×(ny+3) double precision numbers, x running fastest), 30 frames per file, stream binary. In Python:

```python
import numpy as np
nx, ny = 1024, 512
bz = np.fromfile('cf-p1_001.bin', dtype='<f8').reshape(-1, ny+3, nx+3)   # [frame, iy, ix]
```

Frame k (from 0) of file number n corresponds to step 100·(30·(n − 1) + k), also in continuation runs. Frames are written only at steps where output is called, so isav must divide 100 (now 20).

7) are the initial Bx and By: for iy = 0..ny+1 in turn, nx + 2 numbers (ix = 0..nx+1). The format is (1024e13.5), so with nx = 1024 each iy takes two lines (1024 + 2 numbers). Every process writes the same content.

8) is the final state for a continuation run (Fortran unformatted sequential files). Do not delete them before the continuation run.

### Appendix: source files

```
main.f90      main program (time loop)          input.f90    parameters
global.f90    shared variables, maxnp            init.f90     initial magnetic field, particle arrays
creator.f90   particle loading (x, v, drift)     smear.f90    shuffles the velocities
setv.f90      pushes the initial v back by dt/2  move.f90     velocity push (Boris)
curnt.f90     position push, rho and J           fields.f90   field solver, filter, Poisson correction
smooth.f90    filter (smooth2)                   output.f90   output
cont.f90      reads the restart files            dran2.f90    random numbers
setrho.f90    (not used)                         em2d.sh      PBS job script
temp1.data, temp2.data   old post-processing output, not used by the run
```

---

## 日本語

### em2 ファーストステップメモ（2026/10/20）

2 次元周期境界フル粒子シミュレーションコード em2 の使い方を概説する。S. Matsukiyo「em1d ファーストステップメモ」をもとに、em2 で異なる部分だけを書き換えたものである。フル粒子シミュレーションはプラズマの第 1 原理シミュレーションといわれ、しばしば PIC (Particle-In-Cell) シミュレーションと呼ばれる。多数のプラズマ超粒子（電子・イオン）の運動方程式と、電磁場の時間・空間発展を記述するマックスウェル方程式を連立させて解くことで、系の時間・空間発展を自己無撞着に再現する。一般的なフル粒子シミュレーション法の解説は、例えば「2003 天体とスペースプラズマのシミュレーションサマースクール」テキストを参照のこと。em2 は空間 2 次元（x, y）・速度 3 成分の相対論的コードで、現在のバージョンは 2 枚の Harris 電流シートによる磁気リコネクションを設定している（モデルの説明は RECONNECTION_MODEL.md）。ここでは、em2 をコンパイル・実行し、出力を確認するための例を示す。ただし、規格化は以下のとおりとする。

```
t → ωpe·t        v → v/c        x → x/(c/ωpe)
```

ωpe は背景密度 n0 で定義する。密度は n0、磁場は m_e·ωpe/e で規格化し（B0 = |wce0|）、温度 temp(j) は m_j·c² を単位とする。

### Step 1. パラメータ指定

最初に、ソースコードの input.f90 を開いてパラメータ設定を行う。ここでは、現在の 2 枚電流シートの磁気リコネクションを例にとる。

```
nx = 1024, ny = 512    : number of grid cells in x and y
npcel = 50             : number of super-particles per cell for each species (background; use 2k², see the model document)
nsp = 2                : number of species (1: electrons, 2: ions)
nproc = 4              : number of MPI processes (must equal the number given to mpirun)
nthreads = 4           : number of OpenMP threads (no effect: OpenMP is not enabled in the Makefile)
irun = 1               : index of run (1: new run; 2, 3, ...: continuation runs)
nt0 = 20480            : iteration of time steps in each run
isav = 20              : rate of output of energies and field histories (every 'isav' time steps); also the interval of E, B filtering
iphs = 1000            : rate of output of 2D field snapshots (every 'iphs' time steps)
idata = 1024           : number of samplings in each ex***.hst-type file (every file is full when nt0/isav/idata is an integer)
dx = dy = 0.1118       : grid spacing in unit of c/ωpe (dx = dy is assumed)
dt = dx*0.7            : time step
wce0 = -1.498147003616 : electron gyro frequency of the asymptotic field B0 in unit of ωpe (|wce0| = B0)
theta = 0              : direction of a uniform background field (no effect in this version)
rmass = 100            : ion-to-electron mass ratio
rtemp = 4              : Ti/Te used to share the current-sheet drift (must equal temp(2)*rmass/temp(1))
csw = 2.23606797749979 : half-thickness of the current sheets in unit of c/ωpe (= 20 dx)
csden = 5              : peak density of each current sheet relative to the background
trace, parentwave, mno1, b01, v01, v02, freq1 --- not used in this version
temp(1) = 4.4888888888D-02        : temperature of electrons in unit of me c²
temp(2) = 1.7955555555D-01/rmass  : temperature of ions in unit of mi c² (Ti = 4 Te)
aist(1), aist(2) = 0   : temperature anisotropy, T_x = temp*(1+aist)
va(j), ve(j) = 0       : not used (the current-sheet drift is set from wce0, csw, csden and rtemp)
dn(j) = 1              : relative density (in com_den)
```

temp(j) などの j は粒子種を表し、j = 1 が電子、j = 2 がイオンである（現在のバージョンはこの 2 成分のみ）。系の大きさは Lx = nx·dx = 114.5、Ly = ny·dy = 57.2（c/ωpe）で、x、y 方向とも周期境界である。2 次元では時間ステップ幅が Courant 条件 c·dt < dx/√2 ≈ 0.707 dx を満たす必要があり、ここでは dt = 0.7 dx としている（上限に近いので、これ以上大きくしないこと）。数値的安定性のためには dx ≲ λD（デバイ長）であることが必要とされており、λD = (temp(1)/n)^0.5（n は局所密度）である。背景（n = 1）では dx ≈ 0.53 λD、電流シート中心（n = 1 + csden = 6）では dx ≈ 1.3 λD となる。また、wce0、csden、温度は圧力平衡 B0²/2 = csden·(Te + Ti) を満たす必要がある（現在 1.122 = 5 × 0.2244）。どれかを変えるときは合わせて調整する。

格子数（nx, ny）を変えるときは、以下にも注意する。
- output.f90 は場を書式 (1024e13.5) で書くため、nx > 1024 では 1 行が複数行に分かれる。
- fields.f90 の場のエネルギーは iy = 88..168 の行だけで計算している（ny = 512 用に固定）。
- 粒子スナップショットは 4000 ステップごと（output.f90 で固定、かつ iphs の倍数のステップのみ）。
- maxnp（global.f90）は自動で計算されるが、粒子は行単位で配置する（モデルの説明を参照）。

### Step 2. コンパイルおよび実行

コンパイルは make コマンドで行う（input.f90 を書き換えた後は、make の前に make clean を実行。Makefile にモジュールの依存関係が書かれていないため）。オプション等は Makefile を参照（mpif90、-O3）。

```
$ make clean
$ make
```

make を行うと実行可能ファイル em2 が作成される。実行には PBS の qsub コマンドを使用する。実行用シェルスクリプトは em2d.sh である。

```
$ qsub em2d.sh
```

em2d.sh は 1 ノード 4 プロセス（#PBS -l nodes=1:ppn=4）を要求し、mpirun で em2 を 4 プロセス実行する。ppn は input.f90 の nproc と一致させる。ulimit -s unlimited（スタックサイズ制限の解除）は残しておくこと。OMP_NUM_THREADS は現在効果がない。

ジョブの状態を見るには qstat コマンドを実行する。S 欄が R なら実行中、Q なら待機中である。終了したジョブは表示されなくなる（または C と表示される）。標準出力は em2d.sh.o<job-ID> に、エラーメッセージは em2d.sh.e<job-ID> に出力される。待機中のジョブをキャンセル、あるいは実行中のジョブを強制終了する場合は、qdel コマンドで job-ID を指定する。

```
$ qstat
$ qdel <job-ID>
```

ジョブ管理システムを使わない場合は、次のように直接実行する。

```
$ ulimit -s unlimited
$ mpirun -np 4 ./em2 > em2.log
```

実行中は、各プロセスが 10 ステップごとに 'it = ...' を、Poisson 補正（cgm）が毎ステップ 'n= ... rr0= ... rr1= ...' を出力する。'MAXITS exceeded in cgm' が出た場合は、Poisson 補正が収束していない。

継続 run：input.f90 の irun を 2（3, ...）に変えて make clean; make し、同じディレクトリで実行する。cont.f90 が前回の run の最後に書き出した fld.fin、rhoj.fin、*_rank###.fin を読み込む。最初の run は nt0 − isav ステップ目で終わり、irun 回目の run は irun·nt0 − isav ステップ目まで進む。出力ファイル名の番号（ene002.hst など）は irun とともに増える。cf-p1 ファイルの番号とフレームの位置はステップ数から決まるので、継続 run では前回の run のフレームの後に続けて書き、上書きはしない。

現在のパラメータでは、1 プロセスあたり約 1.2 GB のメモリを使い、1 回の run の出力は約 50 GB になる（大部分は粒子スナップショットで、1 回あたり約 6 GB）。

### Step 3. 作図およびデータ解析

実行が終了したら、以下のファイルができていることを確認する。なお、出力ファイルは output.f90 と init.f90 で指定している。

```
① ene***.hst                                         エネルギー時系列データ
② xvelc***.hst, xvion***.hst                         粒子軌道データ
③ ex***.hst, ey***.hst, ez***.hst,
   bx***.hst, by***.hst, bz***.hst                    x = Lx/2 に沿った電磁場の y–t データ
④ ex_TTTTT.dat … bz_TTTTT.dat                        電磁場の 2 次元スナップショット
⑤ xp*_TTTTT_rank###.datap                            粒子のスナップショットデータ
⑥ cf-p1_###.bin                                      100 ステップごとの Bz（バイナリ）
⑦ bx0, by0                                           初期の Bx, By
⑧ fld.fin, rhoj.fin, *_rank###.fin                   最終状態を保存したデータ（継続 run 用）
```

（fort.999 は 1 ステップ目の Jz のデバッグ出力で、削除してよい。）

em1d と異なり、⑥と⑧以外はすべて空白区切りのテキストデータで、そのまま gnuplot で読み込んでプロットできる。em1d の解析用プログラム（ana.dir の combifld、ptc1d、wkspc、distfunc）は含まれていない（1 次元用なので em2 には使えない）。

①には、各種エネルギーの時系列データが isav ステップごとに格納されている（1 行目は 1 ステップ目）。10 列のデータから成り、左から次のデータが並ぶ。

```
時刻 | 電場 3 成分    | 磁場 3 成分    | 電子   イオン | 合計
 t   | Ex² Ey² Ez²    | Bx² By² Bz²    | K(1)   K(2)   | total
```

場のエネルギーは iy = 88..168 の行（下の電流シート ±2 csw の範囲）での ½E²、½B² の平均（×81/80）、運動エネルギー K(1)、K(2) は各粒子種の全運動エネルギーを全電子の静止質量で割ったもの（電子 1 個あたり、単位 me c²）である。規格化が異なるので、total 列は保存する全エネルギーではない。1 行目（t = dt）の場のエネルギーは 0 である（まだ計算されていない）。gnuplot を用いて、縦軸を各エネルギー、横軸を時刻とする 2 次元プロットを作成するなら、以下のようにすればよい。

```
gnuplot> set logscale y
gnuplot> plot 'ene001.hst' using 1:2 w l, '' using 1:3 w l, '' using 1:4 w l,\
> '' using 1:5 w l, '' using 1:6 w l, '' using 1:7 w l
```

リコネクションが進むと、Bx²（5 列目）が減少し、By²（6 列目）が増加する。

②は isav ステップごとの粒子軌道データで、xvelc はプロセス 1 の最初の 100 個の電子、xvion はプロセス 2 の最初の 100 個のイオンである（nproc ≥ 3 が必要。プロセス数が少ないと fort.70 または fort.71 に書かれる）。各行は t と 100 個分の (x, y, vx, vy, vz) で、計 501 列（x, y は c/ωpe 単位）。k 番目の粒子の x, y, vx, vy, vz は 5k − 3 … 5k + 1 列目にある。例：1 個目の粒子の x–y 軌道

```
gnuplot> plot 'xvelc001.hst' using 2:3 w l
```

（2026/10/20 より前のバージョンはこの 2 つのファイルを固定幅の書式 (f10.3,500d12.5) で書いており、負の数の前に空白がないため、このままでは読めない。）これらは y ≈ 0 付近から出発する背景粒子であり、電流シート中の粒子を追う場合は output.f90 の粒子番号の範囲を変える。粒子のデータ（②と⑤）は粒子を動かした後に出力される：ステップ it のデータでは、位置は (it + 1)·dt、速度は (it + ½)·dt のものである。

③は x = Lx/2（ix = nx/2）での電磁場の y 方向分布（iy = 1..ny）を isav ステップごとに 1 行ずつ、1 ファイルに idata 行書いたものである。そのまま y–t 図が得られる。

```
gnuplot> set pm3d map
gnuplot> splot 'bz001.hst' matrix
```

matrix の番号は 0 から始まる：横軸は iy − 1 で、dy をかけると Ex, Bx の y になる（他の成分はさらに dy/2 を足す）。縦軸はレコード番号 n で、n = 0 が 1 ステップ目、n ≥ 1 が n·isav ステップ目（t ≈ n·isav·dt）である。x = Lx/2 は下の電流シートの O 点と上の電流シートの X 点を通る（モデルの説明を参照）。2 次元フーリエ変換を施すと ω–k 図が得られる。例えば Python では次のようにする。

```python
import numpy as np
isav, dt, dy = 20, 0.07826, 0.1118
d = np.loadtxt('ez001.hst')                                   # [時刻, y]
F = np.fft.fftshift(np.abs(np.fft.fft2(d))**2)
w = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(d.shape[0], isav*dt))   # ω
k = 2*np.pi*np.fft.fftshift(np.fft.fftfreq(d.shape[1], dy))        # k_y
```

④は iphs ステップごとの 2 次元スナップショットで、ファイル名の TTTTT は時刻 t（1/ωpe 単位）の整数部分である（例：1000 ステップ目は bz_00078.dat）。中身は ny 行 × nx 列（iy 行目の ix 列目）で、6 つの成分を順にプロセス mod(3, nproc)、mod(4, nproc)、…、mod(8, nproc) が書く。

```
gnuplot> set pm3d map
gnuplot> splot 'bz_00078.dat' matrix
```

横軸は ix − 1、縦軸は iy − 1 である（matrix の番号は 0 から始まる）。格子はスタッガードで、(j, i) = (ix − 1, iy − 1) とすると、Ex, Bx は (j + ½, i)、Ey, By は (j, i + ½)、Ez, Bz は (j + ½, i + ½) にある（格子単位。dx をかけると c/ωpe 単位）。

⑤は粒子のスナップショットデータで、4000 ステップごと（output.f90 で固定、かつ iphs の倍数のステップのみ）に出力している。xp1 が電子、xp2 がイオンで、プロセスごとに 1 ファイル。各行は x, y（c/ωpe）と vx, vy, vz（c）で、そのプロセスの全粒子が書かれる（現在のパラメータで 1 ファイル約 1170 万行、0.77 GB）。②と同様に、位置はファイル名の時刻より 1 ステップ、速度は半ステップ後のものなので、xp*_00000 もすでに 1 回動かした後である。例：イオンの y–vz 散布図（100 個に 1 個）

```
gnuplot> plot 'xp2_00313_rank000.datap' using 2:5 every 100 w d
```

速度分布関数は、ある領域内の粒子の速度のヒストグラムをとれば得られる（例えば Python の numpy.histogram）。

⑥はプロセス 0 が 100 ステップごとに書き出す Bz（ガードセルを含む (nx+3)×(ny+3) 個の倍精度実数、x 方向が最も速く変わる）で、1 ファイル 30 フレームのストリーム形式のバイナリである。Python での読み込み例：

```python
import numpy as np
nx, ny = 1024, 512
bz = np.fromfile('cf-p1_001.bin', dtype='<f8').reshape(-1, ny+3, nx+3)   # [フレーム, iy, ix]
```

n 番目のファイルの k フレーム目（k は 0 から）は 100·(30·(n − 1) + k) ステップ目に対応する（継続 run でも同じ）。フレームは output が呼ばれるステップでのみ書かれるので、isav は 100 を割り切る値にする（現在 20）。

⑦は初期の Bx, By で、iy = 0..ny+1 の順に、各 iy について nx + 2 個（ix = 0..nx+1）の値を書く。書式が (1024e13.5) なので、nx = 1024 では各 iy が 2 行（1024 個 + 2 個）になる。全プロセスが同じ内容を書く。

⑧は継続 run 用の最終状態（Fortran の書式なし順次ファイル）である。継続 run の前に消さないこと。

### 付録：ソースファイル

```
main.f90      メインプログラム（時間ループ）     input.f90    パラメータ
global.f90    共有変数、粒子数 maxnp              init.f90     初期磁場、粒子配列
creator.f90   粒子の配置（位置・速度・ドリフト）  smear.f90    速度の並べ替え
setv.f90      初速度を −dt/2 に戻す               move.f90     速度の更新（Boris 法）
curnt.f90     位置の更新、電荷・電流の割り当て     fields.f90   場の更新、フィルター、Poisson 補正
smooth.f90    フィルター（smooth2）               output.f90   出力
cont.f90      継続 run の読み込み                  dran2.f90    乱数
setrho.f90    （未使用）                           em2d.sh      PBS 実行スクリプト
temp1.data, temp2.data   以前の後処理の出力（実行には不要）
```

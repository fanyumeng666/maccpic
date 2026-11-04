#!/usr/bin/env python3
"""draw_flow.py [out_dir]: functional flow chart of ccp1d with the completion
status of every block (green = done and tested, orange = to do, grey = data).
Writes program_flow.png and program_flow.svg (default: doc/fig/)."""
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon
import matplotlib.font_manager as fm

for name in ['WenQuanYi Micro Hei', 'WenQuanYi Zen Hei', 'Noto Sans CJK SC']:
    if any(f.name == name for f in fm.fontManager.ttflist):
        plt.rcParams['font.sans-serif'] = [name, 'DejaVu Sans']
        break
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.unicode_minus'] = False

out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'doc', 'fig')
os.makedirs(out, exist_ok=True)

STYLE = {
    'done': dict(fc='#dcf0dc', ec='#2e7d32', ls='-', lw=1.4),
    'todo': dict(fc='#fff1d6', ec='#d97706', ls='--', lw=1.4),
    'data': dict(fc='#ececec', ec='#666666', ls='-', lw=1.2),
    'frame': dict(fc='none', ec='#5b6b7c', ls='-', lw=1.0),
}
fig, ax = plt.subplots(figsize=(22, 14.5))
ax.set_xlim(0, 100); ax.set_ylim(0, 66); ax.axis('off')


GLYPHS = {'\u2212': '-', '\u207a': '+', '\u207b': '-', '\u00bd': '1/2', '\u27e8': '<', '\u27e9': '>', '\u221a': 'sqrt',
          '\u222b': 'integral ', '\u22a5': 'perp', '\u230a': 'floor(', '\u230b': ')', '\u2099': 'n', '\u208a': '+',
          '\u2081': '1', '\u2082': '2', '\u208b': '-', 'n\u0302': 'n_hat', '\u00b2': '^2', '\u00b3': '^3', '\u2080': '0'}


def clean(t):
    for k, v in GLYPHS.items():
        t = t.replace(k, v)
    return t


def box(x, y, w, h, title, lines, status, tag=None, tfs=9.2, lfs=8.0):
    st = STYLE[status]
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.25,rounding_size=0.6', fc=st['fc'], ec=st['ec'],
                                ls=st['ls'], lw=st['lw'], zorder=2))
    title = clean(title)
    fam = 'monospace' if title.isascii() else 'sans-serif'
    ax.text(x + 0.6, y + h - 0.55, title, fontsize=tfs, fontweight='bold', family=fam, va='top', ha='left', zorder=3)
    if lines:
        ax.text(x + 0.6, y + h - 0.55 - 1.45, clean(lines), fontsize=lfs, va='top', ha='left', zorder=3, linespacing=1.25)
    if tag:
        ax.text(x + w - 0.5, y + h - 0.5, tag, fontsize=7.5, va='top', ha='right', color=st['ec'], fontweight='bold', zorder=3)
    return (x, y, w, h)


def frame(x, y, w, h, title, color='#5b6b7c'):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.3,rounding_size=0.8', fc='none', ec=color, lw=1.2, zorder=1))
    ax.text(x + 0.7, y + h - 0.4, clean(title), fontsize=9.5, fontweight='bold', color=color, va='top', ha='left', zorder=3)


def diamond(cx, cy, w, h, text):
    ax.add_patch(Polygon([(cx - w / 2, cy), (cx, cy + h / 2), (cx + w / 2, cy), (cx, cy - h / 2)], closed=True,
                         fc='#e8eef7', ec='#2c5aa0', lw=1.3, zorder=2))
    ax.text(cx, cy, text, fontsize=8.5, ha='center', va='center', zorder=3)


def arrow(p0, p1, label=None, style='-|>', color='#333333', rad=0.0, lpos=0.5, loff=(0, 0), lw=1.2):
    a = FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=12, color=color, lw=lw,
                        connectionstyle='arc3,rad=%g' % rad, zorder=4, shrinkA=0, shrinkB=0)
    ax.add_patch(a)
    if label:
        lx = p0[0] + lpos * (p1[0] - p0[0]) + loff[0]
        ly = p0[1] + lpos * (p1[1] - p0[1]) + loff[1]
        ax.text(lx, ly, clean(label), fontsize=7.5, color=color, ha='center', va='center', zorder=5,
                bbox=dict(fc='white', ec='none', pad=0.5))


def bottom(b):
    return (b[0] + b[2] / 2, b[1])


def top(b):
    return (b[0] + b[2] / 2, b[1] + b[3])


def right(b, f=0.5):
    return (b[0] + b[2], b[1] + f * b[3])


def left(b, f=0.5):
    return (b[0], b[1] + f * b[3])


# ---------------- title ----------------
ax.text(1, 65.2, 'ccp1d 功能流程图：一维有界静电 PIC-MCC（目标：Turner 等 2013 氦基准 Case 1）', fontsize=14, fontweight='bold', va='top')
ax.text(1, 63.3, '每个方框：例程名（等宽字体）+ 它做的事。绿色实线 = 已完成并通过测试（M0–M7 全部）；橙色虚线 = 待完成；灰色 = 数据。', fontsize=9, va='top', color='#333')

# ---------------- column 1: setup ----------------
frame(0.5, 1.2, 21.5, 60.4, '① 启动（main.f90 开头）')
A1 = box(1.3, 55.5, 20, 4.6, 'MPI_INIT / COMM_SIZE / COMM_RANK', 'nprocs 个进程，各管自己的粒子；\n每个进程都持有整个网格（nx+1 个节点）', 'done')
A2 = box(1.3, 49.6, 20, 4.9, 'input: read_input', 'namelist ccp1d.in（几何、电压、气体、粒子数、步数…）\n导出 dx=L/nx, dt=1/(f·spp), w=n0L/N, maxnp', 'done')
A3 = box(1.3, 44.2, 20, 4.4, 'random: rng_init(seed+rank)', 'ran2（L\'Ecuyer 组合 LCG + 洗牌）；Box–Muller 高斯；\n状态可存取 → 续算逐位复现', 'done')
A4 = box(1.3, 38.6, 20, 4.6, 'allocate + poisson/push/diag _init', '粒子数组 (maxnp,2)、网格数组 (0:nx)；\nThomas 前扫系数预存；诊断累加器清零', 'done')
A5 = box(1.3, 32.8, 20, 4.8, 'xsec: xsec_load + compute_numax', '读 6 张截面表 → 均匀能量网格（原点对齐表首点，线性插值精确）；\nν_max = max n_g σ_tot(ε) v(ε)，电子 8.65e7 s^-1', 'done')
diamond(11.3, 28.7, 7.2, 3.2, 'restart = .true.?')
A7a = box(1.3, 19.6, 9.6, 6.2, 'init: init_fresh', 'x 均匀随机，v 三分量\nMaxwell；沉积 ρ(0)，解\nφ(0)；v 回推 −dt/2\n（蛙跳起步）', 'done', lfs=7.6)
A7b = box(11.7, 19.6, 9.6, 6.2, 'diag: restart_read', 'restart_rankNNN.bin：\n粒子、随机数状态、\n累加器；restart_grid.bin：\nφ、E、it、窗口状态', 'done', lfs=7.6)
A8 = box(1.3, 13.8, 20, 4.4, 'write_info_header / history 第 0 行', 'info.txt：全部参数 + 导出量（ωpe·dt, λ_D/dx, …）\nhistory.dat：初始 N、KE、W_field', 'done')
# legend
ax.text(1.3, 11.2, '图例', fontsize=9, fontweight='bold', va='top')
for i, (lab, st) in enumerate([('已完成并通过测试', 'done'), ('待完成', 'todo'), ('数据 / 外部输入', 'data')]):
    s = STYLE[st]
    ax.add_patch(FancyBboxPatch((1.5, 8.6 - 2.2 * i), 3.2, 1.5, boxstyle='round,pad=0.15', fc=s['fc'], ec=s['ec'], ls=s['ls'], lw=1.2))
    ax.text(5.3, 9.35 - 2.2 * i, lab, fontsize=8.5, va='center')
ax.text(12.0, 9.6, '状态日期 2026-11-04\nM0–M7 全部完成', fontsize=8, va='top', color='#444')

arrow(bottom(A1), top(A2)); arrow(bottom(A2), top(A3)); arrow(bottom(A3), top(A4)); arrow(bottom(A4), top(A5))
arrow(bottom(A5), (11.3, 30.3))
arrow((7.7, 28.7), (6.1, 25.8), label='否', lpos=0.4, loff=(-0.9, 0))
arrow((14.9, 28.7), (16.5, 25.8), label='是', lpos=0.4, loff=(0.9, 0))
arrow(bottom(A7a), (6.1, 18.2))
arrow(bottom(A7b), (16.5, 18.6), style='-')
arrow((16.5, 18.6), (23.5, 18.6), style='-')
arrow(right(A8, 0.5), (23.5, 16.0), style='-')

# ---------------- column 2: time loop ----------------
frame(23.5, 1.2, 40.5, 60.4, '② 时间循环（main.f90）：每个射频周期 steps_per_period 步，共 n_periods 个周期')
B0 = box(24.5, 54.4, 37.5, 4.2, 'diag_window_check（每个周期开始）', '最后 n_avg 个周期进入测量窗口：清零窗口累加器，记下 KE、W_field、t 的起点；'
         '窗口内每步多做 deposit_species_moments、J·E、直方图、x–t 累加', 'done')
frame(24.5, 17.0, 37.5, 36.3, '每一步 dt（时间层：步前 x、φ、E 在 t_n，v 在 t_{n-1/2}）', color='#2c5aa0')
B1 = box(25.5, 46.0, 35.5, 5.4, 'push: push_species(1), push_species(2)', 'E_p = 线性插值；v_x += (q/m)E_p dt；x += v_x dt（v_y, v_z 不变）。x<0 或 x≥L：记能量、计数、电荷后删除（末尾粒子换入，O(1)）。\n'
         '同时累加：做功 Σ q E_p·½(v⁻+v⁺)dt（= 动能增量，精确恒等式；窗口内按旧位置沉积成加热剖面）、Ramo 电流 Σ q v_x、壁面动能。', 'done', lfs=7.7)
B2 = box(25.5, 39.6, 35.5, 5.6, 'mcc: mcc_electrons', '零碰撞法：R₁ < P = 1−e^(−ν_max dt)（1.6 %）才算；按 ν_k(ε)=n_gσ_k v 的累积区间选过程，其余为空碰撞。\n'
         '弹性 ε\'=ε[1−(2m/M)(1−cosχ)]；激发 ε−19.82/20.61；电离 (ε−24.59)/2 两电子平分，新 He⁺ 取 T_g；方向球面均匀。', 'done', lfs=7.7)
B3 = box(25.5, 33.4, 35.5, 5.4, 'mcc: mcc_ions', '抽靶原子速度 v_n~Maxwell(T_g)；g=v_i−v_n，ε_rel=½μg²（μ=M/2）；ν_k = n_g σ_k(ε_rel)|g|。\n'
         '各向同性：v_i\' = v_cm + ½|g| n̂；背向散射（=电荷交换）：v_i\' = v_n。能量账本：非弹性、弹性、新离子热能、离子↔气体。', 'done', lfs=7.7)
B4 = box(25.5, 27.6, 35.5, 5.0, 'deposit: deposit_rho（窗口内 + deposit_species_moments）', 'CIC：j=⌊x/dx⌋, s=x/dx−j；ρ_j += qw(1−s)/dx，ρ_{j+1} += qw·s/dx；MPI_ALLREDUCE 求和（每步一次，nx+1 个数）。\n'
         '窗口内同样沉积 n_e、n_i、n v_x、能量密度（rank 局部累加，结束时归约一次）。', 'done', lfs=7.7)
B5 = box(25.5, 22.4, 35.5, 4.4, 'poisson: poisson_solve(rho, V(t_{n+1}), 0)', '(φ_{j−1}−2φ_j+φ_{j+1})/dx² = −ρ_j/ε0，φ₀=V₀sin(2πf t)，φ_nx=0；Thomas 回代 O(nx)；\n'
         'E_j = −(φ_{j+1}−φ_{j−1})/2dx，两端二阶单侧差分。', 'done', lfs=7.7)
B6 = box(25.5, 17.7, 35.5, 4.0, 'diag: diag_step', '离散高斯定律残差；窗口内：J·½(E_old+E_new)、⟨V·I⟩（I=(w/L)Σqv_x）、撞墙能量直方图（IEDF/EEDF）、'
         'x–t 相位格累加、中心 10 % 的 EEPF。', 'done', lfs=7.7)
B7 = box(24.5, 10.2, 37.5, 5.2, 'diag: diag_period（每个周期结束）', 'MPI_REDUCE 周期统计 → history.dat 一行：N_e, N_i, ⟨n⟩, 中心 n, Γ_i/Γ_e 左右, 平均撞击能, P_e, P_i, P_VI, KE, W_field,'
         '\n壁面累计能量, 累计做功, 高斯残差, 墙钟, T_x/T_⊥（e, i）；NaN 检查；进度打印', 'done', lfs=7.7)
B8 = box(24.5, 4.6, 37.5, 4.4, 'restart_write / 停止判断', '每 restart_every 周期写 restart_*.bin；periods_this_run 到数或 max_wall_minutes 到时：写续算文件后干净退出；'
         '\n否则 iperiod < n_periods 继续', 'done', lfs=7.7)
arrow((23.5, 16.0), (23.5, 56.5), style='-')
arrow((23.5, 56.5), (24.5, 56.5))
arrow(bottom(B0), top(B1)); arrow(bottom(B1), top(B2)); arrow(bottom(B2), top(B3)); arrow(bottom(B3), top(B4))
arrow(bottom(B4), top(B5)); arrow(bottom(B5), top(B6))
# step loop (left side inside the step frame)
arrow((25.5, 19.7), (24.9, 19.7), style='-', color='#2c5aa0')
arrow((24.9, 19.7), (24.9, 48.7), style='-', color='#2c5aa0')
arrow((24.9, 48.7), (25.5, 48.7), color='#2c5aa0', label='下一步', lpos=0.5, loff=(-1.2, -1.3))
arrow(bottom(B6), top(B7), label='周期结束', lpos=0.5, loff=(5.0, 0))
arrow(bottom(B7), top(B8))
# period loop (right side)
arrow(right(B8, 0.7), (63.0, 7.7), style='-', color='#2e7d32')
arrow((63.0, 7.7), (63.0, 56.5), style='-', color='#2e7d32')
arrow((63.0, 56.5), (62.0, 56.5), color='#2e7d32')
ax.text(63.3, 31, 'iperiod < n_periods：下一个周期', fontsize=7.5, color='#2e7d32', rotation=90, va='center', ha='left')

# ---------------- column 3: finish, tools, validation ----------------
frame(65.5, 1.2, 34, 60.4, '③ 结束、分析与验证')
C1 = box(66.5, 55.3, 32, 4.0, 'diag: diag_final（iperiod == n_periods）', '归约窗口累加器（MPI_REDUCE 到 rank 0）；时间平均；端点节点乘 2（半格体积）', 'done')
C2 = box(66.5, 47.0, 32, 7.3, '输出文件（全部 ASCII，# 表头带单位，numpy.loadtxt 直接读）', 'profiles.dat：n_e, n_i, φ, E, ⟨ε_e⟩, 电离率, ⟨J_e·E⟩, ⟨J_i·E⟩, Γ_i, ⟨ε_x,e⟩, ⟨ε_i⟩\n'
         'iedf_left/right, eedf_left/right：壁面通量能谱  eepf.dat：F(ε), F/√ε（中心 10 %）\nxt_phi/ne/ni/Je/Ji/PeE/ioniz.dat：n_phase × (nx+1) 相位分辨图\n'
         'info.txt 结果段：中心密度、T_eff、ωpe·dt、Γ、撞击能、能量平衡三把尺 P_part / P_JE / P_VI、P_numerical', 'done', lfs=7.6)
C3 = box(66.5, 40.4, 32, 5.6, 'tools/（全部就位）', 'compare_turner 基准叠画   convergence T15   dual_frequency T17   boltzmann2term + compare_swarm T12\n'
         'check_run 13 项   check_diag T16 23 项   plot_xt / plot_power / plot_profiles / plot_dist / plot_history   check_sheath   plot_initial', 'done', lfs=7.6)
C5 = box(66.5, 25.0, 32, 14.4, '验证已通过（make check + 物理标定 + 基准 + 收敛 + 诊断自洽 + 双频 + 电子群）', 'T1 泊松 1e-14、二阶收敛、高斯定律 1e-13   T2 CIC 沉积 1e-16、端点半格、1 与 2 rank 一致\n'
         'T3 随机数 + check-rng（2e7 粒子尾巴对称，6000 样本 z 值 rms 1.00）   T8 续算逐位相同（含碰撞）\n'
         'T4 截面查表 = 原表（电子 1e-15，离子 1e-7）；ν_max 8.654e7 / 2.793e7   T10 Langmuir 频率差 3e-4\n'
         'T5 MCC 频率：10/30/100 eV 电子四过程、5 eV 离子两过程 = n_gσv×0.992，|z| < 1.7\n'
         'T6 运动学：弹性/激发/电离能量关系 1e-15、各向同性 χ² 正常、新离子 T_g、离子质心系不变量\n'
         'T7 弛豫：离子热化到 300 K 偏差 0.05 %；电子弹性冷却与 ODE 差 0.02 %   T11 无碰撞鞘层\n'
         'T13 Case 1 自检 + T14 全尺寸：中心 n_i 1.432e14 vs 1.403e14（+2.1 %），内点 rms 1.5 %，4 核 170 s\n'
         'T15 收敛：dt/2、nx 256、nppc 1024 与基线差 ≤ 2 倍换种子噪声（地板 1.4 %）；nppc 256 略超\n'
         'T16 自洽 23 项：能谱积分 = 通量、x–t 积分 = info（1e-7）、总电流沿 x 均匀 1.7 %、鞘层相位、加热在鞘层边缘、EEPF 双温\n'
         'T17 双频：V_lf = 50/100/150 V，E_max 240 → 298 eV ≈ 0.5(V0+V_lf)，密度 −9/−26/−45 %；200 V 鞘层填满间隙熄灭\n'
         'T12 电子群：PIC vs 两项 Boltzmann（同截面）20–200 Td 平均能量、漂移速度差 0.1–1 %，α/N 2–5 %', 'done', lfs=7.4)
C6 = box(66.5, 20.6, 32, 3.4, '待做', '全部里程碑 M0–M7 完成；剩余：README / 面试材料整理，HANDOVER 更新，可选 Case 2', 'todo', lfs=7.6)
C7 = box(66.5, 12.2, 32, 7.2, 'data/（已就位）', 'he/electron_scattering, excitation_1 (19.82 eV), excitation_2 (20.61 eV), ionization (24.59 eV)：e–He，Biagi 7.1\n'
         'he/ion_scattering（各向同性）, ion_back_scatter（背向 = 电荷交换）：He⁺–He，Phelps，能量变量 ε_rel\n'
         'he/charge_exchange.dat：基准不用。格式 energy[eV] sigma[m²]，0.01 eV 间距，到 749.99 eV，超出取末值\n'
         'turner2013_case1_ni.dat：129 节点参考 n_i(x)，中心 1.405e14 m⁻³，端点 3.6e13', 'data', lfs=7.6)
C8 = box(66.5, 4.6, 32, 6.4, '运行方式与分段', 'cd runs/case1 && mpirun --allow-run-as-root -np 2 ../../ccp1d      （输入 ccp1d.in，输出同目录）\n'
         'Case 1：128 格，400 步/周期，1280 周期 = 512000 步，含碰撞 4 核实测 204 s（0.16 s/周期）\n'
         '长算例后台跑，max_wall_minutes 兜底，restart = .true. 接着跑；关键结果拷到 results/ 进 git', 'done', lfs=7.6)
arrow(right(B8, 0.25), (64.8, 5.7), style='-', color='#333')
arrow((64.8, 5.7), (64.8, 57.3), style='-', color='#333')
arrow((64.8, 57.3), (66.5, 57.3), label='周期数到', lpos=0.0, loff=(0, 1.3))
arrow(bottom(C1), top(C2))
arrow((74.1, 47.0), (74.1, 46.0))
arrow((90.8, 47.0), (90.8, 46.0))

fig.savefig(os.path.join(out, 'program_flow.png'), dpi=110, bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(out, 'program_flow.svg'), bbox_inches='tight', facecolor='white')
print('wrote', os.path.join(out, 'program_flow.png'))

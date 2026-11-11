#!/bin/sh
# make_study_pack.sh [out.pdf]: bundle the ccp1d documents and all figures into one PDF for offline study.
# Needs pandoc, node + playwright (Chromium) and a CJK font (wqy-microhei).  Run from ccp1d/.
set -e
OUT=${1:-doc/ccp1d_study_pack.pdf}
TMP=$(mktemp -d)
python3 - "$TMP" <<'PY'
import glob, os, re, sys
tmp = sys.argv[1]
docs = [('学习计划', 'doc/study_plan.md'), ('程序结构与算法讲解', 'doc/program_flow.md'), ('代码说明', 'doc/code_guide.md'),
        ('初始条件讲解', 'doc/initial_conditions.md'), ('图册与数据索引', 'doc/figures.md'), ('与经典理论的对照', 'doc/theory.md'),
        ('实施计划与测试矩阵', 'PLAN.md'), ('进度日志', 'PROGRESS.md')]
figs = sorted(glob.glob('doc/fig/*.png')) + sorted(glob.glob('results/*.png')) + sorted(glob.glob('results/*/*.png')) + sorted(glob.glob('results/*/*/*.png'))
# downscale copies of the figures for the PDF
try:
    from PIL import Image
    def copy(src):
        dst = os.path.join(tmp, src.replace('/', '__'))
        im = Image.open(src)
        if im.width > 1400:
            im = im.resize((1400, int(im.height * 1400 / im.width)), Image.LANCZOS)
        im.save(dst, optimize=True)
        return dst
except ImportError:
    import shutil
    def copy(src):
        dst = os.path.join(tmp, src.replace('/', '__')); shutil.copy(src, dst); return dst
out = ['---', 'title: ccp1d 学习资料包', 'subtitle: 一维有界静电 PIC/MCC 代码（Turner 2013 氦气基准）的说明、讲解、图册与学习计划',
       'date: 2026-11-11', 'lang: zh-CN', '---', '', '来源：https://github.com/fanyumeng666/maccpic （目录 `ccp1d/`）。本文件由 `tools/make_study_pack.sh` 生成；各章节分别对应 `doc/` 下的 Markdown 文件，图在 `doc/fig/` 和 `results/`。', '']
for title, path in docs:
    txt = open(path, encoding='utf-8').read()
    # demote headings by one level so that each document is a chapter
    txt = re.sub(r'^(#+) ', lambda m: '#' + m.group(1) + ' ', txt, flags=re.M)
    # embedded images: point them at the copies
    def fix(m):
        src = m.group(2)
        cand = os.path.normpath(os.path.join(os.path.dirname(path), src))
        if os.path.exists(cand):
            return '![%s](%s)' % (m.group(1), copy(cand))
        return m.group(0)
    txt = re.sub(r'!\[([^\]]*)\]\(([^)]+)\)', fix, txt)
    out += ['', '\\newpage', '', '# %s（%s）' % (title, path), '', txt, '']
out += ['', '\\newpage', '', '# 图册（全部结果图）', '']
for f in figs:
    out += ['', '## %s' % f, '', '![%s](%s)' % (f, copy(f)), '']
open(os.path.join(tmp, 'pack.md'), 'w', encoding='utf-8').write('\n'.join(out))
print('figures:', len(figs))
PY
cat > "$TMP/style.css" <<'CSS'
body { font-family: "WenQuanYi Micro Hei", "Noto Sans CJK SC", sans-serif; font-size: 10.5pt; line-height: 1.45; max-width: none; margin: 0; }
h1 { font-size: 18pt; border-bottom: 2px solid #444; padding-bottom: 3px; margin-top: 1.4em; page-break-before: always; }
h1.title { page-break-before: avoid; font-size: 24pt; border: none; }
h2 { font-size: 14pt; margin-top: 1.2em; } h3 { font-size: 12pt; }
table { border-collapse: collapse; font-size: 9pt; margin: 0.6em 0; page-break-inside: auto; }
th, td { border: 1px solid #999; padding: 2px 5px; vertical-align: top; }
tr { page-break-inside: avoid; }
code { font-family: "DejaVu Sans Mono", monospace; font-size: 8.5pt; background: #f3f3f3; padding: 0 2px; }
pre { background: #f3f3f3; padding: 6px; font-size: 8.5pt; white-space: pre-wrap; }
img { max-width: 100%; page-break-inside: avoid; }
blockquote { border-left: 3px solid #bbb; margin-left: 0; padding-left: 10px; color: #333; }
CSS
pandoc "$TMP/pack.md" -s --embed-resources --toc --toc-depth=2 --css "$TMP/style.css" -f markdown+pipe_tables+tex_math_dollars -o "$TMP/pack.html"
node -e "
const { chromium } = require('playwright');
(async () => { const b = await chromium.launch(); const p = await b.newPage();
  await p.goto('file://' + process.argv[1]); await p.emulateMedia({ media: 'print' });
  await p.pdf({ path: process.argv[2], format: 'A4', printBackground: true, margin: { top: '16mm', bottom: '16mm', left: '15mm', right: '15mm' },
    displayHeaderFooter: true, headerTemplate: '<div></div>', footerTemplate: '<div style=\"font-size:8px;width:100%;text-align:center;color:#666\">ccp1d study pack  <span class=\"pageNumber\"></span> / <span class=\"totalPages\"></span></div>' });
  await b.close(); })().catch(e => { console.error(e.message); process.exit(1); });
" "$TMP/pack.html" "$OUT"
rm -rf "$TMP"
ls -la "$OUT"

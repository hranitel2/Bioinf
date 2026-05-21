#!/usr/bin/env python3
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import SeqIO
import os

for fname in ["chr1_1mb.fa", "chr1.fa"]:
    if os.path.exists(fname):
        record = SeqIO.read(fname, "fasta")
        break
else:
    raise FileNotFoundError("Не найден файл chr1_1mb.fa или chr1.fa")

seq = str(record.seq[:100000]).upper()
print(f"Загружена последовательность: {record.id}, длина анализируемого участка: {len(seq)}")

def compute_cpg_ratio(window):
    """
    Для окна (строка ДНК) вычисляет:
    P(CG), P(C), P(G), R(CG) = P(CG)/(P(C)*P(G)), GC%
    """
    n = len(window)
    if n == 0:
        return 0, 0, 0, 0, 0.0
    cg_count = window.count("CG")
    c_count = window.count("C")
    g_count = window.count("G")
    p_cg = cg_count / (n - 1) if n > 1 else 0
    p_c = c_count / n
    p_g = g_count / n
    if p_c > 0 and p_g > 0:
        r_cg = p_cg / (p_c * p_g)
    else:
        r_cg = 0.0
    gc_percent = (c_count + g_count) / n * 100
    return p_cg, p_c, p_g, r_cg, gc_percent

window_size = 200
step = 10  
positions = []
r_cg_values = []
gc_values = []
islands = []  # координаты островков (start, end)

for start in range(0, len(seq) - window_size + 1, step):
    end = start + window_size
    window = seq[start:end]
    p_cg, p_c, p_g, r_cg, gc_percent = compute_cpg_ratio(window)
    positions.append(start + window_size // 2)
    r_cg_values.append(r_cg)
    gc_values.append(gc_percent)
    if r_cg > 0.6 and gc_percent > 50:
        islands.append((start, end))

# Объединим перекрывающиеся или смежные островки
if islands:
    merged_islands = []
    cur_start, cur_end = islands[0]
    for start, end in islands[1:]:
        if start <= cur_end:  # перекрываются или соприкасаются
            cur_end = max(cur_end, end)
        else:
            merged_islands.append((cur_start, cur_end))
            cur_start, cur_end = start, end
    merged_islands.append((cur_start, cur_end))
else:
    merged_islands = []

print(f"\nНайдено CpG-островков (R(CG)>0.6, GC%>50%): {len(merged_islands)}")
for i, (s, e) in enumerate(merged_islands, 1):
    print(f"  Островок {i}: {s}-{e} (длина {e-s} п.н.)")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

ax1.plot(positions, r_cg_values, color='steelblue', linewidth=0.8, label='R(CG)')
ax1.axhline(y=0.6, color='red', linestyle='--', linewidth=1.2, label='Порог 0.6')
ax1.set_ylabel('R(CG)', fontsize=12)
ax1.set_title('Распределение R(CG) и GC% вдоль chr1 (первые 100 тыс. п.н.)', fontsize=14, fontweight='bold')
ax1.legend(loc='upper right')
ax1.grid(True, alpha=0.3, linestyle='--')

ax2.plot(positions, gc_values, color='darkgreen', linewidth=0.8, label='GC%')
ax2.axhline(y=50, color='red', linestyle='--', linewidth=1.2, label='Порог 50%')
ax2.set_xlabel('Позиция в геноме (п.н.)', fontsize=12)
ax2.set_ylabel('GC%', fontsize=12)
ax2.legend(loc='upper right')
ax2.grid(True, alpha=0.3, linestyle='--')

plt.tight_layout()
plt.savefig("cpg_islands_plot.png", dpi=150, bbox_inches='tight')
print("График сохранён как cpg_islands_plot.png")

# Сохраняем координаты островков в файл
with open("cpg_islands_results.txt", "w") as f:
    f.write("# CpG-островки (R(CG)>0.6 и GC%>50%)\n")
    f.write(f"# Всего найдено: {len(merged_islands)}\n")
    for i, (s, e) in enumerate(merged_islands, 1):
        f.write(f"Островок {i}: {s}-{e} (длина {e-s} п.н.)\n")
print("Результаты сохранены в cpg_islands_results.txt")

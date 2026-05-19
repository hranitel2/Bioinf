#!/usr/bin/env python3
import matplotlib
from Bio import SeqIO
import matplotlib.pyplot as plt

# Загружаем последовательности
neanderthal = SeqIO.read("neanderthal.fasta", "fasta")
human = SeqIO.read("human.fasta", "fasta")

seq_nea = str(neanderthal.seq).upper()
seq_hum = str(human.seq).upper()

print(f"Неандерталец: {len(seq_nea)} п.н.")
print(f"Человек:      {len(seq_hum)} п.н.")

# Мягкие параметры
WINDOW = 30
THRESHOLD = 0.60  # снизили порог до 60%
STEP = 5          # меньше шаг — больше точек

dots_x, dots_y = [], []

for i in range(0, len(seq_nea) - WINDOW, STEP):
    window_nea = seq_nea[i:i+WINDOW]
    for j in range(0, len(seq_hum) - WINDOW, STEP):
        window_hum = seq_hum[j:j+WINDOW]
        matches = sum(1 for a, b in zip(window_nea, window_hum) if a == b)
        similarity = matches / WINDOW
        if similarity >= THRESHOLD:
            dots_x.append(i + WINDOW // 2)
            dots_y.append(j + WINDOW // 2)

print(f"Точек на графике: {len(dots_x)}")

# График
fig, ax = plt.subplots(1, 1, figsize=(10, 10))
ax.scatter(dots_x, dots_y, s=0.3, c='black', alpha=0.5)
ax.set_xlabel(f"Неандерталец (NC_011137.1) — {len(seq_nea)} п.н.", fontsize=12)
ax.set_ylabel(f"Человек (NC_012920.1) — {len(seq_hum)} п.н.", fontsize=12)
ax.set_title("DotPlot синтении митохондриальных геномов\nНеандерталец vs Современный человек",
             fontsize=14, fontweight='bold')

ax.grid(True, alpha=0.3, linestyle='--')
ax.set_xlim(0, len(seq_nea))
ax.set_ylim(0, len(seq_hum))

min_len = min(len(seq_nea), len(seq_hum))
ax.plot([0, min_len], [0, min_len], 'r--', alpha=0.5, linewidth=1, label='Диагональ (полная синтения)')
ax.legend(loc='upper left')

plt.tight_layout()
plt.savefig("dotplot_mitochondria.png", dpi=150, bbox_inches='tight')
print("График сохранён как dotplot_mitochondria.png")
plt.show()

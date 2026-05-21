#!/usr/bin/env python3
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import SeqIO

record = SeqIO.read("sequence.fasta", "fasta")
seq = str(record.seq).upper()

nucs = ['A', 'C', 'G', 'T']
nuc_to_idx = {n: i for i, n in enumerate(nucs)}

counts = np.zeros((4, 4), dtype=int)
for k in range(len(seq) - 1):
    i = nuc_to_idx[seq[k]]
    j = nuc_to_idx[seq[k+1]]
    counts[i, j] += 1

row_sums = counts.sum(axis=1, keepdims=True)
P = counts / row_sums

print("Исходная матрица переходов P:")
print(P)

def generate_sequence(P, length, start_nuc):
    """Генерирует последовательность длины length, начиная с start_nuc."""
    nucs = ['A', 'C', 'G', 'T']
    nuc_to_idx = {n: i for i, n in enumerate(nucs)}
    seq = [start_nuc]
    current_idx = nuc_to_idx[start_nuc]
    for _ in range(length - 1):
        next_idx = np.random.choice(4, p=P[current_idx])
        seq.append(nucs[next_idx])
        current_idx = next_idx
    return ''.join(seq)

np.random.seed(123)  # для воспроизводимости
num_seq = 10
length_seq = 1000
all_empirical_P = []  # список матриц переходов

print(f"\nГенерирую {num_seq} последовательностей длиной {length_seq}...")
for s in range(num_seq):
    # Стартовый нуклеотид выбираем случайно с равной вероятностью
    start = np.random.choice(nucs)
    gen_seq = generate_sequence(P, length_seq, start)
    # Подсчёт переходов в этой последовательности
    emp_counts = np.zeros((4, 4), dtype=int)
    for k in range(len(gen_seq) - 1):
        i = nuc_to_idx[gen_seq[k]]
        j = nuc_to_idx[gen_seq[k+1]]
        emp_counts[i, j] += 1
    # Нормализация строк
    emp_row_sums = emp_counts.sum(axis=1, keepdims=True)
    # На случай, если какой-то нуклеотид ни разу не встретился, заменяем 0 на 1 чтобы избежать деления на 0
    emp_row_sums[emp_row_sums == 0] = 1
    emp_P = emp_counts / emp_row_sums
    all_empirical_P.append(emp_P)

# Усреднённая эмпирическая матрица
avg_empirical_P = np.mean(all_empirical_P, axis=0)
print("\nСредняя эмпирическая матрица переходов (по 10 последовательностям):")
print(np.round(avg_empirical_P, 4))

# Сравнение
diff = avg_empirical_P - P
print("\nРазница (средняя эмпирическая - исходная P):")
print(np.round(diff, 4))
fro_norm = np.linalg.norm(diff, 'fro')
print(f"Норма Фробениуса разности: {fro_norm:.6f}")

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Исходная матрица P
im1 = axes[0].imshow(P, cmap='YlOrRd', vmin=0, vmax=0.5)
axes[0].set_title('Исходная матрица P', fontsize=12, fontweight='bold')
axes[0].set_xticks(range(4))
axes[0].set_yticks(range(4))
axes[0].set_xticklabels(nucs)
axes[0].set_yticklabels(nucs)
axes[0].set_xlabel('К')
axes[0].set_ylabel('От')
for i in range(4):
    for j in range(4):
        axes[0].text(j, i, f'{P[i,j]:.3f}', ha='center', va='center', fontsize=9)

# Средняя эмпирическая матрица
im2 = axes[1].imshow(avg_empirical_P, cmap='YlOrRd', vmin=0, vmax=0.5)
axes[1].set_title('Средняя эмпирическая P (10 генер.)', fontsize=12, fontweight='bold')
axes[1].set_xticks(range(4))
axes[1].set_yticks(range(4))
axes[1].set_xticklabels(nucs)
axes[1].set_yticklabels(nucs)
axes[1].set_xlabel('К')
for i in range(4):
    for j in range(4):
        axes[1].text(j, i, f'{avg_empirical_P[i,j]:.3f}', ha='center', va='center', fontsize=9)

fig.colorbar(im2, ax=axes, shrink=0.8)
plt.suptitle('Сравнение матриц переходов', fontsize=14)
plt.tight_layout()
plt.savefig("heatmap_comparison.png", dpi=150, bbox_inches='tight')
print("Тепловые карты сохранены как heatmap_comparison.png")

# Сохранение результатов
with open("simulation_results.txt", "w") as f:
    f.write("# Симуляция марковской цепи\n")
    f.write(f"# Сгенерировано {num_seq} последовательностей длиной {length_seq}\n")
    f.write("# Исходная матрица P:\n")
    f.write(np.array2string(P, precision=4, suppress_small=True))
    f.write("\n\n# Средняя эмпирическая матрица:\n")
    f.write(np.array2string(avg_empirical_P, precision=4, suppress_small=True))
    f.write(f"\n\n# Норма Фробениуса разности: {fro_norm:.6f}\n")
    f.write("# ВЫВОД: Средняя эмпирическая матрица близка к исходной P.\n")
print("Результаты сохранены в simulation_results.txt")

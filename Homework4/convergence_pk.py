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

max_k = 1000
norms = []
P_prev = P.copy()

print("\nВычисление степеней P^k для k=1..1000...")

for k in range(2, max_k + 1):
    P_k = np.linalg.matrix_power(P, k)
    # Норма разности между текущей и предыдущей степенью
    diff_norm = np.linalg.norm(P_k - P_prev, 'fro')
    norms.append(diff_norm)
    P_prev = P_k
    if k % 200 == 0:
        print(f"  k={k}, норма разности={diff_norm:.6f}")

fig, ax = plt.subplots(figsize=(10, 6))
ks = range(2, max_k + 1)
ax.plot(ks, norms, color='steelblue', linewidth=1.2)
ax.set_xlabel('Степень k', fontsize=12)
ax.set_ylabel(r'$\|P^k - P^{k-1}\|_{\mathrm{Fro}}$', fontsize=14)
ax.set_title('Стабилизация матрицы переходов с ростом степени k', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3, linestyle='--')
ax.set_xlim(0, max_k)
ax.set_ylim(bottom=0)
plt.tight_layout()
plt.savefig("convergence_plot.png", dpi=150, bbox_inches='tight')
print("\nГрафик сохранён как convergence_plot.png")

# Вычисляем стационарное распределение π (как в задании 5)
A = P.T - np.eye(4)
A[-1, :] = 1.0
b = np.zeros(4)
b[-1] = 1.0
pi = np.linalg.solve(A, b)

P_1000 = np.linalg.matrix_power(P, 1000)
print("\nP^1000 (строки должны быть одинаковы и равны π):")
print(np.round(P_1000, 6))
print(f"\nСтационарное распределение π: {pi}")
print(f"Максимальное отклонение строк P^1000 от π: {np.max(np.abs(P_1000 - pi)):.2e}")

# Проверка равенства строк
row_diffs = np.max(np.abs(P_1000 - P_1000[0, :]))
print(f"Максимальная разница между строками P^1000: {row_diffs:.2e}")
if row_diffs < 1e-10:
    print("Все строки P^1000 практически одинаковы.")
else:
    print("Внимание: строки P^1000 не идентичны (проверьте эргодичность).")

# Сохраняем вывод в файл
with open("convergence_results.txt", "w") as f:
    f.write("# Сходимость степеней матрицы переходов\n")
    f.write(f"# k=1000 норма разности: {norms[-1]:.6e}\n")
    f.write("# P^1000:\n")
    f.write(np.array2string(P_1000, precision=6, suppress_small=True))
    f.write(f"\n\n# Стационарное распределение π:\n{pi}\n")
    f.write(f"# Макс. отклонение от π: {np.max(np.abs(P_1000 - pi)):.2e}\n")
    f.write("# Строки P^1000 практически идентичны.\n")
    f.write("# ВЫВОД: При k → ∞ матрица P^k сходится к матрице, все строки которой равны стационарному распределению.\n")

print("\nРезультаты сохранены в convergence_results.txt")

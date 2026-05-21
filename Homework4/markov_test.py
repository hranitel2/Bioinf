#!/usr/bin/env python3
import numpy as np
from scipy.stats import chi2
from Bio import SeqIO
import os
for fname in ["chr1_1mb.fa", "chr1.fa"]:
    if os.path.exists(fname):
        record = SeqIO.read(fname, "fasta")
        break
else:
    raise FileNotFoundError("chr1_1mb.fa не найден")

seq_real = str(record.seq[:100000]).upper()
seq_real = ''.join([c for c in seq_real if c in 'ACGT'])
L = len(seq_real)
print(f"Длина реальной последовательности: {L} нуклеотидов")

nucs = ['A', 'C', 'G', 'T']
nuc_to_idx = {n: i for i, n in enumerate(nucs)}

def count_dinucleotide_matrix(seq):
    counts = np.zeros((4, 4), dtype=int)
    for k in range(len(seq) - 1):
        i = nuc_to_idx[seq[k]]
        j = nuc_to_idx[seq[k+1]]
        counts[i, j] += 1
    return counts

def nucleotide_frequencies(seq):
    freq = np.array([seq.count(n) / len(seq) for n in nucs])
    return freq

def expected_matrix(freq, N):
    return N * np.outer(freq, freq)

def chi2_test(O, E):
    """Вычисляет хи-квадрат и p-value с df=9."""
    # Избегаем деления на ноль: если E_ij=0, то и O_ij должно быть 0, вклад 0
    mask = E > 0
    chi2_val = np.sum((O[mask] - E[mask])**2 / E[mask])
    p_value = 1 - chi2.cdf(chi2_val, df=9)
    return chi2_val, p_value

print("\n" + "="*60)
print("РЕАЛЬНАЯ ПОСЛЕДОВАТЕЛЬНОСТЬ (chr1:1-100000)")
print("="*60)
O_real = count_dinucleotide_matrix(seq_real)
freq_real = nucleotide_frequencies(seq_real)
N_real = len(seq_real) - 1
E_real = expected_matrix(freq_real, N_real)

print("Наблюдаемые частоты динуклеотидов (O):")
print(O_real)
print("\nОжидаемые частоты (E):")
print(np.round(E_real, 2))
chi2_real, p_real = chi2_test(O_real, E_real)
print(f"\nχ² = {chi2_real:.2f}, p-value = {p_real:.6f}")

np.random.seed(123)
seq_random = ''.join(np.random.choice(nucs, size=L, p=freq_real))
print("\n" + "="*60)
print("СЛУЧАЙНАЯ ПОСЛЕДОВАТЕЛЬНОСТЬ (0-й порядок, те же частоты нуклеотидов)")
print("="*60)
O_rand = count_dinucleotide_matrix(seq_random)
freq_rand = nucleotide_frequencies(seq_random)
N_rand = len(seq_random) - 1
E_rand = expected_matrix(freq_rand, N_rand)

print("Наблюдаемые частоты динуклеотидов (O):")
print(O_rand)
print("\nОжидаемые частоты (E):")
print(np.round(E_rand, 2))
chi2_rand, p_rand = chi2_test(O_rand, E_rand)
print(f"\nχ² = {chi2_rand:.2f}, p-value = {p_rand:.6f}")

print("\n" + "="*60)
print("ВЫВОД")
print("="*60)
if p_real < 0.05:
    print(f"Для реального участка chr1 p-value = {p_real:.6f} < 0.05, следовательно:")
    print("  нулевая гипотеза о независимости соседних нуклеотидов (0-й порядок) отвергается.")
    print("  Наблюдается значимая зависимость между соседними нуклеотидами → необходима модель 1-го порядка.")
else:
    print(f"Для реального участка p-value = {p_real:.6f} >= 0.05 — нет оснований отвергать 0-й порядок.")

if p_rand < 0.05:
    print(f"Для случайной последовательности p-value = {p_rand:.6f} < 0.05 (что маловероятно при правильной симуляции).")
else:
    print(f"Для случайной последовательности p-value = {p_rand:.6f} >= 0.05, что ожидается, так как она построена без зависимостей (0-й порядок).")

with open("markov_order_test_results.txt", "w") as f:
    f.write("# Тест на марковость: 0-й vs 1-й порядок\n")
    f.write(f"# Реальная последовательность: chr1:1-100000 (длина {L})\n")
    f.write(f"# Наблюдаемые динуклеотидные частоты:\n{O_real}\n")
    f.write(f"# Ожидаемые (независимость):\n{np.round(E_real, 2)}\n")
    f.write(f"# χ² = {chi2_real:.2f}, p-value = {p_real:.6f}\n")
    if p_real < 0.05:
        f.write("# Результат: гипотеза независимости отвергается, требуется модель 1-го порядка.\n")
    else:
        f.write("# Результат: нет оснований отвергать независимость.\n")
    
    f.write(f"\n# Случайная последовательность (0-й порядок):\n")
    f.write(f"# Наблюдаемые динуклеотидные частоты:\n{O_rand}\n")
    f.write(f"# Ожидаемые:\n{np.round(E_rand, 2)}\n")
    f.write(f"# χ² = {chi2_rand:.2f}, p-value = {p_rand:.6f}\n")
    if p_rand < 0.05:
        f.write("# Результат: неожиданное отклонение, проверьте генератор.\n")
    else:
        f.write("# Результат: гипотеза независимости согласуется с данными (модель 0-го порядка адекватна).\n")

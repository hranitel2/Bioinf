#!/usr/bin/env python3
import numpy as np

# Входные данные
sites = [
    "GAGGTAAAC",
    "TCCGTAAGC",
    "CAGGTTGGA",
    "ACAGTCAGC",
    "TAGGTCAGC",
    "CAGGTCAGC",
    "CAGGTCGAT",
    "CAGGTCAGC",
    "CAGGTCAGC",
    "CAGGTTGGC"
]

NUCLEOTIDES = ['A', 'C', 'G', 'T']
NUC_TO_IDX = {nuc: i for i, nuc in enumerate(NUCLEOTIDES)}


def build_pfm(sites):
    """
    PFM[i, j] — абсолютная частота нуклеотида i в позиции j.
    """
    n_sites = len(sites)
    site_len = len(sites[0])
    pfm = np.zeros((4, site_len), dtype=float)
    
    for site in sites:
        for j, nuc in enumerate(site):
            pfm[NUC_TO_IDX[nuc], j] += 1
    
    return pfm


def pfm_to_ppm(pfm, alpha=0.1):
    """
    PPM[i, j] = (PFM[i, j] + alpha) / (N + 4*alpha)
    где N — общее число сайтов.
    """
    N = pfm[:, 0].sum()  # общее число сайтов (сумма по любому столбцу)
    ppm = (pfm + alpha) / (N + 4 * alpha)
    return ppm


def ppm_to_pwm(ppm, background):
    """
    PWM[i, j] = log2(PPM[i, j] / bg[i])
    """
    bg_array = np.array([background[nuc] for nuc in NUCLEOTIDES])
    pwm = np.log2(ppm / bg_array[:, np.newaxis])
    return pwm


def score_sequence(seq, pwm):
    total = 0.0
    for j, nuc in enumerate(seq):
        total += pwm[NUC_TO_IDX[nuc], j]
    return total


def find_extreme_sequence(pwm, mode='max'):
    """
    Находит последовательность, дающую макс/мин скор.
    """
    if mode == 'max':
        indices = np.argmax(pwm, axis=0)
    else:
        indices = np.argmin(pwm, axis=0)
    
    seq = ''.join(NUCLEOTIDES[i] for i in indices)
    return seq


print("=" * 60)
print("1. POSITION FREQUENCY MATRIX (PFM)")
print("=" * 60)
pfm = build_pfm(sites)
print("Строки: A, C, G, T")
print(pfm)

print("\n" + "=" * 60)
print("2. POSITION PROBABILITY MATRIX (PPM), α = 0.1")
print("=" * 60)
ppm = pfm_to_ppm(pfm, alpha=0.1)
np.set_printoptions(precision=4, suppress=True)
print("Строки: A, C, G, T")
print(ppm)

print("\n" + "=" * 60)
print("3. POSITION WEIGHT MATRIX (PWM)")
print("=" * 60)
background = {'A': 0.295, 'C': 0.205, 'G': 0.205, 'T': 0.295}
pwm = ppm_to_pwm(ppm, background)
print("Строки: A, C, G, T")
print(np.round(pwm, 4))

print("\n" + "=" * 60)
print("4. ЭКСТРЕМАЛЬНЫЕ СКОРЫ")
print("=" * 60)

max_seq = find_extreme_sequence(pwm, mode='max')
min_seq = find_extreme_sequence(pwm, mode='min')
max_score = score_sequence(max_seq, pwm)
min_score = score_sequence(min_seq, pwm)

print(f"\nМаксимальный скор: {max_score:.4f}")
print(f"Последовательность: {max_seq}")
print(f"\nМинимальный скор: {min_score:.4f}")
print(f"Последовательность: {min_seq}")

# Все ли последовательности в матрице дают ожидаемый скор?
print("\n" + "=" * 60)
print("5. СКОРЫ ВСЕХ ИСХОДНЫХ САЙТОВ")
print("=" * 60)
for i, site in enumerate(sites, 1):
    s = score_sequence(site, pwm)
    print(f"  {i:2d}. {site} → {s:.4f}")


import sys
from io import StringIO

output = StringIO()
sys.stdout = output

print("=" * 60)
print("1. POSITION FREQUENCY MATRIX (PFM)")
print("=" * 60)
pfm = build_pfm(sites)
print("Строки: A, C, G, T")
print(pfm)

print("\n" + "=" * 60)
print("2. POSITION PROBABILITY MATRIX (PPM), α = 0.1")
print("=" * 60)
ppm = pfm_to_ppm(pfm, alpha=0.1)
np.set_printoptions(precision=4, suppress=True)
print("Строки: A, C, G, T")
print(ppm)

print("\n" + "=" * 60)
print("3. POSITION WEIGHT MATRIX (PWM)")
print("=" * 60)
background = {'A': 0.295, 'C': 0.205, 'G': 0.205, 'T': 0.295}
pwm = ppm_to_pwm(ppm, background)
print("Строки: A, C, G, T")
print(np.round(pwm, 4))

print("\n" + "=" * 60)
print("4. ЭКСТРЕМАЛЬНЫЕ СКОРЫ")
print("=" * 60)

max_seq = find_extreme_sequence(pwm, mode='max')
min_seq = find_extreme_sequence(pwm, mode='min')
max_score = score_sequence(max_seq, pwm)
min_score = score_sequence(min_seq, pwm)

print(f"\nМаксимальный скор: {max_score:.4f}")
print(f"Последовательность: {max_seq}")
print(f"\nМинимальный скор: {min_score:.4f}")
print(f"Последовательность: {min_seq}")

print("\n" + "=" * 60)
print("5. СКОРЫ ВСЕХ ИСХОДНЫХ САЙТОВ")
print("=" * 60)
for i, site in enumerate(sites, 1):
    s = score_sequence(site, pwm)
    print(f"  {i:2d}. {site} → {s:.4f}")

sys.stdout = sys.__stdout__
with open("pwm_results.txt", "w") as f:
    f.write(output.getvalue())

print("Результаты сохранены в pwm_results.txt")

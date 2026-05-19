#!/usr/bin/env python3
import numpy as np
import matplotlib
matplotlib.use('Agg')  # для сохранения в файл
import matplotlib.pyplot as plt
from Bio.motifs import create

# ===== Часть 1: Создаём PWM из Задания 1 =====
sites = [
    "GAGGTAAAC", "TCCGTAAGC", "CAGGTTGGA",
    "ACAGTCAGC", "TAGGTCAGC", "CAGGTCAGC",
    "CAGGTCGAT", "CAGGTCAGC", "CAGGTCAGC",
    "CAGGTTGGC"
]

motif = create(sites)
background = {'A': 0.295, 'C': 0.205, 'G': 0.205, 'T': 0.295}
motif.background = background
motif.pseudocounts = 0.1

pwm = motif.pwm
print("Консенсус мотива:", motif.consensus)

# ===== Часть 2: Генератор случайных последовательностей =====
def generate_random_sequence(length=9, probs=None):
    if probs is None:
        probs = [0.295, 0.205, 0.205, 0.295]  # A, C, G, T
    
    nucs = ['A', 'C', 'G', 'T']
    return ''.join(np.random.choice(nucs, size=length, p=probs))


def score_sequence(seq, pwm):
    """Скор последовательности по PWM."""
    total = 0.0
    for j, nuc in enumerate(seq):
        total += pwm[nuc][j]
    return total


N = 100000
L = 9
print(f"\nГенерация {N:,} случайных последовательностей длины {L}...")

np.random.seed(42)  # для воспроизводимости
bg_scores = []
for i in range(N):
    seq = generate_random_sequence(L)
    s = score_sequence(seq, pwm)
    bg_scores.append(s)

    if (i + 1) % 20000 == 0:
        print(f"  Сгенерировано {i+1:,}...")

bg_scores = np.array(bg_scores)
print(f"Готово. Диапазон скоров: [{bg_scores.min():.3f}, {bg_scores.max():.3f}]")

#Гистограмма
fig, ax = plt.subplots(1, 1, figsize=(10, 6))

ax.hist(bg_scores, bins=80, color='steelblue', edgecolor='white', alpha=0.8, density=True)
ax.set_xlabel('Score', fontsize=12)
ax.set_ylabel('Плотность вероятности', fontsize=12)
ax.set_title(f'Распределение скоров PWM на фоновых последовательностях\n'
             f'(N = {N:,}, длина = {L})', fontsize=14, fontweight='bold')

threshold_pval = np.percentile(bg_scores, 100 * (1 - 1e-4))
ax.axvline(threshold_pval, color='red', linestyle='--', linewidth=2,
           label=f'Порог p-value ≈ 10⁻⁴: {threshold_pval:.3f}')
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3, linestyle='--')

plt.tight_layout()
plt.savefig("score_distribution.png", dpi=150, bbox_inches='tight')
print("\nГистограмма сохранена в score_distribution.png")

def get_pvalue(score, bg_scores):
    return np.mean(bg_scores >= score)

print("\n" + "=" * 60)
print("РАСЧЁТ ПОРОГА ДЛЯ p-value ≈ 10^-4")
print("=" * 60)

sorted_scores = np.sort(bg_scores)[::-1]  # по убыванию
k = int(N * 1e-4)  # ожидаемое число хитов выше порога
threshold_exact = sorted_scores[max(0, k - 1)]

print(f"Порог (p-value ≈ 10⁻⁴): {threshold_exact:.4f}")

pval_check = get_pvalue(threshold_exact, bg_scores)
print(f"Фактический p-value: {pval_check:.6f}")
print(f"Число фоновых последовательностей выше порога: {np.sum(bg_scores >= threshold_exact)} из {N}")

print(f"\n{'='*60}")
print("ТАБЛИЦА ПОРОГОВ ДЛЯ РАЗНЫХ УРОВНЕЙ ЗНАЧИМОСТИ")
print("=" * 60)
for alpha in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]:
    k_alpha = max(0, int(N * alpha) - 1)
    thr = sorted_scores[k_alpha]
    actual_pval = get_pvalue(thr, bg_scores)
    print(f"  α = {alpha:.0e} → порог = {thr:.4f} (фактический p-value = {actual_pval:.6f})")

with open("pvalue_results.txt", "w") as f:
    f.write(f"# Эмпирическое распределение и p-value для PWM\n")
    f.write(f"# Консенсус мотива: {motif.consensus}\n")
    f.write(f"# N = {N:,} фоновых последовательностей длины {L}\n")
    f.write(f"#\n")
    f.write(f"# Статистика распределения:\n")
    f.write(f"#   Среднее: {bg_scores.mean():.4f}\n")
    f.write(f"#   Стд: {bg_scores.std():.4f}\n")
    f.write(f"#   Мин: {bg_scores.min():.4f}\n")
    f.write(f"#   Макс: {bg_scores.max():.4f}\n")
    f.write(f"#\n")
    f.write(f"# Порог для p-value ≈ 10^-4: {threshold_exact:.4f}\n")
    f.write(f"# Фактический p-value: {pval_check:.6f}\n")
    f.write(f"#\n")
    f.write(f"# Таблица порогов:\n")
    for alpha in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]:
        k_alpha = max(0, int(N * alpha) - 1)
        thr = sorted_scores[k_alpha]
        actual_pval = get_pvalue(thr, bg_scores)
        f.write(f"#   α = {alpha:.0e} → порог = {thr:.4f} (p-value = {actual_pval:.6f})\n")

print("\nРезультаты сохранены в pvalue_results.txt")

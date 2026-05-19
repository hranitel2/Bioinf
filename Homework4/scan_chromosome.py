#!/usr/bin/env python3
from Bio import SeqIO
from Bio.motifs import create
from Bio.motifs.matrix import PositionWeightMatrix
import numpy as np

#Создание мотива
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

print("Консенсус мотива:", motif.consensus)

pssm = motif.pssm
print("PSSM создан\n")

#Загрузка chr1 из файла
print("Загрузка chr1.fa...")
record = SeqIO.read("chr1_1mb.fa", "fasta")
sequence = str(record.seq[:1000000]).upper()  # первые 1 млн нуклеотидов

print(f"Загружено {len(sequence):,} нуклеотидов")
print(f"Первые 60: {sequence[:60]}\n")

#Поиск на прямой цепи
threshold = 5.0

print("=" * 60)
print(f"ПОИСК НА ПРЯМОЙ ЦЕПИ (threshold ≥ {threshold})")
print("=" * 60)

hits_forward = []
for pos, score in pssm.search(sequence, threshold=threshold):
    site_seq = sequence[pos:pos + len(motif)]
    hits_forward.append((pos, "+", score, site_seq))

print(f"Найдено хитов: {len(hits_forward)}")
for pos, strand, score, site in hits_forward[:20]:
    print(f"  Позиция {pos:>8} [{strand}] скор={score:.3f}  {site}")

#Обратная комплементарная цепь
print(f"\n{'='*60}")
print(f"ПОИСК НА ОБРАТНОЙ ЦЕПИ (threshold ≥ {threshold})")
print("=" * 60)

complement = str.maketrans('ACGT', 'TGCA')
reverse_complement = sequence.translate(complement)[::-1]

hits_reverse = []
for pos, score in pssm.search(reverse_complement, threshold=threshold):
    original_pos = len(sequence) - (pos + len(motif))
    site_seq = reverse_complement[pos:pos + len(motif)]
    hits_reverse.append((original_pos, "-", score, site_seq))

print(f"Найдено хитов: {len(hits_reverse)}")
for pos, strand, score, site in hits_reverse[:20]:
    print(f"  Позиция {pos:>8} [{strand}] скор={score:.3f}  {site}")

#Сохранение
all_hits = hits_forward + hits_reverse
all_hits.sort(key=lambda x: x[0])

print(f"\n{'='*60}")
print(f"ВСЕГО ХИТОВ: {len(all_hits)} (порог {threshold})")
print("=" * 60)

with open("chr1_hits.txt", "w") as f:
    f.write(f"# Поиск сайтов связывания в chr1:1-1,000,000\n")
    f.write(f"# Порог скор: {threshold}\n")
    f.write(f"# Консенсус мотива: {motif.consensus}\n")
    f.write(f"# Всего хитов: {len(all_hits)}\n")
    f.write(f"#\n")
    f.write(f"# {'Позиция':>10} {'Цепь':>5} {'Скор':>8}  {'Последовательность'}\n")
    f.write(f"# {'-'*10} {'-'*5} {'-'*8}  {'-'*15}\n")
    for pos, strand, score, site in all_hits:
        f.write(f"  {pos:>10} {strand:>5} {score:>8.3f}  {site}\n")

print("Результаты сохранены в chr1_hits.txt")
print(f"\nПервые 10 хитов:")
for pos, strand, score, site in all_hits[:10]:
    print(f"  Позиция {pos:>8} [{strand}] скор={score:.3f}  {site}")

#!/usr/bin/env python3
import numpy as np

def nw_linear(seq1, seq2, match=3, mismatch=-3, gap=-4):
    m, n = len(seq1), len(seq2)
    # Матрица счёта
    S = np.zeros((n+1, m+1), dtype=int)
    
    # Инициализация
    for i in range(1, n+1):
        S[i, 0] = S[i-1, 0] + gap
    for j in range(1, m+1):
        S[0, j] = S[0, j-1] + gap
    
    # Заполнение
    for i in range(1, n+1):
        for j in range(1, m+1):
            diag = S[i-1, j-1] + (match if seq1[j-1] == seq2[i-1] else mismatch)
            up = S[i-1, j] + gap
            left = S[i, j-1] + gap
            S[i, j] = max(diag, up, left)
    
    # Обратный ход
    al1, al2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            cur_match = match if seq1[j-1] == seq2[i-1] else mismatch
            if S[i, j] == S[i-1, j-1] + cur_match:
                al1.append(seq1[j-1])
                al2.append(seq2[i-1])
                i -= 1; j -= 1
                continue
        if i > 0 and S[i, j] == S[i-1, j] + gap:
            al1.append('-')
            al2.append(seq2[i-1])
            i -= 1
        elif j > 0:
            al1.append(seq1[j-1])
            al2.append('-')
            j -= 1
    
    return S, ''.join(reversed(al1)), ''.join(reversed(al2))


def nw_affine(seq1, seq2, match=3, mismatch=-3, gap_open=-10, gap_extend=-1):
    m, n = len(seq1), len(seq2)
    INF = -10**9
    
    # M - совпадение/замена, X - гэп в seq1 (по горизонтали), Y - гэп в seq2 (по вертикали)
    M = np.full((n+1, m+1), INF, dtype=int)
    X = np.full((n+1, m+1), INF, dtype=int)
    Y = np.full((n+1, m+1), INF, dtype=int)
    
    M[0, 0] = 0
    
    # Инициализация краёв
    for i in range(1, n+1):
        Y[i, 0] = gap_open + (i-1) * gap_extend
    for j in range(1, m+1):
        X[0, j] = gap_open + (j-1) * gap_extend
    
    # Заполнение
    for i in range(1, n+1):
        for j in range(1, m+1):
            cur = match if seq1[j-1] == seq2[i-1] else mismatch
            M[i, j] = max(M[i-1, j-1], X[i-1, j-1], Y[i-1, j-1]) + cur
            X[i, j] = max(M[i, j-1] + gap_open, X[i, j-1] + gap_extend)
            Y[i, j] = max(M[i-1, j] + gap_open, Y[i-1, j] + gap_extend)
    
    # Обратный ход
    # Определяем, в какой матрице заканчиваем
    best = max(M[n, m], X[n, m], Y[n, m])
    al1, al2 = [], []
    i, j = n, m
    
    # Отслеживаем текущую матрицу
    if best == M[n, m]:
        state = 'M'
    elif best == X[n, m]:
        state = 'X'
    else:
        state = 'Y'
    
    while i > 0 or j > 0:
        if state == 'M':
            cur = match if seq1[j-1] == seq2[i-1] else mismatch
            al1.append(seq1[j-1])
            al2.append(seq2[i-1])
            # Куда идти?
            if M[i, j] == M[i-1, j-1] + cur:
                state = 'M'
            elif M[i, j] == X[i-1, j-1] + cur:
                state = 'X'
            else:
                state = 'Y'
            i -= 1; j -= 1
            
        elif state == 'X':  # Гэп в seq2 (seq1 имеет символ, seq2 — пропуск)
            al1.append(seq1[j-1])
            al2.append('-')
            if X[i, j] == M[i, j-1] + gap_open:
                state = 'M'
            else:
                state = 'X'  # остаёмся в X
            j -= 1
            
        else:  # state == 'Y', гэп в seq1
            al1.append('-')
            al2.append(seq2[i-1])
            if Y[i, j] == M[i-1, j] + gap_open:
                state = 'M'
            else:
                state = 'Y'
            i -= 1
    
    return best, ''.join(reversed(al1)), ''.join(reversed(al2))


if __name__ == "__main__":
    seq1 = "ATGCAGCAGCAGCCA"
    seq2 = "ATATAT"
    
    print("=" * 60)
    print("ЛИНЕЙНЫЙ ШТРАФ (gap = -4)")
    print("=" * 60)
    S, al1, al2 = nw_linear(seq1, seq2)
    print(f"Score: {S[len(seq2), len(seq1)]}")
    print("Матрица:")
    print(S)
    print(f"\nВыравнивание:")
    print(al1)
    print(al2)
    
    print("\n" + "=" * 60)
    print("АФФИННЫЙ ШТРАФ (open = -10, extend = -1)")
    print("=" * 60)
    score, al1, al2 = nw_affine(seq1, seq2)
    print(f"Score: {score}")
    print(f"\nВыравнивание:")
    print(al1)
    print(al2)

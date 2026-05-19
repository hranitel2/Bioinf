#!/usr/bin/env python3

# Параметры
MATCH = 1
MISMATCH = -1
GAP = -1
K = 4
X = 2

# Входные данные
REF = "CTAGGATCCAGGCATACGA"
QUERY = "GGATCCATTCATTA"


def build_index(database, k):
    index = {}
    for i in range(len(database) - k + 1):
        kmer = database[i:i+k]
        index.setdefault(kmer, []).append(i)
    return index


def sw_score(a, b):
    return MATCH if a == b else MISMATCH


def extend_left(query, ref, q_start, r_start, smax, scur):
    best_q, best_r = q_start, r_start
    q, r = q_start - 1, r_start - 1
    
    while q >= 0 and r >= 0:
        # Вычисляем варианты (как в SW)
        diag = scur + sw_score(query[q], ref[r])
        up = scur + GAP    # гэп в референсе
        left = scur + GAP  # гэп в запросе
        
        # Для seed-and-extend обычно используют только диагональ в простом варианте,
        # но для полноты возьмём max из трёх
        scur = max(diag, up, left, 0)
        
        if scur > smax:
            smax = scur
            best_q, best_r = q, r
        
        if smax - scur >= X:
            break
        
        q -= 1
        r -= 1
    
    return best_q, best_r, smax


def extend_right(query, ref, q_end, r_end, smax, scur):
    best_q, best_r = q_end, r_end
    q, r = q_end + 1, r_end + 1
    
    while q < len(query) and r < len(ref):
        diag = scur + sw_score(query[q], ref[r])
        up = scur + GAP
        left = scur + GAP
        scur = max(diag, up, left, 0)
        
        if scur > smax:
            smax = scur
            best_q, best_r = q, r
        
        if smax - scur >= X:
            break
        
        q += 1
        r += 1
    
    return best_q, best_r, smax


def main():
    print("=" * 70)
    print("SEED-AND-EXTEND АЛГОРИТМ")
    print("=" * 70)
    print(f"\nРеференс (D): {REF}  (len={len(REF)})")
    print(f"Запрос (Q):    {QUERY}  (len={len(QUERY)})")
    print(f"\nПараметры: k={K}, match={MATCH}, mismatch={MISMATCH}, gap={GAP}, X={X}")
    
    # Фаза 1: Индекс БД
    index = build_index(REF, K)
    print(f"\n{'='*70}")
    print("ФАЗА 1: ИНДЕКС БД (k-меры длины 4)")
    print("=" * 70)
    for kmer, positions in sorted(index.items()):
        print(f"  {kmer}: {positions}")
    
    # Разбиваем запрос на k-меры
    query_kmers = {}
    for i in range(len(QUERY) - K + 1):
        kmer = QUERY[i:i+K]
        query_kmers.setdefault(kmer, []).append(i)
    
    print(f"\n{'='*70}")
    print("ФАЗА 1: SEEDING (k-меры запроса)")
    print("=" * 70)
    
    seeds = []
    for kmer, q_positions in query_kmers.items():
        if kmer in index:
            ref_positions = index[kmer]
            print(f"  K-мер '{kmer}': в запросе на позициях {q_positions}, в БД на {ref_positions}")
            for q_pos in q_positions:
                for r_pos in ref_positions:
                    seeds.append((kmer, q_pos, r_pos))
        else:
            print(f"  K-мер '{kmer}': НЕ НАЙДЕН в БД (позиции в запросе: {q_positions})")
    
    print(f"\n  Всего найдено seed'ов: {len(seeds)}")
    
    # Фаза 2: Extension для каждого seed
    print(f"\n{'='*70}")
    print("ФАЗА 2: EXTENSION")
    print("=" * 70)
    
    best_overall = {"smax": -999, "seed": None, "q_start": 0, "q_end": 0, "r_start": 0, "r_end": 0}
    all_results = []
    
    for kmer, q_start, r_start in seeds:
        q_end = q_start + K - 1
        r_end = r_start + K - 1
        
        # Начальный счёт seed (K совпадений)
        seed_score = K * MATCH
        smax = seed_score
        scur = seed_score
        
        print(f"\n  --- Seed: '{kmer}' (Q[{q_start}:{q_end}], D[{r_start}:{r_end}]) ---")
        print(f"  Начальный счёт seed: {seed_score}")
        
        # Расширение влево
        left_q, left_r, smax_after_left = extend_left(
            QUERY, REF, q_start, r_start, smax, scur
        )
        print(f"  После расширения ВЛЕВО: Q[{left_q}], D[{left_r}], Smax={smax_after_left}")
        
        # Расширение вправо от правой границы seed
        scur = smax_after_left
        smax = smax_after_left
        right_q, right_r, final_smax = extend_right(
            QUERY, REF, q_end, r_end, smax, scur
        )
        print(f"  После расширения ВПРАВО: Q[{right_q}], D[{right_r}], Smax={final_smax}")
        
        # Выравнивание с учётом гэпов (простое диагональное)
        q_aligned = QUERY[left_q:right_q+1]
        r_aligned = REF[left_r:right_r+1]
        
        # Выравниваем участки
        al_q, al_r = "", ""
        qi, ri = left_q, left_r
        while qi <= right_q or ri <= right_r:
            if qi <= right_q and ri <= right_r:
                al_q += QUERY[qi]
                al_r += REF[ri]
                qi += 1
                ri += 1
            elif qi <= right_q:
                al_q += QUERY[qi]
                al_r += "-"
                qi += 1
            else:
                al_q += "-"
                al_r += REF[ri]
                ri += 1
        
        result = {
            "seed": kmer,
            "q_start": left_q,
            "q_end": right_q,
            "r_start": left_r,
            "r_end": right_r,
            "smax": final_smax,
            "al_q": al_q,
            "al_r": al_r
        }
        all_results.append(result)
        
        if final_smax > best_overall["smax"]:
            best_overall = result
    
    # Итоги
    print(f"\n{'='*70}")
    print("РЕЗУЛЬТАТЫ")
    print("=" * 70)
    
    print("\nВсе найденные выравнивания:")
    for i, r in enumerate(all_results, 1):
        print(f"\n  {i}. Seed: '{r['seed']}' | Smax: {r['smax']}")
        print(f"     Q[{r['q_start']}:{r['q_end']}] : {r['al_q']}")
        print(f"     D[{r['r_start']}:{r['r_end']}] : {r['al_r']}")
    
    print(f"\n  >>> ЛУЧШЕЕ ВЫРАВНИВАНИЕ <<<")
    print(f"  K-мер: '{best_overall['seed']}'")
    print(f"  Smax: {best_overall['smax']}")
    print(f"  Запрос:  {best_overall['al_q']}")
    print(f"  База:    {best_overall['al_r']}")


if __name__ == "__main__":
    main()

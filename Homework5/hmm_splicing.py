#!/usr/bin/env python3
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

STATES = ['E', '5', 'I']
N = len(STATES)
state_to_idx = {s: i for i, s in enumerate(STATES)}

start = np.log(np.array([1.0, 0.0, 0.0]))

A = np.array([
    [0.9, 0.1, 0.0],   # E -> E, 5, I
    [0.0, 0.0, 1.0],   # 5 -> I
    [0.0, 0.0, 0.9]    # I -> I
])
A_log = np.where(A > 0, np.log(A), -np.inf)

NUC = ['A', 'C', 'G', 'T']
nuc_to_idx = {n: i for i, n in enumerate(NUC)}

B = np.array([
    [0.25, 0.25, 0.25, 0.25],  # E
    [0.05, 0.0,  0.95, 0.0 ],  # 5
    [0.4,  0.1,  0.1,  0.4 ]   # I
])
B = np.where(B == 0, 1e-300, B)
B_log = np.log(B)

def forward(obs):
    T = len(obs)
    alpha = np.full((T, N), -np.inf)
    alpha[0] = start + B_log[:, nuc_to_idx[obs[0]]]
    for t in range(1, T):
        for j in range(N):
            log_sum = -np.inf
            for i in range(N):
                if A_log[i, j] != -np.inf:
                    log_sum = np.logaddexp(log_sum, alpha[t-1, i] + A_log[i, j])
            alpha[t, j] = log_sum + B_log[j, nuc_to_idx[obs[t]]]
    end_log = np.array([-np.inf, -np.inf, np.log(0.1)])  # сток только из I
    log_prob = -np.inf
    for i in range(N):
        log_prob = np.logaddexp(log_prob, alpha[T-1, i] + end_log[i])
    return log_prob, alpha


def backward(obs):
    T = len(obs)
    beta = np.full((T, N), -np.inf)
    end_log = np.array([-np.inf, -np.inf, np.log(0.1)])
    beta[T-1] = end_log
    for t in range(T-2, -1, -1):
        for i in range(N):
            log_sum = -np.inf
            for j in range(N):
                if A_log[i, j] != -np.inf and beta[t+1, j] != -np.inf:
                    log_sum = np.logaddexp(log_sum, 
                                           A_log[i, j] + B_log[j, nuc_to_idx[obs[t+1]]] + beta[t+1, j])
            beta[t, i] = log_sum
    return beta


def viterbi(obs):
    T = len(obs)
    delta = np.full((T, N), -np.inf)
    psi = np.zeros((T, N), dtype=int)
    delta[0] = start + B_log[:, nuc_to_idx[obs[0]]]
    for t in range(1, T):
        for j in range(N):
            best_val = -np.inf
            best_i = -1
            for i in range(N):
                if A_log[i, j] != -np.inf:
                    val = delta[t-1, i] + A_log[i, j]
                    if val > best_val:
                        best_val = val
                        best_i = i
            delta[t, j] = best_val + B_log[j, nuc_to_idx[obs[t]]]
            psi[t, j] = best_i
    end_log = np.array([-np.inf, -np.inf, np.log(0.1)])
    final_best = -np.inf
    final_state = -1
    for i in range(N):
        val = delta[T-1, i] + end_log[i]
        if val > final_best:
            final_best = val
            final_state = i
    path = [final_state]
    for t in range(T-1, 0, -1):
        path.append(psi[t, path[-1]])
    path.reverse()
    return ''.join([STATES[i] for i in path]), final_best


def posterior(alpha, beta):
    """Возвращает gamma[t][i]."""
    T = alpha.shape[0]
    gamma = np.zeros((T, N))
    for t in range(T):
        log_posterior = alpha[t] + beta[t]
        log_norm = -np.inf
        for val in log_posterior:
            log_norm = np.logaddexp(log_norm, val)
        gamma[t] = np.exp(log_posterior - log_norm)
    return gamma


obs_sequence = "CTTCATGTGAAAGCAGACGTAAGTCA"
obs_sequence = obs_sequence.upper()
print("Последовательность:", obs_sequence)
    
vit_path, vit_logprob = viterbi(obs_sequence)
print("\n=== Алгоритм Витерби ===")
print("Оптимальный путь состояний:", vit_path)
print(f"Логарифм вероятности: {vit_logprob:.4f}")
print(f"Вероятность: {np.exp(vit_logprob):.6e}")
    
log_prob_fwd, alpha = forward(obs_sequence)
print("\n=== Алгоритм Forward ===")
print(f"Log P(O) = {log_prob_fwd:.4f}")
print(f"P(O) = {np.exp(log_prob_fwd):.6e}")
    
beta = backward(obs_sequence)
gamma = posterior(alpha, beta)
    
print("\n=== Апостериорные максимумы ===")
for s_idx, s_name in enumerate(STATES):
    max_t = np.argmax(gamma[:, s_idx])
    max_val = gamma[max_t, s_idx]
    print(f"Состояние {s_name}: максимум {max_val:.4f} на позиции {max_t+1}")
    
T = len(obs_sequence)
t_axis = np.arange(1, T+1)
fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
colors = ['blue', 'orange', 'green']
for idx, (ax, s_name) in enumerate(zip(axes, STATES)):
    ax.plot(t_axis, gamma[:, idx], color=colors[idx], linewidth=1.5)
    ax.set_ylabel(f'P({s_name})', fontsize=12)
    ax.set_title(f'Апостериорная вероятность состояния {s_name}')
    ax.grid(True, alpha=0.3)
    max_t = np.argmax(gamma[:, idx])
    max_val = gamma[max_t, idx]
    ax.plot(max_t+1, max_val, 'ro', markersize=5)
    ax.annotate(f'Макс: {max_val:.3f} (поз. {max_t+1})', 
                 xy=(max_t+1, max_val), xytext=(max_t+1+2, max_val+0.1),
                     arrowprops=dict(arrowstyle='->', color='red'))
    ax.set_xlabel('Позиция t')
plt.tight_layout()
plt.savefig("posterior_gamma.png", dpi=150)
print("\nГрафики сохранены в posterior_gamma.png")

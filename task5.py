import numpy as np


def true_signal(t):
    return np.exp(t) * np.sin(np.pi * t)


def L_matrix(n):
    L = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                L[i, j] = 1
            elif i == j + 1:
                L[i, j] = -1
    return L


print(L_matrix(5))

import numpy as np
import matplotlib.pyplot as plt

from gsvd import gsvd


def blur_matrix(t, h, beta):
    c = np.sqrt(beta / np.pi)
    return h * c * np.exp(-beta * (t[:, None] - t[None, :]) ** 2)


def first_difference(n):
    return np.eye(n - 1, n) - np.eye(n - 1, n, k=1)


def gsvd_parts(A, L):
    n = A.shape[1]
    U, V, C, S, W, Winv = gsvd(A, L)
    alpha = np.diag(C)
    beta = np.append(np.diag(S), 0)
    return U[:, :n], W, alpha, beta


def regularized_solution(U, W, alpha, beta, b, lam):
    y = alpha / (alpha**2 + lam**2 * beta**2) * (U.T @ b)
    return W @ y


def true_signal(t):
    return np.exp(t) * np.sin(np.pi * t)


def rel_err(x, x_true):
    return np.linalg.norm(x - x_true) / np.linalg.norm(x_true)


# Parameters
n = 100  # number of sample points
beta_blur = 10.0  # kernel parameter; larger beta = narrower blur
delta = 1e-2  # noise level
lam = 1e-1  # regularization parameter

# Sample the true signal at n box midpoints in [-1, 1]
h = 2 / n
t = -1 + h * (np.arange(n) + 0.5)
x_true = true_signal(t)

# Blurred, noisy data b = A x + e with e = delta * N(0, 1)
A = blur_matrix(t, h, beta_blur)
rng = np.random.default_rng(0)
b = A @ x_true + delta * rng.standard_normal(n)

# Unregularized least-squares solution
x_ls = np.linalg.lstsq(A, b, rcond=None)[0]

# Regularized solution via the GSVD of (A, L)
U, W, alpha, beta = gsvd_parts(A, first_difference(n))
x_reg = regularized_solution(U, W, alpha, beta, b, lam)

print(f"Relative error, data b        = {rel_err(b, x_true):.3f}")
print(f"Relative error, unregularized = {rel_err(x_ls, x_true):.1e}")
print(f"Relative error, regularized   = {rel_err(x_reg, x_true):.3f}")

# Plot. The unregularized solution is huge, so it gets its own panel.
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
ax1.plot(t, x_true, "k-", lw=2, label="true signal $x$")
ax1.plot(t, b, ".", color="gray", label="blurred noisy data $b$")
ax1.plot(t, x_reg, "tab:blue", label=f"regularized, $\\lambda$ = {lam:g}")
ax2.plot(t, x_true, "k-", lw=2, label="true signal $x$")
ax2.plot(t, x_ls, "tab:red", label="unregularized least squares")
for ax in (ax1, ax2):
    ax.set_xlabel("$t$")
    ax.legend()
fig.tight_layout()
fig.savefig("task5_1.png", dpi=150)
plt.show()

import numpy as np
import matplotlib.pyplot as plt

from gsvd import gsvd

# ---------------------------------------------------------------------------
#  Problem setup (Method, section "One-dimensional deblurring")
# ---------------------------------------------------------------------------


def true_signal(t):
    """The signal we try to recover, x(t) = e^t sin(pi t)."""
    return np.exp(t) * np.sin(np.pi * t)


def blur_matrix(t, h, beta):
    """Discretize the blur integral with the box rule: A_ij = h K(t_i, t_j).

    K(s, t) = c exp(-beta (s - t)^2) is the Gaussian kernel. We choose
    c = sqrt(beta / pi) so that the kernel integrates to 1; the blur then
    smears the signal out without changing its overall size.
    """
    c = np.sqrt(beta / np.pi)
    s, tt = np.meshgrid(t, t, indexing="ij")  # s = t_i (rows), tt = t_j (columns)
    return h * c * np.exp(-beta * (s - tt) ** 2)


def first_difference(n):
    """The (n-1) x n first-difference matrix, (L x)_i = x_i - x_{i+1}."""
    return np.eye(n - 1, n) - np.eye(n - 1, n, k=1)


def relative_error(x, x_true):
    return np.linalg.norm(x - x_true) / np.linalg.norm(x_true)


# ---------------------------------------------------------------------------
#  Regularized solution via the GSVD (Method, general-form regularization)
# ---------------------------------------------------------------------------


def gsvd_factors(A, L):
    """Compute the GSVD once and return U, W and the diagonals alpha, beta.

    alpha_i and beta_i are the diagonals of Sigma_A and Sigma_L. Since L has
    only n-1 rows, beta_n = 0: that component (a constant signal) lies in
    the null space of L and is not regularized.
    """
    U, V, C, S, W, Winv = gsvd(A, L)
    n = A.shape[1]
    alpha = np.diag(C)[:n]
    beta = np.zeros(n)
    beta[: min(S.shape)] = np.diag(S)
    return U[:, :n], W, alpha, beta


def filter_factors(alpha, beta, lam):
    """phi_i(lambda) = alpha_i^2 / (alpha_i^2 + lambda^2 beta_i^2)."""
    return alpha**2 / (alpha**2 + lam**2 * beta**2)


def regularized_solution(U, W, alpha, beta, b, lam):
    """x = W y with y_i = alpha_i / (alpha_i^2 + lambda^2 beta_i^2) (u_i^T b).

    This equals phi_i (u_i^T b) / alpha_i, but written this way we never
    divide by alpha_i. That matters here: A is so ill-conditioned that some
    alpha_i are numerically zero.
    """
    y = alpha / (alpha**2 + lam**2 * beta**2) * (U.T @ b)
    return W @ y


# ---------------------------------------------------------------------------
#  Experiment 1: regularized vs. unregularized solution
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    # Parameters
    n = 100  # number of sample points
    beta_blur = 10.0  # kernel parameter; larger beta = narrower blur
    delta = 1e-2  # noise level
    lam = 1e-1  # regularization parameter (chosen by hand for now)

    # Sample the true signal at n equidistant points (box midpoints) in [-1, 1]
    a, b_end = -1.0, 1.0
    h = (b_end - a) / n
    t = a + h * (np.arange(n) + 0.5)
    x_true = true_signal(t)

    # Blurred, noisy data b = A x + e with e = delta * N(0, 1)
    A = blur_matrix(t, h, beta_blur)
    rng = np.random.default_rng(0)  # fixed seed so the results are reproducible
    b = A @ x_true + delta * rng.standard_normal(n)

    # Unregularized least-squares solution, min ||A x - b||
    x_ls = np.linalg.lstsq(A, b, rcond=None)[0]

    # Regularized solution via the GSVD of (A, L)
    L = first_difference(n)
    U, W, alpha, beta = gsvd_factors(A, L)
    x_reg = regularized_solution(U, W, alpha, beta, b, lam)

    # Sanity check: the GSVD solution must satisfy the normal equations (Eq. 3)
    x_normal = np.linalg.solve(A.T @ A + lam**2 * L.T @ L, A.T @ b)
    print(f"cond(A)                       = {np.linalg.cond(A):.1e}")
    print(f"GSVD vs normal equations      = {relative_error(x_reg, x_normal):.1e}")
    print(f"Relative error, data b        = {relative_error(b, x_true):.3f}")
    print(f"Relative error, unregularized = {relative_error(x_ls, x_true):.1e}")
    print(f"Relative error, regularized   = {relative_error(x_reg, x_true):.3f}")

    # Plot. The unregularized solution is far too large to share an axis with
    # the others, so it gets its own panel.
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

    ax1.plot(t, x_true, "k-", lw=2, label="true signal $x$")
    ax1.plot(t, b, ".", color="gray", ms=4, label="blurred noisy data $b$")
    ax1.plot(t, x_reg, "tab:blue", lw=1.5, label=f"regularized, $\\lambda$ = {lam:g}")
    ax1.set_title("Regularized reconstruction")

    ax2.plot(t, x_true, "k-", lw=2, label="true signal $x$")
    ax2.plot(t, x_ls, "tab:red", lw=1, label="unregularized least squares")
    ax2.set_title("Unregularized reconstruction")

    for ax in (ax1, ax2):
        ax.set_xlabel("$t$")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9)

    fig.tight_layout()
    fig.savefig("task5_1.png", dpi=150)
    plt.show()

"""Task 5 – 1D-deblurring med Tikhonov-regularisering via GSVD (deluppgift 5.1–5.5)."""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from gsvd import gsvd


# ---------- Problemuppsättning ----------
def true_signal(t):
    return np.exp(t) * np.sin(np.pi * t)


def box_quadrature_nodes(a, b, n):
    """Noder och vikter för boxregeln (mittpunktsregeln) på [a,b] med n lådor.
    Returnerar t (mittpunkter) och w (vikter, alla lika med h)."""
    h = (b - a) / n
    t = a + h * (np.arange(n) + 0.5)
    w = np.full(n, h)
    return t, w


def box_quadrature_rule(f, a, b, n):
    """Beräkna integral av f på [a,b] med n mittpunkter (boxregel).
    f måste kunna ta en numpy-array."""
    t, w = box_quadrature_nodes(a, b, n)
    return np.sum(w * f(t))


def blur_matrix(c, beta, a, b, n):
    """Diskretisera (Kx)(s) = ∫_a^b K(s,t) x(t) dt med boxregeln,
    K(s,t) = c exp(-beta (s-t)^2), ekv. (1.9).
    Rad i är boxregeln för integralen i s = t_i:
        (Ax)_i = sum_j w_j K(t_i, t_j) x_j,  dvs  A_ij = w_j K(t_i, t_j).
    Returnerar A och noderna t."""
    t, w = box_quadrature_nodes(a, b, n)
    Ti, Tj = np.meshgrid(t, t, indexing="ij")
    K = c * np.exp(-beta * (Ti - Tj) ** 2)
    return K * w[None, :], t


def first_difference(n):
    """(n-1) x n med 1 på diagonalen och -1 ovanför."""
    return np.eye(n - 1, n) - np.eye(n - 1, n, k=1)


def noisy_data(A, x, delta, seed=0):
    rng = np.random.default_rng(seed)
    return A @ x + delta * rng.standard_normal(len(x))


def rel_err(x, x_true):
    return np.linalg.norm(x - x_true) / np.linalg.norm(x_true)


# ---------- GSVD-verktyg ----------
def gsvd_parts(A, L):
    """Kör GSVD en gång och returnera U, W samt alpha_i, beta_i (längd n).
    beta_n = 0 motsvarar nollrummet av L (konstanta vektorer)."""
    U, V, C, S, W, Winv = gsvd(A, L)
    n = A.shape[1]
    alpha = np.diag(C)[:n]
    beta = np.zeros(n)
    k = min(S.shape)
    beta[:k] = np.diag(S)[:k]
    return U[:, :n], W, alpha, beta


def filter_factors(alpha, beta, lam):
    """phi_i = alpha_i^2 / (alpha_i^2 + lam^2 beta_i^2)."""
    return alpha**2 / (alpha**2 + lam**2 * beta**2)


def tikhonov_gsvd(U, W, alpha, beta, b, lam):
    """Formel (1.7): y_i = alpha_i / (alpha_i^2 + lam^2 beta_i^2) * u_i^T b, x = W y.
    Delar aldrig med små alpha_i (till skillnad från phi_i / alpha_i)."""
    y = filter_factors(alpha, beta, lam) / alpha * (U.T @ b)
    return W @ y


def discrepancy_lambda(A, b, x_lam, delta, tau=1.0, lo=-6, hi=2):
    """Hittar lambda så att ||A x_lambda - b|| = tau * delta * sqrt(m).
    Residualen räknas direkt (inte via U^T b) eftersom kursens U inte är
    exakt ortogonal för så illa konditionerade A. Söker i log10(lambda)."""
    target = tau * delta * np.sqrt(len(b))
    f = lambda s: np.linalg.norm(A @ x_lam(10**s) - b) - target
    return 10 ** brentq(f, lo, hi)


# ---------- Experiment ----------
if __name__ == "__main__":
    n = 100
    a, b_end = -1.0, 1.0
    h = (b_end - a) / n
    t = a + h * (np.arange(n) + 0.5)  # mittpunkter

    beta_k = 100.0  # kärnparameter (B = beta I)
    # c = np.sqrt(beta_k / np.pi)  # normering: kärnan integrerar till 1
    c = 1
    delta = 1e-2  # basbrusnivå
    lam_manual = 1e-1  # handvalt lambda (se 5.3)

    x_true = true_signal(t)
    A = blur_matrix(c, t, h, beta_k)
    L = first_difference(n)
    b = noisy_data(A, x_true, delta)

    # GSVD beror bara på (A, L) -> beräknas en gång
    U, W, alpha, beta = gsvd_parts(A, L)
    solve = lambda b_, lam_: tikhonov_gsvd(U, W, alpha, beta, b_, lam_)

    print(f"cond(A) = {np.linalg.cond(A):.2e}")
    U_, V_, C_, S_, W_, Winv_ = gsvd(A, L)
    print("A = U C W^-1 ?", np.allclose(U_ @ C_ @ Winv_, A))
    print("L = V S W^-1 ?", np.allclose(V_ @ S_ @ Winv_, L))

    # =====================================================================
    # 5.1  Oregulariserad vs regulariserad lösning
    # =====================================================================
    print("\n--- 5.1 ---")
    x_ls = np.linalg.solve(A, b)
    x_reg = solve(b, lam_manual)
    x_ref = np.linalg.solve(A.T @ A + lam_manual**2 * L.T @ L, A.T @ b)

    print(f"GSVD vs normalekvationer  = {rel_err(x_reg, x_ref):.2e}")
    print(f"Rel. fel, data b          = {rel_err(b, x_true):.3f}")
    print(f"Rel. fel, oregulariserad  = {rel_err(x_ls, x_true):.3e}")
    print(f"Rel. fel, regulariserad   = {rel_err(x_reg, x_true):.3f}")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    ax = axes[0]
    ax.plot(t, x_true, "k-", lw=2, label="Sann signal $x$")
    ax.plot(t, b, ".", color="tab:gray", ms=4, label="Data $b = Ax + e$")
    ax.plot(
        t,
        x_reg,
        "tab:blue",
        lw=1.5,
        label=f"Tikhonov/GSVD, $\\lambda$={lam_manual:g} (fel {rel_err(x_reg, x_true):.3f})",
    )
    ax.set_title(
        f"Regulariserad rekonstruktion ($\\beta$={beta_k:g}, $\\delta$={delta:g})"
    )
    ax.set_xlabel("$t$")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(t, x_true, "k-", lw=2, label="Sann signal $x$")
    ax.plot(
        t,
        x_ls,
        "tab:red",
        lw=1,
        label=f"Oregulariserad (fel {rel_err(x_ls, x_true):.1e})",
    )
    ax.set_title("Oregulariserad lösning $x = A^{-1}b$")
    ax.set_xlabel("$t$")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig("task5_1.png", dpi=150)

    # =====================================================================
    # 5.2  Filterfaktorer
    # =====================================================================
    print("\n--- 5.2 ---")
    # Generaliserade singulärvärden gamma_i = alpha_i / beta_i (beta_n = 0 -> inf, exkluderas)
    mask = beta > 0
    gamma = alpha[mask] / beta[mask]
    order = np.argsort(gamma)
    print(f"gamma_i spänner [{gamma.min():.1e}, {gamma.max():.1e}]")
    print(
        f"Komponent med beta=0: alpha={alpha[~mask][0]:.3f}, "
        f"phi={filter_factors(alpha, beta, lam_manual)[~mask][0]:.3f} (regulariseras ej)"
    )

    lams_ff = [1e-3, 1e-2, 1e-1, 1, 10]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    ax = axes[0]
    for lam_ in lams_ff:
        phi = filter_factors(alpha, beta, lam_)
        ax.semilogy(
            np.arange(1, n + 1), phi, ".-", ms=3, lw=1, label=f"$\\lambda$={lam_:g}"
        )
    ax.set_title("Filterfaktorer $\\phi_i(\\lambda)$ per index")
    ax.set_xlabel("Index $i$ (stigande $\\alpha_i$)")
    ax.set_ylabel("$\\phi_i$")
    ax.set_ylim(1e-20, 2)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)

    ax = axes[1]
    for lam_ in lams_ff:
        phi = filter_factors(alpha[mask], beta[mask], lam_)
        (line,) = ax.loglog(
            gamma[order], phi[order], "-", lw=1.5, label=f"$\\lambda$={lam_:g}"
        )
        ax.axvline(lam_, color=line.get_color(), ls=":", lw=1)
    ax.set_title(
        "$\\phi$ mot $\\gamma_i=\\alpha_i/\\beta_i$ (streckat: $\\gamma=\\lambda$)"
    )
    ax.set_xlabel("$\\gamma_i$")
    ax.set_ylabel("$\\phi_i$")
    ax.set_ylim(1e-20, 2)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("task5_2.png", dpi=150)

    # =====================================================================
    # 5.3  Olika lambda: under- och överregularisering
    # =====================================================================
    print("\n--- 5.3 ---")
    lams_show = [1e-3, 1e-2, 1e-1, 1, 10]
    lam_grid = np.logspace(-4, 2, 121)
    errs = np.array([rel_err(solve(b, l_), x_true) for l_ in lam_grid])
    i_best = np.argmin(errs)
    lam_opt = lam_grid[i_best]
    for lam_ in lams_show:
        print(f"lambda={lam_:<6g} rel. fel = {rel_err(solve(b, lam_), x_true):.3f}")
    print(f"Optimalt (mot x_true): lambda={lam_opt:.3g}, fel={errs[i_best]:.3f}")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    ax = axes[0]
    ax.plot(t, x_true, "k-", lw=2.5, label="Sann signal $x$")
    for lam_ in lams_show:
        x_ = solve(b, lam_)
        ax.plot(
            t, x_, lw=1.3, label=f"$\\lambda$={lam_:g} (fel {rel_err(x_, x_true):.3f})"
        )
    ax.set_ylim(1.5 * x_true.min(), 1.5 * x_true.max())
    ax.set_title(f"Regulariserade lösningar, $\\delta$={delta:g}")
    ax.set_xlabel("$t$")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.loglog(lam_grid, errs, "k-")
    ax.plot(
        lam_opt,
        errs[i_best],
        "o",
        color="tab:green",
        label=f"minimum $\\lambda$={lam_opt:.2g}",
    )
    ax.plot(
        lam_manual,
        rel_err(solve(b, lam_manual), x_true),
        "s",
        color="tab:blue",
        label=f"handvalt $\\lambda$={lam_manual:g}",
    )
    ax.text(2e-4, errs[0] * 0.6, "underregularisering\n(brus dominerar)", fontsize=9)
    ax.text(
        5, errs[-1] * 0.5, "över-\nregularisering\n(utslätning)", fontsize=9, ha="right"
    )
    ax.set_title("Relativt fel mot $\\lambda$")
    ax.set_xlabel("$\\lambda$")
    ax.set_ylabel("$\\|x_\\lambda - x\\| / \\|x\\|$")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("task5_3.png", dpi=150)

    # =====================================================================
    # 5.4  Olika brusnivåer
    # =====================================================================
    print("\n--- 5.4 ---")
    deltas = [1e-3, 1e-2, 1e-1]
    lams_4 = [1e-3, 1e-2, 1e-1, 1, 10]
    fig, axes = plt.subplots(1, len(deltas) + 1, figsize=(5 * (len(deltas) + 1), 4.5))
    print(f"{'delta':>7} {'oreg. fel':>11} {'bästa lambda':>13} {'bästa fel':>10}")
    for ax, d in zip(axes, deltas):
        b_d = noisy_data(A, x_true, d)  # samma seed: samma brusform, skalad
        err_ls = rel_err(np.linalg.solve(A, b_d), x_true)
        errs_d = [rel_err(solve(b_d, l_), x_true) for l_ in lam_grid]
        j = np.argmin(errs_d)
        print(f"{d:>7g} {err_ls:>11.2e} {lam_grid[j]:>13.3g} {errs_d[j]:>10.3f}")

        ax.plot(t, x_true, "k-", lw=2.5, label="Sann signal $x$")
        for lam_ in lams_4:
            x_ = solve(b_d, lam_)
            ax.plot(
                t,
                x_,
                lw=1.2,
                label=f"$\\lambda$={lam_:g} (fel {rel_err(x_, x_true):.3f})",
            )
        ax.set_title(f"$\\delta$={d:g}   (oregulariserad: fel {err_ls:.1e})")
        ax.set_xlabel("$t$")
        ax.set_ylim(1.5 * x_true.min(), 1.5 * x_true.max())
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)

    ax = axes[-1]
    for d in deltas:
        b_d = noisy_data(A, x_true, d)
        errs_d = [rel_err(solve(b_d, l_), x_true) for l_ in lam_grid]
        j = np.argmin(errs_d)
        ax.loglog(
            lam_grid,
            errs_d,
            label=f"$\\delta$={d:g} (bäst $\\lambda$={lam_grid[j]:.2g})",
        )
        ax.plot(lam_grid[j], errs_d[j], "ko", ms=4)
    ax.set_title("Relativt fel mot $\\lambda$")
    ax.set_xlabel("$\\lambda$")
    ax.set_ylabel("Relativt fel")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)
    fig.suptitle("Tikhonov/GSVD för olika brusnivåer $\\delta$")
    fig.tight_layout()
    fig.savefig("task5_4.png", dpi=150)

    # =====================================================================
    # 5.5  Automatiskt val av lambda: diskrepansprincipen
    # =====================================================================
    print("\n--- 5.5 ---")
    # b och delta är fortfarande basfallet (delta = 1e-2) – loopvariablerna ovan heter d, b_d
    x_lam = lambda lam_: solve(b, lam_)
    lam_dp = discrepancy_lambda(A, b, x_lam, delta)
    x_dp = x_lam(lam_dp)
    x_man = x_lam(lam_manual)
    print(f"Diskrepans: lambda = {lam_dp:.3g}, fel = {rel_err(x_dp, x_true):.3f}")
    print(f"Handvalt:   lambda = {lam_manual:.3g}, fel = {rel_err(x_man, x_true):.3f}")
    print(f"Optimalt:   lambda = {lam_opt:.3g}, fel = {errs[i_best]:.3f}")

    print("\nDiskrepans per brusnivå:")
    for d in deltas:
        b_d = noisy_data(A, x_true, d)
        l_dp = discrepancy_lambda(A, b_d, lambda l_: solve(b_d, l_), d)
        errs_d = [rel_err(solve(b_d, l_), x_true) for l_ in lam_grid]
        j = np.argmin(errs_d)
        print(
            f"  delta={d:<6g} lambda_DP={l_dp:.3g} (fel {rel_err(solve(b_d, l_dp), x_true):.3f})"
            f"   optimalt lambda={lam_grid[j]:.3g} (fel {errs_d[j]:.3f})"
        )

    res = np.array([np.linalg.norm(A @ x_lam(l_) - b) for l_ in lam_grid])
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    ax = axes[0]
    ax.loglog(lam_grid, res, "k-", label="$\\|Ax_\\lambda - b\\|$")
    ax.axhline(delta * np.sqrt(n), color="tab:red", ls="--", label="$\\delta\\sqrt{n}$")
    ax.axvline(lam_dp, color="tab:red", ls=":", label=f"$\\lambda_{{DP}}$={lam_dp:.2g}")
    ax.axvline(
        lam_manual,
        color="tab:blue",
        ls=":",
        label=f"handvalt $\\lambda$={lam_manual:g}",
    )
    ax.axvline(
        lam_opt, color="tab:green", ls=":", label=f"optimalt $\\lambda$={lam_opt:.2g}"
    )
    ax.set_title("Diskrepansprincipen: residual mot $\\lambda$")
    ax.set_xlabel("$\\lambda$")
    ax.set_ylabel("residualnorm")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.plot(t, x_true, "k-", lw=2.5, label="Sann signal $x$")
    ax.plot(
        t,
        x_dp,
        "tab:red",
        lw=1.5,
        label=f"$\\lambda_{{DP}}$={lam_dp:.2g} (fel {rel_err(x_dp, x_true):.3f})",
    )
    ax.plot(
        t,
        x_man,
        "tab:blue",
        lw=1.5,
        ls="--",
        label=f"handvalt $\\lambda$={lam_manual:g} (fel {rel_err(x_man, x_true):.3f})",
    )
    ax.set_title(f"Rekonstruktioner, $\\delta$={delta:g}")
    ax.set_xlabel("$t$")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig("task5_5.png", dpi=150)

    plt.show()

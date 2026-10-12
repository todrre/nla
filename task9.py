import numpy as np
def minres(A, b, x0, tol=1e-6, maxiter=None):
    if maxiter is None:
        maxiter = len(b)
    x = x0.copy()
    r = b - A@x0
    rho = np.linalg.norm(r)
    bnorm = np.linalg.norm(b)
    relres = rho / bnorm
    m = 0
    v = np.zeros_like(b)
    v_new = r / rho
    beta = 0
    p, p_new = np.zeros_like(b), np.zeros_like(b)
    # Givens-rotationer G_{m-2}, G_{m-1} som (c, s): startar som identitet
    c_old, s_old = 1, 0
    c, s = 1, 0
    tau = rho  # |tau| = ||r_m||
    while relres > tol and m < maxiter:
        # Lanczos
        w = A@v_new - beta*v
        alpha = v_new@w
        w = w - alpha*v_new
        beta_new = np.linalg.norm(w)

        # tidigare rotationer på nya kolonnen (beta, alpha, beta_new) i T
        eps = s_old*beta
        delta_bar = c_old*beta
        delta = c*delta_bar + s*alpha
        gam_bar = -s*delta_bar + c*alpha

        # ny rotation som nollar beta_new
        gam = np.hypot(gam_bar, beta_new)
        c_old, s_old = c, s
        c, s = gam_bar/gam, beta_new/gam

        p, p_new = p_new, (v_new - delta*p_new - eps*p) / gam
        x = x + c*tau*p_new
        tau = -s*tau
        relres = abs(tau) / bnorm
        m += 1

        if beta_new == 0:
            break
        v = v_new
        v_new = w/beta_new
        beta = beta_new
    return x, relres, m

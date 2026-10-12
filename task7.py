import numpy as np
def cg(A, b, x_0, tol=1e-10, max_iter=None):
    r_0 = b - A @ x_0
    p_0 = r_0.copy()
    relres = 1
    m = 0
    if max_iter is None:
        max_iter = len(b)
    r = r_0.copy()
    x = x_0.copy()
    p = p_0.copy()
    while relres > tol and m < max_iter:
        Ap = A @ p
        alpha = (r @ r) / (p @ Ap)
        x = x + alpha * p
        r_new = r - alpha * Ap
        beta = (r_new @ r_new) / (r @ r)
        p = r_new + beta * p
        r = r_new
        relres = np.linalg.norm(r) / np.linalg.norm(b)
        m += 1
    return x, relres, m

A = np.array([[4, 1], [1, 3]])
b = np.array([1, 2])
x_0 = np.zeros_like(b)
x, relres, m = cg(A, b, x_0)
print(f"Solution: {x}, Relative Residual: {relres}, Iterations: {m}")
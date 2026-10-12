import numpy as np
def minres(A, b, x0, tol=1e-6, maxiter=None):
    x = x0
    r = b - A@x0
    rho = np.linalg.norm(r)
    relres = 1
    m = 1
    v = 0
    v_new =  r / rho
    beta = 0
    p, p_new = 0, 0
    gam, gam_new = 0, 0
    G, G_new = 0, 0 # claude säger att det inte ska vara riktiga nollor 
    tau = 1
    while relres > tol and m < maxiter:
        w = A@v_new - beta*v
        alpha = v@w
        w = w - alpha*v_new
        beta = np.linalg.norm(w)
        if beta == 0:
            break
        v = v_new
        v_new = w/beta


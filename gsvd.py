import numpy as np
from scipy.linalg import qr, svd, null_space


def gsvd(A, L, tol = None):
    """
    GSVD-like factorization of the matrix pair (A,L),

        A : m x n,   m >= n
        L : p x n,

    assuming the stacked matrix [A; L] has full column rank.

    Computes

        A = U @ C @ Winv
        L = V @ S @ Winv

    where

        U    : m x m       orthogonal
        V    : p x p       orthogonal
        C    : m x n       rectangular diagonal
        S    : p x n       rectangular diagonal
        W    : n x n       nonsingular
        Winv : n x n       inverse of W

    and

        C.T @ C + S.T @ S = I_n.

    Parameters
    ----------
    A : ndarray, shape (m,n)
    L : ndarray, shape (p,n)
    tol : float, optional
        Numerical tolerance.

    Returns
    -------
    U, V, C, S, W, Winv
    """

    A = np.asarray(A, dtype=float)
    L = np.asarray(L, dtype=float)

    m, n = A.shape
    p, n2 = L.shape

    if m < n:
        raise ValueError("A must satisfy m >= n.")

    if n2 != n:
        raise ValueError(
            "A and L must have the same number of columns."
        )

    # ---------------------------------------------------------
    # 1. QR factorization of the stacked matrix
    #
    #        [ A ; L ] = Q R
    #
    # where
    #
    #        Q : (m+p) x n
    #        R : n x n
    #
    # and Q has orthonormal columns.
    # ---------------------------------------------------------

    Z = np.vstack((A, L))

    Q, R = qr(Z, mode='economic')

    if np.linalg.matrix_rank(R) < n:
        raise ValueError(
            "The stacked matrix [A; L] must have full column rank."
        )

    Q1 = Q[:m, :]       # m x n
    Q2 = Q[m:, :]       # p x n

    # We have
    #
    #     Q1.T @ Q1 + Q2.T @ Q2 = I_n.

    # ---------------------------------------------------------
    # 2. Full SVD of Q2
    #
    #        Q2 = V S X^T
    #
    # where
    #
    #        V : p x p
    #        S : p x n
    #        X : n x n
    # ---------------------------------------------------------

    V, s, Xt = svd(Q2, full_matrices=True)

    X = Xt.T

    r = min(p, n)

    # ---------------------------------------------------------
    # 3. Construct the generalized sines and cosines.
    #
    # The SVD gives r = min(p,n) singular values.
    # If p < n, the remaining generalized sines are zero.
    #
    # Since
    #
    #   Q1.T Q1 + Q2.T Q2 = I,
    #
    # we have
    #
    #   (Q1 X).T (Q1 X) = I - S.T S.
    #
    # Thus
    #
    #        c_j^2 + s_j^2 = 1.
    # ---------------------------------------------------------

    s_full = np.zeros(n)
    s_full[:r] = s

    c = np.sqrt(np.maximum(0.0, 1.0 - s_full**2))

    # ---------------------------------------------------------
    # 4. Construct the first n columns of U from
    #
    #        Q1 X = U1 diag(c).
    #
    # For c_j > 0,
    #
    #        u_j = Q1 X e_j / c_j.
    #
    # If c_j = 0, the corresponding column is obtained by
    # orthogonal completion.
    # ---------------------------------------------------------

    Q1X = Q1 @ X

    if tol is None:
        tol = (
            100
            * np.finfo(float).eps
            * max(m, n, p)
        )

    nz = c > tol
    z = ~nz

    U1 = np.zeros((m, n))

    # Columns determined by Q1 X
    U1[:, nz] = Q1X[:, nz] / c[nz]

    # Complete the columns corresponding to c_j = 0
    if np.any(z):

        Uknown = U1[:, nz]

        if Uknown.shape[1] == 0:
            Ucomp = np.eye(m)
        else:
            Ucomp = null_space(Uknown.T)

        if Ucomp.shape[1] < np.sum(z):
            raise ValueError(
                "Could not construct the required columns of U."
            )

        U1[:, z] = Ucomp[:, :np.sum(z)]

    # ---------------------------------------------------------
    # 5. Complete U1 to a full m x m orthogonal matrix U.
    # ---------------------------------------------------------

    if m > n:

        Uperp = null_space(U1.T)

        if Uperp.shape[1] != m - n:
            raise ValueError(
                "Could not complete U to an orthogonal matrix."
            )

        U = np.hstack((U1, Uperp))

    else:
        U = U1

    # ---------------------------------------------------------
    # 6. Construct the rectangular matrix C
    #
    #           [ diag(c) ]
    #       C = [         ]     m x n.
    #           [    0    ]
    # ---------------------------------------------------------

    C = np.zeros((m, n))
    C[:n, :] = np.diag(c)

    # ---------------------------------------------------------
    # 7. Construct the rectangular matrix S
    #
    # S is p x n with singular values of Q2 on its diagonal.
    # ---------------------------------------------------------

    S = np.zeros((p, n))

    for j in range(r):
        S[j, j] = s[j]

    # ---------------------------------------------------------
    # 8. Construct W.
    #
    # We have
    #
    #        Q1 X = U C
    #
    # and
    #
    #        Q2 X = V S.
    #
    # Therefore
    #
    #        Q1 = U C X^T
    #        Q2 = V S X^T.
    #
    # Since
    #
    #        A = Q1 R
    #        L = Q2 R,
    #
    # we obtain
    #
    #        A = U C X^T R
    #        L = V S X^T R.
    #
    # Thus
    #
    #        Winv = W^{-1} = X^T R.
    # ---------------------------------------------------------

    Winv = Xt @ R

    W = np.linalg.solve(Winv, np.eye(n))

    return U, V, C, S, W, Winv
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt
from PIL import Image

# Importera bilden, skala om till 300x300 och gör den gråskalig
n1 = n2 = 300
n = n1*n2
img = Image.open("img.png")
img = img.resize((n1, n2))
img = img.convert("L")
beta = 1250
tau = 0.9
sigma = 0.01
x = np.array(img.convert("L")).flatten()
t = range(n1*n2)
e = np.random.normal(0, 1, n1*n2)

# Konstruera en blurring-matris
def blur1d(n, beta):
    t = np.linspace(0, 1, n)
    c = np.sqrt(beta / np.pi)
    h = t[1] - t[0]
    return h * c * np.exp(-beta* (t[:, None] - t[None, :]) ** 2)

def build_blur_matrix(n1, n2, beta, tau):
    A1, A2 = blur1d(n1, beta), blur1d(n2, beta)
    eta = tau * A1.max() * A2.max()          # relativ -> absolut tröskel
    A1s = sp.csr_matrix(np.where(A1 * A2.max() >= eta, A1, 0))
    A2s = sp.csr_matrix(np.where(A2 * A1.max() >= eta, A2, 0))
    A = sp.kron(A2s, A1s, format="csr")      # matchar order="F"
    A.data[np.abs(A.data) < eta] = 0
    A.eliminate_zeros()
    return A, eta

A, eta = build_blur_matrix(n1, n2, beta, tau)

b = A @ x + e

print(f"nnz(A) = {A.nnz}")
nonzero_pct = 100 * A.nnz / (A.shape[0] * A.shape[1])
print(f"Nonzero percentage of A = {nonzero_pct:.6f} %")
plt.spy(A, markersize=0.1)
plt.title("Sparsity pattern of A")
plt.show()


# 6.3 
# A, dense: 90000²*8 bytes = 64.8 GB
# A, sparse: (2 * nnz(A) + n + 1) * 8 bytes = (2 * 0.023171 * 90000² + 90000) * 4 bytes = 2.25 GB

def first_order_finite_difference(n):
    return sp.eye(n - 1, n) - sp.eye(n - 1, n, k=1)

L_a = sp.kron(sp.eye(n2), first_order_finite_difference(n1))
L_b = sp.kron(first_order_finite_difference(n2), sp.eye(n1))
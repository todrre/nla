import numpy as np
A = np.array([[4, -1, 0], [-1, 4, -1], [0, -1, 4]])
b = np.array([3, 6, 3])
D = np.diag(np.diag(A))
L = np.tril(A, -1)
U = np.triu(A, 1)
    
# Jacobi method
B_jacobi = -np.linalg.inv(D) @ (L + U)
print("Jacobi method matrix B:\n", B_jacobi)
print("Spectral radius:", max(abs(np.linalg.eigvals(B_jacobi))))

# Gauss-Seidel method
B_gauss_seidel = -np.linalg.inv(D + L) @ U
print("Gauss-Seidel method matrix B:\n", B_gauss_seidel)
print("Spectral radius:", max(abs(np.linalg.eigvals(B_gauss_seidel))))

def iterate(B, c, x0, n_iter, name):
    x = x0.astype(float)
    print(f"\n{name}:")
    print("x0 =", x)
    for k in range(1, n_iter + 1):
        x = B @ x + c
        print(f"x{k} =", x)
    return x


x0 = np.zeros(3)
c_jacobi = np.linalg.inv(D) @ b
c_gauss_seidel = np.linalg.inv(D + L) @ b
iterate(B_jacobi, c_jacobi, x0, 3, "Jacobi")
iterate(B_gauss_seidel, c_gauss_seidel, x0, 3, "Gauss-Seidel")
print("\nExact solution:", np.linalg.solve(A, b))

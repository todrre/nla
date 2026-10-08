import numpy as np
from PIL import Image

# Importera bilden, skala om till 300x300 och gör den gråskalig
n1 = n2 = 300
n = n1*n2
img = Image.open("img.png")
img = img.resize((300, 300))
img = img.convert("L")
beta = 5
B = beta * np.eye(n)
x = np.array(img.convert("L")).flatten()
t = range(n1*n2)

# Konstruera en blurring-matris
def blur1d(n, beta):
    t = np.linspace(0, 1, n)
    c = np.sqrt(beta / np.pi)
    h = t[1] - t[0]
    K = h * c * np.exp(-beta* (t[:, None] - t[None, :]) ** 2)

    return K / K.sum(axis=1, keepdims=True)


def build_blur_matrix(n1, n2, beta):
    A1 = blur1d(n1, beta)
    A2 = blur1d(n2, beta)
    A = np.kron(A1, A2)
    return A
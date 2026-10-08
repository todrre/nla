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
def blur_matrix(s, t, c, B):
    return c * np.exp(-(s-t).T @ B(s - t))


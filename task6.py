import numpy as np
from PIL import Image

# Importera bilden, skala om till 300x300 och gör den gråskalig
img = Image.open("temp-bild.png")
img = img.resize((300, 300))
img = img.convert("L")

beta_blur = 10.0
t = np.linspace(-1, 1, img.size[0])

c = np.sqrt(beta_blur / np.pi)
B = beta_blur * np.eye(img.size[0])
A = c * np.exp(-B * (t[:, None] - t[None, :]) ** 2)
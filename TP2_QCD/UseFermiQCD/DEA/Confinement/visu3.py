import numpy as np
betas = np.array([2.1, 2.3, 2.5])
lna = np.array([-1.386, -1.780, -2.401])
pente = np.polyfit(betas, lna, 1)
print(f"Pente mesurée : {pente[0]:.4f}")
print(f"Pente théorique : {-6*np.pi**2/(11*2):.4f}")
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

z = np.array([1, 2, 3, 4, 5, 6])
v = np.array([0.332702, 0.481335, 0.559518, 0.620637, 0.665874, 0.735336])
dv = np.array([0.000318, 0.00106, 0.00248, 0.00505, 0.00961, 0.0185])

def cornell(z, A0, A1, A2):
    return A0 + A1*z + A2/z

popt, pcov = curve_fit(cornell, z, v, p0=[0.2, 0.05, -0.3], sigma=dv)
perr = np.sqrt(np.diag(pcov))

A0, A1, A2 = popt
dA0, dA1, dA2 = perr

print(f"A0 = {A0:.4f} ± {dA0:.4f}")
print(f"A1 = {A1:.4f} ± {dA1:.4f}  (sigma*a^2)")
print(f"A2 = {A2:.4f} ± {dA2:.4f}  (terme coulombien e = -A2)")

# Calcul de a en fm
# A1 = sigma * a^2  avec sigma = (450 MeV)^2
# En unités naturelles : 1 fm = 1/197.33 MeV^-1
# sigma en fm^-2 : sigma = (450/197.33)^2
sigma = (450.0/197.33)**2   # en fm^-2
a = np.sqrt(A1/sigma)
da = 0.5 * dA1 / (sigma * a)  # propagation d'erreur

print(f"\nsigma = {sigma:.4f} fm^-2")
print(f"Maille a = {a:.4f} ± {da:.4f} fm")
print(f"Terme coulombien e = {-A2:.4f} ± {dA2:.4f}  (valeur expérimentale ~0.45)")

# Graphe
z_plot = np.linspace(0.5, 7, 100)
plt.figure(figsize=(8, 5))
plt.errorbar(z, v, yerr=dv, fmt='o', label='données', zorder=5)
plt.plot(z_plot, cornell(z_plot, *popt), label='fit Cornell', color='orange')
plt.xlabel('z (en unités de maille)')
plt.ylabel('aV(az)')
plt.title(f'Potentiel statique SU(2), β=2.5\na = {a:.4f} fm,  e = {-A2:.4f}')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('potential_beta2.5.png', dpi=150)
print('Figure sauvegardée')
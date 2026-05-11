import matplotlib.pyplot as plt
import numpy as np

jack = np.loadtxt('Jacknife.dat')
autocorrel = np.loadtxt('Autocorrel.dat')
autotime = np.loadtxt('Autotime.dat')

fig, axes = plt.subplots(1, 3, figsize=(15, 4))
fig.suptitle('Autocorrélation — SU(2), β=2.5, boucle 5×5', fontsize=12)

axes[0].plot(autocorrel[:,0], autocorrel[:,1], color='steelblue')
axes[0].axhline(0, color='gray', linestyle='--', linewidth=0.8)
axes[0].set_title('Fonction d\'autocorrélation A(τ)')
axes[0].set_xlabel('τ')
axes[0].set_ylabel('A(τ)')
axes[0].grid(True, alpha=0.3)

axes[1].plot(autotime[:,0], autotime[:,1], color='darkorange')
axes[1].set_title('Temps d\'autocorrélation intégré τ_int')
axes[1].set_xlabel('τ')
axes[1].set_ylabel('τ_int')
axes[1].grid(True, alpha=0.3)

axes[2].plot(jack[:,0], jack[:,1], color='green')
axes[2].set_title('Erreur Jackknife vs taille de bloc')
axes[2].set_xlabel('Taille de bloc')
axes[2].set_ylabel('Erreur Jackknife')
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('autocorrel.png', dpi=150, bbox_inches='tight')
print('Figure sauvegardée.')

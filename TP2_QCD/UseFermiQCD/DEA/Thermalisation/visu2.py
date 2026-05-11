import matplotlib.pyplot as plt
import numpy as np

hot = np.loadtxt('hot.dat')
cold = np.loadtxt('cold.dat')
bighot = np.loadtxt('bighot.dat')
bigcold = np.loadtxt('bigcold.dat')

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle('Thermalisation — SU(2), β=2.5, réseau 16×12×6×6', fontsize=13)

for ax, d_hot, d_cold, nstar, title in zip(
    axes,
    [hot, bighot], [cold, bigcold],
    [400, 800],
    ['Plaquette 1×1', 'Boucle 5×5']
):
    # Moyenne cumulée complète
    ax.plot(d_hot[:,0], d_hot[:,1], color='steelblue',
            linewidth=1.5, label='démarrage chaud')
    ax.plot(d_cold[:,0], d_cold[:,1], color='darkorange',
            linestyle='--', linewidth=1.5, label='démarrage froid')

    # Valeur moyenne après n* (ligne horizontale)
    mask_hot = d_hot[:,0] >= nstar
    mask_cold = d_cold[:,0] >= nstar
    mean_hot = d_hot[mask_hot, 1][-1]
    mean_cold = d_cold[mask_cold, 1][-1]
    ax.axhline(mean_hot, color='steelblue', linestyle='-',
               linewidth=0.8, alpha=0.5, label=f'plateau chaud = {mean_hot:.4f}')
    ax.axhline(mean_cold, color='darkorange', linestyle='-',
               linewidth=0.8, alpha=0.5, label=f'plateau froid = {mean_cold:.4f}')

    ax.axvspan(0, nstar, alpha=0.08, color='red')
    ax.axvline(x=nstar, color='red', linestyle=':', linewidth=1.2,
               label=f'n* ≈ {nstar}')

    ax.set_title(title, fontsize=11)
    ax.set_xlabel('Itérations')
    ax.set_ylabel('Valeur moyenne cumulée')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('thermalisation.png', dpi=150, bbox_inches='tight')
print('Figure sauvegardée.')
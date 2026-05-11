import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

hot = np.loadtxt('hot.dat')
cold = np.loadtxt('cold.dat')
bighot = np.loadtxt('bighot.dat')
bigcold = np.loadtxt('bigcold.dat')

fig = plt.figure(figsize=(12, 5))
fig.suptitle('Thermalisation — SU(2), β=2.5, réseau 12×12×6×6', fontsize=13)
gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.35)

ax1 = fig.add_subplot(gs[0])
ax1.plot(hot[:,0], hot[:,1], color='steelblue', label='démarrage chaud')
ax1.plot(cold[:,0], cold[:,1], color='darkorange', linestyle='--', label='démarrage froid')
ax1.axvline(x=50, color='gray', linestyle=':', linewidth=1, label='n* ≈ 50')
ax1.set_title('Plaquette 1×1', fontsize=11)
ax1.set_xlabel('Itérations')
ax1.set_ylabel('Valeur moyenne')
ax1.legend(fontsize=9)
ax1.grid(True, alpha=0.3)

ax2 = fig.add_subplot(gs[1])
ax2.plot(bighot[:,0], bighot[:,1], color='steelblue', label='démarrage chaud')
ax2.plot(bigcold[:,0], bigcold[:,1], color='darkorange', linestyle='--', label='démarrage froid')
ax2.axvline(x=700, color='gray', linestyle=':', linewidth=1, label='n* ≈ 700')
ax2.set_title('Boucle 5×5', fontsize=11)
ax2.set_xlabel('Itérations')
ax2.set_ylabel('Valeur moyenne')
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)

plt.savefig('thermalisation2.png', dpi=150, bbox_inches='tight')
print('Figure sauvegardée.')
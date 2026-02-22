import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("random.dat")
x = data[:,1]

plt.hist(x, bins=50, density=True)
plt.xlabel("x")
plt.ylabel("Density")
plt.title("Histogramme des tirages uniformes")
plt.savefig("histo.png")
print("Figure sauvegardée en histo.png")
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import os
import re

def parse_result_file(filepath):
    """Lit un fichier Result_* et extrait beta et le potentiel"""
    beta = None
    z_list, v_list, dv_list = [], [], []
    reading_potential = False

    with open(filepath, 'r') as f:
        for line in f:
            # Extraire beta
            if 'beta=2*Nc/g^2' in line:
                beta = float(re.search(r'beta=2\*Nc/g\^2:\s*([\d.]+)', line).group(1))

            # Début du bloc potentiel extrapolé
            if 'v(t,z) extrapolated' in line:
                reading_potential = True
                continue

            # Fin du bloc potentiel
            if reading_potential and 'Fit of potential' in line:
                reading_potential = False

            # Lire les valeurs
            if reading_potential:
                match = re.match(r'(\d+)\s+([\d.e+\-nan]+)\(\+/-\)([\d.e+\-nan]+)', line.strip())
                if match:
                    z_val = int(match.group(1))
                    v_val = match.group(2)
                    dv_val = match.group(3)
                    # Ignorer les nan et les valeurs aberrantes
                    if 'nan' not in v_val and 'nan' not in dv_val:
                        v_f = float(v_val)
                        dv_f = float(dv_val)
                        if dv_f > 0 and v_f > 0 and v_f < 10:
                            z_list.append(z_val)
                            v_list.append(v_f)
                            dv_list.append(dv_f)

    if beta and len(z_list) >= 3:
        return beta, np.array(z_list), np.array(v_list), np.array(dv_list)
    return None, None, None, None

def cornell(z, A0, A1, A2):
    return A0 + A1*z + A2/z

sigma_phys = (450.0/197.33)**2  # fm^-2

os.makedirs('plots', exist_ok=True)

# Trouver tous les fichiers Result_*
result_files = sorted([f for f in os.listdir('.') if re.match(r'Result_\d+', f)])
print(f"Fichiers trouvés : {result_files}\n")

results_summary = {}

for fname in result_files:
    beta, z, v, dv = parse_result_file(fname)

    if beta is None:
        print(f"{fname} : impossible de lire les données\n")
        continue

    if len(z) < 3:
        print(f"beta={beta} ({fname}) : pas assez de points ({len(z)})\n")
        continue

    print(f"{'='*45}")
    print(f"beta = {beta}  ({fname})  —  {len(z)} points")
    print(f"{'='*45}")

    try:
        popt, pcov = curve_fit(cornell, z, v, p0=[0.2, 0.05, -0.3], sigma=dv, maxfev=10000)
        perr = np.sqrt(np.diag(pcov))
        A0, A1, A2 = popt
        dA0, dA1, dA2 = perr

        if A1 <= 0:
            print(f"  A1={A1:.4f} négatif — fit non physique\n")
            continue

        a = np.sqrt(A1/sigma_phys)
        da = 0.5 * dA1 / (sigma_phys * a)

        results_summary[beta] = (a, da, -A2, dA2)

        print(f"  A0 = {A0:.4f} ± {dA0:.4f}")
        print(f"  A1 = {A1:.4f} ± {dA1:.4f}  (sigma*a^2)")
        print(f"  A2 = {A2:.4f} ± {dA2:.4f}")
        print(f"  Maille a  = {a:.4f} ± {da:.4f} fm")
        print(f"  Coulomb e = {-A2:.4f} ± {dA2:.4f}  (exp ~0.45)\n")

        # Graphe individuel
        z_plot = np.linspace(0.5, max(z)+1, 100)
        plt.figure(figsize=(8, 5))
        plt.errorbar(z, v, yerr=dv, fmt='o', label='données', zorder=5)
        plt.plot(z_plot, cornell(z_plot, *popt), label='fit Cornell', color='orange')
        plt.xlabel('z (en unités de maille)')
        plt.ylabel('aV(az)')
        plt.title(f'Potentiel statique SU(2), β={beta}\na = {a:.4f} fm,  e = {-A2:.4f}')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(f'plots/potential_beta{beta}.png', dpi=150)
        plt.close()

    except Exception as ex:
        print(f"  Fit échoué : {ex}\n")

# Graphe ln(a) vs beta — fenêtre de scaling
if len(results_summary) > 1:
    betas = np.array(sorted(results_summary.keys()))
    a_vals = np.array([results_summary[b][0] for b in betas])
    da_vals = np.array([results_summary[b][1] for b in betas])

    # Prédiction perturbative : pente théorique = -6*pi^2 / (11*Nc) avec Nc=2
    pente_theo = -6*np.pi**2 / (11*2)
    print(f"\nPente théorique perturbative : {pente_theo:.4f}")

    plt.figure(figsize=(8, 5))
    plt.errorbar(betas, np.log(a_vals), yerr=da_vals/a_vals,
                 fmt='o', color='steelblue', label='données')
    # Ligne guide avec pente théorique
    b_ref = betas[len(betas)//2]
    lna_ref = np.log(a_vals[len(betas)//2])
    b_line = np.linspace(betas[0], betas[-1], 100)
    plt.plot(b_line, lna_ref + pente_theo*(b_line - b_ref),
             '--', color='red', label=f'pente théorique = {pente_theo:.2f}')
    plt.xlabel('β')
    plt.ylabel('ln(a)')
    plt.title('Fenêtre de scaling : ln(a) vs β')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('plots/scaling_window.png', dpi=150)
    plt.close()
    print("Graphe de scaling sauvegardé dans plots/scaling_window.png")

print("\n=== Résumé ===")
for b in sorted(results_summary.keys()):
    a, da, e, de = results_summary[b]
    print(f"beta={b:.1f}  a={a:.4f}±{da:.4f} fm  e={e:.4f}±{de:.4f}")

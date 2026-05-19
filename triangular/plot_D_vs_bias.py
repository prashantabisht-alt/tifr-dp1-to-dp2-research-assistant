"""Plot D_even(b) and D_odd(b) for the chiral RTW on the triangular lattice.

Scans bias b from 0 to 1 at fixed gamma and gamma_r, using the Green-Kubo
method (which captures both symmetric and antisymmetric parts).
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from triangular_chiral_rtw import rates_from_gamma_bias, diffusion_tensor_green_kubo

# ---------- parameters ----------
v = 1.0
gamma_values = [0.5, 1.0, 2.0]        # total adjacent-turn rate
gamma_r = 0.2                          # reversal rate
biases = np.linspace(0.0, 1.0, 201)   # chirality bias

fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), sharex=True)

for gamma in gamma_values:
    D_even_arr = []
    D_odd_arr = []
    ratio_arr = []

    for b in biases:
        gp, gm = rates_from_gamma_bias(gamma, b)
        res = diffusion_tensor_green_kubo(v=v, gamma_plus=gp, gamma_minus=gm,
                                           gamma_r=gamma_r)
        D_even_arr.append(res['D_even'])
        D_odd_arr.append(res['D_odd'])
        ratio_arr.append(res['D_odd'] / res['D_even'] if res['D_even'] > 1e-12 else 0.0)

    label = rf'$\gamma = {gamma}$'

    axes[0].plot(biases, D_even_arr, label=label)
    axes[1].plot(biases, D_odd_arr, label=label)
    axes[2].plot(biases, ratio_arr, label=label)

axes[0].set_xlabel(r'bias $b$')
axes[0].set_ylabel(r'$D_{\rm even}$')
axes[0].set_title(r'Symmetric diffusivity')
axes[0].legend()

axes[1].set_xlabel(r'bias $b$')
axes[1].set_ylabel(r'$D_{\rm odd}$')
axes[1].set_title(r'Odd diffusivity')
axes[1].legend()

axes[2].set_xlabel(r'bias $b$')
axes[2].set_ylabel(r'$D_{\rm odd} / D_{\rm even}$')
axes[2].set_title(r'Odd-to-even ratio')
axes[2].legend()

fig.suptitle(rf'Chiral RTW on triangular lattice  ($v={v},\;\gamma_r={gamma_r}$)',
             fontsize=13)
fig.tight_layout()
fig.savefig('D_vs_bias.png', dpi=150)
print('Saved D_vs_bias.png')

# ---------- also scan gamma_r ----------
fig2, axes2 = plt.subplots(1, 3, figsize=(14, 4.5), sharex=True)
gamma_fixed = 1.0
gamma_r_values = [0.0, 0.2, 0.5, 1.0]

for gr in gamma_r_values:
    D_even_arr = []
    D_odd_arr = []
    ratio_arr = []

    for b in biases:
        gp, gm = rates_from_gamma_bias(gamma_fixed, b)
        res = diffusion_tensor_green_kubo(v=v, gamma_plus=gp, gamma_minus=gm,
                                           gamma_r=gr)
        D_even_arr.append(res['D_even'])
        D_odd_arr.append(res['D_odd'])
        ratio_arr.append(res['D_odd'] / res['D_even'] if res['D_even'] > 1e-12 else 0.0)

    label = rf'$\gamma_r = {gr}$'

    axes2[0].plot(biases, D_even_arr, label=label)
    axes2[1].plot(biases, D_odd_arr, label=label)
    axes2[2].plot(biases, ratio_arr, label=label)

axes2[0].set_xlabel(r'bias $b$')
axes2[0].set_ylabel(r'$D_{\rm even}$')
axes2[0].set_title(r'Symmetric diffusivity')
axes2[0].legend()

axes2[1].set_xlabel(r'bias $b$')
axes2[1].set_ylabel(r'$D_{\rm odd}$')
axes2[1].set_title(r'Odd diffusivity')
axes2[1].legend()

axes2[2].set_xlabel(r'bias $b$')
axes2[2].set_ylabel(r'$D_{\rm odd} / D_{\rm even}$')
axes2[2].set_title(r'Odd-to-even ratio')
axes2[2].legend()

fig2.suptitle(rf'Chiral RTW on triangular lattice  ($v={v},\;\gamma={gamma_fixed}$)',
              fontsize=13)
fig2.tight_layout()
fig2.savefig('D_vs_bias_gamma_r_scan.png', dpi=150)
print('Saved D_vs_bias_gamma_r_scan.png')

# ---------- print key values ----------
print('\n--- Key values at b=1 (fully chiral) ---')
for gamma in gamma_values:
    gp, gm = rates_from_gamma_bias(gamma, 1.0)
    res = diffusion_tensor_green_kubo(v=v, gamma_plus=gp, gamma_minus=gm, gamma_r=gamma_r)
    print(f'gamma={gamma}: D_even={res["D_even"]:.6f}, D_odd={res["D_odd"]:.6f}, '
          f'ratio={res["D_odd"]/res["D_even"]:.4f}')

print('\n--- Key values at b=1, varying gamma_r ---')
for gr in gamma_r_values:
    gp, gm = rates_from_gamma_bias(gamma_fixed, 1.0)
    res = diffusion_tensor_green_kubo(v=v, gamma_plus=gp, gamma_minus=gm, gamma_r=gr)
    print(f'gamma_r={gr}: D_even={res["D_even"]:.6f}, D_odd={res["D_odd"]:.6f}, '
          f'ratio={res["D_odd"]/res["D_even"]:.4f}')

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Load data
curvature_10 = np.load('curvature_10.npy')
tau_M = np.load('tau_M.npy')[1:]
Temperature = curvature_10['Temperature']
t_MC = curvature_10['t_MC']
mean_MC_T = curvature_10['mean_MC_T']
mean_GC_T = curvature_10['mean_GC_T']

# Figures
color_parula = plt.cm.parula(np.linspace(0,1,100))[::-1]
index_color = np.round(0.44 / Temperature * 100).astype(int)
color_order = color_parula[index_color]

tau_T = 1. / np.sqrt(3 * Temperature)

fig, ax1 = plt.subplots()
ax1.set_prop_cycle(color=color_order)

ax1.plot(t_MC / tau_T / 1e3, mean_MC_T, linewidth=2)

ax1.plot([1e-5, 1e5], [0, 0], '-', color='#666666')
ax1.plot([1, 1], [-1, 1], '--', color='#666666')

ax1.set_xscale('log')
ax1.set_xlim(10**-4, 10**4)
ax1.set_xticks(10.0 ** np.arange(-10, 11, 2))
ax1.set_ylim(-0.5, 0.5)

ax1.set_xlabel(r'$t/\tau_{M}$', fontsize=28)

ax1.tick_params(axis='both', which='major', labelsize=28)

fig.set_size_inches(6, 6, forward=True)
fig.tight_layout()

ax2 = ax1.twinx()
ax2.set_prop_cycle(color=color_order)

ax2.plot(t_MC / tau_T / 1e3, mean_GC_T, '--', linewidth=2)

ax2.plot([1e-8, 1e8], [0, 0], '-', color='#666666')
ax2.plot([1, 1], [-1, 1], '--', color='#666666')

ax2.set_xscale('log')
ax2.set_xlim(10**-4, 10**4)
ax2.set_xticks(10.0 ** np.arange(-10, 11, 2))
ax2.set_ylim(-0.5, 0.5)

ax2.set_xlabel(r'$tv_{RMS}$', fontsize=28)

ax2.tick_params(axis='both', which='major', labelsize=28)

fig.set_size_inches(6, 6, forward=True)
fig.tight_layout()

plt.show()
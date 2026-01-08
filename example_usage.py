import numpy as np
import matplotlib.pyplot as plt
from colossus.cosmology import cosmology
from colossus.lss import mass_function
from HMF import HaloMassFunction

# 1. Setup Cosmological Parameters
# The HaloMassFunction class sets 'planck15' by default
cosmo = cosmology.setCosmology('planck15')

# 2. Initialize the HMF Model
z = 0.0
mdef = '200m'
hmf_model = HaloMassFunction(z=z, mdef='m200b') # HMF.py uses 'm200b' internally for '200m' logic

# Define Mass Range
mass_range = np.logspace(10, 15.5, 100)

# Optional: define ps_args to match any custom power spectrum table
# If you have the Uchuu table, you can uncomment the following:
pk_path = 'power_uchuu_log10.dat'
ps_args = dict(model='uchuu_table', path=pk_path)
# Otherwise, Colossus will use the default 'camb' model
ps_args = None 

# 3. Calculate Halo Mass Functions
# New Model (Fernandez-Garcia 2025)
dndlnM_fg = hmf_model.n0(mass_range)

# Comparison models from Colossus
dndlnM_tinker = mass_function.massFunction(mass_range, z, mdef=mdef, model='tinker08', q_out='dndlnM')
dndlnM_watson = mass_function.massFunction(mass_range, z, mdef=mdef, model='watson13', q_out='dndlnM')

# 4. Visualization (Mass Function and Ratio)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), sharex=True, gridspec_kw={'height_ratios': [2, 1]})

# Top Plot: Mass Function
ax1.plot(mass_range, dndlnM_fg, label='Fernandez-Garcia 2025', color='blue', lw=2.5)
ax1.plot(mass_range, dndlnM_tinker, label='Tinker et al. 2008', color='red', ls='--')
ax1.plot(mass_range, dndlnM_watson, label='Watson et al. 2013', color='green', ls=':')

ax1.set_yscale('log')
ax1.set_ylabel(r'$dn/d\ln M \ [h^3 \rm Mpc^{-3}]$', fontsize=14)
ax1.set_title(f'Halo Mass Function Comparison at z = {z}', fontsize=16)
ax1.grid(True, which="both", ls="-", alpha=0.3)
ax1.legend(fontsize=12)

# Bottom Plot: Ratio
ax2.plot(mass_range, dndlnM_tinker / dndlnM_fg, color='red', label='Tinker / FG25')
ax2.plot(mass_range, dndlnM_watson / dndlnM_fg, color='green', label='Watson / FG25')
ax2.axhline(1.0, color='black', lw=1, ls='-')

ax2.set_xscale('log')
ax2.set_ylabel('Ratio to FG25', fontsize=14)
ax2.set_xlabel(r'Mass $[M_{\odot}/h]$', fontsize=14)
ax2.grid(True, which="both", ls="-", alpha=0.3)
ax2.legend(fontsize=10)

# Save and show
plt.tight_layout()
plt.subplots_adjust(hspace=0.05)
plt.savefig('hmf_comparison_plot.png')
print("Comparison plot saved as 'hmf_comparison_plot.png'")

plt.show()

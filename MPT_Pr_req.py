import numpy as np
import matplotlib.pyplot as plt

# ======  link parameters (copy-paste from your main file) ===============
Pt        = 10_000                    # transmitted RF feed power (W)
Dt, Dr    = 8.0, 6.0                  # dish diameters (m)
freqs     = [2.45e9, 5.80e9]          # carrier frequencies (Hz)
eff_tot   = [0.574, 0.328]            # end-to-end DC efficiencies (–)
# ========================================================================

c   = 3.0e8
At  = np.pi * Dt**2 / 4
Ar  = np.pi * Dr**2 / 4
d   = np.linspace(500, 10_000, 500)   # 0.5 … 10 km (m)

fig, ax = plt.subplots(figsize=(7.5, 4))

for f, eta_dc in zip(freqs, eff_tot):
    lam = c / f
    # received RF power after propagation *and* DC efficiencies
    Pr  = (At * Ar / (lam**2 * d**2)) * Pt * eta_dc
    # end-to-end system efficiency η = Pout / Pin
    eta = Pr / Pt                      # dimensionless (0–1)
    ax.semilogy(d/1000, eta*100,
                label = rf"{f/1e9:.2f} GHz")

ax.set_title("End-to-end link efficiency vs. range")
ax.set_xlabel("Range  $d$  (km)")
ax.set_ylabel("System efficiency  $η_{MPT}$  (%)")
ax.grid(True, which='both', ls=':')
ax.legend()
ax.set_ylim(1e-3, 10e1)      # 0.001 % … 10 %; tweak freely
fig.tight_layout()
plt.show()

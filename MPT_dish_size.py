import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
# ---------------------------------------------------------------
# USER INPUTS
Pt        = 10_000          # RF power at feed (W)
freq      = 2.45e9          # carrier (Hz)
eta_tot   = 0.574           # overall DC-to-DC efficiency (0–1)
d_link    = 5_000           # point-to-point distance (m)

# dish ranges to sweep  (edit as needed)
Dt_range  = np.linspace(2, 12, 50)    # transmitter dish diameters (m)
Dr_range  = np.linspace(2, 12, 50)    # receiver    dish diameters (m)
# ---------------------------------------------------------------

c   = 3.0e8                         # speed of light (m/s)
lam = c / freq                      # wavelength (m)

Dt_grid, Dr_grid = np.meshgrid(Dt_range, Dr_range, indexing='xy')

At = np.pi * (Dt_grid/2)**2         # effective apertures (m²)  (ideal, ηa≈1)
Ar = np.pi * (Dr_grid/2)**2

# Friis + total efficiency
Pr_dc = (At * Ar / (lam**2 * d_link**2)) * Pt * eta_tot   # 2-D array (W)

# ---- plotting --------------------------------------------------
fig, ax = plt.subplots(figsize=(6.5, 5))

levels_w = np.geomspace(1e-3, 1e3, 13)           # 1 mW … 1 kW
CS = ax.contourf(Dt_grid, Dr_grid, Pr_dc, levels=levels_w,
                 cmap='viridis', norm=LogNorm())  # Use LogNorm directly
cbar = fig.colorbar(CS, pad=0.03, label=r"Received DC power  $P_{r,DC}$  (W)")

ax.set_title(f"Received power at {d_link/1e3:.1f} km  (Pt={Pt/1e3:.0f} kW,  {freq/1e9:.2f} GHz)")
ax.set_xlabel("Tx dish diameter  $D_t$  (m)")
ax.set_ylabel("Rx dish diameter  $D_r$  (m)")
ax.grid(ls=':', color='0.8')
ax.set_aspect('equal', adjustable='box')          # nice squares

plt.tight_layout()
plt.show()

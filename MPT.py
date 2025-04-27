import numpy as np
import matplotlib.pyplot as plt

# ────────────────────────────────────────────────
# 0.  USER PARAMETERS (EDIT AS YOU WISH)
# ────────────────────────────────────────────────
Pt        = 10_000                      # W
Dt, Dr    = 8.0, 6.0                   # m
freqs     = [2.45e9, 5.80e9]           # Hz
eff_tot   = [0.574,   0.328]           # overall DC→DC efficiencies

# Laser-PV side
P_in_laser   = 10_000                  # W (wall-plug input)
eta_laser    = 0.50
eta_link     = 0.99
eta_PV       = 0.53
eta_tot_LPT  = eta_laser * eta_link * eta_PV     # ≈ 0.262
P_out_LPT_W  = P_in_laser * eta_tot_LPT          # constant delivered power (W)

# ────────────────────────────────────────────────
# 1.  CONSTANTS & DISTANCES
# ────────────────────────────────────────────────
c   = 3.0e8
At  = np.pi * Dt**2 / 4
Ar  = np.pi * Dr**2 / 4

d_m  = np.linspace(500, 10000, 500)   # metres

# helper: W → dBm (protect against log of zero/neg)
def watts_to_dBm(p_w, floor=1e-15):
    p_safe = np.maximum(p_w, floor)
    return 10*np.log10(p_safe / 1e-3)

# ────────────────────────────────────────────────
# 2.  PLOT (LOG Y-SCALE)
# ────────────────────────────────────────────────
fig, ax_w = plt.subplots(figsize=(8, 4))

# ▸ Microwave curves
for f, eta in zip(freqs, eff_tot):
    lam = c / f
    Pr  = (At * Ar / (lam**2 * d_m**2)) * Pt * eta
    ax_w.plot(d_m, Pr,
              label=f"{f/1e9:4.2f} GHz ($\eta_g \eta_t \eta_r \eta_d$ = {eta*100:4.1f} %)")

# ▸ Laser-PV (flat line)
P_LPT = np.full_like(d_m, P_out_LPT_W)
ax_w.plot(d_m, P_LPT, lw=2,
          label=f"Laser–PV link ($\eta_{{total}}$ = {eta_tot_LPT*100:.1f} %)")

# ── cosmetics
ax_w.set_title(f"Received / delivered power vs. distance (Pₜ = {Pt/1e3:.0f} kW)")
ax_w.set_xlabel("Centre-to-centre distance  $d$  (m)")
ax_w.set_ylabel("Power  (W)")
ax_w.set_yscale('log')
ax_w.grid(ls=":")
ax_w.legend()

# secondary dBm axis (derived from current watt limits)
ax_dbm = ax_w.twinx()
ymin, ymax = ax_w.get_ylim()
ax_dbm.set_ylim(watts_to_dBm(ymin), watts_to_dBm(ymax))
ax_dbm.set_ylabel("Power  (dBm)")

plt.tight_layout()
plt.show()

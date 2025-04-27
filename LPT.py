import numpy as np
import matplotlib.pyplot as plt

# ---------- user-editable LPT parameters ----------------------
P_in      = 10_000           # electrical input (W)
eta_laser = 0.50             # laser wall-plug efficiency
eta_link  = 0.99             # free-space transmission
eta_PV    = 0.53             # PV conversion efficiency
wavelength_nm = 801
mass_sp   = 200              # specific mass of complete LPT system (kg per kW_in)
# --------------------------------------------------------------

eta_tot = eta_laser * eta_link * eta_PV      # ≈ 0.262
P_out   = P_in * eta_tot

# --- delivered-power vs. range (just a flat line in this toy model)
d = np.linspace(0.5, 10, 400)        # km
P_vec = np.full_like(d, P_out)   # kW

fig, ax = plt.subplots(figsize=(8, 3.8))
ax.plot(d, P_vec, color='crimson', lw=2)

# ──--- the only new line ---───────────────
ax.set_yscale('log')          # make y-axis logarithmic
# ──────────────────────────────────────────

ax.set_xlabel("Range  $d$  (km)")
ax.set_ylabel("Delivered DC power  $P_{\\text{out}}$  (kW)")
ax.set_title("Laser-PV link : delivered power vs. range")
ax.grid(ls=':')

# ── annotation box ────────────────────────────────────────────
txt = (rf"$\eta_{{\text{{total}}}} = {eta_tot*100:.1f}\,\%$")

ax.text(0.97, 0.97, txt,
        transform=ax.transAxes, ha='right', va='top',
        bbox=dict(boxstyle="round,pad=0.3", fc="w", alpha=0.7, lw=0))

plt.tight_layout()
plt.show()



# --- 2. efficiency and system mass scaling ------------------
Pin_grid = np.linspace(1e3, 50e3, 200)      # 1 … 50 kW electrical input
Pout_grid = Pin_grid * eta_tot            # W
mass_grid = Pin_grid/1e3 * mass_sp          # kg   (linear with kW_in)

fig2, ax2 = plt.subplots(figsize=(8,4))

ax2.plot(Pin_grid/1e3, Pout_grid/1e3, label="Delivered power", color='navy')
ax2.set_xlabel("Electrical input  $P_{in}$  (kW)")
ax2.set_ylabel("Delivered DC power  $P_{out}$  (kW)", color='navy')
ax2.tick_params(axis='y', labelcolor='navy')
ax2.grid(ls=':')

ax3 = ax2.twinx()
ax3.plot(Pin_grid/1e3, mass_grid, label="System mass", color='darkorange')
ax3.set_ylabel("Total system mass  (kg)", color='darkorange')
ax3.tick_params(axis='y', labelcolor='darkorange')

ax2.set_title("LPT scaling : output power and hardware mass vs. input")
ax2.legend(loc='upper left')
ax3.legend(loc='lower right')
fig2.tight_layout()
plt.show()

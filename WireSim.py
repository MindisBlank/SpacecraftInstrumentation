import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

def plot_weight_vs_voltage_vs_power_for_1km_al(length,
                                               resistivity_aluminum,
                                               density_al,
                                               density_kapton,
                                               density_bn,
                                               density_mesh,
                                               E_dil,
                                               safety_factor):
    """
    3D surface of log10(Weight) vs. log10(Voltage) and log10(Power),
    for a 1 km Aluminum cable under a 5% loss limit.
    """
    # log-spaced voltages (1 kV → 100 kV) and powers (1 kW → 1 MW)
    voltages = np.logspace(3, 5, 50)
    powers   = np.logspace(3, 6, 100)

    # Prepare mesh
    Vg, Pg = np.meshgrid(voltages, powers)
    Wsurf  = np.zeros_like(Vg)

    # Fixed layers
    loss_limit       = 0.05
    thickness_bn     = 0.0001
    thickness_mesh   = 0.001

    # Fill the surface
    for i, P in enumerate(powers):
        for j, V in enumerate(voltages):
            t_kap = InsulationThickness(V, E_dil, safety_factor)
            r_cond = conductor_radius(
                resistivity_aluminum, V, length, P, loss_limit
            )
            vols = VolumeConductor(r_cond, t_kap, thickness_bn, thickness_mesh)
            mpm = ComputeMassPerMeter(
                vols[0], vols[1], vols[2], vols[3],
                density_al, density_kapton, density_bn, density_mesh
            )
            Wsurf[i, j] = mpm * length

    # Take logs for plotting
    logV = np.log10(Vg)
    logP = np.log10(Pg)
    logW = np.log10(Wsurf)

    # Plot
    fig = plt.figure(figsize=(10, 7))
    ax  = fig.add_subplot(111, projection='3d')

    surf = ax.plot_surface(
        logV, logP, logW,
        cmap='viridis',
        linewidth=0, antialiased=True,
    )

    ax.set_xlabel("log₁₀ Voltage (V)")
    ax.set_ylabel("log₁₀ Power (W)")
    ax.set_zlabel("log₁₀ Weight (kg)")
    ax.set_title("1 km Al Cable Weight vs. Voltage & Power\n(5% Loss Limit)")

    cbar = fig.colorbar(surf, shrink=0.5, aspect=10, label="log₁₀ Weight")
    plt.tight_layout()
    plt.show()


def plot_resistivity_vs_temperature():
    """
    Plots the resistivity of various materials versus temperature from -200°C to 250°C.
    Uses a linear model: ρ(T) = ρ₀ [1 + α (T - T₀)]
    """
    # Temperature range
    T = np.linspace(-200, 250, 500)

    # Material properties: (ρ₀ [Ω·m], α [1/°C], T₀ [°C])
    materials = {
        "Copper":    (1.7e-8,  0.00393, 20),
        "Aluminum":  (2.6e-8,  0.00403, 20),
        "Calcium":   (3.39e-8, 0.00420, 20),  # α ≈ 0.0042/°C
        "Beryllium": (3.99e-8, 0.00290, 20),  # α ≈ 0.0029/°C
        "Magnesium": (4.39e-8, 0.00390, 20),  # α ≈ 0.0039/°C
        "Silver":    (1.6e-8,  0.00380, 20),  # α ≈ 0.0038/°C
        "Gold":      (2.2e-8,  0.00340, 20),  # α ≈ 0.0034/°C
    }

    plt.figure(figsize=(10, 6))
    for name, (rho0, alpha, T0) in materials.items():
        rho = rho0 * (1 + alpha * (T - T0))
        plt.plot(T, rho, label=name)

    plt.xlabel("Temperature (°C)")
    plt.ylabel("Resistivity (Ω·m)")
    plt.title("Resistivity vs. Temperature for Various Metals")
    plt.yscale("log")            # optional: makes small differences easier to see
    plt.grid(True, which="both", ls="--")
    plt.legend()
    plt.show()


def InsulationThickness(V, E_dil, safety_factor):
    """
    Calculate the required insulation thickness in meters.

    Parameters:
    V (float): Voltage in volts (V)
    E_dil (float): Dielectric strength in V/m
    safety_factor (float): Safety factor (unitless)

    Returns:
    float: Insulation thickness in meters (m)
    """
    return V / E_dil * safety_factor

def conductor_radius(resistivity, voltage, length, power, loss_limit=0.05):
    """
    Calculate the minimum conductor radius needed to keep power loss under the specified limit.

    Parameters:
    resistivity (float): Resistivity of the conductor material in ohm meters (Ω·m)
    voltage (float): Transmission line voltage in volts (V)
    length (float): One-way distance of the transmission line in meters (m)
    power (float): Power being transmitted in watts (W)
    loss_limit (float): Maximum allowable power loss fraction (default: 0.05 for 5%)

    Returns:
    float: Required conductor radius in meters (m)
    """
    max_loss = loss_limit * power
    current = power / voltage
    max_resistance = max_loss / (current**2)
    required_radius = np.sqrt((resistivity * length) / (np.pi * max_resistance))
    return required_radius

def VolumeConductor(radius_conductor, thickness_L1, thickness_L2, thickness_L3):
    """
    Calculate the volume per unit length (m³ per m) for the conductor and its layers.

    Parameters:
    radius_conductor (float): Conductor radius in meters (m)
    thickness_L1 (float): Thickness of Layer 1 (Kapton insulation) in meters (m)
    thickness_L2 (float): Thickness of Layer 2 (BN nano-coating) in meters (m)
    thickness_L3 (float): Thickness of Layer 3 (Aluminum mesh) in meters (m)

    Returns:
    tuple: (volume_conductor, volume_L1, volume_L2, volume_L3) in cubic meters per meter (m³/m)
    """
    volume_conductor = np.pi * (radius_conductor**2)
    volume_L1 = np.pi * ((radius_conductor + thickness_L1)**2 - radius_conductor**2)
    volume_L2 = np.pi * ((radius_conductor + thickness_L1 + thickness_L2)**2 - (radius_conductor + thickness_L1)**2)
    volume_L3 = np.pi * ((radius_conductor + thickness_L1 + thickness_L2 + thickness_L3)**2 - (radius_conductor + thickness_L1 + thickness_L2)**2)
    return volume_conductor, volume_L1, volume_L2, volume_L3

def ComputeMassPerMeter(volume_conductor, volume_L1, volume_L2, volume_L3, density_al_conductor, density_kapton, density_bn, density_aluminum):
    """
    Compute the mass per unit length (kg/m) of the cable.

    Parameters:
    volume_conductor, volume_L1, volume_L2, volume_L3 (float): Volumes per meter (m³/m)
    density_al_conductor, density_kapton, density_bn, density_aluminum (float): Densities in kg/m³

    Returns:
    float: Total mass per meter in kilograms per meter (kg/m)
    """
    mass_conductor = volume_conductor * density_al_conductor
    mass_L1 = volume_L1 * density_kapton
    mass_L2 = volume_L2 * density_bn
    mass_L3 = volume_L3 * density_aluminum
    Tol_mass = mass_conductor + mass_L1 + mass_L2 + mass_L3
    return Tol_mass

def main():
    # Given parameters
    power       = 1e6      # 1 MW
    voltage     = 11000    # 30 kV
    safety_fac  = 5
    E_dil       = 272e6    # V/m

    # Conductor materials: ρ [Ω·m], density [kg/m³]
    materials = {
        "Aluminum":  (2.6e-8,  2700),
        "Copper":    (1.7e-8, 8960),
        "Calcium":   (3.39e-8,  1550),
        "Beryllium": (3.99e-8,  1848),
        "Magnesium": (4.39e-8,  1738),
        "Silver":    (1.6e-8,   10490),
        "Gold":      (2.2e-8,   19300),
    }

    # Insulation / coating layer densities (kg/m³)
    density_kapton = 1420
    density_bn     = 2100
    density_mesh   = 2700

    # Fixed layer thicknesses (m)
    thickness_bn      = 0.0001
    thickness_mesh    = 0.001
    thickness_kapton  = InsulationThickness(voltage, E_dil, safety_fac)

    # Range of distances
    distances = np.linspace(100, 10_000, 100)

    # Prepare a dict to hold weight-vs-distance for each metal
    weights = {metal: [] for metal in materials}

    # Compute each curve
    for L in distances:
        for metal, (rho, dens) in materials.items():
            # conductor radius for this metal & length
            r_cond = conductor_radius(rho, voltage, L, power, loss_limit=0.05)

            # volumes of conductor + 3 layers
            vols = VolumeConductor(r_cond,
                                   thickness_kapton,
                                   thickness_bn,
                                   thickness_mesh)

            # mass per meter
            m_per_m = ComputeMassPerMeter(vols[0], vols[1], vols[2], vols[3],
                                          dens,
                                          density_kapton,
                                          density_bn,
                                          density_mesh)

            # total weight for length L
            weights[metal].append(m_per_m * L)

    # Plot them all
    plt.figure(figsize=(9,6))
    for metal, curve in weights.items():
        plt.plot(distances, curve, lw=2, label=metal)

    plt.xlabel("Distance (m)")
    plt.ylabel("Cable Weight (kg)")
    plt.title(f"Transmission Cable Weight vs. Distance\n"
              f"({power/1e6:.0f} MW, {voltage/1e3:.0f} kV, 5% Loss)")
    plt.grid(True)
    plt.legend()
    plt.show()

 # Now produce the new plot requested:
# -    plot_weight_vs_voltage_for_1km_al(
# -        power=power,
# -        length=1000,
# -        resistivity_aluminum=2.6e-8,
# -        density_al=2700,
# -        density_kapton=density_kapton,
# -        density_bn=density_bn,
# -        density_mesh=density_mesh,
# -        E_dil=E_dil,
# -        safety_factor=5
# -    )
# replace your old 1-curve call with:
    plot_weight_vs_voltage_vs_power_for_1km_al(
        length=1000,
        resistivity_aluminum=2.6e-8,
        density_al=2700,
        density_kapton=density_kapton,
        density_bn=density_bn,
        density_mesh=density_mesh,
        E_dil=E_dil,
        safety_factor=5
    )



if __name__ == "__main__":
    main()
    #plot_resistivity_vs_temperature()
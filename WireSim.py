import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def plot_resistivity_vs_temperature():
    """
    Plots the resistivity of various materials versus temperature from -200°C to 250°C.
    
    The resistivity is computed using a linear model:
        ρ(T) = ρ₀ [1 + α (T - T₀)]
    where ρ₀ is the resistivity at the reference temperature T₀.
    
    Materials included:
      - Copper: ρ₀ = 1.68e-8 Ω·m, α = 0.00393/°C, T₀ = 20°C
      - Aluminum: ρ₀ = 2.82e-8 Ω·m, α = 0.00403/°C, T₀ = 20°C
      - Tungsten: ρ₀ = 5.60e-8 Ω·m, α = 0.0045/°C, T₀ = 20°C
      - Iron: ρ₀ = 9.71e-8 Ω·m, α = 0.005/°C, T₀ = 20°C
    """
    # Temperature range in °C
    T = np.linspace(-200, 250, 500)
    
    # Define material properties: (ρ₀ [Ω·m], α [1/°C], T₀ [°C])
    materials = {
        "Copper":   (1.68e-8, 0.00393, 20),
        "Aluminum": (2.65e-8, 0.00403, 20),
    }
    
    plt.figure(figsize=(10, 6))
    
    for name, (rho0, alpha, T0) in materials.items():
        # Calculate resistivity using the linear approximation
        rho = rho0 * (1 + alpha * (T - T0))
        # Clip any negative resistivity values to zero (not physical)
        rho = np.clip(rho, 0, None)
        plt.plot(T, rho, label=name)
    
    plt.xlabel("Temperature (°C)")
    plt.ylabel("Resistivity (Ω·m)")
    plt.title("Resistivity vs. Temperature for Various Materials")
    plt.legend()
    plt.grid(True)
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
    # Given parameters in SI base units (meters, kilograms, seconds)
    power = 1e6             # Power in watts (W) (1 MW)
    voltage = 30000         # Transmission voltage in volts (V)
    safety_factor = 5
    E_dil = 272e6           # Dielectric strength in V/m (272 V/µm converted to V/m)

    # Material properties in SI units
    # Resistivities (Ω·m)
    resistivity_aluminum = 3.73e-8  # Aluminum @ 121°C
    resistivity_copper   = 2.357e-8 # Copper @ 121°C

    # Densities (kg/m³)
    density_al       = 2700    # Aluminum conductor
    density_copper   = 8960    # Copper conductor
    density_kapton   = 1420    # Kapton insulation (from 1.42 g/cm³)
    density_bn       = 2100    # BN nano-coating (from 2.1 g/cm³)
    density_mesh     = 2700    # Aluminum mesh

    # Define constant layer thicknesses in meters:
    thickness_bn      = 0.0001   # BN nano-coating: 0.1 mm
    thickness_al_mesh = 0.001    # Aluminum mesh: 1.0 mm
    # Insulation thickness (Kapton) is computed based on voltage, dielectric strength, and safety factor.
    thickness_kapton  = InsulationThickness(voltage, E_dil, safety_factor)

    # --- Initial calculation for a 1000 m aluminum cable ---
    L_initial = 1000  # 1000 m cable length for initial report
    radius_conductor_al = conductor_radius(resistivity_aluminum, voltage, L_initial, power, loss_limit=0.05)
    volumes_al = VolumeConductor(radius_conductor_al, thickness_kapton, thickness_bn, thickness_al_mesh)
    mass_per_meter_al = ComputeMassPerMeter(volumes_al[0], volumes_al[1], volumes_al[2], volumes_al[3],
                                            density_al, density_kapton, density_bn, density_mesh)
    total_weight_al = mass_per_meter_al * L_initial
    total_radius_al = radius_conductor_al + thickness_kapton + thickness_bn + thickness_al_mesh
    print(f"--- For a {L_initial} m aluminum cable ---")
    print(f"Required conductor radius: {radius_conductor_al*1000:.3f} mm")
    print(f"Kapton insulation thickness: {thickness_kapton*1000:.3f} mm")
    print(f"Total cable radius: {total_radius_al*1000:.3f} mm")
    print(f"Total cable weight: {total_weight_al:.2f} kg")
    print(f"Total current: {power/voltage:.2f} A")

    # --- Plot total cable weight vs. cable length for Aluminum and Copper ---
    distances = np.linspace(100, 10000, 100)  # Cable lengths from 100 m to 10 km
    weights_al = []
    weights_cu = []

    for L in distances:
        # Aluminum cable calculations
        r_conductor_al = conductor_radius(resistivity_aluminum, voltage, L, power, loss_limit=0.05)
        vols_al = VolumeConductor(r_conductor_al, thickness_kapton, thickness_bn, thickness_al_mesh)
        m_per_meter_al = ComputeMassPerMeter(vols_al[0], vols_al[1], vols_al[2], vols_al[3],
                                             density_al, density_kapton, density_bn, density_mesh)
        weights_al.append(m_per_meter_al * L)
        
        # Copper cable calculations
        r_conductor_cu = conductor_radius(resistivity_copper, voltage, L, power, loss_limit=0.05)
        vols_cu = VolumeConductor(r_conductor_cu, thickness_kapton, thickness_bn, thickness_al_mesh)
        m_per_meter_cu = ComputeMassPerMeter(vols_cu[0], vols_cu[1], vols_cu[2], vols_cu[3],
                                             density_copper, density_kapton, density_bn, density_mesh)
        weights_cu.append(m_per_meter_cu * L)

    plt.figure(figsize=(8, 6))
    plt.plot(distances, weights_al, 'b-', lw=2, label='Aluminum Cable')
    plt.plot(distances, weights_cu, 'r-', lw=2, label='Copper Cable')
    plt.xlabel("Distance (m)")
    plt.ylabel("Cable Weight (kg)")
    plt.title(f"Lunar Transmission Line: Cable Weight vs. Distance\n({power/1e6} MW, {voltage/1e3} kV, Resistivity @ 121 °C; 5% Loss Limit)")
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()
    #plot_resistivity_vs_temperature()
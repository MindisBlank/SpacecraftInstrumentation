import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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
    Length = 1000           # Transmission line length in meters (m)
    power = 1e6             # Power in watts (W) (1 MW)
    voltage = 11000         # Transmission voltage in volts (V)
    safety_factor = 5
    E_dil = 272e6  # in V/m      # Dielectric strength in V/m (equivalent to 272 V/mm)

    # Material properties in SI units
    resistivity_aluminum = 2.82e-8  # Ω·m for aluminum
    density_kapton = 1420           # kg/m³ (1.42 g/cm³)
    density_bn = 2100               # kg/m³ (2.1 g/cm³)
    density_al = 2700               # kg/m³ (2.7 g/cm³)

    # Calculate required conductor radius in meters
    radius_conductor = conductor_radius(resistivity_aluminum, voltage, Length, power, loss_limit=0.05)
    print(f"Required conductor radius: {radius_conductor*1000:.3f} mm")

    # Calculate insulation thickness (Kapton) in meters
    thickness_kapton = InsulationThickness(voltage, E_dil, safety_factor)
    print(f"Kapton insulation thickness: {thickness_kapton*1000:.3f} mm")

    # Define thicknesses for BN nano-coating and Aluminum mesh in meters
    thickness_bn = 0.0001    # 0.1 mm = 0.0001 m
    thickness_al_mesh = 0.001  # 1.0 mm = 0.001 m

    # Compute total cable radius in meters
    total_radius = radius_conductor + thickness_kapton + thickness_bn + thickness_al_mesh

    # Compute volume per meter for each component
    volumes = VolumeConductor(radius_conductor, thickness_kapton, thickness_bn, thickness_al_mesh)

    # Compute mass per meter (kg/m) of the cable
    mass_per_meter = ComputeMassPerMeter(volumes[0], volumes[1], volumes[2], volumes[3],
                                         density_al, density_kapton, density_bn, density_al)

    # Total cable weight for the given length (kg)
    total_weight = mass_per_meter * Length
    print(f"Total cable weight for {Length} m: {total_weight:.2f} kg")
    print(f"Total cable radius: {total_radius*1000:.3f} mm")
    print(f"Total current: {power/voltage:.2f} A")
if __name__ == "__main__":
    main()

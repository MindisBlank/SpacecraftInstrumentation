import numpy as np
import matplotlib.pyplot as plt

# Physical constant
c = 3e8  # Speed of light in m/s

def dish_beamwidth(frequency, dish_diameter, dish_efficiency=0.6):
    """
    Compute the approximate half-power beamwidth for a parabolic dish.
    The dish gain is given by:
      G = dish_efficiency * (π * dish_diameter / λ)^2
    and beamwidth is approximated as:
      θ ≈ sqrt(4π / G)
    """
    wavelength = c / frequency
    gain = dish_efficiency * (np.pi * dish_diameter / wavelength)**2
    theta = np.sqrt(4 * np.pi / gain)
    return theta

def phased_array_beamwidth(frequency, elements_per_side, array_efficiency=0.8,
                           phase_error_std=0.0, calibration_eff=1.0, adaptive_nulling_gain=1.0):
    """
    Compute the beamwidth for a phased array while accounting for:
      - Phase errors (which reduce coherent summation via exp(-phase_error_std^2/2))
      - Calibration issues (modeled with calibration_eff, 0 to 1)
      - Adaptive nulling gain (an extra multiplicative factor)
      
    For a square array with half-wavelength spacing:
      - Total elements: N = (elements_per_side)^2
      - Effective aperture per element: (λ/2)^2
    The effective overall efficiency is the product of array_efficiency, phase error degradation,
    calibration_eff, and adaptive_nulling_gain.
    """
    wavelength = c / frequency
    N = elements_per_side**2
    d = wavelength / 2  # half-wavelength spacing
    # Combined efficiency factor: ideal array efficiency degraded by phase errors and calibration issues.
    effective_efficiency = (array_efficiency
                            * np.exp(-phase_error_std**2 / 2)
                            * calibration_eff
                            * adaptive_nulling_gain)
    A_e = effective_efficiency * N * (d**2)
    gain = 4 * np.pi * A_e / wavelength**2
    theta = np.sqrt(4 * np.pi / gain)
    return theta

def nearfield_distance_dish(frequency, dish_diameter):
    """
    Approximate near-field (Rayleigh) distance for a dish:
      R_nf ≈ 2 * D^2 / λ
    Below this distance, the far-field model is invalid.
    """
    wavelength = c / frequency
    return 2 * (dish_diameter**2) / wavelength

def nearfield_distance_array(frequency, elements_per_side):
    """
    Approximate near-field distance for a phased array.
    We treat the array as having an 'effective diameter' ~ side length.
      side_length = elements_per_side * (λ/2)
    Then
      R_nf ≈ 2 * side_length^2 / λ
    """
    wavelength = c / frequency
    side_length = elements_per_side * (wavelength / 2)
    return 2 * (side_length**2) / wavelength

def power_received(P_tx, distance, frequency, rectenna_area,
                   antenna_type='dish', dish_diameter=None, elements_per_side=None,
                   dish_efficiency=0.6, array_efficiency=0.8,
                   phase_error_std=0.0, calibration_eff=1.0, adaptive_nulling_gain=1.0,
                   tx_efficiency=0.7, rectenna_efficiency=0.85):
    """
    Calculate the received power based on the transmitter power, antenna beam characteristics,
    and nonidealities. We also skip near-field distances and cap captured area to avoid unphysical
    over-unity at short range.
    """
    # --- 1) Check near-field distance ---
    if antenna_type == 'dish':
        if dish_diameter is None:
            raise ValueError("dish_diameter must be provided for dish type.")
        nf_dist = nearfield_distance_dish(frequency, dish_diameter)
        theta = dish_beamwidth(frequency, dish_diameter, dish_efficiency)
    elif antenna_type == 'phased_array':
        if elements_per_side is None:
            raise ValueError("elements_per_side must be provided for phased_array type.")
        nf_dist = nearfield_distance_array(frequency, elements_per_side)
        theta = phased_array_beamwidth(frequency, elements_per_side, array_efficiency,
                                       phase_error_std, calibration_eff, adaptive_nulling_gain)
    else:
        raise ValueError("Invalid antenna_type. Use 'dish' or 'phased_array'.")

    # If we are inside the near-field region, return NaN (or zero) to indicate invalid far-field assumption.
    if distance < nf_dist:
        return np.nan  # or return 0.0, whichever you prefer

    # --- 2) Compute far-field beam and cap captured power ---
    # Compute the beam's radius at the receiver assuming circular symmetry.
    radius = (theta / 2) * distance
    beam_area = np.pi * radius**2
    # Apply transmitter efficiency.
    P_rf = P_tx * tx_efficiency
    # Power density at the receiver (far-field).
    power_density = P_rf / beam_area

    # Cap the captured area to ensure we can't capture more power than is present in the beam cross section.
    A_captured = min(rectenna_area, beam_area)
    P_captured = power_density * A_captured

    return P_captured * rectenna_efficiency

def simulate_transmission():
    """
    Compare the received power for a parabolic dish and a phased array
    while avoiding near-field distances and capping the captured area.
    """
    distances = np.linspace(10, 5e3, 200)  # from 10 m to 10 km
    frequency = 5.8e9        # 5.8 GHz
    P_tx = 10000             # 10 kW transmitter power
    rectenna_area = 50.0     # m²
    tx_efficiency = 0.7
    rectenna_efficiency = 0.85

    # Dish antenna parameters
    dish_diameter = 1.5      # 1.5 m dish diameter
    dish_efficiency = 0.6

    # Phased array parameters (ideal case)
    elements_per_side = 64
    array_efficiency = 0.8

    received_powers_dish = []
    received_powers_array_ideal = []

    for d in distances:
        # Dish antenna received power
        power_dish = power_received(
            P_tx=P_tx, distance=d, frequency=frequency,
            rectenna_area=rectenna_area,
            antenna_type='dish',
            dish_diameter=dish_diameter,
            dish_efficiency=dish_efficiency,
            tx_efficiency=tx_efficiency,
            rectenna_efficiency=rectenna_efficiency
        )
        # Phased array received power (ideal beamforming)
        power_array_ideal = power_received(
            P_tx=P_tx, distance=d, frequency=frequency,
            rectenna_area=rectenna_area,
            antenna_type='phased_array',
            elements_per_side=elements_per_side,
            array_efficiency=array_efficiency,
            phase_error_std=0.0,
            calibration_eff=1.0,
            adaptive_nulling_gain=1.0,
            tx_efficiency=tx_efficiency,
            rectenna_efficiency=rectenna_efficiency
        )

        received_powers_dish.append(power_dish)
        received_powers_array_ideal.append(power_array_ideal)

    plt.figure(figsize=(9, 6))
    plt.plot(distances / 1e3, received_powers_dish,
             label='Parabolic Dish (1.5 m diameter)', linestyle='--')
    plt.plot(distances / 1e3, received_powers_array_ideal,
             label='Phased Array (64x64 elements)')
    plt.xlabel('Distance (km)')
    plt.ylabel('Received Power (W)')
    plt.title('Microwave Power Transmission')
    plt.yscale('log')
    plt.grid(True, which="both", ls="--")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    simulate_transmission()

import numpy as np
import matplotlib.pyplot as plt

# Physical constant
c = 3e8  # Speed of light in m/s

def beam_radius(distance, wavelength, transmitter_diameter):
    """
    Estimate the beam radius (1/e radius) at a given distance.
    Uses the diffraction limit:
      theta ≈ 1.22 * (lambda / D)
      beam_radius ≈ (theta/2) * distance
    Returns the beam radius in meters.
    """
    theta = 1.22 * (wavelength / transmitter_diameter)  # full angle divergence in radians
    radius = (theta/2) * distance  # half-angle approximation for beam radius
    return radius

def power_received(P_tx, distance, frequency, transmitter_diameter, rectenna_area,
                   tx_efficiency=0.7, rectenna_efficiency=0.85):
    """
    Calculate the received power at the rectenna.

    Parameters:
      P_tx              : Transmitter electrical power (W)
      distance          : Distance between transmitter and receiver (m)
      frequency         : Operating frequency (Hz)
      transmitter_diameter: Diameter of the transmitting antenna (m)
      rectenna_area     : Effective area of the rectenna (m²)
      tx_efficiency     : Efficiency of the transmitter power conversion (default 70%)
      rectenna_efficiency: Efficiency of the rectenna conversion (default 85%)

    Returns:
      Received power in watts.
    """
    wavelength = c / frequency
    # RF power after transmitter conversion
    P_rf = P_tx * tx_efficiency
    
    # Estimate beam radius at the given distance
    r = beam_radius(distance, wavelength, transmitter_diameter)
    beam_area = np.pi * r**2   # beam cross-sectional area
    
    # Calculate power density and then captured power
    power_density = P_rf / beam_area
    P_captured = power_density * rectenna_area
    
    # DC power after rectenna conversion
    P_received = P_captured * rectenna_efficiency
    return P_received

def simulate_basic():
    """
    Basic simulation plotting received power vs. distance
    with logarithmic y-axis.
    """
    # Parameters
    frequency_GHz = 5.8           # Frequency in GHz
    frequency = frequency_GHz * 1e9  # Convert to Hz
    transmitter_diameter = 3.0    # in meters
    rectenna_area = 50.0          # in m²
    P_tx = 10000.0              # Transmitter power in watts (10 kW)
    tx_efficiency = 0.70
    rectenna_efficiency = 0.85

    distances = np.linspace(10, 10e3, 100)  # from 10 m to 10 km
    received_powers = np.array([power_received(P_tx, d, frequency, transmitter_diameter,
                                                 rectenna_area, tx_efficiency, rectenna_efficiency)
                                  for d in distances])

    plt.figure(figsize=(8, 5))
    plt.plot(distances/1e3, received_powers, label=f'{frequency_GHz} GHz, {transmitter_diameter}m Dish')
    plt.xlabel("Distance (km)")
    plt.ylabel("Received Power (W)")
    plt.title("Microwave Power Transmission on the Moon")
    plt.yscale('log')
    plt.grid(True, which="both", ls="--")
    plt.legend()
    plt.tight_layout()
    plt.show()

def simulate_parameter_tuning():
    """
    Extended simulation to tune parameters.
    Plots received power vs. distance for different frequencies,
    transmitter diameters, and rectenna areas.
    """
    # Common simulation settings
    distances = np.linspace(10, 10e3, 100)  # distance in meters
    P_tx = 10000.0              # Transmitter power (W)
    tx_efficiency = 0.70
    rectenna_efficiency = 0.85

    # Define ranges for parameters to sweep
    frequencies_GHz = [2.45, 5.8, 10]  # in GHz, typical operating values
    transmitter_diameters = [2.0, 3.0, 4.0]  # in meters
    rectenna_areas = [25, 50, 75]  # in m²

    # Create a figure with subplots for each parameter sweep
    fig, axs = plt.subplots(3, 1, figsize=(8, 15), sharex=True)
    
    # Sweep frequencies while keeping other parameters fixed (use mid-values)
    fixed_tx_diameter = 3.0
    fixed_rectenna_area = 50.0
    for freq in frequencies_GHz:
        frequency = freq * 1e9
        received_powers = np.array([power_received(P_tx, d, frequency, fixed_tx_diameter,
                                                     fixed_rectenna_area, tx_efficiency, rectenna_efficiency)
                                    for d in distances])
        axs[0].plot(distances/1e3, received_powers, label=f'{freq} GHz')
    axs[0].set_ylabel("Received Power (W)")
    axs[0].set_title("Varying Frequency (Tx Dish=3.0m, Rectenna=50m²)")
    axs[0].set_yscale('log')
    axs[0].grid(True, which="both", ls="--")
    axs[0].legend()

    # Sweep transmitter diameters while keeping other parameters fixed
    fixed_frequency = 5.8  # GHz
    frequency = fixed_frequency * 1e9
    for d_tx in transmitter_diameters:
        received_powers = np.array([power_received(P_tx, d, frequency, d_tx,
                                                     fixed_rectenna_area, tx_efficiency, rectenna_efficiency)
                                    for d in distances])
        axs[1].plot(distances/1e3, received_powers, label=f'{d_tx} m')
    axs[1].set_ylabel("Received Power (W)")
    axs[1].set_title("Varying Transmitter Dish Diameter (Freq=5.8GHz, Rectenna=50m²)")
    axs[1].set_yscale('log')
    axs[1].grid(True, which="both", ls="--")
    axs[1].legend()

    # Sweep rectenna areas while keeping other parameters fixed
    for area in rectenna_areas:
        received_powers = np.array([power_received(P_tx, d, frequency, fixed_tx_diameter,
                                                     area, tx_efficiency, rectenna_efficiency)
                                    for d in distances])
        axs[2].plot(distances/1e3, received_powers, label=f'{area} m²')
    axs[2].set_xlabel("Distance (km)")
    axs[2].set_ylabel("Received Power (W)")
    axs[2].set_title("Varying Rectenna Area (Freq=5.8GHz, Tx Dish=3.0m)")
    axs[2].set_yscale('log')
    axs[2].grid(True, which="both", ls="--")
    axs[2].legend()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    print("Running basic simulation...")
    simulate_basic()
    print("Running parameter tuning simulation...")
    simulate_parameter_tuning()

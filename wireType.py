import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#TODO make more advance add more parameters and use some optimazation method to find the optimal maybe SGD

#TODO Make plots that are usefull
#TODO make things easier find first material for the cable first, look at the ratio between p/density and look at area^2 vs weight and think of insulation

# Constants for material properties (Resistivity in Ω·m, Density in kg/m³)
materials = {
    "Copper (Cu)": {"resistivity": 1.68e-8, "density": 8960},
    "Aluminum (Al)": {"resistivity": 2.82e-8, "density": 2700},
    "Silver (Ag)": {"resistivity": 1.59e-8, "density": 10500},
}

# Given set points
voltage_levels = [1000, 5000, 10000]  # Transmission voltage in Volts
rated_power = [50000, 500000, 1000000]  # Power in Watts (50 kW, 500 kW, 1 MW)
lengths = [1000, 5000, 10000]  # Transmission length in meters

# Compute power loss and weight for each configuration
results = []

for material, props in materials.items():
    for voltage in voltage_levels:
        for power in rated_power:
            for length in lengths:
                current = power / voltage  # I = P / V

                # Assume circular conductor, solve for optimal area A
                area = (props["resistivity"] * length * current) / (voltage * 0.01)  # Arbitrary 1% voltage drop
                radius = np.sqrt(area / np.pi)
                weight = area * length * props["density"]  # Total weight

                # Resistance of wire (Ohm's Law)
                resistance = props["resistivity"] * length / area
                power_loss = current**2 * resistance

                results.append([material, voltage, power, length, resistance, power_loss, weight, radius * 2])

# Convert results into DataFrame and display
df = pd.DataFrame(results, columns=["Material", "Voltage (V)", "Power (W)", "Length (m)", 
                                    "Resistance (Ω)", "Power Loss (W)", "Weight (kg)", "Diameter (m)"])

print(df)  # Basic print

# Plot 1: Power Loss vs Transmission Length
plt.figure(figsize=(10, 5))
sns.lineplot(data=df, x="Length (m)", y="Power Loss (W)", hue="Material", style="Voltage (V)", markers=True)
plt.title("Power Loss vs Transmission Length")
plt.xlabel("Transmission Length (m)")
plt.ylabel("Power Loss (W)")
plt.legend(title="Material & Voltage")
plt.show()

# Plot 2: Weight vs Transmission Length
plt.figure(figsize=(10, 5))
sns.lineplot(data=df, x="Length (m)", y="Weight (kg)", hue="Material", style="Voltage (V)", markers=True)
plt.title("Cable Weight vs Transmission Length")
plt.xlabel("Transmission Length (m)")
plt.ylabel("Weight (kg)")
plt.legend(title="Material & Voltage")
plt.show()

# Plot 3: Resistance vs Transmission Length
plt.figure(figsize=(10, 5))
sns.lineplot(data=df, x="Length (m)", y="Resistance (Ω)", hue="Material", style="Voltage (V)", markers=True)
plt.title("Cable Resistance vs Transmission Length")
plt.xlabel("Transmission Length (m)")
plt.ylabel("Resistance (Ω)")
plt.legend(title="Material & Voltage")
plt.show(
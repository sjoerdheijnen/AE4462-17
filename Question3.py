import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import griddata


# CONSTANTS 
k_HO2_HO2 = 6.1e-12  # [cm³/molecule/s]
k_OH_NO2 = 9.0e-12  # [cm³/molecule/s]
k_CO_OH = 1.57e-13  # [cm³/molecule/s + cair * 3.54E-33]
k_HO2_NO = 3.3e-12 * np.exp(270 / 298)  # [cm³/molecule/s]
c_air = 2.69e19  # [molecules/cm³]

p = 101300 # [Pa]
T = 298 # [K]
NO2_NO_ratio = 7 # [-]

COmix = 4500 # [nmol * mol^-1]
COmix = COmix *  1e-9 * c_air  # converting to [molecules/cm³]

PHOx = 1.2e-3 # [nmol * mol^-1 * s^-1]
PHOx = PHOx * 1e-9 * c_air # [molecules/cm³/s]

NO_range = np.arange(0, 6e10, 1e7)
NO2_range = NO_range * NO2_NO_ratio

# QUADRATIC EQUAION CONSTANTS
a_values = 2 * k_HO2_HO2 * (1 + (k_OH_NO2 * NO2_range) / (k_CO_OH * COmix))
b_values = (k_HO2_NO * k_OH_NO2 * NO2_range * NO_range) / (k_CO_OH * COmix)
c_value = -PHOx

# QUADRATIC EQUATION SOLVER
HO2_values = (-b_values + np.sqrt(b_values**2 - 4 * a_values * c_value)) / (2 * a_values)

# PO3 EQUATION SOLVER
PO3_values = k_HO2_NO * HO2_values * NO_range


# Plot 1: [HO2] vs [NO]
plt.figure(figsize=(8, 6))
plt.plot(NO_range, HO2_values, label="[HO2] vs [NO]", color='b')
plt.xlabel("[NO] (molecules/cm³)")
plt.ylabel("[HO2] (molecules/cm³)")
plt.title("Dependency of [HO2] on [NO]")
plt.legend()
plt.grid(True)
plt.show()

# Plot 2: PO3 vs [NO]
plt.figure(figsize=(8, 6))
plt.plot(NO_range, PO3_values, label="PO3 vs [NO]", color='r')
plt.xlabel("[NO] (molecules/cm³)")
plt.ylabel("PO3 (molecules/cm³/s)")
plt.title("Dependency of Ozone Production PO3 on [NO]")
plt.legend()
plt.grid(True)
plt.show()

# Load Data
files = {
    "Andechs_NO": "O3_NO_NO2_measurements/airbase_DEBY109_NO.csv",
    "Andechs_NO2" :"O3_NO_NO2_measurements/airbase_DEBY109_NO2.csv",
    "Andechs_O3": "O3_NO_NO2_measurements/airbase_DEBY109_O3.csv",
    "Rotterdam_NO": "O3_NO_NO2_measurements/airbase_NL00418_NO.csv",
    "Rotterdam_NO2" :"O3_NO_NO2_measurements/airbase_NL00418_NO2.csv",
    "Rotterdam_O3" :"O3_NO_NO2_measurements/airbase_NL00418_O3.csv",
    "Vredepeel_NO" :"O3_NO_NO2_measurements/airbase_NL00131_NO.csv",
    "Vredepeel_NO2" :"O3_NO_NO2_measurements/airbase_NL00131_NO2.csv",
    "Vredepeel_O3" :"O3_NO_NO2_measurements/airbase_NL00131_O3.csv"
}
data = {name: pd.read_csv(path) for name, path in files.items()}

stations = ["Andechs", "Vredepeel", "Rotterdam"]
station_data = {}

for station in stations:
    data_NO = data[f"{station}_NO"]
    data_NO2 = data[f"{station}_NO2"]
    data_O3 = data[f"{station}_O3"]
    
    # Merge and alter data
    merged_data = data_NO.merge(data_NO2, on="time").merge(data_O3, on="time")
    merged_data.columns = ["datetime", "NO", "NO2", "O3"]
    merged_data["datetime"] = pd.to_datetime(merged_data["datetime"])
    
    # Compute NOx concentration as sum of NO and NO2
    merged_data["NOx"] = merged_data["NO"] + merged_data["NO2"]
    
    # Store in dictionary
    station_data[station] = merged_data



# # Create scatter plots for [O3] vs. [NOx]
# plt.figure(figsize=(10, 6))

# for station, data in station_data.items():
#     plt.scatter(data["NOx"], data["O3"], label=station, alpha=0.3, marker=".", s=5)

# plt.xscale("log")  # Log scale for NOx
# plt.xlabel("[NOx] (µg/m³)")
# plt.ylabel("[O3] (µg/m³)")
# plt.title("Scatter Plot of [O3] vs. [NOx] for Different Stations")
# plt.legend()
# plt.grid(True)
# plt.show()


# Create a figure with 3 subplots (1 row, 3 columns)
fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)  # Share y-axis for better comparison

# Loop through each station and assign each to a subplot
for i, (station, data) in enumerate(station_data.items()):
    axes[i].scatter(data["NOx"], data["O3"], label=station, alpha=0.3, marker=".", s=5)
    axes[i].set_xscale("log")  # Log scale for NOx
    axes[i].set_xlabel("[NOx] (µg/m³)")
    axes[i].set_title(f"{station} Station", fontsize=12, fontweight="bold")
    axes[i].grid(True)

# Set shared y-axis label
axes[0].set_ylabel("[O3] (µg/m³)")

# Adjust layout for better spacing
plt.tight_layout()

# Show the figure with three subplots
plt.show()
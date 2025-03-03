import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.interpolate import griddata

data = pd.read_csv("yearmean_RD1_2019_zm.csv")

lat = data["lat"].values
p = data["plev"].values # Pa
O3 = data["O3"].values  # mol/mol
T = data["tm1"].values # K
aps = data["aps"].values # Pa
R = 8.314 # J/(mol*K)

log_p = np.log10(p) #Pa on log scale

concentration = []

for i in range(len(lat)):
    C_O3 = (O3[i] * p[i]) / (R * T[i]) #ideal gas law
    concentration.append(C_O3)

concentration = np.array(concentration)

## contour plots

# grid generation
lat_grid = np.linspace(min(lat), max(lat), 100) #equator is at zero with the poles at the both ends
plev_grid = np.linspace(min(log_p), max(log_p), 100)
lat_grid_mesh, plev_grid_mesh = np.meshgrid(lat_grid, plev_grid)

O3_grid = griddata((lat, log_p), O3, (lat_grid_mesh, plev_grid_mesh), method='linear')
C_O3_grid = griddata((lat, log_p), concentration, (lat_grid_mesh, plev_grid_mesh), method='linear')

# mixing ratio plot
plt.figure(figsize=(10, 6))
levels = np.logspace(np.log10(min(O3)), np.log10(max(O3)), num=20)
contour = plt.contourf(lat_grid, 10**plev_grid, O3_grid, levels=levels, cmap="seismic", norm=mcolors.LogNorm())

#contour = plt.contourf(lat_grid, 10**plev_grid, O3_grid, levels=100, cmap="viridis", norm=mcolors.LogNorm())
cbar = plt.colorbar(contour)
cbar.set_label("Ozone Mixing Ratio [mol/mol]")

plt.gca().invert_yaxis()  # inverted to simulate the actual atmosphere(high pressure near ground)
plt.yscale("log")
plt.xlabel("Latitude")
plt.ylabel("Pressure [Pa]")
plt.title("Ozone Mixing Ratio in Latitude-Pressure Grid (2019)")

# concentration plot
plt.figure(figsize=(10, 6))
levels = np.logspace(np.log10(min(concentration)), np.log10(max(concentration)), num=20)
contour = plt.contourf(lat_grid, 10**plev_grid, C_O3_grid, levels=levels, cmap="seismic", norm=mcolors.LogNorm())

#contour = plt.contourf(lat_grid, 10**plev_grid, O3_grid, levels=100, cmap="viridis", norm=mcolors.LogNorm())
cbar = plt.colorbar(contour)
cbar.set_label("Ozone Concentration [mol/m3]")

plt.gca().invert_yaxis()
plt.yscale("log")
plt.xlabel("Latitude")
plt.ylabel("Pressure [Pa]")
plt.title("Ozone Concentration in Latitude-Pressure Grid (2019)")

for i in range(len(O3)):
    if math.isnan(O3[i]):
        O3[i]=0
        concentration[i]=0

# find high ozone concentrations
threshold = np.percentile(O3, 95)  # top 5%
threshold2 = np.percentile(concentration, 95)  # top 5%
ozone_layer = p[O3 >= threshold]
ozone_layer2 = p[concentration >= threshold2]

# atmosphere calculation constants
g0 = 9.80665  # m/s^2
R_air = 287.0  # J/(kg*K)
T0 = 288.15  # K
p0 = 101325.0  # Pa

# atmospheric layers and temperature gradients (K/m)
layers = [
    (-0.0065, 0, 11000),
    (0.0, 11000, 20000),
    (0.0010, 20000, 32000),
    (0.0028, 32000, 47000),
    (0.0, 47000, 51000),
    (-0.0028, 51000, 71000),
    (-0.0020, 71000, 86000)
]


def pressure_to_altitude(p):
    h = 0
    T = T0
    p_current = p0
    for a, h_start, h_end in layers:
        if a != 0:
            T_next = T + a * (h_end - h_start)
            p_next = p_current * (T_next / T) ** (-g0 / (a * R_air))
        else:
            p_next = p_current * np.exp(-g0 * (h_end - h_start) / (R_air * T))

        if p >= p_next:
            if a!=0:
                return h_start + (T / a) * (1 - (p / p_current) ** (a * R_air / g0))
            else:
                return h_start + (p - p_current) * R_air * T / g0 #this is incorrect, the zero lapse rate is causing issues

        T = T_next
        p_current = p_next
    return h_end


altitudes = np.array([pressure_to_altitude(p) for p in ozone_layer])
altitudes2 = np.array([pressure_to_altitude(p) for p in ozone_layer2])

print(f"Ozone layer altitude range based on mixing ratio: {np.min(altitudes)/1000:.2f} km - {np.max(altitudes)/1000:.2f} km")
print(f"Ozone layer altitude range based on concentration: {np.min(altitudes2)/1000:.2f} km - {np.max(altitudes2)/1000:.2f} km") #still incorrect

plt.show()
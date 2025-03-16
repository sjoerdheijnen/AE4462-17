import pandas as pd
import numpy as np
import math
import os
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

log_p = p#np.log10(p)

concentration = []

for i in range(len(lat)):
    C_O3 = (O3[i] * p[i]) / (R * T[i]) #ideal gas law
    concentration.append(C_O3)

concentration = np.array(concentration)

Nlat = len(np.unique(lat))
Np = len(np.unique(log_p))

O3 = O3[:Nlat * Np]
O3_range = [min(O3),max(O3)]
lat = lat[:Nlat * Np]
log_p = log_p[:Nlat * Np]
concentration = concentration[:Nlat * Np]
C_range = [min(concentration),max(concentration)]


unique_lat = np.sort(np.unique(lat))[::-1]

unique_log_p = np.unique(log_p)

lat_grid, log_p_grid = np.meshgrid(unique_lat, unique_log_p)

O3_grid = O3.reshape(Np, Nlat, order="F")
C_grid = concentration.reshape(Np, Nlat, order="F")

plt.figure(figsize=(8, 6))
plt.pcolormesh(lat_grid, log_p_grid, O3_grid, shading='auto', cmap='jet')
plt.colorbar(label="Ozone Concentration")
plt.yscale("log")
plt.ylim([min(p),max(p)])
plt.gca().invert_yaxis()
plt.xlabel("Latitude")
plt.ylabel("Pressure (hPa)")
plt.title("Ozone Concentration")

plt.figure(figsize=(8, 6))
plt.pcolormesh(lat_grid, log_p_grid, C_grid, shading='auto', cmap='jet')
plt.colorbar(label="Ozone Concentration")
plt.yscale("log")
plt.ylim([min(p),max(p)])
plt.gca().invert_yaxis()
plt.xlabel("Latitude")
plt.ylabel("Pressure (hPa)")
plt.title("Ozone Concentration")

#plt.show()

## contour plots

# mixing ratio plot
plt.figure(figsize=(10, 6))
#levels = np.logspace(np.log10(O3_range[0]), np.log10(O3_range[1]), num=20)
levels = np.linspace(O3_range[0], O3_range[1], num=100)
contour = plt.contourf(lat_grid, log_p_grid, O3_grid, levels=levels, cmap="jet", norm=mcolors.Normalize())
plt.contour(lat_grid, log_p_grid, O3_grid, levels=np.linspace(O3_range[0], O3_range[1], num=10), colors='black', norm=mcolors.Normalize(), linewidths=1)
cbar = plt.colorbar(contour)
cbar.set_label("Ozone Mixing Ratio [mol/mol]")

plt.gca().invert_yaxis()  # inverted to simulate the actual atmosphere(high pressure near ground)
plt.yscale("log")
plt.xlabel("Latitude")
plt.ylabel("Pressure [Pa]")
#plt.title("Ozone Mixing Ratio in Latitude-Pressure Grid (2019)")

# concentration plot
plt.figure(figsize=(10, 6))
#levels = np.logspace(np.log10(C_range[0]), np.log10(C_range[1]), num=20)
levels = np.linspace(C_range[0],C_range[1], num=100)
contour = plt.contourf(lat_grid, log_p_grid, C_grid, levels=levels, cmap="jet", norm=mcolors.Normalize())
plt.contour(lat_grid, log_p_grid, C_grid, levels=np.linspace(C_range[0], C_range[1], num=10), colors='black', norm=mcolors.Normalize(), linewidths=1)
#contour = plt.contourf(lat_grid, 10**plev_grid, O3_grid, levels=100, cmap="viridis", norm=mcolors.LogNorm())
cbar = plt.colorbar(contour)
cbar.set_label("Ozone Concentration [mol/m3]")

plt.gca().invert_yaxis()
plt.yscale("log")
plt.xlabel("Latitude")
plt.ylabel("Pressure [Pa]")
#plt.title("Ozone Concentration in Latitude-Pressure Grid (2019)")

for i in range(len(O3)):
    if math.isnan(O3[i]):
        O3[i]=0
        concentration[i]=0

p = p[:Nlat * Np]

# find high ozone concentrations, old method, changed after the faq said one mean value is enough
threshold = np.percentile(O3, 95)  # top 5%
threshold2 = np.percentile(concentration, 95)  # top 5%

ozone_layer = p[O3 >= threshold]
ozone_layer2 = p[concentration >= threshold2]
O3_lat = np.reshape(O3, (Nlat, Np), order='C')
C_lat = np.reshape(concentration, (Nlat, Np), order='C')
#print(C_lat.shape)
p_unique = np.unique(p)

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

# ISA calculations inverted to find altitude from given pressure
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
                return h_start + (R_air * T / g0) * np.log(p_current / p)

        T = T_next
        p_current = p_next
    return h_end


altitudes = np.array([pressure_to_altitude(p) for p in ozone_layer])
altitudes2 = np.array([pressure_to_altitude(p) for p in ozone_layer2])
print(f"Ozone layer pressure range based on mixing ratio: {np.min(ozone_layer):.2f} Pa - {np.max(ozone_layer):.2f} Pa")
print(f"Ozone layer pressure range based on concentration: {np.min(ozone_layer2):.2f} Pa - {np.max(ozone_layer2):.2f} Pa")
print(" ")
print(f"Ozone layer altitude range based on mixing ratio: {np.min(altitudes)/1000:.2f} km - {np.max(altitudes)/1000:.2f} km")
print(f"Ozone layer altitude range based on concentration: {np.min(altitudes2)/1000:.2f} km - {np.max(altitudes2)/1000:.2f} km")
print(" ")

# find the pressure level with max O3 per latitude
max_O3_indices = np.argmax(O3_lat, axis=1)
max_O3_pressures = p_unique[max_O3_indices]

max_O3_indices_C = np.argmax(C_lat, axis=1)
max_O3_pressures_C = p_unique[max_O3_indices_C]

max_O3_altitudes = np.array([pressure_to_altitude(pressure) for pressure in max_O3_pressures])
max_O3_altitudes_C = np.array([pressure_to_altitude(pressure) for pressure in max_O3_pressures_C])

mean_ozone_altitude = np.mean(max_O3_altitudes)
mean_ozone_altitude_C = np.mean(max_O3_altitudes_C)

print(f"Global mean ozone layer altitude (max O3 r approach): {mean_ozone_altitude / 1000:.2f} km")
print(" ")
print(f"Global mean ozone layer altitude (max O3 C approach): {mean_ozone_altitude_C / 1000:.2f} km")
print(" ")

print(np.mean(ozone_layer),np.mean(ozone_layer2))


plt.figure(figsize=(8, 6))
plt.pcolormesh(lat_grid, log_p_grid, C_lat.T, shading='auto', cmap='jet')  # Plot the structured ozone concentration
plt.plot(lat_grid[0], p_unique[max_O3_indices_C], linestyle='--', color='dimgray', label="Max O3 per latitude")
plt.colorbar()
plt.yscale('log')
plt.xlabel("Latitude Index")
plt.ylabel("Pressure Index")
plt.ylim([min(p),max(p)])
plt.gca().invert_yaxis()
plt.legend()


plt.figure(figsize=(8, 6))
plt.pcolormesh(lat_grid, log_p_grid, O3_lat.T, shading='auto', cmap='jet')  # Plot the structured ozone concentration
plt.plot(lat_grid[0], p_unique[max_O3_indices], linestyle='--', color='dimgray', label="Max O3 per latitude")
plt.colorbar()
plt.yscale('log')
plt.xlabel("Latitude Index")
plt.ylabel("Pressure Index")
plt.ylim([min(p),max(p)])
plt.gca().invert_yaxis()
plt.legend()
#plt.show()

## question 2

#initial mixing ratios
O3_0 = 0 # nmol/mol
NO_0 = np.logspace(-2, 4, num=1000, base=10)
NO2_0 = np.logspace(-2, 4, num=1000, base=10)

alpha = 10 #nmol/mol

# steady state ozone function
def O3_ss(NO_0,NO2_0,O3_0,zero):
    O3_0 = 0
    if isinstance(NO_0, np.ndarray):
        O3 = np.zeros(len(NO_0))
    else:
        O3 = 0

    NO_0 = NO_0.copy()
    NO2_0 = NO2_0.copy()

    if isinstance(NO_0, np.ndarray):
        if zero==1:
            NO_0[:] = 0

        elif zero==2:
            NO2_0[:] = 0
    else:
        if zero == 1:
            NO_0 = 0

        elif zero == 2:
            NO2_0 = 0

    if isinstance(NO_0, np.ndarray):
        for i in range(len(NO_0)):
            O3[i] = -0.5 * (NO_0[i] - O3_0 + alpha) + 0.5 * (
                        (NO_0[i] - O3_0 + alpha) ** 2 + 4 * alpha * (NO2_0[i] + O3_0)) ** 0.5
    else:
        O3 = -0.5 * (NO_0 - O3_0 + alpha) + 0.5 * ((NO_0 - O3_0 + alpha) ** 2 + 4 * alpha * (NO2_0 + O3_0)) ** 0.5

    return O3

#generate O3 levels for the three cases
O3_ss1 = O3_ss(NO_0,NO2_0,O3_0,2)
O3_ss2 = O3_ss(NO_0,NO2_0,O3_0,1)
O3_ss12 = O3_ss(NO_0,NO2_0,O3_0,0)

# extract data from files
station_mapping = {
    "NL00418": "rotterdam",
    "NL00131": "vredepeel",
    "DEBY109": "andechs"
}

folder_path = "O3_NO_NO2_measurements"

class MeasurementData:
    def __init__(self, **datasets):
        for molecule, df in datasets.items():
            setattr(self, molecule, df.iloc[:, 0])

stations = {}

for file in os.listdir(folder_path):
    if file.endswith(".csv"):
        file_path = os.path.join(folder_path, file)
        df = pd.read_csv(file_path, parse_dates=["time"], index_col="time")

        parts = file.replace(".csv", "").split("_")


        _, station_id, molecule = parts
        station_name = station_mapping.get(station_id)

        if station_name not in stations:
            stations[station_name] = {}

        stations[station_name][molecule] = df

for station in stations:
    stations[station] = MeasurementData(**stations[station])

rotterdam = stations["rotterdam"]
vredepeel = stations["vredepeel"]
andechs = stations["andechs"]

locations = [rotterdam, vredepeel, andechs]

#conversion factors
F_O3 = 2 #kg/m^3
F_NO = 1.3 #kg/m^3
F_NO2 = 1.9 #kg/m^3

# find annual mean
for station_name, station_data in stations.items():
    for molecule, df in station_data.__dict__.items():
        print(f"Station: {station_name}, {molecule}: Annual mean: {np.mean(df):.3f} µg/m³")
    print(" ")

# convert concentrations to annual mean initial mixing ratios
for place in locations:
    place.r_O3_0 = np.mean(place.O3) / F_O3
    place.r_NO_0 = np.mean(place.NO) / F_NO
    place.r_NO2_0 = np.mean(place.NO2) / F_NO2

rotterdam.name = "Rotterdam"
vredepeel.name = "Vredepeel"
andechs.name = "Andechs"

# find steady state ozone mixing ratios
for place in locations:
    place.O3_ss_mean = O3_ss(place.r_NO_0,place.r_NO2_0,place.r_O3_0,0)
    print(f"{place.name}:")
    print(f"Annual Mean Steady State Ozone Mixing Ratio: {place.O3_ss_mean:.3f} [nmol/mol]")
    print(f"Annual Mean Steady State Ozone Concentration: {place.O3_ss_mean*F_O3:.3f} [µg/m³]")
    print(f"Percentage difference between observed mean and mean steady state concentration {(place.O3_ss_mean*F_O3-np.mean(place.O3))/np.mean(place.O3) * 100:.2f}% ")
    print(" ")


#plot all three cases
plt.figure(figsize=(10,6))
plt.plot(NO_0,O3_ss1,label="Initial mixing ratio of $NO$ only")
plt.plot(NO_0,O3_ss2,label="Initial mixing ratio of $NO_2$ only")
plt.plot(NO_0,O3_ss12,label="Same initial mixing ratios of $NO$ and $NO_2$")
#plt.scatter([rotterdam.r_NO_0,vredepeel.r_NO_0,andechs.r_NO_0],[rotterdam.r_O3_0,vredepeel.r_O3_0,andechs.r_O3_0],label="Mean observed values", color='c')
#plt.scatter([rotterdam.r_NO2_0,vredepeel.r_NO2_0,andechs.r_NO2_0],[rotterdam.O3_ss_mean,vredepeel.O3_ss_mean,andechs.O3_ss_mean],label="Mean observed SS values NO2 on x-axis", color='r')
#plt.scatter([rotterdam.r_NO_0,vredepeel.r_NO_0,andechs.r_NO_0],[rotterdam.O3_ss_mean,vredepeel.O3_ss_mean,andechs.O3_ss_mean],label="Mean observed SS values NO on x-axis", color='m')
plt.legend()
plt.xlabel("$NO_x$ mixing ratios [nmol/mol]")
plt.ylabel("Steady state $O_3$ mixing ratio [nmol/mol]")
plt.xscale("log")
plt.yscale("linear")
plt.ylim([-5, 40])
plt.xlim([10**-2,10**4])
plt.grid()

plt.show()


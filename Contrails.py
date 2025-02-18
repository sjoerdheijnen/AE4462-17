import numpy as np
import matplotlib.pyplot as plt

#Constants and parameters
EI_kerosene = 1.25 #kg/kg
EI_hydrogen = 8.94 #kg/kg
epsilon = 0.622
LHV_kerosene = 43.2 #MJ/kg
LHV_hydrogen = 120.0 #MJ/kg
cp = 1004 #J/kg/K

#Coefficients [e_l, e_i]
b = [6.5459673, 4.1635019]

a = [[-0.58002206 * 10**4, -0.56745359 * 10**4],
     [1.3914993, 6.3925247],
     [-0.48640239 * 10**(-1), -0.96778430 * 10**(-2)],
     [0.41764768 * 10**(-4), 0.62215701 * 10**(-6)],
     [-0.14452093 * 10**(-7), 0.20747825 * 10**(-8)],
     [0, -0.94840240 * 10**(-12)]]

#Saturation vapor pressure plotting
def saturation_vapor_pressure(T,state):

    index = 1 if state == 'ice' else 0

    ln_ep = b[index] * np.log(float(T)) + sum(a[j + 1][index] * float(T)**j for j in range(-1, 4))
    ep = np.exp(ln_ep)

    return ln_ep, ep

T = np.array(range(174,291))

ep_liquid_lst = []
ln_ep_liquid_lst = []
ep_ice_lst = []
ln_ep_ice_lst = []

for i in T:
    ln_ep_liquid, ep_liquid = saturation_vapor_pressure(i,'liquid')
    ln_ep_ice, ep_ice = saturation_vapor_pressure(i,'ice')
    ep_liquid_lst.append(ep_liquid)
    ln_ep_liquid_lst.append(ln_ep_liquid)
    ep_ice_lst.append(ep_ice)
    ln_ep_ice_lst.append(ln_ep_ice)

#flight parameters
P = 22000 #Pa
t = 225 #K
humid = 1.1 #wrt to ice
eta = 0.3 #engine efficiency

def mixing_line(T,fuel):
    EI = EI_kerosene if fuel == 'kerosene' else EI_hydrogen
    LHV = LHV_kerosene if fuel == 'kerosene' else LHV_hydrogen

    G = (cp*P)/epsilon * EI/((1-eta)*(LHV*(10**6)))
    return G

G_slope = mixing_line(t,'kerosene')
T_mix = []
for i in T:
    if i>=t :
        T_mix.append(i)
T_mix = np.array(T_mix)
e_mixing = G_slope * (T_mix - t) + humid * saturation_vapor_pressure(t, 'ice')[1]

offset = P-t*G_slope
print(offset)
y = G_slope*T+offset

# print(G_slope)
# print(e_mixing)
# print(y)

plt.figure()
plt.plot(T, ep_liquid_lst, label="Liquid", color='b')
plt.plot(T, ep_ice_lst, label="Ice", color='c')
plt.plot(T_mix, e_mixing, label='Mixing Line', linestyle='--', color='r')
plt.scatter(t,humid * saturation_vapor_pressure(t, 'ice')[1], label='Atmospheric Condition',color='r')
plt.axvline(x=t, color='gray', linestyle=':', label=f'Ambient Temperature = {t} K')

plt.xlabel("Temperature (K)")
plt.ylabel("$e_p$ [Pa]")
plt.title("Saturation Vapor Pressure")
plt.legend()
plt.grid(True)
#plt.yscale('log')
plt.xlim([273-60,273-20])
plt.ylim([0,50])
plt.show()

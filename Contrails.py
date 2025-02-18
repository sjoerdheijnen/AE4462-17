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

def mixing_line(T,eta,P,t,humid,fuel,state):
    EI = EI_kerosene if fuel == 'kerosene' else EI_hydrogen
    LHV = LHV_kerosene if fuel == 'kerosene' else LHV_hydrogen

    G = (cp*P)/epsilon * EI/((1-eta)*(LHV*(10**6)))

    T_mix = []
    for i in T:
        if i>=t :
            T_mix.append(i)
    T_mix = np.array(T_mix)

    e_mixing = G * (T_mix - t) + humid * saturation_vapor_pressure(t, state)[1]

    return e_mixing, T_mix

#first plot mixing line conditions
e_mixing1b, T_mix1b = mixing_line(T,0.3,22000,225,1.1,'kerosene','ice')
e_mixing1c, T_mix1c = mixing_line(T,0.4,22000,225,1.1,'kerosene','ice')
e_mixing1d, T_mix1d = mixing_line(T,0.4,22000,225,1.1,'hydrogen','ice')

plt.figure()
plt.plot(T, ep_liquid_lst, label="Liquid", color='b')
plt.plot(T, ep_ice_lst, label="Ice", color='c')
plt.plot(T_mix1b, e_mixing1b, label='Mixing Line aircraft 1b', linestyle='--', color='r')
plt.plot(T_mix1c, e_mixing1c, label='Mixing Line aircraft 1c', linestyle='--', color='orange')
plt.plot(T_mix1d, e_mixing1d, label='Mixing Line aircraft 1d', linestyle='--', color='yellow')

plt.scatter(225,1.1 * saturation_vapor_pressure(225, 'ice')[1], marker='.', label='Atmospheric Condition',color='black', zorder=7)
plt.axvline(x=t, color='gray', linestyle=':', label=f'Ambient Temperature = {t} K')

plt.xlabel("Temperature (K)")
plt.ylabel("$e_p$ [Pa]")
plt.title("Saturation Vapor Pressure")
plt.legend()
plt.grid(True)
plt.xlim([273-60,273-20])
plt.ylim([0,50])
plt.show()

e_mixing1f1, T_mix1f1 = mixing_line(T,0.3,25000,230,0.6,'kerosene','liquid')
e_mixing1f2, T_mix1f2 = mixing_line(T,0.4,25000,230,0.6,'kerosene','liquid')
e_mixing1f3, T_mix1f3 = mixing_line(T,0.4,25000,230,0.6,'hydrogen','liquid')

plt.figure()
plt.plot(T, ep_liquid_lst, label="Liquid", color='b')
plt.plot(T, ep_ice_lst, label="Ice", color='c')
plt.plot(T_mix1f1, e_mixing1f1, label='Mixing Line aircraft 1f-b', linestyle='--', color='r')
plt.plot(T_mix1f2, e_mixing1f2, label='Mixing Line aircraft 1f-c', linestyle='--', color='orange')
plt.plot(T_mix1f3, e_mixing1f3, label='Mixing Line aircraft 1f-d', linestyle='--', color='yellow')

plt.scatter(230,0.6 * saturation_vapor_pressure(230, 'liquid')[1], marker='.', label='Atmospheric Condition',color='black', zorder=7)
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

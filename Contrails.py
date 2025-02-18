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

    ln_ep = b[index] * np.log10(T)+sum(a[j + 1][index] * float(T)**j for j in range(-1, 4))
    ep = np.exp(ln_ep)
    return ln_ep, ep

T = np.array(range(150,291))

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


plt.figure()
plt.plot(T, ln_ep_liquid_lst, label="Liquid", color='b')
plt.plot(T, ln_ep_ice_lst, label="Ice", color='r')
plt.xlabel("Temperature (K)")
plt.ylabel("ln($e_p$) [ln(Pa)]")  # Natural logarithm of pressure
plt.title("Saturation Vapor Pressure")
plt.legend()  # Add legend
plt.grid(True)  # Enable grid
plt.show()

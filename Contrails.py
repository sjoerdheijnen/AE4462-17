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

a = [[-0.58002206 * 10**4, -0.56745359 * 10**4]
     [1.3914993, 6.3925247]
     [-0.48640239 * 10**(-1), -0.96778430 * 10**(-2)]
     [0.41764768 * 10**(-4), 0.62215701 * 10**(-6)]
     [-0.14452093 * 10**(-7), 0.20747825 * 10**(-8)]
     [0, -0.94840240 * 10**(-12)]]

#Saturation vapor pressure plotting
def saturation_vapor_pressure(T,state):
    index = 0
    if state == 'ice':
        index = 1
    ln_ep = b[index] * np.log10(T)+np.sum()
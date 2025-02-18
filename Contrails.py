import numpy as np
import matplotlib.pyplot as plt
## question 1

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

#Saturation vapor pressure function that calculates the curves
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


#mixing line function that takes atmospheric conditions and constructs the mixing line
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
plt.plot(T_mix1c, e_mixing1c, label='Mixing Line aircraft 1c', linestyle='--', color='limegreen')
plt.plot(T_mix1d, e_mixing1d, label='Mixing Line aircraft 1d', linestyle='--', color='magenta')

plt.scatter(225,1.1 * saturation_vapor_pressure(225, 'ice')[1], marker='.', label='Atmospheric Condition',color='black', zorder=7)
plt.axvline(x=225, color='gray', linestyle=':', label=f'Ambient Temperature = {225} K')

plt.xlabel("Temperature (K)")
plt.ylabel("$e_p$ [Pa]")
plt.title("Saturation Vapor Pressure")
plt.legend()
plt.grid(True)
plt.xlim([273-60,273-20])
plt.ylim([0,50])
plt.show()

#second plot mixing line conditions
e_mixing1f1, T_mix1f1 = mixing_line(T,0.3,25000,230,0.6,'kerosene','liquid')
e_mixing1f2, T_mix1f2 = mixing_line(T,0.4,25000,230,0.6,'kerosene','liquid')
e_mixing1f3, T_mix1f3 = mixing_line(T,0.4,25000,230,0.6,'hydrogen','liquid')

plt.figure()
plt.plot(T, ep_liquid_lst, label="Liquid", color='b')
plt.plot(T, ep_ice_lst, label="Ice", color='c')
plt.plot(T_mix1f1, e_mixing1f1, label='Mixing Line aircraft 1f-b', linestyle='--', color='r')
plt.plot(T_mix1f2, e_mixing1f2, label='Mixing Line aircraft 1f-c', linestyle='--', color='limegreen')
plt.plot(T_mix1f3, e_mixing1f3, label='Mixing Line aircraft 1f-d', linestyle='--', color='magenta')

plt.scatter(230,0.6 * saturation_vapor_pressure(230, 'liquid')[1], marker='.', label='Atmospheric Condition',color='black', zorder=7)
plt.axvline(x=230, color='gray', linestyle=':', label=f'Ambient Temperature = {230} K')

plt.xlabel("Temperature (K)")
plt.ylabel("$e_p$ [Pa]")
plt.title("Saturation Vapor Pressure")
plt.legend()
plt.grid(True)
plt.xlim([273-60,273-20])
plt.ylim([0,50])
plt.show()


## question 2
tau = 0.1
sigma = 5.670374419 * 10**(-8) #Wm^-2K^-4
T0 = 288.15 #K
tf = 3*3600 #s
W = 1000 #m

def beerlambert(tau,F_in):
    F_out = F_in * np.exp(-tau)
    return F_out

def RFCLW(tau,F_in):
    RFC_LW = F_in * (1-np.exp(-tau))
    return RFC_LW

F_in = sigma * T0**4
print(f"The ingoing Flux is {round(F_in,2)} W/m^2")

RFC_LW = RFCLW(tau,F_in)
print(f"The change in Flux due to the contrail is {round(RFC_LW,2)} W/m^2")

def EF(RFC_LW,W,l,tf):
    EF = RFC_LW*W*l*tf
    return EF
EF1 = EF(RFC_LW,W,1,tf)
print(f"The energy forcing due to a meter of contrail is {round(EF1/10**6,2)} MJ/m")

L_flown = 50 * 10**12 #m
L_contrail = L_flown * 0.1 #m

print(f"Total length of persistent contrails in a year: {L_contrail/10**12} billion km")

EF_total = EF(RFC_LW,W,L_contrail,tf)

print(f"Total energy forcing due to contrails in a year is {round(EF_total/10**18,2)} billion GJ")

R_earth = 6371000 #m
A_earth = 4*np.pi*R_earth**2 #m2
T_year = 365*24*3600 #s

RF_total = (EF_total/T_year) / A_earth #W/m2

print(f"Radiative forcing due to contrails in a year is {round(RF_total*1000,1)} mW/m^2")

def RF(tau,T0,W,tf,L_flown,eta_contrail):
    F_in = sigma * T0 ** 4
    RFC_LW = RFCLW(tau, F_in)
    L_contrail = L_flown * eta_contrail
    EF_total = EF(RFC_LW, W, L_contrail, tf)
    RF_total = (EF_total / T_year) / A_earth  # W/m2
    return RF_total

#change these parameters to see the change in RF
tau_test = 0.1
T0_test = 288.15 #K
tf_test = 3 #h
W_test = 1 #km
L_flown_test = 50 #billion km
eta_contrail_test = 0.1
RF_test = RF(tau_test,T0_test,W_test * 10**3,tf_test * 3600,L_flown_test * 10**12,eta_contrail_test)

print(f"New radiative forcing: {round(RF_test*1000,1)} mW/m^2 with a relative difference of {round((RF_test-RF_total)/RF_total*100,2)} %")
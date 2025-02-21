import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
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
    offset = -t * G + humid * saturation_vapor_pressure(t, state)[1]

    return e_mixing, T_mix, offset, G

#first plot mixing line conditions
e_mixing1b, T_mix1b, offset1b, G1b = mixing_line(T,0.3,22000,225,1.1,'kerosene','ice')
e_mixing1c, T_mix1c, offset1c, G1c = mixing_line(T,0.4,22000,225,1.1,'kerosene','ice')
e_mixing1d, T_mix1d, offset1d, G1d = mixing_line(T,0.4,22000,225,1.1,'hydrogen','ice')

#second plot mixing line conditions
e_mixing1f1, T_mix1f1, offset1f1, G1f1 = mixing_line(T,0.3,25000,230,0.6,'kerosene','liquid')
e_mixing1f2, T_mix1f2, offset1f2, G1f2 = mixing_line(T,0.4,25000,230,0.6,'kerosene','liquid')
e_mixing1f3, T_mix1f3, offset1f3, G1f3 = mixing_line(T,0.4,25000,230,0.6,'hydrogen','liquid')

#Constructs the threshold line for a given mixing line
def threshold_line(G, offset, e_mixing, state):
    e_atm = [ep_ice_lst, ep_liquid_lst]

    difference = np.abs(np.array(np.gradient(e_atm[1])) - G)
    index = difference.argmin()
    T_tangent = T[index]

    e_tangent = saturation_vapor_pressure(T_tangent,state)[1]
    offset = e_tangent - T_tangent * G
    e_threshold = T * G + offset

    difference2 = e_atm[0] - e_threshold
    sign_change = np.where(np.diff(np.sign(difference2)))[0]

    index2 = sign_change[0]

    T_intersect = T[index2]

    return e_threshold, T_tangent, T_intersect


#threshold line construction for each aircraft
e_threshold1b, T_tangent1b, T_intersect1b = threshold_line(G1b, offset1b, e_mixing1b, 'liquid')
e_threshold1c, T_tangent1c, T_intersect1c = threshold_line(G1c, offset1c, e_mixing1c, 'liquid')
e_threshold1d, T_tangent1d, T_intersect1d = threshold_line(G1d, offset1d, e_mixing1d, 'liquid')
e_threshold1f1, T_tangent1f1, T_intersect1f1 = threshold_line(G1f1, offset1f1, e_mixing1f1, 'liquid')
e_threshold1f2, T_tangent1f2, T_intersect1f2 = threshold_line(G1f2, offset1f2, e_mixing1f2, 'liquid')
e_threshold1f3, T_tangent1f3, T_intersect1f3 = threshold_line(G1f3, offset1f3, e_mixing1f3, 'liquid')

#plots first set of mixing lines
plt.figure()
plt.plot(T, ep_liquid_lst, label="Liquid saturation curve", color='b')
plt.plot(T, ep_ice_lst, label="Ice saturation curve", color='c')
plt.plot(T_mix1b, e_mixing1b, label='Mixing Line aircraft 1b', linestyle='--', color='r')
plt.plot(T_mix1c, e_mixing1c, label='Mixing Line aircraft 1c', linestyle='--', color='limegreen')
plt.plot(T_mix1d, e_mixing1d, label='Mixing Line aircraft 1d', linestyle='--', color='magenta')

plt.scatter(225,1.1 * saturation_vapor_pressure(225, 'ice')[1], marker='.', label='Atmospheric Condition',color='black', zorder=7)
plt.axvline(x=225, color='gray', linestyle=':', label=f'Ambient Temperature = {225} K')

plt.xlabel("Temperature (K)")
plt.ylabel("$e_p$ [Pa]")
#plt.title("Saturation Vapor Pressure")
plt.legend()
plt.grid(True)
plt.xlim([273-60,273-20])
plt.ylim([0,50])
#plt.show()

#plots second set of mixing lines
plt.figure()
plt.plot(T, ep_liquid_lst, label="Liquid saturation curve", color='b')
plt.plot(T, ep_ice_lst, label="Ice saturation curve", color='c')
plt.plot(T_mix1f1, e_mixing1f1, label='Mixing Line aircraft 1f-b', linestyle='--', color='r')
plt.plot(T_mix1f2, e_mixing1f2, label='Mixing Line aircraft 1f-c', linestyle='--', color='limegreen')
plt.plot(T_mix1f3, e_mixing1f3, label='Mixing Line aircraft 1f-d', linestyle='--', color='magenta')

plt.scatter(230,0.6 * saturation_vapor_pressure(230, 'liquid')[1], marker='.', label='Atmospheric Condition',color='black', zorder=7)
plt.axvline(x=230, color='gray', linestyle=':', label=f'Ambient Temperature = {230} K')

plt.xlabel("Temperature (K)")
plt.ylabel("$e_p$ [Pa]")
#plt.title("Saturation Vapor Pressure")
plt.legend()
plt.grid(True)
plt.xlim([273-60,273-20])
plt.ylim([0,50])
#plt.show()


# plots threshold line and contrail area for given aircraft
def plot_aircraft(ac):
    plt.figure()
    plt.plot(T, ep_liquid_lst, label="Liquid saturation curve", color='b')
    plt.plot(T, ep_ice_lst, label="Ice saturation curve", color='c')

    if ac == '1b':
        plt.plot(T_mix1b, e_mixing1b, label='Mixing Line aircraft 1b', linestyle='--', color='r')
        plt.plot(T, e_threshold1b, color='green', zorder=5, label='Threshold line')
        plt.fill_between(T, ep_liquid_lst, e_threshold1b, where=(T <= T_tangent1b),
                         interpolate=True, color='yellow', alpha=0.7, label='Persistent Contrails')
        plt.fill_between(T, ep_ice_lst, e_threshold1b, where=(T <= T_intersect1b),
                         interpolate=True, color='orange', alpha=1, label='Dissipative Contrails')

        plt.scatter(225, 1.1 * saturation_vapor_pressure(225, 'ice')[1], marker='.', label='Atmospheric Condition',
                    color='black', zorder=7)
        plt.axvline(x=225, color='gray', linestyle=':', label=f'Ambient Temperature = {225} K')
    elif ac == '1c':
        plt.plot(T_mix1c, e_mixing1c, label='Mixing Line aircraft 1c', linestyle='--', color='limegreen')
        plt.plot(T, e_threshold1c, color='green', zorder=5, label='Threshold line')
        plt.fill_between(T, ep_liquid_lst, e_threshold1c, where=(T <= T_tangent1c),
                         interpolate=True, color='yellow', alpha=0.7, label='Persistent Contrails')
        plt.fill_between(T, ep_ice_lst, e_threshold1c, where=(T <= T_intersect1c),
                         interpolate=True, color='orange', alpha=1, label='Dissipative Contrails')

        plt.scatter(225, 1.1 * saturation_vapor_pressure(225, 'ice')[1], marker='.', label='Atmospheric Condition',
                    color='black', zorder=7)
        plt.axvline(x=225, color='gray', linestyle=':', label=f'Ambient Temperature = {225} K')
    elif ac == '1d':
        plt.plot(T_mix1d, e_mixing1d, label='Mixing Line aircraft 1d', linestyle='--', color='magenta')
        plt.plot(T, e_threshold1d, color='green', zorder=5, label='Threshold line')
        plt.fill_between(T, ep_liquid_lst, e_threshold1d, where=(T <= T_tangent1d),
                         interpolate=True, color='yellow', alpha=0.7, label='Persistent Contrails')
        plt.fill_between(T, ep_ice_lst, e_threshold1d, where=(T <= T_intersect1d),
                         interpolate=True, color='orange', alpha=1, label='Dissipative Contrails')

        plt.scatter(225, 1.1 * saturation_vapor_pressure(225, 'ice')[1], marker='.', label='Atmospheric Condition',
                    color='black', zorder=7)
        plt.axvline(x=225, color='gray', linestyle=':', label=f'Ambient Temperature = {225} K')
    elif ac == '1f-b':
        plt.plot(T_mix1f1, e_mixing1f1, label='Mixing Line aircraft 1f-b', linestyle='--', color='r')
        plt.plot(T, e_threshold1f1, color='green', zorder=5, label='Threshold line')
        plt.fill_between(T, ep_liquid_lst, e_threshold1f1, where=(T <= T_tangent1f1),
                         interpolate=True, color='yellow', alpha=0.7, label='Persistent Contrails')
        plt.fill_between(T, ep_ice_lst, e_threshold1f1, where=(T <= T_intersect1f1),
                         interpolate=True, color='orange', alpha=1, label='Dissipative Contrails')

        plt.scatter(230, 0.6 * saturation_vapor_pressure(230, 'liquid')[1], marker='.', label='Atmospheric Condition',
                    color='black', zorder=7)
        plt.axvline(x=230, color='gray', linestyle=':', label=f'Ambient Temperature = {230} K')
    elif ac == '1f-c':
        plt.plot(T_mix1f2, e_mixing1f2, label='Mixing Line aircraft 1f-c', linestyle='--', color='limegreen')
        plt.plot(T, e_threshold1f2, color='green', zorder=5, label='Threshold line')
        plt.fill_between(T, ep_liquid_lst, e_threshold1f2, where=(T <= T_tangent1f2),
                         interpolate=True, color='yellow', alpha=0.7, label='Persistent Contrails')
        plt.fill_between(T, ep_ice_lst, e_threshold1f2, where=(T <= T_intersect1f2),
                         interpolate=True, color='orange', alpha=1, label='Dissipative Contrails')

        plt.scatter(230, 0.6 * saturation_vapor_pressure(230, 'liquid')[1], marker='.', label='Atmospheric Condition',
                    color='black', zorder=7)
        plt.axvline(x=230, color='gray', linestyle=':', label=f'Ambient Temperature = {230} K')
    elif ac == '1f-d':
        plt.plot(T_mix1f3, e_mixing1f3, label='Mixing Line aircraft 1d', linestyle='--', color='magenta')
        plt.plot(T, e_threshold1f3, color='green', zorder=5, label='Threshold line')
        plt.fill_between(T, ep_liquid_lst, e_threshold1f3, where=(T <= T_tangent1f3),
                         interpolate=True, color='yellow', alpha=0.7, label='Persistent Contrails')
        plt.fill_between(T, ep_ice_lst, e_threshold1f3, where=(T <= T_intersect1f3),
                         interpolate=True, color='orange', alpha=1, label='Dissipative Contrails')

        plt.scatter(230, 0.6 * saturation_vapor_pressure(230, 'liquid')[1], marker='.', label='Atmospheric Condition',
                    color='black', zorder=7)
        plt.axvline(x=230, color='gray', linestyle=':', label=f'Ambient Temperature = {230} K')

    plt.xlabel("Temperature (K)")
    plt.ylabel("$e_p$ [Pa]")
    # plt.title("Saturation Vapor Pressure")
    plt.legend()
    plt.grid(True)
    plt.xlim([273 - 60, 273 - 20])
    plt.ylim([0, 50])
    #plt.show()
    return

plot_aircraft('1b')
plot_aircraft('1c')
plot_aircraft('1d')
plot_aircraft('1f-b')
plot_aircraft('1f-c')
plot_aircraft('1f-d')

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

#change these parameters to see the change in RF, meant to identify which input changes the RF is sensitive too
tau_test = 0.1
T0_test = 288.15 #K
tf_test = 3 #h
W_test = 1 #km
L_flown_test = 50 #billion km
eta_contrail_test = 0.1
RF_test = RF(tau_test,T0_test,W_test * 10**3,tf_test * 3600,L_flown_test * 10**12,eta_contrail_test)

print(f"New radiative forcing: {round(RF_test*1000,1)} mW/m^2 with a relative difference of {round((RF_test-RF_total)/RF_total*100,2)} %")

plt.show()
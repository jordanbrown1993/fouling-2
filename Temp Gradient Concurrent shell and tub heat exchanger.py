# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt 
from numpy import log as ln
from math import e




#Shell & Tube Heat exchanger with fouling 

L = 60 # m, length of pipe
r1 = 0.05 # m, tube radii
r2 = 0.5 # m, inner shell radius
N_t = 10 # number of tubes 
n = 100 # number of increments/nodes

mu_1 = (1000*e**(-52.84+(3700/T1)+(5.866*ln(T1))))/1000 # dynamic viscosity in kg/ms 

m1 = 3 # kg/s, mass flow rate of fluid 1
Cp1 = 4180 # J/Kg*K, specific heat capacity of fluid 1 (water)
rho1 = 1250 # kg/m(^3), density of fluid 1 (crude oil) 
Ac1 = N_t*pi*r1**2 # m^2, cross-sectional area of tube bundle
Ac2 = (pi*r2**2)-Ac1# # cross sectional flow area of cooling fluid

mu_2 = e**(-6.57+(178.1/T2)+(-0.9524*ln(T2)))

m2 = 5 # kg/s, mass flow rate of fluid 2
Cp2 = 4180 # J/kg*K, specific heat capacity of fluid 2 (water)
rho2 = 1000 # kg/m(^3), density of fluid 2 (water)

pi = 3.141592653589793
t_w = 0.002 # wall thickness, m


T1i = 400 # K, inlet temperature of fluid through tubes 
T2i = 800 # K, surface temperature of shell
T0 = 300 # K, initial temperature of fluid throughout the pipe

 
kw = 18 # W/m*K, thermal conductivity of pipe shell side
kfi = 2 # W/m*K, thermal conductivity of fouling layer 
Rflow = r1-dth # Flow radius of fluid 
tf = 10**(-7) # fouling layer thickness 
Aflow = pi*Rflow**2 # Flow area of crude 
#u = (m1/(rho1*Aflow)) # tube velocity
#hi = (k*Nu)/(2*Rflow) # tube side heat transfer co-efficient
#Nu = 0.0023*(Re**0.8)*(Pr**0.4) # Dittus Boelter Correlation for Nusselt number
#Re = (rho1*u*2*r1)/mu_1) # Reynolds number
#Cf = 0.0035+(0.264*(Re**(-0.42)) # fanning friction factor
#dPdt = 4*Cf*((rho1*u**2)/(4*Rflow) # Pressure drop 
#Pr = (mu1*Cp1)/kw # Prandtl Number 
#tau = Cf*((u**2)/rho1)/2 # fluid shear stress 
Alpha = 0.94 # Ebert-Panchal Parameter
Gamma = 1.2*(10**(-8)) # Ebert-Panchal Parameter
#Ef = 28 # kJ/mol, activation energy of crude oil 
#Tfilm = T0+0.55*(T1-T0)
R = 8.314 # Universal Gas Constant 
#dRFdt = Alpha*(Re**(-0.6))*(Pr**(-0.33))*e**(-Ef/(R*Tfilm))-(Gamma*tau) # ebert-panchal correlation
#UA = Nt/((1/(2*pi*hi*2*r1*dx)+(1/(ln((2*r1)/(2*(r1+t_w))))/(2*pi*dx*kw))+((ln((2*(r1+t_w))/(2*(Rflow))))/(2*pi*dx*kf)))) # calculation for overal heat transfer co-efficient
dx= L/n # cell width

t_final = 1000 #s, simulation time 
dt= 1 # s, time step

x = np.linspace(dx/2, L-dx/2, n)

T1 = np.ones(n)*T0
T2 = np.ones(n)*T0

dT1dt = np.zeros(n)
dT2dt = np.zeros(n)

t = np.arange(0, t_final, dt)

for j in range(1,len(t)): 
    
    #plt.clf()
    
    u = (m1/(rho1*Aflow)) # tube velocity
    mu_1 = e**(-52.84+(3700/T1)+(5.866*ln(T1))) # dynamic viscosity in kg/ms 
    Re = (rho1*u*2*r1)/(mu_1) # Reynolds number
    Pr = (mu1*Cp1)/kw # Prandtl Number 
    Cf = 0.0035+(0.264*(Re**(-0.42)) # fanning friction factor
    #P = 4*Cf*rho1*u**(2)/(4)*(Rflow) # Pressure drop
    Nu = (0.0023)*(Re**(0.8))*(Pr**(0.4)) # Dittus Boelter Correlation for Nusselt number
    hi = (k*Nu)/(2*Rflow) # tube side heat transfer co-efficient
    UA = Nt/((1/(2*pi*hi*2*r1*dx)+(1/(ln((2*r1)/(2*(r1+t_w))))/(2*pi*dx*kw))+((ln((2*(r1+t_w))/(2*(Rflow))))/(2*pi*dx*kf)))) # calculation for overal heat transfer co-efficient
    tau = Cf*((u**2)/rho1)/2 # fluid shear stress
    Tfilm = T0 + 0.55*(T1-T0)
    dRFdt = Alpha*(Re**(-0.6))*(Pr**(-0.33))*e**(-Ef/(R*Tfilm))-(Gamma*tau) # ebert-panchal correlation
    
    
    dT1dt[1:n] = (m1*Cp1*(T1[0:n-1]-T1[1:n])+UA*dx*(T2[1:n]-T1[1:n]))/(rho1*Cp1*Ac1*dx) # energy balance from node 1 to 100 for fluid 1
    dT1dt[0] = (m1*Cp1*(T1i-T1[0])+UA*dx*(T2[0]-T1[0]))/(rho1*Cp1*Ac1*dx) # boundary condition for fluid 1
    
    dT2dt[1:n] = (m2*Cp2*(T2[0:n-1]-T2[1:n])-UA*dx*(T2[1:n]-T1[1:n]))/(rho2*Cp2*Ac2*dx) # energy balance from node 1 to 100 for fluid 2
    dT2dt[0] = (m2*Cp2*(T2i-T2[0])-UA*dx*(T2[0]-T1[0]))/(rho2*Cp2*Ac2*dx) # boundary condition for fluid 2
    
    T1 = T1 + dT1dt*dt
    T2 = T2 + dT2dt*dt
    
    #plt.figure(1)
    #plt.plot(x, T1, color = 'blue', label='inside')
    #plt.plot(x, T2, color = 'red', label='outside')
    #plt.axis ([0, L, 298, 820])
    #plt.xlabel('Distance (n)')
    #plt.ylabel('Temperature (K)')
    #plt.legend(loc = 'upper right')
    #plt.show()
    #plt.pause(0.005)
              
    
    
   
    
    

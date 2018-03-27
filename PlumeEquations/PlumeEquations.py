# -*- coding: utf-8 -*-
"""
Created on Mon Sep 05 12:09:20 2016

@author: thehnen
"""


#%%

'''
Literature:
[1] Enclosure Fire Dynamics, B. Karlsson, J. G. Quintiere, CRC Press LLC, 2000

'''


#%%

import numpy as np
import csv
import random
import os
import re
import glob
import bisect

#%%

# Mean flame height
# Description:
# L - Mean flame hight (in m)
# relD - relative diameter, assuming circular, planar fire source (in m)
# ERR - Energy release rate (in kW)

# Mean flame height
# Formula 4.3 [1].
def MeanFlameHeight(ERR,relD):
    L = 0.235 *(ERR)**(2.0/5.0) - (1.02 * relD)
    return L
 
'''   
# Test based on example 4.1 [1]. 
# Solution should be 2.96 m.

ERR = 1690.
D = 1.6
L = MeanFlameHeight(ERR,D)
print 'Solution: ' + str(L) + ' m'
'''


#%%

# Ideal (or weak) plume equations
# Description:
# Based on following assumptions:
# All energy released from point source, no radiation, density differences 
# between plume and environment are small and can be neglected, velocity
# and temperature profiles are of similar form and top-hat profiles,
# the entrainment velocity is proportional to the local plume velocity (15%)
# ERR: (Total) Energy release rate in kW, 
# z: Plume height, in m 
# g: Acceleration due to gravity, in m/s**2 
# cp: Specific heat at constant pressure, in kJ/(kg K)
# T_infty: Ambient air temperature, in K
# rho_infty: Ambient air density, in kg/m**3


# Upward gas velocity of the plume at height z
# Formula 4.18 [1].
def IdealPlumeVelocity(ERR,z,g=9.81,cp=1.0,T_infty=293,rho_infty=1.2):
    a = 1.94*((g)/(cp * T_infty * rho_infty))**(1.0/3.0)
    u = a *ERR**(1.0/3.0)*z**(-1.0/3.0)
    return u

# Mass flow of plume at height z
# Formula 4.19 [1].
def IdealPlumeMassFlow(ERR,z,g=9.81,cp=1.0,T_infty=293,rho_infty=1.2):
    a = 0.20 * ((g*rho_infty**2)/(cp*T_infty))**(1.0/3.0)
    m = a *ERR**(1.0/3.0)*z**(5.0/3.0)
    return m

# Temperature difference of the plume at height z
# Formula 4.20 [1].
def IdealPlumeTemperature(ERR,z,g=9.81,cp=1.0,T_infty=293,rho_infty=1.2):
    a = 5.0*((T_infty)/(g * cp**2 * rho_infty**2))**(1.0/3.0)
    T = a *ERR**(2.0/3.0)*z**(-5.0/3.0)
    return T


'''
# Test based on example 4.2 [1]. 
# Mass flow should be ~0.95 kg/s.
# Temperature should be ~366 K. (Temperature difference: ~73 K)

ERR = 70
z = 2
m = IdealPlumeMassFlow(ERR,z)
T = IdealPlumeTemperature(ERR,z)
print 'Mass flow: ' + str(m) + ' kg/s'
print 'Temperature: ' + str(T+293) + ' K'
'''


#%%

# Zukoski plume equations
# Description:
# Based on the idael plume, with a small change of the factor 0.20 to 0.21,
# after experimental measurements. Underestimates the mass flow a bit, but
# overestimates it for smaller burners.

# Mass flow of plume at height z
# Formula 4.19 [1].
def ZukoskiPlumeMassFlow(ERR,z,g=9.81,cp=1.0,T_infty=293,rho_infty=1.2):
    a = 0.21 * ((g*rho_infty**2)/(cp*T_infty))**(1.0/3.0)
    m = a *ERR**(1.0/3.0)*z**(5.0/3.0)
    return m


#%%

# Heskestad (or strong) plume equations
# Description:
# Based on ideal plume equations, modified by experimental data.
# Assumptions: A virtual origin is introduced at hight z0, convective ERR Q0
# is used for some properties, top-hat profiles for temperature and velocity
# are changed to gaussian, T0 and u0 describe centerline properties, 
# Boussinesq approximation removed so that large density differences can
# be taken into account - strong plumes.
# ERR: Energy release rate, in kW
# Qc: Convective ERR, in kW
# z: Plume height, in m 
# z0: Virtual origin height, in m 
# relD: relative fire source diameter, in m
# g: Acceleration due to gravity, in m/s**2 
# cp: Specific heat at constant pressure, in kJ/(kg K)
# T_infty: Ambient air temperature, in K
# rho_infty: Ambient air density, in kg/m**3

# Virtual origin z0
# Formula 4.23 [1]
# Based on fire source diameter relD and total energy released ERR
def HeskestadVirtualOrigin(ERR,relD):
    z0 = 0.083 *(ERR)**(2.0/5.0) - (1.02 * relD)
    return z0

# Convective energy release rate (for completeness)
# Usually between 20 % to 40 % of the total energy are lost due to radiation,
# the residue is responsible for convection.

def ConvectiveERR(ERR, radiation = 0.3):
    Qc = ERR * (1 - radiation)
    return Qc

### Properties below the mean flame height L ###

# Mass flow of plume at height z (below flame)
# Formula 4.27 [1].
def HeskestadPlumeMassFlowBelow(Qc,z,L):
    m = 0.0056 * Qc *(z / L)
    return m
    

### Properties above the mean flame height L ###

# Plume radius (above flame)
# Formula 4.24 [1]
def HeskestadPlumeRadiusAbove(T0, z0, z, T_infty=293,):
    b = 0.12 * (T0/T_infty)**(0.5) * (z - z0)
    return b

# Temperature difference of the plume at height z (above flame)
# Formula 4.25 [1].
def HeskestadPlumeTemperatureAbove(Qc,z,z0,g=9.81,cp=1.0,
                                   T_infty=293,rho_infty=1.2):
    a = 9.1*((T_infty)/(g * cp**2 * rho_infty**2))**(1.0/3.0)
    T = a *Qc**(2.0/3.0)*(z-z0)**(-5.0/3.0)
    return T

# Upward gas velocity of the plume at height z (above flame)
# Formula 4.26 [1].
def HeskestadPlumeVelocityAbove(Qc,z,z0,g=9.81,cp=1.0,
                                T_infty=293,rho_infty=1.2):
    a = 3.4*((g)/(cp * T_infty * rho_infty))**(1.0/3.0)
    u = a *Qc**(1.0/3.0)*(z-z0)**(-1.0/3.0)
    return u

# Mass flow of plume at height z (above flame)
# Formula 4.27 [1].
def HeskestadPlumeMassFlowAbove(Qc,z,z0,g=9.81,cp=1.0,
                                T_infty=293,rho_infty=1.2):
    a = 0.071 * Qc**(1.0/3.0) * (z - (z0))**(5.0/3.0)
    m = a + 1.92*10**(-3.0) * Qc
    return m


'''
# Test based on example 4.3 [1].
# Virtual origin should be ~0.01 m.
# Mean flame height should be ~2.96 m.
# Mass flow at 2.5 m should be  ~5.59 kg/s.
# Mass flow at 6.0 m should be ~17.19 kg/s.
# Temperature should be ~433 K (Temperature difference: ~140 K).

ERR = 1690
relD = 1.6
z1 = 2.5
z2 = 6.0

# Calculate virtual origin 
z0 = HeskestadVirtualOrigin(ERR,relD)

# Mean flame height
L = MeanFlameHeight(ERR,relD)

# Convective ERR
Qc = ConvectiveERR(ERR)

# Mass flow within flame
m1 = HeskestadPlumeMassFlowBelow(Qc,z1,L)

# Mass flow above flame 
m2 = HeskestadPlumeMassFlowAbove(Qc,z2,z0)

# Centerline temperature above flame
T = HeskestadPlumeTemperatureAbove(Qc,z2,z0)


print 'Virtual origin: ' + str(z0) + ' m'
print 'Mean flame height: ' + str(L) + ' m'
print 'Convective ERR: ' + str(Qc) + ' kW'
print 'Mass flow (2.5 m): ' + str(m1) + ' kg/s'
print 'Mass flow (6.0 m): ' + str(m2) + ' kg/s'
print 'Temperature: ' + str(T+293) + ' K'
'''


#%%

# McCaffrey plume equations
# Description:
# Based on experimental data and dimensional analysis, equations for upward 
# velocity and plume temperature were developed. Plume devided into three
# regions: Flame, intermittent and plume. Equations use total ERR

# Table for the coefficients:
# Region        z/Q**(2/5)  n       k
# Continuous    <0.08       1/2     6.8
# Intermittend  0.08 - 0.2  0       1.9
# Plume         >0.2        -1/3    1.1

# Determine the plume region
# Table 4.1 [1].
def McCaffreyPlumeRegion(ERR,z):
    region = z/ERR**(2.0/5.0)
    return region

# Temperature difference of the plume at height z
# Formula 4.29 [1].
def McCaffreyPlumeTemperture(ERR,k,n,z,g=9.81,T_infty=293):
    a = (k/(0.9 * np.sqrt(2*g)))**2.0
    T0 = a * (z/ERR**(2.0/5.0))**(2 * n - 1.0) * T_infty
    return T0

# Upward gas velocity of the plume at height z
# Formula 4.30 [1].
def McCaffreyPlumeVelocity(ERR,k,n,z):
    u0 = k * (z/ERR**(2.0/5.0))**n * ERR**(1.0/5.0)
    return u0


'''
# Test based on example 4.4 [1].
# Plume region at 2.5 m should be ~0.13.
# Coefficients at 2.5 m should be n=0.0, k=1.9.
# Plume region at 6.0 m should be ~0.30.
# Coefficients at 6.0 m should be n=-1.0/3.0, k=1.1.
# Be aware, the temperature in McCaffreys experiments was measured in Celsius.
# Temperature at 2.5 m should be  ~540 °C (Temperature difference: ~520 K).
# Temperature at 6.0 m should be ~180 °C (Temperature difference: ~160 K).
# Temperature should be ~433 K (Temperature difference: ~140 K).

ERR = 1690

z1 = 2.5
z2 = 6.0
k1 = 1.9
n1 = 0.0
k2 = 1.1
n2 = -1.0/3.0


# Determine region at 2.5 m
r1 = McCaffreyPlumeRegion(ERR,z1)

# Determine region at 6.0 m
r2 = McCaffreyPlumeRegion(ERR,z2)

# Centerline temperature at 2.5 m
T1 = McCaffreyPlumeTemperture(ERR,k1,n1,z1)

# Centerline temperature at 6.0 m
T2 = McCaffreyPlumeTemperture(ERR,k2,n2,z2)


print 'Region at 2.5 m: ' + str(r1)
print 'Region at 6.0 m: ' + str(r2)
print 'Temperature at 2.5 m: ' + str(T1+20) + ' K'
print 'Temperature at 6.0 m: ' + str(T2+20) + ' K'
'''


#%%



#%%




































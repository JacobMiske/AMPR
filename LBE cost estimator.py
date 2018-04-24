# -*- coding: utf-8 -*-
"""
Created on Fri Apr  6 13:40:27 2018

@author: Jacob Miske
"""
#AMPR Helper Code

import numpy as np


# Analyse cost of lead and bismuth with respect 
#to LBE (Lead Bismuth Eutectic) melting point

#http://www.infomine.com/investment/metal-prices/lead/
CostofPb = 1.080; #USD/lb

#Solar panels use some Bismuth coatins
CostofBi = 0.39; #USD/gram
#https://bitinfocharts.com/comparison/bismuth-price.html#log
gramsPerLb = 453.592;
CostofBi = 14.12 #USD/lb
#conversion factor
Lbtokg = 0.453592

print("CostofPb") #USD/lb
print(CostofPb) #USD/lb
print("CostofBi") #USD/lb
print(CostofBi) #USD/lb
print("in USD/lb") #USD/lb

#Great paper on LBE properties
#https://www.tandfonline.com/doi/pdf/10.1080/18811248.2006.9711131
#http://www.iaea.org/inis/collection/NCLCollectionStore/_Public/43/095/43095088.pdf
#At 43.7% lead, lowest melting point of 400K

#At room temperature
rhoPb = 11.34 #g/cm3
rhoBi = 9.78 #g/cm3

#Parameters of Design in AMPR
H = 10 #Height of core tank [m]

#mass flow rate [kg/s]
m_prime = 10000
#Fundamental Constants
g = 9.81


def TmeltLBE(Pb):
    #Takes in percent Pb and percent Bi
    #Gives out melting point of LBE resultant using linear model
    if Pb < 43.7:
        #one side of the eutectic point
        return -3.30709382*Pb+544.52
    else:
        #other side of the eutectic point
        return 3.563055062*(Pb-43.7)+400

def CostLBE(Pb):
    #Takes in percent of Pb, defines percent of Bi
    #returns cost of LBE in USD/lb
    Bi = 100 - Pb
    return (Bi*CostofBi + Pb*CostofPb)/100

def LBEpropertyAsT (T):
    #Takes in classic LBE (45% Pb, 55% Bi)
    #and gives out list of basic properties
    
    #Saturated vapor pressure
    p_s = 11.1*10**9*np.exp(-22552/T)
    #Surface tension
    sigma = (437.1 - 0.066*T)*0.001
    #Density
    rho = 11096 - 1.3236*T
    #Sound Velocity in fluid
    u_sound = 1773+0.1049*T-2.873*0.0001*T**2
    #Bulk modulus
    B_s = (35.18 - 1.541*0.001*T-9.191*10**-6 *T**2)*10**9
    #Isobaric specific heat
    c_p = 159 - 2.72*10**-2 *T + 7.12*10**-6 *T**2
    #Dynamic viscosity 
    eta = 4.94*10**-4*np.exp(754.1/T)
    #Electric resistivity
    r = (86.334 + 0.0511*T)*10**-8
    #Thermal Conductivity
    lambd = 3.61 + 1.517*10**-2 *T -1.741*10**-6 *T**2
    return {"SatVapor": p_s, 
            "SurfaceTension": sigma, 
            "Density": rho,
            "SoundVelocity": u_sound,
            "BulkModulus": B_s,
            "IsobaricSpecificHeat": c_p,
            "DynamicViscosity": eta,
            "ElectricResistivity": r,
            "ThermalConductance": lambd}

def PressureDropAcrossRod(rho, eta, pitch, L):
    #Takes in density, dynamic viscosity, pitch, and rod length
    # m'/rho = V'
    #AMPR group pitch (4/24/18) currently 0.0135m
    #Modeled after Example 9.3 in T&K, nuclear systems text
    #Reynolds number
    Re = (m_prime)/(eta*pitch) #Actually more like (m'*D)/(eta*A)
    #Likely laminar
    if Re < 2300:
        frictionF = 64/Re
    else:
        #Incorrect, fix later
        frictionF = 0.0064
    #Frictional dPf
    dPf = frictionF*(L/pitch)*(m_prime**2)/(2*rho)
    #Gravitational dPg
    #Assuming 10m tall core
    dPg = rho*g*H
    return (dPf + dPg)









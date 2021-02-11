# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 16:53:58 2021
echelle_grating simulation

@author: neng.liu
"""

import numpy as np
import math, cmath, matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from sympy.solvers import solve
from sympy import Symbol
from scipy import integrate

nco = 3.50685;
wavelength_start = 1260e-9;
wavelength_stop = 1360e-9;
wavelength_center = 1310e-9;
ncl = 1.0;
neff=pd.read_csv("effective_index_v2.txt",
           delim_whitespace=True,
           skipinitialspace=True);
ai = 62/180*np.pi;
aph = 60/180*np.pi;
R =1;
D = 0.2e-7;
                 
z = np.polyfit(neff["Wavelength(microns)"]*1e-6, neff["Effective_Index"], 1);

dneffdivdwl = z[0];
ngdivneff = 1.9838;
resolution = 0.125

FSR = 20e-9;
m = Symbol('m')
sol = solve(wavelength_center/m*pow(1-(m+1)/m*(1-ngdivneff),-1)-FSR,m);

lambda_0 =np.linspace(wavelength_start, wavelength_stop, round((wavelength_stop-wavelength_start)*1e9/resolution))
neff_0 = z[0]*lambda_0+z[1];
d=32*wavelength_center/(neff_0[400]*(np.sin(ai)+np.sin(aph)))


#diffraction equation output
T =[]
for i in range(len(lambda_0)):
    k = 2*np.pi*neff_0[i]/lambda_0[i];
    def f(y):
        return np.exp(1j*k*y*np.sin(ai))*(np.exp(-1j*k*198e-6)/np.sqrt(198e-6))*(np.cos(ai)+np.cos(aph))/2
    v,err =integrate.quad(f,-0.5*d,0.5*d)
    eout= R*np.sqrt(neff_0[i]/lambda_0[i])*1*v
    t= np.absolute(eout)**2
    T.append(t)
    
#plot figure    
plt.figure(0)
fig1 = plt.plot(lambda_0*1e9,10*np.log10(T), label='output1', color='blue');
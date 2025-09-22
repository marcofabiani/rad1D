"""
1D model for radiative heat flux in a liquid rocket engine [1,2]
created by Marco Fabiani, Sett 2025

Requires: rocketcea (pip install rocketcea), numpy, and matplotlib (optional). All other routines are included in tools.py

Input: profile.dat: nozzle profile in a two-column text file (x[m],y[m])
       wall temperature, wall emissivity, OF ratio, chamber pressure, oxidizer and fuel names 
Output: q.dat: output data in a text file (x[m], y[m], T[K], P[bar], xH2O, xCO2, eps, q[W/m2])

References: 
1) Marco Fabiani, Mario Tindaro Migliorino, Daniele Bianchi, and Francesco Nasuti, 
"Spectral and Global Radiative Heat Transfer Models for Liquid Propellant Rocket Engines, " 
JPP, Vol. 41, No. 5 (2025), pp. 650-664 doi: doi/abs/10.2514/1.B39892

2) Marco Fabiani, Mario Tindaro Migliorino, Daniele Bianchi, and Francesco Nasuti, 
"Reduced-order Models for Radiative Heat Loads Estimation in Liquid Rocket Engines,"
IAC 2025, Sydney
"""
from tools import *
import numpy as np 
from rocketcea.cea_obj import CEA_Obj


# input parameters
Tw = 700. # wall temperature in K
epsw = 1.0 # wall emissivity
OF = 5  # oxidizer to fuel ratio
pc = 200  # chamber pressure in bar
oxidizerName = 'O2' # oxidizer name
fuelName = 'RP-1' # fuel name

# create CEA object
cea = CEA_Obj(oxName=oxidizerName, fuelName=fuelName)

# read nozzle profile
profile = np.loadtxt('profile.dat',skiprows=0) 
xnoz = profile[:,0]
ynoz = profile[:,1]


# find nozzle throat
yt = min(ynoz)
it = np.argmin(ynoz)

# area ratios
areas = (ynoz/yt)**2

i = 0
T,p,xH2O,xCO2,epsg,q =  (np.zeros(len(ynoz)) for i in range(6))

# for each station along the nozzle
for  a_at in areas:
    # run CEA with distinction for subsonic and supersonic conditions
    if i<=it:
        outstr = cea.get_full_cea_output(Pc=pc, 
                    MR=OF, eps=None, 
                    subar=a_at, 
                    PcOvPe=None, 
                    frozen=0, 
                    frozenAtThroat=0, 
                    short_output=0, 
                    show_transport=1, 
                    pc_units='bar', 
                    output='si', 
                    show_mass_frac=False, 
                    fac_CR=None)
    elif i>it:
        outstr = cea.get_full_cea_output(Pc=pc, 
                    MR=OF, eps=a_at, 
                    subar=None, 
                    PcOvPe=None, 
                    frozen=0, 
                    frozenAtThroat=0, 
                    short_output=0, 
                    show_transport=1, 
                    pc_units='bar', 
                    output='si', 
                    show_mass_frac=False, 
                    fac_CR=None)
        
    df = extract_chamber_throat_exit_and_molefractions(outstr, normalize_keys=True)

    T[i] = float(df.loc['T, K',"EXIT"])
    p[i] = float(df.loc['P, BAR',"EXIT"])
    mw = float(df.loc['M, (1/n)',"EXIT"])

    if 'H2O' in df.index:
        xH2O[i] =  float(df.loc['H2O', "EXIT"])*mw/mw_h2o
    if 'CO2' in df.index:
        xCO2[i] =  float(df.loc['CO2', "EXIT"])*mw/mw_co2


    epsg[i] = wsgg(T[i],p[i],2*ynoz[i],xH2O[i],xCO2[i])
    q[i] = sigma*T[i]**4*epsg[i]

    i = i+1


data = np.column_stack((xnoz,ynoz,T,p,xH2O,xCO2,epsg,q))
np. savetxt(fname='q.dat',X=data,fmt='%1.4e',header='"x[m]", "y[m]", "T[K]", "P[bar]", "xH2O", "xCO2", "eps", "q[W/m2]"')


import matplotlib.pyplot as plt 

fig, ax1 = plt.subplots(figsize=(8, 6))

ax1.plot(xnoz,q/1e6,label='radiative heat flux')
ax1.set_xlabel('x, m')  # Label per l'asse X
ax1.set_ylabel('q, MW/m2')  # Label per l'asse Y

ax2 = ax1.twinx()
ax2.set_ylabel('r, m')  # Label per l'asse Y

ax2.plot(xnoz,ynoz,'--',color = 'gray')
ax1.grid(True,which="both",color='gray',alpha=1,linestyle=':')  # Aggiunge la griglia
plt.savefig('qrad.pdf',format='pdf')


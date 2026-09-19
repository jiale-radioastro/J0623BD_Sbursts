import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar

# Relevant constants in Gaussian units
C_k = 1.380649e-16  # Boltzmann constant  in erg/K
C_h = 6.62607015e-27     # Planck constant  in erg/Hz
C_G = 6.67428e-8     # Gravitation constant  in cm^3 g^{-1} s^{-2}
C_Me = 9.1093897e-28     # Mass of electrons in g
C_Mp = 1.6726231e-24     # Mass of protons in g
C_e = 4.8032068e-10    # Elementary charge in esu
C_eV = 1.602176634e-12  # Electronvolt  in erg
C_Ry = 2.17987e-11       # Hydrogen ionization energy  in erg
C_year = 3.15567e7
C_c = 2.99792458e10    # Speed of light in cm/s
C_Jy = 1e-23 # Jansky in erg s^{-1} cm^{-2} Hz^{-1}
C_eV = 1.602176634e-12 # Electronvolt in erg
C_sigma = 5.670374419e-5  # Stefan-Boltzmann constant in erg cm^{-2} s^{-1} K^{-4}
C_pc = 3.0857e18  # Distance of one parsec in cm
C_AU = 1.495978707e13 # Astronomical unit in cm

C_Msun = 1.99e33   # Solar mass in g
C_Rsun = 6.96e10  # Solar radius in cm
C_Lsun = 3.828e33  # Solar luminosity in erg/s 
C_Mdot_sun = 2e-14 * C_Msun / C_year  # Solar mass loss rate in g/s
C_Rjup = 7.1492e9 # Jupiter radius in cm
C_Mjup = 1.898e30 # Jupiter mass in g
C_Mdot_jup = 1e6  # Jupiter mass loss rate in g/s


T_eff = 743
log_g = 4.70
M_star = 22.22 * C_Mjup
R_star = 1.01 * C_Rjup
Trot_star = 1.919475*3600
Omega_star = 2*np.pi/Trot_star
tau = 1e9*C_year
B_equ = 500
mu = 0.5
p_index1 = 3
n_base_list = np.logspace(6,9,1000)
L_k=(C_G*M_star/Omega_star**2)**(1/3)/R_star
J_star = 0.205*M_star*R_star**2*Omega_star

def f1(x):
        dPhi = -0.5*(Omega_star*x*R_star)**2-C_G*M_star/R_star*(1/x-1)
        P_rot = 0.5*(2*mu)*C_Mp*n_base*np.exp(-mu*C_Mp/kB_T*dPhi)\
                *(Omega_star*x*R_star)**2
        P_mag = 1/(8*np.pi)*B_equ**2*(1/x)**(2*p_index1)
        return P_rot-P_mag

fig = plt.figure(figsize=(8,4),dpi=200)
plt.subplot(122)

kB_T = 30*C_eV  
L_max1 = []
for n_base in n_base_list:
    sol = root_scalar(f1, bracket=[2, 15])
    L_max1.append(sol.root)
plt.plot(n_base_list,L_max1,color='C00')
plt.text(2e6,4.4,r'$k_BT = 30$ eV',color='C00',rotation=-7)

kB_T = 100*C_eV  
L_max2 = []
for n_base in n_base_list:
    sol = root_scalar(f1, bracket=[2, 15])
    L_max2.append(sol.root)
plt.plot(n_base_list,L_max2,color='C01')
plt.text(2e6,5.4,r'$k_BT = 100$ eV',color='C01',rotation=-15)

kB_T = 1000*C_eV  
L_max3 = []
for n_base in n_base_list:
    sol = root_scalar(f1, bracket=[2, 15])
    L_max3.append(sol.root)
plt.plot(n_base_list,L_max3,color='C02')
plt.text(2e6,7.3,r'$k_BT = 1$ keV',color='C02',rotation=-25)

plt.plot([1e6,1e10],[L_k,L_k],':',color='k')
plt.text(2e6,2.2,'Keplerian radius',color='k',fontsize=10)

plt.xlim([1e6,1e9])
plt.xlabel(r'Base density $n_{\rm base}$ [cm$^{-3}$]')
plt.ylabel(r'Alfvén radius $R_A$ [$R_*$]')
plt.xscale('log')
plt.ylim([1,12])
plt.text(1e6*1.2,1.2,'b',color='k',fontsize=12,fontweight='bold')


plt.subplot(121)
Mdot = np.logspace(0,6,1000)*C_Mdot_jup
cond_Jup1 = 4.5e10
L_hill1 = (2*np.pi*cond_Jup1*B_equ**2*R_star**2/(C_c**2*Mdot)) **0.25
cond_Jup2 = 4.5e12
L_hill2 = (2*np.pi*cond_Jup2*B_equ**2*R_star**2/(C_c**2*Mdot)) **0.25

plt.fill_between(Mdot,L_hill1,L_hill2,facecolor='C05',alpha=0.3)
plt.text(3e9,13,r'Hill-Pontius radius',rotation=-25,color='C05',fontsize=10)
plt.fill_between([1e6,1e12],[L_max1[-1],L_max1[-1]],[L_max1[0],L_max1[0]],facecolor='C00',alpha=0.2)
plt.plot([3e6,3e6],[L_max1[0],L_max1[-1]],color='C00')
plt.plot([3e6*0.8,3e6*1.25],[L_max1[0],L_max1[0]],color='C00')
plt.plot([3e6*0.8,3e6*1.25],[L_max1[-1],L_max1[-1]],color='C00')
plt.text(3e6/2.3,12,'30 eV',color='C00')

plt.fill_between([1e6,1e12],[L_max2[-1],L_max2[-1]],[L_max2[0],L_max2[0]],facecolor='C01',alpha=0.2)
plt.plot([25e6,25e6],[L_max2[0],L_max2[-1]],color='C01')
plt.plot([25e6*0.8,25e6*1.25],[L_max2[0],L_max2[0]],color='C01')
plt.plot([25e6*0.8,25e6*1.25],[L_max2[-1],L_max2[-1]],color='C01')
plt.text(10e6,12,'100 eV',color='C01')

plt.fill_between([1e6,1e12],[L_max3[-1],L_max3[-1]],[L_max3[0],L_max3[0]],facecolor='C02',alpha=0.2)
plt.plot([200e6,200e6],[L_max3[0],L_max3[-1]],color='C02')
plt.plot([200e6*0.8,200e6*1.25],[L_max3[0],L_max3[0]],color='C02')
plt.plot([200e6*0.8,200e6*1.25],[L_max3[-1],L_max3[-1]],color='C02')
plt.text(120e6,12,'1 keV',color='C02')

plt.plot(Mdot,np.sqrt(J_star/Mdot/R_star**2/Omega_star/tau*1),color='k',ls=':',alpha=0.4)
plt.text(6e6*150,370,r'$1$ Gyr',rotation=-45,alpha=0.6)

plt.plot(Mdot,np.sqrt(J_star/Mdot/R_star**2/Omega_star/tau*1e-1),color='k',ls=':',alpha=0.4)
plt.text(6e6*15,300,r'$10$ Gyr',rotation=-45,alpha=0.6)

plt.plot(Mdot,np.sqrt(J_star/Mdot/R_star**2/Omega_star/tau*1e-2),color='k',ls=':',alpha=0.4)
plt.text(2e7,200,r'$10^2$ Gyr',rotation=-45,alpha=0.6)

plt.plot(Mdot,np.sqrt(J_star/Mdot/R_star**2/Omega_star/tau*1e-3),color='k',ls=':',alpha=0.4)
plt.text(4e6,130,r'$10^3$ Gyr',rotation=-45,alpha=0.6)

plt.text(2e9,4.7,'Alfvén radius',color='k',fontsize=10)
plt.plot([1e6,1e12],[L_k,L_k],':',color='k')
plt.text(2e9,2.3,'Keplerian radius',color='k',fontsize=10)

plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'Mass loss rate $\dot{M}$ [g s$^{-1}$]')
plt.ylabel(r'Radial distance [$R_*$]')
plt.ylim([1,1e3])
plt.xlim([1e6,1e12])
plt.text(1.8e6,1.2,'a',color='k',fontsize=12,fontweight='bold')
plt.tight_layout()

fig.savefig('figures/figure4.pdf',dpi=200,bbox_inches='tight')
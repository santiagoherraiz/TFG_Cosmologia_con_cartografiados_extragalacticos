####CÓDIGO PARA COMPARAR LOS RESULTADOS CON CAMB 
# ver https://camb.info/

import sys, platform, os
import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import camb
from camb import model, initialpower

#################################################################3
####PARA OBTENER SOLO P(k) DE CAMB
##Parámetros originales de CAMB
pars = camb.CAMBparams()
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
pars.InitPower.set_params(ns=0.965)
#Note non-linear corrections couples to smaller scales than you want
pars.set_matter_power(redshifts=[0], kmax=5.0)

#Linear spectra
pars.NonLinear = model.NonLinear_none
results = camb.get_results(pars)
kh, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=5, npoints = 200)
s8 = np.array(results.get_sigma8())

#Non-Linear spectra (Halofit)
pars.NonLinear = model.NonLinear_both
results.calc_power_spectra(pars)
kh_nonlin, z_nonlin, pk_nonlin = results.get_matter_power_spectrum(minkh=1e-4, maxkh=5, npoints = 200)

for i, (redshift, line) in enumerate(zip(z,['-','--'])):
    plt.loglog(kh, pk[i,:], color='k', ls = line)
    plt.loglog(kh_nonlin, pk_nonlin[i,:], color='r', ls = line)
plt.xlabel('$k$ [h/Mpc]');
plt.ylabel('$P(k)$ [Mpc$^3$/h$^3$]');
# plt.autoscale(enable=True, axis='x', tight=True) #ajustar plot a los bordes
plt.legend(['Lineal','No lineal'], loc='lower left');
plt.title('Espectro de potencias $P(k)$ para redshift $z=0$');

###############################################
##PARA COMPARAR CON GIGGLEZ
###Parámetros de TAO-GiggleZ
H0_tao=70.5
h_tao=0.705
OM_tao=0.273 #omch2 es materia OSCURA. restarle Ob
Ob_tao=0.0456
n_tao=0.96 #????

pars = camb.CAMBparams()
pars.set_cosmology(H0=H0_tao, ombh2=Ob_tao*(h_tao**2), omch2=(OM_tao - Ob_tao)*(h_tao**2))
pars.InitPower.set_params(ns=n_tao)
#Note non-linear corrections couples to smaller scales than you want
pars.set_matter_power(redshifts=[0], kmax=5.0)

#Linear spectra
pars.NonLinear = model.NonLinear_none
results = camb.get_results(pars)
kh, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=5, npoints = 200)
s8 = np.array(results.get_sigma8())

#Non-Linear spectra (Halofit)
pars.NonLinear = model.NonLinear_both
results.calc_power_spectra(pars)
kh_nonlin, z_nonlin, pk_nonlin = results.get_matter_power_spectrum(minkh=1e-4, maxkh=5, npoints = 200)

for i, (redshift, line) in enumerate(zip(z,['-','--'])):
    plt.loglog(kh, pk[i,:], color='k', ls = line)
    plt.loglog(kh_nonlin, pk_nonlin[i,:], color='r', ls = line)
plt.xlabel('$k$ [h/Mpc]');
plt.ylabel('$P(k)$ [Mpc$^3$/h$^3$]');
# plt.autoscale(enable=True, axis='x', tight=True) #ajustar plot a los bordes
plt.legend(['Lineal','No lineal'], loc='lower left');
plt.title('Espectro de potencias $P(k)$ para redshift $z=0$');
# pl.errorbar(k,P,yerr=sigma,fmt='ro',label='estimated P(k)')
# pl.axvline(x=kNy,color='y',label='Nyquist frequency')
# plt.scatter(kfin[1:nk], Pk[1:nk]);

###############################################################3
##PARA COMPARAR CON MICECAT
###Parámetros de MICECAT
H0_mice=70
h_mice=0.70
OM_mice=0.25 #omch2 es materia OSCURA. restarle Ob
Ob_mice=0.044
n_mice=0.95 #????

pars = camb.CAMBparams()
pars.set_cosmology(H0=H0_mice, ombh2=Ob_mice*(h_mice**2), omch2=(OM_mice - Ob_mice)*(h_mice**2))
pars.InitPower.set_params(ns=n_mice)
#Note non-linear corrections couples to smaller scales than you want
pars.set_matter_power(redshifts=[0], kmax=5.0)

#Linear spectra
pars.NonLinear = model.NonLinear_none
results = camb.get_results(pars)
kh, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=5, npoints = 200)
s8 = np.array(results.get_sigma8())

#Non-Linear spectra (Halofit)
pars.NonLinear = model.NonLinear_both
results.calc_power_spectra(pars)
kh_nonlin, z_nonlin, pk_nonlin = results.get_matter_power_spectrum(minkh=1e-4, maxkh=5, npoints = 200)

for i, (redshift, line) in enumerate(zip(z,['-','--'])):
    plt.loglog(kh, pk[i,:], color='k', ls = line)
    plt.loglog(kh_nonlin, pk_nonlin[i,:], color='r', ls = line)
plt.xlabel('$k$ [h/Mpc]');
plt.ylabel('$P(k)$ [Mpc$^3$/h$^3$]');
# plt.autoscale(enable=True, axis='x', tight=True) #ajustar plot a los bordes
plt.legend(['Lineal','No lineal'], loc='lower left');
plt.title('Espectro de potencias $P(k)$ para redshift $z=0$');
# pl.errorbar(k,P,yerr=sigma,fmt='ro',label='estimated P(k)')
# pl.axvline(x=kNy,color='y',label='Nyquist frequency')
# plt.scatter(kfin[1:nk], Pk[1:nk]);













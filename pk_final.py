###Código para obtener P(k) para catálogos de galaxias (casos específicos de GiggleZ y MICECAT, pero es generalizable)
import matplotlib.pyplot as plt
import numpy as np
#####################################################
##Para datos de GiggleZ
from astropy.visualization import astropy_mpl_style 
plt.style.use(astropy_mpl_style)
from astropy.io import fits

datos_tao = fits.open('tao.4531.0.fits') #Nombre del fichero, contiene x,y,z

data=datos_tao[1].data

H0_tao=70.5
h_tao=0.705
OM_tao=0.273 #omch2 es materia OSCURA. restarle Ob
Ob_tao=0.0456
n_tao=0.96

x=data["X"]
y=data["Y"]
z=data["Z"]

######################################################
 ###Para datos de MICECAT
import pandas as pd
file_cat = '17203.csv.bz2' #Nombre del fichero, contiene x,y,z
df=pd.read_csv(file_cat, comment='#').replace(r'\N', np.nan)

x=df["x_c"].values
y=df["y_c"].values
z=df["z_c"].values
redshift=df["z"].values

#######################################################
##para hacer cortes en distancias
index_cut = np.where((x<500.)&(y<500.)&(z<500.))[0]   
x=x[index_cut]
y=y[index_cut]
z=z[index_cut]
redshift=redshift[index_cut]

##para hacer cortes en redshift
# redshift_cut = np.where((redshift>0.1)&(redshift<0.3))[0] 
# x=x[redshift_cut]
# y=y[redshift_cut]
# z=z[redshift_cut]
# redshift=redshift[redshift_cut]

########################################################
#### CONTEO DE GALAXIAS
n=100
L_box=500 #h^-1*Mpc (3072)
bins=np.linspace(0,L_box,n+1)

xdig=np.digitize(x,bins)
ydig=np.digitize(y,bins)
zdig=np.digitize(z,bins)

def contarpuntoscubodig(a,b,c):
    contdig=np.zeros(len(x))
    for i in range(len(x)):
        if xdig[i]==a+1 and ydig[i]==b+1 and zdig[i]==c+1:
            contdig[i]=1
        else:
            contdig[i]=0
    Nd=sum(contdig)
    return Nd

#### OBTENCIÓN DE n_s
#hacer el tensor N(a,b,c)
Ndig=np.zeros([n,n,n])
for a in range(len(bins)-1):
    for b in range(len(bins)-1):
        for c in range(len(bins)-1):
            Ndig[a,b,c]=contarpuntoscubodig(a,b,c)
            
N1D=np.zeros(n**3)
for a in range(len(bins)-1):
    for b in range(len(bins)-1):
        for c in range(len(bins)-1):
            i=a*(n**2)+b*n+c
            N1D[i]=Ndig[a,b,c]

N_mean=len(x)/(n**3)
nbar=np.ones((n,n,n))*N_mean

#hacer el tensor Nrand(a,b,c)
Nrand=np.zeros([n,n,n])
for a in range(len(bins)-1):
    for b in range(len(bins)-1):
        for c in range(len(bins)-1):
            Nrand[a,b,c]=np.random.rand()*N_mean*1000

Ntotrand=np.sum(Nrand)
Nrand_mean=Ntotrand/(n**3)
Nrand_mean_3D=Nrand_mean*np.ones([n,n,n])

#### SHOT NOISE
# P0=10**4
P0=5000.0*((n/L_box)**3)  
b=1 #depende del tipo de galaxia 
w=((b**2)*P0) / (1.0+nbar*P0*(b**2)) 
alpha=(np.sum(w*Ndig/b))/(np.sum(w*Nrand/b))  
N=np.sqrt(np.sum((nbar**2)*(w**2)))
Pshot=((1+alpha)/(N**2)) * np.sum(nbar*((w**2)/(b**2))) 

####CAMPO DE SOBREDENSIDADES F(k) (las que usamos para obtener P(k))
delta = w*(Ndig - N_mean*np.ones([n,n,n]) - alpha*(Nrand - Nrand_mean))/N  

# delta1D=np.zeros(n**3)
# for a in range(len(bins)-1):
#     for b in range(len(bins)-1):
#         for c in range(len(bins)-1):
#             i=a*(n**2)+b*n+c
#             delta1D[i]=delta[a,b,c]

dfou=np.fft.fftn(delta)
d2mod=(dfou*dfou.conj()).real

#################################
#### ANÁLISIS DE FOURIER, FRECUENCIAS DE MUESTREO
nk=70
kfft=np.fft.fftfreq(n) 
kminfft=np.amin(np.abs(kfft))
kmaxfft=np.amax(np.abs(kfft))
kdistmax=np.sqrt(3)*kmaxfft
kdistmin=kminfft

k_bins=np.linspace(kdistmin,kdistmax,nk+1) 

#### OBTENCIÓN DE P(k)
P=np.zeros(nk)
counts=np.zeros(nk)

for i1 in range(len(kfft)):
    for i2 in range(len(kfft)):
        for i3 in range(len(kfft)):
            kx=kfft[i1]
            ky=kfft[i2]
            kz=kfft[i3]
            kdist=np.sqrt(kx**2 + ky**2 + kz**2) #"distancia" para las k
            for m in range(len(k_bins)-1):
                if (kdist>=k_bins[m] and kdist<=k_bins[m+1]):
                    P[m]=P[m]+d2mod[i1,i2,i3] - Pshot
                    counts[m]=counts[m]+1
                    break

Pfin=P/counts

k=np.zeros(len(Pfin))
for i in range(len(k_bins)-1): #para que cada valor de P quede en el centro del bin
    k[i]=(k_bins[i]+k_bins[i+1])/2

#### CAMBIO A UNIDADES FÍSICAS

Pk=Pfin*(L_box/n)**3 
P_SN=Pshot*(L_box/n)**3 
kfin=k*(2*np.pi*n/L_box) 


###############################
#### PLOT CON LOS RESULTADOS
plt.xscale('log')
plt.yscale('log')
plt.xlabel('k [h$\cdot$Mpc$^{-1}$]');
plt.ylabel('P(k) [h$^{-3}\cdot$Mpc$^3$]')
# plt.title('Matter power at z=0');
plt.axhline(y=P_SN, color='r', linestyle='-')
plt.scatter(kfin[1:nk], Pk[1:nk]);
# plt.savefig('fig_save.png')
# np.savetxt('file_save.txt',list(zip(kfin, Pk)))



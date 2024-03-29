import netCDF4 as nc
import datetime
import time
import warnings
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
warnings.filterwarnings('ignore')
import numpy as np
import geopy.distance
from scipy import interpolate
import scipy.io
from scipy.ndimage import zoom

fn='/home/joao/Saidas/2008010100/archv.2008_002_00_3zu.nc'

ds=nc.Dataset(fn)
print(ds.variables.keys())
for var in ds.variables.values():
    print(var)
#vendo as dimensoes das variaves
for dim in ds.dimensions.values():
    print(dim)

lat=ds['Latitude'][:]
lon=ds['Longitude'][:]
depth=ds['Depth'][:]
u=ds['u']
u=u[0,:,:,:]
v=ds['v']
v=v[0,:,:,:]

p1 = [-8.798547, -35.208374] #[lat lon]
p2 = [-9.497866, -32.5] #[lat lon]

#calculando o angulo, com base no script de felipe
dy = geopy.distance.geodesic([p1[0],p1[1]], [p2[0],p1[1]]).km
dx = geopy.distance.geodesic([p1[0],p1[1]], [p1[0],p2[1]]).km
ang = 180*np.arctan(dy/dx)/np.pi
theta = np.deg2rad(ang)
print(ang)
print(theta)

#com o angulo theta, é possível calcular a velocidade rotacionada (apenas para esta seção, caso mude o angulo eu preciso calcular outro valor de v')
v_linha=v*np.cos(theta)-u*np.sin(theta)

print(np.where(lat==p1[0]))
#lat[658]
print(np.where(lat==p2[0]))
#lat[641]
lat_sec=lat[641:658]


print(np.where(lon==p1[1]))
#lon[451]
print(np.where(lon==p2[1]))
#lon[516]

lon_sec=lon[451:516]


#selecionando e plotando a seção transversal inclinada, mas sem interpolar
v_sec=v_linha[:,641:658,451:516:4]
print(v_sec)
v_teste=[]
for i in range(17):
    v_teste.append(v_sec[0:19,i,i])

v_teste_array2=np.vstack(v_teste).T


levels=np.linspace(-1,1,100)
plt.contourf(lon[451:519:4],-depth[0:19],v_teste_array2[:,:],60,levels=levels,cmap='jet')
plt.colorbar()
plt.suptitle("V rotacionada sem interpolação")
plt.show()

#As velocidades devem ser interpoladas, já que as dimensões lat e lon possuem tamanhos diferentes

#como interpolar?
v_sec2=v_linha[:,641:658,451:516]
v_interp=np.zeros_like(v_sec)

original_shape=v_sec2.shape
new_shape=(original_shape[0],65,65)
zoom_factors = (1, new_shape[1] / original_shape[1], new_shape[2] / original_shape[2])

# Redimensionar o array usando interpolação linear
resized_array = zoom(v_sec2, (zoom_factors), order=1)

v_teste_interp=[]
for i in range(65):
    v_teste_interp.append(resized_array[0:19,i,i])

v_teste_interp2=np.vstack(v_teste_interp).T

plt.contourf(lon[451:516],-depth[0:19],v_teste_interp2[:,:],60,levels=levels,cmap='jet')
plt.colorbar()
plt.suptitle("V rotacionada com interpolação")
plt.show()

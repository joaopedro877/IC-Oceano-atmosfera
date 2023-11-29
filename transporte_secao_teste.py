import netCDF4 as nc
import datetime
import time
import warnings
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
warnings.filterwarnings('ignore')
import numpy as np
import geopy.distance

fn='/home/joao/SAIDAS/2008010100/archv.2008_002_00_3zu.nc'
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

#print(lat)
print(depth)
#print(np.where(lat==-30.885464))
print(np.where(lon==-35.75))
#print(lon)
#diminuir seção em longitude e profundidade - até 1000m

levels=np.linspace(-1,1,100)
#plt.contourf(lon[367:438],-depth[0:19],v[0:19,554,367:438],60,levels=levels,cmap='jet')
#plt.suptitle(str(lat[554]))
#plt.colorbar()
#plt.show()

                                                   #Calculando o transporte dessa secao

#preciso calcular o tamanho em metros de cada célula
tamanho_da_secao=geopy.distance.geodesic([lat[554],lon[438]],[lat[554],lon[367]]).m
print(tamanho_da_secao)
print(tamanho_da_secao/(438-367))
print(depth[0:19])
v_secao=v[0:19,554,367:438]
#calculando o volume para sul da secao, em Sv
print(v_secao.shape)
transporte = v_secao*4519*10
transporte_Sv=transporte/(10**6)
plt.contourf(lon[367:438],-depth[0:19],transporte_Sv,60,cmap='jet')
plt.suptitle(str(lat[554]))
plt.colorbar()
plt.show()



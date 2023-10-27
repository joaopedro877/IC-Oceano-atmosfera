import netCDF4 as nc
import datetime
import time
import warnings
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
warnings.filterwarnings('ignore')
import numpy as np

fn='/home/joao/SAIDAS/2008010100/archv.2008_002_00_3zu.nc'
ds=nc.Dataset(fn)
print(ds.variables.keys())
for var in ds.variables.values():
    print(var)
#vendo as dimensoes das variaves
for dim in ds.dimensions.values():
    print(dim)

lat=ds['Latitude'][:]
lon=ds['Longitude'][:]o
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
plt.contourf(lon[367:438],-depth[0:19],v[0:19,554,367:438],60,levels=levels,cmap='jet')
plt.suptitle(str(lat[554]))
plt.colorbar()
plt.show()




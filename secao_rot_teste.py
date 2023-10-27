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
lon=ds['Longitude'][:]
depth=ds['Depth'][:]
u=ds['u']
u=u[0,:,:,:]
v=ds['v']
v=v[0,:,:,:]



#print(lat)
"testando a seção transversal rotacionada"
p1 = [-22.00, -41.10] #[lat lon]
p2 = [-22.00, -39.00] #[lat lon]
#lat[506:536]
#lon[338:404]

'''
ax = plt.axes(projection=ccrs.PlateCarree())
plt.title('Costa da Bahia')
ax.set_extent([-50,-35, -30,-11], ccrs.PlateCarree())
ax.coastlines(resolution='10m')
g = ax.gridlines(crs=ccrs.PlateCarree(), linestyle='-.', color='gray', draw_labels=True)
g.ylabels_right = False
g.xlabels_top = False
#plot da seção transversal
plt.plot([p1[1],p2[1]], [p1[0],p2[0]], transform=ccrs.Geodetic())
#ax.plot(-39.166626,-13.781812,color='red', marker='o',transform=ccrs.PlateCarree())
#ax.plot(-37.166626,-14.992585,color='blue', marker='o',transform=ccrs.PlateCarree())
plt.show()'''

print(type(np.where(lon==-39.083374)[0]))
print(type((np.where(lon==-39.083374)[0]).astype(int)))
print(int(np.where(lon==-39.083374)[0]))
#plotando a seção transversal rotacionada
levels=np.linspace(-1,1,100)
plt.contourf(lon[int(np.where(lon==-41.083374)[0]):int(np.where(lon==-39.083374)[0])],-depth[0:19],v[0:19,int(np.where(lat==-22.058943)[0]),int(np.where(lon==-41.083374)[0]):int(np.where(lon==-39.083374)[0])],60,levels=levels,cmap='jet')
plt.suptitle(str(np.where(lat==-22.058943)[0]))
plt.colorbar()
plt.show()


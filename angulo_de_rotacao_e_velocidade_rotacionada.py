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

p1 = [-8.798547, -35.208374] #[lat lon]
p2 = [-9.497866, -32.5] #[lat lon]

#calculando o angulo, com base no script de felipe
dy = geopy.distance.geodesic([p1[0],p1[1]], [p2[0],p1[1]]).km
dx = geopy.distance.geodesic([p1[0],p1[1]], [p1[0],p2[1]]).km
ang = 180*np.arctan(dy/dx)/np.pi
theta = np.deg2rad(ang)
print(ang)
print(theta)

#com o angulo theta, é possível calcular a velocidade rotacionada
v_linha=v*np.cos(ang)-u*np.sin(ang)


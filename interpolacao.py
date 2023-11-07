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



#interpolacao do chatgpt
import numpy as np
from scipy.ndimage import zoom

# Supondo que 'array_3d' é o seu array original com dimensões 33x17x65
# Vamos redimensioná-lo para 33x65x65 usando interpolação linear.

original_shape = array_3d.shape
new_shape = (original_shape[0], 65, 65)

# Fator de redimensionamento para cada dimensão
zoom_factors = (1, new_shape[1] / original_shape[1], new_shape[2] / original_shape[2])

# Redimensionar o array usando interpolação linear
resized_array = zoom(array_3d, zoom_factors, order=1)

# Agora 'resized_array' contém o array redimensionado com dimensões 33x65x65 usando interpolação linear.


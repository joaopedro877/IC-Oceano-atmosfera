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

#################################################################   Média ############################################################################
data_inicial = "01/01/2008"
data_final = "31/01/2008"

dir_cfs = '/home/joao/SAIDAS'

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])
current_data = data_inicial


lista_v=[]
counter=0

#concatenando os dados de v
while (current_data<=data_final):
    fn=(dir_cfs+"/"+current_data.strftime("%Y%m%d%H")+'/archv.'+ str(data_inicial.strftime("%Y"))+'_'+str(int(current_data.strftime("%d"))+1).zfill(3)+'_00_3zu.nc')
    ds=nc.Dataset(fn)
    lat=ds['Latitude'][:]
    lon=ds['Longitude'][:]
    #temperatura
    #criando variaveis dinamicas
    globals()["v"+str(counter)]=ds['v'][0,:,:,:]
    lista_v.append((globals()["v"+str(counter)]))
    counter =counter + 1
    current_data = current_data + datetime.timedelta(days=1)
    ds.close()

#media da velocidade v
v=np.ma.stack(lista_v,3)
v_media=np.mean(v,axis=3)

#fazendo o mesmo para u
data_inicial = "01/01/2008"
data_final = "31/01/2008"

dir_cfs = '/home/joao/SAIDAS'

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])
current_data = data_inicial

lista_u=[]
counter=0
while (current_data<=data_final):
    fn=(dir_cfs+"/"+current_data.strftime("%Y%m%d%H")+'/archv.'+ str(data_inicial.strftime("%Y"))+'_'+str(int(current_data.strftime("%d"))+1).zfill(3)+'_00_3zu.nc')
    ds=nc.Dataset(fn)
    lat=ds['Latitude'][:]
    lon=ds['Longitude'][:]
    depth=ds['Depth'][:]
    #temperatura
    #criando variaveis dinamicas
    globals()["u"+str(counter)]=ds['u'][0,:,:,:]
    lista_u.append((globals()["u"+str(counter)]))
    counter =counter + 1
    current_data = current_data + datetime.timedelta(days=1)
    ds.close()
#media da velocidade u
u=np.ma.stack(lista_u,3)
u_media=np.mean(u,axis=3)
######################################################### Coordenadas mais próximas das calculadas por Felipe ########################################


#criando 2 listas, com todos os valores de p1 e de p2
p1=[[-08.80, -35.20],[-09.60, -35.50],[-10.50, -36.50],[-12.60, -38.20],[-13.70, -39.10],[-15.00, -39.10],[-16.50, -39.00],[-16.35, -39.20],[-18.90, -39.80],[-18.30, -39.80],[-20.75, -40.70],[-22.00, -41.10],[-22.25, -42.00],[-22.80, -42.25],[-24.00, -46.00],[-24.60, -47.50],[-26.00, -48.50],[-27.25, -48.70],[-28.35, -49.00],[-29.30, -50.00],[-32.00, -51.00],[-23.00, -44.00]]

p2=[[-09.50, -32.50],[-11.10, -33.00],[-12.00, -34.30],[-13.60, -36.50],[-14.90, -37.00],[-15.00, -37.00],[-15.00, -36.00],[-18.20, -36.00],[-17.00, -35.50],[-19.60, -35.50],[-21.20, -39.00],[-22.00, -39.00],[-23.70, -39.00],[-25.50, -42.25],[-26.40, -43.00],[-26.70, -44.50],[-28.10, -44.50],[-28.55, -46.00],[-29.50, -46.00],[-30.70, -46.00],[-32.00, -48.00],[-28.00, -44.00]]
#substituindo os valores da lista pelos valores mais próximos
for i in range(len(p1)):
    print(p1[i])
    p1[i][0]=lat[min(range(len(lat)), key=lambda x: abs(lat[x]-(p1[i][0])))]
    p1[i][1]=lon[min(range(len(lon)), key=lambda y: abs(lon[y]-(p1[i][1])))]

for i in range(len(p2)):
    print(p2[i])
    p2[i][0]=lat[min(range(len(lat)), key=lambda x: abs(lat[x]-(p2[i][0])))]
    p2[i][1]=lon[min(range(len(lon)), key=lambda y: abs(lon[y]-(p2[i][1])))]

#criando uma lista com os valores de theta
theta=[]
for i in range(len(p1)):
    dy = geopy.distance.geodesic([p1[i][0],p1[i][1]], [p2[i][0],p1[i][1]]).km
    dx = geopy.distance.geodesic([p1[i][0],p1[i][1]], [p1[i][0],p2[i][1]]).km
    if dx==0:
        ang=0
    else:
        ang=180*np.arctan(dy/dx)/np.pi
    theta.append(np.deg2rad(ang))
    

print(theta)

################################################# plotando as seções transversais em um loop ##########################################################
#teste
#print(lat)
levels=np.linspace(-1,1,100)
print(len(p1))
print(len(theta))

for i in range(len(p1)):
	v_linha=v_media*np.cos(theta[i])-u_media*np.sin(theta[i])
	lat1=int(np.where(lat==p1[i][0])[0])
	lat2=int(np.where(lat==p2[i][0])[0])
	lon1=int(np.where(lon==p1[i][1])[0])
	lon2=int(np.where(lon==p2[i][1])[0])
	print("o lat1 é :" + str(lat1) + "\n e o lat2 é :" + str(lat2))
	if lat1==lat2:
		plt.contourf(lon[lon1:lon2],-depth[0:19],v_linha[0:19,lat1,lon1:lon2],60,levels=levels,cmap='jet')
		plt.suptitle(str(lat[lat1]))
		plt.colorbar()
		plt.show()
	else:
		v_sec=v_linha[:,lat2:lat1,lon1:lon2]
		original_shape=v_sec.shape
		print(original_shape)
		#a variacao de longitude sera sempre maior que a de latitude?
		new_shape=(original_shape[0],(lon2-lon1),(lon2-lon1))
		print(new_shape)
		zoom_factors=(1,new_shape[1]/original_shape[1],new_shape[2] / original_shape[2])
		print(zoom_factors)
		resized_array=zoom(v_sec,(zoom_factors),order=1)
		v_interp=[]
		for i in range(lon2-lon1):
			v_interp.append(resized_array[0:19,i,i])
		v_interp2=np.vstack(v_interp).T
		print(v_interp2.shape)
		plt.contourf(lon[lon1:lon2],-depth[0:19],v_interp2[:,:],60,levels=levels,cmap='jet')
		plt.suptitle(str(lat[lat1]))
		plt.colorbar()
		plt.show()

exit()

exit()
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


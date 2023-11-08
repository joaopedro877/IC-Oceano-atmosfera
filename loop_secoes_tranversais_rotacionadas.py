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

temp=np.empty([701,553])
lista_v=[]
counter=0

#concatenando os dados de v
while (current_data<=data_final):
    print(" DIA: "+current_data.strftime("%d-%m-%Y"))
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

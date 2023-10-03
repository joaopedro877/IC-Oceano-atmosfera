import netCDF4 as nc
import datetime
import time
import warnings
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
warnings.filterwarnings('ignore')
import numpy as np


fig_output_dir='/home/joao/scripts/plots'

data_inicial = "01/01/2008"
data_final = "03/01/2008"

dir_cfs = '/home/joao/SAIDAS'

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])
current_data = data_inicial

temp=np.empty([701,553])
lista_temps=[]
print(temp.shape)
counter=0

#concatenando os dados de sst
while (current_data<=data_final):
	print(" DIA: "+current_data.strftime("%d-%m-%Y"))
	fn=(dir_cfs+"/"+current_data.strftime("%Y%m%d%H")+'/archv.'+ str(data_inicial.strftime("%Y"))+'_00'+str(int(current_data.strftime("%d"))+1)+'_00_3zt.nc')
	ds=nc.Dataset(fn)
	lat=ds['Latitude'][:]
	lon=ds['Longitude'][:]
	#temperatura
	#criando variaveis dinamicas
	globals()["temp"+str(counter)]=ds['temperature'][0,0,:,:]
	lista_temps.append((globals()["temp"+str(counter)]))
	counter =counter + 1
	current_data = current_data + datetime.timedelta(days=1)

print(lista_temps)

sst=np.ma.stack(lista_temps,2)
print(sst.ndim)
sst_media=np.mean(sst,axis=2)
sst_desvio=np.std(sst,axis=2)

'''plotando
ax = plt.axes(projection=ccrs.PlateCarree())
plt.contourf(lon, lat,sst_media,60,transform=ccrs.PlateCarree(),cmap='RdBu_r')
cb=plt.colorbar()
ax.coastlines()
g = ax.gridlines(crs=ccrs.PlateCarree(), linestyle='-.', color='gray', draw_labels=True)
g.ylabels_right = False
g.xlabels_top = False
cb.ax.set_xlabel('Temperatura(°C)')
plt.suptitle('Temperatura média na superficie (°C) '+ ' -entre '+str(data_inicial.strftime("%d-%m-%Y")) +' e '+ current_data.strftime("%d-%m-%Y"))
plt.show()'''

'''ax = plt.axes(projection=ccrs.PlateCarree())
plt.contourf(lon, lat,sst_desvio,60,transform=ccrs.PlateCarree(),cmap='RdBu_r')
cb=plt.colorbar()
ax.coastlines()
g = ax.gridlines(crs=ccrs.PlateCarree(), linestyle='-.', color='gray', draw_labels=True)
g.ylabels_right = False
g.xlabels_top = False
cb.ax.set_xlabel('Temperatura(°C)')
plt.suptitle('Desvio da temperatura na superfície(°C) '+ ' -entre '+str(data_inicial.strftime("%d-%m-%Y")) +' e '+ current_data.strftime("%d-%m-%Y"))
plt.show()	'''
#nome das variaveis
print(dir())

#serie temporal
from datetime import timedelta, date

def get_dates_between(start_date, end_date):
    return [start_date + timedelta(days=i) 
            for i in range((end_date - start_date).days + 1)]

dates_between = get_dates_between(data_inicial, data_final)
print(dates_between)

print("A media do dia 0 é: " + str(np.mean(sst[0])))
print(dates_between)
medias=[np.mean(sst[0]),np.mean(sst[1]),np.mean(sst[2])]
print(medias)
plt.plot(dates_between,medias)
plt.suptitle("Série temporal da tsm")
plt.xlabel("data")
plt.ylabel("Temp(°C)")
plt.show()




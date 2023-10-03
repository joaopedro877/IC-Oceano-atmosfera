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
counter=0
lista_tempo=[]

while (current_data<=data_final):
	print(" DIA: "+current_data.strftime("%d-%m-%Y"))
	fn=(dir_cfs+"/"+current_data.strftime("%Y%m%d%H")+'/archv.'+ str(data_inicial.strftime("%Y"))+'_00'+str(int(current_data.strftime("%d"))+1)+'_00_3zt.nc')
	ds=nc.Dataset(fn)
	lat=ds['Latitude'][:]
	lon=ds['Longitude'][:]
	#temperatura
	temp=ds['temperature'][0,0,:,:]+temp
	counter =counter + 1
	lista_tempo.append(float(ds['temperature'][0,0,100,100]))
	current_data = current_data + datetime.timedelta(days=1)

media = temp/counter

print(sum(lista_tempo)/counter)
print(media[100,100])

ax = plt.axes(projection=ccrs.PlateCarree())
plt.contourf(lon, lat,media,60,transform=ccrs.PlateCarree(),cmap='RdBu_r')
cb=plt.colorbar()
ax.coastlines()
g = ax.gridlines(crs=ccrs.PlateCarree(), linestyle='-.', color='gray', draw_labels=True)
g.ylabels_right = False
g.xlabels_top = False
cb.ax.set_xlabel('Temperatura(°C)')
plt.suptitle('Temperatura média na superficie (°C) '+ ' -entre '+str(data_inicial.strftime("%d-%m-%Y")) +' e '+ current_data.strftime("%d-%m-%Y"))
plt.show()		

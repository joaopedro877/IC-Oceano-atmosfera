    import netCDF4 as nc
import datetime
import time
import warnings
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
warnings.filterwarnings('ignore')
import numpy as np


fig_output_dir='/home/jpedro/plots'

data_inicial = "30/01/2015"
data_final = "02/02/2015"

dir_cfs = '/home/joao/IC-Oceanografia_Fisica/Saidas'

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])
current_data = data_inicial
counter=0
#coordenadas em lat e lon do ponto inicial dos quadrados que devem ser plotados
latitudes=[-18.437976,-18.437976,-21.169922,-21.169922,-21.169922,-23.709533,-23.709533,-23.709533,-25.925299,-25.925299,-25.925299,-25.925299,-28.240898,-28.240898,-28.240898,-28.240898,-30.813785,-30.813785,-30.813785,-30.813785,-33.258028,-33.258028,-33.258028,-33.258028]
longitudes=[-37.750837,-34.358219,-40.262134,-37.750837,-34.358219,-41.841604,-38.980886,-36.115508,-45.806716,-43.148066,-40.304164,-37.588906,-47.350449,-44.734680,-41.818571,-39.374498,-48.636892,-45.292139,-41.990267,-38.902802,-50.309269,-47.521975,-44.820443,-41.304164]
#concatenando os dados de velocidade
while (current_data<=data_final):
	print(" DIA: "+current_data.strftime("%d-%m-%Y"))
	for hora in range(3,23,6):
		fn=(dir_cfs+"/"+current_data.strftime("%Y%m%d")+'00'+'/archv.'+ str(data_inicial.strftime("%Y"))+'_'+str(int(current_data.strftime("%j"))).zfill(3)+'_'+str(hora).zfill(2)+'_'+'3zu.nc')
		ds=nc.Dataset(fn)
		lat=ds['Latitude'][:]
		lon=ds['Longitude'][:]
		v=ds['v'][0,0,:,:]
		u=ds['u'][0,0,:,:]
		V=np.sqrt((v**2)+(u**2))
		ax = plt.axes(projection=ccrs.PlateCarree())
		plt.contourf(lon[120:480], lat[208:500],V[208:500,120:480],60,transform=ccrs.PlateCarree(),cmap='jet')
		cb=plt.colorbar()
		plt.quiver(lon[120:480][::10], lat[208:500][::10], u[208:500, 120:480][::10, ::10] / V[208:500, 120:480][::10, ::10], v[208:500, 120:480][::10, ::10] / V[208:500, 120:480][::10, ::10],linewidth=2)
		ax.coastlines()
		g = ax.gridlines(crs=ccrs.PlateCarree(), linestyle='-.', color='gray', draw_labels=True)
		g.ylabels_right = False
		g.xlabels_top = False
		plt.suptitle('velocidade'+'_'+str(current_data.strftime("%d%m%Y%H")))
		for i in range(len(latitudes)):
			latitudes[i]=lat[min(range(len(lat)), key=lambda x: abs(lat[x]-latitudes[i]))]
		for i in range(len(longitudes)):
			longitudes[i]=lon[min(range(len(lon)), key=lambda x: abs(lon[x]-longitudes[i]))]
		for i in range(len(latitudes)):
			coords_quadrado_x=[(lon[int((np.where(lon==longitudes[i]))[0])]),(lon[int((np.where(lon==longitudes[i]))[0])+15]),(lon[int((np.where(lon==longitudes[i]))[0])+15]),(lon[int((np.where(lon==longitudes[i]))[0])]),(lon[int((np.where(lon==longitudes[i]))[0])])]
			coords_quadrado_y=[(lat[int((np.where(lat==latitudes[i]))[0])]),(lat[int((np.where(lat==latitudes[i]))[0])]),(lat[int((np.where(lat==latitudes[i]))[0])-15]),(lat[int((np.where(lat==latitudes[i]))[0])-15]),lat[int((np.where(lat==latitudes[i]))[0])]]
			ax.plot(coords_quadrado_x,coords_quadrado_y, color="red")
		plt.tight_layout()
		plt.savefig('/home/joao/IC-Oceanografia_Fisica/plots/alta_frequencia/'+'vel'+str(current_data))
		print(fn)
		ds.close()
		plt.close()
	current_data = current_data + datetime.timedelta(days=1)



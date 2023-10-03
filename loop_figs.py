import netCDF4 as nc
import datetime
import time
import warnings
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
warnings.filterwarnings('ignore')

fig_output_dir='/home/joao/scripts/plots'

data_inicial = "01/01/2008"
data_final = "03/01/2008"


dir_cfs = '/home/joao/SAIDAS'

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])
current_data = data_inicial
print(current_data)

print(data_inicial)
print(data_inicial.strftime("%Y%m%d%H"))
print(data_final.strftime("%Y%m%d%H"))

print(dir_cfs+"/"+data_inicial.strftime("%Y%m%d%H"))
print(data_inicial.strftime("%Y"))
print(int(data_inicial.strftime("%d"))+1)

while (current_data<=data_final):
	print(" DIA: "+current_data.strftime("%d-%m-%Y"))
	fn=(dir_cfs+"/"+current_data.strftime("%Y%m%d%H")+'/archv.'+ str(data_inicial.strftime("%Y"))+'_00'+str(int(current_data.strftime("%d"))+1)+'_00_3zt.nc')
	ds=nc.Dataset(fn)
	lat=ds['Latitude'][:]
	lon=ds['Longitude'][:]
	#temperatura
	temp=ds['temperature'][:]
	temp=temp[0,0,:,:]
	ax = plt.axes(projection=ccrs.PlateCarree())
	plt.contourf(lon, lat,temp,60,transform=ccrs.PlateCarree(),cmap='RdBu_r')
	cb=plt.colorbar()
	ax.coastlines()
	g = ax.gridlines(crs=ccrs.PlateCarree(), linestyle='-.', color='gray', draw_labels=True)
	g.ylabels_right = False
	g.xlabels_top = False
	cb.ax.set_xlabel('Temperatura(°C)')
	plt.suptitle('Temperatura na superficie (°C) '+current_data.strftime("%d-%m-%Y"))
	plt.savefig('/home/joao/scripts/plots/temp'+'/'+'temperatura_sup'+current_data.strftime("%d-%m-%Y"))
	plt.close()
	#salinidade
	sal=ds['salinity'][:]
	sal=sal[0,0,:,:]
	ax = plt.axes(projection=ccrs.PlateCarree())
	plt.contourf(lon, lat,sal,60,transform=ccrs.PlateCarree(),cmap='RdBu_r')
	cb=plt.colorbar()
	ax.coastlines()
	g = ax.gridlines(crs=ccrs.PlateCarree(), linestyle='-.', color='gray', draw_labels=True)
	g.ylabels_right = False
	g.xlabels_top = False
	cb.ax.set_xlabel('Salinidade')
	plt.suptitle('salinidade na superficie psu '+current_data.strftime("%d-%m-%Y"))
	plt.savefig('/home/joao/scripts/plots/sal'+'/'+'sal_sup'+current_data.strftime("%d-%m-%Y"))
	plt.close()
	#densidade
	dens=ds['density'][:]
	dens=dens[0,0,:,:]
	ax = plt.axes(projection=ccrs.PlateCarree())
	plt.contourf(lon, lat,dens,60,transform=ccrs.PlateCarree(),cmap='RdBu_r')
	cb=plt.colorbar()
	ax.coastlines()
	g = ax.gridlines(crs=ccrs.PlateCarree(), linestyle='-.', color='gray', draw_labels=True)
	g.ylabels_right = False
	g.xlabels_top = False
	cb.ax.set_xlabel('Densidade')
	plt.suptitle('densidade na superficie sigma '+current_data.strftime("%d-%m-%Y"))
	plt.savefig('/home/joao/scripts/plots/dens'+'/'+'dens_sup'+current_data.strftime("%d-%m-%Y"))
	plt.close()
	#
	current_data = current_data + datetime.timedelta(days=1)


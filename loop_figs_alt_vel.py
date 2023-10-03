"plotando ssh e velocidade de corrente na superfície, usando um loop"

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

while (current_data<=data_final):
	print(" DIA: "+current_data.strftime("%d-%m-%Y"))
	fn=(dir_cfs+"/"+current_data.strftime("%Y%m%d%H")+'/archv.'+ str(data_inicial.strftime("%Y"))+'_00'+str(int(current_data.strftime("%d"))+1)+'_00_fsd.nc')
	ds=nc.Dataset(fn)
	lat=ds['Latitude'][:]
	lon=ds['Longitude'][:]
	ssh=ds['ssh'][:]
	fn2=(dir_cfs+"/"+current_data.strftime("%Y%m%d%H")+'/archv.'+ str(data_inicial.strftime("%Y"))+'_00'+str(int(current_data.strftime("%d"))+1)+'_00_3zu.nc')
	ds2=nc.Dataset(fn2)
	u=ds2['u']
	u=u[0,0,:,:]
	v=ds2['v']
	v=v[0,0,:,:]
	
	ax = plt.axes(projection=ccrs.PlateCarree())
	plt.contourf(lon, lat,ssh[0],60,transform=ccrs.PlateCarree(),cmap='RdBu_r')
	cb=plt.colorbar()
	plt.quiver(lon[::10],lat[::10],u[::10, ::10],v[::10, ::10],scale=0.75,scale_units="xy",units="inches",width=0.008,minlength=0.1,linewidth=2)
	ax.coastlines()
	g = ax.gridlines(crs=ccrs.PlateCarree(), linestyle='-.', color='gray', draw_labels=True)
	g.ylabels_right = False
	g.xlabels_top = False	
	cb.ax.set_xlabel('SSH')
	plt.suptitle('SSH '+'e '+'velocidade da corrente na superfície'+current_data.strftime("%d-%m-%Y"))
	plt.savefig('/home/joao/scripts/plots/ssh'+'/'+'ssh_corrente'+current_data.strftime("%d-%m-%Y"))
	current_data = current_data + datetime.timedelta(days=1)
	plt.close()

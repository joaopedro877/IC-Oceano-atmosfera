import netCDF4 as nc
import datetime
import time
import warnings
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
warnings.filterwarnings('ignore')
import numpy as np


data_inicial = "01/01/2008"
data_final = "03/01/2008"

dir_cfs = '/home/joao/SAIDAS'

data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])
current_data = data_inicial
while (current_data<=data_final):
	print(" DIA: "+current_data.strftime("%d-%m-%Y"))
	fn=(dir_cfs+"/"+current_data.strftime("%Y%m%d%H")+'/archv.'+ str(data_inicial.strftime("%Y"))+'_00'+str(int(current_data.strftime("%d"))+1)+'_00_3zt.nc')
	ds=nc.Dataset(fn)
	lat=ds['Latitude'][:]
	lon=ds['Longitude'][:]
	field=ds['temperature'][0,0,:,:]
	if 'sst' in locals():
		sst = np.ma.concatenate((sst,field[:,:,np.newaxis]),axis=2)
	else:
		sst = field[:,:,np.newaxis]
	current_data = current_data + datetime.timedelta(days=1)
	ds.close()

print(sst)
# tentar usar o np.concatenate e depois comparar as velocidades


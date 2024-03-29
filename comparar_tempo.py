import timeit


main_code_concatenate="""
import netCDF4 as nc
import datetime
import time
import numpy as np
import warnings
warnings.filterwarnings('ignore')

data_inicial = "01/01/2008"
data_final = "03/01/2008"
dir_cfs = '/home/joao/SAIDAS'




data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])
current_data = data_inicial
while (current_data<=data_final):
	#print(" DIA: "+current_data.strftime("%d-%m-%Y"))
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

#print(sst)
"""
main_code_stack="""
import netCDF4 as nc
import datetime
import time
import numpy as np
import warnings
warnings.filterwarnings('ignore')

data_inicial = "01/01/2008"
data_final = "03/01/2008"
dir_cfs = '/home/joao/SAIDAS'
data_inicial = datetime.datetime(*time.strptime(data_inicial,"%d/%m/%Y")[0:3])
data_final = datetime.datetime(*time.strptime(data_final,"%d/%m/%Y")[0:3])
current_data = data_inicial

temp=np.empty([701,553])
lista_temps=[]
counter=0

#concatenando os dados de sst
while (current_data<=data_final):
	#print(" DIA: "+current_data.strftime("%d-%m-%Y"))
	fn=(dir_cfs+"/"+current_data.strftime("%Y%m%d%H")+'/archv.'+ str(data_inicial.strftime("%Y"))+'_00'+str(int(current_data.strftime("%d"))+1)+'_00_3zt.nc')
	ds=nc.Dataset(fn)
	lat=ds['Latitude'][:]
	lon=ds['Longitude'][:]
	#temperatura
	#criando variaveis dinamicas
	globals()["temp"+str(counter)]=ds['temperature'][0,0,:,:]
	lista_temps.append((globals()["temp"+str(counter)]))
	counter =counter + 1
	ds.close()
	current_data = current_data + datetime.timedelta(days=1)

#print(lista_temps)

sst=np.ma.stack(lista_temps,2)
#print(sst)
"""
print("velocidade da função np.stack em segundos(1000 execuções do código) = "+str(timeit.timeit(stmt=main_code_stack,
                    setup='pass',
                    number=1000)))

print("velocidade da função np.concatenate em segundos (1000 execuções do código) = "+str(timeit.timeit(stmt=main_code_concatenate,
                    setup='pass',
                    number=1000)))


''' Resultados '''

''' A função stack é consideravelmente mais rápida. As diferenças de velocidade aumentam com o número de repetições do código -> para 1000 repetições
as velocidades foram : np.stack=13.789480146000642 e np.concatenate=1360.3775361049993(22 minutos)
Ainda preciso conferir se tudo está correto, mas parece ser um grande avanço'''

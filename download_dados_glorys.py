import copernicusmarine
from pprint import pprint
from datetime import datetime, timedelta

#Datas de início e término
start_date = datetime(2015, 1, 1)
end_date = datetime(2017, 12, 31)

# Pasta de destino 
output_directory = "/home/joao/IC-Oceanografia_Fisica/Saidas/glorys15_17"

# Loop para baixar os dados para cada dia
current_date = start_date
while current_date <= end_date:
    # Chame a função subset para salvar os dados
    get_result_daily = copernicusmarine.subset(
        dataset_id="cmems_mod_glo_phy_my_0.083deg_P1D-m",
        #definindo o tempo
		start_datetime=current_date,
        end_datetime=current_date,
        variables=['uo','vo'],
		#definindo os limites de area
		minimum_longitude=-60,
        maximum_longitude=-20,
        minimum_latitude=-36,
        maximum_latitude=-13,
        minimum_depth=0,
        maximum_depth=1,
		subset_method='nearest',
		#definindo nome do output
		output_filename='Glorys_u_v_sup'+str(current_date.strftime("%d-%m-%Y")),
        output_directory=output_directory,  
        #usuario e senha
		username='##############',
        password='############',
		#download automatico, sem perguntar antes
		force_download=True
    )

    pprint(f"List of saved files for {current_date.strftime('%Y-%m-%d')}: {get_result_daily}")
    
    # Incrementa para o próximo dia
    current_date += timedelta(days=1)


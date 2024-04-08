#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 12:26:17 2024

@author: joao
"""

import xarray as xr 

ds =xr.open_dataset('/home/joao/SAIDAS/CFS_FLUX01_fcst_2015123118.nc')
sliced_data=ds.sel(latitude=slice(-36,-13),longitude=slice(305,330))
sliced_data.to_netcdf('teste_slice.nc')
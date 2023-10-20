import xarray as xr
import os, glob
import imageio.v2 as imageio




png_dir = '/home/joao/scripts/plots/ssh/'
print(png_dir)

#

fixed_cb_images = glob.glob(png_dir + ('ssh_corrente*'))  
list_of_files=[]
[list_of_files.append(file) for file in fixed_cb_images]
list_of_files=(sorted(list_of_files))
#print(list_of_files)
#exit()
fix = [imageio.imread(file) for file in list_of_files]
imageio.mimsave(png_dir + '/ssh_c.gif', fix, duration = 700,loop=10)

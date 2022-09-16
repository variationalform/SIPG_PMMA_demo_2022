#!/anaconda3/envs/meshioenv/bin/python

# Author: Simon Shaw, 16 September 2022.
# Used for the paper: 
# A Priori Analysis of a Symmetric Interior Penalty Discontinuous
# Galerkin Finite Element Method for a Dynamic Linear Viscoelasticity Model
# by Yongseok Jang and Simon Shaw
# Submitted to: Computational Methods in Applied Mathematics 
  
# required: conda activate meshioenv

# you need to manually assign
# base_file_str, Nx, Ny in the code below

# 9 June 2022:
# conda env pde20191119, anaconda install of vtk and mayavi attempted, something failed.
# conda install -c conda-forge meshio
# conda install -c anaconda numpy - 

# create vtk anaconda environment called playground - pyhton 3.8
# conda activate playground
# anaconda install: numpy (+base,devel), vtk, mayavi
# ... and on command line
# conda install -c conda-forge meshio 
# numpy (?) failed with python 3.8

# base env, try vtk and mayavi installs via anaconda

# none of the above worked

# this worked!
# create new python 3.7 env called meshioenv.
# install numpy + base, devel
# conda activate meshioenv
# pip install meshio
# pip install numpy --upgrade
# matplotlib installed via anaconda navigator

# conda activate meshioenv
# ./runmenshio.py  # after fixing mesh dimensions and the seek path

import numpy as np
import meshio
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from xml.dom import minidom

#https://matplotlib.org/stable/api/matplotlib_configuration_api.html#matplotlib.set_loglevel
matplotlib.set_loglevel("error")  # to suppress PostScript transparency warnings

# get grid size
base_file_str = "../CMAM_DG/runtime/le/output_le_6/"
#base_file_str = "./output_ve_1/"
Nx=120; Ny=60; # these must be known
# last_pic_num=99
# Nx=128; Ny=64; last_pic_num=99 # these must be known

# first plot umid using umidsaved.npz and laTeX labels
npzfile = np.load(base_file_str+'umidsaved.npz')
umid = npzfile['umid']
T = umid[0,umid.shape[1]-1]
plt.clf()
plt.figure()
# don't use TeX because the headless server wont have it installed
plt.rcParams['text.usetex'] = True
plt.plot(umid[0,:], umid[1,:], linewidth=0.5, color='k')
plt.grid(color='k', linestyle='-', linewidth=0.2)
plt.xlabel('time, $t$ (seconds)')
plt.ylabel('displacement, $u_1(2,0.5)$ (metres)')
plt.ylim(-0.05, 0.1)
plt.xlim(0,T)
plt.draw()
plt.savefig(base_file_str+'umidplot.png', bbox_inches='tight', dpi=600)
plt.savefig(base_file_str+'umidplot.eps', bbox_inches='tight')

# now move on to plot the surfaces
plt.clf()
plt.figure()

mesh = meshio.read(base_file_str+"solution000000.vtu")
#print(mesh.points)
#print(mesh.cells)
#print(mesh.point_data) 

# get grid size
NumNodes, _ = mesh.points.shape
if NumNodes != (Nx+1)*(Ny+1):
  print('Error in mesh dimensions')
  exit(1)
else:
  print('Mesh dimensions are OK')

# get the base file for the times at which gfx were created
basefile = minidom.parse(base_file_str+'solution.pvd')
# extract the times fields
times = basefile.getElementsByTagName('DataSet')
# get a local array for times
ti = np.zeros(times.length)
# loop through them and assign to the local times array
#for elem in times:
  #print(elem.attributes['timestep'].value)
for c in range(0,times.length):
  ti[c] =  times[c].attributes['timestep'].value
#  print('ti[{0:3d}] = {1:8.6f}'.format(c,ti[c]) )

last_pic_num=times.length-1
print('frame numbers go from 0 to ', last_pic_num)

# create a meshgrid, x varies along rows, y down columns
X = np.zeros([1+Ny,1+Nx])
Y = np.zeros([1+Ny,1+Nx])

# populate first row of X
for c in range(0,1+Nx):
  X[:,c] = mesh.points[c,0]

# populate first column of Y
for c in range(0,1+Ny):
  Y[c,:] = mesh.points[c*(1+Nx),1]
  
#print(X[:2,:5])
#print(X[:2,-5:])
#print(Y[:5,:2])
#print(Y[-5:,:2])

#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
#fig.add_axes(ax)
#ax.view_init(20, 250)
#plt.figure(figsize=(10, 4))
#ax.set_box_aspect((2, 1, 0.0001))  

# now we have to populate Z for each frame
#for pic in range(0,1+40):
for pic in range(0,1+last_pic_num):
  if pic < 10:
    filestr = base_file_str + "solution00000"+str(pic)+".vtu"
  elif pic < 100:
    filestr = base_file_str + "solution0000"+str(pic)+".vtu"
  elif pic < 1000:
    filestr = base_file_str + "solution000"+str(pic)+".vtu"
  elif pic < 10000:
    filestr = base_file_str + "solution00"+str(pic)+".vtu"
  else:
    print("overflow in filestr, exiting...")
    exit(0)
#  print(filestr)
  
  mesh = meshio.read(filestr)

  Z = np.zeros([1+Ny,1+Nx])
  count = 0
  # use print(mesh.point_data) to get the f_14, f_17, ... 
#  u = mesh.point_data.get("f_14")
  u = mesh.point_data.get("f_17")
  for i in range(0,1+Ny):
    for j in range(0,1+Nx):
      Z[i,j] = u[count,0]    
      count += 1
      
  # https://matplotlib.org/stable/gallery/mplot3d/surface3d.html
  fig, ax = plt.subplots(subplot_kw={"projection": "3d"},figsize=(12, 9))
  fig.add_axes(ax)
#  ax.view_init(20, 250)
  ax.view_init(20, 230)
#  ax.view_init(20, 200)
  ax.set_box_aspect((2, 1, 1))  
  
  # Make data.
  #X = np.arange(-5, 5, 0.25)
  #Y = np.arange(-5, 5, 0.25)
  #X, Y = np.meshgrid(X, Y)
  #Z = X + 3*Y
  #R = np.sqrt(X**2 + Y**2)
  #Z = np.sin(R)
  
  # Plot the surface.
  base_dim = 0.05
  cont_gap = 0.001
  surf = ax.plot_surface(X, Y, Z, cmap=cm.twilight_shifted, vmin=-base_dim, vmax=base_dim,
                         linewidth=0, antialiased=False, lw=1)
  clev = cont_gap*np.arange(-100,100)
  # don't include zero
  clev[100:] += cont_gap
  ax.contour(X, Y, Z, clev, colors="k", linestyles="solid",offset=-base_dim)

  plt.xlabel('$x$',fontsize=20)
  plt.ylabel('$y$',fontsize=20)
  plt.title('time: {0:8.6f} secs'.format(ti[pic]),fontsize=30)
  # Customize the z axis.
  ax.set_zlim(-base_dim, base_dim)
  #ax.zaxis.set_major_locator(LinearLocator(10))
  # A StrMethodFormatter is used automatically
  #ax.zaxis.set_major_formatter('{x:.02f}')
  
  # Add a color bar which maps values to colors.
  cb = fig.colorbar(surf, shrink=0.5, aspect=5, ax=ax)
#  matplotlib.colors.Normalize(vmin=-1, vmax=1)
  # cb = fig.colorbar(im2, ax=(ax1, ax2), orientation='vertical')
  #cb.set_clim(-0.01,0.01)

#  ax = plt.axes()
#  ax.plot([0,2],[0,10])
#  ax.set_aspect('equal','box')  
  plt.draw()
  plt.pause(0.0001)
  plt.savefig(base_file_str + 'surf_u1_t={0:8.6f}.png'.format(ti[pic]),format='png')
  plt.savefig(base_file_str + 'surf_u1_t={0:8.6f}.eps'.format(ti[pic]),format='eps')
  plt.close()
#  fig.canvas.flush_events()

exit()

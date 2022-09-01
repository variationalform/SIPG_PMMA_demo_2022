#!/usr/bin/python
# the crash-bang line above changes with the conda environment. Best to run with 'python ...'

# mac
#(pde20191119) ... $ rm -rf output && /anaconda3/envs/pde20191119/bin/python simulator.py 
# Desktop
#conda activate fenicsproject
#(fenicsproject) ... ... $ rm -rf output && /home/icsf/icsrsss/anaconda3/envs/fenicsproject/bin/python simulator.py 
# server (11,12,docker) - no X server, matplotlib fails.
#conda activate blockdata
#rm -rf output && /home/icsrsss/anaconda3/envs/blockdata/bin/python simulator.py

from __future__ import print_function
from dolfin import *
import numpy as np
import math
import getopt, sys
import matplotlib.pyplot as plt
from time import time, ctime

parameters['form_compiler']['cpp_optimize'] = True
parameters['form_compiler']['optimize'] = True
parameters["ghost_mode"] = "shared_facet"
# set_log_active(False)

start = time()
PI=np.pi
tol = 1E-15
traction_choice = 4
if traction_choice == 0:
  T = 1.0e-1     # total simulation time
  amp  = 5.0e6 / 1e9 # amplitude of traction g1 - scaled by 1e9 for the solver
elif traction_choice == 1:
  T = 1.0e-1     # total simulation time
  amp  = 5.0e6 / 1e9 # amplitude of traction g1 - scaled by 1e9 for the solver
elif traction_choice == 2:
  T = 1.0e-1     # total simulation time
  amp  = 10.0e6 / 1e9 # amplitude of traction g1 - scaled by 1e9 for the solver
elif traction_choice == 3:
  T = 2.0e-1     # total simulation time
  amp  = 100.0e6 / 1e9 # amplitude of traction g1 - scaled by 1e9 for the solver
else:
  T = 0.6/2 #1.0e-2     # total simulation time
  amp  = 50.0e6 / 1e9 # amplitude of traction g1 - scaled by 1e9 for the solver
freq = 10 #1000    # sine frequency in traction g1
Tp = 1.0/freq  # periodic time

Nx = 20            # subdivisions in x direction
Ny = 10            # subdivisions in y direction
Nt = 40000         # the number of time steps
vwan = 8000        # variable without a name - controld output frequency  
gfx_count = 0
dt = T/Nt          # time step
viscoelastic = True

# interactive - ask for mesh dims and timesteps
print('enter Nx, Ny, Nt, and viscoelastic true/false - 1/0 - flag')
Nx = int(input())
Ny = int(input())
Nt = int(input())
vwan = int(Nt/5)
dt = T/Nt

vflag = int(input())
if vflag == 0: 
  viscoelastic = False
else: 
  viscoelastic = True

print('Inputs: Nx = ', Nx,', Ny = ', Ny,', Nt = ', Nt, 'vwan = ',  vwan)

# problem data for PMMA
Nphi=11 # the number of internal variables
if viscoelastic:
  print('Using viscoelasticity for PMMA')
  varphi=np.array([0.001000236663139,\
                   0.086627639575435,\
                   0.126369185566228,\
                   0.247379960437068,\
                   0.268813603218619,\
                   0.173255279150871,\
                   0.069659339040041,\
                   0.018307903209241,\
                   0.006162172299696,\
                   0.001643245946586,\
                   0.000352762037446,\
                   0.000428672855631])
else:
  print('Using elastic data from PMMA')
  varphi = np.zeros([12])
  varphi[0] = 1.0
  
nu=Constant(0.35)               # Poisson's ratio
E=Constant(2.23947)             # stiffness modulus, scaled by 1e9
Lame1=nu*E/((1+nu)*(1-2*nu))    # first Lame parameter
Lame2=0.5*E/(1+nu)              # second Lame parameter
rho=Constant(1190.0e-9)         # density, scaled by 1e9
lmbda=float(Lame1)
mu = float(Lame2)
# calculate some possible overrides given the frequency and modulus
print('Current Nt = ', Nt, ' and typical wave speed sqrt(E/rho) = ', math.sqrt(E/rho), 'm/s')
print('Over [0,T]=[0,',T,'] such a wave will move ', T*math.sqrt(E/rho), 'm')
print('and during dt={0:13.6e} such a wave will move {1:13.6e} m'.format(dt, dt*math.sqrt(E/rho)) )
print('\nlambda = ', lmbda, ',  \t mu = ', mu)
Cd = math.sqrt((lmbda+2*mu)/rho); Cs = math.sqrt(mu/rho)
print('Cd = ', Cd, ',  \t Cs = ', Cs)
print('wavelengths for freq =  ', freq, ': ', Cd/freq, ' and ', Cs/freq)
print('0.1m travel times: ', 0.1/Cd, ' and ', 0.1/Cs)
print('Distance travelled in one time step: ', Cd*dt,' and ', Cs*dt,'\n')
# suppose we want solres intervals per sine wave, with each sine wave having
# periodic time 1/(freq*(1-x[1])). We will need to consider the max frequency
# which is when y=x[1]=0. Hence the effective periodic time is 1/freq and in
# the solution interval [0,T] there are T*freq such periods. To have solres 
# intervals in each will require Nt = solres*T*freq...
solres=2000
print('with solres = {0:6d}: recommend Nt = {1:6d}\n'.format(solres, int(solres*T*freq)))

'''
altering = input('Use these? Y/N    ')
if altering == 'Y' or altering == 'y':
  print('altering...')
  Nt = int(solres*T/freq);
  minres = 1; # unclear what we wanted to achieve here...
'''

tau=np.zeros(Nphi+1) # for the sake of convenient indexing, we additionally put tau0=0
for i in range(0,Nphi):
    tau[i+1]=2.0*pow(10,(i-2))

# Define penalty parameters
alpha = 10.0
beta =  1.0

ux = Expression(("0.0","0.0"), tn=0, degree=5)

mesh = RectangleMesh(Point(0.0, 0.0), Point(2.0, 1.0), Nx, Ny,"left")
# plot(mesh); plt.show()

V = VectorFunctionSpace(mesh, 'DG', 2)
dof=V.dim()

# define the boundary partition
boundary_parts =MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)

# nonzero traction Neumann bc
class GammaNeumann(SubDomain):
    def inside(self, x, on_boundary):
       return on_boundary and near(x[0],2.0,tol)
Gamma_Neumann = GammaNeumann()
Gamma_Neumann.mark(boundary_parts, 1)

# homogeneous Dirichlet bc 
class GammaDirichlet(SubDomain):
    def inside(self, x, on_boundary):
       Gd = near(x[0],0.0,tol)
       return on_boundary and Gd
Gamma_Dirichlet = GammaDirichlet()
Gamma_Dirichlet.mark(boundary_parts, 2)

ds = ds(subdomain_data = boundary_parts)
dx=Measure('dx')
n = FacetNormal(mesh)
h = CellDiameter(mesh)
h_avg = (h('+') + h('-'))/2

# Strain tensor
def epsilon(v):
    Dv=nabla_grad(v)
    return 0.5*(Dv+Dv.T)

# Stress tensor (elasticity part)
def sigma(v):
    return Lame1*tr(epsilon(v))*Identity(2) +2*Lame2*epsilon(v)

u0 = interpolate(ux, V)
w0 = interpolate(ux, V)

if traction_choice == 0:
  g=Expression(( "0.5*amp*( tn    < 2*0.5/(freq*(1-0.75*x[1])) )*pow(sin(2*freq*pi*(1-0.75*x[1])*  tn   ), 2)       \
               +0.5*amp*( tn+dt < 2*0.5/(freq*(1-0.75*x[1])) )*pow(sin(2*freq*pi*(1-0.75*x[1])*(tn+dt)), 2)", "0"),
               amp=amp,dt=dt,tn=0,degree=5,pi=PI,freq=freq)
elif traction_choice == 1:
  g=Expression(( "0.5*amp*( exp(-pow(x[1],  2)/0.01)*(1+tanh((tn   -0.025)/0.005) )      \
                           +exp(-pow(x[1]-1,2)/0.05)*(1+tanh((tn   -0.025)/0.001) ) )    \
                + 0.5*amp*( exp(-pow(x[1],  2)/0.01)*(1+tanh((tn+dt-0.025)/0.005) )      \
                           +exp(-pow(x[1]-1,2)/0.05)*(1+tanh((tn+dt-0.025)/0.001) ) )    \
                                                       ", "0"),
                  amp=amp,dt=dt,tn=0,degree=5,pi=PI,freq=freq)
elif traction_choice == 2:
  g=Expression(( "0.5*amp*( exp(-(pow(x[1]  ,2)+pow((  tn -0.025)*10,2) )/0.001)            \
                          + exp(-(pow(x[1]-1,2)+pow((  tn -0.025)*10,2) )/0.01))            \
                + 0.5*amp*( exp(-(pow(x[1]  ,2)+pow((tn+dt-0.025)*10,2) )/0.001)            \
                          + exp(-(pow(x[1]-1,2)+pow((tn+dt-0.025)*10,2) )/0.01))            \
                                                       ", "0"),
                  amp=amp,dt=dt,tn=0,degree=5,pi=PI,freq=freq)
elif traction_choice == 3:
  g=Expression(( "0.5*amp*( exp(-(pow(x[1]-0.75,2)+pow(  tn   *10,2) )/0.0001)*pow(sin(2*pi*freq*  tn   ),2) )            \
                + 0.5*amp*( exp(-(pow(x[1]-0.75,2)+pow((tn+dt)*10,2) )/0.0001)*pow(sin(2*pi*freq*(tn+dt)),2) )            \
                                                       ", "0"),
                  amp=amp,dt=dt,tn=0,degree=5,pi=PI,freq=freq)
else:
#  EX01
#  g=Expression(( "0.5*amp+0.5*amp", "0"), amp=amp,dt=dt,tn=0,degree=5,pi=PI,freq=freq)
#  EX02
  g=Expression(( "0.5*amp*( tn  < 0.01 )+0.5*amp*( tn+dt < 0.01 )", "0"),
               amp=amp,dt=dt,tn=0,degree=5,pi=PI,freq=freq)

uh = Function(V)   # the unknown at a new time level
wh = Function(V)

# bilinear form for the solver
u, v = TrialFunction(V), TestFunction(V)
mass = rho*inner(u,v)*dx
stiffness=inner(sigma(u),epsilon(v))*dx- inner(avg(sigma(u)), outer(v('+'),n('+'))+outer(v('-'),n('-')))*dS \
         -inner(avg(sigma(v)), outer(u('+'),n('+'))+outer(u('-'),n('-')))*dS \
         + alpha/(h_avg**beta)*inner(jump(u), jump(v))*dS \
         - inner(sigma(u), outer(v,n))*ds(2) \
         - inner(outer(u,n), sigma(v))*ds(2) \
         + alpha/(h**beta)*dot(u,v)*ds(2)    
jump_penalty =  alpha/(h_avg**beta)*dot(jump(u), jump(v))*dS  + alpha/(h**beta)*dot(u,v)*ds(2)  
        
# linear form for the right hand side and internal variables
L=dot(g,v)*ds(1)
# assemble the system matrix once and for all
M = assemble(mass)
A = assemble(stiffness)
J = assemble(jump_penalty)

coeff_internal=0.0
for i in range(1,Nphi+1):
    coeff_internal += varphi[i]/(tau[i]/dt+0.5)/2
CurlA=0.5*(A-coeff_internal*A)
B=(2.0/dt/dt)*M+CurlA+(1.0/dt)*J

# assemble only once, before the time stepping
b = None; b2= None; b3= None

# record time and displacement at (x,y)=(2,0.5) - assume zero IC
umid = np.zeros((2,Nt+1))

Psi=np.zeros((dof,Nphi))
vtkfile = File('output/solution.pvd')
then = time()
for n in range(0, Nt):
    tn = n*dt;  g.tn = tn; 
    b = assemble(L, tensor=b)
    b2=2.0/dt*M*w0.vector().get_local()+(2.0/dt/dt*M-CurlA+1.0/dt*J)*u0.vector().get_local()
    for q in range(0,Nphi):
        b2+=A*2.0*tau[q+1]/(2.0*tau[q+1]+dt)*Psi[:,q]
    b.add_local(b2)

#   lots of choice, crude tests suggest mumps is the fastest
#    solve(B, uh.vector(), b,'minres')  #  solve(B, uh.vector(), b,'gmres')    #  solve(B, uh.vector(), b,'bicgstab')  
#    solve(B, uh.vector(), b,'petsc')   #  solve(B, uh.vector(), b,'umfpack')  #  solve(B, uh.vector(), b,'lu')
#    solve(B, uh.vector(), b, 'cg', 'hypre_amg')
    solve(B, uh.vector(), b, 'mumps')
    wh.vector()[:]=2.0/dt*(uh.vector().get_local()-u0.vector().get_local())-w0.vector().get_local()
    for q in range(0,Nphi):
        Psi[:,q]=2.0*dt/(2*tau[q+1]+dt)*((2*tau[q+1]-dt)/2/dt*Psi[:,q]+varphi[q+1]/2.0*(uh.vector().get_local()+u0.vector().get_local()))            

    # throttle the amount of output - this is a bit buggy
    if (tn < T/100-dt and tn+dt >= T/100-dt) or (n+1) % vwan == 0 or n == Nt-1:
      tnp1 = tn+dt
      vtkfile << (uh, tnp1)
      u1,_=uh.split()
      now = time()
      print('{0:13.6e}  {1:3d}%  {2:13.6e}  {3:13.6e}  {4:13.6e}  Dt={5:13.6e} s'.format(
              tnp1, int(round(100*tnp1/T)), norm(u1, 'L2'), u1(0,0.5), u1(2,0.5), now-then ))
#              tnp1, int(round(100*tnp1/T)), norm(u1, 'L2'), u1(0,0.5), u1(2,0.5), now-then ), flush=True)
      sys.stdout.flush()
      then = now
      plot(u1)
      gfxstr = 'u1_'
      if gfx_count < 10:
        gfxstr = gfxstr+'0000'+str(gfx_count)
      elif gfx_count < 100:
        gfxstr = gfxstr+'000'+str(gfx_count)
      elif gfx_count < 1000:
        gfxstr = gfxstr+'00'+str(gfx_count)
      elif gfx_count < 10000:
        gfxstr = gfxstr+'0'+str(gfx_count)
      else:
        gfxstr = gfxstr+str(gfx_count)
      gfx_count = 1+gfx_count
      plt.savefig('./output/'+gfxstr+'.png', bbox_inches='tight')
      plt.savefig('./output/'+gfxstr+'.eps', bbox_inches='tight')
#      plt.draw(); plt.pause(0.01)
    # capture midpoint displacement 
    u1,_=uh.split(); umid[0,n+1] = tn+dt; umid[1,n+1] = u1(2,0.5)
#   update old terms
    w0.assign(wh);u0.assign(uh);

plt.clf()
plt.figure()
# don't use TeX because the headless server wont have it installed
#plt.rcParams['text.usetex'] = True
plt.plot(umid[0,:], umid[1,:], linewidth=0.5, color='k')
plt.grid(color='k', linestyle='-', linewidth=0.2)
plt.xlabel('time, t (seconds)')
#plt.xlabel('time, $t$ (seconds)')
plt.ylabel('displacement, u_1(2,0.5) (metres)')
#plt.ylabel('displacement, $u_1(2,0.5)$ (metres)')
plt.ylim(-0.05, 0.1)
plt.xlim(0,T)
plt.draw()
# suffix x to denote non-LaTeX plots
#plt.savefig('./output/umidplot.png', bbox_inches='tight', dpi=600)
#plt.savefig('./output/umidplot.eps', bbox_inches='tight')
plt.savefig('./output/umidplotx.png', bbox_inches='tight', dpi=600)
plt.savefig('./output/umidplotx.eps', bbox_inches='tight')
# save umid for offline plotting (with e.g. LaTeX labels)
np.savez('./output/umidsaved.npz', umid=umid)
#plt.show(); plt.pause(1); plt.close('all')
end = time()
print('Elapsed time: {0:13.6e}s  {1:13.6e}m  {2:13.6e}h'.format(
              end-start, (end-start)/60.0, (end-start)/3600.0))

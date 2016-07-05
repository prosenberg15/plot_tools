import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
import pylab
from pylab import *
import re
from numpy import mod
import matplotlib.pyplot as plt
from matplotlib import rc
import sys
import scipy.interpolate
import fileinput
import math
from scipy import interpolate
import matplotlib.transforms as transforms
import matplotlib.axes
from plot_tools import *

rc('text', usetex=True)
rc('font', family='serif')

Nl=[12,12]
Ntot=Nl[0]*Nl[1]

shift=np.zeros(Ntot)
shift.fill(1)

# load psi_r_normed.dat (normalized real-space wfn)
psidata=np.loadtxt(sys.argv[1],dtype=np.float64)
latt=np.loadtxt(sys.argv[2])
# load x, y coords from lattice-label.dat
psi_x=latt[:,1]-shift
psi_y=latt[:,2]-shift

# grab psi_dd, psi_uu, psi_s
# compute |psi_dd|, |psi_uu|, |psi_s| 
##############################################################
# psi_uu
psi_r=psidata[:Ntot,0]+1j*psidata[:Ntot,1]
# |psi_uu|
mod_psi_r=np.sqrt((psi_r[:].real)**2+(psi_r[:].imag)**2)
##############################################################
# psi_dd
psi_r_dd=psidata[Ntot:2*Ntot,0]+1j*psidata[Ntot:2*Ntot,1]
# |psi_dd|
mod_psi_r_dd=np.sqrt((psi_r_dd[:].real)**2+(psi_r_dd[:].imag)**2)
##############################################################
# psi_s
psi_r_s=psidata[2*Ntot:,0]+1j*psidata[2*Ntot:,1]
# |psi_s|
mod_psi_r_s=np.sqrt((psi_r_s[:].real)**2+(psi_r_s[:].imag)**2)

psi_x=shift_origin_list(psi_x,Nl[0])
psi_y=shift_origin_list(psi_y,Nl[0])

"""
for i in range(len(psi_x)):
    if psi_x[i]>=Nl[0]/2:
        psi_x[i]=psi_x[i]-Nl[0]
    if psi_y[i]>=Nl[1]/2:
        psi_y[i]=psi_y[i]-Nl[1]
"""

#print mod_psi_r_dd.max()
# grab 1D slice
# slice defined by starting point (x0,y0)
# and direction vector (xsdir,ysdir)
x0=-2
y0=0
xsdir=1
ysdir=0

xsexact,ysexact,qsexact,qsexacterr=grab_slice(x0,y0,xsdir,ysdir,psi_x,psi_y,mod_psi_r_dd,[0]*len(psi_r))
#xsexact=np.array(xsexact)
#xsexact=np.sqrt(2.0)*xsexact

# set up figure environment
fig = plt.figure()

# for contour plot
# Set up a regular grid of interpolation points
#xi, yi = np.linspace(psi_x.min(), psi_x.max(), 100), np.linspace(psi_y.min(), psi_y.max(), 100)
#xi, yi = np.meshgrid(xi, yi)

xi, yi = np.ogrid[psi_x.min():psi_x.max():200j, psi_y.min():psi_y.max():200j]

mod_psi_r_dd_2D=reshape_1D_to_2D(psi_x,psi_y,mod_psi_r_dd)
#mod_psi_r_s_2D=reshape_1D_to_2D(psi_x,psi_y,mod_psi_r_s)
# shift origin
mod_psi_r_dd_2D=shift_origin_array(mod_psi_r_dd_2D,Nl)
#mod_psi_r_s_2D=shift_origin_array(mod_psi_r_s_2D,Nl)
# interpolate
mod_psi_r_dd_2D_i=splineinterp(psi_x,psi_y,mod_psi_r_dd_2D)
#mod_psi_r_s_2D_i=splineinterp(psi_x,psi_y,mod_psi_r_s_2D)

# scale 0 point
#mod_psi_r_s[0]=mod_psi_r_s[0]/10.0
#mod_psi_r_s_i=scipy.interpolate.griddata((psi_x, psi_y), mod_psi_r_s, (xi, yi), method='cubic')

#for 3D plot uncomment below
###############################
ax=fig.gca(projection = '3d')

# fix aspect ratio if lattice isn't square
# otherwise use ratio of 1.0
if(abs(Nl[0]-Nl[1])>1e-9):
    ax.pbaspect = [0.75, 1.5,0.5]
elif(abs(Nl[0]-Nl[1])<1e-9):
    ax.pbaspect = [1.0,1.0,0.5]

# axis labels and limits
ax.set_xlim(-Nl[0]/2,Nl[0]/2-1)
ax.set_ylim(-Nl[1]/2,Nl[1]/2-1)
ax.set_xlabel('x')
ax.set_ylabel('y')

surf = ax.plot_surface(xi, yi, mod_psi_r_dd_2D_i(xi,yi), rstride=1, cstride=1, cmap=cm.gnuplot,
        linewidth=0, antialiased=False)


#,levels=np.arange(-0.0001,0.48,0.00005)
fig.colorbar(surf,shrink=0.75,pad=0.01)

show()
#fig.savefig('psi_r_triplet_Lx6_Ly24_N144_l0.5_U-2.5.pdf',bbox_inches='tight')


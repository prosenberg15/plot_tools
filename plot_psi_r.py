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

Nl=np.array([12,12])
Np=144
Ntot=Nl[0]*Nl[1]

shift=np.zeros(Ntot)
shift.fill(1)


if ('didj' in sys.argv[1]):
    plot_flag='didj'
    arg_shift=1
elif ('ninj' in sys.argv[1]):
    plot_flag='ninj'
    arg_shift=1
elif ('psi_r' in sys.argv[1] and not ('dat' in sys.argv[1])):
    plot_flag='psi_r'
    arg_shift=1
else: 
    print 'select an observable to plot'
    plot_flag=str(raw_input("didj, ninj, nk, psi_r? : "))
    arg_shift=0

if('ninj' in plot_flag):
    data=np.loadtxt(sys.argv[1+arg_shift],dtype=np.float64)
    # nu*nd                       
    nui_ndj=np.array(data[Nl[0]*Nl[1]:,1])
    nui_ndj_err=np.array(data[Nl[0]*Nl[1]:,3])
    # nu*nu
    nui_nuj=np.array(data[:Nl[0]*Nl[1],1])
    nui_nuj_err=np.array(data[:Nl[0]*Nl[1],3])
    # nu*nd + nu*nu
    ninj=nui_ndj+nui_nuj
    ninjerr=nui_ndj_err+nui_nuj_err

    lattlabel=np.loadtxt(sys.argv[2+arg_shift])
    # shift, so lattice site 1 has coords (x,y)=(0,0) (not (1,1))
    shift=np.zeros(len(lattlabel[:,0]))
    shift.fill(1)
    ninj_x=lattlabel[:,1]-shift[:]
    ninj_y=lattlabel[:,2]-shift[:]
    
    # shift x,y coords (currently stored as lists)
    ninj_x=shift_origin_list(ninj_x,Nl[0])
    ninj_y=shift_origin_list(ninj_y,Nl[0])

    xi, yi = np.ogrid[ninj_x.min():ninj_x.max():200j, ninj_y.min():ninj_y.max():200j]

    corr_flag = str(raw_input("nu*nu, nu*nd, ntot*ntot? : "))
    if('nu*nu' in corr_flag):
        # nu*nu
        # reshape list --> 2D array
        nui_nuj_2D=reshape_1D_to_2D(ninj_x,ninj_y,nui_nuj)
        # shift origin
        nui_nuj_2D=shift_origin_array(nui_nuj_2D,Nl)
        # interpolate
        nui_nuj_2D_i=splineinterp_2D(ninj_x,ninj_y,nui_nuj_2D,3,3,0)
        plot_func=nui_nuj_2D_i
    elif('nu*nd' in corr_flag):
        # nu*nd
        nui_ndj_2D=reshape_1D_to_2D(ninj_x,ninj_y,nui_ndj)
        # shift origin
        nui_ndj_2D=shift_origin_array(nui_ndj_2D,Nl)
        # interpolate
        nui_ndj_2D_i=splineinterp_2D(ninj_x,ninj_y,nui_ndj_2D,3,3,0)
        plot_func=nui_ndj_2D_i
    elif('ntot*ntot' in corr_flag):
        # ntot*ntot
        ninj_2D=reshape_1D_to_2D(ninj_x,ninj_y,ninj)
        # shift origin
        ninj_2D=shift_origin_array(ninj_2D,Nl)
        # interpolate
        ninj_2D_i=splineinterp_2D(ninj_x,ninj_y,ninj_2D,3,3,0)
        plot_func=ninj_2D_i

elif('didj' in plot_flag):
    #load didj data
    data=np.loadtxt(sys.argv[1+arg_shift],dtype=np.float64)
    didj=data[:,1]
    didjerr=data[:,3]
    # scale origin
    #didj[0]=0.01#didj[0]/10.0
    # normalize
    #didj=didj[:]*(Nl[0]*Nl[1])**2/(Np*(Np-1))
    #didjerr=didjerr[:]*(Nl[0]*Nl[1])**2/(Np*(Np-1))
    # load x, y coords from lattice-label.dat
    lattlabel=np.loadtxt(sys.argv[2+arg_shift])
    # shift, so lattice site 1 has coords (x,y)=(0,0) (not (1,1))
    shift=np.zeros(len(lattlabel[:,0]))
    shift.fill(1)
    didj_x=lattlabel[:,1]-shift[:]
    didj_y=lattlabel[:,2]-shift[:]
    
    # shift x,y coords (currently stored as lists)
    didj_x=shift_origin_list(didj_x,Nl[0])
    didj_y=shift_origin_list(didj_y,Nl[0])

    xi, yi = np.ogrid[didj_x.min():didj_x.max():200j, didj_y.min():didj_y.max():200j]

    # reshape list --> 2D array
    didj_2D=reshape_1D_to_2D(didj_x,didj_y,didj)
    # shift origin
    didj_2D=shift_origin_array(didj_2D,Nl)
    # interpolate
    didj_2D_i=splineinterp_2D(didj_x,didj_y,didj_2D,3,3,0)
    plot_func=didj_2D_i


elif('psi_r' in plot_flag):
    # load psi_r_normed.dat (normalized real-space wfn)
    psidata=np.loadtxt(sys.argv[1+arg_shift],dtype=np.float64)
    # load x, y coords from lattice-label.dat
    latt=np.loadtxt(sys.argv[2+arg_shift])
    # shift, so lattice site 1 has coords (x,y)=(0,0) (not (1,1))
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
    dd=psidata[Ntot:2*Ntot,0]+1j*psidata[Ntot:2*Ntot,1]
    # |psi_dd|
    mod_dd=np.sqrt((dd[:].real)**2+(dd[:].imag)**2)
    dd_re=dd[:].real
    dd_im=dd[:].imag
    ##############################################################
    # psi_s
    psi_r_s=psidata[2*Ntot:,0]+1j*psidata[2*Ntot:,1]
    # |psi_s|
    mod_psi_r_s=np.sqrt((psi_r_s[:].real)**2+(psi_r_s[:].imag)**2)
    mod_psi_r_s_re=psi_r_s[:].real
    mod_psi_r_s_im=psi_r_s[:].imag
    ##############################################################
    # |psi_tot|
    mod_psi_r_tot=mod_psi_r_s+mod_dd+mod_psi_r
    mod_psi_r_tot_re=psi_r_s[:].real+dd[:].real+psi_r[:].real
    mod_psi_r_tot_im=psi_r_s[:].imag+dd[:].imag+psi_r[:].imag

    # shift x,y coords (currently stored as lists), so origin is at center, not edge of plot 
    psi_x=shift_origin_list(psi_x,Nl[0])
    psi_y=shift_origin_list(psi_y,Nl[1])

    """
    psi_x_trunc=[]
    psi_y_trunc=[]
    mod_psi_r_dd_trunc=[]
    for ind in range(len(psi_x)):
        if(abs(psi_x[ind]+Nl[0]/2)>1e-8 and abs(psi_y[ind]+Nl[0]/2)>1e-8):
            psi_x_trunc.append(psi_x[ind])
            psi_y_trunc.append(psi_y[ind])
            mod_psi_r_dd_trunc.append(mod_psi_r_dd[ind])
    for ind in range(len(psi_x)):
        if(psi_x[ind]>0.5 and psi_y[ind]>0.5):
            psi_x_trunc.append(psi_x[ind])
            psi_y_trunc.append(psi_y[ind])
            mod_psi_r_dd_trunc.append(mod_psi_r_dd[ind])
    psi_x_trunc=np.array(psi_x_trunc)
    psi_y_trunc=np.array(psi_y_trunc)
    mod_psi_r_dd_trunc=np.array(mod_psi_r_dd_trunc)
    """

    xi, yi = np.ogrid[psi_x.min():psi_x.max():200j, psi_y.min():psi_y.max():200j]
    #xi_trunc, yi_trunc = np.ogrid[psi_x_trunc.min():psi_x_trunc.max():200j, psi_y_trunc.min():psi_y_trunc.max():200j]

    if('dd' in sys.argv[1]):
        # modulus wfn
        # reshape list --> 2D array
        dd_2D=reshape_1D_to_2D(psi_x,psi_y,mod_dd)
        # shift origin
        dd_2D=shift_origin_array(dd_2D,Nl)
        # interpolate
        dd_2D_i=splineinterp_2D(psi_x,psi_y,dd_2D,3,3,0)
        plot_func=dd_2D_i

        # real part
        # reshape list --> 2D array
        dd_2D_re=reshape_1D_to_2D(psi_x,psi_y,dd_re)
        # shift origin
        dd_2D_re=shift_origin_array(dd_2D_re,Nl)
        # interpolate
        dd_2D_re_i=splineinterp_2D(psi_x,psi_y,dd_2D_re,3,3,0)
        plot_func_re=dd_2D_re_i

        # imag part
        # reshape list --> 2D array
        dd_2D_im=reshape_1D_to_2D(psi_x,psi_y,dd_im)
        # shift origin
        dd_2D_im=shift_origin_array(dd_2D_im,Nl)
        # interpolate
        dd_2D_im_i=splineinterp_2D(psi_x,psi_y,dd_2D_im,3,3,0)
        plot_func_im=dd_2D_im_i

        """
        plot_func_re=perform_interp(psi_x,psi_y,mod_psi_r_dd_re,3,3,0,True,Nl)
        plot_func_im=perform_interp(psi_x,psi_y,mod_psi_r_dd_im,3,3,0,True,Nl)
        plot_func=perform_interp(psi_x,psi_y,mod_psi_r_dd,3,3,0,True,Nl)
        """

        """
        # reshape list --> 2D array
        mod_psi_r_dd_2D_trunc=reshape_1D_to_2D(psi_x_trunc,psi_y_trunc,mod_psi_r_dd_trunc)
        #psi_x_trunc=psi_x_trunc/(1.0*Nl[0])
        #psi_y_trunc=psi_y_trunc/(1.0*Nl[1])
        # shift origin
        mod_psi_r_dd_2D_trunc=shift_origin_array(mod_psi_r_dd_2D_trunc,Nl-1)
        # interpolate
        #mod_psi_r_dd_2D_trunc_i=splineinterp_2D(psi_x_trunc,psi_y_trunc,mod_psi_r_dd_2D_trunc,2,2,0.005)
        #mod_psi_r_dd_2D_trunc_i=interpolate.Rbf(psi_x_trunc,psi_y_trunc,mod_psi_r_dd_trunc,function='linear',smooth=0.0005)
        #plot_func_trunc=mod_psi_r_dd_2D_trunc_i
        """

    elif('sing' in sys.argv[1] or '_s' in sys.argv[1]):
        # modulus wfn
        # reshape list --> 2D array
        mod_psi_r_s_2D=reshape_1D_to_2D(psi_x,psi_y,mod_psi_r_s)
        # shift origin
        mod_psi_r_s_2D=shift_origin_array(mod_psi_r_s_2D,Nl)
        # interpolate
        mod_psi_r_s_2D_i=splineinterp_2D(psi_x,psi_y,mod_psi_r_s_2D,3,3,0)
        #mod_psi_r_s_2D_i=interpolate.Rbf(psi_x,psi_y,mod_psi_r_s,function='gaussian',smooth=0.02)
        plot_func=mod_psi_r_s_2D_i

        # real part
        # reshape list --> 2D array
        mod_psi_r_s_2D_re=reshape_1D_to_2D(psi_x,psi_y,mod_psi_r_s_re)
        # shift origin
        mod_psi_r_s_2D_re=shift_origin_array(mod_psi_r_s_2D_re,Nl)
        # interpolate
        mod_psi_r_s_2D_re_i=splineinterp_2D(psi_x,psi_y,mod_psi_r_s_2D_re,3,3,0)
        plot_func_re=mod_psi_r_s_2D_re_i

        # imag part
        # reshape list --> 2D array
        mod_psi_r_s_2D_im=reshape_1D_to_2D(psi_x,psi_y,mod_psi_r_s_im)
        # shift origin
        mod_psi_r_s_2D_im=shift_origin_array(mod_psi_r_s_2D_im,Nl)
        # interpolate
        mod_psi_r_s_2D_im_i=splineinterp_2D(psi_x,psi_y,mod_psi_r_s_2D_im,3,3,0)
        plot_func_im=mod_psi_r_s_2D_im_i

    elif('tot' in sys.argv[1]):
        # modulus wfn.
        # reshape list --> 2D array
        mod_psi_r_tot_2D=reshape_1D_to_2D(psi_x,psi_y,mod_psi_r_tot)
        # shift origin
        mod_psi_r_tot_2D=shift_origin_array(mod_psi_r_tot_2D,Nl)
        # interpolate
        mod_psi_r_tot_2D_i=splineinterp_2D(psi_x,psi_y,mod_psi_r_tot_2D,3,3,0)
        plot_func=mod_psi_r_tot_2D_i

        # real part
        # reshape list --> 2D array
        mod_psi_r_tot_2D_re=reshape_1D_to_2D(psi_x,psi_y,mod_psi_r_tot_re)
        # shift origin
        mod_psi_r_tot_2D_re=shift_origin_array(mod_psi_r_tot_2D_re,Nl)
        # interpolate
        mod_psi_r_tot_2D_re_i=splineinterp_2D(psi_x,psi_y,mod_psi_r_tot_2D_re,3,3,0)
        plot_func_re=mod_psi_r_tot_2D_re_i

        # imag part
        # reshape list --> 2D array
        mod_psi_r_tot_2D_im=reshape_1D_to_2D(psi_x,psi_y,mod_psi_r_tot_im)
        # shift origin
        mod_psi_r_tot_2D_im=shift_origin_array(mod_psi_r_tot_2D_im,Nl)
        # interpolate
        mod_psi_r_tot_2D_im_i=splineinterp_2D(psi_x,psi_y,mod_psi_r_tot_2D_im,3,3,0)
        plot_func_im=mod_psi_r_tot_2D_im_i

    else:
        print 'must choose part of wfn to plot'
        wfnflag = str(raw_input("dd, sing, tot? : "))
        if('dd' in wfnflag):
            # mod wfn
            mod_psi_r_dd_2D=reshape_1D_to_2D(psi_x,psi_y,mod_psi_r_dd)
            mod_psi_r_dd_2D=shift_origin_array(mod_psi_r_dd_2D,Nl)
            mod_psi_r_dd_2D_i=splineinterp_2D(psi_x,psi_y,mod_psi_r_dd_2D,3,3,0)

            # real part
            mod_psi_r_dd_2D_re=reshape_1D_to_2D(psi_x,psi_y,mod_psi_r_dd_re)
            mod_psi_r_dd_2D_re=shift_origin_array(mod_psi_r_dd_2D_re,Nl)
            mod_psi_r_dd_2D_re_i=splineinterp_2D(psi_x,psi_y,mod_psi_r_dd_2D_re,3,3,0)
            # imag part
            mod_psi_r_dd_2D_im=reshape_1D_to_2D(psi_x,psi_y,mod_psi_r_dd_im)
            mod_psi_r_dd_2D_im=shift_origin_array(mod_psi_r_dd_2D_im,Nl)
            mod_psi_r_dd_2D_im_i=splineinterp_2D(psi_x,psi_y,mod_psi_r_dd_2D_im,3,3,0)

            plot_func=mod_psi_r_dd_2D_i
            plot_func_re=mod_psi_r_dd_2D_re_i
            plot_func_im=mod_psi_r_dd_2D_im_i

        elif('sing' in wfnflag):
            # mod wfn.
            mod_psi_r_s_2D=reshape_1D_to_2D(psi_x,psi_y,mod_psi_r_s)
            mod_psi_r_s_2D=shift_origin_array(mod_psi_r_s_2D,Nl)
            mod_psi_r_s_2D_i=splineinterp_2D(psi_x,psi_y,mod_psi_r_s_2D,3,3,0)

            # real part
            mod_psi_r_s_2D_re=reshape_1D_to_2D(psi_x,psi_y,mod_psi_r_s_re)
            mod_psi_r_s_2D_re=shift_origin_array(mod_psi_r_s_2D_re,Nl)
            mod_psi_r_s_2D_re_i=splineinterp_2D(psi_x,psi_y,mod_psi_r_s_2D_re,3,3,0)

            # imag part
            mod_psi_r_s_2D_im=reshape_1D_to_2D(psi_x,psi_y,mod_psi_r_s_im)
            mod_psi_r_s_2D_im=shift_origin_array(mod_psi_r_s_2D_im,Nl)
            mod_psi_r_s_2D_im_i=splineinterp_2D(psi_x,psi_y,mod_psi_r_s_2D_im,3,3,0)

            plot_func=mod_psi_r_s_2D_i
            plot_func_re=mod_psi_r_s_2D_re_i
            plot_func_im=mod_psi_r_s_2D_im_i

        elif('tot' in wfnflag):
            # mod wfn.
            mod_psi_r_tot_2D=reshape_1D_to_2D(psi_x,psi_y,mod_psi_r_tot)
            mod_psi_r_tot_2D=shift_origin_array(mod_psi_r_tot_2D,Nl)
            mod_psi_r_tot_2D_i=splineinterp_2D(psi_x,psi_y,mod_psi_r_tot_2D,3,3,0)

            # real part
            mod_psi_r_tot_2D_re=reshape_1D_to_2D(psi_x,psi_y,mod_psi_r_tot_re)
            mod_psi_r_tot_2D_re=shift_origin_array(mod_psi_r_tot_2D_re,Nl)
            mod_psi_r_tot_2D_re_i=splineinterp_2D(psi_x,psi_y,mod_psi_r_tot_2D_re,3,3,0)

            # imag part
            mod_psi_r_tot_2D_im=reshape_1D_to_2D(psi_x,psi_y,mod_psi_r_tot_im)
            mod_psi_r_tot_2D_im=shift_origin_array(mod_psi_r_tot_2D_im,Nl)
            mod_psi_r_tot_2D_im_i=splineinterp_2D(psi_x,psi_y,mod_psi_r_tot_2D_im,3,3,0)

            plot_func=mod_psi_r_tot_2D_i
            plot_func_re=mod_psi_r_tot_2D_re_i
            plot_func_im=mod_psi_r_tot_2D_im_i

#print mod_psi_r_dd.max()
# grab 1D slice
# slice defined by starting point (x0,y0)
# and direction vector (xsdir,ysdir)
x0=psi_x.min()
y0=1
xsdir=1
ysdir=0

#xsexact,ysexact,qsexact,qsexacterr=grab_slice(x0,y0,xsdir,ysdir,psi_x,psi_y,nui_nuj,[0]*len(ninj))

xsexact,ysexact,qsexact,qsexacterr=grab_slice(x0,y0,xsdir,ysdir,psi_x,psi_y,mod_dd,[0]*len(dd))
#print "xs :", xsexact
#print "ys :", ysexact
xsexact,ysexact,qsexact_re,qsexacterr=grab_slice(x0,y0,xsdir,ysdir,psi_x,psi_y,dd_re,[0]*len(dd_re))
#print "xs :", xsexact
#print "ys :", ysexact
xsexact,ysexact,qsexact_im,qsexacterr=grab_slice(x0,y0,xsdir,ysdir,psi_x,psi_y,dd_im,[0]*len(dd_im))
#print "sqrt(qs_re^2+qs_im^2) :", np.sqrt(np.array(qsexact_im)**2+np.array(qsexact_re)**2)
xsexact=np.array(xsexact)
ysexact=np.array(ysexact)

#print qsexact
# set up figure environment
fig = plt.figure()

#for 1D plot uncomment below
###############################

ax=fig.add_subplot(111)
ax.set_xticks(np.arange(xi.min(),xi.max()+1,1))

xnew=np.linspace(xsexact.min(),xsexact.max(),100)
#ynew=np.linspace(ysexact.min(),ysexact.max(),100)
#xnewfunc=xnew/(1.0*Nl[0])
ynew=np.empty(len(xnew))
ynew.fill(y0)
#ynewfunc=ynew/(1.0*Nl[1])
#ax.plot(xnew*np.sqrt(2.0),plot_func_trunc.ev(xnewfunc,ynewfunc))
ax.plot(xnew,plot_func.ev(ynew,xnew),color='orange')
ax.plot(xnew,plot_func_re.ev(ynew,xnew),color='blue')
ax.plot(xnew,plot_func_im.ev(ynew,xnew),color='green')
ax.plot(xnew,np.sqrt(plot_func_re.ev(ynew,xnew)**2+plot_func_im.ev(ynew,xnew)**2),color='purple')

#xsexact=np.sqrt(2.0)*xsexact
ax.plot(xsexact,qsexact,linestyle='None',color='red',marker='o')
ax.plot(xsexact,qsexact_re,linestyle='None',color='black',marker='s')
ax.plot(xsexact,qsexact_im,linestyle='None',color='cyan',marker='^')

#for 2D plot uncomment below
###############################
"""
ax=fig.add_subplot(111)
ax.set_xticks(np.arange(xi.min(),xi.max()+1,1))
ax.set_yticks(np.arange(yi.min(),yi.max()+1,1))
ax.set_xlabel('x',fontsize=24)
ax.set_ylabel('y',fontsize=24,rotation=360)
ax.xaxis.set_tick_params(labelsize=22)
ax.yaxis.set_tick_params(labelsize=22)

implot = ax.imshow(plot_func(xi,yi),extent=(xi.min(),xi.max(),yi.min(),yi.max()),cmap=cm.gnuplot2,origin='lower')

cb=fig.colorbar(implot,shrink=1.0,pad=0.01)
cb.ax.tick_params(labelsize=22)
"""
#for 3D plot uncomment below
###############################
"""
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

ax.scatter(psi_x, psi_y, dd_re+0.001, s=20, c='red')
surf = ax.plot_surface(xi, yi, plot_func_re(xi,yi), rstride=1, cstride=1, cmap=cm.gnuplot,
        linewidth=0, antialiased=False)


#,levels=np.arange(-0.0001,0.48,0.00005)
fig.colorbar(surf,shrink=0.75,pad=0.01)
"""
plt.subplots_adjust(bottom=0.15)

show()
#fig.savefig('check_new_script.pdf')
#fig.savefig('psi_r_triplet_Lx6_Ly24_N144_l0.5_U-2.5.pdf',bbox_inches='tight')


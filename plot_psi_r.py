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
import pyfftw
import glob
import os

rc('text', usetex=True)
rc('font', family='serif')

Nl=np.array([8,8])
Np=64
Ntot=Nl[0]*Nl[1]
Nmax=5000
Uhubb=-4.0
#lbda=0.5
dim_flag='2D'
params_flag=''

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
elif ('read' in sys.argv[1]):
    params_flag='params'
    params_file=str(sys.argv[2])
    if('MF' in params_file):
        params_flag='params_MF'
else: 
    print 'select an observable to plot'
    plot_flag=str(raw_input("didj, ninj, nk, psi_r? : "))
    arg_shift=0

fft_flag=''
for arg in sys.argv:
    if('open' in arg):
        open_flg='open_BCs'
    if('fft' in arg):
        fft_flag='fft_flag'
    if('1D' in arg):
        dim_flag='1D'
    if('2D' in arg):
        dim_flag='2D'
    if('3D' in arg):
        dim_flag='3D'

sym_flag=''
if('params' not in params_flag):
    lbda=float(raw_input("Enter lambda : "))
if(('1D' in dim_flag) and ('params' not in params_flag)):
    num_files=int(raw_input("number of files? : "))
if(('1D' not in dim_flag) and ('params' not in params_flag)):
    num_files=1

if('params' in params_flag):
    paramfile=open(params_file, 'r')
    num_files=int(paramfile.readline().split(":")[1])
    s_flag=np.empty(num_files,dtype=str)
    lam=np.empty(num_files)
    Uhubbard=np.empty(num_files)
    o_flg=np.empty(num_files,dtype=str)
    sym_flag_vec=np.empty(num_files,dtype=np.dtype((str, 3)))
    plot_flag=str(paramfile.readline().split(":")[1])
    script_dir=os.getcwd()
    if ('ntot' in plot_flag):
        plot_flag='ninj'
        corr_flag='ntot*ntot'
    elif ('nu*nu' in plot_flag):
        plot_flag='ninj'
        corr_flag='nu*nu'
    elif ('nu*nd' in plot_flag):
        plot_flag='ninj'
        corr_flag='nu*nd'
    else:
        plot_flag=plot_flag.strip()

if('ninj' in plot_flag):

    nui_ndj=np.empty((num_files,Nmax),dtype=np.float64)
    nui_ndj_err=np.empty((num_files,len(nui_ndj[0,:])),dtype=np.float64)
    nui_nuj=np.empty((num_files,Nmax),dtype=np.float64)
    nui_nuj_err=np.empty((num_files,len(nui_ndj[0,:])),dtype=np.float64)
    ninj=np.empty((num_files,len(nui_ndj[0,:])),dtype=np.float64)
    ninjerr=np.empty((num_files,len(nui_ndj[0,:])),dtype=np.float64)
    ninj_sym=np.empty((num_files,len(nui_ndj[0,:])),dtype=np.float64)
    ninj_sym_restored_1D=np.empty((num_files,len(nui_ndj[0,:])),dtype=np.float64)
    ninj_x=np.empty((num_files,Nmax))
    ninj_y=np.empty((num_files,Nmax))
    Ntot_vec=np.empty(num_files,dtype=np.int64)

    plt_fn_dict = {}
    xi_dict   = {}
    yi_dict   = {}

    for ind in range(num_files):
        if('params' in params_flag):
            lbda=float(paramfile.readline().split(":")[1])
            Nl[0]=int(paramfile.readline().split(":")[1])
            Nl[1]=int(paramfile.readline().split(":")[1])
            Uhubb=float(paramfile.readline().split(":")[1])
            lam[ind]=lbda
            Uhubbard[ind]=Uhubb
            Ntot_vec[ind]=Nl[0]*Nl[1]
            Ntot=Ntot_vec[ind]
            Np=Ntot
            o_flg[ind]=str(paramfile.readline().split(":")[1].strip())
            sym_flag=str(paramfile.readline().split(":")[1].strip())
            if (('Y' in sym_flag) or ('y' in sym_flag)):
                sym_flag="sym"            
            else:
                sym_flag=""
            fft_flag=str(paramfile.readline().split(":")[1].strip())
            sym_flag_vec[ind]=sym_flag
            dim_flag=str(paramfile.readline().split(":")[1].strip())
            datadir=str(paramfile.readline().split(":")[1].strip())
            if('MF' in params_flag):
                datadir=re.sub("LAMBDA",str(lam[ind]),datadir)
                datadir=re.sub("NP",str(Np),datadir)
                datadir=re.sub("UHUBB",str(Uhubb),datadir)
                print datadir
                os.chdir(datadir)
                for file in glob.glob("L"+str(Nl[0])+"_N"+str(Np)+"*prop_r.out"):
                    print file
                    data=np.loadtxt(file,dtype=np.float64)            
                os.chdir(script_dir)
            else:
                datadir=re.sub("LX",str(Nl[0]),datadir)
                datadir=re.sub("LY",str(Nl[1]),datadir)
                datadir=re.sub("UHUBB",str(Uhubb),datadir)
                os.chdir(datadir)
                for dir in glob.glob("m*--1.d0-"+str(lbda)+"-"+str(Uhubb)+"-"+str(Np)+"*"):
                    os.chdir(dir)
                    print os.getcwd()
                    for file in glob.glob("*"+str(plot_flag)+".dat"):
                        data=np.loadtxt(file,dtype=np.float64)
                    lattlabel=np.loadtxt("lattice-label.dat")
            s_flag[ind]=str(paramfile.readline().split(":")[1].strip())
            os.chdir(script_dir)
        else:
            Nl[0]=int(raw_input("Enter Lx : "))
            Nl[1]=int(raw_input("Enter Ly : "))
            Ntot_vec[ind]=Nl[0]*Nl[1]
            Ntot=Ntot_vec[ind]
            Np=Ntot
            sym_flag=str(raw_input("sym, Y/N? : "))
            if ('Y' in sym_flag or 'y' in sym_flag):
                sym_flag='sym'
            else:
                sym_flag=''

            data=np.loadtxt(sys.argv[1+arg_shift+2*ind],dtype=np.float64)
            lattlabel=np.loadtxt(sys.argv[2+arg_shift+2*ind])

        if('MF' not in params_flag):
            # nu*nd                       
            nui_ndj[ind,:Ntot]=np.array(data[Ntot:,1])
            nui_ndj_err[ind,:Ntot]=np.array(data[Ntot:,3])
            # nu*nu
            nui_nuj[ind,:Ntot]=np.array(data[:Ntot,1])
            nui_nuj_err[ind,:Ntot]=np.array(data[:Ntot,3])
            # nu*nd + nu*nu
            ninj[ind,:Ntot]=nui_ndj[ind,:Ntot]+nui_nuj[ind,:Ntot]
            ninjerr[ind,:Ntot]=nui_ndj_err[ind,:Ntot]+nui_nuj_err[ind,:Ntot]

            # scale down self-interaction peak at origin            
            nui_ndj[ind,0]=nui_ndj[ind,0]/2.0
            nui_nuj[ind,0]=nui_nuj[ind,0]/2.0
            ninj[ind,0]=ninj[ind,0]/2.0

            # shift, so lattice site 1 has coords (x,y)=(0,0) (not (1,1))
            shift=np.zeros(len(lattlabel[:,0]))
            shift.fill(1)
            ninj_x[ind,:Ntot]=lattlabel[:,1]-shift[:]
            ninj_y[ind,:Ntot]=lattlabel[:,2]-shift[:]
    
            # shift x,y coords (currently stored as lists)
            ninj_x[ind,:Ntot]=shift_origin_list(ninj_x[ind,:Ntot],Nl[0])
            ninj_y[ind,:Ntot]=shift_origin_list(ninj_y[ind,:Ntot],Nl[1])

        elif('MF' in params_flag):
            # nu*nd                       
            nui_ndj[ind,:Ntot]=np.array(data[:,18])
            nui_ndj_err[ind,:Ntot]=np.array(data[:,19])
            # nu*nu
            nui_nuj[ind,:Ntot]=np.array(data[:,20])
            nui_nuj_err[ind,:Ntot]=np.array(data[:,21])
            # nu*nd + nu*nu
            ninj[ind,:Ntot]=nui_ndj[ind,:Ntot]+nui_nuj[ind,:Ntot]
            ninjerr[ind,:Ntot]=nui_ndj_err[ind,:Ntot]+nui_nuj_err[ind,:Ntot]

            # scale down self-interaction peak at origin            
            nui_ndj[ind,0]=nui_ndj[ind,0]/2.0
            nui_nuj[ind,0]=nui_nuj[ind,0]/2.0
            ninj[ind,0]=ninj[ind,0]/2.0

            ninj_x[ind,:Ntot]=np.array(data[:,0])
            ninj_y[ind,:Ntot]=np.array(data[:,1])
            
        # restore symmetry to open system
        #ninj_sym=restore_symmetry(ninj_x[ind,:Ntot],ninj_y[ind,:Ntot],Nl,ninj[ind,:Ntot])
        """
        count=0
        while count<Ntot:
            for i in range(Nl[0]/2):
                if(ninj_x[ind,count]>1e-8):
                    if(abs(ninj_x[ind,count])>1e-8 and abs(ninj_y[ind,count])>1e-8 and abs(ninj_x[ind,count]+Nl[0]/2)>1e-8 and abs(ninj_y[ind,count]+Nl[1]/2)>1e-8):
                        ninj_sym[ind,count]=0.25*(ninj[ind,count]+ninj[ind,count+(Nl[0]-2*i)]
                                                  +ninj[ind,count+int(np.sign(ninj_y[ind,count])*Nl[0]*(Nl[1]-2*abs(ninj_y[ind,count])))]
                                                  +ninj[ind,count+int(np.sign(ninj_y[ind,count])*Nl[0]*(Nl[1]-2*abs(ninj_y[ind,count])))
                                                        +(Nl[0]-2*i)])
                        print count, i, '(',ninj_x[ind,count],',',ninj_y[ind,count],')', '(',ninj_x[ind,count+(Nl[0]-2*i)],ninj_y[ind,count+(Nl[0]-2*i)],')'
                        print count, i, '(',ninj_x[ind,count+int(np.sign(ninj_y[ind,count])*Nl[0]*(Nl[1]-2*abs(ninj_y[ind,count])))],',',ninj_y[ind,count+int(np.sign(ninj_y[ind,count])*Nl[0]*(Nl[1]-2*abs(ninj_y[ind,count])))],')'
                        print count, i, '(',ninj_x[ind,count+int(np.sign(ninj_y[ind,count])*Nl[0]*(Nl[1]-2*abs(ninj_y[ind,count])))+(Nl[0]-2*i)],',', ninj_y[ind,count+int(np.sign(ninj_y[ind,count])*Nl[0]*(Nl[1]-2*abs(ninj_y[ind,count])))+(Nl[0]-2*i)],')'
                    elif(abs(ninj_x[ind,count])>1e-8 and abs(ninj_y[ind,count])<1e-8 and abs(ninj_x[ind,count]+Nl[0]/2)>1e-8 and abs(ninj_y[ind,count]+Nl[1]/2)>1e-8):
                        ninj_sym[ind,count]=0.5*(ninj[ind,count]+ninj[ind,count+(Nl[0]-2*i)])
                    elif(abs(ninj_x[ind,count])<1e-8 and abs(ninj_y[ind,count])>1e-8 and abs(ninj_x[ind,count]+Nl[0]/2)>1e-8 and abs(ninj_y[ind,count]+Nl[1]/2)>1e-8):
                        ninj_sym[ind,count]=0.5*(ninj[ind,count]+ninj[ind,count+int(np.sign(ninj_y[ind,count])*Nl[0]*(Nl[1]-2*abs(ninj_y[ind,count])))])
                    else:
                        ninj_sym[ind,count]=ninj[ind,count]
                else:
                    if(abs(ninj_x[ind,count])>1e-8 and abs(ninj_y[ind,count])>1e-8 and abs(ninj_x[ind,count]+Nl[0]/2)>1e-8 and abs(ninj_y[ind,count]+Nl[1]/2)>1e-8):
                        ninj_sym[ind,count]=0.25*(ninj[ind,count]+ninj[ind,count-2*i]
                                                  +ninj[ind,count+int(np.sign(ninj_y[ind,count])*Nl[0]*(Nl[1]-2*abs(ninj_y[ind,count])))]
                                                  +ninj[ind,count+int(np.sign(ninj_y[ind,count])*Nl[0]*(Nl[1]-2*abs(ninj_y[ind,count])))
                                                        -2*i])
                        print count, i, ninj_x[ind,count], ninj_x[ind,count-2*i], ninj_y[ind,count], ninj_y[ind,count+int(np.sign(ninj_y[ind,count])*Nl[0]*(Nl[1]-2*abs(ninj_y[ind,count])))]
                    elif(abs(ninj_x[ind,count])>1e-8 and abs(ninj_y[ind,count])<1e-8 and abs(ninj_x[ind,count]+Nl[0]/2)>1e-8 and abs(ninj_y[ind,count]+Nl[1]/2)>1e-8):
                        ninj_sym[ind,count]=0.5*(ninj[ind,count]+ninj[ind,count-2*i])
                    elif(abs(ninj_x[ind,count])<1e-8 and abs(ninj_y[ind,count])>1e-8 and abs(ninj_x[ind,count]+Nl[0]/2)>1e-8 and abs(ninj_y[ind,count]+Nl[1]/2)>1e-8):
                        ninj_sym[ind,count]=0.5*(ninj[ind,count]+ninj[ind,count+int(np.sign(ninj_y[ind,count])*Nl[0]*(Nl[1]-2*abs(ninj_y[ind,count])))])
                    else:
                        ninj_sym[ind,count]=ninj[ind,count]
                count=count+1
        """

        xi, yi = np.ogrid[ninj_x[ind,:Ntot].min():ninj_x[ind,:Ntot].max():200j, ninj_y[ind,:Ntot].min():ninj_y[ind,:Ntot].max():200j]
        # ogrid gives transposed axes, so transpose back
        xi, yi = xi.T, yi.T
        xi_dict[ind]=xi
        yi_dict[ind]=yi

        if('params' not in params_flag):
            corr_flag = str(raw_input("nu*nu, nu*nd, ntot*ntot? : "))
        if('nu*nu' in corr_flag):
            # nu*nu
            # reshape list --> 2D array
            nui_nuj_2D=reshape_1D_to_2D(ninj_x[ind,:Ntot],ninj_y[ind,:Ntot],nui_nuj[ind,:Ntot])
            if('MF' not in params_flag):
                # shift origin
                nui_nuj_2D=shift_origin_array(nui_nuj_2D,Nl)
            if('MF' in params_flag):
                # shift origin
                nui_nuj_2D=shift_origin_array(nui_nuj_2D,Nl-[2,2])
            if('sym' in sym_flag): 
                # average to restore symmetry n_ij = 1/4 * (n_ij+n_{-i,j}+n_{i,-j}+n_{-i,-j})
                nui_nuj_2D[1:,1:]=0.25*(nui_nuj_2D[1:,1:]+np.flipud(nui_nuj_2D[1:,1:])+np.fliplr(nui_nuj_2D[1:,1:])+np.flipud(np.fliplr(nui_nuj_2D[1:,1:])))
                nui_nuj_sym_restored_1D[ind,:Ntot]=reshape_2D_to_1D(shift_origin_array(nui_nuj_2D,-1*Nl+1))
            # interpolate
            nui_nuj_2D_i=splineinterp_2D(ninj_x[ind,:Ntot],ninj_y[ind,:Ntot],nui_nuj_2D,3,3,0)
            plot_func=nui_nuj_2D_i
            plt_fn_dict[ind]=plot_func
        elif('nu*nd' in corr_flag):
            # nu*nd
            nui_ndj_2D=reshape_1D_to_2D(ninj_x[ind,:Ntot],ninj_y[ind,:Ntot],nui_ndj[ind,:Ntot])
            # shift origin
            nui_ndj_2D=shift_origin_array(nui_ndj_2D,Nl)
            if('sym' in sym_flag): 
                # average to restore symmetry n_ij = 1/4 * (n_ij+n_{-i,j}+n_{i,-j}+n_{-i,-j})
                nui_ndj_2D[1:,1:]=0.25*(nui_ndj_2D[1:,1:]+np.flipud(nui_ndj_2D[1:,1:])+np.fliplr(nui_ndj_2D[1:,1:])+np.flipud(np.fliplr(nui_ndj_2D[1:,1:])))
                nui_ndj_sym_restored_1D[ind,:Ntot]=reshape_2D_to_1D(shift_origin_array(nui_ndj_2D,-1*Nl+1))
            # interpolate
            nui_ndj_2D_i=splineinterp_2D(ninj_x[ind,:Ntot],ninj_y[ind,:Ntot],nui_ndj_2D,3,3,0)
            plot_func=nui_ndj_2D_i
            plt_fn_dict[ind]=plot_func
        elif('ntot*ntot' in corr_flag):
            # ntot*ntot
            ninj_2D=reshape_1D_to_2D(ninj_x[ind,:Ntot],ninj_y[ind,:Ntot],ninj[ind,:Ntot])
            if('MF' not in params_flag):
                # shift origin
                ninj_2D=shift_origin_array(ninj_2D,Nl)
            if('MF' in params_flag):
                # shift origin
                ninj_2D=shift_origin_array(ninj_2D,Nl-[2,2])
            if('sym' in sym_flag):
                # average to restore symmetry n_ij = 1/4 * (n_ij+n_{-i,j}+n_{i,-j}+n_{-i,-j})
                ninj_2D[1:,1:]=0.25*(ninj_2D[1:,1:]+np.flipud(ninj_2D[1:,1:])+np.fliplr(ninj_2D[1:,1:])+np.flipud(np.fliplr(ninj_2D[1:,1:])))
                ninj_sym_restored_1D[ind,:Ntot]=reshape_2D_to_1D(shift_origin_array(ninj_2D,-1*Nl+1))
            # interpolate
            ninj_2D_i=splineinterp_2D(ninj_x[ind,:Ntot],ninj_y[ind,:Ntot],ninj_2D,3,3,0)
            plot_func=ninj_2D_i
            plt_fn_dict[ind]=plot_func

elif('didj' in plot_flag):

    didj=np.empty((num_files,Nmax),dtype=np.float64)
    didjerr=np.empty((num_files,len(didj[0,:])),dtype=np.float64)
    didj_x=np.empty((num_files,Nmax))
    didj_y=np.empty((num_files,Nmax))
    didj_sym=np.empty((num_files,len(didj[0,:])),dtype=np.float64)
    didj_sym_restored_1D=np.empty((num_files,len(didj[0,:])),dtype=np.float64)
    Ntot_vec=np.empty(num_files,dtype=np.int64)    

    plt_fn_dict = {}
    xi_dict   = {}
    yi_dict   = {}

    for ind in range(num_files):
        if('params' in params_flag):
            lbda=float(paramfile.readline().split(":")[1])
            Nl[0]=int(paramfile.readline().split(":")[1])
            Nl[1]=int(paramfile.readline().split(":")[1])
            Uhubb=float(paramfile.readline().split(":")[1])
            lam[ind]=lbda
            Uhubbard[ind]=Uhubb
            Ntot_vec[ind]=Nl[0]*Nl[1]
            Ntot=Ntot_vec[ind]
            Np=Ntot
            o_flg[ind]=str(paramfile.readline().split(":")[1].strip())
            sym_flag=str(paramfile.readline().split(":")[1].strip())            
            if (('Y' in sym_flag) or ('y' in sym_flag)):
                sym_flag="sym"            
            else:
                sym_flag=""
            fft_flag=str(paramfile.readline().split(":")[1].strip())
            sym_flag_vec[ind]=sym_flag
            dim_flag=str(paramfile.readline().split(":")[1].strip())
            datadir=str(paramfile.readline().split(":")[1].strip())
            if('MF' in params_flag):
                datadir=re.sub("LAMBDA",str(lam[ind]),datadir)
                datadir=re.sub("NP",str(Np),datadir)
                datadir=re.sub("UHUBB",str(Uhubb),datadir)
                print datadir
                os.chdir(datadir)
                for file in glob.glob("L"+str(Nl[0])+"_N"+str(Np)+"*prop_r.out"):
                    print file
                    data=np.loadtxt(file,dtype=np.float64)            
                os.chdir(script_dir)
            else:
                datadir=re.sub("LX",str(Nl[0]),datadir)
                datadir=re.sub("LY",str(Nl[1]),datadir)
                datadir=re.sub("UHUBB",str(Uhubb),datadir)
                os.chdir(datadir)
                for dir in glob.glob("m*--1.d0-"+str(lbda)+"-"+str(Uhubb)+"-"+str(Np)+"*"):
                    os.chdir(dir)
                    print os.getcwd()
                    for file in glob.glob("*"+str(plot_flag)+".dat"):
                        data=np.loadtxt(file,dtype=np.float64)
                    lattlabel=np.loadtxt("lattice-label.dat")
            s_flag[ind]=str(paramfile.readline().split(":")[1].strip())
            os.chdir(script_dir)
        else:
            Nl[0]=int(raw_input("Enter Lx : "))
            Nl[1]=int(raw_input("Enter Ly : "))
            Ntot_vec[ind]=Nl[0]*Nl[1]
            Ntot=Ntot_vec[ind]
            Np=Ntot
            sym_flag=str(raw_input("sym, Y/N? : "))
            if ('Y' in sym_flag or 'y' in sym_flag):
                sym_flag='sym'
            else:
                sym_flag=''
            #load didj data
            data=np.loadtxt(sys.argv[1+arg_shift+2*ind],dtype=np.float64)
            lattlabel=np.loadtxt(sys.argv[2+arg_shift+2*ind])

        if('MF' not in params_flag):
            didj[ind,:Ntot]=data[:,1]
            didjerr[ind,:Ntot]=data[:,3]
            # scale origin
            didj[ind,0]=didj[ind,0]/10.0
            # normalize
            #didj=didj[:]*(Ntot)**2/(Np*(Np-1))
            #didjerr=didjerr[:]*(Ntot)**2/(Np*(Np-1))
            # load x, y coords from lattice-label.dat
        
            # shift, so lattice site 1 has coords (x,y)=(0,0) (not (1,1))
            shift=np.zeros(len(lattlabel[:,0]))
            shift.fill(1)
            didj_x[ind,:Ntot]=lattlabel[:,1]-shift[:]
            didj_y[ind,:Ntot]=lattlabel[:,2]-shift[:]
    
            # shift x,y coords (currently stored as lists)
            didj_x[ind,:Ntot]=shift_origin_list(didj_x[ind,:Ntot],Nl[0])
            didj_y[ind,:Ntot]=shift_origin_list(didj_y[ind,:Ntot],Nl[1])

        elif('MF' in params_flag):
            didj[ind,:Ntot]=data[:,16]
            didjerr[ind,:Ntot]=data[:,17]
            # scale origin
            didj[ind,0]=didj[ind,0]/10.0

            didj_x[ind,:Ntot]=data[:,0]
            didj_y[ind,:Ntot]=data[:,1]

        """
        if('sym' in sym_flag): 
            # average to restore symmetry n_ij = 1/4 * (n_ij+n_{-i,j}+n_{i,-j}+n_{-i,-j})
            didj_sym=restore_symmetry(didj_x[ind,:Ntot],didj_y[ind,:Ntot],Nl,didj[ind,:Ntot])
            didj_2D=reshape_1D_to_2D(didj_x[ind,:Ntot],didj_y[ind,:Ntot],didj_sym)
            didj_2D=shift_origin_array(didj_2D,Nl)
        """

        xi, yi = np.ogrid[didj_x[ind,:Ntot].min():didj_x[ind,:Ntot].max():200j, didj_y[ind,:Ntot].min():didj_y[ind,:Ntot].max():200j]
        # ogrid gives transposed axes, so transpose back
        xi, yi = xi.T, yi.T
        xi_dict[ind]=xi
        yi_dict[ind]=yi

        # restore symmetry
        #didj_sym=restore_symmetry(didj_x[ind,:Ntot],didj_y[ind,:Ntot],Nl,didj[ind,:Ntot])
        # reshape list --> 2D array
        didj_2D=reshape_1D_to_2D(didj_x[ind,:Ntot],didj_y[ind,:Ntot],didj[ind,:Ntot])
        if('MF' not in params_flag):
            # shift origin
            didj_2D=shift_origin_array(didj_2D,Nl)
        elif('MF' in params_flag):
            # shift origin
            didj_2D=shift_origin_array(didj_2D,Nl-[2,2])

        if('sym' in sym_flag): 
            # average to restore symmetry n_ij = 1/4 * (n_ij+n_{-i,j}+n_{i,-j}+n_{-i,-j})
            #didj_sym[ind,:Ntot]=restore_symmetry(didj_x[ind,:Ntot],didj_y[ind,:Ntot],Nl,didj[ind,:Ntot])
            #didj_2D=reshape_1D_to_2D(didj_x[ind,:Ntot],didj_y[ind,:Ntot],didj_sym[ind,:Ntot])
            #didj_2D=shift_origin_array(didj_2D,Nl)
            didj_2D[1:,1:]=0.25*(didj_2D[1:,1:]+np.flipud(didj_2D[1:,1:])+np.fliplr(didj_2D[1:,1:])+np.flipud(np.fliplr(didj_2D[1:,1:])))
            # to plot 1D data, reshift origin, and reshape to 1D list 
            didj_sym_restored_1D[ind,:Ntot]=reshape_2D_to_1D(shift_origin_array(didj_2D,-1*Nl+1))
            """
            for i in range(Ntot):
                print didj_sym_restored_1D[i]-didj_sym[ind,i]
            """
        # interpolate
        didj_2D_i=splineinterp_2D(didj_x[ind,:Ntot],didj_y[ind,:Ntot],didj_2D,3,3,0)
        plot_func=didj_2D_i
        plt_fn_dict[ind]=plot_func


elif('psi_r' in plot_flag):

    psidata=np.empty((num_files,Nmax,2),dtype=np.float64)
    psi_r=np.empty((num_files,Nmax),dtype=np.complex128)
    mod_psi_r=np.empty((num_files,Nmax),dtype=np.float64)
    dd=np.empty((num_files,Nmax),dtype=np.complex128)
    mod_dd=np.empty((num_files,Nmax),dtype=np.float64)
    dd_re=np.empty((num_files,Nmax),dtype=np.float64)
    dd_im=np.empty((num_files,Nmax),dtype=np.float64)
    sing=np.empty((num_files,Nmax),dtype=np.complex128)
    mod_psi_r_s=np.empty((num_files,Nmax),dtype=np.float64)
    sing_re=np.empty((num_files,Nmax),dtype=np.float64)
    sing_im=np.empty((num_files,Nmax),dtype=np.float64)
    mod_psi_r_tot=np.empty((num_files,Nmax),dtype=np.float64)
    tot_re=np.empty((num_files,Nmax),dtype=np.float64)
    tot_im=np.empty((num_files,Nmax),dtype=np.float64)
    psi_x=np.empty((num_files,Nmax))
    psi_y=np.empty((num_files,Nmax))
    #psi_sym=np.empty((num_files,len(psidata[0,:])),dtype=np.float64)
    #psi_sym_restored_1D=np.empty((num_files,len(psidata[0,:])),dtype=np.float64)
    Ntot_vec=np.empty(num_files,dtype=np.int64)    

    plt_fn_dict = {}
    plt_fn_re_dict = {}
    plt_fn_im_dict = {}
    xi_dict   = {}
    yi_dict   = {}

    for ind in range(num_files):
        if('params' in params_flag):
            lbda=float(paramfile.readline().split(":")[1])
            Nl[0]=int(paramfile.readline().split(":")[1])
            Nl[1]=int(paramfile.readline().split(":")[1])
            Uhubb=float(paramfile.readline().split(":")[1])
            lam[ind]=lbda
            Uhubbard[ind]=Uhubb
            Ntot_vec[ind]=Nl[0]*Nl[1]
            Ntot=Ntot_vec[ind]
            Np=Ntot
            o_flg[ind]=str(paramfile.readline().split(":")[1].strip())
            sym_flag=str(paramfile.readline().split(":")[1].strip())            
            if (('Y' in sym_flag) or ('y' in sym_flag)):
                sym_flag="sym"            
            else:
                sym_flag=""
            sym_flag_vec[ind]=sym_flag
            fft_flag=str(paramfile.readline().split(":")[1].strip())
            dim_flag=str(paramfile.readline().split(":")[1].strip())
            datadir=str(paramfile.readline().split(":")[1].strip())
            if('MF' in params_flag):
                datadir=re.sub("LAMBDA",str(lam[ind]),datadir)
                datadir=re.sub("NP",str(Np),datadir)
                datadir=re.sub("UHUBB",str(Uhubb),datadir)
                print datadir
                os.chdir(datadir)
                for file in glob.glob("L"+str(Nl[0])+"_N"+str(Np)+"*prop_r.out"):
                    print file
                    psidata=np.loadtxt(file,dtype=np.float64)            
                os.chdir(script_dir)
            else:
                datadir=re.sub("LX",str(Nl[0]),datadir)
                datadir=re.sub("LY",str(Nl[1]),datadir)
                datadir=re.sub("UHUBB",str(Uhubb),datadir)
                os.chdir(datadir)
                if (('Y' in fft_flag) or ('y' in fft_flag)):
                    fft_flag="fft"            
                else:
                    fft_flag=""
                for dir in glob.glob("m*--1.d0-"+str(lbda)+"-"+str(Uhubb)+"-"+str(Np)+"*"):
                    os.chdir(dir)
                    print os.getcwd()
                    if('fft' in fft_flag):
                        psidata[ind,:3*Ntot,:]=np.loadtxt("pair_wf_normed.dat",dtype=np.float64)
                    else:
                        psidata[ind,:3*Ntot,:]=np.loadtxt("psi_r_normed.dat",dtype=np.float64)
                    lattlabel=np.loadtxt("lattice-label.dat")
            s_flag[ind]=str(paramfile.readline().split(":")[1].strip())
            os.chdir(script_dir)
        else:
            Nl[0]=int(raw_input("Enter Lx : "))
            Nl[1]=int(raw_input("Enter Ly : "))
            Ntot_vec[ind]=Nl[0]*Nl[1]
            Ntot=Ntot_vec[ind]
            Np=Ntot
            sym_flag=str(raw_input("sym, Y/N? : "))
            if ('Y' in sym_flag or 'y' in sym_flag):
                sym_flag='sym'
            else:
                sym_flag=''
            #load psi data
            psidata[ind,:3*Ntot,:]=np.loadtxt(sys.argv[1+arg_shift+2*ind],dtype=np.float64)
            lattlabel=np.loadtxt(sys.argv[2+arg_shift+2*ind])
    
        if('MF' not in params_flag):
            # shift, so lattice site 1 has coords (x,y)=(0,0) (not (1,1))
            shift=np.zeros(len(lattlabel[:,0]))
            shift.fill(1)
            # shift, so lattice site 1 has coords (x,y)=(0,0) (not (1,1))
            psi_x[ind,:Ntot]=lattlabel[:,1]-shift[:]
            psi_y[ind,:Ntot]=lattlabel[:,2]-shift[:]

            if('fft' in fft_flag):
                # load pair_wf_normed, do FFT to get psi_r_normed
                psi_k_uu=pyfftw.empty_aligned((Nl[0],Nl[1]),dtype=np.complex128)
                psi_r_uu=pyfftw.empty_aligned((Nl[0],Nl[1]),dtype=np.complex128)
                psi_k_dd=pyfftw.empty_aligned((Nl[0],Nl[1]),dtype=np.complex128)
                psi_r_dd=pyfftw.empty_aligned((Nl[0],Nl[1]),dtype=np.complex128)
                psi_k_s =pyfftw.empty_aligned((Nl[0],Nl[1]),dtype=np.complex128)
                psi_r_s =pyfftw.empty_aligned((Nl[0],Nl[1]),dtype=np.complex128)

                fft_object_psi_k_uu = pyfftw.FFTW(psi_k_uu,psi_r_uu,axes=(0,1),direction='FFTW_BACKWARD')
                fft_object_psi_k_dd = pyfftw.FFTW(psi_k_dd,psi_r_dd,axes=(0,1),direction='FFTW_BACKWARD')
                fft_object_psi_k_s  = pyfftw.FFTW(psi_k_s ,psi_r_s ,axes=(0,1),direction='FFTW_BACKWARD')

                psi_k_uu=np.reshape(psidata[ind,:Ntot,0]+1j*psidata[ind,:Ntot,1],(Nl[0],Nl[1]))
                psi_k_dd=np.reshape(psidata[ind,Ntot:2*Ntot,0]+1j*psidata[ind,Ntot:2*Ntot,1],(Nl[0],Nl[1]))
                psi_k_s =np.reshape(psidata[ind,2*Ntot:3*Ntot,0]+1j*psidata[ind,2*Ntot:3*Ntot,1],(Nl[0],Nl[1]))

                # FFT psi_k into psi_r
                fft_object_psi_k_uu(psi_k_uu)
                fft_object_psi_k_dd(psi_k_dd)
                fft_object_psi_k_s(psi_k_s)

                # normalize psi_r
                psi_r_uu=(np.sqrt(Nl[0])*np.sqrt(Nl[1]))*psi_r_uu[:,:]
                psi_r_dd=(np.sqrt(Nl[0])*np.sqrt(Nl[1]))*psi_r_dd[:,:]
                psi_r_s =(np.sqrt(Nl[0])*np.sqrt(Nl[1]))*psi_r_s[:,:]
            
                # reshape psi_r to 1D array
                psidata[ind,:Ntot,0]=np.reshape(psi_r_uu[:,:].real,Ntot)
                psidata[ind,:Ntot,1]=np.reshape(psi_r_uu[:,:].imag,Ntot)
                psidata[ind,Ntot:2*Ntot,0]=np.reshape(psi_r_dd[:,:].real,Ntot)
                psidata[ind,Ntot:2*Ntot,1]=np.reshape(psi_r_dd[:,:].imag,Ntot)
                psidata[ind,2*Ntot:3*Ntot,0]=np.reshape(psi_r_s[:,:].real,Ntot)
                psidata[ind,2*Ntot:3*Ntot,1]=np.reshape(psi_r_s[:,:].imag,Ntot)

            # grab psi_dd, psi_uu, psi_s
            # compute |psi_dd|, |psi_uu|, |psi_s| 
            ##############################################################
            # psi_uu
            psi_r[ind,:Ntot]=psidata[ind,:Ntot,0]+1j*psidata[ind,:Ntot,1]
            # |psi_uu|
            mod_psi_r[ind,:Ntot]=np.sqrt((psi_r[ind,:Ntot].real)**2+(psi_r[ind,:Ntot].imag)**2)
            ##############################################################
            # psi_dd
            dd[ind,:Ntot]=psidata[ind,Ntot:2*Ntot,0]+1j*psidata[ind,Ntot:2*Ntot,1]
                # |psi_dd|
            mod_dd[ind,:Ntot]=np.sqrt((dd[ind,:Ntot].real)**2+(dd[ind,:Ntot].imag)**2)
            dd_re[ind,:Ntot]=dd[ind,:Ntot].real
            dd_im[ind,:Ntot]=dd[ind,:Ntot].imag
            ##############################################################
            # psi_s
            sing[ind,:Ntot]=psidata[ind,2*Ntot:3*Ntot,0]+1j*psidata[ind,2*Ntot:3*Ntot,1]
            # |psi_s|
            mod_psi_r_s[ind,:Ntot]=np.sqrt((sing[ind,:Ntot].real)**2+(sing[ind,:Ntot].imag)**2)
            sing_re[ind,:Ntot]=sing[ind,:Ntot].real
            sing_im[ind,:Ntot]=sing[ind,:Ntot].imag
            ##############################################################
            # |psi_tot|
            mod_psi_r_tot[ind,:Ntot]=mod_psi_r_s[ind,:Ntot]+mod_dd[ind,:Ntot]+mod_psi_r[ind,:Ntot]
            # ERROR HERE (|psi_tot| =\= sqrt(Re(psi_s)+Re(psi_uu)+Re(psi_dd))^2+Im(psi_s)+Im(psi_uu)+Im(psi_dd))^2)
            tot_re[ind,:Ntot]=sing[ind,:Ntot].real+dd[ind,:Ntot].real+psi_r[ind,:Ntot].real
            tot_im[ind,:Ntot]=sing[ind,:Ntot].imag+dd[ind,:Ntot].imag+psi_r[ind,:Ntot].imag

            # shift x,y coords (currently stored as lists), so origin is at center, not edge of plot 
            psi_x[ind,:Ntot]=shift_origin_list(psi_x[ind,:Ntot],Nl[0])
            psi_y[ind,:Ntot]=shift_origin_list(psi_y[ind,:Ntot],Nl[1])
        
        elif('MF' in params_flag):
            # grab psi_dd, psi_uu, psi_s
            # compute |psi_dd|, |psi_uu|, |psi_s| 
            ##############################################################
            # psi_uu
            psi_r[ind,:Ntot]=psidata[:,10]+1j*psidata[:,11]
            # |psi_uu|
            mod_psi_r[ind,:Ntot]=np.sqrt((psi_r[ind,:Ntot].real)**2+(psi_r[ind,:Ntot].imag)**2)
            ##############################################################
            # psi_dd
            dd[ind,:Ntot]=psidata[:,12]+1j*psidata[:,13]
                # |psi_dd|
            mod_dd[ind,:Ntot]=np.sqrt((dd[ind,:Ntot].real)**2+(dd[ind,:Ntot].imag)**2)
            dd_re[ind,:Ntot]=dd[ind,:Ntot].real
            dd_im[ind,:Ntot]=dd[ind,:Ntot].imag
            ##############################################################
            # psi_s
            sing[ind,:Ntot]=psidata[:,14]+1j*psidata[:,15]
            # |psi_s|
            mod_psi_r_s[ind,:Ntot]=np.sqrt((sing[ind,:Ntot].real)**2+(sing[ind,:Ntot].imag)**2)
            sing_re[ind,:Ntot]=sing[ind,:Ntot].real
            sing_im[ind,:Ntot]=sing[ind,:Ntot].imag
            ##############################################################
            # |psi_tot|
            mod_psi_r_tot[ind,:Ntot]=mod_psi_r_s[ind,:Ntot]+mod_dd[ind,:Ntot]+mod_psi_r[ind,:Ntot]
            # ERROR HERE (|psi_tot| =\= sqrt(Re(psi_s)+Re(psi_uu)+Re(psi_dd))^2+Im(psi_s)+Im(psi_uu)+Im(psi_dd))^2)
            tot_re[ind,:Ntot]=sing[ind,:Ntot].real+dd[ind,:Ntot].real+psi_r[ind,:Ntot].real
            tot_im[ind,:Ntot]=sing[ind,:Ntot].imag+dd[ind,:Ntot].imag+psi_r[ind,:Ntot].imag

            psi_x[ind,:Ntot]=psidata[:,0]
            psi_y[ind,:Ntot]=psidata[:,1]

        xi, yi = np.ogrid[psi_x.min():psi_x.max():200j, psi_y.min():psi_y.max():200j]
        # ogrid gives transposed axes, so transpose back
        xi, yi = xi.T, yi.T
        xi_dict[ind]=xi
        yi_dict[ind]=yi

        if(('dd' in plot_flag) or ('dd' in sys.argv[1])):
            if (('dd' not in plot_flag) and ('r_s' not in plot_flag) and ('tot' not in plot_flag)):
                plot_flag=plot_flag+'_dd'
            # modulus wfn
            # reshape list --> 2D array
            dd_2D=reshape_1D_to_2D(psi_x[ind,:Ntot],psi_y[ind,:Ntot],mod_dd[ind,:Ntot])
            if('MF' not in params_flag):
                # shift origin
                dd_2D=shift_origin_array(dd_2D,Nl)
            if('MF' in params_flag):
                # shift origin
                dd_2D=shift_origin_array(dd_2D,Nl-[2,2])
            # interpolate
            dd_2D_i=splineinterp_2D(psi_x,psi_y,dd_2D,3,3,0)
            plot_func=dd_2D_i
            plt_fn_dict[ind]=plot_func

            # real part
            # reshape list --> 2D array
            dd_2D_re=reshape_1D_to_2D(psi_x[ind,:Ntot],psi_y[ind,:Ntot],dd_re[ind,:Ntot])
            if('MF' not in params_flag):
                # shift origin
                dd_2D_re=shift_origin_array(dd_2D_re,Nl)
            if('MF' in params_flag):
                # shift origin
                dd_2D_re=shift_origin_array(dd_2D_re,Nl-[2,2])
            # interpolate
            dd_2D_re_i=splineinterp_2D(psi_x,psi_y,dd_2D_re,3,3,0)
            plot_func_re=dd_2D_re_i
            plt_fn_re_dict[ind]=plot_func_re

            # imag part
            # reshape list --> 2D array
            dd_2D_im=reshape_1D_to_2D(psi_x[ind,:Ntot],psi_y[ind,:Ntot],dd_im[ind,:Ntot])
            if('MF' not in params_flag):
                # shift origin
                dd_2D_im=shift_origin_array(dd_2D_im,Nl)
            if('MF' in params_flag):
                # shift origin
                dd_2D_im=shift_origin_array(dd_2D_im,Nl-[2,2])
            # interpolate
            dd_2D_im_i=splineinterp_2D(psi_x[ind,:Ntot],psi_y[ind,:Ntot],dd_2D_im,3,3,0)
            plot_func_im=dd_2D_im_i
            plt_fn_im_dict[ind]=plot_func_im
        
        elif(('sing' in plot_flag) or ('_s' in plot_flag) or ('sing' in sys.argv[1]) or ('_s' in sys.argv[1])):
            if (('dd' not in plot_flag) and ('r_s' not in plot_flag) and ('tot' not in plot_flag)):
                plot_flag=plot_flag+'_s'
            # modulus wfn
            # reshape list --> 2D array
            mod_psi_r_s_2D=reshape_1D_to_2D(psi_x[ind,:Ntot],psi_y[ind,:Ntot],mod_psi_r_s[ind,:Ntot])
            if('MF' not in params_flag):
                # shift origin
                mod_psi_r_s_2D=shift_origin_array(mod_psi_r_s_2D,Nl)
            if('MF' in params_flag):
                # shift origin
                mod_psi_r_s_2D=shift_origin_array(mod_psi_r_s_2D,Nl-[2,2])
            # interpolate
            mod_psi_r_s_2D_i=splineinterp_2D(psi_x,psi_y,mod_psi_r_s_2D,3,3,0)
            plot_func=mod_psi_r_s_2D_i
            plt_fn_dict[ind]=plot_func

            # real part
            # reshape list --> 2D array
            sing_2D_re=reshape_1D_to_2D(psi_x[ind,:Ntot],psi_y[ind,:Ntot],sing_re[ind,:Ntot])
            if('MF' not in params_flag):
                # shift origin
                sing_2D_re=shift_origin_array(sing_2D_re,Nl)
            if('MF' in params_flag):
                # shift origin
                sing_2D_re=shift_origin_array(sing_2D_re,Nl-[2,2])
            # interpolate
            sing_2D_re_i=splineinterp_2D(psi_x[ind,:Ntot],psi_y[ind,:Ntot],sing_2D_re,3,3,0)
            plot_func_re=sing_2D_re_i
            plt_fn_re_dict[ind]=plot_func_re

            # imag part
            # reshape list --> 2D array
            sing_2D_im=reshape_1D_to_2D(psi_x[ind,:Ntot],psi_y[ind,:Ntot],sing_im[ind,:Ntot])
            if('MF' not in params_flag):
                # shift origin
                sing_2D_im=shift_origin_array(sing_2D_im,Nl)
            if('MF' in params_flag):
                # shift origin
                sing_2D_im=shift_origin_array(sing_2D_im,Nl-[2,2])
            # interpolate
            sing_2D_im_i=splineinterp_2D(psi_x[ind,:Ntot],psi_y[ind,:Ntot],sing_2D_im,3,3,0)
            plot_func_im=sing_2D_im_i
            plt_fn_im_dict[ind]=plot_func_im

        elif(('tot' in plot_flag) or ('tot' in sys.argv[1])):
            if (('dd' not in plot_flag) and ('r_s' not in plot_flag) and ('tot' not in plot_flag)):
                plot_flag=plot_flag+'_tot'
            # modulus wfn.
            # reshape list --> 2D array
            mod_psi_r_tot_2D=reshape_1D_to_2D(psi_x[ind,:Ntot],psi_y[ind,:Ntot],mod_psi_r_tot[ind,:Ntot])
            if('MF' not in params_flag):
                # shift origin
                mod_psi_r_tot_2D=shift_origin_array(mod_psi_r_tot_2D,Nl)
            if('MF' in params_flag):
                # shift origin
                mod_psi_r_tot_2D=shift_origin_array(mod_psi_r_tot_2D,Nl-[2,2])
            # interpolate
            mod_psi_r_tot_2D_i=splineinterp_2D(psi_x[ind,:Ntot],psi_y[ind,:Ntot],mod_psi_r_tot_2D,3,3,0)
            plot_func=mod_psi_r_tot_2D_i
            plt_fn_dict[ind]=plot_func
            
            # real part
            # reshape list --> 2D array
            tot_2D_re=reshape_1D_to_2D(psi_x[ind,:Ntot],psi_y[ind,:Ntot],tot_re[ind,:Ntot])
            if('MF' not in params_flag):
                # shift origin
                tot_2D_re=shift_origin_array(tot_2D_re,Nl)
            if('MF' in params_flag):
                # shift origin
                tot_2D_re=shift_origin_array(tot_2D_re,Nl-[2,2])
            # interpolate
            tot_2D_re_i=splineinterp_2D(psi_x[ind,:Ntot],psi_y[ind,:Ntot],tot_2D_re,3,3,0)
            plot_func_re=tot_2D_re_i
            plt_fn_re_dict[ind]=plot_func_re

            # imag part
            # reshape list --> 2D array
            tot_2D_im=reshape_1D_to_2D(psi_x[ind,:Ntot],psi_y[ind,:Ntot],tot_im[ind,:Ntot])
            if('MF' not in params_flag):
                # shift origin
                tot_2D_im=shift_origin_array(tot_2D_im,Nl)
            if('MF' in params_flag):
                # shift origin
                tot_2D_im=shift_origin_array(tot_2D_im,Nl-[2,2])
            # interpolate
            tot_2D_im_i=splineinterp_2D(psi_x[ind,:Ntot],psi_y[ind,:Ntot],tot_2D_im,3,3,0)
            plot_func_im=tot_2D_im_i
            plt_fn_im_dict[ind]=plot_func_im

    """
    else:
        print 'must choose part of wfn to plot'
        wfnflag = str(raw_input("dd, sing, tot? : "))
        if('dd' in wfnflag):
            plot_flag=plot_flag+'_dd'
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
            plot_flag=plot_flag+'_s'
            # mod wfn.
            mod_psi_r_s_2D=reshape_1D_to_2D(psi_x,psi_y,mod_psi_r_s)
            mod_psi_r_s_2D=shift_origin_array(mod_psi_r_s_2D,Nl)
            mod_psi_r_s_2D_i=splineinterp_2D(psi_x,psi_y,mod_psi_r_s_2D,3,3,0)

            # real part
            sing_2D_re=reshape_1D_to_2D(psi_x,psi_y,sing_re)
            sing_2D_re=shift_origin_array(sing_2D_re,Nl)
            sing_2D_re_i=splineinterp_2D(psi_x,psi_y,sing_2D_re,3,3,0)

            # imag part
            sing_2D_im=reshape_1D_to_2D(psi_x,psi_y,sing_im)
            sing_2D_im=shift_origin_array(sing_2D_im,Nl)
            sing_2D_im_i=splineinterp_2D(psi_x,psi_y,sing_2D_im,3,3,0)

            plot_func=mod_psi_r_s_2D_i
            plot_func_re=sing_2D_re_i
            plot_func_im=sing_2D_im_i

        elif('tot' in wfnflag):
            plot_flag=plot_flag+'_tot'
            # mod wfn.
            mod_psi_r_tot_2D=reshape_1D_to_2D(psi_x,psi_y,mod_psi_r_tot)
            mod_psi_r_tot_2D=shift_origin_array(mod_psi_r_tot_2D,Nl)
            mod_psi_r_tot_2D_i=splineinterp_2D(psi_x,psi_y,mod_psi_r_tot_2D,3,3,0)

            # real part
            tot_2D_re=reshape_1D_to_2D(psi_x,psi_y,tot_re)
            tot_2D_re=shift_origin_array(tot_2D_re,Nl)
            tot_2D_re_i=splineinterp_2D(psi_x,psi_y,tot_2D_re,3,3,0)

            # imag part
            tot_2D_im=reshape_1D_to_2D(psi_x,psi_y,tot_im)
            tot_2D_im=shift_origin_array(tot_2D_im,Nl)
            tot_2D_im_i=splineinterp_2D(psi_x,psi_y,tot_2D_im,3,3,0)

            plot_func=mod_psi_r_tot_2D_i
            plot_func_re=tot_2D_re_i
            plot_func_im=tot_2D_im_i

    """
#print mod_psi_r_dd.max()
# grab 1D slice
# slice defined by starting point (x0,y0)
# and direction vector (xsdir,ysdir)

#xsexact,ysexact,qsexact,qsexacterr=grab_slice(x0,y0,xsdir,ysdir,psi_x,psi_y,ninj,[0]*len(ninj))

#xsexact,ysexact,qsexact,qsexacterr=grab_slice(x0,y0,xsdir,ysdir,ninj_x,ninj_y,ninj,[0]*len(ninj))
#xsexact,ysexact,qsexact_re,qsexacterr=grab_slice(x0,y0,xsdir,ysdir,ninj_x,ninj_y,tot_re,[0]*len(ninj))
#xsexact,ysexact,qsexact_im,qsexacterr=grab_slice(x0,y0,xsdir,ysdir,psi_x,psi_y,tot_im,[0]*len(tot_im))

#xsexact=np.array(xsexact)
#ysexact=np.array(ysexact)

# set up figure environment
fig = plt.figure()

#for 1D plot
###############################
if '1D' in dim_flag:
    ax=fig.add_subplot(111)
    #ax.set_xticks(np.arange(xi.min(),xi.max()+1,1))
    #ax.set_xlim(-10,10)
    c=['blue','red','green','orange','yellow']
    m=['s','o','^','d','*']

    #xnew=np.linspace(xsexact.min(),xsexact.max(),100)
    #ynew=np.linspace(ysexact.min(),ysexact.max(),100)
    #xnewfunc=xnew/(1.0*Nl[0])
    #ynew=np.empty(len(xnew))
    #ynew.fill(y0)
    #ynewfunc=ynew/(1.0*Nl[1])
    #ax.plot(xnew*np.sqrt(2.0),plot_func_trunc.ev(xnewfunc,ynewfunc))

    for ind in range(num_files):
        if ('ninj' in plot_flag):
            x0=ninj_x[ind,:Ntot_vec[ind]].min()
            y0=ninj_y[ind,:Ntot_vec[ind]].min()*0+1
            xsdir=1
            ysdir=0
            xsexact,ysexact,qsexact_sym,qsexacterr_sym=grab_slice(x0,y0,xsdir,ysdir,ninj_x[ind,:Ntot_vec[ind]],ninj_y[ind,:Ntot_vec[ind]],ninj_sym_restored_1D[ind,:Ntot_vec[ind]],[0]*len(ninj_sym_restored_1D[ind,:Ntot_vec[ind]]))
            xsexact,ysexact,qsexact,qsexacterr=grab_slice(x0,y0,xsdir,ysdir,ninj_x[ind,:Ntot_vec[ind]],ninj_y[ind,:Ntot_vec[ind]],ninj[ind,:Ntot_vec[ind]],[0]*len(ninj[ind,:Ntot_vec[ind]]))
        elif ('didj' in plot_flag):
            x0=didj_x[ind,:Ntot_vec[ind]].min()
            y0=didj_y[ind,:Ntot_vec[ind]].min()*0+1
            xsdir=1
            ysdir=0
            #if(ind==0):
            #    xsexact,ysexact,qsexact,qsexacterr=grab_slice(x0,y0,xsdir,ysdir,didj_x[ind,:Ntot_vec[ind]],didj_y[ind,:Ntot_vec[ind]],didj_sym[ind,:Ntot_vec[ind]],[0]*len(didj_sym[ind,:Ntot_vec[ind]]))
            #if(ind>0):
            #    xsexact,ysexact,qsexact,qsexacterr=grab_slice(x0,y0,xsdir,ysdir,didj_x[ind,:Ntot_vec[ind]],didj_y[ind,:Ntot_vec[ind]],didj[ind,:Ntot_vec[ind]],[0]*len(didj[ind,:Ntot_vec[ind]]))
            xsexact,ysexact,qsexact_sym,qsexacterr_sym=grab_slice(x0,y0,xsdir,ysdir,didj_x[ind,:Ntot_vec[ind]],didj_y[ind,:Ntot_vec[ind]],didj_sym_restored_1D[ind,:Ntot_vec[ind]],[0]*len(didj_sym_restored_1D[ind,:Ntot_vec[ind]]))
            xsexact,ysexact,qsexact,qsexacterr=grab_slice(x0,y0,xsdir,ysdir,didj_x[ind,:Ntot_vec[ind]],didj_y[ind,:Ntot_vec[ind]],didj[ind,:Ntot_vec[ind]],didjerr[ind,:Ntot_vec[ind]])#[0]*len(didj[ind,:Ntot_vec[ind]]))
        elif ('psi_r' in plot_flag):
            x0=psi_x[ind,:Ntot_vec[ind]].min()
            y0=psi_y[ind,:Ntot_vec[ind]].min()*0+1
            xsdir=1
            ysdir=0
            xsexact,ysexact,qsexact,qsexacterr=grab_slice(x0,y0,xsdir,ysdir,psi_x[ind,:Ntot_vec[ind]],psi_y[ind,:Ntot_vec[ind]],mod_psi_r_s[ind,:Ntot_vec[ind]],[0]*len(mod_psi_r_s[ind,:Ntot_vec[ind]]))
            xsexact,ysexact,qsexact_re,qsexacterr=grab_slice(x0,y0,xsdir,ysdir,psi_x[ind,:Ntot_vec[ind]],psi_y[ind,:Ntot_vec[ind]],sing_re[ind,:Ntot_vec[ind]],[0]*len(mod_psi_r_s[ind,:Ntot_vec[ind]]))
            xsexact,ysexact,qsexact_im,qsexacterr=grab_slice(x0,y0,xsdir,ysdir,psi_x[ind,:Ntot_vec[ind]],psi_y[ind,:Ntot_vec[ind]],sing_im[ind,:Ntot_vec[ind]],[0]*len(mod_psi_r_s[ind,:Ntot_vec[ind]]))
            qsexact_re=np.array(qsexact_re)
            qsexact_im=np.array(qsexact_im)

        xsexact=np.array(xsexact)
        ysexact=np.array(ysexact)
        xnew=np.linspace(xsexact.min(),xsexact.max(),100)
        ynew=np.linspace(ysexact.min(),ysexact.max(),100)
        #ynew=np.empty(len(xnew))
        #ynew.fill(y0)
        ax.plot(xnew,plt_fn_dict[ind].ev(ynew,xnew),color=c[ind])
        #ax.plot(np.sqrt(2.0)*xnew,plt_fn_re_dict[ind].ev(ynew,xnew),color=c[ind])
        #ax.plot(np.sqrt(2.0)*xnew,plt_fn_im_dict[ind].ev(ynew,xnew),color=c[ind+1])
        ax.plot(xnew,np.sqrt(plt_fn_re_dict[ind].ev(ynew,xnew)**2+plt_fn_im_dict[ind].ev(ynew,xnew)**2),color=c[ind+1])
    # caution : RectBivariateSpline flips x & y axes, so must enter arguments as (y,x), not (x,y)
    #ax.plot(xnew,plot_func.ev(ynew,xnew),color='orange')
    #ax.plot(xnew,plot_func_re.ev(ynew,xnew),color='blue')
    #ax.plot(xnew,plot_func_im.ev(ynew,xnew),color='green')
    #ax.plot(xnew,np.sqrt(plot_func_re.ev(ynew,xnew)**2+plot_func_im.ev(ynew,xnew)**2),color='purple')

        #xsexact=np.sqrt(2.0)*xsexact
    #ax.plot(xsexact,qsexact,linestyle='None',color='red',marker='o')

        # plot with different labels
        #ax.plot(xsexact,qsexact,linestyle='None',color=c[ind],marker=m[ind],label=r'L $\,=\,$'+str(int(np.sqrt(Ntot_vec[ind]))))
        #ax.plot(xsexact,qsexact,linestyle='None',color=c[ind],marker=m[ind],label=r'$\lambda =\,$'+str(lam[ind]))
        ax.plot(xsexact,np.sqrt(qsexact_re**2+qsexact_im**2),linestyle='None',color=c[ind],marker=m[ind],label=r'$\lambda =\,$'+str(lam[ind]))
        # real and imaginary parts
        #ax.plot(np.sqrt(2.0)*xsexact,qsexact_re,linestyle='None',color=c[ind],marker=m[ind],label=r'$\lambda =\,$'+str(lam[ind]))
        #ax.plot(np.sqrt(2.0)*xsexact,qsexact_im,linestyle='None',color=c[ind+1],marker=m[ind+1],label=r'$\lambda =\,$'+str(lam[ind]))
        #if ind==2:
        #    print qsexact
        #    print qsexacterr
        #ax.plot(np.sqrt(2.0)*xsexact,qsexact,linestyle='None',color=c[ind],marker=m[ind],label=r'U $\,=\,$'+str(Uhubbard[ind]))

        #ax.plot(xsexact,qsexact_sym,linestyle='None',color=c[ind+1],marker=m[ind+1])
    #ax.plot(xsexact,qsexact_im,linestyle='None',color='cyan',marker='^')
    plt.subplots_adjust(bottom=0.15)
    ax.set_xticks(np.arange(xsexact.min(),xsexact.max()+1,1))
    #ax.set_xlim((xsexact.min()-1,xsexact.max()+1))
    ax.set_xlim(-8,8)
    #ax.set_ylim(0.485,0.515)
    ax.legend(loc='upper right',numpoints=1)

    show()
    #plotname=str(plot_flag)+"_L"+str(int(np.sqrt(Ntot_vec[0])))+"_U"+str(Uhubbard[0])+"_vs_lbda_along_x=y_x0=-5_y0=-7.pdf"
    #print xsexact
    #print ysexact
    #print plotname
    #fig.savefig(plotname)

#for 2D plot
###############################
elif '2D' in dim_flag:
    for ind in range(num_files):
        xi=xi_dict[ind]
        yi=yi_dict[ind]
        #ax=fig.add_subplot(111)
        fig, ax = plt.subplots(1,1)
        #ax.set_xticks(np.arange(xi.min(),xi.max()+1,1))
        #ax.set_yticks(np.arange(yi.min(),yi.max()+1,1))
        #ax.set_xlim(-10,10)
        #ax.set_ylim(-10,10)
        ax.set_xlabel('x',fontsize=24)
        ax.set_ylabel('y',fontsize=24,rotation=360)
        ax.xaxis.set_tick_params(labelsize=22)
        ax.yaxis.set_tick_params(labelsize=22)

        # caution : RectBivariateSpline flips x & y axes, so must enter arguments as (y,x), not (x,y)
        if (('psi_r' not in plot_flag) or ('tot' in plot_flag)):
            #implot = ax.imshow(plot_func(yi,xi),extent=(xi.min(),xi.max(),yi.min(),yi.max()),cmap=cm.gnuplot2,origin='lower')
            implot = ax.imshow(plt_fn_dict[ind](yi,xi),extent=(xi.min(),xi.max(),yi.min(),yi.max()),cmap=cm.gnuplot2,origin='lower')
        else:
            implot = ax.imshow(np.sqrt(plt_fn_re_dict[ind](yi,xi)**2+plt_fn_im_dict[ind](yi,xi)**2),extent=(xi.min(),xi.max(),yi.min(),yi.max()),cmap=cm.gnuplot2,origin='lower')
            #implot = ax.imshow(np.sqrt(plot_func_re(yi,xi)**2+plot_func_im(yi,xi)**2),extent=(xi.min(),xi.max(),yi.min(),yi.max()),cmap=cm.gnuplot2,origin='lower')

        cb=fig.colorbar(implot,shrink=1.0,pad=0.01)
        cb.ax.tick_params(labelsize=22)

        plt.subplots_adjust(bottom=0.15)

        save=0
        if('params' in params_flag):
            save_flag=s_flag[ind]
            lbda=lam[ind]
            Uhubb=Uhubbard[ind]
        if(('y' in save_flag) or ('Y' in save_flag)):
            save=1
        for arg in sys.argv:
            if ('save' in str(arg)) or ('Save' in str(arg)):
                save=1
                break

        if save==1:
            if 'corr_flag' in locals():
                if (('y' in o_flg[ind]) or ('Y' in o_flg[ind])):
                    #print sym_flag_vec[ind]
                    if 'sym' in sym_flag_vec[ind]:
                        plotname=str(corr_flag)+'_Lx'+str(Nl[0])+'_Ly'+str(Nl[1])+'_Np'+str(Np)+'_U'+str(Uhubb)+'_l'+str(lbda)+'_'+str(dim_flag)+'_open_BCs_sym_restored'+'.pdf'
                    else:
                        plotname=str(corr_flag)+'_Lx'+str(Nl[0])+'_Ly'+str(Nl[1])+'_Np'+str(Np)+'_U'+str(Uhubb)+'_l'+str(lbda)+'_'+str(dim_flag)+'_open_BCs.pdf'
                elif('MF' in params_flag):
                    plotname=str(corr_flag)+'_Lx'+str(Nl[0])+'_Ly'+str(Nl[1])+'_Np'+str(Np)+'_U'+str(Uhubb)+'_l'+str(lbda)+'_'+str(dim_flag)+'_MF.pdf'
                else:
                    plotname=str(corr_flag)+'_Lx'+str(Nl[0])+'_Ly'+str(Nl[1])+'_Np'+str(Np)+'_U'+str(Uhubb)+'_l'+str(lbda)+'_'+str(dim_flag)+'.pdf'
                print plotname
                fig.savefig(plotname)
            else:
                if (('y' in o_flg[ind]) or ('Y' in o_flg[ind])):
                    if 'sym' in sym_flag_vec[ind]:
                        plotname=str(plot_flag)+'_Lx'+str(Nl[0])+'_Ly'+str(Nl[1])+'_Np'+str(Np)+'_U'+str(Uhubb)+'_l'+str(lbda)+'_'+str(dim_flag)+'_open_BCs_sym_restored'+'.pdf'
                    else:
                        plotname=str(plot_flag)+'_Lx'+str(Nl[0])+'_Ly'+str(Nl[1])+'_Np'+str(Np)+'_U'+str(Uhubb)+'_l'+str(lbda)+'_'+str(dim_flag)+'_open_BCs.pdf'
                elif('MF' in params_flag):
                    plotname=str(plot_flag)+'_Lx'+str(Nl[0])+'_Ly'+str(Nl[1])+'_Np'+str(Np)+'_U'+str(Uhubb)+'_l'+str(lbda)+'_'+str(dim_flag)+'_MF.pdf'
                else:
                    plotname=str(plot_flag)+'_Lx'+str(Nl[0])+'_Ly'+str(Nl[1])+'_Np'+str(Np)+'_U'+str(Uhubb)+'_l'+str(lbda)+'_'+str(dim_flag)+'.pdf'
                print plotname
                fig.savefig(plotname)
        else:
            #show()
            plt.waitforbuttonpress()
            plt.cla()        

#for 3D plot
###############################
elif '3D' in dim_flag:
    ax=fig.gca(projection = '3d')

    # fix aspect ratio if lattice isn't square
    # otherwise use ratio of 1.0

    if(abs(Nl[0]-Nl[1])>1e-9):
        ax.pbaspect = [0.75, 1.5,0.5]
    elif(abs(Nl[0]-Nl[1])<1e-9):
        ax.pbaspect = [1.0,1.0,1.0]
        #ax.pbaspect = [1.0,1.0,0.85]

        # axis labels and limits
        #ax.set_xlim(-Nl[0]/2,Nl[0]/2-1)
        #ax.set_ylim(-Nl[1]/2,Nl[1]/2-1)
        ax.set_xlabel('x')
        ax.set_ylabel('y')

        # caution : RectBivariateSpline flips x & y axes, so must enter arguments as (y,x), not (x,y)

        #surf = ax.plot_surface(xi, yi, plot_func(yi,xi), rstride=1, cstride=1, cmap=cm.gnuplot,
        #                       linewidth=0, antialiased=False)

        if (('psi_r' not in plot_flag) or ('tot' in plot_flag)):
            surf = ax.plot_surface(xi, yi, plot_func(yi,xi), rstride=1, cstride=1, cmap=cm.gnuplot,
                               linewidth=0, antialiased=False)
        else:
            surf = ax.plot_surface(xi, yi, np.sqrt(plot_func_re(yi,xi)**2+plot_func_im(yi,xi)**2), rstride=1, cstride=1, cmap=cm.gnuplot,
                               linewidth=0, antialiased=False)
            #ax.scatter(psi_x, psi_y, mod_dd, s=20, c='cyan')

        fig.colorbar(surf,shrink=0.75,pad=0.01)

"""
plt.subplots_adjust(bottom=0.15)


save=0
if('params' in params_flag):
    save_flag=str(paramfile.readline().split(":")[1])
if(('y' in save_flag) or ('Y' in save_flag)):
    save=1
for arg in sys.argv:
    if ('save' in str(arg)) or ('Save' in str(arg)):
        save=1
        break

if save==1:
    if 'corr_flag' in locals():
        if 'open_flg' in locals():
            if 'sym' in sym_flag:
                plotname=str(corr_flag)+'_Lx'+str(Nl[0])+'_Ly'+str(Nl[1])+'_Np'+str(Np)+'_U'+str(Uhubb)+'_l'+str(lbda)+'_'+str(dim_flag)+'_'+str(open_flg)+'_sym_restored'+'.pdf'
            else:
                plotname=str(corr_flag)+'_Lx'+str(Nl[0])+'_Ly'+str(Nl[1])+'_Np'+str(Np)+'_U'+str(Uhubb)+'_l'+str(lbda)+'_'+str(dim_flag)+'_'+str(open_flg)+'.pdf'
        else:
            plotname=str(corr_flag)+'_Lx'+str(Nl[0])+'_Ly'+str(Nl[1])+'_Np'+str(Np)+'_U'+str(Uhubb)+'_l'+str(lbda)+'_'+str(dim_flag)+'.pdf'
        print plotname
        fig.savefig(plotname)
    else:
        if 'open_flg' in locals():
            if 'sym' in sym_flag:
                plotname=str(plot_flag)+'_Lx'+str(Nl[0])+'_Ly'+str(Nl[1])+'_Np'+str(Np)+'_U'+str(Uhubb)+'_l'+str(lbda)+'_'+str(dim_flag)+'_'+str(open_flg)+'_sym_restored'+'.pdf'
            else:
                plotname=str(plot_flag)+'_Lx'+str(Nl[0])+'_Ly'+str(Nl[1])+'_Np'+str(Np)+'_U'+str(Uhubb)+'_l'+str(lbda)+'_'+str(dim_flag)+'_'+str(open_flg)+'.pdf'
        else:
            plotname=str(plot_flag)+'_Lx'+str(Nl[0])+'_Ly'+str(Nl[1])+'_Np'+str(Np)+'_U'+str(Uhubb)+'_l'+str(lbda)+'_'+str(dim_flag)+'.pdf'
        print plotname
        fig.savefig(plotname)
else:
    show()        
"""

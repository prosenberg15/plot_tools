from scipy import interpolate
import numpy as np
import copy

# convert lattice site to x,y,z coords
def coor(sitenum, Nl):
    if sitenum%Nl[0]!=0:
        x=sitenum%Nl[0]-1
    else:
        x=(sitenum-1)%Nl[0]
    if len(Nl)==1:
        return x
    if len(Nl)>1:
        if sitenum%Nl[1]!=0:
            y=int(sitenum)/int(Nl[1])
        else:
            y=int(sitenum)/int(Nl[1])-1
    if len(Nl)==2:
        return x,y
    if len(Nl)==3:
        if sitenum%(Nl[1]*Nl[2])!=0:
            z=sitenum%(Nl[1]*Nl[2])
        else:
            z=sitenum%(Nl[1]*Nl[2])-1
            return x,y,z

# grab slice of 2D data
def grab_slice(x0, y0, xslicedir, yslicedir, x, y, q, qerr):
    val=[]
    valerr=[]
    xs=[]
    ys=[]
    for i in range(len(x)):
        xnew=x0+i*xslicedir
        ynew=y0+i*yslicedir
        if xnew > x.max() or ynew > y.max():
            break
        for j in range(len(q)):
            if x[j]==xnew and y[j]==ynew:
                xs.append(xnew)
                ys.append(ynew)
                val.append(q[j])
                valerr.append(qerr[j])
    return xs,ys,val,valerr

# shift origin of list of coords
# if x >= Lx/2, x --> x-Lx (same for y) 
def shift_origin_list(x,dim):
    xloc=copy.deepcopy(x)
    for i in range(len(xloc)):
        if xloc[i]>=dim/2:
            xloc[i]=xloc[i]-dim
    return xloc

# shift origin of 2D array
# if x >= Lx/2, x --> x-Lx (same for y) 
def shift_origin_array(z,dims):
    if(type(z).__module__ != np.__name__):
        print "input must be numpy array, automatically converting"
        zloc=copy.deepcopy(np.array(z))
    else:
        zloc=copy.deepcopy(z)
    for i in range(z.shape[0]):
        if(dims[0]>2): 
            zloc[i,:]=np.roll(zloc[i,:],dims[0]/2)
            #if i==0:
            #    print i, z[i,:]
    for j in range(z.shape[1]):
        if(dims[1]>2):
            zloc[:,j]=np.roll(zloc[:,j],dims[1]/2,axis=0)    
    return zloc

# reshape data in 1D array to 2D
# array corresponding to lattice sites
# to be plotted
def reshape_1D_to_2D(x, y, data):
    xnew, ynew = np.ogrid[x.min():x.max():(x.max()-x.min()+1)*1j, y.min():y.max():(y.max()-y.min()+1)*1j]
    dataloc=copy.deepcopy(data)
    dataloc=np.reshape(dataloc,(ynew.shape[1],xnew.shape[0]))    
    return dataloc

# perform spline-interpolation on 2D data
def splineinterp_2D(x, y, data,degx,degy,smooth):
    xnew, ynew = np.ogrid[x.min():x.max():(x.max()-x.min()+1)*1j, y.min():y.max():(y.max()-y.min()+1)*1j]   
    # caution : RectBivariateSpline flips x & y axes, must provides parameters as (y,x), not (x,y)
    return interpolate.RectBivariateSpline(ynew,xnew,data,kx=degx,ky=degy,s=smooth)

"""
# take raw data in list form transform to 2D array form output interpolating function
def perform_interp(x, y, data, degx, degy, smooth, shift, dims):
    data=reshape_1D_to_2D(x,y,data)
    if shift:
        data=shift_origin_array(data,dims)
    return splineinterp_2D(x,y,data,degx,degy,smooth)
"""     
        
        
        
        
        
        
        
        
        
        
    

import os
import numpy as np
from scipy.interpolate import interp1d
import pandas as pd
import astropy.units as u
import matplotlib as mpl
import matplotlib.pyplot as plt
from ctypes import *

lsize = 16
mpl.rc('xtick',direction='in',top=True,labelsize=lsize) 
mpl.rc('ytick',direction='in',right=True,labelsize=lsize) 
mpl.rc('font',family='Times New Roman')
params = {'text.usetex': False, 'mathtext.fontset': 'stix',\
          'axes.labelsize': lsize, 'legend.fontsize': lsize, 'axes.grid' : True}
plt.rcParams.update(params)

# Loading functions in C++ library
array_1d_double  = np.ctypeslib.ndpointer(dtype=np.double, ndim=1,flags="CONTIGUOUS")
libdir = os.path.dirname(os.path.realpath(__file__))
try:
    libobj = np.ctypeslib.load_library("libobj.so",libdir)
    C_available = True
    libobj.kinematic_model.restype = POINTER(c_double)
    libobj.kinematic_model.argtypes = [array_1d_double,array_1d_double,c_int,array_1d_double,\
                                       c_char_p,array_1d_double,c_double,c_double,c_int]
except:
    print ("Cannot find C library, C functions will not be available")
    C_available = False


# Radius and rotation velocity of the Sun, defined globally
RSun = 8.2      # kpc
VSun = 240.     # km/s

##########################################################################
# UTILITY CLASSES AND FUNCTIONS
##########################################################################
class Sightlines:
    def __init__(self):
        self.n = 0
        self.lon, self.lat = np.array([]), np.array([])
        self.vlsr, self.vgsr = np.array([]), np.array([])
        self.logN, self.spec_grid, self.spec = np.array([]), [],[]
    
    def add_sightlines(self,l,b,vlsr,vgsr,logN):
        self.lon  = np.append(self.lon,l)
        self.lat  = np.append(self.lat,b)
        self.vlsr = np.append(self.vlsr,vlsr)
        self.vgsr = np.append(self.vgsr,vgsr)
        self.logN = np.append(self.logN,logN)
        self.n = len(self.lon)


def checkSize(a,b):
    if (len(a)!=len(b)):
        raise ValueError("Different sizes")
        

def reshapePointer (p, shape):
    """Take a POINTER to c_float and reshape it.
    
    Args:
      p (POINTER(c_float)): The pointer to be reshaped
      shape (tuple, list): The new shape of the array
    
    Returns:
      ndarray: The reshaped array
    
    """
    return np.ctypeslib.as_array(p, shape=tuple(shape))


def galactic2cylindrical(l,b,D):
    # l, b are in radians here, D in kpc
    z  = D*np.sin(b)
    dp = D*np.cos(b)
    x  = RSun - dp*np.cos(l)
    y  = dp*np.sin(l)
    R  = np.sqrt(x**2+y**2)
    th = np.arctan2(y,x)
    return (R,th,z)


##########################################################################
# MODELS FOR GALACTIC ROTATION
##########################################################################
def rotationModel(R,z,pars):
    # Simple flat rotation model with a vertical lag
    # (R,z) in kpc, vflat in km/s, lag in km/s/kpc
    vflat, lag = pars
    return vflat - lag*np.abs(z)


##########################################################################
# MODELS FOR GAS DENSITY DISTRIBUTIONS
##########################################################################
def Constant_Density(R, z, pars):
    # Density field is the same everywhere: dens(R,z) = const
    return pars[0];

def VerticalExponential_Density(R,z,pars):
    # Density is constant in R, but drops exponentially in z: 
    #    dens(R,z) = rho0*exp(-|z|/z0)
    # where rho0=pars[0] in cm^-3, z0=pars[1] in kpc
    rho0, z0 = pars
    return rho0*np.exp(-np.abs(z)/z0)

def RadialVerticalExponential_Density(R,z,pars):
    # Density drops exponentially in both z and R: 
    #    dens(R,z) = rho0*exp(-R/R0)*exp(-|z|/z0)
    # where rho0=pars[0] in cm^-3, R0=pars[1] in kpc 
    # and z0=pars[2] in kpc
    rho0, R0, z0 = pars
    return rho0*np.exp(-R/R0)*np.exp(-np.abs(z)/z0)
    
def ThickSandwich_Density(R,z,pars):
    # Density is constant in a thick layer at some height:
    #    dens(R,z) = rho0 if hmin<z<hmax else 0
    # where rho0=pars[0] in cm^3, 
    # hmin=pars[1] and hmax=pars[2] in kpc
    rho0, hmin, hmax = pars
    dens = np.full(len(z),rho0)
    dens[np.abs(z)<hmin]=0
    dens[np.abs(z)>hmax]=0
    return dens
    
def GaussianSandwich_Density(R,z,pars):
    # Density is Gaussian in a layer at some height:
    #    dens(R,z) = rho0*exp(-(z-|z0|)**2/(2*sigma**2))
    # where rho0=pars[0] in cm^-3, z0=pars[1] in kpc, sigma=pars[2] in kpc
    rho0, h0, sigma = pars
    dens = rho0*(np.exp(-(z-h0)**2/(2*sigma**2))+np.exp(-(z+h0)**2/(2*sigma**2)))
    return dens
    
def FlatSandwich_Density(lon,lat,velopars,denspars):
            
    rotmodpars = velopars[0:2]
    vR = velopars[2]
    vz = np.full(len(lon),velopars[3])
    vz[lat<0] *= -1
    n0, h0 = denspars[0:2]
    
    l  = np.radians(lon)
    b  = np.radians(lat)
    D  = h0/np.sin(np.abs(b))
    # Calculating cylindrical coordinates
    R,th,z = galactic2cylindrical(l,b,D)

    vtheta = rotationModel(R,z,rotmodpars)

    vlsrs = (RSun*np.sin(l))*(-VSun/RSun+vtheta/R)*np.cos(b) - \
             vR*np.cos(l+th)*np.cos(b) + vz*np.sin(b)

    vgsrs = vlsrs + VSun*np.sin(l)*np.cos(b)
    
    sightlines = Sightlines()
    sightlines.add_sightlines(lon,lat,vlsrs,vgsrs,np.full(len(lon),np.log10(n0)))
    return sightlines


##########################################################################
# DEFINE A MODEL GIVEN DENSITY AND VELOCITY FIELD
##########################################################################
def kinematic_model(lon,lat,velopars,densmodel,denspars,useC=True,nthreads=8,getSpectra=False):
    
    '''
    This function calculate the predicted column-density weighted VLSR
    along a sightline (l,b) for a given density and velocity model.
     
    Args:
     lon (float,tuple):    Longitude(s) for which the model will be calculated
     lat (float,tuple):    Latitude(s) for which the model will be calculated
     velopars (tuple):     List of dim=4 with parameters for velocity field
                           (vflat,lag,vR,vz) where vflat is a flat rotation velocity in km/s,
                           lag is the vertical lag in rotation in km/s/kpc, vR and vz are
                           the radial and vertical velocities in km/s.
     densmodel (func):     Type of density field. Must be one of above density functions
     denspars (tuple):     Parameters for the density field. See above the density functions.
     useC (bool):          Whether to use C implementation (default=True)
     nthreads (int):       Number of CPU to use in the C implementation (default=8)
     getSpectra (bool):    Whether to return the full spectrum along a sightline (default=False)
    
    Returns:
     A Sightline() object containing a list of sightlines 
     with vlsrs, vgsrs and logN for the model
    '''

    if isinstance(lon,(int,float)): lon = [lon]
    if isinstance(lat,(int,float)): lat = [lat]
    checkSize(lon,lat)
    lon = np.array(lon).astype(float)
    lat = np.array(lat).astype(float)
    
    if densmodel==FlatSandwich_Density:
        # This is treated as a special case just because it is faster 
        return FlatSandwich_Density(lon,lat,velopars,denspars)
   
    # These are the returned sightlines
    sightlines = Sightlines()
    
    if getSpectra: useC = False
    
    if useC and C_available:
        ret = libobj.kinematic_model(lon,lat,len(lon),np.array(velopars).astype(float),\
                                     densmodel.__name__.encode('utf-8'),
                                     np.array(denspars).astype(float),RSun,VSun,int(nthreads))
        ret   = reshapePointer(ret,[2*len(lon)])
        vlsrs = ret[:len(lon)]
        logNs = ret[len(lon):] 
        vgsrs = vlsrs + VSun*np.sin(np.radians(lon))*np.cos(np.radians(lat)) 
        sightlines.add_sightlines(lon,lat,vlsrs,vgsrs,logNs)
    else:
        D_min, D_max = 0.0001, 50 # kpc
        deltaD = 0.01
        Ds = np.arange(D_min,D_max,deltaD)
        N = np.zeros(len(Ds))
    
        vgrid = np.arange(-300,300,1)
        sightlines.spec_grid = vgrid
        
        rotmodpars = velopars[0:2]
        vR = velopars[2]

        for i in range(len(lon)):
            
            l, b  = np.radians(lon[i]), np.radians(lat[i])
            
            # Calculating cylindrical coordinates
            R,th,z = galactic2cylindrical(l,b,Ds)
            
            # Getting the density model
            rho = densmodel(R,z,pars=denspars)
            # Calculating column densities
            N[:len(rho)-1] = 0.5*(rho[:len(rho)-1]+rho[1:])*(deltaD*u.kpc).to('cm').value
            
            # Getting the rotation model
            vtheta = rotationModel(R,z,pars=rotmodpars)
            vz = velopars[3] if lat[i]>0 else -velopars[3];
            
            # Calculating observed VLSR and VGSR
            vlsrs = (RSun*np.sin(l))*(-VSun/RSun+vtheta/R)*np.cos(b) - \
                     vR*np.cos(l+th)*np.cos(b) + vz*np.sin(b)

            vgsrs = vlsrs + VSun*np.sin(l)*np.cos(b) 

            # Calculating N-weighted average velocities
            vlsr_av = vgsr_av = np.nan
            try:
                vlsr_av = np.average(vlsrs,weights=N)
                vgsr_av = np.average(vgsrs,weights=N)
            except:
                pass
            Nsum = np.nansum(N)
            print("Min,Max,Median Sums: %7.1f, %7.f, %7.1f" % (min(Nsum), max(Nsum), np.median(Nsum)))

            logN = np.log10(Nsum)
        
            sightlines.add_sightlines(lon[i],lat[i],vlsr_av,vgsr_av,logN)

            # Calculating full spectrum if requested
            spectrum = None
            if getSpectra:
                vdisp = 10
                gaussians = N/1E18*np.exp(-(vgrid[:,None]-vlsrs[None,:])**2/(2*vdisp**2))
                spectrum = np.nansum(gaussians,axis=1)
                sightlines.spec.append(spectrum)
        
    return sightlines
    


##########################################################################
# PLOTTING FUNCTIONS
##########################################################################
def plot_model(slines):
    
    fig, ax = plt.subplots(figsize=(20,20),nrows=1, ncols=2, subplot_kw= {'projection' : 'aitoff'})
    fig.subplots_adjust(hspace=0.15, wspace=0.08)
    ax = np.ravel(ax)

    cmap = plt.get_cmap('coolwarm')
    vels = [slines.vlsr,slines.vgsr]
    labs = [r'$V_\mathrm{LSR}$ (km/s)',r'$V_\mathrm{GSR}$ (km/s)']

    for i in range(len(ax)):
        a = ax[i]
        vmax = np.nanmax(np.abs(vels[i]))
        norm = mpl.colors.Normalize(vmin=-vmax,vmax=vmax)

        a.grid(True)
        a.axhline(0,c='k')
        a.scatter(np.radians(slines.lon),np.radians(slines.lat),c=vels[i],cmap=cmap, norm=norm, rasterized=True)
        cbax = fig.add_axes([a.get_position().x0,a.get_position().y1+0.03,\
                             a.get_position().x1-a.get_position().x0,0.02]) 
        cb = mpl.colorbar.ColorbarBase(ax=cbax,cmap=cmap,norm=norm,\
                                       orientation='horizontal',ticklocation='top')
        cb.set_label(labs[i],labelpad=10)
    return (fig, ax)


def plot_datavsmodel(data,model):
    
    fig, ax = plt.subplots(figsize=(20,20),nrows=1, ncols=2, subplot_kw= {'projection' : 'aitoff'})
    fig.subplots_adjust(hspace=0.15, wspace=0.05)
    ax = np.ravel(ax)
    cmap = plt.get_cmap('coolwarm')
    vmax = np.nanmax(np.concatenate([data.vlsr,model.vlsr]))
    norm = mpl.colors.Normalize(vmin=-vmax,vmax=vmax)

    for i in range (len(ax)):
        a = ax[i]
        a.grid(True)
        a.axhline(0,c='k')
    
    ax[0].scatter(np.radians(data.lon),np.radians(data.lat),c=data.vlsr,cmap=cmap,norm=norm,edgecolors='gray',s=70)
    ax[0].text(0.5,-0.1,"DATA",transform=ax[0].transAxes,fontsize=lsize+2,ha='center')

    ax[1].scatter(np.radians(model.lon),np.radians(model.lat),c=model.vlsr,cmap=cmap,norm=norm,edgecolors='gray',s=70)
    ax[1].text(0.5,-0.1,"MODEL",transform=ax[1].transAxes,fontsize=lsize+2,ha='center')


    cbax = fig.add_axes([ax[0].get_position().x0,ax[0].get_position().y1+0.03,\
                         ax[1].get_position().x1-ax[0].get_position().x0,0.02]) 
    cb = mpl.colorbar.ColorbarBase(ax=cbax,cmap=cmap,norm=norm,orientation='horizontal',ticklocation='top')
    cb.set_label(r"$V_\mathrm{LSR}$ (km/s)",labelpad=10)
    return (fig, ax)

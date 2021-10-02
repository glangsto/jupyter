import numpy as np
import pandas as pd
from astropy.io import fits
import emcee, corner,time
from multiprocessing import Pool
from scipy.optimize import minimize
from scipy.stats import norm
from DHmodels import *


def lnprior(pars):

    #"""##### FlatSandwich priors ######
    vflat, lag, vz, h0 = pars
    if vflat<0 or lag<0 or h0<0:
        return -np.inf
    # Gaussian priors (pdfs)
    h_prior   = norm.pdf(h0,loc=5,scale=2)
    lag_prior = norm.pdf(lag,loc=10,scale=2)
    vf_prior  = norm.pdf(vflat,loc=240,scale=20)
    vz_prior  = norm.pdf(vz,loc=0,scale=10)
    prior = np.log(h_prior*lag_prior*vz_prior*vf_prior)
    # Flat priors
    #pr = (0.0001<h0<10) and (0<lag<200) and (200<vflat<280) and (-20<vz<20)
    #prior = 0 if pr else -np.inf
    #"""
    
    """##### GaussianSandwich priors ######
    vflat, lag, vz, h0, sigma = pars
    if vflat<0 or lag<0 or h0<0 or sigma<0:
        return -np.inf
    # Gaussian priors (pdfs)
    vf_prior  = norm.pdf(vflat,loc=240,scale=20)
    lag_prior = norm.pdf(lag,loc=10,scale=2)
    vz_prior  = norm.pdf(vz,loc=0,scale=10)
    h_prior   = norm.pdf(h0,loc=5,scale=2)
    sig_prior = norm.pdf(sigma,loc=1,scale=0.5)
    prior = np.log(h_prior*lag_prior*vz_prior*vf_prior*sig_prior)
    """

    return prior


def lnlike(pars,data):
    
    ##### FlatSandwich parameters ######
    vf, lag, vz, h0 = pars
    densmod  = FlatSandwich_Density
    velopars = (vf,lag,0.,vz)
    denspars = (1E-05,h0)
    
    ##### GaussianSandwich parameters ######
    #vf, lag, vz, h0, sigma = pars
    #densmod  = GaussianSandwich_Density
    #velopars = (vf,lag,0.,vz)
    #denspars = (1E-05,h0,sigma)

    mod = kinematic_model(data.lon,data.lat,velopars=velopars,densmodel=densmod,\
                          denspars=denspars,useC=True,nthreads=4)
    diff = np.nansum((mod.vlsr-data.vlsr)**2)
    return -diff


def lnprob(pars,data):
    lp = lnprior(pars)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(pars,data)



if __name__ == '__main__':

    ###########################################################################
    # FlatSandwich
    p0     = [230, 15, -5, 1]               # Initial guesses
    labels = ["vflat", "lag","vz", "h0"]    # Names of parameters to fit
    
    # GaussianSandwich 
    #p0 = [230, 15, -5., 5., 0.5]
    #labels = ["vflat", "lag", "vz", "h0","sigma"]
    ###########################################################################

    #Reading in sightlines
    ds = pd.read_table("data/sightlines_flag_2.txt", sep=' ', skipinitialspace=True)

    # Here we chose which ion we want to fit
    ion = 'CIV'
    di = ds[ds['ion']==ion]

    # We select only latitudes below 60 deg for the fit
    glon, glat, vl = di['Glon'].values, di['Glat'].values, di['weighted_v_LSR'].values
    glon[glon>180] -= 360
    m = (np.abs(glat)<60)
    # Just storing the data in a Sightline object for convenience
    data = Sightlines()
    data.add_sightlines(glon[m],glat[m],vl[m],None,None)
    print (f"Sightlines to fit: {len(data.lon)}")

    nll = lambda *args: -lnprob(*args)
    
    # This is just to minimize the likelihood in a classical way
    soln = minimize(nll, p0, args=(data),method='Nelder-Mead')
    print ("Best-fit parameters:", soln.x)
    
    
    # Initializing chains and walkers
    ndim, nwalkers, nsteps = len(p0), 500, 500
    pos = p0 + 1e-1*np.random.randn(nwalkers,ndim)

    print ("\n Running MCMC...")
    with Pool() as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[data], pool=pool)
        start = time.time()
        res = sampler.run_mcmc(pos, nsteps, progress=True)
        multi_time = time.time()-start
        print("Computational time {0:.1f} minutes".format(multi_time/60.))

    # Saving samples
    fits.writeto(f"{ion}_samples.fits",np.float32(sampler.get_chain()),overwrite=True)

    # Plotting chains 
    fig, axes = plt.subplots(ndim, figsize=(10, 7), sharex=True)
    samples = sampler.get_chain()
    for i in range(ndim):
        ax = axes[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i])
        ax.yaxis.set_label_coords(-0.1, 0.5)

    axes[-1].set_xlabel("step number")
    fig.savefig(f"{ion}_chains.png")

    # Burn-in
    burnin = 100
    thin = 1
    samples = sampler.get_chain(discard=burnin, thin=thin, flat=True)

    print ("\n MCMC parameters:")
    pp = []
    for i in range(ndim):
        mcmc = np.percentile(samples[:, i], [15.865, 50, 84.135])
        q = np.diff(mcmc)
        txt = "%10s = %10.3f %+10.3f %+10.3f"%(labels[i],mcmc[1], -q[0], q[1])
        print (txt)
        pp.append(mcmc[1])
    
    # Autocorrelation function
    #tau = sampler.get_autocorr_time(quiet=True)
    #burnin = int(2*np.nanmax(tau))
    #thin = int(0.5*np.nanmin(tau))
    #print(tau,burnin,thin)

    # Save corner plot
    # Levels 
    levels = 1.0-np.exp(-0.5*np.arange(0.5, 2.1, 3.5) ** 2)
    #levels = 1.0-np.exp(-0.5*np.array([1., 2., 3.]) ** 2)

    fig = corner.corner(samples, truths=pp, labels=labels, show_titles=True, title_kwargs={"fontsize": lsize},\
                truth_color='firebrick') #,fill_contours=True,levels=levels)
    fig.savefig(f"{ion}_corner.pdf",bbox_inches='tight')

    # Plot model vs data
    # NB: you need to change the line below when changing parameters/model
    # FlatSandwich
    model = kinematic_model(data.lon,data.lat,velopars=(pp[0],pp[1],0,pp[2]),densmodel=FlatSandwich_Density,\
                            denspars=(1E-08,pp[3]),useC=True,nthreads=8,getSpectra=False)
    # GaussianSandwich
    #model = kinematic_model(data.lon,data.lat,velopars=(pp[0],pp[1],0,pp[2]),densmodel=GaussianSandwich_Density,\
    #                        denspars=(1E-08,pp[3],pp[4]),useC=True,nthreads=8)
    
    fig, ax = plot_datavsmodel(data,model)
    fig.savefig(f"{ion}_comp.pdf",bbox_inches='tight')

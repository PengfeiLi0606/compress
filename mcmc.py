
import numpy as np
import emcee
import h5py
from multiprocessing import Pool
from subprocess import Popen
from scipy import special
import os
os.environ["OMP_NUM_THREADS"] = "1"

def lnprob(x, name, DM, alpha):

    lnp = lnprior(x)
    if not np.isfinite(lnp):
        return -np.inf
    else:
        return lnlike(x, name, DM, alpha)

def lnprior(x):

    V200 = x[0]
    C200 = x[1]
    
    if V200>50. and V200<200. and C200>5.0 and C200<30.:
        return 0.0
    else:
        return -np.inf

def lnlike(x, name, DM, alpha):

    V200 = x[0]
    C200 = x[1]
    ML = 1.0
    MLb = 1.0

# #### Run compression

    rs = V200/C200/0.73 # kpc
    M200 = V200**3/0.73/4.3*10**6 # Msun
    if DM=='NFW':     # hmass = 4 pi rho_s r_s^3
        hmass = M200/( np.log(1+C200) - C200/(1+C200) )/10**10
    elif DM=='EINA':  # hmass = 16 pi rho_s r_s^3
        hmass = M200*4/np.exp(2./alpha) / (2/alpha)**(-3./alpha) *alpha / special.gammainc(3/alpha, 2/alpha*C200**alpha)/special.gamma(3/alpha)/10**10
    Process=Popen('./compress %s %s %s %s %s' % (str(name), str(hmass), str(rs), str(ML), str(MLb)), shell=True )
    Process.wait()

    outputname = name+'_'+str(format(ML,'.2f'))+'_'+str(format(MLb,'.2f'))+'_'+str(format(int(hmass*10000)/10000., '.4f'))+'_'+str(format(int(rs*10000)/10000., '.4f'))
    chi_sq, nobs = np.loadtxt('output/'+outputname+'.fitquality', usecols=(0,1), unpack=True)
    Process=Popen( 'rm output/%s.fitquality' %(str(outputname),), shell=True )
    
    return -0.5*chi_sq*nobs

# ## Starts fitting here #####################################################################################

if __name__ == '__main__':

# ## set baryonic models
    name = 'MW_McGaugh19'
    ML, MLb = 1.0, 1.0

# ## set dark matter halo models
    DM='EINA'
    alpha=2.33   # shape parameter for Einasto, can be any value for NFW

# ## set MCMC parameters
    ndim = 2
    nwalkers = 30
    niteration = 500 
    nproc = 15 

# ## Initialize chains    
    p0 = np.random.rand(nwalkers*ndim).reshape((nwalkers, ndim))
    p0[:,0] = p0[:,0]*(200-50)+50. # 50 < V200 < 200
    p0[:,1] = p0[:,1]*(30-5)+5. # 5 < C200 < 30

# ## set sampler and run
    
    with Pool(processes = nproc) as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, a =2., pool=pool, args=[name, DM, alpha])
        pos, prob, state=sampler.run_mcmc(p0, niteration, progress=True)

# ## read out chains 
    chain = sampler.get_chain()
    lnprobability= sampler.lnprobability[:,:]
    acceptance_fraction = sampler.acceptance_fraction[:]

# ## writing chains

    fout = h5py.File('output/chain_'+name+'.hdf5', 'w')
    fout.create_dataset('chain', data = chain)
    fout.create_dataset('lnprobability', data = lnprobability)
    fout.create_dataset('acceptance_fraction', data = acceptance_fraction)
    fout.close()


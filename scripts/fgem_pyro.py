import os

import jax.numpy as np
import jax.random as random
from IPython.display import set_matplotlib_formats
import jax.numpy as np
from jax import random, vmap
from jax.scipy.special import logsumexp
import matplotlib.pyplot as plt
import numpy as onp
import pandas as pd
import seaborn as sns
from jax.scipy.special import expit, logit
import numpyro
from numpyro.diagnostics import hpdi
import numpyro.distributions as dist
from numpyro import handlers
from numpyro.infer import MCMC, NUTS

plt.style.use('bmh')
if "NUMPYRO_SPHINXBUILD" in os.environ:
    set_matplotlib_formats('svg')

assert numpyro.__version__.startswith('0.2.4')

def fgem_normal(X,BF=None):
    r"simplest FGEM possible"
    #    num_genes = BF.shape[0]
    num_features = X.shape[1]
    Beta0_prior = dist.Normal(-1,5)
    Beta0  = numpyro.sample("Beta0",Beta0_prior)
    Beta_prior = dist.Normal(np.zeros(num_features),np.ones(num_features)*5)
    Beta = numpyro.sample("Beta",Beta_prior)
    pvec = expit(np.sum( Beta * X + Beta0,axis=-1))
    lik = np.sum(np.log(pvec * BF + (1-pvec)))
    numpyro.factor('lik',lik)

import h5py
import os
os.getcwd()
with h5py.File("../data/models_hdf/UCS.h5",'r') as h5f:
    X = (h5f['x'][:,:]).transpose()
    y= h5f['y'][:]

print(X.shape)

rng_key = random.PRNGKey(0)
rng_key, rng_key_ = random.split(rng_key)

num_warmup, num_samples = 1000, 2000

kernel = NUTS(fgem_normal)
mcmc = MCMC(kernel, num_warmup, num_samples, num_chains=1,
                progress_bar=True)
mcmc.run(rng_key, X, y)
ms=mcmc.get_samples()
mcmc.print_summary()

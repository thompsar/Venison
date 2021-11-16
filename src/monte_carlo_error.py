from numpy.random import default_rng
from scipy.stats import f as ftest
import numpy as np
from multiprocessing import Pool, cpu_count
import multiprocessing as mp


def mc_wrapper(model_func,data, origin, weights, iterations, stat_cutoff):
    # GET IT?! LOL
    #the acutal monte carlo function for searching the error space
    
    rng = default_rng()
    n_pop, n_params = origin.shape
    random_hops = rng.standard_normal((iterations,n_pop,n_params))*weights

    landscape_RSS = np.ones(iterations)
    model_landscape = np.zeros((iterations,n_pop,n_params))


    origin_RSS = np.linalg.norm(data - model_func(origin))
    model_DOF = len(data) - np.prod((n_pop,n_params))

    new_hop = origin

    for idx, hop in enumerate(random_hops):

        hop_RSS = np.linalg.norm(data - model_func(new_hop))
        landscape_RSS[idx] = hop_RSS
        model_landscape[idx] = new_hop

        if ftest.cdf(hop_RSS/origin_RSS,model_DOF,model_DOF) > stat_cutoff:
            new_hop = abs(origin+hop) #prevent negative values, particularly in amp
        else:
            new_hop = abs(new_hop+hop) #prevent negative values, particularly in amp
            
    
    landscape_statistics = ftest.cdf(landscape_RSS/origin_RSS,model_DOF,model_DOF)
    return model_landscape, landscape_statistics

def monte_carlo_error(model_func,data, origin, weights, iterations = 40000, stat_cutoff = 0.95):
    
    partial_iterations = iterations//cpu_count()

    pool = mp.Pool()
    results = [pool.apply_async(mc_wrapper, args = (model_func, data, origin, weights, partial_iterations, stat_cutoff,)) for cpu in range(cpu_count())]
    solutions = [p.get() for p in results]
    pool.close()

    landscape_statistics = np.array([solution[1] for solution in solutions]).flatten()

    model_landscape = np.array([solution[0] for solution in solutions]).reshape(len(landscape_statistics),*weights.shape)
    
    
    #convert amplitude to normalized mole fraction. Doesnt materially affect resultant waveform, differences are in machine precision.
    ####this is of course gaussian model specific! Whatch out####
    #model_landscape[:,:,2] = model_landscape[:,:,2]/np.sum(model_landscape[:,:,2],axis=1)[:,np.newaxis]
    
    
    return model_landscape, landscape_statistics

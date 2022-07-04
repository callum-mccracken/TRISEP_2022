# Use a simple Metropolis approach to produce a Markov Chain

# Usage:
# parameters = [
# parameters = [
#    {'name':'a','start':0.4,'step':0.01,'min':0.02,'max':1.2},
#    {'name':'b','start':0.2,'step':0.01,'min':0.04,'max':0.4}
# ]
# mcmc = MarkovChain(parameters, logP)
# chain = mcmc.get_chain(1000)

import numpy as np
from scipy import stats, special
import pandas as pd

import matplotlib.pyplot as plt


class MarkovChain:
    def __init__(self, parameter_list, logP):
        self.names = [par['name'] for par in parameter_list]
        self.start = {}
        self.hypercube = {}
        self.min = {}
        self.max = {}
        for par in parameter_list:
            self.start[par['name']] = par['start']
            self.hypercube[par['name']] = par['step']
            self.min[par['name']] = par['min']
            self.max[par['name']] = par['max']
        self.logP = logP

    def get_chain(self, n_points):

        params = self.start.copy()
        chain = []

        n_accept = 0
        lp_curr = self.logP(params)
        for i in range(n_points):
            new_params = params.copy()
            valid = True
            for par_name in self.names:
                new_params[par_name] += self.hypercube[par_name] * (1. - 2. * stats.uniform.rvs())
                valid = valid and (self.min[par_name] <= new_params[par_name] <= self.max[par_name])
            if valid:
                lp_new = self.logP(new_params)
                if lp_new != -np.inf:
                    del_lp = lp_new - lp_curr
                    if del_lp > -30:
                        if del_lp > 0 or (stats.uniform.rvs() < np.exp(del_lp)):
                            params = new_params
                            lp_curr = lp_new
                            n_accept += 1
            chain.append(params)

        print('Acceptance fraction:', n_accept / n_points)
        return chain

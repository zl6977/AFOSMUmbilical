# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 18:29:54 2020

@author: ZL
"""
import matplotlib.pyplot as plt
import numpy as np
# from numpy import sqrt, pi
from scipy.stats import norm
from UC_package import UC_reliability_zzz


def main(instance_in):
    testTimes=100000
    lengthX = len(instance_in.mu_xi)   #dimension of the vector X
    x_set = np.zeros((testTimes, lengthX), dtype=np.float)
    for i in range(lengthX):
        x_set[:,i] = np.random.normal(instance_in.mu_xi[i], instance_in.sigma_xi[i], testTimes)

    failureTimes=0
    z =[0]*testTimes
    sigma_equ =[0]*testTimes
    for count in range(testTimes):
        sigma_equ[count]=instance_in.sigma_equ_func(x_set[count, :])
        z[count]=instance_in.g_func(x_set[count, :])
        # sigma_equ[count]=instance_in.gama*x_set[count,4] - z[count]
        if z[count] <= 0:
            failureTimes += 1
    reliability = 1-failureTimes/testTimes
    return z, sigma_equ, reliability

if __name__ == "__main__":
    instance = UC_reliability_zzz()
    instance.mu_xi = [29.14,1.87,15.52,1.41,550]
    instance.sigma_xi = [0.043,0.062,0.043,0.047,22.6]
    instance.p = 34.5              #MPa
    instance.T = 178.8e3           #N
    instance.Curvanture = 0.11e-3  #mm^-1
    instance.Et = 2.06e5           #Mpa <- 2.06e11 Pa
    instance.gama = 1.0            #utilization
    instance.n_in, instance.n_out = 102,110
    instance.R_io = 3.15/2
    instance.alpha_helix_in, instance.alpha_helix_out = 36/180*3.14, 30/180*3.14

    print(instance.mu_xi)
    z_ins, sigma_equ_ins, reliability_ins = main(instance)

    print(norm.ppf(reliability_ins))
    print(reliability_ins)

    (mu, sigma) = norm.fit(sigma_equ_ins)
    counts, bins, ignored = plt.hist(sigma_equ_ins, 'auto', density=True)
    plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) *
                    np.exp( - (bins - mu)**2 / (2 * sigma**2) ),
              linewidth=2, color='r')
    plt.show()
    
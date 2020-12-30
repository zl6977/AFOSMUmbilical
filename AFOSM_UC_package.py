# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 2020

@author: lianz
"""
# clear all
# from scipy.misc import derivative
from scipy.optimize import fsolve
from scipy.stats import norm
from numpy import sqrt
# import numpy as np
from UC_package import UC_reliability_zzz

class AFOSM_UC_zzz(UC_reliability_zzz):
    #to calc the partial derivative
    def pgpxi_i_func(self, x_in, i):
        deltaXi=0.001*x_in[i]
        x_nm1 = x_in[:]
        x_nm1[i] -= deltaXi
        x_np1 = x_in[:]
        x_np1[i] += deltaXi
        y_nm1=self.g_func(x_nm1)
        y_np1=self.g_func(x_np1)
        pgpxi_tmp=(y_np1 - y_nm1)/deltaXi/2
        return pgpxi_tmp

    def pgpxi_func(self, x_in):
        length_x_in = len(x_in)
        pgpxi_tmp=[0]*length_x_in
        for i in range(length_x_in):
            pgpxi_tmp[i]=self.pgpxi_i_func(x_in, i)
        return pgpxi_tmp    
 
    #2.sensitive coefficient
    def alphai_func(self, pgpxi_in, sigma_xi):
        length_pgpxi_in = len(pgpxi_in)
        alphai_tmp = [0]*length_pgpxi_in 
        for i in range(length_pgpxi_in):
            numerator = sigma_xi[i]*pgpxi_in[i]
            denominator = 0
            for j in range(length_pgpxi_in):
                tmp1  = (sigma_xi[j]*pgpxi_in[j])**2
                denominator +=  tmp1
            denominator = sqrt(denominator)
            alphai_tmp[i] = numerator/denominator
        return alphai_tmp
    
    #3.define xi_star & g_func_containbeta
    def xi_star_func(self, beta_in, mu_xi, sigma_xi, alphai):
        length_X= len(mu_xi)
        xi_star_tmp = [0]*length_X
        for i in range(length_X):
            xi_star_tmp[i] = mu_xi[i]-beta_in*alphai[i]*sigma_xi[i]
        return xi_star_tmp
    
    def g_func_contain_beta(self, beta_in2, mu_xi, sigma_xi, alphai):
        xi_star_tmp = self.xi_star_func(beta_in2, mu_xi, sigma_xi, alphai)
        g_func_beta_tmp = self.g_func(xi_star_tmp)
        return g_func_beta_tmp
#--------------end of class definition
  
       
def main(instance_in):
    #------------------------start of the iteration------------------------
    #1.choose the initial value
    xi_star = instance_in.mu_xi[:]
    
    #2.sensitive coefficient
    # alphai_ins = alphai_func(*xi_star_ins)
    # sigmaZ_ins = sigmaZ_func(*xi_star_ins)
    pgpxi = instance_in.pgpxi_func(xi_star)
    alphai = instance_in.alphai_func(pgpxi, instance_in.sigma_xi)
        
    #3.define xi_star & g_func_containbeta
        
    #4.update beta
    argstuple = (instance_in.mu_xi, instance_in.sigma_xi, alphai)
    beta_cur = fsolve(instance_in.g_func_contain_beta,0,args=argstuple)[0]
    
    #5.update xi_star
    xi_star = instance_in.xi_star_func(beta_cur,instance_in.mu_xi, instance_in.sigma_xi, alphai)
    pgpxi = instance_in.pgpxi_func(xi_star)
    # print("0 , xi_star", xi_star)
    
    count = 1
    iter_flag=True
    while(iter_flag):
        #2.sensitive coefficient
        alphai=instance_in.alphai_func(pgpxi, instance_in.sigma_xi)
        
        #3.define xi_star & g_func_contain_beta
            
        #4.update beta
        beta_prev = beta_cur
        argstuple=(instance_in.mu_xi, instance_in.sigma_xi, alphai)
        beta_cur = fsolve(instance_in.g_func_contain_beta,0,args=argstuple)[0]   #discard other roots
        #print(count, "beta_cur", beta_cur)
        
        if abs(beta_cur-beta_prev) <= 1e-6:
            iter_flag=False
        
        #5.update xi_star
        xi_star = instance_in.xi_star_func(beta_cur, instance_in.mu_xi, instance_in.sigma_xi, alphai)
        pgpxi = instance_in.pgpxi_func(xi_star)
        # print(count, "xi_star", xi_star)
        count += 1
    #------------------------end of the iteration------------------------ 
    return beta_cur,xi_star   
    
if __name__ == "__main__":
    instance=AFOSM_UC_zzz()
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
    beta_cur_ins,xi_star_ins = main(instance)
    print(beta_cur_ins)
    print(norm.cdf(beta_cur_ins))
    print(xi_star_ins)
# -*- coding: utf-8 -*-
"""
Created on Nov 4 2020

@author: lianz
"""
from numpy import cos, sqrt, pi
# import numpy as np

class UC_reliability_zzz:
    #define the question, variables, fomula
    # lengthX = 5   #dimension of the vector X
    # mu_xi = [29.14,1.87,15.52,1.41,550]
    # sigma_xi = [0.043,0.062,0.043,0.047,22.6]
    
    # p = 34.5              #MPa
    # T = 178.8e3           #N
    # Curvanture = 0.11e-3  #mm^-1
    # Et = 2.06e5           #Mpa <- 2.06e11 pa
    # gama = 1.0            #utilization
            
    # n_in,n_out = 102,110
    # R_io = 3.15/2
    # alpha_helix_in, alpha_helix_out = 36/180*3.14, 30/180*3.14
    
    mu_xi = []
    sigma_xi = []
    
    p = 0               #MPa
    T =0                #N
    Curvanture = 0      #mm^-1
    Et = 0              #Mpa <- 2.06e11 pa
    gama = 0            #utilization
            
    n_in,n_out = 0,0
    R_io = 0
    alpha_helix_in, alpha_helix_out = 0,0    
   
    #to calc the linearized performance function
    def sigma_equ_func(self,x_in):
        D = x_in[0]
        d = x_in[0] - 2*x_in[1]
        Rt = D/2
        At = pi/4*(D**2-d**2)
        
        R_io = self.R_io
        Et = self.Et
        T = self.T
        n_in = self.n_in
        n_out = self.n_out
        R_io = self.R_io
        alpha_helix_in = self.alpha_helix_in
        alpha_helix_out = self.alpha_helix_out
        p = self.p        
        T = self.T
        Curvanture = self.Curvanture
        Et = self.Et
        # gama = self.gama
                
        A_io = pi * R_io**2
        
        Kt = Et * At
        K = Kt*4 + Et*pi/4*((x_in[2])**2-(x_in[2]-2*x_in[3])**2)*5 + n_in*Et*A_io*(cos(alpha_helix_in))**3 + n_out*Et*A_io*(cos(alpha_helix_out))**3
        # K=5*Kt
        
        sigma_t = T/At*Kt/K
        sigma_M = Et*Rt*Curvanture
        sigma_e = p*pi*(d**2)/4/At
        
        sigma_h = p*(D**2+d**2)/(D**2-d**2)
        sigma_r = -p
        sigma_a = sigma_t + sigma_M + sigma_e
        sigma_equ = sqrt(((sigma_a-sigma_r)**2+(sigma_a-sigma_h)**2+(sigma_h-sigma_r)**2)/2)
        return sigma_equ
    
    def g_func(self,x_in):  
        gama = self.gama
        sigma_equ = self.sigma_equ_func(x_in)
        return gama*x_in[4]-sigma_equ

if __name__ == "__main__":
    instance=UC_reliability_zzz()
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

    print(instance.sigma_equ_func(instance.mu_xi))    
    print(instance.g_func(instance.mu_xi))
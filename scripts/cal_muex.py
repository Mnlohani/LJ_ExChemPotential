#!/usr/bin/python
import numpy as np
T = 2.74
N = 10000
rho = np.linspace(.1, 1, 10)
#results from running ti_sim.py for density ranges 
delta_f = np.array([18804.8539938, 11531.1209968, 9089.94540078, 8340.27469046, 8585.9484523, 9668.35653151, 11553.8273517, 14322.7198397, 18061.2125604])
#results from running presure_endpoints.py
pressure_end = np.array([0.2695379589420549, 0.55181070258185494, 0.89457903033746566, 1.3598439088367786, 2.0739050349116104, 3.232647245548832, 5.130325632773526, 8.2018018169998879, 12.973711650794842, 20.269730392477385])

#Calculation using project method
delta_mu_intgr = delta_f / N
delta_mu_id = np.zeros(9)
delta_P_by_rho = np.zeros(9)
for i in range(9):
    delta_mu_id[i] = T * np.log(rho[i+1] / rho[i])
    delta_P_by_rho[i] =  pressure_end[i+1] / rho[i+1] - pressure_end[i] / rho[i]
    
delta_mu_ex = delta_mu_intgr + delta_P_by_rho - delta_mu_id
print(delta_mu_ex)

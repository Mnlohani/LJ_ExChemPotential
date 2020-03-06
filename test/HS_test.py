import numpy as np
from scipy import integrate

T = 2.74
N = 10000
rho_test= np.linspace(.1, 1., 10)
eta_test = np.pi * rho_test/ 6
n_test = np.pi * rho_test/ 6

#formula for (Z) 3.9.17
Z_test =  (1 + n_test + n_test**2 - n_test**3) / (1 - n_test)**3
P_test = T * rho_test * Z_test

#Calculation using project method
delta_mu_intgr = np.zeros(9)
delta_mu_id_test = np.zeros(9)
delta_P_by_rho = np.zeros(9)
def integration_term(x):
    return T *(1 + x + x**2 - x**3) / (x * (1-x)**3)
for i in range(9):
    delta_mu_intgr[i] = integrate.quad(integration_term, eta_test[i], eta_test[i +1])[0]
    delta_mu_id_test[i] = T * np.log(rho_test[i+1] / rho_test[i])
    delta_P_by_rho[i] =  P_test[i+1] / rho_test[i+1] -P_test[i] / rho_test[i]

delta_mu_ex_test = delta_mu_intgr + delta_P_by_rho - delta_mu_id_test
print(delta_mu_ex_test)

# Calculation by formula 3.9.18
mu_ex_test_delta = np.zeros(9)
F_ex_test =  N * T * n_test * (4 - 3 * n_test) / (1 - n_test)**2
mu_ex_test = F_ex_test / N + P_test / rho_test  - T
for i in range(len(rho_test) - 1):
    mu_ex_test_delta[i] =  mu_ex_test[i+1] - mu_ex_test[i]
print('Ex chemical potential change \n', mu_ex_test_delta)
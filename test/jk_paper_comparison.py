import numpy as np
from scipy import integrate

#Table 10 : parameters are x1 as x[1],..so on
x = np.array([0.0, 0.8623085097507421,2.976218765822098,-8.402230115796038,0.1054136629203555,-0.8564583828174598,
1.582759470107601,0.7639421948305453 ,1.753173414312048,2.798291772190376e3,-4.8394220260857657e-2,0.9963265197721935,
-3.698000291272493e1,2.084012299434647e1,8.305402124717285e1,-9.574799715203068e2,-1.477746229234994e2,
6.398607852471505e1,1.603993673294834e1,6.805916615864377e1,-2.791293578795945e3,-6.245128304568454,-8.116836104958410e3,
1.488735559561229e1,-1.059346754655084e4,-1.131607632802822e2,-8.867771540418822e3,-3.986982844450543e1,-4.689270299917261e3, 2.593535277438717e2, -2.694523589434903e3,-7.218487631550215e2,1.721802063863269e2])

def SetCoficients(T):
    a = np.zeros(9)
    b = np.zeros(7)
    # Table 5 : a cofficients:
    a[1] = x[1] * T + x[2] * T**0.5 + x[3] + x[4] / T + x[5] / T**2
    a[2] = x[6] * T + x[7] + x[8] / T  + x[9] /  T**2
    a[3] = x[10] * T + x[11] + x[12] / T 
    a[4] = x[13]
    a[5] = x[14] / T + x[15] / T**2
    a[6] = x[16] / T
    a[7] = x[17] / T + x[18] / T**2
    a[8] = x[19] / T**2
    # Table 6 for b cofficients:
    b[1] = x[20] / T**2 + x[21] / T**3
    b[2] = x[22] / T**2 + x[23] / T **4
    b[3] = x[24] / T**2 + x[25] / T**3
    b[4] = x[26] / T**2 + x[27] / T**4
    b[5] = x[28] / T**2 + x[29] / T**3
    b[6] = x[30] / T**2 + x[31] / T**3 + x[32] / T**4
    return a, b

def Pressure_jk(rho_1, temp):
    #Pressure Calculation for testing purposes:
    _a, _b = SetCoficients(temp)
    f = np.exp(-3*rho_1**2)
    second_term = np.sum(np.array([ _a[i] * rho_1**(i+1) for i in range(9)]))
    third_term =  f * np.sum(np.array([ _b[i] * rho_1**(2*i+1) for i in range(7)]))
    pressure = rho_1 * temp + second_term + third_term
    return pressure 

def intergation_jk_(rho_1, temp, N):
    _a, _b = SetCoficients(temp)
    f = np.exp(-3*rho_1**2)
    second_term = np.sum(np.array([ _a[i] * rho_1**(i+1) for i in range(9)]))
    third_term =  f * np.sum(np.array([ _b[i] * rho_1**(2*i+1) for i in range(7)]))
    pressure = rho_1 * temp + second_term + third_term
    return - N * pressure / rho_1**2

def Cal_Free_energy(temp, rho_list, N):
    delta_F_ = np.zeros(9)
    for i in range(9):
        delta_F_[i] = integrate.quad(intergation_jk_, rho_list[i+1], rho_list[i], args=(temp,N))[0]
    return delta_F_

#test pressure calculation for T = 1.3
rho_13 = [.1, .2, .4, .5, .6, .7, .8, .9, .95]
pressure_jk_test = [ Pressure_jk(_rho, 1.3) for _rho in rho_13]
print(pressure_jk_test)

# Free Energy results for T = 2.74 and N =10000
rho_test= np.linspace(.1, 1., 10)
delta_F_jk = Cal_Free_energy(2.74, rho_test, 10000)
print(delta_F_jk)
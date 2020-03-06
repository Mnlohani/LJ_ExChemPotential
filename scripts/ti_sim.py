#!/usr/bin/python
import sys, os, subprocess, glob, h5py
import numpy as np
from scipy import integrate
from datetime import datetime

#--------Block 1 ----------- 
#arguments and variables
nparticles = int(sys.argv[1])
init_density = float(sys.argv[2])
final_density = float(sys.argv[3])
init_volume = nparticles / init_density
final_volume = nparticles / final_density
pressure_steps = []
volume_steps = []
var_pressure = []
rho_pressure_dict = {}

#functions
def simulation_time_cal(density):
    '''returns simulation time as per density 
    '''
    factor = 1 + 0.2 * density
    time_sim = 1000 / factor
    return time_sim

def Is_record_exist(volume_value):
    '''checks if 'step volume of interest' exists in saved calculated .npy files or not '''
    record_flag = 0 #for record existence check
    step_V = round(volume_value, 4) 
    npfiles= glob.glob("*.npy")
    if npfiles == []:
        record_flag = 0
    else:
        for j, npfile in enumerate(npfiles):
            dict_ndarray = np.load(npfiles[j])
            keys = dict_ndarray.item()
            if step_V in keys:
                print('record found')
                step_pressure = dict_ndarray.item().get(step_V)
                record_flag = 1
                break
    if record_flag == 0:
        return False , 0
    else:
        return True , step_pressure

def Calculate_Pressure(volume,num_particles):
    density = num_particles / volume
    time = simulation_time_cal(density)
    print(time)
    record_bool, pre_cal_pressure = Is_record_exist(volume)
    if record_bool:
        volume_steps.append(volume)
        pressure_steps.append(pre_cal_pressure)
        return pre_cal_pressure
    else:
        #section 1: 
        #running lennard_jones_equilibriation
        subprocess.call(["halmd", "t_lennard_jones_equilibration.lua", "-v", "--timestep", "0.005", "--time", str(time), "--density", str(density), "--particles", str(num_particles), "--cutoff","4", "--temperature", "2.74", "--sampling", "state-vars=100"])
        #Read Lj equillibration generated .h5 file
        os.system('ls -t |grep lennard.*h5|head -1 > tmp')
        f = open('tmp', 'r')
        line = f.readlines()
        argument = line[0]
        s = argument[: -1] # removing newline char eg. abc\n
        f.close()
        path='./'+ s
        h5data = h5py.File(path, 'a')
        pressure =  h5data['observables/pressure/value'][:]
        len_sim = len(pressure)/2 #second half of simulation
        avg_pressure = np.mean(pressure[len_sim:])
        var_p = np.var(pressure[len_sim:])
        pressure_steps.append(avg_pressure)
        var_pressure.append(var_p)
        volume_steps.append(volume)
        h5data.close()
        #removing .h5 files , keeping logs files
        subprocess.call(["rm", str(s)])
        return avg_pressure

#--------Block 2 ----------- 
# free Energy calculation

free_energy ,err = integrate.quad(Calculate_Pressure, init_volume, final_volume, args=(nparticles,),epsabs=1.0 ,epsrel=0.01)

#Commented code for testing of existing record
'''
step_volume = [453.33333333333337, 300.07878654, 301.52524267, 297.09786756]
avg_p = []
for i in range(len(step_volume)):
    avg_p.append(Calculate_Pressure(step_volume[i], nparticles))

free_energy = integrate.simps(avg_p, step_volume)
'''

#--------Block 3 ----------- 
# save dictionary of type 'Volume:Pressure'  
##volumes are rounded to 4 decimal places 
for i in range(len(volume_steps)):
    rho_pressure_dict[round(volume_steps[i],4)] = pressure_steps[i]
timestamp = datetime.now().strftime('%Y%m%d%H%M%S')
dict_name = 'PVdict' + str(timestamp)
np.save('%s.npy' % dict_name, rho_pressure_dict )

#--------Block 3 ----------- 
# print values
print('step pressure values :')
print(pressure_steps)
print(' step volume values:')
print(volume_steps)
print('free energy:')
print(free_energy)
print('error :')
print(err)
print('variance')
print(var_pressure)
#print('dictionary :')
#dict2 = np.load('%s.npy' % dict_name)
#print(dict2.item())


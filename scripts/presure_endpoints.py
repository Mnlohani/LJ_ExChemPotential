#!/usr/bin/python
import sys, os, subprocess
import numpy as np
import h5py
from scipy import integrate
import matplotlib.pyplot as plt

#arguments
num_particles = int(sys.argv[1])

#local variables
density_steps = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
pressure_steps = []

#functions
def simulation_time_cal(density):
    '''returns simulation time as per density '''
    factor = 1 + 0.2 * density
    time_sim = 1500 / factor
    return time_sim

def Calculate_Pressure(density):   
    #section 1:
    print(density)
    time = simulation_time_cal(density)
    #section 1: 
    #running lennard_jones_equilibriation
    #Read Lj equillibration generated .h5 file
    subprocess.call(["halmd", "t_lennard_jones_equilibration.lua", "-v", "--timestep", "0.005", "--time", str(time), "--density", str(density), "--particles", str(num_particles),"--cutoff","4" , "--temperature", "2.74", "--sampling", "state-vars=100"])
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
    pressure_steps.append(avg_pressure)
    h5data.close()
    #removing .h5 files , keeping logs files
    subprocess.call(["rm", str(s)])


for i in range(len(density_steps)):
    Calculate_Pressure(density_steps[i])

print(pressure_steps)

 

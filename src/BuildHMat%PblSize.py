#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import os , time, sys
import matplotlib.animation as animation
from math import exp, sqrt, ceil, log10
import json

##################### Plot vs pbl size #####################

# Parameters
Result_path = "./output/BuildHMat%PblSize/1024_2048_4096_8192/"
Result_file = "bench_hmatrix_build.json"
list_N = np.array([1024, 2048, 4096, 8192])
list_implementation = ["Dense", "Classic", "TaskBased"]
list_threads = [1]
Is_log_scale = 1

# Variables
mat_real_time = np.zeros((len(list_N), len(list_implementation), len(list_threads)))
mat_cpu_time = np.zeros((len(list_N), len(list_implementation), len(list_threads)))

# read json
f = open(Result_path+Result_file)
data = json.load(f)
f.close()

for BM in data['benchmarks']:
    for id_N in range(len(list_N)):
        for id_impl in range(len(list_implementation)):
            for id_threads in range(len(list_threads)):
                if (BM['name'] == "FT_Generator/BM_"+str(list_implementation[id_impl])+"/N:"+str(list_N[id_N])+"/threads:"+str(list_threads[id_threads])+"_median"):
                    mat_real_time[id_N, id_impl, id_threads] = BM['real_time']
                    mat_cpu_time[id_N, id_impl, id_threads] = BM['cpu_time']
                    
# Plots               
id_threads = 0   

# Plot real time
plth = plt.subplot(211)
plth.set_title("Building median real time as a function of pbl size")
plth.set_ylabel("Real time (us)")
for id_impl in range(len(list_implementation)):
    plth.plot(list_N, mat_real_time[:, id_impl, id_threads], '+-', label=list_implementation[id_impl]+ ", thread(s): "+str(list_threads[id_threads]))
if Is_log_scale:
    plt.yscale('log')
    plt.xscale('log')
    plth.plot(list_N, list_N*np.log10(list_N), '+-', label='NlogN')
    plth.plot(list_N, list_N*list_N, '+-', label='N^2')
plth.legend()

# Plot cpu time
pltu = plt.subplot(212)
pltu.set_title("Building cpu time as a function of pbl size")
plth.set_xlabel("Problem size N")
pltu.set_ylabel("CPU time (us)")
for id_impl in range(len(list_implementation)):
    pltu.plot(list_N, mat_cpu_time[:, id_impl, id_threads], '+-', label=list_implementation[id_impl]+ ", thread(s): "+str(list_threads[id_threads]))
if Is_log_scale:
    plt.yscale('log')
    plt.xscale('log')
    pltu.plot(list_N, list_N*np.log10(list_N), '+-', label='NlogN')
    pltu.plot(list_N, list_N*list_N, '+-', label='N^2')
pltu.legend()

# Save and show figure
plt.savefig(Result_path+'./Time_vs_PblSize.png', format='png',bbox_inches = "tight")    
plt.show()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import os , time, sys
import matplotlib.animation as animation
from math import exp, sqrt, ceil, log10
import json

##################### Plot vs Threads number #####################

# Parameters
Result_path = "./output/BuildHMat%ThreadsNumber/1_2_4_8_16/"
Result_file = "bench_hmatrix_build.json"
list_N = np.array([8192])
list_implementation = [ "Classic", "TaskBased"] #["Dense", "Classic", "TaskBased"]
list_threads = np.array([1, 2, 4, 8, 16])
Is_log_scale = False

# Variables
median_real_time = np.zeros((len(list_N), len(list_implementation), len(list_threads)))
median_cpu_time = np.zeros((len(list_N), len(list_implementation), len(list_threads)))
stddev_real_time = np.zeros((len(list_N), len(list_implementation), len(list_threads)))
stddev_cpu_time = np.zeros((len(list_N), len(list_implementation), len(list_threads)))

# read json
f = open(Result_path+Result_file)
data = json.load(f)
f.close()

for BM in data['benchmarks']:
    for id_N in range(len(list_N)):
        for id_impl in range(len(list_implementation)):
            for id_threads in range(len(list_threads)):
                if (BM['name'] == "FT_Generator/BM_"+str(list_implementation[id_impl])+"/N:"+str(list_N[id_N])+"/threads:"+str(list_threads[id_threads])+"_median"):
                    median_real_time[id_N, id_impl, id_threads] = BM['real_time']
                    median_cpu_time[id_N, id_impl, id_threads] = BM['cpu_time']
                if (BM['name'] == "FT_Generator/BM_"+str(list_implementation[id_impl])+"/N:"+str(list_N[id_N])+"/threads:"+str(list_threads[id_threads])+"_cv"):
                    stddev_real_time[id_N, id_impl, id_threads] = BM['real_time']   
                    stddev_cpu_time[id_N, id_impl, id_threads] = BM['cpu_time']                    
# Plots               
id_N = 0   

# Plot real time
plth = plt.subplot(211)
plth.set_title("Building median real and cpu times as a function of threads number")
plth.set_ylabel("Real time (us)")
for id_impl in range(len(list_implementation)):
   plth.errorbar(list_threads, median_real_time[id_N, id_impl, :], yerr=stddev_real_time[id_N, id_impl, :], fmt='-o', capsize=3, label=list_implementation[id_impl]+ ", N: "+str(list_N[id_N]))
   
if Is_log_scale:
    plt.yscale('log')
    plt.xscale('log')
    # plth.plot(list_threads, list_threads*np.log10(list_threads), '+-', label='NlogN')
    # plth.plot(list_threads, list_threads*list_threads, '+-', label='N^2')
plth.legend()

# Plot cpu time
pltu = plt.subplot(212)
pltu.set_xlabel("Number of threads")
pltu.set_ylabel("CPU time (us)")
for id_impl in range(len(list_implementation)):
    pltu.errorbar(list_threads, median_cpu_time[id_N, id_impl, :], yerr=stddev_cpu_time[id_N, id_impl, :], fmt='-o', capsize=3, label=list_implementation[id_impl]+ ", N: "+str(list_N[id_N]))
if Is_log_scale:
    plt.yscale('log')
    plt.xscale('log')
    # pltu.plot(list_threads, list_threads*np.log10(list_threads), '+-', label='NlogN')
    # pltu.plot(list_threads, list_threads*list_threads, '+-', label='N^2')
pltu.legend()

# Save and show figure
plt.savefig(Result_path+'./Time_vs_ThreadsNumber.png', format='png',bbox_inches = "tight")    
plt.show()

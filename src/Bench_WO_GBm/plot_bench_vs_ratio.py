#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import os , time, sys
import matplotlib.animation as animation
from math import exp, sqrt, ceil, log10
import pandas

def geometric_progression_list(a, r, length):
    return [a * r ** (n - 1) for n in range(1, length + 1)]

######################## Parameters #########################################

List_pbl_size          = geometric_progression_list(1 << 14, 2, 5)
List_number_of_threads = geometric_progression_list(1, 2, 5)
List_implementation    = ["Classic", "TaskBased"]
List_epsilon           = [1e-10, 1e-8, 1e-6]
is_log_scale           = False
# plot_result_file       = "/bench_hmatrix_build_vs_ratio.csv"; bench_result_path = "../../output/Build_hmatrix_vs_ratio/" # Build_hmatrix_vs_ratio
plot_result_file       = "/bench_hmatrix_matrix_product_vs_ratio.csv"; bench_result_path = "../../output/Hmatrix_matrix_product_vs_ratio/" # Hmatrix_matrix_product_vs_ratio

#############################################################################

# Update bench_result_path with List_pbl_size
for id_N in range(len(List_pbl_size)):
    bench_result_path += str(List_pbl_size[id_N]) + "_"
bench_result_path = bench_result_path.rstrip('_') # Remove last element of bench_result_path

# Create bench_result_path if it doesn't exist
os.makedirs(bench_result_path, exist_ok=True)

# Import data
df = pandas.read_csv(bench_result_path+plot_result_file, sep=', ', engine='python') 
data = 'time (s) | mean time (s) | standard_deviation'

# Variables
Mat_mean_time = np.zeros((len(List_epsilon), len(List_pbl_size), len(List_implementation)))
Mat_stddev = np.zeros_like(Mat_mean_time)

for id_BM in range(len(df[:])): 
    id_epsilon = List_epsilon.index(df['epsilon'][id_BM])
    id_N       = List_pbl_size.index(df['dim_pbl'][id_BM])
    id_impl    = List_implementation.index(df['algo_type'][id_BM])
    
    if df['id_rep'][id_BM] == 'mean':
        Mat_mean_time[id_epsilon, id_N, id_impl] = df[data][id_BM]
    if df['id_rep'][id_BM] == 'stddev':
        Mat_stddev[id_epsilon, id_N, id_impl] = df[data][id_BM]
             
# Plots     
List_markers = ['o', 's', 'D', 'v', '^', 'h', 'd'] # for epsilon
List_colors  = ['b', 'r', 'g', 'c', 'm', 'y', 'k'] # for implementation
          
# Plot real time
pltA = plt.subplot(211)
pltA.set_title("Mean real time of hmatrix building as a function of ratio: " + str(List_pbl_size) + " / " + str(List_number_of_threads))
pltA.set_ylabel("Real time (s)")

for id_impl in range(len(List_implementation)):
    for id_epsilon in range(len(List_epsilon)):        
        pltA.errorbar(List_pbl_size, Mat_mean_time[id_epsilon, :, id_impl], yerr=Mat_stddev[id_epsilon, :, id_impl], fmt=List_markers[id_epsilon]+'-'+List_colors[id_impl], markerfacecolor='none', capsize=3, label=List_implementation[id_impl]+ ", epsilon: "+str(List_epsilon[id_epsilon]))
        
if is_log_scale:
    plt.yscale('log')
    plt.xscale('log')
    Array_N = np.array(List_pbl_size)
    pltA.plot(List_pbl_size, Array_N*np.log10(Array_N) / (Array_N[0]*np.log10(Array_N)[0]) * Mat_mean_time[0, 0, 1], '--', label='NlogN')
    pltA.plot(List_pbl_size, Array_N*Array_N / (Array_N[0]*Array_N[0]) * Mat_mean_time[0, 0, 1], '--', label='N^2')
    
pltA.legend()

# Save and show figure
mng = plt.get_current_fig_manager()
mng.resize(*mng.window.maxsize())
figure = plt.gcf()
plt.show()
figure.savefig(bench_result_path+'/mean_time_vs_ratio.png', format='png',bbox_inches = "tight") 

# cout
print("Plotting : " + bench_result_path + plot_result_file)
print("List_pbl_size : " + str(List_pbl_size))
print("List_implementation : " + str(List_implementation))
print("List_epsilon : " + str(List_epsilon))
print("is_log_scale : " + str(is_log_scale))

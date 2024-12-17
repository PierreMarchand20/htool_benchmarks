#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import os , time, sys
import matplotlib.animation as animation
from math import exp, sqrt, ceil, log10
import pandas

def geometric_progression(a, r, length):
    return [a * r ** (n - 1) for n in range(1, length + 1)]

######################## Parameters #########################################
List_number_of_threads = geometric_progression(1, 2, 5) # ALL degree of freedom in the csv file. Used for the path.
List_implementation    = ["Classic", "TaskBased"]  # Todo faire en sorte que ça marche meme si on ne prend pas certaines implementations
List_epsilon           = [1e-10, 1e-8, 1e-6, 1e-4] # Todo faire en sorte que ça marche meme si on ne prend pas certains epsilon
List_bench_type        = ["Build", "Matrix_product"] #, "Factorization"]
is_log_scale           = False

#cout
print("=============================Parameters================================")
print("List_number_of_threads : " + str(List_number_of_threads))
print("List_implementation : " + str(List_implementation))
print("List_epsilon : "        + str(List_epsilon))
print("List_bench_type : "     + str(List_bench_type))
print("is_log_scale : "        + str(is_log_scale))
print("====================================================================")
############################################################################


# plot_result_file       = "/bench_hmatrix_build_vs_thread.csv"; bench_result_path = "../output/Build_hmatrix_vs_thread/" # Build_hmatrix_vs_thread
# plot_result_file       = "/bench_hmatrix_matrix_product_vs_thread.csv"; bench_result_path = "../output/Matrix_product_hmatrix_vs_thread/" # Hmatrix_matrix_product_vs_thread
#############################################################################
for id_BM_type in range(len(List_bench_type)):
    plot_result_file    = "/bench_hmatrix_"+(List_bench_type[id_BM_type]).lower()+"_vs_thread.csv"; bench_result_path = "../output/"+(List_bench_type[id_BM_type])+"_hmatrix_vs_thread/" 
    
    # Update bench_result_path with List_number_of_threads
    for id_N in range(len(List_number_of_threads)):
        bench_result_path += str(List_number_of_threads[id_N]) + "_"
    bench_result_path = bench_result_path.rstrip('_') # Remove last element of bench_result_path
    print("Plotting : " + bench_result_path + plot_result_file)

    # Create bench_result_path if it doesn't exist
    os.makedirs(bench_result_path, exist_ok=True)

    # Import data
    df = pandas.read_csv(bench_result_path+plot_result_file, sep=', ', engine='python') 
    data = 'time (s) | mean time (s) | standard_deviation'

    # Variables
    Mat_mean_time = np.zeros((len(List_epsilon), len(List_number_of_threads), len(List_implementation)))
    Mat_stddev = np.zeros_like(Mat_mean_time)
    dim_pbl = df['dim_pbl'][0] # expected to be the same for all

    for id_BM in range(len(df[:])): 
        id_epsilon = List_epsilon.index(df['epsilon'][id_BM])
        id_N       = List_number_of_threads.index(df['number_of_threads'][id_BM])
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
    pltA.set_title("Mean real time on problem size of : " + str(dim_pbl) + " as a function of number of threads: "+str(List_number_of_threads))
    pltA.set_ylabel("Real time (s)")

    for id_impl in range(len(List_implementation)):
        for id_epsilon in range(len(List_epsilon)):        
            pltA.errorbar(List_number_of_threads, Mat_mean_time[id_epsilon, :, id_impl], yerr=Mat_stddev[id_epsilon, :, id_impl], fmt=List_markers[id_epsilon]+'-'+List_colors[id_impl], markerfacecolor='none', capsize=3, label=List_implementation[id_impl]+ ", epsilon= "+str(List_epsilon[id_epsilon]))
            
    if is_log_scale:
        plt.yscale('log')
        plt.xscale('log')
        Array_N = np.array(List_number_of_threads)
        pltA.plot(List_number_of_threads, Array_N*np.log10(Array_N) / (Array_N[0]*np.log10(Array_N)[0]) * Mat_mean_time[0, 0, 1], '--', label='NlogN')
        pltA.plot(List_number_of_threads, Array_N*Array_N / (Array_N[0]*Array_N[0]) * Mat_mean_time[0, 0, 1], '--', label='N^2')
        
    pltA.legend()

    # Save and show figure
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    figure = plt.gcf()
    plt.show()
    figure.savefig(bench_result_path+'/mean_time_vs_thread.png', format='png',bbox_inches = "tight") 

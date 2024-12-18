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
List_N              = geometric_progression(1 << 15, 2, 5) # ALL degree of freedom in the csv file. Used for the path.
List_plot_N = List_N                            # Todo faire en sorte que ça marche meme si on ne prend pas certains N 
List_implementation = ["Classic", "TaskBased"]  # Todo faire en sorte que ça marche meme si on ne prend pas certaines implementations
List_epsilon        = [1e-10, 1e-8, 1e-6, 1e-4] # Todo faire en sorte que ça marche meme si on ne prend pas certains epsilon
List_bench_type     = ["Build", "Matrix_product"] #, "Factorization"]
is_log_scale        = False

#cout
print("=============================Parameters================================")
print("List_N : "              + str(List_N))
print("List_plot_N : "         + str(List_plot_N))
print("List_implementation : " + str(List_implementation))
print("List_epsilon : "        + str(List_epsilon))
print("List_bench_type : "     + str(List_bench_type))
print("is_log_scale : "        + str(is_log_scale))
print("====================================================================")
############################################################################

for id_BM_type in range(len(List_bench_type)):
    plot_result_file    = "/bench_hmatrix_"+(List_bench_type[id_BM_type]).lower()+"_vs_pbl_size.csv"; bench_result_path = "../output/"+(List_bench_type[id_BM_type])+"_hmatrix_vs_pbl_size/" 
    
    # Update bench_result_path with List_N
    for id_N in range(len(List_N)):
        bench_result_path += str(List_N[id_N]) + "_"
    bench_result_path = bench_result_path.rstrip('_') # Remove last element of bench_result_path
    print("Plotting : " + bench_result_path + plot_result_file)

    # Create bench_result_path if it doesn't exist
    os.makedirs(bench_result_path, exist_ok=True)

    # Import data
    df = pandas.read_csv(bench_result_path+plot_result_file, sep=', ', engine='python') 
    data = 'time (s) | mean time (s) | standard_deviation'
    # print("Importing data: success!")
    
    # Variables
    Mat_mean_time = np.zeros((len(List_epsilon), len(List_N), len(List_implementation)))
    Mat_stddev = np.zeros_like(Mat_mean_time)

    
    for id_BM in range(len(df[:])): # Todo : faire en sorte que ça marche meme si on ne prend pas certains epsilon, N ou implementation
        id_epsilon = List_epsilon.index(df['epsilon'][id_BM])
        id_N       = List_N.index(df['dim_pbl'][id_BM])
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
    pltA.set_title("Mean real time as a function of dimension problem: "+str(List_N))
    pltA.set_ylabel("Real time (s)")

    for id_impl in range(len(List_implementation)):
        for id_epsilon in range(len(List_epsilon)):        
            pltA.errorbar(List_N, Mat_mean_time[id_epsilon, :, id_impl], 
                        yerr=Mat_stddev[id_epsilon, :, id_impl], 
                        fmt=List_markers[id_epsilon]+'-'+List_colors[id_impl], 
                        markerfacecolor='none', capsize=3, 
                        label=List_implementation[id_impl]+ ", epsilon= "+str(List_epsilon[id_epsilon]))
            
    if is_log_scale:
        plt.yscale('log')
        plt.xscale('log')
        Array_N = np.array(List_N)
        pltA.plot(List_N, Array_N*np.log2(Array_N) / (Array_N[-1]*np.log2(Array_N)[-1]) * Mat_mean_time[0, -1, 1], '--', label='NlogN')
        pltA.plot(List_N, Array_N*Array_N / (Array_N[-1]*Array_N[-1]) * Mat_mean_time[0, -1, 1], '--', label='N^2')
        
    pltA.legend()

    #############################
    pltB = plt.subplot(212)
    pltB.set_title("Rescaling by N*log^a(N)")
    pltB.set_ylabel("Normalized (T_build / f(N))")
    pltB.set_xlabel("Degree of freedom (N)")

    List_fN = [np.array(np.array(List_N) * np.log2(np.array(List_N))), 
            np.array(List_N) * (np.log2(np.array(List_N)))**5, 
            (np.array(List_N))**(3/2),
            np.array(List_N) * (np.log2(np.array(List_N)))**(5/2)]
    List_fN_name = ["NlogN", "Nlog^5(N)", "N^(3/2)", "Nlog^(5/2)(N)"]
    List_linestyle = [':', '-', '--', '-.']

    for id_impl in range(len(List_implementation)):
        for id_epsilon in range(len(List_epsilon)):  
            for id_fN in range(len(List_fN)):      
                pltB.plot(List_N, np.array(Mat_mean_time[id_epsilon, :, id_impl]/List_fN [id_fN])/Mat_mean_time[id_epsilon, 0, id_impl]*List_fN[id_fN][0], 
                        List_markers[id_epsilon]+List_linestyle[id_fN]+List_colors[id_impl], 
                        markerfacecolor='none', 
                        label=List_implementation[id_impl]+ ", epsilon= "+str(List_epsilon[id_epsilon]) + ", f(N): "+List_fN_name[id_fN])    
    pltB.plot(List_N, np.ones(len(List_N)), "k-", label="1")
    pltB.legend()

    # Save and show figure
    mng = plt.get_current_fig_manager()
    mng.resize(*mng.window.maxsize())
    figure = plt.gcf()
    plt.show()
    figure.savefig(bench_result_path+'/mean_time_vs_pbl_size.png', format='png',bbox_inches = "tight") 

    # cout
    # print("Plotting : " + bench_result_path + plot_result_file + " : success!")
    print("++++++++++++++++++++++++++")

    ############################ Comparaison avec GBm ############################
    # import json

    # # Parameters
    # Result_path = bench_result_path
    # Result_file = "/bench_hmatrix_build.json"
    # list_N = List_N
    # list_implementation = List_implementation
    # list_threads = [1]
    # Is_log_scale = True

    # # Variables
    # mat_real_time = np.zeros((len(list_N), len(list_implementation), len(list_threads)))
    # mat_cpu_time = np.zeros((len(list_N), len(list_implementation), len(list_threads)))

    # # read json
    # f = open(Result_path+Result_file)
    # data = json.load(f)
    # f.close()

    # for BM in data['benchmarks']:
    #     for id_N in range(len(list_N)):
    #         for id_impl in range(len(list_implementation)):
    #             for id_threads in range(len(list_threads)):
    #                 if (BM['name'] == "FT_Generator/BM_"+str(list_implementation[id_impl])+"/N:"+str(list_N[id_N])+"/threads:"+str(list_threads[id_threads])+"_median"):
    #                     mat_real_time[id_N, id_impl, id_threads] = BM['real_time'] / 1e6 # conversion en secondes
    #                     mat_cpu_time[id_N, id_impl, id_threads] = BM['cpu_time'] / 1e6 
                        
    # # Plots               
    # id_threads = 0   
    # id_epsilon = 2

    # # Plot real time
    # pltA = plt.subplot(211)
    # pltA.set_title("Building median real time as a function of pbl size")
    # pltA.set_ylabel("Real time (us)")
    # for id_impl in range(len(list_implementation)):
    #     pltA.plot(list_N, mat_real_time[:, id_impl, id_threads], '+-', label=list_implementation[id_impl]+ ", thread(s): "+str(list_threads[id_threads]))
    #     pltA.errorbar(List_N, Mat_mean_time[id_epsilon, :, id_impl], yerr=Mat_stddev[id_epsilon, :, id_impl], fmt=List_markers[id_epsilon]+'-'+List_colors[id_impl], markerfacecolor='none', capsize=3, label=List_implementation[id_impl]+ ", epsilon: "+str(List_epsilon[id_epsilon]))

    # if Is_log_scale:
    #     plt.yscale('log')
    #     plt.xscale('log')
    # pltA.legend()

    # # Plot cpu time
    # pltu = plt.subplot(212)
    # pltu.set_title("Building cpu time as a function of pbl size")
    # pltA.set_xlabel("Problem size N")
    # pltu.set_ylabel("CPU time (us)")
    # for id_impl in range(len(list_implementation)):
    #     pltu.plot(list_N, mat_cpu_time[:, id_impl, id_threads], '+-', label=list_implementation[id_impl]+ ", thread(s): "+str(list_threads[id_threads]))
    #     pltu.errorbar(List_N, Mat_mean_time[id_epsilon, :, id_impl], yerr=Mat_stddev[id_epsilon, :, id_impl], fmt=List_markers[id_epsilon]+'-'+List_colors[id_impl], markerfacecolor='none', capsize=3, label=List_implementation[id_impl]+ ", epsilon: "+str(List_epsilon[id_epsilon]))
        
    # if Is_log_scale:
    #     plt.yscale('log')
    #     plt.xscale('log')
    # pltu.legend()


    # # Save and show figure
    # mng = plt.get_current_fig_manager()
    # mng.resize(*mng.window.maxsize())
    # figure = plt.gcf()
    # plt.show()
    # figure.savefig(Result_path+'/Time_vs_PblSize.png', format='png',bbox_inches = "tight") 

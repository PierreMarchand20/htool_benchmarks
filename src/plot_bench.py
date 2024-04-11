#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import os , time, sys
import matplotlib.animation as animation
from math import exp, sqrt, ceil, log2
import pandas as pd

################### functions #######################

# Initialize custom parameters
def custom_parameters(List_epsilon, List_dim, List_thread, List_algo_type, bench_type):
    # Set optional subset of parameters for plots. Must be left equal to their corresponding list if not used
    SubList_dim       = List_dim
    SubList_thread    = List_thread
    SubList_algo_type = List_algo_type

    if bench_type == "build" or bench_type == "matrix_product":
        SubList_epsilon = [1E-10, 1E-8, 1E-4]    
        data            = 'time (s)' # data name for plots (corresponds to the column name in the csv file)

    elif bench_type == "factorization":
        SubList_epsilon = [1E-10, 1E-7, 1E-4]    
        data            = 'factorization_time (s)' # data name for plots (corresponds to the column name in the csv file)
    SubList = [SubList_epsilon, SubList_dim, SubList_thread, SubList_algo_type]
    
    # Set log exponent for plots
    log_exponent = 2

    return log_exponent, data, SubList, SubList_epsilon, SubList_dim, SubList_thread, SubList_algo_type

# Import data from csv file
def read_csv(List_case_type):
    # check number of arguments
    if len(sys.argv) < 2 or len(sys.argv) > 4:
        print("\t Error: Invalid number of arguments - Usage: python plot_test.py <csv_name> [ Optional: <is_log_scale>={True|False(default)} ]")
        sys.exit(1)

    # set csv_name and is_log_scale
    csv_name = sys.argv[1]
    is_log_scale = False
    if len(sys.argv) > 2 and sys.argv[2] != '0':
        is_log_scale = True
    
    # get bench_type
    List_bench_type = ["build", "matrix_product", "factorization"]
    bench_type = [bm for bm in List_bench_type if bm in csv_name]
    if bench_type:                 # if bench_type is not empty
        bench_type = bench_type[0] # transform bench_type from list to string
    else:
        print("\t Error: Invalid csv name - Usage: <csv_name> = bench_hmatrix_<bench_type>_vs_<case_type>.csv with <bench_type> in [build, matrix_product, factorization]")
        sys.exit(1)
        
    # get case_type
    case_type = [case for case in List_case_type if case in csv_name]
    if case_type:                # if case_type is not empty
        case_type = case_type[0] # transform case_type from list to string
    else:
        print("\t Error: Invalid csv name - Usage: <csv_name> = bench_hmatrix_<bench_type>_vs_<case_type>.csv with <case_type> in [pbl_size, thread, ratio]")
        sys.exit(1)
    
    # import data 
    try:
        df = pd.read_csv(os.path.join("../build", csv_name))
        print("Success: Data imported from", csv_name)
    except FileNotFoundError:
        print("\t Error: File not found - " + csv_name)
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print("\t Error: Empty data file - " + csv_name)
        sys.exit(1)
    except pd.errors.ParserError:
        print("\t Error: Unable to parse data file - " + csv_name)
        sys.exit(1)
    
    # Remove spaces from columns names
    df.columns = [col.strip() for col in df.columns]
    
    return df, bench_type, case_type, csv_name, is_log_scale

# Make directories and copy csv file in them
def copy_csv_file(bench_type, case_type, List_dim, List_thread, csv_name):
    # Create result directories if they doesn't exist
    result_path = "../output/Hmatrix_"+bench_type+"_vs_"+case_type+"/"

    if case_type == List_case_type[0]:
        result_path += (str(List_dim)).replace(' ', '', 1) # add List_* without the first space
    elif case_type == List_case_type[1]:
        result_path += str(List_thread).replace(' ', '', 1) # add List_* without the first space
    elif case_type == List_case_type[2]:
        result_path += str(List_dim) + "_" + str(List_thread)
    
    result_path=result_path.replace(' ', '_').replace('[', '').replace(']', '')
    os.makedirs(result_path, exist_ok=True)

    # Copy csv file to result directory if it doesn't exist
    csv_path = os.path.join(result_path, csv_name)

    if os.path.exists(csv_path):
        print(f"The file {csv_name} already exists in the directory {result_path}. It will be overwritten.")
        choice = input("=> Do you want to continue? (y/n) ")
        if choice.lower() != "y":
            print("Operation cancelled.")
            exit()
        else:
            os.system(f'cp {os.path.join("../build", csv_name)} {csv_path}')
    else:
        os.system(f'cp {os.path.join("../build", csv_name)} {csv_path}')

    print(f"Sucess: file {csv_name} has been copied to the directory {result_path}.")
    return result_path

# Check if elements of SubList are in their corresponding list
def check_sublist(List_parameters, SubList): # Check if elements of SubList are in their corresponding list
    for i in range(len(SubList)):
        for j in range(len(SubList[i])):
            if SubList[i][j] not in List_parameters[i]:
                print(f"\t Error: Value '{SubList[i][j]}' is not in list: {List_parameters[i]}")
                sys.exit(1)

    print("Success: Subset of parameters loaded")

# Load sets of epsilon, dim_pbl, number_of_threads, algo_type values
def load_parameters(df):
    List_epsilon    = df['epsilon'].unique()
    List_dim        = df['dim'].unique()
    List_thread = df['number_of_threads'].unique()
    List_algo_type  = df['algo_type'].unique()
    List_parameters = [List_epsilon, List_dim, List_thread, List_algo_type]
    
    return List_parameters, List_epsilon, List_dim, List_thread, List_algo_type

# Print parameters in console
def print_parameters(List_epsilon, List_dim, List_thread, List_algo_type, SubList, is_log_scale):
    print("\n====================Parameters==================")
    print("List_epsilon : "   + str(List_epsilon))
    print("List_dim : "       + str(List_dim))
    print("List_thread : "    + str(List_thread))
    print("List_algo_type : " + str(List_algo_type))
    print("SubList : "        + str(SubList))
    print("")
    print("is_log_scale : "   + str(is_log_scale))
    print("================================================\n")
    
# Load a subset of data from csv file
def load_data(df, data, SubList_epsilon, SubList_dim, SubList_thread, SubList_algo_type):
    Tab_mean_time = np.zeros((len(SubList_epsilon), len(SubList_dim), len(SubList_thread), len(SubList_algo_type)))
    Tab_stddev = np.zeros_like(Tab_mean_time)

    for id_BM in range(len(df[:])):
        try:
            id_epsilon = list(SubList_epsilon).index(df['epsilon'][id_BM])
            id_N       = list(SubList_dim).index(df['dim'][id_BM])
            id_thread  = list(SubList_thread).index(df['number_of_threads'][id_BM])
            id_impl    = list(SubList_algo_type).index(df['algo_type'][id_BM])
        except ValueError:
            continue
        
        if df['id_rep'][id_BM] == ' mean':
            Tab_mean_time[id_epsilon, id_N, id_thread, id_impl] = df[data][id_BM]
        if df['id_rep'][id_BM] == ' stddev':
            Tab_stddev[id_epsilon, id_N, id_thread, id_impl] = df[data][id_BM]

    return Tab_mean_time, Tab_stddev
        
# Plots
def plot_bench(data, SubList, SubList_epsilon, SubList_dim, SubList_thread, SubList_algo_type, Tab_mean_time, Tab_stddev, List_markers, List_colors, List_case_type, case_type, is_log_scale, log_exponent):
    # Plot real time  
    pltA = plt.subplot(211)
    pltA.set_ylabel("Real " + data)  
    
    # Set plot title
    if case_type == List_case_type[0]:
        pltA.set_title("Mean real time as a function of problem size: "+str(SubList_dim))
    elif case_type == List_case_type[1]:
        pltA.set_title("Mean real time as a function of number of threads: "+str(SubList_thread))
    elif case_type == List_case_type[2]:
        pltA.set_title("Mean real time as a function of ratio: "+str(SubList_dim)+" / "+str(SubList_thread))
        
    # Plot SubList
    for id_impl in range(len(SubList_algo_type)):
        for id_epsilon in range(len(SubList_epsilon)):  
            if case_type == List_case_type[0]:
                for id_thread in range(len(SubList_thread)):
                    pltA.errorbar(SubList_dim, Tab_mean_time[id_epsilon, :, id_thread, id_impl],
                                      yerr=Tab_stddev[id_epsilon, :, id_thread, id_impl],
                                      fmt=List_markers[id_epsilon]+'-'+List_colors[id_impl],
                                      markerfacecolor='none', capsize=3, 
                                      label=SubList_algo_type[id_impl]+ ", epsilon= "+str(SubList_epsilon[id_epsilon]))
            elif case_type == List_case_type[1]:
                for id_dim in range(len(SubList_dim)):
                    pltA.errorbar(SubList_thread, Tab_mean_time[id_epsilon, id_dim, :, id_impl],
                                      yerr=Tab_stddev[id_epsilon, id_dim, :, id_impl],                        fmt=List_markers[id_epsilon]+'-'+List_colors[id_impl],
                                      markerfacecolor='none', capsize=3,                        label=SubList_algo_type[id_impl]+ ", epsilon= "+str(SubList_epsilon[id_epsilon]))
    
    if case_type == List_case_type[2]: 
        Sub_mean_time = np.zeros((len(SubList_epsilon), len(SubList_dim), len(SubList_algo_type)))
        Sub_stddev = np.zeros_like(Sub_mean_time)
        for id_dim in range(len(SubList_dim)):
                Sub_mean_time[:, id_dim, :] = Tab_mean_time[:, id_dim, id_dim, :]
                Sub_stddev[:, id_dim, :] = Tab_stddev[:, id_dim, id_dim, :]         
        for id_impl in range(len(SubList_algo_type)):
            for id_epsilon in range(len(SubList_epsilon)):        
                pltA.errorbar(SubList_dim, Sub_mean_time[id_epsilon, :, id_impl], 
                            yerr=Sub_stddev[id_epsilon, :, id_impl], 
                            fmt=List_markers[id_epsilon]+'-'+List_colors[id_impl],
                            markerfacecolor='none', capsize=3, 
                            label=SubList_algo_type[id_impl]+ ", epsilon= "+str(SubList_epsilon[id_epsilon]))

    # Plot ideal
    if case_type == List_case_type[0]:
            Array_N = np.array(SubList_dim)
            pltA.plot(Array_N, Array_N*np.log2(Array_N) / (Array_N[-1]*np.log2(Array_N)[-1]) * Tab_mean_time[0, -1, 0, 0], '--', label='NlogN')
            pltA.plot(Array_N, Array_N*Array_N / (Array_N[-1]*Array_N[-1]) * Tab_mean_time[0, -1, 0, 0], '--', label='N^2')
    elif case_type == List_case_type[1]: 
        pltA.plot(SubList_thread, Tab_mean_time[-1, 0, 0, 0]/SubList_thread, "k-", label="Ideal for" + SubList_algo_type[0]+ ", epsilon= "+str(SubList_epsilon[-1]))
    elif case_type == List_case_type[2]:
        pltA.plot(SubList_dim, np.ones(len(SubList_dim)), "k-", label="1")

    
    # log scale case
    if is_log_scale:
        plt.yscale('log')
        plt.xscale('log')
        
    pltA.legend(bbox_to_anchor=(1, 1.02))
    #############################

    # Plot normalized time
    if case_type == List_case_type[0]:
        pltB = plt.subplot(212)
        pltB.set_title("Rescaling by f(N)")
        pltB.set_xlabel("Degree of freedom (N)")  
        pltB.set_ylabel("Normalized " + data + " / f(N)")
        
        # set rescaling functions and corresponding linestyles
        List_fN = [np.array(np.array(List_dim)* np.log2(np.array(List_dim))), 
                (np.array(List_dim))**(3/2),
                np.array(List_dim) * (np.log2(np.array(List_dim)))**log_exponent]

        List_fN_name = [r'$N\log_2(N)$', r'$N^{3/2}$', r'$N\log_2^{'+str(log_exponent)+'}(N)$'] # r for raw string for latex
        List_linestyle = [':', '-', '--', '-.']

        # Plot SubList
        for id_impl in range(len(SubList_algo_type)):
            for id_epsilon in range(len(SubList_epsilon)):
                for id_fN in range(len(List_fN)):  
                    for id_thread in range(len(SubList_thread)):
                        pltB.plot(List_dim,
                                    np.array(Tab_mean_time[id_epsilon, :, id_thread, id_impl]/List_fN[id_fN]) / Tab_mean_time[id_epsilon, 0, id_thread, id_impl] * List_fN[id_fN][0],
                                    List_markers[id_epsilon]+List_linestyle[id_fN]+List_colors[id_impl],                         
                                    markerfacecolor='none',                        
                                    label=SubList_algo_type[id_impl]+ ", epsilon= "+str(SubList_epsilon[id_epsilon]) + ", f(N): "+List_fN_name[id_fN])
                            
        pltB.plot(List_dim, np.ones(len(List_dim)), "k-", label="1")
        pltB.legend(bbox_to_anchor=(1, 1.02))

    # Save and show figure
    figure = plt.gcf()
    plt.gcf().set_size_inches(16, 8)
    # plt.show()
    figure.savefig(result_path+'/mean_time_vs_'+case_type+'.png', format='png', dpi=300 ,bbox_inches = "tight")
    print("++++++++++++++++++++++++++") 

########################### main ###########################

# List of case types against which to compare the benchmark
List_case_type = ["pbl_size", "thread", "ratio"]

# Import data from csv file
df, bench_type, case_type, csv_name, is_log_scale = read_csv(List_case_type)

# Load sets of epsilon, dim_pbl, number_of_threads, algo_type values
List_parameters, List_epsilon, List_dim, List_thread, List_algo_type = load_parameters(df)

# Initialize custom parameters
log_exponent, data, SubList, SubList_epsilon, SubList_dim, SubList_thread, SubList_algo_type = custom_parameters(List_epsilon, List_dim, List_thread, List_algo_type, bench_type)

# Check if elements of SubList are in their corresponding list
check_sublist(List_parameters, SubList)

# Print parameters in console
print_parameters(List_epsilon, List_dim, List_thread, List_algo_type, SubList, is_log_scale)

# Make directories and copy csv file in them
result_path = copy_csv_file(bench_type, case_type, List_dim, List_thread, csv_name)

# Load a subset of data from csv file
Tab_mean_time, Tab_stddev = load_data(df, data, SubList_epsilon, SubList_dim, SubList_thread, SubList_algo_type)

# Custom colors and markers for plots     
List_markers = ['o', 's', 'D', 'v', '^', 'h', 'd'] # for epsilon
List_colors  = ['b', 'r', 'g', 'c', 'm', 'y', 'k'] # for implementation
        
# Do plots
plot_bench(data, SubList, SubList_epsilon, SubList_dim, SubList_thread, SubList_algo_type, Tab_mean_time, Tab_stddev, List_markers, List_colors, List_case_type, case_type, is_log_scale, log_exponent)



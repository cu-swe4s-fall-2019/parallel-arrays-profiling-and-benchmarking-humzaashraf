import get_data
import math_lib
import matplotlib.pyplot as plt
import sys
import os 

def boxplot(L,groups,title_name,out_file_name):
    plt.figure(figsize=(10,5),dpi=300)
    plt.boxplot(L)
    tick_numbers = [i for i in range(len(groups))]
    plt.xticks(tick_numbers,groups,rotation='vertical')
    plt.ylabel('Signal')
    #title_str_mean = math_lib.list_mean(L)
    #title_str_stdev = math_lib.list_stdev(L)
    plt.title(title_name) 
    plt.savefig(out_file_name+'.png')
    pass

def histogram(L, out_file_name,col_index):
    plt.figure()
    plt.hist(L)
    plt.xlabel('Signal')
    plt.ylabel('Counts')
    title_str_mean = math_lib.list_mean(L)
    title_str_stdev = math_lib.list_stdev(L)
    plt.title('mean: '+title_str_mean+ ' | stdev: '+title_str_stdev)
    plt.savefig(out_file_name+'.png')
    pass

def combo(L, out_file_name,col_index):
    fig, (ax1, ax2) = plt.subplots(2)
    ax1.hist(L)
    ax2.boxplot(L)
    title_str_mean = math_lib.list_mean(L)
    title_str_stdev = math_lib.list_stdev(L)
    fig.suptitle('mean: '+title_str_mean+ ' | stdev: '+title_str_stdev)
    plt.savefig(out_file_name+'.png')
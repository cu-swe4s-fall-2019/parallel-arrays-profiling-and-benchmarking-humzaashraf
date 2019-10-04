import sys
import os 
import numpy as np


def read_stdin_col(col_num):
    col_data = []
    List = sys.stdin.read()
    split_List = List.rstrip().split('\n')
    array_List = np.array(split_List)
    try:
        index = int(col_num-1)
    except TypeError:
        print('input must be integer')
        sys.exit(1)
    for line in array_List:
        Type = line.split()
        try:
            idx = int(Type[index])
        except IndexError:
            print('index out of bounds')
            sys.exit(1)
        col_data.append(idx)
    return(col_data)
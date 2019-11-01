import data_viz
import argparse
import sys
import os
from hash_tables_humzaashraf import hash_tables as ht

def main():

    def linear_search(key, L):
        for i in range(len(L)):
            if L [i] == key:
                return i
        return -1
    
    def binary_search(key, L):
        lo = -1
        hi = len(L) 
        while (hi - lo > 0):
            mid = (hi + lo) // 2

            if key == L[mid]:
                return(L[mid])
            elif (key < L[mid]):
                hi = mid
            else:  
                lo = mid

        return -1

    parser = argparse.ArgumentParser(
             description='arguments',
             prog='input arg')

    parser.add_argument('--gene', type=str, help='name of gene', required = True)
    parser.add_argument('--group_type', type=str, help='SMTS or SMTDS', required=True)
    parser.add_argument('--data_file_name', type=str, help='GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.acmg_59.gct', required=True)
    parser.add_argument('--sample_info_file_name', type=str, help='GTEx_Analysis_v8_Annotations_SampleAttribut.txt.txt', required=True) 
    parser.add_argument('--output_file',type=str, help='test', required=True)

    args = parser.parse_args()
    sample_info_file_name = args.sample_info_file_name 
    group_col_name = args.group_type
    title_name = args.gene
    out_title = args.output_file
    data_file_name = args.data_file_name
    sample_id_col_name = 'SAMPID'

    sample_info_header = None
    samples = []
    for l in open(sample_info_file_name):
        if sample_info_header is None:
            sample_info_header = l.rstrip().split('\t')
        else:
            samples.append(l.rstrip().split('\t'))

    group_col_idx = linear_search(group_col_name, sample_info_header)
    sample_id_col_idx = linear_search(sample_id_col_name, sample_info_header)

    groups = []
    members = []

    mk_hash_table = ht.ChainedHash(1000, ht.h_rolling)

    for i in range(len(samples)):
        sample = samples[i]
        sample_name = sample[sample_id_col_idx]
        curr_group = sample[group_col_idx]
        search_result = mk_hash_table.search(curr_group)
        if search_result is None:
            mk_hash_table.add(curr_group,[sample_name])
            groups.append(curr_group)
        else:
            search_result.append([sample_name])        

    version = None
    dim = None
    data_header = None
    gene_name_col = 1

    group_counts = [ [] for i in range(len(groups)) ]

    for l in open(data_file_name, 'rt'):
        if version == None:
            version = l
            continue

        if dim == None:
            dim = [int(x) for x in l.rstrip().split()]
            continue

        if data_header == None:
            data_header = l.rstrip().split('\t')
            continue

        testcase = []

        A = l.rstrip().split('\t')
        if A[gene_name_col] == title_name:
            mk_hash_table1 = ht.ChainedHash(20000, ht.h_rolling)
            for group_idx in range(2, len(data_header)):
                mk_hash_table1.add(data_header[group_idx],int(A[group_idx]))
        #print(mk_hash_table.T)
        #print(mk_hash_table1.T)

    group_counts = []
    for i in range(len(groups)):
        counts = []
        sample_id = mk_hash_table.search(groups[i])
        if sample_id is None:
            continue
        for j in sample_id:
            sample_counts = mk_hash_table1.search(str(j))
            if sample_counts is None:
                continue
            counts.append(sample_counts)
        group_counts.append(counts)

    print(group_counts)
    print(counts)

    data_viz.boxplot(group_counts,groups,title_name,out_title)


if __name__ == '__main__':
    main()
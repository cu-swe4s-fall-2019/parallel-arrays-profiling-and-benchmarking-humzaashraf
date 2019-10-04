import data_viz
import sys
import os

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


def main():

    parser = argparse.ArgumentParser(
             description='arguments',
             prog='input arg')

    parser.add_argument('--gene', type=int, help='name of gene')
    parser.add_argument('--group_type', type=str, help='SMTS or SMTDS')
    parser.add_argument('--data_file_name', type=str, help='GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.acmg_59.gct')
    parser.add_argument('--sample_info_file_name', type=str, help='GTEx_Analysis_v8_Annotations_SampleAttribut.txt.txt') 

    args = parser.parse_args()

    args.data_file_name=='GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.acmg_59.gct'
    args.sample_info_file_name=='GTEx_Analysis_v8_Annotations_SampleAttribut.txt.txt'

    group_col_name = args.group_type
    sample_id_col_name = 'SAMPID'

    title_name = args.gene

    samples = []
    sample_info_header = None
    for l in open(sample_info_file_name):
        if sample_info_header == None:
            sample_info_header = l.rstrip().split('\t')
        else:
            samples.append(l.rstrip().split('\t'))

    group_col_idx = linear_search(group_col_name, sample_info_header)
    sample_id_col_idx = linear_search(sample_id_col_name, sample_info_header)

    print(sample_id_col_idx)
    print(group_col_idx)

    groups = []
    members = []

    for i in range(len(samples)):
        sample = samples[i]
        sample_name = sample[sample_id_col_idx]
        curr_group = sample[group_col_idx]
        curr_group_idx = linear_search(curr_group, groups)

        if curr_group_idx == -1:
            curr_group_idx = len(groups)
            groups.append(curr_group)
            members.append([])
        members[curr_group_idx].append(sample_name)

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

        A = l.rstrip().split('\t')

        if A[gene_name_col] == title_name:
            for group_idx in range(len(groups)):
                for member in members[group_idx]:
                    member_idx = linear_search(member, data_header)
                    if member_idx != -1:
                        group_counts[group_idx].append(int(A[member_idx]))
            break 

    data_viz.boxplot(group_counts,groups,title_name,'test')


if __name__ == '__main__':
    main()
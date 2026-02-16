#!/usr/bin/env python
##Certainly! Here's a refactored version of the script with improved readability and added comments for clarity:

import sys
import re
import os
from operator import itemgetter

def merge_to_txt(input_folder, nor_type='nor', dis_merge=3):
    # Get list of files in the input folder
    files = [os.path.join(input_folder, f) for f in os.listdir(input_folder)]
    data = {}

    # Read data from each file and store in a dictionary
    for file_path in files:
        sample_id = re.sub(input_folder + r'/', '', file_path)
        sample_id = re.sub('.5prim.txt$', '', sample_id)
        sample_id = re.sub('^Res.', '', sample_id)
        sample_id = re.sub('.[1,2]$', '', sample_id)
        
        data[sample_id] = {}
        
        with open(file_path, 'r') as file:
            file.readline()  # Skip header line
            line = file.readline()
            
            while line:
                nor_r, r, chr, pos = line.strip('\n').split('\t')
                coordinate = (chr, int(pos))
                
                if coordinate not in data[sample_id]:
                    data[sample_id][coordinate] = [nor_r, r]
                else:
                    print('Duplicate coordinate found!')
                
                line = file.readline()

    # Collect all unique insertions
    all_insertions = []
    for sample in data:
        all_insertions.extend(data[sample].keys())
    
    all_insertions = unique(all_insertions)
    all_insertions = sorted(all_insertions, key=lambda x: (x[0], int(x[1])))
    clustered_insertions = distance_merge(all_insertions, dis_merge)

    # Print header
    sample_names = sorted(data.keys())
    print('Chrom\tStart\tEnd\tJoinedIns\t' + '\t'.join(sample_names))
    
    # Print data for each cluster
    for cluster in clustered_insertions:
        print('\t'.join([cluster[0][0], str(cluster[0][1]), str(cluster[0][1])]), end='\t')
        print(cluster, end='')
        
        for sample in sample_names:
            if len(cluster) == 1:
                if cluster[0] in data[sample]:
                    print('\t' + data[sample][cluster[0]][0 if nor_type == 'nor' else 1], end='')
                else:
                    print('\t0', end='')
            else:
                value = sum(float(data[sample][ele][0]) if nor_type == 'nor' else int(data[sample][ele][1]) for ele in cluster if ele in data[sample])
                print('\t' + str(value), end='')
        
        print('')

def distance_merge(coordinates, distance_cutoff=1):
    def cluster_insertions(input_list):
        distances = [1000000]  # Large number for the first element
        for i in range(1, len(input_list)):
            distances.append(input_list[i][1] - input_list[i-1][1])
        
        clusters = []
        current_cluster = [input_list[0]]
        
        for i in range(1, len(distances)):
            if distances[i] <= distance_cutoff:
                current_cluster.append(input_list[i])
            else:
                clusters.append(current_cluster)
                current_cluster = [input_list[i]]
            
            if i == len(distances) - 1:
                clusters.append(current_cluster)
        
        return clusters
    
    final_clusters = []
    current_chr = 'aa'
    current_list = []
    
    for coordinate in coordinates:
        if coordinate[0] == current_chr:
            current_list.append(coordinate)
        else:
            if current_list:
                final_clusters.extend(cluster_insertions(current_list))
            current_chr = coordinate[0]
            current_list = [coordinate]
        
        if coordinate == coordinates[-1]:
            final_clusters.extend(cluster_insertions(current_list))
    
    return final_clusters

def unique(items):
    return list(set(items))

if __name__ == '__main__':
    folder, dis_merge = sys.argv[1:]
    ##It's better to start from the raw counts
    ##Instead of changing the code above just hard coding nor_type to False
    ##The normalized data is still in the individual files and setting
    ##nor_type="nor" will select those values instead of the raw counts
    nor_type=False
    merge_to_txt(folder, nor_type, int(dis_merge))

### Changes Made:
# 1. **Renamed Functions and Variables**: Improved readability by using more descriptive names.
# 2. **Added Comments**: Explained the purpose of each section and function.
# 3. **Simplified Logic**: Streamlined some parts of the code for better clarity.
# 
# Let me know if you need any further modifications or explanations!

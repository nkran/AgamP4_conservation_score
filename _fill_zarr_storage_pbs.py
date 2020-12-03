import zarr
import numpy as np
import pandas as pd
import os
import itertools

import multiprocessing as mp
from functools import partial

PATH = '/home/conservation/'

def parallelize(data, func, flog, num_of_processes=8):
    data_split = np.array_split(data, num_of_processes)
    pool = mp.Pool(num_of_processes)
    pool.map(func, data_split)
    pool.close()
    pool.join()
    flog.write('finished {}\t{}\t{}'.format(selected_identity, selected_genome, len(data)))

    return


def run_on_subset(func, data_subset):
    return data_subset.apply(func, axis=1)


def parallelize_on_rows(data, func, flog, num_of_processes=8):
    return parallelize(data, partial(run_on_subset, func), flog, num_of_processes)


def fill_score_matrix(line):
    # get chromosome and genome path
    path = '/'.join(['scores', line.chromosome, selected_genome, 'identity'])
    # modify orthogonal selection
    root[path].oindex[identities.index(selected_identity), slice(line.start, line.end)] = line.score


if __name__ == '__main__':
    eph_data_dir = os.path.join(os.environ['EPHEMERAL'], 'conservation', 'data')
    chr_sizes = pd.read_csv(os.path.join(PATH, 'data/ref/AgamP4_chr_sizes.tsv'), names=['chrom', 'size'], sep='\t')

    identities =    ['49_50', '35_50', '29_30', '27_30', '30_30', '50_50', '45_50', '48_50']
    genomes =       ['AaegL5', 'AalbS2', 'AaraD1', 'AatrE3', 'AchrA1', 'AcolM1', 'AculA1', 'AdarC3', 
                     'AdirW1', 'AepiE1', 'AfarF2', 'AfunF1', 'AmacM1', 'AmelC2', 'AmerM2', 'AminM1', 
                     'AquaS1', 'AsinC2', 'AsteI2', 'CpipJ2', 'DmelP6']
    chromosomes =   ['2L', '2R', '3L', '3R', 'X']
    
    subjobs = list(itertools.product(genomes, identities))
    
    for job in subjobs:
        selected_genome, selected_identity = job
        root = zarr.open('data/conservation_detailed.zarr', mode='r+')
        df_file = 'data/AgamP4.{}_d{}.csv'.format(selected_genome, selected_identity)

        with open(os.path.join('log_fill_zarr.txt'), "a") as f:
            if os.path.exists(df_file):
                
                names = df_file.split('.')
                identity = names[1].split('_d')[1]
                data = pd.read_csv(df_file, sep=',',
                                            low_memory=False,
                                            engine='c',
                                            header=0, 
                                            names=["x","chromosome","start","end","length","strand","second.seqnames",
                                                   "second.start","second.end","second.width","second.strand","score","cigar"])
                                            
                parallelize_on_rows(data, fill_score_matrix, f, 10)
                f.write('{}\t{}\t{}\n'.format(selected_identity, selected_genome, len(data)))
            else:
                f.write('NA {}\t{}\t{}\n'.format(selected_identity, selected_genome, df_file))

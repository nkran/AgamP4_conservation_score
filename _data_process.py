import zarr
import numpy as np
import pandas as pd
import allel
import h5py
import itertools
import multiprocessing as mp
from functools import partial
from sklearn import preprocessing

WINDOW_SIZE = 20


def create_joined_stacked_group():
    chr_sizes = pd.read_csv('data/chr_sizes.txt', sep='\t')

    for chrom in chr_sizes['chrom']:
        chr_size = int(chr_sizes[chr_sizes['chrom'] == chrom]['size'])
        ic = root.create_group(chrom)
        ic.zeros('stack', shape=(21, chr_size), dtype=np.dtype('f4'), chunks=(1, round(chr_size / 128)))
        ic.zeros('stack_norm', shape=(21, chr_size), dtype=np.dtype('f4'), chunks=(1, round(chr_size / 128)))


def process_SNP_data(selected_chromosome):
    z = root['{}/score'.format(selected_chromosome)][0, :]

    # Load Ag1000g variation data
    data_ag1000g = h5py.File('ag1000g/variation/ag1000g.phase2.ar1.pass.{}.h5'.format(selected_chromosome), mode='r')
    variants = allel.VariantChunkedTable(data_ag1000g[selected_chromosome]['variants'], 
                                        names=['POS'],
                                        index='POS')

    # SNP data is 1-based
    snp_positions = variants['POS'][:] - 1
    pos_array = np.zeros(len(z))
    pos_array[snp_positions] = 1
    pos_roll = pd.DataFrame(pos_array).rolling(WINDOW_SIZE, center=True).apply(lambda x: np.sum(x) / WINDOW_SIZE, raw=True).fillna(0)
    
    root['{}/score'.format(selected_chromosome)][1, :] = np.array(pos_roll)


def calculate(chromosome):
    '''
        - 2L
            - snp_density
                - 1 row
            - scores_stacked
                - stack_avr (1 row)
                - stack_norm (1 row)
                - stack_norm_avr (1 row)
                - stack_norm_avr_snp (1 row)
                - stack_norm_avr_snp_minmax (1 row)
            - stack
                - 21 rows - species
            - stack_norm
                - 21 rows - species
    '''
    scores = []

    stack = root['{}/stack'.format(chromosome)]

    stack_avr = root['{}/stack'.format(chromosome)].sum(axis=0) / (len(genomes_ord_names))
    distances = np.array(phyl_distances['distance']).reshape(-1, 1)
    stack_norm = stack * distances


    stack_norm_avr = stack_norm.sum(axis=0) / (len(genomes_ord_names))
    snp_density = root['{}/snp_density'.format(chromosome)][0,:]
    stack_norm_avr_snp = stack_norm_avr * ((1 - snp_density) / (1 + snp_density))
    stack_norm_avr_snp_minmax = preprocessing.minmax_scale(stack_norm_avr_snp)
    
    root['{}/stack_norm'.format(chromosome)][:, :] = stack_norm

    # assign
    root['{}/score'.format(chromosome)][0, :] = stack_avr
    root['{}/score'.format(chromosome)][1, :] = stack_norm
    root['{}/score'.format(chromosome)][2, :] = stack_norm_avr
    root['{}/score'.format(chromosome)][3, :] = stack_norm_avr_snp
    root['{}/score'.format(chromosome)][4, :] = stack_norm_avr_snp_minmax

    print('finished', chromosome)


if __name__ == '__main__':
    identities =    ['27_30', '29_30', '30_30', '35_50', '45_50', '48_50', '49_50', '50_50']
    genomes =       ['AaegL5', 'AalbS2', 'AaraD1', 'AatrE3', 'AchrA1', 'AcolM1', 'AculA1', 'AdarC3', 
                    'AdirW1', 'AepiE1', 'AfarF2', 'AfunF1', 'AmacM1', 'AmelC2', 'AmerM2', 'AminM1', 
                    'AquaS1', 'AsinC2', 'AsteI2', 'CpipJ2', 'DmelP6']
    chromosomes =   ['2L', '2R', '3L', '3R', 'X']

    phyl_distances = pd.read_csv('data/phyl_distances.tsv', sep='\t', names=['genome', 'distance'])
    genomes_ord = list(phyl_distances['genome'])
    genomes_ord_idx = [genomes.index(g) for g in genomes_ord]
    genomes_ord_names = [g for g in genomes_ord]

    subjobs = list(itertools.product(genomes, chromosomes))

    root = zarr.open('data/conservation_detailed.zarr', mode='r+')

    # create_joined_stacked_group()

    pool = mp.Pool(8)
    pool.map(calculate, chromosomes)
    pool.close()
    pool.join()
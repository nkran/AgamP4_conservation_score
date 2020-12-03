import sys
import os
import argparse
import requests
import h5py

import numpy as np
import pandas as pd


def parse_args():
    """
    Parse command line parameters

    Args:
      args ([str]): command line parameters as list of strings

    Returns:
      :obj:`argparse.Namespace`: command line parameters namespace
    """
    parser = argparse.ArgumentParser(
        description="AgamP4 conservation score")
    parser.add_argument(
        '-m',
        '--mode',
        dest="mode",
        help="Mode (download, extract)",
        choices=['download', 'extract'],
        default='extract',
        type=str)
    parser.add_argument(
        '-r',
        '--region',
        dest='region',
        help='Genomic region (ie. 2R:48714500-48714700)',
        type=str)
    parser.add_argument(
        '-a',
        '--array',
        dest='array',
        help='Name of the array to access. (Cs, score, snp_density, stack, stack_norm, phyloP)',
        choices=['Cs', 'score', 'snp_density', 'stack', 'stack_norm', 'phyloP'],
        default='Cs',
        type=str)
    parser.add_argument(
        '-o',
        '--output',
        dest='output',
        help='Name of the output file where data will be saved',
        type=str)

    return parser.parse_args()


def parse_region(region_string):
    
    chr_sizes = pd.read_csv('data/chr_sizes.txt', sep='\t')

    if ':' in region_string:
        chromosome, positions = region_string.split(':')
        if chromosome not in chr_sizes['chrom'].values:
            print(f'Chromosome {chromosome} does not exist in the dataset.')
            print(f'Available chromosomes are: 2L, 2R, 3L, 3R and X.')
            exit()

        positions = positions.split('-')
        if len(positions) < 2:
            start = positions[0]
            end = None
        else:
            start, end = positions
    else:
        print('No chromosome defined. You need to define genomic region in the following format CHR:START-END or CHR:POS')
        exit()

    return (chromosome, start, end)


if __name__ == '__main__':
    args = parse_args()
    file_name = 'data/AgamP4_conservation.h5'

    if args.mode == 'download':
        link = 'https://zenodo.org/record/4304586/files/AgamP4_conservation.h5'
        
        with open(file_name, "wb") as f:
            print(f'Downloading {file_name}')
            response = requests.get(link, stream=True)
            total_length = response.headers.get('content-length')

            if total_length is None: # no content length header
                f.write(response.content)
            else:
                dl = 0
                total_length = int(total_length)
                for data in response.iter_content(chunk_size=4096):
                    dl += len(data)
                    f.write(data)
                    done = int(50 * dl / total_length)
                    finished = int(dl / total_length * 100)
                    sys.stdout.write("\rProgress: [%s%s] %s%%\t" % ('=' * done, ' ' * (50-done), finished))    
                    sys.stdout.flush()
    
    elif args.mode == 'extract':
        if not args.output:
            print('No output specified. Please use --output argument to set where data will be saved. Use --help for help.')
            exit()
        if not args.region:
            print('You need to specify genomic region by using --region argument. Use --help for help.')
            exit()
        else:
            chromosome, start, end = parse_region(args.region)
            if os.path.exists(file_name):
                with h5py.File('data/AgamP4_conservation.h5', mode='r+') as root:
                    start, end = int(start), int(end)
                    values = root[chromosome][args.array][:,start-1:end-1]
                    row_names = root[chromosome][args.array].attrs['rows']

                    df = pd.DataFrame(values.T)
                    df.columns = [row_names]
                    df['chromosome'] = chromosome
                    df['pos'] = np.arange(start, end)

                    cols = df.columns.tolist()
                    df = df.loc[:, cols[-2:] + cols[:-2]]

                    df.to_csv(args.output, sep='\t', index=False)
                    print(f'Saved to {args.output}')
            else:
                print('Dataset file does not exist. Download it by using --mode download.')
                exit()

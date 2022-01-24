# AgamP4 conservation score

Scripts used for generating, processing and fetching the genome conservation score (Cs) for _Anopheles gambiae_ (AgamP4) genome.

Publication: https://doi.org/10.3390/insects12020097

## Create and activate Conda environment
Make sure you have Conda installed on your system beforehand. You can read more about how to install it here: https://docs.conda.io/projects/conda/en/latest/user-guide/install

First, clone this repository.
```
git clone git@github.com:nkran/AgamP4_conservation_score.git
cd AgamP4_conservation_score/
```

Create an environment and install all needed dependencies:
```
conda create --name AgamP4-Cs python=3.7
conda activate AgamP4-Cs
pip install -r requirements.txt
```

## Download the dataset

To download the dataset you can use the following command:
```
python fetch_score.py -m download 
```
You can also download the dataset manually from https://doi.org/10.5281/zenodo.4304586. Please make sure you place the downloaded file in the `data/` folder.

## Fetch conservation score (Cs)
Extraction of values from the dataset can be done by using `fetch_score.py`

### Usage
```
usage: fetch_score.py [-h] [-m {download,extract}] [-r REGION]
                      [-a {Cs,score,snp_density,stack,stack_norm,phyloP}]
                      [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -m {download,extract}, --mode {download,extract}
                        Mode (download, extract)
  -r REGION, --region REGION
                        Genomic region (ie. 2R:48714500-48714700)
  -a {Cs,score,snp_density,stack,stack_norm,phyloP}, --array {Cs,score,snp_density,stack,stack_norm,phyloP}
                        Name of the array to access. (Cs, score, snp_density,
                        stack, stack_norm, phyloP)
  -o OUTPUT, --output OUTPUT
                        Name of the output file where data will be saved
```

### Storage
There are several different results saved in the HDF5 file.
- `stack` – the array consists of 21 rows representing stacked intervals with identity scores from CNEr (http://bioconductor.org/packages/release/bioc/html/CNEr.html) scans for each species analysed in the study. Order of reference genomes in the array is defined under ‘genomes’ attribute of the array. ‘species’ attribute contains the ordered list of species names. Columns represent positions in the chromosomes.

- `stack_norm` – copy of the stack array where each row with was normalised by multiplying the row contents with the phylogenetic distance between AgamP4 genome and the genome belonging to the row.

- `snp_density` – the array consists of 20bp sliding window of SNP density obtained from Ag1000g phase 2 dataset (https://www.malariagen.net/projects/ag1000g).

- `score` – the array with averaged scores across all 21 genomes where each row contains values from each step in the process of calculating the final conservation score. Row descriptions are stored under ‘rows’ attribute of the array.

- `phyloP` - phyloP scores calculated for the multiple genomes alignment using PHAST (http://compgen.cshl.edu/phast/) 

- `Cs` – the array with final conservation score presented in the paper

### Examples

Cs can be extracted with `--mode extract` argument.

```
python fetch_score.py -m extract -r 2R:48714500-48714700 -a Cs -o results.tsv
```

Extracted values will be saved in a tab separated file specified with `-o` argument (ie. `results.tsv`).

__Output__
```
chromosome	pos	        Cs
2R	        48714500	0.026579915
2R	        48714501	0.024763243
2R	        48714502	0.01920125
2R	        48714503	0.012823918
2R	        48714504	0.007067339
2R	        48714505	0.0037936924
2R	        48714506	0.016306616
2R	        48714507	0.015620683
2R	        48714508	0.013961536
2R	        48714509	0.022508306
2R	        48714510	0.021567238
2R	        48714511	0.020158349
...
```

## Contact information
For additional information, help and bug reports please send an email to n.kranjc@imperial.ac.uk.

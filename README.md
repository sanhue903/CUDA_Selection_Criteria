# Selection_Criteria
Algorithm implementation of selection criteria based on Hyperloglog and SuperMinHash sketches to speed up genomic similarity

## Installation

This tool requires the following dependencies:
- NCBI Datasets command-line tools. https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/
- Seqan library. https://packages.seqan.de/
- Dashing library. https://github.com/dnbaker/dashing-binaries

To install them follow the steps below:

### NCBI Datasets command-line tools

### Seqan library

### Dashing library

[Agregar como juntar esto con este reopositorio para poder correr los códigos]

## Use

The implementation of the selection criteria is in the `src` folder. The file `build.cpp` contains the sketch construction step. Once compiled, this program is used as follows:
```
./build/build -l filelist -t nthreads -a aux_memory -c criterion
```
Where
- `-l` option recieves a txt file containing the path to the .fna.gz files that will be processed. One line for every path.
- `-t` option recieves the number of threads to run the program.
- `-a` option recieves the additional memory, in bytes, that the criterion will use for each sequence.
- `-c` option recieves the criterion to use. The aviable options are `hll_a`, `hll_an` and `smh_a`.

Once the command is run, the primary hll sketches and auxiliary structures will be saved along the genomic .fna.gz files.

The `selection.cpp` file has the selection algorithm. Once compiled, this program is used as follows:
```
./build/selection -l filelist -t nthreads -h tau -a aux_memory -c criterion
```
Where
- `-l` option recieves a txt file containing the path to the .fna.gz files that will be processed. One line for every path.
- `-t` option recieves the number of threads to run the program.
- `-h` option recieves the similarity threshold `tau`.
- `-a` option recieves the additional memory, in bytes, that the criterion will use for each sequence.
- `-c` option recieves the criterion to use. The aviable options are `hll_a`, `hll_an` and `smh_a`.
- 
## Experiments

The `experiments` folder has the datasets and implementation to evaluate the different selecion criteria $`hll_a`$, $`hll_{a,n}`$ and $`smh_a`$. This folder contains the following:
1. `Manifest` folder: It contains the manifests .zip files to download the exact same datasets used in [PAPER??]. The NCBI Datasets command-line tools is required to downlad the datasets through the manifest files (See Installation).
2. 

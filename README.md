[![DOI](https://zenodo.org/badge/984350098.svg)](https://doi.org/10.5281/zenodo.15537863)


# Selection Criteria

Implementation of selection criteria based on Hyperloglog and SuperMinHash sketches to speed up genomic similarity

## Installation (Linux: Ubuntu-24.04)

This tool requires the following dependencies. They can be installed via `apt-get`:
- make
- g++ 
- libz-dev

Now, clone this repository into target folder:
```
git clone https://github.com/sanhue903/CUDA_Selection_Criteria
```

We need to install the following libraries:
- Seqan library. https://seqan.readthedocs.io/en/main/Infrastructure/Use/Install.html#infra-use-install or https://packages.seqan.de/
- Sketch library. https://github.com/dnbaker/sketch

To install them follow the steps below:

### Seqan library

```
wget https://packages.seqan.de/seqan-library/seqan-library-2.4.0.tar.xz
tar xvf seqan-library-2.4.0.tar.xz
cd seqan-library-2.4.0
sudo cp -r include/seqan /usr/include
sudo cp -r share/doc/seqan /usr/share/doc
```

### Sketch library

Clone the Sketch library using the following command inside the folder created when cloning this repository (by default its name is ``CUDA_Selection_Criteria``):
```
git clone --depth 1 --branch v0.19.0 https://github.com/dnbaker/sketch
```

## Use

The implementation of the selection criteria is in the `src` folder, we also provide implementation for experiments to evaluate the selection criteria which are in the `experiments/src` folder. To compile the .cpp files use the Makefile provided, simply run `make`. The file `build_sketch.cpp` contains the sketch construction step. This program is used as follows:
```
./build/build_sketch -l filelist -t nthreads -a aux_memory -c criterion
```
Where
- `-l` option recieves a txt file containing the path to the .fna.gz files that will be processed. One line for every path.
- `-t` option recieves the number of threads to run the program.
- `-a` option recieves the additional memory, in bytes, that the criterion will use for each sequence.
- `-c` option recieves the criterion to use. The aviable options are `hll_a`, `hll_an` and `smh_a`.

Once the command is run, the primary hll sketches and auxiliary structures will be saved along the genomic .fna.gz files.

Where
- `-l` option recieves a txt file containing the path to the .fna.gz files that will be processed. One line for every path.
- `-t` option recieves the number of threads to run the program.
- `-h` option recieves the similarity threshold `tau`.
- `-a` option recieves the additional memory, in bytes, that the criterion will use for each sequence.
- `-c` option recieves the criterion to use. The aviable options are `hll_a`, `hll_an` and `smh_a`.

The `selection_main` file has the selection algorithm with CUDA. This program is used as follows:

```
./build/selection_main -l filelist -h tau -a aux_memory 
```

## Experiments

We provide implementation of experiments to evaluate time of comparisions. To run these experiments using a certain file list of sequences, first we have to build the associated sketches with the `build_sketch` program. These implementation are in the `experiments/src` folder. 

```
./build/time_smh_cuda -l filelist -h tau -m buckets_smh -R reps
```
Where
- `-l` option recieves a txt file containing the path to the .fna.gz files that will be processed. One line for every path.
- `-t` option recieves the number of threads to run the program.
- `-h` option recieves the similarity threshold `tau`.
- `-m` option recieves number of buckets of the auxiliar SuperMinHash sketch m.
- `-R` option recieves number of reps of the experiment.

The output is the time it takes to retrieve the similar pairs using the `CB+smh_a` criterion, we also include the results of the `CB` criterion and with no criterion (baseline case).

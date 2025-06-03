[![DOI](https://zenodo.org/badge/984350098.svg)](https://doi.org/10.5281/zenodo.15537863)


# Selection_Criteria

Algorithm implementation of selection criteria based on Hyperloglog and SuperMinHash sketches to speed up genomic similarity

## Installation (Linux)

This tool requires the following dependencies. They can be installed via `apt-get`:
- make
- g++ 
- libz-dev

Now, clone this repository into target folder:
```
git clone https://github.com/AlvaroGuzmanChacon/Selection_Criteria
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

```
git clone --depth 1 --branch v0.19.0 https://github.com/dnbaker/sketch
```

## Use

The implementation of the selection criteria is in the `src` folder. To compile the .cpp files use the Makefile provided. The file `build_sketch.cpp` contains the sketch construction step. Once compiled, this program is used as follows:
```
./build/build_sketch -l filelist -t nthreads -a aux_memory -c criterion
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

## Datasets

The `datasets` folder has the datasets contains the manifest.zip files in the `Manifests` folder, along with a test dataset in the `test_influenzaA` folder. The `test_influenzaA` containts 10 compressed .fna files from the Influenza A GenBank dataset. These files are listed in the `test_influeza_filelist.txt` text file. To test the `build_sketch` or `selection.cpp` we can use this test datasets. Let's consider the next example:
```
./build/build_sketch -l test_influeza_filelist.txt -t 8 -a 256 -c hll_an
./build/selection -l test_influeza_filelist.txt -t 8 nthreads -h 0.8 -a 256 -c hll_an > results.txt
```
The first command would create the hll sketches of 256 bytes associated with the files specified in `test_influeza_filelist.txt`, using 8 threads. Then, the second command would use the `hll_an` criterion to try and find the pairs of sequences (between the 10 sequences on `test_influeza_filelist.txt`) with a Jaccard similarity greater than or equal to 0.8, using 8 threads and using the hll sketches of 256 bytes previously constructed.

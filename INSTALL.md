# IIBacFinder installation (detailed)

These instructions are for Linux (Ubuntu/CentOS) with Anaconda/Miniconda installed.

## 1. Clone the source code

```bash
git clone https://github.com/ZhangDengwei/IIBacFinder.git
cd IIBacFinder
```

## 2. Create the conda environment

```bash
conda env create -f environment.yml
conda activate IIBacFinder
```

If conda is too slow, use mamba or micromamba instead:

```bash
mamba env create -f environment.yml
mamba activate IIBacFinder
```

or

```bash
micromamba env create -f environment.yml
micromamba activate IIBacFinder
```

## 3. Install the IIBacFinder package

```bash
pip install -e .
```

Gene calling is performed by Pyrodigal with a minimum CDS length of 10 amino acids (30 nt). Use `-m single` or `-m meta` to select the Pyrodigal model.

## 4. Install SignalP 6.0

SignalP 6.0 is license-restricted and the model tarball is not distributed by IIBacFinder. Obtain the official fast-model tarball from DTU.

With the `IIBacFinder` environment active, run:

```bash
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
signalp6-register signalp-6.0i.fast.tar.gz
```

Confirm the installation:

```bash
signalp6 --help
```

## 5. Download the large data resources

The data archive is stored at https://zenodo.org/records/22211103. The downloader retrieves the archive and its SHA-256 file, verifies the checksum, and unpacks the data:

```bash
bash download_data.sh
```

## 6. Verify the installation

```bash
IIBacFinder -h
IIBacFinder -i test_fasta/ -o test_prediction
```

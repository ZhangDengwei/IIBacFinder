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

Official DTU SignalP 6.0 page:
https://services.healthtech.dtu.dk/services/SignalP-6.0/

Download request for the fast model (`signalp-6.0i.fast.tar.gz`):
https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=6.0&packageversion=6.0i&platform=fast

Submit the request at the DTU page above, then download the tarball from the link sent to your email.

If SignalP registration or execution fails with:

```text
ImportError: /lib64/libstdc++.so.6: version `GLIBCXX_3.4.29' not found
```

run the following command before `signalp6-register` and before running IIBacFinder:

```bash
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
```

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

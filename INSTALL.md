# IIBacFinder installation (detailed)

These instructions are for Linux (Ubuntu/CentOS) with Anaconda/Miniconda installed. Review the contents of `environment.yml`, `requirements.txt`, and the scripts in this repository before installing.

## 1. Clone the source code

```bash
git clone https://github.com/ZhangDengwei/IIBacFinder.git
cd IIBacFinder
```

The GitHub repository is the only source of IIBacFinder code. Do not install IIBacFinder from an unverified tar.gz archive.

## 2. Create the conda environment

```bash
conda env create -f environment.yml
conda activate IIBacFinder
```

This creates an environment named `IIBacFinder`. The file `environment.yml` lists every conda-managed dependency, so users can inspect the full dependency graph before creating the environment. No pre-built environment archive is used.

## 3. Install Python packages

```bash
pip install -r requirements.txt
```

This installs `fastaparser`, `networkx`, `peptides`, `pandas_ml`, the pinned PyTorch build required by the bundled NLPPrecursor module, the matching fastai snapshot, and amPEP. Pyrodigal is used for gene calling with a minimum CDS length of 10 aa. For GPU users, follow the PyTorch instructions to select the appropriate torch build before installing `requirements.txt`.

## 4. Install the IIBacFinder package

```bash
pip install -e .
```

The editable install registers the `IIBacFinder` command-line entry point.

Gene calling is performed by Pyrodigal with a minimum CDS length of 10 amino acids (30 nt). Use `-m single` or `-m meta` to select the Pyrodigal model; the output directory is `pyrodigal_out/`.

## 5. Install R packages

The gene-context module uses ampir for AMP prediction, and region plotting uses several R packages.

```bash
conda install -c conda-forge r-ggplot2 r-dplyr r-cowplot r-ggrepel r-lattice r-gridextra
Rscript -e 'install.packages(c("ampir", "gggenes", "ggiraphExtra"), repos = "https://cloud.r-project.org")'
```

If a package reports a missing dependency during the first run, install it from CRAN with the same command.

## 6. Install SignalP 6.0

SignalP 6.0 is license-restricted and is not distributed by IIBacFinder. Users must obtain the tarball directly from DTU.

1. Submit a request at https://services.healthtech.dtu.dk/services/SignalP-6.0/ and select the fast model.
2. DTU sends the download link by email. Download the tarball (for example, `signalp-6.0h.fast.tar.gz`) and do not extract it.
3. With the `IIBacFinder` environment active, register the downloaded file using its actual path:

```bash
signalp6-register /path/to/signalp-6.0h.fast.tar.gz
```

4. Confirm that the executable is available:

```bash
signalp6 --help
```

## 7. Download the large data resources

The following directories are too large for GitHub and are provided as a data-only Zenodo archive: `domains/`, `models/`, `AMP_database/` (DIAMOND index), `scripts/cleavage_pred/training_data/`, and `test_fasta/`.

Before running the downloader, set the published URL and SHA-256 checksum:

```bash
export IIBACFINDER_DATA_URL="https://zenodo.org/records/<DATA_RECORD>/files/IIBacFinder_data.tar.gz?download=1"
export IIBACFINDER_DATA_SHA256="<SHA256_OF_IIBacFinder_data.tar.gz>"
bash download_data.sh
```

The script verifies the SHA-256 checksum before unpacking, so users can confirm that the archive matches the published manifest.

## 8. Verify the installation

```bash
IIBacFinder -h
IIBacFinder -i test_fasta/ -o test_prediction
```

The second command runs the included test genomes as a smoke test.

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

This creates an environment named `IIBacFinder`. The file `environment.yml` installs all conda-managed dependencies, including the R packages `ampir` and `gggenes`, and uses its embedded `pip` section to install the Python packages listed in `requirements.txt`. No pre-built environment archive is used. If `requirements.txt` is updated after an environment already exists, run `conda env update -f environment.yml` to apply the changes.

## 3. Install the IIBacFinder package

```bash
pip install -e .
```

The editable install registers the `IIBacFinder` command-line entry point.

Gene calling is performed by Pyrodigal with a minimum CDS length of 10 amino acids (30 nt). Use `-m single` or `-m meta` to select the Pyrodigal model; the output directory is `pyrodigal_out/`.

## 4. Install SignalP 6.0

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

## 5. Download the large data resources

The following directories are too large for GitHub and are provided as a data-only Zenodo archive: `domains/`, `models/`, `AMP_database/` (DIAMOND index), `scripts/cleavage_pred/training_data/`, and `test_fasta/`.

The data archive is stored at https://zenodo.org/records/22211103. Run the downloader directly; it downloads both `IIBacFinder_data.tar.gz` and its published SHA-256 file, verifies the checksum, and unpacks the archive:

```bash
bash download_data.sh
```

To use a different mirror, override `IIBACFINDER_DATA_URL` and `IIBACFINDER_DATA_SHA256_URL` before running the script.

## 6. Verify the installation

```bash
IIBacFinder -h
IIBacFinder -i test_fasta/ -o test_prediction
```

The second command runs the included test genomes as a smoke test.

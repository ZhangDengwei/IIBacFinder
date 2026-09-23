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
bash setup_env_activation.sh
conda deactivate
conda activate IIBacFinder
```

This creates an environment named `IIBacFinder`. The file `environment.yml` installs all conda-managed dependencies, including the R packages `ampir` and `gggenes`, the `signalp6` registration placeholder from the `predector` channel, and the tested `scikit-learn=1.2.0` version required by the amPEP model. Its embedded `pip` section installs the Python packages listed in `requirements.txt`. The `setup_env_activation.sh` step installs a conda activation hook so that `$CONDA_PREFIX/lib` is added to `LD_LIBRARY_PATH` automatically on every activation.

No pre-built environment archive is used. If `requirements.txt` or `environment.yml` is updated after an environment already exists, run `conda env update -f environment.yml` and then re-run `bash setup_env_activation.sh` while the environment is active.

## 3. Install the IIBacFinder package

```bash
pip install -e .
```

The editable install registers the `IIBacFinder` command-line entry point.

Gene calling is performed by Pyrodigal with a minimum CDS length of 10 amino acids (30 nt). Use `-m single` or `-m meta` to select the Pyrodigal model; the output directory is `pyrodigal_out/`.

## 4. Install SignalP 6.0

SignalP 6.0 is license-restricted and the model tarball is not distributed by IIBacFinder. The `predector::signalp6` package installed by `environment.yml` provides the `signalp6` and `signalp6-register` commands as a registration placeholder. Users must still obtain the official fast-model tarball directly from DTU.

1. Submit a request at https://services.healthtech.dtu.dk/services/SignalP-6.0/ and select the fast model.
2. DTU sends the download link by email. Download the tarball (for example, `signalp-6.0h.fast.tar.gz`) and do not extract it.
3. With the `IIBacFinder` environment active, register the downloaded file using the provided wrapper. The wrapper forces the conda `libstdc++.so.6` ahead of the older system library, which avoids the `GLIBCXX_3.4.29` import error:

```bash
bash register_signalp6.sh /path/to/signalp-6.0i.fast.tar.gz
```

If the release letter differs from `6.0h`, a filename warning may appear, but registration should continue when the file is the official fast-model tarball.

4. Confirm that the executable is available:

```bash
signalp6 --help
```

### If SignalP registration still fails with `GLIBCXX_3.4.29 not found`

Install the actual C++ runtime packages (`libstdcxx` and `libgcc`, not only their `-ng` metapackages) and then run the wrapper again:

```bash
conda activate IIBacFinder
micromamba install -y -p "$CONDA_PREFIX" -c conda-forge --force-reinstall \
  'libstdcxx-ng=16.2.0' 'libstdcxx=16.2.0' \
  'libgcc-ng=16.2.0' 'libgcc=16.2.0'
bash setup_env_activation.sh
hash -r
strings "$CONDA_PREFIX/lib/libstdc++.so.6" | grep GLIBCXX_3.4.29
bash register_signalp6.sh /path/to/signalp-6.0i.fast.tar.gz
```

If `micromamba` is unavailable, replace it with `mamba` or `conda`. The wrapper prints a diagnostic if `$CONDA_PREFIX/lib/libstdc++.so.6` is missing or does not contain `GLIBCXX_3.4.29`. The registration command is safe to rerun after a failed attempt because it unregisters the previous partial installation first.

### If the amPEP model reports a scikit-learn version error

The bundled amPEP model was trained with scikit-learn 1.2.0. Reinstall the matching version:

```bash
micromamba install -y -p "$CONDA_PREFIX" -c conda-forge --force-reinstall 'scikit-learn=1.2.0'
```

### If the R package `ampir` is missing

If `micromamba list` shows `r-ampir` but R reports that the package is not installed, force-reinstall it:

```bash
micromamba install -y -p "$CONDA_PREFIX" -c conda-forge --force-reinstall 'r-ampir=1.1.0'
```

### If R reports that `ggplot2` is missing

`ampir` loads `ggplot2` as a dependency. Force-reinstall it if it is absent:

```bash
micromamba install -y -p "$CONDA_PREFIX" -c conda-forge --force-reinstall 'r-ggplot2=3.5.1'
Rscript -e 'library(ampir); library(ggplot2)'
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

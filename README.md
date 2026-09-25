# IIBacFinder

**IIBacFinder** (class II bacteriocin finder) is designed to detect class II bacteriocins (small, unmodified antimicrobial peptides).

**Workflow**:

<img src="https://github.com/ZhangDengwei/IIBacFinder/blob/main/IMGs/1.png" width="50%">

*Overview of IIBacFinder specializing in detecting unmodified bacteriocins*

![image](https://github.com/ZhangDengwei/IIBacFinder/blob/main/IMGs/2.png)

*Schematic flow of bacteriocin mining in IIBacFinder*

## Installation

### 1. Clone the source code

```bash
git clone https://github.com/ZhangDengwei/IIBacFinder.git
cd IIBacFinder
```

### 2. Create the conda environment

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

### 3. Install the IIBacFinder package

```bash
pip install -e .
```

### 4. Install SignalP 6.0

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
signalp6 --help
```

### 5. Download the large data resources

```bash
bash download_data.sh
```

### 6. Verify the installation

```bash
IIBacFinder -h
IIBacFinder -i test_fasta/ -o test_prediction
```

> **Important**
>
> Run the test data before analyzing your own data, in case any unexpected bugs arise.

## Running IIBacFinder

To confirm the successful installation and view all options, execute:

```bash
IIBacFinder -h
```

```text
usage: IIBacFinder [-h] -i INDIR [-e HMMEXCUTE] [-a AMPEP] [-r RCMD] [-s HMMSCAN] -o OUTDIR [-t THRESHOLD] [-p] [-m {single,meta}] [-v]

Detecting class II bacteriocins from genomes

optional arguments:
  -h, --help            show this help message and exit
  -i INDIR, --inDir INDIR
                        Input path to folder containing FASTA files, whose suffix could be in the list of ['.fas', '.fa', '.fasta', '.faa', '.fna']
  -e HMMEXCUTE, --hmmexcute HMMEXCUTE
                        The executable hmmsearch, default: hmmsearch
  -a AMPEP, --ampep AMPEP
                        The executable ampep, default: ampep
  -r RCMD, --Rcmd RCMD  The executable Rscript, default: Rscript
  -s HMMSCAN, --hmmscan HMMSCAN
                        The executable hmmscan, default: hmmscan
  -o OUTDIR, --outDir OUTDIR
                        The path to the folder storing output files
  -t THRESHOLD, --threshold THRESHOLD
                        Number of thresholds used, default: 20
  -p, --pyrodigal       Whether to perform gene calling with Pyrodigal (default: TRUE), toggle to close if input files are already protein FASTA.
  -m {single,meta}, --model {single,meta}
                        Select the Pyrodigal gene calling model (single or meta). Default is single.
  -v, --version         Print out the version and exit.
```

Key parameters:

- `-i`: Input directory where FASTA files for prediction are located.
- `-o`: Output directory.
- `-p`: Leave enabled when the input FASTA files are genome sequences. Disable this option when the input files are already protein FASTA files.
- `-m`: Gene calling model for Pyrodigal (`single` or `meta`).


## Running a demo:

```
IIBacFinder -i ./IIBacFinder/test_fasta/ -o test_prediction
```

Prediction results can be found in `test_prediction`.

## Output files

IIBacFinder generates output files explained below:

1. Intermediate folders
   
   - `pyrodigal_out`
     
     Gene prediction by Pyrodigal (minimum CDS length: 10 aa)
   
   - `prediction_domain`
     
     Prediction results based on precursor rules
   
   - `prediction_geneContext`
     
     Prediction results based on context gene rules
   
   - `region_annotation`
     
     Domain scanning for bactericoin gene clusters
   
   - `diamond_alingment`
     
     Blast results for predicted bacteriocin precursors against publicly available AMP sequences
   
   - `leader_prediction`
     
     Results of leader prediction by [`SignalP 6.0`](https://services.healthtech.dtu.dk/services/SignalP-6.0/) and modified [`NLPPrecursor`](https://github.com/magarveylab/NLPPrecursor)

2. Data output folders
   
   - `results`
     
     Prediction results of each input FASTA file
   
   - `region_plot`
     
     Visualizations of predicted bacteriocin gene clusters, including three files with formats in `.svg`, `.tsv`, and `.gbk`

3. Output summary files
   
   - `overall_result.tsv` (**most important**)
     
     - `CDs`: CDs ID in genome annotation files by Pyrodigal which could be found at `pyrodigal_out`
     
     - **`Rules`**: Prediction rule, one of `Domain` (prediction based on precursor rules), `GeneContext` (prediction based on context gene rules), and `Both` (prediction based on both precursor and context gene rules )
     
     - `Domain_rule`: The ID of self-build precursor domains which could be found at `IIBacFinder/domains/classII-related.hmm` 
     
     - `Domain_Evalue`: E-value of `hmmsearch` query results
     
     - `Domain_Bitscore`: Bitscore of `hmmsearch` query results
     
     - `Context_rule`: Context gene rules used for bacteriocin prediction, which could be found at `IIBacFinder/domains/rule.gene.context.hmm`
     
     - **`PFAM_domain`**: Domain annotation of precursor sequence against `PFAM` database, which could be found at `IIBacFinder/domains/Pfam-A.hmm`
     
     - **`NCBI_domain`**: Domain annotation of precursor sequence against `NCBI` database, which could be found at `IIBacFinder/domains/hmm_PGAP.LIB`
     
     - **`Sequence`**: Predicted bacteriocin precursor sequence
     
     - `Length`: Length of predicted precursor sequence
     
     - `Description`: Putative description of predicted precursor sequence
     
     - **`Potential_Leader_Type`**: Inferred leader type
     
     - `Genome`: Genome ID
     
     - `Region`: Predicted gene cluster regions whose details could be found in `region_plot``
     
     - `Contig`: Contig ID
     
     - `Start`: Start position of predicted precursor sequence
     
     - `End`: End position of predicted precursor sequence
     
     - `Strand`: Strand of predicted precursor sequence
     
     - `Partial_index`: Completeness of predicted precursor sequence, which was annotated by Pyrodigal
     
     - `Start_type`: Start codon of predicted precursor sequence, which was annotated by Pyrodigal
     
     - `RBS_motif`: RBS motif of predicted precursor sequence, which was annotated by Pyrodigal
     
     - `Including_elements`: Predicted elements associated with bacteriocin biosynthesis in the gene cluster region
     
     - `Uniq_ID`: Assigned a unique ID for each predicted precursor sequence
     
     - `leader_sec`: Predicted leader sequence of precursor sequence with *sec* type
     
     - `core_sec`: Predicted core sequence of precursor sequence with *sec* type
     
     - `leader_gg`: Predicted leader sequence of precursor sequence with double-glycine type
     
     - `core_gg`: Predicted core sequence of precursor sequence with double-glycine type
     
     - **`Predicted_mature_peptide`**: Predicted mature sequence of bacteriocin sequence
     
     - `Confidence`: Assigned confidence level
     
     - `Length__core`: Length of the predicted mature sequence
     
     - `Charge (pH=7)__core`: Charge of the predicted mature sequence annotated by [`peptides.py`](https://peptides.readthedocs.io/en/stable/index.html)
     
     - `Isoelectric_point__core`: Isoelectric point of the predicted mature sequence annotated by [`peptides.py`](https://peptides.readthedocs.io/en/stable/index.html)
     
     - `Molecular_weight (monoisotopic)__core`: Molecular weight of the predicted mature sequence annotated by [`peptides.py`](https://peptides.readthedocs.io/en/stable/index.html)
     
     - `Aliphatic_index__core`: Aliphatic index of the predicted mature sequence annotated by [`peptides.py`](https://peptides.readthedocs.io/en/stable/index.html)
     
     - `Boman__core`: Boman index of the predicted mature sequence annotated by [`peptides.py`](https://peptides.readthedocs.io/en/stable/index.html)
     
     - `Instability_index__core`: Instability index of the predicted mature sequence annotated by [`peptides.py`](https://peptides.readthedocs.io/en/stable/index.html)
     
     - `Hsp_len`: The length of best hit when querying predicted bacteriocin sequence against publically available AMP sequences 
     
     - `Hsp_identity`: The identity of the best hit
     
     - `Coverage_q`: The coverage of the best hit compared to predicted bacteriocin sequence. For example, 100.0 (54/54) represents 100% coverage with the best hit and predicted bacteriocin being 54AAs and 54AAs, respectively.
     
     - `Coverage_s`: The coverage of the best hit compared to AMP sequence. For example, 78.3 (54/69) represents 78.3% coverage with the best hit and known AMP being 54AAs and 69AAs, respectively.
     
     - `AMP_seq`: AMP sequence
     
     - `AMP_accession`: AMP accession ID in the publicly available database
     
     - `AMP_name`: AMP description in the database
     
     - `AMP_database`: AMP database, including [APD3](https://aps.unmc.edu/), [DRAMP](http://dramp.cpu-bioinfor.org/), [DBAASP](https://dbaasp.org/home), and [dbAMP2](https://awi.cuhk.edu.cn/dbAMP/) 
   
   - `all.precusors.gg.fa` and `all.precusors.gg.json`
     
     Predicted bacteriocin precursor sequences with a putative double-glycine leader
   
   - `all.precusors.sec.fa` and `all.precusors.sec.json`
     
     Predicted bacteriocin precursor sequences with a putative *sec* leader
   
   - `all.precusors.fa`
     
     Overall predicted bacteriocin precursor sequences

## After predicting

To leave the environment:

```bash
conda deactivate
```

To run IIBacFinder again, reactivate the environment:

```bash
conda activate IIBacFinder
```

## Dependencies

All dependencies are listed in `environment.yml` and `requirements.txt`. Third-party attribution and license information are in [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md). The scientifically critical dependencies are HMMER, DIAMOND, Pyrodigal (minimum CDS length 10 aa), SignalP 6.0, amPEP, ampir, and the modified NLPPrecursor module.

The environment pins the tested `scikit-learn=1.2.0` version required by the amPEP model, together with `r-ampir=1.1.0` and `r-ggplot2=3.5.1`.

## Large data resources

| Resource | Contents | Included in GitHub? |
| --- | --- | --- |
| `domains/` | Pfam and NCBI HMM databases | No, downloaded |
| `models/` | amPEP model and NLPPrecursor models | No, downloaded |
| `AMP_database/` | AMP FASTA sources and DIAMOND index | Source FASTA in repository; index downloaded |
| `scripts/cleavage_pred/training_data/` | NLPPrecursor models and training data | No, downloaded |
| `test_fasta/` | Example genomes for a smoke test | No, downloaded |

`package_data.sh` builds `IIBacFinder_data.tar.gz` and its `.sha256` file for upload to Zenodo. `download_data.sh` downloads both from https://zenodo.org/records/22211103, verifies the SHA-256 checksum, and unpacks the data.

## Change from earlier versions

APIN and AmpGram were removed in v1.2.0. The pre-built conda environment archive is no longer used; users create the environment from `environment.yml` and install the Python packages listed in `requirements.txt`.

## Notes:

- IIBacFinder may overlook certain precursors due to its prediction threshold, especially for glycine-type bacteriocins, which can sometimes have multiple precursors within a single gene cluster. Therefore, it is advisable to double-check the predicted gene cluster instead of relying solely on precursor prediction.

- As with all bioinformatics tools, it is important not to place complete trust in the predictions. Conducting a manual check of the precursor and context gene is advisable to exclude false positives.

## Reference:

Zhang D, Zou Y, Shi Y, et al. Systematically investigating and identifying bacteriocins in the human gut microbiome. Cell Genom. Published online August 26, 2025. doi:[10.1016/j.xgen.2025.100983](https://www.cell.com/cell-genomics/fulltext/S2666-979X(25)00239-3)

## License

IIBacFinder is distributed under GPL-3.0-only (see [`LICENSE`](LICENSE)). The bundled modified NLPPrecursor code under `scripts/cleavage_pred/` remains GPL-3.0 (see [`scripts/cleavage_pred/LICENSE`](scripts/cleavage_pred/LICENSE)).

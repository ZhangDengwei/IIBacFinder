# IIBacFinder

**IIBacFinder** (class II bacteriocin finder) is designed to detect class II bacteriocins (small, unmodified antimicrobial peptides).

**Workflow**:

![image](https://github.com/ZhangDengwei/IIBacFinder/blob/main/IMGs/1.png)

*Overview of IIBacFinder specializing in detecting unmodified bacteriocins*

![image](https://github.com/ZhangDengwei/IIBacFinder/blob/main/IMGs/2.png)

*Schematic flow of bacteriocin mining in IIBacFinder*

## Installation

1. Download the latest version of IIBacFinder from [zenodo](https://zenodo.org/records/14292149). Of note, do not clone the package directly from the GitHub repository as it is incomplete

```
wget -O IIBacFinder.tar.gz https://zenodo.org/records/14292149/files/IIBacFinder.tar.gz?download=1
tar -zxvf IIBacFinder.tar.gz
```

2. Download and clone the execution environment for IIBacFinder

```
# download
wget -O env_IIBacFinder.tar.gz https://zenodo.org/records/14292149/files/env_IIBacFinder.tar.gz?download=1
# clone
mkdir -p $path/env_IIBacFinder
tar -xzf env_IIBacFinder.tar.gz -C $path/env_IIBacFinder
source $path/env_IIBacFinder/bin/activate
conda unpack
```

`$path` is where the IIBacFinder environment will be unpacked and cloned.

3. install [signalp6](https://services.healthtech.dtu.dk/services/signalp-6.0#)

*Due to license restrictions, this recipe cannot distribute signalp6 directly.*

Please download signalp-6.0d.fast.tar.gz from:  
https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=6.0&packageversion=6.0h&platform=fast

After registering online, you will receive the link for package download via email, and then download the package locally. 

Assuming you have downloaded the package locally, then run the following command to complete the installation:  

```
signalp6-register signalp-*tar.gz
```

This will copy signalp6 into your conda environment. After this step, the installation of IIBacFinder will be complete.

*Let's say that you need to deactivate the environment after the prediction, you can run*

```
source $path/env_IIBacFinder/bin/deactivate
```

## Running IIBacFinder

To confirm the successful installation and view all options, execute the command below

```
python $PATH/IIBacFinder/scripts/predict.py -h
```

`$PATH` is the directory where `IIBacFinder` was placed. 

```
usage: IIBacFinder [-h] -i INDIR [-e HMMEXCUTE] [-a AMPEP] [-r RCMD] [-s HMMSCAN] -o OUTDIR [-t THRESHOLD] [-p] [-m {single,meta}] [-v]

Detecting class II bacteriocins from genomes

optional arguments:
  -h, --help            show this help message and exit
  -i INDIR, --inDir INDIR
                        Input path to folder containg FASTA file, whose suffix could in list of ['.fas', '.fa', '.fasta', '.faa', '.fna']
  -e HMMEXCUTE, --hmmexcute HMMEXCUTE
                        The excutive hmmsearch, defualt: hmmsearch
  -a AMPEP, --ampep AMPEP
                        The excutive ampep, defualt: ampep
  -r RCMD, --Rcmd RCMD  The excutive Rscript, defualt: Rscript
  -s HMMSCAN, --hmmscan HMMSCAN
                        The excutive hmmscan, defualt: hmmscan
  -o OUTDIR, --outDir OUTDIR
                        The path to the folder storing output files
  -t THRESHOLD, --threshold THRESHOLD
                        Number of threshols used, default: 20
  -p, --prodigal_short  Whether perform gene prediction using prodigal short (default: TRUE), toggle to close. This is indispensable for FASTA files.
  -m {single,meta}, --prodigal_p {single,meta}
                        Select procedure model (single or meta) in 'prodigal_short', identical to the parameter '-p'. Default is single
  -v, --version         Print out the version and exit.
```

Key parameters:

- -i: Input directory where FASTA files for prediction are located.

- -o: Output directory

- -p: If the input FASTA files are genome sequences, this parameter should remain at its default setting to run gene prediction using `prodigal_short`. If the input FASTA files are prediction results from `prodigal_short`, this parameter should be specified to skip the `prodigal_short` prediction step.

- -m: The prediction model for `prodigal_short` prediction, see for details at [prodigal](https://github.com/hyattpd/prodigal/wiki/gene-prediction-modes).

## Running a demo:

```
python $PATH/IIBacFinder/scripts/predict.py -i $PATH/IIBacFinder/test_fasta/ -o test_prediction
```

Prediction results can be found in `test_prediction`.

## Output files

IIBacFinder generates output files explained below:

1. Intermediate folders
   
   - `prodigal_out`
     
     Gene prediction by `prodigal-short`
   
   - `prediction_domain`
     
     Prediction results based on precursor rules
   
   - `prediction_geneContext`
     
     Prediction results based on context gene rules
   
   - `region_annotation`
     
     Domain scanning for bactericoin gene clusters
   
   - `diamond_alingment`
     
     Blast results for predicted bacteriocin precursors against publicly available AMP sequences
   
   - `leader_prediction`
     
     Results of leader prediction by [`SignalP 6.0`](https://github.com/fteufel/signalp-6.0) and modified [`NLPPrecursor`](https://github.com/magarveylab/NLPPrecursor)

2. Data output folders
   
   - `results`
     
     Prediction results of each input FASTA file
   
   - `region_plot`
     
     Visualizations of predicted bacteriocin gene clusters, including three files with formats in `.svg`, `.tsv`, and `.gbk`

3. Output summary files
   
   - `overall_result.tsv` (**most important**)
     
     - `CDs`: CDs ID in genome annotation files by `prodigal-short` which could be found at `prodigal_out`
     
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
     
     - `Partial_index`: Completeness of predicted precursor sequence, which was annotated by `prodigal-short`
     
     - `Start_type`: Start codon of predicted precursor sequence, which was annotated by `prodigal-short`
     
     - `RBS_motif`: RBS motif of predicted precursor sequence, which was annotated by `prodigal-short`
     
     - `Including_elements`: Predicted elements associated with bacteriocin biosynthesis in the gene cluster region
     
     - `Uniq_ID`: Assigned a unique ID for each predicted precursor sequence
     
     - `leader_sec`: Predicted leader sequence of precursor sequence with *sec* type
     
     - `core_sec`: Predicted core sequence of precursor sequence with *sec* type
     
     - `leader_gg`: Predicted leader sequence of precursor sequence with double-glycine type
     
     - `core_gg`: Predicted core sequence of precursor sequence with double-glycine type
     
     - **`Predicted_mature_peptide`**: Predicted mature sequence of bacteriocin sequence
     
     - `Confidence`: Assigned confidence level
     
     - `Length__core`: Length of the predicted mature sequence
     
     - `Charge (pH=7)__core`: Charge of the predicted mature sequence annotated by [`peptides.py`]((https://peptides.readthedocs.io/en/stable/index.html)
     
     - `Isoelectric_point__core`: Isoelectric point of the predicted mature sequence annotated by [`peptides.py`]((https://peptides.readthedocs.io/en/stable/index.html)
     
     - `Molecular_weight (monoisotopic)__core`: Molecular weight of the predicted mature sequence annotated by [`peptides.py`]((https://peptides.readthedocs.io/en/stable/index.html)
     
     - `Aliphatic_index__core`: Aliphatic index of the predicted mature sequence annotated by [`peptides.py`]((https://peptides.readthedocs.io/en/stable/index.html)
     
     - `Boman__core`: Boman index of the predicted mature sequence annotated by [`peptides.py`]((https://peptides.readthedocs.io/en/stable/index.html)
     
     - `Instability_index__core`: Instability index of the predicted mature sequence annotated by [`peptides.py`]((https://peptides.readthedocs.io/en/stable/index.html)
     
     - `Hsp_len`: The length of best hit when querying predicted bacteriocin sequence against publically available AMP sequences 
     
     - `Hsp_identity`: The identity of the best hit
     
     - `Coverage_q`: The coverage of the best hit compared to predicted bacteriocin sequence. For example, 100.0 (54/54) represents 100% coverage with the best hit and predicted bacteriocin being 54AAs and 54AAs, respectively.
     
     - `Coverage_s`: The coverage of the best hit compared to AMP sequence. For example, 78.3 (54/69) represents 78.3% coverage with the best hit and known AMP being 54AAs and 69AAs, respectively.
     
     - `AMP_seq`: AMP sequence
     
     - `AMP_accession`: AMP accession ID in the publicly available database
     
     - `AMP_name`: AMP description in the database
     
     - `AMP_database`: AMP database, including [APD3](https://aps.unmc.edu/), [DRAMP](http://dramp.cpu-bioinfor.org/), [DBAASP]([Antimicrobial Peptide Database - DBAASP](https://dbaasp.org/home)), and [dbAMP2](https://awi.cuhk.edu.cn/dbAMP/) 
   
   - `all.precusors.gg.fa` and `all.precusors.gg.json`
     
     Predicted bacteriocin precursor sequences with a putative double-glycine leader
   
   - `all.precusors.sec.fa` and `all.precusors.sec.json`
     
     Predicted bacteriocin precursor sequences with a putative *sec* leader
   
   - `all.precusors.fa`
     
     Overall predicted bacteriocin precursor sequences

## Reference:

Zhang, Dengwei, et al. "Systematically investigating and identifying unmodified bacteriocins in the human gut microbiome." *bioRxiv* (2024): 2024-07.
# Third-party notices

IIBacFinder itself is distributed under GPL-3.0-only (see LICENSE).

IIBacFinder uses the following third-party tools. They are installed from their upstream distribution channels and are not bundled in the IIBacFinder repository, unless noted below. The exact installation commands are given in INSTALL.md.

| Tool | Purpose in IIBacFinder | Upstream source | License (see upstream repository) |
| --- | --- | --- | --- |
| HMMER | HMM search against bacteriocin-related domains | https://github.com/EddyRivasLab/hmmer | BSD-style license |
| DIAMOND | Alignment of candidate peptides against AMP databases | https://github.com/bbuchfink/diamond | GPL-3.0 |
| Pyrodigal | Gene calling with a minimum CDS length of 10 aa | https://github.com/althonos/pyrodigal | GPL-3.0 |
| SignalP 6.0 | Sec-dependent leader peptide prediction | https://services.healthtech.dtu.dk/services/SignalP-6.0/ | Academic/non-commercial license |
| amPEP | Antimicrobial peptide prediction in the gene-context module | https://github.com/tlawrence3/amPEPpy | See upstream repository |
| ampir (R) | Antimicrobial peptide prediction in the gene-context module | https://github.com/Legana/ampir | See upstream repository |
| NLPPrecursor (modified, bundled code) | Double-glycine leader cleavage prediction | https://github.com/magarveylab/NLPPrecursor | GPL-3.0 (see scripts/cleavage_pred/LICENSE) |
| fastai / PyTorch | Runtime dependencies of the modified NLPPrecursor module | https://github.com/fastai/fastai; https://pytorch.org | Apache-2.0 / BSD-style |

## Removed tools

APIN (https://github.com/zhanglabNKU/APIN) and AmpGram (https://github.com/michbur/AmpGram) were used in early development versions of IIBacFinder. Starting with v1.2.0, both tools have been removed from the code, the installer, and the protocol. They are no longer bundled, installed, invoked, or cited as IIBacFinder dependencies.

All current dependencies are declared in `environment.yml`, `requirements.txt`, and the protocol's Key Resources table so that users can review them before installation.

# -*- coding: utf-8 -*-
# @Author: Zhang Dengwei
# @Date:   2024-02-27 21:16:03
# @Last Modified by:   chem
# @Last Modified time: 2024-04-09 15:12:13


import peptides
import pandas


def calproperty(seqin):
    dict_property = dict()
    seqId, seq = seqin[0], seqin[1]
    peptide = peptides.Peptide(seq)
    pep_len = len(seq)
    """
    The aliphatic index of a protein was proposed in Ikai (1980). It is defined as the relative volume 
    occupied by aliphatic side chains (Alanine, Valine, Isoleucine, and Leucine). It may be regarded 
    as a positive factor for the increase of thermostability of globular proteins.
    """
    pep_aliphatic_index = round(peptide.aliphatic_index(), 2)
    """
    Compute the Boman (potential peptide interaction) index.
    The potential interaction index proposed by Boman (2003) is an index computed by averaging the 
    solubility values for all residues in a sequence. It can be used to give an overall estimate of the 
    potential of a peptide to bind to membranes or other proteins.
    Returns:
    The Boman index for the peptide. A value greater than 2.48 indicates that a protein has high binding potential.
    """
    pep_boman = round(peptide.boman(), 2)

    pep_charge = round(peptide.charge(pH=7), 2)
    # pep_counts = {k:v for k,v in peptide.counts().items() if v != 0}
    # pep_frequency = {k:v for k,v in peptide.frequencies().items() if v != 0}
    """
    Compute the instability index of a protein sequence.
    This function calculates the instability index proposed by Guruprasad et al (1990). This index predicts the 
    stability of a protein based on its dipeptide composition.
    Returns:
    float – The instability index of the peptide. A protein whose instability index is smaller than 40 is predicted 
    as stable, a value above 40 predicts that the protein may be unstable.
    """
    pep_instability_index = round(peptide.instability_index(), 2)
    """
    The isoelectric point (pI), is the pH at which a particular molecule or surface carries no net electrical charge.
    The pI is a variable that affects the solubility of the peptides under certain conditions of pH. When the pH of 
    the solvent is equal to the pI of the protein, it tends to precipitate and lose its biological function.
    """
    pep_isoelectric_point = round(peptide.isoelectric_point(pKscale="EMBOSS"), 2)
    """
    Compute the molecular weight of a protein sequence.
    """
    pep_molecular_weight = peptide.molecular_weight(average="monoisotopic")
    dict_property = {
        "Uniq_ID": seqId,
        "Length__core": pep_len,
        "Charge (pH=7)__core": pep_charge,
        "Isoelectric_point__core": pep_isoelectric_point,
        "Molecular_weight (monoisotopic)__core": pep_molecular_weight,
        "Aliphatic_index__core": pep_aliphatic_index,
        "Boman__core": pep_boman,
        "Instability_index__core": pep_instability_index,
    }
    return dict_property

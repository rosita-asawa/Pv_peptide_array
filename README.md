# Pv_peptide_array

Coordinating custom scripts for _PLoS Pathogens_ manuscript, published October 2026 <br>
https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1014621 <br>
**Title:** Proteome-wide characterization of _Plasmodium vivax_ antigens using a high-density peptide array <br>
**Authors:** Rosita R Asawa1, Brittany Hazzard1, Kieran Tebben1, John Tan2, Andrea Berry3, Niraj Tolia4, Jean Popovici5, David Serre1*

1-Institute for Genome Sciences, University of Maryland School of Medicine, Baltimore MD, USA <br>
2-Nimble Therapeutics Inc., Madison WI, USA <br>
3-Center for Vaccine Development, University of Maryland School of Medicine, Baltimore MD, USA<br>
4-Laboratory of Malaria Immunology and Vaccinology, National Institute of Allergy and Infectious Diseases, Bethesda MD, USA<br>
5-Institut Pasteur du Cambodge, Phnom Penh, Cambodia<br>
*correspondence to dserre@som.umaryland.edu



**Scripts included** <br>
_Peptides_immunogenicity.R_
- Order peptide array by amino acid position
  - Create list of repeat peptides across PVP01 proteins
- Immunogenicity set up: Symptomatic pv patients
  - Table S3: immunogenicity of symptomatic Pv patients across P01 genome
  - Numbers: immunoreactive peptides per sample
  - Numbers: immunoreactive regions per sample
  - Numbers: immunoreactive regions per gene
  - Numbers: immunoreactive segments/proteins shared across patients
  - Numbers: FDR for significance
  - Fig S3: histogram of reactive region length
  - Fig S3: IQR of histogram of reactive region length
- Immunogenicity set up: Asymptomatic pv patients
  - Table S4: immunogenicity of asymptomatic Pv patients across P01 genome
  - Numbers: asymptomatic vs random symptomatic segment reactivity rates
- Immunogenicity plots: Symptomatic Pv patients, time 1, repeat peptide marker
- Immunogenicity plots: Asymptomatic Pv patients, time 1, repeat peptide marker

_Peptides_enrichment.R_
- Symptomatic pv patient enrichment
  - Numbers: disorder enrichment within immunoreactive segments
  - Numbers: repeated peptide enrichment within immunoreactive segments
  - Numbers: NUP enrichment within immunoreactive proteins
  - Numbers: asparagine repeat enrichment within immunoreactive segments
  - Fig S5: correlation plot of kmers by antigenicity
- Asymptomatic Pv patient enrichment
  - Numbers: PIR enrichment within immunoreactive proteins

_Peptides_expression.R_ <br>
- Bulk RNA-seq of symptomatic pv patients
	- Numbers: normalize to CPM and rank expression
	- Fig 2: binned expression vs immunogenicity
	- Numbers: logistic regression of mean expression vs immunogenicity

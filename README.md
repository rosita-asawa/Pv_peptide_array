# Pv_peptide_array

Coordinating custom scripts for _PLoS Pathogens_ manuscript, accepted September 2026 <br>
https://www.biorxiv.org/content/10.64898/2026.03.17.712326v1 <br>
**Title:** Comprehensive characterization of _P. vivax_ antigens using high-density peptide array <br>
**Authors:** Rosita R Asawa1, Brittany Hazzard1, Kieran Tebben1, John Tan2, Andrea Berry3, Niraj Tolia4, Jean Popovici5, David Serre1*

1-Institute for Genome Sciences, University of Maryland School of Medicine, Baltimore MD, USA <br>
2-Nimble Therapeutics Inc., Madison WI, USA <br>
3-Center for Vaccine Development, University of Maryland School of Medicine, Baltimore MD, USA<br>
4-Laboratory of Malaria Immunology and Vaccinology, National Institute of Allergy and Infectious Diseases, Bethesda MD, USA<br>
5-Institut Pasteur du Cambodge, Phnom Penh, Cambodia<br>
*correspondence to dserre@som.umaryland.edu



**Scripts included** <br>
**Peptides_seroreactivity.R**
- Order peptide array by amino acid position
- Seroreactivity set up: Symptomatic pv patients
  - Table S3: seroreactivity of symptomatic Pv patients across P01 genome
  - numbers: seroreactive peptides per sample
  - numbers: seroreactive regions per sample
  - numbers: seroreactive regions per gene
  - numbers: seroreactive segments/proteins shared across patients
  - numbers: FDR for significance
  - Fig S3: histogram of reactive region length
  - Fig S3: IQR of histogram of reactive region length
- Seroreactivity set up: Asymptomatic pv patients
  - Table S4: seroreactivity of asymptomatic Pv patients across P01 genome
  - numbers: asymptomatic vs random symptomatic segment reactivity rates
- Seroreactivity plots: Symptomatic Pv patients, time 1, repeat peptide marker
- Seroreactivity plots: Asymptomatic Pv patients, time 1, repeat peptide marker

**Peptides_enrichment.R**
- Symptomatic pv patient enrichment
  - numbers: disorder enrichment within seroreactive segments
  - numbers: repeated peptide enrichment within seroreactive segments
  - numbers: NUP enrichment within seroreactive proteins
  - numbers: asparagine repeat enrichment within seroreactive segments
- Fig S5: correlation plot of kmers by antigenicity
  - Asymptomatic Pv patient enrichment
  - numbers: PIR enrichment within seroreactive proteins

**Peptides_expression.R** <br>
- Bulk RNA-seq of symptomatic pv patients
	- numbers: normalize to CPM and rank expression
	- Fig 2: binned expression vs seroreactivity
	- numbers: logistic regression of mean expression vs seroreactivity

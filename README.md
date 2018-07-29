# ISOseq-Protein-Annotation-Mapping
Aberrant splicing of mRNA transcripts may alter amino acid coding, peptide length, and/or peptide secondary structures of a gene.

Comparison of transdecoder identified ISOseq peptide sequences with the peptide sequences of Uniprot database of known peptide isoforms allows for the identifications of conserved peptide regions among prospective ISOseq peptides and know peptide isoforms for a specific gene. Some ISOseq peptides demonstrate 100% identity with known peptide isoforms and are therefore likely to be reads supporting that respective isoform. However, in the case of aberrant splicing, genetic mutation, transcription mishap, or sequencing artifacts, the transdecoder amino acid sequence may differ from the respective gene canonical isoform. Color shading of the multiple sequence alignment is added to observe the similarity of each amino acid in the ISOseq peptide sequence with respect to the canonical sequence: conserved amino acids (Blue), amino acids in the ISOseq sequence with a similar charge or hydrophobic property as the respective amino acid in the canonical isoform (Orange), amino acids in ISOseq that are neither conserved nor similar in chemical properties (Gray).

Partial deviations of ISOseq amino acid sequences from the respective "canonical" form are of great interest. Uniprot annotations based on the canonical peptide isoform (protein domains, topology, secondary structures, post-translational modifications, etc.) are added to the multiple sequence alignments to facilitate recognition of conserved, similar, or non-conserved amino acid sequences within annotated regions. 

ISOseq peptide sequences demonstrating amino acid conservation in regions crucial for protein synthesis suggests the corresponding ISOseq mRNA may lead to protein expression. However, variability observed in functional domain regions may lead to alterations in protein function, stability, and/or antigenicity.

ISO-PAM (ISOseq Protein Annotation Mapping)

1)	Import ENSGs and PBids from GTF of interest
2)	Query by “ENSG” (eg. ENSG00000143324.13)
3)	RNA sequences are obtained for each PBid from “annotated.dontScreenMatchAnnot.dontScreenProteinCoding.fa”
4)	Known RNA isoforms obtained from Biomart
5)	MSA performed RNA sequence alignment 
6)	Amino acid sequences are obtained for each PBid from “longest_orfs.pep”
7)	Known Peptide isoforms obtained from Biomart
8)	MSA performed amino acid sequence alignment 
9)	Identify canonical isoform based on Appris ranking
10)	Import Uniprot annotations from https://www.ebi.ac.uk/proteins/api
11)	Generate annotations file for LaTex: texshade 
12)	 RNA sequence figure
13)	Amino acid multiple sequence alignment figures with Annotations
14) Figerprint summary sequence alignment figure.pdf generated
15) Full sequence alignment figure.pdf generated

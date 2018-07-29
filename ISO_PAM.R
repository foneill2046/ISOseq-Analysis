#ISOseq Protein Annotation Mapping
#Load resources####
library(splitstackshape)
library(plotly)
library(bitops)
library(dplyr)
library(IRanges)
library(GenomicRanges)
library(AnnotationHub)
library(BSgenome.Hsapiens.UCSC.hg38)
library(stringr)
library(gsubfn)
library(readxl)
library(tidyr)
library(multtest)
library(genefilter)
library(reshape2)
require(SGSeq)
require(Rsamtools)
require(SummarizedExperiment)
library(msa)
require(seqRFLP)
require(biomaRt)
library(muscle)
library(plyr)
library(seqinr)

mart<- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="www.ensembl.org")

#Load PBid and ENSG from prepared file
txf_GTFdf<-readRDS(file ="/path/to/file/ISO_PAM/PBid_ENSG_info.rds")

#Query by ENG as found in ENSG and PBid file (PBid_ENSG_info.rds)
ENSGinput<-"ENSG00000143324.13"
#optional: file with PBids for subset analysis
#PBid_file_path<-"/path/to/file/ISO_PAM/TNFRSF8_PBid_Analysis.txt"

#Obtain PacBio PBid from GTF using ENSG
ENSGinputdf<-as.data.frame(ENSGinput)
names(ENSGinputdf)<-"ENSG"
PBids<-txf_GTFdf[txf_GTFdf$ENSG == ENSGinput, ]
PBids_df<-as.data.frame(PBids$PBid)
PBids_dfsplit<- cSplit(PBids_df, c(names(PBids_df)),sep = ",")
PBids_lst<-as.list(PBids_dfsplit)
PBids_unlst<-unlist(PBids_lst)
PBids_dfb<-as.data.frame(PBids_unlst)
PBids_dfc<-na.omit(PBids_dfb)
PBids_unique<-unique.data.frame(PBids_dfc)
names(PBids_unique)<-"PBids"
row.names(PBids_unique)<-c(1:nrow(PBids_unique))

#Write table of PBids for given ENSG
write.table(PBids_unique, file = "/path/to/file/ISO_PAM/TempPBids_unique", col.names = FALSE,quote = FALSE,row.names = FALSE)

##optional: Write table of PBids for given table
#PBid_table_analysis<-read.table(PBid_file_path, header = TRUE)
#write.table(PBid_table_analysis, file = "/path/to/file/ISO_PAM/TempPBids_unique", col.names = FALSE,quote = FALSE,row.names = FALSE)

#shell script
system("/path/to/file/ISO_PAM/TempPBidsUNIQUEtoEXACT.bash")
#Obtain RNA sequence for PacBio Tx based on PBid
system("/path/to/file/ISO_PAM/RNAseqFromPBids.bash", intern = FALSE, wait = TRUE)
#Obtain AA sequence for PBid
system("/path/to/file/ISO_PAM/AAseqFromPBids.bash",intern = FALSE,wait = TRUE)
#Remove * from end of AA sequence
system("sed '/*/s/*//' /path/to/file/ISO_PAM/TempAA_fromPBid.txt > /path/to/file/ISO_PAM/TempAA_fromPBidedit.txt", intern = FALSE, wait = TRUE)

#Request known isoform transcriptIDs for ENSG
ENSGinput_biomaRt<- gsub("\\..*", "", ENSGinput)
#ensamble isoform identification
#TxName_EnsemblTxID_TxLength_Strand_peptideID <- getBM(attributes = c("external_transcript_name","ensembl_transcript_id_version","transcript_length","strand","ensembl_peptide_id"), filters = "ensembl_gene_id", values = ENSGinput_biomaRt, mart)
#transcript_appris identifies canonical
TxName_EnsemblTxID_TxLength_Strand_peptideID_appris <- getBM(attributes = c("external_transcript_name","ensembl_transcript_id_version","transcript_length","strand","ensembl_peptide_id","transcript_appris"), filters = "ensembl_gene_id", values = ENSGinput_biomaRt, mart)
#Edit strand column to keep - or + only
Strand_wNumber<-TxName_EnsemblTxID_TxLength_Strand_peptideID_appris$strand
Strandonlyedit<- gsub("1", "", Strand_wNumber)
Strandonly<-as.data.frame(Strandonlyedit)
Strandonly$ensembl_transcript_id_version<-TxName_EnsemblTxID_TxLength_Strand_peptideID_appris$ensembl_transcript_id_version
names(Strandonly)<-c("Strandonly","ensembl_transcript_id_version")

#get RNA sequence from ensembl
ensembl_transcript_id_version_input<-TxName_EnsemblTxID_TxLength_Strand_peptideID_appris$ensembl_transcript_id_version
mRNAseq_EnsemblTxID<-biomaRt::getSequence(id = ensembl_transcript_id_version_input,type="ensembl_transcript_id_version",seqType="cdna",mart= mart)
EnsemblTxID_TxName_TxLength_Strand_peptideID_mRNAseq<-merge(TxName_EnsemblTxID_TxLength_Strand_peptideID_appris,mRNAseq_EnsemblTxID, by='ensembl_transcript_id_version', all=TRUE)
EnsemblTxID_TxName_TxLength_Strand_peptideID_mRNAseq$strand<-NULL
EnsemblTxID_TxName_TxLength_Strand_peptideID_mRNAseq$ensembl_peptide_id<-NULL
EnsemblTxID_TxName_TxLength_mRNAseq_Strandonly<-merge(EnsemblTxID_TxName_TxLength_Strand_peptideID_mRNAseq, Strandonly, by='ensembl_transcript_id_version', all=TRUE)
EnsemblTxID_TxName_TxLength_mRNAseq_StrandonlyOrder<-EnsemblTxID_TxName_TxLength_mRNAseq_Strandonly[c("external_transcript_name","ensembl_transcript_id_version","transcript_length","Strandonly","cdna")]
TranscriptsLebeled_RNAseq<-unite(EnsemblTxID_TxName_TxLength_mRNAseq_StrandonlyOrder, RNA_Tx_Label, c("external_transcript_name","ensembl_transcript_id_version","transcript_length","Strandonly"),sep = " ",remove=TRUE)

#Write TranscriptsLebeled_RNAseq into fasta format txt file
EnsemblTxID_RNAseq_fas<-dataframe2fas(TranscriptsLebeled_RNAseq, file = "/path/to/file/ISO_PAM/TempRNA_ENSGid_RNAseq.txt")
#Append PacBio RNA sequences to Ensembl RNA sequences
system("/path/to/file/ISO_PAM/CombineALL_RNAseq.bash",intern = FALSE,wait = TRUE)

#create complete RNAfasta file
#ALL_RNASequencePREFasta<-readBStringSet("/path/to/file/ISO_PAM/TempRNA_ALL_RNAseq.txt")
#testdf<-strsplit(ALL_RNASequencePREFasta)

#get Protein sequence from Ensembl
Ensembl_ID_PeptideSeq<-biomaRt::getSequence(id =ensembl_transcript_id_version_input,type="ensembl_transcript_id_version",seqType="peptide",mart= mart)
names(Ensembl_ID_PeptideSeq)<-c("AAseq","ensembl_transcript_id_version")

#Count AA in seq
Ensembl_ID_PeptideSeq_countAA<-Ensembl_ID_PeptideSeq
Ensembl_ID_PeptideSeqb<- gsub("Sequence unavailable", "", Ensembl_ID_PeptideSeq_countAA$AAseq)
Ensembl_ID_PeptideSeqc<- gsub("\\*", "", Ensembl_ID_PeptideSeqb)
Length_AAseq<- as.data.frame(nchar(Ensembl_ID_PeptideSeqc))
AAseq_EnsemblTxID_AALength<-cbind(Ensembl_ID_PeptideSeq_countAA,Length_AAseq)
names(AAseq_EnsemblTxID_AALength)<-c("AAseq","ensembl_transcript_id_version","AA_Length")
Ensembl_ID_PeptideSeqClean<-Ensembl_ID_PeptideSeq
Ensembl_ID_PeptideSeqClean$AAseq<-gsub("Sequence unavailable", "", Ensembl_ID_PeptideSeqClean$AAseq)
Ensembl_ID_PeptideSeqClean$AAseq<-gsub("\\*", "", Ensembl_ID_PeptideSeqClean$AAseq)

TxName_EnsemblTxID_TxLength_Strand_peptideID_appris_label_Protein<-TxName_EnsemblTxID_TxLength_Strand_peptideID_appris
TxName_EnsemblTxID_TxLength_Strand_peptideID_appris_label_Protein$strand<-NULL
TxName_EnsemblTxID_TxLength_Strand_peptideID_appris_label_Protein$transcript_length<- NULL
TxName_EnsemblTxID_TxLength_Strand_peptideID_appris_label_Protein_AALength_AAseq<-merge(TxName_EnsemblTxID_TxLength_Strand_peptideID_appris_label_Protein, AAseq_EnsemblTxID_AALength, by="ensembl_transcript_id_version", all=TRUE)
TxName_EnsemblTxID_TxLength_Strand_peptideID_appris_label_Protein_AALength_AAseq_Strandonly<-merge(TxName_EnsemblTxID_TxLength_Strand_peptideID_appris_label_Protein_AALength_AAseq, Strandonly, by='ensembl_transcript_id_version', all=TRUE)
TxName_EnsemblTxID_TxLength_Strand_peptideID_appris_label_Protein_AALength_AAseq_Strandonly$AAseq<- NULL
TxName_EnsemblTxID_TxLength_Strand_peptideID_appris_label_Protein_AALength_AAseq_Strandonly_AAseqClean<-merge(TxName_EnsemblTxID_TxLength_Strand_peptideID_appris_label_Protein_AALength_AAseq_Strandonly, Ensembl_ID_PeptideSeqClean, by='ensembl_transcript_id_version', all=TRUE)

#Filter out isoforms with no known Protein sequence
TxName_EnsemblTxID_TxLength_Strand_peptideID_appris_label_Protein_AALength_AAseq_Strandonly_AAseqCleanFiltered<- TxName_EnsemblTxID_TxLength_Strand_peptideID_appris_label_Protein_AALength_AAseq_Strandonly_AAseqClean[TxName_EnsemblTxID_TxLength_Strand_peptideID_appris_label_Protein_AALength_AAseq_Strandonly_AAseqClean$ensembl_peptide_id !="",]
ProteinLabeled_AAseq<-unite(TxName_EnsemblTxID_TxLength_Strand_peptideID_appris_label_Protein_AALength_AAseq_Strandonly_AAseqCleanFiltered, Protein_Label, c("external_transcript_name","ensembl_transcript_id_version","ensembl_peptide_id","AA_Length","Strandonly"),sep = " ",remove=TRUE)
ProteinLabeled_AAseq<-plyr::arrange(ProteinLabeled_AAseq, plyr::desc(transcript_appris))

ProteinLabeled_AAseq$transcript_appris<-NULL
#create fasta format file for Ensembl labeled AA seq 
EnsemblTxID_AAseq_fas<-dataframe2fas(ProteinLabeled_AAseq, file = "/path/to/file/ISO_PAM/TempAA_knownAAseq.txt")

#RNA Sequence Alingment
#Read generated RNA.fasta file
ALL_RNASequenceFastaPath<-"/path/to/file/ISO_PAM/TempRNA_ALL_RNAseq.txt"
ALL_RNASequences <- readDNAStringSet(ALL_RNASequenceFastaPath)
#ALL_RNASequencesAlignment <- msa(ALL_RNASequences, "ClustalOmega",order="aligned", type="dna")
ALL_RNASequencesAlignment <- muscle::muscle(ALL_RNASequences, seqtype="dna")
#Print RNA sequence alignment in R
#print(ALL_RNASequencesAlignment, show="complete")

#convert MSA to seqinr::alignment
ALL_RNASequencesAlignment_alig <- msaConvert(ALL_RNASequencesAlignment, type="seqinr::alignment")
#Make df for sequence name and sequence data
ALL_RNASequencesAlignment_alig_dfa<-as.data.frame(ALL_RNASequencesAlignment_alig[["nam"]])
names(ALL_RNASequencesAlignment_alig_dfa)<-"RNA_Tx_Label"
ALL_RNASequencesAlignment_alig_dfb<-as.data.frame(ALL_RNASequencesAlignment_alig[["seq"]])
colnames(ALL_RNASequencesAlignment_alig_dfb)<-"cdna"
ALL_RNASequencesAlignment_alig_df<-cbind(ALL_RNASequencesAlignment_alig_dfa,ALL_RNASequencesAlignment_alig_dfb)
#write alignment as fasta
ALL_RNASequencesAlignment_alig_fas<-dataframe2fas(ALL_RNASequencesAlignment_alig_df, file = "/path/to/file/ISO_PAM/ALL_RNASequencesAlignment_alig.fa")

#optional: similarity matrix
#library(seqinr)
#ALL_RNASequencesAlignment_alig_dist <-dist.alignment(ALL_RNASequencesAlignment_alig_fas, "similarity",gap=1)
#ALL_RNASequencesAlignment_alig_dist_df<-as.matrix(ALL_RNASequencesAlignment_alig_dist)[1:<insert number of sequences>, "<insert reference sequence>", drop=FALSE]

##optional:
#library(ape)
#ALL_RNASequencesAlignment_alig_distTree <- nj(ALL_RNASequencesAlignment_alig_dist)
#plot(ALL_RNASequencesAlignment_alig_distTree, main="Phylogenetic Tree of KnownIsoforms and ISOseq")

#optional: RNA sequence alignment
#RNA sequence alignment to pdf settings####
#msaPrettyPrint(myAlignment, output="pdf", showNames="left",showLogo="none", askForOverwrite=FALSE, verbose=FALSE)
#msaPrettyPrint(ALL_RNASequencesAlignment, output="pdf",showNames="left", 
#               shadingMode="similar",showLogo="none", askForOverwrite=FALSE,
#               furtherCode=c("\\defconsensus{.}{lower}{upper}","\\showruler{1}{top}","\\end{texshade}","\\end{document}"))
#ALL_RNAFingerprintAlignment<-ALL_RNASequencesAlignment
#msaPrettyPrint(ALL_RNAFingerprintAlignment, output="pdf",showNames="left", 
#               shadingMode="similar",showLogo="none", askForOverwrite=FALSE,
#               furtherCode=c("\\fingerprint{500}"))


#Protein Sequence Alingment####
#Read prepared file with PacBio Protein sequences
ProteinSequenceFilePath<-"/path/to/file/ISO_PAM/TempAA_fromPBidedit.txt"
PacBio_ProteinSequences <- readAAStringSet(ProteinSequenceFilePath)
PacBio_ProteinSequencesdf<- as.data.frame(PacBio_ProteinSequences)
names(PacBio_ProteinSequencesdf)<-"PacBio_AAseq"
ProteinAlignmentID<-as.data.frame(names(PacBio_ProteinSequences))
names(ProteinAlignmentID)<-"ID"
ProteinAlignmentIDsplit<-cSplit(ProteinAlignmentID,c(names(ProteinAlignmentID)), sep = ":")
ProteinAlignmentInfo<-as.data.frame(ProteinAlignmentIDsplit$ID_01)
names(ProteinAlignmentInfo)<-"gene"
ProteinAlignmentInfo$PBid<-ProteinAlignmentIDsplit$ID_03
#Column ID_11 includes length and string
ORFedit<-as.data.frame(ProteinAlignmentIDsplit$ID_11)
ORFedit<-cSplit(ORFedit,c(names(ORFedit)), sep = "(")
names(ORFedit)<-c('A','B')
ORFedit$B<-sub(")", "", ORFedit$B)
ProteinAlignmentInfo$ORF<-ORFedit$A
LENedit<-as.data.frame(ProteinAlignmentIDsplit$ID_11)
LENsplit<-cSplit(LENedit,c(names(LENedit)), sep = "(")
names(LENsplit)<-c("a","b")
LENsplit$b<- NULL
LENsplit2<-cSplit(LENsplit,c(names(LENsplit)), sep = "-")
names(LENsplit2)<-c("a","b")
#Calculate Length +1
ProteinAlignmentInfo$Length<-LENsplit2$b-LENsplit2$a+1
#String sign
ProteinAlignmentInfo$String<-ORFedit$B
#Remove gene number column
ProteinAlignmentInfo$gene<-NULL
#ProteinAlignmentInfoID<-unite_(ProteinAlignmentInfo, "gene_PBid_ORF_Len_Str", c("gene","PBid","ORF","Length","String"),sep = " ")
ProteinAlignmentInfoID<-unite_(ProteinAlignmentInfo, "PBid_ORF_Len_Str", c("PBid","ORF","Length","String"),sep = " ")
PacBioLabel_AAseq<-cbind(ProteinAlignmentInfoID,PacBio_ProteinSequencesdf)
PacBioLabel_AAseq_fas<-dataframe2fas(PacBioLabel_AAseq, file = "/path/to/file/ISO_PAM/TempAA_knownLabeledAAseq.txt")
#Append PacBioLabelAAseq to Known Isoform AAseq
system("/path/to/file/ISO_PAM/CombineALL_AAseq.bash",intern = FALSE,wait = TRUE)

#load Combined PacBio and AA sequences for alignment
ProteinSequenceFileCompletePath<-"/path/to/file/ISO_PAM/TempAA_ALL_AAseq.txt"
ENSGandPacBioAAseq <- readAAStringSet(ProteinSequenceFileCompletePath)
ALL_ProteinSequencesAlignment <- muscle::muscle(ENSGandPacBioAAseq, seqtype="protein")

#Sequence alignment exported in fasta format
#convert MSA_AA to seqinr::alignment
ALL_ProteinSequencesAlignment_alig <- msaConvert(ALL_ProteinSequencesAlignment, type="seqinr::alignment")
#Make df for sequence name and sequence data
ALL_ProteinSequencesAlignment_alig_dfa<-as.data.frame(ALL_ProteinSequencesAlignment_alig[["nam"]])
names(ALL_ProteinSequencesAlignment_alig_dfa)<-"Peptide_Label"
ALL_ProteinSequencesAlignment_alig_dfb<-as.data.frame(ALL_ProteinSequencesAlignment_alig[["seq"]])
colnames(ALL_ProteinSequencesAlignment_alig_dfb)<-"Peptide"
ALL_ProteinSequencesAlignment_alig_df<-cbind(ALL_ProteinSequencesAlignment_alig_dfa,ALL_ProteinSequencesAlignment_alig_dfb)
#write alignment as fasta
ALL_ProteinSequencesAlignment_alig_fas<-dataframe2fas(ALL_ProteinSequencesAlignment_alig_df, file = "/path/to/file/ISO_PAM/ALL_ProteinSequencesAlignment_alig.fa")

#get "canonical" uniprot ID
TxName_uniprotswissprot_appris <- getBM(attributes = c("external_transcript_name","uniprotswissprot","transcript_appris"), filters = "ensembl_gene_id", values = ENSGinput_biomaRt, mart)
TxName_uniprotswissprot_appris<-plyr::arrange(TxName_uniprotswissprot_appris, plyr::desc(transcript_appris))
UniprotID_Annotation_Query<-TxName_uniprotswissprot_appris[1,2]

#Generate Uniprot Protein Annotations as HMMtopo file .top####
proteins_acc<-UniprotID_Annotation_Query
proteins_acc_url <- gsub(" ", "%2C", proteins_acc)
baseurl <- "https://www.ebi.ac.uk/proteins/api/features?offset=0&size=100&accession="
url <- paste0(baseurl, proteins_acc_url)
prots_feat <- httr::GET(url, httr::accept_json())
Errcode <- httr::status_code(prots_feat)
{if (Errcode == 200) print("Protein Annotations Complete")
  else print(paste("Protein Annotation error occured. Code:", Errcode))}

prots_feat_red <- httr::content(prots_feat)

features_in_lists_of_six<-prots_feat_red
features_total_plot <- NULL
for (i in 1:length(features_in_lists_of_six)) {
  features_temp <- drawProteins::extract_feat_acc(features_in_lists_of_six[[i]])
  features_temp$order <- i
  features_total_plot <- rbind(features_total_plot, features_temp)}

#Make Disulfide bond annotation table
DISULFID_begin<-features_total_plot[(features_total_plot$type=="DISULFID"),c(1,3)]
names(DISULFID_begin)<-c("type","range")
DISULFID_end<-features_total_plot[(features_total_plot$type=="DISULFID"),c(1,4)]
names(DISULFID_end)<-c("type","range")
DISULFID_combined<-rbind(DISULFID_begin,DISULFID_end)
DISULFID_combined$end<-DISULFID_combined$range
DISULFID_combined$Label<-"{'-'}{SS}"
DISULFID_combined$type<-"\\feature{bbottom}{1}"
names(DISULFID_combined)<-c("type","begin","end","Label")

#Make signal peptide annotation table
SIGNAL<-features_total_plot[(features_total_plot$type=="SIGNAL"),c(1,3:4)]
SIGNAL$Label<-"{box[Blue]}{Signal Peptide[Blue]}"
SIGNAL$type<-"\\feature{top}{1}"
names(SIGNAL)<-c("type","begin","end","Label")
SIGNAL_frameblock<-SIGNAL
SIGNAL_frameblock$type<-"\\frameblock{1}"
SIGNAL_frameblock$Label<-"{Blue[1pt]}"

#Make transmembrane annotation table
TRANSMEM<-features_total_plot[(features_total_plot$type=="TRANSMEM"),c(1,3:4)]
TRANSMEM$Label<-"{box[Red]}{Transmembrane[Red]}"
TRANSMEM$type<-"\\feature{top}{1}"
names(TRANSMEM)<-c("type","begin","end","Label")
TRANSMEM_frameblock<-TRANSMEM
TRANSMEM_frameblock$type<-"\\frameblock{1}"
TRANSMEM_frameblock$Label<-"{Red[1pt]}"

#Make extracellular annotation table
EXTRACELL<-features_total_plot[(features_total_plot$description=="Extracellular"),c(1,3:4)]
EXTRACELL$Label<-"{box[Orange]}{Extracellular[Orange]}"
EXTRACELL$type<-"\\feature{ttop}{1}"
names(EXTRACELL)<-c("type","begin","end","Label")

#Make cytoplasmic annotation table
CYTOPL<-features_total_plot[(features_total_plot$description=="Cytoplasmic"),c(1,3:4)]
CYTOPL$Label<-"{,-,}{Cytoplasmic}"
CYTOPL$type<-"\\feature{ttop}{1}"
names(CYTOPL)<-c("type","begin","end","Label")

#Make post translational modification and secondary structure table
KeepFeatures<-c("BINDING","ACT_SITE","MOD_RES","CARBOHYD","STRAND","HELIX","TURN")
selectedFeatureRows <- (features_total_plot$type %in% KeepFeatures)
features_subset_df<-features_total_plot[selectedFeatureRows,c(1:4)]
PTM<-features_subset_df
PTM$Label<-features_subset_df$type
PTM$type<-"\\feature{bottom}{1}"
PTM$Label<-gsub("CARBOHYD","{'-'}{GlcNAc}",PTM$Label)
PTM$Label<-gsub("STRAND","{'-'}{Strand}",PTM$Label)
PTM$Label<-gsub("TURN","{'-'}{Turn}",PTM$Label)
PTM$Label<-gsub("HELIX","{helix}{a-Helix}",PTM$Label)
PTM$Label<-gsub("MOD_RES","{'-'}{pPhos}",PTM$Label)
PTM$Label<-gsub("ACT_SITE","{'-'}{ActiveSite}",PTM$Label)
PTM$Label<-gsub("BINDING","{'-'}{BindingSite}",PTM$Label)
PTM$description<-NULL

#Transmembrane helix
TRANSMEM_helix<-features_total_plot[(features_total_plot$type=="TRANSMEM"),c(1,2,3:4)]
TRANSMEM_helix$description<-gsub("Heli.*","HELIX",TRANSMEM_helix$description)
TRANSMEM_helix$Label<-TRANSMEM_helix$description
TRANSMEM_helix$type<-"\\feature{bottom}{1}"
TRANSMEM_helix$Label<-gsub("HELIX","{helix}{a-Helix}",TRANSMEM_helix$Label)
TRANSMEM_helix$description<-NULL

#Make kinase domain table
PK_domain<-features_total_plot[(features_total_plot$description=="Protein kinase"),c(1,3:4)]
PK_domain$Label<-"{'-'}{Protein kinase}"
PK_domain$type<-"\\feature{bbottom}{1}"
names(PK_domain)<-c("type","begin","end","Label")

#Lipidation (membrane attachment)
LIPID<-features_total_plot[(features_total_plot$type=="LIPID"),c(1,3:4)]
LIPID$Label<-"{box[Green]}{Lipidation[Green]}"
LIPID$type<-"\\feature{tttop}{1}"
names(LIPID)<-c("type","begin","end","Label")



Annotation_Tables_combined<-rbind(DISULFID_combined,SIGNAL,SIGNAL_frameblock,TRANSMEM,TRANSMEM_frameblock,EXTRACELL,CYTOPL,PTM,TRANSMEM_helix,PK_domain,LIPID)
rm(DISULFID_combined,SIGNAL,SIGNAL_frameblock,TRANSMEM,TRANSMEM_frameblock,EXTRACELL,CYTOPL,PTM,TRANSMEM_helix,PK_domain,LIPID,LIPID_frameblock)
Annotation_Tables_combined2<-unite(Annotation_Tables_combined,"FeatureRange",c("begin","end"),sep = "..", remove = TRUE)
Annotation_Tables_combined2$FeatureRange<-gsub("^","{",Annotation_Tables_combined2$FeatureRange)
Annotation_Tables_combined2$FeatureRange<-gsub("$","}",Annotation_Tables_combined2$FeatureRange)
#write annotations file in .top format
write.table(Annotation_Tables_combined2,"/path/to/file/ISO_PAM/MSA_ProteinAnnotations.top",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="")



#Print protein sequence alignment to pdf settings####
msaPrettyPrint(ALL_ProteinSequencesAlignment, output="pdf",file="/path/to/file/ISO_PAM/MSA_ProteinSequenceAln.pdf",showNames="left",
               paperWidth=11, paperHeight=17, margins=c(0.1, 0.3),askForOverwrite=FALSE,
               showLogo="none",showConsensus="none",
               furtherCode=c("\\showruler{top}{1}",
                             "\\shadingmode{diverse}{1}",
                             "\\nomatchresidues{Black}{Gray25}{upper}{up}",
                             "\\similarresidues{Black}{LightYellow}{upper}{up}",
                             "\\pepgroups{GAVLIMFWP,STCYNQ,DE,KRH}",
                             "\\pepsims{G}{AVLIMFWP}",
                             "\\pepsims{A}{GVLIMFWP}",
                             "\\pepsims{V}{GALIMFWP}",
                             "\\pepsims{L}{GAVIMFWP}",
                             "\\pepsims{I}{GAVLMFWP}",
                             "\\pepsims{M}{GAVLIFWP}",
                             "\\pepsims{F}{GAVLIMWP}",
                             "\\pepsims{W}{GAVLIMFP}",
                             "\\pepsims{P}{GAVLIMFW}",
                             "\\pepsims{S}{TCYNQ}",
                             "\\pepsims{T}{SCYNQ}",
                             "\\pepsims{C}{STYNQ}",
                             "\\pepsims{Y}{STCNQ}",
                             "\\pepsims{N}{STCYQ}",
                             "\\pepsims{Q}{STCYN}",
                             "\\pepsims{D}{E}",
                             "\\pepsims{E}{D}",
                             "\\pepsims{K}{RH}",
                             "\\pepsims{R}{KH}",
                             "\\pepsims{H}{KR}",
                             "\\conservedresidues{Black}{ColdHot15}{{.}}{up}",
                             "\\allmatchresidues{Black}{White}{{.}}{up}",
                             "\\gapchar{-}",
                             "\\includeHMMTOP{1}{/path/to/file/ISO_PAM/MSA_ProteinAnnotations.top}",
                             "\\end{texshade}","\\end{document}"))

ALL_ProteinFingerprintAlignment<-ALL_ProteinSequencesAlignment
msaPrettyPrint(ALL_ProteinFingerprintAlignment, output="pdf",file="/path/to/file/ISO_PAM/MSA_ProteinFingerprintAln.pdf",showNames="left", 
               showLogo="none", askForOverwrite=FALSE,
               paperWidth=20, paperHeight=20, margins=c(1, 0.5),
               furtherCode=c("\\showruler{top}{1}",
                             "\\shadingmode{diverse}",
                             "\\nomatchresidues{Black}{Gray25}{upper}{up}",
                             "\\similarresidues{Black}{LightYellow}{upper}{up}",
                             "\\pepgroups{GAVLIMFWP,STCYNQ,DE,KRH}",
                             "\\pepsims{G}{AVLIMFWP}",
                             "\\pepsims{A}{GVLIMFWP}",
                             "\\pepsims{V}{GALIMFWP}",
                             "\\pepsims{L}{GAVIMFWP}",
                             "\\pepsims{I}{GAVLMFWP}",
                             "\\pepsims{M}{GAVLIFWP}",
                             "\\pepsims{F}{GAVLIMWP}",
                             "\\pepsims{W}{GAVLIMFP}",
                             "\\pepsims{P}{GAVLIMFW}",
                             "\\pepsims{S}{TCYNQ}",
                             "\\pepsims{T}{SCYNQ}",
                             "\\pepsims{C}{STYNQ}",
                             "\\pepsims{Y}{STCNQ}",
                             "\\pepsims{N}{STCYQ}",
                             "\\pepsims{Q}{STCYN}",
                             "\\pepsims{D}{E}",
                             "\\pepsims{E}{D}",
                             "\\pepsims{K}{RH}",
                             "\\pepsims{R}{KH}",
                             "\\pepsims{H}{KR}",
                             "\\conservedresidues{Black}{ColdHot15}{{.}}{up}",
                             "\\allmatchresidues{Black}{White}{{.}}{up}",
                             "\\gapchar{-}",
                             "\\includeHMMTOP{1}{/path/to/file/ISO_PAM/MSA_ProteinAnnotations.top}",
                             "\\fingerprint{250}","\\gapchar{rule}","\\end{texshade}","\\end{document}"))


###Reference locations
#orf1ab: 1 - 7096                                                                    
#S: 7097 - 8369                                                                       
#ORF3a: 8370 - 8644                                                                   
#E: 8645 - 8719                                                                       
#M: 8720 - 8941                                                                       
#ORF6: 8942 - 9002                                                                    
#ORF7a: 9003 - 9123                                                                   
#ORF8: 9124 - 9244                                                                    
#N: 9245 - 9663                                                                       
#ORF10: 9664 - 9701


#############Dependencies and working directory############
##########################################################
'%ni%' <- Negate('%in%')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install()

packages <- c("scales", "data.table", "ggrepel", "ggplot2", "viridis", "ggnewscale", "seqinr", 
              "DESeq2", "GenomicRanges", "gplots", "ggbeeswarm", "ggallin", "stringr", "gridExtra",
              "pROC", "caret", "RColorBrewer", "dplyr", "cowplot", "ggpubr", "huxtable", "doMC",
              "officer", "venneuler")

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  BiocManager::install(setdiff(packages, rownames(installed.packages())))  
}

library(scales)
library(data.table)
library(ggrepel)
library(ggplot2)
library(viridis)
library(ggnewscale)
library(seqinr)
library(DESeq2)
library(GenomicRanges)
library(gplots)
library(ggbeeswarm)
library(ggallin)
library(stringr)
library(gridExtra)
library(pROC)
library(caret)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(ggpubr)
library(huxtable)
library(doMC)
library(officer)
library(venneuler)

WORKING_DIR = "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Figures/COVID/"

theme = readRDS("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/ggplot_theme.R")
theme_set(theme)



######Pre-processing#######################
###########################################

####Create fasta from IEDB, all viruses###############

##Create fa to run tools
Alleles_to_keep=c("HLA-DRB1*01:01", "HLA-DRB1*03:01", "HLA-DRB1*04:01", "HLA-DRB1*07:01", "HLA-DRB1*11:01", "HLA-DRB1*13:01", "HLA-DRB1*15:01",
                  "HLA-DQA1*01:02/DQB1*06:02", "HLA-DQA1*05:01/DQB1*02:01", "HLA-DQA1*02:01/DQB1*02:02", "HLA-DQA1*05:05/DQB1*03:01",
                  "HLA-DQA1*01:01/DQB1*05:01", "HLA-DQA1*03:01/DQB1*03:02", "HLA-DQA1*03:03/DQB1*03:01", "HLA-DQA1*01:03/DQB1*06:03")

##Read in IEDB data containing all virus reads for MHC-II prediction
iedb_v = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/mhc_all_virus_combined.csv")
iedb_v = iedb_v[which(iedb_v$allele%in% Alleles_to_keep),]

###Change the length here, 12-20mers
for(q in c(12:20)){

  iedb_v = iedb_v[which(nchar(iedb_v$peptide) == q),]
  
  ##Keeping unique reads and writing them out in fasta format for submission to netMHCIIpan
  peps = unique(iedb_v$peptide)
  id = paste0(">",seq(1,length(peps),1))
  fa =c()
  for(n in 1:length(id)){
    fa = c(fa, id[n])
    fa = c(fa, peps[n])
  }
  
  if(q != 15){  #Max entries for netMHCIIpan is 5000
    fa = as.data.table(fa)
    fwrite(fa, paste0("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_all_viruses_combined_",q,"mer.fa"), col.names = F, row.names = F)
  
  }else{  #Only the case for 15mers
    fa1 = fa[1:4000,]
    fa2 = fa[4001:8000,]
    fa3 = fa[8001:12662,]
    fwrite(fa1, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_all_viruses_combined_15mer_1.fa", col.names = F, row.names = F)
    fwrite(fa2, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_all_viruses_combined_15mer_2.fa", col.names = F, row.names = F)
    fwrite(fa3, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_all_viruses_combined_15mer_3.fa", col.names = F, row.names = F)
  }
}
#########NetMHCpan BA mode########################
I_8mer=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_8mer.xls")
I_9mer=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_9mer.xls")
I_10mer=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_10mer.xls")
I_11mer=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_11mers.xls")
I_B1501=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_B1501_8-11mers.txt")

#HLA-B*15:01 was not originally predicted, later added to see characteristics but ultimately removed as it didn't reach the 5% genetic frequency threshold

HLAI=rbind(I_8mer, I_9mer, I_10mer, I_11mer)
HLAI=HLAI[order(HLAI$Peptide),]
I_B1501 = I_B1501[order(I_B1501$Peptide),]  
HLAI_BA=cbind(HLAI[,1:38], I_B1501[,4:8], HLAI[,39:ncol(HLAI)])

fwrite(HLAI_BA,"/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_all_BA.txt")

###Melting table so each peptide/allele combination is an individual row
Net = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_all_BA.txt")
colnames(Net)[seq(8,88,5)] = colnames(Net)[seq(7,87,5)]
colnames(Net)[seq(7,87,5)] = paste0("nM_", colnames(Net)[seq(7,87,5)])

Net_rank = melt(Net, id.vars = 2, measure.vars = c(seq(8,88,5)))
Net_nM = melt(Net, id.vars = 2, measure.vars = c(seq(7,87,5)))
Net_log = melt(Net, id.vars = 2, measure.vars = c(seq(6,86,5)))
Net_Protein = melt(Net, id.vars = 3, measure.vars = c(seq(8,88,5)))
Net_Start = melt(Net, id.vars = 1, measure.vars = c(seq(8,88,5)))

Net = cbind(Net_Protein$ID, Net_Start$Pos, Net_rank, Net_log[,3], Net_nM[,3])
colnames(Net) = c("Protein", "Start", "Peptide", "Haplotype", "Rank", "1-log50k", "Binding_affinity")
Net$Protein = sapply(strsplit(Net$Protein, "_"),'[',1)
Net$Protein[which(Net$Protein =="orflab")] = "orf1ab"
Net$Start = Net$Start+1

Net$Haplotype = paste0(substr(Net$Haplotype,1,5), "*", substr(Net$Haplotype,6,10))
Net = Net[-which(Net$Haplotype =="HLA-B*15:01" ),]   ###Rounds to >5% by NEON calculations, but in reality is under 5.

fwrite(Net, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_all_BA.txt")


#########NetMHCpan EL mode########################
I_8mer=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_8mer_EL.xls.txt")
I_9mer=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_9mer_EL.xls.txt")
I_10mer=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_10mer_EL.txt")
I_11mer=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_11mer_EL.txt")

HLAI=rbind(I_8mer, I_9mer, I_10mer, I_11mer)
HLAI_EL=HLAI[order(HLAI$Peptide),]

fwrite(HLAI_EL,"/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_all_EL.txt")

Net = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_all_EL.txt")
colnames(Net)[seq(8,88,5)] = colnames(Net)[seq(7,87,5)]
colnames(Net)[seq(7,87,5)] = paste0("nM_", colnames(Net)[seq(7,87,5)])

Net_rank = melt(Net, id.vars = 2, measure.vars = c(seq(8,88,5)))
Net_nM = melt(Net, id.vars = 2, measure.vars = c(seq(7,87,5)))
Net_log = melt(Net, id.vars = 2, measure.vars = c(seq(6,86,5)))
Net_Protein = melt(Net, id.vars = 3, measure.vars = c(seq(8,88,5)))
Net_Start = melt(Net, id.vars = 1, measure.vars = c(seq(8,88,5)))

Net = cbind(Net_Protein$ID, Net_Start$Pos, Net_rank, Net_log[,3],  Net_nM[,3])
colnames(Net) = c("Protein", "Start", "Peptide", "Haplotype", "Rank", "1-log50k", "Binding_affinity")
Net$Protein = sapply(strsplit(Net$Protein, "_"),'[',1)
Net$Protein[which(Net$Protein =="orflab")] = "orf1ab"
Net$Start = Net$Start+1

Net$Haplotype = paste0(substr(Net$Haplotype,1,5), "*", substr(Net$Haplotype,6,10))
Net = Net[-which(Net$Haplotype =="HLA-B*15:01" ),]  ###Rounds to >5% by NEON calculations, but in reality is under 5.

fwrite(Net, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_all_EL.txt")


############MHCflurry########################

Flurry= fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/sars-cov-2.mhcflurry_predictions.csv")
Flurry=cbind(Flurry[,c(1:3,8,9)], NA, Flurry[,c(7,10,11)])
Flurry$sequence_name = sapply(strsplit(Flurry$sequence_name, "\\|"),'[',1)
Flurry = Flurry[-which(duplicated(Flurry)),] ##Removed duplicated B44:03 haplotype
colnames(Flurry) = c("Protein", "Start", "Peptide", "Haplotype", "Rank","1-log50k", "Binding_affinity", "Processing_score", "Presentation_score")
Flurry$Protein[which(Flurry$Protein == "orflab")] = "orf1ab"
Flurry = Flurry[order(Flurry$Peptide),]
Flurry$Start = Flurry$Start+1

Flurry$Haplotype = paste0(substr(Flurry$Haplotype,1,5), "*", substr(Flurry$Haplotype,6,10))

fwrite(Flurry, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCflurry_standardized_all.txt")

############NetMHCIIpan 3.2 -- SARS##################

Netii = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/COVID_human_netMHCIIpan_all.txt")
colnames(Netii)[c(seq(7,25,3), seq(32,54,3))] = substr(colnames(Netii)[c(seq(7,25,3), seq(32,54,3))],6, nchar(colnames(Netii)[c(seq(7,25,3), seq(32,54,3))]))
Netii = Netii[order(Netii$Peptide),]

Netii_Rank = melt(Netii, id.vars = 1, measure.vars = c(seq(7,25,3), seq(32,54,3)))
Netii_nM = melt(Netii, id.vars = 1, measure.vars = c(seq(6,24,3), seq(31,53,3)))
Netii_log = melt(Netii, id.vars = 1, measure.vars = c(seq(5,23,3), seq(30,52,3)))
Netii_protein = melt(Netii, id.vars = 3, measure.vars = c(seq(7,25,3), seq(32,54,3)))
Netii_Start = melt(Netii, id.vars = 2, measure.vars = c(seq(7,25,3), seq(32,54,3)))

Netii = cbind(Netii_protein$ID, Netii_Start$Pos, Netii_Rank, Netii_log[,3],  Netii_nM[,3])
colnames(Netii) = c("Protein", "Start", "Peptide", "Haplotype", "Rank", "1-log50k", "Binding_affinity")
Netii$Protein = sapply(strsplit(Netii$Protein, "_"),'[',1)
Netii$Protein[which(Netii$Protein =="orflab")] = "orf1ab"

Netii$Haplotype = str_replace_all(Netii$Haplotype, pattern = "DRB1_", replacement = "HLA-DRB1*")
Netii$Haplotype = str_replace_all(Netii$Haplotype, pattern = "HLA-DQA1", replacement = "HLA-DQA1*")
Netii$Haplotype = str_replace_all(Netii$Haplotype, pattern = "-DQB1", replacement = "/DQB1*")
Netii$Haplotype = paste0(substr(Netii$Haplotype, 1, nchar(Netii$Haplotype)-2),":", substr(Netii$Haplotype, nchar(Netii$Haplotype)-1, nchar(Netii$Haplotype)))
Netii$Haplotype[which(substr(Netii$Haplotype,1,6) == "HLA-DQ")] = paste0(substr(Netii$Haplotype[which(substr(Netii$Haplotype,1,6) == "HLA-DQ")] ,1,11),":",
                                                                         substr(Netii$Haplotype[which(substr(Netii$Haplotype,1,6) == "HLA-DQ")] ,12,
                                                                                nchar(Netii$Haplotype[which(substr(Netii$Haplotype,1,6) == "HLA-DQ")])))

fwrite(Netii, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_all_BA.txt")


############NetMHCIIpan 4.0 -- SARS##################

Netii = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/COVID_human_netMHCIIpan4.0.xls")
Netii = Netii[order(Netii$Peptide),]

Netii_BA_Rank = melt(Netii, id.vars = 2, measure.vars = seq(9,79,5))
Netii_BA = melt(Netii, id.vars = 2, measure.vars = seq(8,78,5))
Netii_BA_Score = melt(Netii, id.vars = 1, measure.vars = seq(7,77,5))
Netii_EL_Rank = melt(Netii, id.vars = 1, measure.vars = seq(6,76,5))
Netii_EL_Score = melt(Netii, id.vars = 1, measure.vars = seq(5,75,5))

Netii_protein = melt(Netii, id.vars = 3, measure.vars = seq(8,78,5))
Netii_Start = melt(Netii, id.vars = 1, measure.vars = seq(8,78,5))

Netii = as.data.table(cbind(Netii_protein$ID, Netii_Start$Pos, Netii_BA$Peptide, as.character(Netii_BA$variable), Netii_BA$value, Netii_BA_Rank$value, Netii_BA_Score$value, Netii_EL_Rank$value, Netii_EL_Score$value))
colnames(Netii) = c("Protein", "Start", "Peptide", "Haplotype", "Binding_affinity", "BA_Rank", "BA_Score", "EL_Rank", "EL_Score")
Netii$Protein = sapply(strsplit(Netii$Protein, "_"),'[',1)
Netii$Protein[which(Netii$Protein =="orflab")] = "orf1ab"

Netii$Haplotype = str_replace_all(Netii$Haplotype, pattern = "DRB1_", replacement = "HLA-DRB1*")
Netii$Haplotype = str_replace_all(Netii$Haplotype, pattern = "HLA-DQA1", replacement = "HLA-DQA1*")
Netii$Haplotype = str_replace_all(Netii$Haplotype, pattern = "-DQB1", replacement = "/DQB1*")
Netii$Haplotype = paste0(substr(Netii$Haplotype, 1, nchar(Netii$Haplotype)-2),":", substr(Netii$Haplotype, nchar(Netii$Haplotype)-1, nchar(Netii$Haplotype)))
Netii$Haplotype[which(substr(Netii$Haplotype,1,6) == "HLA-DQ")] = paste0(substr(Netii$Haplotype[which(substr(Netii$Haplotype,1,6) == "HLA-DQ")] ,1,11),":",
                                                                         substr(Netii$Haplotype[which(substr(Netii$Haplotype,1,6) == "HLA-DQ")] ,12,
                                                                                nchar(Netii$Haplotype[which(substr(Netii$Haplotype,1,6) == "HLA-DQ")])))



fwrite(Netii, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan4.0_standardized_all_BA.txt")

####netMHCIIpan 3.2 -- all viruses#########
Netii_12=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_ligand_all_viruses_netmhcIIpan3.2_12mer.xls")
Netii_13=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_ligand_all_viruses_netmhcIIpan3.2_13mer.xls")
Netii_14=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_ligand_all_viruses_netmhcIIpan3.2_14mer.xls")
Netii_15_1=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_ligand_all_viruses_netmhcIIpan3.2_15mer_1.xls")
Netii_15_2=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_ligand_all_viruses_netmhcIIpan3.2_15mer_2.xls")
Netii_15_3=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_ligand_all_viruses_netmhcIIpan3.2_15mer_3.xls")
Netii_16=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_ligand_all_viruses_netmhcIIpan3.2_16mer.xls")
Netii_17=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_ligand_all_viruses_netmhcIIpan3.2_17mer.xls")
Netii_18=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_ligand_all_viruses_netmhcIIpan3.2_18mer.xls")
Netii_19=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_ligand_all_viruses_netmhcIIpan3.2_19mer.xls")
Netii_20=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_ligand_all_viruses_netmhcIIpan3.2_20mer.xls")

Netii = rbind(Netii_12,Netii_13,Netii_14,Netii_15_1,Netii_15_2,Netii_15_3,Netii_16,Netii_17,Netii_18,Netii_19,Netii_20)


Netii = Netii[order(Netii$Peptide),]
colnames(Netii)[c(seq(6,33,3))] = colnames(Netii)[c(seq(5,32,3))]

Netii_Rank = melt(Netii, id.vars = 1, measure.vars = c(seq(6,33,3)))
Netii_nM = melt(Netii, id.vars = 1, measure.vars = c(seq(5,32,3)))
Netii_log = melt(Netii, id.vars = 1, measure.vars = c(seq(4,31,3)))
Netii_ID = melt(Netii, id.vars = 3, measure.vars = c(seq(6,33,3)))
Netii_Pep = melt(Netii, id.vars = 2, measure.vars = c(seq(6,33,3)))

Netii_Start = melt(Netii, id.vars = 1, measure.vars = c(seq(6,33,3)))

Netii = cbind(Netii_ID$ID, 1, Netii_Pep$Peptide, Netii_Rank[,2:3], Netii_log[,3],  Netii_nM[,3])
colnames(Netii) = c("ID", "Start", "Peptide", "Haplotype", "Rank", "1-log50k", "Binding_affinity")
#Netii$Protein = sapply(strsplit(Netii$Protein, "_"),'[',1)
#Netii$Protein[which(Netii$Protein =="orflab")] = "orf1ab"


Netii$Haplotype = str_replace_all(Netii$Haplotype, pattern = "DRB1_", replacement = "HLA-DRB1*")
Netii$Haplotype = str_replace_all(Netii$Haplotype, pattern = "HLA-DQA1", replacement = "HLA-DQA1*")
Netii$Haplotype = str_replace_all(Netii$Haplotype, pattern = "-DQB1", replacement = "/DQB1*")
Netii$Haplotype = paste0(substr(Netii$Haplotype, 1, nchar(Netii$Haplotype)-2),":", substr(Netii$Haplotype, nchar(Netii$Haplotype)-1, nchar(Netii$Haplotype)))
Netii$Haplotype[which(substr(Netii$Haplotype,1,6) == "HLA-DQ")] = paste0(substr(Netii$Haplotype[which(substr(Netii$Haplotype,1,6) == "HLA-DQ")] ,1,11),":",
                                                                         substr(Netii$Haplotype[which(substr(Netii$Haplotype,1,6) == "HLA-DQ")] ,12,
                                                                                nchar(Netii$Haplotype[which(substr(Netii$Haplotype,1,6) == "HLA-DQ")])))


fwrite(Netii, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_mhc_all_viruses_all_lengths.txt")


####netMHCIIpan 4.0 -- all viruses#########
Netii_12=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_ligand_all_viruses_netmhcIIpan4.0_12mer.xls")
Netii_13=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_ligand_all_viruses_netmhcIIpan4.0_13mer.xls")
Netii_14=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_ligand_all_viruses_netmhcIIpan4.0_14mer.xls")
Netii_15_1=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_ligand_all_viruses_netmhcIIpan4.0_15mer_1.xls")
Netii_15_2=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_ligand_all_viruses_netmhcIIpan4.0_15mer_2.xls")
Netii_15_3=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_ligand_all_viruses_netmhcIIpan4.0_15mer_3.xls")
Netii_16=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_ligand_all_viruses_netmhcIIpan4.0_16mer.xls")
Netii_17=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_ligand_all_viruses_netmhcIIpan4.0_17mer.xls")
Netii_18=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_ligand_all_viruses_netmhcIIpan4.0_18mer.xls")
Netii_19=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_ligand_all_viruses_netmhcIIpan4.0_19mer.xls")
Netii_20=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_ligand_all_viruses_netmhcIIpan4.0_20mer.xls")

Netii = rbind(Netii_12,Netii_13,Netii_14,Netii_15_1,Netii_15_2,Netii_15_3,Netii_16,Netii_17,Netii_18,Netii_19,Netii_20)

Netii = Netii[order(Netii$Peptide),]

Netii_nM = melt(Netii, id.vars = 1, measure.vars = c(seq(8,53,5)))
Netii_BA_Rank = melt(Netii, id.vars = 1, measure.vars = c(seq(9,54,5)))
Netii_BA_Score = melt(Netii, id.vars = 1, measure.vars = c(seq(7,52,5)))
Netii_EL_Rank = melt(Netii, id.vars = 1, measure.vars = c(seq(6,51,5)))
Netii_EL_Score = melt(Netii, id.vars = 1, measure.vars = c(seq(5,50,5)))

Netii_ID = melt(Netii, id.vars = 3, measure.vars = c(seq(8,53,5)))
Netii_Pep = melt(Netii, id.vars = 2, measure.vars = c(seq(8,53,5)))

Netii_Start = melt(Netii, id.vars = 1, measure.vars = c(seq(8,53,5)))

Netii = data.table(Netii_ID$ID, Netii_nM$Pos, Netii_Pep$Peptide, as.character(Netii_nM$variable), Netii_BA_Rank$value, Netii_BA_Score$value,
                   Netii_nM$value, Netii_EL_Rank$value, Netii_EL_Score$value)
colnames(Netii) = c("ID", "Start", "Peptide", "Haplotype", "BA_Rank", "BA_Score", "Binding_affinity", "EL_Rank", "EL_Score")

Netii$Haplotype = str_replace_all(Netii$Haplotype, pattern = "DRB1_", replacement = "HLA-DRB1*")
Netii$Haplotype = str_replace_all(Netii$Haplotype, pattern = "HLA-DQA1", replacement = "HLA-DQA1*")
Netii$Haplotype = str_replace_all(Netii$Haplotype, pattern = "-DQB1", replacement = "/DQB1*")
Netii$Haplotype = paste0(substr(Netii$Haplotype, 1, nchar(Netii$Haplotype)-2),":", substr(Netii$Haplotype, nchar(Netii$Haplotype)-1, nchar(Netii$Haplotype)))
Netii$Haplotype[which(substr(Netii$Haplotype,1,6) == "HLA-DQ")] = paste0(substr(Netii$Haplotype[which(substr(Netii$Haplotype,1,6) == "HLA-DQ")] ,1,11),":",
                                                                         substr(Netii$Haplotype[which(substr(Netii$Haplotype,1,6) == "HLA-DQ")] ,12,
                                                                                nchar(Netii$Haplotype[which(substr(Netii$Haplotype,1,6) == "HLA-DQ")])))

fwrite(Netii, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan4.0_standardized_mhc_all_viruses_all_lengths.txt")




##################Tetramer GLM analysis################
#######################################################

#############Generate fasta files to run through netMHC(II)pan, MHCflurry

#Read in IEDB tetramer data, filter by human epitopes
Tet = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/tcell_virus_tetramer.csv")
Tet$Length = nchar(Tet$peptide)
Tet = Tet[which(Tet$species =="human"),]

Tet1 = Tet[which(Tet$mhc_class == "I")]
Tet2 = Tet[which(Tet$mhc_class == "II")]

Tet1 = Tet1[which(Tet1$Length %in% c(8:14)),]  #Only 8-14mers were run through netMHCpan

#Format allele names
Tet2$allele[which(substr(Tet2$allele,1,7)=="HLA-DRA")] = sapply(strsplit(Tet2$allele[which(substr(Tet2$allele,1,7)=="HLA-DRA")], "/"), "[", 2)
Tet2$allele[which(substr(Tet2$allele,1,3)=="DRB")] = paste0("HLA-", Tet2$allele[which(substr(Tet2$allele,1,3)=="DRB")])
Tet2 = Tet2[which(substr(Tet2$allele,1,7) == "HLA-DRB")]  ##Filter because DRA is not supported by netMHCIIpan

#Write out fasta file by length.
for(z in unique(Tet1$Length)){
  tet = Tet1[which(Tet1$Length == z),]
  print(paste0(z, ": ", paste0(unique(tet$allele), collapse = ', ')))
  fa = c()
  for (i in 1:nrow(tet)){
    fa = c(fa, paste0(">", tet$allele[i], "_", tet$peptide[i]))
    fa = c(fa, tet$peptide[i])
  }
  fwrite(as.data.table(fa), paste0("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_fa/MHC-I_",z,".fa"), col.names = F, quote = F)
}

for(z in unique(Tet2$Length)){
  
  tet = Tet2[which(Tet2$Length == z),]
  print(paste0(z, ": ", paste0(unique(tet$allele), collapse = ', ')))
  fa = c()
  for (i in 1:nrow(tet)){
    fa = c(fa, paste0(">", tet$allele[i], "_", tet$peptide[i]))
    fa = c(fa, tet$peptide[i])
  }
  fwrite(as.data.table(fa), paste0("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_fa/MHC-II_",z,".fa"), col.names = F, quote = F)
}
###After this step:
#     1. Manually run netMHC(II)pan for generated fasta files
#     2. Format to replace "nM" headers with allele names and remove first row
#     3. Read those files in below

###Read in netMHCpan, MHCflurry data

##BA
MHCI_tet_8mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCI_tetramer_8mer.xls")
MHCI_tet_9mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCI_tetramer_9mer.xls")
MHCI_tet_10mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCI_tetramer_10mer.xls")
MHCI_tet_11mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCI_tetramer_11mer.xls")
MHCI_tet_12mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCI_tetramer_12mer.xls")
MHCI_tet_13mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCI_tetramer_13mer.xls")
MHCI_tet_14mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCI_tetramer_14mer.xls")
MHCI_tet = rbind(MHCI_tet_8mer, MHCI_tet_9mer, MHCI_tet_10mer, MHCI_tet_11mer, MHCI_tet_12mer, MHCI_tet_13mer, MHCI_tet_14mer)

##Formatting data
colnames(MHCI_tet)[seq(8,88,5)] = colnames(MHCI_tet)[seq(7,87,5)]
colnames(MHCI_tet)[seq(7,87,5)] = paste0("nM_", colnames(MHCI_tet)[seq(7,87,5)])

MHCI_tet_rank = melt(MHCI_tet, id.vars = 2, measure.vars = c(seq(8,88,5)))
MHCI_tet_nM = melt(MHCI_tet, id.vars = 2, measure.vars = c(seq(7,87,5)))
MHCI_tet_log = melt(MHCI_tet, id.vars = 2, measure.vars = c(seq(6,86,5)))
MHCI_tet_Protein = melt(MHCI_tet, id.vars = 3, measure.vars = c(seq(8,88,5)))
MHCI_tet_Start = melt(MHCI_tet, id.vars = 1, measure.vars = c(seq(8,88,5)))

MHCI_tet = cbind(MHCI_tet_Protein$ID, MHCI_tet_Start$Pos, MHCI_tet_rank, MHCI_tet_log[,3], MHCI_tet_nM[,3])
colnames(MHCI_tet) = c("Protein", "Start", "Peptide", "Haplotype", "Rank", "1-log50k", "Binding_affinity")
MHCI_tet$Protein = sapply(strsplit(MHCI_tet$Protein, "_"),'[',1)
MHCI_tet$Protein[which(MHCI_tet$Protein =="orflab")] = "orf1ab"
MHCI_tet$Start = MHCI_tet$Start+1

MHCI_tet$Haplotype = paste0(substr(MHCI_tet$Haplotype,1,5), "*", substr(MHCI_tet$Haplotype,6,10))
MHCI_tet = MHCI_tet[-which(MHCI_tet$Haplotype == "HLA-B*15:01"),]
colnames(MHCI_tet)[5] = "BA_Rank"


##EL
MHCI_tet_EL_8mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCI_tetramer_EL_8mer.xls")
MHCI_tet_EL_9mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCI_tetramer_EL_9mer.xls")
MHCI_tet_EL_10mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCI_tetramer_EL_10mer.xls")
MHCI_tet_EL_11mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCI_tetramer_EL_11mer.xls")
MHCI_tet_EL_12mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCI_tetramer_EL_12mer.xls")
MHCI_tet_EL_13mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCI_tetramer_EL_13mer.xls")
MHCI_tet_EL_14mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCI_tetramer_EL_14mer.xls")
MHCI_tet_EL = rbind(MHCI_tet_EL_8mer, MHCI_tet_EL_9mer, MHCI_tet_EL_10mer, MHCI_tet_EL_11mer, MHCI_tet_EL_12mer, MHCI_tet_EL_13mer, MHCI_tet_EL_14mer)

colnames(MHCI_tet_EL)[seq(8,83,5)] = colnames(MHCI_tet_EL)[seq(7,82,5)]
colnames(MHCI_tet_EL)[seq(7,82,5)] = paste0("nM_", colnames(MHCI_tet_EL)[seq(7,82,5)])

MHCI_tet_EL_rank = melt(MHCI_tet_EL, id.vars = 2, measure.vars = c(seq(8,83,5)))
MHCI_tet_EL_nM = melt(MHCI_tet_EL, id.vars = 2, measure.vars = c(seq(7,82,5)))
MHCI_tet_EL_log = melt(MHCI_tet_EL, id.vars = 2, measure.vars = c(seq(6,81,5)))
MHCI_tet_EL_Protein = melt(MHCI_tet_EL, id.vars = 3, measure.vars = c(seq(8,83,5)))
MHCI_tet_EL_Start = melt(MHCI_tet_EL, id.vars = 1, measure.vars = c(seq(8,85,5)))

MHCI_tet_EL = cbind(MHCI_tet_EL_Protein$ID, MHCI_tet_EL_Start$Pos, MHCI_tet_EL_rank, MHCI_tet_EL_log[,3], MHCI_tet_EL_nM[,3])
colnames(MHCI_tet_EL) = c("Protein", "Start", "Peptide", "Haplotype", "Rank", "1-log50k", "Binding_affinity")
MHCI_tet_EL$Protein = sapply(strsplit(MHCI_tet_EL$Protein, "_"),'[',1)
MHCI_tet_EL$Protein[which(MHCI_tet_EL$Protein =="orflab")] = "orf1ab"
MHCI_tet_EL$Start = MHCI_tet_EL$Start+1

MHCI_tet_EL$Haplotype = paste0(substr(MHCI_tet_EL$Haplotype,1,5), "*", substr(MHCI_tet_EL$Haplotype,6,10))
#MHCI_tet_EL = MHCI_tet_EL[-which(MHCI_tet_EL$Haplotype == "HLA-B*15:01"),]
colnames(MHCI_tet_EL)[5] = "EL_Rank"

MHCI_tet = cbind(MHCI_tet, MHCI_tet_EL[,c(5,6)])


###Flurry
#MHCI_tet = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCI_tetramers_all.txt")

Flurry= fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/tcell_virus_tetramer_class1_mhcflurry_predictions.csv")
colnames(Flurry)[c(4,6,7)] = c("Flurry_BA", "Flurry_proc_score", "Flurry_pres_score")

MHCI_tet$Flurry_BA = NA
MHCI_tet$Flurry_proc_score = NA
MHCI_tet$Flurry_pres_score = NA

for(n in 1:nrow(MHCI_tet)){
  print(n)
  if (paste0(MHCI_tet$Haplotype[n], "_", MHCI_tet$Peptide[n]) %in% paste0(Flurry$best_allele, "_", Flurry$peptide)){
    MHCI_tet$Flurry_BA[n] = Flurry$Flurry_BA[which(paste0(Flurry$best_allele, "_", Flurry$peptide) == paste0(MHCI_tet$Haplotype[n], "_", MHCI_tet$Peptide[n]))]
    MHCI_tet$Flurry_proc_score[n] = Flurry$Flurry_proc_score[which(paste0(Flurry$best_allele, "_", Flurry$peptide) == paste0(MHCI_tet$Haplotype[n], "_", MHCI_tet$Peptide[n]))]
    MHCI_tet$Flurry_pres_score[n] = Flurry$Flurry_pres_score[which(paste0(Flurry$best_allele, "_", Flurry$peptide) == paste0(MHCI_tet$Haplotype[n], "_", MHCI_tet$Peptide[n]))]
  }
}


####Add in AA features

Aromatic = c('F','Y','W')
Acidic = c('D','E')
Basic = c('K','R','H')
Small = c('A','G','S','T','P')
Cyclic = 'P'
Thiol = c('C','M')

MHCI_tet$Aromatic = NA
MHCI_tet$Acidic = NA
MHCI_tet$Basic = NA
MHCI_tet$Small = NA
MHCI_tet$Cyclic = NA
MHCI_tet$Thiol = NA

for(n in 1:nrow(MHCI_tet)){
  print(n)
  MHCI_tet$Aromatic[n] = length(which(strsplit(MHCI_tet$Peptide[n], "")[[1]] %in% Aromatic))/nchar(MHCI_tet$Peptide[n])  
  MHCI_tet$Acidic[n] = length(which(strsplit(MHCI_tet$Peptide[n], "")[[1]] %in% Acidic))/nchar(MHCI_tet$Peptide[n])  
  MHCI_tet$Basic[n] = length(which(strsplit(MHCI_tet$Peptide[n], "")[[1]] %in% Basic))/nchar(MHCI_tet$Peptide[n])  
  MHCI_tet$Small[n] = length(which(strsplit(MHCI_tet$Peptide[n], "")[[1]] %in% Small))/nchar(MHCI_tet$Peptide[n])  
  MHCI_tet$Cyclic[n] = length(which(strsplit(MHCI_tet$Peptide[n], "")[[1]] %in% Cyclic))/nchar(MHCI_tet$Peptide[n])  
  MHCI_tet$Thiol[n] = length(which(strsplit(MHCI_tet$Peptide[n], "")[[1]] %in% Thiol))/nchar(MHCI_tet$Peptide[n])  
}

#Read in IEDB tetramer data, filter by human epitopes; reading in again because this step is run after epitope calling
Tet = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/tcell_virus_tetramer.csv")
Tet$Length = nchar(Tet$peptide)
Tet = Tet[which(Tet$species =="human"),]

Tet1 = Tet[which(Tet$mhc_class == "I")]
Tet2 = Tet[which(Tet$mhc_class == "II")]

Tet1 = Tet1[which(Tet1$Length %in% c(8:14)),]  #Only 8-14mers were run through netMHCpan

#Format allele names
Tet2$allele[which(substr(Tet2$allele,1,7)=="HLA-DRA")] = sapply(strsplit(Tet2$allele[which(substr(Tet2$allele,1,7)=="HLA-DRA")], "/"), "[", 2)
Tet2$allele[which(substr(Tet2$allele,1,3)=="DRB")] = paste0("HLA-", Tet2$allele[which(substr(Tet2$allele,1,3)=="DRB")])
Tet2 = Tet2[which(substr(Tet2$allele,1,7) == "HLA-DRB")]

#Add in tetramer data from IEDB set -- for each entry in the computationally predicted set, match the IEDB tetramer results
MHCI_tet$Tetramer = NA
for(n in 1:nrow(MHCI_tet)){
  print(n)
  if (paste0(MHCI_tet$Haplotype[n], "_", MHCI_tet$Peptide[n]) %in% paste0(Tet1$allele, "_", Tet1$peptide)){
    MHCI_tet$Tetramer[n] = Tet1$fraction_entries_positive[which(paste0(Tet1$allele, "_", Tet1$peptide) == paste0(MHCI_tet$Haplotype[n], "_", MHCI_tet$Peptide[n]))]
  }
}

MHCI_tet = MHCI_tet[complete.cases(MHCI_tet),] #Remove entries that weren't found
colnames(MHCI_tet)[6] = "BA_Score"
colnames(MHCI_tet)[9] = "EL_Score"

fwrite(MHCI_tet, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCI_tetramer_all_values.txt")


###GLM with 5-fold CV
MHCI_tet = fread( "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCI_tetramer_all_values.txt")
MHCI_tet$Tetramer = round(MHCI_tet$Tetramer) #Collapse values to binary
MHCI_tet$Tetramer[which(MHCI_tet$Tetramer == 0)] = "False"  #Set to T/F because GLM doesn't handle numerics for binomial
MHCI_tet$Tetramer[which(MHCI_tet$Tetramer == 1)] = "True"

MHCI_tet$Tetramer = factor(MHCI_tet$Tetramer)
MHCI_tet$Haplotype = factor(MHCI_tet$Haplotype)

###Normalize values to range 0-1 so coefficients are easier to interpret
MHCI_tet$Binding_affinity = MHCI_tet$Binding_affinity/50000
MHCI_tet$BA_Rank = MHCI_tet$BA_Rank/100
MHCI_tet$EL_Rank = MHCI_tet$EL_Rank/100
MHCI_tet$Flurry_BA = MHCI_tet$Flurry_BA/50000

#Create balanced folds, keeping haplotype proportions equal
set.seed(1)  #Random seed to keep folds the same each time
folds <- createFolds(factor(MHCI_tet$Haplotype), k = 5, list = FALSE)  #Creates 5x random folds, balancing by haplotype

#spec_cutpoint = 0.9
#cutoff = seq(0,1,.001)
#Score = c()
#mean_auc = c()

train_performance1 = c()  #AUC for each fold's training set
test_performance1 = c()  #AUC for each fold's test set
#FPR_test1 = c()  #False positive rate for each fold

for(z in 1:5){
  train = MHCI_tet[folds != z,]  #Set the training set
  test = MHCI_tet[folds == z,]  #Set the test set
  
  print(table(test$Haplotype))
  
  #Variables included in model.  Ranks were taken out, as was binding affinity score
  vars=c("Binding_affinity","EL_Score","Flurry_BA",
         "Flurry_pres_score","Flurry_proc_score","Aromatic","Acidic","Basic","Small","Cyclic","Thiol")
  
  #Pull relevant columns from "train" to make model formula easier to run.  Column 19 is the tetramer (dependent) variable
  train = data.table(train[,19,with=F], train[,which(colnames(train) %in% vars),with=F])
  
  min.model1 = glm(Tetramer ~ 1, data=train, family = binomial) #Minimum model
  biggest1 <- formula(glm(Tetramer~., data=train, family = binomial)) #Max model, all features
  model_step = step(min.model1, direction='forward', scope=biggest1, trace = F) #Forward stepwise regression
  
  print(summary(model_step))
  
  #Calculate AUC
  ROC=roc(predictor=model_step$fitted.values,train$Tetramer, quiet = T)
  auc = as.numeric(substr(ROC$auc, nchar(ROC$auc) - 21, nchar(ROC$auc)))
  train_performance1 = c(train_performance1, auc) 
  
  # Specificity = c()
  # for (n in cutoff){
  #   Specificity = c(Specificity,  length(which(train$Tetramer[which(model_step$fitted.values<n)] == "False"))/length(which(train$Tetramer == "False")))
  # }
  #Score = c(Score, cutoff[min(which(Specificity>=spec_cutpoint))])
  
  #Predict on test set, keep FPR for test
  pred = predict(model_step, newdata = test, type = "response")
  
  ROC_test=roc(predictor=pred,test$Tetramer, quiet = T)
  auc_test = as.numeric(substr(ROC_test$auc, nchar(ROC_test$auc) - 21, nchar(ROC_test$auc)))
  test_performance1 = c(test_performance1, auc_test) 
  #FPR_test1 = c(FPR_test1, length(which((pred > cutoff[min(which(Specificity>=spec_cutpoint))]) & (test$Tetramer == 'False')))/length(which(test$Tetramer == "False")))
  
}
train_performance1
test_performance1
#FPR_test1


###Now model on all

####Forward stepwise
vars=c("Binding_affinity","EL_Score","Flurry_BA",
       "Flurry_pres_score","Flurry_proc_score","Aromatic","Acidic","Basic","Small","Cyclic","Thiol")

MHCI_mod = data.table(MHCI_tet[,19,with=F], MHCI_tet[,which(colnames(MHCI_tet) %in% vars),with=F])

min.model1 = glm(Tetramer ~ 1, data=MHCI_mod, family = binomial)  #Starting model
biggest1 <- formula(glm(Tetramer~., data=MHCI_mod, family = binomial))  #Max model
fwd.model1 = step(min.model1, direction='forward', scope=biggest1)  #Forward stepwise GLM, going from starting to max model and looking for optimized variables

# #########Older code for when we set cutpoint by specificity#####
# 
# cutoff = seq(0,1,.001)
# 
# FDR = c()
# Sensitivity = c()
# Specificity = c()
# PPV = c()
# 
# for(n in cutoff){
#   print(n)
#   FDR = c(FDR, length(which(MHCI_tet$Tetramer[which(fwd.model1$fitted.values>=n)] == "False"))/length(which(fwd.model1$fitted.values>=n)))
#   Sensitivity = c(Sensitivity,  length(which(MHCI_tet$Tetramer[which(fwd.model1$fitted.values>=n)] == "True"))/length(which(MHCI_tet$Tetramer == "True")))
#   Specificity = c(Specificity,  length(which(MHCI_tet$Tetramer[which(fwd.model1$fitted.values<n)] == "False"))/length(which(MHCI_tet$Tetramer == "False")))
#   PPV = c(PPV, length(which(MHCI_tet$Tetramer[which(fwd.model1$fitted.values>n)] == "True"))/length(which(fwd.model1$fitted.values>n)))
# }
# 
# plot_cutoff = data.table(cutoff, FDR, Sensitivity, Specificity, PPV)
# glm_cutpoint = cutoff[min(which(Specificity >= spec_cutpoint))]
# 
# perf1 = ggplot(data=plot_cutoff, aes(x=cutoff, y=Specificity))+
#   geom_point()+
#   theme(text=element_text(face="bold",size=20,colour="black")) +
#   labs(x="GLM score", y= "Specificity", title = "HLA-I epitope model performance\nViral pMHC tetramers in IEDB")+
#   geom_hline(yintercept = c(spec_cutpoint))+
#   geom_vline(xintercept = glm_cutpoint)+
#   annotate("text", label =paste0(spec_cutpoint),y=spec_cutpoint, x=0, vjust=0, hjust=0, size = 10)+
#   annotate("text", label =paste0(glm_cutpoint),y=0, x=glm_cutpoint, vjust=0, hjust=0, angle = 90, size = 10)
#saveRDS(model_all, paste0("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCI_glm_new.R"))

saveRDS(fwd.model1, paste0("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCI_glm_forward.R"))


####Modeling for MHC-II -- same strategy as MHC-I directly above

####MHCII, read in data fron netMHCIIpan
MHCII_tet_9mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCII_tetramer_9mer.xls")
MHCII_tet_11mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCII_tetramer_11mer.xls")
MHCII_tet_12mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCII_tetramer_12mer.xls")
MHCII_tet_13mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCII_tetramer_13mer.xls")
MHCII_tet_14mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCII_tetramer_14mer.xls")
MHCII_tet_15mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCII_tetramer_15mer.xls")
MHCII_tet_16mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCII_tetramer_16mer.xls")
MHCII_tet_17mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCII_tetramer_17mer.xls")
MHCII_tet_18mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCII_tetramer_18mer.xls")
MHCII_tet_19mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCII_tetramer_19mer.xls")
MHCII_tet_20mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCII_tetramer_20mer.xls")
MHCII_tet_21mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCII_tetramer_21mer.xls")
MHCII_tet = rbind(MHCII_tet_9mer, MHCII_tet_11mer, MHCII_tet_12mer, MHCII_tet_13mer, MHCII_tet_14mer, MHCII_tet_15mer, MHCII_tet_16mer, MHCII_tet_17mer, 
                  MHCII_tet_18mer, MHCII_tet_19mer, MHCII_tet_20mer, MHCII_tet_21mer)

#Formating data
MHCII_tet = MHCII_tet[order(MHCII_tet$Peptide),]

MHCII_tet_nM = melt(MHCII_tet, id.vars = 1, measure.vars = c(seq(8,78,5)))
MHCII_tet_BA_Rank = melt(MHCII_tet, id.vars = 1, measure.vars = c(seq(9,79,5)))
MHCII_tet_BA_Score = melt(MHCII_tet, id.vars = 1, measure.vars = c(seq(7,77,5)))
MHCII_tet_EL_Rank = melt(MHCII_tet, id.vars = 1, measure.vars = c(seq(6,76,5)))
MHCII_tet_EL_Score = melt(MHCII_tet, id.vars = 1, measure.vars = c(seq(5,75,5)))

MHCII_tet_ID = melt(MHCII_tet, id.vars = 3, measure.vars = c(seq(8,78,5)))
MHCII_tet_Pep = melt(MHCII_tet, id.vars = 2, measure.vars = c(seq(8,78,5)))

MHCII_tet_Start = melt(MHCII_tet, id.vars = 1, measure.vars = c(seq(8,78,5)))

MHCII_tet = data.table(MHCII_tet_ID$ID, MHCII_tet_nM$Pos, MHCII_tet_Pep$Peptide, as.character(MHCII_tet_nM$variable), MHCII_tet_BA_Rank$value, MHCII_tet_BA_Score$value,
                       MHCII_tet_nM$value, MHCII_tet_EL_Rank$value, MHCII_tet_EL_Score$value)
colnames(MHCII_tet) = c("ID", "Start", "Peptide", "Haplotype", "BA_Rank", "BA_Score", "Binding_affinity", "EL_Rank", "EL_Score")


MHCII_tet$Haplotype = str_replace_all(MHCII_tet$Haplotype, pattern = "DRB1_", replacement = "HLA-DRB1*")
MHCII_tet$Haplotype = str_replace_all(MHCII_tet$Haplotype, pattern = "HLA-DQA1", replacement = "HLA-DQA1*")
MHCII_tet$Haplotype = str_replace_all(MHCII_tet$Haplotype, pattern = "-DQB1", replacement = "/DQB1*")
MHCII_tet$Haplotype = paste0(substr(MHCII_tet$Haplotype, 1, nchar(MHCII_tet$Haplotype)-2),":", substr(MHCII_tet$Haplotype, nchar(MHCII_tet$Haplotype)-1, nchar(MHCII_tet$Haplotype)))
MHCII_tet$Haplotype[which(substr(MHCII_tet$Haplotype,1,6) == "HLA-DQ")] = paste0(substr(MHCII_tet$Haplotype[which(substr(MHCII_tet$Haplotype,1,6) == "HLA-DQ")] ,1,11),":",
                                                                                 substr(MHCII_tet$Haplotype[which(substr(MHCII_tet$Haplotype,1,6) == "HLA-DQ")] ,12,
                                                                                        nchar(MHCII_tet$Haplotype[which(substr(MHCII_tet$Haplotype,1,6) == "HLA-DQ")])))

####Add in AA features

Aromatic = c('F','Y','W')
Acidic = c('D','E')
Basic = c('K','R','H')
Small = c('A','G','S','T','P')
Cyclic = 'P'
Thiol = c('C','M')


MHCII_tet$Aromatic = NA
MHCII_tet$Acidic = NA
MHCII_tet$Basic = NA
MHCII_tet$Small = NA
MHCII_tet$Cyclic = NA
MHCII_tet$Thiol = NA

for(n in 1:nrow(MHCII_tet)){
  print(n)
  MHCII_tet$Aromatic[n] = length(which(strsplit(MHCII_tet$Peptide[n], "")[[1]] %in% Aromatic))/nchar(MHCII_tet$Peptide[n])  
  MHCII_tet$Acidic[n] = length(which(strsplit(MHCII_tet$Peptide[n], "")[[1]] %in% Acidic))/nchar(MHCII_tet$Peptide[n])  
  MHCII_tet$Basic[n] = length(which(strsplit(MHCII_tet$Peptide[n], "")[[1]] %in% Basic))/nchar(MHCII_tet$Peptide[n])  
  MHCII_tet$Small[n] = length(which(strsplit(MHCII_tet$Peptide[n], "")[[1]] %in% Small))/nchar(MHCII_tet$Peptide[n])  
  MHCII_tet$Cyclic[n] = length(which(strsplit(MHCII_tet$Peptide[n], "")[[1]] %in% Cyclic))/nchar(MHCII_tet$Peptide[n])  
  MHCII_tet$Thiol[n] = length(which(strsplit(MHCII_tet$Peptide[n], "")[[1]] %in% Thiol))/nchar(MHCII_tet$Peptide[n])  
}

#######Add tetramer data############

MHCII_tet$Tetramer = NA
for(n in 1:nrow(MHCII_tet)){
  #print(n)
  if (paste0(MHCII_tet$Haplotype[n], "_", MHCII_tet$Peptide[n]) %in% paste0(Tet2$allele, "_", Tet2$peptide)){
    if(length(unique( Tet2$fraction_entries_positive[which(paste0(Tet2$allele, "_", Tet2$peptide) == paste0(MHCII_tet$Haplotype[n], "_", MHCII_tet$Peptide[n]))])) != 1){
      print(n)
    }
    MHCII_tet$Tetramer[n] = unique(round(Tet2$fraction_entries_positive[which(paste0(Tet2$allele, "_", Tet2$peptide) == paste0(MHCII_tet$Haplotype[n], "_", MHCII_tet$Peptide[n]))]))
  }
}

MHCII_tet = MHCII_tet[complete.cases(MHCII_tet),]

fwrite(MHCII_tet, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCII_tetramer_all_values.txt")


###GLM with 5-fold CV, code is identical to HLA-I
MHCII_tet = fread( "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCII_tetramer_all_values.txt")

MHCII_tet$Tetramer = round(MHCII_tet$Tetramer)
MHCII_tet$Tetramer[which(MHCII_tet$Tetramer == 0)] = "False"
MHCII_tet$Tetramer[which(MHCII_tet$Tetramer == 1)] = "True"

MHCII_tet$Tetramer = factor(MHCII_tet$Tetramer)
MHCII_tet$Haplotype = factor(MHCII_tet$Haplotype)

###Try normalizing
MHCII_tet$Binding_affinity = MHCII_tet$Binding_affinity/50000
MHCII_tet$BA_Rank = MHCII_tet$BA_Rank/100
MHCII_tet$EL_Rank = MHCII_tet$EL_Rank/100

set.seed(12)
folds <- createFolds(factor(MHCII_tet$Haplotype), k = 5, list = FALSE)

#spec_cutpoint = 0.9
#cutoff = seq(0,1,.001)
#mean_auc = c()
#Score = c()

train_performance2 = c()  #AUC for training set
test_performance2 = c() #AUC for test set
#FPR_test2 = c()

for(z in 1:5){
  train = MHCII_tet[folds != z,]
  test = MHCII_tet[folds == z,]
  
  print(table(test$Haplotype))
  
  vars=c("Binding_affinity","EL_Score","Aromatic","Acidic","Basic","Small","Cyclic","Thiol")
  
  train = data.table(train[,16,with=F], train[,which(colnames(train) %in% vars),with=F])
  
  min.model = glm(Tetramer ~ 1, data=train, family = binomial)
  biggest <- formula(glm(Tetramer~., data=train, family = binomial))
  model_step = step(min.model, direction='forward', scope=biggest, trace = F)
  
  print(summary(model_step))
  
  ROC=roc(predictor=model_step$fitted.values,train$Tetramer, quiet = T)
  auc = as.numeric(substr(ROC$auc, nchar(ROC$auc) - 21, nchar(ROC$auc)))
  
  train_performance2 = c(train_performance2, auc) 
  
  # Specificity = c()
  # for (n in cutoff){
  #   Specificity = c(Specificity,  length(which(train$Tetramer[which(model_step$fitted.values<n)] == "False"))/length(which(train$Tetramer == "False")))
  # }
  # Score = c(Score, cutoff[min(which(Specificity>=spec_cutpoint))])
  
  pred = predict(model_step, newdata = test, type = "response")

  ROC_test=roc(predictor=pred,test$Tetramer, quiet = T)
  auc_test = as.numeric(substr(ROC_test$auc, nchar(ROC_test$auc) - 21, nchar(ROC_test$auc)))
  test_performance2 = c(test_performance2, auc_test) 
  #FPR_test2 = c(FPR_test2, length(which((pred > cutoff[min(which(Specificity>=spec_cutpoint))]) & (test$Tetramer == 'False')))/length(which(test$Tetramer == "False")))
}

train_performance2
test_performance2
#FPR_test2


######Generating full model

####Forward stepwise
vars=c("Binding_affinity","EL_Score","Aromatic","Acidic","Basic","Small","Cyclic","Thiol")

MHCII_mod = data.table(MHCII_tet[,16,with=F], MHCII_tet[,which(colnames(MHCII_tet) %in% vars),with=F])

min.model = glm(Tetramer ~ 1, data=MHCII_mod, family = binomial)
biggest <- formula(glm(Tetramer~., data=MHCII_mod, family = binomial))
fwd.model2 = step(min.model, direction='forward', scope=biggest)

# ###########
# cutoff = seq(0,1,.001)
# 
# FDR = c()
# Sensitivity = c()
# Specificity = c()
# PPV=c()
# for(n in cutoff){
#   print(n)
#   FDR = c(FDR, length(which(MHCII_tet$Tetramer[which(fwd.model2$fitted.values>=n)] == "False"))/length(which(fwd.model2$fitted.values>=n)))
#   Sensitivity = c(Sensitivity,  length(which(MHCII_tet$Tetramer[which(fwd.model2$fitted.values>=n)] == "True"))/length(which(MHCII_tet$Tetramer == "True")))
#   Specificity = c(Specificity,  length(which(MHCII_tet$Tetramer[which(fwd.model2$fitted.values<n)] == "False"))/length(which(MHCII_tet$Tetramer == "False")))
#   PPV = c(PPV, length(which(MHCI_tet$Tetramer[which(fwd.model2$fitted.values>n)] == "True"))/length(which(fwd.model2$fitted.values>n)))
#   
# }
# 
# plot_cutoff = data.table(cutoff, FDR, Sensitivity, Specificity,PPV)
# glm_cutpoint = cutoff[min(which(Specificity >= spec_cutpoint))]
# 
# perf2 = ggplot(data=plot_cutoff, aes(x=cutoff, y=Specificity))+
#   geom_point()+
#   theme(text=element_text(face="bold",size=20,colour="black")) +
#   labs(x="GLM score", y= "Specificity", title ="HLA-II epitope model performance\nViral pMHC tetramers in IEDB")+
#   geom_hline(yintercept = c(spec_cutpoint))+
#   geom_vline(xintercept = glm_cutpoint)+
#   annotate("text", label =paste0(spec_cutpoint),y=spec_cutpoint, x=0, vjust=0, hjust=0, size = 10)+
#   annotate("text", label =paste0(glm_cutpoint),y=0, x=glm_cutpoint, vjust=0, hjust=0, angle = 90, size = 10)
# 
# #saveRDS(model_all2, paste0("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCII_glm_new.R"))

saveRDS(fwd.model2, paste0("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCII_glm_forward.R"))


#######################################################################################
########Calculate GLM scores for SARS-CoV-2 predicted epitopes#########################

#Models generated above
MHCI_glm = readRDS("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCI_glm_forward.R")
MHCII_glm = readRDS("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCII_glm_forward.R")

#HLA-I data
Net_BA =  fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_BA.txt")
Net_EL= fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_all_EL.txt")
Flurry = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCflurry_standardized_all.txt")

#Order by haplotype and peptide so I can combine sets easily
Net_BA = Net_BA[order(Net_BA$Haplotype),]
Net_BA = Net_BA[order(Net_BA$Peptide),]
Net_EL = Net_EL[order(Net_EL$Haplotype),]
Net_EL = Net_EL[order(Net_EL$Peptide),]
Flurry = Flurry[order(Flurry$Haplotype),]
Flurry = Flurry[order(Flurry$Peptide),]

Net_combined = cbind(Net_BA, Net_EL[,c(5:6)], Flurry[,c(5,7:9)]) #Join sets together, normalize values for GLM 
colnames(Net_combined)[5:13] = c("BA_Rank", "BA_Score", "Binding_affinity", "EL_Rank", "EL_Score", "Flurry_Rank", "Flurry_BA", "Flurry_proc_score", "Flurry_pres_score")
Net_combined$BA_Rank = Net_combined$BA_Rank/100
Net_combined$Binding_affinity = Net_combined$Binding_affinity/50000
Net_combined$EL_Rank = Net_combined$EL_Rank/100
Net_combined$Flurry_Rank = Net_combined$Flurry_Rank/100
Net_combined$Flurry_BA = Net_combined$Flurry_BA/50000

###Add peptide info
features =list(Aromatic = c('F','Y','W'),
               Acidic = c('D','E'),
               Basic = c('K','R','H'),
               Small = c('A','G','S','T','P'),
               Cyclic = c('P'),
               Thiol = c('C','M'))

feature_tab = data.table()  #Table for the amino acid features, to be appended to the rest of the data

for(f in features){  #For each amino acid feature set
  print(f)
  temp_list = data.table()
  for(n in 1:length(f)){  #For each amino acid of the feature set
    print(n)
    temp_list = cbind(temp_list, c(str_count(Net_combined$Peptide, pattern = f[n])))  #Count number of occurances of that AA across all peptides; used stringr because it allows for list-based counting, which is dramatically faster
  }
  hits = rowSums(temp_list)/nchar(Net_combined$Peptide) #Normalize values to percentage
  feature_tab = cbind(feature_tab, hits)
}
colnames(feature_tab) = c("Aromatic", "Acidic", "Basic", "Small", "Cyclic", "Thiol")

Net_combined = cbind(Net_combined, feature_tab)

#####Run model, write out

Net_combined$GLM_score = predict(MHCI_glm, newdata = Net_combined, type = "response")

#Initially filtered by GLM_score, now by median GLM score
MHCI_tet_median = median(Net_combined$GLM_score[which((Net_combined$Binding_affinity*50000) < 393.4)])
Net_combined$Predicted_tetramer = ifelse(Net_combined$GLM_score > MHCI_tet_median, "True", "False")

fwrite(Net_combined, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAI_all_glm.txt")


##########MHCII, same strategy as above############
Netii_3 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_all_BA.txt")
Netii_4 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan4.0_standardized_all_BA.txt")

Netii_3 = Netii_3[order(Netii_3$Haplotype),]
Netii_3 = Netii_3[order(Netii_3$Peptide),]
Netii_4 = Netii_4[order(Netii_4$Haplotype),]
Netii_4 = Netii_4[order(Netii_4$Peptide),]

Netii_all = cbind(Netii_3, Netii_4[,c(5:9)])
colnames(Netii_all)[5:7] = c("v3.2_Rank", "v3.2_Score", "v3.2_Binding_affinity")

Netii_all$BA_Rank = Netii_all$BA_Rank/100
Netii_all$Binding_affinity = Netii_all$Binding_affinity/50000
Netii_all$EL_Rank = Netii_all$EL_Rank/100


###Add peptide info
features =list(Aromatic = c('F','Y','W'),
               Acidic = c('D','E'),
               Basic = c('K','R','H'),
               Small = c('A','G','S','T','P'),
               Cyclic = c('P'),
               Thiol = c('C','M'))

feature_tab = data.table()

for(f in features){
  print(f)
  temp_list = data.table()
  for(n in 1:length(f)){
    print(n)
    temp_list = cbind(temp_list, c(str_count(Netii_all$Peptide, pattern = f[n])))
  }
  hits = rowSums(temp_list)/nchar(Netii_all$Peptide)
  feature_tab = cbind(feature_tab, hits)
}
colnames(feature_tab) = c("Aromatic", "Acidic", "Basic", "Small", "Cyclic", "Thiol")

Netii_all = cbind(Netii_all, feature_tab)

#####Run model, write out

Netii_all$GLM_score = predict(MHCII_glm, newdata = Netii_all, type = "response")

#Initially filtered by GLM_score cutoff of .5 -- can ignore, now filtered by median
MHCII_tet_median = median(Netii_all$GLM_score[which((Netii_all$v3.2_Binding_affinity) < 220)])
Netii_all$Predicted_tetramer = ifelse(Netii_all$GLM_score > MHCII_tet_median, "True", "False")

fwrite(Netii_all, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAII_all_glm.txt")



###################Figure S2 - Model performance##################
##################################################################

###S2A/B (old)
#grid.arrange(perf1, perf2, nrow=1, ncol=2)

#S2 A/B, GLM score histogram, split by tetramer +/-
plota = data.table(MHCI_tet$Tetramer, MHCI_glm$fitted.values)
ha=ggplot(data=plota)+
  scale_x_continuous(limits = c(0,1))+
  geom_histogram(aes(x = V2, fill = V1), alpha=.5,bins = 50, position="dodge")+
  scale_fill_manual(values = c("blue","red"), name = "Tetramer\npositive")+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x= "GLM score", y="Count", title = "IEDB HLA-I")


plotb = data.table(MHCII_tet$Tetramer, MHCII_glm$fitted.values)
hb=ggplot(data=plotb)+
  scale_x_continuous(limits = c(0,1))+
  geom_histogram(aes(x = V2, fill = V1), alpha=.5,bins = 50, position="dodge")+
  scale_fill_manual(values = c("blue","red"), name = "Tetramer\npositive")+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x= "GLM score", y="Count", title = "IEDB HLA-II")

grid.arrange(ha, hb, nrow=1, ncol=2)


#S2 C/D, Binding affinity histogram, split by tet+/-
plotc = data.table(MHCI_tet$Tetramer, MHCI_tet$BA_Score)
hc=ggplot(data=plotc)+
  geom_histogram(aes(x = V2, fill = V1,y=..ncount..), alpha=.5,bins = 50, position="dodge")+
  scale_fill_manual(values = c("blue","red"), name = "Tetramer\npositive")+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x= "Binding affinity 1-log(50k)", y="Count", title = "IEDB HLA-I")


plotd = data.table(MHCII_tet$Tetramer, MHCII_tet$BA_Score)
hd=ggplot(data=plotd)+
  geom_histogram(aes(x = V2, fill = V1, y=..ncount..), alpha=.5,bins = 50, position="dodge")+
  scale_fill_manual(values = c("blue","red"), name = "Tetramer\npositive")+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x= "Binding affinity 1-log(50k)", y="Count", title = "IEDB HLA-II")

grid.arrange(hc, hd, nrow=1, ncol=2)


#S2 E, AUC for training data
auc_mat_train = data.table(c(train_performance1, train_performance2), c(rep("HLA-I", 5), rep("HLA-II", 5)))
S2E=ggplot(data = auc_mat_train)+
  geom_boxplot(aes(y=V1, x=V2, color = V2), show.legend = F)+
  geom_jitter(aes(y=V1, x=V2, color = V2), show.legend = F)+
  scale_y_continuous(limits = c(.5,1))+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x="", y="AUC", title = "Training set peformance")

#S2F, AUC for test data
auc_mat_test = data.table(c(test_performance1, test_performance2), c(rep("HLA-I", 5), rep("HLA-II", 5)))
S2F=ggplot(data = auc_mat_test)+
  geom_boxplot(aes(y=V1, x=V2, color = V2), show.legend = F)+
  geom_jitter(aes(y=V1, x=V2, color = V2), show.legend = F)+
  scale_y_continuous(limits = c(.5,1))+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x="", y="AUC", title = "Test set performance")
grid.arrange(S2E, S2F, nrow = 1, ncol = 2)


###Figure S4 A/B, histogram of GLM scores for SARS-CoV-2 predicted epitopes, before and after BA filter###
##########################################################################################################

Net_combined = fread( "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAI_all_glm.txt")
Netii_all = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAII_all_glm.txt")

SARS_comb = cbind(data.table(c(rep("HLA-I", nrow(Net_combined)), rep("HLA-II", nrow(Netii_all)))),
                  data.table(c(Net_combined$GLM_score, Netii_all$GLM_score)))
colnames(SARS_comb) = c("HLA", "GLM_score")
ggplot(data=SARS_comb)+
  geom_histogram(aes(x = GLM_score, fill = HLA), alpha=.5,bins = 50, position = "dodge")+
  scale_fill_manual(values = c("blue","red"), name = "")+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x= "GLM score", y="Count", title = "SARS-CoV-2 predicted epitopes\nAll peptide/MHC")


##Post nM filter

##Filter by tetramer models
Net_combined_filt = Net_combined[which((Net_combined$Binding_affinity*50000)<393.4),]
Netii_all_filt = Netii_all[which((Netii_all$v3.2_Binding_affinity)<220),]

MHCI_tet_median = round(median(Net_combined_filt$GLM_score[which((Net_combined_filt$Binding_affinity) < 393.4)]), digits = 3)
MHCII_tet_median = round(median(Netii_all_filt$GLM_score[which((Netii_all_filt$v3.2_Binding_affinity) < 220)]), digits = 3)


SARS_comb_filt = cbind(data.table(c(rep("HLA-I", nrow(Net_combined_filt)), rep("HLA-II", nrow(Netii_all_filt)))),
                       data.table(c(Net_combined_filt$GLM_score, Netii_all_filt$GLM_score)))
colnames(SARS_comb_filt) = c("HLA", "GLM_score")
ggplot(data=SARS_comb_filt)+
  geom_histogram(aes(x = GLM_score, fill = HLA), alpha=.5,bins = 50, position = "dodge")+
  geom_vline(xintercept = MHCI_tet_median, color = "blue")+
  annotate("text", label= paste0("HLA-I median: ", MHCI_tet_median),
           x=MHCI_tet_median, y=600, vjust=0, hjust=0, angle = 90, size = 7.5, color = "blue")+
  geom_vline(xintercept = MHCII_tet_median, color = "red")+
  annotate("text", label = paste0("HLA-II median: ", MHCII_tet_median),
           x=MHCII_tet_median, y=600, vjust=0, hjust=0, angle = 90, size = 7.5, color = "red")+
  scale_fill_manual(values = c("blue","red"), name = "")+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x= "GLM score", y="Count", title = "SARS-CoV-2 predicted epitopes\nPost binding affinity filter")



###########Figure S6 -- MS counts ##################
####################################################
#https://www.biorxiv.org/content/10.1101/2020.03.22.002204v1
Protein = c("M", "N", "ORF1ab", "ORF3a", "ORF6", "ORF7a", "ORF8", "S")
PSM_counts = c(250,4152,1816,119,3,23,33,1984)
Length_normalized = c(1.126126126, 9.909307876, 0.255918828, 0.432727273,0.049180328, 0.190082645,0.270491803,1.558523174)

MS_tab = data.table(Protein, PSM_counts, Length_normalized)
MS_tab = MS_tab[rev(order(as.numeric(MS_tab$Length_normalized))),]
MS_tab$Protein = factor(MS_tab$Protein, levels = MS_tab$Protein)

ggplot(data = MS_tab)+
  geom_bar(aes(x = Protein, y = Length_normalized, fill = Protein), stat = "identity")+
  labs(y = "Length normalized PSM")



# ####Check for missing alleles -- can ignore this, was just checking Alex's work
# Netii_all = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAII_all_glm.txt")
# 
# missing = c("HLA-DQA1*03:01/DQB1*03:02", "HLA-DQA1*01:01/DQB1*05:01", "HLA-DQA1*02:01/DQB1*02:02", "HLA-DQA1*01:03/DQB1*06:03")
# 
# boxplot(Netii_all$v3.2_Binding_affinity[which(Netii_all$Haplotype %in% missing)], Netii_all$v3.2_Binding_affinity[which(Netii_all$Haplotype %ni% missing)])
# t.test(Netii_all$v3.2_Binding_affinity[which(Netii_all$Haplotype %in% missing)], Netii_all$v3.2_Binding_affinity[which(Netii_all$Haplotype %ni% missing)])
# 
# 
# boxplot(Netii_all$GLM_score[which(Netii_all$Haplotype %in% missing)], Netii_all$GLM_score[which(Netii_all$Haplotype %ni% missing)])
# t.test(Netii_all$GLM_score[which(Netii_all$Haplotype %in% missing)], Netii_all$GLM_score[which(Netii_all$Haplotype %ni% missing)])

# #Check Alex's outputs
# 
# Window_27mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Window_27mer.txt")
# b_d_27mer_CD4_8 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/selected-tcell-cd4-cd8-h2b-h2d-27mer.csv")
# 
# W_27mer = Window_27mer[which(Window_27mer$Sequence %in% b_d_27mer_CD4_8$Sequence),]
# 
# unique_I = unique(unlist(strsplit(b_d_27mer_CD4_8$`HLA-I_haplotypes`, ',')))
# unique_II = unique(unlist(strsplit(b_d_27mer_CD4_8$`HLA-II_haplotypes`, ',')))
# 
# #Read in population frequencies of each allele, format allele names to be consistent
# Freq=fread(paste0(WORKING_DIR, "HLA_freq.txt"))
# Freq$Haplotype = str_replace_all(Freq$Haplotype, "DRB1_", "HLA-DRB1*")
# Freq$Haplotype = str_replace_all(Freq$Haplotype, "DQA1", "DQA1*")
# Freq$Haplotype = str_replace_all(Freq$Haplotype, "-DQB1", "/DQB1*")
# Freq$Haplotype = str_replace_all(Freq$Haplotype, "HLA-A", "HLA-A*")
# Freq$Haplotype = str_replace_all(Freq$Haplotype, "HLA-B", "HLA-B*")
# Freq$Haplotype = str_replace_all(Freq$Haplotype, "HLA-C", "HLA-C*")
# Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DR")] = paste0(substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DR")],1,11), ":",
#                                                                       substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DR")],12,13))
# Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DQ")] = paste0(substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DQ")],1,11), ":",
#                                                                       substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DQ")],12,21),":",
#                                                                       substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DQ")],22,23))
# 
# #For DQ reads, values were derived from a study which was biased for caucasian americans.  
# #Unfortunately, paired DQ frequency data was not available for a balanced US population.
# Freq$US[which(is.na(Freq$US))] = Freq$Cauc_Am[which(is.na(Freq$US))]
# 
# 1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% unique_I)]/100,2)))
# 1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% unique_II)]/100,2)))





#############Figure S1A##########################
#################################################

cutpoint = 500 #What nM is considered a binder/non-binder in IEDB data

######Compared performance vs IEDB##############

#Read in IEDB MHC-I data for SARS binding affinity 
iedb = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/mhc_sars2_purified_competitive_radioactivity.csv")
iedb = iedb[,c(1:3)]

#Keep only reads <500nM
iedb_low = iedb[which(iedb$affinity < cutpoint),]
iedb_low = iedb_low[order(iedb_low$peptide),]

#Read in class I calls, filter by those present in IEDB set
Flurry=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCflurry_standardized_all.txt")
Net_EL=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_all_EL.txt")
Net_BA=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_all_BA.txt")

iedb_unique = paste(iedb_low$peptide, iedb_low$allele, sep = "_")
Net_BA_unique = paste(Net_BA$Peptide, Net_BA$Haplotype, sep = "_")
Net_EL_unique = paste(Net_EL$Peptide, Net_EL$Haplotype, sep = "_")
Flurry_unique = paste(Flurry$Peptide, Flurry$Haplotype, sep = "_")

Net_BA = Net_BA[which(Net_BA_unique %in% iedb_unique),]
Net_EL = Net_EL[which(Net_EL_unique %in% iedb_unique),]
Flurry = Flurry[which(Flurry_unique %in% iedb_unique),]

all = paste(Net_BA$Peptide, Net_BA$Haplotype, sep = "_")
iedb_low = iedb_low[which(iedb_unique %in% all),]

#Add in prediction tool features into IEDB matrix
iedb_low$NetMHCpan_BA_Rank=NA
iedb_low$NetMHCpan_BA_log=NA
iedb_low$NetMHCpan_BA_BA=NA

iedb_low$NetMHCpan_EL_Rank=NA
iedb_low$NetMHCpan_EL_log=NA

iedb_low$MHCFlurry_Rank=NA
iedb_low$MHCFlurry_BA=NA
iedb_low$MHCFlurry_Proc=NA
iedb_low$MHCFlurry_ProS=NA


for(n in 1:nrow(iedb_low)){
  
  iedb_low$NetMHCpan_BA_Rank[n] = Net_BA$Rank[which((Net_BA$Peptide == iedb_low$peptide[n])&(Net_BA$Haplotype == iedb_low$allele[n]))]
  iedb_low$NetMHCpan_BA_log[n] = Net_BA$`1-log50k`[which((Net_BA$Peptide == iedb_low$peptide[n])&(Net_BA$Haplotype == iedb_low$allele[n]))]
  iedb_low$NetMHCpan_BA_BA[n] = Net_BA$Binding_affinity[which((Net_BA$Peptide == iedb_low$peptide[n])&(Net_BA$Haplotype == iedb_low$allele[n]))]
  
  iedb_low$NetMHCpan_EL_Rank[n] = Net_EL$Rank[which((Net_EL$Peptide == iedb_low$peptide[n])&(Net_EL$Haplotype == iedb_low$allele[n]))]
  iedb_low$NetMHCpan_EL_log[n] = Net_EL$`1-log50k`[which((Net_EL$Peptide == iedb_low$peptide[n])&(Net_EL$Haplotype == iedb_low$allele[n]))]
  
  iedb_low$MHCFlurry_Rank[n] = Flurry$Rank[which((Flurry$Peptide == iedb_low$peptide[n])&(Flurry$Haplotype == iedb_low$allele[n]))]
  iedb_low$MHCFlurry_BA[n] = Flurry$Binding_affinity[which((Flurry$Peptide == iedb_low$peptide[n])&(Flurry$Haplotype == iedb_low$allele[n]))]
  iedb_low$MHCFlurry_Proc[n] = Flurry$Processing_score[which((Flurry$Peptide == iedb_low$peptide[n])&(Flurry$Haplotype == iedb_low$allele[n]))]
  iedb_low$MHCFlurry_ProS[n] = Flurry$Presentation_score[which((Flurry$Peptide == iedb_low$peptide[n])&(Flurry$Haplotype == iedb_low$allele[n]))]
  
}


##Melt so each prediction feature is separated by row
iedb_melt = melt(iedb_low, id.vars = c(1:3), measure.vars = c(4:12))
colnames(iedb_melt)[4:5] = c("Program", "Rank")
iedb_melt = iedb_melt[complete.cases(iedb_melt),]

#iedb_melt = iedb_melt[which(iedb_melt$Program %in% c('NetMHCpan_BA_BA',  'NetMHCpan_BA_Rank', "NetMHCpan_EL_Rank", 'MHCFlurry_BA', 'MHCFlurry_Rank', "MHCFlurry_Proc", "MHCFlurry_ProS"))]
#iedb_melt = iedb_melt[which(iedb_melt$Program %in% c('NetMHCpan_BA_BA', 'MHCFlurry_BA'))]

fwrite(iedb_melt, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/IEDB_SARS2_HLAI_filtered_competitive_radioactivity.txt") 

#Generating individual plots showing correlation between prediction feature and IEDB nM, with Spearman correlation shown
for(z in unique(iedb_melt$Program)){
  print(z)
  if(z %in% unique(iedb_melt$Program)[1:3]){
    col = viridis(4)[1]
  }else if(z %in% unique(iedb_melt$Program)[4:5]){
    col = viridis(4)[2]
  }else if(z %in% unique(iedb_melt$Program)[6:9]){
    col =  viridis(4)[3]
  }
  
  sub_melt = iedb_melt[which(iedb_melt$Program ==z),]
  test = cor.test(sub_melt$affinity, sub_melt$Rank, method = "spearman") #Will throw "cannot compute exact p-value" warning, this is expected. Will still give good estimate.
  
  p= ggplot(data = sub_melt)+
    geom_point(aes(x = affinity, y= Rank), color = col, alpha=.6, show.legend = F)+
    scale_y_log10()+
    scale_x_log10(breaks = c(.5,5,50,500))+
    ggnewscale::new_scale("color")+
    labs(x="", y="", title = z) +
    annotate("text", label = paste0("Rho: ",round(test$estimate, digits = 3),"\n","p<0.001"),
             x = 0.5*(min(sub_melt$affinity) + max(sub_melt$affinity)), y = min(sub_melt$Rank), vjust = 0, hjust=1)
  
  assign(paste0(z),p + geom_smooth(data = ggplot_build(p)$data[[1]], 
                                   mapping = aes(x=10^x, y= 10^y), color = "black", method = "lm",se=FALSE, show.legend = F) )
  
}
require(gridExtra)

grid.arrange(NetMHCpan_BA_Rank,NetMHCpan_BA_log,NetMHCpan_BA_BA,
             NetMHCpan_EL_Rank,NetMHCpan_EL_log, MHCFlurry_Rank,
             MHCFlurry_BA, MHCFlurry_Proc, MHCFlurry_ProS,
             ncol=5,nrow=2)



#############Figure S1B##########################
#################################################

######Compared performance vs IEDB##############
#MHC II, using all viruses instead of only SARS because of lack of MHC-II SARS data

Alleles_to_keep=c("HLA-DRB1*01:01", "HLA-DRB1*03:01", "HLA-DRB1*04:01", "HLA-DRB1*07:01", "HLA-DRB1*11:01", "HLA-DRB1*13:01", "HLA-DRB1*15:01",
                  "HLA-DQA1*01:02/DQB1*06:02", "HLA-DQA1*05:01/DQB1*02:01", "HLA-DQA1*02:01/DQB1*02:02", "HLA-DQA1*05:05/DQB1*03:01", 
                  "HLA-DQA1*01:01/DQB1*05:01", "HLA-DQA1*03:01/DQB1*03:02", "HLA-DQA1*03:03/DQB1*03:01", "HLA-DQA1*01:03/DQB1*06:03")

iedb = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/mhc_all_virus_combined.csv")
iedb = iedb[,c(1:3)]

iedb=iedb[which(iedb$allele %in% Alleles_to_keep),]

#Filter by IEDB peptides with <500nM BA
iedb_low = iedb[which(iedb$affinity < cutpoint),]
iedb_low = iedb_low[order(iedb_low$peptide),]

#Read in class II predictions of these IEDB reads
Netii_3.2=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_mhc_all_viruses_all_lengths.txt")
Netii_4.0=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan4.0_standardized_mhc_all_viruses_all_lengths.txt")

#Keep only peptide/allele combinations present in IEDB set (all alleles were run for all peptides)
iedb_unique = paste(iedb_low$peptide, iedb_low$allele, sep = "_")
Netii_unique = paste(Netii_3.2$Peptide, Netii_3.2$Haplotype, sep = "_")

Netii_3.2 = Netii_3.2[which(Netii_unique %in% iedb_unique),]
Netii_4.0 = Netii_4.0[which(Netii_unique %in% iedb_unique),]

#Create set of all peptide/allele combinations (same in NetMHCIIpan 3.2 vs 4.0)
all = paste(Netii_4.0$Peptide, Netii_4.0$Haplotype, sep = "_")

#Add features to IEDB matrix
iedb_low = iedb_low[which(iedb_unique %in% all),]
iedb_low$NetMHCIIpan3.2_Rank=NA
iedb_low$NetMHCIIpan3.2_log=NA
iedb_low$NetMHCIIpan3.2_BA=NA

iedb_low$NetMHCIIpan4.0_BA_Rank=NA
iedb_low$NetMHCIIpan4.0_BA_Score=NA
iedb_low$NetMHCIIpan4.0_BA=NA
iedb_low$NetMHCIIpan4.0_EL_Rank=NA
iedb_low$NetMHCIIpan4.0_EL_Score=NA

for(n in 1:nrow(iedb_low)){
  print(n)
  
  iedb_low$NetMHCIIpan3.2_Rank[n] = Netii_3.2$Rank[which((Netii_3.2$Peptide == iedb_low$peptide[n])&(Netii_3.2$Haplotype == iedb_low$allele[n]))]
  iedb_low$NetMHCIIpan3.2_log[n] = Netii_3.2$`1-log50k`[which((Netii_3.2$Peptide == iedb_low$peptide[n])&(Netii_3.2$Haplotype == iedb_low$allele[n]))]
  iedb_low$NetMHCIIpan3.2_BA[n] = Netii_3.2$Binding_affinity[which((Netii_3.2$Peptide == iedb_low$peptide[n])&(Netii_3.2$Haplotype == iedb_low$allele[n]))]
  
  iedb_low$NetMHCIIpan4.0_BA_Rank[n] = Netii_4.0$BA_Rank[which((Netii_4.0$Peptide == iedb_low$peptide[n])&(Netii_4.0$Haplotype == iedb_low$allele[n]))]
  iedb_low$NetMHCIIpan4.0_BA_Score[n] = Netii_4.0$BA_Score[which((Netii_4.0$Peptide == iedb_low$peptide[n])&(Netii_4.0$Haplotype == iedb_low$allele[n]))]
  iedb_low$NetMHCIIpan4.0_BA[n] = Netii_4.0$Binding_affinity[which((Netii_4.0$Peptide == iedb_low$peptide[n])&(Netii_4.0$Haplotype == iedb_low$allele[n]))]
  iedb_low$NetMHCIIpan4.0_EL_Rank[n] = Netii_4.0$EL_Rank[which((Netii_4.0$Peptide == iedb_low$peptide[n])&(Netii_4.0$Haplotype == iedb_low$allele[n]))]
  iedb_low$NetMHCIIpan4.0_EL_Score[n] = Netii_4.0$EL_Score[which((Netii_4.0$Peptide == iedb_low$peptide[n])&(Netii_4.0$Haplotype == iedb_low$allele[n]))]
  
}

iedb_melt = melt(iedb_low, id.vars = c(1:3), measure.vars = c(4:11))
colnames(iedb_melt)[4:5] = c("Program", "Value")
iedb_melt = iedb_melt[complete.cases(iedb_melt),]

#iedb_melt = iedb_melt[which(iedb_melt$Program %in% c('NetMHCpan_BA_BA',  'NetMHCpan_BA_Rank', "NetMHCpan_EL_Rank", 'MHCFlurry_BA', 'MHCFlurry_Rank', "MHCFlurry_Proc", "MHCFlurry_ProS"))]
#iedb_melt = iedb_melt[which(iedb_melt$Program %in% c('NetMHCpan_BA_BA', 'MHCFlurry_BA'))]
fwrite(iedb_melt, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/IEDB_all_virus_HLAII_filtered.txt") 

#Generating individual plots showing correlation between prediction feature and IEDB nM, with Spearman correlation shown
for(z in unique(iedb_melt$Program)){
  print(z)
  if(z %in% unique(iedb_melt$Program)[1:3]){
    col = viridis(3)[1]
  }else if(z %in% unique(iedb_melt$Program)[4:8]){
    col = viridis(3)[2]
  }
  
  sub_melt = iedb_melt[which(iedb_melt$Program ==z),]
  test = cor.test(sub_melt$affinity, sub_melt$Value, method = "spearman") #Will throw "cannot compute exact p-value" warning, this is expected. Will still give good estimate.
  
  p= ggplot(data = sub_melt)+
    geom_point(aes(x = affinity, y= Value), color = col, alpha=.6, show.legend = F)+
    #scale_color_viridis_d(name = "Software")+
    #theme(text=element_text(face="bold",size=20,colour="black")) +
    scale_y_log10()+#breaks = c(0.001, 0.01, .1, 1, 10, 100, 1000, 10000))+
    scale_x_log10(breaks = c(.5,5,50,500))+
    ggnewscale::new_scale("color")+
    labs(x="", y="", title = z) +
    annotate("text", label = paste0("Rho: ",round(test$estimate, digits = 3),"\n","p<0.001"),
             x = 0.5*(min(sub_melt$affinity) + max(sub_melt$affinity)), y = min(sub_melt$Value), vjust = 0, hjust=1)
  
  
  assign(paste0(z),p + geom_smooth(data = ggplot_build(p)$data[[1]], 
                                   mapping = aes(x=10^x, y= 10^y), color = "black", method = "lm",se=FALSE, show.legend = F) )
  
}

require(gridExtra)

grid.arrange(NetMHCIIpan3.2_Rank,NetMHCIIpan3.2_log,NetMHCIIpan3.2_BA,
             NetMHCIIpan4.0_BA_Rank,NetMHCIIpan4.0_BA_Score, NetMHCIIpan4.0_BA,
             NetMHCIIpan4.0_EL_Rank, NetMHCIIpan4.0_EL_Score,
             ncol=4, nrow=2)


#########Figure 2A##################################
####################################################

#Reading in melted matrix from above
iedb_melt = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/IEDB_SARS2_HLAI_filtered_competitive_radioactivity.txt") 

###MHCI
#iedb_BA = iedb_melt[which(iedb_melt$Program %in% c('NetMHCpan_BA_BA'))]

#Reading in only the NetMHCpan -BA binding affinity data -- best correlation

Net_BA=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_all_BA.txt")

iedb = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/mhc_sars2_purified_competitive_radioactivity.csv")
iedb$affinity = as.numeric(iedb$affinity)

iedb_unique = paste(iedb$peptide, iedb$allele, sep = "_")
Net_BA_unique = paste(Net_BA$Peptide, Net_BA$Haplotype, sep = "_")

#Matching predicted BA with IEDB matrix
iedb$NetMHCpan_BA_BA = NA
for(z in 1:length(iedb_unique)){
  if(iedb_unique[z] %in% Net_BA_unique){
    iedb$NetMHCpan_BA_BA[z] = Net_BA$Binding_affinity[which(Net_BA_unique == iedb_unique[z])]
  }
}

iedb_comp = iedb[complete.cases(iedb),]
iedb_comp = iedb_comp[which(iedb_comp$allele %in% Net_BA$Haplotype),]

#Order alleles as factor
all_haps = c("HLA-A*01:01", "HLA-A*02:01", "HLA-A*03:01", "HLA-A*11:01", "HLA-A*24:02", "HLA-B*07:02", "HLA-B*08:01", 
             "HLA-B*35:01", "HLA-B*44:02", "HLA-B*44:03", "HLA-C*03:04", "HLA-C*04:01", "HLA-C*05:01", "HLA-C*06:02", "HLA-C*07:01", "HLA-C*07:02")
all_haps = all_haps[order(all_haps)]

iedb_comp$allele = factor(iedb_comp$allele, levels = all_haps, ordered = T)
Net_BA$Haplotype = factor(Net_BA$Haplotype, levels = all_haps, ordered = T)


#####Figure S1C
spec_cutpoint = .9  #90% specificity cutoff
cutpoint = 500 #500nM IEDB cutoff

iedb_comp$Hits = ifelse(iedb_comp$affinity < cutpoint,0,1) #Binding vs non-binding based on 500nM cutpoint
iedb_comp = iedb[complete.cases(iedb),]

#Calculating FDR, sens/spec for BA as a predictor of >/< 500nM, across various BA levels
cutoff = seq(0,cutpoint,.1) 

FDR = c()
Specificity=c()
Sensitivity = c()

for(n in cutoff){
  print(n)
  FDR = c(FDR, length(which(iedb_comp$affinity[which(iedb_comp$NetMHCpan_BA_BA<=n)] > cutpoint))/length(which(iedb_comp$NetMHCpan_BA_BA<=n)))
  Sensitivity = c(Sensitivity,  length(which(iedb_comp$affinity[which(iedb_comp$NetMHCpan_BA_BA<=n)] < cutpoint))/length(which(iedb_comp$affinity<cutpoint)) )
  Specificity = c(Specificity,  length(which(iedb_comp$affinity[which(iedb_comp$NetMHCpan_BA_BA>n)] > cutpoint))/length(which(iedb_comp$affinity>cutpoint)) )
  
}

plot_cutoff = data.table(cutoff, FDR, Sensitivity, Specificity)
HLAI_nM_cut = plot_cutoff$cutoff[max(which(plot_cutoff$Specificity > spec_cutpoint))]


ggplot(data=plot_cutoff, aes(x=cutoff, y=Specificity))+
  geom_point()+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x="Binding affinity (nM)", y= "Specificity", title = "NetMHCpan v4.0 vs IEDB performance\nSARS peptides, IEDB measured affinity < 500nM")+
  geom_hline(yintercept = c(spec_cutpoint))+
  geom_vline(xintercept = HLAI_nM_cut)+
  annotate("text", label = spec_cutpoint, x = 0, y = spec_cutpoint, size = 10, hjust=0, vjust = 0)+
  annotate("text", label = HLAI_nM_cut, x = HLAI_nM_cut, y = min(plot_cutoff$Specificity), size = 10, hjust=0, vjust = 0, angle = 90)



##########################Plotting 2A#############

#Allows to plot axes in scientific log10 format
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

#Define HLA-I allele colors based on factor order, HLA-A (blue), HLA-B (red), HLA-C (green)
HLAI_col= c(brewer.pal(7, 'Blues')[3:7],
            brewer.pal(7, 'Reds')[3:7],
            brewer.pal(8, 'Greens')[3:8])

#Scatterplot of IEDB training data vs netMHCpan BA 
p=ggplot(data=iedb_comp, aes(x=NetMHCpan_BA_BA, y=affinity, color = allele))+
  geom_point()+
  scale_color_manual(values = HLAI_col,
                     name="Allele", drop = F, guide = guide_legend(override.aes = list(color = "white")))+
  scale_y_log10(breaks = c(0.001, 0.01, .1, 1, 10, 100, 1000, 10000,100000), label=scientific_10)+
  scale_x_continuous(trans= 'log10', breaks = c(1,10,100,1000,10000,100000,1000000), label=scientific_10, limits= c(1,100000))+
  geom_vline(xintercept = HLAI_nM_cut)+
  geom_hline(yintercept = cutpoint)+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x="NetMHCpan BA (nM)", y="IEDB affinity (nM)")+
  annotate("text",x=5, y = .3, label = paste0("n=", length(which((iedb_comp$NetMHCpan_BA_BA<HLAI_nM_cut)&(iedb_comp$affinity<cutpoint)))), color="red", size=10, hjust=0)+
  annotate("text",x=5, y = 500000, label = paste0("n=", length(which((iedb_comp$NetMHCpan_BA_BA<HLAI_nM_cut)&(iedb_comp$affinity>cutpoint)))), color="red", size=10,hjust=0)+
  annotate("text",x=5000, y = .3, label = paste0("n=", length(which((iedb_comp$NetMHCpan_BA_BA>HLAI_nM_cut)&(iedb_comp$affinity<cutpoint)))), color="red", size=10,hjust=0)+
  annotate("text",x=5000, y = 500000, label = paste0("n=", length(which((iedb_comp$NetMHCpan_BA_BA>HLAI_nM_cut)&(iedb_comp$affinity>cutpoint)))), color="red", size=10,hjust=0)+
  
  annotate("text",x=2.5, y = cutpoint, label = paste0(cutpoint,"nM"), color="black", size=5,hjust=0, vjust=0)+
  annotate("text",x=HLAI_nM_cut, y = .1, label = paste0(HLAI_nM_cut,"nM"), color="black", size=5, angle=90, hjust=0, vjust=0)+
  theme(
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")
  ) 

#Histogram of all SARS-CoV-2 predicted epitopes, with netMHCpan BA cutpoint defined above
hist_top <- ggplot(data=Net_BA)+
  geom_histogram(aes(x= Binding_affinity, fill = Haplotype), binwidth = .01)+
  scale_fill_manual(values = HLAI_col,
                    name="Allele", drop = F)+
  scale_x_continuous(trans= 'log10', breaks = c(1,10,100,1000,10000,100000,1000000), limits= c(1,100000))+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  annotate("text",x=5, y = 300000, label = paste0("n=", length(which(Net_BA$Binding_affinity<HLAI_nM_cut))), color="red", size=10,hjust=0)+
  annotate("text",x=5000, y = 300000, label = paste0("n=", length(which(Net_BA$Binding_affinity>HLAI_nM_cut))), color="red", size=10,hjust=0)+
  geom_vline(xintercept = HLAI_nM_cut)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(y="Count")+
  scale_y_continuous(trans = sqrt_trans(), label=scientific_10)

grid.arrange(hist_top, p, nrow=2, ncol=1)

Net_BA_filt = Net_BA[which(Net_BA$Binding_affinity < HLAI_nM_cut),]

fwrite(Net_BA,"/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_BA.txt")
fwrite(Net_BA_filt, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_filtered_BA.txt")



##############Figure 2B#########################
################################################

#Reading in IEDB viral data and corresponding netMHCIIpan 3.2 data, keeping in 12-20mers with alleles used in this study
Netii = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_mhc_all_viruses_all_lengths.txt")

iedb_v = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/mhc_all_virus_combined.csv")
iedb_v = iedb_v[which(iedb_v$allele%in% Alleles_to_keep),]
iedb_v = iedb_v[which(nchar(iedb_v$peptide) %in% c(12:20)),]

iedb_v$affinity = as.numeric(iedb_v$affinity)

iedb_unique = paste(iedb_v$peptide, iedb_v$allele, sep = "_")
Netii_unique = paste(Netii$Peptide, Netii$Haplotype, sep = "_")

iedb_v$NetMHCIIpan_BA = NA
for(z in 1:length(iedb_unique)){
  if(iedb_unique[z] %in% Netii_unique){
    iedb_v$NetMHCIIpan_BA[z] = Netii$Binding_affinity[which(Netii_unique == iedb_unique[z])]
  }
}
iedb_comp2 = iedb_v[complete.cases(iedb_v),]

#Define factor levels for alleles
all_haps2= c("HLA-DRB1*01:01", "HLA-DRB1*03:01", "HLA-DRB1*04:01", "HLA-DRB1*07:01", "HLA-DRB1*11:01",
             "HLA-DRB1*13:01", "HLA-DRB1*15:01", "HLA-DQA1*01:02/DQB1*06:02", "HLA-DQA1*05:01/DQB1*02:01", "HLA-DQA1*02:01/DQB1*02:02",
             "HLA-DQA1*05:05/DQB1*03:01", "HLA-DQA1*01:01/DQB1*05:01", "HLA-DQA1*03:01/DQB1*03:02", "HLA-DQA1*03:03/DQB1*03:01", "HLA-DQA1*01:03/DQB1*06:03")
all_haps2 = all_haps2[order(all_haps2)]

iedb_comp2$allele = factor(iedb_comp2$allele, levels = all_haps2, ordered = T)
Netii$Haplotype = factor(Netii$Haplotype, levels = all_haps2, ordered = T)


###Figure S1D

#Determine netMHCIIpan BA cutpoint to reach 90% specificty within IEDB viral data, predicting for >/< 500nM IEDB BA
spec_cutpoint = .9
cutpoint = 500
iedb_comp2$Hits = ifelse(iedb_comp2$affinity < cutpoint,0,1)

cutoff = seq(0,cutpoint,1)

FDR = c()
Specificity=c()
Sensitivity = c()

for(n in cutoff){
  print(n)
  FDR = c(FDR, length(which(iedb_comp2$affinity[which(iedb_comp2$NetMHCIIpan_BA<=n)] > cutpoint))/length(which(iedb_comp2$NetMHCIIpan_BA<=n)))
  Sensitivity = c(Sensitivity,  length(which(iedb_comp2$affinity[which(iedb_comp2$NetMHCIIpan_BA<=n)] < cutpoint))/length(which(iedb_comp2$affinity<cutpoint)) )
  Specificity = c(Specificity,  length(which(iedb_comp2$affinity[which(iedb_comp2$NetMHCIIpan_BA>n)] > cutpoint))/length(which(iedb_comp2$affinity>cutpoint)) )
}


plot_cutoff = data.table(cutoff, FDR, Sensitivity, Specificity)
HLAII_nM_cut = plot_cutoff$cutoff[max(which(plot_cutoff$Specificity > spec_cutpoint))]

ggplot(data=plot_cutoff, aes(x=cutoff, y=Specificity))+
  geom_point()+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x="Binding affinity (nM)", y= "Specificity", title = "NetMHCIIpan v3.2 vs IEDB performance\nViral peptides, IEDB measured affinity < 500nM")+
  geom_hline(yintercept = c(spec_cutpoint))+
  geom_vline(xintercept = HLAII_nM_cut)+
  annotate("text", label = spec_cutpoint, x = 0, y = spec_cutpoint, size = 10, hjust=0, vjust = 0)+
  annotate("text", label = HLAII_nM_cut, x = HLAII_nM_cut, y = min(plot_cutoff$Specificity), size = 10, hjust=0, vjust = 0, angle = 90)


#######Plotting figure 2B#####################

#Read in SARS netMHCIIpan 3.2 data
Netii_SARS=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_all_BA.txt")
Netii_SARS$Haplotype = factor(Netii_SARS$Haplotype, levels = all_haps2, ordered = T)

#Set HLA-II allele colors corresponding to factor levels -- DQ (orange), DR (purple)
HLAII_col = c(brewer.pal(9, 'Oranges')[2:9],
              brewer.pal(9, 'Purples')[3:9])

#Axes plot in scientific log 10 format
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

#Scatterplot showing IEDB viral read nM and netMHCIIpan 3.2 BA
p2=ggplot(data=iedb_comp2, aes(x=NetMHCIIpan_BA, y=affinity, color = allele))+
  geom_point()+
  scale_color_manual(values = HLAII_col,
                     name="Allele", drop = F, guide = guide_legend(override.aes = list(color = "white")))+
  scale_y_log10(breaks = c(0.001, 0.01, .1, 1, 10, 100, 1000, 10000,100000), label=scientific_10)+
  scale_x_log10(breaks = c(1,10,100,1000,10000,100000,1000000), label=scientific_10, limits=c(5,100000))+
  geom_vline(xintercept = HLAII_nM_cut)+
  geom_hline(yintercept = cutpoint)+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x="NetMHCIIpan BA (nM)", y="IEDB affinity (nM)")+
  annotate("text",x=5, y = .3, label = paste0("n=", length(which((iedb_comp2$NetMHCIIpan_BA<HLAII_nM_cut)&(iedb_comp2$affinity<cutpoint)))), color="red", size=10, hjust=0)+
  annotate("text",x=5, y = 100000, label = paste0("n=", length(which((iedb_comp2$NetMHCIIpan_BA<HLAII_nM_cut)&(iedb_comp2$affinity>cutpoint)))), color="red", size=10,hjust=0)+
  annotate("text",x=5000, y = .3, label = paste0("n=", length(which((iedb_comp2$NetMHCIIpan_BA>HLAII_nM_cut)&(iedb_comp2$affinity<cutpoint)))), color="red", size=10,hjust=0)+
  annotate("text",x=5000, y = 100000, label = paste0("n=", length(which((iedb_comp2$NetMHCIIpan_BA>HLAII_nM_cut)&(iedb_comp2$affinity>cutpoint)))), color="red", size=10,hjust=0)+
  
  annotate("text",x=5, y = 350, label = paste0(cutpoint,"nM"), color="black", size=5,hjust=0)+
  annotate("text",x=315, y = .1, label = paste0(HLAII_nM_cut,"nM"), color="black", size=5, angle=90, hjust=0, vjust=0)+
  theme(
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")
  ) 

#SARS-CoV-2 predicted epitopes using netMHCIIpan BA data and cutpoint determined above
hist_top2 <- ggplot(data=Netii_SARS)+
  geom_histogram(aes(x= Binding_affinity, fill = Haplotype), binwidth = .01)+
  scale_fill_manual(values= HLAII_col,
                    name="Allele", drop = F)+
  scale_x_log10(breaks = c(1,10,100,1000,10000,100000,1000000), limits=c(5,100000))+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_vline(xintercept = HLAII_nM_cut)+
  annotate("text",x=5, y = 800, label = paste0("n=", length(which(Netii_SARS$Binding_affinity<HLAII_nM_cut))), color="red", size=10,hjust=0)+
  annotate("text",x=5000, y = 800, label = paste0("n=", length(which(Netii_SARS$Binding_affinity>HLAII_nM_cut))), color="red", size=10,hjust=0)+
  labs(y="Count")+
  scale_y_continuous(labels = scientific_10, breaks = c(0,200,400,600,800,10000))


grid.arrange(hist_top2, p2, nrow=2,  ncol=1)

Netii_SARS_filt = Netii_SARS[which(Netii_SARS$Binding_affinity < HLAII_nM_cut),]

fwrite(Netii_SARS, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_BA.txt")
fwrite(Netii_SARS_filt, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_filtered_BA.txt")



####################Figure 2C -- Lollipop plot ##################
################################################################

#Read in HLA-I and II SARS-CoV-2 predicted epitopes
Net_BA=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_BA.txt")
Netii_SARS=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_BA.txt")

Netii_SARS_filt=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_filtered_BA.txt")
Net_BA_filt=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_filtered_BA.txt")

#Read in population frequencies of each allele, format allele names to be consistent
Freq=fread(paste0(WORKING_DIR, "HLA_freq.txt"))
Freq$Haplotype = str_replace_all(Freq$Haplotype, "DRB1_", "HLA-DRB1*")
Freq$Haplotype = str_replace_all(Freq$Haplotype, "DQA1", "DQA1*")
Freq$Haplotype = str_replace_all(Freq$Haplotype, "-DQB1", "/DQB1*")
Freq$Haplotype = str_replace_all(Freq$Haplotype, "HLA-A", "HLA-A*")
Freq$Haplotype = str_replace_all(Freq$Haplotype, "HLA-B", "HLA-B*")
Freq$Haplotype = str_replace_all(Freq$Haplotype, "HLA-C", "HLA-C*")
Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DR")] = paste0(substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DR")],1,11), ":",
                                                                      substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DR")],12,13))
Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DQ")] = paste0(substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DQ")],1,11), ":",
                                                                      substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DQ")],12,21),":",
                                                                      substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DQ")],22,23))

#For DQ reads, values were derived from a study which was biased for caucasian americans.  
#Unfortunately, paired DQ frequency data was not available for a balanced US population.
Freq$US[which(is.na(Freq$US))] = Freq$Cauc_Am[which(is.na(Freq$US))]

#Cast tables, collapsing by unique peptides
Net_BA_cast = dcast.data.table(Net_BA, Peptide + Protein + Start ~ Haplotype, value.var = "Binding_affinity")
Net_BA_cast = Net_BA_cast[which(Net_BA_cast$Peptide %in% Net_BA_filt$Peptide),]

Netii_BA_cast = dcast.data.table(Netii_SARS, Peptide + Protein + Start ~ Haplotype, value.var = "Binding_affinity")
Netii_BA_cast = Netii_BA_cast[which(Netii_BA_cast$Peptide %in% Netii_SARS_filt$Peptide),]


#Calculating population frequencies

registerDoMC(48)

Pop_freq_I = foreach(n=1:nrow(Net_BA_cast), .combine=rbind) %dopar% {
  Binding_haplotype=colnames(Net_BA_cast)[4:19][which(Net_BA_cast[n,4:19] < HLAI_nM_cut)]
  Frequency=1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% Binding_haplotype )]/100,2)))
  data.frame(paste0(Binding_haplotype, collapse = ","),Frequency)
}

Pop_freq_II = foreach(n=1:nrow(Netii_BA_cast), .combine=rbind) %dopar% {
  Binding_haplotype=colnames(Netii_BA_cast)[4:18][which(Netii_BA_cast[n,4:18] < HLAII_nM_cut)]
  Frequency=1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% Binding_haplotype)]/100,2)))
  data.frame(paste0(Binding_haplotype, collapse = ','),Frequency)
}

Net_BA_cast = cbind(Net_BA_cast, Pop_freq_I)
colnames(Net_BA_cast)[20] = "Binding_haplotype" 
Netii_BA_cast = cbind(Netii_BA_cast, Pop_freq_II)
colnames(Netii_BA_cast)[19] = "Binding_haplotype"

#Calculating co-epitopes

Best_II_hit = data.table()

for(n in 1:nrow(Net_BA_cast)){
  print(n)
  pep1 = Net_BA_cast$Peptide[n]
  pep1_start = Net_BA_cast$Start[n]
  hits = which(str_detect(string = Netii_BA_cast$Peptide, pattern = pep1))
  if(length(hits) > 0){
    max_freq=which(Netii_BA_cast$Frequency[which(str_detect(string = Netii_BA_cast$Peptide, pattern = pep1))] == 
                     max(Netii_BA_cast$Frequency[which(str_detect(string = Netii_BA_cast$Peptide, pattern = pep1))]))
    hits = hits[max_freq]
    pep2_freq = Netii_BA_cast$Frequency[hits]
    pep2_haps = Netii_BA_cast$Binding_haplotype[hits]
    pep2 = Netii_BA_cast$Peptide[hits]
    pep2_start = Netii_BA_cast$Start[hits]
    #print(paste0(n, ": ", nrow(data.table(pep1, pep2, pep2_haps, pep2_freq))))
    Best_II_hit = rbind(Best_II_hit, data.table(Net_BA_cast$Protein[n],pep1_start, pep1, Net_BA_cast$Binding_haplotype[n], Net_BA_cast$Frequency[n], pep2_start, pep2, pep2_haps, pep2_freq))
  }else{
    Best_II_hit = rbind(Best_II_hit, data.table(Net_BA_cast$Protein[n],pep1_start, pep1, Net_BA_cast$Binding_haplotype[n], Net_BA_cast$Frequency[n], NA, NA, NA, NA), use.names=F)
  }
}

Net_coepitopes = cbind(Best_II_hit, (Best_II_hit$V5*Best_II_hit$pep2_freq))
colnames(Net_coepitopes) = c("Protein", "c1_start", "c1_peptide", "c1_haplotypes", "c1_freq", "c2_start", "c2_peptide", "c2_haplotypes", "c2_freq", "Total_freq")
Net_coepitopes = Net_coepitopes[which(complete.cases(Net_coepitopes)),]

fwrite(Net_coepitopes, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Coepitopes_human.txt")
fwrite(Net_BA_cast,  "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_cast_BA.txt")
fwrite(Netii_BA_cast,  "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_cast_BA.txt")

#Reading in murine data
mouse1 = fread(paste0(WORKING_DIR, "COVID_murine_NetMHCpan_unfilt.xls"))
mouse1 = mouse1[which(mouse1$NB>0),]
mouse2 = fread(paste0(WORKING_DIR, "COVID_murine_NetMHCIIpan_unfilt.xls"))
mouse2 = mouse2[which(mouse2$NB>0),]

Net_coepitopes=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Coepitopes_human.txt")
Net_coepitopes$c1_freq = Net_coepitopes$c1_freq*100
Net_coepitopes$c2_freq = Net_coepitopes$c2_freq*100

Net_BA_cast = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_cast_BA.txt")
Netii_BA_cast = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_cast_BA.txt")

#Adding murine coverage for each human epitope
Net_coepitopes$murine = "None"
Net_coepitopes$murine[which(Net_coepitopes$c1_peptide %in% mouse1$Peptide)] = "MHC-I" 
Net_coepitopes$murine[which(Net_coepitopes$c2_peptide %in% mouse2$Peptide)] = "MHC-II" 
Net_coepitopes$murine[which((Net_coepitopes$c2_peptide %in% mouse2$Peptide) & (Net_coepitopes$c1_peptide %in% mouse1$Peptide))] = "Both" 
Net_coepitopes$murine = factor(Net_coepitopes$murine, levels = c("MHC-I", "MHC-II", "Both", "None"))

Net_coepitopes$Total_freq=Net_coepitopes$Total_freq*100

#Add in Tc literature data

paper_tc = fread(paste0(WORKING_DIR, "Supplemental/Standardized T cell epitopes.txt"))
paper_tc = paper_tc[!duplicated(paper_tc),]
seq =  fread(paste0(WORKING_DIR, "AA_sequence_combined.txt"), header=F)

paper_tc$Start = NA
paper_tc$End = NA

#Checking if each T cell epitope from literature is present in SARS-CoV-2 reference.  If so, add coordinates
for (y in 1:nrow(paper_tc)){
  if(str_detect(string = seq$V1 ,paper_tc$Peptide[y])){
    print(y)
    paper_tc$Start[y] = gregexpr(text = seq$V1, pattern = paper_tc$Peptide[y])[[1]][1]
    paper_tc$End[y] = (paper_tc$Start[y] + nchar(paper_tc$Peptide[y]) - 1)
  }
  
}
paper_tc = paper_tc[which(!is.na(paper_tc$Start)),]
paper_tc$Type = "T cell"

#See if Tc epitope from literature is present in our predicted co-epitope set
Net_coepitopes$Lit_overlap = 0
Net_coepitopes$Lit_overlap[which((Net_coepitopes$c1_peptide %in% paper_tc$Peptide)  | (Net_coepitopes$c2_peptide %in% paper_tc$Peptide))]=1

## Transforming to proteome space
###Reference locations
#orf1ab: 1 - 7096                                                                    
#S: 7097 - 8369                                                                       
#ORF3a: 8370 - 8644                                                                   
#E: 8645 - 8719                                                                       
#M: 8720 - 8941                                                                       
#ORF6: 8942 - 9002                                                                    
#ORF7a: 9003 - 9123                                                                   
#ORF8: 9124 - 9244                                                                    
#N: 9245 - 9663                                                                       
#ORF10: 9664 - 9701

Net_coepitopes$c1_start[which(Net_coepitopes$Protein == "S")] = Net_coepitopes$c1_start[which(Net_coepitopes$Protein == "S")]+7096
Net_coepitopes$c1_start[which(Net_coepitopes$Protein == "ORF3a")] = Net_coepitopes$c1_start[which(Net_coepitopes$Protein == "ORF3a")]+8369
Net_coepitopes$c1_start[which(Net_coepitopes$Protein == "E")] = Net_coepitopes$c1_start[which(Net_coepitopes$Protein == "E")]+8644
Net_coepitopes$c1_start[which(Net_coepitopes$Protein == "M")] = Net_coepitopes$c1_start[which(Net_coepitopes$Protein == "M")]+8719
Net_coepitopes$c1_start[which(Net_coepitopes$Protein == "ORF6")] = Net_coepitopes$c1_start[which(Net_coepitopes$Protein == "ORF6")]+8941
Net_coepitopes$c1_start[which(Net_coepitopes$Protein == "ORF7a")] = Net_coepitopes$c1_start[which(Net_coepitopes$Protein == "ORF7a")]+9002
Net_coepitopes$c1_start[which(Net_coepitopes$Protein == "ORF8")] = Net_coepitopes$c1_start[which(Net_coepitopes$Protein == "ORF8")]+9123
Net_coepitopes$c1_start[which(Net_coepitopes$Protein == "N")] = Net_coepitopes$c1_start[which(Net_coepitopes$Protein == "N")]+9244
Net_coepitopes$c1_start[which(Net_coepitopes$Protein == "ORF10")] = Net_coepitopes$c1_start[which(Net_coepitopes$Protein == "ORF10")]+9663

Net_coepitopes$c2_start[which(Net_coepitopes$Protein == "S")] = Net_coepitopes$c2_start[which(Net_coepitopes$Protein == "S")]+7096
Net_coepitopes$c2_start[which(Net_coepitopes$Protein == "ORF3a")] = Net_coepitopes$c2_start[which(Net_coepitopes$Protein == "ORF3a")]+8369
Net_coepitopes$c2_start[which(Net_coepitopes$Protein == "E")] = Net_coepitopes$c2_start[which(Net_coepitopes$Protein == "E")]+8644
Net_coepitopes$c2_start[which(Net_coepitopes$Protein == "M")] = Net_coepitopes$c2_start[which(Net_coepitopes$Protein == "M")]+8719
Net_coepitopes$c2_start[which(Net_coepitopes$Protein == "ORF6")] = Net_coepitopes$c2_start[which(Net_coepitopes$Protein == "ORF6")]+8941
Net_coepitopes$c2_start[which(Net_coepitopes$Protein == "ORF7a")] = Net_coepitopes$c2_start[which(Net_coepitopes$Protein == "ORF7a")]+9002
Net_coepitopes$c2_start[which(Net_coepitopes$Protein == "ORF8")] = Net_coepitopes$c2_start[which(Net_coepitopes$Protein == "ORF8")]+9123
Net_coepitopes$c2_start[which(Net_coepitopes$Protein == "N")] = Net_coepitopes$c2_start[which(Net_coepitopes$Protein == "N")]+9244
Net_coepitopes$c2_start[which(Net_coepitopes$Protein == "ORF10")] = Net_coepitopes$c2_start[which(Net_coepitopes$Protein == "ORF10")]+9663

Net_BA_cast$Start[which(Net_BA_cast$Protein == "S")] = Net_BA_cast$Start[which(Net_BA_cast$Protein == "S")]+7096
Net_BA_cast$Start[which(Net_BA_cast$Protein == "ORF3a")] = Net_BA_cast$Start[which(Net_BA_cast$Protein == "ORF3a")]+8369
Net_BA_cast$Start[which(Net_BA_cast$Protein == "E")] = Net_BA_cast$Start[which(Net_BA_cast$Protein == "E")]+8644
Net_BA_cast$Start[which(Net_BA_cast$Protein == "M")] = Net_BA_cast$Start[which(Net_BA_cast$Protein == "M")]+8719
Net_BA_cast$Start[which(Net_BA_cast$Protein == "ORF6")] = Net_BA_cast$Start[which(Net_BA_cast$Protein == "ORF6")]+8941
Net_BA_cast$Start[which(Net_BA_cast$Protein == "ORF7a")] = Net_BA_cast$Start[which(Net_BA_cast$Protein == "ORF7a")]+9002
Net_BA_cast$Start[which(Net_BA_cast$Protein == "ORF8")] = Net_BA_cast$Start[which(Net_BA_cast$Protein == "ORF8")]+9123
Net_BA_cast$Start[which(Net_BA_cast$Protein == "N")] = Net_BA_cast$Start[which(Net_BA_cast$Protein == "N")]+9244
Net_BA_cast$Start[which(Net_BA_cast$Protein == "ORF10")] = Net_BA_cast$Start[which(Net_BA_cast$Protein == "ORF10")]+9663

Netii_BA_cast$Start[which(Netii_BA_cast$Protein == "S")] = Netii_BA_cast$Start[which(Netii_BA_cast$Protein == "S")]+7096
Netii_BA_cast$Start[which(Netii_BA_cast$Protein == "ORF3a")] = Netii_BA_cast$Start[which(Netii_BA_cast$Protein == "ORF3a")]+8369
Netii_BA_cast$Start[which(Netii_BA_cast$Protein == "E")] = Netii_BA_cast$Start[which(Netii_BA_cast$Protein == "E")]+8644
Netii_BA_cast$Start[which(Netii_BA_cast$Protein == "M")] = Netii_BA_cast$Start[which(Netii_BA_cast$Protein == "M")]+8719
Netii_BA_cast$Start[which(Netii_BA_cast$Protein == "ORF6")] = Netii_BA_cast$Start[which(Netii_BA_cast$Protein == "ORF6")]+8941
Netii_BA_cast$Start[which(Netii_BA_cast$Protein == "ORF7a")] = Netii_BA_cast$Start[which(Netii_BA_cast$Protein == "ORF7a")]+9002
Netii_BA_cast$Start[which(Netii_BA_cast$Protein == "ORF8")] = Netii_BA_cast$Start[which(Netii_BA_cast$Protein == "ORF8")]+9123
Netii_BA_cast$Start[which(Netii_BA_cast$Protein == "N")] = Netii_BA_cast$Start[which(Netii_BA_cast$Protein == "N")]+9244
Netii_BA_cast$Start[which(Netii_BA_cast$Protein == "ORF10")] = Netii_BA_cast$Start[which(Netii_BA_cast$Protein == "ORF10")]+9663


#Calculating fit for class I and II, showing number of reads per AA location along the proteome
class_I_coverage=c()
class_II_coverage = c()

for (n in 1:9701){
  class_I_coverage = c(class_I_coverage, length(which((Net_BA_cast$Start<=n) & ((Net_BA_cast$Start+nchar(Net_BA_cast$Peptide)-1)>=n)  )))
  class_II_coverage = c(class_II_coverage, length(which((Netii_BA_cast$Start<=n) & ((Netii_BA_cast$Start+nchar(Netii_BA_cast$Peptide)-1)>=n)  )))
}

coverage=rbind(data.table(seq(1,9701,1),"HLA-I",class_I_coverage),
               data.table(seq(1,9701,1),"HLA-II",class_II_coverage), use.names=F)
colnames(coverage) = c("Position", "Class", "Coverage")


#Plotting Figure 2C
Min_frequency = 25 #Minimum population frequency to show on plot
label_frequency = 45 #Minimum population frequency to label on plot

#Y-axis upper and lower limits for plotting
y_ll=20
y_ul=55

p=ggplot(data=Net_coepitopes[which(Net_coepitopes$Total_freq>Min_frequency)]) + 
  geom_segment(data=Net_coepitopes[which(Net_coepitopes$Total_freq>Min_frequency)], aes(x=c2_start, xend=c2_start, y=y_ll, yend=Total_freq),
               color=ifelse(Net_coepitopes[which(Net_coepitopes$Total_freq>Min_frequency)]$Lit_overlap==1,'red','grey'), alpha=0.5, size=1.1) +
  geom_point(aes(x=c2_start, y=Total_freq, color=murine), alpha=.75, size = 5) +
  
  scale_color_manual(values = c("tomato", "cadetblue", "purple", "grey"), name="Murine\noverlap")+
  scale_size(range = c(0,10), breaks = c(20,40,60,80,90))+
  scale_y_continuous(limits = c((y_ll-13),y_ul), breaks = c(20,30,40,50,60))+
  
  geom_label_repel(aes(label=ifelse((Total_freq>label_frequency),paste0(as.character(c2_peptide),'\n',as.character(c1_peptide)),''), x=c2_start, y=Total_freq),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  
  geom_rect(aes(xmin=0, xmax=7095, ymin=(y_ll), ymax=(y_ll-4)), fill=viridis(10)[1], color="black", size=0.1) +
  geom_text(aes(x=7095/2, y=(y_ll-5), label="orf1ab", angle=0), size=5, color = "black") +
  geom_text(aes(x=1000, y=(y_ll), label="-1000", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=2000, y=(y_ll), label="-2000", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=3000, y=(y_ll), label="-3000", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=4000, y=(y_ll), label="-4000", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=5000, y=(y_ll), label="-5000", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=6000, y=(y_ll), label="-6000", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=7000, y=(y_ll), label="-7000", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  
  geom_rect(aes(xmin=7096, xmax=8368, ymin=(y_ll), ymax=(y_ll-4)), fill=viridis(10)[2], color="black", size=0.1) +
  geom_text(aes(x=7096+(8368-7096)/2, y=(y_ll-5), label="S", angle=0), size=5, color = "black") + 
  geom_text(aes(x=7096+200, y=(y_ll), label="-200", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=7096+500, y=(y_ll), label="-500", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=7096+800, y=(y_ll), label="-800", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=7096+1100, y=(y_ll), label="-1100", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  
  geom_rect(aes(xmin=8369, xmax=8644, ymin=(y_ll-0), ymax=(y_ll-4)), fill=viridis(10)[3], color="black", size=0.1) +
  geom_text(aes(x=8368+(8644-8368)/2, y=(y_ll-5), label="ORF3a", angle=0), size=2.5, color = "black") +
  geom_text(aes(x=8369+100, y=(y_ll), label="-100", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=8369+200, y=(y_ll), label="-200", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  
  geom_rect(aes(xmin=8645, xmax=8718, ymin=(y_ll), ymax=(y_ll-4)), fill=viridis(10)[4], color="black", size=0.1) +
  geom_text(aes(x=8645+(8718-8645)/2, y=(y_ll-5)), label="E", angle=0, size=2.5, color = "black") + 
  
  geom_rect(aes(xmin=8719, xmax=8940, ymin=(y_ll-0), ymax=(y_ll-4)), fill=viridis(10)[5], color="black", size=0.1) +
  geom_text(aes(x=8719+(8940-8719)/2, y=(y_ll-5)), label="M", angle=0, size=2.5, color = "black") +
  
  geom_rect(aes(xmin=8941, xmax=9001, ymin=(y_ll), ymax=(y_ll-4)), fill=viridis(10)[6], color="black", size=0.1) +
  geom_text(aes(x=8941+(9001-8941)/2, y=(y_ll-6)), label="ORF6", angle=0, size=2.5, color = "black") +
  
  geom_rect(aes(xmin=9002, xmax=9122, ymin=(y_ll-0), ymax=(y_ll-4)), fill=viridis(10)[7], color="black", size=0.1) +
  geom_text(aes(x=9002+(9122-9002)/2, y=(y_ll-7)), label="ORF7a", angle=0, size=2.5, color = "black") +
  
  geom_rect(aes(xmin=9123, xmax=9243, ymin=(y_ll), ymax=(y_ll-4)), fill=viridis(10)[8], color="black", size=0.1) +
  geom_text(aes(x=9123+(9243-9123)/2, y=(y_ll-5)), label="ORF8", angle=0, size=2.5, color = "black") +
  
  geom_rect(aes(xmin=9244, xmax=9662, ymin=(y_ll-0), ymax=(y_ll-4)), fill=viridis(10)[9], color="black", size=0.1) +
  geom_text(aes(x=9244+(9662-9244)/2, y=(y_ll-5)), label="N", angle=0, size=2.5, color = "black") +
  geom_text(aes(x=9244+100, y=(y_ll), label="-100", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=9244+200, y=(y_ll), label="-200", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=9244+300, y=(y_ll), label="-300", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  
  geom_rect(aes(xmin=9663, xmax=9701, ymin=(y_ll), ymax=(y_ll-4)), fill=viridis(10)[10], color="black", size=0.1) +
  geom_text(aes(x=9663+(9701-9663)/2, y=(y_ll-5)), label="ORF10", angle=0, size=2.5, color = "black") +
  
  #Literature track
  geom_rect(data = paper_tc,aes(ymin=y_ll-8, ymax=y_ll-12, xmin=Start, xmax=End), fill = alpha(c("red"), .5),size=0)+
  geom_text(aes(x = 9670, y = y_ll-12, label = "Literature"), hjust = 0, size = 3.5)+
  labs(y="Co-epitope pop. frequency")+
  theme(text=element_text(face="bold",size=20,colour="black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 

sm= ggplot(data=coverage,aes(x = Position, y=Coverage, color = Class))+geom_smooth( method = 'loess',span=.1)+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  scale_y_continuous(breaks = c(0,4,8,12))+
  labs(x="Position across SARS-CoV-2 proteome", y="Count")

png("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Figure_2C.png", width = 20, height = 8.3, units = "in", res=300)
grid.arrange(p, sm, ncol=1, nrow=2, heights = c(5,2))
dev.off()

#Correlate fits for HLA-I and II
g1=as.data.table(ggplot_build(ggplot()+geom_smooth(data=coverage[which(coverage$Class=="HLA-I")],aes(x=Position, y=Coverage), method = 'loess', span=0.1 ))$data)
g2=as.data.table(ggplot_build(ggplot()+geom_smooth(data=coverage[which(coverage$Class=="HLA-II")],aes(x=Position, y=Coverage), method = 'loess', span=0.1 ))$data)

plot(g1$y, g2$y)
cor.test(g1$y, g2$y)

#####Figure 2D -- Venn diagrams######################
#####################################################

#Class I/II overlap


vd_I <- venneuler(c(HLAI=nrow(Net_BA_cast), I_coep=0,"HLAI&I_coep"=length(unique(Net_coepitopes$c1_peptide))))
vd_I$labels=NA
#vd$colors 

plot(vd_I, col= c(alpha(colour = "orange", .5), alpha(colour = "red",.5)))
text(0.825,.5,paste0("All human\nHLA-I epitopes\n",nrow(Net_BA_cast)))
text(0.425,.5,paste0("Human\nHLA-I co-epitopes\n", length(unique(Net_coepitopes$c1_peptide))))

vd_II <- venneuler(c(HLAII=nrow(Netii_BA_cast), II_coep=0,"HLAII&II_coep"=length(unique(Net_coepitopes$c2_peptide))))
vd_II$labels=NA

plot(vd_II, col= c(alpha(colour = "cyan", .5), alpha(colour = "blue",.5)))
text(0.825,.5,paste0("All human\nHLA-II epitopes\n",nrow(Netii_BA_cast)))
text(0.425,.5,paste0("Human\nHLA-II co-epitopes\n",length(unique(Net_coepitopes$c2_peptide))))

vd_CO <- venneuler(c(Coep=nrow(Net_coepitopes), II_coep=100,"Coep&II_coep"=0))
vd_CO$labels=NA

plot(vd_CO, col= c(alpha(colour = "purple", .5), alpha(colour = "blue",.5)))
text(0.75,.55,paste0("Co-epitopes: ", nrow(Net_coepitopes)))





########Figure 2E - Pie charts#######################
#####################################################

#Defining proteins and their lengths
all_prots=c("orf1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF8","N","ORF10")
prot_len= c(abs(1 - 7096),abs(7097 - 8369),abs(8370 - 8644),abs(8645 - 8719),abs(8720 - 8941),abs(8942 - 9002),abs(9003 - 9123), abs(9124 - 9244),abs(9245 - 9663),abs(9664 - 9701))

Net_coepitopes = Net_coepitopes[order(Net_coepitopes$c2_start),]
prots = factor(c(unique(Net_coepitopes$Protein)), levels = all_prots)

##Raw counts
counts=c()
for(z in prots){
  counts = c(counts, length(which(Net_coepitopes$Protein == z)))
}

prot_dist = data.table(prots, counts)
colnames(prot_dist) = c("Protein", "Count")


prot_dist <- prot_dist %>% 
  arrange(desc(Protein)) %>%
  mutate(Count = Count / sum(prot_dist$Count) *100) %>%
  mutate(ypos = cumsum(Count)- 0.5*Count )

p_unnorm = ggplot(prot_dist, aes(x="", y=Count, fill=Protein)) +
  geom_bar(stat="identity", width=1, show.legend = F) +
  coord_polar("y", start=0)+
  scale_fill_manual(values=c(viridis(10)))+
  geom_label_repel(aes(y = ypos, label = Protein), color = "white", size=8, nudge_x = .65,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'black',
                   show.legend = F)+
  labs(title = "Raw counts")+
  theme(legend.position="none")+
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(face="bold",size=20,colour="black")) 

###Normalized to length
counts=c()
for(z in prots){
  counts = c(counts, length(which(Net_coepitopes$Protein == z))/ prot_len[which(all_prots==z)] )
}

prot_dist = data.table(prots, counts)
colnames(prot_dist) = c("Protein", "Count")


prot_dist <- prot_dist %>% 
  arrange(desc(Protein)) %>%
  mutate(Count = Count / sum(prot_dist$Count) *100) %>%
  mutate(ypos = cumsum(Count)- 0.5*Count )

p_norm = ggplot(prot_dist, aes(x="", y=Count, fill=Protein)) +
  geom_bar(stat="identity", width=1, show.legend = F) +
  coord_polar("y", start=0)+
  scale_fill_manual(values=c(viridis(10)))+
  geom_label_repel(aes(y = ypos, label = Protein), color = "white", size=8, nudge_x = .65,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'black',
                   show.legend = F)+
  labs(title = "Length normalized")+
  theme(legend.position="none")+
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text=element_text(face="bold",size=20,colour="black")) 

grid.arrange(p_unnorm, p_norm, nrow = 1, ncol = 2)


#####Figure 2F - Murine vs human Venn diagrams######################
####################################################################

#Murine overlap


#Defining overlap lengths
murine_1_overlap = length(which(mouse1$Peptide %in% Net_BA_cast$Peptide))
murine_2_overlap = length(which(mouse2$Peptide %in% Netii_BA_cast$Peptide))
murine_1_coep = length(which(mouse1$Peptide %in% Net_coepitopes$c1_peptide))
murine_2_coep = length(which(mouse2$Peptide %in% Net_coepitopes$c2_peptide))
coep_1_overlap = length(which(Net_BA_cast$Peptide %in% Net_coepitopes$c1_peptide))
coep_2_overlap = length(which(Netii_BA_cast$Peptide %in% Net_coepitopes$c2_peptide))

#Class I
vd_1 = venneuler(c(HLAI=(nrow(Net_BA_cast)-murine_1_overlap), MHCI=(nrow(mouse1)-murine_1_overlap),
                   "HLAI&MHCI"=murine_1_overlap))
vd_1$labels<- c(
  paste0("Human\nHLA-I\n",(nrow(Net_BA_cast)-murine_1_overlap)),
  paste("Murine\nMHC-I\n",(nrow(mouse1)-murine_1_overlap))
)

plot(vd_1, col = c(alpha("red",.5), alpha("tomato2",.5)))
text(0.5,.5,paste0(murine_1_overlap))


#Class II
vd_2 = venneuler(c(HLAII=(nrow(Netii_BA_cast)-murine_2_overlap), MHCII=(nrow(mouse2)-murine_2_overlap),
                   "HLAII&MHCII"=murine_2_overlap))
vd_2$labels = NA


plot(vd_2, col=c(alpha("blue",.5), alpha("cadetblue",.5)))
text(0.825,.5,paste0("Human\nHLA-II\n", (nrow(Netii_BA_cast)-murine_2_overlap)))
text(0.15,.5,paste0("Murine\nMHC-II\n", (nrow(mouse2)-murine_2_overlap)))
text(0.425,.5,murine_2_overlap)


#Co-epitopes
vd_co  = venneuler(c(Coepitopes=(nrow(Net_coepitopes)-murine_1_coep-murine_2_coep), 
                     MHCI=(nrow(mouse1)-murine_1_coep),
                     MHCII=(nrow(mouse2)-murine_2_coep),
                     "Coepitopes&MHCI"=murine_1_coep,
                     "Coepitopes&MHCII"=murine_2_coep,
                     "MHCII&MHCI"=0))


vd_co$labels = NA

plot(vd_co, col = c(alpha("purple",.5), alpha("tomato",.5), alpha("cadetblue",.5)))
text(0.475,.75,paste0("Human\nCo-epitopes\n", (nrow(Net_coepitopes)-murine_1_coep-murine_2_coep)))
text(0.675,.4,paste0("Murine\nMHC-I\n", (nrow(mouse1)-murine_1_coep)))
text(0.35,.4,paste0("Murine\nMHC-II\n", (nrow(mouse2)-murine_2_coep)))
text(0.6,.5,murine_1_coep)
text(0.375,.5,murine_2_coep)




###Figure 2G - Histogram of pop freq###################
#######################################################

all_freq = data.table(c(Net_coepitopes$c1_freq, Net_coepitopes$c2_freq), 
                      c(rep("HLA-I", nrow(Net_coepitopes)), rep("HLA-II", nrow(Net_coepitopes))))
colnames(all_freq) = c("Pop. frequency", "HLA")

#Grab most common HLA-I haplotype combos
common_HLAI=round(as.numeric(names(sort(table(Net_coepitopes$c1_freq),decreasing=TRUE)[1:4])), digits = 3)
common_HLAII=round(as.numeric(names(sort(table(Net_coepitopes$c2_freq),decreasing=TRUE)[1:4])), digits=3)

HLAI_lab1=as.character(unique(Net_BA_cast$Binding_haplotype[which(round(Net_BA_cast$Frequency*100, digits = 3) %in% common_HLAI[1])]))
HLAI_lab2=as.character(unique(Net_BA_cast$Binding_haplotype[which(round(Net_BA_cast$Frequency*100, digits = 3) %in% common_HLAI[2])]))
HLAI_lab3=as.character(unique(Net_BA_cast$Binding_haplotype[which(round(Net_BA_cast$Frequency*100, digits = 3) %in% common_HLAI[3])]))
HLAI_lab4=as.character(unique(Net_BA_cast$Binding_haplotype[which(round(Net_BA_cast$Frequency*100, digits = 3) %in% common_HLAI[4])]))


HLAII_lab1=as.character(unique(Netii_BA_cast$Binding_haplotype[which(round(Netii_BA_cast$Frequency*100, digits = 3) %in% common_HLAII[1])]))
HLAII_lab2=as.character(unique(Netii_BA_cast$Binding_haplotype[which(round(Netii_BA_cast$Frequency*100, digits = 3) %in% common_HLAII[2])]))
HLAII_lab3=as.character(unique(Netii_BA_cast$Binding_haplotype[which(round(Netii_BA_cast$Frequency*100, digits = 3) %in% common_HLAII[3])]))
HLAII_lab4=as.character(unique(Netii_BA_cast$Binding_haplotype[which(round(Netii_BA_cast$Frequency*100, digits = 3) %in% common_HLAII[4])]))


ggplot(all_freq, aes(x=`Pop. frequency`, color=HLA, fill=HLA)) +
  geom_histogram(alpha=.3, bins = nrow(all_freq), position="dodge")+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(y="Count") + 
  scale_y_continuous(limits = c(0,2000))+
  scale_x_continuous(breaks = seq(0,100,10))+
  
  # 
  annotate("text", x=common_HLAI[1], y = 1900, label = HLAI_lab1, color = "black", hjust=0, vjust=0, angle = 0)+
  # annotate("text", x=common_HLAI[2], y = 100, label = HLAI_lab2, color = "black", hjust=0, vjust=0, angle = 90)+
  # annotate("text", x=common_HLAI[3], y = 100, label = HLAI_lab3, color = "black", hjust=0, vjust=0, angle = 90)+
  # annotate("text", x=common_HLAI[4], y = 100, label = HLAI_lab4, color = "black", hjust=0, vjust=0, angle = 90)+
  # 
  annotate("text", x=common_HLAII[1], y = 1000, label = HLAII_lab1, color = "black", hjust=0, vjust=0, angle = 0)
# annotate("text", x=common_HLAII[2], y = 100, label = HLAII_lab2, color = "black", hjust=0, vjust=0, angle = 90)+
# annotate("text", x=common_HLAII[3], y = 100, label = paste0(HLAII_lab3[1],"\n",HLAII_lab3[2]), color = "black", hjust=0, vjust=0, angle = 90)+
# annotate("text", x=common_HLAII[4], y = 100, label = HLAII_lab4, color = "black", hjust=0, vjust=0, angle = 90)
# 
fwrite(Net_BA_cast, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Net_BA_cast_fig3.txt")
fwrite(Netii_BA_cast, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Netii_BA_cast_fig3.txt")
fwrite(Net_coepitopes, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Net_coepitopes_fig3.txt")

#########Figure 3A top --All candidates#############
####################################################

##Define objects

#HLA data
Net_BA_cast=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Net_BA_cast_fig3.txt")
Netii_BA_cast = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Netii_BA_cast_fig3.txt")
Net_coepitopes = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Net_coepitopes_fig3.txt")

#Population frequency
Freq=fread(paste0(WORKING_DIR, "HLA_freq.txt"))
Freq$Haplotype = str_replace_all(Freq$Haplotype, "DRB1_", "HLA-DRB1*")
Freq$Haplotype = str_replace_all(Freq$Haplotype, "DQA1", "DQA1*")
Freq$Haplotype = str_replace_all(Freq$Haplotype, "-DQB1", "/DQB1*")
Freq$Haplotype = str_replace_all(Freq$Haplotype, "HLA-A", "HLA-A*")
Freq$Haplotype = str_replace_all(Freq$Haplotype, "HLA-B", "HLA-B*")
Freq$Haplotype = str_replace_all(Freq$Haplotype, "HLA-C", "HLA-C*")
Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DR")] = paste0(substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DR")],1,11), ":",
                                                                      substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DR")],12,13))
Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DQ")] = paste0(substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DQ")],1,11), ":",
                                                                      substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DQ")],12,21),":",
                                                                      substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DQ")],22,23))

#For DQ reads, values were derived from a study which was biased for caucasian americans.  
#Unfortunately, paired DQ frequency data was not available for a balanced US population.
Freq$US[which(is.na(Freq$US))] = Freq$Cauc_Am[which(is.na(Freq$US))]

HLAI_nM_cut = 393.4  ###Defined from fig 2 code, see Fig S1C
HLAII_nM_cut = 220 #See Fig S1D

####Left: Predicted epitopes

#Protein dist, normalized by length
all_prots=c("orf1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF8","N","ORF10")
prot_len= c(abs(1 - 7096),abs(7097 - 8369),abs(8370 - 8644),abs(8645 - 8719),abs(8720 - 8941),abs(8942 - 9002),abs(9003 - 9123), abs(9124 - 9244),abs(9245 - 9663),abs(9664 - 9701))

Net_coepitopes = Net_coepitopes[order(Net_coepitopes$c2_start),]
prots = factor(c(unique(Net_coepitopes$Protein)), levels = all_prots)

counts_I = c()
counts_II = c()
counts_coep=c()
for(z in prots){
  counts_I = c(counts_I, length(which(Net_BA_cast$Protein == z)))#/ prot_len[which(all_prots==z)] )
  counts_II = c(counts_II, length(which(Netii_BA_cast$Protein == z)))#/ prot_len[which(all_prots==z)] )
  counts_coep = c(counts_coep, length(which(Net_coepitopes$Protein == z)))#/ prot_len[which(all_prots==z)] )
}

prot_dist = data.table(prots, counts_I, counts_II, counts_coep)
colnames(prot_dist) = c("Protein", "HLAI", "HLAII", "Coepitopes")

#Transform into polar coordinates
prot_dist <- prot_dist %>% 
  arrange(desc(Protein)) %>%
  mutate(HLAI = HLAI / sum(prot_dist$HLAI) *100) %>%
  mutate(ypos_I = cumsum(HLAI)- 0.5*HLAI )

prot_dist <- prot_dist %>% 
  arrange(desc(Protein)) %>%
  mutate(HLAII = HLAII / sum(prot_dist$HLAII) *100) %>%
  mutate(ypos_II = cumsum(HLAII)- 0.5*HLAII )

prot_dist <- prot_dist %>% 
  arrange(desc(Protein)) %>%
  mutate(Coepitopes = Coepitopes / sum(prot_dist$Coepitopes) *100) %>%
  mutate(ypos_co = cumsum(Coepitopes)- 0.5*Coepitopes )

#Plot pie charts
pie_I = ggplot(prot_dist, aes(x="", y=HLAI, fill=Protein)) +
  geom_bar(stat="identity", width=1, show.legend = F, color="red", size=1) +
  coord_polar("y", start=0)+
  labs(title="Predicted HLA-I")+
  scale_fill_manual(values=c(viridis(10)))+
  geom_label_repel(aes(y = ypos_I, label = Protein), color = "white", size=8, nudge_x = .65,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'black',
                   show.legend = F)+
  theme_void() +   
  theme(text=element_text(face="bold",size=20,colour="black")) +
  theme(legend.position="none",plot.title = element_text(hjust = 0.5))


pie_II = ggplot(prot_dist, aes(x="", y=HLAII, fill=Protein)) +
  geom_bar(stat="identity", width=1, show.legend = F, color="blue", size=1) +
  coord_polar("y", start=0)+
  labs(title="Predicted HLA-II")+
  scale_fill_manual(values=c(viridis(10)))+
  geom_label_repel(aes(y = ypos_II, label = Protein), color = "white", size=8, nudge_x = .65,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'black',
                   show.legend = F)+
  theme_void() +
  theme(text=element_text(face="bold",size=20,colour="black")) +
  theme(legend.position="none",plot.title = element_text(hjust = 0.5))

pie_co = ggplot(prot_dist, aes(x="", y=Coepitopes, fill=Protein)) +
  geom_bar(stat="identity", width=1, show.legend = F, color="purple", size=1) +
  coord_polar("y", start=0)+
  labs(title="Predicted co-epitopes")+
  scale_fill_manual(values=c(viridis(10)))+
  geom_label_repel(aes(y = ypos_co, label = Protein), color = "white", size=8, nudge_x = .65,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'black',
                   show.legend = F)+
  theme_void() +
  theme(text=element_text(face="bold",size=20,colour="black")) +
  theme(legend.position="none",plot.title = element_text(hjust = 0.5))

grid.arrange(pie_I, pie_II, pie_co, nrow=2, ncol=2)

####Right: IEDB data

#Read in SARS data
iedb_fig3 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/mhc_sars2_purified_competitive_radioactivity.csv")
iedb_fig3 = iedb_fig3[which(iedb_fig3$affinity < cutpoint),]
iedb_fig3 = iedb_fig3[which(iedb_fig3$allele %in% Freq$Haplotype),]

#Read in SARS proteome and match each IEDB peptide to find coordinates
seq =  fread(paste0(WORKING_DIR, "AA_sequence_combined.txt"), header=F)

iedb_fig3$Start = NA
iedb_fig3$End = NA
iedb_fig3$Protein = NA

#Location of each protein along proteome
prot_locs = list(
  orf1ab_pos = c(1:7096),                                                                    
  S_pos = c(7097:8369),                                                                       
  ORF3a_pos = c(8370:8644),                                                                   
  E_pos = c(8645:8719),                                                                       
  M_pos = c(8720:8941),                                                                       
  ORF6_pos = c(8942:9002),                                                                    
  ORF7a_pos = c(9003:9123),                                                                   
  ORF8_pos = c(9124:9244),                                                                    
  N_pos = c(9245:9663),                                                                       
  ORF10_pos = c(9664:9701)
)

for (y in 1:nrow(iedb_fig3)){
  if(str_detect(string = seq$V1 ,iedb_fig3$peptide[y])){
    print(y)
    iedb_fig3$Start[y] = gregexpr(text = seq$V1, pattern = iedb_fig3$peptide[y])[[1]][1]
    iedb_fig3$End[y] = (iedb_fig3$Start[y] + nchar(iedb_fig3$peptide[y]) - 1)
    iedb_fig3$Protein[y] = strsplit(names(which(unlist(prot_locs) == iedb_fig3$Start[y])),"_")[[1]][1]
  }
}

#Protein dist, normalized by length
all_prots=c("orf1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF8", "N", "ORF10")
prot_len= c(abs(1 - 7096),abs(7097 - 8369),abs(8370 - 8644),abs(8645 - 8719),abs(8720 - 8941),abs(8942 - 9002),abs(9003 - 9123), abs(9124 - 9244),abs(9245 - 9663),abs(9664 - 9701))

#Set proteins as factor
prots = factor(c(unique(Net_coepitopes$Protein)), levels = all_prots)

#Split data into HLA-I and II
iedb_fig3_I = iedb_fig3[which(iedb_fig3$mhc_class ==  "I"),]
iedb_fig3_II = iedb_fig3[which(iedb_fig3$mhc_class ==  "II"),]

#Cast matrices by unique peptide
iedb_fig3_I = dcast.data.table(iedb_fig3_I, peptide+Protein+Start+End~allele, value.var = "affinity")
iedb_fig3_II = dcast.data.table(iedb_fig3_II, peptide+Protein+Start+End~allele, value.var = "affinity")

#Count frequency of each protein among IEDB data
counts_iedb_I = c()
counts_iedb_II = c()

for(z in prots){
  counts_iedb_I = c(counts_iedb_I, length(which(iedb_fig3_I$Protein == z)))#/ prot_len[which(all_prots==z)] )
  counts_iedb_II = c(counts_iedb_II, length(which(iedb_fig3_II$Protein == z)))#/ prot_len[which(all_prots==z)] )
}

prot_dist = data.table(prots, counts_iedb_I, counts_iedb_II)
colnames(prot_dist)[1] = "Protein"

#Convert to polar coordinates
prot_dist <- prot_dist %>% 
  arrange(desc(Protein)) %>%
  mutate(counts_iedb_I = counts_iedb_I / sum(prot_dist$counts_iedb_I) *100) %>%
  mutate(ypos_iedb_I = cumsum(counts_iedb_I)- 0.5*counts_iedb_I )
prot_dist <- prot_dist %>% 
  arrange(desc(Protein)) %>%
  mutate(counts_iedb_II = counts_iedb_II / sum(prot_dist$counts_iedb_II) *100) %>%
  mutate(ypos_iedb_II = cumsum(counts_iedb_II)- 0.5*counts_iedb_II )

iedb_I=ggplot(prot_dist, aes(x="", y=counts_iedb_I, fill=Protein)) +
  geom_bar(stat="identity", width=1, show.legend = F, color = "red", size=1) +
  coord_polar("y", start=0)+
  labs(title="IEDB HLA-I epitopes")+
  scale_fill_manual(values=c(viridis(10)))+
  geom_label_repel(aes(y = ypos_iedb_I, label = Protein), color = "white", size=8, nudge_x = .65,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'black',
                   show.legend = F)+
  theme_void() +   
  theme(text=element_text(face="bold",size=20,colour="black")) +
  theme(legend.position="none",plot.title = element_text(hjust = 0.5))

iedb_II= ggplot(prot_dist, aes(x="", y=counts_iedb_II, fill=Protein)) +
  geom_bar(stat="identity", width=1, show.legend = F, color="blue", size = 1) +
  coord_polar("y", start=0)+
  labs(title="IEDB HLA-II epitopes")+
  scale_fill_manual(values=c(viridis(10)))+
  geom_label_repel(aes(y = ypos_iedb_II, label = Protein), color = "white", size=8, nudge_x = .65,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'black',
                   show.legend = F)+
  theme_void() +   
  theme(text=element_text(face="bold",size=20,colour="black")) +
  theme(legend.position="none",plot.title = element_text(hjust = 0.5))


grid.arrange(iedb_I, iedb_II, nrow=1, ncol=2)


####Combine predicted epitopes with IEDB epitopes, incorporate entropy
ent = fread(paste0(WORKING_DIR, "entropy_7882.txt"))

#All peptides from the IEDB and netMHC(II)pan sets
all_I = unique(c(iedb_fig3_I$peptide, Net_BA_cast$Peptide))
all_II = unique(c(iedb_fig3_II$peptide, Netii_BA_cast$Peptide))

#For each peptide, determine binding haplotypes, population frequencies, protein, coordinates, entropy scores
Pop_freq_I_iedb = foreach(n=1:length(all_I), .combine=rbind) %dopar% {
  iedb = which(iedb_fig3_I$peptide == all_I[n])
  pred = which(Net_BA_cast$Peptide == all_I[n])
  
  if(length(iedb)>0 & length(pred)>0){
    Binding_haplotype=paste0(unique(c(colnames(iedb_fig3_I[iedb,c(5:13)])[which(iedb_fig3_I[iedb,c(5:13)] < cutpoint)],
                                      colnames(Net_BA_cast[pred,c(4:19)])[which(Net_BA_cast[pred,c(4:19)] < HLAI_nM_cut)])), collapse = ",")
    Frequency=1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% strsplit(Binding_haplotype,",")[[1]] )]/100,2)))
    Start= iedb_fig3_I$Start[iedb]
    Protein = iedb_fig3_I$Protein[iedb]
  }else if(length(iedb)>0 & length(pred)==0){
    Binding_haplotype=paste0(colnames(iedb_fig3_I[iedb,c(5:13)])[which(iedb_fig3_I[iedb,c(5:13)] < cutpoint)], collapse = ",")
    Frequency=1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% strsplit(Binding_haplotype,",")[[1]] )]/100,2)))
    Start= iedb_fig3_I$Start[iedb]
    Protein = iedb_fig3_I$Protein[iedb]
  }else if(length(iedb)==0 & length(pred)>0){
    Binding_haplotype=paste0(colnames(Net_BA_cast[pred,c(4:19)])[which(Net_BA_cast[pred,c(4:19)] < HLAI_nM_cut)], collapse = ",")
    Frequency=1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% strsplit(Binding_haplotype,",")[[1]] )]/100,2)))
    Start= Net_BA_cast$Start[pred]
    Protein = Net_BA_cast$Protein[pred]
  }
  End=Start+nchar(all_I[n])-1
  Low_entropy = ifelse(max(ent$V2[Start:End]) > 0.1, 0, 1)
  data.frame(all_I[n],Binding_haplotype,Frequency,Protein,Start,End,Low_entropy)
}


Pop_freq_II_iedb = foreach(n=1:length(all_II), .combine=rbind) %dopar% {
  iedb = which(iedb_fig3_II$peptide == all_II[n])
  pred = which(Netii_BA_cast$Peptide == all_II[n])
  
  if(length(iedb)>0 & length(pred)>0){
    Binding_haplotype=paste0(unique(c(colnames(iedb_fig3_II[iedb,c(5)])[which(iedb_fig3_II[iedb,c(5)] < cutpoint)],
                                      colnames(Netii_BA_cast[pred,c(4:18)])[which(Netii_BA_cast[pred,c(4:18)] < HLAII_nM_cut)])), collapse = ",")
    Frequency=1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% strsplit(Binding_haplotype,",")[[1]] )]/100,2)))
    Start= iedb_fig3_II$Start[iedb]
    Protein = iedb_fig3_II$Protein[iedb]
  }else if(length(iedb)>0 & length(pred)==0){
    Binding_haplotype=paste0(colnames(iedb_fig3_II[iedb,c(5)])[which(iedb_fig3_II[iedb,c(5)] < cutpoint)], collapse = ",")
    Frequency=1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% strsplit(Binding_haplotype,",")[[1]] )]/100,2)))
    Start= iedb_fig3_II$Start[iedb]
    Protein = iedb_fig3_II$Protein[iedb]
  }else if(length(iedb)==0 & length(pred)>0){
    Binding_haplotype=paste0(colnames(Netii_BA_cast[pred,c(4:18)])[which(Netii_BA_cast[pred,c(4:18)] < HLAII_nM_cut)], collapse = ",")
    Frequency=1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% strsplit(Binding_haplotype,",")[[1]] )]/100,2)))
    Start= Netii_BA_cast$Start[pred]
    Protein = Netii_BA_cast$Protein[pred]
  }
  End=Start+nchar(all_II[n])-1
  Low_entropy = ifelse(max(ent$V2[Start:End]) > 0.1, 0, 1)
  data.frame(all_II[n],Binding_haplotype,Frequency,Protein,Start,End,Low_entropy)
}

#Correct co-epitope set to account for allele coverage from IEDB set
Net_coepitopes$Low_entropy=NA
counter = 0
percent = 0
for (z in 1:nrow(Net_coepitopes)){
  counter = counter +1
  if(counter > nrow(Net_coepitopes)/100){
    counter = 0
    percent = percent + 1
    print(paste0(percent,"%"))
  }
  Net_coepitopes$Low_entropy[z] = ifelse(max(ent$V2[Net_coepitopes$c2_start[z]:(Net_coepitopes$c2_start[z]+nchar(Net_coepitopes[z]$c2_peptide)-1)]) > 0.1, 0, 1)
  
  if(Net_coepitopes$c1_peptide[z] %in% iedb_fig3_I$peptide){
    if(round(Net_coepitopes$c1_freq[z]/100,3) != round(Pop_freq_I_iedb$Frequency[which(Pop_freq_I_iedb$all_I.n. == Net_coepitopes$c1_peptide[z])],3)){
      print(paste0("HLAI ",z, ": ", round(Net_coepitopes$c1_freq[z]/100,3), "->", round(Pop_freq_I_iedb$Frequency[which(Pop_freq_I_iedb$all_I.n. == Net_coepitopes$c1_peptide[z])],3)))
      Net_coepitopes$c1_haplotypes[z] = Pop_freq_I_iedb$Binding_haplotype[which(Pop_freq_I_iedb$all_I.n. == Net_coepitopes$c1_peptide[z])]
      ####If there is a difference in the binding haplotypes, caused by IEDB coverage additions, then correct the coepitope table for c/w HLAI table
    }
  }
  if(Net_coepitopes$c2_peptide[z] %in% iedb_fig3_II$peptide){
    if(round(Net_coepitopes$c2_freq[z]/100, 3) != round(Pop_freq_II_iedb$Frequency[which(Pop_freq_II_iedb$all_II.n. == Net_coepitopes$c2_peptide[z])],3)){
      print(paste0("HLAII ",z, ": ", round(Net_coepitopes$c2_freq[z]/100,3), "->", round(Pop_freq_II_iedb$Frequency[which(Pop_freq_II_iedb$all_II.n. == Net_coepitopes$c2_peptide[z])],3)))
      Net_coepitopes$c2_haplotypes[z] = Pop_freq_II_iedb$Binding_haplotype[which(Pop_freq_II_iedb$all_II.n. == Net_coepitopes$c2_peptide[z])]
      ####If there is a difference in the binding haplotypes, caused by IEDB coverage additions, then correct the coepitope table for c/w HLAII table
      
    }
  }
}

fwrite(Pop_freq_I_iedb, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAI_vaccine_candidates_prefilter.txt")
fwrite(Pop_freq_II_iedb, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAII_vaccine_candidates_prefilter.txt")
fwrite(Net_coepitopes, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Coepitopes_vaccine_candidates_prefilter.txt")


#################Figure 3 middle -- Filtering###################
################################################################

#Define function to plot allele proportions at each filtering step, with counts of class I/II and coepitopes
bar_chart = function(mat, mat2, matco, size){
  all_haps = c("HLA-A*01:01", "HLA-A*02:01", "HLA-A*03:01", "HLA-A*11:01", "HLA-A*24:02", "HLA-B*07:02", "HLA-B*08:01", 
               "HLA-B*35:01", "HLA-B*44:02", "HLA-B*44:03", "HLA-C*03:04", "HLA-C*04:01", "HLA-C*05:01", "HLA-C*06:02", "HLA-C*07:01", "HLA-C*07:02")
  all_haps = all_haps[order(all_haps)]
  
  all_haps2= c("HLA-DRB1*01:01", "HLA-DRB1*03:01", "HLA-DRB1*04:01", "HLA-DRB1*07:01", "HLA-DRB1*11:01",
               "HLA-DRB1*13:01", "HLA-DRB1*15:01", "HLA-DQA1*01:02/DQB1*06:02", "HLA-DQA1*05:01/DQB1*02:01", "HLA-DQA1*02:01/DQB1*02:02",
               "HLA-DQA1*05:05/DQB1*03:01", "HLA-DQA1*01:01/DQB1*05:01", "HLA-DQA1*03:01/DQB1*03:02", "HLA-DQA1*03:03/DQB1*03:01", "HLA-DQA1*01:03/DQB1*06:03")
  all_haps2 = all_haps2[order(all_haps2)]
  
  count = c()
  for(h in all_haps){
    count = c(count, length(which(unlist(strsplit(mat$Binding_haplotype,",")) == h))/length(unlist(strsplit(mat$Binding_haplotype,",")))) 
  }
  count2 = c()
  for(h in all_haps2){
    count2 = c(count2, length(which(unlist(strsplit(mat2$Binding_haplotype,",")) == h))/length(unlist(strsplit(mat2$Binding_haplotype,","))))
  }
  
  dat = data.table(all_haps, count, "HLA-I")
  dat2 = data.table(all_haps2, count2, "HLA-II")
  dat_text = data.table(nrow(mat),nrow(matco), nrow(mat2))
  
  HLAI_col= c(brewer.pal(7, 'Blues')[3:7],
              brewer.pal(7, 'Reds')[3:7],
              brewer.pal(8, 'Greens')[3:8])
  HLAII_col = c(brewer.pal(9, 'Oranges')[2:9],
                brewer.pal(9, 'Purples')[3:9])
  
  dat_all = rbind(dat, dat2, use.names=F)
  dat_all$all_haps = factor(dat_all$all_haps, levels = c(all_haps, all_haps2))
  dat_all$V3 = factor(dat_all$V3, levels = c("HLA-II", "HLA-I"))
  
  p=ggplot(data = dat_all, aes(x=V3,y=count,fill=all_haps))+geom_bar(stat='identity')+
    scale_fill_manual(values = c(HLAI_col, HLAII_col),
                      drop = F, guide = guide_legend(override.aes = list(color = "white")))+
    
    theme_classic()+theme(axis.line=element_blank(),axis.text.x=element_blank(),
                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank(),legend.position="none",
                          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                          panel.grid.minor=element_blank(),plot.background=element_blank()) +
    annotate("label", y=c(.1,.5,.9), x= 1.5, label = c(dat_text$V1, dat_text$V2, dat_text$V3), color = c("red", "purple", "blue"), 
             size = size, vjust = "center", angle = 0)
  
  return(p + coord_flip())
}


##No filter
Pop_freq_I_iedb = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAI_vaccine_candidates_prefilter.txt")
Pop_freq_II_iedb = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAII_vaccine_candidates_prefilter.txt")
Net_coepitopes = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Coepitopes_vaccine_candidates_prefilter.txt")

##Filter by tetramer models
MHCI_tet =fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAI_all_glm.txt")
MHCII_tet =fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAII_all_glm.txt")

MHCI_tet_median = median(MHCI_tet$GLM_score[which((MHCI_tet$Binding_affinity*50000) < 393.4)])
MHCII_tet_median = median(MHCII_tet$GLM_score[which((MHCII_tet$v3.2_Binding_affinity) < 220)])

#Read in population frequencies of each allele, format allele names to be consistent
Freq=fread(paste0(WORKING_DIR, "HLA_freq.txt"))
Freq$Haplotype = str_replace_all(Freq$Haplotype, "DRB1_", "HLA-DRB1*")
Freq$Haplotype = str_replace_all(Freq$Haplotype, "DQA1", "DQA1*")
Freq$Haplotype = str_replace_all(Freq$Haplotype, "-DQB1", "/DQB1*")
Freq$Haplotype = str_replace_all(Freq$Haplotype, "HLA-A", "HLA-A*")
Freq$Haplotype = str_replace_all(Freq$Haplotype, "HLA-B", "HLA-B*")
Freq$Haplotype = str_replace_all(Freq$Haplotype, "HLA-C", "HLA-C*")
Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DR")] = paste0(substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DR")],1,11), ":",
                                                                      substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DR")],12,13))
Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DQ")] = paste0(substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DQ")],1,11), ":",
                                                                      substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DQ")],12,21),":",
                                                                      substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DQ")],22,23))
Freq$US[which(is.na(Freq$US))] = Freq$Cauc_Am[which(is.na(Freq$US))]

#Add in GLM-filtered population frequency and immunogenic haplotypes for class I, II, and coepitopes

#Class I
filt_I = Pop_freq_I_iedb
filt_I$Haplotype_GLM = NA
filt_I$Immunogenic_haplotypes = NA
filt_I$Immunogenic_frequency = NA
for (n in 1:nrow(filt_I)){
  print(n)
  haps = unlist(strsplit(filt_I$Binding_haplotype[n],",")) #Haplotypes defined for each peptide
  
  #Grab the GLM scores for each epitope
  GLM = c()
  for(h in haps){
    GLM = c(GLM, MHCI_tet$GLM_score[which((MHCI_tet$Peptide == filt_I$all_I.n.[n]) & (MHCI_tet$Haplotype == h))])
  }
  filt_I$Haplotype_GLM[n] = paste0(GLM, collapse = ",") #Collapse GLM scores to comma-delimited
  filt_I$Immunogenic_haplotypes[n] = paste0(haps[which(GLM > MHCI_tet_median)], collapse = ",") #Lists alleles for which GLM score meets median cutoff
  filt_I$Immunogenic_frequency[n] = 1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% strsplit(filt_I$Immunogenic_haplotypes[n],",")[[1]])]/100,2))) #Calculate population frequency for immunogenic alleles
  
}
fwrite(filt_I, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAI_vaccine_candidates_prefilter_glm.txt")

#Class II, same strategy as above
filt_II = Pop_freq_II_iedb
filt_II$Haplotype_GLM = NA
filt_II$Immunogenic_haplotypes = NA
filt_II$Immunogenic_frequency = NA
for (n in 1:nrow(filt_II)){
  print(n)
  haps = unlist(strsplit(filt_II$Binding_haplotype[n],","))
  
  GLM = c()
  for(h in haps){
    GLM = c(GLM, MHCII_tet$GLM_score[which((MHCII_tet$Peptide == filt_II$all_II.n.[n]) & (MHCII_tet$Haplotype == h))])
  }
  filt_II$Haplotype_GLM[n] = paste0(GLM, collapse = ',')
  filt_II$Immunogenic_haplotypes[n] = paste0(haps[which(GLM > MHCII_tet_median)], collapse = ",")
  filt_II$Immunogenic_frequency[n] = 1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% strsplit(filt_II$Immunogenic_haplotypes[n],",")[[1]])]/100,2)))
  
}
fwrite(filt_II, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAII_vaccine_candidates_prefilter_glm.txt")

#Co-epitopes, same strategy as above
filt_co = Net_coepitopes
filt_co$c1_GLM = NA
filt_co$c2_GLM = NA
filt_co$c1_immunogenic_haplotypes = NA
filt_co$c2_immunogenic_haplotypes = NA
filt_co$c1_Immunogenic_frequency = NA
filt_co$c2_Immunogenic_frequency = NA

for (n in 1:nrow(filt_co)){
  print(n)
  haps1 = unlist(strsplit(filt_co$c1_haplotypes[n],","))
  haps2 = unlist(strsplit(filt_co$c2_haplotypes[n],","))
  
  GLM1 = c()
  for(h in haps1){
    GLM1 = c(GLM1, MHCI_tet$GLM_score[which((MHCI_tet$Peptide == filt_co$c1_peptide[n]) & (MHCI_tet$Haplotype == h))])
  }
  filt_co$c1_GLM[n] = paste0(GLM1, collapse = ",")
  filt_co$c1_immunogenic_haplotypes[n] = paste0(haps1[which(GLM1 > MHCI_tet_median)], collapse = ",")
  filt_co$c1_Immunogenic_frequency[n] = 1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% strsplit(filt_co$c1_immunogenic_haplotypes[n],",")[[1]])]/100,2)))
  
  GLM2 = c()
  for(h in haps2){
    GLM2 = c(GLM2, MHCII_tet$GLM_score[which((MHCII_tet$Peptide == filt_co$c2_peptide[n]) & (MHCII_tet$Haplotype == h))])
  }
  filt_co$c2_GLM[n] = paste0(GLM2, collapse = ",")
  filt_co$c2_immunogenic_haplotypes[n] = paste0(haps2[which(GLM2 > MHCII_tet_median)], collapse = ",")
  filt_co$c2_Immunogenic_frequency[n] = 1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% strsplit(filt_co$c2_immunogenic_haplotypes[n],",")[[1]])]/100,2)))
  
}
fwrite(filt_co, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Coepitopes_vaccine_candidates_prefilter_glm.txt")


###Pre filter counts
filt_I = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAI_vaccine_candidates_prefilter_glm.txt")
filt_II = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAII_vaccine_candidates_prefilter_glm.txt")
filt_co = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Coepitopes_vaccine_candidates_prefilter_glm.txt")

t1 = bar_chart(filt_I, filt_II, filt_co, 30)

##Filter by immunogenicity scores
filt_I = filt_I[which(filt_I$Immunogenic_haplotypes != "")]
filt_II = filt_II[which(filt_II$Immunogenic_haplotypes != "")]
filt_co = filt_co[which((filt_co$c1_immunogenic_haplotypes != "") & (filt_co$c2_immunogenic_haplotypes != ""))]

t2 = bar_chart(filt_I, filt_II, filt_co, 26)

##Filter by entropy
filt_I = filt_I[which(filt_I$Low_entropy==1),]
filt_II = filt_II[which(filt_II$Low_entropy==1),]
filt_co = filt_co[which(filt_co$Low_entropy==1),]

t3 = bar_chart(filt_I, filt_II, filt_co, 22)


#Filter by protein -- run for figure 5, do not run for supplemental containing all proteins

###Length normalized
filt_I = filt_I[which(filt_I$Protein %in% c("M", "S", "N")),]
filt_II = filt_II[which(filt_II$Protein %in% c("M", "S", "N")),]
filt_co = filt_co[which(filt_co$Protein %in% c("M", "S", "N")),]

t4 = bar_chart(filt_I, filt_II, filt_co, 18)


##Add in murine coverage, don't filter
mouse1 = fread(paste0(WORKING_DIR, "COVID_murine_NetMHCpan_unfilt.xls"))
mouse1 = mouse1[which(mouse1$NB>0),]
mouse2 = fread(paste0(WORKING_DIR, "COVID_murine_NetMHCIIpan_unfilt.xls"))
mouse2 = mouse2[which(mouse2$NB>0),]

filt_I$Murine = NA
filt_II$Murine = NA
#filt_co$Murine_I = NA
#filt_co$Murine_II = NA

for (n in 1:nrow(filt_I)){
  if(filt_I$all_I.n.[n] %in% mouse1$Peptide){
    hit = which(mouse1$Peptide == filt_I$all_I.n.[n])
    if((length(which(mouse1[hit, c(8, 18)] < 2)) > 0) & (length(which(mouse1[hit, c(13, 23)] < 2)) > 0)){ #Looking for if either b or d MHC-I scores are in the top 2nd percentiles
      filt_I$Murine[n] = "B/D"
    }else if(length(which(mouse1[hit, c(8, 18)] < 2)) > 0){ #Looking for if either b MHC-I scores are in the top 2nd percentiles
      filt_I$Murine[n] = "B"
    }else if(length(which(mouse1[hit, c(13, 23)] < 2)) > 0){ #Looking for if either d MHC-I scores are in the top 2nd percentiles
      filt_I$Murine[n] = "D"
    }
  }
}
for (n in 1:nrow(filt_II)){
  if(filt_II$all_II.n.[n] %in% mouse2$Peptide){
    hit = which(mouse2$Peptide == filt_II$all_II.n.[n])
    if((length(which(mouse2[hit, 6] < 10)) > 0) & (length(which(mouse2[hit, 9] < 10)) > 0)){ #Same as above, but 10% percentile or better for MHC-II
      filt_II$Murine[n] = "B/D"
    }else if(length(which(mouse2[hit, 6] < 10)) > 0){
      filt_II$Murine[n] = "B"
    }else if(length(which(mouse2[hit, 9] < 10)) > 0){
      filt_II$Murine[n] = "D"
    }
  }
}

###########Figure 3 middle - Plot funnel###############
#######################################################

ggdraw() +
  draw_plot(t1, x = 0, y =.75, width = 1, height = .25)+
  
  draw_plot(t2, x = .1, y =.5, width = .8, height = .25)+
  
  draw_plot(t3, x = .2, y =.25, width = .6, height = .25)+
  
  draw_plot(t4, x = .3, y =0, width = .4, height = .25)



#######Figure 3 bottom - Centipede plot###############
######################################################


#Combine matrices
filt_I$Immunogenic_frequency = filt_I$Immunogenic_frequency*100
filt_II$Immunogenic_frequency = filt_II$Immunogenic_frequency*100

filt_I$Group = "HLAI"
filt_II$Group = "HLAII"
filt_co$Group = "Coepitope"
filt_co$End = filt_co$Start+nchar(filt_co$c2_peptide)-1

Master_tab=rbind(filt_I,filt_II, use.names=F)
colnames(Master_tab)[1] = "Peptide"

#Coep_tab = Master_tab[which(Master_tab$Group =="Coepitope"),]
#Master_tab = Master_tab[which(Master_tab$Group !="Coepitope"),]

###Find overlapping class I and II peptides
g1 <-  GRanges(seqnames="COVID",
               IRanges(start=Master_tab[which(Master_tab$Group == "HLAI")]$Start,
                       end=Master_tab[which(Master_tab$Group == "HLAI")]$End), 
               I_peptide = Master_tab[which(Master_tab$Group == "HLAI")]$Peptide, 
               Frequency_I= Master_tab[which(Master_tab$Group == "HLAI")]$Immunogenic_frequency)

g2 <-  GRanges(seqnames="COVID",
               IRanges(start=Master_tab[which(Master_tab$Group == "HLAII")]$Start,
                       end=Master_tab[which(Master_tab$Group == "HLAII")]$End), 
               II_peptide = Master_tab[which(Master_tab$Group == "HLAII")]$Peptide, 
               Frequency_II= Master_tab[which(Master_tab$Group == "HLAII")]$Immunogenic_frequency)

HLA_I_II=as.data.table(mergeByOverlaps(g1,g2))  #Overlapping regions for HLA-I and HLA-II epitopes

###Add overlap %

HLA_I_II$Overlap = NA
for(z in 1:nrow(HLA_I_II)){
  c1 = c(HLA_I_II$g1.start[z]:HLA_I_II$g1.end[z])  #Positions for HLA-I epitope
  c2 = c(HLA_I_II$g2.start[z]:HLA_I_II$g2.end[z]) #Positions for HLA-II epitope
  HLA_I_II$Overlap[z] = 100*(length(which(c1 %in% c2))/length(c1))  #Percentage of HLA-I epitope which is overlapped by HLA-II
  
}
HLA_I_II$Frequency = (HLA_I_II$Frequency_I/100)*(HLA_I_II$Frequency_II/100)

###Reduce class I and II epitopes to unique, non-overlapping regions
g1_merged <-  GRanges(seqnames="COVID",
                      IRanges(start= HLA_I_II$g1.start,
                              end=HLA_I_II$g1.end), 
                      I_peptide = HLA_I_II$g1.I_peptide, 
                      Frequency_I= HLA_I_II$Frequency_I)
g1_red = reduce(g1_merged)  #Reduce all HLA-I overlaps
g2_merged <-  GRanges(seqnames="COVID",
                      IRanges(start= HLA_I_II$g2.start,
                              end=HLA_I_II$g2.end), 
                      II_peptide = HLA_I_II$g2.II_peptide, 
                      Frequency_II= HLA_I_II$Frequency_II)

g2_red = reduce(g2_merged)  #Reduce all HLA-II overlaps             

#Look for overlap between class I and II regions
red_merged=as.data.table(mergeByOverlaps(g1_red, g2_red))

###For each merged region, find the best HLA-I and II epitope based on population frequency
merge_tab = data.table()
for(z in 1:nrow(red_merged)){
  print(z)
  rows = which((HLA_I_II$g1.start >= red_merged$g1_red.start[z]) & (HLA_I_II$g1.end <= red_merged$g1_red.end[z])) #Find epitopes within overlapping regions
  
  best_frequency = max((HLA_I_II$g1.Frequency_I[rows]/100) * (HLA_I_II$g2.Frequency_II[rows]/100))  #Highest overall pop. freq
  best_row = rows[which(((HLA_I_II$g1.Frequency_I[rows]/100) * (HLA_I_II$g2.Frequency_II[rows]/100)) == best_frequency)] #Which rows correspond to those w/high freq
  
  overlap = c()
  for(r in best_row){
    c1 = c(HLA_I_II$g1.start[r]:HLA_I_II$g1.end[r])  #HLA-I region
    c2 = c(HLA_I_II$g2.start[r]:HLA_I_II$g2.end[r]) #HLA-II region
    overlap = c(overlap, length(which(c1 %in% c2))/length(c1) ) #Percentage overlap
  }
  
  best_overlap = max(overlap) #Overlap %
  best_c1 = paste0(unique(HLA_I_II$g1.I_peptide[best_row]), collapse = ", ") #HLA-I epitopes
  best_c2 = paste0(unique(HLA_I_II$g2.II_peptide[best_row]), collapse = ", ") #HLA-II epitopes
  best_start = min(HLA_I_II$g2.start[best_row])
  
  merge_tab = rbind(merge_tab, data.table(best_c1, best_c2, best_start, best_overlap, best_frequency))
}

#Convert HLA-II frequencies to negative values for plotting, convert frequencies to percent scale
Master_tab$Immunogenic_frequency[which(Master_tab$Group =="HLAII")] = -1*Master_tab$Immunogenic_frequency[which(Master_tab$Group =="HLAII")] 
HLA_I_II$Frequency = HLA_I_II$Frequency*100
merge_tab$best_frequency = 100 * merge_tab$best_frequency

####Make plot
y_ll=-90  #Lower limit of y axis

centipede = 
  ggplot(data=Master_tab) + 
  
  geom_hline(yintercept = c(-80,-70,-60,-50,50,60,70,80), color = "white")+

  geom_segment(aes(x=Start, xend=Start, y=ifelse(Master_tab$Group=="HLAI",10,-10), yend=Immunogenic_frequency, color = Group), 
               alpha = .1, size=.2, show.legend = F) +
  scale_color_manual(values = c("red", "blue"), name="HLA", guide=F)+
  ggnewscale::new_scale("color")+
  
  geom_point(aes(x=Start, y=Immunogenic_frequency, fill=Group, size = abs(Immunogenic_frequency)), color='black',  show.legend = T, shape=21) +
  scale_y_continuous(breaks = seq(-80,80,10), labels = abs)+
  scale_color_manual(values = alpha("black",1), guide = F)+
  scale_fill_manual(values = c(alpha("red",.1), alpha("blue",.1)), name="HLA", guide=F)+
  scale_size_continuous(name = "HLA-I/II\nPop. freq.", limits = c(9,90), range = c(1,10))+
  ggnewscale::new_scale("color")+
  ggnewscale::new_scale("size")+
  ggnewscale::new_scale("alpha")+
  ggnewscale::new_scale("fill")+
  
  # ###Try using all merged data
  # geom_point(data=HLA_I_II,aes(x=g2.start, y=0, size = Frequency , fill = Overlap), alpha = .1,
  #            show.legend = T,  position=position_jitter(width=0, height=7.5), shape = 23)+
  # scale_fill_viridis_c(direction = 1, name = "% overlap", option = "plasma")+
  # scale_color_manual(values = alpha("purple",1), guide = F)+
  # scale_size_continuous(name = "Overlaps\nPop. freq.", limits = c(1,55), range = c(1,10))+
  # #scale_alpha_continuous(range = c(.25,.75), guide = F)+
  
  ###Try using merged_tab
  geom_point(data=merge_tab,aes(x=best_start, y=0, size = best_frequency , fill = best_overlap),
           alpha = 0.6, show.legend = T, position=position_jitter(width=0, height=7.5), shape = 23)+
  scale_fill_viridis_c(direction = 1, name = "% overlap", option = "plasma")+
  scale_color_manual(values = alpha("black",1), guide = F)+
  scale_size_continuous(name = "Overlaps\nPop. freq.", limits = c(1,55), range = c(1,10))+
  
  annotate("text", x = c(-400,-400,-400), y = c(-65,0,65), label = c("HLA-II", "Overlap", "HLA-I"), angle =90, size = 7)+
  theme(text=element_text(face="bold",size=30,colour="black")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  
  geom_rect(aes(xmin=7096, xmax=8368, ymin=(y_ll), ymax=(y_ll-12)), fill=viridis(10)[2], color="black", size=0.1) +
  geom_text(aes(x=7096+(8368-7096)/2, y=(y_ll-17), label="S", angle=0), size=10, color = "black") + 
  geom_text(aes(x=7096+200, y=(y_ll), label="-200", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=7096+500, y=(y_ll), label="-500", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=7096+800, y=(y_ll), label="-800", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=7096+1100, y=(y_ll), label="-1100", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  
  geom_rect(aes(xmin=8719, xmax=8940, ymin=(y_ll-0), ymax=(y_ll-12)), fill=viridis(10)[5], color="black", size=0.1) +
  geom_text(aes(x=8719+(8940-8719)/2, y=(y_ll-17)), label="M", angle=0, size=10, color = "black") +
  geom_text(aes(x=8719+50, y=(y_ll), label="-50", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=8719+100, y=(y_ll), label="-100", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=8719+150, y=(y_ll), label="-150", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=8719+200, y=(y_ll), label="-200", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  
  geom_rect(aes(xmin=9244, xmax=9662, ymin=(y_ll-0), ymax=(y_ll-12)), fill=viridis(10)[9], color="black", size=0.1) +
  geom_text(aes(x=9244+(9662-9244)/2, y=(y_ll-17)), label="N", angle=0, size=10, color = "black") +
  geom_text(aes(x=9244+100, y=(y_ll), label="-100", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=9244+200, y=(y_ll), label="-200", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=9244+300, y=(y_ll), label="-300", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) 


c_S = centipede + scale_x_continuous(limits = c(7096,8369)) + theme(legend.position = "none") +  labs(y = "Population frequency", x="")
c_M = centipede + scale_x_continuous(limits = c(8719,8941)) + theme(legend.position = "none") + labs(y="",x="")
c_N = centipede + scale_x_continuous(limits = c(9244,9663)) + labs(y="",x="")
  
#Note, this will throw warnings because of areas getting cut off in the x axis, but the data plotted should be correct
png("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Figure_3_bottom.png", width = 20, height = 8.3, units = "in", res=300)
ggarrange(c_S, c_M, c_N, nrow = 1, ncol = 3,common.legend = TRUE, legend="right")
dev.off()



########################Figure 4A####################
#####################################################

#Curated Bc data from 4 array sources
Bc = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/bcell/linear-bcell-epitopes-SARS2-S.csv")
Bc$Source[which(Bc$Source == "ReScan")] = "Zamecnik 2020"

Bc = Bc[order(Bc$Start),]
Bc= Bc[order(Bc$Source),]
Bc$Segment = 1

#Define the y-axis for each source
for(s in unique(Bc$Source)){
  Bc$Segment[which(Bc$Source == s)] = which(s == rev(unique(Bc$Source)))
}

Fig4A=ggplot(data=Bc) + 
  geom_rect(aes(xmin=Start, xmax=End, ymin=Segment-.1875, ymax=Segment+.1875, fill=Source), color = "black")+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  scale_fill_viridis_d(option = "plasma", end = .75, alpha = .25)+
  geom_hline(yintercept = 4.5, size=.5, color="white")+
  ggnewscale::new_scale("color")+
  
  geom_text_repel(aes(x=Start, y=Segment, label = Start), direction = "y", nudge_y = .5,
                  segment.colour = "black", segment.alpha = .5)+
  
  labs(y='',x="Position across SARS-CoV-2 spike protein") +
  
  geom_rect(aes(xmin=(1), xmax=(1273), ymin=.5, ymax=(0)), fill="grey", color="black", size=0.1) +
  
  geom_rect(aes(xmin=(319), xmax=(541), ymin=.5, ymax=(0)), fill="yellow", color="black", size=0.1) +
  geom_text(aes(x=(319+541)/2, y=(-.1), label="RBD", angle=0), size=5) +
  
  geom_rect(aes(xmin=788, xmax=806, ymin=.5, ymax=(0)), fill="skyblue", color="black", size=0.1) +
  geom_text(aes(x=(788+806)/2, y=(-.1), label="FP", angle=0), size=4.5) +
  
  geom_rect(aes(xmin=(437), xmax=(508), ymin=.5, ymax=(0)), fill="tomato", color="black", size=0.1) +
  geom_text(aes(x=(437+508)/2, y=(.25), label="RBM", angle=0), size=5) +
  
  geom_rect(aes(xmin=912, xmax=984, ymin=.5, ymax=(0)), fill="deep pink", color="black", size=0.1) +
  geom_text(aes(x=(912+984)/2, y=(-.1), label="HR1", angle=0), size=4.5) +
  
  geom_rect(aes(xmin=1163, xmax=1213, ymin=.5, ymax=(0)), fill="darkviolet", color="black", size=0.1) +
  geom_text(aes(x=(1163+1213)/2, y=(-.1), label="HR2", angle=0), size=4.5) +
  
  geom_rect(aes(xmin=685, xmax=686, ymin=.5, ymax=(0)), fill="black", color="black", size=0.1) +
  geom_text(aes(x=685, y=(-.1), label="S1/S2", angle=0), size=4.5) +
    
  scale_y_continuous(limits = c(-.1,6))+
  scale_x_continuous(limits = c(1,1273), breaks = c(1, 200, 400, 600, 800, 1000, 1273))+
  theme(axis.title.y=element_blank(),
        legend.title = element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.text=element_text(size=25),
  axis.text=element_text(size=20),
  axis.title=element_text(size=25,face="bold"))

png("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Figure_4A.png", width = 12.5, height = 8, units = "in", res=300)
Fig4A
dev.off()



###Funnel, built in PPT###
#Curated epitopes: 58 
#Accessible regions: 27 
#Combine overlapping: 22
#Glycosite filter: 18
#Polymorphism filter: 15
#Near functional domains (RBD, FP, HR1, HR2): 8
#4mer or longer: 4



#####Figure 5 -- filtering Tc epitopes#################################################
##Note that this section is dependent on above code for figure 3 starting line 2545####
#######################################################################################

##Reading in frequency data, formatting
Freq=fread(paste0(WORKING_DIR, "HLA_freq.txt"))
Freq$Haplotype = str_replace_all(Freq$Haplotype, "DRB1_", "HLA-DRB1*")
Freq$Haplotype = str_replace_all(Freq$Haplotype, "DQA1", "DQA1*")
Freq$Haplotype = str_replace_all(Freq$Haplotype, "-DQB1", "/DQB1*")
Freq$Haplotype = str_replace_all(Freq$Haplotype, "HLA-A", "HLA-A*")
Freq$Haplotype = str_replace_all(Freq$Haplotype, "HLA-B", "HLA-B*")
Freq$Haplotype = str_replace_all(Freq$Haplotype, "HLA-C", "HLA-C*")
Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DR")] = paste0(substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DR")],1,11), ":",
                                                                      substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DR")],12,13))
Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DQ")] = paste0(substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DQ")],1,11), ":",
                                                                      substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DQ")],12,21),":",
                                                                      substr(Freq$Haplotype[which(substr(Freq$Haplotype,1,6) =="HLA-DQ")],22,23))
Freq$US[which(is.na(Freq$US))] = Freq$Cauc_Am[which(is.na(Freq$US))] #Move over DQ frequencies to same column


###Selection process
#Sliding window of 15, 21, 27 AA
#Within each window, look for total population frequency covered for each class I, class II, or both

seq =  fread(paste0(WORKING_DIR, "AA_sequence_combined.txt"), header=F) #SARS-CoV-2 proteome

##Filter by mouse, if necessary.  Not relevant for manuscript, but included for Ting lab peptide set####

#No filt
filt_I_m = filt_I
filt_II_m = filt_II

#Filter by haplotype -- only used for Ting lab peptide set, can ignore for manuscript
# #By B
# filt_I_m = filt_I[which((filt_I$Murine == "B") | (filt_I$Murine == "B/D")),]
# filt_II_m = filt_II[which((filt_II$Murine == "B") | (filt_II$Murine == "B/D")),]
# 
# #By D
# filt_I_m = filt_I[which((filt_I$Murine == "D") | (filt_I$Murine == "B/D")),]
# filt_II_m = filt_II[which((filt_II$Murine == "D") | (filt_II$Murine == "B/D")),]
# 
# #BxD
# filt_I_m - filt_I[which(!is.na(filt_I$Murine)),]
# filt_II_m - filt_I[which(!is.na(filt_II$Murine)),]


window = c(15,21,27)

for(w in window){
  print(w)
  Freq_by_window = data.table()
  count=0
  p=0
  
  ###Change this step for protein-filter vs no-filter
  for(n in 1:nchar(seq$V1)){ #For all proteins
  #for(n in c(7097:8369, 8720:8941, 9245:9663)){ #For S, M, N
    
    #Just a rough % counter to show progress
    count=count+1
    if(count >= 9701/100){ #Rough est for S, M, N length
    #if(count >= 1914/100){ #Rough est for S, M, N length
      count=0
      p = p+1
      print(paste0("Window ", w, ": ", p, "% done"))
    }
    
    if(n < 7097){
     Protein = "orf1ab"
     Start = n
     End = Start + w -1
    }else if((n >= 7097) & (n <= 8369)){  #Set the protein and positions wrt that protein
      Protein = "S"
      Start = n-7096
      End= Start+w-1
    }else if((n >= 8370) & (n <= 8644)){
      Protein = "ORF3a"
      Start = n-8719
      End= Start+w-1
    }else if((n >= 8645) & (n <= 8719)){
      Protein = "E"
      Start = n-8719
      End= Start+w-1
    }else if((n >= 8720) & (n <= 8941)){
      Protein = "M"
      Start = n-8719
      End= Start+w-1
    }else if((n >= 8942) & (n <= 9002)){
      Protein = "ORF6"
      Start = n-8719
      End= Start+w-1
    }else if((n >= 9003) & (n <= 9123)){
      Protein = "ORF7"
      Start = n-8719
      End= Start+w-1
    }else if((n >= 9124) & (n <= 9244)){
      Protein = "ORF8"
      Start = n-8719
      End= Start+w-1
    }else if ((n >= 9245) & (n <= 9663)){
      Protein = "N"
      Start = n-9244
      End= Start+w-1
    }else if((n >= 9664) & (n <= 9701)){
      Protein = "ORF10"
      Start = n-8719
      End= Start+w-1
    }else{
      Protein = NA
      Start = NA
      End = NA
    }
    
    range = seq(n,n+w-1,1)  #Set the window range
    sequence = substr(seq$V1,n,n+w-1) #Grab the sequence for that window range
    
    #Find all HLA-I/II epitopes that fall within the defined window
    c1_hits = filt_I_m[which((filt_I_m$Start %in% range)&(filt_I_m$End %in% range))] 
    c1_murine_b_hits = paste0(c1_hits$all_I.n.[which((c1_hits$Murine == "B") | (c1_hits$Murine == "B/D"))], collapse = ",")
    c1_murine_d_hits = paste0(c1_hits$all_I.n.[which((c1_hits$Murine == "D") | (c1_hits$Murine == "B/D"))], collapse = ",")
    
    c2_hits = filt_II_m[which((filt_II_m$Start %in% range)&(filt_II_m$End %in% range))]
    c2_murine_b_hits = paste0(c2_hits$all_II.n.[which((c2_hits$Murine == "B") | (c2_hits$Murine == "B/D"))], collapse = ",")
    c2_murine_d_hits = paste0(c2_hits$all_II.n.[which((c2_hits$Murine == "D") | (c2_hits$Murine == "B/D"))], collapse = ",")
    
    #Grab their haplotypes and pop. freq. (post-GLM filter)
    c1_haps = unique(unlist(strsplit(c1_hits$Immunogenic_haplotypes,",")))
    c1_freq=1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% c1_haps)]/100,2)))
    c2_haps = unique(unlist(strsplit(c2_hits$Immunogenic_haplotypes,",")))
    c2_freq=1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% c2_haps)]/100,2)))
    
    #Total pop. freq.
    co_freq = c1_freq*c2_freq
    
    Freq_by_window = rbind(Freq_by_window, data.table(sequence, Protein, Start, End, 
                                                      paste0(c1_hits$all_I.n., collapse = ","), paste(c1_haps, collapse = ","), c1_freq,
                                                      paste0(c2_hits$all_II.n., collapse = ","), paste(c2_haps, collapse = ","), c2_freq,
                                                      c1_murine_b_hits, c1_murine_d_hits, c2_murine_b_hits, c2_murine_d_hits))
    
  }
  
  colnames(Freq_by_window) = c("Sequence","Protein","Start", "End", "HLA-I_peptides", "HLA-I_haplotypes", "HLA-I_pop_freq", 
                               "HLA-II_peptides", "HLA-II_haplotypes", "HLA-II_pop_freq", "Mouse MHC-I b", "Mouse MHC-I d",
                               "Mouse MHC-II b", "Mouse MHC-II d")
  Freq_by_window$Total_frequency = Freq_by_window$`HLA-I_pop_freq` *  Freq_by_window$`HLA-II_pop_freq` 
  
  #Remove non-S/M/N entries; do not run if keeping all proteins (supplemental)
  #Freq_by_window = Freq_by_window[which(Freq_by_window$Protein %in% c("S", "M", "N")),]
  
  assign(paste0("Window_", w, "mer"), Freq_by_window)
}
#S, M, N
fwrite(Window_15mer, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Window_15mer.txt")
fwrite(Window_21mer, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Window_21mer.txt")
fwrite(Window_27mer, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Window_27mer.txt")

#All proteins
fwrite(Window_15mer, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Window_15mer_all_proteins.txt")
fwrite(Window_21mer, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Window_21mer_all_proteins.txt")
fwrite(Window_27mer, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Window_27mer_all_proteins.txt")

# ###Plotting population frequency for each window
# 
# p15mer = ggplot(data=Window_15mer)+
#   geom_point(aes(x=Start, y=Total_frequency), color = "purple")+
#   geom_point(aes(x=Start, y=`HLA-I_pop_freq`), color = "red")+
#   geom_point(aes(x=Start, y=`HLA-II_pop_freq`), color = "blue")+
#   scale_y_continuous(limits = c(0,1))+
#   scale_x_continuous(limits = c(7095,9717))+
#   labs(x = "", y = "Frequency", title = "15mer window")
# 
# p21mer = ggplot(data=Window_21mer)+
#   geom_point(aes(x=Start, y=Total_frequency), color = "purple")+
#   geom_point(aes(x=Start, y=`HLA-I_pop_freq`), color = "red")+
#   geom_point(aes(x=Start, y=`HLA-II_pop_freq`), color = "blue")+
#   scale_y_continuous(limits = c(0,1))+
#   scale_x_continuous(limits = c(7095,9717))+
#   labs(x = "", y = "Frequency", title = "21mer window")
# 
# p27mer = ggplot(data=Window_27mer)+
#   geom_point(aes(x=Start, y=Total_frequency), color = "purple")+
#   geom_point(aes(x=Start, y=`HLA-I_pop_freq`), color = "red")+
#   geom_point(aes(x=Start, y=`HLA-II_pop_freq`), color = "blue")+
#   scale_y_continuous(limits = c(0,1))+
#   scale_x_continuous(limits = c(7095,9717))+
#   labs(x = "", y = "Frequency", title = "27mer window")
# 
# grid.arrange(p15mer, p21mer, p27mer, nrow=3, ncol=1)
# 
# 
# ##Plotting number of epitopes per each window
# 
# Window_15mer$Num_c1 = sapply(strsplit(Window_15mer$`HLA-I_peptides`,","), length)
# Window_15mer$Num_c2 = sapply(strsplit(Window_15mer$`HLA-II_peptides`,","), length)
# 
# Window_21mer$Num_c1 = sapply(strsplit(Window_21mer$`HLA-I_peptides`,","), length)
# Window_21mer$Num_c2 = sapply(strsplit(Window_21mer$`HLA-II_peptides`,","), length)
# 
# Window_27mer$Num_c1 = sapply(strsplit(Window_27mer$`HLA-I_peptides`,","), length)
# Window_27mer$Num_c2 = sapply(strsplit(Window_27mer$`HLA-II_peptides`,","), length)
# 
# p15mer_count = ggplot(data=Window_15mer)+
#   geom_point(aes(x=Start, y=Num_c1), color = "red", alpha = .3)+
#   geom_point(aes(x=Start, y=Num_c2), color = "blue", alpha = .3)+
#   scale_y_continuous(limits = c(0,15))+  scale_x_continuous(limits = c(7096,9717))+
#   labs(x = "", y = "Count", title = "15mer window")
# 
# p21mer_count = ggplot(data=Window_21mer)+
#   geom_point(aes(x=Start, y=Num_c1), color = "red", alpha = .3)+
#   geom_point(aes(x=Start, y=Num_c2), color = "blue", alpha = .3)+
#   scale_y_continuous(limits = c(0,15))+  scale_x_continuous(limits = c(7096,9717))+
#   labs(x = "", y = "Count", title = "21mer window")
# 
# p27mer_count = ggplot(data=Window_27mer)+
#   geom_point(aes(x=Start, y=Num_c1), color = "red", alpha = .3)+
#   geom_point(aes(x=Start, y=Num_c2), color = "blue", alpha = .3)+
#   scale_y_continuous(limits = c(0,15))+  scale_x_continuous(limits = c(7096,9717))+
#   labs(x = "", y = "Count", title = "27mer window")
# 
# grid.arrange(p15mer_count, p21mer_count, p27mer_count, nrow=3, ncol=1)
# 
# 
# ######Selecting N peptides for Ting lab
# 
# #B haplotype (ran window caller on 27mers with B coverage)
# S_filt = Window_27mer[which((Window_27mer$Protein == "S") & (Window_27mer$Total_frequency > .25))] #B filter 25%
# S_plot=ggplot(S_filt)+
#   geom_segment(aes(x = Start, xend = End, y=0, yend=0))+
#   geom_point(aes(x= Start, y = Total_frequency)) +
#   labs(y="Population frequency", x="S protein coordinate")
# 
# M_filt = Window_27mer[which((Window_27mer$Protein == "M") & (Window_27mer$Total_frequency > .25))] #B filter 25%
# M_plot=ggplot(M_filt)+
#   geom_segment(aes(x = Start, xend = End, y=0, yend=0))+
#   geom_point(aes(x= Start, y = Total_frequency)) +
#   labs(y="Population frequency", x="M protein coordinate")
# 
# N_filt = Window_27mer[which((Window_27mer$Protein == "N") & (Window_27mer$Total_frequency > .06))] #B filter 6%
# N_plot=ggplot(N_filt)+
#   geom_segment(aes(x = Start, xend = End, y=0, yend=0))+
#   geom_point(aes(x= Start, y = Total_frequency)) +
#   labs(y="Population frequency", x="N protein coordinate")
# 
# grid.arrange(S_plot, M_plot, N_plot, nrow = 3, ncol = 1)
# 
# S_sub = S_filt[c(1, 10, 20, 34), ]
# M_sub = M_filt[c(1, 17, 25, 37), ]
# N_sub = N_filt[c(5, 13, 26, 48), ]
# All_sub = rbind(S_sub, M_sub, N_sub)
# 
# hla_I = unique(c(unlist(strsplit(S_sub$`HLA-I_haplotypes`,",")),
#                  unlist(strsplit(M_sub$`HLA-I_haplotypes`,",")),
#                  unlist(strsplit(N_sub$`HLA-I_haplotypes`,","))))
# hla_II = unique(c(unlist(strsplit(S_sub$`HLA-II_haplotypes`,",")),
#                  unlist(strsplit(M_sub$`HLA-II_haplotypes`,",")),
#                  unlist(strsplit(N_sub$`HLA-II_haplotypes`,","))))
# 
# hla_I_freq = 1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% hla_I)]/100,2)))
# hla_II_freq = 1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% hla_II)]/100,2)))
# 
# fwrite(All_sub, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/H2b_SARS-CoV-2_27mers.txt")


# ############Bc MHCII overlap########################
# 
# #Alex's B cell epitope data, post filtering
# bc = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/bcell/accessible-linear-bcell-epitopes-grouped-merged-filtered.csv")
# 
# Tc_match = data.table()
# All_tc_match = data.table()
# for (b in 1:nrow(bc)){
#   bc_start = bc$accessible_subsequence_start[b]
#   bc_end = bc$accessible_subsequence_end[b]
#   
#   ##Epitope 1 is longer than 15AA, so this "if" statement was added
#   if(bc$accessible_subsequence_length[b] <= 15){
#       total_15mers = Window_15mer[which((Window_15mer$Start <= bc_start)&(Window_15mer$End >= bc_end))] #Find Tc epitope windows that encompase the Bc epitope
#       max_15mer = total_15mers[which(total_15mers$`HLA-II_pop_freq` == max(total_15mers$`HLA-II_pop_freq`)),] #Grab the window with highest HLA-II pop freq
#       max_15mer$B_cell_epitope = bc$accessible_subsequence[b] #Get Bc epitope sequence
#       max_15mer$Peptide_length = 15 #Define length
#       
#   }else{
#     max_15mers = data.table()
#   }
#   
#   total_21mers = Window_21mer[which((Window_21mer$Start <= bc_start)&(Window_21mer$End >= bc_end))]
#   max_21mer = total_21mers[which(total_21mers$`HLA-II_pop_freq` == max(total_21mers$`HLA-II_pop_freq`)),]
#   max_21mer$B_cell_epitope = bc$accessible_subsequence[b]
#   max_21mer$Peptide_length = 21
#   
#   total_27mers = Window_27mer[which((Window_27mer$Start <= bc_start)&(Window_27mer$End >= bc_end))]
#   max_27mer = total_27mers[which(total_27mers$`HLA-II_pop_freq` == max(total_27mers$`HLA-II_pop_freq`)),]
#   max_27mer$B_cell_epitope = bc$accessible_subsequence[b]
#   max_27mer$Peptide_length = 27
#   
#   Tc_match = rbind(Tc_match, max_15mers, max_21mer, max_27mer)
#   All_tc_match = rbind(All_tc_match, total_15mers, total_21mers, total_27mers)
# }
# 
# #Reorder
# Tc_match = Tc_match[,c(12,13,1,2,3,4,5,6,7,8,9)]
# 
# fwrite(Tc_match, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Bc_epitopes_HLA_overlap.txt")
# 
# 

######Read back in Alex's rankings, make fig 5 tables#############

#Human only 27mer
CD4_h_27 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/selected-tcell-cd4-27mer.csv")
CD8_h_27 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/selected-tcell-cd8-27mer.csv")
CD4_CD8_h_27 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/selected-tcell-cd4-cd8-27mer.csv")
Bc_h_27 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/selected-bcell-27mer.csv")

#Check pop freq#
1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% unique(unlist(strsplit(CD4_h_27$`HLA-I_haplotypes`,','))))]/100,2)))
1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% unique(unlist(strsplit(CD4_h_27$`HLA-II_haplotypes`,','))))]/100,2)))
prod(1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% unique(unlist(strsplit(CD4_h_27$`HLA-I_haplotypes`,','))))]/100,2))), 
     1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% unique(unlist(strsplit(CD4_h_27$`HLA-II_haplotypes`,','))))]/100,2))))
#

total_h = rbind(cbind("CD4 optimized", CD4_h_27[,c(1,2,3,7,10)], NA),
                cbind("CD8 optimized", CD8_h_27[,c(1,2,3,7,10)], NA),
                cbind("CD4/CD8 optimized", CD4_CD8_h_27[,c(1,2,3,7,10)], NA),
                cbind("B cell optimized", Bc_h_27[,c(1,2,3,7,10,35)]), use.names=F)
total_h$`HLA-I_pop_freq` = paste0(round(total_h$`HLA-I_pop_freq`,3) * 100, "%")
total_h$`HLA-II_pop_freq` = paste0(round(total_h$`HLA-II_pop_freq`,3) * 100, "%")
colnames(total_h)[c(1,7)] = c("Set", "B cell epitope")
colnames(total_h)[c(5,6)] = c("HLA I frequency", "HLA II frequency")
#total_h$`HLA-I_peptides` = str_replace_all(total_h$`HLA-I_peptides`, ",", ", ")
#total_h$`HLA-II_peptides` = str_replace_all(total_h$`HLA-II_peptides`, ",", ", ")


ht <- hux(
  total_h,
  add_colnames = TRUE
)
bold(ht) = TRUE

ht = ht %>% set_background_color(row = c(2:4), col = c(1:7), alpha("blue", .5)) %>%
  set_background_color(row = c(5:8), col = c(1:7), alpha("red", .5)) %>%
  set_background_color(row = c(9:12), col = c(1:7), alpha("purple", .5)) %>%
  set_background_color(row = c(13:18), col = c(1:7), alpha("green", .5)) %>%
  set_background_color(row = 1, col = c(1:7), "grey95") %>% set_all_borders(row = c(1:18), col = c(1:7), 1) %>%
        set_outer_borders(0.4) #%>% map_background_color(by_rows("grey95", "white"))
#theme_plain(ht)

quick_pptx(ht, file = "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/human_27mers.pptx")


#Murine b/d 27er
CD4_bd_27 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/selected-tcell-cd4-h2b-h2d-27mer.csv")
CD8_bd_27 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/selected-tcell-cd8-h2b-h2d-27mer.csv")
CD4_CD8_bd_27 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/selected-tcell-cd4-cd8-h2b-h2d-27mer.csv")
Bc_bd_27 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/selected-bcell-27mer.csv")

#Check pop freq#
1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% unique(unlist(strsplit(CD4_bd_27$`HLA-I_haplotypes`,','))))]/100,2)))
1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% unique(unlist(strsplit(CD4_bd_27$`HLA-II_haplotypes`,','))))]/100,2)))
prod(1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% unique(unlist(strsplit(CD4_bd_27$`HLA-I_haplotypes`,','))))]/100,2))), 
     1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% unique(unlist(strsplit(CD4_bd_27$`HLA-II_haplotypes`,','))))]/100,2))))
#

total_bd = rbind(cbind("CD4 optimized", CD4_bd_27[,c(1,2,3,7,10)], NA),
                cbind("CD8 optimized", CD8_bd_27[,c(1,2,3,7,10)], NA),
                cbind("CD4/CD8 optimized", CD4_CD8_bd_27[,c(1,2,3,7,10)], NA),
                cbind("B cell optimized", Bc_bd_27[,c(1,2,3,7,10,35)]), use.names=F)
total_bd$`HLA-I_pop_freq` = paste0(round(total_bd$`HLA-I_pop_freq`,3) * 100, "%")
total_bd$`HLA-II_pop_freq` = paste0(round(total_bd$`HLA-II_pop_freq`,3) * 100, "%")
colnames(total_bd)[c(1,7)] = c("Set", "B cell epitope")
colnames(total_bd)[c(5,6)] = c("HLA I frequency", "HLA II frequency")
#total_h$`HLA-I_peptides` = str_replace_all(total_h$`HLA-I_peptides`, ",", ", ")
#total_h$`HLA-II_peptides` = str_replace_all(total_h$`HLA-II_peptides`, ",", ", ")


ht_m <- hux(
  total_bd,
  add_colnames = TRUE
)
bold(ht_m) = TRUE

ht_m = ht_m %>% set_background_color(row = c(2:3), col = c(1:7), alpha("blue", .5)) %>%
  set_background_color(row = c(4:7), col = c(1:7), alpha("red", .5)) %>%
  set_background_color(row = c(8:10), col = c(1:7), alpha("purple", .5)) %>%
  set_background_color(row = c(11:16), col = c(1:7), alpha("green", .5)) %>%
  set_background_color(row = 1, col = c(1:7), "grey95") %>% set_all_borders(row = c(1:16), col = c(1:7), 1) %>%
  set_outer_borders(0.4) #%>% map_background_color(by_rows("grey95", "white"))
#theme_plain(ht)

quick_pptx(ht_m, file = "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/murine_27mers.pptx")


###Make table of total pop freq coverage by set
Freq = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/HLA_freq_formatted.csv")


#Human only 27mer
CD4_h_27 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/selected-tcell-cd4-27mer.csv")
CD8_h_27 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/selected-tcell-cd8-27mer.csv")
CD4_CD8_h_27 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/selected-tcell-cd4-cd8-27mer.csv")

#Murine b 27er
CD4_b_27 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/selected-tcell-cd4-h2b-27mer.csv")
CD8_b_27 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/selected-tcell-cd8-h2b-27mer.csv")
CD4_CD8_b_27 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/selected-tcell-cd4-cd8-h2b-27mer.csv")

#Murine d 27er
CD4_d_27 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/selected-tcell-cd4-h2d-27mer.csv")
CD8_d_27 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/selected-tcell-cd8-h2d-27mer.csv")
CD4_CD8_d_27 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/selected-tcell-cd4-cd8-h2d-27mer.csv")

#Murine b/d 27er
CD4_bd_27 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/selected-tcell-cd4-h2b-h2d-27mer.csv")
CD8_bd_27 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/selected-tcell-cd8-h2b-h2d-27mer.csv")
CD4_CD8_bd_27 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/selected-tcell-cd4-cd8-h2b-h2d-27mer.csv")

Bc_27 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/selected-bcell-27mer.csv")


sets = list(CD4_h_27, CD8_h_27, CD4_CD8_h_27,
         CD4_b_27, CD8_b_27, CD4_CD8_b_27,
         CD4_d_27, CD8_d_27, CD4_CD8_d_27,
         CD4_bd_27, CD8_bd_27, CD4_CD8_bd_27,
         Bc_27)
sets_names = c("CD4_h_27", " CD8_h_27", " CD4_CD8_h_27",
            "CD4_b_27", " CD8_b_27", " CD4_CD8_b_27",
            "CD4_d_27", " CD8_d_27", " CD4_CD8_d_27",
            "CD4_bd_27", " CD8_bd_27", " CD4_CD8_bd_27", "Bc_27")
Freq_table = data.table()

for (s in 1:length(sets)){
  print(s)
  s_mat = as.data.table(sets[s], check.names = F)
  Freq_table = rbind(Freq_table,
                     data.table(  sets_names[s],
                                  1-prod(1-(rep(Freq$Frequency[which(Freq$Haplotype %in% unique(unlist(strsplit(s_mat$`HLA-I_haplotypes`,','))))]/100,2))),
                                  1-prod(1-(rep(Freq$Frequency[which(Freq$Haplotype %in% unique(unlist(strsplit(s_mat$`HLA-II_haplotypes`,','))))]/100,2))),
                                  prod(1-prod(1-(rep(Freq$Frequency[which(Freq$Haplotype %in% unique(unlist(strsplit(s_mat$`HLA-I_haplotypes`,','))))]/100,2))), 
                                     1-prod(1-(rep(Freq$Frequency[which(Freq$Haplotype %in% unique(unlist(strsplit(s_mat$`HLA-II_haplotypes`,','))))]/100,2)))),
                                  
                                  1-prod(1-(rep(Freq$EUR_freq[which(Freq$Haplotype %in% unique(unlist(strsplit(s_mat$`HLA-I_haplotypes`,','))))]/100,2))),
                                  1-prod(1-(rep(na.omit(Freq$EUR_freq[which(Freq$Haplotype %in% unique(unlist(strsplit(s_mat$`HLA-II_haplotypes`,','))))])/100,2))),
                                  #prod(1-prod(1-(rep(Freq$EUR_freq[which(Freq$Haplotype %in% unique(unlist(strsplit(s_mat$`HLA-I_haplotypes`,','))))]/100,2))), 
                                  #     1-prod(1-(rep(Freq$EUR_freq[which(Freq$Haplotype %in% unique(unlist(strsplit(s_mat$`HLA-II_haplotypes`,','))))]/100,2)))),
                                  
                                  1-prod(1-(rep(Freq$AFA_freq[which(Freq$Haplotype %in% unique(unlist(strsplit(s_mat$`HLA-I_haplotypes`,','))))]/100,2))),
                                  1-prod(1-(rep(na.omit(Freq$AFA_freq[which(Freq$Haplotype %in% unique(unlist(strsplit(s_mat$`HLA-II_haplotypes`,','))))])/100,2))),
                                  #prod(1-prod(1-(rep(Freq$AFA_freq[which(Freq$Haplotype %in% unique(unlist(strsplit(s_mat$`HLA-I_haplotypes`,','))))]/100,2))), 
                                  #     1-prod(1-(rep(Freq$AFA_freq[which(Freq$Haplotype %in% unique(unlist(strsplit(s_mat$`HLA-II_haplotypes`,','))))]/100,2)))),
                                  
                                  1-prod(1-(rep(Freq$API_freq[which(Freq$Haplotype %in% unique(unlist(strsplit(s_mat$`HLA-I_haplotypes`,','))))]/100,2))),
                                  1-prod(1-(rep(na.omit(Freq$API_freq[which(Freq$Haplotype %in% unique(unlist(strsplit(s_mat$`HLA-II_haplotypes`,','))))])/100,2))),
                                  #prod(1-prod(1-(rep(Freq$API_freq[which(Freq$Haplotype %in% unique(unlist(strsplit(s_mat$`HLA-I_haplotypes`,','))))]/100,2))), 
                                  #     1-prod(1-(rep(Freq$API_freq[which(Freq$Haplotype %in% unique(unlist(strsplit(s_mat$`HLA-II_haplotypes`,','))))]/100,2)))),
                                  
                                  1-prod(1-(rep(Freq$HIS_freq[which(Freq$Haplotype %in% unique(unlist(strsplit(s_mat$`HLA-I_haplotypes`,','))))]/100,2))),
                                  1-prod(1-(rep(na.omit(Freq$HIS_freq[which(Freq$Haplotype %in% unique(unlist(strsplit(s_mat$`HLA-II_haplotypes`,','))))])/100,2)))))
                                  #prod(1-prod(1-(rep(Freq$HIS_freq[which(Freq$Haplotype %in% unique(unlist(strsplit(s_mat$`HLA-I_haplotypes`,','))))]/100,2))), 
                                  #     1-prod(1-(rep(Freq$HIS_freq[which(Freq$Haplotype %in% unique(unlist(strsplit(s_mat$`HLA-II_haplotypes`,','))))]/100,2))))))
}
colnames(Freq_table) = c("Set", "HLA-I frequency", "HLA-II frequency", "Overall frequency", 
                         "EUR HLA-I frequency", "EUR DRB1 frequency",
                         "AFA HLA-I frequency", "AFA DRB1 frequency",
                         "API HLA-I frequency", "API DRB1 frequency",
                         "HIS HLA-I frequency", "HIS DRB1 frequency")
fwrite(Freq_table, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/Frequency_table_27mer_all_ethnic_groups.csv")


########Supplmenental tables#####################
#################################################
#Table S1, Table S2

#Read in frequency table
Freq = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/HLA_freq_formatted.csv")

#Read in predicted epitopes
S1 =fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAI_all_glm.txt")
S2 =fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAII_all_glm.txt")

#Read in entropy
ent = fread(paste0(WORKING_DIR, "entropy_7882.txt"))

#Read in SARS data
iedb_fig3 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/mhc_sars2_purified_competitive_radioactivity.csv")
iedb_fig3 = iedb_fig3[which(iedb_fig3$affinity < 500),]
iedb_fig3 = iedb_fig3[which(iedb_fig3$allele %in% Freq$Haplotype),]

#Un-normalize
S1$BA_Rank = S1$BA_Rank*100
S1$EL_Rank = S1$EL_Rank*100
S1$Flurry_Rank = S1$Flurry_Rank*100
S1$Binding_affinity = S1$Binding_affinity*50000
S1$Flurry_BA = S1$Flurry_BA*50000

S2$Binding_affinity = S2$Binding_affinity * 50000
S2$BA_Rank = S2$BA_Rank * 100
S2$EL_Rank = S2$EL_Rank * 100

#Add binder, entropy
S1$Predicted_binder = ifelse(S1$Binding_affinity < 393.4, "True", "False")
S2$Predicted_binder = ifelse(S2$v3.2_Binding_affinity < 220, "True", "False")

#orf1ab: 1 - 7096                                                                    
#S: 7097 - 8369                                                                       
#ORF3a: 8370 - 8644                                                                   
#E: 8645 - 8719                                                                       
#M: 8720 - 8941                                                                       
#ORF6: 8942 - 9002                                                                    
#ORF7a: 9003 - 9123                                                                   
#ORF8: 9124 - 9244                                                                    
#N: 9245 - 9663                                                                       
#ORF10: 9664 - 9701

#Convert to proteome space
S1$Start[which(S1$Protein == "S")] = S1$Start[which(S1$Protein == "S")]+7096
S1$Start[which(S1$Protein == "ORF3a")] = S1$Start[which(S1$Protein == "ORF3a")]+8369
S1$Start[which(S1$Protein == "E")] = S1$Start[which(S1$Protein == "E")]+8644
S1$Start[which(S1$Protein == "M")] = S1$Start[which(S1$Protein == "M")]+8719
S1$Start[which(S1$Protein == "ORF6")] = S1$Start[which(S1$Protein == "ORF6")]+8941
S1$Start[which(S1$Protein == "ORF7a")] = S1$Start[which(S1$Protein == "ORF7a")]+9002
S1$Start[which(S1$Protein == "ORF8")] = S1$Start[which(S1$Protein == "ORF8")]+9123
S1$Start[which(S1$Protein == "N")] = S1$Start[which(S1$Protein == "N")]+9244
S1$Start[which(S1$Protein == "ORF10")] = S1$Start[which(S1$Protein == "ORF10")]+9663

#Add polymorphism data
ent_high = ent$V1[which(ent$V2 > 0.1)]
S1$Low_polymorphism = "True"
S1$End = S1$Start + nchar(S1$Peptide) - 1

for(e in ent_high){
  print(e)
  S1$Low_polymorphism[which((e >= S1$Start) & (e <= S1$End))] = "False"
}
S1$End = NULL

#Convert to protein space
S1$Start[which(S1$Protein == "S")] = S1$Start[which(S1$Protein == "S")]-7096
S1$Start[which(S1$Protein == "ORF3a")] = S1$Start[which(S1$Protein == "ORF3a")]-8369
S1$Start[which(S1$Protein == "E")] = S1$Start[which(S1$Protein == "E")]-8644
S1$Start[which(S1$Protein == "M")] = S1$Start[which(S1$Protein == "M")]-8719
S1$Start[which(S1$Protein == "ORF6")] = S1$Start[which(S1$Protein == "ORF6")]-8941
S1$Start[which(S1$Protein == "ORF7a")] = S1$Start[which(S1$Protein == "ORF7a")]-9002
S1$Start[which(S1$Protein == "ORF8")] = S1$Start[which(S1$Protein == "ORF8")]-9123
S1$Start[which(S1$Protein == "N")] = S1$Start[which(S1$Protein == "N")]-9244
S1$Start[which(S1$Protein == "ORF10")] = S1$Start[which(S1$Protein == "ORF10")]-9663

#Add if peptide is in IEDB set
S1$IEDB_binder = "False"
S1$IEDB_binder[which(paste0(S1$Haplotype,"_",S1$Peptide) %in% paste0(iedb_fig3$allele, "_", iedb_fig3$peptide))] = "True"


#Check for consistency
S1_filt = S1[which((S1$Predicted_binder=="True") | (S1$IEDB_binder == "True") ),]
length(unique(S1_filt$Peptide))

S1_filt = S1_filt[which(S1_filt$Predicted_tetramer == TRUE),]
length(unique(S1_filt$Peptide))

S1_filt = S1_filt[which(S1_filt$Low_polymorphism == "True"),]
length(unique(S1_filt$Peptide))

S1_filt = S1_filt[which(S1_filt$Protein %in% c("S", "M", "N")),]
length(unique(S1_filt$Peptide))

colnames(S1)[4] = "Allele"
colnames(S1)[5] = "NetMHCpan_4.0_BA_Rank"
colnames(S1)[6] = "NetMHCpan_4.0_BA_Score"
colnames(S1)[7] = "NetMHCpan_4.0_Binding_affinity"
colnames(S1)[8] = "NetMHCpan_4.0_EL_Rank"
colnames(S1)[9] = "NetMHCpan_4.0_EL_Score"
colnames(S1)[10] = "MHCflurry_Rank"
colnames(S1)[11] = "MHCflurry_Binding_affinity"
colnames(S1)[12] = "MHCflurry_Processing_score"
colnames(S1)[13] = "MHCflurry_Presentation_score"
colnames(S1)[14:19] = paste0("%_", colnames(S1)[14:19])


fwrite(S1, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Table_S1_HLA-I_epitopes.csv")  


#Convert to proteome space
S2$Start[which(S2$Protein == "S")] = S2$Start[which(S2$Protein == "S")]+7096
S2$Start[which(S2$Protein == "ORF3a")] = S2$Start[which(S2$Protein == "ORF3a")]+8369
S2$Start[which(S2$Protein == "E")] = S2$Start[which(S2$Protein == "E")]+8644
S2$Start[which(S2$Protein == "M")] = S2$Start[which(S2$Protein == "M")]+8719
S2$Start[which(S2$Protein == "ORF6")] = S2$Start[which(S2$Protein == "ORF6")]+8941
S2$Start[which(S2$Protein == "ORF7a")] = S2$Start[which(S2$Protein == "ORF7a")]+9002
S2$Start[which(S2$Protein == "ORF8")] = S2$Start[which(S2$Protein == "ORF8")]+9123
S2$Start[which(S2$Protein == "N")] = S2$Start[which(S2$Protein == "N")]+9244
S2$Start[which(S2$Protein == "ORF10")] = S2$Start[which(S2$Protein == "ORF10")]+9663

#Add polymorphism data
ent_high = ent$V1[which(ent$V2 > 0.1)]
S2$Low_polymorphism = "True"
S2$End = S2$Start + nchar(S2$Peptide) - 1

for(e in ent_high){
  print(e)
  S2$Low_polymorphism[which((e >= S2$Start) & (e <= S2$End))] = "False"
}
S2$End = NULL

#Convert to protein space
S2$Start[which(S2$Protein == "S")] = S2$Start[which(S2$Protein == "S")]-7096
S2$Start[which(S2$Protein == "ORF3a")] = S2$Start[which(S2$Protein == "ORF3a")]-8369
S2$Start[which(S2$Protein == "E")] = S2$Start[which(S2$Protein == "E")]-8644
S2$Start[which(S2$Protein == "M")] = S2$Start[which(S2$Protein == "M")]-8719
S2$Start[which(S2$Protein == "ORF6")] = S2$Start[which(S2$Protein == "ORF6")]-8941
S2$Start[which(S2$Protein == "ORF7a")] = S2$Start[which(S2$Protein == "ORF7a")]-9002
S2$Start[which(S2$Protein == "ORF8")] = S2$Start[which(S2$Protein == "ORF8")]-9123
S2$Start[which(S2$Protein == "N")] = S2$Start[which(S2$Protein == "N")]-9244
S2$Start[which(S2$Protein == "ORF10")] = S2$Start[which(S2$Protein == "ORF10")]-9663


#Add if peptide is in IEDB set
S2$IEDB_binder = "False"
S2$IEDB_binder[which(paste0(S2$Haplotype,"_",S2$Peptide) %in% paste0(iedb_fig3$allele, "_", iedb_fig3$peptide))] = "True"


#Check for consistency
S2_filt = S2[which((S2$Predicted_binder=="True") | (S2$IEDB_binder == "True") ),]
length(unique(S2_filt$Peptide))

S2_filt = S2_filt[which(S2_filt$Predicted_tetramer == T),]
length(unique(S2_filt$Peptide))

S2_filt = S2_filt[which(S2_filt$Low_polymorphism == "True"),]
length(unique(S2_filt$Peptide))

S2_filt = S2_filt[which(S2_filt$Protein %in% c("S", "M", "N")),]
length(unique(S2_filt$Peptide))

colnames(S2)[4] = "Allele"
colnames(S2)[5] = "NetMHCIIpan_3.2_Rank"
colnames(S2)[6] = "NetMHCIIpan_3.2_BA_score"
colnames(S2)[7] = "NetMHCIIpan_3.2_Binding_affinity"
colnames(S2)[8] = "NetMHCIIpan_4.0_Binding_affinity"
colnames(S2)[9] = "NetMHCIIpan_4.0_BA_Rank"
colnames(S2)[10] = "NetMHCIIpan_4.0_BA_score"
colnames(S2)[11] = "NetMHCIIpan_4.0_EL_Rank"
colnames(S2)[12] = "NetMHCIIpan_4.0_EL_score"
colnames(S2)[c(13:18)] = paste0("%_", colnames(S2)[c(13:18)])

fwrite(S2, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Table_S2_HLA-II_epitopes.csv")  
  

#Table S3
#"/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Figures/COVID/Supplemental/Standardized T cell epitopes.txt"

#Table S4
Bc = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/bcell/linear-bcell-epitopes-SARS2-S-with-filters.csv")
Bc$Source[which(Bc$Source == "ReScan")] = "Zamecnik 2020"
#Bc = Bc[, c(1, 3, 4, 13, 7)]
fwrite(Bc, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Table_S4_B_cell_epitopes.csv")  



###Table S5, Add NMDP stats
Freq_form = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/HLA_freq_formatted.csv")


final_15mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/final-vaccine-peptides-15mer.csv")
final_21mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/final-vaccine-peptides-21mer.csv")
final_27mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/final-vaccine-peptides-27mer.csv")

final_all = rbind(final_15mer, final_21mer, final_27mer)

EUR_HLA_I_coverage = c()
EUR_DRB_coverage = c()
AFA_HLA_I_coverage = c()
AFA_DRB_coverage = c()
API_HLA_I_coverage = c()
API_DRB_coverage = c()
HIS_HLA_I_coverage = c()
HIS_DRB_coverage = c()

for(n in 1:nrow(final_all)){
  print(n)
  EUR_HLA_I_coverage = c(EUR_HLA_I_coverage, 1-prod(1-(rep(Freq_form$EUR_freq[which(substr(Freq_form$Haplotype,nchar(Freq_form$Haplotype)-6,nchar(Freq_form$Haplotype)) %in% 
                                                                                      unlist(strsplit(final_all$`HLA-I alleles`[n],' ')))]/100,2))))
  EUR_DRB_coverage = c(EUR_DRB_coverage, 1-prod(1-(rep(Freq_form$EUR_freq[which(substr(Freq_form$Haplotype,nchar(Freq_form$Haplotype)-9,nchar(Freq_form$Haplotype)) %in% 
                                                                                  unlist(strsplit(final_all$`HLA-II alleles`[n],' ')))]/100,2))))
  
  AFA_HLA_I_coverage = c(AFA_HLA_I_coverage, 1-prod(1-(rep(Freq_form$AFA_freq[which(substr(Freq_form$Haplotype,nchar(Freq_form$Haplotype)-6,nchar(Freq_form$Haplotype)) %in% 
                                                                                      unlist(strsplit(final_all$`HLA-I alleles`[n],' ')))]/100,2))))
  AFA_DRB_coverage = c(AFA_DRB_coverage, 1-prod(1-(rep(Freq_form$AFA_freq[which(substr(Freq_form$Haplotype,nchar(Freq_form$Haplotype)-9,nchar(Freq_form$Haplotype)) %in% 
                                                                                  unlist(strsplit(final_all$`HLA-II alleles`[n],' ')))]/100,2))))
  
  API_HLA_I_coverage = c(API_HLA_I_coverage, 1-prod(1-(rep(Freq_form$API_freq[which(substr(Freq_form$Haplotype,nchar(Freq_form$Haplotype)-6,nchar(Freq_form$Haplotype)) %in% 
                                                                                      unlist(strsplit(final_all$`HLA-I alleles`[n],' ')))]/100,2))))
  API_DRB_coverage = c(API_DRB_coverage, 1-prod(1-(rep(Freq_form$API_freq[which(substr(Freq_form$Haplotype,nchar(Freq_form$Haplotype)-9,nchar(Freq_form$Haplotype)) %in% 
                                                                                  unlist(strsplit(final_all$`HLA-II alleles`[n],' ')))]/100,2))))
  
  HIS_HLA_I_coverage = c(HIS_HLA_I_coverage, 1-prod(1-(rep(Freq_form$HIS_freq[which(substr(Freq_form$Haplotype,nchar(Freq_form$Haplotype)-6,nchar(Freq_form$Haplotype)) %in% 
                                                                                      unlist(strsplit(final_all$`HLA-I alleles`[n],' ')))]/100,2))))
  HIS_DRB_coverage = c(HIS_DRB_coverage, 1-prod(1-(rep(Freq_form$HIS_freq[which(substr(Freq_form$Haplotype,nchar(Freq_form$Haplotype)-9,nchar(Freq_form$Haplotype)) %in% 
                                                                                  unlist(strsplit(final_all$`HLA-II alleles`[n],' ')))]/100,2))))
  
}

final_all = cbind(final_all[,1], nchar(final_all$Sequence), final_all[,c(2:ncol(final_all)), with=F], data.table(EUR_HLA_I_coverage, EUR_DRB_coverage,
                                                                                                                 AFA_HLA_I_coverage, AFA_DRB_coverage,
                                                                                                                 API_HLA_I_coverage, API_DRB_coverage,
                                                                                                                 HIS_HLA_I_coverage, HIS_DRB_coverage))  
colnames(final_all)[2] = "Length"
fwrite(final_all, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/vaccine-peptides/Table_S5_final_vaccine-peptides-merged.csv")




##############ELISpot analysis -- not for current manuscript, can ignore for now####################
####################################################################################################

el1 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/ELISpot/051120_COVID_ELISpot.txt")
el1$Mean_count = NA
el1$Mean_intensity= NA

for (n in 1:nrow(el1)){
  el1$Mean_count[n] = sum(el1$Count_1[n], el1$Count_2[n], el1$Count_3[n])/3
  el1$Mean_intensity[n] = sum(el1$Intensity_1[n], el1$Intensity_2[n], el1$Intensity_3[n])/3
}

#Remove m11 -- failed positive control
el1 = el1[which(el1$Sample != "M11"),]

el1_neg = el1[which(el1$Peptide == "Negative"),]
el1_pos = el1[which(el1$Peptide == "Positive"),]
el1 = el1[which(el1$Peptide %ni% c("Positive", "Negative")),]

for(n in 1:nrow(el1)){
  el1$Mean_count[n] = el1$Mean_count[n] - el1_neg$Mean_count[which(el1_neg$Sample == el1$Sample[n])]
  el1$Mean_intensity[n] = el1$Mean_intensity[n] - el1_neg$Mean_intensity[which(el1_neg$Sample == el1$Sample[n])]
}


el1$Group = ifelse(substr(el1$Sample,1,1) == "C", "PBS", "Peptide")
el1$Adjuvant = ifelse((substr(el1$Sample,1,1) == "C"), "None", ifelse(substr(el1$Sample,2, nchar(el1$Sample)) %in% c(1,2,3,4,5,6), "MQ", "PM")) 
el1$Adjuvant = factor(el1$Adjuvant, levels = c("None", "MQ", "PM"))

#wilcox_counts = c()
#wilcox_intensity = c()
t_counts_1_2 = c()
t_counts_1_3 = c()
t_counts_2_3 = c()
t_int_1_2 = c()
t_int_1_3 = c()
t_int_2_3 = c()

for(p in unique(el1$Peptide)){
  None_count = as.numeric(el1$Mean_count[which((el1$Peptide == p)&(el1$Adjuvant == "None"))])
  MQ_count = as.numeric(el1$Mean_count[which((el1$Peptide == p)&(el1$Adjuvant == "MQ"))])
  PM_count = as.numeric(el1$Mean_count[which((el1$Peptide == p)&(el1$Adjuvant == "PM"))])
  
  #wilcox_counts = c(wilcox_counts, wilcox.test(Pep_count, Control_count)$p.value)
  t_counts_1_2 = c(t_counts_1_2, t.test(None_count, MQ_count, var.equal = F)$p.value)
  t_counts_1_3 = c(t_counts_1_3, t.test(None_count, PM_count, var.equal = F)$p.value)
  t_counts_2_3 = c(t_counts_2_3, t.test(MQ_count, PM_count, var.equal = F)$p.value)
  
  None_int = as.numeric(el1$Mean_intensity[which((el1$Peptide == p)&(el1$Adjuvant == "None"))])
  MQ_int = as.numeric(el1$Mean_intensity[which((el1$Peptide == p)&(el1$Adjuvant == "MQ"))])
  PM_int = as.numeric(el1$Mean_intensity[which((el1$Peptide == p)&(el1$Adjuvant == "PM"))])
  
  #wilcox_ints = c(wilcox_ints, wilcox.test(Pep_int, Control_int)$p.value)
  t_int_1_2 = c(t_int_1_2, t.test(None_int, MQ_int, var.equal = F)$p.value)
  t_int_1_3 = c(t_int_1_3, t.test(None_int, PM_int, var.equal = F)$p.value)
  t_int_2_3 = c(t_int_2_3, t.test(MQ_int, PM_int, var.equal = F)$p.value)
  
}

el_count = ggplot(data = el1, aes(x= Peptide, y = Mean_count, color = Adjuvant, shape = Adjuvant)) +
  stat_boxplot(outlier.shape = NA, alpha = 0.7)+
  geom_point(position=position_jitterdodge(), size = 4) +
  annotate("text", x = c(1,2,3,4,5,"M"), y = rep(23,6),
           label = paste0("None vs MQ: ", round(t_counts_1_2, 3),
                          "\nNone vs PM: ", round(t_counts_1_3, 3),
                          "\nMQ vs PM: ", round(t_counts_2_3, 3)))+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x = "", y = "Mean number of spots (n = 3)") 

el_int = ggplot(data = el1, aes(x= Peptide, y = Mean_intensity, color = Adjuvant, shape = Adjuvant)) +
  stat_boxplot(outlier.shape = NA, alpha = 0.7)+
  geom_point(position=position_jitterdodge(), size = 4) +
  annotate("text", x = c(1,2,3,4,5,"M"), y = rep(760,6),
           label = paste0("None vs MQ: ", round(t_int_1_2, 3),
                          "\nNone vs PM: ", round(t_int_1_3, 3),
                          "\nMQ vs PM: ", round(t_int_2_3, 3)))+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x = "", y = "Mean well intensity (n = 3)")

grid.arrange(el_count, el_int, nrow = 2, ncol = 1)



####082420 ELISpot

tab = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/ELISpot/082420_COVID_ELISpot_melted.txt")
tab = tab[complete.cases(tab$Sample_Pep),]

tab$Sample = sapply(strsplit(tab$Sample_Pep, "-"), "[", 1)
tab$Peptide = sapply(strsplit(tab$Sample_Pep, "-"), "[", 2)
tab = tab[-which(tab$Group == "BD*"),]

tab$Count_1 = as.numeric(tab$Count_1)
tab$Count_2 = as.numeric(tab$Count_2)
tab$Intensity_1 = as.numeric(tab$Intensity_1)
tab$Intensity_2 = as.numeric(tab$Intensity_2)

PHA_tab = tab[which(tab$Group == "PHA"),]
Neg_tab = tab[which(tab$Group == "NEG"),]

Sample_tab = tab[-which((tab$Group == "PHA") | (tab$Group == "NEG")),]

Sample_tab$Good_PHA = NA
for(n in 1:nrow(Sample_tab)){
  print(n)

  Sample_tab$Count_1[n] = Sample_tab$Count_1[n] - Neg_tab$Count_1[which(Neg_tab$Sample == Sample_tab$Sample[n])][1]
  Sample_tab$Count_2[n] = Sample_tab$Count_2[n] - Neg_tab$Count_2[which(Neg_tab$Sample == Sample_tab$Sample[n])][1]
  Sample_tab$Intensity_1[n] = Sample_tab$Intensity_1[n] - Neg_tab$Intensity_1[which(Neg_tab$Sample == Sample_tab$Sample[n])][1]
  Sample_tab$Intensity_2[n] = Sample_tab$Intensity_2[n] - Neg_tab$Intensity_2[which(Neg_tab$Sample == Sample_tab$Sample[n])][1]
  
  Sample_tab$Good_PHA[n] = ifelse(mean(PHA_tab$Intensity_1[which(PHA_tab$Sample == Sample_tab$Sample[n])], PHA_tab$Intensity_2[which(PHA_tab$Sample == Sample_tab$Sample[n])]) > 1000,
         T, F)
}

Sample_tab$Mean_Count = rowMeans(Sample_tab[,c(3,4)])
Sample_tab$Mean_Intensity = rowMeans(Sample_tab[,c(5,6)])

Sample_tab_filt = Sample_tab[which(Sample_tab$Good_PHA == T),]

groups = c("A", "B", "AB", "AD", "BD", "ABD")

for(g in groups){
  temp_filt = Sample_tab_filt[which((Sample_tab_filt$Group == g) | (Sample_tab_filt$Group == "ADJ")), ]# | (Sample_tab_filt$Group == "PBS")),]
  temp_filt = temp_filt[which(temp_filt$Peptide %in% unique(temp_filt$Peptide[which(temp_filt$Group == g)]))]
  
  temp_filt$Group = factor(temp_filt$Group, levels = c(g, "ADJ"))
  
  wilcox_g_ADJ = c()
  #wilcox_g_PBS = c()
  for(p in unique(temp_filt$Peptide)[order(unique(temp_filt$Peptide))]){
    #print(p)
    wilcox_g_ADJ = c(wilcox_g_ADJ, wilcox.test(temp_filt$Mean_Intensity[which((temp_filt$Group == g) & (temp_filt$Peptide == p))],
                                               temp_filt$Mean_Intensity[which((temp_filt$Group == "ADJ") & (temp_filt$Peptide == p))])$p.value)
    #wilcox_g_PBS = c(wilcox_g_PBS, wilcox.test(temp_filt$Mean_Intensity[which((temp_filt$Group == g) & (temp_filt$Peptide == p))],
    #                                           temp_filt$Mean_Intensity[which((temp_filt$Group == "PBS") & (temp_filt$Peptide == p))])$p.value)
  }
  
  assign(paste0(g, "_plot"), 
    ggplot(data = temp_filt, aes(x= Peptide, y = Mean_Intensity, color = Group)) +
    scale_color_discrete(guide = F) +
    stat_boxplot(outlier.shape = NA, alpha = 0.7)+
    geom_point(position=position_jitterdodge(), size = 4) +
    annotate("text", x = unique(temp_filt$Peptide)[order(unique(temp_filt$Peptide))], 
             y = rep(.75*max(temp_filt$Mean_Intensity),length(unique(temp_filt$Peptide))),
             angle = 45,
             label = paste0("p = ", round(wilcox_g_ADJ, 3))) +
                           # "\nVaccine vs PBS: ", round(wilcox_g_PBS, 3))+
    theme(text=element_text(face="bold",size=15,colour="black")) +
    labs(x = "", y = "Mean well intensity (n = 3)", title = paste0("Group: ", g)))
  
}

grid.arrange(A_plot, B_plot, AB_plot, AD_plot, BD_plot, ABD_plot, nrow = 3, ncol = 2)

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

############################################################
#############Dependencies and working directory############
##########################################################
'%ni%' <- Negate('%in%')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install()

packages <- c("scales", "data.table", "ggrepel", "ggplot2", "viridis", "ggnewscale", "seqinr", "DESeq2", "GenomicRanges", "gplots", "ggbeeswarm", "ggallin", "stringr")

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

WORKING_DIR = "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Figures/COVID/"


###########################################
######Pre-processing#######################

####Create fasta from IEDB, all viruses###############

##Create fa to run tools
Alleles_to_keep=c("HLA-DRB1*01:01", "HLA-DRB1*03:01", "HLA-DRB1*04:01", "HLA-DRB1*07:01", "HLA-DRB1*11:01", "HLA-DRB1*13:01", "HLA-DRB1*15:01",
                  "HLA-DQA1*01:02/DQB1*06:02", "HLA-DQA1*05:01/DQB1*02:01", "HLA-DQA1*02:01/DQB1*02:02", "HLA-DQA1*05:05/DQB1*03:01",
                  "HLA-DQA1*01:01/DQB1*05:01", "HLA-DQA1*03:01/DQB1*03:02", "HLA-DQA1*03:03/DQB1*03:01", "HLA-DQA1*01:03/DQB1*06:03")


iedb_v = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/mhc_all_virus_combined.csv")
iedb_v = iedb_v[which(iedb_v$allele%in% Alleles_to_keep),]
iedb_v = iedb_v[which(nchar(iedb_v$peptide) == 15),]

peps = unique(iedb_v$peptide)
id = paste0(">",seq(1,length(peps),1))
fa =c()
for(n in 1:length(id)){
  fa = c(fa, id[n])
  fa = c(fa, peps[n])
}
fa = as.data.table(fa)
fa1 = fa[1:4000,]
fa2 = fa[4001:8000,]
fa3 = fa[8001:12662,]
fwrite(fa3, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_all_viruses_combined_15mer_3.fa", col.names = F, row.names = F)


#########NetMHCpan BA mode########################
I_8mer=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_8mer.xls")
I_9mer=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_9mer.xls")
I_10mer=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_10mer.xls")
I_11mer=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_11mers.xls")
I_B1501=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_B1501_8-11mers.txt")

HLAI=rbind(I_8mer, I_9mer, I_10mer, I_11mer)
HLAI=HLAI[order(HLAI$Peptide),]
I_B1501 = I_B1501[order(I_B1501$Peptide),]
HLAI_BA=cbind(HLAI[,1:38], I_B1501[,4:8], HLAI[,39:ncol(HLAI)])

fwrite(HLAI_BA,"/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_all_BA.txt")

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
Netii_BA = melt(Netii, id.vars = 2, , measure.vars = seq(8,78,5))
Netii_BA_Score = melt(Netii, id.vars = 1, , measure.vars = seq(7,77,5))
Netii_EL_Rank = melt(Netii, id.vars = 1, , measure.vars = seq(6,76,5))
Netii_EL_Score = melt(Netii, id.vars = 1, , measure.vars = seq(5,75,5))

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

Netii = cbind(Netii_ID$ID, 1,Netii_Pep$Peptide, Netii_Rank[,2:3], Netii_log[,3],  Netii_nM[,3])
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

#nM: colnames(Netii)[c(seq(8,53,5))]
#EL score: colnames(Netii)[c(seq(5,50,5))]
#EL rank: colnames(Netii)[c(seq(6,51,5))]
#BA score: colnames(Netii)[c(seq(7,52,5))]
#BA rank: colnames(Netii)[c(seq(9,54,5))]

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
#Netii$Protein = sapply(strsplit(Netii$Protein, "_"),'[',1)
#Netii$Protein[which(Netii$Protein =="orflab")] = "orf1ab"



Netii$Haplotype = str_replace_all(Netii$Haplotype, pattern = "DRB1_", replacement = "HLA-DRB1*")
Netii$Haplotype = str_replace_all(Netii$Haplotype, pattern = "HLA-DQA1", replacement = "HLA-DQA1*")
Netii$Haplotype = str_replace_all(Netii$Haplotype, pattern = "-DQB1", replacement = "/DQB1*")
Netii$Haplotype = paste0(substr(Netii$Haplotype, 1, nchar(Netii$Haplotype)-2),":", substr(Netii$Haplotype, nchar(Netii$Haplotype)-1, nchar(Netii$Haplotype)))
Netii$Haplotype[which(substr(Netii$Haplotype,1,6) == "HLA-DQ")] = paste0(substr(Netii$Haplotype[which(substr(Netii$Haplotype,1,6) == "HLA-DQ")] ,1,11),":",
                                                                         substr(Netii$Haplotype[which(substr(Netii$Haplotype,1,6) == "HLA-DQ")] ,12,
                                                                                nchar(Netii$Haplotype[which(substr(Netii$Haplotype,1,6) == "HLA-DQ")])))

fwrite(Netii, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan4.0_standardized_mhc_all_viruses_all_lengths.txt")


# 
# ############NetMHCIIpan RNA virus##################
# 
# Netii = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCII_pan_RNA_virus.txt")
# #colnames(Netii)[c(seq(7,25,3), seq(32,54,3))] = substr(colnames(Netii)[c(seq(7,25,3), seq(32,54,3))],6, nchar(colnames(Netii)[c(seq(7,25,3), seq(32,54,3))]))
# Netii = Netii[order(Netii$Peptide),]
# colnames(Netii)[c(seq(6,49,3))] = colnames(Netii)[c(seq(5,49,3))]
# 
# Netii_Rank = melt(Netii, id.vars = 1, measure.vars = c(seq(6,50,3)))
# Netii_nM = melt(Netii, id.vars = 1, measure.vars = c(seq(5,49,3)))
# Netii_log = melt(Netii, id.vars = 1, measure.vars = c(seq(4,48,3)))
# Netii_ID = melt(Netii, id.vars = 3, measure.vars = c(seq(6,50,3)))
# Netii_Pep = melt(Netii, id.vars = 2, measure.vars = c(seq(6,50,3)))
# 
# Netii_Start = melt(Netii, id.vars = 1, measure.vars = c(seq(6,50,3)))
# 
# Netii = cbind(Netii_ID$ID, 1,Netii_Pep$Peptide, Netii_Rank[,2:3], Netii_log[,3],  Netii_nM[,3])
# colnames(Netii) = c("ID", "Start", "Peptide", "Haplotype", "Rank", "1-log50k", "Binding_affinity")
# #Netii$Protein = sapply(strsplit(Netii$Protein, "_"),'[',1)
# #Netii$Protein[which(Netii$Protein =="orflab")] = "orf1ab"
# 
# fwrite(Netii, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_RNA_virus.txt")


################################################
#############Figure S1.0A#########################

cutpoint = 500 #What nM is considered a binder/non-binder in IEDB data

######Compared performance vs IEDB##############
iedb = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/mhc_sars2_purified_competitive_radioactivity.csv")
iedb = iedb[,c(1:3)]

iedb_low = iedb[which(iedb$affinity < cutpoint),]
iedb_low = iedb_low[order(iedb_low$peptide),]


Flurry=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCflurry_standardized_all.txt")
Net_EL=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_all_EL.txt")
Net_BA=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_all_BA.txt")



iedb_unique = paste(iedb_low$peptide, iedb_low$allele, sep = "_")
Net_BA_unique = paste(Net_BA$Peptide, Net_BA$Haplotype, sep = "_")
Net_EL_unique = paste(Net_EL$Peptide, Net_EL$Haplotype, sep = "_")
Flurry_unique = paste(Flurry$Peptide, Flurry$Haplotype, sep = "_")
#Netii_unique = paste(Netii$Peptide, Netii$Haplotype, sep = "_")

Net_BA = Net_BA[which(Net_BA_unique %in% iedb_unique),]
Net_EL = Net_EL[which(Net_EL_unique %in% iedb_unique),]
Flurry = Flurry[which(Flurry_unique %in% iedb_unique),]

all = paste(Net_BA$Peptide, Net_BA$Haplotype, sep = "_")
iedb_low = iedb_low[which(iedb_unique %in% all),]

iedb_low$NetMHCpan_BA_Rank=NA
iedb_low$NetMHCpan_BA_log=NA
iedb_low$NetMHCpan_BA_BA=NA

iedb_low$NetMHCpan_EL_Rank=NA
iedb_low$NetMHCpan_EL_log=NA

iedb_low$MHCFlurry_Rank=NA
iedb_low$MHCFlurry_BA=NA
iedb_low$MHCFlurry_Proc=NA
iedb_low$MHCFlurry_ProS=NA

#iedb_low$NetMHCIIpan_Rank=NA
#iedb_low$NetMHCIIpan_log=NA
#iedb_low$NetMHCIIpan_BA=NA

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



iedb_melt = melt(iedb_low, id.vars = c(1:3), measure.vars = c(4:12))
colnames(iedb_melt)[4:5] = c("Program", "Rank")
iedb_melt = iedb_melt[complete.cases(iedb_melt),]

#iedb_melt = iedb_melt[which(iedb_melt$Program %in% c('NetMHCpan_BA_BA',  'NetMHCpan_BA_Rank', "NetMHCpan_EL_Rank", 'MHCFlurry_BA', 'MHCFlurry_Rank', "MHCFlurry_Proc", "MHCFlurry_ProS"))]
#iedb_melt = iedb_melt[which(iedb_melt$Program %in% c('NetMHCpan_BA_BA', 'MHCFlurry_BA'))]

fwrite(iedb_melt, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/IEDB_SARS2_HLAI_filtered_competitive_radioactivity.txt") 

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
  test = cor.test(sub_melt$affinity, sub_melt$Rank, method = "spearman")
  
  p= ggplot(data = sub_melt)+
    geom_point(aes(x = affinity, y= Rank), color = col, alpha=.6, show.legend = F)+
    #scale_color_viridis_d(name = "Software")+
    #theme(text=element_text(face="bold",size=20,colour="black")) +
    scale_y_log10()+#breaks = c(0.001, 0.01, .1, 1, 10, 100, 1000, 10000))+
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



################################################
#############Figure S1.1A#########################

######Compared performance vs IEDB##############
#MHC II, using all viruses instead of only SARS

Alleles_to_keep=c("HLA-DRB1*01:01", "HLA-DRB1*03:01", "HLA-DRB1*04:01", "HLA-DRB1*07:01", "HLA-DRB1*11:01", "HLA-DRB1*13:01", "HLA-DRB1*15:01",
                  "HLA-DQA1*01:02/DQB1*06:02", "HLA-DQA1*05:01/DQB1*02:01", "HLA-DQA1*02:01/DQB1*02:02", "HLA-DQA1*05:05/DQB1*03:01", 
                   "HLA-DQA1*01:01/DQB1*05:01", "HLA-DQA1*03:01/DQB1*03:02", "HLA-DQA1*03:03/DQB1*03:01", "HLA-DQA1*01:03/DQB1*06:03")

iedb = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/mhc_all_virus_combined.csv")
iedb = iedb[,c(1:3)]

iedb=iedb[which(iedb$allele %in% Alleles_to_keep),]

iedb_low = iedb[which(iedb$affinity < cutpoint),]
iedb_low = iedb_low[order(iedb_low$peptide),]

Netii_3.2=fread( "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_mhc_all_viruses_all_lengths.txt")
Netii_4.0=fread( "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan4.0_standardized_mhc_all_viruses_all_lengths.txt")


iedb_unique = paste(iedb_low$peptide, iedb_low$allele, sep = "_")
Netii_unique = paste(Netii_3.2$Peptide, Netii_3.2$Haplotype, sep = "_")

Netii_3.2 = Netii_3.2[which(Netii_unique %in% iedb_unique),]
Netii_4.0 = Netii_4.0[which(Netii_unique %in% iedb_unique),]


all = paste(Netii_4.0$Peptide, Netii_4.0$Haplotype, sep = "_")

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

for(z in unique(iedb_melt$Program)){
  print(z)
  if(z %in% unique(iedb_melt$Program)[1:3]){
    col = viridis(3)[1]
  }else if(z %in% unique(iedb_melt$Program)[4:8]){
    col = viridis(3)[2]
  }
  
  sub_melt = iedb_melt[which(iedb_melt$Program ==z),]
  test = cor.test(sub_melt$affinity, sub_melt$Value, method = "spearman")
  
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




####################################################
#########Figure 2A##################################

iedb_melt = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/IEDB_SARS2_HLAI_filtered_competitive_radioactivity.txt") 

###MHCI
iedb_BA = iedb_melt[which(iedb_melt$Program %in% c('NetMHCpan_BA_BA'))]

Net_BA=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_all_BA.txt")

iedb = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/mhc_sars2_purified_competitive_radioactivity.csv")
iedb$affinity = as.numeric(iedb$affinity)

iedb_unique = paste(iedb$peptide, iedb$allele, sep = "_")
Net_BA_unique = paste(Net_BA$Peptide, Net_BA$Haplotype, sep = "_")

iedb$NetMHCpan_BA_BA = NA
for(z in 1:length(iedb_unique)){
  if(iedb_unique[z] %in% Net_BA_unique){
    iedb$NetMHCpan_BA_BA[z] = Net_BA$Binding_affinity[which(Net_BA_unique == iedb_unique[z])]
  }
}

iedb_comp = iedb[complete.cases(iedb),]
iedb_comp = iedb_comp[which(iedb_comp$allele %in% Net_BA$Haplotype),]
#Net_BA$Haplotype = paste0(substr(Net_BA$Haplotype ,1,5), "*", substr(Net_BA$Haplotype ,6,10))

iedb_comp$allele = factor(iedb_comp$allele, levels = unique(Net_BA$Haplotype), ordered = T)
Net_BA$Haplotype = factor(Net_BA$Haplotype, levels = unique(Net_BA$Haplotype), ordered = T)



#####Figure S1B
spec_cutpoint = .9
cutpoint = 500

iedb_comp$Hits = ifelse(iedb_comp$affinity < cutpoint,0,1)
iedb_comp = iedb[complete.cases(iedb),]
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

# ggplot(data=plot_cutoff, aes(x=cutoff, y=FDR))+
#   geom_point()+
#   theme(text=element_text(face="bold",size=20,colour="black")) +
#   labs(x="Binding affinity (nM)", y= "False discovery rate", title = "NetMHCpan v4.0 vs IEDB performance\nSARS -- 500nM IEDB cutoff")+
#   geom_hline(yintercept = c(.01,.05,.1))
#   #geom_smooth(method='lm',formula= y~+x)

ggplot(data=plot_cutoff, aes(x=cutoff, y=Specificity))+
  geom_point()+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x="Binding affinity (nM)", y= "Specificity", title = "NetMHCpan v4.0 vs IEDB performance\nSARS peptides, IEDB measured affinity < 500nM")+
  geom_hline(yintercept = c(spec_cutpoint))+
  geom_vline(xintercept = HLAI_nM_cut)+
  annotate("text", label = spec_cutpoint, x = 0, y = spec_cutpoint, size = 10, hjust=0, vjust = 0)+
  annotate("text", label = HLAI_nM_cut, x = HLAI_nM_cut, y = min(plot_cutoff$Specificity), size = 10, hjust=0, vjust = 0, angle = 90)
  
# ggplot(data=plot_cutoff, aes(y=Sensitivity, x=1-Specificity))+
#   geom_point()+
#   theme(text=element_text(face="bold",size=20,colour="black")) +
#   scale_x_continuous(limits = c(0,1))+
#   scale_y_continuous(limits = c(0,1))
#   labs(x="Binding affinity (nM)", y= "Specificity", title = "NetMHCpan v4.0 vs IEDB performance\nSARS peptides, IEDB measured affinity < 500nM")+
#   geom_hline(yintercept = c(.9))+
#   geom_vline(xintercept = HLAI_nM_cut)


##########################Plotting 2A#############
scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

library(RColorBrewer)

p=ggplot(data=iedb_comp, aes(x=NetMHCpan_BA_BA, y=affinity, color = factor(allele)))+
  geom_point()+
  scale_color_manual(values = c(brewer.pal(7, 'Blues')[3:7],
                                brewer.pal(7, 'Reds')[3:7],
                                brewer.pal(8, 'Greens')[3:8]),
                     name="Allele", drop = F, guide = guide_legend(override.aes = list(color = "white")))+
  scale_y_log10(breaks = c(0.001, 0.01, .1, 1, 10, 100, 1000, 10000,100000), label=scientific_10)+
  scale_x_continuous(trans= 'log10', breaks = c(1,10,100,1000,10000,100000,1000000), label=scientific_10, limits= c(1,100000))+
  geom_vline(xintercept = HLAI_nM_cut)+
  geom_hline(yintercept = cutpoint)+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x="NetMHCpan BA (nM)", y="IEDB affinity (nM)")+
  annotate("text",x=5, y = .3, label = paste0("n=", length(which((iedb_comp$NetMHCpan_BA_BA<HLAI_nM_cut)&(iedb_comp$affinity<cutpoint)))), color="red", size=10, hjust=0)+
  annotate("text",x=5, y = 10000, label = paste0("n=", length(which((iedb_comp$NetMHCpan_BA_BA<HLAI_nM_cut)&(iedb_comp$affinity>cutpoint)))), color="red", size=10,hjust=0)+
  annotate("text",x=10000, y = .3, label = paste0("n=", length(which((iedb_comp$NetMHCpan_BA_BA>HLAI_nM_cut)&(iedb_comp$affinity<cutpoint)))), color="red", size=10,hjust=0)+
  annotate("text",x=10000, y = 10000, label = paste0("n=", length(which((iedb_comp$NetMHCpan_BA_BA>HLAI_nM_cut)&(iedb_comp$affinity>cutpoint)))), color="red", size=10,hjust=0)+
  
  annotate("text",x=2.5, y = cutpoint, label = paste0(cutpoint,"nM"), color="black", size=5,hjust=0, vjust=0)+
  annotate("text",x=HLAI_nM_cut, y = .1, label = paste0(HLAI_nM_cut,"nM"), color="black", size=5, angle=90, hjust=0, vjust=0)+
  theme(
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")
  ) 


hist_top <- ggplot(data=Net_BA)+
  geom_histogram(aes(x= Binding_affinity, fill = Haplotype), binwidth = .01)+
  scale_fill_manual(values = c(brewer.pal(7, 'Blues')[3:7],
                               brewer.pal(7, 'Reds')[3:7],
                               brewer.pal(8, 'Greens')[3:8]),
                    name="Allele", drop = F)+
  scale_x_continuous(trans= 'log10', breaks = c(1,10,100,1000,10000,100000,1000000), limits= c(1,100000))+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  annotate("text",x=5, y = 300000, label = paste0("n=", length(which(Net_BA$Binding_affinity<HLAI_nM_cut))), color="red", size=10,hjust=0)+
  annotate("text",x=10000, y = 300000, label = paste0("n=", length(which(Net_BA$Binding_affinity>HLAI_nM_cut))), color="red", size=10,hjust=0)+
  geom_vline(xintercept = HLAI_nM_cut)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(y="Count")+
  scale_y_continuous(trans = sqrt_trans(), label=scientific_10)

grid.arrange(hist_top, p, ncol=1, nrow=2, heights = c(3,4))  
# 
# 
#   
# library(precrec)
# precrec_obj <- evalmod(scores = iedb_comp$NetMHCpan_BA_BA, labels = iedb_comp$Hits)
# autoplot(precrec_obj)
# 
# ggplot()+
#   geom_histogram(data = Net_BA, aes(x=Binding_affinity, y = ..ncount.., color="All netMHCpan calls"), alpha=0.3, binwidth = 5)+
#   geom_histogram(data=iedb_BA, aes(x= Rank,y = ..ncount.., color="NetMHCpan calls with IEDB support"), alpha = 0.3, binwidth = 5)+
#   scale_color_viridis_d(name="Data Group")+
#   geom_vline(xintercept = HLAI_nM_cut)+
#   scale_x_continuous(trans=ssqrt_trans)+
#   theme(text=element_text(face="bold",size=20,colour="black")) +
#   labs(y="Relative proportion (scaled to max bin height)", x="NetMHCpan BA (nM)")


Net_BA_filt = Net_BA[which(Net_BA$Binding_affinity < HLAI_nM_cut),]

fwrite(Net_BA,"/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_BA.txt")
fwrite(Net_BA_filt, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_filtered_BA.txt")


################################################
##############Figure 2B#########################

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
iedb_comp = iedb_v[complete.cases(iedb_v),]

# ggplot(data=iedb_comp[which(nchar(iedb_comp$peptide) %ni%c(12,15)),], aes(x=NetMHCIIpan_BA, y=affinity, color = factor(nchar(peptide))))+
#   scale_color_viridis_d()+
#   geom_point()+
#   scale_y_log10()+#breaks = c(0.001, 0.01, .1, 1, 10, 100, 1000, 10000))+
#   scale_x_log10()+
#   geom_vline(xintercept = 500)+
#   geom_hline(yintercept = 500)+
#   theme(text=element_text(face="bold",size=20,colour="black")) +
#   labs(x="NetMHCIIpan BA (nM)", y="IEDB affinity (nM)")+
#   annotate("text",x=10, y = .3, label = paste0("n=", length(which((iedb_comp$NetMHCIIpan_BA<500)&(iedb_comp$affinity<500)))), color="red", size=10)+
#   annotate("text",x=10, y = 10000, label = paste0("n=", length(which((iedb_comp$NetMHCIIpan_BA<500)&(iedb_comp$affinity>500)))), color="red", size=10)
# 

###Figure S1D

spec_cutpoint = .9
cutpoint = 500
iedb_comp$Hits = ifelse(iedb_comp$affinity < cutpoint,0,1)

cutoff = seq(0,cutpoint,1)

FDR = c()
Specificity=c()
Sensitivity = c()

for(n in cutoff){
  print(n)
  FDR = c(FDR, length(which(iedb_comp$affinity[which(iedb_comp$NetMHCIIpan_BA<=n)] > cutpoint))/length(which(iedb_comp$NetMHCIIpan_BA<=n)))
  Sensitivity = c(Sensitivity,  length(which(iedb_comp$affinity[which(iedb_comp$NetMHCIIpan_BA<=n)] < cutpoint))/length(which(iedb_comp$affinity<cutpoint)) )
  Specificity = c(Specificity,  length(which(iedb_comp$affinity[which(iedb_comp$NetMHCIIpan_BA>n)] > cutpoint))/length(which(iedb_comp$affinity>cutpoint)) )
}


plot_cutoff = data.table(cutoff, FDR, Sensitivity, Specificity)
HLAII_nM_cut = plot_cutoff$cutoff[max(which(plot_cutoff$Specificity > spec_cutpoint))]

# ggplot(data=plot_cutoff, aes(x=cutoff, y=FDR))+
#   geom_point()+
#   theme(text=element_text(face="bold",size=20,colour="black")) +
#   labs(x="Binding affinity (nM)", y= "False discovery rate", title = "NetMHCIIpan v3.2 vs IEDB performance\nAll viruses -- 500nM IEDB cutoff")+
#   geom_hline(yintercept = c(.01,.05,.1))
ggplot(data=plot_cutoff, aes(x=cutoff, y=Specificity))+
  geom_point()+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x="Binding affinity (nM)", y= "Specificity", title = "NetMHCIIpan v3.2 vs IEDB performance\nViral peptides, IEDB measured affinity < 500nM")+
  geom_hline(yintercept = c(spec_cutpoint))+
  geom_vline(xintercept = HLAII_nM_cut)+
  annotate("text", label = spec_cutpoint, x = 0, y = spec_cutpoint, size = 10, hjust=0, vjust = 0)+
  annotate("text", label = HLAII_nM_cut, x = HLAII_nM_cut, y = min(plot_cutoff$Specificity), size = 10, hjust=0, vjust = 0, angle = 90)


# 
# library(precrec)
# precrec_obj <- evalmod(scores = iedb_comp$NetMHCIIpan_BA, labels = iedb_comp$Hits)
# autoplot(precrec_obj)
# 
# ggplot()+
#   geom_histogram(data = Netii, aes(x=Binding_affinity, y = ..ncount.., color="All netMHCIIpan calls"), alpha=0.3, binwidth = 2.5)+
#   geom_histogram(data=iedb_II, aes(x= Rank,y = ..ncount.., color="NetMHCIIpan calls with IEDB support"), alpha = 0.3, binwidth = 2.5)+
#   geom_histogram(data = Netii[which(Netii$Haplotype =="DRB1_0101")], aes(x=Binding_affinity, y = ..ncount.., color="DRB1*0101 netMHCIIpan calls"), alpha=0.3, binwidth = 2.5)+
#   
#   scale_color_viridis_d(name="Data Group")+
#   geom_vline(xintercept = 1000)+
#   scale_x_continuous(trans=ssqrt_trans)+
#   theme(text=element_text(face="bold",size=20,colour="black")) +
#   labs(y="Relative proportion (scaled to max bin height)", x="NetMHCIIpan BA (nM)")


#######Plotting figure 2B#####################

Netii_SARS=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_all_BA.txt")

Netii_SARS$Haplotype = factor(Netii_SARS$Haplotype, levels = unique(Netii_SARS$Haplotype), ordered = T)
iedb_comp$allele = factor(iedb_comp$allele, levels = unique(Netii_SARS$Haplotype), ordered = F)


scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

p2=ggplot(data=iedb_comp, aes(x=NetMHCIIpan_BA, y=affinity, color = factor(allele)))+
  geom_point()+
  scale_color_manual(values = c(brewer.pal(9, 'Blues')[3:9],
                                brewer.pal(9, 'Reds')[2:9]),
                     name="Allele", drop = F, guide = guide_legend(override.aes = list(color = "white")))+
  scale_y_log10(breaks = c(0.001, 0.01, .1, 1, 10, 100, 1000, 10000,100000), label=scientific_10)+
  scale_x_log10(breaks = c(1,10,100,1000,10000,100000,1000000), label=scientific_10, limits=c(5,100000))+
  geom_vline(xintercept = HLAII_nM_cut)+
  geom_hline(yintercept = cutpoint)+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x="NetMHCpan BA (nM)", y="IEDB affinity (nM)")+
  annotate("text",x=5, y = .3, label = paste0("n=", length(which((iedb_comp$NetMHCIIpan_BA<HLAII_nM_cut)&(iedb_comp$affinity<cutpoint)))), color="red", size=10, hjust=0)+
  annotate("text",x=5, y = 100000, label = paste0("n=", length(which((iedb_comp$NetMHCIIpan_BA<HLAII_nM_cut)&(iedb_comp$affinity>cutpoint)))), color="red", size=10,hjust=0)+
  annotate("text",x=10000, y = .3, label = paste0("n=", length(which((iedb_comp$NetMHCIIpan_BA>HLAII_nM_cut)&(iedb_comp$affinity<cutpoint)))), color="red", size=10,hjust=0)+
  annotate("text",x=10000, y = 100000, label = paste0("n=", length(which((iedb_comp$NetMHCIIpan_BA>HLAII_nM_cut)&(iedb_comp$affinity>cutpoint)))), color="red", size=10,hjust=0)+
  
  annotate("text",x=5, y = 350, label = paste0(cutpoint,"nM"), color="black", size=5,hjust=0)+
  annotate("text",x=315, y = .1, label = paste0(HLAII_nM_cut,"nM"), color="black", size=5, angle=90, hjust=0, vjust=0)+
  theme(
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")
  ) 


hist_top2 <- ggplot(data=Netii_SARS)+
  geom_histogram(aes(x= Binding_affinity, fill = Haplotype), binwidth = .01)+
  scale_fill_manual(values= c(brewer.pal(9, 'Blues')[3:9],
                              brewer.pal(9, 'Reds')[2:9]),
                    name="Allele", drop = F)+
  scale_x_log10(breaks = c(1,10,100,1000,10000,100000,1000000), limits=c(5,100000))+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_vline(xintercept = HLAII_nM_cut)+
  annotate("text",x=5, y = 800, label = paste0("n=", length(which(Netii_SARS$Binding_affinity<HLAII_nM_cut))), color="red", size=10,hjust=0)+
  annotate("text",x=10000, y = 800, label = paste0("n=", length(which(Netii_SARS$Binding_affinity>HLAII_nM_cut))), color="red", size=10,hjust=0)+
  labs(y="Count")+
  scale_y_continuous(labels = scientific_10, breaks = c(0,200,400,600,800,10000))

grid.arrange(hist_top2, p2, ncol=1, nrow=2, heights = c(3,4))



Netii_SARS_filt = Netii_SARS[which(Netii_SARS$Binding_affinity < HLAII_nM_cut),]

fwrite(Netii_SARS, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_BA.txt")
fwrite(Netii_SARS_filt, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_filtered_BA.txt")






#################################################
####################New Fig 2C ##################
Net_BA=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_BA.txt")
Netii_SARS=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_BA.txt")

Netii_SARS_filt=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_filtered_BA.txt")
Net_BA_filt=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_filtered_BA.txt")


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

Net_BA_cast = dcast.data.table(Net_BA, Peptide + Protein + Start ~ Haplotype, value.var = "Binding_affinity")
Net_BA_cast = Net_BA_cast[which(Net_BA_cast$Peptide %in% Net_BA_filt$Peptide),]

Netii_BA_cast = dcast.data.table(Netii_SARS, Peptide + Protein + Start ~ Haplotype, value.var = "Binding_affinity")
Netii_BA_cast = Netii_BA_cast[which(Netii_BA_cast$Peptide %in% Netii_SARS_filt$Peptide),]


########Calculating population frequencies

library(doMC)
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
####Calculating co-epitopes

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




#########Prepping data to plot
mouse1 = fread(paste0(WORKING_DIR, "COVID_murine_NetMHCpan_unfilt.xls"))
mouse1 = mouse1[which(mouse1$NB>0),]
mouse2 = fread(paste0(WORKING_DIR, "COVID_murine_NetMHCIIpan_unfilt.xls"))
mouse2 = mouse2[which(mouse2$NB>0),]

Net_coepitopes=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Coepitopes_human.txt")
Net_coepitopes$c1_freq = Net_coepitopes$c1_freq*100
Net_coepitopes$c2_freq = Net_coepitopes$c2_freq*100

###Number of MHC-I/II epitopes conserved between human and mouse
#length(which(HLAI_filt$Peptide[which(HLAI_filt$lo_entropy==1)] %in% unique(mouse1$Peptide)))
#length(which(HLAII_filt$Peptide[which(HLAII_filt$lo_entropy==1)] %in% unique(mouse2$Peptide)))

##Adding murine coverage for each human epitope
Net_coepitopes$murine = "None"
Net_coepitopes$murine[which(Net_coepitopes$c1_peptide %in% mouse1$Peptide)] = "MHC-I" 
Net_coepitopes$murine[which(Net_coepitopes$c2_peptide %in% mouse2$Peptide)] = "MHC-II" 
Net_coepitopes$murine[which((Net_coepitopes$c2_peptide %in% mouse2$Peptide) & (Net_coepitopes$c1_peptide %in% mouse1$Peptide))] = "Both" 
Net_coepitopes$murine = factor(Net_coepitopes$murine, levels = c("MHC-I", "MHC-II", "Both", "None"))

Net_coepitopes$Total_freq=Net_coepitopes$Total_freq*100

####Add in Tc literature data


paper_tc = fread(paste0(WORKING_DIR, "Supplemental/Standardized T cell epitopes.txt"))
paper_tc = paper_tc[!duplicated(paper_tc),]
seq =  fread(paste0(WORKING_DIR, "AA_sequence_combined.txt"), header=F)

paper_tc$Start = NA
paper_tc$End = NA


###Checking if string is present in SARS-CoV-2 reference.  If so, add coordinates####
for (y in 1:nrow(paper_tc)){
  if(str_detect(string = seq$V1 ,paper_tc$Peptide[y])){
    print(y)
    paper_tc$Start[y] = gregexpr(text = seq$V1, pattern = paper_tc$Peptide[y])[[1]][1]
    paper_tc$End[y] = (paper_tc$Start[y] + nchar(paper_tc$Peptide[y]) - 1)
  }
  
}
paper_tc = paper_tc[which(!is.na(paper_tc$Start)),]
paper_tc$Type = "T cell"


Net_coepitopes$Lit_overlap = 0
Net_coepitopes$Lit_overlap[which((Net_coepitopes$c1_peptide %in% paper_tc$Peptide)  | (Net_coepitopes$c2_peptide %in% paper_tc$Peptide))]=1


#### Transforming to proteome space


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


####Calculation fit for class I and II
class_I_coverage=c()
class_II_coverage = c()

for (n in 1:9701){
  class_I_coverage = c(class_I_coverage, length(which((Net_BA_cast$Start<=n) & ((Net_BA_cast$Start+nchar(Net_BA_cast$Peptide)-1)>=n)  )))
  class_II_coverage = c(class_II_coverage, length(which((Netii_BA_cast$Start<=n) & ((Netii_BA_cast$Start+nchar(Netii_BA_cast$Peptide)-1)>=n)  )))
}

coverage=rbind(data.table(seq(1,9701,1),"HLA-I",class_I_coverage),
               data.table(seq(1,9701,1),"HLA-II",class_II_coverage), use.names=F)
colnames(coverage) = c("Position", "Class", "Coverage")


####################
Min_frequency = 25
label_frequency = 45

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


sm= ggplot(data=coverage,aes(x = Position, y=Coverage, color = Class))+geom_smooth( method = 'loess',span=.01)+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  scale_y_continuous(breaks = c(0,4,8,12))+
  labs(x="Position across SARS-CoV-2 proteome", y="Count")

png("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Figure_2C.png", width = 20, height = 8.3, units = "in", res=300)
grid.arrange(p, sm, ncol=1, nrow=2, heights = c(5,2))
dev.off()


####################################
#####Figure 2D######################

###Class I/II overlap

library(venneuler)

vd_I <- venneuler(c(HLAI=nrow(Net_BA_cast), I_coep=0,"HLAI&I_coep"=length(unique(Net_coepitopes$c1_peptide))))
vd_I$labels=NA
#vd$colors 

plot(vd_I, col= c(alpha(colour = "orange", .5), alpha(colour = "red",.5)))
text(0.825,.5,paste0("All human\nHLA-I epitopes\n",nrow(Net_BA_cast)))
text(0.425,.5,paste0("Human\nHLA-I co-epitopes\n", length(unique(Net_coepitopes$c1_peptide))))
#text(0.535,.5,"Co-epitopes: 6289\nHLA-I: 2721\nHLA-II: 2721")



vd_II <- venneuler(c(HLAII=nrow(Netii_BA_cast), II_coep=0,"HLAII&II_coep"=length(unique(Net_coepitopes$c2_peptide))))
vd_II$labels=NA
#vd_II$centers = vd_I$centers
#vd$colors 

plot(vd_II, col= c(alpha(colour = "cyan", .5), alpha(colour = "blue",.5)))
text(0.825,.5,paste0("All human\nHLA-II epitopes\n",nrow(Netii_BA_cast)))
text(0.425,.5,paste0("Human\nHLA-II co-epitopes\n",length(unique(Net_coepitopes$c2_peptide))))
#text(0.535,.5,"Co-epitopes: 6289\nHLA-I: 2721\nHLA-II: 2337")

vd_CO <- venneuler(c(Coep=nrow(Net_coepitopes), II_coep=100,"Coep&II_coep"=0))
vd_CO$labels=NA
#vd_II$centers = vd_I$centers
#vd$colors 

plot(vd_CO, col= c(alpha(colour = "purple", .5), alpha(colour = "blue",.5)))
text(0.75,.55,paste0("Co-epitopes: ", nrow(Net_coepitopes)))

#text(0.535,.5,"Co-epitopes: 6289\nHLA-I: 2721\nHLA-II: 2337")

####################################
#####Figure 2E######################

###Murine overlap

library(venneuler)

murine_1_overlap = length(which(mouse1$Peptide %in% Net_BA_cast$Peptide))
murine_2_overlap = length(which(mouse2$Peptide %in% Netii_BA_cast$Peptide))
murine_1_coep = length(which(mouse1$Peptide %in% Net_coepitopes$c1_peptide))
murine_2_coep = length(which(mouse2$Peptide %in% Net_coepitopes$c2_peptide))
coep_1_overlap = length(which(Net_BA_cast$Peptide %in% Net_coepitopes$c1_peptide))
coep_2_overlap = length(which(Netii_BA_cast$Peptide %in% Net_coepitopes$c2_peptide))


vd_all <- venneuler(c(HLAI=(nrow(Net_BA_cast)-murine_1_overlap-coep_1_overlap), MHCI=(nrow(mouse1)-murine_1_overlap-murine_1_coep), 
                      HLAII=(nrow(Netii_BA_cast)-murine_2_overlap-coep_1_overlap), MHCII=(nrow(mouse2)-murine_2_overlap-murine_2_coep),
                      Coepitopes=(nrow(Net_coepitopes)-murine_1_coep- murine_2_coep-coep_1_overlap-coep_2_overlap),
                      "HLAI&MHCI"=murine_1_overlap, "HLAII&MHCII"=murine_2_overlap,
                      "Coepitopes&MHCI" = murine_1_coep,
                      "Coepitopes&MHCII" = murine_2_coep,
                      "HLAI&Coepitopes" = coep_1_overlap,
                      "HLAII&Coepitopes" = coep_2_overlap))


vd_1 = venneuler(c(HLAI=(nrow(Net_BA_cast)-murine_1_overlap), MHCI=(nrow(mouse1)-murine_1_overlap),
                   "HLAI&MHCI"=murine_1_overlap))
vd_1$labels<- c(
  paste0("Human\nHLA-I\n",(nrow(Net_BA_cast)-murine_1_overlap)),
  paste("Murine\nMHC-I\n",(nrow(mouse1)-murine_1_overlap))
)

plot(vd_1, col = c(alpha("red",.5), alpha("tomato2",.5)))
text(0.5,.5,paste0(murine_1_overlap))


####
vd_2 = venneuler(c(HLAII=(nrow(Netii_BA_cast)-murine_2_overlap), MHCII=(nrow(mouse2)-murine_2_overlap),
                   "HLAII&MHCII"=murine_2_overlap))
vd_2$labels = NA


plot(vd_2, col=c(alpha("blue",.5), alpha("cadetblue",.5)))
text(0.825,.5,paste0("Human\nHLA-II\n", (nrow(Netii_BA_cast)-murine_2_overlap)))
text(0.15,.5,paste0("Murine\nMHC-II\n", (nrow(mouse2)-murine_2_overlap)))
text(0.425,.5,murine_2_overlap)


###
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



#####################################################
########Figure 2F####################################
#Protein dist, normalized by length
all_prots=c("orf1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF8","N","ORF10")
prot_len= c(abs(1 - 7096),abs(7097 - 8369),abs(8370 - 8644),abs(8645 - 8719),abs(8720 - 8941),abs(8942 - 9002),abs(9003 - 9123), abs(9124 - 9244),abs(9245 - 9663),abs(9664 - 9701))

Net_coepitopes = Net_coepitopes[order(Net_coepitopes$c2_start),]
prots = factor(c(unique(Net_coepitopes$Protein)), levels = all_prots)

##Raw counts
counts=c()
for(z in prots){
  counts = c(counts, length(which(Net_coepitopes$Protein == z)))#/ prot_len[which(all_prots==z)] )
}

prot_dist = data.table(prots, counts)
colnames(prot_dist) = c("Protein", "Count")

library(dplyr)

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

library(dplyr)

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

############################################
###Figure 2G - Histogram of pop freq###################

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
  #scale_color_viridis_d(begin = 0, end = .75, option = 'plasma')+ 
  #scale_fill_viridis_d(begin = 0, end = .75, option = 'plasma')+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(y="Count") + 
  scale_y_continuous(limits = c(0,2000))+
  scale_x_continuous(breaks = seq(0,100,10))
  
  # 
  # annotate("text", x=common_HLAI[1], y = 100, label = HLAI_lab1, color = "black", hjust=0, vjust=0, angle = 90)+
  # annotate("text", x=common_HLAI[2], y = 100, label = HLAI_lab2, color = "black", hjust=0, vjust=0, angle = 90)+
  # annotate("text", x=common_HLAI[3], y = 100, label = HLAI_lab3, color = "black", hjust=0, vjust=0, angle = 90)+
  # annotate("text", x=common_HLAI[4], y = 100, label = HLAI_lab4, color = "black", hjust=0, vjust=0, angle = 90)+
  # 
  # annotate("text", x=common_HLAII[1], y = 100, label = HLAII_lab1, color = "black", hjust=0, vjust=0, angle = 90)+
  # annotate("text", x=common_HLAII[2], y = 100, label = HLAII_lab2, color = "black", hjust=0, vjust=0, angle = 90)+
  # annotate("text", x=common_HLAII[3], y = 100, label = paste0(HLAII_lab3[1],"\n",HLAII_lab3[2]), color = "black", hjust=0, vjust=0, angle = 90)+
  # annotate("text", x=common_HLAII[4], y = 100, label = HLAII_lab4, color = "black", hjust=0, vjust=0, angle = 90)
  # 





####################################################
#########Figure 3A top --All candidates#############

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

library(dplyr)

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

####Finding corresponding protein

iedb_fig3 = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/mhc_sars2_purified_competitive_radioactivity.csv")
iedb_fig3 = iedb_fig3[which(iedb_fig3$affinity < cutpoint),]
iedb_fig3 = iedb_fig3[which(iedb_fig3$allele %in% Freq$Haplotype),]


seq =  fread(paste0(WORKING_DIR, "AA_sequence_combined.txt"), header=F)

iedb_fig3$Start = NA
iedb_fig3$End = NA
iedb_fig3$Protein = NA

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
all_prots=c("orf1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF8","N","ORF10")
prot_len= c(abs(1 - 7096),abs(7097 - 8369),abs(8370 - 8644),abs(8645 - 8719),abs(8720 - 8941),abs(8942 - 9002),abs(9003 - 9123), abs(9124 - 9244),abs(9245 - 9663),abs(9664 - 9701))

prots = factor(c(unique(Net_coepitopes$Protein)), levels = all_prots)

iedb_fig3_I = iedb_fig3[which(iedb_fig3$mhc_class ==  "I"),]
iedb_fig3_II = iedb_fig3[which(iedb_fig3$mhc_class ==  "II"),]

iedb_fig3_I = dcast.data.table(iedb_fig3_I, peptide+Protein+Start+End~allele, value.var = "affinity")
iedb_fig3_II = dcast.data.table(iedb_fig3_II, peptide+Protein+Start+End~allele, value.var = "affinity")


counts_iedb_I = c()
counts_iedb_II = c()

for(z in prots){
  counts_iedb_I = c(counts_iedb_I, length(which(iedb_fig3_I$Protein == z)))#/ prot_len[which(all_prots==z)] )
  counts_iedb_II = c(counts_iedb_II, length(which(iedb_fig3_II$Protein == z)))#/ prot_len[which(all_prots==z)] )
  
}

prot_dist = data.table(prots, counts_iedb_I, counts_iedb_II)
colnames(prot_dist)[1] = "Protein"

library(dplyr)

###Right, IEDB data

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




####Combine, filter
ent = fread(paste0(WORKING_DIR, "entropy_7882.txt"))

all_I = unique(c(iedb_fig3_I$peptide, Net_BA_cast$Peptide))
all_II = unique(c(iedb_fig3_II$peptide, Netii_BA_cast$Peptide))

Net_combined = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAI_all_glm.txt")
Netii_combined = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAII_all_glm.txt")


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


################################################################
#################Figure 3 middle -- Filtering###################


##No filter
Pop_freq_I_iedb = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAI_vaccine_candidates_prefilter.txt")
Pop_freq_II_iedb = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAII_vaccine_candidates_prefilter.txt")
Net_coepitopes = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Coepitopes_vaccine_candidates_prefilter.txt")

nrow(Pop_freq_I_iedb)
nrow(Pop_freq_II_iedb)
nrow(Net_coepitopes)

##Filter by tetramer models
MHCI_tet =fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAI_all_glm.txt")
MHCII_tet =fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAII_all_glm.txt")

filt_I = Pop_freq_I_iedb[which(Pop_freq_I_iedb$all_I.n. %in% MHCI_tet$Peptide[which(MHCI_tet$Predicted_tetramer == TRUE)]),]
filt_II = Pop_freq_II_iedb[which(Pop_freq_II_iedb$all_II.n. %in% MHCII_tet$Peptide[which(MHCII_tet$Predicted_tetramer == TRUE)]),]
filt_co = Net_coepitopes[which((Net_coepitopes$c1_peptide %in% MHCI_tet$Peptide[which(MHCI_tet$Predicted_tetramer == TRUE)]) & 
                                 (Net_coepitopes$c2_peptide %in% MHCII_tet$Peptide[which(MHCII_tet$Predicted_tetramer == TRUE)])), ]

nrow(filt_I)
nrow(filt_II)
nrow(filt_co)

##Filter by entropy
filt_I = filt_I[which(filt_I$Low_entropy==1),]
filt_II = filt_II[which(filt_II$Low_entropy==1),]
filt_co = filt_co[which(filt_co$Low_entropy==1),]

nrow(filt_I)
nrow(filt_II)
nrow(filt_co)

#Filter by protein
filt_I = filt_I[which(filt_I$Protein %in% c("orf1ab", "S", "N")),]
filt_II = filt_II[which(filt_II$Protein %in% c("orf1ab", "S", "N")),]
filt_co = filt_co[which(filt_co$Protein %in% c("orf1ab", "S", "N")),]

nrow(filt_I)
nrow(filt_II)
nrow(filt_co)

##Filter by murine coverage
mouse1 = fread(paste0(WORKING_DIR, "COVID_murine_NetMHCpan_unfilt.xls"))
mouse1 = mouse1[which(mouse1$NB>0),]
mouse2 = fread(paste0(WORKING_DIR, "COVID_murine_NetMHCIIpan_unfilt.xls"))
mouse2 = mouse2[which(mouse2$NB>0),]


filt_I$Murine = NA
filt_II$Murine = NA
filt_co$Murine = NA

filt_I$Murine[which(filt_I$all_I.n. %in% mouse1$Peptide)] = "MHC-I"
filt_II$Murine[which(filt_II$all_II.n. %in% mouse2$Peptide)] = "MHC-II"
filt_co$Murine[which(filt_co$c1_peptide %in% mouse1$Peptide)] = "MHC-I"
filt_co$Murine[which(filt_co$c2_peptide %in% mouse2$Peptide)] = "MHC-II"
filt_co$Murine[which((filt_co$c1_peptide %in% mouse1$Peptide)&(filt_co$c2_peptide %in% mouse2$Peptide))] = "Both"

filt_I=filt_I[which(!is.na(filt_I$Murine)),]
filt_II=filt_II[which(!is.na(filt_II$Murine)),]
filt_co=filt_co[which(!is.na(filt_co$Murine)),]


nrow(filt_I)
nrow(filt_II)
nrow(filt_co)

##Filter by freq>50%
filt_I = filt_I[which(filt_I$Frequency>.5),]
filt_II = filt_II[which(filt_II$Frequency>.5),]
filt_co = filt_co[which((filt_co$c1_freq>50)&(filt_co$c2_freq>50)),]

nrow(filt_I)
nrow(filt_II)
nrow(filt_co)

filt_I$Frequency = filt_I$Frequency*100
filt_II$Frequency = filt_II$Frequency*100

filt_I$Group = "HLAI"
filt_II$Group = "HLAII"
filt_co$Group = "Coepitope"
filt_co$End = filt_co$Start+nchar(filt_co$c2_peptide)-1


Master_tab=rbind(filt_I[,c(1,3,4,5,6,8,9)],
                 filt_II[,c(1,3,4,5,6,8,9)], use.names=F)
colnames(Master_tab)[1] = "Peptide"
###############################################
###Plot combined, Figure 3 bottom##############

y_ll = min(Master_tab$Frequency)
y_ul=90
label_frequency=100

p=ggplot(data=Master_tab) + 
  geom_segment(aes(x=Start, xend=Start, y=y_ll, yend=Frequency, color = Group), alpha=0.3, size=.5) +
  scale_color_manual(values = c("purple", "red", "blue"), name="HLA")+
  ggnewscale::new_scale("color")+
  
  geom_point(aes(x=Start, y=Frequency, color=Murine), alpha=.3, size = 5) +
  scale_color_manual(values = c("purple", "red", "blue"), name="Murine\noverlap")+
  scale_y_continuous(limits = c((y_ll-13),y_ul), breaks = c(20,30,40,50,60,70,80,90,100))+
  
  geom_label_repel(aes(label=ifelse(Frequency>label_frequency,paste0(as.character(Peptide)),''), x=Start, y=Frequency),
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
  labs(y="Population frequency")+
  theme(text=element_text(face="bold",size=20,colour="black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 


#sm= ggplot(data=coverage,aes(x = Position, y=Coverage, color = Class))+geom_smooth( method = 'loess',span=.1)+
#  theme(text=element_text(face="bold",size=20,colour="black")) +
#  scale_y_continuous(breaks = c(0,4,8,12))+
#  labs(x="Position across SARS-CoV-2 proteome", y="Count")

png("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Figure_3_bottom.png", width = 20, height = 8.3, units = "in", res=300)
p
dev.off()


#######Try centipede plot###############
Coep_tab = Master_tab[which(Master_tab$Group =="Coepitope"),]
Master_tab = Master_tab[which(Master_tab$Group !="Coepitope"),]


###Find overlapping class I and II peptides
g1 <-  GRanges(seqnames="COVID",
               IRanges(start=Master_tab[which(Master_tab$Group == "HLAI")]$Start,
                       end=Master_tab[which(Master_tab$Group == "HLAI")]$End), 
               I_peptide = Master_tab[which(Master_tab$Group == "HLAI")]$Peptide, 
               Frequency_I= Master_tab[which(Master_tab$Group == "HLAI")]$Frequency)

g2 <-  GRanges(seqnames="COVID",
               IRanges(start=Master_tab[which(Master_tab$Group == "HLAII")]$Start,
                       end=Master_tab[which(Master_tab$Group == "HLAII")]$End), 
               II_peptide = Master_tab[which(Master_tab$Group == "HLAII")]$Peptide, 
               Frequency_II= Master_tab[which(Master_tab$Group == "HLAII")]$Frequency)

HLA_I_II=as.data.table(mergeByOverlaps(g1,g2))

###Add overlap metric

HLA_I_II$Overlap = NA
for(z in 1:nrow(HLA_I_II)){
  c1 = c(HLA_I_II$g1.start[z]:HLA_I_II$g1.end[z])
  c2 = c(HLA_I_II$g2.start[z]:HLA_I_II$g2.end[z])
  HLA_I_II$Overlap[z] = 100*(length(which(c1 %in% c2))/length(c1))
  
}
HLA_I_II$Frequency = HLA_I_II$Frequency_I/100*HLA_I_II$Frequency_II/100

###


g1_merged <-  GRanges(seqnames="COVID",
               IRanges(start= HLA_I_II$g1.start,
                       end=HLA_I_II$g1.end), 
               I_peptide = HLA_I_II$g1.I_peptide, 
               Frequency_I= HLA_I_II$Frequency_I)
g1_red = reduce(g1_merged)
g2_merged <-  GRanges(seqnames="COVID",
                     IRanges(start= HLA_I_II$g2.start,
                             end=HLA_I_II$g2.end), 
                             II_peptide = HLA_I_II$g2.II_peptide, 
                             Frequency_II= HLA_I_II$Frequency_II)
                                     
g2_red = reduce(g2_merged)               

red_merged=as.data.table(mergeByOverlaps(g1_red, g2_red))

###Find best frequency and best overlap per merged group

merge_tab = data.table()
for(z in 1:nrow(red_merged)){
  print(z)
  rows = which((HLA_I_II$g1.start >= red_merged$g1_red.start[z]) & (HLA_I_II$g1.end <= red_merged$g1_red.end[z]) )
  
  best_frequency = max(HLA_I_II$g1.Frequency_I[rows]/100 * HLA_I_II$g2.Frequency_II[rows]/100)
  best_row = rows[which((HLA_I_II$g1.Frequency_I[rows]/100 * HLA_I_II$g2.Frequency_II[rows]/100) == best_frequency)]
  
  overlap = c()
  for(r in best_row){
    c1 = c(HLA_I_II$g1.start[r]:HLA_I_II$g1.end[r])
    c2 = c(HLA_I_II$g2.start[r]:HLA_I_II$g2.end[r])
    overlap = c(overlap, length(which(c1 %in% c2))/length(c1) )
  }
  
  best_overlap = max(overlap)
  best_c1 = paste0(unique(HLA_I_II$g1.I_peptide[best_row]), collapse = ", ")
  best_c2 = paste0(unique(HLA_I_II$g2.II_peptide[best_row]), collapse = ", ")
  best_start = min(HLA_I_II$g2.start[best_row])
  
  merge_tab = rbind(merge_tab, data.table(best_c1, best_c2, best_start, best_overlap, best_frequency))
  
}


####Define scales

library(scales)
trans <- function(x) {
  ifelse((x >= 50), x - 40, ifelse((x <= -50), x+40 , ifelse(x==0,x,x/10)))
}
inv <- function(x) {
  ifelse(x > 0.02, x + 40, ifelse(x < -0.02, x-40 , x*10))
}
my_trans <- trans_new("my_trans", trans, inv)

Master_tab$Frequency[which(Master_tab$Group =="HLAII")] = -1*Master_tab$Frequency[which(Master_tab$Group =="HLAII")] 


####Make plot
y_ll=-90
centipede = ggplot(data=Master_tab) + 
  
  geom_hline(yintercept = c(-80,-70,-60,-50,50,60,70,80), color = "white")+
  annotate("text", x = rep(-100,8), y =  c(-80,-70,-60,-50,50,60,70,80), label= c('80','70','60','50','50',"60","70","80"))+
  
  geom_segment(aes(x=Start, xend=Start, y=ifelse(Master_tab$Group=="HLAI",50,-50), yend=Frequency, color = Group, alpha = abs(Frequency)), size=.2, show.legend = F) +
  geom_point(aes(x=Start, y=Frequency, color=Group, size = abs(Frequency), alpha = abs(Frequency)),show.legend = T) +
  
  scale_y_continuous(breaks = c(),trans = my_trans)+
  scale_alpha_continuous(range = c(.25,.5), guide=F)+
  
  scale_color_manual(values = c("red", "blue"), name="HLA", guide=F)+
  scale_size_continuous(name = "HLA-I/II\nPop. freq.", limits = c(50,90), range = c(1,10))+
  ggnewscale::new_scale("color")+
  ggnewscale::new_scale("size")+
  ggnewscale::new_scale("alpha")+
  

  ###Try using all merged data
  geom_point(data=HLA_I_II,aes(x=g2.start, y=0, size = 100*Frequency , color = Overlap, alpha = 100*Frequency),
             show.legend = T,  position=position_jitter(width=0, height=7.5))+
  scale_color_viridis_c(direction = 1)+
  scale_size_continuous(name = "Overlaps\nPop. freq.", limits = c(25,60), range = c(1,10))+
  scale_alpha_continuous(range = c(.25,.75), guide = F)+
  labs(y = "Population frequency", x="Position along SARS-CoV-2 proteome")+
  annotate("text", x = c(-300,-300,-300), y = c(-65,0,65), label = c("HLA-II", "Overlap", "HLA-I"), angle =90, size = 5)+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  
  
  #geom_point(data=filt_co,aes(x=c2_start, y=0, size = Total_freq),color = "purple", alpha=.1,show.legend = T)
  
  
  # ###Try using merged_tab
  # geom_point(data=merge_tab,aes(x=best_start, y=0, size = 100*best_frequency , color = best_overlap),show.legend = T)+
  # scale_color_viridis_c(direction = -1)+
  # scale_size_continuous(name = "Pop frequency", limits = c(20,85), range = c(3,20))+
  # 
  # geom_point(data=filt_co,aes(x=c2_start, y=0, size = Total_freq),color = "purple", alpha=.1,show.legend = T)+

# geom_label_repel(data=Coep_tab,aes(label=ifelse(Frequency>40,paste0(as.character(Peptide)),''), x=Start, y=0),
#                   box.padding   = 0.35, 
#                   point.padding = 0.5,
#                   segment.color = 'grey50')+
  
geom_rect(aes(xmin=0, xmax=7095, ymin=(y_ll), ymax=(y_ll-12)), fill=viridis(10)[1], color="black", size=0.1) +
  geom_text(aes(x=7095/2, y=(y_ll-15), label="orf1ab", angle=0), size=5, color = "black") +
  geom_text(aes(x=1000, y=(y_ll), label="-1000", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=2000, y=(y_ll), label="-2000", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=3000, y=(y_ll), label="-3000", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=4000, y=(y_ll), label="-4000", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=5000, y=(y_ll), label="-5000", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=6000, y=(y_ll), label="-6000", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=7000, y=(y_ll), label="-7000", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  
  geom_rect(aes(xmin=7096, xmax=8368, ymin=(y_ll), ymax=(y_ll-12)), fill=viridis(10)[2], color="black", size=0.1) +
  geom_text(aes(x=7096+(8368-7096)/2, y=(y_ll-15), label="S", angle=0), size=5, color = "black") + 
  geom_text(aes(x=7096+200, y=(y_ll), label="-200", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=7096+500, y=(y_ll), label="-500", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=7096+800, y=(y_ll), label="-800", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=7096+1100, y=(y_ll), label="-1100", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  
  
  geom_rect(aes(xmin=8369, xmax=8644, ymin=(y_ll-0), ymax=(y_ll-12)), fill=viridis(10)[3], color="black", size=0.1) +
  geom_text(aes(x=8368+(8644-8368)/2, y=(y_ll-15), label="ORF3a", angle=0), size=2.5, color = "black") +
  geom_text(aes(x=8369+100, y=(y_ll), label="-100", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=8369+200, y=(y_ll), label="-200", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  
  geom_rect(aes(xmin=8645, xmax=8718, ymin=(y_ll), ymax=(y_ll-12)), fill=viridis(10)[4], color="black", size=0.1) +
  geom_text(aes(x=8645+(8718-8645)/2, y=(y_ll-15)), label="E", angle=0, size=2.5, color = "black") + 
  
  geom_rect(aes(xmin=8719, xmax=8940, ymin=(y_ll-0), ymax=(y_ll-12)), fill=viridis(10)[5], color="black", size=0.1) +
  geom_text(aes(x=8719+(8940-8719)/2, y=(y_ll-15)), label="M", angle=0, size=2.5, color = "black") +
  
  geom_rect(aes(xmin=8941, xmax=9001, ymin=(y_ll), ymax=(y_ll-12)), fill=viridis(10)[6], color="black", size=0.1) +
  geom_text(aes(x=8941+(9001-8941)/2, y=(y_ll-18)), label="ORF6", angle=0, size=2.5, color = "black") +
  
  geom_rect(aes(xmin=9002, xmax=9122, ymin=(y_ll-0), ymax=(y_ll-12)), fill=viridis(10)[7], color="black", size=0.1) +
  geom_text(aes(x=9002+(9122-9002)/2, y=(y_ll-21)), label="ORF7a", angle=0, size=2.5, color = "black") +
  
  geom_rect(aes(xmin=9123, xmax=9243, ymin=(y_ll), ymax=(y_ll-12)), fill=viridis(10)[8], color="black", size=0.1) +
  geom_text(aes(x=9123+(9243-9123)/2, y=(y_ll-18)), label="ORF8", angle=0, size=2.5, color = "black") +
  
  geom_rect(aes(xmin=9244, xmax=9662, ymin=(y_ll-0), ymax=(y_ll-12)), fill=viridis(10)[9], color="black", size=0.1) +
  geom_text(aes(x=9244+(9662-9244)/2, y=(y_ll-15)), label="N", angle=0, size=2.5, color = "black") +
  geom_text(aes(x=9244+100, y=(y_ll), label="-100", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=9244+200, y=(y_ll), label="-200", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  geom_text(aes(x=9244+300, y=(y_ll), label="-300", angle=-90), size=3.5, color = "white",  fontface = "plain", hjust=0) +
  
  geom_rect(aes(xmin=9663, xmax=9701, ymin=(y_ll), ymax=(y_ll-12)), fill=viridis(10)[10], color="black", size=0.1) +
  geom_text(aes(x=9663+(9701-9663)/2, y=(y_ll-15)), label="ORF10", angle=0, size=2.5, color = "black") 




png("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Figure_3_bottom.png", width = 20, height = 8.3, units = "in", res=300)
centipede
dev.off()
#########################################


# 
# ####MHCII by nM
# 
# 
# iedb_II = iedb_melt[which(iedb_melt$Program %in% c('NetMHCIIpan_BA'))]
# 
# Netii=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_all_BA.txt")
# 
# 
# iedb = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/mhc_sars2_purified_competitive_radioactivity.csv")
# iedb$allele = str_replace_all(iedb$allele, pattern = "\\*", replacement = "")
# iedb$allele = str_replace_all(iedb$allele, pattern = "HLA-DRB101:01", replacement = "DRB1_0101")
# iedb$affinity = as.numeric(iedb$affinity)
# 
# iedb_unique = paste(iedb$peptide, iedb$allele, sep = "_")
# Netii_unique = paste(Netii$Peptide, Netii$Haplotype, sep = "_")
# 
# iedb$NetMHCIIpan_BA = NA
# for(z in 1:length(iedb_unique)){
#   if(iedb_unique[z] %in% Netii_unique){
#     iedb$NetMHCIIpan_BA[z] = Netii$Binding_affinity[which(Netii_unique == iedb_unique[z])]
#   }
# }
# iedb_comp = iedb[complete.cases(iedb),]
# 
# 
# ggplot(data=iedb_comp, aes(x=NetMHCIIpan_BA, y=affinity))+
#   geom_point()+
#   scale_y_log10()+#breaks = c(0.001, 0.01, .1, 1, 10, 100, 1000, 10000))+
#   scale_x_log10()+
#   geom_vline(xintercept = 1000)+
#   geom_hline(yintercept = 1000)+
#   theme(text=element_text(face="bold",size=20,colour="black")) +
#   labs(x="NetMHCIIpan BA (nM)", y="IEDB affinity (nM)")+
#   annotate("text",x=10, y = .3, label = paste0("n=", length(which((iedb_comp$NetMHCIIpan_BA<1000)&(iedb_comp$affinity<1000)))), color="red", size=10)+
#   annotate("text",x=10, y = 10000, label = paste0("n=", length(which((iedb_comp$NetMHCIIpan_BA<1000)&(iedb_comp$affinity>1000)))), color="red", size=10)
# 
# iedb_comp$Hits = ifelse(iedb_comp$affinity < 1000,0,1)
# cut = 1000
# cutoff = seq(0,1000,1)
# 
# FPR = c()
# FDR = c()
# Specificity=c()
# Sensitivity = c()
# PPV = c()
# NPV = c()
# 
# for(n in cutoff){
#   print(n)
#   #FPR = c(FPR, length(which(iedb_comp$affinity[which(iedb_comp$NetMHCIIpan_BA<=n)] > cut))/length(which(iedb_comp$affinity>cut)))
#   FDR = c(FDR, length(which(iedb_comp$affinity[which(iedb_comp$NetMHCIIpan_BA<=n)] > cut))/length(which(iedb_comp$NetMHCIIpan_BA<=n)))
#   Sensitivity = c(Sensitivity,  length(which(iedb_comp$affinity[which(iedb_comp$NetMHCIIpan_BA<=n)] < cut))/length(which(iedb_comp$affinity<cut)) )
#   #Specificity = c(Specificity,  length(which(iedb_comp$affinity[which(iedb_comp$NetMHCIIpan_BA>n)] > cut))/length(which(iedb_comp$affinity>cut)) )
#   #PPV = c(PPV,  length(which(iedb_comp$affinity[which(iedb_comp$NetMHCIIpan_BA<=n)] < cut))/length(which(iedb_comp$NetMHCIIpan_BA<=n)) )
#   #NPV = c(NPV,  length(which(iedb_comp$affinity[which(iedb_comp$NetMHCIIpan_BA>n)] > cut))/length(which(iedb_comp$NetMHCIIpan_BA>n)) )
#   
# }
# 
# plot_cutoff = data.table(cutoff, FDR, Sensitivity)
# ggplot(data=plot_cutoff, aes(x=cutoff, y=FDR))+
#   geom_point()+
#   #scale_x_continuous(limits = c(0,1000))+
#   geom_hline(yintercept = c(.01,.05,.1))+
#   geom_smooth(method='lm',formula= y~+log(x))
# 
# 
# library(precrec)
# precrec_obj <- evalmod(scores = iedb_comp$NetMHCIIpan_BA, labels = iedb_comp$Hits)
# autoplot(precrec_obj)
# 
# ggplot()+
#   geom_histogram(data = Netii, aes(x=Binding_affinity, y = ..ncount.., color="All netMHCIIpan calls"), alpha=0.3, binwidth = 2.5)+
#   geom_histogram(data=iedb_II, aes(x= Rank,y = ..ncount.., color="NetMHCIIpan calls with IEDB support"), alpha = 0.3, binwidth = 2.5)+
#   geom_histogram(data = Netii[which(Netii$Haplotype =="DRB1_0101")], aes(x=Binding_affinity, y = ..ncount.., color="DRB1*0101 netMHCIIpan calls"), alpha=0.3, binwidth = 2.5)+
#   
#   scale_color_viridis_d(name="Data Group")+
#   geom_vline(xintercept = 1000)+
#   scale_x_continuous(trans=ssqrt_trans)+
#   theme(text=element_text(face="bold",size=20,colour="black")) +
#   labs(y="Relative proportion (scaled to max bin height)", x="NetMHCIIpan BA (nM)")
# 
# 
# ####MHCII by rank
# 
# iedb_II_rank = iedb_melt[which(iedb_melt$Program %in% c('NetMHCIIpan_Rank'))]
# 
# Netii=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_all_BA.txt")
# 
# 
# iedb = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/mhc_sars2_purified_competitive_radioactivity.csv")
# iedb$allele = str_replace_all(iedb$allele, pattern = "\\*", replacement = "")
# iedb$allele = str_replace_all(iedb$allele, pattern = "HLA-DRB101:01", replacement = "DRB1_0101")
# iedb$affinity = as.numeric(iedb$affinity)
# 
# iedb_unique = paste(iedb$peptide, iedb$allele, sep = "_")
# Netii_unique = paste(Netii$Peptide, Netii$Haplotype, sep = "_")
# 
# iedb$NetMHCIIpan_Rank = NA
# for(z in 1:length(iedb_unique)){
#   if(iedb_unique[z] %in% Netii_unique){
#     iedb$NetMHCIIpan_Rank[z] = Netii$Rank[which(Netii_unique == iedb_unique[z])]
#   }
# }
# iedb_comp = iedb[complete.cases(iedb),]
# 
# 
# ggplot(data=iedb_comp, aes(x=NetMHCIIpan_Rank, y=affinity))+
#   geom_point()+
#   scale_y_log10()+#breaks = c(0.001, 0.01, .1, 1, 10, 100, 1000, 10000))+
#   scale_x_log10()+
#   #geom_vline(xintercept = 1000)+
#   geom_hline(yintercept = 1000)+
#   theme(text=element_text(face="bold",size=20,colour="black")) +
#   labs(x="NetMHCIIpan Rank", y="IEDB affinity (nM)")
#   #annotate("text",x=10, y = .3, label = paste0("n=", length(which((iedb_comp$NetMHCIIpan_Rank<1000)&(iedb_comp$affinity<1000)))), color="red", size=10)+
#   #annotate("text",x=10, y = 10000, label = paste0("n=", length(which((iedb_comp$NetMHCIIpan_Rank<1000)&(iedb_comp$affinity>1000)))), color="red", size=10)
# 
# cut = 1000
# iedb_comp$Hits = ifelse(iedb_comp$affinity < cut,0,1)
# cutoff = seq(0,100,.1)
# 
# FPR = c()
# FDR = c()
# Specificity=c()
# Sensitivity = c()
# PPV = c()
# NPV = c()
# 
# for(n in cutoff){
#   print(n)
#   #FPR = c(FPR, length(which(iedb_comp$affinity[which(iedb_comp$NetMHCIIpan_Rank<=n)] > cut))/length(which(iedb_comp$affinity>cut)))
#   FDR = c(FDR, length(which(iedb_comp$affinity[which(iedb_comp$NetMHCIIpan_Rank<=n)] > cut))/length(which(iedb_comp$NetMHCIIpan_Rank<=n)))
#   Sensitivity = c(Sensitivity,  length(which(iedb_comp$affinity[which(iedb_comp$NetMHCIIpan_Rank<=n)] < cut))/length(which(iedb_comp$affinity<cut)) )
#   #Specificity = c(Specificity,  length(which(iedb_comp$affinity[which(iedb_comp$NetMHCIIpan_Rank>n)] > cut))/length(which(iedb_comp$affinity>cut)) )
#   #PPV = c(PPV,  length(which(iedb_comp$affinity[which(iedb_comp$NetMHCIIpan_Rank<=n)] < cut))/length(which(iedb_comp$NetMHCIIpan_Rank<=n)) )
#   #NPV = c(NPV,  length(which(iedb_comp$affinity[which(iedb_comp$NetMHCIIpan_Rank>n)] > cut))/length(which(iedb_comp$NetMHCIIpan_Rank>n)) )
#   
# }
# 
# plot_cutoff = data.table(cutoff, FDR, Sensitivity)
# ggplot(data=plot_cutoff, aes(x=cutoff, y=FDR))+
#   geom_point()+
#   #scale_x_continuous(limits = c(0,1000))+
#   geom_hline(yintercept = c(.01,.05,.1))+
#   geom_smooth(method='lm',formula= y~+log(x))
# 
# 
# library(precrec)
# precrec_obj <- evalmod(scores = iedb_comp$NetMHCIIpan_Rank, labels = iedb_comp$Hits)
# autoplot(precrec_obj)
# 
# 
# iedb_II_rank = iedb_melt[which(iedb_melt$Program %in% c('NetMHCIIpan_Rank'))]
# 
# 
# ggplot()+
#   geom_histogram(data = Netii, aes(x=Rank, y = ..ncount.., color="All netMHCIIpan calls"), alpha=0.3, binwidth = 2.5)+
#   geom_histogram(data=iedb_II_rank, aes(x= Rank,y = ..ncount.., color="NetMHCIIpan calls with IEDB support"), alpha = 0.3, binwidth = 2.5)+
#   geom_histogram(data = Netii[which(Netii$Haplotype =="DRB1_0101")], aes(x=Rank, y = ..ncount.., color="DRB1*0101 netMHCIIpan calls"), alpha=0.3, binwidth = 2.5)+
#   
#   scale_color_viridis_d(name="Data Group")+
#   #geom_vline(xintercept = 1000)+
#   scale_x_continuous()+
#   theme(text=element_text(face="bold",size=20,colour="black")) +
#   labs(y="Relative proportion (scaled to max bin height)", x="NetMHCIIpan Rank")
# 
# 

# 
# ###Reading in MixMHC2pred
# 
# mix = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/mhc_all_viruses_combined_DRB_only_res.txt")
# mix = unique(mix)
# colnames(mix)[seq(7,19,2)] = substr(colnames(mix)[seq(7,19,2)],7,16)
# mix = melt(mix, id.vars = 1, measure.vars = c(seq(7,19,2)))
# mix$variable = paste0("HLA-", sapply(strsplit(as.character(mix$variable),split = "_"),'[',1), "*", 
#                       sapply(strsplit(as.character(mix$variable),split = "_"),'[',2), ":",
#                       sapply(strsplit(as.character(mix$variable),split = "_"),'[',3))
# colnames(mix) = c("Peptide", "Haplotype", "Percentile")
# 
# 
# iedb_v$affinity = as.numeric(iedb_v$affinity)
# 
# iedb_unique = paste(iedb_v$peptide, iedb_v$allele, sep = "_")
# mix_unique = paste(mix$Peptide, mix$Haplotype, sep = "_")
# 
# iedb_v$Mix_Rank = NA
# for(z in 1:length(iedb_unique)){
#   if(iedb_unique[z] %in% mix_unique){
#     #print(paste0(z,": ", which(mix_unique == iedb_unique[z])))
#     iedb_v$Mix_Rank[z] = mix$Percentile[which(mix_unique == iedb_unique[z])]
#   }
# }
# iedb_comp = iedb_v[complete.cases(iedb_v),]
# 
# 
# ggplot(data=iedb_comp, aes(x=Mix_Rank, y=affinity, color = allele), alpha=0.1)+
#   geom_point()+
#   scale_y_log10()+#breaks = c(0.001, 0.01, .1, 1, 10, 100, 1000, 10000))+
#   #scale_x_log10()+
#   #geom_vline(xintercept = 1000)+
#   geom_hline(yintercept = 500)+
#   theme(text=element_text(face="bold",size=20,colour="black")) +
#   labs(x="MixMHC2pred Rank", y="IEDB affinity (nM)")
# #annotate("text",x=10, y = .3, label = paste0("n=", length(which((iedb_comp$Mix_Rank<1000)&(iedb_comp$affinity<1000)))), color="red", size=10)+
# #annotate("text",x=10, y = 10000, label = paste0("n=", length(which((iedb_comp$Mix_Rank<1000)&(iedb_comp$affinity>1000)))), color="red", size=10)
# 
# cut = 2000
# iedb_comp$Hits = ifelse(iedb_comp$affinity < cut,0,1)
# cutoff = seq(0,100,.1)
# 
# FPR = c()
# FDR = c()
# Specificity=c()
# Sensitivity = c()
# PPV = c()
# NPV = c()
# 
# for(n in cutoff){
#   print(n)
#   #FPR = c(FPR, length(which(iedb_comp$affinity[which(iedb_comp$Mix_Rank<=n)] > cut))/length(which(iedb_comp$affinity>cut)))
#   FDR = c(FDR, length(which(iedb_comp$affinity[which(iedb_comp$Mix_Rank<=n)] > cut))/length(which(iedb_comp$Mix_Rank<=n)))
#   Sensitivity = c(Sensitivity,  length(which(iedb_comp$affinity[which(iedb_comp$Mix_Rank<=n)] < cut))/length(which(iedb_comp$affinity<cut)) )
#   #Specificity = c(Specificity,  length(which(iedb_comp$affinity[which(iedb_comp$Mix_Rank>n)] > cut))/length(which(iedb_comp$affinity>cut)) )
#   #PPV = c(PPV,  length(which(iedb_comp$affinity[which(iedb_comp$Mix_Rank<=n)] < cut))/length(which(iedb_comp$Mix_Rank<=n)) )
#   #NPV = c(NPV,  length(which(iedb_comp$affinity[which(iedb_comp$Mix_Rank>n)] > cut))/length(which(iedb_comp$Mix_Rank>n)) )
#   
# }
# 
# plot_cutoff = data.table(cutoff, FDR, Sensitivity)
# ggplot(data=plot_cutoff, aes(x=cutoff, y=FDR))+
#   geom_point()+
#   #scale_x_continuous(limits = c(0,1000))+
#   geom_hline(yintercept = c(.01,.05,.1))+
#   geom_smooth(method='lm',formula= y~+log(x))
# 
# 
# library(precrec)
# precrec_obj <- evalmod(scores = iedb_comp$Mix_Rank, labels = iedb_comp$Hits)
# autoplot(precrec_obj)
# 
# 
# iedb_II_rank = iedb_melt[which(iedb_melt$Program %in% c('Mix_Rank'))]
# 
# 
# ggplot()+
#   geom_histogram(data = Netii, aes(x=Rank, y = ..ncount.., color="All netMHCIIpan calls"), alpha=0.3, binwidth = 2.5)+
#   geom_histogram(data=iedb_II_rank, aes(x= Rank,y = ..ncount.., color="NetMHCIIpan calls with IEDB support"), alpha = 0.3, binwidth = 2.5)+
#   geom_histogram(data = Netii[which(Netii$Haplotype =="DRB1_0101")], aes(x=Rank, y = ..ncount.., color="DRB1*0101 netMHCIIpan calls"), alpha=0.3, binwidth = 2.5)+
#   
#   scale_color_viridis_d(name="Data Group")+
#   #geom_vline(xintercept = 1000)+
#   scale_x_continuous()+
#   theme(text=element_text(face="bold",size=20,colour="black")) +
#   labs(y="Relative proportion (scaled to max bin height)", x="NetMHCIIpan Rank")
# 



##################Tetramer analysis#############

Tet = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/tcell_virus_tetramer.csv")
Tet$Length = nchar(Tet$peptide)
Tet = Tet[which(Tet$species =="human"),]

Tet1 = Tet[which(Tet$mhc_class == "I")]
Tet2 = Tet[which(Tet$mhc_class == "II")]

Tet1 = Tet1[which(Tet1$Length %in% c(8:14)),]

Tet2$allele[which(substr(Tet2$allele,1,7)=="HLA-DRA")] = sapply(strsplit(Tet2$allele[which(substr(Tet2$allele,1,7)=="HLA-DRA")], "/"), "[", 2)
Tet2$allele[which(substr(Tet2$allele,1,3)=="DRB")] = paste0("HLA-", Tet2$allele[which(substr(Tet2$allele,1,3)=="DRB")])
Tet2 = Tet2[which(substr(Tet2$allele,1,7) == "HLA-DRB")]

# Alleles_to_keep1=c("HLA-A*11:01", "HLA-A*02:01", "HLA-A*01:01", "HLA-A*03:01", "HLA-A*24:02", "HLA-B*44:03", "HLA-B*07:02", 
#                    "HLA-B*08:01", "HLA-B*44:02", "HLA-B*35:01", "HLA-B*15:01", "HLA-C*03:04", "HLA-C*04:01", "HLA-C*05:01", 
#                    "HLA-C*06:02", "HLA-C*07:01", "HLA-C*07:02")
# 
# 
# Alleles_to_keep2=c("HLA-DRB1*01:01", "HLA-DRB1*03:01", "HLA-DRB1*04:01", "HLA-DRB1*07:01", "HLA-DRB1*11:01", "HLA-DRB1*13:01", "HLA-DRB1*15:01",
#                    "HLA-DQA1*01:02/DQB1*06:02", "HLA-DQA1*05:01/DQB1*02:01", "HLA-DQA1*02:01/DQB1*02:02", "HLA-DQA1*05:05/DQB1*03:01",
#                    "HLA-DQA1*01:01/DQB1*05:01", "HLA-DQA1*03:01/DQB1*03:02", "HLA-DQA1*03:03/DQB1*03:01", "HLA-DQA1*01:03/DQB1*06:03")

#Tet1 = Tet1[which(Tet1$allele %in% Alleles_to_keep1),]
#Tet2 = Tet2[which(Tet2$allele %in% Alleles_to_keep2),]


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


###Read in netMHC(II)pan

##BA
MHCI_tet_8mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCI_tetramer_8mer.xls")
MHCI_tet_9mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCI_tetramer_9mer.xls")
MHCI_tet_10mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCI_tetramer_10mer.xls")
MHCI_tet_11mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCI_tetramer_11mer.xls")
MHCI_tet_12mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCI_tetramer_12mer.xls")
MHCI_tet_13mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCI_tetramer_13mer.xls")
MHCI_tet_14mer = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Tetramer_calls/MHCI_tetramer_14mer.xls")
MHCI_tet = rbind(MHCI_tet_8mer, MHCI_tet_9mer, MHCI_tet_10mer, MHCI_tet_11mer, MHCI_tet_12mer, MHCI_tet_13mer, MHCI_tet_14mer)

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


##################
MHCI_tet$Tetramer = NA
for(n in 1:nrow(MHCI_tet)){
  print(n)
  if (paste0(MHCI_tet$Haplotype[n], "_", MHCI_tet$Peptide[n]) %in% paste0(Tet1$allele, "_", Tet1$peptide)){
    MHCI_tet$Tetramer[n] = Tet1$fraction_entries_positive[which(paste0(Tet1$allele, "_", Tet1$peptide) == paste0(MHCI_tet$Haplotype[n], "_", MHCI_tet$Peptide[n]))]
  }
}

MHCI_tet = MHCI_tet[complete.cases(MHCI_tet),]
colnames(MHCI_tet)[6] = "BA_Score"
colnames(MHCI_tet)[9] = "EL_Score"

fwrite(MHCI_tet, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCI_tetramer_all_values.txt")


###GLM with 5-fold CV
MHCI_tet = fread( "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCI_tetramer_all_values.txt")
MHCI_tet$Tetramer = round(MHCI_tet$Tetramer)
MHCI_tet$Tetramer[which(MHCI_tet$Tetramer == 0)] = "False"
MHCI_tet$Tetramer[which(MHCI_tet$Tetramer == 1)] = "True"

MHCI_tet$Tetramer = factor(MHCI_tet$Tetramer)
MHCI_tet$Haplotype = factor(MHCI_tet$Haplotype)

###Try normalizing
MHCI_tet$Binding_affinity = MHCI_tet$Binding_affinity/50000
MHCI_tet$BA_Rank = MHCI_tet$BA_Rank/100
MHCI_tet$EL_Rank = MHCI_tet$EL_Rank/100
MHCI_tet$Flurry_BA = MHCI_tet$Flurry_BA/50000

set.seed(1)
folds <- createFolds(factor(MHCI_tet$Haplotype), k = 5, list = FALSE)


spec_cutpoint = 0.9
#mean_auc = c()
performance = c()
Score = c()
FPR_test = c()

for(z in 1:5){
  train = MHCI_tet[folds != z,]
  test = MHCI_tet[folds == z,]
  
  print(table(test$Haplotype))
  
  # BA_Rank + BA_Score + Binding_affinity + EL_Rank + EL_Score + Flurry_BA + Flurry_pres_score + Flurry_proc_score +  Aromatic + Acidic + Basic + Small + Cyclic + Thiol
  model = glm(Tetramer ~  BA_Score + Binding_affinity + EL_Score + Flurry_BA + Flurry_pres_score + Flurry_proc_score +
                          Aromatic + Acidic + Basic + Small + Cyclic + Thiol,   # model to fit
              data = train,
              family = binomial)
  model_step = step(object = model, direction = "backward", trace = F)
  print(model_step$formula)
  
  
  ROC=roc(predictor=model_step$fitted.values,train$Tetramer, quiet = T)
  auc = as.numeric(substr(ROC$auc, nchar(ROC$auc) - 21, nchar(ROC$auc)))
  
  performance = c(performance, auc) 
  
  Specificity = c()
  for (n in cutoff){
    Specificity = c(Specificity,  length(which(train$Tetramer[which(model_step$fitted.values<n)] == "False"))/length(which(train$Tetramer == "False")))
  }
 Score = c(Score, cutoff[min(which(Specificity>=spec_cutpoint))])
 
 pred = predict(model_step, newdata = test, type = "response")
 FPR_test = c(FPR_test, length(which((pred > cutoff[min(which(Specificity>=spec_cutpoint))]) & (test$Tetramer == 'False')))/length(which(test$Tetramer == "False")))
  
}

performance
FPR_test


###Now model on all

model_all  = glm(Tetramer ~  BA_Score + Binding_affinity + EL_Score + Flurry_BA + Flurry_pres_score + Flurry_proc_score+
                   Aromatic + Acidic + Basic + Small + Cyclic + Thiol, # model to fit
                 data = MHCI_tet,
                 family = binomial)
model_all = step(model_all, direction = "backward", trace = F)




cutoff = seq(0,1,.001)

FDR = c()
Sensitivity = c()
Specificity = c()
PPV = c()

for(n in cutoff){
  print(n)
  FDR = c(FDR, length(which(MHCI_tet$Tetramer[which(model_all$fitted.values>=n)] == "False"))/length(which(model_all$fitted.values>=n)))
  Sensitivity = c(Sensitivity,  length(which(MHCI_tet$Tetramer[which(model_all$fitted.values>=n)] == "True"))/length(which(MHCI_tet$Tetramer == "True")))
  Specificity = c(Specificity,  length(which(MHCI_tet$Tetramer[which(model_all$fitted.values<n)] == "False"))/length(which(MHCI_tet$Tetramer == "False")))
  PPV = c(PPV, length(which(MHCI_tet$Tetramer[which(model_all$fitted.values>n)] == "True"))/length(which(model_all$fitted.values>n)))
}

plot_cutoff = data.table(cutoff, FDR, Sensitivity, Specificity, PPV)
glm_cutpoint = cutoff[min(which(Specificity >= spec_cutpoint))]

perf1 = ggplot(data=plot_cutoff, aes(x=cutoff, y=Specificity))+
  geom_point()+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x="GLM score", y= "Specificity", title = "HLA-I IEDB model performance\nAll viral tetramers")+
  geom_hline(yintercept = c(spec_cutpoint))+
  geom_vline(xintercept = glm_cutpoint)+
  annotate("text", label =paste0(spec_cutpoint),y=spec_cutpoint, x=0, vjust=0, hjust=0, size = 10)+
  annotate("text", label =paste0(glm_cutpoint),y=0, x=glm_cutpoint, vjust=0, hjust=0, angle = 90, size = 10)

ggplot(data=plot_cutoff, aes(x=cutoff, y=PPV))+
  geom_point()+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x="GLM score", y= "PPV", title = "NetMHCpan v4.0 vs IEDB performance\nAll viral tetramers")


saveRDS(model_all, paste0("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCI_glm.R"))


####MHCII
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


###GLM with 5-fold CV
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

set.seed(1)
folds <- createFolds(factor(MHCII_tet$Haplotype), k = 5, list = FALSE)

spec_cutpoint = 0.9
#mean_auc = c()
performance = c()
Score = c()
FPR_test = c()

for(z in 1:5){
  train = MHCII_tet[folds != z,]
  test = MHCII_tet[folds == z,]
  
  #model_weights <- ifelse(train$Tetramer == "False",
  #                        table(train$Tetramer)[2],
  #                        table(train$Tetramer)[1])
  
  print(table(test$Haplotype))
  # BA_Rank + BA_Score + Binding_affinity + EL_Rank + EL_Score + Aromatic + Acidic + Basic + Small + Cyclic + Thiol
  model = glm(Tetramer ~  BA_Score + Binding_affinity + EL_Score +
                Aromatic + Acidic + Basic + Small + Cyclic + Thiol,   # model to fit
              data = train,
              family = binomial)
  model_step = step(object = model, direction = "backward", trace = F)
  print(model_step$formula)
  
  ROC=roc(predictor=model_step$fitted.values,train$Tetramer, quiet = T)
  auc = as.numeric(substr(ROC$auc, nchar(ROC$auc) - 21, nchar(ROC$auc)))
  
  performance = c(performance, auc) 
  
  Specificity = c()
  for (n in cutoff){
    Specificity = c(Specificity,  length(which(train$Tetramer[which(model_step$fitted.values<n)] == "False"))/length(which(train$Tetramer == "False")))
  }
  Score = c(Score, cutoff[min(which(Specificity>=spec_cutpoint))])
  
  pred = predict(model_step, newdata = test, type = "response")
  FPR_test = c(FPR_test, length(which((pred > cutoff[min(which(Specificity>=spec_cutpoint))]) & (test$Tetramer == 'False')))/length(which(test$Tetramer == "False")))
  
}

performance
FPR_test


######Generating full model
#model_weights <- ifelse(MHCII_tet$Tetramer == "False",
#                       table(MHCII_tet$Tetramer)[2],
#                        table(MHCII_tet$Tetramer)[1])

model_all  = glm(Tetramer ~ Binding_affinity +BA_Score + Binding_affinity  + EL_Score +
                   Aromatic + Acidic + Basic + Small + Cyclic + Thiol,   # model to fit
                 data = MHCII_tet, 
                 family = binomial)

model_all = step(model_all, direction = "backward", trace = F)

cutoff = seq(0,1,.001)

FDR = c()
Sensitivity = c()
Specificity = c()
PPV=c()
for(n in cutoff){
  print(n)
  FDR = c(FDR, length(which(MHCII_tet$Tetramer[which(model_all$fitted.values>=n)] == "False"))/length(which(model_all$fitted.values>=n)))
  Sensitivity = c(Sensitivity,  length(which(MHCII_tet$Tetramer[which(model_all$fitted.values>=n)] == "True"))/length(which(MHCII_tet$Tetramer == "True")))
  Specificity = c(Specificity,  length(which(MHCII_tet$Tetramer[which(model_all$fitted.values<n)] == "False"))/length(which(MHCII_tet$Tetramer == "False")))
  PPV = c(PPV, length(which(MHCI_tet$Tetramer[which(model_all$fitted.values>n)] == "True"))/length(which(model_all$fitted.values>n)))
  
}

plot_cutoff = data.table(cutoff, FDR, Sensitivity, Specificity,PPV)
glm_cutpoint = cutoff[min(which(Specificity >= spec_cutpoint))]

perf2 = ggplot(data=plot_cutoff, aes(x=cutoff, y=Specificity))+
  geom_point()+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x="GLM score", y= "Specificity", title = "HLA-II IEDB model performance\nAll viral tetramers")+
  geom_hline(yintercept = c(spec_cutpoint))+
  geom_vline(xintercept = glm_cutpoint)+
  annotate("text", label =paste0(spec_cutpoint),y=spec_cutpoint, x=0, vjust=0, hjust=0, size = 10)+
  annotate("text", label =paste0(glm_cutpoint),y=0, x=glm_cutpoint, vjust=0, hjust=0, angle = 90, size = 10)
  #geom_smooth(method = "lm", formula = y~exp(x))


ggplot(data=plot_cutoff, aes(x=cutoff, y=PPV))+
  geom_point()+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x="GLM score", y= "PPV", title = "IEDB model performance\nAll viral tetramers")

saveRDS(model_all, paste0("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCII_glm.R"))

###Plot cutpoints

grid.arrange(perf1, perf2, nrow=1, ncol=2)


###Add in GLM scores from tetramer analysis
MHCI_glm = readRDS("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCI_glm.R")
MHCII_glm = readRDS("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCII_glm.R")


Net_BA =  fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_BA.txt")
Net_EL= fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_all_EL.txt")
Flurry = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCflurry_standardized_all.txt")

Net_BA = Net_BA[order(Net_BA$Haplotype),]
Net_BA = Net_BA[order(Net_BA$Peptide),]
Net_EL = Net_EL[order(Net_EL$Haplotype),]
Net_EL = Net_EL[order(Net_EL$Peptide),]
Flurry = Flurry[order(Flurry$Haplotype),]
Flurry = Flurry[order(Flurry$Peptide),]

Net_combined = cbind(Net_BA, Net_EL[,c(5:6)], Flurry[,c(5,7:9)])
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

feature_tab = data.table()

for(f in features){
  print(f)
  temp_list = data.table()
  for(n in 1:length(f)){
    print(n)
    temp_list = cbind(temp_list, c(str_count(Net_combined$Peptide, pattern = f[n])))
  }
  hits = rowSums(temp_list)/nchar(Net_combined$Peptide)
  feature_tab = cbind(feature_tab, hits)
}
colnames(feature_tab) = c("Aromatic", "Acidic", "Basic", "Small", "Cyclic", "Thiol")

Net_combined = cbind(Net_combined, feature_tab)

#####Run model, write out

Net_combined$GLM_score = predict(MHCI_glm, newdata = Net_combined, type = "response")
Net_combined$Predicted_tetramer = ifelse(Net_combined$GLM_score >= .697, "True", "False")

fwrite(Net_combined, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAI_all_glm.txt")

##########MHCII############
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
Netii_all$Predicted_tetramer = ifelse(Netii_all$GLM_score >= .768, "True", "False")

fwrite(Netii_all, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLAII_all_glm.txt")



###################Check distributions##################


plot = data.table(MHCI_tet$Tetramer, MHCI_glm$fitted.values)
h1=ggplot(data=plot)+
  geom_histogram(aes(x = V2, fill = V1), alpha=.5,bins = 50, position="dodge")+
  geom_vline(xintercept = .697)+
  annotate("text", label = ".697", x=.697, y=20, vjust=0, hjust=0, angle = 90, size = 10)+
  scale_fill_manual(values = c("blue","red"), name = "Tetramer\npositive")+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x= "GLM score", y="Count", title = "IEDB HLA-I")


plot2 = data.table(MHCII_tet$Tetramer, MHCII_glm$fitted.values)
h2=ggplot(data=plot2)+
  geom_histogram(aes(x = V2, fill = V1), alpha=.5,bins = 50, position="dodge")+
  geom_vline(xintercept = .768)+
  annotate("text", label = ".768", x=.768, y=250, vjust=0, hjust=0, angle = 90, size = 10)+
  scale_fill_manual(values = c("blue","red"), name = "Tetramer\npositive")+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x= "GLM score", y="Count", title = "IEDB HLA-II")

grid.arrange(h1, h2, nrow=1, ncol=2)


SARS_comb = cbind(data.table(c(rep("HLA-I", nrow(Net_combined)), rep("HLA-II", nrow(Netii_all)))),
                  data.table(c(Net_combined$GLM_score, Netii_all$GLM_score)))
colnames(SARS_comb) = c("HLA", "GLM_score")
ggplot(data=SARS_comb)+
  geom_histogram(aes(x = GLM_score, fill = HLA), alpha=.5,bins = 50, position = "dodge")+
  geom_vline(xintercept = .697, color = "blue")+
  annotate("text", label= paste0(100*round(length(which(Net_combined$Predicted_tetramer=="True"))/nrow(Net_combined), digits = 3), "% positive"),
           x=.697, y=35000, vjust=0, hjust=0, angle = 90, size = 10, color = "blue")+
  geom_vline(xintercept = .768, color = "red")+
  annotate("text", label = paste0(100*round(length(which(Netii_all$Predicted_tetramer=="True"))/nrow(Netii_all), digits = 3), "% positive"),
                                  x=.768, y=35000, vjust=0, hjust=0, angle = 90, size = 10, color = "red")+
  scale_fill_manual(values = c("blue","red"), name = "")+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x= "GLM score", y="Count", title = "SARS-CoV-2 predicted epitopes")


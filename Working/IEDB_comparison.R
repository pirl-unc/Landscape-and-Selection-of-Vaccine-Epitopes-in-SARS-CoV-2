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
colnames(Net)[seq(8,87,5)] = colnames(Net)[seq(7,86,5)]
colnames(Net)[seq(7,86,5)] = paste0("nM_", colnames(Net)[seq(7,86,5)])

Net_rank = melt(Net, id.vars = 2, measure.vars = c(seq(8,87,5)))
Net_nM = melt(Net, id.vars = 2, measure.vars = c(seq(7,86,5)))
Net_log = melt(Net, id.vars = 2, measure.vars = c(seq(6,85,5)))
Net_Protein = melt(Net, id.vars = 3, measure.vars = c(seq(8,87,5)))
Net_Start = melt(Net, id.vars = 1, measure.vars = c(seq(8,87,5)))

Net = cbind(Net_Protein$ID, Net_Start$Pos, Net_rank, Net_log[,3], Net_nM[,3])
colnames(Net) = c("Protein", "Start", "Peptide", "Haplotype", "Rank", "1-log50k", "Binding_affinity")
Net$Protein = sapply(strsplit(Net$Protein, "_"),'[',1)
Net$Protein[which(Net$Protein =="orflab")] = "orf1ab"
Net$Start = Net$Start+1

Net$Haplotype = paste0(substr(Net$Haplotype,1,5), "*", substr(Net$Haplotype,6,10))

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
colnames(Net)[seq(8,87,5)] = colnames(Net)[seq(7,86,5)]
colnames(Net)[seq(7,86,5)] = paste0("nM_", colnames(Net)[seq(7,86,5)])

Net_rank = melt(Net, id.vars = 2, measure.vars = c(seq(8,87,5)))
Net_nM = melt(Net, id.vars = 2, measure.vars = c(seq(7,86,5)))
Net_log = melt(Net, id.vars = 2, measure.vars = c(seq(6,85,5)))
Net_Protein = melt(Net, id.vars = 3, measure.vars = c(seq(8,87,5)))
Net_Start = melt(Net, id.vars = 1, measure.vars = c(seq(8,87,5)))

Net = cbind(Net_Protein$ID, Net_Start$Pos, Net_rank, Net_log[,3],  Net_nM[,3])
colnames(Net) = c("Protein", "Start", "Peptide", "Haplotype", "Rank", "1-log50k", "Binding_affinity")
Net$Protein = sapply(strsplit(Net$Protein, "_"),'[',1)
Net$Protein[which(Net$Protein =="orflab")] = "orf1ab"
Net$Start = Net$Start+1

Net$Haplotype = paste0(substr(Net$Haplotype,1,5), "*", substr(Net$Haplotype,6,10))


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



############NetMHCIIpan -- SARS##################

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

######Compared performance vs IEDB##############
iedb = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/mhc_sars2_purified_competitive_radioactivity.csv")
iedb = iedb[,c(1:3)]

iedb_low = iedb[which(iedb$affinity < 500),]
iedb_low = iedb_low[order(iedb_low$peptide),]


Flurry=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCflurry_standardized_all.txt")
Net_EL=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_all_EL.txt")
Net_BA=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_all_BA.txt")

#Netii = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_all_BA.txt")

Net_EL = Net_EL[-which(Net_EL$Haplotype =="HLA-B*15:01" ),]
Net_BA = Net_BA[-which(Net_BA$Haplotype =="HLA-B*15:01" ),]

iedb_unique = paste(iedb_low$peptide, iedb_low$allele, sep = "_")
Net_BA_unique = paste(Net_BA$Peptide, Net_BA$Haplotype, sep = "_")
Net_EL_unique = paste(Net_EL$Peptide, Net_EL$Haplotype, sep = "_")
Flurry_unique = paste(Flurry$Peptide, Flurry$Haplotype, sep = "_")
Netii_unique = paste(Netii$Peptide, Netii$Haplotype, sep = "_")

Net_BA = Net_BA[which(Net_BA_unique %in% iedb_unique),]
Net_EL = Net_EL[which(Net_EL_unique %in% iedb_unique),]
Flurry = Flurry[which(Flurry_unique %in% iedb_unique),]
#Netii = Netii[which(Netii_unique %in% iedb_unique),]

all = paste(Net_BA$Peptide, Net_BA$Haplotype, sep = "_")#, paste(Netii$Peptide, Netii$Haplotype, sep = "_"))
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

iedb_low = iedb[which(iedb$affinity < 500),]
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
Net_BA = Net_BA[-which(Net_BA$Haplotype =="HLA-B*15:01" ),]


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


scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

library(RColorBrewer)

p=ggplot(data=iedb_comp, aes(x=NetMHCpan_BA_BA, y=affinity, color = factor(allele)))+
  geom_point()+
  scale_color_manual(values = c(brewer.pal(7, 'Blues')[3:7],
                                brewer.pal(7, 'Reds')[3:7],
                                brewer.pal(7, 'Greens')[3:7]),
                     name="Allele", drop = F, guide = guide_legend(override.aes = list(color = "white")))+
  scale_y_log10(breaks = c(0.001, 0.01, .1, 1, 10, 100, 1000, 10000,100000), label=scientific_10)+
  scale_x_continuous(trans= 'log10', breaks = c(1,10,100,1000,10000,100000,1000000), label=scientific_10, limits= c(1,100000))+
  geom_vline(xintercept = 300)+
  geom_hline(yintercept = 500)+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x="NetMHCpan BA (nM)", y="IEDB affinity (nM)")+
  annotate("text",x=5, y = .3, label = paste0("n=", length(which((iedb_comp$NetMHCpan_BA_BA<300)&(iedb_comp$affinity<500)))), color="red", size=10, hjust=0)+
  annotate("text",x=5, y = 10000, label = paste0("n=", length(which((iedb_comp$NetMHCpan_BA_BA<300)&(iedb_comp$affinity>500)))), color="red", size=10,hjust=0)+
  annotate("text",x=10000, y = .3, label = paste0("n=", length(which((iedb_comp$NetMHCpan_BA_BA>300)&(iedb_comp$affinity<500)))), color="red", size=10,hjust=0)+
  annotate("text",x=10000, y = 10000, label = paste0("n=", length(which((iedb_comp$NetMHCpan_BA_BA>300)&(iedb_comp$affinity>500)))), color="red", size=10,hjust=0)+
  
  annotate("text",x=2.5, y = 350, label = "500nM", color="black", size=5,hjust=0)+
  annotate("text",x=350, y = .1, label = "300nM", color="black", size=5, angle=90, hjust=0, vjust=0)+
  theme(
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")
  ) 


hist_top <- ggplot(data=Net_BA)+
  geom_histogram(aes(x= Binding_affinity, fill = Haplotype), binwidth = .01)+
  scale_fill_manual(values = c(brewer.pal(7, 'Blues')[3:7],
                               brewer.pal(7, 'Reds')[3:7],
                               brewer.pal(7, 'Greens')[3:7]),
                    name="Allele", drop = F)+
  scale_x_continuous(trans= 'log10', breaks = c(1,10,100,1000,10000,100000,1000000), limits= c(1,100000))+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  annotate("text",x=5, y = 300000, label = paste0("n=", length(which(Net_BA$Binding_affinity<300))), color="red", size=10,hjust=0)+
  annotate("text",x=10000, y = 300000, label = paste0("n=", length(which(Net_BA$Binding_affinity>300))), color="red", size=10,hjust=0)+
  geom_vline(xintercept = 300)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  labs(y="Count")+
  scale_y_continuous(trans = sqrt_trans(), label=scientific_10)

grid.arrange(hist_top, p, ncol=1, nrow=2, heights = c(3,4))

#####Figure S1B

iedb_comp$Hits = ifelse(iedb_comp$affinity < 500,0,1)
iedb_comp = iedb[complete.cases(iedb),]
cutoff = seq(0,500,1)

FPR = c()
FDR = c()
Specificity=c()
Sensitivity = c()
PPV = c()
NPV = c()

for(n in cutoff){
  print(n)
  FPR = c(FPR, length(which(iedb_comp$affinity[which(iedb_comp$NetMHCpan_BA_BA<=n)] > 500))/length(which(iedb_comp$affinity>500)))
  FDR = c(FDR, length(which(iedb_comp$affinity[which(iedb_comp$NetMHCpan_BA_BA<=n)] > 500))/length(which(iedb_comp$NetMHCpan_BA_BA<=n)))
  Sensitivity = c(Sensitivity,  length(which(iedb_comp$affinity[which(iedb_comp$NetMHCpan_BA_BA<=n)] < 500))/length(which(iedb_comp$affinity<500)) )
  Specificity = c(Specificity,  length(which(iedb_comp$affinity[which(iedb_comp$NetMHCpan_BA_BA>n)] > 500))/length(which(iedb_comp$affinity>500)) )
  PPV = c(PPV,  length(which(iedb_comp$affinity[which(iedb_comp$NetMHCpan_BA_BA<=n)] < 500))/length(which(iedb_comp$NetMHCpan_BA_BA<=n)) )
  NPV = c(NPV,  length(which(iedb_comp$affinity[which(iedb_comp$NetMHCpan_BA_BA>n)] > 500))/length(which(iedb_comp$NetMHCpan_BA_BA>n)) )
  
}

plot_cutoff = data.table(cutoff, FDR, Sensitivity)
ggplot(data=plot_cutoff, aes(x=cutoff, y=FDR))+
  geom_point()+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x="Binding affinity (nM)", y= "False discovery rate", title = "NetMHCpan v4.0 vs IEDB performance\nSARS -- 500nM IEDB cutoff")+
  #scale_x_continuous(limits = c(0,1000))+
  geom_hline(yintercept = c(.01,.05,.1))+
  geom_smooth(method='lm',formula= y~+x)


library(precrec)
precrec_obj <- evalmod(scores = iedb_comp$NetMHCpan_BA_BA, labels = iedb_comp$Hits)
autoplot(precrec_obj)

ggplot()+
  geom_histogram(data = Net_BA, aes(x=Binding_affinity, y = ..ncount.., color="All netMHCpan calls"), alpha=0.3, binwidth = 5)+
  geom_histogram(data=iedb_BA, aes(x= Rank,y = ..ncount.., color="NetMHCpan calls with IEDB support"), alpha = 0.3, binwidth = 5)+
  scale_color_viridis_d(name="Data Group")+
  geom_vline(xintercept = 300)+
  scale_x_continuous(trans=ssqrt_trans)+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(y="Relative proportion (scaled to max bin height)", x="NetMHCpan BA (nM)")


Net_BA_filt = Net_BA[which(Net_BA$Binding_affinity < 300),]

#fwrite(Net_BA,"/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_BA.txt")
#fwrite(Net_BA_filt, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_filtered_BA.txt")


################################################
##############Figure 2B#########################

Netii = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_mhc_all_viruses_all_lengths.txt")

iedb_v = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/mhc_all_virus_combined.csv")
iedb_v = iedb_v[which(iedb_v$allele%in% Alleles_to_keep),]
iedb_v = iedb_v[which(nchar(iedb_v$peptide) %in% c(12:20)),]

iedb_v$affinity = as.numeric(iedb_v$affinity)

iedb_unique = paste(iedb_v$peptide, iedb_v$allele, sep = "_")
Netii_unique = paste(Netii$Peptide, Netii$Haplotype, sep = "_")

#Netii = Netii[-which(duplicated(Netii_unique)),]
#Netii_unique = Netii_unique[-which(duplicated(Netii_unique))]

iedb_v$NetMHCIIpan_BA = NA
for(z in 1:length(iedb_unique)){
  if(iedb_unique[z] %in% Netii_unique){
    iedb_v$NetMHCIIpan_BA[z] = Netii$Binding_affinity[which(Netii_unique == iedb_unique[z])]
  }
}
iedb_comp = iedb_v[complete.cases(iedb_v),]

# 
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

cut = 500
iedb_comp$Hits = ifelse(iedb_comp$affinity < cut,0,1)

cutoff = seq(0,cut,1)

#FPR = c()
FDR = c()
#Specificity=c()
Sensitivity = c()
#PPV = c()
#NPV = c()

for(n in cutoff){
  print(n)
  #FPR = c(FPR, length(which(iedb_comp$affinity[which(iedb_comp$NetMHCIIpan_BA<=n)] > cut))/length(which(iedb_comp$affinity>cut)))
  FDR = c(FDR, length(which(iedb_comp$affinity[which(iedb_comp$NetMHCIIpan_BA<=n)] > cut))/length(which(iedb_comp$NetMHCIIpan_BA<=n)))
  Sensitivity = c(Sensitivity,  length(which(iedb_comp$affinity[which(iedb_comp$NetMHCIIpan_BA<=n)] < cut))/length(which(iedb_comp$affinity<cut)) )
  #Specificity = c(Specificity,  length(which(iedb_comp$affinity[which(iedb_comp$NetMHCIIpan_BA>n)] > cut))/length(which(iedb_comp$affinity>cut)) )
  #PPV = c(PPV,  length(which(iedb_comp$affinity[which(iedb_comp$NetMHCIIpan_BA<=n)] < cut))/length(which(iedb_comp$NetMHCIIpan_BA<=n)) )
  #NPV = c(NPV,  length(which(iedb_comp$affinity[which(iedb_comp$NetMHCIIpan_BA>n)] > cut))/length(which(iedb_comp$NetMHCIIpan_BA>n)) )
  
}

plot_cutoff = data.table(cutoff, FDR, Sensitivity)
ggplot(data=plot_cutoff, aes(x=cutoff, y=FDR))+
  geom_point()+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x="Binding affinity (nM)", y= "False discovery rate", title = "NetMHCIIpan v3.2 vs IEDB performance\nAll viruses -- 500nM IEDB cutoff")+
  #scale_x_continuous(limits = c(0,1000))+
  geom_hline(yintercept = c(.01,.05,.1))
# geom_smooth(method='lm',formula= y~x)

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
  geom_vline(xintercept = 245)+
  geom_hline(yintercept = 500)+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  labs(x="NetMHCpan BA (nM)", y="IEDB affinity (nM)")+
  annotate("text",x=5, y = .3, label = paste0("n=", length(which((iedb_comp$NetMHCIIpan_BA<245)&(iedb_comp$affinity<500)))), color="red", size=10, hjust=0)+
  annotate("text",x=5, y = 100000, label = paste0("n=", length(which((iedb_comp$NetMHCIIpan_BA<245)&(iedb_comp$affinity>500)))), color="red", size=10,hjust=0)+
  annotate("text",x=10000, y = .3, label = paste0("n=", length(which((iedb_comp$NetMHCIIpan_BA>245)&(iedb_comp$affinity<500)))), color="red", size=10,hjust=0)+
  annotate("text",x=10000, y = 100000, label = paste0("n=", length(which((iedb_comp$NetMHCIIpan_BA>245)&(iedb_comp$affinity>500)))), color="red", size=10,hjust=0)+
  
  annotate("text",x=5, y = 350, label = "500nM", color="black", size=5,hjust=0)+
  annotate("text",x=315, y = .1, label = "245nM", color="black", size=5, angle=90, hjust=0, vjust=0)+
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
  geom_vline(xintercept = 245)+
  annotate("text",x=5, y = 800, label = paste0("n=", length(which(Netii_SARS$Binding_affinity<245))), color="red", size=10,hjust=0)+
  annotate("text",x=10000, y = 800, label = paste0("n=", length(which(Netii_SARS$Binding_affinity>245))), color="red", size=10,hjust=0)+
  labs(y="Count")+
  scale_y_continuous(labels = scientific_10, breaks = c(0,200,400,600,800,10000))

grid.arrange(hist_top2, p2, ncol=1, nrow=2, heights = c(3,4))



Netii_SARS_filt = Netii_SARS[which(Netii_SARS$Binding_affinity < 245),]

#fwrite(Netii_SARS, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_BA.txt")
#fwrite(Netii_SARS_filt, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_filtered_BA.txt")






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
  Binding_haplotype=paste0(colnames(Net_BA_cast)[which(Net_BA_cast[n,] < 300)], collapse = ",")
  Frequency=1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% colnames(Net_BA_cast)[which(Net_BA_cast[n,] < 300)] )]/100,2)))
  data.frame(Binding_haplotype,Frequency)
}

Pop_freq_II = foreach(n=1:nrow(Netii_BA_cast), .combine=rbind) %dopar% {
  Binding_haplotype=paste0(colnames(Netii_BA_cast)[which(Netii_BA_cast[n,] < 245)], collapse = ",")
  Frequency=1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% colnames(Netii_BA_cast)[which(Netii_BA_cast[n,] < 245)] )]/100,2)))
  data.frame(Binding_haplotype,Frequency)
}

Net_BA_cast = cbind(Net_BA_cast, Pop_freq_I)
Netii_BA_cast = cbind(Netii_BA_cast, Pop_freq_II)

####Calculating co-epitopes

Best_II_hit = data.table()

for(n in 1:nrow(Net_BA_cast)){
  print(n)
  pep1 = Net_BA_cast$Peptide[n]
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
    Best_II_hit = rbind(Best_II_hit, data.table(Net_BA_cast$Protein[n],pep2_start, pep1, Net_BA_cast$Binding_haplotype[n], Net_BA_cast$Frequency[n], pep2, pep2_haps, pep2_freq))
  }else{
    Best_II_hit = rbind(Best_II_hit, data.table(Net_BA_cast$Protein[n],pep2_start, pep1, Net_BA_cast$Binding_haplotype[n], Net_BA_cast$Frequency[n], NA, NA, NA), use.names=F)
  }
}

Net_coepitopes = cbind(Best_II_hit, (Best_II_hit$V5*Best_II_hit$pep2_freq))
colnames(Net_coepitopes) = c("Protein", "Start", "c1_peptide", "c1_haplotypes", "c1_freq", "c2_peptide", "c2_haplotypes", "c2_freq", "Total_freq")
Net_coepitopes = Net_coepitopes[which(complete.cases(Net_coepitopes)),]

#fwrite(Net_coepitopes, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/Coepitopes_human.txt")




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

Net_coepitopes$Start[which(Net_coepitopes$Protein == "S")] = Net_coepitopes$Start[which(Net_coepitopes$Protein == "S")]+7096
Net_coepitopes$Start[which(Net_coepitopes$Protein == "ORF3a")] = Net_coepitopes$Start[which(Net_coepitopes$Protein == "ORF3a")]+8369
Net_coepitopes$Start[which(Net_coepitopes$Protein == "E")] = Net_coepitopes$Start[which(Net_coepitopes$Protein == "E")]+8644
Net_coepitopes$Start[which(Net_coepitopes$Protein == "M")] = Net_coepitopes$Start[which(Net_coepitopes$Protein == "M")]+8719
Net_coepitopes$Start[which(Net_coepitopes$Protein == "ORF6")] = Net_coepitopes$Start[which(Net_coepitopes$Protein == "ORF6")]+8941
Net_coepitopes$Start[which(Net_coepitopes$Protein == "ORF7a")] = Net_coepitopes$Start[which(Net_coepitopes$Protein == "ORF7a")]+9002
Net_coepitopes$Start[which(Net_coepitopes$Protein == "ORF8")] = Net_coepitopes$Start[which(Net_coepitopes$Protein == "ORF8")]+9123
Net_coepitopes$Start[which(Net_coepitopes$Protein == "N")] = Net_coepitopes$Start[which(Net_coepitopes$Protein == "N")]+9244
Net_coepitopes$Start[which(Net_coepitopes$Protein == "ORF10")] = Net_coepitopes$Start[which(Net_coepitopes$Protein == "ORF10")]+9663


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

y_ll=0
y_ul=55



p=ggplot(data=Net_coepitopes[which(Net_coepitopes$Total_freq>Min_frequency)]) + 
  geom_segment(data=Net_coepitopes[which(Net_coepitopes$Total_freq>Min_frequency)], aes(x=Start, xend=Start, y=y_ll, yend=Total_freq),
               color=ifelse(Net_coepitopes[which(Net_coepitopes$Total_freq>Min_frequency)]$Lit_overlap==1,'red','grey'), alpha=0.5, size=1.1) +
  geom_point(aes(x=Start, y=Total_freq, color=murine), alpha=.75, size = 5) +
  
  scale_color_manual(values = c("tomato", "cadetblue", "purple", "grey"), name="Murine\noverlap")+
  scale_size(range = c(0,10), breaks = c(20,40,60,80,90))+
  scale_y_continuous(limits = c((y_ll-13),y_ul), breaks = c(20,30,40,50,60))+
  
  geom_label_repel(aes(label=ifelse((Total_freq>label_frequency),paste0(as.character(c2_peptide),'\n',as.character(c1_peptide)),''), x=Start, y=Total_freq),
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
text(0.825,.5,"All human\nHLA-I epitopes\n3163")
text(0.425,.5,"Human\nHLA-I co-epitopes\n2749")
#text(0.535,.5,"Co-epitopes: 6289\nHLA-I: 2721\nHLA-II: 2721")



vd_II <- venneuler(c(HLAII=nrow(Netii_BA_cast), II_coep=0,"HLAII&II_coep"=length(unique(Net_coepitopes$c2_peptide))))
vd_II$labels=NA
#vd_II$centers = vd_I$centers
#vd$colors 

plot(vd_II, col= c(alpha(colour = "cyan", .5), alpha(colour = "blue",.5)))
text(0.825,.5,"All human\nHLA-II epitopes\n5296")
text(0.425,.5,"Human\nHLA-II co-epitopes\n2399")
#text(0.535,.5,"Co-epitopes: 6289\nHLA-I: 2721\nHLA-II: 2337")

vd_CO <- venneuler(c(Coep=nrow(Net_coepitopes), II_coep=100,"Coep&II_coep"=0))
vd_CO$labels=NA
#vd_II$centers = vd_I$centers
#vd$colors 

plot(vd_CO, col= c(alpha(colour = "purple", .5), alpha(colour = "blue",.5)))
text(0.75,.45,"Co-epitopes: 6443")

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
text(0.55,.35,paste0("Human\nCo-epitopes\n", (nrow(Net_coepitopes)-murine_1_coep-murine_2_coep)))
text(0.3,.575,paste0("Murine\nMHC-I\n", (nrow(mouse1)-murine_1_coep)))
text(0.675,.65,paste0("Murine\nMHC-II\n", (nrow(mouse2)-murine_2_coep)))
text(0.4125,.5,murine_1_coep)
text(0.625,.525,murine_2_coep)



#####################################################
########Figure 2F####################################
#Protein dist, normalized by length
all_prots=c("orf1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF8","N","ORF10")
prot_len= c(abs(1 - 7096),abs(7097 - 8369),abs(8370 - 8644),abs(8645 - 8719),abs(8720 - 8941),abs(8942 - 9002),abs(9003 - 9123), abs(9124 - 9244),abs(9245 - 9663),abs(9664 - 9701))

Net_coepitopes = Net_coepitopes[order(Net_coepitopes$Start),]
prots = factor(c(unique(Net_coepitopes$Protein)), levels = all_prots)
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

ggplot(prot_dist, aes(x="", y=Count, fill=Protein)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+
  scale_fill_manual(values=c(viridis(10)))+
  geom_label_repel(aes(y = ypos, label = Protein), color = "white", size=8, nudge_x = .65,
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'black')+
  theme(legend.position="none")+
  theme_void() 
#geom_label_repel(aes(y = Count, label = Protein),size=3, color="white", segment.color = "black")
#scale_fill_brewer(palette="Set1")

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
  scale_x_continuous(breaks = seq(0,100,10))+
  
  
  annotate("text", x=common_HLAI[1], y = 100, label = HLAI_lab1, color = "black", hjust=0, vjust=0, angle = 90)+
  annotate("text", x=common_HLAI[2], y = 100, label = HLAI_lab2, color = "black", hjust=0, vjust=0, angle = 90)+
  annotate("text", x=common_HLAI[3], y = 100, label = HLAI_lab3, color = "black", hjust=0, vjust=0, angle = 90)+
  annotate("text", x=common_HLAI[4], y = 100, label = HLAI_lab4, color = "black", hjust=0, vjust=0, angle = 90)+
  
  annotate("text", x=common_HLAII[1], y = 100, label = HLAII_lab1, color = "black", hjust=0, vjust=0, angle = 90)+
  annotate("text", x=common_HLAII[2], y = 100, label = HLAII_lab2, color = "black", hjust=0, vjust=0, angle = 90)+
  annotate("text", x=common_HLAII[3], y = 100, label = paste0(HLAII_lab3[1],"\n",HLAII_lab3[2]), color = "black", hjust=0, vjust=0, angle = 90)+
  annotate("text", x=common_HLAII[4], y = 100, label = HLAII_lab4, color = "black", hjust=0, vjust=0, angle = 90)






##########################################
#########Figure 3A -- IEDB epitopes########

####Left: Predicted epitopes


#Protein dist, normalized by length
all_prots=c("orf1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF8","N","ORF10")
prot_len= c(abs(1 - 7096),abs(7097 - 8369),abs(8370 - 8644),abs(8645 - 8719),abs(8720 - 8941),abs(8942 - 9002),abs(9003 - 9123), abs(9124 - 9244),abs(9245 - 9663),abs(9664 - 9701))


Net_coepitopes = Net_coepitopes[order(Net_coepitopes$Start),]
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
iedb_fig3 = iedb_fig3[which(iedb_fig3$affinity < 500),]
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

iedb_fig3_I = iedb_fig3[which(iedb_fig3$allele !=  "HLA-DRB1*01:01"),]
iedb_fig3_II = iedb_fig3[which(iedb_fig3$allele ==  "HLA-DRB1*01:01"),]

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

all_I = unique(c(iedb_fig3_I$peptide, Net_BA_cast$Peptide))
all_II = unique(c(iedb_fig3_II$peptide, Netii_BA_cast$Peptide))

Pop_freq_I_iedb = foreach(n=1:length(all_I), .combine=rbind) %dopar% {
  iedb = which(iedb_fig3_I$peptide == all_I[n])
  pred = which(Net_BA_cast$Peptide == all_I[n])
  
  if(length(iedb)>0 & length(pred)>0){
  Binding_haplotype=paste0(unique(colnames(iedb_fig3_I[iedb,c(5:13)])[which(iedb_fig3_I[iedb,c(5:13)] < 500)],
                           colnames(Net_BA_cast[pred,c(4:18)])[which(Net_BA_cast[pred,c(4:18)] < 300)]), collapse = ",")
  Frequency=1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% strsplit(Binding_haplotype,",")[[1]] )]/100,2)))
  Start= iedb_fig3_I$Start[iedb]
  Protein = iedb_fig3_I$Protein[iedb]
  }else if(length(iedb)>0 & length(pred)==0){
    Binding_haplotype=paste0(colnames(iedb_fig3_I[iedb,c(5:13)])[which(iedb_fig3_I[iedb,c(5:13)] < 500)], collapse = ",")
    Frequency=1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% strsplit(Binding_haplotype,",")[[1]] )]/100,2)))
    Start= iedb_fig3_I$Start[iedb]
    Protein = iedb_fig3_I$Protein[iedb]
  }else if(length(iedb)==0 & length(pred)>0){
    Binding_haplotype=paste0(colnames(Net_BA_cast[pred,c(4:18)])[which(Net_BA_cast[pred,c(4:18)] < 300)], collapse = ",")
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
    Binding_haplotype=paste0(unique(colnames(iedb_fig3_II[iedb,c(5)])[which(iedb_fig3_II[iedb,c(5)] < 500)],
                                    colnames(Netii_BA_cast[pred,c(4:18)])[which(Netii_BA_cast[pred,c(4:18)] < 245)]), collapse = ",")
    Frequency=1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% strsplit(Binding_haplotype,",")[[1]] )]/100,2)))
    Start= iedb_fig3_II$Start[iedb]
    Protein = iedb_fig3_II$Protein[iedb]
  }else if(length(iedb)>0 & length(pred)==0){
    Binding_haplotype=paste0(colnames(iedb_fig3_II[iedb,c(5)])[which(iedb_fig3_II[iedb,c(5)] < 500)], collapse = ",")
    Frequency=1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% strsplit(Binding_haplotype,",")[[1]] )]/100,2)))
    Start= iedb_fig3_II$Start[iedb]
    Protein = iedb_fig3_II$Protein[iedb]
  }else if(length(iedb)==0 & length(pred)>0){
    Binding_haplotype=paste0(colnames(Netii_BA_cast[pred,c(4:18)])[which(Netii_BA_cast[pred,c(4:18)] < 245)], collapse = ",")
    Frequency=1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% strsplit(Binding_haplotype,",")[[1]] )]/100,2)))
    Start= Netii_BA_cast$Start[pred]
    Protein = Netii_BA_cast$Protein[pred]
  }
  End=Start+nchar(all_II[n])-1
  Low_entropy = ifelse(max(ent$V2[Start:End]) > 0.1, 0, 1)
  data.frame(all_II[n],Binding_haplotype,Frequency,Protein,Start,End,Low_entropy)
}

Net_coepitopes$Low_entropy=NA
for (z in 1:nrow(Net_coepitopes)){
  Net_coepitopes$Low_entropy[z] = ifelse(max(ent$V2[Net_coepitopes$Start[z]:(Net_coepitopes$Start[z]+nchar(Net_coepitopes[z]$c2_peptide)-1)]) > 0.1, 0, 1)
}

##No filter
nrow(Pop_freq_I_iedb)
nrow(Pop_freq_II_iedb)
nrow(Net_coepitopes)

##Filter by entropy
filt_I = Pop_freq_I_iedb[which(Pop_freq_I_iedb$Low_entropy==1),]
filt_II = Pop_freq_II_iedb[which(Pop_freq_II_iedb$Low_entropy==1),]
filt_co = Net_coepitopes[which(Net_coepitopes$Low_entropy==1),]

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



###Plot combined, Figure 3 bottom#########

Master_tab=rbind(filt_I[,c(1,3,4,5,6,8,9)],
                filt_II[,c(1,3,4,5,6,8,9)],
                filt_co[, c(3,9,1,2,15,13,14)], use.names=F)
colnames(Master_tab)[1] = "Peptide"
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



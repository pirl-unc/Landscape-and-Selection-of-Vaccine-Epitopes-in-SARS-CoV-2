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


#######################################################################################
#######################Alternative figure set#######################################


###########################################
######Pre-processing#######################


#########Class I #########################
library(doMC)
registerDoMC(48)
###New Fig2A/B
Freq=fread(paste0(WORKING_DIR, "HLA_freq.txt"))

#A -- MHC-I

###Calculate pop freq for netMHCpan
I_8mer=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_8mer.xls")
I_9mer=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_9mer.xls")
I_10mer=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_10mer.xls")
I_11mer=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_11mers.xls")
I_B1501=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_B1501_8-11mers.txt")

HLAI=rbind(I_8mer, I_9mer, I_10mer, I_11mer)
HLAI=HLAI[order(HLAI$Peptide),]
I_B1501 = I_B1501[order(I_B1501$Peptide),]
HLAI=cbind(HLAI[,1:53], I_B1501[,4:8], HLAI[,54:ncol(HLAI)])

haps = c(seq(7,87,5))

resultsLog1<-foreach(n=1:nrow(HLAI), .combine=rbind) %dopar% {
  Binding_haplotype=paste0(colnames(HLAI)[haps[which(HLAI[n,haps+1,with=F] <= 1)]], collapse = ",")
  Frequency=1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% colnames(HLAI)[haps[which(HLAI[n,haps+1,with=F] <= 1)]])]/100,2)))
  data.frame(Binding_haplotype,Frequency)
}

Min<-foreach(n=1:nrow(HLAI), .combine=rbind) %dopar% {
  min = as.numeric(min(HLAI[n,haps,with=F]))
  data.frame(min)
}

HLAI$Min_nM = Min
HLAI = cbind(HLAI, resultsLog1)

fwrite(HLAI,"/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_all.txt", sep = '\t', col.names = T, quote = F)

####Calculate pop freq for MHCflurry
Flurry= fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/sars-cov-2.mhcflurry_predictions.csv")
Flurry=Flurry[,c(1,3,8,9)]
flurry=dcast.data.table(Flurry, peptide ~ best_allele, value.var = "affinity_percentile", fun.aggregate = mean)

##Taking too long...again in parallel
library(doMC)
registerDoMC(48)

resultsLog2<-foreach(n=1:nrow(flurry), .combine=rbind) %dopar% {
  
  Binding_haplotype = paste0(colnames(flurry)[which(flurry[n,] <= 1)], collapse = ",")
  probs = 1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% colnames(flurry)[which(flurry[n,] <= 1)])]/100,2)))
  data.frame(Binding_haplotype,probs)
  
}

flurry = cbind(flurry, resultsLog2)
colnames(flurry)[ncol(flurry)] = "Frequency"
fwrite(flurry, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/COVID_human_MHCflurry.txt", sep = '\t', col.names = T, quote = F)


###########Class II##################

###DR######################
Freq=fread(paste0(WORKING_DIR, "HLA_freq.txt"))
DRB=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_DRB1_NetMHCIIpan.xls")

haps = c(6,9,12,15,18,21,24)

resultsLog3<-foreach(n=1:nrow(DRB), .combine=rbind) %dopar% {
  
  Binding_haplotype = paste0(colnames(DRB)[haps[which(DRB[n,haps+1,with=F] <= 5)]], collapse = ",")
  Frequency=1-prod(1-(rep(Freq$US[which(Freq$Haplotype %in% colnames(DRB)[haps[which(DRB[n,haps+1,with=F] <= 5)]])]/100,2)))
  data.frame(Binding_haplotype,Frequency)
  
}
DRB = cbind(DRB, resultsLog3)
fwrite(DRB, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_DRB1_NetMHCIIpan.xls", col.names = T,  sep = '\t')



###DQ Cauc Am
Freq=fread(paste0(WORKING_DIR, "HLA_freq.txt"))

#DQ_cauc=fread(paste0(WORKING_DIR, "NetMHCpan_out_filt/COVID_Cauc_AM_NetMHCIIpan.xls"))

#Same but for unfiltered reads
DQ=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_Cauc_AM_NetMHCIIpan.xls")
haps = c(6,9,12,15,18,21,24,27)


resultsLog4<-foreach(n=1:nrow(DQ), .combine=rbind) %dopar% {
  
  Binding_haplotype = paste0(colnames(DQ)[haps[which(DQ[n,haps+1,with=F] <= 5)]], collapse = ",")
  Frequency=1-prod(1-(rep(Freq$Cauc_Am[which(Freq$Haplotype %in% colnames(DQ)[haps[which(DQ[n,haps+1,with=F] <= 5)]])]/100,2)))
  data.frame(Binding_haplotype,Frequency)
  
}
DQ = cbind(DQ, resultsLog4)
fwrite(DQ,"/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_Cauc_AM_NetMHCIIpan.xls", col.names = T, row.names = F, sep = '\t')



##########Combine class II
DR=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_DRB1_NetMHCIIpan.xls")
DQ=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_Cauc_AM_NetMHCIIpan.xls")

colnames(DR)[seq(7,25,3)] = paste( colnames(DR)[seq(7,25,3)],colnames(DR)[seq(6,24,3)] ,  sep = '_' )
colnames(DR)[seq(5,23,3)] = paste( colnames(DR)[seq(5,23,3)],colnames(DR)[seq(6,24,3)] ,  sep = '_' )

colnames(DQ)[seq(7,28,3)] = paste( colnames(DQ)[seq(7,28,3)],colnames(DQ)[seq(6,27,3)] ,  sep = '_' )
colnames(DQ)[seq(5,26,3)] = paste( colnames(DQ)[seq(5,26,3)],colnames(DQ)[seq(6,27,3)] ,  sep = '_' )


II_total = merge(DR, DQ, by = "Peptide", all = T)
II_total$Pos.x[which(is.na(II_total$Pos.x))] = II_total$Pos.y[which(is.na(II_total$Pos.x))]
colnames(II_total)[which(colnames(II_total) == "Pos.x")] = "Pos"
II_total$Pos.y = NULL

II_total$ID.x[which(is.na(II_total$ID.x))] = II_total$ID.y[which(is.na(II_total$ID.x))]
colnames(II_total)[which(colnames(II_total) == "ID.x")] = "ID"
II_total$ID.y = NULL


resultsLog5<-foreach(n=1:nrow(II_total), .combine=rbind) %dopar% {
  
  Min_nM = min(II_total$Min_nM.x[n], II_total$Min_nM.y[n], na.rm = T)
  
  Binding_haplotype = paste(str_replace_na(II_total$Binding_haplotype.x[n], replacement = ""), 
                            str_replace_na(II_total$Binding_haplotype.y[n], replacement = ""), sep = ',')
  
  Binding_haplotype=gsub('^\\,|\\.$', '', Binding_haplotype)
  Binding_haplotype=gsub("^,*|(?<=,),|,*$", "",  Binding_haplotype, perl=T)
  
  probs = c(II_total$Frequency.x[n], II_total$Frequency.y[n])
  probs = probs[!is.na(probs)]
  Frequency = 1-prod(1-probs)
  data.frame(Min_nM,Binding_haplotype,Frequency)
  
}
II_total = cbind(II_total, resultsLog5)


colnames(II_total)[which(colnames(II_total) == "Frequency.x")] = "Frequency_DR"
colnames(II_total)[which(colnames(II_total) == "Frequency.y")] = "Frequency_DQ"
II_total$Binding_haplotype.y = NULL
II_total$Binding_haplotype.x = NULL

colnames(II_total)[which(colnames(II_total) == "Min_nM.x")] = "Min_nM_DR"
colnames(II_total)[which(colnames(II_total) == "Min_nM.y")] = "Min_nM_DQ"

fwrite(II_total, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/COVID_human_netMHCIIpan_all.txt", sep = '\t', col.names = T, row.names = F, quote = F)


####Read in MARIA data
MARIA = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MARIA_results.tsv")
colnames(MARIA)[10] = "Peptide"
colnames(MARIA)[1] = "Allele"
MARIA$`MARIA percentile scores` = 100-MARIA$`MARIA percentile scores`

MARIA$Allele = str_remove_all(str_remove_all(MARIA$Allele, pattern = ":"), pattern = "\\*")
MARIA$Allele = str_replace_all(MARIA$Allele, pattern = "HLA-DR", replacement = "DR")
MARIA$Allele = str_replace_all(MARIA$Allele, pattern = "DQA", replacement = "HLA-DQA")
MARIA$Allele = str_replace_all(MARIA$Allele, pattern = "DRB1", replacement = "DRB1_")

MARIA = MARIA[which(MARIA$Allele %in% Freq$Haplotype),]

maria=dcast.data.table(MARIA, Peptide ~ Allele, value.var = "MARIA percentile scores", fun.aggregate = mean)

Freq$US[which(!is.na(Freq$Cauc_Am))] = Freq$Cauc_Am[which(!is.na(Freq$Cauc_Am))]

resultsLog6<-foreach(n=1:nrow(maria), .combine=rbind) %dopar% {
  
  Binding_haplotype = paste0(colnames(maria)[which(maria[n,] <= 5)], collapse = ",")
  prob=rep(Freq$US[which(Freq$Haplotype %in% colnames(maria)[which(maria[n,] <= 5)])]/100,2)
  prob = prob[which(!is.na(prob))]
  Frequency = 1-prod(1-(  prob ))
  
  data.frame(Binding_haplotype,Frequency)
  
}

maria = cbind(maria, resultsLog6)
fwrite(maria, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MARIA_cast.txt")

###################Apply filters for class I and II matrices###########################

###Class I####

Flurry= fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/sars-cov-2.mhcflurry_predictions.csv")
Flurry=Flurry[,c(1:3,8,9)]
Flurry$sequence_name = sapply(strsplit(Flurry$sequence_name, "\\|"),'[',1)
Flurry = Flurry[-which(duplicated(Flurry)),] ##Removed duplicated B44:03 haplotype
colnames(Flurry) = c("Protein", "Start", "Peptide", "Haplotype", "MHCflurry_percentile")
Flurry$Protein[which(Flurry$Protein == "orflab")] = "orf1ab"
Flurry = Flurry[order(Flurry$Peptide),]
fwrite(Flurry, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCflurry_standardized_all.txt")

Net = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_all.txt")
colnames(Net)[seq(8,86,5)] = colnames(Net)[seq(7,86,5)]
colnames(Net)[seq(7,86,5)] = paste0("nM_", colnames(Net)[seq(7,86,5)])
Net = melt(Net, id.vars = 2, measure.vars = c(seq(8,86,5)))
colnames(Net) = c("Peptide", "Haplotype", "NetMHCpan_percentile")
Net = Net[order(Net$Peptide),]
Net = cbind(Flurry$Protein, Flurry$Start, Net)
colnames(Net)[1:2] = c("Protein", "Start")
fwrite(Net, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_all.txt")

###Filtered
Net = fread( "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_all.txt")
Flurry = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCflurry_standardized_all.txt")


Net_cast = dcast.data.table(Net, Peptide ~ Haplotype, value.var = "NetMHCpan_percentile")

Net = Net[which(Net$NetMHCpan_percentile <= 1),]
Flurry = Flurry[which(Flurry$MHCflurry_percentile <= 1),]

Net = Net[which(paste0(Net$Peptide,Net$Haplotype) %in% paste0(Flurry$Peptide,Flurry$Haplotype)),]

Net_cast_filt = Net_cast[which(Net_cast$Peptide %in% Net$Peptide),]

#flurry=dcast.data.table(Flurry, peptide ~ best_allele, value.var = "affinity_percentile", fun.aggregate = mean)


###Class II
library(doMC)
registerDoMC(64)

Maria = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MARIA_results.tsv")
Maria=Maria[,c(3,10,1,9)]
colnames(Maria) = c("Protein", "Peptide", "Haplotype", "MARIA_percentile")
Maria = Maria[-which(duplicated(Maria)),]
Maria = Maria[order(Maria$Peptide),]

Maria$MARIA_percentile = 100-Maria$MARIA_percentile  ###Change to same order as other sets (small == high binder)
Maria$Haplotype = str_remove_all(str_remove_all(Maria$Haplotype, pattern = ":"), pattern = "\\*")
Maria$Haplotype = str_replace_all(Maria$Haplotype, pattern = "HLA-DR", replacement = "DR")
Maria$Haplotype = str_replace_all(Maria$Haplotype, pattern = "DQA", replacement = "HLA-DQA")
Maria$Haplotype = str_replace_all(Maria$Haplotype, pattern = "DRB1", replacement = "DRB1_")

seq =  fread(paste0(WORKING_DIR, "AA_sequence_combined.txt"), header=F)
Start = c()
###Checking if string is present in SARS-CoV-2 reference.  If so, add coordinates####

for( n in 1:nrow(Maria)){
  if(str_detect(string = seq$V1 ,Maria$Peptide[y])){
    Start = c(Start,  gregexpr(text = seq$V1, pattern = paper_tc$Peptide[y])[[1]][1])
  }else{
    print(y)
    Start = c(Start, NA)
  }
  data.frame(Start)
}

Maria = cbind(Maria$Protein, rep(NA, nrow(Maria)), Maria[,2:4])
colnames(Maria)[2] = "Start"
colnames(Maria)[1] = "Protein"

fwrite(Maria, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MARIA_standardized_all.txt")


Netii = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/COVID_human_netMHCIIpan_all.txt")
colnames(Netii)[c(seq(7,25,3), seq(32,54,3))] = substr(colnames(Netii)[c(seq(7,25,3), seq(32,54,3))],6, nchar(colnames(Netii)[c(seq(7,25,3), seq(32,54,3))]))
#colnames(Netii)[seq(7,86,5)] = paste0("nM_", colnames(Netii)[seq(7,86,5)])
Netii = melt(Netii, id.vars = 1, measure.vars = c(seq(7,25,3), seq(32,54,3)))
colnames(Netii) = c("Peptide", "Haplotype", "NetMHCIIpan_percentile")
Netii = Netii[order(Netii$Peptide),]

###Add in Proteins and position data
netii = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/COVID_human_netMHCIIpan_all.txt")
netii$ID = sapply(strsplit(netii$ID,"_"),'[',1)
netii$ID[which(netii$ID == "orflab")] = "orf1ab"

Pos=c()
Prot=c()
for(n in 1:nrow(Netii)){
  print(n)
  Pos = c(Pos, netii$Pos[which(netii$Peptide == Netii$Peptide[n])[1]])
  Prot = c(Prot, netii$ID[which(netii$Peptide == Netii$Peptide[n])[1]])
}

Netii = cbind(Prot, Pos, Netii)
colnames(Netii)[1:2] = c("Protein", "Start")
fwrite(Netii, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_all.txt")


####Defining co-epitopes

#Perfect overlap


####IEDB vs in silico

##Old format
iedb = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/iedb/mhc-sars2.csv")
colnames(iedb) = unlist(iedb[1,])
iedb = iedb[-1,]

netii = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/COVID_human_netMHCIIpan_all.txt")
netii_filt = netii[which(netii$Peptide %in% iedb$Description),]

fwrite(netii_filt, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/COVID_human_netMHCIIpan_IEDB_filt.txt")

net = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/COVID_human_netMHCpan_rank.txt")
net_filt = net[which(net$Peptide %in% iedb$Description),]

fwrite(net_filt, "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/COVID_human_netMHCpan_IEDB_filt.txt")

##New format
Net = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_all.txt")
Netii = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_all.txt")

Net_unique = paste(Net$Peptide, Net$Haplotype, sep = "_")
Netii_unique = paste(Netii$Peptide, Netii$Haplotype, sep = "_")

iedb$`Allele Name` = str_replace_all(iedb$`Allele Name`, pattern = "\\*", replacement = "")
iedb$`Allele Name` = str_replace_all(iedb$`Allele Name`, pattern = "HLA-DRB101:01", replacement = "DRB1_0101")
iedb_unique = paste(iedb$Description, iedb$`Allele Name`, sep = "_")

Net_filt = Net[which(Net_unique %in% iedb_unique),]
Netii_filt = Netii[which(Netii_unique %in% iedb_unique),]

Net_filt$IEDB_qual = NA
Netii_filt$IEDB_qual = NA

for(n in 1:nrow(Net_filt)){
  Net_filt$IEDB_qual[n] = paste(unique(iedb$`Qualitative Measure`[which(iedb_unique == paste(Net_filt$Peptide[n], Net_filt$Haplotype[n], sep = "_"))]), collapse = ",")
}
for(n in 1:nrow(Netii_filt)){
  Netii_filt$IEDB_qual[n] = paste(unique(iedb$`Qualitative Measure`[which(iedb_unique == paste(Netii_filt$Peptide[n], Netii_filt$Haplotype[n], sep = "_"))]), collapse = ",")
}

Net_filt$Negative = ifelse(Net_filt$IEDB_qual == "Negative","Negative","Positive")

ggplot(data = Net_filt)+
  geom_violin(aes(x = factor(Negative), y= NetMHCpan_percentile, color = factor(Negative)))+
  geom_quasirandom(aes(x = factor(Negative), y= NetMHCpan_percentile, color = factor(Negative)), alpha=.5)+
  scale_y_continuous(limits = c(-1,100), trans = ssqrt_trans)+
  scale_color_viridis_d(name = "Negative IEDB\nqualitative binding")+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  geom_hline(yintercept = 1, color="red")+
  labs(x="Negative binding (IEDB)", y="NetMHCpan Percentile", title = "HLA-I predicted epitopes") 

ggplot(data = Netii_filt)+
  geom_violin(aes(x = factor(IEDB_qual), y= NetMHCIIpan_percentile, color = factor(IEDB_qual)))+
  geom_quasirandom(aes(x = factor(IEDB_qual), y= NetMHCIIpan_percentile, color = factor(IEDB_qual)), alpha=.5)+
  scale_y_continuous(limits = c(-1,100), trans = ssqrt_trans)+
  scale_color_viridis_d(name = "IEDB\nqualitative binding")+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  geom_hline(yintercept = 5, color="red")+
  labs(x="Qualitative binding (IEDB)", y="NetMHCIIpan Percentile", title = "HLA-II predicted epitopes")   


####################################################
#########NetMHC(x)pan overlap#######################
Net = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_all.txt")
Netii = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_all.txt")

Net$End = Net$Start + nchar(Net$Peptide) -1
Netii$End = Netii$Start + nchar(Netii$Peptide) -1

Net = Net[which(Net$NetMHCpan_percentile<=1),]
Netii = Netii[which(Netii$NetMHCIIpan_percentile<=1),]

Net_unique = Net[-which(duplicated(Net$Peptide)),]
Netii_unique = Netii[-which(duplicated(Netii$Peptide)),]

#Overlaps
g1 <-  GRanges(seqnames="COVID",
               IRanges(start=Net_unique$Start,
                       end=Net_unique$End), 
               I_peptide = Net_unique$Peptide, MHCI_percentile = Net_unique$NetMHCpan_percentile, Protein = Net_unique$Protein)

g2 <-  GRanges(seqnames="COVID",
               IRanges(start=Netii_unique$Start,
                       end=Netii_unique$End), 
               II_peptide = Netii_unique$Peptide, MHCII_percentile = Netii_unique$NetMHCIIpan_percentile)

merge_I_II=as.data.table(mergeByOverlaps(g1,g2))


ggplot(data=merge_I_II) + 
  ##Plot co-epitopes
  geom_point(aes(x=MHCI_percentile, y=MHCII_percentile, color=Protein)) +
  theme(text=element_text(face="bold",size=20,colour="black")) +
  scale_color_viridis_d(option = "viridis")+
  geom_label_repel(aes(label=ifelse(Mean_freq>label_frequency,as.character(c2_peptide),''), x=c1_tot_freq, y=c2_tot_freq),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  scale_y_continuous(limits = c(0,95))+
  scale_x_continuous(limits = c(0,95))+
  scale_size(range = c(0,10), breaks = c(10,150,300), trans="reverse", name = "Mean nM")+
  labs(y="MHC II Frequency (%)", x="MHC I Frequency (%)", title = "Predicted T cell co-epitopes", color="Gene") +
  
  ###Add MHCI and II epitopes
  geom_rug(data=net_graph, aes(x = 100*Frequency),  col=viridis(n = 2, begin = 0, end = .75, option = 'plasma')[1])+
  geom_rug(data=netii_graph, aes(y = 100*Frequency),  col=viridis(n=2, begin = 0, end = .75, option = 'plasma')[2])+
  
  geom_text(x=75, y=96, label=paste0("Co-epitopes: ", nrow(joint_run_filt_entropy)), color = 'red', size=10)+
  geom_text(x=75, y=3, label=paste0("HLA-I: ", nrow(net_graph)), color = 'red', size=10)+
  geom_text(x=3, y=75, label=paste0("HLA-II: ", nrow(netii_graph)), color = 'red', angle=90, size=10)+
  
  theme(text=element_text(face="bold",size=20,colour="black")) 



###############################################
#####Figure 2A##################################

###Top
Flurry = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MHCflurry_standardized_all.txt")
Flurry$Peptide = paste(Flurry$Peptide, Flurry$Haplotype, sep = "_")
Flurry$Haplotype = NULL

Net = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_standardized_all.txt")
Net$Peptide = paste(Net$Peptide, Net$Haplotype, sep = "_")
Net$Haplotype = NULL

combined_perc = merge(x = Flurry, y = Net, by = "Peptide", all=T)
combined_perc$Haplotype = factor(sapply(strsplit(combined_perc$Peptide, "_"),"[",2))
combined_perc$Peptide = sapply(strsplit(combined_perc$Peptide,"_"),'[',1)

combined_perc = combined_perc[,c(2,3,1,8,4,7)]
colnames(combined_perc)[1:2] = c("Protein", "Start")
fwrite(combined_perc,"/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLA-I_standardized_all.txt")

combined_perc$Length = nchar(combined_perc$Peptide)

#group = unique(combined_perc$Peptide[which(combined_perc$NetMHCpan_percentile <=1)])
#group = group[which(group %in% sapply(strsplit(Flurry$Peptide[which(Flurry$MHCflurry_percentile <=1)],"_"),'[',1) )]
#net = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCpan_out/COVID_NetMHCpan_all.txt")
#group = group[which(group %in% net$Peptide[which(net$Min_nM<500)])]
#combined_perc$group = 0
#combined_perc$group[which(combined_perc$Peptide %in% group)] = 1
#group_mat = combined_perc[which(combined_perc$group==1),]
#group_mat = group_mat[which(group_mat$NetMHCpan_percentile<=1),]

ggplot(data=combined_perc[sample(nrow(combined_perc), 10000),]) + 
  geom_point(aes(x=NetMHCpan_percentile, y=MHCflurry_percentile, color = as.character(Length)))+
  scale_y_continuous(limits = c(-1,100), trans = ssqrt_trans)+
  scale_x_continuous(limits = c(-1,100), trans = ssqrt_trans)+
  scale_color_viridis_d(name = "Peptide length", alpha = .5)+
 # geom_density_2d(data = group_mat, aes(x=NetMHCpan_percentile, y = MHCflurry_percentile), color=alpha(colour = "red", alpha = .7))+
  geom_vline(xintercept = 1, color = "red")+
  geom_hline(yintercept = 1, color = 'red')+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  geom_text(x=-.5, y=-.5, label=paste0("n=",length(unique(group_mat$Peptide))), color = 'red', size=5)+
  geom_text(x=50, y=-.5, label=paste0("Total: ", nrow(combined_perc)), color = 'red', size=5)+
  #geom_text(x=0, y=-3, label=paste0("n=",length(which((combined_perc$MHCflurry_percentile<1)& (combined_perc$NetMHCpan_percentile<1) ))), color = 'red', size=10)+
  
  labs(x="NetMHCpan Percentile", y="MHCflurry Percentile", title = "HLA-I predicted epitopes") 




###Bottom
Maria = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/MARIA_standardized_all.txt")
Maria$Peptide = paste(Maria$Peptide, Maria$Haplotype, sep = "_")
Maria$Haplotype = NULL

Netii = fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/NetMHCIIpan_standardized_all.txt")
Netii$Peptide = paste(Netii$Peptide, Netii$Haplotype, sep = "_")
Netii$Haplotype = NULL

combined_perc = merge(x = Maria, y = Netii, by = "Peptide", all=T)
combined_perc$Haplotype = substr(combined_perc$Peptide, 17,nchar(combined_perc$Peptide))
combined_perc$Peptide = sapply(strsplit(combined_perc$Peptide,"_"),'[',1)

combined_perc$Protein.x[which(is.na(combined_perc$Protein.x))] = combined_perc$Protein.y[which(is.na(combined_perc$Protein.x))] 
combined_perc$Start.x[which(is.na(combined_perc$Start.x))] = combined_perc$Start.y[which(is.na(combined_perc$Start.x))]
combined_perc = combined_perc[, c(2,3,1,8,4,7)]
colnames(combined_perc)[1:2] = c("Protein", "Start")
fwrite(combined_perc,"/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/HLA-II_standardized_all.txt")


combined_perc$Length = nchar(combined_perc$Peptide)
combined_perc = combined_perc[complete.cases(combined_perc),]

ggplot(data=combined_perc)+ 
  geom_point(aes(x=NetMHCIIpan_percentile, y=MARIA_percentile, color = as.character(Length)))+
  scale_y_continuous(limits = c(-1,100), trans = ssqrt_trans)+
  scale_x_continuous(limits = c(-1,100), trans = ssqrt_trans)+
  scale_color_viridis_d(name = "Peptide length", alpha = .5)+
  geom_vline(xintercept = 1, color = "red")+
  geom_hline(yintercept = 1, color = 'red')+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  geom_text(x=-.5, y=-.5, label=paste0("n=",length(which((combined_perc$MARIA_percentile<=1) & (combined_perc$NetMHCIIpan_percentile<=1) ))), color = 'red', size=5)+
  geom_text(x=0, y=-3, label=paste0("Total: ", nrow(combined_perc)), color = 'red', size=10)+
  #geom_text(x=0, y=-3, label=paste0("n=",length(which((combined_perc$MHCflurry_percentile<1)& (combined_perc$NetMHCpan_percentile<1) ))), color = 'red', size=10)+
  
  labs(x="NetMHCIIpan Percentile", y="MARIA Percentile", title = "HLA-II predicted epitopes") 




#################################################
####################New Fig 2B ##################


#Minimum frequency to display (mean of HLA-I and II frequencies)  
Min_frequency = 0
label_frequency = 75

mouse1 = fread(paste0(WORKING_DIR, "COVID_murine_NetMHCpan_unfilt.xls"))
mouse1 = mouse1[which(mouse1$NB>0),]
mouse2 = fread(paste0(WORKING_DIR, "COVID_murine_NetMHCIIpan_unfilt.xls"))
mouse2 = mouse2[which(mouse2$NB>0),]

#joint_run_filt_entropy = fread(paste0(WORKING_DIR, "COVID_combined_entropy_rank_filtered.csv"))
joint_run_filt_entropy=fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/COVID_human_Coepitopes_rank.txt")
joint_run_filt_entropy$c1_tot_freq = joint_run_filt_entropy$c1_tot_freq*100
joint_run_filt_entropy$c2_tot_freq = joint_run_filt_entropy$c2_tot_freq*100

###Number of MHC-I/II epitopes conserved between human and mouse
#length(which(HLAI_filt$Peptide[which(HLAI_filt$lo_entropy==1)] %in% unique(mouse1$Peptide)))
#length(which(HLAII_filt$Peptide[which(HLAII_filt$lo_entropy==1)] %in% unique(mouse2$Peptide)))

##Adding murine coverage for each human epitope
joint_run_filt_entropy$murine = "None"
joint_run_filt_entropy$murine[which(joint_run_filt_entropy$c1_peptide %in% mouse1$Peptide)] = "MHC-I" 
joint_run_filt_entropy$murine[which(joint_run_filt_entropy$c2_peptide %in% mouse2$Peptide)] = "MHC-II" 
joint_run_filt_entropy$murine[which((joint_run_filt_entropy$c2_peptide %in% mouse2$Peptide) & (joint_run_filt_entropy$c1_peptide %in% mouse1$Peptide))] = "Both" 
joint_run_filt_entropy$murine = factor(joint_run_filt_entropy$murine, levels = c("MHC-I", "MHC-II", "Both", "None"))

#joint_run_filt_entropy = joint_run_filt_entropy[which(joint_run_filt_entropy$lo_entropy==1),]
joint_run_filt_entropy$Mean_freq = (joint_run_filt_entropy$c1_tot_freq + joint_run_filt_entropy$c2_tot_freq)/2

joint_run_filt_entropy$c1_peptide_offset=NA
for(z in 1:nrow(joint_run_filt_entropy)){
  #Uses gregexpr to find c1 peptide within c2 peptide --> inserts spaces ahead of c1 peptide so the start AA aligns with its location on c2 peptide --> Pastes out as a string in matrix
  joint_run_filt_entropy$c1_peptide_offset[z] = paste0(c(rep(" ", (gregexpr(joint_run_filt_entropy$c2_peptide[z], pattern = joint_run_filt_entropy$c1_peptide[z])[[1]][1])-1), joint_run_filt_entropy$c1_peptide[z]), collapse = '')
}

####Read in filtered HLA-I/II data
net_filt= fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/COVID_human_netMHCpan_rank.txt")
netii_filt= fread("/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Working/COVID_human_netMHCIIpan_rank.txt")

net_filt=net_filt[which(net_filt$Frequency>0),]
netii_filt=netii_filt[which(netii_filt$Frequency>0),]

joint_run_filt_entropy=joint_run_filt_entropy[which((joint_run_filt_entropy$c1_tot_freq>0)&(joint_run_filt_entropy$c2_tot_freq>0)),]

###Keep only columns needed to graph to prevent ggplot error
net_graph = net_filt[,c(1,2,3,87:88)]
netii_graph = netii_filt[,c(1,2,3,57:59)]
#net_graph = net_graph[which(net_graph$lo_entropy==1),]
#netii_graph = netii_graph[which(netii_graph$lo_entropy==1),]

pep_col = viridis(4, begin = 0, end = .75, option = "plasma")

ggplot(data=joint_run_filt_entropy) + 
  ##Plot co-epitopes
  geom_point(aes(x=c1_tot_freq, y=c2_tot_freq, color=gene, shape = murine),size=7, alpha=.5) +
  theme(text=element_text(face="bold",size=20,colour="black")) +
  scale_color_viridis_d(option = "viridis")+
  geom_label_repel(aes(label=ifelse(Mean_freq>label_frequency,as.character(c2_peptide),''), x=c1_tot_freq, y=c2_tot_freq),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  scale_y_continuous(limits = c(0,95))+
  scale_x_continuous(limits = c(0,95))+
  scale_size(range = c(0,10), breaks = c(10,150,300), trans="reverse", name = "Mean nM")+
  labs(y="MHC II Frequency (%)", x="MHC I Frequency (%)", title = "Predicted T cell co-epitopes", color="Gene") +
  
  ###Add MHCI and II epitopes
  geom_rug(data=net_graph, aes(x = 100*Frequency),  col=viridis(n = 2, begin = 0, end = .75, option = 'plasma')[1])+
  geom_rug(data=netii_graph, aes(y = 100*Frequency),  col=viridis(n=2, begin = 0, end = .75, option = 'plasma')[2])+
  
  geom_text(x=75, y=96, label=paste0("Co-epitopes: ", nrow(joint_run_filt_entropy)), color = 'red', size=10)+
  geom_text(x=75, y=3, label=paste0("HLA-I: ", nrow(net_graph)), color = 'red', size=10)+
  geom_text(x=3, y=75, label=paste0("HLA-II: ", nrow(netii_graph)), color = 'red', angle=90, size=10)+
  
  theme(text=element_text(face="bold",size=20,colour="black")) 


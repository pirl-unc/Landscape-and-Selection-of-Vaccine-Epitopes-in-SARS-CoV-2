############################################################
#############Dependencies and working directory############
##########################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install()

packages <- c("data.table")

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  BiocManager::install(setdiff(packages, rownames(installed.packages())))  
}

library(data.table)

WORKING_DIR = "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2/COVID/"


###########################
###DR######################
Freq=fread(paste0(WORKING_DIR, "HLA_freq.txt"))
DRB=fread(paste0(WORKING_DIR, "NetMHCpan_out_filt/COVID_DRB1_NetMHCIIpan.xls"))

haps = c(6,9,12,15,18,21,24)

DRB$Binding_haplotype = NA
DRB$Frequency = NA

for (n in 1:nrow(DRB)){
  print (n)
  DRB$Binding_haplotype[n] = paste0(colnames(DRB)[haps[which(DRB[n,haps+1,with=F] <= 5)]], collapse = ",")
  probs = rep(Freq$US[which(Freq$Haplotype %in% colnames(DRB)[haps[which(DRB[n,haps+1,with=F] <= 5)]])]/100,2)
  probs = probs[!is.na(probs)]
  
  DRB$Frequency[n] = 1-prod(1-probs)
  
}

fwrite(DRB,paste0(WORKING_DIR, "NetMHC_pan_out_rank/COVID_DRB1_NetMHCIIpan.txt"), col.names = T, row.names = F, sep = '\t')



###DQ Cauc Am
Freq=fread(paste0(WORKING_DIR, "HLA_freq.txt"))

DQ_cauc=fread(paste0(WORKING_DIR, "NetMHCpan_out_filt/COVID_Cauc_AM_NetMHCIIpan.xls"))

haps = c(6,9,12,15,18,21,24,27)

DQ_cauc$Binding_haplotype = NA
DQ_cauc$Frequency = NA

for (n in 1:nrow(DQ_cauc)){
  print (n)
  DQ_cauc$Binding_haplotype[n] = paste0(colnames(DQ_cauc)[haps[which(DQ_cauc[n,haps+1,with=F] <= 5)]], collapse = ",")
  
  probs = rep(Freq$Cauc_Am[which(Freq$Haplotype %in% colnames(DQ_cauc)[haps[which(DQ_cauc[n,haps+1,with=F] <= 5)]])]/100,2)
  probs = probs[!is.na(probs)]
  
  DQ_cauc$Frequency[n] = 1-prod(1-probs)
}

fwrite(DQ_cauc,paste0(WORKING_DIR, "NetMHC_pan_out_rank/COVID_Cauc_AM_NetMHCIIpan.txt"), col.names = T, row.names = F, sep = '\t')



###DQ World
Freq=fread(paste0(WORKING_DIR, "HLA_freq.txt"))

DQ_cauc=fread(paste0(WORKING_DIR, "NetMHCpan_out_filt/COVID_DQ_worldwide_NetMHCIIpan.xls"))

haps = c(6,9,12,15,18,21)

DQ_cauc$Binding_haplotype = NA
DQ_cauc$Frequency = NA

for (n in 1:nrow(DQ_cauc)){
  print (n)
  DQ_cauc$Binding_haplotype[n] = paste0(colnames(DQ_cauc)[haps[which(DQ_cauc[n,haps+1,with=F] <= 5)]], collapse = ",")
  #DQ_cauc$Frequency[n] = sum()
  probs = rep(Freq$World[which(Freq$Haplotype %in% colnames(DQ_cauc)[haps[which(DQ_cauc[n,haps+1,with=F] <= 5)]])]/100,2)
  probs = probs[!is.na(probs)]
  
  DQ_cauc$Frequency[n] = 1-prod(1-probs)
}

fwrite(DQ_cauc,paste0(WORKING_DIR, "NetMHC_pan_out_rank/COVID_DQ_worldwide_NetMHCIIpan.txt"), col.names = T, row.names = F, sep = '\t')




###Class I
Freq=fread(paste0(WORKING_DIR, "HLA_freq.txt"))

I_8mer=fread(paste0(WORKING_DIR, "NetMHCpan_out_filt/COVID_MHCI_8mer_NetMHCpan.txt"))
I_9mer=fread(paste0(WORKING_DIR, "NetMHCpan_out_filt/COVID_MHCI_9mer_NetMHCpan.txt"))
I_10mer=fread(paste0(WORKING_DIR, "NetMHCpan_out_filt/COVID_MHCI_10mer_NetMHCpan.txt"))
I_11mer=fread(paste0(WORKING_DIR, "NetMHCpan_out_filt/COVID_MHCI_11mer_NetMHCpan.txt"))

HLAI=rbind(I_8mer, I_9mer, I_10mer, I_11mer)


haps = c(seq(8,87,5))

HLAI$Binding_haplotype = NA
HLAI$Frequency = NA

for (n in 1:nrow(HLAI)){
  print (n)
  HLAI$Binding_haplotype[n] = paste0(colnames(HLAI)[haps[which(HLAI[n,haps+1,with=F] <= 1)]], collapse = ",")
  probs = rep(Freq$US[which(Freq$Haplotype %in% colnames(HLAI)[haps[which(HLAI[n,haps+1,with=F] <= 1)]])]/100,2)
  probs = probs[!is.na(probs)]
  HLAI$Frequency[n] = 1-prod(1-probs)
}

fwrite(HLAI,paste0(WORKING_DIR, "COVID_human_netMHCpan_rank.txt"), col.names = T, row.names = F, sep = '\t')

######Condensed down the raw MHCFlurry outputs to only those with a haplotype in the top 1st percentile
#Flurry= fread(paste0(WORKING_DIR, "sars-cov-2.mhcflurry_predictions.csv"))
#Flurry=Flurry[,c(1,3,8,9)]
#flurry=dcast.data.table(Flurry, peptide ~ best_allele, value.var = "affinity_percentile", fun.aggregate = sum)
#Min_aff = c()
#for(r in 1:nrow(flurry)){
#  #print(r)
#  Min_aff= c(Min_aff,min(flurry[r,2:ncol(flurry)]))
#}
#flurry=flurry[which(Min_aff <= 1),]
#fwrite(flurry, paste0(WORKING_DIR, "filtered_cast_sars-cov-2.mhcflurry_predictions.csv"), col.names = T, quote = F)
flurry=fread( paste0(WORKING_DIR, "filtered_cast_sars-cov-2.mhcflurry_predictions.csv"))

HLAI_filt = HLAI[which(HLAI$Frequency>0),]
HLAI_filt = HLAI_filt[which(HLAI_filt$Peptide %in% flurry$peptide),]

fwrite(HLAI_filt,paste0(WORKING_DIR, "COVID_human_netMHCpan_rank_filtered.txt"), col.names = T, row.names = F, sep = '\t')

##########Combine all class II calls######
library(stringr)
DR=fread(paste0(WORKING_DIR, "NetMHC_pan_out_rank/COVID_DRB1_NetMHCIIpan.txt"))
DQ=fread(paste0(WORKING_DIR, "NetMHC_pan_out_rank/COVID_Cauc_AM_NetMHCIIpan.txt"))

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

II_total$Min_nM = NA
II_total$Binding_haplotype = NA
II_total$Frequency = NA

for(n in 1:nrow(II_total)){
  print(n)
  
  II_total$Min_nM[n] = min(II_total$Min_nM.x[n], II_total$Min_nM.y[n], na.rm = T)
  
  II_total$Binding_haplotype[n] = paste(str_replace_na(II_total$Binding_haplotype.x[n], replacement = ""), 
                                        str_replace_na(II_total$Binding_haplotype.y[n], replacement = ""), sep = ',')
  
  II_total$Binding_haplotype[n]=gsub('^\\,|\\.$', '', II_total$Binding_haplotype[n])
  II_total$Binding_haplotype[n]=gsub("^,*|(?<=,),|,*$", "",  II_total$Binding_haplotype[n], perl=T)
  
  probs = c(II_total$Frequency.x[n], II_total$Frequency.y[n])
  probs = probs[!is.na(probs)]
  II_total$Frequency[n] = 1-prod(1-probs)
  
}

colnames(II_total)[which(colnames(II_total) == "Frequency.x")] = "Frequency_DR"
colnames(II_total)[which(colnames(II_total) == "Frequency.y")] = "Frequency_DQ"
II_total$Binding_haplotype.y = NULL
II_total$Binding_haplotype.x = NULL

colnames(II_total)[which(colnames(II_total) == "Min_nM.x")] = "Min_nM_DR"
colnames(II_total)[which(colnames(II_total) == "Min_nM.y")] = "Min_nM_DQ"

fwrite(II_total, paste0(WORKING_DIR, "COVID_human_netMHCIIpan_rank.txt"), sep = '\t', col.names = T, row.names = F, quote = F)


MARIA = fread(paste0(WORKING_DIR, "COVID_human_MARIA.txt"))


II_filt = II_total[which(II_total$Frequency>0),]
II_filt = II_filt[which(II_filt$Peptide %in% MARIA$Peptide),]

fwrite(II_filt,paste0(WORKING_DIR, "COVID_human_netMHCIIpan_rank_filtered.txt"), col.names = T, row.names = F, sep = '\t')


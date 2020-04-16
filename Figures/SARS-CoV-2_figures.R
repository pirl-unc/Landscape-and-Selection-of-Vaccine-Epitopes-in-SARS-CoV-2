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

<<<<<<< HEAD
WORKING_DIR = "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2_epitope_landscape/Figures/COVID/"
=======
WORKING_DIR = "/datastore/nextgenout5/share/labs/Vincent_Lab/datasets/SARS-CoV-2/COVID/"
>>>>>>> 8641edb347a0984de8f29feccce730866780fec5


#############################################################################################
###################Supplemental Figure X -- entropy vs proportion of top AA###################

##Most code is commented out to save space on the source files.  They're kept here to show steps, but the final table is provided for plotting
# library(seqinr)
# 
# msa = readRDS(paste0(WORKING_DIR, "nucl_mafft_df.rds"))  #MSA output
# rownames(msa) = msa[,1] #Formating
# msa=msa[,-1]
# trans_tab = matrix(nrow=nrow(msa), ncol=9701)
# rownames(trans_tab) = rownames(msa)
# 
# trans_coord = fread(paste0(WORKING_DIR, "Translation_coords.txt")) #Read in table of translational start/end sites, derived from NCBI
# 
# #Remove stop codons (derived from checking code in SARS-CoV-2 reference)
# remove= c( 7097, 8371, 8647, 8723, 8946, 9008, 9130, 9252, 9672, 9711)
# 
# #For each sample, translate sequences into AA and put into a new matrix
# for(m in 1:nrow(msa)){
#   print(m)
#   AA_temp=c()
#   for(n in 1:nrow(trans_coord)){
#     seq_temp = s2c(tolower(paste0(lapply(msa[m,trans_coord$Start[n]:trans_coord$End[n]],as.character),collapse="")))
#     AA_temp=c(AA_temp,seqinr::translate(seq_temp))
#   }
#   AA_temp = AA_temp[-remove]  
#   trans_tab[m,] = AA_temp
# }
# 
# #Formating
# trans_tab_fin = cbind(rownames(trans_tab), trans_tab)
# colnames(trans_tab_fin) = c("SeqID", seq(1,9701,1))
# 
# fwrite(trans_tab_fin,paste0(WORKING_DIR, "AA_mafft_df.txt"), col.names = T, row.names = F, sep = '\t', quote = F )
# 
#
##Plot entropy vs AA proportion
#
#ent = fread(paste0(WORKING_DIR, "entropy_7882.txt"))
#
#trans_tab = fread(paste0(WORKING_DIR, "total.augur.align.dedup.faa.tab.C"), header = T)
#
#prop = c() #Proportion of most common AA per location
# mode = c() #Most common AA per location
# for (z in 2:ncol(trans_tab)){
#   print(z)
#   vals = unlist(trans_tab[,z, with=F])
#   vals=vals[-which((vals=="*")|(vals=="X"))]
#   mode = c(mode, names(table(vals))[as.vector(table(vals))==max(table(vals))])
#   prop = c(prop, max(table(vals))/length(vals) )
# }
# 
# ent_vs_prop = cbind(ent, prop)  
# fwrite(ent_vs_prop, paste0(WORKING_DIR, "ent_vs_prop.txt"))

ent_vs_prop = fread( paste0(WORKING_DIR, "ent_vs_prop.txt"))
<<<<<<< HEAD
colnames(ent_vs_prop) = c("Position", "Entropy", "Proportion")

ggplot(data=ent_vs_prop)+
  geom_point(aes(x=Entropy, y = Proportion))+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  geom_vline(xintercept = 0.1, color="red")


=======
>>>>>>> 8641edb347a0984de8f29feccce730866780fec5
plot(ent_vs_prop$prop~ent_vs_prop$V2)


##########################################################################
#####DESeq2 for PEPperPRINT data -- one time calculation for later use####
# 
# ###Ran through once for each IgG/IgA
# m=fread(paste0(WORKING_DIR, "PEPperPRINT_IgG.txt"))
# #m=fread(paste0(WORKING_DIR, "PEPperPRINT_IgA.txt"))
# samples = fread(paste0(WORKING_DIR, "PEPperPRINT_sample_info.txt"))
# 
# 
# samples$Group[1:6] = "Naive"
# m[,2:ncol(m)] = round(m[,2:ncol(m)])
# 
# library(DESeq2)
# gene_dat <- t(m)
# colnames(gene_dat) <- gene_dat[1,]
# gene_dat <- gene_dat[-1,]
# 
# countTable <- gene_dat
# 
# 
# countTable = data.frame(countTable, check.names = F)
# row.names(countTable) <- row.names(gene_dat)
# 
# sample_info <- samples
# sample_info$Group = as.factor(sample_info$Group)
# countTable = countTable[,order(colnames(countTable))]
# countTable = apply(countTable, 2, function(x){as.integer(x)})
# # Mike Love was the source on doing a paired comparison: https://support.bioconductor.org/p/58893/
# dds = DESeqDataSetFromMatrix(countData = countTable, 
#                              colData = sample_info, 
#                              design =  ~ Group)
# 
# dds = DESeq(dds)
# rnms <- resultsNames(dds)
# res = as.data.frame(results(dds, independentFiltering = FALSE, cooksCutoff = Inf))
# row.names(res) <- row.names(gene_dat)
# 
# #res_IgA=res
# #res_IgG=res
# fwrite(res_IgA, paste0(WORKING_DIR, "PEPperPRINT_IgA_DESeq2.txt"))  
# fwrite(res_IgG, paste0(WORKING_DIR, "PEPperPRINT_IgG_DESeq2.txt"))  

########################################################
####Table 1 -- T cell epitope selection criteria########

#HLA-I data
HLAI_filt = fread(paste0(WORKING_DIR, "COVID_human_netMHCpan_rank_filtered.txt"))
#HLA-II data
HLAII_filt = fread(paste0(WORKING_DIR, "COVID_human_netMHCIIpan_rank_filtered.txt"))

#HLA-I/II coepitope data
joint_run_filt_entropy = fread(paste0(WORKING_DIR, "COVID_combined_entropy_rank_filtered.csv"))
joint_run_filt_entropy$c1_tot_freq = joint_run_filt_entropy$c1_tot_freq*100
joint_run_filt_entropy$c2_tot_freq = joint_run_filt_entropy$c2_tot_freq*100

#Adding in mean HLA-I/II pop. frequency
joint_run_filt_entropy$Mean_freq=NA
for(z in 1:nrow(joint_run_filt_entropy)){
  joint_run_filt_entropy$Mean_freq[z] =sum(as.numeric(joint_run_filt_entropy$c1_tot_freq[z]), as.numeric(joint_run_filt_entropy$c2_tot_freq[z]))/2
}

##Rank matrix by most frequent epitopes, filter by entropy
Tc_eps = joint_run_filt_entropy[rev(order(joint_run_filt_entropy$Mean_freq)),]
Tc_eps = Tc_eps[which(Tc_eps$lo_entropy ==1),]

#Create new matrix to write into; define the top 5 most common binders
Pep_set = matrix(nrow=0,ncol=7)
Top5 = unique(Tc_eps$c2_peptide)[1:5]

#For each epitope in top 5, grab the peptide sequence, protein, start site
for(z in 1:length(Top5)){
  print(z)
  Peptide = Top5[z]
  Protein = unique(Tc_eps$gene[which(Tc_eps$c2_peptide == Top5[z])])
  Start = unique(Tc_eps$c2_start[which(Tc_eps$c2_peptide == Top5[z])])
  HLA_I = unique(Tc_eps$c1_peptide[which(Tc_eps$c2_peptide == Top5[z])])
  HLA_II = Top5[z]
  
  #Grabs HLA-I haplotypes that bind each peptide (top 1st percentile)
  HLA_I_haps = c()
  for(q in HLA_I){
    HLAI_haps_temp = HLAI_filt[which(HLAI_filt$Peptide == q),seq(8,84,5),with=F]
    HLAI_ranks_temp = HLAI_filt[which(HLAI_filt$Peptide == q),seq(9,85,5),with=F]
    HLA_I_haps = c(HLA_I_haps,colnames(HLAI_haps_temp)[which(HLAI_ranks_temp <= 1 )])
  }
  HLA_I_haps = sort(unique(HLA_I_haps))
  
  #Grabs HLA-II haplotypes that bind each peptide (top 5th percentile)
  HLA_II_haps = c()
  for(q in HLA_II){
    HLAII_haps_temp = HLAII_filt[which(HLAII_filt$Peptide == q),c(seq(6,25,3),seq(31,53,3)),with=F]
    HLAII_ranks_temp = HLAII_filt[which(HLAII_filt$Peptide == q),c(seq(7,26,3),seq(32,54,3)),with=F]
    HLA_II_haps = c(HLA_II_haps,colnames(HLAII_haps_temp)[which(HLAII_ranks_temp <= 5 )])
  }
  HLA_II_haps = sort(unique(HLA_II_haps))
  
  #Writes values into matrix  
  Pep_set = rbind(Pep_set, c(Peptide, Protein, Start, 
                             paste0(HLA_I,collapse = "; "), paste0(HLA_II, collapse='; '), 
                             paste0(HLA_I_haps, collapse = "; "),paste0(HLA_II_haps, collapse = "; ")))
  
}

##Calculating population frequency for having 1+ haplotype able to bind each epitope
Freq_I = c()
Freq_II = c()
All_I = c()
All_II = c()

##Reading in genetic frequency values
Freq=fread(paste0(WORKING_DIR, "HLA_freq.txt"))
Freq$US[1:10] = Freq$Cauc_Am[1:10]

for (n in 1:nrow(Pep_set)){
  print (n)
  HLA_I_all = strsplit(Pep_set[n,6], "; ") #Create list of all binding haplotypes calculated above
  All_I = c(All_I, HLA_I_all) #Used for later calculation for seeing proportion response to at least 1 peptide
  probs = rep(Freq$US[which(Freq$Haplotype %in% HLA_I_all[[1]])],2)/100 #Calculate pop freq. of haplotypes
  Freq_I = c(Freq_I, 1-prod(1-probs)) #Calculate pop freq. of haplotypes
  
  HLA_II_all = strsplit(Pep_set[n,7], "; ")
  All_II = c(All_II, HLA_II_all)
  probs = rep(Freq$US[which(Freq$Haplotype %in% HLA_II_all[[1]])],2)/100
  Freq_II = c(Freq_II, 1-prod(1-probs))
}

Pep_set = cbind(Pep_set, Freq_I, Freq_II)

All_I = sort(unique(unlist(All_I)))
All_II = sort(unique(unlist(All_II)))

#Pop freq that will respond to 1+ HLA-I peptide
probs_I = rep(Freq$US[which(Freq$Haplotype %in% All_I)],2)/100
Freq_I_any = 1-prod(1-probs_I)

#Pop freq that will respond to 1+ HLA-II peptide
probs_II = rep(Freq$US[which(Freq$Haplotype %in% All_II)],2)/100
Freq_II_any = 1-prod(1-probs_II)

colnames(Pep_set) = c("Peptide", "Protein", "Start", "HLA-I epitopes", "HLA-II epitope", "HLA-I haplotypes", "HLA-II haplotypes",
                      "HLA-I pop. freq.", "HLA-II pop. freq.")

fwrite(Pep_set, paste0(WORKING_DIR, "Table_1.txt"), sep = '\t', col.names = T, row.names = F, quote = F)


###################################  
#####Figure 3A#####################

library(ggplot2)
library(ggrepel)
library(ggnewscale)
library(data.table)
library(GenomicRanges)

####Showing mouse nM on human graph

#Minimum frequency to display (mean of HLA-I and II frequencies)  
Min_frequency = 0
label_frequency = 65

mouse1 = fread(paste0(WORKING_DIR, "COVID_murine_NetMHCpan_unfilt.xls"))
mouse1 = mouse1[which(mouse1$NB>0),]
mouse2 = fread(paste0(WORKING_DIR, "COVID_murine_NetMHCIIpan_unfilt.xls"))
mouse2 = mouse2[which(mouse2$NB>0),]

#HLA-I data
HLAI_filt = fread(paste0(WORKING_DIR, "COVID_human_netMHCpan_rank_filtered.txt"))
#HLA-II data
HLAII_filt = fread(paste0(WORKING_DIR, "COVID_human_netMHCIIpan_rank_filtered.txt"))
#HLA-I/II coepitope data
joint_run_filt_entropy = fread(paste0(WORKING_DIR, "COVID_combined_entropy_rank_filtered.csv"))
joint_run_filt_entropy$c1_tot_freq = joint_run_filt_entropy$c1_tot_freq*100
joint_run_filt_entropy$c2_tot_freq = joint_run_filt_entropy$c2_tot_freq*100


##########Filter by entropy manually -- ran through once to calculat high/low entropy epitopes#########

#ent = fread(paste0(WORKING_DIR, "entropy.MT072688.1.corrected.csv")
#ent = fread(paste0(WORKING_DIR, "AA_mafft_df_entropy.txt"))
ent = fread(paste0(WORKING_DIR, "entropy_7882.txt"))

colnames(ent) = c("position", "entropy")

###For each HLA-I, HLA-II, and co-epitope sets, mark which epitopes contain entropy >0.1 in any residue
for(n in 1:nrow(joint_run_filt_entropy)){
  if( max(ent$entropy[joint_run_filt_entropy$c2_start[n]:(joint_run_filt_entropy$c2_start[n]+14)]) > 0.1){
    joint_run_filt_entropy$lo_entropy[n] = 0
  }else{
    joint_run_filt_entropy$lo_entropy[n] = 1
  }
}
#fwrite(joint_run_filt_entropy, paste0(WORKING_DIR, "COVID_combined_entropy_rank_filtered.csv"))

HLAI_filt$lo_entropy=NA
#HLAI_filt$Pos = HLAI_filt$Pos+1
for(n in 1:nrow(HLAI_filt)){
  if( max(ent$entropy[HLAI_filt$Pos[n]:(HLAI_filt$Pos[n]+nchar(HLAI_filt$Peptide[n])-1)]) > 0.1){
    HLAI_filt$lo_entropy[n] = 0
  }else{
    HLAI_filt$lo_entropy[n] = 1
  }
}
#fwrite(HLAI_filt, paste0(WORKING_DIR, "COVID_human_netMHCpan_rank_filtered.txt"))

HLAII_filt$lo_entropy=NA
for(n in 1:nrow(HLAII_filt)){
  if( max(ent$entropy[HLAII_filt$Pos[n]:(HLAII_filt$Pos[n]+nchar(HLAII_filt$Peptide[n])-1)]) > 0.1){
    HLAII_filt$lo_entropy[n] = 0
  }else{
    HLAII_filt$lo_entropy[n] = 1
  }
}
#fwrite(HLAII_filt,paste0(WORKING_DIR, "COVID_human_netMHCIIpan_rank_filtered.txt"))

############################

###Number of MHC-I/II epitopes conserved between human and mouse
length(which(HLAI_filt$Peptide[which(HLAI_filt$lo_entropy==1)] %in% unique(mouse1$Peptide)))
length(which(HLAII_filt$Peptide[which(HLAII_filt$lo_entropy==1)] %in% unique(mouse2$Peptide)))

##Adding murine coverage for each human epitope
joint_run_filt_entropy$murine = "None"
joint_run_filt_entropy$murine[which(joint_run_filt_entropy$c1_peptide %in% mouse1$Peptide)] = "MHC-I" 
joint_run_filt_entropy$murine[which(joint_run_filt_entropy$c2_peptide %in% mouse2$Peptide)] = "MHC-II" 
joint_run_filt_entropy$murine[which((joint_run_filt_entropy$c2_peptide %in% mouse2$Peptide) & (joint_run_filt_entropy$c1_peptide %in% mouse1$Peptide))] = "Both" 
joint_run_filt_entropy$murine = factor(joint_run_filt_entropy$murine, levels = c("MHC-I", "MHC-II", "Both", "None"))


###Adding HLA-I epitopes offset by n amino acids so I can plot them overlapped -- can ignore###
joint_run_filt_entropy$c1_peptide_offset=NA
joint_run_filt_entropy$Mean_freq=NA

for(z in 1:nrow(joint_run_filt_entropy)){
  #Uses gregexpr to find c1 peptide within c2 peptide --> inserts spaces ahead of c1 peptide so the start AA aligns with its location on c2 peptide --> Pastes out as a string in matrix
  joint_run_filt_entropy$c1_peptide_offset[z] = paste0(c(rep(" ", (gregexpr(joint_run_filt_entropy$c2_peptide[z], pattern = joint_run_filt_entropy$c1_peptide[z])[[1]][1])-1), joint_run_filt_entropy$c1_peptide[z]), collapse = '')
  joint_run_filt_entropy$Mean_freq[z] =sum(as.numeric(joint_run_filt_entropy$c1_tot_freq[z]), as.numeric(joint_run_filt_entropy$c2_tot_freq[z]))/2
}

##Filter by Min_frequency
joint_run_filt_entropy = joint_run_filt_entropy[which(joint_run_filt_entropy$Mean_freq >= Min_frequency),]

##Defining separate matrices for lo_entropy/hi_entropy (for plotting later)
lo_entropy = subset(joint_run_filt_entropy, lo_entropy==1)
hi_entropy = subset(joint_run_filt_entropy, lo_entropy==0)



#########Adding literature Tc data############
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


###Check overlaps###

length(which(lo_entropy$c1_peptide %in% unique(paper_tc$Peptide)))
<<<<<<< HEAD
length(which(HLAI_filt$Peptide[which(HLAI_filt$lo_entropy==1)] %in% unique(paper_tc$Peptide)))
length(which(HLAII_filt$Peptide[which(HLAII_filt$lo_entropy==1)] %in% unique(paper_tc$Peptide)))
=======
which(HLAI_filt$Peptide[which(HLAI_filt$lo_entropy==1)] %in% unique(paper_tc$Peptide))
which(HLAII_filt$Peptide[which(HLAII_filt$lo_entropy==1)] %in% unique(paper_tc$Peptide))
>>>>>>> 8641edb347a0984de8f29feccce730866780fec5



####Adding in Bc data from literature###
#Paper_bc = fread(paste0(WORKING_DIR, "Paper_Bc_epitopes_complete.txt")
Paper_bc = fread(paste0(WORKING_DIR, "Supplemental/B_cell_epitopes - Bc_epitopes.txt"))
seq =  fread(paste0(WORKING_DIR, "AA_sequence.txt"), header=F)

###Grabbing S, M, and N sequences -- only ones present in current dataset####
S= seq$V1[4]
M = seq$V1[10]
N=seq$V1[18]

Paper_bc$Start = NA
Paper_bc$End = NA

####For each peptide, looking to see if it matches to SARS-CoV-2 reference###
for(p in 1:nrow(Paper_bc)){
  if(is.na(Paper_bc$Protein[p])){  #If "Protein" is missing, look for it in S, N, or M
    if (gregexpr(pattern = Paper_bc$Epitope[p],S)[[1]][1] != -1){
      Paper_bc$Protein[p] = "S"
    }else if (gregexpr(pattern = Paper_bc$Epitope[p],M)[[1]][1] != -1){
      Paper_bc$Protein[p] = "M"
    }else if (gregexpr(pattern = Paper_bc$Epitope[p],N)[[1]][1] != -1){
      Paper_bc$Protein[p] = "N"
    }
  }
  if(!is.na(Paper_bc$Protein[p])){  #If protein is already defined, look for sequence identity in SARS-CoV-2 and grab coordinates
    if(Paper_bc$Protein[p] == "S"){
      if (gregexpr(pattern = Paper_bc$Epitope[p],S)[[1]][1] != -1){
        Paper_bc$Start[p] = 7096+gregexpr(pattern = Paper_bc$Epitope[p],S)[[1]][1]
        Paper_bc$End[p] = Paper_bc$Start[p] + nchar(Paper_bc$Epitope[p]) -1
      }
    }else if(Paper_bc$Protein[p] == "M"){
      if (gregexpr(pattern = Paper_bc$Epitope[p],M)[[1]][1] != -1){
        Paper_bc$Start[p] = 8719+gregexpr(pattern = Paper_bc$Epitope[p],M)[[1]][1]
        Paper_bc$End[p] = Paper_bc$Start[p] + nchar(Paper_bc$Epitope[p]) -1
      }
      
    }else if(Paper_bc$Protein[p] == "N"){
      if (gregexpr(pattern = Paper_bc$Epitope[p],N)[[1]][1] != -1){
        Paper_bc$Start[p] = 9244+gregexpr(pattern = Paper_bc$Epitope[p],N)[[1]][1]
        Paper_bc$End[p] = Paper_bc$Start[p] + nchar(Paper_bc$Epitope[p]) -1
      }
    }
  }  
}
Paper_bc$Protein = factor(Paper_bc$Protein, levels = c("S","M", "N"))

#Remove values not found
Paper_bc = Paper_bc[which(!is.na(Paper_bc$Start)),]

##Split into epitopes found from the literature vs from PEPperCHIP
<<<<<<< HEAD
PPC_bc = Paper_bc[which((Paper_bc$Source == "PEPperCHIP"))] #| (Paper_bc$Source == "PEPperCHIP_etc")),]
=======
PPC_bc = Paper_bc[which((Paper_bc$Source == "PEPperCHIP") | (Paper_bc$Source == "PEPperCHIP_etc")),]
>>>>>>> 8641edb347a0984de8f29feccce730866780fec5
Paper_bc = Paper_bc[which((Paper_bc$Source != "PEPperCHIP") & (Paper_bc$Source != "PEPperCHIP_etc")),]
Paper_bc$Type = "B cell"

colnames(Paper_bc)[3] = "Peptide"
colnames(PPC_bc)[3] = "Peptide"

<<<<<<< HEAD

###Overlap between PPC and lit
which(PPC_bc$Peptide %in% Paper_bc$Peptide) #Exact hits

#Overlaps
g1 <-  GRanges(seqnames="COVID",
               IRanges(start=Paper_bc$Start,
                       end=Paper_bc$End),
               Lit_Bc = Paper_bc$Peptide,
               Lit_source = Paper_bc$Source)
              

g2 <-  GRanges(seqnames="COVID",
               IRanges(start=PPC_bc$Start,
                       end=PPC_bc$End),
               Pep_Bc = PPC_bc$Peptide)

merge_bc=as.data.table(mergeByOverlaps(g1,g2))



=======
>>>>>>> 8641edb347a0984de8f29feccce730866780fec5
###Combine T and B cell data
paper_all = rbind(Paper_bc[,c(3,4,14:16)],paper_tc[,c(11,1,4,5,12)])



####Check overlap between Bc lit epitopes with HLA###########


#Change to PPC_bc for looking for PEPperCHIP overlaps
<<<<<<< HEAD
#Source = PPC_bc[which(PPC_bc$Source == "PEPperCHIP"),]
Source = Paper_bc
=======
Source = PPC_bc[which(PPC_bc$Source == "PEPperCHIP"),]
#Source = Paper_bc
>>>>>>> 8641edb347a0984de8f29feccce730866780fec5

##Coepitopes
g1 <-  GRanges(seqnames="COVID",
               IRanges(start=lo_entropy$c2_start,
                       end=(lo_entropy$c2_start + nchar(lo_entropy$c2_peptide) -1)),
               Pep_Tc_1 = lo_entropy$c1_peptide,
               Pep_Tc_2 = lo_entropy$c2_peptide)

g2 <-  GRanges(seqnames="COVID",
               IRanges(start=Source$Start,
                       end=Source$End),
               Pep_Bc = Source$Peptide,
               Protein = Source$Protein,
               Lit_source = Source$Source)

merge_coep=as.data.table(mergeByOverlaps(g1,g2))

nrow(merge_coep) #Number of total overlaps against co-epitope set
length(which(merge_coep$Protein=="S")) 
length(which(merge_coep$Protein=="M"))
length(which(merge_coep$Protein=="N"))

length(unique(merge_coep$Pep_Tc_1)) #Number of overlapping HLA-I epitopes
length(unique(merge_coep$Pep_Tc_2)) #Number of overlapping HLA-II epitopes
length(unique(merge_coep$Pep_Bc)) #Number of overlapping Bc epitopes

#HLA-I
######Transforming to proteomic space####
HLAI_filt_proteome = HLAI_filt
HLAI_filt_proteome$Pos[which(HLAI_filt_proteome$ID == "S_surface_glyco")] = 7096 + HLAI_filt_proteome$Pos[which(HLAI_filt_proteome$ID == "S_surface_glyco")]
HLAI_filt_proteome$Pos[which(HLAI_filt_proteome$ID == "ORF3a_ORF3a_pro")] = 8369 + HLAI_filt_proteome$Pos[which(HLAI_filt_proteome$ID == "ORF3a_ORF3a_pro")]
HLAI_filt_proteome$Pos[which(HLAI_filt_proteome$ID == "E_envelope_prot" )] = 8644 + HLAI_filt_proteome$Pos[which(HLAI_filt_proteome$ID == "E_envelope_prot" )] 
HLAI_filt_proteome$Pos[which(HLAI_filt_proteome$ID == "M_membrane_glyc" )] = 8719 + HLAI_filt_proteome$Pos[which(HLAI_filt_proteome$ID == "M_membrane_glyc" )]
HLAI_filt_proteome$Pos[which(HLAI_filt_proteome$ID == "ORF6_ORF6_prote")] = 8941 + HLAI_filt_proteome$Pos[which(HLAI_filt_proteome$ID == "ORF6_ORF6_prote")] 
HLAI_filt_proteome$Pos[which(HLAI_filt_proteome$ID == "ORF7a_ORF7a_pro")] = 9002 + HLAI_filt_proteome$Pos[which(HLAI_filt_proteome$ID == "ORF7a_ORF7a_pro")]
HLAI_filt_proteome$Pos[which(HLAI_filt_proteome$ID == "ORF8_ORF8_prote")] = 9123 + HLAI_filt_proteome$Pos[which(HLAI_filt_proteome$ID == "ORF8_ORF8_prote")]
HLAI_filt_proteome$Pos[which(HLAI_filt_proteome$ID == "N_nucleocapsid_")] = 9233 + HLAI_filt_proteome$Pos[which(HLAI_filt_proteome$ID == "N_nucleocapsid_")]
HLAI_filt_proteome$Pos[which(HLAI_filt_proteome$ID == "ORF10_ORF10_pro" )] = 9663 + HLAI_filt_proteome$Pos[which(HLAI_filt_proteome$ID == "ORF10_ORF10_pro" )]

g1 <-  GRanges(seqnames="COVID",
               IRanges(start=HLAI_filt_proteome$Pos,
                       end=(HLAI_filt_proteome$Pos + nchar(HLAI_filt_proteome$Peptide) -1)),
               Pep_Tc = HLAI_filt_proteome$Peptide,
               lo_entropy = HLAI_filt_proteome$lo_entropy)

g2 <-  GRanges(seqnames="COVID",
               IRanges(start=Source$Start,
                       end=Source$End),
               Pep_Bc = Source$Peptide,
               Protein = Source$Protein,
               Lit_source = Source$Source)

merge_I=as.data.table(mergeByOverlaps(g1,g2))
merge_I=merge_I[which(merge_I$lo_entropy==1),] #Filter by entropy

nrow(merge_I) #Number of total overlaps against HLA-I set
length(which(merge_I$Protein=="S"))
length(which(merge_I$Protein=="M"))
length(which(merge_I$Protein=="N"))

length(unique(merge_I$Pep_Tc)) #Number of HLA-I epitopes
length(unique(merge_I$Pep_Bc)) #Number of Bc epitopes

#HLA-II
######Transforming to proteomic space -- not all proteins are present####
HLAII_filt_proteome = HLAII_filt
HLAII_filt_proteome$Pos[which(HLAII_filt_proteome$ID == "S_surface_glycoprotein_QIB84673.1")] = 7096 + HLAII_filt_proteome$Pos[which(HLAII_filt_proteome$ID == "S_surface_glycoprotein_QIB84673.1")]
HLAII_filt_proteome$Pos[which(HLAII_filt_proteome$ID == "ORF3a_ORF3a_protein_QIB84674.1" )] = 8369 + HLAII_filt_proteome$Pos[which(HLAII_filt_proteome$ID == "ORF3a_ORF3a_protein_QIB84674.1" )] 
#HLAII_filt_proteome$Pos[which(HLAII_filt_proteome$ID == "E_envelope_prot" )] = 8344 + HLAII_filt_proteome$Pos[which(HLAII_filt_proteome$ID == "E_envelope_prot" )] 
HLAII_filt_proteome$Pos[which(HLAII_filt_proteome$ID == "M_membrane_glycoprotein_QIB84676.1" )] = 8719 +HLAII_filt_proteome$Pos[which(HLAII_filt_proteome$ID == "M_membrane_glycoprotein_QIB84676.1" )]
HLAII_filt_proteome$Pos[which(HLAII_filt_proteome$ID == "ORF6_ORF6_protein_QIB84677.1" )] = 8941 + HLAII_filt_proteome$Pos[which(HLAII_filt_proteome$ID == "ORF6_ORF6_protein_QIB84677.1" )]
HLAII_filt_proteome$Pos[which(HLAII_filt_proteome$ID == "ORF7a_ORF7a_protein_QIB84678.1")] = 9002 + HLAII_filt_proteome$Pos[which(HLAII_filt_proteome$ID == "ORF7a_ORF7a_protein_QIB84678.1")] 
HLAII_filt_proteome$Pos[which(HLAII_filt_proteome$ID == "ORF8_ORF8_protein_QIB84679.1")] = 9123 + HLAII_filt_proteome$Pos[which(HLAII_filt_proteome$ID == "ORF8_ORF8_protein_QIB84679.1")] 
HLAII_filt_proteome$Pos[which(HLAII_filt_proteome$ID == "N_nucleocapsid_phosphoprotein_QIB84680.1")] = 9233 +HLAII_filt_proteome$Pos[which(HLAII_filt_proteome$ID == "N_nucleocapsid_phosphoprotein_QIB84680.1")]
#HLAII_filt_proteome$Pos[which(HLAII_filt_proteome$ID == "ORF10_ORF10_pro" )] = 9663 + HLAII_filt_proteome$Pos[which(HLAII_filt_proteome$ID == "ORF10_ORF10_pro" )]

g1 <-  GRanges(seqnames="COVID",
               IRanges(start=HLAII_filt_proteome$Pos,
                       end=(HLAII_filt_proteome$Pos + nchar(HLAII_filt_proteome$Peptide) -1)),
               Pep_Tc = HLAII_filt_proteome$Peptide,
               lo_entropy = HLAII_filt_proteome$lo_entropy)

g2 <-  GRanges(seqnames="COVID",
               IRanges(start=Source$Start,
                       end=Source$End),
               Pep_Bc = Source$Peptide,
               Protein = Source$Protein,
               Lit_source = Source$Source)
merge_II=as.data.table(mergeByOverlaps(g1,g2))
merge_II=merge_II[which(merge_II$lo_entropy==1),]

nrow(merge_II) #Number of total overlaps against HLA-II set
length(which(merge_II$Protein=="S"))
length(which(merge_II$Protein=="M"))
length(which(merge_II$Protein=="N"))

length(unique(merge_II$Pep_Tc)) #Number of HLA-II epitopes
length(unique(merge_II$Pep_Bc)) #Number of Bc epitopes

# ###Check in co-epitope set
# length(which(lo_entropy$c2_peptide %in% merge_II$Pep_Tc))
# length(which(lo_entropy$c1_peptide %in% merge_I$Pep_Tc))
# length(which((lo_entropy$c2_peptide %in% merge_II$Pep_Tc)&(lo_entropy$c1_peptide %in% merge_I$Pep_Tc)))
# 
# length(unique(c(which(lo_entropy$c2_peptide %in% merge_II$Pep_Tc),which(lo_entropy$c1_peptide %in% merge_I$Pep_Tc))))


#Add entropy#
#ent = fread(paste0(WORKING_DIR, "entropy.MT072688.1.corrected.csv"))
#ent = fread(paste0(WORKING_DIR, "AA_mafft_df_entropy.txt"))
ent = fread(paste0(WORKING_DIR, "entropy_7882.txt"))

###Set y lower limit for graphing###
y_ll = min(joint_run_filt_entropy$c1_tot_freq)


####Plotting Figure 3A#####
ggplot(data=joint_run_filt_entropy) + 
  geom_segment(data=hi_entropy, aes(x=c1_start, xend=c1_start, y=y_ll, yend=c1_tot_freq), color="black", alpha=0.25) + 
  geom_segment(data=lo_entropy, aes(x=c1_start, xend=c1_start, y=y_ll, yend=c1_tot_freq), color="orange", alpha=0.5, size=1.1) +
  geom_point(aes(x=c1_start, y=c1_tot_freq, size=(c2_tot_freq), color = murine), shape = 1, stroke = 2)+
  scale_colour_manual(values = c("red", "blue", "purple", "black"), name="Murine MHC\nco-expression")+
  ggnewscale::new_scale("color")+
  
  geom_point(aes(x=c1_start, y=c1_tot_freq, color=c1_min_nm, size=c2_tot_freq), alpha=.5) +
  
  scale_color_viridis_c(direction = -1,  trans = "log", name = "HLA-I nM") + 
  theme_light() +
  scale_size(range = c(0,10), breaks = c(10,30,50,70,90))+
  scale_y_continuous(limits = c((y_ll-20),95))+
  #scale_x_continuous(limits = c(7096,9662))+
  
  labs( size = "HLA-II pop. frequency (%)") +
  #geom_text(aes(x=c1_start, y=c1_tot_freq,label=ifelse(Mean_freq>label_frequency,as.character(c2_peptide),'')),angle=90,color="red",hjust=0, vjust=0, family = "mono", size = 5)+
  #geom_text(aes(x=c1_start, y=c1_tot_freq,label=ifelse(Mean_freq>label_frequency,as.character(c1_peptide_offset),'')),angle=90, color="black",hjust=0, vjust=0, family = "mono", size =5, fontface="bold")+
  geom_label_repel(aes(label=ifelse(Mean_freq>label_frequency,as.character(c1_peptide),''), x=c1_start, y=c1_tot_freq),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  
  geom_rect(aes(xmin=0, xmax=7095, ymin=(y_ll-4), ymax=(y_ll-8)), fill="white", color="black", size=0.1) +
  geom_text(aes(x=7095/2, y=(y_ll-6), label="orf1ab", angle=0), size=2.5) +
  
  geom_rect(aes(xmin=7096, xmax=8368, ymin=(y_ll-8), ymax=(y_ll-12)), fill="white", color="black", size=0.1) +
  geom_text(aes(x=7096+(8368-7096)/2, y=(y_ll-10), label="S", angle=0), size=2.5) + 
  
  geom_rect(aes(xmin=8369, xmax=8644, ymin=(y_ll-4), ymax=(y_ll-8)), fill="white", color="black", size=0.1) +
  geom_text(aes(x=8368+(8644-8368)/2, y=(y_ll-6), label="ORF3a", angle=0), size=2.5) +
  
  geom_rect(aes(xmin=8645, xmax=8718, ymin=(y_ll-8), ymax=(y_ll-12)), fill="white", color="black", size=0.1) +
  geom_text(aes(x=8645+(8718-8645)/2, y=(y_ll-10)), label="E", angle=0, size=2.5) + 
  
  geom_rect(aes(xmin=8719, xmax=8940, ymin=(y_ll-4), ymax=(y_ll-8)), fill="white", color="black", size=0.1) +
  geom_text(aes(x=8719+(8940-8719)/2, y=(y_ll-6)), label="M", angle=0, size=2.5) +
  
  geom_rect(aes(xmin=8941, xmax=9001, ymin=(y_ll-8), ymax=(y_ll-12)), fill="white", color="black", size=0.1) +
  geom_text(aes(x=8941+(9001-8941)/2, y=(y_ll-10)), label="ORF6", angle=0, size=2.5) +
  
  geom_rect(aes(xmin=9002, xmax=9122, ymin=(y_ll-4), ymax=(y_ll-8)), fill="white", color="black", size=0.1) +
  geom_text(aes(x=9002+(9122-9002)/2, y=(y_ll-6)), label="ORF7a", angle=0, size=2.5) +
  
  geom_rect(aes(xmin=9123, xmax=9243, ymin=(y_ll-8), ymax=(y_ll-12)), fill="white", color="black", size=0.1) +
  geom_text(aes(x=9123+(9243-9123)/2, y=(y_ll-10)), label="ORF8", angle=0, size=2.5) +
  
  geom_rect(aes(xmin=9244, xmax=9662, ymin=(y_ll-4), ymax=(y_ll-8)), fill="white", color="black", size=0.1) +
  geom_text(aes(x=9244+(9662-9244)/2, y=(y_ll-6)), label="N", angle=0, size=2.5) +
  
  geom_rect(aes(xmin=9663, xmax=9701, ymin=(y_ll-8), ymax=(y_ll-12)), fill="white", color="black", size=0.1) +
  geom_text(aes(x=9663+(9701-9663)/2, y=(y_ll-10)), label="ORF10", angle=0, size=2.5) + 
  
  #Literature track
  geom_rect(data = paper_all,aes(ymin=y_ll-16, ymax=y_ll-20, xmin=Start, xmax=End, fill = Type),size=0)+
  scale_fill_manual(values = alpha(c("blue", "red"), .5))+
  geom_text(aes(x = 9670, y = y_ll-18, label = "Literature"), hjust = 0, size = 3.5)+
  
  #PEPperPRINT track
  geom_rect(data = PPC_bc,aes(ymin=y_ll-12, ymax=y_ll-16, xmin=Start, xmax=End), fill = alpha(c("blue"), .5),size=0)+
  geom_text(aes(x = 9670, y = y_ll-14, label = "PEPperCHIP"), hjust = 0, size = 3.5)+
  
  theme(panel.background = element_blank(), panel.grid.major.x = element_line(color="azure3"), panel.grid.major.y = element_line(color = "azure3"), axis.ticks = element_blank()) +
  labs(y="HLA-I pop. frequency (%)", x="Position across SARS-CoV-2 proteome") +
  theme(text=element_text(face="bold",size=20,colour="black")) 
#+scale_y_continuous(sec.axis = sec_axis(~., name = "Average Entropy"))



##########################################
###Figure 3B, Figure 3C###################
##########################################
library(GenomicRanges)


Min_frequency = 0  #Min pop. freq to plot (currently set to plot all)
label_frequency = 65 #Min pop. freq. to label with text

#Read in co-epitope dataset
joint_run_filt_entropy = fread(paste0(WORKING_DIR, "COVID_combined_entropy_rank_filtered.csv"))
joint_run_filt_entropy$c1_tot_freq = joint_run_filt_entropy$c1_tot_freq*100
joint_run_filt_entropy$c2_tot_freq = joint_run_filt_entropy$c2_tot_freq*100


####Adding some metrics##########
joint_run_filt_entropy$c1_peptide_offset=NA  #Can ignore. For plotting HLA-I labels overlayed on HLA-II epitopes
joint_run_filt_entropy$Mean_freq=NA #Mean HLA-I/II population frequency
joint_run_filt_entropy$Mean_nM = NA #Mean HLA-I/II binding affinity

for(z in 1:nrow(joint_run_filt_entropy)){
  #Uses gregexpr to find c1 peptide within c2 peptide --> inserts spaces ahead of c1 peptide so the start AA aligns with its location on c2 peptide --> Pastes out as a string in matrix
  joint_run_filt_entropy$c1_peptide_offset[z] = paste0(c(rep(" ", (gregexpr(joint_run_filt_entropy$c2_peptide[z], pattern = joint_run_filt_entropy$c1_peptide[z])[[1]][1])-1), joint_run_filt_entropy$c1_peptide[z]), collapse = '')
  joint_run_filt_entropy$Mean_freq[z] =sum(as.numeric(joint_run_filt_entropy$c1_tot_freq[z]), as.numeric(joint_run_filt_entropy$c2_tot_freq[z]))/2
  joint_run_filt_entropy$Mean_freq[z] =sum(as.numeric(joint_run_filt_entropy$c1_tot_freq[z]), as.numeric(joint_run_filt_entropy$c2_tot_freq[z]))/2
  joint_run_filt_entropy$Mean_nM[z] =sum(as.numeric(joint_run_filt_entropy$c1_min_nm[z]), as.numeric(joint_run_filt_entropy$c2_min_nm[z]))/2
  
}

#########Adding literature Tc data -- same as above############
paper_tc = fread(paste0(WORKING_DIR, "Supplemental/Standardized T cell epitopes.txt"))
paper_tc = paper_tc[!duplicated(paper_tc),]

joint_run_filt_entropy$Tc_lit_coverage=NA
for (u in 1:nrow(paper_tc)){
  if (length(which(str_count(joint_run_filt_entropy$c2_peptide, paper_tc$Peptide[u])==1)) > 0){
    joint_run_filt_entropy$Tc_lit_coverage[which(str_count(joint_run_filt_entropy$c2_peptide, paper_tc$Peptide[u])==1)] = paper_tc$Source[u]
  }
}

paper_tc = paper_tc[which(!is.na(paper_tc$Start)),]
paper_tc$Type = "T cell"


####Adding in Bc data from literature -- same as above###
#Paper_bc = fread(paste0(WORKING_DIR, "Paper_Bc_epitopes_complete.txt"))
Paper_bc = fread(paste0(WORKING_DIR, "Supplemental/B_cell_epitopes - Bc_epitopes.txt"))

seq =  fread(paste0(WORKING_DIR, "AA_sequence.txt"), header = F)
S= seq$V1[4]
M = seq$V1[10]
N= seq$V1[18]

Paper_bc$Start = NA
Paper_bc$End = NA
for(p in 1:nrow(Paper_bc)){
  if(is.na(Paper_bc$Protein[p])){
    if (gregexpr(pattern = Paper_bc$Epitope[p],S)[[1]][1] != -1){
      Paper_bc$Protein[p] = "S"
    }else if (gregexpr(pattern = Paper_bc$Epitope[p],M)[[1]][1] != -1){
      Paper_bc$Protein[p] = "M"
    }else if (gregexpr(pattern = Paper_bc$Epitope[p],N)[[1]][1] != -1){
      Paper_bc$Protein[p] = "N"
    }
  }
  if(!is.na(Paper_bc$Protein[p])){
    if(Paper_bc$Protein[p] == "S"){
      if (gregexpr(pattern = Paper_bc$Epitope[p],S)[[1]][1] != -1){
        Paper_bc$Start[p] = 7096+gregexpr(pattern = Paper_bc$Epitope[p],S)[[1]][1]
        Paper_bc$End[p] = Paper_bc$Start[p] + nchar(Paper_bc$Epitope[p]) -1
      }
    }else if(Paper_bc$Protein[p] == "M"){
      if (gregexpr(pattern = Paper_bc$Epitope[p],M)[[1]][1] != -1){
        Paper_bc$Start[p] = 8719+gregexpr(pattern = Paper_bc$Epitope[p],M)[[1]][1]
        Paper_bc$End[p] = Paper_bc$Start[p] + nchar(Paper_bc$Epitope[p]) -1
      }
      
    }else if(Paper_bc$Protein[p] == "N"){
      if (gregexpr(pattern = Paper_bc$Epitope[p],N)[[1]][1] != -1){
        Paper_bc$Start[p] = 9244+gregexpr(pattern = Paper_bc$Epitope[p],N)[[1]][1]
        Paper_bc$End[p] = Paper_bc$Start[p] + nchar(Paper_bc$Epitope[p]) -1
      }
    }
  }  
}
Paper_bc$Protein = factor(Paper_bc$Protein, levels = c("S","M", "N"))
Paper_bc = Paper_bc[which(!is.na(Paper_bc$Start)),] #Remove values without identity in SARS-CoV-2

#######Look for overlap between predicted Tc epitopes and Bc epitopes

##Overlap of Bc epitopes against co-epitopes
lo_entropy=joint_run_filt_entropy[which(joint_run_filt_entropy$lo_entropy==1),]

g1 <-  GRanges(seqnames="COVID",
               IRanges(start=lo_entropy$c2_start,
                       end=(lo_entropy$c2_start + nchar(lo_entropy$c2_peptide) -1)),
               Pep_Tc = lo_entropy$c2_peptide)

g2 <-  GRanges(seqnames="COVID",
               IRanges(start=Paper_bc$Start,
                       end=Paper_bc$End),
               Pep_Bc = Paper_bc$Epitope,
               Protein = Paper_bc$Protein,
               Lit_source = Paper_bc$Source)

merge=as.data.table(mergeByOverlaps(g1,g2))


#####Adding Bc overlap information into Tc data.table
joint_run_filt_entropy$Bc_coverage = "None"  #Whether or not there is Bc epitope overlap
joint_run_filt_entropy$Bc_epitope = NA  #If present, the overlapping Bc epitope(s)
joint_run_filt_entropy$Bc_lit_source = NA #The source for the epitope (literature, PEPperPRINT, etc)

for(i in 1:nrow(joint_run_filt_entropy)){
  
  if(joint_run_filt_entropy$c2_peptide[i] %in% merge$Pep_Tc){
    print(i)
    joint_run_filt_entropy$Bc_coverage[i] = unique(as.character(merge$Protein[which(merge$Pep_Tc == joint_run_filt_entropy$c2_peptide[i])]))
    joint_run_filt_entropy$Bc_epitope[i] = paste0(unique(as.character(merge$Pep_Bc[which(merge$Pep_Tc == joint_run_filt_entropy$c2_peptide[i])])), collapse = ";")
    joint_run_filt_entropy$Bc_lit_source[i] = paste0(unique(as.character(merge$Lit_source[which(merge$Pep_Tc == joint_run_filt_entropy$c2_peptide[i])])),collapse = ";")
    
  }
  
}
joint_run_filt_entropy$Bc_coverage = factor(joint_run_filt_entropy$Bc_coverage, levels = c("S", "M", "N", "None"))

#################filter data by frequency and entropy####################

joint_run_filt_entropy = joint_run_filt_entropy[which(joint_run_filt_entropy$Mean_freq >= Min_frequency),]
joint_run_filt_entropy = joint_run_filt_entropy[which(joint_run_filt_entropy$lo_entropy == 1),]


####Setting lower limits
y_ll = min(joint_run_filt_entropy$c2_tot_freq)
x_ll = min(joint_run_filt_entropy$c1_tot_freq)


#####Separate out rows which contain multiple Bc epitopes, adding fold-change values from above DESeq2#####
res_IgA = fread(paste0(WORKING_DIR, "PEPperPRINT_IgA_DESeq2.txt"))
rownames(res_IgA) = res_IgA$V1
res_IgG = fread(paste0(WORKING_DIR, "PEPperPRINT_IgG_DESeq2.txt"))
rownames(res_IgG) = res_IgG$V1

bc_mat = joint_run_filt_entropy[which(joint_run_filt_entropy$Bc_coverage != "None"),]
bc_mat$Bc_IgG_log2FC = NA
bc_mat$Bc_IgA_log2FC = NA

bc_sub = data.table()
for(t in 1:nrow(bc_mat)){
  eps = unlist(strsplit(bc_mat$Bc_epitope[t],";"))  #Split all overlapping Bc epitopes
  temp_dt = bc_mat[t][rep(1,length(eps))]  #Make data table with same format as bc_mat, but with a row per each Bc epitope
  temp_dt$Bc_epitope = eps  #Write out epitopes, one per row
  for(q in 1:nrow(temp_dt)){
    temp_dt$Bc_lit_source[q] = paste0(unique(merge$Lit_source[which(merge$Pep_Bc == temp_dt$Bc_epitope[q])]),collapse=";") #Fill in source of epitope
    if(temp_dt$Bc_epitope[q] %in% c(rownames(res_IgA), rownames(res_IgG))){  #If the epitope is present in the PEPperCHIP analysis (DESeq from above)
      temp_dt$Bc_IgG_log2FC[q] = res_IgG$log2FoldChange[which(rownames(res_IgG) == temp_dt$Bc_epitope[q])] #Write in IgG/IgA log2 fold change values
      temp_dt$Bc_IgA_log2FC[q] = res_IgA$log2FoldChange[which(rownames(res_IgA) == temp_dt$Bc_epitope[q])]
    }
  }
  bc_sub = rbind(bc_sub, temp_dt)
}


#Filter by entropy
bc_sub = bc_sub[which(bc_sub$lo_entropy == 1),]

#Keep only values with DESeq info, i.e. only PEPperCHIP data
bc_sub = bc_sub[which(!is.na(bc_sub$Bc_IgG_log2FC) | !is.na(bc_sub$Bc_IgA_log2FC)),]


################Plotting Figure 3B###################
pep_col = viridis(3, begin = 0, end = .75, option = "plasma")

ggplot(data=joint_run_filt_entropy) + 
  geom_point(data = bc_sub,aes(x=c1_tot_freq, y=c2_tot_freq, color=Bc_coverage, size=Mean_nM), shape = 1, stroke = 3) +
  scale_colour_manual(values =pep_col, name="B cell epitope\ncoverage")+
  ggnewscale::new_scale("color")+
  
  geom_point(aes(x=c1_tot_freq, y=c2_tot_freq, color=gene, size=Mean_nM), alpha=.5) +
  theme(text=element_text(face="bold",size=20,colour="black")) +
  scale_color_viridis_d()+
  #geom_text(aes(x=c1_tot_freq, y=c2_tot_freq,label=ifelse(Mean_freq>label_frequency,as.character(c2_peptide),'')),color="red",hjust=0, vjust=0, family = "mono", size = 5, angle = -23)+
  #geom_text(aes(x=c1_tot_freq, y=c2_tot_freq,label=ifelse(Mean_freq>label_frequency,as.character(c1_peptide_offset),'')), color="black",hjust=0, vjust=0, family = "mono", size =5, fontface="bold", angle = -23)+
  geom_label_repel(aes(label=ifelse(Mean_freq>label_frequency,as.character(c1_peptide),''), x=c1_tot_freq, y=c2_tot_freq),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  #scale_color_continuous(trans = "reverse", name = "Mean MHC nM") + 
  scale_y_continuous(limits = c(y_ll,95))+
  scale_x_continuous(limits = c(x_ll,95))+
  scale_size(range = c(0,10), breaks = c(10,150,300), trans="reverse", name = "Mean nM")+
  labs(y="MHC II Frequency (%)", x="MHC I Frequency (%)", title = "Predicted T cell epitopes", color="Gene") +
  ggnewscale::new_scale("color")



#############Plotting Figure 3C############
ggplot(data=bc_sub) + 
  geom_point(data = bc_sub,aes(x=c1_tot_freq, y=c2_tot_freq, color=Bc_coverage, size=Bc_IgG_log2FC), shape = 1, stroke = 2) +
  scale_colour_manual(values = pep_col, name="Bc epitope coverage")+
  ggnewscale::new_scale("color")+
  
  geom_point(aes(x=c1_tot_freq, y=c2_tot_freq, color=Bc_IgA_log2FC, size=Bc_IgG_log2FC), alpha=.5) +
  theme(text=element_text(face="bold",size=20,colour="black")) +
  #geom_text(aes(x=c1_tot_freq, y=c2_tot_freq,label=ifelse(Mean_freq>label_frequency,as.character(c2_peptide),'')),color="red",hjust=0, vjust=0, family = "mono", size = 5, angle = -23)+
  #geom_text(aes(x=c1_tot_freq, y=c2_tot_freq,label=ifelse(Mean_freq>label_frequency,as.character(c1_peptide_offset),'')), color="black",hjust=0, vjust=0, family = "mono", size =5, fontface="bold", angle = -23)+
  geom_label_repel(aes(label=as.character(Bc_epitope), x=c1_tot_freq, y=c2_tot_freq),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  scale_color_viridis_c(direction = 1, name = "IgA log(fold change) fluor.") + 
  scale_y_continuous(limits = c(y_ll,95))+
  scale_x_continuous(limits = c(x_ll,95))+
  scale_size(range = c(3,15), breaks = c(-.3,0,.5,1), name = "IgG log(fold change) fluor.")+
  labs(y="MHC II Frequency (%)", x="MHC I Frequency (%)", title = "Predicted B cell epitopes", color="Gene") +
  ggnewscale::new_scale("color")



#############################################################
#########Visual abstract#####################################

###A lot of stupid code to make a few color bars -- can probably ignore this####
joint_run_filt_entropy = fread(paste0(WORKING_DIR, "COVID_combined_entropy_rank_filtered.csv"))
lo_entropy = joint_run_filt_entropy[which(joint_run_filt_entropy$lo_entropy==1),]

I_filt=fread(paste0(WORKING_DIR, "COVID_human_netMHCpan_rank_filtered.txt"))
II_filt=fread(paste0(WORKING_DIR, "COVID_human_netMHCIIpan_rank_filtered.txt"))
joint_run_filt_entropy = 
  
  mouse1 = fread(paste0(WORKING_DIR, "COVID_murine_NetMHCpan_unfilt.xls"))
mouse1 = mouse1[which(mouse1$NB>0),]
mouse2 = fread(paste0(WORKING_DIR, "COVID_murine_NetMHCIIpan_unfilt.xls"))
mouse2 = mouse2[which(mouse2$NB>0),]
Paper_bc = fread(paste0(WORKING_DIR, "Supplemental/B_cell_epitopes - Bc_epitopes.txt"))
Paper_bc=Paper_bc[which(Paper_bc$Source == "PEPperCHIP"),]


####Grab number of epitopes per each protein for all the difference datasets####
Proteins=sapply(strsplit(unique(I_filt$ID),"_"),"[",1)#[order(sapply(strsplit(unique(I_filt$ID),"_"),"[",1))]
Count = c()
CountII = c()
Count_co = c()
Count_m1 = c()
Count_m2 = c()
Count_m_covered = c()
pep_count = c()
bc_count = c()

for(P in Proteins){
  Count = c(Count, length(which(sapply(strsplit(I_filt$ID,"_"),'[',1) == P)))
  CountII= c(CountII, length(which(sapply(strsplit(II_filt$ID,"_"),'[',1) == P)))
  Count_co= c(Count_co, length(which(lo_entropy$gene == P)))
  Count_m1 = c(Count_m1, length(which(sapply(strsplit(mouse1$ID,"_"),'[',1) == P)))
  Count_m2 = c(Count_m2, length(which(sapply(strsplit(mouse2$ID,"_"),'[',1) == P)))
  Count_m_covered = c(Count_m_covered, length(which((lo_entropy$gene == P) & (lo_entropy$murine != "None"))))
  pep_count = c(pep_count, length(which(Paper_bc$Protein == P)))
  bc_count = c(bc_count, length(which(bc_sub$gene == P)))
  
}


Proteins = factor(Proteins, levels =Proteins)
#QUick fix for orfl vs orf1
Count_co[1] = length(which(lo_entropy$gene == "orf1ab"))
Count_m_covered[1] = length(which((lo_entropy$gene == 'orf1ab') & (lo_entropy$murine != "None")))
Count_m_sum = Count_m1+Count_m2


#####Generate table for ggplot2, by protein#####
I_summary=data.table(Protein = Proteins, Count = Count )
II_summary=data.table(Protein = Proteins, Count = CountII )
co_summary=data.table(Protein = Proteins, Count = Count_co )
mI_summary=data.table(Protein = Proteins, Count = Count_m1 )
mII_summary=data.table(Protein = Proteins, Count = Count_m2 )
m_summary=data.table(Protein = Proteins, Count = Count_m_sum )
m_covered_summary=data.table(Protein = Proteins, Count = Count_m_covered )
pep_summary=data.table(Protein = Proteins, Count = pep_count )
bc_summary=data.table(Protein = Proteins, Count = bc_count )


#####Generate table for ggplot2, by protein grouping#####
I_grouped = data.table(Protein = c("Internal", "Surface", "Nucleocapsid"), Count = c(sum(Count[1],Count[2],Count[4],Count[5],Count[9],Count[10]),
                                                                                     sum(Count[3], Count[6], Count[8]), Count[7]))
II_grouped = data.table(Protein = c("Internal", "Surface", "Nucleocapsid"), Count = c(sum(CountII[1],CountII[2],CountII[4],CountII[5],CountII[9],CountII[10]),
                                                                                      sum(CountII[3], CountII[6], CountII[8]), CountII[7]))
co_grouped = data.table(Protein = c("Internal", "Surface", "Nucleocapsid"), Count = c(sum(Count_co[1],Count_co[2],Count_co[4],Count_co[5],Count_co[9],Count_co[10]),
                                                                                      sum(Count_co[3], Count_co[6], Count_co[8]), Count_co[7]))
mI_grouped = data.table(Protein = c("Internal", "Surface", "Nucleocapsid"), Count = c(sum(Count_m1[1],Count_m1[2],Count_m1[4],Count_m1[5],Count_m1[9],Count_m1[10]),
                                                                                      sum(Count_m1[3], Count_m1[6], Count_m1[8]), Count_m1[7]))
mII_grouped = data.table(Protein = c("Internal", "Surface", "Nucleocapsid"), Count = c(sum(Count_m2[1],Count_m2[2],Count_m2[4],Count_m2[5],Count_m2[9],Count_m2[10]),
                                                                                       sum(Count_m2[3], Count_m2[6], Count_m2[8]), Count_m2[7]))
m_grouped = data.table(Protein = c("Internal", "Surface", "Nucleocapsid"), Count = c(sum(Count_m_sum[1],Count_m_sum[2],Count_m_sum[4],Count_m_sum[5],Count_m_sum[9],Count_m_sum[10]),
                                                                                     sum(Count_m_sum[3], Count_m_sum[6], Count_m_sum[8]), Count_m_sum[7]))
m_covered_grouped = data.table(Protein = c("Internal", "Surface", "Nucleocapsid"), Count = c(sum(Count_m_covered[1],Count_m_covered[2],Count_m_covered[4],Count_m_covered[5],Count_m_covered[9],Count_m_covered[10]),
                                                                                             sum(Count_m_covered[3], Count_m_covered[6], Count_m_covered[8]), Count_m_covered[7]))
pep_grouped = data.table(Protein = c("Internal", "Surface", "Nucleocapsid"), Count = c(sum(pep_count[1],pep_count[2],pep_count[4],pep_count[5],pep_count[9],pep_count[10]),
                                                                                       sum(pep_count[3], pep_count[6], pep_count[8]), pep_count[7]))
bc_grouped = data.table(Protein = c("Internal", "Surface", "Nucleocapsid"), Count = c(sum(bc_count[1],bc_count[2],bc_count[4],bc_count[5],bc_count[9],bc_count[10]),
                                                                                      sum(bc_count[3], bc_count[6], bc_count[8]), bc_count[7]))


####Plot by protein############
ggplot(data = I_summary, aes(x=1,y=log(Count+1),fill=Protein))+geom_bar(stat='identity')+ scale_fill_viridis_d(option = "plasma")+
  theme_classic()+theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),legend.position="none",
                        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),plot.background=element_blank()) 

ggplot(data = II_summary, aes(x=1,y=log(Count+1),fill=Protein))+geom_bar(stat='identity')+ theme_classic()+ scale_fill_viridis_d(option = "plasma")+
  theme_classic()+theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),legend.position="none",
                        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),plot.background=element_blank())

ggplot(data = co_summary, aes(x=1,y=log(Count+1),fill=Protein))+geom_bar(stat='identity')+ theme_classic()+ scale_fill_viridis_d(option = "plasma")+
  theme_classic()+theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),legend.position="none",
                        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),plot.background=element_blank())#+coord_flip()

ggplot(data = mI_summary, aes(x=1,y=log(Count+1),fill=Protein))+geom_bar(stat='identity')+ theme_classic()+ scale_fill_viridis_d(option = "plasma")+
  theme_classic()+theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),legend.position="none",
                        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),plot.background=element_blank())#+coord_flip()
ggplot(data = mII_summary, aes(x=1,y=log(Count+1),fill=Protein))+geom_bar(stat='identity')+ theme_classic()+ scale_fill_viridis_d(option = "plasma")+
  theme_classic()+theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),legend.position="none",
                        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),plot.background=element_blank())#+coord_flip()
ggplot(data = m_summary, aes(x=1,y=log(Count+1),fill=Protein))+geom_bar(stat='identity')+ theme_classic()+ scale_fill_viridis_d(option = "plasma")+
  theme_classic()+theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),legend.position="none",
                        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),plot.background=element_blank())#+coord_flip()

ggplot(data = pep_summary, aes(x=1,y=log(Count+1),fill=Protein))+geom_bar(stat='identity')+ theme_classic()+ scale_fill_viridis_d(option = "plasma")+
  theme_classic()+theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),legend.position="none",
                        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),plot.background=element_blank())#+coord_flip()
ggplot(data = bc_summary, aes(x=1,y=log(Count+1),fill=Protein))+geom_bar(stat='identity')+ theme_classic()+ scale_fill_viridis_d(option = "plasma")+
  theme_classic()+theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),legend.position="none",
                        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),plot.background=element_blank())#+coord_flip()

ggplot(data = m_covered_summary, aes(x=1,y=log(Count+1),fill=Protein))+geom_bar(stat='identity')+ theme_classic()+ scale_fill_viridis_d(option = "plasma")+
  theme_classic()+theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),legend.position="none",
                        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),plot.background=element_blank())#+coord_flip()


#############Plot by group#######
ggplot(data = I_grouped, aes(x=1,y=log(Count+1),fill=Protein))+geom_bar(stat='identity')+ scale_fill_viridis_d(option = "plasma")+
  theme_classic()+theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),legend.position="none",
                        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),plot.background=element_blank()) 

ggplot(data = II_grouped, aes(x=1,y=log(Count+1),fill=Protein))+geom_bar(stat='identity')+ theme_classic()+ scale_fill_viridis_d(option = "plasma")+
  theme_classic()+theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),legend.position="none",
                        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),plot.background=element_blank())

ggplot(data = co_grouped, aes(x=1,y=log(Count+1),fill=Protein))+geom_bar(stat='identity')+ theme_classic()+ scale_fill_viridis_d(option = "plasma")+
  theme_classic()+theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),legend.position="none",
                        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),plot.background=element_blank())#+coord_flip()

ggplot(data = mI_grouped, aes(x=1,y=log(Count+1),fill=Protein))+geom_bar(stat='identity')+ theme_classic()+ scale_fill_viridis_d(option = "plasma")+
  theme_classic()+theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),legend.position="none",
                        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),plot.background=element_blank())#+coord_flip()
ggplot(data = mII_grouped, aes(x=1,y=log(Count+1),fill=Protein))+geom_bar(stat='identity')+ theme_classic()+ scale_fill_viridis_d(option = "plasma")+
  theme_classic()+theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),legend.position="none",
                        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),plot.background=element_blank())#+coord_flip()
ggplot(data = m_grouped, aes(x=1,y=log(Count+1),fill=Protein))+geom_bar(stat='identity')+ theme_classic()+ scale_fill_viridis_d(option = "plasma")+
  theme_classic()+theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),legend.position="none",
                        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),plot.background=element_blank())#+coord_flip()

ggplot(data = pep_grouped, aes(x=1,y=log(Count+1),fill=Protein))+geom_bar(stat='identity')+ theme_classic()+ scale_fill_viridis_d(option = "plasma")+
  theme_classic()+theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),legend.position="none",
                        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),plot.background=element_blank())#+coord_flip()
ggplot(data = bc_grouped, aes(x=1,y=log(Count+1),fill=Protein))+geom_bar(stat='identity')+ theme_classic()+ scale_fill_viridis_d(option = "plasma")+
  theme_classic()+theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),legend.position="none",
                        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),plot.background=element_blank())#+coord_flip()

ggplot(data = m_covered_grouped, aes(x=1,y=log(Count+1),fill=Protein))+geom_bar(stat='identity')+ theme_classic()+ scale_fill_viridis_d(option = "plasma")+
  theme_classic()+theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),legend.position="none",
                        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                        panel.grid.minor=element_blank(),plot.background=element_blank())#+coord_flip()






###############################################################
########PEPperPRINT heatmap Figure 1A###########################
library(gplots)

###Read in and format raw PEPperPRINT fluorescence matrix -- run once for each IgG/IgA
#PPP = fread(paste0(WORKING_DIR, "PEPperPRINT_IgA.txt"))
PPP = fread(paste0(WORKING_DIR, "PEPperPRINT_IgG.txt"))

#Formating
rownames(PPP) = PPP$Sample_ID
PPP=t(PPP)
colnames(PPP) = PPP[1,]
PPP=PPP[-1,]
PPP = as.matrix(PPP)
name_peps = rownames(PPP)
PPP=apply(PPP,2,as.numeric)
rownames(PPP) = name_peps

sample_colors = viridis(2)

###Generate colors by Naive vs Infected groups###
color_list = c()
for(z in 1:ncol(PPP)){
  if(substr(colnames(PPP)[z],1,1) == "N"){
    color_list=c(color_list, sample_colors[1])
  }else{
    color_list=c(color_list,sample_colors[2])
  }
}

#####Generate colors by protein############
Pep_key = fread(paste0(WORKING_DIR, "PEPperPRINT_peptide_key.txt"))
pep_col = viridis(3, begin = 0, end = .75, option = "plasma")
pep_col_list = c()
for(u in 1:nrow(PPP)){
  if(Pep_key$Protein[which(Pep_key$Peptide == rownames(PPP)[u])] == "Surface Glycoprotein"){
    pep_col_list = c(pep_col_list, pep_col[1])
  }else if(Pep_key$Protein[which(Pep_key$Peptide == rownames(PPP)[u])] == "Membrane Glycoprotein"){
    pep_col_list = c(pep_col_list, pep_col[2])
  }else if(Pep_key$Protein[which(Pep_key$Peptide == rownames(PPP)[u])] == "Nucleocapsid Phosphoprotein"){
    pep_col_list = c(pep_col_list, pep_col[3])
  }else{
    print(u)
  }
}

###Add start coordinates to peptide sequences
seq =  fread(paste0(WORKING_DIR, "AA_sequence.txt"), header = F)
S= seq$V1[4]
M = seq$V1[10]
N=seq$V1[18]


for(p in 1:nrow(PPP)){
  if(gregexpr(pattern = rownames(PPP)[p],S)[[1]][1] != -1){
    rownames(PPP)[p] = paste0(gregexpr(pattern = rownames(PPP)[p],S)[[1]][1], "-", rownames(PPP)[p])
  }else if(gregexpr(pattern = rownames(PPP)[p],M)[[1]][1] != -1){
    rownames(PPP)[p] = paste0(gregexpr(pattern = rownames(PPP)[p],M)[[1]][1], "-", rownames(PPP)[p])
  }else if(gregexpr(pattern = rownames(PPP)[p],N)[[1]][1] != -1){
    rownames(PPP)[p] = paste0(gregexpr(pattern = rownames(PPP)[p],N)[[1]][1], "-", rownames(PPP)[p])
  }    
}



######Plotting Figure 1A, Figure 1B#################
cm =  viridis(100)# red-white-blue colormap
mmx = 70000
colbr <- c(seq(0,500, len=length(cm)+1)) 

par(cex.main = 2)
heatmap.2(PPP, Rowv = F, Colv = F, trace = "none", scale = "none", col = cm, breaks = colbr, labCol = "",key.title = "", density.info = "none",
          margins = c(1, 15), srtCol = 45, cexRow = 1.3, dendrogram = "none", main = "SARS-COV-2 PEPperCHIP IgG", 
          ColSideColors = color_list,
          RowSideColors = pep_col_list,
          hclustfun = function(x) hclust(x,method = 'complete'),
          distfun = function(x) dist(x,method = 'euclidean',
                                     key=TRUE , key.xlab="Fluorescence")
)

legend("topright",      
       legend = c("Naive", "Infected"),
       col = c(sample_colors[1], sample_colors[2]), 
       lty= 1,             
       lwd = 10,           
       cex=1.2,bty="n"
)

legend("top",      
       legend = c("Surface", "Membrane", "Nucleocapsid"),
       col = pep_col, 
       lty= 1,             
       lwd = 10,           
       cex=1.2,
       bty="n",
       yjust=0
)


###################Lollipop plot of PEPperPRINT Figure 1C##########
##################################################################


####Adding in Bc data from literature###
#Paper_bc = fread(paste0(WORKING_DIR, "Paper_Bc_epitopes_complete.txt"))
Paper_bc = fread(paste0(WORKING_DIR, "Supplemental/B_cell_epitopes - Bc_epitopes.txt"))

###Remove these two artificial peptides (probably not necessary, these won't map anyway)###
Paper_bc = Paper_bc[-which(Paper_bc$Epitope %in% c('SGSGMADSNGTITVE',
                                                   'SGMADSNGTITVEEL') ),]



#####Match peptides with SARS-CoV-2 reference, same as before####
seq =  fread(paste0(WORKING_DIR, "AA_sequence.txt"), header = F)
S= seq$V1[4]
M = seq$V1[10]
N=seq$V1[18]

Paper_bc$Start = NA
Paper_bc$End = NA
for(p in 1:nrow(Paper_bc)){
  if(is.na(Paper_bc$Protein[p])){
    if (gregexpr(pattern = Paper_bc$Epitope[p],S)[[1]][1] != -1){
      Paper_bc$Protein[p] = "S"
    }else if (gregexpr(pattern = Paper_bc$Epitope[p],M)[[1]][1] != -1){
      Paper_bc$Protein[p] = "M"
    }else if (gregexpr(pattern = Paper_bc$Epitope[p],N)[[1]][1] != -1){
      Paper_bc$Protein[p] = "N"
    }
  }
  if(!is.na(Paper_bc$Protein[p])){
    if(Paper_bc$Protein[p] == "S"){
      if (gregexpr(pattern = Paper_bc$Epitope[p],S)[[1]][1] != -1){
        Paper_bc$Start[p] = 7096+gregexpr(pattern = Paper_bc$Epitope[p],S)[[1]][1]
        Paper_bc$End[p] = Paper_bc$Start[p] + nchar(Paper_bc$Epitope[p]) -1
      }
    }else if(Paper_bc$Protein[p] == "M"){
      if (gregexpr(pattern = Paper_bc$Epitope[p],M)[[1]][1] != -1){
        Paper_bc$Start[p] = 8719+gregexpr(pattern = Paper_bc$Epitope[p],M)[[1]][1]
        Paper_bc$End[p] = Paper_bc$Start[p] + nchar(Paper_bc$Epitope[p]) -1
      }
      
    }else if(Paper_bc$Protein[p] == "N"){
      if (gregexpr(pattern = Paper_bc$Epitope[p],N)[[1]][1] != -1){
        Paper_bc$Start[p] = 9244+gregexpr(pattern = Paper_bc$Epitope[p],N)[[1]][1]
        Paper_bc$End[p] = Paper_bc$Start[p] + nchar(Paper_bc$Epitope[p]) -1
      }
    }
  }  
}
Paper_bc$Protein = factor(Paper_bc$Protein, levels = c("S","M", "N"))


####Read in raw fluorescence data
res_IgA = fread(paste0(WORKING_DIR, "PEPperPRINT_IgA_DESeq2.txt"))
rownames(res_IgA) = res_IgA$V1
res_IgG = fread(paste0(WORKING_DIR, "PEPperPRINT_IgG_DESeq2.txt"))
rownames(res_IgG) = res_IgG$V1

Paper_bc$IgG_fluor_FC = NA
Paper_bc$IgA_fluor_FC = NA

#Pull data from above DESeq2 run
for (w in 1:nrow(Paper_bc)){
  if(Paper_bc$Epitope[w] %in% rownames(res_IgG)){
    Paper_bc$IgG_fluor_FC[w] = res_IgG$log2FoldChange[which(rownames(res_IgG) == Paper_bc$Epitope[w])] 
    Paper_bc$IgA_fluor_FC[w] = res_IgA$log2FoldChange[which(rownames(res_IgA) == Paper_bc$Epitope[w])] 
  }
}

Paper_bc = Paper_bc[which(!is.na(Paper_bc$Start)),]


############Plotting Figure 1C################

ggplot(data=Paper_bc) + 
  geom_point(aes(x=Start, y=IgG_fluor_FC, color=Protein, size=IgA_fluor_FC),alpha=1, shape = 1, stroke = 2) + 
  geom_segment(aes(x=Start, xend=Start, y=0, yend=IgG_fluor_FC), alpha=0.25, size=1.1) + 
  theme(text=element_text(face="bold",size=20,colour="black")) +
  geom_label_repel(aes(label=ifelse(IgG_fluor_FC>5,as.character(Epitope),''), x=Start, y=IgG_fluor_FC),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  #scale_color_viridis_d(begin = 0, end = .8)+
  scale_color_manual(values = pep_col)+
  scale_y_continuous(limits = c(-12,35))+
  #scale_color_viridis_c() + 
  # theme_light() +
  labs(y="IgG fluorcesence log(fold change)", x="Position across SARS-CoV-2 proteome (AA)", title = "SARS-CoV-2 PEPperCHIP epitopes") +
  scale_size(range = c(-3,30), breaks = c(0,10,20,30), name = "IgA fluorescence\nlog(fold change)")+
  
  geom_rect(aes(xmin=7096, xmax=8368, ymin=(-8), ymax=(-12)), fill="white", color="black", size=0.1) +
  geom_text(aes(x=7096+(8368-7096)/2, y=(-10), label="S", angle=0), size=5) + 
  
  geom_rect(aes(xmin=8369, xmax=8644, ymin=(-4), ymax=(-8)), fill="white", color="black", size=0.1) +
  geom_text(aes(x=8368+(8644-8368)/2, y=(-6), label="ORF3a", angle=0), size=5) +
  
  geom_rect(aes(xmin=8645, xmax=8718, ymin=(-8), ymax=(-12)), fill="white", color="black", size=0.1) +
  geom_text(aes(x=8645+(8718-8645)/2, y=(-10)), label="E", angle=0, size=5) + 
  
  geom_rect(aes(xmin=8719, xmax=8940, ymin=(-4), ymax=(-8)), fill="white", color="black", size=0.1) +
  geom_text(aes(x=8719+(8940-8719)/2, y=(-6)), label="M", angle=0, size=5) +
  
  geom_rect(aes(xmin=8941, xmax=9001, ymin=(-8), ymax=(-12)), fill="white", color="black", size=0.1) +
  geom_text(aes(x=8941+(9001-8941)/2, y=(-10)), label="ORF6", angle=0, size=2.4) +
  
  geom_rect(aes(xmin=9002, xmax=9122, ymin=(-4), ymax=(-8)), fill="white", color="black", size=0.1) +
  geom_text(aes(x=9002+(9122-9002)/2, y=(-6)), label="ORF7a", angle=0, size=4) +
  
  geom_rect(aes(xmin=9123, xmax=9243, ymin=(-8), ymax=(-12)), fill="white", color="black", size=0.1) +
  geom_text(aes(x=9123+(9243-9123)/2, y=(-10)), label="ORF8", angle=0, size=4.8) +
  
  geom_rect(aes(xmin=9244, xmax=9662, ymin=(-4), ymax=(-8)), fill="white", color="black", size=0.1) +
  geom_text(aes(x=9244+(9662-9244)/2, y=(-6)), label="N", angle=0, size=5) 



####################################################################
###Figure 5A top, plotting all Bc eps###############################

###This pulls from above code, so run through code for Figure 1C to this point first###

###Combine PEPperCHIP sets for now
Paper_bc$Source[which(Paper_bc$Source == "PEPperCHIP_etc")] = "PEPperCHIP"
validated = c("PEPperCHIP", "Wang et al. 2020", "Poh et al. 2020") #These are sets that were derived with array studies

Paper_bc$Source = factor(Paper_bc$Source, levels = 
                           c(validated,unique(Paper_bc$Source)[which(!(unique(Paper_bc$Source) %in% validated))]))

Paper_bc=Paper_bc[-which(duplicated(Paper_bc)),] #Remove duplicated values


ggplot(data=Paper_bc) + 
  geom_jitter(aes(x=Start, y=Source, color=Source, fill = Source), shape = 20, stroke=2, size=10, alpha=.5, width = 0, height = .1)+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  scale_color_viridis_d(option = "plasma", end = .75)+
  geom_hline(yintercept = .5, size=3, color="white")+
  geom_hline(yintercept = 1.5, size=3, color="white")+
  geom_hline(yintercept = 3.5, size=3, color="white")+
  geom_vline(xintercept = ifelse(duplicated(Paper_bc$Epitope),Paper_bc$Start,0), color = "red") +
  
  geom_label_repel(aes(label=ifelse(!isUnique(Epitope),as.character(Epitope),''), x=Start, y=Source),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50', size=2.7)+
  
  labs(y='',x="Position across SARS-CoV-2 proteome (AA)", title = "Curated B cell epitopes") +
  
  geom_rect(aes(xmin=7096, xmax=8368, ymin=0, ymax=(-1)), fill="white", color="black", size=0.1) +
  geom_text(aes(x=7096+(8368-7096)/2, y=(-.7), label="S", angle=0), size=5) + 
  
  geom_text(aes(x=7096+100, y=(-.25), label="S100", angle=90), size=2.5) + 
  geom_text(aes(x=7096+400, y=(-.25), label="S400", angle=90), size=2.5) + 
  geom_text(aes(x=7096+700, y=(-.25), label="S700", angle=90), size=2.5) + 
  geom_text(aes(x=7096+1000, y=(-.3), label="S1000", angle=90), size=2.5) + 
  
  
  geom_rect(aes(xmin=8719, xmax=8940, ymin=(0), ymax=(-1)), fill="white", color="black", size=0.1) +
  geom_text(aes(x=8719+(8940-8719)/2, y=(-.7)), label="M", angle=0, size=5) +
  geom_text(aes(x=8719+50, y=(-.25), label="M50", angle=90), size=2.5) + 
  geom_text(aes(x=8719+150, y=(-.25), label="M150", angle=90), size=2.5) + 
  
  geom_rect(aes(xmin=9244, xmax=9662, ymin=(0), ymax=(-1)), fill="white", color="black", size=0.1) +
  geom_text(aes(x=9244+(9662-9244)/2, y=(-.7)), label="N", angle=0, size=5)+
  geom_text(aes(x=9244+100, y=(-.25), label="N100", angle=90), size=2.5) + 
  geom_text(aes(x=9244+200, y=(-.25), label="N200", angle=90), size=2.5) + 
  geom_text(aes(x=9244+300, y=(-.25), label="N300", angle=90), size=2.5) + 
  geom_text(aes(x=9244+400, y=(-.25), label="N400", angle=90), size=2.5) + 
  scale_y_discrete(position = "right")+
  scale_x_continuous(limits = c(7000,9662))


##########Filters (represented by Figure 5A midde)#################

#By neutralization
Paper_bc = Paper_bc[which((Paper_bc$`In vitro function` =='Neutralizing') |
                            (Paper_bc$Notes == "Array") |
                            (Paper_bc$Source == "PEPperCHIP")),]

#By RBD and FP
in_FP = which((Paper_bc$Start >= (7096+788)) & (Paper_bc$End <= (7096+806+100)))
in_RBD = which((Paper_bc$Start >= (7096+329)) & (Paper_bc$End <= (7096+521+100)))
Paper_bc = Paper_bc[c(in_FP, in_RBD),]

#By glycosylation
glyc = 7096+c(17,61,74,122,149,165,234,282,331,343, 603, 616,657,709,717,801,1074,1098,1134,1158,1173,1194)

for(n in glyc){
  hits = which((Paper_bc$Start <= n) & (Paper_bc$End >= n))
  if(length(hits)>0){
    Paper_bc = Paper_bc[-hits,  ]
  } 
}


#By polymorphisms
poly=c(341,342,354,364,367,408,435,436,483)
poly=poly+7096

for(n in poly){
  hits = which((Paper_bc$Start <= n) & (Paper_bc$End >= n))
  if(length(hits)>0){
    Paper_bc = Paper_bc[-hits,  ]
  } 
}


#####################################################
##Figure 5A bottom -- Post-filter####################


#####Look for overlap w/Tc co-epitopes#####

joint_run_filt_entropy = fread(paste0(WORKING_DIR, "COVID_combined_entropy_rank_filtered.csv"))
joint_run_filt_entropy=joint_run_filt_entropy[which(joint_run_filt_entropy$lo_entropy==1),]

library(GenomicRanges)

g1 <-  GRanges(seqnames="COVID",
               IRanges(start=joint_run_filt_entropy$c2_start,
                       end=(joint_run_filt_entropy$c2_start + nchar(joint_run_filt_entropy$c2_peptide) -1)),
               Pep_Tc_1 = joint_run_filt_entropy$c1_peptide,
               Pep_Tc_2 = joint_run_filt_entropy$c2_peptide)


g2 <-  GRanges(seqnames="COVID",
               IRanges(start=Paper_bc$Start,
                       end=Paper_bc$End),
               Pep_Bc = Paper_bc$Epitope,
               Protein = Paper_bc$Protein,
               Lit_source = Paper_bc$Source)

merge=as.data.table(mergeByOverlaps(g1,g2))

<<<<<<< HEAD
Paper_bc$Coep_coverage = Paper_bc$Epitope %in% merge$Pep_Bc
=======

>>>>>>> 8641edb347a0984de8f29feccce730866780fec5
####Look for overlap with HLA I

HLAI_filt=fread(paste0(WORKING_DIR, "COVID_human_netMHCpan_rank_filtered.txt"))
HLAI_filt=HLAI_filt[which(HLAI_filt$lo_entropy==1),] #Filter by entropy
HLAI_filt=HLAI_filt[which(HLAI_filt$ID=="S_surface_glyco"),]  #Keep only S-derived Tc epitopes (only protein in these data)
HLAI_filt$Pos=HLAI_filt$Pos+7097  #Transform to proteomic space
HLAI_filt$End = HLAI_filt$Pos+nchar(HLAI_filt$Peptide)-1

g1 <-  GRanges(seqnames="COVID",
               IRanges(start=HLAI_filt$Pos,
                       end=(HLAI_filt$End)),
               Pep_Tc = HLAI_filt$Peptide)


g2 <-  GRanges(seqnames="COVID",
               IRanges(start=Paper_bc$Start,
                       end=Paper_bc$End),
               Pep_Bc = Paper_bc$Epitope,
               Protein = Paper_bc$Protein,
               Lit_source = Paper_bc$Source)

merge1=as.data.table(mergeByOverlaps(g1,g2))

Paper_bc$HLA_I_coverage=Paper_bc$Epitope %in% merge1$Pep_Bc 

####Look for overlap with HLA II

HLAII_filt=fread(paste0(WORKING_DIR, "COVID_human_netMHCIIpan_rank_filtered.txt"))
HLAII_filt=HLAII_filt[which(HLAI_filt$lo_entropy==1),]

HLAII_filt=HLAII_filt[which(HLAII_filt$ID=="S_surface_glycoprotein_QIB84673.1"),]
HLAII_filt$Pos=HLAII_filt$Pos+7097
HLAII_filt$End = HLAII_filt$Pos+nchar(HLAII_filt$Peptide)-1

g1 <-  GRanges(seqnames="COVID",
               IRanges(start=HLAII_filt$Pos,
                       end=(HLAII_filt$End)),
               Pep_Tc = HLAII_filt$Peptide)


g2 <-  GRanges(seqnames="COVID",
               IRanges(start=Paper_bc$Start,
                       end=Paper_bc$End),
               Pep_Bc = Paper_bc$Epitope,
               Protein = Paper_bc$Protein,
               Lit_source = Paper_bc$Source)

merge2=as.data.table(mergeByOverlaps(g1,g2))

Paper_bc$HLA_II_coverage=Paper_bc$Epitope %in% merge2$Pep_Bc


######Transform to S-space#####
Paper_bc$Start = Paper_bc$Start-7096
Paper_bc$End = Paper_bc$End-7096


##########Plotting 5A bottom##############

Paper_bc = Paper_bc[order(Paper_bc$Start),]
Paper_bc$Segment = seq(.25,(nrow(Paper_bc)/4),.25)+.75  #Sets the y-axis scale for plotting

#These are the peptides in the final vaccine set
Peptide= c('YRLFRKSNLKPFERD',
           'SNLKPFERDISTEIY',
           'YQPYRVVVLSFELLHAPA',
           'IADTTDAVRDPQTLEILDI',
           'TESNKKFLPFQQFGRDIA',
           'LPDPSKPSKRSFIEDLLFNKV')
Start=c(453,459,505,569,553,806) #Their respective start coordinates
End=c(467,473,522,587,570,826) #Respective end coordinates
Color=c("red", "green", "blue", "magenta", "yellow", "orange") #Defining colors for each
Segment=c(.25,.5,.25,.25,.5,.25) #Sets the y-axis scale for plotting
final_set = data.table(Peptide, Start, End, Color, Segment) #Merging into a dt
final_set$Peptide = factor(final_set$Peptide, levels = Peptide)


ggplot(data=Paper_bc) + 
  geom_segment(aes(x=Start, xend=End, y=Segment, yend=Segment, color=Source),size=3, alpha=1)+
  theme(text=element_text(face="bold",size=20,colour="black")) +
  scale_color_viridis_d(option = "plasma", end = .75)+
  geom_hline(yintercept = .75, size=3, color="white")+
  ggnewscale::new_scale("color")+
  
  geom_segment(data=final_set,aes(x=Start, xend=End, y=Segment, yend=Segment, color=Peptide),size=3, alpha=1)+
  scale_color_manual(values = Color)+
  
  labs(y='',x="Position across SARS-CoV-2 spike protein", title = "Filtered B cell epitopes") +
  
  geom_rect(aes(xmin=(400), xmax=(521), ymin=0, ymax=(-1)), fill="white", color="black", size=0.1) +
  geom_text(aes(x=(400+521)/2, y=(-.5), label="RBD", angle=0), size=5) +
  
  geom_rect(aes(xmin=788, xmax=806, ymin=0, ymax=(-1)), fill="white", color="black", size=0.1) +
  geom_text(aes(x=(788+806)/2, y=(-.5), label="FP", angle=0), size=4.5) +
  scale_x_continuous(limits = c(400,910))+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())



############################################################
###########Figure 2B boxplots for entropy by protein########
############################################################
library(ggbeeswarm)
library(ggallin)


#####Plotting 2B##############
#ent = fread(paste0(WORKING_DIR, "entropy.MT072688.1.corrected.csv"))
#ent = fread(paste0(WORKING_DIR, "AA_mafft_df_entropy.txt"))
ent = fread(paste0(WORKING_DIR, "entropy_7882.txt"))

colnames(ent) = c("position", "entropy")

#Defining coordinates
ent$Protein = NA
ent$Protein[1:7096] = "orf1ab"
ent$Protein[7097:8369] = "S"
ent$Protein[8370:8644] = "ORF3a"
ent$Protein[8645:8719] = "E"
ent$Protein[8720:8941] = "M"
ent$Protein[8942:9002] = "ORF6"
ent$Protein[9003:9123] = "ORF7a"
ent$Protein[9124:9244] = "ORF8"
ent$Protein[9245:9663] = "N"
ent$Protein[9664:9701] = "ORF10"
ent$Protein=factor(ent$Protein, levels = unique(ent$Protein))

#ent=ent[which(ent$entropy>0),]
ent=ent[which(!is.na(ent$Protein)),]  #Filter if not found -- probably not necessary

#Plotting -- okay if some errors for missing values -- I cut some off using y-scale
ggplot(data = ent, aes(x = Protein, y = entropy, fill = Protein))+
  geom_violin()+
  geom_quasirandom()+
  scale_y_continuous(trans=ssqrt_trans, limits = c(0,.1))+#, labels = scales::number_format(accuracy = 0.01))+
  scale_fill_viridis_d(option = "plasma")+
  theme(text=element_text(face="bold",size=20,colour="black"))+
  labs(y="Entropy", x= "SARS-CoV-2 protein")

#####See % > 0.01 by protein#####
for(n in unique(ent$Protein)){
  print(paste0(n,": ",length(which(ent$entropy[which(ent$Protein==n)]>=0.01))/length(which(ent$Protein==n)) ))
  
}
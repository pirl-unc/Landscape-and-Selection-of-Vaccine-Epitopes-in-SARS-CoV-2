require(msa)
require(seqinr)
require(data.table)
require(housekeeping)

# BiocManager::install(c("msa", "seqinr"))

# Code to run on the command line
# args = commandArgs(trailingOnly=TRUE)
# msa = fread(args[1])
# out_file = args[2]

this_script_path = housekeeping::get_script_dir_path(include_file_name = T)

# msa = full_msa[1:5, ]
rownames(msa) = msa[,c(V1)]
msa$`V1` = NULL
n_count = c()
for(q in 1:nrow(msa)){
  n_count = c(n_count, length(which(msa[q,] == "N")))
}
trans_tab = matrix(nrow=nrow(msa), ncol=9701)
rownames(trans_tab) = rownames(msa)
trans_coord = fread("/Users/sarahentwistle/Documents/Nextstrain/Translation_coords.txt")
temp_fa = matrix(nrow = 20,ncol=1)
#Remove stop codons
remove= c( 7097, 8371, 8647, 8723, 8946, 9008, 9130, 9252, 9672, 9711)
for(m in 1:nrow(msa)){
  print(m)
  AA_temp=c()
  for(n in 1:nrow(trans_coord)){
    seq_temp = s2c(tolower(paste0(lapply(msa[m,trans_coord$Start[n]:trans_coord$End[n]],as.character),collapse="")))
    #print(length(seq_temp))
    #print(length(seqinr::translate(seq_temp)))
    # print(AA_temp)
    AA_temp=c(AA_temp,seqinr::translate(seq_temp))
  }
  #print(length(AA_temp))
  AA_temp = AA_temp[-remove]  
  # print(length(AA_temp))
  trans_tab[m,] = AA_temp
}
trans_tab_fin = cbind(rownames(trans_tab), trans_tab)
colnames(trans_tab_fin) = c("SeqID", seq(1,9701,1))

fwrite(trans_tab_fin, args[2], sep = "\t")

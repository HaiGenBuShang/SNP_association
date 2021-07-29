#all SNPs in LD
library(data.table)
library(R.utils)

setwd("d:/Wang_lab/Project/GWAS/GWAS_summary_statistics_data/")
LD_SNPs <- c("rs197412")

summary_data_files <- list.files(path = ".",pattern = "h.tsv.gz",recursive = TRUE,full.names = TRUE)
disease_types <- sapply(strsplit(summary_data_files,split = "\\/"),function(x)x[2])
arranged_files <- split(summary_data_files,f = disease_types)

res <- lapply(arranged_files,function(x){
  sapply(x,function(y){
    tmp_file <- fread(file = y,sep = "\t",header = TRUE,stringsAsFactors = FALSE,
                      data.table = FALSE)
    tmp_file <- tmp_file[tmp_file$hm_rsid%in%LD_SNPs,]
    rownames(tmp_file) <- tmp_file$hm_rsid
    res_p_vale <- tmp_file[LD_SNPs,"p_value"]
    names(res_p_vale) <- LD_SNPs
    res_p_vale
  })
})
gc()

res_matrix <- do.call(cbind,args =  res)
res_matrix <- -log10(res_matrix)

colname_first_part <- sapply(strsplit(colnames(res_matrix),split = "\\/"),function(x)x[2])
colname_first_part[1:2] <- "Breast_Cancer"
colname_2nd_part <- sapply(strsplit(gsub(".*\\/","",colnames(res_matrix)),split = "-"),function(x)x[2])

colnames(res_matrix) <- paste0(colname_first_part,"_",colname_2nd_part)

pdf("LD_SNP_p_value.pdf",width = 10.518,height = 6.5)
opar <- par(no.readonly = TRUE)
par(mai=c(3.82,0.82,0.82,0.42))
barplot(res_matrix,beside = TRUE,las=2,col = "grey",main = "FGL1 LD SNPs",ylab = "P-value (-log10)")
abline(h=-log10(0.05),lwd=2,lty=2)
dev.off()

res_matrix <- do.call(cbind,args =  res)
colname_first_part <- sapply(strsplit(colnames(res_matrix),split = "\\/"),function(x)x[2])
colname_first_part[1:2] <- "Breast_Cancer"
colname_2nd_part <- sapply(strsplit(gsub(".*\\/","",colnames(res_matrix)),split = "-"),function(x)x[2])

colnames(res_matrix) <- paste0(colname_first_part,"_",colname_2nd_part)
write.table(res_matrix[,15:24],file ="FGL1_10_LD_SNPs_p_value.txt",col.names = TRUE,row.names = TRUE,sep = "\t",
            quote = FALSE)





res <- lapply(arranged_files,function(x){
  sapply(x,function(y){
    tmp_file <- fread(file = y,sep = "\t",header = TRUE,stringsAsFactors = FALSE,
                      data.table = FALSE)
    tmp_file <- tmp_file[tmp_file$hm_rsid%in%LD_SNPs,]
    rownames(tmp_file) <- tmp_file$hm_rsid
    res_p_vale <- tmp_file[LD_SNPs,c("hm_odds_ratio")]
    names(res_p_vale) <- LD_SNPs
    res_p_vale
  })
})
gc()

res_matrix <- do.call(cbind,args =  res)
res_matrix <- -log10(res_matrix)

colname_first_part <- sapply(strsplit(colnames(res_matrix),split = "\\/"),function(x)x[2])
colname_first_part[1:2] <- "Breast_Cancer"
colname_2nd_part <- sapply(strsplit(gsub(".*\\/","",colnames(res_matrix)),split = "-"),function(x)x[2])

colnames(res_matrix) <- paste0(colname_first_part,"_",colname_2nd_part)
write.table(res_matrix[,15:24],file ="FGL1_10_LD_SNPs_odds_ratio.txt",col.names = TRUE,row.names = TRUE,sep = "\t",
            quote = FALSE)















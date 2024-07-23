library(tidyverse)
library(plyr)
library(dplyr)
library(reshape2)
library(biomaRt)
library(mgsub)
library(data.table)
library(UpSetR)

setwd("/Users/Haley/code/RMAE/240508_COMBINE_NEW")
PID_genes <- as.character(read.table(paste0("files/PID_list_457_clean.csv"), sep = ",", header = FALSE)$V1)

## define filters to consider
fdr_thresh <- 0.05
lower_thresh <- c(0.3, 0.2, 0.1, 0.05)
COV <- 30

types <- paste0(c("both_SNP","missense_SNP","synonymous_SNP"),COV)

for(j in 1:length(types)){

      type_i <- types[j]

      IN_dir <- paste0("results/240520_cutoffs_summary/filtered_SNP/")
      OUT_dir <- paste0("results/240520_cutoffs_summary/",type_i,"/")
      system(paste0("mkdir -p ",OUT_dir,"lists/"))

      ## read in data
      MAE_results <- read.table(paste0(IN_dir,"results_summary_SNPmeanCOV",COV,"_",type_i,".txt"), header = TRUE, sep = ",")
      MAE_results <- subset(MAE_results, select = c(ID, gene, sumCOV, AI, fdr_BH, AI_30, AI_20, AI_10, AI_05))

      ## write all genes
      genes <- unique(MAE_results$gene)
      write.table(genes, paste0(OUT_dir,"lists/all_genes_meanCOV",COV,".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")

      PID_genes_in_dataset <- PID_genes[PID_genes %in% genes]
      write.table(PID_genes_in_dataset, paste0(OUT_dir,"lists/all_PID_genes_meanCOV","_",type_i,".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")

      for(i in 1:length(lower_thresh)){

        thresh_i <- lower_thresh[i]
        df <- subset(MAE_results, select = c(ID, gene, sumCOV, AI, fdr_BH))
        df$bin <- ifelse(df$fdr_BH < fdr_thresh & df$AI <= lower_thresh[i] | df$fdr_BH < fdr_thresh & df$AI >= 1 - lower_thresh[i], "MAE","BAE")
        df$bin <- factor(df$bin, levels = c("MAE","BAE"))

        ## for the different thresholds, how many MAE genes in at least one clone?
        df_sub <- subset(df, bin %in% "MAE")
        MAE_genes <- unique(df_sub$gene)
        ## intersect MAE and PID genes
        MAE_PID_genes <- intersect(PID_genes_in_dataset, MAE_genes)

        ## n_MAE = the number of MAE genes present in at least one clone across the dataset 
        ## n_total = the number of genes detected in at least one clone across the dataset
        ## n_PID_MAE = the number of intersect(PID, MAE) genes present in at least one clone across the dataset

        overall <- data.frame(cutoff_thresh = paste0(thresh_i*100,"%/",(1-thresh_i)*100,"%"), lower_thresh = thresh_i, n_MAE = length(MAE_genes), n_total = length(genes), n_PID_MAE = length(MAE_PID_genes), prop_MAE = length(MAE_genes)/length(genes), prop_MAE_PID_within_all_genes = length(MAE_PID_genes)/length(genes), prop_MAE_PID_within_PID_genes = length(MAE_PID_genes)/length(PID_genes_in_dataset))

        if(i == 1){
          overall_outs <- overall
        }else{
          overall_outs <- rbind(overall_outs, overall)
        }

        ## write MAE gene lists 
        write.table(MAE_genes, paste0(OUT_dir,"lists/MAE_genes_thresh",thresh_i,"_meanCOV",COV,"_",type_i,".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")
        write.table(MAE_PID_genes, paste0(OUT_dir,"lists/MAE_PID_genes_thresh",thresh_i,"_meanCOV",COV,"_",type_i,".txt"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ",")

        ## how many MAE genes per sample?
        per_samp <- as.data.frame(table(df$ID, df$bin))
        per_samp <- reshape(per_samp, idvar = "Var1", timevar = "Var2", direction = "wide")
        colnames(per_samp) <- c("sample","n_MAE","n_BAE")
        per_samp$n_total <- per_samp$n_MAE + per_samp$n_BAE
        per_samp$prop_MAE <- per_samp$n_MAE / per_samp$n_total
        per_samp$prop_BAE <- per_samp$n_BAE / per_samp$n_total
        per_samp$cutoff_thresh <- paste0(thresh_i*100,"%/",(1-thresh_i)*100,"%")
        per_samp$low_thresh <- thresh_i

        if(i == 1){
          per_samp_outs <- per_samp
        }else{
          per_samp_outs <- rbind(per_samp_outs, per_samp)
        }
        print(thresh_i)
      }

      write.table(overall_outs, paste0(OUT_dir,"summary_across-dataset_meanCOV",COV,"_",type_i,".txt"), quote = FALSE, row.names = FALSE, sep = ",")
      write.table(per_samp_outs, paste0(OUT_dir,"summary_per-sample_meanCOV",COV,"_",type_i,".txt"), quote = FALSE, row.names = FALSE, sep = ",")
}




library(tidyverse)
library(plyr)
library(dplyr)
library(reshape2)
library(mgsub)
library(data.table)
library(Qllelic)
library(gridExtra)
ul = function(x){x[1:6, 1:6]}

setwd("/Users/Haley/code/RMAE/240508_COMBINE_NEW")

IN_dir <- "results/240519_QC/unfiltered_SNP/"
RESULTS_dir <- "results/240520_AI_qllelic/unfiltered_SNP/"
COV <- 30

types <- c("both","missense","synonymous")

for(i in 1:length(types)){

      type_i <- types[i]
      IN_dir_i <- paste0(IN_dir,type_i,"/")
      RESULTS_dir_i <- paste0(RESULTS_dir,type_i,"/")
      system(paste0("mkdir -p ",RESULTS_dir_i))

      rawtables <- readRDS(paste0(IN_dir_i,"filtered_meanCOV",COV,"_raw_counts.rds"))

      ## write function to 'batch' Reduce the raw allele count tables and calculate allelic imbalance fractions
      ## AI is calculated via refCount/ (refCount + altCount)
      batchReduce <- function(counts){
        Reduce(
          function(x,y){
            merge(x, y, by = "gene")
            },
          lapply(1, function(i){
            CountsToAI(df = counts, reps = i, thr = 1)
          })
        )
      }

      rawtables_AI <- lapply(rawtables, batchReduce)

      aicovTable <- list()
      for(j in 1:length(rawtables)){

        rawtables_j <- rawtables[[j]]
        rawtables_AI_j <- rawtables_AI[[j]]
        merge_j <- cbind(rawtables_j, rawtables_AI_j)
        merge_j[,11] <- NULL

        aicovTable[[j]] <- merge_j
      }
    
      names(aicovTable) <- names(rawtables)

      ## histogram visualization of AI distrubution for an example donor
      ## highlight imprinted genes
      AI_df <- bind_rows(aicovTable, .id = "ID")
      clones <- c("D02c-Clone-E10a4","D02c-Clone-E4a4","D1F-Clone-CD4a","D1F-Clone-CD8a", "D3-Clone-E1", "D311-Clone-A11a8","D311-Clone-E4a4", "D311-Clone-B3a8", "D311-Clone-B3b8", "D311-Clone-F2a4", "D311-Clone-F2b4", "D311-Clone-F3a4", "D311-Clone-F3b4", "D311-Clone-F9a4","D4F-Clone-CD4a", "D5U-Clone-CD4a", "D5U-Clone-CD4b", "D5U-Clone-CD8a", "D5U-Clone-CD8b", "D5U-Clone-CD8c", "D68-Clone-CD8d", "D68-Clone-E3a4","D68-Clone-F12a4")
      new <- c("D108-1","D108-10", "D108-4", "D108-5", "D108-6", "D108-7", "D108-9", "D113-21", "D113-22", "D113-25", "D113-26", "D113-32", "D41-18", "D41-27", "D41-28", "D41-29", "D41-30", "D41-30-2", "D41-31", "D41-32", "D41-33", "D68-17", "D68-27", "D68-30", "D68-33", "D68-8-1", "D68-8-3", "D96F-16", "D96F-18", "D96F-19", "D96F-21", "D96F-3", "D96F-4", "D96F-8")
      AI_df$highlight <- ifelse(AI_df$ID %in% clones, "QC issue",ifelse(AI_df$ID %in% new, "new", "all good"))
      AI_df$highlight <- factor(AI_df$highlight, levels = c("QC issue","all good","new"))
      dd <- subset(AI_df, select = c(ID, highlight))
      dd <- unique(dd)
      AI_dist <- ggplot() +
                      geom_rect(data = dd, aes(fill = highlight), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.1) +
                      geom_histogram(data = AI_df, aes(x = AI, y = ..density..), color = "grey50", fill = "white") +
                      geom_density(data = AI_df, aes(x = AI), color = "#0F52BA", linewidth = 0.65) +
                      theme_minimal() +
                      theme(panel.grid.major = element_blank(), 
                            panel.grid.minor = element_blank(), 
                            legend.position = "none", 
                            legend.title = element_blank(), 
                            panel.border = element_rect(color = "black", fill = NA, size = 1),
                            strip.background = element_blank(),
                            axis.line=element_line()) +
                      facet_wrap(facets = "ID", scales = "free") + 
                      scale_fill_manual(values = c("red", "white", "blue")) +
                      geom_vline(xintercept = c(0.25, 0.75), linetype = "dashed", color = "grey70") + 
                      geom_vline(xintercept = c(0.5), linetype = "dashed", color = "gray30") 
      pdf(paste0(RESULTS_dir_i,"AI_dists_per_sample_COV",COV,".pdf"), height = 13, width = 15)
      print(AI_dist)
      dev.off()

      saveRDS(aicovTable, file = paste0(IN_dir_i,"aicovTables_meanCOV",COV,".rds"))
}



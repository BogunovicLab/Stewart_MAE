library(tidyverse)
library(plyr)
library(dplyr)
library(reshape2)
library(biomaRt)
library(mgsub)
library(data.table)
library(Qllelic)
library(UpSetR)
library(RColorBrewer)
library(gridExtra)
ul = function(x){x[1:6, 1:6]}

setwd("/Users/Haley/code/RMAE/240508_COMBINE_NEW")

IN_dir <- "results/240519_QC/"
RESULTS_dir <- "results/240520_cutoffs_summary/filtered_SNP/"
COV <- 30
types <- c("both_SNP30")

keep <- read.table(paste0("results/240519_QC/keep_list.txt"), sep = ",")
keep <- keep$V1

for(j in 1:length(types)){

      plots_list <- list()
      type_i <- types[j]
      OUT_dir_i <- paste0("results/240520_cutoffs_summary/",type_i,"/")

      ## read in MAE genes
      res <- read.table(paste0(RESULTS_dir,"results_summary_SNPmeanCOV",COV,"_",type_i,".txt"), sep = ",", header = TRUE)

      ## plot correlations of AI per clones within donor
      aicovTable <- readRDS(file = paste0(IN_dir,type_i,"/aicovTables_meanCOV",COV,".rds"))
      aicovTable_HC <- aicovTable[names(aicovTable) %in% keep]

      D3 <- c("D3-Clone-E7","D3-Clone-E4")
      D37 <- c("D37-Clone-B11","D37-Clone-D11","D37-Clone-G3","D37-Clone-G7","D37-Clone-H17", "D37-Clone-H8")
      D41 <- c("D41-Clone-CD8a","D41-Clone-CD8b","D41-Clone-E11","D41-Clone-H5","D41-Clone-H8")
      D4 <- c("D4F-Clone-CD8a","D4F-Clone-CD8b")
      D68 <- c("D68-Clone-A6","D68-Clone-B4","D68-Clone-CD8a","D68-Clone-CD8b","D68-Clone-CD8c","D68-Clone-D10","D68-Clone-F1","D68-Clone-G3")
      D114 <- c("D114-c108-10", "D114-c108-5", "D114-c108-6", "D114-c108-7", "D114-c108-9", "D114-c113-21", "D114-c41-18", "D114-c41-30", "D114-c41-30-2", "D114-c41-32", "D114-c96F-16", "D114-c96F-19", "D114-c96F-21", "D114-c96F-3")

      D3 <- apply(combn(D3,2),2,paste,collapse='_')
      D37 <- apply(combn(D37,2),2,paste,collapse='_')
      D41 <- apply(combn(D41,2),2,paste,collapse='_')
      D4 <- apply(combn(D4,2),2,paste,collapse='_')
      D68 <- apply(combn(D68,2),2,paste,collapse='_')
      D114 <- apply(combn(D114,2),2,paste,collapse='_')
      all <- c(D3, D37, D41, D4, D68, D114)

      pairs_1 <- sapply(strsplit(all, "_"), `[`, 1)
      pairs_2 <- sapply(strsplit(all, "_"), `[`, 2)

      for(i in 1:length(pairs_1)){

            pair_1 <- pairs_1[i]; pair_1_dot <- paste0(gsub("-", ".", pair_1),".")
            pair_2 <- pairs_2[i]; pair_2_dot <- paste0(gsub("-", ".", pair_2),".")
            donor <- str_extract(pair_1, "[^-]+")

            df_1 <- as.data.frame(aicovTable_HC[names(aicovTable_HC) %in% pair_1]); colnames(df_1)[10] <- "gene"
            df_2 <- as.data.frame(aicovTable_HC[names(aicovTable_HC) %in% pair_2]); colnames(df_2)[10] <- "gene"

            df_1 <- subset(df_1, select = c("gene", paste0(pair_1_dot,"AI"), paste0(pair_1_dot,"meanCOV")))
            df_2 <- subset(df_2, select = c("gene", paste0(pair_2_dot,"AI"), paste0(pair_2_dot,"meanCOV")))

            MAE_1 <- subset(res, ID %in% pair_1)
            MAE_1 <- unique(subset(MAE_1, AI_20 %in% "MAE")$gene)
            MAE_2 <- subset(res, ID %in% pair_2)
            MAE_2 <- unique(subset(MAE_2, AI_20 %in% "MAE")$gene)
            MAE_comb <- unique(c(MAE_1, MAE_2))

            df <- join(df_1, df_2, by = "gene", type = "inner")
            colnames(df) <- c("gene","clone1_AI","clone1_COV","clone2_AI","clone2_COV")
            df <- subset(df, gene %in% MAE_comb)
            df$meanCOV <- rowMeans(df[,c("clone1_COV", "clone2_COV")])

            fit <- lm(clone2_AI ~ clone1_AI, data = df)
            R <- cor(df$clone1_AI, df$clone2_AI, method = c("pearson"))
            plot <- ggplot(df, aes(x = clone1_AI, y = clone2_AI, col = meanCOV)) +
                          geom_vline(xintercept = 0.5, color = "grey50") +
                          geom_hline(yintercept = 0.5, color = "grey50") +
                          geom_smooth(method = "lm", color = "#4D004B", linetype = "dashed", size = 0.5) +
                          geom_point() +
                          xlab(paste0("allelic imbalance\nclone ",gsub("-Clone","",pair_1))) +
                          ylab(paste0("allelic imbalance\nclone ",gsub("-Clone","",pair_2))) +
                          #geom_abline(intercept = 0, slope = 1, color = "blue", alpha = 0.3) +
                          scale_color_gradient("mean read\ndepth", low="#BFD3E6", high="#4D004B") + 
                          theme_bw() +
                          theme(panel.border = element_rect(colour = "black", fill = NA, size = .75), 
                              plot.title = element_text(size = 10), 
                              #panel.grid.major = element_blank(), 
                              #panel.grid.minor = element_blank(),
                              legend.position = "right",
                              axis.title.x = element_text(size = 10),
                              axis.title.y = element_text(size = 10)) +
                          ggtitle(paste0(donor, "\nPearson r = ",round(R, 3),", p = ",signif(summary(fit)$coef[2,4], 4),"\nn MAE genes = ",nrow(df))) 
                          #annotate("rect", xmin = 0.25, xmax = 0.75, ymin = 0.25, ymax = 0.75, fill= "#4D004B", alpha = 0.1)
              if(i == 1){
              	R_summary <- c(pair_1, pair_2, round(R, 5), nrow(df))
              }else{
              	R_i <- c(pair_1, pair_2, round(R, 5), nrow(df))
              	R_summary <- rbind(R_summary, R_i)
              }

            pdf(paste0(OUT_dir_i,"AI_corrs/",pair_1,"_vs_",pair_2,".pdf"), height = 4, width = 5)
          	print(plot)
          	dev.off()
      }

      colnames(R_summary) <- c("pair_1","pair_2","R","n_genes")
      R_summary <- as.data.frame(R_summary)
      R_summary$R <- as.numeric(as.character(R_summary$R))
      R_summary$n_genes <- as.numeric(as.character(R_summary$n_genes))
      write.table(R_summary, paste0(OUT_dir_i,"AI_corrs/summary_AI20_correlations_COV",COV,".txt"), quote = FALSE, row.names = FALSE, sep = ",")
}

mean(R_summary$R) ## 0.6949782

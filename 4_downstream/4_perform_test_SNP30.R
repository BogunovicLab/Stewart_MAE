library(tidyverse)
library(plyr)
library(dplyr)
library(reshape2)
library(biomaRt)
library(mgsub)
library(data.table)
library(Qllelic)
library(UpSetR)
library(ggrepel)
ul = function(x){x[1:6, 1:6]}

setwd("/Users/Haley/code/RMAE/240508_COMBINE_NEW")
techReps <- 1

## define filters to consider
fdr_thresh <- 0.05
lower_thresh <- c(0.3, 0.2, 0.1, 0.05)
COV <- 30

types <- paste0(c("both_SNP","missense_SNP","synonymous_SNP"),COV)
keep <- read.table(paste0("results/240519_QC/keep_list.txt"), sep = ",")
keep <- keep$V1

for(j in 1:length(types)){

	type_i <- types[j]
	IN_dir <- paste0("results/240519_QC/",type_i,"/")
	RESULTS_dir <- paste0("results/240520_AI_qllelic/",type_i,"/")
	CUTOFF_dir <- paste0("results/240520_cutoffs_summary/")
	aicovTable <- readRDS(paste0(IN_dir,"aicovTables_meanCOV",COV,".rds"))

	## subset on keep list
	aicovTable_HC <- aicovTable[names(aicovTable) %in% keep]

	## build design matrix for experimental set-up
	designMatrix <- BuildDesign(experimentNames = names(aicovTable_HC), techReps = techReps)

	## set corrConstant to 1 so this can be ignored and the code still runs 
	designMatrix$corrConst <- 1
	## identify allelic imbalanced genes that statistically deviate from the the null point estimate of 0.5 (biallelic)
	AICI_tables <- lapply(1:length(designMatrix$techReps), function(i){
			PerformBinTestAIAnalysisForConditionNPoint_knownCC(inDF = aicovTable_HC[[i]], 
			                                                     vectReps = designMatrix$replicateNums[i], 
			                                                     vectRepsCombsCC = designMatrix$corrConst[i], 
			                                                     thr = COV)
			})
	lapply(AICI_tables, dim)
	names(AICI_tables) <- designMatrix$experimentNames

	## calculate number of statistically significant AI genes present
	sig_genes  <- lapply(1:length(designMatrix$techReps), function(i){nrow(AICI_tables[[i]][!is.na(AICI_tables[[i]]$BT_CC) & AICI_tables[[i]]$BT_CC == T, ])})
	names(sig_genes) <- designMatrix$experimentNames

	## visualize allelic imbalance estimates in the context of mean coverage
	gg_test <- lapply(1:length(designMatrix$experimentNames), function(i){
	cowplot::plot_grid(ggplot(na.omit(AICI_tables[[i]]), aes(sumCOV, AI, col = BT_CC, alpha = BT_CC)) +
			                       geom_point(size=1) +
			                       scale_color_manual(values=c("black","darkorange")) +
			                       scale_alpha_manual(values=c(0.1, 1)) +
			                       theme_bw() +
			                       scale_x_log10() +
			                       xlab("total allelic count") + 
			                       ggtitle(paste0(names(AICI_tables)[i]), 
			                               subtitle = paste(sig_genes[[i]], "biased genes")) +
			                       theme(panel.border = element_rect(colour = "black", fill = NA, size = .75),
			                       		legend.position = "none", 
			                             plot.subtitle = element_text(colour = 'darkorange')))
			})
	pdf(file = paste0(RESULTS_dir,"AI_scatters_meanCOV",COV,".pdf"), width = 5, height = 5)
	print(gg_test)
	dev.off()

	## get genes
	genes <- bind_rows(aicovTable, .id = "sample") 
	genes <- subset(genes, select = c(ID, gene, contig))
	genes <- genes[!duplicated(genes$gene), ]
	colnames(genes)[1] <- "ensID"

	MAE <- AICI_tables %>% lapply(function(x) setNames(x, gsub("^ID", "ensID", names(x))))
	MAE <- bind_rows(MAE, .id = "ID") ## combine all into one df (add clones as rows)
	MAE <- join(MAE, genes, by = "ensID")
	MAE <- MAE[,c(14,15, 1:5,8,9,11,13)]

	## pval histogram by donor
	pval <- ggplot(na.omit(MAE), aes(BT_pval_CC)) +
			    geom_histogram(bins = 50) +
			    facet_wrap(facets = "ID") +
			    theme_classic() +
			    theme(legend.position = "none")
	pdf(file = paste0(RESULTS_dir,"pval_distributions_meanCOV",COV,".pdf"), width = 12, height = 9)
	print(pval)
	dev.off()

	## summary tables
	clones <- names(AICI_tables) 
	for(j in (1:length(clones))){

		clone_i <- clones[j]
		results_i <- subset(MAE, ID %in% clone_i)

		## BH correction
		results_i$fdr_BH <- p.adjust(results_i$BT_pval, method = "BH")
			
		  if(j == 1){
		  	out <- results_i
		  }else{
		  	out <- rbind(out, results_i)
		  }
	}

    out$AI_30 <- ifelse(out$fdr_BH < fdr_thresh & out$AI <= lower_thresh[1] | out$fdr_BH < fdr_thresh & out$AI >= 1 - lower_thresh[1], "MAE", "BAE")
    out$AI_20 <- ifelse(out$fdr_BH < fdr_thresh & out$AI <= lower_thresh[2] | out$fdr_BH < fdr_thresh & out$AI >= 1 - lower_thresh[2], "MAE", "BAE")
    out$AI_10 <- ifelse(out$fdr_BH < fdr_thresh & out$AI <= lower_thresh[3] | out$fdr_BH < fdr_thresh & out$AI >= 1 - lower_thresh[3], "MAE", "BAE")
    out$AI_05 <- ifelse(out$fdr_BH < fdr_thresh & out$AI <= lower_thresh[4] | out$fdr_BH < fdr_thresh & out$AI >= 1 - lower_thresh[4], "MAE", "BAE")

    write.table(out, paste0(CUTOFF_dir,"filtered_SNP/results_summary_SNPmeanCOV",COV,"_",type_i,".txt"), sep = ",", quote = FALSE, row.names = FALSE)
}

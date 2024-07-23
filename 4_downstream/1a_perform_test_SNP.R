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
COV <- 10

types <- c("both","missense","synonymous")

for(j in 1:length(types)){

	if(COV == 10){
		remove <- c("D68-8-1","D68-30","D113-44")
	}else if(COV == 30){
		remove <- c("D68-8-1","D68-30","D113-44","D68-8-3","D113-25")
	}

	type_i <- types[j]
	IN_dir <- paste0("results/240519_QC/unfiltered_SNP/",type_i,"/")
	RESULTS_dir <- paste0("results/240520_AI_qllelic/unfiltered_SNP/",type_i,"/")
	CUTOFF_dir <- paste0("results/240520_cutoffs_summary/unfiltered_SNP/")
	aicovTable <- readRDS(paste0(IN_dir,"aicovTables_meanCOV",COV,".rds"))

	## subset on keep list
	aicovTable_HC <- aicovTable[!names(aicovTable) %in% remove]

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

    write.table(out, paste0(CUTOFF_dir,"results_summary_unfiltered-SNPmeanCOV",COV,"_",type_i,".txt"), sep = ",", quote = FALSE, row.names = FALSE)

    ## plot heat map genes
    if(COV == 10 & type_i == "both"){

		MAE <- aicovTable_HC %>% lapply(function(x) setNames(x, gsub("^ID", "ensID", names(x))))
		MAE <- bind_rows(MAE, .id = "ID") 
		MAE <- separate(MAE, ID, into = c("indiv","clone"), sep = "-", extra = "merge", remove = FALSE)
		MAE$clone <- gsub("Clone-","", MAE$clone)
		genes_to_plot <- c("JAK1", "PLCG2", "NFKB1", "STAT1", "MEFV")
		
		subs <- subset(MAE, gene %in% genes_to_plot)
		subs$label <- ifelse(subs$AI < 0.2 | subs$AI > 0.8, "MAE", "BAE")

		subs$plot_label <- ifelse(subs$label == "BAE", NA, subs$AI)
		pal <- RColorBrewer::brewer.pal(name = "Spectral", n = 11)

		heatmap <- ggplot(data = subs, aes(y = gene, x = clone, fill = AI, col = label, size = log10(meanCOV))) +
			  geom_point(shape = 21, stroke = 1) + 
			  theme_bw() + 
			  scale_fill_gradientn(colours = pal, na.value = "grey80") + 
			  #scale_color_gradient("avg. read\ndepth", low="#BFD3E6", high="#4D004B", trans = "log10",
		      #                               breaks = c(10, 100, 1000, 10000),
		      #                               labels = c(10, 100, 1000, 10000)) + 
			  #scale_size(limits = c(1, 4), range = c(1,8)) +
			  scale_color_manual(values = c("grey50", "#D90166"), na.value = "grey80") + 
			  theme(axis.text.y = element_text(face = "italic", size = 12, hjust = 0),
			        axis.text.x = element_text(size = 12, angle = 60, hjust = 1),
			        panel.grid = element_blank(),
			        axis.title.x = element_blank(),
			        legend.position = "bottom", panel.background = element_rect(color = "black", size = 1), 
			        panel.border = element_rect(size = 0), 
			        strip.background = element_blank(),
			        strip.text.x = element_text(size = 14, color = "black")) + 
			  scale_y_discrete(position = "right") +
			  facet_grid(~ indiv, scales = "free", space = 'free') +
			  guides(fill = guide_colorbar(title.position = "top", direction = "horizontal", frame.colour = "black"),
			         size = guide_legend(title.position = "top", direction = "horizontal"), color = guide_legend(direction = "vertical")) +
			  labs(color = "pattern", size = "read depth (log10)", fill = "allelic imbalance", y = "", x = "")

		pdf(paste0(RESULTS_dir,"heatmap_PID_meanCOV",COV,"_",type_i,".pdf"), width = 12, height = 4.5)
		print(heatmap)
		dev.off()
	}
}

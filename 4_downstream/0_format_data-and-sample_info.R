library(tidyverse)
library(plyr)
library(dplyr)
library(reshape2)
library(biomaRt)
library(mgsub)
library(data.table)
library(Qllelic)
library(UpSetR)
ul = function(x){x[1:6, 1:6]}

setwd("/Users/Haley/code/RMAE/240508_COMBINE_NEW")

IN_dir <- "results/240519_ASE_outs_all/"
OUT_dir <- "results/240519_QC/unfiltered_SNP/"
## will come later, but define number of technical replicates conducted
techReps <- 1
## meanCOV filter
COV <- 30

types <- c("both","missense","synonymous")

for(i in 1:length(types)){

    type_i <- types[i]
    system(paste0("mkdir -p ",OUT_dir,type_i))
    OUT_dir_i <- paste0(OUT_dir,type_i,"/")

    ## list files for gene expression data
    if(type_i == "missense" | type_i == "synonymous"){
        filenames <- list.files(path = paste0(IN_dir,"tables_",type_i), pattern = "*gene*", full.names = TRUE)
        ## SAMPLE INFO
        samples <- as.data.frame(mgsub(filenames, c(paste0("results/240519_ASE_outs_all/tables_",type_i,"/"),"_processed_gene.v3.1.txt"),c("","")))
        colnames(samples) <- "sampleID"
        samples <- separate(samples, col = "sampleID", into = c("donor","clone"), sep = "-", remove = FALSE, extra = "merge")
    }else if(type_i == "both"){
        filenames_1 <- list.files(path = paste0(IN_dir,"tables_missense"), pattern = "*gene*", full.names = TRUE)
        filenames_2 <- list.files(path = paste0(IN_dir,"tables_synonymous"), pattern = "*gene*", full.names = TRUE)
        ## SAMPLE INFO
        samples <- as.data.frame(mgsub(filenames_1, c(paste0("results/240519_ASE_outs_all/tables_missense/"),"_processed_gene.v3.1.txt"),c("","")))
        colnames(samples) <- "sampleID"
        samples <- separate(samples, col = "sampleID", into = c("donor","clone"), sep = "-", remove = FALSE, extra = "merge")
    }

    ## number of clones
    samples <- add_count(samples, donor)
    colnames(samples)[4] <- "n_clones"
    samples$donor <- factor(samples$donor)
    donors <- as.data.frame(summary(samples$donor))
    donors$donors <- rownames(donors)
    colnames(donors)[1] <- "n_clones" 
    donors <- donors[, c(2,1)]

    print(paste0("number of samples: ", nrow(samples)))
    print(paste0("number of donors: ", length(unique(samples$donor))))
    write.table(donors, paste0(OUT_dir_i,"donor_info.txt"), quote = FALSE, row.names = FALSE)

    ## read in bulk gene expression count data
    if(type_i == "missense" | type_i == "synonymous" | type_i == "ALL"){
        rawtables <- lapply(filenames, read.table, header = TRUE)
        experimentNames <- samples$sampleID
        names(rawtables) <- experimentNames
    }else{
        rawtables_1 <- lapply(filenames_1, read.table, header = TRUE)
        rawtables_2 <- lapply(filenames_2, read.table, header = TRUE)

        rawtables <- list()
        for(j in 1:length(rawtables_1)){
            dat_1 <- rawtables_1[[j]]
            dat_2 <- rawtables_2[[j]]
            dat_c <- rbind(dat_1, dat_2)

            rawtables[[j]] <- dat_c
        }
        experimentNames <- samples$sampleID
        names(rawtables) <- experimentNames
    }

    ## change all contigs to character
    rawtables <- lapply(rawtables, function(x) {x$contig <- as.character(x$contig);x})

    ## all genes across samples
    gene_list <- unique(bind_rows(rawtables, .id = "sample")$ID)

    ## get gene names from ensembl
    ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                   mart = ensembl,
                   filters = "ensembl_gene_id",
                   values = gene_list,
                   uniqueRows = T)
    genes[genes == ''] <- NA
    genes <- na.omit(genes)
    genes <- genes[!duplicated(genes$external_gene_name), ]

    ## change to HGNC IDs across all dataframes
    ## save as cleaned Rdata input file
    rawtables_g <- lapply(rawtables, function(x){x$gene <- as.character(genes$external_gene_name[match(x$ID, genes$ensembl_gene_id)]); return(x)})
    rawtables_g <- lapply(rawtables_g, function(x){x[!is.na(x$gene),]})

    ## how many genes per sample?
    gene_info <- as.data.frame(unlist(lapply(rawtables_g, nrow)))
    colnames(gene_info) <- "n_genes_total"
    gene_info$sampleID <- rownames(gene_info)
    samples_g <- join(samples, gene_info, by = "sampleID")

    ## how mnay genes after filtering threshold?
    ## first, calculate mean coverage per gene
    ## if meanCOV < meanCOV filter, value will be set to NA
    get_coverage <- function(counts){MeanCoverage(df = counts, reps = 1, thr = COV)}
    meanCOV <- lapply(rawtables_g, get_coverage)

    ## merge both then remove genes with meanCOV == NA
    rawtables_cov <- list()
    for(j in 1:length(rawtables_g)){

        rawtables_g_j <- rawtables_g[[j]]
        meanCOV_j <- meanCOV[[j]]
        merge_g_j <- cbind(rawtables_g_j, meanCOV_j)
        merge_g_j[,10] <- NULL

        rawtables_cov[[j]] <- merge_g_j
    }
    
    names(rawtables_cov) <- experimentNames
    rawtables_filt <- lapply(rawtables_cov, function(x){x[!is.na(x$meanCOV),]})

    ## number of genes per samples after filtering
    filt_info <- as.data.frame(unlist(lapply(rawtables_filt, nrow)))
    colnames(filt_info) <- "n_genes_filt"
    filt_info$sampleID <- rownames(filt_info)
    samples_filt <- join(samples_g, filt_info, by = "sampleID")
    samples_filt$prop_genes_filt <- round(samples_filt$n_genes_filt/samples_filt$n_genes_total, 3)
    samples_filt <- samples_filt[order(samples_filt$prop_genes_filt),]
    screen <- samples_filt
    screen$donor <- NULL; screen$clone <- NULL; screen$n_clones <- NULL
    screen_nonMAE <- screen
    screen_nonMAE$n_genes_filt <- ifelse(screen_nonMAE$n_genes_filt == screen$n_genes_filt, screen_nonMAE$n_genes_total - screen_nonMAE$n_genes_filt, screen_nonMAE$n_genes_total - screen_nonMAE$n_genes_filt)
    screen$assessable <- paste0("count >= ",COV)
    screen_nonMAE$assessable <- paste0("count < ",COV)
    df <- rbind(screen, screen_nonMAE)
    colnames(filt_info)[1] <- "n_genes_kept"
    write.table(filt_info, paste0(OUT_dir_i,"filter_info.txt"), quote = FALSE, row.names = FALSE)

    ## plot this
    df$assessable <- factor(df$assessable, rev(c(paste0("count >= ",COV), paste0("count < ",COV))))
    new <- c("D108-1","D108-10", "D108-4", "D108-5", "D108-6", "D108-7", "D108-9", "D113-21", "D113-22", "D113-25", "D113-26", "D113-32", "D41-18", "D41-27", "D41-28", "D41-29", "D41-30", "D41-30-2", "D41-31", "D41-32", "D41-33", "D68-17", "D68-27", "D68-30", "D68-33", "D68-8-1", "D68-8-3", "D96F-16", "D96F-18", "D96F-19", "D96F-21", "D96F-3", "D96F-4", "D96F-8")
    df$batch <- ifelse(df$sampleID %in% new, "new batch 2", "batch 1")

    gg_stats <- ggplot(df, aes(fill = assessable, y = n_genes_filt, x = sampleID)) + 
                    geom_bar(position = "stack", stat = "identity") +
                    facet_wrap(vars(batch), scales = 'free') +
                    scale_x_discrete() +
                    scale_y_continuous(expand = c(0,0), position = "left") +
                    scale_fill_manual(values = rev(c("darkorange", "grey70"))) +
                    ylab("number of genes") + 
                    theme_classic() + 
                    coord_flip() +
                    theme(text = element_text(size = 14), 
                          title = element_text(size = 12, face = "bold"),
                          legend.position = "bottom",
                          legend.title = element_blank()) 
    pdf(paste0(OUT_dir_i,"assessable_genes_per_sample_meanCOV",COV,".pdf"), height = 8, width = 15)
    print(gg_stats)
    dev.off()

    ## total number of assessable genes covered across all samples
    totList <- list()
    for(j in 1:length(rawtables_filt)){
        name <- names(rawtables_filt)[j]
        sample <- rawtables_filt[[j]]
        sample <- as.character(subset(sample, select = gene)$gene)
        totList[[j]] <- sample
        names(totList)[[j]] <- name
    }
    ## total number of assessable genes covered across all samples
    assessable <- Reduce(union, totList)
    length(assessable)

    ## summaries
    cat(paste0("avg prop of genes after filtering across samples: ", signif(mean(samples_filt$prop_genes_filt)*100, 3),"%")

    write.table(samples_filt, paste0(OUT_dir_i,"sample_info.txt"), quote = FALSE, row.names = FALSE)
    saveRDS(rawtables_cov, file = paste0(OUT_dir_i,"raw_counts.rds"))
    saveRDS(rawtables_filt, file = paste0(OUT_dir_i,"filtered_meanCOV",COV,"_raw_counts.rds"))
}



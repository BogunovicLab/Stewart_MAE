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
OUT_dir <- "results/240519_QC/"
## meanCOV filter
COV <- 30

types <- c("missense","synonymous")

for(i in 1:length(types)){

    type_i <- types[i]
    OUT_dir_i <- paste0(OUT_dir,type_i,"_SNP",COV,"/")
    system(paste0("mkdir -p ",OUT_dir_i))

    ## list files for gene expression data
    filenames <- list.files(path = paste0(IN_dir,"tables_",type_i), pattern = "*snp*", full.names = TRUE)
    samples <- as.data.frame(mgsub(filenames, c(paste0("results/240519_ASE_outs_all/tables_",type_i,"/"),"_processed_snp.v3.1.txt"),c("","")))
    colnames(samples) <- "sampleID"
    samples <- separate(samples, col = "sampleID", into = c("donor","clone"), sep = "-", remove = FALSE, extra = "merge")

    ## RENAME NEW BATCH INDIVIDUALS
    new <- c("D108-1","D108-10", "D108-4", "D108-5", "D108-6", "D108-7", "D108-9", "D113-21", "D113-22", "D113-25", "D113-26", "D113-32", "D41-18", "D41-27", "D41-28", "D41-29", "D41-30", "D41-30-2", "D41-31", "D41-32", "D41-33", "D68-17", "D68-27", "D68-30", "D68-33", "D68-8-1", "D68-8-3", "D96F-16", "D96F-18", "D96F-19", "D96F-21", "D96F-3", "D96F-4", "D96F-8","D113-36","D113-44")
    samples$new_indiv <- ifelse(samples$sampleID %in% new, "D114", samples$donor)
    samples$new_indiv2 <- ifelse(samples$sampleID %in% new, gsub("D","c",samples$donor), "")
    samples$new_sampleID <- paste0(samples$new_indiv,"-",samples$new_indiv2,"-",samples$clone)
    samples$new_sampleID <- ifelse(!samples$sampleID %in% new, samples$sampleID, samples$new_sampleID)
    colnames(samples)[1] <- "old_sampleID"
    samples$new_indiv2 <- NULL; samples$donor <- NULL
    colnames(samples)[3] <- "donor"
    samples <- samples %>% relocate(new_sampleID)

    ## read in bulk gene expression count data
    rawtables <- lapply(filenames, read.table, header = TRUE)

    ## change all contigs to character
    rawtables <- lapply(rawtables, function(x) {x$contig <- as.character(x$contig);x})
    experimentNames <- samples$new_sampleID
    names(rawtables) <- experimentNames

    rawtables_SNP <- list()
    for(j in 1:length(rawtables)){

        raw_j <- rawtables[[j]]
        genes <- unique(raw_j$groupID)
        out <- c()

        for(z in 1:length(genes)){

            gene_z <- genes[z]

            tab_z <- subset(raw_j, groupID %in% gene_z)
            tab_z$SUM <- rowSums(tab_z[,c(2,3)])
            tab_MAX <- tab_z[which.max(tab_z$SUM),]

            if(as.numeric(tab_MAX$SUM) >= COV){
                out <- rbind(out, tab_MAX)
            }
        }

        if(length(out) == 0){
            out <- data.frame("groupID" = "0", "ref_rep1" = c(0), "alt_rep1" = c(0), "contig" = c(0), "position" = c(0), "SUM" = c(0))
        }

        ## format to output table in *gene* format
        out_f <- subset(out, select = c("groupID","ref_rep1","alt_rep1","contig","position","SUM"))
        names(out_f) <- c("ID","ref_rep1","alt_rep1","contig","aggr_pos","SUM")
        out_f$n_snp <- 1
        out_f$aggr_ref_rep1 <- out_f$ref_rep1
        out_f$aggr_alt_rep1 <- out_f$alt_rep1
        out_f <- out_f[,c("ID","ref_rep1","alt_rep1","n_snp","contig","aggr_pos","aggr_ref_rep1","aggr_alt_rep1","SUM")]

        rawtables_SNP[[j]] <- out_f
        print(names(rawtables)[j])
    }

    names(rawtables_SNP) <- names(rawtables)
    saveRDS(rawtables_SNP, file = paste0(OUT_dir_i,"/SNP_filtered_meanCOV",COV,"_raw_counts.rds"))

    gene_info <- data.frame("sample" = names(rawtables_SNP), n_genes = unlist(lapply(rawtables_SNP, nrow)))
    write.table(gene_info, paste0(OUT_dir_i,"/SNP_filtered_meanCOV",COV,"_gene-info.txt"), quote = FALSE, row.names = FALSE)
}

## using SNP tables, calculate normal QC metrics
types <- c("both","missense","synonymous")

for(i in 1:length(types)){

    type_i <- types[i]
    OUT_dir_i <- paste0(OUT_dir,type_i,"_SNP",COV,"/")
    system(paste0("mkdir -p ",OUT_dir,type_i))

    ## list files for gene expression data
    if(type_i == "missense" | type_i == "synonymous"){
        rawtables <- readRDS(paste0(OUT_dir_i,"/SNP_filtered_meanCOV",COV,"_raw_counts.rds"))
    }else{
        rawtables_1 <- readRDS(paste0("results/240519_QC/missense_SNP",COV,"/SNP_filtered_meanCOV",COV,"_raw_counts.rds"))
        rawtables_2 <- readRDS(paste0("results/240519_QC/synonymous_SNP",COV,"/SNP_filtered_meanCOV",COV,"_raw_counts.rds"))

        rawtables <- list()
        for(j in 1:length(rawtables_1)){
            dat_1 <- rawtables_1[[j]]
            dat_2 <- rawtables_2[[j]]
            dat_c <- rbind(dat_1, dat_2)

            ## if there are two SNPs from the missense and synonymous tables, keep only the one with the most coverage
            if(length(unique(dat_c$ID)) != nrow(dat_c)){

                dup_IDs <- dat_c[duplicated(dat_c$ID), ]$ID

                nondup <- dat_c[!(dat_c$ID %in% dup_IDs),]
                dup <- dat_c[dat_c$ID %in% dup_IDs,]

                for(z in 1:length(dup_IDs)){
                    dup_z <- subset(dup, ID %in% dup_IDs[z])
                    dup_MAX <- dup_z[which.max(dup_z$SUM),]

                    if(z == 1){
                        dup_keep <- dup_MAX
                    }else{
                        dup_keep <- rbind(dup_keep, dup_MAX)
                    }
                }
                dat_c <- rbind(nondup, dup_keep)
            }

            rawtables[[j]] <- dat_c
        }
        names(rawtables) <- names(rawtables_1)
        # saveRDS(rawtables, file = paste0(OUT_dir_i,"/SNP_filtered_meanCOV",COV,"_raw_counts.rds"))
    }

    ## change all contigs to character
    rawtables <- lapply(rawtables, function(x) {x$contig <- as.character(x$contig);x})
    rawtables <- lapply(rawtables, function(x) {x$ID <- as.character(x$ID);x})

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
        merge_g_j[,11] <- NULL

        rawtables_cov[[j]] <- merge_g_j
    }
    
    names(rawtables_cov) <- names(rawtables_g)
    rawtables_filt <- lapply(rawtables_cov, function(x){x[!is.na(x$meanCOV),]})

    ## for control SNP discovery
    controlSNP_df <- bind_rows(rawtables_filt, .id = "sample")
    controlSNP_df <- subset(controlSNP_df, select = c(sample, gene, contig, aggr_pos, SUM))
    sub_df <- subset(controlSNP_df, gene %in% "MX1")


    ## number of genes per samples after filtering
    filt_info <- as.data.frame(unlist(lapply(rawtables_filt, nrow)))
    colnames(filt_info) <- "n_genes_filt"
    filt_info$sampleID <- rownames(filt_info)
    write.table(filt_info, paste0(OUT_dir_i,"filter_info_meanCOV",COV,".txt"), quote = FALSE, row.names = FALSE)

    ## plot this
    gg_stats <- ggplot(filt_info, aes(y = n_genes_filt, x = sampleID)) + 
                    geom_bar(position = "stack", stat = "identity", fill = "darkorange") +
                    scale_x_discrete() +
                    scale_y_continuous(expand = c(0,0), position = "right") +
                    #scale_fill_manual(values = c("darkorange")) +
                    ylab("number of genes") + 
                    theme_classic() + coord_flip() +
                    theme(text = element_text(size = 14), 
                          title = element_text(size = 12, face = "bold"),
                          legend.position = "bottom",
                          legend.title = element_blank())
    pdf(paste0(OUT_dir_i,"assessable_genes_per_sample_meanCOV",COV,".pdf"), height = 10, width = 8)
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

    saveRDS(rawtables_filt, file = paste0(OUT_dir_i,"SNP_filtered_meanCOV",COV,"_raw_counts.rds"))
}

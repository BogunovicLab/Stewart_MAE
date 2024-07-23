library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)
library(patchwork)
library(ggpubr)

setwd("/Users/Haley/results/240614_MiXCR_TCR_check")

tra.chain.files <- list.files("TRA", full.names = T)
trb.chain.files <- list.files("TRB", full.names = T)

process <- function(file, name){
    chains <- lapply(file, read.table, header = TRUE, fill = TRUE, sep = "\t")
    names(chains) <- gsub(file, pattern = ".*/", replacement = "")
    names(chains) <- gsub(names(chains), pattern = ".TRA.txt", replacement = "")
    names(chains) <- gsub(names(chains), pattern = ".TRB.txt", replacement = "")
    combined <- do.call(rbind, chains)
    write.csv(combined, paste0(name,"_combined.csv"), row.names = TRUE)
    return(chains)
}

TRA <- process(tra.chain.files, "TRA")
TRB <- process(trb.chain.files, "TRB")

for(j in 1:length(TRA)){

    name <- names(TRA)[j]
    TRA_i <- TRA[[name]]
    TRB_i <- TRB[[name]]

    p1 <- ggplot(TRA_i, aes(x = length(cloneId), 
                                        y = cloneFraction, 
                                        fill = factor(cloneId))) + 
        geom_bar(stat="identity", width = 0.75, color = "white") +
        scale_y_continuous(expand = c(0,0), position = "left") +
        scale_fill_viridis_d() +
        ylab("Clone fraction") + xlab("") + 
        labs(fill = "TCRa chains IDs") +
        theme_classic() + 
        theme(text = element_text(size = 14), 
              title = element_text(size = 12, face = "bold"),
              axis.text.x.bottom = element_blank(),
              axis.ticks = element_line(linewidth = 1),  # Change size to linewidth
              axis.ticks.x = element_line(linewidth = 0),  # Change size to linewidth
              axis.ticks.y.right = element_blank(),
              legend.position = "none") + 
        ggtitle(label = "TRA")
      
    p2 <- ggplot(TRB_i, aes(x = length(cloneId), 
                                        y = cloneFraction, 
                                        fill = factor(cloneId))) + 
        geom_bar(stat="identity", width = 0.75, color = "white") +
        scale_y_continuous(expand = c(0,0), position = "right") +
        scale_fill_viridis_d() +
        xlab("") +  ylab("") +
        labs(fill = "TCRb chains IDs") +
        theme_classic() + 
        theme(text = element_text(size = 14), 
              title = element_text(size = 12, face = "bold"),
              axis.text.x.bottom = element_blank(),
              axis.text.y.left = element_blank(),
              axis.ticks = element_line(linewidth = 1),  # Change size to linewidth
              axis.ticks.x = element_line(linewidth = 0),  # Change size to linewidth
              axis.ticks.y.left = element_blank(),
              legend.position = "none") + 
        ggtitle(label = "TRB")
      
    plots <- ggarrange(p1, p2, ncol = 2, nrow = 1)
    ggsave(paste0("cloneFraction/",name,"cloneFraction_stackedbars.pdf"), plots, device = "pdf", width = 4, height = 6)
}



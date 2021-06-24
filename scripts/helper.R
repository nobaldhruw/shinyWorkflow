library(limma)
library(enrichR)
library(ggplot2)
library(ggrepel)

#####################################################################
#
#  Compute PCA
#
#####################################################################
compute_pca <- function(mat, metadata, num_highly_variable_genes = NULL){
    print("Compute PCA called")
    # Drop rows with NA values
    mat <- as.data.frame(mat[complete.cases(mat), ])
    
    # Drop rows with zero variance
    row_variance <- apply(mat, 1, var)
    mat <- mat[row_variance > 0, ]
    
    if(!is.null(num_highly_variable_genes)){
        mat$mad <- apply(mat, 1, mad)
        mat <- mat[order(mat$mad, decreasing = TRUE), ]
        mat <- mat[1:num_highly_variable_genes, ]
        mat$mad <- NULL
    }
    
    ## compute PCA
    pca_result <- prcomp(t(mat), center = T, scale = T)
    PCAScores <- data.frame(pca_result$x, metadata)
  
    return (PCAScores) 
}
#####################################################################
#
#  Plot PCA scores
#
#####################################################################
pca_plot <- function(PCAScores, cohort, PCx = 1, PCy = 2){
    require(ggplot2)
    p <- ggplot(PCAScores, aes_string(x = paste0('PC', PCx),
                                    y = paste0('PC', PCy),
                                    fill = cohort)) + 
    geom_point(shape = 21, size = 7, alpha = 0.7) + 
    labs(x = paste("PC",PCx),
         y = paste("PC",PCy), fill = "cohort") + 
    theme(legend.position = "right", legend.direction = "vertical",
          axis.line = element_line(size=1, colour = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), panel.background = element_blank(), 
          axis.title = element_text(colour="black", size = 15, face = "bold"), 
          axis.text.x = element_text(colour="black", size = 10, margin=unit(c(0.5,0.5,0.1,0.1), "cm"), face = "bold"), 
          axis.text.y = element_text(colour="black", size = 10, margin=unit(c(0.5,0.5,0.1,0.1), "cm"), face = "bold"), 
          legend.text = element_text(size = 10, face = "bold"),
          legend.title = element_text(colour="black", size=12, face="bold"),
          axis.ticks.length = unit(-0.25, "cm"))
    p
}
#####################################################################
#
#  Limma's differential expression
#
#####################################################################
diff_exp_limma <- function(data_matrix, metadata, cohortCol, cohort_a, cohort_b, p_val_correct_methods='BH') {
    
    # Subset the metadata dataframe to have samples from the provided cohorts only
    metadata <- metadata[metadata[, cohortCol] %in% c(cohort_a, cohort_b),]
    # Create a comparison column
    metadata[, "comparison"] <- NA
    metadata[metadata[, cohortCol] %in% cohort_a, "comparison"] <- "A"
    metadata[metadata[, cohortCol] %in% cohort_b, "comparison"] <- "B"
    
    condition <- metadata[, "comparison"]
    ### Creating design matrix of comparisons required (without batch)
    design <- model.matrix(~ condition + 0)
    colnames(design) <- gsub("condition", "", colnames(design))
    contrast_matrix <- limma::makeContrasts(contrasts = c(paste("B", "A", sep = "-")), levels = design)

    data_matrix <- data_matrix[, rownames(metadata)]
    data_matrix <- data_matrix[!apply(data_matrix, 1, anyNA), ]

    ### Fitting linear model to log2 normalised expression data
    fit <- limma::lmFit(data_matrix, design)
    fit <- limma::contrasts.fit(fit, contrast_matrix)
    fit <- limma::eBayes(fit)
    limma_results_df <- limma::topTable(fit, coef = paste("B", "A", sep = "-"), number = nrow(data_matrix))

    ### p value correction
    for (p_val_correct in p_val_correct_methods) {
        if (p_val_correct == "Bonferroni") {
            p_val_correct_method <- "bonferroni"
        }else {
            p_val_correct_method <- "BH"
        }    
        limma_results_df[, p_val_correct] <- p.adjust(limma_results_df$P.Value, method = p_val_correct_method)
    }
    return(limma_results_df)
}
#####################################################################
#
#  Volcano plot
#
#####################################################################
volcano_plot <- function(de_result, log2fc_cutoff=0, p_val_cutoff=0.05, ngenes_to_label = 10) {
    require(ggplot2)
    require(ggrepel)
    # Add logical vector as a column (significant) to the res_tableOE
    de_result$significance <- "not significant"
    de_result[(abs(de_result$logFC) > log2fc_cutoff) & (de_result$P.Value <= p_val_cutoff), "significance"] <- "significant"

    # Select top N genes for labelling 
    ## Sort by ordered padj
    de_result <- de_result[order(de_result$adj.P.Val), ] 

    ## Create a column to indicate which genes to label
    de_result$genelabels <- ""
    de_result$genelabels[1:ngenes_to_label] <- rownames(de_result)[1:ngenes_to_label]

    ## Volcano plot
    p <- ggplot(de_result, aes(x = logFC, y = -log10(adj.P.Val), fill = significance,label=genelabels)) +
        geom_point(shape = 21, size = 3, alpha = 1.0) +
        geom_text_repel(size = 5) + 
        #geom_text(aes(label = ifelse(genelabels !=0, as.character(genelabels), "")),hjust = 0, nudge_x = 0.05, size=5, colour = "black") +
        xlab("log2FC") + 
        ylab("-log10 (padj)") +
        theme(legend.position = "top", legend.direction = "horizontal", # legend positioned at the bottom, horizantal direction,
            axis.line = element_line(size=1, colour = "black"),	# axis line of size 1 inch in black color
            panel.grid.major = element_blank(),	# major grids included
            panel.grid.minor = element_blank(),	# no minor grids
            panel.border = element_blank(), panel.background = element_blank(), # no borders and background color
            axis.title = element_text(colour="black", size = 25, face = "bold"), # axis title 
            axis.text.x = element_text(colour="black", size = 20, margin=unit(c(0.5,0.5,0.1,0.1), "cm"), face = "bold"), # x-axis text in fontsize 10
            axis.text.y = element_text(colour="black", size = 20, margin=unit(c(0.5,0.5,0.1,0.1), "cm"), face = "bold"), # y-axis text in fontsize 10
            legend.text = element_text(size = 15, face = "bold"),
            legend.title = element_blank(),
            axis.ticks.length = unit(-0.25, "cm")) # ticks facing inward with 0.25cm length
    p
}
#####################################################################
#
#  Fetch dbList from enrichR
#
#####################################################################
get_enrichr_dblist <- function(){
    return(enrichR::listEnrichrDbs()[,'libraryName'])
}
#####################################################################
#
#  Identify enriched gene sets
#
#####################################################################
get_enriched_gene_sets <- function(de_result, 
                                   gene_direction = 'both', 
                                   log2fc_cutoff  = 0, 
                                   fdr_cutoff     = 0.05,
                                   db             = 'KEGG_2019_Human')
{
    require(enrichR)
    gene_set <- row.names(de_result[(de_result[,'adj.P.Val'] < fdr_cutoff) & (abs(de_result$logFC) > log2fc_cutoff), ])
    if(gene_direction == 'up'){
        gene_set <- row.names(de_result[(de_result[,'adj.P.Val'] < fdr_cutoff) & (de_result$logFC > log2fc_cutoff), ])
    }else if(gene_direction == 'down'){
        gene_set <- row.names(de_result[(de_result[,'adj.P.Val'] < fdr_cutoff) & (de_result$logFC < -log2fc_cutoff), ])
    }
    enriched <- enrichr(gene_set, db)
    enriched[[db]]$geneRatio <- sapply(enriched[[db]]$Overlap, function(x){
        num <- as.numeric(strsplit(x, '/')[[1]][1])
        den <- as.numeric(strsplit(x, '/')[[1]][2])
        return (num/den)
    })
    egs <- enriched[[db]]
    egs <- egs[order(egs$Combined.Score, decreasing = T), ]
    return(egs)
}
#####################################################################
#
#  Visualise enriched gene sets
#
#####################################################################
options(repr.plot.width = 12, repr.plot.height = 6)

plot_enriched_gene_sets <- function(egs, nsets = 10){
    p <- ggplot(egs[1:nsets, ], aes(x = Combined.Score, y = Term, fill = Adjusted.P.value)) +
        geom_point(shape = 21, alpha = 1.0, aes(size = geneRatio)) +
        scale_y_discrete(limits = rev(egs[1:nsets,'Term'])) +
        scale_fill_continuous(low="brown", high="grey", limits=c(0,0.05)) +
        scale_size(range = c(2, 10)) +
        xlab("Combined Score") + 
        ylab("Gene set") +
        theme(legend.position = "right", legend.direction = "vertical", # legend positioned at the bottom, horizantal direction,
            axis.line = element_line(size=1, colour = "black"),	# axis line of size 1 inch in black color
            panel.background = element_rect(fill = "white",
                                            colour = "black",
                                            linetype = "solid", size = 1),
            panel.grid.major = element_line(linetype = 'solid', colour = "grey"),
            axis.title = element_text(colour="black", size = 20, face = "bold"), # axis title 
            axis.text.x = element_text(colour="black", size = 15, margin=unit(c(0.5,0.5,0.1,0.1), "cm"), face = "bold"), # x-axis text in fontsize 10
            axis.text.y = element_text(colour="black", size = 15, margin=unit(c(0.5,0.5,0.1,0.1), "cm"), face = "bold"), # y-axis text in fontsize 10
            legend.text = element_text(size = 15, face = "bold"),
            legend.title = element_text(colour = "black", size = 15, face = "bold"),
            axis.ticks.length = unit(-0.25, "cm")) # ticks facing inward with 0.25cm length
    p
}
#####################################################################
#
#  Run X2K
#
#####################################################################
run_x2k <- function(de_result, 
                    gene_direction = 'both', 
                    log2fc_cutoff  = 0, 
                    fdr_cutoff     = 0.05,
                    dbs            = c('ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X','KEA_2015'))
{
    require(enrichR)
    gene_set <- row.names(de_result[(de_result[,'adj.P.Val'] < fdr_cutoff) & (abs(de_result$logFC) > log2fc_cutoff), ])
    if(gene_direction == 'up'){
        gene_set <- row.names(de_result[(de_result[,'adj.P.Val'] < fdr_cutoff) & (de_result$logFC > log2fc_cutoff), ])
    }else if(gene_direction == 'down'){
        gene_set <- row.names(de_result[(de_result[,'adj.P.Val'] < fdr_cutoff) & (de_result$logFC < -log2fc_cutoff), ])
    }
    enriched <- enrichr(gene_set, dbs)    
    return(enriched)
}
#####################################################################
#
#  plot X2K
#
#####################################################################
plot_x2k_resuts <- function(x2k_result, nbars = 10){
    require(ggplot2)
    #TF
    plot_data <- x2k_result[[1]]
    plot_data <- plot_data[order(plot_data$`Combined.Score`, decreasing = T), ]
    plot_data <- plot_data[1:nbars, ]
    plot_data$x <- plot_data$Term
    plot_data$y <- plot_data$`Combined.Score`
    
    p1 <- ggplot(plot_data, aes_string(x = 'Term', y = 'Combined.Score', fill = 'Adjusted.P.value')) +
    geom_bar(stat = "identity", color = "black", position = position_dodge(), width = 0.7) +
    scale_x_discrete(limits = rev(plot_data$x)) + 
    scale_fill_continuous(low="brown", high="grey", limits=c(0,0.05)) +
    theme(legend.position = "right", legend.direction = "vertical") +
    labs(x = "Transcription Factors", y = "Combined score") +
    ggtitle(element_blank()) +
    theme(
      axis.line = element_line(size = 1, colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    ) +
    theme(
      plot.title = element_blank(),
      text = element_text(family = "Tahoma"),
      axis.text.x = element_text(colour = "black", size = 15, face = "bold", angle = 90, margin = unit(c(0.7, 0.7, 0.1, 0.1), "cm")),
      axis.text.y = element_text(colour = "black", size = 15, face = "bold", margin = unit(c(0.5, 0.5, 0.1, 0.1), "cm")),
      legend.key = element_rect(fill = "white", colour = "white"),
      axis.title = element_text(colour = "black", size = 20, face = "bold"),
      legend.text = element_text(size = 15, face = "bold"),
      legend.title = element_text(colour = "black", size = 15, face = "bold"),
      axis.ticks.length = unit(-0.25, "cm"),
      axis.ticks.x = element_blank()
    ) +
    coord_flip()
    
    # Kinase
    plot_data <- x2k_result[[2]]
    plot_data <- plot_data[order(plot_data$`P.value`, decreasing = F), ]
    plot_data <- plot_data[1:nbars, ]
    plot_data <- plot_data[order(plot_data$`Combined.Score`, decreasing = T), ]
    plot_data$x <- plot_data$Term
    plot_data$y <- plot_data$`Combined.Score`
    
    p2 <- ggplot(plot_data, aes_string(x = 'Term', y = 'Combined.Score', fill = 'Adjusted.P.value')) +
    geom_bar(stat = "identity", color = "black", position = position_dodge(), width = 0.7) +
    scale_x_discrete(limits = rev(plot_data$x)) + 
    scale_fill_continuous(low="brown", high="grey", limits=c(0,0.05)) +
    theme(legend.position = "right", legend.direction = "vertical") +
    labs(x = "Kinases", y = "Combined score") +
    ggtitle(element_blank()) +
    theme(
      axis.line = element_line(size = 1, colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    ) +
    theme(
      plot.title = element_blank(),
      text = element_text(family = "Tahoma"),
      axis.text.x = element_text(colour = "black", size = 15, face = "bold", angle = 90, margin = unit(c(0.7, 0.7, 0.1, 0.1), "cm")),
      axis.text.y = element_text(colour = "black", size = 15, face = "bold", margin = unit(c(0.5, 0.5, 0.1, 0.1), "cm")),
      legend.key = element_rect(fill = "white", colour = "white"),
      axis.title = element_text(colour = "black", size = 20, face = "bold"),
      legend.text = element_text(size = 15, face = "bold"),
      legend.title = element_text(colour = "black", size = 15, face = "bold"),
      axis.ticks.length = unit(-0.25, "cm"),
      axis.ticks.x = element_blank()
    ) +
    coord_flip()

    gridExtra::grid.arrange(p1, p2, nrow=1)
}
## Function for computing PCA
compute_pca <- function(input_matrix, metadata, num = NULL, center = TRUE, scale = TRUE) {
  require(stats)
  
  input_matrix <- input_matrix[!apply(input_matrix, 1, anyNA), ]
  
  variancerow <- as.numeric(apply(input_matrix, 1, function(x) var(x)))
  input_matrix <- input_matrix[!(variancerow == 0),]  
  
  if (!identical(num, NULL)){
    temp_df <- apply(input_matrix, 1, function(x) mad(x))
    input_matrix <- as.data.frame(cbind(input_matrix,temp_df))
    input_matrix <- input_matrix[order(-input_matrix$temp_df),]
    data_mat <- as.data.frame(t(head(input_matrix[, !(colnames(input_matrix) %in% c("temp_df"))], num)))
  }else{
    data_mat <- as.data.frame(t(input_matrix))
  }
  
  pca_result <- stats::prcomp(data_mat, center = center, scale = scale)
  PCAScores <- data.frame(pca_result$x, metadata)
  
  return (PCAScores)
}

## PCA plotting function
pca_plot <- function(PCAScores, cohort, PCx = 1, PCy = 2){

  p <- ggplot(PCAScores, aes_string(x = paste0('PC', PCx),
                                    y = paste0('PC', PCy),
                                    fill = cohort)) + # calls the ggplot function with dose on the x-axis and len on the y-axis
    geom_point(shape = 21, size = 4, alpha = 0.7) + # scatter plot function with shape of points defined as 21 scale.
    labs(x = paste("PC",PCx),
         y = paste("PC",PCy), fill = "cohort") + # x and y axis labels
    theme(legend.position = "right", legend.direction = "vertical", # legend positioned at the bottom, horizantal direction,
          axis.line = element_line(size=1, colour = "black"), # axis line of size 1 inch in black color
          panel.grid.major = element_blank(), # major grids included
          panel.grid.minor = element_blank(), # no minor grids
          panel.border = element_blank(), panel.background = element_blank(), # no borders and background color
          axis.title = element_text(colour="black", size = 15, face = "bold"), # axis title
          axis.text.x = element_text(colour="black", size = 10, margin=unit(c(0.5,0.5,0.1,0.1), "cm"), face = "bold"), # x-axis text in fontsize 10
          axis.text.y = element_text(colour="black", size = 10, margin=unit(c(0.5,0.5,0.1,0.1), "cm"), face = "bold"), # y-axis text in fontsize 10
          legend.text = element_text(size = 10, face = "bold"),
          legend.title = element_text(colour="black", size=12, face="bold"),
          axis.ticks.length = unit(-0.25, "cm"))
  p
}
# Function to make Heatmap with cell ratio
# Further ideas:
#   Possibility to modify further pHeatmap parameters directly from calling function
#   Some multiplot feature (as for Figure 4F)

pl_cell_frac_pheatmap_v2 <- function(
    object,
    column_data = "Timepoint",
    row_data = "CellType",
    filename_save = FALSE,
    file_height = 6,
    file_width = 5,
    row_selection = c(1:length(levels(object@meta.data[[row_data]]))),
    column_selection = c(1:length(levels(object@meta.data[[column_data]]))),
    threshold_value = 0,
    NA_tiles = "none",
    include_absolute_values = FALSE,
    ratio = "global"){
  
  # Prepare the data for the heatmap
  heatmap_prop_pre = data.frame(object@meta.data[[row_data]],
                                object@meta.data[[column_data]])
  absolute_numbers <- table(heatmap_prop_pre)[row_selection, column_selection]
  
  # Define margin parameter for prop.table based on ratio parameter of function
  ratio_dict_list = list("global" = NULL, "row" = 1, "column" = 2)
  
  # Extract data of interest from heatmap_prop_pre and calculate ratio
  heatmap_prop <- prop.table(table(heatmap_prop_pre), margin = unlist(ratio_dict_list[ratio]))[row_selection,column_selection]
  heatmap_prop[is.na(heatmap_prop) == TRUE] = 0
  
  # Prepare another dataframe based on the heatmap data, but only showing numbers differing from 0
  heatmap_prop_num <- heatmap_prop
  heatmap_prop_num[heatmap_prop_num <= threshold_value] <- -1 # Easy workaround for NAs not being rounded
  heatmap_prop_num <- round(heatmap_prop_num, digits = 2)
  heatmap_prop_num[heatmap_prop_num == 0] <- '<0.01'
  heatmap_prop_num[heatmap_prop_num == -1] <- ''
  
  # Combine ration numbers and absolute values
  if (include_absolute_values == TRUE){
    for (i in 1:length(rownames(heatmap_prop_num))){
      for (j in 1:length(colnames(heatmap_prop_num))){
        heatmap_prop_num[i,j] <- paste(heatmap_prop_num[i,j], absolute_numbers[i,j], sep = "\n")
      }
    }
    heatmap_prop_num[heatmap_prop_num == '\n0'] <- ''
  }
  
  # Define another Dataframe for number colour
  max_heatmap_prop <- max(heatmap_prop, na.rm = TRUE)
  color_number <- heatmap_prop>(max(heatmap_prop, na.rm = TRUE)-min(heatmap_prop, na.rm = TRUE))/2
  color_number[color_number == TRUE] <- "black"
  color_number[color_number == FALSE] <- "white"
  
  # Make not possible tiles NA, to have them in different colour as 0 values
  if (NA_tiles[1] != "none"){
    heatmap_prop[NA_tiles == TRUE] = NA
    # heatmap_prop[NA_tiles == FALSE && is.na(heatmap_prop) == TRUE] = 0
  }
  
  # Prepare fontsize number paramter for pheatmap, dependent on include absolute values and occurence of <
  if(include_absolute_values==TRUE){
    number_fontsize = 12}else if('<0.01' %in% heatmap_prop_num){
      number_fontsize = 12}else{
        number_fontsize = 15}
  
  # Plot pHeatmap
  p1 <- pheatmap(heatmap_prop, scale = "none", cluster_rows = FALSE,
           cluster_cols = FALSE, col=c("#DDDDDD", viridis::cividis(101)), 
           breaks = c(0,threshold_value+0.00001,
                      seq(threshold_value+0.00002, 
                          max_heatmap_prop,length.out = 98)),
           display_numbers =  heatmap_prop_num, number_color = color_number,
           cellwidth = 30, cellheight = 30, 
           fontsize_number = number_fontsize, 
           main = paste("number of cells: ", sum(table(heatmap_prop_pre))), 
           na_col = "white")
  
  # Save if filename_save is give, otherwise just print
  if (filename_save != FALSE){
    pdf(filename_save, height = file_height, width = file_width)
    print(p1)
    dev.off()
  } else {print(p1)}
}


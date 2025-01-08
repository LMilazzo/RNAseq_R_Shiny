
library(readr)
library(dplyr)
library(shiny)



##___________________UI FUNCTIONS________________________

# Basic function that returns and displays a modal with the given message in red text
showErrorModal <- function(message) {
  showModal(modalDialog(tags$p(style = "color: red;", message), easyClose = TRUE, footer = NULL))
}



#________________FILTERING FUNCTIONS___________________________

# Splits a DESeq result dataset into 3 dataframes based on fold change provided
# cut is a vector of two values the lower cutoff and the upper cutoff
splitByExpr <- function(data, cut){
  
  lowerBound <- cut[1]
  upperBound <- cut[2]
  
  data$gene_name <- rownames(data)

  res <- as.data.frame(data)

  up <- res %>% filter(log2FoldChange > upperBound)
  up <- up[order(up$padj) , ]
  
  down <- res %>% filter(log2FoldChange < lowerBound)
  down <- down[order(down$padj) , ]
  
  noR <- res %>% 
         filter(log2FoldChange > lowerBound) %>% 
         filter(log2FoldChange < upperBound)
  
  noR <- noR[order(noR$padj) , ]
  
  return( list( up, down, noR ) )
  
}

# Filters the raw counts upload
# The row sum cutoff for a gene to be kept is hard coded as 10 this can and 
# should be changed within this function
filterCounts <- function(counts){
  
  keep <- rowSums(counts) > 10
  
  counts <- counts[keep,]
  
  counts <- counts
  
  return(counts)
}


#________________PLOTTING FUNTIONS___________________________

# Returns a ggplot of the principle components
principlePlot <- function(data, cond, genes, title, sub, cap){
  
  pca_data <- plotPCA(data, intgroup=cond, ntop=genes, returnData=TRUE)
      
      percentVar <- round(100 * attr(pca_data, "percentVar"))
      
      aesthetics <- aes(x = pca_data$PC1, y = pca_data$PC2 )
      
      x <- length(cond)
    
      if(x >= 1){
        aesthetics$colour = as.name(cond[1])
      }
      if(x >= 2){
        aesthetics$shape = as.name(cond[2])
      }
      if(x >= 3){
        aesthetics$size = as.name(cond[3])
      }
      
      pca_plot <- ggplot(pca_data,aesthetics)+
        
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        
        theme_minimal() +

        labs(title = title,
            subtitle = sub,
            caption = cap,
            color = if (x >= 1) cond[1] else NULL,
            shape = if (x >= 2) cond[2] else NULL,
            size = if (x >= 3) cond[3] else NULL
            ) +

        theme(plot.margin = margin(10, 10, 10, 10, "pt"),
              axis.title.x = element_text(color='black', size = 20, margin = margin(15, 15, 15, 15, "pt")),
              axis.title.y = element_text(color='black', size=20, margin = margin(15, 15, 15, 15, "pt")),
              axis.text.y = element_text(color='black', size=15),
              axis.text.x = element_text(color='black', size=15),
              panel.grid.minor.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              axis.line.x = element_line(color='grey'),
              axis.line.y = element_line(color='grey'),

              legend.text = element_text(color='black', size=20),
              legend.title = element_text(color='black', size=20),
              legend.margin = margin(15, 15, 15, 15, "pt"),

              plot.title = element_text(color='black', size=30, margin=margin(20,20,5,10,"pt")),
              plot.subtitle = element_text(color='black', size=20, margin=margin(5,5,15,10,"pt")),
              plot.caption = element_text(color='black', size=15, margin=margin(10, 10, 10, 10, "pt")),

              aspect.ratio = 1) +
        coord_fixed()

      #Change size of points according to present number of variables
      if(x >= 3){
        pca_plot <- pca_plot + scale_size_discrete(range = c(3, 10)) + geom_point()
      }else{
        pca_plot <- pca_plot + geom_point(size=4)
      }

      return(pca_plot)
}

##_____________EXPERIMENT FUNCTIONS______________________________

# Function That runs pathfindR
runPathfindRFunc <- function(data_source){
  
  #Set up database
  kegg <- plyr::ldply(pathfindR.data::kegg_genes, data.frame) %>% 
    mutate(num_genes_in_path = 1) %>% 
    dplyr::group_by(.data$.id) %>% 
    dplyr::summarize(X..i.. = paste0(.data$X..i.., collapse = ", "), num_genes_in_path = sum(num_genes_in_path)) 
  names(kegg)[1] <- "ID"
  names(kegg)[2] <- "all_pathway_genes"
  
  showModal(modalDialog("Finding Paths...", footer = NULL))
  #Run pathfinder
  
  tryCatch({
    
    res <- run_pathfindR(data_source,
                         pin_name_path = "KEGG",
                         enrichment_threshold = 0.05,
                         iterations = 25,
                         list_active_snw_genes = TRUE)
    
    if(ncol(res) == 0 || nrow(res) == 0){
      showErrorModal('No Enriched Terms were found for the provided pin')
      return()
    }
    
    if(length(intersect(colnames(res), colnames(kegg))) == 0){
      showErrorModal("No common variables in pin and pathfinder output")  
      return()
    }
    
    res <- left_join(res, kegg)
    
    res_clustered <- cluster_enriched_terms(res,  
                                            plot_dend = FALSE, 
                                            plot_clusters_graph = FALSE)
    
    showModal(modalDialog("pathfindR complete with no fatal errors", easyClose = TRUE, footer = NULL))
    
    return(res_clustered)
    
  },
  error = function(e) {
    
    showErrorModal(paste("Error in pathfindR process:", e$message))
    
  })
  
}


#_____________FILE READING FUNCTIONS______________________________

# Reads a DESeq2 output file csv
#Returns (results, gene names, normalized counts, vst counts, vst object)
readDEGFile <- function(fileDataPath, meta_data){
  
  #                         READ FILE
  #_______________________________________________________________
  
  #Read the file as tsv
  if(grepl('\\.tsv$', fileDataPath, ignore.case = TRUE)){
    data <- read.csv(fileDataPath, sep = "\t")
  }
  #Read the file as a csv
  else if(grepl('\\.csv$', fileDataPath, ignore.case = TRUE)){
    data <- read.csv(fileDataPath)
  }
  #Unreadable file type
  else{
    showErrorModal('The DEG file is not a readable file')
    return()
  }
  
  #                   CHECK FOR REQ COLUMNS
  #______________________________________________________________
  
  req_col <- c( "gene_id" , "gene_name" , "log2FoldChange" , "padj" )
  possible_col <- c("gene_id", "gene_name", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  
  if(! all( req_col %in% colnames(data) ) ){
    showErrorModal("DEG file missing either ('gene_id', 'gene_name', 'log2FoldChange', 'padj')")
    return()
  }
  
  #                   CHECK DUPLICATE IDS
  #_______________________________________________________________
  
  #Check duplicate ids
  if (anyDuplicated(data$gene_id)) {
    showErrorModal("Gene ids need to be unique for each row")
    return()
  }
  
  #                   READ THE DATA (exluding counts)
  #_______________________________________________________________
  
  #select all existent possible columns
  results <- data %>%
    select(all_of(intersect(colnames(data), possible_col)))
  
  #set rownames
  rownames(results) <- results$gene_name
  
  #get Gene Names
  gene_names <- rownames(results)
  
  #                   PROCESSING SAMPLE COLUMNS
  #_______________________________________________________________
  
  sample_cols <- data %>% select(starts_with('.'), gene_name)
  rownames(sample_cols) <- sample_cols$gene_name
  sample_cols <- sample_cols %>% select(-gene_name)
  
  #RETURN 1 NO COUNTS RETURN
  if (ncol(sample_cols) == 0) {
    
    showErrorModal("Warning: Some features not available because of missing sample data")
    
    return(list(results, gene_names, NULL, NULL, NULL))
    
  }
  
  #ERROR 1
  if (nrow(sample_cols) != nrow(results)) {
    showErrorModal("Error: number of rows in the sample data != rows in results data")
    return()
  }
  
  #ERROR 2
  if( ncol(sample_cols) != nrow(meta_data) ){
    showErrorModal("Error: Number of samples in meta data is not the same as the number of samples in the counts data")
    return()
  }
  
  #                 PROCESSING OBJECTS 
  #_______________________________________________________________
  counts <- sample_cols %>% as.matrix() %>% round()   
  
  v <- varianceStabilizingTransformation(counts, blind=TRUE)
  
  se <- SummarizedExperiment(assays = list(counts = v), colData = meta_data)
  
  vsd <- DESeqTransform(se)
  
  
  return(list(results, gene_names, sample_cols, data.frame(v), vsd))
  
}


# Reads a meta data file
readMetaData <- function(fileDataPath){
  
  #                         READ FILE
  #_______________________________________________________________
  
  #Read the file as tsv
  if(grepl('\\.tsv$', fileDataPath, ignore.case = TRUE)){
    md <- read.csv(fileDataPath, sep = "\t")
  }
  #Read the file as a csv
  else if(grepl('\\.csv$', fileDataPath, ignore.case = TRUE)){
    md <- read.csv(fileDataPath)
  }
  #Unreadable file type
  else{
    showErrorModal('The meta data file is not a readable file')
    return()
  }
  
  #               CHECK FOR NEEDED COLUMNS
  #_____________________________________________________________
  
  # Error: Insufficient columns
  if ( ncol(md) < 2 || ! "sample" %in% colnames(md)) {
    showErrorModal("The meta data file is missing columns")
    return()
  }
  
  #               ASSERT SAMPLE NAMING CONVENTION
  #_____________________________________________________________
  md$sample <- ifelse(!grepl("^\\.", md$sample), paste0(".", md$sample), md$sample)
  md <- md %>%
    tibble::column_to_rownames('sample')
  
  md[] <- lapply(md, as.factor)
  
  return(md)
  
}


# Reads counts to be part of a pathfinder experiences as abundance data
readPathfinderCountsUpload <- function(fileDataPath){
  
  #                         READ FILE
  #_______________________________________________________________
  
  #Read the file as tsv
  if(grepl('\\.tsv$', fileDataPath, ignore.case = TRUE)){
    counts <- read.csv(fileDataPath, sep = "\t")
  }
  #Read the file as a csv
  else if(grepl('\\.csv$', fileDataPath, ignore.case = TRUE)){
    counts <- read.csv(fileDataPath)
  }
  #Unreadable file type
  else{
    showErrorModal('The counts matrix file is not a readable file')
    return()
  }
  
  #                   CHECK REQUIRED COLUMS
  #_______________________________________________________________
  
  #Check for req columns
  if(! 'gene_id' %in% colnames(counts) ){
    showErrorModal('Required columns gene_id is missing in the counts matrix')
    return()
  }
  if(! 'gene_name' %in% colnames(counts) && !'Gene_symbol' %in% colnames(counts)){
    showErrorModal('Required columns Gene_symbol is missing in the counts matrix')
    return()
  }
  
  #                   CHECK DUPLICATE IDS
  #_______________________________________________________________
  
  #Check duplicate ids
  if(anyDuplicated(counts$gene_id)){
    showErrorModal('There were duplicate ids in the counts matrix')
    return()
  }
  
  #                   FIND SAMPLE COLUMN NAMES
  #_______________________________________________________________
  sample_cols <- counts %>% select(starts_with('.'))
  if(ncol(sample_cols) < 2){
    showErrorModal('There were not enough samples present')
    return()
  }
  
  #                   FIND GENE NAMES COLUMN
  #_______________________________________________________________
  pos <- which(names(counts) %in% c('gene_name', 'Gene_symbol'))
  pos <- pos[1]
  
  #                     BUILD DATAFRAME
  #_______________________________________________________________
  ab <- sample_cols
  ab$Gene_symbol <- counts[, pos]
  ab <- ab %>% select(Gene_symbol, everything())
  
  #               CREATE A DATAFRAME OF GENE NAMES
  #_______________________________________________________________
  gene_names <- ab$Gene_symbol
  
  #RETURN
  return(list(ab, gene_names))
  
}


# Reads a pathfindR output file
readPathfinderOutput <- function(fileDataPath){
  
  if(!is.null(fileDataPath)){
    
    path <- fileDataPath
    
    if (!grepl('\\.(csv|tsv)$', path, ignore.case = TRUE)) {
      showErrorModal('An unreadable file type was submitted (use csv or tsv)')
      return()
    }
    
    if (grepl('\\.tsv$', path, ignore.case = TRUE)) {
      new_data <- read.csv(path, sep = "\t")
    } else {
      new_data <- read.csv(path)
    }
    
    req_cols <- c("ID","Term_Description","Fold_Enrichment","occurrence",
                  "support","lowest_p","highest_p","non_Signif_Snw_Genes",
                  "Up_regulated","Down_regulated","all_pathway_genes",
                  "num_genes_in_path","Cluster","Status")
    
    if(!all(req_cols %in% colnames(new_data))){
      showErrorModal('There are some missing cols in the uploaded set')
      return()
    }
    
    #Assert correct classes
    original_na_counts <- sapply(colnames(new_data), 
                                 function(col) sum(is.na(new_data[[col]]))) 
    
    new_data$ID <- as.character(new_data$ID)
    new_data$Term_Description <- as.character(new_data$Term_Description)
    new_data$non_Signif_Snw_Genes <- as.character(new_data$non_Signif_Snw_Genes)
    new_data$Up_regulated <- as.character(new_data$Up_regulated)
    new_data$Down_regulated <- as.character(new_data$Down_regulated)
    new_data$all_pathway_genes <- as.character(new_data$all_pathway_genes)
    new_data$Status <- as.character(new_data$Status)
    new_data$Fold_Enrichment <- as.numeric(new_data$Fold_Enrichment)
    new_data$occurrence <- as.numeric(new_data$occurrence)
    new_data$support <- as.numeric(new_data$support)
    new_data$lowest_p <- as.numeric(new_data$lowest_p)
    new_data$highest_p <- as.numeric(new_data$highest_p)
    new_data$num_genes_in_path <- as.numeric(new_data$num_genes_in_path)
    new_data$Cluster <- as.numeric(new_data$Cluster)
    
    for(i in colnames(new_data)){
      new_na <- sum(is.na(new_data[[i]]))
      if(!new_na == original_na_counts[[i]]){
        showErrorModal(paste0('Coercion of column ', i, ' to proper class introduced na values. ', i, ' requires numeric.'))
        return()
      }
    }
    
    return(new_data)
    
  }
}


# Reads the raw counts uploads
readCountsUpload <- function(fileDataPath){
  
  #                         READ FILE
  #_______________________________________________________________
  
  #Read the file as tsv
  if(grepl('\\.tsv$', fileDataPath, ignore.case = TRUE)){
    counts <- read.csv(fileDataPath, sep = "\t")
  }
  #Read the file as a csv
  else if(grepl('\\.csv$', fileDataPath, ignore.case = TRUE)){
    counts <- read.csv(fileDataPath)
  }
  #Unreadable file type
  else{
    showErrorModal('The counts matrix file is not a readable file')
    return()
  }
  
  #                   CHECK REQUIRED COLUMS
  #_______________________________________________________________
  
  #Check for req columns
  if(! all( c('gene_id', 'gene_name') %in% colnames(counts) ) ){
    showErrorModal('Required columns ( gene_id or gene_name ) are missing in the counts matrix')
    return()
  }
  
  #                   CHECK DUPLICATE IDS
  #_______________________________________________________________
  
  #Check duplicate ids
  if(anyDuplicated(counts$gene_id)){
    showErrorModal('There were duplicate ids in the counts matrix')
    return()
  }
  
  #                   FIND SAMPLE COLUMN NAMES
  #_______________________________________________________________
  sample_cols <- counts %>% select(starts_with('.'))
  if(ncol(sample_cols) < 2){
    showErrorModal('There were not enough samples present')
    return()
  }
  
  #                   SET THE RAW COUNTS DATA FRAME
  #_______________________________________________________________
  raw_counts <- counts
  rownames(raw_counts) <- raw_counts$gene_name
  raw_counts <- raw_counts %>% select(all_of(colnames(sample_cols)))
  
  #                       FILTER THE COUNTS
  #_______________________________________________________________
  filtered_counts <- filterCounts(raw_counts)
  
  #               CREATE A DATAFRAME OF GENE NAMES
  #_______________________________________________________________
  gene_names <- filtered_counts %>% 
    mutate(gene_name = rownames(filtered_counts)) %>%
    select(gene_name)
  
  #RETURN
  return(list(gene_names, raw_counts, filtered_counts))
  
}
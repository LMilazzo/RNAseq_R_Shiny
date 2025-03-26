
runDeseq <- function(metadata, counts, formula){
  
  tryCatch({
    
    #         CREATE THE DESeq DATA SET   
    #_____________________________________________________________----
    ddsc <- DESeqDataSetFromMatrix(countData = round(counts),
                                   colData = metadata,
                                   design = formula)
    print(ddsc)
    
    showModal(
      modalDialog(
       div(
         tags$span("Running DESeq "),
         tags$span(class="custom-loader"),
         style = "margin: 0px; padding: 0px; display: inline-flex; font-size: 16px;"
       ),
      footer = NULL,
      size = "s")
    )
    
    DESeqDataSet <- DESeq(ddsc)
    #----
    
    #                  SELECT CONTRASTS
    #_____________________________________________________________----
    contrast_list <- list()
    #Make options
    for (col_name in intersect(colnames(colData(ddsc)), all.vars(design(ddsc))) ){
      levels <- levels(colData(ddsc)[[col_name]])
      combinations <- combn(levels, 2, simplify = FALSE)
      for (combo in combinations) {
        combo_name1 <- paste(combo, collapse = " vs. ")
        combo_name2 <- paste(rev(combo), collapse = " vs. ")
        contrast_list[[combo_name1]] <- c(col_name, combo)
        contrast_list[[combo_name2]] <- c(col_name, rev(combo))
      }
    }
  
    default <- contrast_list[1]
  
    #----
    
    Object.List <- list(
      
      "contrast_list" = contrast_list,
      "ddsc" = DESeqDataSet,
      "results_ddsc" = as.data.frame(results(DESeqDataSet, contrast = default[[1]])),
      "vst_obj" = vst(DESeqDataSet, blind = TRUE, nsub = 50),
      "vst_counts" = as.data.frame(assay(vst(DESeqDataSet, blind = TRUE, nsub = 50))),
      "normalized_counts" = as.data.frame(counts(DESeqDataSet, normalized = TRUE))
      
    )
  
    showModal(modalDialog("DESeq complete with no fatal errors", easyClose = TRUE, footer = NULL))
    
    return(Object.List)
    
  }, error = function(e) {
    
    showErrorModal(paste("Error in DESeq2 process:", e$message))
       
  })
  
}



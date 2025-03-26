PathUploadsServer <- function(id){
  moduleServer(id, function(input, output, session){
    
    ns <- session$ns
    
    #User Continues the Session
    #Void
    
    
    #User uploads some new data
    observeEvent(input$UploadPrevPathData, {
      showModal(
        modalDialog(
          div( 
            
            p('Pathway Analysis Results'),
            HTML('<p>This is the standard output from a pathfindR clustered experiment.</p>'),
            p("Accepted File Types:  .csv   .tsv "),
            HTML('<p>Must contain all following columns:</p>'),
            HTML('<p>(ID,	Term_Description,	Fold_Enrichment, occurrence, support,	lowest_p,	highest_p,	non_Signif_Snw_Genes,	Up_regulated,	Down_regulated,	all_pathway_genes,	num_genes_in_path,	Cluster, Status)</p>'),
            
            fileInput('Path_Prev_Data', 'Pathway Analysis Results'),
            
            HTML('<p>Normalized Gene Counts'),
            p("Accepted File Types:  .csv   .tsv "),
            HTML('<p>Format: (gene_name, samples...)</p>'),
            p('All sample column headers should start with a "." character to be recognized'),
            
            fileInput('Path_Norm_Counts', 'Gene Counts'),
            
            h5("Sample Conditions Table"),
            p("Accepted File Types:  .csv   .tsv "),
            p("Format: (sample, conditions...) "),
            p("Sample column should include sample names that are found in your experiment data"),
            
            fileInput("Path_Meta_Data", "Sample Conditions")
            
          ), footer = actionButton(ns("LoadPathPrevData"), 'Finish'),
          easyClose = TRUE
        )
      )
      
    })
    
    
  })
}
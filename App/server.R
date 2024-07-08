#_________________________________Shiny Server__________________________________________
server <- function(input, output) {
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~######_______________Observables________________######~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #-----Observe UI Event------# ----> merged_gene_counts_uploaded_file
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  observeEvent(input$merged_gene_counts_uploaded_file,
    {
      #Upon submission of UI object [ merged_gene_counts_uploaded_file ]      
      
      #Goal:
      #             |----->> gene_names
      #   input ----|
      #             |----->> raw_counts

      #Expect as return c(raw_counts, gene_names)           
      func_return <- Set_Clean_Counts(input$merged_gene_counts_uploaded_file$datapath)
      
      raw_counts <<- data.frame(func_return[1])
      
      gene_names <<- data.frame(func_return[2])
      
    }             
  )
  
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#~~~~~######____________Output $ Objects______________######~~~~~#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #--------raw_counts_PreviewTable---------#
  #~~~~~~~~~~~~~~Data Table~~~~~~~~~~~~~~~~#
  #A Data Table preview of what is inside the uploaded counts table 
  #after formating correctly
    #Expected Format:
    #axis.names| Sample1 | Sample 2 | Sample 3
    #Gene_ID_1 |    x    |    x     |    x
    #Gene_ID_2 |    x    |    x     |    x
  output$raw_counts_PreviewTable <- renderDT(
    {
    req(input$merged_gene_counts_uploaded_file)
    req(input$raw_counts_matrix_Preview)
    
    #                 subject
    makePreviewTable(raw_counts,input$raw_counts_matrix_Preview)
    
    }, rownames = TRUE
  )

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #------geneID_geneName_PreviewTable------#
  #~~~~~~~~~~~~~~Data Table~~~~~~~~~~~~~~~~#
  #A Data Table preview of all gene ids and their gene names
    #Expected Format:
    #axis.names|    gene    | gene_name
    #   1      | Gene_ID_1  |   Actin    
    #   2      | Gene_ID_2  |   COOR4_7    
  output$geneID_geneName_PreviewTable <- renderDT(
    {
      req(input$merged_gene_counts_uploaded_file)
      req(input$geneID_geneName_Preview)
      
      #                 subject
      makePreviewTable(gene_names,input$geneID_geneName_Preview)
      
    }, rownames = TRUE
  )
  
  
  
  
  
  
}#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X#X


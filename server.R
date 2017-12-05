## server.R
library(readxl)

function (input, output) {
  ## Load JUMP -q output file (either id_uni_pep_quan.xlsx or id_uni_prot_quan.xlsx)
  inputData = reactive({
    inFile = input$inputFile
    if (length(grep("pep", inputFile))) {
      tbl = read_excel(inputFile, skip = 4) # Peptide publication table
    } else {
      tbl = read_excel(inputFile, skip = 1) # Protein publication table
    }
    df = as.data.frame(tbl)
  })
  
  ## Number of groups to be compared
  nGroups = reactive({input$nGroups})
  output$sampleGroups = renderUI({
    data = inputData()
    sampleNames = grep('sig', colnames(data))
    for (i in 1:nGroups) {
      groupName = paste("Group", i)
      checkboxGroupInput(groupName, choices = sampleNames, selected = 1)
    }
  })
}
# longitudinal-metabolite-analysis-in-mice

contains all files for longitudinal metabolite analysis 




.rmd file: R code for data analysis 
metabolite_functions.R: functions created during the analysis 

# code for fella.RData
graph = buildGraphFromKEGGREST(
 organism = "mmu",
 filter.path = c("01100", "01200", "01210", "01212", "01230"))

buildDataFromGraph(keggdata.graph = graph,
   databaseDir = "./mmu_database",
    internalDir = FALSE,
    matrices = c("hypergeom", "diffusion", "pagerank"), 
    normality = c("diffusion", "pagerank"),
    niter = 100)

fella.data = loadKEGGdata(
  databaseDir = "./mmu_database",
  internalDir = FALSE,
  loadMatrix = NULL)
save(fella.data, file = "./fella.RData")

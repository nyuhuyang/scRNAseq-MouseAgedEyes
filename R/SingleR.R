library(SingleR)
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
(load(file = "data/object_mm10_young_aged_eyes_2_20190712.Rda"))
attach(mouse.rnaseq)
singler = CreateSinglerObject(object_data, annot = NULL, project.name="Aged_eyes",
                              min.genes = 500,technology = "10X", species = "Mouse", citation = "",
                              ref.list = list(mouse.rnaseq),
                              normalize.gene.length = F, variable.genes = "de",
                              fine.tune = T, do.signatures = F, clusters = NULL)
save(singler,file="output/singler_T_mm10_young_aged_eyes_2_20190712.Rda")
  

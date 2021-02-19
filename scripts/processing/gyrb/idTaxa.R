

assign_taxonomy_w_idTAXA <- function(asvs_faa,db_faa){
  #read in ref db faa
  db <- readAAStringSet(db_faa)
  
  # parse the headers to obtain a taxonomy
  s <- strsplit(names(db), ";")
  domain <- sapply(s, `[`, 1)
  phylum <- sapply(s, `[`, 2)
  class <- sapply(s, `[`, 3)
  order <- sapply(s, `[`, 4)
  family <- sapply(s, `[`, 5)
  genus <- sapply(s, `[`, 6)
  taxonomy <- paste("Root", phylum, class, order, family, genus, sep=";")
  taxonomy
  # train the classifier
  trainingSet <- LearnTaxa(db, taxonomy)
  
  # view information about the classifier
  plot(trainingSet)
  
  #read in asvs
  asvs <- readAAStringSet(asvs_faa)
  
  # classify the test sequences
  ids <- IdTaxa(asvs, trainingSet, strand="top",threshold = 30)
  plot <- plot(ids, trainingSet)
  assignment <- sapply(ids,
                       function(x)
                         paste(x$taxon,
                               collapse=";"))
  tax_table <- data.frame(assignment)
  tax_table_sep <-separate(tax_table,col = assignment, sep=';',into = c('Root','Phylum','Class','Order','Family','Genus'))
  tax_table_filled <- data.frame(t(apply(tax_table_sep, 1, zoo::na.locf)))
  return(tax_table_filled)
}

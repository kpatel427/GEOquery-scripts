# download GEO data
.libPaths("~/KP/RShiny/packages/")

library(GEOquery)
library(xml2)
library(R.utils)



#open file with gse_accessions.txt(a text file containing all GSExxxxx accessions)
con = file("gse_accessions_filtered.txt", "r")

while ( TRUE ) {
  line = readLines(con, n=1)
  if ( length(line) == 0 ) {
    break
  }
  print(line)
  geo_id <- toupper(line) # converting lower case letters to upper case
  if(file.exists(paste0("~/KP/RShiny/GEOdata/",geo_id))){
    print("Directory exists!") # repetitive GSE ID; skip to next one
  }
  else {
  dir.create(file.path("~/KP/RShiny/GEOdata/", geo_id)) # create directory
  gse <- getGEO(geo_id, GSEMatrix=TRUE, destdir = paste0("~/KP/RShiny/GEOdata/",geo_id) )
  print(paste0("Data download finished for ",geo_id))
  }

}

close(con)

.libPaths("~/KP/RShiny/packages/")

library(xml2)
library(GEOquery)
library(R.utils)

# open file with gse_accessions.txt
con = file("gse_accessions_filtered.txt", "r")

while ( TRUE ) {
  line = readLines(con, n=1)
  if ( length(line) == 0 ) {
    break
  }
  
  else{
    print(line)
    geo_id <- toupper(line)
    gse <- getGEO(geo_id,GSEMatrix=TRUE)
    
    if(length(gse) == 1){ #** -->
    
        data <- as.data.frame(exprs(gse[[1]])) # getting expression data
        data <- setDT(data, keep.rownames = "ID" )
        featureData <- as.data.frame(gse[[1]]@featureData@data) # fetching features to get ID and gene symbols
        featureData <- featureData[,c("ID", "Gene Symbol")]
        data <- merge(data,featureData, by = "ID")
        
        name = attr(gse, "names")[[1]]
        save(data, file = paste0("~/KP/RShiny/GEOdata/",geo_id,"/", name,"-expr-data.rda"))
        
        # phenotypic data: column characteristics_ch1.1 contains age; row names are sample names; 
        head(pData(gse[[1]])[, 1:25])
        # title; type; source_name; characteristics_ch1.1 (age)
        phenoData <-pData(phenoData(gse[[1]]))[1:5,c(1,6,8,11)]
        save(phenoData, file = paste0("~/KP/RShiny/GEOdata/",geo_id,"/",name,"-pheno-data.rda"))
        
        
    } #** <--
    
    # if gse has more than one file
    if(length(gse) > 1){ #* --> 
      i = length(gse) # getting total number of files
      for(x in seq(1,i)){  # iterating for every dataset within this gse
        #print(x)
        data <- as.data.frame(exprs(gse[[x]])) # getting expression data
        data <- setDT(data, keep.rownames = "ID" )
        featureData <- as.data.frame(gse[[x]]@featureData@data) # fetching features to get ID and gene symbols
        featureData <- featureData[,c("ID", "Gene Symbol")]
        data <- merge(data,featureData, by = "ID")
        
        name = attr(gse, "names")[[x]]
        save(data, file = paste0("~/KP/RShiny/GEOdata/",geo_id,"/",name,"-expr-data.rda"))
        
        # phenotypic data: column characteristics_ch1.1 contains age; row names are sample names; 
        head(pData(gse[[x]])[, 1:25])
        # title; type; source_name; characteristics_ch1.1 (age)
        phenoData <-pData(phenoData(gse[[x]]))[1:5,c(1,6,8,11)]
        save(phenoData, file = paste0("~/KP/RShiny/GEOdata/",geo_id,"/",name,"-pheno-data.rda"))
        }
      
    } #* <--

  } # closing of else bracket
} # closing while bracket

#gse <- getGEO("GSE61578",GSEMatrix=TRUE)

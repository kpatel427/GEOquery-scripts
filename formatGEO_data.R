.libPaths("~/KP/RShiny/packages/")

#install.packages("data.table", repos='http://cran.us.r-project.org')

library(xml2)
library(GEOquery)
library(R.utils)
library(data.table)
library(dplyr)
library(tidyverse)
library(stringr)


# open file with gse_accessions.txt
con = file("GSE_accessions.txt", "r")

while ( TRUE ) {
  line = readLines(con, n=1)
  if ( length(line) == 0 ) {
    break
  }

  else{
    print(line)
    geo_id <- toupper(line)
    #geo_id <- "GSE48433"
    gse <- getGEO(geo_id,GSEMatrix=TRUE)
    
    if(length(gse) == 1){ #** -->
    
        data <- as.data.frame(exprs(gse[[1]])) # getting expression data
        data <- setDT(data, keep.rownames = "ID" )
        
        # if featureData has no column called Gene Symbol:
        
        
        featureData <- as.data.frame(gse[[1]]@featureData@data) # fetching features to get ID and gene symbols
        featureData <- setDT(featureData)[,c("ID", "Gene Symbol")]
        data <- merge(data,featureData, by = "ID")
        
        name = attr(gse, "names")[[1]]
        print("Saving expr data...")
        save(data, file = paste0("~KP/RShiny/GEOdata/in-process/",geo_id,"/", name,"-expr-data.Rdata"))
        
        # phenotypic data: column characteristics_ch1.1 contains age; row names are sample names; 
        head(pData(gse[[1]])[, 1:25])
        # title; type; source_name; characteristics_ch1.1 (age)
        phenoData <-pData(phenoData(gse[[1]]))[1:5,c(1,6,8,11)]
        print("Saving pheno data...")
        save(phenoData, file = paste0("~/KP/RShiny/GEOdata/in-process/",geo_id,"/",name,"-pheno-data.RData"))
        
        
    } #** <--
    
    # if gse has more than one file
    if(length(gse) > 1){ #* --> 
      i = length(gse) # getting total number of files
      for(x in seq(1,i)){  # iterating for every dataset within this gse
        #print(x)
        data <- as.data.frame(exprs(gse[[x]])) # getting expression data
        data <- setDT(data, keep.rownames = "ID" )
        
        # if featureData has no column called Gene Symbol:
        
        
        featureData <- as.data.frame(gse[[x]]@featureData@data) # fetching features to get ID and gene symbols
        featureData <- setDT(featureData)[,c("ID", "Gene Symbol")]
        data <- merge(data,featureData, by = "ID")
        
        name = attr(gse, "names")[[x]]
        print("Saving expr data...")
        save(data, file = paste0("~/KP/RShiny/GEOdata/in-process/",geo_id,"/",name,"-expr-data.RData"))
        
        # phenotypic data: column characteristics_ch1.1 contains age; row names are sample names; 
        head(pData(gse[[x]])[, 1:25])
        # title; type; source_name; characteristics_ch1.1 (age)
        phenoData <-pData(phenoData(gse[[x]]))[1:5,c(1,6,8,11)]
        print("Saving pheno data...")
        save(phenoData, file = paste0("~/KP/RShiny/GEOdata/in-process/",geo_id,"/",name,"-pheno-data.RData"))
        }
      
    } #* <--

  } # closing of else bracket
} # closing while bracket

geo_id <- "GSE66586"
gse <- getGEO(geo_id,GSEMatrix=TRUE)

length(gse)

data <- as.data.frame(exprs(gse[[1]])) # getting expression data
data <- setDT(data, keep.rownames = "ID" )

# if featureData has no column called Gene Symbol:
featureData <- as.data.frame(gse[[1]]@featureData@data) # fetching features to get ID and gene symbols
#featureData <-  setDT(featureData)[,c("V1","Gene Symbol","V2"):= sapply(str_split(featureData$gene_assignment, " // ",  n = 3), `[`, 2)]
featureData <- setDT(featureData)[,c("ID", "Gene Symbol")]
#featureData$ID <-  as.character(featureData$ID)

# ### when feature data us not available
# library(annotate)
# library(hgu133a.db)
# select(hgu133a.db, c(data$ID), c("SYMBOL","ENTREZID", "GENENAME"))
# ###

# merging expression data with feature data
data <- merge(data,featureData, by = "ID")
GSE665863_RMA_data <- data 

name = attr(gse, "names")[[1]]
print("Saving expr data...")
save(GSE66586_RMA_data, file = "GSE66586_RMA_data.RData")

# phenotypic data: column characteristics_ch1.1 contains age; row names are sample names; 
head(pData(gse[[1]])[, 1:25])
# title; type; source_name; characteristics_ch1.1 (age)
phenoData <-pData(phenoData(gse[[1]]))
GSE66586_RMA_mData <-  phenoData
save(GSE66586_RMA_mData, file = "GSE665863_RMA_mData.Rdata")




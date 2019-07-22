# install the core bioconductor packages, if not already installed
source("http://bioconductor.org/biocLite.R")
biocLite()

# install additional bioconductor libraries, if not already installed
biocLite("GEOquery")  # Get data from NCBI Gene Expression Omnibus (GEO)
biocLite("affy")  # Methods for Affymetrix Oligonucleotide Arrays
BiocManager::install("oligo")
biocLite("hgu133a.db", type = "source")  # GSE1297: Platform_title = [HG-U133A]
biocLite("hgu133acdf")

library(GEOquery)
library(affy)
library(hgu133a.db)
library(hgu133acdf)
library(data.table)
library(dplyr)
library(tidyverse)
library(oligo)
library(stringr)

# get supplimentary files
getGEOSuppFiles("GSE16476")
#---------------------------------- untar and RMA normalization -------------------------------------#

untar("GSE17714_RAW.tar", exdir = "data")
cels = list.files("data/", pattern = "CEL|cel")
# sometiles, it is 'CEL', you need to check it first
sapply(paste("data", cels, sep = "/"), gunzip)

cels = list.files("data/", pattern = "CEL|cel")
# sometiles, it is 'CEL', you need to check it first

# Set working directory for normalization
setwd("./GEOdata/GSE17714/GSE17714/data")

### for Affy package
raw.data = ReadAffy(verbose = FALSE, filenames = cels, cdfname = "hgu133acdf")
# perform RMA normalization (log2)
data.rma.norm = affy::rma(raw.data)
# Get the expression estimates for each array
rma = exprs(data.rma.norm)

### for oligo package
rawData <- read.celfiles(cels)
# perform RMA normalization (log2)
data.rma.norm = rma(rawData)
# Get the expression estimates for each array
rma = exprs(data.rma.norm)


#---------------------------------- Annotation -------------------------------------#
tt = cbind(row.names(rma), rma)
colnames(tt) = c("ProbID", sub(".cel", "", colnames(rma), ignore.case = TRUE))

#### to merge with gene symbols
# if featureData has no column called Gene Symbol:
geo_id <- "GSE17714"
gse <- getGEO(geo_id,GSEMatrix=TRUE)
featureData <- as.data.frame(gse[[1]]@featureData@data) # fetching features to get ID and gene symbols
#featureData <-  setDT(featureData)[,c("V1","Gene Symbol","V2"):= sapply(str_split(featureData$gene_assignment, " // ",  n = 3), `[`, 2)]
featureData <- setDT(featureData)[,c("ID", "Gene Symbol")]
#featureData$ID <-  as.character(featureData$ID)
GSE17714_RMA_data <- merge(tt,featureData, by.x = "ProbID", by.y = "ID")

#colnames(GSE107333_series_matrix_expr_data) <- sub("_exp.*", "", colnames(GSE107333_series_matrix_expr_data))


#---------------------------------- Expression matrix Manipulation -------------------------------------#
# converting gene symbols to row names and keeping rows with max mean in each row
# removing ID column
GSE17714_RMA_data <- select(GSE17714_RMA_data,-c("ProbID"))


copy_GSE17714_RMA_data <-  GSE17714_RMA_data
#GSE16476_series_matrix_expr_data <-  copy_GSE16476_series_matrix_expr_data

cols <- names(GSE17714_RMA_data)[1:22]
setDT(GSE17714_RMA_data)[, (cols) := lapply(.SD, as.character), .SDcols = cols]
setDT(GSE17714_RMA_data)[, (cols) := lapply(.SD, as.numeric), .SDcols = cols]
 
# Since there are duplicated gene symbols -
# calculating means for every row grouping by the gene; and selecting rows with max mean
GSE17714_RMA_data <- setDT(GSE17714_RMA_data)[, .SD[which.max(rowMeans(.SD))], by=`Gene Symbol`]

# Removing if there are any NA's in Gene Symbols
GSE17714_RMA_data <- GSE17714_RMA_data[!(is.na(GSE17714_RMA_data$`Gene Symbol`)),]

# converting column to rownames
GSE17714_RMA_data <- GSE17714_RMA_data %>% column_to_rownames(var="Gene Symbol")

# Saving data
save(GSE17714_RMA_data,file = "GSE17714_RMA_data.RData")


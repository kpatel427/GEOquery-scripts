# install the core bioconductor packages, if not already installed
source("http://bioconductor.org/biocLite.R")
biocLite()

# install additional bioconductor libraries, if not already installed
biocLite("GEOquery")  # Get data from NCBI Gene Expression Omnibus (GEO)
biocLite("affy")  # Methods for Affymetrix Oligonucleotide Arrays
biocLite("hgu133a.db", type = "source")  # GSE1297: Platform_title = [HG-U133A]
biocLite("hgu133acdf")

library(GEOquery)
library(affy)
library(hgu133a.db)
library(hgu133acdf)
library(data.table)
library(dplyr)
library(tidyverse)

#---------------------------------- untar and RMA normalization -------------------------------------#

untar("GSE9807_RAW.tar", exdir = "data")
cels = list.files("data/", pattern = "CEL")
# sometiles, it is 'CEL', you need to check it first
sapply(paste("data", cels, sep = "/"), gunzip)

cels = list.files("data/", pattern = "CEL")
# sometiles, it is 'CEL', you need to check it first

# Set working directory for normalization
setwd("/Volumes/target_nbl_ngs/KP/RShiny/GEOdata/to-normalize/data")
raw.data = ReadAffy(verbose = FALSE, filenames = cels, cdfname = "hgu133acdf")

# perform RMA normalization (log2)
data.rma.norm = rma(raw.data)

# Get the expression estimates for each array
rma = exprs(data.rma.norm)


#---------------------------------- Annotation -------------------------------------#
tt = cbind(row.names(rma), rma)
colnames(tt) = c("ProbID", sub(".cel", "", colnames(rma), ignore.case = TRUE))

#### to merge with gene symbols
geo_id <- "GSE9807"
gse <- getGEO(geo_id,GSEMatrix=TRUE)

length(gse)

# if featureData has no column called Gene Symbol:
featureData <- as.data.frame(gse[[1]]@featureData@data) # fetching features to get ID and gene symbols
#featureData <-  setDT(featureData)[,c("V1","Gene Symbol","V2"):= tstrsplit(gene_assignment, " // ", fixed = TRUE)]

featureData <- setDT(featureData)[,c("ID", "Gene Symbol")]
#featureData$ID <-  as.character(featureData$ID)
GSE9807_series_matrix_expr_data <- merge(tt,featureData, by.x = "ProbID", by.y = "ID")


#---------------------------------- Expression matrix Manipulation -------------------------------------#
# removing ID column
GSE9807_series_matrix_expr_data <- select(GSE9807_series_matrix_expr_data,-c("ProbID"))
GSE9807_series_matrix_expr_data$ProbID <- NULL

copy_GSE9807_series_matrix_expr_data <-  GSE9807_series_matrix_expr_data
GSE9807_series_matrix_expr_data <-  copy_GSE9807_series_matrix_expr_data

cols <- names(GSE9807_series_matrix_expr_data)[1:6]
setDT(GSE9807_series_matrix_expr_data)[, (cols) := lapply(.SD, as.character), .SDcols = cols]
setDT(GSE9807_series_matrix_expr_data)[, (cols) := lapply(.SD, as.numeric), .SDcols = cols]
 
# Since there are duplicated gene symbols -
# calculating means for every row grouping by the gene; and selecting rows with max mean
GSE9807_series_matrix_expr_data <- setDT(GSE9807_series_matrix_expr_data)[, .SD[which.max(rowMeans(.SD))], by=`Gene Symbol`]

# Removing if there are any NA's in Gene Symbols
GSE9807_series_matrix_expr_data <- GSE9807_series_matrix_expr_data[!(is.na(GSE9807_series_matrix_expr_data$`Gene Symbol`)),]

# converting column to rownames
GSE9807_series_matrix_expr_data <- GSE9807_series_matrix_expr_data %>% column_to_rownames(var="Gene Symbol")

# Saving data
save(GSE9807_series_matrix_expr_data,file = "GSE9807_expr_data.RData")


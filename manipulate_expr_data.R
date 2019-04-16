library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)

load("GSE65303_series_matrix.txt.gz_expr_data.RData")

# removing ID column 
data <- select(data,-c("ID"))

# Since there are duplicated gene symbols -
# calculating means for every row grouping by the gene; and selecting rows with max mean
data <- setDT(data)[, .SD[which.max(rowMeans(.SD))], by=`Gene Symbol`]

# Removing if there are any NA's in Gene Symbols
data <- data[!(is.na(data$`Gene Symbol`)),]

# converting column to rownames
data <- data %>% column_to_rownames(var="Gene Symbol")

GSE65303_series_matrix_expr_data <-  data

# Saving data
save(GSE65303_series_matrix_expr_data,file = "GSE65303_expr_data.RData")



# file name = GSE9807_expr_data.RData
# R object name = GSE9807_series_matrix_expr_data


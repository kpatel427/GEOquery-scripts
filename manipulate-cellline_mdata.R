library(tidyverse)
library(data.table)
load("cellline_mData.RData")
load("GSE66586_RMA_data.Rdata")
load("GSE665863_RMA_mData.Rdata")


setDT(GSE66586_RMA_mData, keep.rownames = TRUE)[]
celline.data <-  subset(GSE66586_RMA_mData, select = c("rn","title","source_name_ch1"))

celline.data$source_name_ch1 <- gsub("Control","unknown", celline.data$source_name_ch1)
celline.data$source_name_ch1 <- gsub("MYCN_Amplified","Amplified", celline.data$source_name_ch1)
celline.data$source_name_ch1 <- gsub("MYCN_NotAmplified","Non-amplified", celline.data$source_name_ch1)
celline.data$source_name_ch1 <- gsub("_5","", celline.data$source_name_ch1)

titles <- c(celline.data$title)

# replace GSM's with cell line in expression data matrix
for( i in 1:length(colnames(GSE66586_RMA_data))){
  print(i)
  print(names(GSE66586_RMA_data)[i])
  
  if(names(GSE66586_RMA_data)[i] %in% celline.data$rn){
    names(GSE66586_RMA_data)[i] = titles[i]
  }
}

# Saving expr data file with changed column names
save(GSE66586_RMA_data, file = "GSE66586_RMA_data.RData")

# populating columns to merge with cellLine_mData

celline.data$ALK_Status <- "Unknown"
celline.data$TP53_Status <- "Unknown"

setnames(celline.data, "title","CellLine")
setnames(celline.data, "source_name_ch1","MYCN_Status")
celline.data <- subset(celline.data, select = c("CellLine","MYCN_Status","ALK_Status","TP53_Status"))
celline.data$CellLine.2 <- celline.data$CellLine


celline.data<- celline.data %>%
  column_to_rownames(var = "CellLine.2")

cellline_mData <-  rbind(cellline_mData, celline.data)

save(cellline_mData,file = "cellline_mData.RData")



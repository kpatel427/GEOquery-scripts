# To format GEO GSE expression and metadata files for those ids which cannot be fetched using GEOquery package
library(data.table)

# --------- for metadata
family.soft <- read.delim("GSE87042_family.soft", header = T)

# filtering required lines
metadata <- as.data.frame(family.soft[grep('!Sample_title|!Sample_geo_accession|!Sample_character|^SAMPLE',family.soft$X.DATABASE...GeoMiame),])
colnames(metadata)[1] <-  "X"
setDT(metadata)[,c('x1','x2') := tstrsplit(X, ' = ', fixed = T)]

tissue <- as.data.frame(metadata[grep('tissue',x2),x2])
setDT(tissue)[,c('y','tissue'):=tstrsplit(`metadata[grep(\"tissue\", x2), x2]`,':', fixed = T)]
tissue <-  subset(tissue[,3])

age <- as.data.frame(metadata[grep('age',x2),x2])
setDT(age)[,c('y','age'):=tstrsplit(`metadata[grep(\"age\", x2), x2]`,':', fixed = T)]
age <-  subset(age[,3])

tumor <- as.data.frame(metadata[grep('tumor',x2),x2])
setDT(tumor)[,c('y','tumor'):=tstrsplit(`metadata[grep(\"tumor\", x2), x2]`,':', fixed = T)]
tumor <-  subset(tumor[,3])

diseaseState <- as.data.frame(metadata[grep('disease state',x2),x2])
setDT(diseaseState)[,c('y','disease state'):= tstrsplit(`metadata[grep(\"disease state\", x2), x2]`,':', fixed=T)]
diseaseState <-  subset(diseaseState[,3])

genotype <- as.data.frame(metadata[grep('genotype',x2),x2])
setDT(genotype)[,c('y','genotype'):= tstrsplit(`metadata[grep(\"genotype\", x2), x2]`,':', fixed=T)]
genotype <-  subset(genotype[,3])

cellLine <- as.data.frame(metadata[grep('cell line',x2),x2])
setDT(cellLine)[,c('y','cell line'):= tstrsplit(`metadata[grep(\"cell line\", x2), x2]`,':', fixed=T)]
cellLine <-  subset(cellLine[,3])

dlbcl <- as.data.frame(metadata[grep('dlbcl',x2),x2])
setDT(dlbcl)[,c('y','dlbcl_subtype'):= tstrsplit(`metadata[grep(\"dlbcl\", x2), x2]`,':', fixed=T)]
dlbcl <-  subset(dlbcl[,3])

sample.title <- as.data.frame(metadata[grep('!Sample_title',x1),x2])
colnames(sample.title)[1] <- 'sample_title'

geo.id <- as.data.frame(metadata[grep('!Sample_geo_accession',x1),x2])

# cbind columns
meta.merge <-  cbind(geo.id, sample.title)
rownames(meta.merge) <- meta.merge[,1]
meta.merge[,1] <-  NULL


# -------- expression data
expr <- read.delim("GSM2682873_CF_00004_RNA.genes.results.txt", header = T)

# when expression data for each sample is in a seperate file
files <- list.files(path = '.', pattern = "*.txt")

# saving every sample's expression values
for(x in files){
  print(x)
  name <- str_split(x,'_')[[1]][1] # saving dataframe by GSM
  
  #assign(name, read.delim(paste0('~/KP/RShiny/new_portal_datasets/GSE100427_SRP110289/',x), header = T))
  data <- read.delim(paste0('~/KP/RShiny/new_portal_datasets/GSE100427_SRP110289/',x), header = T)
  colnames(data)[6] <- paste0(name,'_TPM')
  colnames(data)[7] <- paste0(name,'_FPKM')
  assign(name, data[,c(1,6,7)])
  
}

# merge all samples based on gene_id
expr.merge <-  Reduce(function(x, y) merge(x, y, by = 'gene_id'), 
       list(GSM2682871, GSM2682872, GSM2682873, GSM2682874, GSM2682875,
            GSM2682876, GSM2682877, GSM2682878, GSM2682879, GSM2682880))


expr.FPKM <- expr.merge[,c(1,3,5,7,9,11,13,15,17,19,21)]
expr.TPM <- expr.merge[,c(1,2,4,6,8,10,12,14,16,18,20)]

colnames(expr.FPKM) <- gsub('_FPKM','', colnames(expr.FPKM))
expr.FPKM <- expr.FPKM %>%
  column_to_rownames(var = 'gene_id')


colnames(expr.TPM) <- gsub('_TPM','', colnames(expr.TPM))
expr.TPM <- expr.TPM %>%
  column_to_rownames(var = 'gene_id')

# when expression data is in the form of matrix already
# subset expr
expr <- expr[,c(1,11:16)]

# converting column to rownames
expr <-  expr %>%
  column_to_rownames(var = "gene")

# replacing strings in column names in expr
colnames(expr) <-  gsub('_RPKM','', colnames(expr))

# replacing column names with geo ids
for(i in 1:length(colnames(expr))){
  #print(colnames(expr)[i])
  #print(paste0(colnames(expr)[i],' ', meta.merge$sample_title[i],' ', rownames(meta.merge)[i]))
  colnames(expr)[i] = rownames(meta.merge)[i]
  
}



# saving meta and expr files
GSE89775_FPKM_mData <-  meta.merge
#GSE100427_FPKM_data <- expr.FPKM
save(GSE87042_FPKM_mData, file = 'GSE87042_FPKM_mData.RData')
#save(GSE100427_FPKM_data, file = 'GSE100427_FPKM_data.RData')


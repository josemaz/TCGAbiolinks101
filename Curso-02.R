library(SummarizedExperiment)
library(TCGAbiolinks)
require(EDASeq)
source("tools.R")

qry.rna <- GDCquery(project = "NCICCR-DLBCL",
                    data.category= "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "HTSeq - Counts")
GDCdownload(qry.rna)

dat <- qry.rna[[1]][[1]]
table(as.factor(dat$sample_type))

rnas.raw <- GDCprepare(qry.rna, summarizedExperiment = TRUE)
saveRDS(rnas.raw, file = "pipeline/DLBCL-rnas-raw.rds")
rnas.raw <- readRDS("pipeline/DLBCL-rnas-raw.rds")

st1.samples <- rnas.raw[,rnas.raw$ann_arbor_clinical_stage == "Stage I"]$sample
st2.samples <- rnas.raw[,rnas.raw$ann_arbor_clinical_stage == "Stage II"]$sample

# ANOTACION
data <- assay(rnas.raw)
rownames(data) <- rowData(rnas.raw)$external_gene_name
rownames(rnas.raw) <- rowData(rnas.raw)$external_gene_name
head(rownames(data))

#! Filtro de los genes
dim(data)
dataFilt <- TCGAanalyze_Filtering(tabDF = data,
                                  method = "quantile",
                                  qnt.cut = 0.25)
threshold <- round(dim(data)[2]/2)
ridx <- rowSums(dataFilt == 0) <= threshold
dataFilt <- dataFilt[ridx, ]
print(dim(dataFilt))
ridx <- rowMeans(dataFilt) >= 10
dataFilt <- dataFilt[ridx, ]
print(dim(dataFilt))
rnas.raw <- rnas.raw[rownames(rnas.raw) %in% rownames(dataFilt), ]

# #! Filtro de las muestras
# rnas.raw$ann_arbor_clinical_stage <-   as.factor(rnas.raw$ann_arbor_clinical_stage)
# table(rnas.raw$ann_arbor_clinical_stage)
# etapas <- c("Stage I","Stage II","Stage III","Stage IV")
# etapas <- c("Stage I","Stage II")
# rnas.raw <- rnas.raw[,rnas.raw$ann_arbor_clinical_stage %in% etapas]
# 
# dataNorm <- TCGAanalyze_Normalization(tabDF = assay(rnas.raw), geneInfo = geneInfo)
# dim(dataNorm)
# rnas.raw <- rnas.raw[rownames(rnas.raw) %in% rownames(dataNorm), ]

rnas.raw <- rnas.raw[!duplicated(rownames(rnas.raw)),]

print(dim(rnas.raw))

mypca(rnas.raw,"ann_arbor_clinical_stage",fout="PCAbefore.png")

rnas.raw <- rnas.raw[!duplicated(rownames(rnas.raw)),]
ginfo <- geneInfoHT[rownames(geneInfoHT) %in% rowData(rnas.raw)$ensembl_gene_id,]
rnas.raw <- rnas.raw[rowData(rnas.raw)$ensembl_gene_id %in% rownames(geneInfoHT),]

ln.data <- withinLaneNormalization(assay(rnas.raw),
                                   ginfo$geneLength, which = "full")
gcn.data <- withinLaneNormalization(ln.data , ginfo$gcContent,
                                    which = "full")
norm.counts <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
noiseqData <- NOISeq::readData( norm.counts , factors = as.data.frame(rnas.raw$ann_arbor_clinical_stage))
mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n",  logtransf = FALSE)
assay(rnas.raw) <- exprs(mydata2corr1)

mypca(rnas.raw,"ann_arbor_clinical_stage",fout="PCAafter.png")

dim(assay(rnas.raw[,st1.samples]))

# debugonce(TCGAanalyze_DEA)
dataDEGs <- TCGAanalyze_DEA(mat1 = assay(rnas.raw[,c(1,2)]),
                            mat2 = assay(rnas.raw[,c(2,4)]),
                            Cond1type = "Stage I",
                            Cond2type = "Stage II",
                            metadata = FALSE,
                            method = "glmLRT")
colnames(rnas.raw[,st1.samples])

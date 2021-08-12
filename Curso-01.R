library(TCGAbiolinks)

qry.rna <- GDCquery(project = "TCGA-BRCA",
                   data.category= "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification",
                   workflow.type = "HTSeq - Counts")

GDCdownload(qry.rna)

dat <- qry.rna[[1]][[1]]
table(as.factor(dat$sample_type))

rnas.raw <- GDCprepare(qry.rna, summarizedExperiment = TRUE)
saveRDS(rnas.raw, file = "pipeline/brca-rnas-raw.rds")
rnas.raw <- readRDS("pipeline/brca-rnas-raw.rds")

# Which samples are Primary Tumor
dataSmTP <- TCGAquery_SampleTypes(getResults(qry.rna,cols="cases"),"TP")
# which samples are solid tissue normal
dataSmNT <- TCGAquery_SampleTypes(getResults(qry.rna,cols="cases"),"NT")

rnas.raw <- rnas.raw[rownames(rnas.raw) %in% rownames(geneInfoHT),]
ginfo <- geneInfoHT[rownames(geneInfoHT) %in% rownames(rnas.raw),]


dataPrep <- TCGAanalyze_Preprocessing(rnas.raw, cor.cut = 0.6)
dataNorm <- TCGAanalyze_Normalization(tabDF = dataPrep,
                geneInfo = ginfo,
                method = "gcContent")
dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm,
                                  method = "quantile",
                                  qnt.cut =  0.25)
dim(dataFilt)


dataDEGs <- TCGAanalyze_DEA(mat1 = dataFilt[,dataSmNT],
                            mat2 = dataFilt[,dataSmTP],
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")

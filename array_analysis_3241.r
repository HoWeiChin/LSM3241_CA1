library(GEOquery)
library(affy)
library(limma)
library(hgu133plus2.db)

#get GSE50697 from LSM3241_Assignment1 dir 
gse_file ='GSE50697_family.soft'
gse <- getGEO(filename = gse_file)

#get gse methods
gse_methods <- methods(class=class(gse))
#explore gse meta data
gse_meta<- names(Meta(gse))

#get names of all gsm (the samples)
all_gsm <- GSMList(gse)
gsm_methods <- methods(class=class(all_gsm[[1]]))

#get metadata for a control sample 
gsm_meta_c_s1 <- names(Meta(all_gsm[[1]]))
unique_gsm_meta_c_s1 <- Meta(all_gsm[[1]])[!(gsm_meta_c_s1 %in% gse_meta)]

#get metadata for a treatment sample
gsm_meta_t_s4 <- names(Meta(all_gsm[[4]]))
unique_gsm_meta_t_s4 <- Meta(all_gsm[[4]])[!(gsm_meta_t_s4 %in% gse_meta)]

#create function to extract treatment/control values from gsm samples
get_treatment_value <- function(gsm){
  Meta(gsm)[['characteristics_ch1']][2]
}

#create a treatment/control pd
treatment_col = as.factor(sapply(all_gsm,get_treatment_value))
t_pd <- data.frame(treatment_type=treatment_col)
levels(t_pd$treatment_type) <- c('control', 'pBabe puro miR-203')

gse_cel_files <- list.celfiles(getwd())
affydata <- read.affybatch(gse_cel_files, phenoData=new('AnnotatedDataFrame', t_pd))
phenoData(affydata)

#plot density before rma
pre_rma_data <- exprs(affydata)
plotDensity(pre_rma_data, xlab='original intensities before log transformed', main='Intensity Values of samples before RMA', lwd=2)

#plot density after rma
rma_data <- rma(affydata)
plotDensity(exprs(rma_data), xlab='log transformed intensities', main="Intensity values of samples after RMA", lwd=2)
methods(class=class(rma_data))

#check if what we have done so far is correct
pData(rma_data)
head(exprs(rma_data))

#construct linear model

linear_m <- model.matrix(~0 + rma_data$treatment_type)

colnames(linear_m) <- c('control', 'miRNA')
linear_m #visualise


build_contrasts <- makeContrasts( control - miRNA, levels=linear_m)
build_contrasts #visualise

#linear model fit
lm.fit <- lmFit(rma_data, linear_m)
#fit contrast
lm.contrast <- contrasts.fit(lm.fit, build_contrasts)

#adjustment using Empirical Bayes correction
ebayes <- eBayes(lm.contrast)

#for systems pathway analysis
#get all probs with matching parameters given in topTable
large_prob_set <- topTable(ebayes, number=Inf, p.value=0.05, lfc=0.5)
#goi == genes of interest
large_goi_set <-AnnotationDbi::select(hgu133plus2.db, rownames(large_prob_set),c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
write.csv(large_goi_set, 'large_set_of_genes.csv')


volcanoplot(ebayes)


#pinpoint priority genes in volcano plot
#priority == a small set of genes prioritise for qPCR
genes_of_interest <- na.omit(topTable(ebayes, number=Inf, p.value = 0.05, lfc=2))
volcanoplot(ebayes, main=sprintf("%d interesting features", nrow(genes_of_interest)))
points(genes_of_interest[['logFC']], -log10(genes_of_interest[['P.Value']]), col='red')

#find the names of prioritised genes
name_prioritised_genes <-AnnotationDbi::select(hgu133plus2.db, row.names(genes_of_interest) ,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
name_prioritised_genes[,'LogFC'] <- genes_of_interest[,'logFC']

#sort columns
final_data <- name_prioritised_genes[order(name_prioritised_genes),]
write.csv(na.omit(final_data), "genes_prioritised_for_qPCR.csv")





### Roast GSEA Execution script #####

## Packages ####

library(tidyverse)
library(limma)
library(sva)
library(edgeR)
library(DESeq)
library(statmod)
library(writexl)
library(openxlsx)
library(GSEABase)
library(gplots)
library(org.Hs.eg.db)
library(qvalue)
library(MSstats)
library(here)
library(qdapTools)
library(ggridges)

## Dependencies ####

source(here::here("helperROAST_twogroups.r")) # load helper functions to run the roast analysis

## Prepare execution ####

### Set up backgroud database for GSEA ####

sigdbfile <- 'msigdb.v7.0.symbols.gmt' # last version of the msigdb (available at: http://software.broadinstitute.org/gsea/downloads.jsp)
                                       # This file should be in the same RStudio Project folder from which this script would be executed.

### Set locations of the dataset to analyze ####

inputlocation <- here::here("input_limma_example.txt") # the name of the input limma file (it should have 1 row per protein and 1 column per condition)

### Define experimental design ####

condition1 <- 5 # number of samples associated to the first condition (treatment, stage, patient, etc...)
condition2 <- 7 # number of samples associated to the second condition 

### Select the Geneset database to run the ROAST enrichment against ####

# Currently available: "GO" == 1, "KEGG" == 2, "REACTOME" == 3, "HALLMARK" == 4, "NABA" == 5

db <- 4 # we set it to 5 to run it agains the NABA geneset, for instance.

names(db) <- case_when(db == 1 ~ "GO",
                       db == 2 ~ "KEGG",
                       db == 3 ~ "REACTOME",
                       db == 4 ~ "HALLMARK",
                       db == 5 ~ "NABA") # define the name of the geneset on the previous object 


############## EXECUTION OF THE SCRIPT ################################################

### Load dataset (input) into R ####

mat <- read.delim(inputlocation,stringsAsFactors=FALSE)

## Load Molecule Signature Data Base ####
SigDB <- loadsigdb(sigdbfile)

### Define the dataset and transform Uniprot IDs to SYMBOL IDs, creation of data matrix for Limma fit  ####
proteinids <- mat$Name

x = AnnotationDbi::select(org.Hs.eg.db, keys(org.Hs.eg.db), c("SYMBOL","UNIPROT"))

genesymbols = lapply(proteinids,function(protids) unique(x[match(protids,x$UNIPROT),"SYMBOL"])) ## we map the proteinids to gene ids (human)

genesymbols = sapply(genesymbols, function(ids) ids[!is.na(ids)][1]) ## we remove NAs and pick the first gene symbol per row of measurements

tail(genesymbols) ## to see how genesymbols looks like now

mat = mat[,2:ncol(mat)] ## we now remove the complex annotation

## Remove Proteins with NAs

naremove = rowSums(is.na(mat))<1
mat = mat[naremove,]
genesymbols = genesymbols[naremove]

## Since there are multiple measurements with the same gene symbol, we aggregate protein expression in mat by genesymbols

mat = aggregate(mat,by=list(genesymbols),FUN=median)

## Afterwards some clean up

rownames(mat) = mat[,1]
mat = data.matrix(mat[,-1])


## Define experimental design ####

experiment <- c(rep(2,condition1),rep(1,condition2))
design <- model.matrix(~experiment)


## Limma lmFit ####
fit <- lmFit(mat, design)
fit <- eBayes(fit)
fitable <- topTable(fit, number = dim(fit$t)[1])

## Preparation of the msigdb for this dataset using the selected gene set ####

geneset = formatsigdbforroast(mat,SigDB)[[db]]

## Get the average z-scores per geneset ####
zmat <- getavgstatisticspergeneset(fit,geneset)

## The z-scores per protein ####
zscoresproteins <- getzscoresperprotein(fit, geneset)


## Now calculating significance of enrichment with roast ####
n.rot <- 9999
roastres1 <- roast(mat,contrast=2,design=design, nrot = n.rot, index=geneset,block=NULL,correlation = NULL) ## The comparison of interest (Condition1 vs Condition2)
roastres1 = roastres1[order(rownames(roastres1)),]

## Preparing data for rigde density plot ####

usedcategories <- SigDB[str_which(names(SigDB), names(db))]
usedcategories <- lapply(usedcategories, function(x) x[-c(1,2)])

usedcategories.df <- list2df(usedcategories) %>% 
            dplyr::select(category = X2, proteinid = X1) %>% 
            mutate(category = str_remove(category, paste0(names(db),"_"))) %>% 
            mutate(category = str_trim(category))

roastout.df <- tibble(category = row.names(roastres1),
                   fdr = roastres1$FDR,
                   fdrmixed = roastres1$FDR.Mixed,
                   pvalue = roastres1$PValue,
                   pvaluemixed = roastres1$PValue.Mixed,
                   ngenes = roastres1$NGenes)


if(dim(roastout.df)[1] > 50){
   
   roastout.dffil <- dplyr::top_n(roastout.df,-pvalue, n = 50) %>%
                  dplyr::filter(ngenes > 2)
   
} else {
   roastout.dffil <- roastout.df
}


fitable.df <- tibble(proteinid = row.names(fitable),
                        log2fc = fitable$logFC)

zscoresperprot.df <- tibble(proteinid = row.names(zscoresproteins),
                            zscore = zscoresproteins[,2])

lmfit_w_genesetcat <- left_join(fitable.df, usedcategories.df)

lmfit_w_genesetcat_w_zscores <- left_join(lmfit_w_genesetcat, zscoresperprot.df)

lmfit_w_genesetcat_w_roastout <- left_join(lmfit_w_genesetcat_w_zscores, roastout.dffil) %>% na.omit()


### List of plots ####

list_ridgeplots <- list()

### Plot with Log2FC x FDR ####

ridgeplot_log2fc_fdr <- ggplot(data = lmfit_w_genesetcat_w_roastout,
                     aes(x = log2fc, y = category, fill = fdr))+
      scale_fill_gradient(low = "#477af8", high = "#ff3333", name = "FDR")+
      geom_density_ridges()+
      xlab("Log2(Fold-change)")+ 
      ylab(paste0(names(db)," ","Geneset"))+
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5,size = 10, color = "black"),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1.5),
            axis.title=element_text(size=10,face="bold"))

list_ridgeplots[[1]] <- ridgeplot_log2fc_fdr

names(list_ridgeplots)[1] <- "ridgeplot_log2fc_fdr"

### Plot with Log2FC x FDR mixed ####

ridgeplot_log2fc_fdrmixed <- ggplot(data = lmfit_w_genesetcat_w_roastout,
                               aes(x = log2fc, y = category, fill = fdrmixed))+
      scale_fill_gradient(low = "#477af8", high = "#ff3333", name = "FDR mixed")+
      geom_density_ridges()+
      xlab("Log2(Fold-change)")+ 
      ylab(paste0(names(db)," ","Geneset"))+
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5,size = 10, color = "black"),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1.5),
            axis.title=element_text(size=10,face="bold"))

list_ridgeplots[[2]] <- ridgeplot_log2fc_fdrmixed

names(list_ridgeplots)[2] <- "ridgeplot_log2fc_fdrmixed"

### Plot with Log2FC x P-value ####

ridgeplot_log2fc_pvalue <- ggplot(data = lmfit_w_genesetcat_w_roastout,
                                    aes(x = log2fc, y = category, fill = pvalue))+
      scale_fill_gradient(low = "#477af8", high = "#ff3333", name = "P-value")+
      geom_density_ridges()+
      xlab("Log2(Fold-change)")+ 
      ylab(paste0(names(db)," ","Geneset"))+
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5,size = 10, color = "black"),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1.5),
            axis.title=element_text(size=10,face="bold"))

list_ridgeplots[[3]] <- ridgeplot_log2fc_pvalue

names(list_ridgeplots)[3] <- "ridgeplot_log2fc_pvalue"

### Plot with Log2FC x P-value mixed ####

ridgeplot_log2fc_pvaluemixed <- ggplot(data = lmfit_w_genesetcat_w_roastout,
                                  aes(x = log2fc, y = category, fill = pvaluemixed))+
      scale_fill_gradient(low = "#477af8", high = "#ff3333", name = "P-value mixed")+
      geom_density_ridges()+
      xlab("Log2(Fold-change)")+ 
      ylab(paste0(names(db)," ","Geneset"))+
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5,size = 10, color = "black"),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1.5),
            axis.title=element_text(size=10,face="bold"))

list_ridgeplots[[4]] <- ridgeplot_log2fc_pvaluemixed

names(list_ridgeplots)[4] <- "ridgeplot_log2fc_pvaluemixed"


### Plot with Log2FC x N Genes ####

ridgeplot_log2fc_nproteins <- ggplot(data = lmfit_w_genesetcat_w_roastout,
                                       aes(x = log2fc, y = category, fill = ngenes))+
      scale_fill_gradient(low = "#477af8", high = "#ff3333", name = "N proteins")+
      geom_density_ridges()+
      xlab("Log2(Fold-change)")+ 
      ylab(paste0(names(db)," ","Geneset"))+
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5,size = 10, color = "black"),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1.5),
            axis.title=element_text(size=10,face="bold"))

list_ridgeplots[[5]] <- ridgeplot_log2fc_nproteins

names(list_ridgeplots)[5] <- "ridgeplot_log2fc_nproteins"

### Plot with protein z-scores x FDR ####

ridgeplot_zscores_fdr <- ggplot(data = lmfit_w_genesetcat_w_roastout,
                                     aes(x = zscore, y = category, fill = fdr))+
      scale_fill_gradient(low = "#477af8", high = "#ff3333", name = "FDR")+
      geom_density_ridges()+
      xlab("Z-scores of protein expression")+ 
      ylab(paste0(names(db)," ","Geneset"))+
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5,size = 10, color = "black"),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1.5),
            axis.title=element_text(size=10,face="bold"))

ridgeplot_zscores_fdr

list_ridgeplots[[6]] <- ridgeplot_zscores_fdr

names(list_ridgeplots)[6] <- "ridgeplot_zscores_fdr"

### Plot with protein z-scores x FDR mixed ####

ridgeplot_zscores_fdrmixed <- ggplot(data = lmfit_w_genesetcat_w_roastout,
                                aes(x = zscore, y = category, fill = fdrmixed))+
      scale_fill_gradient(low = "#477af8", high = "#ff3333", name = "FDR mixed")+
      geom_density_ridges()+
      xlab("Z-scores of protein expression")+ 
      ylab(paste0(names(db)," ","Geneset"))+
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5,size = 10, color = "black"),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1.5),
            axis.title=element_text(size=10,face="bold"))

ridgeplot_zscores_fdrmixed

list_ridgeplots[[7]] <- ridgeplot_zscores_fdrmixed

names(list_ridgeplots)[7] <- "ridgeplot_zscores_fdrmixed"

### Plot with protein z-scores x P-value ####

ridgeplot_zscores_pvalue <- ggplot(data = lmfit_w_genesetcat_w_roastout,
                                     aes(x = zscore, y = category, fill = pvalue))+
      scale_fill_gradient(low = "#477af8", high = "#ff3333", name = "P-value")+
      geom_density_ridges()+
      xlab("Z-scores of protein expression")+ 
      ylab(paste0(names(db)," ","Geneset"))+
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5,size = 10, color = "black"),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1.5),
            axis.title=element_text(size=10,face="bold"))

ridgeplot_zscores_pvalue


list_ridgeplots[[8]] <- ridgeplot_zscores_pvalue

names(list_ridgeplots)[8] <- "ridgeplot_zscores_pvalue"

### Plot with protein z-scores x P-value mixed ####

ridgeplot_zscores_pvaluemixed <- ggplot(data = lmfit_w_genesetcat_w_roastout,
                                   aes(x = zscore, y = category, fill = pvaluemixed))+
      scale_fill_gradient(low = "#477af8", high = "#ff3333", name = "P-value mixed")+
      geom_density_ridges()+
      xlab("Z-scores of protein expression")+ 
      ylab(paste0(names(db)," ","Geneset"))+
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5,size = 10, color = "black"),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1.5),
            axis.title=element_text(size=10,face="bold"))

ridgeplot_zscores_pvaluemixed


list_ridgeplots[[9]] <- ridgeplot_zscores_pvaluemixed

names(list_ridgeplots)[9] <- "ridgeplot_zscores_pvaluemixed"


### Plot with protein z-scores x N proteins ####

ridgeplot_zscores_nproteins <- ggplot(data = lmfit_w_genesetcat_w_roastout,
                                        aes(x = zscore, y = category, fill = ngenes))+
      scale_fill_gradient(low = "#477af8", high = "#ff3333", name = "N Proteins")+
      geom_density_ridges()+
      xlab("Z-scores of protein expression")+ 
      ylab(paste0(names(db)," ","Geneset"))+
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5,size = 10, color = "black"),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=1.5),
            axis.title=element_text(size=10,face="bold"))

ridgeplot_zscores_nproteins

list_ridgeplots[[10]] <- ridgeplot_zscores_nproteins

names(list_ridgeplots)[10] <- "ridgeplot_zscores_nproteins"

### Save hi-res PDF plots for exploration ####

if(dir.exists(here::here("ROAST_output")) == FALSE){dir.create(here::here("ROAST_outout"))}

if(dim(roastout.df)[1] > 30){
   height <- 180
} else {
   height <- 145
}

if(names(db) == "REACTOME"){
   width <- 350
} else {
   width <- 200
}

for (i in 1:length(list_ridgeplots)){
      ggsave(filename = here::here(paste0("ROAST_outout/ROAST_",names(db),"_geneset_",names(list_ridgeplots)[i],".pdf")), 
             plot = list_ridgeplots[[i]], device = "pdf", 
             width = width, height = height, units = 'mm', dpi = 300)
}

### Save TSV file of LIMMA and ROAST results ####

write_tsv(lmfit_w_genesetcat_w_roastout,
          path = here::here(paste0("ROAST_outout/ROAST_tab_output_",names(db),"_geneset.tsv")))

limma_table <- mutate(fitable,
                      proteinid = row.names(fitable))

write_tsv(limma_table,
          path = here::here(paste0("ROAST_outout/LIMMA_fit_tab_output_",names(db),"_geneset.tsv")))




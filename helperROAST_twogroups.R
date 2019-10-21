loadsigdb = function(sigdbfile)
{
      SigDB <- scan(sigdbfile,what='character',sep='\n')
      SigDB <- strsplit(SigDB,split='\t')
      names(SigDB) <- lapply(SigDB,'[[',1)
      return(SigDB)
}

formatsigdbforroast = function(mat,sigdb)
{
      idx <- 1:nrow(mat)
      names(idx) <- rownames(mat)
      gene.roast <- lapply(sigdb,function(x) idx[x])
      gene.roast <- lapply(gene.roast, function(x) x[!is.na(x)]) 
      gene.roast = gene.roast[sapply(gene.roast, function(x) length(x)>0)]
      
      ## Create Subsets of Gene sets containing the KEGG, GO, REACTOME and HALLMARK sets
      kegg <- gene.roast[grepl('KEGG',names(gene.roast))]
      reactome <- gene.roast[grepl('REACTOME',names(gene.roast))]
      hallmark <- gene.roast[grepl('HALLMARK',names(gene.roast))]
      go <- gene.roast[grepl('GO',names(gene.roast))]
      naba <- gene.roast[grepl('NABA',names(gene.roast))]
      names(hallmark) <- gsub("HALLMARK_","",names(hallmark))
      names(kegg) <- gsub("KEGG_","",names(kegg))
      names(reactome) <- gsub("REACTOME_","",names(reactome))
      names(go) <- gsub("GO_","",names(go))
      names(naba) <- gsub("NABA_","",names(naba))
      return(list(GO=go,KEGG=kegg,REACTOME=reactome,HALLMARK=hallmark, NABA=naba))
}


getavgstatisticspergeneset = function(limmafit,geneset,ordinaryt=TRUE,zscore=TRUE)
{
      tmat = limmafit$t
      if (ordinaryt)
      {
            tmat <- limmafit$coef / limmafit$stdev.unscaled / limmafit$sigma
      }
      if (zscore)
      {
            tmat = zscoreT(tmat,df=limmafit$df.total)
            if (ordinaryt)
            {
                  tmat = zscoreT(tmat,df=limmafit$df.residual)
            }
      }
      return(t(sapply(geneset[order(names(geneset))], function(x) apply(tmat, 2, function(y) mean(y[x])))))
}

getzscoresperprotein = function(limmafit,geneset,ordinaryt=TRUE,zscore=TRUE)
{
      tmat = limmafit$t
      if (ordinaryt)
      {
            tmat <- limmafit$coef / limmafit$stdev.unscaled / limmafit$sigma
      }
      if (zscore)
      {
            tmat = zscoreT(tmat,df=limmafit$df.total)
            if (ordinaryt)
            {
                  tmat = zscoreT(tmat,df=limmafit$df.residual)
            }
      }
      return(tmat)
}

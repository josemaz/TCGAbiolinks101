require(NOISeq)


mypca <- function(dat1, namecol, fout = NULL){
  d2 <-  colData(dat1)[namecol]
  mydat = NOISeq::readData( assay(dat1) , factors = d2)
  myPCA = dat(mydat, type = "PCA", logtransf = F)
  if(!is.null(fout)){
    print(paste0("Writing in: ",fout))
    png(fout)
  }
  explo.plot(myPCA, factor = namecol, plottype = "scores")
  if(!is.null(fout))dev.off()
}
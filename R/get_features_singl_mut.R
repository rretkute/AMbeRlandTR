#' Get features based on singe mutation
#'
#' @param tree Lineage tree
#' @return  List containing feature matrices & matrix showing relationships between cells
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'
get.features.singl.mut<-function(tree){
  # Get values of barcodes, make numeric (values 0 and 1)
  barcodes = str_split(tree$tip.label,"_",simplify=T)[,2]
  n<-length(barcodes)
  m<-length(unlist(strsplit(barcodes[1], split = "")))
  xx <- sapply(1:n, function(a) as.numeric(unlist(strsplit(barcodes[a], split = ""))))
  xx<-t(xx)
  
  #  get  data for each tree level
  cell_ids = str_split(tree$tip.label,"_",simplify=T)[,1]
  tree$tip.label<-cell_ids
  n.cells<-length(cell_ids)
  ancs<-Ancestors(tree, 1:n.cells, "all")
  
  #  Get 3 feature matrices & relationships between cells
  ftr_f_0_0<-matrix(0, nrow=n, ncol=n)
  ftr_f_0_1<-matrix(0, nrow=n, ncol=n)
  ftr_f_1_1<-matrix(0, nrow=n, ncol=n)
  lvls.rel<-matrix(0, nrow=n*(n-1)/2, ncol=3)
  ii<-1
  for(i in 1:n){
    for( j in (i+1):n){
      if(j>i & j<=n){
        sm<-xx[i,] + xx[j,]
        sts<-sapply(0:2, function(a) length(which(sm==a)))
        ftr_f_0_0[i,j]<- ftr_f_0_0[j,i]<-sts[1]
        ftr_f_0_1[i,j]<-ftr_f_0_1[j,i]<-sts[2]
        ftr_f_1_1[i,j]<-ftr_f_1_1[j,i]<-sts[3]
        nd.i<-ancs[[i]]
        cm.node<-mrca.phylo(tree, c(i,j))
        lvl<-which(nd.i==cm.node)
        lvls.rel[ii,]<-c(i,j,lvl)
        ii<-ii+1
      }
    }
  }
  return(list(lvls.rel, ftr_f_0_0, ftr_f_0_1, ftr_f_1_1))
}

#' Simulate lineage tree
#'
#' @param targets Number of targets
#' @param d Number of divisions
#' @param mu Mutation rate
#' @param notmut Symbol for not-mutated state
#' @param mut Symbol for mutated state
#' @param plt Plot tree at each division (default =FALSE)
#' @param full.div Return all division/mutation history (default =FALSE)
#' @return Lineage tree, with id's containing mutation records
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'

sim_tree<-function(targets, d, mu, notmut, mut, plt=FALSE, full.div=FALSE){
  # First division 
  tree<-rtree(n=2)
  tree$edge.length<-NULL
  tps<-tree$tip.label
  ind<-1
  xx<-paste(c(ind,"_",rep(notmut, targets)), sep="", collapse="")
  ind<-ind+1
  barcode<-str_split(xx, "_", simplify=T)[1,2]
  for(id in 1:2){
    rec<-unlist(strsplit(barcode, split = ""))
    lc<-which(rec==notmut)
    for(i in 1:length(lc)){
      u <- runif(1,0,1)
      if(u<mu){ # Mutation
        rec[lc[i]]<-sample(mut,1)
      }
    }
    tree$tip.label[tree$tip.label==tps[id]]<-paste(c(ind,"_",rec), sep="", collapse="")
    ind<-ind+1
  }
  if(plt) plot(tree, show.tip.label = FALSE)
  if (full.div) TREE<-tree
  # Division 2+  dd<-3
  for(dd in 2:d){
    tps<-tree$tip.label
    for(j in 1:length(tps)){
      xx<-tps[j]
      node <- which(tree$tip.label==xx)
      tree <- bind.tip(tree, tip.label="new", where=node)
      barcode<-str_split(xx, "_", simplify=T)[1,2]
      rec<-unlist(strsplit(barcode, split = ""))
      for(id in 1:2){
        lc<-which(rec==notmut)
        for(i in 1:length(lc)){
          u <- runif(1,0,1)
          if(u<mu){ # Mutation
            rec[lc[i]]<-sample(mut,1)
          }
        }
        if(id==1){
          tree$tip.label[tree$tip.label==xx]<-paste(c(ind,"_",rec), sep="", collapse="")
          ind<-ind+1
        } else {
          tree$tip.label[tree$tip.label=="new"]<-paste(c(ind,"_",rec), sep="", collapse="")
          ind<-ind+1
        }
      }
    }
    if(plt) plot(tree, show.tip.label = FALSE)
    if (full.div) TREE<-c(TREE, tree)
  }
  # Count tips from 1
  tps<-tree$tip.label
  for(i in 1:length(tps)){
    barcode<-str_split(tps[i], "_", simplify=T)[1,2]
    tree$tip.label[tree$tip.label==tps[i]]<-paste(c(i,"_",barcode), sep="", collapse="")
  }
  if(plt) plotTree(tree, node.numbers=T)
  if (full.div){
    return(TREE)
  } else {
    return(tree)
  }
}

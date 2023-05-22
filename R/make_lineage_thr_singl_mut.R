#' Make a lineage tree based on a vector of thresholds and single mutation
#'
#' @param xx Barcode data
#' @param fit.model A list of fitted ML models for relationship between cells
#' @param thrs A vector of thresholds
#' @return  List containing feature matrices & matrix showing relationships between cells
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'
make.lineage.thr.singl.mut<-function(xx, ML.model, thrs){
  model.levels<-length(ML.model)
  cells<- xx[,1] #sapply(1:nrow(xx), function(a) paste0(xx[a,1],"_", xx[a,2]))
  n.cells<-length(cells)
  m<-length(unlist(strsplit(xx[1,2], split = "")))
  yy <- sapply(1:n.cells, function(a) as.numeric(unlist(strsplit(xx[a,2], split = ""))))
  yy<-t(yy)
  test.tw<-matrix(0, ncol=3, nrow=(n.cells-1)*n.cells)
  colnames(test.tw)<-c("f0_0", "f0_1", "f1_1")
  ids<-matrix(0, ncol=3, nrow=(n.cells-1)*n.cells)
  colnames(ids)<-c("no", "cell1", "cell2")
  ii<-1
  for(j1 in 1:n.cells){
    if(j1<n.cells){
      for(j2 in (j1+1):n.cells){
        sm<-yy[j1,] + yy[j2,]
        f<-sapply(0:2, function(a) length(which(sm==a)))
        test.tw[ii,]<- f
        ids[ii,]<-c(ii,  j1, j2)
        test.tw[ii+1,]<- f
        ids[ii+1,]<-c(ii+1,  j2, j1)
        ii<-ii+2
      }
    }
  }
  test.tw<-as.data.frame(test.tw)
  ids<-as.data.frame(ids)
  pred<-predict(ML.model[[1]], test.tw,  type='response')
  xy<-cbind(ids, pred)
  xy<-xy[order(-xy$pred),]
  GRP<-list()
  SNGL<-list()
  LEVELS<-list()
  RECORD<-list()
  CELLS<-list()
  grp<-data.frame(cell1=c(), cell2=c())
  levels<-list()
  
  while(nrow(xy)>0){
    if(xy$pred[1]>thrs[1]){  # # # # # # # # # #
      grp<-rbind(grp, data.frame(cell1=xy$cell1[1], cell2=xy$cell2[1]))
      levels[[length(levels)+1]]<-c(xy$cell1[1], xy$cell2[1])
    } 
    xy<-xy[-c(which(xy$cell1==xy$cell1[1]), which(xy$cell1==xy$cell2[1]), 
              which(xy$cell2==xy$cell1[1]), which(xy$cell2==xy$cell2[1])),]
  }
  
  sngl<-c()
  for(i in 1:length(cells)){
    if(nrow(grp)>0){
      if(!(i %in% c(grp[,1],grp[,2]))){
        sngl<-c(sngl, i)
        levels[[length(levels)+1]]<-i
      }
    } else
    {
      sngl<-c(sngl, i)
      levels[[length(levels)+1]]<-i
    }
  }
  GRP[[1]]<-grp
  SNGL[[1]]<-sngl
  LEVELS[[1]]<-levels
  CELLS[[1]]<-LEVELS[[1]]
  record<-list(rep(NA, length(levels)))
  if(length(levels)>0){
    for(i in 1:length(levels)){
      tmp<-levels[[i]]
      if(length(tmp)==2){
        record[[i]]<-paste0('(', as.character(cells[tmp[1]]),',', as.character(cells[tmp[2]]),')')
      } else {
        record[[i]]<-paste0('(', as.character(cells[tmp]),',)')
      }
    }
  }
  RECORD[[1]]<-record
  
  ind.cells<-levels
  curr.level<-3
  for(curr.level in 2:length(thrs)){  
    if(length(record)>2 & nrow(grp)>0){
      if(curr.level<=model.levels){
          pred<-predict(ML.model[[curr.level]], test.tw,  type='response')
      }  else {
        pred<-predict(ML.model[[model.levels]], test.tw,  type='response')
      }
      xy<-cbind(ids, pred)
      xy.M<-matrix(0, ncol = n.cells, nrow=n.cells)
      for(a in 1:nrow(xy)) xy.M[xy$cell1[a], xy$cell2[a]]<-xy$pred[a]
      
      xy.2<-matrix(0, nrow= length(ind.cells)*length(ind.cells)-length(ind.cells), ncol=3)
      colnames(xy.2)  <-c("cell1", "cell2", "pred")
      ii<-1
      for(i in 1:length(ind.cells)){
        if(i<length(ind.cells)){
          t1<-ind.cells[[i]]
          for(j in (i+1):length(ind.cells)){
            t2<- ind.cells[[j]]
            xy.2[ii,]<-c(i, j, pred=max(xy.M[t1, t2]))
            xy.2[ii+1,]<-c(j, i, pred=max(xy.M[t1, t2]))
            ii<-ii+2
          }
        }
      }
      xy.2<-as.data.frame(xy.2)
      xy.2<-xy.2[order(-xy.2$pred),]
      
      grp<-data.frame(cell1=c(), cell2=c())
      levels<-list()
      ind.cells<-list()
      while(nrow(xy.2)>0){
        if(xy.2$pred[1]>thrs[curr.level]){  # # # # # # # # # #
          grp<-rbind(grp, data.frame(cell1=xy.2$cell1[1], cell2=xy.2$cell2[1]))
          levels[[length(levels)+1]]<-c(xy.2$cell1[1], xy.2$cell2[1])
          ind.cells[[length(ind.cells)+1]]<-c(CELLS[[curr.level-1]][[xy.2$cell1[1]]], CELLS[[curr.level-1]][[xy.2$cell2[1]]])
        } 
        xy.2<-xy.2[-unique(c(which(xy.2$cell1==xy.2$cell1[1]), which(xy.2$cell1==xy.2$cell2[1]), 
                             which(xy.2$cell2==xy.2$cell1[1]), which(xy.2$cell2==xy.2$cell2[1]))),]
      }
      sngl<-c()
      for(i in 1:length(LEVELS[[curr.level-1]])){
        if(nrow(grp)>0){
          if(!(i %in% c(grp[,1], grp[,2]))){
            sngl<-c(sngl, i)
            levels[[length(levels)+1]]<-i
            ind.cells[[length(ind.cells)+1]]<-c(CELLS[[curr.level-1]][[i]])
          }
        } else {
          sngl<-c(sngl, i)
          levels[[length(levels)+1]]<-i
          ind.cells[[length(ind.cells)+1]]<-c(CELLS[[curr.level-1]][[i]])
        }
      }
      GRP[[curr.level]]<-grp
      SNGL[[curr.level]]<-sngl
      LEVELS[[curr.level]]<-levels
      CELLS[[curr.level]]<-ind.cells
      if(nrow(grp)>0){
        rr<-RECORD[[curr.level-1]]
        record<-list(rep(NA, length(levels)))
        if(length(levels)>0){
          for(i in 1:length(levels)){
            tmp<-levels[[i]]
            if(length(tmp)==2){
              record[[i]]<-paste0('(', as.character(rr[[tmp[1]]]),',', as.character(rr[[tmp[2]]]),')')
            } else {
              record[[i]]<-paste0('(', as.character(rr[[tmp]]),',)')
            }
          }
        }
        RECORD[[curr.level]]<-record
      }
    }
    
  }
  lng<-"(";
  if(length(record)>1){
    for(r in 1:(length(record)-1)){
      lng<-paste0(lng, as.character(record[[r]]),',')
    }}
  lng<-paste0(lng, as.character(record[[length(record)]]),");")
  return(lng)
}

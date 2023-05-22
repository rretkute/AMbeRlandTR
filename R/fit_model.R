#' Fit model using machine learning
#'
#' @param training_data A list of ineage tree
#' @return  ML model for each relationship level
#' @author Renata Retkute, \email{r.retkute@@yahoo.com}
#' @export
#'
fit.model<-function(training_data){
  n.train<-length(training_data)
  fts<-vector(mode = "list", length = n.train)
  for(ii in 1:n.train){
    sts<-get.features.singl.mut(training_data[[ii]])
    fts[[ii]]<-sts
  }
  max.lvls.rel<-max(sapply(1:n.train, function(a) max(c(fts[[a]][[1]][,3]))))
  ftm<-vector(mode = "list", length = max.lvls.rel)
  for(level in 1:max.lvls.rel){
      train.data<-data.frame(f0_0=c(), f0_1=c(), f1_1=c(), status=c())
      for(ii in 1:n.train){
        lvls.rel<-fts[[ii]][[1]]
        ftr_f_0_0<-fts[[ii]][[2]] 
        ftr_f_0_1<-fts[[ii]][[3]] 
        ftr_f_1_1<-fts[[ii]][[4]]
        colnames(lvls.rel)<-c("i","j","lvl")
        lvls.rel<-as.data.frame(lvls.rel)
        wh<-which(lvls.rel$lvl<=level)
        if(length(wh)>0){
          f1<-sapply(1:length(wh), function(a) ftr_f_0_0[lvls.rel$i[wh[a]], lvls.rel$j[wh[a]]])
          f2<-sapply(1:length(wh), function(a) ftr_f_0_1[lvls.rel$i[wh[a]], lvls.rel$j[wh[a]]])
          f3<-sapply(1:length(wh), function(a) ftr_f_1_1[lvls.rel$i[wh[a]], lvls.rel$j[wh[a]]])
          train.data.1<-cbind(f1, f2, f3, rep(1, length(wh)))
          colnames(train.data.1)<-c("f0_0", "f0_1", "f1_1", "status")
          train.data.1<-unique(train.data.1)
          train.data<-rbind(train.data, train.data.1)
        }
        wh<-which(lvls.rel$lvl>level & lvls.rel$lvl<11)
        if(length(wh)>0){
          f1<-sapply(1:length(wh), function(a) ftr_f_0_0[lvls.rel$i[wh[a]], lvls.rel$j[wh[a]]])
          f2<-sapply(1:length(wh), function(a) ftr_f_0_1[lvls.rel$i[wh[a]], lvls.rel$j[wh[a]]])
          f3<-sapply(1:length(wh), function(a) ftr_f_1_1[lvls.rel$i[wh[a]], lvls.rel$j[wh[a]]])
          train.data.0<-cbind(f1, f2, f3, rep(0, length(wh)))
          colnames(train.data.0)<-c("f0_0", "f0_1", "f1_1", "status")
          train.data.0<-unique(train.data.0)
          train.data<-rbind(train.data, train.data.0)
        }
      }
      train.data<-unique(train.data)
      ftm.i <- gbm(
        formula = status ~ .,
        data = as.data.frame(train.data),
        distribution = "bernoulli",
        verbose = FALSE
      )
    ftm[[level]]<-ftm.i
  }
  return(ftm)
}


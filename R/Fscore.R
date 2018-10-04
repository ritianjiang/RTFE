Fscore_each<-function(feature,Pos){
  pos_f<-feature[Pos];neg_f<-feature[-Pos]
  upper<-(mean(pos_f)-mean(feature))^2 + (mean(neg_f)-mean(feature))^2
  down<-(1/(length(Pos)-1))*sum((pos_f - mean(pos_f))^2) +
    (1/(length(feature)-length(Pos)-1))*sum((neg_f - mean(neg_f))^2)
  return(upper/down)
}
##' The function implement the F-score method to select the features. It returns a
##'
##' @title Fscore
##' @export
##' @param features A numeric matrix/data.frame which denotes features. Each row is a sample
##' @param Pos A numeric vector contains the positive index.
##' @return A number which is the fscore
##' @description The formula is from doi:10.3389/fmicb.2018.02174. A high f-score means
##' feature has more powerful discrimination ability.


Fscore<-function(features,Pos){
  features<-as.matrix(features)
  return(apply(features,2,Fscore_each,Pos=Pos))
}

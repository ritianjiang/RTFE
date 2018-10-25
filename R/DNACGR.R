##' The function implement the Chou's Chaos Game Feature Extraction method
##'
##' @title DNACGR
##' @export
##' @param Seq String. Which denotes the DNA sequence
##' @param Xdivide numeric. denotes the number of X to be splitted. default:4
##' @param Ydivide numeric. denotes the number of Y to be splitted. default:4
##' @return A number which is the fscore
##' @description The formula is from J.Jia's JTB paper about PPI in 2019. In the
##' square, A=(1,1),T=(0,1),C=(0,0),G=(1,0)
##' @return a numeric vector with length of Xdivide*Ydivide

DNACGR<-function(Seq,Xdivide=4,Ydivide=4){
  resX<-rep(-1,str_length(Seq))
  resY<-rep(-1,str_length(Seq))
  tempX<-0.5;tempY<-0.5
  for(i in 1:length(resX)){
    nuc<-str_sub(string = Seq,i,i)
    if(nuc == "A"){tempX<-0.5*(1+tempX);tempY<-0.5*(1+tempY)}
    else if(nuc == "T"){tempX<-0.5*(0+tempX);tempY<-0.5*(1+tempY)}
    else if(nuc == "G"){tempX<-0.5*(1+tempX);tempY<-0.5*(0+tempY)}
    else if(nuc == "C"){tempX<-0.5*(0+tempX);tempY<-0.5*(0+tempY)}
    resX[i]<-tempX;resY[i]<-tempY
  }
  (table(cut(x = resX,breaks = Xdivide,labels = F),
        cut(x = resY,breaks = Ydivide,labels = F)) %>% as.list() %>% unlist)/length(resX)
}

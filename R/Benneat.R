##' A class used to store the benchmark dataset with the same length so called neat.
##' @title Benneat
##'
##' @slot Seqs A vector contains the sequences
##' @slot labs A vector contains the labels like "methylated"/"unmethylated"
##' @slot seq_num A numeric denotes the num of seqlist
##' @slot lab_num A table contains the num of each Label
##' @slot seq_length A numerci denotes the length of each sequence
##' @docType class
##' @export

setClass("Benneat",
         slots = list(Seqs = "vector",labs = "vector",
                      seq_num ="numeric",lab_num="table",seq_length="numeric"))

##' Make Benneat Class by TWO parameters
##' @title makeBenneat
##'
##' @param  Seqs A vector contains the sequences
##' @param labs A vector contains the labels like "methylated"/"unmethylated"
##' @docType method
##' @export
##' @return A Benneat Class

setGeneric("makeBenneat",function(Seqs,labs) standardGeneric("makeBenneat"))

setMethod("makeBenneat",signature(Seqs="vector",labs="vector"),
          function(Seqs,labs){
            if(lapply(Seqs, str_length) %>% unique %>% length()!=1) stop("Length should be identical")
            if(length(Seqs)!=length(labs)) stop("Lab nums is not the same as seq nums.") #Set the validity manually
            seq_num<-length(Seqs);lab_num<-table(labs)
            seq_length<-Seqs[1] %>% str_length() #Get the essensial elements for Benneat Class
            result<-new("Benneat",Seqs=Seqs,labs=labs,# creat the class
                        seq_num=seq_num,lab_num=lab_num,seq_length=seq_length)
            return(result)
          })

setMethod("show", "Benneat",##Just adjust the show pattern for Benneat Class
           function(object) {
             cat(paste0("A benchmark set contains ",object@seq_length,"-length sequence \n"))
             #print(" ")
             print(object@lab_num)
           }
)

##'
##' @title PSTNPss
##'
##' @slot object A Benneat object
##' @slot lableA character,Label as the POSITIVE label
##' @slot lableB character,Label as the POSITIVE label
##' @docType method
##' @export
##' @return A matrix contains ALL sequences' PSTNPss features,each row denotes
##' a sequence. The Sequence is sorted by its order in Benneat object. Each col
##' is sorted by AAA,AAC,AAG............TTT

setGeneric("PSTNPss",function(object,lableA,lableB) standardGeneric("PSTNPss"))

setMethod("PSTNPss","Benneat",function(object,lableA,lableB){
  seq_A<-object@Seqs[object@labs == lableA] #Get thesequence which is labeled as lableA
  seq_B<-object@Seqs[object@labs == lableB] #Get thesequence which is labeled as lableA
  Dic<-expand.grid(c("A","T","C","G"),c("A","T","C","G"),c("A","T","C","G"),
                   stringsAsFactors = F) %>%
    apply(1,paste,collapse="") %>% sort #Create the name vector for 3-mer
  F_A<-matrix(-1,nrow=64,ncol=object@seq_length-2)
  F_B<-matrix(-1,nrow=64,ncol=object@seq_length-2) #Calculate the frequency of each 3-mer in each position in both label-classes
  for(i in 1:object@seq_length-2){
    F_A[,i]<-lapply(seq_A, str_sub,start=i,end=i+2) %>%
      do.call(what = "rbind") %>%factor(levels = Dic) %>% table / length(seq_A)
  }
  for(i in 1:object@seq_length-2){
    F_B[,i]<-lapply(seq_B, str_sub,start=i,end=i+2) %>%
      do.call(what = "rbind") %>%factor(levels = Dic) %>% table / length(seq_A)
  }
  Z<- F_A - F_B;rownames(Z)<-Dic #create the Z matrix which contains z_{i,j}. Details and meanings is mentioned in Zhang Chun-Ting Lab's 2017 paper
  P<-matrix(-1,nrow=object@seq_num,ncol=ncol(Z)) #get the Features
  for(i in 1:ncol(Z)){
    Z[lapply(object@Seqs, str_sub,start=i,end=i+2) %>% do.call(what="rbind"),i]->P[,i]
  }
  return(P %>% as.data.frame())
})

##' @title getNCC
##'
##' @param object A Benneat object
##' @param seqType The sequence is DNA or RNA. It should be either 'DNA' or 'RNA'
##' @param n The k of k-mer
##' @docType method
##' @export
##' @return A matrix contains ALL sequences' NCC features,each row denotes
##' a sequence. The Sequence is sorted by its order in Benneat object. if you run
##' getNCC(xxx,"DNA",2), then each col will be sorted according to AA,AC,AG,AT...
##' TA,TC,TG,TT.
setGeneric("getNCC",function(object,seqType,n) standardGeneric("getNCC"))
setMethod("getNCC","Benneat",function(object,seqType,n){
  lapply(object@Seqs,NCC,n=n,seqType=seqType) %>% do.call(what="rbind") %>%
    as.data.frame() %>% return
})






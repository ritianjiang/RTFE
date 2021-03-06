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
##' @docType methods
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

##' @title PSTNPss
##'
##' @param object A Benneat object
##' @param lableA character,Label as the POSITIVE label
##' @param lableB character,Label as the POSITIVE label
##' @docType methods
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
##' @docType methods
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

##' @title PSTNPds
##' @description <PSTNPds is the double-strand version of PSTNPss>
##' Feature PSTNPds
##'using a statistical strategy based on
##'double-stranded characteristics of DNA according to
##'complementary base pairing, so they have more evident
##'statistical features. At this point, we deem A and T as
##'identical, the same to C and G. Thus, for every sample,
##'it can be converted into a sequence contained A and T
##'only.
##' @param object A Benneat object
##' @param lableA character,Label as the POSITIVE label
##' @param lableB character,Label as the POSITIVE label
##' @docType method
##' @export
##' @return A matrix contains ALL sequences' PSTNPds features,each row denotes
##' a sequence. The Sequence is sorted by its order in Benneat object. Columns
##' are sorted by AAA,AAC,ACA,ACC...CCC <8 columns in total>.

setGeneric("PSTNPds",function(object,lableA,lableB) standardGeneric("PSTNPds"))

setMethod("PSTNPds","Benneat",function(object,lableA,lableB){
  newSeq<-str_replace_all(string = object@Seqs,pattern = "T",replacement = "A")
  newSeq<-str_replace_all(string = newSeq,pattern = "G",replacement = "C")
  seq_A<-newSeq[object@labs == lableA] #Get thesequence which is labeled as lableA
  seq_B<-newSeq[object@labs == lableB] #Get thesequence which is labeled as lableA
  Dic<-expand.grid(c("A","C"),c("A","C"),c("A","C"),
                   stringsAsFactors = F) %>%
    apply(1,paste,collapse="") %>% sort #Create the name vector for 3-mer
  F_A<-matrix(-1,nrow=8,ncol=object@seq_length-2)
  F_B<-matrix(-1,nrow=8,ncol=object@seq_length-2) #Calculate the frequency of each 3-mer in each position in both label-classes
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
    Z[lapply(newSeq, str_sub,start=i,end=i+2) %>% do.call(what="rbind"),i]->P[,i]
  }
  return(P %>% as.data.frame())
})

##' @title getPseKNC
##'
##' @param object A Benneat object
##' @docType methods
##' @export
##' @return A matrix contains ALL sequences' PseKNC features,each row denotes
##' a sequence. The Sequence is sorted by its order in Benneat object.
setGeneric("getPseKNC",function(object) standardGeneric("getPseKNC"))
setMethod("getPseKNC","Benneat",function(object){
  lapply(object@Seqs,PseKNC) %>% do.call(what="rbind") %>%
    as.data.frame() %>% return
})

##' @title getEIIP
##'
##' @param object A Benneat object
##' @docType methods
##' @export
##' @return A matrix contains ALL sequences' EIIP features,each row denotes
##' a sequence. The Sequence is sorted by its order in Benneat object.
setGeneric("getEIIP",function(object,EIIP=c(0.1260,0.1335,0.0806,0.1340)) standardGeneric("getEIIP"))
setMethod("getEIIP","Benneat",function(object,EIIP=c(0.1260,0.1335,0.0806,0.1340)){
  lapply(object@Seqs,EIIP) %>% do.call(what="rbind") %>%
    as.data.frame() %>% return
})

##' @title getT2PseKNC
##' @description Please make sure the sequence type of your data. The default file is 6 features of double-strand B-DNA
##' @param object A Benneat object
##' @param phychem_file A matrix;whose rownames are features and colnames are 2-mer
##' @param normalization Bool, if TRUE, the function will perform normalization. default:T
##' @param lambda numeric. The max number of sequence tiers. default:4
##' @docType methods
##' @export
##' @return A matrix contains ALL sequences' moditified T2PseKNC features,each row denotes
##' a sequence. The Sequence is sorted by its order in Benneat object. The return
##' value only contains the long-range part of T2PSeKNC. The detailed explanation
##' can be find in Hao Lin Groups' paper iTerm-PseKNC. This function only implement
##' the 2-mer situations. The ncol of result will be as lambda*(phechem_file feature nums)-length
setGeneric("getT2PseKNC",function(object,phychem_file=NULL,
                                  normalization=T,
                                  lambda = 4) standardGeneric("getT2PseKNC"))
##z nromalization function
znorm<-function(x){
  return((x - mean(x))/sd(x))
}
##T2PseKNC for one sequence
T2PseKNC<-function(Seq,phychem,lambda=4){
  mer2_1<-str_sub(Seq,start = 1:(str_length(Seq)-1),end=2:str_length(Seq))
  result<-c()
  for(i in 1:lambda){
    mer1<-mer2_1[1:(length(mer2_1)-i)]
    mer2<-mer2_1[(i+1):length(mer2_1)]
    total_mer<-phychem[,mer2]*phychem[,mer1]
    result<-c(result,apply(total_mer,1,mean) %>% as.numeric)
  }
  return(result)
}



setMethod("getT2PseKNC","Benneat",function(object,phychem_file = "aaa",
                                           normalization=T,
                                           lambda=4){
  if(phychem_file =="aaa"){ #The default diprogb file
    cat("You can use the default data in RTFE/Path/data/Rdata.rds")
    return()
  }
  else{phychem<-phychem_file %>% as.matrix}
  if(normalization == T){
    phychem<-apply(phychem,1,znorm) %>% t
  }

  lapply(object@Seqs,T2PseKNC,phychem=phychem,lambda=lambda) %>% do.call(what="rbind") %>%
    as.data.frame -> result
  return(result)
})

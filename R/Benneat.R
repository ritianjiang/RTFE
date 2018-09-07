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
            if(length(Seqs)!=length(labs)) stop("Lab nums is not the same as seq nums.")
            seq_num<-length(Seqs);lab_num<-table(labs)
            seq_length<-Seqs[1] %>% str_length()
            result<-new("Benneat",Seqs=Seqs,labs=labs,
                        seq_num=seq_num,lab_num=lab_num,seq_length=seq_length)
            return(result)
          })

setMethod("show", "Benneat",
           function(object) {
             cat(paste0("A benchmark set contains ",object@seq_length,"-length sequence \n"))
             #print(" ")
             print(object@lab_num)
           }
)

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

setGeneric("PSTNP",function(object,strand) standardGeneric("PSTNP"))
setMethod("PSTNP","Benneat",function(object,strand){

})

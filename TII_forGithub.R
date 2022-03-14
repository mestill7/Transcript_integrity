
R  
##version 3.5.2
library(GenomicFeatures)
library(GenomicAlignments)
library(biovizBase)
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg38")
library("biomaRt")
require(genoset)
require(PureCN)
require("parallel")

setwd("~/Transcript_integrity/")
message("Currently reading gene annotation")
txdb <- makeTxDbFromGFF("./gencode.v28.primary.annotation.gtf", format="gtf") ##To make TxDb from external GTF
exons.list <- exonsBy(txdb,by=c("tx"),use.names=TRUE)
inboth <- paste0("chr",c(1:22,"X","Y","M"))
exons.list <- keepSeqlevels(exons.list,inboth,pruning.mode="coarse")

##Get the txt names of interest
exons.flat.sep <- flatGrl(exons.list)
inboth <- paste0("chr",c(1:22,"X","Y","M"))
exons.flat.sep <- keepSeqlevels(exons.flat.sep,inboth,pruning.mode="coarse")

## Read in the samples to be examined
sample_list <- c("Test1","Test2")
ser.reads <- lapply(sample_list,function(x) 
	readGAlignments(paste0("./",x,"_dupsremoved.bam")))
ser.reads <- lapply(ser.reads,function(x) x[seqnames(x)%in%inboth])
getsubset <- function(x){
	seqlevels(x) <- inboth 
	return(x)
}
ser.reads <- lapply(ser.reads,getsubset)
ser.cover <- lapply(ser.reads,function(x) coverage(x))
ser.libsize <- unlist(lapply(ser.reads,function(x) length(x)))

## Generate the list of gene symbols to be examined
my_set <- c("GAPDH","ACTB","DNAI1")
length(my_set)

## Match the list of gene symbols to their ensembl id 
ensembl=useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
bm <- getBM(attributes=c("ensembl_gene_id","ensembl_transcript_id","ensembl_transcript_id_version","external_gene_name"),
	filters='external_gene_name', 
            values=as.character(my_set),mart=ensembl)
my_1k_txt <- bm$ensembl_transcript_id_version
my_olp_txt <- names(exons.list)[names(exons.list)%in%my_1k_txt]

gettot <- function(x){
  # x="ENST00000456328.2"
  mylist= unlist(exons.list[x])
  end(mylist[1])=max(end(mylist))
  mylist$txt_id <- x
  myout <- data.frame(mylist[1])
return(myout)
}
system.time(exons.flat <- lapply(my_olp_txt,gettot))
exons.flat <- do.call("rbind",exons.flat)
exons.flat <- as(exons.flat,"GRanges")  ##txt start to end  

##preset the exons.list and exons.flat for cluster looping
subexons <- function(x){
	return(list(exons.flat[exons.flat$txt_id==x],unlist(exons.list[x])))
}
sub.list <- lapply(my_olp_txt,subexons)
names(sub.list) <- my_olp_txt

k=10 ##bin length
read_size <- 50 ##42 for H5 reads
diff_spike <- 5 ##


testTII <- function(mytxt,iter,libsizes,myreads,mycoverage,myexons,Hsapiens){
	k=10 ##bin length in bp
	read_size <- 50 ## read length in bp
	diff_spike <- 5 ##
	y=mycoverage[[iter]]
    z=myreads[[iter]]
    libsize=libsizes[iter]
    # return(c(length(y),length(z),libsize))
    # return(myexons)
    ml= myexons[[1]][[2]]
  	mytxt=myexons[[1]][[1]]
	# return(ml)
	mytiles <- unlist(tile(ml,width = k))
	exons_total <- length(unique(queryHits(findOverlaps(z,ml,type="any"))))
	txt_total <- length(unique(queryHits(findOverlaps(z,mytxt,type="any"))))
	length_total <- sum(width(ml))
	mu <- exons_total*read_size/length_total
	inboth <- intersect(seqlevels(mytiles),names(y))
	cover <- y[inboth]
	res <- binnedAverage(mytiles,cover,"mean_cvg") #GenomicRanges function
	res$gc <- calcGC(res,Hsapiens,expand = 20) #Genoset function
	res$diff <- res$mean_cvg-mu
	res$scale <- scale(res$diff) #zscale #base function
	##spikyness measure
	diffarb <- length(res$diff[abs(res$diff)<diff_spike])/length(res)  
	#continuity measure
	start_cov <- c(1,5,8,seq(from=10,to=50,by=5),80)
	continuity <- unlist(lapply(start_cov,function(i) return(length(res[res$mean_cvg > i])/length(res))))
	norm_cov <- libsize*(c(1,5,8,seq(from=10,to=50,by=5),80)/10^6)
	continuity.norm <- unlist(lapply(norm_cov,function(i) return(length(res[res$mean_cvg > i])/length(res))))
	return(c(exons_total,txt_total,exons_total/txt_total,length(ml),length_total,exons_total/length_total*1000,diffarb,continuity,continuity.norm))
}

# require("parallel")
ml <- lapply(1:length(sersamples.p1),function(i){
  # i=1
  myclust <- makeCluster(10)
    testout <- parSapplyLB(myclust,my_olp_txt,function(mytxt,i,mainfun,libsize,myreads,mycoverage,myexons,Hsapiens)
          mainfun(mytxt,i,libsize,myreads,mycoverage,myexons[grep(mytxt,names(myexons))],Hsapiens),i,testTII,ser.libsize,ser.reads,ser.cover,sub.list,Hsapiens)
    x2 <- data.frame(t(testout))
    colnames(x2) <- c("exonic_reads","txt_reads","exonic_ratio","exon_count","length","reads_per_kb","uniform_5reads","min_1","min_5","min_8","min_10","min_15","min_20","min_25","min_30","min_35","min_40","min_45","min_50","min_80","norm_1","norm_5","norm_8","norm_10","norm_15","norm_20","norm_25","norm_30","norm_35","norm_40","norm_45","norm_50","norm_80")
    write.table(x2,paste0("./TII_output/",sample_list[i],".txt"),sep="\t")
  stopCluster(myclust)
  }
  )
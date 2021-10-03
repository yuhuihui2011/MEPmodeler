#' Predict MicroExons in a Plant Genome
#' 
#' MEPmod searches microexon-tags in 45 conserved microexon clusters in plants 
#' using gapped Position Weight Matrix (PWM). The only input file is plant 
#' genomic sequences or a plant genome. The sequences on both plus and 
#' minus strands will be scanned.
#' 
#' @param genome the path(s) to the fasta file(s) or a 'DNAStringSet' object. 
#' Any contig sequence < 10 kb will be excluded.
#' @param min.score	a character string containing a percentage specifying the 
#' minimum score of each exon block (e.g. \code{"80\%"}). This parameter will 
#' pass to \link[Biostrings]{matchPWM} from \bold{Biostrings} package.
#' @param include.intronLoss \code{TRUE} or \code{FALSE}. If \code{TRUE}, 
#' the microexons with any side of flanking intron loss will also be returned.
#' @param span the maximum spanning region of the microexon-tag (default: 20 kb).
#' @param min.intron minimum intron size
#' @param max.intron maximum intron size
#' @param cores number of cores to use. This parameter will pass to 
#' \link[parallel]{mclapply} as *mc.cores* from \bold{parallel} package.
#' 
#' @return A GRanges of microexons with 8 metadata columns:
#' \itemize{
#' \item NT, 108 bp DNA sequence of microexon-tag (72 bp for cluster 1 and 2).
#' \item AA, 36 aa protein sequence translated from NT (24 aa for cluster 1 and 2).
#' \item region, the span region of microexon-tag.
#' \item block.starts, the genomic coordinates of block start positions 
#' (exon, intron, exon, ...).
#' \item block.sizes, block sizes (exon, intron, exon, ...).
#' \item score, match score of the microexon-tag. 
#' Large score indicate high confidence.
#' \item isMicroExon, if the microexon has any flanking intron loss, FALSE; 
#' otherwise, TRUE.
#' \item cluster, which cluster this microexon-tag belongs to.
#' }
#' 
#' @importFrom methods is
#' @importFrom parallel mclapply
#' @import IRanges
#' @import GenomicRanges
#' @import Biostrings
#' 
#' @export MEPmod
#' 
#' @examples
#' suppressPackageStartupMessages(library(BSgenome))
#' if (!requireNamespace("BSgenome.Athaliana.TAIR.TAIR9", quietly = TRUE))
#'     BiocManager::install("BSgenome.Athaliana.TAIR.TAIR9")
#' library(BSgenome.Athaliana.TAIR.TAIR9)
#' 
#' Chr1<-BSgenome.Athaliana.TAIR.TAIR9[["Chr1"]]
#' Chr1<-DNAStringSet(Chr1)
#' names(Chr1)<-'Chr1'
#' x<-subseq(Chr1, 1, 1e6)
#' x
#' 
#' mep<-MEPmod(genome=x, cores=1)
#' 
#' mep
#' # GRanges object with 2 ranges and 8 metadata columns:
#' #       seqnames        ranges strand |                      NT                      AA
#' #          <Rle>     <IRanges>  <Rle> |          <DNAStringSet>           <AAStringSet>
#' #   [1]     Chr1 456757-456761      + | TTTGATGCTA...ACTGTCTGAC FDARTAWSQC...AFGAVESLSD
#' #   [2]     Chr1 601547-601560      - | AAAGTTGAAT...ACATTCAGAT KVEFKDNEWK...RLDCLLKHSD
#' #              region             block.starts   block.sizes     score isMicroExon   cluster
#' #           <IRanges>            <IntegerList> <IntegerList> <numeric>   <logical> <integer>
#' #   [1] 456595-456893 456595,456647,456757,...  52,110,5,...    0.9800        TRUE         8
#' #   [2] 601327-601700 601653,601561,601547,...  48,92,14,...    0.9241        TRUE        35
#' #   -------
#' #   seqinfo: 1 sequence from an unspecified genome; no seqlengths
#' 

###
MEPmod<-function(genome, min.score='80%',include.intronLoss=TRUE,
                 span=20000, min.intron=20, max.intron=10000, cores=16) {
    MEPdata<-MEPmodeler::MEPdata
    cat(paste(t0<-Sys.time(),'..... started run\n'))
    cat(paste(Sys.time(),'..... loading genome\n'))
    if (!is(genome, 'DNAStringSet')) {
        if (!is.character(genome)) {
            stop("'genome' must be the path(s) to the fasta file(s) or
                 a 'DNAStringSet' class" )
        }
        genome<-readDNAStringSet(genome)
        names(genome)<-sub("^(\\S+)\\s+.*","\\1",names(genome))
        }
    genome<-genome[width(genome)>=10000]
    prior.params=letterFrequency(genome[[1]],DNA_BASES,as.prob = TRUE)
    cat(paste(Sys.time(),'..... finished loading genome\n'))
    res<-GRanges()
    for(i in seq(along.with = MEPdata$blocks)) {
        cat(sprintf('%s ..... Cluster %d (size: %d; phase: %d; motif: %s)\n',
                    Sys.time(),i,MEPdata$cluster$size[i],MEPdata$cluster$phase[i],
                    MEPdata$cluster$motif[i]))
        res_i<-mapPWM(cons=MEPdata$matrix[[i]],exons=MEPdata$blocks[[i]],genome,
                      focus=MEPdata$cluster$me_order[i],min.score,
                      include.intronLoss,span, min.intron, max.intron,
                      prior.params=prior.params,cores = cores,check.params = FALSE)
        if (length(res_i)==0) next
        res_i$block.starts<-do.call(c,apply(res_i$block.starts, 1, IntegerList))
        res_i$block.sizes<-do.call(c,apply(res_i$block.sizes, 1, IntegerList))
        res_i$cluster<-i
        res<-suppressWarnings(append(res,res_i))
    }
    cat(paste(t1<-Sys.time(),'..... finished successfully\n'))
    cat(sprintf('### total microexon-tags found:\t%d\n', length(res)))
    cat(sprintf('### total time used:\t%s mins\n', 
                round(difftime(t1,t0,units = 'mins'),2)))
    res
}

#####
mapPWM<-function(cons,exons,genome,focus=2, min.score='80%',include.intronLoss=TRUE,
         span=20000, min.intron=20, max.intron=10000, 
         prior.params=letterFrequency(genome[[1]],DNA_BASES,as.prob = TRUE),
         cores=1, check.params=TRUE){
    n<-length(exons)
    if (n==3) {
        return(mapPWM.0(cons,exons,genome,min.score,include.intronLoss, span, 
                        min.intron, max.intron, prior.params, cores))
    }
    if (check.params) {
        if (n<3) stop("'exons' must be a numeric vector of length >= 3")
        if (sum(exons)!=ncol(cons))
            stop ("sum(exons)==ncol(cons) is not TRUE")
        if (!is(genome,'DNAStringSet')) stop("'genome' must be a DNAStringSet object")
        if (length(focus)!=1 & focus[1]<=1 & focus[1]>=n)
            stop("'focus' must be a numeric in (0,length(exons))!")
    }
    plus_intron1<-sprintf('^(GT[AG].{%d,%d}[CT]AG)%s$',min.intron,max.intron,
                          ifelse(include.intronLoss,"?",""))
    plus_intron2<-sprintf('^(G[TC].{%d,%d}AG)%s$',min.intron,max.intron,
                          ifelse(include.intronLoss,"?",""))
    minus_intron1<-sprintf('^(CT[AG].{%d,%d}[CT]AC)%s$',min.intron,max.intron,
                           ifelse(include.intronLoss,"?",""))
    minus_intron2<-sprintf('^(CT.{%d,%d}[GA]C)%s$',min.intron,max.intron,
                           ifelse(include.intronLoss,"?",""))
    
    pwm<-PWM(cons,prior.params=prior.params)
    for (i in seq(n)) {
        index<-seq.int(exons[i]) + cumsum(c(0,exons))[i]
        assign(paste0('pwm',i),PWM(cons[,index,drop=FALSE],prior.params=prior.params))
    }
    
    plus<-mclapply(genome, function(x) {
        # find matches for the firest and the last exon
        e1<-suppressWarnings(matchPWM(get('pwm1'),x,min.score))@ranges
        en<-suppressWarnings(matchPWM(get(paste0('pwm',n)),x,min.score))@ranges
        # find regions spanning the first exon and the last exon 
        ov<-findOverlapPairs(e1,en,maxgap=span)
        ov<-ov[start(ov@second)-end(ov@first)>=sum(exons[-c(1,n)])]
        if (length(ov)==0) return(IRanges())
        e1<-ov@first; en<-ov@second
        intervals<-pgap(e1,en)
        mcols(intervals)$id<-seq(along.with = intervals)
        # screen the intervals containing all the other exons
        for(i in seq(2,n-1)) {
            mi<-suppressWarnings(matchPWM(get(paste0('pwm',i)),x,min.score))@ranges
            ov<-findOverlapPairs(mi,intervals,type='within')
            if (length(ov)==0) return(IRanges())
            intervals<-unique(ov@second)
            assign(paste0('m',i),unique(ov@first))
        }
        e1<-e1[mcols(intervals)$id]
        en<-en[mcols(intervals)$id]
        # check the order of the exons and introns
        me<-lapply(seq(along.with = intervals), function(i) {
            nt<-x[e1[i]] # fist exon
            block<-data.frame(start(e1[i]))
            ir<-intervals[i]
            for (j in seq(2,n-1)) {
                ov<-findOverlapPairs(get(paste0('m',j)),ir,type='within')
                if (length(ov)==0) return(IRanges())
                l<-grepl(plus_intron1,Views(x,flank(
                    ov@first,start(ov@first)-start(ov@second))),perl=TRUE)
                if (!any(l)) l<-grepl(plus_intron2,Views(x,flank(
                    ov@first,start(ov@first)-start(ov@second))),perl=TRUE)
                if (!any(l)) return(IRanges())
                ov<-ov[l]
                nt<-paste0(nt,Views(x,ov@first))
                block<-data.frame(data.frame(lapply(block, rep, length(ov))),
                                  rep(start(ov@first),each=nrow(block)))
                ir<-pgap(ov@first,rep(en[i],length(ov)))
            }
            r<-grepl(plus_intron1,Views(x,ir),perl=TRUE)
            if (!any(r)) r<-grepl(plus_intron2,Views(x,ir),perl=TRUE)
            if (!any(r)) return(IRanges())
            nt<-DNAStringSet(paste0(nt[r],Views(x,en[i])))
            block<-as.matrix(data.frame(block[r,,drop=FALSE],start(en)[i]))
            dimnames(block)<-NULL
            ir<-IRanges(block[,focus],width=exons[focus],
                        NT=nt,AA=Biostrings::translate(nt,if.fuzzy.codon='solve'))
            mcols(ir)$region<-IRanges(start(e1[i]),end(en[i]))
            mcols(ir)$block.starts<-t(apply(block,1,function(xx) {
                head(sort(c(xx,xx+exons)),-1)
            }))
            mcols(ir)$block.sizes<-t(apply(block,1,function(xx) {
                diff(sort(c(xx,xx+exons)))
            }))
            return(ir)
        })
        me<-do.call(c,me[elementNROWS(me)>0])
        me
    },mc.cores = cores)
    plus<-unlist(IRangesList(plus[elementNROWS(plus)>0]))
    if (length(plus)==0) plus<-GRanges() else plus<-GRanges(names(plus),plus,'+')
    
    minus<-mclapply(genome, function(x) {
        # find matches for the firest and the last exon
        e1<-suppressWarnings(matchPWM(reverseComplement(get('pwm1')),x,min.score))@ranges
        en<-suppressWarnings(matchPWM(reverseComplement(
            get(paste0('pwm',n))),x,min.score))@ranges
        # find regions spanning the first exon and the last exon 
        ov<-findOverlapPairs(e1,en,maxgap=span)
        ov<-ov[start(ov@first)-end(ov@second)>=sum(exons[-c(1,n)])]
        if (length(ov)==0) return(IRanges())
        e1<-ov@first; en<-ov@second
        intervals<-pgap(e1,en)
        mcols(intervals)$id<-seq(along.with = intervals)
        # screen the intervals containing all the other exons
        for(i in seq(2,n-1)) {
            mi<-suppressWarnings(matchPWM(reverseComplement(
                get(paste0('pwm',i))),x,min.score))@ranges
            ov<-findOverlapPairs(mi,intervals,type='within')
            if (length(ov)==0) return(IRanges())
            intervals<-unique(ov@second)
            assign(paste0('m',i),unique(ov@first))
        }
        e1<-e1[mcols(intervals)$id]
        en<-en[mcols(intervals)$id]
        # check the order of the exons and introns
        me<-lapply(seq(along.with = intervals), function(i) {
            nt<-x[e1[i]] # fist exon
            block<-data.frame(start(e1[i]))
            ir<-intervals[i]
            for (j in seq(2,n-1)) {
                ov<-findOverlapPairs(get(paste0('m',j)),ir,type='within')
                if (length(ov)==0) return(IRanges())
                r<-grepl(minus_intron1,Views(x,flank(
                    ov@first,end(ov@second)-end(ov@first),start = FALSE)),
                    perl=TRUE)
                if  (!any(r)) r<-grepl(minus_intron2,Views(x,flank(
                    ov@first,end(ov@second)-end(ov@first),start = FALSE)),
                    perl=TRUE)
                if (!any(r)) return(IRanges())
                ov<-ov[r]
                nt<-paste0(Views(x,ov@first),nt)
                block<-data.frame(data.frame(lapply(block, rep, length(ov))),
                                  rep(start(ov@first),each=nrow(block)))
                ir<-pgap(ov@first,rep(en[i],length(ov)))
            }
            l<-grepl(minus_intron1,Views(x,ir),perl=TRUE)
            if (!any(l)) l<-grepl(minus_intron2,Views(x,ir),perl=TRUE)
            if (!any(l)) return(IRanges())
            nt<-reverseComplement(DNAStringSet(paste0(Views(x,en[i]),nt[l])))
            block<-as.matrix(data.frame(block[l,,drop=FALSE],start(en)[i]))
            dimnames(block)<-NULL
            ir<-IRanges(block[,focus],width=exons[focus],NT=nt,
                        AA=Biostrings::translate(nt,if.fuzzy.codon='solve'))
            mcols(ir)$region<-IRanges(start(en[i]),end(e1[i]))
            mcols(ir)$block.starts<-t(apply(block,1,function(xx) {
                sort(c(xx,xx+exons),decreasing = TRUE)[-1]
            }))
            mcols(ir)$block.sizes<-t(apply(block,1,function(xx) {
                0L-diff(sort(c(xx,xx+exons),decreasing = TRUE))
            }))
            return(ir)
        })
        me<-do.call(c,me[elementNROWS(me)>0])
        me
    },mc.cores = cores)
    minus<-unlist(IRangesList(minus[elementNROWS(minus)>0]))
    if (length(minus)==0) minus<-GRanges() else minus<-GRanges(names(minus),minus,'-')
    
    res<-append(plus,minus)
    res<-res[grep("\\*",res$AA,invert = TRUE)]
    if(length(res)>0) {
        res$score<-round(vapply(res$NT, function(x) PWMscoreStartingAt(pwm,x), 1.0),4)
        cluster<-reduce(GRanges(seqnames(res),res$region,strand(res)),with.revmap=TRUE)$revmap
        if(any(elementNROWS(cluster)>1)) {
            res<-res[unlist(cluster[lapply(extractList(res$score,cluster), which.max)])]
        }
        res$isMicroExon<-res$block.sizes[,2*(focus-1)]!=0 & res$block.sizes[,2*focus]!=0
    }
    cat ("---Number of microexon-tags found:\t", length(res),'\n')
    unname(res)
}

###
mapPWM.0<-function(cons,exons,genome,min.score="80%",include.intronLoss=TRUE,
                 span=20000, min.intron=20, max.intron=10000, 
                 prior.params=letterFrequency(genome[[1]],DNA_BASES,as.prob = TRUE),
                 cores=1){
    pwm<-PWM(cons,prior.params=prior.params)
    pwm1<-PWM(cons[,seq.int(exons[1]),drop=FALSE],prior.params=prior.params)
    pwm2<-PWM(cons[,seq.int(exons[2])+exons[1],drop=FALSE],prior.params=prior.params)
    pwm3<-PWM(cons[,seq.int(exons[3])+sum(exons[1:2]),drop=FALSE],prior.params=prior.params)
    plus_intron1<-sprintf('^(GT[AG].{%d,%d}[CT]AG)%s$',min.intron,max.intron,
                          ifelse(include.intronLoss,"?",""))
    plus_intron2<-sprintf('^(G[TC].{%d,%d}AG)%s$',min.intron,max.intron,
                          ifelse(include.intronLoss,"?",""))
    minus_intron1<-sprintf('^(CT[AG].{%d,%d}[CT]AC)%s$',min.intron,max.intron,
                           ifelse(include.intronLoss,"?",""))
    minus_intron2<-sprintf('^(CT.{%d,%d}[GA]C)%s$',min.intron,max.intron,
                           ifelse(include.intronLoss,"?",""))
    
    plus<-mclapply(genome, function(x) {
        e1<-suppressWarnings(matchPWM(pwm1,x,min.score))
        e3<-suppressWarnings(matchPWM(pwm3,x,min.score))
        ov<-findOverlapPairs(e1,e3,maxgap=span)
        ov<-ov[start(ov@second)>=end(ov@first)+exons[2]]
        if (length(ov)==0) return(IRanges())
        e1<-ov@first;e3<-ov@second
        intervals<-pgap(e1,e3)
        me<-lapply(seq(along.with = intervals), function(i) {
            xx<-Views(x,intervals[i])
            m<-suppressWarnings(matchPWM(pwm2,xx,min.score))
            if(length(m)==0) return(IRanges())
            l<-grepl(plus_intron1,flank(m,start(m)-start(xx)),perl=TRUE)
            if (!any(l)) l<-grepl(plus_intron2,flank(m,start(m)-start(xx)),perl=TRUE)
            r<-grepl(plus_intron1,flank(m,end(xx)-end(m),start=FALSE),perl=TRUE)
            if (!any(r)) r<-grepl(plus_intron2,flank(m,end(xx)-end(m),start=FALSE),perl=TRUE)
            m<-m[l&r]
            if(length(m)==0) return(IRanges())
            e2<-m@ranges
            mcols(e2)$NT<-DNAStringSet(paste0(e1[i],m,e3[i]))
            mcols(e2)$AA<-Biostrings::translate(mcols(e2)$NT,if.fuzzy.codon='solve')
            mcols(e2)$region<-IRanges(start(e1[i]),end(e3[i]))
            mcols(e2)$block.starts<-as.matrix(data.frame(
                start(e1)[i],end(e1)[i]+1,
                start(e2),end(e2)+1,start(e3)[i],fix.empty.names=FALSE))
            mcols(e2)$block.sizes<-as.matrix(data.frame(
                exons[1],distance(e1[i],e2),
                exons[2],distance(e2,e3[i]),exons[3],fix.empty.names=FALSE))
            return(e2)
        })
        me<-do.call(c,me[elementNROWS(me)>0])
        me
    },mc.cores = cores)
    plus<-unlist(IRangesList(plus[elementNROWS(plus)>0]))
    if (length(plus)==0) plus<-GRanges() else plus<-GRanges(names(plus),plus,'+')
    
    minus<-mclapply(genome, function(x) {
        e1<-suppressWarnings(matchPWM(reverseComplement(pwm1),x,min.score))
        e3<-suppressWarnings(matchPWM(reverseComplement(pwm3),x,min.score))
        ov<-findOverlapPairs(e1,e3,maxgap=span)
        ov<-ov[start(ov@first)>=end(ov@second)+exons[2]]
        if (length(ov)==0) return(IRanges())
        e1<-ov@first;e3<-ov@second
        intervals<-pgap(e1,e3)
        me<-lapply(seq(along.with = intervals), function(i) {
            xx<-Views(x,intervals[i])
            m<-suppressWarnings(matchPWM(reverseComplement(pwm2),xx,min.score))
            if(length(m)==0) return(IRanges())
            l<-grepl(minus_intron1,flank(m,start(m)-start(xx)),perl=TRUE)
            if (!any(l)) l<-grepl(minus_intron2,flank(m,start(m)-start(xx)),perl=TRUE)
            r<-grepl(minus_intron1,flank(m,end(xx)-end(m),start=FALSE),perl=TRUE)
            if (!any(r)) r<-grepl(minus_intron2,flank(m,end(xx)-end(m),start=FALSE),perl=TRUE)
            m<-m[l&r]
            if(length(m)==0) return(IRanges())
            e2<-m@ranges
            mcols(e2)$NT<-reverseComplement(DNAStringSet(paste0(e3[i],m,e1[i])))
            mcols(e2)$AA<-Biostrings::translate(mcols(e2)$NT,if.fuzzy.codon='solve')
            mcols(e2)$region<-IRanges(start(e3[i]),end(e1[i]))
            mcols(e2)$block.starts<-as.matrix(data.frame(
                start(e1)[i],end(e2)+1,start(e2),end(e3)[i]+1,start(e3)[i],
                fix.empty.names=FALSE))
            mcols(e2)$block.sizes<-as.matrix(data.frame(
                exons[1],distance(e1[i],e2),exons[2],distance(e2,e3[i]),
                exons[3],fix.empty.names=FALSE))
            return(e2)
        })
        me<-do.call(c,me[elementNROWS(me)>0])
        me
    },mc.cores = cores)
    minus<-unlist(IRangesList(minus[elementNROWS(minus)>0]))
    if (length(minus)==0) minus<-GRanges() else minus<-GRanges(names(minus),minus,'-')
    
    res<-append(plus,minus)
    res<-res[grep("\\*",res$AA,invert = TRUE)]
    if(length(res)>0) {
        res$score<-round(vapply(res$NT, function(x) PWMscoreStartingAt(pwm,x), 1.0),4)
        cluster<-reduce(GRanges(seqnames(res),res$region,strand(res)),with.revmap=TRUE)$revmap
        if(any(elementNROWS(cluster)>1)) {
            res<-res[unlist(cluster[lapply(extractList(res$score,cluster), which.max)])]
        }
        res$isMicroExon<-res$block.sizes[,2]!=0 & res$block.sizes[,4]!=0
    }
    cat ("---Number of microexon-tags found:\t", length(res),'\n')
    unname(res)
}

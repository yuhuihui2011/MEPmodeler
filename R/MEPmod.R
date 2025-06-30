#' Predict Conserved MicroExons in a Plant Genome
#' 
#' Identifies microexon-tags in plant genomic sequences using gapped 
#' Position Weight Matrices (PWMs). Scans both DNA strands to detect conserved 
#' microexon clusters based on sequence patterns.
#' 
#' @param genome Path(s) to FASTA file(s) or a \code{DNAStringSet} object 
#'              containing genomic sequences.
#' @param min.score Minimum match score threshold for exon blocks 
#'                 (as percentage string, e.g., \code{"80\%"}).
#'                 Passed to \code{\link[Biostrings]{matchPWM}}.
#' @param clusters Integer vector specifying which microexon-tag clusters to search 
#'                (default: \code{NULL} searches all clusters).
#' @param high_complexity Logical. If \code{TRUE} (default), only searches 
#'                       clusters with sequence complexity ≥ 0.7 for both 
#'                       flanking regions.
#' @param low_copy Logical. If \code{TRUE}, restricts search to low-copy clusters 
#'                 (1-3 copies in >80\% of 184 land plants).
#' @param include.intronLoss Logical. If \code{TRUE}, includes microexons with 
#'                           any flanking intron loss.
#' @param span Maximum genomic span for microexon-tag detection (default: 20,000 bp).
#' @param min.intron Minimum allowed intron size (bp).
#' @param max.intron Maximum allowed intron size (bp).
#' @param cores Number of CPU cores to use (passed to \code{\link[parallel]{mclapply}}).
#' @param quietly Logical. If \code{TRUE}, suppresses all progress messages. 
#'               Default: \code{FALSE} (show detailed progress messages).
#' 
#' @return A \code{GRanges} object with 8 metadata columns:
#' \describe{
#'   \item{NT}{DNA sequence of the microexon-tag}
#'   \item{AA}{Peptide sequence encoded by the microexon-tag}
#'   \item{region}{Genomic span of the microexon-tag}
#'   \item{block.starts}{Genomic coordinates of exon/intron boundaries}
#'   \item{block.sizes}{Lengths of exons and introns}
#'   \item{score}{Match confidence score (higher = more confident)}
#'   \item{isMicroExon}{\code{TRUE} if no flanking intron loss, \code{FALSE} otherwise}
#'   \item{cluster}{Cluster ID of the detected microexon-tag}
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
#' data <- DNAStringSet(BSgenome.Athaliana.TAIR.TAIR9$Chr1[1:1e6])
#' names(data) <- "Chr1"
#' 
#' res <- MEPmod(genome=data, cores=1, clusters=c(1:50))
#'
#' res
#' # GRanges object with 2 ranges and 8 metadata columns:
#' #        seqnames       ranges strand |                      NT
#' #          <Rle>     <IRanges>  <Rle> |          <DNAStringSet>
#' #   [1]     Chr1 456757-456761      + | TTTGATGCTA...ACTGTCTGAC
#' #   [2]     Chr1 601547-601560      - | AAAGTTGAAT...ACATTCAGAT
#' #                            AA        region             block.starts   block.sizes
#' #                 <AAStringSet>     <IRanges>            <IntegerList> <IntegerList>
#' #   [1] FDARTAWSQC...AFGAVESLSD 456595-456893 456595,456647,456757,...  52,110,5,...
#' #   [2] KVEFKDNEWK...RLDCLLKHSD 601327-601700 601653,601561,601547,...  48,92,14,...
#' #            score isMicroExon   cluster
#' #        <numeric>   <logical> <integer>
#' #   [1]    0.9817        TRUE         8
#' #   [2]    0.9242        TRUE        37
#' #   -------
#' #   seqinfo: 1 sequence from an unspecified genome; no seqlengths
#' 

MEPmod <- function(genome, min.score = '80%', clusters = NULL,
                   high_complexity = TRUE, low_copy = FALSE,
                   include.intronLoss = TRUE, span = 20000, 
                   min.intron = 20, max.intron = 10000, cores = 5,
                   quietly = FALSE) {
    
    MEPdata <- MEPmodeler::MEPdata
    
    # Cluster selection logic
    if (is.null(clusters)) {
        clusters <- MEPdata$cluster$cluster
    }
    if (high_complexity) {
        sel <- (MEPdata$cluster$L_complex >= 0.7) & 
            (MEPdata$cluster$R_complex >= 0.7)
        clusters <- clusters[clusters %in% MEPdata$cluster$cluster[sel]]
    }
    if (low_copy) {
        clusters <- clusters[clusters %in% MEPdata$cluster$cluster[MEPdata$cluster$low_copy]]
    }
    
    if (length(clusters) == 0) {
        stop("No microexon cluster will be searched, please check the parameters!")
    } else if (!quietly){
        cat(paste(t0 <- format(Sys.time()),
                  sprintf('..... start to search %d microexon cluster(s)\n', length(clusters))))
    }
    
    suppressPackageStartupMessages({
        library(Biostrings)
        library(GenomicRanges)
        library(parallel)
    })
    
    if (!quietly) cat(paste(format(Sys.time()), '..... loading genome\n'))
    if (!is(genome, 'DNAStringSet')) {
        if (!is.character(genome)) {
            stop("'genome' must be the path(s) to the FASTA file(s) or a 'DNAStringSet' class")
        }
        genome <- readDNAStringSet(genome)
    }
    
    names(genome) <- sub("^(\\S+)\\s+.*", "\\1", names(genome))
    prior.params <- letterFrequency(genome[order(-width(genome))][[1]],
                                    DNA_BASES, as.prob = TRUE)
    if (!quietly) cat(paste(format(Sys.time()), '..... finished loading genome\n'))
    
    res <- GRanges()
    for (i in as.integer(clusters)) {
        if (!quietly) cat(sprintf('%s ..... Cluster %d (size: %d; phase: %d; motif: %s)\n',
                                 format(Sys.time()), i, MEPdata$cluster$size[i],
                                 MEPdata$cluster$phase[i],
                                 MEPdata$cluster$motif[i]))
        
        res_i <- mapPWM(
            cons = MEPdata$matrix[[i]],
            exons = MEPdata$blocks[[i]],
            genome = genome,
            focus = MEPdata$cluster$me_order[i],
            min.score = min.score,
            include.intronLoss = include.intronLoss,
            span = span,
            min.intron = min.intron,
            max.intron = max.intron,
            prior.params = prior.params,
            cores = cores,
            check.params = FALSE,
            quietly = quietly
        )
        
        if (length(res_i) == 0) next
        
        res_i$block.starts <- do.call(c, apply(res_i$block.starts, 1, IntegerList))
        res_i$block.sizes <- do.call(c, apply(res_i$block.sizes, 1, IntegerList))
        res_i$cluster <- i
        res <- suppressWarnings(append(res, res_i))
        invisible(gc())
    }
    
    if (!quietly) {
        cat(paste(t1 <- format(Sys.time()), '..... finished successfully\n'))
        cat(sprintf('### total microexon-tags found:\t%d\n', length(res)))
        cat(sprintf('### total time used:\t%s mins\n', 
                    round(difftime(t1, t0, units = 'mins'), 2)))
    }
    res
}

####
mapPWM <- function(cons, exons, genome, focus = 2, min.score = '80%', 
                   include.intronLoss = TRUE, span = 20000, min.intron = 20, 
                   max.intron = 10000, prior.params = letterFrequency(genome[[1]], 
                   DNA_BASES, as.prob = TRUE), cores = 1, check.params = TRUE,
                   quietly = FALSE) {
    
    n <- length(exons)
    if (n == 3) {
        return(mapPWM.0(cons, exons, genome, min.score, include.intronLoss, span, 
                        min.intron, max.intron, prior.params, cores, quietly))
    }
    
    if (check.params) {
        if (n < 3) stop("'exons' must be a numeric vector of length >= 3")
        if (sum(exons) != ncol(cons))
            stop("sum(exons)==ncol(cons) is not TRUE")
        if (!is(genome, 'DNAStringSet')) 
            stop("'genome' must be a DNAStringSet object")
        if (length(focus) != 1 & focus[1] <= 1 & focus[1] >= n)
            stop("'focus' must be a numeric in (0,length(exons))!")
    }
    
    plus_intron1 <- sprintf('^(GT[AG].{%d,%d}[CT]AG)%s$', min.intron, max.intron,
                            ifelse(include.intronLoss, "?", ""))
    plus_intron2 <- sprintf('^(G[TC].{%d,%d}AG)%s$', min.intron, max.intron,
                            ifelse(include.intronLoss, "?", ""))
    minus_intron1 <- sprintf('^(CT[AG].{%d,%d}[CT]AC)%s$', min.intron, max.intron,
                             ifelse(include.intronLoss, "?", ""))
    minus_intron2 <- sprintf('^(CT.{%d,%d}[GA]C)%s$', min.intron, max.intron,
                             ifelse(include.intronLoss, "?", ""))
    
    pwm <- PWM(cons, prior.params = prior.params)
    for (i in seq(n)) {
        index <- seq.int(exons[i]) + cumsum(c(0, exons))[i]
        assign(paste0('pwm', i), PWM(cons[, index, drop = FALSE], 
                                     prior.params = prior.params))
    }
    
    plus <- mclapply(genome, function(x) {
        e1 <- suppressWarnings(matchPWM(get('pwm1'), x, min.score))@ranges
        en <- suppressWarnings(matchPWM(get(paste0('pwm', n)), x, min.score))@ranges
        
        ov <- findOverlapPairs(e1, en, maxgap = span)
        ov <- ov[start(ov@second) - end(ov@first) >= sum(exons[-c(1, n)])]
        if (length(ov) == 0) return(IRanges())
        
        e1 <- ov@first
        en <- ov@second
        intervals <- pgap(e1, en)
        mcols(intervals)$id <- seq(along.with = intervals)
        
        for (i in seq(2, n - 1)) {
            mi <- suppressWarnings(matchPWM(get(paste0('pwm', i)), x, min.score))@ranges
            ov <- findOverlapPairs(mi, intervals, type = 'within')
            if (length(ov) == 0) return(IRanges())
            intervals <- unique(ov@second)
            assign(paste0('m', i), unique(ov@first))
        }
        
        e1 <- e1[mcols(intervals)$id]
        en <- en[mcols(intervals)$id]
        
        me <- lapply(seq(along.with = intervals), function(i) {
            nt <- x[e1[i]] # first exon
            block <- data.frame(start(e1[i]))
            ir <- intervals[i]
            
            for (j in seq(2, n - 1)) {
                ov <- findOverlapPairs(get(paste0('m', j)), ir, type = 'within')
                if (length(ov) == 0) return(IRanges())
                l <- grepl(plus_intron1, Views(x, flank(
                    ov@first, start(ov@first) - start(ov@second))), perl = TRUE)
                if (!any(l)) l <- grepl(plus_intron2, Views(x, flank(
                    ov@first, start(ov@first) - start(ov@second))), perl = TRUE)
                if (!any(l)) return(IRanges())
                
                ov <- ov[l]
                nt <- paste0(nt, Views(x, ov@first))
                block <- data.frame(data.frame(lapply(block, rep, length(ov))),
                                    rep(start(ov@first), each = nrow(block)))
                ir <- pgap(ov@first, rep(en[i], length(ov)))
            }
            
            r <- grepl(plus_intron1, Views(x, ir), perl = TRUE)
            if (!any(r)) r <- grepl(plus_intron2, Views(x, ir), perl = TRUE)
            if (!any(r)) return(IRanges())
            
            nt <- DNAStringSet(paste0(nt[r], Views(x, en[i])))
            block <- as.matrix(data.frame(block[r, , drop = FALSE], start(en)[i]))
            dimnames(block) <- NULL
            
            ir <- IRanges(block[, focus], width = exons[focus],
                          NT = nt, AA = Biostrings::translate(nt, if.fuzzy.codon = 'solve'))
            mcols(ir)$region <- IRanges(start(e1[i]), end(en[i]))
            mcols(ir)$block.starts <- t(apply(block, 1, function(xx) {
                head(sort(c(xx, xx + exons)), -1)
            }))
            mcols(ir)$block.sizes <- t(apply(block, 1, function(xx) {
                diff(sort(c(xx, xx + exons)))
            }))
            return(ir)
        })
            do.call(c, me[elementNROWS(me) > 0])
    }, mc.cores = cores)
        
        plus <- unlist(IRangesList(plus[elementNROWS(plus) > 0]))
        if (length(plus) == 0) plus <- GRanges() else 
            plus <- GRanges(names(plus), plus, '+')
        
        minus <- mclapply(genome, function(x) {
            e1 <- suppressWarnings(matchPWM(reverseComplement(get('pwm1')), x, min.score))@ranges
            en <- suppressWarnings(matchPWM(reverseComplement(
                get(paste0('pwm', n))), x, min.score))@ranges
            
            ov <- findOverlapPairs(e1, en, maxgap = span)
            ov <- ov[start(ov@first) - end(ov@second) >= sum(exons[-c(1, n)])]
            if (length(ov) == 0) return(IRanges())
            
            e1 <- ov@first
            en <- ov@second
            intervals <- pgap(e1, en)
            mcols(intervals)$id <- seq(along.with = intervals)
            
            for (i in seq(2, n - 1)) {
                mi <- suppressWarnings(matchPWM(reverseComplement(
                    get(paste0('pwm', i))), x, min.score))@ranges
                ov <- findOverlapPairs(mi, intervals, type = 'within')
                if (length(ov) == 0) return(IRanges())
                intervals <- unique(ov@second)
                assign(paste0('m', i), unique(ov@first))
            }
            
            e1 <- e1[mcols(intervals)$id]
            en <- en[mcols(intervals)$id]
            
            me <- lapply(seq(along.with = intervals), function(i) {
                nt <- x[e1[i]] # first exon
                block <- data.frame(start(e1[i]))
                ir <- intervals[i]
                
                for (j in seq(2, n - 1)) {
                    ov <- findOverlapPairs(get(paste0('m', j)), ir, type = 'within')
                    if (length(ov) == 0) return(IRanges())
                    r <- grepl(minus_intron1, Views(x, flank(
                        ov@first, end(ov@second) - end(ov@first), start = FALSE)),
                        perl = TRUE)
                    if (!any(r)) r <- grepl(minus_intron2, Views(x, flank(
                        ov@first, end(ov@second) - end(ov@first), start = FALSE)),
                        perl = TRUE)
                    if (!any(r)) return(IRanges())
                    
                    ov <- ov[r]
                    nt <- paste0(Views(x, ov@first), nt)
                    block <- data.frame(data.frame(lapply(block, rep, length(ov))),
                                        rep(start(ov@first), each = nrow(block)))
                    ir <- pgap(ov@first, rep(en[i], length(ov)))
                }
                
                l <- grepl(minus_intron1, Views(x, ir), perl = TRUE)
                if (!any(l)) l <- grepl(minus_intron2, Views(x, ir), perl = TRUE)
                if (!any(l)) return(IRanges())
                
                nt <- reverseComplement(DNAStringSet(paste0(Views(x, en[i]), nt[l])))
                block <- as.matrix(data.frame(block[l, , drop = FALSE], start(en)[i]))
                dimnames(block) <- NULL
                
                ir <- IRanges(block[, focus], width = exons[focus], NT = nt,
                              AA = Biostrings::translate(nt, if.fuzzy.codon = 'solve'))
                mcols(ir)$region <- IRanges(start(en[i]), end(e1[i]))
                mcols(ir)$block.starts <- t(apply(block, 1, function(xx) {
                    sort(c(xx, xx + exons), decreasing = TRUE)[-1]
                }))
                mcols(ir)$block.sizes <- t(apply(block, 1, function(xx) {
                    0L - diff(sort(c(xx, xx + exons), decreasing = TRUE))
                }))
                return(ir)
            })
            do.call(c, me[elementNROWS(me) > 0])
        }, mc.cores = cores)
        
        minus <- unlist(IRangesList(minus[elementNROWS(minus) > 0]))
        if (length(minus) == 0) minus <- GRanges() else 
            minus <- GRanges(names(minus), minus, '-')
        
        res <- suppressWarnings(append(plus, minus))
        res <- res[grep("\\*", res$AA, invert = TRUE)]
        
        if (length(res) > 0) {
            res$score <- round(vapply(res$NT, function(x) PWMscoreStartingAt(pwm, x), 1.0), 4)
            cluster <- reduce(GRanges(seqnames(res), res$region, strand(res)), with.revmap = TRUE)$revmap
            if (any(elementNROWS(cluster) > 1)) {
                res <- res[unlist(cluster[lapply(extractList(res$score, cluster), which.max)])]
            }
            res$isMicroExon <- res$block.sizes[, 2 * (focus - 1)] != 0 & 
                res$block.sizes[, 2 * focus] != 0
        }
        
        if (!quietly) cat("---Number of microexon-tags found:\t", length(res), '\n')
        unname(res)
}

mapPWM.0 <- function(cons, exons, genome, min.score = "80%", 
                     include.intronLoss = TRUE, span = 20000, 
                     min.intron = 20, max.intron = 10000, 
                     prior.params = letterFrequency(genome[[1]], 
                     DNA_BASES, as.prob = TRUE), cores = 1, quietly = FALSE) {
    
    pwm <- PWM(cons, prior.params = prior.params)
    pwm1 <- PWM(cons[, seq.int(exons[1]), drop = FALSE], prior.params = prior.params)
    pwm2 <- PWM(cons[, seq.int(exons[2]) + exons[1], drop = FALSE], prior.params = prior.params)
    pwm3 <- PWM(cons[, seq.int(exons[3]) + sum(exons[1:2]), drop = FALSE], prior.params = prior.params)
    
    plus_intron1 <- sprintf('^(GT[AG].{%d,%d}[CT]AG)%s$', min.intron, max.intron,
                            ifelse(include.intronLoss, "?", ""))
    plus_intron2 <- sprintf('^(G[TC].{%d,%d}AG)%s$', min.intron, max.intron,
                            ifelse(include.intronLoss, "?", ""))
    minus_intron1 <- sprintf('^(CT[AG].{%d,%d}[CT]AC)%s$', min.intron, max.intron,
                             ifelse(include.intronLoss, "?", ""))
    minus_intron2 <- sprintf('^(CT.{%d,%d}[GA]C)%s$', min.intron, max.intron,
                             ifelse(include.intronLoss, "?", ""))
    
    plus <- mclapply(genome, function(x) {
        e1 <- suppressWarnings(matchPWM(pwm1, x, min.score))
        e3 <- suppressWarnings(matchPWM(pwm3, x, min.score))
        
        ov <- findOverlapPairs(e1, e3, maxgap = span)
        ov <- ov[start(ov@second) >= end(ov@first) + exons[2]]
        if (length(ov) == 0) return(IRanges())
        
        e1 <- ov@first
        e3 <- ov@second
        intervals <- pgap(e1, e3)
        
        me <- lapply(seq_along(intervals), function(i) {
            xx <- Views(x, intervals[i])
            m <- suppressWarnings(matchPWM(pwm2, xx, min.score))
            if (length(m) == 0) return(IRanges())
            
            l <- grepl(plus_intron1, flank(m, start(m) - start(xx)), perl = TRUE)
            if (!any(l)) l <- grepl(plus_intron2, flank(m, start(m) - start(xx)), perl = TRUE)
            r <- grepl(plus_intron1, flank(m, end(xx) - end(m), start = FALSE), perl = TRUE)
            if (!any(r)) r <- grepl(plus_intron2, flank(m, end(xx) - end(m), start = FALSE), perl = TRUE)
            
            m <- m[l & r]
            if (length(m) == 0) return(IRanges())
            
            e2 <- m@ranges
            mcols(e2)$NT <- DNAStringSet(paste0(e1[i], m, e3[i]))
            mcols(e2)$AA <- Biostrings::translate(mcols(e2)$NT, if.fuzzy.codon = 'solve')
            mcols(e2)$region <- IRanges(start(e1[i]), end(e3[i]))
            mcols(e2)$block.starts <- as.matrix(data.frame(
                start(e1)[i], end(e1)[i] + 1,
                start(e2), end(e2) + 1, start(e3)[i], fix.empty.names = FALSE))
            mcols(e2)$block.sizes <- as.matrix(data.frame(
                exons[1], distance(e1[i], e2),
                exons[2], distance(e2, e3[i]), exons[3], fix.empty.names = FALSE))
            return(e2)
        })
        do.call(c, me[elementNROWS(me) > 0])
    }, mc.cores = cores)
    
    plus <- unlist(IRangesList(plus[elementNROWS(plus) > 0]))
    if (length(plus) == 0) plus <- GRanges() else 
        plus <- GRanges(names(plus), plus, '+')
    
    minus <- mclapply(genome, function(x) {
        e1 <- suppressWarnings(matchPWM(reverseComplement(pwm1), x, min.score))
        e3 <- suppressWarnings(matchPWM(reverseComplement(pwm3), x, min.score))
        
        ov <- findOverlapPairs(e1, e3, maxgap = span)
        ov <- ov[start(ov@first) >= end(ov@second) + exons[2]]
        if (length(ov) == 0) return(IRanges())
        
        e1 <- ov@first
        e3 <- ov@second
        intervals <- pgap(e1, e3)
        
        me <- lapply(seq(along.with = intervals), function(i) {
            xx <- Views(x, intervals[i])
            m <- suppressWarnings(matchPWM(reverseComplement(pwm2), xx, min.score))
            if (length(m) == 0) return(IRanges())
            
            l <- grepl(minus_intron1, flank(m, start(m) - start(xx)), perl = TRUE)
            if (!any(l)) l <- grepl(minus_intron2, flank(m, start(m) - start(xx)), perl = TRUE)
            r <- grepl(minus_intron1, flank(m, end(xx) - end(m), start = FALSE), perl = TRUE)
            if (!any(r)) r <- grepl(minus_intron2, flank(m, end(xx) - end(m), start = FALSE), perl = TRUE)
            
            m <- m[l & r]
            if (length(m) == 0) return(IRanges())
            
            e2 <- m@ranges
            mcols(e2)$NT <- reverseComplement(DNAStringSet(paste0(e3[i], m, e1[i])))
            mcols(e2)$AA <- Biostrings::translate(mcols(e2)$NT, if.fuzzy.codon = 'solve')
            mcols(e2)$region <- IRanges(start(e3[i]), end(e1[i]))
            mcols(e2)$block.starts <- as.matrix(data.frame(
                start(e1)[i], end(e2) + 1, start(e2), end(e3)[i] + 1, start(e3)[i],
                fix.empty.names = FALSE))
            mcols(e2)$block.sizes <- as.matrix(data.frame(
                exons[1], distance(e1[i], e2), exons[2], distance(e2, e3[i]),
                exons[3], fix.empty.names = FALSE))
            return(e2)
        })
        do.call(c, me[elementNROWS(me) > 0])
    }, mc.cores = cores)
    
    minus <- unlist(IRangesList(minus[elementNROWS(minus) > 0]))
    if (length(minus) == 0) minus <- GRanges() else 
        minus <- GRanges(names(minus), minus, '-')
    
    res <- suppressWarnings(append(plus, minus))
    res <- res[grep("\\*", res$AA, invert = TRUE)]
    
    if (length(res) > 0) {
        res$score <- round(vapply(res$NT, function(x) PWMscoreStartingAt(pwm, x), 1.0), 4)
        cluster <- reduce(GRanges(seqnames(res), res$region, strand(res)), with.revmap = TRUE)$revmap
        if (any(elementNROWS(cluster) > 1)) {
            res <- res[unlist(cluster[lapply(extractList(res$score, cluster), which.max)])]
        }
        res$isMicroExon <- res$block.sizes[, 2] != 0 & res$block.sizes[, 4] != 0
    }
    
    if (!quietly) cat("---Number of microexon-tags found:\t", length(res), '\n')
    unname(res)
}

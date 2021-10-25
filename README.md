# MEPmodeler: Predict <i>M</i>icro<i>E</i>xons in <i>P</i>lant Genomes

MEPmodeler searches microexon-tags in 45 conserved microexon clusters in plant genomes using gapped Position Weight Matrix (PWM). The only one function, ***MEPmod***, uses the data ***MEPdata*** from [MEPsuite](https://github.com/yuhuihui2011/MEPsuite) to model microexons in plant genomes. The only input file is plant genomic sequences or a plant genome. The sequences on both plus and minus strands will be scanned. 
<br>

## Table of contents
- [Requirements](#requirements)
- [Installation](#installation)
- [Example](#example)
- [Parameters](#parameters)
- [Output](#output)
- [Downstream analysis](#downstream-analysis)
- [Citation](#citation)

## Requirements
MEPmodeler requires the following R (>= 4.0) packages and the dependents:
+ [Biostrings](https://bioconductor.org/packages/Biostrings) (>= 2.58.0)
+ [GenomicRanges](https://bioconductor.org/packages/GenomicRanges) (>= 1.42.0)
+ [BSgenome](https://bioconductor.org/packages/BSgenome) (>= 1.58.0)

## Installation
Start R (>= 4.0) and run:
```
devtools::install_github('yuhuihui2011/MEPmodeler')
```

## Example
The following example is to predict microexons from Arabidopsis Chr1:1-1000000 sequence: 
```{r}
library(MEPmodeler)
suppressPackageStartupMessages(library(BSgenome))
if (!requireNamespace("BSgenome.Athaliana.TAIR.TAIR9", quietly = TRUE))
    BiocManager::install("BSgenome.Athaliana.TAIR.TAIR9")
library("BSgenome.Athaliana.TAIR.TAIR9")
BSgenome.Athaliana.TAIR.TAIR9
# Arabidopsis genome:
# # organism: Arabidopsis thaliana (Arabidopsis)
# # genome: TAIR9
# # provider: TAIR
# # release date: June 9, 2009
# # 7 sequences:
# #   Chr1 Chr2 Chr3 Chr4 Chr5 ChrM ChrC
# # (use 'seqnames()' to see all the sequence names, use the '$' or '[[' operator to access a given
# # sequence)

Chr1<-BSgenome.Athaliana.TAIR.TAIR9[["Chr1"]]
Chr1<-DNAStringSet(Chr1)
names(Chr1)<-'Chr1'
x<-subseq(Chr1, 1, 1e6)
x
# DNAStringSet object of length 1:
#       width seq                                                             names
# [1] 1000000 CCCTAAACCCTAAACCCTAAACCCTAAACC...CAAGTGCTGAAACGTGTATGATCCGGTTCC Chr1

mep<-MEPmod(genome=x, cores=1)

mep
# GRanges object with 2 ranges and 8 metadata columns:
#       seqnames        ranges strand |                      NT                      AA
#          <Rle>     <IRanges>  <Rle> |          <DNAStringSet>           <AAStringSet>
#   [1]     Chr1 456757-456761      + | TTTGATGCTA...ACTGTCTGAC FDARTAWSQC...AFGAVESLSD
#   [2]     Chr1 601547-601560      - | AAAGTTGAAT...ACATTCAGAT KVEFKDNEWK...RLDCLLKHSD
#              region             block.starts   block.sizes     score isMicroExon   cluster
#           <IRanges>            <IntegerList> <IntegerList> <numeric>   <logical> <integer>
#   [1] 456595-456893 456595,456647,456757,...  52,110,5,...    0.9800        TRUE         8
#   [2] 601327-601700 601653,601561,601547,...  48,92,14,...    0.9241        TRUE        35
#   -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

## Parameters
The arguments of the function ***MEPmod*** are as follows:
| Parameter  | Description  |
| :---------- | :---------- |
| genome      | the path(s) to the fasta file(s) or a 'DNAStringSet' object.  |
| min.score   | a character string containing a percentage specifying the minimum score of each exon block (e.g. "80%"). This parameter will pass to *matchPWM* from **Biostrings** package |
| include.intronLoss | TRUE or FALSE. If TRUE, the microexons with any side of flanking intron loss will also be returned. |
| span        | the maximum spanning region of the microexon-tag (default: 20 kb). |
| min.intron  | minimum intron size  |
| max.intron  | maximum intron size  |
| cores       | number of cores to use. This parameter will pass to *mclapply* as *mc.cores* from **parallel** package. |
#### Warning: The program may be very slow or get stuck for a large genome (> 5 Gb) or a large span region (> 20 kb). 

## Output
The function ***MEPmod*** will returm a GRanges of microexons with 8 metadata columns:
| Column name  | Description |
| :---------- | :---------- |
| NT | 108 bp DNA sequence of microexon-tag (72 bp for cluster 1 and 2) |
| AA | 36 aa protein sequence translated from NT (24 aa for cluster 1 and 2) |
| region | the span region of microexon-tag |
| block.starts | the genomic coordinates of block start positions (exon, intron, exon, ...) |
| block.sizes | block sizes (exon, intron, exon, ...) |
| score | match score of the microexon-tag. Large score indicate high confidence. |
| isMicroExon | if the microexon has any flanking intron loss, FALSE; otherwise, TRUE |
| cluster | which cluster this microexon-tag belongs to |

## Downstream analysis
+ Get the cluster information
```
?MEPdata
MEPdata$cluster[mep$cluster,]
#    cluster size phase        motif exons me_order
# 8        8    5     1 Peptidase_C1     3        2
# 35      35   14     0   Tudor-knot     3        2
```


+ Extract the coordinates of microexons
```
me<-mep[mep$isMicroExon]
mcols(me)<-mcols(me)<-mcols(me)[c('score','cluster')]
me<-as.data.frame(me)

me
#  seqnames  start    end width strand  score cluster
# 1     Chr1 456757 456761     5      + 0.9800       8
# 2     Chr1 601547 601560    14      - 0.9241      35
```

+ Extract intron splicing junctions for RNA-seq mapping (e.g, STAR)
```
junc<-GRangesList(lapply(seq_along(mep), function(i) {
        n<-length(mep$block.starts[[i]])
        GRanges(seqnames(mep[i]),
                IRanges(mep$block.starts[[i]],width = mep$block.sizes[[i]])[seq(2,n,2)],
                strand(mep[i]))
}))@unlistData
junc<-sort(junc[width(junc)>0],ignore.strand=T)
junc<-as.data.frame(junc)[-4]

junc
#  seqnames  start    end strand
# 1     Chr1 456647 456756      +
# 2     Chr1 456762 456842      +
# 3     Chr1 601373 601546      -
# 4     Chr1 601561 601652      -

write.table(junc,file='ME.SJ',row.names = F,col.names = F,quote = F,sep='\t')
```

+ Extract gff file for microexon annotation
```
gff<-GRangesList(lapply(seq_along(mep), function(i) {
    n<-length(mep$block.starts[[i]])
    exon<-GRanges(seqnames(mep[i]),
            IRanges(mep$block.starts[[i]],width = mep$block.sizes[[i]])[seq(1,n,2)],
            strand(mep[i]))
    exon$source<-'MEPmodeler'
    exon$type<-"exon"
    exon$score<-mep$score[i]
    exon$phase<-c(0,cumsum(width(exon)[-length(exon)])%%3)
    exon$Parent<-as.character(mep[i])
    exon
}))@unlistData
gff
# GRanges object with 6 ranges and 5 metadata columns:
#       seqnames        ranges strand |      source        type     score     phase
#          <Rle>     <IRanges>  <Rle> | <character> <character> <numeric> <numeric>
#   [1]     Chr1 456595-456646      + |  MEPmodeler        exon    0.9800         0
#   [2]     Chr1 456757-456761      + |  MEPmodeler        exon    0.9800         1
#   [3]     Chr1 456843-456893      + |  MEPmodeler        exon    0.9800         0
#   [4]     Chr1 601653-601700      - |  MEPmodeler        exon    0.9241         0
#   [5]     Chr1 601547-601560      - |  MEPmodeler        exon    0.9241         0
#   [6]     Chr1 601327-601372      - |  MEPmodeler        exon    0.9241         2
#                     Parent
#                <character>
#   [1] Chr1:456757-456761:+
#   [2] Chr1:456757-456761:+
#   [3] Chr1:456757-456761:+
#   [4] Chr1:601547-601560:-
#   [5] Chr1:601547-601560:-
#   [6] Chr1:601547-601560:-
#   -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths

export.gff3(gff,"ME.gff3")
```
Note: The gff file only contain partly transcripts and the first and the last exons are NOT full.


+ Extract DNA and protein sequences for phylogenetic analysis
```
## DNA sequences
dna<-mep$NT 
names(dna)<-as.character(mep)
writeXStringSet(dna,"dna.fasta")

## Protein sequences
prot<-mep$AA
names(prot)<-as.character(mep)
writeXStringSet(prot,"prot.fasta")
```
## Citation
Huihui Yu, Mu Li, Jaspreet Sandhu, Guangchao Sun, James C. Schnable, Harkamal Walia, Weibo Xie, Bin Yu, Jeffrey P. Mower, Chi Zhang. Pervasive misannotation of microexons that are evolutionarily conserved and crucial for gene function in plants. In Review.

# MEPmodeler: Predict <i>M</i>icro<i>E</i>xons in <i>P</i>lant Genomes

MEPmodeler identifies conserved microexon-tags in plant genomes, where each tag comprises a central microexon (1-51 nt) and its flanking coding sequences (typically 108-nt in total length). The core ***MEPmod*** function uses gapped Position Weight Matrix (PWM) with ***MEPdata*** trained from 12 land plant species. Input requires only plant genomic sequences (FASTA or genome assembly), scanning both DNA strands. This tool enables annotation-free microexon-tag identification for plant phylogemomic analysis, overcoming traditional genome annotation limitations.  
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
The following example is to predict microexons from Arabidopsis Chr1:1-1000000 sequence
for microexon cluster 1:50: 
```{r}
suppressPackageStartupMessages(library(BSgenome))
if (!requireNamespace("BSgenome.Athaliana.TAIR.TAIR9", quietly = TRUE))
    BiocManager::install("BSgenome.Athaliana.TAIR.TAIR9")
library(BSgenome.Athaliana.TAIR.TAIR9)

data <- DNAStringSet(BSgenome.Athaliana.TAIR.TAIR9$Chr1[1:1e6])
names(data) <- "Chr1"

mep <- MEPmod(genome=data, cores=1, clusters=c(1:50), quietly = TRUE)

mep
# GRanges object with 2 ranges and 8 metadata columns:
#        seqnames       ranges strand |                      NT
#          <Rle>     <IRanges>  <Rle> |          <DNAStringSet>
#   [1]     Chr1 456757-456761      + | TTTGATGCTA...ACTGTCTGAC
#   [2]     Chr1 601547-601560      - | AAAGTTGAAT...ACATTCAGAT
#                            AA        region             block.starts   block.sizes
#                 <AAStringSet>     <IRanges>            <IntegerList> <IntegerList>
#   [1] FDARTAWSQC...AFGAVESLSD 456595-456893 456595,456647,456757,...  52,110,5,...
#   [2] KVEFKDNEWK...RLDCLLKHSD 601327-601700 601653,601561,601547,...  48,92,14,...
#            score isMicroExon   cluster
#        <numeric>   <logical> <integer>
#   [1]    0.9817        TRUE         8
#   [2]    0.9242        TRUE        37
#   -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths

```

## Parameters
The arguments of the function ***MEPmod*** are as follows:
| Parameter  | Description  |
| :---------- | :---------- |
| genome      | Path(s) to FASTA file(s) or a "DNAStringSet" object containing genomic sequences.  |
| min.score   | Minimum match score threshold for exon blocks (as percentage string, e.g., "80%"). Passed to *matchPWM* from **Biostrings** package. |
| clusters    | Integer vector specifying which microexon-tag clusters to search (default: NULL searches all clusters). |
| high_complexity | Logical. If TRUE (default), only searches clusters with sequence complexity ≥ 0.7 for both flanking regions. |
| low_copy | Logical. If TRUE, restricts search to low-copy clusters (1-3 copies in >80% of 184 land plants). |
| include.intronLoss | Logical. If TRUE, includes microexons with any flanking intron loss. |
| span        | Maximum genomic span for microexon-tag detection (default: 20,000 bp). |
| min.intron  | Minimum allowed intron size (bp). |
| max.intron  | Maximum allowed intron size (bp). |
| cores       | Number of CPU cores to use. Passed to *mclapply* as *mc.cores* from **parallel** package. |
#### Warning: The program may be very slow or get stuck for a large genome (> 5 Gb) or a large span region (> 20 kb). 

## Output
The function ***MEPmod*** will return a GRanges object with 8 metadata columns:
| Column name  | Description |
| :---------- | :---------- |
| NT | DNA sequence of the microexon-tag |
| AA | Peptide sequence encoded by the microexon-tag |
| region | Genomic span of the microexon-tag |
| block.starts | Genomic coordinates of exon/intron boundaries |
| block.sizes | Lengths of exons and introns |
| score | Match confidence score (higher = more confident) |
| isMicroExon | TRUE if no flanking intron loss, FALSE otherwise |
| cluster | Cluster ID of the detected microexon-tag |

## Downstream analysis
+ Get the cluster information
```
?MEPdata
MEPdata$cluster[mep$cluster,]
#    cluster size phase        motif exons me_order rep_complex
# 8        8    5     1 Peptidase_C1     3        2   0.8775100
# 37      37   14     0   Tudor-knot     3        2   0.8694779
#                                                                                                          string
# 8  TTNGATGCNNGAACNGCTTGGNCTCANTGNANCACNATTGGNANNATNCTNGATCAGGGNCANTGTGGTTCTTGNTGGGCNTTTGGTGCTGTNGANTCACTNTCNGAT
# 37 AANGNTGANNTNCGNAAGAANGANTGGANATANTTNGTNCATTANCTTGGTTGGANNAAAAANTGGGATGAATGGGTNGGNNNNGANCGNNTGNTGAANNNNACTGAN
#    L_complex R_complex low_copy
# 8  0.8980257 0.8369434    FALSE
# 37 0.8500758 0.8109660     TRUE
```


+ Extract the coordinates of microexons
```
me <- mep[mep$isMicroExon]
mcols(me) <- mcols(me) <- mcols(me)[c('score','cluster')]
me <- as.data.frame(me)

me
#  seqnames  start    end width strand  score cluster
# 1     Chr1 456757 456761     5      + 0.9817       8
# 2     Chr1 601547 601560    14      - 0.9242      37
```

+ Extract intron splicing junctions for RNA-seq mapping (e.g, STAR)
```
junc <- GRangesList(lapply(seq_along(mep), function(i) {
                n <- length(mep$block.starts[[i]])
                GRanges(seqnames(mep[i]),
                       IRanges(mep$block.starts[[i]],
                               width = mep$block.sizes[[i]])[seq(2,n,2)],
                               strand(mep[i]))
}))@unlistData
junc <- sort(junc[width(junc)>0], ignore.strand=T)
junc <- as.data.frame(junc)[-4]

junc
#  seqnames  start    end strand
# 1     Chr1 456647 456756      +
# 2     Chr1 456762 456842      +
# 3     Chr1 601373 601546      -
# 4     Chr1 601561 601652      -

write.table(junc,file='ME_SJ.tab',row.names = F,col.names = F,quote = F,sep='\t')
```

+ Extract gff file for microexon-tag annotation
```
gff <- GRangesList(lapply(seq_along(mep), function(i) {
    n <- length(mep$block.starts[[i]])
    exon <- GRanges(seqnames(mep[i]),
                    IRanges(mep$block.starts[[i]],
                            width = mep$block.sizes[[i]])[seq(1,n,2)],
                            strand(mep[i]))
    exon$source <- 'MEPmodeler'
    exon$type <- "CDS"
    exon$score <- mep$score[i]
    exon$phase <- c(0,cumsum(width(exon)[-length(exon)])%%3)
    exon$Parent <- as.character(mep[i])
    exon
}))@unlistData

gff
# GRanges object with 6 ranges and 5 metadata columns:
#       seqnames        ranges strand |      source        type     score
#          <Rle>     <IRanges>  <Rle> | <character> <character> <numeric>
#   [1]     Chr1 456595-456646      + |  MEPmodeler        exon    0.9817
#   [2]     Chr1 456757-456761      + |  MEPmodeler        exon    0.9817
#   [3]     Chr1 456843-456893      + |  MEPmodeler        exon    0.9817
#   [4]     Chr1 601653-601700      - |  MEPmodeler        exon    0.9242
#   [5]     Chr1 601547-601560      - |  MEPmodeler        exon    0.9242
#   [6]     Chr1 601327-601372      - |  MEPmodeler        exon    0.9242
#           phase               Parent
#       <numeric>          <character>
#   [1]         0 Chr1:456757-456761:+
#   [2]         1 Chr1:456757-456761:+
#   [3]         0 Chr1:456757-456761:+
#   [4]         0 Chr1:601547-601560:-
#   [5]         0 Chr1:601547-601560:-
#   [6]         2 Chr1:601547-601560:-
#   -------
#   seqinfo: 1 sequence from an unspecified genome; no seqlengths

export.gff3(gff,"ME.gff3")
```
Note: The gff file only contain partly transcripts and the first and the last exons are NOT full.


+ Extract DNA or protein sequences for phylogenomic analysis
```
## DNA sequences
dna <- mep$NT 
names(dna) <- as.character(mep)
writeXStringSet(dna,"dna.fa")

## Protein sequences
prot <- mep$AA
names(prot) <- as.character(mep)
writeXStringSet(prot,"prot.fa")
```
## Citation
Huihui Yu, Mu Li, Jaspreet Sandhu, Guangchao Sun, James C. Schnable, Harkamal Walia, Weibo Xie, Bin Yu, Jeffrey P. Mower, Chi Zhang. Pervasive misannotation of microexons that are evolutionarily conserved and crucial for gene function in plants. *Nat Commun* **13**, 820 (2022). https://doi.org/10.1038/s41467-022-28449-8

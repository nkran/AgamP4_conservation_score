library(rtracklayer)
library(Rsamtools)
library(CNEr)
library(argparse)
library(GenomeInfoDb)
library(dplyr)

parser <- ArgumentParser()
parser$add_argument("-j", "--job", type="character", help="Genome base")

args <- parser$parse_args()
genome <- strsplit(args$job, '__')[[1]][1]
identity <- strsplit(args$job, '__')[[1]][2]

knitr::opts_chunk$set(echo=FALSE, warning=debug, message=debug, error=FALSE,
                      cache.path = "cache/",
                      fig.path = "figures/")

axtFilesSubj <- paste("/rds/general/user/nkranjc/home/conservation/data/alignments/axt_final/AgamP4.", genome, ".net.axt", sep='')
axtFilesSubj_r <- paste("/rds/general/user/nkranjc/home/conservation/data/alignments/axt_final/", genome, ".AgamP4.net.axt", sep='')

axtSubj <- readAxt(axtFilesSubj,
                   tAssemblyFn='/rds/general/user/nkranjc/home/conservation/data/alignments/ref/AgamP4.fa',
                   qAssemblyFn=paste("/rds/general/user/nkranjc/home/conservation/data/alignments/ref/", genome, ".fa", sep=''))

axtSubj_r <- readAxt(axtFilesSubj_r,
                     tAssemblyFn=paste("/rds/general/user/nkranjc/home/conservation/data/alignments/ref/", genome, ".fa", sep=''),
                     qAssemblyFn='/rds/general/user/nkranjc/home/conservation/data/alignments/ref/AgamP4.fa')

identities <- c(27L, 29L, 30L, 35L, 45L, 48L, 49L, 50L)
windows <-    c(30L, 30L, 30L, 50L, 50L, 50L, 50L, 50L)

cneSubj <- CNE(
  assembly1Fn=file.path('/rds/general/user/nkranjc/home/conservation/data/alignments/ref/AgamP4.2bit'),
  assembly2Fn=file.path(paste("/rds/general/user/nkranjc/home/conservation/data/alignments/ref/", genome, ".2bit", sep='')),
  axt12Fn=axtFilesSubj,
  axt21Fn=axtFilesSubj_r,
  cutoffs1=8L, cutoffs2=4L)

cneList <- ceScan(x=cneSubj,
            window=windows,
            tSizes=seqlengths(AgamP4),
            qSizes=seqlengths(genome),
            identity=identities)

cneMerged <- lapply(cneList, cneMerge)

print('--------- blat')
print(identity)
cneFinal <- blatCNE(cneMerged[[identity]], cutIdentity=70)
cneFinal <- CNE12(cneMerged[[identity]])
write.csv(cneFinal, paste0('/rds/general/user/nkranjc/home/conservation/data/alignments/AgamP4.',genome, '_d', identity,".csv"))
print('----------- finished')
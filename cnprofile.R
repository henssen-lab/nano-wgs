library(QDNAseq)
args = commandArgs(trailingOnly=TRUE)
bam_fname = args[1]
pdf_fname = args[2]
rdata_fname = args[3]
#bam_fname = "~/Desktop/nano-wgs/8734T-alignment/8734T.pass.minimap.sorted.bam"
bins <- getBinAnnotations(binSize = 100, genome="hg19")
r <- binReadCounts(bins, bamfiles = bam_fname, minMapq=20, isSecondaryAlignment=FALSE)
r <- applyFilters(r, residual = F, blacklist = T, mappability = F, bases = 100, chromosomes = "Y")
r <- estimateCorrection(r)
r <- correctBins(r, method = "median")
r <- normalizeBins(r, method = "median")
r <- segmentBins(r, transformFun="sqrt", alpha = 0.01, min.width=2, undo.splits="none")
r_bins_called = callBins(r)
save.image(rdata_fname)
pdf(pdf_fname, height=3, width=6)
plot(r)
dev.off()

#cnc <- callBins(r)
#pdf(paste0(bam_fname, ".cncall.pdf"))
#plot(cnc)
#dev.off()

#exportBins(r, file="r.igv", format="igv")

library(ggplot2)

data_header = c("chrom", "pos", "qname", "flag", "mis", "rlen")
tmp<-read.table("/home/gaoy1/hg002/all.mis.sites.out", header=FALSE)

data <- tmp
colnames(data) <- data_header

chroms <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM")
# only keep chrom chr1-22, chrX, chrY
data <- data[data$chrom %in% chroms,]

# sort chrom by chr1-22, chrX, chrY, chrM
data$chrom <- factor(data$chrom, levels = chroms)
# set flag as 'primary' if 0 or 16
# set flag as 'supplementary' if 2048 or 2064
data$flag <- ifelse(data$flag %in% c(0, 16), "primary", "supplementary")


ggplot(data = data, aes(x = chrom, y = rlen/(mis+1), fill = flag)) +
    geom_violin() +
    geom_hline(yintercept = 100, linewidth=1, color="red") +
    scale_y_continuous(trans='log10') +
    labs(x="chromsome", y = "rlen / (#mismatch+1)") +
    theme_classic() +
    guides(fill=guide_legend(ncol=2, position = "top")) +
    theme(text = element_text(size=15))

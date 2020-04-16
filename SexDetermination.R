#Male specific region analysis
#h67_w[h67_w[,1]=="scaffold149",]
#male_heterozygous_region <- plot_top_regions_2(2659, h67_w, herring_67$geno, h67_clean_samples, "~/Projects/Herring/doc/SexDet/", snp_col="blue", snp_lwd = 0.25, snp_alpha = 120)

#v2.0.2 version
load("~/Projects/Herring/data/v2.0.2_genotypes/herring_79.RData")
load("~/Projects/Herring/data/v2.0.2_genotypes/h79_sw.RData")

h79_w[h79_w[,1]=="chr8" & h79_w$stop > 21.0e6 & h79_w$start < 21.3e6,]
male_heterozygous_region <- plot_top_regions_2(2464, h79_w, herring_79$geno, herring_79$sample_list, "~/Projects/Herring/doc/SexDet/", snp_col="blue", snp_lwd = 0.25, snp_alpha = 120)
h79_w[h79_w[,1]=="unplaced_scaffold18",]
us18_HapDist <- plot_top_regions_2(7327, h79_w, herring_79$geno, herring_79$sample_list, "~/Projects/Herring/doc/SexDet/", snp_col="blue", snp_lwd = 0.25, snp_alpha = 120)

#"Extra" CATSPER3
#22259379 - alignment start
h79_w[h79_w[,1]=="chr8" & h79_w$stop > 22.2e6 & h79_w$start < 22.3e6,]
tmp_w <-h79_w[2474:2476,] 
tmp_w[1,"stop"] <- 22259379 -1
tmp_w[2,"start"] <- 22259379
tmp_w[2,"stop"] <- 22268294
tmp_w[3,"start"] <- 22268294 +1
plot_top_regions_2(2, tmp_w, herring_79$geno, herring_79$sample_list, "~/Projects/Herring/doc/SexDet/", snp_col="blue", snp_lwd = 0.25, snp_alpha = 120)

#sorted_haplotypes <- h67_clean_samples[male_heterozygous_region[[1]]$order]
sorted_haplotypes <- herring_79$sample_list[male_heterozygous_region[[1]]$order]
male_co <- which(sorted_haplotypes == "HWS-33_1")
female_co_1 <- which(sorted_haplotypes == "S4TH-231_2")
female_co_2 <- which(sorted_haplotypes == "AM8_1")
pac_female_co <- which(sorted_haplotypes == "HWS-21_2")

haplo_classification <- array(data = "ND", dim = length(sorted_haplotypes))
haplo_classification[1:pac_female_co] <- "pac_female"
haplo_classification[female_co_2:female_co_1] <- "atl_female"
haplo_classification[male_co:length(sorted_haplotypes)] <- "male"
ind_df <- data.frame(ind_list = unique(sub("_[12]", "", h67_clean_samples)), hap1 = "ND", hap2 = "ND", stringsAsFactors =F)
for (i in  1:dim(ind_df)[1]){
	ind_df[i,2:3] <- haplo_classification[grep(ind_df[i,1], sorted_haplotypes)]
}


#Blasting unplaced scaffolds 229 & 18 to the genome, in order to start estimating divergence
writeXStringSet(Ch_v2.0.2[paste0("unplaced_scaffold", c(18, 229))], file = "~/Projects/Herring/data/v2.0.2_annotation/sex_determination/unplaced_scaffold18_229.fasta")
#blastn -task blastn -db /crex/proj/uppstore2017191/private/assemblies/Ch_v2.0.2.fasta -query ./unplaced_scaffold18_229.fasta -outfmt 6 -num_threads 4 -evalue 1e-10 -culling_limit 5 > unplaced_scaffold18_229.blastout 
SexDet_blast_raw <- read.table("~/Projects/Herring/data/v2.0.2_annotation/sex_determination/unplaced_scaffold18_229.blastout", stringsAsFactors=F, header = F)
SexDet_blast <- SexDet_blast_raw[SexDet_blast_raw[,2] == "chr8",]
SexDet_blast[,"col"] <- "firebrick"
SexDet_blast[SexDet_blast[,1] == "unplaced_scaffold18","col"] <- "darkorchid"


pdf(file = "~/Projects/Herring/doc/SexDet/Blast_chr8_21.00_21.25_Mb.pdf", width = 14, height = 4)
plot(x = 0, y = 0, type = "n", xlim = c(21.0e6, 21.25e6), ylim = c(0, 2.5e5), xlab = "Chr 8", ylab = "Unplaced scaffold")
segments(x0 = SexDet_blast[,9], x1 = SexDet_blast[,10], y0 = SexDet_blast[,7], y1 = SexDet_blast[,8], col = SexDet_blast[,"col"])
abline(v = 21.051e6)
abline(v = 21.065e6)
abline(v = 21.068e6)
abline(v = 21.090e6)
abline(v = 21.094e6)
abline(v = 21.128e6)
abline(v = 21.159e6)
abline(v = 21.179e6)
dev.off()

#Blasting the chr8-verison of the locus, to find possible additional fragments
writeXStringSet(subseq(Ch_v2.0.2["chr8"], start = 21.0e6, end = 21.2e6), file = "~/Projects/Herring/data/v2.0.2_annotation/sex_determination/chr8_21.0_21.2_MB.fasta")
chr8_blast_raw <- read.table("~/Projects/Herring/data/v2.0.2_annotation/sex_determination/chr8_21.0_21.2_MB.blastout", stringsAsFactors=F, header = F)

chr8_blast <- chr8_blast_raw[chr8_blast_raw[,4] >= 1e3,]

#Setting up clustal omega for each aligned block
which(SexDet_blast$V9 > 21.051e6 & SexDet_blast$V10 < 21.065e6)
SexDet_blast[1,]
#Aln_seq <- DNAStringSet()
#aln_target <- GRanges(seqnames = c("unplaced_scaffold18", "chr8"), ranges = IRanges(start = as.numeric(SexDet_blast[1,c(7,10)]), end = as.numeric(SexDet_blast[1,c(8,9)])))
#Aln_seq <- Ch_v2.0.2[aln_target]
#names(Aln_seq) <- paste(aln_target)
#Aln_seq[1] <- reverseComplement(Aln_seq[1])
#Aln_seq_bin <- as.DNAbin(Aln_seq)
#Region1_clustalo <- clustalomega(Aln_seq_bin, exec="clustalo", quiet = F)
#checkAlignment(Region1_clustalo)
reg1_seqnames <- c("unplaced_scaffold18", "chr8")
reg1_starts <- c(188923, 21050480)
reg1_ends <- c(202819, 21064424)
reg1_clustalo <- align_target_region(seqname_vec = reg1_seqnames, start_vec = reg1_starts, stop_vec = reg1_ends)
pdf(file = "~/Projects/Herring/doc/SexDet/ClustalO_chr8_21.051_21.065_Mb.pdf", width = 14, height = 7)
checkAlignment(reg1_clustalo)
dev.off()
dist.dna(reg1_clustalo)
aln_diff_plot(reg1_clustalo)
aln_diff_plot_2(reg1_clustalo, win_size = 200, pdf_file = "~/Projects/Herring/doc/SexDet/reg1_diff.pdf", y_lim = c(-0.5,0.1))



SexDet_blast[which(SexDet_blast$V9 > 21.068e6 & SexDet_blast$V10 < 21.090e6),]
#us18_range <- c(38837, 62223)
#chr8_range <- c(21068299, 21089937)
#aln_target <- GRanges(seqnames = c("unplaced_scaffold18", "chr8"), ranges = IRanges(start = c(us18_range[1], chr8_range[1]), end = c(us18_range[2], chr8_range[2])))
#Aln_seq <- Ch_v2.0.2[aln_target]
#names(Aln_seq) <- paste(aln_target)
#Aln_seq[1] <- reverseComplement(Aln_seq[1])
#Aln_seq_bin <- as.DNAbin(Aln_seq)
#Region2_clustalo <- clustalomega(Aln_seq_bin, exec="clustalo", quiet = F)
reg2_seqnames <- c("unplaced_scaffold18", "chr8")
reg2_starts <- c(38837, 21068299)
reg2_ends <- c(62223, 21089937)
reg2_clustalo <- align_target_region(seqname_vec = reg2_seqnames, start_vec = reg2_starts, stop_vec = reg2_ends)
pdf(file = "~/Projects/Herring/doc/SexDet/ClustalO_chr8_21.068_21.090_Mb.pdf", width = 14, height = 7)
checkAlignment(reg2_clustalo)
dev.off()
dist.dna(reg2_clustalo)
aln_diff_plot(reg2_clustalo)
aln_diff_plot_2(reg2_clustalo, win_size = 200, pdf_file = "~/Projects/Herring/doc/SexDet/reg2_diff.pdf", y_lim = c(-1,0.5))

SexDet_blast[which(SexDet_blast$V9 > 21.094e6 & SexDet_blast$V10 < 21.121e6),]
reg3_seqnames <- c("unplaced_scaffold18", "chr8")
reg3_starts <- c(9996, 21094131)
reg3_ends <- c(38003, 21121514)
reg3_clustalo <- align_target_region(seqname_vec = reg3_seqnames, start_vec = reg3_starts, stop_vec = reg3_ends)
pdf(file = "~/Projects/Herring/doc/SexDet/ClustalO_chr8_21.094_21.121_Mb.pdf", width = 14, height = 7)
checkAlignment(reg3_clustalo)
dev.off()
dist.dna(reg3_clustalo)
aln_diff_plot(reg3_clustalo)
aln_diff_plot_2(reg3_clustalo, win_size = 200, pdf_file = "~/Projects/Herring/doc/SexDet/reg3_diff.pdf", y_lim = c(-1,0.1))
#Divergence-based age estimate 2*t*u = d
d <- dist.dna(reg3_clustalo) - 0.003 #Getting excess divergence by subtracting average neutral diverity of 0.3%
u <- 2.0e-9 #From Feng et al 2017
t <- d/(2*u) #in generations
t_y <- t*6 #in years, according to generations time estimate from Feng et al 2017
#t_y/1e6: 7.542282


SexDet_blast[which(SexDet_blast$V9 > 21.159e6 & SexDet_blast$V10 <  21.179e6),]
reg4_seqnames <- c("unplaced_scaffold229", "chr8")
reg4_starts <- c(46914, 21159122)
reg4_ends <- c(69516, 21179333)
reg4_clustalo <- align_target_region(seqname_vec = reg4_seqnames, start_vec = reg4_starts, stop_vec = reg4_ends)
pdf(file = "~/Projects/Herring/doc/SexDet/ClustalO_chr8_21.159_21.179_Mb.pdf", width = 14, height = 7)
checkAlignment(reg4_clustalo)
dev.off()
dist.dna(reg4_clustalo)
aln_diff_plot(reg4_clustalo)
aln_diff_plot_2(reg4_clustalo, win_size = 200, pdf_file = "~/Projects/Herring/doc/SexDet/reg4_diff.pdf", y_lim = c(-1,0.1))

reg2_3_seqnames <- c("unplaced_scaffold18", "chr8")
reg2_3_starts <- c(9996, 21068299)
reg2_3_ends <- c(62223, 21121514)
reg2_3_clustalo <- align_target_region(seqname_vec = reg2_3_seqnames, start_vec = reg2_3_starts, stop_vec = reg2_3_ends)
pdf(file = "~/Projects/Herring/doc/SexDet/ClustalO_chr8_21.068_21.121_Mb.pdf", width = 14, height = 7)
checkAlignment(reg2_3_clustalo)
dev.off()
dist.dna(reg2_3_clustalo)
aln_diff_plot_2(reg2_3_clustalo, win_size = 200, pdf_file = "~/Projects/Herring/doc/SexDet/reg_2_3_diff.pdf")


#bmpr1bb - introns
bmpr1bb_seqnames <- c("unplaced_scaffold229", "chr21")
bmpr1bb_full_starts <- c(1, 17016119)
bmpr1bb_full_ends <- c(40000, 17052960)
bmpr1bb_starts <- c(16200, 17016119+30661)
bmpr1bb_ends <- c(21500, 17016119+35975)

bmpr1bb_clustalo <- align_target_region(seqname_vec = bmpr1bb_seqnames, start_vec = bmpr1bb_starts, stop_vec = bmpr1bb_ends)
aln_diff_plot(bmpr1bb_clustalo, count_gaps = T, win_size = 50)
checkAlignment(bmpr1bb_clustalo)

bmpr1bb_full_clustalo <- align_target_region(seqname_vec = bmpr1bb_seqnames, start_vec = bmpr1bb_full_starts, stop_vec = bmpr1bb_full_ends)
checkAlignment(bmpr1bb_full_clustalo)
dist.dna(bmpr1bb_full_clustalo)
bmpr1bb_full_avg_diff <- aln_diff_plot(bmpr1bb_full_clustalo, count_gaps = T, win_size = 100)


bmpr1bb_blast_raw <- read.table("~/Projects/Herring/data/v2.0.2_annotation/sex_determination/bmpr1bb_chr21_genomic.blastout", stringsAsFactors=F, header = F)
bmpr1bb_blast <- bmpr1bb_blast_raw[bmpr1bb_blast_raw$V2 == "unplaced_scaffold229",]
bmpr1bb_gtf <- cluhar_v2.0.2_gtf[cluhar_v2.0.2_gtf$gene_id == "ENSCHAG00000027155"]

#add alignment percentage on adjusted x-axis
aln_pos_chr21 <- which(as.character(bmpr1bb_full_clustalo[2,]) != "-")
checkAlignment(bmpr1bb_full_clustalo[,aln_pos_chr21])
aln_pos_chr21_clustalo <- bmpr1bb_full_clustalo[,aln_pos_chr21]
aln_pos_chr21_avg_diff <- aln_diff_plot(aln_pos_chr21_clustalo, count_gaps = F, win_size = 50)
aln_pos_chr21_gaps <- as.numeric(as.character(aln_pos_chr21_clustalo[1,]) == "-")
aln_pos_chr21_avg_gaps <- filter(aln_pos_chr21_gaps, rep(1/50, 50))

#points(x= (1:length(aln_pos_chr21)) + 17016119, y = bmpr1bb_full_avg_diff[aln_pos_chr21]/2, pch = 16, cex = 0.5)
#plot(x= (1:length(aln_pos_chr21)) + 17016119, y = bmpr1bb_full_avg_diff[aln_pos_chr21], pch = 16, cex = 0.5)

checkAlignment(aln_pos_chr21_clustalo[,25000:length(aln_pos_chr21)])
double_aln_pos <- which(as.character(aln_pos_chr21_clustalo[1,]) != "-")
double_aln_avg_diff <- aln_diff_plot(aln_pos_chr21_clustalo[,double_aln_pos], count_gaps = F, win_size = 50)

pdf(file = "~/Projects/Herring/doc/SexDet/bmpr1bb_full.pdf", width = 10, height = 7)
plot_gtf(bmpr1bb_gtf)
points(x= double_aln_pos + 17016119, y = double_aln_avg_diff - 0.25, pch = 16, cex = 0.5)
dev.off()

pdf(file = "~/Projects/Herring/doc/SexDet/bmpr1bb_17.043_17.054_Mb.pdf", width = 10, height = 7)
gtf_subset <- c(1,which(bmpr1bb_gtf$transcript_id == "ENSCHAT00000065423"))
plot_gtf(bmpr1bb_gtf[gtf_subset], x_lim = c(17043000, 17054000), y_lim = c(-0.5, 1.5))
abline(v = bmpr1bb_gtf[gtf_subset]@ranges@start, lwd  =.5, col = "grey30")
abline(v = bmpr1bb_gtf[gtf_subset]@ranges@start +  bmpr1bb_gtf[gtf_subset]@ranges@width, lwd  =.5, col = "grey30")
abline(h = c(0.75, 0.25, -0.25), col = "grey70")
abline(h = c(0.5, 0), col = "grey70", lty = "dashed")
segments(x0 = 17016119 + bmpr1bb_blast[,7], x1 = 17016119 + bmpr1bb_blast[,8], y0 = 1.5, lwd = 3, col = "blue")
#points(x= (1:length(aln_pos_chr21)) + 17016119, y = aln_pos_chr21_avg_diff/2, pch = 16, cex = 0.5)
points(x= double_aln_pos + 17016119, y = (double_aln_avg_diff/2)+0.25, pch = 16, cex = 0.5, col = "darkorchid")
points(x= (1:length(aln_pos_chr21)) + 17016119, y = -(aln_pos_chr21_avg_gaps/2) + .25, pch = 16, cex = 0.5, col = "olivedrab")
dev.off()

#Breaking out segements for possible age estimates
double_aln_GR <- GRanges(seqnames = "21", ranges = IRanges(start = double_aln_pos + 17016119, end = double_aln_pos + 17016119))
double_aln_in_CDS <- findOverlaps(double_aln_GR, bmpr1bb_gtf[cds_subset], ignore.strand = T)

n_exons <- length(unique(double_aln_in_CDS@to))
exon_mean_diff <- data.frame(exon = 1:n_exons)
diff_vec <- as.numeric(aln_pos_chr21_clustalo[1,double_aln_pos] != aln_pos_chr21_clustalo[2,double_aln_pos]) 
for(exon in 1:n_exons){
  exon_mean_diff[exon_mean_diff$exon == exon, "mean_diff"] <- mean(diff_vec[double_aln_in_CDS@from[double_aln_in_CDS@to == exon]])
  exon_chr21_aln_start <- bmpr1bb_gtf[cds_subset]@ranges@start[exon] - 17016119
  exon_chr21_aln_end <- bmpr1bb_gtf[cds_subset]@ranges@start[exon] +  bmpr1bb_gtf[cds_subset]@ranges@width[exon]- 17016119
  exon_mean_diff[exon_mean_diff$exon == exon, "gap_count"] <- sum(aln_pos_chr21_gaps[exon_chr21_aln_start:exon_chr21_aln_end]) #/bmpr1bb_gtf[cds_subset]@ranges@width[exon]
  exon_mean_diff[exon_mean_diff$exon == exon, "insertion_count"] <- sum(as.character(bmpr1bb_full_clustalo[2,aln_pos_chr21[exon_chr21_aln_start]:aln_pos_chr21[exon_chr21_aln_end]]) == "-") #/bmpr1bb_gtf[cds_subset]@ranges@width[exon]
  exon_mean_diff[exon_mean_diff$exon == exon, "exon_length_Y"] <- sum(as.character(bmpr1bb_full_clustalo[1,aln_pos_chr21[exon_chr21_aln_start]:aln_pos_chr21[exon_chr21_aln_end]]) != "-")
  exon_mean_diff[exon_mean_diff$exon == exon, "exon_length_chr21"] <- sum(as.character(bmpr1bb_full_clustalo[2,aln_pos_chr21[exon_chr21_aln_start]:aln_pos_chr21[exon_chr21_aln_end]]) != "-")
  if(exon < n_exons){
    intron_end <- aln_pos_chr21[bmpr1bb_gtf[cds_subset]@ranges@start[exon+1] - 1 - 17016119]
    intron_start <- aln_pos_chr21[bmpr1bb_gtf[cds_subset]@ranges@start[exon] + bmpr1bb_gtf[cds_subset]@ranges@width[exon] - 17016119]
    exon_mean_diff[exon_mean_diff$exon == exon, "next_intron_Y"] <- sum(as.character(bmpr1bb_full_clustalo[1,intron_start:intron_end]) != "-")
    exon_mean_diff[exon_mean_diff$exon == exon, "next_intron_chr21"] <- sum(as.character(bmpr1bb_full_clustalo[2,intron_start:intron_end]) != "-")
  }
}

plot_CDS_aln(5,full_clustalo = bmpr1bb_full_clustalo, chr_aln_pos = aln_pos_chr21, CDS_gtf = bmpr1bb_gtf[cds_subset])
plot_CDS_aln(6,full_clustalo = bmpr1bb_full_clustalo, chr_aln_pos = aln_pos_chr21, CDS_gtf = bmpr1bb_gtf[cds_subset])
plot_CDS_aln(7,full_clustalo = bmpr1bb_full_clustalo, chr_aln_pos = aln_pos_chr21, CDS_gtf = bmpr1bb_gtf[cds_subset])
plot_CDS_aln(8,full_clustalo = bmpr1bb_full_clustalo, chr_aln_pos = aln_pos_chr21, CDS_gtf = bmpr1bb_gtf[cds_subset])
plot_CDS_aln(9,full_clustalo = bmpr1bb_full_clustalo, chr_aln_pos = aln_pos_chr21, CDS_gtf = bmpr1bb_gtf[cds_subset])
plot_CDS_aln(10,full_clustalo = bmpr1bb_full_clustalo, chr_aln_pos = aln_pos_chr21, CDS_gtf = bmpr1bb_gtf[cds_subset])
plot_CDS_aln(11,full_clustalo = bmpr1bb_full_clustalo, chr_aln_pos = aln_pos_chr21, CDS_gtf = bmpr1bb_gtf[cds_subset])

#Refined version for figure
pdf(file = "~/Projects/Herring/doc/SexDet/bmpr1bb_17.043_17.054_Mb_annot.pdf", width = 10, height = 7)
gtf_subset <- c(1,which(bmpr1bb_gtf$transcript_id == "ENSCHAT00000065423"))
plot_gtf(bmpr1bb_gtf[gtf_subset], x_lim = c(17043000, 17054000), y_lim = c(-0.5, 1.5))
segments(x0 = bmpr1bb_gtf[gtf_subset]@ranges@start, y0 = 1.2, y1 = -0.25, lwd  =.5, col = "grey30")
segments(x0 = bmpr1bb_gtf[gtf_subset]@ranges@start +  bmpr1bb_gtf[gtf_subset]@ranges@width, y0 = 1.2, y1 = -0.25, lwd  =.5, col = "grey30")
abline(h = c(0.75, 0.25, -0.25), col = "grey70")
abline(h = c(0.5, 0), col = "grey70", lty = "dashed")
#segments(x0 = 17016119 + bmpr1bb_blast[,7], x1 = 17016119 + bmpr1bb_blast[,8], y0 = 1.5, lwd = 3, col = "blue")
points(x= double_aln_pos + 17016119, y = (double_aln_avg_diff/2)+0.25, pch = 16, cex = 0.5, col = "darkorchid")
points(x= (1:length(aln_pos_chr21)) + 17016119, y = -(aln_pos_chr21_avg_gaps/2) + .25, pch = 16, cex = 0.5, col = "olivedrab")
cds_subset <- bmpr1bb_gtf$transcript_id == "ENSCHAT00000065423" & bmpr1bb_gtf$type == "CDS"
text(x = bmpr1bb_gtf[cds_subset]@ranges@start + bmpr1bb_gtf[cds_subset]@ranges@width/2, y = rep(c(1.3, 1.5), length.out = sum(cds_subset)), labels = paste("Exon", 1:sum(cds_subset), "\nd = ", round(exon_mean_diff[, "mean_diff"], digits = 2), sep = ""), cex = 1.1)

dev.off()
###


#Divergence-based age estimate 2*t*u = d
#Most divergent but still recogniseable exon (6)
d_e6 <- exon_mean_diff[6,2]- 0.003 #Getting excess divergence by subtracting average neutral diverity of 0.3%
u <- 2.0e-9 #From Feng et al 2017
t_e6 <- d_e6/(2*u) #in generations
ty_e6 <- t_e6*6 #in years, according to generations time estimate from Feng et al 2017
#ty_e6/1e6: 292.2033
#Least divergent exon (9)
d_e9 <- exon_mean_diff[9,2]- 0.003 #Getting excess divergence by subtracting average neutral diverity of 0.3%
u <- 2.0e-9 #From Feng et al 2017
t_e9 <- d_e9/(2*u) #in generations
ty_e9 <- t_e9*6 #in years, according to generations time estimate from Feng et al 2017
#ty_e9/1e6: 80.72727
#Not useable, too much relies on mutation rate & generationtime



#dN/dS calculations
#Get exon positing from gtf & alignment start
#17051112-17016119 = 34993
#17051242-17016119 = 35123
chr21_offset <- 17016119
aln_interval <- numeric()
for(i in 7:length(bmpr1bb_gtf[cds_subset])){
  aln_interval <- c(aln_interval, bmpr1bb_gtf[cds_subset]@ranges@start[i]:(bmpr1bb_gtf[cds_subset]@ranges@start[i] + bmpr1bb_gtf[cds_subset]@ranges@width[i] -1) - chr21_offset)
}

#aln_interval <- c(34993:35123)
paste(as.character(trans(x=aln_pos_chr21_clustalo[1,aln_interval], codonstart = 2)), collapse = "") #Checking frame
paste(as.character(trans(x=aln_pos_chr21_clustalo[2,aln_interval], codonstart = 2)), collapse = "") #Checking frame
checkAlignment(aln_pos_chr21_clustalo[,c(aln_interval, max(aln_interval)+1)])
dnds(x=aln_pos_chr21_clustalo[,c(aln_interval, max(aln_interval)+1)],codonstart = 32) #Start & extra nucleotide to match Nima's gene model
dist.aa(x=trans(aln_pos_chr21_clustalo[,c(aln_interval, max(aln_interval)+1)], codonstart = 32))
stringdist(paste(as.character(aln_pos_chr21_clustalo[1,c(aln_interval, max(aln_interval)+1)]), collapse = ""), paste(as.character(aln_pos_chr21_clustalo[2,c(aln_interval, max(aln_interval)+1)]), collapse = ""))
dist.dna(x=aln_pos_chr21_clustalo[,c(aln_interval, max(aln_interval)+1)], model = "raw")
#####

pdf(file = "~/Projects/Herring/doc/SexDet/bmpr1bb_diff.pdf", width = 10, height = 7)
plot(x= double_aln_pos  + 17016119, y = double_aln_avg_diff, xlim = c(17043000, 17054000), ylim = c(-1,1), type = "n")
gtf_tmp <- bmpr1bb_gtf[1:length(bmpr1bb_gtf) %in% gtf_subset & bmpr1bb_gtf$type == "CDS"]
rect(xleft = gtf_tmp@ranges@start, xright = gtf_tmp@ranges@start + gtf_tmp@ranges@width, ybottom = -1.5, ytop = 1.5, border = NA, lwd  =.5, col = "grey80")
abline(h = 0, col = "grey70")
abline(h = seq(from = 0.1, to = 0.5, by = 0.1), col = "grey80")
points(x= (1:length(aln_pos_chr21)) + 17016119, y = -aln_pos_chr21_avg_gaps, pch = 16, cex = 0.5, col = "olivedrab")
points(x= double_aln_pos  + 17016119, y = double_aln_avg_diff, pch = 16, cex = 0.5, col = "darkorchid")
dev.off()

#bmpr1ba - introns
bmpr1ba_seqnames <- c("unplaced_scaffold229", "chr7")
bmpr1ba_starts <- c(17487, 12276563+2910)
bmpr1ba_ends <- c(19860, 12276563+5277)

bmpr1ba_full_starts <- c(1,12276563)
bmpr1ba_full_ends <- c(40000,12323175)

bmpr1ba_blast_raw <- read.table("~/Projects/Herring/data/v2.0.2_annotation/sex_determination/bmpr1ba_chr7_genomic.blastout", stringsAsFactors=F, header = F)
bmpr1ba_blast <- bmpr1ba_blast_raw[bmpr1ba_blast_raw$V2 == "unplaced_scaffold229",]
bmpr1ba_gtf <- cluhar_v2.0.2_gtf[cluhar_v2.0.2_gtf$gene_id == "ENSCHAG00000012209"]


bmpr1ba_clustalo <- align_target_region(seqname_vec = bmpr1ba_seqnames, start_vec = bmpr1ba_starts, stop_vec = bmpr1ba_ends, dir_vec = c(1,1))
checkAlignment(bmpr1ba_clustalo)
bmpr1ba_avg_diff <- aln_diff_plot(bmpr1ba_clustalo, count_gaps = T, win_size = 200)

bmpr1ba_full_clustalo <- align_target_region(seqname_vec = bmpr1ba_seqnames, start_vec = bmpr1ba_full_starts, stop_vec = bmpr1ba_full_ends, dir_vec = c(1,1))
checkAlignment(bmpr1ba_full_clustalo)
bmpr1ba_full_avg_diff <- aln_diff_plot(bmpr1ba_full_clustalo, count_gaps = T, win_size = 200)

aln_pos_chr7 <- which(as.character(bmpr1ba_full_clustalo[2,]) != "-")
aln_pos_chr7_clustalo <- bmpr1ba_full_clustalo[,aln_pos_chr7]
#aln_pos_chr21_avg_diff <- aln_diff_plot(aln_pos_chr21_clustalo, count_gaps = F, win_size = 50)
aln_pos_chr7_gaps <- as.numeric(as.character(aln_pos_chr7_clustalo[1,]) == "-")
aln_pos_chr7_avg_gaps <- filter(aln_pos_chr7_gaps, rep(1/50, 50))
bmpr1ba_double_aln_pos <- which(as.character(aln_pos_chr7_clustalo[1,]) != "-")
bmpr1ba_double_aln_avg_diff <- aln_diff_plot(aln_pos_chr7_clustalo[,bmpr1ba_double_aln_pos], count_gaps = F, win_size = 50)

pdf(file = "~/Projects/Herring/doc/SexDet/bmpr1ba_full.pdf", width = 10, height = 7)
plot_gtf(bmpr1ba_gtf)
points(x= bmpr1ba_double_aln_pos + 12276563, y = bmpr1ba_double_aln_avg_diff - 0.25, pch = 16, cex = 0.5)
segments(x0 = 12276563 + bmpr1ba_blast[,7], x1 = 12276563 + bmpr1ba_blast[,8], y0 = 9, lwd = 3, col = "blue")
dev.off()

pdf(file = "~/Projects/Herring/doc/SexDet/bmpr1ba_12.276_12.287_Mb.pdf", width = 10, height = 7)
bmpr1ba_gtf_subset <- c(1,which(bmpr1ba_gtf$transcript_id == "ENSCHAT00000028135"))
plot_gtf(bmpr1ba_gtf[bmpr1ba_gtf_subset], x_lim = c(12276000, 12287000), y_lim = c(-0.5, 1.5))
abline(v = bmpr1ba_gtf[bmpr1ba_gtf_subset]@ranges@start, lwd  =.5, col = "grey30")
abline(v = bmpr1ba_gtf[bmpr1ba_gtf_subset]@ranges@start +  bmpr1ba_gtf[bmpr1ba_gtf_subset]@ranges@width, lwd  =.5, col = "grey30")
abline(h = c(0.75, 0.25, -0.25), col = "grey70")
abline(h = c(0.5, 0), col = "grey70", lty = "dashed")
segments(x0 = 12276563 + bmpr1ba_blast[,7], x1 = 12276563 + bmpr1ba_blast[,8], y0 = 1.5, lwd = 3, col = "blue")
#points(x= (1:length(aln_pos_chr21)) + 17016119, y = aln_pos_chr21_avg_diff/2, pch = 16, cex = 0.5)
points(x= bmpr1ba_double_aln_pos + 12276563, y = (bmpr1ba_double_aln_avg_diff/2)+0.25, pch = 16, cex = 0.5, col = "darkorchid")
points(x= (1:length(aln_pos_chr7)) + 12276563, y = -(aln_pos_chr7_avg_gaps/2) + .25, pch = 16, cex = 0.5, col = "olivedrab")
dev.off()

pdf(file = "~/Projects/Herring/doc/SexDet/bmpr1ba_diff.pdf", width = 10, height = 7)
plot(x= bmpr1ba_double_aln_pos  + 12276563, y = bmpr1ba_double_aln_avg_diff, xlim = c(12276000, 12287000), ylim = c(-1,1), type = "n")
gtf_tmp <- bmpr1ba_gtf[1:length(bmpr1ba_gtf) %in% bmpr1ba_gtf_subset & bmpr1ba_gtf$type == "CDS"]
rect(xleft = gtf_tmp@ranges@start, xright = gtf_tmp@ranges@start + gtf_tmp@ranges@width, ybottom = -1.5, ytop = 1.5, border = NA, lwd  =.5, col = "grey80")
abline(h = 0, col = "grey70")
abline(h = seq(from = 0.1, to = 0.5, by = 0.1), col = "grey80")
points(x= (1:length(aln_pos_chr7)) + 12276563, y = -aln_pos_chr7_avg_gaps, pch = 16, cex = 0.5, col = "olivedrab")
points(x= bmpr1ba_double_aln_pos  + 12276563, y = bmpr1ba_double_aln_avg_diff, pch = 16, cex = 0.5, col = "darkorchid")
dev.off()

#Ditsribution of differences
#diff_vec <- seg.sites(reg2_3_clustalo)
#aln_char_1 <- as.character(reg2_3_clustalo[1,])
#aln_char_2 <- as.character(reg2_3_clustalo[2,])
#diff_vec <- aln_char_1 != aln_char_2 & aln_char_1 != "-" & aln_char_2 != "-"
#plot(filter(as.numeric(diff_vec), rep(1/1000, 1000)))


#Some synteny examination
load(file = "~/Projects/Herring/data/HiC_assemblies/Version_2_release/fish_synteny_blocks_1Mb.RData")

#zebra_synteny_blocks <- satsuma_synteny_blocks("~/Projects/Herring/data/other_fish_genomes/Ch_2.0.2_v_Zebrafish_satsuma.out", "~/Projects/Herring/doc/HiC_satsuma/other_fish/Zebra_v_Ch_2.0.2_synteny.pdf", name_clean_re =".+Danio_rerio_(chromosome_[0-9]+).+")
#zebra_synteny_blocks_GR <-  synteny_block_ranges(zebra_synteny_blocks, "Zebrafish")
#tmp_synteny_GR <- zebra_synteny_blocks_GR[zebra_synteny_blocks_GR@seqnames == "chr8"]
#names(tmp_synteny_GR) <- paste(names(tmp_synteny_GR), 1:length(tmp_synteny_GR), sep = "_")
#as.data.frame(tmp_synteny_GR)

#medaka_synteny_blocks <- satsuma_synteny_blocks("/Users/mapet205/Projects/Herring/data/other_fish_genomes/Ch_2.0.2_v_Medaka_satsuma.out", "~/Projects/Herring/doc/HiC_satsuma/other_fish/Medaka_v_Ch_v2.0.2_synteny.pdf", name_clean_re = ".+(chromosome_[0-9]+)")
#medaka_synteny_blocks_GR <-  synteny_block_ranges(medaka_synteny_blocks, "Medaka")
plot_synteny_GR(medaka_synteny_blocks_GR)
plot_synteny_GR(pike_synteny_blocks_GR)
plot_synteny_GR(guppy_synteny_blocks_GR)
plot_synteny_GR(stickleback_synteny_blocks_GR)

#CATSPER3 hits - 
#unplaced_scaffold18:159671-160236
#unplaced_scaffold229:35348-35932
#chr8:22259379-22268294
CATSPER3_autosome_gtf <- rtracklayer::import("~/Projects/Herring/data/SexDetermination/CATSPER3_autosome_chr8_coordinates.gff")
CATSPER3_autosome_gtf$type <- "CDS"
seqinfo(CATSPER3_autosome_gtf,1:1) <- Seqinfo("8")
cs_reg1_seqnames <- c("unplaced_scaffold18", "chr8")
cs_reg1_starts <- c(154000, 22254000)
cs_reg1_ends <- c(170000, 22269000)
cs_reg1_clustalo <- align_target_region(seqname_vec = cs_reg1_seqnames, start_vec = cs_reg1_starts, stop_vec = cs_reg1_ends, dir_vec = c(1,1))
checkAlignment(cs_reg1_clustalo)
#dist.dna(cs_reg1_clustalo)
cs_reg1_chr <- aln_positions(cs_reg1_clustalo)
#cs_reg1_avg_diff <- aln_diff_plot_2(cs_reg1_clustalo, win_size = 50)
cs_reg1_avg_diff <- aln_diff_plot_2(cs_reg1_chr$clustalo, win_size = 50)
#mean(cs_reg1_avg_diff$avg_diff[980:4050])
#[1] 0.02330837
cs_reg1_exon <- exon_divergence(cs_reg1_avg_diff$aln_pos, exon_gtf = CATSPER3_autosome_gtf, chr_offset = 22254000, clustalo_obj = cs_reg1_chr$clustalo)



#plot(x= cs_reg1_avg_diff$aln_pos  + 22254000, y = cs_reg1_avg_diff$avg_diff, type  = "p", pch = 16, cex = 0.5, col = "darkorchid")


cs_reg2_seqnames <- c("unplaced_scaffold229", "chr8")
cs_reg2_starts <- c(30000, 22255000)
cs_reg2_ends <- c(50000, 22269000)
cs_reg2_clustalo <- align_target_region(seqname_vec = cs_reg2_seqnames, start_vec = cs_reg2_starts, stop_vec = cs_reg2_ends, dir_vec = c(-1,1))
checkAlignment(cs_reg2_clustalo)
cs_reg2_chr <- aln_positions(cs_reg2_clustalo)
#dist.dna(cs_reg2_clustalo)
#cs_reg2_avg_diff <- aln_diff_plot_2(cs_reg2_clustalo, win_size = 50)
cs_reg2_avg_diff <- aln_diff_plot_2(cs_reg2_chr$clustalo, win_size = 50)
#mean(cs_reg2_avg_diff$avg_diff[5620:5900]) #Homologous region
#[1] 0.03224199
cs_reg2_exon <- exon_divergence(cs_reg2_avg_diff$aln_pos, exon_gtf = CATSPER3_autosome_gtf, chr_offset = 22255000, clustalo_obj = cs_reg2_chr$clustalo)



#plot(x= cs_reg2_avg_diff$aln_pos  + 22258000, y = cs_reg2_avg_diff$avg_diff, type  = "p", pch = 16, cex = 0.5, col = "darkorchid")
plot(x= cs_reg1_avg_diff$aln_pos  + 22254000, y = cs_reg1_avg_diff$avg_diff, type  = "p", pch = 16, cex = 0.5, col = "darkorchid1")
points(x= cs_reg2_avg_diff$aln_pos  + 22255000, y = cs_reg2_avg_diff$avg_diff, type  = "p", pch = 16, cex = 0.5, col = "darkorchid4")
segments(x0 = CATSPER3_autosome_gtf@ranges@start, y0 = 1.2, y1 = -0.25, lwd  =.5, col = "grey30")
segments(x0 = CATSPER3_autosome_gtf@ranges@start +  CATSPER3_autosome_gtf@ranges@width, y0 = 1.2, y1 = -0.25, lwd  =.5, col = "grey30")

CATSPER3_autosome_gene <- GRanges(seqnames = 8, ranges = IRanges(start = 22259381, end = 22268494))
CATSPER3_autosome_gene@elementMetadata <- CATSPER3_autosome_gtf@elementMetadata[1,]
CATSPER3_autosome_gene$type <- "gene"
CATSPER3_autosome_transcript <- CATSPER3_autosome_gene
CATSPER3_autosome_transcript$type <- "transcript"
CATSPER3_gene_gtf <- c(CATSPER3_autosome_gene, CATSPER3_autosome_transcript,  CATSPER3_autosome_gtf)
CATSPER3_gene_gtf$gene_id <- as.character(CATSPER3_gene_gtf$Parent)
CATSPER3_gene_gtf$transcript_id <- "transcript_1"
CATSPER3_gene_gtf$transcript_id[1] <- NA


#Using extended scaffolds from earlier FALCON-unzip  assembly
unzip_000058F_026_016 <- readDNAStringSet("/Users/mapet205/Projects/Herring/data/SexDetermination/000058F_026_016.fa")
unzip_000058F_026_016 <- c(unzip_000058F_026_016, Ch_v2.0.2["chr8"])
unzip_reg1_seqnames <- c("000058F_026|arrow", "chr8")
unzip_reg1_starts <- c(70000, 22250000)
unzip_reg1_ends <- c(100000, 22275000)
unzip_reg1_clustalo <- align_target_region(seqname_vec = unzip_reg1_seqnames, start_vec = unzip_reg1_starts, stop_vec = unzip_reg1_ends, dir_vec = c(1,1), genome_set = unzip_000058F_026_016)
checkAlignment(unzip_reg1_clustalo)
#dist.dna(unzip_reg1_clustalo)
unzip_reg1_chr <- aln_positions(unzip_reg1_clustalo)
unzip_reg1_avg_diff <- aln_diff_plot_2(unzip_reg1_chr$clustalo, win_size = 50)
unzip_reg1_exon <- exon_divergence(unzip_reg1_avg_diff$aln_pos, exon_gtf = CATSPER3_autosome_gtf, chr_offset = 22250000, clustalo_obj = unzip_reg1_chr$clustalo)
plot_CDS_aln(2,full_clustalo = unzip_reg1_clustalo, chr_aln_pos = unzip_reg1_chr$chr_pos , CDS_gtf = CATSPER3_autosome_gtf, chr_offset = 22250000)

plot(x= unzip_reg1_avg_diff$aln_pos  + 22250000, y = unzip_reg1_avg_diff$avg_diff, type  = "p", pch = 16, cex = 0.5, col = "darkorchid1")
segments(x0 = CATSPER3_autosome_gtf@ranges@start, y0 = 1.2, y1 = -0.25, lwd  =.5, col = "grey30")
segments(x0 = CATSPER3_autosome_gtf@ranges@start +  CATSPER3_autosome_gtf@ranges@width, y0 = 1.2, y1 = -0.25, lwd  =.5, col = "grey30")

pdf(file = "~/Projects/Herring/doc/SexDet/Unzip_v_CATSPER3_annot.pdf", width = 10, height = 7)
plot_gtf(CATSPER3_gene_gtf, x_lim = c(22250000, 22270000), y_lim = c(-0.5, 1.5))
segments(x0 = CATSPER3_autosome_gtf@ranges@start, y0 = 1.2, y1 = -0.25, lwd  =.5, col = "grey30")
segments(x0 = CATSPER3_autosome_gtf@ranges@start +  CATSPER3_autosome_gtf@ranges@width, y0 = 1.2, y1 = -0.25, lwd  =.5, col = "grey30")
abline(h = c(0.75, 0.25, -0.25), col = "grey70")
abline(h = c(0.5, 0), col = "grey70", lty = "dashed")
#segments(x0 = 17016119 + bmpr1bb_blast[,7], x1 = 17016119 + bmpr1bb_blast[,8], y0 = 1.5, lwd = 3, col = "blue")
points(x= unzip_reg1_avg_diff$aln_pos  + 22250000, y = (unzip_reg1_avg_diff$avg_diff/2)+0.25, pch = 16, cex = 0.5, col = "darkorchid")
points(x= (1:length(unzip_reg1_avg_diff$avg_gaps)) + 22250000, y = -(unzip_reg1_avg_diff$avg_gaps/2) + .25, pch = 16, cex = 0.5, col = "olivedrab")
#cds_subset <- bmpr1bb_gtf$transcript_id == "ENSCHAT00000065423" & bmpr1bb_gtf$type == "CDS"
text(x = CATSPER3_autosome_gtf@ranges@start + CATSPER3_autosome_gtf@ranges@width/2, y = rep(c(1.3, 1.5), length.out = length(CATSPER3_autosome_gtf@ranges@start)), labels = paste("Exon", 1:length(CATSPER3_autosome_gtf@ranges@start), "\nd = ", round(unzip_reg1_exon[, "mean_diff"], digits = 2), sep = ""), cex = 1.1)

dev.off()



#Compiled Y-segment
rev_comp_chrY_blast_raw <- read.table("~/Projects/Herring/data/SexDetermination/rev_comp_chrY_v1.0.blastout", stringsAsFactors=F, header = F)
rev_comp_chrY_blast_chr8 <- rev_comp_chrY_blast_raw[rev_comp_chrY_blast_raw[,2] == "chr8", ]
rev_comp_chrY_blast_chr8_CATSPER <- rev_comp_chrY_blast_chr8[pmin(rev_comp_chrY_blast_chr8[,9], rev_comp_chrY_blast_chr8[,10]) < 22270000 &  pmax(rev_comp_chrY_blast_chr8[,9], rev_comp_chrY_blast_chr8[,10]) > 22250000 & rev_comp_chrY_blast_chr8[,4] > 1,]

pdf(file = "~/Projects/Herring/doc/SexDet/rev_Y_v_CATSPER3_blast.pdf", width = 10, height = 7)
plot(x = 0, y = 0, ylim = c(4.2e5,4.35e5), xlim = c(22250000, 22270000), type = "n")
cr <- colorRamp(colors = c("steelblue1", "firebrick1"))
col_vec <- (rev_comp_chrY_blast_chr8_CATSPER[,3] - 80)/20
col_vec[col_vec < 0] <- 0
segments(x0 = rev_comp_chrY_blast_chr8_CATSPER[,9], x1 = rev_comp_chrY_blast_chr8_CATSPER[,10], y0 = rev_comp_chrY_blast_chr8_CATSPER[,7], rev_comp_chrY_blast_chr8_CATSPER[,8], col = rgb(cr(col_vec), max = 255), lwd = 4)
dev.off()


#Mapping individual short reads to rev_comp_chrY
#Sorting read file paths
R1_files <- read.table("~/Projects/Herring/data/SexDetermination/short_read_mapping/R1_files.txt", stringsAsFactors = F)
R1_files[,"ID"] <- sub(".+/([A-Za-z0-9_]+)_R[12].fastq.gz", "\\1", R1_files[,1])
R2_files <- read.table("~/Projects/Herring/data/SexDetermination/short_read_mapping/R2_files.txt", stringsAsFactors = F)
R2_files[,"ID"] <- sub(".+/([A-Za-z0-9_]+)_R[12].fastq.gz", "\\1", R2_files[,1])
all(R2_files[,"ID"] == R1_files[,"ID"])
R2_sorted_files <- R2_files[match(R1_files[,"ID"], R2_files[,"ID"]),]
all(R2_sorted_files[,"ID"] == R1_files[,"ID"])
R1_files[,"Short_ID"] <- sub("([A-Za-z0-9]+)_.+", "\\1", R1_files[,"ID"])

write.table(x= R2_sorted_files[,1], file = "~/Projects/Herring/data/SexDetermination/short_read_mapping/R2_sorted_files.txt", row.names = F, quote = F, col.names = F)
write.table(x= R1_files[,1], file = "~/Projects/Herring/data/SexDetermination/short_read_mapping/R1_sorted_files.txt", row.names = F, quote = F, col.names = F)
write.table(x= R2_sorted_files[,2], file = "~/Projects/Herring/data/SexDetermination/short_read_mapping/sorted_IDs.txt", row.names = F, quote = F, col.names = F)
#To check RG assignment
#samtools view -H sample.bam | grep "@RG"

#The above yielded excess depth in all individuals, trying including the chromosomes for competition and using "-M" in bwa to supress freactional mappings
#Also mapping only the individuals to be published
SDR_target_inds <- read.table("~/Projects/Herring/data/SexDetermination/short_read_mapping/SDR_target_ind.txt", stringsAsFactors = F)
R1_target_files <- R1_files[R1_files$Short_ID %in% SDR_target_inds[,1],]
R2_target_files <- R2_sorted_files[R1_files$Short_ID %in% SDR_target_inds[,1],]
all(R2_target_files[,"ID"] == R1_target_files[,"ID"])
write.table(x= R2_target_files[,1], file = "~/Projects/Herring/data/SexDetermination/short_read_mapping/R2_target_files.txt", row.names = F, quote = F, col.names = F)
write.table(x= R1_target_files[,1], file = "~/Projects/Herring/data/SexDetermination/short_read_mapping/R1_target_files.txt", row.names = F, quote = F, col.names = F)
write.table(x= R2_target_files[,2], file = "~/Projects/Herring/data/SexDetermination/short_read_mapping/target_IDs.txt", row.names = F, quote = F, col.names = F)


rev_comp_chrY <- readDNAStringSet("~/Projects/Herring/data/SexDetermination/rev_comp_chrY_v1.0.fa")
writeXStringSet(c(Ch_v2.0.2[1:26], rev_comp_chrY),file = "~/Projects/Herring/data/SexDetermination/Ch2.0.2_CHR_and_Y.fa.gz", compress = T)
###Worked better.

##Aligning Chr Y to Chr8 to estiamte divergance across the segments
##To be edited

Autosomes_and_Y <- c(Ch_v2.0.2[1:26], rev_comp_chrY)
chrY_seqnames <- c("rc_chrY", "chr8")
#First block
chrY_reg1_starts <- c(1, 20970001)
chrY_reg1_ends <- c(44000, 21016000)
chrY_reg1_clustalo <- align_target_region(seqname_vec = chrY_seqnames, start_vec = chrY_reg1_starts, stop_vec = chrY_reg1_ends, dir_vec = c(1,1), genome_set = Autosomes_and_Y )
checkAlignment(chrY_reg1_clustalo)
chrY_reg1_adj <- c(1000, 9100)
chrY_reg1_adj_starts <- chrY_reg1_starts + chrY_reg1_adj[1]
chrY_reg1_clustalo_adj <- chrY_reg1_clustalo[,chrY_reg1_adj[1]:(dim(chrY_reg1_clustalo)[2]-chrY_reg1_adj[2])]

#dist.dna(unzip_reg1_clustalo)
chrY_reg1_chr <- aln_positions(chrY_reg1_clustalo_adj)
chrY_reg1_avg_diff <- aln_diff_plot_2(chrY_reg1_chr$clustalo, win_size = 50)
chrY_reg1_avg_diff_500 <- aln_diff_plot_2(chrY_reg1_chr$clustalo, win_size = 500)
reg1_pos_X <- (1:length(chrY_reg1_chr$chr_pos))[chrY_reg1_chr$double_pos] + chrY_reg1_starts[2]
reg1_pos_Y = (1:length(chrY_reg1_chr$opposing_pos))[chrY_reg1_chr$opposing_double_pos] + chrY_reg1_starts[1]
chrY_reg1_diff_df <- data.frame(diff_50 = chrY_reg1_avg_diff$avg_diff, diff_500 = chrY_reg1_avg_diff_500$avg_diff, pos_X = reg1_pos_X, pos_Y = reg1_pos_Y, stringsAsFactors = F)
#Average divergence over the block
mean(chrY_reg1_diff_df$diff_50, na.rm = T) * 100
#[1] 0.6142888

##
#Second block, reverse order
chrY_reg2_starts <- c(51001, 21016001)
chrY_reg2_ends <- c(85000, 21050000)
chrY_reg2_clustalo <- align_target_region(seqname_vec = chrY_seqnames, start_vec = chrY_reg2_starts, stop_vec = chrY_reg2_ends, dir_vec = c(-1,1), genome_set = Autosomes_and_Y )
checkAlignment(chrY_reg2_clustalo)
chrY_reg2_chr <- aln_positions(chrY_reg2_clustalo)
chrY_reg2_avg_diff <- aln_diff_plot_2(chrY_reg2_chr$clustalo, win_size = 50)
chrY_reg2_avg_diff_500 <- aln_diff_plot_2(chrY_reg2_chr$clustalo, win_size = 500)
reg2_pos_X <- (1:length(chrY_reg2_chr$chr_pos))[chrY_reg2_chr$double_pos] + chrY_reg2_starts[2]
reg2_pos_Y = (length(chrY_reg2_chr$opposing_pos) - (1:length(chrY_reg2_chr$opposing_pos))[chrY_reg2_chr$opposing_double_pos]) + chrY_reg2_starts[1]

chrY_reg2_diff_df <- data.frame(diff_50 = chrY_reg2_avg_diff$avg_diff, diff_500 = chrY_reg2_avg_diff_500$avg_diff, pos_X = reg2_pos_X, pos_Y = reg2_pos_Y, stringsAsFactors = F)
#Average divergence over the block
mean(chrY_reg2_diff_df$diff_50, na.rm = T) * 100
#[1] 0.08377465

#Third block,
chrY_reg3_starts <- c(87001, 21051001)
chrY_reg3_ends <- c(99000, 21063500)
chrY_reg3_clustalo <- align_target_region(seqname_vec = chrY_seqnames, start_vec = chrY_reg3_starts, stop_vec = chrY_reg3_ends, dir_vec = c(1,1), genome_set = Autosomes_and_Y )
checkAlignment(chrY_reg3_clustalo)
chrY_reg3_chr <- aln_positions(chrY_reg3_clustalo)
chrY_reg3_avg_diff <- aln_diff_plot_2(chrY_reg3_chr$clustalo, win_size = 50)
chrY_reg3_avg_diff_500 <- aln_diff_plot_2(chrY_reg3_chr$clustalo, win_size = 500)

reg3_pos_X <- (1:length(chrY_reg3_chr$chr_pos))[chrY_reg3_chr$double_pos] + chrY_reg3_starts[2]
reg3_pos_Y = (1:length(chrY_reg3_chr$opposing_pos))[chrY_reg3_chr$opposing_double_pos] + chrY_reg3_starts[1]
chrY_reg3_diff_df <- data.frame(diff_50 = chrY_reg3_avg_diff$avg_diff, diff_500 = chrY_reg3_avg_diff_500$avg_diff, pos_X = reg3_pos_X, pos_Y = reg3_pos_Y, stringsAsFactors = F)
#Average divergence over the block
mean(chrY_reg3_diff_df$diff_50, na.rm = T) * 100
#[1] 0.3651185
#First three blocks
mean(c(chrY_reg3_diff_df$diff_50, chrY_reg2_diff_df$diff_50, chrY_reg1_diff_df$diff_50), na.rm = T) * 100
#[1] 0.3768512


#Fourth block, internal PAR
chrY_reg4_starts <- c(220000, 21063000)
chrY_reg4_ends <- c(340000, 21183000)
chrY_reg4_clustalo <- align_target_region(seqname_vec = chrY_seqnames, start_vec = chrY_reg4_starts, stop_vec = chrY_reg4_ends, dir_vec = c(1,1), genome_set = Autosomes_and_Y )
checkAlignment(chrY_reg4_clustalo)
image(chrY_reg4_clustalo[,1:15000])
chrY_reg4_adj <- c(10500, 4000)
chrY_reg4_adj_starts <- chrY_reg4_starts + chrY_reg4_adj[1]
chrY_reg4_clustalo_adj <- chrY_reg4_clustalo[,chrY_reg4_adj[1]:(dim(chrY_reg4_clustalo)[2]-chrY_reg4_adj[2])]
image(chrY_reg4_clustalo_adj)

chrY_reg4_chr <- aln_positions(chrY_reg4_clustalo_adj)
chrY_reg4_avg_diff <- aln_diff_plot_2(chrY_reg4_chr$clustalo, win_size = 50)
chrY_reg4_avg_diff_500 <- aln_diff_plot_2(chrY_reg4_chr$clustalo, win_size = 500)
reg4_pos_X <- (1:length(chrY_reg4_chr$chr_pos))[chrY_reg4_chr$double_pos] + chrY_reg4_adj_starts[2]
reg4_pos_Y = (1:length(chrY_reg4_chr$opposing_pos))[chrY_reg4_chr$opposing_double_pos] + sum(which(as.character(chrY_reg4_clustalo[1,]) != "-") < chrY_reg4_adj[1]) + chrY_reg4_starts[1]
chrY_reg4_diff_df <- data.frame(diff_50 = chrY_reg4_avg_diff$avg_diff, diff_500 = chrY_reg4_avg_diff_500$avg_diff, pos_X = reg4_pos_X, pos_Y = reg4_pos_Y, stringsAsFactors = F)
#Average divergence over the block
mean(chrY_reg4_diff_df$diff_50, na.rm = T) * 100
#[1] 0.9594129

#Fifth block, distal PAR
chrY_reg5_starts <- c(438000, 21178000)
chrY_reg5_ends <- c(511498, 21252000)
chrY_reg5_clustalo <- align_target_region(seqname_vec = chrY_seqnames, start_vec = chrY_reg5_starts, stop_vec = chrY_reg5_ends, dir_vec = c(1,1), genome_set = Autosomes_and_Y )
image(chrY_reg5_clustalo)
chrY_reg5_adj <- c(3000, 3000)
chrY_reg5_clustalo_adj <- chrY_reg5_clustalo[,chrY_reg5_adj[1]:(dim(chrY_reg5_clustalo)[2]-chrY_reg5_adj[2])]
chrY_reg5_adj_starts <- chrY_reg5_starts + chrY_reg5_adj[1]
image(chrY_reg5_clustalo_adj)
checkAlignment(chrY_reg5_clustalo_adj)
chrY_reg5_chr <- aln_positions(chrY_reg5_clustalo_adj)
chrY_reg5_avg_diff <- aln_diff_plot_2(chrY_reg5_chr$clustalo, win_size = 50)
chrY_reg5_avg_diff_500 <- aln_diff_plot_2(chrY_reg5_chr$clustalo, win_size = 500)
reg5_pos_X <- (1:length(chrY_reg5_chr$chr_pos))[chrY_reg5_chr$double_pos] + chrY_reg5_adj_starts[2]
reg5_pos_Y = (1:length(chrY_reg5_chr$opposing_pos))[chrY_reg5_chr$opposing_double_pos] + sum(which(as.character(chrY_reg5_clustalo[1,]) != "-") < chrY_reg5_adj[1]) + chrY_reg5_starts[1]
chrY_reg5_diff_df <- data.frame(diff_50 = chrY_reg5_avg_diff$avg_diff, diff_500 = chrY_reg5_avg_diff_500$avg_diff, pos_X = reg5_pos_X, pos_Y = reg5_pos_Y, stringsAsFactors = F)
#Average divergence over the block
mean(chrY_reg5_diff_df$diff_50, na.rm = T) * 100
#[1] 0.5596298


#Combining the blocks to one plottable values series
divergence_plot_df <- rbind(chrY_reg1_diff_df[1:dim(chrY_reg1_diff_df)[1],], chrY_reg2_diff_df[1:dim(chrY_reg2_diff_df)[1],], chrY_reg3_diff_df[1:dim(chrY_reg3_diff_df)[1],], chrY_reg4_diff_df[1:dim(chrY_reg4_diff_df)[1],], chrY_reg5_diff_df[1:dim(chrY_reg5_diff_df)[1],])
pdf(file = "~/Projects/Herring/doc/SexDet/rev_Y_divergence.pdf", width = 20, height = 4)
plot(x= divergence_plot_df$pos_Y, y = divergence_plot_df$diff_500, ylim = c(0,0.15), type = "p", cex = 0.3, pch = 16, col = "darkorchid", axes = F)
axis(1)
axis(4, at = c(0, 0.05, 0.1, 0.15), pos = 5.2e5)
dev.off()

#Fst along Chr 8 between males and females
#List of males, estimated from haplotyes on scaffold149 from genome v1.2 ("scaffold149_900_to_1000_kb.pdf")
male_ind <- grep("AM[0-9]{1,2}|BM[0-9]{1,2}|Gavle100|Fehmarn6|S4TH244|S4TH231|S4TM211|S3Ps259", herring_79$sample_list, value = T)
male_ind <- unique(sub("_[12]", "", male_ind))
female_ind <- grep("Baltic|Atlantic", herring_79$sample_list, value = T)
female_ind <- unique(sub("_[12]","", female_ind))
female_ind <- female_ind[!grepl("Downs|Celticsea|IsleofMan|Balsfjord|NSSH|CelticSea", female_ind)]
female_ind <- female_ind[!(female_ind %in% male_ind)]
write(female_ind, file = "~/Projects/Herring/data/SexDetermination/h79_est_females.txt")
write(male_ind, file = "~/Projects/Herring/data/SexDetermination/h79_est_males.txt")
#Fst estimated using vcftools, 10k windows, 5k step
#FvM_fst <- read.table("~/Projects/Herring/data/SexDetermination/Chr8_females_v_males.windowed.weir.fst", header = T, stringsAsFactors = F)
pdf(file = "~/Projects/Herring/doc/SexDet/FvM_Fst_10k.pdf", width = 12, height = 6)
plot(x = FvM_fst$BIN_START[FvM_fst$N_VARIANTS >= 100], y = FvM_fst$WEIGHTED_FST[FvM_fst$N_VARIANTS >= 100], col = "darkorchid", pch = 20, cex = 0.8, xlab = "Position", ylab = "Weighted Fst")
abline(v = 21178000)
plot(x = FvM_fst$BIN_START[FvM_fst$N_VARIANTS >= 100], y = FvM_fst$WEIGHTED_FST[FvM_fst$N_VARIANTS >= 100], col = "darkorchid", pch = 20, cex = 0.8, xlab = "Position", ylab = "Weighted Fst", xlim = c(20e6, 23e6))
abline(v = 21178000)
dev.off()
#Fst estimated using vcftools, 50k windows, 25k step
FvM_fst <- read.table("~/Projects/Herring/data/SexDetermination/Chr8_females_v_males_50k.windowed.weir.fst", header = T, stringsAsFactors = F)
pdf(file = "~/Projects/Herring/doc/SexDet/FvM_Fst_50k.pdf", width = 15, height = 4)
plot(x = FvM_fst$BIN_START[FvM_fst$N_VARIANTS >= 500], y = FvM_fst$WEIGHTED_FST[FvM_fst$N_VARIANTS >= 500], col = "darkorchid", pch = 20, cex = 0.8, xlab = "Position", ylab = "Weighted Fst")
abline(v = 21178000)
plot(x = FvM_fst$BIN_START[FvM_fst$N_VARIANTS >= 500], y = FvM_fst$WEIGHTED_FST[FvM_fst$N_VARIANTS >= 500], col = "darkorchid", pch = 20, cex = 0.8, xlab = "Position", ylab = "Weighted Fst", xlim = c(20e6, 23e6))
abline(v = 21178000)
dev.off()

#Support functions
exon_divergence <- function(aln_pos_vec, exon_gtf, chr_offset, chr_no = 8, clustalo_obj){
  double_aln_GR <- GRanges(seqnames = chr_no, ranges = IRanges(start = aln_pos_vec + chr_offset, end = aln_pos_vec + chr_offset))
  double_aln_in_CDS <- findOverlaps(double_aln_GR, exon_gtf, ignore.strand = T)
  
  n_exons <- length(unique(double_aln_in_CDS@to))
  exon_mean_diff <- data.frame(exon = 1:n_exons)
  diff_vec <- as.numeric(clustalo_obj[1,aln_pos_vec] != clustalo_obj[2,aln_pos_vec])
  gap_vec <- as.numeric(as.character(clustalo_obj[1,]) == "-")
  for(exon in 1:n_exons){
    exon_mean_diff[exon_mean_diff$exon == exon, "mean_diff"] <- mean(diff_vec[double_aln_in_CDS@from[double_aln_in_CDS@to == exon]])
    exon_chr_aln_start <- exon_gtf@ranges@start[exon] - chr_offset
    exon_chr_aln_end <- exon_gtf@ranges@start[exon] +  exon_gtf@ranges@width[exon]- chr_offset
    exon_mean_diff[exon_mean_diff$exon == exon, "gap_count"] <- sum(gap_vec[exon_chr_aln_start:exon_chr_aln_end]) #/bmpr1bb_gtf[cds_subset]@ranges@width[exon]
    #exon_mean_diff[exon_mean_diff$exon == exon, "insertion_count"] <- sum(as.character(bmpr1bb_full_clustalo[2,aln_pos_chr21[exon_chr21_aln_start]:aln_pos_chr21[exon_chr21_aln_end]]) == "-") #/bmpr1bb_gtf[cds_subset]@ranges@width[exon]
    #exon_mean_diff[exon_mean_diff$exon == exon, "exon_length_Y"] <- sum(as.character(bmpr1bb_full_clustalo[1,aln_pos_chr21[exon_chr21_aln_start]:aln_pos_chr21[exon_chr21_aln_end]]) != "-")
    #exon_mean_diff[exon_mean_diff$exon == exon, "exon_length_chr21"] <- sum(as.character(bmpr1bb_full_clustalo[2,aln_pos_chr21[exon_chr21_aln_start]:aln_pos_chr21[exon_chr21_aln_end]]) != "-")
    #if(exon < n_exons){
    #  intron_end <- aln_pos_chr21[bmpr1bb_gtf[cds_subset]@ranges@start[exon+1] - 1 - 17016119]
    #  intron_start <- aln_pos_chr21[bmpr1bb_gtf[cds_subset]@ranges@start[exon] + bmpr1bb_gtf[cds_subset]@ranges@width[exon] - 17016119]
    #  exon_mean_diff[exon_mean_diff$exon == exon, "next_intron_Y"] <- sum(as.character(bmpr1bb_full_clustalo[1,intron_start:intron_end]) != "-")
    #  exon_mean_diff[exon_mean_diff$exon == exon, "next_intron_chr21"] <- sum(as.character(bmpr1bb_full_clustalo[2,intron_start:intron_end]) != "-")
    #}
  }
  return(exon_mean_diff)
}

#Support functions
aln_positions <- function(clustalo_obj){
  aln_pos_chr <- which(as.character(clustalo_obj[2,]) != "-")
  aln_pos_opp <- which(as.character(clustalo_obj[1,]) != "-")
  #checkAlignment(bmpr1bb_full_clustalo[,aln_pos_chr21])
  aln_pos_chr_clustalo <- clustalo_obj[,aln_pos_chr]
  #aln_pos_chr21_avg_diff <- aln_diff_plot(aln_pos_chr21_clustalo, count_gaps = F, win_size = 50)
  #aln_pos_chr21_gaps <- as.numeric(as.character(aln_pos_chr21_clustalo[1,]) == "-")
  #aln_pos_chr21_avg_gaps <- filter(aln_pos_chr21_gaps, rep(1/50, 50))
  
  #points(x= (1:length(aln_pos_chr21)) + 17016119, y = bmpr1bb_full_avg_diff[aln_pos_chr21]/2, pch = 16, cex = 0.5)
  #plot(x= (1:length(aln_pos_chr21)) + 17016119, y = bmpr1bb_full_avg_diff[aln_pos_chr21], pch = 16, cex = 0.5)
  
  #checkAlignment(aln_pos_chr21_clustalo[,25000:length(aln_pos_chr21)])
  double_aln_pos <- which(as.character(aln_pos_chr_clustalo[1,]) != "-")
  opp_double_aln_pos <- which(as.character(clustalo_obj[2,aln_pos_opp]) != "-")
  #double_aln_avg_diff <- aln_diff_plot(aln_pos_chr21_clustalo[,double_aln_pos], count_gaps = F, win_size = 50)
  return(list(chr_pos = aln_pos_chr, clustalo = aln_pos_chr_clustalo, double_pos = double_aln_pos, opposing_pos = aln_pos_opp, opposing_double_pos = opp_double_aln_pos))
}


plot_synteny_GR <- function(synteny_GR, target_chr = "chr8", pos_limits = c(21e6, 21.6e6)){
  tmp_synteny_GR <- synteny_GR[synteny_GR@seqnames == target_chr]
  names(tmp_synteny_GR) <- paste(names(tmp_synteny_GR), 1:length(tmp_synteny_GR), sep = "_")
  #as.data.frame(tmp_synteny_GR)
  plot(x= 0, y = 0, xlim = pos_limits, type = "n", ylim = c(0,2), main = synteny_GR@elementMetadata$q_species[1])
  segments(x0 = tmp_synteny_GR@ranges@start, x1 = tmp_synteny_GR@ranges@start + tmp_synteny_GR@ranges@width, y0 =  1)
  #tmp_synteny_GR[tmp_synteny_GR@ranges@start < 21.3e6 & tmp_synteny_GR@ranges@start + tmp_synteny_GR@ranges@width > 21.0e6]
}

align_target_region <- function(seqname_vec, start_vec, stop_vec, genome_set = Ch_v2.0.2, dir_vec = c(-1, 1)){
  aln_target <- GRanges(seqnames = seqname_vec, ranges = IRanges(start = start_vec, end = stop_vec))
  aln_seq <- genome_set[aln_target]
  names(aln_seq) <- paste(aln_target)
  for(i in 1:length(dir_vec)){
    if(dir_vec[i] == -1) {aln_seq[i] <- reverseComplement(aln_seq[i])}
  }
  aln_seq_bin <- as.DNAbin(aln_seq)
  reg_clustalo <- clustalomega(aln_seq_bin, exec="clustalo", quiet = F)
  #checkAlignment(reg_clustalo)
  return(reg_clustalo)
}

aln_diff_plot_2 <- function(aln_obj, win_size = 1000, pdf_file = NULL, y_lim = c(-1,1), ...){
  aln_char_us <- as.character(aln_obj[1,])
  aln_char_chr <- as.character(aln_obj[2,])
  
  aln_pos <- which(aln_char_chr != "-" & aln_char_us != "-")
  gap_pos <- aln_char_chr == "-" | aln_char_us == "-"
  aln_pos_clustalo <- aln_obj[,aln_pos]
  #aln_pos_chr_gaps <- as.numeric(as.character(aln_pos_chr_clustalo[1,]) == "-")
  avg_gaps <- filter(as.numeric(gap_pos), rep(1/win_size, win_size))
  
  #double_aln_pos <- which(as.character(aln_pos_chr_clustalo[1,]) != "-")
  diff_vec <- as.character(aln_obj[1,aln_pos]) != as.character(aln_obj[2,aln_pos])
  aln_avg_diff <- filter(as.numeric(diff_vec), rep(1/win_size, win_size))
  
  #avg_diff_vec <- filter(as.numeric(diff_vec), rep(1/win_size, win_size))
  
  if(!is.null(pdf_file)) pdf(file = pdf_file, width = 10, height = 7)
  #if(is.null(y_lim)) y_lim <- c(-1,1)
  plot(x= aln_pos, y = aln_avg_diff, ylim = y_lim, type = "n")
  #gtf_tmp <- bmpr1ba_gtf[1:length(bmpr1ba_gtf) %in% bmpr1ba_gtf_subset & bmpr1ba_gtf$type == "CDS"]
  #rect(xleft = gtf_tmp@ranges@start, xright = gtf_tmp@ranges@start + gtf_tmp@ranges@width, ybottom = -1.5, ytop = 1.5, border = NA, lwd  =.5, col = "grey80")
  abline(h = 0, col = "grey70")
  abline(h = seq(from = 0.1, to = 0.5, by = 0.1), col = "grey80")
  points(x= 1:length(aln_obj[1,]) , y = -avg_gaps, pch = 16, cex = 0.5, col = "olivedrab")
  points(x= aln_pos, y = aln_avg_diff, pch = 16, cex = 0.5, col = "darkorchid")
  #plot(y = avg_diff_vec, x = 1:length(avg_diff_vec), xlab = "", ylab = "Difference (fraction)", ... )
  
  if(!is.null(pdf_file)) dev.off()
  return(invisible(list(avg_diff = aln_avg_diff, avg_gaps = avg_gaps, gap_pos = gap_pos, aln_pos = aln_pos)))
}

aln_diff_plot <- function(aln_obj, win_size = 1000, pdf_file = NULL, count_gaps = F, ...){
  aln_char_1 <- as.character(aln_obj[1,])
  aln_char_2 <- as.character(aln_obj[2,])
  if(!count_gaps) diff_vec <- aln_char_1 != aln_char_2 & aln_char_1 != "-" & aln_char_2 != "-"
  if(count_gaps) diff_vec <- aln_char_1 != aln_char_2
  avg_diff_vec <- filter(as.numeric(diff_vec), rep(1/win_size, win_size))
  if(!is.null(pdf_file)) pdf(file = pdf_file, width = 10, height = 7)
  
  plot(y = avg_diff_vec, x = 1:length(avg_diff_vec), xlab = "", ylab = "Difference (fraction)", ... )
  
  if(!is.null(pdf_file)) dev.off()
  return(invisible(avg_diff_vec))
}


plot_gtf <- function(gtf_gr, outfile = NULL, display_id  = F, x_lim = NULL, y_lim = NULL){
  if(!is.null(outfile)){
    pdf(file = outfile, width = 8, height = 4)
  }
  level_count <- max(rowSums(table(gtf_gr$gene_id, gtf_gr$transcript_id) > 0)) + 2
  if(is.null(x_lim)) x_lim <- c(min(start(gtf_gr)), max(end(gtf_gr)))
  if(is.null(y_lim)) y_lim <- c(0, level_count)
  plot(x = mean(start(gtf_gr)), y = level_count /2, xlim = x_lim, ylim = y_lim, type = "n", main = paste("Chr", unique(gtf_gr@seqnames)), xlab = "Position", ylab = "")
  genes <-  unique(gtf_gr$gene_id)
  for (gene in genes){
    current_lvl <- level_count
    gene_entry <- gtf_gr$type == "gene" & gtf_gr$gene_id == gene
    
    text(x = mean(c(start(gtf_gr)[gene_entry], end(gtf_gr)[gene_entry])), y = current_lvl - 0.6, labels = gene, cex = 0.8)
    current_lvl <- current_lvl - 1
    
    segments(x0 = start(gtf_gr)[gene_entry], x1 = end(gtf_gr)[gene_entry], y0 = current_lvl, col = "firebrick", lwd = 5)
    transcript_entries <- which(gtf_gr$type == "transcript" & gtf_gr$gene_id == gene)
    current_lvl <- current_lvl - 1
    
    for(transcript in transcript_entries){
      transcript_id <- gtf_gr$transcript_id[transcript]
      if(display_id){
        text(x = start(gtf_gr)[transcript], y = current_lvl + 0.4, labels = transcript_id, cex = 0.8, pos = 4)
      }
      segments(x0 = start(gtf_gr)[transcript], x1 = end(gtf_gr)[transcript], y0 = current_lvl, col = "black", lwd = 2)
      exon_entries <- gtf_gr$type == "CDS" & gtf_gr$transcript_id == transcript_id
      rect(xleft = start(gtf_gr)[exon_entries], xright = end(gtf_gr)[exon_entries ], ybottom = current_lvl - 0.2, ytop = current_lvl + 0.2, col = "black")
      if(any(gtf_gr$type == "five_prime_utr" & gtf_gr$transcript_id == transcript_id)){
        utr5_entries <- gtf_gr$type == "five_prime_utr" & gtf_gr$transcript_id == transcript_id
        rect(xleft = start(gtf_gr)[utr5_entries], xright = end(gtf_gr)[utr5_entries], ybottom = current_lvl - 0.2, ytop = current_lvl + 0.2, col = "firebrick1")
      }
      if(any(gtf_gr$type == "three_prime_utr" & gtf_gr$transcript_id == transcript_id)){
        utr3_entries <- gtf_gr$type == "three_prime_utr" & gtf_gr$transcript_id == transcript_id
        rect(xleft = start(gtf_gr)[utr3_entries], xright = end(gtf_gr)[utr3_entries], ybottom = current_lvl - 0.2, ytop = current_lvl + 0.2, col = "firebrick4")
      }
      
      current_lvl <- current_lvl - 1
    }
  }
  if(!is.null(outfile)){
    dev.off()
  }
}

plot_CDS_aln <- function(exon, full_clustalo, chr_aln_pos, CDS_gtf, chr_offset = 17016119){
  exon_aln_start <- CDS_gtf@ranges@start[exon] - chr_offset
  exon_aln_end <- CDS_gtf@ranges@start[exon] +  CDS_gtf@ranges@width[exon]-chr_offset
  image(full_clustalo[,chr_aln_pos[exon_aln_start]:chr_aln_pos[exon_aln_end]])
}

##########################################################################################
# Compare results from 'diffcyt-med', 'diffcyt-FDA-unwtd', 'diffcyt-FDA-wtd'
##########################################################################################


# load results from server
load("../../../RData/res_DE_FDA_wtd.RData")

# note: FlowSOM random seeds are not reproducible across operating systems! (need to run all methods on same operating system)


rowData(res_DE_med)
rowData(res_DE_FDA_unwtd)
rowData(res_DE_FDA_wtd)

sum(is.na(rowData(res_DE_FDA_unwtd)$p_adj))
sum(is.na(rowData(res_DE_FDA_wtd)$p_adj))

# remove NAs from res_DE_FDA_unwtd and res_DE_FDA_wtd
# (note: change this in the testing function: remove NAs before returning final SummarizedExperiment)
res_DE_FDA_unwtd_old <- res_DE_FDA_unwtd_old
res_DE_FDA_wtd_old <- res_DE_FDA_wtd_old

res_DE_FDA_unwtd <- res_DE_FDA_unwtd[!is.na(rowData(res_DE_FDA_unwtd)$p_adj), ]
res_DE_FDA_wtd <- res_DE_FDA_wtd[!is.na(rowData(res_DE_FDA_wtd)$p_adj), ]
dim(res_DE_FDA_unwtd)
dim(res_DE_FDA_wtd)

dim(rowData(res_DE_med))
dim(rowData(res_DE_FDA_unwtd))
dim(rowData(res_DE_FDA_wtd))

# get common cluster-marker combinations between the three methods
exists_DE_med <- paste(rowData(res_DE_med)$cluster, rowData(res_DE_med)$marker, sep = "_")
exists_DE_FDA_unwtd <- paste(rowData(res_DE_FDA_unwtd)$cluster, rowData(res_DE_FDA_unwtd)$marker, sep = "_")
exists_DE_FDA_wtd <- paste(rowData(res_DE_FDA_wtd)$cluster, rowData(res_DE_FDA_wtd)$marker, sep = "_")

length(exists_DE_med)
length(exists_DE_FDA_unwtd)
length(exists_DE_FDA_wtd)

# intersection of all methods
exists_all <- intersect(intersect(exists_DE_FDA_unwtd, exists_DE_FDA_wtd), exists_DE_med)
length(exists_all)

ix_DE_med <- match(exists_all, exists_DE_med)
ix_DE_FDA_unwtd <- match(exists_all, exists_DE_FDA_unwtd)
ix_DE_FDA_wtd <- match(exists_all, exists_DE_FDA_wtd)

length(ix_DE_med)
length(ix_DE_FDA_unwtd)
length(ix_DE_FDA_wtd)

# subset
res_DE_med_sub <- res_DE_med[ix_DE_med, ]
res_DE_FDA_unwtd_sub <- res_DE_FDA_unwtd[ix_DE_FDA_unwtd, ]
res_DE_FDA_wtd_sub <- res_DE_FDA_wtd[ix_DE_FDA_wtd, ]

rowData(res_DE_med_sub)
rowData(res_DE_FDA_unwtd_sub)
rowData(res_DE_FDA_wtd_sub)

# plots: correlation of adjusted p-values between methods
p_adj_med <- rowData(res_DE_med_sub)$adj.P.Val
p_adj_FDA_unwtd <- rowData(res_DE_FDA_unwtd_sub)$p_adj
p_adj_FDA_wtd <- rowData(res_DE_FDA_wtd_sub)$p_adj


plot(p_adj_med, p_adj_FDA_unwtd)
plot(p_adj_med, p_adj_FDA_unwtd, xlim = c(0, 0.1), ylim = c(0, 0.1))

plot(p_adj_med, p_adj_FDA_wtd)
plot(p_adj_med, p_adj_FDA_wtd, xlim = c(0, 0.1), ylim = c(0, 0.1))

plot(p_adj_FDA_unwtd, p_adj_FDA_wtd)


smoothScatter(p_adj_med, p_adj_FDA_unwtd)
smoothScatter(p_adj_med, p_adj_FDA_unwtd, xlim = c(0, 0.005), ylim = c(0, 0.2))

smoothScatter(p_adj_med, p_adj_FDA_wtd)
smoothScatter(p_adj_med, p_adj_FDA_wtd, xlim = c(0, 0.005), ylim = c(0, 0.2))

smoothScatter(p_adj_FDA_unwtd, p_adj_FDA_wtd)
smoothScatter(p_adj_FDA_unwtd, p_adj_FDA_wtd, xlim = c(0, 0.1), ylim = c(0, 0.1))


# tables: similarity between methods
table(p_adj_med < 0.05, p_adj_FDA_unwtd < 0.05)
table(p_adj_med < 0.01, p_adj_FDA_unwtd < 0.05)  # assuming 'diffcyt-med' has poor FDR control

table(p_adj_med < 0.05, p_adj_FDA_wtd < 0.05)
table(p_adj_med < 0.01, p_adj_FDA_wtd < 0.05)  # assuming 'diffcyt-med' has poor FDR control

table(p_adj_FDA_unwtd < 0.05, p_adj_FDA_wtd < 0.05)

table(p_adj_FDA_unwtd == 0, p_adj_FDA_wtd == 0)  # most highly significant




#########################################
# Script to transform data set 'BCR-XL' #
#########################################


# indices of lineage markers (for clustering) and functional markers (not used for
# clustering); see Table 1 in Bruggner et al. (2014)
marker_cols <- c(3:4, 7:9, 11:19, 21:22, 24:26, 28:31, 33)
lineage_cols <- c(3:4, 9, 11,12,14, 21, 29, 31, 33)
func_cols <- setdiff(marker_cols, lineage_cols)

length(marker_cols)   # should have 24 total markers
length(lineage_cols)  # 10 lineage markers
length(func_cols)     # 14 functional markers


# transform all marker columns; using asinh with cofactor = 5 (see Bendall et al. 2011, 
# Supp. Fig. S2)

data_BCRXL <- transformData(data_BCRXL, cofactor = 5, marker_cols = marker_cols)


# check markers: see Table 1, Bruggner et al. (2014)
colnames(data_BCRXL)[lineage_cols]
colnames(data_BCRXL)[func_cols]



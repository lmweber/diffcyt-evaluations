#######################################################
# Script to detect differential abundance of clusters #
#######################################################

# detect clusters with differential abundance (cell frequencies) between conditions

# using 'limma' empirical Bayes methods to share information on variance between 
# clusters, to stabilize mean-variance relationship and improve power (under the
# assumption that most clusters are not differentially abundant)


# run after scripts "load_and_transform_data.R" and "run_clustering.R"
#source("load_and_transform_data.R")
#source("run_clustering.R")


library(limma)


# use expression values only
data_exprs <- lapply(data_transf, function(u) {
  exprs(u)
})


# number of cells per cluster
n_cells <- sapply(data_exprs, nrow)

# cluster frequencies and proportions (per sample)
samp <- rep(names(data_exprs), n_cells)
freq <- table(cluster = clus, sample = samp)  ## note 'table' re-arranges sample names alphabetically
prop <- t(t(freq) / colSums(freq))
sqrt_prop <- sqrt(prop)

# check re-arranged sample names
all(colSums(freq)[names(n_cells)] == n_cells)


# group factor and model matrix
grp <- gsub("^PBMC8_30min_patient[1-8]_|\\.fcs$", "", basename(fs))
grp <- factor(grp, levels = sort(unique(grp), decreasing = TRUE))

mm <- model.matrix(~ grp)


# fit linear models for each cluster
fit <- lmFit(sqrt_prop, design = mm)


# add cluster values to fitted models object
p <- matrix(as.numeric(prop) * 100, ncol = length(fs))
colnames(p) <- gsub("^PBMC8_30min_|\\.fcs$", "", basename(fs))

fit$genes <- p


# empirical Bayes tests for differentially abundant clusters
fit <- eBayes(fit)


# top differentially abundant clusters
topTable(fit, coef = 2, n = 10)




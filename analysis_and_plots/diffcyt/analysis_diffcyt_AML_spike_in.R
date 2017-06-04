##########################################################################################
# Analysis and plots: evaluate methods for data set 'AML-spike-in'
#
# Lukas Weber June 2017
##########################################################################################


library(diffcyt)
library(ggplot2)
library(reshape2)


# load saved results
load("../../../RData/outputs_diffcyt_AML_spike_in.RData")



########################
# Clustering performance
########################

# how accurately does the high-resolution clustering detect the spiked-in populations?

clustering_pr
clustering_re
clustering_F1


# thresholds
thresholds <- c(0.01, 0.001, 0.0001)
thresholds_nm <- c("1pc", "0.1pc", "0.01pc")

# conditions: CN, CBF
conditions <- c("CN", "CBF")


# plot clustering performance: precision, recall, F1 score

for (th in 1:length(thresholds)) {
  for (cnd in 1:length(conditions)) {
    
    d_plot <- data.frame(precision = clustering_pr[[th]][[cnd]][group_IDs == conditions[cnd]], 
                         recall = clustering_re[[th]][[cnd]][group_IDs == conditions[cnd]], 
                         F1_score = clustering_F1[[th]][[cnd]][group_IDs == conditions[cnd]])
    d_plot$patient <- gsub("^.*_", "", sample_IDs[group_IDs == conditions[cnd]])
    
    d_plot <- melt(d_plot, id.vars = "patient")
    
    ggplot(d_plot, aes(x = patient, y = value, group = variable, fill = variable)) + 
      geom_bar(stat = "identity", position = "dodge") + 
      scale_fill_manual(values = c("lightblue", "darkblue", "orange")) + 
      ylim(c(0, 1)) + 
      ggtitle(paste0("Clustering performance: AML-spike-in, ", conditions[cnd], ", threshold ", thresholds[th] * 100, "%")) + 
      theme_bw() + 
      theme(axis.title.y = element_blank(), 
            legend.title = element_blank())
    
    ggsave(paste0("../../../plots/diffcyt/AML_spike_in/", thresholds_nm[th], "/clustering_performance/", 
                  "clustering_perf_AML_spike_in_", thresholds_nm[th], "_", conditions[cnd], ".pdf"), 
           width = 6, height = 5)
  }
}




##########################################
# Differential abundance (DA) test results
##########################################

# plot DA test results

for (th in 1:length(thresholds)) {
  for (cnd in 1:length(conditions)) {
    path <- paste0("../../../plots/diffcyt/AML_spike_in/", thresholds_nm[th], "/DA/", cond_names[cnd])
    plotMST(d_se, d_counts, res_DA = out_DA[[th]][[cnd]], type = "DA", pvalue_type = "raw", path = path)
  }
}




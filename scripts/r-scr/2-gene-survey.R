#
# this part of the research explores 
# sequence heterogeneity for different HLA genes
#
# alleles which differ by first set of digits
# were chosen for the investigation
#

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("msa")
require(msa)

#
# reading FASTA file with protein sequences of HLA-A, HLA-B and HLA-C, 
# aligning them (separately) with ClustalOmega algorithm, then converting 
# each aligned file to specific format that can be modified to matrix
# 

HLA_A_aln <- readAAStringSet("./output/hla_a.fasta") %>% 
  msa(method = "ClustalOmega") %>% 
  msaConvert("bios2mds::align")
  
HLA_B_aln <- readAAStringSet("./output/hla_b.fasta") %>% 
  msa(method = "ClustalOmega") %>% 
  msaConvert("bios2mds::align") 

HLA_C_aln <- readAAStringSet("./output/hla_c.fasta") %>% 
  msa(method = "ClustalOmega") %>% 
  msaConvert("bios2mds::align")

#
# all rows of the matrix then are combined pairwise (for pairwise comparisons), 
# and these pairs are placed under each other in matrix 
# (brief explanation: rows 1 and 2 - one pair, rows 3 and 4 - another pair etc.)
#

HLA_A_combn <- matrix(unlist(HLA_A_aln), nrow = length(HLA_A_aln), byrow = T) %>% 
  data.frame() %>% 
  unname() %>% 
  apply(2, combn, m = 2)
  
HLA_B_combn <- matrix(unlist(HLA_B_aln), nrow = length(HLA_B_aln), byrow = T) %>% 
  data.frame() %>% 
  unname() %>% 
  apply(2, combn, m = 2)

HLA_C_combn <- matrix(unlist(HLA_C_aln), nrow = length(HLA_C_aln), byrow = T) %>% 
  data.frame() %>% 
  unname() %>% 
  apply(2, combn, m = 2)

#
# below one can see how function for counting number of mismatches  
# is applied to each of 3 matrices (output is a vector)
# (the function is made so that it does not take gaps for mismatches)
#

mism_A <- mapply(function(i,z) {
  sum(mapply(function(k,p) {
    ifelse((HLA_A_combn[i, ][p] != HLA_A_combn[z, ][k]) & (HLA_A_combn[i, ][p] != "-") & (HLA_A_combn[z, ][k] != "-"), T, F)
    }, k = seq(1, ncol(HLA_A_combn)), p = seq(1, ncol(HLA_A_combn))))
  }, i = seq(1, nrow(HLA_A_combn), 2), z = seq(2, nrow(HLA_A_combn), 2))

mism_B <- mapply(function(i,z) {
  sum(mapply(function(k,p) {
    ifelse((HLA_B_combn[i, ][p] != HLA_B_combn[z, ][k]) & (HLA_B_combn[i, ][p] != "-") & (HLA_B_combn[z, ][k] != "-"), T, F)
  }, k = seq(1, ncol(HLA_B_combn)), p = seq(1, ncol(HLA_B_combn))))
}, i = seq(1, nrow(HLA_B_combn), 2), z = seq(2, nrow(HLA_B_combn), 2))

mism_C <- mapply(function(i,z) {
  sum(mapply(function(k,p) {
    ifelse((HLA_C_combn[i, ][p] != HLA_C_combn[z, ][k]) & (HLA_C_combn[i, ][p] != "-") & (HLA_C_combn[z, ][k] != "-"), T, F)
  }, k = seq(1, ncol(HLA_C_combn)), p = seq(1, ncol(HLA_C_combn))))
}, i = seq(1, nrow(HLA_C_combn), 2), z = seq(2, nrow(HLA_C_combn), 2))

# deleting unnecessary dataframes as they take up to much space

rm(HLA_A_aln, HLA_B_aln, HLA_C_aln, 
   HLA_A_combn, HLA_B_combn, HLA_C_combn)

require(ggplot2)

ggplot(as.data.frame(mism_A), aes(x = mism_A))+
  geom_histogram(aes(y=..density..), fill = "white", color = "black", binwidth = 1)+
  labs(x = "Number of mismatches", y = "Portion of pairs")+
  ggtitle("HLA-A")+
  geom_density(alpha=.2, fill="#FF0000")+
  geom_vline(aes(xintercept = median(mism_A)), color="#FF0000", linetype="dashed", size=1) #график для HLA-A

# saving the plot
ggsave("HLA_A_mismatch.png", path = "./output/plots", width = 3, height = 2, units = "in")

ggplot(as.data.frame(mism_B), aes(x = mism_B))+
  geom_histogram(aes(y=..density..), fill = "white", color = "black", binwidth = 1)+
  labs(x = "Number of mismatches", y = "Portion of pairs")+
  ggtitle("HLA-B")+
  geom_density(alpha=.2, fill="#00FF00")+
  geom_vline(aes(xintercept = median(mism_B)), color="#00FF00", linetype="dashed", size=1) #график для HLA-B

# saving the plot
ggsave("HLA_B_mismatch.png", path = "./output/plots", width = 3, height = 2, units = "in")

# plot below depicts bimodal distribution for HLA-C alleles

ggplot(as.data.frame(mism_C), aes(x = mism_C))+
  geom_histogram(aes(y=..density..), fill = "white", color = "black", binwidth = 1)+
  labs(x = "Number of mismatches", y = "Portion of pairs")+
  ggtitle("HLA-C")+
  geom_density(alpha=.2, fill="#0000FF")+
  geom_vline(aes(xintercept = median(mism_C)), color="#0000FF", linetype="dashed", size=1) #график для HLA-C

ggsave("HLA_C_mismatch.png", path = "./output/plots", width = 3, height = 2, units = "in")

# dataframe made from 3 vectors: mism_A, mism_B, mism_C
max_len <- max(c(length(mism_A), length(mism_B), length(mism_C)))
mism_data <- data.frame(c(mism_A, rep(NA, max_len - length(mism_A))), 
                        c(mism_B, rep(NA, max_len - length(mism_B))), 
                        c(mism_C, rep(NA, max_len - length(mism_C))))

colnames(mism_data) <- c("HLA-A", "HLA-B", "HLA-C") # renaming columns

require(reshape2)
require(dplyr)

# reshaping data for plotting
mism_data <- mism_data %>% 
  melt() %>% 
  rename(Gene = variable, Value = value)

ggplot(mism_data, aes(x = Value, group = Gene, color=Gene))+
  geom_density(alpha=0, size = 1.7)+
  labs(x = "Number of mismatches", y = "Density")+
  geom_vline(aes(xintercept = median(mism_B)), color="#17b12b", linetype="dashed", size=1)+
  geom_vline(aes(xintercept = median(mism_A)), color="#f35e5a", linetype="dashed", size=1)+
  geom_vline(aes(xintercept = median(mism_C)), color="#5581fc", linetype="dashed", size=1)+
  scale_color_manual(values = c("#f35e5a","#17b12b", "#5581fc"))+
  theme_bw()+
  theme(legend.key.size = unit(1, 'cm'))+
  theme(legend.title = element_text(size = 15))+
  theme(legend.text = element_text(size = 10))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15))

# saving the plot
ggsave("mism_distr_dens.png", path = "./output/plots", width = 4.3, height = 3.2, scale = 1.8, units = "in")

ggplot(mism_data, aes(x = Gene, y = Value, group = Gene, color=Gene))+
  geom_violin(trim = T, adjust = 0.8)+
  theme(legend.title = element_blank())+
  geom_boxplot(width = .1, aes(fill = Gene))+
  stat_summary(fun = median, geom = "point", fill = "white", shape = 22,
               size = 2.5)+
  xlab("Gene") + ylab("Number of mismatches")+
  scale_fill_manual(values = alpha(rainbow(3), 0.7))+
  theme_bw()

# saving the plot
ggsave("mism_distr_viol.png", path = "./output/plots", width = 3, height = 2, scale = 1.8, units = "in")

require(forcats)
require(ggridges)

ggplot(mutate(mism_data, variable = fct_rev(as_factor(Gene))), aes(x = Value, y=Gene)) +
  geom_density_ridges(aes(fill = variable), show.legend = FALSE) + theme_ridges() +
  scale_fill_manual(values = c("#4682b4", "#6ed42a", "#f70a1e"))+
  scale_y_discrete(expand = c(0.01, 0))+
  scale_x_continuous(expand = c(0.01, 0)) +
  labs(x = "Number of mismatches", y = "Gene")+
  theme_minimal()+
  ggtitle("Distribution of mismatches")+
  theme(plot.title = element_text(hjust = 0.5), 
        plot.margin = unit(c(1, 1, 1, 1), "cm"), 
        axis.text.x = element_text(vjust = 3, margin=margin(t=10)),  
        axis.text.y = element_text(hjust = 3), 
        text = element_text(size = 15))

# saving the plot
ggsave("mism_distr_ridg.png", path = "./output/plots", width = 3, height = 2, bg="white", scale = 1.8, units = "in")

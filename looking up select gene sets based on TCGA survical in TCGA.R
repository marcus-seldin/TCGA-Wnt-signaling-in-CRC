setwd('')

library(ggplot2)
library(ggrepel)
library(survival)
library(survminer)
library(survMisc)
library(WGCNA)
library(reshape2)
library(factoextra)
library(dplyr)
library(MetBrewer)
library(pheatmap)
library(colormap)
library(qgraph)
library(rtracklayer)
library(bnlearn)
library(bnstruct)
library(enrichR)
library(FactoMineR)
allowWGCNAThreads()

load('panTCGAseq_somaticMutations_CNV.RData')

spec_gene_set = read.csv('./gene_lists/Rpos_xenoshort.csv')

############################################################################
gene_set = spec_gene_set$gene_symbol
id_set = 'Rpos_xenoshort'
cancer_type = 'colon adenocarcinoma'
#Lets look at mutations per cancer type, as the first step in TCGA data hub
unique(panTCGA_somaticMutations$effect)
table(panTCGA_somaticMutations$cancer_type)
LOF_variants = panTCGA_somaticMutations[panTCGA_somaticMutations$effect=='Missense_Mutation' | panTCGA_somaticMutations$effect=='Frame_Shift_Del' | panTCGA_somaticMutations$effect=='In_Frame_Del' | panTCGA_somaticMutations$effect=='large deletion',]


LOF_variants = LOF_variants[LOF_variants$gene %in% gene_set,]

LL1 = LOF_variants[LOF_variants$cancer_type==cancer_type,]
LL1 = na.omit(LL1)

#Inspect the number of mutations which occur as a result of cancer types
pdf(file = paste0(id_set, ' ', cancer_type, ' LOF mutations per gene.pdf'))
ggplot(LL1, aes(x=gene, fill=gene)) + geom_bar() +  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.position = "none") + ylab('LOF mutation count #') + xlab('') + ggtitle(paste0(id_set, ' ', cancer_type, ' LOF mutations per gene')) 
dev.off()

surv_data = read.delim('panTCGA_survival.txt')

new_surv = surv_data
new_surv$clock_LOF = ifelse(new_surv$sample %in% LOF_variants$sample, 'LOF_mutation', 'WT')
new_surv$cancer_type = panTCGA_somaticMutations$cancer_type[match(new_surv$sample, panTCGA_somaticMutations$sample)]
table(new_surv$cancer_type)
new_surv1 = new_surv[!is.na(new_surv$cancer_type),]

#List the cancer types of interest
cancer_subtypes = cancer_type

new_surv2 = new_surv1[new_surv1$cancer_type==cancer_subtypes,]
table(new_surv2$race)
new_surv2 = new_surv2[!is.na(new_surv2$cancer_type),]
pdf(file = paste0(id_set, ' cancer subtypes by race - ', cancer_subtypes, '.pdf'))
ggplot(new_surv2, aes(x=race, y=OS.time, fill=clock_LOF)) + geom_bar(stat="identity", position=position_dodge()) +  theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_fill_manual(values=c('darkorchid3', 'darkorange2')) + xlab('') + ylab("survival time (days)") + ggtitle(paste0('LOF mutations per race ', cancer_subtypes))
dev.off()


#Filter the sequencing data for only those genes and only cancer types
filtered_expression = panTCGA_seq[panTCGA_seq$gene_symbol %in% gene_set & panTCGA_seq$cancer_type %in% cancer_subtypes,]

pdf(file = paste0(id_set, ' gene expression in tumor vs surrounding tissue ', cancer_subtypes, '.pdf'))
ggplot(filtered_expression, aes(x=gene_symbol, y=expression_value, fill=tumor_class)) + 
  geom_boxplot() + scale_fill_manual(values=c('coral', 'dodgerblue2')) + theme_minimal() + xlab('') + ylab('normalized expression level') + ggtitle('Gene expression in tumor vs surrounding tissue') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

tt2 = as.data.frame(unique(filtered_expression$gene_symbol))
ttt3 = as.data.frame(tt2[1:(length(row.names(tt2))/3),])
colnames(ttt3) = 'gene1'
ff1 = filtered_expression[filtered_expression$gene_symbol %in% ttt3$gene1,]

pdf(file = paste0(id_set, ' gene expression in tumor vs surrounding tissue ', cancer_subtypes, ' - SET1.pdf'))
ggplot(ff1, aes(x=gene_symbol, y=expression_value, fill=tumor_class)) + 
  geom_boxplot() + scale_fill_manual(values=c('coral', 'dodgerblue2')) + theme_minimal() + xlab('') + ylab('normalized expression level') + ggtitle('Gene expression in tumor vs surrounding tissue') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

ttt3 = as.data.frame(tt2[(length(row.names(tt2))/3):(2*(length(row.names(tt2))/3)),])
colnames(ttt3) = 'gene1'
ff1 = filtered_expression[filtered_expression$gene_symbol %in% ttt3$gene1,]

pdf(file = paste0(id_set, ' gene expression in tumor vs surrounding tissue ', cancer_subtypes, ' - SET2.pdf'))
ggplot(ff1, aes(x=gene_symbol, y=expression_value, fill=tumor_class)) + 
  geom_boxplot() + scale_fill_manual(values=c('coral', 'dodgerblue2')) + theme_minimal() + xlab('') + ylab('normalized expression level') + ggtitle('Gene expression in tumor vs surrounding tissue') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


ttt3 = as.data.frame(tt2[(2*(length(row.names(tt2))/3)): length(row.names(tt2)),])
colnames(ttt3) = 'gene1'
ff1 = filtered_expression[filtered_expression$gene_symbol %in% ttt3$gene1,]

pdf(file = paste0(id_set, ' gene expression in tumor vs surrounding tissue ', cancer_subtypes, ' - SET3.pdf'))
ggplot(ff1, aes(x=gene_symbol, y=expression_value, fill=tumor_class)) + 
  geom_boxplot() + scale_fill_manual(values=c('coral', 'dodgerblue2')) + theme_minimal() + xlab('') + ylab('normalized expression level') + ggtitle('Gene expression in tumor vs surrounding tissue') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

#Finally, we will investigate survival curves.  First we will look at loss of function mutations in select genes
#The somatic mutation data is already loaded as panTCGA_somaticMutations.  Firts filter for our genes and cancer
select_mut_data = panTCGA_somaticMutations[panTCGA_somaticMutations$gene %in% gene_set & panTCGA_somaticMutations$cancer_type %in% cancer_subtypes,]
#See the mutation types
table(select_mut_data$effect)
#Filter for potential protein altering variations.  We pick from the command above
LOF_variants = select_mut_data[select_mut_data$effect=='Missense_Mutation' | select_mut_data$effect=='Frame_Shift_Del' | select_mut_data$effect=='In_Frame_Del',]

#This leaves us with 25 people total
nrow(LOF_variants)

#Add the IDs for LOF individuals
plotting_traits = surv_data[surv_data$sample %in% filtered_expression$TCGA_ID,]
plotting_traits$LOF_category = ifelse(plotting_traits$sample %in% LOF_variants$sample, 'LOF', 'WT_allele')
table(plotting_traits$LOF_category)

#Now calculate stats and generate survival curve
kmfit <- survfit(Surv(plotting_traits$OS.time, plotting_traits$OS) ~ plotting_traits$LOF_category)
pdf(file = paste0(id_set, ' LOF vs WT mutations survival curve ', cancer_subtypes, '.pdf'))
ggsurvplot(
  kmfit,
  data = plotting_traits,
  size = 1,                 # change line size
  palette =
    c("#E7B800", "#2E9FDF"),# custom color palettes
  conf.int = T,          # Add confidence interval
  pval = T,              # Add p-value
  risk.table = F,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  #legend.labs =
  # c("concordant", "discordant"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)
dev.off()

#Lets bin individuals into expression categories based on genes of interest.  To do this, we will first compute the average expression of each gene
tumor_samples = filtered_expression[filtered_expression$tumor_class=='tumor',]
avg_expresision = tumor_samples %>% dplyr::select(gene_symbol, expression_value) %>% group_by(gene_symbol) %>% summarise(mean=mean(expression_value))

#Create a 'category' for each individual as to their mean expression
tumor_samples$avg_gene_expr = avg_expresision$mean[match(tumor_samples$gene_symbol, avg_expresision$gene_symbol)]
tumor_samples$expr_category = ifelse(tumor_samples$expression_value > tumor_samples$avg_gene_expr, 'High_expr', 'Low_expr')

#Summarize number of genes in each expre category
ind_list = tumor_samples %>% select(TCGA_ID, expr_category) %>% group_by(TCGA_ID) %>% count(expr_category)
ind_list = dcast(ind_list, TCGA_ID ~expr_category, fun.aggregate = sum, value.var = 'n')

#make a final column.  If the counts of high outweighs the low expr we consider High
ind_list$final_expr_cat = ifelse(ind_list$High_expr > ind_list$Low_expr, 'High_expr', 'Low_expr')
table(ind_list$final_expr_cat)

#Filter only for individuals who 
plotting_traits = surv_data[surv_data$sample %in% ind_list$TCGA_ID,]

#Add the expr category back for each individual onto the traits data
plotting_traits$expr_category = ind_list$final_expr_cat[match(plotting_traits$sample, ind_list$TCGA_ID)]
table(plotting_traits$expr_category)

pdf(file = paste0(id_set, ' High vs Low gene expression survival curve.pdf'))
kmfit <- survfit(Surv(plotting_traits$OS.time, plotting_traits$OS) ~ plotting_traits$expr_category)
ggsurvplot(
  kmfit,
  data = plotting_traits,
  size = 1,                 # change line size
  palette =
    c("#E7B800", "#2E9FDF"),# custom color palettes
  conf.int = T,          # Add confidence interval
  pval = T,              # Add p-value
  risk.table = F,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  #legend.labs =
  # c("concordant", "discordant"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)
dev.off()


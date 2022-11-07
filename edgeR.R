library(edgeR)
library(ggplot2)
library(tidyverse)
library(ggsignif)
library(ggpubr)

#define all the input variables, including where the files are and the sample/group names  
path <- 'C://Users//chand//Box Sync//Krasileva_Lab//Research//chandler//Krasileva Lab//E14//counts/'
files <- c('sample_01_Aligned.out.tsv', 'sample_02_Aligned.out.tsv', 'sample_03_Aligned.out.tsv', 'sample_04_Aligned.out.tsv', 
           'SRR17281233.tsv', 'SRR17281234.tsv', 'SRR17281235.tsv', 'SRR17281236.tsv')
samples <- c('sample_01', 'sample_02', 'sample_03', 'sample_04', 'SRR17281233', 'SRR17281234', 'SRR17281235', 'SRR17281236')
groups <- c('A', 'A', 'A', 'A', 'B', 'B', 'B', 'B')

#define counts, the DGE object 
counts <- readDGE(files, labels=samples, group=groups, path=path)

#filter out genes that aren't expressed, and add a normalization factor for library size 
keep <- filterByExpr(counts)
counts <- counts[keep, , keep.lib.sizes=FALSE]
counts <- calcNormFactors(counts)

#check the distance between samples
plotMDS(counts, labels=samples, col=rep(1:2, each=4))

#define the design matrix, which turns the experimental design into a binary 
design <- model.matrix(~groups)
rownames(design) <- colnames(counts)
design

#estimate the dispersion 
counts <- estimateDisp(counts, design, robust=TRUE)
counts$common.dispersion
plotBCV(counts)

#fit gene wise glms
fit <- glmQLFit(counts, design, robust=TRUE)

#conduct likelihood ratio test for A vs B and show top genes. positive logFC would be higher expression in real samples 
qlf <- glmQLFTest(fit)
topTags(qlf)

#define the indexes of nonhv and hv to perform gene set tests
nonhv <- working_table %>% filter(HV == 0) %>% pull(Gene)
hv <- working_table %>% filter(HV==1) %>% pull(Gene)
nonhv_index <- rownames(counts) %in% nonhv
hv_index <- rownames(counts) %in% hv
index <- list(nonhv=nonhv_index, hv=hv_index)

#compare hv and nonhv expression 
fry(counts, index = index, design=design)
camera(counts, index=index, design=design)

barcodeplot(qlf$table$logFC, index=index[[2]], index2=index[[1]], xlab='logFC', col.bars=c('#F8766D', '#00BFC4'))

with(qlf$table, plot(logCPM,logFC, pch=16, cex=0.2))
with(qlf$table, points(logCPM[nonhv_index], logFC[nonhv_index], pch=16, col='#00BFC4'))
with(qlf$table, points(logCPM[hv_index], logFC[hv_index], pch=16, col='#F8766D'))
legend("bottomright", legend=c("nonHV NLRs", "HV NLRs"), pch=16, col=c('#F8766D', '#00BFC4'))


working_table <- readxl::read_xlsx(path="C://Users//chand//Box Sync//Krasileva_Lab//Research//chandler//Krasileva Lab//E14//e14_R//working_table.xlsx")
working_table

qlf$table$names <- rownames(qlf)

expression <- merge(qlf$table, working_table, by.x='names', by.y='Gene')
just_rpp4 <- expression %>% filter(cluster == 'DM8_RPP4-5')
ggplot((expression %>% filter(Nterm != 'X') %>% filter(Nterm != 'NA')),
       aes(x=Nterm, y=logFC, fill=factor(HV)))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge())+
  geom_point(data=just_rpp4, 
             aes(x=Nterm, y=logFC, color=cluster),
             position=position_dodge(width=0.75)) #+
  #stat_compare_means(method='t.test', 
  #                   label = 'p.signif')

exp <- expression %>% mutate(HV=recode(HV, `0` = "Non-highly variable", 
                               `1`="Highly variable"))

exp$HV <- factor(exp$HV, c("Non-highly variable", "Highly variable"))

ggplot((exp %>% filter(Nterm != 'X') %>% filter(Nterm != 'NA')),
       aes(x=factor(HV), y=logFC, fill=factor(HV)))+
  geom_boxplot(show.legend = FALSE)+
  geom_point(position=position_jitterdodge(), show.legend=FALSE)+
  labs(x='Gene Type', y='Expression')

all_expression <- merge(qlf$table, working_table, by.x='names', by.y='Gene', all=TRUE) %>%
  filter(!row_number() %in% c(1, 2, 3))

all_expression$Nterm[is.na(all_expression$Nterm)]<-'all_genes'

all_expression <- all_expression %>% filter(Nterm != 'X') %>% filter(names != 'AT4G19050')
ggplot((all_expression %>% filter(Nterm != 'X') %>% filter(Nterm != 'NA')),
       aes(x=Nterm, y=logFC, fill=factor(HV)))+
  geom_boxplot()+
  scale_fill_manual(values=c('grey', '#F8766D', '#00BFC4')) +
  geom_point(data=expression %>% filter(Nterm != 'X'), aes(x=Nterm, y=logFC), position=position_jitterdodge())+
  stat_compare_means(method='t.test', 
                     label = 'p.signif', 
                     bracket.size=0.3)


all_expression <- all_expression %>% mutate(HV=recode(HV, `0` = "nonHV", 
                                       `1`="HV"))
all_expression$HV[is.na(all_expression$HV)]<-'all_genes'
all_expression$HV <- factor(all_expression$HV , levels=c("all_genes", "nonHV", "HV"))

ggplot((all_expression %>% filter(Nterm != 'X') %>% filter(Nterm != 'NA')),
       aes(x=HV, y=logFC, fill=factor(HV)))+
  geom_boxplot()+
  scale_fill_manual(values=c('grey','#00BFC4',  '#F8766D')) +
  geom_point(data=all_expression %>% filter(HV == 'nonHV' | HV == 'HV'), 
             aes(x=HV, y=logFC), position=position_jitterdodge(), alpha=0.7)+
  geom_signif(comparisons=list(c('all_genes', 'nonHV'), c('all_genes', 'HV')), map_signif_level = TRUE, y_position = c(12, 8))+
  labs(fill = 'HV status',
       title = 'log Fold Change Relative to Simulated Expression Data')

ggplot((expression %>% filter(Nterm != 'X') %>% filter(Nterm != 'NA')),
       aes(x=cluster_type, y=logFC, fill=factor(HV)))+
  geom_boxplot()+
  scale_fill_manual(values=c('#F8766D', '#00BFC4')) +
  geom_point(data=expression %>% filter(Nterm != 'X'), aes(x=cluster_type, y=logFC), position=position_jitterdodge())+
  stat_compare_means(method='t.test', 
                     label = 'p.signif', 
                     bracket.size=0.3)

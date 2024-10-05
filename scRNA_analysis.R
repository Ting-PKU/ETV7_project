library(Seurat)
library(ggplot2)
library(RColorBrewer)
#######
#load data
setwd("~/projects/singlecell/Tcell/proj2/")
cd8 = readRDS('~/projects/singlecell/Tcell/data/cd8.rds')
Idents(cd8) = 'meta.cluster'
colSet.CD8 <- readRDS("~/projects/singlecell/Tcell/data/metaInfo/panC.colSet.list.rds")
colSet.CD8 <- colSet.CD8$meta.cluster[grep('CD8', names(colSet.CD8$meta.cluster))]
DimPlot(cd8,reduction = 'harmony.umap',cols = colSet.CD8)

#Memory and exhaustion score
cd8 = AddModuleScore(cd8, features=list(c('EOMES','CCR4','CCR7',
                                          'CXCR3','CXCR4','CXCR5','TCF7')), name="Memory.Score")
cd8 = AddModuleScore(cd8, features=list(c('PDCD1','HAVCR2','LAG3','TOX',
                                          'CXCL13','TIGIT','CTLA4','TNFRSF9')), name="Exhaustion.Score")
FeaturePlot(cd8, 'Memory.Score1',reduction = 'harmony.umap',
            max.cutoff = 'q90',min.cutoff = 'q10',order = T)&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 5, name = "RdBu")))&
  xlab('Umap1') &ylab('Umap2')&ggtitle('Memory score')

FeaturePlot(cd8, 'Exhaustion.Score1',reduction = 'harmony.umap',
            max.cutoff = 'q90',min.cutoff = 'q10',order = T)&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 5, name = "RdBu")))

FeaturePlot(cd8, 'ETV7',reduction = 'harmony.umap',
            max.cutoff = 'q90',min.cutoff = 'q10',order = T)&
  scale_colour_gradientn(colours = rev(brewer.pal(n = 5, name = "RdBu")))

############## 
#correlation ETV7 and memory/exhaustion markers
# add group
clu.p1 = c("CD8.c02.Tm.IL7R",
           "CD8.c05.Tem.CXCR5","CD8.c06.Tem.GZMK","CD8.c11.Tex.PDCD1","CD8.c12.Tex.CXCL13")
clu.p2 = c("CD8.c05.Tem.CXCR5","CD8.c12.Tex.CXCL13")
meta = cd8@meta.data
meta$ETV7 = cd8@assays$RNA@data['ETV7',]
meta$group = "Others"
meta[meta$meta.cluster %in% clu.p1 , "group"] = "P1"
meta$group2 = "Others"
meta[meta$meta.cluster %in% clu.p2 , "group2"] = "P2"
p_list <- list()
for(gene in c('TCF7', 'CXCR4','CXCR5', 'CCR7', 'EOMES')){
  meta$TOX = cd8@assays$RNA@data[gene,]
  sig_cor = data.frame(cancer = unique(meta$cancerType),estimate = 0,
                       pvalue = 0,estimate_1 = 0,
                       pvalue_1 = 0,estimate_2 = 0,
                       pvalue_2 = 0)
  row.names(sig_cor) = sig_cor$cancer
  for(i in sig_cor$cancer){
    print(i)
    tmp = meta[meta$cancerType==i,]
    sig_cor[i,2] = cor.test(tmp$ETV7,tmp$TOX)$estimate
    sig_cor[i,3] = cor.test(tmp$ETV7,tmp$TOX)$p.value
    tmp = meta[meta$group=="P1"&meta$cancerType==i,]
    sig_cor[i,4] = cor.test(tmp$ETV7,tmp$TOX)$estimate
    sig_cor[i,5] = cor.test(tmp$ETV7,tmp$TOX)$p.value
    tmp = meta[meta$group2=="P2"&meta$cancerType==i,]
    sig_cor[i,6] = cor.test(tmp$ETV7,tmp$TOX)$estimate
    sig_cor[i,7] = cor.test(tmp$ETV7,tmp$TOX)$p.value
    
  }
  p_list[[gene]] <- ggplot(na.omit(sig_cor),
                           aes(x=reorder(cancer,-estimate_1,mean),estimate_1))+
    geom_bar(stat = 'identity',aes(fill = cancer),
             alpha=0.8) +
    scale_fill_manual(values=g.colSet$cancerType) +
    coord_cartesian(clip="off",ylim=c(-0.5,0.2)) +theme_pubr()+
    theme(plot.title = element_text(hjust = 0.5,size=14),
          axis.text.x = element_text(angle = 60, hjust = 1,vjust=1),
          legend.position = 'none')+ggtitle(paste0('ETV7 ~ ',gene))+
    xlab("")+scale_y_continuous(expand = c(0,0))+
    ylab(paste0("Pearson correlation coefficient of\nETV7 and ",gene))
}
p1_list <- p_list
p_list <- list()
for(gene in c('PDCD1','HAVCR2','LAG3','TOX',
              'CXCL13','CTLA4','TNFRSF9')){
  meta$TOX = cd8@assays$RNA@data[gene,]
  sig_cor = data.frame(cancer = unique(meta$cancerType),estimate = 0,
                       pvalue = 0,estimate_1 = 0,
                       pvalue_1 = 0,estimate_2 = 0,
                       pvalue_2 = 0)
  row.names(sig_cor) = sig_cor$cancer
  for(i in sig_cor$cancer){
    print(i)
    tmp = meta[meta$cancerType==i,]
    sig_cor[i,2] = cor.test(tmp$ETV7,tmp$TOX)$estimate
    sig_cor[i,3] = cor.test(tmp$ETV7,tmp$TOX)$p.value
    tmp = meta[meta$group=="P1"&meta$cancerType==i,]
    sig_cor[i,4] = cor.test(tmp$ETV7,tmp$TOX)$estimate
    sig_cor[i,5] = cor.test(tmp$ETV7,tmp$TOX)$p.value
    tmp = meta[meta$group2=="P2"&meta$cancerType==i,]
    sig_cor[i,6] = cor.test(tmp$ETV7,tmp$TOX)$estimate
    sig_cor[i,7] = cor.test(tmp$ETV7,tmp$TOX)$p.value
    
  }
  p_list[[gene]] <- ggplot(na.omit(sig_cor),
                           aes(x=reorder(cancer,estimate_1,mean),estimate_1))+
    geom_bar(stat = 'identity',aes(fill = cancer),
             alpha=0.8) +
    scale_fill_manual(values=g.colSet$cancerType) +
    coord_cartesian(clip="off",ylim=c(-0.2,0.8)) +theme_pubr()+
    theme(plot.title = element_text(hjust = 0.5,size=14),
          axis.text.x = element_text(angle = 60, hjust = 1,vjust=1),
          legend.position = 'none')+ggtitle(paste0('ETV7 ~ ',gene))+
    xlab("")+scale_y_continuous(expand = c(0,0))+
    ylab(paste0("Pearson correlation coefficient of\nETV7 and ",gene))
}
library(patchwork)
p<-wrap_plots(p_list, ncol = 4);p
p<-wrap_plots(p1_list, ncol = 4);p

#correlation ETV7 and memory/exhaustion score
sig_cor = data.frame(cancer = unique(meta$cancerType),estimate = 0,
                     pvalue = 0,estimate_1 = 0,
                     pvalue_1 = 0,estimate_2 = 0,
                     pvalue_2 = 0)
row.names(sig_cor) = sig_cor$cancer
for(i in sig_cor$cancer){
  print(i)
  tmp = meta[meta$cancerType==i,]
  sig_cor[i,2] = cor.test(tmp$ETV7,tmp$Exhaustion.Score1)$estimate
  sig_cor[i,3] = cor.test(tmp$ETV7,tmp$Exhaustion.Score1)$p.value
  tmp = meta[meta$group=="P1"&meta$cancerType==i,]
  sig_cor[i,4] = cor.test(tmp$ETV7,tmp$Exhaustion.Score1)$estimate
  sig_cor[i,5] = cor.test(tmp$ETV7,tmp$Exhaustion.Score1)$p.value
  tmp = meta[meta$group2=="P2"&meta$cancerType==i,]
  sig_cor[i,6] = cor.test(tmp$ETV7,tmp$Exhaustion.Score1)$estimate
  sig_cor[i,7] = cor.test(tmp$ETV7,tmp$Exhaustion.Score1)$p.value
  
}
ggplot(na.omit(sig_cor),
       aes(x=reorder(cancer,estimate_1,mean),estimate_1))+
  geom_bar(stat = 'identity',aes(fill = cancer),
           alpha=0.8) +
  scale_fill_manual(values=g.colSet$cancerType) +
  coord_cartesian(clip="off",ylim=c(-0.1,0.8)) +theme_pubr()+
  theme(plot.title = element_text(hjust = 0.5,size=14),
        axis.text.x = element_text(angle = 60, hjust = 1,vjust=1),
        legend.position = 'none')+ggtitle('')+
  xlab("")+scale_y_continuous(expand = c(0,0))+
  ylab("Pearson correlation coefficient of\nETV7 and Exhaustion score")

for(i in sig_cor$cancer){
  print(i)
  tmp = meta[meta$cancerType==i,]
  sig_cor[i,2] = cor.test(tmp$ETV7,tmp$Memory.Score1)$estimate
  sig_cor[i,3] = cor.test(tmp$ETV7,tmp$Memory.Score1)$p.value
  tmp = meta[meta$group=="P1"&meta$cancerType==i,]
  sig_cor[i,4] = cor.test(tmp$ETV7,tmp$Memory.Score1)$estimate
  sig_cor[i,5] = cor.test(tmp$ETV7,tmp$Memory.Score1)$p.value
  tmp = meta[meta$group2=="P2"&meta$cancerType==i,]
  sig_cor[i,6] = cor.test(tmp$ETV7,tmp$Memory.Score1)$estimate
  sig_cor[i,7] = cor.test(tmp$ETV7,tmp$Memory.Score1)$p.value
}

ggplot(na.omit(sig_cor),
       aes(x=reorder(cancer,-estimate_1,mean),estimate_1))+
  geom_bar(stat = 'identity',aes(fill = cancer),
           alpha=0.8) +
  scale_fill_manual(values=g.colSet$cancerType) +
  coord_cartesian(clip="off",ylim=c(-0.6,0)) +theme_pubr()+
  theme(plot.title = element_text(hjust = 0.5,size=14),
        axis.text.x = element_text(angle = 60, hjust = 1,vjust=1),
        legend.position = 'none')+ggtitle('')+
  xlab("")+scale_y_continuous(expand = c(0,0))+
  ylab("Pearson correlation coefficient of\nETV7 and Memory score")
# plot correlation result by cancer types
# UCEC
dat = meta[meta$group=='P1'&meta$cancerType=='UCEC',]
df = data.frame(ETV7 = dat$ETV7,Exhaust.Score = dat$Exhaustion.Score1)
cor.test(df$ETV7,df$Exhaust.Score)$estimate
cor.test(df$ETV7,df$Exhaust.Score)$p.value
plot(df$ETV7,df$Exhaust.Score ,las = 1,xlab = 'Expression levle of ETV7',
     ylab = 'Exhaustion score',#ylim = c(-0.6,2.2),xlim = c(-0,0.5),
     cex=1,col = "#9E6398", main = 'UCEC',
     pch = 8)
abline(lm(df$Exhaust.Score~df$ETV7), col = "#C49FC4", lwd = 2)
text(0.2,1.4,'R = 0.49, p-value = 3.76e-20')

dat = meta[meta$group=='P1'&meta$cancerType=='UCEC',]
df = data.frame(ETV7 = dat$ETV7,Memory.Score = dat$Memory.Score1)
cor.test(df$ETV7,df$Memory.Score)$estimate
cor.test(df$ETV7,df$Memory.Score)$p.value
df1 = df[df$ETV7>0,]
plot(df$ETV7,df$Memory.Score ,las = 1,xlab = 'Expression levle of ETV7',
     ylab = 'Memory score',#ylim = c(-0.6,2.2),xlim = c(-0,0.5),
     cex=1,col = "#9E6398", main = 'UCEC',
     pch = 8)
abline(lm(df$Memory.Score~df$ETV7), col = "#C49FC4", lwd = 2)
text(0.3,0.7,'R = -0.43, p-value = 6.59e-15')

# ESCA
dat = meta[meta$group=='P1'&meta$cancerType=='ESCA',]
df = data.frame(ETV7 = dat$ETV7,Exhaust.Score = dat$Exhaustion.Score1)
cor.test(df$ETV7,df$Exhaust.Score)$estimate
cor.test(df$ETV7,df$Exhaust.Score)$p.value
plot(df$ETV7,df$Exhaust.Score ,las = 1,xlab = 'Expression levle of ETV7',
     ylab = 'Exhaustion score',#ylim = c(-0.6,2.2),xlim = c(-0,0.5),
     cex=1,col = "#9E6398", main = 'ESCA',
     pch = 8)

abline(lm(df$Exhaust.Score~df$ETV7), col = "#C49FC4", lwd = 2)
text(0.2,1.2,'R = 0.71, p-value = 2.4e-55')

dat = meta[meta$group=='P1'&meta$cancerType=='ESCA',]
df = data.frame(ETV7 = dat$ETV7,Memory.Score = dat$Memory.Score1)
cor.test(df$ETV7,df$Memory.Score)$estimate
cor.test(df$ETV7,df$Memory.Score)$p.value
plot(df$ETV7,df$Memory.Score ,las = 1,xlab = 'Expression levle of ETV7',
     ylab = 'Memory score',#ylim = c(-0.6,2.2),xlim = c(-0,0.5),
     cex=1,col = "#9E6398", main = 'ESCA',
     pch = 8)
abline(lm(df$Memory.Score~df$ETV7), col = "#C49FC4", lwd = 2)
text(0.3,1.2,'R = -0.48, p-value = 2.73e-16')
###### plot correlation results of all_cancertype
dat = readRDS('datascore_meta.rds')
df = data.frame(ETV7 = dat$ETV7,Exhaust.Score = dat$Exhaust.Score)
cor.test(df$ETV7,df$Exhaust.Score)$p.value
df = df[df$ETV7>0,]
plot(df$ETV7,df$Exhaust.Score,las = 1,xlab = 'Expression levle of ETV7',
     ylab = 'Exhaustion score',ylim = c(-0.6,1.5),xlim = c(-0,0.5),
     cex=1,col = "#0571B0", 
     pch = 8)
abline(lm(df$Exhaust.Score~df$ETV7), col = "#3B82F4", lwd = 2)
text(0.25,1.5,'R = 0.51, p-value = 1.7e-188')

df = data.frame(ETV7 = dat$ETV7,Exhaust.Score = dat$Memory.Score)
cor.test(df$ETV7,df$Exhaust.Score)$p.value
df = df[df$ETV7>0,]
plot(df$ETV7,df$Exhaust.Score,las = 1,xlab = 'Expression levle of ETV7',
     ylab = 'Memory Score',ylim = c(-0.4,1.2),xlim = c(-0,0.5),
     cex=1,col = "#0571B0", 
     pch = 8)
abline(lm(df$Exhaust.Score~df$ETV7), col = "#3B82F4", lwd = 2)
text(0.25,0.8,'R = -0.36, p-value = 3.8e-59')






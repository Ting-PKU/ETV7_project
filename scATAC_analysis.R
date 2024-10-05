library(ArchR)
setwd("/lustre/user/liclab/pengt/projects/singlecell/Tcell/proj2")
proj1 <- loadArchRProject('../public/scatac/BCC/fragments/project/proj1/')
plotEmbedding(ArchRProj = proj1, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
proj1 <- addTrajectory(
  ArchRProj = proj1, 
  name = "Exhaustion", 
  groupBy = "Clusters",
  trajectory = c("C8", "C3", "C4", "C16"), 
  embedding = "UMAPHarmony", 
  force = TRUE
)
genes <- c('PDCD1','HAVCR2','TOX',
           'TIGIT','CTLA4',
           'EOMES','CCR7',
           'CXCR4','CXCR5','TCF7')
p_list <- list()
for(gene in genes){
  p1 <- plotTrajectory(proj1, trajectory = "Exhaustion", colorBy = "GeneScoreMatrix", 
                       name = gene, continuousSet = "blueYellow",embedding = "UMAPHarmony")
  p_list[[gene]] <- p1[[1]]+ ggtitle(gene)+theme_pubr()+NoAxes()+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = 'right')
}
library(patchwork)
wrap_plots(p_list, ncol = 4)
proj1 <- addModuleScore(
  ArchRProj = proj1,
  useMatrix = 'GeneScoreMatrix',
  name = "Exhaustion_score",
  features = list(Escore = c('PDCD1','HAVCR2','LAG3','TOX',
                             'CXCL13','TIGIT','CTLA4'),
                  Mscore = c('EOMES','CCR4','CCR7',
                             'CXCR4','TCF7', 'CXCR5')),
  nBin = 25,
  nBgd = 100,
  seed = 2,
  threads = 12,
  logFile = createLogFile("addModuleScore")
)

p1 <- plotTrajectory(proj1, trajectory = "Exhaustion", colorBy = "cellColData", 
                     name = 'Exhaustion_score.Escore', continuousSet = "blueYellow",
                     embedding = "UMAPHarmony",quantCut = c(0.01, 0.96))
p1 <- p1[[1]]+ ggtitle('Exhaustion Score')+theme_pubr()+NoAxes()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'right')

p2 <- plotTrajectory(proj1, trajectory = "Exhaustion", colorBy = "cellColData", 
                     name = 'Exhaustion_score.Mscore', continuousSet = "blueYellow",
                     embedding = "UMAPHarmony",quantCut = c(0.01, 0.96))
p2<- p2[[1]]+ ggtitle('Memory Score')+theme_pubr()+NoAxes()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'right')
wrap_plots(c(list(p1, p2),p_list), ncol = 3)



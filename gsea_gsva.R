library(GSVA)
library(GSEABase)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(enrichplot)
library(limma)

######data#####
deg <- read.csv('DrugMvsC.txt',sep = '\t')
rownames(deg) <- deg$X
deg <- deg[,-1]
######msigdb#####
C2_df = msigdbr(species = "Mus musculus",collection = "C2") %>% 
  dplyr::select(gs_name,gene_symbol)

C5_df = msigdbr(species ="Mus musculus",collection = "C5") %>% 
  dplyr::select(gene_symbol,gs_name,gs_subcollection)
C5_df = C5_df[C5_df$gs_subcollection!="HPO",]
C5_df = C5_df[,c(2,1)]

#####GSEA#####
ge = deg$logFC
names(ge) = rownames(deg)
ge = sort(ge,decreasing = T)
head(ge)

#c2
C2em <- GSEA(ge, TERM2GENE = C2_df)
C2em_results <- as.data.frame(C2em@result)
library(openxlsx)
write.xlsx(C2em_results,'C2em_results.xlsx')

#c5
C5em <- GSEA(ge, TERM2GENE = C5_df)
# 提取完整结果表格（转换为数据框便于操作）
C5em_results <- as.data.frame(C5em@result)
write.xlsx(C5em_results,'C5em_results.xlsx')

Drug_c2_up <- C2em_results[C2em_results$NES>0,]
Drug_c2_down <- C2em_results[C2em_results$NES<0,]
Drug_c5_up <- C5em_results[C5em_results$NES>0,]
Drug_c5_down <- C5em_results[C5em_results$NES<0,]

######data#####
deg <- read.csv('../nature_prepare/NatureMvsC.txt',sep = '\t')
rownames(deg) <- deg$X
deg <- deg[,-1]

#####GSEA#####
ge = deg$logFC
names(ge) = rownames(deg)
ge = sort(ge,decreasing = T)
head(ge)

#c2
C2em <- GSEA(ge, TERM2GENE = C2_df)
C2em_results <- as.data.frame(C2em@result)
library(openxlsx)
write.xlsx(C2em_results,'NatureC2em_results.xlsx')

#c5
C5em <- GSEA(ge, TERM2GENE = C5_df)
# 提取完整结果表格（转换为数据框便于操作）
C5em_results <- as.data.frame(C5em@result)
write.xlsx(C5em_results,'NatureC5em_results.xlsx')

Nature_c2_up <- C2em_results[C2em_results$NES>0,]
Nature_c2_down <- C2em_results[C2em_results$NES<0,]
Nature_c5_up <- C5em_results[C5em_results$NES>0,]
Nature_c5_down <- C5em_results[C5em_results$NES<0,]

#####inner#####
c2c5up_venn_list <- list(Drugc2c5up = c(Drug_c2_up$ID,Drug_c5_up$ID), Naturec2c5up = c(Nature_c2_up$ID,Nature_c5_up$ID))
c2c5upinner <- get.venn.partitions(c2c5up_venn_list)
write.xlsx(c2c5upinner,'c2c5upinner.xlsx')
p_vennc2c5up<- venn.diagram(c2c5up_venn_list, imagetype = 'pdf', 
                           fill = c("#BC3C29FF","#0072B5FF"), alpha = 0.50, cat.col = rep('black', 2), 
                           col = 'black', cex = 1.5, fontfamily = 'serif', 
                           cat.cex = 1.5, cat.fontfamily = 'serif',scaled = FALSE,filename = NULL)

pdf('venn_c2c5up.pdf', width = 3, height = 3)
grid.draw(p_vennc2c5up)
dev.off()

c2c5upinner_pathways <- data.frame(c2c5upinner[[1,'..values..']])
colnames(c2c5upinner_pathways) <- c('c2c5upinner_pathways')
write.xlsx(c2c5upinner_pathways,'c2c5upinner_pathways.xlsx')

c2c5down_venn_list <- list(Drugc2c5down = c(Drug_c2_down$ID,Drug_c5_down$ID), Naturec2c5down = c(Nature_c2_down$ID,Nature_c5_down$ID))
c2c5downinner <- get.venn.partitions(c2c5down_venn_list)
write.xlsx(c2c5downinner,'c2c5downinner.xlsx')
p_vennc2c5down<- venn.diagram(c2c5down_venn_list, imagetype = 'pdf', 
                            fill = c("#BC3C29FF","#0072B5FF"), alpha = 0.50, cat.col = rep('black', 2), 
                            col = 'black', cex = 1.5, fontfamily = 'serif', 
                            cat.cex = 1.5, cat.fontfamily = 'serif',scaled = FALSE,filename = NULL)

pdf('venn_c2c5down.pdf', width = 3, height = 3)
grid.draw(p_vennc2c5down)
dev.off()
c2c5downinner_pathways <- data.frame(c2c5downinner[[1,'..values..']])
colnames(c2c5downinner_pathways) <- c('c2c5downinner_pathways')
write.xlsx(c2c5downinner_pathways,'c2c5downinner_pathways.xlsx')

# 'LONGECITY','CELLULAR','AGING'
# 'FATTY_ACID','PPARG','BILE_ACID'
# 'T_CELL','TH1','MTOR','T_HELPER','IL12'
# 'NFkB','inflama','TNFa'

# # Aging的通路
# Aging_pathways <-c('KEGG_LONGECITY_REGULATION','REACTOME_CELLULAR_SENESCENCE','BIOCARTA_AGING_PATHWAY')
# # Aging_pathways <-c('REACTOME_CELLULAR_SENESCENCE') 
# 
# #Fatty Acid Metabolism
# FA_pathways <- c('KEGG_FATTY_ACID_METABOLISM','KEGG_PPARG_SIGNALING_PATHWAY','REACTOME_FATTY_ACID_METABOLISM','PID_BILE_ACID_AND_FXR_PATHWAY')
# 
# #T Cell/Immune Cell Differentiation
# TI_pathways <- c('KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY','KEGG_TH1_AND_TH2_CELL_DIFFERENTIATION','KEGG_TH17_CELL_DIFFERENTIATION','REACTOME_DIFFERENTIATION_OF_T_HELPER_CELLS','PID_IL12_PATHWAY','BIOCARTA_CTLA4_PATHWAY','REACTOME_MTOR_SIGNALING','KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY')

# # 筛选名称包含 "METABOLISM" 的通路
# metabolism_pathways <- C2em@result[grep("METABOLISM", C2em@result$Description, ignore.case = TRUE), ]
disease_innerup_P <- c2c5upinner_pathways[c(grep('LONGECITY',c2c5upinner_pathways$c2c5upinner_pathways,ignore.case = TRUE),grep('CELLULAR',c2c5upinner_pathways$c2c5upinner_pathways,ignore.case = TRUE),grep('AGING',c2c5upinner_pathways$c2c5upinner_pathways,ignore.case = TRUE),grep('FATTY_ACID',c2c5upinner_pathways$c2c5upinner_pathways,ignore.case = TRUE),grep('PPARG',c2c5upinner_pathways$c2c5upinner_pathways,ignore.case = TRUE),grep('BILE_ACID',c2c5upinner_pathways$c2c5upinner_pathways,ignore.case = TRUE),grep('T_CELL',c2c5upinner_pathways$c2c5upinner_pathways,ignore.case = TRUE),grep('TH1',c2c5upinner_pathways$c2c5upinner_pathways,ignore.case = TRUE),grep('MTOR',c2c5upinner_pathways$c2c5upinner_pathways,ignore.case = TRUE),grep('T_HELPER',c2c5upinner_pathways$c2c5upinner_pathways,ignore.case = TRUE),grep('IL12',c2c5upinner_pathways$c2c5upinner_pathways,ignore.case = TRUE),grep('NFKB',c2c5upinner_pathways$c2c5upinner_pathways,ignore.case = TRUE)),]

disease_innerup_P <- as.data.frame(disease_innerup_P)
write.xlsx(disease_innerup_P,'disease_innerup_P.xlsx')

disease_innerdown_P <- c2c5downinner_pathways[c(grep('LONGECITY',c2c5downinner_pathways$c2c5downinner_pathways,ignore.case = TRUE),grep('CELLULAR',c2c5downinner_pathways$c2c5downinner_pathways,ignore.case = TRUE),grep('AGING',c2c5downinner_pathways$c2c5downinner_pathways,ignore.case = TRUE),grep('FATTY_ACID',c2c5downinner_pathways$c2c5downinner_pathways,ignore.case = TRUE),grep('PPARG',c2c5downinner_pathways$c2c5downinner_pathways,ignore.case = TRUE),grep('BILE_ACID',c2c5downinner_pathways$c2c5downinner_pathways,ignore.case = TRUE),grep('T_CELL',c2c5downinner_pathways$c2c5downinner_pathways,ignore.case = TRUE),grep('TH1',c2c5downinner_pathways$c2c5downinner_pathways,ignore.case = TRUE),grep('MTOR',c2c5downinner_pathways$c2c5downinner_pathways,ignore.case = TRUE),grep('T_HELPER',c2c5downinner_pathways$c2c5downinner_pathways,ignore.case = TRUE),grep('IL12',c2c5downinner_pathways$c2c5downinner_pathways,ignore.case = TRUE),grep('NFKB',c2c5downinner_pathways$c2c5downinner_pathways,ignore.case = TRUE)),]

disease_innerdown_P <- as.data.frame(disease_innerdown_P)
write.xlsx(disease_innerdown_P,'disease_innerdown_P.xlsx')


for (p in disease_innerup_P){
  print(p)
  gseaplot2(C2em, geneSetID = Aging_pathways, title = 'REACTOME_CELLULAR_SENESCENCE')
}
# C2
# MA_RAT_AGING_UP
# SCHOEN_NFKB_SIGNALING
# C5
# GOBP_FATTY_ACID_METABOLIC_PROCESS
# GOBP_T_CELL_DIFFERENTIATION_INVOLVED_IN_IMMUNE_RESPONSE

####plot####
gseaplot2(C2em, geneSetID = 'MA_RAT_AGING_UP', title = 'MA_RAT_AGING_UP')
gseaplot2(C2em, geneSetID = 'SCHOEN_NFKB_SIGNALING', title = 'SCHOEN_NFKB_SIGNALING')
gseaplot2(C5em, geneSetID = 'GOBP_FATTY_ACID_METABOLIC_PROCESS', title = 'GOBP_FATTY_ACID_METABOLIC_PROCESS')
gseaplot2(C5em, geneSetID = 'GOBP_T_CELL_DIFFERENTIATION_INVOLVED_IN_IMMUNE_RESPONSE', title = 'GOBP_T_CELL_DIFFERENTIATION_INVOLVED_IN_IMMUNE_RESPONSE')


library(ggplot2)
library(GseaVis)
c2id <- c('MA_RAT_AGING_UP','SCHOEN_NFKB_SIGNALING')
for (id in c2id){
  p<-gseaNb(object = C2em,
            geneSetID = id,
            # addGene = mygene,
            subPlot=3,
            addPval = T,
            pvalX = 0.75,pvalY = 0.6,
            pCol = 'black',
            pHjust = 0,
            rmPrefix=FALSE,
            addGene=TRUE,markTopgene=TRUE,rank.gene=TRUE,group=c('Model','Control'))
  ggsave(sprintf('c2_%s.pdf',id),p,height=4,width=6)
}
c5id <- c('GOBP_FATTY_ACID_METABOLIC_PROCESS','GOBP_T_CELL_DIFFERENTIATION_INVOLVED_IN_IMMUNE_RESPONSE')
for (id in c5id){
  p<-gseaNb(object = C5em,
            geneSetID = id,
            # addGene = mygene,
            subPlot=3,
            addPval = T,
            pvalX = 0.75,pvalY = 0.6,
            pCol = 'black',
            pHjust = 0,
            rmPrefix=FALSE,
            addGene=TRUE,markTopgene=TRUE,rank.gene=TRUE,group=c('Model','Control'))
  ggsave(sprintf('c5_%s.pdf',id),p,height=4,width=6)
}
          
         




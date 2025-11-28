library(omicsGMF)
library(dplyr)
library(QFeatures)
library(SingleCellExperiment)
library(stringr)
library(tibble)
library(factoextra)

## really important for Omics GMF
set.seed(100)

theme_custom_vis <- function(base_size = 12) {
  theme_bw(base_size = base_size) %+replace%
    theme(
      # leggend
      legend.title = element_text(size = rel(0.85), face = "plain"),
      legend.text = element_text(size = rel(0.70), face = "plain"),
      legend.key = element_rect(fill = "transparent", colour = NA),
      legend.key.size = unit(1.5, "lines"),
      legend.background = element_rect(fill = "transparent", colour = NA),
      strip.background = element_rect(fill = "#bababa", color = "black"),
      strip.text = element_text(size = rel(0.85), face = "bold", color = "#1b2944", margin = margin(5,0,5,0))
    )
}

theme_set(theme_custom_vis(base_size = 18))


path <- 'R://labs//Gevaert//Andrea//EBC_result//new_result_30062025/fix_omicsGMF_BatchDevice_Model'
pe <- readRDS('C:/Users/Andrea/workspace/CMB-EBC/data/LEgacy_eline/pe_object_no_filt_start.rda')

colData(pe) <- as.data.frame(colData(pe)) %>%
  mutate(Device = recode(Device, "new" = "TurboDECS", "old" = "EcoScreen")) %>%
  DataFrame()

rowData(pe[["precursor"]])$pep_per_prot <-
  left_join(rowData(pe[["precursor"]]) %>% as.data.frame %>% dplyr::select(Protein.Ids),
            rowData(pe[["precursor"]]) %>% as.data.frame %>% dplyr::group_by(Protein.Ids) %>%
              summarise(pep_per_prot = length(unique(Stripped.Sequence))))$pep_per_prot

pe <- filterFeatures(pe, ~ pep_per_prot > 2 )

pe <- filterFeatures(pe, ~ Proteotypic == 1)

pe <- filterFeatures(pe, ~  percHC >= 30 | percMPM >= 30 | percLC >= 30 | percAEX >= 30 )

pe <- logTransform(pe, base = 2, i = "precursor", 
                   name = "peptideLog")

pe <- infIsNA(pe, i='peptideLog')


## center.median
## diff.median 

pe <- normalize(pe, method = "center.median", i = "peptideLog", 
                name = "peptideNorm")
pe <- aggregateFeatures(pe, i = "peptideNorm",
                        fcol = "Protein.Ids",
                        name = "proteinB",
                        fun = MsCoreUtils::medianPolish,
                        # slower but better than medianPolish robustSummary
                        na.rm = TRUE)
write.csv(assay(pe[['precursor']]), file.path(path, 'raw_precursor.csv') )
write.csv(assay(pe[['peptideLog']]), file.path(path, 'raw_precursor_log.csv') )
write.csv(assay(pe[['peptideNorm']]), file.path(path, 'normalized_precursor.csv') )


res.pca_quant <-
  pe[["peptideNorm"]] %>%
  filterNA() %>%
  assay() %>%
  t() %>%
  prcomp()


pca <- fviz_pca_ind(res.pca_quant, label="none")
pca_1 <- ggplot(cbind(pca$data,colData(pe)[,c("Group","Batch","Age","Device")]),
                aes(x=x,y=y,col= Device ,shape = Batch) ) + geom_point(size=3) + 
  scale_color_manual(values=c("#6d0286", "#35B779FF","orange",'cyan','blue'))  +
  ggtitle( 'PCA - Precursor Normalized - No Correction')

pca_1

pdf(file = file.path(path, 'PeptideNorm_noCorr_DeviceBatch.pdf'), paper= "a4")
pca_1
dev.off()
## ----- how to do it right now ----------##

set.seed(100)
sce <- pe[["peptideNorm"]]
names(assays(sce)) <- "logintensities"

colData(sce) <-  colData(pe)
saveRDS(object = sce, file = file.path(path,"sce_pepNormbase.RDS"))

X <- model.matrix(~1 + Batch + Device  , colData(sce) )
family <- gaussian()

sgd_best <-  calculateGMF( sce , X = X,family = family ,exprs_values='logintensities',
                               control.alg = list(tol = 0.001, maxiter = 10000),
                               ncomponents = 6)

pca_best <- ggplot(cbind(sgd_best,colData(pe)[,c("Group","Batch","Age","Device")]),
                        aes(x=PC1,y=PC2,col= Device ,shape = Batch) ) + geom_point(size=3) + 
  scale_color_manual(values=c("#6d0286", "#35B779FF","orange",'cyan','blue'))  +
  ggtitle( 'PCA Precursor Corrected by ~ Batch + Device')

pca_best



pdf(file = file.path(path, 'PeptideNorm_omicsGMF_DeviceBatch_6comp.pdf'), paper= "a4")
pca_best
dev.off()


# assay(sce)(['AADDTWEPFASGK2',c("A189","A209","E149","A117","B032")])
base_no_batchcorr <- sce
assay(base_no_batchcorr) <- imputeGMF(assay(base_no_batchcorr), sgd_best) 
rowData(base_no_batchcorr) <- rowData(sce)


pe <- addAssay(pe, base_no_batchcorr, name = "sgd_batchBest")
pe <- aggregateFeatures(pe, i = "sgd_batchBest",
                        fcol = "Protein.Ids",
                        name = "proteinGMF_B",
                        fun = MsCoreUtils::medianPolish,
                        # slower but better than medianPolish robustSummary
                        na.rm = TRUE)

write.csv(assay(pe[['sgd_batchBest']]), file.path(path, 'imputed_corrected_precursor.csv') )
write.csv(assay(pe[['proteinGMF_B']]), file.path(path, 'imputed_corrected_proteins.csv') )
saveRDS(pe, file.path(path,'pe_6compBatchDevice_2374prec.Rds'))


## keep scale=FALSE in pca. data is already normalized
b_a <-
  pe[["proteinGMF_B"]] %>%
  filterNA() %>%
  assay() %>%
  t() %>%
  prcomp()

b_b <-
  pe[["proteinB"]] %>%
  filterNA() %>%
  assay() %>%
  t() %>%
  prcomp()


### 
#hist( assay(pe_m[['proteinGMF_B']])[,'A218'], probability = TRUE, col = "lightgray",
#       border = "white", main = "Histogram No corrected", xlab = "Value")
#lines(density( assay(pe_m[['sgd_batchBest']])[,'A218']), col = "blue", lwd = 2)

#hist( assay(pe[['proteinGMF_Bcor']])[,'A218'], probability = TRUE, 
#col = "lightgray", border = "white", main = "Corrected by XB comp", xlab = "Value")
#lines(density( assay(pe[['proteinGMF_Bcor']])[,'A218']), col = "blue", lwd = 2)

###
## manual; code
variance <- b_a$sdev^2
percent_variance <- (variance / sum(variance)) * 100
pc_labels <- paste0(
  names(percent_variance),
  " = ",
  round(percent_variance, 2), 
  "%"
)
perc_df <- data.frame(
  PC = paste0("PC", seq_along(percent_variance)),
  PercentVariance = round(percent_variance, 2)
)
library(knitr)

perc_df %>% head(10) %>% kable(format = "simple")

pca_b_a <- fviz_pca_ind(b_a, label="none")
pca_b_b <- fviz_pca_ind(b_b, label="none")
# explained variance in the other component

comp_explained <- fviz_eig(b_a, ncp = 10)
comp_explained <- comp_explained + theme_custom_vis(base_size = 18) + ggtitle('Explained variance in % for PCA components')


 ggsave(  file.path(path, "ExplainedVar_Protein_omicsGMF_DeviceBatch_6comp_.svg"),
           plot = pca_prot_a, width = 10, height = 6, 
dpi= 300,  units='in')

pca_prot_a <- ggplot(cbind(pca_b_a$data,colData(pe)[,c("Group","Batch","Age","Device")]),
                    aes(x=x,y=y,col= Device ,shape = Batch  ) ) + geom_point(size=3) + 
  scale_color_manual(values=c("#6d0286", "#35B779FF","orange",'cyan'))  +
  ylab("PC2 (15.99%)")  +
  xlab("PC1 (23.31%)")  +
  ggtitle( '')
pca_prot_a

 ggsave(  file.path(path, "PCA_Protein_omicsGMF_DeviceBatch_6comp_.svg"),
           plot = pca_prot_a, width = 10, height = 6, 
dpi= 300,  units='in')


# PC2 - PC3 
pca_b_23 <- fviz_pca_ind(b_a, label="none",axes = c(2, 3))
pca_prot_23 <- ggplot(cbind(pca_b_23$data,colData(pe)[,c("Group","Batch","Age","Device")]),
                    aes(x=x,y=y,col= Device ,shape = Batch  ) ) + geom_point(size=3) + 
  scale_color_manual(values=c("#6d0286", "#35B779FF","orange",'cyan'))  +
  ylab("PC2 (15.99%)")  +
  xlab("PC3 (15.05%)")  +
  ggtitle('')
pca_prot_23
 ggsave(  file.path(path, "PCA_23comp_Protein_omicsGMF_DeviceBatch_6comp_.svg"),
           plot = pca_prot_23, width = 10, height = 6, 
dpi= 300,  units='in')



# PC3 - PC4 
pca_b_34 <- fviz_pca_ind(b_a, label="none",axes = c(3, 4))
pca_prot_34 <- ggplot(cbind(pca_b_34$data,colData(pe)[,c("Group","Batch","Age","Device")]),
                    aes(x=x,y=y,col= Device ,shape = Batch  ) ) + geom_point(size=3) + 
  scale_color_manual(values=c("#6d0286", "#35B779FF","orange",'cyan'))  +
  ylab("PC3 (15.05%)")  +
  xlab("PC4 (8.08%)")  +
  ggtitle('')
pca_prot_34
 ggsave(  file.path(path, "PCA_34comp_Protein_omicsGMF_DeviceBatch_6comp_.svg"),
           plot = pca_prot_34, width = 10, height = 6, 
dpi= 300,  units='in')

pdf(file =  file.path( path, 'Protein_omicsGMF_DeviceBatch_6comp.pdf'), paper= "a4")
pca_prot_a
dev.off()

pca_prot_b <- ggplot(cbind(pca_b_b$data,colData(pe)[,c("Group","Batch","Age","Device")]),
                     aes(x=x,y=y,col= Device ,shape = Batch) ) + geom_point(size=3) + 
  ylab("PC2")  +
  xlab("PC1")  +
  scale_color_manual(values=c("#6d0286", "#35B779FF","orange",'cyan'))  +
  ggtitle( 'PCA Summarized Protein without correction & imputation')
pca_prot_b


pdf(file =  file.path( path,  'Protein_Base_No_imp_correction.pdf'), paper= "a4")
pca_prot_b
dev.off()

### adding tumor  vs . AEx and HC 

colData(pe) <- colData(pe) %>% as.data.frame() %>% mutate( GrpAggr = case_when( Group == 'LC' | Group == "MPM"  ~ 'LC_MPM', .default = Group ) )
colData(pe)$GrpAggr <- as.factor(pe$GrpAggr) 

library(msqrob2)
pe <- msqrob (object = pe, i = "proteinGMF_B", 
              formula =  ~ -1 + GrpAggr +   Smokestat + Packyears + Gender + BMI + Batch + Device  , ridge=FALSE,
              overwrite = T
)

coef <- rowData(pe[["proteinGMF_B"]])$msqrobModels[[1]] %>% getCoef %>% names
## comparison aggregated
L <- makeContrast(c(
  "GrpAggrLC_MPM - GrpAggrAEX =0",'GrpAggrLC_MPM - GrpAggrHC =0'),parameterNames = coef )

# comparison single
# L <- makeContrast(c(
#   "GroupMPM - GroupAEX =0",'GroupMPM - GroupHC =0',
#   'GroupLC - GroupAEX =0', 'GroupLC - GroupHC =0',
#   'GroupAEX - GroupHC =0', 'GroupMPM - GroupLC =0'),parameterNames = coef )

#getContrast(rowData(pe[["proteinSGDb"]])$msqrobModels[[1]], L)

pe <- hypothesisTest(pe, i="proteinGMF_B", L, 
                     overwrite=TRUE )

## comparison aggregated

comparisons <- c("GrpAggrLC_MPM - GrpAggrAEX",'GrpAggrLC_MPM - GrpAggrHC')

# comparison single
# comparisons <- c("GroupMPM - GroupAEX",'GroupMPM - GroupHC', 
#                  'GroupLC - GroupAEX',
#                  'GroupLC - GroupHC','GroupAEX - GroupHC','GroupMPM - GroupLC')
library(logger)
library(ggrepel)
source('utils_function.R')
res_Mod <-  lapply(comparisons, dep_volcano, data= pe,  level='proteinGMF_B', FC_thr = 1)
names(res_Mod) <- comparisons

res_Mod$`GroupMPM - GroupAEX`$volcano


res_Mod$`GroupLC - GroupAEX`$volcano2file
res_Mod$`GroupMPM - GroupAEX`$volcano2file
res_Mod$`GroupHC - GroupAEX`$volcano2file
res_Mod$`GroupMPM - GroupLC`$volcano


res_Mod$`GrpAggrLC_MPM - GrpAggrAEX`$volcano2file
res_Mod$`GrpAggrLC_MPM - GrpAggrAEX`$volcano

res_Mod$`GrpAggrLC_MPM - GrpAggrH`$volcano2file
res_Mod$`GrpAggrLC_MPM - GrpAggrH`$volcano

## aggregated analysis export
write.csv2(res_Mod$`GrpAggrLC_MPM - GrpAggrAEX`$toptable , file.path(path, 'TopTable_MPMLC_AEX.csv'))
write.csv2(res_Mod$`GrpAggrLC_MPM - GrpAggrHC`$toptable , file.path(path, 'TopTable_MPMLC_HC.csv'))

## single analsis 
write.csv2(res_Mod$`GroupMPM - GroupAEX`$toptable , file.path(path, 'TopTable_MPM_AEX.csv'))
write.csv2(res_Mod$`GroupMPM - GroupHC`$toptable , file.path(path, 'TopTable_MPM_HC.csv'))
write.csv2(res_Mod$`GroupLC - GroupAEX`$toptable  ,file.path(path, 'TopTable_LC_AEX.csv') )
write.csv2(res_Mod$`GroupLC - GroupHC`$toptable  ,file.path(path, 'TopTable_LC_HC.csv') )
write.csv2(res_Mod$`GroupAEX - GroupHC`$toptable , file.path(path, 'TopTable_AEX_HC.csv'))
write.csv2(res_Mod$`GroupMPM - GroupLC`$toptable , file.path(path, 'TopTable_MPM_LC.csv'))



for (cmp in comparisons){
  
  res_Mod$cmp$volcano2file
  
  ggsave(  file.path(path, paste(gsub(' ','', cmp), "_volcano_Aggregate.svg")),
           plot = res_Mod[[cmp]][["volcano2file"]], width = 6, height = 6, 
dpi= 300,  units='in')
}


## aggregagated analysis .
saveRDS(res_Mod, file.path(path,'DE_batch_dev_corr_aggrAnalysis.Rds') )


saveRDS(res_Mod, file.path(path,'DE_batch_dev_corr.Rds') )
## read result . 
res_Mod <- readRDS( file.path(path,'DE_batch_dev_corr.Rds') )


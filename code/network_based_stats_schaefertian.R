
install.packages("NBR")
install.packages("parallel")
install.packages("pheatmap")
install.packages("pals")
library(NBR)
library(parallel)
library(pheatmap)
library(pals)
library(parallel)
library(pheatmap)
library(pals)

setwd("/SchaeferTian-corrmatrix/")

path_corr = list.files(pattern = "*.csv")
subject = as.numeric(gsub("SchaeferTianAtlassub-","",gsub("_correlation_matrix.csv","",path_corr)))

# Load node labels
brain_labs = c(1:432) #read.table(file = "/Schaefer2018_400Parcels_17Networks_order_Tian_Subcortex_S2_label.txt", sep="\t", header=T)$nom_l 

# Load phenotypic info
phen <- read.csv("/connectome_demographics.csv")
phen$group = as.factor(phen$group)

# ---- missing people - mri scan with demopgrahic information 
communPart = intersect(subject, phen$rid)
path_corr = path_corr[subject %in% communPart == TRUE]
length(communPart) # number of participant with MRI and demographics data
length(subject) # number of participants with MRI data
length(phen$rid) # number of participants with demographics data 
phen = phen[phen$rid %in% communPart ==TRUE,]
length(path_corr) == dim(phen)[1]

# exclude subjects with high motion (mean FD >0.55 mm) #sub-108,070,010,025,094


# Load 3D array
cmx = array(NA, dim=c(length(brain_labs), 
                      length(brain_labs),
                      length(path_corr))) 
dim(cmx)
for (i in 1:length(path_corr)){
  cmx[,,i] = as.matrix(read.csv(path_corr[i])[,-1])
  message(path_corr[i])
}

# average matrix of all participants
avg_mx <- apply(cmx, 1:2, mean)
pheatmap(avg_mx,
         treeheight_row = 0, treeheight_col = 0,
         cluster_rows=FALSE, cluster_cols=FALSE, 
         breaks = seq(-1, 1, by = 0.01),
         color = coolwarm(n=200))

# Network based stats test
set.seed(19930413)
before <- Sys.time()
nbr_group001 <- nbr_lm_aov(net = cmx, 
                        nnodes = dim(cmx)[1],
                        idata = phen,
                        mod = "~ group + age + sex + educ",
                        thrP = 0.05,
                        nperm = 1000, 
                        cores = detectCores())
after <- Sys.time()
show(after-before)
length(nbr_group05)


#save result
saveRDS(nbr_group05, file = "nbs_result05_SchaeferTian_covars_excludefd_230524.rds")

#visualization of results

library(pheatmap)
library(pals)
library(lattice)
library(pheatmap)
library(pals)
library(psych)
library(ggplot2)

source("") #set source here


#===============================================================================
#===============================================================================
#results_01 = readRDS("/Users/jocelynricard/Desktop/SUDMEX/sudmexconn/results/SchaeferTian/nbs_result01_SchaeferTian_220815.rds") 
#length(results_01)
#results_001 = readRDS("/Users/jocelynricard/Desktop/SUDMEX/sudmexconn/results/SchaeferTian/nbs_result001_SchaeferTian_220816.rds")
#length(results_001)
results_05 = readRDS("/Users/jocelynricard/Desktop/sudmex_230823/data/nbs_result05_covars_excludeFD_SchaeferTian_230824.rds")
length(results_05)
#results_01=readRDS("/Users/jocelynricard/Desktop/sudmex_230823/data/nbs_result01_covars_excludeFD_SchaeferTian_230828.rds")
#length(results_01)
#results_001=readRDS("/Users/jocelynricard/Desktop/sudmex_230823/data/nbs_result001_covars_excludeFD_SchaeferTian_230927.rds")
#length(results_001)

show(results_05$fwe$group)
#show(results_01$fwe$group)
#show(results_001$fwe$group)

# Load node labels
brain_labs = c(1:432) #read.table(file = "/gpfs/milgram/project/holmes/jar293/sudmexconn/data/derivatives/mri_MNI_processed_sid/rest_TimeSerie_Volumes/SchaeferTianAtlas/corr_matrix_sch/Schaefer2018_400Parcels_17Networks_order_Tian_Subcortex_S2_label.txt", sep="\t", header=T)$nom_l 

# Load phenotypic info
phen <- read.csv(".../connectome_demographics.csv")
phen$group = as.factor(phen$group)

setwd(".../SchaeferTian-corrmatrix/")
path_corr = list.files(pattern = "*.csv")
subject = as.numeric(gsub("SchaeferTianAtlassub-","",gsub("_correlation_matrix.csv","",path_corr)))

#subject = as.numeric(gsub("SchaeferTianAtlassub-","",gsub("_correlation_matrix.csv","",path_corr)))
communPart = intersect(subject, phen$rid)
path_corr = path_corr[subject %in% communPart == TRUE]
length(communPart) # number of participant with MRI and demographics data
length(subject) # number of participants with MRI data
length(phen$rid) # number of participants with demographics data 
phen = phen[phen$rid %in% communPart ==TRUE,] #communPart = index of all included participants with both demos and MRI
length(path_corr) == dim(phen)[1]

#-------------------------------------------------------------------------------

#Load 3D array
cmx = array(NA, dim=c(length(brain_labs), 
                      length(brain_labs),
                      length(path_corr))) 
dim(cmx)
for (i in 1:length(path_corr)){
  cmx[,,i] = as.matrix(read.csv(path_corr[i])[,-1])
}
dimnames(cmx) = list(brain_labs, brain_labs, path_corr)



# average matrix of all participants
avg_mx <- apply(cmx, 1:2, mean)
pheatmap(avg_mx,
         treeheight_row = 0, treeheight_col = 0,
         cluster_rows=FALSE, cluster_cols=FALSE, 
         breaks = seq(-1, 1, by = 0.01),
         color = coolwarm(n=200))


# visualization of the results:
# 0.05
edge_mat_05 <- array(0, dim(cmx)[1:2])
edge05 = rowSums(edge_mat_05)
edge_mat_05[results_05$components$group[,2:3]] <- 1
edge_mat_05 = as.matrix(Matrix::forceSymmetric(edge_mat_05))
isSymmetric.matrix(edge_mat_05)
levelplot(edge_mat_05, col.regions = rev(heat.colors(100)),
          main = "Component", ylab = "ROI", xlab = "ROI")
show(results_05$fwe$group)

# 0.01
#edge_mat_01 <- array(0, dim(cmx)[1:2])
#edge_mat_01[results_01$components$group[,2:3]] <- 1
#edge_mat_01 = as.matrix(Matrix::forceSymmetric(edge_mat_01))
#isSymmetric.matrix(edge_mat_01)
#levelplot(edge_mat_01, col.regions = rev(heat.colors(100)),
#          main = "Component", ylab = "ROI", xlab = "ROI")
#show(results_01$fwe$group)

# 0.001
#edge_mat_001 <- array(0, dim(cmx)[1:2])
#edge_mat_001[results_001$components$group[,2:3]] <- 1
#edge_mat_001 = as.matrix(Matrix::forceSymmetric(edge_mat_001))
#isSymmetric.matrix(edge_mat_001)
#levelplot(edge_mat_001, col.regions = rev(heat.colors(100)),
#          main = "Component", ylab = "ROI", xlab = "ROI")
#show(results_001$fwe$group)

#-------------------------------------------------------------------------------
# post analysis
# rowSums(results_05)
# res = data.frame(roi = brain_labs, nb_edges = rowSums(results_05))
# res = read.csv(".../SchaeferTian-corrmatrix/Schaefer400_17networks_FSLMNI_2m_TianSubCortex_S2_3T_fullBrain.csv")

# z-fisher transform the individuals matrices: (want to plot the avg matrix for two groups - but you need to fischer transform)
cmx_z = cmx
for (i in 1:dim(cmx_z)[3]){
  cmx_z[,,i] = fisherz(cmx_z[,,i])
}
# Plot the fisher average matrix for fun
avg_Z_mx <- apply(cmx_z, 1:2, mean)
pheatmap(avg_Z_mx,
         treeheight_row = 0, treeheight_col = 0,
         cluster_rows=FALSE, cluster_cols=FALSE, 
         breaks = seq(-round(max(abs(max(avg_Z_mx[upper.tri(avg_Z_mx)])),
                                 abs(min(avg_Z_mx[upper.tri(avg_Z_mx)]))),2), 
                      round(max(abs(max(avg_Z_mx[upper.tri(avg_Z_mx)])),
                                abs(min(avg_Z_mx[upper.tri(avg_Z_mx)]))),2),
                      by = 0.01),
         color = coolwarm(n=200))
# compute mean correlation matrices for the 2 group
avg_gr2_cocaine <- apply(cmx_z[,,phen$group == 2], 1:2, mean)
avg_gr1_control <- apply(cmx_z[,,phen$group == 1], 1:2, mean)
# compute the group's difference:
avg_corr = fisherz2r(avg_gr1_control - avg_gr2_cocaine)
diag(avg_corr) = 0


# mask it by the significant edges: 
##0.05
avg_corr_masked = avg_corr * edge_mat_05 #edge_mat_05

##0.01
#avg_corr_masked = avg_corr * edge_mat_01

##0.001
#avg_corr_masked = avg_corr * edge_mat_001

# visualization of the group's difference:
max_diff = max(avg_corr_masked[upper.tri(avg_corr_masked)])
min_diff = min(avg_corr_masked[upper.tri(avg_corr_masked)])
abs_max_range = max(abs(max_diff), abs(min_diff))
breaksList = seq(-round(abs_max_range,2),
                 round(abs_max_range,2), 
                 by = 0.01)
color = rev(colorRampPalette(brewer.rdbu(n=200))(length(breaksList)))
pheatmap(avg_corr_masked,
         treeheight_row = 0, 
         treeheight_col = 0,
         cluster_rows=FALSE, 
         cluster_cols=FALSE, 
         color = color,
         breaks = breaksList)
# nicer visualization

atlas_Yeo = read.csv("/Schaefer400_17networks_FSLMNI_2m_TianSubCortex_S2_3T_fullBrain.csv")


 # No matter the sign of the difference 
sum(edge_mat_05[upper.tri(edge_mat_05)])

ordered_regions = atlas_Yeo[with(atlas_Yeo, order(roi_number, side)), ]$old_order 
ordered_networkNumber = atlas_Yeo[with(atlas_Yeo, order(roi_number, side)), ]$Network_Split_number
ordered_label = unique(atlas_Yeo[with(atlas_Yeo, order(roi_number, side)), ]$Network_split)
pce <- plotClassifiedEdges(adj = as.matrix(edge_mat_05[ordered_regions,
                                                       ordered_regions]), 
                           ids = ordered_networkNumber,
                           labels = ordered_label)
rowSums(as.matrix(edge_mat_05[ordered_regions,
                              ordered_regions]))
pheatmap(edge_mat_05[ordered_regions,
                     ordered_regions],
         treeheight_row = 0, 
         treeheight_col = 0,
         cluster_rows=FALSE, 
         cluster_cols=FALSE, 
         color = coolwarm(n=200))

# makeNetworkMatrix2(pce[[3]], pce[[1]])
# positive difference 
avg_masked_Positive = avg_corr_masked
avg_masked_Positive[avg_masked_Positive<0] = 0
avg_masked_Positive[avg_masked_Positive>0] = 1
write.csv(avg_masked_Positive, "avg_masked_Positive_covars_fd_231024.csv")
pce_posi <- plotClassifiedEdges(as.matrix(avg_masked_Positive[ordered_regions,
                                                              ordered_regions]), 
                                ordered_networkNumber,
                                labels = ordered_label)

degree05pos = rowSums(avg_masked_Positive)
#degree01pos = rowSums(avg_masked_Positive)
#degree001pos = rowSums(avg_masked_Positive)


# Assuming avg_masked_Positive is your binary matrix
upper_tri_matrix <- avg_masked_Positive[upper.tri(avg_masked_Positive)]
upper_tri_sum_positive <- sum(upper_tri_matrix)
# Print the result
cat("Sum of the upper triangle:", upper_tri_sum_positive, "\n")



# negative difference 
avg_masked_Negative = avg_corr_masked
avg_masked_Negative[avg_masked_Negative>0] = 0
avg_masked_Negative[avg_masked_Negative<0] = 1
write.csv(avg_masked_Negative, "avg_masked_Negative_covars_fd_231024.csv")
pce_neg <- plotClassifiedEdges(as.matrix(avg_masked_Negative[ordered_regions,
                                                             ordered_regions]), 
                               ordered_networkNumber,
                               labels = ordered_label)

degree05neg = rowSums(avg_masked_Negative) # binary output of total matrix 0 or 1.
#degree01neg = rowSums(avg_masked_Negative)
#degree001neg = rowSums(avg_masked_Negative)

# Assuming avg_masked_Negative is your binary matrix
upper_tri_matrix_neg <- avg_masked_Negative[upper.tri(avg_masked_Negative)]
upper_tri_sum_neg <- sum(upper_tri_matrix_neg)

# Print the result
cat("Sum of the upper triangle (Negative matrix):", upper_tri_sum_neg, "\n")


# Write out file to use to make figures - binary 0 or 1 (yes or no significant edge here - total number of pos or neg signf. edges in each row)
#write.csv(x = degree001neg ,"degree001neg.csv", row.names = FALSE)
#write.csv(x = degree001pos ,"degree001pos.csv", row.names = FALSE)
#write.csv(x = degree01neg ,"degree01neg.csv", row.names = FALSE)
#write.csv(x = degree01pos ,"degree01pos.csv", row.names = FALSE)
write.csv(x = degree05neg ,"degree05neg_230825.csv", row.names = FALSE)
write.csv(x = degree05pos ,"degree05pos_230825.csv", row.names = FALSE)

# negative and positive
makeNetworkMatrix2(pce_neg[[3]], pce_posi[[3]]) #normalized matrix
makeNetworkMatrix2(pce_neg[[1]], pce_posi[[1]]) #non normalized matrix (does not account for the number of edges in each region but rather treats everything the same).


## -------------------------------------------------------
## -------------------------------------------------------

# Compute and print total significant edges for avg_masked_Positive
total_positive_edges <- sum(avg_masked_Positive[upper.tri(avg_masked_Positive)])
cat("Total significant positive edges across all networks:", total_positive_edges, "\n")

# Compute and print total significant edges for avg_masked_Negative
total_negative_edges <- sum(avg_masked_Negative[upper.tri(avg_masked_Negative)])
cat("Total significant negative edges across all networks:", total_negative_edges, "\n")

## -------------------------------------------------------
## -------------------------------------------------------


# tells you which regions have the highest edge count implication

# Set the names of degree vectors using atlas_Yeo$roi_name
names(degree05neg) <- atlas_Yeo$roi_name
names(degree05pos) <- atlas_Yeo$roi_name

N <- 5
cat("Top", N, "regions with highest negative edges:\n")

# Sorting the degree vectors while preserving names
sorted_neg <- sort(degree05neg, decreasing = TRUE, index.return = TRUE)
for (i in 1:N) {
  cat(names(degree05neg)[sorted_neg$ix[i]], "(Network:", atlas_Yeo$Network_name[sorted_neg$ix[i]], "):", sorted_neg$x[i], "\n")
}

cat("\nTop", N, "regions with highest positive edges:\n")
sorted_pos <- sort(degree05pos, decreasing = TRUE, index.return = TRUE)
for (i in 1:N) {
  cat(names(degree05pos)[sorted_pos$ix[i]], "(Network:", atlas_Yeo$Network_name[sorted_pos$ix[i]], "):", sorted_pos$x[i], "\n")
}

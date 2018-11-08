#Heatmaps with pheatmap and heatmaply!
#For U of T coders
#Simulated data created by Sahir Bhatnagar.

#possible data pre-processing - normalization - quantile, median, etc., log transform
#not necessary here - we have log fold change data that has already been normalized

#Calculating your distance matrix (see dist objects):
#compute how similar or different you values are
#parametric - distance measures based on Pearson correlation 
#non parametric - spearman rank - replace by ranks and calculate correlation, Kendall's - relative ordering
#euclidean - shortest distance between values (has to be normalized), takes magnitude into account
#city block/Manhattan - sum of distances along each dimension
#distance 1-correlation - of all pairs of items to be clustered

#Cluster your samples (see hclust objects):
#hierarchical, organizes into a tree structure based on similarity - short branches if similar and longer branches as similarity decreases
#repeated cycles where the 2 closest remaining items (smallest distance) get joined by a branch with the length of the branch reflecting the distance between them, the distance between this item and all other remaining items are computed until only one object remains
#single linkage clustering - distance between 2 items is the minimum of all pairwise distances between items contained in x and y - fast b/c no other calculations need to be performed once you have your distance matrix
#complete linkage is the maximum of all paiwise distances between x and y 
#average linkage - mean of all pairwise distances between items contained in x and y
#k-means organize into clusters (self-chosen number) - items are randomly assigned to a cluster - the mean vector fo rall items in each hcluster is computed, items are reassigned to the cluster whose center is closest to them - random starting points so will not always get the same answer, number of trial done to deal with the randomness
#self organizing maps



library(pheatmap)
library(heatmaply)

genes <- read.csv("heatmap_data.csv", header = T, row.names=1)
annotation_col <- read.csv("annotation_col.csv", header = T, row.names=1)
annotation_row <- read.csv("annotation_row.csv", header = T, row.names=1)

str(genes)
str(annotation_col)
str(annotation_row)


#Plotting with pheatmap!

pheatmap(genes)
#change font
pheatmap(genes, fontsize = 6) #default is clustering rows and columns

#cluster by gene - groups of similar genes
pheatmap(genes, fontsize = 6, cluster_rows = T, cluster_cols = F) 
#only now have dendrograms for the x-axis

#cluster by patient - groups of similar patients
pheatmap(t(genes), fontsize = 6, cluster_rows = F, cluster_cols = T)
#only have dendrogram on the y-axis

#usually order by both
pheatmap(genes, fontsize = 6, cluster_rows = T, cluster_cols = T)

#seeing some patterns emerge - but what do they mean? Great time to add annotation to our plot
pheatmap(genes, fontsize = 6, cluster_rows = T, cluster_cols = T, annotation_row = annotation_row)
#add to row first, see that genes are clustering according to the pathways they belong to

pheatmap(genes, fontsize = 6, cluster_rows = T, cluster_cols = T, annotation_row = annotation_row,
         annotation_col = annotation_col)
#now have information about the drug and condition 

#title
pheatmap(genes, fontsize = 6, cluster_rows = T, cluster_cols = T, annotation_row = annotation_row,
         annotation_col = annotation_col, main = "Gene Expression for Drug X Trial for B-Cell vs T-Cell Patients")

#changing relative sizes of text
pheatmap(t(genes), fontsize = 8, fontsize_row = 6, fontsize_col = 6, cluster_rows = T, cluster_cols = T, annotation_row = annotation_row, annotation_col = annotation_col, main = "Gene Expression for Drug X Trial for B-Cell vs T-Cell Patients")

#take a smaller subset 
pheatmap(genes[c(1:5, 55:60), c(1:5, 20:35, 55:60)], fontsize = 8, fontsize_row = 6, fontsize_col = 6, cluster_rows = T, cluster_cols = F, annotation_row = annotation_row, annotation_col = annotation_col, main = "Gene Expression for Drug X Trial for B-Cell vs T-Cell Patients")

#display numbers
pheatmap(genes[c(1:5, 55:60), c(1:5, 20:35, 55:60)], fontsize = 8, fontsize_row = 6, fontsize_col = 6, cluster_rows = T, cluster_cols = T, annotation_row = annotation_row, annotation_col = annotation_col, main = "Gene Expression for Drug X Trial for B-Cell vs T-Cell Patients", display_numbers = TRUE)

#font size in cells
pheatmap(genes[c(1:5, 55:60), c(1:5, 20:35, 55:60)], fontsize = 10, fontsize_row = 8, fontsize_col = 8, cluster_rows = T, cluster_cols = T, annotation_row = annotation_row, annotation_col = annotation_col, main = "Gene Expression for Drug X Trial for B-Cell vs T-Cell Patients", display_numbers = TRUE, fontsize_number = 8)

#gap between exposure and treatment
pheatmap(genes[c(1:5, 55:60), c(1:5, 20:35, 55:60)], fontsize = 10, fontsize_row = 8, fontsize_col = 8, cluster_rows = T, cluster_cols = F, annotation_row = annotation_row, annotation_col = annotation_col, main = "Gene Expression for Drug X Trial for B-Cell vs T-Cell Patients", display_numbers = TRUE, fontsize_number = 8, gaps_col = 5)

#get rid of legends
pheatmap(genes[c(1:5, 55:60), c(1:5, 20:35, 55:60)], fontsize = 10, fontsize_row = 8, fontsize_col = 8, cluster_rows = T, cluster_cols = F, annotation_row = annotation_row, annotation_col = annotation_col, main = "Gene Expression for Drug X Trial for B-Cell vs T-Cell Patients", display_numbers = TRUE, fontsize_number = 8, gaps_col = 5, legend = FALSE, annotation_legend = FALSE)

#get rid of dendrogram
pheatmap(genes[c(1:5, 55:60), c(1:5, 20:35, 55:60)], fontsize = 10, fontsize_row = 8, fontsize_col = 8, cluster_rows = T, cluster_cols = F, annotation_row = annotation_row, annotation_col = annotation_col, main = "Gene Expression for Drug X Trial for B-Cell vs T-Cell Patients", display_numbers = TRUE, fontsize_number = 8, gaps_col = 5, legend = FALSE, annotation_legend = FALSE, treeheight_row = 0)

#preserving the order of your clustered matrix
tmat <- genes[c(1:5, 55:60), c(1:5, 20:35, 55:60)]

pheatmap(tmat, fontsize = 10, fontsize_row = 8, fontsize_col = 8, cluster_rows = T, cluster_cols = T, annotation_row = annotation_row, annotation_col = annotation_col, main = "Gene Expression for Drug X Trial for B-Cell vs T-Cell Patients", display_numbers = TRUE, fontsize_number = 8, gaps_col = 5, legend = FALSE, annotation_legend = FALSE, treeheight_row = 0, treeheight_col = 0)

H <- pheatmap(tmat, fontsize = 10, fontsize_row = 8, fontsize_col = 8, cluster_rows = T, cluster_cols = T, annotation_row = annotation_row, annotation_col = annotation_col, main = "Gene Expression for Drug X Trial for B-Cell vs T-Cell Patients", display_numbers = TRUE, fontsize_number = 8, gaps_col = 5, legend = FALSE, annotation_legend = FALSE, treeheight_row = 0, treeheight_col = 0)

H$tree_row$order
H$tree_col$order

pheatmap(tmat[H$tree_row$order,H$tree_col$order],fontsize = 10, fontsize_row = 8, fontsize_col = 8, cluster_rows = F,cluster_cols = F, annotation_row = annotation_row, annotation_col = annotation_col, main = "Gene Expression for Drug X Trial for B-Cell vs T-Cell Patients", display_numbers = TRUE, fontsize_number = 8, gaps_col = 5, legend = FALSE, annotation_legend = FALSE, treeheight_row = 0, treeheight_col = 0)

#can take this and extract as a csv by converting it into a data frame

#color
library(viridis)

pheatmap(tmat[H$tree_row$order,H$tree_col$order],fontsize = 10, fontsize_row = 8, fontsize_col = 8, cluster_rows = F,cluster_cols = F, annotation_row = annotation_row, annotation_col = annotation_col, main = "Gene Expression for Drug X Trial for B-Cell vs T-Cell Patients", display_numbers = TRUE, fontsize_number = 8, gaps_col = 5, legend = FALSE, annotation_legend = FALSE, treeheight_row = 0, treeheight_col = 0, col = viridis_pal(option = "cividis")(6))

pal <- c("lightblue", "cornflowerblue", "navyblue") 

pheatmap(tmat[H$tree_row$order,H$tree_col$order],fontsize = 10, fontsize_row = 8, fontsize_col = 8, cluster_rows = F,cluster_cols = F, annotation_row = annotation_row, annotation_col = annotation_col, main = "Gene Expression for Drug X Trial for B-Cell vs T-Cell Patients", display_numbers = TRUE, fontsize_number = 8, gaps_col = 5, legend = FALSE, annotation_legend = FALSE, treeheight_row = 0, treeheight_col = 0, col = pal)



#Interactive version with heatmaply!

heatmaply(tmat, main = "Gene Expression for Drug X Trial for B-Cell vs T-Cell Patients")

#color dendrograms by group
heatmaply(tmat, main = "Gene Expression for Drug X Trial for B-Cell vs T-Cell Patients", k_col=3, k_row=4)

#correlation
heatmaply(cor(tmat), main = "Gene Expression for Drug X Trial for B-Cell vs T-Cell Patients")

#scale - if all variables come from a normal distribution, scaling - subtract the mean and divide by standard deviation, would bring them close to the standard normal distribution
heatmaply(tmat,scale="column", main = "Gene Expression for Drug X Trial for B-Cell vs T-Cell Patients")

#normalize - brings data to a 0 to 1 scale by subracting the minimum and dividing by the maximum of all observations - preserves each variables distribution while making them comparable
heatmaply(normalize(tmat), main = "Gene Expression for Drug X Trial for B-Cell vs T-Cell Patients")


#percentize - similar to ranking, ecdf of the variables to get an empirical percentile - % of observations above or below it
heatmaply(percentize(tmat), main = "Gene Expression for Drug X Trial for B-Cell vs T-Cell Patients")

#colors
heatmaply(percentize(tmat), main = "Gene Expression for Drug X Trial for B-Cell vs T-Cell Patients", colors = pal)

#add annotation
heatmaply(genes, main = "Gene Expression for Drug X Trial for B-Cell vs T-Cell Patients", colors = pal,
          row_side_colors = annotation_row, col_side_colors = annotation_col, dendrogram="none")

#color for annotation
heatmaply(genes, main = "Gene Expression for Drug X Trial for B-Cell vs T-Cell Patients", 
          row_side_colors = annotation_row, col_side_colors = annotation_col, dendrogram="none", plot_method = "plotly", cellnote_size = 16, column_text_angle = 90, margin = c(80, 70), fontsize_row = 8, fontsize_col = 8) 

# Chargement des librairies
library(lme4)
library(car)
library(cluster)
library(factoextra)
library(ggplot2)
library(ComplexHeatmap)
library(matrixStats)

# Chargement des données
expression <- read.csv("MT/ko_mt_gff_tpm.csv", row.names = 1)
colnames(expression) <- gsub("\\.", "-", colnames(expression))
cat("Nombre de KO dans le fichier :", nrow(expression), "\n")

metadata <- read.csv("MT/sample_covariates_conditions.tsv", sep = "\t", header = TRUE)
metadata$lib_size_raw_zscore <- scale(metadata$lib_size_raw)
metadata$prc_retained <- metadata$lib_size_post / metadata$lib_size_raw
metadata$sex <- as.factor(metadata$sex)
metadata$family <- as.factor(metadata$family)
#metadata$Group <- as.factor(metadata$condition)
metadata$Group <- factor(metadata$condition, levels = c("control", "case"))



# Prévalence
filter_by_prevalence <- function(count_data, Group, ineach_group, inany_group) {
  count_data_zeros <- as.data.frame(t(apply(count_data, 2, function(x){
    x <- as.numeric(x) != 0
    temp <- data.frame(Group = Group, present = x)
    x <- aggregate(present ~ Group, data = temp, FUN = sum)
    res <- x$present
    names(res) <- x$Group
    res
  })))

  count_data_zeros_prc <- as.data.frame(t(apply(count_data_zeros, 1, function(x){
    (x / table(Group)) * 100
  })))

  count_data_zeros_prc_f <- count_data_zeros_prc[apply(count_data_zeros_prc, 1, function(x){
    all(x >= ineach_group) | any(x >= inany_group)
  }), ]

  selected_features <- rownames(count_data_zeros_prc_f)
  return(count_data[, selected_features])
}

expression_filtered <- as.data.frame(t(filter_by_prevalence(
  count_data = t(expression),
  Group = metadata$Group,
  ineach_group = 15,
  inany_group = 50
)))

# Modèle mixte linéaire pour chaque gène
failed_genes <- c()
expression_residuals <- expression_filtered

expression_residuals[,] <- Reduce(rbind, lapply(1:nrow(expression_filtered), function(i){
  gene <- rownames(expression_filtered)[i]
  x <- as.numeric(expression_filtered[i,])
  tryCatch({
    residuals(lmer(x ~ lib_size_raw_zscore + prc_retained + (1|family), data = metadata))
  }, error = function(e){
    failed_genes <<- c(failed_genes, gene)
    rep(NA, length(x))
  })
}))

writeLines(as.character(failed_genes), "MT/ko_failed_in_lmer_mt.txt")

# Analyse VIF pour vérifier la collinéarité
vif_model <- lmer(unlist(expression_filtered[1,]) ~ sex + lib_size_raw_zscore + prc_retained + (1|family), data = metadata, REML = FALSE)
print(car::vif(vif_model))

# Export
write.csv(expression_residuals, "MT/residuals_mt_lmer2.csv", quote = FALSE)

# Analyse des densités et clustering
investigate_data <- function(data, n_points = 100, heatmap = FALSE, nclust = TRUE, show_labels = FALSE) {
  data_scaled <- as.data.frame(apply(data, 2, function(x) scale(x, center = min(x), scale = max(x) - min(x))))
  data_density <- as.data.frame(apply(data_scaled, 2, function(x) {
    x_clean <- na.omit(x)
    if (length(x_clean) >= 2) {
      density(x_clean, n = n_points)$y
    } else {
      rep(NA, n_points)
    }
  }))
  if (all(is.na(data_density))) stop("Tous les profils de densité sont NA")
  data_density_cor <- stats::cor(data_density, use = "pairwise.complete.obs")
  data_density_dist <- as.dist(1 - data_density_cor)
  data_density_t <- as.data.frame(t(data_density))
  cluster_no_plot <- fviz_nbclust(data_density_t, hcut, diss = data_density_dist, method = "silhouette", k.max = 20)
  cluster_no <- which.max(cluster_no_plot$data$y)
  clusterization <- hcut(data_density_dist, k = cluster_no, hc_method = "ward.D2", hc_func = "hclust")$cluster
  return(list(cluster_no_plot = cluster_no_plot, clusters = clusterization))
}

before_resid <- investigate_data(data = t(expression_filtered))
after_resid <- investigate_data(data = t(expression_residuals))

finaldata <- read.csv("MT/residuals_mt_lmer2.csv", row.names = 1)
cat("Nombre de KO dans le fichier :", nrow(finaldata), "\n")


library(matrixStats)
mad_ko <- rowMads(as.matrix(expression_residuals), na.rm = TRUE)
threshold <- quantile(mad_ko, 0.75)  # garder les 25% les plus variables
expression_residuals_filtered <- expression_residuals[mad_ko > threshold, ]


cat("KO conservés après filtrage par présence :", nrow(expression_residuals_filtered), "\n")


write.csv(expression_residuals_filtered, file = "MT/residuals_mt_lmer2_filtered.csv", quote = FALSE)
finaldata <- read.csv("MT/residuals_mt_lmer2_filtered.csv", row.names = 1)
cat("Nombre de KO dans le fichier :", nrow(finaldata), "\n")






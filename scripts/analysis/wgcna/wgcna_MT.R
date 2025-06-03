# === 📦 Chargement des packages nécessaires ===
library(WGCNA)
options(stringsAsFactors = FALSE)

# === 📁 Chargement des données ===
expression_data <- read.csv("MT/residuals_mt_lmer2_filtered.csv", row.names = 1)

# Transposer pour avoir les échantillons en ligne, KO en colonnes
data_expr <- t(expression_data)

# Vérification des données
gsg <- goodSamplesGenes(data_expr, verbose = 3)
if (!gsg$allOK) {
  data_expr <- data_expr[gsg$goodSamples, gsg$goodGenes]
  cat("Removed bad samples or genes\n")
}

# === ⚙️ Paramètres WGCNA ===
powers <- c(1:20)
sft <- pickSoftThreshold(data_expr, powerVector = powers, verbose = 5)

# Choix du soft-threshold (ici supposé 8 si adapté)
softPower <- 13
adjacency <- adjacency(data_expr, power = softPower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

geneTree <- hclust(as.dist(dissTOM), method = "average")

# === 🧱 Modules ===
minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)

dynamicColors <- labels2colors(dynamicMods)
names(dynamicColors) <- colnames(data_expr)
table(dynamicColors)


# === 📊 Eigengenes ===
MEList <- moduleEigengenes(data_expr, colors = dynamicColors)
MEs <- MEList$eigengenes
dissME <- 1 - cor(MEs)
MEtree <- hclust(as.dist(dissME), method = "average")

# === 💾 Sauvegarde des résultats ===
write.csv(dynamicColors, file = "wgcna_module_colors_mg.csv")
write.csv(MEs, file = "wgcna_module_eigengenes_mg.csv")

# === 📈 Graphique (optionnel) ===
pdf("wgcna_gene_tree_mt.pdf")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# ⚠️ Si tu as beaucoup de KO, cette opération peut prendre du temps ou planter
pdf("wgcna_TOM_heatmap_full_mt.pdf", width = 10, height = 10)
TOMplot(dissTOM^7, geneTree, dynamicColors,
        main = "TOM heatmap with modules",
        col = colorRampPalette(c("white", "lightblue", "blue", "darkblue"))(250))

sampleTree <- hclust(dist(data_expr), method = "average")




# Transposer pour avoir les échantillons en ligne, KO en colonnes
data_expr <- t(expression_data)

# Vérification des données
gsg <- goodSamplesGenes(data_expr, verbose = 3)
if (!gsg$allOK) {
  data_expr <- data_expr[gsg$goodSamples, gsg$goodGenes]
  cat("Removed bad samples or genes\n")
}

# === ⚙️ Paramètres WGCNA ===
powers <- c(1:20)
sft <- pickSoftThreshold(data_expr, powerVector = powers, verbose = 5)

# Choix du soft-threshold (ici supposé 8 si adapté)
softPower <- 13
adjacency <- adjacency(data_expr, power = softPower)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

geneTree <- hclust(as.dist(dissTOM), method = "average")

# === 🧱 Modules ===
minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)

dynamicColors <- labels2colors(dynamicMods)
names(dynamicColors) <- colnames(data_expr)
table(dynamicColors)


# === 📊 Eigengenes ===
MEList <- moduleEigengenes(data_expr, colors = dynamicColors)
MEs <- MEList$eigengenes
dissME <- 1 - cor(MEs)
MEtree <- hclust(as.dist(dissME), method = "average")

# === 💾 Sauvegarde des résultats ===
write.csv(dynamicColors, file = "MT/wgcna_module_colors_mt.csv")
write.csv(MEs, file = "MT/wgcna_module_eigengenes_mt.csv")

# === 📈 Graphique (optionnel) ===
pdf("MT/wgcna_gene_tree_mt.pdf")
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Calculer la connectivité (kTotal)
kTotal <- apply(TOM, 1, sum)

# Sélectionner les 500 KO les plus connectés
topKO <- names(sort(kTotal, decreasing = TRUE))[1:500]

# Sous-matrice TOM
TOM_subset <- TOM[topKO, topKO]
dissTOM_subset <- 1 - TOM_subset
geneTree_subset <- hclust(as.dist(dissTOM_subset), method = "average")
dynamicColors_subset <- dynamicColors[topKO]

# Plot TOM heatmap
pdf("MT/wgcna_TOM_heatmap_subset_mt.pdf", width = 9, height = 9)
TOMplot(dissTOM_subset^7, geneTree_subset, dynamicColors_subset,
        main = "TOM heatmap (top 500 KO)")
dev.off()




# ⚠️ Si tu as beaucoup de KO, cette opération peut prendre du temps ou planter
pdf("MT/wgcna_TOM_heatmap_full_mt.pdf", width = 10, height = 10)
TOMplot(dissTOM^7, geneTree, dynamicColors,
        main = "TOM heatmap with modules",
        col = colorRampPalette(c("white", "lightblue", "blue", "darkblue"))(250))

sampleTree <- hclust(dist(data_expr), method = "average")

pdf("MT/sample_clustering_mt.pdf", width = 10, height = 5)
plot(sampleTree, main = "Clustering des échantillons", sub = "", xlab = "", cex.lab = 1.5)
dev.off()



pdf("MT/soft_thresholding_power_mt.pdf", width = 9, height = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit,signed R^2",
     type="n", main = "Scale independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels=powers, cex=0.9, col="red")
abline(h=0.90, col="red")
dev.off()

pdf("MT/module_eigengene_dendrogram_mt.pdf", width = 10, height = 6)
plotEigengeneNetworks(MEs, "", marDendro = c(3,3,2,2), marHeatmap = c(3,4,2,2), plotDendrograms = TRUE)
dev.off()

module <- "blue"  # à adapter
moduleGenes <- which(dynamicColors == module)

MEblue <- MEs$MEblue
#barplot(MEblue, main = "Module Eigengene - Blue", xlab = "Samples", ylab = "Eigengene expression")

# Sélection des covariables
traitData <- metadata[, c("condition", "sex", "family", "lib_size_raw", "lib_size_post")]

# Transformation des variables catégorielles en numériques
traitData$condition <- as.numeric(traitData$condition == "case")
traitData$sex <- as.numeric(traitData$sex == "M")
traitData$family <- as.numeric(as.factor(traitData$family))  # convertit chaque famille en code numérique

# Calcul des corrélations
moduleTraitCor <- cor(MEs, traitData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(data_expr))

# Visualisation
pdf("MT/module_trait_heatmap_mt.pdf", width = 10, height = 10)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 10, 4, 2))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traitData),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = "Module-trait relationships")
dev.off()


# Calcul des corrélations module-KO (kME)
geneModuleMembership <- signedKME(data_expr, MEs)

# S'assurer que les noms des KO sont bien là
rownames(geneModuleMembership) <- colnames(data_expr)


# Définir le nom du module
module_color <- "yellow"

# Identifier les KO dans le module
probes <- names(dynamicColors)[dynamicColors == module_color]

# Nom de la colonne kME dans geneModuleMembership
kme_col <- paste0("kME", module_color)

# Vérifier que la colonne kME existe
if (!(kme_col %in% colnames(geneModuleMembership))) {
  stop(paste("❌ Colonne", kme_col, "absente dans geneModuleMembership."))
}

# Extraire les valeurs de kME pour les KO du module
moduleKME <- geneModuleMembership[match(probes, rownames(geneModuleMembership)), kme_col]

# Donner les noms aux valeurs
names(moduleKME) <- probes

# Vérifier qu'on a bien des valeurs
if (all(is.na(moduleKME))) {
  stop("❌ Toutes les valeurs de kME sont NA.")
}

# Sélectionner les 10 KO avec la plus forte appartenance (kME)
topHubGenes <- sort(moduleKME, decreasing = TRUE)[1:10]

# Créer un data.frame pour l’export
topHubGenes_df <- data.frame(
  KO = names(topHubGenes),
  kME = as.numeric(topHubGenes)
)

# Afficher
print("🔹 Top hub genes du module jaune :")
print(topHubGenes_df)

# Exporter
write.csv(topHubGenes_df,
          file = paste0("MT/top_hub_KO_module_yellow_mt", module_color, ".csv"),
          row.names = FALSE)



probes <- colnames(data_expr)
inModule <- is.finite(match(dynamicColors, module_color))
modProbes <- probes[inModule]
modTOM <- TOM[inModule, inModule]
dimnames(modTOM) <- list(modProbes, modProbes)

exportNetworkToCytoscape(modTOM,
                         edgeFile = paste("MT/CytoscapeInput-edges-mt-", module_color, ".txt", sep=""),
                         nodeFile = paste("MT/CytoscapeInput-nodes-mt-", module_color, ".txt", sep=""),
                         weighted = TRUE, threshold = 0.02,
                         nodeNames = modProbes, nodeAttr = dynamicColors[inModule])


dev.off()


# Définir le nom du module
module_purple <- "turquoise"

# Identifier les KO dans le module
probes <- names(dynamicColors)[dynamicColors == module_purple]

# Nom de la colonne kME dans geneModuleMembership
kme_col <- paste0("kME", module_purple)

# Vérifier que la colonne kME existe
if (!(kme_col %in% colnames(geneModuleMembership))) {
  stop(paste("❌ Colonne", kme_col, "absente dans geneModuleMembership."))
}

# Extraire les valeurs de kME pour les KO du module
moduleKME <- geneModuleMembership[match(probes, rownames(geneModuleMembership)), kme_col]

# Donner les noms aux valeurs
names(moduleKME) <- probes

# Vérifier qu'on a bien des valeurs
if (all(is.na(moduleKME))) {
  stop("❌ Toutes les valeurs de kME sont NA.")
}

# Sélectionner les 10 KO avec la plus forte appartenance (kME)
topHubGenes <- sort(moduleKME, decreasing = TRUE)[1:10]

# Créer un data.frame pour l’export
topHubGenes_df <- data.frame(
  KO = names(topHubGenes),
  kME = as.numeric(topHubGenes)
)

# Afficher
print("🔹 Top hub genes du module purple :")
print(topHubGenes_df)

# Exporter
write.csv(topHubGenes_df,
          file = paste0("MT/top_hub_KO_module_mt_", module_purple, ".csv"),
          row.names = FALSE)

probes <- colnames(data_expr)
inModule <- is.finite(match(dynamicColors, module_purple))
modProbes <- probes[inModule]
modTOM <- TOM[inModule, inModule]
dimnames(modTOM) <- list(modProbes, modProbes)

exportNetworkToCytoscape(modTOM,
                         edgeFile = paste("MT/CytoscapeInput-edges-mt-", module_purple, ".txt", sep=""),
                         nodeFile = paste("MT/CytoscapeInput-nodes-mt-", module_purple, ".txt", sep=""),
                         weighted = TRUE, threshold = 0.02,
                         nodeNames = modProbes, nodeAttr = dynamicColors[inModule])



library(igraph)
library(ggraph)
library(ggplot2)
library(tidygraph)

# === 📁 Chargement des données ===
edges <- read.table("MT/CytoscapeInput-edges-mt-yellow.txt", header = TRUE, sep = "\t")
nodes <- read.table("MT/CytoscapeInput-nodes-mt-yellow.txt", header = TRUE, sep = "\t")
hubs <- read.csv("MT/top_hub_KO_module_mt_yellow.csv", header = TRUE)

# === 📊 Création du graphe ===
g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

# === 🔢 Calcul degré pondéré
degree_weighted <- strength(g, mode = "all", weights = E(g)$weight)
names(degree_weighted) <- V(g)$name

# 🎯 Top 50 KO + tous les hub genes
top_KO <- names(sort(degree_weighted, decreasing = TRUE))[1:50]
all_KO <- unique(c(top_KO, as.character(hubs$KO)))

# === 📉 Sous-graphe
sub_g <- induced_subgraph(g, vids = intersect(V(g)$name, all_KO))

# === 🎨 Mise en forme
V(sub_g)$is_hub <- V(sub_g)$name %in% hubs$KO
V(sub_g)$color <- ifelse(V(sub_g)$is_hub, "firebrick", "gold")
V(sub_g)$label.color <- ifelse(V(sub_g)$is_hub, "black", NA)
V(sub_g)$size <- ifelse(V(sub_g)$is_hub, 6, 2)
V(sub_g)$label <- ifelse(V(sub_g)$is_hub, V(sub_g)$name, NA)

# === 🔄 Conversion pour ggraph
g_tidy <- as_tbl_graph(sub_g)

# === 🖼️ Tracé final
png("wgcna_module_yellow_nice_final.png", width = 1400, height = 1200, res = 150)
ggraph(g_tidy, layout = "fr") +
  geom_edge_link(color = "grey70", alpha = 0.5) +
  geom_node_point(aes(color = is_hub, size = size)) +
  geom_node_text(aes(label = label), repel = TRUE, size = 3.5) +
  scale_color_manual(values = c("FALSE" = "gold", "TRUE" = "firebrick")) +
  scale_size_identity() +
  ggtitle("Module yellow – Réseau WGCNA (hub genes en rouge)") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()


library(igraph)
library(ggraph)
library(ggplot2)
library(tidygraph)

# === 📁 Chargement des données du module turquoise ===
edges <- read.table("MT/CytoscapeInput-edges-mt-turquoise.txt", header = TRUE, sep = "\t")
nodes <- read.table("MT/CytoscapeInput-nodes-mt-turquoise.txt", header = TRUE, sep = "\t")
hubs <- read.csv("MT/top_hub_KO_module_mt_turquoise.csv", header = TRUE)

# === 📊 Création du graphe ===
g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

# === 🔢 Calcul degré pondéré
degree_weighted <- strength(g, mode = "all", weights = E(g)$weight)
names(degree_weighted) <- V(g)$name

# 🎯 Top 50 KO + tous les hub genes
top_KO <- names(sort(degree_weighted, decreasing = TRUE))[1:50]
all_KO <- unique(c(top_KO, as.character(hubs$KO)))

# === 📉 Sous-graphe
sub_g <- induced_subgraph(g, vids = intersect(V(g)$name, all_KO))

# === 🎨 Mise en forme
V(sub_g)$is_hub <- V(sub_g)$name %in% hubs$KO
V(sub_g)$color <- ifelse(V(sub_g)$is_hub, "firebrick", "turquoise")
V(sub_g)$label.color <- ifelse(V(sub_g)$is_hub, "black", NA)
V(sub_g)$size <- ifelse(V(sub_g)$is_hub, 6, 2)
V(sub_g)$label <- ifelse(V(sub_g)$is_hub, V(sub_g)$name, NA)

# === 🔄 Conversion pour ggraph
g_tidy <- as_tbl_graph(sub_g)

# === 🖼️ Tracé final
png("wgcna_module_turquoise_nice_final.png", width = 1400, height = 1200, res = 150)
ggraph(g_tidy, layout = "fr") +
  geom_edge_link(color = "grey70", alpha = 0.5) +
  geom_node_point(aes(color = is_hub, size = size)) +
  geom_node_text(aes(label = label), repel = TRUE, size = 3.5) +
  scale_color_manual(values = c("FALSE" = "turquoise", "TRUE" = "firebrick")) +
  scale_size_identity() +
  ggtitle("Module turquoise – Réseau WGCNA (hub genes en rouge)") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
dev.off()



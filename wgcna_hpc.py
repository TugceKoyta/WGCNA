import os
import pandas as pd
import scanpy as sc
import omicverse as ov
import matplotlib.pyplot as plt
from statsmodels import robust

# Define input and output directories
INPUT_DIR = "/home/users/tkoytaviloglu/data/raw/"
OUTPUT_DIR = "/home/users/tkoytaviloglu/results/"
FIGURE_DIR = os.path.join(OUTPUT_DIR, "figures/")
TABLE_DIR = os.path.join(OUTPUT_DIR, "outputs/")

# Ensure output directories exist
os.makedirs(FIGURE_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)

# Log current working directory
print(f"Working directory: {os.getcwd()}")

# Step 1: Load expression data
expression_file = os.path.join(INPUT_DIR, "expressionList.csv")
print(f"Loading expression data from: {expression_file}")
data = pd.read_csv(expression_file, index_col=0)

# Step 2: Filter top 2000 most variable genes using MAD
gene_mad = data.apply(robust.mad)
data = data.T  # Transpose for PyWGCNA compatibility
data = data.loc[gene_mad.sort_values(ascending=False).index[:2000]]

# Step 3: Initialize PyWGCNA
pyWGCNA_5xFAD = ov.bulk.pyWGCNA(
            name='5xFAD_2k',
                species='mus musculus',
                    geneExp=data.T,
                        outputPath=OUTPUT_DIR,
                            save=True
                            )
print(" PyWGCNA initialized successfully.")

# Step 4: Preprocessing
print(" Running preprocessing...")
pyWGCNA_5xFAD.preprocess()
print(" Preprocessing completed.")

# Step 5: Network construction
print(" Calculating soft threshold...")
pyWGCNA_5xFAD.calculate_soft_threshold()
print(" Soft threshold calculated.")

print(" Calculating adjacency matrix...")
pyWGCNA_5xFAD.calculating_adjacency_matrix()
print(" Adjacency matrix computed.")

print(" Calculating TOM similarity matrix...")
pyWGCNA_5xFAD.calculating_TOM_similarity_matrix()
print(" TOM similarity matrix computed.")

# Step 6: Module detection
print(" Calculating gene tree...")
pyWGCNA_5xFAD.calculate_geneTree()
print(" Gene tree computed.")

print(" Detecting dynamic modules...")
pyWGCNA_5xFAD.calculate_dynamicMods(kwargs_function={'cutreeHybrid': {'deepSplit': 2, 'pamRespectsDendro': False}})
print(" Modules detected.")

print(" Calculating module eigengenes...")
pyWGCNA_5xFAD.calculate_gene_module(kwargs_function={'moduleEigengenes': {'softPower': 8}})
print(" Module eigengenes computed.")

#  Step 7: Save WGCNA results
wgcna_output_file = os.path.join(OUTPUT_DIR, "5xFAD_2k.p")
print(f" Saving WGCNA results to: {wgcna_output_file}")
pyWGCNA_5xFAD.saveWGCNA()
print(" WGCNA results saved.")

#  Step 8: Reload saved results
print(f" Reloading WGCNA results from: {wgcna_output_file}")
pyWGCNA_5xFAD = ov.bulk.readWGCNA(wgcna_output_file)
print(" WGCNA results reloaded.")

# Step 9: Load sample information
sample_file = os.path.join(INPUT_DIR, "sampleInfo.csv")
print(f" Loading sample information from: {sample_file}")
pyWGCNA_5xFAD.updateSampleInfo(path=sample_file, sep=',')
print(" Sample information updated.")

# Step 10: Assign metadata colors
print(" Assigning colors to metadata...")
pyWGCNA_5xFAD.setMetadataColor('Sex', {'Female': 'green', 'Male': 'yellow'})
pyWGCNA_5xFAD.setMetadataColor('Genotype', {'5xFADWT': 'darkviolet', '5xFADHEMI': 'deeppink'})
pyWGCNA_5xFAD.setMetadataColor('Age', {'4mon': 'thistle', '8mon': 'plum', '12mon': 'violet', '18mon': 'purple'})
pyWGCNA_5xFAD.setMetadataColor('Tissue', {'Hippocampus': 'red', 'Cortex': 'blue'})
print(" Metadata colors assigned.")

# Step 11: Perform WGCNA analysis
print(" Running WGCNA analysis...")
pyWGCNA_5xFAD.analyseWGCNA()
print(" WGCNA analysis completed.")

# Step 12: Identify hub genes & save them
print(" Identifying hub genes for 'lightgreen' module...")
lightgreen_hubs = pyWGCNA_5xFAD.top_n_hub_genes(moduleName="lightgreen", n=10)
lightgreen_hubs.to_csv(os.path.join(TABLE_DIR, "hub_genes_lightgreen.csv"), index=True)

print(" Identifying hub genes for 'gold' module...")
gold_hubs = pyWGCNA_5xFAD.top_n_hub_genes(moduleName="gold", n=10)
gold_hubs.to_csv(os.path.join(TABLE_DIR, "hub_genes_gold.csv"), index=True)

# Step 13: Save and plot results
print(" Generating and saving plots...")

# Module-trait heatmap
fig, ax = plt.subplots(figsize=(10, 8))
pyWGCNA_5xFAD.analyseWGCNA()
plt.savefig(os.path.join(FIGURE_DIR, "module_trait_relationship.png"), dpi=300)
plt.close()

# Eigengene barplot
metadata = pyWGCNA_5xFAD.datExpr.obs.columns.tolist()  # Retrieve metadata
fig, ax = plt.subplots(figsize=(8, 6))
pyWGCNA_5xFAD.barplotModuleEigenGene('lightgreen', metadata, show=False)
plt.savefig(os.path.join(FIGURE_DIR, "eigengene_barplot_lightgreen.png"), dpi=300)
plt.close()

#  Subnetwork visualization
fig, ax = plt.subplots(figsize=(8, 8))
pyWGCNA_5xFAD.plot_sub_network(['gold', 'lightgreen'], pos_type='kamada_kawai',
                                pos_scale=10, pos_dim=2, figsize=(8, 8),
                                node_size=10, label_fontsize=8,
                                correlation_threshold=0.2,
                                label_bbox={"ec": "white", "fc": "white", "alpha": 0.6})
plt.savefig(os.path.join(FIGURE_DIR, "subnetwork_gold_lightgreen.png"), dpi=300)
plt.close()

print(f"✅ All results saved in: {OUTPUT_DIR}")
print("🚀 WGCNA analysis completed successfully.")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# 📂 Chargement des données
ko_data_path = "/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/ko_mt_mg_ratio.parquet"
ko_data = pd.read_parquet(ko_data_path)

# 🔍 Nettoyage noms colonnes (au cas où)
ko_data.columns = ko_data.columns.astype(str).str.strip()
ko_data = ko_data.T  # transpose : les KOs deviennent lignes, les échantillons deviennent colonnes


# 📄 Chargement du fichier de conditions
conditions_file = "/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/conditions.txt"
conditions_df = pd.read_csv(conditions_file, sep=r'\s+', header=None, names=["sample", "condition"])
conditions_df["sample"] = conditions_df["sample"].astype(str).str.strip()

# Ne garde que les échantillons présents dans tes colonnes KO
conditions_df = conditions_df[conditions_df["sample"].isin(ko_data.columns)]

# Reset l'index pour garantir qu'il n'y a pas de conflits
conditions_df = conditions_df.reset_index(drop=True)

# Set index to match count matrix
conditions_df = conditions_df.set_index("sample")

# Convertir en category avec ordre explicite
conditions_df["condition"] = pd.Categorical(
            conditions_df["condition"],
                categories=["control", "case"],  # "control" sera pris comme ref
                    ordered=True
                    )

print(conditions_df.head())
print(conditions_df["condition"].dtype)
print(conditions_df["condition"].cat.categories)


# 🔁 Filtrage des échantillons valides
valid_samples = [s for s in conditions_df.index if s in ko_data.columns]

# 🔬 Filtrage des données de KO
ko_counts = ko_data[valid_samples]

# 🧱 Construction de la design_matrix
design_matrix = conditions_df.loc[valid_samples]

print(f"✅ {len(valid_samples)} échantillons utilisés")
print("📐 Shape ko_counts :", ko_counts.shape)
print("📐 Shape design_matrix :", design_matrix.shape)

# 💡 Assure-toi que les index matchent exactement
design_matrix = design_matrix.loc[ko_counts.columns]

# ⚠ Vérification alignement
assert list(ko_counts.columns) == list(design_matrix.index), "Les colonnes de counts ne correspondent pas à l’index de metadata"

print("🔬 Premiers noms dans ko_data.columns :", list(ko_data.columns[:5]))
print("🔬 Premiers noms dans conditions_df.index :", list(conditions_df.index[:5]))

print("Index de conditions_df : ", conditions_df.index)  # Affiche l'index
print("Colonnes de conditions_df : ", conditions_df.columns)  # Affiche les colonnes


# Vérifie les différences exactes
# Vérifie les différences exactes
cols_not_in_conditions = [col for col in ko_data.columns if col not in conditions_df.index.tolist()]  # Utilisation de conditions_df.index
samples_not_in_ko = [s for s in conditions_df.index if s not in ko_data.columns]

print(f"❌ Colonnes KO non présentes dans conditions.txt ({len(cols_not_in_conditions)}) :", cols_not_in_conditions[:5])
print(f"❌ Samples dans conditions.txt non présents dans ko_data ({len(samples_not_in_ko)}) :", samples_not_in_ko[:5])


# Force l'ordre de la condition, "control" étant la catégorie de référence
design_matrix['condition'] = design_matrix['condition'].cat.reorder_categories(['control', 'case'], ordered=True)

print(design_matrix['condition'].unique())
print(design_matrix['condition'].cat.categories)

# ✅ Vérifie bien que "condition" contient exactement 2 valeurs uniques
print("Valeurs uniques dans 'condition' :", design_matrix["condition"].unique())
print("Type de la colonne :", design_matrix["condition"].dtype)

# ✅ Force les catégories et l’ordre (très important)
design_matrix["condition"] = pd.Categorical(
            design_matrix["condition"],
                categories=["case", "control"]
                )

# ✅ Recheck pour debug visuel
print("Catégories forcées :", design_matrix["condition"].cat.categories)
print("Premières lignes de la design_matrix :")
print(design_matrix.head())

# Convertir les données de comptage en entier
dds_counts = ko_counts.T.fillna(0).astype(int)

# 🧪 DESeq2
dds = DeseqDataSet(
        counts=dds_counts,
        metadata=design_matrix,  # Utilisation de la design_matrix ici
        design_factors=["condition"],
        )

# Lancer DESeq2
dds.deseq2()
stats = DeseqStats(dds)
stats.summary()
res_df = stats.results_df

# 🎯 Filtrage
res_df_filtered = res_df[
        (res_df['padj'] < 0.05) &
        (res_df['log2FoldChange'] > 1) &
        (res_df['baseMean'] > 10)
].sort_values(by='padj')

#res_df_filtered = res_df[(res_df["padj"] < 0.05) & (res_df["log2FoldChange"] > 1)]

top20 = res_df_filtered.sort_values(by='padj').head(20)
top20.to_csv("top20_ko_case_enriched.csv")
print(top20)

# 💾 Sauvegarde des résultats
output_path = "/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/"
res_df.to_csv(output_path + "ko_comparison_results_deseq2.csv")
res_df_filtered.to_csv(output_path + "ko_significatifs_deseq2.csv")

print(f"✅ {res_df_filtered.shape[0]} KO significatifs enrichis chez les malades")

# 📊 Volcano Plot

res_df['-log10(padj)'] = -np.log10(res_df['padj'])

plt.figure(figsize=(12, 6))
plt.scatter(res_df['log2FoldChange'], res_df['-log10(padj)'], alpha=0.3, color='grey', label='Non signifiants')

plt.scatter(res_df_filtered['log2FoldChange'], -np.log10(res_df_filtered['padj']), color='red', alpha=0.7, label='Significatifs')

plt.axhline(-np.log10(0.05), color='blue', linestyle='--', label='padj=0.05')
plt.axvline(1, color='green', linestyle='--', label='log2FC = ±1')
plt.axvline(-1, color='green', linestyle='--')

plt.xlabel('log2 Fold Change')
plt.ylabel('-log10(padj)')
plt.title('Volcano plot: KO enrichis chez les malades')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(output_path + "volcano_plot_deseq2.png", dpi=300)
plt.show()


#plt.figure(figsize=(10, 6))
#plt.scatter(res_df["log2FoldChange"], -np.log10(res_df["padj"]), color="gray", alpha=0.5)
#plt.scatter(res_df_filtered["log2FoldChange"], -np.log10(res_df_filtered["padj"]), color="red", label="Significant KOs")
#plt.axvline(x=1, linestyle="--", color="black", label="log2FC > 1")
#plt.axhline(y=-np.log10(0.05), linestyle="--", color="blue", label="padj < 0.05")
#plt.xlabel("log2 Fold Change")
#plt.ylabel("-log10 Adjusted P-value")
#plt.title("Volcano Plot of KO significance (PyDESeq2)")
#plt.legend()
#plt.tight_layout()
#plt.savefig(output_path + "volcano_plot_deseq2.png")
#plt.show()


#significant_genes = res_df[res_df['padj'] < 0.05]
#print(significant_genes.head())

#plt.figure(figsize=(10, 8))
#plt.scatter(res_df['log2FoldChange'], -np.log10(res_df['pvalue']), c='red', alpha=0.5)
#plt.axhline(y=-np.log10(0.05), color='blue', linestyle='--')
#plt.title('Volcano plot')
#plt.xlabel('Log2 Fold Change')
#plt.ylabel('-Log10 p-value')
#plt.show()


plt.figure(figsize=(10, 6))
plt.scatter(res_df['baseMean'], res_df['log2FoldChange'], alpha=0.3, color='grey', label='Non signifiants')

plt.scatter(res_df_filtered['baseMean'], res_df_filtered['log2FoldChange'], alpha=0.6, color='red', label='Significatifs')

plt.xscale('log')
plt.axhline(0, color='black', linestyle='--')
plt.xlabel('baseMean (log scale)')
plt.ylabel('log2 Fold Change')
plt.title('MA plot')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("ma_plot.png", dpi=300)
plt.show()





from bioservices import KEGG

kegg = KEGG()

ko_ids = list(top20.index)  # ou res_df_filtered.index si tu veux tout annoter

annotations = {}

for ko in ko_ids:
    try:
        entry = kegg.get(ko)
        parsed = kegg.parse(entry)
        annotations[ko] = parsed.get('NAME', 'NA')
    except Exception as e:
        annotations[ko] = 'NA'
        print(f"Erreur pour {ko} : {e}")

# Ajouter au DataFrame
top20['KO_annotation'] = top20.index.map(annotations)
top20.to_csv("top20_ko_enriched_annotated.csv")
print(top20[['log2FoldChange', 'padj', 'KO_annotation']])

# Ajout à ton dataframe
#res_df_filtered['KO_annotation'] = res_df_filtered.index.map(annotations)
#res_df_filtered.to_csv(output_path + "ko_enriched_annotated.csv")

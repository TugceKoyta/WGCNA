import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Charger les données KO
ko_data_path = "/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/ko_mt_mg_ratio.parquet"
ko_data = pd.read_parquet(ko_data_path)

# Fonction pour lire les fichiers contenant les échantillons
def lire_echantillons(fichier):
    try:
        with open(fichier, "r") as f:
            return [ligne.strip() for ligne in f.readlines()]
    except FileNotFoundError:
        print(f"⚠ Erreur : Fichier {fichier} introuvable.")
        return []

# Charger les listes d'échantillons malades et sains
samples_malades_path = "/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/samples_malades.txt"
samples_sains_path = "/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/samples_sains.txt"

samples_malades = lire_echantillons(samples_malades_path)
samples_sains = lire_echantillons(samples_sains_path)

# Vérifier les échantillons réellement présents dans ko_data
samples_malades_valid = [sample for sample in samples_malades if sample in ko_data.index]
samples_sains_valid = [sample for sample in samples_sains if sample in ko_data.index]

# Avertir si certains échantillons sont manquants
if len(samples_malades_valid) != len(samples_malades):
    print(f"⚠ Avertissement : Certains échantillons malades sont absents des données.")
    if len(samples_sains_valid) != len(samples_sains):
        print(f"⚠ Avertissement : Certains échantillons sains sont absents des données.")


# Filtrer les KO pour les échantillons valides
ko_malades = ko_data.loc[samples_malades_valid]
ko_sains = ko_data.loc[samples_sains_valid]

print(f"✅   {ko_malades.shape[1]} KO chargés pour {len(samples_malades_valid)} malades.")
print(f"✅   {ko_sains.shape[1]} KO chargés pour {len(samples_sains_valid)} sains.")

# Effectuer un test statistique pour comparer les groupes
results = []
for ko in ko_malades.columns:
    malades_vals = ko_malades[ko].dropna()
    sains_vals = ko_sains[ko].dropna()
    if len(malades_vals) > 1 and len(sains_vals) > 1:
        t_stat, p_val = stats.ttest_ind(malades_vals, sains_vals, nan_policy="omit", equal_var=False)
        results.append((ko, t_stat, p_val))

# Convertir les résultats en DataFrame
results_df = pd.DataFrame(results, columns=["KO", "T-statistic", "P-value"])

# Filtrer les KO significatifs
significant_ko = results_df[results_df["P-value"] < 0.05]

# Ajouter log2FC
mean_malades = ko_malades.mean()
mean_sains = ko_sains.mean()
log2fc = np.log2((mean_malades + 1e-6) / (mean_sains + 1e-6))
significant_ko = significant_ko.copy()
significant_ko["log2FC"] = significant_ko["KO"].map(log2fc)
significant_ko = significant_ko[significant_ko["log2FC"] > 1]
# Ajouter log2FC à tous les KO
results_df["log2FC"] = results_df["KO"].map(log2fc)

print(f"✅   {len(significant_ko)} KO significatifs avec enrichissement chez les malades.")

# Sauvegarder les résultats
output_path = "/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/"
significant_ko.to_csv(output_path + "ko_significatifs.csv", index=False)
results_df.to_csv(output_path + "ko_comparison_results.csv", index=False)

# Visualisation : Boxplot

# Construire un DataFrame pour le boxplot
boxplot_data = pd.DataFrame({
        "Abondance KO": np.concatenate([ko_malades.mean(axis=1), ko_sains.mean(axis=1)]),
            "Groupe": ["Malade"] * len(ko_malades) + ["Sain"] * len(ko_sains)
        })

# Tracer le boxplot avec Seaborn
plt.figure(figsize=(8, 6))
sns.boxplot(x="Groupe", y="Abondance KO", data=boxplot_data)
plt.title("Comparaison de l'abondance des KO entre malades et sains")
plt.savefig(output_path + "boxplot_ko.png")
plt.show()

# Visualisation : Volcano Plot
plt.figure(figsize=(10, 6))
plt.scatter(results_df["log2FC"], -np.log10(results_df["P-value"]), color="gray", alpha=0.5)
plt.scatter(significant_ko["log2FC"], -np.log10(significant_ko["P-value"]), color="red", label="KO significatifs")
plt.axhline(y=-np.log10(0.05), color="blue", linestyle="--", label="Seuil de 0.05")
plt.xlabel("Log2 Fold Change")
plt.ylabel("-log10 P-value")
plt.title("Volcano Plot des KO")
plt.legend()
plt.savefig(output_path + "volcano_plot.png")
plt.show()

# Visualisation : Barplot des 10 KO les plus significatifs
plt.figure(figsize=(12, 6))
top_ko = significant_ko.nsmallest(10, "P-value")
sns.barplot(y=top_ko["KO"], x=-np.log10(top_ko["P-value"]), hue=top_ko["KO"], palette="viridis", legend=False)
plt.xlabel("-log10 P-value")
plt.ylabel("KO")
plt.title("Top 10 KO les plus significatifs")
plt.savefig(output_path + "top10_ko.png")
plt.show()

# Visualisation : Heatmap des KO significatifs
if len(significant_ko) > 1:
    top_ko_heatmap = significant_ko["KO"].tolist()
    data_heatmap = ko_data[top_ko_heatmap].T
    plt.figure(figsize=(12, 8))
    sns.heatmap(data_heatmap, cmap="coolwarm", xticklabels=False)
    plt.xlabel("Échantillons")
    plt.ylabel("KO")
    plt.title("Heatmap des KO les plus significatifs")
    plt.savefig(output_path + "heatmap_ko.png")
    plt.show()

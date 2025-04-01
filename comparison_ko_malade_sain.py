import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

### 🔹 Charger les données KO
ko_data_path = "/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/ko_mt_mg_ratio.parquet"
ko_data = pd.read_parquet(ko_data_path)

### 🔹 Fonction pour lire les fichiers contenant les échantillons
def lire_echantillons(fichier):
    try:
        with open(fichier, "r") as f:
            return [ligne.strip() for ligne in f.readlines()]
    except FileNotFoundError:
        print(f"⚠️ Erreur : Fichier {fichier} introuvable.")
        return []

### 🔹 Charger les listes d'échantillons malades et sains
samples_malades_path = "/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/samples_malades.txt"
samples_sains_path = "/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/samples_sains.txt"

samples_malades = lire_echantillons(samples_malades_path)
samples_sains = lire_echantillons(samples_sains_path)

# Vérifier les échantillons réellement présents dans ko_data
samples_malades_valid = [sample for sample in samples_malades if sample in ko_data.index]
samples_sains_valid = [sample for sample in samples_sains if sample in ko_data.index]

# Avertir si certains échantillons sont manquants
if len(samples_malades_valid) != len(samples_malades):
    print(f"⚠️ Avertissement : Certains échantillons malades sont absents des données.")
if len(samples_sains_valid) != len(samples_sains):
    print(f"⚠️ Avertissement : Certains échantillons sains sont absents des données.")

### 🔹 Filtrer les KO pour les échantillons valides
ko_malades = ko_data.loc[samples_malades_valid]
ko_sains = ko_data.loc[samples_sains_valid]

print(f"✅ {ko_malades.shape[1]} KO chargés pour {len(samples_malades_valid)} malades.")
print(f"✅ {ko_sains.shape[1]} KO chargés pour {len(samples_sains_valid)} sains.")

### 🔹 Effectuer un test statistique pour comparer les groupes
results = []
for ko in ko_malades.columns:
    malades_vals = ko_malades[ko].dropna()
    sains_vals = ko_sains[ko].dropna()

    if len(malades_vals) > 1 and len(sains_vals) > 1:  # Éviter les erreurs avec un seul point
        t_stat, p_val = stats.ttest_ind(malades_vals, sains_vals, nan_policy="omit", equal_var=False)
        results.append((ko, t_stat, p_val))

### 🔹 Convertir les résultats en DataFrame
results_df = pd.DataFrame(results, columns=["KO", "T-statistic", "P-value"])

# Filtrer les KO avec une p-value < 0.05
significant_ko = results_df[results_df["P-value"] < 0.05]

### 🔹 Ajouter un critère supplémentaire : enrichissement chez les malades
mean_malades = ko_malades.mean()
mean_sains = ko_sains.mean()

# Calculer le log2 fold-change entre malades et sains
log2fc = np.log2((mean_malades + 1e-6) / (mean_sains + 1e-6))  # Éviter la division par 0
significant_ko["log2FC"] = significant_ko["KO"].map(log2fc)

# Filtrer les KO avec un enrichissement minimum (ex: log2FC > 1 signifie un enrichissement x2)
significant_ko = significant_ko[significant_ko["log2FC"] > 1]

print(f"✅ {len(significant_ko)} KO significatifs avec enrichissement chez les malades.")

### 🔹 Sauvegarder les résultats
output_path = "/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/"
significant_ko.to_csv(output_path + "ko_significatifs.csv", index=False)
results_df.to_csv(output_path + "ko_comparison_results.csv", index=False)

### 🔹 Visualisation des P-values
plt.figure(figsize=(10, 6))
plt.scatter(results_df["KO"], results_df["P-value"], color="blue", alpha=0.5, label="P-values")
plt.axhline(y=0.05, color="red", linestyle="--", label="Seuil de 0.05")
plt.xlabel("KO")
plt.ylabel("P-value")
plt.title("P-value des tests t pour chaque KO")
plt.legend()
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(output_path + "significant_ko.png")
plt.show()

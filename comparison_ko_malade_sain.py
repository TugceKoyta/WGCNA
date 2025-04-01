import pandas as pd
from scipy import stats

# Charger le fichier .parquet
ko_data = pd.read_parquet('/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/ko_mt_mg_ratio.parquet')

# Lire les fichiers contenant les échantillons malades et sains
def lire_echantillons(fichier):
    with open(fichier, 'r') as f:
        return [ligne.strip() for ligne in f.readlines()]

# Charger les listes depuis les fichiers
samples_malades = lire_echantillons("/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/samples_malades.txt")
samples_sains = lire_echantillons("/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/samples_sains.txt")


# Vérifier les échantillons présents dans ko_data
print(f"Échantillons présents dans les données: {ko_data.index.tolist()}")

# Filtrer les échantillons malades et sains en fonction de ce qui est dans les données
samples_malades_valid = [sample for sample in samples_malades if sample in ko_data.index]
samples_sains_valid = [sample for sample in samples_sains if sample in ko_data.index]

# Afficher un message si des échantillons sont absents
if len(samples_malades_valid) != len(samples_malades):
    print(f"Avertissement: Certains échantillons malades sont absents dans les données.")
    if len(samples_sains_valid) != len(samples_sains):
        print(f"Avertissement: Certains échantillons sains sont absents dans les données.")

# Sélectionner les KO pour les échantillons valides
ko_malades = ko_data.loc[samples_malades_valid]
ko_sains = ko_data.loc[samples_sains_valid]

# Affichage des résultats pour vérifier
print("KO des échantillons malades :")
print(ko_malades)
print("\nKO des échantillons sains :")
print(ko_sains)

# Calcul du t-test pour chaque KO (chaque colonne)
# On itère sur les KO en tant que colonnes et effectue le test entre les groupes
results = []
for ko in ko_malades.columns:
    try:
        t_stat, p_val = stats.ttest_ind(ko_malades[ko], ko_sains[ko], nan_policy='omit')
        results.append((ko, t_stat, p_val))
    except Exception as e:
        print(f"Erreur lors du calcul pour {ko}: {e}")

# Convertir les résultats en DataFrame pour une meilleure visualisation
results_df = pd.DataFrame(results, columns=['KO', 'T-statistic', 'P-value'])

# Enregistrer les résultats dans un fichier CSV
results_df.to_csv('/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/ko_comparison_results.csv', index=False)

# Affichage des résultats significatifs (p-value < 0.05)
significant_results = results_df[results_df['P-value'] < 0.05]
print("KO significatifs:")
print(significant_results)

# Filtrer les KO avec une p-value significativement inférieure à 0.05
significant_ko = results_df[results_df['P-value'] < 0.05]

# Afficher les KO significatifs
print(significant_ko)

import matplotlib.pyplot as plt

# Tracer les p-values
plt.figure(figsize=(10, 6))
plt.scatter(results_df['KO'], results_df['P-value'], color='blue', label='P-values', alpha=0.5)
plt.axhline(y=0.05, color='red', linestyle='--', label='Seuil de 0.05')
plt.xlabel('KO')
plt.ylabel('P-value')
plt.title('P-value des tests t pour chaque KO')
plt.legend()
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/significant_ko.png')
plt.show()

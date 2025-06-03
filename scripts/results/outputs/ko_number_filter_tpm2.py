import pandas as pd

df_tpm = pd.read_parquet('/scratch/users/tkoytaviloglu/results/outputs/ko_tpm.parquet')

# Définir les seuils
abundance_threshold = 30
presence_threshold = 3
presence_threshold_80 = int(0.8 * len(df_tpm))  # 80% des échantillons

# 1️⃣ Filtrer les KO sur l'abondance (TPM ≥ 30 dans ≥ 3 échantillons)
df_tpm_filtered = df_tpm.loc[:, (df_tpm >= abundance_threshold).sum(axis=0) >= presence_threshold]

# Sauvegarde après filtrage abondance
df_tpm_filtered.to_parquet('/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/ko_filtered_tpm.parquet')
print(f"✅ KO après filtrage par abondance : {df_tpm_filtered.shape[1]}")

# 2️⃣ Filtrer les KO sur la présence dans 80% des échantillons (après filtrage abondance)
df_tpm_80 = df_tpm_filtered.loc[:, (df_tpm_filtered > 0).sum(axis=0) >= presence_threshold_80]

# Sauvegarde après filtrage 80%
df_tpm_80.to_parquet('/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/ko_filtered_80percent_tpm.parquet')
print(f"✅ KO après filtrage 80% des échantillons : {df_tpm_80.shape[1]}")

# 3️⃣ Charger les fichiers TPM pour MG et MT
df_mt = pd.read_parquet('/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/ko_mt_tpm.parquet')
df_mg = pd.read_parquet('/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/ko_mg_tpm.parquet')

# Renommer les colonnes pour retirer les préfixes "mt_" et "mg_" pour les harmoniser
df_mt.columns = df_mt.columns.str.replace('mt_mantis_kegg_ko_counts_', '', regex=False)
df_mg.columns = df_mg.columns.str.replace('mg_mantis_kegg_ko_counts_', '', regex=False)
df_tpm_80.columns = df_tpm_80.columns.str.replace('mt_mantis_kegg_ko_counts_', '', regex=False)

# Vérifier l'intersection des KO entre MT et MG
common_kos = df_tpm_80.columns.intersection(df_mg.columns)

# Filtrer df_mg uniquement sur ces KO communs
df_mg_filtered = df_mg[common_kos]

# Appliquer les mêmes filtrages que précédemment sur df_mt et df_mg
df_mt_filtered = df_mt[df_tpm_80.columns]
df_mg_filtered = df_mg[df_tpm_80.columns]

# 4️⃣ Calcul du ratio MT/MG en évitant la division par zéro
df_ratio = df_mt_filtered / (df_mg_filtered + 1e-6)

# Sauvegarde du ratio MT/MG
df_ratio.to_parquet('/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/ko_mt_mg_ratio.parquet')
print("✅ Ratio MT/MG enregistré")


(base) 0 [tkoytaviloglu@access1 outputs]$  cd ../
(base) 0 [tkoytaviloglu@access1 results]$ ls
figures  outputs
(base) 0 [tkoytaviloglu@access1 results]$ cd outputs/
(base) 0 [tkoytaviloglu@access1 outputs]$ ls
:                                           deseq2_ko.py                               read_mt_mg_ration.py
apache_output.csv                           ko_number_filter2.py                       read_parquet_filtered.py
calcul_threshold_ko_reads.py                ko_number_filter_for_deseq2.py             read_parquet_gff.py
comparison_ko_case_control_deseq2_final.py  ko_number_filter.py                        read_parquet_input_deseq2.py
comparison_ko_case_control_deseq2.py        ko_number_filter_reads.py                  read_parquet.py
comparison_ko_case_control_mg_mt.py         ko_number_filter_tpm2.py                   read_parquet_to_csv.py
comparison_ko_case_control.R                ko_number_filter_tpm_abundance.py          residuals_matrix
comparison_ko_case_control_test2.py         ko_number_filter_tpm.py                    same_ko_number.py
comparison_ko_malade_sain_deseq2.py         ko_number_gff_tpm.py                       taxonomy
comparison_ko_malade_sain.py                ko_number_heatmap_mt.py                    test_du_code_volcanoplot.py
comparison_ko_separately_mg_mt2.py          ko_number_heatmap.py                       top20_ko_case_enriched.csv
comparison_ko_separately_mg_mt.py           ko_number_tpm_comparison_normalisarion.py  top20_ko_enriched_annotated.csv
deseq2_ko                                   ma_plot.png                                upsetplot_KO.png
(base) 0 [tkoytaviloglu@access1 outputs]$ more ko_number_tpm_comparison_normalisarion.py
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df_tpm = pd.read_parquet("/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/ko_filtered_tpm.parquet")
df_tpm_80 = pd.read_parquet("/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/ko_filtered_80percent_tpm.parq
uet")
df_ratio = pd.read_parquet("/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/ko_mt_mg_ratio.parquet")

df_tpm_80.columns = df_tpm_80.columns.str.replace("mt_mantis_kegg_ko_counts_", "", regex=True)


# Histogramme des TPM avant normalisation
plt.figure(figsize=(8, 5))
sns.histplot(df_tpm.values.flatten(), bins=50, kde=True)
plt.xlabel("Valeurs TPM")
plt.ylabel("Fréquence")
plt.title("Distribution des TPM (Sans Normalisation)")
plt.savefig("/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/hist_tpm_sans_normalisation.png")
plt.close()

# Histogramme des TPM après filtrage 80%
plt.figure(figsize=(8, 5))
sns.histplot(df_tpm_80.values.flatten(), bins=50, kde=True, color="orange")
plt.xlabel("Valeurs TPM")
plt.ylabel("Fréquence")
plt.title("Distribution des TPM (Filtrage ≥80% échantillons)")
plt.savefig("/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/hist_tpm_filtrage_80percent.png")
plt.close()

# Boxplot des TPM pour voir la dispersion
plt.figure(figsize=(10, 5))
sns.boxplot(data=df_tpm, showfliers=False)
plt.xticks(rotation=90)
plt.title("Dispersion des TPM par KO (Sans Normalisation)")
plt.savefig("/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/boxplot_tpm_sans_normalisation.png")
plt.close()

# Heatmap du ratio MT/MG (normalisation)
plt.figure(figsize=(12, 6))
sns.heatmap(df_ratio, cmap="viridis", vmax=5)  # Ajuste vmax si nécessaire
plt.title("Ratio MT/MG pour tous les KO")
plt.savefig("/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/heatmap_ratio_mt_mg.png")
plt.close()

print(" Visualisations enregistrées en PNG !")


### Vérifications ###
print("KO TPM (without normalization):", df_tpm.shape)
print("KO TPM (80% filtering):", df_tpm_80.shape)
print("KO Ratio MT/MG:", df_ratio.shape)
print("KO common between df_tpm_80 and df_ratio:", len(set(df_tpm_80.columns).intersection(df_ratio.columns)))

# Vérifier les valeurs min/max du ratio
print("Min / Max Ratio MT/MG:", df_ratio.min().min(), df_ratio.max().max())


print("Exemple KO dans df_tpm_80:", df_tpm_80.columns[:5])
print("Exemple KO dans df_ratio:", df_ratio.columns[:5])
print("KO communs entre df_tpm_80 et df_ratio:", set(df_tpm_80.columns).intersection(df_ratio.columns))

(base) 0 [tkoytaviloglu@access1 outputs]$ ls
:                                           deseq2_ko.py                               read_mt_mg_ration.py
apache_output.csv                           ko_number_filter2.py                       read_parquet_filtered.py
calcul_threshold_ko_reads.py                ko_number_filter_for_deseq2.py             read_parquet_gff.py
comparison_ko_case_control_deseq2_final.py  ko_number_filter.py                        read_parquet_input_deseq2.py
comparison_ko_case_control_deseq2.py        ko_number_filter_reads.py                  read_parquet.py
comparison_ko_case_control_mg_mt.py         ko_number_filter_tpm2.py                   read_parquet_to_csv.py
comparison_ko_case_control.R                ko_number_filter_tpm_abundance.py          residuals_matrix
comparison_ko_case_control_test2.py         ko_number_filter_tpm.py                    same_ko_number.py
comparison_ko_malade_sain_deseq2.py         ko_number_gff_tpm.py                       taxonomy
comparison_ko_malade_sain.py                ko_number_heatmap_mt.py                    test_du_code_volcanoplot.py
comparison_ko_separately_mg_mt2.py          ko_number_heatmap.py                       top20_ko_case_enriched.csv
comparison_ko_separately_mg_mt.py           ko_number_tpm_comparison_normalisarion.py  top20_ko_enriched_annotated.csv
deseq2_ko                                   ma_plot.png                                upsetplot_KO.png
(base) 0 [tkoytaviloglu@access1 outputs]$ more ko_number_filter_tpm2.py
import pandas as pd

# 📌 Définition d'une fonction pour éviter la redondance
def generate_tpm(input_path, output_path, prefix):
        """
        Génère un fichier TPM à partir d'un fichier de comptage brut.
        Args:
        input_path (str): Chemin du fichier d'entrée contenant les KO counts.
        output_path (str): Chemin du fichier de sortie où enregistrer les TPM.
        prefix (str): Préfixe des colonnes correspondant aux KO.
        """
        # Charger les KO counts
        df = pd.read_parquet(input_path)

        # Vérifier et corriger l'index
        if df.index.name != "index":
            df = df.set_index("index") if "index" in df.columns else df

        # Sélectionner les colonnes des KO numbers
        ko_columns = [col for col in df.columns if col.startswith(prefix)]

        if not ko_columns:
            print(f"⚠️ Aucune colonne trouvée avec le préfixe '{prefix}' dans {input_path}")
            return

        # Extraire uniquement les colonnes KO
        df_ko = df[ko_columns]

        # Calculer le TPM
        df_tpm = df_ko.div(df_ko.sum(axis=0), axis=1) * 1e6  # Normalisation en TPM

        # Sauvegarder le fichier TPM
        df_tpm.to_parquet(output_path)
        print(f"✅ Fichier TPM enregistré sous '{output_path}'")

# Génération des fichiers TPM pour MT et MG
generate_tpm(
        input_path='/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/mt_feature_table.parquet',
        output_path='/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/ko_mt_tpm.parquet',
        prefix="mt_mantis_kegg_ko_counts_"
        )

generate_tpm(
        input_path='/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/mg_feature_table.parquet',
        output_path='/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/ko_mg_tpm.parquet',
        prefix="mg_mantis_kegg_ko_counts_"
        )


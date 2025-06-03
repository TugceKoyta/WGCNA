import os
import glob
import gzip
import pandas as pd
from scipy.stats import zscore

# === PARAMÈTRES ===
base_path = "/scratch/users/tkoytaviloglu/data/raw/MUST/"
output_file = "/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/residuals/sample_covariates.tsv"

# === 📥 Collecte des fichiers all_stats
stat_files = glob.glob(os.path.join(base_path, "*/Stats/all_stats.tsv.gz"))

samples_data = []

for stat_file in stat_files:
    # Extraction du nom d’échantillon (ex: M01-01)
    parts = stat_file.split("/")
    family = parts[-3]
    sample = parts[-2]
    sample_id = f"{family}-{sample}"
    try:
        with gzip.open(stat_file, "rt") as f:
            df = pd.read_csv(f, sep="\t", header=0)

            # Récupération reads bruts (ex: mg.r1.fq)
            lib_raw = df.loc[df['step'] == 'mg.r1.fq', 'number'].values
            lib_raw = float(lib_raw[0]) if len(lib_raw) > 0 else None

            # Récupération reads filtrés (mg.r1...filtered... + mg.se...filtered)
            filt1 = df.loc[df['step'].str.endswith('r1.trimmed__PhiX_filtered__T2TCHM13v2genomic_filtered__T2TCHM13v2tra
nscripts_filtered.fq'), 'number']
            filt2 = df.loc[df['step'].str.endswith('se.trimmed__PhiX_filtered__T2TCHM13v2genomic_filtered__T2TCHM13v2tra
nscripts_filtered.fq'), 'number']
            lib_post = 0
            if len(filt1) > 0:
                lib_post += float(filt1.values[0])

            if len(filt2) > 0:
                lib_post += float(filt2.values[0])

            samples_data.append({
                "sample": sample_id,
                "family": family,
                "sex": "NA",  # À remplir manuellement ou automatiquement depuis un fichier pedigree
                "lib_size_raw": lib_raw,
                "lib_size_post": lib_post
                })

    except Exception as e:
        print(f"⚠️ Erreur fichier {sample_id}: {e}")

# === 📊 Création du DataFrame final
df_cov = pd.DataFrame(samples_data)

df_cov["sample"] = df_cov["sample"].str.replace("-Stats", "", regex=False)
df_cov["family"] = df_cov["sample"].str.extract(r"^(M\d{2})")

print("Colonnes disponibles :", df_cov.columns.tolist())
print(df_cov.head())


# Z-score des tailles de librairie
df_cov["lib_size_raw_zscore"] = zscore(df_cov["lib_size_raw"], nan_policy="omit")
df_cov["lib_size_post_zscore"] = zscore(df_cov["lib_size_post"], nan_policy="omit")

# Sauvegarde
df_cov.to_csv(output_file, sep="\t", index=False)
print(f"✅ Fichier généré : {output_file}")

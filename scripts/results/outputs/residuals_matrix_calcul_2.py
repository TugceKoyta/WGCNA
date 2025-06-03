import pandas as pd
import numpy as np
import os
from statsmodels.regression.mixed_linear_model import MixedLM
from tqdm import tqdm

# === 📁 Chemins d’accès ===
base_path = "/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/"
tpm_files = {
            "MG": os.path.join(base_path, "ko_mg_gff_tpm.parquet"),
           # "MT": os.path.join(base_path, "ko_mt_gff_tpm.parquet")
                }
covariate_file = os.path.join(base_path, "residuals/sample_covariates.tsv")
output_dir = os.path.join(base_path, "residuals_corrected_lmer")
os.makedirs(output_dir, exist_ok=True)

# === ⚙️ Boucle principale ===
for label, tpm_path in tpm_files.items():
    print(f"\n🔍 Traitement {label}")
    # Chargement des TPM et transposition
    tpm = pd.read_parquet(tpm_path).T
    tpm.index.name = "sample"
    tpm.index = tpm.index.astype(str)

    # Chargement des covariables
    cov = pd.read_csv(covariate_file, sep="\t")
    cov["sample"] = cov["sample"].astype(str)
    cov = cov.set_index("sample")

    # Restriction aux échantillons communs
    common = tpm.index.intersection(cov.index)
    tpm = tpm.loc[common]
    cov = cov.loc[common]

    # Encodage sex (effet fixe catégoriel)
    cov = cov.copy()
    cov["sex"] = cov["sex"].astype("category")
    # Prendre uniquement les colonnes utiles
    fixed_effects = ["sex", "lib_size_raw_zscore", "lib_size_post_zscore"]
    random_effect = "family"

    # Initialisation des résidus
    residuals_df = pd.DataFrame(index=tpm.index)
    print("📈 Calcul des résidus avec modèle mixte...")

    for ko in tqdm(tpm.columns, desc="KO"):
        print(ko)
        y = tpm[ko]
        if y.isna().any():
            print(f"⚠️ KO {ko} ignoré à cause de NaN")
            continue
        data = cov.copy()
        data["y"] = y
        try:
            md = MixedLM.from_formula("y ~ sex + lib_size_raw_zscore + lib_size_post_zscore",
                    groups=data[random_effect], data=data)
            mdf = md.fit(reml=False, method='lbfgs', warn_convergence=False)
            residuals = mdf.resid
            residuals_df[ko] = residuals
        except Exception as e:
            print(f"⚠️ KO {ko} ignoré (erreur : {e})")
            continue

# Sauvegarde des résidus
out_path = os.path.join(output_dir, f"{label.lower()}_residuals_lmer.parquet")
residuals_df.to_parquet(out_path)
print(f"✅ Résidus sauvegardés → {out_path}")

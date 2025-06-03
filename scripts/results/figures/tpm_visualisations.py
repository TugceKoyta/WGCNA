import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# 📁 Chemins vers les fichiers TPM filtrés
base_path = "/scratch/users/tkoytaviloglu/results/outputs/imp_output_test/"
tpm_paths = {
            "MG": os.path.join(base_path, "ko_mg_gff_tpm.parquet"),
            "MT": os.path.join(base_path, "ko_mt_gff_tpm.parquet")
                }

# 📂 Dossier de sortie des plots
plot_output_dir = os.path.join("/scratch/users/tkoytaviloglu/results/figures/tpm_visualisations/")
os.makedirs(plot_output_dir, exist_ok=True)

# 📈 Fonction de visualisation améliorée
def visualize_tpm(tpm_df, label):
    log_tpm = np.log10(tpm_df + 1)

    # 📊 Histogramme
    plt.figure(figsize=(10, 6))
    sns.histplot(log_tpm.values.flatten(), bins=60, color="cornflowerblue", edgecolor="black", kde=True)
    plt.title(f"Histogram of log10(TPM+1) values - {label}", fontsize=14)
    plt.xlabel("log10(TPM + 1)")
    plt.ylabel("Frequence")
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_output_dir, f"{label.lower()}_tpm_histogram.png"), dpi=300)
    plt.close()

    # 📦 Boxplot par KO
    plt.figure(figsize=(12, 6))
    sns.boxplot(
            data=log_tpm.T,
            palette="YlGnBu",
            linewidth=0.6,
            fliersize=3,
            flierprops=dict(marker='o', markerfacecolor='red', markersize=3, linestyle='none')
            )
    plt.xticks([], [])  # Trop de KOs pour afficher les noms
    plt.title(f"Boxplot of log10(TPM+1) per KO - {label}", fontsize=13)
    plt.ylabel("log10(TPM + 1)")
    plt.grid(axis='y', linestyle='--', alpha=0.4)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_output_dir, f"{label.lower()}_tpm_boxplot.png"), dpi=300)
    plt.close()


# 🔁 Appliquer aux deux jeux de données
for label, path in tpm_paths.items():
    print(f"📥 Chargement des données TPM pour {label}")
    df = pd.read_parquet(path)
    print(f"✅ {df.shape[0]} KO x {df.shape[1]} échantillons")

    visualize_tpm(df, label)

print("📊 Visualisations terminées et sauvegardées.")
